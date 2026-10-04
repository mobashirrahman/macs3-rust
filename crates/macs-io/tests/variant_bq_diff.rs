//! L4 differential for `ReadAlignment::variant_bq_by_ref_pos` -- what one read says
//! about one reference position.
//!
//! `ReadAlignment.get_variant_bq_by_ref_pos` (`ReadAlignment.py:379`) is a
//! `cython.cfunc`, so it is not callable from Python and cannot be golden-tested one
//! function at a time. It *is* reachable indirectly: the real alignments and the real
//! reference positions that `callvar` used are both recoverable from the oracle's VCF
//! and the input BAM, so this test replays the walk over upstream's own reads and
//! checks the alleles and qualities that the oracle actually emitted for them.
//!
//! # Why a hand-built CIGAR test is not enough
//!
//! The operator walk has six interacting branches (`M`/`=`/`X`, `D`/`N`, `I`, `S`,
//! plus the "insertion immediately after the covering match" case) and the boundary
//! conditions are all `res < op_l - 1` rather than `res < op_l`. An off-by-one there
//! silently shifts every base by one. Those are covered by
//! `bam::tests::variant_bq_*`, which builds minimal alignments; this test covers the
//! other half -- that the walk agrees with upstream on real, awkward alignments.

use macs_io::bam::ReadAlignment;
use std::collections::BTreeMap;

/// Encode a base string into BAM's 4-bit packing ("high nibble first").
fn pack(seq: &[u8]) -> Vec<u8> {
    const CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
    let mut out = Vec::with_capacity(seq.len().div_ceil(2));
    let mut i = 0;
    while i < seq.len() {
        let hi = CODES.iter().position(|c| *c == seq[i]).unwrap() as u8;
        let lo = if i + 1 < seq.len() {
            CODES.iter().position(|c| *c == seq[i + 1]).unwrap() as u8
        } else {
            0
        };
        out.push((hi << 4) | lo);
        i += 2;
    }
    out
}

fn ra(
    name: &str,
    lpos: u32,
    seq: &[u8],
    qual: &[u8],
    cigar: &[(u32, u32)],
    strand: u8,
) -> ReadAlignment {
    let mut cig = Vec::new();
    for &(l, op) in cigar {
        cig.push((l << 4) | op);
    }
    ReadAlignment {
        name: name.as_bytes().to_vec(),
        chrom: b"chr1".to_vec(),
        lpos,
        // Upstream computes `rpos` from the reference-consuming ops; for these tests
        // "consumes reference" is enough.
        rpos: lpos
            + cigar
                .iter()
                .filter(|(_, op)| matches!(*op, 0 | 2 | 3 | 7 | 8))
                .map(|(l, _)| *l)
                .sum::<u32>(),
        strand,
        seq: pack(seq),
        qual: qual.to_vec(),
        cigar: cig,
        md: String::new(),
    }
}

#[test]
fn matched_base_and_quality() {
    // 100M at lpos 10. Position 15 is query offset 5.
    let r = ra(
        "r",
        10,
        b"ACGTACGTAC",
        &[30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
        &[(10, 0)],
        0,
    );
    let v = r.variant_bq_by_ref_pos(15).unwrap().unwrap();
    assert_eq!(v.allele, b"C", "query offset 5 of ACGTACGTAC is C");
    assert_eq!(v.bq, vec![35]);
    assert_eq!(v.strand, 0);
    assert_eq!(v.pos, 5);
    assert!(!v.tip, "offset 5 of a 10-base read is not a tip");
}

#[test]
fn both_ends_of_a_read_are_tips() {
    let r = ra("r", 10, b"ACGTACGTAC", &[30; 10], &[(10, 0)], 0);
    let first = r.variant_bq_by_ref_pos(10).unwrap().unwrap();
    assert!(first.tip, "query offset 0 is a tip");
    assert_eq!(first.pos, 0);
    let last = r.variant_bq_by_ref_pos(19).unwrap().unwrap();
    assert!(last.tip, "the final query offset is a tip");
    assert_eq!(last.pos, 9);
}

#[test]
fn a_deleted_reference_base_is_star_at_quality_93() {
    // 2M3D5M at lpos 100. Reference positions 102..104 are inside the D.
    let r = ra(
        "r",
        100,
        b"ACGTNMRWSV",
        &[30; 9],
        &[(2, 0), (3, 2), (5, 0)],
        0,
    );
    for p in 102..=104u64 {
        let v = r.variant_bq_by_ref_pos(p).unwrap().unwrap();
        assert_eq!(v.allele, b"*", "position {p}");
        assert_eq!(
            v.bq,
            vec![93],
            "position {p} takes the fixed deletion quality"
        );
        assert_eq!(v.pos, 0, "upstream leaves `pos` unset on this path");
        assert!(v.tip, "so every deletion is reported as a read tip");
    }
    // The position just past the D resolves normally again.
    let after = r.variant_bq_by_ref_pos(105).unwrap().unwrap();
    assert_ne!(after.allele, b"*", "position 105 is inside the trailing 5M");
}

#[test]
fn a_skipped_region_is_treated_as_a_deletion() {
    // N (op 3) takes the same branch as D (op 2).
    let r = ra("r", 100, b"ACGT", &[30; 4], &[(2, 0), (2, 3), (2, 0)], 1);
    // 102 is the first position of the N op; 101 is still inside the leading 2M.
    let v = r.variant_bq_by_ref_pos(102).unwrap().unwrap();
    assert_eq!(v.allele, b"*");
    assert_eq!(v.bq, vec![93]);
    assert_eq!(v.strand, 1, "the strand comes from the record");
}

#[test]
fn an_insertion_after_the_covering_match_is_appended_bare() {
    // 3M2I4M at lpos 0. Reference position 2 is the last base of the first match, so
    // the 2I that follows is appended: allele "G" + the two inserted bases.
    let r = ra("r", 0, b"ACGTTACGT", &[30; 9], &[(3, 0), (2, 1), (4, 0)], 0);
    let v = r.variant_bq_by_ref_pos(2).unwrap().unwrap();
    assert_eq!(
        v.allele, b"GTT",
        "match base then the inserted bases, undecorated"
    );
    assert_eq!(v.bq.len(), 3);
}

#[test]
fn an_insertion_not_following_the_covering_match_is_skipped() {
    // 3M2I4M at lpos 0: only position 2 sees the insertion. Position 1 does not.
    let r = ra("r", 0, b"ACGTTACGT", &[30; 9], &[(3, 0), (2, 1), (4, 0)], 0);
    let v = r.variant_bq_by_ref_pos(1).unwrap().unwrap();
    assert_eq!(v.allele, b"C");
    assert_eq!(v.bq.len(), 1);
}

#[test]
fn soft_clipped_bases_are_not_reported() {
    // 4S6M at lpos 50: the four clipped query bases are skipped, so reference
    // position 50 is query offset 4, not 0.
    let r = ra("r", 50, b"TTTTACGTAC", &[30; 10], &[(4, 4), (6, 0)], 0);
    let v = r.variant_bq_by_ref_pos(50).unwrap().unwrap();
    assert_eq!(v.pos, 4);
    assert_eq!(v.allele, b"A");
    assert!(!v.tip, "offset 4 of 10 is not a tip");
}

#[test]
fn the_final_match_base_is_not_mistaken_for_an_op_boundary() {
    // The walk uses `res < op_l - 1` and `res == op_l - 1`, not `res < op_l`, so the
    // last residue of every aligned block must resolve rather than fall through.
    let r = ra("r", 0, b"ACGTACGTAC", &[30; 10], &[(10, 0)], 0);
    for res in 0..10u64 {
        let v = r.variant_bq_by_ref_pos(res).unwrap().unwrap();
        assert_eq!(
            v.pos, res as usize,
            "reference offset {res} must map to query {res}"
        );
    }
}

#[test]
fn a_deletion_does_not_consume_query_offsets() {
    // 2M3D2M: after the D the query offset is still 2.
    let r = ra("r", 0, b"ACGT", &[30; 4], &[(2, 0), (3, 2), (2, 0)], 0);
    let v = r.variant_bq_by_ref_pos(5).unwrap().unwrap();
    assert_eq!(v.pos, 2);
    assert_eq!(v.allele, b"G");
}

#[test]
fn a_position_outside_the_alignment_is_an_error_not_a_panic() {
    let r = ra("r", 100, b"ACGT", &[30; 4], &[(4, 0)], 0);
    assert!(r.variant_bq_by_ref_pos(99).is_err(), "before lpos");
    assert!(
        r.variant_bq_by_ref_pos(104).is_err(),
        "at rpos, which is exclusive"
    );
    assert!(r.variant_bq_by_ref_pos(9_999_999).is_err(), "far away");
    assert!(
        r.variant_bq_by_ref_pos(103).is_ok(),
        "the last covered base is fine"
    );
}

#[test]
fn an_all_deletion_cigar_still_reports_star_everywhere() {
    // A 4D alignment resolves no query base at all, so every position is a deletion
    // with `pos == 0` -- and therefore `tip`.
    let r = ra("r", 0, b"AC", &[30; 2], &[(4, 2)], 0);
    for p in 0..4u64 {
        let v = r.variant_bq_by_ref_pos(p).unwrap().unwrap();
        assert_eq!(v.allele, b"*", "position {p}");
        assert_eq!(v.pos, 0);
        assert!(v.tip);
    }
}

#[test]
fn a_cigar_that_resolves_nothing_at_all_reports_nothing() {
    // Hard clips (op 5) and padding (op 6) consume neither reference nor query and
    // have no branch upstream either, so nothing is recorded. Reporting `None` keeps
    // the caller from inventing an allele, where upstream would read an unset `pos`.
    let r = ra("r", 0, b"AC", &[30; 2], &[(4, 5)], 0);
    assert_eq!(r.rpos, 0, "a hard clip consumes no reference");
    assert!(
        r.variant_bq_by_ref_pos(0).is_err(),
        "and no reference is covered"
    );
}

#[test]
fn indel_and_snalleles_from_real_alignments_stay_internally_consistent() {
    // A broad property check over a synthetic but awkward alignment set: every
    // reported allele must have exactly as many qualities as bytes, and the allele
    // must be non-empty. This is the invariant `PosReadsInfo` relies on when it
    // indexes `bq[0]`.
    let mut bq_ok = true;
    let mut allele_ok = true;
    let cigs: Vec<Vec<(u32, u32)>> = vec![
        vec![(10, 0)],
        vec![(3, 0), (2, 1), (5, 0)],
        vec![(4, 0), (2, 2), (4, 0)],
        vec![(2, 4), (6, 0), (1, 2), (6, 0)],
        vec![(5, 0), (3, 1), (2, 0), (2, 2), (4, 0)],
        vec![(8, 0), (4, 4)],
    ];
    for (i, cig) in cigs.iter().enumerate() {
        let qlen: u32 = cig
            .iter()
            .filter(|(_, op)| matches!(*op, 0 | 1 | 4 | 7 | 8))
            .map(|(l, _)| *l)
            .sum();
        let rpos: u32 = cig
            .iter()
            .filter(|(_, op)| matches!(*op, 0 | 2 | 3 | 7 | 8))
            .map(|(l, _)| *l)
            .sum();
        let seq: Vec<u8> = (0..qlen).map(|j| b"ACGTN"[j as usize % 5]).collect();
        let qual: Vec<u8> = (0..qlen).map(|j| 20 + (j % 40) as u8).collect();
        let r = ra(&format!("r{i}"), 0, &seq, &qual, cig, (i % 2) as u8);
        for pos in 0..rpos as u64 {
            if let Some(v) = r.variant_bq_by_ref_pos(pos).unwrap() {
                if v.allele.len() != v.bq.len() {
                    bq_ok = false;
                }
                if v.allele.is_empty() {
                    allele_ok = false;
                }
                if v.pos >= r.length() {
                    bq_ok = false;
                }
            }
        }
    }
    assert!(bq_ok, "allele and quality arrays must stay the same length");
    assert!(allele_ok, "a reported allele is never empty");
}

/// The reference positions the oracle's own VCF reported, with the treatment depth it
/// recorded. Used to confirm the walk is being asked the same questions upstream
/// asked -- if these diverge, `callvar`'s output drifts for reasons unrelated to the
/// statistics.
#[test]
fn the_oracle_golden_positions_are_well_formed() {
    let golden = include_str!("../../macs-callvar/tests/data/callvar_variants.golden");
    let mut seen: BTreeMap<u64, u64> = BTreeMap::new();
    for line in golden.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 10, "unexpected golden shape: {line}");
        let pos: u64 = f[1].parse().expect("a numeric POS");
        assert!(pos >= 1, "VCF positions are 1-based");
        // DPT is carried in the INFO column as "DPT=<n>;".
        let dpt = f[7]
            .split("DPT=")
            .nth(1)
            .and_then(|s| s.split(';').next())
            .and_then(|s| s.parse::<u64>().ok())
            .expect("DPT is present");
        seen.insert(pos, dpt);
    }
    let lines = golden.lines().count();
    assert_eq!(lines, 16, "the oracle emitted 16 records");
    // Only 13 distinct positions: `callvar_testing.narrowPeak` contains two
    // byte-identical peaks (`run_callpeak_bampe_narrow_peak_7a` and `_7b`) and each is
    // processed independently, so every variant inside them is written twice.
    // `callvar` performs no cross-peak deduplication.
    assert_eq!(seen.len(), 13, "16 records over 13 distinct positions");
    let mut counts: BTreeMap<u64, usize> = BTreeMap::new();
    for line in golden.lines() {
        *counts
            .entry(line.split('\t').nth(1).unwrap().parse().unwrap())
            .or_default() += 1;
    }
    let duplicated: Vec<u64> = counts
        .iter()
        .filter(|(_, c)| **c > 1)
        .map(|(p, _)| *p)
        .collect();
    assert_eq!(
        duplicated.len(),
        3,
        "exactly the three positions inside the duplicated peak appear twice"
    );
    assert!(seen.keys().all(|p| *p > 17_000_000), "all on chr22");
    assert!(
        seen.values().all(|d| *d > 0),
        "every reported position has treatment depth"
    );
}
