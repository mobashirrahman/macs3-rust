//! The Rust BAM reader must reproduce `MACS3.IO.BAM.BAMaccessor` exactly.
//!
//! The golden files in `tests/golden/bam/` were produced by
//! `oracle/dump_bam_region.py`, which calls upstream's own
//! `get_reads_in_region`, and are compared byte-for-byte. Between them they pin
//! every filtering and coordinate rule in `__fw_binary_parse` and
//! `get_reads_in_region`:
//!
//! * which records survive the flag filter (unmapped / QC-fail / secondary /
//!   supplementary) and the paired-read rule (proper pair, mate mapped);
//! * `MAPQ < 1` and `MAPQ == 255` being dropped;
//! * `rightmost` counting only `M`, `D`, `N`, `=`, `X` -- so a soft-clipped or
//!   inserted read has a *shorter* span here than in a pileup, which
//!   `chr1_1300_1500_d1` shows directly (`c_soft` 1400-1420 rather than
//!   1400-1440);
//! * the packed sequence and quality bytes, hex-encoded, plus `n_edits`;
//! * duplicate suppression: `chr1_1990_2060_d1` keeps 1 of 4 identical
//!   alignments and `d4` keeps all 4;
//! * the order, which is file order -- upstream does not sort during the scan.

use std::path::{Path, PathBuf};

use macs_io::bam::BamAccessor;

fn bam() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("workspace root")
        .join("tests/fixtures/bam/reads.bam")
}

fn golden(name: &str) -> String {
    let p = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("workspace root")
        .join("tests/golden/bam")
        .join(name);
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()))
}

/// Render the region exactly as `oracle/dump_bam_region.py` does.
fn dump(path: &Path, chrom: &str, left: u32, right: u32, maxdup: u32) -> String {
    let mut b = BamAccessor::open(path).expect("open BAM");
    let mut s = String::new();
    s.push_str(&format!(
        "### references\t{}\n",
        b.chromosomes()
            .iter()
            .map(|c| String::from_utf8_lossy(c).into_owned())
            .collect::<Vec<_>>()
            .join(",")
    ));
    s.push_str(&format!(
        "### rlengths\t{}\n",
        b.rlengths()
            .iter()
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join(",")
    ));
    let reads = b
        .reads_in_region(chrom.as_bytes(), left, right, maxdup)
        .expect("region query");
    s.push_str(&format!("### n_reads\t{}\n", reads.len()));
    for r in &reads {
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            String::from_utf8_lossy(&r.name),
            String::from_utf8_lossy(&r.chrom),
            r.lpos,
            r.rpos,
            r.strand,
            r.seq.iter().map(|x| format!("{x:02x}")).collect::<String>(),
            r.qual
                .iter()
                .map(|x| format!("{x:02x}"))
                .collect::<String>(),
            r.cigar
                .iter()
                .map(u32::to_string)
                .collect::<Vec<_>>()
                .join(","),
            r.md,
            String::from_utf8_lossy(&r.sequence()),
            r.n_edits(),
            r.length(),
        ));
    }
    s
}

fn check(name: &str, chrom: &str, left: u32, right: u32, maxdup: u32) {
    let got = dump(&bam(), chrom, left, right, maxdup);
    let want = golden(&format!("{name}.tsv"));
    if got != want {
        let gw: Vec<&str> = got.lines().collect();
        let ww: Vec<&str> = want.lines().collect();
        for (i, (a, b)) in gw.iter().zip(ww.iter()).enumerate() {
            assert_eq!(a, b, "{name}: first difference at line {}", i + 1);
        }
        panic!(
            "{name}: line counts differ (got {}, want {})",
            gw.len(),
            ww.len()
        );
    }
}

#[test]
fn chr1_peak_with_maxdup_1() {
    check("chr1_1000_1200_d1", "chr1", 1000, 1200, 1);
}

#[test]
fn chr1_peak_with_maxdup_2() {
    check("chr1_1000_1200_d2", "chr1", 1000, 1200, 2);
}

/// Every CIGAR variant plus every rejected flag class.
#[test]
fn chr1_cigar_and_flag_cases() {
    check("chr1_1300_1500_d1", "chr1", 1300, 1500, 1);
}

#[test]
fn duplicate_suppression_keeps_one() {
    check("chr1_1990_2060_d1", "chr1", 1990, 2060, 1);
}

#[test]
fn duplicate_suppression_can_keep_all_four() {
    check("chr1_1990_2060_d4", "chr1", 1990, 2060, 4);
}

#[test]
fn chr2_peak() {
    check("chr2_3000_3120_d1", "chr2", 3000, 3120, 1);
}

#[test]
fn whole_contig_queries_match() {
    check("chr1_full", "chr1", 0, 20000, 1);
    check("chr2_full", "chr2", 0, 8000, 1);
}

#[test]
fn an_unknown_chromosome_yields_nothing() {
    let mut b = BamAccessor::open(&bam()).expect("open");
    assert!(b
        .reads_in_region(b"chrNope", 0, 1000, 1)
        .unwrap()
        .is_empty());
}

#[test]
fn a_region_with_no_overlapping_chunk_yields_nothing() {
    let mut b = BamAccessor::open(&bam()).expect("open");
    // past the end of chr2 there is no chunk at all
    assert!(b
        .reads_in_region(b"chr2", 7900, 8000, 1)
        .unwrap()
        .is_empty());
}

#[test]
fn an_unsorted_header_is_rejected() {
    // The sorted check looks only at the first header line.
    let dir = std::env::temp_dir().join(format!("macs3rs_bamsort_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for (i, text) in [
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n",
        "@SQ\tSN:chr1\tLN:100\n@HD\tVN:1.6\tSO:coordinate\n",
    ]
    .into_iter()
    .enumerate()
    {
        let p = dir.join(format!("unsorted{i}.bam"));
        // a body-only file is enough: the header check runs before any record
        let mut b = b"BAM\x01".to_vec();
        b.extend_from_slice(&(text.len() as i32).to_le_bytes());
        b.extend_from_slice(text.as_bytes());
        b.extend_from_slice(&0i32.to_le_bytes());
        std::fs::write(&p, macs_io::bam::bgzf_compress(&b, 65536).unwrap()).unwrap();
        let err = macs_io::bam::read_header(&p).unwrap_err().to_string();
        assert!(err.contains("sorted by coordinates"), "case {i}: {err}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_non_bam_file_is_rejected_with_a_clear_error() {
    let dir = std::env::temp_dir().join(format!("macs3rs_notbam_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("x.bam");
    std::fs::write(&p, b"this is not a BAM file at all").unwrap();
    let err = macs_io::bam::read_header(&p).unwrap_err().to_string();
    assert!(err.contains("not a BAM file"), "{err}");
    let _ = std::fs::remove_dir_all(&dir);
}
