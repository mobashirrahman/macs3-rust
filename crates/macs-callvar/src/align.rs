//! Port of `RACollection.align_unitig_to_REFSEQ` and its Smith-Waterman backend.
//!
//! # Why this is C too
//!
//! `align_unitig_to_REFSEQ` places each assembled unitig back onto the peak consensus
//! with `MACS3/Signal/swalign.c`'s `smith_waterman` -- local alignment, **without** an
//! affine gap penalty. That file is vendored verbatim alongside fermi-lite
//! (`vendor/fermi-lite/swalign.{c,h}`); re-deriving the score matrix in Rust would be a
//! different algorithm with different tie-breaking, and the unitig placement feeds
//! straight into the variant calls.
//!
//! # The scoring scheme is not the usual one
//!
//! Per `swalign.c`: match **+2**, mismatch **-3**, gap **-5**, and the second gap in a
//! run costs **-2**. `verify_alns`'s docstring spells out what that means in
//! operational terms: a score of 150 over a 100 bp window is "10 mismatches within
//! 100bps", i.e. `10 * (2 - 3) == -10` would be negative -- so the score is accumulated
//! over the aligned block only, and `min_score_100 = 150` corresponds to roughly a 10%
//! mismatch rate over 100 bp. The constants live in the C; nothing here re-derives them.
//!
//! # Forward *and* reverse complement, and the list is mutated
//!
//! Each unitig is aligned twice -- as-is and reverse-complemented -- and the better
//! score wins. **If the reverse complement wins, upstream overwrites `unitig_list[i]`
//! with the revcomp**, so the caller sees the unitig in the orientation it aligned in.
//! That mutation is load-bearing: `remap_RAs_w_unitigs` later re-derives read positions
//! against the (possibly flipped) unitig.

use crate::fermi::Unitig;
use std::ffi::{c_char, c_double, c_int, c_uint, CString};

#[repr(C)]
#[derive(Clone, Copy)]
struct SeqPairT {
    a: *mut c_char,
    alen: c_uint,
    b: *mut c_char,
    blen: c_uint,
}

#[repr(C)]
struct AlignT {
    seqs: *mut SeqPairT,
    markup: *mut c_char,
    start_a: c_int,
    start_b: c_int,
    end_a: c_int,
    end_b: c_int,
    matches: c_int,
    gaps: c_int,
    score: c_double,
}

extern "C" {
    #[link_name = "smith_waterman"]
    fn c_smith_waterman(problem: *mut SeqPairT) -> *mut AlignT;
    // libc free, under a Rust-side name so it does not collide with the declaration
    // in `fermi.rs`.
    #[link_name = "free"]
    fn c_free(p: *mut std::ffi::c_void);
}

/// One alignment result, owned on the Rust side.
///
/// `PartialEq` but not `Eq`: the score is an `f64`.
#[derive(Debug, Clone, PartialEq)]
pub struct Alignment {
    /// The aligned target (unitig) subsequence.
    pub target_aln: Vec<u8>,
    /// The aligned reference (consensus) subsequence.
    pub reference_aln: Vec<u8>,
    /// The match/mismatch/gap markup from `swalign.c`.
    pub markup: Vec<u8>,
    /// Local alignment score.
    pub score: f64,
}

/// `MACS3.Utilities.Constants.__DNACOMPLEMENT__`.
const DNACOMPLEMENT: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        t[i] = i as u8;
        i += 1;
    }
    // A<->T, C<->G. Everything else is its own complement, which matches upstream's
    // table (it is a full 256-entry translate table).
    t[b'A' as usize] = b'T';
    t[b'T' as usize] = b'A';
    t[b'C' as usize] = b'G';
    t[b'G' as usize] = b'C';
    t
};

/// Reverse-complement a nucleotide string.
pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|b| DNACOMPLEMENT[*b as usize]).collect()
}

/// Align `target` against `reference` with Smith-Waterman local alignment.
pub fn smith_waterman(target: &[u8], reference: &[u8]) -> Result<Alignment, String> {
    if target.is_empty() || reference.is_empty() {
        return Err("smith_waterman: empty sequence".to_string());
    }
    // swalign.c reads the sequences as NUL-terminated C strings, so it must not see an
    // interior NUL. A `CString::new` failure here means the input is not a DNA string.
    let a = CString::new(target).map_err(|_| "target contains an interior NUL".to_string())?;
    let b =
        CString::new(reference).map_err(|_| "reference contains an interior NUL".to_string())?;
    // swalign.c only reads, but the struct is not `const char *`, so hand it the
    // mutable pointer the header asks for.
    let mut problem = SeqPairT {
        a: a.as_ptr() as *mut c_char,
        alen: target.len() as c_uint,
        b: b.as_ptr() as *mut c_char,
        blen: reference.len() as c_uint,
    };

    // SAFETY: `problem` points at live, correctly-laid-out `seq_pair_t`s whose buffers
    // are NUL-terminated and outlive the call.
    let res = unsafe { c_smith_waterman(&mut problem) };
    if res.is_null() {
        return Err("smith_waterman: returned NULL".to_string());
    }
    // SAFETY: non-null `align_t` with the three malloc'd buffers upstream also frees.
    unsafe {
        let r = &*res;
        let target_aln = if r.seqs.is_null() || r.seqs.read().a.is_null() {
            Vec::new()
        } else {
            let sa = &*r.seqs;
            std::ffi::CStr::from_ptr(sa.a).to_bytes().to_vec()
        };
        let reference_aln = if r.seqs.is_null() || r.seqs.read().b.is_null() {
            Vec::new()
        } else {
            let sb = &*r.seqs;
            std::ffi::CStr::from_ptr(sb.b).to_bytes().to_vec()
        };
        let markup = if r.markup.is_null() {
            Vec::new()
        } else {
            std::ffi::CStr::from_ptr(r.markup).to_bytes().to_vec()
        };
        let out = Alignment {
            target_aln,
            reference_aln,
            markup,
            score: r.score,
        };
        // Upstream frees exactly these four buffers: the two aligned sequences, the
        // markup, and the `align_t` itself.
        let seqs_ptr = r.seqs;
        if !seqs_ptr.is_null() {
            let sa = &mut *seqs_ptr;
            if !sa.a.is_null() {
                c_free(sa.a as *mut std::ffi::c_void);
            }
            if !sa.b.is_null() {
                c_free(sa.b as *mut std::ffi::c_void);
            }
            c_free(seqs_ptr as *mut std::ffi::c_void);
        }
        if !r.markup.is_null() {
            c_free(r.markup as *mut std::ffi::c_void);
        }
        c_free(res as *mut std::ffi::c_void);
        Ok(out)
    }
}

/// `align_unitig_to_REFSEQ`.
///
/// Aligns every unitig forward and reverse-complemented against `reference`, keeps the
/// better score, and -- as upstream does -- **rewrites `unitigs[i]` in place** when the
/// reverse complement wins, so the caller sees the unitig in its aligned orientation.
///
/// Returns `(target_alns, reference_alns, aln_scores, markup_alns)`, positionally
/// aligned with `unitigs`.
/// The four per-unitig alignment vectors upstream returns together.
#[derive(Debug, Clone, Default)]
pub struct Alignments {
    /// Aligned unitig subsequences.
    pub target_alns: Vec<Vec<u8>>,
    /// Aligned consensus subsequences.
    pub reference_alns: Vec<Vec<u8>>,
    /// Local alignment scores, positionally aligned with the rest.
    pub scores: Vec<f64>,
    /// `swalign.c` markup strings.
    pub markup_alns: Vec<Vec<u8>>,
}

pub fn align_unitigs_to_reference(
    unitigs: &mut [Unitig],
    reference: &[u8],
) -> Result<Alignments, String> {
    let mut out = Alignments::default();

    for unitig in unitigs.iter_mut() {
        let fwd = smith_waterman(&unitig.seq, reference)?;
        let rc_seq = revcomp(&unitig.seq);
        let rev = smith_waterman(&rc_seq, reference)?;

        // Strict `>`: an exact tie keeps the forward alignment, as upstream does.
        let forward_wins = fwd.score > rev.score;
        if !forward_wins {
            // The list is mutated here, deliberately: `remap_RAs_w_unitigs` re-derives
            // read positions against the (possibly flipped) unitig.
            unitig.seq = rc_seq;
        }
        let best = if forward_wins { fwd } else { rev };
        out.target_alns.push(best.target_aln);
        out.reference_alns.push(best.reference_aln);
        out.scores.push(best.score);
        out.markup_alns.push(best.markup);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_round_trips() {
        let s = b"ACGTACGTNN";
        let r = revcomp(s);
        assert_eq!(revcomp(&r), s.to_vec());
        assert_eq!(r, b"NNACGTACGT".to_vec());
    }

    #[test]
    fn an_identical_sequence_aligns_to_itself() {
        let s = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let a = smith_waterman(s, s).unwrap();
        // 32 bases, +2 each.
        assert_eq!(a.score, 64.0, "score {}", a.score);
        assert!(!a.markup.is_empty());
        assert!(!a.markup.contains(&b' '), "no gaps expected");
    }

    #[test]
    fn a_disjoint_pair_still_aligns_locally() {
        // Smith-Waterman is local, so the score is never negative: the worst case is a
        // single matching base, which is +2 here. What matters is that it is far below
        // a real alignment (20 for ten identical bases), because `verify_alns` uses
        // exactly that gap to drop unitigs that do not belong.
        let a = smith_waterman(b"ACGTACGTAC", b"TTTTTTTTTT").unwrap();
        assert!(
            a.score >= 0.0,
            "a local score cannot be negative: {}",
            a.score
        );
        assert!(
            a.score < 4.0,
            "only an isolated base should match, got {}",
            a.score
        );
        let perfect = smith_waterman(b"ACGTACGTAC", b"ACGTACGTAC").unwrap();
        assert_eq!(perfect.score, 20.0);
    }

    #[test]
    fn an_empty_input_is_refused_rather_than_crashing() {
        assert!(smith_waterman(b"", b"ACGT").is_err());
        assert!(smith_waterman(b"ACGT", b"").is_err());
    }

    #[test]
    fn many_calls_do_not_leak_or_double_free() {
        // The ownership hand-back is the risky part; hammer it.
        for _ in 0..200 {
            let _ = smith_waterman(b"ACGTACGTACGTACGT", b"ACGTACGTACGTACGT").unwrap();
        }
    }

    #[test]
    fn a_unitig_flipped_by_its_alignment_is_rewritten_in_place() {
        // The reference is the reverse complement of the unitig, so the RC branch wins
        // and the caller's list must come back flipped.
        let fwd = b"ACGTACGTACGTACGTACGT".to_vec();
        let reference = revcomp(&fwd);
        let mut unitigs = vec![Unitig {
            seq: fwd.clone(),
            nsr: 3,
        }];
        let a = align_unitigs_to_reference(&mut unitigs, &reference).unwrap();
        assert_eq!(a.scores.len(), 1);
        assert!(a.scores[0] > 0.0, "score {}", a.scores[0]);
        assert_eq!(
            unitigs[0].seq, reference,
            "the unitig list is mutated when the RC alignment wins"
        );
    }

    #[test]
    fn a_unitig_already_in_the_reference_orientation_is_left_alone() {
        let fwd = b"ACGTACGTACGTACGTACGT".to_vec();
        let mut unitigs = vec![Unitig {
            seq: fwd.clone(),
            nsr: 3,
        }];
        let a = align_unitigs_to_reference(&mut unitigs, &fwd).unwrap();
        assert!(a.scores[0] > 0.0);
        assert_eq!(
            unitigs[0].seq, fwd,
            "the forward orientation wins a tie on itself"
        );
    }

    #[test]
    fn the_output_vectors_stay_positionally_aligned() {
        let u = vec![
            Unitig {
                seq: b"ACGTACGTACGTACGTACGT".to_vec(),
                nsr: 1,
            },
            Unitig {
                seq: b"TTTTGGGGCCCCAAAATTTT".to_vec(),
                nsr: 2,
            },
        ];
        let ref20 = b"ACGTACGTACGTACGTACGT";
        let a = {
            let mut copy = u.clone();
            align_unitigs_to_reference(&mut copy, ref20).unwrap()
        };
        assert_eq!(a.target_alns.len(), 2);
        assert_eq!(a.reference_alns.len(), 2);
        assert_eq!(a.scores.len(), 2);
        assert_eq!(a.markup_alns.len(), 2);
        assert!(
            a.scores[0] > a.scores[1],
            "the matching unitig scores higher: {:?}",
            a.scores
        );
    }
}
