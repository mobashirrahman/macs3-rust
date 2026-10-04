//! MACS3-compatible p-scores, q-scores and the AFDR p-to-q table.
//!
//! Port of `MACS3/Signal/ScoreTrack.py:430-544` (`compute_pvalue`,
//! `compute_qvalue`, `make_pq_table`) and the `get_pscore` cache at `:68`.
//!
//! # This is the most compatibility-sensitive module in the project
//!
//! A p-score is a discontinuous function of the local lambda, a q-score is a
//! discontinuous function of the whole genome's p-score histogram, and a peak
//! boundary is a discontinuous function of a q-score. So a one-ULP difference
//! here moves a peak. Every arithmetic step below is therefore transcribed, not
//! improved, and the specific quirks are enumerated:
//!
//! * **p-scores are multiples of 1e-5.** `poisson_cdf(..., log10=True)` rounds
//!   to five decimals, so the score space is a discrete lattice. Nothing
//!   between two lattice points can ever be produced.
//! * **the rank `k` accumulates in `float32`.** Upstream declares
//!   `k: cython.float` and adds a `cython.long` base-pair count to it. Above
//!   2^24 the f32 spacing exceeds 1, so on a human genome the rank is quantised
//!   and `log10(k)` inherits that quantisation.
//! * **`-0.0` and `0.0` share one histogram bucket.** Upstream keys
//!   `pvalue_stat` with a Python `float`, and `hash(-0.0) == hash(0.0)` with
//!   `-0.0 == 0.0`, so the two collapse. `macs-rle::BitKey` deliberately does
//!   *not* collapse them (see its docs), so [`PScoreHistogram`] canonicalises
//!   `-0.0` to `0.0` on insertion. Without this, a genome with many
//!   zero-scoring bases splits into two buckets and every q-value changes.
//! * **chromosomes are summed in lexicographic order** for the histogram, which
//!   matters because the accumulation is not associative in `f32`.
//! * **the `f32` pseudocount is added to both sides**, and the treatment side is
//!   then truncated to an integer.

#![forbid(unsafe_code)]
#![deny(missing_debug_implementations)]

mod pqtable;
mod pscore;
mod scoretrack;

pub use pqtable::{qscore_track, PqTable, QScoreSink};
pub use pscore::{pscore, pscore_track, pseudocounted_inputs, PScoreCache, PScoreHistogram};
pub use scoretrack::{
    log_lr_asym, log_lr_sym, DiffPeak, NormMethod, ScoreMethod, ScoreTrack, ScoreTrack2, TwoScores,
};

use macs_core::Coord;

/// The MACS3 pseudocount, in `-log10` score space.
///
/// Upstream's default is `1.0` and it is added to the treatment count (before
/// truncation to an integer) and to the local lambda (as a float).
pub const DEFAULT_PSEUDOCOUNT: f32 = 1.0;

/// Canonicalise a `-log10` score for use as a histogram key.
///
/// Maps `-0.0` to `+0.0` so that the two share a bucket, matching Python's dict
/// semantics. Every other value is returned unchanged, bit for bit.
#[inline]
pub fn canonical_score_key(v: f32) -> u32 {
    if v == 0.0 {
        0.0f32.to_bits()
    } else {
        v.to_bits()
    }
}

/// Recover a score from [`canonical_score_key`].
#[inline]
pub fn score_from_key(k: u32) -> f32 {
    f32::from_bits(k)
}

/// The number of base pairs upstream counts in `pvalue_stat`.
///
/// The convention is the right-endpoint one (see `F5` in
/// `docs/upstream-findings.md`): a run of length `n` contributes `n`.
#[inline]
pub fn run_length(start: Coord, end: Coord) -> u64 {
    end.saturating_sub(start)
}
