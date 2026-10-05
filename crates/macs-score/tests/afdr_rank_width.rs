//! The AFDR rank `k` of the callpeak p->q table is a C `long`, not a C `float`.
//!
//! `__cal_pvalue_qvalue_table` declares `k: cython.long`
//! (`CallPeakUnit.py:770`, `CallPeakUnit.c:15662`) and `__pre_computes` declares it the
//! same way (`:864`). Above 2^24 a `float` rank has a spacing greater than 1, so adding
//! a base-pair count to it is a no-op; on a human genome the rank passes 2^24 long
//! before the AFDR tail, and every q-score in that tail inherits the quantisation.
//! The near-identical `ScoreTrack.make_pq_table` (`ScoreTrack.py:510`) really does use
//! `cython.float`, which is why the mistake was plausible -- but it is not the function
//! `callpeak` calls.
//!
//! The expectations are the exact `f32` values upstream's walk produces for this
//! histogram, evaluated from the formula at `CallPeakUnit.py:809-835` with an `i64` rank.

use macs_score::{PScoreHistogram, PqTable};

fn histogram(pairs: &[(f32, i64)]) -> PScoreHistogram {
    let mut h = PScoreHistogram::new();
    for &(v, n) in pairs {
        h.add(v, n);
    }
    h
}

/// N = 2.9e9 (a human genome's worth of base pairs), split so that the middle bucket's
/// rank is `1 + 21_000_000` -- past 2^24, where an `f32` rank has a spacing of 2.
fn genome_scale() -> PScoreHistogram {
    let n = 2_900_000_000i64;
    histogram(&[(30.0, 21_000_000), (2.5, n - 21_000_000 - 1), (0.5, 1)])
}

#[test]
fn the_rank_is_exact_above_two_to_the_24() {
    let t = PqTable::from_histogram(&genome_scale());
    // k = 1: q = 30 + (log10(1) - log10(N)), i.e. the f32 nearest 20.537601
    assert_eq!(t.get(30.0).unwrap().to_bits(), 0x41a4_4d02);
    // k = 21_000_001, where an f32 rank would have snapped to 2^24 * 1.25 and lost
    // every subsequent increment. Upstream: the f32 nearest 0.35982173681259155;
    // with an f32 rank: 0x3eb8_3a8d, the f32 nearest 0.35982170701026917.
    assert_eq!(t.get(2.5).unwrap().to_bits(), 0x3eb8_3a8e);
    // F18: the lowest bucket of the histogram is always forced to 0
    assert_eq!(t.get(0.5), Some(0.0));
}

/// The same table with a rank that has been rounded to `f32`, i.e. what the walk would
/// report if `k` were a C float. Kept as a guard: if the two ever collapse, the
/// assertions above have stopped testing anything.
#[test]
fn the_f32_rank_would_differ() {
    let n = 2_900_000_000f64;
    let f: f32 = (-n.log10()) as f32;
    let k = 21_000_001i64;
    let i64_rank = (2.5 + ((k as f64).log10() + f64::from(f))) as f32;
    let f32_rank = (2.5 + ((k as f32).log10() as f64 + f64::from(f))) as f32;
    assert_ne!(
        i64_rank.to_bits(),
        f32_rank.to_bits(),
        "the two rank widths must be distinguishable at this scale"
    );
}
