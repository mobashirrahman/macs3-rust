//! L2 property invariants for the p-score and the p->q table.
//!
//! The acceptance criteria bound the p- and q-score errors numerically (1e-9 and
//! 1e-6). Those bounds are what the recorded tables in `pq_vectors.rs` check. These
//! tests instead check the *shape* of the two functions -- monotonicity and sign
//! conventions -- which is what a wrong `logspace_add` or a flipped comparison
//! produces, and which a fixed corpus can miss whenever the corpus happens not to
//! exercise the offending region.

use macs_core::ChromId;
use macs_rle::SignalTrack;
use macs_score::{pscore, PScoreHistogram, PqTable};
use proptest::prelude::*;

/// Build a histogram whose lengths telescope to `total` (positive, as they always
/// are in practice: the lengths are run extents up to the final paired position).
fn telescoping(parts: &[(f32, i64)]) -> PScoreHistogram {
    let mut h = PScoreHistogram::new();
    let mut v: Vec<(f32, i64)> = parts.to_vec();
    v.sort_by_key(|(_, n)| *n);
    for (s, n) in v {
        h.add(s, n);
    }
    h
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// The p-score never decreases as the observed depth grows: deeper piles are
    /// never less significant.
    #[test]
    fn pscore_is_monotone_in_observation(
        lam in 0.1f32..50.0,
        obs in prop::collection::vec(0u32..500, 2..12),
    ) {
        for pair in obs.windows(2) {
            let (lo, hi) = (pair[0], pair[1]);
            let (a, b) = (pscore(None, lo, lam), pscore(None, hi, lam));
            if hi >= lo {
                prop_assert!(b >= a,
                    "pscore fell from {} to {} for lam={} ({} -> {})", a, b, lam, lo, hi);
            }
        }
    }

    /// A larger local lambda never raises the p-score: the same observation is less
    /// significant against a more competitive background.
    #[test]
    fn pscore_is_antitone_in_expectation(
        obs in 0u32..400,
        lo_lam in 0.05f32..10.0,
        hi_lam in 0.05f32..10.0,
    ) {
        let (a, b) = (lo_lam, hi_lam);
        let (small, big) = if a <= b { (a, b) } else { (b, a) };
        let (ps, pb) = (pscore(None, obs, small), pscore(None, obs, big));
        prop_assert!(pb <= ps,
            "lambda {} gave a higher p-score ({}) than lambda {} ({})", big, pb, small, ps);
    }

    /// `pscore` is a finite, non-negative number for every in-range input, and an
    /// expectation of exactly zero yields the maximum score (the observation is
    /// infinitely unlikely against an empty background).
    #[test]
    fn pscore_is_always_finite_and_nonnegative(
        obs in 0u32..2000,
        lam in 0.0f32..1000.0,
    ) {
        let p = pscore(None, obs, lam);
        prop_assert!(p.is_finite(), "non-finite p-score {} for {} / {}", p, obs, lam);
        prop_assert!(p >= 0.0, "negative p-score {}", p);
    }

    /// The q-table reports the same number of scored base pairs as its histogram.
    ///
    /// The signed total is what the cutoff analysis uses to size the above-cutoff
    /// set, so a lost or double-counted span here silently changes every -log10(q).
    #[test]
    fn qtable_preserves_the_scored_total(
        parts in prop::collection::vec((0.0f32..500.0, 1i64..1000), 1..20),
    ) {
        let h = telescoping(&parts);
        let q = PqTable::from_histogram(&h);
        prop_assert_eq!(q.total(), h.total().unsigned_abs(),
            "q-table total {} != histogram {}", q.total(), h.total());
    }

    /// q is monotone *non-decreasing* in the p-score.
    ///
    /// Sign convention, spelled out because getting it backwards silently inverts
    /// every `-q` cutoff: MACS's `qscore` is already `-log10(q)`, not `q` itself.
    /// `PeakDetect` therefore tests `qscore >= cutoff`, i.e. a *larger* value means
    /// *more* significant. Since the p-score is likewise `-log10(p)`, a more
    /// significant observation must produce both a larger p-score and a larger
    /// q-score -- the two run in the same direction.
    #[test]
    fn qscore_is_monotone_in_pscore(
        parts in prop::collection::vec((10.0f32..400.0, 1i64..500), 2..15),
    ) {
        let h = telescoping(&parts);
        let q = PqTable::from_histogram(&h);
        let mut scores: Vec<f32> = h.iter().map(|(v, _)| v).collect();
        scores.sort_by(|l, r| l.total_cmp(r));
        scores.dedup();
        for w in scores.windows(2) {
            let (lo, hi) = (w[0], w[1]);      // lo < hi, so lo is the *smaller* p-score
            let (qlo, qhi) = (q.qscore_or_zero(lo), q.qscore_or_zero(hi));
            if qlo.is_nan() || qhi.is_nan() {
                continue; // the F177 NaN window; covered separately
            }
            prop_assert!(qhi >= qlo - 1e-4,
                "q not monotone: p {} -> q {}, p {} -> q {}", lo, qlo, hi, qhi);
        }
    }

    /// A q-table built from a track agrees with one built from the same counts added
    /// by hand. `add_track_from` is the only path the parallel code uses, so it must
    /// reproduce the histogram exactly.
    #[test]
    fn add_track_from_matches_manual_counts(
        vals in prop::collection::vec(0.0f32..100.0, 1..25),
        step in 1u64..17,
    ) {
        let mut t = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
        for (i, v) in vals.iter().enumerate() {
            t.push_exact((i as u64 + 1) * step, *v);
        }
        let mut from_track = PScoreHistogram::new();
        from_track.add_track(&t);
        let mut manual = PScoreHistogram::new();
        let mut prev = 0i64;
        for (i, v) in vals.iter().enumerate() {
            let end = ((i as u64 + 1) * step) as i64;
            manual.add(*v, end - prev);
            prev = end;
        }
        prop_assert_eq!(from_track.total(), manual.total());
        let a = PqTable::from_histogram(&from_track);
        let b = PqTable::from_histogram(&manual);
        for v in vals.iter() {
            prop_assert_eq!(from_track.count(*v), manual.count(*v), "count differs at {}", v);
            prop_assert_eq!(a.qscore_or_zero(*v), b.qscore_or_zero(*v));
        }
    }

    /// An empty histogram produces an empty table that answers every lookup with 0
    /// rather than panicking or returning NaN.
    #[test]
    fn empty_histogram_is_a_zero_table(v in 0.0f32..1000.0) {
        let q = PqTable::from_histogram(&PScoreHistogram::new());
        prop_assert_eq!(q.total(), 0);
        prop_assert_eq!(q.qscore_or_zero(v), 0.0);
    }
}
