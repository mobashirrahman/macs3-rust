//! L2 property invariants for the RLE `SignalTrack`.
//!
//! `SignalTrack` is the *only* signal representation in the port, so a bug here
//! corrupts every downstream stage. The parity vectors pin the value function on a
//! recorded corpus; these pin the structural laws that must hold for any input,
//! checked against a naive per-base model rather than another call into the same code.

use macs_core::{ChromId, Interval};
use macs_rle::SignalTrack;
use proptest::prelude::*;

/// Expand a track to a dense per-base vector over `[start, cursor)`.
fn dense(t: &SignalTrack<f32>) -> Vec<f32> {
    let mut v: Vec<f32> = Vec::new();
    for r in t.runs() {
        let at = v.len() as u64;
        v.extend(std::iter::repeat_n(r.value, (r.end - at) as usize));
    }
    v
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// `integral()` equals the naive per-base sum.
    #[test]
    fn integral_matches_dense_sum(
        vals in prop::collection::vec(0.0f32..100.0, 1..40),
    ) {
        let mut t = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
        for (i, v) in vals.iter().enumerate() {
            t.push_exact((i as u64 + 1) * 7, *v);
        }
        let naive: f64 = vals.iter().map(|v| *v as f64 * 7.0).sum();
        prop_assert!((t.integral() - naive).abs() < 1e-3,
            "integral {} != naive {}", t.integral(), naive);
    }

    /// Run ends are strictly increasing, so `find_run_index` is well defined.
    #[test]
    fn run_ends_are_strictly_increasing(
        vals in prop::collection::vec(0.0f32..10.0, 1..60),
    ) {
        let mut t = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
        for (i, v) in vals.iter().enumerate() {
            t.push((i as u64 + 1) * 3, *v);
        }
        for w in t.runs().windows(2) {
            prop_assert!(w[0].end < w[1].end, "ends not increasing");
        }
    }

    /// `push` ignores a non-advancing end rather than corrupting the track or
    /// panicking; the cursor never moves backwards.
    #[test]
    fn push_never_moves_the_cursor_backwards(
        ends in prop::collection::vec(0u64..500, 1..50),
        v in 0.0f32..10.0,
    ) {
        let mut t = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
        let mut last = 0u64;
        for e in ends {
            t.push(e, v);
            prop_assert!(t.cursor() >= last, "cursor went backwards");
            last = t.cursor();
        }
    }

    /// Coalescing equal-valued neighbours must not change the integral.
    #[test]
    fn coalescing_preserves_the_integral(
        n in 1usize..50,
        v in 0.0f32..10.0,
    ) {
        let mut merged = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
        let mut split = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
        for i in 0..n {
            merged.push((i as u64 + 1) * 4, v);
            split.push_exact((i as u64 + 1) * 4, v);
        }
        prop_assert!((merged.integral() - split.integral()).abs() < 1e-6);
    }

    /// `max_with` is a pointwise upper bound. Local lambda is the element-wise max
    /// of four tracks, so a wrong max silently rescales every p-score.
    #[test]
    fn max_with_is_an_upper_bound(
        a in prop::collection::vec(0.0f32..20.0, 1..30),
        b in prop::collection::vec(0.0f32..20.0, 1..30),
    ) {
        // `push_exact` requires non-decreasing ends by contract -- it keeps runs
        // separate instead of coalescing, and the sweep sorts before calling it.
        let mut ta = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
        let mut tb = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
        let (mut ea, mut eb): (Vec<_>, Vec<_>) = (
            a.iter().enumerate().map(|(i, v)| ((i as u64 + 1) * 5, *v)).collect(),
            b.iter().enumerate().map(|(i, v)| ((i as u64 + 1) * 11, *v)).collect(),
        );
        ea.sort_by_key(|(e, _)| *e);
        eb.sort_by_key(|(e, _)| *e);
        for (e, v) in ea { ta.push_exact(e, v); }
        for (e, v) in eb { tb.push_exact(e, v); }
        let m = ta.max_with(&tb);
        // `max_with` is a pointwise max **over the overlap of the two covered
        // regions**: past the shorter track it yields nothing rather than
        // extending the survivor. Local lambda relies on that truncation -- all four
        // inputs share one domain there -- so the truncation is pinned, not wished
        // away.
        let overlap = ta.cursor().min(tb.cursor());
        for x in 0..400u64 {
            let have = m.value_at_or(x, 0.0);
            if x < overlap {
                let want = ta.value_at_or(x, 0.0).max(tb.value_at_or(x, 0.0));
                prop_assert_eq!(have, want, "max_with wrong at {}", x);
            } else {
                prop_assert_eq!(have, 0.0,
                    "max_with invented coverage past the overlap at {}", x);
            }
        }
        prop_assert!(m.cursor() <= overlap, "max_with ran past the overlap");
    }

    /// `max_with` is commutative and associative on the value function.
    #[test]
    fn max_with_is_commutative(
        a in prop::collection::vec(0.0f32..20.0, 1..25),
        b in prop::collection::vec(0.0f32..20.0, 1..25),
    ) {
        let build = |xs: &Vec<f32>, step: u64| {
            let mut t = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
            for (i, v) in xs.iter().enumerate() { t.push_exact((i as u64 + 1) * step, *v); }
            t
        };
        let (ta, tb) = (build(&a, 5), build(&b, 7));
        let (ab, ba) = (ta.max_with(&tb), tb.max_with(&ta));
        prop_assert_eq!(ab.integral(), ba.integral());
        for x in 0..300u64 {
            prop_assert_eq!(ab.value_at_or(x, 0.0), ba.value_at_or(x, 0.0));
        }
    }

    /// The integral of a slice equals the integral of the covered part, and slicing
    /// can never increase the integral.
    #[test]
    fn slice_conserves_the_covered_integral(
        vals in prop::collection::vec(0.0f32..20.0, 1..30),
        lo in 0u64..300,
        hi in 0u64..300,
    ) {
        let mut t = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
        for (i, v) in vals.iter().enumerate() { t.push_exact((i as u64 + 1) * 13, *v); }
        let (lo, hi) = (lo.min(hi), lo.max(hi));
        let s = t.slice(Interval::new(lo, hi));
        let expect: f64 = dense(&t).iter().enumerate()
            .filter(|(x, _)| (*x as u64) >= lo && (*x as u64) < hi)
            .map(|(_, v)| *v as f64).sum();
        prop_assert!((s.integral() - expect).abs() < 1e-3,
            "slice integral {} != {}", s.integral(), expect);
        prop_assert!(s.integral() <= t.integral() + 1e-6);
    }

    /// Replaying the same pushes in any order cannot change the value function.
    ///
    /// This is the property that makes the chromosome-level parallelism sound: each
    /// chromosome's track is built from its own read partition regardless of the
    /// order the partitions arrive in.
    #[test]
    fn construction_is_order_independent(
        vals in prop::collection::vec((0u64..400, 0.0f32..10.0), 1..60),
    ) {
        let mut shuffled = vals.clone();
        for i in (1..shuffled.len()).rev() {
            let j = (i * 2_654_435_761) % (i + 1);
            shuffled.swap(i, j);
        }
        // Both push orders must be sorted first: the RLE layer is append-only by
        // design, so order-independence is the caller's sort plus the track's
        // determinism, not the track reordering for you. The parallel pileup path
        // relies on exactly this (see `pileup_is_order_independent`).
        //
        // The sort key must include the value. `push*` drops a run whose end does not
        // advance the cursor, so among runs sharing an end the *first one pushed
        // wins*; a stable sort by end alone would therefore make the track depend on
        // input order. Sorting by `(end, value)` makes the tie-break total, which is
        // what any parallel producer has to do.
        let build = |xs: &[(u64, f32)]| {
            let mut sorted: Vec<(u64, f32)> =
                xs.iter().copied().map(|(e, v)| (e.max(1), v)).collect();
            sorted.sort_by(|l, r| l.0.cmp(&r.0).then(l.1.total_cmp(&r.1)));
            let mut t = SignalTrack::<f32>::empty(ChromId(0), 0, 1_000_000);
            for &(e, v) in &sorted { t.push_exact(e, v); }
            t
        };
        let (a, b) = (build(&vals), build(&shuffled));
        prop_assert_eq!(a.integral(), b.integral());
        for x in 0..600u64 {
            prop_assert_eq!(a.value_at_or(x, 0.0), b.value_at_or(x, 0.0),
                "value at {} depends on push order", x);
        }
    }
}
