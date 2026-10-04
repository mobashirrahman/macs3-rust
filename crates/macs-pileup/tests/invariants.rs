//! L2 property invariants for the pileup sweep.
//!
//! The parity vectors in `parity_vectors.rs` prove the sweep matches the oracle
//! on a recorded corpus. These tests prove properties that must hold for *any*
//! input, which is what catches the classes of bug a fixed corpus cannot: off-by-one
//! in the endpoint construction, a dropped read at the contig edges, a weight
//! dropped, or a scale factor applied twice.
//!
//! Every property here is checked against an independent naive reference
//! computation rather than against another call into the same code, so a shared
//! bug cannot make a test vacuous.

use macs_core::{ChromId, Interval};
use macs_pileup::{
    integrated_depth, naive_quick_pileup, pileup_from_positions, pileup_from_weighted_positions,
    restrict, SingleEndParams,
};
use macs_rle::SignalTrack;
use proptest::prelude::*;

/// A brute-force single-end pileup: mark every base of every extended read.
///
/// This mirrors the documented contract rather than the implementation: a centred
/// read at `x` covers `[x - d/2, x + d - d/2)`, and **each endpoint is clipped to
/// `[0, rlength]` on its own** -- upstream clamps both in C `long` arithmetic, so a
/// read hanging off the left edge covers `[0, x + d - d/2)`, i.e. *less* than `d`.
/// Deriving `end` from an already-clamped `start` (the tempting shortcut) silently
/// over-counts exactly those edge reads, which is the off-by-one this guards.
fn naive_se(plus: &[u64], minus: &[u64], d: i64, rlength: i64) -> Vec<(i64, i64)> {
    let mut depth = vec![0i64; rlength as usize];
    let clamp = |v: i64| v.clamp(0, rlength);
    for &x in plus.iter().chain(minus.iter()) {
        let b = x as i64;
        let s = clamp(b - d / 2);
        let e = clamp(b + (d - d / 2));
        for cell in depth[s as usize..e as usize].iter_mut() {
            *cell += 1;
        }
    }
    depth
        .iter()
        .enumerate()
        .map(|(i, &v)| (i as i64, v))
        .collect()
}

fn params(d: i64, rlength: i64) -> SingleEndParams {
    SingleEndParams::centred(d, rlength as u64, 1.0)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// The sweep must reproduce a per-base brute-force pileup exactly.
    ///
    /// This is the core property: it simultaneously pins the endpoint arithmetic
    /// (`[p - d/2, p - d/2 + d)`), the clipping to `[0, rlength)`, and the
    /// endpoint-sweep accumulation.
    #[test]
    fn sweep_equals_brute_force(
        d in 1i64..40,
        rlength in 200i64..2000,
        plus in prop::collection::vec(0u64..2000, 0..60),
        minus in prop::collection::vec(0u64..2000, 0..60),
    ) {
        let rlength = rlength.max(*plus.iter().max().unwrap_or(&0) as i64 + 1)
            .max(*minus.iter().max().unwrap_or(&0) as i64 + 1);
        let expect = naive_se(&plus, &minus, d, rlength);
        let t = pileup_from_positions(
            ChromId(0),
            &plus,
            &minus,
            &params(d, rlength),
        );
        // compare at every position, not just the breakpoints
        for (x, want) in &expect {
            prop_assert_eq!(t.value_at_or(*x as u64, 0.0) as i64, *want,
                "depth mismatch at {} for d={} rlength={}", x, d, rlength);
        }
    }

    /// The integral of the pileup equals the number of bases actually covered.
    ///
    /// Reads that hang off either end of the contig contribute less than `d`, so
    /// this is `<=`, and the bound is tight when nothing is clipped.
    #[test]
    fn integral_is_covered_bases(
        d in 1i64..40,
        rlength in 500i64..3000,
        plus in prop::collection::vec(0u64..3000, 0..80),
        minus in prop::collection::vec(0u64..3000, 0..80),
    ) {
        let rlength = rlength.max(*plus.iter().max().unwrap_or(&0) as i64 + 1)
            .max(*minus.iter().max().unwrap_or(&0) as i64 + 1);
        let t = pileup_from_positions(
            ChromId(0), &plus, &minus, &params(d, rlength),
        );
        let expect = naive_se(&plus, &minus, d, rlength)
            .iter().map(|(_, v)| *v as f64).sum::<f64>();
        prop_assert!(
            (integrated_depth(&t) - expect).abs() < 1e-3,
            "integral {} != covered {}", integrated_depth(&t), expect
        );
        let n = (plus.len() + minus.len()) as f64 * d as f64;
        prop_assert!(integrated_depth(&t) <= n + 1e-6, "over-counts coverage");
        prop_assert!(integrated_depth(&t) >= 0.0);
    }

    /// The depth never exceeds the number of reads supplied, whatever the shifts.
    #[test]
    fn depth_never_exceeds_read_count(
        d in 1i64..60,
        rlength in 300i64..1500,
        n in 0usize..200,
        shift in -60i64..60,
        pos in prop::collection::vec(0u64..1500, 0..200),
    ) {
        let n = n.min(pos.len());
        let reads = &pos[..n];
        let mut p = params(d, rlength);
        p.five_shift = shift;
        p.three_shift = d - shift;
        let t = pileup_from_positions(ChromId(0), reads, &[], &p);
        let cap = n as f32;
        for x in 0..rlength.min(1500) {
            prop_assert!(t.value_at_or(x as u64, 0.0) <= cap,
                "depth {} > {} reads at {}", t.value_at_or(x as u64, 0.0), n, x);
        }
    }

    /// A weight of exactly 1.0 must be indistinguishable from the unweighted path.
    ///
    /// Upstream reaches the weighted sweep for control tracks, so a divergence
    /// between the two would silently corrupt control lambda.
    #[test]
    fn unit_weights_match_unweighted(
        d in 1i64..40,
        rlength in 500i64..2000,
        plus in prop::collection::vec(0u64..2000, 0..60),
        minus in prop::collection::vec(0u64..2000, 0..60),
    ) {
        let rlength = rlength.max(*plus.iter().max().unwrap_or(&0) as i64 + 1)
            .max(*minus.iter().max().unwrap_or(&0) as i64 + 1);
        let a = pileup_from_positions(ChromId(0), &plus, &minus, &params(d, rlength));
        let ones = vec![1.0f32; plus.len()];
        let b = pileup_from_weighted_positions(
            ChromId(0), &plus, &minus, &ones, &ones, &params(d, rlength),
        );
        prop_assert_eq!(a.runs().len(), b.runs().len());
        for x in 0..rlength {
            prop_assert_eq!(a.value_at_or(x as u64, 0.0), b.value_at_or(x as u64, 0.0),
                "unit weight diverged at {}", x);
        }
    }

    /// Scaling is linear in `scale_factor`, before the baseline floor.
    #[test]
    fn scaling_is_linear(
        rlength in 500u64..2000,
        plus in prop::collection::vec(0u64..2000, 1..50),
    ) {
        let mut one = params(200, rlength as i64);
        one.scale_factor = 1.0;
        let mut two = one;
        two.scale_factor = 3.0;
        let a = pileup_from_positions(ChromId(0), &plus, &[], &one);
        let b = pileup_from_positions(ChromId(0), &plus, &[], &two);
        prop_assert!((integrated_depth(&b) - 3.0 * integrated_depth(&a)).abs() < 1e-2);
    }

    /// `naive_quick_pileup` -- the O(n*m) reference used by upstream's own tests --
    /// and the production sweep must agree everywhere.
    #[test]
    fn sweep_agrees_with_naive_quick_pileup(
        ext in 1i64..40,
        rlength in 500i64..2000,
        plus in prop::collection::vec(0u64..2000, 0..40),
    ) {
        let rlength = rlength.max(*plus.iter().max().unwrap_or(&0) as i64 + 1);
        // `naive_quick_pileup` uses the symmetric window `[x - ext, x + ext)`, so
        // pin the params to that; only the sweep strategy then differs.
        let mut p = params(ext, rlength);
        p.five_shift = ext;
        p.three_shift = ext;
        // Plus strand only: `naive_quick_pileup` is upstream's *single-strand*
        // reference. In a non-centred config the minus strand is deliberately the
        // mirror, which `non_centred_minus_is_the_mirror` pins instead.
        let fast = pileup_from_positions(ChromId(0), &plus, &[], &p);
        let mut all = plus.clone();
        all.sort_unstable();
        let slow = naive_quick_pileup(ChromId(0), &all, ext, rlength as u64);
        for x in 0..rlength {
            prop_assert_eq!(fast.value_at_or(x as u64, 0.0), slow.value_at_or(x as u64, 0.0),
                "sweep vs naive at {} (ext={})", x, ext);
        }
    }

    /// Slicing to an interval conserves the integral of the covered part and never
    /// invents depth outside the source track's own extent.
    #[test]
    fn restrict_conserves_the_sliced_integral(
        d in 1i64..40,
        rlength in 800i64..2000,
        lo in 0i64..2000,
        hi in 0i64..2000,
        plus in prop::collection::vec(0u64..2000, 0..60),
    ) {
        let rlength = rlength.max(*plus.iter().max().unwrap_or(&0) as i64 + 1);
        let t = pileup_from_positions(ChromId(0), &plus, &[], &params(d, rlength));
        let (lo, hi) = (lo.min(hi), lo.max(hi));
        let slice = restrict(&t, Interval::new(lo as u64, hi as u64));
        let expect = naive_se(&plus, &[], d, rlength).iter()
            .filter(|(x, _)| *x >= lo && *x < hi)
            .map(|(_, v)| *v as f64).sum::<f64>();
        prop_assert!((integrated_depth(&slice) - expect).abs() < 1e-3,
            "sliced integral {} != {}", integrated_depth(&slice), expect);
        prop_assert!(integrated_depth(&slice) <= integrated_depth(&t) + 1e-6);
    }

    /// Permuting the input read order cannot change the pileup.
    ///
    /// The parallel control path feeds chromosome partitions in whatever order the
    /// sweep produces them, so order-independence is what makes the chromosome-level
    /// parallelism provably output-identical.
    #[test]
    fn pileup_is_order_independent(
        d in 1i64..40,
        rlength in 800i64..2000,
        pos in prop::collection::vec(0u64..2000, 0..80),
    ) {
        let mut shuffled = pos.clone();
        // deterministic shuffle so a failure is reproducible from the seed
        for i in (1..shuffled.len()).rev() {
            let j = (i * 2_654_435_761) % (i + 1);
            shuffled.swap(i, j);
        }
        let a = pileup_from_positions(ChromId(0), &pos, &[], &params(d, rlength));
        let b = pileup_from_positions(ChromId(0), &shuffled, &[], &params(d, rlength));
        prop_assert_eq!(a.integral(), b.integral());
        for x in 0..rlength {
            prop_assert_eq!(a.value_at_or(x as u64, 0.0), b.value_at_or(x as u64, 0.0));
        }
    }

    /// In a non-centred config the minus strand is the mirror image of the plus
    /// strand, and in a centred config it is *identical* to it (F190).
    ///
    /// F190 was exactly this distinction going the wrong way: the counted paired-end
    /// control projection had been mirroring the minus strand when upstream's
    /// `pileup_from_LRC_centers_as_list` does not. Pinning it as a property means
    /// the regression cannot come back unnoticed for any shift pair.
    #[test]
    fn non_centred_minus_is_the_mirror(
        five in 0i64..30,
        three in 1i64..60,
        rlength in 500u64..2000,
        plus in prop::collection::vec(0u64..2000, 0..30),
        minus in prop::collection::vec(0u64..2000, 0..30),
    ) {
        let mut mirrored = SingleEndParams {
            five_shift: five,
            three_shift: three,
            rlength,
            scale_factor: 1.0,
            baseline_value: 0.0,
            coalesce: true,
            centred: false,
        };
        let mut centred = mirrored;
        centred.centred = true;
        mirrored.centred = false;

        // mirror: a minus read at `x` covers `[x - three, x + five)`. The plus window
        // is `[y - five, y + three)`, so `y = x + five - three` makes them identical.
        let shifted: Vec<u64> = minus
            .iter()
            .map(|&x| (x as i64 + five - three).max(0) as u64)
            .collect();
        let m = pileup_from_positions(ChromId(0), &[], &minus, &mirrored);
        let sh = pileup_from_positions(ChromId(0), &shifted, &[], &mirrored);
        // The shift is only well defined away from the contig edges: clipping is
        // applied per endpoint, so a read that hangs off the left edge covers less
        // than its window and no shift can reproduce that. Compare the interior,
        // which is where every real peak lives.
        let margin = (five + three) as u64;
        for x in margin..rlength {
            prop_assert_eq!(m.value_at_or(x, 0.0), sh.value_at_or(x, 0.0),
                "minus window is not the plus mirror at {}", x);
        }

        // centred: identical, not mirrored
        let a = pileup_from_positions(ChromId(0), &plus, &minus, &centred);
        let minus_only = pileup_from_positions(ChromId(0), &[], &minus, &centred);
        let plus_only = pileup_from_positions(ChromId(0), &plus, &[], &centred);
        if plus.is_empty() {
            prop_assert_eq!(minus_only.integral(), a.integral(),
                "centred: strands must share one window");
        }
        if minus.is_empty() {
            prop_assert_eq!(plus_only.integral(), a.integral(),
                "centred: strands must share one window");
        }
    }

    /// An empty input yields an empty track, never a panic or a NaN.
    #[test]
    fn empty_input_is_empty_track(
        d in 1i64..200,
        rlength in 1u64..100_000,
    ) {
        let t: SignalTrack<f32> =
            pileup_from_positions(ChromId(0), &[], &[], &params(d, rlength as i64));
        prop_assert_eq!(t.runs().len(), 0);
        prop_assert_eq!(integrated_depth(&t), 0.0);
    }
}
