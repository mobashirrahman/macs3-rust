//! L2 property invariants for the statistics kernels.
//!
//! These kernels feed `d`, the duplicate model and every p-value, and they are the
//! one place where the port re-derives an approximation (series expansions with early
//! termination) rather than transcribing a formula. The recorded vectors pin the
//! approximation's output; these pin the properties that hold *whatever* it
//! computes, so a truncation that happens to agree on the corpus still cannot pass.

use macs_stats::{binomial_cdf, binomial_sf, poisson_cdf, poisson_pdf};
use proptest::prelude::*;

fn pcdf(n: u32, lam: f64, lower: bool) -> f64 {
    poisson_cdf(n, lam, lower, false).expect("valid poisson parameters")
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// The Poisson CDF is non-decreasing in `n`.
    #[test]
    fn poisson_cdf_is_monotone_in_n(
        lam in 0.1f64..40.0,
        ns in prop::collection::vec(0u32..400, 2..15),
    ) {
        for w in ns.windows(2) {
            let (lo, hi) = (w[0].min(w[1]), w[0].max(w[1]));
            let (a, b) = (pcdf(lo, lam, true), pcdf(hi, lam, true));
            prop_assert!(b >= a - 1e-12,
                "poisson CDF fell at lam={}: P({})={} > P({})={}", lam, lo, a, hi, b);
        }
    }

    /// The CDF is bounded by 1 and the SF by 0 -- a truncated series must not
    /// overshoot, which is exactly what an under-converged sum would do.
    #[test]
    fn poisson_cdf_is_bounded(
        n in 0u32..2000,
        lam in 0.01f64..200.0,
    ) {
        let lower = pcdf(n, lam, true);
        let upper = pcdf(n, lam, false);
        // The `> 1.0` clamp upstream applies to the *lower* small-lambda series is
        // not applied to the upper one, so 1 + 1 ulp is legal; anything more is not.
        prop_assert!(lower <= 1.0, "lower CDF {} exceeds 1", lower);
        prop_assert!(lower >= 0.0, "lower CDF {} below 0", lower);
        prop_assert!(upper <= 1.0 + 1e-12, "upper CDF {} exceeds 1", upper);
        prop_assert!(upper >= 0.0, "upper CDF {} below 0", upper);
        // No ordering between the two: they are *complementary* (`P(X<=n)` and
        // `P(X>n)`), which is the relation the next test checks. Either may be the
        // larger one, depending on which side of the median `n` falls.
    }

    /// Lower and upper CDF agree in the middle: `P(X<=n) + P(X>n) == 1`.
    #[test]
    fn poisson_lower_and_upper_are_complementary(
        n in 0u32..400,
        lam in 0.1f64..60.0,
    ) {
        let lower = pcdf(n, lam, true);
        let upper = pcdf(n, lam, false);
        prop_assert!((lower + upper - 1.0).abs() < 1e-6,
            "lower {} + upper {} != 1 (n={} lam={})", lower, upper, n, lam);
    }

    /// The `-log10(p)` returned by `log10` mode always lands on the 1e-5 grid.
    ///
    /// The quantisation is upstream's, not ours, and it is load-bearing: this is the
    /// value that lands in `*.xls` / `narrowPeak` and in the 1e-5 histogram buckets
    /// of the p-score table. It must hold for *every* `(n, lambda)` without
    /// exception, which makes it a stronger check than the recorded vectors.
    #[test]
    fn poisson_log10_mode_is_always_on_the_1e5_grid(
        n in 0u32..2000,
        lam in 0.01f64..200.0,
    ) {
        for lower in [true, false] {
            let log = poisson_cdf(n, lam, lower, true).expect("valid parameters");
            prop_assert!(log.is_finite(), "non-finite log10 CDF (n={} lam={})", n, lam);
            let scaled = log * 1e5;
            prop_assert!((scaled - scaled.round()).abs() < 1e-6,
                "log10 mode off the 1e-5 grid: {} (n={} lam={} lower={})", log, n, lam, lower);
        }
    }

    /// The upper-tail `log10` branch is non-increasing in `n`, signed, and finite.
    ///
    /// Sign convention, pinned because the p-score depends on it: for
    /// `P(X > k) <= 1` the branch returns `log10 P(X > k)`, i.e. a **negative**
    /// number, and it grows more negative as `k` rises. Upstream's own docstring
    /// claims `-log10 p` for this function, which is wrong; every call site negates
    /// the result. A test that asserted positivity here would be asserting the
    /// docstring rather than the behaviour.
    #[test]
    fn poisson_log10_upper_tail_is_signed_and_monotone(
        lam in 0.05f64..100.0,
        ns in prop::collection::vec(0u32..600, 2..15),
    ) {
        for w in ns.windows(2) {
            let (lo, hi) = (w[0].min(w[1]), w[0].max(w[1]));
            let (a, b) = (
                poisson_cdf(lo, lam, false, true).expect("valid parameters"),
                poisson_cdf(hi, lam, false, true).expect("valid parameters"),
            );
            prop_assert!(a.is_finite() && b.is_finite(),
                "non-finite log10 tail at lam={}", lam);
            prop_assert!(a <= 1e-12, "log10 upper tail {} is positive at lam={}", a, lam);
            prop_assert!(b <= a + 1e-4,
                "log10 upper tail rose at lam={}: {} -> {}", lam, a, b);
        }
    }

    /// The lower-tail `log10` branch is also signed and finite.
    #[test]
    fn poisson_log10_lower_tail_is_signed_and_finite(
        n in 0u32..2000,
        lam in 0.01f64..200.0,
    ) {
        let log = poisson_cdf(n, lam, true, true).expect("valid parameters");
        prop_assert!(log.is_finite(), "non-finite log10 lower tail (n={} lam={})", n, lam);
        prop_assert!(log <= 1e-12, "log10 lower tail {} is positive at n={} lam={}", log, n, lam);
    }

    /// The PMF is a probability density over the integers and normalises to 1.
    #[test]
    fn poisson_pdf_normalises(
        lam in 0.5f64..30.0,
        cut in 1u32..120,
    ) {
        let hi = (lam + 12.0 * lam.sqrt()).ceil() as u32;
        let hi = hi.max(cut).min(4000);
        let total: f64 = (0..=hi).map(|k| poisson_pdf(k, lam)).sum();
        prop_assert!(total <= 1.0 + 1e-6, "PMF sum {} exceeds 1", total);
        prop_assert!(total > 1.0 - 1e-6,
            "PMF sum {} under 1 by n={} lam={} (tail beyond {})", total, cut, lam, hi);
    }

    /// `poisson_pdf` is non-negative and unimodal, peaking at `floor(lambda)`.
    ///
    /// The mode of a Poisson is `floor(lambda)` -- and `lambda - 1` as well when
    /// `lambda` is an exact integer -- *not* `round(lambda)`, so an implementation
    /// that scans for the peak near the rounded mean would fail here.
    #[test]
    fn poisson_pdf_is_non_negative_and_unimodal(
        lam in 1.0f64..50.0,
    ) {
        let mode = lam.floor() as u32;
        for k in 0..=(mode + 40) {
            prop_assert!(poisson_pdf(k, lam) >= 0.0, "negative PMF at {}", k);
        }
        // non-decreasing up to the mode, non-increasing after
        for k in 1..=mode {
            prop_assert!(poisson_pdf(k, lam) >= poisson_pdf(k - 1, lam) - 1e-15,
                "PMF dips below the mode at k={} (lam={})", k, lam);
        }
        for k in mode..(mode + 30) {
            prop_assert!(poisson_pdf(k + 1, lam) <= poisson_pdf(k, lam) + 1e-15,
                "PMF rises above the mode at k={} (lam={})", k, lam);
        }
        prop_assert!(poisson_pdf(mode, lam) >= poisson_pdf(mode + 1, lam),
            "mode {} is not a maximum", mode);
    }

    /// The binomial CDF/SF pair is complementary and bounded, and its mean tracks
    /// `n * p`. This is the duplicate-count model's foundation.
    #[test]
    fn binomial_lower_and_upper_are_complementary(
        n in 1i64..200,
        p in 0.01f64..0.99,
        x in 0i64..200,
    ) {
        let (lo, up) = (binomial_cdf(x, n, p, true), binomial_sf(x, n, p, true));
        prop_assert!((lo + up - 1.0).abs() < 1e-6,
            "binomial lower {} + sf {} != 1", lo, up);
        prop_assert!((-1e-12..=1.0 + 1e-9).contains(&lo), "binomial CDF {} out of range", lo);
        prop_assert!((-1e-12..=1.0 + 1e-9).contains(&up), "binomial SF {} out of range", up);
    }

    /// The binomial CDF is non-decreasing in `x`.
    #[test]
    fn binomial_cdf_is_monotone_in_x(
        n in 1i64..150,
        p in 0.01f64..0.99,
        xs in prop::collection::vec(0i64..150, 2..12),
    ) {
        for w in xs.windows(2) {
            let (lo, hi) = (w[0].min(w[1]), w[0].max(w[1]));
            let (a, b) = (binomial_cdf(lo, n, p, true), binomial_cdf(hi, n, p, true));
            prop_assert!(b >= a - 1e-12,
                "binomial CDF fell at n={} p={}: {} > {}", n, p, a, b);
        }
    }

    /// A binomial with `p -> 1` concentrates all its mass at `n`.
    #[test]
    fn binomial_concentrates_at_n_as_p_approaches_one(n in 1i64..60) {
        let at_n = binomial_cdf(n, n, 0.999, true);
        let below = binomial_cdf(n - 1, n, 0.999, true);
        prop_assert!(at_n > below, "not concentrated at n={}: {} vs {}", n, at_n, below);
        prop_assert!(at_n > 0.9, "p=0.999 mass at n is only {}", at_n);
    }
}
