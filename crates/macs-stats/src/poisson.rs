//! Poisson distribution functions.
//!
//! Port of `MACS3/Signal/Prob.py` lines 250-632. Upstream implements its own
//! Poisson tails rather than calling SciPy, with separate small-`lambda` and
//! large-`lambda` code paths and a log-space path that quantises the result to
//! five decimal digits. MACS peak calls go through
//! [`poisson_cdf`], so all of that structure is load-bearing for parity.

use crate::util::{logspace_add, py_round};
use crate::StatsError;

/// Split step used by the large-lambda paths (`Prob.py: LSTEP = 200`).
const LSTEP: i64 = 200;
/// `exp(200)`, exactly: `0x1.73f60ea79f5b9p+288`.
///
/// Written as a hex float literal so the constant is the same bit pattern as
/// upstream's `EXPTHRES: float64_t = exp(LSTEP)` without depending on a
/// const-evaluable `exp`.
const EXPTHRES: f64 = f64::from_bits(0x51F7_3F60_EA79_F5B9);
/// `exp(-200)`, exactly: `0x1.6061812054cfap-289`.
const EXPSTEP: f64 = f64::from_bits(0x2DE6_0618_1205_4CFA);

/// Poisson CDF evaluator.
///
/// Port of `poisson_cdf(n, lam, lower, log10)`.
///
/// * `lower == false` (default, and what MACS uses) is the **upper tail**:
///   `P(X > n)`.
/// * `log10 == true` switches to the log-space paths, which return
///   `-log10(P)` and round to five decimals.
///
/// # Errors
/// Returns [`StatsError::NonPositiveLambda`] for `lam <= 0`, matching upstream's
/// `assert lam > 0.0`.
pub fn poisson_cdf(n: u32, lam: f64, lower: bool, log10: bool) -> Result<f64, StatsError> {
    // `!(x > 0.0)` rather than `x <= 0.0` so that NaN is rejected too, matching
    // upstream's `assert lam > 0.0`
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    if !(lam > 0.0) {
        return Err(StatsError::NonPositiveLambda(lam));
    }

    if log10 {
        return Ok(if lower {
            log10_poisson_cdf_p_large_lambda(n, lam)
        } else {
            log10_poisson_cdf_q_large_lambda(n, lam)
        });
    }

    if lower {
        if lam > 700.0 {
            Ok(poisson_cdf_large_lambda(n, lam))
        } else {
            Ok(poisson_cdf_small_lambda(n, lam))
        }
    } else if lam > 700.0 {
        Ok(poisson_cdf_q_large_lambda(n, lam))
    } else {
        Ok(poisson_cdf_q(n, lam))
    }
}

/// Lower-tail Poisson CDF for small lambda (`__poisson_cdf`).
///
/// Note the upstream clamp `if cdf > 1.0: return 1.0`; it is reproduced.
fn poisson_cdf_small_lambda(k: u32, a: f64) -> f64 {
    let mut nextcdf = (-a).exp();
    let mut cdf = nextcdf;

    for i in 1..=k {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / f64::from(i);
        cdf += nextcdf;
    }
    if cdf > 1.0 {
        1.0
    } else {
        cdf
    }
}

/// Lower-tail Poisson CDF for large lambda (`__poisson_cdf_large_lambda`).
///
/// The running sum is periodically rescaled by `exp(-200)` so that it cannot
/// overflow while iterating up to `k`. The number of rescalings is precomputed
/// from `lambda / 200`, so the final `for i in range(num_parts)` loop applies
/// the *remaining* scalings in one go. If `num_parts` is exhausted early the
/// loop is empty (upstream relies on C semantics for a negative range).
fn poisson_cdf_large_lambda(k: u32, a: f64) -> f64 {
    let mut num_parts: i64 = (a / LSTEP as f64) as i64;
    let mut lastexp = (-(a % LSTEP as f64)).exp();
    let mut nextcdf = EXPSTEP;
    let mut cdf = 0.0_f64;

    num_parts -= 1;

    for i in 1..=k {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / f64::from(i);
        cdf += nextcdf;
        if nextcdf > EXPTHRES || cdf > EXPTHRES {
            if num_parts >= 1 {
                cdf *= EXPSTEP;
                nextcdf *= EXPSTEP;
                num_parts -= 1;
            } else {
                cdf *= lastexp;
                lastexp = 1.0;
            }
        }
    }

    for _ in 0..num_parts {
        cdf *= EXPSTEP;
    }
    cdf *= lastexp;
    cdf
}

/// Upper-tail Poisson CDF for small lambda (`__poisson_cdf_Q`).
///
/// The tail `P(X > k) = sum_{i>k} pmf(i)` is summed **directly** by walking the
/// recurrence forward from `pmf(k)` until it underflows to zero. It is *not*
/// `1 - P(X <= k)`: the two disagree in the last ULPs, and that difference is
/// visible in a q-score.
fn poisson_cdf_q(k: u32, a: f64) -> f64 {
    let mut nextcdf = (-a).exp();
    for i in 1..=k {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / f64::from(i);
    }
    // `nextcdf` is now pmf(k); walk the tail forward until it underflows
    let mut cdf = 0.0_f64;
    let mut i = f64::from(k) + 1.0;
    while nextcdf > 0.0 {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / i;
        cdf += nextcdf;
        i += 1.0;
    }
    cdf
}

/// Upper-tail Poisson CDF for large lambda (`__poisson_cdf_Q_large_lambda`).
///
/// Also a direct tail sum, but with periodic `exp(-200)` rescaling so the
/// running sum cannot overflow while iterating up to `k`. Upstream's unreachable
/// `raise Exception("Unexpected error")` on the first loop is preserved as a
/// silent rescale, because `num_parts` is at least 2 for every `a > 700` and the
/// branch is dead.
fn poisson_cdf_q_large_lambda(k: u32, a: f64) -> f64 {
    let mut num_parts: i64 = (a / LSTEP as f64) as i64;
    let mut lastexp = (-(a % LSTEP as f64)).exp();
    let mut nextcdf = EXPSTEP;
    let mut cdf = 0.0_f64;

    num_parts -= 1;

    for i in 1..=k {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / f64::from(i);
        if nextcdf > EXPTHRES && num_parts >= 1 {
            nextcdf *= EXPSTEP;
            num_parts -= 1;
        }
        // upstream raises "Unexpected error" in the `else`; `num_parts` is at
        // least 2 for every `a > 700`, so the branch is dead and we do nothing
    }

    let mut i = f64::from(k) + 1.0;
    while nextcdf > 0.0 {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / i;
        cdf += nextcdf;
        i += 1.0;
        if nextcdf > EXPTHRES || cdf > EXPTHRES {
            if num_parts >= 1 {
                cdf *= EXPSTEP;
                nextcdf *= EXPSTEP;
                num_parts -= 1;
            } else {
                cdf *= lastexp;
                lastexp = 1.0;
            }
        }
    }

    for _ in 0..num_parts {
        cdf *= EXPSTEP;
    }
    cdf *= lastexp;
    cdf
}

/// `log10` of the lower-tail Poisson probability
/// (`log10_poisson_cdf_P_large_lambda`).
///
/// Returns `log10(P(X <= k))`, rounded to five decimals. See
/// [`log10_poisson_cdf_q_large_lambda`] for the upstream sign convention.
fn log10_poisson_cdf_p_large_lambda(k: u32, lbd: f64) -> f64 {
    let ln_lbd = lbd.ln();
    let mut m = i64::from(k);
    let mut sum_ln_m = 0.0_f64;
    for i in 1..=m {
        sum_ln_m += (i as f64).ln();
    }
    let mut logx = m as f64 * ln_lbd - sum_ln_m;
    let mut residue = logx;

    while m > 1 {
        m -= 1;
        let logy = logx - ln_lbd + (m as f64).ln();
        let pre_residue = residue;
        residue = logspace_add(pre_residue, logy);
        if (pre_residue - residue).abs() < 1e-10 {
            break;
        }
        logx = logy;
    }

    py_round((residue - lbd) / std::f64::consts::LN_10, 5)
}

/// `log10` of the upper-tail Poisson probability
/// (`log10_poisson_cdf_Q_large_lambda`).
///
/// **Sign warning.** Upstream's docstring claims this returns `-log10 p`, but
/// the code computes `log10 p`, which is `<= 0`. Every call site negates the
/// result ([`MACS3/Signal/ScoreTrack.py:83`] and
/// [`MACS3/Signal/CallPeakUnit.py:71`]: `get_pscore = -1 * poisson_cdf(k, lam,
/// False, True)`), so the observable p-score is the positive `-log10 p`. We
/// reproduce the *code*, not the docstring, and this function therefore returns
/// `log10(P(X > k))`.
///
/// The result is rounded to five decimals, so every p-score MACS produces is a
/// multiple of 1e-5. The convergence threshold is `1e-5` here versus `1e-10` in
/// the lower-tail version; that asymmetry is upstream's and is preserved.
///
/// F191: `1e-5` is confirmed correct. The *compiled* extension agrees with this
/// transcription on 1285 of a 1410-point oracle sweep spanning
/// `-log10(p)` from 0 to 1000, and the 125 disagreements are each exactly one
/// unit in the last (5th) decimal -- the rounding boundary, not a different
/// result. See `docs/upstream-findings.md` for what has been ruled out.
fn log10_poisson_cdf_q_large_lambda(k: u32, lbd: f64) -> f64 {
    let ln_lbd = lbd.ln();
    let mut m = i64::from(k) + 1;
    let mut sum_ln_m = 0.0_f64;
    for i in 1..=m {
        sum_ln_m += (i as f64).ln();
    }
    let mut logx = m as f64 * ln_lbd - sum_ln_m;
    let mut residue = logx;

    loop {
        m += 1;
        let logy = logx + ln_lbd - (m as f64).ln();
        let pre_residue = residue;
        residue = logspace_add(pre_residue, logy);
        if (pre_residue - residue).abs() < 1e-5 {
            break;
        }
        logx = logy;
    }

    py_round((residue - lbd) / 10.0_f64.ln(), 5)
}

/// Inverse Poisson CDF (`poisson_cdf_inv`).
///
/// Returns the smallest `i` in `1..=maximum` with
/// `P(X <= i-1) <= cdf <= P(X <= i)`; `0` when `cdf == 0`; `maximum` when no
/// such `i` exists.
///
/// # Errors
/// `lambda` must be `< 740` and `cdf` must lie in `[0, 1]`.
pub fn poisson_cdf_inv(cdf: f64, lam: f64, maximum: i32) -> Result<i32, StatsError> {
    // NaN must also be rejected, as upstream's `assert lam < 740` does
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    if !(lam < 740.0) {
        return Err(StatsError::LambdaTooLarge(lam));
    }
    if !(0.0..=1.0).contains(&cdf) {
        return Err(StatsError::CdfOutOfRange(cdf));
    }
    if cdf == 0.0 {
        return Ok(0);
    }

    let mut newval = (-lam).exp();
    let mut sum2 = newval;

    for i in 1..=maximum {
        let sumold = sum2;
        newval = newval * lam / f64::from(i);
        sum2 += newval;
        if sumold <= cdf && cdf <= sum2 {
            return Ok(i);
        }
    }

    Ok(maximum)
}

/// Upstream's `poisson_cdf_Q_inv`.
///
/// **Upstream bug reproduced deliberately:** the body of
/// `poisson_cdf_Q_inv` (`Prob.py:588-618`) is byte-for-byte the body of
/// `poisson_cdf_inv`; it does not invert the upper tail at all. MACS3 never
/// calls it, so nothing observable depends on the bug, but we expose it under
/// the upstream name so a future differential test can prove the two agree.
pub fn poisson_cdf_q_inv(cdf: f64, lam: f64, maximum: i32) -> Result<i32, StatsError> {
    poisson_cdf_inv(cdf, lam, maximum)
}

/// Poisson PDF (`poisson_pdf`).
///
/// Upstream evaluates `exp(-a) * a**k / k!` directly, which underflows to `0`
/// for moderate `k` and moderate `a` (e.g. `k = 1000, a = 10`). We reproduce
/// that, because MACS3's p-scores are built from [`poisson_cdf`] and a
/// silently "fixed" PDF would be a different program.
pub fn poisson_pdf(k: u32, a: f64) -> f64 {
    if a <= 0.0 {
        return 0.0;
    }
    // `powf`, not `powi`: upstream calls C `pow(a, k)` from libm, which is
    // correctly rounded, whereas `f64::powi` is a repeated-multiplication loop
    // and is 1-2 ULP off for moderate k.
    let r = (-a).exp() * a.powf(f64::from(k)) / factorial(k);
    if r.is_finite() {
        r
    } else {
        // `a**k` and `k!` both overflow for large k, giving inf/inf = NaN.
        // Upstream raises `OverflowError` here; we return 0 so the function is
        // total. Use `poisson_pdf_log` when a stable value is needed.
        0.0
    }
}

/// Poisson PDF computed stably in log space, then exponentiated.
///
/// This is *not* used by MACS-compatible code paths. It exists so that the
/// ported functions can be validated against an independent reference, and for
/// [`poisson_pdf`]'s own unit tests.
pub fn poisson_pdf_log(k: u32, a: f64) -> f64 {
    if a <= 0.0 {
        return 0.0;
    }
    let logp = -a + f64::from(k) * a.ln() - ln_factorial(k);
    logp.exp()
}

/// `n!` as a float (`factorial`).
///
/// Returns `inf` once `n!` overflows `f64`, which is what upstream's `double`
/// accumulator does. Exposed because the oracle vectors it.
pub fn factorial(n: u32) -> f64 {
    let mut fact = 1.0_f64;
    for i in 2..=n {
        fact *= f64::from(i);
    }
    fact
}

/// `log(n!)` via `lgamma`, which does not overflow.
fn ln_factorial(n: u32) -> f64 {
    if n < 2 {
        0.0
    } else {
        ln_gamma_ref(f64::from(n) + 1.0)
    }
}

/// `ln Gamma(x)` for `x > 0` (Lanczos approximation, g=7, n=9).
///
/// Exposed so that the sibling binomial reference implementation in this crate
/// (and the unit tests) can build an *independent* PMF to check the ported
/// ratio-product code against. It is never used by a MACS-compatible path.
pub fn ln_gamma_ref(x: f64) -> f64 {
    const G: f64 = 7.0;
    #[allow(clippy::excessive_precision)]
    const C: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_13,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // reflection formula
        (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma_ref(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = C[0];
        let t = x + G + 0.5;
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Independent reference: exact Poisson tail by direct summation in
    /// arbitrary precision-ish form, for small parameters only.
    fn poisson_upper_tail_ref(k: u32, lam: f64) -> f64 {
        // P(X > k) = 1 - sum_{i=0}^{k} exp(-lam) lam^i / i!
        let mut term = (-lam).exp();
        let mut cdf = term;
        for i in 1..=k {
            term *= lam / f64::from(i);
            cdf += term;
        }
        1.0 - cdf
    }

    fn poisson_lower_tail_ref(k: u32, lam: f64) -> f64 {
        let mut term = (-lam).exp();
        let mut cdf = term;
        for i in 1..=k {
            term *= lam / f64::from(i);
            cdf += term;
        }
        cdf
    }

    #[test]
    fn lower_tail_matches_reference() {
        for lam in [0.001, 0.5, 1.0, 3.7, 50.0, 300.0, 699.0] {
            for k in [0u32, 1, 5, 37, 200] {
                let got = poisson_cdf(k, lam, true, false).unwrap();
                let want = poisson_lower_tail_ref(k, lam);
                assert!(
                    (got - want).abs() < 1e-12,
                    "lower lam={lam} k={k} got={got} want={want}"
                );
            }
        }
    }

    #[test]
    fn upper_tail_small_lambda_matches_reference() {
        for lam in [0.001, 0.5, 1.0, 3.7, 50.0, 300.0, 699.0] {
            for k in [0u32, 1, 5, 37, 200] {
                let got = poisson_cdf(k, lam, false, false).unwrap();
                let want = poisson_upper_tail_ref(k, lam);
                // The reference and the port both compute 1 - cdf, but in a
                // different order, so 1e-12 is the right bound.
                assert!(
                    (got - want).abs() < 1e-12,
                    "upper lam={lam} k={k} got={got} want={want}"
                );
            }
        }
    }

    #[test]
    fn cdf_is_monotone_in_k() {
        let lam = 12.5;
        let mut prev = poisson_cdf(0, lam, true, false).unwrap();
        for k in 1..200 {
            let cur = poisson_cdf(k, lam, true, false).unwrap();
            assert!(cur >= prev, "lower tail decreased at k={k}");
            prev = cur;
        }
        let mut prev = poisson_cdf(0, lam, false, false).unwrap();
        for k in 1..200 {
            let cur = poisson_cdf(k, lam, false, false).unwrap();
            assert!(cur <= prev, "upper tail increased at k={k}");
            prev = cur;
        }
    }

    #[test]
    fn lower_plus_upper_is_one() {
        for lam in [0.1, 5.0, 100.0, 700.0] {
            for k in [0u32, 3, 99] {
                let lo = poisson_cdf(k, lam, true, false).unwrap();
                let up = poisson_cdf(k, lam, false, false).unwrap();
                assert!(
                    (lo + up - 1.0).abs() < 1e-9,
                    "lam={lam} k={k} lo={lo} up={up}"
                );
            }
        }
    }

    /// The p-score MACS actually uses is `-poisson_cdf(k, lam, false, true)`.
    #[test]
    fn pscore_is_minus_the_log10_path_and_is_quantised() {
        for (k, lam) in [(0u32, 1.0f64), (5, 2.0), (10, 3.0), (30, 10.0), (100, 50.0)] {
            let log10p = poisson_cdf(k, lam, false, true).unwrap();
            assert_eq!(log10p, py_round(log10p, 5), "k={k} lam={lam}");
            let pscore = -log10p;
            assert!(pscore >= 0.0, "pscore must be non-negative: {pscore}");
        }
        // P(X > 30 | lambda = 10) = 7.983795e-08, so -log10 = 7.09779.
        // Checked against an independent direct tail summation, upstream's
        // log-space accumulation is accurate to exactly its 5-decimal rounding
        // (docs/upstream-findings.md F1).
        let (k, lam) = (30u32, 10.0f64);
        let pscore = -poisson_cdf(k, lam, false, true).unwrap();
        let true_pscore = -poisson_cdf_q(k, lam).log10();
        assert!((pscore - 7.097_79).abs() < 1e-5, "got {pscore}");
        assert!(
            (pscore - true_pscore).abs() <= 2e-5,
            "pscore={pscore} direct_sum={true_pscore}"
        );
    }

    #[test]
    fn non_positive_lambda_is_an_error_not_a_panic() {
        assert!(matches!(
            poisson_cdf(1, 0.0, false, false),
            Err(StatsError::NonPositiveLambda(_))
        ));
        assert!(matches!(
            poisson_cdf(1, -1.0, false, false),
            Err(StatsError::NonPositiveLambda(_))
        ));
    }

    /// Upstream's `poisson_cdf_inv` looks for the first `i` with
    /// `P(X <= i-1) <= cdf <= P(X <= i)`. Because `sumold` starts at
    /// `P(X = 0) = exp(-lam)`, a `cdf` below that is **unreachable** and the
    /// function returns `maximum`. This is upstream behaviour, not a bug we
    /// introduced, and it is exercised by `--keep-dup` only for large `N`.
    #[test]
    fn inverse_cdf_returns_maximum_for_unreachable_cdf() {
        for lam in [1.0f64, 5.0, 20.0] {
            let cdf: f64 = (-lam).exp() * 0.5;
            assert_eq!(poisson_cdf_inv(cdf, lam, 1000).unwrap(), 1000, "lam={lam}");
        }
    }

    #[test]
    fn inverse_cdf_brackets_the_cdf() {
        for lam in [1.0, 10.0, 100.0, 700.0] {
            for cdf in [0.5, 0.9, 0.99] {
                let i = poisson_cdf_inv(cdf, lam, 1000).unwrap();
                let lo = if i == 0 {
                    0.0
                } else {
                    poisson_cdf((i - 1) as u32, lam, true, false).unwrap()
                };
                let hi = poisson_cdf(i as u32, lam, true, false).unwrap();
                assert!(
                    lo <= cdf && cdf <= hi + 1e-12,
                    "lam={lam} cdf={cdf} i={i} lo={lo} hi={hi}"
                );
            }
        }
        assert_eq!(poisson_cdf_inv(0.0, 1.0, 1000).unwrap(), 0);
        assert!(poisson_cdf_inv(0.5, 800.0, 1000).is_err());
        assert!(poisson_cdf_inv(-0.1, 1.0, 1000).is_err());
    }

    #[test]
    fn pdf_overflows_like_upstream_but_log_version_does_not() {
        // upstream raises OverflowError; we return 0 to stay total
        assert_eq!(poisson_pdf(1000, 10.0), 0.0);
        assert!(!poisson_pdf(200, 30.0).is_nan());
        // the stable path stays accurate deep into the far tail, where the
        // direct one has already lost everything to inf/inf
        let deep = poisson_pdf_log(100, 10.0);
        assert!(deep > 0.0 && deep < 1e-60, "got {deep}");
        // and they agree where neither overflows nor underflows
        for (k, a) in [(0u32, 1.0), (3, 3.0), (10, 10.0), (40, 50.0)] {
            let a1 = poisson_pdf(k, a);
            let a2 = poisson_pdf_log(k, a);
            assert!((a1 - a2).abs() < 1e-15, "k={k} a={a} {a1} vs {a2}");
        }
    }

    #[test]
    fn large_lambda_paths_are_finite() {
        for lam in [700.5, 1000.0, 5000.0, 50_000.0] {
            for k in [0u32, 10, 100] {
                let lo = poisson_cdf(k, lam, true, false).unwrap();
                let up = poisson_cdf(k, lam, false, false).unwrap();
                assert!(lo.is_finite() && up.is_finite(), "lam={lam} k={k}");
                assert!((0.0..=1.0).contains(&lo), "lam={lam} k={k} lo={lo}");
                // Upstream's rescaled accumulation can land a few ULP above 1
                // for k far below lambda. Documented, reproduced, and harmless:
                // the value is compared against a cutoff, not against 1.
                assert!(up <= 1.0 + 1e-9, "lam={lam} k={k} up={up}");
            }
        }
    }

    #[test]
    fn large_lambda_lower_and_upper_agree_where_the_tail_is_visible() {
        // below lambda both tails are O(1) and must sum to 1
        for lam in [700.5f64, 900.0, 1200.0] {
            let k = (lam * 0.5) as u32;
            let lo = poisson_cdf(k, lam, true, false).unwrap();
            let up = poisson_cdf(k, lam, false, false).unwrap();
            assert!(
                (lo + up - 1.0).abs() < 1e-6,
                "lam={lam} k={k} lo={lo} up={up}"
            );
        }
    }

    #[test]
    fn ln_gamma_is_accurate() {
        // known values
        assert!((ln_gamma_ref(1.0) - 0.0).abs() < 1e-12);
        assert!((ln_gamma_ref(5.0) - 24.0_f64.ln()).abs() < 1e-12);
        assert!((ln_gamma_ref(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-12);
    }
}
