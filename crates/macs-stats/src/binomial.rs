//! Binomial distribution functions.
//!
//! Port of `MACS3/Signal/Prob.py` lines 634-880.
//!
//! These matter more than they look: MACS3's `--keep-dup auto` computes the
//! maximum tolerated duplicate count as
//! `binomial_cdf_inv(1 - p, N_total, 1 / effective_genome_size)`, so a
//! one-step difference in [`binomial_cdf_inv`] changes the retained read set
//! and therefore every peak.

use crate::StatsError;

/// Binomial coefficient `C(n, k)` (`binomial_coef`).
///
/// Upstream's very first statement casts `mx + 1` through **`float32`** before
/// widening to `f64`. That is a real precision loss for large `n` and it is
/// reproduced here; "fixing" it would change `binomial_pdf` outputs for large
/// trial counts.
pub fn binomial_coef(n: i64, k: i64) -> f64 {
    let mn = k.min(n - k);
    if mn < 0 {
        return 0.0;
    }
    if mn == 0 {
        return 1.0;
    }
    let mx = k.max(n - k);
    // upstream: cnk: float64_t = <float32_t>(mx + 1)
    let mut cnk = ((mx + 1) as f32) as f64;
    for i in 2..=mn {
        cnk = cnk * ((mx + i) as f64) / (i as f64);
    }
    cnk
}

/// Binomial PDF (`binomial_pdf`, "by H. Gene Shin").
///
/// Computes the coefficient by the ratio product from the mode outwards, with
/// the two rescaling loops that keep the running product inside
/// `[1e-100, 1e+100]`.
pub fn binomial_pdf(x: i64, a: i64, b: f64) -> f64 {
    if a < 1 {
        return 0.0;
    }
    if x < 0 || a < x {
        return 0.0;
    }
    if b == 0.0 {
        return if x == 0 { 1.0 } else { 0.0 };
    }
    if b == 1.0 {
        return if x == a { 1.0 } else { 0.0 };
    }

    let (p, mn, mx) = if x > a - x {
        (1.0 - b, a - x, x)
    } else {
        (b, x, a - x)
    };

    let mut pdf = 1.0_f64;
    let mut t: i64 = 0;

    for q in 1..=mn {
        pdf *= ((a - q + 1) as f64) * p / ((mn - q + 1) as f64);
        if pdf < 1e-100 {
            while pdf < 1e-3 {
                pdf /= 1.0 - p;
                t -= 1;
            }
        }
        if pdf > 1e+100 {
            while pdf > 1e3 && t < mx {
                pdf *= 1.0 - p;
                t += 1;
            }
        }
    }

    for _ in 0..(mx - t) {
        pdf *= 1.0 - p;
    }

    pdf
}

/// Binomial CDF (`binomial_cdf`). `lower = true` is the default and is what
/// MACS uses.
pub fn binomial_cdf(x: i64, a: i64, b: f64, lower: bool) -> f64 {
    if lower {
        binomial_cdf_f(x, a, b)
    } else {
        binomial_cdf_r(x, a, b)
    }
}

/// Binomial survival function (`binomial_sf`) = `1 - CDF`.
pub fn binomial_sf(x: i64, a: i64, b: f64, lower: bool) -> f64 {
    if lower {
        1.0 - binomial_cdf_f(x, a, b)
    } else {
        1.0 - binomial_cdf_r(x, a, b)
    }
}

/// Lower-tail binomial CDF (`_binomial_cdf_f`).
fn binomial_cdf_f(x: i64, a: i64, b: f64) -> f64 {
    // upstream: argmax = <int32_t>(a * b)  -- a * b is a *double* multiply
    // that is then TRUNCATED to int32. The truncation matters.
    let argmax = (a as f64 * b) as i32 as i64;

    if x < 0 {
        return 0.0;
    }
    if a < x {
        return 1.0;
    }
    if b == 0.0 {
        return 1.0;
    }
    if b == 1.0 {
        return 0.0;
    }

    if x > argmax {
        let seedpdf = binomial_pdf(argmax, a, b);
        let mut pdf = seedpdf;
        let mut cdf = pdf;
        for i in (0..argmax).rev() {
            pdf /= ((a - i) as f64) * b / (1.0 - b) / ((i + 1) as f64);
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
        }

        pdf = seedpdf;
        for i in argmax..x {
            pdf *= ((a - i) as f64) * b / (1.0 - b) / ((i + 1) as f64);
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
        }
        cdf.min(1.0)
    } else {
        let mut pdf = binomial_pdf(x, a, b);
        let mut cdf = pdf;
        for i in (0..x).rev() {
            pdf /= ((a - i) as f64) * b / (1.0 - b) / ((i + 1) as f64);
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
        }
        cdf.min(1.0)
    }
}

/// Upper-tail binomial CDF (`_binomial_cdf_r`).
fn binomial_cdf_r(x: i64, a: i64, b: f64) -> f64 {
    let argmax = (a as f64 * b) as i32 as i64;

    if x < 0 {
        return 1.0;
    }
    if a < x {
        return 0.0;
    }
    if b == 0.0 {
        return 0.0;
    }
    if b == 1.0 {
        return 1.0;
    }

    if x < argmax {
        let seedpdf = binomial_pdf(argmax, a, b);
        let mut pdf = seedpdf;
        let mut cdf = pdf;
        // upstream: range(argmax - 1, x, -1)  ==  [argmax-1 ..= x+1]
        for i in (x + 1..argmax).rev() {
            pdf /= ((a - i) as f64) * b / (1.0 - b) / ((i + 1) as f64);
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
        }

        pdf = seedpdf;
        let mut i = argmax;
        loop {
            pdf *= ((a - i) as f64) * b / (1.0 - b) / ((i + 1) as f64);
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
            i += 1;
        }
        cdf.min(1.0)
    } else {
        let mut pdf = binomial_pdf(x + 1, a, b);
        let mut cdf = pdf;
        let mut i = x + 1;
        loop {
            pdf *= ((a - i) as f64) * b / (1.0 - b) / ((i + 1) as f64);
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
            i += 1;
        }
        cdf.min(1.0)
    }
}

/// Inverse lower-tail binomial CDF (`binomial_cdf_inv`).
///
/// Returns the smallest `x` in `0..=a` with `P(X <= x) > cdf`, else `a`.
///
/// This is the function `--keep-dup auto` is built on.
///
/// # Errors
/// `cdf` outside `[0, 1]`.
pub fn binomial_cdf_inv(cdf: f64, a: i64, b: f64) -> Result<i64, StatsError> {
    if !(0.0..=1.0).contains(&cdf) {
        return Err(StatsError::CdfOutOfRange(cdf));
    }
    let mut cdf2 = 0.0_f64;
    for x in 0..=a {
        let pdf = binomial_pdf(x, a, b);
        cdf2 += pdf;
        if cdf < cdf2 {
            return Ok(x);
        }
    }
    Ok(a)
}

/// Probability of a duplicate fragment given a background PMF
/// (`pduplication`).
///
/// Two `float32` behaviours are load-bearing, and both are easy to miss:
///
/// 1. The loop variable is declared `float32_t` upstream
///    ([`MACS3/Signal/Prob.py:691`]: `p: float32_t`), so **every element of the
///    pmf is truncated to `f32` before it reaches `binomial_sf`**, which itself
///    takes an `f64`. The truncation is cubic in `b`, so it moves the result by
///    more than an ULP.
/// 2. The accumulator `sf` is `float32_t`, and C evaluates `f32 += f64` by
///    promoting, adding, and truncating back — so each partial sum is rounded to
///    `f32`.
///
/// The return value is also `float32_t`, divided by `f32(n)`.
pub fn pduplication(pmf: &[f64], n_obs: i64) -> f32 {
    let n = pmf.len();
    let mut sf: f32 = 0.0;
    for &p in pmf {
        // upstream truncates the pmf element to float32 here
        let b = p as f32 as f64;
        sf = ((sf as f64) + binomial_sf(2, n_obs, b, true)) as f32;
    }
    sf / n as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Independent reference for the binomial PMF, via `exp(lgamma)`.
    fn binomial_pdf_ref(x: i64, a: i64, b: f64) -> f64 {
        if x < 0 || x > a {
            return 0.0;
        }
        let ln = crate::poisson::ln_gamma_ref((a + 1) as f64)
            - crate::poisson::ln_gamma_ref((x + 1) as f64)
            - crate::poisson::ln_gamma_ref((a - x + 1) as f64)
            + (x as f64) * b.ln()
            + ((a - x) as f64) * (1.0 - b).ln();
        ln.exp()
    }

    #[test]
    fn pdf_matches_lgamma_reference() {
        for a in [1i64, 2, 5, 20, 100, 1000] {
            for b in [0.01, 0.1, 0.5, 0.9, 0.99] {
                for x in [0i64, 1, a / 2, a - 1, a] {
                    if x < 0 || x > a {
                        continue;
                    }
                    let got = binomial_pdf(x, a, b);
                    let want = binomial_pdf_ref(x, a, b);
                    let tol = 1e-9 * want.abs().max(1e-12) + 1e-300;
                    assert!(
                        (got - want).abs() <= tol,
                        "a={a} b={b} x={x} got={got} want={want}"
                    );
                }
            }
        }
    }

    #[test]
    fn cdf_sums_to_one() {
        for a in [1i64, 3, 10, 50, 200] {
            for b in [0.05, 0.5, 0.8] {
                let total: f64 = (0..=a).map(|x| binomial_pdf(x, a, b)).sum();
                assert!((total - 1.0).abs() < 1e-9, "a={a} b={b} total={total}");
            }
        }
    }

    #[test]
    fn cdf_is_monotone_and_bounded() {
        for a in [1i64, 5, 40, 137] {
            for b in [0.02, 0.3, 0.7, 0.98] {
                let mut prev = -1.0;
                for x in -1..=(a + 1) {
                    let c = binomial_cdf(x, a, b, true);
                    assert!((0.0..=1.0).contains(&c), "a={a} b={b} x={x} c={c}");
                    assert!(c >= prev - 1e-12, "a={a} b={b} x={x} not monotone");
                    prev = c;
                }
            }
        }
    }

    #[test]
    fn cdf_and_sf_are_complementary() {
        for a in [1i64, 7, 33, 100] {
            for b in [0.1, 0.5, 0.9] {
                for x in 0..=a {
                    let lo = binomial_cdf(x, a, b, true);
                    let up = binomial_cdf(x, a, b, false);
                    // `_binomial_cdf_f` returns P(X <= x) and `_binomial_cdf_r`
                    // returns P(X >= x + 1), so they sum to 1.
                    assert!(
                        (lo + up - 1.0).abs() < 1e-9,
                        "a={a} b={b} x={x} lo={lo} up={up}"
                    );
                    let sf = binomial_sf(x, a, b, true);
                    assert!((sf - (1.0 - lo)).abs() < 1e-15);
                }
            }
        }
    }

    #[test]
    fn cdf_matches_a_naive_binomial_sum() {
        // independent reference: sum the PMF directly
        for a in [1i64, 4, 9, 25] {
            for b in [0.2, 0.5, 0.75] {
                for x in 0..=a {
                    let want: f64 = (0..=x).map(|i| binomial_pdf(i, a, b)).sum();
                    let got = binomial_cdf(x, a, b, true);
                    assert!(
                        (got - want).abs() < 1e-12,
                        "a={a} b={b} x={x} got={got} want={want}"
                    );
                }
            }
        }
    }

    #[test]
    fn degenerate_parameters() {
        assert_eq!(binomial_pdf(0, 0, 0.5), 0.0); // a < 1
        assert_eq!(binomial_pdf(0, 5, 0.0), 1.0);
        assert_eq!(binomial_pdf(1, 5, 0.0), 0.0);
        assert_eq!(binomial_pdf(5, 5, 1.0), 1.0);
        assert_eq!(binomial_pdf(4, 5, 1.0), 0.0);
        assert_eq!(binomial_cdf(-1, 5, 0.5, true), 0.0);
        assert_eq!(binomial_cdf(6, 5, 0.5, true), 1.0);
        assert_eq!(binomial_cdf(-1, 5, 0.5, false), 1.0);
        assert_eq!(binomial_cdf(6, 5, 0.5, false), 0.0);
    }

    #[test]
    fn cdf_inv_brackets_the_cdf() {
        for a in [10i64, 100, 1000] {
            for b in [1e-4, 1e-6, 1e-8] {
                for target in [0.5, 0.9, 0.99, 0.999] {
                    let x = binomial_cdf_inv(target, a, b).unwrap();
                    assert!(x >= 0 && x <= a, "a={a} b={b} x={x}");
                    let c_at_x = binomial_cdf(x, a, b, true);
                    assert!(
                        c_at_x > target,
                        "a={a} b={b} target={target} x={x} c={c_at_x}"
                    );
                    if x > 0 {
                        let c_below = binomial_cdf(x - 1, a, b, true);
                        assert!(
                            c_below <= target,
                            "a={a} b={b} target={target} x={x} c_below={c_below}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn cdf_inv_rejects_out_of_range() {
        assert!(binomial_cdf_inv(-0.1, 10, 0.5).is_err());
        assert!(binomial_cdf_inv(1.1, 10, 0.5).is_err());
    }

    #[test]
    fn binomial_coef_small_values() {
        assert_eq!(binomial_coef(5, 0), 1.0);
        assert_eq!(binomial_coef(5, 5), 1.0);
        assert!((binomial_coef(5, 2) - 10.0).abs() < 1e-12);
        assert!((binomial_coef(10, 5) - 252.0).abs() < 1e-10);
        assert_eq!(binomial_coef(3, -1), 0.0);
    }

    #[test]
    fn pduplication_is_in_unit_range() {
        let pmf = vec![1e-6, 1e-6, 2e-6, 1e-7];
        let p = pduplication(&pmf, 100);
        assert!((0.0..=1.0).contains(&p), "got {p}");
        // more observations => more likely to see >= 2 duplicates at a position
        assert!(pduplication(&pmf, 1000) > pduplication(&pmf, 10));
    }
}
