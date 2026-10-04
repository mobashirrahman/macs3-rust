//! Chi-square distribution functions for **even** degrees of freedom.
//!
//! Port of `MACS3/Signal/Prob.py` lines 142-249. Upstream requires an even `df`
//! and gives a wrong answer for odd `df`; [`chisq_pvalue_e`] returns a typed
//! error instead, which is strictly better and cannot change any MACS result
//! because MACS only ever passes even `df`.
//!
//! These are used by `callvar`'s variant model, not by peak calling.

use crate::util::logspace_add;
use crate::StatsError;

/// `bigx` from `Prob.py`: the crossover point between the two summation
/// strategies. It is `20`, not a "big number" in any usual sense — it is where
/// `exp(-a) * a^k` stops being representable in a straightforward way.
const BIGX: f64 = 20.0;

/// `exp(x)`, or `0` when `x < -20` (`ex20`).
#[inline]
fn ex20(x: f64) -> f64 {
    if x < -20.0 {
        0.0
    } else {
        x.exp()
    }
}

/// Chi-square upper tail for even `df`, in linear space (`chisq_pvalue_e`).
///
/// # Errors
/// [`StatsError::BadChiSquareDf`] when `df` is odd or `<= 1`, mirroring the
/// upstream docstring's "do not violate this assumption".
pub fn chisq_pvalue_e(x: f64, df: u32) -> Result<f64, StatsError> {
    if df == 0 || df % 2 == 1 {
        return Err(StatsError::BadChiSquareDf(df));
    }
    if x <= 0.0 {
        return Ok(1.0);
    }

    let mut x = x;
    let a = 0.5 * x;
    let y = ex20(-a);
    let mut s = y;

    if df > 2 {
        x = 0.5 * (df as f64 - 1.0);
        let mut z = 1.0_f64;
        if a > BIGX {
            // log-space accumulation
            let mut e = 0.0_f64;
            let c = a.ln();
            while z <= x {
                e += z.ln();
                s += ex20(c * z - a - e);
                z += 1.0;
            }
            return Ok(s);
        }
        // linear-space accumulation
        let mut e = 1.0_f64;
        let mut c = 0.0_f64;
        while z <= x {
            e *= a / z;
            c += e;
            z += 1.0;
        }
        return Ok(c * y + s);
    }

    Ok(s)
}

/// Chi-square upper tail for even `df`, in log space (`chisq_logp_e`).
///
/// Returns `-ln p`, or `-log10 p` when `log10` is set.
pub fn chisq_logp_e(x: f64, df: u32, log10: bool) -> Result<f64, StatsError> {
    if df == 0 || df % 2 == 1 {
        return Err(StatsError::BadChiSquareDf(df));
    }
    if x <= 0.0 {
        return Ok(0.0);
    }

    let mut x = x;
    let a = 0.5 * x;
    let y = (-a).exp();
    let mut s = -a;

    if df > 2 {
        x = 0.5 * (df as f64 - 1.0);
        let mut z = 1.0_f64;
        if a > BIGX {
            let mut e = 0.0_f64;
            let c = a.ln();
            while z <= x {
                e += z.ln();
                s = logspace_add(s, c * z - a - e);
                z += 1.0;
            }
        } else {
            let mut e = 1.0_f64;
            let mut c = 0.0_f64;
            while z <= x {
                e *= a / z;
                c += e;
                z += 1.0;
            }
            s = (y + c * y).ln();
        }
    }

    Ok(if log10 { -s / 10.0_f64.ln() } else { -s })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Independent reference for the chi-square upper tail with even df, via the
    /// regularised upper incomplete gamma function `Q(df/2, x/2)`.
    ///
    /// For integer `a = df/2` this is the exact finite sum
    /// `Q(a, x) = e^{-x} sum_{i<a} x^i / i!`, a completely different algorithm
    /// from upstream's scaled `exp`-of-`exp` loops.
    fn chisq_q_ref(x: f64, df: u32) -> f64 {
        let a = df as f64 / 2.0;
        let xx = x / 2.0;
        if a < 64.0 {
            // exact finite sum for integer a: P(a,xx) = 1 - e^{-xx} sum_{i<a} xx^i/i!
            let mut term = 1.0_f64;
            let mut acc = 1.0_f64;
            for i in 1..(a as u64) {
                term *= xx / i as f64;
                acc += term;
            }
            return (-xx).exp() * acc;
        }
        if xx < a + 1.0 {
            // lower regularised incomplete gamma P(a,x) by series
            let mut ap = a;
            let mut sum = 1.0 / a;
            let mut del = sum;
            for _ in 0..10_000 {
                ap += 1.0;
                del *= xx / ap;
                sum += del;
                if del.abs() < sum.abs() * 1e-17 {
                    break;
                }
            }
            // P(a, xx) = x^a e^{-x} / Gamma(a) * sum
            sum * (-xx + a * xx.ln() - ln_gamma_series(a)).exp()
        } else {
            // upper regularised incomplete gamma Q(a,x) by continued fraction
            const TINY: f64 = 1e-300;
            let mut b = xx + 1.0 - a;
            let mut c = 1.0 / TINY;
            let mut d = 1.0 / b;
            let mut h = d;
            for i in 1..10_000 {
                let an = -(i as f64) * (i as f64 - a);
                b += 2.0;
                d = an * d + b;
                if d.abs() < TINY {
                    d = TINY;
                }
                c = b + an / c;
                if c.abs() < TINY {
                    c = TINY;
                }
                d = 1.0 / d;
                let del = d * c;
                h *= del;
                if (del - 1.0).abs() < 1e-17 {
                    break;
                }
            }
            (-xx + a * xx.ln() - ln_gamma_series(a)).exp() * h
        }
    }

    /// Series expansion of `ln Gamma(a)` for a > 0.
    fn ln_gamma_series(a: f64) -> f64 {
        // ln Gamma(a) = ln Gamma(a+1) - ln a
        crate::poisson::ln_gamma_ref(a + 1.0) - a.ln()
    }

    #[test]
    fn pvalue_matches_independent_reference() {
        for df in [2u32, 4, 6, 10, 20, 50] {
            for x in [0.5, 1.0, 3.0, 10.0, 25.0, 39.0] {
                let got = chisq_pvalue_e(x, df).unwrap();
                let want = chisq_q_ref(x, df);
                let tol = 1e-9 * want.abs().max(1e-14);
                assert!(
                    (got - want).abs() < tol,
                    "df={df} x={x} got={got} want={want}"
                );
            }
        }
    }

    /// Upstream loses accuracy once `x/2 > bigx = 20`, because `ex20` forces
    /// `exp(t)` to exactly `0` for `t < -20` and the rescaled loops then start
    /// from a truncated seed. We reproduce the loss rather than "fixing" it,
    /// because `callvar` VCF output is defined by upstream's value.
    #[test]
    fn large_x_truncation_is_reproduced_and_bounded() {
        for (df, x) in [(4u32, 41.0), (4, 80.0), (10, 80.0), (50, 200.0)] {
            let got = chisq_pvalue_e(x, df).unwrap();
            let want = chisq_q_ref(x, df);
            let rel = (got - want).abs() / want.max(1e-300);
            // upstream can lose the value entirely (df=4, x=80 returns 0.0
            // where the truth is 1.7e-16); otherwise it is off by a few percent
            assert!(
                got == 0.0 || rel < 0.15,
                "df={df} x={x} got={got} want={want} rel={rel}"
            );
        }
        // and the error grows, but stays bounded, as x grows
        let a = chisq_pvalue_e(80.0, 4).unwrap();
        let b = chisq_pvalue_e(400.0, 4).unwrap();
        let ra = (a - chisq_q_ref(80.0, 4)).abs() / chisq_q_ref(80.0, 4);
        let rb = (b - chisq_q_ref(400.0, 4)).abs() / chisq_q_ref(400.0, 4);
        assert!(rb >= ra, "truncation error should not shrink: {ra} -> {rb}");
    }

    #[test]
    fn df_two_truncates_above_forty() {
        // upstream: chisq_pvalue_e(x, 2) == ex20(-x/2)
        assert_eq!(chisq_pvalue_e(40.0, 2).unwrap(), (-20.0f64).exp());
        assert_eq!(chisq_pvalue_e(40.1, 2).unwrap(), 0.0);
        assert_eq!(chisq_pvalue_e(1000.0, 2).unwrap(), 0.0);
    }

    #[test]
    fn logp_matches_negative_log_of_pvalue() {
        for df in [2u32, 4, 8, 30] {
            for x in [0.5, 5.0, 30.0, 100.0] {
                let p = chisq_pvalue_e(x, df).unwrap();
                let lp = chisq_logp_e(x, df, false).unwrap();
                if p > 0.0 {
                    assert!((lp + p.ln()).abs() < 1e-8, "df={df} x={x} p={p} logp={lp}");
                }
            }
        }
    }

    #[test]
    fn log10_flag_scales_by_log_10() {
        let a = chisq_logp_e(12.0, 4, false).unwrap();
        let b = chisq_logp_e(12.0, 4, true).unwrap();
        assert!((a / 10.0_f64.ln() - b).abs() < 1e-9);
    }

    #[test]
    fn non_positive_x_is_the_identity() {
        assert_eq!(chisq_pvalue_e(0.0, 4).unwrap(), 1.0);
        assert_eq!(chisq_pvalue_e(-1.0, 4).unwrap(), 1.0);
        assert_eq!(chisq_logp_e(0.0, 4, false).unwrap(), 0.0);
    }

    #[test]
    fn odd_df_is_rejected() {
        assert!(matches!(
            chisq_pvalue_e(3.0, 3),
            Err(StatsError::BadChiSquareDf(3))
        ));
        assert!(matches!(
            chisq_logp_e(3.0, 1, false),
            Err(StatsError::BadChiSquareDf(1))
        ));
        assert!(matches!(
            chisq_pvalue_e(3.0, 0),
            Err(StatsError::BadChiSquareDf(0))
        ));
    }

    #[test]
    fn large_df_is_finite_and_in_range() {
        for df in [100u32, 1000, 10_000] {
            for x in [10.0, 500.0, 5000.0] {
                let p = chisq_pvalue_e(x, df).unwrap();
                assert!((0.0..=1.0).contains(&p), "df={df} x={x} p={p}");
                let lp = chisq_logp_e(x, df, false).unwrap();
                assert!(lp.is_finite());
            }
        }
    }
}
