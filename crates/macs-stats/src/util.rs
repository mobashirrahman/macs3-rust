//! Small numeric helpers shared by the distribution functions.

/// Addition in log space: `log(exp(logx) + exp(logy))`.
///
/// Port of upstream `Prob.py:logspace_add`. Two details are load-bearing:
///
/// * The branch is on `>` (not `>=`), matching upstream exactly, which decides
///   which value is returned unchanged when the two arguments are equal.
/// * It uses **`ln_1p`, not `ln(1 + x)`**. Upstream calls `math.log1p`, which
///   is accurate to within an ULP for small arguments; `(1.0 + x).ln()` loses
///   several ULPs there. `chisq_logp_e` sums ~24 of these, which is enough for
///   the 2-ULP divergence the golden vectors caught.
#[inline]
pub fn logspace_add(logx: f64, logy: f64) -> f64 {
    if logx > logy {
        logx + (logy - logx).exp().ln_1p()
    } else {
        logy + (logx - logy).exp().ln_1p()
    }
}

/// Equivalent of CPython's `round(x, ndigits)`: decimal rounding, ties to even,
/// returning a float.
///
/// Formatting and re-parsing (rather than a manual `pow(10, n)` scale) is what
/// makes this bit-identical to CPython: CPython also rounds via decimal
/// conversion, so `round(2.675, 2) == 2.67` there, and a naive
/// `(x * 100.0).round() / 100.0` would give `2.68`.
///
/// Upstream uses this in `log10_poisson_cdf_*`, so the q-score value is
/// quantised to 5 decimal digits before it ever reaches a comparison.
#[inline]
pub fn py_round(x: f64, ndigits: i32) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let n = ndigits.clamp(0, 17) as usize;
    format!("{x:.n$}").parse::<f64>().unwrap_or(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decimal_ties_round_to_even_like_cpython() {
        // The classic binary-representation case CPython also gets "wrong".
        assert_eq!(py_round(2.675, 2), 2.67);
        assert_eq!(py_round(0.125, 2), 0.12);
        assert_eq!(py_round(0.375, 2), 0.38);
    }

    #[test]
    fn nan_and_inf_pass_through() {
        assert!(py_round(f64::NAN, 3).is_nan());
        assert_eq!(py_round(f64::INFINITY, 3), f64::INFINITY);
    }
}
