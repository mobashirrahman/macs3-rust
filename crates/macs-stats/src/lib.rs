//! MACS3-compatible statistical primitives.
//!
//! This is a **faithful line-by-line port of `MACS3/Signal/Prob.py`**
//! (upstream commit `c544319`, MACS3 3.0.5). It is not a "better" implementation:
//! every accumulation order, every early exit, every clamping, and every
//! rounding step is reproduced, because MACS3's peak boundaries are a
//! discontinuous function of these values.
//!
//! Where upstream is mathematically wrong or inconsistent, we reproduce it and
//! say so in a comment citing the upstream line. Fixing it here would silently
//! change peak coordinates.
//!
//! # Deviations from upstream, and why they are safe
//!
//! * Cython's `uint32_t` loop counters become `u32`; upstream's `k < 0` checks
//!   are dead code for unsigned arguments and are dropped.
//! * `assert` becomes [`StatsError`] instead of a panic, so malformed parameter
//!   combinations produce a typed error rather than aborting.
//! * Python's `round(x, 5)` becomes [`py_round`], which routes through decimal
//!   formatting and therefore keeps CPython's round-half-to-even behaviour.
//!
//! # Reference
//!
//! Cross-validation against SciPy is a gate requirement (G2), but SciPy is an
//! *independent* correctness oracle only. MACS3 parity is defined by MACS3.

#![forbid(unsafe_code)]
#![deny(missing_debug_implementations)]

pub mod binomial;
mod chisq;
mod normal;
mod numpy_rng;
mod pcg64;
mod poisson;
mod seed_sequence;
mod util;

pub use binomial::{
    binomial_cdf, binomial_cdf_inv, binomial_coef, binomial_pdf, binomial_sf, pduplication,
};
pub use chisq::{chisq_logp_e, chisq_pvalue_e};
pub use normal::{pnorm, pnorm2, pnorm2_f64, poz};
pub use numpy_rng::NumpyRng;
pub use pcg64::{Pcg64, Pcg64Rng};
pub use poisson::{
    factorial, ln_gamma_ref, poisson_cdf, poisson_cdf_inv, poisson_cdf_q_inv, poisson_pdf,
    poisson_pdf_log,
};
pub use seed_sequence::{randomstate_from_seed_sequence, SeedSequence};
pub use util::{logspace_add, py_round};

/// Errors raised by the statistical primitives.
///
/// These correspond one-to-one to upstream's `assert`/`raise` sites, so that
/// the differential harness can assert identical accept/reject behaviour.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum StatsError {
    /// `lambda` must be strictly positive.
    #[error("lambda must be > 0, however we got {0}")]
    NonPositiveLambda(f64),
    /// `lambda` must be below 740 for the inverse Poisson CDF.
    #[error("lambda must be < 740 for the inverse Poisson CDF, however we got {0}")]
    LambdaTooLarge(f64),
    /// A CDF argument outside `[0, 1]`.
    #[error("CDF must be >= 0 and <= 1, however we got {0}")]
    CdfOutOfRange(f64),
    /// `df` must be even and greater than 1 for the even-degree chi-square.
    #[error("chi-square requires an even degree of freedom greater than 1, however we got {0}")]
    BadChiSquareDf(u32),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn py_round_matches_cpython_half_to_even() {
        // CPython: round(0.5) == 0, round(1.5) == 2, round(2.5) == 2
        assert_eq!(py_round(0.5, 0), 0.0);
        assert_eq!(py_round(1.5, 0), 2.0);
        assert_eq!(py_round(2.5, 0), 2.0);
        assert_eq!(py_round(-0.5, 0), -0.0);
        assert_eq!(py_round(1.234_567_8, 5), 1.234_57);
        assert_eq!(py_round(-1.234_567_8, 5), -1.234_57);
    }

    #[test]
    fn logspace_add_is_symmetric_and_equals_log_sum_exp() {
        // `ln_1p` is *more* accurate than `(1 + x).ln()`, so the reference here
        // is allowed a few ULP of relative slack. The bit-exact requirement is
        // enforced against the oracle in tests/parity_vectors.rs.
        for (a, b) in [
            (-3.25_f64, -7.5_f64),
            (-1.0, -40.0),
            (-20.5, -17.7),
            (-100.0, -3.0),
            (0.0, 0.0),
        ] {
            let got = logspace_add(a, b);
            let want = (a.exp() + b.exp()).ln();
            let rel = (got - want).abs() / want.abs().max(1e-300);
            assert!(rel < 1e-14, "logspace_add({a}, {b}) = {got} vs {want}");
            // symmetric to within ULP
            assert_eq!(
                logspace_add(a, b).to_bits(),
                logspace_add(b, a).to_bits(),
                "logspace_add({a}, {b}) is not symmetric"
            );
        }
        // logspace_add(x, -inf) == x
        assert_eq!(logspace_add(-3.25_f64, f64::NEG_INFINITY), -3.25);
        assert_eq!(
            logspace_add(-3.25_f64, f64::NEG_INFINITY).to_bits(),
            (-3.25f64).to_bits()
        );
    }

    #[test]
    fn logspace_add_uses_ln_1p() {
        // `ln_1p` is not an optimisation here, it is a correctness requirement:
        // upstream calls `math.log1p`, and `(1 + x).ln()` is a different function
        // once `1 + x` has been rounded. `chisq_logp_e` accumulates ~24 of these,
        // which is exactly the 2-ULP divergence the golden vectors caught.
        for (a, b) in [
            (-20.0_f64, -38.0_f64),
            (-20.5_f64, -17.7_f64),
            (-3.25_f64, -7.5_f64),
            (5.0_f64, 5.0_f64),
        ] {
            let want = if a > b {
                a + (b - a).exp().ln_1p()
            } else {
                b + (a - b).exp().ln_1p()
            };
            assert_eq!(
                logspace_add(a, b).to_bits(),
                want.to_bits(),
                "logspace_add({a}, {b})"
            );
        }
    }
}
