//! Normal distribution helpers.
//!
//! Port of `MACS3/Signal/Prob.py` lines 42-138.
//!
//! [`pnorm`] and [`pnorm2`] are used by `predictd`'s strand model, and they
//! **deliberately truncate to `f32` internally**. The sequence of operations is
//! reproduced exactly: upstream casts `v`, `x - u` and `v` to `float32`, does
//! the arithmetic in `float32` (`m**2` is a C `float` multiply, and the whole
//! expression is a C `float`), then returns it widened to `f64`. Changing any
//! of that to `f64` changes the strand profiles that feed the cross-correlation
//! and therefore the predicted `d`.
//!
//! [`poz`] is the classic Abramowitz & Stegun 26.2.17 approximation to
//! `Phi(|z|)`; it saturates at 1.0 for `|z| >= 6`, which is upstream's stated
//! "maximum meaningful z value".

/// `2 * pi`, written exactly as the upstream source writes it.
///
/// This is bit-identical to `f64::consts::TAU` (it is the shortest round-trip
/// decimal for 2*pi), but the literal is kept so that a future "cleanup" to a
/// symbolic constant cannot be mistaken for a change to the numerics. Both
/// `pnorm` and `pnorm2` fold this into an expression whose result is stored at
/// `f32` precision, so any difference at all would show up in `predictd`.
#[allow(clippy::approx_constant)]
const TWO_PI: f64 = 6.283_185_307_179_586;

/// Probability of `X = x` when `X ~ Norm(u, v)` with integer `x`, `u`, `v`.
///
/// The `f32` casts on the *inputs* are load-bearing for large integers (anything
/// above 2^24 is rounded), but the **arithmetic is `f64`**: Cython evaluates
/// `sqrt(...)` and `exp(...)` from double-precision C functions, and `m**2` on
/// a C `float` is promoted to a `double` multiply. Reproducing this as pure
/// `f32` arithmetic loses 9 significant digits and changes `predictd`'s strand
/// profiles.
pub fn pnorm(x: i32, u: i32, v: i32) -> f64 {
    let i = v as f32 as f64;
    let m = (x - u) as f32 as f64;
    let n = v as f32 as f64;
    // 6.283185307179586 is 2 * pi
    1.0 / (TWO_PI * i).sqrt() * (-(m * m) / (2.0 * n)).exp()
}

/// Probability of `X = x` when `X ~ Norm(u, v)` with `f32` parameters.
///
/// **The arithmetic is `f64`; only the result is truncated.** Cython computes
/// `sqrt(6.283185307179586 * v)` and `exp(...)` with double-precision C library
/// functions (the `2.0` literal promotes `v` to `double`), and then stores the
/// product into a `float32_t` return variable. So the observable value is
/// `f32(f64_result)` — visible as a double whose low 32 bits are zero.
///
/// Doing the whole expression in `f32` gives a different answer in the 8th
/// significant digit, which shifts `predictd`'s strand profiles.
pub fn pnorm2(x: f32, u: f32, v: f32) -> f32 {
    // 6.283185307179586 is 2 * pi
    let i = v as f64;
    let d = x as f64 - u as f64;
    (1.0 / (TWO_PI * i).sqrt() * (-(d * d) / (2.0 * i)).exp()) as f32
}

/// The same value widened to `f64`, which is what a Python caller observes.
#[inline]
pub fn pnorm2_f64(x: f32, u: f32, v: f32) -> f64 {
    pnorm2(x, u, v) as f64
}

/// Evaluate `(((c0 * x + c1) * x + c2) ... )` — Horner, in the given order.
///
/// Upstream writes these polynomials as one deeply nested literal. This is the
/// same sequence of multiplications and additions in the same order (Horner
/// evaluation is exactly a nested left-associative polynomial), so the result is
/// bit-identical while being auditable.
#[inline]
fn horner(coeffs: &[f64], x: f64) -> f64 {
    let mut acc = coeffs[0];
    for &c in &coeffs[1..] {
        acc = acc * x + c;
    }
    acc
}

/// A&S 26.2.17 coefficients for the `|z| < 1` branch of `poz`.
const POZ_SMALL: [f64; 9] = [
    0.000_124_818_987,
    -0.001_075_204_047,
    0.005_198_775_019,
    -0.019_198_292_004,
    0.059_054_035_642,
    -0.151_968_751_364,
    0.319_152_932_694,
    -0.531_923_007_300,
    0.797_884_560_593,
];

/// A&S 26.2.17 coefficients for the `1 <= |z| < 6` branch of `poz`, evaluated at
/// `y - 2`.
const POZ_LARGE: [f64; 15] = [
    -0.000_045_255_659,
    0.000_152_529_290,
    -0.000_019_538_132,
    -0.000_676_904_986,
    0.001_390_604_284,
    -0.000_794_620_820,
    -0.002_034_254_874,
    0.006_549_791_214,
    -0.010_557_625_006,
    0.011_630_447_319,
    -0.009_279_453_341,
    0.005_353_579_108,
    -0.002_141_268_741,
    0.000_535_310_849,
    0.999_936_657_524,
];

/// Probability of a normal `z` value: `Phi(z)` (`poz`).
///
/// A&S 26.2.17 with a `Z_MAX = 6.0` saturation. The branch structure and the
/// Horner coefficient order are reproduced literally; reassociating the
/// polynomial changes the last two ULPs and `predictd`'s strand peaks are found
/// by thresholding, so those ULPs are not harmless.
pub fn poz(z: f64) -> f64 {
    const Z_MAX: f64 = 6.0;

    let x = if z == 0.0 {
        0.0
    } else {
        let y = 0.5 * z.abs();
        if y >= (Z_MAX * 0.5) {
            1.0
        } else if y < 1.0 {
            horner(&POZ_SMALL, y * y) * y * 2.0
        } else {
            horner(&POZ_LARGE, y - 2.0)
        }
    };

    if z > 0.0 {
        (x + 1.0) * 0.5
    } else {
        (1.0 - x) * 0.5
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Independent reference: `Phi` by a high-accuracy rational approximation
    /// (Zelen & Severo 26.2.17 in `erf`, evaluated in `f64`).
    fn phi_ref(z: f64) -> f64 {
        // Abramowitz & Stegun 7.1.26, |error| < 1.5e-7
        let negative = z < 0.0;
        let z = z.abs();
        let t = 1.0 / (1.0 + 0.231_641_9 * z);
        let d = 0.398_942_280_401_432_7_f64 * (-z * z / 2.0).exp();
        let p = t
            * (0.319_381_530
                + t * (-0.356_563_782
                    + t * (1.781_477_937 + t * (-1.821_255_978 + t * 1.330_274_429))));
        let tail = d * p;
        if negative {
            tail
        } else {
            1.0 - tail
        }
    }

    #[test]
    fn poz_matches_reference_within_upstreams_stated_accuracy() {
        for z in [
            -8.0, -6.0, -5.5, -3.0, -1.0, -0.5, 0.0, 0.5, 1.0, 3.0, 5.5, 6.0, 8.0,
        ] {
            let got = poz(z);
            let want = phi_ref(z);
            assert!((got - want).abs() < 1.5e-7, "z={z} got={got} want={want}");
        }
    }

    #[test]
    fn poz_saturates_and_is_symmetric() {
        assert_eq!(poz(6.0), 1.0);
        assert_eq!(poz(100.0), 1.0);
        assert_eq!(poz(0.0), 0.5);
        for z in [0.1, 0.9, 1.7, 2.8, 4.2] {
            assert_eq!(poz(z) + poz(-z), 1.0, "z={z}");
        }
    }

    #[test]
    fn poz_is_monotone() {
        let mut prev = -1.0;
        let mut z = -8.0;
        while z <= 8.0 {
            let cur = poz(z);
            assert!(cur >= prev, "not monotone at z={z}");
            prev = cur;
            z += 0.05;
        }
    }

    #[test]
    fn pnorm_is_f64_arithmetic_with_f32_inputs() {
        // u = v = 0 -> division by zero; the result must be NaN, not a panic
        assert!(pnorm(0, 0, 0).is_nan());
        // the mode of a N(10, 4) discrete mass
        let p = pnorm(10, 10, 4);
        assert!((p - 0.199_471_14).abs() < 1e-7, "got {p}");
        // full double precision, NOT the f32 rounding of it
        assert_ne!(p, p as f32 as f64);
        // the f32 casts on the inputs are real: `v` above 2^24 is rounded
        assert_eq!(pnorm(0, 0, 16_777_216), pnorm(0, 0, 16_777_217));
        assert_ne!(pnorm(0, 0, 16_777_216), pnorm(0, 0, 16_777_232));
    }

    #[test]
    fn pnorm2_is_f64_arithmetic_truncated_to_f32() {
        // the computation is f64; only the store is f32
        let f64_answer = 0.398_942_280_401_432_7_f64;
        let p = pnorm2(0.0, 0.0, 1.0);
        assert_eq!(p.to_bits(), (f64_answer as f32).to_bits(), "got {p}");
        // a Python caller sees that f32 widened, so the low 29 bits of the
        // double are zero (f32 and f64 use different exponent biases, so the
        // high words are not comparable bit-for-bit)
        let widened = pnorm2_f64(0.0, 0.0, 1.0);
        assert_eq!(widened.to_bits() & 0x1FFF_FFFF, 0);
        assert_eq!(widened, f64::from(p));
        // and therefore differs from the untruncated f64 answer
        assert_ne!(widened, f64_answer);
    }

    #[test]
    fn pnorm_integrates_to_one() {
        // sum over x of N(x; 0, sigma^2) == 1 for the integer grid
        for sigma in [1i32, 2, 3] {
            let total: f64 = (-20..=20).map(|x| pnorm(x, 0, sigma * sigma)).sum();
            assert!((total - 1.0).abs() < 1e-5, "sigma={sigma} total={total}");
        }
    }
}
