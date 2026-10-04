//! Port of `MACS3/Signal/VariantStat.py` -- the likelihood model behind
//! `callvar`'s allele-specific calls.
//!
//! # Where the numerics come from
//!
//! The Phred conversion is written the way upstream writes it, and the comment there
//! is worth repeating because it is easy to "simplify" into a different number:
//!
//! ```text
//! Phred = -10 log10 E   =>   1-E = 1 - exp(-Phred/10 * ln10)
//! ```
//!
//! `LN10` is spelled out as a literal (`2.3025850929940458`) rather than taken from
//! `std::f64::consts::LN_10`, and `LN10_tenth` as `0.23025850929940458` -- these are
//! upstream's literals and `LN_10/10` is not bit-identical to them. Every
//! `ln` below accumulates in `f64` in the same order as the Cython loop.
//!
//! # The `float` trap
//!
//! `max_allowed_ar` is a Cython `cython.float`, i.e. a **32-bit** float, and the
//! default `0.99` therefore arrives as `0.99f32` widened to `f64`
//! (`0.9900000095367432`), not as the `f64` nearest `0.99`. Worse, in
//! `calculate_ln` the expressions `r < 1 - max_allowed_r` and
//! `log(1 - max_allowed_r)` are evaluated in **C `float`** arithmetic -- Cython
//! emits `1 - max_allowed_r` with both operands float, so the subtraction rounds to
//! `f32` before `log` widens it. Using `f64` for the subtraction changes the last
//! bits of a term that is then compared against `r`, which decides whether the
//! clamped branch is taken at all. [`max_allowed_ar`] keeps the value as `f32` and
//! does the subtraction in `f32` for that reason.
//!
//! The `(int)` casts in [`calculate_gq`] are C casts -- truncation toward zero, not
//! a floor -- and a negative `-4.34294*log(tmp)` therefore rounds *up*. `0.99`
//! clamps are likewise applied before the `1e-110` floor, in that order.
//!
//! # Domain errors are part of the contract
//!
//! These functions run inside CPython, where `math.log1p(-1.0)` and `math.log(0.0)`
//! both raise `ValueError: math domain error`. Rust's `ln_1p`/`ln` return `-inf`
//! instead. That is not a cosmetic difference: a Phred quality of `0` gives
//! `e == 1.0`, so `log1p(-e)` is `log1p(-1.0)` and upstream *aborts the call*,
//! whereas a naive port would emit `-inf` into a VCF and carry on.
//!
//! Hence every logarithm below is checked before it is taken and the error is
//! propagated as `Err`, mirroring the `Exception` that upstream raises. A Phred
//! quality of 0 in a treatment BAM is not exotic -- it is what some aligners write
//! for a masked base.

// Upstream spells both constants out with more digits than an `f64` can hold, and
// uses `LN10_tenth` everywhere the algebra wants `LN10 / 10`. Neither `LN_10` nor
// `LN_10 / 10` is bit-identical to those literals, so the literals are transcribed
// verbatim and both relevant lints are muted with that reason. Upstream's `LN10` is
// defined but never referenced by the functions in this module, so only the tenth is
// carried over.
#[allow(clippy::excessive_precision)]
pub const LN10_TENTH: f64 = 0.23025850929940458;

/// The `-10 * log10` constant hard-coded in `calculate_GQ`.
///
/// It is *not* `-10/LN10`: upstream wrote the five-digit value `-4.34294`, and that
/// five-digit value is what produces the integer truncations in the VCF.
const NEG10_LOG10: f64 = -4.34294;

/// Upstream's `cython.float` default for the allele-ratio clamp.
pub const DEFAULT_MAX_ALLOWED_AR: f32 = 0.99;

/// The error text upstream's `log1p`/`log` domain failures surface as. Cython turns
/// CPython's `ValueError` into a generic `Exception`, so only the type is observable.
const DOMAIN_ERROR: &str = "math domain error";

/// `log1p(-exp(-q * LN10_tenth))`: the log probability that a `top1` base call is
/// correct, for base quality `q`.
///
/// `Err` when `e >= 1.0`, i.e. `q <= 0`, which is upstream's domain error.
///
/// # `ln_1p`, not `ln(1 - x)`
///
/// `log1p(-e)` and `ln(1 - e)` are equal in exact arithmetic and differ in the last
/// bits in `f64`. For a high-quality base `e` is about `1e-3`, so `1 - e` keeps ~13
/// significant digits while `log1p` keeps full precision on the *result* -- and the
/// oracle's golden caught the substitution immediately on the single-base cases. Use
/// `ln_1p`.
#[inline]
fn top1_ln(q: i32) -> Result<f64, String> {
    let e = (-(q as f64) * LN10_TENTH).exp();
    if e >= 1.0 {
        return Err(DOMAIN_ERROR.to_string());
    }
    Ok((-e).ln_1p())
}

/// `log(exp(-q * LN10_tenth))`, kept as its own function because upstream writes
/// it this way for `top2` rather than as a `log1p`. They agree mathematically;
/// the distinction is recorded so a future refactor does not silently swap them.
#[inline]
fn top2_ln(q: i32) -> f64 {
    // `log(exp(x))` is 0 at worst (x -> -inf), so this cannot hit a domain error:
    // the argument is in (0, 1].
    ((-(q as f64) * LN10_TENTH).exp()).ln()
}

/// `CalModel_Homo`: the homozygous model. No free parameters, so `BIC` is just
/// `-2 * lnL`.
pub fn cal_model_homo(
    top1_bq_t: &[i32],
    top1_bq_c: &[i32],
    top2_bq_t: &[i32],
    top2_bq_c: &[i32],
) -> Result<(f64, f64), String> {
    let mut ln_l = 0.0f64;
    for &q in top1_bq_t {
        ln_l += top1_ln(q)?;
    }
    for &q in top1_bq_c {
        ln_l += top1_ln(q)?;
    }
    for &q in top2_bq_t {
        ln_l += top2_ln(q);
    }
    for &q in top2_bq_c {
        ln_l += top2_ln(q);
    }
    let bic = -2.0 * ln_l;
    Ok((ln_l, bic))
}

/// `calculate_ln`: log likelihood of observing `k` top1 bases out of `tn`, given
/// an expected top1 ratio `r`.
///
/// The `max_allowed_ar` clamp mirrors upstream's ordering exactly: the clamp test
/// first, then the combinatorial term, then the per-base error terms.
// Eight parameters, because upstream's has eight. Bundling them into a struct would
// hide the `(m, n, tn)` / `(me, ne)` pairing that makes the function readable.
#[allow(clippy::too_many_arguments)]
fn calculate_ln(
    m: usize,
    n: usize,
    tn: usize,
    me: &[i32],
    ne: &[i32],
    r: f64,
    k: usize,
    max_allowed_ar: f32,
) -> Result<f64, String> {
    debug_assert_eq!(m + n, tn);
    let mut ln_l = 0.0f64;

    // Clamp test. `1 - max_allowed_ar` is evaluated in f32, matching C.
    let one_minus = (1.0f32 - max_allowed_ar) as f64;
    if r > max_allowed_ar as f64 || r < one_minus {
        ln_l += k as f64 * (max_allowed_ar as f64).ln();
        ln_l += (tn - k) as f64 * one_minus.ln();
    } else {
        ln_l += k as f64 * r.ln();
        ln_l += (tn - k) as f64 * (1.0 - r).ln();
    }

    // Entirely biased toward one allele: no combinatorial term.
    if k == 0 || k == tn {
        // pass
    } else if k <= tn / 2 {
        for i in 0..k {
            ln_l += ((tn - i) as f64 / (k - i) as f64).ln();
        }
    } else {
        for i in 0..(tn - k) {
            ln_l += ((tn - i) as f64 / (tn - k - i) as f64).ln();
        }
    }

    // The per-base terms. `e == 1` (a base quality of 0) drives the mixture weight
    // to 0 whenever `k` sits at one end, and `log(0)` is a domain error upstream.
    let kt = k as f64 / tn as f64;
    for &q in me.iter().take(m) {
        let e = (-(q as f64) * LN10_TENTH).exp();
        let arg = (1.0 - e) * kt + e * (1.0 - kt);
        if arg <= 0.0 {
            return Err(DOMAIN_ERROR.to_string());
        }
        ln_l += arg.ln();
    }

    for &q in ne.iter().take(n) {
        let e = (-(q as f64) * LN10_TENTH).exp();
        let arg = (1.0 - e) * (1.0 - kt) + e * kt;
        if arg <= 0.0 {
            return Err(DOMAIN_ERROR.to_string());
        }
        ln_l += arg.ln();
    }

    Ok(ln_l)
}

/// `GreedyMaxFunctionAS`: maximise the likelihood over `k` for a treatment sample
/// that *does* have allele preference, returning `(lnL, k, allele_ratio)`.
///
/// # The redundant `if btemp` branches
///
/// Upstream ends both sweep directions with `if btemp: ... else: ...` where the
/// two arms are character-for-character identical. That looks like an
/// upstream bug, and it probably is one, but it has no effect on the result, so it
/// is collapsed here rather than "fixed". Changing it would change nothing
/// observable and would cost a round trip to prove that.
#[allow(clippy::needless_range_loop)]
fn greedy_max_function_as(
    m: usize,
    n: usize,
    tn: usize,
    me: &[i32],
    ne: &[i32],
    max_allowed_ar: f32,
) -> Result<(f64, usize, f64), String> {
    assert_eq!(m + n, tn);

    if tn == 1 {
        let dl = calculate_ln(m, n, tn, me, ne, 0.0, 0, max_allowed_ar)?;
        let dr = calculate_ln(m, n, tn, me, ne, 1.0, 1, max_allowed_ar)?;
        return if dl > dr {
            Ok((dl, 0, 0.0))
        } else {
            Ok((dr, 1, 1.0))
        };
    }
    if m == 0 {
        return Ok((
            calculate_ln(m, n, tn, me, ne, 0.0, m, max_allowed_ar)?,
            m,
            1.0 - max_allowed_ar as f64,
        ));
    }
    if m == tn {
        return Ok((
            calculate_ln(m, n, tn, me, ne, 1.0, m, max_allowed_ar)?,
            m,
            max_allowed_ar as f64,
        ));
    }

    let k0 = m;
    let d0 = calculate_ln(m, n, tn, me, ne, k0 as f64 / tn as f64, k0, max_allowed_ar)?;
    let d1l = calculate_ln(
        m,
        n,
        tn,
        me,
        ne,
        (k0 - 1) as f64 / tn as f64,
        k0 - 1,
        max_allowed_ar,
    )?;
    let d1r = calculate_ln(
        m,
        n,
        tn,
        me,
        ne,
        (k0 + 1) as f64 / tn as f64,
        k0 + 1,
        max_allowed_ar,
    )?;

    if d0 > d1l - 1e-8 && d0 > d1r - 1e-8 {
        Ok((d0, k0, k0 as f64 / tn as f64))
    } else if d1l > d0 {
        let mut d_old = d1l;
        let mut k_old = k0 - 1;
        let mut r_old = (k0 - 1) as f64 / tn as f64;
        while k_old > 1 {
            let k_new = k_old - 1;
            let r_new = k_new as f64 / tn as f64;
            let d_new = calculate_ln(m, n, tn, me, ne, r_new, k_new, max_allowed_ar)?;
            if d_new - 1e-8 < d_old {
                break;
            }
            k_old = k_new;
            d_old = d_new;
            r_old = r_new;
        }
        Ok((d_old, k_old, r_old))
    } else if d1r > d0 {
        let mut d_old = d1r;
        let mut k_old = k0 + 1;
        let mut r_old = (k0 + 1) as f64 / tn as f64;
        while k_old < tn - 1 {
            let k_new = k_old + 1;
            let r_new = k_new as f64 / tn as f64;
            let d_new = calculate_ln(m, n, tn, me, ne, r_new, k_new, max_allowed_ar)?;
            if d_new - 1e-8 < d_old {
                break;
            }
            k_old = k_new;
            d_old = d_new;
            r_old = r_new;
        }
        Ok((d_old, k_old, r_old))
    } else {
        Err("error in GreedyMaxFunctionAS".to_string())
    }
}

/// `GreedyMaxFunctionNoAS`: maximise over `k` with the ratio pinned at the
/// background 0.5, for a sample assumed to have no allele preference.
#[allow(clippy::needless_range_loop)]
fn greedy_max_function_noas(
    m: usize,
    n: usize,
    tn: usize,
    me: &[i32],
    ne: &[i32],
) -> Result<(f64, usize), String> {
    let bg_r = 0.5f64;

    if tn == 1 {
        let dl = calculate_ln(m, n, tn, me, ne, bg_r, 0, DEFAULT_MAX_ALLOWED_AR)?;
        let dr = calculate_ln(m, n, tn, me, ne, bg_r, 1, DEFAULT_MAX_ALLOWED_AR)?;
        return if dl > dr { Ok((dl, 0)) } else { Ok((dr, 1)) };
    }
    if m == 0 || m == tn {
        return Ok((
            calculate_ln(m, n, tn, me, ne, bg_r, m, DEFAULT_MAX_ALLOWED_AR)?,
            m,
        ));
    }

    let k0 = m;
    let d0 = calculate_ln(m, n, tn, me, ne, bg_r, k0, DEFAULT_MAX_ALLOWED_AR)?;
    let d1l = calculate_ln(m, n, tn, me, ne, bg_r, k0 - 1, DEFAULT_MAX_ALLOWED_AR)?;
    let d1r = calculate_ln(m, n, tn, me, ne, bg_r, k0 + 1, DEFAULT_MAX_ALLOWED_AR)?;

    if d0 > d1l - 1e-8 && d0 > d1r - 1e-8 {
        Ok((d0, k0))
    } else if d1l > d0 {
        let mut d_old = d1l;
        let mut k_old = k0 - 1;
        while k_old >= 1 {
            let k_new = k_old - 1;
            let d_new = calculate_ln(m, n, tn, me, ne, bg_r, k_new, DEFAULT_MAX_ALLOWED_AR)?;
            if d_new - 1e-8 < d_old {
                break;
            }
            k_old = k_new;
            d_old = d_new;
        }
        Ok((d_old, k_old))
    } else if d1r > d0 {
        let mut d_old = d1r;
        let mut k_old = k0 + 1;
        // upstream: `while kold <= tn - 1`, i.e. `kold < tn`. Written with the `- 1`
        // so the correspondence with the source is checkable by eye.
        #[allow(clippy::int_plus_one)]
        while k_old <= tn - 1 {
            let k_new = k_old + 1;
            let d_new = calculate_ln(m, n, tn, me, ne, bg_r, k_new, DEFAULT_MAX_ALLOWED_AR)?;
            if d_new - 1e-8 < d_old {
                break;
            }
            k_old = k_new;
            d_old = d_new;
        }
        Ok((d_old, k_old))
    } else {
        Err("error in GreedyMaxFunctionNoAS".to_string())
    }
}

/// `CalModel_Heter_noAS`: heterogeneous model, no allele specificity on either
/// sample.
pub fn cal_model_heter_noas(
    top1_bq_t: &[i32],
    top1_bq_c: &[i32],
    top2_bq_t: &[i32],
    top2_bq_c: &[i32],
) -> Result<(f64, f64), String> {
    let mut ln_l = 0.0f64;
    let mut bic = 0.0f64;

    let tn_t = top1_bq_t.len() + top2_bq_t.len();
    if tn_t == 0 {
        return Err("Total number of treatment reads is 0!".to_string());
    }
    let (ln_l_t, _) =
        greedy_max_function_noas(top1_bq_t.len(), top2_bq_t.len(), tn_t, top1_bq_t, top2_bq_t)?;
    ln_l += ln_l_t;
    bic += -2.0 * ln_l_t;

    let tn_c = top1_bq_c.len() + top2_bq_c.len();
    if tn_c != 0 {
        let (ln_l_c, _) =
            greedy_max_function_noas(top1_bq_c.len(), top2_bq_c.len(), tn_c, top1_bq_c, top2_bq_c)?;
        ln_l += ln_l_c;
        bic += -2.0 * ln_l_c;
    }

    // `tn_T == 0` is unreachable here (it returned above), but the branch is kept
    // so the port stays a line-for-line image of the cost model.
    if tn_t == 0 {
        bic += (tn_c as f64).ln();
    } else if tn_c == 0 {
        bic += (tn_t as f64).ln();
    } else {
        bic += (tn_t as f64).ln() + (tn_c as f64).ln();
    }

    Ok((ln_l, bic))
}

/// `CalModel_Heter_AS`: heterogeneous model, allele-specific on the treatment and
/// explicitly assumed *not* allele-specific on the control.
pub fn cal_model_heter_as(
    top1_bq_t: &[i32],
    top1_bq_c: &[i32],
    top2_bq_t: &[i32],
    top2_bq_c: &[i32],
    max_allowed_ar: f32,
) -> Result<(f64, f64), String> {
    let mut ln_l = 0.0f64;
    let mut bic = 0.0f64;

    let tn_t = top1_bq_t.len() + top2_bq_t.len();
    if tn_t == 0 {
        return Err("Total number of treatment reads is 0!".to_string());
    }
    let (ln_l_t, _k_t, _ar) = greedy_max_function_as(
        top1_bq_t.len(),
        top2_bq_t.len(),
        tn_t,
        top1_bq_t,
        top2_bq_t,
        max_allowed_ar,
    )?;
    ln_l += ln_l_t;
    bic += -2.0 * ln_l_t;

    let tn_c = top1_bq_c.len() + top2_bq_c.len();
    if tn_c != 0 {
        let (ln_l_c, _k_c) =
            greedy_max_function_noas(top1_bq_c.len(), top2_bq_c.len(), tn_c, top1_bq_c, top2_bq_c)?;
        ln_l += ln_l_c;
        bic += -2.0 * ln_l_c;
    }

    if tn_t == 0 {
        bic += (tn_c as f64).ln();
    } else if tn_c == 0 {
        bic += 2.0 * (tn_t as f64).ln();
    } else {
        bic += 2.0 * (tn_t as f64).ln() + (tn_c as f64).ln();
    }

    Ok((ln_l, bic))
}

/// `calculate_GQ`: genotype quality from three log-likelihoods.
///
/// `(int)` is a C cast, so `-4.34294*log(tmp)` truncates toward zero. The clamps
/// are applied in upstream's order: cap at 1 first, then floor at `1e-110`.
#[allow(clippy::manual_clamp)]
pub fn calculate_gq(ln_l1: f64, ln_l2: f64, ln_l3: f64) -> i32 {
    let mut l2 = (ln_l2 - ln_l1).exp();
    let mut l3 = (ln_l3 - ln_l1).exp();

    if l2 > 1.0 {
        l2 = 1.0;
    }
    if l3 > 1.0 {
        l3 = 1.0;
    }
    if l2 < 1e-110 {
        l2 = 1e-110;
    }
    if l3 < 1e-110 {
        l3 = 1e-110;
    }

    let s = 1.0 + l2 + l3;
    let tmp = (l2 + l3) / s;
    if tmp > 1e-110 {
        (NEG10_LOG10 * tmp.ln()) as i32
    } else {
        255
    }
}

/// `calculate_GQ_heterASsig`: allele-specific significance score.
#[allow(clippy::manual_clamp)]
pub fn calculate_gq_heter_assig(ln_l1: f64, ln_l2: f64) -> i32 {
    let mut l2 = (ln_l2 - ln_l1).exp();
    if l2 > 1.0 {
        l2 = 1.0;
    }
    if l2 < 1e-110 {
        l2 = 1e-110;
    }
    let s = 1.0 + l2;
    let tmp = l2 / s;
    if tmp > 1e-110 {
        (NEG10_LOG10 * tmp.ln()) as i32
    } else {
        255
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn homo_likelihood_matches_closed_form() {
        // One top1 base at q=30 and nothing else. lnL is that base's log
        // probability of being correct.
        let q = [30i32];
        let (ln_l, bic) = cal_model_homo(&q, &[], &[], &[]).unwrap();
        let expect = (-(-(30.0f64) * LN10_TENTH).exp()).ln_1p();
        assert!((ln_l - expect).abs() < 1e-15, "{ln_l} vs {expect}");
        assert!((bic - (-2.0 * ln_l)).abs() < 1e-15);
    }

    #[test]
    fn empty_everything_is_homo_zero() {
        assert_eq!(cal_model_homo(&[], &[], &[], &[]).unwrap(), (0.0, 0.0));
    }

    #[test]
    fn a_zero_base_quality_is_a_domain_error_not_minus_infinity() {
        // Upstream runs in CPython, where `log1p(-1.0)` raises. A port that returned
        // -inf here would write `-inf` into a VCF's GQ calculation instead of
        // refusing the call, so the error is part of the contract.
        assert!(cal_model_homo(&[0], &[], &[], &[]).is_err());
        assert!(cal_model_heter_noas(&[0], &[], &[], &[]).is_err());
        assert!(cal_model_heter_as(&[0], &[], &[], &[], DEFAULT_MAX_ALLOWED_AR).is_err());
        // A zero quality on a *top2* base is fine: `log(exp(0)) == 0`.
        assert!(cal_model_homo(&[], &[], &[0], &[]).unwrap().0 == 0.0);
    }

    #[test]
    fn heter_noas_rejects_empty_treatment() {
        let e = cal_model_heter_noas(&[], &[10], &[], &[]).unwrap_err();
        assert_eq!(e, "Total number of treatment reads is 0!");
    }

    #[test]
    fn heter_as_rejects_empty_treatment() {
        let e = cal_model_heter_as(&[], &[10], &[], &[], DEFAULT_MAX_ALLOWED_AR).unwrap_err();
        assert_eq!(e, "Total number of treatment reads is 0!");
    }

    #[test]
    fn heter_noas_empty_control_is_allowed() {
        let t1 = [20i32, 22, 25];
        let t2 = [18i32, 19];
        let tn_t = (t1.len() + t2.len()) as f64;
        let (ln_l, bic) = cal_model_heter_noas(&t1, &[], &t2, &[]).unwrap();
        assert!(ln_l < 0.0);
        // cost = -2*lnL_T + ln(tn_T); with tn_C == 0 only the treatment is scored.
        let expected = -2.0 * ln_l + tn_t.ln();
        assert!((bic - expected).abs() < 1e-12, "{bic} vs {expected}");
    }

    #[test]
    fn single_read_is_handled_not_panicked() {
        let r = cal_model_heter_noas(&[30], &[], &[], &[]);
        assert!(r.is_ok());
    }

    #[test]
    fn all_top1_and_all_top2_are_symmetric_in_the_greedy_search() {
        let one = [30i32];
        let a = greedy_max_function_noas(1, 0, 1, &one, &[]).unwrap();
        let b = greedy_max_function_noas(0, 1, 1, &[], &one).unwrap();
        assert!((a.0 - b.0).abs() < 1e-15);
    }

    #[test]
    fn greedy_as_picks_the_majority_allele() {
        // 8 top1, 2 top2 -> k should land on the top1 side.
        let me: Vec<i32> = vec![30; 8];
        let ne: Vec<i32> = vec![30; 2];
        let (ln_l, k, ar) =
            greedy_max_function_as(8, 2, 10, &me, &ne, DEFAULT_MAX_ALLOWED_AR).unwrap();
        assert_eq!(k, 8, "should prefer the majority allele, got k={k}");
        assert!((ar - 0.8).abs() < 1e-12);
        assert!(ln_l < 0.0);
    }

    #[test]
    fn calculate_ln_clamps_the_allele_ratio() {
        // r above the clamp uses log(0.99) + (tn-k)*log(0.01), which is finite.
        let me = [30i32, 30, 30, 30];
        let ne = [30i32, 30];
        let low = calculate_ln(4, 2, 6, &me, &ne, 0.999, 4, DEFAULT_MAX_ALLOWED_AR).unwrap();
        assert!(low.is_finite());
    }

    #[test]
    fn gq_is_one_for_identical_likelihoods() {
        // Three equiprobable genotypes: tmp = 2/3, so GQ = (int)(1.761) = 1.
        // Upstream's floor is not a floor on GQ -- 1, not 0, is "uninformative".
        assert_eq!(calculate_gq(-10.0, -10.0, -10.0), 1);
    }

    #[test]
    fn gq_grows_as_alternatives_get_worse() {
        let a = calculate_gq(-10.0, -20.0, -30.0);
        let b = calculate_gq(-10.0, -40.0, -60.0);
        assert!(a < b, "worse alternatives must score higher: {a} vs {b}");
    }

    #[test]
    fn gq_truncates_toward_zero() {
        // A value that lands between two integers proves the cast is `int`, not floor.
        let (l1, l2, l3) = (-1.0f64, -2.0f64, -3.0f64);
        let v = calculate_gq(l1, l2, l3);
        let exact = NEG10_LOG10
            * (((l2 - l1).exp() + (l3 - l1).exp()) / (1.0 + (l2 - l1).exp() + (l3 - l1).exp()))
                .ln();
        assert_eq!(v as f64, exact.trunc(), "must truncate like C's (int)");
    }

    #[test]
    fn gq_assig_is_three_for_identical_likelihoods() {
        // tmp = 1/2, so (int)(-4.34294 * ln 0.5) = (int)3.0104 = 3.
        assert_eq!(calculate_gq_heter_assig(-10.0, -10.0), 3);
    }

    #[test]
    fn gq_assig_grows() {
        assert!(calculate_gq_heter_assig(-1.0, -50.0) > calculate_gq_heter_assig(-1.0, -5.0));
    }
}
