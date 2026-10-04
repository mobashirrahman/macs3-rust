//! F238: the wide-window control factors are built from an **f32** mean fragment
//! length, because upstream's `tp.d` is a C float.
//!
//! `Parser.py:1496` reads `self.d = cython.cast(cython.float, m) / i`. Both operands
//! are C floats, so the division is `divss` and the result is a float32. That value
//! becomes `options.tsize`, and `PeakDetect.py:209` then evaluates
//!
//! ```python
//! tmp_v = float(self.d)/self.lregion*self.ratio_treat2control
//! ```
//!
//! An f64 mean here is off by one ulp in the factor, which changes nothing in the
//! coordinates and flips the fifth decimal of the control lambda.
//!
//! The numbers below are from `sweep/gmini_mfrag_d1200_w180_ctrl_042_B`, where the
//! control depth at the failing row is exactly 2194 (a counted control's depth is an
//! integer sum of counts, so it is exact under either rounding).

/// `float(self.d)` for the fixture's unweighted mean, as upstream stores it.
const MEAN_F32: f32 = 143.555;
/// `ratio_treat2control`, itself the product of two f32-rounded quantities.
const RATIO: f64 = 0.500_000_363_204_480_5;
/// `ctrl_d_s[2]` is `lregion` = 10000.
const LREGION: f64 = 10_000.0;

#[test]
fn control_factor_is_one_ulp_lower_from_the_f32_mean() {
    let f64_mean = 143.555_f64;
    let as_f32 = f64::from(MEAN_F32);

    // the f64 mean puts the factor one ulp high ...
    let wide_f64 = (f64_mean / LREGION * RATIO) as f32;
    // ... and the f32 mean is what upstream computes
    let wide_f32 = (as_f32 / LREGION * RATIO) as f32;

    assert!(
        wide_f32.to_bits() < wide_f64.to_bits(),
        "f32 mean must give the smaller factor: {wide_f32:?} vs {wide_f64:?}"
    );

    // at depth 2194 that ulp is the whole story: 0x417bf7ca prints as 15.74800
    // (ours) and 0x417bf7c9 prints as 15.74799 (upstream).
    let depth = 2194.0_f32;
    assert_eq!((depth * wide_f64).to_bits(), 0x417b_f7ca);
    assert_eq!((depth * wide_f32).to_bits(), 0x417b_f7c9);
    assert_eq!(format!("{:.5}", depth * wide_f32), "15.74799");
}
