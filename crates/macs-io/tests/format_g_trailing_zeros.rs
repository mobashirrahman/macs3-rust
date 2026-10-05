//! `%.6g` strips trailing zeros in **exponent** form too, not only in fixed form.
//!
//! `PeakIO` writes every score column with `%.6g` (`write_to_narrowPeak`,
//! `write_to_broadPeak`, `write_to_summit_bed`), and C's `%g` removes trailing zeros
//! from the fraction *in both* the fixed and the exponent presentation. This port only
//! stripped them in the fixed branch, so a broad peak whose `-log10(qvalue)` is exactly
//! `7.25020e-06` was written `7.25020e-06` where upstream writes `7.2502e-06`. That is
//! not a theoretical difference: 1892 of the 1,132,451 rows of
//! `callpeak --broad -p 0.01` on the 5 M-read fixture carried such a value.
//!
//! Every expectation below is `("%.6g" % v)` from the pinned oracle's Python.

use macs_io::peakout::format_g;

#[test]
fn exponent_form_drops_trailing_zeros() {
    // the value that showed the bug up
    assert_eq!(format_g(7.25020e-06, 6), "7.2502e-06");
    assert_eq!(format_g(1.22680e-05, 6), "1.2268e-05");
    assert_eq!(format_g(4.23300e-06, 6), "4.233e-06");
}

#[test]
fn a_mantissa_that_is_only_zeros_keeps_its_leading_digit() {
    assert_eq!(format_g(1.00000e-16, 6), "1e-16");
    assert_eq!(format_g(1.00000e-05, 6), "1e-05");
}

/// Note: `%.6g` chooses fixed vs exponent form on the *rounded* value, so
/// `9.999995e-05` is `0.0001` and not `1e-04`. This port decides on the raw
/// exponent and gets that one value wrong; it is a separate defect from the trailing
/// zeros fixed here, so the case is left out rather than pinned as correct.
#[test]
fn the_fixed_and_exponent_branches_agree_with_c() {
    for (v, want) in [
        (7.25020e-06, "7.2502e-06"),
        (1.0, "1"),
        (0.0, "0"),
        (1.0 / 3.0, "0.333333"),
        (1234567.0, "1.23457e+06"),
        (1e-05, "1e-05"),
        (1.0000001e-16, "1e-16"),
        (5.0, "5"),
        (1234.5678, "1234.57"),
        (2.5, "2.5"),
        (1e100, "1e+100"),
        (1e-100, "1e-100"),
        (7.19024, "7.19024"),
        (2.34107, "2.34107"),
        (10.0, "10"),
        (100000.0, "100000"),
        (1000000.0, "1e+06"),
    ] {
        assert_eq!(format_g(v, 6), want, "%.6g of {v:?}");
    }
}
