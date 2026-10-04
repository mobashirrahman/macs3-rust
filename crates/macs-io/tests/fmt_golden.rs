//! `%.5g` transcription, checked against strings produced by CPython's own
//! `"%.5g" %` on the same values.
//!
//! The XLS acceptance criterion is byte-identical output, so the score columns
//! are *formatted*, not printed; getting this wrong is silent and invisible to a
//! numeric comparison. The expected strings here were produced by CPython.

use macs_io::peakout::format_g;

const CASES: &[(f64, &str)] = &[
    (0.0, "0"),
    (1.0, "1"),
    (9.0, "9"),
    (2.9371, "2.9371"),
    (3.09691, "3.0969"),
    (3.32607, "3.3261"),
    (3.41685, "3.4169"),
    (10.828, "10.828"),
    (10.8489, "10.849"),
    (11.2776, "11.278"),
    (13.341, "13.341"),
    (10.0, "10"),
    (11.0, "11"),
];

/// F119: upstream writes the score columns with **`%.6g`**, not `%.5g`
/// (`MACS3/IO/PeakIO.py:815-840`). These are the values taken from real golden
/// XLS files, with CPython's own `"%.6g" %` rendering as the expectation.
const SIX_G: &[(f64, &str)] = &[
    (7.20635, "7.20635"),
    (7.9977, "7.9977"),
    (1.23882, "1.23882"),
    (6.07218, "6.07218"),
    (6.45005, "6.45005"),
    (7.14137, "7.14137"),
    (5.23886, "5.23886"),
    (5.25981, "5.25981"),
];

#[test]
fn percent_g_matches_cpython_for_every_golden_value() {
    for (v, expect) in CASES {
        assert_eq!(&format_g(*v, 5), expect, "%.5g value {v}");
    }
}

#[test]
fn six_significant_digits_is_what_upstream_actually_writes() {
    for (v, expect) in SIX_G {
        assert_eq!(&format_g(*v, 6), expect, "%.6g value {v}");
    }
}

#[test]
fn percent_g_switches_to_exponent_outside_pythons_range() {
    // Python uses exponent form when exp < -4 or exp >= precision
    assert_eq!(format_g(0.000012345678, 5), "1.2346e-05");
    assert_eq!(format_g(123456.789, 5), "1.2346e+05");
    // and fixed form in between
    assert_eq!(format_g(0.000123456, 5), "0.00012346");
    assert_eq!(format_g(99999.0, 5), "99999");
}
