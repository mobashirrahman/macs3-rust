//! C `atoi` / `atof` semantics, transcribed.
//!
//! Upstream calls libc's `atoi` and `atof` directly (`Parser.py:26`,
//! `BedGraphIO.py:28`). Those functions are specified by nothing in particular and
//! implemented by everything slightly differently, so this module fixes the
//! behaviour that matters for MACS3 input:
//!
//! * leading whitespace is skipped
//! * an optional `+`/`-` sign is honoured
//! * digits are consumed until the first non-digit
//! * **trailing garbage is ignored**, not an error
//! * an unparseable value is `0`, not an error
//!
//! Saturation at `i32::MIN`/`i32::MAX` matches glibc, which is what upstream
//! runs against. A coordinate that overflows therefore wraps to a sentinel
//! extreme rather than being rejected, and for a *negative* one that sentinel is
//! then discarded by `build_fwtrack`'s `if fpos < 0` check.

/// Parse an integer exactly as C's `atoi` does.
///
/// ```
/// assert_eq!(macs_io::atoi(b"  -42abc"), -42);
/// assert_eq!(macs_io::atoi(b"nope"), 0);
/// assert_eq!(macs_io::atoi(b""), 0);
/// ```
pub fn atoi(s: &[u8]) -> i64 {
    let mut i = 0usize;
    while i < s.len() && is_c_space(s[i]) {
        i += 1;
    }
    let neg = match s.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let mut acc: i64 = 0;
    let mut overflow = false;
    while i < s.len() && s[i].is_ascii_digit() {
        acc = acc.saturating_mul(10).saturating_add((s[i] - b'0') as i64);
        if acc > i32::MAX as i64 + 1 {
            overflow = true;
            acc = i32::MAX as i64 + 1;
        }
        i += 1;
    }
    if !overflow && acc == 0 {
        return 0;
    }
    if neg {
        acc.saturating_neg()
    } else {
        acc.min(i32::MAX as i64)
    }
}

/// Parse a float exactly as C's `atof` does, to the extent that matters here.
///
/// `atof` is a `strtod` that yields `0.0` when no conversion happens, keeps the
/// sign of `-0`, and saturates to infinity on overflow. `NaN` text is accepted,
/// because `strtod` accepts it.
pub fn atof(s: &[u8]) -> f64 {
    let mut i = 0usize;
    while i < s.len() && is_c_space(s[i]) {
        i += 1;
    }
    if i >= s.len() {
        return 0.0;
    }
    // `strtod` also accepts INF and NAN, case-insensitively and with an optional
    // sign. Rust's `parse` accepts them too, but only as whole-string input,
    // whereas `strtod` stops at the first byte it cannot use -- so these are
    // handled before the prefix parse below.
    let rest = &s[i..];
    let (signed_neg, body) = match rest.first() {
        Some(b'-') => (true, &rest[1..]),
        Some(b'+') => (false, &rest[1..]),
        _ => (false, rest),
    };
    let lower: Vec<u8> = body.iter().map(|b| b.to_ascii_lowercase()).collect();
    if lower.starts_with(b"nan") {
        return f64::NAN;
    }
    if lower.starts_with(b"inf") {
        return if signed_neg {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    let text = String::from_utf8_lossy(&s[i..]);
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return 0.0;
    }
    parse_c_float_prefix(trimmed).unwrap_or(0.0)
}

/// `strtod`-style longest-prefix parse; `None` if no digits start a number.
fn parse_c_float_prefix(s: &str) -> Option<f64> {
    let b = s.as_bytes();
    let mut i = 0usize;
    if b[i] == b'+' || b[i] == b'-' {
        i += 1;
    }
    let int_start = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    let had_int = i > int_start;
    let mut had_frac = false;
    if i < b.len() && b[i] == b'.' {
        i += 1;
        let frac_start = i;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        had_frac = i > frac_start;
    }
    if !had_int && !had_frac {
        return None;
    }
    let mut end = i;
    // optional exponent, only consumed if it is complete
    if end < b.len() && (b[end] == b'e' || b[end] == b'E') {
        let mut j = end + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let ds = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > ds {
            end = j;
        }
    }
    s[..end].parse::<f64>().ok()
}

/// `isspace` in the C locale.
#[inline]
fn is_c_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_integers() {
        assert_eq!(atoi(b"0"), 0);
        assert_eq!(atoi(b"7"), 7);
        assert_eq!(atoi(b"-7"), -7);
        assert_eq!(atoi(b"+7"), 7);
    }

    #[test]
    fn leading_whitespace_and_trailing_garbage_are_tolerated() {
        assert_eq!(atoi(b"  \t42xyz"), 42);
        assert_eq!(atoi(b"-3.9"), -3);
        assert_eq!(atoi(b"1000000base"), 1000000);
    }

    #[test]
    fn unparseable_is_zero_not_an_error() {
        assert_eq!(atoi(b""), 0);
        assert_eq!(atoi(b"abc"), 0);
        assert_eq!(atoi(b"-"), 0);
        assert_eq!(atoi(b"   "), 0);
    }

    #[test]
    fn zero_with_a_sign_stays_zero() {
        assert_eq!(atoi(b"-0"), 0);
        assert_eq!(atoi(b"+0"), 0);
    }

    #[test]
    fn overflow_saturates_like_glibc_rather_than_wrapping() {
        assert_eq!(atoi(b"99999999999999999999"), i32::MAX as i64);
        assert_eq!(atoi(b"-99999999999999999999"), i32::MIN as i64);
        assert_eq!(atoi(b"2147483648"), i32::MAX as i64);
    }

    #[test]
    fn negative_coordinates_survive_so_the_caller_can_drop_them() {
        // upstream keeps negative positions through fw_parse_line and drops them
        // in build_fwtrack, so the parse must not clamp
        assert_eq!(atoi(b"-1"), -1);
    }

    #[test]
    fn floats() {
        assert_eq!(atof(b"1.5"), 1.5);
        assert_eq!(atof(b"  -2.25xyz"), -2.25);
        assert_eq!(atof(b"3e2"), 300.0);
        assert_eq!(atof(b"3e"), 3.0);
        assert_eq!(atof(b"1e-3"), 0.001);
        assert_eq!(atof(b".5"), 0.5);
        assert_eq!(atof(b"5."), 5.0);
    }

    #[test]
    fn unparseable_floats_are_zero() {
        assert_eq!(atof(b""), 0.0);
        assert_eq!(atof(b"abc"), 0.0);
        // NaN never compares equal, so this can only be checked with is_nan
        assert!(atof(b"nan").is_nan());
        assert_eq!(atof(b"inf"), f64::INFINITY);
        assert_eq!(atof(b"-Infinity"), f64::NEG_INFINITY);
    }
}
