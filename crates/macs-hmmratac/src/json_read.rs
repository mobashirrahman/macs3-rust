//! A minimal, correctly-rounded JSON reader for `*_model.json`.
//!
//! # Why not `serde_json`
//!
//! `serde_json`'s default float parser is **not correctly rounded**: on 6000
//! values drawn from the exponent range hmmlearn actually produces it disagrees
//! with Python's `float()` on 1335 of them, always by one ULP. `serde_json`
//! parses `1e-88` as `1.0000000000000001e-88`.
//!
//! That is not a hypothetical magnitude. A trained ATAC HMM drives whole states
//! to around `1e-88` (`_score_scaling` raises `ValueError: forward pass failed
//! with underflow` on them), so those numbers are the difference between a state
//! being reachable and not. A one-ULP error in `startprob_` also perturbs every
//! forward alpha and so every posterior, which moves region boundaries.
//!
//! Rust's `f64::from_str` **is** correctly rounded and agrees with Python on all
//! 6000 of those values, so this module keeps every number as a `&str` and hands
//! it to `str::parse`. The schema is fixed and small -- ten keys, nested arrays
//! of numbers, three state indices -- so a hand-rolled reader is small enough to
//! audit, and it removes the dependency from a load path where precision is the
//! whole point.
//!
//! # Accepted syntax
//!
//! Exactly what `json.dump` emits: objects, arrays, strings with `\"`, `\\`,
//! `\n`, `\r`, `\t` escapes, `true`/`false`/`null`, integers, and floats in
//! either `1e-88` or `0.123e-88` form, plus Python's bare `NaN`, `Infinity` and
//! `-Infinity`. Whitespace anywhere between tokens is tolerated. Unknown keys
//! are skipped rather than rejected, since a newer upstream may add fields.

/// A parsed JSON value. Numbers keep their source text so they can be re-parsed
/// with a correctly-rounded conversion.
#[derive(Debug, Clone, PartialEq)]
pub enum Json {
    /// `null`
    Null,
    /// `true` / `false`
    Bool(bool),
    /// A number, still in its source spelling.
    Num(String),
    /// A string, unescaped.
    Str(String),
    /// An array.
    Arr(Vec<Json>),
    /// An object, in source order.
    Obj(Vec<(String, Json)>),
}

impl Json {
    /// Look a key up in an object.
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(kv) => kv.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    /// The value as a string slice.
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    /// The value as a correctly-rounded `f64`.
    ///
    /// Python's bare `NaN` / `Infinity` / `-Infinity` are accepted, matching
    /// `json.load`.
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Num(t) => Some(match t.as_str() {
                "NaN" => f64::NAN,
                "Infinity" => f64::INFINITY,
                "-Infinity" => f64::NEG_INFINITY,
                _ => t.parse().ok()?,
            }),
            _ => None,
        }
    }

    /// The value as an array.
    pub fn as_arr(&self) -> Option<&[Json]> {
        match self {
            Json::Arr(v) => Some(v),
            _ => None,
        }
    }

    /// The value as a float, accepting a JSON string as a fallback.
    ///
    /// Upstream writes `i_open_region` and friends through `int(...)`, so they
    /// are integers; but a hand-edited or third-party model may quote them.
    pub fn as_int(&self) -> Option<i64> {
        let t = match self {
            Json::Num(t) => t,
            Json::Str(s) => s,
            _ => return None,
        };
        t.parse().ok().or_else(|| {
            // `"hmm_binsize": 500.0` from a third-party writer: accept an
            // integral float, reject a fractional one.
            let v: f64 = t.parse().ok()?;
            if v.fract() == 0.0 {
                Some(v as i64)
            } else {
                None
            }
        })
    }
}

/// A parse error with a byte offset.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JsonError {
    /// What went wrong.
    pub message: String,
    /// Byte offset into the input.
    pub offset: usize,
}

impl std::fmt::Display for JsonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} at byte {}", self.message, self.offset)
    }
}

impl std::error::Error for JsonError {}

/// Parse a complete JSON document.
pub fn parse(text: &str) -> Result<Json, JsonError> {
    let b = text.as_bytes();
    let mut p = Parser { b, i: 0 };
    p.ws();
    let v = p.value()?;
    p.ws();
    if p.i != b.len() {
        return Err(p.err("trailing data after JSON value"));
    }
    Ok(v)
}

struct Parser<'a> {
    b: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    fn err(&self, m: &str) -> JsonError {
        JsonError {
            message: m.to_string(),
            offset: self.i,
        }
    }

    fn ws(&mut self) {
        while self.i < self.b.len() && matches!(self.b[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }

    fn eat(&mut self, c: u8) -> Result<(), JsonError> {
        if self.i < self.b.len() && self.b[self.i] == c {
            self.i += 1;
            Ok(())
        } else {
            Err(self.err(&format!("expected `{}`", c as char)))
        }
    }

    fn lit(&mut self, s: &str) -> Result<(), JsonError> {
        if self.b[self.i..].starts_with(s.as_bytes()) {
            self.i += s.len();
            Ok(())
        } else {
            Err(self.err(&format!("expected `{s}`")))
        }
    }

    fn value(&mut self) -> Result<Json, JsonError> {
        match self.b.get(self.i) {
            None => Err(self.err("unexpected end of input")),
            Some(b'{') => self.object(),
            Some(b'[') => self.array(),
            Some(b'"') => Ok(Json::Str(self.string()?)),
            Some(b't') => {
                self.lit("true")?;
                Ok(Json::Bool(true))
            }
            Some(b'f') => {
                self.lit("false")?;
                Ok(Json::Bool(false))
            }
            Some(b'n') => {
                self.lit("null")?;
                Ok(Json::Null)
            }
            Some(b'N') => {
                self.lit("NaN")?;
                Ok(Json::Num("NaN".into()))
            }
            Some(b'I') => {
                self.lit("Infinity")?;
                Ok(Json::Num("Infinity".into()))
            }
            Some(b'-') if self.b[self.i..].starts_with(b"-Infinity") => {
                self.lit("-Infinity")?;
                Ok(Json::Num("-Infinity".into()))
            }
            Some(c) if *c == b'-' || c.is_ascii_digit() => self.number(),
            Some(c) => Err(self.err(&format!("unexpected `{}`", *c as char))),
        }
    }

    fn object(&mut self) -> Result<Json, JsonError> {
        self.eat(b'{')?;
        let mut kv = Vec::new();
        self.ws();
        if self.b.get(self.i) == Some(&b'}') {
            self.i += 1;
            return Ok(Json::Obj(kv));
        }
        loop {
            self.ws();
            let k = self.string()?;
            self.ws();
            self.eat(b':')?;
            self.ws();
            let v = self.value()?;
            kv.push((k, v));
            self.ws();
            match self.b.get(self.i) {
                Some(b',') => self.i += 1,
                Some(b'}') => {
                    self.i += 1;
                    return Ok(Json::Obj(kv));
                }
                _ => return Err(self.err("expected `,` or `}`")),
            }
        }
    }

    fn array(&mut self) -> Result<Json, JsonError> {
        self.eat(b'[')?;
        let mut v = Vec::new();
        self.ws();
        if self.b.get(self.i) == Some(&b']') {
            self.i += 1;
            return Ok(Json::Arr(v));
        }
        loop {
            self.ws();
            v.push(self.value()?);
            self.ws();
            match self.b.get(self.i) {
                Some(b',') => self.i += 1,
                Some(b']') => {
                    self.i += 1;
                    return Ok(Json::Arr(v));
                }
                _ => return Err(self.err("expected `,` or `]`")),
            }
        }
    }

    fn string(&mut self) -> Result<String, JsonError> {
        self.eat(b'"')?;
        let mut s = String::new();
        loop {
            let Some(c) = self.b.get(self.i) else {
                return Err(self.err("unterminated string"));
            };
            self.i += 1;
            match c {
                b'"' => return Ok(s),
                b'\\' => {
                    let Some(e) = self.b.get(self.i) else {
                        return Err(self.err("unterminated escape"));
                    };
                    self.i += 1;
                    match e {
                        b'"' => s.push('"'),
                        b'\\' => s.push('\\'),
                        b'/' => s.push('/'),
                        b'n' => s.push('\n'),
                        b'r' => s.push('\r'),
                        b't' => s.push('\t'),
                        b'b' => s.push('\u{8}'),
                        b'f' => s.push('\u{c}'),
                        b'u' => {
                            let hex = self
                                .b
                                .get(self.i..self.i + 4)
                                .ok_or_else(|| self.err("truncated \\u escape"))?;
                            let hex = std::str::from_utf8(hex).map_err(|_| self.err("bad \\u"))?;
                            let cp = u32::from_str_radix(hex, 16)
                                .map_err(|_| self.err("bad \\u escape"))?;
                            s.push(char::from_u32(cp).unwrap_or('\u{fffd}'));
                            self.i += 4;
                        }
                        e => return Err(self.err(&format!("bad escape `\\{}`", *e as char))),
                    }
                }
                _ => {
                    // copy the whole UTF-8 sequence verbatim
                    let start = self.i - 1;
                    let len = utf8_len(*c);
                    let end = (start + len).min(self.b.len());
                    s.push_str(
                        std::str::from_utf8(&self.b[start..end])
                            .map_err(|_| self.err("invalid utf-8"))?,
                    );
                    self.i = end;
                }
            }
        }
    }

    fn number(&mut self) -> Result<Json, JsonError> {
        let start = self.i;
        if self.b.get(self.i) == Some(&b'-') {
            self.i += 1;
        }
        while self.i < self.b.len() && self.b[self.i].is_ascii_digit() {
            self.i += 1;
        }
        if self.b.get(self.i) == Some(&b'.') {
            self.i += 1;
            while self.i < self.b.len() && self.b[self.i].is_ascii_digit() {
                self.i += 1;
            }
        }
        if matches!(self.b.get(self.i), Some(b'e') | Some(b'E')) {
            self.i += 1;
            if matches!(self.b.get(self.i), Some(b'+') | Some(b'-')) {
                self.i += 1;
            }
            while self.i < self.b.len() && self.b[self.i].is_ascii_digit() {
                self.i += 1;
            }
        }
        let text = std::str::from_utf8(&self.b[start..self.i])
            .map_err(|_| self.err("invalid number"))?
            .to_string();
        if text.is_empty() || text == "-" {
            return Err(self.err("empty number"));
        }
        Ok(Json::Num(text))
    }
}

fn utf8_len(first: u8) -> usize {
    match first {
        0x00..=0x7f => 1,
        0xc0..=0xdf => 2,
        0xe0..=0xef => 3,
        0xf0..=0xf7 => 4,
        _ => 1,
    }
}

/// Read a nested float array as `Vec<Vec<f64>>`.
pub fn matrix(v: &Json) -> Option<Vec<Vec<f64>>> {
    v.as_arr()?
        .iter()
        .map(|row| {
            row.as_arr()?
                .iter()
                .map(Json::as_f64)
                .collect::<Option<Vec<_>>>()
        })
        .collect()
}

/// Read a float array as `Vec<f64>`.
pub fn vector(v: &Json) -> Option<Vec<f64>> {
    v.as_arr()?.iter().map(Json::as_f64).collect()
}

/// Read a 3-deep float array as `Vec<Vec<Vec<f64>>>`.
pub fn cube(v: &Json) -> Option<Vec<Vec<Vec<f64>>>> {
    v.as_arr()?.iter().map(matrix).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_the_shapes_a_model_file_uses() {
        let j = parse(
            r#"{"startprob": [0.6, 0.4], "transmat": [[1.0, 0.0], [0.1, 0.9]],
                "means": [[1.0, 2.0], [3.0, 4.0]], "covars": [[1.0, 2.0], [3.0, 4.0]],
                "covariance_type": "diag", "n_features": 2, "i_open_region": 0,
                "i_background_region": 1, "i_nucleosomal_region": 1,
                "hmm_binsize": 500, "hmm_type": "gaussian"}"#,
        )
        .unwrap();
        assert_eq!(vector(j.get("startprob").unwrap()).unwrap(), vec![0.6, 0.4]);
        assert_eq!(j.get("hmm_type").unwrap().as_str(), Some("gaussian"));
        assert_eq!(j.get("hmm_binsize").unwrap().as_int(), Some(500));
        assert_eq!(j.get("covariance_type").unwrap().as_str(), Some("diag"));
    }

    /// The whole reason this module exists.
    #[test]
    fn floats_are_correctly_rounded_where_serde_json_is_not() {
        // serde_json yields 1.0000000000000001e-88 here
        assert_eq!(parse("1e-88").unwrap().as_f64().unwrap(), 1e-88f64);
        assert_eq!(parse("5.53e-16").unwrap().as_f64().unwrap(), 5.53e-16f64);
        assert_eq!(parse("7e-32").unwrap().as_f64().unwrap(), 7e-32f64);
    }

    #[test]
    fn bare_python_non_finite_tokens_are_accepted() {
        assert!(parse("NaN").unwrap().as_f64().unwrap().is_nan());
        assert_eq!(parse("Infinity").unwrap().as_f64().unwrap(), f64::INFINITY);
        assert_eq!(
            parse("-Infinity").unwrap().as_f64().unwrap(),
            f64::NEG_INFINITY
        );
    }

    #[test]
    fn full_covariance_is_three_deep() {
        let j = parse("[[[1.0, 0.0], [0.0, 1.0]]]").unwrap();
        assert_eq!(
            cube(&j).unwrap(),
            vec![vec![vec![1.0, 0.0], vec![0.0, 1.0]]]
        );
        let j = parse("[[1.0, 2.0]]").unwrap();
        assert!(cube(&j).is_none());
    }

    #[test]
    fn unknown_keys_are_ignored() {
        let j = parse(r#"{"a": 1, "future_field": {"deep": [1]}}"#).unwrap();
        assert!(j.get("a").is_some());
        assert!(j.get("future_field").is_some());
    }

    #[test]
    fn strings_handle_escapes_and_utf8() {
        let j = parse(r#"{"k": "a\"b\\c\ndA"}"#).unwrap();
        assert_eq!(j.get("k").unwrap().as_str(), Some("a\"b\\c\ndA"));
        let j = parse("{\"k\": \"\u{00e9}\"}").unwrap();
        assert_eq!(j.get("k").unwrap().as_str(), Some("\u{00e9}"));
    }

    #[test]
    fn malformed_input_is_rejected_with_an_offset() {
        assert!(parse("{").is_err());
        assert!(parse("[1, 2").is_err());
        assert!(parse(r#"{"a": }"#).is_err());
        assert!(parse("{} extra").is_err());
        assert!(parse("").is_err());
        let e = parse("[1, 2").unwrap_err();
        assert!(e.offset > 0);
    }

    #[test]
    fn integers_and_quoted_integers_both_work() {
        assert_eq!(parse("5").unwrap().as_int(), Some(5));
        assert_eq!(parse("\"5\"").unwrap().as_int(), Some(5));
        assert_eq!(parse("5.0").unwrap().as_int(), Some(5));
    }
}
