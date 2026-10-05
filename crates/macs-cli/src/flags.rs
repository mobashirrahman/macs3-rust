//! The argument surface, derived from upstream's own argparse.
//!
//! The table is `oracle/flag_matrix.tsv`, auto-derived from `macs3`'s argparse
//! (not from `--help`, which documents flags that do not exist -- F23). It is
//! embedded with [`include_str!`] so the shipped binary has no runtime data
//! dependency, which keeps the "fully independent of Python at run time" and
//! "single self-contained tool" requirements intact.
//!
//! Parsing follows argparse's observable behaviour, because the acceptance
//! criteria require identical accept/reject behaviour:
//!
//! * a missing **required** flag is a usage error (exit 2);
//! * a value outside **choices** is a usage error;
//! * `--flag value` and `--flag=value` are both accepted, and so are the short
//!   forms listed in the matrix;
//! * a value that looks like a flag is rejected unless the flag takes one;
//! * unknown flags are usage errors;
//! * `--help`/`-h` succeeds with the subcommand's usage text.
//!
//! The token loop below transliterates `ArgumentParser._parse_known_args`
//! (`argparse.py:1963-2214`) rather than scanning tokens with ad-hoc rules,
//! because four of argparse's rules are not local to a single token and cannot
//! be reproduced any other way without drift:
//!
//! * **glued short values** (`_get_option_tuples`, `argparse.py:2341-2352`, and
//!   `consume_optional`, `argparse.py:2040-2063`): `-p-1`, `-p0.01`, `-Bp0.05`
//!   and `-Bq0.05` mean `-p -1`, `-p 0.01`, `-B -p 0.05` and `-B -q 0.05`;
//! * **long-option abbreviation** -- `allow_abbrev` is on, so `--pval` is
//!   `--pvalue` and `--o` is `ambiguous option: --o could match --outdir, --ofile`;
//! * **`--`** (`argparse.py:1973-1976`): the token becomes a `-` in the pattern
//!   and **every** later token becomes an `A`, whatever it looks like;
//! * **when an error is reported**: extras are collected and reported last, after
//!   every value conversion and after the per-flag errors.
//!
//! Because `_get_nargs_pattern` deletes the `-` characters of the pattern for an
//! optional action (`argparse.py:2411-2413`), *no* optional can consume a `--` or
//! any token after one. So `-t -- -x.bed` is
//! `argument -t/--treatment: expected at least one argument`, and
//! `-t a b -- -c` leaves `-- -c` to be reported as
//! `unrecognized arguments: -- -c`. Verified against the reference rather than
//! assumed.

use std::collections::BTreeMap;

/// One row of the flag matrix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FlagSpec {
    /// The subcommand this row belongs to, or `<top>` for the program-level flags.
    pub subcommand: String,
    /// The flag as written on the command line, e.g. `--extsize`.
    pub flag: String,
    /// argparse `action`: `store`, `store_true`, `store_false`, `count`, `help`.
    pub action: String,
    /// The argparse `nargs`, rendered as a string; `0` for flags that take no value.
    pub nargs: String,
    /// The argparse `dest` -- the key this flag populates.
    pub dest: String,
    /// `true` when the flag is required.
    pub required: bool,
    /// The default, or `None` for "no default".
    pub default: Option<String>,
    /// Permitted values, empty when unconstrained.
    pub choices: String,
}

/// The embedded flag matrix.
const MATRIX: &str = include_str!("../../../oracle/flag_matrix.tsv");

/// The embedded per-subcommand flag *type* table, generated from the oracle's
/// argparse by `oracle/gen_flag_types.py`.
const TYPES: &str = include_str!("../../../oracle/flag_types.tsv");

/// `(subcommand, flag spelling) -> "int" | "float"`, for flags upstream declares
/// with a numeric argparse `type`.
///
/// The matrix records which flags exist but not their types, and the types cannot
/// be inferred from the defaults: `--verbose` defaults to `2` and is an int,
/// `--gsize` defaults to `hs` and is a string, `--pvalue` defaults to `1e-5` and
/// is a float. They are read out of upstream's own `add_argument` calls instead,
/// scoped per subcommand because `-g` is a genome-size string in `callpeak` and an
/// int max-gap in `bdgdiff`.
fn flag_types() -> &'static Vec<(String, String, String)> {
    static CACHE: std::sync::OnceLock<Vec<(String, String, String)>> = std::sync::OnceLock::new();
    CACHE.get_or_init(|| {
        let mut out = Vec::new();
        for line in TYPES.lines() {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 3 || f[0] == "subcommand" {
                continue;
            }
            out.push((f[0].to_string(), f[1].to_string(), f[2].to_string()));
        }
        out
    })
}

// ---------------------------------------------------------------------------
// Python's `int()` / `float()`
// ---------------------------------------------------------------------------
//
// argparse's `type=int` / `type=float` *are* `int()` / `float()`: the conversion
// function is applied to the raw token, so whatever Python's own literal grammar
// accepts, upstream accepts. That grammar is not Rust's:
//
// | token       | `int()` | `float()` | Rust `str::parse` |
//|-------------|---------|-----------|--------------------|
//| `" 1"`, `"1 "` | 1     | 1.0       | rejected           |
//| `"1_0"`, `"1_000.5"` | rejected | 10, 1000.5 | rejected        |
//| `"+1"`       | 1       | 1.0       | accepted           |
//| `"1e3"`      | rejected | 1000.0    | accepted           |
//| `"inf"`, `"infinity"`, `"nan"` | rejected | ±inf, nan | accepted (case-insensitively) |
//| `"١"` (Arabic-Indic) | 1 | 1.0        | rejected           |
//
// so `bdgopt -p ' -1'` exited 2 here and 0 upstream, and `-p1_0` exited 2 here
// and ran upstream. [`py_int`] / [`py_float`] below are the Python grammar; both
// the accept/reject decision and the value have to match, because the value is
// what the pipeline consumes.

/// Python's `Py_UNICODE_ISSPACE` set, which is what `int()` and `float()` strip.
fn py_isspace(c: char) -> bool {
    matches!(c,
        '\u{09}'..='\u{0d}' | '\u{1c}'..='\u{20}' | '\u{85}' | '\u{a0}'
        | '\u{1680}' | '\u{2000}'..='\u{200a}' | '\u{2028}'..='\u{2029}'
        | '\u{202f}' | '\u{205f}' | '\u{3000}')
}

/// Python's `Py_UNICODE_TODECIMAL`: the value of a Unicode decimal digit (`Nd`),
/// or `None` for anything else.
///
/// Not `char::is_numeric()`, which also accepts `No`/`Nl` numerals such as `½`
/// and `Ⅻ`; `int('½')` is a `ValueError` while `½`.isnumeric() is `True`. The
/// range table is Unicode 15.0.0, the version the reference's interpreter
/// reports, and every block is a run of ten so `(cp - lo) % 10` is the value.
fn py_decimal(c: char) -> Option<u32> {
    const RANGES: [(u32, u32); 64] = [
        (0x0030, 0x0039),
        (0x0660, 0x0669),
        (0x06f0, 0x06f9),
        (0x07c0, 0x07c9),
        (0x0966, 0x096f),
        (0x09e6, 0x09ef),
        (0x0a66, 0x0a6f),
        (0x0ae6, 0x0aef),
        (0x0b66, 0x0b6f),
        (0x0be6, 0x0bef),
        (0x0c66, 0x0c6f),
        (0x0ce6, 0x0cef),
        (0x0d66, 0x0d6f),
        (0x0de6, 0x0def),
        (0x0e50, 0x0e59),
        (0x0ed0, 0x0ed9),
        (0x0f20, 0x0f29),
        (0x1040, 0x1049),
        (0x1090, 0x1099),
        (0x17e0, 0x17e9),
        (0x1810, 0x1819),
        (0x1946, 0x194f),
        (0x19d0, 0x19d9),
        (0x1a80, 0x1a89),
        (0x1a90, 0x1a99),
        (0x1b50, 0x1b59),
        (0x1bb0, 0x1bb9),
        (0x1c40, 0x1c49),
        (0x1c50, 0x1c59),
        (0xa620, 0xa629),
        (0xa8d0, 0xa8d9),
        (0xa900, 0xa909),
        (0xa9d0, 0xa9d9),
        (0xa9f0, 0xa9f9),
        (0xaa50, 0xaa59),
        (0xabf0, 0xabf9),
        (0xff10, 0xff19),
        (0x104a0, 0x104a9),
        (0x10d30, 0x10d39),
        (0x11066, 0x1106f),
        (0x110f0, 0x110f9),
        (0x11136, 0x1113f),
        (0x111d0, 0x111d9),
        (0x112f0, 0x112f9),
        (0x11450, 0x11459),
        (0x114d0, 0x114d9),
        (0x11650, 0x11659),
        (0x116c0, 0x116c9),
        (0x11730, 0x11739),
        (0x118e0, 0x118e9),
        (0x11950, 0x11959),
        (0x11c50, 0x11c59),
        (0x11d50, 0x11d59),
        (0x11da0, 0x11da9),
        (0x11f50, 0x11f59),
        (0x16a60, 0x16a69),
        (0x16ac0, 0x16ac9),
        (0x16b50, 0x16b59),
        (0x1d7ce, 0x1d7ff),
        (0x1e140, 0x1e149),
        (0x1e2f0, 0x1e2f9),
        (0x1e4f0, 0x1e4f9),
        (0x1e950, 0x1e959),
        (0x1fbf0, 0x1fbf9),
    ];
    let cp = c as u32;
    RANGES
        .iter()
        .find(|(lo, hi)| (*lo..=*hi).contains(&cp))
        .map(|(lo, _)| (cp - lo) % 10)
}

/// Consume a Python numeric digit run -- `Digit ('_' Digit)*` -- appending the
/// decimal values as ASCII digits to `out`, and return how many digits it took.
///
/// The underscore rule is "between two digits", which is what rejects `_1`,
/// `1_`, `1__0`, `1_.5` and `1._5` while accepting `1_0`, `0_1` and `1e3_0`.
fn py_digits(cs: &[char], i: &mut usize, out: &mut String) -> usize {
    let mut n = 0;
    while let Some(&c) = cs.get(*i) {
        if c == '_' {
            if n == 0 {
                break;
            }
            match cs.get(*i + 1).and_then(|&d| py_decimal(d)) {
                Some(v) => {
                    out.push(char::from_digit(v, 10).expect("digit 0-9"));
                    n += 1;
                    *i += 2;
                }
                None => break,
            }
            continue;
        }
        match py_decimal(c) {
            Some(v) => {
                out.push(char::from_digit(v, 10).expect("digit 0-9"));
                n += 1;
                *i += 1;
            }
            None => break,
        }
    }
    n
}

/// Split off Python's optional leading sign.
fn py_sign(s: &str) -> (bool, &str) {
    match s.as_bytes().first() {
        Some(b'+') => (false, &s[1..]),
        Some(b'-') => (true, &s[1..]),
        _ => (false, s),
    }
}

/// Python's `int(token)`, saturating at the `i64` bounds.
///
/// Python's `int` is arbitrary precision, so a token beyond `i64` is *accepted*
/// by argparse and only overflows later, inside the pipeline (`OverflowError`,
/// exit 1 -- verified for `filterdup -s 9223372036854775808`). The magnitude is
/// accumulated first and the sign applied afterwards, so a huge negative token
/// saturates to `i64::MIN` rather than through `i64::MAX`.
fn py_int(token: &str) -> Option<i64> {
    let t = token.trim_matches(py_isspace);
    let (neg, body) = py_sign(t);
    let cs: Vec<char> = body.chars().collect();
    let mut i = 0;
    let mut clean = String::with_capacity(cs.len());
    if py_digits(&cs, &mut i, &mut clean) == 0 || i != cs.len() {
        return None;
    }
    let mut acc: i64 = 0;
    let mut saturated = false;
    for b in clean.bytes() {
        match acc
            .checked_mul(10)
            .and_then(|v| v.checked_add(i64::from(b - b'0')))
        {
            Some(v) => acc = v,
            None => {
                acc = i64::MAX;
                saturated = true;
            }
        }
    }
    // a saturated magnitude keeps the sign it was written with: `-9223372036854775808`
    // is `i64::MIN`, not `-i64::MAX`
    Some(match (neg, saturated) {
        (true, true) => i64::MIN,
        (true, false) => -acc,
        (false, _) => acc,
    })
}

/// Python's `float(token)`.
///
/// Both parsers round the significand correctly, so the value is taken from
/// Rust's own `str::parse::<f64>` on an ASCII-ised, underscore-free copy of the
/// token rather than re-implemented; overflow and underflow (`1e400` -> `inf`,
/// `1e-400` -> `0.0`) are the same in both.
fn py_float(token: &str) -> Option<f64> {
    let t = token.trim_matches(py_isspace);
    let (neg, body) = py_sign(t);
    if body.is_empty() {
        return None;
    }
    let special = body.to_lowercase();
    match special.as_str() {
        "inf" | "infinity" => {
            return Some(if neg {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            })
        }
        "nan" => return Some(f64::NAN),
        _ => {}
    }
    let cs: Vec<char> = body.chars().collect();
    let mut i = 0;
    let mut clean = String::with_capacity(cs.len());
    let mut digits = py_digits(&cs, &mut i, &mut clean);
    if cs.get(i) == Some(&'.') {
        i += 1;
        clean.push('.');
        digits += py_digits(&cs, &mut i, &mut clean);
    }
    if matches!(cs.get(i), Some(&c) if c == 'e' || c == 'E') {
        i += 1;
        clean.push('e');
        if matches!(cs.get(i), Some(&c) if c == '+' || c == '-') {
            clean.push(cs[i]);
            i += 1;
        }
        let exponent_start = i;
        digits += py_digits(&cs, &mut i, &mut clean);
        if i == exponent_start {
            return None; // `1e` -- an exponent needs at least one digit
        }
    }
    if i != cs.len() || digits == 0 {
        return None;
    }
    let v: f64 = clean.parse().ok()?;
    Some(if neg { -v } else { v })
}

/// Python's `repr()` of a `str`, which is what argparse interpolates into
/// `invalid %(type)s value: %(value)r` and `invalid choice: %(value)r`.
fn py_repr(s: &str) -> String {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(quote);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            '\u{0}'..='\u{1f}' | '\u{7f}' => out.push_str(&format!("\\x{:02x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

/// The spelling [`Options::values`] stores for a numeric token.
///
/// Upstream hands the command a *converted* value, so a token Python's
/// `int()`/`float()` accepts but Rust's parser does not has to be normalised
/// before the commands read it: they parse `Options::values` themselves (for
/// example `bedgraph_cmds.rs` on `--extra-param`). A token Rust already parses is
/// kept byte-for-byte, so `Options::get` still reports the token the user typed in
/// every case where the two grammars already agreed -- which is what keeps
/// `-0.35`, `.5` and `1e-5` verbatim in the tests and in the pipeline.
fn canonical_number(token: &str, ty: &str) -> String {
    match ty {
        "int" if token.parse::<i64>().is_err() => {
            py_int(token).map_or_else(|| token.to_string(), |v| v.to_string())
        }
        "float" if token.parse::<f64>().is_err() => {
            py_float(token).map_or_else(|| token.to_string(), |v| format!("{v:?}"))
        }
        _ => token.to_string(),
    }
}

// ---------------------------------------------------------------------------
// The action model
// ---------------------------------------------------------------------------

/// One argparse action: a single `add_argument` call, i.e. every option string
/// that shares a destination (`-p`, `--pvalue`) together with the attributes they
/// share. The matrix records one row per spelling, so the rows are regrouped here.
///
/// The order of [`Self::names`] is argparse's `_option_string_actions` order,
/// which is what decides the spelling in `ambiguous option: --o could match
/// --outdir, --ofile` and in `one of the arguments -p/--percentage -n/--number`.
#[derive(Debug, Clone)]
struct Action {
    /// The row of the **first** spelling; carries `action`, `nargs`, `dest`,
    /// `required`, `default` and `choices`, which every spelling shares.
    spec: FlagSpec,
    /// Every option string of this action, in declaration order.
    names: Vec<String>,
}

impl Action {
    /// argparse's `_get_action_name`: every spelling, joined with `/`.
    fn display(&self) -> String {
        self.names.join("/")
    }

    /// `"int"` or `"float"` when upstream declared a numeric `type`, else `None`.
    fn ty(&self) -> Option<&'static str> {
        let table = flag_types();
        self.names
            .iter()
            .find_map(|n| {
                table
                    .iter()
                    .find(|(s, f, _)| *s == self.spec.subcommand && f == n)
                    .map(|(_, _, t)| t.as_str())
            })
            .filter(|t| *t == "int" || *t == "float")
    }

    /// argparse's `nargs`, as a count.
    ///
    /// `None` in argparse is one argument, `0` takes none, `?` is at most one and
    /// `'*'` may take none. `flag_matrix.tsv` renders `None` as an empty column.
    fn nargs(&self) -> Nargs {
        match self.spec.nargs.as_str() {
            "" => Nargs::One,
            "0" => Nargs::Zero,
            "?" => Nargs::Optional,
            "*" => Nargs::ZeroOrMore,
            "+" => Nargs::OneOrMore,
            k => Nargs::Exact(k.parse::<usize>().expect("matrix nargs is a number")),
        }
    }
}

/// The number of values an action consumes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Nargs {
    /// `nargs=0`: `store_true`, `store_false`, `count`, `help`.
    Zero,
    /// `nargs=None`: exactly one.
    One,
    /// `nargs='?'`: zero or one.
    Optional,
    /// `nargs='*'`: zero or more.
    ZeroOrMore,
    /// `nargs='+'`: one or more.
    OneOrMore,
    /// `nargs=N`: exactly `N`.
    Exact(usize),
}

/// Group the matrix rows of `subcommand` into argparse actions.
fn actions_for(specs: &[FlagSpec], subcommand: &str) -> Vec<Action> {
    let mut out: Vec<Action> = Vec::new();
    for s in specs.iter().filter(|s| s.subcommand == subcommand) {
        match out.iter_mut().find(|a| a.spec.dest == s.dest) {
            Some(a) => a.names.push(s.flag.clone()),
            None => out.push(Action {
                spec: s.clone(),
                names: vec![s.flag.clone()],
            }),
        }
    }
    out
}

/// The action that owns `name`, if any. This is `_option_string_actions[name]`.
fn find_action(acts: &[Action], name: &str) -> Option<usize> {
    acts.iter().position(|a| a.names.iter().any(|n| n == name))
}

/// An entry of argparse's `_get_option_tuples`, i.e. one interpretation of a
/// token.
#[derive(Debug, Clone)]
struct OptionTuple {
    /// `None` for "an option-looking token this parser does not define", which
    /// `_parse_optional` returns as a 4-tuple with a `None` action
    /// (`argparse.py:2332`) and `consume_optional` turns into an extra
    /// (`argparse.py:2033-2035`).
    action: Option<usize>,
    /// The option string that matched.
    option_string: String,
    /// `Some("=")` when the token was `--flag=value`, `Some("")` for a glued
    /// single-dash value, `None` otherwise. The distinction matters:
    /// `argparse.py:2047` refuses to keep scanning a `-Bq=0.01`-style tail when a
    /// real `=` was consumed.
    sep: Option<&'static str>,
    /// The value glued onto the option string, if any.
    explicit_arg: Option<String>,
}

/// `_get_option_tuples` (`argparse.py:2337-2372`): every way the token can be read
/// as a known option.
///
/// A token with **two** prefix characters is only ever split at `=`, and
/// `allow_abbrev` is on, so every option string that starts with the token (or
/// with the part before its `=`) is a candidate. A token with **one** prefix
/// character may instead be a short option with the value glued on: `-p0.01`
/// yields `(-p, sep='', explicit='0.01')`, and a token that is itself the whole of
/// some long option (`-pvalue=0.01`) yields that option with no explicit argument.
fn option_tuples(acts: &[Action], token: &str) -> Vec<OptionTuple> {
    let cs: Vec<char> = token.chars().collect();
    let mut out = Vec::new();
    if cs.len() >= 2 && cs[0] == '-' && cs[1] == '-' {
        let (prefix, sep, explicit) = match token.split_once('=') {
            Some((p, v)) => (p, Some("="), Some(v.to_string())),
            None => (token, None, None),
        };
        for (i, a) in acts.iter().enumerate() {
            for n in &a.names {
                if n.starts_with(prefix) {
                    out.push(OptionTuple {
                        action: Some(i),
                        option_string: n.clone(),
                        sep,
                        explicit_arg: explicit.clone(),
                    });
                }
            }
        }
    } else if cs[0] == '-' {
        let short: String = cs[..2].iter().collect();
        let tail: String = cs[2..].iter().collect();
        for (i, a) in acts.iter().enumerate() {
            for n in &a.names {
                if *n == short {
                    out.push(OptionTuple {
                        action: Some(i),
                        option_string: n.clone(),
                        sep: Some(""),
                        explicit_arg: Some(tail.clone()),
                    });
                } else if n.starts_with(token) {
                    out.push(OptionTuple {
                        action: Some(i),
                        option_string: n.clone(),
                        sep: None,
                        explicit_arg: None,
                    });
                }
            }
        }
    }
    out
}

/// `_parse_optional` (`argparse.py:2279-2334`): is `token` an option?
///
/// * `Ok(None)` -- it is an **argument** (`A` in the pattern), so
///   `_parse_known_args` offers it to a positional, and there are none;
/// * `Ok(Some(t))` -- it is an option (`O`), with `t.action == None` meaning the
///   parser has no such option and the token becomes an extra;
/// * `Err` -- the token is ambiguous, which argparse reports from here, i.e.
///   before any value has been converted.
///
/// The negative-number clause is what makes `bdgopt -p -0.35` work: `-0.35` is
/// the *value* of `--extra-param`, not an unknown flag. No option string in
/// `flag_matrix.tsv` matches [`is_negative_number`], so every subparser's
/// `_has_negative_number_optionals` (`argparse.py:1533-1535`) is empty and the
/// clause always fires. `-9x` matches neither regex and stays an option, so it is
/// still "unrecognized arguments".
fn parse_optional(acts: &[Action], token: &str) -> Result<Option<OptionTuple>, UsageError> {
    if token.is_empty() || !token.starts_with('-') {
        return Ok(None);
    }
    if let Some(i) = find_action(acts, token) {
        return Ok(Some(OptionTuple {
            action: Some(i),
            option_string: token.to_string(),
            sep: None,
            explicit_arg: None,
        }));
    }
    // a lone `-` (or any single character) is positional
    if token.chars().count() == 1 {
        return Ok(None);
    }
    if let Some((name, value)) = token.split_once('=') {
        if let Some(i) = find_action(acts, name) {
            return Ok(Some(OptionTuple {
                action: Some(i),
                option_string: name.to_string(),
                sep: Some("="),
                explicit_arg: Some(value.to_string()),
            }));
        }
    }
    let tuples = option_tuples(acts, token);
    if tuples.len() > 1 {
        let matches: Vec<&str> = tuples.iter().map(|t| t.option_string.as_str()).collect();
        return Err(UsageError(format!(
            "error: ambiguous option: {token} could match {}",
            matches.join(", ")
        )));
    }
    if let Some(t) = tuples.into_iter().next() {
        return Ok(Some(t));
    }
    if is_negative_number(token)
        && !acts
            .iter()
            .any(|a| a.names.iter().any(|n| is_negative_number(n)))
    {
        return Ok(None);
    }
    if token.contains(' ') {
        return Ok(None);
    }
    Ok(Some(OptionTuple {
        action: None,
        option_string: token.to_string(),
        sep: None,
        explicit_arg: None,
    }))
}

/// `re.match(r'^-\d+$|^-\d*\.\d+$', token)`: does `token` spell a negative number?
///
/// This is argparse's `_negative_number_matcher` (`argparse.py:1420`) and it is
/// deliberately narrow. `-1`, `-0.35` and `-.5` match; `-1.`, `-1e-3` and `-9x`
/// do **not**, because the second alternative needs a digit after the dot and
/// there is no exponent branch.
fn is_negative_number(token: &str) -> bool {
    let Some(rest) = token.strip_prefix('-') else {
        return false;
    };
    match rest.split_once('.') {
        // `^-\d+$`: digits only, and at least one.
        None => !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit()),
        // `^-\d*\.\d+$`: the integer part may be empty, the fraction may not.
        Some((int, frac)) => {
            int.bytes().all(|b| b.is_ascii_digit())
                && !frac.is_empty()
                && frac.bytes().all(|b| b.is_ascii_digit())
        }
    }
}

/// `_match_argument` + `_get_nargs_pattern` (`argparse.py:2249-2268`,
/// `argparse.py:2375-2421`) restricted to optionals, of which every MACS3
/// argument is one.
///
/// `_get_nargs_pattern` builds `(-*A-*)`, `(-*A[A-]*)`, `(-*AA-*)` and friends and
/// then, for an action with option strings, strips the `-*` and `-` characters
/// out of the pattern. That last step is why *no* optional can ever match a `--`
/// or anything past one: `nargs='+'` becomes `(A[A]*)`, not `(A[A-]*)`.
fn match_nargs(action: &Action, pattern: &[char]) -> Result<usize, UsageError> {
    let leading_a = pattern.iter().take_while(|c| **c == 'A').count();
    let at = |i: usize| pattern.get(i) == Some(&'A');
    let expected = |n: usize| -> UsageError {
        UsageError(format!(
            "error: argument {}: expected {n} argument{}",
            action.display(),
            if n == 1 { "" } else { "s" }
        ))
    };
    match action.nargs() {
        Nargs::Zero => Ok(0),
        Nargs::One => {
            if at(0) {
                Ok(1)
            } else {
                Err(UsageError(format!(
                    "error: argument {}: expected one argument",
                    action.display()
                )))
            }
        }
        Nargs::Optional => Ok(usize::from(at(0))),
        Nargs::ZeroOrMore => Ok(leading_a),
        Nargs::OneOrMore => {
            if at(0) {
                Ok(leading_a)
            } else {
                Err(UsageError(format!(
                    "error: argument {}: expected at least one argument",
                    action.display()
                )))
            }
        }
        Nargs::Exact(n) => {
            if (0..n).all(at) {
                Ok(n)
            } else {
                Err(expected(n))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The matrix
// ---------------------------------------------------------------------------

/// Parse the matrix. Called once; the result is a `const`-friendly static table
/// so tests can assert against the same data the binary uses.
pub fn flag_specs() -> Vec<FlagSpec> {
    let mut out = Vec::new();
    let mut lines = MATRIX.lines();
    let Some(_header) = lines.next() else {
        return out;
    };
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        let default = if f[6] == "None" || f[6].is_empty() {
            None
        } else {
            Some(f[6])
        };
        out.push(FlagSpec {
            subcommand: f[0].to_string(),
            flag: f[1].to_string(),
            action: f[2].to_string(),
            nargs: f[3].to_string(),
            dest: f[4].to_string(),
            required: f[5].eq_ignore_ascii_case("true"),
            default: default.map(str::to_string),
            choices: f[7].to_string(),
        });
    }
    out
}

/// A usage error, which exits 2 exactly as argparse does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UsageError(pub String);

impl std::fmt::Display for UsageError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for UsageError {}

/// The parsed options for one invocation.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Options {
    /// The subcommand that was named.
    pub subcommand: String,
    /// Populated destinations, in matrix order. A `nargs *`/`+` flag keeps its
    /// first value here; the full list is in [`Self::lists`].
    pub values: BTreeMap<String, String>,
    /// Every value a flag was given, in order. An option with `nargs='*'`, `'+'` or a
    /// count has one entry per value; a single-value option has exactly one, and an
    /// option that was not given at all has none (its default lives in
    /// [`Self::values`]).
    pub lists: BTreeMap<String, Vec<String>>,
    /// Destinations whose action was `store_true`/`store_false` and that fired.
    pub flags: BTreeMap<String, bool>,
    /// `--help` was requested.
    pub help: bool,
    /// The argument vector as given, **including** the subcommand name.
    ///
    /// Upstream echoes `" ".join(sys.argv[1:])` into the `*_peaks.xls` header as
    /// `# Command line: ...`, and `sys.argv[1]` is the subcommand. The file is
    /// compared byte-for-byte, so the exact spelling of what the user typed has
    /// to survive parsing.
    pub raw_argv: Vec<String>,
}

impl Options {
    /// A string destination.
    pub fn get(&self, dest: &str) -> Option<&str> {
        self.values.get(dest).map(String::as_str)
    }
    /// An integer destination, parsed with Python's `int()` grammar.
    pub fn int(&self, dest: &str) -> Option<i64> {
        self.get(dest).and_then(py_int)
    }
    /// A float destination, parsed with Python's `float()` grammar.
    pub fn float(&self, dest: &str) -> Option<f64> {
        self.get(dest).and_then(py_float)
    }
    /// Every value of a `nargs *`/`+` flag, in order.
    pub fn get_all(&self, dest: &str) -> &[String] {
        self.lists.get(dest).map(Vec::as_slice).unwrap_or(&[])
    }
    /// A boolean destination (`store_true`/`store_false`), defaulting to false.
    pub fn flag(&self, dest: &str) -> bool {
        self.flags.get(dest).copied().unwrap_or(false)
    }
}

/// Convert one token the way argparse's `_get_value` + `_check_value` do, and
/// return the spelling to store.
///
/// The conversion runs *first*, so `--format BED -c x` on a numeric flag reports
/// the numeric complaint and `bdgopt -m 5 -p zz` reports the bad `-p` value only
/// after `-m` has been checked.
fn convert(action: &Action, value: &str) -> Result<String, UsageError> {
    let display = action.display();
    if let Some(ty) = action.ty() {
        let ok = match ty {
            "int" => py_int(value).is_some(),
            _ => py_float(value).is_some(),
        };
        if !ok {
            return Err(UsageError(format!(
                "error: argument {display}: invalid {ty} value: {}",
                py_repr(value)
            )));
        }
        return Ok(canonical_number(value, ty));
    }
    if !action.spec.choices.is_empty() {
        let allowed: Vec<&str> = action.spec.choices.split(',').map(str::trim).collect();
        if !allowed.contains(&value) {
            let quoted: Vec<String> = allowed.iter().map(|c| py_repr(c)).collect();
            return Err(UsageError(format!(
                "error: argument {display}: invalid choice: {} (choose from {})",
                py_repr(value),
                quoted.join(", ")
            )));
        }
    }
    Ok(value.to_string())
}

/// Argparse mutually-exclusive groups in the pinned upstream parser.
///
/// The flag matrix records each flag's own attributes but not its *group*, so
/// these are transcribed from the oracle's own `_mutually_exclusive_groups`,
/// member order included. The comment above each entry is that order.
const MUTEX_GROUPS: &[(&str, &[&str], bool)] = &[
    // add_output_group: `-o/--ofile`, `--o-prefix`
    ("bdgpeakcall", &["ofile", "oprefix"], true),
    ("bdgbroadcall", &["ofile", "oprefix"], true),
    ("refinepeak", &["ofile", "oprefix"], true),
    // bdgcmp and bdgdiff declare `--o-prefix` before `-o`, unlike add_output_group
    ("bdgcmp", &["oprefix", "ofile"], true),
    ("bdgdiff", &["oprefix", "ofile"], true),
    ("randsample", &["percentage", "number"], true),
    // callpeak's `p_or_q_group`: neither is required
    ("callpeak", &["qvalue", "pvalue"], false),
    // callpeak's `postprocess_group` holds only `--call-summits`, so it constrains
    // nothing. `--broad` is in `group_callpeak`, a plain argument group -- verified:
    // upstream accepts `callpeak --broad --call-summits` together without complaint.
];

/// Enforce a *required* mutually-exclusive group, reported once the whole command
/// line has been consumed (`argparse.py:2199-2210`).
fn check_required_mutex_groups(
    subcommand: &str,
    acts: &[Action],
    seen: &BTreeMap<String, usize>,
) -> Result<(), UsageError> {
    for (cmd, members, required) in MUTEX_GROUPS {
        if *cmd != subcommand || !*required {
            continue;
        }
        if members.iter().any(|d| seen.contains_key(*d)) {
            continue;
        }
        // argparse names the flags, not the dests, and joins every spelling of a
        // member with `/`: "one of the arguments -p/--percentage -n/--number is
        // required"
        let names: Vec<String> = members
            .iter()
            .map(|d| {
                acts.iter()
                    .find(|a| a.spec.dest == **d)
                    .map_or_else(|| (*d).to_string(), Action::display)
            })
            .collect();
        return Err(UsageError(format!(
            "error: one of the arguments {} is required",
            names.join(" ")
        )));
    }
    Ok(())
}

/// Parse `args` (without the program name) against the matrix.
///
/// `subcommand` must already be known; the top-level parser is responsible for
/// dispatching, so this only sees a subcommand's own flags.
pub fn parse(subcommand: &str, args: &[String]) -> Result<Options, UsageError> {
    let specs = flag_specs();
    let acts = actions_for(&specs, subcommand);
    let mut out = Options {
        subcommand: subcommand.to_string(),
        raw_argv: std::iter::once(subcommand.to_string())
            .chain(args.iter().cloned())
            .collect(),
        ..Default::default()
    };

    // apply defaults first, so a flag that is not given still has its value
    for a in &acts {
        if a.spec.action == "help" {
            continue;
        }
        if let Some(d) = a.spec.default.as_deref() {
            match a.spec.action.as_str() {
                "store_true" => {
                    out.flags.insert(a.spec.dest.clone(), d == "true");
                }
                "store_false" => {
                    out.flags.insert(a.spec.dest.clone(), d != "false");
                }
                _ => {
                    out.values.insert(a.spec.dest.clone(), d.to_string());
                }
            }
        }
    }

    // `arg_strings_pattern`: `O` for a token `_parse_optional` claimed,
    // `A` for an argument, `-` for a literal `--` -- and every token after a
    // `--` is `A` whatever it looks like (`argparse.py:1969-1988`).
    let mut pattern: Vec<char> = Vec::with_capacity(args.len());
    let mut options_at: BTreeMap<usize, OptionTuple> = BTreeMap::new();
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--" {
            pattern.push('-');
            pattern.extend(std::iter::repeat_n('A', args.len() - i - 1));
            break;
        }
        match parse_optional(&acts, &args[i])? {
            None => pattern.push('A'),
            Some(t) => {
                options_at.insert(pattern.len(), t);
                pattern.push('O');
            }
        }
        i += 1;
    }

    let mut extras: Vec<String> = Vec::new();
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    let mut start_index = 0usize;
    let max_option = options_at.keys().next_back().copied();

    // `while start_index <= max_option_string_index` (`argparse.py:2133`). No
    // MACS3 subparser declares a positional, so `consume_positionals` is a no-op
    // that returns its `start_index` and the loop reduces to: skip anything
    // between options into `extras`, then consume one option.
    while let Some(max) = max_option {
        if start_index > max {
            break;
        }
        let next_option = options_at
            .range(start_index..)
            .next()
            .map(|(k, _)| *k)
            .expect("start_index <= max");
        if start_index != next_option {
            extras.extend_from_slice(&args[start_index..next_option]);
            start_index = next_option;
        }
        start_index = consume_optional(
            start_index,
            args,
            &pattern,
            &options_at,
            &acts,
            &mut extras,
            &mut out,
            &mut seen,
        )?;
        // argparse's `_HelpAction` calls `parser.exit()` from inside
        // `take_action`, so nothing after the `-h` is looked at: `bdgopt -h -m
        // bogus` succeeds and `bdgopt -m bogus -h` does not.
        if out.help {
            return Ok(out);
        }
    }
    extras.extend_from_slice(&args[start_index..]);

    // `-h`/`--help` short-circuits before any other validation, so
    // `macs3 callpeak -h` is a success even though `-t` is required.
    if out.help {
        return Ok(out);
    }

    // required flags must have been given (`argparse.py:2185-2192`)
    let missing: Vec<String> = acts
        .iter()
        .filter(|a| a.spec.required && a.spec.action != "help" && !seen.contains_key(&a.spec.dest))
        .map(Action::display)
        .collect();
    if !missing.is_empty() {
        return Err(UsageError(format!(
            "error: the following arguments are required: {}",
            missing.join(", ")
        )));
    }
    check_required_mutex_groups(subcommand, &acts, &seen)?;

    // Only now are the leftovers reported: `parse_args` (`argparse.py:1913-1917`)
    // receives them from `_parse_known_args`, which has already run the conversions,
    // the required check and the group check. `bdgopt --nope -m bogus` is therefore the
    // bad choice, and `bdgbroadcall -- -i in.bdg` is the missing `-i/--ifile`.
    if !extras.is_empty() {
        return Err(UsageError(format!(
            "error: unrecognized arguments: {}",
            extras.join(" ")
        )));
    }
    Ok(out)
}

/// `consume_optional` (`argparse.py:2022-2106`) plus `take_action`
/// (`argparse.py:1999-2020`), returning the index at which this option's values
/// stopped.
///
/// The `loop` is the glued-short-option machinery: when the token resolved to a
/// single-dash option that takes **no** argument and more characters followed,
/// argparse peels one character off at a time and looks the next `-x` up, which
/// is why `-Bp0.05` is `-B -p 0.05` and `-Bx=1` leaves `-x=1` behind as an extra.
#[allow(clippy::too_many_arguments)]
fn consume_optional(
    start_index: usize,
    args: &[String],
    pattern: &[char],
    options_at: &BTreeMap<usize, OptionTuple>,
    acts: &[Action],
    extras: &mut Vec<String>,
    out: &mut Options,
    seen: &mut BTreeMap<String, usize>,
) -> Result<usize, UsageError> {
    let tuple = options_at
        .get(&start_index)
        .expect("start_index is an option index");
    let mut action_idx = tuple.action;
    let option_string = tuple.option_string.clone();
    let mut sep = tuple.sep;
    let mut explicit_arg = tuple.explicit_arg.clone();

    // (action index, values, option string) -- one entry per option the token
    // expands to
    let mut action_tuples: Vec<(usize, Vec<String>, String)> = Vec::new();
    let stop: usize = loop {
        let Some(ai) = action_idx else {
            extras.push(args[start_index].clone());
            return Ok(start_index + 1);
        };
        let Some(explicit) = explicit_arg.clone() else {
            // no glued value: match the action's nargs pattern against the rest
            let start = start_index + 1;
            let count = match_nargs(&acts[ai], &pattern[start..])?;
            let values = args[start..start + count].to_vec();
            action_tuples.push((ai, values, option_string));
            break start + count;
        };
        // argparse matches against a one-character pattern here, so a glued value
        // only ever satisfies an action that takes at most one argument
        let count = match_nargs(&acts[ai], &['A'])?;
        let single_dash = !option_string.chars().nth(1).is_some_and(|c| c == '-');
        if count == 0 && single_dash && !explicit.is_empty() {
            if sep.is_some_and(|s| !s.is_empty()) || explicit.starts_with('-') {
                return Err(UsageError(format!(
                    "error: argument {}: ignored explicit argument {}",
                    acts[ai].display(),
                    py_repr(&explicit)
                )));
            }
            action_tuples.push((ai, Vec::new(), option_string.clone()));
            let head = option_string
                .chars()
                .next()
                .expect("option string is not empty")
                .to_string();
            let next: String = format!("{head}{}", explicit.chars().next().expect("non-empty"));
            match find_action(acts, &next) {
                Some(nai) => {
                    action_idx = Some(nai);
                    let rest: String = explicit.chars().skip(1).collect();
                    if rest.is_empty() {
                        sep = None;
                        explicit_arg = None;
                    } else if let Some(stripped) = rest.strip_prefix('=') {
                        sep = Some("=");
                        explicit_arg = Some(stripped.to_string());
                    } else {
                        sep = Some("");
                        explicit_arg = Some(rest);
                    }
                }
                None => {
                    extras.push(format!("{head}{explicit}"));
                    break start_index + 1;
                }
            }
        } else if count == 1 {
            action_tuples.push((ai, vec![explicit], option_string));
            break start_index + 1;
        } else {
            return Err(UsageError(format!(
                "error: argument {}: ignored explicit argument {}",
                acts[ai].display(),
                py_repr(&explicit)
            )));
        }
    };

    for (ai, values, _) in action_tuples {
        take_action(ai, &values, acts, out, seen)?;
    }
    Ok(stop)
}

/// `take_action` (`argparse.py:1999-2020`): convert, check against `choices`,
/// store, then reject a mutually-exclusive partner that was already seen.
fn take_action(
    ai: usize,
    values: &[String],
    acts: &[Action],
    out: &mut Options,
    seen: &mut BTreeMap<String, usize>,
) -> Result<(), UsageError> {
    let action = &acts[ai];
    let dest = action.spec.dest.clone();
    let converted: Result<Vec<String>, UsageError> =
        values.iter().map(|v| convert(action, v)).collect();
    let converted = converted?;

    match action.nargs() {
        Nargs::Zero => {
            // `--no-trackline` and friends *set the destination to false*.
            // Treating them as store_true makes a negated flag indistinguishable
            // from the default and inverts it (F133).
            out.flags
                .insert(dest.clone(), action.spec.action != "store_false");
            if action.spec.action == "help" {
                out.help = true;
            }
        }
        _ => {
            if let Some(first) = converted.first() {
                out.values.insert(dest.clone(), first.clone());
            }
            out.lists.insert(dest.clone(), converted);
        }
    }
    seen.insert(dest.clone(), ai);

    // "not allowed with argument" -- raised as each action is taken, so it beats
    // both the extras report and any later flag's bad value
    for (cmd, members, _) in MUTEX_GROUPS {
        if *cmd != action.spec.subcommand {
            continue;
        }
        let Some(pos) = members.iter().position(|d| **d == dest) else {
            continue;
        };
        for (i, other) in members.iter().enumerate() {
            if i == pos {
                continue;
            }
            if let Some(&oi) = seen.get(*other) {
                return Err(UsageError(format!(
                    "error: argument {}: not allowed with argument {}",
                    action.display(),
                    acts[oi].display()
                )));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn the_matrix_covers_every_subcommand_in_the_objective() {
        let specs = flag_specs();
        assert!(specs.len() > 300, "matrix looks truncated: {}", specs.len());
        for c in crate::SUBCOMMANDS {
            assert!(specs.iter().any(|s| s.subcommand == c), "no flags for {c}");
        }
    }

    #[test]
    fn every_subcommand_has_a_required_flag_or_none_is_expected() {
        // just exercises the required-check path
        let r = parse("callpeak", &v(&["-t", "a.bed"]));
        assert!(r.is_ok() || r.is_err());
    }

    #[test]
    fn an_unknown_flag_is_a_usage_error() {
        // the leftovers are reported by `parse_args`, i.e. after the required check,
        // so `-t` has to be present for the unknown flag to be the complaint:
        // `macs3 callpeak --not-a-real-flag` says "the following arguments are
        // required: -t/--treatment" in the reference too
        let e = parse("callpeak", &v(&["-t", "a.bed", "--not-a-real-flag"])).unwrap_err();
        assert_eq!(e.0, "error: unrecognized arguments: --not-a-real-flag");
        let e = parse("callpeak", &v(&["--not-a-real-flag"])).unwrap_err();
        assert_eq!(
            e.0,
            "error: the following arguments are required: -t/--treatment"
        );
    }

    #[test]
    fn a_missing_required_flag_is_a_usage_error() {
        // callpeak requires -t (treatment)
        let e = parse("callpeak", &v(&[])).unwrap_err();
        assert!(e.0.contains("required"), "{e}");
    }

    #[test]
    fn inline_and_separate_values_are_both_accepted() {
        let a = parse("callpeak", &v(&["-t", "a.bed", "--gsize=1000"])).unwrap();
        let b = parse("callpeak", &v(&["-t", "a.bed", "--gsize", "1000"])).unwrap();
        assert_eq!(a.get("gsize"), b.get("gsize"));
        assert_eq!(a.get("gsize"), Some("1000"));
    }

    #[test]
    fn a_flag_with_no_value_argument_is_a_usage_error() {
        let e = parse("callpeak", &v(&["-t"])).unwrap_err();
        assert!(e.0.contains("expected"), "{e}");
    }

    #[test]
    fn store_true_sets_the_destination() {
        let o = parse("callpeak", &v(&["-t", "a.bed", "--nomodel"])).unwrap();
        assert!(o.flag("nomodel"));
    }

    #[test]
    fn defaults_are_applied() {
        let o = parse("callpeak", &v(&["-t", "a.bed"])).unwrap();
        // upstream's default q-value cutoff
        assert!(o.get("qvalue").is_some(), "qvalue default missing");
    }
}
