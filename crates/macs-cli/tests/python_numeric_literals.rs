//! Regression: an `int`/`float` flag's value is checked with **Python's own literal
//! grammar**, because argparse's `type=` *is* `int`/`float`.
//!
//! `_get_value` (`argparse.py:2586-2596`) calls `type_func(arg_string)` and turns a
//! `ValueError` into `invalid <type> value: <repr(token)>`, exit 2. The accept set is
//! therefore whatever `int()`/`float()` accept, which is strictly larger than
//! `str::parse::<i64>()`/`str::parse::<f64>()`:
//!
//! | token                     | `int()`  | `float()`  | Rust `parse` |
//! |---------------------------|----------|------------|--------------|
//! | `" 1"`, `"1 "`, `"\t3"`   | 1        | 1.0        | rejected      |
//! | `"1_0"`, `"1_000.5"`      | error    | 10, 1000.5 | rejected      |
//! | `"+1"`                    | 1        | 1.0        | accepted      |
//! | `"1e3"`                   | error    | 1000.0     | accepted      |
//! | `"inf"`, `"nan"`, `"infinity"` | error | ±inf, nan  | accepted      |
//! | `"١"`, `"１２"`            | 1        | 1.0        | rejected      |
//! | `"1.0"`, `"1e3"` on `type=int` | error | --     | --            |
//!
//! `bdgopt -p ' -1'` exited 2 here and 0 upstream, and because
//! [`macs_cli::Options::int`] / [`macs_cli::Options::float`] return `None` on a parse
//! failure, the old check was an accept/reject divergence *and* a silent wrong answer
//! risk.
//!
//! The tables below are the reference's own classification, obtained by running
//! `.oracle/venv/bin/macs3` with `-i /nonexistent -o /dev/null` so that an accepted
//! value fails at the open with **exit 1** and a rejected one exits **2**. They are
//! embedded here so the test needs no Python and no reference.
//!
//! ```text
//! bdgopt -i /nonexistent -o /dev/null -m multiply -p 1_000.5   1  FileNotFoundError
//! bdgopt -i /nonexistent -o /dev/null -m multiply -p ١٢٣     1  FileNotFoundError
//! bdgopt -i /nonexistent -o /dev/null -m multiply -p 1_0_     2  invalid float value: '1_0_'
//! bdgopt -i /nonexistent -o /dev/null -m multiply -p 1.2.3    2  invalid float value: '1.2.3'
//! bdgopt -i /nonexistent -o /dev/null -m multiply -p ''       2  invalid float value: ''
//! filterdup -i /nonexistent -f BED -o /dev/null -s 1_0        1  FileNotFoundError
//! filterdup -i /nonexistent -f BED -o /dev/null -s ١٢٣        1  FileNotFoundError
//! filterdup -i /nonexistent -f BED -o /dev/null -s 1.0        2  invalid int value: '1.0'
//! filterdup -i /nonexistent -f BED -o /dev/null -s inf        2  invalid int value: 'inf'
//! filterdup -i /nonexistent -f BED -o /dev/null -s "it's"    2  invalid int value: "it's"
//! filterdup -i /nonexistent -f BED -o /dev/null -s "1\t2"    2  invalid int value: '1\t2'
//! filterdup -i /nonexistent -f BED -o /dev/null -s 9223372036854775808
//!                                                                    1  OverflowError
//! ```
//!
//! Note the two spellings that look negative but are **not** values at all:
//! `-1.`, `-1e-3`, `-9x` and `-inf` fail argparse's `_negative_number_matcher`
//! (`^-\d+$|^-\d*\.\d+$`, `argparse.py:1420`) when they stand alone, so they are
//! unknown *options* -- and, after `bdgopt -p`, `-inf` is even `-i nf`, because `-i` is
//! a short option. Glued, they are values. See
//! [`an_option_looking_token_needs_the_glued_spelling`].

use macs_cli::{parse_flags, Options};

fn parse(sub: &str, words: &[&str]) -> Result<Options, String> {
    let args: Vec<String> = words.iter().map(std::string::ToString::to_string).collect();
    parse_flags(sub, &args).map_err(|e| e.0)
}

/// The token set for the reference's `bdgopt -p/--extra-param` (`type=float`,
/// `nargs='*'`), with the value Python's `float()` produces.
const FLOAT_OK: &[(&str, f64)] = &[
    ("0", 0.0),
    ("1", 1.0),
    ("-1", -1.0),
    ("+1", 1.0),
    ("-0", -0.0),
    (" 1", 1.0),
    ("1 ", 1.0),
    ("\t3", 3.0),
    ("2\n", 2.0),
    ("  2  ", 2.0),
    ("1_0", 10.0),
    ("1_000.5", 1000.5),
    ("1_0_0", 100.0),
    ("0_1", 1.0),
    ("1_0e3", 10_000.0),
    ("1e3_0", 1e30),
    ("1_2_3.4_5_6e7_8", 1.23456e80),
    ("+1.5", 1.5),
    (".5", 0.5),
    ("5.", 5.0),
    ("1.", 1.0),
    ("-.5", -0.5),
    ("1e3", 1000.0),
    ("1E3", 1000.0),
    ("1e+3", 1000.0),
    ("1e-3", 0.001),
    ("+1e5", 100_000.0),
    ("1e400", f64::INFINITY),
    ("1e-400", 0.0),
    ("inf", f64::INFINITY),
    ("Inf", f64::INFINITY),
    ("INF", f64::INFINITY),
    ("infinity", f64::INFINITY),
    ("Infinity", f64::INFINITY),
    ("nan", f64::NAN),
    ("NaN", f64::NAN),
    ("\u{0661}", 1.0),          // Arabic-Indic digit one
    ("\u{ff11}", 1.0),          // fullwidth digit one
    ("\u{ff10}.\u{ff15}", 0.5), // fullwidth zero, fullwidth five
    ("\u{0661}\u{0662}\u{0663}", 123.0),
];

/// Tokens the reference rejects with `invalid float value`, quoted as Python would.
const FLOAT_BAD: &[&str] = &[
    "",
    "   ",
    "zz",
    "1_0_",
    "1__0",
    "_1",
    "1_",
    "__1",
    "1_.5",
    "1._5",
    "1.5_",
    ".5_",
    "1_.",
    "._1",
    "1e",
    "1e_3",
    "e3",
    "0x10",
    "1,0",
    "+ 1",
    "- 1",
    "1 e3",
    "inf1",
    "1inf",
    "nan(",
    ".",
    "..",
    "1.2.3",
    "\u{0661}\u{0662}x",
    "nan nan",
];

/// The token set for the reference's `filterdup -s/--tsize` (`type=int`).
const INT_OK: &[(&str, i64)] = &[
    ("0", 0),
    ("1", 1),
    ("-1", -1),
    ("+1", 1),
    ("-0", 0),
    (" 1", 1),
    ("1 ", 1),
    ("\t3", 3),
    ("2\n", 2),
    ("  2  ", 2),
    ("1_0", 10),
    ("1_0_0", 100),
    ("0_1", 1),
    ("1\u{0661}", 11),
    ("\u{0661}\u{0662}\u{0663}", 123),
    ("\u{ff11}\u{ff12}", 12),
];

/// Tokens the reference rejects with `invalid int value`.
const INT_BAD: &[&str] = &[
    "", "   ", "zz", "1.0", "1.5", ".5", "5.", "1e3", "1E3", "1e-3", "inf", "nan", "infinity",
    "0x10", "1,0", "_1", "1_", "1__0", "+ 1", "1 2",
];

#[test]
fn a_float_flag_accepts_everything_python_float_accepts() {
    for (token, want) in FLOAT_OK {
        let words = ["-i", "in.bdg", "-o", "o.bdg", "-m", "multiply", "-p", token];
        let o = parse("bdgopt", &words)
            .unwrap_or_else(|e| panic!("bdgopt -p {token:?} was rejected: {e}"));
        assert_eq!(o.get_all("extraparam").len(), 1, "{token:?}");
        let got = o
            .float("extraparam")
            .unwrap_or_else(|| panic!("{token:?} did not parse to a number: {o:?}"));
        if want.is_nan() {
            assert!(got.is_nan(), "{token:?}: {got}");
        } else {
            assert_eq!(got, *want, "{token:?}");
            // the sign of zero is observable through the pipeline, so compare it
            assert_eq!(got.is_sign_negative(), want.is_sign_negative(), "{token:?}");
        }
    }
}

#[test]
fn a_float_flag_rejects_everything_python_float_rejects() {
    for token in FLOAT_BAD {
        let words = ["-i", "in.bdg", "-o", "o.bdg", "-m", "multiply", "-p", token];
        let Err(e) = parse("bdgopt", &words) else {
            panic!("bdgopt -p {token:?} was accepted")
        };
        assert_eq!(
            e,
            format!(
                "error: argument -p/--extra-param: invalid float value: {}",
                quote(token)
            ),
            "{token:?}"
        );
    }
}

#[test]
fn an_int_flag_accepts_everything_python_int_accepts() {
    for (token, want) in INT_OK {
        let words = ["-i", "in.bed", "-f", "BED", "-o", "o.bed", "-s", token];
        let o = parse("filterdup", &words)
            .unwrap_or_else(|e| panic!("filterdup -s {token:?} was rejected: {e}"));
        assert_eq!(o.int("tsize"), Some(*want), "filterdup -s {token:?}: {o:?}");
    }
}

#[test]
fn an_int_flag_rejects_everything_python_int_rejects() {
    for token in INT_BAD {
        let words = ["-i", "in.bed", "-f", "BED", "-o", "o.bed", "-s", token];
        let Err(e) = parse("filterdup", &words) else {
            panic!("filterdup -s {token:?} was accepted")
        };
        assert_eq!(
            e,
            format!(
                "error: argument -s/--tsize: invalid int value: {}",
                quote(token)
            ),
            "{token:?}"
        );
    }
}

/// Python's `repr()` of a `str`, which is what `%(value)r` interpolates.
fn quote(s: &str) -> String {
    let q = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(q);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == q => {
                out.push('\\');
                out.push(c);
            }
            '\u{0}'..='\u{1f}' | '\u{7f}' => out.push_str(&format!("\\x{:02x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push(q);
    out
}

#[test]
fn an_int_beyond_i64_is_accepted_and_saturated() {
    // Python's `int` is arbitrary precision, so argparse accepts these and the failure
    // is a runtime `OverflowError` inside the pipeline. Saturating keeps the flag from
    // silently reverting to its default, which would be a wrong answer rather than a
    // crash. The magnitude is accumulated first and the sign applied after, so a huge
    // negative token saturates to `i64::MIN` rather than through `i64::MAX`.
    for (token, want) in [
        ("9223372036854775807", i64::MAX),
        ("9223372036854775808", i64::MAX),
        ("99999999999999999999999", i64::MAX),
        ("-9223372036854775807", -9223372036854775807),
        ("-9223372036854775808", i64::MIN),
        ("-99999999999999999999999", i64::MIN),
    ] {
        let words = ["-i", "in.bed", "-f", "BED", "-o", "o.bed", "-s", token];
        let o = parse("filterdup", &words).unwrap_or_else(|e| panic!("-s {token}: {e}"));
        assert_eq!(o.int("tsize"), Some(want), "{token}");
    }
}

#[test]
fn an_option_looking_token_needs_the_glued_spelling() {
    // `-1.`, `-1e-3` and `-9x` do not match `_negative_number_matcher`, so standing
    // alone they are unknown options and the value position is left empty:
    //   bdgopt ... -p -1.     2  macs3: error: unrecognized arguments: -1.
    //   filterdup ... -s -9x  2  macs3 filterdup: error: argument -s/--tsize: expected one argument
    // After `bdgopt -p`, `-inf` is not even unknown: `-i` is a short option there, so
    // `-inf` is `-i nf` and the reference reports
    //   1  ERROR Need EXTRAPARAM for method multiply or add!
    for token in ["-1.", "-1e-3", "-9x"] {
        let Err(e) = parse(
            "bdgopt",
            &["-i", "in.bdg", "-o", "o", "-m", "multiply", "-p", token],
        ) else {
            panic!("bdgopt -p {token} was accepted as a separate token")
        };
        assert_eq!(
            e,
            format!("error: unrecognized arguments: {token}"),
            "{token}"
        );
        let Err(e) = parse(
            "filterdup",
            &["-i", "in.bed", "-f", "BED", "-o", "o", "-s", token],
        ) else {
            panic!("filterdup -s {token} was accepted as a separate token")
        };
        assert_eq!(
            e, "error: argument -s/--tsize: expected one argument",
            "{token}"
        );
    }

    // Glued, they never reach `_parse_optional` at all -- `_get_option_tuples` puts the
    // tail in `explicit_arg` -- so they are values:
    //   bdgopt ... -p-1.  1  (and it writes -5 / -7 for a 5 / 7 track)
    //   filterdup ... -s-inf  2  invalid int value: '-inf'
    let o = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o", "-m", "multiply", "-p-1."],
    )
    .unwrap();
    assert_eq!(o.float("extraparam"), Some(-1.0));
    let o = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o", "-m", "multiply", "-p-1e-3"],
    )
    .unwrap();
    assert_eq!(o.float("extraparam"), Some(-1e-3));
    let o = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o", "-m", "multiply", "-p-1e400"],
    )
    .unwrap();
    assert_eq!(o.float("extraparam"), Some(f64::NEG_INFINITY));
    let o = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o", "-m", "multiply", "-p-inf"],
    )
    .unwrap();
    assert_eq!(o.float("extraparam"), Some(f64::NEG_INFINITY));
    let Err(e) = parse(
        "filterdup",
        &["-i", "in.bed", "-f", "BED", "-o", "o", "-s-inf"],
    ) else {
        panic!("-s-inf was accepted")
    };
    assert_eq!(e, "error: argument -s/--tsize: invalid int value: '-inf'");
    let Err(e) = parse(
        "filterdup",
        &["-i", "in.bed", "-f", "BED", "-o", "o", "-s-inf0"],
    ) else {
        panic!("-s-inf0 was accepted")
    };
    assert_eq!(e, "error: argument -s/--tsize: invalid int value: '-inf0'");
}

#[test]
fn the_stored_value_is_always_parseable_by_the_commands() {
    // The commands read `Options::values` and call `str::parse` themselves (see
    // `bedgraph_cmds.rs` on `--extra-param`), so a token Python accepts and Rust does
    // not has to be normalised before they see it, or `bdgopt -p ' 2 '` would fail at
    // "must be a number" instead of producing upstream's output.
    for (token, _) in FLOAT_OK {
        let words = ["-i", "in.bdg", "-o", "o.bdg", "-m", "multiply", "-p", token];
        let o = parse("bdgopt", &words).unwrap();
        let stored = o.get("extraparam").expect("stored");
        assert!(
            stored.parse::<f64>().is_ok(),
            "{token:?} was stored as {stored:?}, which a command cannot parse"
        );
    }
    for (token, _) in INT_OK {
        let words = ["-i", "in.bed", "-f", "BED", "-o", "o.bed", "-s", token];
        let o = parse("filterdup", &words).unwrap();
        let stored = o.get("tsize").expect("stored");
        assert!(
            stored.parse::<i64>().is_ok(),
            "{token:?} was stored as {stored:?}, which a command cannot parse"
        );
    }
}

#[test]
fn a_token_both_parsers_agree_on_is_stored_verbatim() {
    // Rewriting a token Rust already understands would change what the pipeline and
    // the `# Command line:` header report for no reason, so it must not happen.
    for token in [
        "-0.35", "-1", "-.5", "-0", "0.5", "1e-5", "1E3", "inf", "nan", "50",
    ] {
        let words = ["-i", "in.bdg", "-o", "o.bdg", "-m", "multiply", "-p", token];
        let o = parse("bdgopt", &words).unwrap_or_else(|e| panic!("{token}: {e}"));
        assert_eq!(o.get("extraparam"), Some(token), "{token} was rewritten");
    }
}

#[test]
fn an_exotic_spelling_produces_the_same_bytes_as_the_plain_one() {
    // the accept/reject decision is not the whole point: the *value* has to reach the
    // pipeline. `bdgopt -m multiply -p '1_0'` must write exactly what `-p 10` writes.
    let dir = std::env::temp_dir().join(format!("macs3rs-pynum-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("in.bdg");
    std::fs::write(&input, "chr1\t0\t10\t5\nchr1\t10\t20\t7\n").unwrap();

    let mut outs = Vec::new();
    for (tag, words) in [
        ("plain", vec!["-p", "10"]),
        ("underscore", vec!["-p", "1_0"]),
        ("glued", vec!["-p1_0"]),
        ("abbrev-eq", vec!["--extra-param=1_0"]),
        ("whitespace", vec!["-p", " 1_0 "]),
    ] {
        let out = dir.join(format!("{tag}.bdg"));
        let mut argv: Vec<String> = vec![
            "-i".into(),
            input.to_str().unwrap().into(),
            "-o".into(),
            out.to_str().unwrap().into(),
            "-m".into(),
            "multiply".into(),
        ];
        argv.extend(words.into_iter().map(std::string::ToString::to_string));
        let o = parse_flags("bdgopt", &argv).unwrap_or_else(|e| panic!("{tag}: {e}"));
        macs_cli::commands::bedgraph_cmds::bdgopt(&o).expect("bdgopt runs");
        outs.push((tag, std::fs::read(&out).expect("written")));
    }
    let (_, want) = &outs[0];
    // upstream's bytes for `macs3 bdgopt -m multiply -p 10`
    assert!(String::from_utf8_lossy(want).contains("chr1\t0\t10\t50.00000"));
    for (tag, got) in &outs {
        assert_eq!(got, want, "{tag}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn an_int_flag_reaches_the_pipeline_with_its_python_value() {
    // `bdgpeakcall -l/--min-length` is an int: a 10 bp peak survives `-l 1_0` and is
    // dropped by `-l 100`, which is only true if the value arrived as 10 and 100.
    let dir = std::env::temp_dir().join(format!("macs3rs-pyint-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("in.bdg");
    std::fs::write(&input, "chr1\t0\t10\t5\nchr1\t20\t30\t5\n").unwrap();

    let call = |tag: &str, words: &[&str]| -> String {
        let out = dir.join(format!("{tag}.bed"));
        let mut argv: Vec<String> = vec![
            "-i".into(),
            input.to_str().unwrap().into(),
            "-o".into(),
            out.to_str().unwrap().into(),
        ];
        argv.extend(words.iter().map(std::string::ToString::to_string));
        let o = parse_flags("bdgpeakcall", &argv).unwrap_or_else(|e| panic!("{tag}: {e}"));
        macs_cli::commands::bdgpeakcall::bdgpeakcall(&o).expect("bdgpeakcall runs");
        // the trackline and the peak name both mention the output file, so compare
        // the coordinates only
        std::fs::read_to_string(&out)
            .expect("written")
            .lines()
            .filter(|l| !l.starts_with("track"))
            .map(|l| l.split('\t').take(3).collect::<Vec<_>>().join("\t"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(
        call("ten", &["-l", "10"]),
        call("underscore", &["-l", "1_0"])
    );
    assert_eq!(call("ten", &["-l", "10"]), call("glued", &["-l1_0"]));
    assert_eq!(call("ten", &["-l", "10"]), call("ws", &["-l", " 1_0 "]));
    assert_ne!(call("ten", &["-l", "10"]), call("hundred", &["-l", "100"]));
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn the_error_message_is_argparse_shape_with_pythons_repr() {
    let cases: &[(&str, &str)] = &[
        (
            "zz",
            "error: argument -p/--extra-param: invalid float value: 'zz'",
        ),
        (
            "1_0_",
            "error: argument -p/--extra-param: invalid float value: '1_0_'",
        ),
        (
            "it's",
            "error: argument -p/--extra-param: invalid float value: \"it's\"",
        ),
        (
            "",
            "error: argument -p/--extra-param: invalid float value: ''",
        ),
    ];
    for (token, msg) in cases {
        let words = ["-i", "in.bdg", "-o", "o.bdg", "-m", "multiply", "-p", token];
        let Err(e) = parse("bdgopt", &words) else {
            panic!("{token:?} was accepted")
        };
        assert_eq!(e, *msg, "{token:?}");
    }
    let Err(e) = parse(
        "filterdup",
        &["-i", "x", "-f", "BED", "-o", "o", "-s", "1\t2"],
    ) else {
        panic!("a tab was accepted")
    };
    assert_eq!(e, "error: argument -s/--tsize: invalid int value: '1\\t2'");
}

#[test]
fn an_invalid_choice_is_quoted_the_way_argparse_quotes_it() {
    // `_check_value` (`argparse.py:2611-2617`) interpolates `%(value)r` and
    // `', '.join(map(repr, action.choices))`, so both sides are quoted
    let Err(e) = parse("bdgopt", &["-i", "in.bdg", "-o", "o.bdg", "-m", "nope"]) else {
        panic!("-m nope was accepted")
    };
    assert_eq!(
        e,
        "error: argument -m/--method: invalid choice: 'nope' \
         (choose from 'multiply', 'add', 'p2q', 'max', 'min')"
    );
    let Err(e) = parse("bdgopt", &["-i", "in.bdg", "-o", "o.bdg", "-m", "5"]) else {
        panic!("-m 5 was accepted")
    };
    assert_eq!(
        e,
        "error: argument -m/--method: invalid choice: '5' \
         (choose from 'multiply', 'add', 'p2q', 'max', 'min')"
    );
}

#[test]
fn numeric_defaults_still_parse() {
    // the defaults come out of the matrix as strings and go through the same
    // accessors, so they must not depend on the token spelling
    let o = parse("callpeak", &["-t", "a.bed", "--help"]).unwrap();
    assert_eq!(o.float("qvalue"), Some(0.05));
    assert_eq!(o.int("verbose"), Some(2));
    assert_eq!(o.int("extsize"), Some(200));
    assert_eq!(o.get("mfold"), Some("5,50"));
    assert!(o.get_all("mfold").is_empty());
}

#[test]
fn the_comparison_and_a_glued_value_use_the_same_grammar() {
    // `-p 1_0` and `-p1_0` must agree on both the accept decision and the value
    let sep = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o", "-m", "multiply", "-p", "1_0"],
    )
    .unwrap();
    let glued = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o", "-m", "multiply", "-p1_0"],
    )
    .unwrap();
    assert_eq!(sep.get("extraparam"), glued.get("extraparam"));
    assert_eq!(sep.float("extraparam"), glued.float("extraparam"));
    assert_eq!(glued.float("extraparam"), Some(10.0));
}

/// The two generated oracle tables, so the sweep below can enumerate every numeric
/// flag without the parser having to expose its type table.
const MATRIX: &str = include_str!("../../../oracle/flag_matrix.tsv");
const TYPES: &str = include_str!("../../../oracle/flag_types.tsv");

#[test]
fn every_numeric_flag_uses_the_python_grammar() {
    // `oracle/flag_types.tsv` types 142 spellings, of which 130 name an option the
    // pinned argparse actually registers (see [`the_type_table_has_twelve_stale_rows`]).
    // The check is per (subcommand, spelling), so a mechanical sweep is the only way to
    // be sure none of the live ones was left on Rust's grammar.
    let live: Vec<(String, String)> = MATRIX
        .lines()
        .skip(1)
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f.len() >= 8 && f[0] != "subcommand").then(|| (f[0].to_string(), f[1].to_string()))
        })
        .collect();
    let mut checked = 0;
    for line in TYPES.lines() {
        if line.trim().is_empty() || line.starts_with('#') || line.starts_with("subcommand") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let (sub, flag, ty) = (f[0], f[1], f[2]);
        if !live.iter().any(|(s, g)| s == sub && g == flag) {
            continue;
        }
        checked += 1;

        // Four tokens, because the widest `nargs` in the matrix is 4 (`hmmratac
        // --means`). Spare tokens are harmless: a conversion error is raised while the
        // option is being consumed, i.e. before any leftover token is reported.
        if let Err(e) = parse(sub, &[flag, "1_0", "1_0", "1_0", "1_0"]) {
            assert!(
                !e.contains("invalid"),
                "{sub} {flag} ({ty}) rejected 1_0 with a *conversion* error: {e}"
            );
        }

        // a token no numeric grammar accepts, which must be the conversion error
        let Err(e) = parse(sub, &[flag, "zz", "zz", "zz", "zz"]) else {
            panic!("{sub} {flag} ({ty}) accepted the token 'zz'")
        };
        assert!(
            e.contains(&format!("invalid {ty} value: 'zz'")),
            "{sub} {flag} ({ty}): unexpected error {e}"
        );
    }
    assert_eq!(
        checked, 130,
        "the sweep covered a different number of flags"
    );
}

#[test]
fn the_type_table_has_twelve_stale_rows() {
    // Not a defect: `oracle/flag_types.tsv` is derived by reading `add_argument` calls
    // out of `bin/macs3`, and upstream keeps a dozen of them commented out or behind a
    // disabled branch (`cmbreps -w`, `bin/macs3:465-466`; the `hmmratac` ones come from
    // a module list argparse never registers). `flag_matrix.tsv` is the live authority --
    // it is introspected from the parser objects -- so the type table is only ever
    // consulted for a spelling the matrix has, and a stale row is inert:
    //
    //   macs3 hmmratac -m 3 -i x   ->  2  macs3: error: unrecognized arguments: -m 3
    let live: Vec<(String, String)> = MATRIX
        .lines()
        .skip(1)
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f.len() >= 8 && f[0] != "subcommand").then(|| (f[0].to_string(), f[1].to_string()))
        })
        .collect();
    let stale: Vec<(String, String, String)> = TYPES
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#') && !l.starts_with("subcommand"))
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f.len() >= 3 && !live.iter().any(|(s, g)| s == f[0] && g == f[1]))
                .then(|| (f[0].to_string(), f[1].to_string(), f[2].to_string()))
        })
        .collect();
    assert_eq!(
        stale
            .iter()
            .map(|(s, f, _)| format!("{s} {f}"))
            .collect::<Vec<_>>(),
        [
            "cmbreps -w",
            "hmmratac --minmapq",
            "hmmratac --multiple-processing",
            "hmmratac --states",
            "hmmratac --threshold",
            "hmmratac --trim",
            "hmmratac --window",
            "hmmratac --zscore",
            "hmmratac -m",
            "hmmratac -q",
            "hmmratac -s",
            "hmmratac -z",
        ]
    );
    // and each one really is an unknown option here, not a live flag. (`cmbreps`
    // also requires `-i` and `-o`, and the required check comes first, so the
    // leftover is only named once those are supplied.)
    for (sub, flag, _) in &stale {
        let Err(e) = parse(sub, &[flag, "1"]) else {
            panic!("{sub} {flag} is a live option after all")
        };
        assert!(
            e.contains("unrecognized") || e.contains("required"),
            "{sub} {flag}: {e}"
        );
        // `cmbreps -i` is `nargs='+'` and required; `hmmratac` has no `-o` at all, so
        // one of its stale flags really is read as `-o`'s neighbour
        let base: &[&str] = match sub.as_str() {
            "cmbreps" => &["-i", "a.bed", "b.bed", "-o", "o.bed"],
            "hmmratac" => &["-i", "a.bed"],
            _ => &["-i", "a.bed", "-o", "o.bed"],
        };
        let argv: Vec<&str> = base.iter().copied().chain([flag, "1"]).collect();
        let Err(e) = parse(sub, &argv) else {
            panic!("{sub} {flag} is a live option after all")
        };
        assert_eq!(e, format!("error: unrecognized arguments: {flag} 1"));
    }
}
