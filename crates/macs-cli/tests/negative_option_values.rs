//! Regression: a negative number in a **separate** token is a value, not a flag.
//!
//! argparse decides that with `_parse_optional` (`argparse.py:2276-2334`),
//! which returns `None` -- "this is an argument, mark it `A`" -- for a token that
//!
//! * is empty, or does not start with a prefix character, or is the lone `-`;
//! * matches `_negative_number_matcher`, `^-\d+$|^-\d*\.\d+$`
//!   (`argparse.py:1420`), **and** the parser declares no option string that
//!   matches the same pattern (`argparse.py:1533-1535`); or
//! * contains a space.
//!
//! No flag in `oracle/flag_matrix.tsv` matches that pattern -- not even a
//! program-level one -- so every subparser's `_has_negative_number_optionals` is
//! empty and the clause always fires: `-0.35`, `-1` and `-.5` are values.
//!
//! The port used to stop consuming a `nargs='*'`/`'+'` run at any `-`-prefixed
//! token, so `bdgopt -i x.bdg -m add -p -0.35 -o o.bdg` exited 2 with
//! "unrecognized arguments: -0.35" where argparse accepts it and MACS3 runs.
//! `--extra-param=-0.35` worked because the `=` splits the token first, which is
//! why the bug survived next to the working spelling.
//!
//! The negative-number test is deliberately narrow: `-1.` and `-1e-3` do **not**
//! match it (no digits after the dot / no exponent branch), and `-9x` does not
//! match it at all, so all three stay unknown flags and keep exiting 2.

use macs_cli::{parse_flags, Options};

fn s(x: &str) -> String {
    x.to_string()
}

fn parse(sub: &str, words: &[&str]) -> Result<Options, String> {
    let args: Vec<String> = words.iter().map(|w| s(w)).collect();
    parse_flags(sub, &args).map_err(|e| e.0)
}

/// Upstream's bytes for `macs3 bdgopt -i in.bdg -m add -p -0.35 -o o.bdg`.
const ADD_P_NEG: &str = concat!(
    "track type=bedGraph name=\"ADD_modified_scores\" description=\"Scores calculated by ADD\" visibility=2 alwaysZero=on\n",
    "chr1\t0\t10\t4.65000\n",
    "chr1\t10\t20\t6.65000\n",
);

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-negsum-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn bdgopt_accepts_a_negative_extra_param_as_a_separate_token() {
    let dir = tmpdir("bdgopt");
    let input = dir.join("in.bdg");
    std::fs::write(&input, "chr1\t0\t10\t5\nchr1\t10\t20\t7\n").unwrap();
    let args = vec![
        s("-i"),
        input.to_str().unwrap().to_string(),
        s("-m"),
        s("add"),
        s("-p"),
        s("-0.35"),
        s("-o"),
        s("o.bdg"),
        s("--outdir"),
        dir.to_str().unwrap().to_string(),
    ];
    let o: Options = parse_flags("bdgopt", &args).expect("args parse");
    macs_cli::commands::bedgraph_cmds::bdgopt(&o).expect("bdgopt runs");
    let got = std::fs::read_to_string(dir.join("o.bdg")).expect("output written");
    assert_eq!(got, ADD_P_NEG);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn every_negative_number_spelling_is_a_value_for_every_subcommand() {
    // One representative numeric flag per subcommand, with the destination
    // `oracle/flag_matrix.tsv` records for it. `bdgopt -p` is the `nargs='*'` case
    // that used to fail; `hmmratac --means` is `nargs=4` and `callpeak -m` is
    // `nargs=2`; the rest take a single value.
    let cases: &[(&str, &str, &str, &[&str])] = &[
        ("bdgopt", "-p", "extraparam", &["-1", "-0.35", "-.5", "-0"]),
        (
            "bdgcmp",
            "--pseudocount",
            "pseudocount",
            &["-1", "-0.35", "-.5"],
        ),
        ("bdgpeakcall", "--cutoff", "cutoff", &["-1", "-0.35"]),
        (
            "bdgbroadcall",
            "--cutoff-link",
            "cutofflink",
            &["-1", "-0.35"],
        ),
        ("bdgdiff", "--cutoff", "cutoff", &["-1", "-0.35"]),
        ("callpeak", "--shift", "shift", &["-1", "-75", "-0"]),
        ("callvar", "--max-ar", "maxAR", &["-1", "-0.35"]),
        ("filterdup", "--pvalue", "pvalue", &["-1", "-0.35"]),
        ("hmmratac", "-c", "prescan_cutoff", &["-1", "-0.35"]),
        ("pileup", "--extsize", "extsize", &["-1"]),
        ("predictd", "--tsize", "tsize", &["-1"]),
        ("randsample", "--number", "number", &["-1", "-0.35"]),
        ("refinepeak", "--cutoff", "cutoff", &["-1", "-0.35"]),
    ];
    for &(sub, flag, dest, values) in cases {
        for &v in values {
            // `--help` short-circuits every other check, so a successful parse here
            // means the value was accepted, and `dest` proves it was consumed as
            // this flag's value rather than left over on the command line
            let o = parse(sub, &[flag, v, "--help"])
                .unwrap_or_else(|e| panic!("{sub} {flag} {v}: {e}"));
            assert_eq!(o.get(dest), Some(v), "{sub} {flag} {v}");
            assert!(o.help, "{sub} {flag} {v}");
        }
    }
}

#[test]
fn negative_values_are_accepted_by_fixed_and_variable_nargs_flags() {
    // `nargs=2`, `nargs=4` and `nargs='*'`
    let o = parse("callpeak", &["-m", "-3", "-20", "--help"]).expect("accepted");
    assert_eq!(o.get_all("mfold"), ["-3", "-20"]);
    let o = parse("hmmratac", &["--means", "-1", "-2", "-3", "-4", "--help"]).expect("accepted");
    assert_eq!(o.get_all("em_means"), ["-1", "-2", "-3", "-4"]);
    let o = parse("bdgopt", &["-p", "-0.35", "-1", "--help"]).expect("accepted");
    assert_eq!(o.get_all("extraparam"), ["-0.35", "-1"]);
    // but an option where a value belongs is still "expected N arguments"
    let e = parse("callpeak", &["-m", "-3", "-9x"]).expect_err("rejected");
    assert!(e.contains("expected 2 arguments"), "{e}");
}

#[test]
fn a_negative_number_of_the_wrong_type_is_still_an_invalid_value() {
    // `--shift` is `type=int`; argparse accepts `-0.35` as the *value* and then
    // fails the conversion, which is a different message from "unrecognized".
    let e = parse("callpeak", &["--shift", "-0.35"]).expect_err("rejected");
    assert!(e.contains("invalid int value: '-0.35'"), "{e}");
}

#[test]
fn the_parsed_value_is_the_negative_number() {
    let o = parse("bdgcmp", &["-p", "-1", "--help"]).expect("accepted");
    assert_eq!(o.get("pseudocount"), Some("-1"));
    let o = parse("bdgopt", &["--extra-param=-0.35", "--help"]).expect("accepted");
    assert_eq!(o.get("extraparam"), Some("-0.35"));
    let o = parse("callpeak", &["--shift", "-75", "--help"]).expect("accepted");
    assert_eq!(o.get("shift"), Some("-75"));
}

#[test]
fn a_token_outside_the_negative_number_pattern_is_still_a_usage_error() {
    // `-1.` (no digit after the dot), `-1e-3` (no exponent branch) and `-9x` match
    // neither alternative of `^-\d+$|^-\d*\.\d+$`
    //
    // `-i`/`-o` are given here because `parse_args` reports the leftovers *after*
    // `_parse_known_args` has checked the required actions (`argparse.py:2185-2210`
    // then `argparse.py:1913-1917`), and the reference says the same: bare
    // `macs3 bdgopt -p -9x` is "the following arguments are required", while
    // `macs3 bdgopt -i in.bdg -o o.bdg -p -9x` is
    // "unrecognized arguments: -9x". Both exit 2.
    for bad in ["-9x", "-1e-3", "-1.", "-.e3", "--9x"] {
        let e = parse("bdgopt", &["-i", "in.bdg", "-o", "o.bdg", "-p", bad])
            .expect_err("must be rejected");
        assert_eq!(e, format!("error: unrecognized arguments: {bad}"));
    }
    for bad in ["-9x", "-1e-3", "-1.", "-.e3", "--9x"] {
        let e = parse("bdgopt", &["-p", bad]).expect_err("must be rejected");
        assert!(
            e.contains("unrecognized") || e.contains("expected") || e.contains("required"),
            "{bad}: {e}"
        );
    }
    // a bare unknown flag is still the "unrecognized arguments" usage error
    let e = parse("bdgopt", &["-i", "in.bdg", "-o", "o.bdg", "-9x"]).expect_err("must be rejected");
    assert_eq!(e, "error: unrecognized arguments: -9x");
}
