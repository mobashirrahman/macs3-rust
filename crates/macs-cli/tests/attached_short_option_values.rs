//! Regression: a **short option with its value glued on**, and a short
//! `store_true` option followed by more characters in the same token.
//!
//! This is `_get_option_tuples` (`argparse.py:2341-2352`) plus the "identify
//! additional optionals in the same arg string" loop of `consume_optional`
//! (`argparse.py:2022-2063`). The rule is:
//!
//! * a token with **one** prefix character is first read as its first two
//!   characters, so `-p0.01` is `-p 0.01` and `-p-1` is `-p -1`;
//! * if that short option takes **no** arguments, argparse peels one character
//!   off the tail at a time and looks the next `-x` up, so `-Bp0.05` is
//!   `-B -p 0.05` and `-Bq0.05` is `-B -q 0.05`;
//! * if the peeled `-x` is not an option, the tail is left behind verbatim as an
//!   extra, so `-Bx=1` becomes `-B` plus `unrecognized arguments: -x=1`;
//! * a token with **two** prefix characters is *never* split other than at `=`,
//!   so `--tsize1_0` is not `--tsize 1_0` but an unknown option.
//!
//! Reference bytes (`OPENBLAS_CORETYPE=Haswell .oracle/venv/bin/macs3`, exit
//! status first, then its last stderr line):
//!
//! ```text
//! bdgopt -i in.bdg -o out.bdg -m multiply -p-1      0
//! bdgopt -i in.bdg -o out.bdg -m multiply -p0.01    0
//! bdgopt -i in.bdg -o out.bdg -m multiply -p 0.01   0
//! callpeak -t x.bed -B -q 0.05                     0
//! callpeak -t x.bed -Bq0.05                        0   (-B then -q 0.05)
//! callpeak -t x.bed -Bp0.05                        0   (-B then -p 0.05)
//! callpeak -t x.bed -Bq                            0
//! callpeak -t x.bed -ttreat.bed                    0
//! callpeak -t x.bed -s50                           0
//! callpeak -t x.bed -s-50                          0
//! callpeak -t x.bed -g1e6                          0
//! callpeak -t x.bed -nx                            0
//! bdgopt -i in.bdg -o out.bdg --pvalue=2           2  macs3: error: unrecognized arguments: --pvalue=2
//! callpeak -t x.bed --tsize1_0                     2  macs3: error: unrecognized arguments: --tsize1_0
//! callpeak -t x.bed -B=1                           2  macs3 callpeak: error: argument -B/--bdg: ignored explicit argument '1'
//! callpeak -t x.bed --nomodel=1                    2  macs3 callpeak: error: argument --nomodel: ignored explicit argument '1'
//! callpeak -t x.bed --call-summits=1               2  macs3 callpeak: error: argument --call-summits: ignored explicit argument '1'
//! callpeak -t x.bed -m320                          2  macs3 callpeak: error: argument -m/--mfold: expected 2 arguments
//! ```
//!
//! The last four are the shapes that used to be **accepted** here: `--nomodel=1`
//! and `--call-summits=1` are `store_true` options, and argparse refuses the value
//! outright (`argparse.py:2084-2087`) rather than ignoring it.

use macs_cli::{parse_flags, Options};

fn parse(sub: &str, words: &[&str]) -> Result<Options, String> {
    let args: Vec<String> = words.iter().map(std::string::ToString::to_string).collect();
    parse_flags(sub, &args).map_err(|e| e.0)
}

#[test]
fn a_value_glued_to_a_short_option_is_the_value() {
    // `-p` is `nargs='*'` here, the widest numeric case
    for glued in ["-p0.01", "-p0.5", "-p1e3"] {
        let o = parse(
            "bdgopt",
            &["-i", "in.bdg", "-o", "o.bdg", "-m", "multiply", glued],
        )
        .unwrap_or_else(|e| panic!("{glued}: {e}"));
        assert_eq!(o.get_all("extraparam").len(), 1, "{glued}");
        let want = &glued[2..];
        assert_eq!(o.get("extraparam"), Some(want), "{glued}");
    }
}

#[test]
fn a_negative_number_glued_to_a_short_option_is_still_a_value() {
    // `-p-1` is `-p -1`: `_get_option_tuples` puts the tail in `explicit_arg`, and
    // the negative-number rule of `_parse_optional` is never reached, because the
    // token was already claimed as an option.
    let o = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o.bdg", "-m", "add", "-p-1"],
    )
    .unwrap();
    assert_eq!(o.get_all("extraparam"), ["-1"]);
    assert_eq!(o.float("extraparam"), Some(-1.0));

    let o = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o.bdg", "-m", "add", "-p-0.35"],
    )
    .unwrap();
    assert_eq!(o.get_all("extraparam"), ["-0.35"]);
    assert_eq!(o.float("extraparam"), Some(-0.35));
}

#[test]
fn a_glued_value_equals_the_separate_spelling() {
    // the parsed option must be *identical*, not merely accepted
    let sep = parse(
        "bdgopt",
        &[
            "-i", "in.bdg", "-o", "o.bdg", "-m", "multiply", "-p", "0.01",
        ],
    )
    .unwrap();
    let glued = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o.bdg", "-m", "multiply", "-p0.01"],
    )
    .unwrap();
    assert_eq!(sep.get("extraparam"), glued.get("extraparam"));
    assert_eq!(sep.get("method"), glued.get("method"));
}

#[test]
fn a_short_option_with_no_arguments_splits_into_two_options() {
    let o = parse(
        "callpeak",
        &["-t", "a.bed", "-f", "BED", "--nomodel", "-Bq0.05", "--help"],
    )
    .unwrap();
    assert!(o.flag("store_bdg"), "-B did not fire");
    assert_eq!(o.get("qvalue"), Some("0.05"));
    assert!(o.help);

    let o = parse(
        "callpeak",
        &["-t", "a.bed", "-f", "BED", "--nomodel", "-Bp0.05", "--help"],
    )
    .unwrap();
    assert!(o.flag("store_bdg"), "-B did not fire");
    assert_eq!(o.get("pvalue"), Some("0.05"));

    // the tail may run out exactly, and then the *second* option still wants its
    // own value: `-Bq` is `-B` then `-q`, which is "expected one argument"
    let Err(e) = parse("callpeak", &["-t", "a.bed", "-f", "BED", "-Bq"]) else {
        panic!("-Bq was accepted")
    };
    assert_eq!(e, "error: argument -q/--qvalue: expected one argument");
}

#[test]
fn a_split_stops_at_a_letter_that_is_not_an_option() {
    // `callpeak` has `-B` but no `-x`, so `-Bx=1` runs `-B` and leaves the rest
    let Err(e) = parse("callpeak", &["-t", "a.bed", "-f", "BED", "-Bx=1"]) else {
        panic!("-Bx=1 was accepted")
    };
    assert_eq!(
        e, "error: unrecognized arguments: -x=1",
        "the tail is reported verbatim, not re-parsed"
    );
    // `bdgopt` has no `-B` at all, so the whole token is one unknown option
    let Err(e) = parse("bdgopt", &["-i", "in.bdg", "-o", "o.bdg", "-Bx=1"]) else {
        panic!("bdgopt -Bx=1 was accepted")
    };
    assert_eq!(e, "error: unrecognized arguments: -Bx=1");
}

#[test]
fn a_short_filename_is_glued_to_its_option() {
    let o = parse(
        "callpeak",
        &[
            "-ttreat.bed",
            "-f",
            "BED",
            "--nomodel",
            "--extsize",
            "50",
            "--help",
        ],
    )
    .unwrap();
    assert_eq!(o.get("tfile"), Some("treat.bed"));

    let o = parse("callpeak", &["-ttreat.bed", "-fBED", "--help"]).unwrap();
    assert_eq!(o.get("tfile"), Some("treat.bed"));
    assert_eq!(o.get("format"), Some("BED"));
}

#[test]
fn numeric_short_options_accept_a_glued_value() {
    for (words, dest, want) in [
        (vec!["-t", "a.bed", "-s50", "--help"], "tsize", "50"),
        (vec!["-t", "a.bed", "-s-50", "--help"], "tsize", "-50"),
        (vec!["-t", "a.bed", "-g1e6", "--help"], "gsize", "1e6"),
        (vec!["-t", "a.bed", "-nx", "--help"], "name", "x"),
    ] {
        let o = parse("callpeak", &words).unwrap_or_else(|e| panic!("{words:?}: {e}"));
        assert_eq!(o.get(dest), Some(want), "{words:?}");
    }
}

#[test]
fn a_long_option_never_splits_on_a_glued_value() {
    // `--tsize1_0`: two prefix characters, so the only split is at `=`
    let e = parse("callpeak", &["-t", "a.bed", "--tsize1_0"]).unwrap_err();
    assert_eq!(e, "error: unrecognized arguments: --tsize1_0");

    // `bdgopt`'s flag is `--extra-param`, so `--pvalue` is not an abbreviation
    let e = parse("bdgopt", &["-i", "in.bdg", "-o", "o.bdg", "--pvalue=2"]).unwrap_err();
    assert_eq!(e, "error: unrecognized arguments: --pvalue=2");
}

#[test]
fn an_explicit_argument_on_a_zero_argument_option_is_refused() {
    for (sub, flag, msg) in [
        (
            "callpeak",
            "-B=1",
            "error: argument -B/--bdg: ignored explicit argument '1'",
        ),
        (
            "callpeak",
            "--nomodel=1",
            "error: argument --nomodel: ignored explicit argument '1'",
        ),
        (
            "callpeak",
            "--call-summits=1",
            "error: argument --call-summits: ignored explicit argument '1'",
        ),
        (
            "bdgpeakcall",
            "--call-summits=1",
            "error: argument --call-summits: ignored explicit argument '1'",
        ),
        (
            "bdgpeakcall",
            "--no-trackline=1",
            "error: argument --no-trackline: ignored explicit argument '1'",
        ),
    ] {
        let Err(e) = parse(sub, &[flag]) else {
            panic!("{flag} was accepted")
        };
        assert_eq!(e, msg, "{sub} {flag}");
    }
}

#[test]
fn a_glued_value_cannot_satisfy_a_multi_argument_option() {
    let e = parse("callpeak", &["-t", "a.bed", "-m320"]).unwrap_err();
    assert_eq!(e, "error: argument -m/--mfold: expected 2 arguments");
}

#[test]
fn a_glued_value_runs_the_whole_pipeline() {
    // the accept/reject decision is not the whole point: the *value* has to reach
    // the pipeline, so `bdgopt -m multiply -p2` writes the same file as
    // `-m multiply -p 2`
    let dir = std::env::temp_dir().join(format!("macs3rs-glued-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("in.bdg");
    std::fs::write(&input, "chr1\t0\t10\t5\nchr1\t10\t20\t7\n").unwrap();
    let out = dir.join("o.bdg");

    let args: Vec<String> = [
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-m",
        "multiply",
        "-p2",
    ]
    .iter()
    .map(std::string::ToString::to_string)
    .collect();
    let o = parse_flags("bdgopt", &args).expect("bdgopt parses");
    macs_cli::commands::bedgraph_cmds::bdgopt(&o).expect("bdgopt runs");
    let got = std::fs::read_to_string(&out).expect("output written");
    // upstream's bytes for `macs3 bdgopt -m multiply -p 2`
    assert_eq!(
        got,
        concat!(
            "track type=bedGraph name=\"MULTIPLY_modified_scores\" description=\"Scores calculated by MULTIPLY\" visibility=2 alwaysZero=on\n",
            "chr1\t0\t10\t10.00000\n",
            "chr1\t10\t20\t14.00000\n",
        )
    );
    let _ = std::fs::remove_dir_all(&dir);
}
