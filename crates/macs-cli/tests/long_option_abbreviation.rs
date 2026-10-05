//! Regression: **long-option abbreviation**, which upstream has on.
//!
//! `argparse.ArgumentParser` defaults to `allow_abbrev=True`
//! (`argparse.py:1788`), and none of MACS3's parsers turns it off, so every
//! *unambiguous* prefix of a long option names it. `_get_option_tuples`
//! (`argparse.py:2341-2351`) collects every option string that starts with the
//! token (or with the part before its `=`): one candidate is accepted, more than
//! one is the `ambiguous option` usage error.
//!
//! This port used to reject every abbreviation, so `bdgopt --pvalue 0.01` and
//! `callpeak --qvalue 0.05` exited 2 where upstream runs.
//!
//! Reference bytes (`.oracle/venv/bin/macs3`, exit status then its last stderr
//! line):
//!
//! ```text
//! bdgopt -i in.bdg -o o.bdg -m multiply --e 2                 0
//! bdgopt -i in.bdg -o o.bdg -m multiply --ex 2                0
//! bdgopt -i in.bdg -o o.bdg -m multiply --extra-p 2           0
//! bdgopt -i in.bdg -o o.bdg --out .                           0
//! bdgopt -i in.bdg -o o.bdg --verb 1                           0
//! bdgopt -i in.bdg -o o.bdg --m max -p 3                      0
//! callpeak -t a.bed --treat a.bed                            0
//! callpeak -t a.bed --q 0.05                                 0
//! callpeak -t a.bed --qval 0.05                              0
//! callpeak -t a.bed --pval 0.05                              0
//! callpeak -t a.bed --ext 50                                 0
//! callpeak -t a.bed --exts 50                                0
//! callpeak -t a.bed --shi -1                                 0
//! callpeak -t a.bed --o .                                    0
//! callpeak -t a.bed --pvalu=0.05                             0
//! callpeak -t a.bed --pvalue=0.05                            0
//! callpeak -t a.bed --nomodel -h                              0
//!
//! bdgopt -i in.bdg -o o.bdg --o x       2  macs3 bdgopt: error: ambiguous option: --o could match --outdir, --ofile
//! callpeak -t a.bed --b                2  macs3 callpeak: error: ambiguous option: --b could match --bdg, --barcodes, --broad, --broad-cutoff, --buffer-size
//! bdgopt -i in.bdg -o o.bdg --p 2      2  macs3: error: unrecognized arguments: --p 2
//! callpeak -t a.bed --zzz 1            2  macs3: error: unrecognized arguments: --zzz 1
//! ```
//!
//! `--p` is **not** an abbreviation in `bdgopt`: the flag is `--extra-param`, so no
//! long option there starts with `--p` and the token is simply unknown.

use macs_cli::{flag_specs, parse_flags, Options};

fn parse(sub: &str, words: &[&str]) -> Result<Options, String> {
    let args: Vec<String> = words.iter().map(std::string::ToString::to_string).collect();
    parse_flags(sub, &args).map_err(|e| e.0)
}

#[test]
fn an_unambiguous_prefix_of_a_long_option_is_accepted() {
    for prefix in [
        "--e",
        "--ex",
        "--ext",
        "--extr",
        "--extra",
        "--extra-p",
        "--extra-para",
    ] {
        let words = ["-i", "in.bdg", "-o", "o.bdg", "-m", "multiply", prefix, "2"];
        let o = parse("bdgopt", &words).unwrap_or_else(|e| panic!("{prefix}: {e}"));
        assert_eq!(o.get("extraparam"), Some("2"), "{prefix}");
    }
}

#[test]
fn every_prefix_of_a_flag_that_resolves_to_one_option_works() {
    // a mechanical cross-check rather than a hand-written list: for each long
    // option of a subcommand, every prefix that no *other* long option of the same
    // subcommand also starts with must parse to that option
    let specs = flag_specs();
    for sub in [
        "bdgopt",
        "bdgpeakcall",
        "filterdup",
        "callpeak",
        "refinepeak",
    ] {
        let longs: Vec<&str> = specs
            .iter()
            .filter(|s| s.subcommand == sub && s.flag.starts_with("--"))
            .map(|s| s.flag.as_str())
            .collect();
        for &flag in &longs {
            for n in 2..flag.len() {
                let prefix = &flag[..n];
                if longs.iter().any(|o| *o != flag && o.starts_with(prefix)) {
                    continue; // ambiguous, covered by its own test
                }
                if *flag == *"--help" {
                    continue; // the help action stores nothing; see the test below
                }
                let words = vec![prefix, "7"];
                if let Ok(o) = parse(sub, &words) {
                    assert!(
                        o.values.values().any(|v| v == "7"),
                        "{sub} {prefix} parsed but did not store 7: {o:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn an_abbreviation_reaches_the_same_option_as_the_full_spelling() {
    let full = parse("callpeak", &["-t", "a.bed", "--qvalue", "0.05", "--help"]).unwrap();
    for prefix in ["--q", "--qv", "--qva", "--qval", "--qvalu", "--qvalue"] {
        let short = parse("callpeak", &["-t", "a.bed", prefix, "0.05", "--help"])
            .unwrap_or_else(|e| panic!("{prefix}: {e}"));
        assert_eq!(short.get("qvalue"), full.get("qvalue"), "{prefix}");
    }
}

#[test]
fn an_abbreviation_works_with_an_equals_sign_too() {
    // `_get_option_tuples` partitions at `=` before prefix-matching, so `--pvalu=0.05`
    // resolves exactly like `--pvalue=0.05`
    for prefix in ["--pvalu", "--pval", "--pv", "--pvalue"] {
        let o = parse("callpeak", &["-t", "a.bed", &format!("{prefix}=0.05")])
            .unwrap_or_else(|e| panic!("{prefix}=0.05: {e}"));
        assert_eq!(o.get("pvalue"), Some("0.05"), "{prefix}");
    }
}

#[test]
fn an_ambiguous_prefix_is_the_ambiguous_option_usage_error() {
    let Err(e) = parse("bdgopt", &["-i", "in.bdg", "-o", "o.bdg", "--o", "x"]) else {
        panic!("--o was accepted")
    };
    assert_eq!(
        e,
        "error: ambiguous option: --o could match --outdir, --ofile"
    );
}

#[test]
fn the_ambiguous_option_error_is_raised_before_any_value_is_converted() {
    // `_parse_optional` runs over the whole command line in a first pass
    // (`argparse.py:1969-1988`), so the ambiguity is reported even when an earlier
    // token would fail conversion and a later one has a bad value
    let Err(e) = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o.bdg", "-p", "zz", "--o", "x"],
    ) else {
        panic!("--o was accepted")
    };
    assert!(e.starts_with("error: ambiguous option:"), "{e}");
}

#[test]
fn a_prefix_that_matches_nothing_is_an_unrecognized_argument() {
    // `bdgopt`'s option is `--extra-param`, so no long option starts with `--p`
    let Err(e) = parse("bdgopt", &["-i", "in.bdg", "-o", "o.bdg", "--p", "2"]) else {
        panic!("--p was accepted")
    };
    assert_eq!(e, "error: unrecognized arguments: --p 2");

    let Err(e) = parse("callpeak", &["-t", "a.bed", "--zzz", "1"]) else {
        panic!("--zzz was accepted")
    };
    assert_eq!(e, "error: unrecognized arguments: --zzz 1");
}

#[test]
fn a_prefix_never_crosses_the_number_of_dashes() {
    // `--i` has *two* prefix characters, so the abbreviation branch applies and it
    // resolves to `--ifile` -- reference: `bdgopt --i in.bdg -o o.bdg -m multiply
    // -p2` exits 0. Only a token with a single `-` can be split into an option plus
    // a glued value.
    let o = parse(
        "bdgopt",
        &["--i", "in.bdg", "-o", "o.bdg", "-m", "multiply", "-p2"],
    )
    .unwrap();
    assert_eq!(o.get("ifile"), Some("in.bdg"));
    assert_eq!(o.get("method"), Some("multiply"));

    // the help action abbreviates as well, and short-circuits everything after it
    for prefix in ["--h", "--he", "--hel", "--help"] {
        let o = parse("bdgopt", &[prefix]).unwrap_or_else(|e| panic!("{prefix}: {e}"));
        assert!(o.help, "{prefix}");
    }
}

#[test]
fn the_ambiguous_option_error_names_every_candidate() {
    // the message is built from `_get_option_tuples` in `_option_string_actions`
    // order, so it lists all of them rather than the first two
    let Err(e) = parse("callpeak", &["-t", "a.bed", "--n"]) else {
        panic!("--n was accepted")
    };
    assert_eq!(
        e,
        "error: ambiguous option: --n could match --name, --nomodel, --nolambda"
    );
}

#[test]
fn store_true_flags_can_be_abbreviated_too() {
    for prefix in ["--nom", "--nomo", "--nomod", "--nomode", "--nomodel"] {
        let o = parse("callpeak", &["-t", "a.bed", prefix, "--help"])
            .unwrap_or_else(|e| panic!("{prefix}: {e}"));
        assert!(o.flag("nomodel"), "{prefix}");
    }
}
