//! Regression: what a literal `--` does to the token stream.
//!
//! `_parse_known_args` (`argparse.py:1969-1976`) turns a literal `--` into a `-` in
//! the internal pattern and marks **every** later token `A`, whatever it looks like:
//! no later token can ever be an option, and no later token is an unknown option.
//!
//! The trap is that `_get_nargs_pattern` (`argparse.py:2411-2413`) strips the `-`
//! characters out of the pattern for an action that has option strings, so the `+`
//! and `*` patterns become `(A[A]*)` and `([A-]*)`, and the single-argument pattern
//! becomes `(A)`. **None of them can match a `-`.** No MACS3 subparser declares a
//! positional either, so `consume_positionals` can never take a `A`.
//!
//! The consequences, all verified against `.oracle/venv/bin/macs3` rather than assumed:
//!
//! * `-t -- -x.bed` is **not** `-t -x.bed`. It is
//!   `argument -t/--treatment: expected at least one argument`, because the pattern
//!   right after `-t` starts with `-` and `(A[A]*)` does not match it.
//! * everything from `--` onwards is left over, `--` itself included, and is reported
//!   last -- after every value conversion and after the required checks -- as
//!   `unrecognized arguments: <tokens joined by spaces>`.
//!
//! ```text
//! callpeak -t -- -x.bed              2  macs3 callpeak: error: argument -t/--treatment: expected at least one argument
//! callpeak -t --                     2  macs3 callpeak: error: argument -t/--treatment: expected at least one argument
//! callpeak -t a b -- -c              2  macs3: error: unrecognized arguments: -- -c
//! callpeak -t a.bed --               2  macs3: error: unrecognized arguments: --
//! callpeak -t a.bed -- --            2  macs3: error: unrecognized arguments: -- --
//! callpeak -t a.bed -- -q 0.05       2  macs3: error: unrecognized arguments: -- -q 0.05
//! callpeak -t a.bed -- --nomodel     2  macs3: error: unrecognized arguments: -- --nomodel
//! bdgpeakcall -i in.bdg -o -- out.bed 2  macs3 bdgpeakcall: error: argument -o/--ofile: expected one argument
//! bdgopt -i in.bdg -o o.bdg -m multiply -p 2 --x   2  macs3: error: unrecognized arguments: --x
//! ```
//!
//! This port used to `break` out of its loop at the first `--`, which *silently
//! accepted* the whole tail: `callpeak -t a.bed -- -q 0.05` ran to completion and
//! exited 0 where upstream exits 2.

use macs_cli::{parse_flags, Options};

fn parse(sub: &str, words: &[&str]) -> Result<Options, String> {
    let args: Vec<String> = words.iter().map(std::string::ToString::to_string).collect();
    parse_flags(sub, &args).map_err(|e| e.0)
}

/// `(subcommand, argv, exit status, last stderr line)` straight from the reference.
type Case = (&'static str, &'static [&'static str], i32, &'static str);

const CASES: &[Case] = &[
    (
        "callpeak",
        &["-t", "--", "-x.bed"],
        2,
        "argument -t/--treatment: expected at least one argument",
    ),
    (
        "callpeak",
        &["-t", "--"],
        2,
        "argument -t/--treatment: expected at least one argument",
    ),
    (
        "callpeak",
        &["-t", "a", "b", "--", "-c"],
        2,
        "unrecognized arguments: -- -c",
    ),
    (
        "callpeak",
        &["-t", "a.bed", "--"],
        2,
        "unrecognized arguments: --",
    ),
    (
        "callpeak",
        &["-t", "a.bed", "--", "--"],
        2,
        "unrecognized arguments: -- --",
    ),
    (
        "callpeak",
        &["-t", "a.bed", "--", "-q", "0.05"],
        2,
        "unrecognized arguments: -- -q 0.05",
    ),
    (
        "callpeak",
        &["-t", "a.bed", "--", "--nomodel"],
        2,
        "unrecognized arguments: -- --nomodel",
    ),
    (
        "bdgpeakcall",
        &["-i", "in.bdg", "-o", "--", "out.bed"],
        2,
        "argument -o/--ofile: expected one argument",
    ),
    (
        "bdgopt",
        &[
            "-i", "in.bdg", "-o", "o.bdg", "-m", "multiply", "-p", "2", "--x",
        ],
        2,
        "unrecognized arguments: --x",
    ),
    (
        "bdgbroadcall",
        &["-i", "in.bdg", "-o", "o.bed", "--", "-c", "2"],
        2,
        "unrecognized arguments: -- -c 2",
    ),
    (
        "filterdup",
        &["-i", "in.bed", "-f", "BED", "-o", "o.bed", "--", "-d"],
        2,
        "unrecognized arguments: -- -d",
    ),
];

#[test]
fn the_separator_matches_the_reference_for_every_case() {
    for &(sub, argv, status, message) in CASES {
        let r = parse(sub, argv);
        // every recorded case is a usage error; the exit status is asserted by the
        // harness in `target/agent-scratch/compare.py`, which also runs the reference
        assert_eq!(
            status, 2,
            "{sub} {argv:?} is not a usage error in the table"
        );
        let Err(e) = r else {
            panic!("{sub} {} was accepted", argv.join(" "))
        };
        assert!(
            e.ends_with(message),
            "{sub} {}: got {e:?}, reference says {message:?}",
            argv.join(" ")
        );
    }
}

#[test]
fn the_separator_is_reported_with_the_whole_tail() {
    // argparse's `parse_args` (`argparse.py:1913-1917`) joins the leftovers with a
    // space, in order, `--` first -- not one message per token
    let Err(e) = parse("callpeak", &["-t", "a", "b", "--", "-c", "-d", "extra"]) else {
        panic!("accepted")
    };
    assert_eq!(e, "error: unrecognized arguments: -- -c -d extra");
}

#[test]
fn a_value_position_filled_by_the_separator_is_an_argument_count_error() {
    // `(A[A]*)` cannot match a leading `-`, so `nargs='+'` reports "at least one"
    // and `nargs=None` reports "expected one argument" -- it never silently reads
    // the token after `--` as the value
    for (sub, argv, msg) in [
        (
            "callpeak",
            ["-t", "--", "x.bed"],
            "error: argument -t/--treatment: expected at least one argument",
        ),
        (
            "bdgopt",
            ["-i", "--", "in.bdg"],
            "error: argument -i/--ifile: expected one argument",
        ),
        (
            "bdgcmp",
            ["-t", "--", "a.bdg"],
            "error: argument -t/--tfile: expected one argument",
        ),
        // `filterdup -i` and `hmmratac -i` are `nargs='+'`, so the same shape reads
        // "expected at least one argument" instead
        (
            "filterdup",
            ["-i", "--", "in.bed"],
            "error: argument -i/--ifile: expected at least one argument",
        ),
        (
            "hmmratac",
            ["-i", "--", "in.bed"],
            "error: argument -i/--input: expected at least one argument",
        ),
    ] {
        let Err(e) = parse(sub, &argv) else {
            panic!("{sub} {} was accepted", argv.join(" "))
        };
        assert_eq!(e, msg, "{sub} {}", argv.join(" "));
    }
}

#[test]
fn a_required_group_beats_the_leftover_report() {
    // `bdgbroadcall -i in.bdg --nope` is "one of the arguments ... is required",
    // not "unrecognized arguments": the group check is inside `_parse_known_args`
    let Err(e) = parse("bdgbroadcall", &["-i", "in.bdg", "--nope"]) else {
        panic!("accepted")
    };
    assert_eq!(
        e,
        "error: one of the arguments -o/--ofile --o-prefix is required"
    );
    // and a bad value beats both, because it is raised while the option is consumed
    let Err(e) = parse("bdgopt", &["--nope", "-m", "bogus"]) else {
        panic!("accepted")
    };
    assert_eq!(
        e,
        "error: argument -m/--method: invalid choice: 'bogus' \
         (choose from 'multiply', 'add', 'p2q', 'max', 'min')"
    );
}

#[test]
fn a_separator_before_any_option_leaves_everything_as_an_extra() {
    // `-- -i x.bdg` makes `-i` an argument, so the required check fires first:
    // argparse checks required actions *before* `parse_args` reports extras
    let Err(e) = parse("bdgbroadcall", &["--", "-i", "in.bdg", "-o", "o.bed"]) else {
        panic!("accepted")
    };
    assert_eq!(e, "error: the following arguments are required: -i/--ifile");

    let Err(e) = parse("bdgopt", &["--", "-i", "in.bdg", "-o", "o.bdg"]) else {
        panic!("accepted")
    };
    assert_eq!(
        e,
        "error: the following arguments are required: -i/--ifile, -o/--ofile"
    );
}

#[test]
fn a_value_error_still_beats_the_leftover_report() {
    // the conversion happens while the option is consumed, the extras are reported by
    // `parse_args` afterwards
    let Err(e) = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o.bdg", "-p", "zz", "--x"],
    ) else {
        panic!("accepted")
    };
    assert_eq!(
        e,
        "error: argument -p/--extra-param: invalid float value: 'zz'"
    );
}

#[test]
fn a_flag_before_the_separator_still_takes_its_ordinary_value() {
    // the separator only affects tokens at or after its own position
    let o = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o.bdg", "-m", "max", "-p", "3"],
    )
    .expect("parsed");
    assert_eq!(o.get("method"), Some("max"));
    assert_eq!(o.float("extraparam"), Some(3.0));
}

#[test]
fn the_first_separator_ends_the_options() {
    // only the first `--` is a separator; the rest are ordinary tokens (and this
    // subcommand declares no positionals, so they are all extras)
    let Err(e) = parse(
        "bdgopt",
        &["-i", "in.bdg", "-o", "o.bdg", "--", "--", "-p", "2"],
    ) else {
        panic!("accepted")
    };
    assert_eq!(e, "error: unrecognized arguments: -- -- -p 2");
}

#[test]
fn help_still_wins_over_the_leftover_report() {
    // `-h` exits from inside `take_action`, so a `--` and its tail are never reported
    let o = parse("bdgopt", &["-h", "--", "x"]).expect("help short-circuits");
    assert!(o.help);
}

#[test]
fn the_argument_vector_is_still_echoed_verbatim() {
    // `sys.argv` is what upstream writes into the `# Command line:` header, so the
    // separator and its tail have to survive parsing untouched
    // `-h` first: it exits from inside `take_action`, so the parse succeeds and the
    // whole vector -- separator and tail -- is still observable
    let o = parse_flags(
        "callpeak",
        &[
            "-h".into(),
            "-t".into(),
            "a.bed".into(),
            "--".into(),
            "-q".into(),
            "0.05".into(),
        ],
    )
    .expect("parses");
    assert_eq!(
        o.raw_argv,
        ["callpeak", "-h", "-t", "a.bed", "--", "-q", "0.05"]
    );
}
