//! `randsample`'s percentage/number range checks, and the tag-size limit that
//! `filterdup` and `randsample` inherit from Cython.
//!
//! **Range checks.** `opt_validate_randsample` (`OptValidator.py:375-383`) is
//!
//! ```python
//! if options.percentage:
//!     if options.percentage > 100.0:
//!         logger.error("Percentage can't be bigger than 100.0. ...")
//!         sys.exit(1)
//! elif options.number:
//!     if options.number <= 0:
//!         logger.error("Number of tags can't be smaller than or equal to 0. ...")
//!         sys.exit(1)
//! ```
//!
//! and `-p`/`-n` have **no default** (`bin/macs3:485-488`), so both are `None`
//! unless given. That makes the two tests independent truthiness checks, with two
//! consequences this port got wrong: `-p 1000` was accepted and the tags sampled
//! anyway, and `-p 0` was rejected where upstream succeeds. `-n 0` is *not*
//! rejected here either (`elif options.number:` is false); it dies later, in
//! `randsample_cmd.py:67`, formatting the still-`None` percentage.
//!
//! **Tag size.** Both commands assign the flag straight onto the Cython track --
//! `inputtrack.fw = options.tsize` (`filterdup_cmd.py:58`) and
//! `treat.fw = options.tsize` (`randsample_cmd.py:55`) -- and `FWTrack.fw` is
//! `cython.declare(cython.int)`, so a value above `2**31 - 1` raises
//! `OverflowError: value too large to convert to int`. A non-positive value is a
//! legal C `int` and only trips `assert self.fw > 0` in `print_to_bed`
//! (`FixWidthTrack.py:509`). In both commands the output file is opened *first*
//! (`filterdup_cmd.py:40-43`, `randsample_cmd.py:41-44`), so either failure leaves
//! an **empty output file** behind -- which is why this port's refusal, which
//! created nothing, was also a divergence.
//!
//! Every exit status, file list and byte count here was measured against the pinned
//! oracle (`OPENBLAS_CORETYPE=Haswell .oracle/venv/bin/macs3 ...`); nothing in this
//! file needs Python or the reference at run time.

use std::path::{Path, PathBuf};
use std::process::Command;

const RUST: &str = env!("CARGO_BIN_EXE_macs3-rs");

fn scratch(tag: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/it-randsample-range")
        .join(tag);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn input(dir: &Path) -> String {
    let p = dir.join("in.bed");
    let mut s = String::new();
    for i in 0..10u64 {
        s.push_str(&format!(
            "chr1\t{}\t{}\tr\t0\t+\n",
            100 + i * 7,
            130 + i * 7
        ));
        s.push_str(&format!(
            "chr1\t{}\t{}\tr\t0\t-\n",
            300 + i * 11,
            330 + i * 11
        ));
    }
    std::fs::write(&p, s).unwrap();
    p.display().to_string()
}

struct Run {
    code: i32,
    files: Vec<String>,
    sizes: Vec<(String, u64)>,
}

/// Run one of the two tag-filtering commands inside a scratch directory.
///
/// `-o` is absolute because both commands join it onto `--outdir`
/// (`os.path.join(options.outdir, options.outputfile)`, `filterdup_cmd.py:41`),
/// so a bare `out.bed` would land in `<dir>/<dir>/out.bed`.
fn run(tag: &str, cmd: &str, extra: &[&str]) -> Run {
    let dir = scratch(tag);
    let inp = input(&dir);
    let mut argv: Vec<String> = vec![
        cmd.into(),
        "-i".into(),
        inp,
        "-f".into(),
        "BED".into(),
        "-o".into(),
        dir.join("out.bed").display().to_string(),
        "--outdir".into(),
        dir.display().to_string(),
    ];
    argv.extend(extra.iter().map(|s| (*s).to_string()));
    if cmd == "randsample"
        && !argv
            .iter()
            .any(|a| a == "-p" || a == "-n" || a.starts_with("-p="))
    {
        // `-p`/`-n` are a *required* mutually-exclusive group (`bin/macs3:483`)
        argv.push("-p".into());
        argv.push("50".into());
    }
    let st = Command::new(RUST).args(&argv).output().unwrap();
    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n != "in.bed")
        .collect();
    files.sort();
    let sizes = files
        .iter()
        .map(|f| {
            (
                f.clone(),
                std::fs::metadata(dir.join(f)).map(|m| m.len()).unwrap_or(0),
            )
        })
        .collect();
    Run {
        code: st.status.code().unwrap_or(-1),
        files,
        sizes,
    }
}

#[test]
fn a_percentage_above_one_hundred_is_refused_before_the_writer_opens() {
    let r = run("p1000", "randsample", &["-p", "1000", "-s", "30"]);
    assert_eq!(r.code, 1);
    assert_eq!(r.files, Vec::<String>::new());
    let r = run("p1e30", "randsample", &["-p", "1e30", "-s", "30"]);
    assert_eq!(r.code, 1);
    assert_eq!(r.files, Vec::<String>::new());
}

#[test]
fn a_negative_percentage_keeps_the_empty_output_file() {
    // `sample_percent` allocates `int(total * percentage)`, a negative dimension
    // (`FixWidthTrack.py:449`) -- but only after the writer has been opened.
    let r = run("pneg", "randsample", &["-p", "-1", "-s", "30"]);
    assert_eq!(r.code, 1);
    assert_eq!(r.files, vec!["out.bed".to_string()]);
    assert_eq!(r.sizes, vec![("out.bed".to_string(), 0)]);
}

#[test]
fn a_negative_number_is_refused_and_a_zero_number_is_not() {
    // `-n -5` trips `elif options.number: ... <= 0`, before the writer exists.
    let neg = run("nneg", "randsample", &["-n", "-5", "-s", "30"]);
    assert_eq!(neg.code, 1);
    assert_eq!(neg.files, Vec::<String>::new());
    // `-n 0` slips past that check and dies formatting the unset percentage, after
    // the file has been created.
    let zero = run("nzero", "randsample", &["-n", "0", "-s", "30"]);
    assert_eq!(zero.code, 1);
    assert_eq!(zero.files, vec!["out.bed".to_string()]);
    assert_eq!(zero.sizes, vec![("out.bed".to_string(), 0)]);
}

#[test]
fn asking_for_more_tags_than_the_file_has_leaves_an_empty_file() {
    let r = run("n1e9", "randsample", &["-n", "1000000000", "-s", "30"]);
    assert_eq!(r.code, 1);
    assert_eq!(r.files, vec!["out.bed".to_string()]);
    assert_eq!(r.sizes, vec![("out.bed".to_string(), 0)]);
}

#[test]
fn a_percentage_of_zero_is_a_valid_run_that_keeps_nothing() {
    // `-p 0` is falsy for the range check but still *set*, so `%.2f` formats it and
    // `int(total * 0.0)` is 0: exit 0, every tag dropped.
    let r = run("p0", "randsample", &["-p", "0", "-s", "30"]);
    assert_eq!(r.code, 0);
    assert_eq!(r.files, vec!["out.bed".to_string()]);
    assert_eq!(r.sizes, vec![("out.bed".to_string(), 0)]);
}

#[test]
fn a_percentage_of_one_hundred_keeps_everything() {
    let r = run("p100", "randsample", &["-p", "100", "-s", "30"]);
    assert_eq!(r.code, 0);
    assert_eq!(r.files, vec!["out.bed".to_string()]);
    let bytes = r.sizes[0].1;
    assert!(
        bytes > 0,
        "an empty file would mean the sampler dropped everything"
    );
}

#[test]
fn a_tag_size_above_two_to_the_thirty_one_is_an_overflow() {
    for cmd in ["randsample", "filterdup"] {
        for (tag, extra) in [
            ("big", vec!["-s", "3000000000"]),
            ("huge", vec!["-s", "9223372036854775808"]),
        ] {
            let r = run(&format!("{cmd}-{tag}"), cmd, &extra);
            assert_eq!(r.code, 1, "{cmd} {extra:?} exited {}", r.code);
            // the writer opened at step 0, so the file is there and empty
            assert_eq!(
                r.files,
                vec!["out.bed".to_string()],
                "{cmd} {extra:?} left {:?}",
                r.files
            );
            assert_eq!(r.sizes, vec![("out.bed".to_string(), 0)]);
        }
    }
}

#[test]
fn a_negative_tag_size_is_rejected_by_print_to_bed_not_by_the_validator() {
    for cmd in ["randsample", "filterdup"] {
        for v in ["-50", "-5", "-1"] {
            // `-s -50` is a usage error (argparse reads `-50` as a flag), so the
            // attached spelling is the only way to pass a negative tag size.
            let r = run(&format!("{cmd}-neg"), cmd, &[&format!("-s{v}")]);
            assert_eq!(r.code, 1, "{cmd} -s{v} exited {}", r.code);
            assert_eq!(r.files, vec!["out.bed".to_string()]);
            assert_eq!(r.sizes, vec![("out.bed".to_string(), 0)]);
        }
    }
}

#[test]
fn the_largest_representable_tag_size_is_still_accepted() {
    // 2**31 - 1 fits in a C `int`, so this is a plain (very slow, very wide) run
    // rather than an OverflowError. Only the exit status is asserted: the pileup
    // would be gigabases wide.
    let r = run("max-int", "filterdup", &["-s", "2147483647", "-d"]);
    assert_eq!(r.code, 0);
}

#[test]
fn filterdup_still_opens_its_output_file_for_a_dry_run() {
    // `outfhd = open(...)` is step 0, before the `if not options.dryrun:` branch at
    // the end of `filterdup_cmd.run`, so `--dry-run -o out.bed` still leaves the
    // empty file behind.
    let r = run("dryrun", "filterdup", &["-d"]);
    assert_eq!(r.code, 0);
    assert_eq!(r.files, vec!["out.bed".to_string()]);
    assert_eq!(r.sizes, vec![("out.bed".to_string(), 0)]);
}
