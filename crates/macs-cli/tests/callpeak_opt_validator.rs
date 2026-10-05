//! The option rejections `MACS3/Utilities/OptValidator.py` performs for
//! `callpeak`, and the two `cython.int` conversions that follow it.
//!
//! `opt_validate_callpeak` runs at `callpeak_cmd.py:51`, before a single alignment
//! byte is read, so every rejection it makes is **exit 1 with no output file at
//! all**. Three of its checks are not range checks in the usual sense and were
//! missed by this port in ways that produced *different results*, not just
//! different messages:
//!
//! * `-p` is tested for **truthiness** (`:127`), not for `None`, so `-p 0` is a
//!   valid run that reports a q-value cutoff -- it does not mean "p-value 0";
//! * `-p`, `-q` and `--broad-cutoff` all reach `math.log` (`:130`, `:132`, `:135`),
//!   so a non-positive value is a `ValueError: math domain error` and exit 1;
//! * `--keep-dup` must satisfy `str.isdigit` (`:104-107`), which rejects `-1` even
//!   though it parses as an integer.
//!
//! The rest are ordinary refusals (`--broad` with `--call-summits` at `:123-125`,
//! `--d-min` at `:141-143`, an inverted `--mfold` at `:146-150`, a negative
//! `--max-count` under `-f FRAG` at `:93-96`), and the last three are Cython's
//! `cython.int` narrowing above `2**31 - 1` at
//! `CallerFromAlignments.call_peaks` (`CallPeakUnit.py:1036-1037`) and
//! `CallerFromAlignments.__init__` (`:440`).
//!
//! Every exit status and "no output file" expectation here was measured against the
//! pinned oracle (`OPENBLAS_CORETYPE=Haswell .oracle/venv/bin/macs3 ...`); nothing
//! in this file needs Python or the reference at run time.

use std::path::PathBuf;
use std::process::Command;

const RUST: &str = env!("CARGO_BIN_EXE_macs3-rs");

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures")
        .join(rel)
}

struct Outcome {
    code: i32,
    /// Sorted names of the files the run left in its `--outdir`.
    files: Vec<String>,
}

/// Run `callpeak` on the two-cluster SE fixture with `extra` appended, and report
/// the exit status plus everything the run left behind.
fn run(extra: &[&str]) -> Outcome {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/it-opt-validator-callpeak")
        .join(extra.join("_").replace(['-', '.', '/'], "_"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let se = fixture("se_basic/gauss_two_peaks");
    let mut args: Vec<String> = vec![
        "callpeak".into(),
        "-t".into(),
        se.join("treat_1.bed").display().to_string(),
        "-c".into(),
        se.join("ctrl.bed").display().to_string(),
        "-f".into(),
        "BED".into(),
        "-g".into(),
        "10000".into(),
        "--nomodel".into(),
        "-n".into(),
        "p".into(),
        "--outdir".into(),
        dir.display().to_string(),
    ];
    args.extend(extra.iter().map(|s| (*s).to_string()));
    let st = Command::new(RUST).args(&args).output().unwrap();
    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    files.sort();
    Outcome {
        code: st.status.code().unwrap_or(-1),
        files,
    }
}

/// The three files a successful single-end run writes.
fn three_files() -> Vec<String> {
    vec![
        "p_peaks.narrowPeak".to_string(),
        "p_peaks.xls".to_string(),
        "p_summits.bed".to_string(),
    ]
}

/// Every one of these is exit 1 with nothing written.
#[test]
fn the_validator_rejections_leave_no_output_file() {
    let cases: &[&[&str]] = &[
        &["--broad", "--call-summits"],
        &["-q", "0"],
        &["-q", "-1"],
        &["-q", "-0.5"],
        &["-p", "-1"],
        &["--broad", "--broad-cutoff", "0"],
        &["--broad", "--broad-cutoff", "-0.1"],
        &["--d-min", "-1"],
        &["-m", "100", "5"],
        &["--keep-dup", "-1"],
        &["--keep-dup", "1.5"],
        &["--max-gap", "3000000000"],
        &["--min-length", "3000000000"],
        &["--shift", "3000000000"],
        &["-s", "3000000000"],
    ];
    for case in cases {
        let got = run(case);
        assert_eq!(got.code, 1, "callpeak {case:?} exited {}", got.code);
        assert_eq!(
            got.files,
            Vec::<String>::new(),
            "callpeak {case:?} left {:?} behind; the validator runs before any writer",
            got.files
        );
    }
}

#[test]
fn a_negative_max_count_is_only_refused_for_frag() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/it-opt-validator-callpeak/maxcount-frag");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let frag = fixture("frag_basic/barcode_fragments");
    let st = Command::new(RUST)
        .args([
            "callpeak",
            "-t",
            &frag.join("treat.frag").display().to_string(),
            "-c",
            &frag.join("ctrl.frag").display().to_string(),
            "-f",
            "FRAG",
            "-g",
            "10000",
            "-n",
            "p",
            "--max-count",
            "-1",
            "--outdir",
            &dir.display().to_string(),
        ])
        .output()
        .unwrap();
    assert_eq!(st.status.code(), Some(1));
    assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 0);
    // ... and `-f BED` ignores the option entirely: `OptValidator.py:93-96` only
    // looks at `maxcount` inside the `FRAG` branch, and the header line that reports
    // it (`:189-190`) is likewise FRAG-only.
    assert_eq!(run(&["--max-count", "-1"]).files, three_files());
}

#[test]
fn a_negative_tag_size_is_not_a_rejection() {
    // `-s -50` is a legal C `int`; `callpeak` never assigns it to `FWTrack.fw`
    // (`callpeak_cmd.py:386-403` builds the track from the parser's own estimate),
    // so the run succeeds exactly as it does upstream.
    assert_eq!(run(&["-s-50"]).code, 0);
    assert_eq!(run(&["-s-50"]).files, three_files());
}

#[test]
fn a_zero_pvalue_is_the_same_as_not_giving_one() {
    // `OptValidator.py:127` is `if options.pvalue:` -- a truthiness test -- so
    // `-p 0` takes the q-value branch: `-log10(0.05)` cutoff, q-value columns, and
    // the `# qvalue cutoff = 5.00e-02` header line rather than a p-value one.
    let xls = |extra: &[&str], tag: &str| -> String {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../target/it-opt-validator-callpeak")
            .join(tag);
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let se = fixture("se_basic/gauss_two_peaks");
        let mut args: Vec<String> = vec![
            "callpeak".into(),
            "-t".into(),
            se.join("treat_1.bed").display().to_string(),
            "-c".into(),
            se.join("ctrl.bed").display().to_string(),
            "-f".into(),
            "BED".into(),
            "-g".into(),
            "10000".into(),
            "--nomodel".into(),
            "-n".into(),
            "p".into(),
            "--outdir".into(),
            dir.display().to_string(),
        ];
        args.extend(extra.iter().map(|s| (*s).to_string()));
        let st = Command::new(RUST).args(&args).output().unwrap();
        assert!(
            st.status.success(),
            "{extra:?}: {}",
            String::from_utf8_lossy(&st.stderr)
        );
        std::fs::read_to_string(dir.join("p_peaks.xls")).unwrap()
    };
    let with_zero = xls(&["-p", "0"], "pzero");
    let with_none = xls(&[], "pnone");
    assert!(with_zero.contains("# qvalue cutoff = 5.00e-02\n"));
    assert!(!with_zero.contains("pvalue cutoff"));
    // The recorded `# Command line:` line necessarily differs -- it lists the
    // arguments -- so everything else must match line for line.
    let strip = |x: &str| -> Vec<String> {
        x.lines()
            .filter(|l| !l.starts_with("# Command line:"))
            .map(str::to_string)
            .collect()
    };
    assert_eq!(strip(&with_zero), strip(&with_none));
}

#[test]
fn a_non_positive_window_window_does_not_reject() {
    // `PeakDetect.py:187-190` guards each window with `if self.sregion:` /
    // `if self.lregion:`, so a window set to 0 is skipped, and the *no-control*
    // branch (`:298-410`) has no assertions at all.
    for case in [
        vec!["--slocal", "0"],
        vec!["--llocal", "0"],
        vec!["--slocal", "0", "--llocal", "0"],
        vec!["--slocal", "20000", "--llocal", "0"],
    ] {
        assert_eq!(run(&case).code, 0, "{case:?} should be accepted");
    }
    // `--slocal 20000 --llocal 10000` inverts the ladder and does trip the assert.
    assert_eq!(run(&["--slocal", "20000", "--llocal", "10000"]).code, 1);
}
