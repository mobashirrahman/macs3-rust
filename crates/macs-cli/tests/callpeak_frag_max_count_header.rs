//! `callpeak -f FRAG --max-count N` records the cap in the xls header.
//!
//! `OptValidator.py:189-190` appends
//!
//! ```python
//! if options.format == "FRAG" and options.maxcount:
//!     options.argtxt += "# Maximum count in fragment file is set as %d\n" % (options.maxcount)
//! ```
//!
//! to the argument block, which puts the line between the model-fold line and the
//! cutoff line. Both guards matter: the line is **FRAG-only** (`-f BED --max-count
//! 10` never prints it) and `if options.maxcount:` is a truthiness test, so
//! `--max-count 0` -- documented as "keep all counts" -- does not print it either.
//!
//! This port omitted the line entirely. The peak rows themselves were already
//! correct, which is why a coordinate comparison never saw it: the divergence is
//! one header line inside `*_peaks.xls`.
//!
//! Every expectation is a verbatim copy of the pinned oracle's output for
//! `tests/fixtures/frag_basic/barcode_fragments`; nothing here needs Python or the
//! reference at run time.

use std::path::PathBuf;
use std::process::Command;

const RUST: &str = env!("CARGO_BIN_EXE_macs3-rs");

const MAXCOUNT_LINE: &str = "# Maximum count in fragment file is set as ";

struct Run {
    code: i32,
    header: Vec<String>,
    rows: Vec<String>,
}

/// Run `callpeak` on the counted-fragment fixture with `extra` appended, and return
/// the header lines around the model-fold line plus the peak rows.
fn run(tag: &str, extra: &[&str]) -> Run {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/it-frag-max-count")
        .join(tag);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let frag = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/frag_basic/barcode_fragments");
    let mut args: Vec<String> = vec![
        "callpeak".into(),
        "-t".into(),
        frag.join("treat.frag").display().to_string(),
        "-c".into(),
        frag.join("ctrl.frag").display().to_string(),
        "-f".into(),
        "FRAG".into(),
        "-g".into(),
        "10000".into(),
        "-n".into(),
        "p".into(),
        "--outdir".into(),
        dir.display().to_string(),
    ];
    args.extend(extra.iter().map(|s| (*s).to_string()));
    let st = Command::new(RUST).args(&args).output().unwrap();
    assert!(
        st.status.success(),
        "callpeak {extra:?} failed: {}",
        String::from_utf8_lossy(&st.stderr)
    );
    let xls = std::fs::read_to_string(dir.join("p_peaks.xls")).unwrap_or_default();
    let mut header: Vec<String> = xls
        .lines()
        .take_while(|l| l.starts_with('#'))
        .map(str::to_string)
        .collect();
    // drop the two lines that record the run itself
    header.retain(|l| !l.starts_with("# Command line:"));
    Run {
        code: st.status.code().unwrap_or(-1),
        header,
        rows: xls
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(str::to_string)
            .collect(),
    }
}

#[test]
fn a_positive_max_count_is_recorded_between_the_model_fold_and_cutoff_lines() {
    let r = run("ten", &["--max-count", "10"]);
    assert_eq!(r.code, 0);
    let pos = r
        .header
        .iter()
        .position(|l| l == "# model fold = [5, 50]")
        .expect("model-fold line");
    assert_eq!(
        r.header[pos + 1],
        format!("{MAXCOUNT_LINE}10"),
        "the cap line follows the model-fold line; header was {:?}",
        r.header
    );
    assert_eq!(r.header[pos + 2], "# qvalue cutoff = 5.00e-02");
    // The pileup is count-weighted, so a cap really does change the numbers: the
    // same peak carries depth 7989 uncapped and 1998 at `--max-count 1`.
    assert_eq!(
        r.rows.last().unwrap(),
        "chrM\t101\t200\t100\t121\t7989\t633.836\t1.95522\t632.79\tp_peak_1"
    );
    let capped = run("one", &["--max-count", "1"]);
    assert_eq!(
        capped.rows.last().unwrap(),
        "chrM\t101\t200\t100\t121\t1998\t159.657\t1.9531\t158.611\tp_peak_1"
    );
}

#[test]
fn the_cap_line_reports_the_requested_value() {
    assert!(run("one", &["--max-count", "1"])
        .header
        .contains(&format!("{MAXCOUNT_LINE}1")));
    assert!(run("ten", &["--max-count", "10"])
        .header
        .contains(&format!("{MAXCOUNT_LINE}10")));
}

#[test]
fn max_count_zero_keeps_all_counts_and_prints_no_line() {
    // `if options.maxcount:` -- 0 is the "keep everything" setting and prints nothing.
    let zero = run("zero", &["--max-count", "0"]);
    let unset = run("unset", &[]);
    assert_eq!(zero.code, 0);
    assert!(!zero.header.iter().any(|l| l.starts_with(MAXCOUNT_LINE)));
    assert_eq!(zero.rows, unset.rows, "`--max-count 0` keeps every count");
}

#[test]
fn the_cap_line_is_frag_only() {
    let se = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/se_basic/gauss_two_peaks");
    let dir =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/it-frag-max-count/not-frag");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let st = Command::new(RUST)
        .args([
            "callpeak",
            "-t",
            &se.join("treat_1.bed").display().to_string(),
            "-c",
            &se.join("ctrl.bed").display().to_string(),
            "-f",
            "BED",
            "-g",
            "10000",
            "--nomodel",
            "-n",
            "p",
            "--max-count",
            "10",
            "--outdir",
            &dir.display().to_string(),
        ])
        .output()
        .unwrap();
    assert_eq!(st.status.code(), Some(0));
    let xls = std::fs::read_to_string(dir.join("p_peaks.xls")).unwrap();
    assert!(
        !xls.contains(MAXCOUNT_LINE),
        "`-f BED --max-count 10` is ignored, so no cap line"
    );
}
