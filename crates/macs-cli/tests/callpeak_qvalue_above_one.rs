//! A q-value above 1 makes the `-log10` cutoff **negative**, and every position
//! then clears it.
//!
//! `OptValidator.py:132` computes `log_qvalue = log(qvalue, 10) * -1`, so `-q 2.0`
//! asks for `score > -0.301`. The q-values themselves are floored at zero --
//! `__cal_pvalue_qvalue_table` breaks out of its walk with `q = 0` rather than
//! storing a negative one (`CallPeakUnit.py:990-993`) -- so `0 > -0.301` holds and
//! the above-cutoff set is the **whole paired array**, zero treatment depth
//! included. `__chrom_call_peak_using_certain_criteria` never looks at the
//! treatment pileup (`CallPeakUnit.py:1202`); it is only `apply_multiple_cutoffs`
//! that decides.
//!
//! This port used to add `tpos[i] > 0.0` to that test, which dropped the zero-depth
//! positions at the front of the run. The first surviving chunk then started at
//! `pos[i-1]` instead of at the literal `0` that `above_cutoff[0] == 0` forces
//! (`CallPeakUnit.py:1216-1220`), so the peak start moved: on the fixture below with
//! `-q 2.0`, `chrB 0 1293` upstream against `chrB 1001 1293` here, with the summit
//! offset column 1147 against 146 for the same absolute summit.
//!
//! Every expectation below is a verbatim copy of the pinned oracle's output for this
//! fixture. Nothing here needs Python or the reference at run time.

use std::path::PathBuf;
use std::process::Command;

const RUST: &str = env!("CARGO_BIN_EXE_macs3-rs");

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/sweep/gtiny_mse_d400_w60_ctrl_001")
}

struct Run {
    summits: String,
    narrow: String,
    broad: String,
    gapped: String,
}

/// Run `callpeak` on the fixture and return the peak rows of each output, with the
/// `#`-prefixed xls header lines dropped (they record the command line, which
/// contains this test's own temporary directory).
fn callpeak(extra: &[&str]) -> Run {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target")
        .join("it-qvalue-above-one")
        .join(extra.join("_").replace(['-', '.'], "_"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let f = fixture();
    let mut args: Vec<String> = vec![
        "callpeak".into(),
        "-t".into(),
        f.join("treat_1.bed").display().to_string(),
        "-c".into(),
        f.join("ctrl.bed").display().to_string(),
        "-f".into(),
        "BED".into(),
        "-g".into(),
        "10000".into(),
        "--nomodel".into(),
        "-q".into(),
        "2.0".into(),
        "-n".into(),
        "p".into(),
        "--outdir".into(),
        dir.display().to_string(),
    ];
    args.extend(extra.iter().map(|s| (*s).to_string()));
    let st = Command::new(RUST).args(&args).output().unwrap();
    assert!(
        st.status.success(),
        "callpeak {:?} failed: {}",
        extra,
        String::from_utf8_lossy(&st.stderr)
    );
    let read = |name: &str| match std::fs::read_to_string(dir.join(name)) {
        Ok(text) => text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| format!("{l}\n"))
            .collect::<String>(),
        // a narrow run writes no broadPeak, and a broad run no narrowPeak
        Err(_) => String::new(),
    };
    Run {
        summits: read("p_summits.bed"),
        narrow: read("p_peaks.narrowPeak"),
        broad: read("p_peaks.broadPeak"),
        gapped: read("p_peaks.gappedPeak"),
    }
}

const EXPECTED_NARROW: &str = "chrB\t0\t1293\tp_peak_1\t587\t.\t6.81608\t60.007\t58.7083\t1147\n";
const EXPECTED_SUMMITS: &str = "chrB\t1147\t1148\tp_peak_1\t58.7083\n";
const EXPECTED_BROAD: &str = "chrB\t1033\t1268\tp_peak_1\t447\t.\t5.74758\t45.914\t44.7099\n";
const EXPECTED_GAPPED: &str =
    "chrB\t1033\t1268\tp_peak_1\t447\t.\t0\t0\t0\t2\t1,1\t0,234\t5.74758\t45.914\t44.7099\n";

#[test]
fn a_qvalue_above_one_still_starts_the_peak_at_zero() {
    let r = callpeak(&[]);
    assert_eq!(r.narrow, EXPECTED_NARROW);
    assert_eq!(r.summits, EXPECTED_SUMMITS);
    // The reported start is the xls column, one base above the BED start. The point
    // of the regression is the BED `0`: a zero-depth position opens the run and
    // `above_cutoff[0] == 0` pins its start at coordinate zero.
    assert!(r.narrow.starts_with("chrB\t0\t"), "got {:?}", r.narrow);
}

#[test]
fn a_qvalue_above_one_with_call_summits() {
    let r = callpeak(&["--call-summits"]);
    assert_eq!(
        r.narrow,
        "chrB\t0\t1293\tp_peak_1\t594\t.\t6.94993\t60.9435\t59.4001\t1148\n"
    );
    assert_eq!(r.summits, "chrB\t1148\t1149\tp_peak_1\t59.4001\n");
}

#[test]
fn a_qvalue_above_one_under_broad() {
    let r = callpeak(&["--broad"]);
    assert_eq!(r.broad, EXPECTED_BROAD);
    assert_eq!(r.gapped, EXPECTED_GAPPED);
}

#[test]
fn the_narrow_answer_is_not_the_default_qvalue_answer() {
    // The control for the whole file: at the default `-q 0.05` the same fixture
    // yields `chrB 1033 1267` with summit offset 114. Raising the cutoff above 1
    // therefore *extends* the peak to the whole contig rather than doing nothing --
    // the negative cutoff clears the zero-depth positions too.
    let dir =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../target/it-qvalue-above-one/control");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let f = fixture();
    let st = Command::new(RUST)
        .args([
            "callpeak",
            "-t",
            &f.join("treat_1.bed").display().to_string(),
            "-c",
            &f.join("ctrl.bed").display().to_string(),
            "-f",
            "BED",
            "-g",
            "10000",
            "--nomodel",
            "-q",
            "0.05",
            "-n",
            "p",
            "--outdir",
            &dir.display().to_string(),
        ])
        .output()
        .unwrap();
    assert!(st.status.success());
    let narrow = std::fs::read_to_string(dir.join("p_peaks.narrowPeak")).unwrap_or_default();
    assert_eq!(
        narrow,
        "chrB\t1033\t1267\tp_peak_1\t587\t.\t6.81608\t60.007\t58.7083\t114\n"
    );
    assert_ne!(narrow, EXPECTED_NARROW);
}
