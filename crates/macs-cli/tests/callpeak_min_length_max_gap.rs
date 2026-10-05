//! `--min-length` and `--max-gap` are truthiness tests, and they *change the
//! answer*.
//!
//! `PeakDetect.__init__` (`PeakDetect.py:76-83`):
//!
//! ```python
//! if opt.maxgap: self.maxgap = opt.maxgap
//! else:          self.maxgap = opt.tsize
//! if opt.minlen: self.minlen = opt.minlen
//! else:          self.minlen = self.d
//! ```
//!
//! Both are **truthiness** tests, not `is not None` tests, so `--max-gap 0` and
//! `--min-length 0` fall back to the tag size and `d` -- which is what the xls
//! header says for those values too (`OptValidator.py:206-213`). A *negative*
//! value is truthy and is used as given, and the xls reports it.
//!
//! This port ignored both flags entirely, always using the tag size and `d`. On
//! the fixture below `--min-length -1` is the visible case: a 138-base peak is
//! kept upstream (138 >= -1) and dropped here (138 < d = 400), and the reported
//! start differs with it.
//!
//! A negative `--max-gap` additionally forces the gap comparison to be signed: it
//! is `tl <= max_gap` against a `cython.int`, and every gap is non-negative, so a
//! negative value means **no chunk ever merges**. A `u32` gap cannot express that
//! (it would clamp to 0 and merge the contiguous ones), which is why the merge
//! gap is carried as `i64` through [`macs_peaks::CallParams`].
//!
//! `lvl2_max_gap` for `--broad` is `self.maxgap * 4` on the **float**
//! (`PeakDetect.py:264`), truncated by the `cython.int` parameter at the call, so
//! it comes from the mean fragment length when the option is absent and from the
//! integer when it is given.
//!
//! Every peak row below is a verbatim copy of the pinned oracle's output for
//! `tests/fixtures/sweep/gmini_mse_d400_w60_ctrl_000`; nothing here needs Python
//! or the reference at run time.

use std::path::PathBuf;
use std::process::Command;

const RUST: &str = env!("CARGO_BIN_EXE_macs3-rs");

struct Run {
    code: i32,
    header: Vec<String>,
    rows: Vec<String>,
}

fn callpeak(tag: &str, extra: &[&str]) -> Run {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/it-min-length-max-gap")
        .join(tag);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let f = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/sweep/gmini_mse_d400_w60_ctrl_000");
    let mut argv: Vec<String> = vec![
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
        "-n".into(),
        "p".into(),
        "--outdir".into(),
        dir.display().to_string(),
    ];
    argv.extend(extra.iter().map(|s| (*s).to_string()));
    let st = Command::new(RUST).args(&argv).output().unwrap();
    let narrow = std::fs::read_to_string(dir.join("p_peaks.narrowPeak")).unwrap_or_default();
    let xls = std::fs::read_to_string(dir.join("p_peaks.xls")).unwrap_or_default();
    Run {
        code: st.status.code().unwrap_or(-1),
        header: xls
            .lines()
            .filter(|l| l.starts_with("# The m"))
            .map(str::to_string)
            .collect(),
        rows: narrow
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(str::to_string)
            .collect(),
    }
}

const GAP_LINE: &str =
    "# The maximum gap between significant sites is assigned as the read length/tag size.";
const LEN_LINE: &str =
    "# The minimum length of peaks is assigned as the predicted fragment length \"d\".";

#[test]
fn omitting_the_flags_reports_the_fallback_lines() {
    let r = callpeak("unset", &[]);
    assert_eq!(r.code, 0);
    assert_eq!(r.header, vec![GAP_LINE.to_string(), LEN_LINE.to_string()]);
}

#[test]
fn zero_is_falsy_and_falls_back_like_omitting_the_flag() {
    let zero = callpeak("zero", &["--max-gap", "0", "--min-length", "0"]);
    let unset = callpeak("unset-zero", &[]);
    assert_eq!(
        zero.header, unset.header,
        "`--max-gap 0` is not a value of 0"
    );
    assert_eq!(zero.rows, unset.rows);
}

#[test]
fn a_positive_value_is_reported_and_used() {
    let r = callpeak("positive", &["--max-gap", "200", "--min-length", "500"]);
    assert_eq!(r.code, 0);
    assert_eq!(
        r.header,
        vec![
            "# The maximum gap between significant sites = 200".to_string(),
            "# The minimum length of peaks = 500".to_string(),
        ]
    );
    // d is 400 on this fixture and the whole above-cutoff region is 138 bases, so
    // `--min-length 500` rejects a peak that the default's `min_length = d` also
    // rejects: the default is empty here too. The point of the case is that the
    // value is *read and reported*, not that it changes the answer on this input.
    let default = callpeak("default", &[]);
    assert_eq!(
        default.header,
        vec![GAP_LINE.to_string(), LEN_LINE.to_string()]
    );
    assert_eq!(r.rows, default.rows);
}

#[test]
fn a_negative_value_is_truthy_and_keeps_a_short_peak() {
    // Upstream reports `# The minimum length of peaks = -1` and keeps
    // `chrM 163 300` (138 bases, which is >= -1 but < d = 400).
    let r = callpeak("negative", &["--min-length", "-1"]);
    assert_eq!(r.code, 0);
    assert_eq!(
        r.header,
        vec![
            GAP_LINE.to_string(),
            "# The minimum length of peaks = -1".to_string(),
        ]
    );
    assert_eq!(
        r.rows,
        vec!["chrM\t162\t300\tp_peak_1\t53\t.\t1.87737\t5.76413\t5.34009\t25".to_string()],
        "a 138-base peak survives `--min-length -1`"
    );
    let default = callpeak("default-neg", &[]);
    assert!(
        default.rows.is_empty(),
        "without the flag the same peak is dropped: {:?}",
        default.rows
    );
}

#[test]
fn a_negative_max_gap_keeps_the_region_unsplit() {
    let r = callpeak("neg-gap", &["--max-gap", "-1"]);
    assert_eq!(r.code, 0);
    assert_eq!(
        r.header,
        vec![
            "# The maximum gap between significant sites = -1".to_string(),
            LEN_LINE.to_string(),
        ]
    );
    // every gap is non-negative and therefore exceeds -1, so no chunk ever merges:
    // each one becomes its own region, and at 138 bases each is shorter than
    // `min_length = d = 400`, so none is reported. Combining it with
    // `--min-length -1` lifts that length floor and the 81 split chunks all appear,
    // which is what shows the merge really is disabled rather than merely clamped.
    let default = callpeak("default-gap", &[]);
    assert_eq!(r.rows, default.rows);
    let with_floor = callpeak("neg-gap-floor", &["--max-gap", "-1", "--min-length", "-1"]);
    assert_eq!(
        with_floor.rows.len(),
        81,
        "a negative --max-gap makes every chunk its own region"
    );
}
