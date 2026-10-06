//! `predictd -m/--mfold` is `nargs=2`, and both halves have to reach the model.
//!
//! The upper bound used to be read back as the lower one, so `min_tags` and
//! `max_tags` came out equal. `naive_call_peaks` keeps a summit only when its
//! value is both `> min_v` and `< max_v`, so an equal pair satisfies nothing and
//! every peak was discarded. `macs3-rs predictd` then reported "can only find 0
//! paired peaks" and wrote no file on real CTCF data where upstream finds 3,803
//! paired peaks and an identical `predictd_model.R`.
//!
//! `oracle/check_predictd_randsample.sh` did not catch it because it drives
//! `-m` values through the same reader and therefore agreed with itself. These
//! tests pin the observable behaviour instead, so a reader that drops a value
//! fails here rather than only on someone's ChIP-seq file.
//!
//! The fixture is built so the two bounds are *distinguishable*: with the upper
//! bound intact a summit of ~2000 sits inside `(min, max)`; with the upper bound
//! collapsed onto the lower one the same summit sits above `max` and nothing
//! survives.

use macs_cli::parse_flags;
use std::path::{Path, PathBuf};
use std::process::Command;

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-predictd-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn bin() -> PathBuf {
    let mut p = std::env::current_exe().expect("test exe path");
    p.pop();
    if p.ends_with("deps") {
        p.pop();
    }
    p.join("macs3-rs")
}

/// 150 sites, 3 kb apart, each with 2000 `+` and 2000 `-` reads overlapping in a
/// ~600 bp window. `-g 500000` makes `min_tags = 360 * lmfold` and
/// `max_tags = 360 * umfold`, so `-m 3 50` gives `(1080, 18000)` and a summit of
/// ~2000 is inside it -- while `-m 3 3` would cap at 1080 and find nothing, which
/// is precisely the state the collapsed upper bound used to put every run in.
fn write_model_fixture(dir: &Path) -> PathBuf {
    let bed = dir.join("model.bed");
    let mut s = String::new();
    for site in 0..150i64 {
        let x = 1000 + site * 3000;
        for i in 0..2000i64 {
            let jitter = (i * 7919) % 151; // deterministic, no RNG dependency
            s.push_str(&format!(
                "chr1\t{}\t{}\t0\t0\t+\n",
                x + jitter,
                x + 251 + jitter
            ));
        }
        for i in 0..2000i64 {
            let jitter = (i * 6151) % 151;
            s.push_str(&format!(
                "chr1\t{}\t{}\t0\t0\t-\n",
                x + 200 + jitter,
                x + 451 + jitter
            ));
        }
    }
    std::fs::write(&bed, s).unwrap();
    bed
}

fn run_predictd(bed: &Path, out: &Path, mfold: &str) -> (bool, String) {
    std::fs::create_dir_all(out).unwrap();
    // `-m` is nargs=2, so the two bounds must be two separate argv entries, exactly
    // as a shell would pass them.
    let mut argv: Vec<String> = vec![
        "predictd".into(),
        "-i".into(),
        bed.display().to_string(),
        "-f".into(),
        "BED".into(),
        "-g".into(),
        "500000".into(),
        "-m".into(),
    ];
    argv.extend(mfold.split_whitespace().map(String::from));
    argv.push("--outdir".into());
    argv.push(out.display().to_string());

    let outp = Command::new(bin())
        .args(&argv)
        .output()
        .expect("run macs3-rs predictd");
    let produced = out.join("predictd_model.R").is_file();
    (produced, String::from_utf8_lossy(&outp.stderr).into_owned())
}

#[test]
fn a_wider_mfold_builds_a_model_on_data_where_a_collapsed_one_cannot() {
    let d = tmpdir("mfold-range");
    let bed = write_model_fixture(&d);

    let (wide_ok, _) = run_predictd(&bed, &d.join("wide"), "3 50");
    assert!(
        wide_ok,
        "`-m 3 50` must build a model: the summit value sits inside (min_tags, \
         max_tags) only when the upper bound is actually used"
    );

    // The same run with the bounds collapsed is what a reader that keeps only the
    // first value produces. It must NOT build a model -- if this ever starts
    // succeeding, the fixture no longer distinguishes the two readings and this
    // test has stopped testing anything.
    let (narrow_ok, _) = run_predictd(&bed, &d.join("narrow"), "3 3");
    assert!(
        !narrow_ok,
        "`-m 3 3` makes min_tags == max_tags, so no summit can satisfy both \
         bounds. If this now builds a model the fixture needs retuning."
    );
}

#[test]
fn mfold_parses_as_two_values_not_one() {
    // The mechanical half of the same bug: `Options::get` returns only the first
    // value of an nargs=2 option, which is how the upper bound went missing.
    let o = parse_flags(
        "predictd",
        &[
            "-i".to_string(),
            "x.bed".to_string(),
            "-f".to_string(),
            "BED".to_string(),
            "-m".to_string(),
            "100".to_string(),
            "200".to_string(),
        ],
    )
    .expect("valid invocation");
    assert_eq!(
        o.get_all("mfold"),
        ["100", "200"],
        "`-m 100 200` must parse as two values"
    );
    assert_eq!(
        o.get("mfold"),
        Some("100"),
        "documented behaviour of get(): the first value only, which is why the \
         command must read get_all for an nargs=2 flag"
    );
}

#[test]
fn a_bad_mfold_value_is_rejected() {
    let d = tmpdir("mfold-bad");
    let bed = write_model_fixture(&d);
    // `type=int` with `nargs=2`: upstream refuses a non-integer bound, and a
    // single value or three values cannot satisfy nargs=2. Each case below must be
    // a non-zero exit, whatever layer rejects it. `-m 0 50` is deliberately absent:
    // upstream accepts it (exit 0), so asserting otherwise would pin a divergence
    // rather than a behaviour.
    for bad in ["3 abc", "3", "3 50 7", "3 50.5"] {
        let mut argv: Vec<String> = vec![
            "predictd".into(),
            "-i".into(),
            bed.display().to_string(),
            "-f".into(),
            "BED".into(),
            "-g".into(),
            "500000".into(),
            "-m".into(),
        ];
        argv.extend(bad.split_whitespace().map(String::from));
        argv.push("--outdir".into());
        argv.push(d.join("bad").display().to_string());

        let outp = Command::new(bin())
            .args(&argv)
            .output()
            .expect("run macs3-rs predictd");
        assert_ne!(
            outp.status.code(),
            Some(0),
            "`-m {bad}` must not exit 0; stderr was: {}",
            String::from_utf8_lossy(&outp.stderr)
        );
        assert!(
            !d.join("bad").join("predictd_model.R").exists(),
            "`-m {bad}` must not leave a model behind"
        );
    }
}
