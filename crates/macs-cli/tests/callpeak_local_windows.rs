//! PeakDetect's control-window assertions run after the fragment length is known,
//! including under --nolambda and with fractional paired-end mean lengths.
//! Exit statuses were captured from the pinned MACS3 3.0.5 oracle.

use std::path::Path;
use std::process::Command;

fn run(dir: &Path, format: &str, control: bool, extra: &[&str]) -> std::process::Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_macs3-rs"));
    cmd.args(["callpeak", "-t"]).arg(dir.join("treat"));
    if control {
        cmd.arg("-c").arg(dir.join("control"));
    }
    cmd.args([
        "-f",
        format,
        "-g",
        "100000",
        "--nomodel",
        "--keep-dup",
        "all",
        "-B",
        "-n",
        "target",
        "--outdir",
    ])
    .arg(dir.join("output"))
    .args(extra)
    .output()
    .unwrap()
}

fn single_end_inputs(dir: &Path) {
    std::fs::write(
        dir.join("treat"),
        "chr1\t100\t125\tr1\t0\t+\nchr1\t200\t225\tr2\t0\t+\nchr1\t300\t325\tr3\t0\t-\n",
    )
    .unwrap();
    std::fs::write(
        dir.join("control"),
        "chr1\t20\t45\tc1\t0\t+\nchr1\t250\t275\tc2\t0\t+\nchr1\t600\t625\tc3\t0\t-\n",
    )
    .unwrap();
}

fn assert_rejected(dir: &Path, output: std::process::Output) {
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("can't be smaller"), "{stderr}");
    assert!(
        !dir.join("output").exists()
            || std::fs::read_dir(dir.join("output"))
                .unwrap()
                .next()
                .is_none(),
        "invalid windows must fail before writing peak or signal files"
    );
}

#[test]
fn control_windows_reject_inconsistent_lengths_before_output() {
    let cases: &[&[&str]] = &[
        &["--extsize", "1500"],
        &["--extsize", "1500", "--slocal", "0", "--llocal", "1000"],
        &["--extsize", "200", "--slocal", "2000", "--llocal", "1000"],
        &["--extsize", "200", "--slocal", "-1"],
        &["--extsize", "200", "--llocal", "-1"],
        &["--extsize", "1500", "--nolambda"],
    ];
    for extra in cases {
        let dir = tempfile::tempdir().unwrap();
        single_end_inputs(dir.path());
        assert_rejected(dir.path(), run(dir.path(), "BED", true, extra));
    }
}

#[test]
fn disabled_and_equal_control_windows_remain_valid() {
    for (slocal, llocal) in [("1500", "1500"), ("0", "10000"), ("10000", "0"), ("0", "0")] {
        let dir = tempfile::tempdir().unwrap();
        single_end_inputs(dir.path());
        let output = run(
            dir.path(),
            "BED",
            true,
            &["--extsize", "1500", "--slocal", slocal, "--llocal", llocal],
        );
        assert!(output.status.success(), "{output:?}");
        assert!(dir.path().join("output/target_peaks.xls").exists());
    }
}

#[test]
fn no_control_runs_do_not_apply_control_window_assertions() {
    for extra in [
        &["--extsize", "1500"][..],
        &["--extsize", "200", "--slocal", "2000", "--llocal", "1000"][..],
    ] {
        let dir = tempfile::tempdir().unwrap();
        single_end_inputs(dir.path());
        let output = run(dir.path(), "BED", false, extra);
        assert!(output.status.success(), "{output:?}");
    }
}

#[test]
fn paired_end_windows_use_the_fractional_mean_fragment_length() {
    for format in ["BEDPE", "FRAG"] {
        for slocal in ["1500", "1501"] {
            let dir = tempfile::tempdir().unwrap();
            let suffix = if format == "FRAG" { "\tbarcode\t1" } else { "" };
            std::fs::write(
                dir.path().join("treat"),
                format!("chr1\t100\t1600{suffix}\nchr1\t200\t1701{suffix}\n"),
            )
            .unwrap();
            std::fs::write(
                dir.path().join("control"),
                format!("chr1\t20\t1520{suffix}\nchr1\t250\t1750{suffix}\n"),
            )
            .unwrap();
            let output = run(dir.path(), format, true, &["--slocal", slocal]);
            if slocal == "1500" {
                assert_rejected(dir.path(), output);
            } else {
                assert!(output.status.success(), "{format}: {output:?}");
            }
        }
    }
}
