//! No-control local windows, including upstream's inverted negative window.

use std::process::Command;

#[test]
fn negative_local_window_matches_the_oracle_without_panicking() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.bed");
    std::fs::write(
        &reads,
        "chr1\t100\t125\tr1\t0\t+\nchr1\t120\t145\tr2\t0\t+\nchr1\t150\t175\tr3\t0\t-\n",
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_macs3-rs"))
        .args(["callpeak", "-t"])
        .arg(&reads)
        .args([
            "-f",
            "BED",
            "-g",
            "1000",
            "--nomodel",
            "--extsize",
            "50",
            "--llocal",
            "-3",
            "-B",
            "--keep-dup",
            "all",
            "-n",
            "negative",
            "--outdir",
        ])
        .arg(dir.path())
        .output()
        .unwrap();
    assert!(output.status.success(), "{:?}", output);
    // Pinned MACS3: floor division gives a three-base inverted window;
    // negative depth times the negative scaling factor becomes positive.
    let control = std::fs::read_to_string(dir.path().join("negative_control_lambda.bdg")).unwrap();
    assert_eq!(
        control,
        concat!(
            "chr1\t0\t99\t0.15000\n",
            "chr1\t99\t102\t16.66667\n",
            "chr1\t102\t119\t0.15000\n",
            "chr1\t119\t122\t16.66667\n",
            "chr1\t122\t173\t0.15000\n",
        )
    );
    let peaks = std::fs::read_to_string(dir.path().join("negative_peaks.narrowPeak")).unwrap();
    assert_eq!(
        peaks,
        "chr1\t102\t173\tnegative_peak_1\t24\t.\t3.47826\t4.72783\t2.48978\t35\n"
    );

    // A zero window still fails with the oracle's runtime-error exit status.
    let output = Command::new(env!("CARGO_BIN_EXE_macs3-rs"))
        .args(["callpeak", "-t"])
        .arg(&reads)
        .args(["-g", "1000", "--nomodel", "--llocal", "0"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
}
