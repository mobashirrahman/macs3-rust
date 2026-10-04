//! Layer-5 regression: chromosome-level parallelism must be output-identical to
//! one thread.
//!
//! The porting plan requires "chromosome-level Rayon parallelism that is provably
//! output-identical to 1 thread", and the acceptance criteria require identical
//! output at any thread count. `callpeak` has **no** `--threads` flag -- upstream's
//! argparse does not define one, so adding a CLI flag would break the
//! "identical accept/reject behaviour" criterion -- so the count is an environment
//! override that this harness drives.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The `macs3-rs` binary under test, next to the integration test's own exe.
fn binary() -> PathBuf {
    let mut p = std::env::current_exe().expect("test exe path");
    p.pop(); // deps/
    p.pop(); // debug|release/
    p.join("macs3-rs")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-threads-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A multi-chromosome fixture with reads on both, so more than one chromosome is
/// actually scored and the work can be scheduled.
fn fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let treat = dir.join("treat.bed");
    let ctrl = dir.join("ctrl.bed");
    let mut t = String::new();
    let mut c = String::new();
    for (chrom, base) in [("chrA", 1_000i64), ("chrB", 40_000), ("chrC", 90_000)] {
        for i in 0..40i64 {
            t.push_str(&format!(
                "{chrom}\t{}\t{}\tr\t0\t+\n",
                base + i,
                base + i + 10
            ));
            t.push_str(&format!(
                "{chrom}\t{}\t{}\tr\t0\t-\n",
                base + 500 + i,
                base + 510 + i
            ));
        }
        for i in 0..40i64 {
            c.push_str(&format!(
                "{chrom}\t{}\t{}\tr\t0\t+\n",
                base + 20_000 + i * 60,
                base + 20_060 + i * 60
            ));
        }
    }
    std::fs::write(&treat, &t).unwrap();
    std::fs::write(&ctrl, &c).unwrap();
    (treat, ctrl)
}

fn run(out: &Path, treat: &Path, ctrl: &Path, threads: usize) -> Vec<(String, String)> {
    let _ = std::fs::remove_dir_all(out);
    std::fs::create_dir_all(out).unwrap();
    let status = Command::new(binary())
        .args([
            "callpeak",
            "-n",
            "w",
            "-g",
            "100000",
            "-t",
            treat.to_str().unwrap(),
            "-c",
            ctrl.to_str().unwrap(),
            "-f",
            "BED",
            "--nomodel",
            "--extsize",
            "200",
            "--outdir",
            out.to_str().unwrap(),
        ])
        .env("MACS3_RS_THREADS", threads.to_string())
        .status()
        .expect("macs3-rs runs");
    assert!(status.success(), "callpeak failed at {threads} threads");
    let mut files: Vec<(String, String)> = ["w_peaks.xls", "w_peaks.narrowPeak", "w_summits.bed"]
        .iter()
        .map(|f| {
            (
                (*f).to_string(),
                std::fs::read_to_string(out.join(f)).expect("output written"),
            )
        })
        .collect();
    // the `# Command line:` line embeds --outdir, which differs by construction
    for (_, body) in &mut files {
        if let Some(i) = body.find("--outdir") {
            let end = body[i..].find('\n').map(|n| i + n).unwrap_or(body.len());
            body.replace_range(i..end, "--outdir <DIR>");
        }
    }
    files
}

#[test]
fn output_is_identical_at_1_2_8_and_32_threads() {
    let dir = tmpdir("det");
    let (treat, ctrl) = fixture(&dir);
    let ref_run = run(&dir.join("ref"), &treat, &ctrl, 1);
    assert!(
        ref_run
            .iter()
            .any(|(n, b)| n == "w_peaks.xls" && b.contains("chrA")),
        "the reference run must actually call peaks"
    );
    for threads in [2usize, 8, 32] {
        let got = run(&dir.join(format!("t{threads}")), &treat, &ctrl, threads);
        assert_eq!(
            ref_run, got,
            "output differs at {threads} threads; parallelism must not be observable"
        );
    }
}
