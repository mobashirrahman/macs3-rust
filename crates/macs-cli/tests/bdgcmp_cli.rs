//! Differential harness for `macs3-rs bdgcmp` and `bdgdiff` against the pinned
//! oracle (F137).
//!
//! `bdgcmp` exposes eight scoring methods over an interval-endpoint table, and
//! `bdgdiff` three peak sets from three likelihood-ratio columns. Both are
//! sensitive to float width at every step, so this sweeps method, pseudocount,
//! depth scaling and input shape against upstream and requires byte-identity.
//!
//! Regenerating the oracle side needs `PYTHONPATH` pointing at the pinned
//! source; the comparison itself is pure `diff`.
//!
//! Neither the pinned checkout nor its interpreter is vendored, so both are located
//! through the environment (see `tests/oracle/mod.rs`) and the two oracle-side tests
//! skip, loudly, where there is no reference installed.

use std::path::{Path, PathBuf};
use std::process::Command;

mod oracle;

const RUST: &str = env!("CARGO_BIN_EXE_macs3-rs");

const METHODS: [&str; 8] = [
    "ppois", "qpois", "subtract", "logFE", "FE", "logLR", "slogLR", "max",
];

/// A deterministic LCG, so the corpus needs no extra dependency.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// Write a bedGraph with `nseg` random-length segments per chromosome.
fn write_bedgraph(path: &Path, seed: u64, nchrom: u64, nseg: u64, vmax: f64, zero_frac: f64) {
    let mut rng = Rng(seed);
    let mut s = String::new();
    for ci in 0..nchrom {
        let mut p = 0u64;
        for _ in 0..nseg {
            let len = 40 + rng.below(360);
            p += len;
            let v = if (rng.below(1000) as f64) < zero_frac * 1000.0 {
                0.0
            } else {
                let t = rng.below(1000) as f64 / 1000.0 * vmax;
                (t * 10.0).round() / 10.0
            };
            s.push_str(&format!("chr{}\t{}\t{}\t{:.1}\n", ci + 1, p - len, p, v));
        }
    }
    std::fs::write(path, s).expect("write bedGraph");
}

/// Run the oracle's `macs3` under the pinned interpreter, with the pinned source on
/// `PYTHONPATH`.
///
/// The interpreter has to be the oracle virtualenv's, not whatever `python3` the host
/// has first on PATH: `bdgcmp`'s output is float formatting all the way down, so an
/// unpinned NumPy is a different answer rather than an error.
fn run_oracle(bin: &Path, python: &Path, src: &Path, args: &[&str], outdir: &Path) -> bool {
    Command::new(python)
        .arg(bin)
        .args(args)
        .current_dir(outdir)
        .env("PYTHONPATH", src)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

fn run_rust(args: &[&str]) -> bool {
    Command::new(RUST)
        .args(args)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

fn read(path: &Path) -> Vec<u8> {
    std::fs::read(path).unwrap_or_else(|_| panic!("read {}", path.display()))
}

fn tmpdir(tag: &str, seed: u64) -> PathBuf {
    let d = std::env::temp_dir().join(format!("macs3-rs-bdg-{tag}-{seed}"));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("mkdir");
    d
}

/// `bdgcmp`: every method at two pseudocounts must be byte-identical.
#[test]
fn bdgcmp_matches_oracle_byte_for_byte() {
    let Some(oracle) = oracle::require("bdgcmp_matches_oracle_byte_for_byte") else {
        return;
    };
    let Some((bin, python)) = oracle::entry_point("bdgcmp_matches_oracle_byte_for_byte", &oracle)
    else {
        return;
    };
    for seed in 1..=6u64 {
        let d = tmpdir("cmp", seed);
        write_bedgraph(&d.join("t.bdg"), seed, 2, 30, 20.0, 0.10);
        write_bedgraph(&d.join("c.bdg"), seed + 100, 2, 30, 8.0, 0.10);
        for m in METHODS {
            for pc in ["0.1", "1.0"] {
                let tag = pc.replace('.', "");
                assert!(
                    run_oracle(
                        bin,
                        python,
                        &oracle.src,
                        &[
                            "bdgcmp",
                            "-t",
                            "t.bdg",
                            "-c",
                            "c.bdg",
                            "-m",
                            m,
                            "-p",
                            pc,
                            "--o-prefix",
                            &format!("u{tag}"),
                        ],
                        &d,
                    ),
                    "oracle bdgcmp -m {m} -p {pc} failed (seed {seed})"
                );
                assert!(
                    run_rust(&[
                        "bdgcmp",
                        "-t",
                        d.join("t.bdg").to_str().unwrap(),
                        "-c",
                        d.join("c.bdg").to_str().unwrap(),
                        "-m",
                        m,
                        "-p",
                        pc,
                        "--outdir",
                        d.to_str().unwrap(),
                        "--o-prefix",
                        &format!("m{tag}"),
                    ]),
                    "rust bdgcmp -m {m} -p {pc} failed (seed {seed})"
                );
                let up = read(&d.join(format!("u{tag}_{m}.bdg")));
                let mine = read(&d.join(format!("m{tag}_{m}.bdg")));
                assert_eq!(
                    up,
                    mine,
                    "bdgcmp -m {m} -p {pc} differs (seed {seed})\n--- upstream\n{}\n--- mine\n{}",
                    String::from_utf8_lossy(&up),
                    String::from_utf8_lossy(&mine),
                );
            }
        }
        let _ = std::fs::remove_dir_all(&d);
    }
}

/// `bdgdiff`: cond1/cond2/common BED must be byte-identical at three depth scalings.
#[test]
fn bdgdiff_matches_oracle_byte_for_byte() {
    let Some(oracle) = oracle::require("bdgdiff_matches_oracle_byte_for_byte") else {
        return;
    };
    let Some((bin, python)) = oracle::entry_point("bdgdiff_matches_oracle_byte_for_byte", &oracle)
    else {
        return;
    };
    for seed in 20..=25u64 {
        let d = tmpdir("diff", seed);
        for (name, sd, vmax, zf) in [
            ("t1", seed, 20.0, 0.15),
            ("c1", seed + 50, 5.0, 0.25),
            ("t2", seed + 100, 18.0, 0.15),
            ("c2", seed + 150, 5.0, 0.25),
        ] {
            write_bedgraph(&d.join(format!("{name}.bdg")), sd, 3, 25, vmax, zf);
        }
        for extra in [
            &[][..],
            &["--d1", "2", "--d2", "1"][..],
            &["--d1", "1", "--d2", "4"][..],
        ] {
            let mut base = vec![
                "bdgdiff", "--t1", "t1.bdg", "--c1", "c1.bdg", "--t2", "t2.bdg", "--c2", "c2.bdg",
            ];
            base.extend_from_slice(extra);
            base.extend_from_slice(&["--o-prefix", "X"]);
            assert!(
                run_oracle(bin, python, &oracle.src, &base, &d),
                "oracle bdgdiff failed (seed {seed})"
            );
            // rename so the peak-name prefix matches, since it is user-supplied
            for (src, dst) in [("cond1", "u1"), ("cond2", "u2"), ("common", "u3")] {
                let _ = std::fs::rename(
                    d.join(format!("X_c3.0_{src}.bed")),
                    d.join(format!("{dst}.bed")),
                );
            }
            let dstr = d.to_str().unwrap();
            let t1 = format!("{dstr}/t1.bdg");
            let c1 = format!("{dstr}/c1.bdg");
            let t2 = format!("{dstr}/t2.bdg");
            let c2 = format!("{dstr}/c2.bdg");
            let mut args = vec![
                "bdgdiff", "--t1", &t1, "--c1", &c1, "--t2", &t2, "--c2", &c2,
            ];
            args.extend_from_slice(extra);
            args.extend_from_slice(&["--outdir", dstr, "--o-prefix", "X"]);
            assert!(run_rust(&args), "rust bdgdiff failed (seed {seed})");
            for (src, up) in [("cond1", "u1"), ("cond2", "u2"), ("common", "u3")] {
                let oracle = read(&d.join(format!("{up}.bed")));
                let mine = read(&d.join(format!("X_c3.0_{src}.bed")));
                assert_eq!(
                    oracle, mine,
                    "bdgdiff {src} differs (seed {seed}, extra {extra:?})"
                );
            }
            // clean for the next scaling
            for p in ["X_c3.0_cond1.bed", "X_c3.0_cond2.bed", "X_c3.0_common.bed"] {
                let _ = std::fs::remove_file(d.join(p));
            }
        }
        let _ = std::fs::remove_dir_all(&d);
    }
}

/// An unknown `-m` value is caught by argparse's `choices` check, so it is a
/// usage error (exit 2) rather than a runtime error. Upstream agrees:
/// `macs3 bdgcmp -m bogus` exits 2 with `invalid choice`.
#[test]
fn bdgcmp_rejects_an_unknown_method() {
    let d = tmpdir("bad", 1);
    write_bedgraph(&d.join("t.bdg"), 1, 1, 5, 10.0, 0.0);
    write_bedgraph(&d.join("c.bdg"), 2, 1, 5, 5.0, 0.0);
    let out = Command::new(RUST)
        .args([
            "bdgcmp",
            "-t",
            d.join("t.bdg").to_str().unwrap(),
            "-c",
            d.join("c.bdg").to_str().unwrap(),
            "-m",
            "bogus",
            "--o-prefix",
            "x",
        ])
        .output()
        .expect("run");
    assert_eq!(out.status.code(), Some(2), "usage error must exit 2");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("invalid choice"), "stderr was: {err}");
    assert!(!err.contains("panicked"), "must not panic: {err}");
    let _ = std::fs::remove_dir_all(&d);
}

/// `-o` and `--o-prefix` are mutually exclusive and one is required.
#[test]
fn bdgdiff_requires_an_output_specification() {
    let d = tmpdir("noout", 1);
    for n in ["t1", "c1", "t2", "c2"] {
        write_bedgraph(&d.join(format!("{n}.bdg")), 1, 1, 4, 10.0, 0.0);
    }
    let out = Command::new(RUST)
        .args([
            "bdgdiff",
            "--t1",
            d.join("t1.bdg").to_str().unwrap(),
            "--c1",
            d.join("c1.bdg").to_str().unwrap(),
            "--t2",
            d.join("t2.bdg").to_str().unwrap(),
            "--c2",
            d.join("c2.bdg").to_str().unwrap(),
        ])
        .output()
        .expect("run");
    assert_eq!(out.status.code(), Some(2), "usage error must exit 2");
    let _ = std::fs::remove_dir_all(&d);
}

/// `--max-gap >= --min-len` is *reported and then ignored*.
///
/// `bdgdiff_cmd.py:43` calls `error("MAXGAP should be smaller than MINLEN! ...")`
/// and falls through -- there is no `sys.exit()` -- so the reference prints the
/// complaint and runs anyway, exiting 0 with all three BED files written and
/// byte-identical to a legal run. An earlier port treated it as a hard error.
#[test]
fn bdgdiff_reports_maxgap_ge_minlen_but_still_runs() {
    let d = tmpdir("gap", 1);
    for n in ["t1", "c1", "t2", "c2"] {
        write_bedgraph(&d.join(format!("{n}.bdg")), 1, 1, 4, 10.0, 0.0);
    }
    let out = Command::new(RUST)
        .args([
            "bdgdiff",
            "--t1",
            d.join("t1.bdg").to_str().unwrap(),
            "--c1",
            d.join("c1.bdg").to_str().unwrap(),
            "--t2",
            d.join("t2.bdg").to_str().unwrap(),
            "--c2",
            d.join("c2.bdg").to_str().unwrap(),
            "-g",
            "300",
            "-l",
            "200",
            "--o-prefix",
            "x",
            "--outdir",
            d.to_str().unwrap(),
        ])
        .output()
        .expect("run");
    assert_eq!(out.status.code(), Some(0));
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("MAXGAP"), "stderr was: {err}");
    for f in ["x_c3.0_cond1.bed", "x_c3.0_cond2.bed", "x_c3.0_common.bed"] {
        assert!(d.join(f).exists(), "{f} is written despite the complaint");
    }
    let _ = std::fs::remove_dir_all(&d);
}
