//! `bdgcmp -p 0` on zero-valued inputs must fail exactly like upstream.
//!
//! With a zero pseudocount, `ppois`/`qpois` hit `AssertionError: Lambda must >
//! 0` (`Prob.py:269`) and `logFE` hits `ZeroDivisionError: float division`
//! (`ScoreTrack.py:150-152`); upstream computes the score before opening the
//! output file (`bdgcmp_cmd.py:70-91`), so the failing method exits 1 and
//! leaves NO output file, while files from earlier successful methods stay.
//! Every other method (`subtract`, `FE`, `logLR`, `slogLR`, `max`) succeeds.
//!
//! Expected bytes below were produced by MACS3 3.0.5 on the same fixtures.
//! No Python or oracle is needed at run time.

use macs_cli::{parse_flags, Options};

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-zerodiv-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn s(x: &str) -> String {
    x.to_string()
}

/// Two intervals: `[0,100)` is 0/0 and `[100,200)` is 5/2.
fn write_zero_pair(dir: &std::path::Path) -> (String, String) {
    let t = dir.join("t.bdg");
    let c = dir.join("c.bdg");
    std::fs::write(&t, "chr1\t0\t100\t0\nchr1\t100\t200\t5\n").unwrap();
    std::fs::write(&c, "chr1\t0\t100\t0\nchr1\t100\t200\t2\n").unwrap();
    (
        t.to_string_lossy().into_owned(),
        c.to_string_lossy().into_owned(),
    )
}

fn bdgcmp_args(t: &str, c: &str, out: &str, methods: &[&str], extra: &[String]) -> Vec<String> {
    let mut args = vec![
        s("-t"),
        s(t),
        s("-c"),
        s(c),
        s("-p"),
        s("0"),
        s("--outdir"),
        s(out),
        s("--o-prefix"),
        s("X"),
    ];
    args.push(s("-m"));
    for m in methods {
        args.push(s(m));
    }
    args.extend(extra.iter().cloned());
    args
}

/// `logFE` at `-p 0` on a zero control: exit 1 via `Err`, an error on stderr
/// (checked by the caller convention), and no output file. This exercises the
/// streaming fast path (single distinct method).
#[test]
fn logfe_zero_denominator_is_an_error_with_no_output_file() {
    let dir = tmpdir("logfe");
    let (t, c) = write_zero_pair(&dir);
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = bdgcmp_args(&t, &c, out.to_str().unwrap(), &["logFE"], &[]);
    let o: Options = parse_flags("bdgcmp", &args).expect("parse");
    let err = macs_cli::commands::bdgcmp::bdgcmp(&o).expect_err("logFE 0/0 must fail");
    assert!(
        err.to_string().contains("float division"),
        "unexpected error: {err}"
    );
    assert!(
        !out.join("X_logFE.bdg").exists(),
        "a failing method must leave no output file"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// `ppois` at `-p 0` on a zero control fails the same way (also streaming).
#[test]
fn ppois_zero_lambda_is_an_error_with_no_output_file() {
    let dir = tmpdir("ppois");
    let (t, c) = write_zero_pair(&dir);
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = bdgcmp_args(&t, &c, out.to_str().unwrap(), &["ppois"], &[]);
    let o: Options = parse_flags("bdgcmp", &args).expect("parse");
    let err = macs_cli::commands::bdgcmp::bdgcmp(&o).expect_err("ppois 0/0 must fail");
    assert!(
        err.to_string().contains("Lambda must > 0"),
        "unexpected error: {err}"
    );
    assert!(!out.join("X_ppois.bdg").exists());
    let _ = std::fs::remove_dir_all(&dir);
}

/// `qpois` at `-p 0` on a zero control fails while computing p (materialising
/// path), again with no output file.
#[test]
fn qpois_zero_lambda_is_an_error_with_no_output_file() {
    let dir = tmpdir("qpois");
    let (t, c) = write_zero_pair(&dir);
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = bdgcmp_args(&t, &c, out.to_str().unwrap(), &["qpois"], &[]);
    let o: Options = parse_flags("bdgcmp", &args).expect("parse");
    let err = macs_cli::commands::bdgcmp::bdgcmp(&o).expect_err("qpois 0/0 must fail");
    assert!(
        err.to_string().contains("Lambda must > 0"),
        "unexpected error: {err}"
    );
    assert!(!out.join("X_qpois.bdg").exists());
    let _ = std::fs::remove_dir_all(&dir);
}

/// The methods that do not divide by the control succeed on the same zeros,
/// byte-identical to upstream.
#[test]
fn other_methods_succeed_on_zeros() {
    let dir = tmpdir("others");
    let (t, c) = write_zero_pair(&dir);
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let want: &[(&str, &str)] = &[
        (
            "subtract",
            "chr1\t0\t100\t0.00000\nchr1\t100\t200\t3.00000\n",
        ),
        ("FE", "chr1\t0\t200\tnan\n"),
        ("logLR", "chr1\t0\t100\t0.00000\nchr1\t100\t200\t0.68682\n"),
        ("slogLR", "chr1\t0\t100\t0.00000\nchr1\t100\t200\t0.68682\n"),
        ("max", "chr1\t0\t100\t0.00000\nchr1\t100\t200\t5.00000\n"),
    ];
    for (m, bytes) in want {
        let args = bdgcmp_args(&t, &c, out.to_str().unwrap(), &[m], &[]);
        let o: Options = parse_flags("bdgcmp", &args).expect("parse");
        macs_cli::commands::bdgcmp::bdgcmp(&o).expect("method {m} must succeed");
        let got = std::fs::read_to_string(out.join(format!("X_{m}.bdg"))).unwrap();
        assert_eq!(got, *bytes, "bdgcmp -m {m} differs from upstream");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// With several methods in one invocation, a later failure keeps the earlier
/// files, creates nothing for the failing method, and never runs the methods
/// after it -- matching `bdgcmp_cmd.py:56-91`.
#[test]
fn later_method_failure_keeps_earlier_files_only() {
    let dir = tmpdir("multi");
    let (t, c) = write_zero_pair(&dir);
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = bdgcmp_args(
        &t,
        &c,
        out.to_str().unwrap(),
        &["subtract", "logFE", "FE"],
        &[],
    );
    let o: Options = parse_flags("bdgcmp", &args).expect("parse");
    macs_cli::commands::bdgcmp::bdgcmp(&o).expect_err("logFE must fail the run");
    assert_eq!(
        std::fs::read_to_string(out.join("X_subtract.bdg")).unwrap(),
        "chr1\t0\t100\t0.00000\nchr1\t100\t200\t3.00000\n"
    );
    assert!(
        !out.join("X_logFE.bdg").exists(),
        "the failing method must leave no file"
    );
    assert!(
        !out.join("X_FE.bdg").exists(),
        "methods after the failure must not run"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// A failing method given twice still goes through the streaming fast path
/// (one distinct method) and leaves no file.
#[test]
fn duplicate_failing_method_leaves_no_file() {
    let dir = tmpdir("dup");
    let (t, c) = write_zero_pair(&dir);
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = bdgcmp_args(&t, &c, out.to_str().unwrap(), &["logFE", "logFE"], &[]);
    let o: Options = parse_flags("bdgcmp", &args).expect("parse");
    macs_cli::commands::bdgcmp::bdgcmp(&o).expect_err("logFE 0/0 must fail");
    assert!(!out.join("X_logFE.bdg").exists());
    let _ = std::fs::remove_dir_all(&dir);
}
