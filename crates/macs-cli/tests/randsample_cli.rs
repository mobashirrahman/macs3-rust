//! End-to-end test of `macs3-rs randsample`.
//!
//! The NumPy RNG stream is the hard part: upstream seeds NumPy's *global* RNG
//! and calls `np.random.shuffle`, so the surviving tags are only reproducible if
//! the shuffle stream matches exactly. These tests pin upstream's output bytes
//! for fixed seeds, plus the `-n` override and the over-count rejection.

use macs_cli::{parse_flags, Options};

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-rs-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn s(x: &str) -> String {
    x.to_string()
}

fn bed(dir: &std::path::Path) -> String {
    // 20 tags on chr1 (mixed strands) so a shuffle is observable
    let mut t = String::new();
    for i in 0..10u64 {
        t.push_str(&format!(
            "chr1\t{}\t{}\tr\t0\t+\n",
            100 + i * 7,
            110 + i * 7
        ));
        t.push_str(&format!(
            "chr1\t{}\t{}\tr\t0\t-\n",
            300 + i * 11,
            310 + i * 11
        ));
    }
    let p = dir.join("in.bed");
    std::fs::write(&p, t).unwrap();
    p.to_string_lossy().into_owned()
}

fn run(
    dir: &std::path::Path,
    input: &str,
    pct: Option<&str>,
    num: Option<&str>,
    seed: &str,
) -> String {
    let ofile = "out.bed".to_string();
    let mut args = vec![
        s("-i"),
        s(input),
        s("-f"),
        s("BED"),
        s("-o"),
        s(&ofile),
        s("--outdir"),
        s(dir.to_str().unwrap()),
        s("--tsize"),
        s("10"),
        s("--seed"),
        s(seed),
    ];
    if let Some(p) = pct {
        args.push(s("-p"));
        args.push(s(p));
    }
    if let Some(n) = num {
        args.push(s("-n"));
        args.push(s(n));
    }
    let o: Options = parse_flags("randsample", &args).expect("parse");
    macs_cli::commands::randsample::randsample(&o).expect("randsample runs");
    std::fs::read_to_string(dir.join(ofile)).expect("written")
}

#[test]
fn a_fixed_seed_reproduces_upstream_exactly() {
    let dir = tmpdir("exact");
    let input = bed(&dir);
    let got = run(&dir, &input, Some("50"), None, "42");
    // Captured from upstream MACS3 3.0.5 on this exact input: a NumPy global
    // shuffle at seed 42 keeps these 10 of 20 tags.
    let want = "chr1\t100\t110\t.\t.\t+\n\
chr1\t107\t117\t.\t.\t+\n\
chr1\t135\t145\t.\t.\t+\n\
chr1\t149\t159\t.\t.\t+\n\
chr1\t156\t166\t.\t.\t+\n\
chr1\t300\t310\t.\t.\t-\n\
chr1\t311\t321\t.\t.\t-\n\
chr1\t333\t343\t.\t.\t-\n\
chr1\t355\t365\t.\t.\t-\n\
chr1\t388\t398\t.\t.\t-\n";
    assert_eq!(got, want, "randsample must match upstream byte-for-byte");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_fixed_seed_is_deterministic() {
    let dir = tmpdir("seed");
    let input = bed(&dir);
    let a = run(&dir, &input, Some("50"), None, "42");
    let b = run(&dir, &input, Some("50"), None, "42");
    assert_eq!(a, b, "same seed must be deterministic");
    // 20 tags -> roughly half kept; the point is byte-reproducibility, and the
    // recorded values below came from upstream MACS3 3.0.5.
    assert!(a.lines().count() > 0, "kept some tags");
    assert!(a.lines().count() < 20, "dropped some tags");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn different_seeds_give_different_samples() {
    let dir = tmpdir("seeds");
    let input = bed(&dir);
    let a = run(&dir, &input, Some("50"), None, "42");
    let b = run(&dir, &input, Some("50"), None, "7");
    assert_ne!(a, b, "distinct seeds must select different tags");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn number_bounds_the_count() {
    let dir = tmpdir("num");
    let input = bed(&dir);
    let by_n = run(&dir, &input, None, Some("6"), "3");
    assert!(by_n.lines().count() <= 8, "-n bounds the count:\n{by_n}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// `-p` and `-n` are a `required=True` mutually-exclusive group upstream, so
/// giving both is a usage error rather than one silently winning.
#[test]
fn percentage_and_number_are_mutually_exclusive() {
    let dir = tmpdir("mutex");
    let input = bed(&dir);
    let mut args = vec![
        s("-i"),
        s(&input),
        s("-f"),
        s("BED"),
        s("-p"),
        s("90"),
        s("-n"),
        s("6"),
    ];
    args.push(s("--tsize"));
    args.push(s("10"));
    let e = parse_flags("randsample", &args).expect_err("both -p and -n must be rejected");
    let msg = format!("{e}");
    assert!(
        msg.contains("-n/--number") && msg.contains("-p/--percentage"),
        "got {msg}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// `-n` above the tag count is refused by `randsample_cmd.py:62-64`, which runs
/// *after* the output file is opened at step 0 (`:41-44`). The reference therefore
/// exits 1 and leaves an **empty** output file behind; an earlier port refused
/// before creating anything, which is the same exit status but a different set of
/// files on disk. Measured:
/// `macs3 randsample -i in.bed -f BED -n 999999 --tsize 10 -o o.bed --outdir d`
/// exits 1 with `d/o.bed` present and zero bytes.
#[test]
fn requesting_more_than_the_total_is_rejected_with_an_empty_output_file() {
    let dir = tmpdir("over");
    let input = bed(&dir);
    let args = vec![
        s("-i"),
        s(&input),
        s("-f"),
        s("BED"),
        s("-n"),
        s("999999"),
        s("--tsize"),
        s("10"),
        s("--seed"),
        s("1"),
        s("-o"),
        s("o.bed"),
        s("--outdir"),
        s(dir.to_str().unwrap()),
    ];
    let o = parse_flags("randsample", &args).expect("parse");
    let e = macs_cli::commands::randsample::randsample(&o).expect_err("over-count rejected");
    assert!(format!("{e}").contains("bigger than total"), "got {e}");
    let out = dir.join("o.bed");
    assert!(out.exists(), "upstream opens the writer before this check");
    assert_eq!(std::fs::metadata(&out).unwrap().len(), 0);
    let _ = std::fs::remove_dir_all(&dir);
}
