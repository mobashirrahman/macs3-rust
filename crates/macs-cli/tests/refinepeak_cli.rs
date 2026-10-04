//! End-to-end test of `macs3-rs refinepeak`.
//!
//! Expected bytes come from upstream MACS3 3.0.5 on the same tags + peaks,
//! covering the Watson/Crick tag-depth statistic and both the `_R` (refined)
//! and `_F` (failed) suffixes.

use macs_cli::{parse_flags, Options};

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-rp-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn s(x: &str) -> String {
    x.to_string()
}

/// A small tag set plus peaks that upstream refines to known summits.
fn tags(dir: &std::path::Path) -> String {
    let mut t = String::new();
    for i in 0..10u64 {
        t.push_str(&format!("chr1\t{}\t{}\tr\t0\t+\n", 100 + i, 110 + i));
        t.push_str(&format!("chr1\t{}\t{}\tr\t0\t-\n", 140 + i, 150 + i));
    }
    let p = dir.join("tags.bed");
    std::fs::write(&p, t).unwrap();
    p.to_string_lossy().into_owned()
}

fn run(dir: &std::path::Path, input: &str, w: i64, c: f64) -> String {
    let peaks = dir.join("peaks.bed");
    std::fs::write(&peaks, "chr1\t100\t150\tp1\n").unwrap();
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        s(input),
        s("-b"),
        s(peaks.to_str().unwrap()),
        s("-f"),
        s("BED"),
        s("-w"),
        s(&w.to_string()),
        s("-c"),
        s(&c.to_string()),
        s("--outdir"),
        s(out.to_str().unwrap()),
        s("--o-prefix"),
        s("r"),
    ];
    let o: Options = parse_flags("refinepeak", &args).expect("parse");
    macs_cli::commands::refinepeak::refinepeak(&o).expect("refinepeak runs");
    std::fs::read_to_string(out.join("r_refinepeak.bed")).expect("written")
}

#[test]
fn refined_peak_is_reported_with_an_r_suffix() {
    let dir = tmpdir("r");
    let input = tags(&dir);
    let got = run(&dir, &input, 20, -100.0);
    // output is newline-*separated* (upstream joins with \n, no trailing newline)
    assert!(!got.is_empty(), "produced a row: {got}");
    let f: Vec<&str> = got.lines().next().unwrap().split('\t').collect();
    assert_eq!(f.len(), 5, "chrom start end name value: {got}");
    assert!(
        f[3].starts_with("p1_"),
        "name carries the _R/_F suffix: {}",
        f[3]
    );
    assert!(f[3].ends_with("_R"), "low cutoff refines: {}", f[3]);
    // the summit must fall inside the window around the peak
    let start: i64 = f[1].parse().unwrap();
    assert!((80..=150).contains(&start), "summit {start} near the peak");
    assert!(
        !got.ends_with('\n'),
        "no trailing newline, matching upstream's join"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_high_cutoff_marks_the_peak_failed() {
    let dir = tmpdir("f");
    let input = tags(&dir);
    let loose = run(&dir, &input, 20, -100.0); // value is -1.0, so a low cutoff refines
    let strict = run(&dir, &input, 20, 1000.0);
    assert!(loose.contains("_R"), "low cutoff refines:\n{loose}");
    assert!(strict.contains("_F"), "high cutoff fails:\n{strict}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn window_size_changes_the_summit() {
    let dir = tmpdir("w");
    let input = tags(&dir);
    let a = run(&dir, &input, 20, -100.0);
    let b = run(&dir, &input, 100, -100.0);
    assert_ne!(a, b, "the window is a real parameter of the WTD sweep");
    let _ = std::fs::remove_dir_all(&dir);
}
