//! End-to-end test of `macs3-rs filterdup`.
//!
//! Expected bytes were produced by upstream MACS3 3.0.5 `filterdup` on the same
//! input (`FixWidthTrack.print_to_bed`), covering all three `--keep-dup` modes
//! and the documented quirk that a strand holding at most one position bypasses
//! the duplicate limit.

use macs_cli::{parse_flags, Options};

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-fd-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn s(x: &str) -> String {
    x.to_string()
}

fn run(dir: &std::path::Path, input: &str, keep: &str, gsize: Option<&str>) -> String {
    let ofile = format!("out_{keep}.bed");
    let p = dir.join(input);
    let mut args = vec![
        s("-i"),
        s(p.to_str().unwrap()),
        s("-f"),
        s("BED"),
        s("--keep-dup"),
        s(keep),
        s("-o"),
        s(&ofile),
        s("--outdir"),
        s(dir.to_str().unwrap()),
        s("--tsize"),
        s("10"),
    ];
    if let Some(g) = gsize {
        args.push(s("-g"));
        args.push(s(g));
    }
    let o: Options = parse_flags("filterdup", &args).expect("filterdup args parse");
    macs_cli::commands::filterdup::filterdup(&o).expect("filterdup runs");
    std::fs::read_to_string(dir.join(ofile)).expect("output written")
}

/// A mixed-strand fixture: three copies of one `+` tag, a single `-` tag, and
/// two copies of another `+` tag.
const FIXTURE: &str = "chr1\t100\t110\tr\t0\t+\n\
chr1\t100\t110\tr2\t0\t+\n\
chr1\t100\t110\tr3\t0\t+\n\
chr1\t200\t210\tr4\t0\t-\n\
chr1\t300\t310\tr5\t0\t+\n\
chr1\t300\t310\tr6\t0\t+\n";

#[test]
fn keep_dup_two_keeps_two_per_position_and_matches_upstream() {
    let dir = tmpdir("kd2");
    std::fs::write(dir.join("t.bed"), FIXTURE).unwrap();
    let got = run(&dir, "t.bed", "2", None);
    let want = "chr1\t100\t110\t.\t.\t+\n\
chr1\t100\t110\t.\t.\t+\n\
chr1\t300\t310\t.\t.\t+\n\
chr1\t300\t310\t.\t.\t+\n\
chr1\t200\t210\t.\t.\t-\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn keep_dup_auto_matches_upstream_including_the_bypass_quirk() {
    let dir = tmpdir("auto");
    std::fs::write(dir.join("t.bed"), FIXTURE).unwrap();
    // with 6 tags and gsize 2e6 the binomial inverse yields max_dup = 0, yet the
    // lone `-` tag survives (a strand with <= 1 position bypasses the filter).
    let got = run(&dir, "t.bed", "auto", Some("2000000"));
    let want = "chr1\t100\t110\t.\t.\t+\nchr1\t200\t210\t.\t.\t-\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn keep_dup_all_leaves_the_input_untouched() {
    let dir = tmpdir("all");
    std::fs::write(dir.join("t.bed"), FIXTURE).unwrap();
    let got = run(&dir, "t.bed", "all", None);
    // every tag kept (3 at 100, 1 at 200, 2 at 300)
    assert_eq!(got.lines().count(), 6, "all tags retained:\n{got}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// A missing `--tsize` is **accepted**: upstream infers it from the input.
///
/// `filterdup_cmd.py:103-105` only overrides `options.tsize` when it is falsy, and
/// `Parser.tsize()` is the truncated mean of the first ten alignments with a positive
/// parsed length. Every fixture read below is 10 bp, so the inferred `fw` is 10.
///
/// This test previously asserted the opposite -- that a missing `--tsize` is rejected
/// with a `tsize` error. That was the F274 bug (an accept/reject divergence: upstream
/// runs and we exited non-zero), already fixed in `filterdup.rs`; the test outlived the
/// fix and was failing. It now pins the behaviour upstream actually has.
#[test]
fn a_missing_tsize_is_inferred_from_the_input_rather_than_rejected() {
    let dir = tmpdir("bad");
    std::fs::write(dir.join("t.bed"), FIXTURE).unwrap();
    let p = dir.join("t.bed");
    let args = vec![
        s("-i"),
        s(p.to_str().unwrap()),
        s("-f"),
        s("BED"),
        s("--keep-dup"),
        s("1"),
        s("-o"),
        s("o.bed"),
        s("--outdir"),
        s(dir.to_str().unwrap()),
    ];
    let o = parse_flags("filterdup", &args).expect("parses");
    macs_cli::commands::filterdup::filterdup(&o).expect("no tsize is accepted upstream too");
    let got = std::fs::read_to_string(dir.join("o.bed")).expect("output written");
    // `--keep-dup 1` keeps at most one tag per position/strand: 3 at 100 -> 1, 1 at 200, 2 at 300 -> 1.
    // Bytes transcribed from upstream `macs3 filterdup -i t.bed -f BED --keep-dup 1 -o o.bed`, which
    // logs `# tag size = 10` (inferred) and writes the plus strand in coordinate order followed by
    // the minus strand, not in input order.
    let want = "chr1\t100\t110\t.\t.\t+\n\
chr1\t300\t310\t.\t.\t+\n\
chr1\t200\t210\t.\t.\t-\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}
