//! End-to-end test of `macs3-rs pileup`.
//!
//! Expected bytes come from upstream MACS3 3.0.5 `PileupV2.pileup_and_write_se`
//! on the same input, covering both the directional and `--both-direction`
//! extension modes (which select different five/three shifts).

use macs_cli::{parse_flags, Options};

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-pl-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn s(x: &str) -> String {
    x.to_string()
}

fn bed(dir: &std::path::Path) -> String {
    let mut t = String::new();
    t.push_str("chr1\t100\t110\tr\t0\t+\n");
    t.push_str("chr1\t105\t115\tr2\t0\t+\n");
    t.push_str("chr1\t200\t210\tr3\t0\t-\n");
    t.push_str("chr1\t300\t310\tr4\t0\t+\n");
    let p = dir.join("in.bed");
    std::fs::write(&p, t).unwrap();
    p.to_string_lossy().into_owned()
}

fn run(dir: &std::path::Path, input: &str, ext: i64, both: bool) -> String {
    let mut args = vec![
        s("-i"),
        s(input),
        s("-f"),
        s("BED"),
        s("-o"),
        s("o.bedGraph"),
        s("--outdir"),
        s(dir.to_str().unwrap()),
        s("--extsize"),
        s(&ext.to_string()),
    ];
    if both {
        args.push(s("-B"));
    }
    let o: Options = parse_flags("pileup", &args).expect("parse");
    macs_cli::commands::pileup::pileup(&o).expect("pileup runs");
    std::fs::read_to_string(dir.join("o.bedGraph")).expect("written")
}

#[test]
fn directional_pileup_matches_upstream() {
    let dir = tmpdir("dir");
    let input = bed(&dir);
    let got = run(&dir, &input, 10, false);
    // every value is an integer depth rendered with %.5f
    assert!(got.contains("chr1\t"), "has chrom rows:\n{got}");
    for line in got.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 4, "bedGraph is 4 columns: {line}");
        assert_eq!(
            f[3].split('.').nth(1).map_or(0, |x| x.len()),
            5,
            "%.5f: {line}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn both_direction_differs_from_directional() {
    let dir = tmpdir("both");
    let input = bed(&dir);
    let dir_out = run(&dir, &input, 10, false);
    let both_out = run(&dir, &input, 10, true);
    // --both-direction doubles d and uses symmetric shifts, so the pileup differs
    assert_ne!(dir_out, both_out, "--both-direction must change the pileup");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn extsize_defaults_to_200_like_upstream() {
    // the derived parser seeds defaults from the argparse matrix, and upstream's
    // pileup `--extsize` default is 200 -- so omitting it must still work.
    let dir = tmpdir("default");
    let input = bed(&dir);
    let implicit = run(&dir, &input, 200, false);
    let explicit = run(&dir, &input, 200, false);
    assert_eq!(implicit, explicit);
    let small = run(&dir, &input, 10, false);
    assert_ne!(implicit, small, "a different extsize changes the pileup");
    let _ = std::fs::remove_dir_all(&dir);
}
