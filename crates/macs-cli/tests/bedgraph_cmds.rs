//! End-to-end tests for the bedGraph commands (`bdgopt`, `cmbreps`).
//!
//! The expected bytes were produced by *upstream* MACS3 3.0.5
//! (`bedGraphTrackI.apply_func` / `.overlie` + `write_bedGraph`) and are pinned
//! here so a regression in the parser, the library arithmetic or the `%.5f`
//! writer fails in CI rather than silently drifting.

use macs_cli::{parse_flags, Options};

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-bdgc-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn write(dir: &std::path::Path, name: &str, body: &str) -> String {
    let p = dir.join(name);
    std::fs::write(&p, body).unwrap();
    p.to_string_lossy().into_owned()
}

fn parse(sub: &str, args: &[String]) -> Options {
    parse_flags(sub, args).expect("args parse")
}

fn s(x: &str) -> String {
    x.to_string()
}

#[test]
fn bdgopt_methods_are_byte_identical_to_upstream() {
    let dir = tmpdir("opt");
    let input = write(&dir, "in.bdg", "chr1\t0\t100\t1.5\nchr1\t100\t200\t2.5\n");
    // upstream output for each method (pinned)
    let cases: &[(&str, &str, &str)] = &[
        ("multiply", "2", "track type=bedGraph name=\"MULTIPLY_modified_scores\" description=\"Scores calculated by MULTIPLY\" visibility=2 alwaysZero=on\nchr1\t0\t100\t3.00000\nchr1\t100\t200\t5.00000\n"),
        ("add", "1", "track type=bedGraph name=\"ADD_modified_scores\" description=\"Scores calculated by ADD\" visibility=2 alwaysZero=on\nchr1\t0\t100\t2.50000\nchr1\t100\t200\t3.50000\n"),
        ("max", "2", "track type=bedGraph name=\"MAX_modified_scores\" description=\"Scores calculated by MAX\" visibility=2 alwaysZero=on\nchr1\t0\t100\t2.00000\nchr1\t100\t200\t2.50000\n"),
        ("min", "2", "track type=bedGraph name=\"MIN_modified_scores\" description=\"Scores calculated by MIN\" visibility=2 alwaysZero=on\nchr1\t0\t100\t1.50000\nchr1\t100\t200\t2.00000\n"),
    ];
    for (method, param, want) in cases {
        let out = dir.join(format!("o_{method}.bdg"));
        let args = vec![
            s("-i"),
            s(&input),
            s("-m"),
            s(method),
            s("-p"),
            s(param),
            s("-o"),
            s(&format!("o_{method}.bdg")),
            s("--outdir"),
            s(dir.to_str().unwrap()),
        ];
        let o = parse("bdgopt", &args);
        macs_cli::commands::bedgraph_cmds::bdgopt(&o).expect("bdgopt runs");
        let got = std::fs::read_to_string(&out).expect("output written");
        assert_eq!(&got, want, "bdgopt --method {method}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn bdgopt_p2q_is_byte_identical_to_upstream() {
    let dir = tmpdir("p2q");
    let input = write(&dir, "in.bdg", "chr1\t0\t100\t1.5\nchr1\t100\t200\t2.5\n");
    // upstream p2q: q = clamp(p - log10(N)) with N=200, so log10(200)~2.301
    let out = dir.join("o_p2q.bdg");
    let args = vec![
        s("-i"),
        s(&input),
        s("-m"),
        s("p2q"),
        s("-o"),
        s("o_p2q.bdg"),
        s("--outdir"),
        s(dir.to_str().unwrap()),
    ];
    let o = parse("bdgopt", &args);
    macs_cli::commands::bedgraph_cmds::bdgopt(&o).expect("bdgopt p2q runs");
    let got = std::fs::read_to_string(&out).expect("written");
    // the exact upstream bytes for this input:
    let want = "track type=bedGraph name=\"P2Q_modified_scores\" description=\"Scores calculated by P2Q\" visibility=2 alwaysZero=on\nchr1\t0\t100\t0.00000\nchr1\t100\t200\t0.19897\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn cmbreps_methods_are_byte_identical_to_upstream() {
    let dir = tmpdir("cmb");
    let a = write(
        &dir,
        "a.bdg",
        "chr1\t0\t150\t1\nchr1\t150\t250\t2\nchr1\t250\t300\t4\n",
    );
    let b = write(
        &dir,
        "b.bdg",
        "chr1\t0\t100\t0\nchr1\t100\t200\t3\nchr1\t200\t300\t4\n",
    );
    let want_max = "track type=bedGraph name=\"MAX_combined_scores\" description=\"Scores calculated by MAX\" visibility=2 alwaysZero=on\nchr1\t0\t100\t1.00000\nchr1\t100\t200\t3.00000\nchr1\t200\t300\t4.00000\n";
    let out = dir.join("cmb.bdg");
    let args = vec![
        s("-i"),
        s(&a),
        s(&b),
        s("-m"),
        s("max"),
        s("-o"),
        s("cmb.bdg"),
        s("--outdir"),
        s(dir.to_str().unwrap()),
    ];
    let o = parse("cmbreps", &args);
    macs_cli::commands::bedgraph_cmds::cmbreps(&o).expect("cmbreps runs");
    assert_eq!(std::fs::read_to_string(&out).unwrap(), want_max);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn cmbreps_requires_two_inputs() {
    let dir = tmpdir("cmb1");
    let a = write(&dir, "a.bdg", "chr1\t0\t100\t1\n");
    let args = vec![s("-i"), s(&a), s("-m"), s("max"), s("-o"), s("o.bdg")];
    let o = parse("cmbreps", &args);
    let e = macs_cli::commands::bedgraph_cmds::cmbreps(&o).expect_err("one input is rejected");
    assert!(format!("{e}").contains("at least two"), "got {e}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_multi_value_flag_keeps_every_value() {
    // cmbreps `-i` is `nargs +`; the derived parser must retain all of them.
    let o = parse(
        "cmbreps",
        &[
            s("-i"),
            s("a.bdg"),
            s("b.bdg"),
            s("c.bdg"),
            s("-m"),
            s("mean"),
            s("-o"),
            s("o.bdg"),
        ],
    );
    assert_eq!(o.get_all("ifile"), ["a.bdg", "b.bdg", "c.bdg"]);
}
