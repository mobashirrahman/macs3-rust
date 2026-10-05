//! Regression: `bdgopt`'s `-p` operand is a Python `float`, applied in `f64`.
//!
//! `bdgopt_cmd.py:47-53` builds `extraparam = float(options.extraparam[0])` -- a
//! `double` -- and hands it to `bedGraphTrackI.apply_func`
//! (`BedGraph.py:888-901`), which calls `func(s[i])` on each element of a
//! `array.array('f')` (`BedGraph.py:138-139`, `182-183`). An `array.array`
//! element indexes as a **Python float**, not a numpy scalar, so `x * extraparam`
//! is a `double` multiply and only the store back into the `f` array narrows.
//!
//! Doing the multiply in `f32` instead rounds the operand to `f32` *first*,
//! which is a different number and not merely a wider final rounding: with
//! `x = 48.0` and `-p 1.7` upstream prints `81.60000` (`f32(48.0 * 1.7)`) where
//! the `f32` product `48.0f32 * 1.7000000476837158f32` rounds to `81.60001`.
//! 227360 of 8990989 rows of `macs3 callpeak --bdg` differ that way.
//!
//! The expected bytes are upstream's for this exact input
//! (`macs3 bdgopt -i in.bdg -m multiply -p 1.7 -o o.bdg`).

use macs_cli::{parse_flags, Options};

fn s(x: &str) -> String {
    x.to_string()
}

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-bdgopt-f64-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn bdgopt(dir: &std::path::Path, input: &str, method: &str, param: &str) -> String {
    let out = dir.join(format!("o_{method}.bdg"));
    let args = vec![
        s("-i"),
        s(input),
        s("-m"),
        s(method),
        s("-p"),
        s(param),
        s("-o"),
        s(&format!("o_{method}.bdg")),
        s("--outdir"),
        s(dir.to_str().unwrap()),
    ];
    let o: Options = parse_flags("bdgopt", &args).expect("args parse");
    macs_cli::commands::bedgraph_cmds::bdgopt(&o).expect("bdgopt runs");
    std::fs::read_to_string(&out).expect("output written")
}

#[test]
fn bdgopt_multiply_applies_p_in_double() {
    let dir = tmpdir("mult");
    let input = dir.join("in.bdg");
    std::fs::write(
        &input,
        "chr1\t0\t10\t48.00000\nchr1\t10\t20\t53.00000\nchr1\t20\t30\t58.00000\n",
    )
    .unwrap();
    let got = bdgopt(&dir, input.to_str().unwrap(), "multiply", "1.7");
    let want = "track type=bedGraph name=\"MULTIPLY_modified_scores\" description=\"Scores calculated by MULTIPLY\" visibility=2 alwaysZero=on\n\
                chr1\t0\t10\t81.60000\n\
                chr1\t10\t20\t90.10000\n\
                chr1\t20\t30\t98.60000\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn bdgopt_add_applies_p_in_double() {
    // `f32` addition of two `f32`s is itself correctly rounded, so this case
    // cannot separate the two readings; it is pinned because `-p` is parsed as a
    // `f64` for every method and a future `f32` parse of `-p` would show up here
    // first. Upstream: `macs3 bdgopt -i in.bdg -m add -p 0.1 -o o.bdg`.
    let dir = tmpdir("add");
    let input = dir.join("in.bdg");
    std::fs::write(&input, "chr1\t0\t10\t2.50000\nchr1\t10\t20\t0.00000\n").unwrap();
    let got = bdgopt(&dir, input.to_str().unwrap(), "add", "0.1");
    let want = "track type=bedGraph name=\"ADD_modified_scores\" description=\"Scores calculated by ADD\" visibility=2 alwaysZero=on\n\
                chr1\t0\t10\t2.60000\n\
                chr1\t10\t20\t0.10000\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn bdgopt_max_and_min_still_narrow_p() {
    // `max`/`min` do no arithmetic: upstream compares in `f64` and stores either
    // the untouched `f32` element or the `f64` operand narrowed once. `-p 1.7`
    // is not representable in `f32`, so this pins that the comparison happens
    // before the narrowing. Upstream: `-m max`/`-m min` with `-p 1.7`.
    let dir = tmpdir("clamp");
    let input = dir.join("in.bdg");
    std::fs::write(
        &input,
        "chr1\t0\t10\t1.69999\nchr1\t10\t20\t1.70000\nchr1\t20\t30\t1.70001\n",
    )
    .unwrap();
    let got_max = bdgopt(&dir, input.to_str().unwrap(), "max", "1.7");
    let want_max = "track type=bedGraph name=\"MAX_modified_scores\" description=\"Scores calculated by MAX\" visibility=2 alwaysZero=on\n\
                   chr1\t0\t10\t1.70000\n\
                   chr1\t10\t20\t1.70000\n\
                   chr1\t20\t30\t1.70001\n";
    assert_eq!(got_max, want_max);
    let got_min = bdgopt(&dir, input.to_str().unwrap(), "min", "1.7");
    let want_min = "track type=bedGraph name=\"MIN_modified_scores\" description=\"Scores calculated by MIN\" visibility=2 alwaysZero=on\n\
                   chr1\t0\t10\t1.69999\n\
                   chr1\t10\t20\t1.70000\n\
                   chr1\t20\t30\t1.70000\n";
    assert_eq!(got_min, want_min);
    let _ = std::fs::remove_dir_all(&dir);
}
