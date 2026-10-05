//! `hmmratac`'s `*_cutoff_analysis.tsv`, byte for byte, on the fixture that exposed the
//! f32/f64 ladder bug.
//!
//! `hmmratac_cmd.py:307` calls `bedGraphTrackI.cutoff_analysis`, whose score column is
//! printed with `"%.2f"` applied to a `cython.float` (`BedGraph.py:1281`), i.e. to the
//! **f32** image of each `round(value, 3)` ladder entry. On this fixture the ladder step
//! is `0.020308969542384148` and four entries land on a `.xx5` decimal boundary, so the
//! f32 and f64 renderings disagree: the reference prints `1.00 0.98 0.32 0.31` where an
//! f64 ladder prints `0.99 0.97 0.33 0.30`.
//!
//! The expected text below is the reference's own file for
//! `tests/fixtures/frag_basic/barcode_fragments/treat.frag`. The fixture is committed,
//! so this test needs neither Python nor the reference binary.
//!
//! Both implementations stop with "Not enough training regions!" -- this fragment file
//! calls no training region at the default `-l 10 -u 20` -- and both write the report
//! before stopping, which is what makes the comparison possible at all.

use macs_cli::{parse_flags, Options};

/// Upstream's `barcode_fragments` report: header plus 73 rows, 100 ladder entries of
/// which the 27 highest call nothing.
const REFERENCE: &str = concat!(
    "score	npeaks	lpeaks	avelpeak\n",
    "1.46	1	98	98.00\n",
    "1.44	1	98	98.00\n",
    "1.42	1	98	98.00\n",
    "1.40	1	98	98.00\n",
    "1.38	1	98	98.00\n",
    "1.36	1	98	98.00\n",
    "1.34	1	98	98.00\n",
    "1.32	1	98	98.00\n",
    "1.30	1	99	99.00\n",
    "1.28	1	99	99.00\n",
    "1.26	1	99	99.00\n",
    "1.24	1	99	99.00\n",
    "1.22	1	99	99.00\n",
    "1.20	1	99	99.00\n",
    "1.18	1	99	99.00\n",
    "1.16	1	99	99.00\n",
    "1.14	1	100	100.00\n",
    "1.12	1	100	100.00\n",
    "1.10	1	100	100.00\n",
    "1.08	1	100	100.00\n",
    "1.06	1	100	100.00\n",
    "1.04	1	100	100.00\n",
    "1.01	1	101	101.00\n",
    "1.00	1	101	101.00\n",
    "0.98	1	101	101.00\n",
    "0.95	1	101	101.00\n",
    "0.93	1	101	101.00\n",
    "0.91	1	101	101.00\n",
    "0.89	1	101	101.00\n",
    "0.87	1	101	101.00\n",
    "0.85	1	102	102.00\n",
    "0.83	1	102	102.00\n",
    "0.81	1	102	102.00\n",
    "0.79	1	102	102.00\n",
    "0.77	1	102	102.00\n",
    "0.75	1	102	102.00\n",
    "0.73	1	102	102.00\n",
    "0.71	1	103	103.00\n",
    "0.69	1	103	103.00\n",
    "0.67	1	103	103.00\n",
    "0.65	1	103	103.00\n",
    "0.63	1	103	103.00\n",
    "0.61	1	103	103.00\n",
    "0.59	1	103	103.00\n",
    "0.57	1	104	104.00\n",
    "0.55	1	104	104.00\n",
    "0.53	1	104	104.00\n",
    "0.51	1	104	104.00\n",
    "0.49	1	104	104.00\n",
    "0.47	1	104	104.00\n",
    "0.45	1	104	104.00\n",
    "0.43	1	105	105.00\n",
    "0.41	1	105	105.00\n",
    "0.39	1	105	105.00\n",
    "0.37	1	105	105.00\n",
    "0.34	1	105	105.00\n",
    "0.32	1	106	106.00\n",
    "0.31	1	106	106.00\n",
    "0.28	1	106	106.00\n",
    "0.26	1	106	106.00\n",
    "0.24	1	106	106.00\n",
    "0.22	1	107	107.00\n",
    "0.20	1	107	107.00\n",
    "0.18	1	107	107.00\n",
    "0.16	1	107	107.00\n",
    "0.14	1	108	108.00\n",
    "0.12	1	108	108.00\n",
    "0.10	1	109	109.00\n",
    "0.08	1	109	109.00\n",
    "0.06	1	110	110.00\n",
    "0.04	1	111	111.00\n",
    "0.02	1	113	113.00\n",
    "0.00	1	118	118.00\n",
);

fn fixture() -> String {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/frag_basic/barcode_fragments/treat.frag")
        .to_str()
        .unwrap()
        .to_string()
}

fn s(value: &str) -> String {
    value.to_string()
}

#[test]
fn the_cutoff_analysis_report_matches_the_reference_byte_for_byte() {
    let out = std::env::temp_dir().join(format!("macs3rs-hmmratac-cutoff-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        fixture(),
        s("-f"),
        s("FRAG"),
        s("-n"),
        s("p"),
        s("--outdir"),
        s(out.to_str().unwrap()),
    ];
    let options: Options = parse_flags("hmmratac", &args).expect("valid options");
    let err = macs_cli::commands::hmmratac::hmmratac(&options)
        .expect_err("no training region is called at the default cutoffs");
    assert!(
        err.to_string().contains("Not enough training regions"),
        "unexpected error: {err}"
    );

    let got = std::fs::read_to_string(out.join("p_cutoff_analysis.tsv"))
        .expect("the report is written before the error is raised");
    assert_eq!(got, REFERENCE, "cutoff analysis report differs");
}

/// The rows the bug moved, spelled out so a failure names them instead of dumping a
/// 74-line diff.
#[test]
fn the_moved_rows_read_as_the_f32_image() {
    let out = std::env::temp_dir().join(format!("macs3rs-hmmratac-cutoff2-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        fixture(),
        s("-f"),
        s("FRAG"),
        s("-n"),
        s("p"),
        s("--outdir"),
        s(out.to_str().unwrap()),
    ];
    let options: Options = parse_flags("hmmratac", &args).expect("valid options");
    let _ = macs_cli::commands::hmmratac::hmmratac(&options);
    let got = std::fs::read_to_string(out.join("p_cutoff_analysis.tsv")).expect("report written");
    for want in [
        "1.00\t1\t101\t101.00",
        "0.98\t1\t101\t101.00",
        "0.32\t1\t106\t106.00",
        "0.31\t1\t106\t106.00",
    ] {
        assert!(got.lines().any(|l| l == want), "missing row {want:?}");
    }
    for wrong in [
        "0.99\t1\t101\t101.00",
        "0.97\t1\t101\t101.00",
        "0.33\t1\t106\t106.00",
        "0.30\t1\t106\t106.00",
    ] {
        assert!(
            !got.lines().any(|l| l == wrong),
            "row {wrong:?} is the f64 rendering, not the f32 one"
        );
    }
}
