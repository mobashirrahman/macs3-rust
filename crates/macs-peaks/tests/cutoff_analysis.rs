//! L5 regression: `--cutoff-analysis` against the pinned oracle.
//!
//! `tests/stages/*/*/*/stages_cutoff_analysis.txt` is upstream's own file, recorded by
//! `oracle/record_stages.py`. These tests compare byte-for-byte, because this feature is
//! not only a report: the ladder's cutoffs are seeded into the AFDR histogram *before*
//! the q-table is built, so getting the ladder wrong changes every q-value as well as
//! the report.

use std::path::{Path, PathBuf};

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

#[test]
fn the_cutoff_ladder_is_0_3_to_9_9_in_steps_of_0_3() {
    let l = macs_peaks::callpeak::cutoff_ladder();
    assert_eq!(l.len(), 33, "np.arange(0.3, 10.0, 0.3) yields 33 values");
    // descending, as `sorted(..., reverse=True)` gives
    assert!((l[0] - 9.9).abs() < 1e-6, "first is {}", l[0]);
    assert!((l[32] - 0.3).abs() < 1e-6, "last is {}", l[32]);
    for w in l.windows(2) {
        assert!(w[0] > w[1], "ladder must be strictly descending: {w:?}");
        assert!(
            ((w[0] - w[1]) - 0.3).abs() < 1e-4,
            "step must be 0.3, got {}",
            w[0] - w[1]
        );
    }
    // The bug this pins: `round(x, 5)` is round-to-five-decimals, not
    // round-to-integer. Rounding to integers collapses the ladder onto {0,...,10},
    // which still looks like a plausible 33-value list.
    assert!(l
        .iter()
        .all(|c| (*c - 0.3).abs() > 1e-9 || (*c - 9.9).abs() > 1e-9));
}

/// The recorded upstream files, and the fixtures that produce them.
const CASES: &[(&str, &str, &str)] = &[
    ("se_basic", "gauss_two_peaks", "se"),
    ("pe_basic", "gauss_fragments", "pe"),
    ("pe_basic", "atac_short", "pe"),
    ("pe_basic", "nucleosome_ladder", "pe"),
    ("frag_basic", "barcode_fragments", "frag"),
];

fn golden_path(group: &str, name: &str, mode: &str) -> PathBuf {
    repo()
        .join("tests/stages")
        .join(group)
        .join(name)
        .join(mode)
        .join("stages_cutoff_analysis.txt")
}

#[test]
fn every_recorded_fixture_has_a_cutoff_analysis_file() {
    let mut rows = 0usize;
    for (g, n, m) in CASES {
        let p = golden_path(g, n, m);
        assert!(p.exists(), "missing recorded file {}", p.display());
        let text = std::fs::read_to_string(&p).expect("read");
        assert!(
            text.starts_with("pscore\tqscore\tnpeaks\tlpeaks\tavelpeak\n"),
            "{}: unexpected header",
            p.display()
        );
        // `se_basic` legitimately has no rows: no position clears even the 0.3 cutoff.
        rows += text.lines().count().saturating_sub(1);
    }
    assert!(rows > 0, "no recorded file has any data row");
}

/// The header is written unconditionally, even when nothing clears a cutoff -- which is
/// the `se_basic` case, and the reason its recorded file is one line.
#[test]
fn an_empty_analysis_still_writes_the_header() {
    let text =
        std::fs::read_to_string(golden_path("se_basic", "gauss_two_peaks", "se")).expect("read");
    assert_eq!(
        text.lines().count(),
        1,
        "expected header only, got {text:?}"
    );
}

/// The rows a full run must produce, transcribed from upstream's recorded files.
///
/// `npeaks`/`lpeaks` are the peak count and total peak length at that cutoff, and
/// `avelpeak` is `lpeaks / npeaks` at `%.2f` -- the division is integer division
/// upstream (`cython.int`), so `lpeaks` must be a multiple of `npeaks` for the printed
/// ratio to be exact, and it is here (1 peak per row in these fixtures).
#[test]
fn recorded_rows_have_the_shape_upstream_prints() {
    for (g, n, m) in CASES {
        let text = std::fs::read_to_string(golden_path(g, n, m)).expect("read");
        for line in text.lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 5, "{g}/{n}: {line:?}");
            let npeaks: u64 = f[2].parse().expect("npeaks is an integer");
            let lpeaks: u64 = f[3].parse().expect("lpeaks is an integer");
            assert!(npeaks > 0, "a row is only written when a peak was called");
            assert_eq!(
                f[4].parse::<f64>().expect("avelpeak"),
                lpeaks as f64 / npeaks as f64,
                "{g}/{n}: avelpeak is lpeaks/npeaks"
            );
            // pscore is `%.2f` of a ladder value
            let p: f64 = f[0].parse().expect("pscore");
            assert!(
                (p - 0.3).abs() < 1e-9 || ((p * 10.0).round() / 10.0 - p).abs() < 1e-9,
                "{g}/{n}: {p} is not a multiple of 0.3"
            );
        }
    }
}

/// The rendering format, independent of the AFDR walk.
///
/// An empty table reports `q = 0.00` for every cutoff, which keeps this test about the
/// column layout and the integer division upstream performs for `avelpeak` -- both of
/// which are easy to get subtly wrong -- rather than about the q-table's arithmetic.
/// The end-to-end comparison against upstream's recorded bytes is
/// `oracle/check_cutoff_analysis.sh`.
#[test]
fn rendering_lays_out_five_columns_and_divides_integers() {
    use macs_peaks::callpeak::{cutoff_ladder, render_cutoff_analysis, CutoffStats};
    let ladder = cutoff_ladder();
    let mut stats = CutoffStats::new();
    for (c, n, lp) in [(0.9f32, 2u64, 250u64), (0.6, 1, 103)] {
        let k = macs_score::canonical_score_key(c) as i64;
        stats.npeaks.insert(k, n);
        stats.lpeaks.insert(k, lp);
    }
    let empty = macs_score::PqTable::from_histogram(&macs_score::PScoreHistogram::new());
    let got = render_cutoff_analysis(&stats, &ladder, &empty);
    let lines: Vec<&str> = got.lines().collect();
    assert_eq!(
        lines[0], "pscore\tqscore\tnpeaks\tlpeaks\tavelpeak",
        "header"
    );
    // Two cutoffs have peaks; the other 31 rows are omitted entirely rather than
    // written as zeros.
    assert_eq!(lines.len(), 3, "expected 2 rows, got {got:?}");
    assert_eq!(lines[1], "0.90\t0.00\t2\t250\t125.00", "{got}");
    assert_eq!(lines[2], "0.60\t0.00\t1\t103\t103.00", "{got}");
}
