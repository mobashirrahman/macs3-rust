//! The `cutoff_analysis` score column is read through a C float.
//!
//! `bedGraphTrackI.cutoff_analysis` builds the ladder as
//! `[round(value, 3) for value in np.arange(minv, maxv, s)]` -- f64 entries -- but
//! its loop variable is declared `cutoff: cython.float` (`BedGraph.py:1281`), so both
//! `cutoff = cutoff_list[n]` statements (`BedGraph.py:1323`, `BedGraph.py:1375`)
//! assign into a C float (`BedGraph.c:29325`, `BedGraph.c:30006`). The comparison
//! `score_array > cutoff` and the `"%.2f"` of the report therefore both see the **f32**
//! image of the rounded cutoff, never the f64 one.
//!
//! That is load-bearing, because the ladder step is not a round decimal. With a track
//! spanning `0 .. 2.0308969020843506` (the fold-change track of
//! `tests/fixtures/frag_basic/barcode_fragments/treat.frag`) and `steps = 100`, the
//! step is `0.020308969542384148` and the 3-decimal rounding lands seven entries on a
//! `.xx5` boundary. Six of them sit *below* the decimal midpoint as f64 and *above* it
//! as f32, so `%.2f` disagrees:
//!
//! | ladder entry | `%.2f` of the f64 | `%.2f` of the f32 |
//! |--------------|-------------------|-------------------|
//! | 1.015        | `1.01`            | `1.01`            |
//! | 0.995        | `0.99`            | `1.00`            |
//! | 0.975        | `0.97`            | `0.98`            |
//! | 0.955        | `0.95`            | `0.95`            |
//! | 0.345        | `0.34`            | `0.34`            |
//! | 0.325        | `0.33`            | `0.32`            |
//! | 0.305        | `0.30`            | `0.31`            |
//!
//! Printing the f64 ladder turns upstream's `1.00 0.98 0.32 0.31` into
//! `0.99 0.97 0.33 0.30`, which is the whole of the byte difference against the
//! reference `*_cutoff_analysis.tsv` on that fixture.
//!
//! The step is narrowed as well: `s = float(maxv - minv)/steps` is stored back into the
//! `cython.float s` of `BedGraph.py:1309` (`BedGraph.c:28941`), so `np.arange` is handed
//! `f32(f64(span) / steps)` rather than the f64 quotient.
//!
//! No Python and no reference binary is involved: the track below reproduces upstream's
//! ladder from its own min and max, and the expected rows are the ones upstream prints.

use macs_bedgraph::BedGraph;

/// The fold-change maximum of `barcode_fragments/treat.frag`, as the reference reports
/// it: a C float, so this f64 literal is the exact value of `np.float32(7989 /
/// 3933.72998046875)` and `MAXV as f32` is the track's own `maxvalue`.
const MAXV: f64 = 2.030_896_902_084_350_6;

/// One chromosome whose value extremes are the fixture's, so that `minv = 0`,
/// `maxv = 2.0308969020843506` and `steps = 100` gives bit-for-bit the ladder
/// `hmmratac` builds for `barcode_fragments/treat.frag`.
///
/// Only the `[1, 102)` run clears the cutoffs between 0.3 and 1.02, and the `[150, 151)`
/// run -- the track's maximum -- clears everything below 2.03. It is bracketed by a
/// below-cutoff run, so neither above-cutoff run is the track's *first* run and the
/// `pos[-1]` wrap-around in `above_cutoff_startpos` (`BedGraph.py:1326`) does not fire.
fn ladder_track() -> BedGraph {
    let mut bg = BedGraph::new(0.0);
    let runs: &[(u32, f32)] = &[(1, 0.0), (102, 1.05), (150, 0.0), (151, MAXV as f32)];
    let mut pre: u32 = 0;
    for &(end, value) in runs {
        bg.add_loc(b"chr1", pre, end, value);
        pre = end;
    }
    bg
}

/// `cutoff_analysis(min_length=98, max_gap=1000, steps=100, max_score=100)`: the
/// arguments `hmmratac_cmd.py:307` passes, with `minlen` from the fixture's average
/// fragment length (98) and `max_gap` from `--training-flanking`'s default (1000).
fn report() -> String {
    macs_bedgraph::peakcall::cutoff_analysis(&ladder_track(), 1000, 98, 100, 0.0, 100.0)
}

fn row(text: &str, score: &str) -> String {
    text.lines()
        .skip(1)
        .find(|l| l.starts_with(&format!("{score}\t")))
        .unwrap_or_else(|| panic!("no row for score {score} in {text:?}"))
        .to_string()
}

#[test]
fn a_half_centimal_ladder_entry_prints_as_its_f32_image() {
    let text = report();
    let scores: Vec<&str> = text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    // The rows the bug moved, and the neighbours that straddle them. A ladder that is
    // off by a whole step would lose the neighbours too, so both are asserted.
    for want in [
        "1.01", "1.00", "0.98", "0.95", "0.34", "0.32", "0.31", "0.28",
    ] {
        assert!(
            scores.contains(&want),
            "missing {want} in the score column: {scores:?}"
        );
    }
    // Rows upstream prints for this track: the two above-cutoff runs merge (the 48 bp
    // gap is inside `max_gap`) into one 150 bp peak.
    assert_eq!(row(&text, "1.00"), "1.00\t1\t150\t150.00");
    assert_eq!(row(&text, "0.98"), "0.98\t1\t150\t150.00");
    assert_eq!(row(&text, "0.32"), "0.32\t1\t150\t150.00");
    assert_eq!(row(&text, "0.31"), "0.31\t1\t150\t150.00");
}

/// The score column must descend and must stay inside `(minv, maxv]`: `np.arange` is
/// half-open, so `maxv` itself is never a cutoff, and rows are emitted top-down.
#[test]
fn the_score_column_descends_and_never_reaches_maxv() {
    let text = report();
    let scores: Vec<f32> = text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap().parse::<f32>().unwrap())
        .collect();
    assert!(
        scores.len() > 50,
        "expected a full ladder, got {} rows",
        scores.len()
    );
    for w in scores.windows(2) {
        assert!(w[0] > w[1], "the report must descend: {w:?}");
        // The ladder step is 0.0203 and each entry is printed to two decimals, so a
        // gap between printed rows is 0.02 or 0.03 (and 0.01 across a `xx5` rounding).
        assert!(
            (0.005..=0.0305).contains(&(w[0] - w[1])),
            "gap of {} between {w:?}",
            w[0] - w[1]
        );
    }
    assert!(
        scores[0] < MAXV as f32,
        "top row {} reached maxv",
        scores[0]
    );
}

/// The ladder itself: `ceil((maxv - minv) / s)` entries, `s` stored back into a C float,
/// and every entry narrowed to f32 before it is used.
#[test]
fn the_ladder_is_half_open_and_stored_as_f32() {
    let minv = 0.0f32;
    let maxv = MAXV as f32;
    let step = (f64::from(maxv - minv) / 100.0) as f32;
    assert_ne!(
        f64::from(step),
        f64::from(maxv - minv) / 100.0,
        "the f32 step must differ from the f64 quotient, else this test proves nothing"
    );
    let n = ((f64::from(maxv) - f64::from(minv)) / f64::from(step)).ceil() as usize;
    assert_eq!(n, 100, "one entry per step, with maxv excluded");

    let mut last: Option<f32> = None;
    for i in 0..n {
        // `round(value, 3)`, then `cutoff = cutoff_list[n]` into a C float.
        let entry: f32 =
            macs_stats::py_round(f64::from(minv) + i as f64 * f64::from(step), 3) as f32;
        if i + 1 < n {
            assert!(entry < maxv, "entry {i} = {entry} reaches maxv = {maxv}");
        }
        if let Some(p) = last {
            assert!(entry > p, "the ladder must increase: {entry} after {p}");
        }
        last = Some(entry);
    }
}

/// The seven entries whose f64 and f32 images disagree at two decimals. `1.015` is the
/// control: it agrees, so it is not a probe, and only the other six can move the file.
#[test]
fn six_of_the_seven_boundary_entries_disagree_between_f32_and_f64() {
    for v in [0.995f64, 0.975, 0.325, 0.305] {
        assert_ne!(
            format!("{v:.2}"),
            format!("{:.2}", v as f32),
            "{v}: the two widths would agree, this test would be stale"
        );
    }
    for v in [1.015f64, 0.955, 0.345] {
        assert_eq!(
            format!("{v:.2}"),
            format!("{:.2}", v as f32),
            "{v}: expected to agree in both widths"
        );
    }
}
