//! **F180/F168** — the paired treatment/control arrays are
//! `__chrom_pair_treat_ctrl`'s **union** walk.
//!
//! The p-score track and the AFDR histogram are both built on that array. Building
//! them on the positions where both tracks happen to change instead loses rows and
//! removes the negative first span, which is what makes the p->q walk produce its
//! `NaN` window (F179). See `docs/upstream-findings.md`.

use macs_core::ChromId;
use macs_rle::SignalTrack;

const C: ChromId = ChromId(0);

fn track(runs: &[(u64, f32)]) -> SignalTrack<f32> {
    let mut t = SignalTrack::empty(C, 0, runs.last().map(|r| r.0).unwrap_or(0));
    for &(end, v) in runs {
        t.push(end, v);
    }
    t
}

/// **F180** — the paired arrays are `__chrom_pair_treat_ctrl`'s **union** walk: one
/// row per union position, carrying whatever each cursor currently holds.
///
/// The treatment and control tracks here break at disjoint coordinates, so the
/// union has strictly more rows than the set of coincident boundaries.
#[test]
fn the_union_walk_emits_one_row_per_union_position() {
    let treat = track(&[(10, 100.0), (20, 5.0)]);
    let ctrl = track(&[(15, 1.0), (25, 2.0)]);
    let (pos, t, c) = macs_peaks::callpeak::paired_union(&treat, &ctrl);

    // upstream: 10 < 15 -> (10, 100, 1); 20 > 15 -> (15, 5, 1); 20 < 25 -> (20, 5, 2)
    assert_eq!(pos, vec![10, 15, 20]);
    assert_eq!(t, vec![100.0, 5.0, 5.0]);
    assert_eq!(c, vec![1.0, 1.0, 2.0]);
}

/// The treatment's value is carried forward across every position the *control*
/// contributes, and vice versa -- that repetition is what the histogram weights.
#[test]
fn the_union_walk_repeats_a_cursor_value_across_the_other_tracks_boundaries() {
    let treat = track(&[(10, 7.0), (12, 9.0)]);
    let ctrl = track(&[(11, 3.0), (14, 4.0)]);
    let (pos, t, c) = macs_peaks::callpeak::paired_union(&treat, &ctrl);
    // position 14 never appears: the treatment is exhausted first
    assert_eq!(pos, vec![10, 11, 12]);
    // the treatment's 9 is carried across the control boundary at 11; at 12 the
    // *control* cursor has advanced, so that row carries the control's 4
    assert_eq!(t, vec![7.0, 9.0, 9.0]);
    assert_eq!(c, vec![3.0, 3.0, 4.0]);
}

/// The walk stops as soon as either input is exhausted (F168), so the longer
/// array's tail never appears.
#[test]
fn the_union_walk_stops_with_the_shorter_input() {
    let treat = track(&[(10, 1.0)]);
    let ctrl = track(&[(5, 2.0), (20, 3.0)]);
    let (pos, _, _) = macs_peaks::callpeak::paired_union(&treat, &ctrl);
    assert_eq!(pos, vec![5, 10]);
}
