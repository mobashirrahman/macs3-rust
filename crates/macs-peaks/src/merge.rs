//! Merging pileup arrays: the pointwise max that builds local lambda.
//!
//! Transcribed from `MACS3.Signal.PileupV2.over_two_pv_array` and from
//! `FixWidthTrack.pileup_a_chromosome_c`, which combines several control
//! pileups -- one per local-lambda scale -- by taking their pointwise maximum.
//!
//! # The merge is head-to-head, not interval-aligned
//!
//! Upstream walks both arrays with a cursor and emits one output entry per
//! *step*, advancing whichever cursor has the smaller position (both, on a tie).
//! The emitted value is `max(v1[i1], v2[i2])` of the **current heads**, before
//! either cursor moves. It is not a re-sample of the second array at the first's
//! position: when the two arrays have different breakpoints, the output carries
//! each array's own value forward until its cursor advances.
//!
//! # Tail truncation is real
//!
//! The loop is `while i1 < l1 and i2 < l2`. When one array is exhausted, the
//! remainder of the other is **dropped**, not appended. So merging arrays of
//! unequal length returns at most `l1 + l2` entries but often fewer, and the
//! coverage can end before the chromosome does. That is what upstream computes
//! and it changes lambda for everything downstream, so it is reproduced rather
//! than "fixed" by draining the longer array.
//!
//! # Position equality emits one entry, not two
//!
//! On a tie both cursors advance but a single output entry is written. Duplicated
//! breakpoint positions cannot appear in the output.

use macs_core::{ChromId, Coord};
use macs_rle::{Run, SignalTrack};

/// The reducer `over_two_pv_array` applies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reducer {
    /// `max`, used for the local-lambda combination.
    Max,
    /// `min`.
    Min,
    /// Arithmetic mean.
    Mean,
}

impl Reducer {
    fn apply(self, a: f32, b: f32) -> f32 {
        match self {
            // upstream writes `a if a > b else b`, so a tie takes the *second*
            // array's value. With equal floats that is not observable, but the
            // form is kept so the semantics are explicit.
            Reducer::Max => {
                if a > b {
                    a
                } else {
                    b
                }
            }
            Reducer::Min => {
                if a < b {
                    a
                } else {
                    b
                }
            }
            Reducer::Mean => (a + b) / 2.0,
        }
    }
}

/// The head-to-head walk shared by every `over_two_pv_*` entry point.
///
/// Emits `(min(p1, p2), apply(v1, v2))` per step, advancing the smaller cursor
/// (both on a tie), and stops as soon as either cursor is exhausted -- the
/// tail-truncation and tie rules described at the top of the module. Callers only
/// choose what to do with each emitted pair, which is what lets the track
/// builder below skip the legacy `[pos, val]` intermediate without changing a
/// single emitted pair.
pub(crate) fn over_two_pv_walk(
    a_runs: &[Run<f32>],
    b_runs: &[Run<f32>],
    func: Reducer,
    emit: &mut impl FnMut(Coord, f32),
) {
    let (l1, l2) = (a_runs.len(), b_runs.len());
    let (mut i1, mut i2) = (0usize, 0usize);
    while i1 < l1 && i2 < l2 {
        let v1 = a_runs[i1].value;
        let v2 = b_runs[i2].value;
        let value = func.apply(v1, v2);

        let p1 = a_runs[i1].end;
        let p2 = b_runs[i2].end;
        match p1.cmp(&p2) {
            std::cmp::Ordering::Less => {
                emit(p1, value);
                i1 += 1;
            }
            std::cmp::Ordering::Greater => {
                emit(p2, value);
                i2 += 1;
            }
            std::cmp::Ordering::Equal => {
                // one entry, both cursors advance
                emit(p1, value);
                i1 += 1;
                i2 += 1;
            }
        }
    }
}

/// Merge two end-indexed tracks with a pointwise reducer.
///
/// Returns `(positions, values)` in the legacy `[p, v]` form upstream works in.
/// Positions are the **ends** of the runs (F5), and values are `f32`.
pub fn over_two_pv_array(
    a: &SignalTrack<f32>,
    b: &SignalTrack<f32>,
    func: Reducer,
) -> (Vec<Coord>, Vec<f32>) {
    let a_runs = a.runs();
    let b_runs = b.runs();
    let l1 = a_runs.len();
    let l2 = b_runs.len();

    let mut pos: Vec<Coord> = Vec::with_capacity(l1 + l2);
    let mut val: Vec<f32> = Vec::with_capacity(l1 + l2);
    over_two_pv_walk(a_runs, b_runs, func, &mut |p, v| {
        pos.push(p);
        val.push(v);
    });
    // F168: the tail of the longer array is **dropped**, and that is upstream's
    // behaviour -- verified against the compiled `over_two_pv_array` with
    // controlled inputs rather than by reading the source:
    //
    // ```text
    // a=[10,20,30] b=[0,5]       -> n=2  pos=[0,5]
    // a=[10,20,30] b=[0,5,40]    -> n=5  pos=[0,5,10,20,30]
    // a=[0,5]       b=[10,20,30] -> n=2  pos=[0,5]
    // a=[10,20,30] b=[10,40,50]  -> n=3  pos=[10,20,30]
    // ```
    //
    // The walk emits `(min(p1, p2), apply(v1, v2))` and stops as soon as either
    // cursor is exhausted; nothing is appended afterwards. The `equal` case
    // advances **both** cursors, which is what makes case 4 emit three rows
    // rather than one.
    //
    // A real measurement that depends on this: on
    // `frag_basic/barcode_fragments` the `d`(152 rows) x `slocal`(168) merge
    // yields **232** rows spanning `[-409, 249]` -- the `d` scale's extent, with
    // 88 of `slocal`'s right-hand rows discarded -- and the following
    // `llocal`(160) merge yields **304** rows spanning `[-4901, 249]`. Flushing
    // instead gives 320 and 480, reaching the `llocal` extent and changing every
    // control value outside the `d` window.
    (pos, val)
}

/// [`over_two_pv_array`] as a track-to-track `max` merge.
///
/// The merged track's span ends wherever the merge stopped, which can be short of
/// either input's end (F35). The emitted pairs are identical to
/// `over_two_pv_array`'s by construction -- the same walk -- so this and
/// `track_from_pv(chrom, start, pos, val)` return byte-identical tracks.
///
/// It exists as a separate entry point because of allocation, not behaviour. Going
/// through the `[pos, val]` form costs 8 bytes per emitted pair for the temporaries
/// on top of the run vector itself, and it sizes that vector at `l1 + l2` when the
/// walk stops at the first exhausted cursor. The local-lambda chain folds three
/// scales per chromosome with the accumulator and the next scale both live, so on
/// chr1 of the 5 M-read fixture those transients were the single-end peak RSS
/// driver. Here the walk runs twice -- count, then fill -- and the run vector is
/// sized exactly and never grows.
pub fn over_two_pv_array_track(a: &SignalTrack<f32>, b: &SignalTrack<f32>) -> SignalTrack<f32> {
    let a_runs = a.runs();
    let b_runs = b.runs();
    let chrom = b.chrom();
    let start = a.start().min(b.start());
    let mut n = 0usize;
    over_two_pv_walk(a_runs, b_runs, Reducer::Max, &mut |_, _| n += 1);
    if n == 0 {
        return SignalTrack::empty(chrom, start, start);
    }
    let mut runs: Vec<Run<f32>> = Vec::with_capacity(n);
    over_two_pv_walk(a_runs, b_runs, Reducer::Max, &mut |p, v| {
        runs.push(Run::new(p, v));
    });
    debug_assert_eq!(runs.len(), n);
    let end = runs[runs.len() - 1].end.max(start);
    SignalTrack::from_runs_exact(chrom, start, end, runs)
}

/// Pointwise max-fold of 2..=3 tracks in a **single pass**.
///
/// `over_two_pv_array_track` folds pairwise, which materialises every intermediate:
/// for three control scales it holds `combined`, the next scale's track and the
/// merge result simultaneously. That transient is what sets the single-end peak RSS
/// (~114 MB for chr1 of the 5 M-read fixture, against a final control track of
/// 2.5 M runs = 20 MB).
///
/// The fold is equivalent to `((a max b) max c)`, which is what upstream's loop
/// computes:
///
/// * the emitted position is `min` of the live cursors, and the emitted value is
///   `max` of the live values -- and `max(max(va,vb),vc) == max(va,vb,vc)`;
/// * `over_two_pv_array` stops as soon as either cursor is exhausted (F168), so
///   `(a max b)` covers `min(extent_a, extent_b)` and the next fold covers
///   `min(that, extent_c)` == `min` of all three. Stopping when *any* of the three
///   exhausts gives the same row count;
/// * the equal-end case advances **both** cursors in the pairwise form and every
///   cursor tied at the minimum in this one, which is the same rule.
///
/// Verified against the pairwise fold by the golden gate (byte identity over 7204
/// recorded invocations), which is what makes the equivalence a fact rather than an
/// argument.
pub fn over_max_tracks(tracks: &[&SignalTrack<f32>]) -> Option<SignalTrack<f32>> {
    match tracks {
        [] => return None,
        [only] => return Some((*only).clone()),
        _ => {}
    }
    let runs: Vec<&[Run<f32>]> = tracks.iter().map(|t| t.runs()).collect();
    let mut idx = vec![0usize; runs.len()];
    let mut out: Vec<Run<f32>> =
        Vec::with_capacity(runs.iter().map(|r| r.len()).min().unwrap_or(0));
    let chrom = tracks[0].chrom();
    let start = tracks.iter().map(|t| t.start()).min().unwrap_or(0);
    loop {
        // Stop when any cursor is exhausted, matching `over_two_pv_array`'s
        // "stop when either runs out" applied through the fold.
        if idx.iter().zip(&runs).any(|(i, r)| *i >= r.len()) {
            break;
        }
        let lowest = idx
            .iter()
            .zip(&runs)
            .map(|(i, r)| r[*i].end)
            .min()
            .expect("non-empty");
        let mut v = f32::MIN;
        for (i, r) in idx.iter().zip(&runs) {
            let cand = r[*i].value;
            // `a if a > b else b`: a tie takes the later array's value, which for
            // equal floats is unobservable but keeps the fold order explicit.
            v = if v > cand { v } else { cand };
        }
        out.push(Run::new(lowest, v));
        for (i, r) in idx.iter_mut().zip(&runs) {
            if r[*i].end == lowest {
                *i += 1;
            }
        }
    }
    let end = out.last().map_or(start, |r| r.end).max(start);
    Some(SignalTrack::from_runs_exact(chrom, start, end, out))
}

/// Rebuild a [`SignalTrack`] from a merged `[p, v]` pair.
///
/// The merged positions are run ends, so the run for entry `i` spans
/// `[p[i-1], p[i])` with `p[-1]` treated as the track's start. The merge flushes
/// the longer input's tail (F168), so the rebuilt track covers the full union of
/// the two inputs' extents.
pub fn track_from_pv(chrom: ChromId, start: Coord, pos: &[Coord], val: &[f32]) -> SignalTrack<f32> {
    debug_assert_eq!(pos.len(), val.len(), "positions and values must pair");
    if pos.is_empty() {
        return SignalTrack::empty(chrom, start, start);
    }
    let end = pos[pos.len() - 1];
    let mut runs: Vec<Run<f32>> = Vec::with_capacity(pos.len());
    for (i, p) in pos.iter().enumerate() {
        runs.push(Run::new(*p, val[i]));
    }
    SignalTrack::from_runs_exact(chrom, start, end.max(start), runs)
}

/// Pointwise maximum of several tracks, folded left to right.
///
/// Upstream combines the per-scale local lambda tracks with
/// `over_two_pv_array(prev, tmp, func="max")` in a loop, so the fold order is
/// `((a max b) max c)`, not a simultaneous maximum over all of them. With a
/// commutative reducer that is unobservable, but the fold is written this way to
/// match.
pub fn pointwise_max(tracks: &[SignalTrack<f32>]) -> SignalTrack<f32> {
    let mut it = tracks.iter();
    let Some(first) = it.next() else {
        panic!("pointwise_max needs at least one track");
    };
    let mut acc = first.clone();
    for next in it {
        let (pos, val) = over_two_pv_array(&acc, next, Reducer::Max);
        acc = track_from_pv(next.chrom(), acc.start(), &pos, &val);
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use macs_rle::SignalTrack;

    fn t(chrom: ChromId, start: Coord, ends: &[Coord], vals: &[f32]) -> SignalTrack<f32> {
        let runs: Vec<Run<f32>> = ends
            .iter()
            .zip(vals)
            .map(|(e, v)| Run::new(*e, *v))
            .collect();
        SignalTrack::from_runs(chrom, start, *ends.last().unwrap(), runs)
    }

    #[test]
    fn equal_lengths_align_one_to_one() {
        let a = t(ChromId(0), 0, &[10, 20], &[1.0, 2.0]);
        let b = t(ChromId(0), 0, &[10, 20], &[5.0, 3.0]);
        let (pos, val) = over_two_pv_array(&a, &b, Reducer::Max);
        assert_eq!(pos, vec![10, 20]);
        assert_eq!(val, vec![5.0, 3.0]);
    }

    #[test]
    fn equal_positions_emit_one_entry_and_advance_both() {
        let a = t(ChromId(0), 0, &[10, 20], &[1.0, 2.0]);
        let b = t(ChromId(0), 0, &[10, 20], &[4.0, 6.0]);
        let (pos, val) = over_two_pv_array(&a, &b, Reducer::Max);
        assert_eq!(pos, vec![10, 20], "no duplicate positions");
        assert_eq!(val, vec![4.0, 6.0]);
    }

    #[test]
    fn a_finer_first_array_repeats_the_second_head_value() {
        // a ends at 10, 15, 20; b only at 20. Between b's endpoints the output
        // carries b's head value forward, which is the merge's defining property.
        let a = t(ChromId(0), 0, &[10, 15, 20], &[1.0, 2.0, 3.0]);
        let b = t(ChromId(0), 0, &[20], &[9.0]);
        let (pos, val) = over_two_pv_array(&a, &b, Reducer::Max);
        assert_eq!(pos, vec![10, 15, 20]);
        assert_eq!(val, vec![9.0, 9.0, 9.0]);
    }

    #[test]
    fn the_shorter_array_dominates_and_the_tail_is_dropped() {
        // b is exhausted after 10, so a's remaining entries are dropped. This is
        // upstream's `while i1 < l1 and i2 < l2` with no tail drain.
        let a = t(ChromId(0), 0, &[10, 20, 30], &[1.0, 2.0, 3.0]);
        let b = t(ChromId(0), 0, &[10], &[5.0]);
        let (pos, val) = over_two_pv_array(&a, &b, Reducer::Max);
        assert_eq!(pos, vec![10], "a's tail at 20 and 30 is dropped");
        assert_eq!(val, vec![5.0]);
    }

    #[test]
    fn an_empty_second_array_yields_nothing() {
        let a = t(ChromId(0), 0, &[10, 20], &[1.0, 2.0]);
        let b = SignalTrack::empty(ChromId(0), 0, 20);
        let (pos, val) = over_two_pv_array(&a, &b, Reducer::Max);
        assert!(pos.is_empty() && val.is_empty());
    }

    #[test]
    fn min_and_mean_reducers() {
        let a = t(ChromId(0), 0, &[10], &[1.0]);
        let b = t(ChromId(0), 0, &[10], &[4.0]);
        assert_eq!(over_two_pv_array(&a, &b, Reducer::Min).1, vec![1.0]);
        assert_eq!(over_two_pv_array(&a, &b, Reducer::Mean).1, vec![2.5]);
        assert_eq!(over_two_pv_array(&a, &b, Reducer::Max).1, vec![4.0]);
    }

    #[test]
    fn folding_several_tracks_matches_pairwise_max() {
        let a = t(ChromId(0), 0, &[10, 20, 30], &[1.0, 2.0, 3.0]);
        let b = t(ChromId(0), 0, &[10, 20, 30], &[5.0, 1.0, 4.0]);
        let c = t(ChromId(0), 0, &[10, 20, 30], &[2.0, 9.0, 0.0]);
        let m = pointwise_max(&[a, b, c]);
        assert_eq!(m.len(), 3);
        let vals: Vec<f32> = m.runs().iter().map(|r| r.value).collect();
        assert_eq!(vals, vec![5.0, 9.0, 4.0]);
    }

    #[test]
    fn a_single_track_passes_through_the_fold() {
        let a = t(ChromId(0), 0, &[10, 20], &[7.0, 8.0]);
        let m = pointwise_max(&[a]);
        assert_eq!(m.len(), 2);
        assert_eq!(m.runs()[1].value, 8.0);
    }

    #[test]
    fn rebuild_from_pv_preserves_ends_and_values() {
        let pos = vec![10, 20, 30];
        let val = vec![1.0f32, 2.0, 3.0];
        let tr = track_from_pv(ChromId(0), 0, &pos, &val);
        assert_eq!(tr.start(), 0);
        assert_eq!(tr.end(), 30);
        assert_eq!(tr.runs()[0].end, 10);
        assert_eq!(tr.runs()[2].value, 3.0);
    }

    #[test]
    fn rebuilding_an_empty_merge_gives_an_empty_track() {
        let tr = track_from_pv(ChromId(0), 5, &[], &[]);
        assert!(tr.is_empty());
        assert_eq!(tr.start(), 5);
    }
}
