//! The `CallerFromAlignments` driver: per-chromosome treatment/control pairing.
//!
//! Transcribed from `MACS3.Signal.CallPeakUnit.pileup_treat_ctrl_a_chromosome`
//! and `__chrom_pair_treat_ctrl`. This is the stage that turns two read tracks
//! into the `(pos, treat, ctrl)` triple every downstream stage consumes.
//!
//! # Pairing is a cursor merge, and it drops a tail
//!
//! `__chrom_pair_treat_ctrl` walks the treatment and control `[p, v]` arrays
//! together and emits one entry per step, advancing whichever cursor has the
//! smaller end position (both on a tie). The emitted pair is the **current head**
//! of each array — again not a resample — and the loop is
//! `while it < lt and ic < lc`, so the remainder of the longer array is
//! discarded.
//!
//! Consequence worth stating plainly: the paired arrays are shorter than either
//! input, and they stop wherever the *shorter* input stops. Everything downstream
//! — p-scores, cutoffs, peak boundaries — is computed only over that prefix. A
//! port that drained the longer array would call peaks that upstream never sees.
//!
//! # The three lambda scales
//!
//! `ctrl_d_s` is built by `PeakDetect` in a fixed order, and the order is part of
//! the contract because the merge is a left fold:
//!
//! | index | scale |
//! |---|---|
//! | 0 | `d` (the fragment length) |
//! | 1 | `--slocal`, default 1000 |
//! | 2 | `--llocal`, default 10000 |
//!
//! Each gets its own bidirectional pileup with its own scaling factor, and the
//! three are combined with `max`. Because `over_two_pv_array` truncates, and
//! larger windows produce *more* breakpoints, the widest scale is usually the one
//! that decides where coverage ends.
//!
//! # `--nolambda` replaces control with a single baseline point
//!
//! Upstream does not compute a control pileup at all. It substitutes
//!
//! ```text
//! ctrl_pv = [treat_pv[0][-1:], np.array([lambda_bg])]
//! ```
//!
//! — a one-entry array holding just the global lambda, positioned at the
//! treatment's *last* breakpoint.
//!
//! Pairing that against the treatment therefore does **not** yield one entry. The
//! merge advances on the smaller of the two ends, so it walks the treatment to its
//! final breakpoint carrying `lambda_bg` as the control value the whole way, and
//! stops when both cursors land on that last position together. For a treatment of
//! `n` runs the result is `n` entries, every control value `lambda_bg`.
//!
//! Verified against MACS3 3.0.5: `callpeak --nolambda --bdg` on the
//! `se_model/realistic` fixture writes 376 pileup lines spanning
//! `chr20 0..93003`, i.e. the full treatment extent, not a single point. So the
//! observable effect of `--nolambda` is that the *control* is a constant, and the
//! treatment is untouched — not that the arrays collapse.

use crate::merge::pointwise_max;

use macs_core::Coord;
use macs_rle::SignalTrack;

/// One scale of the local lambda, i.e. one entry of `ctrl_d_s`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LambdaScale {
    /// Extension size in bp. `d` for index 0, `--slocal` for 1, `--llocal` for 2.
    pub d: i64,
    /// Multiplier applied to that scale's control pileup.
    pub scale_factor: f32,
}

/// The three local-lambda scales, in upstream's order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LambdaScales {
    /// Fragment length `d`.
    pub d: i64,
    /// `--slocal`, the small local window.
    pub slocal: i64,
    /// `--llocal`, the large local window.
    pub llocal: i64,
    /// Multiplier for the `d` scale.
    pub d_factor: f32,
    /// Multiplier for the `slocal` scale.
    pub slocal_factor: f32,
    /// Multiplier for the `llocal` scale.
    pub llocal_factor: f32,
}

impl LambdaScales {
    /// The scales as upstream's `ctrl_d_s` / `ctrl_scaling_factor_s` pair.
    ///
    /// Upstream appends `d` first, then `sregion`, then `lregion`, and only when
    /// `--to-small` scaling is in play. With `nocontrol` the lists are empty and
    /// the global lambda is used instead.
    pub fn as_pairs(&self) -> Vec<LambdaScale> {
        let mut out = vec![LambdaScale {
            d: self.d,
            scale_factor: self.d_factor,
        }];
        // `slocal` and `llocal` are conditional upstream: `PeakDetect` appends the
        // small window only when `sregion` is set, and the large window only when
        // `lregion > sregion`. A zero here means "not configured", which is not the
        // same as a window of width 0.
        if self.slocal > 0 {
            out.push(LambdaScale {
                d: self.slocal,
                scale_factor: self.slocal_factor,
            });
            if self.llocal > self.slocal {
                out.push(LambdaScale {
                    d: self.llocal,
                    scale_factor: self.llocal_factor,
                });
            }
        }
        out
    }
}

/// The paired per-chromosome arrays, in the end-indexed `[p, v]` form.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PairedSignal {
    /// Run ends; entry `i` covers `[pos[i-1], pos[i])` (F5).
    pub pos: Vec<Coord>,
    /// Treatment pileup, possibly scaled.
    pub treat: Vec<f32>,
    /// Control lambda, the pointwise max over the scales.
    pub ctrl: Vec<f32>,
}

impl PairedSignal {
    /// Number of runs.
    pub fn len(&self) -> usize {
        self.pos.len()
    }

    /// `true` when no runs were produced.
    pub fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }

    /// The treatment and control values at run `i`, if present.
    pub fn at(&self, i: usize) -> Option<(f32, f32)> {
        match (self.treat.get(i), self.ctrl.get(i)) {
            (Some(t), Some(c)) => Some((*t, *c)),
            _ => None,
        }
    }

    /// Position at the *start* of run `i`, following F5.
    ///
    /// This is what `__chrom_call_peak_using_certain_criteria` reads as
    /// `above_cutoff_startpos`; index 0 starts at 0 by definition.
    pub fn start_of(&self, i: usize) -> Coord {
        if i == 0 {
            0
        } else {
            self.pos[i - 1]
        }
    }
}

/// Pair a treatment and a control `[p, v]` pair, dropping the longer tail.
///
/// Transcribed from `__chrom_pair_treat_ctrl`.
pub fn pair_treat_ctrl(
    treat_pos: &[Coord],
    treat_val: &[f32],
    ctrl_pos: &[Coord],
    ctrl_val: &[f32],
) -> PairedSignal {
    let lt = treat_pos.len().min(treat_val.len());
    let lc = ctrl_pos.len().min(ctrl_val.len());
    let mut out = PairedSignal::default();
    out.pos.reserve(lt + lc);

    let (mut it, mut ic) = (0usize, 0usize);
    while it < lt && ic < lc {
        // head-to-head: each array's own current value, not a resample
        out.treat.push(treat_val[it]);
        out.ctrl.push(ctrl_val[ic]);
        let p1 = treat_pos[it];
        let p2 = ctrl_pos[ic];
        match p1.cmp(&p2) {
            std::cmp::Ordering::Less => {
                out.pos.push(p1);
                it += 1;
            }
            std::cmp::Ordering::Greater => {
                out.pos.push(p2);
                ic += 1;
            }
            std::cmp::Ordering::Equal => {
                out.pos.push(p1);
                it += 1;
                ic += 1;
            }
        }
    }
    out
}

/// [`pair_treat_ctrl`] over two run-length tracks.
///
/// The walk needs the raw `(end, value)` lists rather than the merged arrays, so
/// the caller's chunk boundaries and its control values come from the *same*
/// pairing -- see the F150 note at the single-end call site.
pub fn pair_treat_ctrl_runs(
    treat: &SignalTrack<f32>,
    ctrl: &SignalTrack<f32>,
) -> (Vec<Coord>, Vec<f32>, Vec<f32>) {
    let t_runs = treat.runs();
    let c_runs = ctrl.runs();
    let mut out = PairedSignal::default();
    out.pos.reserve(t_runs.len() + c_runs.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < t_runs.len() && j < c_runs.len() {
        out.treat.push(t_runs[i].value);
        out.ctrl.push(c_runs[j].value);
        let p1 = t_runs[i].end;
        let p2 = c_runs[j].end;
        match p1.cmp(&p2) {
            std::cmp::Ordering::Less => {
                out.pos.push(p1);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.pos.push(p2);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                out.pos.push(p1);
                i += 1;
                j += 1;
            }
        }
    }
    (out.pos, out.treat, out.ctrl)
}

/// Combine several control pileups into one local lambda.
///
/// Folds left with `max`, matching upstream's loop over `ctrl_d_s`.
pub fn combine_control_scales(tracks: &[SignalTrack<f32>]) -> SignalTrack<f32> {
    if tracks.is_empty() {
        panic!("combine_control_scales needs at least one track");
    }
    if tracks.len() == 1 {
        return tracks[0].clone();
    }
    pointwise_max(tracks)
}

/// The `--nolambda` control substitute: one entry at the treatment's last end.
///
/// Upstream builds `[treat_pv[0][-1:], [lambda_bg]]`, so merging it against the
/// treatment yields a single entry. Reproduced literally.
pub fn nolambda_control(treat_pos: &[Coord], lambda_bg: f32) -> (Vec<Coord>, Vec<f32>) {
    match treat_pos.last() {
        Some(last) => (vec![*last], vec![lambda_bg]),
        None => (Vec::new(), Vec::new()),
    }
}

/// Merge a treatment `[p, v]` with a control `[p, v]`, honouring `--nolambda`.
///
/// This is the `pileup_treat_ctrl_a_chromosome` path minus the pileup computation:
/// the caller supplies the already-built arrays.
pub fn pair_for_chromosome(
    treat_pos: &[Coord],
    treat_val: &[f32],
    ctrl_pos: &[Coord],
    ctrl_val: &[f32],
) -> PairedSignal {
    pair_treat_ctrl(treat_pos, treat_val, ctrl_pos, ctrl_val)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pairing_equal_length_arrays_aligns_one_to_one() {
        let p = pair_treat_ctrl(&[10, 20], &[1.0, 2.0], &[10, 20], &[5.0, 6.0]);
        assert_eq!(p.pos, vec![10, 20]);
        assert_eq!(p.treat, vec![1.0, 2.0]);
        assert_eq!(p.ctrl, vec![5.0, 6.0]);
    }

    #[test]
    fn pairing_repeats_the_other_head_while_a_cursor_waits() {
        // treatment at 10,15,20 and control only at 20: control's head value is
        // carried forward, and the treatment's own value advances underneath it
        let p = pair_treat_ctrl(&[10, 15, 20], &[1.0, 2.0, 3.0], &[20], &[9.0]);
        assert_eq!(p.pos, vec![10, 15, 20]);
        assert_eq!(p.treat, vec![1.0, 2.0, 3.0]);
        assert_eq!(
            p.ctrl,
            vec![9.0, 9.0, 9.0],
            "control's head is carried forward"
        );
    }

    #[test]
    fn pairing_drops_the_longer_tail() {
        // upstream's `while it < lt and ic < lc` with no drain
        let p = pair_treat_ctrl(&[10, 20, 30], &[1.0, 2.0, 3.0], &[10], &[5.0]);
        assert_eq!(
            p.pos,
            vec![10],
            "entries past the shorter array are dropped"
        );
        assert_eq!(p.len(), 1);
    }

    #[test]
    fn pairing_an_empty_control_yields_nothing() {
        let p = pair_treat_ctrl(&[10, 20], &[1.0, 2.0], &[], &[]);
        assert!(p.is_empty());
    }

    #[test]
    fn run_starts_follow_the_end_indexed_convention() {
        let p = pair_treat_ctrl(
            &[10, 20, 30],
            &[1.0, 2.0, 3.0],
            &[10, 20, 30],
            &[1.0, 2.0, 3.0],
        );
        assert_eq!(p.start_of(0), 0, "the first run starts at 0 by definition");
        assert_eq!(p.start_of(1), 10);
        assert_eq!(p.start_of(2), 20);
        assert_eq!(p.pos[2], 30, "and each run ends at its own position");
    }

    #[test]
    fn nolambda_control_is_a_single_baseline_entry() {
        let (pos, val) = nolambda_control(&[10, 20, 30], 0.75);
        assert_eq!(pos, vec![30], "only the last treatment breakpoint");
        assert_eq!(val, vec![0.75]);
    }

    #[test]
    fn nolambda_then_pairing_covers_the_whole_treatment_with_a_constant_control() {
        // The one-entry control sits at the treatment's *last* breakpoint, so the
        // merge walks the treatment to that point carrying lambda_bg throughout.
        // It does not collapse to a single entry.
        let treat_pos = [10u32, 20, 30];
        let treat_val = [1.0f32, 5.0, 2.0];
        let (cpos, cval) = nolambda_control(&treat_pos, 0.75);
        assert_eq!(cpos, vec![30], "control is anchored at the last breakpoint");
        let p = pair_for_chromosome(&treat_pos, &treat_val, &cpos, &cval);
        assert_eq!(p.pos, vec![10, 20, 30], "the treatment extent is preserved");
        assert_eq!(
            p.treat,
            vec![1.0, 5.0, 2.0],
            "treatment values are untouched"
        );
        assert_eq!(
            p.ctrl,
            vec![0.75, 0.75, 0.75],
            "control is constant lambda_bg"
        );
    }

    #[test]
    fn nolambda_stops_where_both_cursors_meet() {
        // treatment longer than the control's anchor: the merge ends at the anchor
        let p = pair_treat_ctrl(&[10, 20, 30, 40], &[1.0, 2.0, 3.0, 4.0], &[20], &[0.5]);
        assert_eq!(p.pos, vec![10, 20], "stops at the shared endpoint 20");
        assert_eq!(p.ctrl, vec![0.5, 0.5]);
    }

    #[test]
    fn nolambda_on_an_empty_treatment_yields_nothing() {
        let (pos, val) = nolambda_control(&[], 1.0);
        assert!(pos.is_empty() && val.is_empty());
    }

    #[test]
    fn at_returns_the_pair_or_none() {
        let p = pair_treat_ctrl(&[10, 20], &[1.0, 2.0], &[10, 20], &[3.0, 4.0]);
        assert_eq!(p.at(0), Some((1.0, 3.0)));
        assert_eq!(p.at(1), Some((2.0, 4.0)));
        assert_eq!(p.at(2), None);
    }

    #[test]
    fn the_scale_order_is_d_then_slocal_then_llocal() {
        let s = LambdaScales {
            d: 200,
            slocal: 1000,
            llocal: 10_000,
            d_factor: 1.0,
            slocal_factor: 0.5,
            llocal_factor: 0.25,
        };
        let pairs = s.as_pairs();
        assert_eq!(pairs.len(), 3, "all three configured");
        assert_eq!(pairs[0].d, 200, "index 0 is d");
        assert_eq!(pairs[1].d, 1000, "index 1 is slocal");
        assert_eq!(pairs[2].d, 10_000, "index 2 is llocal");
        assert_eq!(pairs[0].scale_factor, 1.0);
        assert_eq!(pairs[2].scale_factor, 0.25);
    }

    #[test]
    fn a_disabled_window_drops_its_scale() {
        // upstream omits slocal when it is unset and llocal unless it exceeds slocal
        let s = LambdaScales {
            d: 200,
            slocal: 0,
            llocal: 10_000,
            d_factor: 1.0,
            slocal_factor: 0.0,
            llocal_factor: 0.0,
        };
        assert_eq!(s.as_pairs().len(), 1, "only the d scale remains");

        let s = LambdaScales {
            d: 200,
            slocal: 1000,
            llocal: 500,
            ..s
        };
        assert_eq!(
            s.as_pairs().len(),
            2,
            "llocal <= slocal is dropped, not clamped"
        );
    }

    #[test]
    fn combining_scales_takes_the_pointwise_max() {
        use macs_core::ChromId;
        use macs_rle::Run;
        let mk = |vals: &[f32]| {
            let ends: Vec<Coord> = (1..=vals.len() as u32 * 10).step_by(10).collect();
            let runs: Vec<Run<f32>> = ends
                .iter()
                .zip(vals)
                .map(|(e, v)| Run::new(*e, *v))
                .collect();
            SignalTrack::from_runs(ChromId(0), 0, *ends.last().unwrap(), runs)
        };
        let a = mk(&[1.0, 2.0, 3.0]);
        let b = mk(&[5.0, 1.0, 4.0]);
        let c = mk(&[2.0, 9.0, 0.0]);
        let m = combine_control_scales(&[a, b, c]);
        let vals: Vec<f32> = m.runs().iter().map(|r| r.value).collect();
        assert_eq!(vals, vec![5.0, 9.0, 4.0]);
    }

    #[test]
    fn a_single_scale_passes_through() {
        use macs_core::ChromId;
        use macs_rle::Run;
        let runs = vec![Run::new(10u32, 7.0f32), Run::new(20, 8.0)];
        let a = SignalTrack::from_runs(ChromId(0), 0, 20, runs);
        let m = combine_control_scales(&[a]);
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn the_merge_helper_round_trips_through_track_from_pv() {
        use crate::merge::{over_two_pv_array, track_from_pv, Reducer};
        use macs_core::ChromId;
        let (pos, val) = over_two_pv_array(
            &track_from_pv(ChromId(0), 0, &[10, 20], &[1.0, 2.0]),
            &track_from_pv(ChromId(0), 0, &[10, 20], &[3.0, 4.0]),
            Reducer::Max,
        );
        assert_eq!(pos, vec![10, 20]);
        assert_eq!(val, vec![3.0, 4.0]);
    }
}
