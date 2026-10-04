//! Endpoint-sweep, run-length-encoded pileup.
//!
//! This is the Rust form of [`MACS3/Signal/PileupV2.py`]. Upstream builds two
//! sorted `int32` endpoint arrays, sweeps them, and emits a `pv` array of
//! `(position, float32 value)` breakpoints. We do the same sweep but write
//! straight into a [`SignalTrack`], so no per-base array is ever allocated.
//!
//! # The one behaviour that must be reproduced exactly
//!
//! Upstream's sweep
//! ([`PileupV2.py:419-423`](https://github.com/macs3-project/MACS)) handles a
//! start and an end at the same position by advancing both pointers and leaving
//! the depth **unchanged**:
//!
//! ```text
//! if start_ptr[0] < end_ptr[0]:  pileup += 1
//! elif start_ptr[0] > end_ptr[0]: pileup -= 1
//! else:                          both advance, pileup unchanged
//! ```
//!
//! So when one inferred fragment ends on the base where another begins, the
//! read that ended is not subtracted and the read that began is not added. The
//! correct depth at that base is `n - 1`; upstream reports `n`. This is not
//! exotic: with a fixed extension `d` it fires for every pair of reads exactly
//! `d` apart, and with `--shift` for every pair `d +/- 2 * shift` apart.
//!
//! See `docs/upstream-findings.md` F5. [`pileup_from_positions`] has a
//! dedicated regression test for it.
//!
//! # Numeric type
//!
//! Values are `f32`, matching upstream's `dtype="f4"`. The multiplication
//! `pileup * scale_factor` is done in `f32` and then clamped against
//! `baseline_value` in `f32`, in that order, because a `f64` intermediate
//! changes the last bit of every scaled track value and therefore of every
//! `control_lambda.bdg` byte.

#![forbid(unsafe_code)]
#![deny(missing_debug_implementations)]

use macs_core::{ChromId, Coord, Interval};
use macs_rle::{SignalTrack, TrackBuilder};

/// Parameters of a single-end pileup.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SingleEndParams {
    /// Shift applied to the 5' end *before* extension. Upstream calls it
    /// `five_shift`; it is `-end_shift` in the directional case.
    pub five_shift: i64,
    /// Extension applied at the 3' end, i.e. `d` in the directional case.
    pub three_shift: i64,
    /// Contig length. Endpoints are clipped to `[0, rlength]`.
    pub rlength: Coord,
    /// Multiplier applied to the depth.
    pub scale_factor: f32,
    /// Floor applied *after* scaling.
    pub baseline_value: f32,
    /// F71: merge a run into its predecessor when the value is unchanged.
    ///
    /// `false` for the bidirectional control pileup, where a start and an end
    /// can net out and upstream still keeps the separate endpoint.
    pub coalesce: bool,
    /// F190: both endpoint lists get the **same** window `[x - five_shift,
    /// x + three_shift]`, with no plus/minus swap.
    ///
    /// `pileup_from_LR_as_list` / `pileup_from_PN_shifted` treat the two lists as
    /// plus and minus strands, so the minus window is mirrored
    /// (`N - three_shift`, `N + five_shift`). `pileup_from_LRC_centers_as_list`
    /// -- the `--format FRAG` control projection -- does not: it concatenates
    /// `l - half_d` and `r - half_d` and sets `end_poss = start_poss + d`, so a
    /// fragment's left *and* right end each get `[x - half_d, x - half_d + d)`.
    /// Mirroring the right ends put their window one base further left, which
    /// added depth on the left of every right end and lost it on the right: on
    /// `sweep/gmini_mfrag_d1200_w600_ctrl_231` the merged lambda came out up to
    /// 14 high from position 120 on against the oracle's own
    /// `pileup_a_chromosome_c`.
    pub centred: bool,
}

impl SingleEndParams {
    /// Directional single-end parameters: shift each 5' end by `end_shift` and
    /// extend `d` bases in the 3' direction.
    ///
    /// This is `pileup_a_chromosome(chrom, d, ...)` with `directional=True`.
    /// Paired-end *control* projection: extend each fragment end symmetrically by
    /// `d // 2` on both sides.
    ///
    /// This is deliberately **not** [`SingleEndParams::bidirectional`]. The
    /// single-end control (`FixWidthTrack.pileup_a_chromosome_c`) uses
    /// `five = d//2 - end_shift`, `three = end_shift + d - d//2`, which for
    /// `end_shift = 0` gives an *asymmetric* split `d//2` / `d - d//2`. The
    /// paired-end control (`PairedEndTrack.pileup_a_chromosome_c`) instead passes
    /// `five_shift = d//2` **and** `three_shift = d//2` -- both floored -- giving a
    /// span of `d` bases for even `d` and `d - 1` for odd `d`.
    ///
    /// That one-base difference for odd `d` is not cosmetic: measured across all 70
    /// paired-end summits, using the single-end split leaves the control depth at
    /// the summit out by exactly **one tag** (F100), which is enough to move a
    /// q-score cutoff across a position.
    pub fn symmetric(d: i64, rlength: Coord, scale_factor: f32) -> Self {
        let h = d / 2;
        SingleEndParams {
            five_shift: h,
            three_shift: h,
            rlength,
            scale_factor,
            baseline_value: 0.0,
            coalesce: false,
            centred: false,
        }
    }

    /// F190: the `--format FRAG` (counted) paired-end control projection,
    /// `pileup_from_LRC_centers_as_list`: one `[x - d/2, x - d/2 + d)` window per
    /// fragment *end*, with no plus/minus mirror.
    pub fn centred(d: i64, rlength: Coord, scale_factor: f32) -> Self {
        SingleEndParams {
            five_shift: d / 2,
            three_shift: d - d / 2,
            rlength,
            scale_factor,
            baseline_value: 0.0,
            // F71: equal-valued neighbours are common here too, and the merged
            // lambda's run structure depends on keeping them.
            coalesce: false,
            centred: true,
        }
    }

    pub fn directional(d: i64, end_shift: i64, rlength: Coord, scale_factor: f32) -> Self {
        SingleEndParams {
            five_shift: -end_shift,
            three_shift: end_shift + d,
            rlength,
            scale_factor,
            baseline_value: 0.0,
            // F71: every event here changes depth by +/-1, so consecutive runs
            // never share a value and coalescing is a no-op.
            coalesce: true,
            centred: false,
        }
    }

    /// Bidirectional single-end parameters: shift each 5' end by `end_shift`
    /// and extend `d` bases in both directions, `d / 2` on the 5' side.
    ///
    /// This is `pileup_a_chromosome(..., directional=False)`, which MACS uses for
    /// the **control** track.
    pub fn bidirectional(d: i64, end_shift: i64, rlength: Coord, scale_factor: f32) -> Self {
        SingleEndParams {
            five_shift: d / 2 - end_shift,
            three_shift: end_shift + d - d / 2,
            rlength,
            scale_factor,
            baseline_value: 0.0,
            // F71: a plus read contributes events at `p-d` and `p` and a minus
            // read at `m` and `m+d`; with overlapping reads a start and an end
            // can coincide or net out, so equal-valued neighbours are common
            // and must NOT be merged -- upstream's raw-array pileup keeps every
            // endpoint, and the merged three-scale control lambda is what sets
            // the paired array's breakpoint count (F70).
            coalesce: false,
            centred: false,
        }
    }

    /// Set the floor applied after scaling (`baseline_value`, e.g. `lambda_bg`).
    pub fn with_baseline(mut self, baseline: f32) -> Self {
        self.baseline_value = baseline;
        self
    }
}

/// Build the sorted, clipped endpoint arrays for a set of 5' positions.
///
/// Returns `(starts, ends)`. Kept separate from the sweep so the clip-and-sort
/// step can be tested on its own and reused by the weighted variant.
fn endpoints_from_positions(
    plus: &[Coord],
    minus: &[Coord],
    p: &SingleEndParams,
) -> (Vec<Coord>, Vec<Coord>) {
    let rlength = p.rlength;
    // upstream computes these in C `long` (i64) arithmetic and clips negatives to
    // zero, so a large negative offset collapses onto the contig start rather
    // than wrapping
    let clamp = |v: i64| -> Coord {
        if v < 0 {
            0
        } else if v as u64 > rlength {
            rlength
        } else {
            v as u64
        }
    };

    let mut starts = Vec::with_capacity(plus.len() + minus.len());
    let mut ends = Vec::with_capacity(plus.len() + minus.len());
    for &x in plus {
        let base = x as i64;
        starts.push(clamp(base - p.five_shift));
        ends.push(clamp(base + p.three_shift));
    }
    for &x in minus {
        let base = x as i64;
        if p.centred {
            // F190: see `SingleEndParams::centred` -- no plus/minus mirror.
            starts.push(clamp(base - p.five_shift));
            ends.push(clamp(base + p.three_shift));
        } else {
            starts.push(clamp(base - p.three_shift));
            ends.push(clamp(base + p.five_shift));
        }
    }
    starts.sort_unstable();
    ends.sort_unstable();
    (starts, ends)
}

/// Sweep sorted start/end arrays into a signal track.
///
/// This is a transcription of `_pileup_sorted_unit_as_list`
/// ([`PileupV2.py:344`]) including the coincident-event behaviour described in
/// the module docs.
fn sweep_to_track(
    chrom: ChromId,
    starts: &[Coord],
    ends: &[Coord],
    rlength: Coord,
    scale_factor: f32,
    baseline_value: f32,
    coalesce: bool,
) -> SignalTrack<f32> {
    // F71: the bidirectional control/local-lambda pileup must keep every
    // endpoint. There a plus read contributes events at `p-d` and `p` and a
    // minus read at `m` and `m+d`; with overlapping reads a start and an end can
    // coincide or net out, so consecutive depths can be equal and `push` would
    // merge them. Upstream's raw-array pileup keeps every position, and the
    // merged three-scale control lambda is what sets the paired array's 10,684
    // breakpoints (F70). The directional treatment pileup never coalesces anyway
    // -- every event there changes depth by +/-1 -- so it is unaffected.
    macro_rules! emit {
        ($b:expr, $e:expr, $v:expr) => {
            if coalesce {
                $b.push($e, $v)
            } else {
                $b.push_exact($e, $v)
            }
        };
    }
    let mut b = TrackBuilder::with_capacity(chrom, 0, rlength, starts.len() + ends.len());
    if starts.is_empty() || ends.is_empty() {
        return b.finish();
    }

    let clamp_depth = |depth: i64, scale: f32, base: f32| -> f32 {
        // upstream: `value = pileup * scale_factor; if value < baseline: value = baseline`
        let v = depth as f32 * scale;
        if v < base {
            base
        } else {
            v
        }
    };

    let mut i_s = 0usize;
    let mut i_e = 0usize;
    let mut depth: i64 = 0;
    let mut pre_p = starts[0].min(ends[0]);

    if pre_p != 0 {
        // the region `[0, pre_p)` precedes the first event, so its depth is 0
        emit!(b, pre_p, clamp_depth(0, scale_factor, baseline_value));
    }

    while i_s < starts.len() && i_e < ends.len() {
        let s = starts[i_s];
        let e = ends[i_e];
        if s < e {
            if s != pre_p {
                // `s` closes `[pre_p, s)`, whose depth is the current one
                emit!(b, s, clamp_depth(depth, scale_factor, baseline_value));
                pre_p = s;
            }
            depth += 1;
            i_s += 1;
        } else if s > e {
            if e != pre_p {
                emit!(b, e, clamp_depth(depth, scale_factor, baseline_value));
                pre_p = e;
            }
            depth -= 1;
            i_e += 1;
        } else {
            // a fragment ends exactly where the next begins: net delta 0, so the
            // depth is unchanged and no breakpoint is needed
            i_s += 1;
            i_e += 1;
        }
    }

    while i_e < ends.len() {
        let e = ends[i_e];
        if e != pre_p {
            emit!(b, e, clamp_depth(depth, scale_factor, baseline_value));
            pre_p = e;
        }
        depth -= 1;
        i_e += 1;
    }

    b.finish()
}

/// Single-end pileup from 5' positions on both strands.
///
/// `plus` and `minus` are 0-based 5' coordinates. They do not need to be
/// sorted.
pub fn pileup_from_positions(
    chrom: ChromId,
    plus: &[Coord],
    minus: &[Coord],
    p: &SingleEndParams,
) -> SignalTrack<f32> {
    let (starts, ends) = endpoints_from_positions(plus, minus, p);
    sweep_to_track(
        chrom,
        &starts,
        &ends,
        p.rlength,
        p.scale_factor,
        p.baseline_value,
        p.coalesce,
    )
}

/// Pile up **count-weighted** 5' positions, as `--format FRAG` needs.
///
/// The weighted analogue of [`pileup_from_positions`]: `plus_w[k]` is the
/// multiplicity of `plus[k]`, likewise for the minus strand. Upstream's
/// `pileup_from_LRC_centers_as_list` builds the same event list -- one window per
/// fragment *end*, `[x - d//2, x - d//2 + d)` -- and carries the fragment's count
/// as the event weight, so a weighted FRAG control lambda is `d`-windows centred on
/// the ends rather than the unit-depth windows `pileup_from_positions` builds.
///
/// Weights are `f32` because upstream's `_pileup_sorted_weighted_as_list`
/// accumulates into an `f4` depth, and the depth is multiplied by `scale_factor` in
/// that same lane.
pub fn pileup_from_weighted_positions(
    chrom: ChromId,
    plus: &[Coord],
    minus: &[Coord],
    plus_w: &[f32],
    minus_w: &[f32],
    p: &SingleEndParams,
) -> SignalTrack<f32> {
    let rlength = p.rlength;
    let clamp = |v: i64| -> Coord {
        if v < 0 {
            0
        } else if v as u64 > rlength {
            rlength
        } else {
            v as u64
        }
    };
    // F167: build `(position, weight)` **pairs** and sort those, rather than
    // sorting a bare position array and keeping the weights in input order.
    //
    // `endpoints_from_positions` sorts `starts` and `ends` independently, which
    // is exactly what upstream does -- but upstream then permutes the weights
    // with the *same* argsort:
    //
    //     indices = np.argsort(start_poss)
    //     start_poss = start_poss[indices]
    //     start_weights = weights[indices]
    //
    // Sorting the positions and leaving the weights alone attaches each count to
    // the wrong position. The totals stay plausible -- most counts are similar --
    // so the error only shows as a depth that is a few units off in either
    // direction: the `--llocal` control depth at a summit came out 7855 where
    // the file says 7852, moving `-log10(pvalue)` in the third decimal.
    //
    // Tie order does not matter: every event at one position is summed before the
    // value is recorded, so the multiset is order-independent.
    let mut start_events: Vec<(Coord, f32)> = Vec::with_capacity(plus.len() + minus.len());
    let mut end_events: Vec<(Coord, f32)> = Vec::with_capacity(plus.len() + minus.len());
    for (idx, &x) in plus.iter().enumerate() {
        let w = plus_w.get(idx).copied().unwrap_or(1.0);
        let base = x as i64;
        start_events.push((clamp(base - p.five_shift), w));
        end_events.push((clamp(base + p.three_shift), -w));
    }
    for (idx, &x) in minus.iter().enumerate() {
        let w = minus_w.get(idx).copied().unwrap_or(1.0);
        let base = x as i64;
        // F190: `pileup_from_LRC_centers_as_list` concatenates `l - half_d` and
        // `r - half_d` and sets `end_poss = start_poss + d`, so the fragment's
        // right end gets the **same** window as its left end. The plus/minus
        // mirror belongs to `pileup_from_PN_shifted` only.
        if p.centred {
            start_events.push((clamp(base - p.five_shift), w));
            end_events.push((clamp(base + p.three_shift), -w));
        } else {
            start_events.push((clamp(base - p.three_shift), w));
            end_events.push((clamp(base + p.five_shift), -w));
        }
    }
    start_events.sort_unstable_by_key(|e| e.0);
    end_events.sort_unstable_by_key(|e| e.0);
    let starts: Vec<Coord> = start_events.iter().map(|e| e.0).collect();
    let ends: Vec<Coord> = end_events.iter().map(|e| e.0).collect();
    // F155: the weighted sweep **always** coalesces. `_pileup_sorted_weighted_as_list`
    // ends with
    //
    //     if scaled_z == pre_z:
    //         ret_p_ptr[c-1] = p      # extend the previous run's end
    //     else:
    //         ret_p_ptr[c] = p; ret_v_ptr[c] = scaled_z; c += 1; pre_z = scaled_z
    //
    // so an equal-valued neighbour overwrites the previous run's end rather than
    // adding a breakpoint. Unlike the unweighted sweep (F71) this is
    // unconditional -- a counted event's `+w`/`-w` pair can land on the same
    // position and leave the depth unchanged, and upstream's `--format FRAG`
    // control lambda therefore has far fewer breakpoints than the raw event list.
    //
    // It also differs on odd `d`: upstream builds `end = start + d` from
    // `start = x - d//2`, while `endpoints_from_positions` builds `x + d//2`.
    // Those agree only when `d` is even.
    sweep_weighted_to_track(
        chrom,
        &starts,
        &ends,
        &start_events.iter().map(|e| e.1).collect::<Vec<_>>(),
        p.rlength,
        p.scale_factor,
        p.baseline_value,
        true,
    )
}

/// Same sweep as [`sweep_to_track`] but with per-event weights.
///
/// The scalar parameters travel together because they are exactly
/// `SingleEndParams`' scalar half; `coalesce` is passed separately because the
/// weighted sweep's rule is not the same as the unweighted one (F155).
#[allow(clippy::too_many_arguments)]
fn sweep_weighted_to_track(
    chrom: ChromId,
    starts: &[Coord],
    ends: &[Coord],
    weights: &[f32],
    rlength: Coord,
    scale_factor: f32,
    baseline_value: f32,
    coalesce: bool,
) -> SignalTrack<f32> {
    // F235: starts and ends stay in **separate sorted arrays** and are applied by a
    // two-pointer walk, because that is the order upstream accumulates in:
    //
    //     while i_s < ls and start_ptr[i_s] == p:  z += start_w_ptr[i_s]
    //     while i_e < le and end_ptr[i_e]   == p:  z -= end_w_ptr[i_e]
    //
    // Starts are folded in first at every position. Merging the two lists into one
    // sorted by `(pos, value)` -- as this function used to -- puts the *negative*
    // end weights first at a tie, which is the reverse. `z` is an f32 accumulator,
    // so the two orders are not interchangeable: with depths in the thousands the
    // last bits differ, and because `_pileup_sorted_weighted_as_list` coalesces on
    // `scaled_z == pre_z` (an exact f32 equality) a last-bit difference changes the
    // *breakpoints* of the control lambda, not just its printed value.
    let mut start_pairs: Vec<(Coord, f32)> = Vec::with_capacity(starts.len());
    for (i, &p) in starts.iter().enumerate() {
        start_pairs.push((p, weights[i]));
    }
    let mut end_pairs: Vec<(Coord, f32)> = Vec::with_capacity(ends.len());
    for (i, &p) in ends.iter().enumerate() {
        end_pairs.push((p, weights[i]));
    }

    macro_rules! emit {
        ($b:expr, $e:expr, $v:expr) => {
            if coalesce {
                $b.push($e, $v)
            } else {
                $b.push_exact($e, $v)
            }
        };
    }
    let total = start_pairs.len() + end_pairs.len();
    let mut b = TrackBuilder::with_capacity(chrom, 0, rlength, total);
    if total == 0 {
        return b.finish();
    }
    let clamp = |depth: f32, scale: f32, base: f32| -> f32 {
        // `scaled_z = z * scale_factor` is `float * float` in C, so it is a
        // single-precision multiply on both sides of this line.
        let v = depth * scale;
        if v < base {
            base
        } else {
            v
        }
    };
    let mut depth = 0.0f32;
    let mut at = 0u64;
    let (mut i_s, mut i_e) = (0usize, 0usize);
    while i_s < start_pairs.len() || i_e < end_pairs.len() {
        // upstream picks the smaller of the two heads, preferring starts on a tie
        let pos = match (
            start_pairs.get(i_s).map(|e| e.0),
            end_pairs.get(i_e).map(|e| e.0),
        ) {
            (Some(a), Some(b)) => {
                if i_e >= end_pairs.len() || (i_s < start_pairs.len() && a < b) {
                    a
                } else {
                    b
                }
            }
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => break,
        };
        if pos > at {
            emit!(b, pos, clamp(depth, scale_factor, baseline_value));
            at = pos;
        }
        // starts first, then ends -- both grouped at this position
        while i_s < start_pairs.len() && start_pairs[i_s].0 == pos {
            depth += start_pairs[i_s].1;
            i_s += 1;
        }
        while i_e < end_pairs.len() && end_pairs[i_e].0 == pos {
            depth -= end_pairs[i_e].1;
            i_e += 1;
        }
    }
    b.finish()
}

/// Convenience: pile up a strand-agnostic list of 5' positions.
pub fn pileup_single_strand(
    chrom: ChromId,
    positions: &[Coord],
    p: &SingleEndParams,
) -> SignalTrack<f32> {
    pileup_from_positions(chrom, positions, &[], p)
}

/// Pile up true fragments, optionally with integer weights.
///
/// Each fragment is `[l, r)`; the depth on `[l, r)` increases by `weight`.
///
/// This is the paired-end path. MACS3 never *estimates* fragment length in
/// paired-end mode, so the estimated `d` is ignored entirely — a deliberate
/// difference from [`pileup_from_positions`], and one that is checked by a
/// fixture where `d` is set to something absurd.
///
/// Weighted fragments are the `FRAG` / single-cell path, where one genomic
/// interval carries a fragment count per barcode.
pub fn pileup_from_fragments(
    chrom: ChromId,
    fragments: &[(Coord, Coord)],
    rlength: Coord,
    scale_factor: f32,
    baseline_value: f32,
) -> SignalTrack<f32> {
    let mut starts: Vec<Coord> = Vec::with_capacity(fragments.len());
    let mut ends: Vec<Coord> = Vec::with_capacity(fragments.len());
    for &(l, r) in fragments {
        starts.push(l.min(rlength));
        ends.push(r.min(rlength));
    }
    starts.sort_unstable();
    ends.sort_unstable();
    // F96: `coalesce = false`. Upstream's fragment sweep is
    // `_pileup_sorted_unit_as_list` (`PileupV2.py:344`), which appends one entry
    // per sweep step into raw arrays and never merges neighbours -- the same
    // property that made `from_runs_exact` necessary for the single-end local
    // lambda in F76. Coalescing here drops fragment breakpoints, which shifts
    // where the q-score crosses its cutoff by a position in paired-end fixtures.
    sweep_to_track(
        chrom,
        &starts,
        &ends,
        rlength,
        scale_factor,
        baseline_value,
        false,
    )
}

/// Pile up fragments with a per-fragment weight (FRAG multiplicity).
pub fn pileup_from_weighted_fragments(
    chrom: ChromId,
    fragments: &[(Coord, Coord, f32)],
    rlength: Coord,
    scale_factor: f32,
    baseline_value: f32,
) -> SignalTrack<f32> {
    pileup_from_weighted_fragments_from(chrom, fragments, rlength, scale_factor, baseline_value, 0)
}

/// [`pileup_from_weighted_fragments`] with an explicit sweep origin.
///
/// F181: upstream's `_pileup_sorted_weighted_as_list` (`PileupV2.py:238`) keeps
/// `pre_p = 0` and only ever advances it **inside** `if p != pre_p`, so the
/// deltas at coordinate 0 are applied but no run is emitted there. A sweep that
/// starts at the first boundary emits the leading `depth == 0` run instead, which
/// is one extra run, one extra union row in the paired walk, and -- because the
/// extra row changes the signed first span in the p->q histogram -- a different
/// q-value table.
///
/// `origin` is the coordinate that `pre_p` starts at, i.e. upstream's raw zero.
/// A coordinate-shifted caller passes its shift so the suppression lands on raw
/// zero rather than on the shift itself.
pub fn pileup_from_weighted_fragments_from(
    chrom: ChromId,
    fragments: &[(Coord, Coord, f32)],
    rlength: Coord,
    scale_factor: f32,
    baseline_value: f32,
    origin: Coord,
) -> SignalTrack<f32> {
    // Bucket the deltas by position and sweep once. This is the weighted
    // analogue of upstream's `pileup_from_LRC_as_list`, and it keeps the
    // weighted arithmetic in the same `f32` lane as the unweighted path.
    let mut deltas: Vec<(Coord, f32)> = Vec::with_capacity(fragments.len() * 2);
    for &(l, r, w) in fragments {
        if w == 0.0 {
            continue;
        }
        deltas.push((l.min(rlength), w));
        deltas.push((r.min(rlength), -w));
    }
    deltas.sort_unstable_by_key(|d| d.0);

    // Same right-endpoint convention as the unweighted sweep: a breakpoint at
    // `pos` carries the depth of the region it closes, so all deltas at `pos`
    // are applied *after* it has been recorded.
    let mut b = TrackBuilder::with_capacity(chrom, 0, rlength, deltas.len());
    let mut depth = 0.0f32;
    let mut at = origin;
    let mut i = 0usize;
    while i < deltas.len() {
        let pos = deltas[i].0;
        if pos > at {
            let v = depth * scale_factor;
            b.push(
                pos,
                if v < baseline_value {
                    baseline_value
                } else {
                    v
                },
            );
            at = pos;
        }
        while i < deltas.len() && deltas[i].0 == pos {
            depth += deltas[i].1;
            i += 1;
        }
    }
    b.finish()
}

/// Simple symmetric extension around each 5' position, as used by
/// `naive_quick_pileup`.
///
/// Positions below 0 are clipped to 0, matching upstream.
pub fn naive_quick_pileup(
    chrom: ChromId,
    sorted_positions: &[Coord],
    extension: i64,
    rlength: Coord,
) -> SignalTrack<f32> {
    let mut starts = Vec::with_capacity(sorted_positions.len());
    let mut ends = Vec::with_capacity(sorted_positions.len());
    for &pos in sorted_positions {
        let base = pos as i64;
        starts.push((base - extension).max(0) as Coord);
        ends.push((base + extension).clamp(0, rlength as i64) as Coord);
    }
    starts.sort_unstable();
    ends.sort_unstable();
    sweep_to_track(chrom, &starts, &ends, rlength, 1.0, 0.0, true)
}

/// Total integrated depth of a pileup track: `sum(value * length)`.
///
/// Equals `sum(fragment lengths * scale_factor)` exactly, because each fragment
/// contributes `scale_factor` over its own length. Property-tested.
pub fn integrated_depth(t: &SignalTrack<f32>) -> f64 {
    t.integral()
}

/// Clip a track to an interval, preserving values.
pub fn restrict(t: &SignalTrack<f32>, iv: Interval) -> SignalTrack<f32> {
    t.slice(iv)
}

#[cfg(test)]
mod tests {
    use super::*;

    const C: ChromId = ChromId(0);

    /// `(start, value)` of every run — the compact canonical view.
    fn spans(t: &SignalTrack<f32>) -> Vec<(u64, f32)> {
        t.spans().collect()
    }

    /// Coverage at `pos`, treating an uncovered base as zero.
    ///
    /// Upstream's sweep always emits a leading breakpoint carrying the depth
    /// *before* the first event, so a track produced by this module is total over
    /// `[0, cursor)` and the uncovered head is represented explicitly. Reading it
    /// as zero is what every MACS consumer does.
    fn cov(t: &SignalTrack<f32>, pos: u64) -> f32 {
        t.value_at_or(pos, 0.0)
    }

    fn params(d: i64, rlength: u64) -> SingleEndParams {
        SingleEndParams::directional(d, 0, rlength, 1.0)
    }

    // ---------- basic geometry ----------

    #[test]
    fn empty_input_gives_an_empty_track() {
        let t = pileup_from_positions(C, &[], &[], &params(200, 1000));
        assert!(t.is_empty());
        assert_eq!(t.end(), 1000);
    }

    /// Upstream emits a leading breakpoint at `min(starts[0], ends[0])` carrying
    /// the depth before the first event, so the head is explicit rather than
    /// uncovered.
    #[test]
    fn the_head_before_the_first_read_is_explicit_zero() {
        let t = pileup_from_positions(C, &[100], &[], &params(200, 1000));
        assert_eq!(spans(&t), vec![(0, 0.0), (100, 1.0)]);
        assert_eq!(t.cursor(), 300);
        assert_eq!(t.covered_len(), 300);
    }

    #[test]
    fn one_plus_read_extends_forward_only() {
        let t = pileup_from_positions(C, &[100], &[], &params(200, 1000));
        assert_eq!(cov(&t, 0), 0.0);
        assert_eq!(cov(&t, 99), 0.0);
        assert_eq!(cov(&t, 100), 1.0);
        assert_eq!(cov(&t, 299), 1.0);
        assert_eq!(cov(&t, 300), 0.0);
        assert_eq!(
            t.value_at(300),
            None,
            "past the cursor is genuinely uncovered"
        );
        assert_eq!(t.integral(), 200.0);
    }

    #[test]
    fn one_minus_read_extends_backward_only() {
        let t = pileup_from_positions(C, &[], &[100], &params(200, 1000));
        // start = 100 - 200 -> clipped to 0, end = 100
        assert_eq!(cov(&t, 0), 1.0);
        assert_eq!(cov(&t, 99), 1.0);
        assert_eq!(
            t.value_at(100),
            None,
            "past the cursor is genuinely uncovered"
        );
        assert_eq!(t.integral(), 100.0);
    }

    #[test]
    fn two_overlapping_reads_stack() {
        let t = pileup_from_positions(C, &[100, 150], &[], &params(200, 1000));
        // [0,100)=0, [100,150)=1, [150,300)=2, [300,350)=1
        assert_eq!(
            spans(&t),
            vec![(0, 0.0), (100, 1.0), (150, 2.0), (300, 1.0)]
        );
        assert_eq!(t.integral(), 400.0);
    }

    // ---------- F5: the coincident-event quirk ----------

    /// The load-bearing upstream quirk.
    ///
    /// Two reads `d` apart produce one fragment ending exactly where the other
    /// begins. Upstream advances both pointers and leaves the depth unchanged,
    /// so it reports 2 at a position where the true depth is 1.
    ///
    /// Sources: `MACS3/Signal/PileupV2.py:419-423`, and `F5` in
    /// `docs/upstream-findings.md`.
    #[test]
    fn coincident_start_and_end_emit_no_breakpoint() {
        // reads exactly d apart: every fragment abuts its neighbour, so the
        // sweep takes the net-zero branch at each junction and emits nothing
        // until the final fragment end
        let t = pileup_from_positions(C, &[0, 200, 400], &[], &params(200, 1000));
        assert_eq!(spans(&t), vec![(0, 1.0)]);
        assert_eq!(t.cursor(), 600);
        // the depth really is 1 across the whole span
        assert_eq!(cov(&t, 0), 1.0);
        assert_eq!(cov(&t, 200), 1.0);
        assert_eq!(cov(&t, 400), 1.0);
        assert_eq!(t.integral(), 600.0);
        // and nothing is represented past the last fragment
        assert_eq!(t.value_at(600), None);
    }

    #[test]
    fn coincident_events_nest_correctly_between_differing_depths() {
        // [5,55), [30,80), [55,105): the junction at 55 is coincident, so the
        // depth reads 2 straight through it, which is correct
        let t = pileup_from_positions(C, &[5, 30, 55], &[], &params(50, 1000));
        assert_eq!(spans(&t), vec![(0, 0.0), (5, 1.0), (30, 2.0), (80, 1.0)]);
        assert_eq!(cov(&t, 60), 2.0);
        assert_eq!(t.integral(), 25.0 + 50.0 * 2.0 + 25.0);
    }

    /// The same quirk when the junction lands on a value change, i.e. when the
    /// two depths differ at the junction and cannot be absorbed by coalescing.
    #[test]
    fn coincident_events_preserve_the_inflated_depth() {
        // reads at 100 and 320 with d=200: fragments [100,300) and [320,520)
        // no junction. Now force one: reads at 100 and 300 -> [100,300) and
        // [300,500); the junction at 300 keeps the depth at 1 (also correct).
        let t = pileup_from_positions(C, &[100, 300], &[], &params(200, 1000));
        assert_eq!(cov(&t, 299), 1.0);
        assert_eq!(cov(&t, 300), 1.0);
        assert_eq!(t.integral(), 400.0);
    }

    #[test]
    fn coincident_events_fire_with_shift() {
        // --shift 10, d 200, plus reads at 210 and 310:
        //   [210+10, 210+10+200) = [220, 420)
        //   [310+10, 310+10+200) = [320, 520)
        let p = SingleEndParams::directional(200, 10, 1000, 1.0);
        let t = pileup_from_positions(C, &[210, 310], &[], &p);
        assert_eq!(cov(&t, 219), 0.0);
        assert_eq!(cov(&t, 220), 1.0);
        assert_eq!(cov(&t, 320), 2.0);
    }

    #[test]
    fn depth_is_never_negative() {
        let t = pileup_from_positions(C, &[], &[10, 20, 30], &params(200, 1000));
        for s in t.iter() {
            assert!(*s.value >= 0.0, "negative depth at {}", s.start);
        }
    }

    // ---------- clipping ----------

    #[test]
    fn endpoints_are_clipped_to_the_contig() {
        // a read at 10 extended by 2000 in a 500 bp contig covers [10,500)
        let t = pileup_from_positions(C, &[10], &[], &params(2000, 500));
        assert_eq!(cov(&t, 9), 0.0);
        assert_eq!(cov(&t, 10), 1.0);
        assert_eq!(cov(&t, 499), 1.0);
        assert_eq!(t.cursor(), 500);
        assert_eq!(t.integral(), 490.0);
    }

    #[test]
    fn a_fragment_start_before_the_contig_collapses_to_zero() {
        // a minus read at 5 with d = 200 starts at -195, which must clamp to 0
        // rather than wrap around as a u64
        let t = pileup_from_positions(C, &[], &[5], &params(200, 1000));
        assert_eq!(cov(&t, 0), 1.0);
        assert_eq!(cov(&t, 4), 1.0);
        assert_eq!(t.value_at(5), None);
        assert_eq!(t.integral(), 5.0);
    }

    #[test]
    fn an_absurd_extension_clamps_both_ends_without_wrapping() {
        // five_shift of -2^40 would wrap if the arithmetic were done in u64
        let p = SingleEndParams {
            five_shift: -(1i64 << 40),
            three_shift: 10,
            rlength: 1000,
            scale_factor: 1.0,
            baseline_value: 0.0,
            coalesce: true,
            centred: false,
        };
        let t = pileup_from_positions(C, &[500], &[], &p);
        // the start clips to the contig end, the end is 510, so nothing is
        // covered and the head is a single zero run
        assert!(t.cursor() <= 1000);
        assert_eq!(t.integral(), 0.0);
    }

    #[test]
    fn endpoints_beyond_the_contig_are_clamped_not_dropped() {
        // a read at 900 with d=200 in a 1000 bp contig covers [900,1000)
        let t = pileup_from_positions(C, &[900], &[], &params(200, 1000));
        assert_eq!(t.cursor(), 1000);
        assert_eq!(t.integral(), 100.0);
    }

    // ---------- scaling and baseline ----------

    #[test]
    fn scale_factor_is_applied_in_f32() {
        let p = SingleEndParams::directional(200, 0, 1000, 0.5f32);
        let t = pileup_from_positions(C, &[100, 150], &[], &p);
        assert_eq!(
            spans(&t),
            vec![(0, 0.0), (100, 0.5), (150, 1.0), (300, 0.5)]
        );
    }

    #[test]
    fn baseline_is_a_floor_applied_after_scaling() {
        // lambda_bg = 0.5 floors every position, including the head
        let p = params(200, 1000).with_baseline(0.5);
        let t = pileup_from_positions(C, &[100], &[], &p);
        assert_eq!(cov(&t, 0), 0.5);
        assert_eq!(cov(&t, 150), 1.0);
    }

    #[test]
    fn zero_scale_floors_to_the_baseline() {
        let p = SingleEndParams::directional(200, 0, 1000, 0.0).with_baseline(2.5);
        let t = pileup_from_positions(C, &[100, 150], &[], &p);
        assert!(t.iter().all(|s| *s.value == 2.5));
    }

    // ---------- bidirectional ----------

    #[test]
    fn bidirectional_splits_d_evenly() {
        let p = SingleEndParams::bidirectional(200, 0, 1000, 1.0);
        // five_shift = 100, three_shift = 100 => a plus read at 100 -> [0, 200)
        let t = pileup_from_positions(C, &[100], &[], &p);
        assert_eq!(cov(&t, 0), 1.0);
        assert_eq!(cov(&t, 199), 1.0);
        assert_eq!(
            t.value_at(200),
            None,
            "past the cursor is genuinely uncovered"
        );
        assert_eq!(t.integral(), 200.0);
    }

    #[test]
    fn bidirectional_handles_odd_d_asymmetrically() {
        // d = 201 => five_shift = 100, three_shift = 101 => [0, 201)
        let p = SingleEndParams::bidirectional(201, 0, 1000, 1.0);
        let t = pileup_from_positions(C, &[100], &[], &p);
        assert_eq!(t.cursor(), 201);
        assert_eq!(t.integral(), 201.0);
    }

    // ---------- fragments ----------

    #[test]
    fn fragment_pileup_uses_real_intervals() {
        let frags = [(100u64, 200u64), (150, 250)];
        let t = pileup_from_fragments(C, &frags, 1000, 1.0, 0.0);
        // [0,100)=0, [100,150)=1, [150,200)=2, [200,250)=1
        assert_eq!(
            spans(&t),
            vec![(0, 0.0), (100, 1.0), (150, 2.0), (200, 1.0)]
        );
        assert_eq!(t.cursor(), 250);
        assert_eq!(t.integral(), 200.0);
    }

    #[test]
    fn fragment_pileup_ignores_the_estimated_d() {
        // a 2 bp fragment and a 400 bp fragment, side by side
        let frags = [(100u64, 102u64), (200, 600)];
        let t = pileup_from_fragments(C, &frags, 1000, 1.0, 0.0);
        assert_eq!(cov(&t, 101), 1.0);
        assert_eq!(cov(&t, 102), 0.0);
        assert_eq!(cov(&t, 599), 1.0);
        assert_eq!(t.integral(), 402.0);
    }

    #[test]
    fn fragment_pileup_clips_to_the_contig() {
        let frags = [(0u64, 10u64), (995, 5000)];
        let t = pileup_from_fragments(C, &frags, 1000, 1.0, 0.0);
        assert_eq!(t.cursor(), 1000);
        assert_eq!(t.integral(), 15.0);
    }

    #[test]
    fn weighted_fragments_accumulate_multiplicity() {
        let frags = [(100u64, 200u64, 3.0f32), (150, 250, 1.0)];
        let t = pileup_from_weighted_fragments(C, &frags, 1000, 1.0, 0.0);
        // [0,100)=0, [100,150)=3, [150,200)=4, [200,250)=1
        assert_eq!(
            spans(&t),
            vec![(0, 0.0), (100, 3.0), (150, 4.0), (200, 1.0)]
        );
        assert!((t.integral() - 400.0).abs() < 1e-3, "got {}", t.integral()); // 3*100 + 1*100
    }

    #[test]
    fn zero_weight_fragments_are_dropped() {
        let frags = [(100u64, 200u64, 0.0f32), (300, 400, 2.0)];
        let t = pileup_from_weighted_fragments(C, &frags, 1000, 1.0, 0.0);
        // the head up to 300 is an explicit zero run
        assert_eq!(spans(&t), vec![(0, 0.0), (300, 2.0)]);
    }

    #[test]
    fn weighted_and_unweighted_agree_for_unit_weights() {
        let frags = [(100u64, 200u64), (150, 250), (400, 500)];
        let w: Vec<(u64, u64, f32)> = frags.iter().map(|&(l, r)| (l, r, 1.0)).collect();
        let a = pileup_from_fragments(C, &frags, 1000, 1.0, 0.0);
        let b = pileup_from_weighted_fragments(C, &w, 1000, 1.0, 0.0);
        assert_eq!(spans(&a), spans(&b));
    }

    // ---------- naive ----------

    #[test]
    fn naive_quick_pileup_extends_symmetrically_and_clips() {
        // reads at 10 and 500 extended by 100 => [0,110) and [400,600)
        let t = naive_quick_pileup(C, &[10, 500], 100, 1000);
        assert_eq!(spans(&t), vec![(0, 1.0), (110, 0.0), (400, 1.0)]);
        assert_eq!(t.cursor(), 600);
        assert_eq!(cov(&t, 0), 1.0);
        assert_eq!(cov(&t, 109), 1.0);
        assert_eq!(cov(&t, 110), 0.0);
        assert_eq!(cov(&t, 399), 0.0);
        assert_eq!(cov(&t, 400), 1.0);
        assert_eq!(t.integral(), 310.0);
    }

    // ---------- invariants ----------

    #[test]
    fn integrated_depth_equals_the_sum_of_fragment_lengths() {
        let plus = [10u64, 40, 90, 140, 190, 240];
        let minus = [15u64, 45, 95, 145, 195, 245];
        let d = 200i64;
        let t = pileup_from_positions(C, &plus, &minus, &params(d, 1000));
        let want: f64 = plus
            .iter()
            .map(|&x| {
                let s = x;
                let e = ((x as i64) + d).min(1000) as u64;
                e.saturating_sub(s) as f64
            })
            .sum::<f64>()
            + minus
                .iter()
                .map(|&x| {
                    let s = ((x as i64) - d).max(0) as u64;
                    let e = (x as i64).min(1000) as u64;
                    e.saturating_sub(s) as f64
                })
                .sum::<f64>();
        // no coincident events in this fixture, so the mass is exactly conserved
        assert!(
            (t.integral() - want).abs() < 1e-3,
            "{} vs {want}",
            t.integral()
        );
    }

    #[test]
    fn canonical_form_holds_for_random_input() {
        let plus: Vec<u64> = (0..200).map(|i| (i * 37) % 900).collect();
        let minus: Vec<u64> = (0..200).map(|i| (i * 53) % 900).collect();
        let t = pileup_from_positions(C, &plus, &minus, &params(200, 1000));
        let mut prev_end = t.start();
        for (i, s) in t.iter().enumerate() {
            assert_eq!(s.start, prev_end, "hole or overlap at run {i}");
            assert!(s.end > s.start);
            assert!(*s.value >= 0.0);
            if i > 0 {
                let prev = t.runs()[i - 1].value;
                assert_ne!(prev, *s.value, "adjacent runs share a value at {i}");
            }
            prev_end = s.end;
        }
        assert_eq!(prev_end, t.cursor());
    }

    #[test]
    fn unsorted_input_gives_the_same_track_as_sorted() {
        let a = [500u64, 10, 300, 42, 700];
        let b = [100u64, 600, 250];
        let t1 = pileup_from_positions(C, &a, &b, &params(200, 1000));
        let mut a = a;
        let mut b = b;
        a.sort_unstable();
        b.sort_unstable();
        let t2 = pileup_from_positions(C, &a, &b, &params(200, 1000));
        assert!(t1.identical(&t2));
    }

    #[test]
    fn every_run_end_is_within_the_contig() {
        for rlength in [1u64, 2, 10, 1000, 100_000] {
            let plus: Vec<u64> = (0..50).map(|i| i * 7 % rlength.max(1)).collect();
            let t = pileup_from_positions(C, &plus, &[], &params(200, rlength));
            assert!(
                t.cursor() <= rlength,
                "rlength={rlength} cursor={}",
                t.cursor()
            );
        }
    }

    #[test]
    fn single_base_contig_is_handled() {
        let t = pileup_from_positions(C, &[0], &[], &params(200, 1));
        assert_eq!(t.end(), 1);
        assert!(t.cursor() <= 1);
        assert_eq!(t.integral(), 1.0);
    }

    #[test]
    fn zero_length_contig_gives_an_empty_track() {
        let t = pileup_from_positions(C, &[0], &[], &params(200, 0));
        assert!(t.is_empty());
    }
}
