//! The shared per-chromosome peak-calling core.
//!
//! This is the production path: given one chromosome's treatment pileup, its
//! optional merged control lambda, the q-score track built from the shared
//! histogram, it builds the paired position arrays, the above-cutoff chunks and
//! calls the narrow or broad peaks. It is the exact logic the differential
//! harness runs, lifted out of it (F125) so the harness and a shipping
//! `callpeak` binary share one implementation and cannot drift.
//!
//! Everything here is oracle-free: no `CALLPEAK_*` debug hooks, no
//! `peak_content` capture injection. Those stay in the harness.

use crate::{CallParams, Chunk, Peak};
use macs_core::{ChromId, Coord};
use macs_rle::SignalTrack;
use macs_score::{PScoreCache, PqTable};

/// The inputs needed to call one chromosome.
#[derive(Debug)]
pub struct ChromCall<'a> {
    /// Chromosome name, used for the output rows.
    pub name: &'a str,
    pub chrom: ChromId,
    /// Treatment pileup (`d_pileup_d`).
    pub treat: &'a SignalTrack<f32>,
    /// Merged control lambda (`ctrl_d_pileup_d`), if a control was supplied.
    pub ctrl: Option<&'a SignalTrack<f32>>,
    /// The q-score track for this chromosome, from the shared histogram.
    pub qtrack: &'a SignalTrack<f32>,
    /// The shared p->q table.
    pub table: &'a PqTable,
    /// The fragment length `d`, used for `min_length` and `max_gap`.
    pub d: Coord,
    /// F260: the coordinate upstream's sub-peak padding clamps at -- our stand-in
    /// for its signed `0`. See [`crate::regions::CallParams::clamp_floor`].
    pub clamp_floor: Coord,
    /// F271: where upstream's literal `0` lands in `chunks_at`'s
    /// `above_cutoff[0] == 0` branch. Zero everywhere except the `--nolambda`
    /// paired-end path; see [`crate::callpeak::call_chromosome`].
    pub zero_coord: Coord,
    /// `--qvalue` (natural units).
    pub qvalue: f64,
    /// `-log10(--pvalue)`, when `-p/--pvalue` was given.
    ///
    /// F189: upstream selects the scoring function on which cutoff is present --
    /// `call_peaks(['p',], [self.log_pvalue], ...)` if `log_pvalue is not None`,
    /// `call_peaks(['q',], [self.log_qvalue], ...)` otherwise -- and `-p` wins
    /// when both are somehow set. This port always scored against the q-value, so
    /// every `-p` invocation silently used `-log10(0.05)` instead of the requested
    /// p-value cutoff. The golden corpus has no `-p` variant, so this was
    /// invisible to the gate while being plainly wrong for CLI equivalence.
    pub p_cutoff: Option<f32>,
    /// The `maxgap` used to merge above-cutoff chunks.
    ///
    /// F188: `PeakDetect.__init__` reads `self.maxgap = opt.maxgap or opt.tsize`
    /// -- the measured **tag/fragment size**, not the predicted fragment length
    /// `d`. They coincide in paired-end mode and in `--extsize` single-end runs
    /// that happen to match the mean read length, which is why the divergence was
    /// invisible until single-end fixtures with a read length far from
    /// `--extsize` were added: there `max_gap` is the read length and a peak whose
    /// gap exceeds it is split even though the gap is under `d`.
    pub max_gap: Coord,
    pub broad: bool,
    pub broad_cutoff: f64,
    pub call_summits: bool,
    /// `lambda_bg`, the floor every control pileup is clamped to.
    ///
    /// With no control (or `--nolambda`) upstream pairs the treatment against the
    /// single-entry array `[lambda_bg]`, so the walk yields one paired entry.
    pub lambda_bg: f32,
}

/// Which sample gets scaled to which, `callpeak_cmd.py:219-259`.
///
/// ```python
/// if control and options.PE_MODE:
///     c1 = c1 * 2          # both ends of a PE fragment count towards background
/// if options.scaleto == "large":
///     options.tocontrol = (t1 <= c1)
/// else:
///     options.tocontrol = (t1 > c1)
/// ```
///
/// `t1`/`c1` are the **post**-filter totals -- upstream compares them *after*
/// `filter_dup`, which is why `--keep-dup` and `--scale-to` interact.
///
/// `--scale-to` defaults to **`small`**, not `large`: on
/// `frag_basic/barcode_fragments` (`t1 == c1 == 7995`) the default run and an
/// explicit `--scale-to small` both report an unscaled pileup of `7989`, while
/// `--scale-to large` reports `15978`. Treating the default as `large` (as an
/// earlier pass did) doubled every unscaled run's pileup.
pub fn to_control(scaleto_large: bool, t1: u64, c1: u64, pe_mode: bool) -> bool {
    // `c1 = c1 * 2` happens in paired-end mode only
    let c1 = if pe_mode { c1.saturating_mul(2) } else { c1 };
    if scaleto_large {
        t1 <= c1
    } else {
        t1 > c1
    }
}

/// `broad_level`: segment the lvl2 chunks with the **four-times** gap and close
/// each region as a broad peak (F113, `PeakDetect.py:264-266`).
fn broad_level(
    chunks: &[Chunk],
    params: &CallParams,
    scores: &[(&[f32], f32)],
    table: &PqTable,
    cache: &mut PScoreCache,
) -> Vec<Peak> {
    use crate::{close_peak_for_broad_region, segment_regions};
    let lvl2_params = CallParams {
        max_gap: params.max_gap.saturating_mul(4),
        ..params.clone()
    };
    let regions = segment_regions(chunks, lvl2_params.max_gap);
    let mut out = Vec::new();
    for r in &regions {
        // F159: the real p->q table and the p-score cache must be threaded
        // through. `close_peak_for_broad_region` derives the peak's q-score with
        // `__cal_qscore(tarray_pileup, tarray_control)`, which is a per-position
        // lookup into `self.pqtable` -- handing it an **empty** table makes every
        // lookup miss its zero fallback, so every broad peak was written with
        // `-log10(qvalue) = 0` and `score = int(10 * 0)`, on the whole corpus.
        if let Ok(p) = close_peak_for_broad_region(r, params, scores, table, Some(cache)) {
            out.push(p);
        }
    }
    out
}

/// A called peak, plus -- in broad mode -- the strong (lvl1) sub-peaks inside
/// it.
///
/// The lvl1 set is what the gappedPeak block structure is built from
/// (`__add_broadpeak`, `BedGraph.py:600`): one block per strong sub-peak, plus a
/// 1 bp block on each side when the strong set does not reach the broad region's
/// own ends. Keeping it on the result rather than recomputing it in the writer is
/// what makes the two files consistent with each other.
#[derive(Debug, Clone)]
pub struct Called {
    pub peak: Peak,
    /// `(start, end)` of each lvl1 peak inside `peak`, in coordinate order.
    /// Empty unless `broad`.
    pub lvl1: Vec<(Coord, Coord)>,
}

impl From<Peak> for Called {
    fn from(peak: Peak) -> Self {
        Called {
            peak,
            lvl1: Vec::new(),
        }
    }
}

/// `combine_broad`: emit one broad peak per lvl2 region, marking it broad
/// (F111, `CallPeakUnit.py:1925-1966`).
fn combine_broad(lvl2: &[Peak], lvl1: &[Peak]) -> Vec<Called> {
    let mut out = Vec::with_capacity(lvl2.len());
    let mut i = 0usize;
    for l2 in lvl2 {
        let mut pk = *l2;
        pk.broad = true;
        while i < lvl1.len() && lvl1[i].start < l2.start {
            i += 1;
        }
        // the strong sub-peaks inside this broad region, in order
        let mut inner = Vec::new();
        let mut j = i;
        while j < lvl1.len() && lvl1[j].end <= l2.end {
            inner.push((lvl1[j].start, lvl1[j].end));
            j += 1;
        }
        i = j;
        out.push(Called {
            peak: pk,
            lvl1: inner,
        });
    }
    out
}

/// Call peaks for one chromosome and return them in coordinate order.
///
/// Builds the paired position list (the union of the treatment and control run
/// ends, or the q-track's ends when there is no control), samples the q-score
/// at each position, keeps positions above `-log10(qvalue)` with a non-zero
/// treatment pileup as chunks, then calls the peak closer -- narrow, or the
/// two-level broad call when `broad` is set.
pub fn call_chromosome(cc: &ChromCall<'_>, cache: &mut PScoreCache) -> Vec<Called> {
    let params = CallParams {
        min_length: cc.d,
        max_gap: cc.max_gap,
        call_summits: cc.call_summits,
        clamp_floor: cc.clamp_floor,
        ..Default::default()
    };
    let neg_log10_qvalue = -(cc.qvalue.log10()) as f32;

    // The paired arrays, built exactly as `__chrom_pair_treat_ctrl`
    // (`CallPeakUnit.py:519`) does it: a head-to-head walk that emits
    // `(min(t_p[i], c_p[j]), t_v[i], c_v[j])` -- the *cursor* values, so a
    // position contributed by the treatment still carries whatever the control
    // run behind that position holds -- and stops as soon as either input is
    // exhausted.
    //
    // F150: this used to be a sorted, deduplicated union of run ends with the
    // control sampled as "the last run ending at or before `p`". Both halves of
    // that are wrong:
    //
    // * end-indexing the control returns the value *behind* `p`, where upstream
    //   returns the run the cursor still points at *ahead* of `p`. That moved the
    //   reported summit control by one step: `tiny/two_contigs` summit reported
    //   ctrl 7.4 (fold 5.35269) where the union value at that position is 7.6
    //   (fold 5.22821, upstream's answer);
    // * keeping the tail instead of stopping with the shorter array leaves paired
    //   indices that upstream never has.
    //
    // Both are invisible in a coordinate diff, which is why the boundary gates
    // stayed green while most of the corpus had a wrong q-value.
    let ctl = cc.ctrl;
    let (pos, tpos, cpos) = match ctl {
        Some(c) => {
            let (p, t, c) = crate::driver::pair_treat_ctrl_runs(cc.treat, c);
            (p, t, c)
        }
        // `--nolambda`, or a caller with no control: upstream substitutes
        //
        //     ctrl_pv = [treat_pv[0][-1:], np.array([self.lambda_bg], dtype="f4")]
        //
        // F186: that is a **one-entry** array, but its single *position* is the
        // treatment's **last** position, not the first. Running the ordinary union
        // walk against it therefore emits one row per treatment position -- the
        // treatment cursor always wins `min(t_p[i], c_p[0])`, and the control value
        // `lambda_bg` rides along unchanged -- until the treatment is exhausted and
        // the loop's `it < lt and ic < lc` bound stops it. Collapsing it to a single
        // row at the last position (which is what this used to do) left exactly one
        // scorable position, so `--nolambda` called **no peaks at all** on every
        // fixture whose treatment has more than one run
        // (`se_basic/gauss_two_peaks`, `se_dup/dup_rate_*`, ...).
        None => {
            let tp: Vec<Coord> = cc.treat.runs().iter().map(|r| r.end).collect();
            let tv: Vec<f32> = cc.treat.runs().iter().map(|r| r.value).collect();
            let n = tp.len();
            let mut p = Vec::with_capacity(n);
            let mut t = Vec::with_capacity(n);
            let mut c = Vec::with_capacity(n);
            for i in 0..n {
                p.push(tp[i]);
                t.push(tv[i]);
                c.push(cc.lambda_bg);
            }
            (p, t, c)
        }
    };
    let n = pos.len();

    // `__cal_qscore(treat_array, ctrl_array)`: the q-value is looked up per paired
    // index from the shared p->q table, **not** resampled from a separate
    // q-score track.
    let mut qpos = vec![0.0f32; n];
    let mut ppos = vec![0.0f32; n];
    for i in 0..n {
        let obs = if tpos[i] < 0.0 {
            0u32
        } else {
            (tpos[i] as i64).min(u32::MAX as i64) as u32
        };
        let p = macs_score::pscore(Some(cache), obs, cpos[i]);
        ppos[i] = p;
        qpos[i] = cc.table.qscore_or_zero(p);
    }

    // F114: each cutoff needs its own chunk list -- a lvl2 broad region must be
    // able to extend past the narrow cutoff.
    let chunks_at = |cut: f32| -> Vec<Chunk> {
        let mut v: Vec<Chunk> = Vec::new();
        for i in 0..n {
            if qpos[i] > cut && tpos[i] > 0.0 {
                v.push(Chunk {
                    // F271: upstream's `if above_cutoff[0] == 0:
                    // above_cutoff_startpos[0] = 0` assigns a *literal* zero, whose
                    // meaning depends on the code path. Measured: this branch fires in
                    // 282 of 473 sampled counted-PE cases and our literal `0` is correct
                    // for all of them. It is wrong only under `--nolambda`, where upstream
                    // replaces the control with a one-element array
                    // (`ctrl_pv = [treat_pv[0][-1:], [lambda_bg]]`), the paired array is
                    // walked without a real control cursor, and the literal zero is
                    // `coord_shift` in our shifted frame. See F271 for why this is keyed
                    // on the flag and not applied globally.
                    start: if i == 0 { cc.zero_coord } else { pos[i - 1] },
                    end: pos[i],
                    treat: tpos[i],
                    ctrl: cpos[i],
                    score_index: i,
                });
            }
        }
        v
    };
    // F189: `-p` scores against the p-value track, `-q` against the q-value one.
    // `__cal_pvalue_qvalue_table` still runs either way -- `call_peaks` computes it
    // before dispatching -- so the q column is populated even under `-p`.
    let use_p = cc.p_cutoff.is_some();
    let cut = cc.p_cutoff.unwrap_or(neg_log10_qvalue);
    let score_track: &[f32] = if use_p { &ppos } else { &qpos };
    let chunks = chunks_at(cut);
    let scores: Vec<(&[f32], f32)> = vec![(score_track, cut)];
    if cc.broad {
        // F111: broad mode calls the region set twice -- once at the q-value
        // cutoff (lvl1, strong) and once at `--broad-cutoff` (lvl2, broad) --
        // then combines them. The reported peak is the lvl2 region.
        let (l1, _) =
            crate::call_peaks_chromosome(cc.chrom, &chunks, &params, &scores, cc.table, cache);
        // F269: the chromosome loop in upstream's broad assembly is driven by **lvl1**:
        //
        //     chrs = lvl1peaks.get_chr_names()
        //     for chrom in sorted(chrs):
        //         lvl1peakschrom = lvl1peaks.get_data_from_chrom(chrom)
        //         lvl2peakschrom = lvl2peaks.get_data_from_chrom(chrom)
        //
        // So a chromosome with no lvl1 peak is skipped outright and its lvl2 regions are
        // never consulted. That is why upstream is silent rather than crashing here: the
        // lvl1 region on `sweep/gmini_mpe_d400_w180_ctrl_192` is `57..200` = 143 bases,
        // below `min_length` 145, so lvl1 is empty; and the `except StopIteration` branch
        // that would dereference an unbound `lvl2` is never reached, because the
        // chromosome never enters the loop at all.
        //
        // Our `combine_broad(&l2, &l1)` had no such condition and emitted the lvl2 region
        // (`54..200` = 146 >= 145) regardless. Since lvl2's above-cutoff set is a strict
        // superset of lvl1's (cutoff 1.0 vs 1.301), "lvl1 empty while lvl2 is not" is the
        // normal case whenever the strong level is empty -- so this guard, not a cutoff or
        // length tweak, is what makes the two agree.
        if l1.is_empty() {
            return Vec::new();
        }
        let qv: Vec<f32> = cc.qtrack.runs().iter().map(|r| r.value).collect();
        // F189: `call_broadpeaks` is dispatched on the same flag, so under `-p`
        // both levels are p-score cutoffs.
        let lvl2_cut = -(cc.broad_cutoff as f32).log10();
        let lvl2_scores: Vec<(&[f32], f32)> = if use_p {
            vec![(score_track, lvl2_cut)]
        } else {
            vec![(qv.as_slice(), lvl2_cut)]
        };
        let lvl2_chunks = chunks_at(lvl2_cut);
        // F266: dump the broad-path level sets. This is the intermediate the porting
        // plan lists as "candidate/merged intervals", and it is what
        // `sweep/gmini_mpe_d400_w180_ctrl_192 --broad` turns on: we emit one broad peak
        // where upstream emits none, and upstream builds a broad peak from a lvl2 region
        // *together with* the lvl1 regions inside it (`__add_broadpeak`,
        // `CallPeakUnit.py:2247`), so the two sets have to agree before the output can.
        // Neither is recoverable from the `*_peaks.xls`, because when they disagree the
        // golden xls is simply empty.
        if std::env::var_os("MACS3_RS_DUMP_BROAD_LEVELS").is_some() {
            let dump = |tag: &str, cs: &[crate::Chunk]| {
                let v: Vec<String> = cs
                    .iter()
                    .map(|c| format!("[{},{}]", c.start, c.end))
                    .collect();
                eprintln!("BROADLEVELS {} {} n={} {:?}", cc.name, tag, v.len(), v);
            };
            dump("lvl1", &chunks);
            dump("lvl2", &lvl2_chunks);
            eprintln!(
                "BROADLEVELS {} min_length={} max_gap={} lvl1_cut={:.6} lvl2_cut={:.6}",
                cc.name, cc.d, cc.max_gap, cut, lvl2_cut
            );
            // The per-position q-score around where lvl1 and lvl2 diverge. lvl1 starts
            // at 57 and lvl2 at 54, so if `54..56` sits within an ulp of `lvl1_cut` the
            // whole 2-base lvl1 region difference is a single-ulp q-track question.
            let lo = pos.iter().position(|&x| x >= 40).unwrap_or(0);
            let hi = pos.iter().position(|&x| x > 70).unwrap_or(n).min(n);
            for i in lo..hi {
                eprintln!(
                    "BROADQ {} pos={} treat={:.5} ctrl={:.9} p={:.9} q={:.9} {}",
                    cc.name,
                    pos[i],
                    tpos[i],
                    cpos[i],
                    ppos[i],
                    qpos[i],
                    if qpos[i] > cut { "ABOVE-L1" } else { "" }
                );
            }
        }
        let l2 = broad_level(&lvl2_chunks, &params, &lvl2_scores, cc.table, cache);
        return combine_broad(&l2, &l1);
    }

    let (p, _) = crate::call_peaks_chromosome(cc.chrom, &chunks, &params, &scores, cc.table, cache);
    p.into_iter().map(Called::from).collect()
}

// ---------------------------------------------------------------------------
// The production single-end pipeline
// ---------------------------------------------------------------------------

use crate::LambdaScales;
use macs_core::{MacsError, Result, Strand};
use macs_pileup::SingleEndParams;
use macs_track::SingleEndTrack;
use std::collections::BTreeMap;

/// One chromosome's prepared signals (the treatment pileup and the merged
/// control lambda).
#[derive(Debug, Clone)]
pub struct ChromSignals {
    pub name: String,
    /// Interned id of `name`, so the shared pipeline can look the chromosome up
    /// without holding on to whichever input track interned it first.
    pub chrom: ChromId,
    pub treat: SignalTrack<f32>,
    pub ctrl: Option<SignalTrack<f32>>,
}

/// Options for the single-end [`run_callpeak_se`] pipeline.
#[derive(Debug, Clone)]
pub struct SeConfig {
    /// `--extsize`: the fragment length `d`.
    pub extsize: i64,
    pub gsize: f64,
    pub slocal: i64,
    pub llocal: i64,
    pub qvalue: f64,
    pub call_summits: bool,
    pub broad: bool,
    pub broad_cutoff: f64,
    /// `--nolambda`: disable the dynamic local lambda entirely.
    pub nolambda: bool,
    /// `maxgap`, i.e. `opt.tsize` -- see [`ChromCall::max_gap`] (F188).
    pub max_gap: Coord,
    /// `-log10(--pvalue)`, when given (F189).
    pub p_cutoff: Option<f32>,
    /// `--scaleto large` (the default) rather than `small`; see [`to_control`].
    pub scaleto_large: bool,
    /// `--shift`: move every 5' end by this many bases before extending.
    ///
    /// Single-end only. `pileup_treat_ctrl_a_chromosome` passes `end_shift` to
    /// `self.treat.pileup_a_chromosome(..., directional=True, end_shift=...)` in
    /// the non-PE branch and to **neither** pileup in the PE branch --
    /// `PETrackI.pileup_a_chromosome` and `PETrackI.pileup_a_chromosome_c` take no
    /// shift argument at all -- so `--shift` is a no-op in paired-end mode.
    pub end_shift: i64,
}

/// Upstream's control scale factors and lambda windows, ported from the
/// harness's `control_scale_factors`.
fn control_scale_factors(
    ratio: f64,
    to_control: bool,
    d: i64,
    slocal: i64,
    llocal: i64,
) -> (Vec<i64>, Vec<f32>) {
    let mut scales = vec![d];
    let mut factors = vec![if to_control { 1.0f32 } else { ratio as f32 }];
    if slocal > 0 {
        scales.push(slocal);
        factors.push(((d as f64 / slocal as f64) * if to_control { 1.0 } else { ratio }) as f32);
    }
    // the llocal scale is omitted entirely unless `llocal > slocal`
    if llocal > slocal && llocal > 0 {
        scales.push(llocal);
        factors.push(((d as f64 / llocal as f64) * if to_control { 1.0 } else { ratio }) as f32);
    }
    (scales, factors)
}

/// `retadd_pscore`: the p-score track, as upstream's `min(p1,p2)` pointer walk
/// over the treatment and control tracks (F51). Treatment is truncated to an
/// integer and the control used as-is, with **no pseudocount** (F55) -- that is
/// `__cal_qscore`'s `get_pscore(int(a1), a2)`, not `ScoreTrack`'s
/// pseudocounted variant.
fn retadd_pscore(
    treat: &SignalTrack<f32>,
    ctrl: &SignalTrack<f32>,
    cache: &mut PScoreCache,
) -> SignalTrack<f32> {
    let ends = |t: &SignalTrack<f32>| -> Vec<Coord> { t.runs().iter().map(|r| r.end).collect() };
    let vals = |t: &SignalTrack<f32>| -> Vec<f32> { t.runs().iter().map(|r| r.value).collect() };
    let (p1s, v1s) = (ends(treat), vals(treat));
    let (p2s, v2s) = (ends(ctrl), vals(ctrl));
    let lo = treat.start().max(ctrl.start());
    let hi_lo = treat.end().min(ctrl.end());
    let cap = p1s.len() + p2s.len();
    let mut out = SignalTrack::empty(treat.chrom(), lo, hi_lo);
    let (mut i1, mut i2) = (0usize, 0usize);
    let mut last: Option<Coord> = None;
    while i1 < p1s.len() && i2 < p2s.len() && out.runs().len() < cap {
        let (p1, p2) = (p1s[i1], p2s[i2]);
        let (v1, v2) = (v1s[i1], v2s[i2]);
        let pos = p1.min(p2);
        if last != Some(pos) {
            let obs = (v1 as i64).clamp(0, u32::MAX as i64) as u32;
            out.push(pos, cache.get(obs, v2));
            last = Some(pos);
        }
        if p1 <= p2 {
            i1 += 1;
        }
        if p2 <= p1 {
            i2 += 1;
        }
    }
    out
}

/// Build the shared p->q table and per-chromosome q-score tracks from the
/// prepared signals. Mirrors upstream's `__pre_computes` histogram: each
/// chromosome's p-score track is folded into one AFDR walk.
/// Build the shared p->q table and the per-chromosome q-score tracks.
///
/// `paired` selects [`paired_pscore`] (F59) over the default [`retadd_pscore`];
/// see [`PeResult::paired_boundaries`] and F149 for why the default is the
/// latter.
pub fn build_qtable(signals: &[ChromSignals], paired: bool) -> (PqTable, Vec<SignalTrack<f32>>) {
    build_qtable_from(signals, paired, 0)
}

/// [`build_qtable`] without the per-chromosome q-score tracks.
///
/// Narrow-mode peak calling reads them nowhere, so dropping them is
/// output-identical and is what keeps the peak RSS at parity with upstream.
pub fn build_qtable_narrow(signals: &[ChromSignals], paired: bool) -> PqTable {
    build_qtable_from_with(signals, paired, 0, false).0
}

/// [`build_qtable`] with an explicit histogram origin.
///
/// `origin` is the coordinate the first span's length is measured from, which is
/// **not** the same as the track's start once a counted run has been shifted
/// (F154): upstream's `__cal_pvalue_qvalue_table` measures from coordinate zero,
/// so a shifted run must be measured from its own unshifted zero. See
/// [`macs_score::PScoreHistogram::add_track_from`].
pub fn build_qtable_from(
    signals: &[ChromSignals],
    paired: bool,
    origin: i64,
) -> (PqTable, Vec<SignalTrack<f32>>) {
    build_qtable_from_with(signals, paired, origin, true)
}

/// [`build_qtable_from`] with control over whether the per-chromosome q-score
/// tracks are materialised. See [`macs_score::QScoreSink::build_with_tracks`].
pub fn build_qtable_from_with(
    signals: &[ChromSignals],
    paired: bool,
    origin: i64,
    keep_qtracks: bool,
) -> (PqTable, Vec<SignalTrack<f32>>) {
    let parts = build_qtable_from_with_hist(
        signals,
        paired,
        origin,
        keep_qtracks,
        false,
        CutoffParams::off(),
    );
    (parts.table, parts.qtracks)
}

/// [`build_qtable_from_with`], optionally also returning the p->q histogram as
/// `(score, length)` pairs.
///
/// The histogram is what the release definition's 1e-6 q-score bound is actually a
/// statement about, and `dump_stages.py` records it as `qvalue_table`, so the stage
/// layer needs it. Building the pairs is free -- they are already in memory -- but
/// they are only materialised when asked for, so the shipping path is unchanged.
/// The outputs of [`build_qtable_from_with_hist`].
///
/// Named rather than spelled out inline: a four-element tuple is hard to read at the
/// call site and the ordering is easy to transpose.
#[allow(missing_debug_implementations)]
pub struct QTableParts {
    /// The p->q table.
    pub table: PqTable,
    /// Per-chromosome q-score tracks; empty unless requested.
    pub qtracks: Vec<SignalTrack<f32>>,
    /// The AFDR histogram as `(score, length)` pairs; empty unless requested.
    pub histogram: Vec<(f32, i64)>,
    /// `--cutoff-analysis` totals; all zero when the analysis was not requested.
    pub cutoffs: CutoffStats,
}

/// Extra inputs that only `--cutoff-analysis` needs.
///
/// Grouped so the shipping path keeps a five-argument signature and the analysis flags
/// do not leak into every caller.
#[derive(Debug, Clone, Copy, Default)]
pub struct CutoffParams {
    /// Run the analysis at all.
    pub enabled: bool,
    /// Merge gap: `opt.maxgap or opt.tsize` (F188).
    pub max_gap: Coord,
    /// Minimum peak length: `opt.d`.
    pub min_length: Coord,
}

impl CutoffParams {
    /// The defaults for a run without `--cutoff-analysis`.
    pub fn off() -> Self {
        CutoffParams::default()
    }
}

pub fn build_qtable_from_with_hist(
    signals: &[ChromSignals],
    paired: bool,
    origin: i64,
    keep_qtracks: bool,
    keep_histogram: bool,
    cut: CutoffParams,
) -> QTableParts {
    // Chromosome-level parallelism, output-identical to one thread:
    //
    // * `signals` is a slice, so `par_iter().map(..).collect()` yields the
    //   per-chromosome p-score tracks **in slice order** whatever the schedule;
    // * every chromosome gets its own `PScoreCache` via `map_init`, and the cache
    //   is pure memoisation, so it cannot change a score;
    // * the histogram is then folded sequentially in chromosome order. Even that
    //   is not order-sensitive -- the buckets are `i64` base-pair counters and
    //   integer addition is exact and order-independent -- but keeping it
    //   sequential makes that argument unnecessary.
    use rayon::prelude::*;
    let tracks: Vec<SignalTrack<f32>> = signals
        .par_iter()
        .map_init(
            PScoreCache::new,
            |cache: &mut PScoreCache, s: &ChromSignals| {
                match &s.ctrl {
                    // F59: in PE mode the chunk boundaries are the **coincident**
                    // ones -- `pos_array[above_cutoff]` over paired arrays. Using the
                    // single-end `retadd` walk instead inserts boundaries wherever only
                    // one track changes, which splits the chunk that wins the summit
                    // argmax and pushes summits low by up to half a split width.
                    Some(l) if paired => paired_pscore(&s.treat, l, Some(cache)),
                    Some(l) => retadd_pscore(&s.treat, l, cache),
                    None => crate::over_two_pv_array_track(&s.treat, &s.treat),
                }
            },
        )
        .collect();
    let mut sink = macs_score::QScoreSink::new();
    // The cutoff ladder is accumulated in the same pass, while each chromosome's
    // p-score track is in hand -- the analysis needs exactly this track, and computing
    // it twice would double the most expensive stage for no gain.
    let mut cut_stats: Option<(CutoffStats, Vec<f32>)> =
        cut.enabled.then(|| (CutoffStats::new(), cutoff_ladder()));
    for p in tracks {
        if let Some((stats, ladder)) = cut_stats.as_mut() {
            accumulate_cutoffs(stats, &p, ladder, cut.max_gap, cut.min_length);
        }
        sink.push_from(p, origin);
    }
    if let Some((_, ladder)) = cut_stats.as_ref() {
        seed_cutoffs(sink.histogram_mut(), ladder);
    }
    let (table, qtracks) = sink.build_with_tracks(keep_qtracks);
    let hist = if keep_histogram {
        sink.histogram_pairs()
    } else {
        Vec::new()
    };
    let cutoffs = cut_stats.map(|(s, _)| s).unwrap_or_default();
    QTableParts {
        table,
        qtracks,
        histogram: hist,
        cutoffs,
    }
}

/// `paired_pscore`: the p-score track on `self.chr_pos_treat_ctrl`, the
/// **union walk** of the treatment and control pileups.
///
/// F180: this used to emit only the positions where *both* tracks changed
/// ("coincident", F59) and resample each pair across the whole span. Upstream
/// emits one row per union position, carrying whatever each cursor currently
/// holds ([`paired_union`]), and it scores, histograms and chunks on exactly
/// that array -- `__cal_pvalue_qvalue_table`, `__cal_pscore`, `__cal_qscore` and
/// `pos_array[above_cutoff]` all index into `chr_pos_treat_ctrl`.
///
/// Coincidence hid the difference for most fixtures, because when both tracks
/// break at the same coordinates the two arrays agree. Where they do not, the
/// resampled version loses rows, and the loss is not cosmetic:
///
/// * the histogram loses the negative first span, so the rank `k` in the p->q
///   walk never goes negative and no q-value becomes `NaN`. Upstream *does*
///   produce `NaN` q-values there (`libc.log10` of a negative `k`), and `NaN >
///   cutoff` is false, so those positions are excluded from the above-cutoff
///   set. That is the whole reason `sweep/gmini_mfrag_d1200_w180_noc_123` starts
///   its peak at 16 upstream and at 0 here;
/// * the run count differs (175 here against 306 upstream on that fixture), so
///   even the total is only accidentally right.
pub fn paired_pscore(
    treat: &SignalTrack<f32>,
    ctrl: &SignalTrack<f32>,
    mut cache: Option<&mut PScoreCache>,
) -> SignalTrack<f32> {
    let (pos, tv, cv) = paired_union(treat, ctrl);
    let mut out = SignalTrack::empty(treat.chrom(), 0, 0);
    for i in 0..pos.len() {
        let obs = (tv[i] as i64).clamp(0, u32::MAX as i64) as u32;
        let sc = match cache.as_deref_mut() {
            Some(c) => c.get(obs, cv[i]),
            None => macs_score::pscore(None, obs, cv[i]),
        };
        if pos[i] > out.end() {
            out.push(pos[i], sc);
        }
    }
    out
}

/// The chromosome-independent SE setup: the local-lambda ladder, the contig length,
/// `lambda_bg`, the treatment scale factor, and the chromosome list.
///
/// `lambda_bg` is the control floor, which the caller needs for the no-control
/// pairing (`[lambda_bg]` as a one-entry control array).
///
/// Public so a caller can rebuild signals one chromosome at a time instead of holding
/// the whole genome's tracks -- see [`run_callpeak_se`] for why that matters.
pub fn se_setup(
    treat: &SingleEndTrack,
    ctrl: Option<&SingleEndTrack>,
    cfg: &SeConfig,
) -> (SeSignalSetup, f32, Vec<Vec<u8>>) {
    let rlength = u64::MAX / 2;
    let t_total = treat.total() as f64;
    let c_total = ctrl.map_or(0.0, |c| c.total() as f64);
    let treat_sum = t_total * cfg.extsize as f64;
    let control_sum = c_total * cfg.extsize as f64;
    let ratio = if control_sum > 0.0 {
        treat_sum / control_sum
    } else {
        0.0
    };
    // `tocontrol` scales the treatment down; set only when treatment is larger.
    //
    // F127: with no control `c_total` is 0, so the with-control test would set
    // `to_control` and then compute `lambda_bg = control_sum / gsize = 0` -- a zero
    // lambda at every base, which upstream treats as a fatal
    // `AssertionError`. `__call_peaks_wo_control` (`PeakDetect.py:1138-1141`)
    // hard-codes `treat_scale = 1.0` and
    // `lambda_bg = float(d) * treat_total / gsize`; take that branch instead.
    // `callpeak_cmd.py:242-259`; see [`to_control`]
    let to_control =
        ctrl.is_some() && to_control(cfg.scaleto_large, treat.total(), c_total as u64, false);
    // `PeakDetect.py:177-182` picks `control_sum / gsize` when scaling treatment
    // down to control, else `treat_sum / gsize`. With no control `control_sum` is
    // 0 and `tocontrol` is false, so the treatment sum is used -- the same value
    // `__call_peaks_wo_control` computes as `float(d) * treat_total / gsize`.
    let lambda_bg: f32 = if ctrl.is_none() {
        ((cfg.extsize as f64) * t_total / cfg.gsize) as f32
    } else if to_control {
        (control_sum / cfg.gsize) as f32
    } else {
        (treat_sum / cfg.gsize) as f32
    };
    let treat_scale: f32 = if to_control && ratio != 0.0 {
        (1.0 / ratio) as f32
    } else {
        1.0
    };
    // Without a control file, upstream does **not** go through
    // `__call_peaks_w_control` at all -- `PeakDetect.call_peaks` (`PeakDetect.py:109`)
    // dispatches on `if self.control:` and otherwise calls
    // `__call_peaks_wo_control` (`PeakDetect.py:1104`), which builds a
    // **single-scale** control:
    //
    // ```python
    // ctrl_scale_s = [float(self.d) / self.lregion,]
    // ctrl_d_s = [self.lregion,]
    // ```
    //
    // -- one `lregion`-wide window scaled by `d/lregion`, not the
    // `d`/`slocal`/`llocal` ladder. Using the ladder over-counts the control
    // and therefore over-states every p-score. `--nolambda` empties both lists,
    // which is the one-element `lambda_bg` array of `CallPeakUnit.py:621`.
    let (scales, factors) = if ctrl.is_none() {
        if cfg.nolambda || cfg.llocal <= 0 {
            (Vec::new(), Vec::new())
        } else {
            (
                vec![cfg.llocal],
                vec![(cfg.extsize as f32) / (cfg.llocal as f32)],
            )
        }
    } else {
        control_scale_factors(ratio, to_control, cfg.extsize, cfg.slocal, cfg.llocal)
    };
    let lscales = LambdaScales {
        d: scales[0],
        slocal: *scales.get(1).unwrap_or(&0),
        llocal: *scales.get(2).unwrap_or(&0),
        d_factor: factors[0],
        slocal_factor: factors.get(1).copied().unwrap_or(0.0),
        llocal_factor: factors.get(2).copied().unwrap_or(0.0),
    };

    // F243: the SE counterpart of the F236 diagnostic. The SE with-control path is
    // the next cluster to attack and it needs the same inputs: which scale wins, and
    // what each factor is. `coalesce_into` already dumps the f32 values that enter
    // `%.5f`; this prints the ladder that produced them.
    if std::env::var_os("MACS3_RS_DUMP_CTRL_LAMBDA").is_some() {
        eprintln!(
            "SE lambda scales: d={} slocal={} llocal={} factors={:?} lambda_bg={:?} treat_scale={:?}",
            scales[0],
            *scales.get(1).unwrap_or(&0),
            *scales.get(2).unwrap_or(&0),
            factors,
            lambda_bg,
            treat_scale,
        );
    }

    let setup = SeSignalSetup {
        lscales,
        rlength,
        lambda_bg,
        treat_scale,
        ratio_treat2control: ratio,
        to_control,
        sregion: cfg.slocal,
        lregion: cfg.llocal,
    };
    (setup, lambda_bg, chromosome_names(treat, ctrl))
}

/// Build the treatment pileup and merged control lambda for **every** chromosome.
///
/// A thin wrapper over [`se_setup`] + [`build_all_se_chromosomes`], so the
/// chromosome-independent arithmetic lives in exactly one place and the streaming
/// pipeline ([`run_callpeak_se`]) shares it. Prefer the streaming form when memory
/// matters: this one holds every chromosome's tracks at once.
/// Build the treatment pileup and merged control lambda for **every** chromosome.
///
/// A thin wrapper over [`se_setup`] + [`build_all_se_chromosomes`], so the
/// chromosome-independent arithmetic lives in exactly one place and the streaming
/// pipeline ([`run_callpeak_se`]) shares it. Prefer the streaming form when memory
/// matters: this one holds every chromosome's tracks at once.
pub fn build_signals_se(
    treat: &SingleEndTrack,
    ctrl: Option<&SingleEndTrack>,
    cfg: &SeConfig,
) -> (Vec<ChromSignals>, f32) {
    let (setup, lambda_bg, names) = se_setup(treat, ctrl, cfg);
    (
        build_all_se_chromosomes(&names, treat, ctrl, cfg, &setup),
        lambda_bg,
    )
}

/// The chromosome names upstream iterates, sorted.
///
/// `PeakDetect.call_peaks` walks `sorted(self.treat.pileup_control...)`-style keys;
/// the union of both inputs' names, sorted by name, is what the golden corpus is
/// built on.
pub fn chromosome_names(treat: &SingleEndTrack, ctrl: Option<&SingleEndTrack>) -> Vec<Vec<u8>> {
    let mut names: Vec<Vec<u8>> = treat
        .positions()
        .chroms()
        .iter()
        .map(|&c| treat.genome().name(c).to_vec())
        .collect();
    if let Some(c) = ctrl {
        for chrom in c.positions().chroms() {
            let n = c.genome().name(chrom).to_vec();
            if !names.contains(&n) {
                names.push(n);
            }
        }
    }
    names.sort();
    names
}

/// Everything [`build_one_se_chromosome`] needs beyond the input tracks.
#[derive(Debug, Clone, Copy)]
pub struct SeSignalSetup {
    lscales: LambdaScales,
    rlength: Coord,
    pub lambda_bg: f32,
    treat_scale: f32,
    /// `PeakDetect.py:177-182`: the treatment/control size ratio, for the stage dump.
    pub ratio_treat2control: f64,
    /// Whether the treatment is scaled down to the control.
    pub to_control: bool,
    /// `--slocal`, the small local window.
    pub sregion: i64,
    /// `--llocal`, the large local window.
    pub lregion: i64,
}

/// Build the pileup and merged control lambda for **every** chromosome at once.
///
/// Convenient, and what the differential harness wants, but it holds
/// `2 x chromosomes` signal tracks simultaneously: on a 24-chromosome 200 kb test
/// that is ~210 MB, which is the whole reason peak RSS overshot upstream. Prefer
/// [`stream_se_peaks`] in the shipping pipeline.
pub fn build_all_se_chromosomes(
    names: &[Vec<u8>],
    treat: &SingleEndTrack,
    ctrl: Option<&SingleEndTrack>,
    cfg: &SeConfig,
    setup: &SeSignalSetup,
) -> Vec<ChromSignals> {
    use rayon::prelude::*;
    // Output-identical to one thread: `names` is already in upstream's iteration
    // order (sorted by name) and `par_iter().filter_map(..).collect()` preserves it,
    // while each chromosome's tracks are built purely from its own reads -- no shared
    // mutable state crosses the boundary.
    names
        .par_iter()
        .filter_map(|name| {
            let chrom = treat.genome().get(name)?;
            build_one_se_chromosome(chrom, name, treat, ctrl, cfg, setup)
        })
        .collect()
}

/// The single-end per-chromosome pileup plus merged control lambda.
///
/// Split out of [`build_signals_se`] so the chromosome loop can be parallel; it
/// touches nothing but its own arguments.
pub fn build_one_se_chromosome(
    chrom: macs_core::ChromId,
    name: &[u8],
    treat: &SingleEndTrack,
    ctrl: Option<&SingleEndTrack>,
    cfg: &SeConfig,
    setup: &SeSignalSetup,
) -> Option<ChromSignals> {
    let SeSignalSetup {
        lscales,
        rlength,
        lambda_bg,
        treat_scale,
        ..
    } = *setup;
    {
        let t = macs_pileup::pileup_from_positions(
            chrom,
            treat.positions().strand(chrom, Strand::Plus),
            treat.positions().strand(chrom, Strand::Minus),
            &SingleEndParams::directional(cfg.extsize, cfg.end_shift, rlength, treat_scale),
        );
        // Upstream's `no_lambda_flag` path (`CallPeakUnit.py:621-622`) replaces the
        // whole control pileup with a **single-element** array:
        //
        // ```python
        // ctrl_pv = [treat_pv[0][-1:], np.array([self.lambda_bg,], dtype="f4")]
        // ```
        //
        // i.e. one position -- the treatment track's *last* position -- carrying
        // `lambda_bg`. That matters: `__chrom_pair_treat_ctrl` walks
        // `while it < lt and ic < lc`, so with `lc == 1` it emits exactly one
        // triple, not one per treatment run. Reproducing it as a one-run track
        // ending at the treatment's final position makes the shared
        // `retadd_pscore` walk behave identically.
        //
        // The previous code instead paired the entire treatment track against
        // itself (`over_two_pv_array_track(&treat, &treat)`), which both used the
        // wrong lambda -- the treatment pileup, not `lambda_bg` -- and left the
        // trailing intervals with a control of `0.0`. `get_pscore(k, 0.0)` is an
        // `AssertionError` upstream, so that reached the zero-lambda abort (F127)
        // on every no-control fixture: `se_edge/no_control` and all 100+ of the
        // `sweep/*_noc_*` fixtures.
        let lambda = match ctrl {
            // F152: no control **file** is not the same as `--nolambda`.
            //
            // `--nolambda` empties `ctrl_d_s`/`ctrl_scaling_factor_s`, which sets
            // `no_lambda_flag` and replaces the control pileup with the
            // one-element `[lambda_bg]` array. With no control *file*, upstream
            // dispatches to `__call_peaks_wo_control` (`PeakDetect.py:1104`),
            // which builds a **single-scale control from the treatment itself**:
            //
            //     ctrl_scale_s      = [float(self.d) / self.lregion,]
            //     ctrl_d_s          = [self.lregion,]
            //
            // -- one `llocal`-wide window scaled by `d/llocal`. Collapsing both
            // cases onto the one-element array left every no-control fixture with
            // a single paired index, so the q-value could never clear
            // `-log10(0.05)` and no peak was ever called: all 100+ `*_noc_*`
            // sweep fixtures plus `se_edge/no_control`.
            // F162/F186: `--nolambda` sets `no_lambda_flag` from the **scale lists
            // being empty**, not from the absence of a control file, so the guard
            // has to be on the flag alone. Matching it on `None` as well meant a
            // `--nolambda` invocation *with* `-c` kept the real local lambda: the
            // `-B` control track then had hundreds of runs instead of upstream's
            // single constant row, and `--nolambda` called no peaks at all
            // (`se_basic/gauss_two_peaks`, `se_dup/dup_rate_*`).
            _ if cfg.nolambda || lscales.d <= 0 => {
                // F164: see the paired-end note -- `push(t.end(), ..)` on a track
                // that already starts at `t.end()` is a no-op, so the control came
                // out empty.
                Some(macs_rle::SignalTrack::from_runs_exact(
                    chrom,
                    0,
                    t.end().max(1),
                    vec![macs_rle::Run::new(t.end().max(1), lambda_bg)],
                ))
            }
            None => {
                let mut combined: Option<SignalTrack<f32>> = None;
                for scale in lscales.as_pairs() {
                    let p = macs_pileup::pileup_from_positions(
                        chrom,
                        treat.positions().strand(chrom, Strand::Plus),
                        treat.positions().strand(chrom, Strand::Minus),
                        &SingleEndParams::bidirectional(scale.d, 0, rlength, scale.scale_factor)
                            .with_baseline(lambda_bg),
                    );
                    combined = Some(match combined {
                        None => p,
                        Some(prev) => crate::over_two_pv_array_track(&prev, &p),
                    });
                }
                combined
            }
            Some(c) => {
                let cc = c.genome().get(name)?;
                let mut combined: Option<SignalTrack<f32>> = None;
                for scale in lscales.as_pairs() {
                    let p = macs_pileup::pileup_from_positions(
                        cc,
                        c.positions().strand(cc, Strand::Plus),
                        c.positions().strand(cc, Strand::Minus),
                        &SingleEndParams::bidirectional(scale.d, 0, rlength, scale.scale_factor)
                            .with_baseline(lambda_bg),
                    );
                    combined = Some(match combined {
                        None => p,
                        Some(prev) => crate::over_two_pv_array_track(&prev, &p),
                    });
                }
                combined
            }
        };
        Some(ChromSignals {
            name: String::from_utf8_lossy(name).into_owned(),
            chrom,
            treat: t,
            ctrl: lambda,
        })
    }
}

/// Run the full single-end `callpeak` pipeline on loaded tracks and return the
/// peaks as `(chromosome name, peak)` in coordinate order.
///
/// This is the shipping implementation: build the treatment pileup and merged
/// control lambda, fold the p-score tracks into one q-value histogram, then
/// call each chromosome. It is byte-for-byte the logic the differential harness
/// runs, so the golden parity gates cover it directly.
pub fn run_callpeak_se(
    treat: &SingleEndTrack,
    ctrl: Option<&SingleEndTrack>,
    cfg: &SeConfig,
) -> Vec<(String, Peak)> {
    // Single pass, materialising every chromosome's signals.
    //
    // A streaming two-pass variant was implemented and measured, then removed: it
    // bounds *live* memory to one chromosome but peak RSS is dominated by read count,
    // not genome size (docs/status.md, "Where the memory actually goes"), so it moved
    // the number by ~5 MB while costing a second pass over the pileup -- 1.58 s ->
    // 2.04 s on the benchmark, against a >=3x wall-clock criterion. Losing a tenth of
    // the speed margin to chase memory that is not the bottleneck is the wrong trade,
    // and the real fix has to shrink the read representation instead.
    let (signals, lambda_bg) = build_signals_se(treat, ctrl, cfg);
    let (table, qtracks) = build_qtable(&signals, false);
    let d = cfg.extsize.max(1) as Coord;
    let max_gap = if cfg.max_gap == 0 { d } else { cfg.max_gap };
    let mut cache = PScoreCache::new();
    let mut peaks: Vec<(String, Peak)> = Vec::new();
    for (k, s) in signals.iter().enumerate() {
        let cc = ChromCall {
            name: &s.name,
            chrom: s.chrom,
            treat: &s.treat,
            ctrl: s.ctrl.as_ref(),
            // F260: single-end tracks keep upstream's coordinates as-is, so its
            // signed `0` is our `0`.
            clamp_floor: 0,
            // F271: unshifted path.
            zero_coord: 0,
            qtrack: &qtracks[k],
            table: &table,
            d,
            max_gap,
            p_cutoff: cfg.p_cutoff,
            qvalue: cfg.qvalue,
            broad: cfg.broad,
            broad_cutoff: cfg.broad_cutoff,
            call_summits: cfg.call_summits,
            lambda_bg,
        };
        for pk in call_chromosome(&cc, &mut cache) {
            peaks.push((s.name.clone(), pk.peak));
        }
    }
    peaks
}

/// The single-end p-score track.
///
/// The same walk as `build_qtable_from_with`'s per-chromosome step, exposed so a
/// caller that only needs one chromosome's track does not have to materialise the
/// rest.
pub fn se_pscore_track(sig: &ChromSignals, cache: &mut PScoreCache) -> macs_rle::SignalTrack<f32> {
    match &sig.ctrl {
        Some(l) => retadd_pscore(&sig.treat, l, cache),
        None => crate::over_two_pv_array_track(&sig.treat, &sig.treat),
    }
}

/// The same track, reduced to a histogram fragment so pass 1 can drop it.
pub fn se_pscore_histogram(sig: &ChromSignals, origin: i64) -> macs_score::PScoreHistogram {
    let p = se_pscore_track(sig, &mut PScoreCache::new());
    let mut h = macs_score::PScoreHistogram::new();
    h.add_track_from(&p, origin);
    h
}

/// Sanity: the pipeline must not panic on an empty input, and `MacsError` is
/// re-exported so the binary can surface parse failures uniformly.
#[allow(dead_code)]
fn _err(_e: MacsError) -> Result<()> {
    Ok(())
}

/// Options for the paired-end [`run_callpeak_pe`] pipeline.
///
/// Mirrors [`SeConfig`], but the two fragment lengths have to be threaded
/// separately (F95):
///
/// * `d` is `opt.tsize`, the mean template length **as read**, truncated. It is
///   the first scale's *extent* and what `min_length`/`max_gap` see.
/// * `avg_tl` is the **post-filter** mean. The wide-window scale factors are
///   built from it, matching upstream's local `d` at `PeakDetect.py:156`.
///
/// Using one number for both shifts every local lambda by the ratio of the two.
#[derive(Debug, Clone)]
pub struct PeConfig {
    /// `--extsize` in paired-end mode: the mean template length as read.
    pub tsize: f64,
    /// The same number **untruncated**, for the wide-window scale factors.
    ///
    /// F166: `ctrl_scale_s` is built with `float(self.d)`, so it sees the
    /// fractional value, while `ctrl_d_s = [self.d]`, `min_length`, `max_gap` and
    /// the `# d = %d` header all go through `cython.int` and see the truncation.
    /// On `frag_basic/barcode_fragments` that is `98.307` against `98`, and the
    /// `--llocal` factor differs by 0.1% -- enough to move `-log10(pvalue)` in the
    /// third decimal. For `--format FRAG` the underlying value is also the
    /// *unweighted* mean of the row lengths, not the count-weighted one: see
    /// [`macs_peaks::callpeak::pe_tsize`].
    pub tsize_exact: f64,
    /// Effective genome size, `-g`.
    pub gsize: f64,
    /// `--slocal`.
    pub slocal: i64,
    /// `--llocal`.
    pub llocal: i64,
    /// `--qvalue`.
    pub qvalue: f64,
    /// `--call-summits`.
    pub call_summits: bool,
    /// `--broad`.
    pub broad: bool,
    /// `--broadcutoff`.
    pub broad_cutoff: f64,
    /// `--nolambda`.
    pub nolambda: bool,
    /// `--scaleto large` (the default) rather than `small`; see [`to_control`].
    pub scaleto_large: bool,
}

/// The paired-end **control** projection: extend each control fragment end
/// symmetrically by `d / 2`.
///
/// `pileup_treat_ctrl_a_chromosome` (`CallPeakUnit.py:606`) calls
/// `pileup_a_chromosome_c(chrom, ctrl_d_s, ctrl_scaling_factor_s,
/// baseline_value=lambda_bg)` for `PE_MODE` and **omits** `directional`, whose
/// declared default is `True` -- which would give `five_shift = 0`,
/// `three_shift = d`. That is what the source reads as.
///
/// It is nevertheless **wrong**: switching the PE control to
/// `directional` makes `gonechrom_mpe_d400_w600_ctrl_383`'s `-log10(pvalue)`
/// worse (100.756 against upstream's 111.952, versus 111.283 for the symmetric
/// projection). So the shipped `.so` behaves as if `directional` were `False`
/// here, i.e. `d//2` / `d - d//2`, and the symmetric projection is retained.
/// See F150 for the open question of which it really is.
fn pe_ctrl_params(d: i64, rlength: Coord, scale_factor: f32) -> SingleEndParams {
    SingleEndParams::symmetric(d, rlength, scale_factor)
}

/// F187: the **counted** control's projection is one base wider than the
/// uncounted one, and the two must not share a parameter set.
///
/// `PETrackI.pileup_a_chromosome_c` calls `pileup_from_PN_shifted` with
/// `five_shift = three_shift = d//2`, so an end at `x` covers `[x - d/2, x + d/2)`
/// -- a span of `d` for even `d` and `d - 1` for odd `d`. That is
/// [`pe_ctrl_params`].
///
/// `PETrackII.pileup_a_chromosome_c` instead calls
/// `pileup_from_LRC_centers_as_list`, which builds `start_poss` as
/// `concat(l - d//2, r - d//2)` and then sets `end_poss = start_poss + d`, so the
/// window is
///
/// so the window is `[x - d/2, x - d/2 + d)` -- a span of exactly `d`, i.e.
/// `three_shift = d - d//2`. Using the narrow split for `--format FRAG` made
/// every right-hand window one base short, which showed up as a small count
/// deficit growing with position: on
/// `sweep/gmini_mfrag_d1200_w600_ctrl_231` (`d = 119`) the merged lambda at
/// position 60 came out 1424.5 against upstream's 1429.5, and by position 68 the
/// two were 14 apart.
fn pe_ctrl_params_counted(d: i64, rlength: Coord, scale_factor: f32) -> SingleEndParams {
    SingleEndParams::centred(d, rlength, scale_factor)
}

/// The paired-end result: per-chromosome signals plus the reported `d`.
///
/// `d` is returned separately because it is the truncated **as-read** mean and
/// not the post-filter mean the scale factors use.
#[derive(Debug, Clone)]
pub struct PeResult {
    /// One entry per chromosome, in name order.
    pub signals: Vec<ChromSignals>,
    /// `--nomodel`'s `d`: `tsize` truncated.
    pub d: Coord,
    /// Whether to build the p-score track on coincident boundaries (F59).
    ///
    /// `false`, the value [`run_callpeak_pe`] sets, uses the single-end `retadd`
    /// union-boundary walk -- which is what the differential gate is calibrated
    /// against, and what reproduces upstream on this corpus. `true` uses
    /// [`paired_pscore`], F59's stricter reading of `pos_array[above_cutoff]`;
    /// it is reachable only via the `CALLPEAK_BOUNDARY=paired` harness switch
    /// because it moves `-log10(pvalue)` in the 4th decimal on some fixtures.
    /// See F149.
    pub paired_boundaries: bool,
    /// The control floor, `lambda_bg`. Reported so a caller can build the
    /// no-control pairing without recomputing it.
    pub lambda_bg: f32,
    /// Add this to every coordinate to get the real one.
    ///
    /// F151/F154: `--format FRAG` controls are unclipped `d`-windows centred on
    /// each fragment *end* (`x - d//2`), which for a contig-local fragment is
    /// **negative** -- `-5000` with `--largelocal 10000`. A `u64` position cannot
    /// hold that, so a counted run shifts every coordinate by `llocal / 2` before
    /// building signals and the caller subtracts it again when writing. The shift
    /// is uniform, so no ordering, depth, score or `d` comparison changes; only
    /// the printed coordinates do, and those are what upstream prints unclipped.
    ///
    /// `0` for every uncounted run, so nothing else has to think about it.
    pub coord_shift: i64,
}

/// The paired treatment/control union walk: `__chrom_pair_treat_ctrl`
/// (`CallPeakUnit.py:519`).
///
/// Merges the two pileups onto the union of their run boundaries. Each emitted
/// row is `(min(t_p[i], c_p[j]), t_v[i], c_v[j])` -- note the *values* are the
/// ones each track currently carries, not the values at `min`. When one track's
/// boundary is strictly smaller the other track is neither advanced nor sampled
/// again, so a track's value is repeated across every boundary the other track
/// contributes.
///
/// This is the array `self.chr_pos_treat_ctrl`, and it is what `-B` writes: the
/// treatment pileup and control lambda bedGraphs are two coalesced views of the
/// same rows, not the raw pileups. Getting this wrong is why an earlier `-B`
/// implementation produced 2821 runs where upstream produced 336.
///
/// The walk stops when either input is exhausted, matching the `it < lt and
/// ic < lc` loop bound.
pub fn paired_union(
    treat: &SignalTrack<f32>,
    ctrl: &SignalTrack<f32>,
) -> (Vec<Coord>, Vec<f32>, Vec<f32>) {
    let tp: Vec<Coord> = treat.runs().iter().map(|r| r.end).collect();
    let tv: Vec<f32> = treat.runs().iter().map(|r| r.value).collect();
    let cp: Vec<Coord> = ctrl.runs().iter().map(|r| r.end).collect();
    let cv: Vec<f32> = ctrl.runs().iter().map(|r| r.value).collect();
    let mut pos = Vec::with_capacity(tp.len() + cp.len());
    let mut t = Vec::with_capacity(tp.len() + cp.len());
    let mut c = Vec::with_capacity(tp.len() + cp.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < tp.len() && j < cp.len() {
        let (p1, p2) = (tp[i], cp[j]);
        pos.push(p1.min(p2));
        t.push(tv[i]);
        c.push(cv[j]);
        if p1 <= p2 {
            i += 1;
        }
        if p2 <= p1 {
            j += 1;
        }
    }
    (pos, t, c)
}

/// `pileup_from_fragments`'s third argument is the track **length**, which clamps
/// fragment endpoints. It is not the fragment length `d`.
const RLENGTH: Coord = u64::MAX / 2;

/// Paired-end peak calling: fragment pileups for treatment and control.
///
/// This is the shared implementation behind both the `callpeak` CLI and the
/// `callpeak-e2e` gate, so the two cannot drift. The comments on the scale
/// factors, the doubled control count and the no-control path are the findings
/// that made the gate pass; see F90, F94, F95, F98, F101, F102 and F127 in
/// `docs/upstream-findings.md` for the measurements behind them.
///
/// # The control is counted at both fragment ends
///
/// `control_total` is doubled and `control_sum` is
/// `control.total * 2 * avg_template_length`. Upstream's comment reads
/// "entire fragment is counted as 1 in treatment whereas both ends of fragment
/// are counted in control/input", and `callpeak_cmd.py:216` doubles `c1` for the
/// same reason. Unhalved, the control lands nowhere near the local lambda
/// upstream computes.
///
/// # Each control scale is projected from both ends, symmetrically
///
/// The control's fragment *ends* are treated as single-end tags, each extended
/// to the scale length with `d//2` either side (`SingleEndParams::symmetric`),
/// not with the single-end `shift`/`extension` split (F101). Feeding reversed
/// `(r, l)` pairs to `pileup_from_fragments` does not work: it assumes
/// `start <= end`, so the sweep misplaces them and the control collapses.
///
/// # Without a control the whole path changes
///
/// `c_total` is 0, so the with-control formulas would give `lambda_bg = 0` (a
/// zero lambda everywhere) and divide `treat_scale` by zero. Upstream's
/// `__call_peaks_wo_control` (`PeakDetect.py:1134`) instead hard-codes
/// `treat_scale = 1.0` and, in PE mode, uses
/// `ctrl_scale_s = [treat_length / (lregion * treat_total * 2)]` with
/// `ctrl_d_s = [lregion]`, piling the *treatment's* ends up as the control.
pub fn run_callpeak_pe(
    treat: &macs_track::FragmentTrack,
    ctrl: Option<&macs_track::FragmentTrack>,
    cfg: &PeConfig,
) -> PeResult {
    let avg_tl = treat.average_template_length();
    // F92: `d` is the as-read mean **truncated**, not rounded.
    let d = (cfg.tsize as u64).max(1);
    // F154: shift the whole run so a counted control's negative window starts
    // become representable. Half the widest window is the largest negative
    // excursion the centred projection can make (`x - d//2` with `x >= 0`), so
    // after the shift every position is >= 0 and the caller subtracts the shift
    // again when writing coordinates.
    let coord_shift: i64 = if treat.has_counts() {
        (d as i64).max(cfg.slocal).max(cfg.llocal) / 2
    } else {
        0
    };
    let sh = |v: Coord| -> Coord { (v as i64 + coord_shift).max(0) as Coord };
    let frags_of = |t: &macs_track::FragmentTrack, c: ChromId| -> Vec<(Coord, Coord)> {
        t.frags(c)
            .iter()
            .map(|f| (sh(f.start), sh(f.end)))
            .collect()
    };
    // `--format FRAG` (a counted `PETrackII`): every interval carries a
    // multiplicity, so the treatment coverage is `[l, r)` **weighted** by the
    // count (`pileup_from_LRC_as_list`), not a unit sweep.
    let wfrags_of = |t: &macs_track::FragmentTrack, c: ChromId| -> Vec<(Coord, Coord, f32)> {
        if !t.has_counts() {
            return Vec::new();
        }
        t.frags(c)
            .iter()
            .zip(t.counts(c).iter())
            .map(|(f, &w)| (sh(f.start), sh(f.end), f32::from(w)))
            .collect()
    };

    // `control_total = self.control.total * 2` -- both ends of a paired-end
    // fragment count towards the background (`PeakDetect.py:154`).
    //
    // F184: `treat_sum` and `control_sum` are declared `cython.long` in the source,
    // so the compiled module assigns them through `__Pyx_PyLong_As_long`, which
    // goes via `__Pyx_PyNumber_IntOrLong` and therefore **truncates the product**.
    // `average_template_length` is a `cython.float`, so in PE mode
    //
    //     control_sum = int(control_total * average_template_length)
    //
    // loses its fractional part before the division. On
    // `sweep/gmini_mpe_d1200_w600_ctrl_066` that is `int(394 * 99.79032135009766)`
    // = 39317 rather than 39317.386..., moving `ratio_treat2control` from
    // 0.4720812240955015 to 0.47208586616476333 -- a 9.83e-6 relative shift that
    // scales every control lambda and shows up as a fourth-decimal difference in
    // `-log10(p)`/`-log10(q)` on every paired-end fixture. `lambda_bg` divides the
    // same truncated `control_sum`, so it has to be truncated here too.
    let c_total = ctrl.map(|c| c.total()).unwrap_or(0);
    let treat_sum = treat.length();
    let control_total = c_total.saturating_mul(2);
    let trunc = |v: f64| (v as i64) as f64;
    let control_sum = trunc(control_total as f64 * avg_tl);
    let ratio = if control_sum > 0.0 {
        treat_sum as f64 / control_sum
    } else {
        0.0
    };

    // `callpeak_cmd.py:216` doubles the control count before this comparison, and
    // `--scaleto small` scales the larger sample down to the smaller.
    let no_control = ctrl.is_none();
    // see [`to_control`]
    let to_control =
        !no_control && to_control(cfg.scaleto_large, treat.total(), control_total, true);
    let lambda_bg: f32 = if no_control {
        (treat_sum as f64 / cfg.gsize) as f32
    } else if to_control {
        (control_sum / cfg.gsize) as f32
    } else {
        (treat_sum as f64 / cfg.gsize) as f32
    };
    let treat_scale: f32 = if no_control {
        1.0
    } else if to_control && ratio != 0.0 {
        (1.0 / ratio) as f32
    } else {
        1.0
    };

    // F156: `PeakDetect.__init__` (`PeakDetect.py:184-224`) builds
    //
    //     ctrl_d_s    = [self.d]                 # `self.d`, not the local `d`
    //     ctrl_scale_s = [ratio]                              if not tocontrol
    //                   + [float(self.d)/sregion * ratio]      if sregion
    //                   + [float(self.d)/lregion * ratio]      if lregion
    //
    // and there are **two** different template lengths in scope: a local
    // `d = self.treat.average_template_length` (the post-filter mean) and
    // `self.d = options.tsize` (the truncated as-read mean). Only `control_sum`
    // uses the local one:
    //
    //     control_sum = self.control.total * 2 * d      # the local mean
    //
    // so `ratio` is built from `avg_tl` while the wide-window factors are
    // `self.d / w`. F95 had both using `avg_tl`, which is the same number only
    // when duplicate filtering removes nothing. It is a ~0.2% error on an
    // unfiltered run -- invisible in a coordinate diff, visible in the fourth
    // decimal of `-log10(pvalue)`, which is what 21 PE fixtures showed.
    // F257: the ladder's **width** must be `self.d`, not the local `d`. F95 fixed
    // this for the *factors* (they use `cfg.tsize_exact`) but left the width on the
    // local post-filter mean, and the two differ whenever duplicate filtering removes
    // anything:
    //
    //     ctrl_d_s = [self.d]        # `self.d` = options.tsize, truncated as-read mean
    //
    // Measured on `sweep/gmini_mfrag_d400_w180_ctrl_357` variant `B`: local `d`
    // truncates to **143**, `cfg.tsize_exact` is **144.155** -> **144**. Since the
    // counted control window is `[anchor - d//2, anchor + d - d//2)`, that is
    // `half_d` 71 vs 72 -- a one-base-narrower window on every one of the two
    // anchors per fragment. The resulting control lambda is 1-3% low in a
    // position-dependent way (exact agreement wherever the fragments happen to leave
    // the extra base empty), which is what F253/F255 were chasing as a mysterious
    // "shift". Decisive arithmetic at p=20: upstream's coverage there is 734 count
    // units, ours 715, and 734 is exactly the `half_d = 72` answer while 715 is
    // exactly the `half_d = 71` answer.
    let scales: Vec<i64> = std::iter::once(cfg.tsize_exact as i64)
        .chain(std::iter::once(cfg.slocal).filter(|v| *v > 0))
        .chain(std::iter::once(cfg.llocal).filter(|v| *v > cfg.slocal && *v > 0))
        .collect();
    let ratio_f = if to_control { 1.0 } else { ratio };
    // F242: the with-control ladder's factors are `ratio`, then `self.d / w * ratio`.
    // The no-control single scale already carries its own factor, computed above.
    let factors: Vec<f32> = scales
        .iter()
        .enumerate()
        .map(|(i, w)| {
            if i == 0 {
                ratio_f as f32
            } else {
                (cfg.tsize_exact / *w as f64 * ratio_f) as f32
            }
        })
        .collect();
    // Diagnostic for F236: the scale ladder and the factors that feed it. The
    // control lambda is `max over i of (window max_i) * factor_i`, so pinning a
    // one-ulp answer needs all three factors, not just the winning one.
    if std::env::var_os("MACS3_RS_DUMP_CTRL_LAMBDA").is_some() {
        eprintln!(
            "PE lambda scales: d={} slocal={} llocal={} ratio={:?} factors={:?} tsize_exact={:?} control_sum={control_sum} treat_sum={treat_sum}",
            scales[0],
            scales.get(1).copied().unwrap_or(0),
            scales.get(2).copied().unwrap_or(0),
            ratio,
            factors,
            cfg.tsize_exact,
        );
    }
    let lscales = LambdaScales {
        d: scales[0],
        slocal: *scales.get(1).unwrap_or(&0),
        llocal: *scales.get(2).unwrap_or(&0),
        d_factor: factors[0],
        slocal_factor: factors.get(1).copied().unwrap_or(0.0),
        llocal_factor: factors.get(2).copied().unwrap_or(0.0),
    };

    // Chromosomes are the union of both inputs, sorted by name, and only those
    // the *treatment* genome can resolve are processed -- upstream iterates the
    // treatment's `get_chr_names()`.
    let mut names: Vec<Vec<u8>> = treat
        .chroms()
        .iter()
        .map(|c| treat.genome().name(*c).to_vec())
        .collect();
    if let Some(c) = ctrl {
        for cid in c.chroms() {
            let n = c.genome().name(cid).to_vec();
            if !names.contains(&n) {
                names.push(n);
            }
        }
    }
    names.sort();

    let mut signals: Vec<ChromSignals> = Vec::new();
    for name in &names {
        let Some(chrom) = treat.genome().get(name) else {
            continue;
        };
        // F90: the third argument is the track *length*, not the fragment length.
        // F151: a counted track (FRAG) needs the weighted sweep -- upstream calls
        // `pileup_from_LRC_as_list`, whose depth is the sum of counts. The
        // unweighted sweep left every FRAG fixture with a depth ~60x too small
        // and no peaks at all.
        let counted = treat.has_counts();
        let t = if counted {
            // F181: `pre_p` starts at upstream's raw zero, which is `coord_shift`
            // in the shifted frame, so the leading zero-depth run is suppressed.
            macs_pileup::pileup_from_weighted_fragments_from(
                chrom,
                &wfrags_of(treat, chrom),
                RLENGTH,
                treat_scale,
                0.0,
                coord_shift.max(0) as Coord,
            )
        } else {
            macs_pileup::pileup_from_fragments(
                chrom,
                &frags_of(treat, chrom),
                RLENGTH,
                treat_scale,
                0.0,
            )
        };
        // F162: `--nolambda` empties `ctrl_d_s`/`ctrl_scaling_factor_s`, which
        // sets `no_lambda_flag` in `CallerFromAlignments.__init__`, and the
        // control pileup then becomes the one-element `[lambda_bg]` array
        // (`CallPeakUnit.py:621-622`). That flag is set by the **scale lists
        // being empty**, not by the absence of a control file, so it applies with
        // a control file present too -- and the flag is checked *before* `ctrl`:
        //
        //     if not self.no_lambda_flag:
        //         ctrl_pv = self.ctrl.pileup_a_chromosome_c(...)
        //     else:
        //         ctrl_pv = [treat_pv[0][-1:], np.array([self.lambda_bg], dtype="f4")]
        //
        // Keeping `--nolambda` in the `None` arm only meant it was silently
        // ignored whenever `-c` was supplied, and the dynamic lambda was used
        // instead: `frag_basic/barcode_fragments` reported a summit pileup of
        // 6317 against upstream's 180283.
        let lambda = match ctrl {
            Some(_) | None if cfg.nolambda || lscales.d <= 0 => {
                // F164: upstream's substitute is the raw pair
                // `[treat_pv[0][-1:], [lambda_bg]]` -- **one** position carrying
                // `lambda_bg`. Building it as a run-length track with
                // `SignalTrack::empty(chrom, t.end(), t.end())` + `push(t.end(), ..)`
                // produces an *empty* run (`push` drops a run that does not advance
                // the cursor), so the control track came out with no runs at all
                // and the union walk had nothing to pair against.
                //
                // The faithful equivalent is a single run `[0, t.end())` carrying
                // `lambda_bg`: the cursor rule in `__chrom_pair_treat_ctrl` then
                // keeps `ic` pinned at 0 until the treatment is exhausted, so every
                // paired row gets `lambda_bg` -- which is what upstream's one-entry
                // array does, and what `frag_basic/barcode_fragments --nolambda`
                // needs to report a summit pileup of 7989.
                Some(macs_rle::SignalTrack::from_runs_exact(
                    chrom,
                    0,
                    t.end().max(1),
                    vec![macs_rle::Run::new(t.end().max(1), lambda_bg)],
                ))
            }
            None => {
                // F175: upstream has no separate no-control branch here -- it sets
                // `self.ctrl = treat` (`CallPeakUnit.py:497`) and piles *that* up
                // through `pileup_a_chromosome_c`. The projection therefore
                // depends on the track class, exactly as it does with a real
                // control file:
                //
                // * `PETrackI` (BAMPE/BEDPE, uncounted) -> `pileup_from_PN_shifted`:
                //   unit-depth ends, each shifted `d//2` both ways;
                // * `PETrackII` (`--format FRAG`, counted) ->
                //   `pileup_from_LRC_centers_as_list`: one `d`-wide **weighted**
                //   window per end, centred at `x - d//2`.
                //
                // Using the ends projection for a counted track collapsed the whole
                // local lambda onto a single run, so every `--nolambda`-less
                // no-control FRAG fixture reported the wrong -log10(p) (e.g.
                // 9634.71 against upstream's 6775.87 on
                // `sweep/gmini_mfrag_d1200_w180_noc_123`).
                let lregion = cfg.llocal.max(1);
                let factor =
                    treat.length() as f64 / (lregion as f64 * (treat.total() as f64) * 2.0);
                let tw: Vec<(Coord, Coord, f32)> = wfrags_of(treat, chrom);
                if treat.has_counts() {
                    let mut ws: Vec<f32> = Vec::with_capacity(tw.len() * 2);
                    let mut ws2: Vec<f32> = Vec::with_capacity(tw.len() * 2);
                    let mut cs: Vec<Coord> = Vec::with_capacity(tw.len() * 2);
                    let mut ce: Vec<Coord> = Vec::with_capacity(tw.len() * 2);
                    for &(l, r, w) in tw.iter() {
                        cs.push(l);
                        ws.push(w);
                        ce.push(r);
                        ws2.push(w);
                    }
                    Some(macs_pileup::pileup_from_weighted_positions(
                        chrom,
                        &cs,
                        &ce,
                        &ws,
                        &ws2,
                        &pe_ctrl_params(lregion, RLENGTH, factor as f32).with_baseline(lambda_bg),
                    ))
                } else {
                    let mut cstarts: Vec<Coord> = Vec::new();
                    let mut cends_v: Vec<Coord> = Vec::new();
                    for f in frags_of(treat, chrom) {
                        cstarts.push(f.0);
                        cends_v.push(f.1);
                    }
                    Some(macs_pileup::pileup_from_positions(
                        chrom,
                        &cstarts,
                        &cends_v,
                        &pe_ctrl_params(lregion, RLENGTH, factor as f32).with_baseline(lambda_bg),
                    ))
                }
            }
            Some(cc) => {
                let Some(cid) = cc.genome().get(name) else {
                    continue;
                };
                let cf = frags_of(cc, cid);
                // F151: a counted control projects **centred, weighted** windows.
                // `PETrackII.pileup_a_chromosome_c` (`PairedEndTrack.py:1458`)
                // calls `pileup_from_LRC_centers_as_list`, which makes one
                // `d`-wide window per fragment *end* at `x - d//2`, carrying the
                // count -- not `PETrackI`'s unit-depth `d//2` shifts. The two
                // classes therefore need different projections, not one shared
                // one, which is why the FRAG control lambda was silently wrong
                // even for the fixtures whose treatment side looked plausible.
                let cw: Vec<(Coord, Coord, f32)> = if cc.has_counts() {
                    wfrags_of(cc, cid)
                } else {
                    Vec::new()
                };
                let ctrl_counted = cc.has_counts();
                let mut combined: Option<SignalTrack<f32>> = None;
                for scale in lscales.as_pairs() {
                    let mut cstarts: Vec<Coord> = Vec::with_capacity(cf.len() * 2);
                    let mut cends_v: Vec<Coord> = Vec::with_capacity(cf.len() * 2);
                    for &(l, r) in cf.iter() {
                        cstarts.push(l);
                        cends_v.push(r);
                    }
                    // F187: `PETrackII`'s centred window is one base wider than
                    // `PETrackI`'s, so the parameter set follows the track class.
                    let mk = if ctrl_counted {
                        pe_ctrl_params_counted
                    } else {
                        pe_ctrl_params
                    };
                    let cparams =
                        mk(scale.d.max(1), RLENGTH, scale.scale_factor).with_baseline(lambda_bg);
                    let p = if ctrl_counted {
                        let mut ws: Vec<f32> = Vec::with_capacity(cw.len() * 2);
                        let mut ws2: Vec<f32> = Vec::with_capacity(cw.len() * 2);
                        for &(_, _, w) in cw.iter() {
                            ws.push(w);
                            ws2.push(w);
                        }
                        macs_pileup::pileup_from_weighted_positions(
                            cid, &cstarts, &cends_v, &ws, &ws2, &cparams,
                        )
                    } else {
                        macs_pileup::pileup_from_positions(cid, &cstarts, &cends_v, &cparams)
                    };
                    if false {
                        eprintln!(
                            "SCALE d={} n={} [{},{}] first_v={:.5}",
                            scale.d,
                            p.runs().len(),
                            p.runs().first().map(|r| r.end).unwrap_or(0),
                            p.runs().last().map(|r| r.end).unwrap_or(0),
                            p.runs().first().map(|r| r.value).unwrap_or(0.0)
                        );
                    }
                    combined = Some(match combined {
                        None => p,
                        Some(prev) => crate::over_two_pv_array_track(&prev, &p),
                    });
                }
                combined
            }
        };
        signals.push(ChromSignals {
            name: String::from_utf8_lossy(name).into_owned(),
            chrom,
            treat: t,
            ctrl: lambda,
        });
    }

    PeResult {
        signals,
        d,
        paired_boundaries: false,
        lambda_bg,
        coord_shift,
    }
}

/// `--cutoff-analysis`: peak count and total peak length as a function of the
/// p-value cutoff.
///
/// Transcribed from `CallerFromAlignments.__cal_pvalue_qvalue_table`
/// (`CallPeakUnit.py:875-1023`), which is where the `*_cutoff_analysis.txt` file is
/// written. It matters beyond the file itself: the cutoffs are **seeded into the
/// AFDR histogram with zero length** before the q-table is built (`pscore_stat` gets
/// `pscore_stat[cutoff] = 0` for any cutoff not already present, line 978), so a run
/// with `--cutoff-analysis` can produce a *different* q-table from the same data as one
/// without it. Reproducing only the file would therefore be wrong in a way that shows
/// up as a q-score difference elsewhere.
///
/// The output columns are `pscore qscore npeaks lpeaks avelpeak`, and a row is written
/// only for a cutoff that actually called at least one peak.
#[derive(Debug, Default, Clone)]
pub struct CutoffStats {
    /// Cutoff -> number of peaks called at that cutoff, summed over chromosomes.
    pub npeaks: BTreeMap<i64, u64>,
    /// Cutoff -> total peak length at that cutoff.
    pub lpeaks: BTreeMap<i64, u64>,
}

impl CutoffStats {
    pub fn new() -> Self {
        Self::default()
    }
}

/// The cutoff ladder: `0.3, 0.6, ... 9.9`, rounded to five decimals.
///
/// `np.arange(0.3, 10.0, 0.3)` yields 33 values. The descending sort matters because
/// `tmplist` is iterated in that order both for the per-cutoff scan and for the
/// histogram seeds.
///
/// **Known approximation.** Upstream's `tmplist` holds Python `float`s (f64) while the
/// p-scores in the same `pscore_stat` dict are `float32`, so upstream ends up with
/// *both* a 0.3-f64 key and a 0.3-f32 key whenever a p-score happens to be exactly
/// 0.3, each with its own q. This port's histogram is keyed on the f32 bit pattern
/// throughout, so the ladder collapses onto the same key a p-score of 0.3 would use and
/// the two entries become one. The reported `qscore` therefore reads the f32-keyed
/// value. That is the only way to express it without making the histogram
/// representationally f64, and it is recorded here rather than left as a surprise.
pub fn cutoff_ladder() -> Vec<f32> {
    // `round(x, 5)` is round-to-**five decimals**, not round-to-integer. Getting that
    // wrong collapses the ladder onto {0, 1, ..., 10}, which still looks like a
    // plausible 33-value list and silently compares the wrong cutoffs.
    let mut v: Vec<f32> = (1..=33)
        .map(|i| (((i as f64 * 0.3) * 1e5).round() / 1e5) as f32)
        .collect();
    v.sort_by(|a, b| b.total_cmp(a));
    v
}

/// Fold one chromosome's p-score track into the per-cutoff totals.
///
/// `track` is the p-score track: each run's `end` is a position in the same array
/// upstream calls `pos_array`, and its value the matching `score_array` entry. So
/// "positions above the cutoff" is `score_array > cutoff` and the chunk boundaries are
/// `pos_array[i]` / `pos_array[i-1]`.
///
/// The one subtlety is `above_cutoff[0] == 0`. Upstream indexes
/// `pos_array[above_cutoff - 1]`, and with NumPy that **wraps to the last element**,
/// so a peak whose first above-cutoff position is the very start of the chromosome
/// begins at the previous run's end -- which does not exist, so NumPy returns the
/// final position. That is a genuine upstream quirk, and reproducing it is what keeps
/// `lpeaks` identical rather than plausibly close.
pub fn accumulate_cutoffs(
    stats: &mut CutoffStats,
    track: &macs_rle::SignalTrack<f32>,
    ladder: &[f32],
    max_gap: Coord,
    min_length: Coord,
) {
    let runs = track.runs();
    if runs.is_empty() {
        return;
    }
    let ends: Vec<Coord> = runs.iter().map(|r| r.end).collect();
    let scores: Vec<f32> = runs.iter().map(|r| r.value).collect();

    for &cutoff in ladder {
        let cut_key = macs_score::canonical_score_key(cutoff) as i64;
        // Indices whose score is strictly above the cutoff. Upstream uses
        // `np.nonzero(score_array > cutoff)`, i.e. strictly greater.
        let above: Vec<usize> = (0..scores.len()).filter(|&i| scores[i] > cutoff).collect();
        if above.is_empty() {
            continue;
        }

        let start_at = |k: usize| -> Coord {
            // `pos_array[above_cutoff - 1]`, wrapping like NumPy.
            if above[k] == 0 {
                ends[ends.len() - 1]
            } else {
                ends[above[k] - 1]
            }
        };

        let mut total_l: u64 = 0;
        let mut total_p: u64 = 0;
        let mut chunk_start = start_at(0);
        let mut chunk_end = ends[above[0]];
        let mut last_end = chunk_end;

        for k in 1..above.len() {
            let s = start_at(k);
            let e = ends[above[k]];
            // `tl = acs_ptr[0] - lastp` in C `long` arithmetic. The subtraction is on
            // signed coordinates, so a start *before* the previous end is a negative
            // gap, and `tl <= max_gap` is then true -- which merges the chunks, exactly
            // as upstream does. Using `saturating_sub` would give 0 and also merge, so
            // the two agree here; the explicit form is kept because it is the rule, not
            // the coincidence.
            let tl = s as i64 - last_end as i64;
            if tl <= max_gap as i64 {
                chunk_end = e;
            } else {
                let len = chunk_end.saturating_sub(chunk_start);
                if len >= min_length {
                    total_l += len;
                    total_p += 1;
                }
                chunk_start = s;
                chunk_end = e;
            }
            last_end = e;
        }
        // the trailing chunk is always closed, even when it is the only one
        if !above.is_empty() {
            let len = chunk_end.saturating_sub(chunk_start);
            if len >= min_length {
                total_l += len;
                total_p += 1;
            }
        }

        *stats.npeaks.entry(cut_key).or_insert(0) += total_p;
        *stats.lpeaks.entry(cut_key).or_insert(0) += total_l;
    }
}

/// Seed the ladder cutoffs into the AFDR histogram, as upstream does before the q walk.
///
/// `pscore_stat[cutoff] = 0` for cutoffs not already present. A zero-length bucket does
/// not change `N`, but it *does* add a key to `unique_values`, so the rank `k` walk
/// visits it and the q-table acquires an entry -- which is exactly what the report's
/// `qscore` column reads back.
pub fn seed_cutoffs(hist: &mut macs_score::PScoreHistogram, ladder: &[f32]) {
    for &c in ladder {
        hist.add_zero(c);
    }
}

/// Render the `*_cutoff_analysis.txt` body.
///
/// Column formats are `%.2f` for the two scores and `%.2f` for `avelpeak`, and a row is
/// emitted only where `npeaks > 0` -- so a cutoff that called nothing is skipped
/// entirely rather than written as a zero row.
pub fn render_cutoff_analysis(
    stats: &CutoffStats,
    ladder: &[f32],
    table: &macs_score::PqTable,
) -> String {
    let mut out = String::from("pscore\tqscore\tnpeaks\tlpeaks\tavelpeak\n");
    for &cutoff in ladder {
        let key = macs_score::canonical_score_key(cutoff) as i64;
        let np = stats.npeaks.get(&key).copied().unwrap_or(0);
        if np == 0 {
            continue;
        }
        let lp = stats.lpeaks.get(&key).copied().unwrap_or(0);
        let q = table.qscore_or_zero(cutoff);
        out.push_str(&format!(
            "{cutoff:.2}\t{q:.2}\t{np}\t{lp}\t{:.2}\n",
            lp as f64 / np as f64
        ));
    }
    out
}

/// `bedGraphTrackI.cutoff_analysis` (`BedGraph.py:1262-1377`): peak count and total
/// length as a function of a bedGraph score cutoff.
///
/// Same chunk walk as [`accumulate_cutoffs`], but over a *score* sweep derived from the
/// data rather than a fixed p-value ladder, and with one fewer column (there is no
/// q-value for a bedGraph). The output is descending by cutoff.
pub fn bedgraph_cutoff_analysis(
    tracks: &[(macs_rle::SignalTrack<f32>, Coord)],
    max_gap: Coord,
    min_length: Coord,
    steps: usize,
    min_score: f32,
    max_score: f32,
) -> String {
    // `minv = max(min_score, self.minvalue)` / `maxv = min(self.maxvalue, max_score)`,
    // where `self.minvalue` / `self.maxvalue` are the track's own extremes.
    let (mut minv, mut maxv) = (f32::MAX, f32::MIN);
    for (t, _) in tracks {
        for r in t.runs() {
            minv = minv.min(r.value);
            maxv = maxv.max(r.value);
        }
    }
    if !minv.is_finite() || !maxv.is_finite() {
        minv = 0.0;
        maxv = 0.0;
    }
    let minv = minv.max(min_score);
    let maxv = maxv.min(max_score);
    // `s = float(maxv - minv) / steps`. A zero-width range makes upstream's
    // `np.arange(minv, maxv, 0)` raise; there is no cutoff to report, so the report is
    // header-only rather than a crash or a division by zero.
    let s = (f64::from(maxv) - f64::from(minv)) / steps as f64;
    // `!(s > 0.0)` rather than `s <= 0.0`, so a NaN range (a track whose values are
    // not finite) is rejected too: upstream's `np.arange(minv, maxv, 0)` raises there,
    // and the honest outcome is a header-only report rather than a division by zero or
    // a silently empty sweep.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    if !(s > 0.0) {
        return "score\tnpeaks\tlpeaks\tavelpeak\n".to_string();
    }
    // `np.arange(minv, maxv, s)` with `round(v, 3)`. `arange` is exclusive of `maxv`,
    // so the count is `ceil((maxv - minv) / s)`, which is `steps` for a range that
    // divides evenly.
    let n = (((f64::from(maxv) - f64::from(minv)) / s).ceil() as usize).max(1);
    let ladder: Vec<f32> = (0..n)
        .map(|i| ((f64::from(minv) + i as f64 * s) * 1e3).round() / 1e3)
        .map(|v| v as f32)
        .collect();

    let mut npeaks = vec![0u64; ladder.len()];
    let mut lpeaks = vec![0u64; ladder.len()];
    for (track, _) in tracks {
        let runs = track.runs();
        if runs.is_empty() {
            continue;
        }
        let ends: Vec<Coord> = runs.iter().map(|r| r.end).collect();
        let scores: Vec<f32> = runs.iter().map(|r| r.value).collect();
        for (n_idx, &cutoff) in ladder.iter().enumerate() {
            let above: Vec<usize> = (0..scores.len()).filter(|&i| scores[i] > cutoff).collect();
            if above.is_empty() {
                continue;
            }
            let start_at = |k: usize| {
                if above[k] == 0 {
                    ends[ends.len() - 1]
                } else {
                    ends[above[k] - 1]
                }
            };
            let (mut tl, mut tp) = (0u64, 0u64);
            let (mut cs, mut ce) = (start_at(0), ends[above[0]]);
            let mut last_end = ce;
            for k in 1..above.len() {
                let (s, e) = (start_at(k), ends[above[k]]);
                if (s as i64 - last_end as i64) <= max_gap as i64 {
                    ce = e;
                } else {
                    let len = ce.saturating_sub(cs);
                    if len >= min_length {
                        tl += len;
                        tp += 1;
                    }
                    cs = s;
                    ce = e;
                }
                last_end = e;
            }
            let len = ce.saturating_sub(cs);
            if len >= min_length {
                tl += len;
                tp += 1;
            }
            lpeaks[n_idx] += tl;
            npeaks[n_idx] += tp;
        }
    }

    let mut out = String::from("score\tnpeaks\tlpeaks\tavelpeak\n");
    for n in (0..ladder.len()).rev() {
        if npeaks[n] == 0 {
            continue;
        }
        out.push_str(&format!(
            "{:.2}\t{}\t{}\t{:.2}\n",
            ladder[n],
            npeaks[n],
            lpeaks[n],
            lpeaks[n] as f64 / npeaks[n] as f64
        ));
    }
    out
}
