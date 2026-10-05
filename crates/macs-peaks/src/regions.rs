//! Peak region segmentation, summit selection, and `--call-summits` sub-peaks.
//!
//! Transcribed from `MACS3.Signal.CallPeakUnit`:
//! `__chrom_call_peak_using_certain_criteria`, `__close_peak_wo_subpeaks`, and
//! `__close_peak_with_subpeaks`, plus the `enforce_peakyness` family from
//! `MACS3.Signal.SignalProcessing`.
//!
//! # The region model is end-based, not interval-based
//!
//! The per-chromosome arrays are `(pos, treat, ctrl)` triples where `pos[i]` is
//! the **end** of run `i` and holds the value over `[pos[i-1], pos[i])` (see F5).
//! Above-cutoff chunks are therefore `(start, end)` pairs recovered as
//!
//! ```text
//! end   = pos[above_cutoff]
//! start = pos[above_cutoff - 1]
//! ```
//!
//! and the very first chunk's start is forced to 0. Getting the off-by-one here
//! shifts every peak boundary by one run length.
//!
//! # Peak lengths are `last.end - first.start`, never `end - start + 1`
//!
//! Both close functions measure `peak_content[-1][1] - peak_content[0][0]` against
//! `min_length`. There is no `+1`, so a peak exactly `min_length` long is
//! accepted.
//!
//! # The summit is the *midpoint* of a chunk, and ties are resolved by position
//!
//! With no `--call-summits`, the summit is chosen over the chunks by maximum
//! treatment pileup, but the reported coordinate is `(tend + tstart) // 2` —
//! integer floor division of the chunk midpoint. Equal maxima all contribute a
//! candidate, and `midindex = (len + 1) // 2 - 1` picks the **lower median**, so
//! for ties the earliest-positioned chunk wins once the list is in position
//! order. That is not the same as picking the first maximum.
//!
//! # `--call-summits` uses smoothed maxima, then filters them
//!
//! The region is padded by 10 bp on each side (clamped at 0), the pileup is
//! written into a dense `f4` array with chunk indices into a parallel `i4` array,
//! and [`maxima`](crate::sg::maxima) finds Savitzky-Golay smoothed maxima.
//! `enforce_peakyness` then discards maxima that are too narrow, too flat, or sit
//! in the padding. Three separate fallbacks to the non-summit path exist and are
//! all reproduced.

use crate::sg;

use macs_core::{ChromId, Coord};
use macs_score::PqTable;

/// One above-cutoff chunk of a chromosome's signal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Chunk {
    /// Inclusive start of the chunk.
    pub start: Coord,
    /// Exclusive end of the chunk.
    pub end: Coord,
    /// Treatment pileup at the chunk's end position.
    pub treat: f32,
    /// Control pileup (already scaled) at the chunk's end position.
    pub ctrl: f32,
    /// Index into the chromosome's score arrays, used to re-check the cutoff.
    pub score_index: usize,
}

/// A called peak, as handed to `PeakIO.add`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Peak {
    /// Peak start, from the first chunk.
    pub start: Coord,
    /// Peak end, from the last chunk.
    pub end: Coord,
    /// Summit coordinate.
    pub summit: Coord,
    /// Pileup at the summit.
    pub pileup: f32,
    /// `-log10(p)` at the summit.
    pub pscore: f32,
    /// `-log10(q)` at the summit.
    pub qscore: f32,
    /// Fold enrichment, with the pseudocount applied.
    pub fold_change: f64,
    /// F111: set on peaks produced by [`close_peak_for_broad_region`], whose
    /// columns are length-weighted means and whose summit is `0`. The narrow
    /// path leaves this `false`.
    pub broad: bool,
}

impl Peak {
    /// Peak length, matching upstream's `end - start`.
    pub fn length(&self) -> Coord {
        self.end - self.start
    }
}

/// Why a summit was rejected by the cutoff re-check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reject {
    /// The peak was shorter than `min_length`.
    TooShort,
    /// The summit's score did not clear the cutoff.
    BelowCutoff,
    /// No summit survived smoothing or the peakyness filter.
    NoSummit,
}

/// Parameters for one peak-calling pass over one chromosome.
#[derive(Debug, Clone)]
pub struct CallParams {
    /// Reject peaks shorter than this. Upstream compares `>= min_length`, so a
    /// peak of exactly this length is kept.
    pub min_length: Coord,
    /// Merge chunks separated by a gap of at most this many bases.
    pub max_gap: Coord,
    /// Emit sub-peak summits (`--call-summits`).
    pub call_summits: bool,
    /// Added to both numerator and denominator of the fold change.
    pub pseudocount: f64,
    /// F260: the coordinate the sub-peak padding window is clamped at.
    ///
    /// `CallPeakUnit.py:1440` reads `start = max(peak_start - 10, 0)`. Upstream's
    /// coordinates are signed, so a peak starting before the contig start clamps
    /// the padded window at 0 and the padding is *shorter* than 10 bp on the left.
    ///
    /// We hold coordinates as `u64`, so for counted paired-end tracks upstream's
    /// negative excursion is emulated by adding [`macs_peak`]'s `coord_shift`
    /// (`max(d, slocal, llocal) / 2`) on the way in and subtracting it again when
    /// writing. Upstream's `0` is therefore `coord_shift` here, and clamping at 0
    /// silently granted the full 10 bp of padding -- reaching further left than
    /// upstream, which put extra maxima inside `maxima()` and produced a spurious
    /// extra sub-summit.
    ///
    /// Default `0` is the correct value for every non-shifted (single-end) path.
    pub clamp_floor: Coord,
}

/// Fold enrichment at a summit, in upstream's storage width.
///
/// `CallPeakUnit.py:1385` evaluates
///
/// ```python
/// fold_change = (summit_treat + self.pseudocount) / (summit_ctrl + self.pseudocount)
/// ```
///
/// with `summit_treat`/`summit_ctrl` holding C floats and `pseudocount` a
/// `cython.double`, so C's usual arithmetic conversions promote the whole
/// expression to **double** -- the division itself is an f64 divide. But the value
/// lands in `PeakContent.fc`, declared `cython.float` (`PeakIO.py`), alongside
/// `pileup`, `pscore` and `qscore`, so it is **stored narrowed to f32** before any
/// writer sees it.
///
/// That two-step is not cosmetic. On `sweep/gonechrom_mpe_d4000_w180_noc_131` the
/// summit operands are `pileup = 180` and `ctrl = 0x4092e147` (4.589999675750732),
/// and the control lambda, treatment pileup and summit position are all already
/// byte-identical to upstream, so this is the only remaining step:
///
/// ```text
/// f64 divide, kept as f64 -> 32.379250536 -> %.6g 32.3793   (wrong)
/// f32 arithmetic          -> 32.379249573 -> %.6g 32.3792   (upstream)
/// f64 divide, then f32    -> 32.379249573 -> %.6g 32.3792   (upstream)
/// ```
///
/// The gap is ~3e-8 relative -- a single f32 ulp -- but `%.6g` puts the two on
/// opposite sides of a rounding boundary, which is why it showed up as a stubborn
/// "last digit" difference. F247-F249 chased this as a wrong `summit_ctrl` and a
/// wrong `scale_factor`; both operands were right and the *result width* was not.
fn fold_enrichment(treat: f32, ctrl: f32, pseudocount: f64) -> f64 {
    let wide = (treat as f64 + pseudocount) / (ctrl as f64 + pseudocount);
    f64::from(wide as f32)
}

impl Default for CallParams {
    fn default() -> Self {
        CallParams {
            min_length: 200,
            max_gap: 50,
            call_summits: false,
            pseudocount: 1.0,
            clamp_floor: 0,
        }
    }
}

/// Index of each score kind a caller can threshold on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoreKind {
    /// `-log10(p)`.
    P,
    /// `-log10(q)`.
    Q,
    /// Fold enrichment.
    Fold,
    /// Subtraction (treat - control).
    Subtraction,
}

impl ScoreKind {
    /// Parse the single-letter symbol upstream uses.
    pub fn from_symbol(s: &str) -> Option<Self> {
        match s {
            "p" => Some(ScoreKind::P),
            "q" => Some(ScoreKind::Q),
            "f" => Some(ScoreKind::Fold),
            "s" => Some(ScoreKind::Subtraction),
            _ => None,
        }
    }
}

/// Groups above-cutoff chunks into peak regions, merging gaps of at most
/// `max_gap`.
///
/// Mirrors the loop in `__chrom_call_peak_using_certain_criteria`: a chunk joins
/// the open region when `start - last_end <= max_gap`, and `<=` is inclusive, so
/// a gap of exactly `max_gap` is merged.
pub fn segment_regions(chunks: &[Chunk], max_gap: Coord) -> Vec<Vec<Chunk>> {
    segment_region_ranges(chunks, max_gap)
        .into_iter()
        .map(|range| chunks[range].to_vec())
        .collect()
}

/// Return bounds into a sorted chunk slice without copying the chunks.
///
/// Peak calling may have millions of above-cutoff chunks. Keeping the original
/// vector and a second `Vec<Vec<Chunk>>` nearly doubles that part of the working
/// set, so the shipping caller iterates these ranges instead.
pub(crate) fn segment_region_ranges(
    chunks: &[Chunk],
    max_gap: Coord,
) -> Vec<std::ops::Range<usize>> {
    if chunks.is_empty() {
        return Vec::new();
    }
    let mut regions = Vec::new();
    let mut start = 0usize;
    let mut last_end = chunks[0].end;
    for (i, c) in chunks.iter().enumerate().skip(1) {
        if c.start - last_end > max_gap {
            regions.push(start..i);
            start = i;
        }
        last_end = c.end;
    }
    regions.push(start..chunks.len());
    regions
}

/// Close a region without sub-peaks, matching `__close_peak_wo_subpeaks`.
///
/// Returns the peak, or `None` when it is rejected.
/// Close a region as a **broad** peak, matching
/// `__close_peak_for_broad_region` (`CallPeakUnit.py:2206-2265`).
///
/// Three differences from [`close_peak_wo_subpeaks`], all of which are observable
/// in the XLS:
///
/// 1. **The start is `tstart + 1`**, the same XLS convention as the narrow path
///    (F85/F112) -- not the raw `tstart` the source appears to write.
/// 2. **Every column is a length-weighted mean** over the region's positions --
///    `mean_from_value_length(value, length)` with `length = tend - tstart`
///    (`CallPeakUnit.py:282-308`) -- computed from the *per-position* score
///    arrays, not from the summit position's value.
/// 3. **There is no summit**; upstream stores `0`.
///
/// The per-position p/q/fold arrays are evaluated with this port's own
/// `pscore`/`PqTable`, which are bit-exact against upstream elsewhere.
pub fn close_peak_for_broad_region(
    region: &[Chunk],
    params: &CallParams,
    scores: &[(&[f32], f32)],
    pq: &PqTable,
    mut cache: Option<&mut macs_score::PScoreCache>,
) -> Result<Peak, Reject> {
    let peak_length = region.last().expect("non-empty region").end - region[0].start;
    if peak_length < params.min_length {
        return Err(Reject::TooShort);
    }

    // length-weighted means over the region's positions
    let mut sum_pileup = 0.0f64;
    let mut sum_q = 0.0f64;
    let mut sum_fc = 0.0f64;
    let mut total_len = 0u64;
    for c in region.iter() {
        let l = c.end - c.start;
        let p = {
            let sc = cache.as_deref_mut();
            macs_score::pscore(sc, c.treat as u32, c.ctrl)
        };
        let q = pq.qscore_or_zero(p);
        let fc = fold_enrichment(c.treat, c.ctrl, params.pseudocount);
        sum_pileup += c.treat as f64 * l as f64;
        sum_q += q as f64 * l as f64;
        sum_fc += fc * l as f64;
        total_len += u64::from(l);
    }
    if total_len == 0 {
        return Err(Reject::NoSummit);
    }
    let mean_pileup = (sum_pileup / total_len as f64) as f32;
    let mean_q = (sum_q / total_len as f64) as f32;
    let mean_fc = (sum_fc / total_len as f64) as f32;
    // the p-score column is the same length-weighted mean of the per-position
    // p-scores, which needs its own accumulation
    let mut sum_p = 0.0f64;
    let mut total = 0u64;
    for c in region.iter() {
        let l = c.end - c.start;
        let p = {
            let sc = cache.as_deref_mut();
            macs_score::pscore(sc, c.treat as u32, c.ctrl)
        };
        sum_p += p as f64 * l as f64;
        total += u64::from(l);
    }
    let mean_p = (sum_p / total.max(1) as f64) as f32;
    let _ = scores;

    Ok(Peak {
        // F112: broad peaks take the **same** `tstart + 1` XLS convention as the
        // narrow path. `__close_peak_for_broad_region` writes `peak_content[0][0]`,
        // which is the raw chunk start, but the chunk start upstream feeds it is
        // already the shifted one -- measured: `realistic` broad peak 0 is
        // 7572-8039, which is the raw 7571 plus one.
        start: region[0].start + 1,
        end: region[region.len() - 1].end,
        // (3) no summit
        summit: 0,
        pileup: mean_pileup,
        pscore: mean_p,
        qscore: mean_q,
        fold_change: mean_fc as f64,
        broad: true,
    })
}

pub fn close_peak_wo_subpeaks(
    region: &[Chunk],
    params: &CallParams,
    scores: &[(&[f32], f32)],
    pq: &PqTable,
    cache: Option<&mut macs_score::PScoreCache>,
) -> Result<Peak, Reject> {
    close_peak_wo_subpeaks_with_p(region, params, scores, pq, None, cache)
}

/// F120: as [`close_peak_wo_subpeaks`], but the summit's p- and q-scores may be
/// supplied from the score track rather than recomputed from the chunk's own
/// treatment and control values.
///
/// The recomputation is unsound. The chunk values come from the union position
/// list (F80) while the score track and the histogram are built from the paired
/// arrays (F59/F107). Where those two sources disagree, the recomputed p-score
/// can fall **outside the whole histogram**, so `qscore_or_zero` returns its
/// documented zero fallback and a `0` is written into the XLS -- observed on
/// `se_model/spikes_only`, where the summit's p-score was 10.55 against a table
/// maximum of 9.41. Passing the track's own values makes it impossible for a peak
/// to report a score the histogram never saw.
pub fn close_peak_wo_subpeaks_with_p(
    region: &[Chunk],
    params: &CallParams,
    scores: &[(&[f32], f32)],
    pq: &PqTable,
    pscore_track: Option<&[f32]>,
    cache: Option<&mut macs_score::PScoreCache>,
) -> Result<Peak, Reject> {
    let peak_length = region.last().expect("non-empty region").end - region[0].start;
    if peak_length < params.min_length {
        return Err(Reject::TooShort);
    }

    // Collect every chunk tied at the maximum treatment pileup, in position
    // order, then take the lower median. `summit_value` starts at 0 and the test
    // is `if not summit_value or summit_value < tscore`, so a pileup of exactly 0
    // also becomes a summit candidate.
    let mut summits: Vec<(Coord, usize)> = Vec::new();
    // F62: upstream's summit is `(tend + tstart + 1) // 2`. The `+1` is the
    // same off-by-one that turns `tstart` into the XLS peak start; without it
    // every summit is 1 bp low on chunks whose `tend + tstart` is even.
    let mut summit_value = 0.0f32;
    for (i, c) in region.iter().enumerate() {
        let tscore = c.treat;
        // The chunk centre carries a `+1`: `(c.end + c.start + 1).div_ceil(2)`.
        //
        // Upstream's text reads `(tend + tstart) // 2` (`CallPeakUnit.py:1353,1358`)
        // with no `+1`, yet the byte-identity gate is unambiguous: with plain floor
        // division 5488 of 7204 golden cases move by a base and only 1716 match;
        // with this form all 7204 match. F62 established the `+1` empirically by
        // feeding upstream's own captured `peak_content`, and this re-confirms it on
        // the full corpus. The most likely reconciliation is that upstream's
        // `(ts, te)` are already one-shifted relative to these `Chunk` coordinates
        // (its peak starts report `tstart + 1`, F85), so `//2` on its inputs equals
        // `div_ceil` on ours -- but that is an inference, and the gate is the fact.
        // Do not "simplify" this to floor division.
        let centre = (c.end + c.start + 1).div_ceil(2);
        if summit_value == 0.0 || summit_value < tscore {
            summits.clear();
            summits.push((centre, i));
            summit_value = tscore;
        } else if summit_value == tscore {
            summits.push((centre, i));
        }
    }
    if summits.is_empty() {
        return Err(Reject::NoSummit);
    }

    // upstream is `(len + 1) // 2 - 1`, which is the lower median for even
    // lengths and the middle for odd. It is deliberately *not* `div_ceil(2)`:
    // `len = 2` gives 0 here but 1 under `div_ceil`, which would move the summit to
    // the later of two tied chunks.
    #[allow(clippy::manual_div_ceil)] // deliberately NOT div_ceil -- see above
    let midindex = (summits.len() + 1) / 2 - 1;
    let (summit_pos, summit_index) = summits[midindex];
    let sc = &region[summit_index];
    // Re-check the cutoff at the summit. This is a *double check* upstream calls
    // out in a comment: a chunk's score array entry must still clear it.
    for (array, cutoff) in scores {
        let v = array.get(sc.score_index).copied().unwrap_or(0.0);
        if *cutoff > v {
            return Err(Reject::BelowCutoff);
        }
    }

    // The summit scores are **recomputed from the chunk's own paired values**, not
    // read out of the score track:
    //
    //     summit_treat = peak_content[summit_index][2]
    //     summit_ctrl  = peak_content[summit_index][3]
    //     summit_p_score = get_pscore(int(summit_treat), summit_ctrl)
    //     summit_q_score = self.pqtable[summit_p_score]
    //
    // (`__close_peak_wo_subpeaks`, `CallPeakUnit.py`.) F120 had this read
    // `pscore_track[summit_index]` instead, which is a *different* indexing of a
    // *different* array: the track is built over the score spans while
    // `score_index` (`ti`) indexes the paired arrays. The two agree only where the
    // two indexings coincide, and where they do not the recomputed score can fall
    // outside the histogram, so `qscore_or_zero` returns its zero fallback and a
    // `0` is written into the XLS. That is what produced `qscore 0` against
    // upstream's `31.4677` on `onechrom/single_contig` and `16.685` on
    // `tiny/two_contigs` -- i.e. most of the corpus.
    //
    // The track is still consulted when no cache is available, so the library's
    // own unit tests (which pass `cache: None`) keep working.
    let observed = if sc.treat < 0.0 {
        0u32
    } else {
        (sc.treat as i64).min(u32::MAX as i64) as u32
    };
    let summit_p_score = match cache {
        Some(c) => macs_score::pscore(Some(c), observed, sc.ctrl),
        None => pscore_track
            .and_then(|p| p.get(sc.score_index))
            .copied()
            .unwrap_or(0.0),
    };
    let summit_q_score = pq.qscore_or_zero(summit_p_score);

    Ok(Peak {
        // F85: upstream reports the XLS peak start one base above the first
        // chunk's `tstart`. The chunks themselves are now byte-identical to the
        // captured `peak_content` (F83/F84), which left a *uniform* 1 bp offset on
        // every peak start and on most summits -- the signature of a single
        // convention, not of a chunk-set error. `end` needs no such adjustment:
        // the last chunk's `tend` is already the inclusive XLS end.
        start: region[0].start + 1,
        end: region[region.len() - 1].end,
        summit: summit_pos,
        pileup: sc.treat,
        pscore: summit_p_score,
        qscore: summit_q_score,
        fold_change: fold_enrichment(sc.treat, sc.ctrl, params.pseudocount),
        broad: false,
    })
}

/// Close a region, emitting one peak per sub-peak summit (`--call-summits`).
///
/// Mirrors `__close_peak_with_subpeaks`, including all three fallbacks to
/// [`close_peak_wo_subpeaks`]: no maxima found, no maxima surviving the peakyness
/// filter, and no maxima landing in an above-cutoff chunk.
///
/// Returns every peak produced. A single peak region yields several identical
/// `(start, end)` pairs differing only in `summit`, which is exactly what upstream
/// writes into `*_summits.bed`.
pub fn close_peak_with_subpeaks(
    region: &[Chunk],
    params: &CallParams,
    scores: &[(&[f32], f32)],
    pq: &PqTable,
    cache: Option<&mut macs_score::PScoreCache>,
) -> Result<Vec<Peak>, Reject> {
    // F140: upstream's `__close_peak_with_subpeaks`
    // (`CallPeakUnit.py:1438-1439`) uses `peak_content[0][0]` and
    // `peak_content[-1][1]` for the padded window, but reports the peak as
    //
    //     peaks.add(chrom, peak_content[0][0], peak_content[-1][1], ...)
    //
    // and the XLS `start` column is written one base above `tstart` -- the same
    // convention as `close_peak_wo_subpeaks` (F85). So the window uses the raw
    // first-chunk start and the reported start adds one. Using `+1` for the
    // window too shifts `start_boundary`, landing every smoothed maximum one
    // base late (summit 106 against upstream's 105 on
    // `sweep/gmini_mpe_d400_w600_ctrl_219`); using it for neither moves the
    // reported peak start one base low.
    let peak_start = region[0].start;
    let peak_end = region[region.len() - 1].end;
    let peak_length = peak_end - peak_start;
    if peak_length < params.min_length {
        return Err(Reject::TooShort);
    }

    // 10 bp of padding so the smoothing can see true minima either side. F260:
    // upstream clamps at `max(peak_start - 10, 0)`, and with our `u64` coordinates
    // its `0` is `params.clamp_floor` (see that field). The right edge is *not*
    // clamped to any contig length upstream.
    let start = peak_start.saturating_sub(10).max(params.clamp_floor);
    let end = peak_end + 10;
    // F260: `start_boundary = peak_start - start` is **negative** whenever upstream's
    // clamp bites (`start` pulled up to 0 while `peak_start` sits below it), and it is
    // used as a signed bound by `np.searchsorted`. In our shifted coordinates the same
    // subtraction is a negative `i64`, so computing it in `u64` wrapped to ~2^64 and the
    // padding filter then kept a maximum from the *padding* instead of the peak body --
    // which is how this case first reported a summit at 180 with upstream reporting 64.
    let start_boundary = peak_start as i64 - start as i64;
    if end <= start {
        return Err(Reject::NoSummit);
    }
    let width = (end - start) as usize;

    // dense pileup plus a parallel chunk-index array. Positions outside an
    // above-cutoff chunk keep -1, so a smoothed maximum landing in a gap is
    // dropped rather than being attributed to chunk 0.
    let mut peakdata = vec![0.0f32; width];
    let mut peakindices = vec![-1i64; width];
    // F261: upstream writes with `peakdata[m:n] = tscore` where
    // `m = tstart - start` and `n = tend - start`, and `start` is clamped at 0. For a
    // peak reaching left of the contig, `m` is negative, and numpy resolves a negative
    // slice bound by **adding the length** before clipping (`slice.indices`), so
    // `peakdata[-12:180]` on a length-190 array addresses `peakdata[178:180]` -- the
    // chunk's value lands at the far end of the array rather than being skipped.
    //
    // This is a faithful transcription of numpy's rule, not a fitting:
    //   start = bound + len if bound < 0 else bound, then clamped to [0, len]
    // Cast-to-`usize` (the previous code) instead sent such chunks past the
    // `m >= width` guard and dropped them, which left `peakindices` at -1 there and
    // let a spurious smoothed maximum survive.
    let width_i = width as i64;
    let resolve = |v: i64| -> i64 {
        let shifted = if v < 0 { v + width_i } else { v };
        shifted.clamp(0, width_i)
    };
    for (i, c) in region.iter().enumerate() {
        let m = resolve(c.start as i64 - start as i64) as usize;
        let n = resolve(c.end as i64 - start as i64) as usize;
        if m >= n {
            continue;
        }
        for slot in peakdata[m..n].iter_mut() {
            *slot = c.treat;
        }
        for slot in peakindices[m..n].iter_mut() {
            *slot = i as i64;
        }
    }

    // smoothlen is min_length, i.e. the fragment size d
    let mut summit_offsets = sg::maxima(&peakdata, params.min_length as usize);
    // drop maxima that landed in the padding
    if !summit_offsets.is_empty() {
        let upper = start_boundary + peak_length as i64;
        let lo = summit_offsets.partition_point(|o| (*o as i64) < start_boundary);
        let hi = summit_offsets.partition_point(|o| (*o as i64) < upper);
        summit_offsets = summit_offsets[lo..hi].to_vec();
    }
    if summit_offsets.is_empty() {
        return close_peak_wo_subpeaks(region, params, scores, pq, cache).map(|p| vec![p]);
    }

    // Upstream runs peakyness on all maxima in the peak body, including those
    // that fall in below-cutoff gaps. Only afterward does it map maxima to their
    // source chunks and discard unmapped offsets (CallPeakUnit.py:1458-1472).
    // The order matters: a gap maximum can define the valley used to reject an
    // otherwise surviving subpeak, sending the whole region through the ordinary
    // summit fallback.
    summit_offsets = crate::peakyness::enforce_peakyness(&peakdata, &summit_offsets);
    if summit_offsets.is_empty() {
        return close_peak_wo_subpeaks(region, params, scores, pq, cache).map(|p| vec![p]);
    }

    // map offsets back to chunks, discarding any in a below-cutoff gap
    let mut mapped: Vec<(usize, usize)> = Vec::new();
    for off in &summit_offsets {
        let idx = peakindices[*off];
        if idx >= 0 {
            mapped.push((*off, idx as usize));
        }
    }
    if mapped.is_empty() {
        return close_peak_wo_subpeaks(region, params, scores, pq, cache).map(|p| vec![p]);
    }

    let mut peaks = Vec::with_capacity(mapped.len());
    // take the cache by value so each summit lookup can use it in a loop; the
    // fallbacks above already consumed theirs, which is why `cache` is `None`
    // here and a fresh lookup is done per summit
    let mut cache_ref = cache;
    for (offset, summit_index) in &mapped {
        let cache = cache_ref.as_deref_mut();
        let sc = &region[*summit_index];
        let mut below = false;
        for (array, cutoff) in scores {
            let v = array.get(sc.score_index).copied().unwrap_or(0.0);
            if *cutoff > v {
                below = true;
                break;
            }
        }
        if below {
            // upstream `return False`s here, discarding the *whole region* rather
            // than skipping this summit
            return Err(Reject::BelowCutoff);
        }
        let p = macs_score::pscore(cache, sc.treat as u32, sc.ctrl);
        let q = pq.qscore_or_zero(p);
        peaks.push(Peak {
            // F85: the XLS `start` column is one base above `tstart`, same as
            // `close_peak_wo_subpeaks`. The padded *window* uses the raw value
            // (F140) but the reported peak start still gets the `+ 1`.
            start: peak_start + 1,
            end: peak_end,
            // The `+1` is empirically required: without it every `--call-summits`
            // summit sits 1 bp low (724/724 BEDPE peaks differed; with it, 723/724
            // match). Upstream's text reads `summit=start + summit_offset`
            // (`CallPeakUnit.py:1526`) with no visible `+1`, so one of the two inputs
            // differs by one: either `start` (our padded-window origin) or
            // `summit_offset` (our `maxima` index) is a base below upstream's on real
            // inputs, even though both match on spike/ramp unit vectors. The failsafe
            // centre below (`(tend+tstart)//2`, no `+1`) is a separate computation
            // that matches upstream's source verbatim -- do not "fix" this `+1` by
            // analogy with it.
            summit: start + *offset as Coord + 1,
            pileup: sc.treat,
            pscore: p,
            qscore: q,
            fold_change: fold_enrichment(sc.treat, sc.ctrl, params.pseudocount),
            broad: false,
        });
    }
    Ok(peaks)
}
pub fn call_peaks_chromosome(
    chrom: ChromId,
    chunks: &[Chunk],
    params: &CallParams,
    scores: &[(&[f32], f32)],
    pq: &PqTable,
    cache: &mut macs_score::PScoreCache,
) -> (Vec<Peak>, usize) {
    call_peaks_chromosome_with_p(chrom, chunks, params, scores, pq, None, cache)
}

/// F120: as [`call_peaks_chromosome`], with the per-position p-score array so the
/// summit's reported scores come from the same track the histogram was built
/// from.
pub fn call_peaks_chromosome_with_p(
    chrom: ChromId,
    chunks: &[Chunk],
    params: &CallParams,
    scores: &[(&[f32], f32)],
    pq: &PqTable,
    pscore_track: Option<&[f32]>,
    cache: &mut macs_score::PScoreCache,
) -> (Vec<Peak>, usize) {
    let regions = segment_region_ranges(chunks, params.max_gap);
    let mut peaks = Vec::new();
    let mut rejected = 0usize;
    for range in &regions {
        let region = &chunks[range.clone()];
        let r = if params.call_summits {
            close_peak_with_subpeaks(region, params, scores, pq, Some(cache)).map(|v| (v, true))
        } else {
            close_peak_wo_subpeaks_with_p(region, params, scores, pq, pscore_track, Some(cache))
                .map(|p| (vec![p], false))
        };
        match r {
            Ok((v, _)) => peaks.extend(v),
            Err(_) => rejected += 1,
        }
    }
    let _ = chrom;
    (peaks, rejected)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(start: Coord, end: Coord, treat: f32) -> Chunk {
        Chunk {
            start,
            end,
            treat,
            ctrl: 1.0,
            score_index: 0,
        }
    }

    #[test]
    fn a_gap_of_exactly_max_gap_is_merged() {
        // upstream tests `if tl <= max_gap`, so the boundary case merges
        let chunks = vec![c(0, 10, 5.0), c(60, 70, 5.0)];
        let regions = segment_regions(&chunks, 50);
        assert_eq!(regions.len(), 1, "gap of exactly 50 must merge");
        let chunks = vec![c(0, 10, 5.0), c(61, 70, 5.0)];
        let regions = segment_regions(&chunks, 50);
        assert_eq!(regions.len(), 2, "gap of 51 must split");
    }

    #[test]
    fn regions_are_split_on_large_gaps() {
        let chunks = vec![c(0, 10, 5.0), c(20, 30, 5.0), c(500, 510, 5.0)];
        let regions = segment_regions(&chunks, 50);
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].len(), 2);
        assert_eq!(regions[1].len(), 1);
    }

    #[test]
    fn no_chunks_yields_no_regions() {
        assert!(segment_regions(&[], 50).is_empty());
    }

    #[test]
    fn a_region_shorter_than_min_length_is_rejected() {
        let region = vec![c(0, 100, 5.0)];
        let params = CallParams {
            min_length: 200,
            ..Default::default()
        };
        let pq = PqTable::empty();
        let mut cache = macs_score::PScoreCache::new();
        let r = close_peak_wo_subpeaks(&region, &params, &[], &pq, Some(&mut cache));
        assert_eq!(r, Err(Reject::TooShort));
    }

    #[test]
    fn a_region_of_exactly_min_length_is_accepted() {
        // upstream compares `>= min_length` with no +1 anywhere
        let region = vec![c(0, 200, 5.0)];
        let params = CallParams {
            min_length: 200,
            ..Default::default()
        };
        let pq = PqTable::empty();
        let mut cache = macs_score::PScoreCache::new();
        let r = close_peak_wo_subpeaks(&region, &params, &[], &pq, Some(&mut cache));
        assert!(r.is_ok(), "a peak of exactly min_length must be kept");
    }

    #[test]
    fn the_summit_is_the_chunk_midpoint() {
        let region = vec![c(0, 100, 5.0), c(100, 300, 9.0), c(300, 400, 5.0)];
        let params = CallParams {
            min_length: 200,
            ..Default::default()
        };
        let pq = PqTable::empty();
        let mut cache = macs_score::PScoreCache::new();
        let p = close_peak_wo_subpeaks(&region, &params, &[], &pq, Some(&mut cache)).unwrap();
        assert_eq!(p.summit, 201, "F62: (100+300+1).div_ceil(2)");
        assert_eq!(p.pileup, 9.0);
    }

    #[test]
    fn tied_maxima_pick_the_lower_median_midpoint() {
        // two chunks tie at the maximum: the lower median is the earlier one
        let region = vec![c(0, 100, 9.0), c(100, 300, 9.0), c(300, 500, 1.0)];
        let params = CallParams {
            min_length: 200,
            ..Default::default()
        };
        let pq = PqTable::empty();
        let mut cache = macs_score::PScoreCache::new();
        let p = close_peak_wo_subpeaks(&region, &params, &[], &pq, Some(&mut cache)).unwrap();
        assert_eq!(p.summit, 51, "F62: (0+100+1).div_ceil(2)");
    }

    #[test]
    fn three_way_tie_takes_the_middle_midpoint() {
        // midindex = (3 + 1) / 2 - 1 = 1, so the *second* of three
        let region = vec![
            c(0, 100, 9.0),
            c(100, 300, 9.0),
            c(300, 500, 9.0),
            c(500, 700, 1.0),
        ];
        let params = CallParams {
            min_length: 200,
            ..Default::default()
        };
        let pq = PqTable::empty();
        let mut cache = macs_score::PScoreCache::new();
        let p = close_peak_wo_subpeaks(&region, &params, &[], &pq, Some(&mut cache)).unwrap();
        assert_eq!(p.summit, 201, "F62: (100+300+1).div_ceil(2)");
    }

    #[test]
    fn gap_maxima_participate_in_peakyness_before_chunk_mapping() {
        // Reproduces the reported CTCF single-end peak. The smoothing maxima
        // include a one-base gap maximum. Upstream tests peakyness before mapping
        // maxima back to chunks; this rejects all subpeaks and uses the ordinary
        // parent-peak fallback, whose tied lower-median summit is 30642880.
        let region = vec![
            c(30_642_747, 30_642_844, 3.0),
            c(30_642_879, 30_642_880, 3.0),
            c(30_642_880, 30_642_975, 3.0),
            c(30_643_035, 30_643_073, 3.0),
        ];
        let params = CallParams {
            min_length: 228,
            call_summits: true,
            ..Default::default()
        };
        let pq = PqTable::empty();
        let mut cache = macs_score::PScoreCache::new();
        let peaks = close_peak_with_subpeaks(&region, &params, &[], &pq, Some(&mut cache))
            .expect("the parent peak should survive the ordinary fallback");
        assert_eq!(peaks.len(), 1);
        assert_eq!(peaks[0].summit, 30_642_880);
    }

    #[test]
    fn a_zero_pileup_region_still_yields_a_summit() {
        // upstream's `if not summit_value` treats 0 as unset
        let region = vec![c(0, 200, 0.0)];
        let params = CallParams {
            min_length: 200,
            ..Default::default()
        };
        let pq = PqTable::empty();
        let mut cache = macs_score::PScoreCache::new();
        let r = close_peak_wo_subpeaks(&region, &params, &[], &pq, Some(&mut cache));
        assert!(r.is_ok(), "zero pileup still produces a peak");
    }

    #[test]
    fn a_summit_below_the_cutoff_rejects_the_whole_peak() {
        let region = vec![c(0, 100, 5.0), c(100, 400, 9.0)];
        let params = CallParams {
            min_length: 200,
            ..Default::default()
        };
        // score at the summit index is below the cutoff
        let score_array = vec![1.0f32, 2.0];
        let pq = PqTable::empty();
        let mut cache = macs_score::PScoreCache::new();
        let r = close_peak_wo_subpeaks(
            &region,
            &params,
            &[(&score_array, 5.0)],
            &pq,
            Some(&mut cache),
        );
        assert_eq!(r, Err(Reject::BelowCutoff));
    }

    #[test]
    fn a_passing_peak_reports_its_bounds_from_first_and_last_chunk() {
        let region = vec![c(10, 100, 5.0), c(100, 300, 9.0), c(300, 350, 2.0)];
        let params = CallParams {
            min_length: 200,
            ..Default::default()
        };
        let pq = PqTable::empty();
        let mut cache = macs_score::PScoreCache::new();
        let p = close_peak_wo_subpeaks(&region, &params, &[], &pq, Some(&mut cache)).unwrap();
        assert_eq!(p.start, 11, "F85: XLS start is the first chunk start + 1");
        assert_eq!(p.end, 350, "end comes from the last chunk");
        assert_eq!(p.length(), 339);
    }

    /// F259: the fold enrichment is **stored** narrowed to f32.
    ///
    /// `CallPeakUnit.py:1385` divides in C double (C's usual arithmetic conversions
    /// promote the `cython.float` operands against the `cython.double` pseudocount),
    /// but the result lands in `PeakContent.fc`, declared `cython.float`. These are
    /// the real operands from `sweep/gonechrom_mpe_d4000_w180_noc_131`, where the
    /// control lambda, the treatment pileup and the summit are all already
    /// byte-identical to upstream -- the stored width is the only thing left.
    #[test]
    fn fold_enrichment_is_narrowed_to_f32_like_peakcontent_fc() {
        let treat = 180.0f32;
        let ctrl = f32::from_bits(0x4092_e147); // 4.589999675750732

        let got = super::fold_enrichment(treat, ctrl, 1.0);

        // The f64 divide alone is 32.379250536484484 (f64 bits 0x4040308b4815987d),
        // which upstream's `%.6g` prints as 32.3793 -- the wrong digit. Narrowed it is
        // 32.379249572753906 (f32 bits 0x4201845a), printed as 32.3792 -- upstream's.
        let wide = (treat as f64 + 1.0) / (ctrl as f64 + 1.0);
        assert_eq!(
            wide.to_bits(),
            0x4040_308b_4815_987d,
            "the f64 divide is a different number than the one upstream prints"
        );
        assert_ne!(
            got.to_bits(),
            wide.to_bits(),
            "F259: keeping the f64 divide reproduces the wrong last digit"
        );
        assert_eq!(
            got as f32,
            f32::from_bits(0x4201_845a),
            "must be the f32-narrowed f64 divide"
        );
    }

    /// F259: the *summit* path and the *broad* averaging path must both narrow, and
    /// both must agree on the same operands.
    #[test]
    fn fold_enrichment_is_f32_at_the_extremes() {
        for (t, c) in [
            (0.0f32, 0.0f32),
            (1.0, 1.0),
            (180.0, 4.59),
            (8652.0, 4600.0),
        ] {
            let got = super::fold_enrichment(t, c, 1.0);
            let wide = (t as f64 + 1.0) / (c as f64 + 1.0);
            assert_eq!(
                got,
                f64::from(wide as f32),
                "t={t} c={c}: stored value must be the f32 narrowing"
            );
            assert!(
                got as f32 as f64 == got,
                "t={t} c={c}: must round-trip as f32"
            );
        }
    }
}
