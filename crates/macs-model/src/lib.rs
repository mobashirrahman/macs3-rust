//! Peak model building -- a port of `MACS3/Signal/PeakModel.py`.
//!
//! The model estimates the fragment length `d` from ChIP-seq tags by finding
//! paired `+`/`-` strand peaks and cross-correlating the two strand
//! distributions. `predictd` is the only consumer, but `callpeak` uses the same
//! `d` when the model path is enabled, so this is the G8 gate.
//!
//! The algorithm, in order:
//!
//! 1. `naive_quick_pileup` -- extend every tag by `peaksize/2` on each side and
//!    pile up, keeping **every** breakpoint (not a merged RLE).
//! 2. `naive_call_peaks` -- call summits on that pileup with the
//!    `lmfold`/`umfold` tag-count bounds and upstream's default `max_gap=50`,
//!    `min_length=200`.
//! 3. `find_paired_peaks` / `find_pair_center` -- pair `+` and `-` summits whose
//!    positions overlap within `peaksize`.
//! 4. `model_add_line` / `count` -- accumulate each strand's tags into a
//!    `window_size = 1 + 2*peaksize + tag_expansion_size` profile.
//! 5. Cross-correlate the two normalised profiles, smooth, take local maxima
//!    above `d_min` as the candidate `d`s, and the highest-correlation one as
//!    `d`.

use std::collections::HashMap;

use macs_core::{Coord, MacsError, Result, Strand};
use macs_track::SingleEndTrack;

/// Upstream's fixed tag-expansion constant (`PeakModel.py:84`).
pub const TAG_EXPANSION_SIZE: i64 = 10;

/// Not enough paired peaks to fit a model.
#[derive(Debug, thiserror::Error)]
#[error("MACS needs at least 100 paired peaks at + and - strand to build the model")]
pub struct NotEnoughPairs;

/// Parameters for [`PeakModel`].
#[derive(Debug, Clone)]
pub struct ModelOptions {
    /// `--gsize`: effective genome size.
    pub gsize: f64,
    /// `--mfold` as `(lmfold, umfold)`.
    pub mfold: (f64, f64),
    /// `--bw`: model band width; `peaksize = 2 * bw`.
    pub bw: i64,
    /// `--d-min`: smallest acceptable lag.
    pub d_min: i64,
}

impl Default for ModelOptions {
    fn default() -> Self {
        Self {
            gsize: 2.65e9,
            mfold: (5.0, 50.0),
            bw: 300,
            d_min: 20,
        }
    }
}

/// The fitted peak model.
/// `(plus_line, minus_line, ycorr, xcorr, d, alternative_d)`
type ModelParts = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, f64, Vec<i64>);

#[derive(Debug, Clone)]
pub struct PeakModel {
    pub d: f64,
    pub alternative_d: Vec<i64>,
    pub plus_line: Vec<f64>,
    pub minus_line: Vec<f64>,
    pub ycorr: Vec<f64>,
    pub xcorr: Vec<f64>,
    pub min_tags: f64,
    pub max_tags: f64,
    pub peaksize: i64,
    pub gsize: f64,
    pub scan_window: i64,
    /// Per-strand and paired peak counts from `find_paired_peaks`.
    ///
    /// Upstream prints these under `--verbose 3`
    /// (`PeakModel.py:121,177-192`) and nothing else exposes them, so they are the
    /// only way to tell *where* a model build diverges: a `min_tags`/`max_tags`
    /// disagreement, a strand that yielded fewer summits, or a pairing step that
    /// dropped candidates. Kept on the model so a caller can report them without
    /// re-deriving anything.
    pub paired: PairedCounts,
}

/// `(chromosome, plus summits, minus summits, paired centres)`.
type ChromCounts = Vec<(Vec<u8>, usize, usize, usize)>;

/// The pairing stage's output: paired centres per chromosome, plus the counts.
pub type PairedPeakResult = (HashMap<Vec<u8>, Vec<Coord>>, ChromCounts);

/// Counts observed while pairing `+`/`-` strand summits.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PairedCounts {
    /// Total paired peak centres found, summed over chromosomes.
    pub total_pairs: usize,
    /// `(chromosome, plus summits, minus summits, paired centres)`, in the order
    /// `find_paired_peaks` visits them.
    pub per_chrom: ChromCounts,
}

/// `naive_quick_pileup`: extend each sorted tag by `extension` both ways
/// (clamped at 0) and pile up, emitting **every** breakpoint.
pub fn naive_quick_pileup(sorted_poss: &[Coord], extension: i64) -> (Vec<Coord>, Vec<f32>) {
    if sorted_poss.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let e = extension.max(0) as u64;
    let mut starts: Vec<Coord> = sorted_poss.iter().map(|&p| p.saturating_sub(e)).collect();
    let mut ends: Vec<Coord> = sorted_poss.iter().map(|&p| p + e).collect();
    starts.sort_unstable();
    ends.sort_unstable();
    sweep_raw(&starts, &ends)
}

/// The endpoint sweep behind `naive_quick_pileup`, keeping every breakpoint.
fn sweep_raw(starts: &[Coord], ends: &[Coord]) -> (Vec<Coord>, Vec<f32>) {
    let mut p_out = Vec::with_capacity(starts.len() + ends.len());
    let mut v_out = Vec::with_capacity(starts.len() + ends.len());
    let mut pileup: i64 = 0;
    let mut pre_p: Coord = starts[0].min(ends[0]);
    if pre_p != 0 {
        p_out.push(pre_p);
        v_out.push(0.0);
    }
    let (mut i_s, mut i_e) = (0usize, 0usize);
    while i_s < starts.len() && i_e < ends.len() {
        if starts[i_s] < ends[i_e] {
            let p = starts[i_s];
            if p != pre_p {
                p_out.push(p);
                v_out.push(pileup as f32);
                pre_p = p;
            }
            pileup += 1;
            i_s += 1;
        } else if starts[i_s] > ends[i_e] {
            let p = ends[i_e];
            if p != pre_p {
                p_out.push(p);
                v_out.push(pileup as f32);
                pre_p = p;
            }
            pileup -= 1;
            i_e += 1;
        } else {
            // Coincident start and end: upstream consumes **both** events without
            // emitting a breakpoint and without touching `pileup`
            // (`PileupV2.py:419-423`).
            //
            // This branch used to emit a breakpoint and do `pileup -= 1`. Both are
            // wrong, and together they cost one unit of depth at *every* coincident
            // position. With `extension = peaksize/2 = 300` on a real ChIP-seq file
            // that is most positions, so the strand pileup was systematically too
            // shallow and `naive_call_peaks` shattered long above-cutoff runs into
            // fragments shorter than `min_length`. Measured on upstream's own
            // `CTCF_SE_ChIP_chr22_50k.bed.gz`: 235 plus-strand summits instead of
            // 1722, so pairing found 66 centres instead of 464, below the 100
            // upstream's `NotEnoughPairsException` needs -- and the fragment model
            // could not be fitted at all. Pinned by
            // `tests/paired_peaks_vs_upstream.rs`.
            i_s += 1;
            i_e += 1;
        }
    }
    // drain remaining ends, then remaining starts (mirrors upstream's two tails)
    while i_e < ends.len() {
        let p = ends[i_e];
        if p != pre_p {
            p_out.push(p);
            v_out.push(pileup as f32);
            pre_p = p;
        }
        pileup -= 1;
        i_e += 1;
    }
    while i_s < starts.len() {
        let p = starts[i_s];
        if p != pre_p {
            p_out.push(p);
            v_out.push(pileup as f32);
            pre_p = p;
        }
        pileup += 1;
        i_s += 1;
    }
    (p_out, v_out)
}

/// Close one peak: the lower-median midpoint of the maximal-value runs, kept
/// only if its height is below `max_v`.
fn close_peak(content: &[(Coord, Coord, f32)], max_v: f32, out: &mut Vec<(Coord, f32)>) {
    let mut tsummit: Vec<Coord> = Vec::new();
    let mut summit_value = 0.0f32;
    for &(ts, te, tv) in content {
        if summit_value == 0.0 || summit_value < tv {
            tsummit.clear();
            tsummit.push((te + ts) / 2);
            summit_value = tv;
        } else if summit_value == tv {
            tsummit.push((te + ts) / 2);
        }
    }
    if tsummit.is_empty() {
        return;
    }
    // upstream is `(len + 1) / 2 - 1`, the lower median -- deliberately not
    // `div_ceil(2)`, which picks the later of two tied runs (cf. F47).
    #[allow(clippy::manual_div_ceil)]
    let mid = (tsummit.len() + 1) / 2 - 1;
    if summit_value < max_v {
        out.push((tsummit[mid], summit_value));
    }
}

/// `naive_call_peaks` with upstream's defaults (`max_gap=50`, `min_length=200`).
pub fn naive_call_peaks(pv: &(Vec<Coord>, Vec<f32>), min_v: f32, max_v: f32) -> Vec<(Coord, f32)> {
    const MAX_GAP: Coord = 50;
    const MIN_LENGTH: Coord = 200;
    let (ps, vs) = pv;
    let mut ret: Vec<(Coord, f32)> = Vec::new();
    let mut content: Vec<(Coord, Coord, f32)> = Vec::new();
    let mut pre_p: Coord = 0;
    let mut i = 0usize;
    while i < ps.len() {
        let (p, v) = (ps[i], vs[i]);
        i += 1;
        if v > min_v {
            content.push((pre_p, p, v));
            pre_p = p;
            break;
        }
        pre_p = p;
    }
    while i < ps.len() {
        let (p, v) = (ps[i], vs[i]);
        i += 1;
        if v <= min_v {
            pre_p = p;
            continue;
        }
        let last_end = content.last().map(|c| c.1).unwrap_or(0);
        if pre_p.saturating_sub(last_end) <= MAX_GAP {
            content.push((pre_p, p, v));
        } else {
            let len = content.last().unwrap().1.saturating_sub(content[0].0);
            if len >= MIN_LENGTH {
                close_peak(&content, max_v, &mut ret);
            }
            content = vec![(pre_p, p, v)];
        }
        pre_p = p;
    }
    if !content.is_empty() {
        let len = content.last().unwrap().1.saturating_sub(content[0].0);
        if len >= MIN_LENGTH {
            close_peak(&content, max_v, &mut ret);
        }
    }
    ret
}

impl PeakModel {
    /// Build the model from single-end tags.
    pub fn build(track: &SingleEndTrack, options: ModelOptions) -> Result<Self> {
        let peaksize = 2 * options.bw;
        let total = track.total() as f64;
        let (lmfold, umfold) = options.mfold;
        let min_tags = (total * lmfold * peaksize as f64 / options.gsize / 2.0).round();
        let max_tags = (total * umfold * peaksize as f64 / options.gsize / 2.0).round();

        let (paired, counts) = find_paired_peaks(track, peaksize, min_tags, max_tags);
        let num: usize = paired.values().map(Vec::len).sum();
        let counts = PairedCounts {
            total_pairs: num,
            per_chrom: counts,
        };
        if num < 100 {
            // F200: `predictd_cmd.py:82-83` catches upstream's `NotEnoughPairsException`
            // and merely warns, so the command still exits 0. Surfacing a hard error
            // here made `macs3-rs predictd` exit 1 on an input upstream accepts.
            return Err(MacsError::NotEnoughPairs {
                found: num as u64,
                needed: 100,
            });
        }

        let (plus_line, minus_line, ycorr, xcorr, d, alternative_d) =
            paired_peak_model(&paired, track, peaksize, options.d_min);

        let scan_window = (d as i64).max(TAG_EXPANSION_SIZE) * 2;
        Ok(Self {
            d,
            alternative_d,
            plus_line,
            minus_line,
            ycorr,
            xcorr,
            min_tags,
            max_tags,
            peaksize,
            gsize: options.gsize,
            scan_window,
            paired: counts,
        })
    }
}

/// `__find_paired_peaks` + `__naive_find_peaks`.
///
/// Also returns the per-strand summit counts, which upstream reports under
/// `--verbose 3` and which localise a divergence to either the strand pileup or
/// the pairing step.
fn find_paired_peaks(
    track: &SingleEndTrack,
    peaksize: i64,
    min_tags: f64,
    max_tags: f64,
) -> PairedPeakResult {
    let mut out: HashMap<Vec<u8>, Vec<Coord>> = HashMap::new();
    let mut counts: ChromCounts = Vec::new();
    for chrom in track.positions().chroms_sorted() {
        let plus = naive_call_peaks(
            &naive_quick_pileup(track.positions().strand(chrom, Strand::Plus), peaksize / 2),
            min_tags as f32,
            max_tags as f32,
        );
        let minus = naive_call_peaks(
            &naive_quick_pileup(track.positions().strand(chrom, Strand::Minus), peaksize / 2),
            min_tags as f32,
            max_tags as f32,
        );
        if plus.is_empty() || minus.is_empty() {
            // upstream discards the chromosome entirely (`PeakModel.py:187-188`)
            counts.push((
                track.genome().name(chrom).to_vec(),
                plus.len(),
                minus.len(),
                0,
            ));
            continue;
        }
        let centers = find_pair_center(&plus, &minus, peaksize);
        counts.push((
            track.genome().name(chrom).to_vec(),
            plus.len(),
            minus.len(),
            centers.len(),
        ));
        out.insert(track.genome().name(chrom).to_vec(), centers);
    }
    (out, counts)
}

/// `__find_pair_center`: pair `+`/`-` summits overlapping within `peaksize`.
fn find_pair_center(
    plus_peaks: &[(Coord, f32)],
    minus_peaks: &[(Coord, f32)],
    peaksize: i64,
) -> Vec<Coord> {
    let (mut ip, mut im, mut im_prev) = (0usize, 0usize, 0usize);
    let mut centers = Vec::new();
    let mut overlap = false;
    while ip < plus_peaks.len() && im < minus_peaks.len() {
        let pp = plus_peaks[ip].0;
        let mp = minus_peaks[im].0;
        let (pp_i, mp_i) = (pp as i64, mp as i64);
        if pp_i - peaksize > mp_i {
            im += 1;
        } else if pp_i + peaksize < mp_i {
            ip += 1;
            im = im_prev;
            overlap = false;
        } else {
            if !overlap {
                overlap = true;
                im_prev = im;
            }
            // only pair peaks whose tag counts are comparable, and only when the
            // plus summit lies left of the minus one (upstream `PeakModel.py:440-445`)
            let (pn, mn) = (f64::from(plus_peaks[ip].1), f64::from(minus_peaks[im].1));
            if pn / mn < 2.0 && pn / mn > 0.5 && pp < mp {
                centers.push((pp + mp) / 2);
            }
            im += 1;
        }
    }
    centers
}

/// `__model_add_line`: place tags into the `start`/`end` profile arrays.
fn model_add_line(
    paired_centers: &[Coord],
    tags: &[Coord],
    start: &mut [i64],
    end: &mut [i64],
    peaksize: i64,
) {
    let max_index = start.len() - 1;
    let pad = peaksize + TAG_EXPANSION_SIZE / 2;
    let half = TAG_EXPANSION_SIZE / 2;
    let mut i1 = 0usize;
    let mut i2 = 0usize;
    let mut i2_prev = 0usize;
    let mut overlap = false;
    while i1 < paired_centers.len() && i2 < tags.len() {
        let p1 = paired_centers[i1];
        let p2 = tags[i2];
        let (p1_i, p2_i) = (p1 as i64, p2 as i64);
        if p1_i - pad > p2_i {
            i2 += 1;
        } else if p1_i + pad < p2_i {
            i1 += 1;
            i2 = i2_prev;
            overlap = false;
        } else {
            if !overlap {
                overlap = true;
                i2_prev = i2;
            }
            let s = ((p2_i - half - p1_i + pad).max(0) as usize).min(max_index);
            start[s] += 1;
            let e = (p2_i + half - p1_i + pad).clamp(0, max_index as i64) as usize;
            end[e] -= 1;
            i2 += 1;
        }
    }
}

/// `__count`: prefix-sum the `+1`/`-1` events into a profile line.
fn count(start: &[i64], end: &[i64], line: &mut [f64]) {
    let mut acc = 0i64;
    for i in 0..line.len() {
        acc += start[i] + end[i];
        line[i] = acc as f64;
    }
}

/// `smooth(x, window_len=11, window='flat')` -- reflected-edge moving average.
fn smooth_flat(x: &[f64], window_len: usize) -> Vec<f64> {
    let n = x.len();
    if n < window_len || window_len < 3 {
        return x.to_vec();
    }
    // reflect: x[w-1:0:-1] ++ x ++ x[-1:-w:-1]
    let mut s: Vec<f64> = Vec::with_capacity(n + 2 * (window_len - 1));
    for i in (1..window_len).rev() {
        s.push(x[i]);
    }
    s.extend_from_slice(x);
    for i in 1..window_len {
        s.push(x[n - i]);
    }
    // F201: upstream is `np.convolve(w/w.sum(), s, mode='valid')`
    // (`PeakModel.py:505-506`). That is **not** "sum the window then divide" -- the
    // window is divided first, so every product is individually scaled and rounded:
    //
    //     w = np.ones(window_len); kernel = w / w.sum()      # each 1/window_len
    //     y = np.convolve(kernel, s, 'valid')
    //
    // Scaling each term before summing is not the same operation as summing and
    // scaling once: on this window the two disagree in 31 of 50 elements, which is
    // exactly enough to put different last digits into the `ycorr <- c(...)` line of
    // `*_model.r`, and `*_model.r` is an acceptance-gate byte-identical file.
    //
    // The indices are visited in **forward** `j` order; the flat kernel is symmetric so
    // the index *set* is orientation-independent, but the accumulation order is
    // observable. Empirically `sum_j (1/n) * s[k+j]` matches `np.convolve` exactly and
    // `sum_j (1/n) * s[k+n-1-j]` does not.
    let kernel = 1.0f64 / window_len as f64;
    let ylen = s.len() - window_len + 1;
    let mut y = vec![0.0f64; ylen];
    for (k, yk) in y.iter_mut().enumerate() {
        let mut acc = 0.0;
        for &v in &s[k..k + window_len] {
            acc += kernel * v;
        }
        *yk = acc;
    }
    let half = window_len / 2;
    y[half..ylen - half].to_vec()
}

/// `__paired_peak_model`: build the profiles and cross-correlate them.
fn paired_peak_model(
    paired: &HashMap<Vec<u8>, Vec<Coord>>,
    track: &SingleEndTrack,
    peaksize: i64,
    d_min: i64,
) -> ModelParts {
    let window_size = (1 + 2 * peaksize + TAG_EXPANSION_SIZE) as usize;
    let mut plus_start = vec![0i64; window_size];
    let mut plus_end = vec![0i64; window_size];
    let mut minus_start = vec![0i64; window_size];
    let mut minus_end = vec![0i64; window_size];
    let mut chroms: Vec<&Vec<u8>> = paired.keys().collect();
    chroms.sort();
    for chrom in chroms {
        let Some(cid) = track.genome().get(chrom) else {
            continue;
        };
        let centers = &paired[chrom];
        model_add_line(
            centers,
            track.positions().strand(cid, Strand::Plus),
            &mut plus_start,
            &mut plus_end,
            peaksize,
        );
        model_add_line(
            centers,
            track.positions().strand(cid, Strand::Minus),
            &mut minus_start,
            &mut minus_end,
            peaksize,
        );
    }
    let mut plus_line = vec![0.0f64; window_size];
    let mut minus_line = vec![0.0f64; window_size];
    count(&plus_start, &plus_end, &mut plus_line);
    count(&minus_start, &minus_end, &mut minus_line);

    // normalise, then full cross-correlation truncated to +/- peaksize
    let norm = |line: &[f64]| -> Vec<f64> {
        let n = line.len() as f64;
        let mean = np_mean(line);
        let sd = np_std(line, mean);
        line.iter().map(|&v| (v - mean) / (sd * n)).collect()
    };
    let plus_data = norm(&plus_line);
    let minus_data = norm(&minus_line);
    let full = correlate(&minus_data, &plus_data);
    let lo = window_size - peaksize as usize;
    let hi = window_size + peaksize as usize;
    let mut ycorr: Vec<f64> = full[lo..hi].to_vec();
    // upstream: `np.linspace(len(ycorr)//2*-1, len(ycorr)//2, num=len(ycorr))`.
    // Note it is built on the UNSMOOTHED length and spans the same count of
    // points, so the lag axis is *not* the integers -half..=half (which would be
    // one element longer and shift every candidate d).
    let xcorr: Vec<f64> = linspace(
        -((ycorr.len() / 2) as f64),
        (ycorr.len() / 2) as f64,
        ycorr.len(),
    );
    ycorr = smooth_flat(&ycorr, 11);

    // local maxima above d_min, strongest first
    let mut idx: Vec<usize> = Vec::new();
    for i in 1..ycorr.len().saturating_sub(1) {
        if ycorr[i] > ycorr[i - 1] && ycorr[i] > ycorr[i + 1] && xcorr[i] > d_min as f64 {
            idx.push(i);
        }
    }
    idx.sort_by(|&a, &b| ycorr[b].partial_cmp(&ycorr[a]).expect("finite correlation"));
    let d = if let Some(&i) = idx.first() {
        xcorr[i]
    } else {
        0.0
    };
    let mut alternative_d: Vec<i64> = idx.iter().map(|&i| xcorr[i] as i64).collect();
    alternative_d.sort_unstable();
    (plus_line, minus_line, ycorr, xcorr, d, alternative_d)
}

/// `np.linspace(start, stop, num)`.
fn linspace(start: f64, stop: f64, num: usize) -> Vec<f64> {
    if num <= 1 {
        return vec![start];
    }
    let step = (stop - start) / (num as f64 - 1.0);
    (0..num).map(|i| start + step * i as f64).collect()
}

/// NumPy's **pairwise** summation (`PW_BLOCKSIZE = 128`, 8-way unrolled).
///
/// This is not an optimisation detail: float addition is not associative, and
/// NumPy's `mean`/`std` use this routine while a naive `iter().sum()` is
/// sequential. The two disagree in the last bits, and that difference survives
/// normalisation, correlation and smoothing -- enough to change the digits
/// written into `*_model.r`. Reproducing the exact traversal is what makes the
/// acceptance criterion (byte-identical `*_model.r`) achievable.
fn numpy_pairwise_sum(a: &[f64]) -> f64 {
    const PW_BLOCKSIZE: usize = 128;
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &v in a {
            res += v;
        }
        return res;
    }
    if n <= PW_BLOCKSIZE {
        let mut r = [0.0f64; 8];
        r.copy_from_slice(&a[..8]);
        let mut i = 8usize;
        while i < n - (n % 8) {
            for k in 0..8 {
                r[k] += a[i + k];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        return res;
    }
    let mut n2 = n / 2;
    n2 -= n2 % 8;
    numpy_pairwise_sum(&a[..n2]) + numpy_pairwise_sum(&a[n2..])
}

/// NumPy `ndarray.mean()`: pairwise sum divided by `n`.
fn np_mean(a: &[f64]) -> f64 {
    numpy_pairwise_sum(a) / a.len() as f64
}

/// NumPy `ndarray.std()` (ddof = 0): `sqrt(mean(|x - mean|^2))`, with the
/// squared deviations themselves pairwise-summed.
fn np_std(a: &[f64], mean: f64) -> f64 {
    let devs: Vec<f64> = a.iter().map(|&v| (v - mean) * (v - mean)).collect();
    (numpy_pairwise_sum(&devs) / a.len() as f64).sqrt()
}

/// `np.correlate(a, b, 'full')`.
///
/// `out[k] = sum_j a[j + lag] * b[j]` with `lag = k - (len(b) - 1)` -- i.e. the
/// pair index difference is `k - (len(b)-1)`, **not** `i + j == k`. Getting that
/// wrong mirrors the lag axis and moves the correlation peak.
fn correlate(a: &[f64], b: &[f64]) -> Vec<f64> {
    let lb = b.len() as i64;
    let n = a.len() + b.len() - 1;
    let mut out = vec![0.0f64; n];
    for (k, slot) in out.iter_mut().enumerate() {
        let lag = k as i64 - (lb - 1);
        let mut acc = 0.0;
        for (j, &bv) in b.iter().enumerate() {
            let i = j as i64 + lag;
            if i >= 0 && (i as usize) < a.len() {
                acc += a[i as usize] * bv;
            }
        }
        *slot = acc;
    }
    out
}

/// `model2r_script` (`OutputWriter.py:227-278`): the `*_model.r` acceptance file.
pub fn model2r_script(model: &PeakModel, filename: &std::path::Path, name: &str) -> Result<()> {
    let w = model.plus_line.len();
    let sum_p: f64 = model.plus_line.iter().sum();
    let sum_m: f64 = model.minus_line.iter().sum();
    let norm_p: Vec<String> = model
        .plus_line
        .iter()
        .map(|&v| fmt_num(v * 100.0 / sum_p))
        .collect();
    let norm_m: Vec<String> = model
        .minus_line
        .iter()
        .map(|&v| fmt_num(v * 100.0 / sum_m))
        .collect();
    let ycorr: Vec<String> = model.ycorr.iter().map(|&v| fmt_num(v)).collect();
    let xcorr: Vec<String> = model.xcorr.iter().map(|&v| fmt_num(v)).collect();
    let alt_str: Vec<String> = model.alternative_d.iter().map(|&v| v.to_string()).collect();
    let alt = alt_str.join(", ");
    let mut out = String::new();
    out.push_str("# R script for Peak Model\n");
    out.push_str("#  -- generated by MACS\n");
    out.push_str(&format!("p <- c({})\n", norm_p.join(",")));
    out.push_str(&format!("m <- c({})\n", norm_m.join(",")));
    out.push_str(&format!("ycorr <- c({})\n", ycorr.join(",")));
    out.push_str(&format!("xcorr <- c({})\n", xcorr.join(",")));
    out.push_str(&format!("altd  <- c({})\n", alt));
    out.push_str("x <- seq.int((length(p)-1)/2*-1,(length(p)-1)/2)\n");
    out.push_str(&format!("pdf('{name}_model.pdf',height=6,width=6)\n"));
    out.push_str("plot(x,p,type='l',col=c('red'),main='Peak Model',xlab='Distance to the middle',ylab='Percentage')\n");
    out.push_str("lines(x,m,col=c('blue'))\n");
    out.push_str(
        "legend('topleft',c('forward tags','reverse tags'),lty=c(1,1,1),col=c('red','blue'))\n",
    );
    out.push_str("plot(xcorr,ycorr,type='l',col=c('black'),main='Cross-Correlation',xlab='Lag between + and - tags',ylab='Correlation')\n");
    out.push_str("abline(v=altd,lty=2,col=c('red'))\n");
    out.push_str("legend('topleft','alternative lag(s)',lty=2,col='red')\n");
    out.push_str(&format!("legend('right','alt lag(s) : {alt}',bty='n')\n"));
    out.push_str("dev.off()\n");
    let _ = w;
    std::fs::write(filename, out).map_err(MacsError::Io)?;
    Ok(())
}

/// Python `str(float)` (== `repr`): the shortest round-tripping digit string,
/// in scientific notation when the decimal exponent is `< -4` or `>= 16`, else
/// plain fixed-point.
///
/// The `*_model.r` acceptance criterion is byte-identity, and upstream writes
/// these values with Python's `%s` on a float -- e.g. `3.098604902685399e-05`,
/// not `0.000030986049026853796`. Both name the same number; only the rendering
/// differs, so this must be reproduced rather than reformatted.
fn fmt_num(v: f64) -> String {
    if v == 0.0 {
        return "0.0".to_string();
    }
    if v.is_nan() {
        return "nan".to_string();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf".into() } else { "-inf".into() };
    }
    // shortest round-trip mantissa + exponent via Rust's LowerExp
    let sci = format!("{:e}", v); // e.g. "3.098604902685399e-5"
    let (mant, exp_str) = sci.split_once('e').expect("LowerExp has an exponent");
    let exp: i32 = exp_str.parse().expect("LowerExp exponent parses");
    if !(-4..16).contains(&exp) {
        // scientific: mantissa with a trailing ".0" if integral, 2-digit exponent
        let mant = if mant.contains('.') {
            mant.to_string()
        } else {
            format!("{mant}.0")
        };
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{mant}e{sign}{:02}", exp.abs())
    } else {
        // fixed point; Python always shows a decimal point on a float
        let t = format!("{v}");
        if t.contains('.') {
            t
        } else {
            format!("{t}.0")
        }
    }
}
/// The pairing stage on its own, for differential diagnosis.
///
/// [`PeakModel::build`] returns `NotEnoughPairs` before producing a model, which
/// discards the very counts needed to explain *why* pairing fell short. This
/// exposes the stage so a caller (or a differential test) can compare per-strand
/// summits and paired centres against upstream's `--verbose 3` numbers.
///
/// Prefer [`PeakModel::build`]; this is for localisation.
pub mod probe {
    use super::*;

    /// `(paired centres by chromosome, per-chromosome (name, plus, minus, paired))`.
    pub fn find_paired_peaks_pub(
        track: &SingleEndTrack,
        peaksize: i64,
        min_tags: f64,
        max_tags: f64,
    ) -> PairedPeakResult {
        super::find_paired_peaks(track, peaksize, min_tags, max_tags)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use macs_core::Strand;
    use macs_track::SingleEndTrackBuilder;

    /// Build a track with `n` plus clusters each followed by a minus cluster
    /// `d` bp downstream, so the model should recover `d`.
    fn synth(n: usize, d: i64) -> SingleEndTrack {
        let mut b = SingleEndTrackBuilder::new();
        for i in 0..n {
            let c = (1000 + i * 1000) as u64;
            for k in 0..12u64 {
                let p = c + (k * 7 % 80);
                b.push(b"chr1", p, Strand::Plus);
                b.push(b"chr1", p + d as u64, Strand::Minus);
            }
        }
        b.finalize();
        b.build()
    }

    #[test]
    fn recover_a_known_fragment_length() {
        let track = synth(300, 200);
        let opts = ModelOptions {
            gsize: 2e6,
            ..Default::default()
        };
        let m = PeakModel::build(&track, opts).expect("model fits");
        // the lag axis is fractional and the clusters overlap, so the recovered
        // lag sits a few bases off the true separation; require it to be close
        assert!(
            (m.d - 200.0).abs() <= 15.0,
            "recovered d={} should be near the true 200",
            m.d
        );
        assert!(
            m.alternative_d.iter().any(|&a| (a - 200).abs() <= 15),
            "an alternative lag should be near 200: {:?}",
            m.alternative_d
        );
    }

    #[test]
    fn correlate_matches_numpy_full() {
        // np.correlate([1,2,3],[10,20],'full') == [20,50,80,30]
        assert_eq!(
            correlate(&[1.0, 2.0, 3.0], &[10.0, 20.0]),
            vec![20.0, 50.0, 80.0, 30.0]
        );
    }

    /// F201: `smooth_flat` must scale each window element by `1/n` *before* summing.
    ///
    /// The vector below is the one used to demonstrate the divergence against
    /// NumPy: summing then dividing once disagrees with `np.convolve(w/w.sum(), ...)`
    /// in 31 of these 50 outputs, while per-term scaling matches all 50.
    #[test]
    fn smooth_flat_scales_each_term_before_summing() {
        let x: Vec<f64> = (0..50)
            .map(|i| ((i as f64 * 37.0) % 23.0) - 11.0 + (i as f64) * 0.37)
            .collect();
        let got = smooth_flat(&x, 11);
        assert_eq!(got.len(), x.len());

        let mut s: Vec<f64> = Vec::new();
        for i in (1..11).rev() {
            s.push(x[i]);
        }
        s.extend_from_slice(&x);
        for i in 1..11 {
            s.push(x[x.len() - i]);
        }
        let kernel = 1.0f64 / 11.0f64;

        // Exact reference: per-term scaling, forward traversal.
        let mut expect: Vec<f64> = Vec::with_capacity(50);
        for k in 0..50usize {
            let start = k + 5;
            let mut acc = 0.0f64;
            for &v in &s[start..start + 11] {
                acc += kernel * v;
            }
            expect.push(acc);
        }
        assert_eq!(got, expect);

        // And explicitly NOT "sum then divide", which is the bug this test guards.
        let mut wrong_count = 0;
        for (k, &g) in got.iter().enumerate() {
            let start = k + 5;
            let sum: f64 = s[start..start + 11].iter().sum();
            if sum / 11.0 != g {
                wrong_count += 1;
            }
        }
        assert!(
            wrong_count > 0,
            "the two orderings converged on this input; the test no longer discriminates"
        );
    }

    #[test]
    fn linspace_matches_numpy() {
        let v = linspace(-600.0, 600.0, 1200);
        assert_eq!(v.len(), 1200);
        assert_eq!(v[0], -600.0);
        assert!((v[1199] - 600.0).abs() < 1.5);
    }
}
