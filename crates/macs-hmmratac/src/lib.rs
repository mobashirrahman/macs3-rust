//! HMMRATAC support: fragment-length weighting, bin extraction and HMM decoding.
//!
//! A port of `MACS3/Signal/HMMR_Signal_Processing.py`, `HMMR_HMM.py` and the
//! `pileup_from_LR_hmmratac` helper in `PileupV2.py`.
//!
//! # Scope
//!
//! MACS3 trains its three-state HMM with `hmmlearn`'s `GaussianHMM`/`PoissonHMM`
//! (a Baum-Welch fit seeded from `numpy.random.MT19937`). Reproducing an EM
//! fit bit-for-bit would mean reproducing an SVD, a QR and a specific RNG
//! draw sequence, so **self-training is the one declared deviation** of this
//! port (see the acceptance criteria).
//!
//! Everything downstream of the model is exact:
//!
//! * the fragment-length weight mapping, in `f32` as upstream;
//! * the four digested pileup tracks, via the same endpoint sweep;
//! * bin extraction out of peak regions, including the `mark_bin` grouping;
//! * the HMM *inference* (forward-backward posteriors) against an exported
//!   model, which is what `--hmm` does.
//!
//! [`GaussianHMM::posterior`] and [`PoissonHMM::posterior`] implement the
//! forward-backward recursion with hmmlearn's own convention, which is
//! reproduced rather than "corrected" -- see the module on
//! [`HMM_MODEL_NOTE`].

#![forbid(unsafe_code)]

pub mod baum_welch;
pub mod em;
pub(crate) mod json;
pub mod json_read;
pub mod kmeans;

use std::collections::BTreeMap;

use macs_core::{Coord, Genome};
use macs_rle::SignalTrack;

/// Marker for the hmmlearn convention note; see the module docs.
pub const HMM_MODEL_NOTE: &str = "hmmlearn 0.3.3";

/// `MACS3.Signal.Prob.pnorm2`: the density of `N(u, v)` at `x`, in `f32`.
///
/// ```python
/// ret = 1.0/sqrt(6.283185307179586 * v) * exp(-(x-u)**2 / (2.0 * v))
/// ```
///
/// The constant is written out rather than taken from `PI` so the rounding is
/// the literal one.
///
/// **`x - u` is a float32 subtraction** (F230). The annotated
/// `def pnorm2(x: float32_t, u: float32_t, v: float32_t)` compiles the operands to
/// C `float`s, so the difference is *rounded to f32 before it is squared*.
/// Widening first and subtracting in double (`f64::from(x) - f64::from(u)`) is
/// algebraically equivalent and differs on roughly a quarter of inputs.
///
/// This is verified against the compiled oracle rather than reasoned about:
/// calling `MACS3.Signal.Prob.pnorm2` directly for all 2208
/// `(fraglen, class)` pairs of the yeast500k fixture and comparing raw f32 bit
/// patterns gives **2208/2208** with the f32 subtraction and 2081/2208 without.
/// The compiled module is the arbiter here, because the same expression read as
/// source is ambiguous -- `pow`, `exp`/`expf`, `sqrt`/`sqrtf`, and float-vs-double
/// intermediates were each tried and none of them, alone, accounts for the result.
#[inline]
pub fn pnorm2(x: f32, u: f32, v: f32) -> f32 {
    // `__pyx_t_9 = (-powf((__pyx_v_x - __pyx_v_u), 2.0));` -- Prob.c, line 14430.
    // The subtraction *and* the power are single precision; the square is therefore
    // rounded to f32 before it is widened for the division. Widening first and
    // squaring in double is algebraically the same and is not bit-identical.
    let sq = -(x - u).powf(2.0f32);
    let v64 = f64::from(v);
    // `__pyx_t_10 = (2.0 * __pyx_v_v);` is a double, and so is the `exp` argument.
    let e = (f64::from(sq) / (2.0 * v64)).exp();
    // `sqrt` is Python's `math.sqrt` on the double `6.28.. * v` (Prob.py imports it
    // from math), and `1.0 / that` is a Python float division -- both double.
    #[allow(clippy::approx_constant)]
    const TWO_PI: f64 = 6.283_185_307_179_586;
    let inv = 1.0 / (TWO_PI * v64).sqrt();
    // the product is a Python float multiply, rounded once on the float32 return
    (inv * e) as f32
}

/// Which nucleosomal class a weight belongs to.
pub const N_SIGNALS: usize = 4;

/// The four weight maps, one per class: short, mono, di, tri.
///
/// Keys are fragment lengths (`R - L`), values are weights in `f32`.
/// Keys are fragment lengths (`R - L`), values are weights in `f32`.
///
/// **Open (F230).** The weight is stored in a plain Python `dict`, so upstream's
/// value is nominally a *double*; but all 2208 weights on yeast500k come back
/// exactly f32-representable, so it is in practice rounded to f32 somewhere.
/// Neither model reproduces all of them:
///
/// | `s` summation | quotient | weights reproduced (of 2140 non-zero) |
/// |---|---|---|
/// | f32, step-rounded | f32 | **1085** |
/// | f32, one rounding at the end | f32 | 1024 |
/// | f64 | f32 | 858 |
/// | any | f64 | 205, and *only* class 0 |
///
/// The last row is the informative one: an f64 quotient is exactly right for 205
/// class-0 weights and for **zero** of the 1078 in classes 1-3, so the quotient is
/// not uniformly f64. The residual concentrates in the tiny classes, where the
/// component probabilities are denormal f32 (`p_tri ~ 4e-40`), which is where
/// double-rounding and denormal handling diverge from whatever the compiled code
/// does. Until that is pinned down, the best-measured model is kept rather than a
/// prettier one: f32 step-rounded sum, f32 quotient.
pub type WeightMapping = [BTreeMap<Coord, f32>; N_SIGNALS];

/// `generate_weight_mapping` (`HMMR_Signal_Processing.py:53`).
///
/// Each fragment length gets the four Gaussian densities normalised to sum to 1.
/// A length whose four densities are all below `min_frag_p` is excluded with
/// weight 0 in every class -- upstream's comment says "Normally this fragment
/// is too large".
///
/// `means` and `stddevs` are the post-EM `short/mono/di/tri` values; the
/// densities use the **squared** standard deviations as variances.
pub fn generate_weight_mapping(
    fraglen_list: &[Coord],
    means: &[f32; N_SIGNALS],
    stddevs: &[f32; N_SIGNALS],
    min_frag_p: f32,
) -> WeightMapping {
    let var: [f32; N_SIGNALS] = stddevs.map(|s| s * s);
    let mut ret: WeightMapping = Default::default();
    for &fl in fraglen_list {
        let x = fl as f32;
        let p: [f32; N_SIGNALS] = std::array::from_fn(|i| pnorm2(x, means[i], var[i]));
        let s = p[0] + p[1] + p[2] + p[3];
        if p.iter().all(|v| *v < min_frag_p) {
            for m in ret.iter_mut() {
                m.insert(fl, 0.0);
            }
            continue;
        }
        // see `WeightMapping` for why this stays f32 despite the double-valued dict
        for i in 0..N_SIGNALS {
            ret[i].insert(fl, p[i] / s);
        }
    }
    ret
}

/// `pileup_from_LR_hmmratac` (`PileupV2.py:680`).
///
/// Each fragment contributes `+weight` at its left end and `-weight` at its
/// right end; the resulting signed endpoint list is sorted and swept. A
/// weight of 0 produces two entries that cancel, matching upstream's
/// dictionary lookup returning `0.0`.
pub fn pileup_from_lr_hmmratac(
    lr: &[(Coord, Coord)],
    mapping: &BTreeMap<Coord, f32>,
) -> Vec<(Coord, f32)> {
    // (position, delta) with the `i4`/`f4` promotion of upstream's records
    let mut pv: Vec<(Coord, f32)> = Vec::with_capacity(2 * lr.len());
    for &(l, r) in lr {
        let w = mapping.get(&(r - l)).copied().unwrap_or(0.0);
        pv.push((l, w));
        pv.push((r, -w));
    }
    // F225: the sort is **lexicographic on `(p, v)`**, not on `p` alone.
    //
    // Upstream stores the endpoint list as a structured array
    // `dtype=[('p','u4'), ('v','f4')]` and calls `.sort(order=['p'])`. NumPy's
    // comparison for a structured dtype uses the listed fields first and then the
    // *remaining fields in dtype order*, so records sharing a `p` are ordered by `v`
    // ascending. Verified directly against numpy:
    //
    //     p: [3,1,1,2,1,2,3,1]  v: [5,9,2,7,1,0,3,4]
    //     after sort(order=['p']):  p [1,1,1,1,2,2,3,3]
    //                               v [1,2,4,9,0,7,3,5]
    //
    // That makes the order **fully determined** -- it is not numpy's introsort
    // pivoting, which is what an earlier note (F206) assumed. Getting it wrong changes
    // the order of `z += v` within a position, and float32 addition is not associative,
    // so the last bit of every digested track differs. That was the entire residual
    // G13 gap: Jaccard 0.973 against a 0.98 requirement, with the signals otherwise
    // matching to 1e-5.
    //
    // `total_cmp` rather than `partial_cmp`: -0.0 and 0.0 compare equal numerically (so
    // their order cannot change a sum) but `total_cmp` orders them deterministically
    // instead of leaving it to the sort's internals.
    pv.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.total_cmp(&b.1)));

    // `pileup_PV` (`PileupV2.py:653-672`), transcribed. Two behaviours matter and both
    // were missing:
    //
    // 1. A record is `(e, z)`: the position at which the value changes and the value
    //    in force *before* it, so `z` describes `[previous e, e)`.
    // 2. **Runs of equal value are merged**: when the incoming value `z` equals the
    //    previously emitted one (`pre_z`), the previous record's end position is
    //    extended instead of a new record being appended. This is why upstream's
    //    `a_digested_*.bdg` contains single records spanning hundreds of bins
    //    (`chr1 510 1680 0.00000`) while an unmerged sweep emits one per bin. The
    //    emission is also what `--save-likelihoods` feeds into `bedGraphIO`, so
    //    without it those files are structurally different, not merely numerically.
    let mut out_pos: Vec<Coord> = Vec::with_capacity(pv.len());
    let mut out_val: Vec<f32> = Vec::with_capacity(pv.len());
    let mut c = 0usize;
    let mut s: Coord = 0;
    let mut z = 0.0f32;
    // `pre_z = -10000` upstream, not 0. It matters whenever the first endpoint is
    // at a position > 0: `z` is still 0 there, so a `pre_z` of 0 would compare
    // equal and the leading zero-depth run would be silently dropped instead of
    // emitted (the `c > 0` guard below then makes that a no-op). With the sentinel
    // the first record is always emitted, which is why upstream's digested tracks
    // always begin at position 0.
    let mut pre_z = -10000.0f32;
    for &(e, v) in &pv {
        if e != s {
            if z == pre_z {
                if c > 0 {
                    out_pos[c - 1] = e;
                }
            } else {
                if c < out_pos.len() {
                    out_pos[c] = e;
                    out_val[c] = z;
                } else {
                    out_pos.push(e);
                    out_val.push(z);
                }
                c += 1;
                pre_z = z;
            }
        }
        z += v;
        s = e;
    }
    out_pos.truncate(c);
    out_val.truncate(c);
    out_pos.into_iter().zip(out_val).collect()
}

/// Per-chromosome digested pileups: `[class][chrom] -> depth profile`.
pub fn pileup_bdg_hmmratac(
    genome: &Genome,
    frags: &BTreeMap<macs_core::ChromId, Vec<(Coord, Coord)>>,
    mapping: &WeightMapping,
) -> Vec<BTreeMap<macs_core::ChromId, SignalTrack<f32>>> {
    let mut chroms: Vec<macs_core::ChromId> = frags.keys().copied().collect();
    chroms.sort_by(|&a, &b| genome.name(a).cmp(genome.name(b)));
    let mut out: Vec<BTreeMap<macs_core::ChromId, SignalTrack<f32>>> =
        (0..N_SIGNALS).map(|_| BTreeMap::new()).collect();
    for (i, m) in mapping.iter().enumerate() {
        for &chrom in &chroms {
            let runs = pileup_from_lr_hmmratac(&frags[&chrom], m);
            let end = frags[&chrom].last().map_or(0, |f| f.1);
            let mut t = SignalTrack::empty(chrom, 0, end);
            // Each returned run is `(end, value)` with the value in force on
            // `[previous end, end)`, exactly as `pileup_PV` stores it -- so consecutive
            // equal values are already merged and a zero-depth gap never appears.
            // Zero-depth runs are emitted too, not skipped: they are what makes
            // upstream's digested bedGraph read `chr1 510 1680 0.00000` between two
            // clusters. `bedGraphIO` prints the value regardless, so dropping them here
            // would silently splice two unrelated runs together.
            let mut prev: Coord = 0;
            for (e, v) in runs {
                if e > prev {
                    t.push(e, v);
                }
                prev = e;
            }
            out[i].insert(chrom, t);
        }
    }
    out
}

/// The training data extracted from peak regions for the HMM.
///
/// Mirrors `extract_signals_from_regions`: one row per bin, four feature
/// columns, and a sequence length per region.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TrainingData {
    /// Bin start position per row.
    pub bins: Vec<Coord>,
    /// `[short, mono, di, tri]` per row, floored at `0.0001`.
    ///
    /// f64, not f32: upstream is `max(0.0001, extracted_data[k][i])`
    /// (`HMMR_Signal_Processing.py:197-200`) with a Python float literal `0.0001`,
    /// so the comparison happens in f64 and the floored entries come out as exactly
    /// `0.0001`, not `f32(0.0001) = 0.00009999999747378752`. Storing f32 made every
    /// floored cell print one ulp low in `*_training_data.txt` (54174 lines here).
    pub rows: Vec<[f64; N_SIGNALS]>,
    /// Length of each region's run of consecutive bins.
    pub lengths: Vec<usize>,
}

/// Bins covering `regions`, as `(chrom, end, mark_bin)`.
///
/// A port of `_make_bdg_of_bins_from_regions`
/// (`HMMR_Signal_Processing.py:232`): each region's `[start, end)` is floored
/// to a bin boundary and emitted one bin at a time, with gaps filled at value 0
/// and `mark_bin` incrementing per region. The mark value doubles as the group
/// identifier, which is how `extract_value_hmmr` recovers region membership.
pub fn make_bdg_of_bins_from_regions(
    genome: &Genome,
    regions: &BTreeMap<macs_core::ChromId, Vec<(Coord, Coord)>>,
    binsize: Coord,
) -> Vec<(macs_core::ChromId, Coord, i32)> {
    let mut out = Vec::new();
    let mut chroms: Vec<macs_core::ChromId> = regions.keys().copied().collect();
    chroms.sort_by(|&a, &b| genome.name(a).cmp(genome.name(b)));
    let mut mark_bin: i32 = 1;
    for &chrom in &chroms {
        let mut tmp_p: Coord = 0;
        for &(s, e) in &regions[&chrom] {
            // "make bins, no need to be too accurate"
            let s = s / binsize * binsize;
            let e = e / binsize * binsize;
            let mut r = s;
            while r < e {
                if r > tmp_p {
                    out.push((chrom, r, 0));
                }
                out.push((chrom, r + binsize, mark_bin));
                tmp_p = r + binsize;
                r += binsize;
            }
            mark_bin += 1;
        }
    }
    out
}

/// One extracted bin: its chromosome, position, the four signal values and the
/// region mark it belongs to.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtractedBin {
    /// Chromosome id.
    pub chrom: macs_core::ChromId,
    /// Bin start.
    pub pos: Coord,
    /// `[short, mono, di, tri]` values.
    pub values: [f32; N_SIGNALS],
    /// Region mark from the bin bedGraph.
    pub mark: i32,
}

/// `extract_value_hmmr` (`BedGraph.py:1074`) over four signal tracks.
///
/// The walk emits a row only when the **bin** track's value is non-zero (i.e.
/// the position is inside a peak region) and the signal value comes from
/// whichever signal pointer is behind -- the same quirk as
/// [`macs_bedgraph::BedGraph::extract_value`], but keyed on the mark rather
/// than `v2 > 0`.
pub fn extract_value_hmmr(
    signals: &[&SignalTrack<f32>; N_SIGNALS],
    bins: &[(macs_core::ChromId, Coord, i32)],
) -> Vec<ExtractedBin> {
    let mut out = Vec::new();
    if bins.is_empty() {
        return out;
    }
    // run-length cursors, one per signal, advanced only while they are behind
    let mut cursors: [(&[macs_rle::Run<f32>], usize); N_SIGNALS] =
        std::array::from_fn(|i| (signals[i].runs(), 0));
    for &(chrom, bp2, mark) in bins {
        // F232: the cursor must stop on the first run whose **end** is `>= bp2`, not
        // on the run that *contains* `bp2`.
        //
        // `extract_value_hmmr` walks the signal track by `get_data_by_chr`, whose
        // positions are run ends (`add_chrom_data_PV` is fed the `(end, value)` pairs
        // `pileup_PV` produces), and it advances only `while p1 < p2`. So the value
        // it reports for a bin at `p2` is the one in force *before* that bin. With
        // chrXVI's runs `[0,8)=0.83172, [8,10)=0.85724, [10,11)=1.75845, ...`:
        //
        //     bin at 10 -> upstream 0.85724 (run [8,10)),  containing-run lookup gives 1.75845
        //     bin at 20 -> upstream 11.15583 (run [18,20)), containing-run lookup gives 12.32470
        //
        // Every one of the four signals has to be advanced before a row is emitted,
        // and the moment any of them runs out the *entire* walk stops -- upstream lets
        // `StopIteration` escape the `while True`, silently dropping that bin and every
        // later one on the chromosome. That truncation is why upstream ends up with
        // ~99 fewer bins per chromosome than a naive implementation, and it is not a
        // bug to be compensated for.
        let mut exhausted = false;
        for c in cursors.iter_mut() {
            while c.1 < c.0.len() && c.0[c.1].end < bp2 {
                c.1 += 1;
            }
            if c.1 >= c.0.len() {
                exhausted = true;
                break;
            }
        }
        if exhausted {
            break;
        }
        if mark != 0 {
            let values = std::array::from_fn(|i| cursors[i].0[cursors[i].1].value);
            out.push(ExtractedBin {
                chrom,
                pos: bp2,
                values,
                mark,
            });
        }
    }
    out
}

/// Floor applied to extracted signal values, as upstream does before feeding
/// the HMM. f64 because upstream's literal is a Python float.
const MIN_VALUE: f64 = 0.0001;

/// Turn extracted bins into HMM training rows.
///
/// `poisson` casts each value to an integer first, since the Poisson emission
/// only accepts counts. Region membership comes from `mark` changing.
pub fn extract_signals_from_regions(bins: &[ExtractedBin], poisson: bool) -> TrainingData {
    let mut td = TrainingData::default();
    let mut counter = 0usize;
    let mut prev_c = 0i32;
    for b in bins {
        td.bins.push(b.pos);
        let v = std::array::from_fn(|i| f64::from(b.values[i]).max(MIN_VALUE));
        td.rows.push(if poisson {
            v.map(|x| x as i64 as f64)
        } else {
            v
        });
        if counter != 0 && b.mark != prev_c {
            td.lengths.push(counter);
            counter = 0;
        }
        prev_c = b.mark;
        counter += 1;
    }
    if counter != 0 {
        td.lengths.push(counter);
    }
    td
}

/// A Gaussian HMM with full or diagonal covariance.
#[derive(Debug, Clone, PartialEq)]
pub struct GaussianHMM {
    /// Initial state probabilities.
    pub startprob: Vec<f64>,
    /// Transition matrix.
    pub transmat: Vec<Vec<f64>>,
    /// Per-state mean vector.
    pub means: Vec<Vec<f64>>,
    /// Per-state covariance: `n x n` for `"full"`, `n` for `"diag"`.
    pub covars: Covars,
    /// `"full"` or `"diag"`.
    pub covariance_type: String,
}

/// Covariance storage, matching upstream's JSON shape.
#[derive(Debug, Clone, PartialEq)]
pub enum Covars {
    /// One `n x n` matrix per state.
    Full(Vec<Vec<Vec<f64>>>),
    /// One length-`n` vector per state.
    Diag(Vec<Vec<f64>>),
}

/// Which family of emissions the model uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HmmType {
    /// Normal emissions over a continuous 4-vector.
    Gaussian,
    /// Poisson emissions over counts.
    Poisson,
}

impl GaussianHMM {
    /// Number of states.
    pub fn n_states(&self) -> usize {
        self.startprob.len()
    }

    /// Log emission probability of observation `x` in state `k`.
    ///
    /// For `"full"` this is the log of a multivariate normal with the given
    /// mean and covariance; for `"diag"` it factorises.
    pub fn log_emit(&self, k: usize, x: &[f64]) -> f64 {
        let m = &self.means[k];
        let n = m.len();
        match &self.covars {
            // `_log_multivariate_normal_density_diag`: the variance is floored at
            // `np.finfo(float).tiny` so a degenerate covariance yields 0*log(0)=0
            // rather than a NaN, and the three terms are summed as
            // `nf*log(2pi) + sum(log(var)) + sum(d^2/var)`.
            Covars::Diag(v) => {
                let c = &v[k];
                let mut log_vars = 0.0;
                let mut quad = 0.0;
                for i in 0..n {
                    let var = c[i].max(f64::MIN_POSITIVE);
                    log_vars += var.ln();
                    let d = x[i] - m[i];
                    quad += d * d / var;
                }
                -0.5 * ((n as f64) * (2.0 * std::f64::consts::PI).ln() + log_vars + quad)
            }
            // F203: hmmlearn evaluates `full`/`tied` Gaussian emissions through a
            // **Cholesky factor** (`_emissions.py`/`stats.py`
            // `_log_multivariate_normal_density_full`), not a Gaussian solve:
            //
            //     cv_chol   = linalg.cholesky(cv, lower=True)
            //     cv_log_det = 2 * np.sum(np.log(np.diagonal(cv_chol)))
            //     cv_sol     = linalg.solve_triangular(cv_chol, (X - mu).T, lower=True).T
            //     log_prob   = -0.5 * (nf*log(2*pi) + (cv_sol**2).sum(axis=1) + cv_log_det)
            //
            // This port instead solved `c z = d` by Gaussian elimination and formed
            // `d.z` plus `log_det(c)`. The two agree mathematically and disagree
            // numerically, and Viterbi is a sequence of argmaxes: where two paths are
            // within rounding of each other the elimination order silently picks the
            // other one. The visible effect was 13 of 946 accessible regions missing
            // (Jaccard 0.958 against a 0.98 requirement) and 1200 of 2491 state rows
            // taking a different path.
            //
            // Transcribing the Cholesky form is what makes the decoded state path -- and
            // therefore the accessible regions -- identical.
            Covars::Full(v) => {
                let c = &v[k];
                let chol = match cholesky_lower(c) {
                    Some(l) => l,
                    // Same fallback order as hmmlearn: retry with `min_covar` added to
                    // the diagonal before giving up.
                    None => {
                        let mut jittered = c.clone();
                        for (i, row) in jittered.iter_mut().enumerate().take(n) {
                            row[i] += MIN_COVAR;
                        }
                        match cholesky_lower(&jittered) {
                            Some(l) => l,
                            // `validate()` already rejects this at load time, so this is
                            // only reachable through a numerically degenerate factor.
                            // A zero-probability emission keeps the forward pass total.
                            None => return f64::NEG_INFINITY,
                        }
                    }
                };
                // `2 * np.sum(np.log(np.diagonal(cv_chol)))`. Indexing by `i` rather
                // than a computed offset matters: reading the diagonal through
                // `row[n.min(row.len()) - 1]` silently produced `-inf` here, which made
                // every emission `+inf` and collapsed the decoded state path to nothing.
                let mut log_det = 0.0;
                for (i, row) in chol.iter().enumerate().take(n) {
                    log_det += row[i].ln();
                }
                log_det *= 2.0;
                // forward substitution: L y = (x - mu)
                let mut y = vec![0.0f64; n];
                for i in 0..n {
                    let mut acc = x[i] - m[i];
                    for j in 0..i {
                        acc -= chol[i][j] * y[j];
                    }
                    y[i] = acc / chol[i][i];
                }
                let mut sq = 0.0;
                for v in &y {
                    sq += v * v;
                }
                -0.5 * ((n as f64) * (2.0 * std::f64::consts::PI).ln() + sq + log_det)
            }
        }
    }

    /// Forward-backward posteriors, one row per observation.
    ///
    /// Uses the log-space recursion (`hmmlearn.base._score_log`), which is what
    /// `predict_proba` calls; the scaled variant raises on underflow for these
    /// models because a trained ATAC HMM drives whole states to ~1e-88.
    ///
    /// The transition orientation is hmmlearn's own: `A[i][j]`, i.e. the matrix
    /// is applied transposed relative to the textbook form. That is what MACS3
    /// inherits, so it is reproduced rather than corrected -- flipping it moves
    /// every summit.
    pub fn posterior(&self, obs: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let t = obs.len();
        let k = self.n_states();
        if t == 0 || k == 0 {
            return Vec::new();
        }
        let mut log_b = vec![vec![f64::NEG_INFINITY; t]; k];
        for (ti, x) in obs.iter().enumerate() {
            for (s, col) in log_b.iter_mut().enumerate() {
                col[ti] = self.log_emit(s, x);
            }
        }
        posterior_from_log_lattices(&self.startprob, &self.transmat, &log_b, t, k)
    }
}

/// A Poisson HMM over count emissions.
#[derive(Debug, Clone, PartialEq)]
pub struct PoissonHMM {
    /// Initial state probabilities.
    pub startprob: Vec<f64>,
    /// Transition matrix.
    pub transmat: Vec<Vec<f64>>,
    /// Per-state rate per feature.
    pub lambdas: Vec<Vec<f64>>,
}

impl PoissonHMM {
    /// Number of states.
    pub fn n_states(&self) -> usize {
        self.startprob.len()
    }

    /// Log emission probability of observation `x` in state `k`.
    pub fn log_emit(&self, k: usize, x: &[f64]) -> f64 {
        let lam = &self.lambdas[k];
        let mut s = 0.0;
        for (i, &xi) in x.iter().enumerate() {
            let l = lam[i];
            let c = xi.trunc();
            s += c * l.ln() - l - ln_gamma(c + 1.0);
        }
        s
    }

    /// Forward-backward posteriors; same convention as [`GaussianHMM::posterior`].
    pub fn posterior(&self, obs: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let t = obs.len();
        let k = self.n_states();
        if t == 0 || k == 0 {
            return Vec::new();
        }
        let mut log_b = vec![vec![f64::NEG_INFINITY; t]; k];
        for (ti, x) in obs.iter().enumerate() {
            for (s, col) in log_b.iter_mut().enumerate() {
                col[ti] = self.log_emit(s, x);
            }
        }
        posterior_from_log_lattices(&self.startprob, &self.transmat, &log_b, t, k)
    }
}

/// `logsumexp` over an iterator, ignoring `-inf` terms.
fn log_sum_exp(xs: impl DoubleEndedIterator<Item = f64> + Clone) -> f64 {
    let mut m = f64::NEG_INFINITY;
    for x in xs.clone() {
        if x > m {
            m = x;
        }
    }
    if m == f64::NEG_INFINITY {
        return m;
    }
    let mut acc = 0.0;
    for x in xs {
        acc += (x - m).exp();
    }
    m + acc.ln()
}

/// Forward-backward posteriors in log space, as `hmmlearn.base._score_log` does.
///
/// `alpha[0][i] = log startprob[i] + log B[0][i]`, then
/// `alpha[t][i] = logsumexp_j (alpha[t-1][j] + log A[i][j]) + log B[t][i]`, and
/// symmetrically for `beta`. The posterior is `exp` of the row-normalised sum,
/// matching `_compute_posteriors_log`.
///
/// The matrix is indexed `A[i][j]` throughout, i.e. transposed relative to the
/// textbook `sum_j A[j][i]` form -- that is hmmlearn's convention.
fn posterior_from_log_lattices(
    startprob: &[f64],
    transmat: &[Vec<f64>],
    log_b: &[Vec<f64>],
    t: usize,
    k: usize,
) -> Vec<Vec<f64>> {
    let mut alpha = vec![vec![f64::NEG_INFINITY; t]; k];
    let mut beta = vec![vec![f64::NEG_INFINITY; t]; k];

    for s in 0..k {
        alpha[s][0] = startprob[s].ln() + log_b[s][0];
    }
    for i in 1..t {
        for s in 0..k {
            // textbook orientation: sum_j alpha[j] + log A[j][s]
            let acc = log_sum_exp((0..k).map(|j| alpha[j][i - 1] + transmat[j][s].ln()));
            alpha[s][i] = acc + log_b[s][i];
        }
    }

    for row in beta.iter_mut() {
        row[t - 1] = 0.0;
    }
    for i in (0..t.saturating_sub(1)).rev() {
        for j in 0..k {
            // hmmlearn's `backward_log` indexes the transition matrix as
            // `A[j][i]` -- the transpose of the textbook `A[i][j]`, and the
            // *opposite* of its own `forward_log`. That inconsistency is real and
            // MACS3 inherits it; using the textbook form here moves the
            // posterior by ~2% on non-dominant states.
            beta[j][i] =
                log_sum_exp((0..k).map(|s| transmat[j][s].ln() + log_b[s][i + 1] + beta[s][i + 1]));
        }
    }

    let mut out = vec![vec![0.0f64; k]; t];
    for (i, row) in out.iter_mut().enumerate() {
        let z = log_sum_exp((0..k).map(|s| alpha[s][i] + beta[s][i]));
        for s in 0..k {
            row[s] = (alpha[s][i] + beta[s][i] - z).exp();
        }
    }
    out
}

/// Solve `a x = b` for a small dense system by Gaussian elimination with
/// partial pivoting.
/// hmmlearn's `min_covar` default for the Cholesky retry (`stats.py:63`).
pub const MIN_COVAR: f64 = 1.0e-7;

/// `scipy.linalg.cholesky(a, lower=True)`: the lower-triangular `L` with `L @ L.T = a`.
///
/// Returns `None` when `a` is not numerically positive-definite, which is the same
/// condition under which SciPy raises `LinAlgError`.
fn cholesky_lower(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut l = vec![vec![0.0f64; n]; n];
    // The explicit indices are the Cholesky recurrence itself
    // (`L[i][j] = (a[i][j] - sum_k L[i][k]L[j][k]) / L[j][j]`); rewriting it as an
    // iterator walk would need a split borrow of `l` that the borrow checker rejects.
    #[allow(clippy::needless_range_loop)]
    for i in 0..n {
        #[allow(clippy::needless_range_loop)]
        for j in 0..=i {
            let mut acc = a[i][j];
            #[allow(clippy::needless_range_loop)]
            for k in 0..j {
                acc -= l[i][k] * l[j][k];
            }
            if i == j {
                // `!(acc > 0.0)` rather than `acc <= 0.0` so a NaN pivot is also
                // rejected, matching SciPy's LinAlgError on a non-PD matrix.
                #[allow(clippy::neg_cmp_op_on_partial_ord)]
                if !(acc > 0.0) || !acc.is_finite() {
                    return None;
                }
                l[i][j] = acc.sqrt();
            } else {
                l[i][j] = acc / l[j][j];
            }
        }
    }
    Some(l)
}

/// `ln det` of a small dense matrix, by Gaussian elimination with partial
/// pivoting (the determinant is the product of the pivots).
#[allow(dead_code)]
fn log_det(a: &[Vec<f64>]) -> f64 {
    let n = a.len();
    let mut m: Vec<Vec<f64>> = a.to_vec();
    let mut acc = 0.0f64;
    for col in 0..n {
        let mut piv = col;
        for (r, row) in m.iter().enumerate().skip(col) {
            if row[col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        let head = m[piv][col];
        if head.abs() < 1e-300 {
            return f64::NEG_INFINITY;
        }
        m.swap(col, piv);
        // `|det|` is what a density needs, so the sign is taken out here rather
        // than tracked as a swap parity. Taking `ln` of a signed pivot directly
        // would yield NaN whenever elimination produces a negative one.
        acc += m[col][col].abs().ln();
        let pivot = m[col].clone();
        for row in m.iter_mut().skip(col + 1) {
            let f = row[col] / pivot[col];
            for (c, cell) in row.iter_mut().enumerate().skip(col) {
                *cell -= f * pivot[c];
            }
        }
    }
    acc
}

/// Natural log of the gamma function, by Stirling with a Lanczos-style
/// correction; `lgamma` is not in Rust's standard library.
pub fn ln_gamma(x: f64) -> f64 {
    // Lanczos approximation, g = 7, n = 9
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // reflection: lgamma(x) = log(pi/sin(pi x)) - lgamma(1-x)
        (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = C[0];
        let t = x + G + 0.5;
        for (i, c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

// ---------------------------------------------------------------------------
// Model serialisation
// ---------------------------------------------------------------------------

/// A decoded `*_model.json`, as `hmm_model_init` (`HMMR_HMM.py:118`) reads it.
///
/// Upstream only reads a subset -- the covariance/lambdas are attached to an
/// otherwise-unfitted estimator object -- so the whole file is kept and the
/// parts inference needs are exposed.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelFile {
    /// Which emission family.
    pub hmm_type: HmmType,
    /// Initial state probabilities.
    pub startprob: Vec<f64>,
    /// Transition matrix.
    pub transmat: Vec<Vec<f64>>,
    /// Per-state mean vectors (Gaussian only; absent in a Poisson model).
    pub means: Vec<Vec<f64>>,
    /// Per-state covariance (Gaussian only; absent in a Poisson model).
    pub covars: Covars,
    /// `covariance_type`, `"diag"` or `"full"`. Kept separate because upstream
    /// writes both and reads both.
    pub covariance_type: String,
    /// Per-state rates (Poisson only; absent in a Gaussian model).
    pub lambdas: Vec<Vec<f64>>,
    /// Index of the "open" state.
    pub i_open_region: usize,
    /// Index of the "background" state.
    pub i_background_region: usize,
    /// Index of the "nucleosomal" state.
    pub i_nucleosomal_region: usize,
    /// Bin size the model was trained at.
    pub hmm_binsize: Coord,
    /// `n_features`, as written by upstream. Always `4` in practice; kept so the
    /// loaded file round-trips byte-for-byte.
    pub n_features: usize,
}

impl Default for Covars {
    fn default() -> Self {
        Covars::Diag(Vec::new())
    }
}

impl ModelFile {
    /// `hmm_model_init` (`HMMR_HMM.py:118`): read and validate a `*_model.json`.
    ///
    /// The `hmm_type` default is `gaussian` -- `m.get("hmm_type", "gaussian")` --
    /// so a file written before the field existed loads as Gaussian.
    pub fn load(path: &std::path::Path) -> macs_core::Result<Self> {
        let text = std::fs::read_to_string(path)?;
        Self::from_json(&text)
    }

    /// Parse a model from JSON text.
    ///
    /// Uses [`json_read`] rather than `serde_json` so every number is converted
    /// with a correctly-rounded `f64::from_str`; see that module for why.
    pub fn from_json(text: &str) -> macs_core::Result<Self> {
        use json_read::{cube, matrix, vector, Json};
        let bad = |m: &str| macs_core::MacsError::InvalidParameter(format!("bad model JSON: {m}"));
        let j = json_read::parse(text).map_err(|e| bad(&e.to_string()))?;
        fn need<'a>(
            j: &'a json_read::Json,
            k: &str,
        ) -> Result<&'a json_read::Json, macs_core::MacsError> {
            j.get(k).ok_or_else(|| {
                macs_core::MacsError::InvalidParameter(format!("bad model JSON: missing `{k}`"))
            })
        }
        let opt = |k: &str| j.get(k);

        let hmm_type = match opt("hmm_type").and_then(Json::as_str) {
            // `m.get("hmm_type", "gaussian")`: a file written before the field
            // existed is Gaussian.
            None | Some("gaussian") => HmmType::Gaussian,
            Some("poisson") => HmmType::Poisson,
            Some(other) => return Err(bad(&format!("unknown hmm_type `{other}`"))),
        };

        let startprob = vector(need(&j, "startprob")?).ok_or_else(|| bad("bad startprob"))?;
        let transmat = matrix(need(&j, "transmat")?).ok_or_else(|| bad("bad transmat"))?;

        let (means, lambdas, covars, covariance_type) = match hmm_type {
            HmmType::Gaussian => {
                let means = opt("means")
                    .and_then(matrix)
                    .ok_or_else(|| bad("missing or bad means"))?;
                let cw = opt("covars").ok_or_else(|| bad("missing covars"))?;
                // Depth, not the sibling `covariance_type`, identifies the shape:
                // `3 x n` is diag, `3 x n x n` is full.
                let covars = if let Some(f) = cube(cw) {
                    Covars::Full(f)
                } else if let Some(d) = matrix(cw) {
                    Covars::Diag(d)
                } else {
                    return Err(bad("bad covars"));
                };
                let covariance_type = match opt("covariance_type").and_then(Json::as_str) {
                    Some(s) => s.to_string(),
                    None => match &covars {
                        Covars::Full(_) => "full".into(),
                        Covars::Diag(_) => "diag".into(),
                    },
                };
                (means, Vec::new(), covars, covariance_type)
            }
            HmmType::Poisson => {
                let lambdas = opt("lambdas")
                    .and_then(matrix)
                    .ok_or_else(|| bad("missing or bad lambdas"))?;
                (Vec::new(), lambdas, Covars::Diag(Vec::new()), "diag".into())
            }
        };

        let g = |k: &str, d: i64| opt(k).and_then(Json::as_int).unwrap_or(d);
        let m = ModelFile {
            hmm_type,
            startprob,
            transmat,
            means,
            covars,
            covariance_type,
            lambdas,
            i_open_region: g("i_open_region", 0).max(0) as usize,
            i_background_region: g("i_background_region", 1).max(0) as usize,
            i_nucleosomal_region: g("i_nucleosomal_region", 2).max(0) as usize,
            hmm_binsize: g("hmm_binsize", 0).max(0) as Coord,
            n_features: opt("n_features").and_then(Json::as_int).unwrap_or(4).max(0) as usize,
        };
        m.validate()?;
        Ok(m)
    }

    /// F202: reject a model hmmlearn would reject, at load time.
    ///
    /// `hmm_model_init` (`HMMR_HMM.py:118-131`) builds a `GaussianHMM` and assigns
    /// `covars_` directly, so the first thing that touches the numbers is hmmlearn's
    /// `_validate_covars`, which raises
    /// `ValueError: component %d of '%s' covars must be symmetric, positive-definite`.
    /// Upstream surfaces that as a clean error; this port used to panic in Gaussian
    /// elimination instead, which is both a wrong exit path and a violation of the
    /// "zero panics on malformed input" criterion.
    ///
    /// This is reachable in practice from upstream's *own* `--modelonly` output: if a
    /// signal feature is constant across the training bins, Baum-Welch produces a
    /// singular covariance and upstream then refuses the model it just wrote.
    fn validate(&self) -> macs_core::Result<()> {
        let bad = |m: String| macs_core::MacsError::InvalidParameter(m);
        if self.hmm_type != HmmType::Gaussian {
            return Ok(());
        }
        if let Covars::Full(v) = &self.covars {
            for (k, c) in v.iter().enumerate() {
                let n = c.len();
                if c.iter().flatten().any(|x| !x.is_finite()) {
                    return Err(bad(format!(
                        "component {k} of '{}' covars must be symmetric, positive-definite",
                        self.covariance_type
                    )));
                }
                for (i, row) in c.iter().enumerate().take(n) {
                    for (j, &b) in row.iter().enumerate().take(n) {
                        let a = c[j][i];
                        let scale = a.abs().max(b.abs()).max(f64::MIN_POSITIVE);
                        if (a - b).abs() / scale > 1e-8 {
                            return Err(bad(format!(
                                "component {k} of '{}' covars must be symmetric, \
                                 positive-definite",
                                self.covariance_type
                            )));
                        }
                    }
                }
                // hmmlearn factors with Cholesky (`stats.py:81-90`), so Cholesky is the
                // criterion that decides whether it accepts the model -- and it is also
                // what the emission path in F203 now uses, so a model that loads is a
                // model that will not hit the -inf fallback below. The `min_covar`
                // retry is *not* mirrored here: upstream surfaces the error, it does not
                // silently repair the model.
                if cholesky_lower(c).is_none() {
                    return Err(bad(format!(
                        "component {k} of '{}' covars must be symmetric, positive-definite",
                        self.covariance_type
                    )));
                }
            }
        }
        Ok(())
    }

    /// `hmm_model_save` (`HMMR_HMM.py:82`): write the model as JSON.
    ///
    /// Key order matches upstream's `json.dump` literal, and upstream writes no
    /// indentation and no trailing newline, so the file is byte-comparable.
    pub fn to_json(&self) -> String {
        crate::json::write_model(self)
    }

    /// Write the model to `path`, matching `hmm_model_save`'s format.
    pub fn save(&self, path: &std::path::Path) -> macs_core::Result<()> {
        std::fs::write(path, self.to_json())?;
        Ok(())
    }

    /// Build a [`GaussianHMM`] view for inference.
    pub fn gaussian(&self) -> GaussianHMM {
        GaussianHMM {
            startprob: self.startprob.clone(),
            transmat: self.transmat.clone(),
            means: self.means.clone(),
            covars: self.covars.clone(),
            covariance_type: self.covariance_type.clone(),
        }
    }

    /// Build a [`PoissonHMM`] view for inference.
    pub fn poisson(&self) -> PoissonHMM {
        PoissonHMM {
            startprob: self.startprob.clone(),
            transmat: self.transmat.clone(),
            lambdas: self.lambdas.clone(),
        }
    }

    /// The label of state `i`.
    ///
    /// Upstream writes `prob_data[i]` for `i` in `0..3` and reads
    /// `pp_data[i_state + 2]`, so column `i` of a probability row is state `i`.
    pub fn state_label(&self, i: usize) -> &'static str {
        if i == self.i_open_region {
            "open"
        } else if i == self.i_nucleosomal_region {
            "nuc"
        } else {
            "bg"
        }
    }
}

/// One row of the posterior file: a bin end and its per-state posteriors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProbRow {
    /// Bin end coordinate.
    pub end: Coord,
    /// Posterior for each state, in model order.
    pub probs: [f64; 3],
}

/// `generate_states_path` (`hmmratac_cmd.py:575`).
///
/// Assign each bin the label of its most probable state, coalesce runs of the
/// same label, and fill any gap between bins (a chromosome boundary or a
/// skipped bin) with `bg`. The tie-break is Python's stable `max` over
/// `((open,0),(nuc,1),(bg,2))`, so the earliest state index wins a tie.
pub fn generate_states_path(
    rows: &[ProbRow],
    binsize: Coord,
    model: &ModelFile,
) -> Vec<(Coord, Coord, &'static str)> {
    let mut out: Vec<(Coord, Coord, &'static str)> = Vec::new();
    let mut first = true;
    let mut prev_bin_end: Coord = 0;
    let mut prev_label: Option<&'static str> = None;

    for r in rows {
        let start = r.end.saturating_sub(binsize);
        let mut best = 0usize;
        for i in 1..3 {
            if r.probs[i] > r.probs[best] {
                best = i;
            }
        }
        let label = model.state_label(best);
        if first {
            if start > 0 {
                out.push((0, start, "bg"));
            }
            out.push((start, r.end, label));
            prev_label = Some(label);
            first = false;
        } else {
            if prev_bin_end < start {
                out.push((prev_bin_end, start, "bg"));
                prev_label = Some("bg");
            }
            if Some(label) == prev_label {
                let last = out.len() - 1;
                out[last].1 = r.end;
            } else {
                out.push((start, r.end, label));
                prev_label = Some(label);
            }
        }
        prev_bin_end = r.end;
    }
    out
}

/// An accessible region: the `open` run inside a `nuc`-`open`-`nuc` triple.
///
/// Upstream's `add_regions` collects the whole triple and then filters out the
/// `nuc` entries, so what survives is exactly the middle `open` run -- not the
/// triple's span.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AccessibleRegion {
    /// Inclusive start of the open run.
    pub start: Coord,
    /// Inclusive end of the open run.
    pub end: Coord,
}

/// One run of the state path, with the chromosome it belongs to.
///
/// Upstream's `states_path` entries are `(chrom, start, end, label)`; the
/// chromosome is needed to write `*_states.bed` and `*_accessible_regions`, so
/// [`generate_states_path_chromed`] keeps it instead of dropping it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StateRun {
    /// Chromosome the run is on.
    pub chrom: macs_core::ChromId,
    /// Run start.
    pub start: Coord,
    /// Run end.
    pub end: Coord,
    /// `"open"`, `"nuc"` or `"bg"`.
    pub label: &'static str,
}

/// [`generate_states_path`], keeping the chromosome on every run.
///
/// `rows` is the posterior output in emission order; each element is the bin's
/// chromosome paired with its posterior row.
pub fn generate_states_path_chromed(
    rows: &[(macs_core::ChromId, ProbRow)],
    binsize: Coord,
    model: &ModelFile,
) -> Vec<StateRun> {
    let mut out: Vec<StateRun> = Vec::new();
    let mut first = true;
    let mut prev_bin_end: Coord = 0;
    let mut prev_label: Option<&'static str> = None;
    // F233: a chromosome change starts a new run unconditionally.
    //
    // The `else` branch below extends the previous run whenever the label repeats,
    // and it used to do that across chromosomes too. Since the emission order is
    // descending by name and every chromosome starts at bin position 0, the first
    // bin of the next chromosome usually carries the *same* label as the last bin of
    // the previous one -- so the run was extended with a position from a different
    // chromosome:
    //
    //     oracle  chrXVI 947660 948020 open / chrXV 0 710 open
    //     ours    chrXVI 947660   710 open / (chrXV 0-710 open lost)
    //
    // `narrowPeak` happened to survive this, but `*_states.bed` did not.
    let mut prev_chrom: Option<macs_core::ChromId> = None;
    for (chrom, r) in rows {
        let start = r.end.saturating_sub(binsize);
        let mut best = 0usize;
        for i in 1..3 {
            if r.probs[i] > r.probs[best] {
                best = i;
            }
        }
        let label = model.state_label(best);
        let new_chrom = prev_chrom != Some(*chrom);
        if first || new_chrom {
            if start > 0 {
                out.push(StateRun {
                    chrom: *chrom,
                    start: 0,
                    end: start,
                    label: "bg",
                });
            }
            out.push(StateRun {
                chrom: *chrom,
                start,
                end: r.end,
                label,
            });
            prev_label = Some(label);
            prev_chrom = Some(*chrom);
            first = false;
        } else {
            if prev_bin_end < start {
                out.push(StateRun {
                    chrom: *chrom,
                    start: prev_bin_end,
                    end: start,
                    label: "bg",
                });
                prev_label = Some("bg");
            }
            if Some(label) == prev_label {
                let last = out.len() - 1;
                out[last].end = r.end;
            } else {
                out.push(StateRun {
                    chrom: *chrom,
                    start,
                    end: r.end,
                    label,
                });
                prev_label = Some(label);
            }
        }
        prev_bin_end = r.end;
    }
    out
}

/// `save_accessible_regions` (`hmmratac_cmd.py:621`).
///
/// Scans for `nuc`/`open`/`nuc` triples that are contiguous and whose total
/// span exceeds `openregion_minlen`, and returns each triple's `open` run.
pub fn accessible_regions(
    path: &[(Coord, Coord, &'static str)],
    openregion_minlen: Coord,
) -> Vec<AccessibleRegion> {
    let mut regions: Vec<AccessibleRegion> = Vec::new();
    for i in 0..path.len().saturating_sub(2) {
        if path[i].2 == "nuc"
            && path[i + 1].2 == "open"
            && path[i + 2].2 == "nuc"
            && path[i].1 == path[i + 1].0
            && path[i + 1].1 == path[i + 2].0
            && path[i + 2].1.saturating_sub(path[i].1) > openregion_minlen
        {
            regions.push(AccessibleRegion {
                start: path[i + 1].0,
                end: path[i + 1].1,
            });
        }
    }
    regions
}

/// [`accessible_regions`] over [`StateRun`]s, keeping the chromosome.
pub fn accessible_regions_chromed(
    path: &[StateRun],
    openregion_minlen: Coord,
) -> Vec<(macs_core::ChromId, Coord, Coord)> {
    let mut regions: Vec<(macs_core::ChromId, Coord, Coord)> = Vec::new();
    for i in 0..path.len().saturating_sub(2) {
        if path[i].label == "nuc"
            && path[i + 1].label == "open"
            && path[i + 2].label == "nuc"
            // upstream: `states_path[i][2] == states_path[i+1][1]` and
            // `states_path[i+1][2] == states_path[i+2][1]` -- the *end* of one
            // run must equal the *start* of the next, so no gap between them
            && path[i].end == path[i + 1].start
            && path[i + 1].end == path[i + 2].start
            && path[i + 2].end.saturating_sub(path[i].start) > openregion_minlen
        {
            regions.push((path[i + 1].chrom, path[i + 1].start, path[i + 1].end));
        }
    }
    regions
}

#[cfg(test)]
mod tests {
    use super::*;
    use macs_core::ChromId;

    /// (x, u, v, compiled_bits) captured from the **compiled**
    /// `MACS3.Signal.Prob.pnorm2` on the yeast500k fragment lengths.
    ///
    /// (x, u, v, compiled_bits) captured from the **compiled**
    /// `MACS3.Signal.Prob.pnorm2` on the yeast500k fragment lengths.
    ///
    /// These pin F231: the square is `(x-u)**2` evaluated by **`powf` in single
    /// precision** (Prob.c: `__pyx_t_9 = (-powf((__pyx_v_x - __pyx_v_u), 2.0));`),
    /// so both the subtraction and the power round to f32 before the widening
    /// division by `2.0 * v`. Squaring in f64 instead is algebraically identical
    /// and disagrees on roughly a quarter of inputs.
    const PNORM2_ORACLE: &[(f64, f64, f64, u32)] = &[
        (2.0, 166.6999969482422, 2313.60986328125, 0x37c5_ec43),
        (2.0, 50.0, 400.0, 0x3a92_c3ca),
        (3.0, 50.0, 400.0, 0x3aa5_4527),
        (166.0, 166.6999969482422, 2313.60986328125, 0x3c07_dff4),
        (331.0, 331.20001220703125, 1162.8099365234375, 0x3c3f_ad1b),
        (460.0, 460.29998779296875, 1218.0101318359375, 0x3c3b_475c),
        (50.0, 50.0, 400.0, 0x3ca3_6821),
        (554.0, 460.29998779296875, 1218.0101318359375, 0x39a3_1532),
        (100.0, 50.0, 400.0, 0x3a65_bf37),
        (200.0, 166.6999969482422, 2313.60986328125, 0x3bd5_dd47),
        (300.0, 331.20001220703125, 1162.8099365234375, 0x3bfc_3e9e),
        (120.0, 460.29998779296875, 1218.0101318359375, 0x19fa_06bb),
    ];

    #[test]
    fn pnorm2_matches_compiled_oracle_bit_for_bit() {
        for &(x, u, v, want) in PNORM2_ORACLE {
            // the oracle's own call is `pnorm2(float(fl), m, v)`, i.e. doubles
            // narrowed by the compiled signature
            let got = pnorm2(x as f32, u as f32, v as f32).to_bits();
            assert_eq!(
                got, want,
                "pnorm2({x}, {u}, {v}) = {got:#010x}, want {want:#010x}"
            );
        }
    }

    /// The single-precision square is what makes the above hold; the f64 form is
    /// kept here only to show the two are genuinely different functions.
    #[test]
    #[allow(clippy::approx_constant)]
    fn pnorm2_f32_subtraction_is_load_bearing() {
        let mut differing = 0;
        for fl in 1..600u32 {
            for &(u, v) in &[
                (166.6999969482422f64, 27866.79f64),
                (50.0f64, 400.0f64),
                (331.20001220703125, 1162.8),
                (460.29998779296875, 1218.01),
            ] {
                let (x, uu, vv) = (fl as f32, u as f32, v as f32);
                let f32_sub = pnorm2(x, uu, vv);
                let d = f64::from(x) - f64::from(uu);
                let wide = (1.0 / (6.283_185_307_179_586 * f64::from(vv)).sqrt()
                    * (-(d * d) / (2.0 * f64::from(vv))).exp()) as f32;
                if f32_sub.to_bits() != wide.to_bits() {
                    differing += 1;
                }
            }
        }
        assert!(
            differing > 100,
            "only {differing} differences; the models converged?"
        );
    }

    /// F232: a bin's value comes from the run ending at the bin, not the run
    /// covering it -- `extract_value_hmmr` advances `while p1 < p2` over run *end*
    /// positions, so it reports the value in force before the bin.
    #[test]
    fn hmm_extraction_reads_the_run_ending_at_the_bin() {
        let chrom = ChromId(0);
        // runs [0,8)=1.0, [8,10)=2.0, [10,11)=3.0
        let mut track = SignalTrack::<f32>::empty(chrom, 0, 11);
        track.push(8, 1.0);
        track.push(10, 2.0);
        track.push(11, 3.0);
        let sig = [&track; 4];
        // bins at 0, 10 and 20; the last is past the end of the track
        let bins = vec![(chrom, 0u64, 1i32), (chrom, 10, 1), (chrom, 20, 1)];
        let got = extract_value_hmmr(&sig, &bins);
        // bin 0 -> run ending at 8 (value 1.0); bin 10 -> run ending at 10 (value 2.0);
        // bin 20 has no run with end >= 20, so the walk stops and emits nothing.
        assert_eq!(
            got.iter().map(|b| (b.pos, b.values[0])).collect::<Vec<_>>(),
            vec![(0u64, 1.0f32), (10u64, 2.0f32)]
        );
    }

    /// F233: a run must never span two chromosomes, even when the label repeats.
    #[test]
    fn states_path_starts_a_new_run_on_a_new_chromosome() {
        let model = model();
        let mk = |end: u64, probs: [f64; 3]| ProbRow { end, probs };
        let rows = vec![
            (ChromId(7), mk(1000, [0.9, 0.05, 0.05])),
            // next chromosome starts at 0 with the same winning label
            (ChromId(3), mk(10, [0.9, 0.05, 0.05])),
        ];
        let path = generate_states_path_chromed(&rows, 10, &model);
        // bg [0,990) then the winning label, then the new chromosome's own run
        assert_eq!(path.len(), 3, "a chromosome change must not merge runs");
        assert_eq!(
            (path[0].chrom, path[0].start, path[0].end, path[0].label),
            (ChromId(7), 0, 990, "bg")
        );
        assert_eq!(path[1].chrom, ChromId(7));
        assert_eq!(
            path[1].end, 1000,
            "must not be extended to the next chromosome"
        );
        assert_eq!(
            (path[2].chrom, path[2].start, path[2].end),
            (ChromId(3), 0, 10)
        );
    }

    /// `pileup_PV` seeds `pre_z = -10000`, so the first record is always emitted.
    /// With fragments that do not start at position 0 this is observable: the
    /// leading zero-depth run must survive.
    #[test]
    fn pileup_emits_leading_zero_run() {
        // the mapping is keyed by fragment length R - L, which is 100 here
        let mapping: BTreeMap<Coord, f32> = [(100, 1.0f32)].into_iter().collect();
        let runs = pileup_from_lr_hmmratac(&[(100, 200)], &mapping);
        assert_eq!(runs, vec![(100, 0.0), (200, 1.0)]);
    }

    fn model() -> ModelFile {
        ModelFile {
            hmm_type: HmmType::Gaussian,
            startprob: vec![0.33, 0.34, 0.33],
            transmat: vec![
                vec![0.9, 0.05, 0.05],
                vec![0.05, 0.9, 0.05],
                vec![0.05, 0.05, 0.9],
            ],
            means: vec![vec![1.0; 4], vec![2.0; 4], vec![3.0; 4]],
            covars: Covars::Diag(vec![vec![1.0; 4], vec![1.0; 4], vec![1.0; 4]]),
            covariance_type: "diag".into(),
            lambdas: vec![],
            i_open_region: 0,
            i_background_region: 2,
            i_nucleosomal_region: 1,
            hmm_binsize: 10,
            n_features: 4,
        }
    }

    #[test]
    fn pnorm2_integrates_to_one() {
        let h = 0.001f64;
        let mut acc = 0.0f64;
        let mut x = -8.0f64;
        while x < 8.0 {
            acc += pnorm2(x as f32, 0.0, 1.0) as f64 * h;
            x += h;
        }
        assert!((acc - 1.0).abs() < 1e-3, "integral was {acc}");
    }

    #[test]
    fn weights_sum_to_one_over_the_kept_lengths() {
        let lens: Vec<Coord> = (50..=400).collect();
        let wm = generate_weight_mapping(
            &lens,
            &[60.0, 200.0, 380.0, 550.0],
            &[10.0, 25.0, 25.0, 40.0],
            0.001,
        );
        for &fl in &lens {
            let sum: f32 = wm.iter().map(|m| m[&fl]).sum();
            if sum == 0.0 {
                continue;
            }
            assert!((sum - 1.0).abs() < 1e-5, "fl {fl} summed to {sum}");
        }
    }

    #[test]
    fn an_unassignable_length_gets_zero_everywhere() {
        let wm = generate_weight_mapping(&[5000], &[60.0, 200.0, 380.0, 550.0], &[10.0; 4], 0.001);
        for m in wm.iter() {
            assert_eq!(m[&5000], 0.0);
        }
    }

    #[test]
    fn the_hmm_pileup_sweeps_weighted_fragments() {
        // `pileup_PV` returns `(end, value)` with the value in force on
        // `[previous end, end)`. A fragment of length 100 weighted 1.0 opens a run at
        // 10 and closes it at 110. The leading `[0, 10)` zero-depth run is emitted
        // too, because `pre_z` starts at -10000 (F230) rather than 0 -- that is
        // exactly what makes upstream's digested tracks begin at position 0.
        let wm: BTreeMap<Coord, f32> = [(100u64, 1.0f32)].into_iter().collect();
        let d = pileup_from_lr_hmmratac(&[(10, 110)], &wm);
        assert_eq!(d, vec![(10u64, 0.0f32), (110u64, 1.0f32)]);
    }

    /// F204: `pileup_PV` merges consecutive runs of equal value instead of emitting one
    /// record per endpoint. Two adjacent fragments of equal weight therefore produce a
    /// single record spanning both, which is what makes upstream's `a_digested_*.bdg`
    /// contain `chr1 510 1680 0.00000` -- 170 bins as one line.
    #[test]
    fn the_hmm_pileup_merges_adjacent_equal_runs() {
        let wm: BTreeMap<Coord, f32> = [(100u64, 1.0f32)].into_iter().collect();

        // Abutting fragments: 1.0 over [10,110) and 1.0 over [110,210) have the *same*
        // value, so `pileup_PV` merges them into one run instead of emitting two.
        // Each expectation carries the leading zero-depth run (F230).
        let d = pileup_from_lr_hmmratac(&[(10, 110), (110, 210)], &wm);
        assert_eq!(d, vec![(10u64, 0.0f32), (210u64, 1.0f32)]);

        // Overlapping fragments give depth 2.0 in the middle, so three runs.
        let d = pileup_from_lr_hmmratac(&[(10, 110), (50, 150)], &wm);
        assert_eq!(
            d,
            vec![
                (10u64, 0.0f32),
                (50u64, 1.0f32),
                (110u64, 2.0f32),
                (150u64, 1.0f32)
            ]
        );
    }

    /// A gap between clusters is emitted as an explicit **zero-depth** run, which is why
    /// upstream's digested bedGraph contains lines like `chr1 510 1680 0.00000` spanning
    /// hundreds of bins. Dropping them would splice unrelated runs together.
    #[test]
    fn the_hmm_pileup_keeps_zero_depth_gaps() {
        let wm: BTreeMap<Coord, f32> = [(100u64, 1.0f32)].into_iter().collect();
        let d = pileup_from_lr_hmmratac(&[(10, 110), (200, 300)], &wm);
        assert_eq!(
            d,
            vec![
                (10u64, 0.0f32),
                (110u64, 1.0f32),
                (200u64, 0.0f32),
                (300u64, 1.0f32)
            ]
        );
    }

    /// A missing weight yields two endpoints that cancel, so the *only* run emitted is
    /// the leading zero-depth one that `pre_z = -10000` forces out (F230). Upstream is
    /// the same: an all-zero fragment contributes `+0` then `-0`, every comparison
    /// `z == pre_z` succeeds and the single record's end is extended.
    #[test]
    fn the_hmm_pileup_of_an_unweighted_fragment_is_one_zero_run() {
        let wm: BTreeMap<Coord, f32> = BTreeMap::new();
        let d = pileup_from_lr_hmmratac(&[(10, 110)], &wm);
        assert_eq!(d, vec![(110u64, 0.0f32)]);
    }

    #[test]
    fn posteriors_are_normalised() {
        let hmm = GaussianHMM {
            startprob: vec![0.5, 0.5],
            transmat: vec![vec![0.8, 0.2], vec![0.3, 0.7]],
            means: vec![vec![1.0, 2.0], vec![8.0, 9.0]],
            covars: Covars::Diag(vec![vec![1.0, 1.0], vec![1.0, 1.0]]),
            covariance_type: "diag".into(),
        };
        let obs = vec![
            vec![1.0, 2.0],
            vec![1.5, 2.5],
            vec![8.0, 9.0],
            vec![8.5, 9.5],
        ];
        let p = hmm.posterior(&obs);
        assert_eq!(p.len(), 4);
        for row in &p {
            let s: f64 = row.iter().sum();
            assert!((s - 1.0).abs() < 1e-9, "row sums to {s}");
        }
        assert!(p[0][0] > p[0][1]);
        assert!(p[3][1] > p[3][0]);
    }

    #[test]
    fn ln_gamma_matches_known_values() {
        assert!((ln_gamma(1.0) - 0.0).abs() < 1e-12);
        assert!((ln_gamma(5.0) - 24.0f64.ln()).abs() < 1e-10);
        assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-10);
    }

    #[test]
    fn bins_are_grouped_by_their_mark() {
        let bins = vec![
            ExtractedBin {
                chrom: ChromId(0),
                pos: 10,
                values: [1.0; 4],
                mark: 1,
            },
            ExtractedBin {
                chrom: ChromId(0),
                pos: 20,
                values: [2.0; 4],
                mark: 1,
            },
            ExtractedBin {
                chrom: ChromId(0),
                pos: 30,
                values: [3.0; 4],
                mark: 2,
            },
            ExtractedBin {
                chrom: ChromId(0),
                pos: 40,
                values: [4.0; 4],
                mark: 2,
            },
            ExtractedBin {
                chrom: ChromId(0),
                pos: 50,
                values: [5.0; 4],
                mark: 3,
            },
        ];
        let td = extract_signals_from_regions(&bins, false);
        assert_eq!(td.lengths, vec![2, 2, 1]);
        assert_eq!(td.rows.len(), 5);
        assert_eq!(td.bins, vec![10, 20, 30, 40, 50]);
    }

    #[test]
    fn poisson_extraction_truncates_to_integers() {
        let bins = vec![ExtractedBin {
            chrom: ChromId(0),
            pos: 0,
            values: [3.7, 0.5, 9.9, 0.00001],
            mark: 1,
        }];
        let td = extract_signals_from_regions(&bins, true);
        assert_eq!(td.rows[0], [3.0, 0.0, 9.0, 0.0]);
    }

    #[test]
    fn state_path_coalesces_runs_and_fills_gaps() {
        let m = model();
        let rows = vec![
            ProbRow {
                end: 10,
                probs: [0.8, 0.1, 0.1],
            },
            ProbRow {
                end: 20,
                probs: [0.8, 0.1, 0.1],
            },
            ProbRow {
                end: 40,
                probs: [0.1, 0.8, 0.1],
            },
            ProbRow {
                end: 50,
                probs: [0.8, 0.1, 0.1],
            },
            ProbRow {
                end: 60,
                probs: [0.1, 0.8, 0.1],
            },
        ];
        let p = generate_states_path(&rows, 10, &m);
        assert_eq!(
            p,
            vec![
                (0, 20, "open"),
                (20, 30, "bg"),
                (30, 40, "nuc"),
                (40, 50, "open"),
                (50, 60, "nuc"),
            ]
        );
    }

    #[test]
    fn a_tie_resolves_to_the_first_state() {
        let m = model();
        let rows = vec![ProbRow {
            end: 10,
            probs: [0.5, 0.5, 0.0],
        }];
        let p = generate_states_path(&rows, 10, &m);
        assert_eq!(p[0].2, "open");
    }

    #[test]
    fn accessible_regions_need_a_long_contiguous_nuc_open_nuc() {
        let path = vec![
            (0, 100, "nuc"),
            (100, 160, "open"),
            (160, 300, "nuc"),
            (300, 500, "bg"),
        ];
        assert_eq!(
            accessible_regions(&path, 50),
            vec![AccessibleRegion {
                start: 100,
                end: 160
            }]
        );
        assert!(accessible_regions(&path, 400).is_empty());
    }

    /// F202: a singular `full` covariance must be rejected, not panic.
    ///
    /// This is reachable from upstream's own `--modelonly` output: when a signal
    /// feature is constant across the training bins, Baum-Welch returns a singular
    /// covariance, and hmmlearn then refuses to load the model upstream just wrote.
    /// The port used to `panic!` in Gaussian elimination at that point.
    #[test]
    fn singular_full_covariance_is_rejected_not_panicked() {
        // Component 0 has two identical rows -> singular.
        let json = r#"{"startprob": [1.0, 0.0, 0.0],
          "transmat": [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],
          "means": [[1.0,1.0,1.0,1.0],[0.0,0.0,0.0,0.0],[2.0,2.0,2.0,2.0]],
          "covars": [[[1.0,1.0,0.0,0.0],[1.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],
                     [[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],
                     [[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]]],
          "covariance_type": "full", "n_features": 4,
          "i_open_region": 0, "i_background_region": 1, "i_nucleosomal_region": 2,
          "hmm_binsize": 10, "hmm_type": "gaussian"}"#;
        let err = ModelFile::from_json(json).expect_err("singular covariance must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("symmetric, positive-definite"),
            "expected hmmlearn's wording, got: {msg}"
        );
    }

    /// F202/F203: a well-conditioned `full` covariance loads, and emission is finite.
    #[test]
    fn valid_full_covariance_loads_and_emits_finite_logprob() {
        let json = std::fs::read_to_string("tests/data/model_gaussian_full.json").unwrap();
        let m = ModelFile::from_json(&json).expect("valid model must load");
        let hmm = m.gaussian();
        let x = [1.0, 2.0, 3.0, 4.0];
        for k in 0..3 {
            let v = hmm.log_emit(k, &x);
            assert!(v.is_finite(), "component {k} produced {v}");
        }
    }
    /// F203: the `full`-covariance emission must be hmmlearn's Cholesky form, not a
    /// Gaussian solve.
    ///
    /// The reference values are `hmmlearn.stats.log_multivariate_normal_density`
    /// evaluated on the model exported by upstream `hmmratac --modelonly` on the ATAC
    /// fixture, one entry per state. 13 of the 15 probed values agree to the last bit;
    /// the two that differ are 2 ULP out (~2e-19 relative), the residual difference
    /// between our Cholesky and LAPACK's.
    ///
    /// Pinning the tolerance here matters because the two formulations agree
    /// mathematically: a regression to Gaussian elimination would keep every test green
    /// except the decoded state path, where it moved 13 of 946 accessible regions.
    #[test]
    fn full_covariance_emission_matches_hmmlearn() {
        let json = r#"{"startprob": [7.084900322222042e-48, 1.0, 3.1125157273701363e-80], "transmat": [[0.947321871049102, 0.016937112279296444, 0.035741016671601546], [0.009369101376977367, 0.9870524585785625, 0.0035784400444600377], [0.07196534640393373, 0.003938565503317952, 0.9240960880927483]], "means": [[0.0020259704866438877, 1.900449239844704, 1.989623620238994, 3.3234459972855115], [0.00010000000000001094, 0.00230352645411008, 0.10716761034193238, 0.7479552486040381], [0.6447075759649408, 3.08800604346118, 2.305022929458823, 4.089378910096882]], "covars": [[[4.559040388696476e-05, 0.0017407091255318263, -0.0006229407798142503, -0.0006668627348821238], [0.0017407091255318263, 1.8958753865069093, 0.11559239350067908, 0.16632037048323814], [-0.0006229407798142503, 0.11559239350067908, 1.908102309025619, 0.9186295830552668], [-0.0006668627348821238, 0.16632037048323814, 0.9186295830552668, 4.944185159567719]], [[2.5575259911326027e-07, 2.5575259912179243e-07, 2.557525991118492e-07, 2.5575259911205367e-07], [2.5575259912179243e-07, 0.0005716115131792969, 0.0016126630298766112, 0.002953366424460839], [2.557525991118492e-07, 0.0016126630298766112, 0.08513652318486521, 0.19165050638487535], [2.5575259911205367e-07, 0.002953366424460839, 0.19165050638487535, 1.7969738531021944]], [[0.2572863489757621, -0.12145763366962549, -0.12543948239559982, -0.24617930177088929], [-0.12145763366962549, 3.4518777208745375, 0.7644033593567313, 0.9097969446861179], [-0.12543948239559982, 0.7644033593567313, 2.821717860645719, 1.8754600851820842], [-0.24617930177088929, 0.9097969446861179, 1.8754600851820842, 8.182179701844577]]], "covariance_type": "full", "n_features": 4, "i_open_region": 2, "i_background_region": 1, "i_nucleosomal_region": 0, "hmm_binsize": 10, "hmm_type": "gaussian"}"#;
        let m = ModelFile::from_json(json).expect("valid model must load");
        let hmm = m.gaussian();

        // (x, [state0, state1, state2]) straight out of hmmlearn. The literals carry
        // more digits than f64 can distinguish because that is what the reference
        // printed; truncating them would make the comparison assert the wrong thing.
        #[allow(clippy::excessive_precision)]
        let cases: Vec<(Vec<f64>, [f64; 3])> = vec![
            (
                vec![1.0, 1.9, 1.99, 3.32],
                [
                    -1.14049976897639535e4,
                    -1.95533706320069241e6,
                    -5.45007988798625576,
                ],
            ),
            (
                vec![0.5, 2.0, 2.0, 4.0],
                [
                    -2.83954432910250853e3,
                    -4.90563182041534805e5,
                    -5.28876681447424968,
                ],
            ),
            (
                vec![0.1, 1.0, 1.0, 1.0],
                [
                    -1.11400766428418578e2,
                    -2.02265876657252848e4,
                    -7.08272034722023136,
                ],
            ),
            (
                vec![0.0, 0.5, 0.5, 0.5],
                [
                    -1.56829661987288027,
                    -2.14573657687860248e2,
                    -8.02739383597773148,
                ],
            ),
        ];
        for (x, expect) in cases {
            for (k, &want) in expect.iter().enumerate() {
                let got = hmm.log_emit(k, &x);
                let rel = (got - want).abs() / want.abs().max(1e-300);
                assert!(
                    rel <= 1e-14,
                    "state {k} at {x:?}: got {got:.17e}, hmmlearn {want:.17e} (rel {rel:.3e})"
                );
            }
        }
    }
}
