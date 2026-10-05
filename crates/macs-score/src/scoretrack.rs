//! `ScoreTrackII` and `TwoConditionScores`, a port of `MACS3/Signal/ScoreTrack.py`.
//!
//! These back `bdgcmp` and `bdgdiff`. Unlike the peak-caller's p-score track,
//! they work on **bedGraph interval endpoints**, not per-base pileup: the
//! treatment and control tracks are merged into a common set of transition
//! points, each carrying `(pos, treat, ctrl, score)` rows.
//!
//! # float32 is load-bearing here
//!
//! Every numeric column in upstream is a `numpy` array with an explicit
//! `dtype`: `int32` for positions, `float32` for pileups and scores. The
//! likelihood ratios are therefore computed from `f32` inputs in `f32` and
//! stored as `f32`. This is *not* an implementation detail -- the asymmetry in
//! `logLR_asym` and the `f32` rounding both feed directly into the `%.5f`
//! bedGraph output, so promoting to `f64` would visibly change results.
//! Each arithmetic operation below is annotated with the upstream line it
//! mirrors.

use std::collections::BTreeMap;
use std::path::Path;

use macs_core::{ChromId, Coord, Genome, MacsError};
use macs_rle::SignalTrack;
use macs_stats::poisson_cdf;

/// `LOG10_E` (`ScoreTrack.py:62`), `1/ln(10)`.
///
/// The annotation reads `LOG10_E: cython.float = 0.43429448190325176`, which
/// looks like an `f32` constant, but the binding is at **module scope**: Cython
/// stores a module-level typed variable as a Python attribute, so the value is a
/// full `double` and the annotation is not enforced there. Verified against
/// upstream over all 76 rows of both `-m logLR` and `-m slogLR`: using the
/// `f32` narrowing of this constant mismatches 21/76 and 20/76 rows, while the
/// `f64` value matches every one.
const LOG10_E: f64 = 0.434_294_481_903_251_76;

/// Which score a [`ScoreTrack`] column holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoreMethod {
    /// `-log10(p)`.
    P,
    /// `-log10(q)`; implies [`ScoreMethod::P`] first.
    Q,
    /// asymmetric `log10` likelihood ratio.
    LogLR,
    /// symmetric `log10` likelihood ratio.
    SymLogLR,
    /// `log10` fold enrichment.
    LogFE,
    /// linear fold enrichment.
    FE,
    /// treatment minus control.
    Subtract,
    /// treatment pileup per million reads.
    SPMR,
    /// element-wise `max(treat, ctrl)`.
    Max,
    /// not yet computed.
    None,
}

/// Normalisation of the treatment/control columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormMethod {
    /// scale to the treatment depth.
    Treat,
    /// scale to the control depth.
    Ctrl,
    /// scale to one million reads.
    Million,
    /// raw pileup.
    Raw,
}

/// One chromosome's four-column score table: `(end position, treat, ctrl, score)`.
///
/// `pos` is a right endpoint: the row describes `[prev_pos, pos)`.
#[derive(Debug, Clone, Default)]
pub struct ScoreTrack {
    /// End position of each interval.
    pub pos: Vec<Coord>,
    /// Treatment pileup.
    pub treat: Vec<f32>,
    /// Control pileup.
    pub ctrl: Vec<f32>,
    /// The score column, in whatever unit `scoring_method` says.
    pub score: Vec<f32>,
}

/// `ScoreTrackII` (`ScoreTrack.py:166-760`): treatment/control interval scores.
///
/// Ported for `bdgcmp`. `TwoConditionScores` ([`TwoScores`]) is the `bdgdiff`
/// analogue.
#[derive(Debug, Clone)]
pub struct ScoreTrack2 {
    /// Chromosome -> score table. A `BTreeMap` keyed by id; iteration that
    /// affects output order goes through `chroms_sorted`.
    data: BTreeMap<ChromId, ScoreTrack>,
    /// The genome shared with the inputs, for name resolution and ordering.
    pub(crate) genome: Genome,
    treat_edm: f32,
    ctrl_edm: f32,
    scoring_method: ScoreMethod,
    normalization_method: NormMethod,
    pseudocount: f32,
    /// `-log10 p` -> number of scored base pairs, for the p->q map.
    pvalue_stat: BTreeMap<u32, u64>,
    /// Whether to emit a UCSC `track` line.
    trackline: bool,
}

impl ScoreTrack2 {
    /// A score track for treatments at `treat_depth` million reads and controls
    /// at `ctrl_depth` million reads.
    pub fn new(genome: Genome, treat_depth: f32, ctrl_depth: f32) -> Self {
        Self {
            data: BTreeMap::new(),
            genome,
            treat_edm: treat_depth,
            ctrl_edm: ctrl_depth,
            scoring_method: ScoreMethod::None,
            normalization_method: NormMethod::Raw,
            pseudocount: 1.0,
            pvalue_stat: BTreeMap::new(),
            trackline: false,
        }
    }

    /// The genome the input chromosomes were interned into.
    pub fn genome(&self) -> &Genome {
        &self.genome
    }

    /// `set_pseudocount`.
    pub fn set_pseudocount(&mut self, pseudocount: f32) {
        self.pseudocount = pseudocount;
    }

    /// `enable_trackline`.
    pub fn enable_trackline(&mut self) {
        self.trackline = true;
    }

    /// Append one interval row (`add`, `ScoreTrack.py:263-282`).
    pub fn add(&mut self, chrom: ChromId, endpos: Coord, treat: f32, ctrl: f32) {
        let e = self.data.entry(chrom).or_default();
        e.pos.push(endpos);
        e.treat.push(treat);
        e.ctrl.push(ctrl);
        // score is allocated lazily; push a placeholder so it stays parallel
        e.score.push(0.0);
    }

    /// Chromosome ids in byte-wise name order.
    pub fn chroms_sorted(&self) -> Vec<ChromId> {
        let mut ids: Vec<ChromId> = self.data.keys().copied().collect();
        ids.sort_by(|&a, &b| self.genome.name(a).cmp(self.genome.name(b)));
        ids
    }

    /// The score table for a chromosome.
    pub fn get(&self, chrom: ChromId) -> Option<&ScoreTrack> {
        self.data.get(&chrom)
    }

    /// `change_normalization_method` (`ScoreTrack.py:317-374`).
    ///
    /// The upstream comment warns that flipping between methods loses precision
    /// ("I only keep two digits"). That warning is real and the transition
    /// table below is reproduced verbatim, including the asymmetry where
    /// switching *to* `Treat` from `Raw` normalises only the control side.
    pub fn change_normalization_method(&mut self, method: NormMethod) {
        let (treat_edm, ctrl_edm) = (self.treat_edm, self.ctrl_edm);
        match method {
            NormMethod::Treat => match self.normalization_method {
                NormMethod::Treat => {}
                NormMethod::Ctrl => {
                    let s = treat_edm / ctrl_edm;
                    self.normalize(s, s)
                }
                NormMethod::Million => self.normalize(treat_edm, treat_edm),
                NormMethod::Raw => self.normalize(1.0, treat_edm / ctrl_edm),
            },
            NormMethod::Ctrl => match self.normalization_method {
                NormMethod::Treat => {
                    let s = ctrl_edm / treat_edm;
                    self.normalize(s, s)
                }
                NormMethod::Ctrl => {}
                NormMethod::Million => self.normalize(ctrl_edm, ctrl_edm),
                NormMethod::Raw => self.normalize(ctrl_edm / treat_edm, 1.0),
            },
            NormMethod::Million => match self.normalization_method {
                NormMethod::Treat => self.normalize(1.0 / treat_edm, 1.0 / treat_edm),
                NormMethod::Ctrl => self.normalize(1.0 / ctrl_edm, 1.0 / ctrl_edm),
                NormMethod::Million => {}
                NormMethod::Raw => self.normalize(1.0 / treat_edm, 1.0 / ctrl_edm),
            },
            NormMethod::Raw => match self.normalization_method {
                NormMethod::Treat => self.normalize(treat_edm, treat_edm),
                NormMethod::Ctrl => self.normalize(ctrl_edm, ctrl_edm),
                NormMethod::Million => self.normalize(treat_edm, ctrl_edm),
                NormMethod::Raw => {}
            },
        }
        self.normalization_method = method;
    }

    /// `normalize` (`ScoreTrack.py:376-390`): scale both columns in place.
    fn normalize(&mut self, treat_scale: f32, control_scale: f32) {
        for e in self.data.values_mut() {
            for v in e.treat.iter_mut() {
                *v *= treat_scale;
            }
            for v in e.ctrl.iter_mut() {
                *v *= control_scale;
            }
        }
    }

    /// `change_score_method` (`ScoreTrack.py:391-423`).
    pub fn change_score_method(&mut self, method: ScoreMethod) {
        match method {
            ScoreMethod::P => self.compute_pvalue(),
            ScoreMethod::Q => {
                // p must exist to build the pq table
                if self.scoring_method != ScoreMethod::P {
                    self.compute_pvalue();
                }
                self.compute_qvalue();
            }
            ScoreMethod::LogLR => self.compute_likelihood(),
            ScoreMethod::SymLogLR => self.compute_sym_likelihood(),
            ScoreMethod::LogFE => self.compute_logfe(),
            ScoreMethod::FE => self.compute_foldenrichment(),
            ScoreMethod::Subtract => self.compute_subtraction(),
            ScoreMethod::SPMR => self.compute_spmr(),
            ScoreMethod::Max => self.compute_max(),
            ScoreMethod::None => {}
        }
    }

    /// `compute_pvalue` (`ScoreTrack.py:425-452`).
    ///
    /// Note the treatment side is **truncated to an integer** after the `f32`
    /// pseudocount add (`cython.cast(cython.int, ...)`), while the lambda side
    /// stays `f32` -- the same asymmetry as the peak caller's p-score.
    fn compute_pvalue(&mut self) {
        let pc = self.pseudocount;
        for chrom in self.chroms_sorted() {
            let e = &self.data[&chrom];
            let mut prev_pos: Coord = 0;
            let mut score = Vec::with_capacity(e.pos.len());
            for i in 0..e.pos.len() {
                let v = pscore_interval(e.treat[i], e.ctrl[i], pc);
                score.push(v);
                let len = e.pos[i] - prev_pos;
                *self.pvalue_stat.entry(macs_score_key(v)).or_insert(0) += u64::from(len);
                prev_pos = e.pos[i];
            }
            self.data.get_mut(&chrom).expect("present").score = score;
        }
        self.scoring_method = ScoreMethod::P;
    }

    /// `compute_qvalue` (`ScoreTrack.py:454-476`): map `-log10 p` through the
    /// AFDR table built by [`Self::pq_table`].
    fn compute_qvalue(&mut self) {
        let table = self.pq_table();
        for chrom in self.chroms_sorted() {
            let src = self.data[&chrom].score.clone();
            let mapped: Vec<f32> = src
                .iter()
                .map(|&v| table.get(&macs_score_key(v)).copied().unwrap_or(0.0))
                .collect();
            self.data.get_mut(&chrom).expect("present").score = mapped;
        }
        self.scoring_method = ScoreMethod::Q;
    }

    /// `make_pq_table` (`ScoreTrack.py:478-521`).
    ///
    /// Mirrors `BedGraph.p2q`: the rank `k` starts at 1 and is incremented by
    /// the block length, so `log10(k)` is **not** dropped here (unlike the
    /// peak-caller variant) -- it genuinely ranks by cumulative bases.
    fn pq_table(&self) -> BTreeMap<u32, f32> {
        let total: u64 = self.pvalue_stat.values().sum();
        let mut out = BTreeMap::new();
        if total == 0 {
            return out;
        }
        // `f = -log10(N)` where N is an integer count. The `log10` runs in f64
        // (C promotes the `long` operand) and the result is stored into a
        // `cython.float`, so `f` is an f32 *rounded from* the f64 logarithm --
        // not `log10` of an f32 N. Verified: 43 of 165 `qpois` rows mismatch
        // when the sum is evaluated in f32, 0 when evaluated in f64.
        let f = (-((total as f64).log10())) as f32;
        let mut k = 1.0f32;
        // upstream initialises `pre_q = 2147483647`, an `int` literal, not
        // FLT_MAX -- both work as a clamp ceiling, but the literal documents
        // that this is the same sentinel `BedGraph.p2q` uses.
        let mut pre_q = 2_147_483_647.0f32;
        // distinct scores, descending
        let mut uniq: Vec<f32> = self
            .pvalue_stat
            .keys()
            .map(|&b| f32::from_bits(b))
            .collect();
        uniq.sort_by(|a, b| b.partial_cmp(a).expect("not NaN"));
        // The `q` expression is evaluated in f64 and narrowed once, even though
        // `v`, `k` and `f` are all f32: C promotes the operands and the local
        // `q` is declared `cython.float`, so the rounding happens at the store.
        //
        // Upstream's post-loop back-fill starts at `i`, the loop counter, and
        // `i` is left at `len - 1` when the loop finishes without breaking. So
        // the **last** score is unconditionally overwritten with 0, whether or
        // not the `q <= 0` break fired. Verified against upstream on the
        // 76-row `qpois` fixture: p = 0.00639 keeps q = 0.00156 while the
        // immediately following p = 0.00537 is forced to 0 even though its
        // computed q is a healthy 0.0053 > 0.
        let mut backfill_from = uniq.len();
        for (i, &v) in uniq.iter().enumerate() {
            let ln = self.pvalue_stat[&macs_score_key(v)] as f32;
            let q = ((f64::from(v) + (f64::from(k).log10() + f64::from(f))) as f32).min(pre_q);
            backfill_from = i;
            if q <= 0.0 {
                break;
            }
            out.insert(v.to_bits(), q);
            pre_q = q;
            k += ln;
        }
        // bottom-ranked pscores all get qscore 0
        for &v in &uniq[backfill_from.min(uniq.len())..] {
            out.insert(v.to_bits(), 0.0);
        }
        out
    }

    /// `compute_likelihood` (`ScoreTrack.py:523-549`).
    fn compute_likelihood(&mut self) {
        let pc = self.pseudocount;
        for chrom in self.chroms_sorted() {
            let e = self.data.get_mut(&chrom).expect("present");
            for i in 0..e.pos.len() {
                e.score[i] = log_lr_asym(e.treat[i] + pc, e.ctrl[i] + pc);
            }
        }
        self.scoring_method = ScoreMethod::LogLR;
    }

    /// `compute_sym_likelihood` (`ScoreTrack.py:551-573`).
    fn compute_sym_likelihood(&mut self) {
        let pc = self.pseudocount;
        for chrom in self.chroms_sorted() {
            let e = self.data.get_mut(&chrom).expect("present");
            for i in 0..e.pos.len() {
                e.score[i] = log_lr_sym(e.treat[i] + pc, e.ctrl[i] + pc);
            }
        }
        self.scoring_method = ScoreMethod::SymLogLR;
    }

    /// `compute_logFE` (`ScoreTrack.py:575-597`).
    /// `compute_logFE` (`ScoreTrack.py:575-597`).
    ///
    /// The division happens in `f32` -- it is the argument to `get_logFE`, whose
    /// parameter is `cython.float` -- but the `log10` runs in `f64`, since C
    /// promotes the `float` argument to `double` on the way into `log10`. The
    /// `f32` rounding of the ratio is therefore observable in the result: with
    /// `treat=15.2`, `ctrl=0.0`, `pseudocount=1.0`, upstream scores
    /// `1.2095149755477905`, which is `f32(log10(f64(f32(16.2/1.0))))`. Doing
    /// the division in `f64` as well gives `1.209515095`, a different `f32`
    /// printing as `1.20952` instead of `1.20951`.
    fn compute_logfe(&mut self) {
        let pc = self.pseudocount;
        for chrom in self.chroms_sorted() {
            let e = self.data.get_mut(&chrom).expect("present");
            for i in 0..e.pos.len() {
                let ratio = (e.treat[i] + pc) / (e.ctrl[i] + pc);
                e.score[i] = f64::from(ratio).log10() as f32;
            }
        }
        self.scoring_method = ScoreMethod::LogFE;
    }

    /// `compute_foldenrichment` (`ScoreTrack.py:599-621`).
    ///
    /// Unlike every other score here, the pseudocount is added in `f64` and the
    /// division runs in `f64`, with a single narrowing to `f32` at the end.
    /// Upstream writes `(p[i] + pseudocount)/(c[i] + pseudocount)` with numpy
    /// `float32` array elements and a `cython.float` pseudocount, and that
    /// mixed-type expression is evaluated in double before being stored back
    /// into the `float32` score array. Verified against upstream: with
    /// `treat=13.4`, `ctrl=6.3`, `pseudocount=0.1`, upstream scores
    /// `2.109374761581421`, which is `f32((f64(13.4) + f64(0.1)) / (f64(6.3) + f64(0.1)))`;
    /// adding in `f32` first gives `13.5/6.4 = 2.109375`, a different bit
    /// pattern that prints as `2.10938` instead of `2.10937`.
    fn compute_foldenrichment(&mut self) {
        let pc = f64::from(self.pseudocount);
        for chrom in self.chroms_sorted() {
            let e = self.data.get_mut(&chrom).expect("present");
            for i in 0..e.pos.len() {
                e.score[i] = ((f64::from(e.treat[i]) + pc) / (f64::from(e.ctrl[i]) + pc)) as f32;
            }
        }
        self.scoring_method = ScoreMethod::FE;
    }

    /// `compute_subtraction` (`ScoreTrack.py:623-643`).
    fn compute_subtraction(&mut self) {
        for chrom in self.chroms_sorted() {
            let e = self.data.get_mut(&chrom).expect("present");
            for i in 0..e.pos.len() {
                e.score[i] = e.treat[i] - e.ctrl[i];
            }
        }
        self.scoring_method = ScoreMethod::Subtract;
    }

    /// `compute_SPMR` (`ScoreTrack.py:645-670`).
    fn compute_spmr(&mut self) {
        // `score` is in the f32 array, so the division rounds there
        let scale = match self.normalization_method {
            NormMethod::Treat | NormMethod::Raw => self.treat_edm,
            NormMethod::Ctrl => self.ctrl_edm,
            NormMethod::Million => 1.0,
        };
        for chrom in self.chroms_sorted() {
            let e = self.data.get_mut(&chrom).expect("present");
            for i in 0..e.pos.len() {
                e.score[i] = e.treat[i] / scale;
            }
        }
        self.scoring_method = ScoreMethod::SPMR;
    }

    /// `compute_max` (`ScoreTrack.py:672-694`).
    fn compute_max(&mut self) {
        for chrom in self.chroms_sorted() {
            let e = self.data.get_mut(&chrom).expect("present");
            for i in 0..e.pos.len() {
                e.score[i] = e.treat[i].max(e.ctrl[i]);
            }
        }
        self.scoring_method = ScoreMethod::Max;
    }

    /// `write_bedGraph` (`ScoreTrack.py:697-770`).
    ///
    /// Two upstream quirks matter here:
    /// - the change-detection threshold is `abs(pre_v - v) > 1e-5` for
    ///   `ScoreTrackII` but `>= 1e-6` for `TwoConditionScores`;
    /// - a run is emitted as `[pre, pos[i-1])` using the **previous** row's
    ///   position, then the final row is written unconditionally.
    pub fn write_bedgraph(
        &self,
        path: &Path,
        name: &str,
        description: &str,
        column: usize,
    ) -> macs_core::Result<()> {
        assert!((1..=3).contains(&column), "column should be 1, 2 or 3");
        let mut out = String::new();
        if self.trackline {
            out.push_str(&format!(
                "track type=bedGraph name=\"{}\" description=\"{}\"\n",
                name, description
            ));
        }
        for chrom in self.chroms_sorted() {
            let e = &self.data[&chrom];
            let values: &Vec<f32> = match column {
                1 => &e.treat,
                2 => &e.ctrl,
                _ => &e.score,
            };
            if e.pos.is_empty() {
                continue;
            }
            let mut pre: Coord = 0;
            let mut pre_v = values[0];
            for (i, &v) in values.iter().enumerate().skip(1) {
                let p = e.pos[i - 1];
                if (pre_v - v).abs() > 1e-5 {
                    out.push_str(&row(&self.genome, chrom, pre, p, pre_v));
                    pre_v = v;
                    pre = p;
                }
            }
            out.push_str(&row(
                &self.genome,
                chrom,
                pre,
                *e.pos.last().expect("non-empty"),
                pre_v,
            ));
        }
        std::fs::write(path, out).map_err(MacsError::Io)
    }
}

/// One bedGraph row in upstream's format.
fn row(genome: &Genome, chrom: ChromId, start: Coord, end: Coord, value: f32) -> String {
    format!(
        "{}\t{}\t{}\t{:.5}\n",
        String::from_utf8_lossy(genome.name(chrom)),
        start,
        end,
        value
    )
}

/// `-log10 P(X > observed | expectation)` for one interval.
///
/// Takes the raw treatment pileup and control pileup plus the pseudocount, and
/// reproduces `get_pscore(int(p[i] + pseudocount), c[i] + pseudocount)`
/// (`ScoreTrack.py:451-453`). The treatment add is widened to `f64` before the
/// integer truncation -- see [`pseudocounted_inputs`](crate::pseudocounted_inputs) for why the
/// `f32` rounding of that sum gives the wrong count.
///
/// # Panics
/// If the control pileup plus pseudocount is not positive. With the default
/// `bdgcmp` pseudocount of `0.0` and a control track containing explicit zeros
/// upstream raises `AssertionError: Lambda must > 0` in exactly the same place,
/// so this is a faithful failure rather than a divergence -- but it does mean
/// `bdgcmp -m ppois -p 0` crashes on sparse control bedGraphs in both.
fn pscore_interval(treat: f32, ctrl: f32, pseudocount: f32) -> f32 {
    let (observed, expectation) = crate::pseudocounted_inputs(treat, ctrl, pseudocount);
    let v = poisson_cdf(observed, f64::from(expectation), false, true)
        .expect("pseudocounted control is positive");
    (-v) as f32
}

/// The two likelihood-ratio formulas, evaluated the way C does.
///
/// # Why the arithmetic is `f64` throughout, not `f32`
///
/// Upstream's signature is `logLR_asym(x: cython.float, y: cython.float)`, and
/// in a C expression every `float` operand is promoted to `double`. So the
/// whole body -- both logarithms, the products, the sum and the `log10(e)`
/// scale -- runs in `f64`, and the result is narrowed to `f32` **once**, on
/// assignment to the local `s` which is also declared `cython.float`.
///
/// Computing in `f32` is *not* equivalent, and fails twice over:
///
/// 1. `ln` is not correctly rounded, so `f32::ln` and `f64::ln` narrowed to
///    `f32` disagree on about 4% of inputs -- 15 of 49999 sampled values.
/// 2. Each `f32` rounding of an intermediate lands in a different digit than a
///    single rounding at the end.
///
/// Verified against upstream on `treat=12.3`, `ctrl=1.9`, `pseudocount=0.1`,
/// which upstream scores as `5.308994770050049`:
///
/// ```text
/// f32 throughout          -> 5.308995723724365   (wrong, prints 5.30900)
/// f32 log, f64 arithmetic -> 5.308994293212891   (wrong, prints 5.30899 -- last digit)
/// f64 throughout          -> 5.308994770050049   (exact)
/// ```
///
/// The inputs are still `f32` -- they come from `float32` arrays, and the
/// pseudocount is added in `f32` to match `(p[i] + self.pseudocount)` -- so only
/// the *evaluation* is widened.
fn log_lr(x: f32, y: f32) -> [f64; 2] {
    let xd = f64::from(x);
    let yd = f64::from(y);
    let lx = xd.ln();
    let ly = yd.ln();
    [
        (xd * (lx - ly) + yd - xd) * LOG10_E,
        (xd * (-lx + ly) - yd + xd) * LOG10_E,
    ]
}

/// `logLR_asym` (`ScoreTrack.py:110-127`).
///
/// The two branches are algebraically different expressions that only coincide
/// at `x == y`, so both are kept rather than simplified:
/// ```text
/// x > y:  (x*(ln x - ln y) + y - x) * log10(e)
/// x < y:  (x*(-ln x + ln y) - y + x) * log10(e)
/// ```
/// Upstream memoises on the `f32` bit patterns; the arithmetic is identical
/// either way (a cache hit returns the same value), so no cache is kept.
pub fn log_lr_asym(x: f32, y: f32) -> f32 {
    let [gt, lt] = log_lr(x, y);
    if x > y {
        gt as f32
    } else if x < y {
        lt as f32
    } else {
        0.0
    }
}

/// `logLR_sym` (`ScoreTrack.py:138-153`).
///
/// Note the `y > x` branch reuses `log(x) - log(y)`, i.e. the *negative* of
/// what the name suggests -- it is not a symmetric rewrite of the `x > y`
/// branch. Reproduced verbatim; "fixing" it would change `bdgcmp -m slogLR`
/// output.
pub fn log_lr_sym(x: f32, y: f32) -> f32 {
    let [gt, _] = log_lr(x, y);
    if x > y {
        gt as f32
    } else if y > x {
        // upstream: `(y*(log(x) - log(y)) + y - x) * log10(e)`
        let xd = f64::from(x);
        let yd = f64::from(y);
        ((yd * (xd.ln() - yd.ln()) + yd - xd) * LOG10_E) as f32
    } else {
        0.0
    }
}

/// Cache key for a score value.
///
/// Two `-0.0` and `0.0` are the same score for lookup purposes, and `-0.0` must
/// sort below `0.0` to keep the descending walk monotone, so the sign of a zero
/// is cleared here rather than relying on float comparison.
fn macs_score_key(v: f32) -> u32 {
    if v == 0.0 {
        0
    } else {
        v.to_bits()
    }
}

/// One interval of a [`TwoScores`] table: end position and three LLR columns.
#[derive(Debug, Clone, Default)]
pub struct TwoScoreRow {
    /// End position of the interval.
    pub end: Coord,
    /// `logLR` of condition 1 treatment vs condition 1 control.
    pub t1_vs_c1: f32,
    /// `logLR` of condition 2 treatment vs condition 2 control.
    pub t2_vs_c2: f32,
    /// symmetric `logLR` of condition 1 treatment vs condition 2 treatment.
    pub t1_vs_t2: f32,
}

/// A peak reported by [`TwoScores::call_peaks`].
#[derive(Debug, Clone)]
pub struct DiffPeak {
    /// Chromosome.
    pub chrom: ChromId,
    /// Inclusive start coordinate, XLS convention.
    pub start: Coord,
    /// Inclusive end coordinate, XLS convention.
    pub end: Coord,
    /// Base-length-weighted mean LLR across the merged region.
    pub score: f32,
}

/// `TwoConditionScores` (`ScoreTrack.py:1378-1960`): differential scores for
/// `bdgdiff`.
///
/// Holds three likelihood-ratio columns per interval -- condition 1 vs its own
/// control, condition 2 vs its own control, and condition 1 vs condition 2 --
/// and calls three peak sets from them.
#[derive(Debug, Clone)]
pub struct TwoScores {
    data: BTreeMap<ChromId, Vec<TwoScoreRow>>,
    genome: Genome,
    cond1_factor: f32,
    cond2_factor: f32,
    pseudocount: f32,
}

impl TwoScores {
    /// Build from four bedGraph tracks: `(t1, c1, t2, c2)`, each mapping
    /// chromosome to a run-length track, together with the genome they share.
    ///
    /// Only chromosomes present in **all four** tracks are scored, matching
    /// `get_common_chrs`.
    /// Four bedGraphs plus the two scaling factors and pseudocount; the arity
    /// mirrors upstream's constructor.
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        genome: Genome,
        t1: &BTreeMap<ChromId, SignalTrack<f32>>,
        c1: &BTreeMap<ChromId, SignalTrack<f32>>,
        t2: &BTreeMap<ChromId, SignalTrack<f32>>,
        c2: &BTreeMap<ChromId, SignalTrack<f32>>,
        cond1_factor: f32,
        cond2_factor: f32,
        pseudocount: f32,
    ) -> Self {
        let mut me = TwoScores {
            data: BTreeMap::new(),
            genome,
            cond1_factor,
            cond2_factor,
            pseudocount,
        };
        // four-way chromosome intersection, in name order
        let mut common: Vec<ChromId> = t1
            .keys()
            .copied()
            .filter(|k| c1.contains_key(k) && t2.contains_key(k) && c2.contains_key(k))
            .collect();
        common.sort_by(|&a, &b| me.genome.name(a).cmp(me.genome.name(b)));

        for chrom in common {
            let rows =
                me.build_chromosome(chrom, &t1[&chrom], &c1[&chrom], &t2[&chrom], &c2[&chrom]);
            me.data.insert(chrom, rows);
        }
        me
    }

    /// `build_chromosome` (`ScoreTrack.py:1468-1538`).
    ///
    /// A four-way merge: each iteration records the *current* four values
    /// against the **previous** minimum, then advances every track sitting at
    /// that minimum. So the row is `[pre_p, minp)` and the loop stops as soon as
    /// any track is exhausted -- the tail of the longer tracks is dropped, which
    /// is why a single unmatched row can change the last interval.
    fn build_chromosome(
        &self,
        chrom: ChromId,
        t1: &SignalTrack<f32>,
        c1: &SignalTrack<f32>,
        t2: &SignalTrack<f32>,
        c2: &SignalTrack<f32>,
    ) -> Vec<TwoScoreRow> {
        let r1 = t1.runs();
        let r2 = c1.runs();
        let r3 = t2.runs();
        let r4 = c2.runs();
        let mut out = Vec::with_capacity(r1.len() + r2.len() + r3.len() + r4.len());
        let (mut i1, mut i2, mut i3, mut i4) = (0usize, 0usize, 0usize, 0usize);
        let mut pre_p: Coord = 0;
        while let (Some(&p1), Some(&p2), Some(&p3), Some(&p4)) = (
            r1.get(i1).map(|x| &x.end),
            r2.get(i2).map(|x| &x.end),
            r3.get(i3).map(|x| &x.end),
            r4.get(i4).map(|x| &x.end),
        ) {
            let minp = p1.min(p2).min(p3).min(p4);
            out.push(self.score_row(
                pre_p,
                r1[i1].value,
                r2[i2].value,
                r3[i3].value,
                r4[i4].value,
            ));
            pre_p = minp;
            if p1 == minp {
                i1 += 1;
            }
            if p2 == minp {
                i2 += 1;
            }
            if p3 == minp {
                i3 += 1;
            }
            if p4 == minp {
                i4 += 1;
            }
        }
        let _ = chrom;
        out
    }

    /// `add` (`ScoreTrack.py:1571-1600`).
    fn score_row(&self, end: Coord, t1: f32, c1: f32, t2: f32, c2: f32) -> TwoScoreRow {
        let pc = self.pseudocount;
        TwoScoreRow {
            end,
            t1_vs_c1: log_lr_asym((t1 + pc) * self.cond1_factor, (c1 + pc) * self.cond1_factor),
            t2_vs_c2: log_lr_asym((t2 + pc) * self.cond2_factor, (c2 + pc) * self.cond2_factor),
            t1_vs_t2: log_lr_sym((t1 + pc) * self.cond1_factor, (t2 + pc) * self.cond2_factor),
        }
    }

    /// The genome shared with the inputs.
    pub fn genome(&self) -> &Genome {
        &self.genome
    }

    /// Chromosome ids in byte-wise name order.
    pub fn chroms_sorted(&self) -> Vec<ChromId> {
        let mut ids: Vec<ChromId> = self.data.keys().copied().collect();
        ids.sort_by(|&a, &b| self.genome.name(a).cmp(self.genome.name(b)));
        ids
    }

    /// `write_bedGraph` (`ScoreTrack.py:1641-1694`).
    ///
    /// As with `ScoreTrackII`, rows are `[pre, pos[i-1])` and the last row is
    /// always written -- but the change threshold is `>= 1e-6`, not `> 1e-5`.
    pub fn write_bedgraph(&self, path: &Path, column: usize) -> macs_core::Result<()> {
        assert!((1..=3).contains(&column), "column should be 1, 2 or 3");
        let mut out = String::new();
        for chrom in self.chroms_sorted() {
            let rows = &self.data[&chrom];
            if rows.is_empty() {
                continue;
            }
            let value_at = |r: &TwoScoreRow| -> f32 {
                match column {
                    1 => r.t1_vs_c1,
                    2 => r.t2_vs_c2,
                    _ => r.t1_vs_t2,
                }
            };
            let mut pre: Coord = 0;
            let mut pre_v = value_at(&rows[0]);
            for i in 1..rows.len() {
                let v = value_at(&rows[i]);
                let p = rows[i - 1].end;
                if (pre_v - v).abs() >= 1e-6 {
                    out.push_str(&row(&self.genome, chrom, pre, p, pre_v));
                    pre_v = v;
                    pre = p;
                }
            }
            out.push_str(&row(
                &self.genome,
                chrom,
                pre,
                rows.last().expect("non-empty").end,
                pre_v,
            ));
        }
        std::fs::write(path, out).map_err(MacsError::Io)
    }

    /// `write_matrix` (`ScoreTrack.py:1719-1751`): all three columns per row.
    pub fn write_matrix(&self, path: &Path) -> macs_core::Result<()> {
        let mut out = String::new();
        for chrom in self.chroms_sorted() {
            let rows = &self.data[&chrom];
            let mut pre: Coord = 0;
            for r in rows {
                out.push_str(&format!(
                    "{}:{}_{}\t{:.5}\t{:.5}\t{:.5}\n",
                    String::from_utf8_lossy(self.genome.name(chrom)),
                    pre,
                    r.end,
                    r.t1_vs_c1,
                    r.t2_vs_c2,
                    r.t1_vs_t2
                ));
                pre = r.end;
            }
        }
        std::fs::write(path, out).map_err(MacsError::Io)
    }

    /// `call_peaks` (`ScoreTrack.py:1749-1860`): the three peak sets.
    ///
    /// Returns `(cond1-only, cond2-only, common)`. Peaks are regions where the
    /// score stays above `cutoff`, merged across gaps of at most `max_gap`, and
    /// kept only if at least `min_length` bases long.
    pub fn call_peaks(
        &self,
        cutoff: f32,
        min_length: Coord,
        max_gap: Coord,
    ) -> (Vec<DiffPeak>, Vec<DiffPeak>, Vec<DiffPeak>) {
        let mut cat1 = Vec::new();
        let mut cat2 = Vec::new();
        let mut cat3 = Vec::new();
        for chrom in self.chroms_sorted() {
            let rows = &self.data[&chrom];
            let n = rows.len();
            let mut idx1 = Vec::new();
            let mut idx2 = Vec::new();
            let mut idx3 = Vec::new();
            for (i, r) in rows.iter().enumerate() {
                let cond1_over_cond2 = r.t1_vs_t2 >= cutoff;
                let cond2_over_cond1 = r.t1_vs_t2 <= -cutoff;
                let cond1_eq_cond2 = r.t1_vs_t2 >= -cutoff && r.t1_vs_t2 <= cutoff;
                let cond1_sig = r.t1_vs_c1 >= cutoff;
                let cond2_sig = r.t2_vs_c2 >= cutoff;
                if cond1_sig && cond1_over_cond2 {
                    idx1.push(i);
                }
                if cond2_over_cond1 && cond2_sig {
                    idx2.push(i);
                }
                if cond1_sig && cond2_sig && cond1_eq_cond2 {
                    idx3.push(i);
                }
            }
            for (idx, sign, sink) in [
                (idx1, 1.0f32, &mut cat1),
                (idx2, -1.0f32, &mut cat2),
                (idx3, 0.0f32, &mut cat3),
            ] {
                // signed/abs score used for the weighted mean
                let sc: Vec<f32> = rows
                    .iter()
                    .map(|r| match sign {
                        0.0 => r.t1_vs_t2.abs(),
                        s => r.t1_vs_t2 * s,
                    })
                    .collect();
                self.add_a_peak(sink, chrom, &idx, rows, &sc, max_gap, min_length);
            }
            let _ = n;
        }
        (cat1, cat2, cat3)
    }

    /// `__add_a_peak` (`ScoreTrack.py:1862-1920`).
    ///
    /// Note the length test uses `>=` against `min_length`, that the first
    /// region is snapped to start at 0 when it is the very first row (otherwise
    /// `pos[-1]` would wrap around), and that the score is a base-length-weighted
    /// mean, not a maximum.
    #[allow(clippy::too_many_arguments)]
    fn add_a_peak(
        &self,
        sink: &mut Vec<DiffPeak>,
        chrom: ChromId,
        indices: &[usize],
        rows: &[TwoScoreRow],
        score: &[f32],
        max_gap: Coord,
        min_length: Coord,
    ) {
        if indices.is_empty() {
            return;
        }
        // (start, end, score) for each significant interval
        let mut content: Vec<(Coord, Coord, f32)> = Vec::with_capacity(indices.len());
        for (n, &i) in indices.iter().enumerate() {
            // start is the previous row's end; row 0 has no previous, so 0
            let mut start = if i == 0 { 0 } else { rows[i - 1].end };
            if n == 0 && indices[0] == 0 {
                start = 0;
            }
            content.push((start, rows[i].end, score[i]));
        }
        let mut group: Vec<(Coord, Coord, f32)> = vec![content[0]];
        for &seg in &content[1..] {
            if seg.0.saturating_sub(group.last().expect("non-empty").1) <= max_gap {
                group.push(seg);
            } else {
                if let Some(p) = self.finish_group(chrom, &group, min_length) {
                    sink.push(p);
                }
                group = vec![seg];
            }
        }
        if let Some(p) = self.finish_group(chrom, &group, min_length) {
            sink.push(p);
        }
    }

    /// Emit one merged region if it is long enough, with its weighted mean score.
    fn finish_group(
        &self,
        chrom: ChromId,
        group: &[(Coord, Coord, f32)],
        min_length: Coord,
    ) -> Option<DiffPeak> {
        let last = group.last().expect("non-empty");
        let first = group.first().expect("non-empty");
        if last.1 - first.0 < min_length {
            return None;
        }
        // `mean_from_peakcontent` (`ScoreTrack.py:1922-1948`): accumulate the
        // weighted sum in f64 for precision, then narrow to f32 for storage.
        let mut ln: u64 = 0;
        let mut sum_v: f64 = 0.0;
        for &(s, e, v) in group {
            sum_v += f64::from(v) * (e - s) as f64;
            ln += u64::from(e - s);
        }
        Some(DiffPeak {
            chrom,
            start: first.0,
            end: last.1,
            score: (sum_v / ln as f64) as f32,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_lr_is_zero_on_the_diagonal() {
        assert_eq!(log_lr_asym(3.0, 3.0), 0.0);
        assert_eq!(log_lr_sym(3.0, 3.0), 0.0);
    }

    /// (x, y, expected f32 bit pattern) triples computed by upstream's own
    /// formulas with `cython.float` arithmetic (see `docs/upstream-findings.md`
    /// F137).
    ///
    /// Pinned as bit patterns rather than decimal literals because the `f32`
    /// narrowing is load-bearing and the shortest decimal that round-trips is
    /// fragile to write by hand.
    const ASYM: [(f32, f32, u32); 5] = [
        (10.0, 2.0, 0x4060_fb66),
        (2.0, 10.0, 0xc004_e3ff),
        (5.0, 5.0, 0x0000_0000),
        (100.0, 3.0, 0x42dc_5297),
        (1.0, 50.0, 0xc19c_a6d4),
    ];

    #[test]
    fn log_lr_asym_matches_upstream() {
        for (x, y, want) in ASYM {
            assert_eq!(log_lr_asym(x, y).to_bits(), want, "asym({x}, {y})");
        }
    }

    #[test]
    fn log_lr_sym_matches_upstream() {
        // note `sym(2, 10) = -3.5153444` but `asym(2, 10) = -2.0764158`: the two
        // are different functions, and upstream's `y > x` branch reuses
        // `log(x) - log(y)` rather than swapping the arguments.
        const SYM: [(f32, f32, u32); 5] = [
            (10.0, 2.0, 0x4060_fb66),
            (2.0, 10.0, 0xc060_fb66),
            (5.0, 5.0, 0x0000_0000),
            (100.0, 3.0, 0x42dc_5297),
            (1.0, 50.0, 0xc27e_ac1b),
        ];
        for (x, y, want) in SYM {
            assert_eq!(log_lr_sym(x, y).to_bits(), want, "sym({x}, {y})");
        }
    }

    #[test]
    fn log_lr_asym_is_not_antisymmetric() {
        // upstream's two branches are different expressions, so swapping the
        // arguments does not merely negate the result. Documented so a future
        // "simplification" to a single symmetric expression is caught.
        let a = log_lr_asym(10.0, 2.0);
        let b = log_lr_asym(2.0, 10.0);
        assert!((a + b).abs() > 1e-3, "expected {a} and {b} to differ");
    }

    #[test]
    fn negative_zero_collapses_to_one_key() {
        assert_eq!(macs_score_key(0.0), macs_score_key(-0.0));
    }

    #[test]
    fn pscore_truncates_treatment_to_int() {
        // 5.9 + 0.0 -> 5 observed, not 6
        assert_eq!(
            pscore_interval(5.9, 1.0, 0.0),
            pscore_interval(5.0, 1.0, 0.0)
        );
    }

    #[test]
    fn pscore_counts_a_sum_that_rounds_up_in_f32() {
        // 19.899999618530273 + 0.10000000149011612 == 20.0 in f32 but
        // 19.999999620020389 in f64, and upstream counts the latter
        let (observed, _) = crate::pseudocounted_inputs(19.9, 5.0, 0.1);
        assert_eq!(observed, 19, "must not round the sum up to 20");
        assert_eq!((19.9f32 + 0.1f32) as i64, 20, "f32 addition would give 20");
    }
}
