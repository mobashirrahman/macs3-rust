//! Port of `MACS3/Signal/PosReadsInfo.py` -- the per-position evidence record and
//! the genotype call that `callvar` writes.
//!
//! # The thing that will bite you: dict insertion order is load-bearing
//!
//! Upstream picks the top two alleles with
//!
//! ```python
//! [self.top1allele, self.top2allele] = sorted(self.n_reads, key=self.n_reads_T.get, reverse=True)[:2]
//! ```
//!
//! `sorted` is **stable**, and `reverse=True` does not reverse equal elements -- it
//! only inverts the comparison. So ties in `n_reads_T` are broken by the *insertion
//! order of the `n_reads` dict*, and that order is fixed by the initialiser:
//!
//! ```python
//! {ref_allele: 0, b'A': 0, b'C': 0, b'G': 0, b'T': 0, b'N': 0, b'*': 0}
//! ```
//!
//! i.e. the reference base comes **first**, then `A C G T N *`, then any further
//! allele in the order `add_T`/`add_C` first saw it. A reference base of `G` gives
//! the order `G A C T N *` -- `G` is a key *before* `A`, not after it. This is why
//! [`AlleleTable`] keeps `Vec` insertion order rather than a `HashMap` iteration
//! order, and why [`AlleleTable::top_two`] sorts a vector of indices with a stable
//! sort instead of collecting from a map.
//!
//! # `cython.float` again
//!
//! `min_top12alleles_ratio`, `max_allowed_ar` and `top12alleles_ratio` are all
//! `cython.float`, so the ratio is computed in `f64` and then **stored as `f32`**
//! before being compared against the `f32` threshold. That is a real rounding step
//! on the value being tested, so it is done rather than elided -- see F213 for the
//! same trap in `VariantStat`.
//!
//! # Thresholds are comparisons, not casts
//!
//! `deltaBIC` is a `double` but `min_delta_BIC` is a `cython.float`, and the
//! `(int)` casts that render `GQ` and the `PL` triple truncate toward zero like C.

use crate::variant_stat::{
    cal_model_heter_as, cal_model_heter_noas, cal_model_homo, DEFAULT_MAX_ALLOWED_AR,
};
use macs_stats::binomial::binomial_cdf;

/// `LN10` as upstream spells it -- not `std::f64::consts::LN_10`, and not
/// `LN10_TENTH * 10`. Both mutations shift every `PL` value in the VCF.
#[allow(clippy::excessive_precision, clippy::approx_constant)]
const LN10: f64 = 2.3025850929940458;

/// An observed allele: `A`, `C`, `G`, `T`, `N`, `*`, or a multi-byte indel string.
///
/// Kept as a byte string rather than an enum because indels really are multi-byte
/// (`*` is a single byte but an insertion is longer) and upstream compares
/// `len(allele) == 1` against these values.
pub type Allele = Vec<u8>;

/// Per-allele counters, in insertion order.
#[derive(Debug, Clone, Default)]
struct Entry {
    allele: Allele,
    bq_t: Vec<i32>,
    bq_c: Vec<i32>,
    n_t: i32,
    n_c: i32,
    n_all: i32,
    n_strand_plus: i32,
    n_strand_minus: i32,
    n_tips: i32,
}

/// An insertion-ordered allele table.
///
/// The order is the contract -- see the module docs. Lookup goes through a side
/// index so the hot `add_*` path stays O(1) without giving up the ordering.
#[derive(Debug, Clone, Default)]
pub struct AlleleTable {
    entries: Vec<Entry>,
}

impl AlleleTable {
    fn find(&self, a: &[u8]) -> Option<usize> {
        self.entries.iter().position(|e| e.allele == a)
    }

    fn get_or_create(&mut self, a: &[u8]) -> usize {
        match self.find(a) {
            Some(i) => i,
            None => {
                self.entries.push(Entry {
                    allele: a.to_vec(),
                    ..Default::default()
                });
                self.entries.len() - 1
            }
        }
    }

    /// `sorted(n_reads, key=n_reads_T.get, reverse=True)[:2]`.
    ///
    /// Rust's `sort_by` is stable, and Python's `sorted(reverse=True)` is stable in
    /// the same sense, so equal treatment counts keep insertion order in both.
    fn top_two(&self) -> Option<(usize, usize)> {
        let mut order: Vec<usize> = (0..self.entries.len()).collect();
        order.sort_by(|&a, &b| self.entries[b].n_t.cmp(&self.entries[a].n_t));
        match order.len() {
            0 => None,
            1 => Some((order[0], order[0])),
            _ => Some((order[0], order[1])),
        }
    }

    fn sum_n_all(&self) -> i32 {
        self.entries.iter().map(|e| e.n_all).sum()
    }

    fn sum_n_t(&self) -> i32 {
        self.entries.iter().map(|e| e.n_t).sum()
    }

    fn sum_n_c(&self) -> i32 {
        self.entries.iter().map(|e| e.n_c).sum()
    }

    /// Zero an allele out entirely, exactly as upstream's `= []` / `= 0` run does.
    fn zero(&mut self, i: usize) {
        let e = &mut self.entries[i];
        e.bq_t.clear();
        e.bq_c.clear();
        e.n_t = 0;
        e.n_c = 0;
        e.n_all = 0;
        e.n_tips = 0;
    }
}

/// Which sample a read came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sample {
    Treatment,
    Control,
}

/// The read depth breakdown `raw_read_depth` can ask for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DepthOpt {
    All,
    Treatment,
    Control,
}

/// The genotype strings upstream can assign.
pub const GT_UNSURE: &str = "unsure";

/// A fully populated `PosReadsInfo`.
#[derive(Debug, Clone)]
pub struct PosReadsInfo {
    pub ref_pos: u64,
    pub ref_allele: Allele,
    pub alt_allele: Vec<u8>,
    pub filterout: bool,
    bq_set_t: AlleleTable,
    bq_set_c: AlleleTable,
    n_reads: AlleleTable,
    top1: usize,
    top2: usize,
    /// `f32` upstream; see the module docs.
    pub top12alleles_ratio: f32,
    pub ln_l_homo_major: f64,
    pub ln_l_heter_as: f64,
    pub ln_l_heter_noas: f64,
    pub ln_l_homo_minor: f64,
    pub bic_homo_major: f64,
    pub bic_heter_as: f64,
    pub bic_heter_noas: f64,
    pub bic_homo_minor: f64,
    pub pl_00: f64,
    pub pl_01: f64,
    pub pl_11: f64,
    pub delta_bic: f64,
    pub gq: f64,
    pub gt: String,
    pub var_type: String,
    pub mutation_type: String,
}

/// Default quality cutoff: "only consider Q20 or read_bq > 20", i.e. `read_bq <= 20`
/// is dropped.
const DEFAULT_Q_CUTOFF: i32 = 20;

impl PosReadsInfo {
    /// Upstream's `__init__`.
    pub fn new(ref_pos: u64, ref_allele: &[u8]) -> Self {
        let mut n_reads = AlleleTable::default();
        // Insertion order is the contract: ref_allele first, then A C G T N *.
        n_reads.get_or_create(ref_allele);
        for a in [b"A", b"C", b"G", b"T", b"N", b"*"] {
            n_reads.get_or_create(a);
        }
        Self {
            ref_pos,
            ref_allele: ref_allele.to_vec(),
            alt_allele: b".".to_vec(),
            filterout: false,
            bq_set_t: n_reads.clone(),
            bq_set_c: n_reads.clone(),
            n_reads,
            // Placeholders; `update_top_alleles` sets the real ones before use.
            top1: 0,
            top2: 0,
            top12alleles_ratio: 0.0,
            ln_l_homo_major: 0.0,
            ln_l_heter_as: 0.0,
            ln_l_heter_noas: 0.0,
            ln_l_homo_minor: 0.0,
            bic_homo_major: 0.0,
            bic_heter_as: 0.0,
            bic_heter_noas: 0.0,
            bic_homo_minor: 0.0,
            pl_00: 0.0,
            pl_01: 0.0,
            pl_11: 0.0,
            delta_bic: 0.0,
            gq: 0.0,
            gt: GT_UNSURE.to_string(),
            var_type: String::new(),
            mutation_type: String::new(),
        }
    }

    /// `add_T`. `strand` is 0 for plus, 1 for minus.
    ///
    /// Note that only `add_T` touches `n_strand`; the control contributes to
    /// `n_reads_C` and `n_reads` but never to the strand tally. That is why `SB` in
    /// the VCF can exceed neither `DPT`'s two alleles nor be non-zero where `DPC` is.
    pub fn add_t(&mut self, read_allele: &[u8], read_bq: i32, strand: u8, tip: bool) {
        self.add_read(Sample::Treatment, read_allele, read_bq, strand, tip);
    }

    /// `add_C`.
    pub fn add_c(&mut self, read_allele: &[u8], read_bq: i32) {
        self.add_read(Sample::Control, read_allele, read_bq, 0, false);
    }

    fn add_read(&mut self, which: Sample, read_allele: &[u8], read_bq: i32, strand: u8, tip: bool) {
        if read_bq <= DEFAULT_Q_CUTOFF {
            return;
        }
        // `if not self.n_reads.has_key(read_allele)` creates the allele in all five
        // tables at once.
        let i = self.n_reads.get_or_create(read_allele);
        if self.bq_set_t.find(read_allele).is_none() {
            self.bq_set_t.get_or_create(read_allele);
            self.bq_set_c.get_or_create(read_allele);
        }

        match which {
            Sample::Treatment => {
                self.bq_set_t.entries[i].bq_t.push(read_bq);
                self.bq_set_t.entries[i].n_t += 1;
                self.n_reads.entries[i].n_t += 1;
                self.n_reads.entries[i].n_all += 1;
                if strand == 0 {
                    self.n_reads.entries[i].n_strand_plus += 1;
                } else {
                    self.n_reads.entries[i].n_strand_minus += 1;
                }
                if tip {
                    self.n_reads.entries[i].n_tips += 1;
                }
            }
            Sample::Control => {
                self.bq_set_c.entries[i].bq_c.push(read_bq);
                self.bq_set_c.entries[i].n_c += 1;
                self.n_reads.entries[i].n_c += 1;
                self.n_reads.entries[i].n_all += 1;
            }
        }
    }

    /// `raw_read_depth`.
    pub fn raw_read_depth(&self, opt: DepthOpt) -> i32 {
        match opt {
            DepthOpt::All => self.n_reads.sum_n_all(),
            DepthOpt::Treatment => self.n_reads.sum_n_t(),
            DepthOpt::Control => self.n_reads.sum_n_c(),
        }
    }

    fn ref_index(&self) -> usize {
        self.n_reads.find(&self.ref_allele).unwrap_or(0)
    }

    fn allele(&self, i: usize) -> &[u8] {
        &self.n_reads.entries[i].allele
    }

    pub fn top1_allele(&self) -> &[u8] {
        self.allele(self.top1)
    }

    pub fn top2_allele(&self) -> &[u8] {
        self.allele(self.top2)
    }

    pub fn n_t(&self, allele: &[u8]) -> i32 {
        self.n_reads
            .find(allele)
            .map_or(0, |i| self.n_reads.entries[i].n_t)
    }

    pub fn n_c(&self, allele: &[u8]) -> i32 {
        self.n_reads
            .find(allele)
            .map_or(0, |i| self.n_reads.entries[i].n_c)
    }

    pub fn strand_counts(&self, allele: &[u8]) -> (i32, i32) {
        self.n_reads.find(allele).map_or((0, 0), |i| {
            (
                self.n_reads.entries[i].n_strand_plus,
                self.n_reads.entries[i].n_strand_minus,
            )
        })
    }

    /// `update_top_alleles`. Returns `true` if the position is still live.
    ///
    /// `min_top12alleles_ratio` and `max_allowed_ar` are `cython.float`.
    #[allow(clippy::excessive_precision)]
    pub fn update_top_alleles(
        &mut self,
        min_top12alleles_ratio: f32,
        min_altallele_count: i32,
        max_allowed_ar: f32,
    ) -> bool {
        let Some((t1, t2)) = self.n_reads.top_two() else {
            self.filterout = true;
            return false;
        };
        self.top1 = t1;
        self.top2 = t2;

        if self.n_reads.entries[t1].n_t + self.n_reads.entries[t2].n_t == 0 {
            self.filterout = true;
            return false;
        }

        // Upstream:
        //   (len(top1)==1 and len(top2)==1) and (
        //       (top2 != ref and (n_t[top2]-tips[top2]) < min_alt) or ar > max_ar )
        let both_single = self.allele(t1).len() == 1 && self.allele(t2).len() == 1;
        if both_single {
            let a = {
                let e1 = &self.n_reads.entries[t1];
                let e2 = &self.n_reads.entries[t2];
                let ar = e1.n_t as f64 / (e1.n_t + e2.n_t) as f64;
                (self.allele(t2) != self.ref_allele.as_slice()
                    && (e2.n_t - e2.n_tips) < min_altallele_count)
                    || ar > max_allowed_ar as f64
            };
            if a {
                self.n_reads.zero(t2);
                self.bq_set_t.zero(t2);
                self.bq_set_c.zero(t2);
                let e1 = &self.n_reads.entries[t1];
                if self.allele(t1) != self.ref_allele.as_slice()
                    && (e1.n_t - e1.n_tips) < min_altallele_count
                {
                    self.n_reads.zero(t1);
                    self.bq_set_t.zero(t1);
                    self.bq_set_c.zero(t1);
                }
            }
        }

        if self.n_reads.entries[t1].n_t + self.n_reads.entries[t2].n_t == 0 {
            self.filterout = true;
            return false;
        }

        let total = self.n_reads.sum_n_all();
        // Computed in f64, stored to f32, compared against an f32 threshold.
        self.top12alleles_ratio = ((self.n_reads.entries[t1].n_all + self.n_reads.entries[t2].n_all)
            as f64
            / total as f64) as f32;
        if self.top12alleles_ratio < min_top12alleles_ratio {
            self.filterout = true;
            return false;
        }

        if self.allele(t1) == self.ref_allele.as_slice() && self.n_reads.entries[t2].n_all == 0 {
            self.var_type = "homo_ref".to_string();
            self.filterout = true;
            return false;
        }
        true
    }

    /// `call_GT`. Requires `update_top_alleles` first.
    ///
    /// `Err` propagates the domain error from `VariantStat` -- upstream raises there
    /// and aborts the run.
    #[allow(clippy::excessive_precision)]
    pub fn call_gt(&mut self, max_allowed_ar: f32) -> Result<(), String> {
        if self.filterout {
            return Ok(());
        }
        let (t1, t2) = (self.top1, self.top2);
        let top1_bq_t = std::mem::take(&mut self.bq_set_t.entries[t1].bq_t);
        let top2_bq_t = std::mem::take(&mut self.bq_set_t.entries[t2].bq_t);
        let top1_bq_c = std::mem::take(&mut self.bq_set_c.entries[t1].bq_c);
        let top2_bq_c = std::mem::take(&mut self.bq_set_c.entries[t2].bq_c);

        let r = self.score_models(
            &top1_bq_t,
            &top1_bq_c,
            &top2_bq_t,
            &top2_bq_c,
            max_allowed_ar,
        );

        self.bq_set_t.entries[t1].bq_t = top1_bq_t;
        self.bq_set_t.entries[t2].bq_t = top2_bq_t;
        self.bq_set_c.entries[t1].bq_c = top1_bq_c;
        self.bq_set_c.entries[t2].bq_c = top2_bq_c;

        r?;

        let ref_allele = self.ref_allele.clone();
        let top1_allele = self.allele(t1).to_vec();
        let top2_allele = self.allele(t2).to_vec();
        let n_top2_all = self.n_reads.entries[t2].n_all;

        if top1_allele != ref_allele && n_top2_all == 0 {
            // No minor allele anywhere: assume 1/1, keep it only if BIC agrees.
            self.delta_bic = self
                .bic_heter_noas
                .min(self.bic_heter_as)
                .min(self.bic_homo_minor)
                - self.bic_homo_major;
            if self.delta_bic < 2.0 {
                self.filterout = true;
                return Ok(());
            }
            self.var_type = "homo".to_string();
            self.gt = "1/1".to_string();
            self.pl_00 = -10.0 * self.ln_l_homo_minor / LN10;
            self.pl_01 = -10.0 * self.ln_l_heter_noas.max(self.ln_l_heter_as) / LN10;
            self.pl_11 = -10.0 * self.ln_l_homo_major / LN10;
            // This branch -- and only this one -- clamps the two losers at 0.
            self.pl_00 = (self.pl_00 - self.pl_11).max(0.0);
            self.pl_01 = (self.pl_01 - self.pl_11).max(0.0);
            self.pl_11 = 0.0;
            self.gq = self.pl_00.min(self.pl_01);
            self.alt_allele = top1_allele;
        } else {
            self.assign_genotype(&ref_allele, &top1_allele, &top2_allele);
        }

        self.set_mutation_type();
        Ok(())
    }

    /// The five-model comparison ladder. Each branch is guarded by `BIC + 2 <=`,
    /// i.e. a 2-nat unit preference, and `deltaBIC` is always measured against the
    /// best alternative.
    fn assign_genotype(&mut self, ref_allele: &[u8], top1: &[u8], top2: &[u8]) {
        let (hm, hn, na, as_) = (
            self.bic_homo_major,
            self.bic_homo_minor,
            self.bic_heter_noas,
            self.bic_heter_as,
        );

        if ref_allele != top1 && hm + 2.0 <= hn && hm + 2.0 <= na && hm + 2.0 <= as_ {
            self.var_type = "homo".to_string();
            self.delta_bic = na.min(as_).min(hn) - hm;
            self.gt = "1/1".to_string();
            self.alt_allele = top1.to_vec();
            self.pl_00 = -10.0 * self.ln_l_homo_minor / LN10;
            self.pl_01 = -10.0 * fmax(self.ln_l_heter_noas, self.ln_l_heter_as) / LN10;
            self.pl_11 = -10.0 * self.ln_l_homo_major / LN10;
            self.pl_00 -= self.pl_11;
            self.pl_01 -= self.pl_11;
            self.pl_11 = 0.0;
            self.gq = self.pl_00.min(self.pl_01);
        } else if na + 2.0 <= hm && na + 2.0 <= hn && na + 2.0 <= as_ {
            self.var_type = "heter_noAS".to_string();
            self.delta_bic = hm.min(hn) - na;
            self.pl_00 = -10.0 * self.ln_l_homo_minor / LN10;
            self.pl_01 = -10.0 * self.ln_l_heter_noas / LN10;
            self.pl_11 = -10.0 * self.ln_l_homo_major / LN10;
            self.pl_00 -= self.pl_01;
            self.pl_11 -= self.pl_01;
            self.pl_01 = 0.0;
            self.gq = self.pl_00.min(self.pl_11);
        } else if as_ + 2.0 <= hm && as_ + 2.0 <= hn && as_ + 2.0 <= na {
            self.var_type = "heter_AS".to_string();
            self.delta_bic = hm.min(hn) - as_;
            self.pl_00 = -10.0 * self.ln_l_homo_minor / LN10;
            self.pl_01 = -10.0 * self.ln_l_heter_as / LN10;
            self.pl_11 = -10.0 * self.ln_l_homo_major / LN10;
            self.pl_00 -= self.pl_01;
            self.pl_11 -= self.pl_01;
            self.pl_01 = 0.0;
            self.gq = self.pl_00.min(self.pl_11);
        } else if as_ + 2.0 <= hm && as_ + 2.0 <= hn {
            self.var_type = "heter_unsure".to_string();
            // Measured against whichever heter model scored *better* (lower BIC).
            self.delta_bic = hm.min(hn) - fmax(as_, na);
            self.pl_00 = -10.0 * self.ln_l_homo_minor / LN10;
            self.pl_01 = -10.0 * fmax(self.ln_l_heter_noas, self.ln_l_heter_as) / LN10;
            self.pl_11 = -10.0 * self.ln_l_homo_major / LN10;
            self.pl_00 -= self.pl_01;
            self.pl_11 -= self.pl_01;
            self.pl_01 = 0.0;
            self.gq = self.pl_00.min(self.pl_11);
        } else if ref_allele == top1 && hm < hn && hm < na && hm < as_ {
            self.var_type = "homo_ref".to_string();
            self.gt = "0/0".to_string();
            self.filterout = true;
        } else {
            self.var_type = "unsure".to_string();
            self.filterout = true;
        }

        // Evaluated *after* the ladder, at the same level as it: upstream asks
        // `self.type.startswith("heter")` once a type has actually been chosen. Reading
        // it beforehand sees the empty initialiser and silently drops every
        // heterozygous call.
        if self.var_type.starts_with("heter") {
            if ref_allele == top1 {
                self.alt_allele = top2.to_vec();
                self.gt = "0/1".to_string();
            } else if ref_allele == top2 {
                self.alt_allele = top1.to_vec();
                self.gt = "0/1".to_string();
            } else {
                let mut joined = top1.to_vec();
                joined.push(b',');
                joined.extend_from_slice(top2);
                self.alt_allele = joined;
                self.gt = "1/2".to_string();
            }
        }
    }

    fn score_models(
        &mut self,
        top1_bq_t: &[i32],
        top1_bq_c: &[i32],
        top2_bq_t: &[i32],
        top2_bq_c: &[i32],
        max_allowed_ar: f32,
    ) -> Result<(), String> {
        let (hm, hbm) = cal_model_homo(top1_bq_t, top1_bq_c, top2_bq_t, top2_bq_c)?;
        self.ln_l_homo_major = hm;
        self.bic_homo_major = hbm;
        let (hmin, hbmin) = cal_model_homo(top2_bq_t, top2_bq_c, top1_bq_t, top1_bq_c)?;
        self.ln_l_homo_minor = hmin;
        self.bic_homo_minor = hbmin;
        let (n, bn) = cal_model_heter_noas(top1_bq_t, top1_bq_c, top2_bq_t, top2_bq_c)?;
        self.ln_l_heter_noas = n;
        self.bic_heter_noas = bn;
        let (a, ba) =
            cal_model_heter_as(top1_bq_t, top1_bq_c, top2_bq_t, top2_bq_c, max_allowed_ar)?;
        self.ln_l_heter_as = a;
        self.bic_heter_as = ba;
        Ok(())
    }

    fn set_mutation_type(&mut self) {
        let mut parts: Vec<&str> = Vec::new();
        for alt in self.alt_allele.split(|&b| b == b',') {
            parts.push(if alt == b"*" {
                "Deletion"
            } else if alt.len() > 1 {
                "Insertion"
            } else {
                "SNV"
            });
        }
        self.mutation_type = parts.join(",");
    }

    /// `apply_deltaBIC_cutoff`. `min_delta_bic` is a `cython.float`.
    #[allow(clippy::excessive_precision)]
    pub fn apply_delta_bic_cutoff(&mut self, min_delta_bic: f32) {
        if self.filterout {
            return;
        }
        if self.delta_bic < min_delta_bic as f64 {
            self.filterout = true;
        }
    }

    /// `apply_GQ_cutoff`.
    ///
    /// The comparison is `double < int`, so Python promotes the threshold to `f64` and
    /// compares **as floats**. Casting `gq` to `i32` first would be wrong: a `GQ` of
    /// 99.9 passes `>= 100` in Python and fails it after truncation.
    #[allow(clippy::if_same_then_else)]
    pub fn apply_gq_cutoff(&mut self, min_homo_gq: i32, min_heter_gq: i32) {
        if self.filterout {
            return;
        }
        let homo = self.var_type.starts_with("homo");
        let heter = self.var_type.starts_with("heter");
        if homo && self.gq < min_homo_gq as f64 {
            self.filterout = true;
        } else if heter && self.gq < min_heter_gq as f64 {
            self.filterout = true;
        }
    }

    /// `SB_score_ChIP`: the strand-bias score. `>= 1` means filter out.
    ///
    /// `a`/`b` are the major/minor allele on the plus strand, `c`/`d` on the minus.
    pub fn sb_score_chip(&self, a: i32, b: i32, c: i32, d: i32) -> f32 {
        if a + b == 0 || c + d == 0 {
            return 0.0;
        }
        self.sb_from_cdfs(a, b, c, d)
    }

    /// `SB_score_ATAC`. Byte-for-byte the same body as [`Self::sb_score_chip`]:
    /// upstream keeps two verbatim copies of it. Kept as two functions so that if the
    /// two ever diverge upstream, this port has somewhere to diverge too.
    #[allow(clippy::if_same_then_else)]
    pub fn sb_score_atac(&self, a: i32, b: i32, c: i32, d: i32) -> f32 {
        if a + b == 0 || c + d == 0 {
            return 0.0;
        }
        self.sb_from_cdfs(a, b, c, d)
    }

    fn sb_from_cdfs(&self, a: i32, b: i32, c: i32, d: i32) -> f32 {
        // `lower=True` everywhere: upstream passes it explicitly in all four calls.
        let p1_l = binomial_cdf(a as i64, a as i64 + c as i64, 0.5, true);
        let p1_r = binomial_cdf(c as i64, a as i64 + c as i64, 0.5, true);
        let p2_l = binomial_cdf(b as i64, b as i64 + d as i64, 0.5, true);
        let p2_r = binomial_cdf(d as i64, b as i64 + d as i64, 0.5, true);
        if (p1_l < 0.05 && p2_r < 0.05) || (p1_r < 0.05 && p2_l < 0.05) {
            1.0
        } else {
            0.0
        }
    }

    /// The `INFO` and sample columns of `to_vcf`.
    ///
    /// `%.2f` fields and the `%d` truncations are reproduced with `format!`, which
    /// rounds `%f` half-to-even exactly like C's `printf` on the same target.
    pub fn to_vcf(&self) -> String {
        let t1 = self.top1_allele().to_vec();
        let t2 = self.top2_allele().to_vec();
        let (s1p, s1m) = self.strand_counts(&t1);
        let (s2p, s2m) = self.strand_counts(&t2);
        let dpt = self.n_reads.sum_n_t();
        let dpc = self.n_reads.sum_n_c();
        let n1t = self.n_t(&t1);
        let n2t = self.n_t(&t2);
        let n1c = self.n_c(&t1);
        let n2c = self.n_c(&t2);
        let ar = n1t as f64 / (n1t + n2t) as f64;

        let info = format!(
            "M={};MT={};DPT={};DPC={};DP1T={}{};DP2T={}{};DP1C={}{};DP2C={}{};\
             SB={},{},{},{};DBIC={:.2};BICHOMOMAJOR={:.2};BICHOMOMINOR={:.2};\
             BICHETERNOAS={:.2};BICHETERAS={:.2};AR={:.2}",
            self.var_type,
            self.mutation_type,
            dpt,
            dpc,
            n1t,
            String::from_utf8_lossy(&t1),
            n2t,
            String::from_utf8_lossy(&t2),
            n1c,
            String::from_utf8_lossy(&t1),
            n2c,
            String::from_utf8_lossy(&t2),
            s1p,
            s2p,
            s1m,
            s2m,
            self.delta_bic,
            self.bic_homo_major,
            self.bic_homo_minor,
            self.bic_heter_noas,
            self.bic_heter_as,
            ar
        );

        let sample = format!(
            "{}:{}:{}:{},{},{}",
            self.gt,
            self.raw_read_depth(DepthOpt::All),
            self.gq as i64,
            self.pl_00 as i64,
            self.pl_01 as i64,
            self.pl_11 as i64
        );

        [
            String::from_utf8_lossy(&self.ref_allele).into_owned(),
            String::from_utf8_lossy(&self.alt_allele).into_owned(),
            format!("{}", self.gq as i64),
            ".".to_string(),
            info,
            "GT:DP:GQ:PL".to_string(),
            sample,
        ]
        .join("\t")
    }

    /// Convenience for tests and for the orchestration layer: the reference index.
    pub fn ref_allele_index(&self) -> usize {
        self.ref_index()
    }

    /// The default allele-ratio clamp `call_GT` is given by the command line.
    pub fn default_max_allowed_ar() -> f32 {
        DEFAULT_MAX_ALLOWED_AR
    }
}

/// `max` written out because `f64::max` propagates the *NaN* payload differently
/// from Python's `max`, which returns the other operand. These values are never NaN
/// in practice (a NaN would mean `VariantStat` returned one, which cannot happen),
/// but the comparison direction is kept so a future NaN behaves like upstream's.
#[inline]
fn fmax(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    /// Feed `n` reads of one allele at a fixed quality.
    fn feed(p: &mut PosReadsInfo, allele: &[u8], n: i32, bq: i32, plus: bool) {
        for _ in 0..n {
            p.add_t(allele, bq, u8::from(!plus), false);
        }
    }

    #[test]
    fn insertion_order_puts_the_reference_allele_first() {
        // Upstream's dict literal is {ref: 0, A, C, G, T, N, *}. With ref = G the
        // order must be G A C T N *, NOT A C G T N * -- the key appears where the
        // literal puts it.
        let p = PosReadsInfo::new(100, b"G");
        let order: Vec<&[u8]> = (0..p.n_reads.entries.len()).map(|i| p.allele(i)).collect();
        assert_eq!(order, vec![b"G".as_slice(), b"A", b"C", b"T", b"N", b"*"]);
    }

    #[test]
    fn top_two_breaks_ties_by_insertion_order() {
        // Two alleles, equal treatment counts. `sorted(reverse=True)` is stable, so
        // the first in dict order wins -- here the reference base.
        let mut p = PosReadsInfo::new(1, b"G");
        feed(&mut p, b"G", 3, 30, true);
        feed(&mut p, b"A", 3, 30, true);
        assert!(p.update_top_alleles(0.0, 2, 0.95));
        assert_eq!(p.top1_allele(), b"G");
        assert_eq!(p.top2_allele(), b"A");
    }

    #[test]
    fn sorted_covers_every_allele_so_top2_is_often_a_zero_count_one() {
        // `sorted` runs over `n_reads`, which holds all seven seeded alleles, and keys
        // on `n_reads_T`. So the 50 control G reads rank *below* A, C, G, T, N and
        // `*` all at zero treatment count -- and among those ties the first in dict
        // order (C, since A is the reference here) becomes top2.
        //
        // That is upstream's behaviour, not a porting slip, and it is why
        // `n_reads[top2] == 0` is tested for at all: with no second allele, a
        // reference-major position is `homo_ref` and gets filtered.
        let mut p = PosReadsInfo::new(1, b"A");
        feed(&mut p, b"A", 2, 30, true);
        for _ in 0..50 {
            p.add_c(b"G", 30);
        }
        assert!(!p.update_top_alleles(0.0, 2, 0.95));
        assert_eq!(p.var_type, "homo_ref");
        assert!(p.filterout);
    }

    #[test]
    fn quality_cutoff_drops_q20_and_below() {
        let mut p = PosReadsInfo::new(1, b"A");
        p.add_t(b"A", 20, 0, false);
        p.add_t(b"A", 21, 0, false);
        assert_eq!(p.n_t(b"A"), 1, "Q20 is dropped; only Q21 survives");
    }

    #[test]
    fn strand_counts_are_treatment_only() {
        let mut p = PosReadsInfo::new(1, b"A");
        p.add_t(b"A", 30, 0, false);
        p.add_t(b"A", 30, 1, false);
        for _ in 0..5 {
            p.add_c(b"A", 30);
        }
        assert_eq!(
            p.strand_counts(b"A"),
            (1, 1),
            "control must not touch n_strand"
        );
        assert_eq!(p.raw_read_depth(DepthOpt::All), 7);
        assert_eq!(p.raw_read_depth(DepthOpt::Treatment), 2);
        assert_eq!(p.raw_read_depth(DepthOpt::Control), 5);
    }

    #[test]
    fn tip_reads_are_tracked_and_can_zero_the_major_allele() {
        // top1 = G with 2 reads, one of them a tip; top2 = the reference A with none.
        // The allele ratio 2/2 = 1.0 exceeds `max_allowed_ar`, so top2 is zeroed; and
        // because G is not the reference and its non-tip count is 1 < 2, G is zeroed
        // too. Nothing is left, so the position is filtered.
        let mut p = PosReadsInfo::new(1, b"A");
        p.add_t(b"G", 30, 0, false); // real
        p.add_t(b"G", 30, 0, true); // tip
        assert!(!p.update_top_alleles(0.0, 2, 0.95));
        assert_eq!(p.n_t(b"G"), 0);
        assert_eq!(p.n_t(b"A"), 0);
        assert!(p.filterout);
    }

    #[test]
    fn a_tip_only_matters_once_the_allele_is_already_being_zeroed() {
        // The tip rule is *nested inside* the `ar > max_allowed_ar` / low-alt-count
        // condition, so a tip on its own never zeroes anything. Here G is 3 reads with
        // 1 tip and A has none: the ratio 3/3 trips the clamp, A is zeroed, and G then
        // survives because 3 - 1 = 2 is not `< min_altallele_count` (2).
        let mut p = PosReadsInfo::new(1, b"A");
        for _ in 0..2 {
            p.add_t(b"G", 30, 0, false);
        }
        p.add_t(b"G", 30, 0, true); // tip
        assert!(p.update_top_alleles(0.0, 2, 0.95));
        assert_eq!(p.top1_allele(), b"G");
        assert_eq!(p.n_t(b"G"), 3, "one tip does not clear a count of 3");
    }

    #[test]
    fn top12_ratio_is_computed_over_all_reads() {
        let mut p = PosReadsInfo::new(1, b"A");
        feed(&mut p, b"A", 6, 30, true);
        feed(&mut p, b"G", 2, 30, true);
        assert!(p.update_top_alleles(0.0, 2, 0.95));
        // 8 of 8 reads are on the two top alleles.
        assert!((p.top12alleles_ratio - 1.0f32).abs() < 1e-6);
    }

    #[test]
    fn ratio_threshold_rejects_a_third_allele() {
        let mut p = PosReadsInfo::new(1, b"A");
        feed(&mut p, b"A", 5, 30, true);
        feed(&mut p, b"G", 3, 30, true);
        for _ in 0..10 {
            p.add_c(b"C", 30);
        }
        // Only 8 of 18 reads are on the top two alleles -> 0.444 < 0.8.
        assert!(!p.update_top_alleles(0.8, 2, 0.95));
        assert!(p.filterout);
    }

    #[test]
    fn reference_only_position_is_homo_ref_and_filtered() {
        let mut p = PosReadsInfo::new(1, b"A");
        feed(&mut p, b"A", 5, 30, true);
        assert!(!p.update_top_alleles(0.0, 2, 0.95));
        assert_eq!(p.var_type, "homo_ref");
        assert!(p.filterout);
    }

    #[test]
    fn a_strong_single_alt_allele_calls_homo_1_1() {
        let mut p = PosReadsInfo::new(4242, b"A");
        for _ in 0..20 {
            p.add_t(b"G", 35, 0, false);
        }
        assert!(p.update_top_alleles(0.8, 2, 0.99));
        p.call_gt(0.99).unwrap();
        assert_eq!(p.var_type, "homo");
        assert_eq!(p.gt, "1/1");
        assert_eq!(p.alt_allele, b"G".to_vec());
        assert_eq!(p.mutation_type, "SNV");
        // The 1/1 branch is the only one that clamps the losers at zero.
        assert_eq!(p.pl_11, 0.0);
        assert!(p.pl_00 >= 0.0 && p.pl_01 >= 0.0);
    }

    #[test]
    fn a_balanced_mixture_calls_heterozygous() {
        let mut p = PosReadsInfo::new(7, b"A");
        feed(&mut p, b"G", 10, 30, true);
        feed(&mut p, b"A", 10, 30, false);
        assert!(p.update_top_alleles(0.8, 2, 0.99));
        p.call_gt(0.99).unwrap();
        assert!(
            p.var_type.starts_with("heter"),
            "expected a heter call, got {}",
            p.var_type
        );
        assert_eq!(p.gt, "0/1");
        assert_eq!(p.alt_allele, b"G".to_vec());
        assert_eq!(p.pl_01, 0.0, "the heterozygote likelihood is the reference");
    }

    #[test]
    fn neither_ref_nor_top1_calls_the_multiallelic_genotype() {
        let mut p = PosReadsInfo::new(9, b"T");
        feed(&mut p, b"A", 10, 30, true);
        feed(&mut p, b"G", 10, 30, false);
        assert!(p.update_top_alleles(0.8, 2, 0.99));
        p.call_gt(0.99).unwrap();
        assert_eq!(p.gt, "1/2");
        assert_eq!(p.alt_allele, b"A,G".to_vec());
        assert_eq!(p.mutation_type, "SNV,SNV");
    }

    #[test]
    fn insertion_and_deletion_mutation_types() {
        // A balanced mixture with a reference that is neither of the two observed
        // alleles: that is what produces the multiallelic `1/2` genotype, and it is
        // the only way to get past the 2-unit BIC margin without a control sample.
        let mut p = PosReadsInfo::new(1, b"T");
        // 11, not 10: an exact tie would be broken by dict insertion order, and `A`
        // is seeded while `AGGG` is appended, so `A` would win and the ALT order (and
        // therefore the mutation-type order) would come out reversed.
        for _ in 0..11 {
            p.add_t(b"AGGG", 30, 0, false);
        }
        for _ in 0..10 {
            p.add_t(b"A", 30, 1, false);
        }
        assert!(p.update_top_alleles(0.8, 2, 0.99));
        p.call_gt(0.99).unwrap();
        assert!(p.var_type.starts_with("heter"), "got {}", p.var_type);
        assert_eq!(p.gt, "1/2");
        assert_eq!(p.alt_allele, b"AGGG,A".to_vec());
        assert_eq!(p.mutation_type, "Insertion,SNV");

        let mut d = PosReadsInfo::new(2, b"T");
        for _ in 0..11 {
            d.add_t(b"*", 30, 0, false);
        }
        for _ in 0..10 {
            d.add_t(b"A", 30, 1, false);
        }
        assert!(d.update_top_alleles(0.8, 2, 0.99));
        d.call_gt(0.99).unwrap();
        assert_eq!(d.gt, "1/2");
        assert_eq!(d.alt_allele, b"*,A".to_vec());
        assert_eq!(d.mutation_type, "Deletion,SNV");
    }

    #[test]
    fn vcf_line_has_every_field_upstream_writes() {
        let mut p = PosReadsInfo::new(4242, b"A");
        feed(&mut p, b"G", 12, 33, true);
        feed(&mut p, b"G", 8, 33, false);
        feed(&mut p, b"A", 5, 30, true);
        assert!(p.update_top_alleles(0.8, 2, 0.99));
        p.call_gt(0.99).unwrap();
        let line = p.to_vcf();
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 7, "{line}");
        assert_eq!(f[0], "A");
        assert_eq!(f[1], "G");
        assert_eq!(f[3], ".");
        assert_eq!(f[5], "GT:DP:GQ:PL");
        for key in [
            "M=",
            "MT=",
            "DPT=",
            "DPC=",
            "DP1T=",
            "DP2T=",
            "DP1C=",
            "DP2C=",
            "SB=",
            "DBIC=",
            "BICHOMOMAJOR=",
            "BICHOMOMINOR=",
            "BICHETERNOAS=",
            "BICHETERAS=",
            "AR=",
        ] {
            assert!(f[4].contains(key), "INFO missing {key}: {}", f[4]);
        }
        // DPT counts only treatment reads (12+8 G plus 5 A = 25); DPC only control.
        assert!(f[4].contains("DPT=25;"), "{}", f[4]);
        assert!(f[4].contains("DPC=0;"), "{}", f[4]);
        // The four SB counts must partition the treatment depth.
        let sb: Vec<i64> = f[4]
            .split("SB=")
            .nth(1)
            .unwrap()
            .split(';')
            .next()
            .unwrap()
            .split(',')
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(sb.iter().sum::<i64>(), 25, "SB must sum to DPT: {}", f[4]);
        // PL must have exactly three entries and end with the reference likelihood 0.
        let sample = f[6];
        let (gt, rest) = sample.split_once(':').unwrap();
        let parts: Vec<&str> = rest.split(':').collect();
        assert_eq!(gt, p.gt);
        // rest is "DP:GQ:PL" -- three fields, the genotype having been split off.
        assert_eq!(parts.len(), 3, "{sample}");
        assert_eq!(parts[1].parse::<i64>().unwrap(), p.gq as i64);
        assert_eq!(parts[2].split(',').count(), 3);
        // Whichever model won, its own likelihood is the reference and lands at 0 --
        // the middle entry for every `heter_*` call, the last entry for a `homo` 1/1.
        // All 16 records of the oracle golden read like `159,0,58`, which is this rule.
        let pl: Vec<&str> = parts[2].split(',').collect();
        let zeroed = if p.var_type.starts_with("heter") {
            pl[1]
        } else {
            pl[2]
        };
        assert_eq!(
            zeroed, "0",
            "PL should be zeroed at the winning model: {sample}"
        );
    }

    #[test]
    fn ar_field_is_the_treatment_allele_ratio() {
        let mut p = PosReadsInfo::new(1, b"A");
        feed(&mut p, b"G", 7, 30, true);
        feed(&mut p, b"A", 3, 30, false);
        assert!(p.update_top_alleles(0.8, 2, 0.99));
        p.call_gt(0.99).unwrap();
        assert!(p.to_vcf().contains("AR=0.70"), "{}", p.to_vcf());
    }

    #[test]
    fn strand_bias_score_is_zero_when_one_strand_is_empty() {
        let p = PosReadsInfo::new(1, b"A");
        assert_eq!(p.sb_score_chip(5, 2, 0, 0), 0.0);
        assert_eq!(p.sb_score_chip(0, 0, 5, 2), 0.0);
        assert_eq!(p.sb_score_atac(5, 2, 0, 0), 0.0);
    }

    #[test]
    fn strand_bias_score_rejects_opposing_biases() {
        // top1 on the minus strand, top2 on the plus strand, both strongly -- exactly
        // the shape upstream calls a false positive.
        let p = PosReadsInfo::new(1, b"A");
        assert_eq!(p.sb_score_chip(0, 20, 20, 0), 1.0);
        // Consistent bias (both alleles favour the same strand) is allowed.
        assert_eq!(p.sb_score_chip(20, 18, 0, 0), 0.0);
    }

    #[test]
    fn gq_cutoff_compares_as_floats_not_as_truncated_ints() {
        let mut p = PosReadsInfo::new(1, b"A");
        feed(&mut p, b"G", 5, 30, true);
        assert!(p.update_top_alleles(0.0, 2, 0.99));
        p.call_gt(0.99).unwrap();
        p.gq = 99.9;
        p.apply_gq_cutoff(100, 100);
        assert!(p.filterout, "99.9 < 100 in Python, so it must be filtered");
        let mut q = PosReadsInfo::new(1, b"A");
        feed(&mut q, b"G", 5, 30, true);
        assert!(q.update_top_alleles(0.0, 2, 0.99));
        q.call_gt(0.99).unwrap();
        q.gq = 100.5;
        q.apply_gq_cutoff(100, 100);
        assert!(!q.filterout, "100.5 >= 100 must survive");
    }

    #[test]
    fn delta_bic_cutoff_uses_the_given_threshold() {
        let mut p = PosReadsInfo::new(1, b"A");
        p.delta_bic = 9.99;
        p.apply_delta_bic_cutoff(10.0);
        assert!(p.filterout);
        let mut q = PosReadsInfo::new(1, b"A");
        q.delta_bic = 10.01;
        q.apply_delta_bic_cutoff(10.0);
        assert!(!q.filterout);
    }

    #[test]
    fn the_quality_cutoff_keeps_a_zero_quality_base_out_of_the_model() {
        // `add_T` drops `read_bq <= Q` (Q = 20) before anything else, so a Q0 base
        // never reaches `VariantStat`. That matters because `log1p(-1.0)` is a domain
        // error there -- the error path itself is pinned by the 410-case golden in
        // tests/variant_stat_golden.rs, and this test pins the guard that keeps
        // `PosReadsInfo` off it.
        let mut p = PosReadsInfo::new(1, b"A");
        p.add_t(b"G", 0, 0, false);
        assert_eq!(p.raw_read_depth(DepthOpt::Treatment), 0);
        assert_eq!(p.n_t(b"G"), 0);
    }

    #[test]
    fn filtered_positions_are_never_scored() {
        let mut p = PosReadsInfo::new(1, b"A");
        assert!(!p.update_top_alleles(0.0, 2, 0.99)); // no reads at all
        p.call_gt(0.99).unwrap();
        assert_eq!(p.gt, GT_UNSURE);
    }
}
