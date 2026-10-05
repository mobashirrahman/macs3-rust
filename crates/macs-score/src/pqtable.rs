//! The AFDR p-to-q table.
//!
//! Port of `MACS3/Signal/CallPeakUnit.py:752-835` (`__cal_pvalue_qvalue_table`) --
//! the table `callpeak` builds, since this is the only caller of it. The
//! near-identical `ScoreTrack.make_pq_table` (`ScoreTrack.py:490-544`) differs in
//! one declaration, `k`: see the note on `k` below.
//!
//! # The algorithm, exactly
//!
//! Given the histogram `pvalue_stat` of `-log10 p` -> base pairs:
//!
//! ```python
//! N = sum(pvalue_stat.values())
//! k = 1
//! f = -log10(N)
//! pre_q = 2147483647
//! unique_values = sorted(pvalue_stat.keys(), reverse=True)
//! for i in range(len(unique_values)):
//!     v = unique_values[i]
//!     ln = pvalue_stat[v]
//!     q = v + (log10(k) + f)
//!     if q > pre_q: q = pre_q
//!     if q <= 0:
//!         q = 0
//!         break
//!     pvalue2qvalue[v] = q
//!     pre_q = q
//!     k += ln
//! for j in range(i, len(unique_values)):
//!     pvalue2qvalue[unique_values[j]] = 0
//! ```
//!
//! Three details decide whether this is bit-compatible:
//!
//! * `k` is declared `cython.long` here (`CallPeakUnit.py:770`, and likewise in
//!   `__pre_computes`, `:864`), so the rank accumulates **exactly**; it is *not*
//!   the `cython.float` of `ScoreTrack.make_pq_table` (`ScoreTrack.py:510`). The
//!   difference only shows once the rank passes 2^24, which on a human genome it
//!   does long before the AFDR tail: with `k` in `f32` the spacing above 2^24
//!   exceeds 1, `log10(k)` inherits the quantisation, and every q-score in the
//!   tail moves. Measured on the 5 M-read fixture with `callpeak -p 0.01`: 1126 of
//!   69,159 peaks carried a `-log10(qvalue)` that was wrong in its last few
//!   digits, e.g. `5.77227e-05` upstream against `5.48469e-05`.
//! * `q` is `cython.float`, so `v + (log10(k) + f)` is evaluated in `f64`
//!   (`log10` is a Python function, `v` and `k` are widened) and then
//!   **truncated** to `f32`. Reproduced.
//! * The `break` at `q <= 0` leaves the breaking p-score **out** of the table in
//!   the first loop; the second loop then assigns it `0` along with everything
//!   below.
//! * The second loop is `range(i, len)`, and `i` is `len - 1` when the walk
//!   completes. So the **lowest** p-score in the histogram is always forced to
//!   `0` even when the walk produced a positive q for it — see `F18` in
//!   `docs/upstream-findings.md` and the test
//!   `pq_tables_lowest_bucket_is_forced_to_zero`.

use macs_rle::SignalTrack;
use std::collections::HashMap;

use crate::pscore::PScoreHistogram;
use crate::{canonical_score_key, score_from_key, PScoreCache};

/// The initial `pre_q`, as upstream declares it.
///
/// Upstream writes `2147483647` into a `cython.float`, which is not representable
/// and rounds to `2^31`. It only has to be larger than any real q-score, so the
/// exact value is irrelevant — but it is kept as written.
const PRE_Q_INIT: f32 = 2_147_483_647.0;

/// The p-score to q-score map.
#[derive(Debug, Clone)]
pub struct PqTable {
    /// Canonical `f32` bit pattern of the p-score -> q-score.
    ///
    /// The key is [`canonical_score_key`], not the raw bits, so that `-0.0` and
    /// `+0.0` are one entry. Upstream stores its `pvalue2qvalue` in a Python
    /// dict keyed by float, where `-0.0 == 0.0` and `hash(-0.0) == hash(0.0)`,
    /// so a lookup with either sign finds the same q. Keying on raw bits would
    /// make a `+0.0` p-score miss a table built from a `-0.0` one.
    map: HashMap<u32, f32>,
    /// Number of base pairs that fed the table.
    total: u64,
    /// The p-score at which the AFDR walk first produced `q <= 0`, if any.
    cut_at: Option<f32>,
}

impl Default for PqTable {
    /// An empty table: every p-score maps to q-score 0.
    fn default() -> Self {
        Self::empty()
    }
}

impl PqTable {
    /// An empty table: every p-score maps to q-score 0.
    ///
    /// Useful for exercising the peak-calling geometry, where the q-score is not
    /// under test, and as the neutral element when a caller has no p-value track.
    pub fn empty() -> Self {
        PqTable {
            map: HashMap::new(),
            total: 0,
            cut_at: None,
        }
    }

    /// Build the table from a p-score histogram.
    ///
    /// The AFDR rank `k` accumulates in `i64`, as `CallPeakUnit.py:770` declares
    /// it -- see the module doc for why that is not the `f32` of
    /// `ScoreTrack.make_pq_table`.
    ///
    /// An empty histogram yields an empty table (upstream would raise a
    /// `NameError` here, because `i` is never bound; we return an empty table
    /// and let the caller decide, which is strictly safer).
    pub fn from_histogram(h: &PScoreHistogram) -> Self {
        let mut table = PqTable {
            map: HashMap::with_capacity(h.len()),
            total: h.total().max(0) as u64,
            cut_at: None,
        };
        if h.is_empty() {
            return table;
        }

        let n = h.total();
        // upstream: `f: cython.float = -log10(N)`. `log10` is a Python function so
        // the *computation* is f64, but the declared type truncates the result to
        // float32, and that truncation is worth up to a ULP in every q below.
        // Caught by the golden pq tables.
        let f: f32 = (-(n as f64).log10()) as f32;
        let mut k: i64 = 1;
        let mut pre_q: f32 = PRE_Q_INIT;
        let mut cut_at: Option<f32> = None;
        let ordered = h.sorted_descending();
        let mut last_index = 0usize;

        for (i, (v, ln)) in ordered.iter().copied().enumerate() {
            last_index = i;
            // upstream: q = v + (log10(k) + f); `v` and `f` are C floats and `k` a C
            // long, all widened to f64 for the addition (C promotes, and `log10` is a
            // Python function), then the result is truncated to float32 on
            // assignment to `q: cython.float`
            let q = (f64::from(v) + ((k as f64).log10() + f64::from(f))) as f32;
            let q = if q > pre_q { pre_q } else { q };
            if q <= 0.0 {
                cut_at = Some(v);
                break;
            }
            table.map.insert(canonical_score_key(v), q);
            pre_q = q;
            k += ln;
        }

        // F18: upstream's second loop runs from the loop variable `i`, which is
        // `len - 1` when the walk completes rather than only after a break. So
        // the lowest p-score in the histogram is *always* overwritten with 0,
        // discarding whatever the walk computed for it. Reproduced.
        for (v, _n) in ordered.iter().skip(last_index) {
            table.map.insert(canonical_score_key(*v), 0.0);
        }

        table.cut_at = cut_at;
        table
    }

    /// The q-score for a p-score.
    ///
    /// Returns `None` when the p-score is not in the table, which happens for
    /// every score at or below the AFDR cut. Upstream raises `KeyError` there;
    /// callers treat it as `0.0`, which is what the second loop of
    /// `make_pq_table` writes into the table anyway.
    pub fn get(&self, pscore: f32) -> Option<f32> {
        self.map.get(&canonical_score_key(pscore)).copied()
    }

    /// The q-score for a p-score, with `0.0` for anything not in the table.
    pub fn qscore_or_zero(&self, pscore: f32) -> f32 {
        self.map
            .get(&canonical_score_key(pscore))
            .copied()
            .unwrap_or(0.0)
    }

    /// Number of base pairs that fed the table.
    pub fn total(&self) -> u64 {
        self.total
    }

    /// Number of distinct p-scores in the table.
    ///
    /// Smaller than the histogram's `len` when the AFDR walk stopped early.
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// True when the table maps nothing.
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// The p-score at which the AFDR walk first produced a non-positive
    /// q-score, if it stopped early.
    pub fn cut_at(&self) -> Option<f32> {
        self.cut_at
    }

    /// All `(p, q)` pairs, sorted by descending p-score.
    pub fn entries(&self) -> Vec<(f32, f32)> {
        let mut out: Vec<(f32, f32)> = self
            .map
            .iter()
            .map(|(&p, &q)| (score_from_key(p), q))
            .collect();
        out.sort_unstable_by(|a, b| b.0.partial_cmp(&a.0).expect("NaN p-score"));
        out
    }
}

/// Convert a p-score track into a q-score track.
///
/// Bases whose p-score is not in `table` get `0.0`, which is what upstream's
/// table stores for them.
///
/// `cache` is threaded through so the p-scores are not recomputed for every
/// base of a large track; it does not affect the result.
pub fn qscore_track(
    chrom: macs_core::ChromId,
    pscore: &SignalTrack<f32>,
    table: &PqTable,
    _cache: Option<&mut PScoreCache>,
) -> SignalTrack<f32> {
    let mut out = SignalTrack::empty(chrom, pscore.start(), pscore.end());
    for span in pscore.iter() {
        out.push(span.end, table.qscore_or_zero(*span.value));
    }
    out
}

/// Sink that accumulates p-scores across chromosomes into one histogram, and
/// then produces the q-scores.
///
/// The chromosome iteration order matters: upstream's `compute_pvalue` walks
/// `sorted(self.data.keys())`, and the histogram accumulation is not
/// associative, so the order is part of the contract. [`QScoreSink::build`]
/// takes the tracks already in ascending chromosome-name order.
#[derive(Debug, Default)]
pub struct QScoreSink {
    tracks: Vec<SignalTrack<f32>>,
    hist: PScoreHistogram,
}

impl QScoreSink {
    /// An empty sink.
    pub fn new() -> Self {
        QScoreSink::default()
    }

    /// Accumulate one chromosome's p-score track.
    ///
    /// Must be called in ascending chromosome-name order to match upstream.
    pub fn push(&mut self, t: SignalTrack<f32>) {
        self.hist.add_track(&t);
        self.tracks.push(t);
    }

    /// [`Self::push`] with the first span's length measured from `origin`.
    ///
    /// See [`PScoreHistogram::add_track_from`] -- upstream measures it from
    /// coordinate zero, which a coordinate-shifted run no longer starts at.
    pub fn push_from(&mut self, t: SignalTrack<f32>, origin: i64) {
        self.hist.add_track_from(&t, origin);
        self.tracks.push(t);
    }

    /// Merge a histogram fragment computed elsewhere.
    ///
    /// This is what lets the calling pipeline be *streaming*: it can reduce each
    /// chromosome's p-score track to a sparse score -> length map, drop the track, and
    /// still accumulate the genome-wide AFDR histogram. Base-pair counts are `i64`, so
    /// the merge is exact and order-independent; the caller still folds in chromosome
    /// order to keep the argument simple rather than necessary.
    ///
    /// Unlike [`Self::push`] this records no track, so [`Self::build_with_tracks`]
    /// cannot produce per-chromosome q-score tracks from it -- the caller maps those
    /// itself, per chromosome, once the table exists.
    pub fn push_histogram(&mut self, h: &PScoreHistogram) {
        self.hist.merge(h);
    }

    /// The accumulated histogram.
    pub fn histogram(&self) -> &PScoreHistogram {
        &self.hist
    }

    /// Mutable access to the accumulated histogram.
    ///
    /// `--cutoff-analysis` seeds its ladder into the histogram *after* every
    /// chromosome has been folded in and *before* the q-table is built, which is not
    /// expressible as a `push`.
    pub fn histogram_mut(&mut self) -> &mut PScoreHistogram {
        &mut self.hist
    }

    /// The accumulated histogram as `(score, length)` pairs, ascending by score.
    ///
    /// This is the shape `oracle/dump_stages.py` records as `qvalue_table`, and the
    /// release definition's 1e-6 q-score bound is a statement about this structure --
    /// so the stage comparator needs the raw pairs, not just a q-score lookup.
    /// `PqTable` keeps only the q-value map, so this lives here, where the histogram
    /// still does.
    pub fn histogram_pairs(&self) -> Vec<(f32, i64)> {
        self.hist.pairs()
    }

    /// Number of chromosomes accumulated.
    pub fn len(&self) -> usize {
        self.tracks.len()
    }

    /// True when nothing has been accumulated.
    pub fn is_empty(&self) -> bool {
        self.tracks.is_empty()
    }

    /// Build the p-to-q table and convert every accumulated track.
    pub fn build(&self) -> (PqTable, Vec<SignalTrack<f32>>) {
        self.build_with_tracks(true)
    }

    /// [`Self::build`], optionally **discarding** the per-chromosome q-score
    /// tracks.
    ///
    /// The tracks are the largest single allocation in a callpeak run on a real
    /// genome: `SignalTrack` is a `Vec<Run<f32>>`, i.e. 16 bytes per breakpoint,
    /// so a 5 Mbp genome scored at every position is ~80 MB. They are read in
    /// exactly one place -- the `--broad` level-2 cutoff
    /// (`__close_peak_for_broad_region`, which needs the q track and not the
    /// p track) -- so a narrow-mode run can drop them entirely and still produce
    /// byte-identical output. That is the difference between 2x and parity on the
    /// peak-RSS criterion.
    /// [`Self::build_with_tracks`].
    pub fn build_with_tracks(&self, keep: bool) -> (PqTable, Vec<SignalTrack<f32>>) {
        let table = PqTable::from_histogram(&self.hist);
        let out = if keep {
            self.tracks
                .iter()
                .map(|t| qscore_track(t.chrom(), t, &table, None))
                .collect()
        } else {
            Vec::new()
        };
        (table, out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use macs_core::ChromId;

    const C: ChromId = ChromId(0);

    fn hist(pairs: &[(f32, u64)]) -> PScoreHistogram {
        let mut h = PScoreHistogram::new();
        for &(v, n) in pairs {
            h.add(v, n as i64);
        }
        h
    }

    // ---------- the AFDR walk ----------

    #[test]
    fn q_is_pscore_plus_log10_rank_minus_log10_total() {
        // N = 1000, so f = -3; k starts at 1
        let h = hist(&[(5.0, 900), (2.0, 99), (0.5, 1)]);
        let t = PqTable::from_histogram(&h);
        // v = 5.0, k = 1 => q = 5 + (log10(1) - 3) = 2
        assert!((t.get(5.0).unwrap() - 2.0).abs() < 1e-6, "{:?}", t.get(5.0));
        // k becomes 1 + 900 = 901 => q = 2 + (log10(901) - 3)
        let want = 2.0 + (901f64).log10() - 3.0;
        assert!(
            (t.get(2.0).unwrap() - want as f32).abs() < 1e-5,
            "got {:?} want {want}",
            t.get(2.0)
        );
        // 0.5 is the lowest bucket, so F18 forces it to 0
        assert_eq!(t.get(0.5), Some(0.0));
    }

    #[test]
    fn q_is_clamped_to_be_non_increasing() {
        // a p-score that would produce a larger q than the previous one is
        // capped, which is the monotonicity enforcement
        let h = hist(&[(10.0, 1), (9.99, 1), (0.00001, 10_000_000)]);
        let t = PqTable::from_histogram(&h);
        let entries = t.entries();
        assert!(entries.len() >= 2);
        for w in entries.windows(2) {
            assert!(
                w[1].1 <= w[0].1,
                "q increased while p decreased: {:?} then {:?}",
                w[0],
                w[1]
            );
        }
    }

    #[test]
    fn the_walk_stops_at_the_first_non_positive_q_and_the_rest_are_zero() {
        // N = 2e9, so f = -9.3. A p-score of 20 survives (q = 10.7) and 1 does
        // not (q = -8.3), so the walk breaks at the second bucket.
        let h = hist(&[(20.0, 1), (1.0, 1), (0.00001, 2_000_000_000)]);
        let t = PqTable::from_histogram(&h);
        assert_eq!(t.cut_at(), Some(1.0), "the walk should stop at p = 1");
        // everything from the cut downwards is 0
        assert_eq!(t.get(1.0), Some(0.0));
        assert_eq!(t.qscore_or_zero(0.00001), 0.0);
        // ...and the bucket above the cut keeps its computed q
        let want = 20.0f64 - (2_000_000_002f64).log10();
        assert!(
            (f64::from(t.get(20.0).unwrap()) - want).abs() < 1e-3,
            "got {:?} want {want}",
            t.get(20.0)
        );
        // the table holds the surviving bucket plus the forced-zero tail
        assert_eq!(t.len(), h.len(), "every bucket is present, the tail as 0");
    }

    /// The AFDR rank `k` is a C `long` in `__cal_pvalue_qvalue_table`
    /// (`CallPeakUnit.py:770`), so it accumulates exactly. `f32` is exact below
    /// 2^24 and decisive above it, so this is only observable on a genome-scale
    /// histogram: with the rank rounded to `f32` the middle bucket of this one
    /// would report `0.35982170701026917` rather than the
    /// `0.35982173681259155` asserted here.
    #[test]
    fn the_rank_accumulates_as_a_c_long_not_as_a_float() {
        // a human genome's worth of base pairs, split so that the middle bucket's
        // rank is `1 + 21_000_000` -- past 2^24, where f32 spacing is 2
        let n = 2_900_000_000u64;
        let h = hist(&[(30.0, 21_000_000), (2.5, n - 21_000_000 - 1), (0.5, 1)]);
        let t = PqTable::from_histogram(&h);
        // the top bucket is scored at k = 1, so its q is `v - log10(N)` alone:
        // 0x41a4_4d02, the f32 nearest 20.537601
        assert_eq!(t.get(30.0).unwrap().to_bits(), 0x41a4_4d02);
        // the middle bucket is the only one that depends on a rank above 2^24:
        // 0x3eb8_3a8e, the f32 nearest 0.35982173681259155
        assert_eq!(t.get(2.5).unwrap().to_bits(), 0x3eb8_3a8e);
        // the lowest bucket is forced to 0 by F18
        assert_eq!(t.get(0.5), Some(0.0));
    }

    /// With a single bucket, the walk computes `q = v - log10(N)` and F18 then
    /// throws it away, because the only bucket is also the lowest.
    #[test]
    fn a_single_bucket_is_also_the_lowest_and_so_maps_to_zero() {
        let h = hist(&[(4.0, 1000)]);
        let t = PqTable::from_histogram(&h);
        assert_eq!(t.get(4.0), Some(0.0));
        assert!(t.cut_at().is_none(), "the walk did not break; F18 did this");
    }

    /// With two buckets the higher one keeps its computed q, which is the
    /// non-degenerate case.
    #[test]
    fn a_higher_bucket_keeps_its_computed_q() {
        let h = hist(&[(4.0, 1000), (7.0, 1)]);
        let t = PqTable::from_histogram(&h);
        // N = 1001, f = -3.0004; v = 7, k = 1 => q = 7 - 3.0004 = 3.9996
        let want = 7.0 - (1001f64).log10();
        assert!(
            (t.get(7.0).unwrap() - want as f32).abs() < 1e-5,
            "got {:?} want {want}",
            t.get(7.0)
        );
        assert_eq!(t.get(4.0), Some(0.0));
    }

    /// F18: the lowest p-score in the histogram always maps to 0, because
    /// upstream's `for j in range(i, len)` loop starts at the loop variable, which
    /// is `len - 1` when the AFDR walk completes rather than only after a break.
    #[test]
    fn pq_tables_lowest_bucket_is_forced_to_zero() {
        // N = 1000, f = -3. The walk gives q(5) = 2, then k = 901 and
        // q(2) = 2 + (log10(901) - 3) = 1.9547 -- which upstream then throws away.
        let h = hist(&[(5.0, 900), (2.0, 100)]);
        let t = PqTable::from_histogram(&h);
        assert!(
            (t.get(5.0).unwrap() - 2.0).abs() < 1e-6,
            "q(5) = {:?}",
            t.get(5.0)
        );
        assert_eq!(
            t.get(2.0),
            Some(0.0),
            "the lowest bucket must be forced to 0, not the computed 1.9547"
        );
        assert!(t.cut_at().is_none(), "the walk did not break here");
    }

    /// The force-to-zero applies to exactly one bucket when the walk completes,
    /// and to every bucket from the break point down when it does not.
    #[test]
    fn force_to_zero_covers_one_bucket_without_a_break_and_a_suffix_with_one() {
        // no break: only the lowest bucket is forced
        let h = hist(&[(9.0, 10), (6.0, 10), (3.0, 10)]);
        let t = PqTable::from_histogram(&h);
        assert!(t.get(9.0).unwrap() > 0.0);
        assert!(t.get(6.0).unwrap() > 0.0);
        assert_eq!(t.get(3.0), Some(0.0));

        // a break: everything from the cut down is 0
        let h2 = hist(&[(6.0, 1), (3.0, 1), (1.0, 4_000_000_000)]);
        let t2 = PqTable::from_histogram(&h2);
        assert!(t2.cut_at().is_some());
        assert_eq!(t2.qscore_or_zero(1.0), 0.0);
    }

    #[test]
    fn an_empty_histogram_gives_an_empty_table() {
        let t = PqTable::from_histogram(&PScoreHistogram::new());
        assert!(t.is_empty());
        assert_eq!(t.qscore_or_zero(1.0), 0.0);
    }

    #[test]
    fn a_one_base_pair_genome_maps_to_zero_via_f18() {
        // N = 1 => f = 0, k = 1 => the walk computes q = v, and F18 forces the
        // single (and therefore lowest) bucket to 0
        let h = hist(&[(3.0, 1)]);
        let t = PqTable::from_histogram(&h);
        assert_eq!(t.get(3.0), Some(0.0));
    }

    #[test]
    fn q_never_exceeds_p() {
        // -log10(q) <= -log10(p) always, since q >= p
        let h = hist(&[(6.0, 500), (4.0, 2000), (2.0, 40_000), (0.5, 500_000)]);
        let t = PqTable::from_histogram(&h);
        for (p, q) in t.entries() {
            assert!(q <= p + 1e-6, "q {q} > p {p}");
        }
        // and the lowest bucket is 0, which trivially satisfies it
        assert_eq!(t.get(0.5), Some(0.0));
    }

    // ---------- q-score track ----------

    #[test]
    fn qscore_track_preserves_run_structure() {
        let p = SignalTrack::from_runs(
            C,
            0,
            100,
            vec![macs_rle::Run::new(20, 3.0f32), macs_rle::Run::new(60, 1.0)],
        );
        let h = hist(&[(3.0, 20), (1.0, 40)]);
        let t = PqTable::from_histogram(&h);
        let q = qscore_track(C, &p, &t, None);
        assert_eq!(q.value_at(10), Some(t.qscore_or_zero(3.0)));
        assert_eq!(q.value_at(40), Some(t.qscore_or_zero(1.0)));
    }

    #[test]
    fn qscore_track_maps_cut_scores_to_zero() {
        let p = SignalTrack::from_runs(C, 0, 10, vec![macs_rle::Run::new(10, 0.00001f32)]);
        let h = hist(&[(0.00001, 1_000_000_000)]);
        let t = PqTable::from_histogram(&h);
        assert_eq!(qscore_track(C, &p, &t, None).value_at(5), Some(0.0));
    }

    // ---------- the sink ----------

    #[test]
    fn sink_sums_across_chromosomes() {
        // 150 base pairs at p-score 4 => N = 150, f = -2.176, q = 4 - 2.176 > 0
        let a = SignalTrack::from_runs(ChromId(0), 0, 100, vec![macs_rle::Run::new(100, 4.0f32)]);
        let b = SignalTrack::from_runs(ChromId(1), 0, 50, vec![macs_rle::Run::new(50, 4.0f32)]);
        let mut sink = QScoreSink::new();
        sink.push(a);
        sink.push(b);
        assert_eq!(sink.len(), 2);
        assert_eq!(sink.histogram().count(4.0), 150);
        let (table, tracks) = sink.build();
        assert_eq!(tracks.len(), 2);
        assert!(table.get(4.0).is_some(), "q = 4 - log10(150) = 1.82 > 0");
        // both chromosomes get the same q, because the table is global
        assert_eq!(tracks[0].value_at(10), tracks[1].value_at(10));
    }
}
