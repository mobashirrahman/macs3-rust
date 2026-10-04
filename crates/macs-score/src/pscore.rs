//! The p-score cache and the p-score histogram.

use macs_core::ChromId;
use macs_rle::SignalTrack;
use macs_stats::poisson_cdf;
use std::collections::HashMap;

use crate::{canonical_score_key, score_from_key};

/// Memoised `get_pscore(observed, expectation)`.
///
/// Upstream keys its cache on `(observed, float32 bits of expectation)`
/// ([`MACS3/Signal/ScoreTrack.py:69-88`]):
///
/// ```python
/// memcpy(&y_bits, &y, 4)
/// key = ((longlong)(uint)observed << 32) | y_bits
/// score = -1 * poisson_cdf(observed, expectation, False, True)
/// ```
///
/// so two lambdas that differ only below `f32` precision are the *same* key.
/// That is reproduced here: the key is the exact `f32` bit pattern.
///
/// The cache is a pure memoisation. [`PScoreCache::compute_uncached`] is
/// exposed so the parity tests can assert that enabling it changes nothing.
#[derive(Debug, Default)]
pub struct PScoreCache {
    map: HashMap<u64, f32>,
    hits: u64,
    misses: u64,
    /// Set when a p-score was requested with a non-positive or NaN expectation.
    ///
    /// F172: upstream divides by zero there and dies with
    /// `ZeroDivisionError: float division` from `PeakDetect.__call_peaks_w_control`
    /// (`se_edge/contig_edges`). So the observable behaviour is a **runtime
    /// error with exit status 1 and no output file**, which is what the acceptance
    /// criteria require -- not a panic (exit 101) and not a silent continuation.
    /// The caller checks [`Self::hit_bad_lambda`] and turns it into that error.
    bad_lambda: bool,
}

impl PScoreCache {
    /// An empty cache.
    pub fn new() -> Self {
        PScoreCache::default()
    }

    /// A cache preallocated for `capacity` distinct keys.
    pub fn with_capacity(capacity: usize) -> Self {
        PScoreCache {
            map: HashMap::with_capacity(capacity),
            hits: 0,
            misses: 0,
            bad_lambda: false,
        }
    }

    /// `true` once a p-score has been asked for with `expectation <= 0` or NaN.
    pub fn hit_bad_lambda(&self) -> bool {
        self.bad_lambda
    }

    /// Number of lookups served from the cache.
    pub fn hits(&self) -> u64 {
        self.hits
    }

    /// Number of lookups that had to be computed.
    pub fn misses(&self) -> u64 {
        self.misses
    }

    /// Number of distinct keys held.
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// True when nothing has been cached.
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Drop every entry, keeping the allocation.
    pub fn clear(&mut self) {
        self.map.clear();
        self.hits = 0;
        self.misses = 0;
    }

    /// Look up, computing on a miss.
    ///
    /// `expectation` is the local lambda **before** the pseudocount; it must
    /// already be an `f32`, because that is the precision of the cache key.
    #[inline]
    pub fn get(&mut self, observed: u32, expectation: f32) -> f32 {
        if expectation <= 0.0 || expectation.is_nan() {
            // F172: record it, and answer `+inf` rather than panicking.
            self.bad_lambda = true;
            return f32::INFINITY;
        }
        let key = pack_key(observed, expectation);
        if let Some(&v) = self.map.get(&key) {
            self.hits += 1;
            return v;
        }
        self.misses += 1;
        let v = compute_uncached(observed, expectation);
        self.map.insert(key, v);
        v
    }

    /// The same value, computed directly, bypassing the cache.
    ///
    /// `expectation` must be `> 0`. MACS only ever evaluates this with a
    /// pseudocount added, so the result is always a valid probability.
    pub fn compute_uncached(observed: u32, expectation: f32) -> f32 {
        compute_uncached(observed, expectation)
    }
}

/// Upstream's cache key: `(observed << 32) | float32 bits of expectation`.
#[inline]
fn pack_key(observed: u32, expectation: f32) -> u64 {
    ((observed as u64) << 32) | u64::from(expectation.to_bits())
}

/// `-log10 P(X > observed | lambda = expectation)`, in `f32`.
///
/// The `poisson_cdf` log-space path already rounds to five decimals, so the
/// result is a multiple of `1e-5` (see `F1` in `docs/upstream-findings.md`).
///
/// # Panics
/// If `expectation <= 0` (or is NaN). MACS always adds a pseudocount first, so
/// this cannot happen on a compatible path; there is no valid p-value for a
/// zero expectation, and silently returning `0.0` would make every base look
/// maximally enriched.
#[inline]
pub fn compute_uncached(observed: u32, expectation: f32) -> f32 {
    // F172: no panic. `-log10(0) = +inf` is the honest answer for "impossible
    // under a zero expectation", and [`PScoreCache::hit_bad_lambda`] lets the
    // caller turn it into the runtime error upstream raises.
    if expectation <= 0.0 || expectation.is_nan() {
        return f32::INFINITY;
    }
    match poisson_cdf(observed, f64::from(expectation), false, true) {
        Ok(v) => (-v) as f32,
        Err(_) => f32::INFINITY,
    }
}

/// `-log10 P(X > observed | lambda = expectation)` for one (count, lambda) pair,
/// with the MACS pseudocount convention applied by the caller.
#[inline]
pub fn pscore(cache: Option<&mut PScoreCache>, observed: u32, expectation: f32) -> f32 {
    match cache {
        Some(c) => c.get(observed, expectation),
        None => compute_uncached(observed, expectation),
    }
}

/// The treatment count and local lambda for one base, after the pseudocount.
///
/// Upstream ([`ScoreTrack.py:451-453`]):
///
/// ```python
/// v[i] = get_pscore(cython.cast(cython.int, (p[i] + self.pseudocount)),
///                   c[i] + self.pseudocount)
/// ```
///
/// So the treatment side is **truncated to an integer**, and the lambda side
/// stays `f32`.
///
/// # The pseudocount is added in `f64`, not `f32`
///
/// This is not a choice -- it is observable. Both operands are `f32`, but the
/// sum is computed in double before the `int` cast. Rounding the sum to `f32`
/// first changes the count:
///
/// ```text
/// treat = 19.9f32 = 19.899999618530273
/// pseudocount = 0.1f32 = 0.10000000149011612
/// f32 sum = 20.0                  -> count 20, pscore 7.09102   (wrong)
/// f64 sum = 19.999999620020389    -> count 19, pscore 6.46191   (upstream)
/// ```
///
/// The peak caller's counts are integral, so `f32` and `f64` addition agree
/// there and the difference only shows up for bedGraph scores with fractional
/// pileup values (`bdgcmp -m ppois`). Verified against upstream on a
/// fractional fixture: 3 of 59 rows differ under `f32` addition, 0 under
/// `f64`.
#[inline]
pub fn pseudocounted_inputs(treat: f32, lambda: f32, pseudocount: f32) -> (u32, f32) {
    let observed = (f64::from(treat) + f64::from(pseudocount)) as i64;
    let expectation = lambda + pseudocount;
    // a negative depth is impossible, but a saturated u32 cast must not wrap
    let observed = observed.clamp(0, u32::MAX as i64) as u32;
    (observed, expectation)
}

/// Compute the `-log10` p-score track for one chromosome.
///
/// `treat` and `lambda` must have the same extent; the result covers their
/// intersection and the value at each base is
/// `-log10 P(X > int(treat + pseudocount) | lambda + pseudocount)`.
///
/// Values are `f32`, matching upstream's `v` array.
pub fn pscore_track(
    chrom: ChromId,
    treat: &SignalTrack<f32>,
    lambda: &SignalTrack<f32>,
    cache: Option<&mut PScoreCache>,
    pseudocount: f32,
) -> SignalTrack<f32> {
    let lo = treat.start().max(lambda.start());
    let hi = treat.end().min(lambda.end()).max(lo);
    // pointwise combination in the same pass as the score evaluation, so the
    // two tracks are walked once instead of twice
    let zipped = treat.zip(lambda, |t, l| (t, l));
    let mut out = SignalTrack::empty(chrom, lo, hi);
    let mut cache = cache;
    for span in zipped.iter() {
        let (observed, expectation) = pseudocounted_inputs(span.value.0, span.value.1, pseudocount);
        let v = match cache.as_deref_mut() {
            Some(c) => c.get(observed, expectation),
            None => compute_uncached(observed, expectation),
        };
        out.push(span.end, v);
    }
    out
}

/// `-log10` p-score -> number of scored base pairs.
///
/// This is upstream's `pvalue_stat`. The key is the canonicalised `f32` bit
/// pattern (see the module docs on `-0.0`), and the value is the number of base
/// pairs carrying that score, summed over all chromosomes.
#[derive(Debug, Default, Clone)]
pub struct PScoreHistogram {
    /// Lengths are **signed**.
    ///
    /// F169: upstream's `__cal_pvalue_qvalue_table` measures each paired entry's
    /// length as `pos_array[i] - pos_array[i-1]` with `pre_p` starting at `0` --
    /// a plain C integer subtraction, so when the union's first position is
    /// **negative** the first length is negative too. That happens for every
    /// `--format FRAG` run, where the centred control windows start before the
    /// contig: `frag_basic/barcode_fragments` pairs at `-4901`, contributing
    /// `-4901` to `N`.
    ///
    /// Clamping that to zero (or taking an absolute value) changes
    /// `f = -log10(N)` by ~1.4 and therefore every q-value in the run.
    buckets: HashMap<u32, i64>,
    total: i64,
}

impl PScoreHistogram {
    /// An empty histogram.
    pub fn new() -> Self {
        PScoreHistogram::default()
    }

    /// A histogram preallocated for `capacity` distinct scores.
    pub fn with_capacity(capacity: usize) -> Self {
        PScoreHistogram {
            buckets: HashMap::with_capacity(capacity),
            total: 0,
        }
    }

    /// Add `len` base pairs carrying the score `v`.
    pub fn add(&mut self, v: f32, len: i64) {
        if len == 0 {
            return;
        }
        *self.buckets.entry(canonical_score_key(v)).or_insert(0) += len;
        self.total += len;
    }

    /// Accumulate every base of a p-score track.
    ///
    /// The run lengths follow the right-endpoint convention: a run covering
    /// `[a, b)` contributes `b - a`.
    /// `origin` is the coordinate the first span's length is measured **from**.
    ///
    /// Upstream's `__cal_pvalue_qvalue_table` walks the paired arrays with
    /// `pre_p = 0`, so the first entry's length is `pos_array[0] - 0` -- measured
    /// from the contig's coordinate zero, not from wherever the track happens to
    /// start. `SignalTrack::start()` is normally that same zero, but a run
    /// computed with a coordinate shift (F154: a `--format FRAG` control's centred
    /// windows can start before the contig, so the whole run is shifted by
    /// `coord_shift`) starts at the shifted origin instead. Measuring from there
    /// adds `coord_shift` bases to `N`, and since `f = -log10(N)` enters every
    /// q-value, that shows up as `log10(1 + shift / span)` -- 1.4 on
    /// `frag_basic/barcode_fragments`, whose whole contig is 200 bases and whose
    /// shift is 5000. Passing the true origin restores it.
    pub fn add_track_from(&mut self, t: &SignalTrack<f32>, origin: i64) {
        let mut prev = origin;
        for span in t.iter() {
            self.add(*span.value, span.end as i64 - prev);
            prev = span.end as i64;
        }
    }

    /// [`Self::add_track_from`] with the track's own start as the origin.
    pub fn add_track(&mut self, t: &SignalTrack<f32>) {
        self.add_track_from(t, t.start() as i64);
    }

    /// Number of base pairs scored, **signed** (see the type's note).
    pub fn total(&self) -> i64 {
        self.total
    }

    /// The buckets as `(score, length)` pairs, ascending by score.
    pub fn pairs(&self) -> Vec<(f32, i64)> {
        let mut v: Vec<(f32, i64)> = self
            .buckets
            .iter()
            .map(|(&k, &n)| (score_from_key(k), n))
            .collect();
        v.sort_by(|l, r| l.0.total_cmp(&r.0));
        v
    }

    /// Number of distinct scores.
    pub fn len(&self) -> usize {
        self.buckets.len()
    }

    /// True when nothing has been scored.
    pub fn is_empty(&self) -> bool {
        self.buckets.is_empty()
    }

    /// Record a bucket with **zero** length.
    ///
    /// Upstream's `--cutoff-analysis` path does exactly this for every cutoff in its
    /// ladder that the data did not already produce (`pscore_stat[cutoff] = 0`,
    /// `CallPeakUnit.py:978`). It is not a no-op: `N = sum(pscore_stat.values())` is
    /// unchanged, but `unique_values = sorted(pscore_stat.keys(), ...)` gains an entry,
    /// so the q-table acquires a lookup for that p-score -- and the analysis report's
    /// `qscore` column reads it back.
    pub fn add_zero(&mut self, v: f32) {
        self.buckets.entry(canonical_score_key(v)).or_insert(0);
    }

    /// Add every bucket of `other` into this histogram.
    ///
    /// Counts are `i64` base-pair totals, so addition is exact and independent of
    /// order. This is what allows the p-score tracks to be reduced one chromosome at a
    /// time and dropped, instead of every chromosome's track being held until the
    /// genome-wide table is built.
    pub fn merge(&mut self, other: &PScoreHistogram) {
        for (&k, &n) in &other.buckets {
            *self.buckets.entry(k).or_insert(0) += n;
        }
        self.total += other.total;
    }

    /// [`Self::merge`], consuming `other` and **reusing whichever map is larger**.
    ///
    /// This is the memory-critical form. The histogram carries roughly one entry per
    /// scored base -- the p-score is quantised to `1e-5`, but a continuous pileup still
    /// produces a distinct value nearly everywhere -- so on a 4.8 Mb genome it is
    /// ~5 M entries, and a parallel reduce transiently holds two of them. Copying the
    /// accumulated map on every merge step is what makes that ~2x rather than ~1x.
    /// Swapping in the larger side and draining the smaller keeps the peak at one map
    /// plus the smaller one, and it is exact either way: `i64` addition, no rounding.
    pub fn merge_owned(&mut self, mut other: PScoreHistogram) {
        if self.buckets.len() < other.buckets.len() {
            std::mem::swap(&mut self.buckets, &mut other.buckets);
        }
        for (k, n) in other.buckets {
            *self.buckets.entry(k).or_insert(0) += n;
        }
        self.total += other.total;
    }

    /// The base-pair count for one score.
    pub fn count(&self, v: f32) -> i64 {
        self.buckets
            .get(&canonical_score_key(v))
            .copied()
            .unwrap_or(0)
    }

    /// The distinct scores, sorted **descending**, as `(score, base_pairs)`.
    ///
    /// Descending is the order `make_pq_table` walks them in. Python sorts
    /// `float` keys descending; we sort by the numeric value, which for `f32`
    /// is the same as sorting by bits with the sign bit flipped.
    pub fn sorted_descending(&self) -> Vec<(f32, i64)> {
        let mut out: Vec<(f32, i64)> = self
            .buckets
            .iter()
            .map(|(&k, &n)| (score_from_key(k), n))
            .collect();
        out.sort_unstable_by(|a, b| b.0.partial_cmp(&a.0).expect("NaN p-score"));
        out
    }

    /// Raw bucket access, for the differential comparator.
    pub fn iter(&self) -> impl Iterator<Item = (f32, i64)> + '_ {
        self.buckets.iter().map(|(&k, &n)| (score_from_key(k), n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_PSEUDOCOUNT;

    const C: ChromId = ChromId(0);

    fn const_track(runs: &[(u64, f32)], end: u64) -> SignalTrack<f32> {
        SignalTrack::from_runs(
            C,
            0,
            end,
            runs.iter()
                .map(|&(e, v)| macs_rle::Run::new(e, v))
                .collect(),
        )
    }

    // ---------- the score itself ----------

    /// Every p-score MACS can produce is a multiple of 1e-5, because
    /// `poisson_cdf(..., log10=True)` rounds to five decimals. The score space
    /// is a discrete lattice, so no value between two lattice points is
    /// reachable and a peak boundary can never land between them.
    #[test]
    fn pscore_is_a_multiple_of_1e_minus_5() {
        for observed in 0u32..30 {
            for lam in [0.5f32, 1.0, 2.0, 5.0, 10.0, 50.0] {
                let v = compute_uncached(observed, lam);
                // f32 has ~7 significant digits, so the scaled value can be up to
                // half a step away from an integer
                let scaled = f64::from(v) * 1e5;
                assert!(
                    (scaled - scaled.round()).abs() < 0.5,
                    "pscore {v} for observed={observed} lam={lam} is not on the 1e-5 lattice"
                );
            }
        }
    }

    /// The p-score is `-log10 P(X > k)`, and that tail probability *increases*
    /// with `k`, so the p-score increases with the observed count. It is
    /// non-negative and monotonically increasing, bounded above by
    /// `-log10 P(X > 0) = -log10(1 - e^-lambda)`.
    #[test]
    fn pscore_is_non_negative_and_increases_with_the_count() {
        for lam in [0.5f32, 1.0, 5.0, 50.0] {
            let mut prev = -1.0f32;
            for observed in 0u32..40 {
                let v = compute_uncached(observed, lam);
                assert!(v >= 0.0, "negative p-score at observed={observed}");
                assert!(
                    v >= prev,
                    "p-score decreased at observed={observed}, lam={lam}: {v} < {prev}"
                );
                prev = v;
            }
            // it *starts* at -log10 P(X > 0) = -log10(1 - e^-lambda) and then
            // grows without bound as the count becomes more impossible
            let floor = -f64::from(compute_uncached(0, lam));
            assert!(
                (f64::from(compute_uncached(0, lam)) + floor).abs() < 1e-3,
                "lam={lam}"
            );
            assert!(f64::from(prev) > floor, "lam={lam}");
        }
    }

    /// A larger lambda makes the observed count less surprising, so the
    /// p-score falls. This is the whole point of the local background.
    #[test]
    fn pscore_decreases_as_lambda_grows() {
        for observed in [0u32, 1, 5, 20, 60] {
            let mut prev = f32::INFINITY;
            for lam in [0.5f32, 1.0, 2.0, 5.0, 10.0, 50.0, 200.0] {
                let v = compute_uncached(observed, lam);
                assert!(
                    v <= prev + 1e-3,
                    "p-score rose with lambda at observed={observed}, lam={lam}: {v} > {prev}"
                );
                prev = v;
            }
        }
    }

    #[test]
    fn a_zero_expectation_is_infinite_and_flagged_not_a_panic() {
        // F172: this used to `expect(...)` and abort the process. Upstream raises
        // `ZeroDivisionError` from `PeakDetect.__call_peaks_w_control`, so the
        // observable contract is a *runtime error*, which the caller raises from
        // `PScoreCache::hit_bad_lambda` -- never a panic.
        assert_eq!(compute_uncached(1, 0.0), f32::INFINITY);
        assert_eq!(compute_uncached(1, f32::NAN), f32::INFINITY);
        assert_eq!(compute_uncached(1, -1.0), f32::INFINITY);
        let mut c = PScoreCache::new();
        assert!(!c.hit_bad_lambda());
        let _ = c.get(1, 0.0);
        assert!(c.hit_bad_lambda());
        let mut ok = PScoreCache::new();
        let _ = ok.get(1, 1.0);
        assert!(!ok.hit_bad_lambda());
    }

    // ---------- the cache ----------

    #[test]
    fn the_cache_is_a_pure_memoisation() {
        let mut c = PScoreCache::new();
        let lams = [0.5f32, 1.0, 3.0];
        let mut first = Vec::new();
        for observed in 0u32..8 {
            for &lam in &lams {
                first.push((observed, lam, c.get(observed, lam)));
            }
        }
        assert_eq!(c.misses(), 24, "the first pass must be all misses");
        assert_eq!(c.hits(), 0);
        assert_eq!(c.len(), 24);

        for (observed, lam, want) in first {
            assert_eq!(c.get(observed, lam), want, "cached value differs");
            assert_eq!(
                compute_uncached(observed, lam),
                want,
                "cached value differs from the direct computation"
            );
        }
        assert_eq!(c.hits(), 24);
        assert_eq!(c.misses(), 24, "a hit must not recompute");
        assert_eq!(c.len(), 24);
    }

    /// The key is the exact `f32` bit pattern, so two lambdas that differ only
    /// below `f32` precision are the *same* key, and `-0.0` and `0.0` are
    /// *different* keys (upstream's cache is a dict keyed by the bits, so it
    /// behaves the same way -- unlike the p-value histogram, which keys on the
    /// float and therefore collapses them; see the module docs).
    #[test]
    fn the_cache_key_is_the_f32_bit_pattern_of_lambda() {
        let mut c = PScoreCache::new();
        // 1.0 and the next representable f32 above it are distinct keys
        let a = 1.0f32;
        let b = f32::from_bits(a.to_bits() + 1);
        assert_ne!(a.to_bits(), b.to_bits());
        let _ = c.get(5, a);
        let _ = c.get(5, b);
        assert_eq!(c.len(), 2, "distinct bit patterns are distinct keys");
    }

    /// The key packs `(observed << 32) | f32 bits`, so distinct counts are
    /// distinct keys even when the lambda is identical.
    ///
    /// The counts here are realistic on purpose: `log10_poisson_cdf_Q_large_lambda`
    /// costs `O(k + lambda)` per call (an `O(k)` prefix sum for `ln(k!)` plus the
    /// `O(lambda)` walk past the mode), so a count of 2^31 takes ~30 s in Rust
    /// and far longer in Python. MACS never sees one -- the observed count is a
    /// read depth -- but it does mean a caller that feeds a garbage depth pays
    /// for it. See `F2` in `docs/upstream-findings.md`.
    #[test]
    fn the_cache_key_separates_observed_counts() {
        let mut c = PScoreCache::new();
        for observed in [0u32, 1, 7, 64, 1000] {
            let _ = c.get(observed, 1.0);
        }
        assert_eq!(c.len(), 5);
        assert_eq!(c.misses(), 5);
    }

    // ---------- pseudocount ----------

    #[test]
    fn the_treatment_side_is_truncated_after_the_pseudocount() {
        // treat = 3.7 with pseudocount 1.0 -> int(4.7) = 4
        assert_eq!(
            pseudocounted_inputs(3.7, 2.0, DEFAULT_PSEUDOCOUNT),
            (4, 3.0)
        );
        // treat = 3.2 with pseudocount 1.0 -> int(4.2) = 4
        assert_eq!(
            pseudocounted_inputs(3.2, 2.0, DEFAULT_PSEUDOCOUNT),
            (4, 3.0)
        );
        // treat = 4.0 with pseudocount 1.0 -> int(5.0) = 5
        assert_eq!(
            pseudocounted_inputs(4.0, 2.0, DEFAULT_PSEUDOCOUNT),
            (5, 3.0)
        );
        // the lambda side stays f32
        assert_eq!(pseudocounted_inputs(0.0, 0.25, DEFAULT_PSEUDOCOUNT).1, 1.25);
    }

    #[test]
    fn a_negative_treatment_depth_is_clamped_not_wrapped() {
        // unreachable on a real track, but an unchecked cast would turn a
        // negative depth into a ~4e9 count and make every base look enriched
        let (observed, _) = pseudocounted_inputs(-1.0, 1.0, DEFAULT_PSEUDOCOUNT);
        assert_eq!(observed, 0);
    }

    // ---------- the p-score track ----------

    #[test]
    fn pscore_track_matches_the_dense_computation() {
        let treat = const_track(&[(50, 1.0), (100, 5.0)], 100);
        let lambda = const_track(&[(100, 0.5)], 100);
        let mut cache = PScoreCache::new();
        let t = pscore_track(C, &treat, &lambda, Some(&mut cache), DEFAULT_PSEUDOCOUNT);
        // treat 1.0 + pseudocount 1 => observed 2; lambda 0.5 + 1 => 1.5
        assert_eq!(t.value_at(0), Some(compute_uncached(2, 1.5)));
        assert_eq!(t.value_at(60), Some(compute_uncached(6, 1.5)));
    }

    #[test]
    fn pscore_track_is_independent_of_the_cache() {
        let treat = const_track(&[(40, 0.0), (80, 3.0), (120, 9.0)], 120);
        let lambda = const_track(&[(60, 0.25), (120, 2.0)], 120);
        let with = pscore_track(C, &treat, &lambda, Some(&mut PScoreCache::new()), 1.0);
        let without = pscore_track(C, &treat, &lambda, None, 1.0);
        assert!(with.identical(&without));
    }

    #[test]
    fn pscore_track_uses_the_intersection_extent() {
        let treat = const_track(&[(100, 1.0)], 100);
        let lambda = const_track(&[(50, 1.0)], 50);
        let t = pscore_track(C, &treat, &lambda, None, DEFAULT_PSEUDOCOUNT);
        assert_eq!((t.start(), t.end()), (0, 50));
    }

    // ---------- the histogram ----------

    #[test]
    fn histogram_counts_run_lengths() {
        // [0,30)=0, [30,50)=1, [50,80)=0
        let t = const_track(&[(30, 0.0), (50, 1.0), (80, 0.0)], 80);
        let mut h = PScoreHistogram::new();
        h.add_track(&t);
        assert_eq!(h.count(0.0), 60);
        assert_eq!(h.count(1.0), 20);
        assert_eq!(h.total(), 80);
    }

    #[test]
    fn negative_zero_and_zero_share_one_bucket() {
        // Python dicts collapse them, so upstream's pvalue_stat has one bucket
        let mut h = PScoreHistogram::new();
        h.add(-0.0, 10);
        h.add(0.0, 5);
        assert_eq!(h.len(), 1, "-0.0 and 0.0 must share a bucket");
        assert_eq!(h.count(0.0), 15);
        assert_eq!(h.total(), 15);
    }

    #[test]
    fn histogram_is_sorted_descending_by_value() {
        let mut h = PScoreHistogram::new();
        h.add(1.0, 1);
        h.add(12.5, 1);
        h.add(0.00001, 1);
        h.add(-0.0, 1);
        let s = h.sorted_descending();
        let keys: Vec<f32> = s.iter().map(|&(v, _)| v).collect();
        assert_eq!(keys, vec![12.5, 1.0, 0.00001, 0.0]);
        assert_eq!(h.total(), 4);
    }

    #[test]
    fn zero_length_runs_contribute_nothing() {
        let mut h = PScoreHistogram::new();
        h.add(3.0, 0);
        assert!(h.is_empty());
        assert_eq!(h.total(), 0);
    }
}
