//! The p-value histogram and the cutoff-analysis peak counts.
//!
//! Transcribed from `MACS3.Signal.CallPeakUnit.__pre_computes`, the path taken
//! when `--cutoff-analysis` is requested. It differs from the ordinary path
//! (`__cal_pvalue_qvalue_table`) in two ways that matter for output:
//!
//! * the p-value histogram is built over **run lengths**, not base counts, so the
//!   histogram and the table agree with each other;
//! * it additionally counts, per cutoff, how many peaks would be called and their
//!   total length, which is what `--cutoff-analysis` writes out.
//!
//! # The histogram is keyed by p-score and weighted by run length
//!
//! ```text
//! pre_p = 0
//! for each run i:
//!     this_l = pos[i] - pre_p
//!     pscore_stat[score[i]] += this_l
//!     pre_p = pos[i]
//! ```
//!
//! The first run's length is therefore `pos[0] - 0`, i.e. the run is credited from
//! the contig start. Combined with F5 (run `i` covers `[pos[i-1], pos[i])`) that is
//! consistent, and the lengths sum to `pos[-1]`, not to the contig length — the
//! paired arrays stop where the shorter input stopped (F35), so the contig tail is
//! not counted at all.
//!
//! # `q` is forced monotonically non-increasing, and the tail is zeroed
//!
//! The walk visits p-scores from high to low, clamping `q` to the previous value
//! and stopping at the first `q <= 0`. The loop variable is then reused for a
//! second loop that overwrites everything from the stopping point onward with `0`.
//! See `macs_score::PqTable::from_histogram`, which implements the same walk and
//! the F18 quirk that follows from reusing `i`.
//!
//! # `tmplist` is a fixed 0.3-step ladder, descending
//!
//! `sorted(np.arange(0.3, 10.0, 0.3), reverse=True)`, each value rounded to 5
//! decimals. It is *not* derived from the data, so the q column of
//! `*_cutoff_analysis.txt` has the same 33 rows for every input (or fewer, since
//! rows with `npeaks == 0` are skipped when writing).
//!
//! # The ladder is f64 while the score track is f32
//!
//! `np.arange(0.3, 10.0, 0.3)` produces `numpy.float64`, so every ladder element is
//! an f64, and the comparison `score_array > cutoff` widens the f32 score to f64
//! before comparing. The ladder is therefore *not* rounded to f32: holding
//! `9.9` as `f32` would give `9.899999618530273`, which is below the true `9.9`
//! and would admit a score that upstream rejects. The types here are deliberate.

use macs_core::Coord;
use macs_rle::SignalTrack;
use macs_score::{PScoreHistogram, PqTable};

/// The cutoff ladder `--cutoff-analysis` walks, descending.
///
/// `sorted(np.arange(0.3, 10.0, 0.3), reverse=True)`, each element
/// `round(x, 5)`. The step is floating point, so the values are not exact
/// multiples of 0.3 and the rounding matters: `round(0.6000000000000001, 5)` is
/// `0.6` while the raw value is not.
///
/// Kept as `f64` because that is what upstream holds; see the module docs.
pub fn cutoff_ladder() -> Vec<f64> {
    let mut v: Vec<f64> = Vec::new();
    let mut x = 0.3f64;
    while x < 10.0 {
        v.push(macs_stats::py_round(x, 5));
        x += 0.3;
    }
    v.reverse();
    v
}

/// Per-cutoff peak accounting, accumulated across chromosomes.
#[derive(Debug, Clone, Default)]
pub struct CutoffStats {
    /// `pvalue_length[c]`: total length of regions above cutoff `c`.
    pub pvalue_length: Vec<(f64, u64)>,
    /// `pvalue_npeaks[c]`: number of regions above cutoff `c`.
    pub pvalue_npeaks: Vec<(f64, u64)>,
}

impl CutoffStats {
    /// Add one chromosome's contribution.
    pub fn accumulate(&mut self, other: &CutoffStats) {
        for (c, l) in &other.pvalue_length {
            upsert(&mut self.pvalue_length, *c, *l);
        }
        for (c, n) in &other.pvalue_npeaks {
            upsert(&mut self.pvalue_npeaks, *c, *n);
        }
    }

    /// Total length for a cutoff.
    pub fn length_at(&self, cutoff: f64) -> u64 {
        self.pvalue_length
            .iter()
            .find(|(c, _)| *c == cutoff)
            .map(|(_, l)| *l)
            .unwrap_or(0)
    }

    /// Peak count for a cutoff.
    pub fn npeaks_at(&self, cutoff: f64) -> u64 {
        self.pvalue_npeaks
            .iter()
            .find(|(c, _)| *c == cutoff)
            .map(|(_, n)| *n)
            .unwrap_or(0)
    }

    /// The rows `--cutoff-analysis` writes: p, q, npeaks, lpeaks, avelpeak.
    ///
    /// Cutoffs with no peaks are **skipped**, which is why real files have fewer
    /// than 33 rows.
    pub fn rows(&self, pq: &PqTable) -> Vec<(f32, f32, u64, u64, f64)> {
        let mut out = Vec::new();
        for cutoff in cutoff_ladder() {
            let n = self.npeaks_at(cutoff);
            if n == 0 {
                continue;
            }
            let l = self.length_at(cutoff);
            out.push((
                cutoff as f32,
                pq.qscore_or_zero(cutoff as f32),
                n,
                l,
                l as f64 / n as f64,
            ));
        }
        out
    }
}

fn upsert(v: &mut Vec<(f64, u64)>, cutoff: f64, delta: u64) {
    if let Some(slot) = v.iter_mut().find(|(c, _)| *c == cutoff) {
        slot.1 += delta;
    } else {
        v.push((cutoff, delta));
    }
}

/// The per-cutoff peak count and length for one chromosome.
///
/// Uses the same `tl <= max_gap` merge and the same `>= min_length` acceptance as
/// peak calling, but only counts — it does not build peaks.
pub fn chromosome_cutoff_stats(
    pos: &[Coord],
    score: &[f32],
    max_gap: Coord,
    min_length: Coord,
) -> CutoffStats {
    let mut stats = CutoffStats::default();
    if pos.is_empty() || pos.len() != score.len() {
        return stats;
    }
    let ladder = cutoff_ladder();

    for cutoff in ladder {
        let mut total_l: u64 = 0;
        let mut total_p: u64 = 0;

        // `score_array > cutoff` is strict, and numpy widens the f32 score to f64
        // to compare against the f64 ladder -- so widen here too rather than
        // narrowing the cutoff, which would admit runs upstream rejects.
        let above: Vec<usize> = (0..score.len())
            .filter(|i| f64::from(score[*i]) > cutoff)
            .collect();
        if above.is_empty() {
            upsert(&mut stats.pvalue_length, cutoff, 0);
            upsert(&mut stats.pvalue_npeaks, cutoff, 0);
            continue;
        }

        // end positions are the runs' own ends; start positions are the previous
        // end, with the first forced to 0 (F5 + the upstream `if above_cutoff[0]==0`
        // fix-up)
        let ends: Vec<Coord> = above.iter().map(|i| pos[*i]).collect();
        let starts: Vec<Coord> = above
            .iter()
            .map(|i| if *i == 0 { 0 } else { pos[i - 1] })
            .collect();

        let mut region: Vec<(Coord, Coord)> = vec![(starts[0], ends[0])];
        let mut lastp = ends[0];
        for i in 1..starts.len() {
            let tl = starts[i] - lastp;
            if tl <= max_gap {
                region.push((starts[i], ends[i]));
            } else {
                let len = region[region.len() - 1].1 - region[0].0;
                if len >= min_length {
                    total_l += len;
                    total_p += 1;
                }
                region = vec![(starts[i], ends[i])];
            }
            lastp = ends[i];
        }
        let len = region[region.len() - 1].1 - region[0].0;
        if len >= min_length {
            total_l += len;
            total_p += 1;
        }

        upsert(&mut stats.pvalue_length, cutoff, total_l);
        upsert(&mut stats.pvalue_npeaks, cutoff, total_p);
    }
    stats
}

/// Build the p-value histogram from a chromosome's p-score track.
///
/// Each run contributes its **length**, `pos[i] - pos[i-1]` with `pos[-1]` taken
/// as 0, so the total equals the last breakpoint rather than the contig length.
pub fn accumulate_histogram(hist: &mut PScoreHistogram, pos: &[Coord], score: &[f32]) {
    // F169: `this_l = pos[i] - pre_p` is a plain C subtraction with `pre_p`
    // starting at `0`, so a union that begins **before** the contig contributes a
    // negative length. `saturating_sub` silently clamped it to 0.
    let mut pre_p: i64 = 0;
    for i in 0..pos.len().min(score.len()) {
        let this_l = pos[i] as i64 - pre_p;
        hist.add(score[i], this_l);
        pre_p = pos[i] as i64;
    }
}

/// The full `--cutoff-analysis` pre-computation for a set of chromosomes.
///
/// Returns the p-value histogram and the per-cutoff peak accounting. The q table
/// itself comes from [`macs_score::PqTable::from_histogram`], which performs the
/// same AFDR walk.
pub fn pre_computes(
    chromosomes: &[(&[Coord], &[f32])],
    max_gap: Coord,
    min_length: Coord,
) -> (PScoreHistogram, CutoffStats) {
    let mut hist = PScoreHistogram::new();
    let mut stats = CutoffStats::default();
    for (pos, score) in chromosomes {
        accumulate_histogram(&mut hist, pos, score);
        stats.accumulate(&chromosome_cutoff_stats(pos, score, max_gap, min_length));
    }
    (hist, stats)
}

/// A convenience wrapper for callers holding run-length tracks.
pub fn pre_computes_tracks(
    tracks: &[(&SignalTrack<f32>, &SignalTrack<f32>)],
    max_gap: Coord,
    min_length: Coord,
) -> (PScoreHistogram, CutoffStats) {
    // `SignalTrack` exposes run ends and values; there is no `values()` accessor
    // that yields a contiguous slice, so collect them.
    let mut owned: Vec<(Vec<Coord>, Vec<f32>)> = Vec::with_capacity(tracks.len());
    for (pos, score) in tracks {
        let ends: Vec<Coord> = pos.runs().iter().map(|r| r.end).collect();
        let score_values: Vec<f32> = score.runs().iter().map(|r| r.value).collect();
        owned.push((ends, score_values));
    }
    let borrowed: Vec<(&[Coord], &[f32])> = owned
        .iter()
        .map(|(p, s)| (p.as_slice(), s.as_slice()))
        .collect();
    pre_computes(&borrowed, max_gap, min_length)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_ladder_is_descending_and_starts_just_above_a_third() {
        let l = cutoff_ladder();
        assert!(l.len() >= 32, "expected ~33 entries, got {}", l.len());
        assert!((l[0] - 9.9).abs() < 1e-12, "starts at 9.9, got {}", l[0]);
        assert!(l[0] > l[1], "descending");
        assert!((l[l.len() - 1] - 0.3).abs() < 1e-12, "ends at 0.3");
        for w in l.windows(2) {
            assert!(w[0] > w[1], "strictly descending: {w:?}");
        }
    }

    #[test]
    fn the_ladder_matches_upstreams_arange() {
        // np.arange(0.3, 10.0, 0.3) then round(_,5), reversed
        let expect: Vec<f64> = {
            let mut v = Vec::new();
            let mut x = 0.3f64;
            while x < 10.0 {
                v.push(round5(x));
                x += 0.3;
            }
            v.reverse();
            v
        };
        let got: Vec<f64> = cutoff_ladder();
        assert_eq!(got.len(), expect.len());
        for (g, e) in got.iter().zip(expect.iter()) {
            assert_eq!(*g, *e, "the ladder must be bit-identical to arange+round");
        }
    }

    fn round5(v: f64) -> f64 {
        let s = format!("{v:.5}");
        s.parse().unwrap_or(v)
    }

    #[test]
    fn the_histogram_weights_by_run_length() {
        let mut h = PScoreHistogram::new();
        // three runs ending at 10, 20, 30 => lengths 10, 10, 10
        accumulate_histogram(&mut h, &[10, 20, 30], &[1.0, 1.0, 2.0]);
        assert_eq!(h.total(), 30, "sums to the last breakpoint");
        assert_eq!(h.count(1.0), 20, "two runs of length 10 at p=1");
        assert_eq!(h.count(2.0), 10);
    }

    #[test]
    fn the_first_run_is_credited_from_zero() {
        let mut h = PScoreHistogram::new();
        accumulate_histogram(&mut h, &[50], &[3.0]);
        assert_eq!(h.count(3.0), 50, "the first run is pos[0] - 0");
    }

    #[test]
    fn mismatched_lengths_are_truncated_not_panicked() {
        let mut h = PScoreHistogram::new();
        accumulate_histogram(&mut h, &[10, 20, 30], &[1.0]);
        assert_eq!(h.total(), 10);
    }

    #[test]
    fn an_f32_score_of_zero_point_three_clears_the_f64_cutoff_of_zero_point_three() {
        // The comparison is `score_array > cutoff` where `score_array` is f4 and
        // `cutoff` is a numpy.float64 from np.arange. numpy promotes the *array*
        // to f64 for an f64 scalar, and f64(float32(0.3)) is
        // 0.30000001192092896, which IS greater than 0.3.
        //
        // This is worth pinning because the opposite reading is very natural:
        // a *Python* float scalar is weak and would be cast down to f32, giving
        // equality and `False`. Upstream holds numpy.float64, so the run is
        // counted. Narrowing the cutoff to f32 would silently drop it.
        let pos = vec![0u64, 1000, 2000];
        let score = vec![0.3f32, 0.3, 0.3];
        assert!(
            f64::from(0.3f32) > 0.3,
            "f32(0.3) widened must exceed the f64 cutoff"
        );
        let s = chromosome_cutoff_stats(&pos, &score, 50, 200);
        assert_eq!(
            s.npeaks_at(0.3),
            1,
            "the run clears the cutoff under upstream's promotion rules"
        );
    }

    #[test]
    fn a_score_exactly_equal_in_f64_does_not_clear_the_cutoff() {
        // the comparison really is strict once both sides are the same type
        let pos = vec![0u64, 1000, 2000];
        let score = vec![9.9f32, 9.9, 9.9]; // f32(9.9) widens to 9.899999618...
        let s = chromosome_cutoff_stats(&pos, &score, 50, 200);
        assert_eq!(
            s.npeaks_at(9.9),
            0,
            "f32(9.9) is below the f64 cutoff 9.9, so strict > is false"
        );
    }

    #[test]
    fn a_single_above_cutoff_region_is_counted_once() {
        // one run well above 0.3, well over min_length
        let pos = vec![0u64, 500, 1000, 1500, 2000, 2500];
        let score = vec![0.0f32, 0.0, 5.0, 0.0, 0.0, 0.0];
        let s = chromosome_cutoff_stats(&pos, &score, 50, 200);
        assert_eq!(s.npeaks_at(0.3), 1);
        assert_eq!(
            s.length_at(0.3),
            500,
            "run length, and 500 >= min_length 200"
        );
    }

    #[test]
    fn regions_below_min_length_are_not_counted() {
        let pos = vec![0u64, 100];
        let score = vec![0.0f32, 5.0];
        let s = chromosome_cutoff_stats(&pos, &score, 50, 200);
        assert_eq!(s.npeaks_at(0.3), 0, "100 < 200");
        assert_eq!(s.length_at(0.3), 0);
    }

    #[test]
    fn nearby_runs_merge_and_far_runs_split() {
        // two runs 40bp apart (merge, max_gap 50) and one 400bp away (split)
        let pos = vec![0u64, 200, 400, 1200, 1400];
        let score = vec![0.0f32, 5.0, 5.0, 0.0, 5.0];
        let s = chromosome_cutoff_stats(&pos, &score, 50, 200);
        assert_eq!(
            s.npeaks_at(0.3),
            2,
            "the merged pair is one region, then another"
        );
    }

    #[test]
    fn accumulate_sums_across_chromosomes() {
        let mut a = CutoffStats::default();
        upsert(&mut a.pvalue_length, 0.3, 100);
        upsert(&mut a.pvalue_npeaks, 0.3, 1);
        let mut b = CutoffStats::default();
        upsert(&mut b.pvalue_length, 0.3, 250);
        upsert(&mut b.pvalue_npeaks, 0.3, 2);
        a.accumulate(&b);
        assert_eq!(a.length_at(0.3), 350);
        assert_eq!(a.npeaks_at(0.3), 3);
    }

    #[test]
    fn an_empty_chromosome_contributes_nothing() {
        let s = chromosome_cutoff_stats(&[], &[], 50, 200);
        assert!(s.pvalue_npeaks.is_empty());
    }

    #[test]
    fn rows_skip_cutoffs_with_no_peaks() {
        let mut s = CutoffStats::default();
        upsert(&mut s.pvalue_npeaks, 0.3, 0); // explicitly zero
        upsert(&mut s.pvalue_npeaks, 0.6, 4);
        upsert(&mut s.pvalue_length, 0.6, 1000);
        let pq = PqTable::empty();
        let rows = s.rows(&pq);
        assert_eq!(rows.len(), 1, "the zero-peak cutoff is not written");
        let (c, _q, n, l, a) = rows[0];
        assert!((c - 0.6f32).abs() < 1e-6, "cutoff {c}");
        assert_eq!((n, l), (4, 1000));
        assert_eq!(a, 250.0);
    }

    #[test]
    fn pre_computes_aggregates_histogram_and_stats() {
        let a_pos = [0u64, 500, 1000];
        let a_score = [0.0f32, 5.0, 0.0];
        let b_pos = [0u64, 400, 800];
        let b_score = [0.0f32, 6.0, 0.0];
        let (hist, stats) = pre_computes(&[(&a_pos, &a_score), (&b_pos, &b_score)], 50, 200);
        assert_eq!(hist.total(), 1000 + 800);
        assert_eq!(stats.npeaks_at(0.3), 2, "one region per chromosome");
        assert_eq!(stats.length_at(0.3), 500 + 400);
    }
}
