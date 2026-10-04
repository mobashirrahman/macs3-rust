//! The per-chromosome peak caller, tying the stages together.
//!
//! This is `MACS3.Signal.CallPeakUnit.call_peaks` and
//! `__chrom_call_peak_using_certain_criteria` assembled into one pass:
//!
//! ```text
//! for each chromosome:
//!     paired  = pair_treat_ctrl(treat_pileup, ctrl_pileup)
//!     scores  = [score_fn(paired.treat, paired.ctrl) for each requested kind]
//!     above   = indices where every score clears its cutoff
//!     regions = merge chunks separated by > max_gap
//!     for each region: close it (with or without sub-peaks)
//! ```
//!
//! # Scores are computed on the paired arrays, before segmentation
//!
//! `__cal_pscore` / `__cal_qscore` are element-wise over the paired arrays, so a
//! run's score is a function of the pileup at its **end** position. Chunks are
//! then recovered as `(pos[i-1], pos[i])` (F5), which is why the score attached to
//! a chunk is the score at its end.
//!
//! # A p-score is `get_pscore(int(treat), ctrl)` — the cast is a truncation
//!
//! Upstream writes `cython.cast(cython.int, a1_ptr[0])`, so a treatment depth of
//! `7.9` becomes observed count `7`, not `8`. `macs_score::pseudocounted_inputs`
//! does the same thing for the pseudocount path; here the pseudocount has already
//! been folded into `ctrl` upstream, so no additional add is applied.
//!
//! # `--call-summits` emits several peaks with identical bounds
//!
//! One peak region yields one entry per sub-peak summit, all sharing
//! `(start, end)`. That is what upstream writes to `*_summits.bed`, and the XLS
//! records each as its own peak, so the *count* of peaks differs between `--broad`
//! / narrow / `--call-summits` modes while the boundaries do not.

use crate::driver::{pair_treat_ctrl, PairedSignal};
use crate::regions::{segment_regions, CallParams, Chunk, Peak, Reject, ScoreKind};

use macs_core::Coord;
use macs_score::{PScoreCache, PqTable};

/// Which score a caller thresholds on, and the cutoff for it.
#[derive(Debug, Clone)]
pub struct Criterion {
    /// The score kind.
    pub kind: ScoreKind,
    /// `-log10` cutoff, or a fold-enrichment / subtraction threshold.
    pub cutoff: f32,
}

// The score kinds live in `regions::ScoreKind`, shared with the close functions
// that re-check the cutoff. Declaring a second enum here would let a `Criterion`
// name one kind and the close logic test another.

/// The p-score track for one chromosome's paired arrays.
///
/// `observed` is `int(treat)` — a truncating cast, matching
/// `cython.cast(cython.int, a1_ptr[0])`.
pub fn cal_pscore(paired: &PairedSignal, cache: &mut PScoreCache) -> Vec<f32> {
    paired
        .treat
        .iter()
        .zip(paired.ctrl.iter())
        .map(|(&t, &c)| {
            // a negative or saturated depth must not wrap; a pileup is non-negative
            let observed = if t < 0.0 {
                0u32
            } else {
                (t as i64).min(u32::MAX as i64) as u32
            };
            cache.get(observed, c)
        })
        .collect()
}

/// The q-score track for one chromosome's paired arrays.
///
/// Upstream indexes `self.pqtable[pscore]` directly, which raises `KeyError` for a
/// p-score below the AFDR cut. The table's tail is already zeroed by
/// [`macs_score::PqTable::from_histogram`], so a miss and a stored zero are the
/// same answer, and a miss is returned as `0.0`.
pub fn cal_qscore(paired: &PairedSignal, pq: &PqTable, cache: &mut PScoreCache) -> Vec<f32> {
    cal_pscore(paired, cache)
        .iter()
        .map(|p| pq.qscore_or_zero(*p))
        .collect()
}

/// Fold enrichment, `(treat + pseudocount) / (ctrl + pseudocount)`.
pub fn cal_fold(paired: &PairedSignal, pseudocount: f32) -> Vec<f32> {
    paired
        .treat
        .iter()
        .zip(paired.ctrl.iter())
        .map(|(&t, &c)| (t + pseudocount) / (c + pseudocount))
        .collect()
}

/// `treat - ctrl`.
pub fn cal_subtraction(paired: &PairedSignal) -> Vec<f32> {
    paired
        .treat
        .iter()
        .zip(paired.ctrl.iter())
        .map(|(&t, &c)| t - c)
        .collect()
}

/// The mutable state a calling pass needs.
///
/// Bundled because it is threaded through every stage; upstream keeps the same
/// things as attributes of `CallerFromAlignments`.
#[derive(Debug)]
pub struct CallerState<'a> {
    /// The pre-computed p -> q table.
    pub pq: &'a PqTable,
    /// The p-score memo, keyed by `(observed << 32) | f32 bits of expectation`.
    pub cache: &'a mut PScoreCache,
    /// Fold-enrichment pseudocount.
    pub pseudocount: f32,
}

/// Compute one score track.
pub fn cal_score(kind: ScoreKind, paired: &PairedSignal, state: &mut CallerState<'_>) -> Vec<f32> {
    let pq = state.pq;
    let cache = &mut *state.cache;
    let pseudocount = state.pseudocount;
    match kind {
        ScoreKind::P => cal_pscore(paired, cache),
        ScoreKind::Q => cal_qscore(paired, pq, cache),
        ScoreKind::Fold => cal_fold(paired, pseudocount),
        ScoreKind::Subtraction => cal_subtraction(paired),
    }
}

/// The chunks whose scores clear every criterion.
///
/// Upstream allows a **combination** of criteria, and the combination is an AND:
/// `apply_multiple_cutoffs` keeps an index only if every score clears its own
/// cutoff. An empty criterion list would keep everything, which cannot happen
/// because `call_peaks` asserts the two lists are the same length.
pub fn above_cutoff(
    paired: &PairedSignal,
    scores: &[Vec<f32>],
    criteria: &[Criterion],
) -> Vec<usize> {
    let n = paired.pos.len();
    (0..n)
        .filter(|i| {
            criteria.iter().all(|c| {
                let v = scores
                    .iter()
                    .find(|s| s.len() == n)
                    .map(|s| s[*i])
                    .unwrap_or(0.0);
                // which score array belongs to which criterion is positional
                v > c.cutoff
            })
        })
        .collect()
}

/// Call peaks for one chromosome from its paired treatment/control arrays.
///
/// `treat` and `ctrl` are the already-piled-up `[pos, val]` pairs; see
/// [`pair_treat_ctrl`] for how they are merged into the paired form.
pub fn call_chromosome(
    treat_pos: &[Coord],
    treat_val: &[f32],
    ctrl_pos: &[Coord],
    ctrl_val: &[f32],
    criteria: &[Criterion],
    params: &CallParams,
    state: &mut CallerState<'_>,
) -> (Vec<Peak>, usize) {
    state.pseudocount = params.pseudocount as f32;
    let paired = pair_treat_ctrl(treat_pos, treat_val, ctrl_pos, ctrl_val);
    call_chromosome_paired(&paired, criteria, params, state)
}

/// Call peaks from an already-paired chromosome.
pub fn call_chromosome_paired(
    paired: &PairedSignal,
    criteria: &[Criterion],
    params: &CallParams,
    state: &mut CallerState<'_>,
) -> (Vec<Peak>, usize) {
    let pq = state.pq;
    let n = paired.pos.len();
    if n == 0 || criteria.is_empty() {
        return (Vec::new(), 0);
    }

    // one score array per criterion, in the order given
    let scores: Vec<Vec<f32>> = criteria
        .iter()
        .map(|c| cal_score(c.kind, paired, state))
        .collect();

    // AND across criteria, matching apply_multiple_cutoffs
    let above: Vec<usize> = (0..n)
        .filter(|i| {
            criteria
                .iter()
                .enumerate()
                .all(|(k, c)| scores[k][*i] > c.cutoff)
        })
        .collect();
    if above.is_empty() {
        return (Vec::new(), 0);
    }

    // chunks: end is the run's own position, start is the previous run's, with
    // the first forced to 0 (F5 + the upstream `if above_cutoff[0] == 0` fix-up)
    let chunks: Vec<Chunk> = above
        .iter()
        .map(|i| Chunk {
            start: if *i == 0 { 0 } else { paired.pos[i - 1] },
            end: paired.pos[*i],
            treat: paired.treat[*i],
            ctrl: paired.ctrl[*i],
            score_index: *i,
        })
        .collect();

    if std::env::var("CALLPEAK_CHUNKS").is_ok() {
        eprintln!("CHUNKS {} REGIONS_PENDING", chunks.len());
    }
    let regions = segment_regions(&chunks, params.max_gap);
    if std::env::var("CALLPEAK_CHUNKS").is_ok() {
        eprintln!("REGIONS {}", regions.len());
    }
    let mut peaks = Vec::new();
    let mut rejected = 0usize;
    // the cutoff re-check inside the close functions reads the *combined* score
    // arrays positionally, so hand them over in criterion order
    let combined: Vec<(&[f32], f32)> = criteria
        .iter()
        .enumerate()
        .map(|(k, c)| (scores[k].as_slice(), c.cutoff))
        .collect();

    for region in &regions {
        let r = if params.call_summits {
            crate::regions::close_peak_with_subpeaks(
                region,
                params,
                &combined,
                pq,
                Some(state.cache),
            )
        } else {
            crate::regions::close_peak_wo_subpeaks(region, params, &combined, pq, Some(state.cache))
                .map(|p| vec![p])
        };
        match r {
            Ok(v) => peaks.extend(v),
            Err(Reject::BelowCutoff) => {
                // upstream returns False from close, which discards this region
                rejected += 1;
            }
            Err(_) => rejected += 1,
        }
    }
    (peaks, rejected)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Borrow two caches mutably in one call, which a tuple of references cannot
    /// express at the call site.
    fn state<'a>(pq: &'a PqTable, cache: &'a mut PScoreCache) -> CallerState<'a> {
        CallerState {
            pq,
            cache,
            pseudocount: 1.0,
        }
    }

    fn paired() -> PairedSignal {
        // five runs, a clear hump in the middle
        PairedSignal {
            pos: vec![200, 400, 600, 800, 1000],
            treat: vec![1.0, 9.0, 12.0, 9.0, 1.0],
            ctrl: vec![1.0, 1.0, 1.0, 1.0, 1.0],
        }
    }

    #[test]
    fn a_p_score_track_has_one_value_per_run() {
        let p = paired();
        let mut cache = PScoreCache::new();
        let s = cal_pscore(&p, &mut cache);
        assert_eq!(s.len(), p.pos.len());
        // the deepest run has the largest -log10 p
        assert!(s[2] > s[1], "deeper treatment -> larger pscore: {s:?}");
        assert!(s[1] > s[0]);
    }

    #[test]
    fn the_observed_count_truncates_the_treatment_depth() {
        // 7.9 must be observed count 7, not 8
        let p = PairedSignal {
            pos: vec![10],
            treat: vec![7.9],
            ctrl: vec![1.0],
        };
        let mut a = PScoreCache::new();
        let got = cal_pscore(&p, &mut a);
        let mut b = PScoreCache::new();
        let want = b.get(7, 1.0);
        assert_eq!(got[0], want, "int(7.9) == 7");
    }

    #[test]
    fn a_negative_depth_does_not_wrap() {
        let p = PairedSignal {
            pos: vec![10],
            treat: vec![-3.0],
            ctrl: vec![1.0],
        };
        let mut cache = PScoreCache::new();
        let got = cal_pscore(&p, &mut cache);
        let mut b = PScoreCache::new();
        assert_eq!(got[0], b.get(0, 1.0), "clamped to 0, not 2^32-3");
    }

    #[test]
    fn fold_and_subtraction_are_element_wise() {
        let p = PairedSignal {
            pos: vec![10, 20],
            treat: vec![4.0, 6.0],
            ctrl: vec![2.0, 3.0],
        };
        let f = cal_fold(&p, 1.0);
        assert!((f[0] - 5.0 / 3.0).abs() < 1e-6, "{f:?}");
        assert!((f[1] - 7.0 / 4.0).abs() < 1e-6);
        let s = cal_subtraction(&p);
        assert_eq!(s, vec![2.0, 3.0]);
    }

    #[test]
    fn q_scores_are_the_table_lookup_of_the_p_scores() {
        let p = paired();
        let mut cache = PScoreCache::new();
        let hist = macs_score::PScoreHistogram::new();
        let pq = PqTable::from_histogram(&hist);
        let ps = cal_pscore(&p, &mut cache);
        let qs = cal_qscore(&p, &pq, &mut cache);
        assert_eq!(qs.len(), ps.len());
        for (p_, q_) in ps.iter().zip(qs.iter()) {
            assert_eq!(*q_, pq.qscore_or_zero(*p_));
        }
    }

    #[test]
    fn a_cutoff_below_every_run_makes_the_whole_chromosome_one_peak() {
        // every run's p-score exceeds 0, so all five chunks join one region whose
        // start is forced to 0 by the `above_cutoff[0] == 0` rule (F5)
        let p = paired();
        let params = CallParams {
            // F260: unit-test call sites; the shipped path takes this
            // from `ChromCall::clamp_floor`.
            clamp_floor: 0,
            min_length: 200,
            max_gap: 50,
            call_summits: false,
            pseudocount: 1.0,
        };
        let pq = PqTable::empty();
        let mut cache = PScoreCache::new();
        let crit = vec![Criterion {
            kind: ScoreKind::P,
            cutoff: 0.0,
        }];
        let (peaks, _) = call_chromosome_paired(&p, &crit, &params, &mut state(&pq, &mut cache));
        assert_eq!(peaks.len(), 1);
        assert_eq!(
            peaks[0].start, 1,
            "XLS start is the first chunk start + 1 (F85)"
        );
        assert_eq!(peaks[0].end, 1000, "the last chunk's own end");
    }

    #[test]
    fn a_cutoff_above_the_edges_isolates_the_hump() {
        // p-scores are [0.578, 6.95, 10.20, 6.95, 0.578], so a cutoff of 1 keeps
        // only the middle three runs
        let p = paired();
        let params = CallParams {
            // F260: unit-test call sites; the shipped path takes this
            // from `ChromCall::clamp_floor`.
            clamp_floor: 0,
            min_length: 200,
            max_gap: 50,
            call_summits: false,
            pseudocount: 1.0,
        };
        let pq = PqTable::empty();
        let mut cache = PScoreCache::new();
        let crit = vec![Criterion {
            kind: ScoreKind::P,
            cutoff: 1.0,
        }];
        let (peaks, _) = call_chromosome_paired(&p, &crit, &params, &mut state(&pq, &mut cache));
        assert_eq!(peaks.len(), 1);
        assert_eq!(peaks[0].start, 201, "XLS start is pos[i-1] + 1 (F85)");
        assert_eq!(peaks[0].end, 800, "end is the last kept run's own end");
        // the deepest chunk is index 2, i.e. [400, 600); F62: (600+400+1).div_ceil(2).
        assert_eq!(peaks[0].summit, 501);
    }

    #[test]
    fn a_high_cutoff_rejects_everything() {
        let p = paired();
        let params = CallParams::default();
        let pq = PqTable::empty();
        let mut cache = PScoreCache::new();
        let crit = vec![Criterion {
            kind: ScoreKind::Subtraction,
            cutoff: 1_000.0,
        }];
        let (peaks, rejected) =
            call_chromosome_paired(&p, &crit, &params, &mut state(&pq, &mut cache));
        assert!(peaks.is_empty());
        assert_eq!(rejected, 0, "nothing above cutoff is not a rejection");
    }

    #[test]
    fn criteria_combine_with_and_not_or() {
        let p = paired();
        let params = CallParams::default();
        let pq = PqTable::empty();
        let mut cache = PScoreCache::new();
        // subtraction > 5 AND subtraction > 50: the second is unsatisfiable
        let crit = vec![
            Criterion {
                kind: ScoreKind::Subtraction,
                cutoff: 5.0,
            },
            Criterion {
                kind: ScoreKind::Subtraction,
                cutoff: 50.0,
            },
        ];
        let (peaks, _) = call_chromosome_paired(&p, &crit, &params, &mut state(&pq, &mut cache));
        assert!(peaks.is_empty(), "an unsatisfiable second criterion wins");
    }

    #[test]
    fn call_summits_yields_several_peaks_with_identical_bounds() {
        // Two broad humps inside one above-cutoff region, over a smoothly rising
        // background. Two properties are needed for `--call-summits` to find two
        // summits:
        //
        // * the humps must be wider than the smoothing window, which is
        //   `min_length`; a narrower hump is smoothed flat;
        // * the signal must not be too flat. `too_flat` rejects any region with
        //   fewer than 6 distinct values, so a constant background discards every
        //   maximum -- which is correct upstream behaviour, not a bug.
        //
        // A sawtooth background does *not* work: it makes the smoothed derivative
        // cross zero every few samples, so `maxima` returns dozens of spurious
        // candidates and each is then bounded by a nearby "minimum" only a few
        // bases away, failing the 50 bp width test. Real pileup profiles are
        // smooth for this reason.
        let mut pos = Vec::new();
        let mut treat = Vec::new();
        let mut ctrl = Vec::new();
        for i in 0..600u64 {
            pos.push((i + 1) * 10);
            let x = i as f32;
            let background = 5.0 + x * 0.01;
            let t1 = (x - 175.0).abs();
            let t2 = (x - 455.0).abs();
            let humps = 30.0 * (1.0 - t1 / 75.0).max(0.0) + 25.0 * (1.0 - t2 / 75.0).max(0.0);
            treat.push(background + humps);
            ctrl.push(1.0);
        }
        let p = PairedSignal { pos, treat, ctrl };
        let mut cache = PScoreCache::new();
        let pq = PqTable::empty();
        let crit = vec![Criterion {
            kind: ScoreKind::Subtraction,
            cutoff: 1.0,
        }];
        let narrow = CallParams {
            // F260: unit-test call sites; the shipped path takes this
            // from `ChromCall::clamp_floor`.
            clamp_floor: 0,
            call_summits: false,
            min_length: 200,
            max_gap: 2000,
            pseudocount: 1.0,
        };
        let summits = CallParams {
            // F260: unit-test call sites; the shipped path takes this
            // from `ChromCall::clamp_floor`.
            clamp_floor: 0,
            call_summits: true,
            ..narrow.clone()
        };
        let (one, _) = call_chromosome_paired(&p, &crit, &narrow, &mut state(&pq, &mut cache));
        let (many, _) = call_chromosome_paired(&p, &crit, &summits, &mut state(&pq, &mut cache));
        assert_eq!(one.len(), 1, "narrow mode calls one peak per region");
        assert!(
            many.len() > one.len(),
            "--call-summits should find more summits: {} vs {}",
            many.len(),
            one.len()
        );
        for pk in &many {
            assert_eq!((pk.start, pk.end), (one[0].start, one[0].end));
        }
    }

    #[test]
    fn an_empty_paired_array_calls_nothing() {
        let p = PairedSignal::default();
        let params = CallParams::default();
        let pq = PqTable::empty();
        let mut cache = PScoreCache::new();
        let crit = vec![Criterion {
            kind: ScoreKind::P,
            cutoff: 0.0,
        }];
        let (peaks, rejected) =
            call_chromosome_paired(&p, &crit, &params, &mut state(&pq, &mut cache));
        assert!(peaks.is_empty() && rejected == 0);
    }

    #[test]
    fn call_chromosome_pairs_its_inputs_first() {
        let tp = [0u64, 200, 400, 600, 800];
        let tv = [1.0f32, 9.0, 12.0, 9.0, 1.0];
        let cp = [0u64, 200, 400, 600, 800];
        let cv = [1.0f32, 1.0, 1.0, 1.0, 1.0];
        let params = CallParams {
            // F260: unit-test call sites; the shipped path takes this
            // from `ChromCall::clamp_floor`.
            clamp_floor: 0,
            min_length: 200,
            max_gap: 50,
            ..Default::default()
        };
        let pq = PqTable::empty();
        let mut cache = PScoreCache::new();
        let crit = vec![Criterion {
            kind: ScoreKind::P,
            cutoff: 0.0,
        }];
        let mut st = state(&pq, &mut cache);
        let (peaks, _) = call_chromosome(&tp, &tv, &cp, &cv, &crit, &params, &mut st);
        assert!(!peaks.is_empty());
    }
}
