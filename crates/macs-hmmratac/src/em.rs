//! `MACS3.Signal.HMMR_EM.HMMR_EM` -- the fragment-length EM that supplies the
//! mono/di/tri means and standard deviations.
//!
//! A port of `HMMR_EM.py:87`. The whole class is deterministic given the
//! downsample, and the downsample is seeded from `--randomSeed`, so this is
//! reproducible -- unlike the HMM Baum-Welch fit that follows it.
//!
//! # The EM fits three modes, not four
//!
//! The **short** mode's mean and standard deviation are taken from the user's
//! `--means`/`--stddevs` and never touched. EM only refines indices 1..3 (mono,
//! di, tri), which is why `hmmratac_cmd.py` passes `options.em_means[1:4]` and
//! re-prepends index 0 afterwards.
//!
//! # Assignment is strict-max; ties discard the fragment
//!
//! `return_greater` finds the indices equal to the array maximum and returns
//! `-1` if there is more than one. A fragment that ties between two modes is
//! **dropped entirely** -- it is not assigned to the lower-indexed mode. It then
//! contributes to neither `total` nor any `c[j]`, so it also does not affect the
//! weight renormalisation.
//!
//! # The update is a damped online (Welford) update
//!
//! Per mode, `c/mean/s` are updated by Welford's recurrence -- `mean += (x -
//! old_mean)/c`, `s += (x - new_mean) * (x - old_mean)` -- but the *committed*
//! values are damped by `jump`:
//!
//! ```text
//! means_new = means_old + jump * (mean_running - means_old)
//! var_new   = vars_old   + jump * (var_running   - vars_old)
//! ```
//!
//! with `jump` defaulting to `0.5`. `fragStddevs` is recomputed as
//! `sqrt(fragVars)` after each mode's commit, so it is always consistent with the
//! variance rather than damped independently.
//!
//! # `s[0] = 0.0` on the first observation, so `var = s/c = 0.0`
//!
//! `online_update` returns `float(0.0)` for `s` when `c == 1`, so the first
//! fragment of a mode gives that mode zero variance, and the damped commit
//! halves the old variance toward 0.
//!
//! # Weights are renormalised only over assigned fragments
//!
//! `weights[j] = c[j] / total` where `total` counts assigned fragments. A mode
//! with `c[j] == 0` is **skipped entirely** -- neither its mean, variance nor
//! weight is updated -- which leaves its weight at the previous value while the
//! others are renormalised against a smaller `total`. The three weights
//! therefore need not sum to 1 after an iteration.
//!
//! # Convergence needs all three modes to settle at once
//!
//! The counter increments per mode whose mean, weight **and** standard deviation
//! each moved by less than `epsilon` (`0.05`); convergence requires all three.
//! The loop is also hard-capped at `maxIter` (20) iterations.

use macs_core::Coord;

/// `online_update` (`HMMR_EM.py:35`): Welford's recurrence, returning
/// `(count, mean, m2)`.
fn online_update(x: f32, c: i64, m: f32, s: f32) -> (i64, f32, f32) {
    let c = c + 1;
    if c == 1 {
        return (1, x, 0.0);
    }
    let delta = x - m;
    let m = m + delta / c as f32;
    let s = s + delta * (x - m);
    (c, m, s)
}

/// `get_weighted_density` (`HMMR_EM.py:56`): `w * pnorm2(x, m, v)` where the third
/// argument is a **variance**, not a standard deviation.
fn weighted_density(x: f32, m: f32, v: f32, w: f32) -> f32 {
    w * crate::pnorm2(x, m, v)
}

/// `return_greater` (`HMMR_EM.py:76`): the index of the unique maximum, or `-1`
/// when the maximum is tied.
fn return_greater(d: [f32; 3]) -> i32 {
    let amax = d.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let idxs: Vec<usize> = (0..3).filter(|&i| d[i] == amax).collect();
    if idxs.len() > 1 {
        -1
    } else {
        idxs[0] as i32
    }
}

/// The EM result: per-mode means, standard deviations, and the weights.
#[derive(Debug, Clone, PartialEq)]
pub struct EmResult {
    /// Fitted means for `[mono, di, tri]`.
    pub means: [f32; 3],
    /// Fitted standard deviations for `[mono, di, tri]`.
    pub stddevs: [f32; 3],
    /// Fitted mode weights, normalised over assigned fragments.
    pub weights: [f32; 3],
    /// Whether the loop reached convergence before `max_iter`.
    pub converged: bool,
    /// Iterations actually run.
    pub iterations: usize,
}

/// Tunables for [`train_em`], all with upstream's defaults.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EmParams {
    /// Shortest fragment length retained for training.
    pub min_fraglen: Coord,
    /// Longest fragment length retained for training.
    pub max_fraglen: Coord,
    /// Fraction of fragments sampled for training.
    pub sample_percentage: f32,
    /// Convergence threshold on mean, weight and standard deviation.
    pub epsilon: f32,
    /// Iteration cap.
    pub max_iter: usize,
    /// Damping factor for the parameter update.
    pub jump: f32,
    /// Seed for the downsample.
    pub seed: i64,
}

impl Default for EmParams {
    fn default() -> Self {
        Self {
            min_fraglen: 100,
            max_fraglen: 1000,
            sample_percentage: 10.0,
            epsilon: 0.05,
            max_iter: 20,
            jump: 0.5,
            seed: 12345,
        }
    }
}

/// One E-step/M-step sweep over the retained lengths.
fn iterate(
    data: &[Coord],
    means: &mut [f32; 3],
    vars: &mut [f32; 3],
    stddevs: &mut [f32; 3],
    weights: &mut [f32; 3],
    jump: f32,
) -> u64 {
    let mut temp = [0.0f32; 3];
    let mut total: u64 = 0;
    let mut run_mean = [0.0f32; 3];
    let mut run_s = [0.0f32; 3];
    let mut c = [0i64; 3];
    for slot in temp.iter_mut() {
        *slot = 0.0;
    }
    for &len in data {
        let x = len as f32;
        for j in 0..3 {
            temp[j] = weighted_density(x, means[j], vars[j], weights[j]);
        }
        let index = return_greater(temp);
        if index >= 0 {
            let j = index as usize;
            let (nc, nm, ns) = online_update(x, c[j], run_mean[j], run_s[j]);
            c[j] = nc;
            run_mean[j] = nm;
            run_s[j] = ns;
            total += 1;
        }
    }
    for j in 0..3 {
        if c[j] == 0 {
            // no fragment of this size: leave mean, variance, stddev and weight
            // untouched
            continue;
        }
        let variance = run_s[j] / c[j] as f32;
        means[j] += jump * (run_mean[j] - means[j]);
        vars[j] += jump * (variance - vars[j]);
        stddevs[j] = vars[j].sqrt();
        weights[j] = c[j] as f32 / total as f32;
    }
    total
}

/// `HMMR_EM.__init__` + `__learn` + `__iterate`.
///
/// `data` is the already-downsampled, already-windowed list of fragment lengths
/// (see [`EmParams::min_fraglen`]/[`max_fraglen`]); the caller performs the
/// downsample so this function stays free of RNG concerns.
///
/// `init_means`/`init_stddevs` are the **`[mono, di, tri]`** slices, matching
/// what `hmmratac_cmd.py` passes (`options.em_means[1:4]`).
#[allow(clippy::needless_range_loop)]
pub fn train_em(
    data: &[Coord],
    init_means: [f32; 3],
    init_stddevs: [f32; 3],
    p: &EmParams,
) -> EmResult {
    let mut means = init_means;
    let mut stddevs = init_stddevs;
    let mut vars = [
        stddevs[0] * stddevs[0],
        stddevs[1] * stddevs[1],
        stddevs[2] * stddevs[2],
    ];
    let mut converged = false;

    // initial weights: split by the midpoints between consecutive init means
    let cutoff1 = (init_means[1] - init_means[0]) / 2.0 + init_means[0];
    let cutoff2 = (init_means[2] - init_means[1]) / 2.0 + init_means[1];
    let sum3 = data.len() as f32;
    let sum2 = data.iter().filter(|&&d| (d as f32) < cutoff2).count() as f32;
    let sum1 = data.iter().filter(|&&d| (d as f32) < cutoff1).count() as f32;
    let counter1 = sum1;
    let counter2 = sum2 - sum1;
    let counter3 = sum3 - sum2;
    let mut weights = [counter1 / sum3, counter2 / sum3, counter3 / sum3];

    let mut iterations = 0usize;
    while !converged {
        let mut old_means = [0.0f32; 3];
        let mut old_stddevs = [0.0f32; 3];
        let mut old_weights = [0.0f32; 3];
        old_means.copy_from_slice(&means);
        old_stddevs.copy_from_slice(&stddevs);
        old_weights.copy_from_slice(&weights);
        iterate(
            data,
            &mut means,
            &mut vars,
            &mut stddevs,
            &mut weights,
            p.jump,
        );
        iterations += 1;
        let mut counter = 0;
        for i in 0..3 {
            if (old_means[i] - means[i]).abs() < p.epsilon
                && (old_weights[i] - weights[i]).abs() < p.epsilon
                && (old_stddevs[i] - stddevs[i]).abs() < p.epsilon
            {
                counter += 1;
            }
        }
        if counter == 3 {
            converged = true;
        }
        if iterations >= p.max_iter {
            break;
        }
    }

    EmResult {
        means,
        stddevs,
        weights,
        converged,
        iterations,
    }
}

/// The `short/mono/di/tri` vectors `hmmratac_cmd.py` prints and feeds to
/// `generate_weight_mapping`: index 0 is the user's value, 1..4 the EM output.
pub fn prepend_short(short_mean: f32, short_stddev: f32, em: &EmResult) -> ([f32; 4], [f32; 4]) {
    let mut means = [0.0f32; 4];
    means[0] = short_mean;
    means[1..].copy_from_slice(&em.means);
    let mut stddevs = [0.0f32; 4];
    stddevs[0] = short_stddev;
    stddevs[1..].copy_from_slice(&em.stddevs);
    (means, stddevs)
}

/// Round every entry to one decimal, as `hmmratac_cmd.py` does before building
/// the weight mapping.
///
/// Python's `round(x, 1)` is banker's rounding, so `0.25 -> 0.2` and `0.35 ->
/// 0.4`; `f32::round` rounds half away from zero and would give `0.3` and `0.4`.
/// [`macs_stats::py_round`] implements the banker's rule, which is what is
/// wanted here.
pub fn round_to_one_decimal(v: &mut [f32; 4]) {
    for x in v.iter_mut() {
        *x = macs_stats::py_round(f64::from(*x), 1) as f32;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn return_greater_picks_the_unique_max() {
        assert_eq!(return_greater([1.0, 3.0, 2.0]), 1);
        assert_eq!(return_greater([5.0, 1.0, 2.0]), 0);
        assert_eq!(return_greater([1.0, 2.0, 9.0]), 2);
    }

    #[test]
    fn return_greater_discards_a_tie() {
        assert_eq!(return_greater([3.0, 3.0, 1.0]), -1);
        assert_eq!(return_greater([1.0, 1.0, 1.0]), -1);
    }

    #[test]
    fn online_update_first_observation_zeroes_m2() {
        let (c, m, s) = online_update(100.0, 0, 0.0, 0.0);
        assert_eq!((c, m, s), (1, 100.0, 0.0));
    }

    #[test]
    fn online_update_is_welford() {
        let (c, m, s) = online_update(4.0, 0, 0.0, 0.0);
        let (c, m, s) = online_update(2.0, c, m, s);
        assert_eq!(c, 2);
        assert_eq!(m, 3.0);
        assert_eq!(s, 2.0);
    }

    #[test]
    fn em_recovers_a_well_separated_bimodal_mixture() {
        // 200 fragments near 200 (mono) and 200 near 400 (tri)
        let mut data = Vec::new();
        for i in 0..100 {
            data.push(200 + (i % 7) as u64);
            data.push(400 + (i % 11) as u64);
        }
        let p = EmParams {
            max_iter: 60,
            ..Default::default()
        };
        let r = train_em(&data, [200.0, 400.0, 550.0], [20.0, 30.0, 40.0], &p);
        // Each cluster is captured by the mode it was initialised nearest, so
        // the 200 cluster lands in mode 0 and the 400 cluster in mode 1. The
        // third mode is initialised at 550 with nothing near it, keeps
        // `c[2] == 0`, and is therefore never updated at all.
        assert!(
            (r.means[0] - 203.0).abs() < 5.0,
            "mono mean {:?}",
            r.means[0]
        );
        assert!((r.means[1] - 405.0).abs() < 5.0, "di mean {:?}", r.means[1]);
        assert_eq!(r.means[2], 550.0, "empty mode untouched");
        assert_eq!(r.stddevs[2], 40.0, "empty mode stddev untouched");
        assert_eq!(r.weights[2], 0.0);
        // The within-cluster spread is 0..6 and 0..10, so the fitted stddevs must
        // be near that, not near the 20/30 initial values.
        assert!(r.stddevs[0] < 5.0 && r.stddevs[1] < 6.0, "{:?}", r.stddevs);
        assert!(r.converged);
    }

    #[test]
    fn em_is_deterministic() {
        let data: Vec<Coord> = (100..=1000).map(|x| x + (x % 13)).collect();
        let p = EmParams::default();
        let a = train_em(&data, [150.0, 250.0, 400.0], [20.0, 30.0, 40.0], &p);
        let b = train_em(&data, [150.0, 250.0, 400.0], [20.0, 30.0, 40.0], &p);
        assert_eq!(a, b);
    }

    #[test]
    fn em_stops_at_max_iter() {
        let data: Vec<Coord> = (100..=1000).collect();
        let p = EmParams {
            max_iter: 3,
            ..Default::default()
        };
        let r = train_em(&data, [150.0, 250.0, 400.0], [20.0, 30.0, 40.0], &p);
        assert_eq!(r.iterations, 3);
    }

    #[test]
    fn a_mode_with_no_fragments_keeps_its_weight() {
        // every length is far from the 3rd mode, so c[2] stays 0 and its weight
        // must not be renormalised against the other two
        let data: Vec<Coord> = (100..200).collect();
        let p = EmParams::default();
        let r = train_em(&data, [110.0, 150.0, 900.0], [10.0, 10.0, 10.0], &p);
        assert_eq!(r.weights[2], 0.0, "third mode was never assigned");
        assert!(r.weights[0] > 0.0 && r.weights[1] > 0.0);
    }

    #[test]
    fn prepend_short_keeps_the_user_short_value() {
        let em = EmResult {
            means: [1.0, 2.0, 3.0],
            stddevs: [4.0, 5.0, 6.0],
            weights: [0.0; 3],
            converged: true,
            iterations: 1,
        };
        let (m, s) = prepend_short(99.0, 98.0, &em);
        assert_eq!(m, [99.0, 1.0, 2.0, 3.0]);
        assert_eq!(s, [98.0, 4.0, 5.0, 6.0]);
    }

    #[test]
    fn round_to_one_decimal_is_bankers_rounding() {
        // 0.25 -> 0.2 under Python's round-half-to-even; f32::round would give 0.3
        let mut v = [0.25f32, 1.05, 2.35, 3.15];
        round_to_one_decimal(&mut v);
        assert_eq!(v[0], 0.2);
    }
}
