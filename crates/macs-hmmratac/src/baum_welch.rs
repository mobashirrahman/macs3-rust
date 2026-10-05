//! hmmlearn `GaussianHMM`/`PoissonHMM` Baum-Welch, functionally.
//!
//! `HMMR_HMM.hmm_training` fits with `random_state=RandomState(MT19937(SeedSequence))`
//! and hmmlearn defaults (`n_iter=10`, `tol=1e-2`, full covariance). Initialization
//! is KMeans means ([`crate::kmeans`]) plus global covariance and uniform/random
//! startprob/transmat from the same stream.
//!
//! This implements the EM algorithm correctly (forward-backward E-step, closed-form
//! M-step, log-likelihood convergence) to reach the same optimum hmmlearn reaches
//! from the same init. It does not replicate hmmlearn's C forward-backward
//! summation order or sklearn's exact KMeans draws bit-for-bit; validated by final
//! model agreement and accessible-region identity, not intermediate exactness.

/// Multivariate Gaussian log-density with full covariance.
///
/// `log N(x | mean, cov) = -0.5*(d*ln2pi + ln|det| + (x-m)' inv(cov) (x-m))`.
/// Uses a Cholesky factorization and forward substitution.
pub fn gaussian_logpdf(x: &[f64], mean: &[f64], cov: &[Vec<f64>]) -> f64 {
    let d = x.len();
    // Cholesky: C = L L', stable for positive-definite covariances.
    // log|det| = 2*sum(ln(diag(L))), and the quadratic form solves via forward/back
    // substitution rather than an explicit inverse. Gauss-Jordan was tried first and
    // produced NaN on the yeast global covariance (cond ~2e2) from pivot growth;
    // Cholesky needs no pivoting and degrades gracefully.
    let n = d;
    let mut l = vec![vec![0.0; n]; n];
    let mut det_ln = 0.0;
    for i in 0..n {
        for j in 0..=i {
            let mut s = cov[i][j];
            #[allow(clippy::needless_range_loop)]
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                if s <= 0.0 {
                    return f64::NEG_INFINITY;
                }
                l[i][j] = s.sqrt();
                det_ln += l[i][j].ln();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    det_ln *= 2.0;
    let mut diff = vec![0.0; d];
    for (i, (xi, mi)) in x.iter().zip(mean.iter()).enumerate() {
        diff[i] = xi - mi;
    }
    // Solve L y = diff, then q = y'y.
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut s = diff[i];
        for k in 0..i {
            s -= l[i][k] * y[k];
        }
        y[i] = s / l[i][i];
    }
    let q: f64 = y.iter().map(|v| v * v).sum();
    -0.5 * (d as f64 * (2.0 * std::f64::consts::PI).ln() + det_ln + q)
}

/// NumPy RandomState's legacy Gaussian/Gamma sampler. The polar normal caches
/// its second variate, which is load-bearing for exact gamma draw consumption.
/// Transcribed from NumPy's `legacy-distributions.c` RandomState compatibility
/// path, whose outputs NumPy preserves across Generator changes.
struct LegacySampler<'a> {
    rng: &'a mut macs_stats::NumpyRng,
    cached_gauss: Option<f64>,
}

impl LegacySampler<'_> {
    fn standard_exponential(&mut self) -> f64 {
        -(1.0 - self.rng.random_sample()).ln()
    }

    fn gauss(&mut self) -> f64 {
        if let Some(cached) = self.cached_gauss.take() {
            return cached;
        }
        loop {
            let x1 = 2.0 * self.rng.random_sample() - 1.0;
            let x2 = 2.0 * self.rng.random_sample() - 1.0;
            let r2 = x1 * x1 + x2 * x2;
            if r2 >= 1.0 || r2 == 0.0 {
                continue;
            }
            let f = (-2.0 * r2.ln() / r2).sqrt();
            self.cached_gauss = Some(f * x1);
            return f * x2;
        }
    }

    fn standard_gamma(&mut self, shape: f64) -> f64 {
        if shape == 1.0 {
            return self.standard_exponential();
        }
        if shape == 0.0 {
            return 0.0;
        }
        if shape < 1.0 {
            loop {
                let u = self.rng.random_sample();
                let v = self.standard_exponential();
                if u <= 1.0 - shape {
                    let x = u.powf(1.0 / shape);
                    if x <= v {
                        return x;
                    }
                } else {
                    let y = -((1.0 - u) / shape).ln();
                    let x = (1.0 - shape + shape * y).powf(1.0 / shape);
                    if x <= v + y {
                        return x;
                    }
                }
            }
        }
        let b = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * b).sqrt();
        loop {
            let (x, v0) = loop {
                let x = self.gauss();
                let v0 = 1.0 + c * x;
                if v0 > 0.0 {
                    break (x, v0);
                }
            };
            let v = v0 * v0 * v0;
            let u = self.rng.random_sample();
            if u < 1.0 - 0.0331 * x * x * x * x || u.ln() < 0.5 * x * x + b * (1.0 - v + v.ln()) {
                return b * v;
            }
        }
    }

    fn dirichlet(&mut self, n: usize) -> Vec<f64> {
        // BaseHMM._init passes alpha=1/n_components for every entry.
        let alpha = 1.0 / n as f64;
        let mut values: Vec<f64> = (0..n).map(|_| self.standard_gamma(alpha)).collect();
        let sum: f64 = values.iter().sum();
        for value in &mut values {
            *value /= sum;
        }
        values
    }
}

/// Log-sum-exp over a slice.
fn logsumexp(xs: &[f64]) -> f64 {
    let m = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return m;
    }
    m + xs.iter().map(|x| (x - m).exp()).sum::<f64>().ln()
}

/// Forward-backward in log space for one sequence.
///
/// Returns `(log_alpha, log_beta, log_likelihood)`. `log_em[i][k]` is the log
/// emission density of frame `i` in state `k`.
fn forward_backward(
    log_em: &[Vec<f64>],
    log_start: &[f64],
    log_trans: &[Vec<f64>],
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, f64) {
    let t = log_em.len();
    let k = log_start.len();
    let mut alpha = vec![vec![0.0; k]; t];
    for s in 0..k {
        alpha[0][s] = log_start[s] + log_em[0][s];
    }
    for i in 1..t {
        for s in 0..k {
            let mut v = vec![0.0; k];
            for p in 0..k {
                v[p] = alpha[i - 1][p] + log_trans[p][s];
            }
            alpha[i][s] = logsumexp(&v) + log_em[i][s];
        }
    }
    let mut beta = vec![vec![0.0; k]; t];
    beta[t - 1].fill(0.0);
    for i in (0..t - 1).rev() {
        for s in 0..k {
            let mut v = vec![0.0; k];
            for n in 0..k {
                v[n] = log_trans[s][n] + log_em[i + 1][n] + beta[i + 1][n];
            }
            beta[i][s] = logsumexp(&v);
        }
    }
    let ll = logsumexp(&alpha[t - 1]);
    (alpha, beta, ll)
}

/// Fitted Gaussian HMM parameters.
#[derive(Debug, Clone)]
pub struct GaussianHmm {
    /// Initial state log-probabilities.
    pub log_start: Vec<f64>,
    /// Transition log-probabilities `[from][to]`.
    pub log_trans: Vec<Vec<f64>>,
    /// Emission means `[state][feature]`.
    pub means: Vec<Vec<f64>>,
    /// Emission covariances `[state][i][j]`, full.
    pub covars: Vec<Vec<Vec<f64>>>,
}

/// Baum-Welch for a Gaussian HMM with full covariances.
///
/// `x` is the full feature matrix, `lengths` the sequence break points.
/// `init_means` comes from KMeans; startprob/transmat are initialized from
/// hmmlearn's Dirichlet(1/n_states) priors and covariance starts from the global estimate.
/// Runs at most `n_iter` iterations, stopping when
/// the log-likelihood gain falls below `tol`.
pub fn fit_gaussian(
    x: &[Vec<f64>],
    lengths: &[usize],
    init_means: &[Vec<f64>],
    rng: &mut macs_stats::NumpyRng,
    n_iter: usize,
    tol: f64,
    min_covar: f64,
) -> GaussianHmm {
    let n_states = init_means.len();
    let n_feat = init_means[0].len();
    assert!(n_states > 0 && n_feat > 0 && !x.is_empty());
    // hmmlearn BaseHMM uses Dirichlet(alpha=1/n_states) for starts and rows.
    let mut sampler = LegacySampler {
        rng,
        cached_gauss: None,
    };
    let log_start: Vec<f64> = sampler
        .dirichlet(n_states)
        .into_iter()
        .map(f64::ln)
        .collect();
    let log_trans: Vec<Vec<f64>> = (0..n_states)
        .map(|_| {
            sampler
                .dirichlet(n_states)
                .into_iter()
                .map(f64::ln)
                .collect()
        })
        .collect();
    fit_gaussian_from_init(
        x, lengths, init_means, log_start, log_trans, n_iter, tol, min_covar,
    )
}

/// Train a full-covariance Gaussian HMM from scratch with hmmlearn's default
/// initialization order: Dirichlet start/transitions, then ten sklearn KMeans
/// restarts, then the global covariance estimate.
pub fn train_gaussian(
    x: &[Vec<f64>],
    lengths: &[usize],
    n_states: usize,
    rng: &mut macs_stats::NumpyRng,
    n_iter: usize,
    tol: f64,
    min_covar: f64,
) -> GaussianHmm {
    assert!(n_states > 0 && !x.is_empty() && !x[0].is_empty());
    let (log_start, log_trans) = {
        let mut sampler = LegacySampler {
            rng,
            cached_gauss: None,
        };
        let start = sampler
            .dirichlet(n_states)
            .into_iter()
            .map(f64::ln)
            .collect();
        let trans = (0..n_states)
            .map(|_| {
                sampler
                    .dirichlet(n_states)
                    .into_iter()
                    .map(f64::ln)
                    .collect()
            })
            .collect();
        (start, trans)
    };
    let (means, _) = crate::kmeans::kmeans(x, n_states, 10, 300, 1e-4, rng);
    fit_gaussian_from_init(
        x, lengths, &means, log_start, log_trans, n_iter, tol, min_covar,
    )
}

#[allow(clippy::too_many_arguments)]
fn fit_gaussian_from_init(
    x: &[Vec<f64>],
    lengths: &[usize],
    init_means: &[Vec<f64>],
    mut log_start: Vec<f64>,
    mut log_trans: Vec<Vec<f64>>,
    n_iter: usize,
    tol: f64,
    min_covar: f64,
) -> GaussianHmm {
    let n_states = init_means.len();
    let n_feat = init_means[0].len();
    // Global covariance + min_covar*I, replicated per state.
    let mut global_cov = vec![vec![0.0; n_feat]; n_feat];
    let mut mean_all = vec![0.0; n_feat];
    for row in x.iter() {
        for (i, v) in row.iter().enumerate() {
            mean_all[i] += v;
        }
    }
    for v in mean_all.iter_mut() {
        *v /= x.len() as f64;
    }
    for row in x.iter() {
        for i in 0..n_feat {
            for j in 0..n_feat {
                global_cov[i][j] += (row[i] - mean_all[i]) * (row[j] - mean_all[j]);
            }
        }
    }
    for row in global_cov.iter_mut() {
        for v in row.iter_mut() {
            *v /= (x.len() - 1) as f64;
        }
    }
    for (i, row) in global_cov.iter_mut().enumerate() {
        row[i] += min_covar;
    }
    // `min_covar` is used ONLY here, at init (`hmm.py:316-320`); the M-step adds
    // `covars_prior` instead. Adding `min_covar` every iteration (an earlier version)
    // inflated the variances and changed the EM trajectory.
    let mut covars = vec![global_cov; n_states];
    let mut means = init_means.to_vec();

    // Sequence offsets.
    let offsets = sequence_offsets(x.len(), lengths);

    let mut prev_ll = f64::NEG_INFINITY;
    for _ in 0..n_iter {
        // E-step accumulators.
        let mut start_acc = vec![0.0; n_states];
        let mut trans_acc = vec![vec![0.0; n_states]; n_states];
        let mut w_sum = vec![0.0; n_states];
        let mut wx_sum = vec![vec![0.0; n_feat]; n_states];
        // Raw second moments `sum_t gamma * x x'`. hmmlearn forms the covariance as
        // `obs*obs.T - outer(obs, mu) - outer(mu, obs) + outer(mu, mu)*post`
        // (`hmm.py:382-388`), which is algebraically `sum gamma (x-mu)(x-mu)'` but
        // with a different floating-point order; accumulating the raw moments matches
        // its arithmetic.
        let mut wxx_sum = vec![vec![vec![0.0; n_feat]; n_feat]; n_states];
        let mut total_ll = 0.0;
        for s in 0..lengths.len() {
            let (a, b) = (offsets[s], offsets[s + 1]);
            let seq = &x[a..b];
            let mut log_em = vec![vec![0.0; n_states]; seq.len()];
            for (i, row) in seq.iter().enumerate() {
                for k in 0..n_states {
                    log_em[i][k] = gaussian_logpdf(row, &means[k], &covars[k]);
                }
            }

            let (alpha, beta, ll) = forward_backward(&log_em, &log_start, &log_trans);
            total_ll += ll;
            // Gamma and xi.
            for (i, row) in seq.iter().enumerate() {
                let mut gamma = vec![0.0; n_states];
                for k in 0..n_states {
                    gamma[k] = (alpha[i][k] + beta[i][k] - ll).exp();
                }
                if i == 0 {
                    for k in 0..n_states {
                        start_acc[k] += gamma[k];
                    }
                }
                if i + 1 < seq.len() {
                    for p in 0..n_states {
                        for n in 0..n_states {
                            let xi =
                                (alpha[i][p] + log_trans[p][n] + log_em[i + 1][n] + beta[i + 1][n]
                                    - ll)
                                    .exp();
                            trans_acc[p][n] += xi;
                        }
                    }
                }
                for k in 0..n_states {
                    w_sum[k] += gamma[k];
                    for f in 0..n_feat {
                        wx_sum[k][f] += gamma[k] * row[f];
                    }
                    for fi in 0..n_feat {
                        for fj in 0..n_feat {
                            wxx_sum[k][fi][fj] += gamma[k] * row[fi] * row[fj];
                        }
                    }
                }
            }
        }
        // M-step.
        normalize_logs(&mut log_start, &start_acc);
        for p in 0..n_states {
            normalize_logs(&mut log_trans[p], &trans_acc[p]);
        }
        for k in 0..n_states {
            if w_sum[k] > 0.0 {
                for f in 0..n_feat {
                    means[k][f] = wx_sum[k][f] / w_sum[k];
                }
                for fi in 0..n_feat {
                    for fj in 0..n_feat {
                        // hmmlearn defaults to covars_prior=0.01 and
                        // covars_weight=1. For full covariance the divisor is
                        // post[k] + max(weight - n_features, 0) == post[k].
                        covars[k][fi][fj] = (wxx_sum[k][fi][fj] / w_sum[k]
                            - means[k][fi] * means[k][fj])
                            + 0.01 / w_sum[k];
                    }
                }
            }
        }

        if total_ll - prev_ll < tol {
            break;
        }
        prev_ll = total_ll;
    }
    GaussianHmm {
        log_start,
        log_trans,
        means,
        covars,
    }
}

/// Fitted Poisson HMM parameters.
#[derive(Debug, Clone)]
pub struct PoissonHmm {
    /// Initial state log-probabilities.
    pub log_start: Vec<f64>,
    /// Transition log-probabilities `[from][to]`.
    pub log_trans: Vec<Vec<f64>>,
    /// Poisson rates `[state][feature]`.
    pub lambdas: Vec<Vec<f64>>,
}

/// Fit a Poisson HMM using hmmlearn's default initialization and M-step.
/// `lengths` separates independent sequences within the flattened feature matrix.
pub fn fit_poisson(
    x: &[Vec<f64>],
    lengths: &[usize],
    rng: &mut macs_stats::NumpyRng,
    n_iter: usize,
    tol: f64,
) -> PoissonHmm {
    assert!(!x.is_empty() && !x[0].is_empty());
    let n_states = 3;
    let n_feat = x[0].len();
    let mut sampler = LegacySampler {
        rng,
        cached_gauss: None,
    };
    let mut log_start: Vec<f64> = sampler
        .dirichlet(n_states)
        .into_iter()
        .map(f64::ln)
        .collect();
    let mut log_trans: Vec<Vec<f64>> = (0..n_states)
        .map(|_| {
            sampler
                .dirichlet(n_states)
                .into_iter()
                .map(f64::ln)
                .collect()
        })
        .collect();

    // PoissonHMM._init uses a method-of-moments Gamma draw for every rate,
    // based on the population moments of all entries of X.
    let count = (x.len() * n_feat) as f64;
    let mean = x.iter().flatten().sum::<f64>() / count;
    let variance = x.iter().flatten().map(|v| (v - mean).powi(2)).sum::<f64>() / count;
    let (shape, scale) = if mean > 0.0 && variance > 0.0 {
        (mean * mean / variance, variance / mean)
    } else {
        (1.0, mean.max(f64::MIN_POSITIVE))
    };
    let mut lambdas = vec![vec![0.0; n_feat]; n_states];
    for row in &mut lambdas {
        for rate in row {
            *rate = sampler.standard_gamma(shape) * scale;
        }
    }

    let offsets = sequence_offsets(x.len(), lengths);
    let mut prev_ll = f64::NEG_INFINITY;
    for _ in 0..n_iter {
        let mut start_acc = vec![0.0; n_states];
        let mut trans_acc = vec![vec![0.0; n_states]; n_states];
        let mut post = vec![0.0; n_states];
        let mut obs = vec![vec![0.0; n_feat]; n_states];
        let mut total_ll = 0.0;
        for pair in offsets.windows(2) {
            let seq = &x[pair[0]..pair[1]];
            let log_em: Vec<Vec<f64>> = seq
                .iter()
                .map(|row| {
                    lambdas
                        .iter()
                        .map(|rates| {
                            row.iter()
                                .zip(rates)
                                .map(|(&v, &rate)| {
                                    if rate == 0.0 && v == 0.0 {
                                        0.0
                                    } else if rate <= 0.0 {
                                        f64::NEG_INFINITY
                                    } else {
                                        v * rate.ln() - rate - macs_stats::ln_gamma_ref(v + 1.0)
                                    }
                                })
                                .sum()
                        })
                        .collect()
                })
                .collect();
            let (alpha, beta, ll) = forward_backward(&log_em, &log_start, &log_trans);
            total_ll += ll;
            for t in 0..seq.len() {
                let gamma: Vec<f64> = (0..n_states)
                    .map(|s| (alpha[t][s] + beta[t][s] - ll).exp())
                    .collect();
                if t == 0 {
                    for s in 0..n_states {
                        start_acc[s] += gamma[s];
                    }
                }
                for s in 0..n_states {
                    post[s] += gamma[s];
                    for f in 0..n_feat {
                        obs[s][f] += gamma[s] * seq[t][f];
                    }
                }
                if t + 1 < seq.len() {
                    for i in 0..n_states {
                        for j in 0..n_states {
                            trans_acc[i][j] +=
                                (alpha[t][i] + log_trans[i][j] + log_em[t + 1][j] + beta[t + 1][j]
                                    - ll)
                                    .exp();
                        }
                    }
                }
            }
        }
        normalize_logs(&mut log_start, &start_acc);
        for s in 0..n_states {
            normalize_logs(&mut log_trans[s], &trans_acc[s]);
        }
        for s in 0..n_states {
            if post[s] > 0.0 {
                for f in 0..n_feat {
                    lambdas[s][f] = obs[s][f] / post[s];
                }
            }
        }
        if total_ll - prev_ll < tol {
            break;
        }
        prev_ll = total_ll;
    }
    PoissonHmm {
        log_start,
        log_trans,
        lambdas,
    }
}

fn sequence_offsets(total: usize, lengths: &[usize]) -> Vec<usize> {
    if lengths.is_empty() {
        return vec![0, total];
    }
    let mut offsets = Vec::with_capacity(lengths.len() + 1);
    offsets.push(0);
    for &length in lengths {
        offsets.push(offsets.last().copied().unwrap() + length);
    }
    assert_eq!(
        *offsets.last().unwrap(),
        total,
        "sequence lengths must sum to samples"
    );
    offsets
}

fn normalize_logs(log_values: &mut [f64], counts: &[f64]) {
    let sum: f64 = counts.iter().sum();
    if sum > 0.0 {
        for (log_value, &count) in log_values.iter_mut().zip(counts) {
            *log_value = if count > 0.0 {
                (count / sum).ln()
            } else {
                f64::NEG_INFINITY
            };
        }
    }
}

#[cfg(test)]
mod sampler_tests {
    use super::LegacySampler;

    #[test]
    fn legacy_dirichlet_and_gamma_match_numpy_randomstate() {
        let mut rng = macs_stats::randomstate_from_seed_sequence(12345);
        let mut sampler = LegacySampler {
            rng: &mut rng,
            cached_gauss: None,
        };
        let start = sampler.dirichlet(3);
        let expected_start = [
            0.03382386122172712,
            0.47787368483958914,
            0.48830245393868377,
        ];
        for (got, expected) in start.iter().zip(expected_start) {
            assert!((got - expected).abs() < 2e-15);
        }
        let expected_trans = [
            [0.002063143671635695, 0.7033087878421972, 0.2946280684861671],
            [0.7815487110901412, 0.035926403105419245, 0.1825248858044394],
            [
                0.7497777964230784,
                0.036611391197625014,
                0.21361081237929663,
            ],
        ];
        for expected_row in expected_trans {
            let got = sampler.dirichlet(3);
            for (got, expected) in got.iter().zip(expected_row) {
                assert!((got - expected).abs() < 2e-15);
            }
        }
        let expected_gamma = [
            0.6146697609669675,
            3.103981184016897,
            14.959772659815386,
            0.7188776138920241,
            5.03668016605001,
            0.9764026291891866,
        ];
        for expected in expected_gamma {
            let got = sampler.standard_gamma(2.4) * 1.3;
            assert!((got - expected).abs() < 2e-14);
        }
    }
}
