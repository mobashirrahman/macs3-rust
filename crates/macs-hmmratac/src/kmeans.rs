//! `sklearn.cluster.KMeans` as hmmlearn invokes it for HMM initialisation.
//!
//! `HMMR_HMM.hmm_training` fits `GaussianHMM`/`PoissonHMM` with
//! `random_state=RandomState(MT19937(SeedSequence(seed)))`. hmmlearn's
//! `GaussianHMM._init` runs `KMeans(n_clusters=3, random_state=rs, n_init=10)`
//! (sklearn defaults: `k-means++`, `max_iter=300`, `tol=1e-4`, Lloyd's) and takes
//! the cluster centers as initial means. Covariances start from the global
//! `np.cov(X.T)`, not per-cluster.
//!
//! This implements the algorithm functionally (k-means++ seeding per sklearn's
//! greedy variant, Lloyd iterations to `tol`, best-of-`n_init` by inertia). It is
//! validated by converging to the same clusters hmmlearn reaches on the training
//! data, not by bit-identical intermediate draws.

/// Squared Euclidean distance between two vectors.
fn sq_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| (x - y) * (x - y)).sum()
}

/// k-means++ seeding (sklearn greedy variant).
///
/// First center uniform; each subsequent center is the best of `2+ln(k)` candidates
/// sampled proportional to squared distance, where "best" maximises the potential
/// (total distance) reduction. Returns the centers and their row indices.
fn kmeans_plusplus(
    x: &[Vec<f64>],
    n_clusters: usize,
    rng: &mut macs_stats::NumpyRng,
) -> Vec<Vec<f64>> {
    kmeans_plusplus_indexed(x, n_clusters, rng).0
}

/// As above, also returning chosen row indices (for differential debugging).
pub fn kmeans_plusplus_indexed(
    x: &[Vec<f64>],
    n_clusters: usize,
    rng: &mut macs_stats::NumpyRng,
) -> (Vec<Vec<f64>>, Vec<usize>) {
    let n = x.len();
    let n_trials = 2 + (n_clusters as f64).ln() as usize;
    let mut centers: Vec<Vec<f64>> = Vec::with_capacity(n_clusters);
    let mut closest_sq = vec![f64::INFINITY; n];

    // `random_state.choice(n, p=uniform)`: draw a uniform float and take the first
    // index whose cumulative weight reaches it (`searchsorted`). NOT `randint` --
    // that consumes different RNG words and picks a different first center (2781 vs
    // 16912 here), flipping the whole downstream basin.
    let r0 = rng.random_sample() * n as f64;
    // Uniform weights: cumsum is 1,2,...,n; first index with cumsum > r0, i.e.
    // floor(r0), clamped. (`searchsorted` side='left' gives first `cumsum >= r`;
    // with unit weights that is `ceil(r)` for non-integer r, i.e. `floor(r)` as an
    // index, except r integer-valued which never happens for a continuous draw.)
    let first_idx = (r0.floor() as usize).min(n - 1);
    let mut indices = vec![first_idx];
    centers.push(x[first_idx].clone());
    for (i, row) in x.iter().enumerate() {
        closest_sq[i] = sq_dist(row, &centers[0]);
    }

    for _ in 1..n_clusters {
        let total: f64 = closest_sq.iter().sum();
        let mut best_idx = 0usize;
        let mut best_gain = f64::NEG_INFINITY;
        for _ in 0..n_trials {
            // `uniform(size) * current_pot` then `searchsorted(cumsum, ...)`.
            let r = rng.random_sample() * total;
            let mut cum = 0.0;
            let mut cand = n - 1;
            for (i, &w) in closest_sq.iter().enumerate() {
                cum += w;
                if cum >= r {
                    cand = i;
                    break;
                }
            }
            let mut gain = 0.0;
            for (i, row) in x.iter().enumerate() {
                let d = sq_dist(row, &x[cand]);
                if d < closest_sq[i] {
                    gain += closest_sq[i] - d;
                }
            }
            if gain > best_gain {
                best_gain = gain;
                best_idx = cand;
            }
        }
        centers.push(x[best_idx].clone());
        indices.push(best_idx);
        for (i, row) in x.iter().enumerate() {
            let d = sq_dist(row, &x[best_idx]);
            if d < closest_sq[i] {
                closest_sq[i] = d;
            }
        }
    }
    (centers, indices)
}

/// One Lloyd run from `init` centers. Returns `(centers, inertia, n_iter)`.
///
/// Iterates assign-then-update until sklearn's scaled center shift is below
/// tolerance, or labels stop changing. Empty clusters are reseeded from the
/// samples farthest from their assigned center, as in sklearn's Lloyd kernel.
fn lloyd(
    x: &[Vec<f64>],
    mut centers: Vec<Vec<f64>>,
    max_iter: usize,
    tol: f64,
) -> (Vec<Vec<f64>>, f64, usize) {
    let k = centers.len();
    let dim = x[0].len();
    let mut labels = vec![0usize; x.len()];
    // sklearn compares squared center movement against tol * mean(var(X)).
    let mut feature_means = vec![0.0; dim];
    for row in x {
        for (mean, value) in feature_means.iter_mut().zip(row) {
            *mean += value;
        }
    }
    for mean in &mut feature_means {
        *mean /= x.len() as f64;
    }
    let mut mean_variance = 0.0;
    for (j, mean) in feature_means.iter().enumerate() {
        mean_variance += x.iter().map(|row| (row[j] - mean).powi(2)).sum::<f64>() / x.len() as f64;
    }
    mean_variance /= dim as f64;
    let scaled_tol = tol * mean_variance;
    let mut old_labels = vec![usize::MAX; x.len()];
    for it in 0..max_iter {
        // Assign.
        for (i, row) in x.iter().enumerate() {
            let (mut best, mut bd) = (0usize, f64::INFINITY);
            for (j, c) in centers.iter().enumerate() {
                let d = sq_dist(row, c);
                if d < bd {
                    bd = d;
                    best = j;
                }
            }
            labels[i] = best;
        }
        let labels_unchanged = labels == old_labels;
        old_labels.clone_from_slice(&labels);
        // Update.
        let mut new_centers = vec![vec![0.0; dim]; k];
        let mut counts = vec![0usize; k];
        for (i, row) in x.iter().enumerate() {
            let j = labels[i];
            counts[j] += 1;
            for (d, v) in new_centers[j].iter_mut().zip(row.iter()) {
                *d += v;
            }
        }
        // Relocate empty clusters using distinct samples with the largest
        // residuals. Moving each sample's mass out of its donor preserves the
        // actual mean of every nonempty cluster.
        let mut empty: Vec<usize> = counts
            .iter()
            .enumerate()
            .filter_map(|(j, &count)| (count == 0).then_some(j))
            .collect();
        if !empty.is_empty() {
            let mut farthest: Vec<(f64, usize)> = x
                .iter()
                .enumerate()
                .map(|(i, row)| (sq_dist(row, &centers[labels[i]]), i))
                .collect();
            farthest.sort_by(|a, b| b.0.total_cmp(&a.0));
            for (rank, cluster) in empty.drain(..).enumerate() {
                let Some(&(_, sample)) = farthest.get(rank) else {
                    break;
                };
                let donor = labels[sample];
                if counts[donor] <= 1 {
                    continue;
                }
                counts[donor] -= 1;
                counts[cluster] = 1;
                for d in 0..dim {
                    new_centers[donor][d] -= x[sample][d];
                    new_centers[cluster][d] = x[sample][d];
                }
            }
        }
        let mut shift = 0.0;
        for j in 0..k {
            if counts[j] == 0 {
                // More clusters than samples is invalid for sklearn too; keep
                // a finite center in this defensive path.
                new_centers[j].clone_from(&centers[j]);
                continue;
            }
            for v in new_centers[j].iter_mut() {
                *v /= counts[j] as f64;
            }
            shift += sq_dist(&centers[j], &new_centers[j]);
        }
        centers = new_centers;
        if labels_unchanged || shift <= scaled_tol {
            // Recompute inertia at the new centers for the reported value.
            let mut final_inertia = 0.0;
            for row in x.iter() {
                let mut bd = f64::INFINITY;
                for c in centers.iter() {
                    let d = sq_dist(row, c);
                    if d < bd {
                        bd = d;
                    }
                }
                final_inertia += bd;
            }
            return (centers, final_inertia, it + 1);
        }
    }
    // Max iter: report last inertia.
    let mut inertia = 0.0;
    for row in x {
        inertia += centers
            .iter()
            .map(|c| sq_dist(row, c))
            .fold(f64::INFINITY, f64::min);
    }
    (centers, inertia, max_iter)
}

/// Full KMeans: `n_init` restarts, keep lowest inertia.
///
/// `rng` must already be the `RandomState(MT19937(SeedSequence(seed)))` stream;
/// each restart draws k-means++ seeds sequentially from it, as sklearn does.
pub fn kmeans(
    x: &[Vec<f64>],
    n_clusters: usize,
    n_init: usize,
    max_iter: usize,
    tol: f64,
    rng: &mut macs_stats::NumpyRng,
) -> (Vec<Vec<f64>>, f64) {
    let mut best_centers = Vec::new();
    let mut best_inertia = f64::INFINITY;
    for _ in 0..n_init {
        let init = kmeans_plusplus(x, n_clusters, rng);
        let (centers, inertia, _) = lloyd(x, init, max_iter, tol);
        if inertia < best_inertia {
            best_inertia = inertia;
            best_centers = centers;
        }
    }
    (best_centers, best_inertia)
}
