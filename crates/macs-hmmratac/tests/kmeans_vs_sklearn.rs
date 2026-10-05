//! KMeans against a fixed oracle result from the pinned sklearn version.

use macs_hmmratac::kmeans::kmeans;

/// Recorded from sklearn 1.9.1 KMeans(n_clusters=3,n_init=10,max_iter=300,
/// tol=1e-4) with RandomState(MT19937(SeedSequence(12345))).
#[test]
fn kmeans_matches_sklearn_seeded_oracle() {
    let x = vec![
        vec![0.0, 0.0],
        vec![0.2, 0.0],
        vec![0.0, 0.2],
        vec![10.0, 10.0],
        vec![10.2, 10.0],
        vec![10.0, 10.2],
        vec![20.0, 0.0],
        vec![20.2, 0.0],
        vec![20.0, 0.2],
        vec![10.3, 10.1],
        vec![0.1, 0.2],
        vec![20.1, 0.1],
    ];
    let mut rng = macs_stats::randomstate_from_seed_sequence(12345);
    let (mut centers, inertia) = kmeans(&x, 3, 10, 300, 1e-4, &mut rng);
    assert_eq!(centers.len(), 3);
    centers.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let expected = [[0.075, 0.1], [10.125, 10.075], [20.075, 0.075]];
    for (center, expected) in centers.iter().zip(expected) {
        assert!((center[0] - expected[0]).abs() < 1e-12);
        assert!((center[1] - expected[1]).abs() < 1e-12);
    }
    assert!((inertia - 0.2175).abs() < 1e-12);
}
