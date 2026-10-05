use macs_hmmratac::baum_welch::{fit_gaussian, fit_poisson, gaussian_logpdf, train_gaussian};
use macs_hmmratac::kmeans::kmeans;

#[test]
fn gaussian_log_density_has_full_normalization_constant() {
    let got = gaussian_logpdf(&[0.0], &[0.0], &[vec![1.0]]);
    assert!((got + 0.5 * (2.0 * std::f64::consts::PI).ln()).abs() < 1e-14);
}

#[test]
fn gaussian_log_density_uses_correlations_in_cholesky_factor() {
    let got = gaussian_logpdf(&[1.0, 2.0], &[0.0, 0.0], &[vec![2.0, 1.0], vec![1.0, 2.0]]);
    let expected = -0.5 * (2.0 * (2.0 * std::f64::consts::PI).ln() + 3.0_f64.ln() + 2.0);
    assert!((got - expected).abs() < 1e-14);
}

#[test]
fn gaussian_m_step_centers_covariance_on_updated_mean() {
    let x = vec![vec![9.0], vec![10.0], vec![11.0]];
    let means = vec![vec![0.0]];
    let mut rng = macs_stats::randomstate_from_seed_sequence(1);
    let fitted = fit_gaussian(&x, &[3], &means, &mut rng, 4, 1e-8, 1e-3);
    assert!((fitted.means[0][0] - 10.0).abs() < 1e-10);
    // ML variance is 2/3; hmmlearn's default full-covariance prior contributes
    // 0.01 / occupancy to this diagonal.
    assert!((fitted.covars[0][0][0] - (2.0 / 3.0 + 0.01 / 3.0)).abs() < 1e-8);
}

#[test]
fn poisson_baum_welch_returns_normalized_probabilities_and_finite_rates() {
    let x = vec![
        vec![1.0, 2.0],
        vec![1.0, 1.0],
        vec![2.0, 1.0],
        vec![14.0, 12.0],
        vec![12.0, 15.0],
        vec![13.0, 14.0],
        vec![3.0, 4.0],
        vec![4.0, 3.0],
        vec![3.0, 3.0],
    ];
    let mut rng = macs_stats::randomstate_from_seed_sequence(12345);
    let fitted = fit_poisson(&x, &[3, 3, 3], &mut rng, 25, 1e-5);
    let start_sum: f64 = fitted.log_start.iter().map(|p| p.exp()).sum();
    assert!((start_sum - 1.0).abs() < 1e-12);
    for row in &fitted.log_trans {
        assert!((row.iter().map(|p| p.exp()).sum::<f64>() - 1.0).abs() < 1e-12);
    }
    assert!(fitted
        .lambdas
        .iter()
        .flatten()
        .all(|v| v.is_finite() && *v > 0.0));
    let (min_rate, max_rate) = fitted
        .lambdas
        .iter()
        .map(|row| row[0])
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), rate| {
            (lo.min(rate), hi.max(rate))
        });
    assert!(max_rate - min_rate > 1.0);
}

#[test]
fn kmeans_empty_cluster_relocation_keeps_centers_on_data() {
    let x: Vec<Vec<f64>> = (0..20)
        .map(|i| vec![if i < 10 { 100.0 } else { 200.0 }])
        .collect();
    let mut rng = macs_stats::randomstate_from_seed_sequence(17);
    let (centers, _) = kmeans(&x, 3, 1, 10, 1e-4, &mut rng);
    assert!(centers
        .iter()
        .all(|center| center[0] >= 100.0 && center[0] <= 200.0));
}

fn training_fixture() -> Vec<Vec<f64>> {
    vec![
        vec![2., 4., 3., 1.],
        vec![3., 3., 1., 2.],
        vec![4., 4., 2., 1.],
        vec![2., 3., 3., 2.],
        vec![3., 4., 1., 1.],
        vec![4., 3., 2., 2.],
        vec![2., 4., 3., 1.],
        vec![3., 3., 1., 2.],
        vec![8., 6., 6., 3.],
        vec![9., 5., 4., 4.],
        vec![10., 6., 5., 3.],
        vec![8., 5., 6., 4.],
        vec![9., 6., 4., 3.],
        vec![10., 5., 5., 4.],
        vec![8., 6., 6., 3.],
        vec![9., 5., 4., 4.],
        vec![16., 10., 9., 5.],
        vec![17., 9., 7., 6.],
        vec![18., 10., 8., 5.],
        vec![16., 9., 9., 6.],
        vec![17., 10., 7., 5.],
        vec![18., 9., 8., 6.],
        vec![16., 10., 9., 5.],
        vec![17., 9., 7., 6.],
    ]
}

#[test]
fn gaussian_training_matches_pinned_hmmlearn_fixture() {
    // Captured from hmmlearn 0.3.3 / sklearn 1.9.1 on training_fixture with
    // RandomState(MT19937(SeedSequence(12345))), lengths [8, 8, 8], n_iter=3,
    // tol=0, covariance_type=full.
    let x = training_fixture();
    let mut rng = macs_stats::randomstate_from_seed_sequence(12345);
    let fitted = train_gaussian(&x, &[8, 8, 8], 3, &mut rng, 3, 0.0, 1e-3);
    let expected_means = [
        [
            6.262400596479795,
            4.1116811818253005,
            3.35471875499051,
            3.0670087090721054,
        ],
        [
            5.749999981838459,
            4.999999992697806,
            3.7499999859350943,
            1.999999992697805,
        ],
        [
            16.85459434007683,
            9.511690328525738,
            8.001442804346734,
            5.488309625823207,
        ],
    ];
    for (got, expected) in fitted.means.iter().zip(expected_means) {
        for (got, expected) in got.iter().zip(expected) {
            assert!((got - expected).abs() < 1e-5, "means {got} != {expected}");
        }
    }
    let expected_covars = [
        [
            [
                12.309841529199469,
                4.216911144669883,
                5.358636975800619,
                3.7038320348800773,
            ],
            [
                4.216911144669883,
                1.5248190662439205,
                1.9796163859129212,
                1.306445776991586,
            ],
            [
                5.358636975800619,
                1.9796163859129212,
                3.357225823115789,
                1.7748569215591925,
            ],
            [
                3.7038320348800773,
                1.306445776991586,
                1.7748569215591925,
                1.175421803500741,
            ],
        ],
        [
            [
                9.688750027460515,
                3.0012500050430337,
                4.063749998348536,
                3.0012500050430373,
            ],
            [
                3.0012500050430337,
                1.0012500004324112,
                1.50124999738433,
                1.0012500004324236,
            ],
            [
                4.063749998348536,
                1.50124999738433,
                2.938749991309191,
                1.5012499973843336,
            ],
            [
                3.0012500050430373,
                1.0012500004324236,
                1.5012499973843336,
                1.0012500004324258,
            ],
        ],
        [
            [
                0.5992970168858146,
                -0.0522404921368258,
                -0.38079363125904375,
                0.05479933680326821,
            ],
            [
                -0.0522404921368258,
                0.25114265356029475,
                0.1284635982276283,
                -0.2485840413704457,
            ],
            [
                -0.38079363125904375,
                0.1284635982276283,
                0.7643999729180997,
                -0.1259050395168978,
            ],
            [
                0.05479933680326821,
                -0.2485840413704457,
                -0.1259050395168978,
                0.25114260660074167,
            ],
        ],
    ];
    for (got, expected) in fitted.covars.iter().zip(expected_covars) {
        for (got_row, expected_row) in got.iter().zip(expected) {
            for (got, expected) in got_row.iter().zip(expected_row) {
                assert!((got - expected).abs() < 1e-5, "covar {got} != {expected}");
            }
        }
    }
    let expected_start = [
        6.198983784522477e-18,
        0.6666666554675527,
        0.33333334453244734,
    ];
    for (got, expected) in fitted.log_start.iter().zip(expected_start) {
        assert!((got.exp() - expected).abs() < 1e-6);
    }
    let expected_trans = [
        [4.428744080511974e-17, 0.97324465978355, 0.02675534021644989],
        [0.9999999999523642, 0.0, 4.763578569348484e-11],
        [
            0.02674044314303273,
            3.302356089610097e-11,
            0.9732595568239437,
        ],
    ];
    for (got_row, expected_row) in fitted.log_trans.iter().zip(expected_trans) {
        for (got, expected) in got_row.iter().zip(expected_row) {
            assert!((got.exp() - expected).abs() < 1e-6);
        }
    }
}

#[test]
fn poisson_training_matches_pinned_hmmlearn_fixture() {
    // Captured from hmmlearn 0.3.3 PoissonHMM with the same RandomState/X;
    // n_iter=10 and tol=1e-2.
    let x = training_fixture();
    let mut rng = macs_stats::randomstate_from_seed_sequence(12345);
    let poisson = fit_poisson(&x, &[8, 8, 8], &mut rng, 10, 1e-2);
    let expected_lambdas = [
        [
            2.9708121435541925,
            3.4339350188033837,
            1.8887174707878345,
            1.566064981668229,
        ],
        [
            12.875000000435369,
            7.500000000242076,
            6.500000000217246,
            4.500000000076585,
        ],
        [
            2.874999934693468,
            3.500000046277479,
            2.000000077853836,
            1.499999954583119,
        ],
    ];
    for (got, expected) in poisson.lambdas.iter().zip(expected_lambdas) {
        for (got, expected) in got.iter().zip(expected) {
            assert!((got - expected).abs() < 1e-5, "lambda {got} != {expected}");
        }
    }
}
