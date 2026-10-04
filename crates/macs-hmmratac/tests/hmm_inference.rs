//! Differential test of the HMM inference against hmmlearn.
//!
//! The acceptance criteria make HMMRATAC's *self-training* the one declared
//! deviation (reproducing an EM fit means reproducing an SVD, a QR and a
//! specific RNG draw sequence). Inference must still be exact, and this test is
//! what pins that: the reference posteriors here were produced by
//! `hmmlearn.GaussianHMM.predict_proba` on a model fitted from real
//! `*_model.json` parameters, and `macs_hmmratac` must reproduce them.
//!
//! The fixtures are checked in rather than generated, so this test needs no
//! Python at run time.

use macs_hmmratac::{Covars, GaussianHMM};

/// A fitted 3-state, 4-feature model: start probabilities, transition matrix,
/// means and full covariances, exactly as `hmm_model_save` writes them.
fn model() -> GaussianHMM {
    GaussianHMM {
        startprob: vec![
            0.00011037684098707946,
            0.33333333333333337,
            0.6665562898256796,
        ],
        transmat: vec![
            vec![6.92650648e-05, 7.37424048e-32, 9.99930735e-01],
            vec![5.53187914e-16, 1.00000000e+00, 9.04520895e-59],
            vec![5.07831454e-02, 8.46276847e-88, 9.49216855e-01],
        ],
        means: vec![
            vec![0.061352, 0.073463, 0.027539, 0.059180],
            vec![0.300755, 0.253242, 0.118908, 0.057968],
            vec![0.066147, 0.039851, 0.034064, 0.024259],
        ],
        covars: Covars::Full(vec![
            vec![
                vec![0.00253679, 0.00059579, 0.00120496, 0.00103350],
                vec![0.00059579, 0.00291602, 0.00120239, 0.00039461],
                vec![0.00120496, 0.00120239, 0.00116707, 0.00088898],
                vec![0.00103350, 0.00039461, 0.00088898, 0.00181530],
            ],
            vec![
                vec![0.00092037, 0.00007274, 0.00008648, 0.00008952],
                vec![0.00007274, 0.00098045, 0.00002475, -0.00012804],
                vec![0.00008648, 0.00002475, 0.00081038, 0.00009909],
                vec![0.00008952, -0.00012804, 0.00009909, 0.00095274],
            ],
            vec![
                vec![0.00116812, 0.00032385, 0.00025317, 0.00009402],
                vec![0.00032385, 0.00095332, 0.00033852, 0.00011061],
                vec![0.00025317, 0.00033852, 0.00094496, 0.00009610],
                vec![0.00009402, 0.00011061, 0.00009610, 0.00055649],
            ],
        ]),
        covariance_type: "full".into(),
    }
}

/// hmmlearn's `predict_proba` on a 10-bin observation sequence.
const REFERENCE: &[[f64; 3]] = &[
    [1.7846270e-06, 9.9823923e-46, 9.9999822e-01],
    [7.7469041e-05, 9.1366107e-46, 9.9992253e-01],
    [9.9999862e-01, 1.4442218e-22, 1.3818065e-06],
    [9.9893158e-01, 1.4442218e-22, 1.0684236e-03],
    [1.3818244e-03, 4.8660399e-40, 9.9861818e-01],
    [4.4440239e-02, 4.0300326e-66, 9.5555976e-01],
    [7.3195811e-06, 7.2035194e-71, 9.9999268e-01],
    [9.9999998e-01, 8.9604590e-43, 1.8280984e-08],
    [7.1720659e-06, 1.8041708e-63, 9.9999283e-01],
    [1.0627671e-02, 5.8566363e-66, 9.8937233e-01],
];

/// Observations that produced [`REFERENCE`], printed by the generator.
const OBS: &[[f64; 4]] = &[
    [0.05832779, 0.0001, 0.02371375, 0.0001],
    [0.02397781, 0.04657704, 0.0001, 0.0001],
    [0.2865066, 0.2993998, 0.1684696, 0.04175152],
    [0.3299857, 0.2673675, 0.1278294, 0.07063165],
    [0.05821674, 0.08400132, 0.08179904, 0.1109973],
    [0.1074559, 0.0001, 0.01378695, 0.03990314],
    [0.05195207, 0.0001, 0.005681001, 0.0001],
    [0.2983716, 0.246182, 0.1474714, 0.1037851],
    [0.1003591, 0.0265411, 0.02395469, 0.009991651],
    [0.04817165, 0.007063554, 0.0001, 0.0001],
];

#[test]
fn posteriors_track_hmmlearn_to_the_precision_that_matters() {
    let obs: Vec<Vec<f64>> = OBS.iter().map(|r| r.to_vec()).collect();
    let got = model().posterior(&obs);
    assert_eq!(got.len(), REFERENCE.len());
    let mut worst_tv: f64 = 0.0;
    for (i, (mine, want)) in got.iter().zip(REFERENCE.iter()).enumerate() {
        let row: f64 = mine.iter().sum();
        assert!((row - 1.0).abs() < 1e-9, "row {i} sums to {row}");
        // total-variation distance: what actually decides the state path, and
        // therefore every hmmratac output
        let tv: f64 = mine
            .iter()
            .zip(want.iter())
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            / 2.0;
        worst_tv = worst_tv.max(tv);
        assert!(
            tv < 1e-5,
            "frame {i}: total-variation distance {tv:e}\n  mine {mine:?}\n  ref  {:?}",
            want.to_vec()
        );
    }
    println!("worst total-variation distance: {worst_tv:e}");
}

#[test]
fn every_row_picks_the_same_state_as_hmmlearn() {
    let obs: Vec<Vec<f64>> = OBS.iter().map(|r| r.to_vec()).collect();
    let got = model().posterior(&obs);
    for (i, (mine, want)) in got.iter().zip(REFERENCE.iter()).enumerate() {
        let a = mine
            .iter()
            .enumerate()
            .max_by(|x, y| x.1.partial_cmp(y.1).unwrap())
            .map(|(k, _)| k)
            .unwrap();
        let b = want
            .iter()
            .enumerate()
            .max_by(|x, y| x.1.partial_cmp(y.1).unwrap())
            .map(|(k, _)| k)
            .unwrap();
        assert_eq!(a, b, "frame {i} picks state {a}, hmmlearn picks {b}");
    }
}
