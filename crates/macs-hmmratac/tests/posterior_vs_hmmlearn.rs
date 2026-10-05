//! Differential harness for [`macs_hmmratac::GaussianHMM::posterior`] against hmmlearn.
//!
//! Upstream's `hmm_predict` is `hmm_model.predict_proba(signals, lens)`, so the only
//! honest check is against hmmlearn's own output on the same observations. `lens`
//! matters: each length is an independent sequence, and `posterior` is called per
//! group exactly as `hmmratac.rs` does it per region mark.
//!
//! The committed capture contains non-symmetric transitions and full correlated
//! covariances generated with the pinned hmmlearn public API. Environment
//! overrides allow checking a separate capture:
//!
//! ```text
//! MACS_HMM_CAP=/tmp/hmmcap MACS_HMM_MODEL=/tmp/atr/a_model.json \
//!   cargo test --release -p macs-hmmratac --test posterior_vs_hmmlearn -- --ignored
//! ```

use macs_hmmratac::ModelFile;

fn rows(path: &std::path::Path, want: usize) -> Vec<Vec<f64>> {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    text.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| {
            let v: Vec<f64> = l.split_whitespace().map(|t| t.parse().unwrap()).collect();
            assert_eq!(v.len(), want, "bad width in {}", path.display());
            assert!(
                v.iter().all(|x| x.is_finite()),
                "non-finite value in {}",
                path.display()
            );
            v
        })
        .collect()
}

#[test]
fn posterior_matches_hmmlearn() {
    let default_cap =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/hmmlearn_correlated");
    let cap = std::env::var("MACS_HMM_CAP")
        .map(std::path::PathBuf::from)
        .unwrap_or(default_cap);
    let model_path = std::env::var("MACS_HMM_MODEL")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| cap.join("model.json"));

    let model = ModelFile::load(&model_path).expect("model");
    let hmm = model.gaussian();
    assert_eq!(
        hmm.covariance_type, "full",
        "fixture exercises full covariance"
    );
    assert_ne!(
        hmm.transmat[0][1], hmm.transmat[1][0],
        "fixture transitions are non-symmetric"
    );
    match &hmm.covars {
        macs_hmmratac::Covars::Full(covars) => {
            assert_ne!(covars[2][0][1], 0.0, "fixture includes correlated features");
        }
        macs_hmmratac::Covars::Diag(_) => panic!("fixture must have full covariance"),
    }

    let obs = rows(&cap.join("obs.txt"), 4);
    let want = rows(&cap.join("probs.txt"), 3);
    let lens: Vec<usize> = std::fs::read_to_string(cap.join("lens.txt"))
        .expect("lens")
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.trim().parse().unwrap())
        .collect();

    assert!(!obs.is_empty(), "posterior observation capture is empty");
    assert!(!want.is_empty(), "hmmlearn posterior capture is empty");
    assert!(!lens.is_empty(), "sequence lengths capture is empty");
    assert_eq!(
        obs.len(),
        want.len(),
        "one expected posterior per observation"
    );
    assert!(
        lens.iter().all(|&n| n > 0),
        "sequence lengths must be positive"
    );
    assert_eq!(
        lens.iter().sum::<usize>(),
        obs.len(),
        "lengths must cover the observations"
    );

    let mut ours: Vec<Vec<f64>> = Vec::with_capacity(obs.len());
    let mut at = 0usize;
    for &n in &lens {
        ours.extend(hmm.posterior(&obs[at..at + n]));
        at += n;
    }

    let mut worst = 0.0f64;
    let mut worst_at = 0usize;
    let mut argmax_disagree = 0usize;
    for (i, (a, b)) in ours.iter().zip(&want).enumerate() {
        for (j, (&u, &w)) in a.iter().zip(b.iter()).enumerate() {
            let d = (u - w).abs();
            if d > worst {
                worst = d;
                worst_at = i * 3 + j;
            }
        }
        let am = |v: &[f64]| {
            (0..3)
                .max_by(|&x, &y| v[x].partial_cmp(&v[y]).unwrap())
                .unwrap()
        };
        if am(a) != am(b) {
            argmax_disagree += 1;
        }
    }
    eprintln!(
        "rows={} max|ours-hmmlearn|={worst:.3e} at flat index {worst_at}",
        ours.len()
    );
    eprintln!("argmax disagreements: {argmax_disagree}/{}", want.len());
    assert!(
        worst <= 1e-9,
        "posterior max abs error {worst:.3e} exceeds 1e-9"
    );
}
