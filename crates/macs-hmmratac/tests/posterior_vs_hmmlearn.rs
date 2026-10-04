//! Differential harness for [`macs_hmmratac::GaussianHMM::posterior`] against hmmlearn.
//!
//! Upstream's `hmm_predict` is `hmm_model.predict_proba(signals, lens)`, so the only
//! honest check is against hmmlearn's own output on the same observations. `lens`
//! matters: each length is an independent sequence, and `posterior` is called per
//! group exactly as `hmmratac.rs` does it per region mark.
//!
//! Needs a capture from the pinned oracle (`oracle/grab_hmm_probs.py`), so it is
//! ignored by default:
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
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let v: Vec<f64> = l.split_whitespace().map(|t| t.parse().unwrap()).collect();
            assert_eq!(v.len(), want, "bad width in {}", path.display());
            v
        })
        .collect()
}

#[test]
#[ignore = "needs a capture from the pinned oracle"]
fn posterior_matches_hmmlearn() {
    let cap = std::env::var("MACS_HMM_CAP").expect("MACS_HMM_CAP");
    let cap = std::path::PathBuf::from(cap);
    let model_path = std::env::var("MACS_HMM_MODEL").expect("MACS_HMM_MODEL");

    let model = ModelFile::load(std::path::Path::new(&model_path)).expect("model");
    let hmm = model.gaussian();

    let obs = rows(&cap.join("obs.txt"), 4);
    let want = rows(&cap.join("probs.txt"), 3);
    let lens: Vec<usize> = std::fs::read_to_string(cap.join("lens.txt"))
        .expect("lens")
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.trim().parse().unwrap())
        .collect();

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
    let mut dest = String::new();
    for (i, row) in ours.iter().enumerate() {
        if i > 0 {
            dest.push('\n');
        }
        for (j, v) in row.iter().enumerate() {
            if j > 0 {
                dest.push('\t');
            }
            dest.push_str(&format!("{v:.17e}"));
        }
    }
    let out = cap.join("ours_posteriors.tsv");
    std::fs::write(&out, dest).expect("write");

    eprintln!(
        "rows={} max|ours-hmmlearn|={worst:.3e} at flat index {worst_at}",
        ours.len()
    );
    eprintln!("argmax disagreements: {argmax_disagree}/{}", want.len());
    eprintln!("wrote {}", out.display());
    assert!(
        worst <= 1e-9,
        "posterior max abs error {worst:.3e} exceeds 1e-9"
    );
}
