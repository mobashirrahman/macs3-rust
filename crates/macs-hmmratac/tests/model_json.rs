//! `*_model.json` read/write parity with `MACS3.Signal.HMMR_HMM`.
//!
//! The reference files in `tests/data/` were produced by CPython's
//! `json.dumps` over the exact key order of `hmm_model_save`'s literal, so this
//! pins the byte layout -- `", "` item separators, `": "` key separators, no
//! indentation, no trailing newline, and `repr()` float formatting -- as well as
//! the numeric round trip.
//!
//! The round trip matters more than it looks. `serde_json`'s float parser is not
//! correctly rounded: it parses `1e-88` as `1.0000000000000001e-88`, and across
//! 6000 values in the exponent range hmmlearn produces it disagrees with Python's
//! `float()` on 1335 of them. A trained ATAC HMM really does put entries that
//! small in `startprob_`, so `macs-hmmratac`'s reader parses with Rust's
//! correctly-rounded `f64::from_str` instead.

use std::path::Path;

use macs_hmmratac::{Covars, HmmType, ModelFile};

fn data(name: &str) -> String {
    let p = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name);
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()))
}

#[test]
fn gaussian_diag_file_rewrites_byte_identically() {
    let m = ModelFile::from_json(&data("model_gaussian_diag.json")).unwrap();
    assert_eq!(m.hmm_type, HmmType::Gaussian);
    assert!(matches!(m.covars, Covars::Diag(_)));
    assert_eq!(m.covariance_type, "diag");
    assert_eq!(m.to_json(), data("model_gaussian_diag.json"));
}

#[test]
fn gaussian_full_file_rewrites_byte_identically() {
    let m = ModelFile::from_json(&data("model_gaussian_full.json")).unwrap();
    assert!(matches!(m.covars, Covars::Full(_)));
    assert_eq!(m.covariance_type, "full");
    assert_eq!(m.to_json(), data("model_gaussian_full.json"));
}

#[test]
fn poisson_file_rewrites_byte_identically() {
    let m = ModelFile::from_json(&data("model_poisson.json")).unwrap();
    assert_eq!(m.hmm_type, HmmType::Poisson);
    assert_eq!(m.to_json(), data("model_poisson.json"));
}

#[test]
fn a_one_ulp_misparse_would_be_caught() {
    // If the reader used serde_json, this assertion would fail: serde_json maps
    // `1e-88` to 1.0000000000000001e-88.
    let m = ModelFile::from_json(&data("model_gaussian_diag.json")).unwrap();
    assert_eq!(m.startprob[1], 1e-88f64);
    assert_eq!(m.startprob[1].to_bits(), 1e-88f64.to_bits());
    assert_eq!(m.transmat[1][0], 5.53e-16f64);
    assert_eq!(m.transmat[2][1], 7e-32f64);
}

#[test]
fn state_indices_and_binsize_survive() {
    let m = ModelFile::from_json(&data("model_gaussian_diag.json")).unwrap();
    assert_eq!(m.i_open_region, 2);
    assert_eq!(m.i_background_region, 0);
    assert_eq!(m.i_nucleosomal_region, 1);
    assert_eq!(m.hmm_binsize, 500);
    assert_eq!(m.n_features, 4);
}

#[test]
fn state_labels_follow_the_recorded_indices() {
    let m = ModelFile::from_json(&data("model_gaussian_diag.json")).unwrap();
    assert_eq!(m.state_label(m.i_open_region), "open");
    assert_eq!(m.state_label(m.i_nucleosomal_region), "nuc");
    assert_eq!(m.state_label(m.i_background_region), "bg");
}

#[test]
fn a_missing_field_is_an_error_not_a_panic() {
    for bad in [
        "{}",
        r#"{"startprob":[1.0]}"#,
        r#"{"startprob":[1.0],"transmat":[[1.0]]}"#,
        r#"{"startprob":[1.0],"transmat":[[1.0]],"hmm_type":"quantum"}"#,
    ] {
        assert!(
            ModelFile::from_json(bad).is_err(),
            "should have been rejected: {bad}"
        );
    }
}
