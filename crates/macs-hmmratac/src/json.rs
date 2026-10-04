//! Byte-compatible JSON serialisation of a [`ModelFile`].
//!
//! Upstream uses `json.dump` with default arguments, which means:
//!
//! * **no** indentation and **no** trailing newline,
//! * `", "` between items -- one space, unusual but that is `json`'s default
//!   separator pair (`item_separator=", "`),
//! * `repr()` for floats, i.e. the shortest string that round-trips, and
//!   `NaN` / `Infinity` written as bare tokens.
//!
//! Getting this wrong does not change what the model *means* -- `hmm_model_init`
//! reads it back with `json.load` -- but it does change the bytes, and
//! `*_model.json` is one of the files a golden comparison looks at.

use crate::{Covars, HmmType, ModelFile};

/// Serialise one `f64` the way Python's `json` encoder does.
///
/// Python uses `float.__repr__`, the shortest decimal that round-trips, and
/// always keeps a `.0` on integral values so the token reads back as a float.
/// Rust's `{:?}` is also shortest-round-trip, so the two agree on digits; Rust
/// only differs in dropping the `.0`.
fn num(s: &mut String, v: f64) {
    if v.is_nan() {
        s.push_str("NaN");
    } else if v.is_infinite() {
        s.push_str(if v > 0.0 { "Infinity" } else { "-Infinity" });
    } else {
        let d = format!("{v:?}");
        s.push_str(&d);
        if !d.contains(['.', 'e', 'E', 'n', 'i']) {
            // `3` -> `3.0`, matching `repr(3.0)`
            s.push_str(".0");
        }
    }
}

fn vec(s: &mut String, v: &[f64]) {
    s.push('[');
    for (i, x) in v.iter().enumerate() {
        if i > 0 {
            s.push_str(", ");
        }
        num(s, *x);
    }
    s.push(']');
}

fn mat(s: &mut String, m: &[Vec<f64>]) {
    s.push('[');
    for (i, row) in m.iter().enumerate() {
        if i > 0 {
            s.push_str(", ");
        }
        vec(s, row);
    }
    s.push(']');
}

fn cube(s: &mut String, m: &[Vec<Vec<f64>>]) {
    s.push('[');
    for (i, slab) in m.iter().enumerate() {
        if i > 0 {
            s.push_str(", ");
        }
        mat(s, slab);
    }
    s.push(']');
}

fn string(s: &mut String, v: &str) {
    s.push('"');
    for c in v.chars() {
        match c {
            '"' => s.push_str("\\\""),
            '\\' => s.push_str("\\\\"),
            '\n' => s.push_str("\\n"),
            '\r' => s.push_str("\\r"),
            '\t' => s.push_str("\\t"),
            c => s.push(c),
        }
    }
    s.push('"');
}

fn key(s: &mut String, k: &str) {
    string(s, k);
    // Python's default `item_separator`/`key_separator` pair is `(", ", ": ")`,
    // so there is a space after the colon too.
    s.push_str(": ");
}

/// Serialise `m` in upstream's exact byte layout.
pub fn write_model(m: &ModelFile) -> String {
    let mut s = String::new();
    s.push('{');
    key(&mut s, "startprob");
    vec(&mut s, &m.startprob);
    s.push_str(", ");
    key(&mut s, "transmat");
    mat(&mut s, &m.transmat);
    match m.hmm_type {
        HmmType::Gaussian => {
            s.push_str(", ");
            key(&mut s, "means");
            mat(&mut s, &m.means);
            s.push_str(", ");
            key(&mut s, "covars");
            match &m.covars {
                Covars::Diag(d) => mat(&mut s, d),
                Covars::Full(f) => cube(&mut s, f),
            }
            s.push_str(", ");
            key(&mut s, "covariance_type");
            string(&mut s, &m.covariance_type);
        }
        HmmType::Poisson => {
            s.push_str(", ");
            key(&mut s, "lambdas");
            mat(&mut s, &m.lambdas);
        }
    }
    s.push_str(&format!(", \"n_features\": {}", m.n_features));
    s.push_str(&format!(", \"i_open_region\": {}", m.i_open_region));
    s.push_str(&format!(
        ", \"i_background_region\": {}",
        m.i_background_region
    ));
    s.push_str(&format!(
        ", \"i_nucleosomal_region\": {}",
        m.i_nucleosomal_region
    ));
    s.push_str(&format!(", \"hmm_binsize\": {}", m.hmm_binsize));
    s.push_str(", ");
    key(&mut s, "hmm_type");
    string(
        &mut s,
        match m.hmm_type {
            HmmType::Gaussian => "gaussian",
            HmmType::Poisson => "poisson",
        },
    );
    s.push('}');
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A Gaussian model with the magnitudes a real trained HMM has.
    fn model() -> ModelFile {
        ModelFile {
            hmm_type: HmmType::Gaussian,
            startprob: vec![0.6666666666666666, 1e-88, 0.3333333333333333],
            transmat: vec![
                vec![1.0, 0.0, 0.0],
                vec![5.53e-16, 1.0, 0.0],
                vec![0.05344, 7.0e-32, 0.94656],
            ],
            means: vec![vec![1.0, 2.0, 3.0, 4.0]; 3],
            covars: Covars::Diag(vec![vec![1.0, 1.0, 1.0, 1.0]; 3]),
            covariance_type: "diag".into(),
            lambdas: Vec::new(),
            i_open_region: 2,
            i_background_region: 0,
            i_nucleosomal_region: 1,
            hmm_binsize: 500,
            n_features: 4,
        }
    }

    #[test]
    fn integral_floats_keep_their_decimal_point() {
        let mut s = String::new();
        num(&mut s, 1.0);
        assert_eq!(s, "1.0");
        s.clear();
        num(&mut s, 3.0);
        assert_eq!(s, "3.0");
    }

    #[test]
    fn exponent_notation_is_kept() {
        let mut s = String::new();
        num(&mut s, 5.53e-16);
        assert!(s.contains('e'), "got {s}");
    }

    #[test]
    fn nan_and_infinity_use_pythons_bare_tokens() {
        let mut s = String::new();
        num(&mut s, f64::NAN);
        assert_eq!(s, "NaN");
        s.clear();
        num(&mut s, f64::INFINITY);
        assert_eq!(s, "Infinity");
        s.clear();
        num(&mut s, f64::NEG_INFINITY);
        assert_eq!(s, "-Infinity");
    }

    #[test]
    fn separators_match_pythons_json_default() {
        let j = write_model(&model());
        assert!(j.contains("[1.0, 1.0, 1.0, 1.0]"), "comma-space: {j}");
        assert!(!j.contains('\n'), "no pretty printing");
    }

    #[test]
    fn key_order_matches_the_upstream_literal() {
        let j = write_model(&model());
        let order: Vec<usize> = [
            "\"startprob\"",
            "\"transmat\"",
            "\"means\"",
            "\"covars\"",
            "\"covariance_type\"",
            "\"n_features\"",
            "\"i_open_region\"",
            "\"i_background_region\"",
            "\"i_nucleosomal_region\"",
            "\"hmm_binsize\"",
            "\"hmm_type\"",
        ]
        .iter()
        .map(|k| j.find(k).unwrap_or_else(|| panic!("missing {k}")))
        .collect();
        assert!(order.windows(2).all(|w| w[0] < w[1]), "{order:?}");
    }

    #[test]
    fn round_trips_through_from_json() {
        let m = model();
        let back = ModelFile::from_json(&m.to_json()).unwrap();
        assert_eq!(back, m);
    }

    fn poisson_model() -> ModelFile {
        let mut m = model();
        m.hmm_type = HmmType::Poisson;
        m.lambdas = vec![vec![1.0, 2.0, 3.0, 4.0]; 3];
        m.means = Vec::new();
        m.covars = Covars::Diag(Vec::new());
        m
    }

    #[test]
    fn a_poisson_model_writes_lambdas_and_no_covars() {
        let m = poisson_model();
        let j = m.to_json();
        assert!(j.contains("\"lambdas\""));
        assert!(!j.contains("\"covars\""));
        assert_eq!(ModelFile::from_json(&j).unwrap(), m);
    }

    #[test]
    fn a_file_without_hmm_type_defaults_to_gaussian() {
        let m = model();
        let full = m.to_json();
        let j = full.replace(", \"hmm_type\": \"gaussian\"", "");
        assert!(!j.contains("\"hmm_type\""), "{j}");
        let back = ModelFile::from_json(&j).unwrap();
        assert_eq!(back.hmm_type, HmmType::Gaussian);
    }

    #[test]
    fn full_covariance_round_trips() {
        let mut m = model();
        m.covars = Covars::Full(vec![vec![vec![1.0, 0.0], vec![0.0, 1.0]]]);
        m.covariance_type = "full".into();
        let back = ModelFile::from_json(&m.to_json()).unwrap();
        assert_eq!(back, m);
        assert!(matches!(back.covars, Covars::Full(_)));
    }

    #[test]
    fn covariance_depth_disambiguates_when_covariance_type_is_absent() {
        let diag = r#"{"startprob":[1.0],"transmat":[[1.0]],
            "means":[[1.0,2.0]],"covars":[[1.0,1.0]],
            "n_features":2,"i_open_region":0,"i_background_region":0,
            "i_nucleosomal_region":0,"hmm_binsize":500,"hmm_type":"gaussian"}"#;
        assert!(matches!(
            ModelFile::from_json(diag).unwrap().covars,
            Covars::Diag(_)
        ));
        let full = r#"{"startprob":[1.0],"transmat":[[1.0]],
            "means":[[1.0,2.0]],"covars":[[[1.0,0.0],[0.0,1.0]]],
            "n_features":2,"i_open_region":0,"i_background_region":0,
            "i_nucleosomal_region":0,"hmm_binsize":500,"hmm_type":"gaussian"}"#;
        assert!(matches!(
            ModelFile::from_json(full).unwrap().covars,
            Covars::Full(_)
        ));
    }
}
