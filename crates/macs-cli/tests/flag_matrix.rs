//! The Rust subcommand list must agree with the auto-derived upstream flag matrix.
//!
//! `oracle/flag_matrix.tsv` is generated from upstream's argparse by
//! `oracle/gen_flag_matrix.py`, so it is the authority on which subcommands exist.
//! This test fails if the port gains, loses, or renames one without the matrix
//! being regenerated.

use std::collections::{BTreeMap, BTreeSet};

/// `oracle/flag_matrix.tsv`, relative to the repository root.
fn matrix_path() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("oracle")
        .join("flag_matrix.tsv")
}

#[test]
fn the_matrix_exists_and_has_a_header() {
    let path = matrix_path();
    assert!(
        path.is_file(),
        "{} missing; run oracle/gen_flag_matrix.py",
        path.display()
    );
    let text = std::fs::read_to_string(&path).expect("readable");
    let mut lines = text.lines();
    assert_eq!(
        lines.next(),
        Some("subcommand\tflag\taction\tnargs\tdest\trequired\tdefault\tchoices\thelp"),
        "unexpected header"
    );
}

#[test]
fn the_rust_subcommand_list_matches_the_matrix() {
    let text = std::fs::read_to_string(matrix_path()).expect("readable");
    let from_matrix: BTreeSet<&str> = text
        .lines()
        .skip(1)
        .filter_map(|line| line.split('\t').next())
        .filter(|sub| *sub != "<top>")
        .collect();

    let from_rust: BTreeSet<&str> = macs_cli::SUBCOMMANDS.iter().copied().collect();

    assert_eq!(
        from_rust, from_matrix,
        "subcommand set differs from oracle/flag_matrix.tsv"
    );
    assert_eq!(from_rust.len(), 14, "expected 14 subcommands");
}

#[test]
fn every_subcommand_carries_at_least_one_option() {
    let text = std::fs::read_to_string(matrix_path()).expect("readable");
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for line in text.lines().skip(1) {
        if let Some(sub) = line.split('\t').next() {
            if sub != "<top>" {
                *counts.entry(sub).or_default() += 1;
            }
        }
    }
    for sub in macs_cli::SUBCOMMANDS {
        let n = counts.get(sub).copied().unwrap_or(0);
        assert!(n >= 5, "{sub} has only {n} options, which cannot be right");
    }
}

#[test]
fn callpeak_exposes_the_options_the_corpus_exercises() {
    // A sanity check that the matrix is the real thing and not an empty shell:
    // these flags are all used by oracle/run_oracle.py's variant matrix.
    let text = std::fs::read_to_string(matrix_path()).expect("readable");
    let callpeak: BTreeSet<&str> = text
        .lines()
        .skip(1)
        .filter_map(|line| {
            let mut f = line.split('\t');
            match (f.next(), f.next()) {
                (Some("callpeak"), Some(flag)) => Some(flag),
                _ => None,
            }
        })
        .collect();
    for flag in [
        "--treatment",
        "--control",
        "--format",
        "--gsize",
        "--keep-dup",
        "--nomodel",
        "--extsize",
        "--shift",
        "--call-summits",
        "--broad",
        "--broad-cutoff",
        "--nolambda",
        "--slocal",
        "--llocal",
        "--bw",
        "--mfold",
        "--SPMR",
        "--scale-to",
        "--pvalue",
        "--qvalue",
        "--cutoff-analysis",
        "--bdg",
    ] {
        assert!(
            callpeak.contains(flag),
            "callpeak {flag} missing from the matrix"
        );
    }
}
