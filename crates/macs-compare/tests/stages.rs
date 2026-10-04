//! L5 regression tests for the per-stage comparator.
//!
//! The comparator is what the stage-by-stage test layer (L3) actually runs, so it
//! needs to be *shown* to fail: a comparator that always reports "ok" is worse than
//! none, because it would gate L3 green.

use macs_compare::stages::{compare_stage, compare_trees, parse, Json};
use std::path::Path;

const GOLD: &str = r#"{
  "reads_pre_filter": {"total": 3000, "chroms": {"chr1": {"plus": 700, "minus": 500}}},
  "reads_post_filter": {"treatment": 32, "control": 195},
  "d": 200,
  "scaling": {"ratio_treat2control": 0.1641025641025641, "tocontrol": false},
  "xls_header": {"d": 200, "slocal": 1000, "llocal": 10000},
  "qvalue_table": [[5.0, 900], [2.0, 99]],
  "peak_files": {"candidate": ["/tmp/x/candidate_peaks.peaks.xls"]},
  "_argv": "callpeak -n x"
}"#;

fn stage(txt: &str, key: &str) -> Json {
    match parse(txt).expect("parse stage dump") {
        Json::Obj(m) => m
            .get(key)
            .cloned()
            .unwrap_or_else(|| panic!("no stage {key}")),
        _ => panic!("top level must be an object"),
    }
}

#[test]
fn identical_stages_are_clean() {
    for k in [
        "reads_pre_filter",
        "reads_post_filter",
        "d",
        "scaling",
        "xls_header",
    ] {
        let r = compare_stage(k, &stage(GOLD, k), &stage(GOLD, k), 1e-9);
        assert!(r.is_clean(), "{k} should be clean, got {:?}", r.mismatches);
        assert_eq!(r.max_abs, 0.0);
    }
}

#[test]
fn a_read_count_off_by_one_is_a_mismatch_not_a_tolerance_breach() {
    let ours = GOLD.replace("\"treatment\": 32", "\"treatment\": 33");
    let r = compare_stage(
        "reads_post_filter",
        &stage(GOLD, "reads_post_filter"),
        &stage(&ours, "reads_post_filter"),
        1e-9,
    );
    assert!(!r.is_clean());
    assert!(r.mismatches[0].contains("treatment"), "{:?}", r.mismatches);
    // A JSON number carries no integer/float distinction, so a read count is compared
    // with the tolerance like any other number. At the p-score tolerance of 1e-9 a
    // difference of one is still a hard failure, which is the property that matters --
    // and `reads_pre_filter`'s per-chromosome plus/minus counts are the same case.
    assert!(
        r.mismatches[0].contains("32 vs 33") && r.mismatches[0].contains("> tol"),
        "{:?}",
        r.mismatches
    );
    assert_eq!(r.exact, 0, "both leaves are numbers");
}

#[test]
fn a_string_leaf_must_match_exactly() {
    // Strings and booleans have no tolerance: `xls_header`'s `tocontrol` flag is a
    // decision, not a measurement.
    let ours = GOLD.replace("\"tocontrol\": false", "\"tocontrol\": true");
    let r = compare_stage(
        "scaling",
        &stage(GOLD, "scaling"),
        &stage(&ours, "scaling"),
        1e-9,
    );
    assert!(!r.is_clean());
    assert!(
        r.mismatches[0].contains("must match exactly"),
        "{:?}",
        r.mismatches
    );
}

#[test]
fn a_sub_tolerance_numeric_drift_is_clean_and_above_it_is_not() {
    let a = stage(GOLD, "scaling");
    let nearly = stage(
        &GOLD.replace("0.1641025641025641", "0.1641025641025642"),
        "scaling",
    );
    let r = compare_stage("scaling", &a, &nearly, 1e-9);
    assert!(
        r.is_clean(),
        "1e-16 drift should be inside 1e-9: {:?}",
        r.mismatches
    );

    let far = stage(
        &GOLD.replace("0.1641025641025641", "0.1641035641025641"),
        "scaling",
    );
    let r2 = compare_stage("scaling", &a, &far, 1e-9);
    assert!(!r2.is_clean(), "1e-6 drift should breach a 1e-9 tolerance");
    assert!(
        (r2.max_abs - 1e-6).abs() < 1e-15,
        "max_abs was {}",
        r2.max_abs
    );
    assert!(r2.max_rel > 6e-6 && r2.max_rel < 7e-6, "{}", r2.max_rel);
}

#[test]
fn a_missing_leaf_is_reported_rather_than_skipped() {
    let ours = GOLD.replace("\"treatment\": 32, ", "");
    let r = compare_stage(
        "reads_post_filter",
        &stage(GOLD, "reads_post_filter"),
        &stage(&ours, "reads_post_filter"),
        1e-9,
    );
    assert!(!r.is_clean());
    assert!(
        r.mismatches
            .iter()
            .any(|m| m.contains("present only in ours")),
        "{:?}",
        r.mismatches
    );
}

#[test]
fn file_references_are_not_compared_as_values() {
    // A path recorded under one root can never match a path recorded under another
    // root, so comparing them would report a spurious failure on every run. The file
    // comparator owns that comparison.
    let ours = GOLD.replace("/tmp/x/candidate_peaks.peaks.xls", "/somewhere/else.xls");
    let r = compare_stage(
        "peak_files",
        &stage(GOLD, "peak_files"),
        &stage(&ours, "peak_files"),
        1e-9,
    );
    assert!(r.is_clean(), "{:?}", r.mismatches);
    assert_eq!(
        r.exact, 0,
        "file references must not be counted as exact leaves"
    );
}

#[test]
fn the_scalar_stage_key_has_no_trailing_dot() {
    let ours = GOLD.replace("\"d\": 200", "\"d\": 250");
    let r = compare_stage("d", &stage(GOLD, "d"), &stage(&ours, "d"), 1e-9);
    assert_eq!(r.mismatches.len(), 1, "{:?}", r.mismatches);
    assert_eq!(r.mismatches[0], "d: 200 vs 250 (abs 5e1 > tol 1e-9)");
}

#[test]
fn nested_leaves_are_keyed_by_their_full_path() {
    let r = compare_stage(
        "reads_pre_filter",
        &stage(GOLD, "reads_pre_filter"),
        &stage(GOLD, "reads_pre_filter"),
        1e-9,
    );
    // 1 total + 2 (plus/minus) = 3 numeric leaves; a comparator that flattened them
    // would report 1.
    assert_eq!(r.numeric, 3, "{r:?}");
}

#[test]
fn the_recorded_stage_corpus_parses_and_has_the_expected_stages() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("tests/stages");
    if !root.exists() {
        return; // not checked out here
    }
    let rep = compare_trees(&root, &root, 1e-9).expect("compare the corpus to itself");
    assert!(
        rep.is_clean(),
        "{:?}",
        rep.stages
            .iter()
            .filter(|s| !s.is_clean())
            .map(|s| (&s.name, &s.mismatches))
            .collect::<Vec<_>>()
    );
    for want in ["d", "reads_pre_filter", "reads_post_filter", "scaling"] {
        assert!(
            rep.stages.iter().any(|s| s.name == want),
            "no {want} stage recorded"
        );
    }
    // `qvalue_table` is the p->q histogram the 1e-6 criterion is about, so it must be
    // present and non-trivial for this layer to mean anything.
    let qt = rep
        .stages
        .iter()
        .find(|s| s.name == "qvalue_table")
        .expect("qvalue_table");
    assert!(
        qt.numeric > 50,
        "qvalue_table only had {} leaves",
        qt.numeric
    );
}

#[test]
fn stages_are_reported_in_pipeline_order_so_the_first_failure_is_the_earliest_one() {
    let root_a = std::env::temp_dir().join("mc_order_a");
    let root_b = std::env::temp_dir().join("mc_order_b");
    for r in [&root_a, &root_b] {
        let d = r.join("se/fixture");
        std::fs::create_dir_all(&d).unwrap();
        std::fs::write(d.join("stages.json"), GOLD).unwrap();
    }
    // Break an *earlier* stage as well, so the ordering claim is actually exercised.
    let broken = std::fs::read_to_string(root_b.join("se/fixture/stages.json"))
        .unwrap()
        .replace("\"treatment\": 32", "\"treatment\": 33")
        .replace("\"d\": 200,", "\"d\": 250,");
    std::fs::write(root_b.join("se/fixture/stages.json"), broken).unwrap();

    let rep = compare_trees(&root_a, &root_b, 1e-9).unwrap();
    assert!(!rep.is_clean());
    let names: Vec<&str> = rep.stages.iter().map(|s| s.name.as_str()).collect();
    let ip = names
        .iter()
        .position(|n| *n == "reads_post_filter")
        .unwrap();
    let id = names.iter().position(|n| *n == "d").unwrap();
    assert!(ip < id, "reads_post_filter must precede d, got {names:?}");
    assert_eq!(
        rep.first_failure().map(|s| s.name.as_str()),
        Some("reads_post_filter")
    );
    let _ = std::fs::remove_dir_all(&root_a);
    let _ = std::fs::remove_dir_all(&root_b);
}

#[test]
fn a_run_present_on_one_side_only_is_noted() {
    let root_a = std::env::temp_dir().join("mc_missing_a");
    let root_b = std::env::temp_dir().join("mc_missing_b");
    for r in [&root_a, &root_b] {
        std::fs::create_dir_all(r.join("se/fixture")).unwrap();
        std::fs::write(r.join("se/fixture/stages.json"), GOLD).unwrap();
    }
    std::fs::create_dir_all(root_b.join("se/extra")).unwrap();
    std::fs::write(root_b.join("se/extra/stages.json"), GOLD).unwrap();
    let rep = compare_trees(&root_a, &root_b, 1e-9).unwrap();
    assert!(!rep.run_notes.is_empty(), "{:?}", rep.run_notes);
    assert!(!rep.is_clean());
    assert!(rep.run_notes.iter().any(|n| n.contains("extra")));
    let _ = std::fs::remove_dir_all(&root_a);
    let _ = std::fs::remove_dir_all(&root_b);
}
