//! Self-training CLI regression coverage. The checked-in fragment fixture keeps
//! this focused on training selection, model serialization, and --modelonly.

use macs_cli::{parse_flags, Options};

fn s(value: &str) -> String {
    value.to_string()
}

#[test]
fn modelonly_self_trains_and_saves_training_artifacts() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/pe_basic/gauss_fragments/treat.bedpe");
    let out = std::env::temp_dir().join(format!("macs3rs-hmmratac-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        s(root.to_str().unwrap()),
        s("-f"),
        s("BEDPE"),
        s("--no-fragem"),
        s("--modelonly"),
        s("--save-training-data"),
        s("--lower"),
        s("0"),
        s("--upper"),
        s("100000"),
        s("--training-flanking"),
        s("20"),
        s("--binsize"),
        s("10"),
        s("--name"),
        s("selftrain"),
        s("--outdir"),
        s(out.to_str().unwrap()),
    ];
    let options: Options = parse_flags("hmmratac", &args).expect("valid options");
    macs_cli::commands::hmmratac::hmmratac(&options).expect("self-training succeeds");

    for filename in [
        "selftrain_model.json",
        "selftrain_training_regions.bed",
        "selftrain_training_data.txt",
        "selftrain_training_lengths.txt",
    ] {
        let path = out.join(filename);
        assert!(path.is_file(), "missing training artifact {filename}");
        assert!(
            !std::fs::read(path).unwrap().is_empty(),
            "empty training artifact {filename}"
        );
    }
    assert!(
        !out.join("selftrain_accessible_regions.narrowPeak").exists(),
        "--modelonly should return before decoding"
    );
    let model = std::fs::read_to_string(out.join("selftrain_model.json")).unwrap();
    assert!(model.contains("\"hmm_type\": \"gaussian\""));
    assert!(model.contains("\"hmm_binsize\": 10"));
    assert!(model.contains("\"n_features\": 4"));
    let _ = std::fs::remove_dir_all(out);
}

#[test]
fn poisson_training_accepts_a_supplied_bed_without_rewriting_it() {
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/pe_basic/gauss_fragments/treat.bedpe");
    let out = std::env::temp_dir().join(format!("macs3rs-hmmratac-poisson-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    let bed = out.join("training.bed");
    std::fs::write(&bed, "chrM\t0\t220\n").unwrap();
    let args = vec![
        s("-i"),
        s(fixture.to_str().unwrap()),
        s("-f"),
        s("BEDPE"),
        s("--no-fragem"),
        s("--hmm-type"),
        s("poisson"),
        s("--modelonly"),
        s("--save-training-data"),
        s("--training"),
        s(bed.to_str().unwrap()),
        s("--binsize"),
        s("10"),
        s("--name"),
        s("poisson"),
        s("--outdir"),
        s(out.to_str().unwrap()),
    ];
    let options: Options = parse_flags("hmmratac", &args).expect("valid options");
    macs_cli::commands::hmmratac::hmmratac(&options).expect("Poisson self-training succeeds");
    let model = std::fs::read_to_string(out.join("poisson_model.json")).unwrap();
    assert!(model.contains("\"hmm_type\": \"poisson\""));
    assert!(out.join("poisson_training_data.txt").is_file());
    assert!(!out.join("poisson_training_regions.bed").exists());
    let _ = std::fs::remove_dir_all(out);
}

#[test]
fn empty_training_bed_returns_an_error_without_panicking() {
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/pe_basic/gauss_fragments/treat.bedpe");
    let out = std::env::temp_dir().join(format!("macs3rs-hmmratac-empty-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    let bed = out.join("empty.bed");
    std::fs::write(&bed, "# no intervals\n").unwrap();
    let args = vec![
        s("-i"),
        s(fixture.to_str().unwrap()),
        s("-f"),
        s("BEDPE"),
        s("--no-fragem"),
        s("--modelonly"),
        s("--training"),
        s(bed.to_str().unwrap()),
        s("--name"),
        s("empty"),
        s("--outdir"),
        s(out.to_str().unwrap()),
    ];
    let options: Options = parse_flags("hmmratac", &args).expect("valid options");
    assert!(macs_cli::commands::hmmratac::hmmratac(&options).is_err());
    assert!(!out.join("empty_model.json").exists());
    let _ = std::fs::remove_dir_all(out);
}
