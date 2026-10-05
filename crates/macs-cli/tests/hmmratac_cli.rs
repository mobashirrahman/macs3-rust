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

/// `*_training_data.txt` is written with CPython's `repr`, not Rust's `Display`.
///
/// `hmmratac_cmd.py:350` interpolates each cell as `f"{v[k]}"`, and a Gaussian
/// `v[k]` is a `numpy.float32`, whose empty format spec widens to `double` and
/// defers to `float.__repr__`. Rust's `Display` is shortest-round-trip as well but
/// rounds an exact halfway digit away from zero where CPython rounds to even, so
/// `89.1022f32` came out as `89.10220336914063` instead of `89.10220336914062`.
/// A Poisson `v[k]` is `int(max(...))` -- a Python `int`, hence a bare integer.
#[test]
fn saved_training_data_is_byte_identical_to_the_reference() {
    const GAUSSIAN: &str = "\
b'chrM'\t10\t0.0001\t283.0\t0.0001\t0.0001
b'chrM'\t20\t0.0001\t398.0\t0.0001\t0.0001
b'chrM'\t30\t0.0001\t510.0\t0.0001\t0.0001
b'chrM'\t40\t0.0001\t631.0\t0.0001\t0.0001
b'chrM'\t50\t0.00021014314552303404\t757.0\t0.0001\t0.0001
b'chrM'\t60\t0.0027049786876887083\t895.9983520507812\t0.0001\t0.0001
b'chrM'\t70\t0.029203087091445923\t1028.971923828125\t0.0001\t0.0001
b'chrM'\t80\t0.2977858781814575\t1162.70166015625\t0.0001\t0.0001
b'chrM'\t90\t2.5513854026794434\t1295.44677734375\t0.0001\t0.0001
b'chrM'\t100\t14.736224174499512\t1391.2623291015625\t0.0001\t0.0001
b'chrM'\t110\t45.056800842285156\t1441.940673828125\t0.0001\t0.0001
b'chrM'\t120\t89.10220336914062\t1451.8955078125\t0.0001\t0.0001
b'chrM'\t130\t115.05833435058594\t1452.9395751953125\t0.0001\t0.0001
b'chrM'\t140\t115.05833435058594\t1452.9395751953125\t0.0001\t0.0001
b'chrM'\t150\t121.04974365234375\t1452.9481201171875\t0.0001\t0.0001
b'chrM'\t160\t136.04336547851562\t1452.954345703125\t0.0001\t0.0001
b'chrM'\t170\t144.04209899902344\t1452.955322265625\t0.0001\t0.0001
b'chrM'\t180\t144.04209899902344\t1452.955322265625\t0.0001\t0.0001
b'chrM'\t190\t144.04209899902344\t1169.955322265625\t0.0001\t0.0001
b'chrM'\t200\t144.04209899902344\t1054.955322265625\t0.0001\t0.0001
";
    const POISSON: &str = "\
b'chrM'\t10\t0\t283\t0\t0
b'chrM'\t20\t0\t398\t0\t0
b'chrM'\t30\t0\t510\t0\t0
b'chrM'\t40\t0\t631\t0\t0
b'chrM'\t50\t0\t757\t0\t0
b'chrM'\t60\t0\t895\t0\t0
b'chrM'\t70\t0\t1028\t0\t0
b'chrM'\t80\t0\t1162\t0\t0
b'chrM'\t90\t2\t1295\t0\t0
b'chrM'\t100\t14\t1391\t0\t0
b'chrM'\t110\t45\t1441\t0\t0
b'chrM'\t120\t89\t1451\t0\t0
b'chrM'\t130\t115\t1452\t0\t0
b'chrM'\t140\t115\t1452\t0\t0
b'chrM'\t150\t121\t1452\t0\t0
b'chrM'\t160\t136\t1452\t0\t0
b'chrM'\t170\t144\t1452\t0\t0
b'chrM'\t180\t144\t1452\t0\t0
b'chrM'\t190\t144\t1169\t0\t0
b'chrM'\t200\t144\t1054\t0\t0
";

    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/sweep/gmini_mfrag_d400_w180_ctrl_033/treat.frag");
    for (hmm_type, want) in [("gaussian", GAUSSIAN), ("poisson", POISSON)] {
        let out = std::env::temp_dir().join(format!(
            "macs3rs-hmmratac-repr-{hmm_type}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&out);
        std::fs::create_dir_all(&out).unwrap();
        let args = vec![
            s("-i"),
            s(fixture.to_str().unwrap()),
            s("-f"),
            s("FRAG"),
            s("--hmm-type"),
            s(hmm_type),
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
            s("t"),
            s("--outdir"),
            s(out.to_str().unwrap()),
        ];
        let options: Options = parse_flags("hmmratac", &args).expect("valid options");
        macs_cli::commands::hmmratac::hmmratac(&options).expect("self-training succeeds");
        let got = std::fs::read_to_string(out.join("t_training_data.txt")).unwrap();
        assert_eq!(got, want, "{hmm_type} training data");
        let _ = std::fs::remove_dir_all(out);
    }
}
