//! Regressions for commands whose `-i` option accepts multiple treatment files.

use macs_cli::{commands, parse_flags, Options};
use std::path::PathBuf;

fn tempdir(tag: &str) -> PathBuf {
    let path =
        std::env::temp_dir().join(format!("macs3rs-multi-input-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&path);
    std::fs::create_dir_all(&path).unwrap();
    path
}

fn write_tags(path: &std::path::Path, chrom: &str, base: u64, count: usize) {
    let mut text = String::new();
    for i in 0..count {
        let start = base + i as u64 * 10;
        text.push_str(&format!("{chrom}\t{start}\t{}\tr{i}\t0\t+\n", start + 5));
    }
    std::fs::write(path, text).unwrap();
}

fn parse(command: &str, args: &[String]) -> Options {
    parse_flags(command, args).unwrap_or_else(|e| panic!("{command} args parse: {e}"))
}

fn write_fragments(path: &std::path::Path, base: u64, count: usize) {
    let mut text = String::new();
    for i in 0..count {
        let start = base + i as u64 * 100;
        text.push_str(&format!("chr1\t{start}\t{}\n", start + 50));
    }
    std::fs::write(path, text).unwrap();
}

#[test]
fn filterdup_and_randsample_pool_every_treatment_file() {
    let dir = tempdir("se");
    let first = dir.join("first.bed");
    let second = dir.join("second.bed");
    write_tags(&first, "chr1", 100, 10);
    write_tags(&second, "chr1", 1000, 10);
    let first = first.to_string_lossy().into_owned();
    let second = second.to_string_lossy().into_owned();

    let filter_out = dir.join("filterdup.bed");
    let filter_args = vec![
        "-i".into(),
        first.clone(),
        second.clone(),
        "-f".into(),
        "BED".into(),
        "--keep-dup".into(),
        "all".into(),
        "-o".into(),
        filter_out.to_string_lossy().into_owned(),
    ];
    commands::filterdup::filterdup(&parse("filterdup", &filter_args)).unwrap();
    let filtered = std::fs::read_to_string(&filter_out).unwrap();
    assert_eq!(
        filtered.lines().count(),
        20,
        "both input files are retained"
    );
    assert!(filtered.contains("chr1\t1000\t1005"));

    let sample_out = dir.join("sample.bed");
    let sample_args = vec![
        "-i".into(),
        first,
        second,
        "-f".into(),
        "BED".into(),
        "-n".into(),
        "10".into(),
        "--seed".into(),
        "17".into(),
        "-o".into(),
        sample_out.to_string_lossy().into_owned(),
    ];
    commands::randsample::randsample(&parse("randsample", &sample_args)).unwrap();
    let sampled = std::fs::read_to_string(&sample_out).unwrap();
    assert_eq!(
        sampled,
        concat!(
            "chr1\t110\t115\t.\t.\t+\n",
            "chr1\t120\t125\t.\t.\t+\n",
            "chr1\t130\t135\t.\t.\t+\n",
            "chr1\t150\t155\t.\t.\t+\n",
            "chr1\t180\t185\t.\t.\t+\n",
            "chr1\t1000\t1005\t.\t.\t+\n",
            "chr1\t1010\t1015\t.\t.\t+\n",
            "chr1\t1020\t1025\t.\t.\t+\n",
            "chr1\t1030\t1035\t.\t.\t+\n",
            "chr1\t1040\t1045\t.\t.\t+\n",
        ),
        "seeded sampling must include both files and match the upstream stream",
    );

    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn explicit_tsize_skips_bowtie_inference_for_short_inputs() {
    let dir = tempdir("bowtie-explicit-tsize");
    let input = dir.join("short.bowtie");
    std::fs::write(&input, "r1\t+\tchr1\t10\tACGT\tIIII\t0\t\n").unwrap();
    let output = dir.join("filtered.bed");
    let args = vec![
        "-i".into(),
        input.to_string_lossy().into_owned(),
        "-f".into(),
        "BOWTIE".into(),
        "-s".into(),
        "4".into(),
        "--keep-dup".into(),
        "all".into(),
        "-o".into(),
        output.to_string_lossy().into_owned(),
    ];
    commands::filterdup::filterdup(&parse("filterdup", &args))
        .expect("explicit --tsize avoids BowtieParser.tsize() on short input");
    assert_eq!(
        std::fs::read_to_string(&output).unwrap(),
        "chr1\t10\t14\t.\t.\t+\n"
    );
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn paired_filterdup_and_randsample_pool_every_treatment_file() {
    let dir = tempdir("pe");
    let first = dir.join("first.bedpe");
    let second = dir.join("second.bedpe");
    write_fragments(&first, 100, 4);
    write_fragments(&second, 1000, 4);
    let first = first.to_string_lossy().into_owned();
    let second = second.to_string_lossy().into_owned();

    let filter_out = dir.join("filterdup.bedpe");
    let filter_args = vec![
        "-i".into(),
        first.clone(),
        second.clone(),
        "-f".into(),
        "BEDPE".into(),
        "--keep-dup".into(),
        "all".into(),
        "-o".into(),
        filter_out.to_string_lossy().into_owned(),
    ];
    commands::filterdup::filterdup(&parse("filterdup", &filter_args)).unwrap();
    let filtered = std::fs::read_to_string(&filter_out).unwrap();
    assert_eq!(
        filtered,
        concat!(
            "chr1\t100\t150\n",
            "chr1\t200\t250\n",
            "chr1\t300\t350\n",
            "chr1\t400\t450\n",
            "chr1\t1000\t1050\n",
            "chr1\t1100\t1150\n",
            "chr1\t1200\t1250\n",
            "chr1\t1300\t1350\n",
        ),
        "filterdup must append all fragments from both files",
    );

    let sample_out = dir.join("sample.bedpe");
    let sample_args = vec![
        "-i".into(),
        first,
        second,
        "-f".into(),
        "BEDPE".into(),
        "-n".into(),
        "4".into(),
        "--seed".into(),
        "17".into(),
        "-o".into(),
        sample_out.to_string_lossy().into_owned(),
    ];
    commands::randsample::randsample(&parse("randsample", &sample_args)).unwrap();
    let sampled = std::fs::read_to_string(&sample_out).unwrap();
    assert_eq!(
        sampled,
        concat!(
            "chr1\t300\t350\n",
            "chr1\t400\t450\n",
            "chr1\t1000\t1050\n",
            "chr1\t1300\t1350\n",
        ),
        "seeded paired sampling matches the pinned oracle over both inputs",
    );
    let _ = std::fs::remove_dir_all(dir);
}
