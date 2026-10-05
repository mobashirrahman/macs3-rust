//! End-to-end tests for `bdgpeakcall` / `bdgbroadcall`.
//!
//! Byte expectations come from upstream MACS3 3.0.5 on the same bedGraph,
//! covering the summit selection, the track line, the `%.1f` default filename and
//! the `--no-trackline` negation.

use macs_cli::{parse_flags, Options};

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-bpc-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn s(x: &str) -> String {
    x.to_string()
}

fn fixture(dir: &std::path::Path) -> String {
    let p = dir.join("in.bdg");
    std::fs::write(
        &p,
        "chr1\t0\t100\t2\n\
         chr1\t100\t200\t10\n\
         chr1\t200\t300\t12\n\
         chr1\t300\t400\t1\n\
         chr1\t400\t600\t8\n\
         chr1\t600\t700\t3\n",
    )
    .unwrap();
    p.to_string_lossy().into_owned()
}

#[test]
fn bdgpeakcall_matches_upstream_bytes() {
    let dir = tmpdir("main");
    let input = fixture(&dir);
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        s(&input),
        s("-c"),
        s("5"),
        s("-l"),
        s("100"),
        s("-g"),
        s("30"),
        s("--outdir"),
        s(out.to_str().unwrap()),
        s("--o-prefix"),
        s("t"),
    ];
    let o: Options = parse_flags("bdgpeakcall", &args).expect("parse");
    macs_cli::commands::bdgpeakcall::bdgpeakcall(&o).expect("bdgpeakcall runs");
    let got = std::fs::read_to_string(out.join("t_c5.0_l100_g30_peaks.narrowPeak")).unwrap();
    // upstream bytes: track line uses `name` (t), peaks use peakprefix (t_narrowPeak)
    let want = "track type=narrowPeak name=\"t\" description=\"t\" nextItemButton=on\n\
chr1\t100\t300\tt_narrowPeak1\t120\t.\t0\t0\t0\t150\n\
chr1\t400\t600\tt_narrowPeak2\t80\t.\t0\t0\t0\t100\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn no_trackline_suppresses_only_the_track_line() {
    let dir = tmpdir("notrack");
    let input = fixture(&dir);
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        s(&input),
        s("-c"),
        s("5"),
        s("-l"),
        s("100"),
        s("-g"),
        s("30"),
        s("--no-trackline"),
        s("--outdir"),
        s(out.to_str().unwrap()),
        s("--o-prefix"),
        s("t"),
    ];
    let o: Options = parse_flags("bdgpeakcall", &args).expect("parse");
    macs_cli::commands::bdgpeakcall::bdgpeakcall(&o).expect("bdgpeakcall runs");
    let got = std::fs::read_to_string(out.join("t_c5.0_l100_g30_peaks.narrowPeak")).unwrap();
    assert!(
        !got.starts_with("track "),
        "--no-trackline must drop the track line"
    );
    assert!(got.starts_with("chr1\t"), "but keep the data rows:\n{got}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn the_default_filename_formats_the_cutoff_with_one_decimal() {
    let dir = tmpdir("fname");
    let input = fixture(&dir);
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        s(&input),
        s("-c"),
        s("5"),
        s("--outdir"),
        s(out.to_str().unwrap()),
        s("--o-prefix"),
        s("t"),
    ];
    let o: Options = parse_flags("bdgpeakcall", &args).expect("parse");
    macs_cli::commands::bdgpeakcall::bdgpeakcall(&o).expect("runs");
    // upstream default name is %s_c%.1f_l%d_g%d_peaks.narrowPeak
    let name = std::fs::read_dir(&out)
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .file_name();
    let name = name.to_string_lossy();
    assert!(name.contains("_c5.0_"), "cutoff formatted %.1f: {name}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn bdgbroadcall_matches_upstream_gappedpeak_bytes() {
    let dir = tmpdir("broad");
    let input = fixture(&dir);
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        s(&input),
        s("-c"),
        s("5"),
        s("-C"),
        s("2"),
        s("-l"),
        s("100"),
        s("-g"),
        s("30"),
        s("-G"),
        s("800"),
        s("--outdir"),
        s(out.to_str().unwrap()),
        s("--o-prefix"),
        s("t"),
    ];
    let o: Options = parse_flags("bdgbroadcall", &args).expect("parse");
    macs_cli::commands::bdgpeakcall::bdgbroadcall(&o).expect("bdgbroadcall runs");
    let got = std::fs::read_to_string(out.join("t_c5.0_C2.00_l100_g30_G800_broad.bed12")).unwrap();
    // upstream: track line name/description default to "peak"; one gapped broad
    // region 0..700 with 4 blocks (the two 1bp complements + the two lvl1 peaks).
    let want = "track name=\"peak\" description=\"peak\" type=gappedPeak nextItemButton=on\n\
chr1\t0\t700\tt_broadRegion1\t120\t.\t0\t0\t0\t4\t1,200,200,1\t0,100,400,699\t0\t0\t0\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}

/// Peak numbers run continuously across chromosomes, matching upstream's
/// `write_to_narrowPeak` (`PeakIO.py:719-735`), where `n_peak` is initialised
/// once before the chromosome loop.
///
/// Bytes below are MACS3 3.0.5's for `bdgpeakcall -c 5 -l 50 -g 30
/// --o-prefix t` on this two-chromosome bedGraph: the chr2 peak is
/// `t_narrowPeak2`, not `t_narrowPeak1`.
#[test]
fn peak_numbers_run_continuously_across_chromosomes() {
    let dir = tmpdir("continuous");
    let input = dir.join("in.bdg");
    std::fs::write(
        &input,
        "chr1\t0\t100\t10\n\
         chr1\t100\t300\t1\n\
         chr2\t0\t150\t10\n\
         chr2\t150\t400\t1\n",
    )
    .unwrap();
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        s(input.to_str().unwrap()),
        s("-c"),
        s("5"),
        s("-l"),
        s("50"),
        s("-g"),
        s("30"),
        s("--outdir"),
        s(out.to_str().unwrap()),
        s("--o-prefix"),
        s("t"),
    ];
    let o: Options = parse_flags("bdgpeakcall", &args).expect("parse");
    macs_cli::commands::bdgpeakcall::bdgpeakcall(&o).expect("bdgpeakcall runs");
    let got = std::fs::read_to_string(out.join("t_c5.0_l50_g30_peaks.narrowPeak")).unwrap();
    let want = "track type=narrowPeak name=\"t\" description=\"t\" nextItemButton=on\n\
chr1\t0\t100\tt_narrowPeak1\t100\t.\t0\t0\t0\t50\n\
chr2\t0\t150\tt_narrowPeak2\t100\t.\t0\t0\t0\t75\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}

/// `call_broadpeaks` (`BedGraph.py:564-567`) iterates only the chromosomes of
/// the lvl1 peaks, so a chromosome that only reaches the linking cutoff
/// contributes no broad peak. Here chr2 peaks at 3, above the linking cutoff
/// 2 but below the peak cutoff 5.
///
/// Bytes below are MACS3 3.0.5's for `bdgbroadcall -c 5 -C 2 -l 50 -g 30 -G
/// 800 --o-prefix t` on this bedGraph.
#[test]
fn bdgbroadcall_drops_chromosomes_without_strong_peaks() {
    let dir = tmpdir("broad-drop");
    let input = dir.join("in.bdg");
    std::fs::write(
        &input,
        "chr1\t0\t100\t10\n\
         chr1\t100\t300\t1\n\
         chr2\t0\t150\t3\n\
         chr2\t150\t400\t1\n",
    )
    .unwrap();
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        s(input.to_str().unwrap()),
        s("-c"),
        s("5"),
        s("-C"),
        s("2"),
        s("-l"),
        s("50"),
        s("-g"),
        s("30"),
        s("-G"),
        s("800"),
        s("--outdir"),
        s(out.to_str().unwrap()),
        s("--o-prefix"),
        s("t"),
    ];
    let o: Options = parse_flags("bdgbroadcall", &args).expect("parse");
    macs_cli::commands::bdgpeakcall::bdgbroadcall(&o).expect("bdgbroadcall runs");
    let got = std::fs::read_to_string(out.join("t_c5.0_C2.00_l50_g30_G800_broad.bed12")).unwrap();
    let want = "track name=\"peak\" description=\"peak\" type=gappedPeak nextItemButton=on\n\
chr1\t0\t100\tt_broadRegion1\t100\t.\t0\t0\t0\t1\t100\t0\t0\t0\t0\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}

/// With `-o`, upstream renames the o-prefix to the ofile name
/// (`bdgbroadcall_cmd.py:53-55`), so the peak names use the ofile too.
/// Bytes below are MACS3 3.0.5's for the same bedGraph as above with `-o
/// custom.bed12`.
#[test]
fn bdgbroadcall_names_peaks_after_ofile_when_o_given() {
    let dir = tmpdir("broad-ofile");
    let input = dir.join("in.bdg");
    std::fs::write(
        &input,
        "chr1\t0\t100\t10\n\
         chr1\t100\t300\t1\n\
         chr2\t0\t150\t3\n\
         chr2\t150\t400\t1\n",
    )
    .unwrap();
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        s(input.to_str().unwrap()),
        s("-c"),
        s("5"),
        s("-C"),
        s("2"),
        s("-l"),
        s("50"),
        s("-g"),
        s("30"),
        s("-G"),
        s("800"),
        s("--outdir"),
        s(out.to_str().unwrap()),
        s("-o"),
        s("custom.bed12"),
    ];
    let o: Options = parse_flags("bdgbroadcall", &args).expect("parse");
    macs_cli::commands::bdgpeakcall::bdgbroadcall(&o).expect("bdgbroadcall runs");
    let got = std::fs::read_to_string(out.join("custom.bed12")).unwrap();
    let want = "track name=\"peak\" description=\"peak\" type=gappedPeak nextItemButton=on\n\
chr1\t0\t100\tcustom.bed12_broadRegion1\t100\t.\t0\t0\t0\t1\t100\t0\t0\t0\t0\n";
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&dir);
}
