//! Localisation probe for `predictd`'s paired-peak search.
//!
//! Prints the quantities upstream reports under `--verbose 3` -- total tags,
//! `min_tags`/`max_tags`, and the per-chromosome plus/minus/paired counts -- so a
//! divergence can be attributed to the strand pileup, the per-strand summit
//! filter, or the pairing step, rather than surfacing as a bare
//! "can only find 0".
//!
//! Usage: probe_predictd <file> [format] [gsize] [lmfold] [umfold] [bw]

use macs_cli::commands::input::{load_single_end_files, single_end_format};
use macs_model::probe::find_paired_peaks_pub;
use std::path::Path;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let path = a
        .first()
        .expect("usage: probe_predictd <file> [format] ...");
    let format = a.get(1).cloned().unwrap_or_else(|| "BED".into());
    let gsize: f64 = a
        .get(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(2_913_022_398.0);
    let lmfold: f64 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(100.0);
    let umfold: f64 = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(200.0);
    let bw: i64 = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(300);

    let fmt = single_end_format(Path::new(path), &format).expect("format");
    let (track, tsize) = load_single_end_files(std::slice::from_ref(path), &fmt).expect("load");

    let peaksize = 2 * bw;
    let total = track.total() as f64;
    let min_tags = (total * lmfold * peaksize as f64 / gsize / 2.0).round();
    let max_tags = (total * umfold * peaksize as f64 / gsize / 2.0).round();
    println!(
        "tsize={tsize} total={total} peaksize={peaksize} min_tags={min_tags} max_tags={max_tags}"
    );

    let (_paired, counts) = find_paired_peaks_pub(&track, peaksize, min_tags, max_tags);
    let tp: usize = counts.iter().map(|c| c.1).sum();
    let tm: usize = counts.iter().map(|c| c.2).sum();
    let tt: usize = counts.iter().map(|c| c.3).sum();
    println!(
        "chroms={} plus_summits={tp} minus_summits={tm} paired_centres={tt}",
        counts.len()
    );
    for (n, p, m, pr) in counts.iter().take(8) {
        println!(
            "  {} plus={p} minus={m} paired={pr}",
            String::from_utf8_lossy(n)
        );
    }
}
