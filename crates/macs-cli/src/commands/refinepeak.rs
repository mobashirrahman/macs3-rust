//! `macs3-rs refinepeak`: refine peak summits using Watson/Crick tag-depth.
//!
//! Port of `refinepeak_cmd.run` + `find_summit`
//! (`refinepeak_cmd.py:29-102`). For each input peak, tags within
//! `[peak.start - window, peak.end + window]` are collected per strand and a
//! running **Watson-Crick Tag Depth** statistic
//!
//! ```text
//! wtd(j) = 2 * sqrt(watson_left(j) * crick_right(j)) - watson_right(j) - crick_left(j)
//! ```
//!
//! is evaluated at every position; the argmax is the refined summit. The peak is
//! tagged `_R` (refined) when the best value exceeds `--cutoff`, else `_F`.

use std::collections::HashMap;
use std::path::Path;

use macs_core::{MacsError, Result, Strand};
use macs_track::SingleEndTrack;

use crate::Options;

fn load_bed(path: &Path) -> Result<SingleEndTrack> {
    super::input::load_single_end_bed(path)
}

/// A peak from the `--bedfile`: `chrom start end name`.
#[derive(Debug, Clone, PartialEq)]
struct BedPeak {
    chrom: String,
    start: i64,
    end: i64,
    name: String,
}

fn read_peaks(path: &Path) -> Result<Vec<BedPeak>> {
    let text = std::fs::read_to_string(path).map_err(MacsError::Io)?;
    let mut out = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 4 {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[1].parse::<i64>(), f[2].parse::<i64>()) else {
            continue;
        };
        out.push(BedPeak {
            chrom: f[0].to_string(),
            start: s,
            end: e,
            name: f[3].to_string(),
        });
    }
    // upstream sorts the peaks before processing
    out.sort_by(|a, b| (a.chrom.as_str(), a.start, a.end).cmp(&(b.chrom.as_str(), b.start, b.end)));
    Ok(out)
}

fn counter(pos: &[u64], lo: i64, hi: i64) -> HashMap<i64, i64> {
    let mut m: HashMap<i64, i64> = HashMap::new();
    for &p in pos {
        let p = p as i64;
        if p >= lo && p <= hi {
            *m.entry(p).or_insert(0) += 1;
        }
    }
    m
}

/// `find_summit`: running WTD over the window, returning `(best_pos, best_val)`.
fn find_summit(plus: &[u64], minus: &[u64], peak_start: i64, peak_end: i64, w: i64) -> (i64, f64) {
    // The counters hold tags in `[peak_start, peak_end]` only -- the same window
    // upstream's `compute_region_tags_from_peaks` collects into `rt_plus`/`rt_minus`.
    // The `w`-expansion happens inside `sum_le`/`sum_ge`, not here. Counting
    // `[peak_start-w, peak_end+w]` instead pulled in tags outside the window, which
    // inflated the edge-bin sums and moved the WTD maximum (2.83 vs 0.00 on a peak
    // with no plus-strand tags in range).
    let watson = counter(plus, peak_start, peak_end);
    let crick = counter(minus, peak_start, peak_end);
    let sum_le = |c: &HashMap<i64, i64>, pos: i64| -> i64 {
        c.iter()
            .filter(|&(x, _)| *x <= pos && *x >= pos - w)
            .map(|(_, &v)| v)
            .sum()
    };
    let sum_ge = |c: &HashMap<i64, i64>, pos: i64| -> i64 {
        c.iter()
            .filter(|&(x, _)| *x >= pos && *x <= pos + w)
            .map(|(_, &v)| v)
            .sum()
    };
    let get = |c: &HashMap<i64, i64>, pos: i64| -> i64 { *c.get(&pos).unwrap_or(&0) };

    let mut wl = sum_le(&watson, peak_start);
    let mut cl = sum_le(&crick, peak_start);
    let mut wr = sum_ge(&watson, peak_start);
    let mut cr = sum_ge(&crick, peak_start);

    let mut best_val = f64::MIN;
    let mut best_pos = peak_start;
    for j in peak_start..=peak_end {
        let v = 2.0 * ((wl as f64) * (cr as f64)).sqrt() - (wr as f64) - (cl as f64);
        if v > best_val {
            best_val = v;
            best_pos = j;
        }
        // slide the window one base right
        wl += get(&watson, j) - get(&watson, j - w);
        wr += get(&watson, j + w) - get(&watson, j);
        cl += get(&crick, j) - get(&crick, j - w);
        cr += get(&crick, j + w) - get(&crick, j);
    }
    (best_pos, best_val)
}

pub fn refinepeak(o: &Options) -> Result<()> {
    let bedfile = o
        .get("bedfile")
        .ok_or_else(|| MacsError::InvalidParameter("-b/--bedfile is required".into()))?;
    let ifile = o
        .get("ifile")
        .ok_or_else(|| MacsError::InvalidParameter("-i/--ifile is required".into()))?;
    let format = o.get("format").unwrap_or("AUTO").to_uppercase();
    if format == "BAMPE" || format == "BEDPE" {
        return Err(MacsError::InvalidParameter(
            "paired-end input is not yet supported".into(),
        ));
    }
    let track = if format == "BAM" {
        super::input::load_single_end_bam(Path::new(ifile))?.0
    } else if format == "SAM" {
        super::input::load_single_end_sam(Path::new(ifile))?.0
    } else {
        load_bed(Path::new(ifile))?
    };
    let peaks = read_peaks(Path::new(bedfile))?;
    let w = o.int("windowsize").unwrap_or(200);
    let cutoff = o.float("cutoff").unwrap_or(5.0);

    // Group tags by chromosome once.
    //
    // The vectors hold **every** tag on the chromosome, unfiltered. Filtering to a
    // peak's window happens inside `counter`, per peak. Pre-filtering here is wrong
    // for two independent reasons:
    //
    // 1. The window `find_summit` needs is `[pk.start-2w, pk.end+2w]`, not
    //    `[pk.start-w, pk.end+w]`, because it scans an already-expanded window with
    //    an internal half-width-`w` slide. A narrower pre-filter truncated the edge
    //    bins (2 of 730 scores wrong on CTCF data).
    // 2. Accumulating per-peak windows into one per-chromosome vector double-counts
    //    every tag in the overlap of two peaks' windows. Widening the window to fix
    //    (1) made (2) worse instead (different peaks broke). There is no correct
    //    pre-filter width; the filter belongs per peak.
    //
    // Upstream keeps the full arrays (`Counter(plus)` over everything) and filters
    // inside `left_sum`/`right_sum`, so this matches it exactly.
    let mut by_chrom: HashMap<String, (Vec<u64>, Vec<u64>)> = HashMap::new();
    for pk in &peaks {
        let Some(chrom) = track.genome().get(pk.chrom.as_bytes()) else {
            continue;
        };
        by_chrom.entry(pk.chrom.clone()).or_insert_with(|| {
            (
                track.positions().strand(chrom, Strand::Plus).to_vec(),
                track.positions().strand(chrom, Strand::Minus).to_vec(),
            )
        });
    }

    let mut out = String::new();
    for pk in &peaks {
        let (plus, minus) = match by_chrom.get(&pk.chrom) {
            Some(v) => v,
            None => continue,
        };
        // the window over which find_summit runs
        let ps = pk.start - w;
        let pe = pk.end + w;
        let (best_pos, best_val) = find_summit(plus, minus, ps, pe, w);
        let tag = if best_val > cutoff { "R" } else { "F" };
        if !out.is_empty() {
            out.push('\n');
        }
        // upstream: `b"%s\t%d\t%d\t%s\t%.2f"`, name suffixed `_R`/`_F`
        out.push_str(&format!(
            "{}\t{}\t{}\t{}_{}\t{:.2}",
            pk.chrom,
            best_pos,
            best_pos + 1,
            pk.name,
            tag,
            best_val
        ));
    }
    let ofile = o.get("ofile").map(str::to_string).unwrap_or_else(|| {
        format!(
            "{}_refinepeak.bed",
            o.get("oprefix").unwrap_or("macs3rscore")
        )
    });
    let outdir = o.get("outdir").unwrap_or(".");
    std::fs::create_dir_all(outdir)?;
    std::fs::write(std::path::Path::new(outdir).join(ofile), out)?;
    Ok(())
}
