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

use crate::Options;

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

fn counter(pos: &[u32]) -> HashMap<i64, i64> {
    let mut m: HashMap<i64, i64> = HashMap::new();
    for &p in pos {
        *m.entry(i64::from(p)).or_insert(0) += 1;
    }
    m
}

/// The tags of one peak's window, as `compute_region_tags_from_peaks` hands them
/// to `find_summit` (`FixWidthTrack.py:663-691`).
///
/// This is **not** "every tag in `[startpos, endpos]`". Upstream walks a single
/// forward cursor per chromosome and strand, keeps it across peaks, and rewinds
/// it by less than `window_size` afterwards (`FixWidthTrack.py:695-704`):
///
/// ```text
/// for i in range(prev_i, plus.shape[0]):
///     pos = plus[i]
///     if pos < startpos:   continue
///     elif pos > endpos:   prev_i = i; break
///     else:                temp.append(pos)
/// ...
/// for i in range(prev_i, 0, -1):
///     if plus[prev_i] - plus[i] >= window_size: break
/// prev_i = i
/// ```
///
/// Two consequences have to be reproduced verbatim, because both decide which
/// tags a later peak sees:
///
/// 1. `range(prev_i, 0, -1)` never visits index `0`, so a rewind that does not
///    `break` leaves `prev_i == 1` -- index 0 of that strand is dropped for
///    every later peak on the chromosome.
/// 2. When the rewind loop body never runs (`prev_i == 0`) the assignment reads
///    the collection loop's own loop variable, which is the `break` index, or
///    the last index if the collection loop ran off the end. It is not a
///    no-op, so `prev_i` can jump forward past the tags a later peak needs.
///
/// On the 5 M-read CTCF fixture this is worth 24 of 36769 refined summits: peak
/// 1436 (`chr1 86968634 86968834`) gets 1 plus and 3 minus tags instead of 5 and
/// 7, so MACS3 reports `2.8284271247461903` (which fails the `--cutoff`) where a
/// stateless window reports `5.745966692414834`.
///
/// `cursor` is the per-chromosome, per-strand index carried between peaks.
fn collect_region_tags(
    tags: &[u32],
    startpos: i64,
    endpos: i64,
    w: i64,
    cursor: &mut usize,
) -> Vec<u32> {
    let mut out = Vec::new();
    let mut i = *cursor;
    while i < tags.len() {
        let pos = i64::from(tags[i]);
        if pos < startpos {
            i += 1;
        } else if pos > endpos {
            *cursor = i;
            break;
        } else {
            out.push(tags[i]);
            i += 1;
        }
    }
    if *cursor == 0 {
        // the rewind loop is `range(0, 0, -1)`: empty, so it re-assigns the
        // index the collection loop stopped on. `i` is already that index when
        // it broke, and one past the end when it exhausted the array (upstream
        // leaves it on the last index).
        *cursor = i.min(tags.len().saturating_sub(1));
    } else {
        // walk back to the first index within `w` of the cursor; upstream's
        // `range` stops at 1, never 0.
        let mut j = *cursor;
        loop {
            if i64::from(tags[*cursor]) - i64::from(tags[j]) >= w {
                *cursor = j;
                break;
            }
            if j <= 1 {
                *cursor = 1;
                break;
            }
            j -= 1;
        }
    }
    out
}

/// `find_summit`: running WTD over the window, returning `(best_pos, best_val)`.
fn find_summit(plus: &[u32], minus: &[u32], peak_start: i64, peak_end: i64, w: i64) -> (i64, f64) {
    // The counters hold the tags `collect_region_tags` returned, i.e. those
    // inside `[peak_start, peak_end]`. The `w`-expansion happens inside
    // `sum_le`/`sum_ge`, not here.
    let watson = counter(plus);
    let crick = counter(minus);
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
        // `(watson_left * crick_right)**0.5` is `pow(x, 0.5)`, not `sqrt`
        // (`refinepeak_cmd.py:88`), and `powf` is the same libm call.
        let v = 2.0 * ((wl as f64) * (cr as f64)).powf(0.5) - (wr as f64) - (cl as f64);
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
    let ifiles = super::input::input_files(o)?;
    let format = o.get("format").unwrap_or("AUTO").to_uppercase();
    if format == "BAMPE" || format == "BEDPE" {
        return Err(MacsError::InvalidParameter(
            "paired-end input is not yet supported".into(),
        ));
    }
    let track = super::input::load_single_end_files(&ifiles, &format)?.0;
    let peaks = read_peaks(Path::new(bedfile))?;
    let w = o.int("windowsize").unwrap_or(200);
    let cutoff = o.float("cutoff").unwrap_or(5.0);

    // Group tags by chromosome once.
    //
    // The vectors hold **every** tag on the chromosome, unfiltered: upstream
    // keeps the full `plus`/`minus` arrays and only ever *indexes* into them
    // (`FixWidthTrack.py:661,671-690`). The per-peak window is applied by
    // `collect_region_tags` below, and that function has to see the same
    // indices upstream does -- the cursor it carries is an index into these
    // arrays, so any pre-filtering would move it.
    let mut by_chrom: HashMap<String, (Vec<u32>, Vec<u32>)> = HashMap::new();
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

    // `compute_region_tags_from_peaks` walks chromosomes in name order and peaks
    // in ascending start order (`FixWidthTrack.py:659,665`), which is the order
    // `peaks` is already in, and carries `prev_i`/`prev_j` across the peaks of a
    // chromosome -- so both cursors reset when the chromosome changes.
    let mut cursors: HashMap<&str, (usize, usize)> = HashMap::new();
    let mut out = String::new();
    for pk in &peaks {
        let Some(tags) = by_chrom.get(&pk.chrom) else {
            continue;
        };
        let (prev_i, prev_j) = cursors.entry(pk.chrom.as_str()).or_insert((0, 0));
        // the window over which find_summit runs
        let ps = pk.start - w;
        let pe = pk.end + w;
        let rt_plus = collect_region_tags(&tags.0, ps, pe, w, prev_i);
        let rt_minus = collect_region_tags(&tags.1, ps, pe, w, prev_j);
        let (best_pos, best_val) = find_summit(&rt_plus, &rt_minus, ps, pe, w);
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
