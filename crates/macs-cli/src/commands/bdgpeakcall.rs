//! `macs3-rs bdgpeakcall` and `bdgbroadcall`: call peaks directly from a
//! bedGraph, with no treatment/control model.
//!
//! Both are ports of `bedGraphTrackI.call_peaks` / `call_broadpeaks`
//! (`BedGraph.py:405-527`): scan a bedGraph for runs at or above a cutoff, merge
//! runs separated by no more than `max_gap`, drop regions shorter than
//! `min_length`, and report the region with its summit (the lower median of the
//! midpoints of the runs holding the maximum value).
//!
//! Note the end-indexed convention that the callpeak pipeline shares: a peak's
//! **start is the end of the run immediately before the first above-cutoff run**
//! (`peak_content[0][0]`), and its end is the last above-cutoff run's end
//! (`peak_content[-1][1]`). The summit midpoint is `(tend + tstart) / 2` --
//! integer division, no `+1`.

use macs_bedgraph::BedGraph;
use macs_core::{Coord, Result};

use crate::Options;

/// One called peak from a bedGraph, in `PeakIO` field terms.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BdgPeak {
    pub start: Coord,
    pub end: Coord,
    pub summit: Coord,
    /// `peak_score`: the maximum bedGraph value inside the peak.
    pub score: f32,
}

/// `bedGraphTrackI.call_peaks`: regions of the track at or above `cutoff`.
///
/// Returns peaks per chromosome in name order, each with its own run list so the
/// broad level can re-segment at a second cutoff.
pub fn call_peaks_from_bedgraph(
    bg: &BedGraph,
    cutoff: f32,
    min_length: Coord,
    max_gap: Coord,
) -> Vec<(Vec<u8>, Vec<BdgPeak>)> {
    let mut out = Vec::new();
    for (chrom, t) in bg.iter_sorted() {
        let runs = t.runs();
        // peak_content entries are (start, end, value) in upstream's end-indexed
        // coordinates; the run before an above-cutoff run supplies its start.
        let mut regions: Vec<Vec<(Coord, Coord, f32)>> = Vec::new();
        let mut pre_p: Coord = 0;
        for r in runs {
            let p = r.end;
            let v = r.value;
            if v < cutoff {
                pre_p = p;
                continue;
            }
            match regions.last_mut() {
                Some(cur) if pre_p.saturating_sub(cur[cur.len() - 1].1) <= max_gap => {
                    cur.push((pre_p, p, v));
                }
                _ => regions.push(vec![(pre_p, p, v)]),
            }
            pre_p = p;
        }
        let mut peaks = Vec::new();
        for content in regions {
            if let Some(pk) = close_peak(&content, min_length) {
                peaks.push(pk);
            }
        }
        let name = bg.genome().name(chrom).to_vec();
        out.push((name, peaks));
    }
    out
}

/// `__close_peak`: keep the region if long enough, and compute its summit.
fn close_peak(content: &[(Coord, Coord, f32)], min_length: Coord) -> Option<BdgPeak> {
    let peak_length = content[content.len() - 1].1.saturating_sub(content[0].0);
    if peak_length < min_length {
        return None;
    }
    // summit = lower median of the midpoints of the maximal-value runs
    let mut tsummit: Vec<Coord> = Vec::new();
    let mut summit_value = 0.0f32;
    for &(tstart, tend, tvalue) in content {
        if summit_value == 0.0 || summit_value < tvalue {
            tsummit.clear();
            tsummit.push((tend + tstart) / 2);
            summit_value = tvalue;
        } else if summit_value == tvalue {
            tsummit.push((tend + tstart) / 2);
        }
    }
    // upstream is `(len + 1) / 2 - 1`, the lower median -- deliberately NOT
    // `div_ceil(2)`, which would pick the later of two tied runs (cf. F47).
    #[allow(clippy::manual_div_ceil)]
    let mid = (tsummit.len() + 1) / 2 - 1;
    Some(BdgPeak {
        start: content[0].0,
        end: content[content.len() - 1].1,
        summit: tsummit[mid],
        score: summit_value,
    })
}

/// Render peaks as an ENCODE narrowPeak file, matching upstream's
/// `write_to_narrowPeak` with `score_column="score"`.
/// `name` labels the track line (upstream uses `options.oprefix`); `peakprefix`
/// prefixes each peak name (`options.oprefix + "_narrowPeak"`). They are
/// different strings, and upstream uses each in exactly one place.
pub fn write_narrowpeak(
    peaks: &[(Vec<u8>, Vec<BdgPeak>)],
    name: &str,
    peakprefix: &str,
    trackline: bool,
) -> String {
    let mut out = String::new();
    if trackline {
        out.push_str(&format!(
            "track type=narrowPeak name=\"{name}\" description=\"{name}\" nextItemButton=on\n"
        ));
    }
    // upstream's `write_to_narrowPeak` (`PeakIO.py:719-735`) initialises
    // `n_peak = 0` once, before the chromosome loop, and increments it per
    // peak group -- so numbering is continuous across chromosomes.
    let mut n_peak = 0usize;
    for (chrom, list) in peaks {
        for p in list {
            n_peak += 1;
            let s = (p.summit as i64 - p.start as i64).max(-1);
            out.push_str(&format!(
                "{}\t{}\t{}\t{}{}\t{}\t.\t0\t0\t0\t{}\n",
                String::from_utf8_lossy(chrom),
                p.start,
                p.end,
                peakprefix,
                n_peak,
                (10.0 * f64::from(p.score)) as i64,
                s
            ));
        }
    }
    out
}

fn read_bedgraph(o: &Options) -> Result<BedGraph> {
    let ifile = o
        .get("ifile")
        .ok_or_else(|| macs_core::MacsError::InvalidParameter("-i/--ifile is required".into()))?;
    BedGraph::read(std::path::Path::new(ifile), 0.0)
}

/// `bdgpeakcall`.
pub fn bdgpeakcall(o: &Options) -> Result<()> {
    let bg = read_bedgraph(o)?;
    let cutoff = o.float("cutoff").unwrap_or(5.0) as f32;
    let minlen = o.int("minlen").unwrap_or(200) as Coord;
    let maxgap = o.int("maxgap").unwrap_or(30) as Coord;
    let oprefix = o.get("oprefix").unwrap_or("bdgpeakcall").to_string();
    let ofile = o.get("ofile").map(str::to_string);
    let outdir = o.get("outdir").unwrap_or(".");
    std::fs::create_dir_all(outdir)?;

    // `--cutoff-analysis` **replaces** peak calling, it does not accompany it:
    // `bdgpeakcall_cmd.py:44-56` is an `if/else`, and the peak call lives in the `else`.
    // So a run with the flag produces a report and no narrowPeak at all.
    //
    // Where the report lands is also not obvious: upstream writes it to
    // `options.ofile` when `-o` was given -- *overwriting* the path the peaks would
    // have used -- and only falls back to `<oprefix>_l<minlen>_g<maxgap>_cutoff_analysis.txt`
    // when `-o` is absent. Getting that backwards produces a file full of peaks where a
    // report belongs.
    if o.flag("cutoff_analysis") {
        let steps = o.int("cutoff_analysis_steps").unwrap_or(100).max(1) as usize;
        let max_score = o.int("cutoff_analysis_max").unwrap_or(100) as f32;
        // upstream passes `min_score=btrack.minvalue`, so the lower bound is the
        // track's own minimum rather than 0
        let min_value = bg
            .iter_sorted()
            .flat_map(|(_, t)| t.runs().iter().map(|r| r.value))
            .fold(f32::MAX, f32::min);
        let tracks: Vec<(macs_rle::SignalTrack<f32>, macs_core::Coord)> =
            bg.iter_sorted().map(|(_, t)| (t.clone(), 0u32)).collect();
        let report = macs_peaks::callpeak::bedgraph_cutoff_analysis(
            &tracks,
            maxgap.into(),
            minlen,
            steps,
            min_value,
            max_score,
        );
        let path = ofile
            .clone()
            .unwrap_or_else(|| format!("{oprefix}_l{minlen}_g{maxgap}_cutoff_analysis.txt"));
        std::fs::write(std::path::Path::new(outdir).join(&path), report)?;
        eprintln!("bdgpeakcall: cutoff analysis written to {path}");
        return Ok(());
    }

    let peaks = call_peaks_from_bedgraph(&bg, cutoff, minlen, maxgap);
    let total: usize = peaks.iter().map(|(_, v)| v.len()).sum();
    // F198: `bdgpeakcall_cmd.py:60` does `options.oprefix = options.ofile` before
    // writing, so `-o NAME` renames the *track* and every peak too, not just the file.
    // Left alone, `-o peaks.narrowPeak` writes a track line called `bdgpeakcall` into
    // `peaks.narrowPeak` -- same coordinates, wrong names.
    let (ofile, oprefix) = match ofile {
        Some(f) => {
            let name = f.clone();
            (f, name)
        }
        None => {
            let f = format!(
                "{}_c{:.1}_l{}_g{}_peaks.narrowPeak",
                oprefix, cutoff, minlen, maxgap
            );
            (f, oprefix)
        }
    };
    let trackline = o.flag("trackline"); // default true; --no-trackline sets false
    let body = write_narrowpeak(
        &peaks,
        &oprefix,
        &format!("{oprefix}_narrowPeak"),
        trackline,
    );
    std::fs::write(std::path::Path::new(outdir).join(&ofile), body)?;
    eprintln!("bdgpeakcall: {total} peaks");
    Ok(())
}

/// `bdgbroadcall`: build broad peaks by nesting strict (lvl1) regions inside
/// loose (lvl2) regions, then write the gappedPeak/bed12 format
/// (`__add_broadpeak`, `BedGraph.py:600-668`; `write_to_gappedPeak`,
/// `PeakIO.py:1370-1395`).
pub fn bdgbroadcall(o: &Options) -> Result<()> {
    let bg = read_bedgraph(o)?;
    let lvl1_cutoff = o.float("cutoffpeak").unwrap_or(2.0) as f32;
    let lvl2_cutoff = o.float("cutofflink").unwrap_or(1.0) as f32;
    if lvl1_cutoff <= lvl2_cutoff {
        return Err(macs_core::MacsError::InvalidParameter(
            "level 1 cutoff should be larger than level 2".into(),
        ));
    }
    let minlen = o.int("minlen").unwrap_or(200) as Coord;
    let lvl1_max_gap = o.int("lvl1maxgap").unwrap_or(30) as Coord;
    let lvl2_max_gap = o.int("lvl2maxgap").unwrap_or(800) as Coord;
    if lvl1_max_gap >= lvl2_max_gap {
        return Err(macs_core::MacsError::InvalidParameter(
            "level 2 maximum gap should be larger than level 1".into(),
        ));
    }

    let lvl2_all = call_peaks_from_bedgraph(&bg, lvl2_cutoff, minlen, lvl2_max_gap);
    let lvl1_all = call_peaks_from_bedgraph(&bg, lvl1_cutoff, minlen, lvl1_max_gap);

    let mut out = String::new();
    let trackline = o.flag("trackline"); // default true; --no-trackline sets false
                                         // `bdgbroadcall_cmd.py:53-58`: when `-o` is given, `options.oprefix` is
                                         // overwritten with the ofile name before writing, so the peak prefix uses
                                         // the ofile -- mirroring `bdgpeakcall`'s F198 handling above.
    let ofile_opt = o.get("ofile").map(str::to_string);
    let oprefix = match &ofile_opt {
        Some(f) => f.clone(),
        None => o.get("oprefix").unwrap_or("bdgbroadcall").to_string(),
    };
    if trackline {
        // upstream's write_to_gappedPeak defaults name/description to "peak"
        // (the command does not pass them), so the track line does not use the
        // o-prefix.
        out.push_str(
            "track name=\"peak\" description=\"peak\" type=gappedPeak nextItemButton=on\n",
        );
    }
    let peakprefix = format!("{oprefix}_broadRegion");
    let mut n_peak = 0usize;
    let mut total = 0usize;
    for (chrom, lvl2_list) in &lvl2_all {
        // `call_broadpeaks` (`BedGraph.py:564-567`) iterates only the
        // chromosomes of the lvl1 peaks, so a chromosome with lvl2 regions but
        // no lvl1 peak contributes nothing upstream. Both lists come from the
        // same track, so the chromosome is always present here; only an empty
        // lvl1 list (no peaks on this chromosome) skips output.
        let (_, lvl1_list) = lvl1_all
            .iter()
            .find(|(c, _)| c == chrom)
            .expect("lvl1 and lvl2 share chromosomes");
        // `call_peaks_from_bedgraph` returns an entry per input chromosome even
        // when it holds no peaks; upstream's `PeakIO` only contains chromosomes
        // with at least one peak, and `call_broadpeaks` iterates those -- so an
        // empty lvl1 list also means "skip this chromosome".
        if lvl1_list.is_empty() {
            continue;
        }
        let inner: Vec<&BdgPeak> = lvl1_list.iter().collect();
        for p2 in lvl2_list {
            let set: Vec<&&BdgPeak> = inner
                .iter()
                .filter(|p| p.start >= p2.start && p.end <= p2.end)
                .collect();
            let (start, end) = (p2.start, p2.end);
            let (_thick_start, block_num, block_sizes, block_starts) = if set.is_empty() {
                // no strong peaks inside: complement with 1bp blocks
                (
                    start,
                    2usize,
                    "1,1".to_string(),
                    format!("0,{}", end - start - 1),
                )
            } else {
                let mut ts = set[0].start;
                let te = set[set.len() - 1].end;
                let mut bn = set.len();
                let mut sizes: Vec<String> =
                    set.iter().map(|p| (p.end - p.start).to_string()).collect();
                let mut starts: Vec<String> =
                    set.iter().map(|p| (p.start - start).to_string()).collect();
                if ts != start {
                    ts = start;
                    bn += 1;
                    sizes.insert(0, "1".into());
                    starts.insert(0, "0".into());
                }
                if te != end {
                    bn += 1;
                    sizes.push("1".into());
                    starts.push(format!("{}", end - start - 1));
                }
                (ts, bn, sizes.join(","), starts.join(","))
            };
            n_peak += 1;
            total += 1;
            out.push_str(&format!(
                "{}\t{}\t{}\t{}{}\t{}\t.\t0\t0\t0\t{}\t{}\t{}\t0\t0\t0\n",
                String::from_utf8_lossy(chrom),
                start,
                end,
                peakprefix,
                n_peak,
                (10.0 * f64::from(p2.score)) as i64,
                block_num,
                block_sizes,
                block_starts
            ));
        }
    }

    let ofile = ofile_opt.unwrap_or_else(|| {
        format!(
            "{}_c{:.1}_C{:.2}_l{}_g{}_G{}_broad.bed12",
            oprefix, lvl1_cutoff, lvl2_cutoff, minlen, lvl1_max_gap, lvl2_max_gap
        )
    });
    let outdir = o.get("outdir").unwrap_or(".");
    std::fs::create_dir_all(outdir)?;
    std::fs::write(std::path::Path::new(outdir).join(ofile), out)?;
    eprintln!("bdgbroadcall: {total} broad peaks");
    Ok(())
}
