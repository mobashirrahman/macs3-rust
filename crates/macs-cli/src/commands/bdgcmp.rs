//! `macs3-rs bdgcmp` and `bdgdiff`: the two differential bedGraph commands.
//!
//! `bdgcmp` compares one treatment against one control and emits one bedGraph
//! per scoring method. `bdgdiff` compares two conditions (each with its own
//! control) and emits three BED files. Both are thin drivers over
//! `macs_score`'s `ScoreTrack2`/`TwoScores` and `macs_bedgraph`'s merge, which
//! are ports of `MACS3/Signal/ScoreTrack.py`.

use std::path::{Path, PathBuf};

use macs_bedgraph::BedGraph;
use macs_core::{ChromId, MacsError, Result};
use macs_io::format_g;
use macs_score::{NormMethod, ScoreMethod};

use crate::Options;

/// The eight score methods `bdgcmp -m` accepts, in upstream's order.
const BDCMP_METHODS: [(&str, ScoreMethod); 8] = [
    ("ppois", ScoreMethod::P),
    ("qpois", ScoreMethod::Q),
    ("subtract", ScoreMethod::Subtract),
    ("logFE", ScoreMethod::LogFE),
    ("FE", ScoreMethod::FE),
    ("logLR", ScoreMethod::LogLR),
    ("slogLR", ScoreMethod::SymLogLR),
    ("max", ScoreMethod::Max),
];

fn read_bed(path: &Path) -> Result<BedGraph> {
    BedGraph::read(path, 0.0)
}

/// The output directory, created if missing.
///
/// Upstream's `--outdir` has no `type=`, no validator and no `os.makedirs` in
/// `add_outdir_option` (`bin/macs3:168-170`), yet `bdgdiff --outdir nodir` exits 0
/// and leaves `nodir/` behind -- something further down creates it. This port
/// failed the write instead (exit 1), so it is a divergence in the stricter
/// direction: the same invocation succeeded upstream and failed here.
fn outdir(o: &Options) -> Result<PathBuf> {
    let d = PathBuf::from(o.get("outdir").unwrap_or(""));
    if !d.as_os_str().is_empty() {
        std::fs::create_dir_all(&d)?;
    }
    Ok(d)
}

fn required<'a>(o: &'a Options, key: &str, flag: &str) -> Result<&'a str> {
    o.get(key)
        .ok_or_else(|| MacsError::InvalidParameter(format!("{flag} is required")))
}

/// `bdgcmp`: score one treatment against one control.
///
/// Follows `bdgcmp_cmd.run` exactly, including two upstream details:
///
/// 1. `-S/--scaling-factor` is inverted (`pseudo_depth = 1/S`) to reach the
///    `Million` normalisation, and the override is only applied when
///    `|S - 1| > 1e-6`.
/// 2. Repeated methods are **skipped**, not recomputed: only the first
///    occurrence of each method in `-m` order produces an output file.
pub fn bdgcmp(o: &Options) -> Result<()> {
    let tfile = required(o, "tfile", "-t/--tfile")?;
    let cfile = required(o, "cfile", "-c/--cfile")?;
    let sfactor: f64 = o.float("sfactor").unwrap_or(1.0);
    // upstream: `pseudo_depth = 1.0/scaling_factor`, "a trick to override SPMR"
    let pseudo_depth = (1.0 / sfactor) as f32;

    let tbtrack = read_bed(Path::new(tfile))?;
    let cbtrack = read_bed(Path::new(cfile))?;

    let mut sbtrack = tbtrack.make_score_track_for_macs(&cbtrack, pseudo_depth, pseudo_depth);
    if (sfactor - 1.0).abs() > 1e-6 {
        eprintln!(
            "Values in your input bedGraph files will be multiplied by {:.6} ...",
            sfactor
        );
        sbtrack.change_normalization_method(NormMethod::Million);
    }
    sbtrack.set_pseudocount(o.float("pseudocount").unwrap_or(0.0) as f32);

    let methods: Vec<String> = o.get_all("method").to_vec();
    let methods = if methods.is_empty() {
        vec!["ppois".to_string()]
    } else {
        methods
    };
    for m in &methods {
        if !BDCMP_METHODS.iter().any(|(n, _)| n == m) {
            return Err(MacsError::InvalidParameter(format!("Invalid method: {m}")));
        }
    }
    let ofiles: Vec<String> = o.get_all("ofile").to_vec();
    if !ofiles.is_empty() && ofiles.len() != methods.len() {
        return Err(MacsError::InvalidParameter(
            "The number and the order of arguments for --ofile must be the same as for -m.".into(),
        ));
    }
    let oprefix = o.get("oprefix").unwrap_or("bdgcmp").to_string();

    let dir = outdir(o)?;
    let mut done: Vec<String> = Vec::new();
    for (i, method) in methods.iter().enumerate() {
        if done.contains(method) {
            continue;
        }
        done.push(method.clone());
        let ofile = if ofiles.is_empty() {
            format!("{oprefix}_{method}.bdg")
        } else {
            ofiles[i].clone()
        };
        let sm = BDCMP_METHODS
            .iter()
            .find(|(n, _)| n == method)
            .expect("validated")
            .1;
        sbtrack.change_score_method(sm);
        let upper = method.to_uppercase();
        sbtrack.write_bedgraph(
            &dir.join(&ofile),
            &format!("{upper}_Scores"),
            &format!("Scores calculated by {upper}"),
            3,
        )?;
    }
    Ok(())
}

/// `bdgdiff`: differential peaks between two conditions.
///
/// Follows `bdgdiff_cmd.run`: the two depths are normalised so the larger
/// sample is scaled down to the smaller one (leaving the other at 1.0), the
/// three peak sets come from `TwoScores::call_peaks`, and the filenames use a
/// `%.1f` cutoff when `--o-prefix` is given.
pub fn bdgdiff(o: &Options) -> Result<()> {
    let t1 = required(o, "t1bdg", "--t1")?;
    let c1 = required(o, "c1bdg", "--c1")?;
    let t2 = required(o, "t2bdg", "--t2")?;
    let c2 = required(o, "c2bdg", "--c2")?;
    let cutoff = o.float("cutoff").unwrap_or(3.0) as f32;
    let minlen = o.int("minlen").unwrap_or(200).max(0) as u64;
    let maxgap = o.int("maxgap").unwrap_or(100).max(0) as u64;
    if maxgap >= minlen {
        return Err(MacsError::InvalidParameter(format!(
            "MAXGAP should be smaller than MINLEN! Your input is MAXGAP = {maxgap} and MINLEN = {minlen}"
        )));
    }
    let depth1 = o.float("depth1").unwrap_or(1.0);
    let depth2 = o.float("depth2").unwrap_or(1.0);
    // scale the larger condition down to the level of the smaller one
    let (f1, f2) = if depth1 > depth2 {
        (depth2 / depth1, 1.0)
    } else if depth1 < depth2 {
        (1.0, depth1 / depth2)
    } else {
        (1.0, 1.0)
    };

    let t1bdg = read_bed(Path::new(t1))?;
    let c1bdg = read_bed(Path::new(c1))?;
    let t2bdg = read_bed(Path::new(t2))?;
    let c2bdg = read_bed(Path::new(c2))?;
    // all four share one genome so chromosome ids agree
    let mut genome = t1bdg.genome().clone();
    for bg in [&c1bdg, &t2bdg, &c2bdg] {
        for chrom in bg.chroms_sorted() {
            let _ = genome.intern(bg.genome().name(chrom));
        }
    }

    let scores = macs_score::TwoScores::build(
        genome,
        t1bdg.tracks(),
        c1bdg.tracks(),
        t2bdg.tracks(),
        c2bdg.tracks(),
        f1 as f32,
        f2 as f32,
        0.01,
    );
    let (cat1, cat2, cat3) = scores.call_peaks(cutoff, minlen, maxgap);

    let dir = outdir(o)?;
    let ofiles: Vec<String> = o.get_all("ofile").to_vec();
    let oprefix = o.get("oprefix").unwrap_or("bdgdiff").to_string();
    let (paths, names): (Vec<PathBuf>, Vec<(String, &str, &str)>) = if !ofiles.is_empty() {
        (
            vec![
                dir.join(&ofiles[0]),
                dir.join(&ofiles[1]),
                dir.join(&ofiles[2]),
            ],
            vec![
                (
                    ofiles[0].clone(),
                    "condition 1",
                    "unique regions in condition 1",
                ),
                (
                    ofiles[1].clone(),
                    "condition 2",
                    "unique regions in condition 2",
                ),
                (
                    ofiles[2].clone(),
                    "common",
                    "common regions in both conditions",
                ),
            ],
        )
    } else {
        let c = format!("{cutoff:.1}");
        (
            vec![
                dir.join(format!("{oprefix}_c{c}_cond1.bed")),
                dir.join(format!("{oprefix}_c{c}_cond2.bed")),
                dir.join(format!("{oprefix}_c{c}_common.bed")),
            ],
            vec![
                (
                    format!("{oprefix}_cond1_"),
                    "condition 1",
                    "unique regions in condition 1",
                ),
                (
                    format!("{oprefix}_cond2_"),
                    "condition 2",
                    "unique regions in condition 2",
                ),
                (
                    format!("{oprefix}_common_"),
                    "common",
                    "common regions in both conditions",
                ),
            ],
        )
    };
    for (i, peaks) in [&cat1, &cat2, &cat3].into_iter().enumerate() {
        write_diff_bed(&paths[i], &names[i], scores.genome(), peaks)?;
    }
    Ok(())
}

/// Write one BED5 file of differential peaks, matching `PeakIO._to_bed`.
///
/// Peak names are `PREFIX_N` with `N` counting across chromosomes in name
/// order; the score column is `%.6g` (significant digits, not fixed decimals),
/// so `format_g` rather than a fixed-precision formatter. Upstream groups peaks by `end` and only
/// adds a subpeak letter when two peaks share an end -- `DiffPeak`s from
/// `call_peaks` are already disjoint and unique-ended, so no letters appear
/// here, but the grouping is preserved so a future caller that can produce ties
/// stays consistent.
fn write_diff_bed(
    path: &Path,
    spec: &(String, &str, &str),
    genome: &macs_core::Genome,
    peaks: &[macs_score::DiffPeak],
) -> Result<()> {
    let (prefix, name, description) = spec;
    let mut out = String::new();
    out.push_str(&format!(
        "track name=\"{} (peaks)\" description=\"{}\" visibility=1\n",
        name.replace('"', "\\\""),
        description.replace('"', "\\\"")
    ));
    // group by (chrom, end) so tie handling matches upstream's `groupby`
    let mut by_chrom: std::collections::BTreeMap<ChromId, Vec<&macs_score::DiffPeak>> =
        std::collections::BTreeMap::new();
    for p in peaks {
        by_chrom.entry(p.chrom).or_default().push(p);
    }
    let mut chroms: Vec<ChromId> = by_chrom.keys().copied().collect();
    chroms.sort_by(|&a, &b| genome.name(a).cmp(genome.name(b)));
    let mut n_peak = 0usize;
    for chrom in chroms {
        let list = &by_chrom[&chrom];
        let mut ends: Vec<u64> = list.iter().map(|p| p.end).collect();
        ends.sort_unstable();
        ends.dedup();
        for end in ends {
            let group: Vec<&&macs_score::DiffPeak> = list.iter().filter(|p| p.end == end).collect();
            n_peak += 1;
            let cname = String::from_utf8_lossy(genome.name(chrom));
            for (i, p) in group.iter().enumerate() {
                if group.len() > 1 {
                    out.push_str(&format!(
                        "{}\t{}\t{}\t{}{}{}\t{}\n",
                        cname,
                        p.start,
                        p.end,
                        prefix,
                        n_peak,
                        subpeak_letter(i),
                        format_g(f64::from(p.score), 6)
                    ));
                } else {
                    out.push_str(&format!(
                        "{}\t{}\t{}\t{}{}\t{}\n",
                        cname,
                        p.start,
                        p.end,
                        prefix,
                        n_peak,
                        format_g(f64::from(p.score), 6)
                    ));
                }
            }
        }
    }
    std::fs::write(path, out).map_err(MacsError::Io)
}

/// `subpeak_letters` (`PeakIO.py`): `a`, `b`, ... `z`, `A`, `B`, ...
fn subpeak_letter(i: usize) -> char {
    const LOWER: usize = 26;
    if i < LOWER {
        (b'a' + i as u8) as char
    } else {
        (b'A' + (i - LOWER) as u8) as char
    }
}
