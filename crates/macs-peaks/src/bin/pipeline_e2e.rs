//! End-to-end: fixture BED -> pileup -> paired arrays, against upstream's own
//! `--bdg` output.
//!
//! This is the check that the *stages agree*, not just that each stage matches
//! upstream in isolation. It runs the real pipeline -- read, dedup, pile up at
//! every local-lambda scale, merge, pair -- and compares the resulting treatment
//! pileup and control lambda against the bedGraph files upstream writes for the
//! same invocation.
//!
//! Usage:
//! ```text
//! macs-io-pipeline-e2e <fixture-dir> --extsize 200 --treat-pileup F \
//!                                        --control-lambda F [--summits F]
//! ```
//!
//! `--bdg` is what makes the comparison possible: `callpeak --bdg` writes
//! `<name>_treat_pileup.bdg` and `<name>_control_lambda.bdg`, which are exactly the
//! merged arrays `CallerFromAlignments` hands to the scoring stage, rendered with
//! `%.5f`.

use macs_core::{MacsError, Result, Strand};
use macs_peaks::LambdaScales;
use macs_pileup::SingleEndParams;
use macs_track::{SingleEndTrack, SingleEndTrackBuilder};

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Parse a BED file into a deduplicated single-end track.
///
/// Mirrors `callpeak`'s `--nomodel --extsize D` path: parse with the BED rules,
/// sort, then `filter_dup(1)`.
fn load_bed(path: &Path) -> Result<SingleEndTrack> {
    let file = std::fs::File::open(path).map_err(|e| {
        MacsError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })?;
    let mut reader = BufReader::new(file);
    let mut line = Vec::new();
    let mut b = SingleEndTrackBuilder::new();
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        let Some(rec) = macs_io::parse_bed_line(&line)? else {
            continue;
        };
        if rec.pos < 0 || rec.chrom.is_empty() {
            continue;
        }
        b.push(&rec.chrom, rec.pos as u64, rec.strand);
    }
    b.finalize();
    let mut t = b.build();
    // `--keep-dup` defaults to 1, so the corpus is deduplicated to one copy per
    // (position, strand)
    t.filter_dup(1)?;
    Ok(t)
}

/// A parsed bedGraph: `(chrom, start, end, value)` runs in file order.
type BedGraph = BTreeMap<String, Vec<(u64, u64, f64)>>;

fn load_bedgraph(path: &Path) -> Result<BedGraph> {
    let file = std::fs::File::open(path).map_err(|e| {
        MacsError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })?;
    let mut out: BTreeMap<String, Vec<(u64, u64, f64)>> = BTreeMap::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with("track") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            continue;
        }
        let (Ok(start), Ok(end), Ok(v)) = (
            f[1].parse::<u64>(),
            f[2].parse::<u64>(),
            f[3].parse::<f64>(),
        ) else {
            continue;
        };
        out.entry(f[0].to_string())
            .or_default()
            .push((start, end, v));
    }
    Ok(out)
}

fn main() -> std::process::ExitCode {
    match run() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("macs-io-pipeline-e2e: {e}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run() -> Result<()> {
    let mut fixture = String::new();
    let mut treat_pileup = String::new();
    let mut control_lambda = String::new();
    let mut extsize: i64 = 200;
    let mut gsize: f64 = 1_000_000.0;
    let mut slocal: i64 = 1000;
    let mut llocal: i64 = 10_000;

    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--extsize" => {
                extsize = args[i + 1].parse().expect("extsize");
                i += 2;
            }
            "--gsize" => {
                gsize = args[i + 1].parse().expect("gsize");
                i += 2;
            }
            "--slocal" => {
                slocal = args[i + 1].parse().expect("slocal");
                i += 2;
            }
            "--llocal" => {
                llocal = args[i + 1].parse().expect("llocal");
                i += 2;
            }
            "--treat-pileup" => {
                treat_pileup = args[i + 1].clone();
                i += 2;
            }
            "--control-lambda" => {
                control_lambda = args[i + 1].clone();
                i += 2;
            }
            other => {
                fixture = other.to_string();
                i += 1;
            }
        }
    }
    if fixture.is_empty() {
        return Err(MacsError::InvalidParameter(
            "usage: macs-io-pipeline-e2e <fixture-dir> [--extsize N] --treat-pileup F \
             --control-lambda F"
                .into(),
        ));
    }

    let dir = Path::new(&fixture);
    let treat = load_bed(&dir.join("treat.bed"))?;
    let ctrl_path = dir.join("ctrl.bed");
    let ctrl = if ctrl_path.exists() {
        Some(load_bed(&ctrl_path)?)
    } else {
        None
    };

    // `rlength`: upstream uses INT_MAX when no chrom size file is given, which
    // means endpoints are only clipped at zero. The fixtures carry a
    // `genome.txt`, but `--gsize` is what `callpeak` is given here, so mirror
    // upstream and leave the right edge unclipped.
    let rlength = u64::MAX / 2;

    // lambda_bg = treat_sum / gsize when scaling treatment to control, which is
    // what `callpeak` does when the control is larger (the default for
    // `--to-small`).
    let t_total = treat.total() as f64;
    let c_total = ctrl.as_ref().map(|c| c.total()).unwrap_or(0) as f64;
    // `callpeak --to-small` scales the *larger* sample down to the smaller one.
    // `PeakDetect` names this `tocontrol`, which reads as "scale to control" and is
    // true only when the **treatment** is the larger one. So:
    //
    //   treat > ctrl  ->  tocontrol = true   (treatment scaled down, lambda_bg from control_sum)
    //   treat < ctrl  ->  tocontrol = false  (control scaled down,    lambda_bg from treat_sum)
    //
    // Getting this backwards scales the wrong sample and sets lambda_bg from the
    // wrong total, which shifts every control-lambda value.
    let to_control = t_total > c_total;
    let treat_sum = t_total * extsize as f64;
    let control_sum = c_total * extsize as f64;
    let ratio = if control_sum > 0.0 {
        treat_sum / control_sum
    } else {
        0.0
    };
    let lambda_bg = if to_control {
        (control_sum / gsize) as f32
    } else {
        (treat_sum / gsize) as f32
    };
    let treat_scale: f32 = if to_control && ratio != 0.0 {
        (1.0 / ratio) as f32
    } else {
        1.0
    };
    let _ = treat_scale;

    // The control scaling factors for the three local-lambda scales.
    let (scale_ds, scale_factors) =
        control_scale_factors(ratio, to_control, extsize, slocal, llocal);
    let mut scales = LambdaScales {
        d: scale_ds[0],
        slocal: scale_ds.get(1).copied().unwrap_or(0),
        llocal: scale_ds.get(2).copied().unwrap_or(0),
        d_factor: scale_factors[0],
        slocal_factor: scale_factors.get(1).copied().unwrap_or(0.0),
        llocal_factor: scale_factors.get(2).copied().unwrap_or(0.0),
    };
    // `as_pairs` always emits three entries; drop any the caller disabled
    if scale_ds.len() < 3 {
        scales.llocal = 0;
    }

    // ---- compare the pileups ----
    let mut problems: Vec<String> = Vec::new();
    let mut notes: Vec<String> = Vec::new();
    let mut worst_treat: f64 = 0.0;
    let mut worst_ctrl: f64 = 0.0;
    let mut n_treat = 0usize;
    let mut n_ctrl = 0usize;

    if !treat_pileup.is_empty() {
        let want = load_bedgraph(Path::new(&treat_pileup))?;
        let got = build_treat_bedgraph(&treat, rlength, extsize, treat_scale);
        let (n, worst, probs, nts) = compare_graph("treat_pileup", &want, &got);
        n_treat = n;
        worst_treat = worst;
        problems.extend(probs);
        notes.extend(nts);
    }
    if !control_lambda.is_empty() {
        let want = load_bedgraph(Path::new(&control_lambda))?;
        // Upstream's `--bdg` output comes from `chr_pos_treat_ctrl`, i.e. the
        // **paired** arrays, and the pairing stops wherever the shorter of
        // treatment and control runs out (F35). The treatment ends at 93003 here,
        // so upstream's control bedGraph is truncated there even though the raw
        // control lambda continues. Truncate the same way before comparing,
        // otherwise the tail reads as a mismatch that is really the pairing.
        let treat_end = treat
            .positions()
            .chroms()
            .iter()
            .filter_map(|c| {
                macs_pileup::pileup_from_positions(
                    *c,
                    treat.positions().strand(*c, Strand::Plus),
                    treat.positions().strand(*c, Strand::Minus),
                    &SingleEndParams::directional(extsize, 0, rlength, treat_scale),
                )
                .runs()
                .last()
                .map(|r| r.end)
            })
            .max()
            .unwrap_or(0);
        let got = build_ctrl_bedgraph(
            ctrl.as_ref(),
            rlength,
            &scales,
            lambda_bg,
            treat_sum,
            control_sum,
            gsize,
            treat_end,
        );
        let (n, worst, probs, nts) = compare_graph("control_lambda", &want, &got);
        n_ctrl = n;
        worst_ctrl = worst;
        problems.extend(probs);
        notes.extend(nts);
    }

    println!(
        "chromosomes: treat={} ctrl={}",
        treat.total(),
        ctrl.as_ref().map(|c| c.total()).unwrap_or(0)
    );
    println!("scaling: tocontrol={to_control} ratio={ratio:.6} lambda_bg={lambda_bg} treat_scale={treat_scale}");
    println!(
        "lambda scales: d_s={:?} factors={:?}",
        scale_ds, scale_factors
    );
    println!("treat pileup runs compared: {n_treat}, max abs diff {worst_treat:.3e}");
    println!("control lambda runs compared: {n_ctrl}, max abs diff {worst_ctrl:.3e}");
    for p in notes.iter().take(20) {
        println!("NOTE {p}");
    }
    for p in problems.iter().take(20) {
        println!("MISMATCH {p}");
    }
    println!("total mismatches: {}", problems.len());
    if problems.is_empty() {
        println!("E2E PILEUP MATCH");
        Ok(())
    } else {
        Err(MacsError::Internal(format!(
            "{} pileup mismatches",
            problems.len()
        )))
    }
}

/// The three control scaling factors, transcribed from
/// `PeakDetect.__call_peaks_w_control`.
///
/// ```text
/// ctrl_d_s = [d]
/// # d scale
/// tmp_v = ratio_treat2control  if not tocontrol   else 1.0
/// ctrl_scale_s.append(tmp_v)
/// # slocal
/// ctrl_d_s.append(sregion)
/// tmp_v = d/sregion*ratio       if not tocontrol   else d/sregion
/// # llocal, only when lregion > sregion
/// ctrl_d_s.append(lregion)
/// tmp_v = d/lregion*ratio       if not tocontrol   else d/lregion
/// ```
///
/// Two details that are easy to miss: the local-window factors are the window
/// size **divided into** the region (`d / sregion`), not the other way round; and
/// the `llocal` scale is omitted entirely when `lregion <= sregion`, which is why
/// `treat_scale` does not appear here at all — upstream folds it into `ratio`
/// before this point.
#[allow(clippy::too_many_arguments)]
fn control_scale_factors(
    ratio: f64,
    to_control: bool,
    d: i64,
    slocal: i64,
    llocal: i64,
) -> (Vec<i64>, Vec<f32>) {
    let mut ds = vec![d];
    let mut factors: Vec<f32> = Vec::new();

    factors.push(if to_control { 1.0 } else { ratio as f32 });

    if slocal > 0 {
        ds.push(slocal);
        factors.push(
            (if to_control {
                d as f64 / slocal as f64
            } else {
                d as f64 / slocal as f64 * ratio
            }) as f32,
        );
    }

    if llocal > 0 && llocal > slocal {
        ds.push(llocal);
        factors.push(
            (if to_control {
                d as f64 / llocal as f64
            } else {
                d as f64 / llocal as f64 * ratio
            }) as f32,
        );
    }

    (ds, factors)
}

/// Render the treatment pileup as a bedGraph, for comparison.
/// Render the treatment pileup as a bedGraph, for comparison.
///
/// Note there is **no baseline argument**: upstream builds the treatment pileup
/// with `baseline_value=0.0` and gives `lambda_bg` only to the control. Passing it
/// to treatment as well raises every zero-depth base, and zero-depth is most of a
/// sparse ChIP pileup.
fn build_treat_bedgraph(
    track: &SingleEndTrack,
    rlength: u64,
    extsize: i64,
    scale: f32,
) -> BedGraph {
    let mut out: BedGraph = BTreeMap::new();
    for chrom in track.positions().chroms() {
        let name = String::from_utf8_lossy(track.genome().name(chrom)).into_owned();
        let t = macs_pileup::pileup_from_positions(
            chrom,
            track.positions().strand(chrom, Strand::Plus),
            track.positions().strand(chrom, Strand::Minus),
            &SingleEndParams::directional(extsize, 0, rlength, scale),
        );
        let mut runs = Vec::new();
        let mut prev = 0u64;
        for r in t.runs() {
            runs.push((prev, r.end, f64::from(r.value)));
            prev = r.end;
        }
        out.insert(name, runs);
    }
    out
}

/// Render the merged control lambda as a bedGraph, for comparison.
#[allow(clippy::too_many_arguments)]
fn build_ctrl_bedgraph(
    ctrl: Option<&SingleEndTrack>,
    rlength: u64,
    scales: &LambdaScales,
    lambda_bg: f32,
    treat_sum: f64,
    control_sum: f64,
    gsize: f64,
    limit: u64,
) -> BedGraph {
    let mut out: BedGraph = BTreeMap::new();
    let Some(c) = ctrl else {
        return out;
    };
    let _ = (treat_sum, control_sum, gsize);
    for chrom in c.positions().chroms() {
        let name = String::from_utf8_lossy(c.genome().name(chrom)).into_owned();
        let mut combined = None::<macs_rle::SignalTrack<f32>>;
        for scale in scales.as_pairs() {
            let p = macs_pileup::pileup_from_positions(
                chrom,
                c.positions().strand(chrom, Strand::Plus),
                c.positions().strand(chrom, Strand::Minus),
                &SingleEndParams::bidirectional(scale.d, 0, rlength, scale.scale_factor)
                    .with_baseline(lambda_bg),
            );
            combined = Some(match combined {
                None => p,
                Some(prev) => macs_peaks::over_two_pv_array_track(&prev, &p),
            });
        }
        let Some(t) = combined else { continue };
        let mut runs = Vec::new();
        let mut prev = 0u64;
        for r in t.runs() {
            // stop where the pairing would have stopped
            if r.end > limit {
                break;
            }
            runs.push((prev, r.end, f64::from(r.value)));
            prev = r.end;
        }
        out.insert(name, runs);
    }
    out
}

/// Compare two bedGraphs as the piecewise-constant functions they describe.
///
/// `__write_bedGraph_for_a_chromosome` only emits a line when the value changes by
/// more than `1e-5`, so the file is a *coalesced* view of the underlying array: the
/// number and position of run boundaries is an artefact of that threshold, not
/// data. Comparing run-by-run therefore reports spurious differences — a boundary
/// that lands on one side of the threshold but not the other.
///
/// The right comparison is of the function each file defines. Both are step
/// functions over intervals, so the union of all interval starts is a refinement of
/// both, and on that refinement the two must agree.
///
/// Returns `(intervals compared, worst absolute deviation, problems)`.
fn compare_graph(
    label: &str,
    want: &BedGraph,
    got: &BedGraph,
) -> (usize, f64, Vec<String>, Vec<String>) {
    let mut problems = Vec::new();
    let mut notes = Vec::new();
    let mut worst = 0.0f64;
    let mut n = 0usize;
    let mut keys: Vec<&String> = want.keys().chain(got.keys()).collect();
    keys.sort();
    keys.dedup();
    for k in keys {
        let w = want.get(k).map(|v| v.as_slice()).unwrap_or(&[]);
        let g = got.get(k).map(|v| v.as_slice()).unwrap_or(&[]);
        if w.is_empty() || g.is_empty() {
            problems.push(format!(
                "{label}/{k}: present in only one side (upstream {} runs, rust {})",
                w.len(),
                g.len()
            ));
            continue;
        }
        // Each side's domain is `[first_start, last_end)`. Upstream writes the
        // *paired* arrays (F35), so its file stops where the shorter of the
        // treatment and control tracks stops, while an untruncated port keeps
        // going. Comparing the union of the start positions would then evaluate
        // upstream's final value at positions past the end of its own file --
        // which reads as "one extra read" but is purely a support difference.
        //
        // So: compare over the intersection of the two domains, and report the
        // support difference once, separately, instead of per-position.
        let w_end = w.last().map(|(_, e, _)| *e).unwrap_or(0);
        let g_end = g.last().map(|(_, e, _)| *e).unwrap_or(0);
        let common_end = w_end.min(g_end);
        if w_end != g_end {
            notes.push(format!(
                "{label}/{k}: domain differs -- upstream ends at {w_end}, rust ends at {g_end} \
                 (comparing [0,{common_end}))"
            ));
        }

        // every start from either side refines both step functions
        let mut bounds: Vec<u64> = w.iter().map(|(s, _, _)| *s).collect();
        bounds.extend(g.iter().map(|(s, _, _)| *s));
        bounds.retain(|b| *b < common_end);
        bounds.sort_unstable();
        bounds.dedup();

        for (i, st) in bounds.iter().enumerate() {
            let next = bounds.get(i + 1).copied().unwrap_or(*st);
            if next == *st {
                continue;
            }
            let wv = value_at(w, *st);
            let gv = value_at(g, *st);
            let (Some(wv), Some(gv)) = (wv, gv) else {
                problems.push(format!("{label}/{k} @{st}: covered by one side only"));
                continue;
            };
            n += 1;
            let d = (wv - gv).abs();
            if d > worst {
                worst = d;
            }
            if d > 1e-4 {
                problems.push(format!(
                    "{label}/{k} [{st},{next}): value {wv:.5} vs {gv:.5} (diff {d:.2e})"
                ));
            }
        }
    }
    (n, worst, problems, notes)
}

/// The value a run list reports at `pos`, if `pos` is covered by it.
fn value_at(runs: &[(u64, u64, f64)], pos: u64) -> Option<f64> {
    runs.iter()
        .find(|(s, e, _)| pos >= *s && pos < *e)
        .map(|(_, _, v)| *v)
        .or_else(|| {
            // a zero-length query at the very end still has a value upstream
            runs.last()
                .filter(|(_, e, _)| pos == *e)
                .map(|(_, _, v)| *v)
        })
}
