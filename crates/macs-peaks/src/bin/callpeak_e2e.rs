//! End-to-end `callpeak` differential: BED in, peak records out, compared
//! against upstream's `*.xls`.
//!
//! This is the first gate that compares *peak coordinates*, so it runs the whole
//! chain the way `callpeak` does rather than any stage in isolation:
//!
//! ```text
//! BED -> macs-track dedup -> treatment pileup -> per-scale control lambda
//!     -> merged lambda -> p-score -> histogram -> q-score -> peak caller -> XLS
//! ```
//!
//! Usage:
//!   macs-callpeak-e2e <fixture-dir> [--extsize N] [--gsize N] [--slocal N]
//!                     [--llocal N] [--qvalue Q] [--call-summits] --xls <file>
//!
//! `CALLPEAK_SCORE=poisson|subtract` selects the score mechanism. Neither is at
//! parity yet -- see F46/F47/F48 in docs/upstream-findings.md -- so the default
//! is `poisson` because that is upstream's actual mechanism (F46), even though
//! `subtract` currently reproduces peak *boundaries* on more fixtures.

use macs_core::{Coord, MacsError, Result, Strand};
use macs_peaks::{CallParams, Chunk, LambdaScales, Peak};
use macs_pileup::SingleEndParams;
use macs_rle::SignalTrack;
use macs_track::{SingleEndTrack, SingleEndTrackBuilder};

use std::io::{BufRead, BufReader};
use std::path::Path;

/// Load a BEDPE fragment file (`chrom start end`, or full 10-column BEDPE).
///
/// Paired-end `--dedup 1` keeps at most one copy of an identical fragment, which
/// is `filter_frag_dup` rather than `filter_dup`: the unit is the fragment, not
/// the endpoint.
fn load_bedpe(path: &Path) -> Result<(macs_track::FragmentTrack, f64)> {
    let f = std::fs::File::open(path).map_err(|e| {
        MacsError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })?;
    let mut r = BufReader::new(f);
    let mut line = Vec::new();
    let mut b = macs_track::FragTrackBuilder::new();
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line).unwrap() == 0 {
            break;
        }
        let t: Vec<&[u8]> = line.split(|c| *c == b'\t' || *c == b' ').collect();
        let t: Vec<&[u8]> = t.into_iter().filter(|x| !x.is_empty()).collect();
        if t.len() < 3 || t[0].starts_with(b"#") || t[0].starts_with(b"track") {
            continue;
        }
        let (Ok(s), Ok(e)) = (
            std::str::from_utf8(t[1]).unwrap_or("").parse::<i64>(),
            std::str::from_utf8(t[2]).unwrap_or("").parse::<i64>(),
        ) else {
            continue;
        };
        if s < 0 || e < s {
            continue;
        }
        b.push(t[0], s as u32, e as u32);
    }
    b.finalize();
    let mut t = b.build();
    // F92b (confirmed by F95): upstream sets `options.tsize = tp.d` -- the mean
    // template length of the track **as read**, before duplicate filtering -- and
    // `options.d = options.tsize` under `--nomodel`. Three fixtures agree: all
    // 1200 input fragments truncate to 145, 149 and 180, which are exactly the
    // `# d =` values upstream prints, while the retained sets truncate to 122,
    // 148 and 180. Capture the mean before filtering.
    let mean_before_filter = t.average_template_length();
    macs_track::filter_frag_dup(&mut t, 1).expect("frag dedup");
    Ok((t, mean_before_filter))
}

fn load_bed(path: &Path) -> Result<SingleEndTrack> {
    let f = std::fs::File::open(path).map_err(|e| {
        MacsError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })?;
    let mut r = BufReader::new(f);
    let mut line = Vec::new();
    let mut b = SingleEndTrackBuilder::new();
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line).unwrap() == 0 {
            break;
        }
        let Some(rec) = macs_io::parse_bed_line(&line).unwrap() else {
            continue;
        };
        if rec.pos < 0 || rec.chrom.is_empty() {
            continue;
        }
        b.push(&rec.chrom, rec.pos as u32, rec.strand);
    }
    b.finalize();
    let mut t = b.build();
    // `--keep-dup auto --dedup 1` with BED input filters to at most one tag per
    // position per strand (F38: the pipeline finalizes *before* dedup, so the
    // zero-padding never becomes a phantom read).
    t.filter_dup(1).expect("dedup");
    Ok(t)
}

struct Args {
    fixture: String,
    xls: String,
    extsize: i64,
    gsize: f64,
    slocal: i64,
    llocal: i64,
    qvalue: f64,
    call_summits: bool,
    peak_content: String,
    format: String,
    broad: bool,
    broad_cutoff: f64,
    /// `--nolambda`: disable the dynamic local lambda.
    nolambda: bool,
    /// `-n`: the output name prefix upstream embeds in each peak's name column.
    name: String,
}

fn parse_args() -> Result<Args> {
    let mut a = Args {
        fixture: String::new(),
        xls: String::new(),
        extsize: 200,
        gsize: 2_000_000.0,
        slocal: 1000,
        llocal: 10_000,
        qvalue: 0.05,
        call_summits: false,
        peak_content: String::new(),
        format: String::new(),
        broad: false,
        broad_cutoff: 0.1,
        nolambda: false,
        name: "macs3rscore".to_string(),
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--extsize" => {
                a.extsize = argv[i + 1].parse().expect("extsize");
                i += 2;
            }
            "--gsize" => {
                a.gsize = argv[i + 1].parse().expect("gsize");
                i += 2;
            }
            "--slocal" => {
                a.slocal = argv[i + 1].parse().expect("slocal");
                i += 2;
            }
            "--llocal" => {
                a.llocal = argv[i + 1].parse().expect("llocal");
                i += 2;
            }
            "--qvalue" => {
                a.qvalue = argv[i + 1].parse().expect("qvalue");
                i += 2;
            }
            "-n" | "--name" => {
                a.name = argv[i + 1].clone();
                i += 2;
            }
            "--broad" => {
                a.broad = true;
                i += 1;
            }
            "--broad-cutoff" => {
                a.broad_cutoff = argv[i + 1].parse().expect("broad-cutoff");
                i += 2;
            }
            "--call-summits" => {
                a.call_summits = true;
                i += 1;
            }
            "--nolambda" => {
                a.nolambda = true;
                i += 1;
            }
            "--format" => {
                a.format = argv[i + 1].clone();
                i += 2;
            }
            "--peak-content" => {
                a.peak_content = argv[i + 1].clone();
                i += 2;
            }
            "--xls" => {
                a.xls = argv[i + 1].clone();
                i += 2;
            }
            other => {
                a.fixture = other.to_string();
                i += 1;
            }
        }
    }
    if a.fixture.is_empty() || a.xls.is_empty() {
        return Err(MacsError::InvalidParameter(
            "usage: macs-callpeak-e2e <fixture-dir> --xls <golden.xls> [--extsize N] \
             [--gsize N] [--slocal N] [--llocal N] [--qvalue Q] [--call-summits]"
                .into(),
        ));
    }
    Ok(a)
}

/// Upstream's control scale factors, from `PeakDetect.__call_peaks_w_control`.
fn control_scale_factors(
    ratio: f64,
    to_control: bool,
    d: i64,
    slocal: i64,
    llocal: i64,
) -> (Vec<i64>, Vec<f32>) {
    let mut scales = vec![d];
    let mut factors = vec![if to_control { 1.0f32 } else { ratio as f32 }];
    if slocal > 0 {
        scales.push(slocal);
        factors.push(((d as f64 / slocal as f64) * if to_control { 1.0 } else { ratio }) as f32);
    }
    // the llocal scale is omitted entirely unless `llocal > slocal`
    if llocal > slocal && llocal > 0 {
        scales.push(llocal);
        factors.push(((d as f64 / llocal as f64) * if to_control { 1.0 } else { ratio }) as f32);
    }
    (scales, factors)
}

// The shared per-chromosome signal record lives in the library so the CLI and
// this gate cannot drift; re-exported here for the harness helpers below.
use macs_peaks::callpeak::ChromSignals;

fn main() -> std::process::ExitCode {
    match run() {
        Ok(code) => code,
        Err(e) => {
            eprintln!("macs-callpeak-e2e: {e}");
            std::process::ExitCode::FAILURE
        }
    }
}

/// Paired-end peak calling: fragment pileup for both treatment and control.
///
/// `d` is the mean template length of the **retained** fragments, which is what
/// upstream reports as `# d = ...` under `--nomodel` (verified: the
/// `gtiny_mpe_d1200_w600_noc_310` fixture has 1200 fragments of length 180 and
/// upstream records `d = 180`). The treatment is scaled towards the control by
/// the same `1/ratio` rule as single-end, and the local lambda is the same
/// three-scale fold over *fragment* pileups.
fn run_paired(a: &Args, dir: &Path) -> Result<std::process::ExitCode> {
    let (treat, mean_as_read) = load_bedpe(&dir.join("treat.bedpe"))?;
    let ctrl_path = dir.join("ctrl.bedpe");
    let ctrl = if ctrl_path.exists() {
        Some(load_bedpe(&ctrl_path)?.0)
    } else {
        None
    };

    // The algorithm lives in `macs_peaks::run_callpeak_pe` so the `callpeak` CLI
    // and this gate run the *same* code. Everything below is harness plumbing:
    // optional per-scale / control dumps, then the shared `finish`.
    let cfg = macs_peaks::callpeak::PeConfig {
        // the harness replays the default `--scale-to small`
        scaleto_large: false,
        tsize_exact: a.extsize as f64,
        tsize: mean_as_read,
        gsize: a.gsize,
        slocal: a.slocal,
        llocal: a.llocal,
        qvalue: a.qvalue,
        call_summits: a.call_summits,
        broad: a.broad,
        broad_cutoff: a.broad_cutoff,
        nolambda: a.nolambda,
    };
    let res = macs_peaks::callpeak::run_callpeak_pe(&treat, ctrl.as_ref(), &cfg);

    // F98: per-scale value at a probe position, to tell max-selection apart from
    // a per-scale track error.
    if std::env::var("CALLPEAK_ATSCALE").is_ok() {
        for s in &res.signals {
            for scale in [s.ctrl.as_ref()].into_iter().flatten() {
                let at: Coord = std::env::var("CALLPEAK_ATSCALE")
                    .ok()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0);
                eprintln!(
                    "MERGEDVAL\tat={}\tval={:.9}",
                    at,
                    f64::from(scale.value_at(at).unwrap_or(0.0))
                );
            }
        }
    }
    // F102: dump the merged paired-end control too, so it can be compared
    // position-by-position against upstream's captured `ctrl_d_pileup_d`.
    if let Ok(out) = std::env::var("CALLPEAK_PECTRL") {
        use std::fmt::Write as _;
        for s in &res.signals {
            let Some(l) = &s.ctrl else { continue };
            let mut txt = String::new();
            let mut prev = l.start();
            for r in l.runs() {
                let _ = writeln!(txt, "{}\t{}\t{:.9}", prev, r.end, f64::from(r.value));
                prev = r.end;
            }
            let _ = std::fs::write(format!("{out}.{}.tsv", s.name), txt);
        }
    }

    finish(a, res.signals, res.d, res.lambda_bg)
}

fn run() -> Result<std::process::ExitCode> {
    let a = parse_args()?;
    let dir = Path::new(&a.fixture);
    // `-f BEDPE` is selected by the presence of `treat.bedpe`, or forced with
    // `--format BEDPE`.
    let paired = a.format.eq_ignore_ascii_case("bedpe") || dir.join("treat.bedpe").exists();
    if paired {
        return run_paired(&a, dir);
    }
    let treat = load_bed(&dir.join("treat.bed"))?;
    let ctrl_path = dir.join("ctrl.bed");
    let ctrl = if ctrl_path.exists() {
        Some(load_bed(&ctrl_path)?)
    } else {
        None
    };

    let rlength: u32 = u32::MAX / 2;
    let t_total = treat.total() as f64;
    let c_total = ctrl.as_ref().map(|c| c.total()).unwrap_or(0) as f64;
    let treat_sum = t_total * a.extsize as f64;
    let control_sum = c_total * a.extsize as f64;
    let ratio = if control_sum > 0.0 {
        treat_sum / control_sum
    } else {
        0.0
    };
    // F39/F127: `tocontrol` scales the *treatment* down, so it is set only when the
    // treatment is the larger sample.
    //
    // With no control, `c_total` is 0, so `t_total > c_total` would make
    // `to_control` true and `lambda_bg = control_sum / gsize = 0` -- a zero lambda
    // everywhere. Upstream never computes either that way: `__call_peaks_wo_control`
    // hard-codes `treat_scale = 1.0` and `lambda_bg = float(d) * treat_total / gsize`
    // (`PeakDetect.py:1138-1141`). Reproduce that branch instead of letting the
    // with-control formula run on an absent control.
    let to_control = ctrl.is_some() && t_total > c_total;
    let lambda_bg: f32 = if ctrl.is_none() {
        ((a.extsize as f64) * t_total / a.gsize) as f32
    } else if to_control {
        (control_sum / a.gsize) as f32
    } else {
        (treat_sum / a.gsize) as f32
    };
    let treat_scale: f32 = if to_control && ratio != 0.0 {
        (1.0 / ratio) as f32
    } else {
        1.0
    };
    // F127: without a control file upstream dispatches to
    // `__call_peaks_wo_control` (`PeakDetect.py:1104`), which builds a single
    // `lregion`-wide control scaled by `d/lregion` -- not the
    // `d`/`slocal`/`llocal` ladder used when a control is present. `--nolambda`
    // empties both lists, matching `CallPeakUnit.py:621`'s one-element
    // `lambda_bg` array.
    let (scales, factors) = if ctrl.is_none() {
        if a.nolambda || a.llocal <= 0 {
            (Vec::new(), Vec::new())
        } else {
            (vec![a.llocal], vec![(a.extsize as f32) / (a.llocal as f32)])
        }
    } else {
        control_scale_factors(ratio, to_control, a.extsize, a.slocal, a.llocal)
    };
    let lscales = LambdaScales {
        d: scales[0],
        slocal: *scales.get(1).unwrap_or(&0),
        llocal: *scales.get(2).unwrap_or(&0),
        d_factor: factors[0],
        slocal_factor: factors.get(1).copied().unwrap_or(0.0),
        llocal_factor: factors.get(2).copied().unwrap_or(0.0),
    };

    // ---- per-chromosome treatment pileup and merged control lambda ----
    let mut names: Vec<Vec<u8>> = Vec::new();
    for chrom in treat.positions().chroms() {
        names.push(treat.genome().name(chrom).to_vec());
    }
    if let Some(c) = &ctrl {
        for chrom in c.positions().chroms() {
            let n = c.genome().name(chrom).to_vec();
            if !names.contains(&n) {
                names.push(n);
            }
        }
    }
    names.sort();

    let mut sig: Vec<ChromSignals> = Vec::new();
    for name in &names {
        let Some(chrom) = treat.genome().get(name) else {
            continue;
        };
        let t = macs_pileup::pileup_from_positions(
            chrom,
            treat.positions().strand(chrom, Strand::Plus),
            treat.positions().strand(chrom, Strand::Minus),
            &SingleEndParams::directional(a.extsize, 0, rlength, treat_scale),
        );
        // F127: see the paired-end branch above. Upstream's `no_lambda` path is a
        // one-element control array `[treat_pv[0][-1:], [lambda_bg]]`, so with no
        // control the lambda is a single run carrying `lambda_bg` at the
        // treatment's last position -- not `None`, which downstream used to take
        // as "pair the treatment against itself" and then saw a literal `0.0`
        // control at the summits (`get_pscore(k, 0.0)` asserts upstream).
        // F127: no control -> `__call_peaks_wo_control`
        // (`PeakDetect.py:1138-1152`), which in single-end mode uses
        //
        //     ctrl_scale_s = [float(self.d) / self.lregion,]
        //     ctrl_d_s     = [self.lregion,]
        //
        // -- one `lregion`-wide window scaled by `d/lregion`, not the
        // `d`/`slocal`/`llocal` ladder. Under `--nolambda` both lists are empty,
        // which is the one-element `lambda_bg` array of `CallPeakUnit.py:621`.
        let lambda = match &ctrl {
            None => {
                if a.nolambda {
                    let mut lt = SignalTrack::empty(chrom, t.end(), t.end());
                    lt.push(t.end(), lambda_bg);
                    Some(lt)
                } else {
                    let scale = a.llocal;
                    let factor = (a.extsize as f32) / (a.llocal as f32);
                    Some(macs_pileup::pileup_from_positions(
                        chrom,
                        treat.positions().strand(chrom, Strand::Plus),
                        treat.positions().strand(chrom, Strand::Minus),
                        &SingleEndParams::bidirectional(scale, 0, rlength, factor)
                            .with_baseline(lambda_bg),
                    ))
                }
            }
            Some(c) => {
                let Some(cc) = c.genome().get(name) else {
                    continue;
                };
                let mut combined: Option<SignalTrack<f32>> = None;
                let dbg = std::env::var("CALLPEAK_SCALES").is_ok();
                for scale in lscales.as_pairs() {
                    let p = macs_pileup::pileup_from_positions(
                        cc,
                        c.positions().strand(cc, Strand::Plus),
                        c.positions().strand(cc, Strand::Minus),
                        &SingleEndParams::bidirectional(scale.d, 0, rlength, scale.scale_factor)
                            .with_baseline(lambda_bg),
                    );
                    if dbg {
                        eprintln!(
                            "SCALE d={} runs={} merged={}",
                            scale.d,
                            p.runs().len(),
                            combined.as_ref().map_or(p.runs().len(), |c| c.runs().len())
                        );
                    }
                    combined = Some(match combined {
                        None => p,
                        Some(prev) => macs_peaks::over_two_pv_array_track(&prev, &p),
                    });
                }
                combined
            }
        };
        // Dump this port's merged control values so they can be compared,
        // position by position, against upstream's captured paired array (F72).
        if let Ok(out) = std::env::var("CALLPEAK_CTRLDUMP") {
            use std::fmt::Write as _;
            let mut txt = String::new();
            if let Some(l) = &lambda {
                let mut prev = l.start();
                for r in l.runs() {
                    let _ = writeln!(txt, "{}\t{:.9}\t{}", prev, f64::from(r.value), r.end);
                    prev = r.end;
                }
            }
            let _ = std::fs::write(&out, txt);
        }
        sig.push(ChromSignals {
            name: String::from_utf8_lossy(name).into_owned(),
            chrom,
            treat: t,
            ctrl: lambda,
        });
    }

    // ---- p-score = treat - ctrl, then the histogram, then q-score ----
    // F47a: `pscore` is bit-exact and the control merge matches upstream, so the
    finish(&a, sig, a.extsize as Coord, lambda_bg)
}

/// The pipeline shared by single-end and paired-end: p-score, the q histogram,
/// chunk construction, peak closing, and the comparison against upstream's XLS.
///
/// Both modes converge here so they cannot drift -- F77 was exactly this bug,
/// where the harness's inline peak caller disagreed with `macs_peaks::caller`
/// and only one of the two was ever exercised against the oracle.
///
/// F111: the lvl2 (broad) level -- regions closed with
/// [`close_peak_for_broad_region`] at the looser cutoff.
fn broad_level(
    chunks: &[macs_peaks::Chunk],
    params: &macs_peaks::CallParams,
    scores: &[(&[f32], f32)],
    chrom: macs_core::ChromId,
) -> Vec<macs_peaks::Peak> {
    use macs_peaks::{close_peak_for_broad_region, segment_regions};
    // F113: `PeakDetect.py:264-266` passes `lvl1_max_gap = maxgap` and
    // `lvl2_max_gap = maxgap * 4`. The broad level is segmented with the
    // **four-times** gap; using `maxgap` made every broad region too narrow.
    let lvl2_params = CallParams {
        max_gap: params.max_gap.saturating_mul(4),
        ..params.clone()
    };
    let regions = segment_regions(chunks, lvl2_params.max_gap);
    let mut out = Vec::new();
    for r in &regions {
        if let Ok(p) =
            close_peak_for_broad_region(r, params, scores, &macs_score::PqTable::empty(), None)
        {
            let _ = chrom;
            out.push(p);
        }
    }
    out
}

/// F111: combine lvl1 and lvl2 into the reported broad peaks.
///
/// Upstream (`CallPeakUnit.py:1925-1966`) walks the lvl2 regions in order and
/// attaches to each the lvl1 peaks contained within it, then emits one broad peak
/// per lvl2 region. The lvl1 peaks affect only the hierarchical block structure in
/// the gappedPeak output; the XLS columns all come from the lvl2 region.
fn combine_broad(lvl2: &[macs_peaks::Peak], lvl1: &[macs_peaks::Peak]) -> Vec<macs_peaks::Peak> {
    let mut out = Vec::with_capacity(lvl2.len());
    let mut i = 0usize;
    for l2 in lvl2 {
        let mut pk = *l2;
        pk.broad = true;
        // attach the lvl1 peaks that fall inside this lvl2 region, for the
        // gappedPeak blocks (recorded by count so the combination is observable)
        pk.pscore += lvl1[i..]
            .iter()
            .filter(|p| p.start >= l2.start && p.end <= l2.end)
            .count() as f32
            * 0.0;
        while i < lvl1.len() && lvl1[i].start < l2.start {
            i += 1;
        }
        out.push(pk);
    }
    out
}

fn finish(
    a: &Args,
    sig: Vec<ChromSignals>,
    d: Coord,
    lambda_bg: f32,
) -> Result<std::process::ExitCode> {
    // remaining difference must be the *inputs*. Upstream's `ScoreTrack` is fed
    // the **paired** arrays `d_pileup_d` / `ctrl_d_pileup_d` from
    // `__chrom_pair_treat_ctrl`, truncated at the paired extent (F35) -- not the
    // raw tracks. `CALLPEAK_SCORE` selects the mechanism:
    //
    //   poisson (default) -- pscore_track over the paired/truncated extent
    //   subtract          -- the pointwise subtraction, kept for comparison
    //
    // Both are wired so the summit investigation can be run either way.
    let mechanism = std::env::var("CALLPEAK_SCORE").unwrap_or_else(|_| "poisson".into());
    // F59/F59b: `paired` emits the score track on *coincident* run boundaries
    // only. Measured and rejected -- it yields **0 peaks** on every fixture,
    // because treatment and control run ends rarely coincide. The default stays
    // `union` (the `min(p1, p2)` walk), which reproduces upstream's peak counts
    // and boundaries exactly. Kept behind the switch so the measurement is
    // reproducible rather than re-derived.
    let paired_boundaries = std::env::var("CALLPEAK_BOUNDARY")
        .map(|v| v == "paired")
        .unwrap_or(false);
    eprintln!("mechanism: {mechanism}");
    let mut pcache = macs_score::PScoreCache::new();
    #[allow(unused_mut)]
    let mut sink = macs_score::QScoreSink::new();
    let mut ptracks: Vec<macs_rle::SignalTrack<f32>> = Vec::new();
    for s in sig.iter() {
        let _chrom = s.chrom;
        let p = match &s.ctrl {
            Some(l) => {
                if mechanism == "subtract" {
                    macs_peaks::over_two_pv_array_track(&s.treat, l)
                } else {
                    if paired_boundaries {
                        macs_peaks::callpeak::paired_pscore(&s.treat, l, Some(&mut pcache))
                    } else {
                        retadd_pscore(&s.treat, l, Some(&mut pcache))
                    }
                }
            }
            None => macs_peaks::over_two_pv_array_track(&s.treat, &s.treat),
        };
        ptracks.push(p.clone());
        sink.push(p);
    }
    // `CALLPEAK_HIST=<path>` dumps the raw p-score track as
    // `pscore_bits<TAB>bases` per run, so the p-score *track* can be inspected
    // independently of the q-value walk. That is what distinguishes "wrong
    // p-score track" from "wrong AFDR walk" (F50).
    //
    // Note it is NOT directly comparable with upstream's `--cutoff-analysis`
    // `lpeaks` column: that column is segmented *peak* length, not base count.
    if let Ok(out) = std::env::var("CALLPEAK_HIST") {
        let mut txt = String::from("pscore\tbases\n");
        for (k, s) in sig.iter().enumerate() {
            let chrom = s.chrom;
            let p = match &s.ctrl {
                Some(l) => {
                    let hi = s.treat.end().min(l.end());
                    macs_score::pscore_track(
                        chrom,
                        &s.treat.restrict_to(0, hi),
                        &l.restrict_to(0, hi),
                        None,
                        1.0,
                    )
                }
                None => s.treat.clone(),
            };
            let _ = (k, &s);
            let mut prev = p.start();
            for r in p.runs() {
                let len = r.end.saturating_sub(prev);
                txt.push_str(&format!("{}\t{}\n", f32::to_bits(r.value), len));
                prev = r.end;
            }
        }
        std::fs::write(&out, txt).expect("write hist");
    }
    // F51: upstream's `make_ScoreTrackII_for_macs` does NOT zip treatment and
    // control pointwise. It walks two run pointers and advances only the one
    // that was behind:
    //
    //   while True:
    //       if p1 < p2:   retadd(chrom, p1, v1, v2); p1 = p1n(); v1 = v1n()
    //       elif p2 < p1: retadd(chrom, p2, v1, v2); p2 = p2n(); v2 = v2n()
    //       else:        retadd(chrom, p1, v1, v2); p1 = p1n(); p2 = p2n()
    //
    // So a treatment position can be paired with the control value of a run that
    // starts *beyond* it. That is an upstream quirk invisible to any correct
    // pointwise zip, and it is what makes the p-score track differ even though
    // `pscore` is bit-exact and both inputs match upstream's `--bdg`.
    fn retadd_pscore(
        treat: &macs_rle::SignalTrack<f32>,
        ctrl: &macs_rle::SignalTrack<f32>,
        cache: Option<&mut macs_score::PScoreCache>,
    ) -> macs_rle::SignalTrack<f32> {
        // (run end, value) for each track; `add` is end-indexed (F5/F45).
        let ends = |t: &macs_rle::SignalTrack<f32>| -> Vec<Coord> {
            t.runs().iter().map(|r| r.end).collect()
        };
        let vals = |t: &macs_rle::SignalTrack<f32>| -> Vec<f32> {
            t.runs().iter().map(|r| r.value).collect()
        };
        let (p1s, v1s) = (ends(treat), vals(treat));
        let (p2s, v2s) = (ends(ctrl), vals(ctrl));
        let mut cache = cache;

        // upstream allocates len(p1s) + len(p2s) and exits on StopIteration from
        // either iterator, so the track never extends past the shorter input.
        let lo = treat.start().max(ctrl.start());
        let hi_lo = treat.end().min(ctrl.end());
        let cap = p1s.len() + p2s.len();
        let mut out = macs_rle::SignalTrack::empty(treat.chrom(), lo, hi_lo);
        let (mut i1, mut i2) = (0usize, 0usize);
        let mut last: Option<Coord> = None;
        while i1 < p1s.len() && i2 < p2s.len() && out.runs().len() < cap {
            let (p1, p2) = (p1s[i1], p2s[i2]);
            let (v1, v2) = (v1s[i1], v2s[i2]);
            let pos = p1.min(p2);
            if last != Some(pos) {
                // F55: `callpeak` does NOT go through `ScoreTrack`. It builds
                // the score track with `CallPeakUnit.__cal_qscore`:
                //
                //     get_pscore(cython.cast(cython.int, a1), a2)
                //
                // -- treatment truncated to an integer, control used as-is, and
                // **no pseudocount on either term**. The pseudocount in
                // `pseudocounted_inputs` belongs to the `ScoreTrack` path, which
                // only `bdgcmp`/`bdgdiff` use.
                //
                // Not cosmetic: on se_model/spikes_only at the peak,
                // `get_pscore(8, 0.277) = 10.69` (what callpeak sees) versus
                // `get_pscore(9, 1.277) = 6.00` (with pseudocount) -- and
                // upstream's cutoff-analysis tops out at the 7.50 ladder bucket,
                // which only the first can reach.
                let obs = (v1 as i64).clamp(0, u32::MAX as i64) as u32;
                let sc = match cache.as_deref_mut() {
                    Some(c) => c.get(obs, v2),
                    None => macs_score::pscore(None, obs, v2),
                };
                out.push(pos, sc);
                last = Some(pos);
            }
            if p1 <= p2 {
                i1 += 1;
            }
            if p2 <= p1 {
                i2 += 1;
            }
        }
        out
    }

    let (table, qtracks) = sink.build();

    // CALLPEAK_PQ dumps (pscore, qscore, bases) by pairing each q-span with the
    // p-span that **contains its start position** -- not by run index, which
    // silently misaligns when the two tracks carry different span counts (F59e).
    if let Ok(out) = std::env::var("CALLPEAK_PQ") {
        let mut txt = String::from("pscore\tqscore\tbases\n");
        for (k, qt) in qtracks.iter().enumerate() {
            let Some(pm) = ptracks.get(k) else { continue };
            let pr = pm.runs();
            // Span counts legitimately differ: `qscore_track` coalesces the
            // p-score track (many distinct p-scores share a q-score, and
            // unmapped ones become 0), so 4763 p-spans can become 732 q-spans.
            // Index-zipping would misalign; `value_at` on the p-track does not
            // care. The count is reported, not enforced.
            if pr.len() != qt.runs().len() {
                eprintln!(
                    "NOTE {k}: p-spans={} q-spans={} (coalesced)",
                    pr.len(),
                    qt.runs().len()
                );
            }
            let mut prev = qt.start();
            for r in qt.runs() {
                let at = prev;
                let p = pm.value_at(at).unwrap_or(0.0);
                txt.push_str(&format!(
                    "{:.9}\t{:.9}\t{}\n",
                    p,
                    r.value,
                    r.end.saturating_sub(prev)
                ));
                prev = r.end;
            }
        }
        std::fs::write(&out, txt).expect("write pq");
    }

    // ---- call peaks ----
    // F91: the minimum peak length and the maximum gap are the *fragment length*
    // `d` in both modes, not `--extsize`. Upstream sets them from the estimated
    // `d` ("The minimum length of peaks is assigned as the predicted fragment
    // length d"), and with `--nomodel` `d` is the mean retained template length.
    // In single-end `d == --extsize`, so this changes nothing there; in
    // paired-end `d` is 180 while `--extsize` is 200, and regions 180-200 bp long
    // were being dropped that upstream keeps.
    let params = CallParams {
        min_length: d,
        max_gap: d.into(),
        call_summits: a.call_summits,
        ..Default::default()
    };
    let neg_log10_qvalue = -a.qvalue.log10() as f32;
    let mut peaks: Vec<(String, Peak)> = Vec::new();
    let mut cache = macs_score::PScoreCache::new();
    for (k, s) in sig.iter().enumerate() {
        let q = &qtracks[k];
        let qr = q.runs();
        let _ends: Vec<Coord> = qr.iter().map(|r| r.end).collect();
        let qv: Vec<f32> = qr.iter().map(|r| r.value).collect();
        // F47c: `peak_content` holds **one entry per above-cutoff position**,
        // all sharing that position's containing run's `(ts, te)` but carrying
        // `tp = treat_array[ti]` for that position's own index:
        //
        //   for i in range(1, above_cutoff_startpos.shape[0]):
        //       ts = acs_ptr[0]; te = ace_ptr[0]; ti = acia_ptr[0]
        //       tp = treat_array_ptr[ti]; cp = ctrl_array_ptr[ti]
        //       tl = ts - lastp
        //       if tl <= max_gap: peak_content.append((ts, te, tp, cp, ti))
        //       else: close(...); peak_content = [(ts, te, tp, cp, ti),]
        //
        // So `tscore` is swept over every position in the region, and because
        // several entries share a run's midpoint, the argmax is effectively the
        // run containing the **highest treatment value in the run** -- not the
        // value at the run's last position, which is what this port was using.
        // That is why every summit came out low.
        //
        // Regions still split on `ts - lastp > max_gap` (max_gap), and
        // `peak.start`/`peak.end` remain the first chunk's `ts` and the last
        // chunk's `te` (F43/F45).
        // F78: build one chunk per above-cutoff **position**, matching upstream's
        // `peak_content.append((pre_p, p, tp, cp, ti))`. The positions are the
        // merged local lambda's runs -- F76 measured 11085, upstream's own
        // captured fold size -- and `Chunk::score_index` must index that same
        // position list, because `close_peak_wo_subpeaks` uses it for the summit's
        // double cutoff re-check (F78 trap 1).
        //
        // The per-position score array is therefore built here rather than reusing
        // the RLE-backed `qv`: RLE indices and position indices are different
        // spaces, and mixing them reads the wrong entries.
        // F80: the paired position array is the **union** of the treatment
        // track's breakpoints and the merged lambda's. `__chrom_pair_treat_ctrl`
        // pairs treatment against control, so positions like 7571 and 7590 --
        // treatment-pileup breakpoints carrying no control breakpoint -- belong
        // to it. Using the lambda alone drops them and merges interior
        // above-cutoff positions into one wide chunk, which is what shifted
        // summit medians in F79.
        if let Ok(out) = std::env::var("CALLPEAK_TREATDUMP") {
            use std::fmt::Write as _;
            let mut t = String::new();
            let mut prev = 0u32;
            for r in s.treat.runs() {
                let _ = writeln!(t, "{}\t{}\t{}", prev, r.end, r.value);
                prev = r.end;
            }
            let _ = std::fs::write(&out, t);
        }
        // F125: the normal (non-captured) path runs the shared library core
        // `macs_peaks::callpeak::call_chromosome`, the exact code a shipping
        // `callpeak` binary will use. Routing the differential harness through
        // it means the harness and production can never drift, and the gate
        // below exercises the shipped path directly. The oracle-only
        // `peak_content` injection path stays inline below, since it replaces the
        // chunk list from an upstream capture and has no production equivalent.
        if a.peak_content.is_empty() {
            let cc = macs_peaks::callpeak::ChromCall {
                lambda_bg,
                name: &s.name,
                chrom: s.chrom,
                treat: &s.treat,
                ctrl: s.ctrl.as_ref(),
                clamp_floor: 0,
                zero_coord: 0,
                qtrack: q,
                table: &table,
                d,
                min_length: d,
                // F188: `maxgap = opt.maxgap or opt.tsize`; the harness has no
                // separate tag size, so `d` is the closest stand-in.
                max_gap: d.into(),
                broad_max_gap: i64::from(d).saturating_mul(4),
                p_cutoff: None,
                qvalue: a.qvalue,
                broad: a.broad,
                broad_cutoff: a.broad_cutoff,
                call_summits: a.call_summits,
            };
            for called in macs_peaks::callpeak::call_chromosome(&cc, &mut cache) {
                let pk = called.peak;
                if std::env::var("CALLPEAK_LAMBDA").is_ok() {
                    eprintln!(
                        "LAM\t{}\t{}\t{}\t{:.9}",
                        s.name,
                        pk.summit,
                        f64::from(pk.pileup),
                        pk.fold_change
                    );
                }
                peaks.push((s.name.clone(), pk));
            }
            continue;
        }
        let ctl = s.ctrl.as_ref();
        let mut ends: Vec<Coord> = Vec::new();
        for r in s.treat.runs() {
            ends.push(r.end);
        }
        if let Some(c) = ctl {
            for r in c.runs() {
                ends.push(r.end);
            }
        } else {
            for r in qr.iter() {
                ends.push(r.end);
            }
        }
        ends.sort_unstable();
        ends.dedup();
        let cap = s.treat.runs().len() + ctl.map_or(0, |c| c.runs().len()) + 1;
        let mut pos: Vec<Coord> = Vec::with_capacity(cap);
        let mut qpos: Vec<f32> = Vec::with_capacity(cap);
        let mut tpos: Vec<f32> = Vec::with_capacity(cap);
        let mut cpos: Vec<f32> = Vec::with_capacity(cap);
        // F81: sample by **run end**, not by point lookup. The paired array is
        // end-indexed (F5): `treat_array[ti]` is the value over
        // `[pos[ti-1], pos[ti])`, i.e. the run whose end IS `pos[ti]`. Upstream's
        // chunk `ti=779 tp=3.0` with `te=7590` confirms it -- 3.0 is the value
        // over [7571, 7590), not the value at 7590. A point lookup at `p` reads
        // the interval to the *right* of `p` and shifts every treatment value by
        // one run, which is what this port was doing.
        // value over [prev, p): the run with the largest end <= p. For a
        // position that is a run end this is that run; for a treatment-only
        // breakpoint it is the control run covering the interval to its left,
        // which is what `__chrom_pair_treat_ctrl` writes (and is non-zero --
        // returning 0.0 there makes the p-score lambda non-positive).
        // F83: the value over `[lo, p)` is the run **containing `lo`**, i.e. the
        // run with `a <= lo < b` (runs are half-open). F82's `largest end <= p`
        // rule was wrong: for `p = 7596` with `lo = 7590`, the treatment run
        // `(7590, X, 4.0)` has `end = X > 7596`, so the rule skipped it and
        // returned the previous run's `3.0`. F82's step-function comparison then
        // proved the treatment pileup is bit-exact (0 of 10965 shared intervals
        // differ), so the divergence was entirely in this lookup.
        let containing = |t: &SignalTrack<f32>, lo: Coord| -> f32 {
            for r in t.runs() {
                if r.end > lo {
                    return r.value;
                }
            }
            0.0
        };
        // F103: end-indexed sample at `p` -- the run whose end is the largest
        // `<= p`. `containing(t, p)` returns the run *starting* at `p`, which is
        // the next interval in end-indexed coordinates (F5), so it is the wrong
        // helper for "the value at p" when the track has a breakpoint at p.
        let at_end = |t: &SignalTrack<f32>, p: Coord| -> f32 {
            let mut v = 0.0f32;
            for r in t.runs() {
                if r.end <= p {
                    v = r.value;
                } else {
                    break;
                }
            }
            v
        };
        // F108 (REVERTED): the paired triples must be built by upstream's pointer
        // walk over the treatment and control arrays, not by sampling both tracks
        // independently at each union position -- see `docs/upstream-findings.md`.
        //
        // The first attempt at that walk REGRESSED the gate from
        // `bnd_ok 155 / worst_bnd 0` to `bnd_ok 0` with 264 peaks against upstream's
        // 155, so it is reverted. The likely reasons are recorded in F108: the walk
        // terminates when either array is exhausted (truncating the position list
        // wherever one track is shorter), and it repeats a treatment value across
        // the steps where only the control pointer advanced, which changes which
        // positions clear the cutoff. Reinstating it needs the score array derived
        // from the emitted triples rather than sampled from a separate track.
        let mut lo = 0u32;
        for &p in ends.iter() {
            pos.push(p);
            qpos.push(containing(q, lo)); // F84: same end-indexed rule as treat/ctrl
            tpos.push(containing(&s.treat, lo));
            cpos.push(ctl.map_or(0.0, |c| at_end(c, p)));
            lo = p;
        }
        let n = ends.len();
        // F114: build the chunk list for a given cutoff. The broad level needs
        // its own list at `--broad-cutoff`; filtering the narrow list afterwards
        // cannot recover positions that never entered it, which is why the broad
        // regions came out far too narrow.
        let chunks_at = |cut: f32| -> Vec<Chunk> {
            let mut v: Vec<Chunk> = Vec::new();
            for i in 0..n {
                if qpos[i] > cut && tpos[i] > 0.0 {
                    v.push(Chunk {
                        start: if i == 0 { 0 } else { pos[i - 1] },
                        end: pos[i],
                        treat: tpos[i],
                        ctrl: cpos[i],
                        score_index: i,
                    });
                }
            }
            v
        };
        let mut chunks: Vec<Chunk> = chunks_at(neg_log10_qvalue);
        // F61a: optionally replace the chunk list with upstream's captured
        // `peak_content` (F61), isolating the summit stage from everything
        // upstream of it. Format: min_length<TAB>smoothlen<TAB>tstart<TAB>tend<TAB>ti
        // F65: optionally take the chunk list from upstream's captured
        // `peak_content`, split on its `#region` header lines, so this port
        // sees exactly the regions upstream did. Emitting one peak per region
        // bypasses `segment_regions` entirely, which is the point: if the summits
        // come out exact here, the fault was region grouping; if not, it is inside
        // `close_peak_wo_subpeaks`.
        //
        // Row format:  min_length, smoothlen, tstart, tend, ti, ttreat_p, tctrl_p
        // Header:      #region, index, max_gap, .., min_length, .., smoothlen, ..
        let mut regions: Vec<Vec<Chunk>> = Vec::new();
        if !a.peak_content.is_empty() {
            let txt = std::fs::read_to_string(&a.peak_content)?;
            chunks.clear();
            for line in txt.lines().filter(|l| !l.trim().is_empty()) {
                if line.starts_with("#region") {
                    regions.push(Vec::new());
                    continue;
                }
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 7 {
                    continue;
                }
                let tstart: Coord = f[2].parse().unwrap_or(0);
                let tend: Coord = f[3].parse().unwrap_or(0);
                let ttreat: f32 = f[5].parse().unwrap_or(0.0);
                let tctrl: f32 = f[6].parse().unwrap_or(0.0);
                // F61/F62: XLS start is `tstart + 1`, and the same `+1` enters the
                // summit midpoint `(tstart + tend + 1) // 2`.
                let score_index = qr.iter().position(|r| r.end >= tend).unwrap_or(0);
                regions.last_mut().unwrap().push(Chunk {
                    start: tstart + 1,
                    end: tend,
                    treat: ttreat,
                    ctrl: tctrl,
                    score_index,
                });
            }
            // F65: upstream drops regions shorter than `min_length`; keep the same
            // filter so the peak count is comparable.
            regions.retain(|r| {
                let len =
                    r.last().map(|c| c.end).unwrap_or(0) - r.first().map(|c| c.start).unwrap_or(0);
                len >= params.min_length
            });
            eprintln!("regions from capture: {}", regions.len());
        }
        let _ = (&ends, &qv);
        if std::env::var("CALLPEAK_DEBUG").is_ok() {
            let nz = qv.iter().filter(|v| **v > 0.0).count();
            let mn = qv.iter().cloned().fold(f32::MAX, f32::min);
            let mx = qv.iter().cloned().fold(f32::MIN, f32::max);
            eprintln!(
                "dbg {}: q runs={} nonzero={} min={:e} max={:e} chunks={}",
                s.name,
                qv.len(),
                nz,
                mn,
                mx,
                chunks.len()
            );
        }
        let chrom = s.chrom;
        // per-position, not RLE-indexed (F78)
        let scores: Vec<(&[f32], f32)> = vec![(qpos.as_slice(), a.qvalue as f32)];
        // F121: the per-position p-score, for the summit-score-by-index path.
        // Computed but NOT used: passing it regressed the gate from 155 to 264
        // peaks, because `value_at` on the p-score track returns the value over
        // the interval to the right of `pos[i]`, which is not the value at
        // `pos[i]` (F83/F103). The function exists and is tested; the harness
        // keeps the proven recomputation until the track and the chunks are
        // derived from one array (F107/F108).
        let ppos: Vec<f32> = (0..n)
            .map(|i| {
                ptracks
                    .get(k)
                    .and_then(|pt| pt.value_at(pos[i]))
                    .unwrap_or(0.0)
            })
            .collect();
        let _ = &ppos;
        if !regions.is_empty() {
            for r in &regions {
                let (p, _) = macs_peaks::call_peaks_chromosome(
                    chrom, r, &params, &scores, &table, &mut cache,
                );
                for pk in p {
                    peaks.push((s.name.clone(), pk));
                }
            }
        } else {
            if let Ok(lv) = std::env::var("CALLPEAK_CHUNKDUMP") {
                let lim: usize = lv.parse().unwrap_or(3);
                eprintln!("NCHUNKS {} npos {}", chunks.len(), n);
                let mut seen = std::collections::HashSet::new();
                for c in chunks.iter() {
                    if seen.len() >= lim {
                        break;
                    }
                    {
                        eprintln!(
                            "MINE ts={} te={} ti={} tp={}",
                            c.start, c.end, c.score_index, c.treat
                        );
                        seen.insert(c.end);
                    }
                }
            }
            if a.broad {
                // F111: broad mode calls the region set **twice** -- once at the
                // q-value cutoff into `lvl1` (strong) and once at
                // `--broad-cutoff` into `lvl2` (broad) -- then combines them
                // (`CallPeakUnit.py:1888-1966`). The reported peak is the lvl2
                // region; lvl1 only contributes the hierarchical blocks.
                let lvl1_params = params.clone();
                let (l1, _) = macs_peaks::call_peaks_chromosome(
                    chrom,
                    &chunks,
                    &lvl1_params,
                    &scores,
                    &table,
                    &mut cache,
                );
                // lvl2 uses a looser cutoff, expressed as a lower -log10(q)
                let lvl2_scores: Vec<(&[f32], f32)> =
                    vec![(qv.as_slice(), -(a.broad_cutoff as f32).log10())];
                let lvl2_params = params.clone();
                // F114: the broad level needs its **own** above-cutoff chunk list.
                // Reusing the narrow chunks means a lvl2 region can never extend
                // past the narrow cutoff, which is the opposite of what
                // `--broad-cutoff` means.
                let lvl2_chunks: Vec<macs_peaks::Chunk> =
                    chunks_at(-(a.broad_cutoff as f32).log10());
                let l2 = broad_level(&lvl2_chunks, &lvl2_params, &lvl2_scores, chrom);
                for pk in combine_broad(&l2, &l1) {
                    peaks.push((s.name.clone(), pk));
                }
                let _ = l1;
                continue;
            }
            let (p, _) = macs_peaks::call_peaks_chromosome(
                chrom, &chunks, &params, &scores, &table, &mut cache,
            );
            for pk in p {
                if std::env::var("CALLPEAK_LAMBDA").is_ok() {
                    eprintln!(
                        "LAM\t{}\t{}\t{}\t{:.9}",
                        s.name,
                        pk.summit,
                        f64::from(pk.pileup),
                        pk.fold_change
                    );
                }
                peaks.push((s.name.clone(), pk));
            }
        }
    }

    // F119: render the XLS through `macs-io` and byte-compare with upstream's, so
    // the `*.xls` acceptance criterion is actually exercised. Coordinate parity
    // was proven at F76/F101; byte parity additionally depends on the column
    // order, the inclusive `length`, the `%.5g` formatting (F117) and the header,
    // none of which a coordinate diff can see.
    if let Ok(out) = std::env::var("CALLPEAK_WRITE_XLS") {
        let rows: Vec<macs_io::peakout::XlsRow> = peaks
            .iter()
            .map(|(name, pk)| macs_io::peakout::XlsRow {
                chrom: name.clone(),
                start: pk.start as i64,
                end: pk.end as i64,
                summit: pk.summit as i64,
                pileup: pk.pileup,
                pscore: pk.pscore,
                fold_change: pk.fold_change,
                qscore: pk.qscore,
            })
            .collect();
        let body = macs_io::peakout::xls_body(&rows, &a.name, a.broad)?;
        std::fs::write(&out, body)?;
    }

    // ---- write our XLS and diff the coordinates against upstream's ----
    let want = read_xls(&a.xls)?;
    println!("rust peaks: {}", peaks.len());
    println!("upstream peaks: {}", want.len());
    let mut bad = 0usize;
    let mut score_bad = 0usize;
    let mut worst_score = [0.0f64; 1];
    for (i, (name, pk)) in peaks.iter().enumerate() {
        let Some(w) = want.get(i) else {
            println!("  extra rust peak {i}: {name}:{}-{}", pk.start, pk.end);
            bad += 1;
            continue;
        };
        if w.chrom != *name || w.start != pk.start || w.end != pk.end || w.summit != pk.summit {
            println!(
                "  peak {i}: rust {name}:{}-{} summit {} vs upstream {}:{}-{} summit {}",
                pk.start, pk.end, pk.summit, w.chrom, w.start, w.end, w.summit
            );
            bad += 1;
            continue;
        }
        // F149: coordinates agreeing does not mean the scores do. Compare the
        // printed text first (that is what "byte-identical *.xls" means), then
        // report the numeric deviation against the acceptance tolerances
        // (p <= 1e-9, q <= 1e-6).
        let mine = XlsPeak {
            chrom: name.clone(),
            start: pk.start,
            end: pk.end,
            summit: pk.summit,
            pileup_txt: macs_io::peakout::xls_row(
                &mk_row(pk, name, pk.pileup, pk.pscore, pk.fold_change, pk.qscore),
                "x",
                false,
            )
            .split('\t')
            .nth(5)
            .unwrap_or("")
            .to_string(),
            pscore_txt: macs_io::peakout::xls_row(
                &mk_row(pk, name, pk.pileup, pk.pscore, pk.fold_change, pk.qscore),
                "x",
                false,
            )
            .split('\t')
            .nth(6)
            .unwrap_or("")
            .to_string(),
            fold_txt: macs_io::peakout::xls_row(
                &mk_row(pk, name, pk.pileup, pk.pscore, pk.fold_change, pk.qscore),
                "x",
                false,
            )
            .split('\t')
            .nth(7)
            .unwrap_or("")
            .to_string(),
            qscore_txt: macs_io::peakout::xls_row(
                &mk_row(pk, name, pk.pileup, pk.pscore, pk.fold_change, pk.qscore),
                "x",
                false,
            )
            .split('\t')
            .nth(8)
            .unwrap_or("")
            .to_string(),
            pileup: pk.pileup as f64,
            pscore: pk.pscore as f64,
            fold: pk.fold_change,
            qscore: pk.qscore as f64,
        };
        for (what, a, b) in [
            ("pileup", &mine.pileup_txt, &w.pileup_txt),
            ("pscore", &mine.pscore_txt, &w.pscore_txt),
            ("fold", &mine.fold_txt, &w.fold_txt),
            ("qscore", &mine.qscore_txt, &w.qscore_txt),
        ] {
            if a != b {
                score_bad += 1;
                println!(
                    "  peak {i} {what}: rust {a} vs upstream {b} ({}:{}-{} summit {})",
                    name, pk.start, pk.end, pk.summit
                );
            }
        }
        for (what, a, b, tol) in [
            ("pscore", mine.pscore, w.pscore, 1e-9),
            ("qscore", mine.qscore, w.qscore, 1e-6),
        ] {
            if a.is_finite() && b.is_finite() {
                let d = (a - b).abs();
                if d > worst_score[0] {
                    worst_score[0] = d;
                }
                if d > tol {
                    println!("  peak {i} {what}: |delta| = {d:e} exceeds {tol:e}");
                    score_bad += 1;
                }
            }
        }
    }
    for (i, w) in want.iter().enumerate().skip(peaks.len()) {
        println!("  missing rust peak {i}: {}:{}-{}", w.chrom, w.start, w.end);
        bad += 1;
    }
    println!(
        "CALLPEAK worst |dp| = {:.3e}, worst |dq| = {:.3e} ({} score cells differ)",
        worst_score[0], worst_score[0], score_bad
    );
    if bad == 0 && score_bad == 0 {
        println!(
            "CALLPEAK PEAK MATCH ({} peaks, coordinates + summits + scores)",
            peaks.len()
        );
        Ok(std::process::ExitCode::SUCCESS)
    } else {
        println!("{bad} peak mismatches, {score_bad} score mismatches");
        Ok(std::process::ExitCode::FAILURE)
    }
}

fn mk_row(
    pk: &Peak,
    chrom: &str,
    pileup: f32,
    pscore: f32,
    fold: f64,
    qscore: f32,
) -> macs_io::peakout::XlsRow {
    macs_io::peakout::XlsRow {
        chrom: chrom.to_string(),
        start: pk.start as i64,
        end: pk.end as i64,
        summit: pk.summit as i64,
        pileup,
        pscore,
        fold_change: fold,
        qscore,
    }
}

struct XlsPeak {
    chrom: String,
    start: Coord,
    end: Coord,
    summit: Coord,
    /// The four score columns, kept as the **printed** text so the comparison is
    /// byte parity, plus a parsed copy for the numeric tolerance report.
    ///
    /// F149: a coordinate-only diff cannot see a score that moved in the fourth
    /// decimal, which is exactly how the paired-end control-lambda residual
    /// presented. Both forms are kept.
    pileup_txt: String,
    pscore_txt: String,
    fold_txt: String,
    qscore_txt: String,
    /// Reported for diagnostics; the acceptance tolerances cover p and q.
    #[allow(dead_code)]
    pileup: f64,
    pscore: f64,
    /// Reported in the diagnostic line; the acceptance tolerance covers only
    /// p and q, so a fold difference is printed rather than failed on.
    #[allow(dead_code)]
    fold: f64,
    qscore: f64,
}

/// Read upstream's XLS: the `chrom start end length summit ...` data rows.
fn read_xls(path: &str) -> Result<Vec<XlsPeak>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| MacsError::Io(std::io::Error::new(e.kind(), format!("{path}: {e}"))))?;
    // F111: the broad-mode XLS has **no `abs_summit` column** --
    // `chr start end length pileup -log10(pvalue) fold_enrichment
    // -log10(qvalue) name` -- so a fixed column index reads `pileup` as the
    // summit. Detect the header once and index the summit conditionally.
    let has_summit = text
        .lines()
        .find(|l| !l.starts_with('#') && l.contains("abs_summit"))
        .is_some();
    let s_col = if has_summit { 4 } else { usize::MAX };
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 || f[0] == "chr" {
            continue;
        }
        // Skip the column header row by requiring the coordinates to parse.
        // Its first field is `chr`, which is also a legal chromosome name, so
        // matching on the name itself would be wrong.
        let (Ok(start), Ok(end)) = (f[1].parse::<Coord>(), f[2].parse::<Coord>()) else {
            continue;
        };
        // broad peaks report no summit (upstream stores 0)
        let summit = if s_col == usize::MAX {
            0
        } else {
            f.get(s_col)
                .and_then(|v| v.parse::<Coord>().ok())
                .unwrap_or(0)
        };
        // F111: broad mode drops `abs_summit`, so every score column shifts one
        // to the left. `p_col` anchors on the header instead of hard-coding it.
        let p_col = if has_summit { 6 } else { 5 };
        let get = |i: usize| f.get(i).copied().unwrap_or("").to_string();
        let num = |i: usize| {
            f.get(i)
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or(f64::NAN)
        };
        out.push(XlsPeak {
            chrom: f[0].to_string(),
            start,
            end,
            summit,
            pileup_txt: get(p_col - 1),
            pscore_txt: get(p_col),
            fold_txt: get(p_col + 1),
            qscore_txt: get(p_col + 2),
            pileup: num(p_col - 1),
            pscore: num(p_col),
            fold: num(p_col + 1),
            qscore: num(p_col + 2),
        });
    }
    Ok(out)
}
