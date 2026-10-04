//! `macs3-rs randsample`: down-sample a BED alignment by percentage or count.
//!
//! Delegates the sampling to `macs_track::sample_percent`, which reproduces
//! upstream's `np.random.shuffle`-based `FWTrack.sample_percent` (including
//! NumPy's global-seed stream via `macs_stats::NumpyRng`), and writes the
//! surviving tags in upstream's `print_to_bed` layout -- the same writer
//! `filterdup` uses.
//!
//! Port notes (`randsample_cmd.py:28-78`):
//! * `-n/--number` overrides `-p/--percentage`, and upstream rejects
//!   `number > total` before writing anything.
//! * `--seed` >= 0 seeds NumPy's global RNG; negative uses the global state.

use std::path::{Path, PathBuf};

use macs_core::{MacsError, Result, Strand};
use macs_track::SingleEndTrack;

use crate::Options;

fn load_bed(path: &Path) -> Result<SingleEndTrack> {
    super::input::load_single_end_bed(path)
}

pub fn randsample(o: &Options) -> Result<()> {
    let ifile = o
        .get("ifile")
        .ok_or_else(|| MacsError::InvalidParameter("-i/--ifile is required".into()))?;
    let format = o.get("format").unwrap_or("BED").to_uppercase();
    // Paired-end modes sample fragments (`randsample_cmd.py:45-72`). FRAG is still
    // rejected: it needs barcode/count handling on top of fragment sampling.
    if format == "BAMPE" || format == "BEDPE" {
        return randsample_pe(o, Path::new(ifile), &format);
    }
    if format == "FRAG" {
        return Err(MacsError::InvalidParameter(
            "paired-end input is not yet supported".into(),
        ));
    }
    let mut track = if format == "BAM" {
        super::input::load_single_end_bam(Path::new(ifile))?.0
    } else if format == "SAM" {
        super::input::load_single_end_sam(Path::new(ifile))?.0
    } else {
        load_bed(Path::new(ifile))?
    };
    let t0 = track.total();

    // percentage (default upstream) and optional -n override
    let mut percent = o.float("percentage").unwrap_or(100.0);
    if let Some(num) = o.int("number") {
        if num < 0 {
            return Err(MacsError::InvalidParameter("-n must be >= 0".into()));
        }
        if num as u64 > t0 {
            return Err(MacsError::InvalidParameter(format!(
                "number requested ({num}) is bigger than total tags ({t0})"
            )));
        }
        percent = (num as f64) / (t0 as f64) * 100.0;
    }
    // F199: `-s/--tsize` is optional. `randsample_cmd.py:87-89` falls back to the
    // parser's own estimate whenever the user omitted it, and that is the normal path
    // -- rejecting the invocation instead made every `randsample` without an explicit
    // `-s` a hard error where upstream succeeds.
    let fw = match o.int("tsize") {
        Some(v) if v > 0 => v as i32,
        Some(v) => {
            return Err(MacsError::InvalidParameter(format!(
                "--tsize must be > 0 (got {v})"
            )))
        }
        None => match macs_io::detect_tsize(Path::new(ifile), &format)? {
            Some(v) => v,
            None => {
                return Err(MacsError::InvalidParameter(
                    "could not estimate tag size from the input; pass -s/--tsize".into(),
                ))
            }
        },
    };
    eprintln!("randsample: tag size is determined as {fw} bps");
    let seed = o.int("seed").unwrap_or(-1);

    let kept = macs_track::sample_percent(&mut track, percent / 100.0, seed)?;
    eprintln!("randsample: {t0} tags, {kept} kept ({percent:.2}%)");

    // print_to_bed
    let ofile = o.get("outputfile").unwrap_or("stdout");
    let mut out = String::new();
    let pos = track.positions();
    for chrom in pos.chroms_sorted() {
        let name = String::from_utf8_lossy(track.genome().name(chrom));
        for &p in pos.strand(chrom, Strand::Plus) {
            out.push_str(&format!("{name}\t{p}\t{}\t.\t.\t+\n", p + fw as u64));
        }
        for &p in pos.strand(chrom, Strand::Minus) {
            // Upstream writes negative starts verbatim (`chrIV -11 39`); it does
            // not clamp at zero. `saturating_sub` hid every minus-strand tag within
            // `fw` of the contig start behind a `0`, diverging on real data.
            let lo = p as i64 - i64::from(fw);
            out.push_str(&format!("{name}\t{lo}\t{p}\t.\t.\t-\n"));
        }
    }
    if ofile == "stdout" {
        print!("{out}");
    } else {
        let outdir = PathBuf::from(o.get("outdir").unwrap_or("."));
        std::fs::create_dir_all(&outdir)?;
        std::fs::write(outdir.join(ofile), out)?;
    }
    Ok(())
}

/// Down-sample paired-end fragments (BEDPE/BAMPE).
///
/// `randsample_cmd.py:45-72` via `PETrack.sample_percent` (`PairedEndTrack.py:550-580`):
/// `-n` becomes a fraction of the total, the quota per chromosome is
/// `uint(round(len * f32(percent), 5))`, and the shuffle uses
/// `RandomState(MT19937(SeedSequence(seed)))` in sorted chromosome order.
/// Single-chromosome output is byte-identical; multi-chromosome inherits the
/// single-chromosome stream per chromosome in sorted order (upstream iterates a
/// set, so only sorted content can match there).
fn randsample_pe(o: &crate::Options, path: &Path, format: &str) -> Result<()> {
    use std::collections::BTreeMap;

    let mut frags: BTreeMap<Vec<u8>, Vec<(u64, u64)>> = BTreeMap::new();
    if format == "BAMPE" {
        let (records, _) = macs_io::bam::bampe_fragments(path)?;
        for fr in &records {
            frags
                .entry(fr.chrom.clone())
                .or_default()
                .push((u64::from(fr.start), u64::from(fr.start + fr.len)));
        }
    } else {
        use std::io::BufRead;
        let mut r = std::io::BufReader::new(macs_io::open_maybe_gzip(path)?);
        let mut line = Vec::new();
        loop {
            line.clear();
            if r.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            let Some(rec) = macs_io::parse_bedpe_line(&line)? else {
                continue;
            };
            if rec.chrom.is_empty() || rec.left < 0 || rec.right < rec.left {
                continue;
            }
            frags
                .entry(rec.chrom.clone())
                .or_default()
                .push((rec.left as u64, rec.right as u64));
        }
    }
    let t0: u64 = frags.values().map(|v| v.len() as u64).sum();

    // `-n` becomes a fraction; upstream keeps it 0-1 for the quota math.
    let fraction = if let Some(num) = o.int("number") {
        if num < 0 {
            return Err(MacsError::InvalidParameter("-n must be >= 0".into()));
        }
        if num as u64 > t0 {
            return Err(MacsError::InvalidParameter(format!(
                "number {num} exceeds total fragments {t0}"
            )));
        }
        num as f64 / t0 as f64
    } else {
        o.float("percentage").unwrap_or(100.0) / 100.0
    };
    let percent = fraction * 100.0;
    let seed = o.int("seed").unwrap_or(-1);

    // `rs = RandomState(MT19937(SeedSequence(seed)))`, sorted chromosomes, shuffle
    // each independently, retain the quota, restore coordinate order.
    //
    // Fragments are sorted by `(start, end)` before shuffling: upstream's
    // `finalize()` sorts, so the shuffle input is sorted regardless of file order.
    // A coordinate-sorted BAM streams fragments out of `(start, end)` order (sorted
    // by leftmost position only), and shuffling stream order gave a different subset
    // than shuffling sorted order here.
    for spans in frags.values_mut() {
        spans.sort_unstable();
    }
    let mut rng = macs_stats::randomstate_from_seed_sequence(seed as u64);
    type Spans = Vec<(u64, u64)>;
    let mut sampled: Vec<(Vec<u8>, Spans)> = Vec::new();
    let mut kept_total = 0u64;
    for (chrom, spans) in &frags {
        // `uint(round(len * f32(percent), 5))`: the percent narrows to f32 on entry
        // to `sample_percent`, but the multiply promotes back to f64.
        let pct_f32 = (fraction) as f32;
        let num = macs_track::retained_count(spans.len() as f64 * f64::from(pct_f32));
        let mut idx: Vec<usize> = (0..spans.len()).collect();
        rng.shuffle(&mut idx);
        let mut kept: Vec<(u64, u64)> = idx.into_iter().take(num).map(|i| spans[i]).collect();
        kept.sort_unstable();
        kept_total += kept.len() as u64;
        sampled.push((chrom.clone(), kept));
    }
    eprintln!("randsample: {t0} fragments, {kept_total} kept ({percent:.2}%)");

    let ofile = o.get("outputfile").unwrap_or("stdout");
    let mut out = String::new();
    for (chrom, spans) in &sampled {
        let name = String::from_utf8_lossy(chrom);
        for (s, e) in spans {
            out.push_str(&format!("{name}\t{s}\t{e}\n"));
        }
    }
    if ofile == "stdout" {
        print!("{out}");
    } else {
        let outdir = std::path::PathBuf::from(o.get("outdir").unwrap_or("."));
        std::fs::create_dir_all(&outdir)?;
        std::fs::write(outdir.join(ofile), out)?;
    }
    Ok(())
}
