//! `macs3-rs pileup`: pileup a BED alignment into a bedGraph.
//!
//! Delegates to `macs_pileup::pileup_from_positions`, whose
//! `SingleEndParams::directional` / `bidirectional` shifts already reproduce
//! upstream's `PileupV2.pileup_and_write_se` five/three-shift cases exactly
//! (F-era G5 work, 1776 bit-exact vectors), then writes the bedGraph with
//! upstream's `%.5f` format and post-scale baseline floor.
//!
//! Port notes (`pileup_cmd.py:26-89`, `PileupV2.py:1134-1187`):
//! * `--both-direction` doubles `--extsize` and switches to the symmetric
//!   (`five = d/2`, `three = d - d/2`) extension; otherwise the directional
//!   (`five = 0`, `three = d`) one.
//! * `halfextension` is **False** here, which is what selects those shift pairs
//!   -- the `--half-extension` variant is a different code path.
//! * Values are floored at the baseline *after* scaling, and written with
//!   `%.5f` (`_write_pv_to_bedGraph`).
//! * Single-end chromosomes are written in file order: upstream's
//!   `pileup_and_write_se` iterates `chrlengths.keys()` (insertion order).
//!   Paired-end chromosomes are written **sorted**: `pileup_and_write_pe` iterates
//!   `sorted(chrlengths.keys())`. The two functions disagree, so the order is
//!   per-branch, not global.

use std::path::{Path, PathBuf};

use macs_core::{MacsError, Result, Strand};
use macs_pileup::SingleEndParams;

use crate::Options;

pub fn pileup(o: &Options) -> Result<()> {
    let ifiles = super::input::input_files(o)?;
    let format = o.get("format").unwrap_or("BED").to_uppercase();
    // Paired-end modes pile up fragments, not 5' ends (`pileup_cmd.py:50-68`,
    // `PileupV2.pileup_and_write_pe`). FRAG is still rejected: it needs the barcode
    // subset and count weighting, which live in the callpeak/hmmratac FRAG paths.
    if format == "BAMPE" || format == "BEDPE" {
        return pileup_pe(o, &ifiles, &format);
    }
    // FRAG piles up count-weighted fragments, optionally subset by `--barcodes`
    // (`pileup_cmd.py:54-68`). Sorted chromosome order, like the other PE path.
    if format == "FRAG" {
        return pileup_frag(o, &ifiles);
    }
    // SE BAM streams without an index; everything else text goes through the BED
    // loader (which sniffs gzip). Upstream selects the parser by `-f` directly.
    let track = super::input::load_tag_files(&ifiles, &format, true)?.0;

    let extsize = o.int("extsize").unwrap_or(0);
    if extsize <= 0 {
        return Err(MacsError::InvalidParameter("--extsize must be > 0".into()));
    }
    let both = o.flag("bothdirection");
    // upstream passes `extsize * 2` and symmetric shifts for --both-direction
    let d = if both { extsize * 2 } else { extsize };
    let rlength: u32 = u32::MAX / 2;
    let baseline = 0.0f32;

    let ofile = o.get("outputfile").unwrap_or("pileup.bedGraph");
    let outdir = PathBuf::from(o.get("outdir").unwrap_or("."));
    std::fs::create_dir_all(&outdir)?;
    // Stream to a buffered writer rather than accumulating the whole bedGraph in a
    // `String`. The genome-wide text is hundreds of MB on a real run (5 M reads ->
    // ~500 MB here) and the whole thing was held resident alongside the read track,
    // which made `pileup` the one command using *more* memory than upstream
    // (313 MB against 138 MB). Row order and formatting are unchanged.
    use std::io::Write as _;
    let mut out = std::io::BufWriter::new(std::fs::File::create(outdir.join(ofile))?);

    let pos = track.positions();
    // File order, not sorted order: upstream's `pileup_and_write_se` iterates
    // `chrlengths.keys()`, which follows the track's insertion (file) order, so a
    // multi-contig input like `contigs50k.bed.gz` emits `unitig0, unitig1, unitig2, ...`.
    // Sorting lexicographically emitted `unitig0, unitig1, unitig10, ...` instead --
    // same rows, wrong order. (`callpeak` is different: it sorts explicitly, so it
    // keeps `chroms_sorted`.)
    for chrom in pos.chroms() {
        let params = if both {
            SingleEndParams::bidirectional(d, 0, rlength, 1.0).with_baseline(baseline)
        } else {
            SingleEndParams::directional(d, 0, rlength, 1.0).with_baseline(baseline)
        };
        let t = macs_pileup::pileup_from_positions(
            chrom,
            pos.strand(chrom, Strand::Plus),
            pos.strand(chrom, Strand::Minus),
            &params,
        );
        let name = String::from_utf8_lossy(track.genome().name(chrom));
        let mut pre = 0u32;
        for r in t.runs() {
            let v = r.value.max(baseline);
            writeln!(out, "{name}\t{pre}\t{}\t{v:.5}", r.end)?;
            pre = r.end;
        }
    }
    out.flush()?;
    Ok(())
}

/// Pile up paired-end fragments (BEDPE/BAMPE) into a bedGraph.
///
/// `pileup_and_write_pe` (`PileupV2.py:1190-1217`): stream the fragments, pile up
/// per chromosome with `pileup_from_LR` (unweighted; BEDPE/BAMPE carry no counts),
/// and write in **sorted** chromosome order -- unlike the single-end path, which
/// uses file order.
fn pileup_pe(o: &Options, paths: &[String], format: &str) -> Result<()> {
    use std::collections::BTreeMap;

    // Collect (start, end) per chromosome name, preserving input order within each.
    let mut frags: BTreeMap<Vec<u8>, Vec<(u32, u32)>> = BTreeMap::new();
    let track = super::input::load_fragment_files(paths, format)?;
    for chrom in track.chroms() {
        let name = track.genome().name(chrom).to_vec();
        for fragment in track.frags(chrom) {
            frags
                .entry(name.clone())
                .or_default()
                .push((fragment.start, fragment.end));
        }
    }

    let ofile = o.get("outputfile").unwrap_or("pileup.bedGraph");
    let outdir = PathBuf::from(o.get("outdir").unwrap_or("."));
    std::fs::create_dir_all(&outdir)?;
    // Buffered writer, not a whole-file `String`: see the SE path's note.
    use std::io::Write as _;
    let mut out = std::io::BufWriter::new(std::fs::File::create(outdir.join(ofile))?);
    // BTreeMap iterates keys in byte-lexicographic order, which is upstream's
    // `sorted(chrlengths.keys())`.
    for (name, spans) in &frags {
        let cname = String::from_utf8_lossy(name);
        // `rlength` clamps fragment ends; upstream passes the chromosome length.
        // Without a genome file the length is unknown, so use u64::MAX (no clamp),
        // which matches for all fragments in practice.
        let t = macs_pileup::pileup_from_fragments(
            // `pileup_from_fragments` needs a ChromId, but we only have names here.
            // Use a dummy id; the track is only read for its runs, never interned.
            macs_core::ChromId(0),
            spans,
            u32::MAX / 2,
            1.0,
            0.0,
        );
        let mut pre = 0u32;
        for r in t.runs() {
            writeln!(out, "{cname}\t{pre}\t{}\t{:.5}", r.end, r.value)?;
            pre = r.end;
        }
    }
    out.flush()?;
    Ok(())
}

/// Pile up FRAG fragments (single-cell), count-weighted.
///
/// `pileup_cmd.py:54-68` via `pileup_and_write_pe`: subset by `--barcodes` when given,
/// cap counts at `--max-count` when positive, pile up with per-fragment weights
/// (`pileup_from_LRC`), sorted chromosomes. A `--max-count` of 0 means "keep all
/// counts" (upstream: "If this is set as 0, MACS3 will behave as the default setting
/// to keep all counts").
fn pileup_frag(o: &Options, paths: &[String]) -> Result<()> {
    use std::collections::BTreeMap;

    // `--barcodes` allow-list, read once.
    let barcodes: Option<std::collections::HashSet<Vec<u8>>> = match o.get("barcodefile") {
        Some(bf) if !bf.is_empty() => {
            use std::io::BufRead;
            let mut set = std::collections::HashSet::new();
            for line in macs_io::open_maybe_gzip(Path::new(bf))?.split(b'\n') {
                let line = line?;
                let mut end = line.len();
                while end > 0 && matches!(line[end - 1], b'\n' | b'\r' | b' ' | b'\t') {
                    end -= 1;
                }
                if end > 0 {
                    set.insert(line[..end].to_vec());
                }
            }
            Some(set)
        }
        _ => None,
    };
    let max_count = o.int("maxcount").unwrap_or(0).max(0) as u32;

    let mut frags: BTreeMap<Vec<u8>, Vec<(u32, u32, u32)>> = BTreeMap::new();
    for path in paths {
        use std::io::BufRead;
        let mut r = std::io::BufReader::new(macs_io::open_maybe_gzip(Path::new(path))?);
        let mut line = Vec::new();
        loop {
            line.clear();
            if r.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            let Some(rec) = macs_io::parse_frag_line(&line)? else {
                continue;
            };
            if rec.chrom.is_empty() || rec.left < 0 || rec.right < rec.left {
                continue;
            }
            if let Some(allow) = &barcodes {
                match &rec.barcode {
                    Some(bc) if allow.contains(bc) => {}
                    _ => continue,
                }
            }
            let mut count = rec.count.unwrap_or(1);
            if max_count > 0 {
                count = count.min(max_count);
            }
            if count == 0 {
                continue;
            }
            frags.entry(rec.chrom.clone()).or_default().push((
                rec.left as u32,
                rec.right as u32,
                count,
            ));
        }
    }

    let ofile = o.get("outputfile").unwrap_or("pileup.bedGraph");
    let outdir = PathBuf::from(o.get("outdir").unwrap_or("."));
    std::fs::create_dir_all(&outdir)?;
    use std::io::Write as _;
    let mut out = std::io::BufWriter::new(std::fs::File::create(outdir.join(ofile))?);
    for (name, spans) in &frags {
        let cname = String::from_utf8_lossy(name);
        let weighted: Vec<(u32, u32, f32)> =
            spans.iter().map(|(s, e, c)| (*s, *e, *c as f32)).collect();
        let t = macs_pileup::pileup_from_weighted_fragments(
            macs_core::ChromId(0),
            &weighted,
            u32::MAX / 2,
            1.0,
            0.0,
        );
        let mut pre = 0u32;
        for r in t.runs() {
            writeln!(out, "{cname}\t{pre}\t{}\t{:.5}", r.end, r.value)?;
            pre = r.end;
        }
    }
    out.flush()?;
    Ok(())
}
