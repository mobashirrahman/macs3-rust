//! `macs3-rs predictd`: estimate the fragment length `d` and emit the
//! `*_model.r` plotting script.
//!
//! Delegates to `macs_model::PeakModel` (the `PeakModel.py` port) and writes the
//! R script via `model2r_script`. Upstream also emits the R script from
//! `callpeak`'s model path, so both share this.

use std::path::{Path, PathBuf};

use macs_core::{MacsError, Result};
use macs_model::{model2r_script, ModelOptions, PeakModel};
use macs_track::SingleEndTrack;

use crate::Options;

fn load_bed(path: &Path) -> Result<SingleEndTrack> {
    super::input::load_single_end_bed(path)
}

/// Load a BEDPE/BAMPE fragment track, for the paired-end branch of `predictd`.
fn load_frag_track(path: &Path) -> Result<macs_track::FragmentTrack> {
    use std::io::BufRead;
    // Sniffing reader, so a gzipped BEDPE is decompressed: upstream's `Parser._open`
    // does this for every format, and `CTCF_PE_*_chr22_50k.bedpe.gz` is gzipped.
    let mut r = macs_io::open_maybe_gzip(path)?;
    let mut line = Vec::new();
    let mut b = macs_track::FragTrackBuilder::new();
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let Some(rec) = macs_io::parse_bedpe_line(&line)? else {
            continue;
        };
        if rec.left < 0 || rec.right <= rec.left || rec.chrom.is_empty() {
            continue;
        }
        if let Some(n) = rec.count {
            b.push_with_count(&rec.chrom, rec.left as u64, rec.right as u64, n);
        } else {
            b.push(&rec.chrom, rec.left as u64, rec.right as u64);
        }
    }
    b.finalize();
    Ok(b.build())
}

pub fn predictd(o: &Options) -> Result<()> {
    let ifile = o
        .get("ifile")
        .ok_or_else(|| MacsError::InvalidParameter("-i/--ifile is required".into()))?;
    let format = o.get("format").unwrap_or("AUTO").to_uppercase();
    if format == "BAMPE" || format == "BEDPE" {
        // `predictd_cmd.py:54-60`: in paired-end mode predictd does not fit a model at
        // all. It prints the average insertion length and returns, so there is no output
        // file and the exit status is 0. Reporting this as an unsupported input put a
        // spurious error on a completely ordinary `--format BAMPE` run.
        // `-f BAMPE` is a BAM, not a BEDPE text file. Routing it through the BEDPE
        // parser read the BGZF header as text and failed with "alignment record too
        // large or malformed"; upstream's `BAMPEParser` streams the BAM and reports
        // the mean template length, so `predictd -f BAMPE` printed nothing at all here.
        let t = if format == "BAMPE" {
            let (frags, _summary) = macs_io::bam::bampe_fragments(Path::new(ifile))?;
            let mut b = macs_track::FragTrackBuilder::new();
            for fr in &frags {
                b.push(&fr.chrom, u64::from(fr.start), u64::from(fr.start + fr.len));
            }
            b.finalize();
            b.build()
        } else {
            load_frag_track(Path::new(ifile))?
        };
        // `%d` in upstream (`predictd_cmd.py:56`) truncates the float, so a mean of
        // 98.307 is reported as "98 bps".
        let d = t.average_template_length() as i64;
        eprintln!(
            "predictd: total fragments/pairs in alignment file: {}",
            t.total()
        );
        // Trailing space: upstream's log format is
        // `'%(levelname)-5s @ %(asctime)s: %(message)s '` (`MACS3/Utilities/Logger.py:28`),
        // so every line it emits ends in one, and `cmdlinetest`'s
        // `perl -pe 's/^.*\s(\d+)\sbps/$1/'` preserves it -- `253 `, not `253`.
        eprintln!("# Average insertion length of all pairs is {d} bps ");
        return Ok(());
    }
    let track = if format == "BAM" {
        super::input::load_single_end_bam(Path::new(ifile))?.0
    } else if format == "SAM" {
        super::input::load_single_end_sam(Path::new(ifile))?.0
    } else {
        load_bed(Path::new(ifile))?
    };

    let mfold: (f64, f64) = match o.get("mfold") {
        Some(s) => {
            let mut it = s.split(',');
            let lo = it
                .next()
                .and_then(|x| x.parse().ok())
                .ok_or_else(|| MacsError::InvalidParameter("--mfold must be `lo,hi`".into()))?;
            let hi = it.next().and_then(|x| x.parse().ok()).unwrap_or(lo);
            (lo, hi)
        }
        None => (5.0, 50.0),
    };
    // `-g/--gsize` is a string that is either a shortcut or a number, and upstream
    // rejects anything else with exit 1 before doing any work
    // (`OptValidator.py:274-281`). Reading it with `Options::float` made a bad value
    // fall through to the default genome size, so `predictd -g nonsense` ran to
    // completion here while upstream refused.
    let gsize = match o.get("gsize") {
        Some(spec) => {
            macs_core::genomesize::resolve_gsize(spec).map_err(MacsError::InvalidParameter)?
        }
        None => 2.65e9,
    };
    let options = ModelOptions {
        gsize,
        mfold,
        bw: o.int("bw").unwrap_or(300),
        d_min: o.int("d_min").unwrap_or(20),
    };

    // F200: upstream catches `NotEnoughPairsException` (`predictd_cmd.py:82-83`), warns,
    // and returns normally -- exit status 0, no output file. An error here would flip
    // the exit status on inputs upstream treats as a successful, empty run.
    let model = match PeakModel::build(&track, options) {
        Ok(m) => m,
        Err(MacsError::NotEnoughPairs { found, needed }) => {
            eprintln!("predictd: MACS3 needs at least {needed} paired peaks at + and - strand to build the model, but can only find {found}! Please make your MFOLD range broader and try again.");
            eprintln!("predictd: Can't find enough pairs of symmetric peaks to build model!");
            return Ok(());
        }
        Err(e) => return Err(e),
    };
    let outdir = PathBuf::from(o.get("outdir").unwrap_or("."));
    std::fs::create_dir_all(&outdir)?;
    let rfile = o.get("rfile").unwrap_or("predictd_model.R");
    let name = rfile; // upstream: `pdf("<rfile>_model.pdf")`

    // Upstream's wording and formatting are load-bearing here, not cosmetic: its own
    // `test/cmdlinetest` derives `standard_results_predictd/run_predictd.txt` by
    // grepping these two lines out of the log --
    //
    //     grep "predicted fragment length"   ... | perl -pe 's/^.*\s(\d+)\sbps/$1/'
    //     grep "alternative fragment length" ... | perl -pe 's/^.*\s(\S+)\sbps/$1/'
    //
    // so the phrases and the `%d` are what the expected value is read from. `%d`
    // applied to a Python float truncates toward zero, hence `d.trunc()` and not
    // rounding: a fitted d of 229.69 is reported as 229.
    //
    // The trailing space is not ours: upstream's log format is
    // `'%(levelname)-5s @ %(asctime)s: %(message)s '` (`MACS3/Utilities/Logger.py:28`),
    // so every line it emits ends in one, and the perl above preserves it as
    // `229 ` rather than `229`. Matching it makes the extracted file byte-identical
    // to `test/standard_results_predictd/run_predictd.txt`.
    eprintln!(
        "# predicted fragment length is {} bps ",
        model.d.trunc() as i64
    );
    let alt = model
        .alternative_d
        .iter()
        .map(|d| d.to_string())
        .collect::<Vec<_>>()
        .join(",");
    eprintln!("# alternative fragment length(s) may be {alt} bps ");
    eprintln!(
        "# Generate R script for model : {} ",
        outdir.join(rfile).display()
    );
    model2r_script(&model, &outdir.join(rfile), name)?;
    Ok(())
}
