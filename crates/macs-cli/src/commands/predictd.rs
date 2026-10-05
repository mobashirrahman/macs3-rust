//! `macs3-rs predictd`: estimate the fragment length `d` and emit the
//! `*_model.r` plotting script.
//!
//! Delegates to `macs_model::PeakModel` (the `PeakModel.py` port) and writes the
//! R script via `model2r_script`. Upstream also emits the R script from
//! `callpeak`'s model path, so both share this.

use std::path::PathBuf;

use macs_core::{MacsError, Result};
use macs_model::{model2r_script, ModelOptions, PeakModel};

use crate::Options;

pub fn predictd(o: &Options) -> Result<()> {
    let ifiles = super::input::input_files(o)?;
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
        let t = super::input::load_fragment_files(&ifiles, &format)?;
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
    let track = super::input::load_tag_files(&ifiles, &format, true)?.0;

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
