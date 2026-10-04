//! `macs3-rs callvar`: call variants in given peak regions from alignment BAM files.
//!
//! The flag surface comes from `oracle/flag_matrix.tsv` via [`crate::parse_flags`], so
//! this layer never re-implements argparse -- which matters here, because `callvar`'s
//! help text and defaults are unusual (`-Q` defaults to 20, `--max-ar` to 0.95,
//! `--top2alleles-mratio` to 0.8) and are all recorded in the matrix.
//!
//! Port notes (`callvar_cmd.py`, `bin/macs3:635-700`):
//! * `opt_validate_callvar` (`OptValidator.py:663-678`) validates essentially
//!   nothing: its whole body is `if options.np <= 0: options.np = 1`. So there is no
//!   filesystem validation to reproduce, and adding any here would *break* the
//!   accept/reject criterion by rejecting invocations upstream accepts.
//! * `--outdir` exists for argparse parity but is never used: the VCF path comes
//!   from `-o`.
//! * The VCF header is written before any peak is processed, and
//!   `##Program_Args` is a *reconstruction* (the parsed options appended to what the
//!   user typed), not an echo of `argv`. See [`macs_callvar::program_args`].
//!
//! The variant-calling kernel is not implemented yet; see `macs_callvar`'s module
//! docs. It refuses **before** creating the VCF, which keeps "errors are raised
//! before any output file is created" true rather than violated by a header-only file.

use std::path::PathBuf;

use macs_callvar::{CallvarOptions, MACS_VERSION};
use macs_core::{MacsError, Result};

use crate::Options;

fn want(o: &Options, dest: &str) -> Result<String> {
    o.values
        .get(dest)
        .cloned()
        .ok_or_else(|| MacsError::Rejected(format!("callvar: missing required option {dest}")))
}

pub fn callvar(o: &Options) -> Result<()> {
    let np = o
        .values
        .get("np")
        .and_then(|v| v.parse::<i64>().ok())
        .unwrap_or(1);

    let opts = CallvarOptions {
        peak_bed: PathBuf::from(want(o, "peakbed")?),
        tfile: PathBuf::from(want(o, "tfile")?),
        cfile: o.values.get("cfile").map(PathBuf::from),
        ofile: PathBuf::from(want(o, "ofile")?),
        outdir: o.values.get("outdir").map(PathBuf::from),
        gq_cutoff_hetero: num(o, "GQCutoffHetero", 0.0),
        gq_cutoff_homo: num(o, "GQCutoffHomo", 0.0),
        q: int(o, "Q", 20),
        max_duplicate: int(o, "maxDuplicate", 1),
        fermi: o
            .values
            .get("fermi")
            .cloned()
            .unwrap_or_else(|| "auto".into()),
        fermi_min_overlap: int(o, "fermiMinOverlap", 30),
        top2_alleles_min_ratio: num(o, "top2allelesMinRatio", 0.8),
        alt_allele_min_count: int(o, "altalleleMinCount", 2),
        max_ar: num(o, "maxAR", 0.95),
        np,
        verbose: int(o, "verbose", 2),
    };

    // Upstream clamps `np <= 0` to 1 in the validator rather than rejecting.
    if opts.np <= 0 {
        eprintln!("callvar: -m must be positive; using 1");
    }

    eprintln!("callvar: macs3-rs callvar (targets MACS {MACS_VERSION})");
    macs_callvar::run(&opts, &o.raw_argv)
}

fn num(o: &Options, dest: &str, default: f64) -> f64 {
    o.values
        .get(dest)
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn int(o: &Options, dest: &str, default: i64) -> i64 {
    o.values
        .get(dest)
        .and_then(|v| v.parse::<i64>().ok())
        .unwrap_or(default)
}
