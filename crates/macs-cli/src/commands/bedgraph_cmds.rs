//! `macs3-rs bdgopt` and `cmbreps`: the two bedGraph-only commands.
//!
//! Both are thin -- read bedGraph(s), transform/combine them with the
//! `macs-bedgraph` library (whose `overlie` is already verified bit-exact
//! against upstream, F124), write a bedGraph. The arithmetic lives in the
//! library; this layer only maps the derived flags onto it.

use std::path::Path;

use macs_bedgraph::{BedGraph, Op};
use macs_core::{MacsError, Result};

use crate::Options;

fn read_bed(path: &Path) -> Result<BedGraph> {
    BedGraph::read(path, 0.0)
}

/// `bdgopt`: modify a single bedGraph (`multiply`/`add`/`max`/`min`/`p2q`).
pub fn bdgopt(o: &Options) -> Result<()> {
    let ifile = o
        .get("ifile")
        .ok_or_else(|| MacsError::InvalidParameter("-i/--ifile is required".into()))?;
    let method = o.get("method").unwrap_or("multiply").to_lowercase();
    let track = read_bed(Path::new(ifile))?;

    let modified = match method.as_str() {
        "p2q" => track.p2q(),
        "multiply" | "add" | "max" | "min" => {
            let p = o
                .get("extraparam")
                .ok_or_else(|| {
                    MacsError::InvalidParameter(format!("--method {method} needs -p/--extra-param"))
                })?
                .parse::<f32>()
                .map_err(|_| MacsError::InvalidParameter("-p must be a number".into()))?;
            let m = method.clone();
            track.apply_func(move |x| match m.as_str() {
                "multiply" => x * p,
                "add" => x + p,
                "max" => {
                    if x > p {
                        x
                    } else {
                        p
                    }
                }
                "min" => {
                    if x < p {
                        x
                    } else {
                        p
                    }
                }
                _ => x,
            })
        }
        other => {
            return Err(MacsError::InvalidParameter(format!(
                "unknown --method {other}"
            )));
        }
    };

    write_out(
        o,
        &modified,
        &format!("{}_modified_scores", method.to_uppercase()),
    )
}

/// `cmbreps`: combine several bedGraphs with `--method` (`max`/`mean`/`fisher`).
pub fn cmbreps(o: &Options) -> Result<()> {
    let ifiles = o.get_all("ifile");
    if ifiles.len() < 2 {
        return Err(MacsError::InvalidParameter(
            "cmbreps needs at least two -i/--ifile".into(),
        ));
    }
    let method = o.get("method").unwrap_or("max").to_lowercase();
    let op = Op::parse(method.as_bytes()).ok_or_else(|| {
        MacsError::InvalidParameter(format!("unsupported cmbreps --method {method}"))
    })?;
    let mut tracks: Vec<BedGraph> = Vec::new();
    for f in ifiles {
        tracks.push(read_bed(Path::new(f))?);
    }
    let first = &tracks[0];
    let others: Vec<&BedGraph> = tracks[1..].iter().collect();
    // Streamed: the combined track is only ever written, so building it as a
    // `BedGraph` first costs a whole extra copy of the genome's runs.
    let ofile = o
        .get("ofile")
        .ok_or_else(|| MacsError::InvalidParameter("-o/--ofile is required".into()))?;
    let outdir = o.get("outdir").unwrap_or(".");
    std::fs::create_dir_all(outdir)?;
    let desc = format!("{}_combined_scores", method.to_uppercase());
    let name = desc.trim_end_matches("_combined_scores");
    first.overlie_write(
        &others,
        op,
        &Path::new(outdir).join(ofile),
        true,
        &desc,
        &format!("Scores calculated by {}", name.to_uppercase()),
    )
}

fn write_out(o: &Options, bg: &BedGraph, desc: &str) -> Result<()> {
    let ofile = o
        .get("ofile")
        .ok_or_else(|| MacsError::InvalidParameter("-o/--ofile is required".into()))?;
    let outdir = o.get("outdir").unwrap_or(".");
    std::fs::create_dir_all(outdir)?;
    // upstream's write_bedGraph writes a `track` line with the method name as
    // the track name and "<method> ..." as the description.
    let name = desc
        .trim_end_matches("_modified_scores")
        .trim_end_matches("_combined_scores");
    bg.write(
        &Path::new(outdir).join(ofile),
        true,
        desc,
        &format!("Scores calculated by {}", name.to_uppercase()),
    )
}
