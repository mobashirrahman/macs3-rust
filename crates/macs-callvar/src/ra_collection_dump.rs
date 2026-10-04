//! A small driver that reproduces `RACollection`'s consensus construction for a peak
//! file plus BAMs, in exactly the field order `oracle/dump_callvar_refseq.py` emits.
//!
//! It exists so a consensus divergence is reported per peak and per field, rather than
//! inferred from a VCF that is wrong everywhere downstream of the first shifted byte.

use crate::{CallParams, RACollection};
use macs_io::bam::BamAccessor;
use std::path::Path;

/// One peak's reconstruction.
#[derive(Debug, Clone)]
pub struct PeakRecord {
    pub chrom: Vec<u8>,
    pub left: i64,
    pub right: i64,
    pub read_start: i64,
    pub read_end: i64,
    pub count_t: usize,
    pub count_c: usize,
    pub ext_len: usize,
    /// `peak_refseq`, lower-case hex.
    pub seq_hex: String,
    /// `peak_refseq_ext`, lower-case hex.
    pub ext_hex: String,
}

/// Inputs to the dump.
#[derive(Debug, Clone)]
pub struct DumpOptions<'a> {
    pub peaks: &'a Path,
    pub tbam: &'a Path,
    pub cbam: Option<&'a Path>,
    pub max_duplicate: u32,
}

fn hex(b: &[u8]) -> String {
    let mut s = String::with_capacity(b.len() * 2);
    for x in b {
        s.push_str(&format!("{x:02x}"));
    }
    s
}

/// Read a peak file as `callvar_cmd.run` does: first three whitespace fields.
fn read_peaks(path: &Path) -> std::io::Result<Vec<(Vec<u8>, i64, i64)>> {
    let text = std::fs::read_to_string(path)?;
    let mut out = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 3 {
            continue;
        }
        let start: i64 = f[1].parse().map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("bad start in {line:?}"),
            )
        })?;
        let end: i64 = f[2].parse().map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("bad end in {line:?}"),
            )
        })?;
        out.push((f[0].as_bytes().to_vec(), start, end));
    }
    out.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    Ok(out)
}

/// Build a record for every read-covered peak. Peaks with no treatment reads are
/// skipped, matching upstream's "No reads found in this peak. Skipped".
pub fn dump_peak_records(o: &DumpOptions<'_>) -> Result<Vec<PeakRecord>, String> {
    let peaks = read_peaks(o.peaks).map_err(|e| e.to_string())?;
    let mut tbam = BamAccessor::open(o.tbam).map_err(|e| e.to_string())?;
    let mut cbam = match o.cbam {
        Some(p) => Some(BamAccessor::open(p).map_err(|e| e.to_string())?),
        None => None,
    };
    let mut out = Vec::new();
    for (chrom, start, end) in peaks {
        let t = tbam
            .reads_in_region(
                &chrom,
                start.max(0) as u32,
                end.max(0) as u32,
                o.max_duplicate,
            )
            .map_err(|e| e.to_string())?;
        if t.is_empty() {
            continue;
        }
        let c = match cbam.as_mut() {
            Some(cb) => cb
                .reads_in_region(
                    &chrom,
                    start.max(0) as u32,
                    end.max(0) as u32,
                    o.max_duplicate,
                )
                .map_err(|e| e.to_string())?,
            None => Vec::new(),
        };
        let mut coll = match RACollection::new(&chrom, start, end, t, c) {
            Ok(x) => x,
            Err(e) if e == crate::NO_READS => continue,
            Err(e) => return Err(e),
        };
        coll.remove_outliers_in_place(5);
        out.push(PeakRecord {
            chrom,
            left: coll.left,
            right: coll.right,
            read_start: coll.read_start,
            read_end: coll.read_end,
            count_t: coll.treatment_reads().len(),
            count_c: coll.control_reads().len(),
            ext_len: coll.peak_refseq_ext.len(),
            seq_hex: hex(&coll.peak_refseq),
            ext_hex: hex(&coll.peak_refseq_ext),
        });
    }
    let _ = CallParams::default();
    Ok(out)
}
