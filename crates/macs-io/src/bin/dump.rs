//! Dumps parsed records in a stable text form so the differential harness can
//! compare Rust against the oracle record-by-record.
//!
//! ```text
//! macs-io-dump <format> <file>
//! ```
//!
//! Output is one record per line, tab-separated, chromosome first, so a diff
//! against upstream's own dump localises a disagreement to a single record
//! instead of failing the whole file. Order is input order on purpose: upstream
//! does not sort during parsing, so a sorted comparison would hide a reordering
//! bug that changes peak widths.

use macs_core::Result;
use macs_io::{parse_bedgraph_line, parse_bedpe_line, parse_frag_line};
use std::io::{BufRead, BufWriter, Write};
use std::path::PathBuf;
use std::process::ExitCode;

/// `macs-io-dump bam <bam> <chrom> <left> <right> <maxdup>`
///
/// Emits `BAMaccessor.get_reads_in_region` in the same column order as
/// `oracle/dump_bam_region.py`, so the two can be diffed directly. The packed
/// sequence and quality are hex-encoded for the same reason: a reader that
/// decodes correctly but packs incorrectly still shows up.
fn dump_bam(args: &[String]) -> Result<()> {
    use macs_io::bam::BamAccessor;
    let [bam, chrom, left, right, maxdup] = args else {
        eprintln!("usage: macs-io-dump bam <bam> <chrom> <left> <right> <maxdup>");
        std::process::exit(2);
    };
    let mut b = BamAccessor::open(std::path::Path::new(bam))?;
    let stdout = std::io::stdout();
    let mut w = std::io::BufWriter::new(stdout.lock());
    writeln!(
        w,
        "### references\t{}",
        b.chromosomes()
            .iter()
            .map(|c| String::from_utf8_lossy(c).into_owned())
            .collect::<Vec<_>>()
            .join(",")
    )?;
    writeln!(
        w,
        "### rlengths\t{}",
        b.rlengths()
            .iter()
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )?;
    let num = |s: &str| -> Result<u32> {
        s.parse()
            .map_err(|e| macs_core::MacsError::InvalidParameter(format!("{s}: {e}")))
    };
    let reads = b.reads_in_region(chrom.as_bytes(), num(left)?, num(right)?, num(maxdup)?)?;
    writeln!(w, "### n_reads\t{}", reads.len())?;
    for r in &reads {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            String::from_utf8_lossy(&r.name),
            String::from_utf8_lossy(&r.chrom),
            r.lpos,
            r.rpos,
            r.strand,
            hex(&r.seq),
            hex(&r.qual),
            r.cigar
                .iter()
                .map(u32::to_string)
                .collect::<Vec<_>>()
                .join(","),
            r.md,
            String::from_utf8_lossy(&r.sequence()),
            r.n_edits(),
            r.length(),
        )?;
    }
    w.flush()?;
    Ok(())
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}

fn main() -> ExitCode {
    if std::env::args().nth(1).as_deref() == Some("bam") {
        return match dump_bam(&std::env::args().skip(2).collect::<Vec<_>>()) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("macs-io-dump bam: {e}");
                ExitCode::FAILURE
            }
        };
    }
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("macs-io-dump: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let format = args
        .next()
        .map(|s| s.to_string_lossy().into_owned())
        .ok_or_else(|| {
            macs_core::MacsError::InvalidParameter(
                "usage: macs-io-dump <bed|bedpe|frag|bedgraph> <file>".into(),
            )
        })?;
    let path =
        PathBuf::from(args.next().ok_or_else(|| {
            macs_core::MacsError::InvalidParameter("missing file argument".into())
        })?);

    let mut r = macs_io::open_maybe_gzip(&path)?;
    macs_io::skip_leading_headers(&mut r);

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    match format.as_str() {
        "bed" => {
            let reader = macs_io::SingleEndReader::new(r);
            for rec in reader {
                let rec = rec?;
                writeln!(
                    out,
                    "{}\t{}\t{}",
                    String::from_utf8_lossy(&rec.chrom),
                    rec.pos,
                    strand_char(rec.strand)
                )?;
            }
        }
        // BEDPE emits 3 columns and FRAG emits 4 (the count). The FRAG *barcode*
        // is deliberately omitted: upstream does not retain it on the track, so
        // the oracle side cannot report it and emitting it here would compare
        // something the differential cannot check. The barcode is still parsed
        // and unit-tested; the count is what the track actually keeps.
        "bedpe" | "frag" => {
            let mut line = Vec::new();
            loop {
                line.clear();
                let n = r.read_until(b'\n', &mut line)?;
                if n == 0 {
                    break;
                }
                if format == "frag" {
                    if let Some(f) = parse_frag_line(&line)? {
                        writeln!(
                            out,
                            "{}\t{}\t{}\t{}",
                            String::from_utf8_lossy(&f.chrom),
                            f.left,
                            f.right,
                            f.count.unwrap_or(1)
                        )?;
                    }
                } else if let Some(f) = parse_bedpe_line(&line)? {
                    writeln!(
                        out,
                        "{}\t{}\t{}",
                        String::from_utf8_lossy(&f.chrom),
                        f.left,
                        f.right
                    )?;
                }
            }
        }
        "bedgraph" => {
            let mut line = Vec::new();
            loop {
                line.clear();
                let n = r.read_until(b'\n', &mut line)?;
                if n == 0 {
                    break;
                }
                if let Some(g) = parse_bedgraph_line(&line)? {
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{:?}",
                        String::from_utf8_lossy(&g.chrom),
                        g.start,
                        g.end,
                        g.value
                    )?;
                }
            }
        }
        other => {
            return Err(macs_core::MacsError::InvalidParameter(format!(
                "unknown format {other}"
            )))
        }
    }
    out.flush()?;
    Ok(())
}

fn strand_char(s: macs_core::Strand) -> char {
    match s {
        macs_core::Strand::Plus => '+',
        macs_core::Strand::Minus => '-',
        macs_core::Strand::Unknown => '.',
    }
}
