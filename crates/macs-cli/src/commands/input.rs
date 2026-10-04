//! Shared input loading for the subcommands.
//!
//! # Why this exists
//!
//! Five subcommands -- `filterdup`, `pileup`, `randsample`, `refinepeak` and
//! `predictd` -- each carried their own private `load_bed`, and all five were
//! byte-for-byte identical copies that opened the file with
//! `std::fs::File::open` and read it as UTF-8 text. That is what made gzipped
//! input fail:
//!
//! ```text
//! $ macs3-rs filterdup -i CTCF_SE_ChIP_chr22_50k.bed.gz -f BED -o o.bed
//! macs3-rs filterdup: alignment record too large or malformed: <binary>
//! ```
//!
//! Upstream runs that command, because `Parser._open` sniffs the gzip magic
//! itself and wraps the stream in `gzip.open` when it is present. `callpeak`
//! already went through [`macs_io::open_maybe_gzip`], which is why `callpeak
//! -t x.bed.gz` worked while these five did not -- an inconsistency that also
//! silently took out `bdgcmp`, `bdgpeakcall`, `bdgbroadcall` and `cmbreps` in
//! any pipeline where the bedGraph came from `pileup`.
//!
//! One copy, going through the same sniffing reader `callpeak` uses, is the fix.

use macs_core::Result;
use macs_track::SingleEndTrack;
use std::path::Path;

/// Read a single-end BED/ELAND/bowtie file into a track, transparently
/// decompressing gzip.
///
/// Detection is by content, not by file extension: `macs_io::open_maybe_gzip`
/// probes the `1f 8b` magic, so a gzipped stream named `.bed` is still
/// decompressed and a plain file named `.bed.gz` is still read as text. That
/// matches upstream, which also sniffs rather than trusting the suffix.
pub fn load_single_end_bed(path: &Path) -> Result<SingleEndTrack> {
    use std::io::BufRead;
    let mut r = macs_io::open_maybe_gzip(path)?;
    let mut line = Vec::new();
    let mut b = macs_track::SingleEndTrackBuilder::new();
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
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
    Ok(b.build())
}

/// Read a single-end BAM into a track, without requiring an index.
///
/// Upstream's `BAMParser` streams the file and needs no `.bai`; the indexed
/// `BamAccessor` would reject index-less inputs upstream accepts. Returns the track
/// and the mean query length (`BAMParser.tsize`, first-10 average) for the header.
pub fn load_single_end_bam(path: &Path) -> Result<(SingleEndTrack, f64)> {
    let (tags, mean_qlen) = macs_io::bam::se_bam_tags(path)?;
    let mut b = macs_track::SingleEndTrackBuilder::new();
    for t in &tags {
        let strand = if t.strand == 1 {
            macs_core::Strand::Minus
        } else {
            macs_core::Strand::Plus
        };
        b.push(&t.chrom, u64::from(t.pos), strand);
    }
    b.finalize();
    Ok((b.build(), mean_qlen))
}

/// Read a single-end SAM file into a track, transparently decompressing gzip.
///
/// `SAMParser` skips `@` headers and parses FLAG/RNAME/POS/CIGAR/SEQ
/// (`Parser.py:908`). The filter, 5' convention and `.fa` stripping match
/// `parse_sam_line`; returns the track and the first-10 mean SEQ length for the
/// `# tag size` header (`SAMParser.tlen_parse_line`).
pub fn load_single_end_sam(path: &Path) -> Result<(SingleEndTrack, f64)> {
    use std::io::BufRead;
    let mut r = macs_io::open_maybe_gzip(path)?;
    let mut line = Vec::new();
    let mut b = macs_track::SingleEndTrackBuilder::new();
    let (mut qsum, mut qn) = (0u64, 0u64);
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let Some(rec) = macs_io::parse_sam_line(&line)? else {
            continue;
        };
        if rec.pos < 0 || rec.chrom.is_empty() {
            continue;
        }
        if qn < 10 && rec.seqlen > 0 {
            qsum += rec.seqlen as u64;
            qn += 1;
        }
        b.push(&rec.chrom, rec.pos as u64, rec.strand);
    }
    b.finalize();
    let mean_qlen = if qn > 0 { qsum as f64 / qn as f64 } else { 0.0 };
    Ok((b.build(), mean_qlen))
}
