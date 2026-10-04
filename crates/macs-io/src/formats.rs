//! Per-line parsers, transcribed from `MACS3.IO.Parser`.
//!
//! These are pure functions over one raw line so the differential harness can
//! compare them record-by-record without buffering.

use crate::atoi::{atof, atoi};
use crate::{RawBedGraph, RawFragment, RawSingleEnd};
use macs_core::{MacsError, Result, Strand};

/// Which format a line is being read as.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordKind {
    /// BED / ELAND / bowtie, single-end. One 5' position per record.
    SingleEnd,
    /// BEDPE / BAMPE, paired-end. One fragment per record.
    Fragment,
    /// bedGraph, for coverage tracks.
    BedGraph,
    /// FRAG, the single-cell fragment format.
    Frag,
}

/// How a line's trailing whitespace is treated before splitting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineSplit {
    /// `rstrip()` then `split(b"\t")`: BED, BEDPE, FRAG.
    TabAfterRstrip,
    /// `split()` with no argument: bedGraph. Runs of whitespace collapse, so
    /// this differs from the tab case in two ways — any whitespace separates, and
    /// no empty field is produced for adjacent separators.
    AnyWhitespace,
}

/// A bounded, printable rendering of a raw input line, for error messages.
///
/// Parse errors used to carry `String::from_utf8_lossy(line)`, which is a
/// problem whenever the "line" is not text: point a BED parser at a BGZF/BAM
/// file, or at a gzipped BED that was not decompressed, and the first "line" is
/// a few kilobytes of binary. The resulting diagnostic dumped the gzip/BAM magic
/// and raw record bytes to stderr -- reproduced with
/// `macs3-rs pileup -i reads.bam -f BAM`, which emitted ~2.3 KB of binary.
///
/// So an error message reports *where* the failure was and shows only a short
/// ASCII-safe prefix. Non-printable bytes are escaped rather than passed through,
/// so the message stays on one line and cannot smuggle terminal control
/// sequences into a log.
pub(crate) fn line_preview(line: &[u8]) -> String {
    const MAX: usize = 60;
    let mut s = String::with_capacity(MAX + 8);
    for &b in line.iter().take(MAX) {
        match b {
            b'\t' => s.push_str("\\t"),
            b'\n' => s.push_str("\\n"),
            b'\r' => s.push_str("\\r"),
            0x20..=0x7e => s.push(b as char),
            // Escape C0 controls (including NUL) and DEL so a binary record
            // cannot inject control characters into the message.
            _ => s.push_str(&format!("\\x{b:02x}")),
        }
    }
    if line.len() > MAX {
        s.push_str("...");
    }
    if s.is_empty() {
        s.push_str("<empty line>");
    }
    s
}

/// Strips the trailing bytes of `line` that C's `str.rstrip()` would.///
/// Python's `bytes.rstrip()` removes whitespace and, notably, also treats some
/// non-ASCII bytes as whitespace. In practice BED files are ASCII, so ASCII
/// whitespace is the whole of it, but `\x0b`/`\x0c` must be included or a file
/// carrying them changes meaning.
fn rstrip_bytes(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && matches!(line[end - 1], b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c) {
        end -= 1;
    }
    &line[..end]
}

/// ASCII whitespace, in Python's `str.split()` sense.
#[inline]
fn is_ws(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

/// Splits a line the way the given format does.
pub fn split_line(line: &[u8], how: LineSplit) -> Vec<&[u8]> {
    match how {
        LineSplit::TabAfterRstrip => rstrip_bytes(line).split(|&b| b == b'\t').collect(),
        // Rust's `split` with a predicate emits an empty field between adjacent
        // separators; Python's argument-less `split()` does not. Collapse runs
        // by hand so bedGraph with aligned columns parses like upstream's.
        LineSplit::AnyWhitespace => {
            let mut out = Vec::new();
            let mut i = 0usize;
            while i < line.len() {
                while i < line.len() && is_ws(line[i]) {
                    i += 1;
                }
                if i >= line.len() {
                    break;
                }
                let start = i;
                while i < line.len() && !is_ws(line[i]) {
                    i += 1;
                }
                out.push(&line[start..i]);
            }
            out
        }
    }
}

/// Parses one BED line into a single 5' end.
///
/// Transcribed from `BEDParser.fw_parse_line` ([`Parser.py:441`]). The strand
/// column drives the coordinate: `+` reads use `fields[1]` and `-` reads use
/// `fields[2]`. `Ok(None)` means "no record" and is never produced by this
/// parser; every line yields either a record or an error, as upstream does.
///
/// # Errors
///
/// A strand byte that is neither `+` nor `-` raises a parse error, matching
/// upstream's `StrandFormatError`. `b'.'` and `b'?'` are *not* accepted here even
/// though `Strand::from_byte` tolerates them, because upstream rejects them.
pub fn parse_bed_line(line: &[u8]) -> Result<Option<RawSingleEnd>> {
    let f = split_line(line, LineSplit::TabAfterRstrip);
    let chrom = f[0].to_vec();
    // Indexing `f[0]` and `f[1]` deliberately propagates a panic-free error for a
    // one-field line; upstream's `except IndexError` handler itself indexes
    // fields[1], so such a line raises out of the parser.
    let start = *f
        .get(1)
        .ok_or_else(|| MacsError::BadAlignment(line_preview(line)))?;
    let bad = |b: &[u8]| {
        MacsError::BadAlignment(format!(
            "invalid strand byte {:?} in BED line: {}",
            b.first().map(|c| *c as char).unwrap_or('?'),
            line_preview(line)
        ))
    };
    match f.get(5).copied() {
        // No strand column: upstream's `IndexError` path, plus strand at the start.
        None | Some(b"+") => {
            let end = *f.get(2).unwrap_or(&start);
            Ok(Some(RawSingleEnd {
                chrom,
                pos: atoi(start),
                strand: Strand::Plus,
                end: atoi(end),
            }))
        }
        Some(b"-") => {
            let end = *f
                .get(2)
                .ok_or_else(|| MacsError::BadAlignment(line_preview(line)))?;
            Ok(Some(RawSingleEnd {
                chrom,
                pos: atoi(end),
                strand: Strand::Minus,
                end: atoi(start),
            }))
        }
        Some(other) => Err(bad(other)),
    }
}

/// Parses one BEDPE line into a fragment.
///
/// `BEDPEParser.pe_parse_line` ([`Parser.py:507`]) is `fields[0..3]` and nothing
/// else: the middle column is ignored, and there is no strand handling.
pub fn parse_bedpe_line(line: &[u8]) -> Result<Option<RawFragment>> {
    let f = split_line(line, LineSplit::TabAfterRstrip);
    let chrom = *f
        .first()
        .ok_or_else(|| MacsError::BadAlignment(line_preview(line)))?;
    let left = *f
        .get(1)
        .ok_or_else(|| MacsError::BadAlignment(line_preview(line)))?;
    let right = *f
        .get(2)
        .ok_or_else(|| MacsError::BadAlignment(line_preview(line)))?;
    Ok(Some(RawFragment {
        chrom: chrom.to_vec(),
        left: atoi(left),
        right: atoi(right),
        barcode: None,
        count: None,
    }))
}

/// Parses one FRAG line into a fragment with its barcode.
///
/// Transcribed from `FragParser.pe_parse_line` ([`Parser.py:585`]). The column
/// order is `chrom start end barcode count` — **the barcode is column 4 and the
/// count is column 5**, which is the reverse of the order in the MACS3
/// documentation's prose and is an easy thing to get wrong. The count is stored
/// in an `unsigned short`, so anything above 65535 saturates rather than
/// wrapping.
pub fn parse_frag_line(line: &[u8]) -> Result<Option<RawFragment>> {
    let f = split_line(line, LineSplit::TabAfterRstrip);
    let get = |i: usize| -> Result<&[u8]> {
        f.get(i)
            .copied()
            .ok_or_else(|| MacsError::BadAlignment(line_preview(line)))
    };
    Ok(Some(RawFragment {
        chrom: get(0)?.to_vec(),
        left: atoi(get(1)?),
        right: atoi(get(2)?),
        barcode: Some(get(3)?.to_vec()),
        // The count lands in an `unsigned short` via a plain C cast, which
        // **truncates** rather than clamping. Verified against MACS3 3.0.5:
        //
        //   count 70000 -> stored 4464   (= 70000 & 0xFFFF)
        //   count     -5 -> stored 65531  (= -5 as u16)
        //
        // The `except OverflowError: thiscount = 65535` branch in
        // `FragParser.pe_parse_line` therefore never fires -- Cython emits a
        // truncating cast, not a range check -- so upstream's "will be capped at
        // 65535" warning is dead code and clamping here would be wrong.
        count: Some((atoi(get(4)?) as i32 as u16) as u32),
    }))
}

/// Parses one bedGraph line.
///
/// From `BedGraphIO.read_bedGraph` ([`BedGraphIO.py:84`](https://github.com/macs3-project/MACS/blob/c544319/MACS3/IO/BedGraphIO.py)):
/// whitespace-split, four fields, `atoi` for the coordinates and `atof` for the
/// value.
pub fn parse_bedgraph_line(line: &[u8]) -> Result<Option<RawBedGraph>> {
    let f = split_line(line, LineSplit::AnyWhitespace);
    let get = |i: usize| -> Result<&[u8]> {
        f.get(i)
            .copied()
            .ok_or_else(|| MacsError::BadAlignment(line_preview(line)))
    };
    // A blank or all-whitespace line yields no fields at all once runs are
    // collapsed, and is not a record.
    if f.is_empty() {
        return Ok(None);
    }
    Ok(Some(RawBedGraph {
        chrom: get(0)?.to_vec(),
        start: atoi(get(1)?),
        end: atoi(get(2)?),
        value: atof(get(3)?),
    }))
}

/// One parsed SAM alignment: 5' end plus query length.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawSam {
    /// Reference name, with a trailing `.fa...` suffix stripped.
    pub chrom: Vec<u8>,
    /// 5' end, 0-based: `POS-1` for `+`, `POS-1` plus reference length for `-`.
    pub pos: i64,
    /// Strand from flag bit 4.
    pub strand: Strand,
    /// `SEQ` length, for `SAMParser.tsize` (first-10 mean).
    pub seqlen: i64,
}

/// Parse one SAM alignment line (`SAMParser.fw_parse_line`, `Parser.py:908`).
///
/// Returns `Ok(None)` for `@` headers, blank lines, and records the filter rejects.
/// The filter is the BAMPE-style one (`Parser.py:929-931`): reject unmapped (4),
/// secondary (256), QC-fail (512), supplementary (2048); for paired reads reject
/// not-proper-pair, mate-unmapped (8) and second mate (128). There is deliberately
/// **no** MAPQ gate -- upstream never checks it here.
///
/// Reference-consuming CIGAR ops for the `-` end are M/D/N/=/X (`Parser.py:946`),
/// matching the BAM decoder. `SEQ` of `*` yields length 1, which upstream counts
/// (it takes `len(fields[9])` unconditionally).
pub fn parse_sam_line(line: &[u8]) -> Result<Option<RawSam>> {
    use crate::atoi::atoi;
    if line.first() == Some(&b'@') {
        return Ok(None);
    }
    // rstrip, then tab-split. Upstream splits the rstripped line.
    let mut end = line.len();
    while end > 0 && matches!(line[end - 1], b'\n' | b'\r' | b' ' | b'\t') {
        end -= 1;
    }
    if end == 0 {
        return Ok(None);
    }
    let line = &line[..end];
    let mut fields = line.split(|&b| b == b'\t');
    let mut get = || fields.next().unwrap_or(b"");
    let (_qname, flag_s, rname, pos_s, _mapq, cigar, _rnext, _pnext, _tlen, seq) = (
        get(),
        get(),
        get(),
        get(),
        get(),
        get(),
        get(),
        get(),
        get(),
        get(),
    );
    if rname.is_empty() || rname == b"*" {
        return Ok(None);
    }
    let flag: i64 = atoi(flag_s);
    if (flag & 2820) != 0 || (flag & 1 != 0 && ((flag & 136) != 0 || (flag & 2) == 0)) {
        return Ok(None);
    }
    let pos1 = atoi(pos_s);
    if pos1 <= 0 {
        return Ok(None);
    }
    let pos0 = pos1 - 1;
    let (pos, strand) = if flag & 16 != 0 {
        // Minus strand: shift by reference-consuming CIGAR length.
        let mut ref_len = 0i64;
        let mut num = 0i64;
        for &b in cigar {
            if b.is_ascii_digit() {
                num = num * 10 + (b - b'0') as i64;
            } else {
                match b {
                    b'M' | b'D' | b'N' | b'=' | b'X' => ref_len += num,
                    _ => {}
                }
                num = 0;
            }
        }
        (pos0 + ref_len, Strand::Minus)
    } else {
        (pos0, Strand::Plus)
    };
    // Strip a `.fa...` suffix (`thisref[:rindex(".fa")]`).
    let chrom = match rname.windows(3).position(|w| w == b".fa") {
        Some(i) => rname[..i].to_vec(),
        None => rname.to_vec(),
    };
    if chrom.is_empty() {
        return Ok(None);
    }
    Ok(Some(RawSam {
        chrom,
        pos,
        strand,
        seqlen: seq.len() as i64,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bed_strand_selects_the_coordinate() {
        let plus = parse_bed_line(b"chr1\t10\t20\tn1\t0\t+").unwrap().unwrap();
        assert_eq!((plus.pos, plus.strand), (10, Strand::Plus));
        let minus = parse_bed_line(b"chr1\t10\t20\tn1\t0\t-").unwrap().unwrap();
        assert_eq!((minus.pos, minus.strand), (20, Strand::Minus));
    }

    #[test]
    fn bed_without_a_strand_column_defaults_to_plus() {
        let r = parse_bed_line(b"chr1\t10\t20\tn1\t0").unwrap().unwrap();
        assert_eq!(r.pos, 10);
        assert_eq!(r.strand, Strand::Plus);
    }

    #[test]
    fn bed_rejects_an_unknown_strand_byte() {
        assert!(parse_bed_line(b"chr1\t10\t20\tn1\t0\t.").is_err());
        assert!(parse_bed_line(b"chr1\t10\t20\tn1\t0\tx").is_err());
    }

    #[test]
    fn bed_splits_on_tabs_only() {
        // a space-delimited line is one field, which is an error upstream
        assert!(parse_bed_line(b"chr1 10 20 n1 0 +").is_err());
    }

    #[test]
    fn bedpe_uses_the_first_three_columns_only() {
        let r = parse_bedpe_line(b"chr1\t100\t200\tchr1\t500\t600\tname")
            .unwrap()
            .unwrap();
        assert_eq!(
            (r.chrom.as_slice(), r.left, r.right),
            (&b"chr1"[..], 100, 200)
        );
    }

    #[test]
    fn frag_puts_the_barcode_in_column_four_and_the_count_in_five() {
        let r = parse_frag_line(b"chr1\t100\t200\tAAAC\t3")
            .unwrap()
            .unwrap();
        assert_eq!(r.barcode.as_deref(), Some(&b"AAAC"[..]));
        assert_eq!(r.count, Some(3));
        assert_eq!((r.left, r.right), (100, 200));
    }

    #[test]
    fn frag_count_truncates_to_sixteen_bits_rather_than_clamping() {
        // both values below are what MACS3 3.0.5 actually stores, measured with
        // FragParser.build_petrack rather than inferred from the source
        let over = parse_frag_line(b"chr1\t1\t2\tA\t70000").unwrap().unwrap();
        assert_eq!(
            over.count,
            Some(4464),
            "70000 & 0xFFFF, not a clamp to 65535"
        );

        let neg = parse_frag_line(b"chr1\t1\t2\tA\t-5").unwrap().unwrap();
        assert_eq!(neg.count, Some(65531), "two's complement, not a clamp to 0");

        let exact = parse_frag_line(b"chr1\t1\t2\tA\t65535").unwrap().unwrap();
        assert_eq!(exact.count, Some(65535));
    }

    #[test]
    fn bedgraph_splits_on_any_whitespace_and_skips_blanks() {
        let r = parse_bedgraph_line(b"chr1  100 200  1.5").unwrap().unwrap();
        assert_eq!((r.start, r.end, r.value), (100, 200, 1.5));
        assert!(parse_bedgraph_line(b"\n").unwrap().is_none());
        assert!(parse_bedgraph_line(b"   \t \n").unwrap().is_none());
    }
}
