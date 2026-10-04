//! Tag-size estimation, transcribed from upstream `Parser.tsize()`.
//!
//! `randsample_cmd.py:87-89` only falls back to the parser's own estimate when the
//! user did *not* pass `-s/--tsize`:
//!
//! ```python
//! if not options.tsize:           # override tsize if user specified --tsize
//!     ttsize = tp.tsize()
//!     options.tsize = ttsize
//! ```
//!
//! So `--tsize` is **optional**, and omitting it is the normal path. Rejecting the
//! invocation when it is absent is not a drop-in-compatible behaviour.
//!
//! Two details of the upstream implementation are load-bearing and are reproduced
//! exactly here:
//!
//! * The text path accumulates into a `cython.int`, so `s/n` is **C integer division**
//!   -- a floor, not a rounding. A mean of 149.1 read lengths reports 149.
//! * The BAM path accumulates into a `cython.double`, so its `s/n` is **float**
//!   division truncated on the final cast. The same input therefore truncates in two
//!   different directions depending on the container format.
//!
//! The scan also *reads the stream twice*: `tsize()` seeks back to 0 afterwards and
//! re-skips the header lines, so detection never consumes a record.

use std::io::BufRead;
use std::path::Path;

use crate::{open_maybe_gzip, MacsError, Result};

/// Formats for which `tlen_parse_line` is `end - start`.
///
/// `Parser.py:432-439` is the BED implementation and `BEDPEParser` reuses it, so a
/// single-ended read and a whole fragment are both measured the same way.
fn is_interval_len_format(format: &str) -> bool {
    matches!(
        format.to_ascii_uppercase().as_str(),
        "BED" | "BEDPE" | "BED3" | "BED4" | "BED5" | "BED6" | "ELAND" | "ELANDMULTI" | "ELANDPE"
    )
}

/// `Parser.tlen_parse_line` for the tab-separated formats: field 3 minus field 2.
///
/// Mirrors the C `atoi` helper upstream uses, including its zero-on-garbage
/// behaviour, so a malformed line contributes `0` and is skipped rather than
/// aborting the scan.
fn interval_len(line: &[u8]) -> i32 {
    let line = trim_rstrip(line);
    if line.is_empty() {
        return 0;
    }
    let mut fields = line.split(|&b| b == b'\t');
    let Some(_chrom) = fields.next() else {
        return 0;
    };
    let Some(start) = fields.next() else {
        return 0;
    };
    let Some(end) = fields.next() else {
        return 0;
    };
    (crate::atoi(end) - crate::atoi(start)) as i32
}

fn trim_rstrip(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && matches!(line[end - 1], b'\n' | b'\r' | b' ' | b'\t') {
        end -= 1;
    }
    &line[..end]
}

/// `Parser.tsize()` for the text formats: mean length of the first 10 usable lines,
/// by C integer division.
///
/// `Parser.py:279-293` bounds the scan two ways -- at most 10 *successful* lines and
/// at most 10000 lines read -- because a file can begin with arbitrarily many
/// non-alignment lines. Returns `None` for the upstream `-1` sentinel when no line
/// yielded a positive length.
pub fn detect_tsize_text(path: &Path) -> Result<Option<i32>> {
    let mut r = open_maybe_gzip(path)?;
    let mut s: i64 = 0;
    let mut n: i64 = 0;
    let mut m: i64 = 0;
    let mut line = Vec::new();
    while n < 10 && m < 10000 {
        m += 1;
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let tlen = interval_len(&line);
        if tlen > 0 {
            s += tlen as i64;
            n += 1;
        }
    }
    Ok(if n != 0 { Some((s / n) as i32) } else { None })
}

/// `Parser.tsize()` for BAM: the mean `l_seq` of the first 10 alignment records.
///
/// `Parser.py:1036-1070` walks the raw BGZF container by hand rather than going
/// through the aligner, reading `l_seq` out of each fixed-size record at byte 16.
/// It averages in `double` and truncates on the cast, and -- unlike the text path --
/// it applies no validity filter, so it reads exactly 10 records.
pub fn detect_tsize_bam(path: &Path) -> Result<Option<i32>> {
    let raw = crate::bam::read_all_maybe_gzip(path)?;
    let mut s: f64 = 0.0;
    let mut n: i64 = 0;
    // `header_len` lives at offset 4; the reference dictionary follows it.
    if raw.len() < 12 {
        return Ok(None);
    }
    let header_len = i32::from_le_bytes([raw[4], raw[5], raw[6], raw[7]]) as usize;
    let mut pos = 8 + header_len;
    if pos + 4 > raw.len() {
        return Ok(None);
    }
    let nc = i32::from_le_bytes([raw[pos], raw[pos + 1], raw[pos + 2], raw[pos + 3]]);
    pos += 4;
    for _ in 0..nc {
        if pos + 4 > raw.len() {
            return Ok(None);
        }
        let nlength = i32::from_le_bytes([raw[pos], raw[pos + 1], raw[pos + 2], raw[pos + 3]]);
        pos += 4 + nlength as usize + 4;
    }
    while n < 10 {
        if pos + 4 > raw.len() {
            break;
        }
        let entrylength =
            i32::from_le_bytes([raw[pos], raw[pos + 1], raw[pos + 2], raw[pos + 3]]) as usize;
        pos += 4;
        if pos + 20 > raw.len() {
            break;
        }
        s +=
            i32::from_le_bytes([raw[pos + 16], raw[pos + 17], raw[pos + 18], raw[pos + 19]]) as f64;
        n += 1;
        pos += entrylength;
    }
    Ok(if n != 0 {
        Some((s / n as f64) as i32)
    } else {
        None
    })
}

/// Dispatch on the resolved input format.
///
/// `AUTO` sniffs the content, as upstream's `guess_parser` does: a BAM magic means
/// BAM, otherwise the interval-length text path (BED/ELAND/bowtie all carry
/// `chrom start end` as the first three fields, which is all `detect_tsize_text`
/// reads). Upstream tries BAM first, then BED, then the others, so a BAM is never
/// misread as text and text is never misread as BAM.
pub fn detect_tsize(path: &Path, format: &str) -> Result<Option<i32>> {
    let f = format.to_ascii_uppercase();
    if f == "BAM" {
        return detect_tsize_bam(path);
    }
    if f == "AUTO" && is_bam_file(path) {
        return detect_tsize_bam(path);
    }
    if f == "AUTO" || is_interval_len_format(&f) {
        return detect_tsize_text(path);
    }
    Err(MacsError::InvalidParameter(format!(
        "cannot estimate tag size for input format {format}"
    )))
}

/// True when the file (after gzip sniffing) starts with the BAM magic.
fn is_bam_file(path: &Path) -> bool {
    use std::io::Read;
    let Ok(f) = std::fs::File::open(path) else {
        return false;
    };
    let mut r = std::io::BufReader::new(f);
    let mut magic = [0u8; 4];
    if r.read_exact(&mut magic).is_err() {
        return false;
    }
    if magic[..2] != [0x1f, 0x8b] {
        // Not gzipped: compare the raw magic.
        return magic == *b"BAM\x01";
    }
    // Gzipped: decode the first block and check its magic. BGZF's first block holds
    // the BAM header, which starts with `BAM\1`.
    drop(r);
    let Ok(f) = std::fs::File::open(path) else {
        return false;
    };
    let mut dec = flate2::bufread::MultiGzDecoder::new(std::io::BufReader::new(f));
    let mut inner = [0u8; 4];
    dec.read_exact(&mut inner).is_ok() && &inner == b"BAM\x01"
}
