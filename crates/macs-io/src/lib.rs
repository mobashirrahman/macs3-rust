//! Input parsers for `macs3-rs`.
//!
//! Every parser here is a transcription of `MACS3.IO.Parser` from MACS3 3.0.5
//! (`c544319`). The transcription is literal, including the parts that look like
//! mistakes, because a "sensible" fix changes every peak:
//!
//! # C `atoi` / `atof` semantics
//!
//! Upstream calls libc's `atoi`/`atof` directly ([`Parser.py:26`]). They are
//! *lenient*: leading whitespace and a sign are consumed, digits are accumulated,
//! and **anything unparseable yields `0` rather than an error**. A BED line with
//! `start=abc` is therefore read as position 0. [`atoi`] and [`atof`] reproduce
//! that, and the differential harness relies on it.
//!
//! # BED single-end: the 5' coordinate depends on the strand column
//!
//! [`BEDParser.fw_parse_line`](Parser.py:441):
//!
//! ```text
//! strand == b"+"  ->  position = atoi(fields[1])   # the start
//! strand == b"-"  ->  position = atoi(fields[2])   # the END, not end-1
//! anything else    ->  StrandFormatError
//! fewer than 6 columns -> position = atoi(fields[1]), strand = plus
//! ```
//!
//! So a BED interval `[100, 200)` contributes a plus read at 100 and a minus read
//! at 200. The minus read is one base *past* the last covered base, and that is
//! what the pileup then extends leftwards from. Reproduced exactly.
//!
//! # Line splitting
//!
//! `rstrip()` then `split(b"\t")` for BED/BEDPE/FRAG — trailing whitespace of any
//! kind is stripped, and **tab is the only separator**, so a space-delimited file
//! parses as one field per line and produces a position of 0. bedGraph is the
//! exception: it uses `i.split()`, i.e. any whitespace run.
//!
//! # Header handling
//!
//! `skip_first_commentlines` runs once, when the parser is opened, and only
//! skips lines at the *top* of the file whose first five bytes are `track`, or
//! first seven are `browser`, or first byte is `#`. A `#` line in the middle of
//! the file is *not* skipped by BED; it is parsed as a record.

#![forbid(unsafe_code)]
#![deny(missing_debug_implementations)]

mod atoi;
pub mod bam;
mod formats;
pub mod peakout;
pub mod tsize;

pub use atoi::{atof, atoi};
pub use formats::{
    parse_bed_line, parse_bedgraph_line, parse_bedpe_line, parse_frag_line, parse_sam_line,
    split_line, LineSplit, RawSam, RecordKind,
};
pub use peakout::{format_g, narrowpeak_row, python_repr, summit_row, xls_body, xls_row, XlsRow};
pub use tsize::{detect_tsize, detect_tsize_bam, detect_tsize_text};

use macs_core::{ChromId, Coord, MacsError, Result, Strand};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek};
use std::path::Path;

/// The read buffer size upstream uses when streaming a file
/// (`MACS3.IO.Parser.READ_BUFFER_SIZE`). Only affects I/O granularity, not
/// results, but matching it keeps syscall profiles comparable in benchmarks.
pub const READ_BUFFER_SIZE: usize = 100_000;

/// A 5' end on a strand, before it is interned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FivePrime {
    /// 0-based position, already adjusted for the strand by the parser.
    pub pos: i64,
    /// The read's strand.
    pub strand: Strand,
}

/// A parsed single-end record, with the chromosome still as raw bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawSingleEnd {
    /// Chromosome name, borrowed from the input line.
    pub chrom: Vec<u8>,
    /// The 5' position, as `atoi` produced it.
    pub pos: i64,
    /// Strand: `+` unless the file said `-`.
    pub strand: Strand,
    /// The tag's **end**, as column 3.
    ///
    /// Not needed to build a track -- only the 5' position is piled up -- but
    /// upstream's `BEDParser.d` is the mean of `end - pos` over every tag read,
    /// and `callpeak_cmd.py:80` prints it in the xls header. A minus-strand line
    /// carries `end` as the 5' coordinate, so the length is `pos - end` there.
    pub end: i64,
}

/// A parsed fragment record (BEDPE, FRAG, BAMPE).
#[derive(Debug, Clone, PartialEq)]
pub struct RawFragment {
    /// Chromosome name.
    pub chrom: Vec<u8>,
    /// Left end.
    pub left: i64,
    /// Right end.
    pub right: i64,
    /// Barcode, for FRAG.
    pub barcode: Option<Vec<u8>>,
    /// Multiplicity, for FRAG. Capped at 65535 by upstream.
    pub count: Option<u32>,
}

/// A bedGraph record.
#[derive(Debug, Clone, PartialEq)]
pub struct RawBedGraph {
    /// Chromosome name.
    pub chrom: Vec<u8>,
    /// Start.
    pub start: i64,
    /// End.
    pub end: i64,
    /// Value.
    pub value: f64,
}

/// Opens a possibly gzipped file for buffered line reading.
///
/// Upstream sniffs the gzip magic itself and wraps `gzip.open` when present, so
/// a `.gz` file does not need a `.gz` name. Reproduced: the extension is
/// irrelevant.
pub fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|e| {
        MacsError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })?;
    let mut magic = [0u8; 2];
    let mut probe = file;
    let gz = match probe.read_exact(&mut magic) {
        Ok(()) => magic == [0x1f, 0x8b],
        // fewer than 2 bytes: not gzip
        Err(_) => false,
    };
    if gz {
        // rewind: the magic probe consumed two bytes, and the decoder must see
        // the whole stream or it reports a truncated-member error
        probe
            .seek(std::io::SeekFrom::Start(0))
            .map_err(MacsError::Io)?;
        Ok(Box::new(BufReader::with_capacity(
            READ_BUFFER_SIZE,
            flate2::read::MultiGzDecoder::new(probe),
        )))
    } else {
        probe
            .seek(std::io::SeekFrom::Start(0))
            .map_err(MacsError::Io)?;
        Ok(Box::new(BufReader::with_capacity(READ_BUFFER_SIZE, probe)))
    }
}

/// `True` when the stream starts with the gzip magic.
pub fn is_gzipped(path: &Path) -> bool {
    let mut f = match File::open(path) {
        Ok(f) => f,
        Err(_) => return false,
    };
    let mut magic = [0u8; 2];
    matches!(f.read_exact(&mut magic), Ok(())) && magic == [0x1f, 0x8b]
}

/// Skips the leading `track` / `browser` / `#` lines of a stream.
///
/// Mirrors `BEDParser.skip_first_commentlines`
/// ([`Parser.py:415`](https://github.com/macs3-project/MACS/blob/c544319/MACS3/IO/Parser.py)):
/// the test is on the *first bytes* of the line, not on a whole-line match, and
/// it only applies at the top of the file.
pub fn skip_leading_headers<R: BufRead>(r: &mut R) {
    loop {
        let buf = match r.fill_buf() {
            Ok(b) => b,
            Err(_) => return,
        };
        if buf.is_empty() {
            return;
        }
        let line: &[u8] = match buf.iter().position(|&b| b == b'\n') {
            Some(i) => &buf[..=i],
            // A first line longer than the buffer is data, not a header.
            None => return,
        };
        let skip =
            line.starts_with(b"track") || line.starts_with(b"browser") || line.starts_with(b"#");
        if !skip {
            return;
        }
        let n = line.len();
        r.consume(n);
    }
}

/// A record source that yields already-parsed single-end records.
///
/// This is the boundary the differential harness tests: everything downstream
/// sees only these values, so if they match upstream's track contents, the rest
/// of the pipeline is being compared on equal terms.
#[derive(Debug)]
pub struct SingleEndReader<R: BufRead> {
    inner: R,
    line_no: u64,
    /// Records rejected because the position was negative or the name empty.
    pub skipped: u64,
    last: Vec<u8>,
}

impl<R: BufRead> SingleEndReader<R> {
    /// Wrap a stream. Call [`skip_leading_headers`] first, as upstream does at
    /// open time.
    pub fn new(inner: R) -> Self {
        SingleEndReader {
            inner,
            line_no: 0,
            skipped: 0,
            last: Vec::new(),
        }
    }

    /// The most recent raw line, for error messages.
    pub fn last_line(&self) -> &[u8] {
        &self.last
    }

    /// The 1-based line number of the most recent record.
    pub fn line_no(&self) -> u64 {
        self.line_no
    }
}

impl<R: BufRead> Iterator for SingleEndReader<R> {
    type Item = Result<RawSingleEnd>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let mut line = std::mem::take(&mut self.last);
            line.clear();
            let n = match self.inner.read_until(b'\n', &mut line) {
                Ok(n) => n,
                Err(e) => return Some(Err(MacsError::Io(e))),
            };
            if n == 0 {
                return None;
            }
            self.line_no += 1;
            match parse_bed_line(&line) {
                Err(e) => return Some(Err(e)),
                Ok(None) => continue,
                Ok(Some(rec)) => {
                    if rec.pos < 0 || rec.chrom.is_empty() {
                        // upstream: `if fpos < 0 or not chromosome: continue`
                        self.skipped += 1;
                        continue;
                    }
                    return Some(Ok(rec));
                }
            }
        }
    }
}

/// Interns parsed records onto a genome, producing dense per-chromosome arrays.
///
/// This is the `Chromosome { plus: Vec<Coord>, minus: Vec<Coord> }` layout the
/// plan calls for: no `String` per read, no per-read allocation after this point.
#[derive(Debug, Default)]
pub struct SingleEndStore {
    /// Interned chromosome ids, in first-appearance order.
    genome: macs_core::Genome,
    plus: Vec<Vec<Coord>>,
    minus: Vec<Vec<Coord>>,
    total: u64,
}

impl SingleEndStore {
    /// An empty store.
    pub fn new() -> Self {
        SingleEndStore {
            genome: macs_core::Genome::new(),
            plus: Vec::new(),
            minus: Vec::new(),
            total: 0,
        }
    }

    /// The genome dictionary.
    pub fn genome(&self) -> &macs_core::Genome {
        &self.genome
    }

    /// Total reads stored.
    pub fn total(&self) -> u64 {
        self.total
    }

    /// Add one record, interning its chromosome if new.
    pub fn push(&mut self, rec: &RawSingleEnd) {
        let id = self.intern(&rec.chrom);
        let pos = rec.pos.max(0) as Coord;
        match rec.strand {
            Strand::Plus => self.plus[id.0 as usize].push(pos),
            Strand::Minus => self.minus[id.0 as usize].push(pos),
            Strand::Unknown => self.plus[id.0 as usize].push(pos),
        }
        self.total += 1;
    }

    fn intern(&mut self, name: &[u8]) -> ChromId {
        if let Some(id) = self.genome.get(name) {
            return id;
        }
        let id = self.genome.intern(name);
        self.plus.push(Vec::new());
        self.minus.push(Vec::new());
        id
    }

    /// Sort each per-chromosome array and return `(chrom, plus, minus)`.
    ///
    /// Upstream sorts its `locations` arrays in `FWTrack.finalize`; peak calling
    /// depends on it because the pileup sweep assumes sorted endpoints.
    pub fn finalize(&mut self) {
        for v in &mut self.plus {
            v.sort_unstable();
        }
        for v in &mut self.minus {
            v.sort_unstable();
        }
    }

    /// The 5' positions on `chrom` for one strand.
    pub fn positions(&self, chrom: ChromId, strand: Strand) -> &[Coord] {
        match strand {
            Strand::Minus => &self.minus[chrom.0 as usize],
            _ => &self.plus[chrom.0 as usize],
        }
    }

    /// Reads per chromosome, for diagnostics and the differential report.
    pub fn counts(&self, chrom: ChromId) -> (usize, usize) {
        (
            self.plus[chrom.0 as usize].len(),
            self.minus[chrom.0 as usize].len(),
        )
    }

    /// Every chromosome id, in interning order.
    pub fn chroms(&self) -> Vec<ChromId> {
        self.genome.ids_file_order()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    // `write_all` on the gzip encoder; the non-test build does not need `Write`
    use std::io::Write as _;

    fn parse_all(text: &[u8]) -> Vec<RawSingleEnd> {
        let mut r = Cursor::new(text.to_vec());
        skip_leading_headers(&mut r);
        SingleEndReader::new(r)
            .map(|x| x.expect("parse error"))
            .collect()
    }

    #[test]
    fn plus_reads_use_the_start_and_minus_reads_use_the_end() {
        // the single most load-bearing rule in the BED parser
        let text = b"chr1\t100\t200\tr1\t0\t+\nchr1\t100\t200\tr2\t0\t-\n";
        let recs = parse_all(text);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].pos, 100, "a + read sits at the BED start");
        assert_eq!(recs[1].pos, 200, "a - read sits at the BED end, not end-1");
        assert_eq!(recs[0].strand, Strand::Plus);
        assert_eq!(recs[1].strand, Strand::Minus);
    }

    #[test]
    fn a_short_line_defaults_to_plus_at_the_start() {
        // BED3: only 3 columns, so the IndexError path fires
        let recs = parse_all(b"chr1\t100\t200\n");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].pos, 100);
        assert_eq!(recs[0].strand, Strand::Plus);
    }

    #[test]
    fn an_unknown_strand_column_is_rejected() {
        let mut r = Cursor::new(b"chr1\t100\t200\tr\t0\t.\n".to_vec());
        skip_leading_headers(&mut r);
        let first = SingleEndReader::new(r).next().expect("one item");
        assert!(first.is_err(), "a '.' strand must be an error, not a guess");
    }

    #[test]
    fn leading_headers_are_skipped_but_later_ones_are_not() {
        let text = b"# a comment\ntrack name=x\nbrowser hide all\nchr1\t10\t20\tr\t0\t+\n";
        let recs = parse_all(text);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].pos, 10);
    }

    #[test]
    fn trailing_whitespace_is_stripped_and_crlf_is_tolerated() {
        let recs = parse_all(b"chr1\t10\t20\tr\t0\t+\r\n");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].strand, Strand::Plus);
    }

    #[test]
    fn a_missing_trailing_newline_still_yields_a_record() {
        let recs = parse_all(b"chr1\t10\t20\tr\t0\t+");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].pos, 10);
    }

    #[test]
    fn unparseable_coordinates_become_zero_like_c_atoi() {
        let recs = parse_all(b"chr1\tabc\txyz\tr\t0\t+\n");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].pos, 0, "atoi returns 0 for garbage, not an error");
    }

    #[test]
    fn a_line_with_one_field_is_an_error_not_a_silent_zero() {
        // upstream's `except IndexError` handler itself indexes fields[1], so a
        // one-field line raises IndexError out of the parser
        let mut r = Cursor::new(b"chr1\n".to_vec());
        skip_leading_headers(&mut r);
        let first = SingleEndReader::new(r).next().expect("one item");
        assert!(first.is_err());
    }

    #[test]
    fn the_store_interns_chromosomes_and_sorts() {
        let mut s = SingleEndStore::new();
        for (c, p, st) in [
            ("chr2", 300i64, Strand::Plus),
            ("chr1", 100, Strand::Plus),
            ("chr1", 50, Strand::Minus),
            ("chr1", 900, Strand::Plus),
        ] {
            s.push(&RawSingleEnd {
                chrom: c.as_bytes().to_vec(),
                pos: p,
                strand: st,
                end: p + 50,
            });
        }
        s.finalize();
        assert_eq!(s.total(), 4);
        assert_eq!(s.genome().name_string(macs_core::ChromId(0)), "chr2");
        assert_eq!(s.genome().name_string(macs_core::ChromId(1)), "chr1");
        let c1 = s.genome().get(b"chr1").unwrap();
        assert_eq!(s.positions(c1, Strand::Plus), &[100, 900]);
        assert_eq!(s.positions(c1, Strand::Minus), &[50]);
        assert_eq!(s.counts(c1), (2, 1));
    }

    #[test]
    fn a_plain_file_reads_from_byte_zero() {
        // regression: the gzip-magic probe used to leave the cursor at 2, which
        // silently truncated the first chromosome name of every real file
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.bed");
        std::fs::write(&p, b"chrM\t1\t2\tr\t0\t+\nchrM\t3\t4\tr\t0\t+\n").unwrap();
        let mut r = open_maybe_gzip(&p).unwrap();
        skip_leading_headers(&mut r);
        let recs: Vec<_> = SingleEndReader::new(r).map(|x| x.unwrap()).collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].chrom, b"chrM");
        assert_eq!(recs[1].chrom, b"chrM");
    }

    #[test]
    fn a_gzipped_file_round_trips_and_does_not_lose_its_magic() {
        use flate2::write::GzEncoder;
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.bed.gz");
        let raw = b"chr1\t1\t2\tr\t0\t+\nchr1\t5\t6\tr\t0\t+\n";
        let mut enc = GzEncoder::new(Vec::new(), flate2::Compression::fast());
        enc.write_all(raw).unwrap();
        std::fs::write(&p, enc.finish().unwrap()).unwrap();

        assert!(
            is_gzipped(&p),
            "the magic must be detected by content, not name"
        );
        let mut r = open_maybe_gzip(&p).unwrap();
        skip_leading_headers(&mut r);
        let recs: Vec<_> = SingleEndReader::new(r).map(|x| x.unwrap()).collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].chrom, b"chr1");
        assert_eq!(recs[0].pos, 1);
    }

    #[test]
    fn a_gzipped_file_named_without_the_extension_is_still_read() {
        use flate2::write::GzEncoder;
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.bed");
        let mut enc = GzEncoder::new(Vec::new(), flate2::Compression::fast());
        enc.write_all(b"chr1\t7\t8\tr\t0\t+\n").unwrap();
        std::fs::write(&p, enc.finish().unwrap()).unwrap();
        let mut r = open_maybe_gzip(&p).unwrap();
        skip_leading_headers(&mut r);
        let recs: Vec<_> = SingleEndReader::new(r).map(|x| x.unwrap()).collect();
        assert_eq!(recs[0].pos, 7);
    }

    #[test]
    fn reader_skips_negative_and_nameless_records_like_upstream() {
        let mut r = Cursor::new(
            b"chr1\t-5\t20\tr\t0\t+\n\t10\t20\tr\t0\t+\nchr1\t7\t8\tr\t0\t+\n".to_vec(),
        );
        skip_leading_headers(&mut r);
        let mut it = SingleEndReader::new(r);
        let recs: Vec<_> = it.by_ref().map(|x| x.unwrap()).collect();
        assert_eq!(recs.len(), 1, "only the valid record survives");
        assert_eq!(recs[0].pos, 7);
        assert_eq!(it.skipped, 2);
    }
}
