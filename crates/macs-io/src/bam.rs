//! BAM reading: BGZF block decompression, the BAM header, alignment records,
//! and the BAI region index.
//!
//! A port of `MACS3.IO.BAM` (`BAM.py`). It is written directly against the
//! BGZF and BAM byte layouts rather than delegating to a BAM crate, because
//! several of upstream's *filtering* rules are load-bearing and easier to keep
//! faithful when the bytes are parsed here.
//!
//! # BGZF is gzip with a block-size field
//!
//! Each block is a gzip member carrying a `BC` extra subfield whose two-byte
//! payload is the total block size. Upstream reads it with
//!
//! ```python
//! self.bamfile.seek(10, 1)                  # past ID1 ID2 CM FLG MTIME XLEN(2)
//! xlen = unpack("H", self.bamfile.read(2))[0]
//! extra = self.bamfile.read(xlen)
//! bsize = unpack("H", extra[extra.index(b'BC\x02\x00') + 4:])[0]
//! cdata = self.bamfile.read(bsize - xlen - 19)
//! block = decompress(cdata, -MAX_WBITS)     # raw deflate, no wrapper
//! self.bamfile.seek(8, 1)                   # CRC32 + ISIZE trailer
//! ```
//!
//! so `bsize - xlen - 19` is the exact deflate payload length: 12 bytes of fixed
//! header, `xlen` of extra, and 8 of trailer. Reproduced with `flate2`'s raw
//! deflate decoder.
//!
//! # The header is read through plain gzip, not BGZF
//!
//! `__parse_header` opens the file with Python's `gzip`, which concatenates the
//! members, so it sees the uncompressed BAM stream. That matters only in that a
//! truncated final member would raise there and not here.
//!
//! # Coordinate sorting is *checked*, not enforced
//!
//! `__check_sorted` looks at the **first line** of the header text and accepts it
//! only if it starts with `@HD` and the five characters at `SO:` are `coord`. A
//! BAM without an `@HD` line, or with `SO:unsorted`, is rejected outright -- even
//! though the record scan itself would work.
//!
//! # The record filter is four flags, then MAPQ, then a required MD tag
//!
//! In order:
//!
//! 1. skip if `flag` has any of `4` (unmapped), `512` (QC fail), `256`
//!    (secondary) or `2048` (supplementary);
//! 2. if `flag & 1` (paired), additionally require `flag & 2` (proper pair) and
//!    `!(flag & 8)` (mate mapped);
//! 3. skip if `MAPQ < 1` or `MAPQ == 255`;
//! 4. raise `MDTagMissingError` if no `MD` tag is present -- this is the one
//!    hard failure on a well-formed record.
//!
//! Note `callvar` and `callpeak` share this filter, and it is *stricter* than
//! `filterdup`'s, which keeps both reads of a pair.
//!
//! # `rightmost` counts only reference-consuming CIGAR operations
//!
//! ```text
//! for op in cigar: if op & 15 in {0, 2, 3, 7, 8}: rightmost += op >> 4
//! ```
//!
//! i.e. `M`, `D`, `N`, `=`, `X`. `I`, `S`, `H` and `P` advance the query, not the
//! reference, so they are excluded -- which is why a read with a long soft clip
//! gets a *shorter* span here than in a pileup.
//!
//! # Duplicate suppression compares four fields and keeps the first `maxDuplicate`
//!
//! Consecutive reads are duplicates when `(lpos, rpos, strand, cigar)` all match.
//! The counter is reset on any non-duplicate, and reads are kept while
//! `cur_duplicates <= maxDuplicate`, so `--keep-dup N` keeps exactly `N` copies of
//! each distinct alignment in file order.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use flate2::write::DeflateEncoder;
use flate2::Compression;

use macs_core::{MacsError, Result};

/// One parsed alignment, as `MACS3.Signal.ReadAlignment.ReadAlignment`.
///
/// `seq` and `qual` hold the raw 4-bit-packed sequence and the raw Phred scores;
/// [`ReadAlignment::sequence`] decodes the former and
/// [`ReadAlignment::sequence_string`] spells it as ASCII, which is what
/// `RACollection` needs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadAlignment {
    /// QNAME, without the NUL.
    pub name: Vec<u8>,
    /// Reference name, resolved through the header.
    pub chrom: Vec<u8>,
    /// Leftmost aligned position, 0-based inclusive.
    pub lpos: u32,
    /// Rightmost aligned position, exclusive. See the module note.
    pub rpos: u32,
    /// `0` for the forward strand, `1` for reverse.
    pub strand: u8,
    /// 4-bit-packed sequence.
    pub seq: Vec<u8>,
    /// Per-base Phred scores, raw (no `+33`).
    pub qual: Vec<u8>,
    /// Raw CIGAR words, `len << 4 | op`.
    pub cigar: Vec<u32>,
    /// The `MD` tag value.
    pub md: String,
}

impl ReadAlignment {
    /// Decode the 4-bit-packed sequence into a base per position.
    ///
    /// BAM packs two bases per byte, "high nibble first", with `=ACMGRSVTWYHKDBN`
    /// as the code. For an odd `l_seq` the final nibble is padding, encoded as
    /// `=`, and upstream trims it -- so the result is always exactly `l_seq`
    /// bases. That length is load-bearing: `__get_SEQ_QUAL` asserts
    /// `len(seq) == len(qual)`, so a port that emits the padding base trips the
    /// assertion for every odd-length read.
    pub fn sequence(&self) -> Vec<u8> {
        const CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
        let mut out = Vec::with_capacity(self.seq.len() * 2);
        for &b in &self.seq {
            out.push(CODES[(b >> 4) as usize]);
            out.push(CODES[(b & 0x0f) as usize]);
        }
        if out.last() == Some(&b'=') {
            out.pop();
        }
        out
    }

    /// `get_n_edits` (`ReadAlignment.py:99`): CIGAR `I` and `S` lengths, plus
    /// every byte of the `MD` tag in `'A'..='Z'`.
    ///
    /// The MD rule is a raw byte-range test (`c > 64 and c < 91`), so `^NN`
    /// deletions contribute **two** edits and `^` itself contributes none. That
    /// is not what a biologist would write, and it is reproduced because
    /// `RACollection.remove_outliers` filters on it.
    pub fn n_edits(&self) -> u32 {
        let mut n = 0u32;
        for op in &self.cigar {
            if matches!(op & 15, 1 | 4) {
                n += op >> 4;
            }
        }
        for b in self.md.as_bytes() {
            if *b > 64 && *b < 91 {
                n += 1;
            }
        }
        n
    }

    /// The read length, as `len(binaryqual)`.
    pub fn length(&self) -> usize {
        self.qual.len()
    }

    /// The sequence as a printable ASCII string.
    pub fn sequence_string(&self) -> Vec<u8> {
        self.sequence()
    }

    /// One base out of the 4-bit-packed sequence.
    #[inline]
    fn base_at(&self, p: usize) -> u8 {
        const CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
        let byte = self.seq[p / 2];
        let nibble = if p % 2 == 0 { byte >> 4 } else { byte & 0x0f };
        CODES[nibble as usize]
    }

    /// `get_REFSEQ` (`ReadAlignment.py:246`): the read's own bases, edited so that
    /// the result is one base per reference position across `[lpos, rpos)`.
    ///
    /// Two steps, as upstream:
    ///
    /// 1. drop `I` and `S` from the sequence, using the CIGAR;
    /// 2. walk the `MD` tag -- digits accumulate a match run, letters are mismatches
    ///    (overwrite the base, or *insert* when inside a `^NN` deletion), and `^`
    ///    begins a deletion whose bases come from the tag itself.
    ///
    /// The `MD` rules are raw byte-range tests, so a space or any other byte that is
    /// not a digit, `^` or `A..=Z` is an error. Upstream raises there; returning
    /// `Err` keeps the "zero panics on malformed input" contract while still refusing
    /// the read.
    pub fn refseq(&self) -> std::result::Result<Vec<u8>, String> {
        let mut seq = self.sequence();

        // Step 1: remove insertions and soft clips.
        let mut ind = 0usize;
        for &word in &self.cigar {
            let op = word & 15;
            let op_l = (word >> 4) as usize;
            match op {
                // D, H, P consume nothing here; the deletion bases are restored in
                // step 2 from the MD tag.
                2 | 5 | 6 => {}
                0 | 7 | 8 => ind += op_l,
                // I, S -- Python's `seq[ind:ind+l] = b''`, i.e. a splice that
                // *shortens* the sequence.
                1 | 4 => splice(&mut seq, ind, ind + op_l, &[]),
                _ => {}
            }
        }

        // Step 2: walk the MD tag.
        ind = 0;
        let mut run: Vec<u8> = Vec::new();
        let mut in_del = false;
        for &c in self.md.as_bytes() {
            if (48..=57).contains(&c) {
                in_del = false;
                run.push(c);
            } else if (65..=90).contains(&c) && !in_del {
                // A mismatch: skip the preceding match run, overwrite one base.
                ind += digits(&run)?;
                run.clear();
                if ind >= seq.len() {
                    return Err(format!(
                        "MD mismatch at {} is past the {} bases left after step 1",
                        ind,
                        seq.len()
                    ));
                }
                seq[ind] = c;
                ind += 1;
            } else if (65..=90).contains(&c) && in_del {
                // A mismatch *inside* a deletion: insert, shifting the rest right.
                splice(&mut seq, ind, ind, &[c]);
                ind += 1;
            } else if c == b'^' {
                // Deletion: the bases come from the tag, so just skip the run.
                in_del = true;
                ind += digits(&run)?;
                run.clear();
            } else {
                return Err(format!(
                    "Don't understand this operator in MD: {}",
                    c as char
                ));
            }
        }
        Ok(seq)
    }

    /// What one read says about one reference position --
    /// `ReadAlignment.get_variant_bq_by_ref_pos` (`ReadAlignment.py:379`).
    ///
    /// Returns the observed allele, its raw Phred quality, the strand, and whether
    /// the base sits at either end of the read (`tip`).
    ///
    /// # Allele shapes
    ///
    /// * a single matched base (`A`/`C`/`G`/`T`/`N`/...), for `M`/`=`/`X` CIGAR ops;
    /// * `*` with quality **93**, for a position inside a `D` or `N` op -- the
    ///   reference base is deleted;
    /// * the matched base *followed by the inserted bases*, for an `I` op
    ///   immediately after the match that covers the position. Upstream's docstring
    ///   describes this as `^<bases>+` but the code appends them bare, so `AGGG`
    ///   really is the allele string.
    ///
    /// # `pos` and `tip`
    ///
    /// `pos` is the *query* offset of the matched base, and `tip` is
    /// `pos == 0 or pos == len - 1`. `PosReadsInfo` uses `tip` to discount read-end
    /// evidence when counting an alternate allele.
    ///
    /// # `Ok(None)`
    ///
    /// Upstream declares `pos` without initialising it, so a CIGAR that never
    /// matches -- a position in a `D`/`N` that runs off the end of the operator
    /// list, say -- leaves `pos` holding whatever it held before. Rather than
    /// reproduce an uninitialised read, this returns `Ok(None)` and the caller skips
    /// the read, which is what a correct alignment walk would have concluded anyway.
    pub fn variant_bq_by_ref_pos(
        &self,
        ref_pos: u64,
    ) -> std::result::Result<Option<VariantAtPos>, String> {
        if !(self.lpos as u64 <= ref_pos && ref_pos < self.rpos as u64) {
            return Err(format!(
                "position {ref_pos} is outside alignment {} [{}, {})",
                String::from_utf8_lossy(&self.name),
                self.lpos,
                self.rpos
            ));
        }
        let mut res = (ref_pos - self.lpos as u64) as i64;
        let mut p: usize = 0;
        let mut allele: Vec<u8> = Vec::new();
        let mut bq: Vec<u8> = Vec::new();
        let mut pos: Option<usize> = None;

        let mut m = 0usize;
        while m < self.cigar.len() {
            let word = self.cigar[m];
            let op = word & 15;
            let op_l = (word >> 4) as i64;

            match op {
                // M, =, X -- aligned to the reference
                0 | 7 | 8 => {
                    if res < op_l - 1 {
                        p += res as usize;
                        push_base(self, p, &mut allele, &mut bq);
                        pos = Some(p);
                        break;
                    } else if res == op_l - 1 {
                        p += res as usize;
                        push_base(self, p, &mut allele, &mut bq);
                        pos = Some(p);
                        // An `I` op immediately after contributes the inserted bases.
                        if let Some(&next) = self.cigar.get(m + 1) {
                            let n_op = next & 15;
                            let n_len = (next >> 4) as usize;
                            if n_op == 1 {
                                for _ in 0..n_len {
                                    p += 1;
                                    push_base(self, p, &mut allele, &mut bq);
                                }
                            }
                        }
                        break;
                    } else {
                        p += op_l as usize;
                        res -= op_l;
                    }
                }
                // D, N -- present in the reference, absent from the query
                2 | 3 => {
                    if res < op_l {
                        allele.push(b'*');
                        bq.push(93);
                        // Upstream never assigns `pos` on this path. `pos` is only
                        // written inside the M branch, which `break`s, so reaching
                        // here means no match was ever recorded and `pos` reads back
                        // as 0. Every deletion is therefore reported with `pos == 0`
                        // and so is always a `tip`.
                        //
                        // That is very likely an upstream bug, and it is not
                        // cosmetic: `tip` feeds
                        // `PosReadsInfo.update_top_alleles`'s
                        // `n_t[allele] - n_tips[allele]` alt-allele count, so a
                        // deletion allele is discounted as if every read supporting it
                        // were a read end.
                        pos = Some(0);
                        break;
                    } else {
                        res -= op_l;
                    }
                }
                // I -- consumes query only; already handled on the preceding op
                1 => p += op_l as usize,
                // S -- soft clip
                4 => p += op_l as usize,
                // H (5) and P (6) consume neither; upstream has no branch for them.
                _ => {}
            }
            m += 1;
        }

        let Some(pos) = pos else {
            return Ok(None);
        };
        let tip = pos == 0 || pos == self.length().saturating_sub(1);
        Ok(Some(VariantAtPos {
            allele,
            bq,
            strand: self.strand,
            tip,
            pos,
        }))
    }
}

/// One read's contribution at one reference position.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantAtPos {
    /// The observed allele: a base, `*`, or a base followed by inserted bases.
    pub allele: Vec<u8>,
    /// Raw Phred qualities, one per byte of `allele`.
    pub bq: Vec<u8>,
    /// `0` forward, `1` reverse.
    pub strand: u8,
    /// Whether the base is at either end of the read.
    pub tip: bool,
    /// The query offset of the matched base.
    pub pos: usize,
}

#[inline]
fn push_base(ra: &ReadAlignment, p: usize, allele: &mut Vec<u8>, bq: &mut Vec<u8>) {
    allele.push(ra.base_at(p));
    bq.push(ra.qual[p]);
}

/// `int(bytearray_of_ascii_digits)` -- Python's `int()` over a `bytearray`, so a
/// malformed run like `b"1a"` would raise there and must raise here.
fn digits(run: &[u8]) -> std::result::Result<usize, String> {
    if run.is_empty() {
        return Ok(0);
    }
    let s = std::str::from_utf8(run).map_err(|_| format!("bad MD run {run:?}"))?;
    s.parse::<usize>()
        .map_err(|_| format!("Don't understand this operator in MD: {s}"))
}

/// Python's `bytearray` slice assignment: `buf[a:b] = value`.
///
/// When `value` is a different length from `b - a` the sequence **changes length** --
/// it is not a fixed-capacity overwrite. `RACollection.__fill_refseq` relies on that
/// (probably unknowingly), so the port has to reproduce it rather than pad.
fn splice(buf: &mut Vec<u8>, a: usize, b: usize, value: &[u8]) {
    let a = a.min(buf.len());
    let b = b.min(buf.len()).max(a);
    buf.splice(a..b, value.iter().copied());
}

/// The BAM header text and reference list.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BamHeader {
    /// The raw `@`-prefixed header text, newline separated.
    pub text: String,
    /// Reference names, in header order. Index `i` is refID `i`.
    pub references: Vec<Vec<u8>>,
    /// Reference lengths, parallel to `references`.
    pub rlengths: Vec<u32>,
}

impl BamHeader {
    /// Reference index for `name`.
    pub fn ref_index(&self, name: &[u8]) -> Option<usize> {
        self.references.iter().position(|r| r.as_slice() == name)
    }

    /// Reference length for `name`.
    pub fn ref_length(&self, name: &[u8]) -> Option<u32> {
        self.ref_index(name).map(|i| self.rlengths[i])
    }
}

/// `__check_sorted`: the first header line must be `@HD` with `SO:coord`.
pub fn header_is_coordinate_sorted(text: &[u8]) -> bool {
    let first = match text.split(|b| *b == b'\n').next() {
        Some(f) => f,
        None => return false,
    };
    if !first.starts_with(b"@HD") {
        return false;
    }
    match find_sub(first, b"SO:") {
        Some(i) => first.len() >= i + 8 && &first[i + 3..i + 8] == b"coord",
        None => false,
    }
}

fn find_sub(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || hay.len() < needle.len() {
        return None;
    }
    hay.windows(needle.len()).position(|w| w == needle)
}

/// `MDTagMissingError`: a well-formed record with no `MD` tag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MdTagMissingError {
    /// QNAME of the offending record.
    pub read_name: Vec<u8>,
}

impl std::fmt::Display for MdTagMissingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "MD tag missing for read {}",
            String::from_utf8_lossy(&self.read_name)
        )
    }
}

impl std::error::Error for MdTagMissingError {}

/// Parse one BAM alignment block.
///
/// `data` is the record body *after* the 4-byte `block_size` prefix, i.e. starting
/// at `refID`. `min_mapq` defaults to 1, as upstream's `min_MAPQ=1`.
///
/// Returns `Ok(None)` for a record the filter rejects, and an error only for the
/// missing-`MD` case.
#[allow(clippy::too_many_arguments)]
pub fn parse_alignment(
    data: &[u8],
    header: &BamHeader,
    min_mapq: u8,
) -> std::result::Result<Option<ReadAlignment>, MdTagMissingError> {
    parse_alignment_inner(data, header, min_mapq, true)
}

/// As [`parse_alignment`], but without the MAPQ gate.
///
/// Upstream's `bam_fw_binary_parse` (SE BAM 5' ends) never inspects MAPQ, so a
/// `min_mapq` filter here drops reads upstream keeps (167k on yeast, all MAPQ 0/255).
/// `BamAccessor` (callvar) keeps the gated version; tag loaders use this one.
pub fn parse_alignment_nomapq(
    data: &[u8],
    header: &BamHeader,
) -> std::result::Result<Option<ReadAlignment>, MdTagMissingError> {
    parse_alignment_inner(data, header, 0, false)
}

fn parse_alignment_inner(
    data: &[u8],
    header: &BamHeader,
    min_mapq: u8,
    gate_mapq: bool,
) -> std::result::Result<Option<ReadAlignment>, MdTagMissingError> {
    if data.is_empty() {
        return Ok(None);
    }
    if data.len() < 32 {
        // A record too short to hold its own fixed fields cannot be a valid
        // alignment; upstream would raise on the unpack, so treat it as rejected
        // rather than panicking (the zero-panic acceptance criterion).
        return Ok(None);
    }
    let u16at = |o: usize| u16::from_le_bytes([data[o], data[o + 1]]);
    let i32at = |o: usize| i32::from_le_bytes([data[o], data[o + 1], data[o + 2], data[o + 3]]);
    let u32at = |o: usize| u32::from_le_bytes([data[o], data[o + 1], data[o + 2], data[o + 3]]);

    let n_cigar_op = u16at(12) as usize;
    let bwflag = u16at(14);
    let l_read_name = data[8] as usize;
    let mapq = data[9];

    // 1. the four "problematic" flags
    if bwflag & (4 | 512 | 256 | 2048) != 0 {
        return Ok(None);
    }
    // 2. paired reads must be a proper pair with a mapped mate
    if bwflag & 1 != 0 {
        if bwflag & 2 == 0 {
            return Ok(None);
        }
        if bwflag & 8 != 0 {
            return Ok(None);
        }
    }
    // 3. mapping quality (skipped for SE tag loading; upstream never checks it there)
    if gate_mapq && (mapq < min_mapq || mapq == 255) {
        return Ok(None);
    }

    let ref_id = i32at(0);
    let leftmost = u32at(4);
    let l_seq = i32at(16).max(0) as usize;

    let name_end = (32 + l_read_name).min(data.len());
    let mut read_name = data[32..name_end].to_vec();
    if read_name.last() == Some(&0) {
        read_name.pop();
    }

    let mut i = 32 + l_read_name;
    let mut cigar = Vec::with_capacity(n_cigar_op);
    for k in 0..n_cigar_op {
        let o = i + k * 4;
        if o + 4 > data.len() {
            break;
        }
        cigar.push(u32at(o));
    }
    i += n_cigar_op * 4;

    let packed = l_seq.div_ceil(2);
    let seq_end = (i + packed).min(data.len());
    let seq = data[i..seq_end].to_vec();
    i += packed;
    let qual_end = (i + l_seq).min(data.len());
    let qual = data[i.min(data.len())..qual_end].to_vec();
    i += l_seq;

    // `rightmost` sums only the reference-consuming CIGAR ops: M D N = X
    let mut rightmost = leftmost;
    for &op in &cigar {
        if matches!(op & 15, 0 | 2 | 3 | 7 | 8) {
            rightmost = rightmost.wrapping_add(op >> 4);
        }
    }
    let strand = if bwflag & 16 != 0 { 1u8 } else { 0u8 };

    // The MD tag: upstream searches the raw aux bytes for the literal "MDZ",
    // i.e. tag `MD` with type `Z`, then reads to the NUL. Searching for the
    // three-byte run means a value that happens to contain "MDZ" can be
    // mistaken for the tag -- reproduced.
    let tag = &data[i.min(data.len())..];
    let md = match find_sub(tag, b"MDZ") {
        Some(j) => {
            let rest = &tag[j..];
            let end = rest.iter().position(|b| *b == 0).unwrap_or(rest.len());
            String::from_utf8_lossy(&rest[3..end]).into_owned()
        }
        None => {
            return Err(MdTagMissingError {
                read_name: read_name.clone(),
            })
        }
    };

    let chrom = ref_id
        .try_into()
        .ok()
        .and_then(|i: usize| header.references.get(i).cloned())
        .unwrap_or_default();

    Ok(Some(ReadAlignment {
        name: read_name,
        chrom,
        lpos: leftmost,
        rpos: rightmost,
        strand,
        seq,
        qual,
        cigar,
        md,
    }))
}

/// A BGZF block: its uncompressed payload and the compressed offset that follows
/// it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BgzfBlock {
    /// Offset of this block in the file.
    pub coffset: u64,
    /// The decompressed bytes.
    pub data: Vec<u8>,
}

/// The decoded payload of one BGZF block, read from an already-positioned stream.
///
/// Returns the block's total size so the caller can advance. A short read yields
/// `None` (end of file), which is how the record loop terminates.
pub fn read_bgzf_block_at(f: &mut BufReader<File>, coffset: u64) -> Result<Option<BgzfBlock>> {
    f.seek(SeekFrom::Start(coffset))?;
    let mut fixed = [0u8; 12];
    if f.read_exact(&mut fixed).is_err() {
        return Ok(None);
    }
    if fixed[0] != 0x1f || fixed[1] != 0x8b {
        return Err(MacsError::InvalidParameter(format!(
            "not a BGZF block at offset {coffset}: bad gzip magic"
        )));
    }
    let xlen = u16::from_le_bytes([fixed[10], fixed[11]]) as usize;
    let mut extra = vec![0u8; xlen];
    f.read_exact(&mut extra)?;
    let bsize = find_bsize(&extra).ok_or_else(|| {
        MacsError::InvalidParameter(format!("BGZF block at offset {coffset} has no BC subfield"))
    })? as usize;
    // 12 fixed + xlen extra + 8 trailer = the whole member
    let cdata_len = bsize.checked_sub(xlen + 19).ok_or_else(|| {
        MacsError::InvalidParameter(format!("BGZF block at {coffset} has a short BSIZE"))
    })?;
    let mut cdata = vec![0u8; cdata_len];
    f.read_exact(&mut cdata)?;
    // the 8-byte trailer (CRC32 + ISIZE) is not read, only skipped over by the
    // next seek, but it must be consumed for a sequential walk
    let mut trailer = [0u8; 8];
    let _ = f.read_exact(&mut trailer);
    let data = inflate_raw(&cdata)?;
    Ok(Some(BgzfBlock { coffset, data }))
}

/// The `BC` subfield's two-byte block size.
fn find_bsize(extra: &[u8]) -> Option<u16> {
    let mut i = 0usize;
    while i + 4 <= extra.len() {
        let si1 = extra[i];
        let si2 = extra[i + 1];
        let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        if si1 == b'B' && si2 == b'C' && slen == 2 && i + 4 + 2 <= extra.len() {
            return Some(u16::from_le_bytes([extra[i + 4], extra[i + 5]]));
        }
        i += 4 + slen;
    }
    None
}

/// Raw DEFLATE, no zlib or gzip wrapper -- upstream's `decompress(c, -MAX_WBITS)`.
pub fn inflate_raw(data: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    let mut dec = flate2::bufread::DeflateDecoder::new(data);
    std::io::Read::read_to_end(&mut dec, &mut out)
        .map_err(|e| MacsError::InvalidParameter(format!("BGZF inflate failed: {e}")))?;
    Ok(out)
}

/// Read the BAM header (and therefore every reference name and length).
///
/// Uses plain concatenated-gzip decoding, matching upstream's `gzip.open` path.
pub fn read_header(path: &Path) -> Result<BamHeader> {
    let raw = read_all_maybe_gzip(path)?;
    let mut c: &[u8] = &raw;
    if c.len() < 8 || &c[0..4] != b"BAM\x01" {
        return Err(MacsError::InvalidParameter(format!(
            "{} is not a BAM file: the magic string is not \"BAM\\1\"",
            path.display()
        )));
    }
    c = &c[4..];
    let header_len = i32::from_le_bytes([c[0], c[1], c[2], c[3]]) as usize;
    c = &c[4..];
    if c.len() < header_len + 4 {
        return Err(MacsError::InvalidParameter(format!(
            "{}: truncated BAM header",
            path.display()
        )));
    }
    let text = String::from_utf8_lossy(&c[..header_len]).into_owned();
    c = &c[header_len..];
    let nc = i32::from_le_bytes([c[0], c[1], c[2], c[3]]) as usize;
    c = &c[4..];
    let mut references = Vec::with_capacity(nc);
    let mut rlengths = Vec::with_capacity(nc);
    for _ in 0..nc {
        if c.len() < 4 {
            return Err(MacsError::InvalidParameter(format!(
                "{}: truncated BAM reference list",
                path.display()
            )));
        }
        let nlength = i32::from_le_bytes([c[0], c[1], c[2], c[3]]) as usize;
        c = &c[4..];
        if c.len() < nlength + 4 {
            return Err(MacsError::InvalidParameter(format!(
                "{}: truncated BAM reference name",
                path.display()
            )));
        }
        // the trailing NUL is not part of the name
        references.push(c[..nlength.saturating_sub(1)].to_vec());
        c = &c[nlength..];
        rlengths.push(i32::from_le_bytes([c[0], c[1], c[2], c[3]]) as u32);
        c = &c[4..];
    }
    if !header_is_coordinate_sorted(text.as_bytes()) {
        return Err(MacsError::InvalidParameter(
            "BAM should be sorted by coordinates!".into(),
        ));
    }
    Ok(BamHeader {
        text,
        references,
        rlengths,
    })
}

/// Read a whole file, transparently handling gzip and BGZF.
pub fn read_all_maybe_gzip(path: &Path) -> Result<Vec<u8>> {
    use std::io::BufReader as BR;
    let f = File::open(path).map_err(|e| {
        MacsError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })?;
    let mut magic = [0u8; 2];
    let mut probe = BR::new(&f);
    let gz = probe.read_exact(&mut magic).is_ok() && magic == [0x1f, 0x8b];
    // The probe shares the file handle, so the cursor has to go back before the
    // decoder sees the stream -- otherwise it starts two bytes into the first
    // gzip member and reports "unexpected end of file".
    probe.into_inner().seek(SeekFrom::Start(0))?;
    if !gz {
        let mut out = Vec::new();
        let mut r = std::io::BufReader::new(f);
        r.read_to_end(&mut out)?;
        return Ok(out);
    }
    // BGZF is a valid multi-member gzip, so a plain gzip decoder sees the whole
    // concatenated stream.
    let mut out = Vec::new();
    let mut dec = flate2::bufread::MultiGzDecoder::new(BR::new(f));
    dec.read_to_end(&mut out)
        .map_err(|e| MacsError::InvalidParameter(format!("{}: gzip: {e}", path.display())))?;
    Ok(out)
}

/// `reg2bins` (`BAM.py:79`): the BAI bin ids intersecting `[rbeg, rend)`.
///
/// The `assert` upstream is a hard failure; here an out-of-range interval is a
/// usage error, so nothing is silently clamped.
pub fn reg2bins(rbeg: u32, rend: u32) -> Result<Vec<u32>> {
    const MAX_RNG: u32 = (1 << 29) - 1;
    // The first bin id at each level, and its shift, from the BAI spec.
    const BIN_LEVEL_1_START: u32 = 1;
    const BIN_LEVEL_5_START: u32 = 4681;
    if !(rbeg <= rend && rend <= MAX_RNG) {
        return Err(MacsError::InvalidParameter(format!(
            "Invalid region {rbeg}, {rend}"
        )));
    }
    let starts = [0, BIN_LEVEL_1_START, 9, 73, 585, BIN_LEVEL_5_START];
    let shifts = [29u32, 26, 23, 20, 17, 14];
    let mut out = Vec::new();
    for (&start, &shift) in starts.iter().zip(shifts.iter()) {
        let i = if rbeg > 0 { rbeg >> shift } else { 0 };
        let j = if rend < MAX_RNG {
            rend >> shift
        } else {
            MAX_RNG >> shift
        };
        for off in i..=j {
            out.push(start + off);
        }
    }
    Ok(out)
}

/// One BAI chunk: a virtual offset range.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chunk {
    /// Virtual offset of the chunk's first record.
    pub beg: u64,
    /// Virtual offset just past its last record.
    pub end: u64,
}

/// A parsed `.bai` index.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BaiIndex {
    /// Per reference: bin -> chunks, in file order.
    bins: Vec<BTreeMap<u32, Vec<Chunk>>>,
    /// Per reference: linear index (16 kb windows) -> minimum virtual offset.
    linear: Vec<Vec<u64>>,
}

impl BaiIndex {
    /// Parse a BAI file.
    pub fn parse(path: &Path) -> Result<Self> {
        let raw = read_all_maybe_gzip(path)?;
        let mut c: &[u8] = &raw;
        if c.len() < 4 || &c[0..4] != b"BAI\x01" {
            return Err(MacsError::InvalidParameter(format!(
                "{} is not a BAI file: the magic string is not \"BAI\\1\"",
                path.display()
            )));
        }
        c = &c[4..];
        let n_ref = i32::from_le_bytes([c[0], c[1], c[2], c[3]]) as usize;
        c = &c[4..];
        let mut bins = Vec::with_capacity(n_ref);
        let mut linear = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            let n_bin = i32::from_le_bytes([c[0], c[1], c[2], c[3]]) as usize;
            c = &c[4..];
            let mut m = BTreeMap::new();
            for _ in 0..n_bin {
                let bin = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
                c = &c[4..];
                let n_chunk = i32::from_le_bytes([c[0], c[1], c[2], c[3]]) as usize;
                c = &c[4..];
                let mut chunks = Vec::with_capacity(n_chunk);
                for _ in 0..n_chunk {
                    if c.len() < 16 {
                        return Err(MacsError::InvalidParameter(
                            "truncated BAI chunk list".into(),
                        ));
                    }
                    let beg = u64::from_le_bytes(c[0..8].try_into().unwrap());
                    let end = u64::from_le_bytes(c[8..16].try_into().unwrap());
                    c = &c[16..];
                    chunks.push(Chunk { beg, end });
                }
                m.insert(bin, chunks);
            }
            let n_intv = i32::from_le_bytes([c[0], c[1], c[2], c[3]]) as usize;
            c = &c[4..];
            let mut iv = Vec::with_capacity(n_intv);
            for _ in 0..n_intv {
                if c.len() < 8 {
                    return Err(MacsError::InvalidParameter(
                        "truncated BAI linear index".into(),
                    ));
                }
                iv.push(u64::from_le_bytes(c[0..8].try_into().unwrap()));
                c = &c[8..];
            }
            bins.push(m);
            linear.push(iv);
        }
        Ok(Self { bins, linear })
    }

    /// Number of references the index covers.
    pub fn n_refs(&self) -> usize {
        self.bins.len()
    }

    /// Every chunk overlapping `[beg, end)` on reference `ref_n`, in file order.
    pub fn chunks_by_region(&self, ref_n: usize, beg: u32, end: u32) -> Result<Vec<Chunk>> {
        let bins = self
            .bins
            .get(ref_n)
            .ok_or_else(|| MacsError::InvalidParameter("BAI: reference out of range".into()))?;
        let want = reg2bins(beg, end)?;
        // Upstream iterates the reference's bin *dictionary* in insertion order
        // and picks the ones that match, so the result is index order. A BTreeMap
        // gives bin-id order instead, which is the same set in a different
        // order -- harmless here, because only `min(coffset)` is taken.
        let mut out: Vec<Chunk> = Vec::new();
        for bin in want {
            if let Some(cs) = bins.get(&bin) {
                out.extend(cs.iter().copied());
            }
        }
        Ok(out)
    }

    /// `get_coffset_by_region` (`BAM.py:414`): the smallest compressed offset
    /// among the overlapping chunks, or `0` when there are none.
    pub fn coffset_by_region(&self, ref_n: usize, beg: u32, end: u32) -> Result<u64> {
        let chunks = self.chunks_by_region(ref_n, beg, end)?;
        Ok(chunks.iter().map(|c| c.beg >> 16).min().unwrap_or(0))
    }
}

/// A BAM file opened for region queries.
#[derive(Debug)]
pub struct BamAccessor {
    header: BamHeader,
    bai: BaiIndex,
    bam: BufReader<File>,
    coffset_cache: u64,
    block_cache: Vec<u8>,
    coffset_next: u64,
    header_path: String,
}

impl BamAccessor {
    /// Open a BAM and its `.bai`, which is found by appending `.bai` to the BAM
    /// path -- the samtools convention, so `reads.bam` pairs with `reads.bam.bai`.
    pub fn open(path: &Path) -> Result<Self> {
        let header = read_header(path)?;
        let bai_path = {
            let mut s = path.as_os_str().to_os_string();
            s.push(".bai");
            PathBuf::from(s)
        };
        let bai = BaiIndex::parse(&bai_path)?;
        let f = File::open(path).map_err(|e| {
            MacsError::Io(std::io::Error::new(
                e.kind(),
                format!("{}: {e}", path.display()),
            ))
        })?;
        Ok(Self {
            header,
            bai,
            bam: BufReader::new(f),
            coffset_cache: 0,
            block_cache: Vec::new(),
            coffset_next: 0,
            header_path: path.display().to_string(),
        })
    }

    /// The parsed header.
    pub fn header(&self) -> &BamHeader {
        &self.header
    }

    /// Reference names, in header order.
    pub fn chromosomes(&self) -> &[Vec<u8>] {
        &self.header.references
    }

    /// Reference lengths, in header order.
    pub fn rlengths(&self) -> &[u32] {
        &self.header.rlengths
    }

    /// Decompress the block at `coffset`, replacing the cache.
    fn retrieve_block(&mut self, coffset: u64) -> Result<bool> {
        match read_bgzf_block_at(&mut self.bam, coffset)? {
            Some(b) => {
                self.block_cache = b.data;
                Ok(true)
            }
            None => {
                self.block_cache.clear();
                Ok(false)
            }
        }
    }

    /// `get_reads_in_region` (`BAM.py:681`).
    ///
    /// Starts at the BAI-derived block and walks forward, stopping at the first
    /// record whose `lpos > right`. Records are filtered exactly as
    /// [`parse_alignment`] does, and consecutive duplicates beyond
    /// `max_duplicate` are dropped.
    pub fn reads_in_region(
        &mut self,
        chrom: &[u8],
        left: u32,
        right: u32,
        max_duplicate: u32,
    ) -> Result<Vec<ReadAlignment>> {
        let mut out: Vec<ReadAlignment> = Vec::new();
        let Some(ref_n) = self.header.ref_index(chrom) else {
            return Ok(out);
        };
        let coffset = self.bai.coffset_by_region(ref_n, left, right)?;
        if coffset == 0 {
            return Ok(out);
        }
        // Either a new start block, or a cached offset whose payload was
        // consumed by a previous call's block walk.
        if coffset != self.coffset_cache || self.block_cache.is_empty() {
            self.coffset_cache = coffset;
            if !self.retrieve_block(coffset)? {
                return Ok(out);
            }
        }

        let mut cur_duplicates: u32 = 0;
        let mut previous: Option<(u32, u32, u8, Vec<u32>)> = None;
        loop {
            // Scoped borrow: the block is cloned-free (records are copied out)
            // but the cache has to be replaced at the bottom of the loop.
            let block_len = self.block_cache.len();
            let mut i_bytes = 0usize;
            let mut end_searching = false;
            while i_bytes + 4 <= block_len {
                let entrylength = {
                    let b = &self.block_cache;
                    u32::from_le_bytes(b[i_bytes..i_bytes + 4].try_into().unwrap()) as usize
                };
                i_bytes += 4;
                if i_bytes + entrylength > block_len {
                    // a record straddling the block boundary: upstream slices
                    // short and the unpack raises; stopping here keeps the
                    // zero-panic criterion while still yielding the records that
                    // did parse
                    break;
                }
                // A missing MD tag is a hard failure upstream; surfacing it as an
                // error keeps that, and keeps the zero-panic criterion.
                let parsed = {
                    let (rec, hdr) = (
                        &self.block_cache[i_bytes..i_bytes + entrylength],
                        &self.header,
                    );
                    parse_alignment(rec, hdr, 1)
                        .map_err(|e| MacsError::InvalidParameter(e.to_string()))?
                };
                let Some(read) = parsed else {
                    i_bytes += entrylength;
                    continue;
                };
                i_bytes += entrylength;
                if read.lpos > right {
                    end_searching = true;
                    break;
                }
                if read.rpos <= left {
                    continue;
                }
                let sig = (read.lpos, read.rpos, read.strand, read.cigar.clone());
                if previous.as_ref() == Some(&sig) {
                    cur_duplicates += 1;
                } else {
                    cur_duplicates = 1;
                }
                if cur_duplicates <= max_duplicate {
                    out.push(read);
                }
                previous = Some(sig);
            }
            if end_searching {
                break;
            }
            self.coffset_next = self.coffset_cache;
            // advance past the block we just consumed
            let next = next_block_offset(&mut self.bam, self.coffset_cache)?;
            let Some(next) = next else { break };
            self.coffset_cache = next;
            if !self.retrieve_block(next)? {
                break;
            }
        }
        let _ = &self.header_path;
        Ok(out)
    }
}

/// The compressed offset of the block following `coffset`.
///
/// Upstream records `self.bamfile.tell()` *after* decompressing a block -- i.e.
/// the start of the next one -- and then re-reads from there. Reproduced by
/// decompressing the block at `coffset` and taking the resulting stream position,
/// which is what `tell()` would have returned.
fn next_block_offset(f: &mut BufReader<File>, coffset: u64) -> Result<Option<u64>> {
    match read_bgzf_block_at(f, coffset)? {
        Some(_) => Ok(Some(f.stream_position()?)),
        None => Ok(None),
    }
}

/// One fragment as upstream's `BAMPEParser` infers it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BampeFragment {
    /// Reference name, resolved through the header.
    pub chrom: Vec<u8>,
    /// `min(pos, next_pos)`: the leftmost position of the pair.
    pub start: u32,
    /// `abs(tlen)`: the template length, taken from the record's TLEN field
    /// rather than derived from a CIGAR.
    pub len: u32,
}

/// The fragment count and mean template length upstream reports for a BAMPE file.
///
/// `BAMPEParser` sets `self.n` and `self.d = m / i` (`Parser.py:1234-1235`) and
/// `callpeak`/`predictd` print them, so they are part of the observable contract.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BampeSummary {
    /// Number of accepted records (`i`).
    pub n: u64,
    /// Mean of `abs(tlen)` over accepted records (`m / i`).
    pub d: f64,
}

/// Read a BAMPE file as upstream's `BAMPEParser` does.
///
/// # Why this is not `BamAccessor::reads_in_region`
///
/// Upstream does not run a BAM decoder over `-f BAMPE` input. `build_petrack`
/// (`Parser.py:1189-1238`) streams the BGZF stream itself and reads four fixed
/// offsets out of each record body:
///
/// * `flag` at `[14:16]` -- the rejection test,
/// * `refID` at `[0:4]` and `pos` at `[4:8]`,
/// * `next_pos` at `[24:28]` and `tlen` at `[28:32]`.
///
/// Three consequences follow, and all three are observable:
///
/// 1. **The fragment length is `abs(tlen)`**, not the alignment's reference span.
///    A real decoder given the same file and asked for `lpos`/`rpos` reports the
///    *read* length, so a BAMPE run measured `d = 97` where upstream reports
///    `d = 253`, and produced 1040 peaks instead of 724.
/// 2. **Only one record per pair survives.** The filter rejects any paired record
///    whose mate is on the reverse strand (`flag & 128`), so the leftmost,
///    forward record of each proper pair is the one kept -- "we keep only the
///    leftmost position which means this must be at + strand" (`Parser.py:186-188`).
///    A generic decoder yields both mates and must not double-count them.
/// 3. **No MAPQ filter and no `.bai`.** `bampe_pe_binary_parse` never looks at
///    `mapq`, and the sequential stream needs no index.
///
/// # The rejection test
///
/// ```text
/// if (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2)):
///     return (-1, -1, -1)
/// ```
///
/// `2820 = 4 | 256 | 512 | 2048` (unmapped, secondary, QC-fail, supplementary)
/// and `136 = 8 | 128` (mate unmapped, mate reverse).
pub fn bampe_fragments(path: &Path) -> Result<(Vec<BampeFragment>, BampeSummary)> {
    use std::io::Read as _;

    let f = File::open(path).map_err(|e| {
        MacsError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })?;
    // BGZF is multi-member gzip, so a plain gzip decoder sees the whole stream.
    let mut r = flate2::bufread::MultiGzDecoder::new(std::io::BufReader::new(f));

    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)
        .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?;
    if &magic != b"BAM\x01" {
        return Err(MacsError::InvalidParameter(format!(
            "{}: not a BAM file (bad magic)",
            path.display()
        )));
    }
    let mut i32buf = [0u8; 4];
    let mut read_i32 = |r: &mut dyn Read| -> Result<i32> {
        r.read_exact(&mut i32buf)
            .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?;
        Ok(i32::from_le_bytes(i32buf))
    };
    let l_text = read_i32(&mut r)?;
    let mut text = vec![0u8; l_text.max(0) as usize];
    r.read_exact(&mut text)
        .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?;
    // Drop the NUL terminator so the text parser sees the same bytes upstream's does.
    if text.last() == Some(&0) {
        text.pop();
    }
    let n_ref = read_i32(&mut r)?;
    let mut references: Vec<Vec<u8>> = Vec::with_capacity(n_ref.max(0) as usize);
    for _ in 0..n_ref.max(0) {
        let l_name = read_i32(&mut r)?;
        let mut name = vec![0u8; l_name.max(0) as usize];
        r.read_exact(&mut name)
            .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?;
        if name.last() == Some(&0) {
            name.pop();
        }
        references.push(name);
        let _l_ref = read_i32(&mut r)?;
    }

    let mut out: Vec<BampeFragment> = Vec::new();
    let mut m: u64 = 0;
    loop {
        let mut bs = [0u8; 4];
        match r.read_exact(&mut bs) {
            Ok(()) => {}
            Err(_) => break, // clean EOF, or a truncated tail: upstream breaks too
        }
        let block_size = i32::from_le_bytes(bs);
        if block_size <= 0 {
            break;
        }
        let mut body = vec![0u8; block_size as usize];
        r.read_exact(&mut body)
            .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?;
        if body.len() < 32 {
            continue;
        }
        let u16at = |o: usize| u16::from_le_bytes([body[o], body[o + 1]]);
        let i32at = |o: usize| i32::from_le_bytes([body[o], body[o + 1], body[o + 2], body[o + 3]]);
        let bwflag = u16at(14);
        // `bwflag & 2820` or `bwflag & 1 and (bwflag & 136 or not bwflag & 2)`
        if (bwflag & 2820) != 0 || (bwflag & 1 != 0 && ((bwflag & 136) != 0 || (bwflag & 2) == 0)) {
            continue;
        }
        let thisref = i32at(0);
        let pos = i32at(4);
        let nextpos = i32at(24);
        let thistlen = i32at(28);
        let thisstart = pos.min(nextpos).max(0) as u32;
        let thistlen = (thistlen as i64).unsigned_abs() as u32;
        let Some(name) = usize::try_from(thisref)
            .ok()
            .and_then(|i| references.get(i))
            .cloned()
        else {
            continue;
        };
        m += u64::from(thistlen);
        out.push(BampeFragment {
            chrom: name,
            start: thisstart,
            len: thistlen,
        });
    }
    let n = out.len() as u64;
    let d = if n > 0 { m as f64 / n as f64 } else { 0.0 };
    Ok((out, BampeSummary { n, d }))
}

/// One single-end tag from a BAM: chromosome, 5' position, strand, query length.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeBamTag {
    /// Reference name.
    pub chrom: Vec<u8>,
    /// 5' end: `lpos` for `+`, `rpos - 1` for `-` (upstream `ReadAlignment.__str__`
    /// convention, as `callpeak -f BAM` uses it).
    pub pos: u32,
    /// `0` for forward, `1` for reverse.
    pub strand: u8,
    /// Query length (`l_seq`), for `BAMParser.tsize`.
    pub qlen: u32,
}

/// Read a BAM as single-end tags, sequentially, without requiring an index.
///
/// Upstream's `BAMParser` streams the file (`build_fwtrack`) and needs no `.bai`;
/// requiring one here would reject inputs upstream accepts. The record filter is
/// `parse_alignment`'s (unmapped/secondary/QC/supplementary/proper-pair/mate/mapq),
/// **not** the BAMPE leftmost-only filter -- every kept alignment contributes one 5'
/// end, both mates included.
///
/// Returns the tags in stream order plus the mean query length over the first 10
/// valid records (`BAMParser.tsize`, `Parser.py:1024-1050`).
pub fn se_bam_tags(path: &Path) -> Result<(Vec<SeBamTag>, f64)> {
    use std::io::Read as _;

    let f = File::open(path).map_err(|e| {
        MacsError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {e}", path.display()),
        ))
    })?;
    let mut r = flate2::bufread::MultiGzDecoder::new(std::io::BufReader::new(f));

    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)
        .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?;
    if &magic != b"BAM\x01" {
        return Err(MacsError::InvalidParameter(format!(
            "{}: not a BAM file (bad magic)",
            path.display()
        )));
    }
    let mut i32buf = [0u8; 4];
    let mut read_i32 = |r: &mut dyn Read| -> Result<i32> {
        r.read_exact(&mut i32buf)
            .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?;
        Ok(i32::from_le_bytes(i32buf))
    };
    let l_text = read_i32(&mut r)?;
    let mut text = vec![0u8; l_text.max(0) as usize];
    r.read_exact(&mut text)
        .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?;
    if text.last() == Some(&0) {
        text.pop();
    }
    let n_ref = read_i32(&mut r)?;
    let mut references: Vec<Vec<u8>> = Vec::with_capacity(n_ref.max(0) as usize);
    let mut rlengths: Vec<u32> = Vec::with_capacity(n_ref.max(0) as usize);
    for _ in 0..n_ref.max(0) {
        let l_name = read_i32(&mut r)?;
        let mut name = vec![0u8; l_name.max(0) as usize];
        r.read_exact(&mut name)
            .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?;
        if name.last() == Some(&0) {
            name.pop();
        }
        references.push(name);
        rlengths.push(read_i32(&mut r)?.max(0) as u32);
    }
    let header = BamHeader {
        text: String::from_utf8_lossy(&text).into_owned(),
        references,
        rlengths,
    };

    let mut out: Vec<SeBamTag> = Vec::new();
    let (mut qsum, mut qn) = (0u64, 0u64);
    loop {
        let mut bs = [0u8; 4];
        match r.read_exact(&mut bs) {
            Ok(()) => {}
            Err(_) => break,
        }
        let block_size = i32::from_le_bytes(bs);
        if block_size <= 0 {
            break;
        }
        let mut body = vec![0u8; block_size as usize];
        r.read_exact(&mut body)
            .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?;
        if body.len() < 32 {
            continue;
        }
        // Upstream's `bam_fw_binary_parse` keeps only the first mate of a proper
        // pair (`bwflag & 128` rejects the second). `parse_alignment` does not check
        // this -- it serves `BamAccessor` (callvar), where both mates are needed --
        // so filter here. Without it a paired BAM yields ~2x tags (601576 vs 467908
        // on yeast), and the fragment model fits the wrong `d`.
        if body.len() >= 16 {
            let flag = u16::from_le_bytes([body[14], body[15]]);
            if flag & 1 != 0 && flag & 128 != 0 {
                continue;
            }
        }
        let Some(rec) = parse_alignment_nomapq(&body, &header)
            .map_err(|e| MacsError::InvalidParameter(format!("{}: {e}", path.display())))?
        else {
            continue;
        };
        // Query length from CIGAR (M/I/S/=/X consume query). `parse_alignment`
        // does not return `l_seq`, so derive it here for `tsize`.
        let mut qlen = 0u32;
        for &w in &rec.cigar {
            match w & 0xF {
                0 | 1 | 4 | 7 | 8 => qlen += w >> 4,
                _ => {}
            }
        }
        if qn < 10 && qlen > 0 {
            qsum += u64::from(qlen);
            qn += 1;
        }
        // 5' end: `lpos` for `+`, `rpos` (exclusive end) for `-`.
        //
        // Upstream's `bam_fw_binary_parse` returns `pos + sum(M/D/N/=/X)` for minus
        // strand -- the *exclusive* rightmost, not inclusive. Subtracting one here
        // put every minus-strand 5' end a base low (68285 vs 68286), shifting SE
        // pileup breakpoints. `callpeak -f BAM` masked it behind a 200 bp extension.
        let (pos, strand) = if rec.strand == 1 {
            (rec.rpos, 1u8)
        } else {
            (rec.lpos, 0u8)
        };
        out.push(SeBamTag {
            chrom: rec.chrom,
            pos,
            strand,
            qlen,
        });
    }
    let mean_qlen = if qn > 0 { qsum as f64 / qn as f64 } else { 0.0 };
    Ok((out, mean_qlen))
}

/// Compress `data` into BGZF blocks, for writing test fixtures.
///
/// This is a fixture-authoring helper, not part of the read path; it exists so
/// the differential harness can produce BAMs without shelling out to samtools.
pub fn bgzf_compress(data: &[u8], block_size: usize) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    for chunk in data.chunks(block_size.max(1)) {
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
        std::io::Write::write_all(&mut enc, chunk)?;
        let cdata = enc.finish()?;
        // The `BC` field holds the member size **minus one**: 12 bytes of fixed
        // header + 6 of `BC` extra + the payload + an 8-byte CRC/ISIZE trailer,
        // less one. That off-by-one is why the reader computes the payload as
        // `bsize - xlen - 19` and can then simply skip 8 bytes to land on the
        // next block's start.
        let bsize = cdata.len() + 25;
        if bsize > 65536 {
            return Err(MacsError::InvalidParameter(
                "BGZF block payload too large".into(),
            ));
        }
        out.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04]); // ID1 ID2 CM FLG
        out.extend_from_slice(&[0, 0, 0, 0]); // MTIME
        out.extend_from_slice(&[0x00, 0xff]); // XFL, OS
        out.extend_from_slice(&6u16.to_le_bytes()); // XLEN
        out.extend_from_slice(b"BC");
        out.extend_from_slice(&2u16.to_le_bytes());
        out.extend_from_slice(&(bsize as u16).to_le_bytes());
        out.extend_from_slice(&cdata);
        out.extend_from_slice(&crc32(chunk).to_le_bytes());
        out.extend_from_slice(&(chunk.len() as u32).to_le_bytes());
    }
    Ok(out)
}

fn crc32(data: &[u8]) -> u32 {
    // BGZF's trailer carries the standard zlib CRC-32 of the *uncompressed* data.
    let mut crc = flate2::Crc::new();
    crc.update(data);
    crc.sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(refs: &[(&[u8], u32)]) -> BamHeader {
        BamHeader {
            text: String::from("@HD\tVN:1.6\tSO:coordinate\n"),
            references: refs.iter().map(|(n, _)| n.to_vec()).collect(),
            rlengths: refs.iter().map(|(_, l)| *l).collect(),
        }
    }

    /// Build a BAM alignment block body (everything after `block_size`).
    #[allow(clippy::too_many_arguments)]
    fn record(
        ref_id: i32,
        lpos: u32,
        flag: u16,
        mapq: u8,
        name: &[u8],
        cigar: &[u32],
        seq: &[u8],
        qual: &[u8],
        md: Option<&[u8]>,
    ) -> Vec<u8> {
        let mut b = Vec::new();
        b.extend_from_slice(&ref_id.to_le_bytes());
        b.extend_from_slice(&lpos.to_le_bytes());
        let mut n = name.to_vec();
        n.push(0);
        b.push(n.len() as u8);
        b.push(mapq);
        b.extend_from_slice(&0u16.to_le_bytes()); // bin
        b.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
        b.extend_from_slice(&flag.to_le_bytes());
        b.extend_from_slice(&(seq.len() as i32).to_le_bytes());
        b.extend_from_slice(&(-1i32).to_le_bytes()); // next_refID
        b.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos
        b.extend_from_slice(&0i32.to_le_bytes()); // tlen
        b.extend_from_slice(&n);
        for op in cigar {
            b.extend_from_slice(&op.to_le_bytes());
        }
        let packed: Vec<u8> = seq
            .chunks(2)
            .map(|p| {
                let hi = code(p[0]) << 4;
                let lo = if p.len() > 1 { code(p[1]) } else { 0 };
                hi | lo
            })
            .collect();
        b.extend_from_slice(&packed);
        b.extend_from_slice(qual);
        if let Some(m) = md {
            b.extend_from_slice(b"MDZ");
            b.extend_from_slice(m);
            b.push(0);
        }
        b.extend_from_slice(b"NMC");
        b.push(b'i');
        b.extend_from_slice(&1i32.to_le_bytes());
        b
    }

    fn code(b: u8) -> u8 {
        b"=ACMGRSVTWYHKDBN"
            .iter()
            .position(|c| *c == b)
            .unwrap_or(15) as u8
    }

    #[test]
    fn sorted_header_is_detected_only_from_the_first_line() {
        assert!(header_is_coordinate_sorted(
            b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\n"
        ));
        assert!(!header_is_coordinate_sorted(b"@HD\tVN:1.6\tSO:queryname\n"));
        // not the first line -> rejected
        assert!(!header_is_coordinate_sorted(
            b"@SQ\tSN:chr1\n@HD\tSO:coordinate\n"
        ));
        // no @HD at all -> rejected
        assert!(!header_is_coordinate_sorted(b"@SQ\tSN:chr1\n"));
        // only the first five characters after `SO:` are compared, so a longer
        // sort order that starts with "coord" is accepted -- upstream compares
        // `tl[i+3:i+8]` and nothing else
        assert!(header_is_coordinate_sorted(b"@HD\tSO:coordinateX\n"));
        assert!(!header_is_coordinate_sorted(b"@HD\tSO:co\n"));
    }

    #[test]
    fn a_plain_forward_read_parses() {
        let h = header(&[(b"chr1", 1000)]);
        // cigar 10M, seq ACGT, mapq 30
        let r = record(
            0,
            100,
            0,
            30,
            b"r1",
            &[10 << 4],
            b"ACGT",
            &[30; 4],
            Some(b"10"),
        );
        let got = parse_alignment(&r, &h, 1).unwrap().unwrap();
        assert_eq!(got.name, b"r1");
        assert_eq!(got.chrom, b"chr1");
        assert_eq!(got.lpos, 100);
        assert_eq!(got.rpos, 110);
        assert_eq!(got.strand, 0);
        assert_eq!(got.sequence_string(), b"ACGT");
        assert_eq!(got.qual, vec![30, 30, 30, 30]);
        assert_eq!(got.md, "10");
    }

    #[test]
    fn the_four_problematic_flags_are_rejected() {
        let h = header(&[(b"chr1", 1000)]);
        for flag in [4u16, 512, 256, 2048] {
            let r = record(
                0,
                100,
                flag,
                30,
                b"r",
                &[10 << 4],
                b"ACGT",
                &[30; 4],
                Some(b"10"),
            );
            assert!(
                parse_alignment(&r, &h, 1).unwrap().is_none(),
                "flag {flag} should be rejected"
            );
        }
    }

    #[test]
    fn paired_reads_must_be_a_proper_pair_with_a_mapped_mate() {
        let h = header(&[(b"chr1", 1000)]);
        let mk = |flag: u16| {
            record(
                0,
                100,
                flag,
                30,
                b"r",
                &[10 << 4],
                b"ACGT",
                &[30; 4],
                Some(b"10"),
            )
        };
        // 1|4: paired but mate unmapped -> reject
        assert!(parse_alignment(&mk(1 | 4), &h, 1).unwrap().is_none());
        // 1: paired, not proper -> reject
        assert!(parse_alignment(&mk(1), &h, 1).unwrap().is_none());
        // 1|2|64: proper pair, forward first -> keep
        assert!(parse_alignment(&mk(1 | 2 | 64), &h, 1).unwrap().is_some());
        // 1|2|16: proper pair, reverse -> keep, strand 1
        assert_eq!(
            parse_alignment(&mk(1 | 2 | 16), &h, 1)
                .unwrap()
                .unwrap()
                .strand,
            1
        );
    }

    #[test]
    fn mapq_zero_and_255_are_rejected() {
        let h = header(&[(b"chr1", 1000)]);
        for mapq in [0u8, 255] {
            let r = record(
                0,
                100,
                0,
                mapq,
                b"r",
                &[10 << 4],
                b"ACGT",
                &[30; 4],
                Some(b"10"),
            );
            assert!(parse_alignment(&r, &h, 1).unwrap().is_none());
        }
        let r = record(
            0,
            100,
            0,
            1,
            b"r",
            &[10 << 4],
            b"ACGT",
            &[30; 4],
            Some(b"10"),
        );
        assert!(
            parse_alignment(&r, &h, 1).unwrap().is_some(),
            "MAPQ 1 is kept"
        );
    }

    #[test]
    fn rightmost_counts_only_reference_consuming_cigar_ops() {
        let h = header(&[(b"chr1", 1000)]);
        // 4M2I4M: the I advances the query only -> rpos = 108, not 110
        let r = record(
            0,
            100,
            0,
            30,
            b"r",
            &[(4 << 4), (2 << 4 | 1), (4 << 4)],
            b"ACGTACGT",
            &[30; 8],
            Some(b"8"),
        );
        let got = parse_alignment(&r, &h, 1).unwrap().unwrap();
        assert_eq!(got.rpos, 108);
        // 4M2D4M: D consumes the reference -> 110
        let r = record(
            0,
            100,
            0,
            30,
            b"r",
            &[(4 << 4), (2 << 4 | 2), (4 << 4)],
            b"ACGT",
            &[30; 4],
            Some(b"10"),
        );
        assert_eq!(parse_alignment(&r, &h, 1).unwrap().unwrap().rpos, 110);
    }

    #[test]
    fn a_missing_md_tag_is_an_error_not_a_silent_drop() {
        let h = header(&[(b"chr1", 1000)]);
        let r = record(0, 100, 0, 30, b"r", &[4 << 4], b"ACGT", &[30; 4], None);
        match parse_alignment(&r, &h, 1) {
            Err(MdTagMissingError { read_name }) => assert_eq!(read_name, b"r"),
            other => panic!("expected MdTagMissingError, got {other:?}"),
        }
    }

    #[test]
    fn a_truncated_record_is_rejected_rather_than_panicking() {
        let h = header(&[(b"chr1", 1000)]);
        for n in 0..20usize {
            let r = vec![0u8; n];
            assert!(parse_alignment(&r, &h, 1).is_ok(), "len {n} must not panic");
        }
    }

    #[test]
    fn reg2bins_covers_the_region() {
        let bins = reg2bins(0, 1).unwrap();
        assert!(bins.contains(&0));
        let bins = reg2bins(1000, 2000).unwrap();
        // the deepest level starts at bin 4681 (shift 14); the shift-26 level
        // starts at bin 1. 1000 >> 14 and 1000 >> 26 are both 0, so these
        // assert the level-0 ids are present.
        assert!(bins.contains(&4681u32));
        assert!(bins.contains(&1u32));
        assert!(reg2bins(10, 5).is_err());
        assert!(reg2bins(0, 1 << 29).is_err());
    }

    #[test]
    fn bgzf_round_trips() {
        let data: Vec<u8> = (0..5000u32).map(|i| (i % 251) as u8).collect();
        let comp = bgzf_compress(&data, 1024).unwrap();
        let dir = std::env::temp_dir().join(format!("macs3rs_bgzf_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.bgzf");
        std::fs::write(&path, &comp).unwrap();
        let f = File::open(&path).unwrap();
        let mut r = BufReader::new(f);
        let mut all = Vec::new();
        let mut off = 0u64;
        while let Some(b) = read_bgzf_block_at(&mut r, off).unwrap() {
            all.extend_from_slice(&b.data);
            off = r.stream_position().unwrap();
        }
        assert_eq!(all, data);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn bsize_is_found_among_other_extra_subfields() {
        let mut extra = vec![b'X', b'Y', 2, 0, 0xaa, 0xbb];
        extra.extend_from_slice(b"BC");
        extra.extend_from_slice(&2u16.to_le_bytes());
        extra.extend_from_slice(&1234u16.to_le_bytes());
        assert_eq!(find_bsize(&extra), Some(1234));
        assert_eq!(
            find_bsize(b"BC\x03\x00\x00\x00"),
            None,
            "wrong subfield length"
        );
        assert_eq!(find_bsize(b"XX"), None);
    }
}
