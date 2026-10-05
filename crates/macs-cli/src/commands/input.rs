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

use macs_core::{MacsError, Result, Strand};
use macs_track::{FragmentTrack, SingleEndTrack};
use std::path::Path;

/// Return all `-i/--ifile` values in their command-line order.
///
/// `Options::get` exposes only the first value for a `nargs +` option; callers
/// that load treatment tracks must use this list so repeated input files are not
/// silently ignored. The fallback keeps manually constructed `Options` values
/// and older test helpers working.
pub fn input_files(options: &crate::Options) -> Result<Vec<String>> {
    let files = options.get_all("ifile");
    let files: Vec<String> = if files.is_empty() {
        options
            .get("ifile")
            .map(str::to_owned)
            .into_iter()
            .collect()
    } else {
        files.to_vec()
    };
    if files.is_empty() {
        return Err(MacsError::InvalidParameter("-i/--ifile is required".into()));
    }
    Ok(files)
}

/// Resolve AUTO from decompressed content. Paired-end formats require an explicit flag.
pub fn single_end_format(path: &Path, format: &str) -> Result<String> {
    use std::io::BufRead;
    if format != "AUTO" {
        return Ok(format.to_string());
    }
    let mut reader = macs_io::open_maybe_gzip(path)?;
    if reader.fill_buf()?.starts_with(b"BAM\x01") {
        return Ok("BAM".into());
    }
    let mut line = Vec::new();
    for _ in 0..10000 {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        while line.last().is_some_and(|b| b.is_ascii_whitespace()) {
            line.pop();
        }
        if line.starts_with(b"@HD\t") || line.starts_with(b"@SQ\t") {
            return Ok("SAM".into());
        }
        if line.is_empty()
            || line.starts_with(b"#")
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
        {
            continue;
        }
        let fields: Vec<_> = line.split(|b| *b == b'\t').collect();
        let integer = |v: &[u8]| {
            std::str::from_utf8(v)
                .ok()
                .and_then(|v| v.trim().parse::<i64>().ok())
                .is_some()
        };
        let format = if fields.len() >= 3 && integer(fields[1]) && integer(fields[2]) {
            "BED"
        } else if fields.len() >= 9
            && matches!(
                fields[2],
                b"U0" | b"U1" | b"U2" | b"R0" | b"R1" | b"R2" | b"NM"
            )
        {
            "ELAND"
        } else if fields.len() >= 4 && fields[2].contains(&b':') {
            "ELANDMULTI"
        } else if fields.len() >= 14 && matches!(fields[13], b"F" | b"R") {
            "ELANDEXPORT"
        } else if fields.len() >= 11 && integer(fields[1]) && integer(fields[3]) {
            "SAM"
        } else if fields.len() >= 5 && matches!(fields[1], b"+" | b"-") {
            "BOWTIE"
        } else {
            continue;
        };
        return Ok(format.into());
    }
    Err(MacsError::Rejected(format!(
        "{}: cannot detect alignment format",
        path.display()
    )))
}

/// Load every supported single-end format and its parser-reported tag length.
pub fn load_single_end(path: &Path, format: &str) -> Result<(SingleEndTrack, f64)> {
    let format = single_end_format(path, format)?;
    match format.as_str() {
        "BAM" => load_single_end_bam(path),
        "SAM" => load_single_end_sam(path),
        "BED" => Ok((
            load_single_end_bed(path)?,
            macs_io::detect_tsize_text(path)?.map_or(-1.0, f64::from),
        )),
        "BOWTIE" | "ELAND" | "ELANDMULTI" | "ELANDEXPORT" => load_legacy(path, &format),
        _ => Err(MacsError::InvalidParameter(format!(
            "unsupported single-end format {format}"
        ))),
    }
}

/// Load multiple single-end files into one track, retaining the first parser's
/// tag-size estimate just like upstream `load_tag_files_options`. With AUTO,
/// upstream invokes `guess_parser` independently for each input path.
pub fn load_single_end_files(paths: &[String], format: &str) -> Result<(SingleEndTrack, f64)> {
    let first = paths
        .first()
        .ok_or_else(|| MacsError::InvalidParameter("-i/--ifile is required".into()))?;
    if paths.len() == 1 {
        return load_single_end(Path::new(first), format);
    }
    let mut builder = macs_track::SingleEndTrackBuilder::new();
    let mut first_size = -1.0;
    for (i, path) in paths.iter().enumerate() {
        let (track, size) = load_single_end(Path::new(path), format)?;
        if i == 0 {
            first_size = size;
        }
        for chrom in track.genome().ids_file_order() {
            let name = track.genome().name(chrom);
            for strand in [Strand::Plus, Strand::Minus] {
                for &pos in track.positions().strand(chrom, strand) {
                    builder.push(name, pos, strand);
                }
            }
        }
    }
    builder.finalize();
    Ok((builder.build(), first_size))
}

/// Upstream commands that call `load_tag_files_options` first ask the first
/// parser for a tag-size estimate before appending later input files. `infer_tsize`
/// distinguishes commands that always infer from filterdup/randsample, which skip
/// this parser call when an explicit `--tsize` is provided.
pub fn load_tag_files(
    paths: &[String],
    format: &str,
    infer_tsize: bool,
) -> Result<(SingleEndTrack, f64)> {
    let first = paths
        .first()
        .ok_or_else(|| MacsError::InvalidParameter("-i/--ifile is required".into()))?;
    let resolved = if infer_tsize {
        Some(single_end_format(Path::new(first), format)?)
    } else {
        None
    };
    let result = load_single_end_files(paths, format)?;
    if infer_tsize && resolved.as_deref() == Some("BOWTIE") {
        validate_bowtie_tsize(Path::new(first))?;
    }
    Ok(result)
}

pub(super) fn validate_bowtie_tsize(path: &Path) -> Result<()> {
    use std::io::BufRead;
    let mut reader = macs_io::open_maybe_gzip(path)?;
    let mut line = Vec::new();
    let (mut valid, mut attempts) = (0usize, 0usize);
    while valid < 10 && attempts < 10_000 {
        line.clear();
        attempts += 1;
        if reader.read_until(b'\n', &mut line)? == 0 {
            // BowtieParser.tlen_parse_line is declared to return cython.int but
            // returns a tuple for EOF, so Cython raises "an integer is required"
            // before build_fwtrack gets a chance to consume short files.
            return Err(MacsError::Rejected(
                "BowtieParser.tlen_parse_line: an integer is required".into(),
            ));
        }
        let fields: Vec<_> = line.trim_ascii_end().split(|b| *b == b'\t').collect();
        if fields.len() <= 4 {
            return Err(MacsError::Rejected(
                "BowtieParser.tlen_parse_line: list index out of range".into(),
            ));
        }
        if !fields[4].is_empty() {
            valid += 1;
        }
    }
    Ok(())
}

/// Load all BEDPE/BAMPE files into one unweighted fragment track.
pub fn load_fragment_files(paths: &[String], format: &str) -> Result<FragmentTrack> {
    if paths.is_empty() {
        return Err(MacsError::InvalidParameter("-i/--ifile is required".into()));
    }
    if !matches!(format, "BEDPE" | "BAMPE") {
        return Err(MacsError::InvalidParameter(format!(
            "unsupported paired-end format {format}"
        )));
    }
    let mut builder = macs_track::FragTrackBuilder::new();
    for path in paths {
        if format == "BAMPE" {
            let (records, _) = macs_io::bam::bampe_fragments(Path::new(path))?;
            for fragment in &records {
                builder.push(
                    &fragment.chrom,
                    fragment.start,
                    fragment.start + fragment.len,
                );
            }
        } else {
            use std::io::BufRead;
            let mut reader = std::io::BufReader::new(macs_io::open_maybe_gzip(Path::new(path))?);
            let mut line = Vec::new();
            loop {
                line.clear();
                if reader.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                let Some(record) = macs_io::parse_bedpe_line(&line)? else {
                    continue;
                };
                if record.chrom.is_empty() || record.left < 0 || record.right < record.left {
                    continue;
                }
                builder.push(&record.chrom, record.left as u32, record.right as u32);
            }
        }
    }
    builder.finalize();
    Ok(builder.build())
}

fn load_legacy(path: &Path, format: &str) -> Result<(SingleEndTrack, f64)> {
    use std::io::BufRead;
    let mut reader = macs_io::open_maybe_gzip(path)?;
    let mut builder = macs_track::SingleEndTrackBuilder::new();
    let mut line = Vec::new();
    let (mut sum, mut n, mut lines) = (0u64, 0u64, 0u64);
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        lines += 1;
        while line.last().is_some_and(|b| b.is_ascii_whitespace()) {
            line.pop();
        }
        if line.is_empty() || line.starts_with(b"#") {
            continue;
        }
        let fields: Vec<_> = line.split(|b| *b == b'\t').collect();
        if format == "ELANDEXPORT" && fields.len() <= 12 {
            continue;
        }
        let get = |i: usize| {
            fields.get(i).copied().ok_or_else(|| {
                MacsError::Rejected(format!("{}: malformed {format} record", path.display()))
            })
        };
        let sequence = get(match format {
            "BOWTIE" => 4,
            "ELANDEXPORT" => 8,
            _ => 1,
        })?;
        let length = if matches!(format, "ELAND" | "ELANDMULTI")
            && sequence.iter().all(u8::is_ascii_digit)
        {
            0
        } else {
            sequence.len() as i64
        };
        let successful = format != "ELANDEXPORT" || fields.get(12).is_some_and(|v| !v.is_empty());
        if n < 10 && lines <= 10000 && length > 0 && successful {
            sum += length as u64;
            n += 1;
        }
        let (chrom, pos, strand) = match format {
            "BOWTIE" => (get(2)?, macs_io::atoi(get(3)?), get(1)?),
            "ELAND" => {
                if fields.len() <= 6 {
                    continue;
                }
                if !matches!(get(2)?, b"U0" | b"U1" | b"U2") {
                    continue;
                }
                (get(6)?, macs_io::atoi(get(7)?) - 1, get(8)?)
            }
            "ELANDMULTI" => {
                if fields.len() < 4 {
                    continue;
                }
                // Pinned MACS3 3.0.5 casts bytes tokens from `fields[2].split(b':')`
                // directly to cython.int before checking whether the row is a
                // multi-hit. Cython therefore raises TypeError("an integer is
                // required") for every otherwise valid ELANDMULTI row.
                return Err(MacsError::Rejected(
                    "ELANDMultiParser.fw_parse_line: an integer is required".into(),
                ));
            }
            _ => {
                if !successful {
                    continue;
                }
                (get(10)?, macs_io::atoi(get(12)?) - 1, get(13)?)
            }
        };
        let strand = match strand {
            b"+" | b"F" => Strand::Plus,
            b"-" | b"R" => Strand::Minus,
            _ => {
                return Err(MacsError::Rejected(format!(
                    "{}: invalid strand",
                    path.display()
                )));
            }
        };
        let chrom = if format == "ELANDEXPORT" {
            chrom
        } else {
            chrom
                .windows(3)
                .rposition(|s| s == b".fa")
                .map_or(chrom, |i| &chrom[..i])
        };
        let pos = pos + if strand == Strand::Minus { length } else { 0 };
        if pos >= 0 && !chrom.is_empty() {
            builder.push(chrom, pos as u32, strand);
        }
    }
    builder.finalize();
    Ok((
        builder.build(),
        sum.checked_div(n).map_or(-1.0, |v| v as f64),
    ))
}

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
        b.push(&rec.chrom, rec.pos as u32, rec.strand);
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
        b.push(&t.chrom, t.pos, strand);
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
        b.push(&rec.chrom, rec.pos as u32, rec.strand);
    }
    b.finalize();
    let mean_qlen = if qn > 0 { qsum as f64 / qn as f64 } else { 0.0 };
    Ok((b.build(), mean_qlen))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_and_auto_formats_preserve_strand_coordinates_and_tag_lengths() {
        let dir = std::env::temp_dir().join(format!("macs-input-formats-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let cases = [
            (
                "BOWTIE",
                "r1\t+\tchr1\t10\tACGT\tIIII\t0\nr2\t-\tchr1\t30\tACGT\tIIII\t0\n",
            ),
            (
                "ELAND",
                "r1\tACGT\tU0\t0\t0\t0\tchr1.fa\t11\tF\nr2\tACGT\tU1\t0\t0\t0\tchr1.fa\t31\tR\n",
            ),
            (
                "ELANDEXPORT",
                "a\tb\tc\td\te\tf\tg\th\tACGT\tIIII\tchr1\tx\t11\tF\na\tb\tc\td\te\tf\tg\th\tACGT\tIIII\tchr1\tx\t31\tR\n",
            ),
        ];
        for (format, text) in cases {
            let file = dir.join(format);
            std::fs::write(&file, text).unwrap();
            for requested in [format, "AUTO"] {
                let (track, size) = load_single_end(&file, requested).unwrap();
                let chrom = track.genome().get(b"chr1").unwrap();
                assert_eq!(track.total(), 2, "{format}/{requested}");
                assert_eq!(size, 4.0);
                assert_eq!(track.positions().strand(chrom, Strand::Plus), &[10]);
                assert_eq!(track.positions().strand(chrom, Strand::Minus), &[34]);
            }
        }
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn unique_elandmulti_rows_match_pinned_cython_type_error() {
        let dir = std::env::temp_dir().join(format!("macs-elandmulti-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("unique.txt");
        std::fs::write(&file, "r1\tACGT\t1:0:0\tchr1.fa:11F0\n").unwrap();
        let error = load_single_end(&file, "ELANDMULTI").expect_err("upstream Cython fails");
        assert!(error.to_string().contains("an integer is required"));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn multihit_elandmulti_rows_match_pinned_cython_type_error() {
        let dir =
            std::env::temp_dir().join(format!("macs-elandmulti-multi-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("multi.txt");
        std::fs::write(&file, "r1\tACGT\t2:0:0\tchr1.fa:11F0\n").unwrap();
        let error = load_single_end(&file, "ELANDMULTI").expect_err("upstream Cython fails");
        assert!(error.to_string().contains("an integer is required"));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn bowtie_tsize_reproduces_the_upstream_short_file_failure_only_when_inferred() {
        let dir = std::env::temp_dir().join(format!("macs-bowtie-size-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("short.bowtie");
        std::fs::write(&file, "r1\t+\tchr1\t10\tACGT\tIIII\t0\t\n").unwrap();
        let paths = vec![file.to_string_lossy().into_owned()];
        let error = load_tag_files(&paths, "BOWTIE", true).expect_err("upstream EOF tuple fails");
        assert!(error.to_string().contains("an integer is required"));
        let (track, _) = load_tag_files(&paths, "BOWTIE", false)
            .expect("explicit --tsize bypasses the upstream tsize() failure");
        assert_eq!(track.total(), 1);
        let long = dir.join("long.bowtie");
        let text = (0..10)
            .map(|i| format!("r{i}\t+\tchr1\t{}\tACGT\tIIII\t0\t\n", i * 10))
            .collect::<String>();
        std::fs::write(&long, text).unwrap();
        let paths = vec![long.to_string_lossy().into_owned()];
        let (track, size) = load_tag_files(&paths, "BOWTIE", true).unwrap();
        assert_eq!(
            track.total(),
            10,
            "the supported 10-tag branch remains enabled"
        );
        assert_eq!(size, 4.0);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn multiple_files_merge_tags_and_keep_first_files_size_estimate() {
        let dir = std::env::temp_dir().join(format!("macs-input-merge-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let first = dir.join("first.bed");
        let second = dir.join("second.bed");
        std::fs::write(&first, "chr1\t10\t20\ta\t0\t+\nchr1\t30\t40\tb\t0\t-\n").unwrap();
        std::fs::write(&second, "chr1\t50\t65\tc\t0\t+\nchr2\t70\t82\td\t0\t-\n").unwrap();
        let paths = vec![
            first.to_string_lossy().into_owned(),
            second.to_string_lossy().into_owned(),
        ];
        let (track, size) = load_tag_files(&paths, "BED", true).unwrap();
        assert_eq!(track.total(), 4);
        assert_eq!(size, 10.0, "size is estimated from the first input only");
        let chr1 = track.genome().get(b"chr1").unwrap();
        let chr2 = track.genome().get(b"chr2").unwrap();
        assert_eq!(track.positions().strand(chr1, Strand::Plus), &[10, 50]);
        assert_eq!(track.positions().strand(chr1, Strand::Minus), &[40]);
        assert_eq!(track.positions().strand(chr2, Strand::Minus), &[82]);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn auto_format_detects_each_file_and_pools_mixed_bed_and_sam() {
        let dir = std::env::temp_dir().join(format!("macs-input-auto-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let first = dir.join("first.bed");
        let second = dir.join("second.sam");
        let bed = (0..20)
            .map(|i| format!("chr1\t{}\t{}\tb{i}\t0\t+\n", 10 + i * 50, 50 + i * 50))
            .collect::<String>();
        let sam = format!(
            "@HD\tVN:1.0\tSO:unsorted\n{}",
            (0..20)
                .map(|i| format!(
                    "s{i}\t0\tchr1\t{}\t60\t40M\t*\t0\t0\t{}\t{}\n",
                    2000 + i * 50,
                    "A".repeat(40),
                    "I".repeat(40)
                ))
                .collect::<String>()
        );
        std::fs::write(&first, bed).unwrap();
        std::fs::write(&second, sam).unwrap();
        let paths = vec![
            first.to_string_lossy().into_owned(),
            second.to_string_lossy().into_owned(),
        ];
        let (track, size) = load_single_end_files(&paths, "AUTO").unwrap();
        assert_eq!(track.total(), 40);
        assert_eq!(size, 40.0, "tag size comes from the first file only");
        std::fs::remove_dir_all(dir).unwrap();
    }
}
