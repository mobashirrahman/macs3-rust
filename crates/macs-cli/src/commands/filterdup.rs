//! `macs3-rs filterdup`: duplicate filtering over BED tags.
//!
//! Reads a BED alignment, optionally limits the number of tags per
//! (position, strand) with `--keep-dup` (`auto` uses the binomial inverse from
//! `macs-stats`), and writes the surviving tags as BED in upstream's
//! `FixWidthTrack.print_to_bed` layout.
//!
//! Port notes (upstream `filterdup_cmd.py:28-96`, `FixWidthTrack.py:490-531`):
//!
//! * `print_to_bed` writes `[pos, pos+fw)` for `+` reads and `[pos-fw, pos)` for
//!   `-` reads, where `pos` is the *stored* coordinate: BED parsing keeps the
//!   start for `+` and the **end** for `-` (`parse_bed_line`), so a minus read
//!   reproduces its input interval.
//! * Per chromosome it emits every `+` record, then every `-` record. Upstream
//!   iterates `get_chr_names()`, a `set`, so its multi-chromosome order is
//!   arbitrary; this port writes chromosomes in name order (deterministic, and
//!   identical for the single-chromosome and small-n cases).

use std::path::PathBuf;

use macs_core::{MacsError, Result, Strand};

use crate::Options;

/// Upstream `cal_max_dup_tags`: `binomial_cdf_inv(1 - p, tags, 1 / gsize)`.
fn cal_max_dup_tags(gsize: f64, tags: u64, p: f64) -> i64 {
    let n = tags as i64;
    let pr = 1.0 / gsize;
    macs_stats::binomial_cdf_inv(1.0 - p, n, pr)
        .unwrap_or(1)
        .max(0)
}

pub fn filterdup(o: &Options) -> Result<()> {
    let ifiles = super::input::input_files(o)?;
    let format = o.get("format").unwrap_or("BED").to_uppercase();
    // Paired-end modes filter duplicate *fragments* (`filterdup_cmd.py:48-74`).
    // FRAG forces `--keep-dup all` upstream (`OptValidator.py:110-112`), so there is
    // nothing to filter and it is rejected here rather than silently accepted.
    if format == "BAMPE" || format == "BEDPE" {
        return filterdup_pe(o, &ifiles, &format);
    }
    if format == "FRAG" {
        return Err(MacsError::InvalidParameter(
            "paired-end input is not yet supported".into(),
        ));
    }
    // SE BAM streams without an index; SAM/BED/ELAND go through text.
    // (`filterdup` rejects BAMPE/BEDPE/FRAG above; BAM here is single-end.)
    // The BAM mean query length is captured for `--tsize` inference below.
    let infer_tsize = matches!(o.int("tsize"), None | Some(0));
    let (mut track, inferred_size) = super::input::load_tag_files(&ifiles, &format, infer_tsize)?;
    let t0 = track.total();

    // tag size fw (upstream sets `inputtrack.fw = options.tsize`)
    //
    // F274: with no `--tsize`, upstream infers it from the first 10 valid reads rather
    // than rejecting. Only an inference that yields nothing (`-1`) is an error.
    let fw = match o.int("tsize") {
        Some(v) if v != 0 => v,
        _ => {
            // BAM infers from query lengths (`BAMParser.tsize`); text infers from
            // the selected text parser. Both sample the first 10 valid records.
            let inferred = inferred_size.trunc() as i64;
            if inferred <= 0 {
                return Err(MacsError::InvalidParameter(
                    "could not infer --tsize from the input: no valid alignments".into(),
                ));
            }
            inferred
        }
    };
    if fw <= 0 {
        return Err(MacsError::InvalidParameter("--tsize must be > 0".into()));
    }

    // F275: upstream's `--keep-dup` defaults to **`auto`**, not `1`, and `-g/--gsize`
    // defaults to `hs` (2,913,022,398) rather than being required:
    //
    //     argparser_filterdup.add_argument("-g", "--gsize", dest="gsize", type=str,
    //                                      default="hs", ...)          # bin/macs3:328
    //     --keep-dup ... default=auto                                     # flag_matrix.tsv
    //
    // so `macs3 filterdup -i reads.bed -f BED -o out.txt` is a valid upstream invocation
    // and it filters with a binomial-derived cap. We rejected it for a missing `-g` and
    // defaulted the cap to `1`, which is a different *result* rather than only a different
    // accept/reject outcome. Verified against upstream: with no `-g` it reports
    // `max_dup_tags based on binomal = 1`, with `-g 24000000` it reports `2`, so the
    // shortcut genuinely participates in the computation.
    let keep = o.get("keepduplicates").unwrap_or("auto").to_lowercase();
    if keep != "all" {
        let max_dup: i64 = if keep == "auto" {
            let spec = o.get("gsize").unwrap_or("hs");
            let gsize =
                macs_core::genomesize::resolve_gsize(spec).map_err(MacsError::InvalidParameter)?;
            let pval = o.float("pvalue").unwrap_or(1e-5);
            cal_max_dup_tags(gsize, t0, pval)
        } else {
            keep.parse()
                .map_err(|_| MacsError::InvalidParameter(format!("invalid --keep-dup {keep}")))?
        };
        let t1 = track.filter_dup(max_dup)?;
        eprintln!("filterdup: {t0} tags, {t1} after filtering (max {max_dup} per position/strand)");
    }

    if o.flag("dryrun") {
        return Ok(());
    }

    // print_to_bed
    // Stream to a buffered writer; the whole-file `String` is hundreds of MB on a
    // real run and made `filterdup` exceed upstream's peak RSS. Bytes unchanged.
    use std::io::Write as _;
    let ofile = o.get("outputfile").unwrap_or("stdout");
    let outdir = PathBuf::from(o.get("outdir").unwrap_or("."));
    let mut file: Option<std::fs::File> = None;
    if ofile != "stdout" {
        std::fs::create_dir_all(&outdir)?;
        file = Some(std::fs::File::create(outdir.join(ofile))?);
    }
    let mut out = std::io::BufWriter::new(match file {
        Some(f) => Box::new(f) as Box<dyn std::io::Write>,
        None => Box::new(std::io::stdout()) as Box<dyn std::io::Write>,
    });
    let pos = track.positions();
    for chrom in pos.chroms_sorted() {
        let name = String::from_utf8_lossy(track.genome().name(chrom));
        for &p in pos.strand(chrom, Strand::Plus) {
            writeln!(out, "{name}\t{p}\t{}\t.\t.\t+", p + fw as u32)?;
        }
        for &p in pos.strand(chrom, Strand::Minus) {
            // Upstream writes negative starts verbatim (`chrIV -11 39`); it does
            // not clamp at zero. `saturating_sub` hid every minus-strand tag within
            // `fw` of the contig start behind a `0`, diverging on real data.
            let lo = i64::from(p) - fw;
            writeln!(out, "{name}\t{lo}\t{p}\t.\t.\t-")?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Filter duplicate paired-end fragments (BEDPE/BAMPE).
///
/// `filterdup_cmd.py:48-74`: load as `PETrackI`, resolve `--keep-dup` (number or
/// binomial `auto`), `filter_dup`, then `print_to_bed` (3-column `chrom start end`).
/// Chromosome order follows upstream's `get_chr_names` set iteration, which is
/// hash-randomized and therefore non-deterministic across runs; within a chromosome
/// `finalize()` sorts by position. This port emits chromosomes sorted and positions
/// in order, which is deterministic and matches upstream on single-chromosome inputs.
/// Multi-chromosome byte-identity is impossible by upstream's nature (F29).
fn filterdup_pe(o: &crate::Options, paths: &[String], format: &str) -> Result<()> {
    use std::collections::BTreeMap;

    let mut frags: BTreeMap<Vec<u8>, Vec<(u32, u32)>> = BTreeMap::new();
    let track = super::input::load_fragment_files(paths, format)?;
    for chrom in track.chroms() {
        let name = track.genome().name(chrom).to_vec();
        for fragment in track.frags(chrom) {
            frags
                .entry(name.clone())
                .or_default()
                .push((fragment.start, fragment.end));
        }
    }

    // Build a FragmentTrack for duplicate filtering.
    let mut b = macs_track::FragTrackBuilder::new();
    for (chrom, spans) in &frags {
        for &(s, e) in spans {
            b.push(chrom, s, e);
        }
    }
    b.finalize();
    let mut track = b.build();
    let t0 = track.total();

    // `--keep-dup`: number, `all` (no-op), or `auto` (binomial).
    let keep = o.get("keepduplicates").unwrap_or("auto");
    let max_dup = if keep == "all" {
        i64::MAX
    } else if keep == "auto" {
        let gsize = match o.get("gsize") {
            Some(spec) => {
                macs_core::genomesize::resolve_gsize(spec).map_err(MacsError::InvalidParameter)?
            }
            None => 2.65e9,
        };
        cal_max_dup_tags(gsize, t0, 1e-5)
    } else {
        keep.parse::<i64>()
            .map_err(|_| MacsError::InvalidParameter(format!("invalid --keep-dup {keep}")))?
    };
    if max_dup != i64::MAX {
        macs_track::filter_frag_dup(&mut track, max_dup)?;
    }
    let t1 = track.total();

    eprintln!(
        "filterdup: {} fragments, {} after filtering (max {} per position)",
        t0,
        t1,
        if max_dup == i64::MAX {
            "unlimited".to_string()
        } else {
            max_dup.to_string()
        }
    );

    // `--dry-run` reports counts without writing (`filterdup_cmd.py:80-85`).
    if o.flag("dryrun") {
        return Ok(());
    }
    let ofile = o.get("outputfile").unwrap_or("filterdup.bedpe");
    let mut out = String::new();
    // Sorted for determinism (see above); `FragmentTrack::chroms` is file order.
    let mut chroms = track.chroms();
    chroms.sort_by(|a, b| track.genome().name(*a).cmp(track.genome().name(*b)));
    for chrom in chroms {
        let name = String::from_utf8_lossy(track.genome().name(chrom));
        let mut spans: Vec<(u32, u32)> = track
            .frags(chrom)
            .iter()
            .map(|f| (f.start, f.end))
            .collect();
        spans.sort_unstable();
        for (s, e) in spans {
            out.push_str(&format!("{name}\t{s}\t{e}\n"));
        }
    }
    if ofile == "stdout" {
        print!("{out}");
    } else {
        let outdir = std::path::PathBuf::from(o.get("outdir").unwrap_or("."));
        std::fs::create_dir_all(&outdir)?;
        std::fs::write(outdir.join(ofile), out)?;
    }
    Ok(())
}
