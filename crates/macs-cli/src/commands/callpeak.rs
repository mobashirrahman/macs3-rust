//! `macs3-rs callpeak`: the thin CLI shim that feeds the shared library
//! pipeline. The flag surface is validated by the derived matrix parser
//! (`crate::parse_flags`); this module only maps parsed destinations onto
//! [`macs_peaks::callpeak`] and writes the output files upstream writes.
//!
//! # Input formats and which pipeline they select
//!
//! | `-f` | reader | pipeline |
//! |---|---|---|
//! | `AUTO`, `BED`, `BOWTIE`, `ELAND*`, `SAM` | BED / SAM text | single-end |
//! | `BAM` | BGZF + BAM | single-end |
//! | `BAMPE` | BGZF + BAM | paired-end |
//! | `BEDPE` | three-column BEDPE | paired-end |
//! | `FRAG` | five-column FRAG | paired-end, counted |
//!
//! Paired-end mode piles whole fragments for the treatment and both fragment
//! *ends* for the control, and reports `d` as the truncated mean template length
//! **as read**. See [`macs_peaks::callpeak::run_callpeak_pe`] for the scale
//! factors and the doubled control count.
//!
//! # Duplicate filtering
//!
//! `--keep-dup 1` (the default) keeps one tag per position per strand for
//! single-end, and one per `(start, end)` pair for paired-end. `--keep-dup all`
//! and `--keep-dup auto` do no filtering here: `auto` needs MACS's own binomial
//! model and is rejected rather than approximated, so the caller is told instead
//! of getting a subtly different answer.
//!
//! # `d` is captured before filtering
//!
//! `options.tsize` is the mean template length of the track as read, and
//! `options.d` is that truncated -- *not* the post-filter mean. The `gmini_mpe`
//! fixtures have 1200 input fragments truncating to 145, 149 and 180, which are
//! exactly the `# d =` values upstream prints, while the retained sets truncate
//! to 122, 148 and 180. So the mean is captured before `filter_dup`.

use std::collections::{BTreeMap, BTreeSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};

use macs_core::genome::ChromId;
use macs_core::{MacsError, Result, Strand};
use macs_io::peakout::{narrowpeak_body, summits_body, ScoreCol, XlsRow};
use macs_peaks::callpeak::{call_chromosome, ChromCall, PeConfig, SeConfig};
use macs_score::PScoreCache;
use macs_track::{
    filter_frag_dup, FragTrackBuilder, FragmentTrack, SingleEndTrack, SingleEndTrackBuilder,
};

use crate::commands::stagedump::{self, StageDump};
use crate::Options;

const VERSION: &str = "3.0.5";

/// The `# ARGUMENTS LIST:` block `opt_validate_callpeak` builds
/// (`OptValidator.py:174-256`).
///
/// Every line is conditional on the flags in force, and the order is
/// load-bearing because the whole block is byte-compared. Two details are easy
/// to get wrong:
///
/// * the list fields are Python `list` reprs, so a single `-t` prints
///   `['path']` with single quotes;
/// * `--pvalue` replaces the `qvalue cutoff` line *and* adds a
///   "qvalue will not be calculated and reported as -1" line;
/// * the lambda line has **one** window when there is no control.
fn argtxt(o: &Options, name: &str, gsize: f64, pe: bool) -> String {
    let mut a = String::new();
    a.push_str(&format!(
        "# Command line: {}
",
        o.raw_argv.join(" ")
    ));
    a.push_str(
        "# ARGUMENTS LIST:
",
    );
    a.push_str(&format!(
        "# name = {name}
"
    ));
    a.push_str(&format!(
        "# format = {}\n",
        o.get("format").unwrap_or("AUTO")
    ));
    a.push_str(&format!(
        "# ChIP-seq file = {}\n",
        py_list(o.get_all("tfile"))
    ));
    // `-c` is `nargs="*"` with **no** default, so an absent control
    // interpolates as the Python object `None`, not an empty list.
    a.push_str(&format!(
        "# control file = {}\n",
        if o.get_all("cfile").is_empty() {
            "None".to_string()
        } else {
            py_list(o.get_all("cfile"))
        }
    ));
    a.push_str(&format!("# effective genome size = {}\n", py_e(gsize, 2)));
    a.push_str(&format!("# band width = {}\n", o.int("bw").unwrap_or(300)));
    // `mfold` is `nargs=2` with `type=int`, so it interpolates as a Python list
    // repr: `[5, 50]`.
    let mfold: Vec<i64> = o
        .get_all("mfold")
        .iter()
        .filter_map(|x| x.parse().ok())
        .collect();
    // `# model fold = [3, 20]` -- a Python list repr of the two bounds. The
    // `{:?}` of a `Vec<f64>` is `[3.0, 20.0]`, so both the element formatting and
    // the bracket nesting have to be corrected: printing the `Vec` *inside*
    // brackets gives `[[3, 20]]`, which is a list of one list.
    let mfold = if mfold.is_empty() {
        "[5, 50]".to_string()
    } else {
        format!(
            "[{}]",
            mfold
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )
    };
    a.push_str(&format!("# model fold = {mfold}\n"));
    // the trailing newline of the model-fold line is the "\n" in its own
    // format string, which is why the block already ends in a blank line
    let pvalue = o.float("pvalue");
    let broad = o.flag("broad");
    match pvalue {
        Some(p) if broad => {
            a.push_str(&format!(
                "# pvalue cutoff for narrow/strong regions = {}\n",
                py_e(p, 2)
            ));
            a.push_str(&format!(
                "# pvalue cutoff for broad/weak regions = {}\n",
                py_e(o.float("broadcutoff").unwrap_or(0.1), 2)
            ));
            a.push_str("# qvalue will not be calculated and reported as -1 in the final output.\n");
        }
        Some(p) => {
            a.push_str(&format!("# pvalue cutoff = {}\n", py_e(p, 2)));
            a.push_str("# qvalue will not be calculated and reported as -1 in the final output.\n");
        }
        None if broad => {
            a.push_str(&format!(
                "# qvalue cutoff for narrow/strong regions = {}\n",
                py_e(o.float("qvalue").unwrap_or(0.05), 2)
            ));
            a.push_str(&format!(
                "# qvalue cutoff for broad/weak regions = {}\n",
                py_e(o.float("broadcutoff").unwrap_or(0.1), 2)
            ));
        }
        None => {
            a.push_str(&format!(
                "# qvalue cutoff = {}\n",
                py_e(o.float("qvalue").unwrap_or(0.05), 2)
            ));
        }
    }
    match o.int("maxgap") {
        Some(g) if g > 0 => a.push_str(&format!("# The maximum gap between significant sites = {g}\n")),
        _ => a.push_str(
            "# The maximum gap between significant sites is assigned as the read length/tag size.\n",
        ),
    }
    match o.int("minlen") {
        Some(l) if l > 0 => a.push_str(&format!("# The minimum length of peaks = {l}\n")),
        _ => a.push_str(
            "# The minimum length of peaks is assigned as the predicted fragment length \"d\".\n",
        ),
    }
    if o.flag("downsample") {
        a.push_str("# Larger dataset will be randomly sampled towards smaller dataset.\n");
        if let Some(seed) = o.int("seed") {
            if seed >= 0 {
                a.push_str(&format!("# Random seed has been set as: {seed}\n"));
            }
        }
    } else if o.get("scaleto") == Some("large") {
        a.push_str("# Smaller dataset will be scaled towards larger dataset.\n");
    } else {
        a.push_str("# Larger dataset will be scaled towards smaller dataset.\n");
    }
    let ratio = o.float("ratio").unwrap_or(1.0);
    if ratio != 1.0 {
        a.push_str(&format!(
            "# Using a custom scaling factor: {}\n",
            py_e(ratio, 2)
        ));
    }
    if o.get_all("cfile").is_empty() {
        a.push_str(&format!(
            "# Range for calculating regional lambda is: {} bps\n",
            o.int("largelocal").unwrap_or(10000)
        ));
    } else {
        a.push_str(&format!(
            "# Range for calculating regional lambda is: {} bps and {} bps\n",
            o.int("smalllocal").unwrap_or(1000),
            o.int("largelocal").unwrap_or(10000)
        ));
    }
    a.push_str(if broad {
        "# Broad region calling is on\n"
    } else {
        "# Broad region calling is off\n"
    });
    let fecutoff = o.float("fecutoff").unwrap_or(1.0);
    if fecutoff != 1.0 {
        a.push_str(&format!(
            "# Additional cutoff on fold-enrichment is: {fecutoff:.2}\n"
        ));
    }
    if pe {
        // `opt_validate_callpeak` also sets `options.shift = 0` here, which is
        // why the xls has no "shifted by" line in PE mode even with `--shift`.
        a.push_str("# Paired-End mode is on\n");
    } else {
        a.push_str("# Paired-End mode is off\n");
    }
    if o.flag("call_summits") {
        a.push_str("# Searching for subpeak summits is on\n");
    }
    a
}

/// Python's `"%.Ne"`: a signed two-digit exponent.
///
/// Rust's `{:.2e}` writes `4.82e4` where Python writes `4.82e+04` -- Python
/// always emits the sign and pads the exponent to at least two digits. The
/// xls header carries three such fields (`effective genome size`, the q/p value
/// cutoffs), so the difference shows up on every run.
/// The per-chromosome treatment pileups, for the stage dump.
fn signal_tracks(
    sigs: &[macs_peaks::callpeak::ChromSignals],
) -> Vec<(String, &macs_rle::SignalTrack<f32>)> {
    sigs.iter().map(|c| (c.name.clone(), &c.treat)).collect()
}

/// Resident set size in kB, for `MACS3_RS_RSS_TRACE` stage-by-stage attribution.
pub fn rss_kb() -> u64 {
    std::fs::read_to_string("/proc/self/statm")
        .ok()
        .and_then(|t| {
            t.split_whitespace()
                .nth(1)
                .and_then(|p| p.parse::<u64>().ok())
        })
        .map(|pages| pages * 4)
        .unwrap_or(0)
}

fn py_e(v: f64, prec: usize) -> String {
    let s = format!("{v:.*e}", prec);
    // Rust: "4.82e4" / "1e-9"; Python: "4.82e+04" / "1.00e-09"
    match s.split_once('e') {
        None => s,
        Some((mant, exp)) => {
            let (sign, digits) = match exp.strip_prefix('-') {
                Some(d) => ('-', d),
                None => ('+', exp),
            };
            if digits.len() < 2 {
                format!("{mant}e{sign}0{digits}")
            } else {
                format!("{mant}e{sign}{digits}")
            }
        }
    }
}

/// A Python `list` repr of strings, as `options.tfile` interpolates.
fn py_list(v: &[String]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|x| format!("'{x}'"))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

/// The per-track counts block (`callpeak_cmd.py:79-154`).
///
/// `tag` is `"fragment"` in PE mode and `"tag"` otherwise, and the
/// "maximum duplicate" line differs too. `Redundant rate` is `(t0 - t1) / t0`
/// printed with `%.2f`, so a track with zero input reads `0.00` rather than
/// dividing by zero.
fn tagsinfo(
    pe: bool,
    tsize: i64,
    t0: u64,
    t1: u64,
    max_dup: Option<i64>,
    control: Option<(u64, u64)>,
) -> String {
    let tag = if pe { "fragment" } else { "tag" };
    let mut s = format!("# {tag} size is determined as {tsize} bps\n");
    s.push_str(&format!("# total {tag}s in treatment: {t0}\n"));
    if let Some(md) = max_dup {
        s.push_str(&format!("# {tag}s after filtering in treatment: {t1}\n"));
        if pe {
            s.push_str(&format!(
                "# maximum duplicate fragments in treatment = {md}\n"
            ));
        } else {
            s.push_str(&format!(
                "# maximum duplicate tags at the same position in treatment = {md}\n"
            ));
        }
        let rate = if t0 == 0 {
            0.0
        } else {
            (t0 - t1) as f64 / t0 as f64
        };
        s.push_str(&format!("# Redundant rate in treatment: {rate:.2}\n"));
    }
    if let Some((c0, c1)) = control {
        s.push_str(&format!("# total {tag}s in control: {c0}\n"));
        if let Some(md) = max_dup {
            s.push_str(&format!("# {tag}s after filtering in control: {c1}\n"));
            if pe {
                s.push_str(&format!(
                    "# maximum duplicate fragments in control = {md}\n"
                ));
            } else {
                s.push_str(&format!(
                    "# maximum duplicate tags at the same position in control = {md}\n"
                ));
            }
            let rate = if c0 == 0 {
                0.0
            } else {
                (c0 - c1) as f64 / c0 as f64
            };
            s.push_str(&format!("# Redundant rate in control: {rate:.2}\n"));
        }
    }
    s
}

/// Formats that select the paired-end pipeline.
const PE_FORMATS: [&str; 3] = ["BAMPE", "BEDPE", "FRAG"];
/// Formats that need the BGZF/BAM reader.
const BAM_FORMATS: [&str; 2] = ["BAM", "BAMPE"];

/// Read a BED/SAM-ish single-end file.
///
/// Tab-separated only, matching upstream: a space-separated BED is *not*
/// accepted, because `split(b"\t")` is the whole parser and a file with spaces
/// yields one field whose `pos` does not parse.
/// The measured read length, as upstream's `GenericParser.get_tag_size` computes it.
///
/// F153: upstream does **not** average the whole file. `get_tag_size`
/// (`Parser.py:275-293`) reads forward until it has ten alignments with a
/// positive length:
///
/// ```python
/// while n < 10 and m < 10000:
///     m += 1
///     thisline = self.fhd.readline()
///     this_taglength = self.tlen_parse_line(thisline)
///     if this_taglength > 0:
///         s += this_taglength
///         n += 1
/// self.tag_size = cython.cast(cython.int, (s/n))
/// ```
///
/// then rewinds. So `options.tsize` -- and therefore the xls header's
/// `# tag size is determined as N bps` -- is the **truncated mean of the first ten
/// tags**: 190 on `gmini_mse_d1200_w60_ctrl_009`, whose whole-file mean is 151.
/// Averaging everything, or reporting `--extsize`, gets both wrong. Note this is
/// the *measured* size: `--extsize` still overrides `d` for the pileup extension
/// afterwards, which is why a 190 bp report pairs with a 200 bp extension.
///
/// Returned alongside the track so the caller can put it in the header.
fn load_se_text(path: &Path) -> Result<(SingleEndTrack, f64)> {
    use std::io::BufReader;
    let f = macs_io::open_maybe_gzip(path)?;
    let mut r = BufReader::new(f);
    let mut line = Vec::new();
    let mut b = SingleEndTrackBuilder::new();
    // the first ten positive-length tags, exactly as `get_tag_size` samples them
    let (mut sum, mut n) = (0f64, 0usize);
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
        // F192: `BEDParser.tlen_parse_line` is `atoi(col3) - atoi(col2)` -- a *signed*
        // difference between the two columns exactly as written, with no strand
        // handling at all -- and `GenericParser.d` only accumulates it
        // `if this_taglength > 0`. So a record with `col3 < col2` contributes
        // nothing. Taking `|end - pos|` instead counted those records with a positive
        // magnitude: on `3 x (col3-col2 = 100)` plus one reversed `col3 - col2 = -800`,
        // upstream reported `tag size = 100` and this reported `275`.
        //
        // `rec.pos` is column 3 for a minus-strand record and column 2 otherwise, so
        // `col3 - col2` is `pos - end` on the minus strand and `end - pos` on the plus.
        if n < 10 {
            let tlen = match rec.strand {
                Strand::Minus => rec.pos - rec.end,
                _ => rec.end - rec.pos,
            };
            if tlen > 0 {
                sum += tlen as f64;
                n += 1;
            }
        }
    }
    b.finalize();
    // `cython.cast(cython.int, s/n)`: truncated, and `-1` when nothing parsed
    let mean = if n != 0 {
        (sum / n as f64).trunc()
    } else {
        -1.0
    };
    Ok((b.build(), mean))
}

/// What one paired-end input contributes to `options.tsize`.
///
/// `mean_before_filter` is the count-weighted mean
/// ([`FragmentTrack::average_template_length`]); `mean_row` is upstream's
/// `self.d` -- the **unweighted** mean of the row lengths, which is a different
/// number as soon as a row carries a multiplicity (F166).
struct PeLoad {
    track: FragmentTrack,
    /// Total fragments *before* `--barcodes` filtering.
    ///
    /// Upstream reports `t0 = treat.total` in the xls *before* applying
    /// `treat.subset(barcodes_set)` (`callpeak_cmd.py:80-92`), so with `--barcodes`
    /// the header shows the unfiltered count (67421) while peak calling uses the
    /// subset (31359). Filtering at load time never sees the unfiltered total, so
    /// it is counted separately here.
    pre_subset_total: u64,
    /// Count-weighted mean (`FragmentTrack::average_template_length`).
    ///
    /// Retained because it is the value `# fragment size is determined as %d bps`
    /// truncates, and a caller that wants `options.tsize` needs it; the scale
    /// factors themselves use `mean_row` (F166).
    #[allow(dead_code)]
    mean_before_filter: f64,
    mean_row: f64,
    rows: u64,
}

/// Read a paired-end file, returning the track and the two pre-filter means.
///
/// The BEDPE parser is deliberately looser than upstream's `FragmentTrackParser`
/// here: the fixtures are whitespace-separated, and upstream's own
/// `pe_parse_line` splits on tab and then falls back to space. Both are accepted;
/// the *coordinates* are still validated (`0 <= start <= end`).
/// Read the `--barcodes` allow-list.
///
/// `callpeak_cmd.py:86-92`: when `-f FRAG` and `--barcodes` are both given, upstream
/// reads the file line by line, strips trailing whitespace, and keeps only fragments
/// whose barcode is in the set (`treat.subset(barcodes_set)`). The set is matched as
/// raw bytes; there is no normalization beyond the rstrip.
fn read_barcodes(path: &Path) -> Result<std::collections::HashSet<Vec<u8>>> {
    use std::io::BufRead;
    let r = macs_io::open_maybe_gzip(path)?;
    let mut out = std::collections::HashSet::new();
    for line in r.split(b'\n') {
        let line = line?;
        // `l.rstrip()`: strip trailing ASCII whitespace (newline, CR, spaces, tabs).
        let mut end = line.len();
        while end > 0 && matches!(line[end - 1], b'\n' | b'\r' | b' ' | b'\t') {
            end -= 1;
        }
        if end > 0 {
            out.insert(line[..end].to_vec());
        }
    }
    Ok(out)
}

fn load_pe_text(
    path: &Path,
    frag: bool,
    barcodes: Option<&std::collections::HashSet<Vec<u8>>>,
) -> Result<PeLoad> {
    let f = macs_io::open_maybe_gzip(path)?;
    let mut r = std::io::BufReader::new(f);
    let mut line = Vec::new();
    let mut b = if frag {
        FragTrackBuilder::with_barcodes()
    } else {
        FragTrackBuilder::new()
    };
    let (mut sum_len, mut n_rows) = (0f64, 0u64);
    let mut pre_subset_total = 0u64;
    // Unweighted row mean over *all* rows, before `--barcodes` filtering. Upstream's
    // `FragParser.d` is set during parsing (`m / i` over every row) and `subset()`
    // preserves it, so with `--barcodes` the reported fragment size (199) is the
    // full-set mean (199.77), not the subset mean (198.09 -> 198).
    let (mut sum_len_all, mut n_rows_all) = (0f64, 0u64);
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if frag {
            let Some(rec) = macs_io::parse_frag_line(&line)? else {
                continue;
            };
            if rec.chrom.is_empty() || rec.left < 0 || rec.right < rec.left {
                continue;
            }
            // Counted before the `--barcodes` filter: pre-subset total (`t0`) and the
            // pre-subset row mean (`d`) upstream reports.
            sum_len_all += (rec.right - rec.left) as f64;
            n_rows_all += 1;
            if barcodes.is_some() {
                pre_subset_total += rec.count.unwrap_or(1) as u64;
            }
            // `--barcodes` subsets the fragments (`callpeak_cmd.py:86-92`). It has no
            // effect when the format is not FRAG, and when it is absent every row is
            // kept. A row with no barcode column can never match a non-empty set.
            if let Some(allow) = barcodes {
                match &rec.barcode {
                    Some(bc) if allow.contains(bc) => {}
                    _ => continue,
                }
            }
            sum_len += (rec.right - rec.left) as f64;
            n_rows += 1;
            b.push_with_count(
                &rec.chrom,
                rec.left as u64,
                rec.right as u64,
                rec.count.unwrap_or(1),
            );
        } else {
            let Some(rec) = macs_io::parse_bedpe_line(&line)? else {
                continue;
            };
            if rec.chrom.is_empty() || rec.left < 0 || rec.right < rec.left {
                continue;
            }
            sum_len += (rec.right - rec.left) as f64;
            n_rows += 1;
            b.push(&rec.chrom, rec.left as u64, rec.right as u64);
        }
    }
    b.finalize();
    let t = b.build();
    // F92: `options.tsize` is the mean template length **as read**, before
    // duplicate filtering, so it has to be captured here.
    //
    // F166: upstream's parser accumulates `m += right_pos - left_pos` and sets
    // `self.d = m / i` over **rows**, with no weight, for every paired-end format.
    // `FragmentTrack::average_template_length` is the *count-weighted* mean
    // (`length / total`), which is a different number as soon as a row carries a
    // multiplicity -- 98.307 against 98.4048 on
    // `frag_basic/barcode_fragments`, and it lands in `float(self.d)` for the
    // wide-window control scale factors. So the unweighted row mean is
    // accumulated here, while the parser is still at the rows.
    let mean_before_filter = t.average_template_length();
    // With `--barcodes`, upstream's `d` is the full-set row mean, not the subset's.
    let (sum_len, n_rows) = if barcodes.is_some() {
        (sum_len_all, n_rows_all)
    } else {
        (sum_len, n_rows)
    };
    let mean_row = if n_rows > 0 {
        sum_len / n_rows as f64
    } else {
        0.0
    };
    let pre_subset_total = if barcodes.is_some() {
        pre_subset_total
    } else {
        t.total()
    };
    Ok(PeLoad {
        track: t,
        pre_subset_total,
        mean_before_filter,
        mean_row,
        rows: n_rows,
    })
}

/// `--keep-dup`: the integer form, or an error for the two symbolic values.
///
/// `all` means no filtering and is handled by returning `i64::MAX` as the cap;
/// `auto` needs MACS's binomial model and is rejected explicitly rather than
/// approximated.
/// `--format FRAG` forces duplicate filtering off.
///
/// `OptValidator.validate_options` (`MACS3/Utilities/OptValidator.py:110-112`):
///
/// ```python
/// if options.format == 'FRAG' and options.keepduplicates != "all":
///     info("Since the FRAG format specified, all duplicates in reads will be kept.")
///     options.keepduplicates = "all"
/// ```
///
/// FRAG rows already carry a multiplicity, so deduplicating them would discard
/// the weights that make the format meaningful -- and, because the xls count
/// block is gated on `keepduplicates != "all"`, it would also add four
/// "after filtering" lines upstream never writes for FRAG.
fn keep_dup_for(format: &str, o: &Options) -> Result<DupPolicy> {
    if format == "FRAG" {
        return Ok(DupPolicy::Fixed(i64::MAX));
    }
    match o.get("keepduplicates") {
        Some("auto") => Ok(DupPolicy::Auto),
        _ => keep_dup(o).map(DupPolicy::Fixed),
    }
}

/// `--keep-dup`, before it is resolved against the loaded data.
#[derive(Debug, Clone, Copy, PartialEq)]
enum DupPolicy {
    Fixed(i64),
    /// Resolve with [`cal_max_dup_tags`].
    Auto,
}

/// `cal_max_dup_tags(gsize, tags_number, p=1e-5)` (`callpeak_cmd.py:331`):
///
/// ```python
/// return binomial_cdf_inv(1-p, tags_number, 1.0/genome_size)
/// ```
///
/// i.e. the smallest `x` whose binomial CDF reaches `1 - p` for
/// `Binomial(tags_number, 1/gsize)`. It needs the **retained** treatment count, so
/// it can only be resolved after loading; the same number is then used for the
/// control (`callpeak_cmd.py:141` filters the control with
/// `treatment_max_dup_tags`, not a control-specific one).
fn cal_max_dup_tags(gsize: f64, tags: u64) -> Result<i64> {
    if gsize <= 0.0 {
        return Err(MacsError::InvalidParameter(
            "--keep-dup auto needs a positive effective genome size".into(),
        ));
    }
    macs_stats::binomial::binomial_cdf_inv(1.0 - 1e-5, tags as i64, 1.0 / gsize)
        .map_err(|e| MacsError::InvalidParameter(format!("--keep-dup auto: {e}")))
}

fn keep_dup(o: &Options) -> Result<i64> {
    match o.get("keepduplicates") {
        None => Ok(1),
        Some("all") => Ok(i64::MAX),
        Some("auto") => Err(MacsError::InvalidParameter(
            "--keep-dup auto needs MACS's own binomial model and is not reproduced; \
             use an integer or `all`"
                .into(),
        )),
        Some(v) => v.parse::<i64>().map_err(|_| {
            MacsError::InvalidParameter(format!(
                "--keep-dup must be an integer or `all`, got `{v}`"
            ))
        }),
    }
}

/// Pool several input files into one track, as `-t A B C` and `-c A B` do.
///
/// Duplicate filtering is applied by the caller on the pooled track, not per
/// file, so `-t A B` where A and B share a position keeps one copy -- which is
/// what upstream does, since it appends to a single `FWTrack`.
/// F27 — `check_names` (`callpeak_cmd.py:32`): refuse to run when the treatment and
/// the control share no chromosome name at all.
///
/// Upstream prints the two diagnostic lines and then dies with status 1 (its
/// third `error_stream` call raises `TypeError` on `bytes`, so `sys.exit()` is
/// never reached — see `crates/macs-cli/src/error.rs`). The corpus records that
/// as exit 1 with no output file, which is what this reproduces. The check has to
/// happen while the files are being read, before any writer is opened.
fn ensure_shared_chroms<'a, I, J>(treat: I, ctrl: J) -> Result<()>
where
    I: IntoIterator<Item = &'a [u8]>,
    J: IntoIterator<Item = &'a [u8]>,
{
    let cnames: BTreeSet<Vec<u8>> = ctrl.into_iter().map(|n| n.to_vec()).collect();
    if treat.into_iter().any(|n| cnames.contains(n)) {
        Ok(())
    } else {
        Err(MacsError::Rejected(
            crate::error::NO_COMMON_CHROMOSOMES.to_string(),
        ))
    }
}

/// Stand-in for a chromosome's q-score track when `--broad` is off and the
/// tracks were therefore not materialised. [`ChromCall::qtrack`] is read only
/// inside the `broad` branch, so an empty track is never indexed.
static EMPTY_TRACK: std::sync::LazyLock<macs_rle::SignalTrack<f32>> =
    std::sync::LazyLock::new(|| macs_rle::SignalTrack::empty(macs_core::ChromId(0), 0, 0));

/// The rayon pool the per-chromosome work runs on.
///
/// callpeak has **no** `--threads` flag -- upstream's argparse does not define
/// one, so adding a CLI flag would break the "identical accept/reject behaviour"
/// criterion (a real `macs3` rejects `--threads` with exit 2). Thread count is
/// therefore an environment override, which the determinism harness sets.
///
/// `MACS3_RS_THREADS` unset means "rayon's default", i.e. one thread per core.
fn configure_pool() -> rayon::ThreadPool {
    match std::env::var("MACS3_RS_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
    {
        Some(n) if n >= 1 => rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .expect("a positive thread count is always a valid pool"),
        _ => rayon::ThreadPoolBuilder::new()
            .build()
            .expect("rayon's default pool is always constructible"),
    }
}

fn pool_se(paths: &[String], format: &str, infer_tsize: bool) -> Result<(SingleEndTrack, f64)> {
    let first = paths
        .first()
        .ok_or_else(|| MacsError::InvalidParameter("-i/--ifile is required".into()))?;
    let first_format = super::input::single_end_format(Path::new(first), format)?;
    // `load_tag_files_options` asks only the first treatment parser for tsize,
    // and skips that call when `--tsize` was supplied. Bowtie's parser can fail
    // while estimating tsize on short files, but control parsers never estimate it.
    if infer_tsize && first_format == "BOWTIE" {
        super::input::validate_bowtie_tsize(Path::new(first))?;
    }
    let mut b = SingleEndTrackBuilder::new();
    let mut first_mean = None;
    for p in paths {
        let path = Path::new(p);
        // OptValidator assigns `guess_parser` for AUTO; `load_tag_files_options`
        // invokes it for each file, so mixed-format pools must resolve each path.
        let detected = super::input::single_end_format(path, format)?;
        let (t, mean) = match detected.as_str() {
            "BED" => load_se_text(path)?,
            _ => super::input::load_single_end(path, &detected)?,
        };
        if paths.len() == 1 {
            return Ok((t, mean));
        }
        first_mean.get_or_insert(mean);
        let pos = t.positions();
        for c in pos.genome().ids_file_order() {
            let name = pos.genome().name(c).to_vec();
            for &p in pos.strand(c, macs_core::Strand::Plus) {
                b.push(&name, p, macs_core::Strand::Plus);
            }
            for &p in pos.strand(c, macs_core::Strand::Minus) {
                b.push(&name, p, macs_core::Strand::Minus);
            }
        }
    }
    b.finalize();
    Ok((b.build(), first_mean.unwrap_or(0.0)))
}

/// Select the temporary directory used by callpeak's intermediate files.
/// The flag matrix carries the upstream machine's `/scratch/.../tmp` default,
/// so an omitted flag must use the current platform's standard temp directory.
fn callpeak_tempdir(o: &Options) -> PathBuf {
    let explicitly_set = o
        .raw_argv
        .iter()
        .any(|arg| arg == "--tempdir" || arg.starts_with("--tempdir="));
    if explicitly_set {
        PathBuf::from(o.get("tempdir").unwrap_or_default())
    } else {
        std::env::temp_dir()
    }
}

fn pool_pe(
    paths: &[String],
    frag: bool,
    barcodes: Option<&std::collections::HashSet<Vec<u8>>>,
) -> Result<(FragmentTrack, f64, f64, u64)> {
    let mut b = if frag {
        FragTrackBuilder::with_barcodes()
    } else {
        FragTrackBuilder::new()
    };
    // `(pooled_count_weighted_mean, unweighted_row_mean)`
    //
    // F166: the two differ as soon as a row carries a multiplicity. The first is
    // what `FragmentTrack::average_template_length` reports and is only used for
    // `options.tsize`'s reported value; the second is upstream's `self.d`, which
    // `PeakDetect.__init__` feeds to `control_sum` and, unrounded, to
    // `float(self.d)/lregion*ratio`.
    // `--barcodes` is read once and shared across inputs. Upstream reads it inside
    // `load_frag_files_options` per file, but the file does not change between inputs
    // so the set is identical; reading it once avoids re-parsing.

    let (mut sum, mut n, mut sum_row, mut n_row) = (0f64, 0f64, 0f64, 0f64);
    let mut pre_subset_total = 0u64;
    for p in paths {
        let PeLoad {
            track: t,
            pre_subset_total: pre,
            mean_before_filter: _,
            mean_row,
            rows,
        } = load_pe_text(Path::new(p), frag, barcodes)?;
        pre_subset_total += pre;
        sum += t.length() as f64;
        n += t.total() as f64;
        sum_row += mean_row * rows as f64;
        n_row += rows as f64;
        for c in t.chroms() {
            let name = t.genome().name(c).to_vec();
            if t.has_counts() {
                // F151: pooling must carry the multiplicity. Dropping it here made
                // every `--format FRAG` run unit-weighted: the depth at a summit
                // came out 947 where upstream reports 3783, i.e. divided by the
                // mean count.
                for (f, &w) in t.frags(c).iter().zip(t.counts(c).iter()) {
                    b.push_with_count(&name, f.start, f.end, u32::from(w));
                }
            } else {
                for f in t.frags(c) {
                    b.push(&name, f.start, f.end);
                }
            }
        }
    }
    b.finalize();
    let mean = if n > 0.0 { sum / n } else { 0.0 };
    let mean_row = if n_row > 0.0 { sum_row / n_row } else { 0.0 };
    Ok((b.build(), mean, mean_row, pre_subset_total))
}

/// The chromosome name list the paired-end pipeline needs, from both tracks.
fn names_of(treat: &FragmentTrack, ctrl: Option<&FragmentTrack>) -> Vec<ChromId> {
    let mut seen: BTreeMap<Vec<u8>, ChromId> = BTreeMap::new();
    for c in treat.chroms() {
        seen.insert(treat.genome().name(c).to_vec(), c);
    }
    if let Some(c) = ctrl {
        for cid in c.chroms() {
            seen.entry(c.genome().name(cid).to_vec()).or_insert(cid);
        }
    }
    seen.into_values().collect()
}

/// Run `callpeak` with the already-validated options.
pub fn run(o: &Options) -> Result<()> {
    let treat_paths = o.get_all("tfile").to_vec();
    if treat_paths.is_empty() {
        return Err(MacsError::InvalidParameter("-t/--tfile is required".into()));
    }
    let ctrl_paths = o.get_all("cfile").to_vec();
    let format = o.get("format").unwrap_or("AUTO").to_uppercase();
    let is_pe = PE_FORMATS.contains(&format.as_str());
    let is_bam = BAM_FORMATS.contains(&format.as_str());
    let dup_policy = keep_dup_for(&format, o)?;
    let name = o.get("name").unwrap_or("MACS").to_string();
    let outdir = PathBuf::from(o.get("outdir").unwrap_or("."));
    let gsize = match o.get("gsize") {
        Some(spec) => {
            macs_core::genomesize::resolve_gsize(spec).map_err(MacsError::InvalidParameter)?
        }
        None => return Err(MacsError::InvalidParameter("-g/--gsize is required".into())),
    };
    let cfg_common = PeLike {
        gsize,
        slocal: o.int("smalllocal").unwrap_or(1000),
        llocal: o.int("largelocal").unwrap_or(10000),
        qvalue: o.float("qvalue").unwrap_or(0.05),
        call_summits: o.flag("call_summits"),
        broad: o.flag("broad"),
        broad_cutoff: o.float("broadcutoff").unwrap_or(0.1),
        nolambda: o.flag("nolambda"),
        // F189: `log_pvalue = -log10(pvalue)`; `call_peaks` dispatches on it
        // being set, and `-p` takes precedence over `-q`.
        p_cutoff: o
            .float("pvalue")
            .filter(|v| *v > 0.0)
            .map(|v| -(v as f32).log10()),
    };
    let extsize = o.int("extsize").unwrap_or(0);

    // `# alternative fragment length(s) may be ...` is emitted into the xls header only
    // when a `PeakModel` was actually fitted (`PeakIO.py` writes it from
    // `opt.alternative_d`, which only the model path sets), so it is carried out of the
    // load branch rather than recomputed at header time.
    let mut alt_d_line: Option<String> = None;

    // `--mfold` is `nargs=2`; upstream reads `lmfold = mfold[0]`, `umfold = mfold[1]`
    // (`OptValidator.py:146-147`) and defaults to `[5, 50]`.
    let mfold_pair = {
        let v = o.get_all("mfold");
        let lo = v.first().and_then(|x| x.parse::<f64>().ok()).unwrap_or(5.0);
        let hi = v.get(1).and_then(|x| x.parse::<f64>().ok()).unwrap_or(50.0);
        (lo, hi)
    };

    // ---- load -------------------------------------------------------------
    // `tsize`/`t0`/`t1`/`c0`/`c1` are retained for the xls header's per-track
    // count block; each load branch fills them.
    #[allow(unused_assignments)]
    let mut tsize: i64 = 0;
    #[allow(unused_assignments)]
    let mut t0 = 0u64;
    // The per-chromosome read counts **before** duplicate filtering, for
    // `reads_pre_filter`. Captured here rather than inside `record_reads` because by the
    // time the pipeline calls that, `filter_dup` has already rewritten the track --
    // reporting post-filter counts under a `pre_filter` name would manufacture a
    // divergence on every paired-end fixture, which is exactly what a stage dump must
    // not do.
    let mut sdump_pre: Option<String> = None;
    #[allow(unused_assignments)]
    let mut t1 = 0u64;
    let mut c0 = 0u64;
    let mut c1 = 0u64;
    // F166: `float(self.d)` as the paired-end text loader reports it, and whether
    // any loader set it at all (the BAM path has no row-length mean).
    let mut pe_d_exact = 0.0f64;
    let mut pe_d_exact_set = false;
    #[allow(unused_mut)]
    let mut max_dup_line: Option<i64> = match dup_policy {
        DupPolicy::Fixed(i64::MAX) | DupPolicy::Auto => None,
        DupPolicy::Fixed(v) => Some(v),
    };
    // Harness-only stage capture (`MACS3_RS_DUMP_STAGES`). An env var, not a flag:
    // a flag upstream does not have would be an accept/reject divergence. Inert unless
    // the variable is set.
    let rss_trace = std::env::var("MACS3_RS_RSS_TRACE").is_ok();
    if rss_trace {
        eprintln!("rss[entry] {} kB", rss_kb());
    }
    let mut sdump = StageDump::from_env();
    let sdump_on = sdump.enabled();

    // `SeSignalSetup` when the SE path runs, so `scaling` can be recorded; `None` on the
    // paired-end paths, which compute their own ratio.
    let mut scaling_setup: Option<macs_peaks::callpeak::SeSignalSetup> = None;
    let mut lambda_ladder: Option<String> = None;
    let mut se_source: Option<SeSignalSource> = None;
    let mut pe_source: Option<PeSignalSource> = None;
    let mut stream_bdg: Option<(tempfile::NamedTempFile, tempfile::NamedTempFile)> = None;
    let mut se_signal_spool: Option<tempfile::NamedTempFile> = None;
    let tempdir = callpeak_tempdir(o);
    let (signals, d, paired_boundaries, lambda_bg, coord_shift) = if is_pe {
        let pe_cfg = PeConfig {
            tsize: 0.0,       // replaced below from the loaded track
            tsize_exact: 0.0, // ditto
            gsize: cfg_common.gsize,
            slocal: cfg_common.slocal,
            llocal: cfg_common.llocal,
            qvalue: cfg_common.qvalue,
            call_summits: cfg_common.call_summits,
            broad: cfg_common.broad,
            broad_cutoff: cfg_common.broad_cutoff,
            nolambda: cfg_common.nolambda,
            scaleto_large: o.get("scaleto").unwrap_or("small") == "large",
        };
        if is_bam {
            // BAMPE: one file at a time, pooled like the text readers
            let mut b = FragTrackBuilder::new();
            let mut sum = 0f64;
            let mut n = 0f64;
            // F161: `--keep-dup auto` resolves against the treatment's retained
            // count, so BAM loads need it up front; `--keep-dup all` is a plain
            // pass-through.
            let mut dup_bam = match dup_policy {
                DupPolicy::Fixed(v) => v,
                DupPolicy::Auto => i64::MAX,
            };
            for p in &treat_paths {
                let (frags, summary) = macs_io::bam::bampe_fragments(Path::new(p))?;
                sum += summary.d * summary.n as f64;
                n += summary.n as f64;
                for fr in &frags {
                    b.push(&fr.chrom, u64::from(fr.start), u64::from(fr.start + fr.len));
                }
            }
            b.finalize();
            let mut treat = b.build();
            let mean = if n > 0.0 { sum / n } else { 0.0 };
            t0 = treat.total();
            let pre_treatment = sdump_on.then(|| stagedump::pre_json_frag(&treat, "treatment"));
            if dup_policy == DupPolicy::Auto {
                dup_bam = cal_max_dup_tags(gsize, t0)?;
                max_dup_line = if dup_bam == i64::MAX {
                    None
                } else {
                    Some(dup_bam)
                };
                if dup_bam != i64::MAX {
                    filter_frag_dup(&mut treat, dup_bam)?;
                }
            } else if dup_bam != i64::MAX {
                filter_frag_dup(&mut treat, dup_bam)?;
            }
            t1 = treat.total();
            let ctrl = if ctrl_paths.is_empty() {
                None
            } else {
                let c = {
                    let mut cb = FragTrackBuilder::new();
                    for path in &ctrl_paths {
                        let (frags, _) = macs_io::bam::bampe_fragments(Path::new(path))?;
                        for fr in &frags {
                            cb.push(
                                &fr.chrom,
                                u64::from(fr.start),
                                u64::from(fr.start) + u64::from(fr.len),
                            );
                        }
                    }
                    cb.finalize();
                    let mut c = cb.build();
                    c0 = c.total();
                    if sdump_on {
                        let pre = pre_treatment
                            .as_deref()
                            .expect("the treatment pre-filter snapshot was recorded");
                        sdump_pre = Some(stagedump::add_pre_control_frag(pre, Some(&c)));
                    }
                    // The control gets the same duplicate filter as the treatment
                    // (`callpeak_cmd.py:130-137`). Omitting this emptied the control
                    // track entirely and reported `Redundant rate in control: 1.00`.
                    if dup_bam != i64::MAX {
                        filter_frag_dup(&mut c, dup_bam)?;
                    }
                    c1 = c.total();
                    c
                };
                Some(c)
            };
            if sdump_on && ctrl.is_none() {
                let pre = pre_treatment
                    .as_deref()
                    .expect("the treatment pre-filter snapshot was recorded");
                sdump_pre = Some(stagedump::add_pre_control_frag(pre, None));
            }
            tsize = (o.int("tsize").map_or(mean, |v| v as f64)) as i64;
            let mut cfg = pe_cfg;
            // F148: with `--nomodel` in PE mode, `options.d = options.tsize`
            // (`callpeak_cmd.py:165`) and `--extsize` is **ignored** -- there is
            // no `info("#2 Use %d as fragment length")` line for PE, and the
            // measured mean is used even when `--extsize 200` is passed.
            // `--tsize` (`-s`) is the only override, and it is applied earlier
            // (`if not options.tsize: options.tsize = ttsize`).
            cfg.tsize = o.int("tsize").map_or(mean, |v| v as f64);
            // F238: upstream's `tp.d` is **f32** -- `Parser.py:1496` is
            // `self.d = cython.cast(cython.float, m) / i`, and both operands are
            // C floats -- so `options.tsize` is a float32 that `float(self.d)`
            // then widens without loss. Keeping our mean in f64 shifts the wide
            // window factors by one ulp, which is invisible in a coordinate diff
            // and visible in the fifth decimal of the control lambda:
            //
            //     f64 mean 143.555                    -> factor 0.007177755236625671
            //                                            2194 * factor = 0x417bf7ca (15.74800)
            //     f32 mean 143.55499267578125          -> factor 0.007177754770964384
            //                                            2194 * factor = 0x417bf7c9 (15.74799) <- upstream
            //
            // `depth = 2194` is exact either way (a counted control's depth is an
            // integer sum of counts), so the factor is the only thing that moved.
            cfg.tsize_exact = if pe_d_exact_set {
                pe_d_exact as f32 as f64
            } else {
                cfg.tsize as f32 as f64
            };
            if sdump_on {
                let r = macs_peaks::callpeak::run_callpeak_pe(&treat, ctrl.as_ref(), &cfg);
                record_pe_scaling_stage(&mut sdump, &treat, ctrl.as_ref(), &cfg);
                lambda_ladder = Some(pe_lambda_ladder_json(
                    &treat,
                    ctrl.as_ref(),
                    &cfg,
                    r.lambda_bg,
                ));
                record_lambda_merged_stage(
                    &mut sdump,
                    &r.signals,
                    r.coord_shift,
                    spmr_denominator(o, t1, c1),
                );
                stagedump::record_reads(
                    &mut sdump,
                    sdump_pre.as_deref(),
                    ctrl.is_some(),
                    tsize as f64,
                    t1,
                    c1,
                    cfg_common.slocal,
                    cfg_common.llocal,
                );
                stagedump::record_tracks(
                    &mut sdump,
                    &signal_tracks(&r.signals),
                    &[],
                    r.coord_shift,
                );
                (
                    r.signals,
                    r.d,
                    r.paired_boundaries,
                    r.lambda_bg,
                    r.coord_shift,
                )
            } else {
                let (source, d, boundaries, lambda, shift) = PeSignalSource::new(treat, ctrl, cfg);
                pe_source = Some(source);
                (Vec::new(), d, boundaries, lambda, shift)
            }
        } else {
            let frag = format == "FRAG";
            // `--barcodes` subsets `-f FRAG` fragments to an allow-list
            // (`callpeak_cmd.py:86-92`). It is read once here and shared by the
            // treatment and control pools; it has no effect for other formats.
            let barcode_set = if frag {
                match o.get("barcodefile") {
                    Some(bf) if !bf.is_empty() => Some(read_barcodes(Path::new(bf))?),
                    _ => None,
                }
            } else {
                None
            };
            // F258: the count-weighted mean is no longer read here; upstream's
            // `self.d` is the unweighted row mean (`mean_row`).
            let (mut treat, _count_weighted_mean, mean_row, pre_subset) =
                pool_pe(&treat_paths, frag, barcode_set.as_ref())?;
            let pre_treatment = sdump_on.then(|| stagedump::pre_json_frag(&treat, "treatment"));
            // With `--barcodes`, `t0` is the pooled total *before* subsetting
            // (`callpeak_cmd.py:80`), not the filtered track's total. Without a
            // barcode filter the two are the same.
            t0 = if barcode_set.is_some() {
                pre_subset
            } else {
                treat.total()
            };
            let dup = match dup_policy {
                DupPolicy::Fixed(v) => v,
                DupPolicy::Auto => cal_max_dup_tags(gsize, t0)?,
            };
            // F185: the "after filtering" header block is gated on
            // `options.keepduplicates != "all"` (`callpeak_cmd.py:96`), not on the
            // value being a literal. `--keep-dup auto` resolves to a number and
            // therefore *does* print all four lines, with the resolved number in
            // the "maximum duplicate" line. Suppressing them for `auto` left every
            // `keepdup_auto` fixture's `*.xls` three lines short.
            max_dup_line = (dup != i64::MAX).then_some(dup);
            if dup != i64::MAX {
                filter_frag_dup(&mut treat, dup)?;
            }
            t1 = treat.total();
            let ctrl = if ctrl_paths.is_empty() {
                None
            } else {
                let (mut c, _, _, _) = pool_pe(&ctrl_paths, frag, barcode_set.as_ref())?;
                ensure_shared_chroms(
                    treat.chroms().iter().map(|&cid| treat.genome().name(cid)),
                    c.chroms().iter().map(|&cid| c.genome().name(cid)),
                )?;
                c0 = c.total();
                if sdump_on {
                    let pre = pre_treatment
                        .as_deref()
                        .expect("the treatment pre-filter snapshot was recorded");
                    sdump_pre = Some(stagedump::add_pre_control_frag(pre, Some(&c)));
                }
                if dup != i64::MAX {
                    filter_frag_dup(&mut c, dup)?;
                }
                c1 = c.total();
                Some(c)
            };
            if sdump_on && ctrl.is_none() {
                let pre = pre_treatment
                    .as_deref()
                    .expect("the treatment pre-filter snapshot was recorded");
                sdump_pre = Some(stagedump::add_pre_control_frag(pre, None));
            }
            // F258: the reported tag/fragment size is upstream's `self.d`, which
            // `Parser.py:1496` computes **unweighted**:
            //
            //     self.d = cython.cast(cython.float, m) / i   # m += right - left
            //
            // `m` accumulates fragment lengths once per *row*, with no multiplicity,
            // and `callpeak_cmd.py:382/415` then prints `options.tsize`. We were
            // passing the **count-weighted** mean here, which is a different number
            // as soon as a row carries a count -- exactly the situation every
            // `--format FRAG` fixture with counts is in. Measured on
            // `sweep/gmini_mfrag_d400_w180_ctrl_357`:
            //
            //     count-weighted mean -> 143.xx -> "# fragment size is determined as 143 bps"
            //     unweighted row mean  -> 144.155 -> 144, which is what upstream prints
            //
            // so this was wrong in the xls header (`# d = ...` and the
            // `# fragment size is determined as ... bps` line) for every counted
            // paired-end fixture. `mean_row` is the same quantity F166 already routes
            // to `cfg.tsize_exact`, so the header and the scale factors now agree on
            // one `d` instead of two.
            tsize = o.int("tsize").map_or(mean_row, |v| v as f64) as i64;
            let mut cfg = pe_cfg;
            // F148: `--extsize` is ignored in PE mode; `-s/--tsize` is the only
            // override of the measured mean.
            cfg.tsize = tsize as f64;
            // F166: `float(self.d)` in the wide-window control factors sees the
            // parser's **unweighted** row mean, unrounded; `self.d` itself is
            // truncated where it reaches `cython.int`.
            pe_d_exact = mean_row;
            pe_d_exact_set = true;
            // F238: upstream's `tp.d` is **f32** -- `Parser.py:1496` is
            // `self.d = cython.cast(cython.float, m) / i`, and both operands are
            // C floats -- so `options.tsize` is a float32 that `float(self.d)`
            // then widens without loss. Keeping our mean in f64 shifts the wide
            // window factors by one ulp, which is invisible in a coordinate diff
            // and visible in the fifth decimal of the control lambda:
            //
            //     f64 mean 143.555                    -> factor 0.007177755236625671
            //                                            2194 * factor = 0x417bf7ca (15.74800)
            //     f32 mean 143.55499267578125          -> factor 0.007177754770964384
            //                                            2194 * factor = 0x417bf7c9 (15.74799) <- upstream
            //
            // `depth = 2194` is exact either way (a counted control's depth is an
            // integer sum of counts), so the factor is the only thing that moved.
            cfg.tsize_exact = if pe_d_exact_set {
                pe_d_exact as f32 as f64
            } else {
                cfg.tsize as f32 as f64
            };
            if sdump_on {
                let r = macs_peaks::callpeak::run_callpeak_pe(&treat, ctrl.as_ref(), &cfg);
                record_pe_scaling_stage(&mut sdump, &treat, ctrl.as_ref(), &cfg);
                lambda_ladder = Some(pe_lambda_ladder_json(
                    &treat,
                    ctrl.as_ref(),
                    &cfg,
                    r.lambda_bg,
                ));
                record_lambda_merged_stage(
                    &mut sdump,
                    &r.signals,
                    r.coord_shift,
                    spmr_denominator(o, t1, c1),
                );
                stagedump::record_reads(
                    &mut sdump,
                    sdump_pre.as_deref(),
                    ctrl.is_some(),
                    tsize as f64,
                    t1,
                    c1,
                    cfg_common.slocal,
                    cfg_common.llocal,
                );
                stagedump::record_tracks(
                    &mut sdump,
                    &signal_tracks(&r.signals),
                    &[],
                    r.coord_shift,
                );
                (
                    r.signals,
                    r.d,
                    r.paired_boundaries,
                    r.lambda_bg,
                    r.coord_shift,
                )
            } else {
                let (source, d, boundaries, lambda, shift) = PeSignalSource::new(treat, ctrl, cfg);
                pe_source = Some(source);
                (Vec::new(), d, boundaries, lambda, shift)
            }
        }
    } else {
        let infer_tsize = o.int("tsize").unwrap_or(0) == 0;
        let (mut treat, mean_treat) = pool_se(&treat_paths, &format, infer_tsize)?;
        if sdump_on {
            sdump_pre = Some(stagedump::pre_json_se(&treat));
        }
        t0 = treat.total();
        // F161: `--keep-dup auto` can only be resolved once the treatment's
        // retained count is known.
        let dup = match dup_policy {
            DupPolicy::Fixed(v) => v,
            DupPolicy::Auto => cal_max_dup_tags(gsize, t0)?,
        };
        if dup != i64::MAX {
            treat.filter_dup(dup)?;
        }
        t1 = treat.total();
        // F185: see the paired-end site -- `auto` prints the filtering block too.
        max_dup_line = (dup != i64::MAX).then_some(dup);
        // F153: the header prints the **measured** mean tag length
        // (`callpeak_cmd.py:80`, `%d` so truncated), not `--extsize`. `d` itself is
        // the extsize, which is why a 190 bp fixture with `--extsize 200` reports
        // 190 in the header yet extends reads by 200.
        tsize = o
            .int("tsize")
            .filter(|v| *v != 0)
            .unwrap_or(mean_treat as i64);
        let ctrl = if ctrl_paths.is_empty() {
            None
        } else {
            let (mut c, _mean_c) = pool_se(&ctrl_paths, &format, false)?;
            ensure_shared_chroms(
                treat
                    .positions()
                    .genome()
                    .ids_file_order()
                    .into_iter()
                    .map(|cid| treat.genome().name(cid)),
                c.positions()
                    .genome()
                    .ids_file_order()
                    .into_iter()
                    .map(|cid| c.genome().name(cid)),
            )?;
            c0 = c.total();
            if sdump_on {
                let pre = sdump_pre
                    .take()
                    .expect("the treatment pre-filter snapshot was recorded");
                sdump_pre = Some(stagedump::add_pre_control_se(&pre, Some(&c)));
            }
            if dup != i64::MAX {
                c.filter_dup(dup)?;
            }
            c1 = c.total();
            Some(c)
        };
        if sdump_on && ctrl.is_none() {
            let pre = sdump_pre
                .take()
                .expect("the treatment pre-filter snapshot was recorded");
            sdump_pre = Some(stagedump::add_pre_control_se(&pre, None));
        }
        if extsize <= 0 && !o.flag("nomodel") {
            return Err(MacsError::InvalidParameter(
                "--extsize is required when the model cannot be fitted".into(),
            ));
        }
        // #2 Build Peak Model (`callpeak_cmd.py:159-210`).
        //
        // This is the *default* path -- `callpeak` without `--nomodel` -- and it used to
        // be unimplemented: `d` was silently taken from `--extsize` (or, when that was
        // absent, from the tag size), so the run completed with exit 0 and wrote three
        // plausible-looking files built on a guessed fragment length. On upstream's own
        // `CTCF_12878_5M` fixture that reported `d = 200` where upstream infers
        // `d = 101`, and 36770 peaks where upstream calls 40656. Silent wrong numbers are
        // worse than a refusal, so the model is now actually fitted.
        let ext = if o.flag("nomodel") {
            if extsize <= 0 {
                return Err(MacsError::InvalidParameter(
                    "--extsize is required with --nomodel".into(),
                ));
            }
            extsize as f64
        } else {
            let options = macs_model::ModelOptions {
                gsize,
                mfold: mfold_pair,
                bw: o.int("bw").unwrap_or(300),
                d_min: o.int("d_min").unwrap_or(20),
            };
            match macs_model::PeakModel::build(&treat, options) {
                Ok(model) => {
                    // `options.d = peakmodel.d` keeps upstream's *float*; only the
                    // `# d = %d` header line truncates for display.
                    let mut d = model.d;
                    // `if options.d <= 2*options.tsize` (`callpeak_cmd.py:192`): a
                    // fragment shorter than two read lengths is treated as suspect.
                    if d <= 2.0 * tsize as f64 {
                        eprintln!(
                            "#2 Since the d ({d:.0}) calculated from paired-peaks are smaller \
                             than 2*tag length, it may be influenced by unknown sequencing problem!",
                        );
                        if o.flag("onauto") {
                            d = extsize as f64;
                            eprintln!(
                                "#2 MACS will use {d:.0} as EXTSIZE/fragment length d. NOTE: if \
                                 the d calculated is still acceptable, please do not use \
                                 --fix-bimodal option!",
                            );
                        } else {
                            let alt = model
                                .alternative_d
                                .iter()
                                .map(|x| x.to_string())
                                .collect::<Vec<_>>()
                                .join(",");
                            eprintln!(
                                "#2 You may need to consider one of the other alternative d(s): {alt}"
                            );
                            eprintln!(
                                "#2 You can restart the process with --nomodel --extsize XXX with \
                                 your choice or an arbitrary number. Nontheless, MACS will \
                                 continute computing."
                            );
                        }
                    }
                    // `#2.2 Generate R script for model : <outdir>/<name>_model.r`
                    // (`OptValidator.py:174`, `callpeak_cmd.py:188-189`).
                    let model_r = outdir.join(format!("{name}_model.r"));
                    macs_model::model2r_script(&model, &model_r, &name)?;
                    alt_d_line = Some(format!(
                        "# alternative fragment length(s) may be {} bps\n",
                        model
                            .alternative_d
                            .iter()
                            .map(|x| x.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    ));
                    d
                }
                Err(MacsError::NotEnoughPairs { found, needed }) => {
                    // `except NotEnoughPairsException: if not options.onauto: sys.exit(1)`
                    // (`callpeak_cmd.py:202-208`). Upstream writes **no** output file here.
                    if !o.flag("onauto") {
                        return Err(MacsError::NotEnoughPairs { found, needed });
                    }
                    eprintln!("#2 Skipped...");
                    let d = extsize as f64;
                    eprintln!(
                        "#2 Since --fix-bimodal is set, MACS will use {d:.0} as fragment length"
                    );
                    d
                }
                Err(e) => return Err(e),
            }
        };
        let cfg = SeConfig {
            // `d` is the fragment length: `--extsize` under `--nomodel`, otherwise the
            // value fitted by the `PeakModel`. Upstream keeps it a float here
            // (`options.d = peakmodel.d`) and truncates only when printing `# d = %d`,
            // so the pileup extension is the *truncated* integer, exactly as upstream's
            // Cython `int` parameters receive it.
            extsize: ext.trunc() as i64,
            scaleto_large: o.get("scaleto").unwrap_or("small") == "large",
            gsize: cfg_common.gsize,
            slocal: cfg_common.slocal,
            llocal: cfg_common.llocal,
            qvalue: cfg_common.qvalue,
            call_summits: cfg_common.call_summits,
            broad: cfg_common.broad,
            broad_cutoff: cfg_common.broad_cutoff,
            nolambda: cfg_common.nolambda,
            // F188: `maxgap = opt.maxgap or opt.tsize`; callpeak has no
            // `--max-gap`, so it is the measured tag size.
            max_gap: tsize.max(1) as macs_core::Coord,
            // F189: see `ChromCall::p_cutoff`.
            p_cutoff: cfg_common.p_cutoff,
            // F160: `--shift` moves every 5' end before the single-end extension.
            // It is ignored in paired-end mode because neither `PETrackI`
            // pileup takes a shift argument.
            end_shift: if is_pe {
                0
            } else {
                o.int("shift").unwrap_or(0)
            },
        };
        let (setup_se, lambda_bg, _names) =
            macs_peaks::callpeak::se_setup(&treat, ctrl.as_ref(), &cfg);
        lambda_ladder = Some(se_lambda_ladder_json(&setup_se, &cfg, ctrl.is_some()));
        scaling_setup = Some(setup_se);
        let has_ctrl = ctrl.is_some();
        let sigs = if sdump_on {
            // The stage harness intentionally records full tracks. Keep its
            // materialized path isolated from normal peak calling.
            macs_peaks::callpeak::build_all_se_chromosomes(
                &_names,
                &treat,
                ctrl.as_ref(),
                &cfg,
                &setup_se,
            )
        } else {
            se_source = Some(SeSignalSource {
                treat,
                ctrl,
                cfg,
                setup: setup_se,
                names: _names,
            });
            Vec::new()
        };
        if sdump_on {
            stagedump::record_reads(
                &mut sdump,
                sdump_pre.as_deref(),
                has_ctrl,
                ext,
                t1,
                c1,
                cfg_common.slocal,
                cfg_common.llocal,
            );
            // `-B` makes upstream persist these two as `<name>_treat_pileup.bdg` and
            // `<name>_control_lambda.bdg`, so dumping them here gives the stage
            // differential something real to compare against the oracle's own bytes.
            // The SE branch never shifts (`coord_shift` is 0 there: it is only
            // non-zero for a counted FRAG run), so zero is the correct value here and is
            // not merely a placeholder.
            //
            // The control pileup is deliberately left empty: upstream's `-B` writes the
            // control *lambda*, not the control pileup, so there is no oracle counterpart
            // to compare against and emitting it would register as "present only in ours"
            // on every fixture forever.
            stagedump::record_tracks(&mut sdump, &signal_tracks(&sigs), &[], 0);
            record_lambda_merged_stage(&mut sdump, &sigs, 0, spmr_denominator(o, t1, c1));
        }
        (sigs, ext.trunc().max(1.0) as u64, false, lambda_bg, 0i64)
    };
    let _ = names_of;
    stagedump::record_duplicates(
        &mut sdump,
        t0,
        t1,
        c0,
        c1,
        !ctrl_paths.is_empty(),
        format == "FRAG",
    );

    // ---- call -------------------------------------------------------------
    // F165: a counted (`--format FRAG`) run is computed with its coordinates
    // shifted, but upstream's AFDR histogram measures the first span from
    // coordinate **zero**. Feeding the shifted start inflates `N` by the shift
    // and moves every q-value by `log10(1 + shift / span)`.
    // The per-chromosome q-score tracks are read in exactly one place, the
    // `--broad` level-2 cutoff, so a narrow-mode run drops them.
    //
    // For single-end runs the signal source rebuilds one chromosome at a time.
    // `hist_pairs` is the p->q histogram, and is only materialised when a stage dump
    // was requested; see `build_qtable_from_with_hist`.
    // `--cutoff-analysis` is not just a report: the ladder's cutoffs are seeded into
    // the AFDR histogram before the q-table is built, so the q-values differ from a run
    // without the flag. That is why the flag is threaded into the table pass rather than
    // being a post-hoc dump.
    let cut_analysis = o.flag("cutoff_analysis");
    let (table, qtracks, cut_stats) = if let Some(source) = &se_source {
        use std::io::{BufWriter, Write as _};
        let mut histogram = macs_score::PScoreHistogram::new();
        let ladder = cut_analysis.then(macs_peaks::callpeak::cutoff_ladder);
        let mut cut_stats = macs_peaks::callpeak::CutoffStats::new();
        let mut cache = PScoreCache::new();
        let spool = tempfile::NamedTempFile::new_in(&tempdir)?;
        let mut spool_writer = BufWriter::new(spool.as_file());
        let mut stream_bodies = if o.flag("store_bdg") {
            Some((
                tempfile::NamedTempFile::new_in(&tempdir)?,
                tempfile::NamedTempFile::new_in(&tempdir)?,
            ))
        } else {
            None
        };
        let denominator: f64 = if o.flag("do_SPMR") {
            if t1 as f64 <= c1 as f64 * 2.0 {
                t1 as f64 / 1e6
            } else {
                c1 as f64 / 1e6
            }
        } else {
            1.0
        };
        for name_bytes in &source.names {
            let Some(signal) = source.chromosome(name_bytes) else {
                spool_writer.write_all(&[0])?;
                continue;
            };
            spool_writer.write_all(&[1])?;
            write_signal_track(&mut spool_writer, &signal.treat)?;
            if let Some(ctrl) = &signal.ctrl {
                spool_writer.write_all(&[1])?;
                write_signal_track(&mut spool_writer, ctrl)?;
            } else {
                spool_writer.write_all(&[0])?;
            }
            cache.clear();
            let ptrack = macs_peaks::callpeak::se_pscore_track(&signal, &mut cache);
            histogram.add_track_from(&ptrack, 0);
            if let Some(ladder) = &ladder {
                macs_peaks::callpeak::accumulate_cutoffs(
                    &mut cut_stats,
                    &ptrack,
                    ladder,
                    tsize.max(1) as macs_core::Coord,
                    d,
                );
            }
            if let (Some((tfile, cfile)), Some(_)) = (&mut stream_bodies, &signal.ctrl) {
                let mut t_body = String::new();
                let mut c_body = String::new();
                append_paired_bdg_signal(
                    &signal,
                    coord_shift,
                    denominator,
                    &mut t_body,
                    &mut c_body,
                );
                tfile.write_all(t_body.as_bytes())?;
                cfile.write_all(c_body.as_bytes())?;
            }
        }
        if let Some(ladder) = &ladder {
            macs_peaks::callpeak::seed_cutoffs(&mut histogram, ladder);
        }
        if let Some((tfile, cfile)) = &mut stream_bodies {
            tfile.as_file_mut().sync_all()?;
            cfile.as_file_mut().sync_all()?;
        }
        spool_writer.flush()?;
        drop(spool_writer);
        spool.as_file().sync_all()?;
        let table = macs_score::PqTable::from_histogram(&histogram);
        let cut_stats = if cut_analysis {
            cut_stats
        } else {
            macs_peaks::callpeak::CutoffStats::default()
        };
        stream_bdg = stream_bodies;
        se_signal_spool = Some(spool);
        (table, Vec::new(), cut_stats)
    } else if let Some(source) = &mut pe_source {
        use std::io::{BufWriter, Write as _};
        let mut histogram = macs_score::PScoreHistogram::new();
        let ladder = cut_analysis.then(macs_peaks::callpeak::cutoff_ladder);
        let mut cut_stats = macs_peaks::callpeak::CutoffStats::new();
        let mut cache = PScoreCache::new();
        let spool = tempfile::NamedTempFile::new_in(&tempdir)?;
        let mut spool_writer = BufWriter::new(spool.as_file());
        let mut stream_bodies = if o.flag("store_bdg") {
            Some((
                tempfile::NamedTempFile::new_in(&tempdir)?,
                tempfile::NamedTempFile::new_in(&tempdir)?,
            ))
        } else {
            None
        };
        let denominator: f64 = if o.flag("do_SPMR") {
            if t1 as f64 <= c1 as f64 * 2.0 {
                t1 as f64 / 1e6
            } else {
                c1 as f64 / 1e6
            }
        } else {
            1.0
        };
        let names = source.names.clone();
        for name_bytes in &names {
            let result = source.chromosome(name_bytes);
            let Some(signal) = result.signals.into_iter().next() else {
                spool_writer.write_all(&[0])?;
                continue;
            };
            spool_writer.write_all(&[1])?;
            write_signal_track(&mut spool_writer, &signal.treat)?;
            if let Some(ctrl) = &signal.ctrl {
                spool_writer.write_all(&[1])?;
                write_signal_track(&mut spool_writer, ctrl)?;
            } else {
                spool_writer.write_all(&[0])?;
            }
            cache.clear();
            let ptrack = macs_peaks::callpeak::se_pscore_track(&signal, &mut cache);
            histogram.add_track_from(&ptrack, coord_shift);
            if let Some(ladder) = &ladder {
                macs_peaks::callpeak::accumulate_cutoffs(
                    &mut cut_stats,
                    &ptrack,
                    ladder,
                    tsize.max(1) as macs_core::Coord,
                    d,
                );
            }
            if let (Some((tfile, cfile)), Some(_)) = (&mut stream_bodies, &signal.ctrl) {
                let mut t_body = String::new();
                let mut c_body = String::new();
                append_paired_bdg_signal(
                    &signal,
                    coord_shift,
                    denominator,
                    &mut t_body,
                    &mut c_body,
                );
                tfile.write_all(t_body.as_bytes())?;
                cfile.write_all(c_body.as_bytes())?;
            }
        }
        if let Some(ladder) = &ladder {
            macs_peaks::callpeak::seed_cutoffs(&mut histogram, ladder);
        }
        if let Some((tfile, cfile)) = &mut stream_bodies {
            tfile.as_file_mut().sync_all()?;
            cfile.as_file_mut().sync_all()?;
        }
        spool_writer.flush()?;
        drop(spool_writer);
        spool.as_file().sync_all()?;
        let table = macs_score::PqTable::from_histogram(&histogram);
        let cut_stats = if cut_analysis {
            cut_stats
        } else {
            macs_peaks::callpeak::CutoffStats::default()
        };
        stream_bdg = stream_bodies;
        se_signal_spool = Some(spool);
        (table, Vec::new(), cut_stats)
    } else {
        let qparts = macs_peaks::callpeak::build_qtable_from_with_hist(
            &signals,
            paired_boundaries,
            coord_shift,
            cfg_common.broad,
            sdump_on,
            macs_peaks::callpeak::CutoffParams {
                enabled: cut_analysis,
                // F188: `maxgap = opt.maxgap or opt.tsize`, and callpeak defines no
                // `--max-gap`, so the merge gap is the measured tag size -- the same
                // value the peak caller uses. `min_length` is `opt.d`.
                max_gap: tsize.max(1) as macs_core::Coord,
                min_length: d,
            },
        );
        (qparts.table, qparts.qtracks, qparts.cutoffs)
    };
    let (table, qtracks) = (&table, &qtracks);
    let cut_stats = &cut_stats;

    if rss_trace {
        eprintln!("rss[signals built] {} kB", rss_kb());
        eprintln!("rss[after qtable] {} kB", rss_kb());
    }
    if sdump_on {
        let cutoff_text = macs_peaks::callpeak::render_cutoff_analysis(
            cut_stats,
            &macs_peaks::callpeak::cutoff_ladder(),
            table,
        );
        sdump.put("qvalue_table", &stagedump::cutoff_table_json(&cutoff_text));
    }
    if sdump_on {
        if let Some(ladder) = lambda_ladder {
            sdump.put("lambda_ladder", &ladder);
        }
        // `scaling` records the treatment/control ratio and the two local windows --
        // the inputs to `lambda_bg`, so a p-score divergence can be traced to here.
        if let Some(setup) = scaling_setup {
            // Exactly upstream's four keys: an extra key here reads as a spurious
            // mismatch rather than as extra information.
            sdump.put(
                "scaling",
                &stagedump::mixed_json(&[
                    (
                        "ratio_treat2control",
                        stagedump::num(setup.ratio_treat2control),
                    ),
                    ("sregion", stagedump::num(setup.sregion as f64)),
                    ("lregion", stagedump::num(setup.lregion as f64)),
                    ("tocontrol", setup.to_control.to_string()),
                ]),
            );
        }
        sdump.finish(&[]);
    }

    let mut peaks: Vec<(String, macs_peaks::Peak)> = Vec::new();
    // the strong (lvl1) sub-peaks inside each broad peak, in the same order --
    // the gappedPeak block structure is derived from them
    let mut broad_lvl1: Vec<Vec<(macs_core::Coord, macs_core::Coord)>> = Vec::new();
    // Chromosome-level parallelism. The porting plan requires it, and it is
    // **output-identical to one thread by construction**:
    //
    // * `signals` is a slice, so `par_iter().map(..).collect::<Vec<_>>()` returns
    //   the per-chromosome results in slice order regardless of how the work was
    //   scheduled, and the flatten below preserves that order;
    // * each chromosome gets its **own** `PScoreCache` via `map_init`. The cache
    //   is pure memoisation (`HashMap` keyed on the f32 bit pattern of lambda), so
    //   a fresh one per chromosome cannot change any score;
    // * the only cross-chromosome state, the AFDR histogram, is still accumulated
    //   sequentially below in chromosome order -- and its buckets are `i64`
    //   counters, so even a parallel reduce would be exact and order-independent.
    let per_chrom: Vec<(String, Vec<macs_peaks::callpeak::Called>, bool)> = if let Some(source) =
        &se_source
    {
        let mut out = Vec::with_capacity(source.names.len());
        let mut cache = PScoreCache::new();
        let spool = se_signal_spool
            .as_ref()
            .expect("streamed SE signals are spooled during qtable construction");
        let mut spool_reader = std::io::BufReader::new(std::fs::File::open(spool.path())?);
        for name_bytes in &source.names {
            let chrom = source
                .treat
                .genome()
                .get(name_bytes)
                .expect("signal source chromosome is interned in its genome");
            let Some(signal) = read_spooled_signal(&mut spool_reader, name_bytes, chrom)? else {
                continue;
            };
            cache.clear();
            let qtrack = if cfg_common.broad {
                let p = macs_peaks::callpeak::se_pscore_track(&signal, &mut cache);
                macs_score::qscore_track(signal.chrom, &p, table, None)
            } else {
                EMPTY_TRACK.clone()
            };
            cache.clear();
            let cc = ChromCall {
                name: &signal.name,
                chrom: signal.chrom,
                treat: &signal.treat,
                ctrl: signal.ctrl.as_ref(),
                clamp_floor: 0,
                zero_coord: 0,
                qtrack: &qtrack,
                table,
                d,
                max_gap: tsize.max(1) as macs_core::Coord,
                p_cutoff: cfg_common.p_cutoff,
                qvalue: cfg_common.qvalue,
                broad: cfg_common.broad,
                broad_cutoff: cfg_common.broad_cutoff,
                call_summits: cfg_common.call_summits,
                lambda_bg,
            };
            let called = call_chromosome(&cc, &mut cache);
            let bad = cache.hit_bad_lambda();
            out.push((signal.name, called, bad));
        }
        out
    } else if let Some(source) = &pe_source {
        let mut out = Vec::with_capacity(source.names.len());
        let mut cache = PScoreCache::new();
        let spool = se_signal_spool
            .as_ref()
            .expect("streamed PE signals are spooled during qtable construction");
        let mut spool_reader = std::io::BufReader::new(std::fs::File::open(spool.path())?);
        for name_bytes in &source.names {
            let chrom = source
                .treat
                .genome()
                .get(name_bytes)
                .expect("paired signal chromosome is interned in its treatment genome");
            let Some(signal) = read_spooled_signal(&mut spool_reader, name_bytes, chrom)? else {
                continue;
            };
            cache.clear();
            let qtrack = if cfg_common.broad {
                let p = macs_peaks::callpeak::se_pscore_track(&signal, &mut cache);
                macs_score::qscore_track(signal.chrom, &p, table, None)
            } else {
                EMPTY_TRACK.clone()
            };
            cache.clear();
            let cc = ChromCall {
                name: &signal.name,
                chrom: signal.chrom,
                treat: &signal.treat,
                ctrl: signal.ctrl.as_ref(),
                clamp_floor: coord_shift.max(0) as macs_core::Coord,
                zero_coord: if cfg_common.nolambda {
                    coord_shift.max(0) as macs_core::Coord
                } else {
                    0
                },
                qtrack: &qtrack,
                table,
                d,
                max_gap: tsize.max(1) as macs_core::Coord,
                p_cutoff: cfg_common.p_cutoff,
                qvalue: cfg_common.qvalue,
                broad: cfg_common.broad,
                broad_cutoff: cfg_common.broad_cutoff,
                call_summits: cfg_common.call_summits,
                lambda_bg,
            };
            let called = call_chromosome(&cc, &mut cache);
            let bad = cache.hit_bad_lambda();
            out.push((signal.name, called, bad));
        }
        out
    } else {
        configure_pool().install(|| {
            use rayon::prelude::*;
            signals
                .par_iter()
                .enumerate()
                .map_init(PScoreCache::new, |cache: &mut PScoreCache, (k, s)| {
                    let cc = ChromCall {
                        name: &s.name,
                        chrom: s.chrom,
                        treat: &s.treat,
                        ctrl: s.ctrl.as_ref(),
                        // F260: counted paired-end tracks add `coord_shift` on the way
                        // in so `u64` coordinates can stand in for upstream's signed
                        // ones, so upstream's clamp-at-0 is a clamp at `coord_shift`.
                        clamp_floor: coord_shift.max(0) as macs_core::Coord,
                        // F271: upstream's literal 0 in the first-above-cutoff branch is
                        // `coord_shift` in our frame only under `--nolambda`, where the
                        // control side is a one-element array and the paired array has no
                        // real control cursor. Zero on every other path.
                        zero_coord: if cfg_common.nolambda {
                            coord_shift.max(0) as macs_core::Coord
                        } else {
                            0
                        },
                        qtrack: qtracks.get(k).unwrap_or(&EMPTY_TRACK),
                        table,
                        d,
                        max_gap: tsize.max(1) as macs_core::Coord,
                        p_cutoff: cfg_common.p_cutoff,
                        qvalue: cfg_common.qvalue,
                        broad: cfg_common.broad,
                        broad_cutoff: cfg_common.broad_cutoff,
                        call_summits: cfg_common.call_summits,
                        lambda_bg,
                    };
                    let called = call_chromosome(&cc, cache);
                    let bad = cache.hit_bad_lambda();
                    (s.name.clone(), called, bad)
                })
                .collect()
        })
    };
    for (name, called, bad) in per_chrom {
        // F172: upstream raises `ZeroDivisionError: float division` from
        // `PeakDetect.__call_peaks_w_control` when the control lambda is zero at
        // some position (`se_edge/contig_edges` is the recorded example). That is
        // a runtime error with exit status 1 and **no output file**.
        if bad {
            return Err(MacsError::Rejected(
                "callpeak: local lambda is zero at some position, which makes the \
                 p-value undefined (upstream raises ZeroDivisionError here)"
                    .into(),
            ));
        }
        for called in called {
            broad_lvl1.push(called.lvl1);
            peaks.push((name.clone(), called.peak));
        }
    }

    let rows: Vec<XlsRow> = peaks
        .iter()
        .map(|(chrom, pk)| XlsRow {
            chrom: chrom.clone(),
            // F154: a counted (`--format FRAG`) run is computed with its
            // coordinates shifted by `coord_shift`, because upstream's centred
            // control windows can start before the contig. Undo it here, in
            // `i64`, so `gmini_mfrag_d1200_w600_ctrl_069` reports its xls `start`
            // of `-13` instead of wrapping.
            start: pk.start as i64 - coord_shift,
            end: pk.end as i64 - coord_shift,
            summit: pk.summit as i64 - coord_shift,
            pileup: pk.pileup,
            pscore: pk.pscore,
            fold_change: pk.fold_change,
            qscore: pk.qscore,
        })
        .collect();

    std::fs::create_dir_all(&outdir)?;
    // `callpeak_cmd.py:273-290`: version line, the argument block, a blank
    // line, the per-track counts, the shift line, `d`, the alternative-d line
    // (absent under `--nomodel`), the `--nolambda` line, then the table.
    let mut body = String::new();
    body.push_str(&format!(
        "# This file is generated by MACS version {VERSION}\n"
    ));
    body.push_str(&argtxt(o, &name, gsize, is_pe));
    body.push('\n');
    // F163: the `# total %ss in control` line is emitted whenever a control
    // exists; only the four lines *after* it ("after filtering", "maximum
    // duplicate", "Redundant rate") are gated on `keepduplicates != "all"`
    // (`callpeak_cmd.py:126-155`). Tying the whole block to `max_dup_line`
    // dropped the total line for every `--keep-dup all` and every FRAG run --
    // which is exactly where FRAG forces `all`.
    let ctrl_counts = if ctrl_paths.is_empty() {
        None
    } else {
        Some((c0, c1))
    };
    body.push_str(&tagsinfo(is_pe, tsize, t0, t1, max_dup_line, ctrl_counts));
    if !is_pe {
        // F148: PE mode forces `shift = 0`, so the line never appears there.
        let shift = o.int("shift").unwrap_or(0);
        if shift > 0 {
            body.push_str(&format!(
                "# Sequencing ends will be shifted towards 3' by {shift} bp(s)\n"
            ));
        } else if shift < 0 {
            body.push_str(&format!(
                "# Sequencing ends will be shifted towards 5' by {} bp(s)\n",
                -shift
            ));
        }
    }
    body.push_str(&format!("# d = {d}\n"));
    if let Some(line) = &alt_d_line {
        body.push_str(line);
    }
    if cfg_common.nolambda {
        body.push_str("# local lambda is disabled!\n");
    }
    body.push_str(&macs_io::peakout::xls_body(&rows, &name, cfg_common.broad)?);
    // The report is written alongside the peaks, in the same order upstream writes it
    // (`cutpeak_cmd.py` writes the xls first, then the analysis). `min_length` for the
    // analysis is `opt.d` and the merge gap is `opt.maxgap or opt.tsize` -- the same two
    // the peak caller uses -- which is why the help text warns that `minlen` and
    // `maxgap` affect the results.
    if cut_analysis {
        let ladder = macs_peaks::callpeak::cutoff_ladder();
        let body = macs_peaks::callpeak::render_cutoff_analysis(cut_stats, &ladder, table);
        std::fs::write(outdir.join(format!("{name}_cutoff_analysis.txt")), body)?;
    }

    std::fs::write(outdir.join(format!("{name}_peaks.xls")), body)?;

    // F189: `callpeak_cmd.py:298` picks the score column the same way it picks the
    // scoring function: `score_column = "pscore"` under `-p`, `"qscore"` otherwise,
    // and every BED-family writer emits `int(10 * peak[score_column])`.
    let score_col = if cfg_common.p_cutoff.is_some() {
        ScoreCol::P
    } else {
        ScoreCol::Q
    };
    if cfg_common.broad {
        std::fs::write(
            outdir.join(format!("{name}_peaks.broadPeak")),
            macs_io::peakout::broadpeak_body(&rows, &name, score_col),
        )?;
        // the gappedPeak block spans, i.e. `start = xls start - 1`, `end = xls end`
        let blocks: Vec<macs_io::peakout::BroadBlocks> = rows
            .iter()
            .zip(broad_lvl1.iter())
            .map(|(r, lvl1)| {
                macs_io::peakout::broad_blocks(
                    r.start - 1,
                    r.end,
                    &lvl1
                        .iter()
                        // F171: the lvl1 sub-peaks are in shifted coordinates
                        // too, so they need the same correction the peak rows get.
                        // Without it a counted broad peak reports phantom 1 bp
                        // flank blocks and a block start of 5000.
                        .map(|(s, e)| (*s as i64 - coord_shift - 1, *e as i64 - coord_shift))
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        std::fs::write(
            outdir.join(format!("{name}_peaks.gappedPeak")),
            macs_io::peakout::gappedpeak_body(&rows, &blocks, &name, score_col),
        )?;
    } else {
        // F189: `callpeak_cmd.py:298` picks the score column the same way it picks
        // the scoring function: `score_column = "pscore"` under `-p`, `"qscore"`
        // otherwise, and both writers emit `int(10 * peak[score_column])`.
        std::fs::write(
            outdir.join(format!("{name}_peaks.narrowPeak")),
            narrowpeak_body(&rows, &name, score_col),
        )?;
        std::fs::write(
            outdir.join(format!("{name}_summits.bed")),
            summits_body(&rows, &name, score_col),
        )?;
    }

    // `-B`: the treatment pileup and control lambda bedGraphs. Both are written
    // from the paired union walk (`chr_pos_treat_ctrl`), not from the raw
    // pileups -- see `write_paired_bdg`.
    if o.flag("store_bdg") {
        if let Some((treat_temp, ctrl_temp)) = stream_bdg.take() {
            // Keep the full bedGraph bodies on disk until all chromosome calls
            // have passed validation. A zero-lambda error therefore leaves no
            // final bedGraph behind, and the streaming path never buffers the
            // genome-wide text in RAM.
            std::fs::copy(
                treat_temp.path(),
                outdir.join(format!("{name}_treat_pileup.bdg")),
            )?;
            std::fs::copy(
                ctrl_temp.path(),
                outdir.join(format!("{name}_control_lambda.bdg")),
            )?;
        } else {
            write_paired_bdg(
                o,
                &signals,
                &outdir,
                &name,
                SpmrDenom {
                    treat_total: t1,
                    ctrl_total: c1,
                    // `--scaleto small` scales treatment down exactly when its
                    // (post-filter) read count exceeds double the control's
                    treat_scale_is_one: t1 as f64 <= c1 as f64 * 2.0,
                },
                coord_shift,
            )?;
        }
    }

    eprintln!(
        "callpeak: {} peaks written to {}",
        rows.len(),
        outdir.display()
    );

    Ok(())
}

/// `__write_bedGraph_for_a_chromosome` (`CallPeakUnit.py:1717`).
///
/// The treatment pileup and the control lambda are **two coalesced views of the
/// same paired union walk**, not serialisations of the two pileups. Four details
/// are byte-visible:
///
/// * there is **no `track` line** -- the file is written with bare `fprintf`s to
///   a file handle, unlike `bedGraphIO.write_bedGraph` which does emit one;
/// * the run start is a running `pre` initialised to `0`, not the track's own
///   start, so the first row spans `[0, pos[0])`;
/// * a new row is emitted only when the value changed by **more than `1e-5`**
///   (the comment calls this "precision is 5 digits") -- runs within that
///   tolerance are merged, which is what turns ~2800 union rows into ~336;
/// * **the position lags the values by one row.** The treatment and control
///   pointers are advanced once *before* the loop and again inside it, while
///   `pos_array_ptr` is only advanced inside it, so step `i` reads
///   `treat_array[i]`/`ctrl_array[i]` against `pos_array[i-1]`:
///
///   ```c
///   pre_p_t = 0;  pre_v_t = treat_array[0];  treat_array_ptr += 1;  ctrl_array_ptr += 1;
///   for (i = 1; i < l; i++) {
///       v_t = treat_array[i];  v_c = ctrl_array[i];  p = pos_array[i-1];
///       if (fabsf(pre_v_t - v_t) > 1e-5) { fprintf(ft, pre_p_t, p, pre_v_t); ... }
///   }
///   fprintf(ft, pre_p_t, pos_array[l-1], pre_v_t);
///   ```
///
///   Reading `pos_array[i]` instead makes the first row span `[0, pos[1])`
///   rather than `[0, pos[0])` -- for `gonechrom_mpe_d400_w600_ctrl_383` that is
///   `[0, 7)` where upstream writes `[0, 2)`, and for
///   `gmini_mse_d1200_w180_noc_117` `[0, 2)` against `[0, 1)`. Every row after
///   the first is shifted the same way, which is why the run count matched while
///   no byte did;
/// * values are `"%.5f"`, trailing zeros kept, and the final row is always
///   emitted unconditionally.
///
/// `--SPMR` divides by the deeper sample's million-normalised total
/// (`denominator` below); otherwise `denominator` is `1.0`.
/// The flags shared by the single-end and paired-end peak-calling configs.
///
/// `callpeak` builds one of these from the parser and then fills in the
/// mode-specific half (`extsize` for single-end, `tsize` for paired-end).
struct PeLike {
    gsize: f64,
    slocal: i64,
    llocal: i64,
    qvalue: f64,
    call_summits: bool,
    broad: bool,
    broad_cutoff: f64,
    nolambda: bool,
    /// `-log10(--pvalue)` when `-p` was given (F189).
    p_cutoff: Option<f32>,
}

/// The inputs `--SPMR` needs for its denominator (`CallPeakUnit.py:1743`).
///
/// When the treatment was not scaled (`treat_scaling_factor == 1`) the *control*
/// is expressed per million treatment reads, and otherwise the *treatment* is
/// expressed per million control reads -- i.e. both bedGraphs are normalised by
/// the deeper of the two samples.
struct SpmrDenom {
    treat_total: u64,
    ctrl_total: u64,
    treat_scale_is_one: bool,
}

/// Rebuild single-end chromosome signals on demand after pooling reads.
///
/// The inputs stay resident, but the genome-wide pileup and lambda tracks do
/// not. The first pass reduces each chromosome to the q-value histogram; the
/// second pass calls peaks and drops each chromosome before moving on.
struct SeSignalSource {
    treat: SingleEndTrack,
    ctrl: Option<SingleEndTrack>,
    cfg: SeConfig,
    setup: macs_peaks::callpeak::SeSignalSetup,
    names: Vec<Vec<u8>>,
}

impl SeSignalSource {
    fn chromosome(&self, name: &[u8]) -> Option<macs_peaks::callpeak::ChromSignals> {
        let chrom = self.treat.genome().get(name)?;
        macs_peaks::callpeak::build_one_se_chromosome(
            chrom,
            name,
            &self.treat,
            self.ctrl.as_ref(),
            &self.cfg,
            &self.setup,
        )
    }
}

/// Paired-end input tracks remain resident, while only one chromosome's
/// treatment/lambda signals are built at a time and spooled between passes.
struct PeSignalSource {
    treat: FragmentTrack,
    ctrl: Option<FragmentTrack>,
    cfg: PeConfig,
    names: Vec<Vec<u8>>,
    pending: Option<(Vec<u8>, macs_peaks::callpeak::PeResult)>,
}

impl PeSignalSource {
    fn new(
        treat: FragmentTrack,
        ctrl: Option<FragmentTrack>,
        cfg: PeConfig,
    ) -> (Self, macs_core::Coord, bool, f32, i64) {
        let mut names: Vec<Vec<u8>> = treat
            .chroms()
            .iter()
            .map(|&chrom| treat.genome().name(chrom).to_vec())
            .collect();
        names.sort();
        let pending = names.first().map(|name| {
            (
                name.clone(),
                macs_peaks::callpeak::run_callpeak_pe_chromosome(&treat, ctrl.as_ref(), &cfg, name),
            )
        });
        let metadata = if let Some((_, result)) = &pending {
            (
                result.d,
                result.paired_boundaries,
                result.lambda_bg,
                result.coord_shift,
            )
        } else {
            let result = macs_peaks::callpeak::run_callpeak_pe(&treat, ctrl.as_ref(), &cfg);
            (
                result.d,
                result.paired_boundaries,
                result.lambda_bg,
                result.coord_shift,
            )
        };
        (
            Self {
                treat,
                ctrl,
                cfg,
                names,
                pending,
            },
            metadata.0,
            metadata.1,
            metadata.2,
            metadata.3,
        )
    }

    fn chromosome(&mut self, name: &[u8]) -> macs_peaks::callpeak::PeResult {
        if self
            .pending
            .as_ref()
            .is_some_and(|(pending_name, _)| pending_name.as_slice() == name)
        {
            return self.pending.take().expect("pending signal checked").1;
        }
        macs_peaks::callpeak::run_callpeak_pe_chromosome(
            &self.treat,
            self.ctrl.as_ref(),
            &self.cfg,
            name,
        )
    }
}

fn record_pe_scaling_stage(
    sdump: &mut StageDump,
    treat: &FragmentTrack,
    ctrl: Option<&FragmentTrack>,
    cfg: &PeConfig,
) {
    if !sdump.enabled() {
        return;
    }
    let control_total = ctrl.map_or(0, FragmentTrack::total).saturating_mul(2);
    let control_sum = if ctrl.is_some() {
        (control_total as f64 * treat.average_template_length()) as i64 as f64
    } else {
        0.0
    };
    let ratio = if control_sum > 0.0 {
        treat.length() as f64 / control_sum
    } else {
        0.0
    };
    let tocontrol = ctrl.is_some()
        && macs_peaks::callpeak::to_control(cfg.scaleto_large, treat.total(), control_total, true);
    sdump.put(
        "scaling",
        &stagedump::mixed_json(&[
            ("ratio_treat2control", stagedump::num(ratio)),
            ("sregion", stagedump::num(cfg.slocal as f64)),
            ("lregion", stagedump::num(cfg.llocal as f64)),
            ("tocontrol", tocontrol.to_string()),
        ]),
    );
}

fn se_lambda_ladder_json(
    setup: &macs_peaks::callpeak::SeSignalSetup,
    cfg: &macs_peaks::callpeak::SeConfig,
    has_control: bool,
) -> String {
    let mut d_s = Vec::new();
    let mut factors = Vec::new();
    if !cfg.nolambda {
        if has_control {
            d_s.push(cfg.extsize);
            let first = if setup.to_control {
                1.0
            } else {
                setup.ratio_treat2control
            };
            factors.push(first);
            if cfg.slocal != 0 {
                d_s.push(cfg.slocal);
                factors.push(if setup.to_control {
                    cfg.extsize as f64 / cfg.slocal as f64
                } else {
                    cfg.extsize as f64 / cfg.slocal as f64 * setup.ratio_treat2control
                });
            }
            if cfg.llocal > cfg.slocal {
                d_s.push(cfg.llocal);
                factors.push(if setup.to_control {
                    cfg.extsize as f64 / cfg.llocal as f64
                } else {
                    cfg.extsize as f64 / cfg.llocal as f64 * setup.ratio_treat2control
                });
            }
        } else if cfg.llocal > 0 {
            d_s.push(cfg.llocal);
            factors.push(cfg.extsize as f64 / cfg.llocal as f64);
        }
    }
    let treat_scale = if setup.to_control && setup.ratio_treat2control != 0.0 {
        (1.0 / setup.ratio_treat2control) as f32 as f64
    } else {
        1.0
    };
    stagedump::lambda_ladder_json(&d_s, &factors, setup.lambda_bg as f64, treat_scale)
}

fn pe_lambda_ladder_json(
    treat: &FragmentTrack,
    ctrl: Option<&FragmentTrack>,
    cfg: &PeConfig,
    lambda_bg: f32,
) -> String {
    let control_total = ctrl.map_or(0, FragmentTrack::total).saturating_mul(2);
    let control_sum = if ctrl.is_some() {
        (control_total as f64 * treat.average_template_length()) as i64 as f64
    } else {
        0.0
    };
    let ratio = if control_sum > 0.0 {
        treat.length() as f64 / control_sum
    } else {
        0.0
    };
    let tocontrol = ctrl.is_some()
        && macs_peaks::callpeak::to_control(cfg.scaleto_large, treat.total(), control_total, true);
    let mut d_s = Vec::new();
    let mut factors = Vec::new();
    if !cfg.nolambda {
        if ctrl.is_some() {
            let d = cfg.tsize_exact as i64;
            d_s.push(d);
            let ratio_factor = if tocontrol { 1.0 } else { ratio };
            factors.push(ratio_factor);
            if cfg.slocal > 0 {
                d_s.push(cfg.slocal);
                factors.push(cfg.tsize_exact / cfg.slocal as f64 * ratio_factor);
            }
            if cfg.llocal > cfg.slocal && cfg.llocal > 0 {
                d_s.push(cfg.llocal);
                factors.push(cfg.tsize_exact / cfg.llocal as f64 * ratio_factor);
            }
        } else if cfg.llocal > 0 {
            d_s.push(cfg.llocal);
            factors.push(treat.length() as f64 / (cfg.llocal as f64 * treat.total() as f64 * 2.0));
        }
    }
    let treat_scale = if tocontrol && ratio != 0.0 {
        (1.0 / ratio) as f32 as f64
    } else {
        1.0
    };
    stagedump::lambda_ladder_json(&d_s, &factors, lambda_bg as f64, treat_scale)
}

fn record_lambda_merged_stage(
    sdump: &mut StageDump,
    signals: &[macs_peaks::callpeak::ChromSignals],
    coord_shift: i64,
    denominator: f64,
) {
    if !sdump.enabled() {
        return;
    }
    let (mut treat_body, mut ctrl_body) = (String::new(), String::new());
    for signal in signals {
        append_paired_bdg_signal(
            signal,
            coord_shift,
            denominator,
            &mut treat_body,
            &mut ctrl_body,
        );
    }
    sdump.put("lambda_merged", &stagedump::bedgraph_body_json(&ctrl_body));
}

fn spmr_denominator(o: &Options, treat_total: u64, ctrl_total: u64) -> f64 {
    if !o.flag("do_SPMR") {
        return 1.0;
    }
    if treat_total as f64 <= ctrl_total as f64 * 2.0 {
        treat_total as f64 / 1e6
    } else {
        ctrl_total as f64 / 1e6
    }
}

/// Store raw RLE runs so the peak pass can reuse the first pass's pileups
/// without holding every chromosome in memory or recomputing the pileup.
/// Coordinates are u64 and values are written as their exact f32 bit pattern.
fn write_signal_track<W: std::io::Write>(
    writer: &mut W,
    track: &macs_rle::SignalTrack<f32>,
) -> std::io::Result<()> {
    writer.write_all(&track.start().to_le_bytes())?;
    writer.write_all(&track.end().to_le_bytes())?;
    writer.write_all(&(track.runs().len() as u64).to_le_bytes())?;
    for run in track.runs() {
        writer.write_all(&run.end.to_le_bytes())?;
        writer.write_all(&run.value.to_bits().to_le_bytes())?;
    }
    Ok(())
}

fn read_signal_track<R: std::io::Read>(
    reader: &mut R,
    chrom: ChromId,
) -> std::io::Result<macs_rle::SignalTrack<f32>> {
    let mut coord_buf = [0u8; 8];
    reader.read_exact(&mut coord_buf)?;
    let start = u64::from_le_bytes(coord_buf);
    reader.read_exact(&mut coord_buf)?;
    let end = u64::from_le_bytes(coord_buf);
    reader.read_exact(&mut coord_buf)?;
    let run_count = usize::try_from(u64::from_le_bytes(coord_buf)).map_err(|_| {
        std::io::Error::new(std::io::ErrorKind::InvalidData, "too many signal runs")
    })?;
    let mut runs = Vec::with_capacity(run_count);
    for _ in 0..run_count {
        reader.read_exact(&mut coord_buf)?;
        let run_end = u64::from_le_bytes(coord_buf);
        let mut value_buf = [0u8; 4];
        reader.read_exact(&mut value_buf)?;
        runs.push(macs_rle::Run::new(
            run_end,
            f32::from_bits(u32::from_le_bytes(value_buf)),
        ));
    }
    Ok(macs_rle::SignalTrack::from_runs_exact(
        chrom, start, end, runs,
    ))
}

fn read_spooled_signal<R: std::io::Read>(
    reader: &mut R,
    name: &[u8],
    chrom: ChromId,
) -> std::io::Result<Option<macs_peaks::callpeak::ChromSignals>> {
    let mut marker = [0u8; 1];
    reader.read_exact(&mut marker)?;
    if marker[0] == 0 {
        return Ok(None);
    }
    if marker[0] != 1 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "invalid signal spool record marker",
        ));
    }
    let treat = read_signal_track(reader, chrom)?;
    reader.read_exact(&mut marker)?;
    let ctrl = match marker[0] {
        0 => None,
        1 => Some(read_signal_track(reader, chrom)?),
        _ => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "invalid control spool marker",
            ));
        }
    };
    Ok(Some(macs_peaks::callpeak::ChromSignals {
        name: String::from_utf8_lossy(name).into_owned(),
        chrom,
        treat,
        ctrl,
    }))
}

#[cfg(test)]
mod signal_spool_tests {
    use super::*;

    #[test]
    fn binary_spool_preserves_raw_runs_and_optional_control_alignment() {
        let chrom = ChromId(0);
        let treat = macs_rle::SignalTrack::from_runs_exact(
            chrom,
            3,
            30,
            vec![
                macs_rle::Run::new(8, f32::from_bits(0x8000_0000)),
                macs_rle::Run::new(17, 2.5),
                macs_rle::Run::new(30, 2.5),
            ],
        );
        let ctrl = macs_rle::SignalTrack::from_runs_exact(
            chrom,
            0,
            30,
            vec![macs_rle::Run::new(11, 0.25), macs_rle::Run::new(30, 0.5)],
        );
        let mut bytes = Vec::new();
        bytes.push(1); // present
        write_signal_track(&mut bytes, &treat).unwrap();
        bytes.push(1); // control present
        write_signal_track(&mut bytes, &ctrl).unwrap();
        bytes.push(0); // next chromosome absent

        let mut reader = std::io::Cursor::new(bytes);
        let signal = read_spooled_signal(&mut reader, b"chr1", chrom)
            .unwrap()
            .unwrap();
        assert_eq!((signal.treat.start(), signal.treat.end()), (3, 30));
        assert_eq!(
            signal
                .treat
                .runs()
                .iter()
                .map(|run| (run.end, run.value.to_bits()))
                .collect::<Vec<_>>(),
            treat
                .runs()
                .iter()
                .map(|run| (run.end, run.value.to_bits()))
                .collect::<Vec<_>>()
        );
        let decoded_ctrl = signal.ctrl.unwrap();
        assert_eq!((decoded_ctrl.start(), decoded_ctrl.end()), (0, 30));
        assert_eq!(decoded_ctrl.runs(), ctrl.runs());
        assert!(read_spooled_signal(&mut reader, b"chr2", chrom)
            .unwrap()
            .is_none());
    }
}

fn write_paired_bdg(
    o: &Options,
    signals: &[macs_peaks::callpeak::ChromSignals],
    // ^ materialised by the caller; `-B` writes a bedGraph per chromosome, so it
    //   genuinely needs them all at once. That is the caller's decision, not this
    //   function's, and `-B` is already the memory-heavy path.
    outdir: &std::path::Path,
    name: &str,
    sp: SpmrDenom,
    coord_shift: i64,
) -> Result<()> {
    // F191: `denominator` is a Python float in upstream -- `self.treat.total/1e6`
    // is a **double** division -- and it is used against an f32 array element in a
    // Python `/`, so the result is a double that is rounded to f32 once, on the
    // assignment to `float __pyx_v_pre_v_c`. Narrowing the denominator to f32
    // first, and dividing in f32, each lose a bit independently.
    let denominator: f64 = if o.flag("do_SPMR") {
        if sp.treat_scale_is_one {
            sp.treat_total as f64 / 1e6
        } else {
            sp.ctrl_total as f64 / 1e6
        }
    } else {
        1.0
    };

    let mut treat_bdg = String::new();
    let mut ctrl_bdg = String::new();
    for s in signals {
        let Some(c) = &s.ctrl else { continue };
        // F170: the bedGraphs carry the **same** coordinates as the peak rows, so
        // a counted run's shift has to come off here too. Without it
        // `frag_basic/barcode_fragments -B` reports a first breakpoint of `99`
        // where upstream reports `-4901`.
        let shift = |p: macs_core::Coord| -> i64 { p as i64 - coord_shift };
        let (pos, t, ctl) = macs_peaks::callpeak::paired_union(&s.treat, c);
        if pos.is_empty() {
            continue;
        }
        let pos: Vec<i64> = pos.iter().map(|p| shift(*p)).collect();
        if std::env::var_os("MACS3_DEBUG_UNION").is_some() {
            let at: Option<i64> = std::env::var("MACS3_DEBUG_AT")
                .ok()
                .and_then(|v| v.parse().ok());
            for i in 0..pos.len() {
                let keep = match at {
                    Some(a) => (a - 40..=a + 40).contains(&pos[i]),
                    None => i < 10,
                };
                if keep {
                    eprintln!("UNION {}\t{}\t{}\t{}", s.name, pos[i], t[i], ctl[i]);
                }
            }
        }
        coalesce_into(&mut treat_bdg, &s.name, &pos, &t, denominator);
        coalesce_into(&mut ctrl_bdg, &s.name, &pos, &ctl, denominator);
    }
    std::fs::write(outdir.join(format!("{name}_treat_pileup.bdg")), treat_bdg)?;
    std::fs::write(outdir.join(format!("{name}_control_lambda.bdg")), ctrl_bdg)?;
    Ok(())
}

fn append_paired_bdg_signal(
    signal: &macs_peaks::callpeak::ChromSignals,
    coord_shift: i64,
    denominator: f64,
    treat_body: &mut String,
    ctrl_body: &mut String,
) {
    let Some(ctrl) = &signal.ctrl else { return };
    let shift = |p: macs_core::Coord| p as i64 - coord_shift;
    let (raw_pos, treat, control) = macs_peaks::callpeak::paired_union(&signal.treat, ctrl);
    if raw_pos.is_empty() {
        return;
    }
    let pos: Vec<i64> = raw_pos.iter().map(|p| shift(*p)).collect();
    coalesce_into(treat_body, &signal.name, &pos, &treat, denominator);
    coalesce_into(ctrl_body, &signal.name, &pos, &control, denominator);
}

fn coalesce_into(out: &mut String, chrom: &str, pos: &[i64], vals: &[f32], denom: f64) {
    let last = pos.len() - 1;
    let mut pre_p: i64 = 0;
    // `vals[i] / denom` upstream is `float / Python float` -> a **double**, narrowed
    // to f32 exactly once on assignment. Dividing in f32 double-rounds.
    let mut pre_v = (f64::from(vals[0]) / denom) as f32;
    // step i reads vals[i] against pos[i-1], for i in 1..=last
    // (`for i in range(1, l)` upstream: the last step reads `pos[last-1]` and
    // `vals[last]`, which is how the final value change is caught.)
    for (i, v) in vals.iter().enumerate().take(last + 1).skip(1) {
        let v = (f64::from(*v) / denom) as f32;
        if (pre_v - v).abs() > 1e-5 {
            out.push_str(&format!("{chrom}\t{pre_p}\t{}\t{pre_v:.5}\n", pos[i - 1]));
            pre_v = v;
            pre_p = pos[i - 1];
        }
    }
    // The closing row is **unconditional**. When the last value change landed on
    // step `last`, `pre_p` is already `pos[last-1]`, so this writes a second row
    // `[pos[last-1], pos[last])` carrying the same value as its predecessor --
    // a redundant split that a reader cannot detect but that is in the bytes.
    // Merging it instead loses one line: `gmini_mse_d1200_w180_noc_117` ended
    // `[359, 378)` where upstream writes `[359, 361)` and `[361, 378)`.
    out.push_str(&format!("{chrom}\t{pre_p}\t{}\t{pre_v:.5}\n", pos[last]));

    // Diagnostic: the exact f32 values entering `%.5f`, which is otherwise
    // unobservable -- the file only ever records five decimals, so a one-ulp
    // difference is invisible until it flips a printed digit (F191/F236).
    if let Some(path) = std::env::var_os("MACS3_RS_DUMP_CTRL_LAMBDA") {
        use std::io::Write as _;
        let path = std::path::Path::new(&path).join(format!("{chrom}.tsv"));
        if let Ok(mut f) = std::fs::File::create(&path) {
            for i in 0..=last {
                let v = (f64::from(vals[i]) / denom) as f32;
                let _ = writeln!(
                    f,
                    "{}\t{}\t{}\t{:.9e}\t{:08x}",
                    chrom,
                    pos[i],
                    i,
                    v,
                    v.to_bits()
                );
            }
        }
    }
}

#[cfg(test)]
mod tsize_tests {
    use super::load_se_text;
    use std::io::Write;

    fn write_tmp(name: &str, body: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(name);
        let mut f = std::fs::File::create(&p).expect("create fixture");
        f.write_all(body.as_bytes()).expect("write fixture");
        p
    }

    /// F192: `tlen_parse_line` is the **signed** `col3 - col2`, and only lengths `> 0`
    /// are accumulated.
    ///
    /// Taking `|end - pos|` instead counted records with `col3 < col2` using a
    /// positive magnitude. On this input upstream reports `tag size = 100` while the
    /// absolute-value version reported `275`.
    #[test]
    fn reversed_records_do_not_contribute_to_the_tag_size() {
        let p = write_tmp(
            "f192_reversed.bed",
            "chr1\t100\t200\tr1\t0\t+\n\
             chr1\t300\t400\tr2\t0\t+\n\
             chr1\t500\t600\tr3\t0\t+\n\
             chr1\t900\t100\tbad\t0\t+\n",
        );
        let (_track, mean) = load_se_text(&p).expect("load");
        assert_eq!(mean, 100.0, "reversed records must not shift the mean");
        let _ = std::fs::remove_file(p);
    }

    /// The same rule applies on the minus strand, where `col3` is the 5' coordinate:
    /// the length is still `col3 - col2` as written, not a swapped pair.
    #[test]
    fn minus_strand_uses_the_columns_as_written() {
        let p = write_tmp(
            "f192_minus.bed",
            "chr1\t100\t200\tr1\t0\t-\n\
             chr1\t300\t400\tr2\t0\t-\n",
        );
        let (_track, mean) = load_se_text(&p).expect("load");
        assert_eq!(mean, 100.0, "two 100 bp minus tags average to 100");
        let _ = std::fs::remove_file(p);
    }

    /// A minus-strand record with `col3 < col2` is likewise rejected, because the
    /// sign comes from the columns, not from the strand.
    #[test]
    fn reversed_minus_records_are_also_ignored() {
        let p = write_tmp(
            "f192_minus_rev.bed",
            "chr1\t100\t200\tr1\t0\t-\n\
             chr1\t900\t300\tbad\t0\t-\n",
        );
        let (_track, mean) = load_se_text(&p).expect("load");
        assert_eq!(mean, 100.0);
        let _ = std::fs::remove_file(p);
    }

    /// Only the first ten qualifying tags are sampled (`while n < 10 and m < 10000`),
    /// and the mean is truncated. Both properties are load-bearing for the header.
    #[test]
    fn the_mean_is_truncated_and_capped_at_ten_tags() {
        // 12 tags of length 100..111. The first ten give sum 1045, and upstream's
        // `cython.cast(cython.int, s/n)` is integer division, so the report is 104 --
        // not the 12-tag mean of 105.5 and not a rounded 105.
        let mut body = String::new();
        for i in 0..12i64 {
            body.push_str(&format!(
                "chr1\t{}\t{}\tr{}\t0\t+\n",
                i * 1000,
                i * 1000 + 100 + i,
                i
            ));
        }
        let p = write_tmp("f192_ten.bed", &body);
        let (_track, mean) = load_se_text(&p).expect("load");
        assert_eq!(mean, 104.0, "truncated mean of the first ten tags only");
        let _ = std::fs::remove_file(p);
    }

    /// A file with no qualifying tag reports `-1`, not `0` -- upstream's sentinel,
    /// which `callpeak` turns into an error rather than a zero-width pileup.
    #[test]
    fn no_qualifying_tag_reports_minus_one() {
        let p = write_tmp("f192_none.bed", "chr1\t900\t100\tr1\t0\t+\n");
        let (_track, mean) = load_se_text(&p).expect("load");
        assert_eq!(mean, -1.0);
        let _ = std::fs::remove_file(p);
    }
}
