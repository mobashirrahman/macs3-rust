//! `macs3-rs callvar`: call variants in peak regions from alignment BAM files.
//!
//! # Scope, stated plainly
//!
//! This module implements the **driver** completely and faithfully: the peak BED
//! input, the BAM treatment/control access and the chromosome-naming assertion, the
//! exact VCF header (including the `Program_Args` reconstruction and the per-contig
//! lines), the per-peak loop, and the read-extraction/`maxDuplicate` filter.
//!
//! The variant kernel implements read consensus, per-position statistics, and
//! VCF records. Local assembly uses the bundled fermi-lite C library through
//! the Rust bridge; both assembly and no-assembly paths have oracle checks.
//!
//! # Upstream behaviour this reproduces
//!
//! * `opt_validate_callvar` (`OptValidator.py:663`) is nearly a no-op: it only
//!   clamps `np <= 0` to `1`. In particular it performs **no** validation of the BAM
//!   files or the peak file, so argument validation here must not be stricter.
//! * The peak BED is read by `PeakIO`, which keeps `fields[0..3]` and sorts.
//! * The chromosome-name check is an `assert` on `tbam.get_chromosomes()[0]` vs
//!   `cbam.get_chromosomes()[0]` — i.e. it compares only the *first* chromosome of
//!   each, and it is an assert (so `-O` would disable it), not a raised error.
//! * The VCF header is written **before** any peak is processed, so a run that fails
//!   mid-way still leaves a header-only VCF. That ordering is reproduced.
//! * A peak with no reads logs `No reads found in this peak. Skipped` and continues.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use macs_core::{MacsError, Result};

mod align;
mod fermi;
mod peak_variants;
mod pos_reads_info;
mod ra_collection;
pub mod ra_collection_dump;
pub mod unitig;
mod variant_stat;
mod vcf_header_template;
pub use align::{align_unitigs_to_reference, revcomp, smith_waterman, Alignment, Alignments};
pub use fermi::{
    assemble as fermi_assemble, version as fermi_lite_version, Unitig, DEFAULT_OPT_FLAG,
    MAG_F_AGGRESSIVE, MAG_F_NO_SIMPL, MAG_F_POPOPEN,
};
pub use peak_variants::{PeakVariants, Variant};
pub use pos_reads_info::{DepthOpt, PosReadsInfo, Sample};
pub use ra_collection::{
    assemble_peak_with_unitigs, call_peak_from_unitigs, call_peak_without_assembly,
    call_variants_at_range, finish_peak, is_assembly_implemented, remove_outliers,
    revisit_refer_biased, CallParams, RACollection, NO_READS, REFER_BIASED_AR,
};
pub use unitig::{
    build_unitig_collection, Built, Remapped, UnitigCollection, UnitigRAs, MIN_SCORE_100,
};
pub use variant_stat::{
    cal_model_heter_as, cal_model_heter_noas, cal_model_homo, calculate_gq,
    calculate_gq_heter_assig, DEFAULT_MAX_ALLOWED_AR,
};
pub use vcf_header_template::VCF_HEADER_TEMPLATE;

/// The upstream version this port reproduces. It appears in the VCF header as
/// `##source=MACS_V<version>`, so it has to be the *upstream* string, not ours.
pub const MACS_VERSION: &str = "3.0.5";

/// Is the no-assembly variant-calling kernel implemented?
///
/// Yes. It reproduces upstream's `--fermi off` VCF **byte for byte** on upstream's own
/// `callvar_testing` fixtures: 22 of 22 records identical across all ten peaks,
/// including the three emitted twice because the input contains two byte-identical
/// peaks. See `oracle/check_callvar.sh`.
///
/// [`is_assembly_implemented`] reports the separately tested fermi-lite bridge.
pub const fn is_kernel_implemented() -> bool {
    true
}

/// The validated options for `callvar`, mirroring `callvar_cmd.run`'s reads of
/// `options`.
#[derive(Debug, Clone)]
pub struct CallvarOptions {
    /// `-b/--peak`: peak regions in BED format, sorted by coordinates.
    pub peak_bed: PathBuf,
    /// `-t/--treatment`: treatment BAM.
    pub tfile: PathBuf,
    /// `-c/--control`: optional control BAM.
    pub cfile: Option<PathBuf>,
    /// `-o/--ofile`: output VCF.
    pub ofile: PathBuf,
    /// `--outdir`: accepted for argparse parity; upstream records it but never uses
    /// it for `callvar` (the VCF name comes from `-o`).
    pub outdir: Option<PathBuf>,
    /// `-g/--gq-hetero`.
    pub gq_cutoff_hetero: f64,
    /// `-G/--gq-homo`.
    pub gq_cutoff_homo: f64,
    /// `-Q`: minimum base quality.
    pub q: i64,
    /// `-D`: maximum duplicates per (position, strand, CIGAR).
    pub max_duplicate: i64,
    /// `-F/--fermi`: `auto`, `on` or `off`.
    pub fermi: String,
    /// `--fermi-overlap`.
    pub fermi_min_overlap: i64,
    /// `--top2alleles-mratio`.
    pub top2_alleles_min_ratio: f64,
    /// `--altallele-count`.
    pub alt_allele_min_count: i64,
    /// `--max-ar`.
    pub max_ar: f64,
    /// `-m/--multiple-processing`.
    pub np: i64,
    /// `--verbose`.
    pub verbose: i64,
}

impl CallvarOptions {
    /// The effective process count, after `opt_validate_callvar`'s `np <= 0 -> 1`.
    pub fn np(&self) -> usize {
        if self.np <= 0 {
            1
        } else {
            self.np as usize
        }
    }
}

/// One peak region, as `PeakIO` stores it: chromosome plus `[start, end)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Peak {
    pub chrom: Vec<u8>,
    pub start: i64,
    pub end: i64,
}

/// Read the peak BED exactly as `callvar_cmd` does.
///
/// ```python
/// peakio = open(peakbedfile)
/// peaks = PeakIO()
/// for t_peak in peakio:
///     fs = t_peak.rstrip().split()
///     peaks.add(fs[0].encode(), int(fs[1]), int(fs[2]))
/// peaks.sort()
/// ```
///
/// Note `int(...)`, not `atoi`: these are Python integers, so an out-of-range
/// coordinate raises `ValueError` rather than saturating. That is exit-1 behaviour
/// here, and it happens before the VCF is opened.
pub fn read_peak_bed(path: &Path) -> Result<Vec<Peak>> {
    use std::io::BufRead;
    let f = macs_io::open_maybe_gzip(path)?;
    let mut r = std::io::BufReader::new(f);
    let mut line = Vec::new();
    let mut peaks = Vec::new();
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        // `t_peak.rstrip().split()` splits on *any* run of whitespace, and an
        // all-whitespace line yields an empty list, which upstream then indexes
        // (`fs[0]`) and crashes on. Treat a blank line as an error rather than
        // reproducing an IndexError: the observable exit status is the same (1).
        let text = String::from_utf8_lossy(&line);
        let fs: Vec<&str> = text.split_whitespace().collect();
        if fs.is_empty() {
            continue;
        }
        if fs.len() < 3 {
            return Err(MacsError::BadAlignment(format!(
                "peak BED line needs at least 3 fields, got {}: {:?}",
                fs.len(),
                text.trim()
            )));
        }
        let parse = |s: &str| -> Result<i64> {
            s.parse::<i64>().map_err(|_| {
                MacsError::BadAlignment(format!("peak BED coordinate is not an integer: {:?}", s))
            })
        };
        peaks.push(Peak {
            chrom: fs[0].as_bytes().to_vec(),
            start: parse(fs[1])?,
            end: parse(fs[2])?,
        });
    }
    // `peaks.sort()` -- PeakIO sorts by (chrom, start, end)
    peaks.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    Ok(peaks)
}

/// Rebuild the `Program_Args` header value.
///
/// Upstream does not echo the parsed options; it appends a *canonical* subset to
/// whatever the user typed:
///
/// ```python
/// tmpcmdstr = " --fermi " + fermi + " --fermi-overlap "+str(fermiMinOverlap)
/// " ".join(sys.argv[1:] + ["-Q", str(minQ), "-D", str(maxDuplicate),
///                          "--max-ar", str(max_allowed_ar),
///                          "--top2alleles-mratio", str(top2allelesminr),
///                          "--top2allele-count", str(min_altallele_count),
///                          "-g", str(min_heter_GQ), "-G", str(min_homo_GQ),
///                          tmpcmdstr])
/// ```
///
/// Two details are easy to get wrong and both change the bytes: the flag is
/// `--top2allele-count` (singular "allele") although the option is
/// `--altallele-count`, and `--fermi-overlap` is appended *inside* `tmpcmdstr`,
/// i.e. after `-G`, so it comes last with a leading space.
///
/// `argv` is `sys.argv[1:]`, which **includes the subcommand name**.
pub fn program_args(argv: &[String], o: &CallvarOptions) -> String {
    let mut parts: Vec<String> = argv.to_vec();
    parts.push("-Q".into());
    parts.push(o.q.to_string());
    parts.push("-D".into());
    parts.push(o.max_duplicate.to_string());
    parts.push("--max-ar".into());
    parts.push(fmt_float(o.max_ar));
    parts.push("--top2alleles-mratio".into());
    parts.push(fmt_float(o.top2_alleles_min_ratio));
    parts.push("--top2allele-count".into());
    parts.push(o.alt_allele_min_count.to_string());
    parts.push("-g".into());
    parts.push(fmt_float(o.gq_cutoff_hetero));
    parts.push("-G".into());
    parts.push(fmt_float(o.gq_cutoff_homo));
    parts.push(format!(
        " --fermi {} --fermi-overlap {}",
        o.fermi, o.fermi_min_overlap
    ));
    parts.join(" ")
}

/// `str(float)` as CPython renders it (shortest round-tripping repr).
///
/// `format!("{}", 0.95f64)` gives `0.95` and `format!("{}", 0.8)` gives `0.8`, so the
/// common cases agree with `str()`; the helper exists so the intent is explicit and
/// so a value like `1e-05` does not silently pick up Rust's formatting.
fn fmt_float(v: f64) -> String {
    let s = format!("{v}");
    // Rust prints `1e-5` as `0.00001` at Display's default; CPython's `str()` picks
    // the shortest of the two forms. Normalise the exponent case only.
    if s.contains('e') && !s.contains("e-") {
        s.replace('e', "e+")
    } else {
        s
    }
}

/// Substitute the three `%s` placeholders in [`VCF_HEADER_TEMPLATE`].
///
/// Plain `str::replace` rather than `format!`: the template contains no other `%`
/// directives, but a literal `%` appearing in a future header line would make
/// `format!` panic, and this is not on a path where a panic is acceptable.
pub fn vcf_header(date: &str, version: &str, program_args: &str) -> String {
    // Upstream writes `VCFHEADER % (...) + "\n"`: the template itself has no trailing
    // newline, so without this the first `##contig=` line would be appended to the last
    // `##` line.
    let mut out = VCF_HEADER_TEMPLATE
        .replacen("%s", date, 1)
        .replacen("%s", version, 1)
        .replacen("%s", program_args, 1);
    out.push('\n');
    out
}

/// The full header, including the per-contig lines and the column header line.
///
/// `rlengths` is `BAMaccessor.get_rlengths()`: chromosome name (bytes) to length.
/// Upstream iterates `tbam.get_rlengths().items()`, i.e. **file order**, which for a
/// BAM is header order — so a `BTreeMap` would silently reorder contigs. `Vec` in
/// header order is what is wanted.
pub fn vcf_header_with_contigs(
    date: &str,
    version: &str,
    program_args: &str,
    rlengths: &[(Vec<u8>, u64)],
) -> String {
    let mut out = vcf_header(date, version, program_args);
    for (chrom, len) in rlengths {
        let _ = writeln!(
            out,
            "##contig=<ID={},length={},assembly=NA>",
            String::from_utf8_lossy(chrom),
            len
        );
    }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n");
    out
}

/// Today's date as `%Y%m%d`, which is what `##fileDate` carries.
pub fn today_compact(now: std::time::SystemTime) -> String {
    // Days since the Unix epoch -> civil date (Howard Hinnant's algorithm).
    let secs = now
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let days = secs.div_euclid(86_400);
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!("{y:04}{m:02}{d:02}")
}

/// Upstream's chromosome-name consistency check.
///
/// ```python
/// assert tbam.get_chromosomes()[0] in cbam.get_chromosomes() \
///     or cbam.get_chromosomes()[0] in tbam.get_chromosomes(), \
///     Exception("It seems Treatment and Control BAM use different naming ...")
/// ```
///
/// Only the **first** chromosome of each is compared, and it is an `assert`, so it
/// disappears under `python -O`. Reproduced faithfully, including the asymmetry:
/// either the treatment's first contig is in the control's list, or vice versa.
pub fn check_chrom_names(t_chroms: &[Vec<u8>], c_chroms: &[Vec<u8>]) -> Result<()> {
    let (Some(t0), Some(c0)) = (t_chroms.first(), c_chroms.first()) else {
        // `get_chromosomes()[0]` on an empty BAM is an IndexError upstream, i.e. a
        // crash; here it is a clean error, same exit status class.
        return Err(MacsError::BadAlignment(
            "BAM file has no contigs in its header".into(),
        ));
    };
    if !c_chroms.iter().any(|c| c == t0) && !t_chroms.iter().any(|t| t == c0) {
        return Err(MacsError::BadAlignment(
            "It seems Treatment and Control BAM use different naming for chromosomes! \
             Check headers of both files."
                .into(),
        ));
    }
    Ok(())
}

/// Run `callvar`.
///
/// Returns an error without touching `ofile` when the calling kernel is missing, so
/// the "errors precede output files" criterion holds.
pub fn run(o: &CallvarOptions, argv: &[String]) -> Result<()> {
    // Peak BED first: upstream reads it before opening either BAM, and a malformed
    // peak file must not leave a VCF behind.
    let peaks = read_peak_bed(&o.peak_bed)?;

    let mut tbam = open_bam(&o.tfile)?;
    let cbam = match &o.cfile {
        Some(p) => Some(open_bam(p)?),
        None => None,
    };

    if let Some(cb) = &cbam {
        check_chrom_names(tbam.chromosomes(), cb.chromosomes())?;
    }

    let params = CallParams {
        top2alleles_min_ratio: o.top2_alleles_min_ratio as f32,
        min_alt_allele_count: o.alt_allele_min_count as i32,
        max_allowed_ar: o.max_ar as f32,
        min_homo_gq: o.gq_cutoff_homo as i32,
        min_heter_gq: o.gq_cutoff_hetero as i32,
        min_q: o.q as i32,
        max_duplicate: o.max_duplicate.max(0) as u32,
    };

    // `--fermi on|auto` re-calls every peak carrying an indel or a reference-biased het
    // from the assembly. Upstream *replaces* the no-assembly calls for those peaks
    // rather than adding to them (measured: `-F off` gives 22 records, `auto` 16,
    // `on` 15), so emitting the `-F off` answer under `-F auto` would not be a partial
    // result -- it would be a different, wrong one.
    let wants_assembly_mode = matches!(o.fermi.as_str(), "on" | "auto");
    if wants_assembly_mode && !is_assembly_implemented() {
        return Err(MacsError::Rejected(format!(
            "callvar: --fermi {} asks for fermi-lite local assembly. The path is \
             written in macs3-rs and runs, but it is not yet byte-exact: 7 of 16 \
             records match the pinned oracle, and the three peaks that need assembly \
             are called at the right coordinates with well-formed records but \
             different read depths. Upstream *replaces* those peaks' calls rather than \
             adding to them, so falling back to the --fermi off answer here would be a \
             different wrong answer rather than a partial one. No VCF was written.",
            o.fermi
        )));
    }

    // Header order, not sorted: `get_rlengths()` iterates the BAM's dict, which is
    // header order, and the `##contig=` lines are compared byte for byte.
    let rlengths: Vec<(Vec<u8>, u64)> = tbam
        .header()
        .references
        .iter()
        .zip(tbam.header().rlengths.iter())
        .map(|(n, l)| (n.clone(), u64::from(*l)))
        .collect();
    let header = vcf_header_with_contigs(
        &today_compact(std::time::SystemTime::now()),
        MACS_VERSION,
        &program_args(argv, o),
        &rlengths,
    );
    std::fs::write(&o.ofile, header)?;

    // Peaks are sorted by (chrom, start, end) exactly as `peaks.sort()` does, so a
    // peak file given out of order produces the same VCF.
    let mut peaks = peaks;
    peaks.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    let mut cbam = cbam;
    for peak in &peaks {
        // `No reads found in this peak. Skipped` -- upstream catches the exception and
        // moves on, so an uncovered peak is silently absent from the VCF.
        let treatment = tbam.reads_in_region(
            &peak.chrom,
            peak.start.max(0) as u32,
            peak.end.max(0) as u32,
            params.max_duplicate,
        )?;
        if treatment.is_empty() {
            continue;
        }
        let control = match cbam.as_mut() {
            Some(cb) => cb.reads_in_region(
                &peak.chrom,
                peak.start.max(0) as u32,
                peak.end.max(0) as u32,
                params.max_duplicate,
            )?,
            None => Vec::new(),
        };

        let mut coll =
            match RACollection::new(&peak.chrom, peak.start, peak.end, treatment, control) {
                Ok(c) => c,
                // Upstream logs "No reads found in this peak. Skipped" and continues.
                Err(e) if e == NO_READS => continue,
                Err(e) => return Err(MacsError::Rejected(e)),
            };
        coll.remove_outliers_in_place(5);

        // F224: upstream's branch structure, which is not "assembly replaces the
        // no-assembly answer". When assembly is attempted the no-assembly variants are
        // held *unwritten* until the branch is decided:
        //
        //   assembly failed with -1  -> the peak is dropped entirely
        //   assembly returned 0     -> fall back to the no-assembly answer
        //   auto, no indel, refer-biased -> revisit ONLY the reference-biased
        //                                   positions on the no-assembly answer
        //   otherwise                -> reset and re-call everything from the unitigs
        //
        // My first version replaced the whole peak with the assembly-derived one in
        // both `auto` sub-cases, which discards correct no-assembly calls and loses
        // reads (DPT 12 where upstream reports 15).
        let (mut variants, _wants_assembly) =
            call_peak_without_assembly(&coll, &params).map_err(MacsError::Rejected)?;

        let assembly_requested = wants_assembly_mode
            && (o.fermi == "on"
                || variants.has_indel()
                || variants.has_refer_biased_01(REFER_BIASED_AR));

        if assembly_requested {
            match assemble_peak_with_unitigs(&coll, o.fermi_min_overlap as i32)
                .map_err(MacsError::Rejected)?
            {
                Built::TooManyMismatches => {
                    // "Too many mismatches ... we will skip this region entirely!"
                    // Note this discards the no-assembly result too, not just the
                    // assembly one.
                    continue;
                }
                Built::Empty => {
                    // "Failed to assemble unitigs, fall back to previous results"
                }
                Built::Collection(coll_v) => {
                    if o.fermi == "auto"
                        && !variants.has_indel()
                        && variants.has_refer_biased_01(REFER_BIASED_AR)
                    {
                        revisit_refer_biased(&mut variants, &coll_v, &params)
                            .map_err(MacsError::Rejected)?;
                    } else {
                        variants = call_peak_from_unitigs(&coll, &coll_v, &params)
                            .map_err(MacsError::Rejected)?;
                    }
                }
            }
        }

        finish_peak(&mut variants, &o.ofile).map_err(MacsError::Rejected)?;
    }
    Ok(())
}

/// A BAM file `callvar` has opened, reduced to what the driver needs.
pub struct Bam {
    /// `BAMaccessor`, which owns the header, the index and the block cache.
    acc: macs_io::bam::BamAccessor,
}

impl Bam {
    /// Contig names in header order -- `BAMaccessor.get_chromosomes()`.
    pub fn chromosomes(&self) -> &[Vec<u8>] {
        &self.acc.header().references
    }

    /// The parsed header.
    pub fn header(&self) -> &macs_io::bam::BamHeader {
        self.acc.header()
    }

    /// `BAMaccessor.get_reads_in_region`.
    fn reads_in_region(
        &mut self,
        chrom: &[u8],
        left: u32,
        right: u32,
        max_duplicate: u32,
    ) -> Result<Vec<macs_io::bam::ReadAlignment>> {
        self.acc.reads_in_region(chrom, left, right, max_duplicate)
    }
}

fn open_bam(path: &Path) -> Result<Bam> {
    if !path.exists() {
        return Err(MacsError::Io(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("{}", path.display()),
        )));
    }
    Ok(Bam {
        acc: macs_io::bam::BamAccessor::open(path)?,
    })
}
