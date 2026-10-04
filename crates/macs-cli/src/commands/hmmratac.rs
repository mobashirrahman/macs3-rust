//! `macs3-rs hmmratac`: decompose an ATAC-seq fragment set into short /
//! mono- / di- / tri-nucleosomal signal and decode it with a three-state HMM.
//!
//! Port of `hmmratac_cmd.run` (`Commands/hmmratac_cmd.py:50`):
//!
//! ```text
//! fragments ─► fragment-length EM ─► weight mapping ─► 4 digested pileups
//!                                                          │
//!                       ┌──────────────────────────────────┘
//!                       ▼
//!   fold-change pileup ─► candidate regions (prescan, expanded, merged)
//!                                     │
//!                                     ▼
//!                        bin extraction ─► HMM posteriors
//!                                                       │
//!                     ┌─────────────────────────────────┤
//!                     ▼                                 ▼
//!           *_open/nuc/bg.bdg                 state path ─┬─► *_states.bed
//!                                                       └─► *_accessible_regions.narrowPeak
//! ```
//!
//! Everything deterministic is ported. **Inference is exact** against hmmlearn --
//! see `crates/macs-hmmratac/tests/hmm_inference.rs` (worst total-variation
//! distance 1.4e-6) and F143 for why that needed hmmlearn's *inconsistent*
//! matrix orientation. The single declared deviation is **self-training**, which
//! needs hmmlearn's Baum-Welch (an SVD, a QR, and a specific `MT19937` draw);
//! `--model` supplies a trained model instead.
//!
//! # Bins come out in *reverse* chromosome order
//!
//! `extract_value_hmmr` (`BedGraph.py:1105`) builds `common_chr = sorted(...)`
//! and then iterates `for i in range(len(common_chr)): chrom = common_chr.pop()`
//! -- it pops, so chromosomes are visited from last to first. `pos` rows are
//! appended in that order, and `hmmratac_cmd` writes them to the posterior
//! temporary file in the same order, which is the order `generate_states_path`
//! and `save_accessible_regions` then read. So every hmmratac output is ordered by
//! *descending* chromosome name. Reproduced here, because the output files are
//! compared byte-for-byte.
//!
//! # `save_accessible_regions` drops the last region
//!
//! `hmmratac_cmd.py:660` writes `accessible_regions[:-1]`, so the final
//! `nuc`-`open`-`nuc` triple in the state path is never reported. That is an
//! upstream off-by-one, reproduced rather than fixed.
//!
//! # Batching does not change the answer
//!
//! Upstream pops `decoding_steps` (10000) candidate regions at a time purely to
//! bound memory. Each region's bins are its own HMM sequence, so decoding is
//! independent per region; this port batches the same way and concatenates.
//!
//! # `--maxTrain` truncation is deterministic, not `randomly_pick`
//!
//! When more training regions survive than `--maxTrain`, upstream calls
//! `PeakIO.randomly_pick`, which shuffles with **Python's `random` module** -- a
//! different stream from the NumPy one the EM downsample uses. That needs a port
//! of `random.shuffle`; until then the cap keeps the first `--maxTrain` regions in
//! (chromosome, position) order.

use std::collections::{BTreeMap, HashSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};

use macs_bedgraph::peakcall;
use macs_bedgraph::BedGraph;
use macs_core::genome::ChromId;
use macs_core::{Coord, Genome, MacsError, Result};
use macs_hmmratac::em::{prepend_short, round_to_one_decimal, train_em, EmParams};
use macs_hmmratac::{
    accessible_regions_chromed, generate_states_path_chromed, generate_weight_mapping,
    make_bdg_of_bins_from_regions, pileup_bdg_hmmratac, ExtractedBin, HmmType, ModelFile, ProbRow,
    N_SIGNALS,
};
use macs_peaks::hmm_regions::RegionSet;
use macs_rle::SignalTrack;
use macs_track::{filter_frag_dup, sample_frag_percent, FragTrackBuilder, FragmentTrack};

use crate::Options;

/// Output file names, from `hmmratac_cmd.run`'s section 0.
#[derive(Debug)]
struct Outputs {
    outdir: PathBuf,
    name: String,
}

impl Outputs {
    fn new(outdir: &Path, name: &str) -> Self {
        Self {
            outdir: outdir.to_path_buf(),
            name: name.to_string(),
        }
    }
    fn digested(&self, k: &str) -> PathBuf {
        self.outdir
            .join(format!("{}_digested_{}.bdg", self.name, k))
    }
    fn training_regions(&self) -> PathBuf {
        self.outdir
            .join(format!("{}_training_regions.bed", self.name))
    }
    fn state_bdg(&self, k: &str) -> PathBuf {
        self.outdir.join(format!("{}_{}.bdg", self.name, k))
    }
    fn states_bed(&self) -> PathBuf {
        self.outdir.join(format!("{}_states.bed", self.name))
    }
    fn accessible(&self) -> PathBuf {
        self.outdir
            .join(format!("{}_accessible_regions.narrowPeak", self.name))
    }
    fn cutoff_analysis(&self) -> PathBuf {
        self.outdir
            .join(format!("{}_cutoff_analysis.tsv", self.name))
    }
}

/// Read every `-i`/`--input` file into one fragment track.
///
/// `--format FRAG` with `--barcodes` restricts to the listed barcodes, and
/// `--max-count` keeps one row per `(chromosome, barcode)` with summed counts --
/// the counted track that `build_petrack(max_count=...)` builds.
fn load_fragments(
    inputs: &[String],
    format: &str,
    barcodes: Option<&HashSet<Vec<u8>>>,
    max_count: bool,
    cap: u32,
) -> Result<FragmentTrack> {
    // `-f BAMPE` is a BAM, not text. Upstream's `BAMPEParser` streams the BAM and
    // keeps one fragment per proper pair (`tlen` at byte offset 28); routing it
    // through the BEDPE text parser read the BGZF header as text and failed. The
    // shared `bampe_fragments` reader reproduces upstream's filter exactly.
    if format == "BAMPE" {
        let mut b = FragTrackBuilder::new();
        for path in inputs {
            let (frags, _) = macs_io::bam::bampe_fragments(Path::new(path))?;
            for fr in &frags {
                b.push(&fr.chrom, u64::from(fr.start), u64::from(fr.start + fr.len));
            }
        }
        b.finalize();
        return Ok(b.build());
    }
    // `FragParser.build_petrack` always builds a `PETrackII`, so a FRAG input is
    // a counted track whether or not `--max-count` was given; only the count
    // value is capped. BEDPE is a `PETrackI`.
    let mut b = if format == "FRAG" {
        FragTrackBuilder::with_barcodes()
    } else {
        FragTrackBuilder::new()
    };
    for path in inputs {
        let rdr = macs_io::open_maybe_gzip(Path::new(path))?;
        let mut r = std::io::BufReader::new(rdr);
        let mut line: Vec<u8> = Vec::new();
        loop {
            line.clear();
            if r.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            let rec = match format {
                "FRAG" => macs_io::parse_frag_line(&line)?,
                _ => macs_io::parse_bedpe_line(&line)?,
            };
            let Some(rec) = rec else { continue };
            if rec.chrom.is_empty() || rec.left < 0 || rec.right < 0 || rec.right < rec.left {
                continue;
            }
            if let Some(want) = barcodes {
                match rec.barcode.as_ref() {
                    Some(bc) if want.contains(bc) => {}
                    _ => continue,
                }
            }
            if format == "FRAG" {
                // `FragParser.build_petrack` always builds a `PETrackII`, so the
                // fifth column is the multiplicity whether or not `--max-count`
                // is given. `--max-count` only caps it:
                // `if max_count: count = min(count, max_count)` -- it does *not*
                // merge barcode rows.
                let raw = rec.count.unwrap_or(1);
                let count = if max_count { raw.min(cap) } else { raw };
                b.push_with_count(&rec.chrom, rec.left as Coord, rec.right as Coord, count);
            } else {
                // BEDPE is a `PETrackI`: no counts, depth 1 per row.
                b.push(&rec.chrom, rec.left as Coord, rec.right as Coord);
            }
        }
    }
    b.finalize();
    Ok(b.build())
}

/// `PETrackII.pileup_bdg(scale_factor=1, baseline_value=0)` (`PairedEndTrack.py:1513`).
///
/// The pileup is **count-weighted**: `pileup_from_LRC_as_list` passes the `c`
/// column as the weight for both endpoints, so a FRAG row with count 7 contributes
/// depth 7. Treating every row as depth 1 changes the fold-change track, and
/// through it both the candidate regions and the HMM input.
fn pileup_fragments(track: &FragmentTrack) -> BedGraph {
    // F226: the tracks are keyed by `track`'s own ChromIds, so the BedGraph must
    // resolve names through that same dictionary.
    let mut bg = BedGraph::with_genome(track.genome(), 0.0);
    for chrom in track.chroms() {
        let frags = track.frags(chrom);
        if frags.is_empty() {
            continue;
        }
        let counts = track.counts(chrom);
        let end = frags[frags.len() - 1].end;
        let mut t = SignalTrack::empty(chrom, 0, end);
        let mut pv: Vec<(Coord, f32)> = Vec::with_capacity(2 * frags.len());
        for (f, k) in frags.iter().zip(counts.iter()) {
            let w = *k as f32;
            pv.push((f.start, w));
            pv.push((f.end, -w));
        }
        pv.sort_by_key(|x| x.0);
        let mut depth = 0.0f32;
        let mut prev: Coord = 0;
        for &(p, d) in &pv {
            if p > prev {
                t.push(p, depth);
            }
            depth += d;
            prev = p;
        }
        bg.insert_track(chrom, t);
    }
    bg
}

/// `fc_bdg.apply_func(lambda x: x/mean_v)`.
fn apply_fold_change(bg: &BedGraph, mean: f32) -> BedGraph {
    // Same genome as the source, for the same reason.
    let mut out = BedGraph::with_genome(bg.genome(), bg.baseline());
    for (chrom, t) in bg.iter_sorted() {
        let runs = t.runs();
        let mut n = SignalTrack::empty(chrom, 0, runs.last().map_or(0, |r| r.end));
        for r in runs {
            n.push(r.end, r.value / mean);
        }
        out.insert_track(chrom, n);
    }
    out
}

/// Group a flat region list into a `RegionSet` with chromosomes in name order.
fn regions_of(genome: &Genome, v: &[(ChromId, Coord, Coord)]) -> RegionSet {
    let mut by: BTreeMap<usize, Vec<(Coord, Coord)>> = BTreeMap::new();
    for &(c, s, e) in v {
        by.entry(c.0 as usize).or_default().push((s, e));
    }
    let mut ids: Vec<ChromId> = by.keys().map(|&i| ChromId(i as u32)).collect();
    ids.sort_by(|&a, &b| genome.name(a).cmp(genome.name(b)));
    let mut rs = RegionSet::new(genome.clone());
    for c in ids {
        let mut list = by.remove(&(c.0 as usize)).unwrap_or_default();
        list.sort();
        rs.extend_chrom(c, list);
    }
    rs
}

/// Bins and extracted rows for one chromosome's regions.
///
/// `extract_value_hmmr` visits the shared chromosomes by **popping** a sorted
/// list, so the rows come out in descending name order; see the module note.
fn extract_one_chrom(
    signals: &[&SignalTrack<f32>; N_SIGNALS],
    genome: &Genome,
    chrom: ChromId,
    regions: &[(Coord, Coord)],
    binsize: Coord,
    poisson: bool,
) -> Vec<ExtractedBin> {
    let bins = make_bdg_of_bins_from_regions(
        genome,
        &BTreeMap::from([(chrom, regions.to_vec())]),
        binsize,
    );
    let rows = macs_hmmratac::extract_value_hmmr(signals, &bins);
    // `extract_signals_from_regions` floors each value at 1e-4 and, for a
    // Poisson model, truncates it to an integer count first.
    if poisson {
        rows.into_iter()
            .map(|mut r| {
                r.values = r.values.map(|v| (v.max(0.0001) as i64) as f32);
                r
            })
            .collect()
    } else {
        rows.into_iter()
            .map(|mut r| {
                r.values = r.values.map(|v| v.max(0.0001));
                r
            })
            .collect()
    }
}

/// Decode every region's bin sequence and return `(chrom, bin_end, posteriors)`
/// in upstream's emission order.
fn decode(
    model: &ModelFile,
    digested: &[BTreeMap<ChromId, SignalTrack<f32>>],
    genome: &Genome,
    regions: &[(ChromId, Coord, Coord)],
    binsize: Coord,
) -> Vec<(ChromId, ProbRow)> {
    let poisson = model.hmm_type == HmmType::Poisson;
    let hmm = if poisson {
        None
    } else {
        Some(model.gaussian())
    };
    let phmm = if poisson { Some(model.poisson()) } else { None };
    let mut out: Vec<(ChromId, ProbRow)> = Vec::new();
    // group by chromosome, then visit descending by name (the pop-order quirk)
    let mut by: BTreeMap<usize, Vec<(Coord, Coord)>> = BTreeMap::new();
    for &(c, s, e) in regions {
        by.entry(c.0 as usize).or_default().push((s, e));
    }
    let mut ids: Vec<ChromId> = by.keys().map(|&i| ChromId(i as u32)).collect();
    ids.sort_by(|&a, &b| genome.name(b).cmp(genome.name(a)));
    for c in ids {
        let regions = by.remove(&(c.0 as usize)).unwrap_or_default();
        // `extract_value_hmmr` only visits the chromosomes present in *both* the
        // signal tracks and the bin track, so a chromosome missing from any
        // digested signal contributes no rows.
        let present: Option<Vec<&SignalTrack<f32>>> =
            digested.iter().map(|signals| signals.get(&c)).collect();
        let Some(present) = present else { continue };
        let sig: [&SignalTrack<f32>; N_SIGNALS] = present.try_into().ok().unwrap();
        let rows = extract_one_chrom(&sig, genome, c, &regions, binsize, poisson);
        // Diagnostic dump of the exact bins fed to the HMM, for differencing against
        // the oracle's `cr_bins`/`cr_data`. Off unless MACS3_RS_DUMP_HMM_BINS is set.
        if let Some(dir) = std::env::var_os("MACS3_RS_DUMP_HMM_BINS") {
            use std::io::Write as _;
            let path = std::path::Path::new(&dir).join(format!(
                "bins_{}.tsv",
                String::from_utf8_lossy(genome.name(c)).replace('/', "_")
            ));
            if let Ok(mut f) = std::fs::File::create(&path) {
                for r in &rows {
                    let _ = writeln!(
                        f,
                        "{}\t{}\t{:.9e}\t{:.9e}\t{:.9e}\t{:.9e}",
                        r.pos, r.mark, r.values[0], r.values[1], r.values[2], r.values[3]
                    );
                }
            }
        }
        // split into region sequences: `mark` identifies the region
        let mut i = 0usize;
        while i < rows.len() {
            let mark = rows[i].mark;
            let mut j = i;
            while j < rows.len() && rows[j].mark == mark {
                j += 1;
            }
            let obs: Vec<Vec<f64>> = rows[i..j]
                .iter()
                .map(|r| r.values.iter().map(|v| f64::from(*v)).collect())
                .collect();
            let post = if poisson {
                phmm.as_ref().unwrap().posterior(&obs)
            } else {
                hmm.as_ref().unwrap().posterior(&obs)
            };
            for (k, p) in post.into_iter().enumerate() {
                let end = rows[i + k].pos;
                let probs = [
                    p.first().copied().unwrap_or(0.0),
                    p.get(1).copied().unwrap_or(0.0),
                    p.get(2).copied().unwrap_or(0.0),
                ];
                out.push((c, ProbRow { end, probs }));
            }
            i = j;
        }
    }
    out
}

/// `write_bedGraph`, the format `bedGraphIO.write_bedGraph` emits.
/// F204: `bedGraphIO.write_bedGraph` (`BedGraphIO.py:98-143`).
///
/// Three things were wrong here, and together they made `--save-likelihoods` and
/// `--save-digested` unreadable as bedGraph:
///
/// * **The start column was missing.** Upstream keeps a running `pre` and writes
///   `chrom, pre, pos, value` -- four fields. This port wrote only the end, producing
///   `chr1 10 0`, which is not a valid bedGraph line at all.
/// * **Wrong value format.** Upstream is unconditionally `"%.5f"`. The general `%g`
///   helper `fmt_bdg` was used instead, emitting 20 significant digits, so values like
///   `4.03e-13` did not collapse to `0.00000` the way upstream's do.
/// * **Wrong track line.** Upstream writes
///   `track type=bedGraph name=".." description=".." visibility=2 alwaysZero=on`
///   and only when `trackline=True`. `--save-likelihoods` passes `trackline=False`, so
///   those three files have no track line at all (`hmmratac_cmd.py:551`).
fn bedgraph_text(
    rows: &[BdgRow],
    genome: &Genome,
    name: &str,
    desc: &str,
    trackline: bool,
) -> String {
    let mut s = String::new();
    if trackline {
        s.push_str(&format!(
            "track type=bedGraph name=\"{name}\" description=\"{desc}\" \
             visibility=2 alwaysZero=on\n"
        ));
    }
    for &(c, s0, e0, v) in rows {
        s.push_str(&format!(
            "{}\t{}\t{}\t{:.5}\n",
            String::from_utf8_lossy(genome.name(c)),
            s0,
            e0,
            f64::from(v)
        ));
    }
    s
}

fn parse_vec4(s: Option<&str>, default: &[f32; 4], what: &str) -> Result<[f32; 4]> {
    let Some(s) = s else { return Ok(*default) };
    let mut v: Vec<f32> = Vec::new();
    for part in s.split(',') {
        v.push(part.trim().parse::<f32>().map_err(|_| {
            MacsError::InvalidParameter(format!("{what} must be four comma-separated numbers"))
        })?);
    }
    if v.len() != 4 {
        return Err(MacsError::InvalidParameter(format!(
            "{what} must be four comma-separated numbers, got {}",
            v.len()
        )));
    }
    Ok([v[0], v[1], v[2], v[3]])
}

fn need(o: &Options, dest: &str) -> Result<f32> {
    o.float(dest)
        .map(|v| v as f32)
        .ok_or_else(|| MacsError::InvalidParameter(format!("--{dest} is required")))
}

/// Run `hmmratac`.
pub fn hmmratac(o: &Options) -> Result<()> {
    let inputs = o.get_all("input_file").to_vec();
    if inputs.is_empty() {
        return Err(MacsError::InvalidParameter("-i/--input is required".into()));
    }
    let format = o.get("format").unwrap_or("BAMPE").to_uppercase();
    if !matches!(format.as_str(), "BAMPE" | "BEDPE" | "FRAG") {
        return Err(MacsError::InvalidParameter(format!(
            "wrong format `{format}`: expected BAMPE, BEDPE or FRAG"
        )));
    }
    // `opt_validate_hmmratac`: every usage error it raises before touching a file.
    let min_frag_p = need(o, "min_frag_p")?;
    if !(min_frag_p > 0.0 && min_frag_p < 1.0) {
        return Err(MacsError::InvalidParameter(
            "`--min-frag-p` should be larger than 0 and smaller than 1!".into(),
        ));
    }
    // `opt_validate_hmmratac`: `if options.hmm_binsize <= 0: error`. The flag is
    // a signed int upstream, so a negative `--binsize` is a usage error rather
    // than something to clamp.
    let binsize_raw = o.int("hmm_binsize").unwrap_or(10);
    if binsize_raw <= 0 {
        return Err(MacsError::InvalidParameter(
            "`--binsize` must be larger than 0.".into(),
        ));
    }
    let binsize = binsize_raw as Coord;
    let lower = o.float("hmm_lower").unwrap_or(10.0);
    let upper = o.float("hmm_upper").unwrap_or(20.0);
    if lower < 0.0 {
        return Err(MacsError::InvalidParameter(
            "`-l` or `--lower` should not be negative!".into(),
        ));
    }
    if upper < 0.0 {
        return Err(MacsError::InvalidParameter(
            "`-u` or `--upper` should not be negative!".into(),
        ));
    }
    if lower > upper {
        return Err(MacsError::InvalidParameter(
            "Upper limit of fold change range should be greater than lower limit!".into(),
        ));
    }
    let max_train = o.int("hmm_maxTrain").unwrap_or(1000);
    if max_train <= 0 {
        return Err(MacsError::InvalidParameter(
            "`--maxTrain` should be larger than 0!".into(),
        ));
    }
    let prescan_cutoff = o.float("prescan_cutoff").unwrap_or(1.2);
    if prescan_cutoff <= 1.0 {
        return Err(MacsError::InvalidParameter(
            "In order to use -c or --prescan-cutoff, the cutoff must be larger than 1.".into(),
        ));
    }
    let minlen = o.int("openregion_minlen").unwrap_or(100);
    if minlen < 0 {
        return Err(MacsError::InvalidParameter(
            "In order to use --minlen, the length should not be negative.".into(),
        ));
    }
    let em_means = parse_vec4(o.get("em_means"), &[50.0, 200.0, 400.0, 600.0], "means")?;
    let em_stddevs = parse_vec4(o.get("em_stddevs"), &[20.0, 20.0, 20.0, 20.0], "stddevs")?;
    if em_means.iter().any(|v| *v < 0.0) {
        return Err(MacsError::InvalidParameter(
            "`--means` should not be negative!".into(),
        ));
    }
    if em_stddevs.iter().any(|v| *v < 0.0) {
        return Err(MacsError::InvalidParameter(
            "`--stddev` should not be negative!".into(),
        ));
    }
    let hmm_type = match o
        .get("hmm_type")
        .unwrap_or("gaussian")
        .to_lowercase()
        .as_str()
    {
        "gaussian" => HmmType::Gaussian,
        "poisson" => HmmType::Poisson,
        other => {
            return Err(MacsError::InvalidParameter(format!(
                "unknown `--hmm-type` `{other}`: expected gaussian or poisson"
            )))
        }
    };
    let _ = hmm_type;
    if o.flag("modelonly") && o.get("hmm_file").is_none() {
        return Err(MacsError::InvalidParameter(SELF_TRAINING_DEVIATION.into()));
    }

    let name = o.get("name").unwrap_or("hmmratac").to_string();
    let outdir = PathBuf::from(o.get("outdir").unwrap_or("."));
    let out = Outputs::new(&outdir, &name);
    let flanking = o.int("hmm_training_flanking").unwrap_or(1000).max(0) as Coord;

    // ---- 1. read fragments ------------------------------------------------
    let barcodes: Option<HashSet<Vec<u8>>> = o.get("barcodefile").map(|p| {
        std::fs::read_to_string(p)
            .unwrap_or_default()
            .lines()
            .map(|l| l.trim_end().as_bytes().to_vec())
            .collect()
    });
    let mut petrack = load_fragments(
        &inputs,
        &format,
        barcodes.as_ref(),
        o.flag("maxcount"),
        o.int("maxcount").unwrap_or(0).max(0) as u32,
    )?;
    if o.flag("misc_remove_duplicates") {
        filter_frag_dup(&mut petrack, 1)?;
    }
    if petrack.total() == 0 {
        return Err(MacsError::InvalidParameter(
            "no fragments were read from the input".into(),
        ));
    }
    let genome = petrack.genome().clone();

    // ---- 2. fragment-length EM -------------------------------------------
    let (mut means4, mut stddevs4) = if o.flag("em_skip") {
        (em_means, em_stddevs)
    } else {
        let p = EmParams {
            min_fraglen: o.int("min_fraglen").unwrap_or(100).max(0) as Coord,
            max_fraglen: o.int("max_fraglen").unwrap_or(1000).max(0) as Coord,
            sample_percentage: o.float("sample_percentage").unwrap_or(10.0) as f32,
            seed: o.int("hmm_randomSeed").unwrap_or(10151),
            jump: o.float("em_jump").unwrap_or(0.5) as f32,
            ..Default::default()
        };
        let sampled = sample_frag_percent(
            &petrack,
            f64::from(p.sample_percentage) / 100.0,
            p.seed,
            p.min_fraglen,
            p.max_fraglen,
        )?;
        // `HMMR_EM.__init__` logs the post-filter sample size
        // (`HMMR_EM.py:168`); emitting it makes the RNG/sampling step directly
        // comparable against the oracle instead of inferring it from the means.
        eprintln!(
            "hmmratac: # Downsampled {} fragments will be used for EM training...",
            sampled.len()
        );
        let r = train_em(
            &sampled,
            [em_means[1], em_means[2], em_means[3]],
            [em_stddevs[1], em_stddevs[2], em_stddevs[3]],
            &p,
        );
        prepend_short(em_means[0], em_stddevs[0], &r)
    };
    round_to_one_decimal(&mut means4);
    round_to_one_decimal(&mut stddevs4);

    // Upstream logs these two lines (`hmmratac_cmd.py:198-199`) in a fixed-width
    // `{:>10.4g}` layout. Emitting them makes the EM result directly comparable
    // against the oracle's log without having to guess at internal state.
    {
        let fm = |v: &[f32; 4]| -> String {
            v.iter()
                .map(|x| format!("{x:>10.4}"))
                .collect::<Vec<_>>()
                .join(" ")
        };
        eprintln!(
            "hmmratac: #             means: {}",
            fm(&[means4[0], means4[1], means4[2], means4[3]])
        );
        eprintln!(
            "hmmratac: #           stddevs: {}",
            fm(&[stddevs4[0], stddevs4[1], stddevs4[2], stddevs4[3]])
        );
    }

    // ---- digested signals --------------------------------------------------
    let fl_list: Vec<Coord> = petrack.count_fraglengths().into_keys().collect();
    let weight_mapping = generate_weight_mapping(&fl_list, &means4, &stddevs4, min_frag_p);
    // `pileup_bdg_hmmr` (`PairedEndTrack.py:2078`) "explodes" the counts with
    // `np.repeat(arange(n), counts)` *before* applying the weight mapping, so a
    // row with count 7 adds the class weight seven times. That is not the same
    // as adding `7 * w` once in f32, so the list is expanded rather than scaled.
    let mut frags: BTreeMap<ChromId, Vec<(Coord, Coord)>> = BTreeMap::new();
    for c in petrack.chroms() {
        let fs = petrack.frags(c);
        let ks = petrack.counts(c);
        let mut v: Vec<(Coord, Coord)> = Vec::new();
        for (f, k) in fs.iter().zip(ks.iter()) {
            for _ in 0..*k {
                v.push((f.start, f.end));
            }
        }
        if !v.is_empty() {
            frags.insert(c, v);
        }
    }
    let digested = pileup_bdg_hmmratac(&genome, &frags, &weight_mapping);

    // ---- fold-change track -------------------------------------------------
    let min_template = petrack.average_template_length().max(0.0) as Coord;
    let fc_source = if o.flag("pileup_short") {
        // F226: the tracks are keyed by `petrack`'s ChromIds, so the BedGraph must
        // resolve names through the same dictionary -- `BedGraph::new` starts empty.
        let mut bg = BedGraph::with_genome(&genome, 0.0);
        for (c, t) in digested[0].iter() {
            bg.insert_track(*c, t.clone());
        }
        bg
    } else {
        pileup_fragments(&petrack)
    };
    let mean_v = peakcall::summary(&fc_source).mean;
    if !mean_v.is_finite() || mean_v == 0.0 {
        return Err(MacsError::InvalidParameter(
            "the fragment pileup has zero mean; nothing to normalise against".into(),
        ));
    }
    let fc_bdg = apply_fold_change(&fc_source, mean_v);

    if o.flag("cutoff_analysis_only") {
        std::fs::create_dir_all(&outdir)?;
        let t = peakcall::cutoff_analysis(
            &fc_bdg,
            flanking,
            min_template,
            o.int("cutoff_analysis_steps").unwrap_or(100),
            0.0,
            o.float("cutoff_analysis_max").unwrap_or(100.0) as f32,
        );
        std::fs::write(out.cutoff_analysis(), t)?;
        return Ok(());
    }

    // ---- 4. the model ------------------------------------------------------
    let model = match o.get("hmm_file") {
        Some(p) => ModelFile::load(Path::new(p))?,
        None => return Err(MacsError::InvalidParameter(SELF_TRAINING_DEVIATION.into())),
    };
    let binsize = model.hmm_binsize.max(binsize);
    // ---- 5. decode ---------------------------------------------------------
    let mut candidates: Vec<(ChromId, Coord, Coord)> = Vec::new();
    for (c, peaks) in peakcall::call_peaks(&fc_bdg, prescan_cutoff as f32, min_template, flanking) {
        for pk in peaks {
            candidates.push((c, pk.start, pk.end));
        }
    }
    eprintln!("hmmratac: #5  Total candidate peaks : {}", candidates.len());
    eprintln!("hmmratac: #  We expand the candidate regions with {flanking} and merge overlap");
    let mut cr = regions_of(&genome, &candidates);
    cr.expand(flanking);
    cr.merge_overlap();
    eprintln!(
        "hmmratac: #   after expanding and merging, we have {} candidate regions",
        cr.total()
    );
    // Diagnostic: the exact candidate regions handed to the decoder, so they can be
    // differenced against the oracle's `Regions` object. Off unless
    // MACS3_RS_DUMP_HMM_REGIONS is set.
    if let Some(path) = std::env::var_os("MACS3_RS_DUMP_HMM_REGIONS") {
        use std::io::Write as _;
        if let Ok(mut f) = std::fs::File::create(std::path::Path::new(&path)) {
            for c in cr.chroms_sorted() {
                for (s, e) in cr.chrom(c) {
                    let _ = writeln!(f, "{}\t{s}\t{e}", String::from_utf8_lossy(genome.name(c)));
                }
            }
        }
    }

    let decoding_steps = o.int("decoding_steps").unwrap_or(5000).max(1) as usize;
    let mut decoded: Vec<(ChromId, ProbRow)> = Vec::new();
    while !cr.is_empty() {
        let batch = cr.pop(decoding_steps).expect("non-empty");
        let flat = batch.flatten();
        decoded.extend(decode(&model, &digested, &genome, &flat, binsize));
    }

    // ---- 3. training regions (informational outputs only) ------------------
    // Upstream writes these before training; with `--model` step 3 is skipped
    // entirely. They are reproduced here so `--save-digested`/`--save-train`
    // users still see the same files, but they cannot change the decode.
    if o.get("hmm_file").is_none()
        && o.get("hmm_training_regions").is_none()
        && o.flag("save_train")
    {
        let mut kept: Vec<(ChromId, Coord, Coord)> = Vec::new();
        for (c, peaks) in peakcall::call_peaks(&fc_bdg, lower as f32, min_template, flanking) {
            for pk in peaks {
                if pk.score >= lower as f32 && pk.score < upper as f32 {
                    kept.push((c, pk.start, pk.end));
                }
            }
        }
        if kept.is_empty() {
            return Err(MacsError::InvalidParameter(
                "Not enough training regions! Please adjust the lower or upper cutoff.".into(),
            ));
        }
        kept.truncate(max_train.max(0) as usize);
        let mut rs = regions_of(&genome, &kept);
        rs.expand(flanking);
        rs.merge_overlap();
        std::fs::create_dir_all(&outdir)?;
        std::fs::write(out.training_regions(), rs.to_bed_string())?;
    }

    // ---- 6. outputs --------------------------------------------------------
    std::fs::create_dir_all(&outdir)?;
    if o.flag("save_digested") {
        for (i, k) in ["short", "mono", "di", "tri"].iter().enumerate() {
            std::fs::write(
                out.digested(k),
                // `--save-digested` uses the default `trackline=True`.
                bedgraph_text(&rows_of(&digested[i]), &genome, k, k, true),
            )?;
        }
    }

    let decoded: Vec<(ChromId, ProbRow)> = decoded;
    if o.flag("save_likelihoods") {
        let (open, nuc, bg) = likelihood_rows(&decoded, binsize, &model);
        for (i, (label, desc, src)) in [
            ("open", "Likelihoods of being Open States", &open),
            ("nuc", "Likelihoods of being Nucleosomal States", &nuc),
            ("bg", "Likelihoods of being Background States", &bg),
        ]
        .into_iter()
        .enumerate()
        {
            let _ = i;
            std::fs::write(
                out.state_bdg(label),
                // `hmmratac_cmd.py:551-553`: `write_bedGraph(..., trackline=False)`.
                bedgraph_text(src, &genome, label, desc, false),
            )?;
        }
    }

    let states_path = generate_states_path_chromed(&decoded, binsize, &model);
    if o.flag("save_states") {
        std::fs::write(out.states_bed(), states_bed_text(&states_path, &genome))?;
    }

    // `accessible_regions[:-1]` upstream: the final triple is dropped.
    let mut open = accessible_regions_chromed(&states_path, minlen as Coord);
    if !open.is_empty() {
        open.pop();
    }
    std::fs::write(
        out.accessible(),
        narrowpeak_text(&peakcall::refine_peaks(&fc_bdg, &open), &genome, &name),
    )?;
    Ok(())
}

/// The declared deviation, reported identically wherever it applies.
const SELF_TRAINING_DEVIATION: &str =
    "HMMRATAC self-training (hmmlearn Baum-Welch) is this port's one declared deviation and is \
     not reproduced; supply a trained model with --model. Inference against a supplied model is \
     exact.";

fn rows_of(t: &BTreeMap<ChromId, SignalTrack<f32>>) -> Vec<BdgRow> {
    let mut out = Vec::new();
    for (c, t) in t {
        let mut prev = 0;
        for r in t.runs() {
            out.push((*c, prev, r.end, r.value));
            prev = r.end;
        }
    }
    out
}

/// One emitted bedGraph record: `(chrom, start, end, value)`.
type BdgRow = (ChromId, Coord, Coord, f32);
/// The `open` / `nuc` / `bg` likelihood tracks, in that order.
type LikelihoodTracks = (Vec<BdgRow>, Vec<BdgRow>, Vec<BdgRow>);

/// `save_proba_to_bedGraph` (`hmmratac_cmd.py:517`): three bedGraphs, one per
/// state.
///
/// Gaps are filled with `bg = 1.0` and `open = nuc = 0.0` -- the value in force
/// outside any decoded bin. `add_loc` coalesces an equal-valued predecessor, so
/// a run of consecutive bins with identical posteriors collapses to one record.
fn likelihood_rows(
    rows: &[(ChromId, ProbRow)],
    binsize: Coord,
    model: &ModelFile,
) -> LikelihoodTracks {
    let mut open: Vec<BdgRow> = Vec::new();
    let mut nuc: Vec<BdgRow> = Vec::new();
    let mut bg: Vec<BdgRow> = Vec::new();
    let mut prev_chrom: Option<ChromId> = None;
    let mut prev_end: Coord = 0;
    for (chrom, r) in rows {
        let start = r.end.saturating_sub(binsize);
        if prev_chrom != Some(*chrom) {
            if start > 0 {
                for v in [&mut open, &mut nuc] {
                    v.push((*chrom, 0, start, 0.0));
                }
                bg.push((*chrom, 0, start, 1.0));
            }
            prev_chrom = Some(*chrom);
        } else if prev_end < start {
            for v in [&mut open, &mut nuc] {
                v.push((*chrom, prev_end, start, 0.0));
            }
            bg.push((*chrom, prev_end, start, 1.0));
        }
        open.push((*chrom, start, r.end, r.probs[model.i_open_region] as f32));
        nuc.push((
            *chrom,
            start,
            r.end,
            r.probs[model.i_nucleosomal_region] as f32,
        ));
        bg.push((
            *chrom,
            start,
            r.end,
            r.probs[model.i_background_region] as f32,
        ));
        prev_end = r.end;
    }
    (open, nuc, bg)
}

/// `save_states_bed`: `chrom\tstart\tend\tlabel`, skipping `bg`.
///
/// The `chr` byte slice is written with `%s` after `.decode()`, so a non-UTF-8
/// chromosome name would raise upstream; here it is lossily replaced rather than
/// aborting the run.
fn states_bed_text(path: &[macs_hmmratac::StateRun], genome: &Genome) -> String {
    let mut out = String::new();
    for r in path {
        if r.label != "bg" {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\n",
                String::from_utf8_lossy(genome.name(r.chrom)),
                r.start,
                r.end,
                r.label
            ));
        }
    }
    out
}

/// `PeakIO.write_to_narrowPeak` over the refined accessible regions.
///
/// The format is ENCODE narrowPeak (BED6+4), and **`PeakIO` writes columns 7-9 as
/// three `%.6g` values, not coordinates**:
///
/// ```text
/// "%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\t%d\n"
/// ```
///
/// * column 5 -- `int(10 * peak[score_column])`: an **integer**, truncating toward zero.
/// * column 7 -- `peak['fc']`, fold enrichment.
/// * column 8 -- `peak['pscore']`, `-log10(p)`.
/// * column 9 -- `peak['qscore']`, `-log10(q)`.
/// * column 10 -- the **summit offset from the peak start**, or `-1` when unset.
///
/// F227: this port wrote `start`, `end` and the absolute summit into columns 7-9 and
/// the fold change into column 5. Column 10 was already right, which is why the summit
/// *offsets* looked plausible while the rest of the line was structurally wrong -- and
/// `_accessible_regions.narrowPeak` is covered by the "byte-identical
/// `*_peaks.narrowPeak`" criterion.
///
/// Peaks are emitted grouped by equal `end` within a chromosome, and a group with more
/// than one member (only reachable via `--call-summits`) gets the `a`/`b`/`c` suffixes.
/// `sorted(chrs)` is Python's byte-order sort of the names.
///
/// `name` is the dataset label, and is deliberately **unused**: `save_accessible_regions`
/// calls `write_to_narrowPeak(fhd)` with its defaults, so the peak names are
/// `MACS_peak_<n>` regardless of `--name`. The parameter is kept so the call site reads
/// the way the signature does.
fn narrowpeak_text(peaks: &[(ChromId, peakcall::BdgPeak)], genome: &Genome, _name: &str) -> String {
    let mut by_chrom: BTreeMap<Vec<u8>, Vec<(u64, usize, &peakcall::BdgPeak)>> = BTreeMap::new();
    for (i, (c, pk)) in peaks.iter().enumerate() {
        by_chrom
            .entry(genome.name(*c).to_vec())
            .or_default()
            .push((pk.end, i, pk));
    }

    let mut s = String::new();
    let mut n_peak = 0usize;
    // `sorted(chrs)` sorts the decoded names, which is byte order.
    let mut chroms: Vec<Vec<u8>> = by_chrom.keys().cloned().collect();
    chroms.sort();
    for chrom in chroms {
        let mut rows = by_chrom.remove(&chrom).unwrap_or_default();
        rows.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut i = 0;
        while i < rows.len() {
            let end = rows[i].0;
            let mut j = i;
            while j < rows.len() && rows[j].0 == end {
                j += 1;
            }
            let group = &rows[i..j];
            n_peak += 1;
            // `save_accessible_regions` calls `write_to_narrowPeak(fhd)` with **no**
            // `name_prefix`/`name`, so the prefix is the default `b"MACS_peak_"` and
            // the dataset name is the default `b"MACS"` -- the user's `--name` is *not*
            // interpolated here. Getting that wrong is visible on every line.
            if group.len() > 1 {
                for (k, (_, _, pk)) in group.iter().enumerate() {
                    s.push_str(&narrowpeak_line(
                        &chrom,
                        pk,
                        &format!("MACS_peak_{}{}", n_peak, subpeak_letters(k)),
                    ));
                }
            } else {
                let (_, _, pk) = group[0];
                s.push_str(&narrowpeak_line(
                    &chrom,
                    pk,
                    &format!("MACS_peak_{}", n_peak),
                ));
            }
            i = j;
        }
    }
    s
}

/// `subpeak_letters`: `a`, `b`, ... `z`, `aa`, `ab`, ... for grouped sub-peaks.
fn subpeak_letters(mut i: usize) -> String {
    let mut out = Vec::new();
    loop {
        out.push(b'a' + (i % 26) as u8);
        if i < 26 {
            break;
        }
        i = i / 26 - 1;
    }
    out.reverse();
    String::from_utf8_lossy(&out).into_owned()
}

/// One narrowPeak record, in `PeakIO.write_to_narrowPeak`'s exact format.
///
/// Columns 7-9 are the peak's `fc`, `pscore` and `qscore` attributes with `%.6g`.
/// `save_accessible_regions` builds its `PeakIO` with `openpeak.add(chromosome=...,
/// start=..., end=...)` and never sets any of those three, so all three are the
/// default **0** -- *not* the fold change the digest tracks carry. The fold change is
/// only in column 5, as `int(10 * score)`.
///
/// `summit` is `u64` here so upstream's `-1` sentinel for an unset summit cannot be
/// represented; `refine_peaks` always sets it on this path.
fn narrowpeak_line(chrom: &[u8], pk: &peakcall::BdgPeak, name: &str) -> String {
    let summit_off = pk.summit.saturating_sub(pk.start) as i64;
    format!(
        "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\t{}\n",
        String::from_utf8_lossy(chrom),
        pk.start,
        pk.end,
        name,
        (10.0f32 * pk.score) as i64,
        macs_io::peakout::format_g(0.0, 6),
        macs_io::peakout::format_g(0.0, 6),
        macs_io::peakout::format_g(0.0, 6),
        summit_off,
    )
}
