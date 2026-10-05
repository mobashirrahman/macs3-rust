//! BedGraph track model and arithmetic -- a port of
//! `MACS3/Signal/BedGraph.py` (`bedGraphTrackI`).
//!
//! A bedGraph is a per-chromosome list of `(end, value)` pairs: `value[i]` is
//! the signal over `[end[i-1], end[i])`, with the first run starting at the
//! chromosome's left edge (0, or the `baseline` when the first record starts
//! later). That is exactly [`macs_rle::SignalTrack`]'s run representation, so a
//! `BedGraph` is a genome plus one `SignalTrack<f32>` per chromosome and the
//! upstream pointer walks port directly onto the run slices.
//!
//! Upstream fidelity notes:
//!
//! * [`BedGraph::add_loc`] reproduces `add_loc`: it clamps a negative start to
//!   0, drops `end <= 0`, back-fills `[0, start)` with the baseline for a new
//!   chromosome, and **coalesces an equal-valued predecessor** by extending its
//!   end. Like upstream, it assumes the input is sorted and non-overlapping --
//!   `read_bedGraph` does not sort, so an unsorted file is reproduced faithfully
//!   (and incorrectly) rather than silently repaired.
//! * [`BedGraph::overlie`] is `overlie`: a pointer walk that emits the union of
//!   all breakpoints and applies a combining function to the values covering
//!   each interval. It stops as soon as *any* input is exhausted, exactly like
//!   upstream's `StopIteration`.
//! * [`BedGraph::extract_value`] is `extract_value`: for every interval in which
//!   the *second* track's value is `> 0`, report the *first* track's value.

pub mod peakcall;

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use macs_core::genome::ChromId;
use macs_core::{Coord, Genome, MacsError, Result};
use macs_rle::SignalTrack;

/// `log10(e)`, the constant upstream's `fisher_func` divides by
/// (`BedGraph.py:58-61`).
use std::f64::consts::LOG10_E;

/// The combining functions `overlie` accepts.
///
/// The `subtract` and `divide` variants index the value tuple exactly as
/// upstream does: `subtract` is `v[1] - v[0]` and `divide` is `v[1] / v[2]`,
/// where `v[0]` is the receiving track and `v[1..]` are the overlaid ones.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    Max,
    Sum,
    Subtract,
    Product,
    Divide,
    Mean,
    Fisher,
}

impl Op {
    /// Parse the name used by `bdgopt --transform` / `cmbreps --method`.
    pub fn parse(name: &[u8]) -> Option<Self> {
        match name {
            b"max" => Some(Self::Max),
            b"sum" => Some(Self::Sum),
            b"subtract" => Some(Self::Subtract),
            b"product" => Some(Self::Product),
            b"divide" => Some(Self::Divide),
            b"mean" => Some(Self::Mean),
            b"fisher" => Some(Self::Fisher),
            _ => None,
        }
    }

    /// Combine the values covering one interval.
    ///
    /// `vals[0]` is the receiving track, `vals[1..]` the overlaid ones. The
    /// `subtract`/`divide` index patterns are upstream's and are deliberately
    /// not "corrected" to `v[0] - v[1]`.
    pub fn apply(self, vals: &[f32]) -> f32 {
        // `sum`, `mean` and `product` accumulate in `f64`, not `f32`.
        //
        // Upstream's `mean_func` is `sum(x)/len(x)`, `product_func` is `math.prod(x)`,
        // and the `sum` combiner is `sum(x)`, all over a tuple of Python floats. So
        // the reduction runs in f64, left to right, and narrows to f32 only when the
        // result is stored into the f32 track. Summing in f32 instead loses a unit in
        // the last place on long enough inputs -- observed as `cmbreps -m mean`
        // reporting 28.74796 where upstream reports 28.74795 on 1200 of 420948 rows
        // of a real `callpeak --bdg` track.
        match self {
            Self::Max => vals.iter().copied().fold(f32::MIN, f32::max),
            Self::Sum => vals.iter().map(|&v| f64::from(v)).sum::<f64>() as f32,
            Self::Subtract => vals[1] - vals[0],
            Self::Product => vals.iter().map(|&v| f64::from(v)).product::<f64>() as f32,
            Self::Divide => {
                // upstream indexes x[1]/x[2] and *raises* IndexError with only
                // two overlaid values, so `divide` is not reachable for the
                // 2-track case in practice. Return 0.0 instead of panicking,
                // per the zero-panic criterion.
                if vals.len() < 3 {
                    0.0
                } else {
                    vals[1] / vals[2]
                }
            }
            Self::Mean => {
                (vals.iter().map(|&v| f64::from(v)).sum::<f64>() / vals.len() as f64) as f32
            }
            Self::Fisher => {
                // combine -log10 p-values with Fisher's method:
                // chisq_logp_e(2*sum(x)/log10(e), 2*len(x), log10=True)
                let s: f64 = vals.iter().map(|&v| f64::from(v)).sum();
                let df = 2 * vals.len() as u32;
                let x = 2.0 * s / LOG10_E;
                // `df` is even and >= 2 for any non-empty tuple, so this cannot
                // fail; fall back to 0 rather than propagate an impossible error.
                macs_stats::chisq_logp_e(x, df, true).unwrap_or(0.0) as f32
            }
        }
    }
}

/// A parsed bedGraph: one run-length track per chromosome.
#[derive(Debug, Clone, Default)]
pub struct BedGraph {
    genome: Genome,
    /// Keyed by `ChromId`; a `BTreeMap` keeps iteration in id order, but the
    /// operations that care about output order sort by *name* explicitly.
    tracks: BTreeMap<ChromId, SignalTrack<f32>>,
    baseline: f32,
}

impl BedGraph {
    /// An empty track set with the given baseline value.
    pub fn new(baseline: f32) -> Self {
        Self {
            genome: Genome::new(),
            tracks: BTreeMap::new(),
            baseline,
        }
    }

    /// A `BedGraph` that resolves chromosome names through `genome`.
    ///
    /// F226: [`BedGraph::new`] starts with an **empty** `Genome`, so inserting tracks
    /// whose `ChromId`s came from a different dictionary leaves `chroms_sorted` to index
    /// an empty `names` vector -- a panic, not a wrong answer. That is what
    /// `hmmratac --save-digested` did on a valid BEDPE run: `BedGraph::new(0.0)` plus
    /// tracks keyed by `petrack`'s ids, and every digest write aborted in
    /// `Genome::name`.
    ///
    /// Ids are positional, so the genome must be the *same* dictionary the tracks were
    /// keyed with. Passing a clone of it is enough; interning order is identity and
    /// `Genome::clone` reproduces ids exactly.
    pub fn with_genome(genome: &Genome, baseline: f32) -> Self {
        Self {
            genome: genome.clone(),
            tracks: BTreeMap::new(),
            baseline,
        }
    }

    /// The baseline value used to back-fill gaps before a chromosome's first record.
    pub fn baseline(&self) -> f32 {
        self.baseline
    }

    /// The genome, for resolving chromosome names to ids.
    pub fn genome(&self) -> &Genome {
        &self.genome
    }

    /// Chromosome ids present, in byte-wise name order (upstream's
    /// `sorted(common_chr)`).
    pub fn chroms_sorted(&self) -> Vec<ChromId> {
        let mut ids: Vec<ChromId> = self.tracks.keys().copied().collect();
        ids.sort_by(|&a, &b| self.genome.name(a).cmp(self.genome.name(b)));
        ids
    }

    /// The track for a chromosome, if present.
    pub fn track(&self, chrom: ChromId) -> Option<&SignalTrack<f32>> {
        self.tracks.get(&chrom)
    }

    /// All per-chromosome tracks, keyed by `ChromId`.
    ///
    /// `bdgdiff` hands these to `macs-score`'s `TwoScores::build`, which does
    /// its own four-way merge over the raw runs.
    pub fn tracks(&self) -> &BTreeMap<ChromId, SignalTrack<f32>> {
        &self.tracks
    }

    /// `make_ScoreTrackII_for_macs` (`BedGraph.py:1170-1259`): merge treatment
    /// and control into one interval table carrying both values.
    ///
    /// This is `bdgcmp`'s first step. The walk is the same `p1 < p2` / `p1 == p2`
    /// pointer advance as upstream, and -- as there -- it stops at the end of
    /// *either* track, so the tail of the longer one is dropped. Chromosomes are
    /// restricted to the intersection of the two inputs.
    pub fn make_score_track_for_macs(
        &self,
        other: &BedGraph,
        depth1: f32,
        depth2: f32,
    ) -> macs_score::ScoreTrack2 {
        let mut st = macs_score::ScoreTrack2::new(self.genome.clone(), depth1, depth2);
        // shared chromosome names, in name order
        let mut common: Vec<ChromId> = self
            .tracks
            .keys()
            .copied()
            .filter(|c| other.tracks.contains_key(c))
            .collect();
        common.sort_by(|&a, &b| self.genome.name(a).cmp(self.genome.name(b)));
        for chrom in common {
            let r1 = self.tracks[&chrom].runs();
            let r2 = other.tracks[&chrom].runs();
            let (mut i1, mut i2) = (0usize, 0usize);
            while let (Some(p1), Some(p2)) = (r1.get(i1).map(|x| x.end), r2.get(i2).map(|x| x.end))
            {
                // `p1 < p2` and `p1 == p2` emit at `p1`; `p2 < p1` emits at `p2`
                st.add(chrom, p1.min(p2), r1[i1].value, r2[i2].value);
                if p1 <= p2 {
                    i1 += 1;
                }
                if p2 <= p1 {
                    i2 += 1;
                }
            }
        }
        st
    }

    /// Iterate every `(chrom, track)` in name order.
    pub fn iter_sorted(&self) -> impl Iterator<Item = (ChromId, &SignalTrack<f32>)> {
        self.chroms_sorted()
            .into_iter()
            .filter_map(move |c| self.tracks.get(&c).map(|t| (c, t)))
    }

    /// Append `[startpos, endpos)` with `value` on `chromosome`.
    ///
    /// A port of `add_loc`: negative starts clamp to 0, `end <= 0` is dropped,
    /// a new chromosome back-fills `[0, start)` with the baseline, and an
    /// equal-valued predecessor is extended rather than duplicated.
    pub fn add_loc(&mut self, chromosome: &[u8], startpos: Coord, endpos: Coord, value: f32) {
        // upstream: `if endpos <= 0: return` and `if startpos < 0: startpos = 0`.
        // `Coord` is u64, so both reduce to `endpos == 0` and an unchanged start.
        if endpos == 0 {
            return;
        }
        let chrom = self.genome.intern(chromosome);
        let t = self
            .tracks
            .entry(chrom)
            .or_insert_with(|| SignalTrack::empty(chrom, 0, 0));
        if startpos > 0 && t.runs().is_empty() {
            // first record for this chromosome: baseline up to `startpos`
            t.push(startpos, self.baseline);
        }
        t.push(endpos, value);
    }

    /// Attach a prebuilt run-length track, taking ownership.
    ///
    /// `add_loc` cannot express a track whose leading value is not the baseline
    /// and whose first run boundary is not `0` -- the digested pileups and the
    /// fold-change track both start at their first fragment endpoint -- so this
    /// inserts the track directly.
    pub fn insert_track(&mut self, chrom: ChromId, track: SignalTrack<f32>) {
        // The caller must have interned the name already; re-interning is a
        // no-op and would only cost an allocation.
        self.tracks.insert(chrom, track);
    }

    /// Read a bedGraph file, reproducing `BedGraphIO.read_bedGraph`.
    ///
    /// Lines beginning with `track`, `#` or `browse` are skipped; the rest are
    /// `chrom start end value`. Records are appended with [`Self::add_loc`], so
    /// like upstream the file is assumed sorted and non-overlapping.
    pub fn read(path: &Path, baseline: f32) -> Result<Self> {
        let mut bg = Self::new(baseline);
        let file = std::fs::File::open(path).map_err(MacsError::Io)?;
        let reader = BufReader::new(file);
        for (i_line, line) in reader.lines().enumerate() {
            let lineno = (i_line + 1) as u64;
            let line = line.map_err(MacsError::Io)?;
            let t = line.trim_start();
            if t.starts_with("track")
                || t.starts_with('#')
                || t.starts_with("browse")
                || t.is_empty()
            {
                continue;
            }
            let f: Vec<&str> = t.split_whitespace().collect();
            if f.len() < 4 {
                return Err(MacsError::InvalidFormat {
                    format: "bedGraph",
                    line: lineno,
                    detail: format!("expected 4 columns, got {}", f.len()),
                });
            }
            let start = f[1].parse::<Coord>().map_err(|_| MacsError::ParseInt {
                format: "bedGraph",
                line: lineno,
                value: f[1].to_string(),
                expected: "start",
            })?;
            let end = f[2].parse::<Coord>().map_err(|_| MacsError::ParseInt {
                format: "bedGraph",
                line: lineno,
                value: f[2].to_string(),
                expected: "end",
            })?;
            let value = f[3].parse::<f32>().map_err(|_| MacsError::ParseInt {
                format: "bedGraph",
                line: lineno,
                value: f[3].to_string(),
                expected: "float",
            })?;
            bg.add_loc(f[0].as_bytes(), start, end, value);
        }
        Ok(bg)
    }

    /// `overlie`: combine `self` with `others` using `op`, a port of the
    /// upstream pointer walk.
    ///
    /// Breakpoints from every input are unioned; for each interval the values
    /// covering it are combined by `op`. The walk halts as soon as any input is
    /// exhausted (`StopIteration` upstream), so the result covers only the
    /// common prefix. Chromosomes are processed in name order and only those
    /// present in *every* input are emitted.
    pub fn overlie(&self, others: &[&BedGraph], op: Op) -> Self {
        assert!(!others.is_empty(), "overlie needs at least one other track");
        let mut out = BedGraph::new(self.baseline);
        // `ChromId` is genome-local, so the common-chromosome set is computed by
        // *name*: intersect the name sets, then resolve each name to an id in
        // every participating track.
        let mut common_names: Vec<Vec<u8>> = self
            .chroms_sorted()
            .iter()
            .map(|&c| self.genome.name(c).to_vec())
            .collect();
        for o in others {
            let on: Vec<Vec<u8>> = o
                .chroms_sorted()
                .iter()
                .map(|&c| o.genome.name(c).to_vec())
                .collect();
            common_names.retain(|n| on.contains(n));
        }
        for name in common_names {
            let self_chrom = self.genome.get(&name).expect("common chrom in self");
            let mut tracks: Vec<&SignalTrack<f32>> = Vec::new();
            tracks.push(self.tracks.get(&self_chrom).expect("common chrom in self"));
            for o in others {
                let oc = o.genome.get(&name).expect("common chrom in other");
                tracks.push(o.tracks.get(&oc).expect("common chrom in other"));
            }
            let merged = overlie_chrom(&tracks, op);
            // intern under the *output* genome and key the track by that id
            let out_chrom = out.genome.intern(&name);
            let mut t = SignalTrack::empty(out_chrom, 0, 0);
            for (end, value) in merged {
                t.push(end, value);
            }
            out.tracks.insert(out_chrom, t);
        }
        out
    }

    /// `extract_value`: for each interval where `other`'s value is `> 0`,
    /// report `self`'s value over that interval as `(start, end, value)`.
    ///
    /// The walk stops when either track is exhausted, matching upstream's
    /// `StopIteration`.
    pub fn extract_value(
        &self,
        chrom: ChromId,
        other: &SignalTrack<f32>,
    ) -> Vec<(Coord, Coord, f32)> {
        let Some(t1) = self.tracks.get(&chrom) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        let r1 = t1.runs();
        let r2 = other.runs();
        let (mut i1, mut i2) = (0usize, 0usize);
        let mut pre_p: Coord = 0;
        while i1 < r1.len() && i2 < r2.len() {
            let (p1, v1) = (r1[i1].end, r1[i1].value);
            let (p2, v2) = (r2[i2].end, r2[i2].value);
            let end = p1.min(p2);
            if v2 > 0.0 && end > pre_p {
                out.push((pre_p, end, v1));
            }
            pre_p = end;
            if p1 <= p2 {
                i1 += 1;
            }
            if p2 <= p1 {
                i2 += 1;
            }
        }
        out
    }

    /// `apply_func`: map every run's value through `f`, coalescing runs that
    /// become equal (the representation is a run-length track, so this is the
    /// natural, lossless form of upstream's per-value map).
    /// `p2q`: convert `-log10(p)` scores to `-log10(q)` in place, a port of
    /// `bedGraphTrackI.p2q` (`BedGraph.py:904-984`).
    ///
    /// It is a base-weighted AFDR walk over the whole track, but -- note the
    /// upstream quirk -- the rank `k` is **never incremented** inside the loop
    /// (`nhcal` is, but it is unused). So `log10(k)` is always 0 and the mapping
    /// reduces to a monotonic clamp of `v - log10(N)`, walked over the distinct
    /// scores in descending order. That is deliberately reproduced rather than
    /// "fixed" to a textbook BH, which would diverge from the oracle.
    pub fn p2q(&self) -> Self {
        // base-weighted frequency of each score, accumulated per chromosome in
        // name order (the histogram's order is part of the contract).
        let mut stat: std::collections::BTreeMap<u32, u64> = std::collections::BTreeMap::new();
        let mut total: u64 = 0;
        for (chrom, t) in self.iter_sorted() {
            let _ = chrom;
            let mut pre: Coord = 0;
            for r in t.runs() {
                *stat.entry(r.value.to_bits()).or_insert(0) += u64::from(r.end.saturating_sub(pre));
                total += u64::from(r.end.saturating_sub(pre));
                pre = r.end;
            }
        }
        if total == 0 {
            return self.clone();
        }
        let f = -(total as f64).log10();
        let mut table: std::collections::HashMap<u32, f32> = std::collections::HashMap::new();
        let mut pre_q = f32::MAX;
        // descending distinct scores
        let mut uniq: Vec<f32> = stat.keys().map(|&b| f32::from_bits(b)).collect();
        uniq.sort_by(|a, b| b.partial_cmp(a).expect("NaN score"));
        for v in uniq {
            let q = v as f64 + (0.0 + f); // log10(k)=0, k==1
            let q = q.max(0.0).min(pre_q as f64) as f32;
            table.insert(v.to_bits(), q);
            pre_q = q;
        }
        let mut out = BedGraph::new(self.baseline);
        out.genome = self.genome.clone();
        for (chrom, t) in &self.tracks {
            let mut nt = SignalTrack::empty(*chrom, t.start(), t.end());
            for r in t.runs() {
                let q = table.get(&r.value.to_bits()).copied().unwrap_or(0.0);
                nt.push(r.end, q);
            }
            out.tracks.insert(*chrom, nt);
        }
        out.merge_regions();
        out
    }

    /// `merge_regions`: coalesce adjacent runs carrying the same value within
    /// each chromosome (`BedGraph.py:269-309`). The run-length representation
    /// coalesces equal neighbours on push already, so this only needs to sort
    /// chromosomes by name and normalise the map.
    pub fn merge_regions(&self) -> Self {
        let mut out = BedGraph::new(self.baseline);
        out.genome = self.genome.clone();
        for (chrom, t) in &self.tracks {
            let mut nt = SignalTrack::empty(*chrom, t.start(), t.end());
            for r in t.runs() {
                nt.push(r.end, r.value);
            }
            out.tracks.insert(*chrom, nt);
        }
        out
    }

    /// Map every value, preserving interval boundaries.
    ///
    /// `apply_func` (`BedGraph.py:888-901`) rewrites each run's value in place, and
    /// its docstring is explicit about the consequence: *"Two adjacent regions with
    /// same value after applying func will not be merged."*
    ///
    /// So this must use the non-coalescing push. Using the ordinary `push` merged
    /// every pair of neighbours that mapped to the same value, which for
    /// `-m max`/`-m min` is most of the file: `bdgopt -i pileup.bdg -m max -p 100`
    /// emitted 74 rows where upstream emits one per input row (91983 on a real
    /// `callpeak --bdg` track), because long stretches of sub-threshold coverage all
    /// clamp to the same constant.
    pub fn apply_func<F: Fn(f32) -> f32>(&self, f: F) -> Self {
        let mut out = BedGraph::new(self.baseline);
        out.genome = self.genome.clone();
        for (chrom, t) in &self.tracks {
            let mut nt = SignalTrack::empty(*chrom, t.start(), t.end());
            for r in t.runs() {
                nt.push_exact(r.end, f(r.value));
            }
            out.tracks.insert(*chrom, nt);
        }
        out
    }

    /// Serialise to a bedGraph file, matching `write_bedGraph`'s body: a
    /// `track` line (when `trackline`), then `chrom start end value` per run.
    /// Serialise to a bedGraph file, matching `write_bedGraph`
    /// (`BedGraphIO.py:97-143`): an optional UCSC `track` line, then one
    /// `chrom start end value` row per run, chromosomes in name order and each
    /// row's value formatted **`%.5f`** (fixed five decimals -- *not* the
    /// `%.6g` the peak writers use).
    /// Serialise to a bedGraph file, streaming.
    ///
    /// A whole-file `String` was held here, which on real data is hundreds of MB
    /// (the `callpeak --bdg` treatment pileup of a 5 M-read run is ~290 MB and the
    /// control lambda ~510 MB) and dominated every `bdg*` command's peak. Measured
    /// on those files: `cmbreps -m max` 1096 MB against upstream's 401 MB and
    /// `bdgcmp -m ppois` 1383 MB against 626 MB -- both *worse* than upstream purely
    /// because the output was buffered. The bytes are unchanged: same track line,
    /// same row order, same `%.5f`.
    pub fn write(&self, path: &Path, trackline: bool, name: &str, description: &str) -> Result<()> {
        use std::io::Write as _;
        let file = std::fs::File::create(path).map_err(MacsError::Io)?;
        let mut out = std::io::BufWriter::new(file);
        if trackline {
            writeln!(
                out,
                "track type=bedGraph name=\"{}\" description=\"{}\" visibility=2 alwaysZero=on",
                name.replace('"', "\\\""),
                description.replace('"', "\\\"")
            )
            .map_err(MacsError::Io)?;
        }
        for (chrom, t) in self.iter_sorted() {
            let cname = String::from_utf8_lossy(self.genome.name(chrom));
            let mut pre = 0u32;
            for r in t.runs() {
                writeln!(out, "{}\t{}\t{}\t{:.5}", cname, pre, r.end, r.value)
                    .map_err(MacsError::Io)?;
                pre = r.end;
            }
        }
        out.flush().map_err(MacsError::Io)?;
        Ok(())
    }
}

/// `overlie` for one chromosome: union the breakpoints of `tracks` and apply
/// `op` to the values covering each interval.
fn overlie_chrom(tracks: &[&SignalTrack<f32>], op: Op) -> Vec<(Coord, f32)> {
    let runs: Vec<&[macs_rle::Run<f32>]> = tracks.iter().map(|t| t.runs()).collect();
    let mut idx = vec![0usize; tracks.len()];
    let mut pre_p: Coord = 0;
    let mut out: Vec<(Coord, f32)> = Vec::new();
    loop {
        // current end + value of each track; stop if any is exhausted
        let mut lowest = Coord::MAX;
        let mut vals = Vec::with_capacity(tracks.len());
        for (i, rr) in runs.iter().enumerate() {
            let Some(r) = rr.get(idx[i]) else {
                return out; // StopIteration: halt, matching upstream
            };
            lowest = lowest.min(r.end);
            vals.push(r.value);
        }
        let v = op.apply(&vals);
        // emit [pre_p, lowest) with value v, merging equal neighbours
        if lowest > pre_p {
            match out.last_mut() {
                Some(last) if last.1 == v => last.0 = lowest,
                _ => out.push((lowest, v)),
            }
        }
        pre_p = lowest;
        // advance every track sitting at `lowest`
        for (i, rr) in runs.iter().enumerate() {
            if let Some(r) = rr.get(idx[i]) {
                if r.end == lowest {
                    idx[i] += 1;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bg(pairs: &[(&str, Coord, Coord, f32)], baseline: f32) -> BedGraph {
        let mut b = BedGraph::new(baseline);
        for (c, s, e, v) in pairs {
            b.add_loc(c.as_bytes(), *s, *e, *v);
        }
        b
    }

    #[test]
    fn add_loc_backfills_baseline_and_merges() {
        // first record starts at 100 -> baseline [0,100), then [100,200)
        let b = bg(&[("chr1", 100, 200, 5.0)], 0.0);
        let c = b.genome().get(b"chr1").unwrap();
        let t = b.track(c).unwrap();
        assert_eq!(t.runs().len(), 2);
        assert_eq!((t.runs()[0].end, t.runs()[0].value), (100, 0.0));
        assert_eq!((t.runs()[1].end, t.runs()[1].value), (200, 5.0));
        // equal-valued neighbour extends instead of duplicating
        let b = bg(&[("chr1", 0, 100, 5.0), ("chr1", 100, 200, 5.0)], 0.0);
        let c = b.genome().get(b"chr1").unwrap();
        assert_eq!(b.track(c).unwrap().runs().len(), 1);
    }

    #[test]
    fn add_loc_drops_nonpositive_end() {
        // `add_loc` returns before interning when end <= 0, so nothing is recorded
        let b = bg(&[("chr1", 0, 0, 1.0)], 0.0);
        assert!(b.chroms_sorted().is_empty());
    }

    #[test]
    fn overlie_max_matches_the_documented_example() {
        // upstream docstring example:
        //  #1: (0,100,0) (100,200,3) (200,300,4)
        //  #2: (0,150,1) (150,250,2) (250,300,4)
        //  max -> 1,3,3,4,4 over the five intervals
        let a = bg(
            &[
                ("c", 0, 100, 0.0),
                ("c", 100, 200, 3.0),
                ("c", 200, 300, 4.0),
            ],
            0.0,
        );
        let b2 = bg(
            &[
                ("c", 0, 150, 1.0),
                ("c", 150, 250, 2.0),
                ("c", 250, 300, 4.0),
            ],
            0.0,
        );
        let m = a.overlie(&[&b2], Op::Max);
        let c = m.genome().get(b"c").unwrap();
        let got: Vec<(Coord, f32)> = m
            .track(c)
            .unwrap()
            .runs()
            .iter()
            .map(|r| (r.end, r.value))
            .collect();
        // upstream's docstring lists five intervals 1,3,3,4,4; as an RLE bedGraph the
        // two adjacent 3s (and the two 4s) coalesce, so the runs are 1,3,4.
        assert_eq!(got, vec![(100, 1.0), (200, 3.0), (300, 4.0)]);
    }

    #[test]
    fn overlie_subtract_is_second_minus_first() {
        let a = bg(&[("c", 0, 100, 1.0)], 0.0);
        let b2 = bg(&[("c", 0, 100, 4.0)], 0.0);
        let m = a.overlie(&[&b2], Op::Subtract);
        let c = m.genome().get(b"c").unwrap();
        assert_eq!(m.track(c).unwrap().runs()[0].value, 3.0); // 4 - 1
    }

    #[test]
    fn overlie_emits_only_common_prefix_and_common_chroms() {
        // b is shorter: the walk stops when b is exhausted -> only [0,100)
        let a = bg(&[("c", 0, 200, 1.0)], 0.0);
        let b2 = bg(&[("c", 0, 100, 1.0)], 0.0);
        let m = a.overlie(&[&b2], Op::Max);
        let c = m.genome().get(b"c").unwrap();
        let got: Vec<(Coord, f32)> = m
            .track(c)
            .unwrap()
            .runs()
            .iter()
            .map(|r| (r.end, r.value))
            .collect();
        assert_eq!(got, vec![(100, 1.0)]);
        // disjoint chromosomes -> nothing
        let a = bg(&[("x", 0, 100, 1.0)], 0.0);
        let b2 = bg(&[("y", 0, 100, 1.0)], 0.0);
        let m = a.overlie(&[&b2], Op::Max);
        assert!(m.chroms_sorted().is_empty());
    }

    #[test]
    fn extract_value_reports_self_where_other_is_positive() {
        let a = bg(&[("c", 0, 100, 7.0), ("c", 100, 200, 9.0)], 0.0);
        // other positive only on [100,200)
        let b2 = bg(&[("c", 0, 100, 0.0), ("c", 100, 200, 1.0)], 0.0);
        let c = a.genome().get(b"c").unwrap();
        let other = b2.track(b2.genome().get(b"c").unwrap()).unwrap();
        let got = a.extract_value(c, other);
        assert_eq!(got, vec![(100, 200, 9.0)]);
    }

    /// `apply_func` preserves interval boundaries even when neighbours map to the
    /// same value.
    ///
    /// This test previously asserted the opposite -- that two equal-valued runs
    /// coalesce into one -- which contradicted `BedGraph.py:888-901` and its own
    /// docstring, "Two adjacent regions with same value after applying func will not
    /// be merged". It passed only because the implementation used the coalescing
    /// push; the pair encoded the bug rather than catching it.
    ///
    /// The input values differ (1.0 and 1.5) because `add_loc` merges equal-valued
    /// neighbours, as upstream's does, so two already-equal runs cannot be built
    /// through it. Both still map to 1.0, which is the case that matters.
    #[test]
    fn apply_func_maps_values_without_merging_neighbours() {
        let a = bg(&[("c", 0, 100, 1.0), ("c", 100, 200, 1.5)], 0.0);
        let c0 = a.genome().get(b"c").unwrap();
        assert_eq!(a.track(c0).unwrap().runs().len(), 2, "input has two runs");
        let b = a.apply_func(|v| if v < 2.0 { 1.0 } else { 2.0 });
        let c = b.genome().get(b"c").unwrap();
        let runs = b.track(c).unwrap().runs();
        assert_eq!(runs.len(), 2, "both intervals survive: {runs:?}");
        assert!(runs.iter().all(|r| r.value == 1.0));
        // `Run` carries only `end`; the start of run *n* is the end of run *n-1*.
        assert_eq!(b.track(c).unwrap().start(), 0);
        assert_eq!(runs[0].end, 100);
        assert_eq!(runs[1].end, 200);
    }

    /// The `-m max` shape: every value below the threshold clamps to one constant,
    /// and the row count must not change.
    #[test]
    fn apply_func_max_keeps_one_row_per_input_row() {
        let a = bg(
            &[
                ("c", 0, 100, 1.0),
                ("c", 100, 200, 2.0),
                ("c", 200, 300, 3.0),
                ("c", 300, 400, 4.0),
            ],
            0.0,
        );
        let b = a.apply_func(|v| if v > 100.0 { v } else { 100.0 });
        let c = b.genome().get(b"c").unwrap();
        let runs = b.track(c).unwrap().runs();
        assert_eq!(runs.len(), 4, "{runs:?}");
        assert!(runs.iter().all(|r| r.value == 100.0));
    }

    #[test]
    fn op_parses_upstream_names() {
        assert_eq!(Op::parse(b"max"), Some(Op::Max));
        assert_eq!(Op::parse(b"subtract"), Some(Op::Subtract));
        assert_eq!(Op::parse(b"nope"), None);
    }

    /// Every value below was produced by *upstream* MACS3 3.0.5
    /// (`bedGraphTrackI.overlie`) on this input pair, dumped with `%.6g`:
    ///
    /// ```text
    /// a.bdg: chr1 0 100 0 / 100 200 3 / 200 300 4
    /// b.bdg: chr1 0 150 1 / 150 250 2 / 250 300 4
    /// ```
    fn upstream_pair() -> (BedGraph, BedGraph) {
        let a = bg(
            &[
                ("chr1", 0, 100, 0.0),
                ("chr1", 100, 200, 3.0),
                ("chr1", 200, 300, 4.0),
            ],
            0.0,
        );
        let b = bg(
            &[
                ("chr1", 0, 150, 1.0),
                ("chr1", 150, 250, 2.0),
                ("chr1", 250, 300, 4.0),
            ],
            0.0,
        );
        (a, b)
    }

    fn overlie_runs(a: &BedGraph, b: &BedGraph, op: Op) -> Vec<(Coord, f32)> {
        let m = a.overlie(&[b], op);
        let c = m.genome().get(b"chr1").expect("chr1 in result");
        m.track(c)
            .unwrap()
            .runs()
            .iter()
            .map(|r| (r.end, r.value))
            .collect()
    }

    #[test]
    fn overlie_matches_upstream_for_every_function() {
        let (a, b) = upstream_pair();
        // max: adjacent equal values coalesce, exactly as upstream's add_loc does
        assert_eq!(
            overlie_runs(&a, &b, Op::Max),
            vec![(100, 1.0), (200, 3.0), (300, 4.0)]
        );
        assert_eq!(
            overlie_runs(&a, &b, Op::Sum),
            vec![(100, 1.0), (150, 4.0), (200, 5.0), (250, 6.0), (300, 8.0)]
        );
        // subtract is b - a (upstream's x[1] - x[0])
        assert_eq!(
            overlie_runs(&a, &b, Op::Subtract),
            vec![
                (100, 1.0),
                (150, -2.0),
                (200, -1.0),
                (250, -2.0),
                (300, 0.0)
            ]
        );
        assert_eq!(
            overlie_runs(&a, &b, Op::Product),
            vec![(100, 0.0), (150, 3.0), (200, 6.0), (250, 8.0), (300, 16.0)]
        );
        assert_eq!(
            overlie_runs(&a, &b, Op::Mean),
            vec![(100, 0.5), (150, 2.0), (200, 2.5), (250, 3.0), (300, 4.0)]
        );
    }

    #[test]
    fn overlie_fisher_matches_upstream() {
        let (a, b) = upstream_pair();
        let got = overlie_runs(&a, &b, Op::Fisher);
        let want = [0.481146f32, 2.99096, 3.90264, 4.82928, 6.71174];
        assert_eq!(got.len(), want.len());
        for (i, (&(_, g), &w)) in got.iter().zip(want.iter()).enumerate() {
            assert!(
                (g - w).abs() <= 1e-4 * w.abs().max(1.0),
                "fisher[{i}]: got {g} want {w}"
            );
        }
    }
}
