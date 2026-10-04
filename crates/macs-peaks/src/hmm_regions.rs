//! `MACS3.Signal.Region.Regions` -- the interval container `hmmratac` uses for
//! its training and candidate regions.
//!
//! A port of `Region.py`. Unlike the peak machinery elsewhere in this crate,
//! `Regions` stores real `(start, end)` intervals rather than end-based runs,
//! because `hmmratac` expands, merges and pops them before binning.
//!
//! # `merge_overlap` runs once, then it is skipped
//!
//! `self.__merged` latches: `expand` clears it, `merge_overlap` sets it, and a
//! second call returns immediately without doing anything. Since `expand` sorts
//! and `merge_overlap` sorts, calling either twice is *not* idempotent in the
//! naive sense -- the latch is what makes the sequence
//! `expand; merge_overlap; merge_overlap` equal to `expand; merge_overlap`.
//! [`RegionSet::merge_overlap`] reproduces the latch because `hmmratac` calls
//! it exactly once per region set, but `pop` does not clear it (upstream's
//! `pop` does not either), so a popped-from set stays merged.
//!
//! # `pop` takes a chromosome-major prefix, not a genome-wide one
//!
//! `pop(n)` walks chromosomes in sorted order and gives each one a slice of the
//! *remaining* budget: it slices `[:n_taken]` off the front and decrements
//! `n_taken` by however many it actually took (including 0). So `pop(10000)` on
//! a set holding 3 regions on `chr1` and 5 on `chr2` takes all 8 and leaves the
//! set empty -- it never runs past `chr1`.
//!
//! # Expansion is not clipped to the contig end
//!
//! `expand` clamps the start at 0 but adds `flanking` to the end unconditionally,
//! so a region near the contig end is pushed past it. The bins generated from it
//! are then clipped by the binning loop's `while r < e`, not by the region.

use std::collections::BTreeMap;

use macs_core::genome::ChromId;
use macs_core::{Coord, Genome};

/// A set of `(start, end)` intervals, keyed by chromosome.
#[derive(Debug, Clone, Default)]
pub struct RegionSet {
    genome: Genome,
    regions: BTreeMap<ChromId, Vec<(Coord, Coord)>>,
    total: usize,
    merged: bool,
}

impl RegionSet {
    /// An empty set.
    pub fn new(genome: Genome) -> Self {
        Self {
            genome,
            regions: BTreeMap::new(),
            total: 0,
            merged: false,
        }
    }

    /// The genome the regions are indexed against.
    pub fn genome(&self) -> &Genome {
        &self.genome
    }

    /// Number of regions.
    pub fn total(&self) -> usize {
        self.total
    }

    /// True when there is nothing left to iterate.
    pub fn is_empty(&self) -> bool {
        self.total == 0
    }

    /// `add_loc`.
    pub fn add(&mut self, chrom: ChromId, start: Coord, end: Coord) {
        self.regions.entry(chrom).or_default().push((start, end));
        self.total += 1;
    }

    /// Add a whole chromosome's list at once, as `init_from_PeakIO` does.
    pub fn extend_chrom(&mut self, chrom: ChromId, v: Vec<(Coord, Coord)>) {
        self.total += v.len();
        self.regions.insert(chrom, v);
    }

    /// Regions on one chromosome, ascending.
    pub fn chrom(&self, chrom: ChromId) -> &[(Coord, Coord)] {
        self.regions.get(&chrom).map_or(&[], |v| v.as_slice())
    }

    /// Chromosome ids in name order.
    pub fn chroms_sorted(&self) -> Vec<ChromId> {
        let mut ids: Vec<ChromId> = self.regions.keys().copied().collect();
        ids.sort_by(|&a, &b| self.genome.name(a).cmp(self.genome.name(b)));
        ids
    }

    /// `sort`: per-chromosome ascending by `(start, end)`.
    pub fn sort(&mut self) {
        for v in self.regions.values_mut() {
            v.sort();
        }
    }

    /// `expand`: grow every region by `flanking` on both sides, clamping the
    /// start at 0 but **not** the end.
    pub fn expand(&mut self, flanking: Coord) {
        self.sort();
        for v in self.regions.values_mut() {
            for r in v.iter_mut() {
                r.0 = r.0.saturating_sub(flanking);
                r.1 += flanking;
            }
            v.sort();
        }
        self.merged = false;
    }

    /// `merge_overlap`: coalesce regions that overlap (`start <= prev.end`).
    ///
    /// Note `start <= prev.end`, i.e. **abutting** regions merge too.
    pub fn merge_overlap(&mut self) {
        if self.merged {
            return;
        }
        self.sort();
        let mut out: BTreeMap<ChromId, Vec<(Coord, Coord)>> = BTreeMap::new();
        self.total = 0;
        for (&chrom, rs) in &self.regions {
            let mut merged: Vec<(Coord, Coord)> = Vec::new();
            let mut prev: Option<(Coord, Coord)> = None;
            for r in rs {
                match prev {
                    None => prev = Some(*r),
                    Some(p) => {
                        if r.0 <= p.1 {
                            prev = Some((p.0, r.1));
                        } else {
                            merged.push(p);
                            prev = Some(*r);
                        }
                    }
                }
            }
            if let Some(p) = prev {
                merged.push(p);
            }
            self.total += merged.len();
            out.insert(chrom, merged);
        }
        self.regions = out;
        self.sort();
        self.merged = true;
    }

    /// `pop(n)`: remove and return a chromosome-major prefix of at most `n`
    /// regions.
    ///
    /// Each chromosome takes `min(remaining, len)`, so a set with fewer than `n`
    /// regions in total is drained completely even if the first chromosome is
    /// small. `total` is recomputed from what remains.
    pub fn pop(&mut self, n: usize) -> Option<RegionSet> {
        if self.total == 0 {
            return None;
        }
        let mut ret = RegionSet::new(self.genome.clone());
        let mut n_taken = n;
        for chrom in self.chroms_sorted() {
            let src = self.regions.get(&chrom).cloned().unwrap_or_default();
            let take = n_taken.min(src.len());
            ret.regions.insert(chrom, src[..take].to_vec());
            self.regions.insert(chrom, src[take..].to_vec());
            if self.regions[&chrom].is_empty() {
                self.regions.remove(&chrom);
            }
            n_taken -= take;
        }
        self.total = self.regions.values().map(Vec::len).sum();
        ret.total = ret.regions.values().map(Vec::len).sum();
        Some(ret)
    }

    /// `write_to_bed`: `chrom\tstart\tend`, chromosomes in name order.
    pub fn to_bed_string(&self) -> String {
        let mut s = String::new();
        for chrom in self.chroms_sorted() {
            for &(rs, re) in self.chrom(chrom) {
                let _ = re;
                s.push_str(&format!(
                    "{}\t{}\t{}\n",
                    String::from_utf8_lossy(self.genome.name(chrom)),
                    rs,
                    re
                ));
            }
        }
        s
    }

    /// Flatten to `(chrom, start, end)` in name order -- the form the hmmratac
    /// binner wants.
    pub fn flatten(&self) -> Vec<(ChromId, Coord, Coord)> {
        let mut out = Vec::with_capacity(self.total);
        for chrom in self.chroms_sorted() {
            for &(s, e) in self.chrom(chrom) {
                out.push((chrom, s, e));
            }
        }
        out
    }
}

/// `PeakIO.filter_score`: keep peaks whose `fold_change` lies in
/// `[low_cutoff, high_cutoff)`.
///
/// The comparison is `low <= score < high` -- the lower bound is inclusive and
/// the upper is exclusive. `hmmratac` uses this to keep training regions whose
/// fold-change sits inside the user's `[--lower, --upper)` window.
pub fn filter_fc<I: Iterator<Item = (ChromId, Coord, Coord, f32)>>(
    peaks: I,
    low: f32,
    high: f32,
) -> Vec<(ChromId, Coord, Coord, f32)> {
    peaks.filter(|p| p.3 >= low && p.3 < high).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(name: &[u8]) -> Genome {
        let mut g = Genome::new();
        let _ = g.intern(name);
        g
    }

    #[test]
    fn expand_clamps_start_but_not_end() {
        let gm = g(b"chr1");
        let c = gm.get(b"chr1").unwrap();
        let mut r = RegionSet::new(gm);
        r.add(c, 50, 150);
        r.add(c, 200_000, 200_010);
        r.expand(100);
        assert_eq!(r.chrom(c), &[(0, 250), (199_900, 200_110)]);
    }

    #[test]
    fn merge_overlap_joins_abutting_regions() {
        let gm = g(b"chr1");
        let c = gm.get(b"chr1").unwrap();
        let mut r = RegionSet::new(gm);
        r.add(c, 100, 200);
        r.add(c, 200, 300); // abuts: start == prev end
        r.add(c, 400, 500);
        r.merge_overlap();
        assert_eq!(r.total(), 2);
        assert_eq!(r.chrom(c), &[(100, 300), (400, 500)]);
    }

    #[test]
    fn merge_overlap_latches_unless_expanded_again() {
        let gm = g(b"chr1");
        let c = gm.get(b"chr1").unwrap();
        let mut r = RegionSet::new(gm);
        r.add(c, 100, 200);
        r.add(c, 150, 300);
        r.merge_overlap();
        r.add(c, 250, 400);
        r.merge_overlap(); // latched -> 400 stays separate
        assert_eq!(r.chrom(c), &[(100, 300), (250, 400)]);
        r.expand(0); // clears the latch
        r.merge_overlap();
        assert_eq!(r.chrom(c), &[(100, 400)]);
    }

    #[test]
    fn pop_is_chromosome_major_and_drains_when_short() {
        let mut gm = g(b"chr1");
        let c1 = gm.intern(b"chr1");
        let c2 = gm.intern(b"chr2");
        let mut r = RegionSet::new(gm);
        for i in 0..3 {
            r.add(c1, i * 10, i * 10 + 5);
        }
        for i in 0..5 {
            r.add(c2, i * 10, i * 10 + 5);
        }
        let taken = r.pop(10_000).unwrap();
        assert_eq!(taken.total(), 8, "pop(10000) drains all 8");
        assert_eq!(r.total(), 0);
        assert!(r.pop(1).is_none());
        // chromosome-major: chr1 first, then chr2
        assert_eq!(taken.chroms_sorted(), vec![c1, c2]);
    }

    #[test]
    fn pop_budget_is_shared_across_chromosomes() {
        let mut gm = g(b"chr1");
        let c1 = gm.intern(b"chr1");
        let c2 = gm.intern(b"chr2");
        let mut r = RegionSet::new(gm);
        for i in 0..3 {
            r.add(c1, i * 10, i * 10 + 5);
        }
        for i in 0..5 {
            r.add(c2, i * 10, i * 10 + 5);
        }
        let taken = r.pop(5).unwrap();
        assert_eq!(taken.total(), 5);
        assert_eq!(taken.chrom(c1).len(), 3);
        assert_eq!(taken.chrom(c2).len(), 2);
        assert_eq!(r.total(), 3);
    }

    #[test]
    fn filter_score_bounds_are_half_open() {
        let gm = g(b"chr1");
        let c = gm.get(b"chr1").unwrap();
        let all = vec![
            (c, 0, 10, 0.9),
            (c, 10, 20, 1.0),
            (c, 20, 30, 1.5),
            (c, 30, 40, 2.0),
        ];
        let kept = filter_fc(all.into_iter(), 1.0, 2.0);
        assert_eq!(kept.len(), 2, "1.0 included, 2.0 excluded");
    }
}
