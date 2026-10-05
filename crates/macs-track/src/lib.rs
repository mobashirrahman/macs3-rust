//! Read and fragment tracks: storage, duplicate filtering, and down-sampling.
//!
//! Transcribed from `MACS3.Signal.FixWidthTrack` (`FWTrack`) and
//! `MACS3.Signal.PairedEndTrack` (`PETrackI` / `PETrackII`). The load-bearing
//! details are all things that look like oversights upstream but are load-bearing
//! for output equality:
//!
//! # Duplicate filtering keeps the *first* `maxnum` copies, not a sample
//!
//! `FWTrack.filter_dup` walks the sorted positions and keeps a position while the
//! running run-length `n <= maxnum`. Because the array is sorted, that is a prefix
//! of each run — so `--keep-dup 1` keeps exactly one copy of every distinct
//! position and the survivor's value is irrelevant (they are all equal). This is
//! a *filter*, not a random subsample, so it is deterministic and needs no seed.
//!
//! # A chromosome with a single position is passed through untouched
//!
//! Upstream tests `if len(plus) <= 1: new_plus = plus` and skips the filtering
//! body *entirely* — including skipping the `self.total` and `self.pointer`
//! updates for that chromosome. So a chromosome holding exactly one 5' end
//! reports its position unconditionally, even when `maxnum` is 0, and it is not
//! counted in `total`. Reproduced, including the counter, because a divergence in
//! `total` propagates into lambda and from there into every peak coordinate.
//!
//! # `PETrack.filter_dup` compares both endpoints
//!
//! Paired-end duplicates are identified by the `(start, end)` *pair*, not the start
//! alone, and the filtered track's `average_template_length` is recomputed from
//! the survivors.
//!
//! # Sampling uses NumPy's global RNG
//!
//! `sample_percent` / `sample_num` call `np.random.shuffle` after seeding the
//! *global* NumPy state, then re-sort. Reproducing that exactly requires NumPy's
//! Mersenne-Twister stream, which is implemented in `macs-stats`; see
//! [`sample_percent`].

#![forbid(unsafe_code)]
#![deny(missing_debug_implementations)]

mod frag;
mod sample;
mod single;

pub use frag::{filter_frag_dup, FragTrackBuilder, Fragment, FragmentTrack, PETrackI, PETrackII};
pub use sample::{retained_count, sample_frag_percent, sample_num, sample_percent};

use macs_core::{ChromId, Genome, Result};

/// Sorted per-strand positions, interned onto a genome.
#[derive(Debug, Default, Clone)]
pub struct SortedPositions {
    genome: Genome,
    plus: Vec<Vec<macs_core::Coord>>,
    minus: Vec<Vec<macs_core::Coord>>,
    total: u64,
    /// Whether [`finalize`](Self::finalize) has run. Mirrors upstream's flag so the
    /// samplers can skip a redundant sort exactly where upstream does.
    sorted: bool,
}

impl SortedPositions {
    /// An empty track.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append one 5' end, interning the chromosome if new.
    pub fn push(&mut self, chrom: &[u8], pos: macs_core::Coord, strand: macs_core::Strand) {
        let id = self.intern(chrom);
        let bucket = match strand {
            macs_core::Strand::Minus => &mut self.minus[id.0 as usize],
            _ => &mut self.plus[id.0 as usize],
        };
        bucket.push(pos);
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

    /// Number of positions appended.
    pub fn total(&self) -> u64 {
        self.total
    }

    /// The genome dictionary.
    pub fn genome(&self) -> &Genome {
        &self.genome
    }

    /// One chromosome's positions on one strand, ascending.
    pub fn strand(&self, chrom: ChromId, strand: macs_core::Strand) -> &[macs_core::Coord] {
        match strand {
            macs_core::Strand::Minus => &self.minus[chrom.0 as usize],
            _ => &self.plus[chrom.0 as usize],
        }
    }

    /// Every chromosome id, in **sorted byte order** of the chromosome name.
    ///
    /// Not file order. Upstream's `get_chr_names` returns `set(sorted(keys))` — the
    /// `sorted()` is discarded by the `set()` — and every stage iterates that set,
    /// so upstream's chromosome processing order is CPython set order, which is
    /// randomised per process for `bytes` keys (F29).
    ///
    /// Sorted byte order is the only stable choice, and the one upstream's own
    /// `sorted()` was reaching for. It matters only for the downsamplers, where
    /// chromosome visit order decides which array draws from the shared RNG stream
    /// first; counts are unaffected either way.
    pub fn chroms_sorted(&self) -> Vec<ChromId> {
        let mut ids = self.genome.ids_file_order();
        ids.sort_by(|a, b| self.genome.name(*a).cmp(self.genome.name(*b)));
        ids
    }

    /// Every chromosome id, in first-appearance order.
    ///
    /// This is the order upstream's *track storage* uses, and is what the parser
    /// and the pileup must use so that a chromosome's id is stable and diagnostics
    /// read in input order. Prefer [`chroms_sorted`](Self::chroms_sorted) for the
    /// samplers.
    pub fn chroms(&self) -> Vec<ChromId> {
        self.genome.ids_file_order()
    }

    /// Mutable access to one strand, for the samplers.
    pub(crate) fn strand_mut(
        &mut self,
        chrom: ChromId,
        strand: macs_core::Strand,
    ) -> &mut Vec<macs_core::Coord> {
        match strand {
            macs_core::Strand::Minus => &mut self.minus[chrom.0 as usize],
            _ => &mut self.plus[chrom.0 as usize],
        }
    }

    /// Sort every strand ascending.
    ///
    /// Upstream's `finalize` also shrinks the backing arrays to the fill count and
    /// resets `total`; this crate never zero-pads (see F19), so only the sort is
    /// observable.
    /// F210: release the spare capacity a growing `Vec` keeps.
    ///
    /// `Vec` doubles its allocation, so a track that grew to `n` elements can hold
    /// nearly `2n` of them. `sort_unstable` does not give that back, and neither does
    /// duplicate filtering, which only shortens the vector. Since the reads are held in
    /// memory for the whole of the pileup stage -- the peak-RSS moment -- that slack is
    /// charged directly against the peak.
    ///
    /// Measured on the 24-chromosome x 200 kbp benchmark (3 M + 1.5 M reads), this takes
    /// peak RSS from 232 MB to 224 MB. Real, but smaller than hoped: the slack is not
    /// where the memory actually goes. The tracks are -- see F211.
    pub fn finalize(&mut self) {
        for v in self.plus.iter_mut().chain(self.minus.iter_mut()) {
            v.sort_unstable();
            v.shrink_to_fit();
        }
        self.sorted = true;
    }

    /// `True` once [`finalize`](Self::finalize) has run.
    pub fn is_sorted(&self) -> bool {
        self.sorted
    }
}

/// Builder that accumulates positions and sorts them at the end.
///
/// Kept separate so the "insert order" state (which upstream models with a
/// `is_sorted` flag and a pointer array) is explicit rather than implicit.
#[derive(Debug, Default)]
pub struct SingleEndTrackBuilder {
    positions: SortedPositions,
    sorted: bool,
}

impl SingleEndTrackBuilder {
    /// An empty builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record one 5' end. Does not sort.
    pub fn push(&mut self, chrom: &[u8], pos: macs_core::Coord, strand: macs_core::Strand) {
        self.positions.push(chrom, pos, strand);
        self.sorted = false;
    }

    /// Sort, if not already sorted.
    pub fn finalize(&mut self) {
        if !self.sorted {
            self.positions.finalize();
            self.sorted = true;
        }
    }

    /// `true` once [`finalize`](Self::finalize) has run.
    pub fn is_sorted(&self) -> bool {
        self.sorted
    }

    /// Consume into a track.
    pub fn build(self) -> SingleEndTrack {
        SingleEndTrack::from_positions(self.positions)
    }

    /// The accumulated positions, for callers that want the raw view.
    pub fn positions(&self) -> &SortedPositions {
        &self.positions
    }
}

/// A single-end read track: sorted per-strand positions.
///
/// Upstream's `FWTrack` stores two arrays per chromosome (plus and minus 5'
/// positions) pre-sized to `buffer_size`. This crate never zero-pads (see F19), so
/// the arrays are exactly their fill length and `finalize` only has to sort.
#[derive(Debug, Clone, Default)]
pub struct SingleEndTrack {
    pub(crate) positions: SortedPositions,
    pub(crate) sorted: bool,
}

impl SingleEndTrack {
    /// An empty track.
    pub fn new() -> Self {
        Self::default()
    }

    /// Build a track directly from accumulated positions.
    pub fn from_positions(positions: SortedPositions) -> Self {
        let sorted = positions.is_sorted();
        SingleEndTrack { positions, sorted }
    }

    /// Number of retained 5' ends.
    pub fn total(&self) -> u64 {
        self.positions.total
    }

    /// The genome.
    pub fn genome(&self) -> &Genome {
        self.positions.genome()
    }

    /// The position storage.
    pub fn positions(&self) -> &SortedPositions {
        &self.positions
    }

    /// Mutable position storage, for in-place transforms.
    pub(crate) fn positions_mut(&mut self) -> &mut SortedPositions {
        &mut self.positions
    }

    /// `true` when the backing arrays are known to be ascending.
    pub fn is_sorted(&self) -> bool {
        self.sorted
    }

    /// Limit duplicates per position and strand, returning upstream's total.
    ///
    /// See the module docs for the two upstream quirks this reproduces: a strand
    /// holding at most one position bypasses the filter *and* the total update,
    /// and a run is truncated to its first `maxnum` entries rather than sampled.
    pub fn filter_dup(&mut self, maxnum: i64) -> Result<u64> {
        if maxnum < 0 {
            return Ok(self.positions.total);
        }
        if !self.sorted {
            self.positions.finalize();
            self.sorted = true;
        }
        Ok(self.positions.filter_dup(maxnum))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use macs_core::Strand;

    fn track(rows: &[(&[u8], u32, Strand)]) -> SingleEndTrack {
        let mut b = SingleEndTrackBuilder::new();
        for (c, p, st) in rows {
            b.push(c, *p, *st);
        }
        b.finalize();
        b.build()
    }

    #[test]
    fn keep_dup_one_keeps_exactly_one_copy_per_position() {
        // upstream truncates each sorted run to its first `maxnum` entries
        let mut t = track(&[
            (b"chr1", 10, Strand::Plus),
            (b"chr1", 10, Strand::Plus),
            (b"chr1", 10, Strand::Plus),
            (b"chr1", 20, Strand::Plus),
            (b"chr1", 20, Strand::Plus),
        ]);
        assert_eq!(t.filter_dup(1).unwrap(), 2);
    }

    #[test]
    fn strands_are_filtered_independently() {
        let mut t = track(&[
            (b"chr1", 10, Strand::Plus),
            (b"chr1", 10, Strand::Plus),
            (b"chr1", 10, Strand::Minus),
        ]);
        // upstream's `if len <= 1: continue` skips the `total` update, so the
        // surviving minus-strand position is retained but not counted. total is 1,
        // not 2, and both facts are load-bearing: `total` feeds the sampling
        // fraction and then lambda.
        assert_eq!(
            t.filter_dup(1).unwrap(),
            1,
            "the 1-element minus strand is uncounted"
        );
        let p = t.positions();
        let g = p.genome();
        assert_eq!(
            p.strand(g.get(b"chr1").unwrap(), Strand::Minus),
            &[10],
            "but it is still retained"
        );
    }

    #[test]
    fn keep_dup_two_keeps_two() {
        let mut t = track(&[
            (b"chr1", 10, Strand::Plus),
            (b"chr1", 10, Strand::Plus),
            (b"chr1", 10, Strand::Plus),
        ]);
        assert_eq!(t.filter_dup(2).unwrap(), 2);
    }

    #[test]
    fn negative_max_dup_is_a_no_op() {
        let mut t = track(&[(b"chr1", 10, Strand::Plus), (b"chr1", 10, Strand::Plus)]);
        assert_eq!(t.filter_dup(-1).unwrap(), 2);
    }

    #[test]
    fn a_single_position_bypasses_the_filter_and_the_total_update() {
        // upstream's `if len(plus) <= 1` branch skips the whole body, so the
        // position survives even at maxnum 0 and is not counted
        let mut t = track(&[(b"chr1", 42, Strand::Plus)]);
        assert_eq!(t.total(), 1);
        assert_eq!(t.filter_dup(0).unwrap(), 0, "not counted in total");
        let p = t.positions();
        let g = p.genome();
        assert_eq!(
            p.strand(g.get(b"chr1").unwrap(), Strand::Plus),
            &[42],
            "but the position itself is untouched"
        );
    }

    #[test]
    fn chromosomes_are_interned_in_first_appearance_order() {
        let t = track(&[
            (b"chrZ", 1, Strand::Plus),
            (b"chrA", 1, Strand::Plus),
            (b"chrZ", 2, Strand::Plus),
        ]);
        let g = t.genome();
        assert_eq!(g.name_string(ChromId(0)), "chrZ");
        assert_eq!(g.name_string(ChromId(1)), "chrA");
    }

    #[test]
    fn filtering_sorts_first_when_the_track_is_dirty() {
        let mut b = SingleEndTrackBuilder::new();
        // inserted out of order, and never finalized
        b.push(b"chr1", 30, Strand::Plus);
        b.push(b"chr1", 10, Strand::Plus);
        b.push(b"chr1", 30, Strand::Plus);
        let mut t = b.build();
        assert!(!t.is_sorted());
        assert_eq!(t.filter_dup(1).unwrap(), 2);
        assert!(t.is_sorted(), "filter_dup must sort, as upstream does");
        let p = t.positions();
        let g = p.genome();
        assert_eq!(p.strand(g.get(b"chr1").unwrap(), Strand::Plus), &[10, 30]);
    }
}
