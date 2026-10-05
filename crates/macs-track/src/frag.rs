//! Paired-end storage. Filled in below.

use macs_core::{ChromId, Genome, Result};

/// A single fragment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Fragment {
    /// Left (5') end, 0-based inclusive.
    pub start: macs_core::Coord,
    /// Right (3') end, 0-based exclusive.
    pub end: macs_core::Coord,
}

/// A paired-end track: fragments per chromosome, plus mean template length.
///
/// `PETrackI` stores only `(start, end)`; `PETrackII` additionally carries the
/// FRAG barcode count in a `u16` field. Both share the deduplication and
/// sampling behaviour, so one type covers both with an optional count.
#[derive(Debug, Default)]
pub struct FragmentTrack {
    genome: Genome,
    frags: Vec<Vec<Fragment>>,
    counts: Vec<Vec<u16>>,
    /// `PETrackII` when true: every fragment carries a multiplicity.
    with_counts: bool,
    total: u64,
    /// Sum of `end - start` over all fragments, as upstream accumulates it.
    length: u64,
}

/// `Genome` is deliberately not `Clone` — interning order is identity — so this is
/// written out rather than derived.
impl Clone for FragmentTrack {
    fn clone(&self) -> Self {
        FragmentTrack {
            genome: self.genome.clone(),
            frags: self.frags.clone(),
            counts: self.counts.clone(),
            with_counts: self.with_counts,
            total: self.total,
            length: self.length,
        }
    }
}

impl FragmentTrack {
    /// An empty `PETrackI`-style track.
    pub fn new() -> Self {
        Self::default()
    }

    /// An empty `PETrackII`-style track, whose fragments carry counts.
    pub fn with_barcodes() -> Self {
        FragmentTrack {
            with_counts: true,
            ..Default::default()
        }
    }

    /// `True` when fragments carry a multiplicity.
    pub fn has_counts(&self) -> bool {
        self.with_counts
    }

    fn intern(&mut self, name: &[u8]) -> ChromId {
        if let Some(id) = self.genome.get(name) {
            return id;
        }
        let id = self.genome.intern(name);
        self.frags.push(Vec::new());
        self.counts.push(Vec::new());
        id
    }

    /// Append one fragment, interning the chromosome if new.
    pub fn push(&mut self, chrom: &[u8], frag: Fragment, count: u16) {
        let id = self.intern(chrom);
        self.frags[id.0 as usize].push(frag);
        self.counts[id.0 as usize].push(count);
        // `PETrackII.add_loc`: `self.length += (end - start) * count`. The
        // weighting only bites for a counted track, but it changes
        // `average_template_length`, which hmmratac uses as `min_length`.
        //
        // F208: the product must be computed *wrapping*. `length` is `cython.long` in
        // C, so upstream silently wraps on overflow; Rust's `*` panics under the
        // overflow checks that the fuzz and debug builds enable, and a FRAG record can
        // reach here with any coordinate and count the input supplies. That is a direct
        // violation of the "zero panics on malformed input" criterion, and it was found
        // by `fuzz/fuzz_targets/parse_frag.rs` within seconds of being run for the first
        // time -- the target had never been built, because the manifest pointed at a
        // file that did not exist.
        self.length = self.length.wrapping_add(
            u64::from(frag.end.saturating_sub(frag.start)).wrapping_mul(u64::from(count)),
        );
        self.total += 1;
    }

    /// Number of fragments.
    pub fn total(&self) -> u64 {
        self.total
    }

    /// Sum of template lengths, as upstream's `length`.
    pub fn length(&self) -> u64 {
        self.length
    }

    /// Upstream's `average_template_length`; `0.0` when empty.
    ///
    /// Upstream computes this in `finalize` as a Python float division, so it is
    /// `f64` and not an integer ratio.
    pub fn average_template_length(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        // F183/F184: upstream declares `average_template_length` as a
        // `cython.float`, so the mean is rounded to f32. That is not cosmetic here:
        // `PeakDetect` declares `control_sum` as a `cython.long` too, so
        //
        //     control_sum = int(control_total * average_template_length)
        //
        // truncates, and the two roundings decide which side of the integer the
        // product falls on. On `sweep/gmini_mpe_d1200_w600_ctrl_228` the f64 mean is
        // 119.05166666666667, giving 2400 * mean = 285724.00000000006 ->
        // `control_sum` 285724 and a `ratio_treat2control` of exactly 0.5; the f32
        // mean 119.05166625976562 gives 285723.9990234375 -> 285723 and
        // 0.500001749947, which is what the oracle reports.
        (self.length as f32 / self.total as f32) as f64
    }

    /// The genome.
    pub fn genome(&self) -> &Genome {
        &self.genome
    }

    /// Fragments on one chromosome, ascending by `(start, end)`.
    pub fn frags(&self, chrom: ChromId) -> &[Fragment] {
        &self.frags[chrom.0 as usize]
    }

    /// Counts on one chromosome, parallel to [`frags`](Self::frags).
    pub fn counts(&self, chrom: ChromId) -> &[u16] {
        &self.counts[chrom.0 as usize]
    }

    /// Every chromosome id, in first-appearance order.
    pub fn chroms(&self) -> Vec<ChromId> {
        self.genome.ids_file_order()
    }

    /// `fraglengths` (`PairedEndTrack.py:325`), concatenated in **first-appearance
    /// chromosome order** -- not sorted order.
    ///
    /// The order matters: `HMMR_EM` indexes this array positionally while
    /// downsampling, so any reordering changes which fragments land in the
    /// 100-1000bp window.
    pub fn fraglengths(&self) -> Vec<macs_core::Coord> {
        let mut out = Vec::new();
        for c in self.genome.ids_file_order() {
            let fs = &self.frags[c.0 as usize];
            let ks = &self.counts[c.0 as usize];
            if self.with_counts {
                // `PETrackII.fraglengths`: `np.repeat(sizes, counts)`
                for (f, k) in fs.iter().zip(ks.iter()) {
                    for _ in 0..*k {
                        out.push(f.end - f.start);
                    }
                }
            } else {
                for f in fs {
                    out.push(f.end - f.start);
                }
            }
        }
        out
    }

    /// `count_fraglengths` (`PairedEndTrack.py:300`): a histogram of fragment
    /// lengths, weighted by count when the track carries counts.
    ///
    /// `HMMR_EM.generate_weight_mapping` only iterates the keys, so chromosome
    /// order does not matter here, but the counts do: a `PETrackII` built from a
    /// barcode-collapsed FRAG file has every length weighted by its barcode
    /// count.
    pub fn count_fraglengths(&self) -> std::collections::BTreeMap<macs_core::Coord, u64> {
        let mut out = std::collections::BTreeMap::new();
        for c in self.genome.ids_file_order() {
            for (f, k) in self.frags[c.0 as usize]
                .iter()
                .zip(self.counts[c.0 as usize].iter())
            {
                *out.entry(f.end - f.start).or_insert(0) += u64::from(*k);
            }
        }
        out
    }

    /// `frags` for one chromosome as `(start, end)` pairs, ascending.
    pub fn pairs(&self, chrom: ChromId) -> Vec<(macs_core::Coord, macs_core::Coord)> {
        self.frags[chrom.0 as usize]
            .iter()
            .map(|f| (f.start, f.end))
            .collect()
    }

    /// Sort each chromosome's fragments by `(start, end)`.
    pub fn finalize(&mut self) {
        for (i, v) in self.frags.iter_mut().enumerate() {
            // counts are parallel, so they must be permuted with the fragments
            let mut pairs: Vec<(Fragment, u16)> = v
                .iter()
                .copied()
                .zip(self.counts[i].iter().copied())
                .collect();
            pairs.sort_unstable_by_key(|(f, _)| *f);
            v.clear();
            self.counts[i].clear();
            v.extend(pairs.iter().map(|(f, _)| *f));
            self.counts[i].extend(pairs.iter().map(|(_, c)| *c));
        }
        // `PETrackII.finalize`: `self.total += np.sum(locations[c]['c'])` --
        // the total is the **count** sum for a counted track, not the row count.
        // Since `length` is count-weighted too, `average_template_length` comes
        // out as the count-weighted mean fragment length either way; what changes
        // is the divisor, so using rows here would inflate it by the mean count.
        self.total = if self.with_counts {
            self.counts.iter().flatten().map(|c| u64::from(*c)).sum()
        } else {
            self.frags.iter().map(|v| v.len() as u64).sum::<u64>()
        };
    }
}

/// Builder for [`FragmentTrack`].
#[derive(Debug, Default)]
pub struct FragTrackBuilder {
    track: FragmentTrack,
}

impl FragTrackBuilder {
    /// An empty `PETrackI`-style builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// An empty builder for a counted (`PETrackII`-style) track.
    pub fn with_barcodes() -> Self {
        Self {
            track: FragmentTrack::with_barcodes(),
        }
    }

    /// Append a fragment.
    pub fn push(&mut self, chrom: &[u8], start: macs_core::Coord, end: macs_core::Coord) {
        self.track.push(chrom, Fragment { start, end }, 1);
    }

    /// Append a fragment with a multiplicity.
    pub fn push_with_count(
        &mut self,
        chrom: &[u8],
        start: macs_core::Coord,
        end: macs_core::Coord,
        count: u32,
    ) {
        let count = count.min(u32::from(u16::MAX)) as u16;
        self.track.push(chrom, Fragment { start, end }, count);
    }

    /// Sort by `(start, end)`.
    pub fn finalize(&mut self) {
        self.track.finalize();
    }

    /// Consume into a track.
    pub fn build(self) -> FragmentTrack {
        self.track
    }
}

/// `PETrackI`: fragments without barcodes.
pub type PETrackI = FragmentTrack;
/// `PETrackII`: fragments with barcode counts.
pub type PETrackII = FragmentTrack;

/// Duplicate filtering for paired-end tracks, matching `PETrack.filter_dup`.
///
/// Upstream identifies a duplicate by the full `(start, end)` pair, keeps the
/// first `maxnum` occurrences of each pair, subtracts the removed fragments' spans
/// from `length`, and recomputes `average_template_length`.
pub fn filter_frag_dup(track: &mut FragmentTrack, maxnum: i64) -> Result<u64> {
    if maxnum < 0 {
        return Ok(track.total);
    }
    let cap = maxnum.max(0) as usize;
    let n = track.frags.len();
    for idx in 0..n {
        // take the chromosome out so the loop can mutate in place without
        // holding a second borrow of `track`
        let mut frags = std::mem::take(&mut track.frags[idx]);
        if frags.len() <= 1 {
            // upstream: `if locs_size == 1: total += locs_size; continue`
            track.frags[idx] = frags;
            continue;
        }
        let mut counts = std::mem::take(&mut track.counts[idx]);
        // upstream initialises `selected_idx` to all-True, so the first fragment of
        // a run is kept unconditionally and only entries 1.. are candidates for
        // removal. Writing index 0 up front is what makes that true.
        let mut write = 0usize;
        write += 1; // frags[0] is already in place
        let mut run = 1usize;
        let mut removed_len: u64 = 0;
        for read in 1..frags.len() {
            let f = frags[read];
            if f == frags[read - 1] {
                run += 1;
                if run > cap.max(1) {
                    removed_len += u64::from(f.end.saturating_sub(f.start));
                    continue;
                }
            } else {
                run = 1;
            }
            frags[write] = f;
            counts[write] = counts[read];
            write += 1;
        }
        frags.truncate(write);
        counts.truncate(write);
        track.frags[idx] = frags;
        track.counts[idx] = counts;
        track.length = track.length.saturating_sub(removed_len);
    }
    let total: u64 = track.frags.iter().map(|v| v.len() as u64).sum();
    track.total = total;
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn f(start: u32, end: u32) -> Fragment {
        Fragment { start, end }
    }

    #[test]
    fn fragments_sort_by_start_then_end() {
        let mut b = FragTrackBuilder::new();
        b.push(b"chr1", 100, 200);
        b.push(b"chr1", 10, 20);
        b.push(b"chr1", 10, 30);
        b.finalize();
        let t = b.build();
        assert_eq!(
            t.frags(t.genome().get(b"chr1").unwrap()),
            &[f(10, 20), f(10, 30), f(100, 200)]
        );
    }

    #[test]
    fn average_template_length_matches_upstreams_float_division() {
        let mut b = FragTrackBuilder::new();
        b.push(b"chr1", 0, 100);
        b.push(b"chr1", 0, 101);
        let t = b.build();
        assert_eq!(t.length(), 201);
        assert_eq!(t.average_template_length(), 100.5);
    }

    #[test]
    fn empty_track_has_zero_average_length_not_nan() {
        let t = FragmentTrack::new();
        assert_eq!(t.average_template_length(), 0.0);
    }

    #[test]
    fn dup_filtering_compares_both_endpoints() {
        // same start, different end -> not duplicates
        let mut b = FragTrackBuilder::new();
        for _ in 0..3 {
            b.push(b"chr1", 10, 20);
        }
        b.push(b"chr1", 10, 21);
        b.finalize();
        let mut t = b.build();
        let before = t.total();
        assert_eq!(before, 4);
        let after = filter_frag_dup(&mut t, 1).unwrap();
        assert_eq!(after, 2, "one (10,20) survives, (10,21) survives");
    }

    #[test]
    fn dup_filtering_subtracts_removed_spans_from_length() {
        let mut b = FragTrackBuilder::new();
        for _ in 0..3 {
            b.push(b"chr1", 100, 200);
        }
        b.finalize();
        let mut t = b.build();
        assert_eq!(t.length(), 300);
        filter_frag_dup(&mut t, 1).unwrap();
        assert_eq!(t.length(), 100, "two 100bp fragments removed");
        assert_eq!(t.average_template_length(), 100.0);
    }

    #[test]
    fn negative_max_dup_is_a_no_op() {
        let mut b = FragTrackBuilder::new();
        b.push(b"chr1", 10, 20);
        b.push(b"chr1", 10, 20);
        b.finalize();
        let mut t = b.build();
        assert_eq!(filter_frag_dup(&mut t, -1).unwrap(), 2);
    }

    #[test]
    fn counts_stay_parallel_to_fragments_after_sorting() {
        let mut b = FragTrackBuilder::new();
        b.push_with_count(b"chr1", 100, 200, 5);
        b.push_with_count(b"chr1", 10, 20, 7);
        b.finalize();
        let t = b.build();
        let c = t.genome().get(b"chr1").unwrap();
        assert_eq!(t.frags(c), &[f(10, 20), f(100, 200)]);
        assert_eq!(
            t.counts(c),
            &[7, 5],
            "count 7 belongs to the (10,20) fragment"
        );
    }

    /// F208: a huge `count` on a long fragment must not overflow-panic.
    ///
    /// `PETrackII.add_loc` accumulates `length += (end - start) * count` into a C
    /// `long`, so upstream wraps; Rust's `*` panics under the overflow checks that debug
    /// and fuzz builds enable. Found by `fuzz/fuzz_targets/parse_frag.rs` the first time
    /// that target was ever run -- the manifest had pointed at a file that did not exist,
    /// so the whole fuzz layer had never been built.
    #[test]
    fn a_huge_count_on_a_long_fragment_does_not_overflow() {
        let mut b = FragTrackBuilder::new();
        // `push_with_count` clamps to u16::MAX, and u32::MAX * 65535 is past u64::MAX's
        // 2^64 by two orders of magnitude only when the fragment length is near 2^32;
        // with `Coord = u32` the widest possible fragment is u32::MAX, which is the
        // value used here -- the accumulation is still exact in u64, so this pins the
        // no-panic contract rather than the old wrapping one.
        let count = u32::from(u16::MAX);
        let long = u32::MAX;
        b.push_with_count(b"chr1", 0, long, count);
        b.finalize();
        let t = b.build();
        // Wrapping is the documented behaviour (it is what the C original does); what
        // matters is that asking for the average length does not panic.
        let _ = t.average_template_length();
        assert_eq!(t.total(), 1);
        assert_eq!(t.length(), u64::from(long).wrapping_mul(u64::from(count)));
    }
}
