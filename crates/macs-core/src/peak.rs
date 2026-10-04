//! Peak representation shared by every caller and writer.

use crate::{Coord, Interval, Score, Strand};

/// A called peak.
///
/// Coordinates are 0-based half-open, matching BED. `summit` is an absolute
/// genomic coordinate of the summit base; it is always inside
/// `[start, end)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Peak {
    /// Inclusive start.
    pub start: Coord,
    /// Exclusive end.
    pub end: Coord,
    /// Absolute summit coordinate.
    pub summit: Coord,
    /// Sum of pileup across the peak.
    pub pileup: Score,
    /// Peak's -log10 p-value.
    pub pscore: Score,
    /// Peak's -log10 q-value.
    pub qscore: Score,
    /// Fold enrichment against local lambda.
    pub fold_enrichment: Score,
    /// Width of the summit region, in bases.
    pub summit_width: Coord,
    /// Number of reads whose 5' end falls in the summit region.
    pub summit_reads: Count,
    /// Minimum q-score base pair inside the peak.
    pub min_qscore: Score,
    /// Number of base pairs in the peak.
    pub length: Coord,
    /// Number of tags in the peak.
    pub tags: Count,
    /// Number of tags in the summit region.
    pub summit_tags: Count,
    /// Local lambda at the summit, before smoothing.
    pub local_lambda: Score,
    /// Local lambda used for the call, after smoothing.
    pub lambda_input: Score,
    /// Strands represented inside the peak, when tracked.
    pub strand: Strand,
}

/// Alias kept short for read sites.
pub type Count = u64;

/// A peak together with its chromosome.
#[derive(Debug, Clone, PartialEq)]
pub struct PeakRecord {
    /// Interned chromosome.
    pub chrom: crate::ChromId,
    /// The peak itself.
    pub peak: Peak,
}

impl Peak {
    /// A peak with every score zeroed; used when a caller only knows geometry.
    pub fn geometry(start: Coord, end: Coord, summit: Coord) -> Self {
        Peak {
            start,
            end,
            summit,
            pileup: 0.0,
            pscore: 0.0,
            qscore: 0.0,
            fold_enrichment: 0.0,
            summit_width: 0,
            summit_reads: 0,
            min_qscore: 0.0,
            length: end.saturating_sub(start),
            tags: 0,
            summit_tags: 0,
            local_lambda: 0.0,
            lambda_input: 0.0,
            strand: Strand::Unknown,
        }
    }

    /// The peak as an [`Interval`].
    pub fn interval(&self) -> Interval {
        Interval::new(self.start, self.end)
    }

    /// Number of bases covered.
    pub fn len(&self) -> Coord {
        self.end.saturating_sub(self.start)
    }

    /// Always false: a peak has non-negative width by construction.
    pub fn is_empty(&self) -> bool {
        false
    }
}

/// Ordering used everywhere peaks are emitted: chromosome (per caller-chosen
/// order), then start, then end, then summit.
///
/// This is a *total* order with no ties beyond genuinely identical peaks, so
/// output is stable regardless of parallelism.
pub fn sort_peaks_by_position(peaks: &mut [(crate::ChromId, Peak)]) {
    peaks.sort_by(|(ca, pa), (cb, pb)| {
        ca.0.cmp(&cb.0)
            .then(pa.start.cmp(&pb.start))
            .then(pa.end.cmp(&pb.end))
            .then(pa.summit.cmp(&pb.summit))
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geometry_defaults() {
        let p = Peak::geometry(10, 20, 15);
        assert_eq!(p.len(), 10);
        assert!(!p.is_empty());
        assert_eq!(p.interval(), Interval::new(10, 20));
    }

    #[test]
    fn peak_sort_is_total() {
        let mut v = vec![
            (crate::ChromId(1), Peak::geometry(5, 10, 7)),
            (crate::ChromId(0), Peak::geometry(9, 12, 10)),
            (crate::ChromId(0), Peak::geometry(1, 4, 2)),
        ];
        sort_peaks_by_position(&mut v);
        assert_eq!(v[0].0, crate::ChromId(0));
        assert_eq!(v[0].1.start, 1);
        assert_eq!(v[1].1.start, 9);
        assert_eq!(v[2].0, crate::ChromId(1));
    }
}
