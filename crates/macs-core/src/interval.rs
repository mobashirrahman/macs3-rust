//! Half-open genomic intervals, `[start, end)`.
//!
//! Every interval in `macs3-rs` is 0-based half-open, matching the BED
//! convention that MACS3 uses for all BED-family output. XLS output uses a
//! different convention and is handled by the writer, not here.

use crate::{Coord, Len};
use std::cmp::{max, min};

/// A 0-based, half-open genomic interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Interval {
    start: Coord,
    end: Coord,
}

impl Interval {
    /// Construct an interval, clamping `end >= start`.
    ///
    /// Upstream occasionally produces an empty interval (e.g. a zero-length
    /// fragment); we normalise it to an empty interval rather than panicking.
    #[inline]
    pub fn new(start: Coord, end: Coord) -> Self {
        Interval {
            start,
            end: max(start, end),
        }
    }

    /// Construct a single-base interval at `pos` (i.e. `[pos, pos+1)`).
    #[inline]
    pub fn point(pos: Coord) -> Self {
        Interval {
            start: pos,
            end: pos + 1,
        }
    }

    /// Inclusive-start coordinate.
    #[inline]
    pub fn start(&self) -> Coord {
        self.start
    }

    /// Exclusive-end coordinate.
    #[inline]
    pub fn end(&self) -> Coord {
        self.end
    }

    /// Number of bases covered. Saturating, so a malformed input cannot wrap.
    #[inline]
    pub fn len(&self) -> Len {
        self.end.saturating_sub(self.start).into()
    }

    /// True when the interval covers no base.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }

    /// True when `pos` lies inside the half-open interval.
    #[inline]
    pub fn contains(&self, pos: Coord) -> bool {
        pos >= self.start && pos < self.end
    }

    /// True when the two intervals share at least one base.
    ///
    /// An **empty** interval shares no bases with anything, so it never overlaps --
    /// including another empty interval and including one that contains it. The
    /// emptiness guard is what makes this consistent with [`Self::overlap_len`], which
    /// is zero for exactly those cases; without it the two methods disagree, and a
    /// caller that used `overlaps` to mean "merge me" (the peak walk uses
    /// `gap == 0`, which is a *different*, deliberately looser test) would treat a
    /// zero-length fragment as a real overlap.
    #[inline]
    pub fn overlaps(&self, other: &Interval) -> bool {
        !self.is_empty() && !other.is_empty() && self.start < other.end && other.start < self.end
    }

    /// Length of the intersection, zero if disjoint.
    #[inline]
    pub fn overlap_len(&self, other: &Interval) -> Len {
        min(self.end, other.end)
            .saturating_sub(max(self.start, other.start))
            .into()
    }

    /// Smallest interval containing both.
    #[inline]
    pub fn union(&self, other: &Interval) -> Interval {
        Interval {
            start: min(self.start, other.start),
            end: max(self.end, other.end),
        }
    }

    /// Gap between two intervals, zero when they touch or overlap.
    ///
    /// Touching intervals (`self.end == other.start`) have gap 0; disjoint
    /// intervals have gap `other.start - self.end`. This definition is the one
    /// used by the peak mergers, so the `max_gap` boundary behaviour of
    /// `bdgpeakcall` is fully determined by it.
    ///
    /// Touching intervals give 0, so `max_gap = 0` merges only abutting or
    /// overlapping intervals. `bdgpeakcall`'s `max_gap` boundary behaviour is
    /// fully determined by this function.
    #[inline]
    pub fn gap(&self, other: &Interval) -> Coord {
        // `saturating_sub` rather than `-`: the branches already guarantee no
        // underflow, but `gap` feeds the peak mergers, where a wrong value is a
        // wrong peak, and saturating is a no-op on well-formed intervals.
        if other.start >= self.end {
            other.start.saturating_sub(self.end)
        } else if self.start >= other.end {
            self.start.saturating_sub(other.end)
        } else {
            0
        }
    }

    /// `self.end + pad`, saturating.
    #[inline]
    pub fn pad_right(&self, pad: Coord) -> Interval {
        Interval {
            start: self.start,
            end: self.end.saturating_add(pad),
        }
    }
}

/// Sort by start, then by end; ties keep input order because the sort is stable.
///
/// Stable ordering is mandatory: anywhere output order is observable we must
/// reproduce upstream's order, and upstream's order is the input order for ties.
pub fn sort_intervals(v: &mut [Interval]) {
    v.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));
}

/// Merge overlapping/adjacent intervals.
///
/// `max_gap` is the number of *uncovered* bases that may separate two intervals
/// and still be merged. `max_gap == 0` merges only touching or overlapping
/// intervals, because a gap of 0 is the degenerate touching case.
pub fn merge_intervals(v: &[Interval], max_gap: Coord) -> Vec<Interval> {
    if v.is_empty() {
        return Vec::new();
    }
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].start.cmp(&v[b].start).then(v[a].end.cmp(&v[b].end)));

    let mut out: Vec<Interval> = Vec::with_capacity(idx.len());
    let mut cur = v[idx[0]];
    for &i in &idx[1..] {
        let next = v[i];
        if next.gap(&cur) <= max_gap {
            cur = cur.union(&next);
        } else {
            out.push(cur);
            cur = next;
        }
    }
    out.push(cur);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn len_and_containment() {
        let iv = Interval::new(10, 20);
        assert_eq!(iv.len(), 10);
        assert!(iv.contains(10));
        assert!(iv.contains(19));
        assert!(!iv.contains(20));
        assert!(!iv.is_empty());
    }

    #[test]
    fn end_is_clamped_to_start() {
        let iv = Interval::new(20, 10);
        assert!(iv.is_empty());
        assert_eq!(iv.len(), 0);
    }

    #[test]
    fn gap_semantics_drive_max_gap_boundary() {
        let a = Interval::new(0, 100);
        // gap 49, 50, 51
        let b49 = Interval::new(149, 200);
        let b50 = Interval::new(150, 200);
        let b51 = Interval::new(151, 200);
        assert_eq!(a.gap(&b49), 49);
        assert_eq!(a.gap(&b50), 50);
        assert_eq!(a.gap(&b51), 51);

        // touching => gap 0, merged even with max_gap 0
        let touching = Interval::new(100, 200);
        assert_eq!(a.gap(&touching), 0);
        assert_eq!(merge_intervals(&[a, touching], 0).len(), 1);

        // The three-way boundary test required by gate G6.
        for (g, expect_merged) in [(49u32, true), (50, true), (51, false)] {
            let ivs = [a, Interval::new(100 + g, 200)];
            let merged = merge_intervals(&ivs, 50);
            assert_eq!(merged.len() == 1, expect_merged, "gap={g}");
        }
    }

    #[test]
    fn merge_is_order_independent_and_covers_input() {
        let ivs = vec![
            Interval::new(300, 400),
            Interval::new(0, 100),
            Interval::new(90, 200),
        ];
        let m = merge_intervals(&ivs, 10);
        assert_eq!(m, vec![Interval::new(0, 200), Interval::new(300, 400)]);
    }

    #[test]
    fn pad_right_saturates() {
        let iv = Interval::new(u32::MAX - 2, u32::MAX);
        assert_eq!(iv.pad_right(10).end(), u32::MAX);
    }
}
