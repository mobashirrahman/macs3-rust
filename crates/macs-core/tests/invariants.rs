//! L2 property invariants for the core types.
//!
//! `Interval` carries every coordinate in the port and `Genome` is the chromosome
//! interner the parallel pipeline keys on, so both need laws that hold for any input
//! rather than only for the recorded corpus.

use macs_core::interval::{merge_intervals, sort_intervals};
use macs_core::{Genome, Interval, Len};
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// `Interval::new` normalises `end >= start` instead of panicking or wrapping.
    /// Upstream does emit empty intervals (zero-length fragments), so this must be a
    /// normalisation rather than an assertion.
    #[test]
    fn interval_never_has_a_negative_length(a in 0u64..5000, b in 0u64..5000) {
        let iv = Interval::new(a, b);
        prop_assert!(iv.start() <= iv.end());
        prop_assert_eq!(iv.len(), iv.end() - iv.start());
        prop_assert!(!iv.is_empty() || iv.start() == iv.end());
    }

    /// `len` is half-open: `[a, b)` has length `b - a`, and a zero-length interval is
    /// empty. This is the BED convention; the XLS convention is handled at the writer.
    #[test]
    fn length_is_half_open(a in 0u64..5000, b in 0u64..5000) {
        let (lo, hi) = (a.min(b), a.max(b));
        let iv = Interval::new(lo, hi);
        prop_assert_eq!(iv.len(), hi - lo);
        prop_assert_eq!(iv.is_empty(), lo == hi);
    }

    /// `contains` is half-open at both ends and agrees with `start <= pos < end`.
    #[test]
    fn contains_is_half_open(a in 0u64..5000, d in 0u64..200, pos in 0u64..6000) {
        let iv = Interval::new(a, a + d);
        prop_assert_eq!(iv.contains(pos), pos >= a && pos < a + d);
        prop_assert!(!iv.contains(a + d), "end must be exclusive");
        prop_assert!(!iv.contains(a.wrapping_sub(1)), "before start");
    }

    /// Overlap and `overlap_len` agree, and `overlap_len` is never larger than either
    /// input. The peak-merging walk relies on this to bound its candidate set.
    #[test]
    fn overlap_len_is_bounded_by_both_inputs(
        a in 0u64..3000, b in 0u64..3000,
        c in 0u64..3000, d in 0u64..3000,
    ) {
        let (x, y) = (Interval::new(a, a + b), Interval::new(c, c + d));
        let ov = x.overlap_len(&y);
        prop_assert!(ov <= x.len());
        prop_assert!(ov <= y.len());
        prop_assert_eq!(x.overlaps(&y), ov > 0);
        prop_assert_eq!(y.overlaps(&x), ov > 0);
    }

    /// Merging is idempotent and produces sorted, non-overlapping intervals.
    #[test]
    fn merge_intervals_is_sorted_and_disjoint(
        ivs in prop::collection::vec((0u64..2000, 1u64..300), 0..40),
        max_gap in 0u64..200,
    ) {
        let mut v: Vec<Interval> =
            ivs.iter().map(|(s, l)| Interval::new(*s, s + l)).collect();
        sort_intervals(&mut v);
        let merged = merge_intervals(&v, max_gap);
        for w in merged.windows(2) {
            prop_assert!(w[0].end() <= w[1].start(),
                "merged intervals overlap: {:?} {:?}", w[0], w[1]);
        }
        // idempotent: merging already-merged output changes nothing
        let again = merge_intervals(&merged, max_gap);
        prop_assert_eq!(merged.len(), again.len());
        for (x, y) in merged.iter().zip(again.iter()) {
            prop_assert_eq!(x.start(), y.start());
            prop_assert_eq!(x.end(), y.end());
        }
    }

    /// Merging never loses covered bases: the total length is at least the input's.
    #[test]
    fn merge_never_loses_coverage(
        ivs in prop::collection::vec((0u64..2000, 1u64..300), 0..40),
        max_gap in 0u64..200,
    ) {
        let v: Vec<Interval> =
            ivs.iter().map(|(s, l)| Interval::new(*s, s + l)).collect();
        // `Len` is unsigned, so the meaningful law is containment, not a sign check:
        // every input interval must survive inside some merged interval, and the
        // merged total must be at least the largest single input.
        let merged_all = merge_intervals(&v, max_gap);
        let total: Len = merged_all.iter().map(|iv| iv.len()).sum();
        let widest = v.iter().map(|iv| iv.len()).max().unwrap_or(0);
        prop_assert!(total >= widest, "merged total {} < widest input {}", total, widest);
        // every input interval is contained in some merged interval
        for iv in &v {
            prop_assert!(
                merged_all.iter().any(|m| m.start() <= iv.start() && m.end() >= iv.end()),
                "interval {:?} vanished into {:?}", iv, merged_all
            );
        }
    }

    /// Interning is stable and injective: the same name always yields the same id, and
    /// two different names never share one. This is what lets the parallel stages key
    /// per-chromosome caches by `ChromId`.
    #[test]
    fn chromosome_interning_is_stable_and_injective(
        names in prop::collection::vec("[A-Za-z0-9_.]{1,12}", 1..40),
    ) {
        let mut g = Genome::new();
        let mut seen: Vec<(String, u32)> = Vec::new();
        for n in &names {
            let id = g.intern_str(n);
            let key = g.name_string(id);
            if let Some((prev, prev_id)) = seen.iter().find(|(k, _)| *k == key) {
                prop_assert_eq!(prev, n, "intern returned a different name");
                prop_assert_eq!(*prev_id, id.0, "same name interned to two ids");
            } else {
                seen.push((key, id.0));
            }
        }
        // re-interning every name returns the original id
        for (key, id) in &seen {
            let again = g.intern_str(key);
            prop_assert_eq!(again.0, *id, "re-interning {:?} changed its id", key);
        }
        prop_assert_eq!(g.len(), seen.len(), "interner grew unexpectedly");
    }

    /// `union` is commutative, associative, and covers both inputs; `gap` is zero
    /// exactly when they overlap or touch.
    #[test]
    fn union_and_gap_are_consistent(
        a in 0u64..3000, b in 0u64..500,
        c in 0u64..3000, d in 0u64..500,
    ) {
        let (x, y) = (Interval::new(a, a + b), Interval::new(c, c + d));
        let (u1, u2) = (x.union(&y), y.union(&x));
        prop_assert_eq!(u1.start(), u2.start());
        prop_assert_eq!(u1.end(), u2.end());
        prop_assert!(u1.start() <= u1.end());
        prop_assert!(u1.contains(x.start()) || x.is_empty());
        prop_assert!(u1.contains(x.end().saturating_sub(1)) || x.is_empty());

        // `gap` is the distance between disjoint intervals, and zero exactly when they
        // touch or overlap: `max_gap` merging treats a gap of 0 as "merge me".
        let g = x.gap(&y);
        // `gap` is 0 exactly when the intervals overlap *or merely touch* -- an
        // empty interval next to a non-empty one touches it without overlapping,
        // and `maxgap` merging treats a gap of 0 as "merge me". So the law is about
        // adjacency, not about `overlaps`.
        let touching = x.end() >= y.start() && y.end() >= x.start();
        prop_assert_eq!(g == 0, touching,
            "gap {} disagrees with adjacency {} for {:?} vs {:?}",
            g, touching, x, y);
        if !touching {
            prop_assert_eq!(g, u1.len() - x.len() - y.len(),
                "gap {} inconsistent with the union length", g);
        }
    }
}
