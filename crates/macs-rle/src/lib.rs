//! Run-length-encoded genomic signal tracks.
//!
//! # Why this exists
//!
//! Upstream MACS3 materialises per-base NumPy arrays (`int32` positions +
//! `float32` values) and, for a 3 Gbp human chromosome, then pickles them to
//! disk per chromosome ([`MACS3/Signal/CallPeakUnit.py:554`]). The *values* it
//! actually needs, however, are the local lambda, the p-score and the q-score,
//! all of which are piecewise constant with a small number of breakpoints. This
//! crate is that representation.
//!
//! # Representation and invariants
//!
//! ```text
//! start = 0, end = 1000
//! runs  = [ Run { end: 120, value: 1.0 },
//!           Run { end: 150, value: 2.0 },
//!           Run { end: 180, value: 1.0 } ]
//! ```
//!
//! * `runs[i]` covers `[prev_end, runs[i].end)`, where `prev_end` is `start` for
//!   `i == 0`.
//! * **Contiguity.** Runs are contiguous: a track covers exactly
//!   `[start, last_end)` with no interior gaps. Only the tail
//!   `[last_end, end)` is uncovered. This is exactly upstream's structure: a
//!   `pv` array's value at `p[i]` applies to `[p[i], p[i+1])` and the final
//!   value applies to everything after it.
//! * **Canonical form.** Adjacent runs never carry equal values. Every
//!   constructor coalesces. This is a representation invariant only: coalescing
//!   does not change the pointwise value function, and everything downstream
//!   (histograms, thresholds, integral conservation, peak calling) is a function
//!   of the pointwise values. Upstream's `over_two_pv_array` does *not*
//!   coalesce; we do, and [`SignalTrack::zip`] documents why the difference is
//!   invisible.
//! * **Extents are explicit.** `end` is metadata, not part of the value
//!   function. `zip` restricts to the intersection of extents but the value of
//!   the *last* run persists to that extent, matching upstream's merge-walk.
//!
//! # Numeric type
//!
//! Upstream stores track values as `float32` (`dtype="f4"` throughout
//! `PileupV2.py`). Instantiating `SignalTrack<f32>` reproduces that exactly;
//! `SignalTrack<f64>` is for the sites where upstream widens (bedGraph writers
//! that format as `f64`, HMMRATAC) and each use must be justified per call site.

#![forbid(unsafe_code)]
#![deny(missing_debug_implementations)]

use std::collections::HashMap;
use std::fmt::Debug;

use macs_core::{ChromId, Coord, Interval, Len};

/// One constant-valued run, ending at `end` (exclusive).
#[derive(Debug, Clone, Copy, PartialEq, Hash)]
pub struct Run<T> {
    /// Exclusive end coordinate of this run.
    pub end: Coord,
    /// The value on the run.
    pub value: T,
}

impl<T> Run<T> {
    /// A run ending at `end` with `value`.
    pub fn new(end: Coord, value: T) -> Self {
        Run { end, value }
    }
}

/// A borrowed run with its resolved start.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Span<'a, T> {
    /// Inclusive start.
    pub start: Coord,
    /// Exclusive end.
    pub end: Coord,
    /// Value on `[start, end)`.
    pub value: &'a T,
}

impl<T> Span<'_, T> {
    /// Number of bases covered.
    pub fn len(&self) -> Len {
        (self.end - self.start).into()
    }

    /// Always false: spans are non-empty by construction.
    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }
}

/// A value usable as an exact, collision-free map key.
///
/// `f32` and `f64` deliberately do not implement `Eq + Hash` because NaN is not
/// reflexive, and a numeric-keyed `HashMap` would also merge `-0.0` with `0.0`.
/// MACS3's score histograms must not do either, so we key on the IEEE-754 bit
/// pattern.
///
/// The mapping is **order-preserving**: `a < b` implies `a.to_key() <
/// b.to_key()`, so sorting by key is sorting by value. For positive floats the
/// raw bits already increase with the value; for negative floats the sign bit
/// is inverted, which we achieve by flipping the top bit of the sign-extended
/// 64-bit key.
pub trait BitKey: Copy {
    /// The order-preserving bit-pattern key.
    fn to_key(self) -> u64;
    /// Recover a value from a key. Total for the float and integer types.
    fn from_key(k: u64) -> Self;
}

/// The classic order-preserving float-to-uint transform.
///
/// For a sign-magnitude IEEE-754 word, `mask` is all-ones for negatives (so
/// XOR inverts the whole word, and larger magnitudes give smaller keys) and just
/// the sign bit for non-negatives. The result is a `uN` whose unsigned ordering
/// matches the float ordering exactly, including across zero and for NaN.
#[inline]
const fn order_key_i32(bits: u32) -> u64 {
    let negative = bits >> 31; // 0 or 1
    let mask = 0u32.wrapping_sub(negative) | 0x8000_0000; // sign bit, or all ones
    (bits ^ mask) as u64
}

#[inline]
const fn order_key_i64(bits: u64) -> u64 {
    let negative = bits >> 63; // 0 or 1
    let mask = 0u64.wrapping_sub(negative) | 0x8000_0000_0000_0000; // sign bit, or all ones
    bits ^ mask
}

/// Inverse of [`order_key_i32`].
///
/// The mask is recovered from the key's own sign bit, which is the complement
/// of the original's: a non-negative float maps to a key with the top bit set,
/// a negative float to a key with the top bit clear.
#[inline]
const fn unorder_key_i32(k: u32) -> u32 {
    let mask = 0u32.wrapping_sub((!k) >> 31) | 0x8000_0000;
    k ^ mask
}

/// Inverse of [`order_key_i64`].
#[inline]
const fn unorder_key_i64(k: u64) -> u64 {
    let mask = 0u64.wrapping_sub((!k) >> 63) | 0x8000_0000_0000_0000;
    k ^ mask
}

impl BitKey for f32 {
    #[inline]
    fn to_key(self) -> u64 {
        order_key_i32(self.to_bits())
    }
    #[inline]
    fn from_key(k: u64) -> Self {
        f32::from_bits(unorder_key_i32(k as u32))
    }
}

impl BitKey for f64 {
    #[inline]
    fn to_key(self) -> u64 {
        order_key_i64(self.to_bits())
    }
    #[inline]
    fn from_key(k: u64) -> Self {
        f64::from_bits(unorder_key_i64(k))
    }
}

impl BitKey for i32 {
    #[inline]
    fn to_key(self) -> u64 {
        order_key_i32(self as u32)
    }
    #[inline]
    fn from_key(k: u64) -> Self {
        unorder_key_i32(k as u32) as i32
    }
}

impl BitKey for u64 {
    #[inline]
    fn to_key(self) -> u64 {
        self
    }
    #[inline]
    fn from_key(k: u64) -> Self {
        k
    }
}

/// A run-length-encoded signal over one chromosome.
///
/// See the module docs for the invariants this type maintains.
#[derive(Debug, Clone, PartialEq)]
pub struct SignalTrack<T> {
    chrom: ChromId,
    start: Coord,
    end: Coord,
    runs: Vec<Run<T>>,
}

impl<T> AsRef<SignalTrack<T>> for SignalTrack<T> {
    fn as_ref(&self) -> &SignalTrack<T> {
        self
    }
}

impl<T: Copy + PartialEq> SignalTrack<T> {
    /// An empty (fully uncovered) track over `[start, end)`.
    ///
    /// A reversed range is normalised to an empty range at the left edge.
    pub fn empty(chrom: ChromId, start: Coord, end: Coord) -> Self {
        SignalTrack {
            chrom,
            start,
            end: end.max(start),
            runs: Vec::new(),
        }
    }

    /// A track covering `[start, end)` with the constant value `v`.
    pub fn constant(chrom: ChromId, start: Coord, end: Coord, v: T) -> Self {
        let mut t = SignalTrack::empty(chrom, start, end);
        if t.end > t.start {
            t.runs.push(Run {
                end: t.end,
                value: v,
            });
        }
        t
    }

    /// Build from runs, coalescing equal adjacent values, clamping to the
    /// extent, and dropping zero-width and backwards runs.
    ///
    /// # Panics (debug only)
    /// If runs are supplied out of order, this is a caller bug. In release the
    /// offending run is dropped instead, so no malformed input can produce a
    /// corrupt value function.
    pub fn from_runs(chrom: ChromId, start: Coord, end: Coord, runs: Vec<Run<T>>) -> Self {
        Self::from_runs_impl(chrom, start, end, runs, true)
    }

    /// Build from runs, **preserving every breakpoint** the caller supplied.
    ///
    /// [`SignalTrack::from_runs`] coalesces runs whose values compare equal.
    /// That is correct for the bedGraph and local-lambda reporting paths, whose
    /// upstream counterparts also coalesce, but it is wrong for the peak
    /// pipeline's local-lambda fold: upstream's `over_two_pv_array` appends one
    /// entry per step of its lagging-pointer union walk into a raw numpy array
    /// and never merges neighbours, so the breakpoint set it hands downstream is
    /// part of the observable contract.
    ///
    /// The breakpoints are load-bearing well beyond the track itself. Peak
    /// selection indexes positions in the lambda track, so dropping equal-valued
    /// neighbours shrinks the index space and changes which regions are selected
    /// at all -- not merely how they are drawn. See `docs/upstream-findings.md`
    /// F75.
    ///
    /// Clamping to the extent and dropping zero-width and backwards runs still
    /// applies; only value-based coalescing is omitted.
    pub fn from_runs_exact(chrom: ChromId, start: Coord, end: Coord, runs: Vec<Run<T>>) -> Self {
        Self::from_runs_impl(chrom, start, end, runs, false)
    }

    fn from_runs_impl(
        chrom: ChromId,
        start: Coord,
        end: Coord,
        runs: Vec<Run<T>>,
        coalesce: bool,
    ) -> Self {
        let mut t = SignalTrack::empty(chrom, start, end);
        let mut prev_end = start;
        for run in runs {
            // Defensive: a run that does not advance the cursor carries no
            // value and is dropped. This is reachable from untrusted parsers,
            // so it must not be a panic.
            if run.end <= prev_end {
                continue;
            }
            let run_end = run.end.min(t.end);
            if run_end <= prev_end {
                continue;
            }
            if coalesce {
                if let Some(last) = t.runs.last_mut() {
                    if last.value == run.value {
                        last.end = run_end;
                        prev_end = run_end;
                        continue;
                    }
                }
            }
            t.runs.push(Run {
                end: run_end,
                value: run.value,
            });
            prev_end = run_end;
        }
        t
    }

    /// Build from upstream's legacy `[positions, values]` "pv array".
    ///
    /// # Value convention
    ///
    /// The pv array is **right-endpoint indexed**: `values[i]` is the value on
    /// `[positions[i-1], positions[i])`, with `positions[-1] := start`. This is
    /// fixed by [`MACS3/Signal/PileupV2.py::_write_pv_to_bedGraph`], which
    /// writes `fprintf("%s\t%d\t%d\t%.5f", chrom, pre, pos, value)` where
    /// `pre` is the *previous* position.
    ///
    /// So the first entry describes the region before the first read, and the
    /// last entry closes the covered region. Nothing is represented after
    /// `positions[n-1]`, and that tail is left uncovered.
    ///
    /// Getting this backwards silently shifts the whole signal by one interval,
    /// which is exactly the kind of error a value-only spot check misses.
    ///
    /// # Panics
    /// If `positions` is not strictly increasing, or the lengths disagree. Both
    /// would mean the oracle produced a malformed array.
    pub fn from_breakpoints(
        chrom: ChromId,
        start: Coord,
        end: Coord,
        positions: &[Coord],
        values: &[T],
    ) -> Self {
        assert_eq!(
            positions.len(),
            values.len(),
            "pv array: {} positions but {} values",
            positions.len(),
            values.len()
        );
        let mut runs = Vec::with_capacity(values.len());
        let mut prev_pos: Option<Coord> = None;
        for (i, (&p, &v)) in positions.iter().zip(values.iter()).enumerate() {
            if let Some(pp) = prev_pos {
                assert!(
                    p > pp,
                    "pv array positions must strictly increase at {i}: {p} after {pp}"
                );
            }
            prev_pos = Some(p);
            // `v` is the value on the region this breakpoint *closes*
            if p > start {
                runs.push(Run { end: p, value: v });
            }
        }
        SignalTrack::from_runs(chrom, start, end, runs)
    }

    /// Convert back to upstream's `[positions, values]` pv array.
    ///
    /// The exact inverse of [`SignalTrack::from_breakpoints`]: the run *ends* are
    /// the breakpoint positions, and each run's value is the value on the region
    /// that run closes.
    pub fn to_breakpoints(&self) -> (Vec<Coord>, Vec<T>) {
        let mut pos = Vec::with_capacity(self.runs.len());
        let mut val = Vec::with_capacity(self.runs.len());
        for run in &self.runs {
            if run.end <= self.start {
                continue;
            }
            pos.push(run.end);
            val.push(run.value);
        }
        (pos, val)
    }

    /// Append a run ending at `end`, coalescing with the last run when values
    /// match. A run at or before the cursor is ignored.
    pub fn push_run(&mut self, run: Run<T>) {
        if run.end <= self.cursor() {
            return;
        }
        if let Some(last) = self.runs.last_mut() {
            if last.value == run.value {
                last.end = run.end;
                return;
            }
        }
        self.runs.push(run);
    }

    /// Append `value` on `[cursor, end)`.
    #[inline]
    /// Append a run **without** coalescing an equal-valued predecessor.
    ///
    /// F69: upstream's pileup arrays keep an entry at *every* read start and
    /// read end, even when the depth does not change -- the captured paired
    /// array shows consecutive entries carrying identical depths. The ordinary
    /// [`Self::push`] coalesces them, which is right for bedGraph output and for
    /// the local-lambda merge but wrong for the track that feeds the score and
    /// chunk stages: dropping redundant endpoints shrinks the paired index space
    /// (10,684 upstream indices become 4,763 here) and therefore shrinks the
    /// `above_cutoff` candidate set (915 upstream chunks become 522).
    ///
    /// This is a distinct representation on purpose. Callers pick the one their
    /// consumer needs rather than the track guessing.
    pub fn push_exact(&mut self, end: Coord, value: T) {
        if end <= self.cursor() {
            return;
        }
        self.runs.push(Run { end, value });
    }

    pub fn push(&mut self, end: Coord, value: T) {
        self.push_run(Run { end, value });
    }

    /// Exclusive end of the last run, or `start` when empty. This is where the
    /// covered region stops.
    #[inline]
    pub fn cursor(&self) -> Coord {
        self.runs.last().map_or(self.start, |r| r.end)
    }

    /// The chromosome this track belongs to.
    pub fn chrom(&self) -> ChromId {
        self.chrom
    }

    /// Left edge of the track.
    pub fn start(&self) -> Coord {
        self.start
    }

    /// Right edge (exclusive) of the track's extent.
    pub fn end(&self) -> Coord {
        self.end
    }

    /// Number of runs.
    pub fn len(&self) -> usize {
        self.runs.len()
    }

    /// True when no base is covered.
    pub fn is_empty(&self) -> bool {
        self.runs.is_empty()
    }

    /// The raw runs.
    pub fn runs(&self) -> &[Run<T>] {
        &self.runs
    }

    /// Consume the track and yield its runs.
    pub fn into_runs(self) -> Vec<Run<T>> {
        self.runs
    }

    /// Iterate over all spans, resolving each run's start.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = Span<'_, T>> + '_ {
        let mut at = self.start;
        self.runs.iter().map(move |r| {
            let s = Span {
                start: at,
                end: r.end,
                value: &r.value,
            };
            at = r.end;
            s
        })
    }

    /// `(start, value)` of every run. The compact, canonical view.
    pub fn spans(&self) -> impl ExactSizeIterator<Item = (Coord, T)> + '_ {
        self.iter().map(|s| (s.start, *s.value))
    }

    /// Value at `pos`, or `None` if `pos` is outside the covered region.
    pub fn value_at(&self, pos: Coord) -> Option<T> {
        if pos < self.start || pos >= self.cursor() {
            return None;
        }
        self.runs.get(self.find_run_index(pos)).map(|r| r.value)
    }

    /// Value at `pos`, defaulting to `default` when `pos` is uncovered.
    pub fn value_at_or(&self, pos: Coord, default: T) -> T {
        self.value_at(pos).unwrap_or(default)
    }

    /// Index of the run containing `pos`. `pos` must be covered.
    #[inline]
    fn find_run_index(&self, pos: Coord) -> usize {
        // partition_point gives the first run whose end is > pos
        let i = self.runs.partition_point(|r| r.end <= pos);
        debug_assert!(i < self.runs.len(), "pos {pos} is not covered");
        i.min(self.runs.len().saturating_sub(1))
    }

    /// Index of the first run whose end is strictly greater than `pos`, or
    /// `runs.len()` when `pos` is past the covered region.
    #[inline]
    fn first_index_after(&self, pos: Coord) -> usize {
        self.runs.partition_point(|r| r.end <= pos)
    }

    /// The sub-track over `iv`, clamped to the track's extent.
    ///
    /// Runs are split at the interval edges, so the result has the same pointwise
    /// values restricted to `iv`.
    pub fn slice(&self, iv: Interval) -> SignalTrack<T> {
        let lo = iv.start().max(self.start);
        let hi = iv.end().min(self.end).max(lo);
        let mut out = SignalTrack::empty(self.chrom, lo, hi);
        let mut idx = self.first_index_after(lo);
        let mut at = lo;
        while at < hi && idx < self.runs.len() {
            let r = self.runs[idx];
            if r.end <= at {
                idx += 1;
                continue;
            }
            let stop = r.end.min(hi);
            if r.end > at {
                out.runs.push(Run {
                    end: stop,
                    value: r.value,
                });
            }
            at = stop;
            idx += 1;
        }
        out
    }

    /// Number of covered bases: `cursor() - start()`.
    pub fn covered_len(&self) -> Len {
        (self.cursor() - self.start).into()
    }

    /// `sum(f(value) * length)` over the covered region.
    pub fn integral_of<F: Fn(&T) -> f64>(&self, f: F) -> f64 {
        let mut at = self.start;
        let mut acc = 0.0_f64;
        for r in &self.runs {
            acc += f(&r.value) * (r.end - at) as f64;
            at = r.end;
        }
        acc
    }

    /// `sum(value * length)` over the covered region.
    pub fn integral(&self) -> f64
    where
        T: Copy + Into<f64>,
    {
        self.integral_of(|v| (*v).into())
    }

    /// Minimum value over all runs.
    pub fn min(&self) -> Option<T>
    where
        T: PartialOrd,
    {
        self.runs
            .iter()
            .map(|r| r.value)
            .min_by(|a, b| a.partial_cmp(b).expect("NaN in track value"))
    }

    /// Maximum value over all runs.
    pub fn max(&self) -> Option<T>
    where
        T: PartialOrd,
    {
        self.runs
            .iter()
            .map(|r| r.value)
            .max_by(|a, b| a.partial_cmp(b).expect("NaN in track value"))
    }

    /// Apply `f` to every value.
    pub fn map<U: Copy + PartialEq>(&self, f: impl Fn(T) -> U) -> SignalTrack<U> {
        SignalTrack {
            chrom: self.chrom,
            start: self.start,
            end: self.end,
            runs: self
                .runs
                .iter()
                .map(|r| Run {
                    end: r.end,
                    value: f(r.value),
                })
                .collect(),
        }
    }

    /// Multiply every value by `k`, in `f64`.
    ///
    /// The conversion happens *before* the multiply, so widening then scaling
    /// is what a caller gets; callers that need `f32` arithmetic should use
    /// [`SignalTrack::map`] with an explicit `f32` closure instead.
    pub fn scale(&self, k: f64) -> SignalTrack<f64>
    where
        T: Into<f64> + Copy,
    {
        let k = k as f32 as f64;
        self.map(|v| v.into() * k)
    }

    /// Histogram of value bits -> number of covered bases.
    ///
    /// This is exactly the quantity MACS3's p-score to q-score map is built
    /// from: for each distinct p-score, how many base pairs carry it. Keying by
    /// bits means `+0.0` and `-0.0` stay distinct and NaN payloads do not
    /// collapse, both of which a numeric key would get wrong.
    pub fn value_histogram(&self) -> HashMap<u64, u64>
    where
        T: BitKey,
    {
        let mut h: HashMap<u64, u64> = HashMap::with_capacity(self.runs.len());
        let mut at = self.start;
        for r in &self.runs {
            *h.entry(r.value.to_key()).or_insert(0) += u64::from(r.end - at);
            at = r.end;
        }
        h
    }

    /// The same histogram with the value recovered for each key, sorted by key
    /// (and therefore by value).
    pub fn value_histogram_sorted(&self) -> Vec<(T, u64)>
    where
        T: BitKey,
    {
        let mut out: Vec<(T, u64)> = self
            .value_histogram()
            .into_iter()
            .map(|(k, n)| (T::from_key(k), n))
            .collect();
        out.sort_by_key(|(v, _)| v.to_key());
        out
    }

    /// Pointwise combination of two tracks.
    ///
    /// This is the Rust form of upstream's `over_two_pv_array`
    /// ([`MACS3/Signal/PileupV2.py:950`]). The merge walks both breakpoint
    /// lists and **stops as soon as either is exhausted**, which is what
    /// upstream's `while i1 < l1 and i2 < l2` does. The last emitted value
    /// therefore persists to the end of the intersection of the two extents.
    ///
    /// This truncation is load-bearing: the `d`-window control track is strictly
    /// shorter than the `llocal`-window one, so `lambda_local` is only defined
    /// over the shorter of the two. Reproducing it is required for parity.
    ///
    /// Coalescing (which upstream does not do here) is safe: it merges adjacent
    /// runs with equal reduced values, leaving the pointwise value function
    /// identical.
    pub fn zip<U: Copy + PartialEq, V: Copy + PartialEq>(
        &self,
        other: &SignalTrack<U>,
        f: impl Fn(T, U) -> V,
    ) -> SignalTrack<V> {
        let lo = self.start.max(other.start);
        let hi = self.end.min(other.end).max(lo);
        let mut out = SignalTrack::empty(self.chrom, lo, hi);
        if lo >= hi {
            return out;
        }
        let mut i = self.first_index_after(lo);
        let mut j = other.first_index_after(lo);
        while i < self.runs.len() && j < other.runs.len() {
            let ra = self.runs[i];
            let rb = other.runs[j];
            let stop = ra.end.min(rb.end);
            if stop <= lo {
                // both breakpoints precede the window; advance and retry
                if ra.end <= rb.end {
                    i += 1;
                }
                if rb.end <= ra.end {
                    j += 1;
                }
                continue;
            }
            out.push(stop, f(ra.value, rb.value));
            if ra.end <= rb.end {
                i += 1;
            }
            if rb.end <= ra.end {
                j += 1;
            }
        }
        out
    }

    /// Pointwise maximum of an iterator of tracks.
    ///
    /// `tracks` must be non-empty. Because [`SignalTrack::zip`] stops at the
    /// shorter track, folding a set of tracks left-to-right is **not** the same
    /// as taking the pointwise max of the whole set. Callers that need the
    /// true pointwise max must first extend every track to a common coverage.
    /// Use [`SignalTrack::pad_to`] for that.
    ///
    /// # Panics
    /// If `tracks` is empty.
    pub fn pointwise_max<'a, I, Tr>(tracks: I) -> SignalTrack<T>
    where
        I: IntoIterator<Item = &'a Tr>,
        Tr: AsRef<SignalTrack<T>> + 'a,
        T: PartialOrd + 'a,
    {
        let mut it = tracks.into_iter();
        let Some(first) = it.next() else {
            panic!("pointwise_max requires at least one track");
        };
        let mut acc = first.as_ref().clone();
        for t in it {
            acc = acc.max_with(t.as_ref());
        }
        acc
    }

    /// Pointwise maximum of two tracks.
    pub fn max_with(&self, other: &SignalTrack<T>) -> SignalTrack<T>
    where
        T: PartialOrd,
    {
        self.zip(other, |a, b| if a > b { a } else { b })
    }

    /// A copy whose extent is exactly the covered region.
    ///
    /// [`SignalTrack::to_breakpoints`] is lossy about where coverage stops (the
    /// last pv value simply persists), so a pv round trip must be compared
    /// against `track.covered()`.
    pub fn covered(&self) -> SignalTrack<T> {
        SignalTrack {
            chrom: self.chrom,
            start: self.start,
            end: self.cursor(),
            runs: self.runs.clone(),
        }
    }

    /// Extend the covered region to `end` by repeating the last value.
    ///
    /// The right way to make a short track comparable with a long one.
    pub fn pad_to(&self, end: Coord) -> SignalTrack<T> {
        if end <= self.cursor() {
            return self.clone();
        }
        let mut out = self.clone();
        out.end = out.end.max(end);
        if let Some(last) = out.runs.last() {
            out.runs.push(Run {
                end,
                value: last.value,
            });
        }
        out
    }

    /// Restrict the extent to `[start, end)` without changing covered values.
    pub fn restrict_to(&self, start: Coord, end: Coord) -> SignalTrack<T> {
        let mut s = self.slice(Interval::new(start, end));
        s.end = end.max(start);
        s
    }

    /// Pointwise difference `self - other` over the overlap.
    pub fn sub<U: Copy + Into<f64> + PartialEq>(&self, other: &SignalTrack<U>) -> SignalTrack<f64>
    where
        T: Into<f64> + Copy,
    {
        self.zip(other, |a, b| a.into() - b.into())
    }

    /// The maximal intervals on which the predicate holds.
    pub fn intervals_where(&self, pred: impl Fn(T) -> bool) -> Vec<Interval> {
        self.iter()
            .filter(|s| pred(*s.value))
            .map(|s| Interval::new(s.start, s.end))
            .collect()
    }

    /// Map the predicate to two values, keeping the extent.
    pub fn select(&self, keep_value: T, drop_value: T, pred: impl Fn(T) -> bool) -> SignalTrack<T> {
        self.map(|v| if pred(v) { keep_value } else { drop_value })
    }

    /// Concatenate `other` after `self`, extending the extent.
    ///
    /// # Panics
    /// If the two tracks do not abut or overlap, or disagree on the overlap.
    pub fn concat(mut self, other: &SignalTrack<T>) -> SignalTrack<T> {
        assert!(
            other.start <= self.cursor(),
            "tracks must abut or overlap: {} > {}",
            other.start,
            self.cursor()
        );
        if other.start < self.cursor() {
            let overlap = self.slice(Interval::new(other.start, self.cursor()));
            for s in overlap.iter() {
                assert!(
                    *s.value == other.value_at_or(s.start, *s.value),
                    "tracks disagree on the overlap at {}",
                    s.start
                );
            }
            self.runs.retain(|r| r.end > other.start);
            if let Some(last) = self.runs.last_mut() {
                last.end = other.start;
            }
        }
        self.end = self.end.max(other.end);
        for r in &other.runs {
            self.push_run(*r);
        }
        self
    }

    /// True when both tracks have the same extent and identical runs.
    pub fn identical(&self, other: &SignalTrack<T>) -> bool {
        self.start == other.start && self.end == other.end && self.runs == other.runs
    }

    /// True when both tracks have the same extent and pointwise-equal values
    /// within `tol`.
    ///
    /// Merge-walks the two run lists, so differing run boundaries are not
    /// reported as a mismatch.
    pub fn approx_identical(&self, other: &SignalTrack<T>, tol: f64) -> bool
    where
        T: Copy + Into<f64> + PartialEq,
    {
        if self.start != other.start || self.end != other.end {
            return false;
        }
        let mut i = 0;
        let mut j = 0;
        let mut at = self.start;
        while i < self.runs.len() && j < other.runs.len() {
            let stop = self.runs[i].end.min(other.runs[j].end);
            let (va, vb) = (self.runs[i].value.into(), other.runs[j].value.into());
            if (va - vb).abs() > tol {
                return false;
            }
            at = stop;
            if self.runs[i].end == stop {
                i += 1;
            }
            if other.runs[j].end == stop {
                j += 1;
            }
        }
        at == self.cursor() && i == self.runs.len() && j == other.runs.len()
    }

    /// Sample every integer position in the covered region.
    ///
    /// For tests and small synthetic genomes only; a real chromosome is never
    /// materialised this way.
    pub fn dense(&self) -> impl Iterator<Item = (Coord, T)> + '_ {
        self.iter()
            .flat_map(|s| (s.start..s.end).map(move |p| (p, *s.value)))
    }
}

impl<T: Copy + PartialEq> Default for SignalTrack<T> {
    fn default() -> Self {
        SignalTrack {
            chrom: ChromId(0),
            start: 0,
            end: 0,
            runs: Vec::new(),
        }
    }
}

/// Incremental builder for a sweep-produced track.
///
/// The pileup engine writes into this: it never seeks, never sorts, and reuses
/// one allocation for a whole chromosome.
#[derive(Debug, Clone)]
pub struct TrackBuilder<T> {
    chrom: ChromId,
    start: Coord,
    end: Coord,
    runs: Vec<Run<T>>,
    cursor: Coord,
}

impl<T: Copy + PartialEq> TrackBuilder<T> {
    /// A builder over `[start, end)`.
    pub fn new(chrom: ChromId, start: Coord, end: Coord) -> Self {
        TrackBuilder::with_capacity(chrom, start, end, 1024)
    }

    /// A builder with a preallocated capacity.
    pub fn with_capacity(chrom: ChromId, start: Coord, end: Coord, cap: usize) -> Self {
        let end = end.max(start);
        TrackBuilder {
            chrom,
            start,
            end,
            runs: Vec::with_capacity(cap),
            cursor: start,
        }
    }

    /// Write `value` on `[cursor, end)`.
    ///
    /// # Panics (debug only)
    /// A non-monotone `end` is a caller bug: the sweep must produce a monotone
    /// event stream.
    pub fn push(&mut self, end: Coord, value: T) {
        debug_assert!(
            end >= self.cursor,
            "non-monotone push: {end} after {}",
            self.cursor
        );
        if end <= self.cursor {
            return;
        }
        let end = end.min(self.end);
        if end <= self.cursor {
            return;
        }
        if let Some(last) = self.runs.last_mut() {
            if last.value == value {
                last.end = end;
                self.cursor = end;
                return;
            }
        }
        self.runs.push(Run { end, value });
        self.cursor = end;
    }

    /// Append without coalescing an equal-valued predecessor. See
    /// [`SignalTrack::push_exact`] and F69.
    pub fn push_exact(&mut self, end: Coord, value: T) {
        if end <= self.cursor {
            return;
        }
        let end = end.min(self.end);
        if end <= self.cursor {
            return;
        }
        self.runs.push(Run { end, value });
        self.cursor = end;
    }

    /// Current write position.
    pub fn cursor(&self) -> Coord {
        self.cursor
    }

    /// Number of runs written so far.
    pub fn len(&self) -> usize {
        self.runs.len()
    }

    /// True when nothing has been written.
    pub fn is_empty(&self) -> bool {
        self.runs.is_empty()
    }

    /// Reuse the allocation for another run list.
    pub fn reset(&mut self) {
        self.runs.clear();
        self.cursor = self.start;
    }

    /// Finish and return the track. Any unwritten tail is left uncovered, which
    /// callers that need a total function combine with [`SignalTrack::pad_to`] or
    /// [`SignalTrack::constant`].
    pub fn finish(self) -> SignalTrack<T> {
        SignalTrack {
            chrom: self.chrom,
            start: self.start,
            end: self.end,
            runs: self.runs,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const C: ChromId = ChromId(0);

    fn track(runs: &[(u32, f32)], end: u32) -> SignalTrack<f32> {
        SignalTrack::from_runs(
            C,
            0,
            end,
            runs.iter().map(|&(e, v)| Run::new(e, v)).collect(),
        )
    }

    /// `(start, value)` of every run — the compact, canonical view.
    fn spans(t: &SignalTrack<f32>) -> Vec<(u32, f32)> {
        t.spans().collect()
    }

    /// `(position, value)` at every base — the dense view.
    fn dense(t: &SignalTrack<f32>) -> Vec<(u32, f32)> {
        t.dense().collect()
    }

    // ---------- construction and canonical form ----------

    #[test]
    fn empty_track_covers_nothing() {
        let t = SignalTrack::<f32>::empty(C, 100, 200);
        assert!(t.is_empty());
        assert_eq!(t.len(), 0);
        assert_eq!(t.cursor(), 100);
        assert_eq!(t.covered_len(), 0);
        assert_eq!(t.value_at(150), None);
        assert_eq!(t.value_at_or(150, -1.0), -1.0);
        assert_eq!(t.integral(), 0.0);
    }

    #[test]
    fn reversed_extent_is_normalised() {
        let t = SignalTrack::<f32>::empty(C, 200, 100);
        assert_eq!(t.start(), 200);
        assert_eq!(t.end(), 200);
    }

    #[test]
    fn adjacent_equal_values_coalesce() {
        // 1.0 on [0,20), 2.0 on [20,40)
        let t = track(&[(10, 1.0), (20, 1.0), (30, 2.0), (40, 2.0)], 40);
        assert_eq!(t.len(), 2);
        assert_eq!(spans(&t), vec![(0, 1.0), (20, 2.0)]);
    }

    #[test]
    fn runs_are_clamped_to_the_extent() {
        let t = SignalTrack::from_runs(C, 0, 30, vec![Run::new(10, 1.0), Run::new(99, 5.0)]);
        assert_eq!(t.end(), 30);
        assert_eq!(t.cursor(), 30);
        assert_eq!(spans(&t), vec![(0, 1.0), (10, 5.0)]);
    }

    #[test]
    fn zero_width_runs_are_dropped() {
        let t = SignalTrack::from_runs(
            C,
            0,
            100,
            vec![Run::new(10, 1.0), Run::new(10, 2.0), Run::new(20, 4.0)],
        );
        assert_eq!(spans(&t), vec![(0, 1.0), (10, 4.0)]);
    }

    #[test]
    fn backwards_runs_are_dropped_defensively() {
        // untrusted parsers must not be able to corrupt the value function
        let t = SignalTrack::from_runs(
            C,
            0,
            100,
            vec![Run::new(10, 1.0), Run::new(5, 3.0), Run::new(20, 4.0)],
        );
        assert_eq!(spans(&t), vec![(0, 1.0), (10, 4.0)]);
    }

    #[test]
    fn constant_track() {
        let t = SignalTrack::constant(C, 5, 15, 3.0);
        assert_eq!(t.len(), 1);
        assert_eq!(t.covered_len(), 10);
        assert_eq!(spans(&t), vec![(5, 3.0)]);
        assert!(SignalTrack::constant(C, 5, 5, 3.0).is_empty());
    }

    #[test]
    fn contiguity_holds() {
        // every covered base has a value, and runs never leave a hole
        let t = track(&[(3, 1.0), (9, 2.0), (10, 3.0), (25, 1.0)], 100);
        let mut at = t.start();
        for s in t.iter() {
            assert_eq!(s.start, at, "hole before run at {}", s.start);
            at = s.end;
        }
        assert_eq!(at, t.cursor());
        assert_eq!(t.cursor(), 25);
        // positions past the cursor are uncovered, not zero
        assert_eq!(t.value_at(24), Some(1.0));
        assert_eq!(t.value_at(25), None);
        assert_eq!(t.value_at(99), None);
    }

    // ---------- lookup ----------

    #[test]
    fn basic_lookup() {
        let t = track(&[(120, 1.0), (150, 2.0), (180, 1.0)], 180);
        assert_eq!(t.value_at(0), Some(1.0));
        assert_eq!(t.value_at(119), Some(1.0));
        assert_eq!(t.value_at(120), Some(2.0));
        assert_eq!(t.value_at(149), Some(2.0));
        assert_eq!(t.value_at(150), Some(1.0));
        assert_eq!(t.value_at(179), Some(1.0));
        assert_eq!(t.value_at(180), None);
        assert_eq!(t.value_at(u32::MAX), None);
    }

    #[test]
    fn value_at_or_falls_back_outside() {
        let t = track(&[(10, 1.0)], 10); // covers [0,10)
        assert_eq!(t.value_at_or(0, -1.0), 1.0);
        assert_eq!(t.value_at_or(9, -1.0), 1.0);
        assert_eq!(t.value_at_or(10, -1.0), -1.0);
        assert_eq!(t.value_at_or(500, -1.0), -1.0);
    }

    #[test]
    fn lookup_is_correct_for_a_long_track() {
        // a track with many runs must binary-search correctly at every boundary
        let runs: Vec<Run<f32>> = (1..500u32).map(|i| Run::new(i * 3, i as f32)).collect();
        let t = SignalTrack::from_runs(C, 0, 1500, runs);
        for i in 1..500u32 {
            let run_start = (i - 1) * 3;
            assert_eq!(t.value_at(run_start), Some(i as f32), "at {run_start}");
            assert_eq!(t.value_at(run_start + 1), Some(i as f32));
            assert_eq!(t.value_at(run_start + 2), Some(i as f32));
        }
    }

    // ---------- mass conservation ----------

    #[test]
    fn integral_conserves_mass() {
        // [0,120) at 1, [120,150) at 2, [150,180) at 1 => 120 + 60 + 30
        let t = track(&[(120, 1.0), (150, 2.0), (180, 1.0)], 180);
        assert_eq!(t.integral(), 210.0);
        assert_eq!(t.covered_len(), 180);
    }

    #[test]
    fn integral_matches_dense_sum() {
        let t = track(&[(7, 1.5), (13, 0.25), (14, 9.0)], 20);
        let d: f64 = dense(&t).iter().map(|&(_, v)| f64::from(v)).sum();
        assert!((t.integral() - d).abs() < 1e-9, "{} vs {d}", t.integral());
    }

    // ---------- upstream pv array interop ----------

    #[test]
    fn from_breakpoints_round_trips() {
        // runs [0,120)=1, [120,150)=2, [150,180)=1
        let t = track(&[(120, 1.0), (150, 2.0), (180, 1.0)], 200);
        let (pos, val) = t.to_breakpoints();
        // right-endpoint indexing: the run *ends* are the breakpoints
        assert_eq!(pos, vec![120, 150, 180]);
        assert_eq!(val, vec![1.0, 2.0, 1.0]);
        let back = SignalTrack::from_breakpoints(C, 0, 200, &pos, &val);
        assert!(back.identical(&t));
    }

    /// The pv format's first entry describes the region *before* the first read.
    #[test]
    fn from_breakpoints_reads_the_value_as_closing_the_previous_region() {
        // positions [5, 155], values [0.0, 1.0] means 0.0 on [0,5) and 1.0 on
        // [5,155)
        let t = SignalTrack::from_breakpoints(C, 0, 1000, &[5, 155], &[0.0f32, 1.0]);
        assert_eq!(t.value_at(4), Some(0.0));
        assert_eq!(t.value_at(5), Some(1.0));
        assert_eq!(t.value_at(154), Some(1.0));
        assert_eq!(
            t.value_at(155),
            None,
            "nothing is represented past the last bp"
        );
        assert_eq!(t.integral(), 150.0);
    }

    #[test]
    fn from_breakpoints_stops_at_the_last_breakpoint() {
        // there is no "persisting tail" in the pv format: the last entry closes
        // the covered region and nothing beyond it is represented
        let t = SignalTrack::from_breakpoints(C, 0, 1000, &[0, 10, 20], &[1.0f32, 2.0, 3.0]);
        assert_eq!(t.cursor(), 20);
        assert_eq!(t.value_at(19), Some(3.0));
        assert_eq!(t.value_at(20), None);
        assert_eq!(spans(&t), vec![(0, 2.0), (10, 3.0)]);
    }

    #[test]
    fn from_breakpoints_ignores_breakpoints_at_or_before_the_start() {
        let t = SignalTrack::from_breakpoints(C, 50, 1000, &[10, 60, 80], &[1.0f32, 2.0, 3.0]);
        assert_eq!(spans(&t), vec![(50, 2.0), (60, 3.0)]);
        assert_eq!(t.value_at(55), Some(2.0));
    }

    #[test]
    #[should_panic(expected = "strictly increase")]
    fn from_breakpoints_rejects_unsorted_input() {
        SignalTrack::from_breakpoints(C, 0, 100, &[0, 50, 50], &[1.0, 2.0, 3.0]);
    }

    #[test]
    #[should_panic(expected = "positions but")]
    fn from_breakpoints_rejects_length_mismatch() {
        SignalTrack::from_breakpoints(C, 0, 100, &[0, 50], &[1.0]);
    }

    // ---------- slicing ----------

    #[test]
    fn slice_splits_runs_and_preserves_values() {
        let t = track(&[(120, 1.0), (150, 2.0), (180, 1.0)], 200);
        let s = t.slice(Interval::new(100, 160));
        assert_eq!(s.start(), 100);
        assert_eq!(s.end(), 160);
        assert_eq!(spans(&s), vec![(100, 1.0), (120, 2.0), (150, 1.0)]);
        // the original is untouched
        assert!(t.identical(&track(&[(120, 1.0), (150, 2.0), (180, 1.0)], 200)));
    }

    #[test]
    fn slice_clamps_to_the_track() {
        let t = track(&[(50, 1.0)], 100);
        let s = t.slice(Interval::new(0, 1000));
        assert_eq!((s.start(), s.end()), (0, 100));
        let e = t.slice(Interval::new(500, 600));
        assert!(e.is_empty());
        assert_eq!((e.start(), e.end()), (500, 500));
    }

    // ---------- algebra ----------

    #[test]
    fn zip_is_pointwise() {
        let a = track(&[(100, 5.0), (300, 10.0)], 300);
        let b = track(&[(150, 3.0), (300, 7.0)], 300);
        let c = a.zip(&b, |x, y| x + y);
        assert_eq!(spans(&c), vec![(0, 8.0), (100, 13.0), (150, 17.0)]);
        assert_eq!(c.cursor(), 300);
    }

    #[test]
    fn zip_restricts_the_extent_to_the_overlap() {
        let a = track(&[(100, 5.0), (300, 10.0)], 300);
        let b = track(&[(150, 3.0), (1000, 7.0)], 1000);
        let c = a.zip(&b, |x, y| x + y);
        assert_eq!((c.start(), c.end()), (0, 300));
        // [0,100) 5+3, [100,150) 10+3, [150,300) 10+7
        assert_eq!(spans(&c), vec![(0, 8.0), (100, 13.0), (150, 17.0)]);
    }

    /// The load-bearing upstream quirk: `over_two_pv_array` stops as soon as
    /// either breakpoint list runs out, so the last value persists over the
    /// region where the shorter track has no breakpoints.
    #[test]
    fn zip_stops_at_the_shorter_breakpoint_list() {
        // a covers [0,300) with breaks at 100; b covers [0,200) with breaks at 50
        let a = track(&[(100, 5.0), (300, 10.0)], 300);
        let b = track(&[(50, 3.0), (200, 7.0)], 300);
        let c = a.zip(&b, |x, y| x + y);
        // b's breakpoint list runs out after 200, so the walk stops there
        assert_eq!(c.cursor(), 200);
        assert_eq!(spans(&c), vec![(0, 8.0), (50, 12.0), (100, 17.0)]);
        // the extent still covers [0,300), but the value function is only
        // *defined* up to the last breakpoint; `pad_to` makes the persisting
        // tail explicit
        assert_eq!(c.value_at(250), None);
        assert_eq!(c.end(), 300);
        assert_eq!(c.pad_to(300).value_at(250), Some(17.0));
    }

    #[test]
    fn zip_of_an_uncovered_track_is_empty() {
        let a = SignalTrack::<f32>::empty(C, 0, 10);
        let b = SignalTrack::constant(C, 0, 20, 1.0);
        assert!(a.zip(&b, |x, y| x + y).is_empty());
    }

    #[test]
    fn zip_of_disjoint_extents_is_empty() {
        let a = SignalTrack::constant(C, 0, 10, 1.0);
        let b = SignalTrack::constant(C, 50, 100, 1.0);
        assert!(a.zip(&b, |x, y| x + y).is_empty());
    }

    #[test]
    fn pad_to_makes_tracks_comparable() {
        let a = track(&[(100, 5.0), (300, 10.0)], 300);
        let b = track(&[(50, 3.0), (200, 7.0)], 300).pad_to(300);
        let m = a.max_with(&b);
        assert_eq!(m.cursor(), 300);
        // b padded to 300 is 3 on [0,50) and 7 on [50,300)
        // a is 5 on [0,100) and 10 on [100,300)
        assert_eq!(spans(&m), vec![(0, 5.0), (50, 7.0), (100, 10.0)]);
    }

    #[test]
    fn pointwise_max_over_many_tracks() {
        let a = track(&[(100, 5.0), (300, 1.0)], 300);
        let b = track(&[(50, 3.0), (300, 7.0)], 300);
        let c = track(&[(150, 4.0), (300, 2.0)], 300);
        let m = SignalTrack::pointwise_max([&a, &b, &c]);
        // a: 5 on [0,100), 1 after.  b: 3 on [0,50), 7 after.  c: 4 on [0,150), 2 after.
        // max: 5 on [0,50), 7 on [50,300)
        assert_eq!(m.cursor(), 300);
        assert_eq!(spans(&m), vec![(0, 5.0), (50, 7.0)]);
    }

    #[test]
    fn pointwise_max_equals_iterated_max_with() {
        let a = track(&[(37, 5.0), (300, 1.0)], 300);
        let b = track(&[(11, 3.0), (300, 7.0)], 300);
        let c = track(&[(150, 4.0), (300, 2.0)], 300);
        let iter = a.max_with(&b).max_with(&c);
        let one_shot = SignalTrack::pointwise_max([&a, &b, &c]);
        assert!(iter.identical(&one_shot));
    }

    #[test]
    fn sub_and_scale() {
        let a = track(&[(10, 5.0), (20, 3.0)], 20);
        let b = track(&[(5, 1.0), (20, 1.0)], 20);
        let d = a.sub(&b);
        assert_eq!(
            d.spans().collect::<Vec<_>>(),
            vec![(0u32, 4.0f64), (10, 2.0)],
            "a-b is 4 on [0,10) and 2 on [10,20)"
        );
        let s = a.scale(2.0);
        assert_eq!(
            s.spans().collect::<Vec<_>>(),
            vec![(0u32, 10.0f64), (10, 6.0)]
        );
    }

    // ---------- histogram ----------

    #[test]
    fn histogram_counts_covered_bases() {
        let t = track(&[(120, 1.0), (150, 2.0), (180, 1.0)], 180);
        let h = t.value_histogram();
        assert_eq!(h[&1.0f32.to_key()], 150); // [0,120) + [150,180)
        assert_eq!(h[&2.0f32.to_key()], 30);
        assert_eq!(h.values().sum::<u64>(), 180);
        assert_eq!(
            t.value_histogram_sorted(),
            vec![(1.0f32, 150), (2.0f32, 30)]
        );
    }

    #[test]
    fn histogram_ignores_the_uncovered_tail() {
        let t = track(&[(10, 1.0)], 1000);
        assert_eq!(t.value_histogram_sorted(), vec![(1.0f32, 10)]);
    }

    #[test]
    fn bit_key_distinguishes_signed_zero() {
        assert_ne!(0.0f32.to_key(), (-0.0f32).to_key());
        assert_ne!(0.0f64.to_key(), (-0.0f64).to_key());
        assert_eq!(
            f32::from_key((-0.0f32).to_key()).to_bits(),
            (-0.0f32).to_bits()
        );
        assert_eq!(
            f64::from_key(f64::NAN.to_key()).to_bits(),
            f64::NAN.to_bits()
        );
        assert_eq!(i32::from_key((-5i32).to_key()), -5);
        assert_eq!(u64::from_key(7u64.to_key()), 7);
    }

    #[test]
    fn bit_key_is_order_preserving() {
        // values are listed in ascending numeric order; the keys must be too
        let vals: Vec<f32> = vec![
            -1e30, -100.0, -1.0, -0.5, -1e-30, -0.0, 0.0, 1e-30, 0.5, 1.0, 100.0, 1e30,
        ];
        let keys: Vec<u64> = vals.iter().map(|v| v.to_key()).collect();
        let mut sorted = keys.clone();
        sorted.sort_unstable();
        assert_eq!(keys, sorted, "f32 keys must sort like the values");

        let vs: Vec<f64> = vec![-1e300, -1.0, -0.0, 0.0, f64::MIN_POSITIVE, 1.0, 1e300];
        let ks: Vec<u64> = vs.iter().map(|v| v.to_key()).collect();
        let mut ss = ks.clone();
        ss.sort_unstable();
        assert_eq!(ks, ss, "f64 keys must sort like the values");

        // and a shuffle-then-sort-by-key round trip restores numeric order
        let mut shuffled = vals.clone();
        shuffled.reverse();
        shuffled.sort_by_key(|v| v.to_key());
        assert_eq!(shuffled, vals);
        let mut shuffled64 = vs.clone();
        shuffled64.reverse();
        shuffled64.sort_by_key(|v| v.to_key());
        assert_eq!(shuffled64, vs);
    }

    // ---------- builder ----------

    #[test]
    fn builder_coalesces_and_ignores_stalls() {
        let mut b = TrackBuilder::new(C, 0, 100);
        b.push(10, 1.0);
        b.push(20, 1.0); // coalesced: 1.0 now covers [0,20)
        b.push(20, 9.0); // zero width, ignored
        b.push(30, 2.0);
        let t = b.finish();
        assert_eq!(t.len(), 2);
        assert_eq!(spans(&t), vec![(0, 1.0), (20, 2.0)]);
    }

    #[test]
    fn builder_clamps_to_extent() {
        let mut b = TrackBuilder::new(C, 0, 50);
        b.push(100, 7.0);
        let t = b.finish();
        assert_eq!(t.end(), 50);
        assert_eq!(spans(&t), vec![(0, 7.0)]);
        assert_eq!(t.covered_len(), 50);
    }

    #[test]
    fn builder_reset_reuses_the_allocation() {
        let mut b = TrackBuilder::with_capacity(C, 0, 10, 16);
        b.push(5, 1.0);
        b.reset();
        assert!(b.is_empty());
        assert_eq!(b.cursor(), 0);
        b.push(5, 2.0);
        assert_eq!(spans(&b.finish()), vec![(0, 2.0)]);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "non-monotone")]
    fn builder_rejects_backwards_push() {
        let mut b = TrackBuilder::new(C, 0, 100);
        b.push(30, 1.0);
        b.push(25, 1.0);
    }

    // ---------- concat ----------

    #[test]
    fn concat_extends() {
        let a = SignalTrack::constant(C, 0, 100, 2.0);
        let b = SignalTrack::from_runs(C, 0, 50, vec![Run::new(10, 2.0)]);
        let c = a.concat(&b);
        assert_eq!(c.end(), 100);
        assert_eq!(spans(&c), vec![(0, 2.0)]);
    }

    #[test]
    #[should_panic(expected = "tracks disagree")]
    fn concat_rejects_inconsistent_overlap() {
        let a = SignalTrack::constant(C, 0, 100, 2.0);
        let b = SignalTrack::constant(C, 0, 50, 5.0);
        let _ = a.concat(&b);
    }

    #[test]
    #[should_panic(expected = "must abut")]
    fn concat_rejects_a_hole() {
        let a = SignalTrack::constant(C, 0, 10, 1.0);
        let b = SignalTrack::constant(C, 50, 100, 1.0);
        let _ = a.concat(&b);
    }

    // ---------- comparison ----------

    #[test]
    fn approx_identical_reports_value_differences() {
        let a = track(&[(100, 1.0), (200, 2.0)], 200);
        let b = track(&[(100, 1.0), (200, 2.5)], 200);
        assert!(!a.identical(&b));
        assert!(!a.approx_identical(&b, 0.0));
        assert!(!a.approx_identical(&b, 0.4));
        assert!(a.approx_identical(&b, 0.5));
    }

    #[test]
    fn approx_identical_respects_extent() {
        let a = track(&[(100, 1.0)], 100);
        let b = track(&[(100, 1.0)], 200);
        assert!(!a.approx_identical(&b, 0.0));
    }

    #[test]
    fn identical_and_approx_identical_agree_on_equals() {
        let a = track(&[(10, 1.0), (20, 2.0)], 20);
        let b = a.clone();
        assert!(a.identical(&b));
        assert!(a.approx_identical(&b, 0.0));
    }

    // ---------- selection ----------

    #[test]
    fn intervals_where_finds_blocks() {
        let t = track(&[(10, 0.0), (20, 5.0), (30, 0.0), (50, 7.0)], 60);
        let ivs = t.intervals_where(|v| v >= 5.0);
        assert_eq!(ivs, vec![Interval::new(10, 20), Interval::new(30, 50)]);
    }

    #[test]
    fn select_maps_predicate_to_two_values() {
        let t = track(&[(10, 0.0), (20, 5.0)], 20);
        let s = t.select(1.0, 0.0, |v| v >= 5.0);
        assert_eq!(spans(&s), vec![(0, 0.0), (10, 1.0)]);
    }

    #[test]
    fn min_and_max() {
        let t = track(&[(10, 3.0), (20, 1.0), (30, 7.0)], 30);
        assert_eq!(t.min(), Some(1.0));
        assert_eq!(t.max(), Some(7.0));
        assert_eq!(SignalTrack::<f32>::empty(C, 0, 1).min(), None);
    }
}
