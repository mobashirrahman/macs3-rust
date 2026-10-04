//! Down-sampling of single-end positions.
//!
//! Upstream's `FWTrack.sample_percent` and `sample_num` seed NumPy's **global**
//! Mersenne-Twister RNG, shuffle each per-chromosome strand array in place,
//! truncate, and re-sort. Reproducing the retained *set* exactly therefore
//! requires NumPy's `shuffle` stream, not merely a shuffle with the same
//! distribution.
//!
//! NumPy's legacy `RandomState.shuffle` is Fisher-Yates from the end, consuming
//! one bounded integer per step via Lemire-style masked rejection:
//!
//! ```text
//! for i = n-1 down to 1:
//!     j = rng_interval(i + 1)          # uniform in [0, i]
//!     swap(a[i], a[j])
//! ```
//!
//! [`macs_stats::NumpyRng`] implements that stream. The samplers here are kept
//! separate from the storage types so the RNG dependency is explicit and testable.

use macs_core::Result;

use crate::frag::FragmentTrack;
use crate::SingleEndTrack;

/// Retain `percent` of positions on each strand, per chromosome.
///
/// `num = round(len * percent, 5)`, matching Python's two-argument `round`, which
/// is banker's rounding — `round(0.5)` is `0`, not `1`. Getting that wrong shifts
/// the retained count by one on exact halves.
///
/// Returns the new total.
pub fn sample_percent(track: &mut SingleEndTrack, percent: f64, seed: i64) -> Result<u64> {
    if !track.is_sorted() {
        track.positions_mut().finalize();
    }
    // sorted byte order, not file order: see `chroms_sorted` and F29
    let chroms = track.positions().chroms_sorted();
    // Legacy MT19937 on purpose. `FWTrack.sample_percent` (`FWTrack.py:632-649`) goes
    // through `np.random.seed(seed)` / `np.random.shuffle` -- the *old* `RandomState`
    // API, not `default_rng`. `randsample` and `filterdup` both come through here and
    // are byte-identical against the oracle with this stream; switching it to PCG64
    // breaks them. `sample_frag_percent` below is the one that needs PCG64 (F206).
    let mut rng = macs_stats::NumpyRng::seeded(seed);
    let mut total: u64 = 0;
    for chrom in chroms {
        for strand in [macs_core::Strand::Plus, macs_core::Strand::Minus] {
            let num = {
                let slice = track.positions().strand(chrom, strand);
                retained_count(slice.len() as f64 * percent)
            };
            // copy out, shuffle, truncate, sort, write back: a borrow per phase
            // rather than one long &mut
            let mut copy = track.positions().strand(chrom, strand).to_vec();
            // The shuffle happens **unconditionally**, including when `num` is 0.
            // It consumes a variable number of words from the stream, so skipping
            // it for an empty or tiny quota desynchronises every later draw and
            // the whole run diverges. And the write-back is unconditional too:
            // upstream resizes to `num` every time, so percent 0 empties the track
            // rather than leaving it untouched.
            rng.shuffle(&mut copy);
            // `resize(num)` semantics, not `truncate(num)`: a num larger than the
            // array zero-fills rather than clamping (F30)
            if copy.len() < num {
                copy.resize(num, 0);
            } else {
                copy.truncate(num);
            }
            copy.sort_unstable();
            *track.positions_mut().strand_mut(chrom, strand) = copy;
            total += num as u64;
        }
    }
    Ok(total)
}

/// Retain approximately `samplesize` positions.
///
/// Upstream's `sample_num` is not a separate algorithm: it converts the requested
/// count into a fraction and delegates.
///
/// ```text
/// percent = float(samplesize) / self.total
/// sample_percent(percent, seed)
/// ```
///
/// # A request above the track size zero-pads instead of clamping (F30)
///
/// When `samplesize > total`, `percent > 1` and the retained length computed from
/// it exceeds the array. NumPy's `ndarray.resize` to a larger size does not clamp,
/// it **zero-fills**, so the track ends up holding fabricated positions at 0 and
/// `total` reports roughly `samplesize`. Measured on MACS3 3.0.5:
///
/// ```text
/// before total: 37 ; sample_num(1000) -> after total: 999
/// plus: [0, 0, 0, 0, 0, 0, 0, 0] ...
/// ```
///
/// Those zeros are then pileup'd as real reads. [`sample_percent`] reproduces it
/// with `resize`-style zero-fill rather than `truncate`.
pub fn sample_num(track: &mut SingleEndTrack, samplesize: u64, seed: i64) -> Result<u64> {
    let total = track.total();
    if total == 0 {
        return Ok(0);
    }
    // `percent` is declared `cython.float` upstream -- a C `float`, so the
    // quotient is computed in f64 and then **truncated to f32** before use.
    //
    // That is observable. With 15 positions and a target of 1000, f64 gives
    // 1000/15 = 66.66666666666667 and the two strands truncate to 600 and 400,
    // totalling 1000. As an f32 the quotient is 66.66666412353516, which yields
    // 599.99997... and 399.99998..., truncating to 599 and 399 -- total 998.
    // Measured against MACS3 3.0.5, which reports 998.
    let percent = (samplesize as f64 / total as f64) as f32 as f64;
    sample_percent(track, percent, seed)
}

/// The retained count for one array: `int(round(len * percent, 5))`.
///
/// Upstream writes this as `cython.cast(cython.int, round(n * percent, 5))`, which
/// is **two** steps and they do different things:
///
/// * `round(x, 5)` is Python's two-argument round — round-half-to-even on the
///   5-decimal string;
/// * the `cython.int` cast then *truncates* toward zero.
///
/// So it is not a rounding step at the end. `round(0.5, 5)` is `0.5` and the cast
/// makes it `0`, not `1`; and `round(2.7, 5)` is `2.7` and the cast makes it `2`,
/// not `3`. Rounding instead of truncating would keep one extra position per array
/// and desynchronise the whole sampling stream that follows.
/// `uint(round(v, 5))`: round to five decimals, then truncate toward zero.
///
/// Upstream sizes every sampling quota this way (`PairedEndTrack.py:605`,
/// `FWTrack.py` likewise).
pub fn retained_count(v: f64) -> usize {
    if !v.is_finite() || v <= 0.0 {
        return 0;
    }
    let rounded = macs_stats::py_round(v, 5);
    if rounded <= 0.0 {
        return 0;
    }
    rounded.trunc() as usize
}

/// `PETrackII.sample_percent_copy(percent, seed)` (`PairedEndTrack.py:1958`),
/// then `fraglengths()` and the `[min_fraglen, max_fraglen]` window that
/// `HMMR_EM.__init__` applies.
///
/// # A counted track is sampled by **count**, not by row
///
/// The target is `n_sample = round(sum(counts) * percent)` *fragments*, and the
/// sampling unit is the count-expanded row index:
///
/// ```text
/// idx_flat = np.repeat(np.arange(len(loc)), counts)
/// shuffle(idx_flat)
/// idx_flat = idx_flat[:n_sample]
/// unique_idx, new_counts = np.unique(idx_flat, return_counts=True)
/// ```
///
/// So a row with count 7 is entered 7 times and may be drawn several times; the
/// result keeps one row per surviving index with the drawn multiplicity. Sampling
/// rows instead would keep 10% of the *rows*, which is a different number of
/// fragments whenever counts vary -- on the `mfrag` fixture it keeps 125 rows
/// carrying 160 fragments, not 40 rows.
///
/// # The RNG is a `SeedSequence`-keyed MT19937, and this port does not match it
///
/// Upstream uses `RandomState(MT19937(SeedSequence(seed)))`. This port's
/// [`macs_stats::randomstate_from_seed_sequence`] reproduces the
/// `SeedSequence` key -- verified against the installed NumPy -- but its
/// `MT19937` array-seeding stage does not yet reproduce NumPy's, so the
/// **permutation differs**. That is why this sampler is not claimed exact: with
/// `--no-fragem` the EM is skipped and `means`/`stddevs` come straight from the
/// user, which *is* exact; without it the EM sees a different 10% down-sample.
/// See F144 in `docs/upstream-findings.md`.
///
/// # The two track classes sample from different generators
///
/// `PETrackI.sample_percent_copy` uses `np.random.default_rng` (PCG64), while
/// `PETrackII` uses `SeedSequence`-keyed MT19937. Same seed, different stream, so
/// the two cannot share an implementation.
///
/// # A chromosome that samples nothing does not advance the stream
///
/// The `n == 0 or n_sample == 0` guard comes *before* the shuffle, so a
/// chromosome whose target rounds to zero is skipped without consuming any
/// randomness. `np.unique` returns sorted indices, so the surviving rows come
/// back in original order regardless of the shuffle.
///
/// Returns the filtered fragment lengths, i.e. what EM consumes.
pub fn sample_frag_percent(
    track: &FragmentTrack,
    percent: f64,
    seed: i64,
    min_fraglen: macs_core::Coord,
    max_fraglen: macs_core::Coord,
) -> Result<Vec<macs_core::Coord>> {
    let genome = track.genome();
    let mut chroms = track.chroms();
    chroms.sort_by(|&a, &b| genome.name(a).cmp(genome.name(b)));
    // F206: `PETrackI.sample_percent_copy` (`PairedEndTrack.py:653-659`) seeds
    // `np.random.default_rng(seed)` -- **PCG64** -- so this is the one sampler that must
    // not use the legacy MT19937 stream. They are not interchangeable: with the legacy
    // stream the EM down-sample retained 1041 fragments where upstream retains 1040,
    // which shifts HMMRATAC's fitted nucleosome means by a few bp and so moves every
    // accessible region derived from the digested signals.
    //
    // This is *not* `np.random.seed(seed)`; the two APIs are unrelated even though both
    // are "seeded by an integer". See `macs_stats::pcg64` for the three details that
    // make the permutation match.
    let mut rng = macs_stats::Pcg64Rng::seeded(seed as u128);
    let mut out: Vec<macs_core::Coord> = Vec::new();
    for c in chroms {
        let frags = track.frags(c);
        let counts = track.counts(c);
        let n: u64 = counts.iter().map(|k| u64::from(*k)).sum();
        // F206: `PETrackI.sample_percent_copy` computes
        //
        //     num = cython.cast(cython.uint, round(loc.shape[0] * percent, 5))
        //
        // which is **not** `round(n * percent)`. `round(x, 5)` keeps a *float* --
        // it rounds to five decimal places -- and the subsequent Cython cast to
        // `uint` truncates toward zero. So `11859 * 0.1` becomes 1185.9 and then
        // 1185, not the 1186 that integer rounding would give.
        //
        // That single fragment is the whole residual of the HMMRATAC divergence: the
        // EM then saw 1040 fragments instead of 1041, which moved the fitted
        // nucleosome means (mono 182.6 -> 186.3 against upstream's 186.3) and, through
        // the weight mapping and digested signals, every accessible region.
        let n_sample = if n == 0 {
            0
        } else {
            let rounded5 = macs_stats::py_round(n as f64 * percent, 5);
            if rounded5 <= 0.0 {
                0
            } else {
                rounded5 as u64
            }
        };
        if n == 0 || n_sample == 0 {
            continue;
        }
        // expand to one entry per counted fragment, then shuffle and truncate
        let mut idx_flat: Vec<usize> = Vec::with_capacity(n as usize);
        for (i, k) in counts.iter().enumerate() {
            for _ in 0..*k {
                idx_flat.push(i);
            }
        }
        rng.shuffle(&mut idx_flat);
        idx_flat.truncate(n_sample as usize);
        // `np.unique(..., return_counts=True)`: sorted unique values and their
        // multiplicities.
        idx_flat.sort_unstable();
        let mut kept: Vec<(usize, u32)> = Vec::new();
        for i in idx_flat {
            match kept.last_mut() {
                Some((last, c)) if *last == i => *c += 1,
                _ => kept.push((i, 1)),
            }
        }
        // `fraglengths()` expands each surviving row by its *new* count, so a
        // row drawn three times contributes its length three times.
        for (i, drawn) in kept {
            let f = frags[i];
            out.extend(std::iter::repeat_n(
                f.end.saturating_sub(f.start),
                drawn as usize,
            ));
        }
    }
    out.retain(|&l| l >= min_fraglen && l <= max_fraglen);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn retained_count_truncates_after_rounding_to_five_decimals() {
        // round(x, 5) then a truncating int cast -- so these are truncations,
        // not round-half-even outcomes
        assert_eq!(retained_count(0.5), 0, "round(0.5,5)=0.5, cast to int is 0");
        assert_eq!(retained_count(1.5), 1);
        assert_eq!(retained_count(2.7), 2, "round(2.7,5)=2.7, cast to int is 2");
        assert_eq!(retained_count(2.700001), 2, "still truncates");
        assert_eq!(retained_count(2.9), 2);
        assert_eq!(retained_count(3.0), 3);
    }

    #[test]
    fn negative_and_non_finite_are_zero() {
        assert_eq!(retained_count(-1.0), 0);
        assert_eq!(retained_count(f64::NAN), 0);
    }
}
