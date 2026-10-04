//! `numpy.random.default_rng` -- the PCG64 **Generator** API.
//!
//! Why this exists: `PETrackI.sample_percent_copy` (`PairedEndTrack.py:655-659`)
//! seeds a *modern* NumPy generator, not the legacy one:
//!
//! ```python
//! if seed >= 0:
//!     rs = np.random.default_rng(seed)
//! else:
//!     rs = np.random.default_rng()
//! rs_shuffle = rs.shuffle
//! ```
//!
//! [`crate::numpy_rng::NumpyRng`] is the legacy `RandomState` (MT19937 + Gauss) stream
//! that upstream's `randsample`/`filterdup` shuffle uses, and it is bit-exact against
//! NumPy. But `default_rng` is PCG64, seeded through `SeedSequence`, and its
//! `Generator.shuffle` draws through a *different* bounded-integer routine
//! (`bounded_lemire_uint64`, not `RandomState`'s masked rejection). The two are
//! unrelated algorithms, so a legacy stream cannot stand in for it.
//!
//! The symptom was small and easy to misread: upstream's fragment down-sample for
//! HMMRATAC's EM retained **1040** fragments where this port retained **1041**. One
//! fragment moves the fitted nucleosome means by several bp, which moves the digested
//! signals, which moves every decoded posterior -- 13 of 946 accessible regions
//! missing. Everything downstream was faithful; the shuffle was not.
//!
//! This is the *only* place upstream uses `default_rng` (2 call sites, both this
//! sampling path), so it is a bounded piece of work.
//!
//! # Scope
//!
//! [`Pcg64`] is complete and verified bit-for-bit against the installed NumPy: the
//! `SeedSequence` expansion, the 128-bit seeding constants, and the XSL-RR output
//! permutation.
//!
//! [`Pcg64Rng::interval`] and [`Pcg64Rng::shuffle`] transcribe
//! `numpy/random/src/distributions/distributions.c` and `_generator.pyx` directly.
//! Two details are load-bearing and neither is guessable:
//!
//! * `random_interval` draws from **`next_uint32`** whenever `max <= 0xffffffff`, which
//!   is every case a shuffle hits, and the bound is **inclusive**.
//! * `next_uint32` is the **high** half of the 64-bit XSL-RR output, not a second LCG
//!   step and not the low half.

/// PCG64 (XSL-RR, 128-bit state, 128-bit increment) as NumPy instantiates it.
///
/// Note this is plain PCG64, **not** `PCG64DXSM`: `default_rng` selects `PCG64` and the
/// DXSM variant is only reachable by asking for it explicitly.
#[derive(Debug, Clone)]
pub struct Pcg64 {
    state: u128,
    inc: u128,
}

/// `PCG_DEFAULT_MULTIPLIER_128` from `pcg64.h`.
const MULTIPLIER_128: u128 = 0x2360_ED05_1FC6_5DA4_4385_DF64_9FCC_F645;

/// `pcg64_step_r`: one LCG advance, `state = state * MULT + inc`.
fn step_r(state: u128, inc: u128) -> u128 {
    state.wrapping_mul(MULTIPLIER_128).wrapping_add(inc)
}

/// `pcg64_set_seed(state, seed, inc)` from `pcg64.h`.
///
/// Note `inc` arrives already shifted and OR'd with 1 by NumPy before it gets here;
/// `numpy/random/_pcg64.pyx` computes `inc = (seq << 1) | 1` and hands the raw `seq`
/// here, so that shift is applied below.
fn pcg64_set_seed(initstate: u128, initseq: u128) -> (u128, u128) {
    let inc = (initseq << 1) | 1;
    let mut state = 0u128;
    state = step_r(state, inc);
    state = state.wrapping_add(initstate);
    state = step_r(state, inc);
    (state, inc)
}

impl Pcg64 {
    /// `np.random.default_rng(seed)`: `SeedSequence(seed)` -> 4 u64 words -> state/inc.
    pub fn from_seed(seed: u128) -> Self {
        let seq = crate::seed_sequence::SeedSequence::new(seed as u64);
        // `SeedSequence.generate_state(4, uint64)`: eight u32 words, packed
        // little-endian into four u64s.
        let w32 = seq.generate_state(8);
        let w: Vec<u64> = w32
            .chunks(2)
            .map(|c| (c[0] as u64) | ((c[1] as u64) << 32))
            .collect();
        // PCG's `pcg128_t` takes the **first** word as the high half, so the pairing is
        // big-endian: `initstate = w[0]<<64 | w[1]`, `initseq = w[2]<<64 | w[3]`.
        // Pairing little-endian produces a plausible-looking but different 128-bit
        // constant, which then yields a completely different stream.
        let initstate = ((w[0] as u128) << 64) | (w[1] as u128);
        let initseq = ((w[2] as u128) << 64) | (w[3] as u128);
        let (state, inc) = pcg64_set_seed(initstate, initseq);
        Self { state, inc }
    }

    /// The current 128-bit state, for tests.
    pub fn state(&self) -> u128 {
        self.state
    }

    /// The current 128-bit increment, for tests.
    pub fn inc(&self) -> u128 {
        self.inc
    }

    /// `next_uint64`: `pcg_output_xsl_rr_128_64` applied *after* the LCG step.
    ///
    /// The order matters. `pcg_setseq_128_xsl_rr_64_random_r` steps first and permutes
    /// the **new** state; emitting the permutation of the pre-step state shifts the whole
    /// stream by one.
    pub fn next_u64(&mut self) -> u64 {
        self.state = step_r(self.state, self.inc);
        let hi = (self.state >> 64) as u64;
        let lo = self.state as u64;
        let xorshifted = hi ^ lo;
        let rot = (self.state >> 122) as u32;
        xorshifted.rotate_right(rot)
    }
}

/// `gen_mask` from `numpy/random/src/distributions/distributions.c`: the smallest
/// bitmask that is `>= max`.
fn gen_mask(mut max: u64) -> u64 {
    max |= max >> 1;
    max |= max >> 2;
    max |= max >> 4;
    max |= max >> 8;
    max |= max >> 16;
    max |= max >> 32;
    max
}

/// `numpy.random.Generator`: PCG64 plus the bounded draw its methods share.
#[derive(Debug, Clone)]
pub struct Pcg64Rng {
    bitgen: Pcg64,
    /// The high half of the last 64-bit draw, pending as the next 32-bit value.
    u32_buf: Option<u32>,
}

impl Pcg64Rng {
    /// `np.random.default_rng(seed)`.
    pub fn seeded(seed: u128) -> Self {
        Self {
            bitgen: Pcg64::from_seed(seed),
            u32_buf: None,
        }
    }

    /// `bit_generator.random_raw()`.
    pub fn next_u64(&mut self) -> u64 {
        self.bitgen.next_u64()
    }

    /// `random_interval(bitgen, max)` verbatim from `distributions.c`.
    ///
    /// ```c
    /// uint64_t random_interval(bitgen_t *bitgen_state, uint64_t max) {
    ///   uint64_t mask, value;
    ///   if (max == 0) return 0;
    ///   mask = max;
    ///   mask |= mask >> 1;  mask |= mask >> 2;  mask |= mask >> 4;
    ///   mask |= mask >> 8;  mask |= mask >> 16; mask |= mask >> 32;
    ///   if (max <= 0xffffffffUL) {
    ///     while ((value = (next_uint32(bitgen_state) & mask)) > max) ;
    ///   } else {
    ///     while ((value = (next_uint64(bitgen_state) & mask)) > max) ;
    ///   }
    ///   return value;
    /// }
    /// ```
    ///
    /// The `max <= 0xffffffff` branch is the one every shuffle iteration takes: the
    /// loop bounds are element indices, so `max` is at most `n-1`. Getting that wrong
    /// -- using the 64-bit draw here, or Lemire's `bounded_lemire_uint64` that
    /// `random_bounded_uint64` uses elsewhere -- silently yields a different
    /// permutation while still looking like a correct shuffle.
    pub fn interval(&mut self, max: u64) -> u64 {
        if max == 0 {
            return 0;
        }
        let mask = gen_mask(max);
        if max <= 0xffff_ffff {
            // The C code masks the 32-bit draw with the *64-bit* mask; for
            // `max <= 0xffffffff` the mask fits in 32 bits, so the results agree.
            let mask32 = mask as u32;
            loop {
                let v = (self.next_uint32() & mask32) as u64;
                if v <= max {
                    return v;
                }
            }
        }
        loop {
            let v = self.next_u64() & mask;
            if v <= max {
                return v;
            }
        }
    }

    /// `Generator.shuffle` on a 1-D slice, via `_shuffle_raw`
    /// (`_generator.pyx`): backward Fisher-Yates from `first = 1`, so the whole slice
    /// is shuffled, with each swap partner drawn by [`Pcg64Rng::interval`].
    pub fn shuffle<T>(&mut self, data: &mut [T]) {
        let n = data.len();
        if n < 2 {
            return;
        }
        for i in (1..n).rev() {
            let j = self.interval(i as u64) as usize;
            data.swap(i, j);
        }
    }

    /// `bitgen.next_uint32` for PCG64: a **32-bit buffer** over the 64-bit stream,
    /// low half first, then high half, then a fresh draw.
    ///
    /// This is the detail that cost the most to pin down, because every plausible
    /// alternative looks reasonable:
    ///
    /// * a second LCG step through PCG's `pcg_output_xsl_rr_128_32` -- wrong;
    /// * the low half only, or the high half only -- wrong;
    /// * one step per 32-bit value -- wrong, that skips a draw and desynchronises.
    ///
    /// Recoverable exactly from Python, because `rng.integers(0, 2**32, dtype=uint32)`
    /// hits `random_bounded_uint32`'s `rng == 0xFFFFFFFF` fast path, which returns
    /// `next_uint32` unfiltered. Over seeds 0, 1, 7, 999, 10151 and 12345 that
    /// sequence is reproduced exactly by
    ///
    /// ```text
    /// step; yield lo32; yield hi32; step; yield lo32; yield hi32; ...
    /// ```
    ///
    /// Note `random(dtype=float32)` is *not* a usable probe: its values are
    /// `(next_uint32 >> 8) * 2**-24`, so only the top 24 bits are observable, and
    /// several wrong variants happen to agree there.
    fn next_uint32(&mut self) -> u32 {
        if let Some(v) = self.u32_buf.take() {
            return v;
        }
        let v = self.bitgen.next_u64();
        let lo = v as u32;
        self.u32_buf = Some((v >> 32) as u32);
        lo
    }

    /// The raw 32-bit stream, for the test that pins it against NumPy.
    #[cfg(test)]
    pub(crate) fn next_uint32_for_test(&mut self) -> u32 {
        self.next_uint32()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// States captured from `numpy.random.default_rng(seed).bit_generator.state`.
    #[test]
    fn seeding_matches_numpy_default_rng() {
        // (seed, state, inc)
        let cases: [(u128, u128, u128); 4] = [
            (
                0,
                35399562948360463058890781895381311971,
                87136372517582989555478159403783844777,
            ),
            (
                1,
                207833532711051698738587646355624148094,
                194290289479364712180083596243593368443,
            ),
            (
                10151,
                56718792556032530308105358412900135578,
                8662422394028795114193864146514913943,
            ),
            (
                12345,
                33261208707367790463622745601869196757,
                268209174141567072605526753992732310247,
            ),
        ];
        for (seed, state, inc) in cases {
            let g = Pcg64::from_seed(seed);
            assert_eq!(g.state(), state, "state mismatch for seed {seed}");
            assert_eq!(g.inc(), inc, "inc mismatch for seed {seed}");
        }
    }

    /// `random_raw()` six times from `default_rng(7)`.
    #[test]
    fn next_u64_matches_numpy() {
        let mut g = Pcg64::from_seed(7);
        let want: [u64; 6] = [
            11530976094092348043,
            16550673365885938325,
            14308875409591826786,
            4154339397315733314,
            5537090637313560901,
            16114216841932056372,
        ];
        for (i, &w) in want.iter().enumerate() {
            assert_eq!(g.next_u64(), w, "draw {i}");
        }
    }
    #[test]
    fn shuffle_matches_numpy() {
        let cases: [(u128, usize, Vec<usize>); 3] = [
            (0, 5, vec![2, 4, 3, 0, 1]),
            (1, 12, vec![8, 11, 4, 7, 5, 0, 1, 9, 2, 10, 6, 3]),
            (10151, 12, vec![11, 2, 10, 5, 6, 1, 8, 0, 9, 4, 7, 3]),
        ];
        for (seed, n, want) in cases {
            let mut v: Vec<usize> = (0..n).collect();
            Pcg64Rng::seeded(seed).shuffle(&mut v);
            assert_eq!(v, want, "shuffle mismatch seed {seed} n {n}");
        }
    }

    /// The first 14 entries of `default_rng(10151).shuffle(range(1186))` -- the exact
    /// size HMMRATAC's EM down-sample uses at 10% of a ~12k-fragment track.
    #[test]
    fn shuffle_matches_numpy_at_em_sample_size() {
        let mut v: Vec<usize> = (0..1186).collect();
        Pcg64Rng::seeded(10151).shuffle(&mut v);
        assert_eq!(
            &v[..14],
            &[373, 369, 52, 977, 616, 835, 269, 88, 811, 942, 103, 184, 281, 210]
        );
    }

    /// `rng.integers(0, 2**32, dtype=uint32)` reaches `next_uint32` unfiltered
    /// (`random_bounded_uint32`'s `rng == 0xFFFFFFFF` fast path), so it is an exact
    /// probe of the 32-bit stream.
    #[test]
    fn next_uint32_matches_numpy() {
        let cases: [(u128, [u32; 4]); 3] = [
            (0, [3653403231, 2735729615, 2195314465, 1158725112]),
            (1, [2032329983, 2198257139, 3243419750, 4082210491]),
            (10151, [3144192899, 1974701719, 3953748964, 2371281220]),
        ];
        for (seed, want) in cases {
            let mut g = Pcg64Rng::seeded(seed);
            for (i, &w) in want.iter().enumerate() {
                assert_eq!(g.next_uint32_for_test(), w, "seed {seed} draw {i}");
            }
        }
    }

    /// The permutation at the size HMMRATAC actually shuffles: 11859 fragments,
    /// truncated to 1186. Checked because `random_interval`'s rejection rate depends on
    /// the bound, so a match at n=1186 is not evidence for n=11859.
    #[test]
    fn shuffle_matches_numpy_at_full_em_track_size() {
        const N: usize = 11859;
        let mut v: Vec<usize> = (0..N).collect();
        Pcg64Rng::seeded(10151).shuffle(&mut v);
        let want_head: Vec<usize> = vec![
            2300, 8679, 10978, 3270, 346, 8821, 4881, 9797, 5013, 8463, 6194, 4799, 7909, 5454,
        ];
        assert_eq!(&v[..14], want_head.as_slice());
    }
}
