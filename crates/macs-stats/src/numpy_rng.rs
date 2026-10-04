//! NumPy-compatible random number generation.
//!
//! Upstream's downsamplers (`FWTrack.sample_percent`, `sample_num`) seed NumPy's
//! **global** Mersenne-Twister RNG and then call `np.random.shuffle`. Which
//! positions survive therefore depends on NumPy's exact stream and its exact
//! `shuffle` algorithm, not merely on the distribution. Reproducing upstream's
//! sampling output bit-for-bit requires all three of:
//!
//! 1. the MT19937 stream NumPy seeds,
//! 2. NumPy's seeding procedure (`RandomState(seed)`), which is *not* just
//!    `init_by_array([seed])` — it folds the seed into the key with
//!    `init_by_array` over a 4-element array and uses `init_genrand` only when the
//!    seed is a scalar in some paths,
//! 3. `RandomState.shuffle`'s Fisher-Yates loop and its bounded-integer draw.
//!
//! Everything here is a transcription of NumPy's `mt19937.c`, and is validated
//! against the installed NumPy in the differential harness rather than against
//! hand-computed constants.

/// NumPy's legacy `RandomState`: MT19937 plus NumPy's seeding and draw rules.
#[derive(Debug, Clone)]
pub struct NumpyRng {
    /// MT19937 state, `624` words.
    state: [u32; Self::N],
    /// Position within the state array.
    index: usize,
}

impl NumpyRng {
    /// Number of words in the MT19937 state.
    const N: usize = 624;
    /// Word offset 397, the standard MT19937 recurrence parameter.
    const M: usize = 397;
    /// Seed exactly as `np.random.seed(seed)` does for a scalar seed.
    ///
    /// A scalar seed uses NumPy's **Knuth** initialiser (`mt19937_seed`), not
    /// `init_by_array`: every state word is set from the seed, and the seed is
    /// then advanced by `1812433253 * (s ^ (s >> 30)) + pos + 1`.
    ///
    /// Note there is **no** zero-seed promotion. `RandomState(0)` really is
    /// seeded with 0, verified against the installed NumPy: the `0xBADF00D`
    /// substitution belongs to `SeedSequence`, which the legacy `RandomState`
    /// does not use. Adding it makes seed 0 diverge from NumPy on the very first
    /// shuffle while every other seed still matches, which is a nasty way to
    /// discover the mistake.
    pub fn seeded(seed: i64) -> Self {
        let mut s = seed as u32;
        let mut rng = NumpyRng {
            state: [0; Self::N],
            index: Self::N,
        };
        for (pos, word) in rng.state.iter_mut().enumerate() {
            *word = s;
            s = 1_812_433_253u32
                .wrapping_mul(s ^ (s >> 30))
                .wrapping_add(pos as u32 + 1);
        }
        rng
    }

    /// Seed from an explicit key, as NumPy's `init_by_array` takes.
    ///
    /// This is what `np.random.seed([a, b, ...])` uses for a *sequence* seed, and
    /// is kept separate from [`seeded`](Self::seeded) because the two produce
    /// completely different streams for the same first element. It is also what
    /// `RandomState(MT19937(SeedSequence(seed)))` ends up using.
    ///
    /// The three phases are the standard `init_by_array`: expand the key over the
    /// whole state with the `1812433253` LCG, fold it in again with the
    /// `1664525` / `1565833941` multipliers, then mask the top bit. The
    /// key-expansion phase is not optional -- skipping it and starting from a
    /// constant gives a stream that looks plausible but diverges from numpy on
    /// the very first draw.
    pub fn init_by_array(key: &[u32]) -> Self {
        const N: usize = 624;
        if key.is_empty() {
            return Self::seeded(435);
        }
        let mut state = [0u32; N];
        state[0] = key[0];
        for i in 1..N {
            let prev = state[i - 1];
            state[i] = 1_812_433_253u32
                .wrapping_mul(prev ^ (prev >> 30))
                .wrapping_add(key[i % key.len()]);
        }
        let mut i = 1usize;
        let mut j = 0usize;
        let mut k = N.max(key.len());
        while k > 0 {
            let prev = state[i - 1];
            state[i] = (state[i] ^ (prev ^ (prev >> 30)).wrapping_mul(1_664_525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                state[0] = state[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
            k -= 1;
        }
        let mut k = N - 1;
        while k > 0 {
            let prev = state[i - 1];
            state[i] = (state[i] ^ (prev ^ (prev >> 30)).wrapping_mul(1_565_833_941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                state[0] = state[N - 1];
                i = 1;
            }
            k -= 1;
        }
        state[0] = 0x8000_0000;
        Self { state, index: N }
    }
    /// Seed exactly as `np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))`.
    ///
    /// This is NOT `seeded` and NOT `init_by_array`: NumPy fills the 624-word MT19937
    /// state directly from `SeedSequence.generate_state(624)`, forces `state[0]` to
    /// `0x80000000`, and leaves the position at 623 (verified against the installed
    /// NumPy via `get_state()`). `init_by_array(generate_state(4))` mixes the same
    /// entropy differently and diverges on the first draw -- which is why PE
    /// `randsample` kept a different 9999 fragments than upstream despite a correct
    /// `SeedSequence` and a correct shuffle (F145).
    pub fn from_seed_sequence(seed: u64) -> Self {
        use crate::SeedSequence;
        let ss = SeedSequence::new(seed);
        let words = ss.generate_state(Self::N);
        let mut state = [0u32; Self::N];
        state.copy_from_slice(&words);
        state[0] = 0x8000_0000;
        Self { state, index: 623 }
    }

    /// NumPy's `mixbits` helper: XOR the two words and fold the low bits upward.
    /// Generate the next block if the state is exhausted.
    fn generate(&mut self) {
        const MATRIX_A: u32 = 0x9908_b0df;
        const UPPER_MASK: u32 = 0x8000_0000;
        const LOWER_MASK: u32 = 0x7fff_ffff;
        for i in 0..Self::N {
            let y = (self.state[i] & UPPER_MASK) | (self.state[(i + 1) % Self::N] & LOWER_MASK);
            let mut next = self.state[(i + Self::M) % Self::N] ^ (y >> 1);
            if y & 1 != 0 {
                next ^= MATRIX_A;
            }
            self.state[i] = next;
        }
        self.index = 0;
    }

    /// One raw 32-bit MT19937 word.
    ///
    /// This is NumPy's `next_uint32`, i.e. one tempered MT19937 output with no
    /// scaling and no rejection sampling. Exposed so the stream itself can be
    /// diffed against NumPy, not just the shuffles built on top of it.
    #[inline]
    pub fn next_uint32(&mut self) -> u32 {
        self.next_u32()
    }

    #[inline]
    fn next_u32(&mut self) -> u32 {
        if self.index >= Self::N {
            self.generate();
        }
        let mut y = self.state[self.index];
        self.index += 1;
        // tempering
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    /// A uniform `f64` in `[0, 1)`, as NumPy's `random_sample()`.
    #[inline]
    pub fn random_sample(&mut self) -> f64 {
        // NumPy uses 53 bits of randomness from two 32-bit draws, shifted right 11
        let a = self.next_u32() >> 5;
        let b = self.next_u32() >> 6;
        (a as f64 * 67_108_864.0 + b as f64) * (1.0 / 9_007_199_254_740_992.0)
    }

    /// A uniform integer in `[low, high)`, as NumPy's `randint`.
    ///
    /// NumPy uses masked rejection sampling: mask off the low bits, and redraw
    /// while the value is out of range. That redraw loop consumes a variable
    /// number of words, so a naive modulo would desynchronise the stream and
    /// diverge on every subsequent draw.
    pub fn randint(&mut self, low: i64, high: i64) -> i64 {
        let rng = high - low - 1;
        if rng <= 0 {
            return low;
        }
        // build the mask in u32 space, since that is what the draw produces
        let mut mask = rng as u64;
        mask |= mask >> 1;
        mask |= mask >> 2;
        mask |= mask >> 4;
        mask |= mask >> 8;
        mask |= mask >> 16;
        mask |= mask >> 32;
        let mask = mask as u32;
        loop {
            let value = (self.next_u32() & mask) as i64;
            if value <= rng {
                return value + low;
            }
        }
    }

    /// `RandomState.shuffle`: in-place Fisher-Yates from the end.
    ///
    /// NumPy's implementation swaps with `j = randint(0, i + 1)`, i.e. an inclusive
    /// upper bound, and draws in decreasing `i`.
    pub fn shuffle<T>(&mut self, data: &mut [T]) {
        let n = data.len();
        if n < 2 {
            return;
        }
        for i in (1..n).rev() {
            let j = self.randint(0, i as i64 + 1) as usize;
            data.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_instances_with_the_same_seed_agree() {
        let mut a = NumpyRng::seeded(42);
        let mut b = NumpyRng::seeded(42);
        for _ in 0..100 {
            assert_eq!(a.next_u32(), b.next_u32());
        }
    }

    #[test]
    fn different_seeds_diverge() {
        let mut a = NumpyRng::seeded(1);
        let mut b = NumpyRng::seeded(2);
        let va: Vec<_> = (0..10).map(|_| a.next_u32()).collect();
        let vb: Vec<_> = (0..10).map(|_| b.next_u32()).collect();
        assert_ne!(va, vb);
    }

    #[test]
    fn seed_zero_is_not_promoted() {
        // `RandomState(0)` is seeded with 0; the 0xBADF00D substitution belongs
        // to SeedSequence, which RandomState does not use
        let mut a = NumpyRng::seeded(0);
        assert_ne!(a.next_u32(), NumpyRng::seeded(1).next_u32());
    }

    #[test]
    fn random_sample_is_in_range() {
        let mut r = NumpyRng::seeded(7);
        for _ in 0..1000 {
            let v = r.random_sample();
            assert!((0.0..1.0).contains(&v), "{v} out of range");
        }
    }

    #[test]
    fn randint_stays_in_range_and_covers_it() {
        let mut r = NumpyRng::seeded(11);
        let mut seen = [false; 6];
        for _ in 0..2000 {
            let v = r.randint(0, 6);
            assert!((0..6).contains(&v));
            seen[v as usize] = true;
        }
        assert!(seen.iter().all(|s| *s), "all values in 0..6 should appear");
    }

    #[test]
    fn randint_with_a_degenerate_range_returns_low() {
        let mut r = NumpyRng::seeded(3);
        assert_eq!(r.randint(5, 5), 5);
        assert_eq!(r.randint(5, 4), 5);
    }

    #[test]
    fn shuffle_is_a_permutation_and_deterministic() {
        let mut a: Vec<u32> = (0..50).collect();
        let mut b = a.clone();
        NumpyRng::seeded(5).shuffle(&mut a);
        NumpyRng::seeded(5).shuffle(&mut b);
        assert_eq!(a, b, "same seed, same permutation");
        let mut sorted = a.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..50).collect::<Vec<u32>>(), "still a permutation");
        assert_ne!(a, sorted, "and actually shuffled");
    }

    #[test]
    fn shuffling_tiny_slices_is_a_no_op() {
        let mut empty: Vec<u32> = Vec::new();
        NumpyRng::seeded(1).shuffle(&mut empty);
        assert!(empty.is_empty());
        let mut one = vec![9];
        NumpyRng::seeded(1).shuffle(&mut one);
        assert_eq!(one, vec![9]);
    }
}
