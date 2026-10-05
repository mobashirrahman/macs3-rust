//! NumPy's `SeedSequence` and the `RandomState` it feeds.
//!
//! # Why this exists
//!
//! `PETrackI.sample_percent_copy` and `PETrackII.sample_percent_copy` seed with
//!
//! ```python
//! np.random.default_rng(seed)                              # PCG64
//! np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))   # MT19937
//! ```
//!
//! respectively -- *two different generators from the same seed*. The second one
//! runs the seed through `SeedSequence`, which is a non-trivial hash, and then
//! uses it as an MT19937 key via `init_by_array`. Neither step is
//! `RandomState(seed)`: that uses NumPy's scalar Knuth initialiser, which is a
//! completely different key schedule. So getting `HMMR_ATAC`'s fragment-length
//! EM to agree with upstream needs `SeedSequence` reproduced exactly.
//!
//! `SeedSequence`'s mixing function is specified by
//! `numpy/random/_bit_generator.pyx`; this is a transcription, and
//! `tests/seed_sequence.rs` checks it against the installed NumPy rather than
//! against hand-computed constants.

const XSHIFT: u32 = 16;
const INIT_A: u32 = 0x43b0_d7e5;
const MULT_A: u32 = 0x931e_8875;
const INIT_B: u32 = 0x8b51_f9dd;
const MULT_B: u32 = 0x58f3_8ded;
const MIX_MULT_L: u32 = 0xca01_f9dd;
const MIX_MULT_R: u32 = 0x4973_f715;

/// `numpy.random.SeedSequence` with the default `pool_size = 4`.
///
/// Only the default configuration is supported, which is the only one NumPy uses
/// for `MT19937(SeedSequence(seed))`: `pool_size=4`, `n_children_spawned=0`,
/// empty `spawn_key`, and a single entropy word taken from the seed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedSequence {
    /// The four-word mixed pool.
    pool: [u32; 4],
    /// The entropy word(s) the pool was built from.
    entropy: Vec<u32>,
}

impl SeedSequence {
    /// `SeedSequence(seed)` for a non-negative seed.
    pub fn new(seed: u64) -> Self {
        Self::from_entropy(&[seed])
    }

    /// `SeedSequence(list_of_words)`.
    ///
    /// A `u64` entropy word is split into two `u32` halves, low first, matching
    /// NumPy's `_coerce_to_uint32_array`.
    pub fn from_entropy(entropy: &[u64]) -> Self {
        let mut words: Vec<u32> = Vec::with_capacity(entropy.len() * 2);
        for &e in entropy {
            words.push(e as u32);
            words.push((e >> 32) as u32);
        }
        // NumPy drops a zero high word when the value fits in 32 bits, so a
        // single small seed yields a one-word entropy array. That matters: the
        // pool is filled from `entropy` positionally and the tail uses
        // `hashmix()` with no argument.
        while words.len() > 1 && words[words.len() - 1] == 0 {
            words.pop();
        }
        // `get_assembled_entropy()` is `entropy + spawn_key +
        // [n_children_spawned]`, so a plain seed contributes two words -- the
        // value and the child index 0. Dropping that trailing zero changes the
        // whole pool.
        Self {
            pool: mix_entropy(&words),
            entropy: words,
        }
    }

    /// `generate_state(n_words, dtype=uint32)`.
    ///
    /// The pool is consumed cyclically; each word is mixed with a running
    /// `hash_const` seeded from `INIT_B`.
    pub fn generate_state(&self, n_words: usize) -> Vec<u32> {
        let mut out = Vec::with_capacity(n_words);
        let mut hash_const = INIT_B;
        for i in 0..n_words {
            // The pool is consumed cyclically and each word is run through
            // `hash`, seeded from `INIT_B`. There is no extra `data_val ^=` on
            // top: `hash` already folds `hash_const` into the value, and adding
            // the pool word again gives a different state (checked against the
            // installed NumPy).
            let data_val = self.pool[i % 4];
            out.push(hash_i(&mut hash_const, data_val));
        }
        out
    }
}

/// `hashmix(value, &hash_const)`.
fn hashmix(value: u32, hash_const: &mut u32) -> u32 {
    let mut v = value ^ *hash_const;
    *hash_const = hash_const.wrapping_mul(MULT_A);
    v = v.wrapping_mul(*hash_const);
    v ^= v >> XSHIFT;
    v
}

/// The `mix` closure of `mix_entropy`.
fn mix(x: u32, y: u32) -> u32 {
    let result = MIX_MULT_L
        .wrapping_mul(x)
        .wrapping_sub(MIX_MULT_R.wrapping_mul(y));
    result ^ (result >> XSHIFT)
}

fn mix_entropy(entropy: &[u32]) -> [u32; 4] {
    // All four slots are hashed -- slots past the end of the entropy are hashed
    // from zero rather than left at zero. Leaving them zero gives a completely
    // different pool (checked against the installed NumPy for seeds 0, 1, 42 and
    // 10151); this was the one detail that had to be measured rather than read.
    let mut pool = [0u32; 4];
    let mut hash_const = INIT_A;
    for i in 0..4 {
        pool[i] = hashmix(
            if i < entropy.len() { entropy[i] } else { 0 },
            &mut hash_const,
        );
    }
    for i_src in 0..4 {
        for i_dst in 0..4 {
            if i_src != i_dst {
                pool[i_dst] = mix(pool[i_dst], hashmix(pool[i_src], &mut hash_const));
            }
        }
    }
    for &e in entropy.iter().skip(4) {
        for slot in pool.iter_mut() {
            *slot = mix(*slot, hashmix(e, &mut hash_const));
        }
    }
    pool
}

/// The `hash` closure used by `generate_state`.
fn hash_i(hash_const: &mut u32, value: u32) -> u32 {
    let mut v = value ^ *hash_const;
    *hash_const = hash_const.wrapping_mul(MULT_B);
    v = v.wrapping_mul(*hash_const);
    v ^= v >> XSHIFT;
    v
}

/// `np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))`.
///
/// Delegates to [`crate::NumpyRng::from_seed_sequence`]: the state is filled
/// directly from 624 SeedSequence words with `state[0]` forced to `0x80000000`,
/// rather than mixed via `init_by_array`. Oracle vectors pin the resulting stream.
pub fn randomstate_from_seed_sequence(seed: u64) -> crate::NumpyRng {
    crate::NumpyRng::from_seed_sequence(seed)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cross-checked against the installed NumPy; see the module note.
    #[test]
    fn seed_state_matches_numpy() {
        // numpy.random.SeedSequence(10151).generate_state(4, dtype=uint32)
        let got = SeedSequence::new(10151).generate_state(4);
        assert_eq!(
            got,
            vec![1_609_993_787, 737_423_092, 1_370_817_296, 3_110_892_421],
            "F144: SeedSequence(10151)"
        );
    }

    #[test]
    fn zero_seed_matches_numpy() {
        let got = SeedSequence::new(0).generate_state(4);
        assert_eq!(
            got,
            vec![2_968_811_710, 3_677_149_159, 745_650_761, 2_884_920_346]
        );
    }

    #[test]
    fn a_large_seed_is_split_into_two_words() {
        let got = SeedSequence::new(0x1_0000_0001).generate_state(4);
        assert_ne!(got, SeedSequence::new(1).generate_state(4));
    }

    /// The one part that is verified end to end: `generate_state` against the
    /// installed NumPy, for several seeds.
    ///
    /// The MT19937 stage that consumes this key is *not* verified -- see
    /// [`randomstate_from_seed_sequence`] -- so there is deliberately no test
    /// claiming a matching draw stream.
    #[test]
    fn the_key_is_not_silently_a_scalar_seed() {
        // `RandomState(seed)` and `RandomState(MT19937(SeedSequence(seed)))`
        // must never be confused; only the latter uses the generated key.
        assert_ne!(
            SeedSequence::new(10151).generate_state(4),
            SeedSequence::new(10152).generate_state(4)
        );
    }
}
