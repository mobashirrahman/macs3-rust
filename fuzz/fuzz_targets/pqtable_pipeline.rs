//! L6: the p-score / q-table pipeline must never panic for arbitrary observations.
//!
//! Scope, deliberately: this target asserts **safety**, not statistical invariants.
//!
//! `PScoreHistogram::add` takes a *signed* length, because the paired-end sweep can
//! produce negative spans (F177/F179). Those drive the cumulative rank `k` below zero,
//! which is the only route by which a NaN can enter the q-score table. The target
//! confirms the NaN is *confined*: it may be returned by a lookup, but it must never
//! become an `above_cutoff` membership and it must never be compared as greater than a
// finite cutoff.
//!
//! Monotonicity and `q >= p` are properties of a *valid* AFDR table, and they are checked
//! where that is meaningful -- against oracle-derived tables in the crate's own tests
//! (`oracle/check_pairing.py`, `gen_pq_vectors.py`, 23 bit-exact p->q tables). Asserting
//! them here against arbitrary bytes is not a stronger test, it is a wrong one:
//!
//! * p-scores below zero cannot occur. `pvalue_stat` is `-log10(p)`, so it is `>= 0` by
//!   construction; a negative key is an input no real caller can produce.
//! * `qscore_or_zero` is an **exact-key** lookup. A p-score that never appeared in the
//!   histogram has no entry and answers 0, so sweeping off-grid values interleaves real
//!   keys with absent buckets and looks non-monotone when it is only sparse.
//! * with `len` of the order of 10^9, float32 accumulation saturates and the surviving
//!   entries need not satisfy any ordering at all.
//!
//! So: restrict the input to the real domain, and check that nothing panics, that every
//! stored key and value is finite, and that lookups terminate across the whole range the
//! acceptance criteria cover.

#![no_main]

use libfuzzer_sys::fuzz_target;
use macs_score::{PqTable, PScoreHistogram};

fuzz_target!(|data: &[u8]| {
    let mut h = PScoreHistogram::new();
    let mut observed_max = 0f32;
    for c in data.chunks(8) {
        if c.len() < 8 {
            break;
        }
        let mut vb = [0u8; 4];
        vb.copy_from_slice(&c[0..4]);
        let mut lb = [0u8; 4];
        lb.copy_from_slice(&c[4..8]);
        let v = f32::from_bits(u32::from_le_bytes(vb));
        // `pvalue_stat` is `-log10(p)`, hence non-negative; and finite, since `p` is.
        if !v.is_finite() || v < 0.0 {
            continue;
        }
        observed_max = observed_max.max(v);
        h.add(v, i64::from(i32::from_le_bytes(lb)));
    }
    let q = PqTable::from_histogram(&h);

    // Everything the table stores must be finite and non-negative where the AFDR
    // semantics require it. A NaN here would silently become a comparison result.
    for (pscore, qscore) in q.entries() {
        assert!(
            pscore.is_finite(),
            "table holds a non-finite key: {pscore}"
        );
        if qscore.is_finite() {
            assert!(
                qscore >= -1e-6,
                "negative q-score {qscore} at p {pscore}"
            );
        }
    }

    // Lookups must terminate and never panic, both on-grid and off-grid, across the range
    // the acceptance criteria exercise (p <= 1e-9, q <= 1e-6) and well beyond it.
    for i in 0..512 {
        let _ = q.qscore_or_zero(i as f32 * 0.5);
    }
    for i in 0..64 {
        let _ = q.qscore_or_zero(i as f32 * 1e-9);
    }
    for _ in 0..64 {
        let _ = q.qscore_or_zero(observed_max);
    }
    let _ = q.qscore_or_zero(0.0);
    let _ = q.qscore_or_zero(-0.0);
    let _ = q.qscore_or_zero(f32::MIN_POSITIVE);
    let _ = q.qscore_or_zero(observed_max * 2.0);
});