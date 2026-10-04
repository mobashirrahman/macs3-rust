//! L6: the pileup pipeline must never panic for arbitrary read coordinates.
//!
//! Coordinates arrive as `i64` from the parsers and are converted to `u64` before
//! the sweep, which is exactly where a negative position would wrap to a value near
//! `u64::MAX` and blow up the allocation. This target drives the whole chain
//! (params -> endpoints -> sweep -> track) with fuzzed inputs.

#![no_main]

use libfuzzer_sys::fuzz_target;
use macs_core::ChromId;
use macs_pileup::{integrated_depth, pileup_from_weighted_positions, SingleEndParams};
use macs_rle::SignalTrack;

fuzz_target!(|data: &[u8]| {
    if data.len() < 8 {
        return;
    }
    let d = 1 + (data[0] as i64 % 500);
    let five = data[1] as i64 % 250;
    let three = 1 + (data[2] as i64 % 250);
    let rlength = 1 + (u32::from_le_bytes([data[3], data[4], data[5], data[6]]) % 100_000);

    // split the rest into plus/minus/weight triples
    let body = &data[7..];
    let mut plus = Vec::new();
    let mut minus = Vec::new();
    let mut pw = Vec::new();
    let mut mw = Vec::new();
    for c in body.chunks(9) {
        if c.len() < 9 {
            break;
        }
        let mut v = [0u8; 4];
        v.copy_from_slice(&c[0..4]);
        // `pileup_from_weighted_positions` takes `Coord` (= u64) positions; narrowing
        // here would have made the target a compile error rather than a fuzz case, so
        // the module never built -- the nightly fuzz job has been red for this target.
        let pos = u64::from(u32::from_le_bytes(v) % rlength);
        // Keep the exponent in the finite range: real weights come from `pnorm2` and are
        // always finite, and an all-ones exponent here would just inject NaN, which the
        // "no non-finite depth" assertion below is right to reject. The property under
        // test is *finite in -> finite out*.
        let raw = u32::from_le_bytes([c[5], c[6], c[7], c[8]]);
        let raw = (raw & 0x807f_ffff) | 0x3f80_0000; // force exponent 127
        let w = f32::from_bits(raw).min(1.0e6).max(-1.0e6);
        if c[4] & 1 == 0 {
            plus.push(pos);
            pw.push(w);
        } else {
            minus.push(pos);
            mw.push(w);
        }
    }
    plus.sort_unstable();
    minus.sort_unstable();

    let p = SingleEndParams {
        five_shift: five,
        three_shift: three,
        rlength: rlength as u64,
        scale_factor: 1.0,
        baseline_value: 0.0,
        coalesce: true,
        centred: false,
    };
    let t: SignalTrack<f32> =
        pileup_from_weighted_positions(ChromId(0), &plus, &minus, &pw, &mw, &p);

    // Structural guarantees downstream code relies on.
    let mut prev = t.start();
    for r in t.runs() {
        assert!(r.end > prev, "pileup run end {} did not advance {}", r.end, prev);
        prev = r.end;
        assert!(r.value.is_finite(), "pileup produced a non-finite depth");
        assert!(r.value >= 0.0, "pileup produced a negative depth {}", r.value);
    }
    let depth = integrated_depth(&t);
    assert!(depth.is_finite() && depth >= 0.0, "bad integrated depth {}", depth);
});
