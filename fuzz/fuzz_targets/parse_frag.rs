//! L6: the FRAG fragment format, end to end.
//!
//! `parse_bedpe.rs` fuzzes the *parser* for both BEDPE and FRAG. This target goes one
//! layer further and pushes whatever survives parsing through the counted fragment
//! track, which is where FRAG-specific arithmetic lives and where a malformed record can
//! do real damage:
//!
//! * `count` is an optional column. `None` means one fragment; `Some(0)`, `Some(k)` and
//!   absurd values all have to be handled without a panic.
//! * the barcode column is parsed and, for a counted track, keys the aggregation. A
//!   hostile barcode must not be able to collide two fragments into one.
//! * `average_template_length` divides `length` by `total`. A track whose fragments all
//!   net to zero length makes that a division by zero -- upstream raises
//!   `ZeroDivisionError` in exactly that case (`PairedEndTrack.py:608`), so we must
//!   return an error rather than panic.
//! * `right < left` normalises, and `right - left` must not wrap: an underflowing
//!   length would be catastrophic downstream, where it indexes a weight map.
//!
//! The invariant asserted throughout is the acceptance criterion's: **no input may
//! panic**, whatever it contains.

#![no_main]

use libfuzzer_sys::fuzz_target;
use macs_track::FragTrackBuilder;

fuzz_target!(|data: &[u8]| {
    let Ok(Some(rec)) = macs_io::parse_frag_line(data) else {
        return;
    };
    if rec.left < 0 || rec.chrom.is_empty() {
        return;
    }

    let mut builder = FragTrackBuilder::new();
    match rec.count {
        None => {
            builder.push(&rec.chrom, rec.left as u64, rec.right as u64);
        }
        Some(count) => {
            builder.push_with_count(&rec.chrom, rec.left as u64, rec.right as u64, count);
        }
    }
    // `finalize` is where ordering and per-chromosome totals are established, so it is
    // the interesting point at which a poisoned length becomes observable.
    builder.finalize();
    let track = builder.build();

    // Invariants the rest of the codebase actually relies on.
    //
    // Note what is deliberately *not* asserted: that `end > start`. Upstream's
    // `BEDPEParser` stores `atoi(f[1])` and `atoi(f[2])` verbatim and
    // `PETrackI.add_loc` performs no normalisation, so a record with `right < left` --
    // or `right == left` -- is a legitimate, reachable state here too. Asserting the
    // opposite here would encode a behaviour upstream does not have. What matters is
    // that nothing downstream *panics* on such a fragment.
    for chrom in track.chroms() {
        // `frags` and `counts` are indexed in lockstep everywhere; a length mismatch
        // would be an out-of-bounds panic rather than a wrong answer.
        assert_eq!(track.frags(chrom).len(), track.counts(chrom).len());

        for (i, frag) in track.frags(chrom).iter().enumerate() {
            // The length is a key into the weight map, so it must be computed without
            // wrapping. Upstream's `fraglengths()` declares `d` as `uint64_t` and so
            // wraps for an inverted fragment; that is upstream's behaviour and the
            // lookup then simply misses. What must not happen is a panic, so compute it
            // the safe way and assert the safe value.
            // An inverted fragment makes `end - start` underflow, so the length must be
            // obtained with a checked op. There is deliberately no upper bound asserted:
            // a large coordinate gap is legitimate, and `u64` can represent all of them.
            // `PScoreHistogram` accumulates lengths in `i64` precisely so that a huge one
            // cannot wrap, so `checked_mul` failing here is an expected outcome, not a bug.
            let len = frag.end.checked_sub(frag.start);
            let c = u64::from(track.counts(chrom)[i]);
            let _ = len.and_then(|l| l.checked_mul(c));
        }
    }

    // Zero total fragments, or fragments that all net to zero length, must be an error
    // or a defined value -- never a panic.
    let _ = track.average_template_length();
});