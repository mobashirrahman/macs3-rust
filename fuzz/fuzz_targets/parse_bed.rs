//! L6: the BED line parser must never panic on arbitrary bytes.
//!
//! The acceptance criterion is "zero panics on malformed input". Upstream raises on
//! bad input; we return `MacsError`, which the CLI turns into exit 1. This target
//! asserts only that no *panic* escapes -- any `Err` is an acceptable, expected
//! outcome -- plus the consistency laws that must hold for whatever does parse.

#![no_main]

use libfuzzer_sys::fuzz_target;
use macs_core::Strand;
use macs_io::{parse_bed_line, split_line, LineSplit};

fuzz_target!(|data: &[u8]| {
    // splitting must not panic on any byte string, including empty and NUL-heavy
    let _ = split_line(data, LineSplit::TabAfterRstrip);
    let _ = split_line(data, LineSplit::AnyWhitespace);

    if let Ok(Some(r)) = parse_bed_line(data) {
        assert!(
            matches!(r.strand, Strand::Plus | Strand::Minus | Strand::Unknown),
            "parser invented a strand"
        );
        // A minus-strand tag carries `end` as the 5' coordinate, so its length is
        // `pos - end`.
        //
        // That difference is *allowed* to come out negative: `BEDParser` stores
        // `atoi(f[1])`/`atoi(f[2])` verbatim and performs no normalisation, so a record
        // with `end < start` is a reachable, legitimate state that upstream also carries
        // into `BEDParser::d`. Asserting non-negativity here would be asserting a
        // behaviour upstream does not have, and it fires within seconds on arbitrary
        // bytes -- so it tests the target, not the parser.
        //
        // What the acceptance criteria do require is that obtaining the length cannot
        // panic, which is why it is computed at all: on the i64 coordinates the parser
        // returns, `+`/`-` cannot overflow, so a wrap would mean a sign mistake.
        let len = match r.strand {
            Strand::Minus => r.pos - r.end,
            _ => r.end - r.pos,
        };
        assert!(
            len.checked_abs().is_some(),
            "template length {len} is not representable"
        );
    }
});
