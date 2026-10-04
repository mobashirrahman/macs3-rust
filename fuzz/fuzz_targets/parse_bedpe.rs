//! L6: the BEDPE and FRAG fragment parsers must never panic.
//!
//! Fragment coordinates are the one place where malformed input can yield an
//! *inverted* interval (`right < left`) rather than an error, because upstream
//! accepts it. However we normalise it, it must not be a panic, and the resulting
//! interval length must not wrap -- an underflowing `len` would blow up the peak
//! union walk.

#![no_main]

use libfuzzer_sys::fuzz_target;
use macs_core::Interval;
use macs_io::{parse_bedpe_line, parse_frag_line, split_line, LineSplit};

fuzz_target!(|data: &[u8]| {
    let _ = split_line(data, LineSplit::TabAfterRstrip);

    for rec in [parse_bedpe_line(data), parse_frag_line(data)] {
        if let Ok(Some(f)) = rec {
            // `Interval::new` normalises `end >= start`; that is the guarantee the
            // rest of the port relies on, so assert it directly.
            let iv = Interval::new(f.left as u64, f.right as u64);
            assert!(iv.end() >= iv.start(), "Interval::new produced an inverted range");
            if let Some(c) = f.count {
                assert!(c <= 65535, "multiplicity {} exceeds the upstream cap", c);
            }
        }
    }
});
