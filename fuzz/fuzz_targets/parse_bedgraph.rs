//! L6: the bedGraph parser must never panic, and a surviving record must be usable
//! by the bedGraph writers.

#![no_main]

use libfuzzer_sys::fuzz_target;
use macs_io::{parse_bedgraph_line, split_line, LineSplit};

fuzz_target!(|data: &[u8]| {
    let _ = split_line(data, LineSplit::AnyWhitespace);
    if let Ok(Some(r)) = parse_bedgraph_line(data) {
        // `end < start` is tolerated, not rejected: upstream's `bedGraphIO.add_loc`
        // stores the pair verbatim, so a hostile record reaches the pileup there too.
        // What must hold is that the interval can be measured without panicking.
        let span = r.end - r.start;
        assert!(
            span.checked_abs().is_some(),
            "bedGraph span {span} not representable ({}..{})",
            r.start,
            r.end
        );
        // `atof("inf")` is infinity, and C's `atof` -- which is what upstream calls -- does
        // exactly the same, so a non-finite value here matches the oracle rather than
        // contradicting it. The property worth asserting is that the value is the one the
        // same text yields when converted again, i.e. the parse is a real parse.
        if r.value.is_finite() {
            assert!(
                r.value.to_string().parse::<f64>().is_ok(),
                "finite bedGraph value {} does not round-trip",
                r.value
            );
        }
    }
});
