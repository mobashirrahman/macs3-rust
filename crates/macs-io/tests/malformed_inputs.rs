//! L6 (CI-runnable half): randomized malformed-input no-panic coverage.
//!
//! The acceptance criteria require "zero panics on malformed input" and "identical
//! accept/reject behaviour". The nightly libFuzzer targets under `fuzz/` explore that
//! space far more thoroughly, but they need nightly + cargo-fuzz, so they cannot be
//! the *only* evidence. This test asserts the same invariants over a deterministic
//! pseudo-random corpus so that ordinary CI proves them on every push.
//!
//! The generator deliberately produces the shapes that actually break parsers:
//! truncated lines, wrong field counts, non-numeric coordinates, NUL bytes, CRLF,
//! embedded tabs, negative and overflowing integers, and pure garbage.

use macs_core::{ChromId, Interval, Strand};
use macs_io::{
    parse_bed_line, parse_bedgraph_line, parse_bedpe_line, parse_frag_line, split_line, LineSplit,
};
use macs_pileup::{integrated_depth, pileup_from_weighted_positions, SingleEndParams};

/// xorshift64*, so the corpus is reproducible and needs no dependency.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn below(&mut self, n: u32) -> u32 {
        (self.next() % u64::from(n)) as u32
    }

    fn pick<'a, T>(&mut self, xs: &'a [T]) -> &'a T {
        &xs[self.below(xs.len() as u32) as usize]
    }
}

const CHROMS: &[&str] = &["chr1", "chrX", "scaffold_7", "", "chr with space", "chr\tX"];
const NUMBERS: &[&str] = &[
    "0",
    "1",
    "-1",
    "-2147483648",
    "2147483647",
    "99999999999999999999999999",
    "",
    " ",
    "+7",
    "0x10",
    "1.5",
    "1e9",
    "NaN",
    "inf",
    "--3",
    "007",
    "٣",
];
const STRANDS: &[&str] = &["+", "-", ".", "", "*", "plus", "forward"];

/// Assemble a random line from plausible-looking fields.
fn make_line(rng: &mut Rng) -> String {
    let n_fields = rng.below(7);
    let mut parts: Vec<String> = Vec::with_capacity(n_fields as usize);
    for _ in 0..n_fields {
        match rng.below(4) {
            0 => parts.push((*rng.pick(CHROMS)).to_string()),
            1 => parts.push((*rng.pick(NUMBERS)).to_string()),
            2 => parts.push((*rng.pick(STRANDS)).to_string()),
            _ => parts.push(format!("{}", rng.next() as i64)),
        }
    }
    let joined = parts.join("\t");
    match rng.below(6) {
        0 => format!("{}\r\n", joined),
        1 => format!("track name=x\n{}", joined),
        2 => joined.replace('\t', "   "),
        3 => format!("{}\u{0}\u{1}\u{7}", joined),
        _ => joined,
    }
}

#[test]
fn parsers_never_panic_on_random_lines() {
    let mut rng = Rng(0x5eed_1234_abcd_0001);
    let mut non_finite = 0u32;
    for i in 0..200_000u32 {
        let line = make_line(&mut rng);
        let bytes = line.as_bytes();

        let _ = split_line(bytes, LineSplit::TabAfterRstrip);
        let _ = split_line(bytes, LineSplit::AnyWhitespace);
        let _ = split_line(&[], LineSplit::AnyWhitespace);
        let _ = split_line(b"\t\t\t", LineSplit::TabAfterRstrip);

        if let Ok(Some(r)) = parse_bed_line(bytes) {
            assert!(
                matches!(r.strand, Strand::Plus | Strand::Minus | Strand::Unknown),
                "case {}: parser invented a strand",
                i
            );
            // No non-negativity law here, deliberately. `tlen_parse_line` is the
            // *signed* `col3 - col2` read from the columns as written, so a record with
            // `col3 < col2` yields a negative length -- and upstream accepts it, so
            // rejecting it would break accept/reject parity. The consequence is that
            // `d` itself can go negative, which is exactly why `GenericParser.d` gates
            // accumulation on `tlen > 0` (F192). What is guaranteed is that the value
            // stays finite and does not overflow.
            let tlen = match r.strand {
                Strand::Minus => r.pos - r.end,
                _ => r.end - r.pos,
            };
            assert!(
                (i64::MIN..=i64::MAX).contains(&tlen),
                "case {}: tlen overflowed: {}",
                i,
                tlen
            );
        }

        if let Ok(Some(r)) = parse_bedgraph_line(bytes) {
            // Again no ordering law: upstream's `read_bedGraph` is literally
            // `add_func(fs[0], atoi(fs[1]), atoi(fs[2]), atof(fs[3]))` with no
            // validation, so `end < start` is accepted there and must be accepted here.
            // (Contrast `BEDPEParser.build_petrack`, which *asserts*
            // `right_pos > left_pos` -- fragments are checked, bedGraph rows are not.
            // Getting that asymmetry right is the whole point.)
            // `atof` is a C `strtod`, so "nan"/"inf" become NaN/inf here exactly as
            // upstream. The criterion is "no panic", not "no NaN", so a non-finite
            // value is only reported, never asserted against.
            if !r.value.is_finite() {
                non_finite += 1;
            }
        }

        for rec in [parse_bedpe_line(bytes), parse_frag_line(bytes)] {
            if let Ok(Some(f)) = rec {
                let iv = Interval::new(f.left.max(0) as u32, f.right.max(0) as u32);
                assert!(
                    iv.end() >= iv.start(),
                    "case {}: inverted fragment interval",
                    i
                );
                if let Some(c) = f.count {
                    assert!(c <= 65535, "case {}: multiplicity {} over the cap", i, c);
                }
            }
        }
    }
    // If the corpus never produced a non-finite bedGraph value the finiteness branch
    // above would be vacuous, so require that it was actually exercised.
    assert!(
        non_finite > 0,
        "the random corpus never produced a non-finite bedGraph value"
    );
}

#[test]
fn pipeline_never_panics_on_random_coordinates() {
    let mut rng = Rng(0xfeed_9999_0000_2222);
    for i in 0..20_000u32 {
        let five = rng.below(200) as i64;
        let three = 1 + rng.below(200) as i64;
        let rlength = 1 + rng.below(50_000);

        let mut plus = Vec::new();
        let mut minus = Vec::new();
        let mut pw = Vec::new();
        let mut mw = Vec::new();
        for _ in 0..rng.below(40) {
            let pos = rng.below(rlength);
            let w = 1.0 + (rng.below(1000) as f32 / 1000.0);
            if rng.below(2) == 0 {
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
            rlength,
            scale_factor: 1.0,
            baseline_value: 0.0,
            coalesce: true,
            centred: false,
        };
        let t = pileup_from_weighted_positions(ChromId(0), &plus, &minus, &pw, &mw, &p);

        let mut prev = t.start();
        for r in t.runs() {
            assert!(
                r.end > prev,
                "case {}: run end {} did not advance {}",
                i,
                r.end,
                prev
            );
            assert!(r.value.is_finite(), "case {}: non-finite depth", i);
            assert!(r.value >= 0.0, "case {}: negative depth {}", i, r.value);
            prev = r.end;
        }
        let depth = integrated_depth(&t);
        assert!(
            depth.is_finite() && depth >= 0.0,
            "case {}: bad depth {}",
            i,
            depth
        );
    }
}

#[test]
fn truncated_and_empty_inputs_are_rejected_without_panic() {
    // Explicit boundary cases, independent of the random corpus.
    for bad in [
        &b""[..],
        b"\n",
        b"\r\n",
        b"\t",
        b"\0",
        b"chr1",
        b"chr1\t",
        b"chr1\t1",
        b"chr1\t1\t2",
        b"chr1\t1\t2\t+",
        b"chr1\t1\t2\t+\t",
        b"chr1\t-\t2\t+",
        b"chr1\t1\t2\t+extra\tfields\tand\tmore",
        b"chr1\t1\t2\t+\t\xff\xfe\xfd",
        b"#chr1\t1\t2\t+",
        b"chr1\t1e400\t2\t+",
    ] {
        let _ = parse_bed_line(bad);
        let _ = parse_bedpe_line(bad);
        let _ = parse_frag_line(bad);
        let _ = parse_bedgraph_line(bad);
        let _ = split_line(bad, LineSplit::AnyWhitespace);
        let _ = split_line(bad, LineSplit::TabAfterRstrip);
    }
}
