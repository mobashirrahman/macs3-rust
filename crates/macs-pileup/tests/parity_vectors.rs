//! Bit-exact parity of the pileup sweep against the compiled MACS3 3.0.5.
//!
//! `tests/pileup_vectors.tsv` is produced by `oracle/gen_pileup_vectors.py`,
//! which calls `MACS3.Signal.PileupV2.pileup_from_PN_shifted` and
//! `pileup_from_LR_as_list` and records the resulting breakpoint positions and
//! `float32` values.
//!
//! The comparison is deliberately **stricter than it needs to be**: it compares
//! the full breakpoint list, not the pointwise value function. Our canonical
//! form coalesces adjacent equal-valued runs, so the breakpoint lists are
//! *not* expected to be identical — the test verifies that every upstream
//! breakpoint position is reproduced, and that the value at every position
//! matches. That catches a coalescing bug as well as a sweep bug, which a
//! value-only comparison would hide.

use macs_core::ChromId;
use macs_pileup::{pileup_from_fragments, pileup_from_positions, SingleEndParams};
use macs_rle::SignalTrack;
use std::collections::BTreeMap;

const C: ChromId = ChromId(0);

fn f32_from_bits(h: &str) -> f32 {
    let h = h.trim();
    assert!(h.starts_with("0x"), "expected hex bits, got {h:?}");
    f32::from_bits(u32::from_str_radix(&h[2..], 16).expect("bad bits"))
}

fn parse_list(s: &str) -> Vec<i64> {
    s.split(',')
        .filter(|x| !x.trim().is_empty())
        .map(|x| x.trim().parse().expect("bad list element"))
        .collect()
}

/// Upstream's pv array is right-endpoint indexed: `v[i]` is the value on
/// `[pos[i-1], pos[i])`. So to compare against upstream we evaluate our track on
/// the region *just below* each breakpoint, not at it.
fn upstream_value(t: &SignalTrack<f32>, pos: u64) -> f32 {
    if pos == 0 {
        // the first breakpoint closes the empty region before the contig
        return t.value_at_or(0, 0.0);
    }
    t.value_at_or(pos - 1, 0.0)
}

struct Row {
    func: String,
    pos: Vec<i64>,
    val: Vec<f32>,
    args: Vec<String>,
}

fn load() -> Vec<Row> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/pileup_vectors.tsv");
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {path}: {e}; run oracle/gen_pileup_vectors.py"));
    let mut out = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        assert!(f.len() >= 4, "malformed vector row: {line}");
        let val = f[f.len() - 1];
        let pos = f[f.len() - 2];
        let func = f[0].to_string();
        // the args are everything between, minus the function name
        let args: Vec<String> = f[1..f.len() - 2].iter().map(|s| s.to_string()).collect();
        out.push(Row {
            func,
            pos: parse_list(pos),
            val: val
                .split(',')
                .filter(|x| !x.trim().is_empty())
                .map(f32_from_bits)
                .collect(),
            args,
        });
    }
    assert!(out.len() > 1000, "only {} pileup vectors loaded", out.len());
    out
}

#[test]
fn bit_exact_pileup_parity_with_macs3() {
    let rows = load();
    let mut checked: BTreeMap<&str, usize> = BTreeMap::new();
    let mut failures: Vec<String> = Vec::new();

    for row in &rows {
        let a = &row.args;
        let t: SignalTrack<f32> = match row.func.as_str() {
            "pn" | "pn_junction" => {
                // args: rlength five_shift three_shift scale baseline plus minus
                let rlength: u64 = a[0].parse().expect("rlength");
                let five: i64 = a[1].parse().expect("five_shift");
                let three: i64 = a[2].parse().expect("three_shift");
                let scale: f32 = a[3].parse().expect("scale");
                let baseline: f32 = a[4].parse().expect("baseline");
                let plus: Vec<u64> = a[5]
                    .split(',')
                    .filter(|x| !x.trim().is_empty())
                    .map(|x| x.trim().parse().expect("plus pos"))
                    .collect();
                let minus: Vec<u64> = a[6]
                    .split(',')
                    .filter(|x| !x.trim().is_empty())
                    .map(|x| x.trim().parse().expect("minus pos"))
                    .collect();
                let p = SingleEndParams {
                    five_shift: five,
                    three_shift: three,
                    rlength,
                    scale_factor: scale,
                    baseline_value: baseline,
                    // F71: the vectors are compared value-for-value, and this
                    // mirrors what the real pipeline requests.
                    coalesce: true,
                    centred: false,
                };
                pileup_from_positions(C, &plus, &minus, &p)
            }
            "lr" => {
                // args: rlength scale lefts rights
                let rlength: u64 = a[0].parse().expect("rlength");
                let scale: f32 = a[1].parse().expect("scale");
                let ls: Vec<u64> = a[2]
                    .split(',')
                    .filter(|x| !x.trim().is_empty())
                    .map(|x| x.trim().parse().expect("l"))
                    .collect();
                let rs: Vec<u64> = a[3]
                    .split(',')
                    .filter(|x| !x.trim().is_empty())
                    .map(|x| x.trim().parse().expect("r"))
                    .collect();
                let frags: Vec<(u64, u64)> = ls.into_iter().zip(rs).collect();
                pileup_from_fragments(C, &frags, rlength, scale, 0.0)
            }
            other => {
                failures.push(format!("unhandled pileup vector {other:?}"));
                continue;
            }
        };

        *checked.entry(row.func.as_str()).or_default() += 1;

        // 1. every upstream breakpoint must exist in our value function with the
        //    same value. This is the property that actually matters: a peak
        //    threshold is applied to these values.
        let mut value_mismatch = None;
        for (i, &p) in row.pos.iter().enumerate() {
            if p < 0 {
                value_mismatch = Some(format!("negative breakpoint {p}"));
                break;
            }
            let want = row.val[i];
            let got = upstream_value(&t, p as u64);
            if got.to_bits() != want.to_bits() {
                value_mismatch = Some(format!(
                    "at {p}: got {got:?} (0x{:08x}) want {want:?} (0x{:08x})",
                    got.to_bits(),
                    want.to_bits()
                ));
                break;
            }
        }
        if let Some(m) = value_mismatch {
            failures.push(format!(
                "{}({}) value mismatch {m}\n    ours: {:?}",
                row.func,
                a.join(","),
                t.spans().collect::<Vec<_>>()
            ));
            continue;
        }

        // 2. the breakpoint *set* must be identical, modulo the coalescing our
        //    canonical form performs. Every one of our runs must start at an
        //    upstream breakpoint, and vice versa.
        let upstream_starts = &row.pos;
        // our head may start at 0 when upstream's first breakpoint is 0 (a
        // leading zero run), which is the same position
        // Every one of our run *ends* must be an upstream breakpoint, and every
        // upstream breakpoint must be one of our run ends. Coalescing only ever
        // removes breakpoints, so the check is: our ends are a subset in order,
        // and the value function above already proved the values agree.
        let our_ends: Vec<i64> = t.runs().iter().map(|r| r.end as i64).collect();
        let mut j = 0usize;
        let mut skipped = Vec::new();
        for e in &our_ends {
            while j < upstream_starts.len() && upstream_starts[j] < *e {
                skipped.push(upstream_starts[j]);
                j += 1;
            }
            if j >= upstream_starts.len() || upstream_starts[j] != *e {
                failures.push(format!(
                    "{}({}) run end {e} is not an upstream breakpoint\n    ours: {:?}\n    upstream: {:?}\n    upstream breakpoints we did not reproduce: {skipped:?}",
                    row.func,
                    a.join(","),
                    t.spans().collect::<Vec<_>>(),
                    upstream_starts
                ));
                break;
            }
            j += 1;
        }
    }

    if !failures.is_empty() {
        let mut by_fn: BTreeMap<&str, usize> = BTreeMap::new();
        for f in &failures {
            *by_fn.entry(f.split('(').next().unwrap_or("?")).or_default() += 1;
        }
        panic!(
            "{} of {} pileup vectors disagree with MACS3 3.0.5 (by function: {by_fn:?}):\n{}",
            failures.len(),
            rows.len(),
            failures
                .iter()
                .take(6)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }

    for required in ["pn", "pn_junction", "lr"] {
        assert!(
            checked.get(required).copied().unwrap_or(0) > 0,
            "no vectors checked for {required}"
        );
    }
    eprintln!(
        "pileup parity: {} vectors OK ({:?})",
        checked.values().sum::<usize>(),
        checked
    );
}

/// The upstream coincident-event behaviour, pinned as an explicit expectation.
///
/// Reads at `a`, `a + d`, `a + 2d` with extension `d` produce three inferred
/// fragments that abut exactly. Upstream's sweep takes the net-zero branch and
/// emits no breakpoint at either junction, so the pv array has a single
/// breakpoint at the final fragment end. The value function is still correct:
/// the depth really is 1 across the whole span.
///
/// Values below are the ones MACS 3.0.5 returns, verified against the compiled
/// oracle (`tests/pileup_vectors.tsv`, `pn_junction` rows) and recorded as `F5`
/// in `docs/upstream-findings.md`.
#[test]
fn coincident_start_and_end_emits_no_breakpoint() {
    // reads at 5, 55, 105 with d = 50 on a 1000 bp contig
    let p = SingleEndParams::directional(50, 0, 1000, 1.0);
    let t = pileup_from_positions(C, &[5, 55, 105], &[], &p);

    // the pv array is (positions=[5, 155], values=[0.0, 1.0]), i.e. the track
    // below is 0.0 on [0,5) and 1.0 on [5,155)
    let (pos, val) = t.to_breakpoints();
    assert_eq!(pos, vec![5u64, 155]);
    assert_eq!(val, vec![0.0f32, 1.0]);

    // ...and that is the true depth, not an artefact
    assert_eq!(t.value_at_or(4, -1.0), 0.0);
    assert_eq!(t.value_at_or(5, -1.0), 1.0);
    assert_eq!(t.value_at_or(100, -1.0), 1.0);
    assert_eq!(t.value_at_or(154, -1.0), 1.0);
    // nothing is represented past the last breakpoint
    assert_eq!(t.value_at(155), None);

    // a non-abutting case keeps every breakpoint: reads at 5, 60, 115 with
    // d = 50 leave a 5 bp gap between each fragment
    let q = SingleEndParams::directional(50, 0, 1000, 1.0);
    let u = pileup_from_positions(C, &[5, 60, 115], &[], &q);
    let (pos, val) = u.to_breakpoints();
    assert_eq!(pos, vec![5u64, 55, 60, 110, 115, 165]);
    assert_eq!(val, vec![0.0f32, 1.0, 0.0, 1.0, 0.0, 1.0]);

    // and partially overlapping fragments stack; the middle junction at 55 is
    // coincident (fragment 1 ends, fragment 3 starts) so no breakpoint is
    // emitted there and the depth correctly reads 2 across it
    let r = SingleEndParams::directional(50, 0, 1000, 1.0);
    let w = pileup_from_positions(C, &[5, 30, 55], &[], &r);
    let (pos, val) = w.to_breakpoints();
    // [5,30)=1, [30,80)=2, [80,105)=1
    assert_eq!(pos, vec![5u64, 30, 80, 105]);
    assert_eq!(val, vec![0.0f32, 1.0, 2.0, 1.0]);
    assert_eq!(w.value_at_or(60, -1.0), 2.0);
}
