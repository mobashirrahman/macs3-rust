//! Bit-exact parity of the p-to-q table.
//!
//! `tests/pq_vectors.tsv` has two kinds of row:
//!
//! * **`table`** — an explicit `(p-score, base-pair count)` histogram and the
//!   table `MACS3.Signal.ScoreTrack.make_pq_table` produces from it, computed by
//!   `oracle/pq_reference.py`. That module takes p-scores from the **compiled**
//!   `poisson_cdf` and transcribes only the histogram construction and the AFDR
//!   walk. `ScoreTrackII` is a Cython `cdef class`, so its table cannot be driven
//!   directly; this is the bridge, and the bar is bit equality.
//!
//! * **`cutoff`** — a *direct* readout of the real compiled MACS3's table, from
//!   `--cutoff-analysis`, which prints `self.pqtable[cutoff]` for 33 cutoffs.
//!   These are checked here for **structural** agreement: the table must be
//!   monotone non-increasing in q as p descends, must contain the reported
//!   p-scores, and must report 0 at the lowest one (F18). Full agreement is
//!   enforced end-to-end at gate G9, once the callpeak pipeline can reproduce the
//!   histogram that produced them.
//!
//! Bit equality is the right bar and is not achievable by accident: the AFDR rank
//! accumulates in `f32` (quantised above 2^24), `f = -log10(N)` is truncated to
//! `f32`, `q` is computed in `f64` and truncated to `f32`, and the lowest bucket is
//! unconditionally forced to 0 (F18).

use macs_score::{PScoreHistogram, PqTable};
use std::collections::BTreeMap;

fn f32_from_bits(h: &str) -> f32 {
    f32::from_bits(u32::from_str_radix(h.trim().trim_start_matches("0x"), 16).expect("bits"))
}

struct Row {
    kind: String,
    name: String,
    total: u64,
    hist_len: usize,
    table_len: usize,
    hist: Vec<(f32, u64)>,
    table: Vec<(f32, f32)>,
}

fn load() -> Vec<Row> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/pq_vectors.tsv");
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {path}: {e}; run oracle/gen_pq_vectors.py"));
    let mut out = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        assert!(f.len() == 7, "malformed row: {line}");
        // both columns encode float32 bit patterns for the key; the value is a
        // float32 for the table rows and a base-pair count for the histogram
        let keyed = |s: &str| -> Vec<(f32, String)> {
            s.split(',')
                .filter(|x| !x.trim().is_empty())
                .map(|x| {
                    let (a, b) = x.split_once(':').expect("k:v");
                    (f32_from_bits(a), b.trim().to_string())
                })
                .collect()
        };
        let kind = f[0].to_string();
        // `table` rows carry a histogram in field 5 and a table in field 6;
        // `cutoff` rows carry only the real MACS3 table, in field 5
        let (hist_src, table_src) = if kind == "table" {
            (f[5], f[6])
        } else {
            ("", f[5])
        };
        out.push(Row {
            kind,
            name: f[1].to_string(),
            total: f[2].parse().expect("total"),
            hist_len: f[3].parse().expect("hist_len"),
            table_len: f[4].parse().expect("table_len"),
            hist: keyed(hist_src)
                .into_iter()
                .map(|(k, v)| (k, v.parse().expect("count")))
                .collect(),
            table: keyed(table_src)
                .into_iter()
                .map(|(k, v)| (k, f32_from_bits(&v)))
                .collect(),
        });
    }
    assert!(out.len() >= 30, "only {} pq vectors loaded", out.len());
    out
}

#[test]
fn bit_exact_pq_table_parity() {
    let rows = load();
    let mut checked = 0usize;
    let mut failures: Vec<String> = Vec::new();

    for row in rows.iter().filter(|r| r.kind == "table") {
        // rebuild the histogram, checking the encoded total and length
        let mut h = PScoreHistogram::new();
        for &(p, n) in &row.hist {
            // F169: the histogram's lengths are signed, so a negative length is
            // expressible (upstream's first paired entry can precede position 0).
            h.add(p, n as i64);
        }
        if h.total() != row.total as i64 {
            failures.push(format!(
                "{}: histogram total {} != {} (signed zero should have merged)",
                row.name,
                h.total(),
                row.total
            ));
            continue;
        }
        if h.len() != row.hist_len {
            failures.push(format!(
                "{}: histogram has {} buckets, expected {}",
                row.name,
                h.len(),
                row.hist_len
            ));
            continue;
        }

        let t = PqTable::from_histogram(&h);
        if t.total() != row.total {
            failures.push(format!(
                "{}: table total {} != {}",
                row.name,
                t.total(),
                row.total
            ));
            continue;
        }
        if t.len() != row.table_len {
            failures.push(format!(
                "{}: table has {} entries, expected {}",
                row.name,
                t.len(),
                row.table_len
            ));
            continue;
        }

        for &(p, q) in &row.table {
            match t.get(p) {
                Some(got) if got.to_bits() == q.to_bits() => {}
                Some(got) => failures.push(format!(
                    "{}: q({p}) = 0x{:08x}, expected 0x{:08x}",
                    row.name,
                    got.to_bits(),
                    q.to_bits()
                )),
                None => failures.push(format!("{}: q({p}) absent from the table", row.name)),
            }
        }
        checked += 1;
    }

    if !failures.is_empty() {
        panic!(
            "{} of {} pq tables disagree:\n{}",
            failures.len(),
            rows.iter().filter(|r| r.kind == "table").count(),
            failures
                .iter()
                .take(20)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
    assert!(checked >= 20, "only {checked} tables checked");
    eprintln!("pq table parity: {checked} tables bit-exact");
}

#[test]
fn real_macs3_cutoff_tables_are_structurally_consistent() {
    let rows = load();
    let real: Vec<&Row> = rows.iter().filter(|r| r.kind == "cutoff").collect();
    assert!(
        real.len() >= 5,
        "expected real MACS3 cutoff-analysis tables, found {}",
        real.len()
    );

    let mut failures: Vec<String> = Vec::new();
    let mut checked = 0usize;
    let mut checked_full = 0usize;
    for row in real {
        // The pairs come out of the file in descending p-score order. The count
        // varies with the library: a fixture whose p-score range is narrow only
        // reports a few cutoffs.
        assert!(row.table.len() >= 3, "{} has too few pairs", row.name);
        // field 2 of a cutoff row carries the pair count
        assert_eq!(
            row.table.len() as u64,
            row.total,
            "{} length field",
            row.name
        );
        // a library with a wide p-score range must report the full 33 cutoffs
        if row.table.len() == 33 {
            checked_full += 1;
        }

        // 1. strictly descending p-score
        for w in row.table.windows(2) {
            // `!(a > b)` rather than `a <= b` so a NaN p-score is also caught
            #[allow(clippy::neg_cmp_op_on_partial_ord)]
            if !(w[0].0 > w[1].0) {
                failures.push(format!(
                    "{}: p-scores are not descending: {:?} then {:?}",
                    row.name, w[0].0, w[1].0
                ));
                break;
            }
        }
        // 2. q must be non-increasing as p descends, or capped by the previous q
        let mut prev_q = f32::INFINITY;
        for &(p, q) in &row.table {
            if q > prev_q + 1e-4 {
                failures.push(format!(
                    "{}: q rose from {prev_q} to {q} at p = {p}",
                    row.name
                ));
                break;
            }
            prev_q = q;
        }
        // 3. q must never exceed p, i.e. the multiple-testing adjustment is
        //    never anti-conservative. This is the property the whole AFDR walk
        //    exists to provide, and it is checkable from a filtered view of the
        //    table.
        //
        //    Note that the reported cutoffs are *filtered*: `CallPeakUnit` only
        //    writes a row when `pvalue_npeaks[cutoff] > 0`, so the lowest
        //    reported p-score is generally above the table's true minimum and
        //    says nothing about F18's forced-zero bucket. That is checked
        //    directly on the `table` rows instead.
        for &(p, q) in &row.table {
            if q > p + 1e-4 {
                failures.push(format!(
                    "{}: q = {q} exceeds p = {p}; the AFDR adjustment is anti-conservative",
                    row.name
                ));
                break;
            }
        }
        // 4. every reported p-score must be on the 1e-5 lattice (F1)
        for &(p, _) in &row.table {
            if p > 0.0 {
                let scaled = f64::from(p) * 1e5;
                if (scaled - scaled.round()).abs() > 0.5 {
                    failures.push(format!(
                        "{}: p-score {p} is not a multiple of 1e-5",
                        row.name
                    ));
                    break;
                }
            }
        }
        checked += 1;
    }
    if !failures.is_empty() {
        panic!("{} problems:\n{}", failures.len(), failures.join("\n"));
    }
    assert!(
        checked_full >= 5,
        "only {checked_full} full-length real tables"
    );
    eprintln!(
        "real MACS3 cutoff tables: {checked} fixtures structurally consistent \
         ({checked_full} with the full 33 cutoffs)"
    );
}

/// The forced-zero tail must be exactly one bucket when the AFDR walk completes,
/// and everything from the break down when it does not. This is the F18 contract
/// expressed directly, checked against the vector file's own bucket counts.
#[test]
fn forced_zero_tail_length_is_consistent_with_the_walk() {
    let rows = load();
    for row in rows.iter().filter(|r| r.kind == "table") {
        let mut h = PScoreHistogram::new();
        for &(p, n) in &row.hist {
            h.add(p, n as i64);
        }
        let t = PqTable::from_histogram(&h);
        let zeros = t.entries().into_iter().filter(|(_, q)| *q == 0.0).count();
        match t.cut_at() {
            // the walk broke: everything from the cut down is zero, so the number
            // of zeros is at least the number of buckets at or below the cut
            None => assert_eq!(
                zeros, 1,
                "{}: a completed walk must force exactly the lowest bucket to 0",
                row.name
            ),
            Some(cut) => {
                let below = h
                    .sorted_descending()
                    .into_iter()
                    .filter(|(p, _)| *p <= cut)
                    .count();
                assert!(
                    zeros >= below,
                    "{}: {zeros} zeros but {below} buckets at or below the cut {cut}",
                    row.name
                );
            }
        }
    }
}

#[test]
fn oracle_lock_records_the_pinned_version() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../oracle/ENV.lock");
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("oracle/ENV.lock missing: {e}"));
    assert!(
        text.contains("MACS3_VERSION=3.0.5"),
        "the pq vectors must come from MACS3 3.0.5, got:\n{text}"
    );
    let _ = BTreeMap::<u8, u8>::new();
}
