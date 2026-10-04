//! Differential harness for [`macs_stats::poisson_cdf`] against the pinned oracle.
//!
//! Upstream's `poisson_cdf(k, lam, False, True)` is what produces every MACS
//! p-score, and it ends in `round(..., 5)` -- so the p-score is quantised to 1e-5
//! before any comparison. A disagreement of one unit in the 5th decimal is
//! therefore a byte difference in `*_peaks.xls` even though the pre-rounding value
//! may agree to 1e-12 (F191).
//!
//! Needs a capture from the oracle, so it is ignored by default:
//!
//! ```text
//! MACS_POIS_TSV=/tmp/pois.tsv cargo test --release -p macs-stats --test poisson_vs_oracle -- --ignored --nocapture
//! ```

#[test]
#[ignore = "needs a capture from the pinned oracle"]
fn poisson_upper_tail_log10_matches_oracle() {
    let path = std::env::var("MACS_POIS_TSV").expect("MACS_POIS_TSV");
    let text = std::fs::read_to_string(path).expect("read capture");
    let mut n = 0usize;
    let mut bad = 0usize;
    let mut first: Option<String> = None;
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        let (k, lam, want, pre) = (
            f[0].parse::<u32>().unwrap(),
            f[1].parse::<f64>().unwrap(),
            f[2].parse::<f64>().unwrap(),
            f[3].parse::<f64>().unwrap(),
        );
        let got = macs_stats::poisson_cdf(k, lam, false, true).expect("oracle lambda is valid");
        n += 1;
        if got != want {
            bad += 1;
            if first.is_none() {
                first = Some(format!(
                    "k={k} lam={lam}\n  oracle   {want:?}\n  ours     {got:?}\n  \
                     pre-round value (oracle side) {pre:?}\n  \
                     py_round(ours, 5) = {:?}",
                    macs_stats::py_round(got, 5)
                ));
            }
        }
    }
    eprintln!("points {n}, mismatches {bad}");
    if let Some(f) = first {
        eprintln!("first mismatch:\n{f}");
    }
    assert_eq!(bad, 0, "{bad}/{n} points disagree with the oracle");
}
