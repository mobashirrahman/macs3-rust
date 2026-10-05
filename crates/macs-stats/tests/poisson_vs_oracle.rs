//! Differential harness for [`macs_stats::poisson_cdf`] against the pinned oracle.
//!
//! Upstream's `poisson_cdf(k, lam, False, True)` is what produces every MACS
//! p-score, and it ends in `round(..., 5)` -- so the p-score is quantised to 1e-5
//! before any comparison. A disagreement of one unit in the 5th decimal is
//! therefore a byte difference in `*_peaks.xls` even though the pre-rounding value
//! may agree to 1e-12 (F191).
//!
//! The committed capture was generated with the pinned MACS3 3.0.5 public API.
//! Set `MACS_POIS_TSV` to check another capture from the same API.
//!
//! ```text
//! MACS_POIS_TSV=/tmp/pois.tsv cargo test --release -p macs-stats --test poisson_vs_oracle -- --ignored --nocapture
//! ```

#[test]
fn poisson_upper_tail_log10_matches_oracle() {
    let path = std::env::var("MACS_POIS_TSV").unwrap_or_else(|_| {
        format!(
            "{}/tests/data/poisson_oracle.tsv",
            env!("CARGO_MANIFEST_DIR")
        )
    });
    let text = std::fs::read_to_string(path).expect("read capture");
    let mut n = 0usize;
    let mut bad = 0usize;
    let mut first: Option<String> = None;
    for (line_number, line) in text
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 3, "bad row at capture line {}", line_number + 1);
        let k = f[0].parse::<u32>().expect("valid integer observation");
        let lam = f[1].parse::<f64>().expect("valid lambda");
        let want = f[2].parse::<f64>().expect("valid oracle result");
        assert!(lam.is_finite() && lam > 0.0 && want.is_finite());
        let got = macs_stats::poisson_cdf(k, lam, false, true).expect("oracle lambda is valid");
        n += 1;
        if got != want {
            bad += 1;
            if first.is_none() {
                first = Some(format!(
                    "k={k} lam={lam}\n  oracle   {want:?}\n  ours     {got:?}\n  \
                     py_round(ours, 5) = {:?}",
                    macs_stats::py_round(got, 5)
                ));
            }
        }
    }
    assert!(
        n > 0,
        "the Poisson oracle capture must contain observations"
    );
    eprintln!("points {n}, mismatches {bad}");
    if let Some(f) = first {
        eprintln!("first mismatch:\n{f}");
    }
    assert_eq!(bad, 0, "{bad}/{n} points disagree with the oracle");
}
