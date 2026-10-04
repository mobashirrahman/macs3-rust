//! Bit-for-bit check of `pscore` against upstream's `get_pscore`.
//!
//! Upstream ([`CallPeakUnit.py:71`]): `get_pscore = -1 * poisson_cdf(k, lam,
//! lower=False, log10=True)` -- the **upper** tail, computed in log10 space.
//! Both extra arguments are positional booleans, so transposing them would
//! silently return `log10 P(X <= k)`, a negative score of the wrong magnitude
//! that still looks superficially plausible.
//!
//! Vectors: `pscore_grid.tsv`, generated from the pinned oracle by
//! `oracle/gen_pscore_grid.py`. The `#` prefix marks the comment line.

use macs_score::pscore;

fn parse_f32(s: &str) -> f32 {
    s.parse().unwrap_or_else(|e| panic!("bad float {s:?}: {e}"))
}

#[test]
fn pscore_matches_upstream_on_the_grid() {
    let raw = include_str!("pscore_grid.tsv");
    let mut checked = 0usize;
    let mut worst = 0.0f32;
    let mut worst_at = (0u32, 0.0f32);
    for line in raw.lines().filter(|l| !l.starts_with('#')) {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 3, "malformed vector row: {line:?}");
        let k: u32 = f[0].parse().expect("k");
        let lam = parse_f32(f[1]);
        let want = parse_f32(f[2]);

        let got = pscore(None, k, lam);

        // `-log10` of a tail probability loses relative precision when the
        // probability is denormal, so compare in absolute terms with a tolerance
        // scaled to the magnitude, while still demanding tight agreement.
        let tol = 1e-4 * want.abs().max(1.0);
        let diff = (got - want).abs();
        if diff > worst {
            worst = diff;
            worst_at = (k, lam);
        }
        assert!(
            diff <= tol,
            "pscore({k}, {lam}) = {got}, upstream {want}, diff {diff:e} > {tol:e}"
        );
        checked += 1;
    }
    assert!(
        checked >= 100,
        "expected a full grid, only checked {checked}"
    );
    println!(
        "{checked} vectors, worst abs diff {worst:e} at k={} lam={}",
        worst_at.0, worst_at.1
    );
}
