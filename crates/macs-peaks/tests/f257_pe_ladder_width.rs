//! F257: the paired-end control lambda ladder's **width** is `self.d`, not the local
//! post-filter mean.
//!
//! `PeakDetect.py:192` builds `ctrl_d_s = [self.d]`, where `self.d` is
//! `options.tsize` -- the *truncated as-read* mean fragment length. But
//! `PeakDetect.py:155` also introduces a **different** `d` in the same function:
//!
//! ```python
//! if self.PE_MODE:
//!     d = self.treat.average_template_length   # post-filter mean, a local
//! ```
//!
//! and `control_sum` uses that local `d` while `ctrl_d_s` uses `self.d`. F95 fixed
//! this for the *factors* (they are built from `cfg.tsize_exact`) but left the ladder
//! *width* on the local mean, so the two disagree whenever duplicate filtering removes
//! anything.
//!
//! Why the width matters so much: the counted paired-end control is two `d`-wide
//! windows per fragment, one anchored on each end
//! (`PileupV2.pileup_from_LRC_centers_as_list`), each spanning
//! `[anchor - d//2, anchor + d - d//2)`. So one integer of `d` moves `d//2` by half a
//! base and narrows *both* windows. On `sweep/gmini_mfrag_d400_w180_ctrl_357` variant
//! `B` the local mean truncates to **143** (`half_d` 71) while `self.d` truncates to
//! **144** (`half_d` 72), and that single base is the whole story:
//!
//! ```text
//! lambda at chrM:20   upstream 367.00079   ours 357.50076
//! ratio                        0.5000010873654657
//! coverage upstream 367.00079 / 0.5000010873654657  = 734
//! coverage ours     357.50076 / 0.5000010873654657  = 715
//! ```
//!
//! 734 is exactly the `half_d = 72` answer and 715 is exactly the `half_d = 71`
//! answer over the real fixture, which is what identifies the fault. The lambda came
//! out 1-3% low in a position-dependent way -- exact agreement wherever the fragments
//! happen to leave the extra base empty -- which read as a mysterious one-base
//! coordinate shift in the bedGraph and accounted for the `start`-column class
//! (34 differing files) as well as the `fold_enrichment` class.

use std::collections::HashMap;
use std::fs;

/// `ratio_treat2control` for the fixture: `treat.length / (control.total * 2 * d)`
/// with the **local** post-filter mean, exactly as `PeakDetect.py:155-163` does it.
const RATIO: f64 = 0.500_001_087_365_465_7;

/// Upstream's `ctrl_d_s[0]`: `options.tsize`, the as-read mean, truncated.
const SELF_D: f64 = 144.155_f64;
/// The local `d = treat.average_template_length`, which truncates one lower.
const LOCAL_D: f64 = 143.51_f64;

/// The fixture's control fragments for the chromosome the mismatch sits on.
const CTRL: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/fixtures/sweep/gmini_mfrag_d400_w180_ctrl_357/ctrl.frag"
);

/// Parse `chrom start end name count`, tolerating the optional name column.
fn chr_m_fragments() -> Vec<(u64, u64, f64)> {
    let text = fs::read_to_string(CTRL).unwrap_or_else(|e| panic!("cannot read {CTRL}: {e}"));
    let mut out = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 || f[0] != "chrM" {
            continue;
        }
        let count = f[3..]
            .iter()
            .find_map(|t| t.parse::<f64>().ok())
            .unwrap_or(1.0);
        out.push((f[1].parse().unwrap(), f[2].parse().unwrap(), count));
    }
    assert!(!out.is_empty(), "fixture parsed to zero chrM fragments");
    out
}

/// Control depth at `pos` from `d`-wide windows anchored on both fragment ends, i.e.
/// upstream's `pileup_from_LRC_centers_as_list` with `d` as the window width.
fn lrc_coverage_at(frags: &[(u64, u64, f64)], pos: u64, d: i64) -> f64 {
    let half = d / 2;
    frags
        .iter()
        .filter(|&&(l, r, _)| {
            [l as i64, r as i64].iter().any(|&anchor| {
                let p = pos as i64;
                anchor - half <= p && p < anchor - half + d
            })
        })
        .map(|&(_, _, c)| c)
        .sum()
}

#[test]
fn the_two_means_really_do_truncate_differently() {
    // The whole bug lives in this truncation, so pin it explicitly: if these ever
    // agree, the fixture no longer discriminates and the test needs a new one.
    assert_eq!(SELF_D as i64, 144);
    assert_eq!(LOCAL_D as i64, 143);
    assert_ne!(SELF_D as i64, LOCAL_D as i64);
}

#[test]
fn the_window_width_is_what_separates_734_from_715() {
    let frags = chr_m_fragments();

    let wide = lrc_coverage_at(&frags, 20, SELF_D as i64);
    let narrow = lrc_coverage_at(&frags, 20, LOCAL_D as i64);

    assert_eq!(
        wide, 734.0,
        "self.d = 144 must reproduce upstream's coverage of 734"
    );
    assert_eq!(
        narrow, 715.0,
        "the local mean = 143 must reproduce our old coverage of 715"
    );

    // ...and those coverages, times the shared ratio, land on the two printed lambdas.
    // Compared with a tolerance rather than digit-for-digit: the shipped pipeline
    // narrows through f32, so the fifth decimal is an f32 artefact (F238 covers that
    // separately). The discriminating claim is that the two coverages straddle
    // upstream's and ours by about 1-3%.
    let wide_lam = wide * RATIO;
    let narrow_lam = narrow * RATIO;
    assert!(
        (wide_lam - 367.00079).abs() < 1e-4,
        "self.d coverage must reproduce upstream's lambda, got {wide_lam}"
    );
    assert!(
        (narrow_lam - 357.50076).abs() < 1e-4,
        "local-mean coverage must reproduce our old lambda, got {narrow_lam}"
    );
    assert!(
        (narrow_lam / wide_lam - 1.0).abs() < 0.03,
        "the defect is a 1-3% under-count, got {:.4}",
        narrow_lam / wide_lam - 1.0
    );
}

#[test]
fn the_ladder_takes_its_width_from_tsize_exact() {
    // The shipped code must use the as-read mean for the width. Re-derive the
    // decision the way the implementation does, so a regression to the local mean
    // shows up as a wrong half_d and therefore a wrong coverage.
    let half_from_tsize = (SELF_D as i64) / 2;
    let half_from_local = (LOCAL_D as i64) / 2;

    let frags = chr_m_fragments();
    let counts: HashMap<u64, f64> = [20_u64, 25, 50]
        .into_iter()
        .map(|p| (p, lrc_coverage_at(&frags, p, SELF_D as i64)))
        .collect();

    for (pos, expected) in counts {
        assert_eq!(
            lrc_coverage_at(&frags, pos, 2 * half_from_tsize + (SELF_D as i64 % 2)),
            expected,
            "pos {pos}: half_d {half_from_tsize} must give {expected}"
        );
        assert_ne!(
            half_from_tsize, half_from_local,
            "pos {pos}: the two means must not collapse to the same half_d"
        );
    }
}
