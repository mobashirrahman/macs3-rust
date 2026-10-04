//! Regressions for the paired-end peak-calling defects found against MACS3 3.0.5.
//!
//! Each case names the finding in `docs/upstream-findings.md` and pins the
//! behaviour that was wrong, so a future refactor of the p->q walk cannot silently
//! reintroduce it.
//!
//! The union pairing these interact with is covered by
//! `crates/macs-peaks/tests/paired_union_regressions.rs`.

use macs_score::{PScoreHistogram, PqTable};

/// **F179** — the p->q walk's rank `k` accumulates *signed* base-pair counts and
/// goes negative as soon as a negative span is added. `log10` of a negative is
/// `NaN`, `NaN > cutoff` is false, and upstream stores that `NaN` in its table.
///
/// The fixture is `sweep/gmini_mfrag_d1200_w180_noc_123` reduced to its shape: a
/// histogram whose highest score carries `-4959` base pairs.
#[test]
fn a_negative_rank_produces_nan_qscores_that_never_pass_a_cutoff() {
    let mut h = PScoreHistogram::new();
    // The lengths telescope to the last paired position upstream, so `N` is always
    // positive even though the *first* span is negative. 6000 - 4959 = 1041 here.
    h.add(6775.87, 6000);
    h.add(765.3187, -4959);
    let q = PqTable::from_histogram(&h);

    assert!(
        !q.qscore_or_zero(6775.87).is_nan(),
        "the top bucket is finite: k is still 1 there"
    );
    assert!(
        q.qscore_or_zero(6775.87) > 0.0,
        "summit q must be a real number"
    );
    // The negative span's own q is still finite -- it is the bins *below* it that
    // inherit the NaN, because `k` only goes negative after `k += l`.
    assert!(!q.qscore_or_zero(765.3187).is_nan());
}

/// The NaN window is the observable consequence: a lookup there must not compare
/// greater than any cutoff, which is what removes those positions from the
/// above-cutoff set (F177).
///
/// The shape mirrors `sweep/gmini_mfrag_d1200_w180_noc_123`: a large positive
/// span at the summit, then the negative first span that drives `k` below zero,
/// then enough further spans to carry it back above zero.
#[test]
fn nan_qscores_never_exceed_a_cutoff() {
    let mut h = PScoreHistogram::new();
    h.add(6775.87, 100);
    h.add(765.3187, -4959); // k: 1 -> 101 -> -4858
    h.add(700.0, 1); // k = -4858 -> NaN, then k = -4857
    h.add(650.0, 6000); // k = -4857 -> NaN, then k = 1143
    h.add(600.0, 1); // k = 1143 -> finite
    let q = PqTable::from_histogram(&h);

    assert!(
        q.qscore_or_zero(700.0).is_nan(),
        "the bin below the negative span must inherit the NaN"
    );
    assert!(
        q.qscore_or_zero(650.0).is_nan(),
        "the window extends while k stays negative"
    );
    // and the walk really did recover: the bin after the NaN window is a real
    // number again (its *value* is 0, which is F18 -- upstream always forces the
    // lowest p-score bucket to q = 0).
    assert!(!q.qscore_or_zero(600.0).is_nan());
    for cutoff in [-1.0f32, 0.0, 1.0, 1.301_03] {
        let score = q.qscore_or_zero(700.0);
        assert!(
            score.partial_cmp(&cutoff) != Some(std::cmp::Ordering::Greater),
            "a NaN q-score must never pass a cutoff (cutoff {cutoff})"
        );
    }
}
