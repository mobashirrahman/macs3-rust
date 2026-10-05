use macs_core::ChromId;
use macs_peaks::callpeak::{build_qtable_from_with_hist, ChromSignals, CutoffParams};
use macs_rle::{Run, SignalTrack};

fn track(chrom: ChromId, runs: &[(u64, f32)]) -> SignalTrack<f32> {
    SignalTrack::from_runs_exact(
        chrom,
        0,
        runs.last().map_or(0, |(end, _)| *end),
        runs.iter()
            .map(|(end, value)| Run::new(*end, *value))
            .collect(),
    )
}

#[test]
fn narrow_histogram_reduction_matches_materialized_qtracks() {
    let chrom = ChromId(0);
    let signals = vec![ChromSignals {
        name: "chr1".into(),
        chrom,
        treat: track(chrom, &[(20, 1.0), (40, 2.0), (60, 1.0)]),
        ctrl: Some(track(chrom, &[(10, 0.5), (30, 1.0), (60, 0.5)])),
    }];

    let full = build_qtable_from_with_hist(&signals, false, 0, true, true, CutoffParams::off());
    let narrow = build_qtable_from_with_hist(&signals, false, 0, false, true, CutoffParams::off());

    assert_eq!(narrow.table.entries(), full.table.entries());
    assert_eq!(narrow.histogram, full.histogram);
    assert!(narrow.qtracks.is_empty());
    assert_eq!(full.qtracks.len(), 1);
}
