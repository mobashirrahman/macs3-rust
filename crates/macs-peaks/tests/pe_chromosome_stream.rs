use macs_peaks::callpeak::{run_callpeak_pe, run_callpeak_pe_chromosome, PeConfig};
use macs_track::FragTrackBuilder;

fn tracks() -> (macs_track::FragmentTrack, macs_track::FragmentTrack) {
    let mut treat = FragTrackBuilder::new();
    let mut ctrl = FragTrackBuilder::new();
    for (chrom, offset) in [(b"chr1".as_slice(), 0u64), (b"chr2".as_slice(), 10_000)] {
        for i in 0..80u64 {
            let start = offset + 100 + i * 11;
            treat.push(chrom, start, start + 140 + i % 7);
            let cstart = offset + 90 + i * 13;
            ctrl.push(chrom, cstart, cstart + 135 + i % 9);
        }
    }
    treat.finalize();
    ctrl.finalize();
    (treat.build(), ctrl.build())
}

fn counted_tracks() -> (macs_track::FragmentTrack, macs_track::FragmentTrack) {
    let mut treat = FragTrackBuilder::with_barcodes();
    let mut ctrl = FragTrackBuilder::with_barcodes();
    for (chrom, offset) in [(b"chr1".as_slice(), 0u64), (b"chr2".as_slice(), 10_000)] {
        for i in 0..80u64 {
            let start = offset + 100 + i * 11;
            treat.push_with_count(chrom, start, start + 140 + i % 7, 1 + (i % 4) as u32);
            let cstart = offset + 90 + i * 13;
            ctrl.push_with_count(chrom, cstart, cstart + 135 + i % 9, 1 + (i % 3) as u32);
        }
    }
    treat.finalize();
    ctrl.finalize();
    (treat.build(), ctrl.build())
}

fn assert_same_track(left: &macs_rle::SignalTrack<f32>, right: &macs_rle::SignalTrack<f32>) {
    assert_eq!((left.start(), left.end()), (right.start(), right.end()));
    assert_eq!(left.runs().len(), right.runs().len());
    for (a, b) in left.runs().iter().zip(right.runs()) {
        assert_eq!(a.end, b.end);
        assert_eq!(a.value.to_bits(), b.value.to_bits());
    }
}

#[test]
fn selected_chromosome_signals_match_the_full_paired_pipeline() {
    let (treat, ctrl) = tracks();
    let cfg = PeConfig {
        tsize: 143.0,
        tsize_exact: 143.75,
        gsize: 25_000.0,
        slocal: 500,
        llocal: 1_000,
        qvalue: 0.05,
        call_summits: false,
        broad: false,
        broad_cutoff: 0.1,
        nolambda: false,
        scaleto_large: false,
    };
    let all = run_callpeak_pe(&treat, Some(&ctrl), &cfg);
    assert_eq!(all.signals.len(), 2);

    for expected in &all.signals {
        let one = run_callpeak_pe_chromosome(&treat, Some(&ctrl), &cfg, expected.name.as_bytes());
        assert_eq!(one.d, all.d);
        assert_eq!(one.paired_boundaries, all.paired_boundaries);
        assert_eq!(one.lambda_bg.to_bits(), all.lambda_bg.to_bits());
        assert_eq!(one.coord_shift, all.coord_shift);
        assert_eq!(one.signals.len(), 1);
        let actual = &one.signals[0];
        assert_eq!(actual.name, expected.name);
        assert_same_track(&actual.treat, &expected.treat);
        assert_same_track(
            actual.ctrl.as_ref().unwrap(),
            expected.ctrl.as_ref().unwrap(),
        );
    }
}

#[test]
fn selected_chromosome_signals_match_for_counted_fragments() {
    let (treat, ctrl) = counted_tracks();
    let cfg = PeConfig {
        tsize: 143.0,
        tsize_exact: 143.75,
        gsize: 25_000.0,
        slocal: 500,
        llocal: 1_000,
        qvalue: 0.05,
        call_summits: false,
        broad: false,
        broad_cutoff: 0.1,
        nolambda: false,
        scaleto_large: false,
    };
    let all = run_callpeak_pe(&treat, Some(&ctrl), &cfg);
    assert_eq!(all.signals.len(), 2);

    for expected in &all.signals {
        let one = run_callpeak_pe_chromosome(&treat, Some(&ctrl), &cfg, expected.name.as_bytes());
        assert_eq!(one.d, all.d);
        assert_eq!(one.paired_boundaries, all.paired_boundaries);
        assert_eq!(one.lambda_bg.to_bits(), all.lambda_bg.to_bits());
        assert_eq!(one.coord_shift, all.coord_shift);
        assert_eq!(one.signals.len(), 1);
        let actual = &one.signals[0];
        assert_eq!(actual.name, expected.name);
        assert_same_track(&actual.treat, &expected.treat);
        assert_same_track(
            actual.ctrl.as_ref().unwrap(),
            expected.ctrl.as_ref().unwrap(),
        );
    }
}
