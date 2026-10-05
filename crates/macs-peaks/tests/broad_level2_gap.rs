//! `--broad`'s level-2 merge gap is `int(opt.tsize * 4)`, not `int(opt.tsize) * 4`.
//!
//! `PeakDetect.maxgap` is `opt.maxgap or opt.tsize`
//! (`PeakDetect.py:71-74`), and `opt.tsize` is a Python **float** -- the mean
//! fragment length, `Parser.py:1496`. `call_broadpeaks` is then called with
//! `lvl1_max_gap=self.maxgap, lvl2_max_gap=self.maxgap*4`
//! (`PeakDetect.py:259-262`), and only there does Cython truncate, at the
//! `cython.int` parameters of `call_broadpeaks` (`CallPeakUnit.py:1793-1794`).
//! With a mean fragment length of `253.5` that is `int(253.5) = 253` for level 1 and
//! `int(1014.0) = 1014` for level 2 -- not `253 * 4 = 1012`. Two bases of difference,
//! which is one merged-versus-split broad region: on
//! `CTCF_PE_ChIP_chr22_50k.bedpe.gz` (`-f BEDPE`, mean fragment size 253.5) the region
//! at `chr22:18187588-18189977` has an internal gap of `1013` bases, so upstream's level-2
//! walk keeps it whole and a `max_gap * 4` gap splits it in two.
//!
//! This drives `call_chromosome` with that geometry directly: two 1000 bp plateaus of
//! treatment pileup separated by a 1013 bp dip, which level 1 (gap 253) splits and
//! level 2 must merge.

use macs_core::ChromId;
use macs_peaks::callpeak::{call_chromosome, ChromCall};
use macs_rle::{Run, SignalTrack};
use macs_score::{PScoreCache, PqTable};

const C: ChromId = ChromId(0);

/// Plateau [1000, 2000), dip [2000, 3013), plateau [3013, 4013).
fn treatment() -> SignalTrack<f32> {
    SignalTrack::from_runs(
        C,
        0,
        5000,
        vec![
            Run::new(1000, 0.0f32),
            Run::new(2000, 40.0),
            Run::new(3013, 0.0),
            Run::new(4013, 40.0),
            Run::new(5000, 0.0),
        ],
    )
}

/// A flat local lambda over the whole contig, so every plateau base clears the cutoff.
fn control() -> SignalTrack<f32> {
    SignalTrack::from_runs(C, 0, 5000, vec![Run::new(5000, 1.0f32)])
}

/// `int(tsize) = 253` and `int(tsize * 4) = 1014` for the mean fragment length 253.5.
fn call(broad_max_gap: macs_core::Coord) -> Vec<macs_peaks::callpeak::Called> {
    let treat = treatment();
    let ctrl = control();
    let empty = SignalTrack::empty(C, 0, 5000);
    let table = PqTable::empty();
    let cc = ChromCall {
        name: "chr1",
        chrom: C,
        treat: &treat,
        ctrl: Some(&ctrl),
        clamp_floor: 0,
        zero_coord: 0,
        qtrack: &empty,
        table: &table,
        d: 500,
        max_gap: 253,
        broad_max_gap,
        p_cutoff: Some(1.0),
        qvalue: 1.0,
        broad: true,
        broad_cutoff: 0.1,
        call_summits: false,
        lambda_bg: 1.0,
    };
    let mut cache = PScoreCache::new();
    call_chromosome(&cc, &mut cache)
}

#[test]
fn level_two_merges_a_gap_of_int_tsize_times_four() {
    let peaks = call(1014);
    assert_eq!(
        peaks.len(),
        1,
        "the 1013 bp dip is <= int(253.5 * 4) = 1014, so it stays one region"
    );
    assert_eq!(peaks[0].peak.start, 1001);
    assert_eq!(peaks[0].peak.end, 4013);
    // both strong (level-1) sub-peaks are reported inside the broad region
    assert_eq!(peaks[0].lvl1, vec![(1001, 2000), (3014, 4013)]);
}

#[test]
fn level_two_splits_a_gap_one_base_over_int_tsize_times_four() {
    // `max_gap * 4 = 253 * 4 = 1012` is the value the port used to compute, and it is
    // two bases short of upstream's 1014: 1013 > 1012 splits where upstream merges.
    let peaks = call(1012);
    assert_eq!(peaks.len(), 2, "1013 > 1012 splits the broad region in two");
    assert_eq!((peaks[0].peak.start, peaks[0].peak.end), (1001, 2000));
    assert_eq!((peaks[1].peak.start, peaks[1].peak.end), (3014, 4013));
}
