//! Peak calling.
//!
//! Segments a chromosome's above-cutoff signal into peak regions, picks a summit
//! for each, and — under `--call-summits` — emits one peak per sub-peak summit.
//!
//! The transcription is from `MACS3.Signal.CallPeakUnit`
//! (`__chrom_call_peak_using_certain_criteria`, `__close_peak_wo_subpeaks`,
//! `__close_peak_with_subpeaks`) and `MACS3.Signal.SignalProcessing`
//! (`maxima`, `enforce_peakyness`, `savitzky_golay_order2_deriv1`). See
//! [`regions`] for the conventions that decide peak coordinates — notably that the
//! signal arrays are end-indexed (F5) and that a tied summit resolves to the
//! *lower median* of the candidate midpoints, not the first maximum.

#![forbid(unsafe_code)]
#![deny(missing_debug_implementations)]

pub mod caller;
pub mod callpeak;
pub mod driver;
pub mod hmm_regions;
mod merge;
mod peakyness;
mod precomputes;
mod regions;
pub mod sg;

pub use caller::{
    above_cutoff, cal_fold, cal_pscore, cal_qscore, cal_score, cal_subtraction, call_chromosome,
    call_chromosome_paired, Criterion,
};
pub use driver::{
    combine_control_scales, nolambda_control, pair_for_chromosome, pair_treat_ctrl, LambdaScale,
    LambdaScales, PairedSignal,
};
pub use merge::{
    over_max_tracks, over_two_pv_array, over_two_pv_array_track, pointwise_max, track_from_pv,
    Reducer,
};
pub use peakyness::{enforce_peakyness, hard_clip, internal_minima, is_valid_peak, too_flat};
pub use precomputes::{
    accumulate_histogram, chromosome_cutoff_stats, cutoff_ladder, pre_computes,
    pre_computes_tracks, CutoffStats,
};
pub use regions::{
    call_peaks_chromosome, call_peaks_chromosome_with_p, close_peak_for_broad_region,
    close_peak_with_subpeaks, close_peak_wo_subpeaks, close_peak_wo_subpeaks_with_p,
    segment_regions, CallParams, Chunk, Peak, Reject, ScoreKind,
};
pub use sg::{deriv1_coefficients, maxima, savitzky_golay_order2_deriv1};
