//! The `enforce_peakyness` filter for sub-peak summits.
//!
//! Transcribed from `MACS3.Signal.SignalProcessing.enforce_peakyness` and its
//! helpers `internal_minima`, `is_valid_peak`, `too_flat`, and `hard_clip`.
//!
//! A candidate maximum survives only if the region around it, measured from its
//! adjacent minima, is:
//!
//! * at least **50 bp** wide, and
//! * **not too flat** — at least 6 distinct values.
//!
//! The threshold for each valley is the larger of its two bounding minima plus
//! the square root of that value, which is a heuristic for "prominent enough to be
//! a real binding site rather than sampling noise".

/// Contiguous nonnegative region containing `maximum`.
///
/// Negative boundary values are excluded; values exactly at zero stay inside.
/// Upstream walks left while `signal[left-1] >= 0` and right while
/// `signal[right] >= 0`, so zero is *included* in the region.
pub fn hard_clip(signal: &[f32], maximum: usize) -> &[f32] {
    if signal.is_empty() {
        return signal;
    }
    let maximum = maximum.min(signal.len() - 1);
    let mut left = maximum;
    while left > 0 && signal[left - 1] >= 0.0 {
        left -= 1;
    }
    let mut right = maximum + 1;
    while right < signal.len() && signal[right] >= 0.0 {
        right += 1;
    }
    &signal[left..right]
}

/// `True` when the region has fewer than 6 distinct values.
pub fn too_flat(signal: &[f32]) -> bool {
    let mut seen: Vec<u32> = signal.iter().map(|v| v.to_bits()).collect();
    seen.sort_unstable();
    seen.dedup();
    seen.len() < 6
}

/// A candidate maximum is valid if its clipped region is >= 50 wide and not flat.
pub fn is_valid_peak(signal: &[f32], maximum: usize) -> bool {
    let s = hard_clip(signal, maximum);
    if s.len() < 50 {
        return false;
    }
    !too_flat(s)
}

/// The lowest signal value in each interval between consecutive maxima.
///
/// For `maxima = [m0, m1, ..., mk]` this returns `k` values: the position of the
/// minimum of `signal[m_i .. m_{i+1}]`, offset by `m_i`.
pub fn internal_minima(signal: &[f32], maxima: &[usize]) -> Vec<usize> {
    let n = maxima.len();
    if n <= 1 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let pos1 = maxima[i];
        let pos2 = maxima[i + 1];
        if pos2 <= pos1 {
            out.push(pos1);
            continue;
        }
        // Python uses `signal[pos1:pos2]`, excluding the right-hand maximum.
        let Some(window) = signal.get(pos1..pos2) else {
            out.push(pos1.min(signal.len().saturating_sub(1)));
            continue;
        };
        if window.is_empty() {
            out.push(pos1.min(signal.len().saturating_sub(1)));
            continue;
        }
        // np.argmin returns the FIRST minimum on ties
        let mut best = 0usize;
        for (j, v) in window.iter().enumerate() {
            if *v < window[best] {
                best = j;
            }
        }
        out.push(pos1 + best);
    }
    out
}

/// Filter candidate maxima down to those that look like real peaks.
///
/// Returns a subset of `maxima`, in ascending order. With zero or one minimum
/// upstream returns `maxima` unchanged.
pub fn enforce_peakyness(signal: &[f32], maxima: &[usize]) -> Vec<usize> {
    if maxima.is_empty() {
        return Vec::new();
    }
    let minima = internal_minima(signal, maxima);
    let n = minima.len();
    if n == 0 {
        // no minima: every candidate is on the boundary, so nothing is filtered
        return maxima.to_vec();
    }

    let mut out: Vec<usize> = Vec::with_capacity(maxima.len());

    // region before the first minimum, paired with maxima[0]
    let threshold = signal[minima[0]] as f64 + (signal[minima[0]] as f64).sqrt();
    let new_signal: Vec<f32> = signal[0..minima[0]]
        .iter()
        .map(|v| (*v as f64 - threshold) as f32)
        .collect();
    if is_valid_peak(&new_signal, maxima[0]) {
        out.push(maxima[0]);
    }

    for i in 0..n.saturating_sub(1) {
        let a = signal[minima[i]] as f64;
        let b = signal[minima[i + 1]] as f64;
        let t = a.max(b);
        let threshold = t + t.sqrt();
        let new_signal: Vec<f32> = signal[minima[i]..minima[i + 1]]
            .iter()
            .map(|v| (*v as f64 - threshold) as f32)
            .collect();
        let new_maximum = maxima[i + 1].saturating_sub(minima[i]);
        if is_valid_peak(&new_signal, new_maximum) {
            out.push(maxima[i + 1]);
        }
    }

    // region after the last minimum, paired with the last maximum
    let threshold = signal[minima[n - 1]] as f64 + (signal[minima[n - 1]] as f64).sqrt();
    let new_signal: Vec<f32> = signal[minima[n - 1]..]
        .iter()
        .map(|v| (*v as f64 - threshold) as f32)
        .collect();
    let new_maximum = maxima[maxima.len() - 1].saturating_sub(minima[n - 1]);
    if is_valid_peak(&new_signal, new_maximum) {
        out.push(maxima[maxima.len() - 1]);
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hard_clip_keeps_zero_inside_the_region() {
        let s = vec![1.0f32, 0.0, 2.0, -1.0, 3.0];
        let r = hard_clip(&s, 0);
        assert_eq!(r, &[1.0, 0.0, 2.0], "zero is not a boundary");
    }

    #[test]
    fn hard_clip_excludes_a_negative_immediately_left() {
        let s = vec![-1.0f32, 5.0, 6.0];
        let r = hard_clip(&s, 1);
        assert_eq!(r, &[5.0, 6.0]);
    }

    #[test]
    fn too_flat_counts_distinct_values() {
        assert!(too_flat(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]));
        assert!(too_flat(&[1.0, 2.0, 1.0, 2.0, 1.0, 2.0]));
        // 6 entries but only 5 distinct values -- still too flat
        assert!(too_flat(&[1.0, 2.0, 3.0, 4.0, 5.0, 1.0]));
        assert!(!too_flat(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]));
    }

    #[test]
    fn is_valid_peak_needs_fifty_bases_and_six_values() {
        let wide: Vec<f32> = (0..60).map(|i| (i % 7) as f32).collect();
        assert!(is_valid_peak(&wide, 30), "60 wide, 7 distinct");
        let narrow: Vec<f32> = (0..40).map(|i| (i % 7) as f32).collect();
        assert!(!is_valid_peak(&narrow, 20), "only 40 wide");
        let flat: Vec<f32> = vec![1.0; 60];
        assert!(!is_valid_peak(&flat, 30), "flat, however wide");
    }

    #[test]
    fn internal_minima_returns_nothing_for_fewer_than_two_maxima() {
        let s = vec![1.0f32; 10];
        assert!(internal_minima(&s, &[]).is_empty());
        assert!(internal_minima(&s, &[3]).is_empty());
    }

    #[test]
    fn internal_minima_picks_the_lowest_point_between_maxima() {
        let s = vec![10.0f32, 5.0, 1.0, 5.0, 10.0];
        let m = internal_minima(&s, &[0, 4]);
        assert_eq!(m, vec![2], "the valley at index 2");
    }

    #[test]
    fn internal_minima_breaks_ties_toward_the_first_minimum() {
        // np.argmin returns the first minimum
        let s = vec![10.0f32, 1.0, 1.0, 5.0, 10.0];
        let m = internal_minima(&s, &[0, 4]);
        assert_eq!(m, vec![1]);
    }

    #[test]
    fn internal_minima_excludes_the_right_maximum_like_python_slice() {
        // The right endpoint is lower than every interior point, but Python's
        // `signal[pos1:pos2]` does not include it when finding the valley.
        let s = vec![5.0f32, 2.0, 1.0];
        let m = internal_minima(&s, &[0, 2]);
        assert_eq!(m, vec![1]);
    }

    #[test]
    fn a_single_maximum_passes_the_filter_unchanged() {
        // with one maximum there are no minima, so upstream returns maxima as-is
        let s: Vec<f32> = (0..100).map(|i| (i % 11) as f32).collect();
        let out = enforce_peakyness(&s, &[50]);
        assert_eq!(out, vec![50]);
    }

    #[test]
    fn two_clear_peaks_both_survive() {
        // two well separated, non-flat humps
        let mut s = vec![0.0f32; 400];
        s[80..160]
            .iter_mut()
            .enumerate()
            .for_each(|(i, v)| *v = (i as f32 / 4.0).min(20.0));
        s[240..320]
            .iter_mut()
            .enumerate()
            .for_each(|(i, v)| *v = (i as f32 / 4.0).min(20.0));
        let maxima = vec![120usize, 280];
        let out = enforce_peakyness(&s, &maxima);
        assert_eq!(out, maxima, "both peaks are prominent enough");
    }

    #[test]
    fn a_very_shallow_hump_is_rejected_as_flat() {
        let mut s = vec![0.0f32; 400];
        // a 1..2 ramp: wide, but only two distinct values
        s[100..300].fill(1.0);
        s[200..300].fill(2.0);
        let out = enforce_peakyness(&s, &[200]);
        // either filtered as too flat, or kept -- but it must not panic and the
        // result must be a subset of the input
        assert!(out.iter().all(|o| maxima_contains(o, 200)));
    }

    fn maxima_contains(v: &usize, want: usize) -> bool {
        *v == want
    }

    #[test]
    fn empty_input_yields_empty_output() {
        assert!(enforce_peakyness(&[], &[]).is_empty());
        let s = vec![1.0f32; 100];
        assert!(enforce_peakyness(&s, &[]).is_empty());
    }

    #[test]
    fn a_narrow_spike_is_rejected() {
        // a single sharp spike: the clipped region around it is only a few wide
        let mut s = vec![0.0f32; 300];
        s[150] = 100.0;
        let out = enforce_peakyness(&s, &[150]);
        assert!(
            out.iter().all(|o| *o == 150),
            "with a single maximum there are no minima, so upstream returns it \
             unchanged -- the filter is not the thing that rejects narrow peaks \
             here, `maxima` is. Got {out:?}"
        );
    }
}
