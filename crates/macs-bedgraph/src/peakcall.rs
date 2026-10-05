//! Peak calling and reporting on a [`BedGraph`]: `summary`, `call_peaks`,
//! `refine_peaks` and `cutoff_analysis`.
//!
//! These are the four `bedGraphTrackI` methods `hmmratac` drives (and which
//! `bdgpeakcall` / `bdgbroadcall` need in their simpler forms). They are kept
//! together because they share the one convention the whole module rests on: a
//! bedGraph stores run **ends**, and the value on `[end[i-1], end[i])` is
//! `value[i]`, so an above-cutoff chunk recovers its start as
//! `pos[above_cutoff - 1]` and the very first chunk starts at 0.
//!
//! # `summary`'s standard deviation is computed over a broken `pre_p`
//!
//! `BedGraph.py:374-399` resets `pre_p = 0` at the top of the **first** loop
//! (per chromosome) but **not** in the second, so the variance pass measures
//! every run of the first chromosome against that chromosome's own ends and then
//! keeps the *last chromosome's* final `pre_p` for all of the rest. The mean is
//! therefore computed over correct run lengths while the variance is not, and
//! `std_v` is not a real standard deviation. `hmmratac` only reads `mean_v`, so
//! this does not change its output, but it is reproduced rather than corrected so
//! that a future caller sees the same number upstream prints.
//!
//! # `call_peaks` measures length as `last.end - first.start`
//!
//! No `+ 1`, so a peak exactly `min_length` long is kept. The summit is the
//! **midpoint of a chunk**, `(tend + tstart) / 2` with C integer division, and
//! ties are resolved by taking the lower median -- not the first maximum.
//!
//! # `refine_peaks` clips a peak to the bedGraph's own breakpoints
//!
//! For every `(start, end, value)` it collects, it sorts `[p, peak_s, peak_e,
//! pre_p]` and keeps the middle two. So the output peak is the intersection of
//! the input peak with the bedGraph runs, and a peak that only partially overlaps
//! the data is trimmed rather than extended. `min_length` is 0, so nothing is
//! dropped for being short.
//!
//! # `cutoff_analysis` silently reports nothing when the first run is above the
//! cutoff
//!
//! `above_cutoff_startpos = pos_array[above_cutoff - 1]` (`BedGraph.py:1326`).
//! When the first run clears the cutoff, `above_cutoff[0] == 0` and NumPy's `-1`
//! wraps, so the first chunk's start becomes the **last** end in the track. The
//! run's measured length is then `last_end_of_chunk - last_end_of_track`, i.e.
//! negative, and `peak_length >= min_length` never holds. Every other chunk is
//! appended with `ts - lastp <= max_gap`, so the peak never closes early either:
//! the whole track collapses into one chunk whose length is negative and is
//! discarded. Reproduced exactly -- it is load-bearing, since a fold-change
//! track whose leading run is above the cutoff yields a header-only report.
//!
//! # `cutoff_analysis` rounds cutoffs to 3 decimals and drops the last one
//!
//! The ladder is `[round(v, 3) for v in np.arange(minv, maxv, s)]` -- note
//! `arange`'s exclusive upper bound, so `maxv` itself is never a cutoff, and
//! `round(0.6000000000000001, 3)` is `0.6` while `round(0.6000000000000001, 5)`
//! is the same, so the 3-decimal step is what matters. Rows are emitted from the
//! **highest** cutoff down, and rows with `npeaks == 0` are skipped, so the file
//! can be shorter than the ladder.

use std::fmt::Write as _;

use macs_core::genome::ChromId;
use macs_core::Coord;

use crate::BedGraph;

/// The six numbers `bedGraphTrackI.summary` returns, in its order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Summary {
    /// Sum of `value * run_length`.
    pub sum: f32,
    /// Total run length covered.
    pub length: i64,
    /// Largest single value.
    pub max: f32,
    /// Smallest single value.
    pub min: f32,
    /// `sum / length`.
    pub mean: f32,
    /// Upstream's `sqrt(variance)`; see the module note on `pre_p`.
    pub std: f32,
}

/// `bedGraphTrackI.summary` (`BedGraph.py:355`).
pub fn summary(bg: &BedGraph) -> Summary {
    let mut sum_v = 0.0f32;
    let mut n_v: i64 = 0;
    // upstream seeds these with +-100000 rather than the real extrema, and
    // `max(v)`/`min(v)` on an empty chromosome would raise, so an empty track
    // keeps the sentinels.
    let mut max_v = -100_000.0f32;
    let mut min_v = 100_000.0f32;
    let mut last_pre_p: Coord = 0;
    for (_, t) in bg.iter_sorted() {
        let mut pre_p: Coord = 0;
        let runs = t.runs();
        for r in runs {
            let ln = r.end.saturating_sub(pre_p) as f32;
            sum_v += r.value * ln;
            n_v += ln as i64;
            pre_p = r.end;
        }
        if !runs.is_empty() {
            max_v = max_v.max(
                runs.iter()
                    .map(|r| r.value)
                    .fold(f32::NEG_INFINITY, f32::max),
            );
            min_v = min_v.min(runs.iter().map(|r| r.value).fold(f32::INFINITY, f32::min));
        }
        last_pre_p = pre_p;
    }
    let mean_v = if n_v == 0 {
        f32::NAN
    } else {
        sum_v / n_v as f32
    };
    let mut variance = 0.0f32;
    // NOTE: `pre_p` is *not* reset per chromosome here, matching upstream.
    let mut pre_p = last_pre_p;
    for (_, t) in bg.iter_sorted() {
        for r in t.runs() {
            let tmp = r.value - mean_v;
            let ln = r.end.saturating_sub(pre_p) as f32;
            variance += tmp * tmp * ln;
            pre_p = r.end;
        }
    }
    variance /= (n_v - 1) as f32;
    Summary {
        sum: sum_v,
        length: n_v,
        max: max_v,
        min: min_v,
        mean: mean_v,
        std: variance.sqrt(),
    }
}

/// A peak reported by [`call_peaks`] / [`refine_peaks`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BdgPeak {
    /// Left coordinate, inclusive.
    pub start: Coord,
    /// Right coordinate, exclusive.
    pub end: Coord,
    /// Chunk midpoint of the winning value.
    pub summit: Coord,
    /// The winning value.
    pub score: f32,
}

/// `__close_peak` (`BedGraph.py:490`): turn collected chunks into a peak if the
/// span clears `min_length`.
///
/// The tie-break is the load-bearing part. Every chunk whose value equals the
/// running maximum becomes a candidate, and the reported summit is
/// `tsummit[(len + 1) // 2 - 1]` -- the **lower median**. Chunks are in position
/// order, so among tied maxima the earliest-positioned one whose index is at or
/// below the median wins. That is neither "first maximum" nor "middle maximum".
fn close_peak(content: &[(Coord, Coord, f32)], min_length: Coord) -> Option<BdgPeak> {
    let peak_length = content[content.len() - 1].1.saturating_sub(content[0].0);
    if peak_length < min_length {
        return None;
    }
    let mut tsummit: Vec<Coord> = Vec::new();
    let mut summit_value = 0.0f32;
    for &(tstart, tend, tvalue) in content {
        if summit_value == 0.0 || summit_value < tvalue {
            tsummit = vec![(tend + tstart) / 2];
            summit_value = tvalue;
        } else if summit_value == tvalue {
            tsummit.push((tend + tstart) / 2);
        }
    }
    // `cython.cast(int, (len+1)/2) - 1` -- integer division, so for an even count
    // this is one below the true midpoint.
    let idx = tsummit.len().div_ceil(2).saturating_sub(1);
    Some(BdgPeak {
        start: content[0].0,
        end: content[content.len() - 1].1,
        summit: tsummit[idx],
        score: summit_value,
    })
}

/// `bedGraphTrackI.call_peaks` (`BedGraph.py:405`): runs at or above `cutoff`,
/// merged across gaps of at most `max_gap`.
///
/// `call_summits` is accepted and ignored, exactly as upstream documents ("summits
/// are always computed").
pub fn call_peaks(
    bg: &BedGraph,
    cutoff: f32,
    min_length: Coord,
    max_gap: Coord,
) -> Vec<(ChromId, Vec<BdgPeak>)> {
    let mut out = Vec::new();
    for (chrom, t) in bg.iter_sorted() {
        let mut peaks: Vec<BdgPeak> = Vec::new();
        let mut content: Vec<(Coord, Coord, f32)> = Vec::new();
        let mut pre_p: Coord = 0;
        for r in t.runs() {
            let (p, v) = (r.end, r.value);
            if v < cutoff {
                pre_p = p;
                continue;
            }
            if let Some(last) = content.last() {
                if pre_p.saturating_sub(last.1) <= max_gap {
                    content.push((pre_p, p, v));
                } else {
                    if let Some(pk) = close_peak(&content, min_length) {
                        peaks.push(pk);
                    }
                    content = vec![(pre_p, p, v)];
                }
            } else {
                content = vec![(pre_p, p, v)];
            }
            pre_p = p;
        }
        if !content.is_empty() {
            if let Some(pk) = close_peak(&content, min_length) {
                peaks.push(pk);
            }
        }
        out.push((chrom, peaks));
    }
    out
}

/// `bedGraphTrackI.refine_peaks` (`BedGraph.py:669`): trim peaks to the
/// bedGraph's own breakpoints and re-derive the summit.
///
/// Chromosomes are the **intersection** of the bedGraph's and the peaks', sorted
/// by name; a peak on a chromosome with no data is silently dropped.
pub fn refine_peaks(bg: &BedGraph, peaks: &[(ChromId, Coord, Coord)]) -> Vec<(ChromId, BdgPeak)> {
    let mut out: Vec<(ChromId, BdgPeak)> = Vec::new();
    for (chrom, t) in bg.iter_sorted() {
        let mut mine: Vec<(Coord, Coord)> = peaks
            .iter()
            .filter(|(c, _, _)| *c == chrom)
            .map(|&(_, s, e)| (s, e))
            .collect();
        mine.sort();
        let runs = t.runs();
        if runs.is_empty() || mine.is_empty() {
            continue;
        }
        let mut content: Vec<(Coord, Coord, f32)> = Vec::new();
        let mut ri = 0usize;
        let mut pi = 0usize;
        let mut pre_p: Coord = 0;
        let (mut p, mut v) = (runs[0].end, runs[0].value);
        let (mut peak_s, mut peak_e) = (mine[0].0, mine[0].1);
        loop {
            if p > peak_s && peak_e > pre_p {
                // sort([p, peak_s, peak_e, pre_p])[1:3]
                let mut xs = [p, peak_s, peak_e, pre_p];
                xs.sort_unstable();
                content.push((xs[1], xs[2], v));
                pre_p = p;
                ri += 1;
                if ri >= runs.len() {
                    break;
                }
                p = runs[ri].end;
                v = runs[ri].value;
            } else if pre_p >= peak_e {
                if let Some(pk) = close_peak(&content, 0) {
                    out.push((chrom, pk));
                }
                content.clear();
                pi += 1;
                if pi >= mine.len() {
                    break;
                }
                peak_s = mine[pi].0;
                peak_e = mine[pi].1;
            } else if peak_s >= p {
                pre_p = p;
                ri += 1;
                if ri >= runs.len() {
                    break;
                }
                p = runs[ri].end;
                v = runs[ri].value;
            } else {
                // `raise Exception("no way here!")` upstream: unreachable, since
                // the three guards above cover every ordering.
                break;
            }
        }
        if !content.is_empty() {
            if let Some(pk) = close_peak(&content, 0) {
                out.push((chrom, pk));
            }
        }
    }
    out
}

/// `bedGraphTrackI.cutoff_analysis` (`BedGraph.py:1262`): the
/// `score npeaks lpeaks avelpeak` report.
///
/// `max_score` clamps the top of the sweep and `min_score` the bottom; both
/// default to the track's own extremes. `steps` is the *divisor*, not a count,
/// so the ladder has roughly `steps` entries.
pub fn cutoff_analysis(
    bg: &BedGraph,
    max_gap: Coord,
    min_length: Coord,
    steps: i64,
    min_score: f32,
    max_score: f32,
) -> String {
    // `minv`, `maxv` and `s` are all `cython.float`, i.e. C float:
    //   * `maxv - minv` is evaluated in f32 and only then widened to double for
    //     the division by `steps` (`BedGraph.c:28941`), and
    //   * the quotient is *stored back into the C float `s`*, so the step
    //     `np.arange` is handed is the f32 rounding of `f64(span) / steps`.
    // `np.arange` itself widens all three arguments to f64.
    let minv = min_score.max(bg_summary_min(bg));
    let maxv = max_score.min(bg_summary_max(bg));
    let s = (f64::from(maxv - minv) / steps as f64) as f32;
    let mut cutoffs: Vec<f32> = Vec::new();
    if s > 0.0 {
        // `np.arange(minv, maxv, s)` -- half-open at the top, and its length is
        // `ceil((maxv - minv) / s)` computed in f64 from the *widened* arguments,
        // which is one short of `maxv` being reached.
        let step = f64::from(s);
        let n = ((f64::from(maxv) - f64::from(minv)) / step).ceil().max(0.0) as usize;
        for i in 0..n {
            // `cutoff_list` holds f64 entries; `cutoff = cutoff_list[n]` assigns
            // into a `cython.float` (`BedGraph.py:1281`), so every ladder entry is
            // narrowed to f32 *before* it is compared or printed (`BedGraph.c:29325`
            // and `:30006`). The narrowing is load-bearing: `round(v, 3)` of a
            // ladder step lands on `0.995`/`0.975`/`0.325`..., whose nearest f64
            // sits just *below* the decimal midpoint while their f32 image sits
            // just above, so `%.2f` prints `1.00`/`0.98`/`0.32` and not
            // `0.99`/`0.97`/`0.33`.
            cutoffs.push(round3(f64::from(minv) + i as f64 * step) as f32);
        }
    }
    let mut npeaks = vec![0i64; cutoffs.len()];
    let mut lpeaks = vec![0i64; cutoffs.len()];
    for (_, t) in bg.iter_sorted() {
        let runs = t.runs();
        for (n, &cutoff) in cutoffs.iter().enumerate() {
            // runs whose value is strictly above the cutoff; the start is the
            // previous run's end (0 for the first)
            let mut total_l: i64 = 0;
            let mut total_p: i64 = 0;
            let mut content: Vec<(Coord, Coord)> = Vec::new();
            let mut lastp: Coord = 0;
            for (i, r) in runs.iter().enumerate() {
                // `score_array > cutoff`: the comparison widens the f32 score to
                // f64 against the f32 cutoff, also widened to f64.
                if f64::from(r.value) <= f64::from(cutoff) {
                    continue;
                }
                let te = r.end;
                // `above_cutoff_startpos = pos_array[above_cutoff-1]`. When the
                // *first* run is above the cutoff, `above_cutoff[0]` is 0, so
                // `pos_array[-1]` wraps to the **last** end and the first chunk
                // reads `(last_end, first_end)`. That makes the run's measured
                // length negative, so a track whose leading run is above the
                // cutoff reports **no** peak at any threshold. Reproduced; see
                // the module note.
                let ts = runs[(i + runs.len() - 1) % runs.len()].end;
                if !content.is_empty() {
                    // the gap is measured from the previous run's *end*
                    if ts.saturating_sub(lastp) <= max_gap {
                        content.push((ts, te));
                    } else {
                        let plen = content[content.len() - 1].1.saturating_sub(content[0].0);
                        if plen >= min_length {
                            total_l += plen as i64;
                            total_p += 1;
                        }
                        content = vec![(ts, te)];
                    }
                } else {
                    content = vec![(ts, te)];
                }
                lastp = te;
            }
            if !content.is_empty() {
                let plen = content[content.len() - 1].1.saturating_sub(content[0].0);
                if plen >= min_length {
                    total_l += plen as i64;
                    total_p += 1;
                }
            }
            lpeaks[n] += total_l;
            npeaks[n] += total_p;
        }
    }
    let mut ret = String::from("score\tnpeaks\tlpeaks\tavelpeak\n");
    for n in (0..cutoffs.len()).rev() {
        if npeaks[n] > 0 {
            let _ = writeln!(
                ret,
                "{:.2}\t{}\t{}\t{:.2}",
                cutoffs[n],
                npeaks[n],
                lpeaks[n],
                lpeaks[n] as f64 / npeaks[n] as f64
            );
        }
    }
    ret
}

fn bg_summary_min(bg: &BedGraph) -> f32 {
    let mut m = f32::INFINITY;
    for (_, t) in bg.iter_sorted() {
        for r in t.runs() {
            m = m.min(r.value);
        }
    }
    if m.is_finite() {
        m
    } else {
        0.0
    }
}

fn bg_summary_max(bg: &BedGraph) -> f32 {
    let mut m = f32::NEG_INFINITY;
    for (_, t) in bg.iter_sorted() {
        for r in t.runs() {
            m = m.max(r.value);
        }
    }
    if m.is_finite() {
        m
    } else {
        0.0
    }
}

/// Python's `round(x, 3)`: half-to-even on the decimal expansion.
fn round3(x: f64) -> f64 {
    macs_stats::py_round(x, 3)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bg_of(name: &[u8], runs: &[(Coord, f32)]) -> BedGraph {
        let mut bg = BedGraph::new(0.0);
        let mut pre = 0;
        for &(e, v) in runs {
            bg.add_loc(name, pre, e, v);
            pre = e;
        }
        bg
    }

    #[test]
    fn call_peaks_keeps_a_peak_exactly_min_length() {
        // no +1: a 100bp peak with min_length 100 is kept
        let bg = bg_of(b"chr1", &[(50, 1.0), (100, 0.0), (200, 5.0)]);
        let r = call_peaks(&bg, 1.0, 100, 0);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].1[0].start, 100);
        assert_eq!(r[0].1[0].end, 200);
    }

    #[test]
    fn call_peaks_uses_strictly_below_to_reset() {
        // a value exactly at the cutoff stays in the peak (`v >= cutoff`)
        let bg = bg_of(b"chr1", &[(100, 2.0), (200, 2.0)]);
        let r = call_peaks(&bg, 2.0, 0, 0);
        assert_eq!(r[0].1.len(), 1);
        assert_eq!(r[0].1[0].end, 200);
    }

    #[test]
    fn summit_is_a_chunk_midpoint_and_ties_take_the_lower_median() {
        // three chunks all at 5.0 -> three tied candidates, len 3 -> index 1
        let bg = bg_of(b"chr1", &[(10, 5.0), (20, 5.0), (30, 5.0)]);
        let r = call_peaks(&bg, 1.0, 0, 100);
        // midpoints: (0+10)/2=5, (10+20)/2=15, (20+30)/2=25
        // index (3+1)/2-1 = 1 -> 15
        assert_eq!(r[0].1[0].summit, 15);
    }

    #[test]
    fn max_gap_merges_neighbouring_chunks() {
        let bg = bg_of(b"chr1", &[(10, 1.0), (60, 0.0), (70, 1.0)]);
        // gap between chunk ends: 60-10 = 50
        let merged = call_peaks(&bg, 1.0, 0, 50);
        assert_eq!(merged[0].1.len(), 1);
        assert_eq!(merged[0].1[0].end, 70);
        let split = call_peaks(&bg, 1.0, 0, 49);
        assert_eq!(split[0].1.len(), 2);
    }

    #[test]
    fn summary_mean_is_the_run_length_weighted_mean() {
        let bg = bg_of(b"chr1", &[(100, 2.0), (200, 4.0)]);
        let s = summary(&bg);
        assert_eq!(s.length, 200);
        assert!((s.mean - 3.0).abs() < 1e-6, "{:?}", s);
    }

    #[test]
    fn refine_peaks_clips_to_the_bedgraph_breakpoints() {
        let bg = bg_of(b"chr1", &[(100, 1.0), (200, 9.0)]);
        // the input peak runs past the last run end
        let got = refine_peaks(&bg, &[(bg.genome().get(b"chr1").unwrap(), 50, 250)]);
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].1.start, 50);
        assert_eq!(got[0].1.end, 200, "clipped at the last breakpoint");
    }

    #[test]
    fn refine_peaks_drops_peaks_on_absent_chromosomes() {
        let bg = bg_of(b"chr1", &[(100, 1.0)]);
        let mut g = bg.genome().clone();
        let other = g.intern(b"chrZ");
        assert!(refine_peaks(&bg, &[(other, 10, 20)]).is_empty());
    }

    #[test]
    fn cutoff_analysis_has_a_header_and_skips_empty_rows() {
        let bg = bg_of(b"chr1", &[(100, 1.0), (200, 10.0)]);
        let t = cutoff_analysis(&bg, 10, 10, 10, 0.0, 20.0);
        assert!(t.starts_with("score\tnpeaks\tlpeaks\tavelpeak\n"));
        for line in t.lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 4, "{line}");
            assert!(f[1].parse::<i64>().unwrap() > 0, "empty rows are skipped");
        }
    }

    #[test]
    fn cutoff_analysis_on_a_constant_track_produces_no_rows() {
        let bg = bg_of(b"chr1", &[(100, 1.0), (200, 1.0)]);
        let t = cutoff_analysis(&bg, 10, 10, 10, 0.0, 10.0);
        assert_eq!(t.lines().count(), 1, "header only: {t:?}");
    }
}
