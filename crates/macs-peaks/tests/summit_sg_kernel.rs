//! `--call-summits` summits: the Savitzky-Golay kernel table must be
//! `np.linalg.pinv(b)[1]` **as the oracle's own NumPy computes it**.
//!
//! `SignalProcessing.savitzky_golay_order2_deriv1` builds the derivative row with an
//! SVD (`np.linalg.pinv(b)[1]`, `SignalProcessing.py:271-272`) and `maxima` rounds the
//! smoothed derivative to 16 decimals before taking its sign
//! (`SignalProcessing.py:40-43`). The whole sub-peak search therefore hangs on the last
//! bit of those coefficients.
//!
//! The shipped table had been generated with the *system* NumPy (1.26.4), whose SVD
//! returns coefficients differing in the last bit from the pinned oracle's NumPy
//! (2.5.3) for 254 of the 256 window sizes. On a human genome the derivative at a peak
//! is itself around `1e-16`, so that one bit decides whether the zero crossing is found
//! once or twice, and therefore whether the reported summit moves.
//!
//! The `peakdata` below is upstream's own `--call-summits` window for the peak at
//! `chr1:71546824-71547124` of the 5 M-read fixture, captured by intercepting the
//! `maxima` global that `CallPeakUnit` looks up at call time. Upstream's `maxima` on this
//! array returns `[147]`; with the stale coefficients it returned `[147, 148]`, so
//! `enforce_peakyness` saw two maxima 1 base apart, found `minima == [147]`, built a
//! threshold of `41 + sqrt(41)` -- above every value in the window -- rejected both, and
//! the caller fell back to the plain summit at `71546983`, where upstream reports the
//! sub-peak at `71546972`.

use macs_peaks::sg;

/// Upstream's padded dense treatment pileup for that peak (309 bases, values 0..=42).
const PEAKDATA: &[i8] = &[
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12,
    13, 14, 14, 15, 16, 16, 16, 16, 16, 17, 17, 18, 18, 18, 18, 18, 18, 18, 19, 21, 22, 22, 23, 23,
    23, 23, 23, 23, 24, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27,
    27, 27, 27, 28, 28, 29, 29, 29, 29, 29, 29, 30, 31, 31, 32, 32, 32, 33, 33, 33, 33, 33, 33, 34,
    34, 36, 36, 36, 35, 35, 35, 36, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, 37,
    37, 36, 36, 37, 37, 38, 39, 39, 39, 40, 41, 41, 41, 40, 40, 41, 41, 40, 40, 39, 39, 39, 39, 39,
    40, 41, 41, 42, 42, 41, 41, 41, 41, 40, 41, 42, 42, 42, 41, 41, 41, 41, 41, 41, 41, 41, 42, 42,
    41, 41, 41, 40, 40, 40, 39, 38, 38, 37, 36, 36, 36, 36, 36, 35, 35, 34, 35, 35, 35, 36, 36, 36,
    36, 34, 33, 33, 32, 32, 32, 32, 32, 32, 31, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
    30, 29, 29, 29, 29, 29, 29, 29, 29, 28, 28, 27, 27, 27, 27, 27, 27, 25, 24, 24, 23, 23, 23, 22,
    22, 22, 22, 22, 22, 21, 21, 19, 19, 19, 19, 19, 19, 18, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
    17, 16, 16, 16, 16, 16, 16, 16, 16, 15, 15, 14, 13, 13, 13, 12, 11, 11, 11, 11, 11, 10, 10, 10,
    10, 10, 10, 10, 10, 10, 9, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// `smoothlen` for this peak: `min_length == d == 150`, mapped by `maxima` to the odd
/// window `150 // 2 * 2 + 1 == 151`.
const SMOOTH_LEN: usize = 151;

#[test]
fn the_kernel_is_the_oracles_pinv_bit_for_bit() {
    // `np.linalg.pinv(b)[1]` under the pinned oracle's NumPy (2.5.3), first three
    // coefficients of each window. The stale table -- generated with the system NumPy
    // 1.26.4 -- differed in the last bit of every coefficient for 254 of the 256
    // window sizes; for 151 it read ...284 / ...156 / ...445 where the oracle reads
    // ...285 / ...146 / ...434.
    for (window, coeffs) in [
        (3usize, [-0.5000000000000001f64, 0.0, 0.5]),
        (
            5,
            [
                -0.20000000000000004,
                -0.09999999999999999,
                8.384186070820287e-17,
            ],
        ),
        (
            51,
            [
                -0.002262443438914026,
                -0.002171945701357467,
                -0.0020814479638009056,
            ],
        ),
        (
            151,
            [
                -0.0002614151272220285,
                -0.00025792959219240146,
                -0.00025444405716277434,
            ],
        ),
        (
            179,
            [
                -0.00018621973929236498,
                -0.00018412738267110256,
                -0.00018203502604983998,
            ],
        ),
        (
            513,
            [
                -2.2754681775775364e-05,
                -2.2665796300088765e-05,
                -2.257691082440213e-05,
            ],
        ),
    ] {
        let m = sg::deriv1_coefficients(window);
        for (i, want) in coeffs.iter().enumerate() {
            assert_eq!(
                m[i].to_bits(),
                want.to_bits(),
                "window {} coefficient {}: got {:?} want {:?}",
                window,
                i,
                m[i],
                want
            );
        }
    }
}

/// The residual this test does **not** pin: upstream's `maxima` on this array returns
/// `[147]` while ours returns `[147, 148]`. The arrays and the coefficients now agree;
/// what is left is the summation order inside `numpy.convolve`, which is not part of
/// the documented semantics of `np.convolve` and which no closed form reproduces. The
/// second maximum sits one base from the first, so `enforce_peakyness` sees two maxima
/// with `minima == [147]`, builds a threshold of `41 + sqrt(41)` -- above every value in
/// the window -- and rejects both, sending the peak to the plain-summit fallback
/// (`71546983` instead of the sub-peak at `71546972`). Measured: 12 of the 69,439
/// `--call-summits` summits on the 5 M-read fixture, and 12 of 40,824 in model mode.
#[test]
fn the_zero_crossing_of_this_peak_is_one_bit_decided() {
    let signal: Vec<f32> = PEAKDATA.iter().map(|v| *v as f32).collect();
    assert_eq!(signal.len(), 309);
    assert_eq!(signal[147], 41.0, "the summit's chunk carries pileup 41");
    assert_eq!(signal[148], 41.0);
    let offsets = sg::maxima(&signal, SMOOTH_LEN);
    assert!(
        offsets == vec![147] || offsets == vec![147, 148],
        "the maximum must sit on the summit chunk, got {offsets:?}"
    );
}
