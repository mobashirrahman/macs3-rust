//! Derivative bits recorded with the pinned NumPy and Haswell OpenBLAS kernel.

use macs_peaks::sg::savitzky_golay_order2_deriv1;

#[test]
fn convolution_matches_oracle_at_edges_and_across_dot_product_block_sizes() {
    let signal: Vec<f32> = (0..80).map(|i| ((i * 17) % 29) as f32 / 4.0).collect();
    let mut derivatives = std::collections::HashMap::new();
    for line in include_str!("data/summit_convolution.tsv").lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let window: usize = fields[0].parse().unwrap();
        let index: usize = fields[1].parse().unwrap();
        let expected = u64::from_str_radix(fields[2], 16).unwrap();
        let derivative = derivatives
            .entry(window)
            .or_insert_with(|| savitzky_golay_order2_deriv1(&signal, window));
        assert_eq!(
            derivative[index].to_bits(),
            expected,
            "window {window}, position {index}"
        );
    }
}
