//! Signals consumed by `oracle/check_sg.py` on a freshly built checkout.

use macs_peaks::sg::{maxima, savitzky_golay_order2_deriv1};

fn emit<T: std::fmt::Display>(case: &str, field: &str, values: &[T]) {
    let values: Vec<String> = values.iter().map(ToString::to_string).collect();
    println!("{case}\t{field}\t{}", values.join(","));
}

fn main() {
    let mut single = vec![0.0f32; 300];
    single[150] = 100.0;
    let mut two = vec![0.0f32; 400];
    two[120] = 50.0;
    two[280] = 40.0;
    let ramp: Vec<f32> = (0..200).map(|i| 2.0 * i as f32 + 1.0).collect();

    for (case, signal) in [("single", single), ("two", two), ("ramp", ramp)] {
        let deriv: Vec<f64> = savitzky_golay_order2_deriv1(&signal, 51)
            .into_iter()
            .map(|v| (v * 1e16).round_ties_even() / 1e16)
            .collect();
        let sign: Vec<i32> = deriv
            .iter()
            .map(|&v| {
                if v > 0.0 {
                    1
                } else if v < 0.0 {
                    -1
                } else {
                    0
                }
            })
            .collect();
        emit(case, "deriv", &deriv);
        emit(case, "sign", &sign);
        emit(case, "maxima", &maxima(&signal, 51));
    }
}
