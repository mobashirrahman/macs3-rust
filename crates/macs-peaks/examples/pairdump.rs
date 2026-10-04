//! Prints the treatment/control pairing for the cases `oracle/check_pairing.py`
//! checks, one JSON object per line, prefixed with `case `.
use macs_peaks::pair_treat_ctrl;

/// One pairing case: treatment and control, each a position/value pair.
type Case = (Vec<u64>, Vec<f32>, Vec<u64>, Vec<f32>);

fn main() {
    // (treat_pos, treat_val, ctrl_pos, ctrl_val) -- kept in step with CASES in
    // oracle/check_pairing.py
    let cases: [Case; 6] = [
        (vec![10, 20], vec![1.0, 2.0], vec![10, 20], vec![5.0, 6.0]),
        (vec![10, 15, 20], vec![1.0, 2.0, 3.0], vec![20], vec![9.0]),
        (vec![10, 20, 30], vec![1.0, 2.0, 3.0], vec![10], vec![5.0]),
        (
            vec![10, 20],
            vec![1.0, 2.0],
            vec![5, 15, 25],
            vec![0.5, 1.5, 2.5],
        ),
        (
            vec![10, 20, 30, 40],
            vec![1.0, 2.0, 3.0, 4.0],
            vec![20],
            vec![0.5],
        ),
        (vec![10], vec![1.0], vec![10, 20, 30], vec![1.0, 2.0, 3.0]),
    ];
    for (tp, tv, cp, cv) in cases {
        let p = pair_treat_ctrl(&tp, &tv, &cp, &cv);
        // hand-rolled JSON: no serde dependency in this crate
        let arr = |v: &Vec<u64>| {
            v.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        };
        let arrf = |v: &Vec<f32>| {
            v.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        };
        println!(
            "case {{\"pos\":\"{}\",\"treat\":\"{}\",\"ctrl\":\"{}\"}}",
            arr(&p.pos),
            arrf(&p.treat),
            arrf(&p.ctrl)
        );
    }
}
