//! Prints `--cutoff-analysis` peak counts for synthetic signals, so
//! `oracle/check_precomputes.py` can diff them against a transcription.
use macs_peaks::{chromosome_cutoff_stats, cutoff_ladder};

fn main() {
    // (name, pos, score) triples exercising the merge/split paths
    let cases: Vec<(&str, Vec<u32>, Vec<f32>)> = vec![
        ("single_run", vec![0, 500], vec![0.0, 5.0]),
        ("two_merged", vec![0, 200, 400], vec![0.0, 5.0, 5.0]),
        (
            "two_split",
            vec![0, 200, 900, 1100],
            vec![0.0, 5.0, 0.0, 5.0],
        ),
        ("below_min", vec![0, 100], vec![0.0, 5.0]),
        ("exactly_min", vec![0, 200], vec![0.0, 5.0]),
        ("empty", vec![], vec![]),
        ("flat_at_03", vec![0, 1000, 2000], vec![0.3, 0.3, 0.3]),
        ("flat_at_99", vec![0, 1000, 2000], vec![9.9, 9.9, 9.9]),
        (
            "all_above",
            vec![0, 500, 1000, 1500],
            vec![5.0, 6.0, 7.0, 8.0],
        ),
        ("mismatch", vec![0, 200], vec![5.0]),
    ];

    let ladder = cutoff_ladder();
    let mut first = true;
    for c in &ladder {
        if !first {
            print!(",");
        }
        first = false;
        print!("{c}");
    }
    println!();

    for (name, pos, score) in cases {
        let s = chromosome_cutoff_stats(&pos, &score, 50, 200);
        let npeaks: Vec<String> = ladder.iter().map(|c| s.npeaks_at(*c).to_string()).collect();
        let lens: Vec<String> = ladder.iter().map(|c| s.length_at(*c).to_string()).collect();
        println!("{name}\tnpeaks\t{}", npeaks.join(","));
        println!("{name}\tlength\t{}", lens.join(","));
    }
}
