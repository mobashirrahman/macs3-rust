//! Prints the RNG's outputs so `oracle/check_rng.py` can diff them against NumPy.
use macs_stats::NumpyRng;

fn join<T: std::fmt::Display>(v: &[T]) -> String {
    v.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn main() {
    for seed in [0i64, 1, 42, 12345] {
        let mut r = NumpyRng::seeded(seed);
        let words: Vec<u32> = (0..20).map(|_| r.next_uint32()).collect();
        println!("words_{seed}\t{}", join(&words));
    }
    // one 50-element permutation per seed, matching check_rng.py
    for seed in [0i64, 1, 7, 42, 12345] {
        let mut r = NumpyRng::seeded(seed);
        let mut a: Vec<i64> = (0..50).collect();
        r.shuffle(&mut a);
        println!("shuffle_{seed}\t{}", join(&a));
    }
    for seed in [0i64, 1, 42] {
        let mut r = NumpyRng::seeded(seed);
        let mut a: Vec<i64> = (0..100).collect();
        r.shuffle(&mut a);
        println!("shuffle100_{seed}\t{}", join(&a[..10]));
    }
}
