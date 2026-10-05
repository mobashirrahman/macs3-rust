//! Regression: `-m FE` is evaluated in `f32` end to end, in both score paths.
//!
//! `compute_foldenrichment` (`MACS3/Signal/ScoreTrack.py:624-645`) reads
//!
//! ```text
//! v[i] = (p[i] + pseudocount)/(c[i] + pseudocount)
//! ```
//!
//! The generated C makes that look like double arithmetic -- `p[i]` is boxed,
//! `PyNumber_Add`/`PyNumber_TrueDivide` run on boxed doubles, and only the store
//! back into the `float32` score array narrows. It is not: `p` and `c` are
//! `float32` numpy arrays (`add_chromosome`, `ScoreTrack.py:236-247`) and
//! numpy >= 2 returns a **numpy scalar** from `arr[i]`, so
//! `__Pyx__PyNumber_Add_object_float`'s `PyFloat_CheckExact` fast path is
//! skipped and the add lands in `numpy.float32.__add__`. Under NEP 50 a Python
//! float is a "weak" operand, so the sum is computed and returned in `float32`
//! -- and so is the quotient of the two sums.
//!
//! The inputs below are chosen so the two readings differ in the fifth decimal
//! that `%.5f` prints. The expected bytes are upstream's
//! (`macs3 bdgcmp -t t.bdg -c c.bdg -m FE -p 1 -o fe.bdg`); with the `f64`
//! reading the first row is `1.70680`.

use macs_core::Genome;
use macs_score::{bedgraph_value, score_value, ScoreMethod, ScoreTrack2};

/// Upstream's rows for `t.bdg` = `2.0`/`1.0` over `c.bdg` = `0.75767`/`0.94709`.
const UPSTREAM_ROWS: &str = "chr1\t0\t100\t1.70681\nchr1\t100\t200\t1.02717\n";

fn fe(method_rows: &[(f32, f32)]) -> Vec<String> {
    method_rows
        .iter()
        .map(|&(t, c)| {
            bedgraph_value(
                score_value(ScoreMethod::FE, t, c, 1.0)
                    .expect("FE cannot fail")
                    .expect("FE is row-local"),
            )
        })
        .collect()
}

#[test]
fn streaming_fe_keeps_the_pseudocount_add_in_f32() {
    // 3 / 1.75767 in f32 is 1.7068051 and prints 1.70681; the same quotient in
    // f64 is 1.7068050459, which prints 1.70680.
    assert_eq!(
        fe(&[(2.0, 0.75767), (1.0, 0.94709)]),
        ["1.70681", "1.02717"]
    );
    // A pseudocount small enough to vanish into the treatment's f32 grid is the
    // other half of the bug: with `-p 0.1` and treat=79539840.0 the f32 add keeps
    // 79539840.0 while an f64 add keeps 79539840.1, and the two quotients are
    // whole ulps apart (upstream prints 16528447.00000, f64 would print
    // 16528446.00000).
    let v = score_value(ScoreMethod::FE, 79539840.0, 4.7123, 0.1)
        .expect("FE cannot fail")
        .expect("FE is row-local");
    assert_eq!(bedgraph_value(v), "16528447.00000");
}

#[test]
fn materialising_fe_matches_the_same_bytes() {
    let mut genome = Genome::new();
    let chrom = genome.intern_str("chr1");
    let mut st = ScoreTrack2::new(genome, 1.0, 1.0);
    st.add(chrom, 100, 2.0, 0.75767);
    st.add(chrom, 200, 1.0, 0.94709);
    st.set_pseudocount(1.0);
    st.change_score_method(ScoreMethod::FE).expect("FE runs");
    let dir = std::env::temp_dir().join(format!("macs3rs-fe-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("fe.bdg");
    st.write_bedgraph(&out, "FE_Scores", "Scores calculated by FE", 3)
        .expect("written");
    let got = std::fs::read_to_string(&out).expect("read back");
    assert_eq!(got, UPSTREAM_ROWS, "compute_foldenrichment bytes");
    let _ = std::fs::remove_dir_all(&dir);
}
