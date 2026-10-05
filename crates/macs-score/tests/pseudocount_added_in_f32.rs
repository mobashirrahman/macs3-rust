//! Regression: the p-score pseudocount is added in `float32`, not `float64`.
//!
//! `compute_pvalue` (`MACS3/Signal/ScoreTrack.py:451-453`) reads
//!
//! ```text
//! v[i] = get_pscore(cython.cast(cython.int, (p[i] + self.pseudocount)),
//!                   c[i] + self.pseudocount)
//! ```
//!
//! and the generated C for that line (`ScoreTrack.c:14012-14030`) is
//!
//! ```c
//! __pyx_t_2 = __Pyx_GetItemInt(p, i, ...);            /* np.float32 scalar */
//! __pyx_t_4 = PyFloat_FromDouble(self->pseudocount);   /* python float      */
//! __pyx_t_8 = PyNumber_Add(__pyx_t_2, __pyx_t_4);
//! __pyx_t_9 = __Pyx_PyLong_As_int(__pyx_t_8);          /* truncate          */
//! ```
//!
//! so the sum is `np.float32 + python float`, which under NEP 50 (numpy 2.5.3, the
//! pinned version) keeps the `float32` width -- and only then is it truncated to
//! `int`. Adding the pseudocount in `f64` first is a different computation:
//!
//! ```text
//! treat = 9.9f32  = 9.899999618530273
//! pc    = 0.1f32  = 0.10000000149011612
//! f32 sum = 10.0               -> count 10 -> pscore 18.64095   (upstream)
//! f64 sum = 9.999999620020389  -> count  9 -> pscore 16.59923   (wrong)
//! ```
//!
//! The expected bytes below are upstream's, from
//! `macs3 bdgcmp -t t6.bdg -c c6.bdg -m ppois -p 0.1 --o-prefix r`. The control
//! side (`c[i] + pseudocount`) also retains float32 width.

use macs_core::Genome;
use macs_score::{
    bedgraph_value, pseudocounted_inputs, score_value, ScoreMethod, ScoreTrack2,
    DEFAULT_PSEUDOCOUNT,
};

/// `t6.bdg` / `c6.bdg`: `(treat, ctrl)` per 10 bp interval, with four of the six
/// treatment values chosen so the `f32` and `f64` pseudocount sums truncate
/// differently.
const ROWS: &[(f32, f32)] = &[
    (9.9, 0.0),
    (19.9, 5.0),
    (0.0, 0.0),
    (3.7, 1.1),
    (4.7, 2.0),
    (2.0, 0.25),
];

/// Upstream's `-m ppois -p 0.1` rows for [`ROWS`].
const PPOIS_P01: &str = concat!(
    "chr1\t0\t10\t18.64095\n",
    "chr1\t10\t20\t6.95136\n",
    "chr1\t20\t30\t1.02153\n",
    "chr1\t30\t40\t1.47148\n",
    "chr1\t40\t50\t1.20673\n",
    "chr1\t50\t60\t2.25893\n",
);

/// Upstream's `-m ppois -p 1.0` rows for [`ROWS`].
const PPOIS_P1: &str = concat!(
    "chr1\t0\t10\t7.99793\n",
    "chr1\t10\t20\t5.83711\n",
    "chr1\t20\t30\t0.57800\n",
    "chr1\t30\t40\t1.20673\n",
    "chr1\t40\t50\t1.07615\n",
    "chr1\t50\t60\t1.41715\n",
);

/// Upstream's `-m qpois -p 0.1` rows for [`ROWS`]; `qpois` cannot take the
/// streaming path (`Q` needs the genome-wide p-value table), so this is the
/// materialising `ScoreTrack2::compute_pvalue` + `compute_qvalue`.
const QPOIS_P01: &str = concat!(
    "chr1\t0\t10\t16.86280\n",
    "chr1\t10\t20\t6.21460\n",
    "chr1\t20\t30\t0.00000\n",
    "chr1\t30\t40\t1.18469\n",
    "chr1\t40\t50\t1.04136\n",
    "chr1\t50\t60\t1.80300\n",
);

fn streamed(pseudocount: f32) -> Vec<String> {
    ROWS.iter()
        .map(|&(t, c)| {
            bedgraph_value(
                score_value(ScoreMethod::P, t, c, pseudocount)
                    .expect("lambda is positive on every row")
                    .expect("P is row-local"),
            )
        })
        .collect()
}

fn materialised(pseudocount: f32, method: ScoreMethod, tag: &str) -> String {
    let mut genome = Genome::new();
    let chrom = genome.intern_str("chr1");
    let mut st = ScoreTrack2::new(genome, 1.0, 1.0);
    for (i, &(t, c)) in ROWS.iter().enumerate() {
        st.add(chrom, 10 * (i as u32 + 1), t, c);
    }
    st.set_pseudocount(pseudocount);
    st.change_score_method(method).expect("scoring runs");
    let dir = std::env::temp_dir().join(format!("macs3rs-pc-f32-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let out = dir.join("score.bdg");
    st.write_bedgraph(&out, "P_Scores", "Scores calculated by ppois", 3)
        .expect("written");
    let got = std::fs::read_to_string(&out).expect("read back");
    let _ = std::fs::remove_dir_all(&dir);
    got
}

#[test]
fn the_treatment_sum_is_truncated_after_rounding_to_f32() {
    // 9.899999618530273 + 0.10000000149011612 is exactly 10.0 in f32 and
    // 9.999999620020389 in f64, so the count is 10, not 9.
    assert_eq!(pseudocounted_inputs(9.9, 0.0, 0.1).0, 10);
    assert_eq!(pseudocounted_inputs(19.9, 5.0, 0.1).0, 20);
    // the lambda side stays f32 and is unaffected by the width of the treatment add
    assert_eq!(pseudocounted_inputs(9.9, 0.0, 0.1).1, 0.1f32 + 0.0f32);
    // rows where both readings agree must stay put
    assert_eq!(pseudocounted_inputs(0.0, 0.0, DEFAULT_PSEUDOCOUNT).0, 1);
    assert_eq!(pseudocounted_inputs(3.7, 1.1, 0.1).0, 3);
    assert_eq!(pseudocounted_inputs(4.7, 2.0, 0.3).0, 5);
}

#[test]
fn streaming_ppois_counts_the_f32_sum() {
    let want = |rows: &str| -> Vec<String> {
        rows.lines()
            .map(|l| l.split('\t').nth(3).unwrap().to_string())
            .collect()
    };
    assert_eq!(streamed(0.1), want(PPOIS_P01));
    assert_eq!(streamed(1.0), want(PPOIS_P1));
}

#[test]
fn materialising_ppois_and_qpois_count_the_f32_sum() {
    assert_eq!(materialised(0.1, ScoreMethod::P, "p"), PPOIS_P01);
    assert_eq!(materialised(0.1, ScoreMethod::Q, "q"), QPOIS_P01);
}
