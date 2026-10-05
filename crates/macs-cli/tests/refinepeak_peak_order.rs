//! Regression: `refinepeak` emits peaks in `PeakIO.sort` order, not `(chrom,
//! start, end)` order.
//!
//! `refinepeak_cmd.run` calls `peaks.sort()` (`refinepeak_cmd.py:42`) before
//! handing the peaks to `compute_region_tags_from_peaks`, and that is
//! [`PeakIO.sort`](MACS3/IO/PeakIO.py:357-369):
//!
//! ```text
//! chrs = sorted(list(self.peaks.keys()))
//! for chrom in sorted(chrs):
//!     self.peaks[chrom].sort(key=lambda x: x['start'])
//! ```
//!
//! The per-chromosome key is **`start` alone**, and Python's `list.sort` is
//! stable, so peaks sharing a start keep the order they had in the `--bedfile`.
//! Adding `end` to the key -- which is what this port used to do -- permutes them,
//! and the permutation is visible in the output: `compute_region_tags_from_peaks`
//! carries a per-chromosome forward cursor across the peaks
//! (`FixWidthTrack.py:665-700`), so the peak that is scanned first also gets the
//! tags the cursor has not yet consumed.
//!
//! The fixture is two peaks on `chr1` that both start at 1000, listed
//! `b_second` (end 1100) before `a_first` (end 1050) -- the order the argument
//! says to reverse. The expected bytes are upstream's, from
//! `macs3 refinepeak -b peaks.bed -i reads.bed -f BED --o-prefix r`.

use macs_cli::{parse_flags, Options};

/// Two peaks with the same start, the longer one listed first.
const PEAKS: &str = concat!(
    "chr1\t1000\t1100\tb_second\n",
    "chr1\t1000\t1050\ta_first\n",
);

/// One plus and one minus read, both reaching both peak windows.
const TAGS: &str = concat!(
    "chr1\t1050\t1100\tread1\t0\t+\n",
    "chr1\t1040\t1100\tread2\t0\t-\n",
);

/// Two chromosomes and a same-start tie, to pin the chromosome order too.
/// `sorted(list(self.peaks.keys()))` compares the raw `bytes` names, so `chr1`
/// precedes `chr10` precedes `chr2`.
const MULTI_PEAKS: &str = concat!(
    "chr10\t500\t600\tz1\n",
    "chr2\t100\t200\ta1\n",
    "chr2\t100\t150\ta0\n",
    "chr1\t900\t950\tm1\n",
);

const MULTI_TAGS: &str = concat!(
    "chr10\t520\t600\tr1\t0\t+\n",
    "chr2\t110\t200\tr2\t0\t+\n",
    "chr2\t110\t190\tr3\t0\t-\n",
    "chr1\t910\t950\tr4\t0\t+\n",
    "chr1\t905\t940\tr5\t0\t-\n",
);

fn s(x: &str) -> String {
    x.to_string()
}

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-rp-order-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn refinepeak(dir: &std::path::Path, peaks: &str, tags: &str, tag: &str) -> String {
    let bed = dir.join("peaks.bed");
    std::fs::write(&bed, peaks).unwrap();
    let reads = dir.join("reads.bed");
    std::fs::write(&reads, tags).unwrap();
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        reads.to_str().unwrap().to_string(),
        s("-b"),
        bed.to_str().unwrap().to_string(),
        s("-f"),
        s("BED"),
        s("--outdir"),
        s(out.to_str().unwrap()),
        s("--o-prefix"),
        s(tag),
    ];
    let o: Options = parse_flags("refinepeak", &args).expect("args parse");
    macs_cli::commands::refinepeak::refinepeak(&o).expect("refinepeak runs");
    std::fs::read_to_string(out.join(format!("{tag}_refinepeak.bed"))).expect("written")
}

/// Upstream's rows for the tie fixture. `b_second` comes first even though
/// `a_first` ends earlier: `end` is not part of the sort key.
const EXPECTED: &str = concat!(
    "chr1\t1051\t1052\tb_second_F\t2.00\n",
    "chr1\t1051\t1052\ta_first_F\t2.00",
);

#[test]
fn peaks_sharing_a_start_keep_their_bedfile_order() {
    let dir = tmpdir("tie");
    let got = refinepeak(&dir, PEAKS, TAGS, "r");
    assert_eq!(got, EXPECTED);
    // and the reverse listing really does reverse the output, so the assertion
    // above is not vacuous
    let swapped = PEAKS
        .lines()
        .rev()
        .map(|l| format!("{l}\n"))
        .collect::<String>();
    let got = refinepeak(&dir, &swapped, TAGS, "s");
    assert_eq!(
        got,
        concat!(
            "chr1\t1051\t1052\ta_first_F\t2.00\n",
            "chr1\t1051\t1052\tb_second_F\t2.00",
        )
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn chromosomes_are_walked_in_byte_wise_name_order() {
    let dir = tmpdir("chroms");
    let got = refinepeak(&dir, MULTI_PEAKS, MULTI_TAGS, "m");
    // Upstream's bytes: chr1 < chr10 < chr2, and the a1/a0 tie keeps file order.
    assert_eq!(
        got,
        concat!(
            "chr1\t911\t912\tm1_F\t2.00\n",
            "chr10\t300\t301\tz1_F\t0.00\n",
            "chr2\t111\t112\ta1_F\t2.00\n",
            "chr2\t111\t112\ta0_F\t2.00",
        )
    );
    let _ = std::fs::remove_dir_all(&dir);
}
