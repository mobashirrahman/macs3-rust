//! Regression: `refinepeak` must carry upstream's per-chromosome tag cursor.
//!
//! `find_summit` (`refinepeak_cmd.py:67-102`) is pure, but the tags it is handed
//! are not: `FWTrack.compute_region_tags_from_peaks`
//! (`FixWidthTrack.py:659-704`) walks **one forward cursor per chromosome and
//! strand**, keeps it across peaks, and rewinds it by less than `window_size`
//! afterwards:
//!
//! ```text
//! for i in range(prev_i, plus.shape[0]):
//!     pos = plus[i]
//!     if pos < startpos:   continue
//!     elif pos > endpos:   prev_i = i; break
//!     else:                temp.append(pos)
//! ...
//! for i in range(prev_i, 0, -1):
//!     if plus[prev_i] - plus[i] >= window_size: break
//! prev_i = i
//! ```
//!
//! Two details of that rewind decide which tags a *later* peak sees, and both
//! have to be reproduced verbatim:
//!
//! 1. `range(prev_i, 0, -1)` never visits index `0`, so a rewind that does not
//!    `break` parks the cursor on index `1` and index `0` of that strand is
//!    dropped for the rest of the chromosome.
//! 2. When `prev_i` was `0` the loop body never runs, so `prev_i = i` reads the
//!    **collection loop's own loop variable** -- the `break` index, or the last
//!    index when it ran off the end. It is not a no-op.
//!
//! Collecting "every tag in `[startpos, endpos]`" per peak instead (what this
//! port used to do) is correct arithmetic over the wrong tag set. On the 5 M-read
//! CTCF fixture it is worth 24 of 36769 summits; peak 1436 below sees 1 plus and
//! 3 minus tags rather than 5 and 7, so upstream reports `2.8284271247461903` --
//! which fails the default `--cutoff 5`, hence the `_F` suffix and the summit at
//! 86968778 -- where the stateless window reports `5.745966692414834` (`_R`) at
//! 86968506.
//!
//! The tags are the 60 BED reads of `chr1:86967000-86969500` from
//! `test/CTCF_12878_5M.bed.gz` that reach the three peaks below; the expected
//! bytes are upstream's
//! (`macs3 refinepeak -b peaks.narrowPeak -i reads.bed --o-prefix r`).

use macs_cli::{parse_flags, Options};

/// Three peaks of that run, verbatim from `callpeak --outdir` output.
const PEAKS: &str = concat!(
    "chr1\t86967794\t86968154\tprep_peak_1434\t305\t.\t15.4076\t33.5298\t30.5438\t171\n",
    "chr1\t86968358\t86968577\tprep_peak_1435\t75\t.\t6.82696\t10.1584\t7.59483\t166\n",
    "chr1\t86968634\t86968834\tprep_peak_1436\t22\t.\t3.79275\t4.67973\t2.27537\t111\n",
);

/// The 60 reads, verbatim from that slice of `CTCF_12878_5M.bed.gz`.
const TAGS: &str = concat!(
    "chr1	86967216	86967252	.	.	+\n",
    "chr1	86967242	86967278	.	.	+\n",
    "chr1	86967340	86967376	.	.	+\n",
    "chr1	86967521	86967557	.	.	+\n",
    "chr1	86967735	86967771	.	.	+\n",
    "chr1	86967872	86967908	.	.	+\n",
    "chr1	86967889	86967925	.	.	+\n",
    "chr1	86967890	86967926	.	.	+\n",
    "chr1	86967901	86967937	.	.	+\n",
    "chr1	86967901	86967937	.	.	+\n",
    "chr1	86967916	86967952	.	.	+\n",
    "chr1	86967919	86967955	.	.	+\n",
    "chr1	86967923	86967959	.	.	+\n",
    "chr1	86967923	86967959	.	.	+\n",
    "chr1	86967928	86967964	.	.	+\n",
    "chr1	86967934	86967970	.	.	+\n",
    "chr1	86967936	86967972	.	.	+\n",
    "chr1	86967938	86967974	.	.	+\n",
    "chr1	86967940	86967976	.	.	+\n",
    "chr1	86967940	86967976	.	.	+\n",
    "chr1	86967942	86967978	.	.	+\n",
    "chr1	86967943	86967979	.	.	+\n",
    "chr1	86967949	86967985	.	.	+\n",
    "chr1	86967951	86967987	.	.	+\n",
    "chr1	86967951	86967987	.	.	+\n",
    "chr1	86967952	86967988	.	.	+\n",
    "chr1	86967952	86967988	.	.	+\n",
    "chr1	86967954	86967990	.	.	+\n",
    "chr1	86967958	86967994	.	.	+\n",
    "chr1	86967963	86967999	.	.	+\n",
    "chr1	86968051	86968087	.	.	+\n",
    "chr1	86968191	86968227	.	.	+\n",
    "chr1	86968197	86968233	.	.	+\n",
    "chr1	86968467	86968503	.	.	+\n",
    "chr1	86968471	86968507	.	.	+\n",
    "chr1	86968505	86968541	.	.	+\n",
    "chr1	86968634	86968670	.	.	+\n",
    "chr1	86968731	86968767	.	.	+\n",
    "chr1	86969214	86969250	.	.	+\n",
    "chr1	86969469	86969505	.	.	+\n",
    "chr1	86967436	86967472	.	.	-\n",
    "chr1	86967548	86967584	.	.	-\n",
    "chr1	86967925	86967961	.	.	-\n",
    "chr1	86967949	86967985	.	.	-\n",
    "chr1	86967955	86967991	.	.	-\n",
    "chr1	86967958	86967994	.	.	-\n",
    "chr1	86967963	86967999	.	.	-\n",
    "chr1	86968063	86968099	.	.	-\n",
    "chr1	86968066	86968102	.	.	-\n",
    "chr1	86968072	86968108	.	.	-\n",
    "chr1	86968081	86968117	.	.	-\n",
    "chr1	86968111	86968147	.	.	-\n",
    "chr1	86968507	86968543	.	.	-\n",
    "chr1	86968511	86968547	.	.	-\n",
    "chr1	86968514	86968550	.	.	-\n",
    "chr1	86968522	86968558	.	.	-\n",
    "chr1	86968541	86968577	.	.	-\n",
    "chr1	86968841	86968877	.	.	-\n",
    "chr1	86968863	86968899	.	.	-\n",
    "chr1	86969043	86969079	.	.	-\n",
);

fn s(x: &str) -> String {
    x.to_string()
}

fn tmpdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("macs3rs-rp-cursor-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn refinepeak(dir: &std::path::Path) -> String {
    let tags = dir.join("reads.bed");
    std::fs::write(&tags, TAGS).unwrap();
    let peaks = dir.join("peaks.narrowPeak");
    std::fs::write(&peaks, PEAKS).unwrap();
    let out = dir.join("o");
    std::fs::create_dir_all(&out).unwrap();
    let args = vec![
        s("-i"),
        tags.to_str().unwrap().to_string(),
        s("-b"),
        peaks.to_str().unwrap().to_string(),
        s("-f"),
        s("BED"),
        s("--outdir"),
        s(out.to_str().unwrap()),
        s("--o-prefix"),
        s("r"),
    ];
    let o: Options = parse_flags("refinepeak", &args).expect("args parse");
    macs_cli::commands::refinepeak::refinepeak(&o).expect("refinepeak runs");
    std::fs::read_to_string(out.join("r_refinepeak.bed")).expect("written")
}

#[test]
fn a_later_peak_sees_only_the_tags_the_cursor_left_it() {
    let dir = tmpdir("cursor");
    let got = refinepeak(&dir);
    // output is newline-*separated* (upstream joins with \n, no trailing newline)
    assert_eq!(
        got,
        concat!(
            "chr1\t86967959\t86967960\tprep_peak_1434_R\t28.98\n",
            "chr1\t86968506\t86968507\tprep_peak_1435_R\t6.75\n",
            "chr1\t86968778\t86968779\tprep_peak_1436_F\t2.83",
        )
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn the_third_peak_is_not_refined_from_the_full_window() {
    let dir = tmpdir("tagset");
    let got = refinepeak(&dir);
    let last = got.lines().last().expect("three rows");
    let f: Vec<&str> = last.split('\t').collect();
    // The stateless window finds the Watson/Crick optimum at 86968506 with
    // 5.745966692414834, which clears the cutoff. Upstream's cursor hands the
    // peak one plus and three minus tags instead, the optimum drops to 2.83 and
    // the peak is reported as *failed* at a different summit.
    assert_eq!(f[1], "86968778", "summit: {last}");
    assert!(f[3].ends_with("_F"), "fails the default cutoff: {last}");
    assert_eq!(f[4], "2.83", "wtd: {last}");
    let _ = std::fs::remove_dir_all(&dir);
}
