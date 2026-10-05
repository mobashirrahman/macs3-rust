//! Paired-peak discovery against upstream's own CTCF ChIP-seq data.
//!
//! The generated corpus in `tests/fixtures` is too small for the fragment model:
//! every fixture has `model_capable=0`, so `predictd`/`callpeak` model mode was
//! never gated against the oracle on real data. That gap hid a real divergence --
//! on `CTCF_SE_ChIP_chr22_50k.bed.gz` upstream finds **464** paired strand peaks
//! and fits `d = 229`, while this port found 66 and refused to build a model at
//! all. Upstream's own `test/cmdlinetest` asserts the 229.
//!
//! This test drives the same input through `PeakModel` and pins the counts
//! upstream reports under `--verbose 3`
//! (`PeakModel.py:121,177-192`, captured from the pinned oracle):
//!
//! ```text
//! #2 min_tags: 1; max_tags:14;
//! Number of unique tags on + strand: 24867
//! Number of peaks in + strand: 1722
//! Number of unique tags on - strand: 24755
//! Number of peaks in - strand: 1682
//! Number of paired peaks in this chromosome: 464
//! #2 Total number of paired peaks: 464
//! # predicted fragment length is 229 bps
//! ```
//!
//! The oracle tree is not vendored (see `oracle/ENV.lock`), so the test locates
//! it through the environment, skips with a clear message when it is absent, and
//! never silently passes.

use macs_core::{Result, Strand};
use macs_model::{ModelOptions, PeakModel};
use macs_track::SingleEndTrackBuilder;

mod oracle;

fn load_se(path: &std::path::Path) -> Result<SingleEndTrackBuilder> {
    let mut r = macs_io::open_maybe_gzip(path).expect("open");
    let mut line = Vec::new();
    let mut b = SingleEndTrackBuilder::new();
    loop {
        line.clear();
        use std::io::BufRead;
        if r.read_until(b'\n', &mut line).expect("read") == 0 {
            break;
        }
        let Some(rec) = macs_io::parse_bed_line(&line).expect("parse") else {
            continue;
        };
        if rec.pos < 0 || rec.chrom.is_empty() {
            continue;
        }
        b.push(&rec.chrom, rec.pos as u32, rec.strand);
    }
    b.finalize();
    Ok(b)
}

/// Upstream's `min_tags`/`max_tags` for this input: `round(total * fold * 2*bw / gz / 2)`
/// with `total = 49622`, `bw = 300`, `gsize = 52e6`, `--mfold 5 50`.
#[test]
fn paired_peak_counts_match_upstream_on_ctcf() {
    let Some(oracle) = oracle::require("paired_peak_counts_match_upstream_on_ctcf") else {
        return;
    };
    let testdir = oracle.test_dir();
    let chip = testdir.join("CTCF_SE_ChIP_chr22_50k.bed.gz");
    if !chip.exists() {
        eprintln!("skipping: {} not present", chip.display());
        return;
    }

    let builder = load_se(&chip).expect("load CTCF ChIP");
    let track = builder.build();
    assert_eq!(
        track.total(),
        49622,
        "upstream reports 24867 + strand and 24755 - strand tags, total 49622"
    );

    let options = ModelOptions {
        gsize: 52_000_000.0,
        mfold: (5.0, 50.0),
        bw: 300,
        d_min: 10,
    };

    let model = PeakModel::build(&track, options).expect(
        "upstream finds 464 paired peaks on this input and fits d=229; \
         a NotEnoughPairs error here means find_paired_peaks regressed",
    );

    assert_eq!(
        (model.min_tags, model.max_tags),
        (1.0, 14.0),
        "min_tags/max_tags must match upstream's `#2 min_tags: 1; max_tags:14`"
    );
    assert_eq!(
        model.paired.total_pairs, 464,
        "upstream: `#2 Total number of paired peaks: 464`"
    );
    assert_eq!(
        model.d.trunc(),
        229.0,
        "upstream: `# predicted fragment length is %d bps` \
         -- `%d` on the float truncates, so the reported value is d.trunc()"
    );
    assert_eq!(
        model.alternative_d,
        vec![229],
        "upstream: `# alternative fragment length(s) may be 229 bps`"
    );

    // Per-strand summits, from upstream's `--verbose 3` output.
    let chr22 = model
        .paired
        .per_chrom
        .iter()
        .find(|(c, ..)| c.as_slice() == b"chr22")
        .expect("chr22");
    assert_eq!(
        chr22.1, 1722,
        "upstream: `Number of peaks in + strand: 1722`"
    );
    assert_eq!(
        chr22.2, 1682,
        "upstream: `Number of peaks in - strand: 1682`"
    );
    assert_eq!(
        chr22.3, 464,
        "upstream: `Number of paired peaks in this chromosome: 464`"
    );

    // Strand tag counts are a cheap precondition: if the pileup itself is wrong,
    // every downstream count is wrong for the same reason.
    let pos = track.positions();
    let plus: usize = pos
        .chroms_sorted()
        .iter()
        .map(|c| pos.strand(*c, Strand::Plus).len())
        .sum();
    let minus: usize = pos
        .chroms_sorted()
        .iter()
        .map(|c| pos.strand(*c, Strand::Minus).len())
        .sum();
    assert_eq!(plus, 24867, "upstream: `unique tags on + strand: 24867`");
    assert_eq!(minus, 24755, "upstream: `unique tags on - strand: 24755`");
}
