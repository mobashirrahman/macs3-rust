//! Diagnostic: per-unitig geometry and read assignment for one peak.
//!
//! `--fermi auto` is gated off because the differential stands at 10/16 records. The
//! divergence is localised to treatment read assignment (F222): at the three failing
//! positions the *control* depth matches upstream exactly while the treatment depth is
//! short. This example prints the intermediate that decides that, because the VCF alone
//! cannot show it:
//!
//! ```text
//! cargo run --release -p macs-callvar --example unitig_stats -- \
//!   <treat.bam> <ctrl.bam> <chrom-start> <chrom-end>
//! ```
//!
//! It prints, per unitig, the `[lpos, rpos)` span, the alignment length and how many
//! treatment and control reads were assigned to it, then the totals. Reads assigned to a
//! unitig that does not cover a position of interest never reach `PosReadsInfo`, so the
//! per-unitig numbers say whether a depth shortfall is an assignment problem or a
//! per-position filtering problem.
//!
//! On upstream's own `callvar_testing` peak 7 it reports:
//!
//! ```text
//! reads after outlier removal: T=147 C=17
//! unitigs: 7
//!   lpos=17652311 rpos=17652552 span=241 readsT=4  readsC=3
//!   lpos=17652312 rpos=17652553 span=241 readsT=6  readsC=0
//!   lpos=17652391 rpos=17652553 span=162 readsT=7  readsC=0
//!   lpos=17652416 rpos=17653027 span=611 readsT=87 readsC=7
//!   ...
//! ```
//!
//! Three unitigs overlap position 17652406 and hold 4 + 6 + 7 = 17 candidate treatment
//! reads between them, of which 12 survive the per-position `ra_pos` range filter -- the
//! 15 upstream reports. So the reads exist and are assigned; the question the next
//! attempt has to answer is which unitig each read lands on, since `remap` takes the
//! **first** containing unitig and several overlap.
use macs_callvar::unitig::build_unitig_collection;
use macs_callvar::RACollection;
use macs_io::bam::BamAccessor;
use std::path::Path;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut tbam = BamAccessor::open(Path::new(&a[1])).unwrap();
    let mut cbam = BamAccessor::open(Path::new(&a[2])).unwrap();
    let (chrom, start, end): (&[u8], i64, i64) =
        (b"chr22", a[3].parse().unwrap(), a[4].parse().unwrap());
    let t = tbam
        .reads_in_region(chrom, start as u32, end as u32, 1)
        .unwrap();
    let c = cbam
        .reads_in_region(chrom, start as u32, end as u32, 1)
        .unwrap();
    let mut coll = RACollection::new(chrom, start, end, t, c).unwrap();
    coll.remove_outliers_in_place(5);
    println!(
        "reads after outlier removal: T={} C={}",
        coll.treatment_reads().len(),
        coll.control_reads().len()
    );
    let built = build_unitig_collection(
        chrom,
        coll.left,
        coll.right,
        coll.read_start,
        coll.read_end,
        &coll.peak_refseq_ext,
        coll.treatment_reads(),
        coll.control_reads(),
        30,
    )
    .unwrap();
    let coll_v = match built {
        macs_callvar::Built::Collection(c) => *c,
        _ => {
            println!("no collection");
            return;
        }
    };
    println!("unitigs: {}", coll_v.ura_list.len());
    for u in &coll_v.ura_list {
        println!(
            "  lpos={} rpos={} span={} aln_len={} readsT={} readsC={}",
            u.lpos,
            u.rpos,
            u.rpos - u.lpos,
            u.aln_length(),
            u.reads[0].len(),
            u.reads[1].len()
        );
    }
    let total_t: usize = coll_v.ura_list.iter().map(|u| u.reads[0].len()).sum();
    let total_c: usize = coll_v.ura_list.iter().map(|u| u.reads[1].len()).sum();
    println!("reads assigned to some unitig: T={total_t} C={total_c}");
}
