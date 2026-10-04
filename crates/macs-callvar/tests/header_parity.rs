//! L4 golden: the `callvar` VCF header must match the pinned oracle byte for byte.
//!
//! `tests/data/header.golden` is the first 25 lines of a real upstream run
//! (`macs3 callvar -b peaks.bed -t reads.bam -o out.vcf`), captured with
//! `oracle/ENV.lock`'s pinned MACS3 3.0.5. It pins:
//!
//! * the 21 `##` header lines, transcribed verbatim;
//! * the `Program_Args` reconstruction, including `--top2allele-count` (singular
//!   "allele") for an option spelled `--altallele-count`, and the **double space**
//!   before `--fermi`, which falls out of `tmpcmdstr` carrying a leading space;
//! * the `##contig=` lines, in BAM *header* order rather than sorted order;
//! * the `#CHROM` column line.
//!
//! The date and the version are substituted, so the golden stores the recorded date
//! and the test passes the same values.

use macs_callvar::{
    check_chrom_names, program_args, read_peak_bed, today_compact, vcf_header_with_contigs,
    CallvarOptions,
};
use std::path::{Path, PathBuf};

fn opts() -> CallvarOptions {
    CallvarOptions {
        peak_bed: PathBuf::from("peaks.bed"),
        tfile: PathBuf::from("reads.bam"),
        cfile: None,
        ofile: PathBuf::from("out.vcf"),
        outdir: None,
        gq_cutoff_hetero: 0.0,
        gq_cutoff_homo: 0.0,
        q: 20,
        max_duplicate: 1,
        fermi: "auto".into(),
        fermi_min_overlap: 30,
        top2_alleles_min_ratio: 0.8,
        alt_allele_min_count: 2,
        max_ar: 0.95,
        np: 1,
        verbose: 2,
    }
}

/// The invocation the golden was recorded with, i.e. `sys.argv[1:]`.
const ARGV: &[&str] = &[
    "callvar",
    "-b",
    "peaks.bed",
    "-t",
    "reads.bam",
    "-o",
    "uv/out.vcf",
];

#[test]
fn vcf_header_is_byte_identical_to_upstream() {
    let golden = include_str!("data/header.golden");
    // The recorded run's two contigs, in BAM header order. `tests/fixtures/bam/reads.bam`
    // is a 20 kb / 8 kb two-contig BAM; the lengths are what the golden recorded, so
    // they also pin that `##contig=` uses the BAM's own lengths rather than anything
    // derived from the peak file.
    let rlengths = vec![(b"chr1".to_vec(), 20_000u64), (b"chr2".to_vec(), 8_000)];

    let got = vcf_header_with_contigs(
        "20261002",
        "3.0.5",
        &program_args(
            &ARGV.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
            &opts(),
        ),
        &rlengths,
    );

    if got != golden {
        let g: Vec<&str> = got.lines().collect();
        let w: Vec<&str> = golden.lines().collect();
        for i in 0..g.len().max(w.len()) {
            if g.get(i) != w.get(i) {
                panic!(
                    "header line {} differs\n  ours: {:?}\n  gold: {:?}",
                    i + 1,
                    g.get(i),
                    w.get(i)
                );
            }
        }
        panic!(
            "headers differ in length: ours {} lines, golden {} lines",
            g.len(),
            w.len()
        );
    }
}

#[test]
fn program_args_keeps_the_double_space_before_fermi() {
    let s = program_args(
        &ARGV.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        &opts(),
    );
    assert_eq!(
        s,
        "callvar -b peaks.bed -t reads.bam -o uv/out.vcf -Q 20 -D 1 --max-ar 0.95 \
         --top2alleles-mratio 0.8 --top2allele-count 2 -g 0 -G 0  --fermi auto \
         --fermi-overlap 30"
    );
    // `tmpcmdstr` starts with a space and is joined with another one.
    assert!(
        s.contains("-G 0  --fermi"),
        "expected the double space: {:?}",
        s
    );
}

/// The template's shape: 21 `##` lines, the three substitutions, and a trailing
/// newline that upstream supplies with `+ "\n"` rather than baking into the literal.
#[test]
fn header_template_has_the_expected_shape() {
    let h = macs_callvar::vcf_header("20261002", "3.0.5", "args");
    assert_eq!(
        h.matches('\n').count(),
        22,
        "unexpected template line count"
    );
    assert!(h.starts_with("##fileformat=VCFv4.1\n##fileDate=20261002\n"));
    assert!(h.contains("\n##source=MACS_V3.0.5\n"));
    assert!(h.contains("\n##Program_Args=args\n"));
    // The last template line is the PL FORMAT declaration, and the template carries no
    // trailing newline of its own -- upstream appends one with `+ "\\n"`.
    let last = h.lines().last().expect("non-empty");
    assert!(
        last.starts_with("##FORMAT=<ID=PL,"),
        "last line is {:?}",
        last
    );
    assert!(
        last.ends_with("\">"),
        "template line is not newline-terminated: {:?}",
        last
    );
    // Only the three `%s` placeholders are substituted, so no stray `%` survives.
    assert!(!h.contains("%s"));
    // `VCFHEADER` (not `VCFHEADER_0`) declares DBIC; picking the wrong template
    // silently drops it.
    assert!(h.contains("##INFO=<ID=DBIC,"));
    assert_eq!(h.matches("%s").count(), 0);
}

#[test]
fn date_formatting_is_yyyymmdd() {
    assert_eq!(today_compact(std::time::UNIX_EPOCH), "19700101");
    // 2001-09-09T01:46:40Z, the classic overflow boundary.
    let t = std::time::UNIX_EPOCH + std::time::Duration::from_secs(1_000_000_000);
    assert_eq!(today_compact(t), "20010909");
}

#[test]
fn chromosome_assert_is_asymmetric_and_first_only() {
    let a = vec![b"chr1".to_vec(), b"chr2".to_vec()];
    // treatment's first contag present in control -> ok
    assert!(check_chrom_names(&a, &[b"chrX".to_vec(), b"chr1".to_vec()]).is_ok());
    // control's first contig present in treatment -> also ok (the `or` branch)
    assert!(check_chrom_names(&a, &[b"chr1".to_vec()]).is_ok());
    // neither first contig shared -> rejected
    assert!(check_chrom_names(&a, &[b"chrX".to_vec()]).is_err());
    // Only the FIRST contig is compared: a control whose second contig is chr1 but
    // whose first is not still passes, because upstream indexes [0] on both sides.
    assert!(check_chrom_names(&a, &[b"scaffold".to_vec(), b"chr1".to_vec()]).is_ok());
    assert!(check_chrom_names(&[], &[b"chr1".to_vec()]).is_err());
}

#[test]
fn peak_bed_sorts_and_parses_like_peakio() {
    let dir = std::env::temp_dir();
    let p = dir.join("callvar_peaks.bed");
    std::fs::write(
        &p,
        "chr2\t100\t200\nchr1\t300\t400\nchr1\t100\t200\nchr1\t300\t350\n",
    )
    .expect("write");
    let peaks = read_peak_bed(Path::new(&p)).expect("read");
    assert_eq!(peaks.len(), 4);
    assert_eq!(peaks[0].chrom, b"chr1");
    assert_eq!((peaks[0].start, peaks[0].end), (100, 200));
    assert_eq!((peaks[1].start, peaks[1].end), (300, 350));
    assert_eq!((peaks[2].start, peaks[2].end), (300, 400));
    assert_eq!(peaks[3].chrom, b"chr2");

    // A line with fewer than 3 fields is an error, not a panic.
    let bad = dir.join("callvar_bad.bed");
    std::fs::write(&bad, "chr1\t100\n").expect("write");
    assert!(read_peak_bed(Path::new(&bad)).is_err());
    // A non-integer coordinate is an error too (`int()` raises upstream).
    let bad2 = dir.join("callvar_bad2.bed");
    std::fs::write(&bad2, "chr1\tnot_a_number\t200\n").expect("write");
    assert!(read_peak_bed(Path::new(&bad2)).is_err());
    let _ = std::fs::remove_file(p);
    let _ = std::fs::remove_file(bad);
    let _ = std::fs::remove_file(bad2);
}
