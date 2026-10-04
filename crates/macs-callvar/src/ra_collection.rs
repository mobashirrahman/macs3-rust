//! Port of `MACS3/Signal/RACollection.py` and the `callvar` calling loop --
//! everything `callvar` does without local assembly.
//!
//! # The reference sequence is a read consensus, not the genome
//!
//! `callvar` never reads a genome FASTA. `__get_peak_REFSEQ` builds the peak's
//! reference **from the reads themselves**: every gap in the treatment pileup is
//! back-filled with the previous read's reference bases, then the control pileup
//! back-fills over the top. That is what `s[i - left]` indexes in the calling loop,
//! and it is why the whole thing is self-contained: no FASTA, no `--genome`.
//!
//! # `__fill_refseq` writes the previous read's bases into the current read's span
//!
//! ```python
//! for i in range(len(ralist)):
//!     read = ralist[i]
//!     if read["lpos"] > prev_r:
//!         read = ralist[i - 1]                  # <-- previous read
//!         read_refseq = read.get_REFSEQ()       # <-- its bases
//!         ind   = read["lpos"] - start          # <-- current read's lpos
//!         ind_r = ind + read["rpos"] - read["lpos"]
//!         seq[ind: ind_r] = read_refseq
//!         prev_r = read["rpos"]
//! ```
//!
//! The span is sized by the **current** read and filled with the **previous** one's
//! bases, so the two lengths routinely disagree -- and Python's `bytearray` slice
//! assignment *changes the buffer length* when they do, shifting every later
//! coordinate. [`py_splice`] reproduces the resize; a fixed-capacity overwrite would
//! not, and would silently produce different reference bases downstream. This is not
//! hypothetical: it fires on upstream's own `callvar_testing` BAM on the first peak,
//! where a 100-byte slice receives 101 bytes of consensus.
//!
//! # `remove_outliers` is a global percentile across both samples
//!
//! The edit-count percentile is computed over **all** treatment *and* control reads
//! together, then applied to both lists. So a pile of noisy control reads can raise
//! the threshold and evict clean treatment reads. It runs *after* the collection is
//! built, so the consensus still comes from the unfiltered reads.
//!
//! # What is not here
//!
//! `--fermi on|auto` sends peaks carrying an indel or a reference-biased het through
//! fermi-lite assembly, which needs a de novo assembler -- the last item on the
//! porting plan. [`is_assembly_implemented`] reports that, and the driver refuses
//! rather than quietly emitting a VCF with records missing.

use crate::peak_variants::{PeakVariants, Variant};
use crate::pos_reads_info::{DepthOpt, PosReadsInfo};
use macs_io::bam::ReadAlignment;

/// Thresholds, mirroring the CLI defaults in `callvar_cmd.run`.
#[derive(Debug, Clone)]
pub struct CallParams {
    /// `top2allelesMinRatio`, passed to `PosReadsInfo::update_top_alleles`.
    #[allow(clippy::excessive_precision)]
    pub top2alleles_min_ratio: f32,
    /// `altalleleMinCount`.
    pub min_alt_allele_count: i32,
    /// `maxAR`.
    #[allow(clippy::excessive_precision)]
    pub max_allowed_ar: f32,
    /// `GQCutoffHomo`, the `-G` default.
    pub min_homo_gq: i32,
    /// `GQCutoffHetero`, the `-g` default.
    pub min_heter_gq: i32,
    /// `Q`, the base-quality cutoff handed to `add_T`/`add_C`.
    pub min_q: i32,
    /// `maxDuplicate`.
    pub max_duplicate: u32,
}

impl Default for CallParams {
    /// Upstream's argparse defaults, read off `macs3 callvar --help`.
    fn default() -> Self {
        Self {
            top2alleles_min_ratio: 0.8,
            min_alt_allele_count: 2,
            max_allowed_ar: 0.99,
            min_homo_gq: 50,
            min_heter_gq: 100,
            min_q: 20,
            max_duplicate: 1,
        }
    }
}

/// One peak's reads, plus the reference consensus built from them.
#[derive(Debug, Clone)]
pub struct RACollection {
    pub chrom: Vec<u8>,
    /// Peak start.
    pub left: i64,
    /// Peak end.
    pub right: i64,
    /// Leftmost coordinate covered by any read (or the peak, whichever is smaller).
    pub read_start: i64,
    /// Rightmost coordinate covered by any read (or the peak, whichever is larger).
    pub read_end: i64,
    treatment: Vec<ReadAlignment>,
    control: Vec<ReadAlignment>,
    /// `peak_refseq`: the consensus over `[left, right)`.
    pub peak_refseq: Vec<u8>,
    /// `peak_refseq_ext`: the consensus over `[read_start, read_end)`.
    pub peak_refseq_ext: Vec<u8>,
}

/// Upstream raises `Exception("No reads from ChIP sample to construct RAcollection!")`.
pub const NO_READS: &str = "No reads from ChIP sample to construct RAcollection!";

impl RACollection {
    /// `RACollection.__init__`.
    ///
    /// `Err` when the treatment list is empty -- upstream's one hard precondition
    /// here -- or when a read's `MD` cannot be interpreted.
    pub fn new(
        chrom: &[u8],
        peak_start: i64,
        peak_end: i64,
        treatment: Vec<ReadAlignment>,
        control: Vec<ReadAlignment>,
    ) -> Result<Self, String> {
        if treatment.is_empty() {
            return Err(NO_READS.to_string());
        }
        let mut treatment = treatment;
        let mut control = control;
        // `sort(key=itemgetter("lpos"))` -- Python's list sort is stable and so is
        // Rust's, so equal `lpos` keeps BAM order in both.
        treatment.sort_by_key(|r| r.lpos);
        control.sort_by_key(|r| r.lpos);

        // Upstream seeds from the treatment list and then widens over both.
        let mut read_start = treatment.iter().map(|r| r.lpos as i64).min().unwrap();
        let mut read_end = treatment.iter().map(|r| r.rpos as i64).max().unwrap();
        for r in treatment.iter().chain(control.iter()) {
            read_start = read_start.min(r.lpos as i64);
            read_end = read_end.max(r.rpos as i64);
        }
        read_start = read_start.min(peak_start);
        read_end = read_end.max(peak_end);

        let mut ext = vec![b'N'; (read_end - read_start).max(0) as usize];
        fill_refseq(&mut ext, read_start, &treatment)?;
        fill_refseq(&mut ext, read_start, &control)?;
        // `peak_refseq_ext[self.left - start : self.right - start]` -- a Python slice,
        // so both ends clamp to the (possibly resized) buffer. A shrunk buffer yields a
        // *short* `peak_refseq`, and the calling loop's `peak_refseq[i - left]` simply
        // stops being answerable past its end.
        let lo = clamp_slice(peak_start - read_start, ext.len());
        let hi = clamp_slice(peak_end - read_start, ext.len()).max(lo);
        let peak_refseq = ext[lo..hi].to_vec();

        Ok(Self {
            chrom: chrom.to_vec(),
            left: peak_start,
            right: peak_end,
            read_start,
            read_end,
            treatment,
            control,
            peak_refseq,
            peak_refseq_ext: ext,
        })
    }

    /// `n_edits_sum`.
    pub fn n_edits_sum(&self) -> u32 {
        self.treatment
            .iter()
            .chain(self.control.iter())
            .map(|r| r.n_edits())
            .sum()
    }

    /// `remove_outliers(percent=5)`, applied in place.
    ///
    /// Upstream builds the collection first and *then* removes outliers, so the
    /// consensus is built from the unfiltered reads. Re-deriving it afterwards would
    /// change the reference bases the calls are made against.
    pub fn remove_outliers_in_place(&mut self, percent: i64) {
        remove_outliers(&mut self.treatment, &mut self.control, percent);
    }

    /// Treatment reads, for tests and diagnostics.
    pub fn treatment_reads(&self) -> &[ReadAlignment] {
        &self.treatment
    }

    /// Control reads, for tests and diagnostics.
    pub fn control_reads(&self) -> &[ReadAlignment] {
        &self.control
    }

    /// `get_FASTQ`, for debugging exactly which reads the collection holds.
    pub fn fastq(&self) -> Vec<u8> {
        let mut out = Vec::new();
        for r in self.treatment.iter().chain(self.control.iter()) {
            out.extend_from_slice(b"@");
            out.extend_from_slice(&r.name);
            out.push(b'\n');
            out.extend_from_slice(&r.sequence());
            out.extend_from_slice(b"\n+\n");
            for &q in &r.qual {
                out.push(q + 33);
            }
            out.push(b'\n');
        }
        out
    }

    /// `get_PosReadsInfo_ref_pos`: what every read says about one position.
    ///
    /// The `Q` comparison is `read_bq <= Q` -> drop, applied here rather than inside
    /// `add_T` because upstream threads `Q` through as a parameter and the default is
    /// 20; `PosReadsInfo` hard-codes that default, so the caller has to filter first
    /// to honour a non-default `-Q`.
    pub fn pos_reads_info(&self, ref_pos: i64, ref_nt: &[u8], q: i32) -> PosReadsInfo {
        let mut p = PosReadsInfo::new(ref_pos as u64, ref_nt);
        for ra in &self.treatment {
            if ra.lpos as i64 <= ref_pos && ref_pos < ra.rpos as i64 {
                if let Ok(Some(v)) = ra.variant_bq_by_ref_pos(ref_pos as u64) {
                    if let Some(&bq) = v.bq.first() {
                        if bq as i32 > q {
                            p.add_t(&v.allele, bq as i32, v.strand, v.tip);
                        }
                    }
                }
            }
        }
        for ra in &self.control {
            if ra.lpos as i64 <= ref_pos && ref_pos < ra.rpos as i64 {
                if let Ok(Some(v)) = ra.variant_bq_by_ref_pos(ref_pos as u64) {
                    if let Some(&bq) = v.bq.first() {
                        if bq as i32 > q {
                            p.add_c(&v.allele, bq as i32);
                        }
                    }
                }
            }
        }
        p
    }
}

/// `__fill_refseq`, one sample's pass.
fn fill_refseq(seq: &mut Vec<u8>, start: i64, ralist: &[ReadAlignment]) -> Result<(), String> {
    if ralist.is_empty() {
        return Ok(());
    }
    let mut prev_r = ralist[0].lpos as i64;
    for i in 0..ralist.len() {
        let read = &ralist[i];
        if read.lpos as i64 > prev_r {
            // Upstream rebinds `read = ralist[i - 1]` and then computes the span from
            // *that* read, so the gap-filling write places read i-1's own bases at
            // read i-1's own coordinates. Reading it as "the previous read's bases
            // into the current read's span" is wrong, and it is off by the distance
            // between the two reads -- which is what left our consensus full of `N`
            // at the start of every peak.
            let prev = &ralist[i - 1];
            let refseq = prev.refseq()?;
            let ind = prev.lpos as i64 - start;
            let ind_r = ind + prev.rpos as i64 - prev.lpos as i64;
            py_splice(seq, ind, ind_r, &refseq)?;
            // And `prev_r = read["rpos"]` likewise reads the *previous* read. It is
            // deliberately not advanced for overlapping reads, so once a gap is found
            // the condition stays true for a while -- reproduced by only assigning it
            // here.
            prev_r = prev.rpos as i64;
        }
    }
    let last = ralist[ralist.len() - 1].clone();
    let refseq = last.refseq()?;
    let ind = last.lpos as i64 - start;
    let ind_r = ind + last.rpos as i64 - last.lpos as i64;
    py_splice(seq, ind, ind_r, &refseq)?;
    Ok(())
}

/// Python's slice-index clamping for a non-negative index.
fn clamp_slice(i: i64, len: usize) -> usize {
    (i.max(0) as u64).min(len as u64) as usize
}

/// Python's `bytearray` slice assignment, **including the resize**.
fn py_splice(buf: &mut Vec<u8>, a: i64, b: i64, value: &[u8]) -> Result<(), String> {
    let n = buf.len() as i64;
    // Python clamps both ends into the buffer.
    let a = if a < 0 { (n + a).max(0) } else { a.min(n) };
    let b = if b < 0 { (n + b).max(0) } else { b.min(n) }.max(a);
    if a > b {
        return Err(format!("fill_refseq: bad slice {a}..{b}"));
    }
    buf.splice(a as usize..b as usize, value.iter().copied());
    // A mis-sized replacement can in principle walk the buffer far away from its
    // intended size. Cap the drift so a pathological read set fails loudly instead of
    // allocating without bound.
    let expected = n as usize + value.len();
    let expected = expected.saturating_sub((b - a) as usize);
    if buf.len() > expected.saturating_mul(2).saturating_add(64) {
        return Err(format!(
            "fill_refseq: consensus grew to {} bytes, more than twice the {n} it \
             started at -- upstream's bytearray would have done this too, but the \\
             calling loop indexes it as peak_refseq[i - left] and cannot follow",
            buf.len()
        ));
    }
    Ok(())
}

/// Is the fermi-lite assembly path implemented?
///
/// Is the fermi-lite assembly path implemented *and* exact?
///
/// **Yes.** `--fermi auto` -- upstream's default -- is byte-identical to the pinned
/// oracle on its own `callvar_testing` fixtures: 16 of 16 records, across all ten peaks,
/// including the three emitted twice because the input holds two byte-identical peaks.
/// Checked by `oracle/check_callvar.sh --fermi`.
///
/// Reaching that took getting the *branch structure* right rather than the arithmetic:
/// when assembly is attempted, upstream holds the no-assembly variants **unwritten** and
/// then decides between dropping the peak, falling back, revisiting only the
/// reference-biased positions, or re-calling the whole peak from the unitigs. Getting
/// that last one wrong looked like a read-assignment bug -- DPT 12 where upstream reports
/// 15 -- and was chased through `remap_reads_with_unitigs` and the per-position
/// `ra_pos` filter before the real cause appeared. See F224.
pub const fn is_assembly_implemented() -> bool {
    true
}

/// `remove_outliers(percent=5)`: drop reads whose `n_edits` exceeds the
/// `(1 - percent/100)` quantile of the **combined** list.
pub fn remove_outliers(
    treatment: &mut Vec<ReadAlignment>,
    control: &mut Vec<ReadAlignment>,
    percent: i64,
) {
    let mut all: Vec<u32> = treatment
        .iter()
        .chain(control.iter())
        .map(|r| r.n_edits())
        .collect();
    if all.is_empty() {
        return;
    }
    all.sort_unstable();
    // `n_edits_list[int(len * (1 - percent*.01))]` -- a plain index, so `percent` of
    // 100 indexes past the end and upstream raises there. Clamping keeps the "zero
    // panics" contract.
    let idx = ((all.len() as f64) * (1.0 - percent as f64 * 0.01)) as usize;
    let cutoff = all[idx.min(all.len() - 1)];
    treatment.retain(|r| r.n_edits() <= cutoff);
    control.retain(|r| r.n_edits() <= cutoff);
}

/// `call_variants_at_range`: scan `[from, to)` and return `(pos, Variant)`.
///
/// # The order of the three filters
///
/// `update_top_alleles` may itself set `filterout`, in which case `call_GT` returns
/// immediately and `apply_GQ_cutoff` still runs (and no-ops). The GQ cutoff is applied
/// **only** on the branch where `call_GT` ran: a position filtered by
/// `update_top_alleles` never reaches it.
pub fn call_variants_at_range(
    coll: &RACollection,
    from: i64,
    to: i64,
    p: &CallParams,
) -> Result<Vec<(i64, Variant)>, String> {
    let mut out = Vec::new();
    for i in from..to {
        let Some(&nt) = coll.peak_refseq.get((i - coll.left) as usize) else {
            // `s[i - left]` past the end: upstream raises IndexError. The peak loop
            // already clamps to `[left, right)`, so this is only reachable when the
            // bytearray resize shrank the consensus below the peak; skipping keeps the
            // malformed-input contract.
            continue;
        };
        if nt == b'N' {
            continue;
        }
        let mut pri = coll.pos_reads_info(i, &[nt], p.min_q);
        if pri.raw_read_depth(DepthOpt::Treatment) == 0 {
            continue;
        }
        let live = pri.update_top_alleles(
            p.top2alleles_min_ratio,
            p.min_alt_allele_count,
            p.max_allowed_ar,
        );
        if live && !pri.filterout {
            pri.call_gt(p.max_allowed_ar)?;
            pri.apply_gq_cutoff(p.min_homo_gq, p.min_heter_gq);
        }
        if !pri.filterout {
            out.push((i, Variant::from_pos_reads_info(&pri)));
        }
    }
    Ok(out)
}

/// The `--fermi off` path for one peak: call, then repair indels.
///
/// The returned `bool` is whether the peak *would* have gone to fermi-lite assembly
/// under `--fermi auto`, i.e. `has_indel() or has_refer_biased_01()`.
///
/// # `fix_indels` is NOT applied here
///
/// Upstream applies it at *write* time, once, after deciding which branch it is in.
/// Applying it eagerly breaks `--fermi auto`: the no-assembly variants must still be
/// mutable when the refer-biased revisit replaces entries in them, and
/// `fix_indels` merges and removes positions. Use [`finish_peak`] instead, which calls
/// the repair and writes.
pub fn call_peak_without_assembly(
    coll: &RACollection,
    p: &CallParams,
) -> Result<(PeakVariants, bool), String> {
    let mut variants = PeakVariants::new(
        String::from_utf8_lossy(&coll.chrom).into_owned(),
        coll.left,
        coll.right,
        coll.peak_refseq.clone(),
    );
    for (pos, v) in call_variants_at_range(coll, coll.left, coll.right, p)? {
        variants.add_variant(pos, v);
    }
    #[allow(clippy::excessive_precision)]
    let wants_assembly = variants.has_indel() || variants.has_refer_biased_01(0.85f32);
    Ok((variants, wants_assembly))
}

/// Re-call every position in `[left, right)` from the unitig consensus.
///
/// This is upstream's `fermi == "on"`, or `--fermi auto` on a peak that carries an
/// indel: `peak_variants` is **reset** and rebuilt from the assembly, so the
/// no-assembly calls for that peak are discarded entirely.
pub fn call_peak_from_unitigs(
    coll: &RACollection,
    collection: &crate::unitig::UnitigCollection,
    p: &CallParams,
) -> Result<PeakVariants, String> {
    let mut variants = PeakVariants::new(
        String::from_utf8_lossy(&coll.chrom).into_owned(),
        coll.left,
        coll.right,
        coll.peak_refseq.clone(),
    );
    for pos in collection.left..collection.right {
        let Some(&nt) = coll.peak_refseq.get((pos - coll.left) as usize) else {
            continue;
        };
        if nt == b'N' {
            continue;
        }
        let mut pri = collection.pos_reads_info(pos, &[nt], p.min_q);
        if pri.raw_read_depth(crate::pos_reads_info::DepthOpt::Treatment) == 0 {
            continue;
        }
        if !pri.update_top_alleles(
            p.top2alleles_min_ratio,
            p.min_alt_allele_count,
            p.max_allowed_ar,
        ) {
            continue;
        }
        pri.call_gt(p.max_allowed_ar)?;
        pri.apply_gq_cutoff(p.min_homo_gq, p.min_heter_gq);
        if !pri.filterout {
            variants.add_variant(
                pos,
                crate::peak_variants::Variant::from_pos_reads_info(&pri),
            );
        }
    }
    Ok(variants)
}

/// `fix_indels` then write, exactly as every one of upstream's write sites does.
pub fn finish_peak(variants: &mut PeakVariants, ofile: &std::path::Path) -> Result<(), String> {
    if variants.n_variants() == 0 {
        return Ok(());
    }
    variants.fix_indels()?;
    use std::io::Write;
    let mut f = std::fs::OpenOptions::new()
        .append(true)
        .open(ofile)
        .map_err(|e| e.to_string())?;
    f.write_all(variants.to_vcf().as_bytes())
        .map_err(|e| e.to_string())
}

/// `fermi == "auto"`'s refer-biased revisit (`callvar_cmd.run`).
///
/// When the assembly ran but the peak has **no** indel, upstream re-calls every
/// reference-biased `0/1` from the unitig consensus rather than leaving the no-assembly
/// call in place -- a reference-biased het is exactly the artefact local assembly is
/// supposed to correct, so the no-assembly answer is known to be biased toward the
/// reference.
///
/// `Ok(false)` means the driver must drop the variant; `Ok(true)` means keep (and it may
/// have been replaced in place).
///
/// Note the order: `update_top_alleles` -> `call_GT` -> `apply_GQ_cutoff`, with **no**
/// `raw_read_depth == 0` skip. That check belongs to `call_variants_at_range`; the
/// revisit path goes straight to the model.
pub fn revisit_refer_biased(
    variants: &mut PeakVariants,
    collection: &crate::unitig::UnitigCollection,
    p: &CallParams,
) -> Result<(), String> {
    // A peak with an indel is handled by the assembly path instead.
    if variants.has_indel() {
        return Ok(());
    }
    for pos in variants.get_refer_biased_01s(REFER_BIASED_AR) {
        let Some(v) = variants.get(pos) else {
            continue;
        };
        let ref_nt = v.ref_allele.clone();
        if ref_nt == "N" {
            variants.remove_variant(pos)?;
            continue;
        }
        let mut pri = collection.pos_reads_info(pos, ref_nt.as_bytes(), p.min_q);
        if pri.raw_read_depth(crate::pos_reads_info::DepthOpt::Treatment) == 0 {
            variants.remove_variant(pos)?;
            continue;
        }
        pri.update_top_alleles(
            p.top2alleles_min_ratio,
            p.min_alt_allele_count,
            p.max_allowed_ar,
        );
        pri.call_gt(p.max_allowed_ar)?;
        pri.apply_gq_cutoff(p.min_homo_gq, p.min_heter_gq);
        if pri.filterout {
            variants.remove_variant(pos)?;
        } else {
            variants.replace_variant(
                pos,
                crate::peak_variants::Variant::from_pos_reads_info(&pri),
            )?;
        }
    }
    Ok(())
}

/// `PeakVariants.has_refer_biased_01()` takes no argument upstream, so
/// `is_refer_biased_01`'s own default of 0.85 is what applies.
pub const REFER_BIASED_AR: f32 = 0.85;

/// Build the unitig collection for one peak.
///
/// This is `build_unitig_collection` and nothing else -- it deliberately does **not**
/// re-call. Which variant set results is decided by the caller, because `--fermi auto`
/// on a peak with no indel keeps the no-assembly calls and only revisits the
/// reference-biased ones.
pub fn assemble_peak_with_unitigs(
    coll: &RACollection,
    fermi_min_overlap: i32,
) -> Result<crate::unitig::Built, String> {
    crate::unitig::build_unitig_collection(
        &coll.chrom,
        coll.left,
        coll.right,
        coll.read_start,
        coll.read_end,
        &coll.peak_refseq_ext,
        coll.treatment_reads(),
        coll.control_reads(),
        fermi_min_overlap,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An all-`M` alignment, so `refseq()` is exactly the read's own bases and its
    /// length always equals `ref_len`.
    fn ra_span(lpos: u32, ref_len: u32) -> ReadAlignment {
        // BAM packs two bases per byte, high nibble first, using indices into
        // "=ACMGRSVTWYHKDBN" -- so 'A' is nibble 1, not 0x41.
        const CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
        let qlen = ref_len as usize;
        let nibble = |b: u8| CODES.iter().position(|c| *c == b).unwrap() as u8;
        let mut seq = vec![0u8; qlen.div_ceil(2)];
        for i in 0..qlen {
            let n = nibble(b"ACGT"[i % 4]);
            if i % 2 == 0 {
                seq[i / 2] = n << 4;
            } else {
                seq[i / 2] |= n;
            }
        }
        ReadAlignment {
            name: b"r".to_vec(),
            chrom: b"chr1".to_vec(),
            lpos,
            rpos: lpos + ref_len,
            strand: 0,
            seq,
            qual: vec![30; qlen],
            // CIGAR word: `len << 4 | op`, op 0 = M.
            cigar: vec![ref_len << 4],
            md: format!("{ref_len}"),
        }
    }

    #[test]
    fn an_empty_treatment_list_is_refused() {
        let e = RACollection::new(b"chr1", 0, 10, vec![], vec![]).unwrap_err();
        assert_eq!(e, NO_READS);
    }

    #[test]
    fn the_read_span_widens_past_the_peak() {
        let c = RACollection::new(b"chr1", 100, 200, vec![ra_span(90, 10)], vec![]).unwrap();
        assert_eq!(c.read_start, 90);
        assert_eq!(c.read_end, 200, "the peak's own right end always widens it");
        assert_eq!(c.peak_refseq.len(), 100);
        assert_eq!(c.peak_refseq_ext.len(), 110);
    }

    #[test]
    fn the_reference_is_the_read_consensus_with_n_elsewhere() {
        let c = RACollection::new(b"chr1", 90, 96, vec![ra_span(90, 10)], vec![]).unwrap();
        assert_eq!(c.peak_refseq, b"ACGTAC".to_vec());
        assert!(!c.peak_refseq.contains(&b'N'));

        // A peak running past the read leaves `N` where nothing was filled.
        let d = RACollection::new(b"chr1", 90, 106, vec![ra_span(90, 10)], vec![]).unwrap();
        assert_eq!(&d.peak_refseq[10..], b"NNNNNN");
    }

    #[test]
    fn a_mis_sized_replacement_resizes_the_consensus_like_python() {
        // `__fill_refseq` writes the *previous* read's bases into the *current* read's
        // span. With a 10-byte span and an 11-byte value, Python's bytearray grows by
        // one and every later coordinate shifts. The port has to do the same or the
        // calling loop's `peak_refseq[i - left]` means something different.
        let mut buf = vec![b'N'; 100];
        py_splice(&mut buf, 10, 20, &[b'A'; 11]).unwrap();
        assert_eq!(
            buf.len(),
            101,
            "the buffer must resize, not overwrite in place"
        );
        assert_eq!(buf[10..21], [b'A'; 11]);
        assert_eq!(buf[21], b'N', "the tail shifted right by one");

        // Shrinking works the other way.
        py_splice(&mut buf, 10, 21, &[b'C'; 5]).unwrap();
        assert_eq!(buf.len(), 95);
    }

    #[test]
    fn py_splice_clamps_like_a_python_slice() {
        let mut buf = vec![b'N'; 10];
        // A slice past the end is clamped, not an error.
        py_splice(&mut buf, 8, 999, b"AA").unwrap();
        assert_eq!(buf.len(), 10);
        assert_eq!(&buf[8..], b"AA");
        // Python treats a reversed slice () as *empty*, so assigning to it
        // changes nothing. Clamping both ends to the same index reproduces that.
        let mut b2 = vec![b'Z'; 4];
        py_splice(&mut b2, 3, 1, b"").unwrap();
        assert_eq!(b2, b"ZZZZ".to_vec(), "an empty slice replaces nothing");
    }

    #[test]
    fn control_reads_fill_after_treatment() {
        let t = ra_span(100, 10);
        let mut ctl = ra_span(100, 10);
        ctl.seq = vec![0x88; 5]; // five bytes, high nibble first -> ten 'T'
        let coll = RACollection::new(b"chr1", 100, 110, vec![t], vec![ctl]).unwrap();
        assert_eq!(
            coll.peak_refseq,
            b"TTTTTTTTTT".to_vec(),
            "the control pass is applied second and overwrites treatment"
        );
    }

    #[test]
    fn remove_outliers_drops_only_the_extreme_tail() {
        let clean: Vec<ReadAlignment> = (0..20).map(|i| ra_span(i * 20, 10)).collect();
        let mut noisy = ra_span(1000, 10);
        noisy.md = "^ACGTACGTAC".to_string();
        let mut t = clean.clone();
        let mut c = clean.clone();
        c.push(noisy);
        remove_outliers(&mut t, &mut c, 5);
        assert_eq!(t.len(), 20);
        assert_eq!(c.len(), 20, "the noisy control read is the only casualty");
    }

    #[test]
    fn the_percentile_is_taken_over_both_samples_not_per_sample() {
        // 39 clean reads and exactly one with 5 edits.
        //
        // Combined (40 reads): idx = int(40 * 0.95) = 38, sorted[38] = 0, so the
        // single edited read is evicted. Per-sample treatment alone (20 reads):
        // idx = int(20 * 0.95) = 19, sorted_t[19] = 5, so nothing is evicted. That
        // difference is why this test exists.
        let mut edited = ra_span(1000, 10);
        edited.md = "^ACGTA".to_string();
        assert_eq!(edited.n_edits(), 5);

        let clean: Vec<ReadAlignment> = (0..39).map(|i| ra_span(i * 20, 10)).collect();
        let mut t = clean;
        t.push(edited);
        let mut c = Vec::new();
        remove_outliers(&mut t, &mut c, 5);
        assert_eq!(t.len(), 39, "the combined percentile evicts the outlier");
    }

    #[test]
    fn remove_outliers_never_panics() {
        let (mut t, mut c): (Vec<ReadAlignment>, Vec<ReadAlignment>) = (vec![], vec![]);
        remove_outliers(&mut t, &mut c, 5);
        assert!(t.is_empty() && c.is_empty());
        // A percent of 100 indexes past the end upstream; here it must not panic.
        let mut one = vec![ra_span(0, 10)];
        remove_outliers(&mut one, &mut Vec::new(), 100);
        assert_eq!(one.len(), 1);
    }

    #[test]
    fn fastq_covers_both_samples() {
        let coll = RACollection::new(b"chr1", 100, 110, vec![ra_span(100, 10)], vec![]).unwrap();
        let text = String::from_utf8(coll.fastq()).unwrap();
        assert_eq!(text.lines().count(), 4, "one read -> 4 FASTQ lines");
        assert!(text.starts_with("@r\n"));
        assert!(text.lines().nth(1).unwrap().starts_with("ACGT"));
    }

    #[test]
    fn a_unanimous_single_allele_peak_calls_nothing() {
        let coll = RACollection::new(b"chr1", 100, 120, vec![ra_span(100, 10)], vec![]).unwrap();
        let out =
            call_variants_at_range(&coll, coll.left, coll.right, &CallParams::default()).unwrap();
        assert!(
            out.is_empty(),
            "every covered position is homo_ref and filtered"
        );
    }

    #[test]
    fn an_n_base_is_never_called() {
        let coll = RACollection::new(
            b"chr1",
            100,
            115,
            vec![ra_span(100, 10), ra_span(105, 10)],
            vec![],
        )
        .unwrap();
        for (pos, _) in
            call_variants_at_range(&coll, coll.left, coll.right, &CallParams::default()).unwrap()
        {
            assert_ne!(
                coll.peak_refseq[(pos - coll.left) as usize],
                b'N',
                "at {pos}"
            );
        }
    }

    #[test]
    fn a_bad_md_tag_is_refused_rather_than_panicking() {
        let mut bad = ra_span(100, 10);
        bad.md = "3 3".to_string();
        let e = RACollection::new(b"chr1", 100, 110, vec![bad], vec![]).unwrap_err();
        assert!(e.contains("Don't understand this operator in MD"), "{e}");
    }

    #[test]
    fn n_edits_sum_covers_both_samples() {
        let coll = RACollection::new(
            b"chr1",
            100,
            110,
            vec![ra_span(100, 10)],
            vec![ra_span(100, 10)],
        )
        .unwrap();
        assert_eq!(coll.n_edits_sum(), 0);
        assert_eq!(coll.treatment_reads().len(), 1);
        assert_eq!(coll.control_reads().len(), 1);
    }
}
