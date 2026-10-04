//! Port of the rest of `MACS3/Signal/RACollection.py`'s assembly path and
//! `MACS3/Signal/UnitigRACollection.py`: the unitig bookkeeping between "fermi has
//! produced some sequences" and "here is a `PosReadsInfo` for reference position `i`".
//!
//! # The shape of the algorithm
//!
//! 1. assemble the peak's reads into unitigs (`fermi::assemble`);
//! 2. align each unitig back onto the peak consensus (`align::align_unitigs_to_reference`);
//! 3. **drop** unitigs whose alignment is too poor (`verify_alns`);
//! 4. assign each read to the first unitig that **contains its whole sequence**
//!    (`remap_RAs_w_unitigs`);
//! 5. if any treatment read is unmapped, assemble the unmapped ones and repeat --
//!    once at the full overlap and once at half (`build_unitig_collection`);
//! 6. convert each surviving unitig's alignment into an absolute `[lpos, rpos)` span.
//!
//! # Two rules that decide everything downstream
//!
//! **`verify_alns` normalises by alignment length, not unitig length.** The test is
//!
//! ```python
//! if aln_scores[i] * 100 / len(markup_alns[i]) < min_score_100:
//! ```
//!
//! i.e. score per 100 aligned columns, including the gap columns in the denominator.
//! `min_score_100` defaults to 150. Getting the denominator wrong changes which unitigs
//! survive, which changes the calls.
//!
//! **A read is assigned to the first unitig containing it, and only on an exact
//! substring match** (`tmp_ra_seq in unitig`). No alignment, no partial match. A read
//! that matches nothing goes to the unmapped list and is re-assembled in the next round.
//! Note `break` on first match: order matters, and unitig order is the order fermi
//! produced them in.

use crate::align::{align_unitigs_to_reference, revcomp, Alignments};
use crate::fermi::{assemble, Unitig, DEFAULT_OPT_FLAG};
use crate::pos_reads_info::PosReadsInfo;
use macs_io::bam::ReadAlignment;

/// `'-'` in the alignment strings.
const GAP: u8 = b'-';

/// Default `verify_alns` threshold: score per 100 aligned columns.
pub const MIN_SCORE_100: f64 = 150.0;

/// `verify_alns`: drop unitigs whose alignment scores below `min_score_100` per 100
/// aligned columns.
///
/// Upstream iterates **backwards** (`range(len-1, -1, -1)`) so that popping from a list
/// is safe. The result is order-preserving either way, but the walk order is kept so a
/// future side effect cannot depend on it.
pub fn verify_alns(unitigs: &mut Vec<Unitig>, alns: &mut Alignments, min_score_100: f64) {
    for i in (0..unitigs.len()).rev() {
        let Some(&score) = alns.scores.get(i) else {
            continue;
        };
        let Some(markup) = alns.markup_alns.get(i) else {
            continue;
        };
        if markup.is_empty() {
            // Upstream divides by `len(markup_alns[i])` unguarded, so an empty markup is
            // a ZeroDivisionError there. Refusing to keep it is the safe reading: a
            // zero-length alignment carries no evidence.
            unitigs.remove(i);
            alns.target_alns.remove(i);
            alns.reference_alns.remove(i);
            alns.scores.remove(i);
            alns.markup_alns.remove(i);
            continue;
        }
        let per_100 = score * 100.0 / markup.len() as f64;
        if per_100 < min_score_100 {
            unitigs.remove(i);
            alns.target_alns.remove(i);
            alns.reference_alns.remove(i);
            alns.scores.remove(i);
            alns.markup_alns.remove(i);
        }
    }
}

/// Reads split across unitigs, per `remap_RAs_w_unitigs`.
#[derive(Debug, Clone, Default)]
pub struct Remapped {
    /// Treatment reads per unitig.
    pub treatment: Vec<Vec<ReadAlignment>>,
    /// Control reads per unitig.
    pub control: Vec<Vec<ReadAlignment>>,
    /// Treatment reads that matched no unitig.
    pub unmapped_treatment: Vec<ReadAlignment>,
    /// Control reads that matched no unitig.
    pub unmapped_control: Vec<ReadAlignment>,
}

/// `remap_RAs_w_unitigs`: assign each read to the first unitig containing its sequence.
///
/// Exact substring containment, first match wins, no alignment. `ra.find(ra_seq)` in
/// upstream's later use of the mapping is the same containment test.
pub fn remap_reads_with_unitigs(
    unitigs: &[Unitig],
    treatment: &[ReadAlignment],
    control: &[ReadAlignment],
) -> Remapped {
    let mut out = Remapped {
        treatment: vec![Vec::new(); unitigs.len()],
        control: vec![Vec::new(); unitigs.len()],
        ..Default::default()
    };
    for ra in treatment {
        let seq = ra.sequence();
        match unitigs.iter().position(|u| contains(&u.seq, &seq)) {
            Some(i) => out.treatment[i].push(ra.clone()),
            None => out.unmapped_treatment.push(ra.clone()),
        }
    }
    for ra in control {
        let seq = ra.sequence();
        match unitigs.iter().position(|u| contains(&u.seq, &seq)) {
            Some(i) => out.control[i].push(ra.clone()),
            None => out.unmapped_control.push(ra.clone()),
        }
    }
    out
}

/// Python's `haystack.find(needle) != -1`, i.e. `needle in haystack`.
fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() {
        return true;
    }
    if needle.len() > haystack.len() {
        return false;
    }
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// `add_to_unitig_list`: append the second-round unitigs that are not already covered.
///
/// A second-round unitig `u0` is dropped when it is a substring of any existing unitig
/// `u1`, **or of `u1`'s reverse complement**. New unitigs are placed *before* the
/// existing ones (`new_unitig_list.extend(unitig_list)`), which changes the
/// first-match-wins order used by [`remap_reads_with_unitigs`].
pub fn add_to_unitig_list(unitig_list: &[Unitig], unitigs_2nd: &[Unitig]) -> Vec<Unitig> {
    let mut new_list: Vec<Unitig> = Vec::new();
    for u0 in unitigs_2nd {
        let mut covered = false;
        for u1 in unitig_list {
            if contains(&u1.seq, &u0.seq) {
                covered = true;
                break;
            }
            let rc = revcomp(&u1.seq);
            if contains(&rc, &u0.seq) {
                covered = true;
                break;
            }
        }
        if !covered {
            new_list.push(u0.clone());
        }
    }
    new_list.extend_from_slice(unitig_list);
    new_list
}

/// `UnitigRAs`: one unitig's placement on the consensus, plus the reads mapped to it.
#[derive(Debug, Clone)]
pub struct UnitigRAs {
    pub chrom: Vec<u8>,
    pub lpos: i64,
    pub rpos: i64,
    /// The alignment, with `-` for gaps. Same length as `reference_aln`.
    pub unitig_aln: Vec<u8>,
    /// The consensus side of the alignment.
    pub reference_aln: Vec<u8>,
    /// `unitig_aln` with gaps removed: the unitig's own sequence.
    pub seq: Vec<u8>,
    /// `[treatment, control]`.
    pub reads: [Vec<ReadAlignment>; 2],
}

impl UnitigRAs {
    /// `UnitigRAs.__init__`.
    ///
    /// `Err` when the two alignment strings differ in length -- upstream `assert`s it,
    /// and it would otherwise silently mis-pair the columns.
    pub fn new(
        chrom: &[u8],
        lpos: i64,
        rpos: i64,
        unitig_aln: Vec<u8>,
        reference_aln: Vec<u8>,
        reads: [Vec<ReadAlignment>; 2],
    ) -> Result<Self, String> {
        if unitig_aln.len() != reference_aln.len() {
            return Err(format!(
                "aln on unitig and reference should be the same length! \
                 ({} vs {})",
                unitig_aln.len(),
                reference_aln.len()
            ));
        }
        let seq: Vec<u8> = unitig_aln.iter().copied().filter(|b| *b != GAP).collect();
        Ok(Self {
            chrom: chrom.to_vec(),
            lpos,
            rpos,
            unitig_aln,
            reference_aln,
            seq,
            reads,
        })
    }

    /// Alignment length: the number of columns.
    pub fn aln_length(&self) -> usize {
        self.unitig_aln.len()
    }

    /// `get_variant_bq_by_ref_pos`: what the mapped reads say about `ref_pos`.
    ///
    /// Returns `(allele, bq_t, bq_c, strand_t, strand_c, tip_t, pos_t, pos_c)`.
    ///
    /// # The residue walk counts reference columns only
    ///
    /// ```python
    /// residue = ref_pos - self.lpos + 1
    /// for i in range(self.aln_length):
    ///     if self.reference_aln[i] != 45: residue -= 1
    ///     if residue == 0: break
    ///     index_aln += 1
    /// ```
    ///
    /// Gap columns in the *reference* do not consume a residue, so `index_aln` can point
    /// at a `-` in the unitig -- which is the deletion case below.
    #[allow(clippy::type_complexity)]
    pub fn variant_bq_by_ref_pos(
        &self,
        ref_pos: i64,
    ) -> (
        Vec<u8>,
        Vec<i32>,
        Vec<i32>,
        Vec<u8>,
        Vec<u8>,
        Vec<bool>,
        Vec<usize>,
        Vec<usize>,
    ) {
        let mut bq_t = Vec::new();
        let mut bq_c = Vec::new();
        let mut strand_t = Vec::new();
        let mut strand_c = Vec::new();
        let mut tip_t = Vec::new();
        let mut pos_t = Vec::new();
        let mut pos_c = Vec::new();

        let mut residue = ref_pos - self.lpos + 1;
        let mut index_aln = 0usize;
        for i in 0..self.aln_length() {
            if self.reference_aln[i] != GAP {
                residue -= 1;
            }
            if residue == 0 {
                break;
            }
            index_aln += 1;
        }
        // Upstream indexes `unitig_aln[index_aln:index_aln+1]`, which is empty rather
        // than an error if the walk ran off the end.
        let mut s: Vec<u8> = self
            .unitig_aln
            .get(index_aln..index_aln + 1)
            .map_or_else(Vec::new, |x| x.to_vec());

        // Unitig offset of the aligned column: gaps before it do not count.
        let index_unitig = self.unitig_aln[..(index_aln + 1).min(self.unitig_aln.len())]
            .iter()
            .filter(|b| **b != GAP)
            .count();

        // Deletion: every mapped read supports it, at fixed quality 93, and every read
        // is reported at its position in the unitig.
        if s.first() == Some(&GAP) {
            for ra in &self.reads[0] {
                let seq = ra.sequence();
                let found = find(&self.seq, &seq);
                // `self.seq.find(ra_seq)` is -1 when absent, which upstream subtracts as
                // -1 -- kept, because the subsequent `ra_pos < l_read and >= 0` filter
                // drops those reads anyway.
                let ra_pos = index_unitig as i64 - found - 1;
                let l_read = ra.length();
                tip_t.push(ra_pos == 0 || ra_pos == l_read as i64 - 1);
                bq_t.push(93);
                strand_t.push(ra.strand);
                pos_t.push(ra_pos.max(0) as usize);
            }
            for ra in &self.reads[1] {
                let seq = ra.sequence();
                let found = find(&self.seq, &seq);
                let ra_pos = index_unitig as i64 - found - 1;
                bq_c.push(93);
                strand_c.push(ra.strand);
                pos_c.push(ra_pos.max(0) as usize);
            }
            return (
                vec![b'*'],
                bq_t,
                bq_c,
                strand_t,
                strand_c,
                tip_t,
                pos_t,
                pos_c,
            );
        }

        // An insertion is the reference running out while the unitig has bases left.
        if index_aln < self.aln_length().saturating_sub(1) {
            for i in (index_aln + 1)..self.aln_length() {
                if self.reference_aln[i] == GAP {
                    s.extend_from_slice(&self.unitig_aln[i..i + 1]);
                } else {
                    break;
                }
            }
        }

        for ra in &self.reads[0] {
            let seq = ra.sequence();
            let l_read = ra.length();
            let ra_pos = index_unitig as i64 - find(&self.seq, &seq) - 1;
            if ra_pos < l_read as i64 && ra_pos >= 0 {
                pos_t.push(ra_pos as usize);
                tip_t.push(ra_pos == 0 || ra_pos == l_read as i64 - 1);
                bq_t.push(ra.qual[ra_pos as usize] as i32);
                strand_t.push(ra.strand);
            }
        }
        for ra in &self.reads[1] {
            let seq = ra.sequence();
            let l_read = ra.length();
            let ra_pos = index_unitig as i64 - find(&self.seq, &seq) - 1;
            if ra_pos < l_read as i64 && ra_pos >= 0 {
                pos_c.push(ra_pos as usize);
                bq_c.push(ra.qual[ra_pos as usize] as i32);
                strand_c.push(ra.strand);
            }
        }
        (s, bq_t, bq_c, strand_t, strand_c, tip_t, pos_t, pos_c)
    }
}

/// Python's `bytes.find`, returning `-1` when absent.
fn find(haystack: &[u8], needle: &[u8]) -> i64 {
    if needle.is_empty() {
        return 0;
    }
    if needle.len() > haystack.len() {
        return -1;
    }
    match haystack.windows(needle.len()).position(|w| w == needle) {
        Some(i) => i as i64,
        None => -1,
    }
}

/// `UnitigCollection`: the peak's unitigs, sorted by `lpos`.
#[derive(Debug, Clone)]
pub struct UnitigCollection {
    pub chrom: Vec<u8>,
    pub left: i64,
    pub right: i64,
    pub ura_list: Vec<UnitigRAs>,
}

impl UnitigCollection {
    pub fn new(chrom: &[u8], peak_start: i64, peak_end: i64, mut ura_list: Vec<UnitigRAs>) -> Self {
        // `URAs_list.sort(key=itemgetter("lpos"))` -- Python's sort is stable.
        ura_list.sort_by_key(|u| u.lpos);
        Self {
            chrom: chrom.to_vec(),
            left: peak_start,
            right: peak_end,
            ura_list,
        }
    }

    /// `URAs_left`.
    pub fn uras_left(&self) -> i64 {
        self.ura_list.iter().map(|u| u.lpos).min().unwrap_or(0)
    }

    /// `URAs_right`.
    pub fn uras_right(&self) -> i64 {
        self.ura_list.iter().map(|u| u.rpos).max().unwrap_or(0)
    }

    /// `get_PosReadsInfo_ref_pos`.
    pub fn pos_reads_info(&self, ref_pos: i64, ref_nt: &[u8], q: i32) -> PosReadsInfo {
        let mut p = PosReadsInfo::new(ref_pos as u64, ref_nt);
        for ura in &self.ura_list {
            if ura.lpos <= ref_pos && ref_pos < ura.rpos {
                let (s, bq_t, bq_c, strand_t, _strand_c, tip_t, _pos_t, _pos_c) =
                    ura.variant_bq_by_ref_pos(ref_pos);
                for j in 0..bq_t.len() {
                    if bq_t[j] > q {
                        p.add_t(&s, bq_t[j], strand_t[j], tip_t[j]);
                    }
                }
                // `add_C` takes no strand: `n_strand` is treatment-only, which is why
                // the control's strand list is discarded rather than stored.
                for &bq in &bq_c {
                    if bq > q {
                        p.add_c(&s, bq);
                    }
                }
            }
        }
        p
    }
}

/// What [`build_unitig_collection`] can report back.
#[derive(Debug)]
pub enum Built {
    /// `return 0`: fermi produced nothing at all.
    Empty,
    /// `return -1`: every unitig failed `verify_alns`.
    TooManyMismatches,
    /// A collection.
    Collection(Box<UnitigCollection>),
}

/// `build_unitig_collection`: the whole iterative assembly.
#[allow(clippy::too_many_arguments)]
pub fn build_unitig_collection(
    chrom: &[u8],
    peak_start: i64,
    peak_end: i64,
    read_start: i64,
    read_end: i64,
    peak_refseq_ext: &[u8],
    treatment: &[ReadAlignment],
    control: &[ReadAlignment],
    fermi_min_overlap: i32,
) -> Result<Built, String> {
    let reads_for = |t: &[ReadAlignment], c: &[ReadAlignment]| -> Vec<(Vec<u8>, Vec<u8>)> {
        t.iter()
            .chain(c.iter())
            .map(|r| (r.sequence(), r.qual.clone()))
            .collect()
    };

    let mut unitigs = assemble(
        &reads_for(treatment, control),
        fermi_min_overlap,
        DEFAULT_OPT_FLAG,
    )?;
    if unitigs.is_empty() {
        return Ok(Built::Empty);
    }
    let mut alns = align_unitigs_to_reference(&mut unitigs, peak_refseq_ext)?;
    verify_alns(&mut unitigs, &mut alns, MIN_SCORE_100);
    if unitigs.is_empty() {
        return Ok(Built::TooManyMismatches);
    }

    let mut mapped = remap_reads_with_unitigs(&unitigs, treatment, control);
    // The loop's exit condition is "no unmapped treatment reads AND the unitig count
    // stopped growing", so both are tracked exactly.
    let mut prev_n = usize::MAX;

    while !mapped.unmapped_treatment.is_empty() && prev_n != unitigs.len() {
        prev_n = unitigs.len();

        // Round 2a: the unmapped reads at the full overlap.
        let second = assemble(
            &reads_for(&mapped.unmapped_treatment, &mapped.unmapped_control),
            fermi_min_overlap,
            DEFAULT_OPT_FLAG,
        )?;
        if !second.is_empty() {
            unitigs = add_to_unitig_list(&unitigs, &second);
            alns = align_unitigs_to_reference(&mut unitigs, peak_refseq_ext)?;
            verify_alns(&mut unitigs, &mut alns, MIN_SCORE_100);
            mapped = remap_reads_with_unitigs(&unitigs, treatment, control);
        }

        // Round 2b: and again at *half* the overlap, unconditionally -- upstream runs
        // this second block whether or not the first produced anything.
        let half = fermi_min_overlap / 2;
        let second = assemble(
            &reads_for(&mapped.unmapped_treatment, &mapped.unmapped_control),
            half,
            DEFAULT_OPT_FLAG,
        )?;
        if !second.is_empty() {
            unitigs = add_to_unitig_list(&unitigs, &second);
            alns = align_unitigs_to_reference(&mut unitigs, peak_refseq_ext)?;
            verify_alns(&mut unitigs, &mut alns, MIN_SCORE_100);
            mapped = remap_reads_with_unitigs(&unitigs, treatment, control);
        }
    }

    if unitigs.is_empty() {
        // Upstream returns `None` here; the driver treats that as "reset to the
        // previous result", so `Empty` is the faithful mapping.
        return Ok(Built::Empty);
    }

    let start = peak_start.min(read_start);
    let end = peak_end.max(read_end);

    let mut ura_list = Vec::with_capacity(unitigs.len());
    for i in 0..unitigs.len() {
        let unitig_aln = alns.target_alns.get(i).cloned().unwrap_or_default();
        let reference_aln = alns.reference_alns.get(i).cloned().unwrap_or_default();

        // The reference side with its gaps dropped is what gets located in the
        // extended consensus.
        let ref_seq: Vec<u8> = reference_aln
            .iter()
            .copied()
            .filter(|b| *b != GAP)
            .collect();
        let left_padding_ref = match find(peak_refseq_ext, &ref_seq) {
            i if i >= 0 => i,
            // `bytes.find` returning -1 is upstream's most dangerous case: it makes
            // `left_padding_ref` -1 and the lpos arithmetic below goes wrong by one.
            // Refused rather than propagated.
            _ => {
                return Err(format!(
                    "build_unitig_collection: unitig's aligned reference not found in \
                     peak_refseq_ext (unitig length {})",
                    ref_seq.len()
                ))
            }
        };
        let right_padding_ref =
            peak_refseq_ext.len() as i64 - left_padding_ref - ref_seq.len() as i64;

        let left_padding_unitig = leading_gaps(&unitig_aln);
        let right_padding_unitig = trailing_gaps(&unitig_aln);

        let mut tmp_lpos = start + left_padding_ref;
        let mut tmp_rpos = end - right_padding_ref;

        for j in 0..left_padding_unitig {
            if reference_aln.get(j).copied().unwrap_or(GAP) != GAP {
                tmp_lpos += 1;
            }
        }
        for j in 1..=right_padding_unitig {
            let idx = reference_aln.len() as i64 - j as i64;
            if idx >= 0 && reference_aln[idx as usize] != GAP {
                tmp_rpos -= 1;
            }
        }

        let cut_unitig = &unitig_aln[left_padding_unitig..unitig_aln.len() - right_padding_unitig];
        let cut_ref =
            &reference_aln[left_padding_unitig..reference_aln.len() - right_padding_unitig];

        ura_list.push(UnitigRAs::new(
            chrom,
            tmp_lpos,
            tmp_rpos,
            cut_unitig.to_vec(),
            cut_ref.to_vec(),
            [mapped.treatment[i].clone(), mapped.control[i].clone()],
        )?);
    }

    Ok(Built::Collection(Box::new(UnitigCollection::new(
        chrom, peak_start, peak_end, ura_list,
    ))))
}

fn leading_gaps(s: &[u8]) -> usize {
    s.iter().take_while(|b| **b == GAP).count()
}

fn trailing_gaps(s: &[u8]) -> usize {
    s.iter().rev().take_while(|b| **b == GAP).count()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ra(lpos: u32, seq: &[u8], qual: u8, strand: u8) -> ReadAlignment {
        const CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
        let qlen = seq.len();
        let nibble = |b: u8| CODES.iter().position(|c| *c == b).unwrap() as u8;
        let mut packed = vec![0u8; qlen.div_ceil(2)];
        for (i, &b) in seq.iter().enumerate() {
            let n = nibble(b);
            if i % 2 == 0 {
                packed[i / 2] = n << 4;
            } else {
                packed[i / 2] |= n;
            }
        }
        ReadAlignment {
            name: b"r".to_vec(),
            chrom: b"chr1".to_vec(),
            lpos,
            rpos: lpos + qlen as u32,
            strand,
            seq: packed,
            qual: vec![qual; qlen],
            cigar: vec![(qlen as u32) << 4],
            md: format!("{qlen}"),
        }
    }

    fn unitig(seq: &[u8]) -> Unitig {
        Unitig {
            seq: seq.to_vec(),
            nsr: 1,
        }
    }

    #[test]
    fn verify_alns_normalises_by_alignment_length() {
        // A perfect 10-column alignment scores 20 -> 200 per 100, above the 150
        // threshold, so it survives.
        let mut u = vec![unitig(b"ACGTACGTAC")];
        let mut a = Alignments {
            target_alns: vec![b"ACGTACGTAC".to_vec()],
            reference_alns: vec![b"ACGTACGTAC".to_vec()],
            scores: vec![20.0],
            markup_alns: vec![b"||||||||||".to_vec()],
        };
        verify_alns(&mut u, &mut a, MIN_SCORE_100);
        assert_eq!(u.len(), 1, "a perfect short alignment must survive");

        // A poor one is dropped, and all five parallel vectors stay the same length.
        let mut u = vec![unitig(b"ACGTACGTAC"), unitig(b"TTTT")];
        let mut a = Alignments {
            target_alns: vec![b"ACGTACGTAC".to_vec(), b"TTTT".to_vec()],
            reference_alns: vec![b"ACGTACGTAC".to_vec(), b"GGGG".to_vec()],
            scores: vec![20.0, 2.0],
            markup_alns: vec![b"||||||||||".to_vec(), b"||||".to_vec()],
        };
        verify_alns(&mut u, &mut a, MIN_SCORE_100);
        assert_eq!(u.len(), 1);
        assert_eq!(a.target_alns.len(), 1);
        assert_eq!(a.reference_alns.len(), 1);
        assert_eq!(a.scores.len(), 1);
        assert_eq!(a.markup_alns.len(), 1);
        assert_eq!(u[0].seq, b"ACGTACGTAC");
    }

    #[test]
    fn the_denominator_is_columns_not_unitig_bases() {
        // 10 columns, 5 of them gaps, score 12: 12*100/10 = 120 < 150 -> dropped. Had
        // the denominator been the 5 non-gap bases it would have been 240 and kept.
        let mut u = vec![unitig(b"ACGTA")];
        let mut a = Alignments {
            target_alns: vec![b"ACGTA-----".to_vec()],
            reference_alns: vec![b"ACGTA-----".to_vec()],
            scores: vec![12.0],
            markup_alns: vec![b"|||||     ".to_vec()],
        };
        verify_alns(&mut u, &mut a, MIN_SCORE_100);
        assert!(u.is_empty(), "gapped columns count toward the denominator");
    }

    #[test]
    fn an_empty_markup_drops_the_unitig_rather_than_dividing_by_zero() {
        let mut u = vec![unitig(b"ACGT")];
        let mut a = Alignments {
            target_alns: vec![b"ACGT".to_vec()],
            reference_alns: vec![b"ACGT".to_vec()],
            scores: vec![8.0],
            markup_alns: vec![Vec::new()],
        };
        verify_alns(&mut u, &mut a, MIN_SCORE_100);
        assert!(u.is_empty());
    }

    #[test]
    fn reads_are_assigned_to_the_first_containing_unitig() {
        let us = vec![unitig(b"TTTTACGTACGTTTTT"), unitig(b"GGGGCCCC")];
        let (t, c) = (
            vec![ra(0, b"ACGTACGT", 30, 0), ra(10, b"CCCC", 30, 1)],
            vec![ra(20, b"GGGG", 30, 0)],
        );
        let m = remap_reads_with_unitigs(&us, &t, &c);
        assert_eq!(m.treatment.len(), 2);
        assert_eq!(m.treatment[0].len(), 1, "ACGTACGT is in unitig 0");
        assert_eq!(m.treatment[1].len(), 1, "CCCC is only in unitig 1");
        assert_eq!(m.control[1].len(), 1);
        assert!(m.unmapped_treatment.is_empty());
    }

    #[test]
    fn an_unmatched_read_goes_to_the_unmapped_list() {
        let us = vec![unitig(b"ACGTACGT")];
        let t = vec![ra(0, b"TTTTTTTT", 30, 0)];
        let m = remap_reads_with_unitigs(&us, &t, &[]);
        assert!(m.treatment[0].is_empty());
        assert_eq!(m.unmapped_treatment.len(), 1);
    }

    #[test]
    fn first_match_wins_even_when_a_later_unitig_also_contains_the_read() {
        // The read must be a substring of *both* unitigs for the order to matter;
        // "ACGT" does not contain "ACGTTT", so ordering would be unobservable there.
        let us = vec![unitig(b"ACGTTT"), unitig(b"TTACGTTTAA")];
        let t = vec![ra(0, b"ACGTTT", 30, 0)];
        let m = remap_reads_with_unitigs(&us, &t, &[]);
        assert_eq!(m.treatment[0].len(), 1, "the earlier unitig takes it");
        assert!(m.treatment[1].is_empty(), "and the later one gets nothing");
    }

    #[test]
    fn add_to_unitig_list_puts_new_unitigs_first_and_skips_covered_ones() {
        // revcomp("AAAACCCC") == "GGGGTTTT", so "GGGG" is covered by the first unitig
        // *only* through its reverse complement -- which is the case that a naive
        // "is it already in the list" check would miss.
        let existing = vec![unitig(b"AAAACCCC"), unitig(b"TTTTGGGG")];
        let second = vec![unitig(b"CCCC"), unitig(b"GGGG"), unitig(b"NNNN")];
        let out = add_to_unitig_list(&existing, &second);
        // "CCCC" is inside unitig 0 -> dropped.
        // "GGGG" is inside the *reverse complement* of unitig 0 -> dropped.
        // "NNNN" is new, and goes to the FRONT.
        assert_eq!(out.len(), 3, "got {} unitigs", out.len());
        assert_eq!(out[0].seq, b"NNNN");
        assert_eq!(out[1].seq, b"AAAACCCC");
        assert_eq!(out[2].seq, b"TTTTGGGG");
    }

    #[test]
    fn unitig_ras_requires_equal_alignment_lengths() {
        let e = UnitigRAs::new(
            b"chr1",
            0,
            10,
            b"ACGT".to_vec(),
            b"ACG".to_vec(),
            [vec![], vec![]],
        )
        .unwrap_err();
        assert!(e.contains("same length"), "{e}");
    }

    #[test]
    fn a_ura_strips_gaps_into_seq() {
        let u = UnitigRAs::new(
            b"chr1",
            0,
            10,
            b"AC-GT".to_vec(),
            b"ACGGT".to_vec(),
            [vec![], vec![]],
        )
        .unwrap();
        assert_eq!(u.seq, b"ACGT");
        assert_eq!(u.aln_length(), 5);
    }

    #[test]
    fn a_deletion_column_reports_star_at_93_for_every_mapped_read() {
        // unitig_aln has a gap at column 2; reference_aln does not, so the residue walk
        // lands on the gap.
        let reads = vec![ra(0, b"ACGT", 30, 0)];
        let u = UnitigRAs::new(
            b"chr1",
            0,
            4,
            b"AC-GT".to_vec(),
            b"ACAGT".to_vec(),
            [reads, vec![]],
        )
        .unwrap();
        let (s, bq_t, _, _, _, tip_t, _, _) = u.variant_bq_by_ref_pos(2);
        assert_eq!(s, b"*");
        assert_eq!(bq_t, vec![93], "a deletion takes the fixed quality");
        assert_eq!(tip_t.len(), 1);
    }

    #[test]
    fn an_insertion_appends_the_gapped_unitig_bases() {
        // An insertion is a gap in the *reference* (`reference_aln[i] == '-'`), not in
        // the unitig. Getting that backwards was my first attempt, and it silently
        // produced the matched base alone.
        let u = UnitigRAs::new(
            b"chr1",
            0,
            4,
            b"ACGT-".to_vec(),
            b"ACG-T".to_vec(),
            [vec![], vec![]],
        )
        .unwrap();
        let (s, _, _, _, _, _, _, _) = u.variant_bq_by_ref_pos(2);
        assert_eq!(s, b"GT", "the inserted base follows the match");
    }

    #[test]
    fn a_unitig_collection_sorts_by_lpos() {
        let mk = |l: i64, r: i64| {
            UnitigRAs::new(
                b"chr1",
                l,
                r,
                b"AC".to_vec(),
                b"AC".to_vec(),
                [vec![], vec![]],
            )
            .expect("equal-length alignments")
        };
        let c = UnitigCollection::new(b"chr1", 0, 100, vec![mk(50, 60), mk(10, 20), mk(30, 40)]);
        assert_eq!(c.uras_left(), 10);
        assert_eq!(c.uras_right(), 60);
        assert_eq!(c.ura_list[0].lpos, 10);
        assert_eq!(c.ura_list[2].lpos, 50);
    }

    #[test]
    fn a_position_outside_every_ura_yields_no_reads() {
        let c = UnitigCollection::new(
            b"chr1",
            0,
            100,
            vec![UnitigRAs::new(
                b"chr1",
                10,
                20,
                b"ACGT".to_vec(),
                b"ACGT".to_vec(),
                [vec![], vec![]],
            )
            .unwrap()],
        );
        let p = c.pos_reads_info(50, b"A", 20);
        assert_eq!(p.raw_read_depth(crate::pos_reads_info::DepthOpt::All), 0);
    }
}
