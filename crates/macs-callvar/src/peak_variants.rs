//! Port of `MACS3/Signal/PeakVariants.py` -- the variant record `callvar` writes,
//! and the per-peak container that repairs indels.
//!
//! # The narrowing that actually reaches the VCF
//!
//! `PosReadsInfo::toVariant` feeds this module, and its constructor signature is
//! where two lossy conversions happen:
//!
//! | field | Python type at the call site | stored as |
//! |---|---|---|
//! | `GQ`, `PL_00`, `PL_01`, `PL_11` | `double` | `cython.int` -> **truncated** |
//! | `deltaBIC`, the four BICs, `AR` | `double` | `cython.float` -> **`f32`** |
//!
//! The truncation is harmless (it happens once, at the same place either way). The
//! `f32` narrowing is not: `toVCF` prints those six fields with `%.2f`, so a `f64`
//! value sitting a hair above a `x.xx5` boundary can round *down* after being stored
//! as `f32`. [`Variant`] therefore keeps them as `f32`, and [`Variant::to_vcf`]
//! widens them only for the format call, exactly as Cython does when it passes a
//! `float` to `printf("%f", ...)`.
//!
//! # `is_indel` is a substring test; `is_only_del` is equality
//!
//! ```python
//! def is_indel(self):    return self.v_mutation_type.find("Insertion") != -1 or ... == "Deletion"
//! def is_only_del(self): return self.v_mutation_type == "Deletion"
//! ```
//!
//! So a multiallelic `Deletion,SNV` counts as an indel but is neither "only a
//! deletion" nor "only an insertion". That asymmetry is load-bearing:
//! [`PeakVariants::fix_indels`] keys its three repair passes off `is_only_del` and
//! `is_only_insertion`, so a multiallelic indel is deliberately left alone.
//!
//! # `toVCF` is 1-based
//!
//! `str(p + 1)`, while `fix_indels` indexes `refseq` as `rs[p - start]` (0-based
//! within the peak). Both conventions appear in one module.

use crate::pos_reads_info::{DepthOpt, PosReadsInfo};
use std::collections::BTreeMap;

/// One variant call.
///
/// The six `f32` fields are not an oversight -- see the module docs.
#[derive(Debug, Clone, PartialEq)]
pub struct Variant {
    pub ref_allele: String,
    pub alt_allele: String,
    /// Truncated from `PosReadsInfo::gq` (`double` -> C `int`).
    pub gq: i32,
    pub filter: String,
    pub var_type: String,
    pub mutation_type: String,
    pub top1_allele: String,
    pub top2_allele: String,
    pub dpt: i32,
    pub dpc: i32,
    pub dp1t: i32,
    pub dp2t: i32,
    pub dp1c: i32,
    pub dp2c: i32,
    pub plus1t: i32,
    pub plus2t: i32,
    pub minus1t: i32,
    pub minus2t: i32,
    pub delta_bic: f32,
    pub bic_homo_major: f32,
    pub bic_homo_minor: f32,
    pub bic_heter_noas: f32,
    pub bic_heter_as: f32,
    pub ar: f32,
    pub gt: String,
    pub dp: i32,
    pub pl_00: i32,
    pub pl_01: i32,
    pub pl_11: i32,
}

/// `is_indel`: a **substring** search, so `Deletion,SNV` counts.
pub fn is_indel(mutation_type: &str) -> bool {
    mutation_type.contains("Insertion") || mutation_type.contains("Deletion")
}

/// `is_only_del`: exact equality.
pub fn is_only_del(mutation_type: &str) -> bool {
    mutation_type == "Deletion"
}

/// `is_only_insertion`: exact equality.
pub fn is_only_insertion(mutation_type: &str) -> bool {
    mutation_type == "Insertion"
}

impl Variant {
    /// `PosReadsInfo.toVariant`, including the two narrowing conversions.
    ///
    /// `cython.int` conversion truncates toward zero (C semantics), which is what
    /// `d as i32` does for the non-negative values GQ and PL always are.
    #[allow(clippy::excessive_precision)]
    pub fn from_pos_reads_info(p: &PosReadsInfo) -> Self {
        let t1 = p.top1_allele().to_vec();
        let t2 = p.top2_allele().to_vec();
        let (plus1, minus1) = p.strand_counts(&t1);
        let (plus2, minus2) = p.strand_counts(&t2);
        let dp1t = p.n_t(&t1);
        let dp2t = p.n_t(&t2);
        let ar = if dp1t + dp2t == 0 {
            f32::NAN
        } else {
            (dp1t as f64 / (dp1t + dp2t) as f64) as f32
        };
        Self {
            ref_allele: String::from_utf8_lossy(&p.ref_allele).into_owned(),
            alt_allele: String::from_utf8_lossy(&p.alt_allele).into_owned(),
            gq: p.gq as i32,
            filter: ".".to_string(),
            var_type: p.var_type.clone(),
            mutation_type: p.mutation_type.clone(),
            top1_allele: String::from_utf8_lossy(&t1).into_owned(),
            top2_allele: String::from_utf8_lossy(&t2).into_owned(),
            dpt: p.raw_read_depth(DepthOpt::Treatment),
            dpc: p.raw_read_depth(DepthOpt::Control),
            dp1t,
            dp2t,
            dp1c: p.n_c(&t1),
            dp2c: p.n_c(&t2),
            plus1t: plus1,
            plus2t: plus2,
            minus1t: minus1,
            minus2t: minus2,
            delta_bic: p.delta_bic as f32,
            bic_homo_major: p.bic_homo_major as f32,
            bic_homo_minor: p.bic_homo_minor as f32,
            bic_heter_noas: p.bic_heter_noas as f32,
            bic_heter_as: p.bic_heter_as as f32,
            ar,
            gt: p.gt.clone(),
            dp: p.raw_read_depth(DepthOpt::All),
            pl_00: p.pl_00 as i32,
            pl_01: p.pl_01 as i32,
            pl_11: p.pl_11 as i32,
        }
    }

    pub fn is_indel(&self) -> bool {
        is_indel(&self.mutation_type)
    }

    pub fn is_only_del(&self) -> bool {
        is_only_del(&self.mutation_type)
    }

    pub fn is_only_insertion(&self) -> bool {
        is_only_insertion(&self.mutation_type)
    }

    /// `top1isreference`.
    pub fn top1_is_reference(&self) -> bool {
        self.ref_allele == self.top1_allele
    }

    /// `top2isreference`.
    pub fn top2_is_reference(&self) -> bool {
        self.ref_allele == self.top2_allele
    }

    /// `is_refer_biased_01`. `ar` is a `cython.float`.
    #[allow(clippy::excessive_precision)]
    pub fn is_refer_biased_01(&self, ar: f32) -> bool {
        self.ar >= ar && self.ref_allele == self.top1_allele
    }

    /// `toVCF`: the seven body columns, no chromosome or position.
    ///
    /// The six `f32` fields are widened to `f64` for `%.2f`, matching a C `float`
    /// passed to `printf`'s `%f`.
    pub fn to_vcf(&self) -> String {
        let info = format!(
            "M={};MT={};DPT={};DPC={};DP1T={}{};DP2T={}{};DP1C={}{};DP2C={}{};\
             SB={},{},{},{};DBIC={:.2};BICHOMOMAJOR={:.2};BICHOMOMINOR={:.2};\
             BICHETERNOAS={:.2};BICHETERAS={:.2};AR={:.2}",
            self.var_type,
            self.mutation_type,
            self.dpt,
            self.dpc,
            self.dp1t,
            self.top1_allele,
            self.dp2t,
            self.top2_allele,
            self.dp1c,
            self.top1_allele,
            self.dp2c,
            self.top2_allele,
            self.plus1t,
            self.plus2t,
            self.minus1t,
            self.minus2t,
            self.delta_bic as f64,
            self.bic_homo_major as f64,
            self.bic_homo_minor as f64,
            self.bic_heter_noas as f64,
            self.bic_heter_as as f64,
            self.ar as f64,
        );
        let sample = format!(
            "{}:{}:{}:{},{},{}",
            self.gt, self.dp, self.gq, self.pl_00, self.pl_01, self.pl_11
        );
        [
            self.ref_allele.clone(),
            self.alt_allele.clone(),
            format!("{}", self.gq),
            self.filter.clone(),
            info,
            "GT:DP:GQ:PL".to_string(),
            sample,
        ]
        .join("\t")
    }
}

/// All variants found in one peak, keyed by 0-based reference position.
///
/// A [`BTreeMap`] is the equivalent of Python's `dict` here because every loop in
/// upstream is `for p in sorted(self.d_Variants.keys())` -- and crucially, that
/// `sorted(...)` is a **snapshot**, so keys inserted or removed during a pass are not
/// visited in that same pass. The code below takes explicit snapshots for the same
/// reason.
#[derive(Debug, Clone)]
pub struct PeakVariants {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    /// The peak's reference sequence, indexed `rs[p - start]` (0-based).
    pub refseq: Vec<u8>,
    variants: BTreeMap<i64, Variant>,
}

impl PeakVariants {
    pub fn new(chrom: impl Into<String>, start: i64, end: i64, s: Vec<u8>) -> Self {
        Self {
            chrom: chrom.into(),
            start,
            end,
            refseq: s,
            variants: BTreeMap::new(),
        }
    }

    pub fn n_variants(&self) -> usize {
        self.variants.len()
    }

    /// `add_variant`: `d_Variants[p] = v`, so an existing position is replaced.
    pub fn add_variant(&mut self, p: i64, v: Variant) {
        self.variants.insert(p, v);
    }

    pub fn get(&self, p: i64) -> Option<&Variant> {
        self.variants.get(&p)
    }

    pub fn contains(&self, p: i64) -> bool {
        self.variants.contains_key(&p)
    }

    /// `remove_variant`. Upstream `assert`s the key is present.
    pub fn remove_variant(&mut self, p: i64) -> Result<(), String> {
        self.variants
            .remove(&p)
            .map(|_| ())
            .ok_or_else(|| format!("remove_variant({p}) on a position with no variant"))
    }

    /// `replace_variant`. Upstream `assert`s the key is present.
    pub fn replace_variant(&mut self, p: i64, v: Variant) -> Result<(), String> {
        if !self.variants.contains_key(&p) {
            return Err(format!(
                "replace_variant({p}) on a position with no variant"
            ));
        }
        self.variants.insert(p, v);
        Ok(())
    }

    /// Positions in ascending order, matching `sorted(d_Variants.keys())`.
    pub fn positions(&self) -> Vec<i64> {
        self.variants.keys().copied().collect()
    }

    /// `has_indel`.
    pub fn has_indel(&self) -> bool {
        self.variants.values().any(Variant::is_indel)
    }

    /// `has_refer_biased_01`.
    #[allow(clippy::excessive_precision)]
    pub fn has_refer_biased_01(&self, ar: f32) -> bool {
        self.variants.values().any(|v| v.is_refer_biased_01(ar))
    }

    /// `get_refer_biased_01s`.
    #[allow(clippy::excessive_precision)]
    pub fn get_refer_biased_01s(&self, ar: f32) -> Vec<i64> {
        self.variants
            .iter()
            .filter(|(_, v)| v.is_refer_biased_01(ar))
            .map(|(p, _)| *p)
            .collect()
    }

    /// `fix_indels`: the three indel-repair passes.
    ///
    /// Returns `Err` where upstream raises. The one reachable case is a pure
    /// deletion at absolute position 0: pass 1 evaluates `d_Variants[p0]` with
    /// `p0 == -1` and `p == p1 + 1 == 0`, which is a `KeyError`. Surfacing it as an
    /// error rather than a panic keeps the "zero panics" contract while still
    /// refusing to invent output that upstream never produces.
    pub fn fix_indels(&mut self) -> Result<(), String> {
        self.merge_contiguous_deletions()?;
        self.anchor_leading_deletions()?;
        self.drop_deletion_after_insertion();
        Ok(())
    }

    /// Pass 1: fold a run of adjacent pure deletions into its first position.
    ///
    /// `p1` tracks the last *consumed* position (including one that was just
    /// removed), which is what makes a deleted position still extend the run.
    fn merge_contiguous_deletions(&mut self) -> Result<(), String> {
        let mut p0: i64 = -1;
        let mut p1: i64 = -1;
        for p in self.positions() {
            let continues = p == p1 + 1 && self.variants.get(&p).is_some_and(Variant::is_only_del);
            let head_ok = if continues {
                let head = self.variants.get(&p0).ok_or_else(|| {
                    format!("fix_indels: no variant at {p0} while merging at {p}")
                })?;
                head.is_only_del()
            } else {
                false
            };
            if continues && head_ok {
                let (t1, t2, ref_b) = {
                    let head = &self.variants[&p0];
                    let tail = &self.variants[&p];
                    let t1 = if head.top1_allele == "*" {
                        String::new()
                    } else {
                        head.top1_allele.clone()
                    };
                    let t2 = if head.top2_allele == "*" {
                        String::new()
                    } else {
                        head.top2_allele.clone()
                    };
                    (t1, t2, tail.ref_allele.clone())
                };
                let head_ref = self.variants[&p].ref_allele.clone();
                let head = self.variants.get_mut(&p0).expect("checked above");
                // Only the *reference* allele is extended, and only on whichever of
                // top1/top2 is the reference.
                if head.top1_is_reference() {
                    head.top1_allele = t1 + &head_ref;
                } else if head.top2_is_reference() {
                    head.top2_allele = t2 + &head_ref;
                }
                head.ref_allele += &ref_b;
                self.variants.remove(&p);
                p1 = p;
            } else {
                p0 = p;
                p1 = p;
            }
        }
        Ok(())
    }

    /// Pass 2: a deletion with no variant at `p-1` borrows the preceding reference
    /// base, so the VCF anchor exists.
    fn anchor_leading_deletions(&mut self) -> Result<(), String> {
        for p in self.positions() {
            let is_del = self.variants.get(&p).is_some_and(Variant::is_only_del);
            if !is_del || self.variants.contains_key(&(p - 1)) || p <= self.start {
                continue;
            }
            let offset = (p - self.start) as usize;
            let base = *self
                .refseq
                .get(offset)
                .ok_or_else(|| format!("fix_indels: refseq too short for position {p}"))?;
            let base = (base as char).to_string();

            // The branch is chosen from the **original** variant at `p`, before its
            // REF/ALT are rewritten. Testing the rewritten copy instead flips the
            // answer: prepending the base to REF almost always makes `ref != top1`,
            // so top1/top2 would be left stale (still holding the old `*`).
            let top1_ref = self.variants[&p].top1_is_reference();
            let top2_ref = self.variants[&p].top2_is_reference();

            let mut moved = self.variants[&p].clone();
            moved.ref_allele = format!("{}{}", base, moved.ref_allele);
            moved.alt_allele = base.clone();
            if top1_ref {
                moved.top1_allele = moved.ref_allele.clone();
                moved.top2_allele = moved.alt_allele.clone();
            } else if top2_ref {
                moved.top1_allele = moved.alt_allele.clone();
                moved.top2_allele = moved.ref_allele.clone();
            }
            self.variants.insert(p - 1, moved);
            self.variants.remove(&p);
        }
        Ok(())
    }

    /// Pass 3: a deletion immediately after an insertion means either a third
    /// genotype or a simple repeat; upstream drops both.
    fn drop_deletion_after_insertion(&mut self) {
        for p in self.positions() {
            let is_del = self.variants.get(&p).is_some_and(Variant::is_only_del);
            if !is_del {
                continue;
            }
            let prev_is_ins = self
                .variants
                .get(&(p - 1))
                .is_some_and(Variant::is_only_insertion);
            if prev_is_ins {
                self.variants.remove(&p);
                self.variants.remove(&(p - 1));
            }
        }
    }

    /// `toVCF`: one line per variant, `\t`-joined with a trailing newline.
    ///
    /// The position is `p + 1` -- VCF is 1-based even though the keys are 0-based.
    pub fn to_vcf(&self) -> String {
        let mut out = String::new();
        for (p, v) in &self.variants {
            out.push_str(&self.chrom);
            out.push('\t');
            out.push_str(&(p + 1).to_string());
            out.push_str("\t.\t");
            out.push_str(&v.to_vcf());
            out.push('\n');
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn variant(ref_allele: &str, alt: &str, mt: &str) -> Variant {
        Variant {
            ref_allele: ref_allele.to_string(),
            alt_allele: alt.to_string(),
            gq: 58,
            filter: ".".to_string(),
            var_type: "heter_unsure".to_string(),
            mutation_type: mt.to_string(),
            top1_allele: ref_allele.to_string(),
            top2_allele: alt.to_string(),
            dpt: 7,
            dpc: 0,
            dp1t: 5,
            dp2t: 2,
            dp1c: 0,
            dp2c: 0,
            plus1t: 0,
            plus2t: 0,
            minus1t: 5,
            minus2t: 2,
            delta_bic: 23.21,
            bic_homo_major: 37.77,
            bic_homo_minor: 84.27,
            bic_heter_noas: 13.53,
            bic_heter_as: 14.56,
            ar: 0.714,
            gt: "0/1".to_string(),
            dp: 7,
            pl_00: 159,
            pl_01: 0,
            pl_11: 58,
        }
    }

    #[test]
    fn is_indel_is_substring_but_only_del_is_equality() {
        assert!(Variant::is_indel_mt("Deletion"));
        assert!(Variant::is_indel_mt("Deletion,SNV"));
        assert!(Variant::is_indel_mt("Insertion,SNV"));
        assert!(!Variant::is_indel_mt("SNV"));
        // A multiallelic indel is an indel but neither "only" kind, which is why
        // fix_indels leaves it alone.
        assert!(!is_only_del("Deletion,SNV"));
        assert!(is_only_del("Deletion"));
        assert!(!is_only_insertion("Deletion,SNV"));
        assert!(is_only_insertion("Insertion"));
    }

    /// Helper so the substring/equality asymmetry is stated once in a test name.
    impl Variant {
        #[allow(non_snake_case)]
        fn is_indel_mt(mt: &str) -> bool {
            is_indel(mt)
        }
    }

    #[test]
    fn vcf_body_has_seven_columns_and_all_info_keys() {
        let v = variant("A", "G", "SNV");
        let line = v.to_vcf();
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 7, "{line}");
        assert_eq!(f[0], "A");
        assert_eq!(f[1], "G");
        assert_eq!(f[2], "58");
        assert_eq!(f[3], ".");
        assert_eq!(f[5], "GT:DP:GQ:PL");
        assert_eq!(f[6], "0/1:7:58:159,0,58");
        assert!(f[4].starts_with("M=heter_unsure;MT=SNV;DPT=7;DPC=0;DP1T=5A;DP2T=2G;"));
        assert!(f[4].contains("SB=0,0,5,2;"));
        assert!(f[4].contains("DBIC=23.21;"));
        assert!(f[4].ends_with("AR=0.71"), "{}", f[4]);
    }

    #[test]
    fn floats_are_stored_as_f32_and_widened_for_printing() {
        // 0.714 is not representable in f32; the %.2f must still be 0.71.
        let mut v = variant("A", "G", "SNV");
        v.ar = 5.0f32 / 7.0f32;
        assert!(v.to_vcf().contains("AR=0.71"), "{}", v.to_vcf());
        // The narrowing is not cosmetic. These are values where the f64 and the f32
        // renderings of the *same* number differ in the second decimal, so keeping
        // them as f64 would produce a visibly different INFO column:
        //
        //     value   as f64    as f32
        //     23.215  23.21     23.22
        //     13.535  13.54     13.53
        //    138.185 138.19    138.18
        for (value, f32_text) in [
            (23.215f32, "DBIC=23.22"),
            (13.535, "DBIC=13.53"),
            (138.185, "DBIC=138.18"),
        ] {
            v.delta_bic = value;
            assert!(v.to_vcf().contains(f32_text), "{}", v.to_vcf());
            // And the f64 rendering is the one we must *not* produce.
            let f64_text = format!("DBIC={:.2}", value as f64);
            if f64_text != f32_text {
                assert!(
                    !v.to_vcf().contains(f64_text.as_str()),
                    "f32 storage must not fall back to the f64 rendering {f64_text}"
                );
            }
        }
    }

    #[test]
    fn peak_vcf_is_one_based_and_sorted() {
        let mut pv = PeakVariants::new("chr22", 100, 200, vec![b'A'; 100]);
        pv.add_variant(105, variant("A", "G", "SNV"));
        pv.add_variant(101, variant("C", "T", "SNV"));
        let out = pv.to_vcf();
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("chr22\t102\t.\tC\tT"), "{}", lines[0]);
        assert!(lines[1].starts_with("chr22\t106\t.\tA\tG"), "{}", lines[1]);
        assert!(out.ends_with('\n'));
    }

    #[test]
    fn adjacent_deletions_merge_into_the_first_position() {
        let mut pv = PeakVariants::new("chr22", 100, 200, vec![b'A'; 100]);
        pv.add_variant(110, variant("A", "*", "Deletion"));
        pv.add_variant(111, variant("C", "*", "Deletion"));
        pv.add_variant(112, variant("G", "*", "Deletion"));
        pv.fix_indels().unwrap();
        // Pass 1 folds 110/111/112 into 110 (REF "ACG"); pass 2 then anchors that
        // deletion at 109 using refseq[10] == 'A'. So the survivor is 109, not 110.
        assert_eq!(pv.n_variants(), 1);
        assert!(!pv.contains(110));
        let v = pv.get(109).unwrap();
        assert_eq!(
            v.ref_allele, "AACG",
            "merged REF then prefixed with the anchor base"
        );
        assert_eq!(v.alt_allele, "A");
        assert_eq!(v.mutation_type, "Deletion");
        // The branch came from the ORIGINAL variant, whose top1 was its REF, so
        // top1/top2 were rewritten rather than left holding the stale `*`.
        assert_eq!(v.top1_allele, "AACG");
        assert_eq!(v.top2_allele, "A");
    }

    #[test]
    fn non_adjacent_deletions_do_not_merge() {
        let mut pv = PeakVariants::new("chr22", 100, 200, vec![b'A'; 100]);
        pv.add_variant(110, variant("A", "*", "Deletion"));
        pv.add_variant(112, variant("G", "*", "Deletion"));
        pv.fix_indels().unwrap();
        assert_eq!(pv.n_variants(), 2, "a gap of one base breaks the run");
    }

    #[test]
    fn a_deletion_gets_a_reference_base_anchored_before_it() {
        // Only position 130 exists, so pass 2 synthesises 129 from refseq.
        let mut pv = PeakVariants::new("chr22", 100, 200, vec![b'A'; 100]);
        pv.add_variant(130, variant("A", "*", "Deletion"));
        pv.fix_indels().unwrap();
        assert_eq!(pv.n_variants(), 1);
        let v = pv.get(129).unwrap();
        assert_eq!(v.ref_allele, "AA");
        assert_eq!(v.alt_allele, "A");
        assert_eq!(v.mutation_type, "Deletion");
    }

    #[test]
    fn a_deletion_at_the_peak_start_is_left_alone() {
        // `p > self.start` guards pass 2, so nothing is synthesised before the peak.
        let mut pv = PeakVariants::new("chr22", 100, 200, vec![b'A'; 100]);
        pv.add_variant(100, variant("A", "*", "Deletion"));
        pv.fix_indels().unwrap();
        assert!(pv.contains(100));
        assert!(!pv.contains(99));
    }

    #[test]
    fn a_deletion_following_an_insertion_removes_both() {
        let mut pv = PeakVariants::new("chr22", 100, 200, vec![b'A'; 100]);
        pv.add_variant(140, variant("A", "ACGT", "Insertion"));
        pv.add_variant(141, variant("C", "*", "Deletion"));
        pv.fix_indels().unwrap();
        assert_eq!(pv.n_variants(), 0, "both are dropped together");
    }

    #[test]
    fn multiallelic_indels_are_never_merged_or_dropped() {
        let mut pv = PeakVariants::new("chr22", 100, 200, vec![b'A'; 100]);
        pv.add_variant(150, variant("A", "*,G", "Deletion,SNV"));
        pv.add_variant(151, variant("C", "*,T", "Deletion,SNV"));
        pv.fix_indels().unwrap();
        assert_eq!(
            pv.n_variants(),
            2,
            "is_only_del is false for a multiallelic type"
        );
    }

    #[test]
    fn refer_biased_het_calls_are_reported_in_position_order() {
        let mut pv = PeakVariants::new("chr22", 100, 200, vec![b'A'; 100]);
        let mut v1 = variant("A", "G", "SNV");
        v1.ar = 0.90;
        let mut v2 = variant("C", "T", "SNV");
        v2.ar = 0.10;
        pv.add_variant(120, v1);
        pv.add_variant(115, v2);
        assert!(pv.has_refer_biased_01(0.85));
        assert_eq!(pv.get_refer_biased_01s(0.85), vec![120]);
        assert!(!pv.has_refer_biased_01(0.95), "0.90 does not clear 0.95");
    }

    #[test]
    fn has_indel_sees_the_multiallelic_case() {
        let mut pv = PeakVariants::new("chr22", 100, 200, vec![b'A'; 100]);
        pv.add_variant(120, variant("A", "G", "SNV"));
        assert!(!pv.has_indel());
        pv.add_variant(121, variant("A", "*,G", "Deletion,SNV"));
        assert!(pv.has_indel());
    }

    #[test]
    fn removing_a_missing_variant_is_an_error_not_a_panic() {
        let mut pv = PeakVariants::new("chr22", 100, 200, vec![b'A'; 100]);
        assert!(pv.remove_variant(5).is_err());
        assert!(pv.replace_variant(5, variant("A", "G", "SNV")).is_err());
    }

    #[test]
    fn a_short_refseq_is_an_error_not_a_panic() {
        let mut pv = PeakVariants::new("chr22", 100, 105, vec![b'A'; 5]);
        pv.add_variant(130, variant("A", "*", "Deletion"));
        assert!(pv.fix_indels().is_err());
    }
}
