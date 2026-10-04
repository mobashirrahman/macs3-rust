//! Core genomic types for `macs3-rs`.
//!
//! # Numerical contract (see `docs/numerical-contract.md`)
//!
//! * [`Coord`] is `u64`. Genomic offsets are non-negative; BED-family output is
//!   0-based half-open `[start, end)`. `u64` makes overflow of contigs larger
//!   than 4.29 Gbp and of offsets beyond 2^32 impossible.
//! * [`Count`] is `u64`. Counting in 64 bits removes the 32-bit overflow class of
//!   bug that upstream MACS3 has had to patch.
//! * Floating point values are `f64` internally, but upstream `f32` behaviour is
//!   emulated at the specific, enumerated call sites listed in the numerical
//!   contract. Silently "improving" precision is a compatibility regression,
//!   because a threshold crossing changes a genomic interval and therefore an
//!   entire peak.
//!
//! This crate performs no I/O and has no dependency beyond `thiserror`.

#![forbid(unsafe_code)]
#![deny(missing_debug_implementations)]

pub mod error;
pub mod genome;
pub mod interval;
pub mod peak;
pub mod strand;

pub use error::{MacsError, Result};
pub use genome::{ChromId, Genome};
pub use interval::Interval;
pub use peak::Peak;
pub use strand::Strand;

/// 0-based genomic offset. See module docs for the BED convention.
pub type Coord = u64;

/// Length of a genomic span.
pub type Len = u64;

/// A read / fragment / pileup count.
pub type Count = u64;

/// An integer weight applied to a fragment (usually `1`).
pub type Weight = i64;

/// A floating point score (p-value, q-value, fold enrichment, ...).
pub type Score = f64;

/// Default effective genome size table, mirroring `MACS3.Utilities.Genome`.
pub mod genomesize {
    /// A curated subset of upstream's `Genome.sizes` table.
    ///
    /// This is intentionally *not* the whole upstream table in this first
    /// revision; the full table is transcribed in
    /// `docs/numerical-contract.md` and added in bulk. Values are exact copies
    /// of upstream's `hg19.sizes` / `mm9.sizes` / `ce10.sizes` etc.
    pub const TABLE: &[(&str, u64)] = &[
        ("ce10", 121718998),
        ("ce11", 121685502),
        ("dm2", 162367812),
        ("dm3", 162367812),
        ("dm6", 142490431),
        ("hg17", 268551150),
        ("hg18", 274330980),
        ("hg19", 293128983),
        ("hg20", 293128983),
        ("hg38", 293128983),
        ("mm10", 265278350),
        ("mm9", 272115926),
        ("sacCer3", 121571434),
    ];

    /// Look up an effective genome size by name. Unknown names return `None`,
    /// which the CLI turns into the "must be numeric" validation error.
    pub fn lookup(name: &str) -> Option<u64> {
        TABLE
            .iter()
            .find(|(n, _)| n.eq_ignore_ascii_case(name))
            .map(|(_, s)| *s)
    }

    /// The four `-g` shortcuts, from `MACS3.Utilities.Constants.EFFECTIVEGS`.
    ///
    /// These are the *only* names `-g` accepts verbatim; anything else must be a
    /// number. They differ from the [`TABLE`] entries on purpose -- `hs` is
    /// GRCh38 (2 913 022 398), not the `hg38` row above, because upstream's
    /// table is keyed by *species* and is a different list. Getting these mixed
    /// up changes every lambda by a few percent.
    pub const EFFECTIVE_GS: &[(&str, u64)] = &[
        ("hs", 2_913_022_398), // GRCh38
        ("mm", 2_652_783_500), // GRCm38
        ("ce", 100_286_401),   // WBcel235
        ("dm", 142_573_017),   // dm6
    ];

    /// `opt_validate_callpeak`'s gsize resolution: try the shortcut table, then
    /// fall back to `float()`.
    ///
    /// Upstream's failure message names the available shortcuts, so the error
    /// text here lists the same four in the same order.
    pub fn resolve_gsize(spec: &str) -> Result<f64, String> {
        if let Some((_, v)) = EFFECTIVE_GS.iter().find(|(n, _)| *n == spec) {
            return Ok(*v as f64);
        }
        match spec.parse::<f64>() {
            Ok(v) => Ok(v),
            Err(_) => Err(format!(
                "Error when interpreting --gsize option: {spec}\nAvailable shortcuts of \
                 effective genome sizes are {}",
                EFFECTIVE_GS
                    .iter()
                    .map(|(n, _)| *n)
                    .collect::<Vec<_>>()
                    .join(",")
            )),
        }
    }
}

#[cfg(test)]
mod gsize_tests {
    use super::genomesize::*;

    #[test]
    fn the_four_shortcuts_match_upstream_constants() {
        assert_eq!(resolve_gsize("hs").unwrap(), 2_913_022_398.0);
        assert_eq!(resolve_gsize("mm").unwrap(), 2_652_783_500.0);
        assert_eq!(resolve_gsize("ce").unwrap(), 100_286_401.0);
        assert_eq!(resolve_gsize("dm").unwrap(), 142_573_017.0);
    }

    #[test]
    fn hs_is_grch38_not_the_hg38_table_row() {
        // The species shortcut and the assembly row are different lists; mixing
        // them shifts every lambda.
        assert_ne!(resolve_gsize("hs").unwrap(), lookup("hg38").unwrap() as f64);
    }

    #[test]
    fn numbers_pass_through() {
        assert_eq!(resolve_gsize("1e9").unwrap(), 1e9);
        assert_eq!(resolve_gsize("1000000000").unwrap(), 1e9);
        assert_eq!(resolve_gsize("0").unwrap(), 0.0);
    }

    #[test]
    fn a_shortcut_is_matched_exactly_not_case_insensitively() {
        // `efgsize[options.gsize]` is a plain dict lookup: "HS" is a KeyError.
        assert!(resolve_gsize("HS").is_err());
    }

    #[test]
    fn the_error_lists_the_available_shortcuts_in_order() {
        let e = resolve_gsize("nonsense").unwrap_err();
        assert!(
            e.starts_with("Error when interpreting --gsize option: nonsense"),
            "{e}"
        );
        assert!(e.ends_with("hs,mm,ce,dm"), "{e}");
    }
}
