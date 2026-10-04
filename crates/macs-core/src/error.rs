//! Error type shared by every `macs3-rs` crate.
//!
//! The contract is deliberately about *error conditions*, not error text: the
//! same malformed input must fail, and it must fail **before** any output file
//! is created. Upstream's human-readable messages are not reproduced verbatim;
//! see `docs/compatibility.md` §1.

use std::fmt;

/// Convenient result alias.
pub type Result<T> = std::result::Result<T, MacsError>;

/// Every failure mode that a user can trigger.
///
/// Variants map one-to-one onto upstream validation sites so that the
/// differential CLI harness (`tests/differential/cli_matrix`) can assert
/// accept/reject parity without string matching.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum MacsError {
    /// Malformed input record. `context` names the file/format.
    #[error("invalid {format} record at line {line}: {detail}")]
    InvalidFormat {
        format: &'static str,
        line: u64,
        detail: String,
    },

    /// Field could not be parsed as the required type.
    #[error("cannot parse {value:?} as {expected} in {format} at line {line}")]
    ParseInt {
        format: &'static str,
        line: u64,
        value: String,
        expected: &'static str,
    },

    /// Effective genome size neither numeric nor a known name.
    #[error("genome size {0:?} is neither a known genome name nor a positive integer")]
    InvalidGenomeSize(String),

    /// `--slocal` smaller than `d`, `--mfold` not length 3, bad cutoff order, ...
    #[error("invalid parameter: {0}")]
    InvalidParameter(String),

    /// Treatment and control share no chromosomes.
    #[error("no common chromosomes between treatment and control: {treat:?} vs {control:?}")]
    NoCommonChromosomes {
        treat: Vec<String>,
        control: Vec<String>,
    },

    /// `--keep-dup` value out of range.
    #[error("invalid --keep-dup {0}; expected a positive integer, 'all' or 'auto'")]
    InvalidKeepDup(String),

    /// A BAM/SAM/CRAM file was requested but its index is missing.
    #[error("index file not found for {path}; run samtools index first")]
    MissingIndex { path: String },

    /// BAM/SAM record exceeded what we can represent.
    #[error("alignment record too large or malformed: {0}")]
    BadAlignment(String),

    /// `predictd` could not build a model (too few paired strand peaks).
    #[error("failed to build fragment model: {0}")]
    ModelBuildFailed(String),

    /// The model could not be fitted because there were too few paired strand peaks.
    ///
    /// Distinct from [`MacsError::ModelBuildFailed`] because upstream does **not**
    /// treat it as a failure: `predictd_cmd.py:82-83` catches
    /// `NotEnoughPairsException`, emits a warning, and falls out of the function
    /// normally. So `macs3 predictd` exits **0** and writes nothing, and a caller
    /// that turns this into an error changes the exit status of a command upstream
    /// considers successful.
    #[error("can only find {found} paired peaks; MACS needs at least {needed}")]
    NotEnoughPairs { found: u64, needed: u64 },

    /// Paired-end mode was requested with a single-end input format, or vice versa.
    #[error("paired-end mode requires one of BAMPE, BEDPE, FRAG; got {0}")]
    PairedEndFormatRequired(&'static str),

    /// I/O failure.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),

    /// Output file already exists and would be silently clobbered.
    #[error("refusing to overwrite existing file: {0}")]
    OutputExists(String),

    /// Internal invariant violation. Reaching this is a bug in `macs3-rs`.
    #[error("internal error: {0}")]
    Internal(String),

    /// A rejection upstream performs with a bare `error_stream(...)` diagnostic,
    /// rendered here without an added prefix so the explanation matches.
    ///
    /// F27 and F172: `check_names` and `PeakDetect.__call_peaks_w_control` both
    /// print straight to stderr and then die with status 1. Exit status 1 is what
    /// the corpus asserts; the message is reproduced so a user sees the same
    /// explanation `macs3` gives.
    #[error("{0}")]
    Rejected(String),
}

impl MacsError {
    /// Convenience constructor for parse failures of integer fields.
    pub fn parse_int(format: &'static str, line: u64, value: &str, expected: &'static str) -> Self {
        MacsError::ParseInt {
            format,
            line,
            value: value.to_string(),
            expected,
        }
    }

    /// Convenience constructor for malformed records.
    pub fn invalid_format(format: &'static str, line: u64, detail: impl Into<String>) -> Self {
        MacsError::InvalidFormat {
            format,
            line,
            detail: detail.into(),
        }
    }

    /// Short machine-readable tag, used by the differential CLI harness to
    /// compare *conditions* between implementations without matching prose.
    pub fn tag(&self) -> &'static str {
        match self {
            MacsError::InvalidFormat { .. } => "invalid_format",
            MacsError::ParseInt { .. } => "parse_int",
            MacsError::InvalidGenomeSize(_) => "invalid_genome_size",
            MacsError::InvalidParameter(_) => "invalid_parameter",
            MacsError::NoCommonChromosomes { .. } => "no_common_chromosomes",
            MacsError::InvalidKeepDup(_) => "invalid_keep_dup",
            MacsError::MissingIndex { .. } => "missing_index",
            MacsError::BadAlignment(_) => "bad_alignment",
            MacsError::ModelBuildFailed(_) => "model_build_failed",
            MacsError::NotEnoughPairs { .. } => "not_enough_pairs",
            MacsError::PairedEndFormatRequired(_) => "paired_end_format_required",
            MacsError::Io(_) => "io",
            MacsError::OutputExists(_) => "output_exists",
            MacsError::Internal(_) => "internal",
            MacsError::Rejected(_) => "rejected",
        }
    }
}

/// `std::error::Error` for a bare message, used by `anyhow`-style contexts.
#[derive(Debug)]
pub struct Msg(pub String);

impl fmt::Display for Msg {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Msg {}
