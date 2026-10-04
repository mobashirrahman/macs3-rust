//! Implementations of the `macs3` subcommands.
//!
//! Each command maps its already-validated [`crate::Options`] onto the library
//! and writes the files upstream writes. The flag surface is *derived* from the
//! auto-derived argparse matrix (`oracle/flag_matrix.tsv`) by
//! [`crate::parse_flags`], so this layer never re-implements argparse.

pub mod bdgcmp;
pub mod bdgpeakcall;
pub mod bedgraph_cmds;
pub mod callpeak;
pub mod callvar;
pub mod filterdup;
pub mod hmmratac;
pub mod input;
pub mod pileup;
pub mod predictd;
pub mod randsample;

pub mod refinepeak;
pub(crate) mod stagedump;
