//! Per-subcommand `--help` text, captured from the pinned oracle.
//!
//! # Why embedded text instead of generated help
//!
//! Every `<subcommand> --help` used to print a one-line stub
//! (`usage: macs3-rs <cmd> [options]`). Upstream prints full argparse help -- usage,
//! every flag, every default -- and users and scripts rely on it. A stub is a
//! drop-in incompatibility.
//!
//! These files are byte captures of `macs3 <cmd> --help` from the pinned oracle
//! (MACS3 3.0.5 at `c544319`, see `oracle/ENV.lock`), stored under `help/`.
//! They never go stale because the oracle never changes; if the pin moves, they
//! are re-captured.
//!
//! # The program name
//!
//! Upstream's text says `usage: macs3 ...`. This binary is `macs3-rs`, so the
//! usage line is rewritten at print time to the actual program name
//! ([`crate::PROGRAM`]). The option documentation is reproduced verbatim --
//! rewording it would risk documenting a flag this port parses differently from
//! the text.

use crate::PROGRAM;

/// Captured `--help` for `subcommand`, with the usage program name rewritten.
///
/// Returns `None` for unknown subcommands. The rewrite targets only the leading
/// `usage: macs3 ` prefix; `macs3` appearing inside prose (e.g. "MACS3 will ...")
/// is left alone, as is the `macs3 ` inside example command lines, which document
/// upstream's invocation rather than ours.
pub fn subcommand_help(subcommand: &str) -> Option<String> {
    let text: &str = match subcommand {
        "callpeak" => include_str!("../help/callpeak.txt"),
        "bdgpeakcall" => include_str!("../help/bdgpeakcall.txt"),
        "bdgbroadcall" => include_str!("../help/bdgbroadcall.txt"),
        "bdgcmp" => include_str!("../help/bdgcmp.txt"),
        "bdgopt" => include_str!("../help/bdgopt.txt"),
        "cmbreps" => include_str!("../help/cmbreps.txt"),
        "bdgdiff" => include_str!("../help/bdgdiff.txt"),
        "filterdup" => include_str!("../help/filterdup.txt"),
        "predictd" => include_str!("../help/predictd.txt"),
        "pileup" => include_str!("../help/pileup.txt"),
        "randsample" => include_str!("../help/randsample.txt"),
        "refinepeak" => include_str!("../help/refinepeak.txt"),
        "hmmratac" => include_str!("../help/hmmratac.txt"),
        "callvar" => include_str!("../help/callvar.txt"),
        _ => return None,
    };
    // Rewrite `usage: macs3 <subcommand>` to `usage: <PROGRAM> <subcommand>`.
    // Only the usage prefix is touched; the body keeps upstream's wording.
    let prefix = format!("usage: macs3 {subcommand}");
    if let Some(rest) = text.strip_prefix(&prefix) {
        Some(format!("usage: {PROGRAM} {subcommand}{rest}"))
    } else {
        Some(text.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every subcommand has captured help, and the usage line names this binary.
    #[test]
    fn every_subcommand_has_help_naming_this_binary() {
        for cmd in crate::SUBCOMMANDS {
            let h = subcommand_help(cmd).unwrap_or_else(|| panic!("no help for {cmd}"));
            assert!(
                h.starts_with(&format!("usage: {PROGRAM} {cmd}")),
                "help for {cmd} does not start with the usage line"
            );
            // The option documentation survived capture.
            assert!(h.contains("-h"), "help for {cmd} has no options");
        }
    }

    /// Unknown subcommands have no help text.
    #[test]
    fn unknown_subcommand_has_no_help() {
        assert!(subcommand_help("nosuchcmd").is_none());
    }
}
