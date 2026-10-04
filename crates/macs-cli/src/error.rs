//! Error types and exit-code policy for the `macs3` compatible CLI.
//!
//! The compatibility target is not just "the right numbers on the happy path" but
//! also *how the program fails*. Concretely, three properties are part of the
//! contract and are checked by the golden corpus:
//!
//! 1. the **exit code** — `0` on success, `1` for an upstream `sys.exit()` with no
//!    argument, `2` for an argparse usage error;
//! 2. **no output file is created** when a run fails. Upstream validates before
//!    it writes, and several of its failure paths are crashes that happen before
//!    the first `open()` (F27, F28), so a port that buffers lazily must not create
//!    a partial file where upstream created nothing;
//! 3. the **diagnostic text on stderr**, which is user-visible.
//!
//! # Tracebacks are not reproduced
//!
//! Two upstream failure paths are unhandled Python exceptions on perfectly valid
//! input:
//!
//! * **F27** — `check_names` formats a set of `bytes` with `",".join(...)`, raising
//!   `TypeError` after having already printed the two useful diagnostic lines.
//! * **F28** — `ratio_treat2control` divides by a zero control sum, raising
//!   `ZeroDivisionError`.
//!
//! For both, the port matches the **exit code**, the **printed lines**, and the
//! **absence of output files**, but emits a clean diagnostic instead of a Python
//! traceback. That was confirmed with the user. A traceback is not part of any
//! documented `macs3` interface, reproducing it would mean a Rust binary printing
//! Python frames it never ran, and the objective requires zero panics on
//! malformed input — an unhandled crash is the behaviour being retired, not
//! preserved.
//!
//! The mapping from variant to exit code is what the golden corpus asserts: 7685
//! runs over 425 fixtures, of which 7204 exit 0, 425 exit 2 (the deliberate
//! `--mfold 2 4 8` arity error) and 56 exit 1 (F27 and F28).

use std::fmt;

/// Why a run failed, in terms of what the user did.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FailureKind {
    /// A required argument was missing, or a flag was unrecognised or had the
    /// wrong arity. Upstream gets these from argparse, which exits 2.
    Usage,
    /// Upstream called `sys.exit()` after printing a diagnostic. Exit 1.
    Rejected,
    /// An input file could not be read or parsed. Upstream raises; exit 1.
    BadInput,
    /// A stage of the pipeline failed. Upstream raises; exit 1.
    Internal,
}

impl FailureKind {
    /// The process exit code upstream produces for this class of failure.
    pub const fn exit_code(self) -> i32 {
        match self {
            // argparse: "error: ...", then usage, then exit 2
            FailureKind::Usage => 2,
            // sys.exit() with no argument, or an unhandled exception
            FailureKind::Rejected | FailureKind::BadInput | FailureKind::Internal => 1,
        }
    }
}

/// A failure, with the diagnostic text that goes to stderr.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CliError {
    /// What class of failure this is, which fixes the exit code.
    pub kind: FailureKind,
    /// The diagnostic, without a trailing newline. Several lines are allowed.
    pub message: String,
}

impl CliError {
    /// A usage error: exit 2, as argparse produces.
    pub fn usage(message: impl Into<String>) -> Self {
        CliError {
            kind: FailureKind::Usage,
            message: message.into(),
        }
    }

    /// An upstream `sys.exit()` after a diagnostic: exit 1.
    pub fn rejected(message: impl Into<String>) -> Self {
        CliError {
            kind: FailureKind::Rejected,
            message: message.into(),
        }
    }

    /// An unreadable or malformed input: exit 1.
    pub fn bad_input(message: impl Into<String>) -> Self {
        CliError {
            kind: FailureKind::BadInput,
            message: message.into(),
        }
    }

    /// A pipeline stage failure: exit 1.
    pub fn internal(message: impl Into<String>) -> Self {
        CliError {
            kind: FailureKind::Internal,
            message: message.into(),
        }
    }

    /// The process exit code for this failure.
    pub fn exit_code(&self) -> i32 {
        self.kind.exit_code()
    }
}

impl fmt::Display for CliError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for CliError {}

/// The exact diagnostic upstream's `check_names` reaches the user with, before it
/// crashes on `",".join` of a set of `bytes` (F27).
///
/// Only the first two lines are ever seen by a user: the third line is where the
/// `TypeError` happens. Reproduced so a `macs3` invocation that finds no shared
/// chromosome produces the same explanation, then the same exit code.
pub const NO_COMMON_CHROMOSOMES: &str = "No common chromosome names can be found from treatment and control!\n\
Please make sure that the treatment and control alignment files were generated by using the same genome assembly!";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn usage_errors_exit_two_and_others_exit_one() {
        assert_eq!(CliError::usage("bad flag").exit_code(), 2);
        assert_eq!(CliError::rejected("no").exit_code(), 1);
        assert_eq!(CliError::bad_input("junk").exit_code(), 1);
        assert_eq!(CliError::internal("boom").exit_code(), 1);
    }

    #[test]
    fn the_no_common_chromosomes_diagnostic_is_two_lines() {
        let lines: Vec<&str> = NO_COMMON_CHROMOSOMES.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("No common chromosome names"));
        assert!(lines[1].starts_with("Please make sure"));
        // the crash site must not leak into the diagnostic
        assert!(!NO_COMMON_CHROMOSOMES.contains("TypeError"));
        assert!(!NO_COMMON_CHROMOSOMES.contains("Traceback"));
    }

    #[test]
    fn display_is_the_bare_message_with_no_python_framing() {
        let e = CliError::internal("stage failed");
        assert_eq!(e.to_string(), "stage failed");
        assert!(!e.to_string().contains("Traceback"));
    }
}
