//! Finding the pinned MACS3 oracle a differential test needs.
//!
//! The oracle is deliberately not vendored: it is a git checkout of MACS3 3.0.5
//! plus a Python environment with NumPy and the compiled extensions, which no
//! clone of this repository can carry. `oracle/ENV.lock` records what the pin *is*
//! (version, commit, dependency versions) and never where it lives;
//! `oracle/provision_oracle.sh` records the location of a local provisioning in
//! `oracle/ENV.provisioned` and exports it to later steps on CI.
//!
//! So a test that needs the reference looks it up through the environment only:
//!
//! | variable               | meaning                                            |
//! |------------------------|----------------------------------------------------|
//! | `MACS3_SRC`            | the pinned checkout, e.g. `.oracle/macs3-src`      |
//! | `MACS3_ORACLE_BIN`     | the oracle's `macs3` entry point                   |
//! | `MACS3_ORACLE_PYTHON`  | the interpreter that can import it                  |
//! | `MACS3_VENV`           | the virtualenv holding that interpreter and NumPy  |
//! | `MACS3_RS_NO_ORACLE`   | skip even when one of the above resolves           |
//!
//! Each falls back to `oracle/ENV.provisioned` and then to the location
//! `provision_oracle.sh` creates by default (`<repo>/.oracle`), so a provisioned
//! checkout works without exporting anything.
//!
//! With no reference installed, [`require`] prints a skip line and the caller
//! returns successfully: CI's unit job runs on a stock Rust image with no Python
//! at all, and the release definition requires the shipped tool not to need one.
//! That is a *missing* reference, which is not a divergence. A reference that is
//! present and wrong is still a hard failure: nothing here swallows the
//! comparison itself.

// Each crate's test binary uses a different subset of this.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// A provisioned oracle: the pinned source tree, and what can run it.
pub struct Oracle {
    /// The pinned MACS3 checkout.
    pub src: PathBuf,
    /// The oracle's `macs3` entry point, if one was found.
    pub bin: Option<PathBuf>,
    /// An interpreter that can `import MACS3` -- the provisioned virtualenv is the
    /// only place NumPy and the compiled extensions are installed.
    pub python: Option<PathBuf>,
}

impl Oracle {
    /// Upstream's own test data, a sibling of the `MACS3` package inside the pinned
    /// checkout. Reading it from there is what keeps multi-megabyte BAMs out of this
    /// repository.
    pub fn test_dir(&self) -> PathBuf {
        self.src.join("test")
    }
}

/// This crate's repository root.
fn repo_root() -> PathBuf {
    // `<repo>/crates/<crate>`; the ancestors one and two up are `crates` and `<repo>`.
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("every crate lives at <repo>/crates/<name>")
        .to_path_buf()
}

/// `KEY=value` from a shell-style record, ignoring comments and empty values.
fn record_value(record: &Path, key: &str) -> Option<PathBuf> {
    let prefix = format!("{key}=");
    std::fs::read_to_string(record)
        .ok()?
        .lines()
        .find_map(|line| line.strip_prefix(&prefix))
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(PathBuf::from)
}

/// The checkout holding `MACS3/`, given a recorded `MACS3_PATH`.
fn checkout_of_init(path: &Path) -> PathBuf {
    // MACS3_PATH points at `<checkout>/MACS3/__init__.py`.
    path.parent()
        .and_then(Path::parent)
        .unwrap_or(path)
        .to_path_buf()
}

/// The first candidate location that is really there.
fn first_existing<I, F>(candidates: I, is_there: F) -> Option<PathBuf>
where
    I: IntoIterator<Item = PathBuf>,
    F: Fn(&Path) -> bool,
{
    candidates.into_iter().find(|candidate| is_there(candidate))
}

/// The provisioned reference, or `None` when this machine has none.
pub fn oracle() -> Option<Oracle> {
    if std::env::var_os("MACS3_RS_NO_ORACLE").is_some() {
        return None;
    }
    let root = repo_root();
    let provisioned = root.join("oracle").join("ENV.provisioned");
    let src = first_existing(
        [
            std::env::var_os("MACS3_SRC").map(PathBuf::from),
            record_value(&provisioned, "MACS3_SRC"),
            record_value(&provisioned, "MACS3_PATH").map(|p| checkout_of_init(&p)),
            Some(root.join(".oracle").join("macs3-src")),
        ]
        .into_iter()
        .flatten(),
        |candidate| candidate.join("MACS3").is_dir(),
    )?;
    let venv = first_existing(
        [
            std::env::var_os("MACS3_VENV").map(PathBuf::from),
            record_value(&provisioned, "MACS3_VENV"),
            Some(root.join(".oracle").join("venv")),
        ]
        .into_iter()
        .flatten(),
        |candidate| candidate.join("bin").join("python").is_file(),
    );
    let bin = first_existing(
        [
            std::env::var_os("MACS3_ORACLE_BIN").map(PathBuf::from),
            venv.as_ref().map(|v| v.join("bin").join("macs3")),
            Some(src.join("bin").join("macs3")),
        ]
        .into_iter()
        .flatten(),
        |candidate| candidate.is_file(),
    );
    let python = first_existing(
        [
            std::env::var_os("MACS3_ORACLE_PYTHON").map(PathBuf::from),
            venv.map(|v| v.join("bin").join("python")),
        ]
        .into_iter()
        .flatten(),
        |candidate| candidate.is_file(),
    );
    Some(Oracle { src, bin, python })
}

/// The reference for a test that cannot do anything without it.
///
/// Prints the skip line and returns `None` when there is no reference installed,
/// so the caller reads as `let Some(oracle) = require("...") else { return };`.
pub fn require(test: &str) -> Option<Oracle> {
    match oracle() {
        Some(oracle) => Some(oracle),
        None => {
            skip(test);
            None
        }
    }
}

/// The oracle's entry point and the interpreter that can import it, for a test that
/// runs the oracle as a subprocess.
pub fn entry_point<'a>(test: &str, oracle: &'a Oracle) -> Option<(&'a Path, &'a Path)> {
    match (&oracle.bin, &oracle.python) {
        (Some(bin), Some(python)) => Some((bin.as_path(), python.as_path())),
        _ => {
            skip(test);
            None
        }
    }
}

/// Say, on stderr, that a differential did not run and why.
pub fn skip(test: &str) {
    eprintln!(
        "skipped: no MACS3 oracle for {test}; provision one with \
         `bash oracle/provision_oracle.sh`, or point MACS3_SRC/MACS3_ORACLE_BIN at an \
         existing one"
    );
}
