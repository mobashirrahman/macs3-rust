//! `macs3-rs` entry point.
//!
//! Dispatch is deliberately thin: the subcommand surface is parsed by
//! [`macs_cli::flags::parse`], which is driven by `oracle/flag_matrix.tsv`, so the
//! CLI's flag surface is *derived* from upstream argparse rather than
//! hand-maintained. Exit codes follow argparse and upstream:
//!
//! * `0` -- success
//! * `2` -- a usage error (unknown flag, missing required flag, bad choice)
//! * `1` -- everything else

use std::process::ExitCode;

use macs_cli::{is_subcommand, parse_flags, PROGRAM, SUBCOMMANDS, TARGET_VERSION};

/// Tune glibc's allocator, then hand off to the real entry.
///
/// The pipeline is chromosome-parallel, so every rayon worker gets its own malloc
/// arena and freed blocks are retained per-arena instead of being returned to the
/// OS. That is a *fragmentation* cost, not a working-set one, and it is large:
/// measured on the 5 M-read fixture (SE `--SPMR`, release, 24 chromosomes), peak
/// RSS is 237 MB with the default arena count and 177 MB with one arena, at
/// unchanged wall clock (7.88 s vs 7.81 s).
///
/// Returning the freed blocks matters as much as not fragmenting them. The chunked
/// passes allocate and drop one chromosome's signal (tens of MB) at a time, and
/// glibc's main arena keeps those blocks in its free list rather than trimming the
/// heap, so RSS creeps up across the pass even though live memory is bounded to a
/// window: with one arena the SE `--SPMR` peak still reached 170 MB. Trimming on
/// every free (`MALLOC_TRIM_THRESHOLD_=0`) drops that to 142 MB -- **0.47x** of
/// upstream, inside the 50% target -- for about 4% wall clock (8.3 s -> 8.7 s,
/// still 3.6x). The same setting helps every command with a chunked pass.
///
/// The settings cannot be applied in-process: glibc reads them through its
/// tunables machinery during `ptmalloc_init`, which runs before `main` -- the Rust
/// runtime itself allocates first, so a `set_var` here is always too late
/// (verified: it left RSS at 237 MB). Re-executing the same binary with the
/// variables set is the standard way around that, and `CommandExt::exec` *replaces*
/// the process image rather than forking, so there is no extra process, no
/// double-wait, and the exit status, argv and stdio all pass through unchanged.
///
/// A user-supplied value wins for each variable individually, and the presence of
/// the private sentinel on the re-exec'd image is what terminates the recursion.
#[cfg(unix)]
fn cap_malloc_arenas() {
    use std::os::unix::process::CommandExt as _;
    if std::env::var_os("MACS3_RS_ALLOC_TUNED").is_some() {
        return;
    }
    let Ok(exe) = std::env::current_exe() else {
        return;
    };
    let mut child = std::process::Command::new(exe);
    child.args(std::env::args_os().skip(1));
    child.env("MACS3_RS_ALLOC_TUNED", "1");
    for (name, value) in [("MALLOC_ARENA_MAX", "1"), ("MALLOC_TRIM_THRESHOLD_", "0")] {
        if std::env::var_os(name).is_none() {
            child.env(name, value);
        }
    }
    // `exec` only returns on failure; fall through and run normally in that case
    // rather than refusing to start.
    let _ = child.exec();
}

#[cfg(not(unix))]
fn cap_malloc_arenas() {}

fn main() -> ExitCode {
    cap_malloc_arenas();
    let argv: Vec<String> = std::env::args().skip(1).collect();

    if argv.is_empty() {
        eprintln!("usage: {PROGRAM} [-h] [--version] <command> [<args>]");
        eprintln!();
        eprintln!("commands: {}", SUBCOMMANDS.join(", "));
        return ExitCode::from(2);
    }

    let head = argv[0].as_str();
    match head {
        "-h" | "--help" => {
            println!("usage: {PROGRAM} [-h] [--version] <command> [<args>]");
            println!();
            println!("commands: {}", SUBCOMMANDS.join(", "));
            println!();
            println!("compatibility target: macs3 {TARGET_VERSION}");
            ExitCode::SUCCESS
        }
        "--version" => {
            println!("{PROGRAM} {TARGET_VERSION}");
            ExitCode::SUCCESS
        }
        _ if !is_subcommand(head) => {
            eprintln!("{PROGRAM}: error: argument command: invalid choice: {head}");
            eprintln!("(choose from {})", SUBCOMMANDS.join(", "));
            ExitCode::from(2)
        }
        _ => {
            match parse_flags(head, &argv[1..]) {
                Err(e) => {
                    eprintln!("usage: {PROGRAM} {head} [-h] ...");
                    eprintln!("{PROGRAM} {head}: {e}");
                    ExitCode::from(2)
                }
                Ok(o) => {
                    if o.help {
                        match macs_cli::help::subcommand_help(head) {
                            Some(text) => print!("{text}"),
                            // Unreachable: `head` is a validated subcommand and every
                            // subcommand has captured help. The stub stays as a
                            // failsafe so a missing file can never become a panic.
                            None => println!("usage: {PROGRAM} {head} [options]"),
                        }
                        return ExitCode::SUCCESS;
                    }
                    // Dispatch. Each command consumes its own validated
                    // `Options` (flag surface derived from the argparse matrix)
                    // and calls the shared library pipeline.
                    let result = match head {
                        "callpeak" => macs_cli::commands::callpeak::run(&o),
                        "bdgopt" => macs_cli::commands::bedgraph_cmds::bdgopt(&o),
                        "cmbreps" => macs_cli::commands::bedgraph_cmds::cmbreps(&o),
                        "filterdup" => macs_cli::commands::filterdup::filterdup(&o),
                        "randsample" => macs_cli::commands::randsample::randsample(&o),
                        "pileup" => macs_cli::commands::pileup::pileup(&o),
                        "bdgpeakcall" => macs_cli::commands::bdgpeakcall::bdgpeakcall(&o),
                        "bdgbroadcall" => macs_cli::commands::bdgpeakcall::bdgbroadcall(&o),
                        "refinepeak" => macs_cli::commands::refinepeak::refinepeak(&o),
                        "predictd" => macs_cli::commands::predictd::predictd(&o),
                        "bdgcmp" => macs_cli::commands::bdgcmp::bdgcmp(&o),
                        "bdgdiff" => macs_cli::commands::bdgcmp::bdgdiff(&o),
                        "hmmratac" => macs_cli::commands::hmmratac::hmmratac(&o),
                        "callvar" => macs_cli::commands::callvar::callvar(&o),
                        _ => {
                            eprintln!("{PROGRAM} {head}: not yet implemented");
                            return ExitCode::FAILURE;
                        }
                    };
                    match result {
                        Ok(()) => ExitCode::SUCCESS,
                        Err(e) => {
                            eprintln!("{PROGRAM} {head}: {e}");
                            ExitCode::FAILURE
                        }
                    }
                }
            }
        }
    }
}
