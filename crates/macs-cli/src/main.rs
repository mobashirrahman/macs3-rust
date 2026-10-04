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

fn main() -> ExitCode {
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
