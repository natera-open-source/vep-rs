// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Entry point for the VEP Rust CLI binary.
//!
//! Exit status: 0 on success, and when the command line asks for help or the
//! version or names no argument at all (the usage, as `vep` prints it, vep:29);
//! 1 on a refusal or a failed run, after one line `ERROR: [<component>]
//! <message>` on standard error; 2 on an argument syntax error, in clap's words.

use std::ffi::OsString;

use vep_cli::args;
use vep_cli::config::{self, CliError, ParseFailure};

/// The annotation path allocates and frees per (variant x transcript
/// consequence), and glibc malloc dominates its profile; mimalloc's
/// thread-local heaps take that cost off the rayon workers.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    std::process::exit(run(std::env::args_os().collect()));
}

fn run(argv: Vec<OsString>) -> i32 {
    if argv.len() < 2 {
        let mut command = args::Args::command_with_error_table();
        if command.print_long_help().is_err() {
            return 1;
        }
        return 0;
    }

    let parsed = match config::parse_command_line(argv) {
        Ok(parsed) => parsed,
        Err(ParseFailure::Usage(error)) => error.exit(),
        Err(ParseFailure::Refused(error)) => {
            eprintln!("{error}");
            return 1;
        }
    };

    // vep:32 warns before anything else and whatever --quiet says.
    if let [_, .., last] = parsed.args.output_file.as_slice() {
        eprintln!(
            "WARNING: More than one output file command-line flag provided. Using the last value provided: '{last}'"
        );
    }

    let level = if parsed.args.verbose {
        tracing::Level::DEBUG
    } else if parsed.args.quiet {
        tracing::Level::ERROR
    } else {
        tracing::Level::INFO
    };
    let quiet = parsed.args.quiet;

    let subscriber = tracing_subscriber::fmt()
        .with_max_level(level)
        .with_writer(std::io::stderr)
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("failed to set tracing subscriber");

    let config = match config::Config::from_parsed(parsed) {
        Ok(config) => config,
        Err(error) => {
            eprintln!("{error}");
            return 1;
        }
    };

    if !quiet {
        eprintln!("-------------------- -------------------- --------------------");
        eprintln!(
            "Ensembl VEP (Rust)   v{}.{}",
            vep_core::VEP_VERSION,
            vep_core::VEP_SUB_VERSION,
        );
        eprintln!("-------------------- -------------------- --------------------");
    }

    match vep_cli::runner::Runner::new(config).and_then(|runner| runner.run()) {
        Ok(()) => 0,
        Err(error) => {
            eprintln!("{}", CliError::from_run(&error));
            1
        }
    }
}
