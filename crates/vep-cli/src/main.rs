// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Entry point for the VEP Rust CLI binary.
//!
//! Exit status: 0 on success, and when the command line asks for help or the
//! version or names no argument at all (the usage, as `vep` prints it, vep:29);
//! 1 on a refusal or a failed run, after one line `ERROR: [<component>]
//! <message>` on standard error; 2 on an argument syntax error, in clap's words.
//!
//! The command line is read through [`argv::prepass`] first, which adds the
//! flags of the `--config` file, `<dir>/vep.ini` and the `VEP_*` variables
//! behind the command line's own and spells every option as the parser
//! declares it.

use std::ffi::OsString;

use tracing_subscriber::filter::filter_fn;
use tracing_subscriber::prelude::*;

use vep_cli::config::{self, CliError, ParseFailure};
use vep_cli::{args, argv, warnings};

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

    let prepassed = match argv::prepass_with_sources(argv, std::env::vars_os()) {
        Ok(prepassed) => prepassed,
        Err(error) => {
            eprintln!("{error}");
            return 1;
        }
    };
    let parsed = match config::parse_command_line(prepassed.argv) {
        Ok(parsed) => parsed.with_defined_off(
            prepassed
                .defined_off
                .iter()
                .map(|(flag, value)| (flag.as_str(), value.as_str())),
        ),
        Err(ParseFailure::Usage(error)) => error.exit(),
        Err(ParseFailure::Refused(error)) => {
            eprintln!("{error}");
            return 1;
        }
    };

    // vep:32 warns before anything else and whatever --quiet says; its line
    // ends with a space before the newline.
    if let [_, .., last] = parsed.args.output_file.as_slice() {
        eprintln!(
            "WARNING: More than one output file command-line flag provided. Using the last value provided: '{last}' "
        );
    }

    // Every warning goes to the sink instead of the status layer, whatever the
    // level, so --quiet keeps it in the warnings file (BaseVEP.pm:821 drops the
    // standard-error copy alone). --quiet silences every status message
    // (BaseVEP.pm:782), --verbose beside it included: the pair passes under --safe.
    let status_level = if parsed.args.quiet {
        tracing::Level::ERROR
    } else if parsed.args.verbose {
        tracing::Level::DEBUG
    } else {
        tracing::Level::INFO
    };
    let quiet = parsed.args.quiet;
    let status = tracing_subscriber::fmt::layer()
        .with_writer(std::io::stderr)
        .with_filter(filter_fn(move |metadata| {
            *metadata.level() != tracing::Level::WARN && *metadata.level() <= status_level
        }));
    let warnings_layer = warnings::WarningLayer.with_filter(filter_fn(|metadata| {
        *metadata.level() == tracing::Level::WARN
    }));
    tracing_subscriber::registry()
        .with(status)
        .with(warnings_layer)
        .init();
    // Configured before the resolution, whose own warnings are the run's too. A
    // skip names no line number once --input_data is given, whatever its value
    // (BaseVEP.pm:854 tests `defined`), while the value decides what is read.
    warnings::global().configure(
        warnings::target_for(&config::warning_file_for(&parsed.args)),
        quiet,
        parsed.args.input_data.is_some(),
    );

    let config = match config::Config::from_parsed(parsed) {
        Ok(config) => config,
        Err(error) => {
            eprintln!("{error}");
            return 1;
        }
    };

    // The resolved knobs are a property of the configuration and the input, so
    // they print before any cache is resolved or opened.
    if config.resolved_config {
        let written =
            vep_cli::runner::write_resolved_config(&config, &mut std::io::stdout().lock());
        return exit_status(written);
    }

    // Cache shards parse on rayon's global pool, whose threads would otherwise
    // inherit the name of whichever thread first touches it.
    if let Err(error) = rayon::ThreadPoolBuilder::new()
        .thread_name(|i| format!("vep-parse-{i}"))
        .build_global()
    {
        eprintln!("{}", CliError::new("run", error.to_string()));
        return 1;
    }

    if !quiet {
        eprintln!("-------------------- -------------------- --------------------");
        eprintln!(
            "Ensembl VEP (Rust)   v{}.{}",
            vep_core::VEP_VERSION,
            vep_core::VEP_SUB_VERSION,
        );
        eprintln!("-------------------- -------------------- --------------------");
    }

    exit_status(vep_cli::runner::Runner::new(config).and_then(|runner| runner.run()))
}

/// The exit status of a run or a `--resolved_config` probe: 0 when it succeeded,
/// else 1 after its `ERROR: [<component>] <message>` line on standard error.
fn exit_status(result: anyhow::Result<()>) -> i32 {
    match result {
        Ok(()) => 0,
        Err(error) => {
            eprintln!("{}", CliError::from_run(&error));
            1
        }
    }
}
