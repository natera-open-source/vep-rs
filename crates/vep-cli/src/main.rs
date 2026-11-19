// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Entry point for the VEP Rust CLI binary.
//!
//! Exit status: 0 on success, and when the command line asks for help or the
//! version or names no argument at all (the usage, as `vep` prints it, vep:29);
//! 1 on a refusal or a failed run, after one line `ERROR: [<component>]
//! <message>` on standard error; 2 on an argument syntax error, in clap's words.
//! A line standard error refuses is dropped ([`vep_cli::stderr`]), and the exit
//! status is the run's either way.
//!
//! The command line is read through [`argv::prepass`] first, which adds the
//! flags of the `--config` file, `<dir>/vep.ini` and the `VEP_*` variables
//! behind the command line's own and spells every option as the parser
//! declares it. `vep filter ...` and `vep recode ...` are the companion tools
//! ([`vep_cli::filter_cli`], [`vep_cli::recode_cli`]), dispatched on their
//! first word before any of that: each reads its own options, none of the
//! sources above.
//!
//! The process ends as soon as its outputs are written and closed: the loaded
//! cache is left to the operating system, which reclaims the heap in one step
//! ([`Teardown::LeaveToExit`]), and the exit skips the runtime's and the
//! allocator's own teardown. An allocator asked to report at exit
//! ([`allocator_reports_at_exit`]) gets the full teardown instead, so the
//! report prints.

use std::ffi::OsString;
use std::io::Write;

use tracing_subscriber::filter::filter_fn;
use tracing_subscriber::prelude::*;

use vep_cli::config::{self, CliError, ParseFailure};
use vep_cli::runner::Teardown;
use vep_cli::{args, argv, stderr, warnings};

/// The annotation path allocates and frees per (variant x transcript
/// consequence), and glibc malloc dominates its profile; mimalloc's
/// thread-local heaps take that cost off the rayon workers.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    stderr::drop_refused_writes();
    let teardown = if allocator_reports_at_exit() {
        Teardown::Release
    } else {
        Teardown::LeaveToExit
    };
    let code = run(std::env::args_os().collect(), teardown);
    end_process(code, teardown)
}

/// Whether the allocator prints a report when the process ends, which the
/// immediate exit of [`end_process`] would skip: mimalloc's `show_stats` and
/// `verbose` options, the second and third members of its `mi_option_e`, set
/// by `MIMALLOC_SHOW_STATS`, `MIMALLOC_VERBOSE` or a build's own defaults.
#[cfg(feature = "mimalloc")]
fn allocator_reports_at_exit() -> bool {
    const MI_OPTION_SHOW_STATS: std::ffi::c_int = 1;
    const MI_OPTION_VERBOSE: std::ffi::c_int = 2;
    extern "C" {
        fn mi_option_is_enabled(option: std::ffi::c_int) -> bool;
    }
    // SAFETY: a read of two options of the allocator this binary links; the
    // function takes a plain integer and touches nothing of ours.
    unsafe { mi_option_is_enabled(MI_OPTION_SHOW_STATS) || mi_option_is_enabled(MI_OPTION_VERBOSE) }
}

#[cfg(not(feature = "mimalloc"))]
fn allocator_reports_at_exit() -> bool {
    false
}

/// Ends the process with `code` once every output is written and closed.
/// Standard output is flushed first. Under [`Teardown::LeaveToExit`] the exit
/// skips the runtime's and the allocator's teardown, which would walk a heap the
/// kernel is about to reclaim; under [`Teardown::Release`] it is the ordinary
/// one, with the allocator's exit hook and its report.
fn end_process(code: i32, teardown: Teardown) -> ! {
    let _ = std::io::stdout().flush();
    match teardown {
        Teardown::Release => std::process::exit(code),
        Teardown::LeaveToExit => exit_immediately(code),
    }
}

#[cfg(unix)]
fn exit_immediately(code: i32) -> ! {
    extern "C" {
        fn _exit(status: std::ffi::c_int) -> !;
    }
    // SAFETY: `_exit` takes one plain integer, never returns and runs no
    // handler of ours; every output is written and closed before it is reached.
    unsafe { _exit(code) }
}

#[cfg(not(unix))]
fn exit_immediately(code: i32) -> ! {
    std::process::exit(code)
}

fn run(argv: Vec<OsString>, teardown: Teardown) -> i32 {
    if argv.len() < 2 {
        let mut command = args::Args::command_with_error_table();
        if command.print_long_help().is_err() {
            return 1;
        }
        return 0;
    }
    if argv[1] == "filter" {
        return vep_cli::filter_cli::run(argv[2..].to_vec());
    }
    if argv[1] == "recode" {
        return vep_cli::recode_cli::run(argv[2..].to_vec());
    }

    let prepassed = match argv::prepass_with_sources(argv, std::env::vars_os()) {
        Ok(prepassed) => prepassed,
        Err(error) => {
            stderr::write(format_args!("{error}\n"));
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
            stderr::write(format_args!("{error}\n"));
            return 1;
        }
    };

    // vep:32 warns before anything else and whatever --quiet says; its line
    // ends with a space before the newline.
    if let [_, .., last] = parsed.args.output_file.as_slice() {
        stderr::write(format_args!(
            "WARNING: More than one output file command-line flag provided. Using the last value provided: '{last}' \n"
        ));
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
        .with_writer(|| stderr::Stream)
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
            stderr::write(format_args!("{error}\n"));
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
    // inherit the name of whichever thread first touches it; --parse_threads
    // sets its width, the pool's own default (the available parallelism) when
    // the flag is 0.
    if let Err(error) = rayon::ThreadPoolBuilder::new()
        .thread_name(|i| format!("vep-parse-{i}"))
        .num_threads(config.parse_threads.unwrap_or(0))
        .build_global()
    {
        stderr::write(format_args!(
            "{}\n",
            CliError::new("run", error.to_string())
        ));
        return 1;
    }

    if !quiet {
        stderr::write(format_args!(
            "-------------------- -------------------- --------------------\n"
        ));
        stderr::write(format_args!(
            "Ensembl VEP (Rust)   v{}.{}\n",
            vep_core::VEP_VERSION,
            vep_core::VEP_SUB_VERSION,
        ));
        stderr::write(format_args!(
            "-------------------- -------------------- --------------------\n"
        ));
    }

    exit_status(vep_cli::runner::Runner::new(config).and_then(|runner| runner.run_with(teardown)))
}

/// The exit status of a run or a `--resolved_config` probe: 0 when it succeeded,
/// else 1 after its `ERROR: [<component>] <message>` line on standard error.
fn exit_status(result: anyhow::Result<()>) -> i32 {
    match result {
        Ok(()) => 0,
        Err(error) => {
            stderr::write(format_args!("{}\n", CliError::from_run(&error)));
            1
        }
    }
}
