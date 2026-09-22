// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Entry point for the VEP Rust CLI binary.

use clap::Parser;

use vep_cli::args;
use vep_cli::config;

/// The annotation path allocates and frees per (variant x transcript
/// consequence), and glibc malloc dominates its profile; mimalloc's
/// thread-local heaps take that cost off the rayon workers.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> anyhow::Result<()> {
    let args = args::Args::parse();

    let level = if args.verbose {
        tracing::Level::DEBUG
    } else if args.quiet {
        tracing::Level::ERROR
    } else {
        tracing::Level::INFO
    };

    let subscriber = tracing_subscriber::fmt()
        .with_max_level(level)
        .with_writer(std::io::stderr)
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("failed to set tracing subscriber");

    if !args.quiet {
        eprintln!("-------------------- -------------------- --------------------");
        eprintln!(
            "Ensembl VEP (Rust)   v{}.{}",
            vep_core::VEP_VERSION,
            vep_core::VEP_SUB_VERSION,
        );
        eprintln!("-------------------- -------------------- --------------------");
    }

    let config = config::Config::from_args(args);
    let runner = vep_cli::runner::Runner::new(config);
    runner.run()
}
