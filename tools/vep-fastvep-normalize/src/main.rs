// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Rewrites a fastVEP tab output so that its intergenic rows carry the
//! Feature_type Ensembl VEP writes.
//!
//! fastVEP through v0.4.0 writes `Transcript` in the Feature_type cell of a row whose Gene
//! and Feature cells are `-`, where Ensembl VEP writes `-`; from v0.5.0 it writes `-` there
//! itself, and no row of such an output is rewritten. Feature_type is part of the key the
//! tuple comparison pairs rows on, so a fastVEP output is rewritten before it is scored. The
//! command line is two file names, the output to read and the file to write; the count of
//! rows rewritten goes to stderr.

mod rewrite;
mod stream;
mod threads;

use std::ffi::{OsStr, OsString};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use std::process::ExitCode;

const USAGE: &str = "usage: vep-fastvep-normalize <in> <out>";

fn named(path: &OsStr, detail: impl std::fmt::Display) -> String {
    format!("{}: {detail}", Path::new(path).display())
}

fn run(input: &OsStr, output: &OsStr) -> Result<rewrite::Counts, String> {
    let threads = threads::thread_count()?;
    // The input is opened before the output is made: a run whose input cannot be read
    // leaves no output file, and a run whose two names are one file reads it empty.
    let source = File::open(input).map_err(|error| named(input, error))?;
    let is_directory = source
        .metadata()
        .map_err(|error| named(input, error))?
        .is_dir();
    if is_directory {
        return Err(named(input, "is a directory"));
    }
    let sink = File::create(output).map_err(|error| named(output, error))?;
    stream::normalise(source, sink, threads, stream::RUN_BYTES).map_err(|failure| match failure {
        stream::Failure::Read(error) => named(input, error),
        stream::Failure::Write(error) => named(output, error),
        stream::Failure::NotUtf8(offset) => named(
            input,
            format_args!("the byte at offset {offset} is not UTF-8"),
        ),
    })
}

fn main() -> ExitCode {
    let arguments: Vec<OsString> = std::env::args_os().skip(1).collect();
    // Not locked: a thread that panics writes to stderr while this one waits on it.
    let mut stderr = io::stderr();
    let [input, output] = arguments.as_slice() else {
        let _ = writeln!(stderr, "{USAGE}");
        return ExitCode::FAILURE;
    };
    match run(input, output) {
        Ok(counts) => {
            let _ = writeln!(
                stderr,
                "normalize: rewrote {}/{} intergenic Feature_type=Transcript \u{2192} -",
                counts.rewritten, counts.rows
            );
            ExitCode::SUCCESS
        }
        Err(detail) => {
            let _ = writeln!(stderr, "ERROR: [vep-fastvep-normalize] {detail}");
            ExitCode::FAILURE
        }
    }
}
