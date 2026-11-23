// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep-storable-dump <cache_dir> <ids.txt> <out.txt>`: the `Data::Dumper` text of named
//! transcripts of an Ensembl VEP Storable cache.

use std::io::Write;
use std::os::unix::ffi::OsStrExt;
use std::process::ExitCode;

use vep_storable_dump::run::{run, Request, DIE};

/// The environment variable that sets the number of worker threads.
const THREADS_ENV: &str = "VEP_STORABLE_DUMP_THREADS";

fn threads() -> Result<usize, String> {
    match std::env::var_os(THREADS_ENV) {
        None => Ok(std::thread::available_parallelism().map_or(1, usize::from)),
        Some(value) => value
            .to_str()
            .and_then(|text| text.parse::<usize>().ok())
            .filter(|count| *count >= 1)
            .ok_or_else(|| {
                format!(
                    "{THREADS_ENV}={}: expected a whole number of threads, 1 or more",
                    value.to_string_lossy()
                )
            }),
    }
}

fn main() -> ExitCode {
    let mut arguments = std::env::args_os();
    let program = arguments.next().unwrap_or_default();
    let arguments: Vec<Vec<u8>> = arguments
        .map(|argument| argument.as_bytes().to_vec())
        .collect();
    let mut stderr = std::io::stderr();
    // Three arguments are needed; what follows them is not read.
    let (Some(cache_dir), Some(ids), Some(out)) =
        (arguments.first(), arguments.get(1), arguments.get(2))
    else {
        let usage = [
            b"usage: ",
            program.as_bytes(),
            b" <cache_dir> <ids.txt> <out.txt>\n",
        ]
        .concat();
        stderr.write_all(&usage).ok();
        return ExitCode::from(DIE as u8);
    };
    let threads = match threads() {
        Ok(threads) => threads,
        Err(message) => {
            eprintln!("ERROR: [vep-storable-dump] {message}");
            return ExitCode::from(2);
        }
    };
    let request = Request {
        cache_dir: cache_dir.clone(),
        ids: ids.clone(),
        out: out.clone(),
        threads,
    };
    ExitCode::from(run(&request, &mut stderr) as u8)
}
