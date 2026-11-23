// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep-storable-to-json cache_dir output_dir [chr1,chr2,...]`: converts an Ensembl VEP
//! Storable cache into the JSON cache vep-rs annotates from, with `sequence_regions.tsv`
//! beside its `info.json`: the length and circular flag of every chromosome that has a
//! transcript.

use std::io::{IsTerminal, Write};
use std::os::unix::ffi::OsStrExt;
use std::process::ExitCode;

use vep_storable_to_json::run::{run, Request};

/// The environment variable that sets the number of worker threads.
const THREADS_ENV: &str = "VEP_STORABLE_TO_JSON_THREADS";

/// A Perl string is false when it is empty or `0`; such an argument counts as missing.
fn is_true(argument: &[u8]) -> bool {
    !(argument.is_empty() || argument == b"0")
}

/// `split(/,/, $list)`: the fields between commas, without the empty ones at the end.
fn split_commas(list: &[u8]) -> Vec<Vec<u8>> {
    let mut fields: Vec<Vec<u8>> = list
        .split(|byte| *byte == b',')
        .map(<[u8]>::to_vec)
        .collect();
    while fields.last().is_some_and(Vec::is_empty) {
        fields.pop();
    }
    fields
}

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
    let (Some(cache_dir), Some(output_dir)) = (
        arguments.first().filter(|argument| is_true(argument)),
        arguments.get(1).filter(|argument| is_true(argument)),
    ) else {
        let usage = [
            b"Usage: ",
            program.as_bytes(),
            b" cache_dir output_dir [chr1,chr2,...]\n",
        ]
        .concat();
        stderr.write_all(&usage).ok();
        return ExitCode::from(USAGE_EXIT);
    };
    let threads = match threads() {
        Ok(threads) => threads,
        Err(message) => {
            eprintln!("ERROR: [storable_to_json] {message}");
            return ExitCode::from(2);
        }
    };
    let request = Request {
        cache_dir: cache_dir.clone(),
        output_dir: output_dir.clone(),
        chromosomes: arguments
            .get(2)
            .filter(|argument| is_true(argument))
            .map(|argument| split_commas(argument))
            .unwrap_or_default(),
        threads,
    };
    // Perl holds standard output in a buffer unless it is a terminal, where it writes
    // each line as it is printed.
    let stdout = std::io::stdout();
    let code = if stdout.is_terminal() {
        let mut out = std::io::LineWriter::new(stdout.lock());
        run(&request, &mut out, &mut stderr)
    } else {
        let mut out = std::io::BufWriter::with_capacity(8192, stdout.lock());
        let code = run(&request, &mut out, &mut stderr);
        out.flush().ok();
        code
    };
    ExitCode::from(code as u8)
}

/// The exit code of a Perl `die` when no error number is set: the usage line's.
const USAGE_EXIT: u8 = 255;
