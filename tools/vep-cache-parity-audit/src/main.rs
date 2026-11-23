// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep-cache-parity-audit <storable_cache_dir> <json_cache_dir> [N] [SEED]`: audits the
//! conversion of an Ensembl VEP Storable cache into the JSON cache vep-rs annotates from.

use std::io::Write;
use std::os::unix::ffi::OsStrExt;
use std::process::ExitCode;

use vep_cache_parity_audit::audit::{run, Request, DIE};

/// The environment variable that sets the number of worker threads.
const THREADS_ENV: &str = "VEP_CACHE_PARITY_AUDIT_THREADS";

/// A Perl string is false when it is empty or `0`; such an argument counts as missing.
fn is_true(argument: &[u8]) -> bool {
    !(argument.is_empty() || argument == b"0")
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
    let (Some(storable_root), Some(json_root)) = (
        arguments.first().filter(|argument| is_true(argument)),
        arguments.get(1).filter(|argument| is_true(argument)),
    ) else {
        let usage = [
            b"Usage: ",
            program.as_bytes(),
            b" storable_cache_dir json_cache_dir [N]\n",
        ]
        .concat();
        stderr.write_all(&usage).ok();
        return ExitCode::from(DIE as u8);
    };
    let threads = match threads() {
        Ok(threads) => threads,
        Err(message) => {
            eprintln!("ERROR: [vep-cache-parity-audit] {message}");
            return ExitCode::from(2);
        }
    };
    let request = Request {
        storable_root: storable_root.clone(),
        json_root: json_root.clone(),
        sample: arguments.get(2).cloned(),
        seed: arguments.get(3).cloned(),
        threads,
    };
    let mut stdout = std::io::stdout().lock();
    ExitCode::from(run(&request, &mut stdout, &mut stderr) as u8)
}
