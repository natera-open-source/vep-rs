// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep-storable-transcript-classes`: classes the field differences of two cache
//! builds from the dumps of their differing transcripts.

use std::io::Write;
use std::os::unix::ffi::OsStrExt;
use std::process::ExitCode;

use vep_storable_transcript_classes::argparse::{parse, Opt, Parsed};
use vep_storable_transcript_classes::classify::{run, Request};
use vep_storable_transcript_classes::py::path_str;

const PROGRAM: &str = "vep-storable-transcript-classes";

/// The environment variable that sets the number of worker threads.
const THREADS_ENV: &str = "VEP_STORABLE_TRANSCRIPT_CLASSES_THREADS";

const OPTIONS: [Opt; 7] = [
    Opt {
        name: "--version-a",
        required: true,
    },
    Opt {
        name: "--dump-a",
        required: true,
    },
    Opt {
        name: "--version-b",
        required: true,
    },
    Opt {
        name: "--dump-b",
        required: true,
    },
    Opt {
        name: "--out-dir",
        required: true,
    },
    Opt {
        name: "--json-cache",
        required: false,
    },
    Opt {
        name: "--json-cache-label",
        required: false,
    },
];

const USAGE: &str = "usage: vep-storable-transcript-classes [-h] --version-a VERSION_A --dump-a DUMP_A --version-b VERSION_B --dump-b DUMP_B --out-dir OUT_DIR [--json-cache JSON_CACHE] [--json-cache-label JSON_CACHE_LABEL]";

const HELP: &str = "
Classes how two Ensembl VEP Storable cache builds differ on the transcripts whose content
differs, from two dumps of vep-storable-dump for the same transcript ids: peptide,
predictions, protein_features, other. Writes one row per transcript and a summary under
--out-dir. With --json-cache the stored peptide of each transcript in a vep-rs JSON cache
is matched against both builds.

options:
  -h, --help            show this help message and exit
  --version-a VERSION_A
  --dump-a DUMP_A
  --version-b VERSION_B
  --dump-b DUMP_B
  --out-dir OUT_DIR
  --json-cache JSON_CACHE
  --json-cache-label JSON_CACHE_LABEL
                        name recorded for the JSON cache (default: the directory's name)
";

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
    let argv: Vec<Vec<u8>> = std::env::args_os()
        .skip(1)
        .map(|argument| argument.as_bytes().to_vec())
        .collect();
    let values = match parse(&OPTIONS, &argv) {
        Ok(Parsed::Values(values)) => values,
        Ok(Parsed::Help) => {
            println!("{USAGE}\n{HELP}");
            return ExitCode::SUCCESS;
        }
        Err(message) => {
            eprintln!("{USAGE}\n{PROGRAM}: error: {message}");
            return ExitCode::from(2);
        }
    };
    let text = |name: &str| -> Result<Option<String>, String> {
        values
            .get(name)
            .map(|value| {
                String::from_utf8(value.clone())
                    .map_err(|_| format!("the value of {name} is not UTF-8"))
            })
            .transpose()
    };
    let path = |name: &str| values.get(name).map(|value| path_str(value));
    let request = (|| -> Result<Request, String> {
        Ok(Request {
            version_a: text("--version-a")?.unwrap_or_default(),
            dump_a: path("--dump-a").unwrap_or_default(),
            version_b: text("--version-b")?.unwrap_or_default(),
            dump_b: path("--dump-b").unwrap_or_default(),
            out_dir: path("--out-dir").unwrap_or_default(),
            json_cache: path("--json-cache"),
            json_cache_label: text("--json-cache-label")?,
            threads: threads()?,
        })
    })();
    let request = match request {
        Ok(request) => request,
        Err(detail) => {
            eprintln!("ERROR: [{PROGRAM}] {detail}");
            return ExitCode::from(1);
        }
    };
    let mut stdout = std::io::stdout().lock();
    match run(&request, &mut stdout) {
        Ok(code) => ExitCode::from(code as u8),
        Err(failure) => {
            stdout.flush().ok();
            eprintln!("{}", failure.line);
            ExitCode::from(failure.code as u8)
        }
    }
}
