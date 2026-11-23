// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep-storable-digest-join`: joins two digest files of `vep-storable-digest` on the
//! transcript id.

use std::io::Write;
use std::os::unix::ffi::OsStrExt;
use std::process::ExitCode;

use vep_storable_digest_join::argparse::{parse, Opt, Parsed};
use vep_storable_digest_join::join::{run, Request};
use vep_storable_digest_join::py::path_str;

const PROGRAM: &str = "vep-storable-digest-join";

const OPTIONS: [Opt; 8] = [
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
        name: "--parity-csv",
        required: false,
    },
    Opt {
        name: "--sweep-id",
        required: false,
    },
    Opt {
        name: "--assembly",
        required: false,
    },
];

const USAGE: &str = "usage: vep-storable-digest-join [-h] --version-a VERSION_A --dump-a DUMP_A --version-b VERSION_B --dump-b DUMP_B --out-dir OUT_DIR [--parity-csv PARITY_CSV] [--sweep-id SWEEP_ID] [--assembly ASSEMBLY]";

const HELP: &str = "
Joins two digest files of vep-storable-digest on the transcript id and writes, under
--out-dir, the ids in one cache only, those identical in content alone and those that
differ, with a summary. With --parity-csv and --assembly it also extends that assembly's
row of the parity table with the counts.

options:
  -h, --help            show this help message and exit
  --version-a VERSION_A
  --dump-a DUMP_A
  --version-b VERSION_B
  --dump-b DUMP_B
  --out-dir OUT_DIR
  --parity-csv PARITY_CSV
  --sweep-id SWEEP_ID   id of the run that produced the two dumps (sweep-<UTC>), written to the parity row
  --assembly ASSEMBLY
";

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
            parity_csv: path("--parity-csv"),
            sweep_id: text("--sweep-id")?,
            assembly: text("--assembly")?,
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
