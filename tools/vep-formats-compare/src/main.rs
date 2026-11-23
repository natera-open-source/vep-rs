// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use std::ffi::OsString;
use std::process::ExitCode;

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().skip(1).collect();
    let code = vep_formats_compare::handover::run(&argv, vep_formats_compare::cli::main);
    ExitCode::from(u8::try_from(code).unwrap_or(1))
}
