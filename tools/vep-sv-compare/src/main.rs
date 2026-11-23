// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `vep-sv-compare` binary: the command line of the structural-variant comparator.

use std::process::ExitCode;

fn main() -> ExitCode {
    vep_sv_compare::cli::main()
}
