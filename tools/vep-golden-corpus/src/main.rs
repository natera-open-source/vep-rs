// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `vep-golden-corpus` binary: the command line goes to the library and the run's exit
//! code comes back.

fn main() {
    let arguments = std::env::args_os().skip(1).collect();
    std::process::exit(vep_golden_corpus::run(arguments));
}
