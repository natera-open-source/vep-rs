// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the binary ends a run that met an input no reading of the comparison expects.

mod support;

use std::ffi::OsStr;
use std::os::unix::ffi::OsStrExt;
use support::{run, scratch};

/// The run of a comparison under `floors`, which is read before any file.
fn run_under(scratch_name: &str, floors: &OsStr) -> support::Ran {
    let dir = scratch(scratch_name);
    run(
        &dir,
        &[
            OsStr::new("--format"),
            OsStr::new("default"),
            OsStr::new("--expected"),
            OsStr::new("expected.txt"),
            OsStr::new("--actual"),
            OsStr::new("actual.txt"),
            OsStr::new("--floors"),
            floors,
        ],
    )
}

#[test]
fn an_unexpected_ending_is_one_error_line_exit_code_1_and_no_standard_output() {
    let ran = run_under("handover-plain", OsStr::new("IMPACT=9223372036854775808"));
    assert_eq!(ran.code, 1);
    assert!(ran.stdout.is_empty());
    assert_eq!(
        ran.err(),
        "ERROR: [compare_vep_formats] --floors IMPACT=9223372036854775808: \
         a minimum of 2**63 or more is not read\n"
    );
}

#[test]
fn an_unexpected_ending_spells_an_argument_byte_that_is_not_utf8_as_its_escape() {
    let ran = run_under(
        "handover-held-byte",
        OsStr::from_bytes(b"IM\xe9=9223372036854775808"),
    );
    assert_eq!(ran.code, 1);
    assert!(ran.stdout.is_empty());
    assert_eq!(
        ran.stderr,
        b"ERROR: [compare_vep_formats] --floors IM\\udce9=9223372036854775808: \
          a minimum of 2**63 or more is not read\n"
    );
}
