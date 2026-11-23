// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cache mode declares `--json-dir` and `--output-dir` alone. A name the
//! parser does not declare, `--region-size` among them, ends the run with
//! clap's usage error, exit 2, the name in the message, before the input is
//! looked for or the output directory created.

use std::path::Path;
use std::process::{Command, Output};

/// A cache conversion from a flat cache that does not exist under `scratch`
/// into `scratch/out`, with `extra` after the two declared flags.
fn convert(scratch: &Path, extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .env("NO_COLOR", "1")
        .arg("--json-dir")
        .arg(scratch.join("flat"))
        .arg("--output-dir")
        .arg(scratch.join("out"))
        .args(extra)
        .output()
        .expect("vep-cache-converter runs")
}

/// `--region-size` is refused as a name no flag has is: the same exit status,
/// clap's line naming the flag, nothing on standard output and no output
/// directory, whether its value follows a blank or `=`.
#[test]
fn region_size_is_refused_as_any_undeclared_flag_is() {
    let scratch = tempfile::tempdir().unwrap();
    let undeclared = convert(scratch.path(), &["--no-such-flag", "1000000"]);
    assert_eq!(undeclared.status.code(), Some(2));
    assert!(
        String::from_utf8_lossy(&undeclared.stderr)
            .starts_with("error: unexpected argument '--no-such-flag' found\n"),
        "{}",
        String::from_utf8_lossy(&undeclared.stderr)
    );

    for spelling in [
        &["--region-size", "1000000"][..],
        &["--region-size=1000000"][..],
    ] {
        let refused = convert(scratch.path(), spelling);
        let stderr = String::from_utf8_lossy(&refused.stderr);
        assert_eq!(
            refused.status.code(),
            undeclared.status.code(),
            "exit status of {spelling:?}: {stderr}"
        );
        assert!(
            stderr.starts_with("error: unexpected argument '--region-size' found\n"),
            "standard error of {spelling:?}: {stderr}"
        );
        assert!(
            refused.stdout.is_empty(),
            "standard output of {spelling:?}: {}",
            String::from_utf8_lossy(&refused.stdout)
        );
        assert!(
            !scratch.path().join("out").exists(),
            "{spelling:?} created the output directory"
        );
    }
}
