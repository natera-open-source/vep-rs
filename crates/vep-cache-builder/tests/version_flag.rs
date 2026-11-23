// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep-cache-builder --version` and `-V` print the binary's name and the
//! crate version on one line of standard output and exit 0, with none of the
//! flags a build requires on the line and with them.

use std::process::Command;

/// Runs the builder with `args` and asserts the version line, exit 0 and an
/// empty standard error.
fn assert_prints_version(args: &[&str]) {
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(args)
        .output()
        .unwrap();
    assert_eq!(
        String::from_utf8(output.stderr).unwrap(),
        "",
        "standard error of {args:?}"
    );
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("vep-cache-builder {}\n", env!("CARGO_PKG_VERSION")),
        "standard output of {args:?}"
    );
    assert_eq!(output.status.code(), Some(0), "exit code of {args:?}");
}

#[test]
fn version_flag_prints_the_name_and_crate_version_and_exits_0() {
    assert_prints_version(&["--version"]);
    assert_prints_version(&["-V"]);
}

#[test]
fn version_flag_wins_over_the_flags_of_a_build() {
    let parent = tempfile::tempdir().unwrap();
    let output_dir = parent.path().join("cache");
    assert_prints_version(&[
        "--release",
        "116",
        "--assembly",
        "GRCh38",
        "--output-dir",
        output_dir.to_str().unwrap(),
        "--version",
    ]);
    assert!(!output_dir.exists());
}
