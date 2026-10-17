// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The exit-status contract of the `vep` binary: 0 on success, 1 with one
//! `ERROR: [<component>] <message>` line on a refusal or a failed run, 2 with
//! clap's `error:` on an argument syntax error.

use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn vep() -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
    cmd.stdin(Stdio::null());
    cmd
}

fn corpus_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-release116")
}

fn corpus_args(out: &Path) -> Vec<String> {
    let dir = corpus_dir();
    vec![
        "-i".to_string(),
        dir.join("variants.vcf").display().to_string(),
        "-o".to_string(),
        out.display().to_string(),
        "--offline".to_string(),
        "--json_cache".to_string(),
        dir.join("json_cache").display().to_string(),
        "--assembly".to_string(),
        "GRCh38".to_string(),
        "--force_overwrite".to_string(),
        "--no_stats".to_string(),
        "--max_sv_size".to_string(),
        "-1".to_string(),
        "--quiet".to_string(),
    ]
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

#[test]
fn success_is_0() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = vep().args(corpus_args(&out)).output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
    assert!(out.is_file());
    assert!(!stderr(&output).contains("ERROR:"));
}

/// A run that fails after parsing exits 1 with one contract line: here the
/// input file cannot be opened.
#[test]
fn run_error_is_1() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let mut args = corpus_args(&out);
    args.extend(["-i".to_string(), "/nonexistent/input.vcf".to_string()]);
    let output = vep().args(args).output().unwrap();
    assert_eq!(output.status.code(), Some(1));
    let err = stderr(&output);
    let lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(lines.len(), 1, "{err}");
    assert!(
        lines[0].starts_with("ERROR: [run] Failed to open input file '/nonexistent/input.vcf'"),
        "{err}"
    );
    assert!(!out.exists(), "no output is written for a failed run");
}

/// An argument syntax error is clap's: `error:` prefix, exit 2.
#[test]
fn usage_error_is_2() {
    for argv in [
        vec!["--no_such_flag"],
        vec!["--fork"],
        vec!["--output_format", "xml"],
        vec!["--buffer_size", "many"],
    ] {
        let output = vep().args(&argv).output().unwrap();
        assert_eq!(output.status.code(), Some(2), "{argv:?}");
        let err = stderr(&output);
        assert!(err.starts_with("error: "), "{argv:?}: {err}");
        assert!(!err.contains("ERROR: ["), "{argv:?}: {err}");
    }
}

/// A flag the error table refuses exits 1 with its one line, before the input is
/// read; the clap usage path (exit 2) is not taken.
#[test]
fn refused_flag_is_1() {
    let output = vep()
        .args(["--database", "-i", "/nonexistent/input.vcf", "--quiet"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        stderr(&output),
        "ERROR: [cli] --database is not supported: vep-rs reads a local cache only and never connects to an Ensembl database; name the cache with --json_cache, or --cache with --dir_cache\n"
    );
    assert!(output.stdout.is_empty());
}
