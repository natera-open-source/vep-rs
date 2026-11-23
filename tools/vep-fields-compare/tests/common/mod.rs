// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the test files share: a scratch directory removed on drop, the text of a golden, and
//! a run of the binary from the package root.

#![allow(dead_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub const BIN: &str = env!("CARGO_BIN_EXE_vep-fields-compare");
pub const ROOT: &str = env!("CARGO_MANIFEST_DIR");

pub struct Scratch(pub PathBuf);

impl Scratch {
    /// A directory of its own for one test; `tag` is unique within a test file.
    pub fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-fields-compare-test-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    /// Where `run` asks for the report.
    pub fn report(&self) -> PathBuf {
        self.0.join("report").join("fields.json")
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

pub fn fixture(name: &str) -> PathBuf {
    Path::new(ROOT).join("tests/fixtures").join(name)
}

pub fn golden(name: &str) -> String {
    fs::read_to_string(fixture(name)).unwrap()
}

/// Runs the binary from the package root with exactly `args`.
pub fn run_args<S: AsRef<std::ffi::OsStr>>(args: &[S]) -> Output {
    run_in(Path::new(ROOT), args)
}

/// Runs the binary from the package root with exactly `args` and `path` as its `PATH`, so
/// the `sort` it starts is the one under `path` or none.
pub fn run_with_path<S: AsRef<std::ffi::OsStr>>(path: &Path, args: &[S]) -> Output {
    Command::new(BIN)
        .current_dir(ROOT)
        .env("PATH", path)
        .args(args)
        .output()
        .unwrap()
}

/// Runs the binary from `dir` with exactly `args`.
pub fn run_in<S: AsRef<std::ffi::OsStr>>(dir: &Path, args: &[S]) -> Output {
    Command::new(BIN)
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap()
}

/// Runs the binary from the package root with the fixture paths relative to it, so the
/// report's `inputs` block carries the strings the golden holds. Returns the run and the
/// report's text, empty when none was written.
pub fn run(scratch: &Scratch, reference: &str, scored: &str, extra: &[&str]) -> (Output, String) {
    let out = scratch.report();
    let output = Command::new(BIN)
        .current_dir(ROOT)
        .args(["--perl", reference, "--rust", scored, "--out"])
        .arg(&out)
        .arg("--tmp-dir")
        .arg(&scratch.0)
        .args(extra)
        .output()
        .unwrap();
    let json = fs::read_to_string(&out).unwrap_or_default();
    (output, json)
}

pub fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

pub fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

pub fn assert_clean_exit(output: &Output) {
    assert!(
        output.status.success(),
        "exit {:?}: {}",
        output.status.code(),
        stderr(output)
    );
    assert!(output.stderr.is_empty(), "{}", stderr(output));
}
