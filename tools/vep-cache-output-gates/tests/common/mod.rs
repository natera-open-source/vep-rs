// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the test files share: a scratch directory removed on drop, the path and the bytes
//! of a fixture, and a run of the binary from a directory.

#![allow(dead_code)]

use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub const BIN: &str = env!("CARGO_BIN_EXE_vep-cache-output-gates");
pub const ROOT: &str = env!("CARGO_MANIFEST_DIR");

/// The five thresholds `compare_cache_outputs.sh` passes when none is given to it.
pub const DEFAULTS: [&str; 5] = ["0.999", "0.001", "0.01", "0.01", "0.01"];

pub struct Scratch(pub PathBuf);

impl Scratch {
    /// A directory of its own for one test; `tag` is unique within a test file.
    pub fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-cache-output-gates-test-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
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

/// A report of `tests/fixtures/reports/`.
pub fn report(name: &str) -> PathBuf {
    fixture("reports").join(name)
}

pub fn golden(name: &str) -> Vec<u8> {
    fs::read(fixture("expected").join(name)).unwrap()
}

/// Runs the binary in `directory` with exactly `args`.
pub fn run_in<S: AsRef<OsStr>>(directory: &Path, args: &[S]) -> Output {
    Command::new(BIN)
        .current_dir(directory)
        .args(args)
        .output()
        .unwrap()
}

/// The one line a run that reaches no verdict writes on standard error.
pub fn error_line(output: &Output) -> String {
    let stderr = String::from_utf8(output.stderr.clone()).unwrap();
    assert!(
        stderr.starts_with("ERROR: [vep-cache-output-gates] "),
        "{stderr}"
    );
    assert!(
        stderr.ends_with('\n') && stderr.matches('\n').count() == 1,
        "{stderr}"
    );
    stderr
}
