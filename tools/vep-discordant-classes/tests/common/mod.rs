// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the test files share: a scratch directory removed on drop, the bytes of a fixture
//! or a golden, and a run of the binary.

#![allow(dead_code)]

use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

pub const BIN: &str = env!("CARGO_BIN_EXE_vep-discordant-classes");
pub const ROOT: &str = env!("CARGO_MANIFEST_DIR");
pub const THREADS_ENV: &str = "VEP_DISCORDANT_CLASSES_THREADS";

/// The thread bounds every golden is held under; `None` is the host's cores.
pub const BOUNDS: [Option<&str>; 5] = [Some("1"), Some("2"), Some("3"), Some("8"), None];

pub struct Scratch(pub PathBuf);

impl Scratch {
    /// A directory of its own for one test; `tag` is unique within a test file.
    pub fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-discordant-classes-test-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    pub fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
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

/// The bytes of a fixture or golden; empty where the file does not exist, as the golden of
/// a run that writes nothing.
pub fn bytes(name: &str) -> Vec<u8> {
    fs::read(fixture(name)).unwrap_or_default()
}

/// Runs the binary from `dir` with exactly `args` under a thread bound.
pub fn run_bound<S: AsRef<OsStr>>(dir: &Path, args: &[S], bound: Option<&str>) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(dir).args(args).stdin(Stdio::null());
    match bound {
        Some(bound) => command.env(THREADS_ENV, bound),
        None => command.env_remove(THREADS_ENV),
    };
    command.output().unwrap()
}

/// Runs the binary from `dir` with exactly `args` on the host's cores.
pub fn run_in<S: AsRef<OsStr>>(dir: &Path, args: &[S]) -> Output {
    run_bound(dir, args, None)
}

pub fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// Holds a run to exit code 0 and an empty stderr.
pub fn assert_clean_exit(output: &Output) {
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(output.stderr.is_empty(), "{}", text(&output.stderr));
}

/// Holds a run to exit code 1 and the one error line that carries `detail`.
pub fn assert_refused(output: &Output, detail: &str) {
    assert_eq!(output.status.code(), Some(1), "{}", text(&output.stderr));
    assert_eq!(
        text(&output.stderr),
        format!("ERROR: [vep-discordant-classes] {detail}\n")
    );
}

/// Holds a run to exit code 2, the usage line and the error line that carries `message`,
/// with nothing on stdout.
pub fn assert_usage_error(output: &Output, message: &str) {
    assert_eq!(output.status.code(), Some(2), "{}", text(&output.stderr));
    assert_eq!(
        text(&output.stderr),
        format!(
            "usage: vep-discordant-classes [-h] [-o OUTPUT] input\n\
             vep-discordant-classes: error: {message}\n"
        )
    );
    assert!(output.stdout.is_empty());
}
