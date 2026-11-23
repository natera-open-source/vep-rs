// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the binary share: where it and the fixtures are, and a directory
//! of their own to run in.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub const BINARY: &str = env!("CARGO_BIN_EXE_vep-storable-digest");

pub fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

/// A fresh directory under the system's temporary one, removed when dropped.
pub struct Scratch(pub PathBuf);

impl Scratch {
    pub fn new(name: &str) -> Self {
        let path = std::env::temp_dir().join(format!(
            "vep-storable-digest-{}-{}",
            std::process::id(),
            name.replace('/', "_")
        ));
        std::fs::remove_dir_all(&path).ok();
        std::fs::create_dir_all(&path).expect("scratch directory");
        Scratch(path)
    }

    pub fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).ok();
    }
}

pub fn copy_tree(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).expect("create directory");
    for entry in
        std::fs::read_dir(from).unwrap_or_else(|error| panic!("{}: {error}", from.display()))
    {
        let entry = entry.expect("directory entry");
        let target = to.join(entry.file_name());
        if entry.file_type().expect("file type").is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), &target).expect("copy file");
        }
    }
}

/// Runs the binary in `directory` on `threads` worker threads.
pub fn run_in(directory: &Path, arguments: &[&str], threads: &str) -> Output {
    Command::new(BINARY)
        .args(arguments)
        .current_dir(directory)
        .env("VEP_STORABLE_DIGEST_THREADS", threads)
        .output()
        .expect("run the binary")
}

/// A recorded standard error with the program's own path in the place of `<program>`.
pub fn with_program(recorded: Vec<u8>) -> Vec<u8> {
    let marker = b"<program>";
    match recorded
        .windows(marker.len())
        .position(|window| window == marker)
    {
        Some(at) => [
            &recorded[..at],
            BINARY.as_bytes(),
            &recorded[at + marker.len()..],
        ]
        .concat(),
        None => recorded,
    }
}
