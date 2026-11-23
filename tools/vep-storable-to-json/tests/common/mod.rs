// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the binary share: where it and the fixtures are, and a directory
//! of their own to run in.

#![allow(dead_code)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub const BINARY: &str = env!("CARGO_BIN_EXE_vep-storable-to-json");

pub fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

/// A fresh directory under the system's temporary one, removed when dropped.
pub struct Scratch(pub PathBuf);

impl Scratch {
    pub fn new(name: &str) -> Self {
        let path = std::env::temp_dir().join(format!(
            "vep-storable-to-json-{}-{}",
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

/// Every file under `root`, by its path relative to it; empty for a root that is absent.
pub fn tree(root: &Path) -> BTreeMap<String, Vec<u8>> {
    fn walk(root: &Path, directory: &Path, files: &mut BTreeMap<String, Vec<u8>>) {
        let Ok(entries) = std::fs::read_dir(directory) else {
            return;
        };
        for entry in entries {
            let path = entry.expect("directory entry").path();
            if path.is_dir() {
                walk(root, &path, files);
            } else {
                let relative = path.strip_prefix(root).expect("under the root");
                files.insert(
                    relative.to_string_lossy().into_owned(),
                    std::fs::read(&path).expect("read file"),
                );
            }
        }
    }
    let mut files = BTreeMap::new();
    walk(root, root, &mut files);
    files
}

/// Runs the binary in `directory` on four worker threads.
pub fn run_in(directory: &Path, arguments: &[&str]) -> Output {
    Command::new(BINARY)
        .args(arguments)
        .current_dir(directory)
        .env("VEP_STORABLE_TO_JSON_THREADS", "4")
        .output()
        .expect("run the binary")
}
