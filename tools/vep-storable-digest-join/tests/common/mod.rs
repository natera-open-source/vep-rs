// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the binary share: where it and the fixtures are, and a directory
//! of their own to run in.

#![allow(dead_code)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub const BINARY: &str = env!("CARGO_BIN_EXE_vep-storable-digest-join");

pub fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

/// A fresh directory under the system's temporary one, removed when dropped.
pub struct Scratch(pub PathBuf);

impl Scratch {
    pub fn new(name: &str) -> Self {
        let path = std::env::temp_dir().join(format!(
            "vep-storable-digest-join-{}-{}",
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

/// Every file a reference run wrote or changed, by its path in the directory of the run:
/// the files under `tree/` of `expected` and, for a case with `files.tsv`, the files of
/// `files/` under the paths that table spells.
///
/// Windows refuses a file name with a double quote or a backslash, and a file system
/// that normalises names can change one with a character outside ASCII, so a file of
/// such a name lies in `files/` under a plain one. The first line of `files.tsv` states
/// its grammar; each line after it holds a file of `files/`, a tab and the path, and
/// every file of `files/` has one line.
pub fn written(expected: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut files = tree(&expected.join("tree"));
    let mut plain = tree(&expected.join("files"));
    let table = expected.join("files.tsv");
    if table.exists() {
        let text = std::fs::read_to_string(&table).expect("files.tsv");
        let mut rows = text.lines();
        let grammar = rows.next().unwrap_or_default();
        assert!(grammar.starts_with('#'), "files.tsv: no grammar line");
        for row in rows {
            let (file, path) = row
                .split_once('\t')
                .unwrap_or_else(|| panic!("files.tsv: no tab in {row:?}"));
            let bytes = plain
                .remove(file)
                .unwrap_or_else(|| panic!("files.tsv: files/{file} is absent or on two lines"));
            let other = files.insert(unescape(path), bytes);
            assert!(other.is_none(), "files.tsv: {path} is held twice");
        }
    }
    assert!(
        plain.is_empty(),
        "files/ without a line of files.tsv: {:?}",
        plain.keys()
    );
    files
}

/// The path a spelling of `files.tsv` stands for: `\uHHHH`, `H` a hexadecimal digit, is
/// the character of that code point, `\"` a double quote, `\\` a backslash, and any other
/// printable ASCII character itself. A spelling outside this grammar is refused.
fn unescape(spelling: &str) -> String {
    let mut path = String::new();
    let mut characters = spelling.chars();
    while let Some(character) = characters.next() {
        match character {
            '\\' => match characters.next() {
                Some('u') => {
                    let digits: String = characters.by_ref().take(4).collect();
                    let code = u32::from_str_radix(&digits, 16)
                        .ok()
                        .filter(|_| {
                            digits.len() == 4 && digits.bytes().all(|b| b.is_ascii_hexdigit())
                        })
                        .and_then(char::from_u32);
                    path.push(code.unwrap_or_else(|| {
                        panic!("files.tsv: {spelling}: \\u{digits} is no character")
                    }));
                }
                Some(escaped @ ('"' | '\\')) => path.push(escaped),
                other => panic!("files.tsv: {spelling}: a backslash before {other:?}"),
            },
            ' '..='~' if character != '"' => path.push(character),
            other => panic!("files.tsv: {spelling}: {other:?} stands unescaped"),
        }
    }
    path
}

/// Runs the binary in `directory`.
pub fn run_in(directory: &Path, arguments: &[&str]) -> Output {
    Command::new(BINARY)
        .args(arguments)
        .current_dir(directory)
        .output()
        .expect("run the binary")
}
