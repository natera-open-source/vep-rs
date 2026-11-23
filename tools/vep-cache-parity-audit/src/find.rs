// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The files under a directory in the order `File::Find::find` reports them: a
//! directory's own entries in the order the system lists them, then each of its
//! subdirectories in that order, whole, before the next.

use std::ffi::OsString;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::PathBuf;

fn os_path(bytes: &[u8]) -> PathBuf {
    PathBuf::from(OsString::from_vec(bytes.to_vec()))
}

/// The text of an operating system error without Rust's ` (os error <n>)`.
pub fn os_text(error: &std::io::Error) -> String {
    let text = error.to_string();
    match text.rfind(" (os error ") {
        Some(at) if text.ends_with(')') => text[..at].to_string(),
        _ => text,
    }
}

/// What the walk found: the names, and the lines `File::Find` warns with.
#[derive(Debug, Default)]
pub struct Found {
    /// `$File::Find::name` of every entry the test took, in the walk's order.
    pub names: Vec<Vec<u8>>,
    pub warnings: Vec<Vec<u8>>,
}

/// Walks `root` and keeps every entry that is a file, through a symbolic link too, and
/// whose own name passes `wanted`. A root that is a file is its own only entry; a root
/// that is a symbolic link to a directory is not walked unless its name ends in `/`.
pub fn find(root: &[u8], wanted: impl Fn(&[u8]) -> bool) -> Found {
    let mut found = Found::default();
    let is_file = |path: &[u8]| std::fs::metadata(os_path(path)).is_ok_and(|meta| meta.is_file());
    let top_meta = match std::fs::symlink_metadata(os_path(root)) {
        Ok(meta) => meta,
        Err(error) => {
            let mut line = b"Can't stat ".to_vec();
            line.extend_from_slice(top_item(root));
            line.extend_from_slice(format!(": {}", os_text(&error)).as_bytes());
            found.warnings.push(line);
            return found;
        }
    };
    let top = top_item(root);
    if !top_meta.is_dir() {
        // `fileparse`: the directory with its slash, `./` where the name has none.
        let split = top
            .iter()
            .rposition(|byte| *byte == b'/')
            .map_or(0, |at| at + 1);
        let name = &top[split..];
        let mut full = if split == 0 {
            b"./".to_vec()
        } else {
            top[..split].to_vec()
        };
        full.extend_from_slice(name);
        if is_file(top) && wanted(name) {
            found.names.push(full);
        }
        return found;
    }
    // Directories still to walk; the last one is walked next.
    let mut stack: Vec<Vec<u8>> = vec![top.to_vec()];
    while let Some(directory) = stack.pop() {
        let listing = match std::fs::read_dir(os_path(&directory)) {
            Ok(listing) => listing,
            Err(error) => {
                let mut line = b"Can't opendir(".to_vec();
                line.extend_from_slice(&directory);
                line.extend_from_slice(format!("): {}", os_text(&error)).as_bytes());
                found.warnings.push(line);
                continue;
            }
        };
        let prefix = if directory == b"/" {
            b"/".to_vec()
        } else {
            [directory.as_slice(), b"/"].concat()
        };
        let mut subdirectories = Vec::new();
        for entry in listing.filter_map(Result::ok) {
            let name = entry.file_name();
            let name = name.as_bytes();
            let full = [prefix.as_slice(), name].concat();
            // A symbolic link is no directory to the walk, whatever it points at.
            if entry.file_type().is_ok_and(|kind| kind.is_dir()) {
                subdirectories.push(full);
            } else if is_file(&full) && wanted(name) {
                found.names.push(full);
            }
        }
        // The first subdirectory listed is walked first.
        stack.extend(subdirectories.into_iter().rev());
    }
    found
}

/// The root without one final slash, unless it is `/` itself.
fn top_item(root: &[u8]) -> &[u8] {
    match root.strip_suffix(b"/") {
        Some(rest) if !rest.is_empty() => rest,
        _ => root,
    }
}

/// A name ends in `suffix`, where Perl's `$` also stands before one final line feed.
pub fn name_ends(name: &[u8], suffix: &[u8]) -> bool {
    name.ends_with(suffix)
        || name
            .strip_suffix(b"\n")
            .is_some_and(|name| name.ends_with(suffix))
}

#[cfg(test)]
#[path = "find_tests.rs"]
mod tests;
