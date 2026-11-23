// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};

use super::{find, name_ends};

/// A directory of the test's own under the system's temporary one.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let path = std::env::temp_dir().join(format!(
            "vep-cache-parity-audit-find-{}-{name}",
            std::process::id()
        ));
        std::fs::remove_dir_all(&path).ok();
        std::fs::create_dir_all(&path).expect("scratch directory");
        Scratch(path)
    }

    fn file(&self, relative: &str) {
        let path = self.0.join(relative);
        std::fs::create_dir_all(path.parent().expect("a parent")).expect("directories");
        std::fs::write(path, b"x").expect("a file");
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).ok();
    }
}

fn listed(directory: &Path) -> Vec<String> {
    std::fs::read_dir(directory)
        .expect("a listing")
        .map(|entry| {
            entry
                .expect("an entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect()
}

/// What the walk must give for a tree: a directory's files as the system lists them,
/// then each subdirectory as the system lists it, whole.
fn expected(directory: &Path, prefix: &str, wanted: &dyn Fn(&str) -> bool, out: &mut Vec<String>) {
    let names = listed(directory);
    for name in &names {
        if directory.join(name).is_file() && wanted(name) {
            out.push(format!("{prefix}/{name}"));
        }
    }
    for name in &names {
        let path = directory.join(name);
        if path.is_dir() && !path.is_symlink() {
            expected(&path, &format!("{prefix}/{name}"), wanted, out);
        }
    }
}

#[test]
fn files_come_before_subdirectories_and_each_subdirectory_whole() {
    let scratch = Scratch::new("order");
    for file in [
        "top.gz",
        "b/2-3.gz",
        "b/1-2.gz",
        "b/deep/er/9-9.gz",
        "b/deep/0-1.gz",
        "a/5-6.gz",
        "a/1-1_var.gz",
        "a/1-1_reg.gz",
        "a/notes.txt",
        "c/x.gz",
        "empty/.keep",
    ] {
        scratch.file(file);
    }
    std::os::unix::fs::symlink(scratch.0.join("b"), scratch.0.join("link-to-b")).expect("a link");
    std::os::unix::fs::symlink(scratch.0.join("top.gz"), scratch.0.join("a/link.gz"))
        .expect("a link");
    let root = scratch.0.to_string_lossy().into_owned();
    let wanted = |name: &str| {
        name.ends_with(".gz") && !name.ends_with("_var.gz") && !name.ends_with("_reg.gz")
    };
    let found = find(root.as_bytes(), |name| {
        wanted(&String::from_utf8_lossy(name))
    });
    assert!(found.warnings.is_empty());
    let names: Vec<String> = found
        .names
        .iter()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect();
    let mut want = Vec::new();
    expected(&scratch.0, &root, &wanted, &mut want);
    assert_eq!(names, want);
    // A link to a file is a file; a link to a directory is not walked.
    assert!(names.iter().any(|name| name.ends_with("/a/link.gz")));
    assert!(!names.iter().any(|name| name.contains("link-to-b")));
    assert_eq!(names.len(), 8);
    // One final slash of the root is dropped.
    let slashed = find(format!("{root}/").as_bytes(), |name| {
        wanted(&String::from_utf8_lossy(name))
    });
    assert_eq!(slashed.names, found.names);
}

#[test]
fn a_root_that_is_a_file_is_its_own_entry() {
    let scratch = Scratch::new("file-root");
    scratch.file("one.gz");
    let file = scratch.0.join("one.gz");
    let found = find(file.as_os_str().as_bytes(), |name| name == b"one.gz");
    assert_eq!(found.names, vec![file.as_os_str().as_bytes().to_vec()]);
    assert!(find(file.as_os_str().as_bytes(), |_| false)
        .names
        .is_empty());
}

#[test]
fn a_root_that_is_absent_is_warned_of() {
    let found = find(b"/no/such/root/", |_| true);
    assert!(found.names.is_empty());
    assert_eq!(
        found.warnings,
        vec![b"Can't stat /no/such/root: No such file or directory".to_vec()]
    );
}

#[test]
fn a_name_ends_where_perl_dollar_stands() {
    assert!(name_ends(b"1-2.gz", b".gz"));
    assert!(name_ends(b"1-2.gz\n", b".gz"));
    assert!(!name_ends(b"1-2.gz\n\n", b".gz"));
    assert!(!name_ends(b"1-2.gzip", b".gz"));
}
