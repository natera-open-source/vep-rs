// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A scratch directory, of a test or of the questions `classify` puts to the format
//! comparator: created under a parent with a prefix and a suffix of its own, removed with
//! its contents when it is dropped.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static CREATED: AtomicU64 = AtomicU64::new(0);

/// A directory of mode 0700 that is removed with its contents on drop.
#[derive(Debug)]
pub struct TempDir {
    pub path: PathBuf,
}

impl TempDir {
    /// Creates `<parent>/<prefix><suffix>`, the suffix made of the process, the clock and a
    /// count so that no two calls share a directory.
    pub fn create(parent: &Path, prefix: &str) -> std::io::Result<TempDir> {
        use std::os::unix::fs::DirBuilderExt;
        loop {
            let nanos = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |elapsed| elapsed.subsec_nanos());
            let count = CREATED.fetch_add(1, Ordering::Relaxed);
            let suffix = format!("{}_{nanos:09}_{count}", std::process::id());
            let path = parent.join(format!("{prefix}{suffix}"));
            match std::fs::DirBuilder::new().mode(0o700).create(&path) {
                Ok(()) => return Ok(TempDir { path }),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            }
        }
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_scratch_directory_is_its_own_and_leaves_with_its_contents() {
        let parent = std::env::temp_dir();
        let first = TempDir::create(&parent, "golden_scratch_test_").unwrap();
        let second = TempDir::create(&parent, "golden_scratch_test_").unwrap();
        assert_ne!(first.path, second.path);
        std::fs::write(first.path.join("file"), "x").unwrap();
        let kept = first.path.clone();
        drop(first);
        assert!(!kept.exists());
        assert!(second.path.is_dir());
    }
}
