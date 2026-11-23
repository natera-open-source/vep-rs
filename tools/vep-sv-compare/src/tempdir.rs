// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A scratch directory removed with its contents when dropped.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT: AtomicU64 = AtomicU64::new(0);

/// A directory of mode 0700 under `parent`, removed with its contents on drop.
pub struct TempDir {
    pub path: PathBuf,
}

impl TempDir {
    pub fn create(parent: &Path, prefix: &str) -> std::io::Result<TempDir> {
        use std::os::unix::fs::DirBuilderExt;
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos() as u64);
        loop {
            let serial = NEXT.fetch_add(1, Ordering::Relaxed);
            let name = format!("{prefix}{:x}_{:x}_{serial:x}", std::process::id(), nanos);
            let path = parent.join(name);
            match std::fs::DirBuilder::new().mode(0o700).create(&path) {
                Ok(()) => return Ok(TempDir { path }),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e),
            }
        }
    }

    /// A scratch directory under the system's temporary directory.
    pub fn new(prefix: &str) -> std::io::Result<TempDir> {
        TempDir::create(&std::env::temp_dir(), prefix)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}
