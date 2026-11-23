// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A temporary directory as `tempfile.mkdtemp(prefix=...)` makes one: under the first of
//! `TMPDIR`, `TEMP`, `TMP`, `/tmp`, `/var/tmp` and `/usr/tmp` that takes it, named with the
//! prefix and eight characters that differ from run to run, readable by its owner alone.

use crate::pypath;
use std::io::Read;
use std::os::unix::fs::DirBuilderExt;

const ALPHABET: &[u8] = b"abcdefghijklmnopqrstuvwxyz0123456789_";

fn random_name(prefix: &str) -> String {
    let mut bytes = [0u8; 8];
    let read = std::fs::File::open("/dev/urandom").and_then(|mut f| f.read_exact(&mut bytes));
    if read.is_err() {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos() as u64);
        bytes = (nanos ^ (u64::from(std::process::id()) << 32)).to_le_bytes();
    }
    let suffix: String = bytes
        .iter()
        .map(|b| ALPHABET[usize::from(*b) % ALPHABET.len()] as char)
        .collect();
    format!("{prefix}{suffix}")
}

/// The path of a directory made for this run; the error is one line.
pub fn make_temp_dir(prefix: &str) -> Result<String, String> {
    let mut candidates: Vec<String> = ["TMPDIR", "TEMP", "TMP"]
        .iter()
        .filter_map(|name| std::env::var(name).ok())
        .filter(|directory| !directory.is_empty())
        .collect();
    candidates.extend(
        ["/tmp", "/var/tmp", "/usr/tmp"]
            .iter()
            .map(|d| d.to_string()),
    );
    let mut last_error = "no usable temporary directory".to_string();
    for directory in candidates {
        let directory = pypath::resolve(&directory);
        for _ in 0..100 {
            let path = pypath::join(&directory, &random_name(prefix));
            match std::fs::DirBuilder::new().mode(0o700).create(&path) {
                Ok(()) => return Ok(path),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => {
                    last_error = pypath::os_error(&error, &path);
                    break;
                }
            }
        }
    }
    Err(last_error)
}
