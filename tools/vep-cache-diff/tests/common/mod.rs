// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the test files share: a scratch directory removed on drop, the fixtures' directory,
//! a run of the binary from that directory under a thread count and a buffer bound, and
//! the text `json.dumps` writes for a path.

#![allow(dead_code)]

use std::ffi::OsStr;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub const BIN: &str = env!("CARGO_BIN_EXE_vep-cache-diff");
pub const ROOT: &str = env!("CARGO_MANIFEST_DIR");
pub const THREADS_ENV: &str = "VEP_CACHE_DIFF_THREADS";
pub const BUFFER_ENV: &str = "VEP_CACHE_DIFF_BUFFER_BYTES";

pub struct Scratch(pub PathBuf);

impl Scratch {
    /// A directory of its own for one test, by its resolved path; `tag` is unique within a
    /// test file.
    pub fn new(tag: &str) -> Scratch {
        let path =
            std::env::temp_dir().join(format!("vep-cache-diff-test-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(fs::canonicalize(path).unwrap())
    }

    pub fn path(&self, relative: &str) -> PathBuf {
        self.0.join(relative)
    }

    /// Writes `bytes` under the scratch directory, making the directories above the file.
    pub fn write(&self, relative: &str, bytes: impl AsRef<[u8]>) -> PathBuf {
        let path = self.path(relative);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(&path, bytes).unwrap();
        path
    }

    /// Writes `text` gzipped under the scratch directory.
    pub fn write_gzipped(&self, relative: &str, text: impl AsRef<[u8]>) -> PathBuf {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(text.as_ref()).unwrap();
        self.write(relative, encoder.finish().unwrap())
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// The fixtures' directory by its resolved path, the path the binary reports.
pub fn fixtures() -> PathBuf {
    fs::canonicalize(Path::new(ROOT).join("tests/fixtures")).unwrap()
}

/// How one run is bounded: the thread count and the bytes of variation records a
/// chromosome holds in memory, each the binary's default where `None`.
#[derive(Clone, Copy, Debug)]
pub struct Bounds {
    pub threads: Option<&'static str>,
    pub buffer: Option<&'static str>,
}

impl Bounds {
    pub const DEFAULT: Bounds = Bounds {
        threads: None,
        buffer: None,
    };

    /// The bounds every golden is held under: the defaults, one thread, two threads with
    /// records written out every few lines, and a scratch file for every record.
    pub const ALL: [Bounds; 4] = [
        Bounds::DEFAULT,
        Bounds {
            threads: Some("1"),
            buffer: None,
        },
        Bounds {
            threads: Some("2"),
            buffer: Some("300"),
        },
        Bounds {
            threads: Some("3"),
            buffer: Some("1"),
        },
    ];
}

/// Runs the binary from the fixtures' directory with exactly `args`.
pub fn run_args<S: AsRef<OsStr>>(args: &[S], bounds: Bounds) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(fixtures()).args(args);
    for (name, value) in [(THREADS_ENV, bounds.threads), (BUFFER_ENV, bounds.buffer)] {
        match value {
            Some(value) => command.env(name, value),
            None => command.env_remove(name),
        };
    }
    command.output().unwrap()
}

/// Runs the binary over two cache directories with `report` as its report.
pub fn run(baseline: &Path, builder: &Path, report: &Path, bounds: Bounds) -> Output {
    let args: [&OsStr; 6] = [
        "--baseline-dir".as_ref(),
        baseline.as_os_str(),
        "--builder-dir".as_ref(),
        builder.as_os_str(),
        "--output-report".as_ref(),
        report.as_os_str(),
    ];
    run_args(&args, bounds)
}

/// The text of `path` inside a JSON string as `json.dumps` writes it.
pub fn in_json(path: &Path) -> String {
    let mut out = String::new();
    for c in path.to_str().expect("a UTF-8 path").chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            ' '..='~' => out.push(c),
            _ => {
                let mut units = [0u16; 2];
                for unit in c.encode_utf16(&mut units) {
                    out.push_str(&format!("\\u{unit:04x}"));
                }
            }
        }
    }
    out
}

pub fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}
