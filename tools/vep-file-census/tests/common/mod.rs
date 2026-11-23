// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the test files share: a case of a golden test, a run of the binary inside the
//! fixture directory, a scratch directory removed on drop, and the count and sum of some
//! lines added up a byte at a time.

#![allow(dead_code)]

use sha2::{Digest, Sha256};
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

pub const BIN: &str = env!("CARGO_BIN_EXE_vep-file-census");
pub const ROOT: &str = env!("CARGO_MANIFEST_DIR");

/// The thread counts every golden is held under; `None` leaves the host's cores.
pub const THREADS: [Option<&str>; 3] = [Some("1"), Some("2"), None];

/// One golden: the shell pipeline that wrote `expected/<name>.out` and left with
/// `status`, run inside the fixture directory, and the command line that stands for it.
pub struct Case {
    pub name: &'static str,
    pub rule: &'static str,
    pub stdin: Option<&'static str>,
    pub args: &'static [&'static str],
    pub status: i32,
}

pub fn fixtures() -> PathBuf {
    Path::new(ROOT).join("tests/fixtures")
}

pub fn expected(name: &str) -> Vec<u8> {
    let path = fixtures().join("expected").join(name);
    fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// Runs the binary inside the fixture directory; `stdin` names a fixture.
pub fn run<S: AsRef<std::ffi::OsStr>>(
    args: &[S],
    stdin: Option<&str>,
    threads: Option<&str>,
) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(fixtures()).args(args);
    match threads {
        Some(count) => command.env("VEP_FILE_CENSUS_THREADS", count),
        None => command.env_remove("VEP_FILE_CENSUS_THREADS"),
    };
    match stdin {
        Some(name) => command.stdin(File::open(fixtures().join(name)).unwrap()),
        None => command.stdin(Stdio::null()),
    };
    command.output().unwrap()
}

/// Holds the binary to the case's golden and status under every thread count.
pub fn hold(case: &Case) {
    let golden = expected(&format!("{}.out", case.name));
    for threads in THREADS {
        let output = run(case.args, case.stdin, threads);
        assert_eq!(
            output.stdout,
            golden,
            "{} on {threads:?} thread(s): standard output is not that of `{}`; stderr: {}",
            case.name,
            case.rule,
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            output.status.code(),
            Some(case.status),
            "{} on {threads:?} thread(s): the exit status is not that of `{}`",
            case.name,
            case.rule
        );
    }
}

/// The line `line-sum` prints of `lines`: their count, a tab, and the sum of their
/// SHA256 digests modulo 2**256 in 64 hexadecimal digits, added up a byte at a time from
/// the last.
pub fn line_sum<'a>(lines: impl Iterator<Item = &'a [u8]>) -> String {
    let mut sum = [0u8; 32];
    let mut count = 0u64;
    for line in lines {
        let mut carry = 0u16;
        for (byte, more) in sum.iter_mut().zip(Sha256::digest(line)).rev() {
            let total = u16::from(*byte) + u16::from(more) + carry;
            *byte = total.to_le_bytes()[0];
            carry = total >> 8;
        }
        count += 1;
    }
    let hex: String = sum.iter().map(|byte| format!("{byte:02x}")).collect();
    format!("{count}\t{hex}\n")
}

pub struct Scratch(pub PathBuf);

impl Scratch {
    /// A directory of its own for one test; `tag` is unique within a test file.
    pub fn new(tag: &str) -> Scratch {
        let path =
            std::env::temp_dir().join(format!("vep-file-census-test-{}-{tag}", std::process::id()));
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
