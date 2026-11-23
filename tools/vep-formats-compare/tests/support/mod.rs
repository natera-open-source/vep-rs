// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the command-line tests share: the binary, a scratch directory, the golden corpora
//! and a report read back.

#![allow(dead_code)]

use std::ffi::{OsStr, OsString};
use std::io::Read;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;
use vep_formats_compare::pyjson::{self, Value};

/// The variables that would change what a run writes; a run names the ones it sets.
const CLEARED: &[&str] = &[
    "COLUMNS",
    "LINES",
    "NO_COLOR",
    "FORCE_COLOR",
    "PYTHON_COLORS",
    "VEP_FORMATS_COMPARE_THREADS",
];

/// How one run of the binary ended.
#[derive(Debug)]
pub struct Ran {
    pub code: i32,
    pub stdout: Vec<u8>,
    pub stderr: Vec<u8>,
}

impl Ran {
    pub fn out(&self) -> String {
        String::from_utf8_lossy(&self.stdout).into_owned()
    }

    pub fn err(&self) -> String {
        String::from_utf8_lossy(&self.stderr).into_owned()
    }
}

/// The binary on `args`, in the directory `cwd`, with `env` set.
pub fn run_in<A: AsRef<OsStr>>(cwd: &Path, env: &[(OsString, OsString)], args: &[A]) -> Ran {
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-formats-compare"));
    command.current_dir(cwd).args(args);
    for name in CLEARED {
        command.env_remove(name);
    }
    for (name, value) in env {
        command.env(name, value);
    }
    let output = command.output().expect("the binary runs");
    Ran {
        code: output.status.code().expect("the binary exits"),
        stdout: output.stdout,
        stderr: output.stderr,
    }
}

/// The binary on `args` in `cwd`.
pub fn run<A: AsRef<OsStr>>(cwd: &Path, args: &[A]) -> Ran {
    run_in(cwd, &[], args)
}

/// An empty directory of this test's own; what an earlier run left there is removed.
pub fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("vep-formats-compare-{name}"));
    if dir.exists() {
        make_removable(&dir);
        std::fs::remove_dir_all(&dir).expect("the old scratch directory is removed");
    }
    std::fs::create_dir_all(&dir).expect("the scratch directory is made");
    dir.canonicalize()
        .expect("the scratch directory has a path")
}

fn make_removable(dir: &Path) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() && !path.is_symlink() {
            let _ = std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755));
            make_removable(&path);
        }
    }
}

/// `tests/golden/116` of the repository.
pub fn golden() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .canonicalize()
        .expect("the golden corpora are two directories above the tool")
}

/// The text of a file, inflated where it is gzip.
pub fn text_of(path: &Path) -> String {
    let bytes = std::fs::read(path).expect("the file is read");
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut text)
            .expect("the file inflates to UTF-8");
        text
    } else {
        String::from_utf8(bytes).expect("the file is UTF-8")
    }
}

/// The lines of an output that open with `#`, and the others.
pub fn header_and_data_lines(text: &str) -> (Vec<&str>, Vec<&str>) {
    text.lines().partition(|line| line.starts_with('#'))
}

/// The header lines and `rows`, one a line.
pub fn rebuilt(header: &[&str], rows: &[&str]) -> String {
    let mut text = String::new();
    for line in header.iter().chain(rows) {
        text.push_str(line);
        text.push('\n');
    }
    text
}

/// The report at `path` as the value it holds.
pub fn report(path: &Path) -> Value {
    let text = std::fs::read_to_string(path).expect("the report is read");
    pyjson::loads(&text).expect("the report is JSON")
}

/// The member of `value` at `path`, a name a level.
pub fn member<'a>(value: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(value, |value, name| {
        value
            .get(name)
            .unwrap_or_else(|| panic!("the report has no member {name} of {path:?}"))
    })
}

/// `value` as compact JSON, for an assertion on a part of a report.
pub fn compact(value: &Value) -> String {
    pyjson::dumps(
        value,
        &pyjson::DumpOptions {
            indent: None,
            sort_keys: false,
            ensure_ascii: false,
        },
    )
}

pub fn args_of(args: &[&str]) -> Vec<OsString> {
    args.iter().map(OsString::from).collect()
}
