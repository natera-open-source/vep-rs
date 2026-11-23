// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the test files share: a scratch directory removed on drop, the fixtures, a run of
//! the binary, and the placeholders a golden holds for what depends on the checkout and the
//! run.

#![allow(dead_code)]

use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub const BIN: &str = env!("CARGO_BIN_EXE_vep-prepare-benchmark-vcfs");
pub const ROOT: &str = env!("CARGO_MANIFEST_DIR");
pub const THREADS_ENV: &str = "VEP_PREPARE_BENCHMARK_VCFS_THREADS";

/// The usage block every refused command line starts with on stderr.
pub const USAGE: &str = "\
usage: vep-prepare-benchmark-vcfs [-h] --input-dir INPUT_DIR
                                  --output-dir OUTPUT_DIR
                                  [--chrom-mode {auto,keep,strip,add}]
                                  [--cache-dir CACHE_DIR]
                                  [--max-variants MAX_VARIANTS]
                                  [--assembly ASSEMBLY]
                                  [--canonical-contigs [FILE_OR_LIST]]
                                  [--no-canonical-contigs]
";

pub struct Scratch(pub PathBuf);

impl Scratch {
    /// A directory of its own for one test, named without a link on the way so that the
    /// manifest names it as given; `tag` is unique within a test file.
    pub fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().canonicalize().unwrap().join(format!(
            "vep-prepare-benchmark-vcfs-test-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    pub fn at(&self, rel: &str) -> PathBuf {
        self.0.join(rel)
    }

    /// Where the runs of a test ask for their output.
    pub fn out(&self) -> PathBuf {
        self.0.join("out")
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// The fixtures directory, named without a link on the way.
pub fn fixtures() -> PathBuf {
    Path::new(ROOT)
        .join("tests/fixtures")
        .canonicalize()
        .unwrap()
}

pub fn fixture(name: &str) -> PathBuf {
    fixtures().join(name)
}

/// Runs the binary with exactly `args` and the thread bound `threads`, the host's cores
/// where it is None.
pub fn run_with<S: AsRef<OsStr>>(args: &[S], threads: Option<&str>) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(ROOT).args(args).env_remove(THREADS_ENV);
    if let Some(threads) = threads {
        command.env(THREADS_ENV, threads);
    }
    command.output().unwrap()
}

/// Runs the binary with exactly `args`.
pub fn run<S: AsRef<OsStr>>(args: &[S]) -> Output {
    run_with(args, None)
}

/// Runs the binary over `input` into `output` with `options` after the two directories.
pub fn prepare(input: &Path, output: &Path, options: &[&str], threads: Option<&str>) -> Output {
    let mut args: Vec<&OsStr> = vec![
        OsStr::new("--input-dir"),
        input.as_os_str(),
        OsStr::new("--output-dir"),
        output.as_os_str(),
    ];
    args.extend(options.iter().map(OsStr::new));
    run_with(&args, threads)
}

pub fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

pub fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

pub fn assert_clean_exit(output: &Output, what: &str) {
    assert!(
        output.status.success(),
        "{what}: exit {:?}: {}",
        output.status.code(),
        stderr(output)
    );
    assert!(output.stderr.is_empty(), "{what}: {}", stderr(output));
}

/// The names in a directory, sorted.
pub fn names_in(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// A path as `json.dumps` writes it.
pub fn json_escaped(path: &Path) -> String {
    let mut out = String::new();
    for c in path.to_str().expect("a UTF-8 test path").chars() {
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

fn replace_all(text: &[u8], from: &[u8], to: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.windows(from.len()).position(|window| window == from) {
        out.extend_from_slice(&rest[..at]);
        out.extend_from_slice(to);
        rest = &rest[at + from.len()..];
    }
    out.extend_from_slice(rest);
    out
}

/// A golden with its placeholders filled in: each `{{NAME}}` of `values` by its value.
pub fn filled(golden: &[u8], values: &[(&str, &str)]) -> Vec<u8> {
    values.iter().fold(golden.to_vec(), |text, (name, value)| {
        replace_all(
            &text,
            format!("{{{{{name}}}}}").as_bytes(),
            value.as_bytes(),
        )
    })
}

/// The manifest a run wrote, its `generated_at_utc` value replaced by `{{GENERATED_AT}}`
/// once it is seen to be a UTC time as `isoformat()` writes one: the date, `T`, the time,
/// six decimals or none, and `+00:00`.
pub fn manifest_without_its_time(output_dir: &Path) -> Vec<u8> {
    let text = fs::read(output_dir.join("manifest.json")).unwrap();
    let key = b"  \"generated_at_utc\": \"";
    let start = text
        .windows(key.len())
        .position(|window| window == key)
        .expect("the manifest holds its time")
        + key.len();
    let length = text[start..]
        .iter()
        .position(|byte| *byte == b'"')
        .expect("the time is a string");
    let stamp = std::str::from_utf8(&text[start..start + length]).unwrap();
    let shape: String = stamp
        .chars()
        .map(|c| if c.is_ascii_digit() { '9' } else { c })
        .collect();
    assert!(
        shape == "9999-99-99T99:99:99+99:99" || shape == "9999-99-99T99:99:99.999999+99:99",
        "generated_at_utc is {stamp}"
    );
    assert!(stamp.ends_with("+00:00"), "generated_at_utc is {stamp}");
    [&text[..start], b"{{GENERATED_AT}}", &text[start + length..]].concat()
}
