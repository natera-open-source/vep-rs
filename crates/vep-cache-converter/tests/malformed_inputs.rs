// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A flat cache whose region and variation files are named for the parent
//! directory (`.._1-1000000.json`, `.._all_vars.gz` and its index): the
//! converter skips them with its unexpected-name warning and writes nothing
//! outside the output directory, whose entries are `info.json`, the
//! `transcripts/` directory and the checksum manifest written over them.
//! Standard error that refuses every write changes neither a conversion nor
//! its exit status.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/parent_named_shards")
}

/// A scratch directory under the system temporary directory, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "vep-cache-converter-malformed-{}-{nanos}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The names of every entry directly under `dir`, sorted.
fn entries(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

#[test]
fn files_named_for_the_parent_directory_are_skipped_and_write_nothing_beside_the_output() {
    let scratch = Scratch::new();
    let out = scratch.0.join("out");
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .arg("--json-dir")
        .arg(fixture())
        .arg("--output-dir")
        .arg(&out)
        .output()
        .expect("vep-cache-converter runs");
    let err = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{err}");
    assert_eq!(
        err.matches("Skipping file with unexpected name").count(),
        3,
        "{err}"
    );
    assert_eq!(entries(&scratch.0), ["out"]);
    assert_eq!(
        entries(&out),
        ["CHECKSUMS.sha256", "info.json", "transcripts"]
    );
    assert_eq!(entries(&out.join("transcripts")), ["21"]);
    assert_eq!(entries(&out.join("transcripts/21")), ["1-1000000.json"]);
}

/// The converter's command on the flat cache `json_dir` into `out`.
fn command(json_dir: &Path, out: &Path) -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"));
    cmd.arg("--json-dir")
        .arg(json_dir)
        .arg("--output-dir")
        .arg(out);
    cmd
}

/// The device that refuses every write, where the platform has one.
const FULL_DEVICE: &str = "/dev/full";

/// The streams that take no bytes, by name: the device that refuses every
/// write, where the platform has one, and a pipe whose reader has gone.
fn unwritable_streams() -> &'static [&'static str] {
    if Path::new(FULL_DEVICE).exists() {
        &[FULL_DEVICE, "a pipe whose reader has gone"]
    } else {
        &["a pipe whose reader has gone"]
    }
}

/// The stream [`unwritable_streams`] names as `stream`, opened for a child
/// process to write.
fn unwritable(stream: &str) -> Stdio {
    if stream == FULL_DEVICE {
        return std::fs::File::options()
            .write(true)
            .open(stream)
            .unwrap()
            .into();
    }
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    writer.into()
}

/// Every file below `dir`, by its path from `dir`, with its bytes.
fn files(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut found = BTreeMap::new();
    let mut pending = vec![dir.to_path_buf()];
    while let Some(next) = pending.pop() {
        for entry in std::fs::read_dir(&next).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                pending.push(path);
            } else {
                let name = path.strip_prefix(dir).unwrap().display().to_string();
                found.insert(name, std::fs::read(&path).unwrap());
            }
        }
    }
    found
}

/// A conversion whose every log line standard error refuses, on a device that
/// refuses every write or a pipe whose reader has gone, completes: it exits 0
/// and writes the files of the same conversion on a readable stream, byte for
/// byte.
#[test]
fn a_conversion_completes_on_an_unwritable_error_stream() {
    let scratch = Scratch::new();
    let readable = scratch.0.join("readable");
    let output = command(&fixture(), &readable)
        .output()
        .expect("vep-cache-converter runs");
    let err = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{err}");
    assert!(err.contains("Conversion complete: "), "{err}");
    let expected = files(&readable);
    assert!(
        expected.contains_key("transcripts/21/1-1000000.json"),
        "{:?}",
        expected.keys()
    );
    for (n, stream) in unwritable_streams().iter().enumerate() {
        let out = scratch.0.join(format!("unwritable{n}"));
        let status = command(&fixture(), &out)
            .stdout(Stdio::null())
            .stderr(unwritable(stream))
            .status()
            .expect("vep-cache-converter runs");
        assert!(status.success(), "{stream}: {status}");
        assert!(files(&out) == expected, "{stream}: the files differ");
    }
}

/// A refused conversion keeps its exit status of 1 on such a stream, though
/// its log lines and its error are lost: a flat cache whose `info.json` is an
/// empty object is refused after the transcript files are logged, and the
/// output directory is left without the `info.json` that marks a cache.
#[test]
fn a_refused_conversion_exits_1_on_an_unwritable_error_stream() {
    let scratch = Scratch::new();
    let flat = scratch.0.join("flat");
    std::fs::create_dir_all(flat.join("transcripts")).unwrap();
    std::fs::write(flat.join("info.json"), "{}").unwrap();
    std::fs::write(flat.join("transcripts/21_1-1000000.json"), "[]\n").unwrap();

    let output = command(&flat, &scratch.0.join("readable"))
        .output()
        .expect("vep-cache-converter runs");
    let err = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let logged = err.find("Total transcripts: ");
    let refused = err.find("info.json lacks a non-empty species string");
    assert!(
        matches!((logged, refused), (Some(logged), Some(refused)) if logged < refused),
        "a log line precedes the error:\n{err}"
    );
    for (n, stream) in unwritable_streams().iter().enumerate() {
        let out = scratch.0.join(format!("unwritable{n}"));
        let status = command(&flat, &out)
            .stdout(Stdio::null())
            .stderr(unwritable(stream))
            .status()
            .expect("vep-cache-converter runs");
        assert_eq!(status.code(), Some(1), "{stream}: {status}");
        assert!(!out.join("info.json").exists(), "{stream}");
    }
}
