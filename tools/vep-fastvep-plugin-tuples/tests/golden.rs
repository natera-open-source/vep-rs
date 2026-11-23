// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over the cases of `tests/fixtures/cases.tsv` and holds its exit code,
//! its stderr and its tuple file to the bytes of `tests/fixtures/expected/`, on one thread
//! and on several.

mod common;

use common::{Scratch, BIN};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const ROOT: &str = env!("CARGO_MANIFEST_DIR");

/// The cases `cases.tsv` holds; a table cut short fails here, not silently.
const CASES: usize = 58;

fn fixtures() -> PathBuf {
    Path::new(ROOT).join("tests/fixtures")
}

/// One line of `cases.tsv`: the case, its exit code, what is held to the goldens (`bytes`:
/// stderr and the tuple file; `exit`: the exit code alone, for a usage error or the help)
/// and the arguments, separated by blanks, `{tuples}` standing for a tuple file of a
/// scratch directory.
struct Case {
    name: String,
    exit: i32,
    bytes: bool,
    arguments: Vec<String>,
}

fn cases() -> Vec<Case> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    let mut lines = table.lines();
    assert_eq!(lines.next(), Some("case\texit\tcheck\targuments"));
    lines
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let [name, exit, check, arguments] = fields.as_slice() else {
                panic!("four fields a case, not {line:?}");
            };
            Case {
                name: name.to_string(),
                exit: exit.parse().unwrap(),
                bytes: match *check {
                    "bytes" => true,
                    "exit" => false,
                    other => panic!("{name}: unknown check {other:?}"),
                },
                arguments: arguments
                    .split(' ')
                    .filter(|argument| !argument.is_empty())
                    .map(str::to_string)
                    .collect(),
            }
        })
        .collect()
}

/// Runs the binary from the package root, so a fixture path is the relative one the
/// goldens carry.
fn run(arguments: &[String], threads: &str) -> Output {
    Command::new(BIN)
        .current_dir(ROOT)
        .args(arguments)
        .env("VEP_FASTVEP_PLUGIN_TUPLES_THREADS", threads)
        .output()
        .unwrap()
}

fn golden(name: &str, kind: &str) -> Option<Vec<u8>> {
    fs::read(fixtures().join("expected").join(format!("{name}.{kind}"))).ok()
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// Every fault of the cases run on `threads` threads.
fn faults_on(threads: &str) -> Vec<String> {
    let scratch = Scratch::new();
    let cases = cases();
    assert_eq!(cases.len(), CASES);
    let mut faults = Vec::new();
    for case in &cases {
        let tuples = scratch.path(&format!("{}.tuples.tsv", case.name));
        let arguments: Vec<String> = case
            .arguments
            .iter()
            .map(|argument| argument.replace("{tuples}", &common::text(&tuples)))
            .collect();
        let output = run(&arguments, threads);
        if output.status.code() != Some(case.exit) {
            faults.push(format!(
                "{}: exit {:?}, expected {}: {}",
                case.name,
                output.status.code(),
                case.exit,
                text(&output.stderr)
            ));
        }
        if !case.bytes {
            continue;
        }
        if !output.stdout.is_empty() {
            faults.push(format!("{}: stdout {}", case.name, text(&output.stdout)));
        }
        let expected = golden(&case.name, "stderr.txt").expect("a golden stderr a case");
        if output.stderr != expected {
            faults.push(format!(
                "{}: stderr\n  wrote    {}\n  expected {}",
                case.name,
                text(&output.stderr),
                text(&expected)
            ));
        }
        // A refused run has no golden tuple file and writes none.
        let written = fs::read(&tuples).ok();
        if written != golden(&case.name, "tuples.tsv") {
            faults.push(format!(
                "{}: the tuple file\n{}",
                case.name,
                written.map_or("is absent".to_string(), |bytes| text(&bytes))
            ));
        }
    }
    faults
}

#[test]
fn every_case_writes_the_expected_bytes() {
    let faults = faults_on("1");
    assert!(faults.is_empty(), "{}", faults.join("\n"));
}

#[test]
fn every_case_writes_the_expected_bytes_on_several_threads() {
    let faults = faults_on("7");
    assert!(faults.is_empty(), "{}", faults.join("\n"));
}

#[test]
fn every_golden_belongs_to_a_case_that_holds_bytes() {
    let cases = cases();
    let mut names: Vec<String> = fs::read_dir(fixtures().join("expected"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    for name in &names {
        let case = name
            .strip_suffix(".stderr.txt")
            .or_else(|| name.strip_suffix(".tuples.tsv"))
            .unwrap_or_else(|| panic!("{name}: neither a stderr nor a tuple file"));
        assert!(
            cases.iter().any(|known| known.name == case && known.bytes),
            "{name}: no case holds it"
        );
    }
    let held = cases.iter().filter(|case| case.bytes).count();
    let stderr = names
        .iter()
        .filter(|name| name.ends_with(".stderr.txt"))
        .count();
    assert_eq!(stderr, held);
}
