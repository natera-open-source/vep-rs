// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary against stored runs of its command line.
//!
//! `tests/data/cli/cases.tsv` holds one command line a row with the exit code, the
//! standard output and the standard error the run must end with, and
//! `tests/data/cli/reports` the report of each that wrote one. A field is ASCII: `\\`,
//! `\t`, `\n`, `\r` and `\xNN` for every other byte outside 0x20 to 0x7e; the arguments
//! are joined by one space, an empty argument written `\e`. `<<G>>` stands for the
//! `GRCh38-release116` corpus, `<<T>>` for the directory the planted files are in
//! (`<<T->>` without its leading `/`), `<<usage>>` for the usage block, `<<help>>` for the
//! whole help and `<<prog>>` for the program's name.
//!
//! Three things are the binary's own and are written so in the rows: the program's name;
//! the usage block, which is the `argparse` layout of CPython at 80 columns under that
//! name; and the one line that ends a run on an input no reading of the comparison
//! expects. No row holds a run at another width, in colour or on a terminal.

mod support;

use std::ffi::OsString;
use std::os::unix::ffi::OsStringExt;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use vep_formats_compare::pyjson;

const DATA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/cli");
const USAGE: &str = include_str!("../data/usage.txt");
const HELP: &str = include_str!("../data/help.txt");

struct Case {
    name: String,
    group: String,
    code: i32,
    env: Vec<(OsString, OsString)>,
    argv: Vec<OsString>,
    stdout: Vec<u8>,
    stderr: Vec<u8>,
    report: Option<PathBuf>,
}

fn unescape(field: &str) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(field.len());
    let mut rest = field.as_bytes();
    while let Some((&byte, tail)) = rest.split_first() {
        if byte != b'\\' {
            bytes.push(byte);
            rest = tail;
            continue;
        }
        let (escape, tail) = tail.split_first().expect("an escape has a letter");
        rest = tail;
        match escape {
            b'\\' => bytes.push(b'\\'),
            b't' => bytes.push(b'\t'),
            b'n' => bytes.push(b'\n'),
            b'r' => bytes.push(b'\r'),
            b'e' => {}
            b'x' => {
                let digits = std::str::from_utf8(&rest[..2]).expect("two hex digits");
                bytes.push(u8::from_str_radix(digits, 16).expect("two hex digits"));
                rest = &rest[2..];
            }
            other => panic!("unknown escape \\{}", *other as char),
        }
    }
    bytes
}

fn replaced(bytes: &[u8], from: &str, to: &[u8]) -> Vec<u8> {
    let from = from.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut at = 0;
    while at < bytes.len() {
        if bytes[at..].starts_with(from) {
            out.extend_from_slice(to);
            at += from.len();
        } else {
            out.push(bytes[at]);
            at += 1;
        }
    }
    out
}

/// What stands for the placeholders of a row: the corpus, the planted directory and the
/// tool's own texts.
struct Places {
    corpus: PathBuf,
    planted: PathBuf,
}

impl Places {
    fn fill(&self, bytes: &[u8]) -> Vec<u8> {
        let corpus = self.corpus.to_str().expect("the corpus path is UTF-8");
        let planted = self.planted.to_str().expect("the scratch path is UTF-8");
        let mut bytes = bytes.to_vec();
        for (from, to) in [
            ("<<G>>", corpus),
            ("<<T>>", planted),
            ("<<T->>", planted.trim_start_matches('/')),
            ("<<usage>>", USAGE),
            ("<<help>>", HELP),
            ("<<prog>>", "vep-formats-compare"),
        ] {
            bytes = replaced(&bytes, from, to.as_bytes());
        }
        bytes
    }

    /// A stored report, the two directories spelt as a JSON text spells them.
    fn fill_report(&self, bytes: &[u8]) -> Vec<u8> {
        let spelt = |path: &Path| {
            let quoted = pyjson::dumps_str(path.to_str().expect("the path is UTF-8"), true);
            quoted[1..quoted.len() - 1].to_string()
        };
        let bytes = replaced(bytes, "<<G>>", spelt(&self.corpus).as_bytes());
        replaced(&bytes, "<<T>>", spelt(&self.planted).as_bytes())
    }

    fn word(&self, word: &str) -> OsString {
        OsString::from_vec(self.fill(&unescape(word)))
    }
}

fn cases(group: &str, places: &Places) -> Vec<Case> {
    let table = std::fs::read_to_string(format!("{DATA}/cases.tsv")).expect("the cases are read");
    let mut cases = Vec::new();
    for row in table.lines().skip(1) {
        let fields: Vec<&str> = row.split('\t').collect();
        let [name, row_group, code, env, argv, stdout, stderr, report] = fields[..] else {
            panic!("a row of cases.tsv has eight fields: {row}");
        };
        if row_group != group {
            continue;
        }
        let words = |field: &str| -> Vec<OsString> {
            field
                .split(' ')
                .filter(|word| !word.is_empty())
                .map(|word| places.word(word))
                .collect()
        };
        let env = env
            .split(' ')
            .filter(|word| !word.is_empty())
            .map(|word| {
                let (name, value) = word.split_once('=').expect("NAME=value");
                (places.word(name), places.word(value))
            })
            .collect();
        cases.push(Case {
            name: name.to_string(),
            group: row_group.to_string(),
            code: code.parse().expect("an exit code"),
            env,
            argv: words(argv),
            stdout: places.fill(&unescape(stdout)),
            stderr: places.fill(&unescape(stderr)),
            report: (report != "-").then(|| PathBuf::from(places.word(report))),
        });
    }
    assert!(!cases.is_empty(), "no case of the group {group}");
    cases
}

/// The files the command lines read, in a directory of this test's own: three warnings
/// files, one of them under a name that is not ASCII and under one that is not UTF-8, a
/// JSON list, a directory and a file that cannot be written, and the corpus's default
/// output without its last row.
fn planted(name: &str, corpus: &Path) -> PathBuf {
    let dir = support::scratch(name);
    for fixture in ["w.txt", "w-expected.txt", "w-actual.txt", "list.json"] {
        std::fs::copy(format!("{DATA}/fixtures/{fixture}"), dir.join(fixture))
            .expect("a fixture is copied");
    }
    let warnings = std::fs::read(dir.join("w.txt")).expect("the warnings file is read");
    std::fs::write(dir.join("w\u{e4}rn.txt"), &warnings).expect("a file is planted");
    let not_utf8 = OsString::from_vec(b"w\xff.txt".to_vec());
    std::fs::write(dir.join(not_utf8), &warnings).expect("a file is planted");
    std::fs::create_dir(dir.join("readonly")).expect("a directory is planted");
    std::fs::write(dir.join("readonly-file.json"), "{}\n").expect("a file is planted");
    for (path, mode) in [("readonly", 0o555), ("readonly-file.json", 0o444)] {
        std::fs::set_permissions(dir.join(path), std::fs::Permissions::from_mode(mode))
            .expect("the permissions are set");
    }
    let default = support::text_of(&corpus.join("expected/default.txt.gz"));
    let (header, rows) = support::header_and_data_lines(&default);
    let shorter = support::rebuilt(&header, &rows[..rows.len() - 1]);
    std::fs::write(dir.join("default-one-row-less.txt"), shorter).expect("a file is planted");
    dir
}

/// Whether the permissions of the planted directory keep this user from writing there;
/// they do not keep the superuser.
fn permissions_bind(planted: &Path) -> bool {
    let probe = planted.join("readonly/probe");
    match std::fs::write(&probe, "") {
        Ok(()) => {
            let _ = std::fs::remove_file(&probe);
            false
        }
        Err(_) => true,
    }
}

fn shown(bytes: &[u8]) -> String {
    let text = String::from_utf8_lossy(bytes);
    match text.char_indices().nth(600) {
        Some((at, _)) => format!("{}...", &text[..at]),
        None => text.into_owned(),
    }
}

/// Runs every case of `group` and holds the exit code, both streams and the report to
/// the stored ones.
fn hold(group: &str) {
    let corpus = support::golden().join("GRCh38-release116");
    let planted = planted(&format!("captures-{group}"), &corpus);
    let places = Places { corpus, planted };
    let binds = permissions_bind(&places.planted);
    let mut failures = Vec::new();
    let mut ran = 0;
    for case in cases(group, &places) {
        if !binds && case.name.starts_with("own-report-unwritable") {
            continue;
        }
        if let Some(report) = &case.report {
            let _ = std::fs::remove_file(report);
        }
        let run = support::run_in(&places.planted, &case.env, &case.argv);
        ran += 1;
        let mut differs = Vec::new();
        if run.code != case.code {
            differs.push(format!("exit code {}, expected {}", run.code, case.code));
        }
        for (stream, got, expected) in [
            ("standard output", &run.stdout, &case.stdout),
            ("standard error", &run.stderr, &case.stderr),
        ] {
            if got != expected {
                differs.push(format!(
                    "{stream}:\n{}\nexpected:\n{}",
                    shown(got),
                    shown(expected)
                ));
            }
        }
        if let Some(report) = &case.report {
            let template = std::fs::read(format!("{DATA}/reports/{}.json", case.name))
                .expect("the case's report is read");
            match std::fs::read(report) {
                Ok(written) if written == places.fill_report(&template) => {}
                Ok(written) => differs.push(format!("report:\n{}", shown(&written))),
                Err(error) => differs.push(format!("report: {error}")),
            }
        }
        if !differs.is_empty() {
            failures.push(format!(
                "{} ({}):\n{}",
                case.name,
                case.group,
                differs.join("\n")
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {ran} command lines differ:\n\n{}",
        failures.len(),
        failures.join("\n\n")
    );
}

#[test]
fn usage_errors_help_and_the_checks_before_any_read_end_as_their_stored_runs() {
    hold("parser");
}

#[test]
fn comparison_mode_runs_end_as_their_stored_runs() {
    hold("compare");
}

#[test]
fn warnings_mode_runs_end_as_their_stored_runs() {
    hold("warnings");
}

#[test]
fn self_parity_mode_runs_end_as_their_stored_runs() {
    hold("self-parity");
}
