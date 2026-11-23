// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over the cases of `tests/fixtures/cases.tsv` and holds its exit code,
//! stdout, stderr and report file to the bytes of `tests/fixtures/expected/`, then over the
//! command lines that name files outside the fixtures.

use std::ffi::{OsStr, OsString};
use std::fs;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-stats-compare");
const ROOT: &str = env!("CARGO_MANIFEST_DIR");

/// The cases `cases.tsv` holds; a table cut short fails here, not silently.
const CASES: usize = 77;

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-stats-compare-golden-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn fixtures() -> PathBuf {
    Path::new(ROOT).join("tests/fixtures")
}

/// One line of `cases.tsv`: the case, its exit code, what is held to the goldens (`bytes`:
/// stdout, stderr and the report; `exit`: the exit code alone, for a usage error or the
/// help) and the arguments, `{report}` standing for a report path under a scratch directory.
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
            let mut fields = line.split('\t');
            let name = fields.next().unwrap().to_string();
            let exit = fields.next().unwrap().parse().unwrap();
            let bytes = match fields.next().unwrap() {
                "bytes" => true,
                "exit" => false,
                other => panic!("{name}: unknown check {other:?}"),
            };
            // An argument writes a tab as `\t` and a backslash as `\\`.
            let arguments = fields
                .map(|field| {
                    field
                        .replace("\\\\", "\0")
                        .replace("\\t", "\t")
                        .replace('\0', "\\")
                })
                .collect();
            Case {
                name,
                exit,
                bytes,
                arguments,
            }
        })
        .collect()
}

/// Runs the binary from the package root, so a fixture path is the relative one the
/// goldens carry.
fn run(arguments: impl IntoIterator<Item = impl AsRef<OsStr>>) -> Output {
    Command::new(BIN)
        .current_dir(ROOT)
        .args(arguments)
        .output()
        .unwrap()
}

fn golden(name: &str, kind: &str) -> Option<Vec<u8>> {
    fs::read(fixtures().join("expected").join(format!("{name}.{kind}"))).ok()
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

#[test]
fn every_case_writes_the_expected_bytes() {
    let scratch = Scratch::new("cases");
    let cases = cases();
    assert_eq!(cases.len(), CASES);
    let mut faults = Vec::new();
    for case in &cases {
        let report = scratch.0.join(&case.name).join("out").join("stats.json");
        let arguments = case
            .arguments
            .iter()
            .map(|argument| match argument.as_str() {
                "{report}" => report.clone().into_os_string(),
                other => OsString::from(other),
            });
        let output = run(arguments);
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
        let written = fs::read(&report).ok();
        for (kind, actual, expected) in [
            (
                "stdout",
                Some(output.stdout),
                golden(&case.name, "stdout").or(Some(Vec::new())),
            ),
            (
                "stderr",
                Some(output.stderr),
                golden(&case.name, "stderr").or(Some(Vec::new())),
            ),
            ("json", written, golden(&case.name, "json")),
        ] {
            if actual != expected {
                faults.push(format!(
                    "{}: {kind} differs\n--- expected\n{}\n--- actual\n{}",
                    case.name,
                    expected.map_or("(none)".to_string(), |bytes| text(&bytes)),
                    actual.map_or("(none)".to_string(), |bytes| text(&bytes)),
                ));
            }
        }
    }
    assert!(faults.is_empty(), "{}", faults.join("\n"));
}

#[test]
fn every_golden_belongs_to_a_case() {
    let cases = cases();
    let mut names: Vec<String> = fs::read_dir(fixtures().join("expected"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    let mut expected = Vec::new();
    for case in cases.iter().filter(|case| case.bytes) {
        for kind in ["json", "stderr", "stdout"] {
            if golden(&case.name, kind).is_some() {
                expected.push(format!("{}.{kind}", case.name));
            }
        }
    }
    expected.sort();
    assert_eq!(names, expected);
    // A case held to its bytes writes to stdout or to stderr, never to neither.
    for case in cases.iter().filter(|case| case.bytes) {
        let quiet =
            golden(&case.name, "stdout").is_none() && golden(&case.name, "stderr").is_none();
        assert!(!quiet, "{} has no golden", case.name);
    }
}

/// The command line on files named by absolute paths: the paths under `expected` and
/// `actual`, the JSON report written, `stats: <verdict>` on stdout, exit 0 on green, 1 on
/// red, 2 on an error.
#[test]
fn command_line_writes_the_report_the_verdict_and_the_exit_code() {
    let scratch = Scratch::new("command-line");
    let reference = fixtures().join("reference.summary.log");
    let reference = reference.to_str().unwrap();
    let actual = scratch.0.join("candidate_summary.txt");
    let actual_arg = actual.to_str().unwrap();
    fs::copy(fixtures().join("candidate.summary.log"), &actual).unwrap();
    let written = scratch.0.join("out").join("stats.json");
    let output = run([
        "--expected",
        reference,
        "--actual",
        actual_arg,
        "--report",
        written.to_str().unwrap(),
    ]);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(text(&output.stdout), "stats: green\n");
    assert!(output.stderr.is_empty());
    let on_disk = fs::read_to_string(&written).unwrap();
    assert!(on_disk.starts_with("{\n \"mode\": \"stats\",\n \"format\": \"stats\",\n \"verdict\": \"green\",\n \"reasons\": [],\n"));
    assert!(on_disk.ends_with(&format!(
        " \"differences\": {{}},\n \"expected\": \"{reference}\",\n \"actual\": \"{actual_arg}\"\n}}\n"
    )));
    fs::copy(fixtures().join("one_value.summary.log"), &actual).unwrap();
    let output = run(["--expected", reference, "--actual", actual_arg]);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        text(&output.stdout),
        "stats: red\n  [General statistics]: 1 difference\n"
    );
    let output = run([
        "--expected",
        reference,
        "--actual",
        actual_arg,
        "--mask-row",
        "Variants processed",
    ]);
    assert_eq!(output.status.code(), Some(0));
    let missing = scratch.0.join("none.txt");
    let missing = missing.to_str().unwrap();
    let output = run(["--expected", reference, "--actual", missing]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        text(&output.stderr),
        format!("ERROR: [compare_vep_stats] cannot read {missing}: [Errno 2] No such file or directory: '{missing}'\n")
    );
    assert!(output.stdout.is_empty());
    let output = run([
        "--expected",
        reference,
        "--actual",
        actual_arg,
        "--tolerance",
        "-1",
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        text(&output.stderr),
        "ERROR: [compare_vep_stats] --tolerance -1 is negative\n"
    );
    // A gzipped report reads like a plain one; a corrupt one is an error line.
    let zipped = fixtures().join("reference.summary.log.gz");
    let output = run([
        "--expected",
        zipped.to_str().unwrap(),
        "--actual",
        reference,
    ]);
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(text(&output.stdout), "stats: green\n");
    let corrupt = fixtures().join("corrupt.summary.log.gz");
    let corrupt = corrupt.to_str().unwrap();
    let output = run(["--expected", corrupt, "--actual", reference]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        text(&output.stderr),
        format!("ERROR: [compare_vep_stats] cannot read {corrupt}: Error -3 while decompressing data: invalid bit length repeat\n")
    );
}

#[test]
fn a_report_is_written_only_where_one_is_asked_for_and_the_run_succeeds() {
    let scratch = Scratch::new("report");
    let reference = fixtures().join("reference.summary.log");
    let report = scratch.0.join("deep").join("er").join("stats.json");
    let output = run([
        OsStr::new("--expected"),
        reference.as_os_str(),
        OsStr::new("--actual"),
        fixtures().join("not_a_report.txt").as_os_str(),
        OsStr::new("--report"),
        report.as_os_str(),
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert!(!scratch.0.join("deep").exists());
    // A red verdict writes the report as a green one does, over a report already there.
    fs::create_dir_all(report.parent().unwrap()).unwrap();
    fs::write(
        &report,
        "an older report, longer than nothing\n".repeat(400),
    )
    .unwrap();
    let output = run([
        OsStr::new("--expected"),
        reference.as_os_str(),
        OsStr::new("--actual"),
        fixtures().join("one_value.summary.log").as_os_str(),
        OsStr::new("--report"),
        report.as_os_str(),
    ]);
    assert_eq!(output.status.code(), Some(1));
    let on_disk = fs::read_to_string(&report).unwrap();
    assert!(on_disk.starts_with("{\n \"mode\": \"stats\","));
    assert!(on_disk.contains("\n \"verdict\": \"red\",\n"));
    assert!(on_disk.ends_with("\"\n}\n"));
}

/// An argument that is not UTF-8 is written as Python writes the `str` it decodes it to:
/// `\udcxx` for each undecodable byte, in the report and on stderr alike.
#[test]
fn arguments_that_are_not_utf8_are_written_with_their_escapes() {
    let scratch = Scratch::new("bytes");
    let reference = fixtures().join("reference.summary.log");
    let report = scratch.0.join("stats.json");
    let output = run([
        OsStr::new("--expected"),
        reference.as_os_str(),
        OsStr::new("--actual"),
        fixtures().join("candidate.summary.log").as_os_str(),
        OsStr::new("--mask-row"),
        OsStr::from_bytes(b"\xffrow \xc3\xa9"),
        OsStr::new("--report"),
        report.as_os_str(),
    ]);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let on_disk = fs::read_to_string(&report).unwrap();
    assert!(on_disk.contains(
        "\n  \"Run time\",\n  \"Start time\",\n  \"\\udcffrow \\u00e9\"\n ],\n \"sections\": {"
    ));
    let mut missing = scratch.0.clone().into_os_string().into_vec();
    missing.extend_from_slice(b"/none\xff\xc3\xa9.txt");
    let missing = OsString::from_vec(missing);
    let output = run([
        OsStr::new("--expected"),
        reference.as_os_str(),
        OsStr::new("--actual"),
        missing.as_os_str(),
    ]);
    assert_eq!(output.status.code(), Some(2));
    let shown = format!("{}/none\\udcff\u{e9}.txt", scratch.0.display());
    // The kernel's refusal of the name differs by system; the name's spelling does not.
    let stderr = text(&output.stderr);
    assert!(
        stderr.starts_with(&format!(
            "ERROR: [compare_vep_stats] cannot read {shown}: [Errno "
        )),
        "{stderr}"
    );
    assert!(stderr.ends_with(&format!(": '{shown}'\n")), "{stderr}");
}
