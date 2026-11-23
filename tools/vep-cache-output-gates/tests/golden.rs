// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary held to the goldens: every case of `tests/fixtures/cases.tsv` is run in an
//! empty directory over the reports it names, and its standard output, its verdict report
//! and its exit code are compared with `tests/fixtures/expected/` byte for byte.

mod common;

use common::{error_line, fixture, golden, report, run_in, Scratch};
use std::collections::BTreeSet;
use std::ffi::OsString;
use std::fs;

struct Case {
    name: String,
    /// Baseline metrics, builder metrics, cache comparison report, concordance summary.
    reports: [String; 4],
    thresholds: [String; 5],
    status: i32,
}

fn cases() -> Vec<Case> {
    let table = fs::read_to_string(fixture("cases.tsv")).unwrap();
    let mut lines = table.lines();
    assert_eq!(
        lines.next(),
        Some(
            "case\tbaseline\tbuilder\tcache_diff\tconcordance\tf1\ttotal\tcolocated\tfrequency\tconsequence\tstatus"
        )
    );
    lines
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            assert_eq!(cells.len(), 11, "{line}");
            Case {
                name: cells[0].to_string(),
                reports: std::array::from_fn(|i| cells[1 + i].to_string()),
                // A threshold's line feed stands in its cell as a backslash and an `n`.
                thresholds: std::array::from_fn(|i| cells[5 + i].replace("\\n", "\n")),
                status: cells[10].parse().unwrap(),
            }
        })
        .collect()
}

/// The operand of one report cell: `-` is the empty operand, `absent` a name in the run's
/// directory that names nothing, anything else a file of `tests/fixtures/reports/`.
fn operand(cell: &str, scratch: &Scratch) -> OsString {
    match cell {
        "-" => OsString::new(),
        "absent" => scratch.0.join("absent.json").into_os_string(),
        name => {
            let path = report(name);
            assert!(path.is_file(), "{name}");
            path.into_os_string()
        }
    }
}

#[test]
fn every_case_matches_its_goldens() {
    let cases = cases();
    let mut failures = Vec::new();
    for case in &cases {
        let scratch = Scratch::new(&format!("golden-{}", case.name));
        let mut args: Vec<OsString> = case
            .reports
            .iter()
            .map(|cell| operand(cell, &scratch))
            .collect();
        args.push("verdict.json".into());
        args.extend(case.thresholds.iter().map(OsString::from));
        let output = run_in(&scratch.0, &args);

        let mut wrong = Vec::new();
        if output.status.code() != Some(case.status) {
            wrong.push(format!(
                "exit code {:?}, expected {}",
                output.status.code(),
                case.status
            ));
        }
        if output.stdout != golden(&format!("{}.stdout", case.name)) {
            wrong.push(format!(
                "standard output:\n{}",
                String::from_utf8_lossy(&output.stdout)
            ));
        }
        let written = fs::read(scratch.0.join("verdict.json")).ok();
        let expected =
            fs::read(fixture("expected").join(format!("{}.verdict.json", case.name))).ok();
        if written != expected {
            wrong.push(format!(
                "verdict report:\n{}",
                written.map_or("none written".to_string(), |bytes| {
                    String::from_utf8_lossy(&bytes).into_owned()
                })
            ));
        }
        // A verdict is on standard error only where none was reached, as one line.
        if case.status == 2 {
            error_line(&output);
            assert!(expected.is_none(), "{}", case.name);
        } else {
            assert!(expected.is_some(), "{}", case.name);
            if !output.stderr.is_empty() {
                wrong.push(format!(
                    "standard error:\n{}",
                    String::from_utf8_lossy(&output.stderr)
                ));
            }
        }
        if !wrong.is_empty() {
            failures.push(format!("{}: {}", case.name, wrong.join("\n")));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n\n"));
}

#[test]
fn the_cases_reach_each_exit_code_and_own_every_golden() {
    let cases = cases();
    for status in [0, 1, 2] {
        assert!(cases.iter().any(|case| case.status == status), "{status}");
    }
    let names: BTreeSet<&str> = cases.iter().map(|case| case.name.as_str()).collect();
    assert_eq!(names.len(), cases.len());
    let mut goldens = 0;
    for entry in fs::read_dir(fixture("expected")).unwrap() {
        let file = entry.unwrap().file_name().into_string().unwrap();
        let case = file
            .strip_suffix(".stdout")
            .or_else(|| file.strip_suffix(".verdict.json"))
            .unwrap_or_else(|| panic!("{file} is no golden of a case"));
        assert!(names.contains(case), "{file}");
        goldens += 1;
    }
    let verdicts = cases.iter().filter(|case| case.status != 2).count();
    assert_eq!(goldens, cases.len() + verdicts);
}
