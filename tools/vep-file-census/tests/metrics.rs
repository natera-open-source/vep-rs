// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `metrics` held to its goldens under `tests/fixtures/expected/`: `<name>.json` is the
//! report and `<name>.out` the summary of each fixture, in the bytes of the rule the
//! README states: the file read as UTF-8 text under universal newlines, the cells
//! stripped with `str.strip()`, the terms in `Counter.most_common()` order and the report
//! as `json.dump(indent=2)` writes it, then a newline.

mod common;

use common::{expected, run, Scratch, THREADS};

/// A fixture, the SECONDS word and the status; a status of 1 writes no report and no
/// summary.
const CASES: &[(&str, &str, &str, i32)] = &[
    ("metrics_output", "output.txt", "12", 0),
    ("metrics_no_final_newline", "output_no_newline.txt", "0", 0),
    // `int()` reads blanks, a sign and an underscore between digits.
    ("metrics_line_ends", "output_crlf.txt", " +0_07 ", 0),
    ("metrics_late_header", "output_late_header.txt", "3", 0),
    ("metrics_no_header", "output_no_header.txt", "1", 0),
    (
        "metrics_odd_terms",
        "output_odd_terms.txt",
        "99999999999999999999999",
        0,
    ),
    ("metrics_empty", "empty.vcf", "5", 0),
    ("metrics_not_utf8", "output_not_utf8.txt", "5", 1),
    ("metrics_missing", "missing.txt", "5", 1),
    ("metrics_seconds_no_number", "output.txt", "1.5", 1),
    ("metrics_spaced_name", "spaced name.vcf", "2", 0),
    // A gzip file is no UTF-8 text.
    ("metrics_gzip", "records.vcf.gz", "2", 1),
];

#[test]
fn the_report_and_the_summary_are_fixed_to_the_byte() {
    for (name, output, seconds, status) in CASES {
        for (at, threads) in THREADS.into_iter().enumerate() {
            let scratch = Scratch::new(&format!("{name}-{at}"));
            let report = scratch.path("baseline_metrics.json");
            let run = run(
                &["metrics", output, report.to_str().unwrap(), seconds],
                None,
                threads,
            );
            assert_eq!(
                run.status.code(),
                Some(*status),
                "{name} on {threads:?} thread(s)"
            );
            assert_eq!(
                run.stdout,
                expected(&format!("{name}.out")),
                "{name} on {threads:?} thread(s): the summary"
            );
            if *status == 0 {
                assert_eq!(
                    std::fs::read(&report).unwrap(),
                    expected(&format!("{name}.json")),
                    "{name} on {threads:?} thread(s): the report"
                );
            } else {
                assert!(!report.exists(), "{name}: a failed run writes no report");
            }
        }
    }
}

#[test]
fn a_report_that_cannot_be_written_leaves_no_summary() {
    let run = run(
        &[
            "metrics",
            "output.txt",
            "no_such_directory/report.json",
            "1",
        ],
        None,
        Some("2"),
    );
    assert_eq!(run.status.code(), Some(1));
    assert!(run.stdout.is_empty());
}
