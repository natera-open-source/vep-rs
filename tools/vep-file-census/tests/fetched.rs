// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `fetched` held to its goldens under `tests/fixtures/expected/`: `<name>.out` is the
//! standard output and `<name>.report` the report over `tests/fixtures/fetched/`. The
//! ETags of `fetched/list.tsv` are those of the files themselves:
//!
//! ```text
//! md5sum root/data/one.txt                                      # data/one.txt, the md5
//! md5sum root/data/two.txt | cut -c1-32 | xxd -r -p | md5sum    # data/two.txt, one part: <that>-1
//! sha256sum root/three.txt                                      # three.txt, the sha256
//! ```
//!
//! and the rows after them state a count of parts no part size gives, another file's
//! digest, another size, a missing file, a directory and an ETag of no shape.

mod common;

use common::{expected, run, Scratch};
use std::io::Write;

const CASES: &[(&str, &str, &str, i32)] = &[
    ("fetched_list", "fetched/list.tsv", "fetched/root", 3),
    ("fetched_held", "fetched/held.tsv", "fetched/root", 0),
    // A key with a space, a line ended by `\r\n` and a last line without a newline.
    (
        "fetched_no_final_newline",
        "fetched/no_final_newline.tsv",
        "fetched/root",
        0,
    ),
    ("fetched_no_rows", "names_none.txt", "fetched/root", 0),
    (
        "fetched_two_columns",
        "fetched/two_columns.tsv",
        "fetched/root",
        2,
    ),
    (
        "fetched_no_number",
        "fetched/no_number.tsv",
        "fetched/root",
        2,
    ),
    ("fetched_no_root", "fetched/held.tsv", "fetched/missing", 2),
];

#[test]
fn the_verdicts_the_report_and_the_status_are_fixed() {
    for (name, list, root, status) in CASES {
        for jobs in ["1", "2", "7"] {
            let scratch = Scratch::new(&format!("{name}-{jobs}"));
            let report = scratch.path("report.tsv");
            let run = run(
                &[
                    "fetched",
                    "--list",
                    list,
                    "--root",
                    root,
                    "--jobs",
                    jobs,
                    "--report",
                    report.to_str().unwrap(),
                ],
                None,
                None,
            );
            assert_eq!(
                run.status.code(),
                Some(*status),
                "{name} with {jobs} job(s)"
            );
            assert_eq!(
                run.stdout,
                expected(&format!("{name}.out")),
                "{name}: standard output"
            );
            if *status == 2 {
                assert!(!report.exists(), "{name}: a refused run writes no report");
                assert_eq!(
                    String::from_utf8_lossy(&run.stderr)
                        .lines()
                        .filter(|line| line.starts_with("ERROR: [vep-file-census] "))
                        .count(),
                    1
                );
            } else {
                assert_eq!(
                    std::fs::read(&report).unwrap(),
                    expected(&format!("{name}.report")),
                    "{name}: the report"
                );
            }
        }
    }
}

#[test]
fn fewer_than_one_job_is_refused() {
    let run = run(
        &[
            "fetched",
            "--list",
            "fetched/held.tsv",
            "--root",
            "fetched/root",
            "--jobs",
            "0",
        ],
        None,
        None,
    );
    assert_eq!(run.status.code(), Some(2));
    assert!(run.stdout.is_empty());
    assert!(String::from_utf8_lossy(&run.stderr).contains("--jobs takes 1 or more, not 0"));
}

/// The md5 of `bytes` as its sixteen bytes.
fn md5(bytes: &[u8]) -> [u8; 16] {
    md5::compute(bytes).0
}

// { head -c 5242880 F | md5sum; tail -c +5242881 F | md5sum; } | cut -c1-32 | xxd -r -p | md5sum
#[test]
fn an_object_uploaded_in_two_parts_is_held_by_the_part_size_that_gives_its_etag() {
    let scratch = Scratch::new("two-parts");
    std::fs::create_dir(scratch.path("root")).unwrap();
    // Six mebibytes of bytes that repeat at no part's length.
    let bytes: Vec<u8> = (0..6u32 << 20)
        .map(|at| (at.wrapping_mul(2_654_435_761) >> 24) as u8)
        .collect();
    std::fs::write(scratch.path("root/object.bin"), &bytes).unwrap();
    let mut list = std::fs::File::create(scratch.path("list.tsv")).unwrap();
    for (part_mib, key) in [(5usize, "object.bin"), (8, "object.bin")] {
        let digests: Vec<u8> = bytes.chunks(part_mib << 20).flat_map(md5).collect();
        let parts = bytes.len().div_ceil(part_mib << 20);
        let etag = md5(&digests)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>();
        writeln!(list, "{key}\t{}\t{etag}-{parts}", bytes.len()).unwrap();
    }
    // Two parts that no part size gives: the md5 of the whole file with a count of 2.
    let whole: String = md5(&bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    writeln!(list, "object.bin\t{}\t{whole}-2", bytes.len()).unwrap();
    drop(list);
    let report = scratch.path("report.tsv");
    let run = run(
        &[
            "fetched",
            "--list",
            scratch.path("list.tsv").to_str().unwrap(),
            "--root",
            scratch.path("root").to_str().unwrap(),
            "--report",
            report.to_str().unwrap(),
        ],
        None,
        Some("2"),
    );
    assert_eq!(run.status.code(), Some(0));
    assert_eq!(
        String::from_utf8(run.stdout).unwrap(),
        "object.bin\tbytes-alone\nheld: 0 by md5, 2 by multipart md5, 1 by bytes alone, 0 differ of 3\n"
    );
    assert_eq!(
        std::fs::read_to_string(&report).unwrap(),
        "object.bin\t6291456\tmultipart-md5 5\nobject.bin\t6291456\tmultipart-md5 8\nobject.bin\t6291456\tbytes-alone\n"
    );
}
