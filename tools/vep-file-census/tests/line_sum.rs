// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `line-sum` held to its goldens under `tests/fixtures/expected/`: one line a file, the
//! count of its lines, a tab and their sum, in the bytes of the rule the README states.
//! Each golden is that rule in CPython over the fixture read as UTF-8 text with
//! `newline=""`: every line the file object yields without its `\r\n` end, an empty line
//! dropped, `hashlib.sha256` of the line's bytes, the digests added as integers and the
//! total printed as `format(total % 2**256, "064x")`.

mod common;

use common::{expected, fixtures, hold, line_sum, run, Case, THREADS};

const CASES: &[Case] = &[
    Case {
        name: "line_sum_empty",
        rule: "the count and sum of no line",
        stdin: None,
        args: &["line-sum", "empty.vcf"],
        status: 0,
    },
    Case {
        name: "line_sum_no_final_newline",
        rule: "the count and sum of the lines of records_no_newline.vcf",
        stdin: None,
        args: &["line-sum", "records_no_newline.vcf"],
        status: 0,
    },
    // Lines that end in `\r\n`, in `\n` and in a lone `\r`.
    Case {
        name: "line_sum_line_ends",
        rule: "the count and sum of the lines of output_crlf.txt",
        stdin: None,
        args: &["line-sum", "output_crlf.txt"],
        status: 0,
    },
    Case {
        name: "line_sum_rows_line_ends",
        rule: "the count and sum of the lines of output_crlf.txt that do not start with #",
        stdin: None,
        args: &[
            "line-sum",
            "--without-lines-starting",
            "#",
            "output_crlf.txt",
        ],
        status: 0,
    },
    Case {
        name: "line_sum_rows",
        rule: "the count and sum of the lines of output.txt that do not start with #",
        stdin: None,
        args: &["line-sum", "--without-lines-starting=#", "output.txt"],
        status: 0,
    },
    Case {
        name: "line_sum_records",
        rule: "the count and sum of the lines of records.vcf that do not start with #",
        stdin: None,
        args: &["line-sum", "--without-lines-starting", "#", "records.vcf"],
        status: 0,
    },
    Case {
        name: "line_sum_gzip",
        rule: "the count and sum of the lines records.vcf.gz inflates to",
        stdin: None,
        args: &["line-sum", "--gunzip", "records.vcf.gz"],
        status: 0,
    },
    Case {
        name: "line_sum_gzip_rows",
        rule: "the count and sum of the lines records.vcf.gz inflates to that do not start with #",
        stdin: None,
        args: &[
            "line-sum",
            "--without-lines-starting",
            "#",
            "--gunzip",
            "records.vcf.gz",
        ],
        status: 0,
    },
    Case {
        name: "line_sum_two_members",
        rule: "the count and sum of the lines two_members.vcf.gz inflates to",
        stdin: None,
        args: &["line-sum", "--gunzip", "two_members.vcf.gz"],
        status: 0,
    },
    Case {
        name: "line_sum_block_gzip",
        rule: "the count and sum of the lines chr_blocks.vcf.gz inflates to",
        stdin: None,
        args: &["line-sum", "--gunzip", "chr_blocks.vcf.gz"],
        status: 0,
    },
    Case {
        name: "line_sum_zero_padded",
        rule: "the count and sum of the lines zero_padded.vcf.gz inflates to",
        stdin: None,
        args: &["line-sum", "--gunzip", "zero_padded.vcf.gz"],
        status: 0,
    },
    // One line three times and one twice.
    Case {
        name: "line_sum_repeated",
        rule: "the count and sum of the lines of repeated_lines.txt",
        stdin: None,
        args: &["line-sum", "repeated_lines.txt"],
        status: 0,
    },
    // The same lines in another order, each ended by `\r\n` but the last, which has no end.
    Case {
        name: "line_sum_other_order",
        rule: "the count and sum of the lines of repeated_lines_other_order.txt",
        stdin: None,
        args: &["line-sum", "repeated_lines_other_order.txt"],
        status: 0,
    },
    Case {
        name: "line_sum_spaced_name",
        rule: "the count and sum of the lines of 'spaced name.vcf'",
        stdin: None,
        args: &["line-sum", "spaced name.vcf"],
        status: 0,
    },
    Case {
        name: "line_sum_several",
        rule: "one line a file, in the order of the operands",
        stdin: None,
        args: &[
            "line-sum",
            "records.vcf",
            "output.txt",
            "empty.vcf",
            "repeated_lines.txt",
        ],
        status: 0,
    },
    // gzip warns of the bytes after the member and leaves with 2; the member is summed.
    Case {
        name: "line_sum_gzip",
        rule: "the count and sum of the lines the member of garbage.vcf.gz inflates to",
        stdin: None,
        args: &["line-sum", "--gunzip", "garbage.vcf.gz"],
        status: 2,
    },
];

#[test]
fn the_count_and_the_sum_are_fixed_to_the_byte() {
    for case in CASES {
        hold(case);
    }
}

#[test]
fn the_same_lines_in_another_order_and_under_other_line_ends_sum_the_same() {
    assert_eq!(
        expected("line_sum_repeated.out"),
        expected("line_sum_other_order.out")
    );
    // A gzip file sums as the text it inflates to, in one member or in two.
    assert_eq!(
        expected("line_sum_gzip.out"),
        expected("line_sum_no_final_newline.out")
    );
    assert_eq!(
        expected("line_sum_gzip.out"),
        expected("line_sum_two_members.out")
    );
}

#[test]
fn a_file_that_is_no_utf8_text_is_summed_by_its_bytes() {
    let bytes = std::fs::read(fixtures().join("zero_bytes.bin")).unwrap();
    assert!(std::str::from_utf8(&bytes).is_err());
    let lines = bytes
        .split(|byte| matches!(byte, b'\n' | b'\r'))
        .filter(|line| !line.is_empty());
    let sum = line_sum(lines);
    for threads in THREADS {
        let output = run(&["line-sum", "zero_bytes.bin"], None, threads);
        assert_eq!(output.status.code(), Some(0));
        assert_eq!(String::from_utf8(output.stdout).unwrap(), sum);
    }
}

#[test]
fn a_file_that_cannot_be_read_has_no_line_and_ends_the_run() {
    let first = expected("line_sum_no_final_newline.out");
    for threads in THREADS {
        // The file before it has its line; the file after it is not read.
        let output = run(
            &["line-sum", "records.vcf", "missing.vcf", "output.txt"],
            None,
            threads,
        );
        assert_eq!(output.status.code(), Some(1));
        assert_eq!(output.stdout, first);
        assert_eq!(
            String::from_utf8(output.stderr).unwrap(),
            "ERROR: [vep-file-census] missing.vcf: No such file or directory (os error 2)\n"
        );
        for args in [
            &["line-sum", "tree"][..],
            &["line-sum", "--gunzip", "truncated.vcf.gz"],
            &["line-sum", "--gunzip", "truncated_blocks.vcf.gz"],
            &["line-sum", "--gunzip", "not_gzip.vcf.gz"],
            &["line-sum", "--gunzip", "missing.vcf.gz"],
        ] {
            let output = run(args, None, threads);
            assert_eq!(output.status.code(), Some(1), "{args:?}");
            assert!(output.stdout.is_empty(), "{args:?}");
            assert!(
                output.stderr.starts_with(b"ERROR: [vep-file-census] "),
                "{args:?}"
            );
        }
    }
}

#[test]
fn a_command_line_the_subcommand_cannot_take_is_refused() {
    for args in [
        &["line-sum"][..],
        &["line-sum", "--gunzip"],
        &["line-sum", "--gunzip=yes", "records.vcf"],
        &["line-sum", "--without-lines-starting"],
        &[
            "line-sum",
            "--without-lines-starting",
            "a\nb",
            "records.vcf",
        ],
        &["line-sum", "--algorithm", "sha256", "records.vcf"],
    ] {
        let output = run(args, None, None);
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert!(output.stdout.is_empty(), "{args:?}");
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(
            stderr.starts_with("ERROR: [vep-file-census] ") && stderr.contains("\nusage: "),
            "{args:?}: {stderr}"
        );
    }
    // After `--` a name that starts with a dash is a file.
    let output = run(&["line-sum", "--", "--gunzip"], None, None);
    assert_eq!(output.status.code(), Some(1));
}
