// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line: what is refused with exit code 2, the two ways an option takes its
//! value, names that the digest programs escape, and a gzip stream that ends early.

mod common;

use common::{expected, run, Scratch, BIN};
use sha2::{Digest, Sha256};
use std::process::Command;

fn refused(args: &[&str]) {
    let output = run(args, None, Some("2"));
    assert_eq!(output.status.code(), Some(2), "{args:?}");
    assert!(output.stdout.is_empty(), "{args:?}");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.starts_with("ERROR: [vep-file-census] "),
        "{args:?}: {stderr}"
    );
    assert!(
        stderr.contains("\nusage: vep-file-census count "),
        "{args:?}: {stderr}"
    );
}

#[test]
fn a_command_line_the_tool_cannot_take_is_refused_with_the_usage() {
    for args in [
        &[][..],
        &["tally"],
        &["count", "records.vcf"],
        &["count", "--not-starting", "#"],
        &["count", "--not-starting", "#", "records.vcf", "empty.vcf"],
        &["count", "--not-starting", "#\n", "records.vcf"],
        &["count", "--not-starting"],
        &[
            "count",
            "--gunzip=yes",
            "--not-starting",
            "#",
            "records.vcf",
        ],
        &["count", "--matching", "#", "records.vcf"],
        &["count", "-v", "#", "records.vcf"],
        &["digest", "records.vcf"],
        &["digest", "--algorithm", "sha1", "records.vcf"],
        &["digest", "--algorithm", "sha256"],
        &[
            "digest",
            "--algorithm",
            "sha256",
            "--names",
            "lines",
            "records.vcf",
        ],
        &["digest", "--algorithm", "sha256", "--names", "list"],
        &[
            "digest",
            "--algorithm",
            "sha256",
            "--chunk-bytes",
            "0",
            "records.vcf",
        ],
        &[
            "digest",
            "--algorithm",
            "sha256",
            "--chunk-bytes",
            "many",
            "records.vcf",
        ],
        &[
            "digest",
            "--algorithm",
            "sha256",
            "--gunzip",
            "--chunk-bytes",
            "64",
            "records.vcf.gz",
        ],
        &[
            "digest",
            "--algorithm",
            "sha256",
            "--without-lines-starting",
            "##",
            "--after-leading-lines-starting",
            "##",
            "output.txt",
        ],
        &["strip-chr"],
        &["strip-chr", "--gunzip", "chr.vcf.gz"],
        &["metrics", "output.txt", "report.json"],
        &["contigs"],
        &["slice", "--records", "3", "chr.vcf"],
        &[
            "slice",
            "--records",
            "three",
            "--tally",
            "chr.tally.tsv",
            "chr.vcf",
        ],
        &["fetched", "--list", "fetched/held.tsv"],
        &[
            "fetched",
            "--list",
            "fetched/held.tsv",
            "--root",
            "fetched/root",
            "--jobs",
            "two",
        ],
    ] {
        refused(args);
    }
}

#[test]
fn a_thread_count_that_is_no_whole_number_of_one_or_more_is_refused() {
    for count in ["0", "-1", "four", " 2"] {
        let output = run(
            &["count", "--not-starting", "#", "records.vcf"],
            None,
            Some(count),
        );
        assert_eq!(output.status.code(), Some(2), "{count:?}");
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8_lossy(&output.stderr)
            .contains("VEP_FILE_CENSUS_THREADS takes a whole number of 1 or more"));
    }
}

#[test]
fn help_prints_the_usage_and_leaves_with_zero() {
    for flag in ["--help", "-h"] {
        let output = run(&[flag], None, None);
        assert_eq!(output.status.code(), Some(0));
        let usage = String::from_utf8(output.stdout).unwrap();
        for subcommand in [
            "count",
            "digest",
            "strip-chr",
            "metrics",
            "contigs",
            "slice",
            "fetched",
        ] {
            assert!(
                usage.contains(&format!("vep-file-census {subcommand} ")),
                "{subcommand}"
            );
        }
    }
}

#[test]
fn an_option_takes_its_value_after_an_equals_sign_as_well() {
    let output = run(
        &["count", "--not-starting=#", "records.vcf"],
        None,
        Some("2"),
    );
    assert_eq!(output.stdout, expected("count_records.out"));
    let output = run(
        &["digest", "--algorithm=md5", "--names=read"],
        Some("names_json.txt"),
        Some("2"),
    );
    assert_eq!(output.stdout, expected("digest_names_read.out"));
}

fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

// sha256sum -- 'back\slash' $'new\nline' -dash; xxh128sum 'back\slash'
#[test]
fn a_name_with_a_backslash_or_a_newline_is_escaped_as_sha256sum_and_md5sum_escape_it() {
    let scratch = Scratch::new("escaped-names");
    for name in ["back\\slash", "new\nline", "-dash"] {
        std::fs::write(scratch.path(name), "x\n").unwrap();
    }
    let run_in = |args: &[&str]| {
        let output = Command::new(BIN)
            .current_dir(&scratch.0)
            .env("VEP_FILE_CENSUS_THREADS", "2")
            .args(args)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(0), "{args:?}");
        String::from_utf8(output.stdout).unwrap()
    };
    let hex = sha256_hex(b"x\n");
    assert_eq!(
        run_in(&[
            "digest",
            "--algorithm",
            "sha256",
            "--",
            "back\\slash",
            "new\nline",
            "-dash"
        ]),
        format!("\\{hex}  back\\\\slash\n\\{hex}  new\\nline\n{hex}  -dash\n")
    );
    assert_eq!(
        run_in(&["digest", "--algorithm", "md5", "back\\slash"]),
        "\\401b30e3b8b5d629635a5c613cdb7919  back\\\\slash\n"
    );
    assert_eq!(
        run_in(&["digest", "--algorithm", "xxh128", "back\\slash"]),
        "2528250d92619b8e1bf67208bddf1edb  back\\slash\n"
    );
}

// gzip -dc truncated.vcf.gz | sha256sum: gzip leaves with 1 and sha256sum digests what it wrote.
#[test]
fn a_gzip_stream_that_ends_early_is_digested_to_its_end_and_leaves_with_one() {
    for file in ["truncated.vcf.gz", "truncated_blocks.vcf.gz"] {
        let output = run(
            &["digest", "--algorithm", "sha256", "--gunzip", file],
            None,
            Some("2"),
        );
        assert_eq!(output.status.code(), Some(1), "{file}");
        let line = String::from_utf8(output.stdout).unwrap();
        let (hex, name) = line.split_once("  ").unwrap();
        assert_eq!(name, "-\n");
        assert!(hex.len() == 64 && hex.bytes().all(|b| b.is_ascii_hexdigit()));
        assert!(String::from_utf8_lossy(&output.stderr).starts_with("ERROR: [vep-file-census] "));
    }
}
