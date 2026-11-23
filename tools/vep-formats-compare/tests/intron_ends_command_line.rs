// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two intron classes through the command line: each needs the candidate's JSON cache.

mod compare_common;
mod intron_ends_common;
mod support;

use std::ffi::OsString;
use std::path::Path;

use intron_ends_common::input_vcf;
use intron_ends_common::intron_end_cache;
use intron_ends_common::option_help;
use intron_ends_common::row;
use intron_ends_common::Row;
use intron_ends_common::INTRON_END;
use intron_ends_common::INTRON_START;
use support::compact;
use support::member;
use support::report;
use support::run;
use support::scratch;
use vep_formats_compare::cli;

fn text(path: &Path) -> String {
    path.to_str().expect("a path of text").to_string()
}

/// The row's default-format and VCF pair and the VCF's input under `directory`: the
/// arguments that compare the default-format pair, and those that compare the VCF pair
/// against the input.
fn planted_pairs(directory: &Path, row: &Row) -> (Vec<String>, Vec<String>) {
    let mut compared = Vec::new();
    for fmt in ["default", "vcf"] {
        let (reference, candidate, inputs) = row.text(fmt);
        let reference_path = directory.join(format!("reference.{fmt}"));
        let candidate_path = directory.join(format!("candidate.{fmt}"));
        std::fs::write(&reference_path, reference).unwrap();
        std::fs::write(&candidate_path, candidate).unwrap();
        let mut argv = vec![
            "--format".to_string(),
            fmt.to_string(),
            "--expected".to_string(),
            text(&reference_path),
            "--actual".to_string(),
            text(&candidate_path),
        ];
        if let Some(inputs) = inputs {
            let input = input_vcf(&directory.join("in.vcf"), &inputs);
            argv.extend(["--input".to_string(), text(&input)]);
        }
        compared.push(argv);
    }
    let vcf = compared.pop().unwrap();
    (compared.pop().unwrap(), vcf)
}

fn with(argv: &[String], more: &[&str]) -> Vec<OsString> {
    argv.iter()
        .map(String::as_str)
        .chain(more.iter().copied())
        .map(OsString::from)
        .collect()
}

/// What the report at `path` counts as applied, as compact JSON.
fn applied(path: &Path) -> String {
    compact(member(&report(path), &["masks", "applied"]))
}

#[test]
fn the_intron_end_class_needs_the_candidate_json_cache() {
    let directory = scratch("intron-end-class-needs-the-cache");
    let cache = text(&intron_end_cache(&directory.join("cache")));
    let none = text(&directory.join("none"));
    let (default, vcf) = planted_pairs(
        &directory,
        &row("ENST00000000001", 11199, "G", "G", "c.51-1dup"),
    );
    let argv = with(&default, &["--masks", INTRON_END]);
    let args = cli::parse_args(&cli::build_parser(), &argv).expect("the command line parses");
    let refused = cli::run(&args, None).unwrap_err().to_string();
    assert!(refused.contains("needs --json-cache"), "{refused}");
    let ran = run(&directory, &argv);
    assert_eq!(ran.code, 2);
    assert_eq!(
        ran.err(),
        "ERROR: [compare_vep_formats] --masks hgvsp_intron_end_insertion needs --json-cache, the candidate JSON cache\n"
    );
    let ran = run(
        &directory,
        &with(&default, &["--masks", INTRON_END, "--json-cache", &none]),
    );
    assert_eq!(ran.code, 2);
    assert_eq!(
        ran.err(),
        format!(
            "ERROR: [compare_vep_formats] --json-cache: {none} has no transcripts/ directory\n"
        )
    );
    let default_report = directory.join("default.json");
    let ran = run(
        &directory,
        &with(
            &default,
            &[
                "--masks",
                INTRON_END,
                "--json-cache",
                &cache,
                "--report",
                &text(&default_report),
            ],
        ),
    );
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(
        applied(&default_report),
        format!(r#"{{"{INTRON_END}": 1}}"#)
    );
    let vcf_report = directory.join("vcf.json");
    let ran = run(
        &directory,
        &with(
            &vcf,
            &[
                "--masks",
                INTRON_END,
                "--json-cache",
                &cache,
                "--report",
                &text(&vcf_report),
            ],
        ),
    );
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(applied(&vcf_report), format!(r#"{{"{INTRON_END}": 1}}"#));
    assert!(option_help(cli::build_parser().format_help(), "--json-cache").contains(INTRON_END));
}

#[test]
fn the_intron_start_class_needs_the_candidate_json_cache() {
    let directory = scratch("intron-start-class-needs-the-cache");
    let cache = text(&intron_end_cache(&directory.join("cache")));
    let (default, vcf) = planted_pairs(
        &directory,
        &row("ENST00000000001", 11097, "G", "TAT", "c.50_50+1insTAT").both("HGVS_OFFSET", "2"),
    );
    let ran = run(&directory, &with(&default, &["--masks", INTRON_START]));
    assert_eq!(ran.code, 2);
    assert_eq!(
        ran.err(),
        "ERROR: [compare_vep_formats] --masks hgvsp_intron_start_insertion needs --json-cache, the candidate JSON cache\n"
    );
    let both = format!("{INTRON_START},{INTRON_END}");
    let ran = run(&directory, &with(&default, &["--masks", &both]));
    assert_eq!(ran.code, 2);
    assert_eq!(
        ran.err(),
        "ERROR: [compare_vep_formats] --masks hgvsp_intron_end_insertion,hgvsp_intron_start_insertion needs --json-cache, the candidate JSON cache\n"
    );
    let default_report = directory.join("default.json");
    let ran = run(
        &directory,
        &with(
            &default,
            &[
                "--masks",
                INTRON_START,
                "--json-cache",
                &cache,
                "--report",
                &text(&default_report),
            ],
        ),
    );
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(
        applied(&default_report),
        format!(r#"{{"{INTRON_START}": 1}}"#)
    );
    let vcf_report = directory.join("vcf.json");
    let ran = run(
        &directory,
        &with(
            &vcf,
            &[
                "--masks",
                INTRON_START,
                "--json-cache",
                &cache,
                "--report",
                &text(&vcf_report),
            ],
        ),
    );
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(applied(&vcf_report), format!(r#"{{"{INTRON_START}": 1}}"#));
    assert!(option_help(cli::build_parser().format_help(), "--json-cache").contains(INTRON_START));
}

/// Without the cache a command that names one of the two classes is refused before any row
/// is compared, and a directory without `transcripts/` is refused by its name.
#[test]
fn a_command_that_names_an_intron_class_without_its_cache_is_refused() {
    let directory = scratch("intron-classes-refused-without-the-cache");
    let none = text(&directory.join("none"));
    let (default, _) = planted_pairs(
        &directory,
        &row("ENST00000000001", 11199, "G", "G", "c.51-1dup"),
    );
    for (masks, named) in [
        (INTRON_END, INTRON_END.to_string()),
        (INTRON_START, INTRON_START.to_string()),
        (
            "hgvsp_intron_start_insertion,hgvsp_intron_end_insertion",
            format!("{INTRON_END},{INTRON_START}"),
        ),
    ] {
        let ran = run(&directory, &with(&default, &["--masks", masks]));
        assert_eq!(ran.code, 2);
        assert_eq!(
            ran.err(),
            format!(
                "ERROR: [compare_vep_formats] --masks {named} needs --json-cache, the candidate JSON cache\n"
            )
        );
        let ran = run(
            &directory,
            &with(&default, &["--masks", masks, "--json-cache", &none]),
        );
        assert_eq!(ran.code, 2);
        assert_eq!(
            ran.err(),
            format!(
                "ERROR: [compare_vep_formats] --json-cache: {none} has no transcripts/ directory\n"
            )
        );
    }
    let help = option_help(cli::build_parser().format_help(), "--json-cache");
    assert!(help.contains(INTRON_END) && help.contains(INTRON_START));
}
