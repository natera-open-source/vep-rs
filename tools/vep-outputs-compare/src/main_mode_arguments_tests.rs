// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line: which mode an argument list selects and which arguments each mode needs.

use super::*;

fn parsed(args: &[&str]) -> Args {
    match parse_args(args.iter().map(OsString::from).collect()) {
        Ok(args) => args,
        Err(message) => panic!("{message}"),
    }
}

/// The usage error an argument list is refused with.
fn refused(args: &[&str]) -> String {
    match parse_args(args.iter().map(OsString::from).collect()) {
        Ok(_) => panic!("accepted: {args:?}"),
        Err(message) => message,
    }
}

/// The message naming the arguments the mode of `args` lacks.
fn lacking(args: &[&str]) -> String {
    require_mode_args(&parsed(args)).expect_err("a missing argument")
}

#[test]
fn every_option_stands_once_and_the_cache_option_is_read_by_any_spelling() {
    let mut names: Vec<&str> = OPTIONS.iter().map(|(name, _)| *name).collect();
    names.sort_unstable();
    names.dedup();
    assert_eq!(names.len(), OPTIONS.len(), "an option is declared twice");
    let scratch = TempDir::create(&std::env::temp_dir(), "cache_option_test_").unwrap();
    let cache = scratch.path.to_str().unwrap();
    for spelling in [
        vec!["--json-cache", cache],
        vec!["--json", cache],
        vec![&format!("--json-cache={cache}")[..]],
    ] {
        assert_eq!(
            parsed(&spelling).json_cache.as_deref(),
            Some(&*scratch.path)
        );
    }
    assert_eq!(parsed(&[]).json_cache, None);
}

#[test]
fn an_argument_list_without_a_mode_flag_is_a_compare_run() {
    let args = parsed(&["--perl-dir", "p", "--rust-dir", "r", "--report-dir", "o"]);
    assert_eq!(args.mode.name(), "compare");
    assert_eq!(require_mode_args(&args), Ok(()));
}

#[test]
fn the_by_consequence_class_flag_selects_its_mode() {
    let args = parsed(&[
        "--by-consequence-class",
        "--perl-dir",
        "p",
        "--discordant-in",
        "d.tsv",
        "--dataset",
        "s01",
        "--class-csv",
        "c.csv",
        "--sweep-id",
        "sweep-20260101T000000Z",
        "--binary-md5",
        "deadbeef",
    ]);
    assert_eq!(args.mode.name(), "by-consequence-class");
    assert_eq!(require_mode_args(&args), Ok(()));
}

#[test]
fn the_per_class_mode_needs_the_sweep_and_the_binary_pin() {
    let message = lacking(&[
        "--by-consequence-class",
        "--perl-dir",
        "p",
        "--discordant-in",
        "d.tsv",
        "--dataset",
        "s01",
        "--class-csv",
        "c.csv",
    ]);
    assert!(message.contains("--sweep-id"), "{message}");
    assert!(message.contains("--binary-md5"), "{message}");
}

#[test]
fn the_per_class_mode_never_asks_for_the_vep_rs_output_directory() {
    let message = lacking(&["--by-consequence-class"]);
    assert_eq!(
        message,
        "error: mode by-consequence-class requires --perl-dir, --discordant-in, --dataset, \
         --class-csv, --sweep-id, --binary-md5"
    );
    assert!(!message.contains("--rust-dir"));
}

#[test]
fn the_compare_mode_needs_its_three_directories() {
    let message = lacking(&["--perl-dir", "p"]);
    assert!(message.contains("--rust-dir"), "{message}");
}

#[test]
fn the_pool_mode_needs_the_two_csv_paths_alone() {
    let args = parsed(&[
        "--pool-consequence-class",
        "--class-csv",
        "c.csv",
        "--pooled-csv",
        "p.csv",
        "--per-term",
        "shapes/per_term.tsv",
    ]);
    assert_eq!(args.mode.name(), "pool-consequence-class");
    assert_eq!(args.per_term, Some(PathBuf::from("shapes/per_term.tsv")));
    assert_eq!(args.per_term_engine, "vep-rs");
    assert_eq!(require_mode_args(&args), Ok(()));
    assert_eq!(
        lacking(&["--pool-consequence-class", "--class-csv", "c.csv"]),
        "error: mode pool-consequence-class requires --pooled-csv"
    );
    assert_eq!(
        lacking(&["--pool-consequence-class"]),
        "error: mode pool-consequence-class requires --class-csv, --pooled-csv"
    );
}

#[test]
fn the_two_mode_flags_are_not_taken_together() {
    assert_eq!(
        refused(&["--by-consequence-class", "--pool-consequence-class"]),
        "argument --pool-consequence-class: not allowed with argument --by-consequence-class"
    );
}

#[test]
fn the_engine_label_is_one_of_the_two_spellings_the_csvs_carry() {
    let message = refused(&["--engine", "vep_rs"]);
    assert!(
        message.starts_with("argument --engine: invalid choice: 'vep_rs'"),
        "{message}"
    );
    assert_eq!(parsed(&["--engine", "vep-rs"]).engine, "vep-rs");
}

#[test]
fn the_tuple_threshold_of_an_empty_argument_list_is_the_fixed_one() {
    assert_eq!(parsed(&[]).min_tuples, classes::PER_CLASS_MIN_TUPLES);
}
