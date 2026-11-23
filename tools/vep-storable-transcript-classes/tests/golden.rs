// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary against reference bytes.
//!
//! `tests/fixtures/inputs/` holds pairs of dumps and JSON caches, and
//! `tests/fixtures/expected/<case>/` one reference run over them, made by CPython 3.14:
//! the arguments (`args`, one a line), standard output, standard error, the exit code,
//! and under `tree/` every file the run wrote. A file the run names with a character
//! outside printable ASCII, a double quote or a backslash lies under `files/` by a
//! plain name, and `files.tsv` beside it spells the path the run wrote it to. Each test
//! runs the binary with the same arguments in a copy of the inputs, on one worker
//! thread and on four, and asks for the same bytes: the two streams, the exit code, and
//! the directory as the run leaves it. Two things are held as the tool gives them (the
//! README's "Outside the contract"): the order of the summary's `per_class` keys, and
//! standard error where the reference is `argparse`'s own text or a traceback.

mod common;

use common::{copy_tree, fixtures, run_in, tree, written, Scratch};

fn lines(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn case(id: &str) {
    let expected = fixtures().join("expected").join(id);
    let arguments: Vec<String> = std::fs::read_to_string(expected.join("args"))
        .expect("args")
        .lines()
        .map(str::to_string)
        .collect();
    let argument_refs: Vec<&str> = arguments.iter().map(String::as_str).collect();
    let exit: i32 = std::fs::read_to_string(expected.join("exit"))
        .expect("exit")
        .trim()
        .parse()
        .expect("an exit code");
    let stdout = std::fs::read(expected.join("stdout")).expect("stdout");
    let stderr = std::fs::read(expected.join("stderr")).expect("stderr");
    // The directory after the run: the inputs and what the reference run wrote.
    let mut wanted = tree(&fixtures().join("inputs"));
    wanted.extend(written(&expected));
    for threads in ["1", "4"] {
        let scratch = Scratch::new(&format!("{id}-{threads}"));
        copy_tree(&fixtures().join("inputs"), scratch.path());
        let run = run_in(scratch.path(), &argument_refs, threads);
        assert_eq!(
            run.status.code(),
            Some(exit),
            "{id}: exit code; stderr: {}",
            lines(&run.stderr)
        );
        assert_eq!(lines(&run.stdout), lines(&stdout), "{id}: standard output");
        assert_eq!(run.stdout, stdout, "{id}: standard output bytes");
        assert_eq!(lines(&run.stderr), lines(&stderr), "{id}: standard error");
        assert_eq!(run.stderr, stderr, "{id}: standard error bytes");
        let left = tree(scratch.path());
        let names = |files: &std::collections::BTreeMap<String, Vec<u8>>| {
            files.keys().cloned().collect::<Vec<_>>()
        };
        assert_eq!(names(&left), names(&wanted), "{id}: files after the run");
        for (name, bytes) in &wanted {
            assert_eq!(
                lines(&left[name]),
                lines(bytes),
                "{id}: {name} after the run"
            );
            assert!(left[name] == *bytes, "{id}: the bytes of {name}");
        }
    }
}

macro_rules! cases {
    ($($name:ident => $id:literal,)*) => {
        $(
            #[test]
            fn $name() {
                case($id);
            }
        )*

        #[test]
        fn every_expected_directory_has_its_test() {
            let mut listed: Vec<&str> = vec![$($id),*];
            listed.sort_unstable();
            let mut present: Vec<String> = std::fs::read_dir(fixtures().join("expected"))
                .expect("expected directory")
                .map(|entry| entry.expect("entry").file_name().to_string_lossy().into_owned())
                .collect();
            present.sort();
            assert_eq!(present, listed);
        }
    };
}

cases! {
    // The classes of two dumps.
    every_class_and_combination_of_two_dumps => "plain",
    the_output_directory_is_made_with_its_parents => "nested-out-dir",
    one_dump_against_itself_is_identical_throughout => "same-dump",
    the_two_dumps_the_other_way_round => "swapped",
    bytes_that_are_not_utf8_and_blocks_of_odd_shape => "odd-bytes",
    dumps_without_a_block_class_nothing => "empty-dumps",
    ids_in_one_dump_only_are_named_and_fail_the_run => "one-sided",
    // The stored translation of a JSON cache.
    each_transcript_is_matched_against_a_json_cache => "json-cache",
    a_label_names_the_json_cache => "json-cache-label",
    an_empty_label_is_the_directory_name => "json-cache-empty-label",
    a_json_cache_named_by_a_dot_has_an_empty_name => "json-cache-here",
    a_json_cache_without_transcripts_holds_no_id => "json-cache-without-transcripts",
    a_json_cache_that_is_absent_holds_no_id => "json-cache-absent",
    a_cache_no_object_under_an_id_read_before_is_passed_over => "json-cache-late-bad-cache",
    versions_named_like_the_fixed_matches_share_their_counts => "versions-named-like-matches",
    // The command line.
    options_shortened_and_with_an_equals_sign => "prefixes",
    options_that_are_required_are_named => "missing-required",
    an_option_that_is_unknown_is_named => "unknown-option",
    a_prefix_of_two_options_is_named => "ambiguous-option",
    an_option_without_its_value_is_named => "no-value",
    // Refusals in the tool's own words.
    a_dump_that_is_absent => "tb-missing-dump",
    a_block_whose_first_line_has_no_end => "tb-no-line-end",
    an_output_directory_that_is_a_file => "tb-out-dir-is-a-file",
    a_shard_that_is_no_json => "tb-bad-json",
    a_stable_id_that_is_a_number => "tb-bad-id",
    a_cache_that_is_text_under_a_wanted_id => "tb-bad-cache",
}

#[test]
fn help_goes_to_standard_output() {
    let scratch = Scratch::new("help");
    for arguments in [&["-h"][..], &["--help"], &["--version-a", "1", "--he"]] {
        let run = run_in(scratch.path(), arguments, "1");
        assert_eq!(run.status.code(), Some(0), "{arguments:?}");
        assert!(lines(&run.stdout)
            .starts_with("usage: vep-storable-transcript-classes [-h] --version-a"));
        assert!(run.stderr.is_empty());
    }
}
