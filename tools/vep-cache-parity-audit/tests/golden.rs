// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary against reference bytes.
//!
//! `tests/fixtures/storable/` holds small Storable caches written by perl 5.34.0 with
//! Storable 3.23, `tests/fixtures/json/` the JSON caches `tools/vep-storable-to-json`
//! writes from them and variants of those (`tests/fixtures/README.md` says what each
//! holds), and `tests/fixtures/expected/<case>/` the reference run over a pair, made on
//! the same perl with JSON::PP 4.06 and File::Find 1.39: the arguments (`args`, one a
//! line), standard output, standard error and the exit code. The three are held as the
//! tool gives them, under the rules of the README's "Where the tool does not write what
//! Perl writes": the manifest with the keys of every object in the order of their
//! names, standard error without Perl's own lines, and exit code 255 for a run that
//! died. The manifest of a run over several shards is held with the number of the ids
//! of each gap in the place of the ids. Each test runs the binary with the same
//! arguments in a directory that holds the same caches and asks for the same bytes, on
//! one worker thread and on four.

mod common;

use common::{copy_tree, fixtures, run_in, with_program, Scratch};

fn lines(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// A fresh directory with the Storable and the JSON fixture caches, and the cache of
/// odd shard files the digest tool's fixtures hold.
fn scratch(name: &str) -> Scratch {
    let scratch = Scratch::new(name);
    copy_tree(
        &fixtures().join("storable"),
        &scratch.path().join("storable"),
    );
    copy_tree(&fixtures().join("json"), &scratch.path().join("json"));
    copy_tree(
        &fixtures().join("../../../vep-storable-digest/tests/fixtures/caches/files"),
        &scratch.path().join("storable/files"),
    );
    copy_tree(
        &fixtures().join("../../../vep-storable-digest/tests/fixtures/caches/values"),
        &scratch.path().join("storable/values"),
    );
    scratch
}

/// The manifest with the number of ids of each gap in the place of the ids. The ids a
/// gap is exemplified by are the first three in the sample's order, and over several
/// shards that order follows the order the system lists the cache's directories in.
fn without_example_ids(manifest: &[u8]) -> Vec<u8> {
    let open = b"\"parity_gap_examples\":{";
    let Some(start) = manifest
        .windows(open.len())
        .position(|window| window == open)
    else {
        return manifest.to_vec();
    };
    let mut out = manifest[..start + open.len()].to_vec();
    let mut at = start + open.len();
    let mut in_string = false;
    while at < manifest.len() {
        let byte = manifest[at];
        if in_string {
            out.push(byte);
            if byte == b'\\' {
                out.push(manifest[at + 1]);
                at += 1;
            } else if byte == b'"' {
                in_string = false;
            }
            at += 1;
            continue;
        }
        match byte {
            b'"' => {
                in_string = true;
                out.push(byte);
                at += 1;
            }
            b'[' => {
                // A list of strings and numbers: its items are one more than its commas
                // outside the strings.
                let mut ids = 0;
                let mut quoted = false;
                at += 1;
                if manifest[at] != b']' {
                    ids = 1;
                }
                while manifest[at] != b']' || quoted {
                    match manifest[at] {
                        b'\\' if quoted => at += 1,
                        b'"' => quoted = !quoted,
                        b',' if !quoted => ids += 1,
                        _ => {}
                    }
                    at += 1;
                }
                at += 1;
                out.extend_from_slice(ids.to_string().as_bytes());
            }
            b'}' => {
                out.extend_from_slice(&manifest[at..]);
                break;
            }
            _ => {
                out.push(byte);
                at += 1;
            }
        }
    }
    out
}

fn case(id: &str, several_shards: bool) {
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
    let stderr = with_program(std::fs::read(expected.join("stderr")).expect("stderr"));
    for threads in ["1", "4"] {
        let scratch = scratch(&format!("{id}-{threads}"));
        let run = run_in(scratch.path(), &argument_refs, threads);
        assert_eq!(
            run.status.code(),
            Some(exit),
            "{id}: exit code; stderr: {}",
            lines(&run.stderr)
        );
        let manifest = if several_shards {
            without_example_ids(&run.stdout)
        } else {
            run.stdout.clone()
        };
        assert_eq!(lines(&manifest), lines(&stdout), "{id}: standard output");
        assert_eq!(manifest, stdout, "{id}: standard output bytes");
        assert_eq!(lines(&run.stderr), lines(&stderr), "{id}: standard error");
        assert_eq!(run.stderr, stderr, "{id}: standard error bytes");
    }
}

macro_rules! cases {
    ($($name:ident => $id:literal, $several_shards:literal;)*) => {
        $(
            #[test]
            fn $name() {
                case($id, $several_shards);
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

// Each case with whether its cache holds several shards.
cases! {
    // Every class of path, on the converter's own output and on edited copies of it.
    a_converted_cache_audits_with_its_fixture_fields_as_gaps => "clean", false;
    fields_taken_out_added_and_changed_in_kind => "gaps", false;
    files_that_are_no_shards_and_ids_of_every_kind => "junk", false;
    // The walk of a whole cache: shards, files that are none, chromosomes of both classes.
    a_cache_of_several_shards_from_its_parent_directory => "walk", true;
    a_cache_of_several_shards_from_its_own_directory => "walk-root", true;
    a_root_that_is_one_shard_file => "file-root", false;
    a_json_root_that_is_one_shard_file => "json-file-root", false;
    a_json_root_that_is_absent_leaves_every_transcript_unmatched => "json-missing", false;
    // The sample and its seed.
    the_seed_of_a_run_without_one_is_the_number_one => "seed-default", false;
    a_seed_that_is_given_is_a_string => "seed-given", false;
    a_seed_that_is_no_number_is_a_string => "seed-text", false;
    two_transcripts_are_sampled => "sample-two", false;
    no_transcript_is_sampled => "sample-zero", false;
    a_sample_that_is_no_number_is_none => "sample-text", false;
    a_sample_of_a_fraction_is_rounded_up => "sample-fraction", false;
    a_negative_sample_is_none => "sample-negative", false;
    arguments_past_the_fourth_are_not_read => "extra-arguments", false;
    // Runs that die.
    a_storable_root_that_is_absent_stops_the_run => "storable-missing", false;
    a_storable_root_without_a_shard_stops_the_run => "storable-without-shards", false;
    a_list_item_that_is_an_array_stops_the_run => "array-item", false;
    no_argument_prints_the_usage_line => "usage-none", false;
    one_argument_prints_the_usage_line => "usage-one", false;
    a_first_argument_of_zero_prints_the_usage_line => "usage-zero", false;
    an_empty_second_argument_prints_the_usage_line => "usage-empty", false;
}

/// The same walk in Perl has no end on an array that holds itself.
#[test]
fn an_array_that_holds_itself_is_refused() {
    let scratch = scratch("array-cycle");
    let run = run_in(scratch.path(), &["storable/values", "json/tiny"], "2");
    assert_eq!(run.status.code(), Some(255));
    assert_eq!(
        lines(&run.stderr),
        "a transcript holds an array that holds itself: Perl walks it without end\n"
    );
    assert!(run.stdout.is_empty());
}

#[test]
fn the_ids_of_a_gap_give_way_to_their_number() {
    let manifest = br#"{"matched_in_json":2,"parity_gap_examples":{"a[]":["x","y,z","q\"]"],"b":[12345],"c":[]},"parity_gaps":{"a[]":3}}"#;
    assert_eq!(
        lines(&without_example_ids(manifest)),
        r#"{"matched_in_json":2,"parity_gap_examples":{"a[]":3,"b":1,"c":0},"parity_gaps":{"a[]":3}}"#
    );
}
