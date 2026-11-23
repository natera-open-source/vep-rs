// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary against reference bytes.
//!
//! The caches are small Storable caches written by perl 5.34.0 with Storable 3.23: the
//! ones of `tools/vep-storable-digest/tests/fixtures/caches/` and three of this tool's
//! own under `tests/fixtures/caches/` (each directory's `README.md` says what its
//! caches hold). `tests/fixtures/lists/` holds lists of ids, and
//! `tests/fixtures/expected/<case>/` the reference run over a cache and a list, made on
//! the same perl with Data::Dumper 2.179 and IO::Uncompress::Gunzip 2.102: the
//! arguments (`args`, one a line), standard output, standard error, the exit code and
//! the file written (`out.txt`). Standard error and the exit code are held as the tool
//! gives them, under the rules of the README's "Where the tool does not write what Perl
//! writes". Each test runs the binary with the same arguments in a directory that holds
//! the same cache and lists and asks for the same bytes, on one worker thread and on
//! four.

mod common;

use common::{copy_tree, fixtures, run_in, with_program, Scratch};

fn lines(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// A fresh directory with the lists and the cache `top`, from this tool's fixtures or
/// from the digest tool's.
fn scratch_with(name: &str, top: Option<&str>) -> Scratch {
    let scratch = Scratch::new(name);
    copy_tree(&fixtures().join("lists"), &scratch.path().join("lists"));
    if let Some(top) = top {
        let own = fixtures().join("caches").join(top);
        let shared = fixtures()
            .join("../../../vep-storable-digest/tests/fixtures/caches")
            .join(top);
        if own.is_dir() {
            copy_tree(&own, &scratch.path().join(top));
        } else if shared.is_dir() {
            copy_tree(&shared, &scratch.path().join(top));
        }
    }
    scratch
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
    let stderr = with_program(std::fs::read(expected.join("stderr")).expect("stderr"));
    let written = std::fs::read(expected.join("out.txt")).ok();
    // The cache is the first component of the first argument.
    let top = arguments
        .first()
        .map(|cache| cache.split('/').next().expect("a path"));
    for threads in ["1", "4"] {
        let scratch = scratch_with(&format!("{id}-{threads}"), top);
        let run = run_in(scratch.path(), &argument_refs, threads);
        assert_eq!(
            run.status.code(),
            Some(exit),
            "{id}: exit code; stderr: {}",
            lines(&run.stderr)
        );
        assert_eq!(lines(&run.stdout), lines(&stdout), "{id}: standard output");
        assert_eq!(lines(&run.stderr), lines(&stderr), "{id}: standard error");
        assert_eq!(run.stderr, stderr, "{id}: standard error bytes");
        let out = std::fs::read(scratch.path().join("out.txt")).ok();
        assert_eq!(
            out.as_deref().map(lines),
            written.as_deref().map(lines),
            "{id}: the file written"
        );
        assert_eq!(out, written, "{id}: the bytes of the file written");
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
    // The walk and the list of ids.
    chromosomes_in_their_own_order_and_an_id_once => "walk",
    an_empty_list_dumps_nothing => "empty-list",
    a_list_without_an_id_dumps_nothing => "no-id-in-list",
    a_transcript_in_two_shards_is_dumped_once => "conflict",
    a_cache_directory_that_is_absent_writes_an_empty_file => "no-cache",
    arguments_past_the_third_are_not_read => "extra-arguments",
    // What a block holds.
    every_kind_of_value_at_one_indentation => "values",
    prediction_matrices_are_summarised => "predictions",
    a_matrix_that_is_a_list_stops_the_run => "predictions-list",
    ids_of_every_kind_against_a_list_of_bytes => "ids",
    // Items that are no transcripts.
    a_string_among_the_transcripts_stops_the_run => "plain",
    an_undef_item_is_passed_over_and_a_string_stops_the_run => "items",
    // Shard files as gunzip and fd_retrieve take them.
    an_image_that_is_not_compressed_is_read_as_it_stands => "files-plain-image",
    the_first_gzip_member_alone_is_read => "files-two-members",
    bytes_after_the_gzip_member_are_not_read => "files-trailing",
    a_member_without_its_trailer_is_read => "files-no-trailer",
    a_member_with_a_wrong_checksum_is_read => "files-bad-checksum",
    an_image_that_ends_early_holds_no_transcript => "files-short-image",
    a_file_that_holds_no_image_stops_the_run => "files-no-image",
    an_empty_file_stops_the_run => "files-empty-file",
    an_empty_gzip_member_stops_the_run => "files-empty-member",
    a_deflate_stream_that_ends_early_stops_the_run => "files-short-gzip",
    a_transcript_that_is_an_array_stops_the_run => "files-array-item",
    a_list_that_is_a_hash_stops_the_run => "files-hash-list",
    a_list_that_is_a_string_stops_the_run => "files-string-list",
    a_blessed_hash_at_the_root_stops_the_run => "files-blessed-root",
    a_scalar_at_the_root_stops_the_run => "files-scalar-root",
    a_stable_id_that_is_a_reference_names_no_transcript => "files-id-reference",
    // The command line and the two files.
    no_argument_prints_the_usage_line => "usage-none",
    two_arguments_print_the_usage_line => "usage-two",
    a_list_that_cannot_be_opened_is_named => "unreadable-list",
    an_output_that_cannot_be_opened_is_named => "unwritable",
}

/// Perl would write the length of the reference's address in memory.
#[test]
fn a_matrix_that_is_a_reference_where_a_length_is_asked_is_refused() {
    let scratch = scratch_with("matrix-reference", Some("predictions"));
    let run = run_in(
        scratch.path(),
        &["predictions", "lists/predictions-reference.txt", "out.txt"],
        "2",
    );
    assert_eq!(run.status.code(), Some(255));
    assert_eq!(
        lines(&run.stderr),
        "a reference is used as a string: Perl writes its address in memory, which no other run repeats\n"
    );
    assert_eq!(
        std::fs::read(scratch.path().join("out.txt")).expect("the file"),
        b""
    );
}

/// The number of worker threads is a whole number, 1 or more.
#[test]
fn a_thread_count_that_is_no_number_is_refused() {
    let scratch = scratch_with("threads", None);
    let run = run_in(
        scratch.path(),
        &["no-such-cache", "lists/one.txt", "out.txt"],
        "0",
    );
    assert_eq!(run.status.code(), Some(2));
    assert_eq!(
        lines(&run.stderr),
        "ERROR: [vep-storable-dump] VEP_STORABLE_DUMP_THREADS=0: expected a whole number of threads, 1 or more\n"
    );
}
