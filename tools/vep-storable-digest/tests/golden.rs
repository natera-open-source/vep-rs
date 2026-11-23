// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary against reference bytes.
//!
//! `tests/fixtures/caches/` holds small Storable caches written by perl 5.34.0 with
//! Storable 3.23, and `tests/fixtures/expected/<case>/` the reference run over one of
//! them (`tests/fixtures/caches/README.md` says what each holds), made on the
//! same perl with Data::Dumper 2.179, Digest::MD5 2.58 and IO::Uncompress::Gunzip 2.102:
//! the arguments (`args`, one a line), standard output, standard error, the exit code
//! and the file written (`out.tsv`). Standard error and the exit code are held as the
//! tool gives them, under the rules of the README's "Where the tool does not write what
//! Perl writes". Each test runs the binary with the same arguments in a directory that
//! holds the same cache and asks for the same bytes, on one worker thread and on four.

mod common;

use common::{copy_tree, fixtures, run_in, with_program, Scratch};

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
    let stderr = with_program(std::fs::read(expected.join("stderr")).expect("stderr"));
    let written = std::fs::read(expected.join("out.tsv")).ok();
    for threads in ["1", "4"] {
        let scratch = Scratch::new(&format!("{id}-{threads}"));
        if let Some(cache) = arguments.first() {
            // The cache is the first component of the first argument.
            let top = cache.split('/').next().expect("a path");
            let source = fixtures().join("caches").join(top);
            if source.is_dir() {
                copy_tree(&source, &scratch.path().join(top));
            }
        }
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
        let out = std::fs::read(scratch.path().join("out.tsv")).ok();
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
    // The walk: chromosomes, shard names, a transcript in two shards.
    chromosomes_in_string_order_and_shards_by_name => "plain",
    a_transcript_whose_two_copies_differ_is_warned_of_and_fails_the_run => "conflict",
    a_cache_directory_that_is_absent_writes_an_empty_file => "no-cache",
    arguments_past_the_second_are_not_read => "extra-arguments",
    // The two digests over every kind of value, and how a transcript is named.
    every_kind_of_value_under_both_digests => "values",
    stable_ids_versions_and_coordinates_of_every_kind => "ids",
    // Shard files as gunzip and fd_retrieve take them.
    an_image_that_is_not_compressed_is_read_as_it_stands => "files-plain-image",
    the_first_gzip_member_alone_is_read => "files-two-members",
    bytes_after_the_gzip_member_are_not_read => "files-trailing",
    a_member_without_its_trailer_is_read => "files-no-trailer",
    a_member_with_a_wrong_checksum_is_read => "files-bad-checksum",
    an_image_that_ends_early_counts_and_holds_no_transcript => "files-short-image",
    a_file_that_holds_no_image_stops_the_run => "files-no-image",
    an_empty_file_stops_the_run => "files-empty-file",
    an_empty_gzip_member_stops_the_run => "files-empty-member",
    a_deflate_stream_that_ends_early_stops_the_run => "files-short-gzip",
    // Roots and items of other kinds.
    a_transcript_that_is_an_array_stops_the_run => "files-array-item",
    a_list_that_is_a_hash_stops_the_run => "files-hash-list",
    a_list_that_is_a_string_stops_the_run => "files-string-list",
    a_blessed_hash_at_the_root_stops_the_run => "files-blessed-root",
    a_scalar_at_the_root_stops_the_run => "files-scalar-root",
    // The command line and the output file.
    no_argument_prints_the_usage_line => "usage-none",
    one_argument_prints_the_usage_line => "usage-one",
    an_output_that_cannot_be_opened_is_named => "unwritable",
}

/// Perl would write the reference's address in memory, which no other run repeats.
#[test]
fn a_stable_id_that_is_a_reference_is_refused() {
    let scratch = Scratch::new("id-reference");
    copy_tree(
        &fixtures().join("caches/files"),
        &scratch.path().join("files"),
    );
    let run = run_in(scratch.path(), &["files/id-reference", "out.tsv"], "2");
    assert_eq!(run.status.code(), Some(255));
    assert_eq!(
        lines(&run.stderr),
        "a reference is used as a string: Perl writes its address in memory, which no other run repeats\n"
    );
    assert_eq!(
        std::fs::read(scratch.path().join("out.tsv")).expect("the file"),
        b""
    );
}

/// The number of worker threads is a whole number, 1 or more.
#[test]
fn a_thread_count_that_is_no_number_is_refused() {
    let scratch = Scratch::new("threads");
    let run = run_in(scratch.path(), &["no-such-cache", "out.tsv"], "many");
    assert_eq!(run.status.code(), Some(2));
    assert_eq!(
        lines(&run.stderr),
        "ERROR: [vep-storable-digest] VEP_STORABLE_DIGEST_THREADS=many: expected a whole number of threads, 1 or more\n"
    );
}
