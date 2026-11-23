// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `sequence_regions.tsv` as the binary writes it for caches of `tests/fixtures/caches/`:
//! its rows and their order, the header alone for a cache without slices, and no file
//! after a run that fails.

mod common;

use std::process::Output;

use common::{copy_tree, fixtures, run_in, Scratch};

const HEADER: &[u8] = b"name\tlength\tcircular\n";

/// Converts a copy of the fixture cache `cache` into `<cache>.out` in a scratch directory
/// named for the test, where `stale`, when given, stands as the `sequence_regions.tsv` of
/// a run before this one.
fn convert(test: &str, cache: &str, stale: Option<&[u8]>) -> (Scratch, Output) {
    let scratch = Scratch::new(test);
    copy_tree(
        &fixtures().join("caches").join(cache),
        &scratch.path().join(cache),
    );
    let output = format!("{cache}.out");
    std::fs::create_dir_all(scratch.path().join(&output)).expect("output directory");
    if let Some(stale) = stale {
        std::fs::write(
            scratch.path().join(&output).join("sequence_regions.tsv"),
            stale,
        )
        .expect("the stale file");
    }
    let run = run_in(scratch.path(), &[cache, &output]);
    (scratch, run)
}

fn count(haystack: &[u8], needle: &str) -> usize {
    haystack
        .windows(needle.len())
        .filter(|window| *window == needle.as_bytes())
        .count()
}

#[test]
fn a_linear_a_circular_and_an_unstated_sequence_are_listed_in_contig_order() {
    // Chromosomes 1, 2, MT and scaffold_a; the transcripts of 2 carry no slice.
    let (scratch, run) = convert("rows", "regions", None);
    assert_eq!(run.status.code(), Some(0));
    let written = std::fs::read(scratch.path().join("regions.out/sequence_regions.tsv"))
        .expect("sequence_regions.tsv");
    assert_eq!(
        written,
        [HEADER, b"1\t5000\t0\nMT\t16000\t1\nscaffold_a\t700\t-\n"].concat()
    );
    assert_eq!(
        count(
            &run.stdout,
            "Wrote regions.out/sequence_regions.tsv (3 sequences)\n"
        ),
        1
    );
}

#[test]
fn rows_follow_contig_order_where_the_directory_names_sort_otherwise() {
    // Chromosome 1 stands a second time as 10, and MT as 9 and as X: the names sort 1, 10,
    // 2, 9, MT, X, scaffold_a, and the contig order is 1, 2, 9, 10, X, MT, scaffold_a.
    let scratch = Scratch::new("order");
    let cache = scratch.path().join("regions");
    copy_tree(&fixtures().join("caches/regions"), &cache);
    for (from, to) in [("1", "10"), ("MT", "9"), ("MT", "X")] {
        copy_tree(&cache.join(from), &cache.join(to));
    }
    std::fs::create_dir_all(scratch.path().join("regions.out")).expect("output directory");
    let run = run_in(scratch.path(), &["regions", "regions.out"]);
    assert_eq!(run.status.code(), Some(0));
    let written = std::fs::read(scratch.path().join("regions.out/sequence_regions.tsv"))
        .expect("sequence_regions.tsv");
    assert_eq!(
        written,
        [
            HEADER,
            b"1\t5000\t0\n9\t16000\t1\n10\t5000\t0\nX\t16000\t1\nMT\t16000\t1\nscaffold_a\t700\t-\n"
        ]
        .concat()
    );
    assert_eq!(
        count(
            &run.stdout,
            "Wrote regions.out/sequence_regions.tsv (6 sequences)\n"
        ),
        1
    );
}

#[test]
fn contig_order_keeps_its_names_as_strings_and_no_transcript_gains_a_slice() {
    let (scratch, run) = convert("records", "regions", None);
    assert_eq!(run.status.code(), Some(0));
    let info = std::fs::read(scratch.path().join("regions.out/info.json")).expect("info.json");
    assert_eq!(
        count(
            &info,
            "\"contig_order\" : [\n      \"1\",\n      \"2\",\n      \"MT\",\n      \"scaffold_a\"\n   ]"
        ),
        1
    );
    let shard = std::fs::read(
        scratch
            .path()
            .join("regions.out/transcripts/MT/1-1000000.json"),
    )
    .expect("the shard of MT");
    assert_eq!(count(&shard, "\"stable_id\":"), 2);
    assert_eq!(count(&shard, "slice"), 0);
}

#[test]
fn a_cache_whose_transcripts_carry_no_slice_gets_the_header_line_alone() {
    let (scratch, run) = convert("header", "both-versions", None);
    assert_eq!(run.status.code(), Some(0));
    let written = std::fs::read(
        scratch
            .path()
            .join("both-versions.out/sequence_regions.tsv"),
    )
    .expect("sequence_regions.tsv");
    assert_eq!(written, HEADER);
}

#[test]
fn a_failed_run_leaves_no_file_that_stood_before_it() {
    // One shard of chromosome 1 is no Storable image.
    let stale = [HEADER, b"stale\t1\t0\n"].concat();
    let (scratch, run) = convert("failed", "regions-corrupt", Some(&stale));
    assert_eq!(run.status.code(), Some(1));
    assert!(!scratch
        .path()
        .join("regions-corrupt.out/sequence_regions.tsv")
        .exists());
}
