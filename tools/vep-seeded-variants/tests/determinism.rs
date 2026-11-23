// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The determinism contract: the same cache and seed give the same bytes, whatever
//! order a directory lists its files in, and the light records of each pilot slice for
//! seed 1 are the committed ones.

mod common;

use std::fs;
use std::path::Path;

use common::{copy_tree, pilot_cache, run, run_to_bytes, Scratch, PILOTS};
use vep_seeded_variants::cache::Cache;
use vep_seeded_variants::GENERATOR_VERSION;

#[test]
fn the_same_cache_and_seed_give_the_same_bytes() {
    let cache = pilot_cache("BDGP6.54-fly");
    let scratch = Scratch::new("same-bytes");
    let first = run_to_bytes(&cache, &scratch, &["--seed", "1"]);
    let second = run_to_bytes(&cache, &scratch, &["--seed", "1"]);
    assert!(first == second, "two runs of one cache and seed differ");

    let other_seed = run_to_bytes(&cache, &scratch, &["--seed", "2"]);
    assert!(
        first.0 != other_seed.0,
        "seeds 1 and 2 give the same records"
    );

    let rebuilt = scratch.path().join("rebuilt");
    copy_tree(&cache, &rebuilt, true);
    let reversed = run_to_bytes(&rebuilt, &scratch, &["--seed", "1"]);
    assert!(
        first == reversed,
        "a copy of the cache built from the last file to the first gives other bytes"
    );

    let wide = run_to_bytes(&cache, &scratch, &["--seed", "1", "--records", "wide"]);
    let wide_again = run_to_bytes(&rebuilt, &scratch, &["--seed", "1", "--records", "wide"]);
    assert!(wide == wide_again, "two wide runs differ");
    assert!(wide.0 != first.0, "the wide records are the light ones");
}

#[test]
fn the_shards_of_a_sequence_are_read_by_their_start() {
    let cache = Cache::open(&pilot_cache("BDGP6.54-fly")).expect("a cache");
    let shards = cache.shards("2L").expect("the shards of 2L");
    let starts: Vec<u64> = shards
        .iter()
        .map(|shard| {
            let name = shard.file_name().expect("a name").to_string_lossy();
            name.split('-')
                .next()
                .and_then(|start| start.parse().ok())
                .expect("a start")
        })
        .collect();
    assert!(starts.len() > 2, "2L has {} shards", starts.len());
    assert!(
        starts.windows(2).all(|pair| pair[0] < pair[1]),
        "{starts:?}"
    );
}

#[test]
fn the_output_names_the_generator_the_seed_and_the_mode() {
    let cache = pilot_cache("R64-1-1-yeast");
    let scratch = Scratch::new("header");
    let (records, summary) = run_to_bytes(&cache, &scratch, &["--seed", "41"]);
    let text = String::from_utf8(records).expect("text");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(
        &lines[..6],
        [
            "##fileformat=VCFv4.2",
            &format!("##source=vep-seeded-variants generator_version={GENERATOR_VERSION}"),
            "##seed=41",
            "##records=light",
            "##reference_alleles=cache",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        ]
    );
    let order: Vec<String> = Cache::open(&cache).expect("a cache").sequences().to_vec();
    let mut last = (0, 0u64);
    let mut ids = Vec::new();
    for line in &lines[6..] {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 8, "{line}");
        assert_eq!(&fields[5..], [".", ".", "."]);
        ids.push(fields[2].to_owned());
        // A synonym record names its sequence another way and sorts with the sequence.
        if let Some(index) = order.iter().position(|name| name == fields[0]) {
            let place = (index, fields[1].parse::<u64>().expect("a position"));
            assert!(last <= place, "{line} stands before its turn");
            last = place;
        }
    }
    let count = ids.len();
    ids.sort();
    ids.dedup();
    assert_eq!(ids.len(), count, "an ID stands twice");
    let summary = String::from_utf8(summary).expect("text");
    assert!(summary.starts_with("category\twanted\twritten\tpassed_over\treason\n"));
    let version = run(&["--version"]);
    assert_eq!(version.code, Some(0));
    assert_eq!(
        version.stdout,
        format!("vep-seeded-variants 0.1.0 (generator version {GENERATOR_VERSION})\n")
    );
}

#[test]
fn a_run_that_cannot_start_ends_on_one_error_line() {
    let scratch = Scratch::new("errors");
    let output = scratch.path().join("out.vcf");
    let output = output.to_str().expect("a path");
    let missing = run(&[
        "--cache",
        "/no/such/cache",
        "--seed",
        "1",
        "--output",
        output,
    ]);
    assert_eq!(missing.code, Some(1));
    assert!(missing
        .stderr
        .starts_with("ERROR: [vep-seeded-variants] /no/such/cache/info.json: "));
    assert_eq!(missing.stderr.lines().count(), 1);

    let cache = pilot_cache("R64-1-1-yeast");
    let cache = cache.to_str().expect("a path");
    for bad in [
        vec!["--cache", cache, "--output", output],
        vec!["--cache", cache, "--seed", "x", "--output", output],
        vec![
            "--cache", cache, "--seed", "1", "--seed", "2", "--output", output,
        ],
        vec![
            "--cache",
            cache,
            "--seed",
            "1",
            "--output",
            output,
            "--records",
            "deep",
        ],
        vec![
            "--cache",
            cache,
            "--seed",
            "1",
            "--output",
            output,
            "--per-category",
            "3",
        ],
        vec![
            "--cache", cache, "--seed", "1", "--output", output, "--depth", "3",
        ],
        vec!["--cache", cache, "--seed", "1", "--output", output, "stray"],
        vec![
            "--cache",
            cache,
            "--seed",
            "1",
            "--output",
            output,
            "--summary",
            output,
        ],
    ] {
        let ended = run(&bad);
        assert_eq!(ended.code, Some(2), "{bad:?}");
        assert!(
            ended.stderr.starts_with("ERROR: [vep-seeded-variants] "),
            "{bad:?}: {}",
            ended.stderr
        );
    }
    let unknown = run(&[
        "--cache",
        cache,
        "--seed",
        "1",
        "--output",
        output,
        "--sequence",
        "nowhere",
    ]);
    assert_eq!(unknown.code, Some(1));
    assert!(unknown
        .stderr
        .contains("the cache holds no sequence nowhere"));
}

#[test]
fn the_light_records_of_each_pilot_are_pinned() {
    let scratch = Scratch::new("pinned");
    let update = std::env::var_os("VEP_SEEDED_VARIANTS_UPDATE").is_some();
    for (corpus, _) in PILOTS {
        let (records, _) = run_to_bytes(&pilot_cache(corpus), &scratch, &["--seed", "1"]);
        let pinned = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/expected")
            .join(format!("{corpus}.vcf"));
        if update {
            fs::write(&pinned, &records).expect("a written file");
            continue;
        }
        let expected = fs::read(&pinned).unwrap_or_default();
        assert!(
            records == expected,
            "the light records of {corpus} for seed 1 are not those of {}. A deliberate \
             change of the algorithm bumps GENERATOR_VERSION in src/lib.rs and regenerates \
             the files under tests/expected: run the tests once with \
             VEP_SEEDED_VARIANTS_UPDATE=1.",
            pinned.display()
        );
    }
}
