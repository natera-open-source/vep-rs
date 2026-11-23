// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The gathered table: one row a genome, written out by hand from the caches the test
//! reads.

mod common;

use std::fs;
use std::path::Path;

use common::{coding, edit_info, species_cache, transcript, write_shard, Scratch};
use serde_json::json;
use vep_genome_features::extract::{extract, Request};
use vep_genome_features::features::Features;
use vep_genome_features::gather::{gather, COLUMNS};

/// Writes `features` as `<dir>/<name>/features.json`.
fn put(dir: &Path, name: &str, features: &Features) {
    let holder = dir.join(name);
    fs::create_dir_all(&holder).unwrap();
    fs::write(
        holder.join("features.json"),
        serde_json::to_vec_pretty(features).unwrap(),
    )
    .unwrap();
}

/// Two genomes, each a copy of the committed cache with a second coding transcript on
/// the standard code in a second shard of sequence `1`. The first copy's file names no
/// genome and so takes the key of its directory, `b__made_up`. The second copy has
/// SIFT and regulatory data, two synonym pairs, a side file of two rows (`1` of 5,000
/// bases, `MT` of 16,000 bases and circular) and a sequence `MT` of two coding
/// transcripts on code 2, one of them with `cds_end_NF`: one of its four coding
/// transcripts is cut off. It names its own key, `a__made_up`, one directory deeper.
/// The rows stand in key order under the header.
#[test]
fn one_row_a_genome_in_key_order_with_the_counts_as_columns() {
    let scratch = Scratch::new("gather");
    let gathered = scratch.path().join("gathered");
    let second_coding = coding(
        transcript(
            "SECOND_CODING",
            1,
            "protein_coding",
            &[(1_500_001, 1_500_009)],
        ),
        Some(1),
        "ATGAAATAA",
        Some("MK"),
    );
    let plain = Scratch::new("gather");
    let plain_cache = plain.species_copy();
    write_shard(
        &plain_cache,
        "1",
        "1000001-2000000.json",
        std::slice::from_ref(&second_coding),
    );
    put(
        &gathered,
        "b__made_up",
        &extract(&plain_cache, &Request::default()).unwrap(),
    );

    let cache = scratch.species_copy();
    write_shard(&cache, "1", "1000001-2000000.json", &[second_coding]);
    edit_info(&cache, |info| {
        info["sift"] = json!("b");
        info["regulatory"] = json!(true);
    });
    fs::write(cache.join("chr_synonyms.txt"), "1\tchr1\nMT\tchrM\n").unwrap();
    fs::write(
        cache.join("sequence_regions.tsv"),
        "name\tlength\tcircular\n1\t5000\t0\nMT\t16000\t1\n",
    )
    .unwrap();
    let mut cut = coding(
        transcript("MT_REVERSE", -1, "protein_coding", &[(201, 209)]),
        Some(2),
        "ATGAAATAA",
        Some("MK"),
    );
    cut["attributes"] = json!([{"code": "cds_end_NF"}]);
    write_shard(
        &cache,
        "MT",
        "1-1000000.json",
        &[
            coding(
                transcript("MT_FORWARD", 1, "protein_coding", &[(101, 109)]),
                Some(2),
                "ATGAAATAA",
                Some("MK"),
            ),
            cut,
        ],
    );
    let request = Request {
        genome: Some("a__made_up"),
        ..Request::default()
    };
    put(
        &gathered.join("deeper"),
        "another_name",
        &extract(&cache, &request).unwrap(),
    );

    let table = gather(&gathered).unwrap();
    let expected = [
        COLUMNS.join("\t"),
        [
            "a__made_up",
            "synthetic_species",
            "SYN1",
            "5",
            "4",
            "1:2,2:2",
            "2",
            "1",
            "0.250000",
            "2",
            "2",
            "2",
            "MT",
            "16000",
            "1",
            "2",
            "1",
            "1",
            "1",
            "0",
        ]
        .join("\t"),
        [
            "b__made_up",
            "synthetic_species",
            "SYN1",
            "3",
            "2",
            "1:2",
            "2",
            "0",
            "0.000000",
            "1",
            "1",
            "-",
            "-",
            "-",
            "-",
            "0",
            "1",
            "0",
            "0",
            "0",
        ]
        .join("\t"),
    ];
    assert_eq!(table, expected.join("\n") + "\n");
    assert_eq!(COLUMNS[0], "genome");
    assert_eq!(COLUMNS.len(), 20);
}

/// A directory without a features file, a key that stands twice, a file of another
/// shape number, a higher one or a lower one, and a file that is not JSON are each
/// refused with the file's path. A file of another shape number that lacks a key this
/// build reads is refused by its number.
#[test]
fn what_cannot_be_gathered_is_refused_by_name() {
    let scratch = Scratch::new("refused");
    let empty = scratch.path().join("empty");
    fs::create_dir_all(&empty).unwrap();
    let error = gather(&empty).unwrap_err();
    assert_eq!(
        error,
        format!("{}: no features.json under it", empty.display())
    );

    let features = extract(&species_cache(), &Request::default()).unwrap();
    let twice = scratch.path().join("twice");
    put(&twice.join("first"), "same__key", &features);
    put(&twice.join("second"), "same__key", &features);
    let error = gather(&twice).unwrap_err();
    assert_eq!(
        error,
        format!(
            "{}: the genome same__key also stands in {}",
            twice.join("second/same__key/features.json").display(),
            twice.join("first/same__key/features.json").display()
        )
    );

    let other = scratch.path().join("other");
    let mut newer = features.clone();
    newer.features_version = 3;
    put(&other, "newer", &newer);
    let error = gather(&other).unwrap_err();
    assert_eq!(
        error,
        format!(
            "{}: features_version is 3, and this build reads 2",
            other.join("newer/features.json").display()
        )
    );

    let earlier = scratch.path().join("earlier");
    let mut older = features.clone();
    older.features_version = 1;
    put(&earlier, "older", &older);
    let error = gather(&earlier).unwrap_err();
    assert_eq!(
        error,
        format!(
            "{}: features_version is 1, and this build reads 2",
            earlier.join("older/features.json").display()
        )
    );

    let reshaped = scratch.path().join("reshaped");
    let mut without_a_key = serde_json::to_value(&newer).unwrap();
    without_a_key.as_object_mut().unwrap().remove("transcripts");
    fs::create_dir_all(reshaped.join("newer")).unwrap();
    fs::write(
        reshaped.join("newer/features.json"),
        serde_json::to_vec(&without_a_key).unwrap(),
    )
    .unwrap();
    let error = gather(&reshaped).unwrap_err();
    assert_eq!(
        error,
        format!(
            "{}: features_version is 3, and this build reads 2",
            reshaped.join("newer/features.json").display()
        )
    );

    let broken = scratch.path().join("broken");
    fs::create_dir_all(broken.join("genome")).unwrap();
    fs::write(broken.join("genome/features.json"), "{").unwrap();
    let error = gather(&broken).unwrap_err();
    assert!(
        error.starts_with(&broken.join("genome/features.json").display().to_string()),
        "the error is {error}"
    );
}
