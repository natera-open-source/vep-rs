// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The built binary: what it writes, its exit codes and its error lines.

mod common;

use std::fs;

use common::{coding, run, species_cache, transcript, write_shard, Scratch};
use serde_json::{json, Value};
use vep_genome_features::extract::{extract, Request};
use vep_genome_features::features::Features;

/// `extract` writes the features of the committed cache as JSON that reads back to the
/// library's own value, with the genome key the command line names; `gather` over the
/// directory that holds the file writes the header and that genome's row.
#[test]
fn extract_writes_the_features_and_gather_their_table() {
    let scratch = Scratch::new("command");
    let holder = scratch.path().join("made_up__genome");
    fs::create_dir_all(&holder).unwrap();
    let out = holder.join("features.json");
    let cache = species_cache();
    let ended = run(&[
        "extract",
        cache.to_str().unwrap(),
        "--genome",
        "made_up__genome",
        "--threads=2",
        "--out",
        out.to_str().unwrap(),
    ]);
    assert_eq!(ended.code, Some(0), "stderr: {}", ended.stderr);
    assert_eq!(ended.stdout, "");
    assert_eq!(ended.stderr, "");
    let text = fs::read_to_string(&out).unwrap();
    assert!(text.ends_with("}\n"));
    let written: Features = serde_json::from_str(&text).unwrap();
    let request = Request {
        genome: Some("made_up__genome"),
        ..Request::default()
    };
    assert_eq!(written, extract(&cache, &request).unwrap());

    let table = scratch.path().join("features.tsv");
    let ended = run(&[
        "gather",
        scratch.path().to_str().unwrap(),
        "--out",
        table.to_str().unwrap(),
    ]);
    assert_eq!(ended.code, Some(0), "stderr: {}", ended.stderr);
    let lines: Vec<String> = fs::read_to_string(&table)
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[0].starts_with("genome\tspecies\tassembly\t"));
    assert!(lines[1].starts_with("made_up__genome\tsynthetic_species\tSYN1\t2\t1\t"));
}

/// The file `extract` writes for `a_stable_id_with_several_records`, whole.
const SEVERAL_RECORDS_FEATURES: &str = r#"{
  "features_version": 2,
  "genome": "made_up__genome",
  "species": "synthetic_species",
  "assembly": "SYN1",
  "transcripts": 6,
  "coding_transcripts": 4,
  "genetic_codes": [
    {
      "code": 1,
      "coding_transcripts": 4,
      "sequences": [
        {
          "name": "1",
          "coding_transcripts": 1
        },
        {
          "name": "2",
          "coding_transcripts": 3
        }
      ],
      "opening_codons": [
        {
          "codon": "ATG",
          "transcripts": 3,
          "peptide_opens_with_methionine": 3,
          "peptide_opens_otherwise": 0
        },
        {
          "codon": "GTG",
          "transcripts": 1,
          "peptide_opens_with_methionine": 0,
          "peptide_opens_otherwise": 1
        }
      ]
    }
  ],
  "biotypes": [
    {
      "biotype": "lncRNA",
      "transcripts": 2
    },
    {
      "biotype": "protein_coding",
      "transcripts": 4
    }
  ],
  "incomplete_coding": {
    "start": 0,
    "end": 0,
    "transcripts": 0,
    "share": 0.0
  },
  "sequences": {
    "with_features": 2,
    "with_transcripts": 2,
    "side_file": null
  },
  "synonyms": {
    "present": false,
    "pairs": 0
  },
  "variant_data": true,
  "sift_data": false,
  "regulatory_data": false,
  "exons_against_position": 1
}
"#;

/// The one shard of sequence `2` of a made-up cache. The stable id `PLACED_TWICE` stands
/// at two places, a coding transcript on code 1 and a lncRNA. The stable id `SHARED`
/// stands at one place under two `dbID`s: 7, listed first, whose exons run with position
/// and whose peptide opens with valine, and 3, whose exons are listed against position
/// and whose peptide opens with methionine.
fn a_stable_id_with_several_records() -> Vec<Value> {
    let shared = |db_id: u64, exons: &[(u64, u64)], sequence: &str, peptide: &str| {
        let mut record = coding(
            transcript("SHARED", 1, "protein_coding", exons),
            Some(1),
            sequence,
            Some(peptide),
        );
        record["dbID"] = json!(db_id);
        record
    };
    vec![
        coding(
            transcript("PLACED_TWICE", 1, "protein_coding", &[(101, 109)]),
            Some(1),
            "ATGAAATAA",
            Some("MK"),
        ),
        transcript("PLACED_TWICE", 1, "lncRNA", &[(5001, 5100)]),
        shared(7, &[(10_101, 10_150), (10_301, 10_400)], "GTGAAATAA", "VK"),
        shared(3, &[(10_301, 10_400), (10_101, 10_150)], "ATGAAATAA", "MK"),
    ]
}

/// A copy of the committed cache gains sequence `2` with the records of
/// `a_stable_id_with_several_records`. Each record is a transcript of its own, so
/// `extract` writes six transcripts with the committed two, four of them coding on code
/// 1: three open with `ATG` and a peptide that opens with methionine, the committed one
/// among them, and one opens with `GTG` and a peptide that opens with another residue.
/// One transcript lists its exons against position.
#[test]
fn extract_writes_every_record_under_a_stable_id_as_a_transcript_of_its_own() {
    let scratch = Scratch::new("several-records");
    let cache = scratch.species_copy();
    write_shard(
        &cache,
        "2",
        "1-1000000.json",
        &a_stable_id_with_several_records(),
    );
    let out = scratch.path().join("features.json");
    let ended = run(&[
        "extract",
        cache.to_str().unwrap(),
        "--genome",
        "made_up__genome",
        "--out",
        out.to_str().unwrap(),
    ]);
    assert_eq!(ended.code, Some(0), "stderr: {}", ended.stderr);
    assert_eq!(ended.stdout, "");
    assert_eq!(ended.stderr, "");
    assert_eq!(fs::read_to_string(&out).unwrap(), SEVERAL_RECORDS_FEATURES);
}

/// `--version` prints the tool's version and the number of the features file's shape.
#[test]
fn version_names_the_tool_and_the_features_version() {
    let ended = run(&["--version"]);
    assert_eq!(ended.code, Some(0));
    assert_eq!(
        ended.stdout,
        "vep-genome-features 0.1.0 (features version 2)\n"
    );
}

/// A command line the tool cannot read ends with exit 2, one error line and the usage.
#[test]
fn a_command_line_that_cannot_be_read_exits_two() {
    let cases: [(&[&str], &str); 9] = [
        (&[], "no subcommand"),
        (&["list"], "no subcommand list"),
        (
            &["extract", "--out", "f.json"],
            "the directory to read is missing",
        ),
        (&["extract", "cache"], "--out is required"),
        (
            &["extract", "cache", "more", "--out", "f.json"],
            "more is one argument too many",
        ),
        (
            &["extract", "cache", "--out", "f.json", "--seed", "1"],
            "no option --seed",
        ),
        (
            &["extract", "cache", "--out", "f.json", "--threads", "0"],
            "--threads takes a whole number above zero",
        ),
        (
            &["gather", "dir", "--out", "a.tsv", "--out", "b.tsv"],
            "--out is given twice",
        ),
        (&["gather", "dir", "--out"], "--out takes a value"),
    ];
    for (args, what) in cases {
        let ended = run(args);
        assert_eq!(ended.code, Some(2), "{args:?}");
        assert_eq!(ended.stdout, "", "{args:?}");
        let mut lines = ended.stderr.lines();
        assert_eq!(
            lines.next(),
            Some(format!("ERROR: [vep-genome-features] {what}").as_str()),
            "{args:?}"
        );
        assert_eq!(lines.next(), Some("usage:"), "{args:?}");
    }
}

/// A cache that cannot be read ends with exit 1 and one error line, and no output file.
#[test]
fn a_failed_run_exits_one_on_one_error_line() {
    let scratch = Scratch::new("failed");
    let out = scratch.path().join("features.json");
    let ended = run(&[
        "extract",
        scratch.path().to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
    ]);
    assert_eq!(ended.code, Some(1));
    assert_eq!(ended.stdout, "");
    assert_eq!(ended.stderr.lines().count(), 1);
    assert!(ended.stderr.starts_with(&format!(
        "ERROR: [vep-genome-features] {}",
        scratch.path().join("info.json").display()
    )));
    assert!(!out.exists());
}
