// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the corpus cache: its directory under either layout, and the transcript
//! records of a made-up cache in the order, and with the skips, of the reader they follow.

use super::*;
use crate::json::Object;
use crate::TempDir;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;
use std::path::Path;

fn manifest(text: &str) -> Object {
    match loads(text) {
        Ok(Value::Object(object)) => object,
        other => panic!("no JSON object: {other:?}"),
    }
}

fn gzip(path: &Path, text: &str) {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(text.as_bytes()).unwrap();
    std::fs::write(path, encoder.finish().unwrap()).unwrap();
}

fn scratch(tag: &str) -> TempDir {
    TempDir::create(&std::env::temp_dir(), &format!("golden_cache_{tag}_")).unwrap()
}

/// The chromosome and the `stable_id` of every record the reader gives.
fn read(cache: &Path) -> Vec<(String, String)> {
    iter_cache_transcripts(cache.to_str().unwrap())
        .map(|item| {
            let (chromosome, transcript) = item.unwrap();
            (chromosome, py_str(transcript.get("stable_id").unwrap()))
        })
        .collect()
}

#[test]
fn the_cache_directory_follows_the_layout_the_species_and_the_flags() {
    let yeast = manifest(
        r#"{"release": "116", "assembly": "R64-1-1", "species": "saccharomyces_cerevisiae"}"#,
    );
    assert_eq!(
        corpus_cache_dir("tmp/yeast", &yeast),
        "tmp/yeast/json_cache"
    );
    let mut laid_out = yeast.clone();
    laid_out.insert("cache_layout", "ensembl");
    assert_eq!(
        corpus_cache_dir("tmp/yeast", &laid_out),
        "tmp/yeast/cache_root/saccharomyces_cerevisiae/116_R64-1-1"
    );
    let human = manifest(r#"{"release": "116", "assembly": "GRCh38", "cache_layout": "ensembl"}"#);
    assert_eq!(
        corpus_cache_dir("tmp/human", &human),
        "tmp/human/cache_root/homo_sapiens/116_GRCh38"
    );
    for (flags, directory) in [
        (r#"["--refseq"]"#, "homo_sapiens_refseq"),
        (r#"["--merged", "--hgvs"]"#, "homo_sapiens_merged"),
        (r#"["--merged", "--refseq"]"#, "homo_sapiens_refseq_merged"),
        (r#"["--refseq=1", 7]"#, "homo_sapiens"),
        ("[]", "homo_sapiens"),
        ("null", "homo_sapiens"),
    ] {
        let mut flagged = human.clone();
        flagged.insert("flags", loads(flags).unwrap());
        assert_eq!(
            corpus_cache_dir("tmp/human/", &flagged),
            format!("tmp/human/cache_root/{directory}/116_GRCh38"),
            "{flags}"
        );
    }
    // A layout other than the Ensembl one is the JSON cache, whatever the flags.
    let mut flagged = manifest(r#"{"flags": ["--refseq"]}"#);
    assert_eq!(corpus_cache_dir(".", &flagged), "json_cache");
    flagged.insert("cache_layout", "other");
    assert_eq!(corpus_cache_dir("c", &flagged), "c/json_cache");
}

#[test]
fn transcripts_come_in_the_order_of_directories_shards_and_records() {
    let scratch = scratch("order");
    let transcripts = scratch.path.join("cache/transcripts");
    for directory in ["1", "10", "2", "X", "empty"] {
        std::fs::create_dir_all(transcripts.join(directory)).unwrap();
    }
    let write = |name: &str, text: &str| std::fs::write(transcripts.join(name), text).unwrap();
    write("notes.txt", "no directory\n");
    // A shard of the converted layout: the records of each chromosome under its name, the
    // empty places and the values that hold no list left out.
    gzip(
        &transcripts.join("1/1-1000000.json.gz"),
        concat!(
            r#"{"1": [{"stable_id": "T1", "start": 5}, null, {"stable_id": "", "start": 6}, "#,
            r#"{"start": 7}, {"stable_id": "T2"}], "other": null, "text": "ab", "#,
            r#""map": {"k": 1}, "zero": 0}"#
        ),
    );
    // A shard of the native layout: a list, of which only an object with a `stable_id` that
    // is not empty is a record.
    write(
        "1/1000001-2000000.json",
        r#"[{"stable_id": "T3"}, null, 5, "x", [1], {"stable_id": 0}, {"stable_id": 7}]"#,
    );
    write("1/bad.json", "{not json");
    write("1/not-gzip.json.gz", "[]");
    write("1/skipped.txt", r#"[{"stable_id": "NO"}]"#);
    write("1/scalar.json", "5");
    std::fs::write(
        transcripts.join("1/latin.json"),
        b"[{\"stable_id\": \"\xe9\"}]",
    )
    .unwrap();
    std::fs::write(
        transcripts.join("1/bom.json"),
        b"\xef\xbb\xbf[{\"stable_id\": \"BOM\"}]",
    )
    .unwrap();
    write("1/.json", r#"[{"stable_id": "HIDDEN"}]"#);
    std::fs::create_dir(transcripts.join("1/dir.json")).unwrap();
    gzip(
        &transcripts.join("10/a.json.gz"),
        r#"[{"stable_id": "T10"}]"#,
    );
    write(
        "2/B.json",
        r#"{"2": [{"stable_id": "TB"}], "3": [{"stable_id": "TB3"}]}"#,
    );
    write("2/a.json", r#"[{"stable_id": "Ta"}]"#);
    // A name written twice holds its last value.
    write(
        "X/x.json",
        r#"{"X": [{"stable_id": "TX", "n": 1.5}], "X": [{"stable_id": "TX2"}]}"#,
    );
    let expected = [
        ("1", "HIDDEN"),
        ("1", "T1"),
        ("1", "T2"),
        ("1", "T3"),
        ("1", "7"),
        ("10", "T10"),
        ("2", "TB"),
        ("2", "TB3"),
        ("2", "Ta"),
        ("X", "TX2"),
    ];
    let expected: Vec<(String, String)> = expected
        .iter()
        .map(|(chromosome, id)| (chromosome.to_string(), id.to_string()))
        .collect();
    assert_eq!(read(&scratch.path.join("cache")), expected);
    // A record is given whole.
    let second = iter_cache_transcripts(scratch.path.join("cache").to_str().unwrap())
        .nth(1)
        .unwrap()
        .unwrap();
    assert_eq!(second.0, "1");
    assert_eq!(second.1, manifest(r#"{"stable_id": "T1", "start": 5}"#));
}

#[test]
fn a_cache_without_transcripts_holds_no_record() {
    let scratch = scratch("absent");
    assert_eq!(read(&scratch.path.join("none")), []);
    std::fs::create_dir(scratch.path.join("cache")).unwrap();
    assert_eq!(read(&scratch.path.join("cache")), []);
    std::fs::write(scratch.path.join("cache/transcripts"), "a file").unwrap();
    assert_eq!(read(&scratch.path.join("cache")), []);
}

#[test]
fn a_shard_the_reader_cannot_hold_ends_the_reading_in_the_tools_own_words() {
    let scratch = scratch("refused");
    let cache = scratch.path.join("cache");
    let directory = cache.join("transcripts/1");
    std::fs::create_dir_all(&directory).unwrap();
    std::fs::write(directory.join("a.json"), r#"[{"stable_id": "T1"}]"#).unwrap();
    std::fs::write(directory.join("b.json"), r#"{"1": 5}"#).unwrap();
    std::fs::write(directory.join("c.json"), r#"[{"stable_id": "T3"}]"#).unwrap();
    let mut records = iter_cache_transcripts(cache.to_str().unwrap());
    assert_eq!(records.next().unwrap().unwrap().0, "1");
    let refusal = records.next().unwrap().unwrap_err();
    let shard = directory.join("b.json");
    assert_eq!(
        refusal.line,
        format!(
            "ERROR: [vep-golden-corpus] {} holds a number where the records of 1 stand",
            shard.to_str().unwrap()
        )
    );
    assert_eq!(refusal.code, 1);
    assert!(records.next().is_none());
    // A lone surrogate escape is refused; an integer past the reader's bound skips the
    // shard.
    let lone = format!("[{{\"stable_id\": \"{}ud83d\"}}]", '\\');
    std::fs::write(directory.join("b.json"), lone).unwrap();
    let outcomes: Vec<bool> = iter_cache_transcripts(cache.to_str().unwrap())
        .map(|item| item.is_ok())
        .collect();
    assert_eq!(outcomes, [true, false]);
    let long = format!("[{{\"stable_id\": {}}}]", "9".repeat(4301));
    std::fs::write(directory.join("b.json"), long).unwrap();
    assert_eq!(
        read(&cache),
        [
            ("1".to_string(), "T1".to_string()),
            ("1".to_string(), "T3".to_string())
        ]
    );
}
