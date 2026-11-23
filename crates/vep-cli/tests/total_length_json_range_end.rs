// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--total_length` in the JSON format: the end key of a position range holds the range's
//! last position with its total (`"487/2674"`), where Ensembl VEP writes the first position
//! under both keys and neither key of a range whose first position is unknown
//! (`OutputFactory/JSON.pm:337` at release/116.3). The two records are those of the
//! reproduction of the documented difference `total_length_json_range_end`; the first is a
//! record of the committed corpus `GRCh37-output-flags`, whose stored Ensembl VEP outputs the
//! expected objects are read from.

mod common;

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use serde_json::{json, Value};

const TRANSCRIPT: &str = "ENST00000298232";
const CORPUS_RECORD: &str = "rs11184139";
const MADE_UP_RECORD: &str = "del_over_exon_start";

fn repository() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn corpus() -> PathBuf {
    repository().join("tests/golden/116/GRCh37-output-flags")
}

/// The lines of the reproduction's inline input, the corpus's deletion first.
fn reproduction_records() -> Vec<String> {
    let path =
        repository().join("scripts/divergences/reproductions/total_length_json_range_end.toml");
    let text = fs::read_to_string(path).unwrap();
    let (_, after) = text.split_once("text = \"\"\"\n").expect("an inline input");
    let (records, _) = after.split_once("\"\"\"").expect("the input's end");
    records.lines().map(str::to_string).collect()
}

/// A stored Ensembl VEP output of the corpus, decompressed.
fn stored(expected_dir: &str, name: &str) -> String {
    let bytes = fs::read(corpus().join(expected_dir).join(name)).unwrap();
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(&bytes[..])
        .read_to_string(&mut text)
        .unwrap();
    text
}

/// The object of [`TRANSCRIPT`] in the JSON record whose `input` opens with the first five
/// columns of `record`.
fn transcript_object(output: &str, record: &str) -> Value {
    let key: Vec<&str> = record.split('\t').take(5).collect();
    let key = key.join("\t");
    let object: Value = output
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|object| {
            object["input"]
                .as_str()
                .is_some_and(|input| input.starts_with(&key))
        })
        .expect("the record's object");
    object["transcript_consequences"]
        .as_array()
        .expect("transcript rows")
        .iter()
        .find(|row| row["transcript_id"] == TRANSCRIPT)
        .expect("the transcript's row")
        .clone()
}

/// The JSON output of `vep` on `records` over the corpus's cache with the flags of the
/// corpus's own runs and `flags`.
fn annotate(name: &str, records: &[String], flags: &[&str]) -> String {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.vcf");
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    fs::write(&input, format!("{header}{}\n", records.join("\n"))).unwrap();
    let out = dir.path().join(format!("{name}.json"));
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["--offline", "--json_cache"])
        .arg(corpus().join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh37"])
        .args([
            "--buffer_size",
            "5000",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .arg("--fasta")
        .arg(repository().join("tests/golden/116/GRCh37-hgvs/reference.fa.gz"))
        .arg("--json")
        .args(flags)
        .output()
        .expect("vep binary runs");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read_to_string(out).unwrap()
}

/// The corpus's deletion of four bases inside one exon, cDNA 484 to 487 of 2,674, coding
/// bases 116 to 119 of 1,602 and residues 39 to 40 of 533. The whole object equals the one
/// Ensembl VEP 116.3 stored for the flag set but for the three end keys, which hold the last
/// position with its total where Ensembl VEP repeats the start.
#[test]
fn a_range_end_under_total_length_holds_the_last_position_and_the_total() {
    let records = reproduction_records();
    let record = &records[0];
    assert!(
        fs::read_to_string(corpus().join("variants.vcf"))
            .unwrap()
            .lines()
            .any(|line| line == record),
        "the reproduction's first record is a record of the corpus"
    );
    assert!(record.contains(CORPUS_RECORD));
    let reference = transcript_object(&stored("expected_total-length", "json.jsonl.gz"), record);
    for (start, end, first) in [
        ("cdna_start", "cdna_end", 484),
        ("cds_start", "cds_end", 116),
        ("protein_start", "protein_end", 39),
    ] {
        assert_eq!(reference[start], json!(first));
        assert_eq!(reference[end], json!(first), "Ensembl VEP's {end}");
    }
    let mut expected = reference;
    expected["cdna_end"] = json!("487/2674");
    expected["cds_end"] = json!("119/1602");
    expected["protein_end"] = json!("40/533");
    let written = annotate("range", &records[..1], &["--total_length"]);
    assert_eq!(transcript_object(&written, record), expected);
}

/// Without the flag the object is Ensembl VEP's stored one, whose end keys are the numbers
/// 487, 119 and 40: the last positions the flag's end keys carry.
#[test]
fn without_the_flag_the_end_keys_are_the_last_positions_as_numbers() {
    let records = reproduction_records();
    let record = &records[0];
    let reference = transcript_object(&stored("expected", "json.jsonl.gz"), record);
    assert_eq!(reference["cdna_end"], json!(487));
    assert_eq!(reference["cds_end"], json!(119));
    assert_eq!(reference["protein_end"], json!(40));
    let written = annotate("plain", &records[..1], &[]);
    assert_eq!(transcript_object(&written, record), reference);
}

/// A made-up deletion of two intron bases and the first two bases of the same exon, whose
/// first base is cDNA 434, coding base 66: the first position is unknown and the last is
/// cDNA 435, coding base 67, residue 23. The end keys stand alone, as `?-435` writes
/// `cdna_end` alone without the flag.
#[test]
fn a_range_of_unknown_start_writes_its_end_key_alone_under_total_length() {
    let records = reproduction_records();
    let record = &records[1];
    assert!(record.contains(MADE_UP_RECORD));
    let row = transcript_object(
        &annotate("open", &records[1..], &["--total_length"]),
        record,
    );
    let ends: Vec<(&str, &Value)> = row
        .as_object()
        .unwrap()
        .iter()
        .filter(|(key, _)| key.ends_with("_start") || key.ends_with("_end"))
        .map(|(key, value)| (key.as_str(), value))
        .collect();
    assert_eq!(
        ends,
        [
            ("cdna_end", &json!("435/2674")),
            ("cds_end", &json!("67/1602")),
            ("protein_end", &json!("23/533")),
        ]
    );
    let plain = transcript_object(&annotate("open-plain", &records[1..], &[]), record);
    assert_eq!(
        (&plain["cdna_end"], &plain["cds_end"], &plain["protein_end"]),
        (&json!(435), &json!(67), &json!(23))
    );
    assert_eq!(plain.get("cdna_start"), None);
}

/// The reproduction records the values this file pins for vep-rs.
#[test]
fn the_reproduction_records_the_end_keys_vep_rs_writes() {
    let path =
        repository().join("scripts/divergences/reproductions/total_length_json_range_end.toml");
    let text = fs::read_to_string(path).unwrap();
    assert!(text.contains("ensembl_vep = [\"484\"]\n"));
    assert!(text.contains("vep_rs = [\"435/2674\", \"487/2674\"]\n"));
}
