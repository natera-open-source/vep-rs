// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--plugin NearestExonJB` on a transcript whose exons have no stable id,
//! through the `vep` binary: the item is written with the empty id in every
//! format.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

mod common;

/// The GRCh38-merged corpus, whose cache carries RefSeq aligned transcripts
/// (`compmerge.*`) with no exon stable ids beside the Ensembl and RefSeq
/// models.
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-merged")
}

/// A substitution 23 bases into the first exon of `compmerge.69.pooled.chr21`
/// (10521527-10521694), which `--all_refseq` admits to the output.
const RECORD: &str = "21\t10521550\tidless_exon1\tG\tA";

/// The 18 aligned transcripts the record lies in an exon of, with the item
/// Ensembl VEP 116.3 writes on each under `--merged --all_refseq`: the exon's
/// id is empty, the distance is to the exon's start and the length is the
/// exon's.
const EXPECTED: &[(&str, &str)] = &[
    ("compmerge.61.Testes.chr21", "|23|start|168"),
    ("compmerge.62.Testes.chr21", "|23|start|168"),
    ("compmerge.63.Testes.chr21", "|23|start|168"),
    ("compmerge.68.pooled.chr21", "|23|start|168"),
    ("compmerge.69.pooled.chr21", "|23|start|168"),
    ("compmerge.70.pooled.chr21", "|23|start|168"),
    ("compmerge.56.Testes.chr21", "|6|start|151"),
    ("compmerge.57.Testes.chr21", "|6|start|151"),
    ("compmerge.58.Testes.chr21", "|6|start|151"),
    ("compmerge.59.Testes.chr21", "|6|start|151"),
    ("compmerge.63.pooled.chr21", "|6|start|151"),
    ("compmerge.64.pooled.chr21", "|6|start|151"),
    ("compmerge.66.pooled.chr21", "|6|start|151"),
    ("compmerge.67.pooled.chr21", "|6|start|151"),
    ("compmerge.183.Testes.chr21", "|4|start|149"),
    ("compmerge.383.pooled.chr21", "|4|start|149"),
    ("compmerge.185.Testes.chr21", "|15|start|160"),
    ("compmerge.385.pooled.chr21", "|15|start|160"),
];

/// Runs the binary on the record in `format` and returns the output's lines,
/// the `##` header lines included.
fn annotate(dir: &Path, format: &str) -> Vec<String> {
    let input = dir.join("in.vcf");
    std::fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{RECORD}\t.\t.\t.\n"),
    )
    .unwrap();
    let output = dir.join(format!("out.{format}"));
    let mut command = common::isolated_vep_anywhere();
    command
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .args([
            "--json_cache",
            corpus().join("json_cache").to_str().unwrap(),
        ])
        .args([
            "--fasta",
            corpus().join("reference.fa.gz").to_str().unwrap(),
        ])
        .args([
            "--offline",
            "--assembly",
            "GRCh38",
            "--merged",
            "--all_refseq",
        ])
        .args(["--plugin", "NearestExonJB"])
        .args(["--force_overwrite", "--no_stats", "--quiet"]);
    if format != "default" {
        command.arg(format!("--{format}"));
    }
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "--{format}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

/// The `NearestExonJB` value of every aligned-transcript row that carries one,
/// by transcript, in the default or tab format.
fn text_items(rows: &[String], format: &str) -> BTreeMap<String, String> {
    let column = (format == "tab").then(|| {
        let header = rows
            .iter()
            .find(|r| r.starts_with('#') && !r.starts_with("##"))
            .unwrap();
        header
            .split('\t')
            .position(|c| c == "NearestExonJB")
            .unwrap()
    });
    rows.iter()
        .filter(|r| !r.starts_with('#'))
        .filter_map(|row| {
            let columns: Vec<&str> = row.split('\t').collect();
            let feature = columns[4];
            if !feature.starts_with("compmerge.") {
                return None;
            }
            let item = match column {
                Some(index) => Some(columns[index]).filter(|value| *value != "-"),
                None => columns[13]
                    .split(';')
                    .find_map(|kv| kv.strip_prefix("NearestExonJB=")),
            };
            item.map(|item| (feature.to_string(), item.to_string()))
        })
        .collect()
}

/// The `NearestExonJB` value of every aligned-transcript CSQ entry that carries
/// one, by transcript, the `+`-joined parts as the entry prints them.
fn vcf_items(rows: &[String]) -> BTreeMap<String, String> {
    let csq_format = rows
        .iter()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .map(|fields| fields.trim_end_matches("\">"))
        .unwrap_or_default();
    let field = |name: &str| {
        csq_format
            .split('|')
            .position(|f| f == name)
            .unwrap_or_else(|| panic!("CSQ format lacks {name}: {csq_format}"))
    };
    let (feature_field, item_field) = (field("Feature"), field("NearestExonJB"));
    let line = rows
        .iter()
        .find(|r| !r.starts_with('#'))
        .expect("--vcf: the record's line");
    let csq = line
        .split('\t')
        .nth(7)
        .unwrap()
        .split(';')
        .find_map(|kv| kv.strip_prefix("CSQ="))
        .expect("--vcf: a CSQ field");
    csq.split(',')
        .filter_map(|entry| {
            let fields: Vec<&str> = entry.split('|').collect();
            let feature = fields[feature_field];
            (feature.starts_with("compmerge.") && !fields[item_field].is_empty())
                .then(|| (feature.to_string(), fields[item_field].to_string()))
        })
        .collect()
}

/// The `nearestexonjb` value of every aligned-transcript consequence that
/// carries one, by transcript.
fn json_items(rows: &[String]) -> BTreeMap<String, String> {
    let record: serde_json::Value = serde_json::from_str(&rows[0]).unwrap();
    record["transcript_consequences"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|tc| {
            let feature = tc["transcript_id"].as_str().unwrap();
            let item = tc.get("nearestexonjb")?.as_str()?;
            feature
                .starts_with("compmerge.")
                .then(|| (feature.to_string(), item.to_string()))
        })
        .collect()
}

/// Ensembl VEP 116.3 writes the item of an exon that has no stable id with the
/// empty id (`NearestExonJB.pm:220, :269`): `|23|start|168` on each of the 18
/// aligned transcripts of `21 10521550 G A` whose exon holds it, `+23+start+168`
/// under `--vcf`, and the row carries the field in every format; the aligned
/// transcripts the record lies upstream of carry none.
#[test]
fn an_exon_without_a_stable_id_is_written_with_the_empty_id_in_every_format() {
    let dir = tempfile::tempdir().unwrap();
    let expected: BTreeMap<String, String> = EXPECTED
        .iter()
        .map(|(feature, item)| (feature.to_string(), item.to_string()))
        .collect();
    for format in ["default", "tab"] {
        let rows = annotate(dir.path(), format);
        assert_eq!(text_items(&rows, format), expected, "--{format}");
    }
    let rows = annotate(dir.path(), "vcf");
    let plus_joined: BTreeMap<String, String> = expected
        .iter()
        .map(|(feature, item)| (feature.clone(), item.replace('|', "+")))
        .collect();
    assert_eq!(vcf_items(&rows), plus_joined, "--vcf");
    let rows = annotate(dir.path(), "json");
    assert_eq!(json_items(&rows), expected, "--json");
}
