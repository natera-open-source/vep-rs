// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The documented difference `transcript_across_circular_origin`, through the
//! `vep` binary: a transcript an Ensembl cache stores across the origin of a
//! circular sequence, with a start below 1 and a mapper pair that keeps the
//! database coordinates. The cache and the records are the reproduction's own,
//! `scripts/divergences/reproductions/transcript_across_circular_origin.toml`
//! and the `json_cache` beside its cache builder.
//!
//! Ensembl VEP annotates such a transcript on neither side of the origin. A
//! transcript takes the records between its stored start less the flank and
//! its stored end plus the flank (`AnnotationType/Transcript.pm:103-109`), so
//! a record among the sequence's last bases gets no row. A record among the
//! sequence's first bases maps to no transcript position, the interval query
//! over the mapper pairs finding none (ensembl core `Mapper.pm:311-316`,
//! `:431-435`), no predicate holds, and the row carries the default term
//! (`BaseVariationFeatureOverlapAllele.pm:285`).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

use vep_core::cache_manifest::CacheManifest;

mod common;

const REPRODUCTION: &str =
    "scripts/divergences/reproductions/transcript_across_circular_origin.toml";
const CACHE: &str =
    "scripts/divergences/reproductions/inputs/transcript_across_circular_origin/json_cache";

/// The row of each record of the reproduction, in the order of its input: a
/// missense change at coding base 121, the kept stop codon at the gene's last
/// base, the upstream row 911 bases before the gene's first base and a
/// missense change at coding base 40. Each follows from the record's position
/// on a sequence of 12,000 bases whose gene runs from base 11,911 through the
/// origin to base 60.
const ROWS: [&str; 4] = [
    "low_side_31\tring:31\tG\tRINGG1\tRINGT1\tTranscript\tmissense_variant\t121\t121\t41\tT/A\tAct/Gct\t-\tIMPACT=MODERATE;STRAND=1",
    "last_gene_base_60\tring:60\tG\tRINGG1\tRINGT1\tTranscript\tstop_retained_variant\t150\t150\t50\t*\ttaA/taG\t-\tIMPACT=LOW;STRAND=1",
    "before_the_gene_11000\tring:11000\tG\tRINGG1\tRINGT1\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=911;STRAND=1",
    "high_side_11950\tring:11950\tG\tRINGG1\tRINGT1\tTranscript\tmissense_variant\t40\t40\t14\tI/V\tAtt/Gtt\t-\tIMPACT=MODERATE;STRAND=1",
];

/// The columns the reproduction compares, as indexes of the default output.
const COMPARED: [usize; 9] = [0, 1, 4, 5, 6, 7, 8, 9, 10];

fn repository() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The lines of the TOML basic multi-line string or array that opens on the
/// line starting with `opening`, up to the line that is `closing`.
fn block<'a>(text: &'a str, opening: &str, closing: &str) -> Vec<&'a str> {
    text.lines()
        .skip_while(|line| !line.starts_with(opening))
        .skip(1)
        .take_while(|line| *line != closing)
        .collect()
}

/// The reproduction's input records and the rows it records for vep-rs.
fn reproduction() -> (Vec<String>, Vec<String>) {
    let text = fs::read_to_string(repository().join(REPRODUCTION)).unwrap();
    let records = block(&text, "text = \"\"\"", "\"\"\"")
        .into_iter()
        .map(str::to_string)
        .collect();
    let recorded = block(&text, "vep_rs = [", "]")
        .into_iter()
        .map(|line| {
            line.trim()
                .trim_end_matches(',')
                .trim_matches('"')
                .to_string()
        })
        .collect();
    (records, recorded)
}

/// Runs `vep` on `records` over a copy of the committed cache that carries its
/// checksum manifest, and returns the rows it writes.
fn annotate(records: &[String]) -> Vec<String> {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("cache");
    let shard = Path::new("transcripts/ring/1-1000000.json");
    fs::create_dir_all(cache.join(shard).parent().unwrap()).unwrap();
    for file in [Path::new("info.json"), shard] {
        fs::copy(repository().join(CACHE).join(file), cache.join(file)).unwrap();
    }
    CacheManifest::write(&cache).unwrap();
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for record in records {
        text.push_str(record);
        text.push('\n');
    }
    let input = dir.path().join("records.vcf");
    fs::write(&input, text).unwrap();
    let out = dir.path().join("out.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(&cache)
        .args(["--species", "circular_example", "--assembly", "RING1"])
        .args(["--format", "vcf"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let warnings = dir.path().join("out.txt_warnings.txt");
    assert!(
        !warnings.exists(),
        "{}",
        fs::read_to_string(&warnings).unwrap_or_default()
    );
    fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

#[test]
fn a_transcript_stored_across_the_origin_is_annotated_on_both_sides_of_it() {
    let (records, _) = reproduction();
    assert_eq!(records.len(), ROWS.len());
    assert_eq!(annotate(&records), ROWS);
}

#[test]
fn the_reproduction_records_the_rows_vep_rs_writes() {
    let (_, recorded) = reproduction();
    let mut compared: Vec<String> = ROWS
        .iter()
        .map(|row| {
            let columns: Vec<&str> = row.split('\t').collect();
            COMPARED.map(|index| columns[index]).join(" | ")
        })
        .collect();
    compared.sort();
    assert_eq!(recorded, compared);
}
