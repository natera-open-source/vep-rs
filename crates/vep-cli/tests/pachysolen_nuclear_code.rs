// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The documented difference `pachysolen_nuclear_code_read_as_standard`,
//! through the `vep` binary: a transcript whose cache record names genetic
//! code 26, the nuclear code of Pachysolen tannophilus, under which `CTG`
//! codes alanine. The cache and the records are the reproduction's own,
//! `scripts/divergences/reproductions/pachysolen_nuclear_code_read_as_standard.toml`
//! and the `json_cache` beside its cache builder.
//!
//! Ensembl VEP translates an allele's codon under the transcript's table id
//! (`TranscriptVariationAllele.pm:712` and `:736`) through
//! `Bio::Tools::CodonTable`, which holds no table 26: its `id` warns and
//! falls to id 0 (`CodonTable.pm:328-338`), and `translate` then reads the
//! module's last table, the standard residues (`CodonTable.pm:412-460`). It
//! writes leucine for `CTG`, and `synonymous_variant` where the other codon
//! of the row is a leucine codon.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

use vep_core::cache_manifest::CacheManifest;

mod common;

const REPRODUCTION: &str =
    "scripts/divergences/reproductions/pachysolen_nuclear_code_read_as_standard.toml";
const CACHE: &str =
    "scripts/divergences/reproductions/inputs/pachysolen_nuclear_code_read_as_standard/json_cache";

/// The row of each record of the reproduction, in the order of its input. The
/// coding sequence is `ATG GCT CTG CTA GAA TAA` from base 1,001, and genetic
/// code 26 reads `CTG` as alanine: base 1,008 leaves `CCG` in place of `CTG`
/// (alanine to proline), base 1,009 leaves `CTA` (alanine to leucine) and
/// base 1,012 leaves `CTG` in place of `CTA` (leucine to alanine).
const ROWS: [&str; 3] = [
    "alanine_codon_to_proline\tcontig:1008\tC\tCODEG1\tCODET1\tTranscript\tmissense_variant\t8\t8\t3\tA/P\tcTg/cCg\t-\tIMPACT=MODERATE;STRAND=1",
    "alanine_codon_to_leucine\tcontig:1009\tA\tCODEG1\tCODET1\tTranscript\tmissense_variant\t9\t9\t3\tA/L\tctG/ctA\t-\tIMPACT=MODERATE;STRAND=1",
    "leucine_codon_to_alanine\tcontig:1012\tG\tCODEG1\tCODET1\tTranscript\tmissense_variant\t12\t12\t4\tL/A\tctA/ctG\t-\tIMPACT=MODERATE;STRAND=1",
];

/// The columns the reproduction compares, as indexes of the default output;
/// `IMPACT`, an item of the last column, follows them.
const COMPARED: [usize; 9] = [0, 1, 4, 6, 7, 8, 9, 10, 11];

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
    let shard = Path::new("transcripts/contig/1-1000000.json");
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
        .args([
            "--species",
            "nuclear_code_example",
            "--assembly",
            "EXAMPLE1",
        ])
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
fn a_transcript_on_the_pachysolen_nuclear_code_reads_ctg_as_alanine() {
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
            let impact = columns[13]
                .split(';')
                .find_map(|item| item.strip_prefix("IMPACT="))
                .unwrap();
            let mut cells = COMPARED.map(|index| columns[index]).to_vec();
            cells.push(impact);
            cells.join(" | ")
        })
        .collect();
    compared.sort();
    assert_eq!(recorded, compared);
}
