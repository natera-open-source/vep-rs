// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` key that spells a column or a flag field that is a list, on
//! the rows the binary writes for the records of the documented difference
//! `custom_key_list_shared_across_rows`: lines 9 and 10 of the custom corpus's
//! input with a twin of line 9 whose alternate alleles stand in the other
//! order, and line 42, one allele with sixteen rows. The corpus's
//! `custom/sites.vcf.gz` holds one record at each site and is read under the
//! key's name, over the cache of the `GRCh38-everything` corpus.
//!
//! Each row lists the source's values and the row's own items, each once.
//! Ensembl VEP 116.3 gives every row of a record the one list its record's
//! hash holds under the key (OutputFactory.pm:391): the values of an `overlap`
//! source's record (:1003-1014) or the nearest names (:1017). Each row pushes
//! onto it: the co-located identifiers (:1086), the synonyms (:1090), which
//! the row then joins into a string of its own (:1141), `MANE_Select` (:1601)
//! and the records of the row's allele (:1426-1435). The lines are written
//! once every row is built (:266-272). So Ensembl VEP lists an identifier once
//! for each row of its record, repeats the synonyms once more on each later
//! row, writes `MANE_Select` on the rows of every transcript and lists an
//! allele's record once for each row. The reproduction record of the entry
//! holds Ensembl VEP's cells beside the cells these tests hold.

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The start of the two corpus records of the first run, on lines 9 and 10 of
/// the custom corpus's `variants.vcf`.
const TWO_ALLELE_LINES: [(usize, &str); 2] = [(9, "21\t5380280\t"), (10, "21\t8987421\t")];

/// The twin of line 9: its alternate alleles in the other order.
const TWIN: &str = "21\t5380280\talleles_reversed\tA\tG,C\t.\t.\t.";

/// The start of the record of sixteen rows, on line 42 of the same file.
const SIXTEEN_ROW_LINE: (usize, &str) = (42, "21\t33549625\t");

/// The sixteen transcripts of that record's rows, in the order of the rows.
const SIXTEEN_TRANSCRIPTS: [&str; 16] = [
    "ENST00000290239",
    "ENST00000300278",
    "ENST00000356577",
    "ENST00000381679",
    "ENST00000381692",
    "ENST00000436227",
    "ENST00000455528",
    "ENST00000475072",
    "ENST00000492229",
    "ENST00000695559",
    "ENST00000704334",
    "ENST00001036731",
    "ENST00001036732",
    "ENST00001129604",
    "ENST00001129605",
    "ENST00001131645",
];

/// The transcript of the MANE Select set among them.
const MANE_SELECT_TRANSCRIPT: &str = "ENST00000356577";

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The record on the 1-based `line` of the custom corpus's input, which opens
/// with `start`.
fn corpus_record(input: &str, (line, start): (usize, &str)) -> String {
    let record = input.lines().nth(line - 1).expect("the record's line");
    assert!(
        record.starts_with(start),
        "line {line} of the custom corpus's input is the record"
    );
    record.to_string()
}

/// A VCF file under `dir` that holds the header lines of the custom corpus's
/// input and `records`.
fn records_file(dir: &Path, records: &[String]) -> PathBuf {
    let input = fs::read_to_string(corpus("GRCh38-custom").join("variants.vcf")).unwrap();
    let mut text: String = input
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    for record in records {
        text.push_str(record);
        text.push('\n');
    }
    let path = dir.join("records.vcf");
    fs::write(&path, text).unwrap();
    path
}

/// The three records of two alternate alleles: line 9, its twin and line 10.
fn two_allele_records() -> Vec<String> {
    let input = fs::read_to_string(corpus("GRCh38-custom").join("variants.vcf")).unwrap();
    let [first, second] = TWO_ALLELE_LINES.map(|line| corpus_record(&input, line));
    vec![first, TWIN.to_string(), second]
}

/// The one record of sixteen rows: line 42.
fn sixteen_row_record() -> Vec<String> {
    let input = fs::read_to_string(corpus("GRCh38-custom").join("variants.vcf")).unwrap();
    vec![corpus_record(&input, SIXTEEN_ROW_LINE)]
}

/// The cells `(Uploaded_variation, Allele, Feature, key)` of every row of a
/// tab run over `records` under `flags`, with the custom corpus's site file
/// as a source of `kind` (`overlap` or `exact`) named `key` that asks for its
/// `AF` and `CLASS` fields.
fn key_cells(
    records: &[String],
    key: &str,
    kind: &str,
    flags: &[&str],
) -> Vec<(String, String, String, String)> {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let out = dir.join("tab.out");
    let source = format!(
        "file={},short_name={key},format=vcf,type={kind},fields=AF%CLASS",
        corpus("GRCh38-custom")
            .join("custom/sites.vcf.gz")
            .display()
    );
    let run = common::isolated_vep(dir)
        .arg("-i")
        .arg(records_file(dir, records))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus("GRCh38-everything").join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
        .args(["--tab", "--custom", &source])
        .args(flags)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "vep under {key}:\n{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let text = fs::read_to_string(&out).unwrap();
    let mut lines = text.lines().filter(|line| !line.starts_with("##"));
    let columns: Vec<&str> = lines
        .next()
        .expect("the column line")
        .trim_start_matches('#')
        .split('\t')
        .collect();
    let place = |name: &str| {
        let places: Vec<usize> = (0..columns.len())
            .filter(|&at| columns[at] == name)
            .collect();
        assert_eq!(places.len(), 1, "the column {name} is listed once");
        places[0]
    };
    let places = ["Uploaded_variation", "Allele", "Feature", key].map(place);
    lines
        .map(|row| {
            let cells: Vec<&str> = row.split('\t').collect();
            let [record, allele, feature, cell] = places.map(|at| cells[at].to_string());
            (record, allele, feature, cell)
        })
        .collect()
}

/// The cells of the record of sixteen rows, as `(Feature, key)` in the order
/// of the rows, after a check that every row is the record's one allele.
fn sixteen_cells(key: &str, kind: &str, flags: &[&str]) -> Vec<(String, String)> {
    key_cells(&sixteen_row_record(), key, kind, flags)
        .into_iter()
        .map(|(record, allele, feature, cell)| {
            assert_eq!((record.as_str(), allele.as_str()), ("rs2145814066", "T"));
            (feature, cell)
        })
        .collect()
}

/// `cell` on the row of each of the sixteen transcripts.
fn on_every_transcript(cell: impl Fn(&str) -> String) -> Vec<(String, String)> {
    SIXTEEN_TRANSCRIPTS
        .iter()
        .map(|transcript| (transcript.to_string(), cell(transcript)))
        .collect()
}

/// Under `--check_existing` an `overlap` source named `Existing_variation`
/// lists its record's name and then each co-located identifier once, on every
/// row: on the two rows of a record of two alleles and on the four rows of a
/// record of two alleles and two transcripts.
#[test]
fn a_whole_record_custom_key_named_existing_variation_lists_each_identifier_once() {
    let cells = key_cells(
        &two_allele_records(),
        "Existing_variation",
        "overlap",
        &["--check_existing"],
    );
    let row = |record: &str, allele: &str, feature: &str, cell: &str| {
        (
            record.to_string(),
            allele.to_string(),
            feature.to_string(),
            cell.to_string(),
        )
    };
    assert_eq!(
        cells,
        [
            row("rs1444031311", "C", "-", "exact_2,rs1444031311"),
            row("rs1444031311", "G", "-", "exact_2,rs1444031311"),
            row("alleles_reversed", "G", "-", "exact_2,rs1444031311"),
            row("alleles_reversed", "C", "-", "exact_2,rs1444031311"),
            row(
                "rs766964238",
                "C",
                "ENST00000577708",
                "overlap_2,rs766964238"
            ),
            row(
                "rs766964238",
                "T",
                "ENST00000577708",
                "overlap_2,rs766964238"
            ),
            row(
                "rs766964238",
                "C",
                "ENST00000581792",
                "overlap_2,rs766964238"
            ),
            row(
                "rs766964238",
                "T",
                "ENST00000581792",
                "overlap_2,rs766964238"
            ),
        ]
    );
}

/// Under `--mane_select` an `overlap` source named `MANE` lists its record's
/// name on every row and `MANE_Select` after it on the row of the one
/// transcript of that set.
#[test]
fn a_whole_record_custom_key_named_mane_names_the_set_on_its_own_transcripts_row() {
    assert_eq!(
        sixteen_cells("MANE", "overlap", &["--mane_select"]),
        on_every_transcript(|transcript| {
            if transcript == MANE_SELECT_TRANSCRIPT {
                "multi_7,MANE_Select".to_string()
            } else {
                "multi_7".to_string()
            }
        })
    );
}

/// Under `--check_existing --var_synonyms` an `overlap` source named
/// `VAR_SYNONYMS` lists its record's name and the co-located variant's
/// synonyms once on every row, joined by two hyphens.
#[test]
fn a_whole_record_custom_key_named_var_synonyms_lists_the_synonyms_once_on_every_row() {
    assert_eq!(
        sixteen_cells(
            "VAR_SYNONYMS",
            "overlap",
            &["--check_existing", "--var_synonyms"]
        ),
        on_every_transcript(|_| "multi_7--ClinVar::RCV002273474,VCV001699617".to_string())
    );
}

/// Under `--nearest symbol` an `exact` source named `NEAREST` lists the
/// nearest gene's symbol and then the record of the row's allele once, on
/// every row.
#[test]
fn an_exact_custom_key_named_nearest_lists_the_alleles_record_once_on_every_row() {
    assert_eq!(
        sixteen_cells("NEAREST", "exact", &["--nearest", "symbol"]),
        on_every_transcript(|_| "SON,multi_7".to_string())
    );
}
