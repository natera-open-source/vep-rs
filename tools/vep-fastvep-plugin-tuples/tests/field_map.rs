// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The field map: every malformed map the binary refuses, with its line, and the shipped
//! map, `data/fastvep_field_map.tsv`, against the references and the input records of the
//! repository's golden corpora.

mod common;

use common::{one_cell, Scratch, ERROR, FIELD_MAP, MAP_HEADER};
use std::collections::BTreeSet;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

const GOOD_ROW: &str = "FV_X\t1\tX_ONE\tX\ttext\tfabricated\t-\t-";

fn a_malformed_field_map_is_refused_with_its_line(map: &str, message: &str) {
    let scratch = Scratch::new();
    let custom = scratch.write("map.tsv", format!("{map}\n"));
    let run = one_cell(&scratch, "FV_X", "G|a", &["--field-map", &custom]);
    assert_eq!(run.code, 1);
    assert!(
        run.stderr.starts_with(&format!("{ERROR}{custom}{message}")),
        "{}",
        run.stderr
    );
    assert_eq!(run.stderr.matches('\n').count(), 1);
    assert!(!run.wrote());
}

/// The header line, then `row`.
fn with_header(row: &str) -> String {
    format!("{MAP_HEADER}\n{row}")
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_blank_line() {
    a_malformed_field_map_is_refused_with_its_line("", ":1: the header is not fastvep_column part");
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_column_count() {
    a_malformed_field_map_is_refused_with_its_line(
        &with_header("FV_X\t1\tX_ONE\tX\ttext\tfabricated\t-"),
        ":2: 7 columns, not 8",
    );
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_empty_reference_set() {
    a_malformed_field_map_is_refused_with_its_line(
        &with_header("FV_X\t1\tX_ONE\tX\ttext\tfabricated\t-\t"),
        ":2: an empty fastvep_column, vep_field, list_separator or reference_set",
    );
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_empty_separator() {
    a_malformed_field_map_is_refused_with_its_line(
        &with_header("FV_X\t1\tX_ONE\tX\ttext\tfabricated\t\t-"),
        ":2: an empty fastvep_column, vep_field, list_separator or reference_set",
    );
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_falling_range() {
    a_malformed_field_map_is_refused_with_its_line(
        &with_header(&GOOD_ROW.replace("\t1\t", "\t3-2\t")),
        ":2: part '3-2' must count from 1, and a range must rise",
    );
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_field_twice() {
    a_malformed_field_map_is_refused_with_its_line(
        &format!(
            "{}\n{}",
            with_header(GOOD_ROW),
            GOOD_ROW.replace("\t1\t", "\t2\t")
        ),
        ":3: X_ONE is mapped twice",
    );
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_header() {
    a_malformed_field_map_is_refused_with_its_line(
        &format!("fastvep_column\tpart\n{GOOD_ROW}"),
        ":1: the header is not fastvep_column part",
    );
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_kind() {
    a_malformed_field_map_is_refused_with_its_line(
        &with_header(&GOOD_ROW.replace("\ttext\t", "\tfloat\t")),
        ":2: kind 'float' is not text or number",
    );
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_part() {
    a_malformed_field_map_is_refused_with_its_line(
        &with_header(&GOOD_ROW.replace("\t1\t", "\tfirst\t")),
        ":2: part 'first' is not N, N-M, N.word or N.score",
    );
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_part_zero() {
    a_malformed_field_map_is_refused_with_its_line(
        &with_header(&GOOD_ROW.replace("\t1\t", "\t0\t")),
        ":2: part '0' must count from 1",
    );
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_separator_on_a_half() {
    a_malformed_field_map_is_refused_with_its_line(
        &with_header("FV_X\t1.word\tX_ONE\tX\ttext\tfabricated\t|\t-"),
        ":2: a list separator needs a part written N, not '1.word'",
    );
}

#[test]
fn a_malformed_field_map_is_refused_with_its_line_separator_on_a_range() {
    a_malformed_field_map_is_refused_with_its_line(
        &with_header("FV_X\t1-2\tX_ONE\tX\ttext\tfabricated\t|\t-"),
        ":2: a list separator needs a part written N, not '1-2'",
    );
}

#[test]
fn a_missing_field_map_is_refused() {
    let scratch = Scratch::new();
    let absent = common::text(&scratch.path("absent.tsv"));
    let run = one_cell(&scratch, "FV_X", "G|a", &["--field-map", &absent]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.stderr,
        format!("{ERROR}cannot read the field map {absent}: No such file or directory\n")
    );
}

#[test]
fn a_field_map_that_is_not_utf8_text_is_refused() {
    let scratch = Scratch::new();
    let custom = scratch.write("map.tsv", b"\xff\xfe\n");
    let run = one_cell(&scratch, "FV_X", "G|a", &["--field-map", &custom]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.stderr,
        format!("{ERROR}the field map {custom} is not UTF-8 text\n")
    );
}

#[test]
fn an_empty_field_map_is_refused() {
    let scratch = Scratch::new();
    let custom = scratch.write("map.tsv", b"");
    let run = one_cell(&scratch, "FV_X", "G|a", &["--field-map", &custom]);
    assert_eq!(run.code, 1);
    assert!(run.stderr.starts_with(&format!(
        "{ERROR}{custom}:1: the header is not fastvep_column part"
    )));
}

/// One row of the shipped map, by column name.
struct MapRow {
    cells: Vec<String>,
}

impl MapRow {
    fn get(&self, name: &str) -> &str {
        let at = MAP_HEADER
            .split('\t')
            .position(|column| column == name)
            .expect("a column of the map");
        &self.cells[at]
    }
}

fn map_rows() -> Vec<MapRow> {
    let text = fs::read_to_string(FIELD_MAP).unwrap();
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some(MAP_HEADER));
    lines
        .map(|line| MapRow {
            cells: line.split('\t').map(str::to_string).collect(),
        })
        .inspect(|row| assert_eq!(row.cells.len(), 8))
        .collect()
}

/// The golden corpora of the repository this tool is part of.
fn golden() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116")
}

/// The entries of `directory` whose name starts with `prefix`, in name order.
fn starting_with(directory: &Path, prefix: &str) -> Vec<PathBuf> {
    let mut found: Vec<PathBuf> = fs::read_dir(directory)
        .unwrap_or_else(|error| panic!("{}: {error}", directory.display()))
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.starts_with(prefix))
        })
        .collect();
    found.sort();
    found
}

/// The field names the header of a gzipped default-format reference declares.
fn header_fields(reference: &Path) -> BTreeSet<String> {
    let file = fs::File::open(reference).unwrap();
    let mut names = BTreeSet::new();
    for line in BufReader::new(flate2::read::MultiGzDecoder::new(file)).lines() {
        let line = line.unwrap();
        if !line.starts_with("##") {
            break;
        }
        if line.contains(" : ") {
            let declared = line.get(3..).unwrap_or_default();
            let name = declared.split(" : ").next().unwrap_or_default();
            names.insert(name.to_string());
        }
    }
    names
}

/// The fields every GRCh38 reference of a source's corpora declares, over each of their
/// expected sets.
fn reference_fields(source: &str) -> BTreeSet<String> {
    let corpus = format!("GRCh38-plugin-{}", source.to_lowercase());
    starting_with(&golden(), &corpus)
        .iter()
        .flat_map(|corpus| starting_with(corpus, "expected"))
        .map(|set| set.join("default.txt.gz"))
        .filter(|reference| reference.is_file())
        .flat_map(|reference| header_fields(&reference))
        .collect()
}

#[test]
fn the_shipped_map_covers_each_source_column_and_no_constraint_z_score() {
    let rows = map_rows();
    let columns: BTreeSet<&str> = rows.iter().map(|row| row.get("fastvep_column")).collect();
    assert_eq!(
        columns,
        BTreeSet::from([
            "FV_ALPHAMISSENSE",
            "FV_REVEL",
            "SpliceAI",
            "FV_PRIMATEAI",
            "FV_DBNSFP",
            "FV_GNOMAD_GENE",
            "FV_CLINVAR",
            "FV_GNOMAD",
        ])
    );
    let kinds: BTreeSet<&str> = rows.iter().map(|row| row.get("kind")).collect();
    assert_eq!(kinds, BTreeSet::from(["text", "number"]));
    let constraint: Vec<&str> = rows
        .iter()
        .filter(|row| row.get("fastvep_column") == "FV_GNOMAD_GENE")
        .map(|row| row.get("part"))
        .collect();
    assert_eq!(constraint, ["1", "2"]);
}

#[test]
fn every_reference_set_of_the_map_declares_its_field_on_both_assemblies() {
    let rows: Vec<MapRow> = map_rows()
        .into_iter()
        .filter(|row| row.get("reference_set") != "-")
        .collect();
    let sources: BTreeSet<&str> = rows.iter().map(|row| row.get("source")).collect();
    assert_eq!(sources, BTreeSet::from(["SpliceAI", "PrimateAI", "dbNSFP"]));
    for row in &rows {
        for assembly in ["GRCh37", "GRCh38"] {
            let corpus = format!("{assembly}-plugin-{}", row.get("source").to_lowercase());
            let references: Vec<PathBuf> = starting_with(&golden(), &corpus)
                .iter()
                .map(|corpus| corpus.join(row.get("reference_set")).join("default.txt.gz"))
                .filter(|reference| reference.is_file())
                .collect();
            let [reference] = references.as_slice() else {
                panic!("{corpus}*: one reference, not {references:?}");
            };
            assert!(
                header_fields(reference).contains(row.get("vep_field")),
                "{assembly} {}",
                row.get("vep_field")
            );
        }
    }
}

#[test]
fn every_plugin_field_of_the_map_is_a_field_of_a_committed_reference() {
    let rows: Vec<MapRow> = map_rows()
        .into_iter()
        .filter(|row| !["ClinVar", "gnomAD"].contains(&row.get("source")))
        .collect();
    let sources: BTreeSet<&str> = rows.iter().map(|row| row.get("source")).collect();
    assert_eq!(sources.len(), 7);
    for row in &rows {
        assert!(
            reference_fields(row.get("source")).contains(row.get("vep_field")),
            "{}",
            row.get("vep_field")
        );
    }
}

#[test]
fn every_custom_field_of_the_map_is_its_short_name_and_a_committed_info_key() {
    let mut keys = BTreeSet::new();
    for assembly in ["GRCh37", "GRCh38"] {
        let records = fs::read_to_string(golden().join(assembly).join("variants.vcf")).unwrap();
        for line in records.lines().filter(|line| !line.starts_with('#')) {
            let info = line.split('\t').nth(7).expect("an INFO column");
            for item in info.split(';') {
                keys.insert(item.split('=').next().unwrap().to_string());
            }
        }
    }
    let rows: Vec<MapRow> = map_rows()
        .into_iter()
        .filter(|row| ["ClinVar", "gnomAD"].contains(&row.get("source")))
        .collect();
    assert_eq!(rows.len(), 19);
    for row in &rows {
        let (short_name, key) = row
            .get("vep_field")
            .split_once('_')
            .expect("a short name, `_` and an INFO key");
        assert_eq!(short_name, row.get("source"), "{}", row.get("vep_field"));
        assert!(keys.contains(key), "{}", row.get("vep_field"));
    }
}
