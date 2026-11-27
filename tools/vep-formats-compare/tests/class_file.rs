// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `data/classes.tsv` as a program outside the tool reads it: by its path, by the names of
//! its header, a cell empty where the row's class is not in the column's table and else the
//! class's place in that table.

use std::collections::BTreeSet;
use std::path::Path;

use vep_formats_compare::classes::{tables, ARCHIVE_DBSNP_SYNONYMS};
use vep_formats_compare::parse::recoder::IDENTIFIER_INDEL_NOTATIONS;

/// The column of the classes a recoder corpus's `masks` may carry.
const RECODER_CORPUS_MASKS: &str = "RECODER_CORPUS_MASKS";

/// The header and the rows of the file, each row with as many cells as the header.
fn rows() -> (Vec<String>, Vec<Vec<String>>) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("data/classes.tsv");
    let text = std::fs::read_to_string(path).expect("the file is read");
    let mut lines = text.lines().map(|line| {
        line.split('\t')
            .map(str::to_string)
            .collect::<Vec<String>>()
    });
    let header = lines.next().expect("a header");
    let rows: Vec<Vec<String>> = lines.collect();
    assert_eq!(header[0], "class");
    for row in &rows {
        assert_eq!(row.len(), header.len(), "{}", row[0]);
    }
    let names: BTreeSet<&str> = rows.iter().map(|row| row[0].as_str()).collect();
    assert_eq!(names.len(), rows.len(), "a class has one row");
    (header, rows)
}

/// The classes of the column `table`, in the order of their places, which count from 0
/// without a gap.
fn column(header: &[String], rows: &[Vec<String>], table: &str) -> Vec<String> {
    let at = header
        .iter()
        .position(|name| name == table)
        .unwrap_or_else(|| panic!("no column {table}"));
    let mut placed: Vec<(usize, String)> = rows
        .iter()
        .filter(|row| !row[at].is_empty())
        .map(|row| (row[at].parse().expect("a place"), row[0].clone()))
        .collect();
    placed.sort();
    let places: Vec<usize> = placed.iter().map(|(place, _)| *place).collect();
    assert_eq!(places, (0..placed.len()).collect::<Vec<_>>(), "{table}");
    placed.into_iter().map(|(_, class)| class).collect()
}

#[test]
fn the_four_tables_of_corrected_behaviours_are_columns_of_the_file() {
    let (header, rows) = rows();
    let all = tables();
    let fields: Vec<&str> = all
        .field_classes()
        .iter()
        .map(|(class, _)| *class)
        .collect();
    let warnings: Vec<&str> = all
        .plugin_warning_classes()
        .iter()
        .map(|(class, _)| *class)
        .collect();
    assert_eq!(column(&header, &rows, "FIELD_CLASSES"), fields);
    assert_eq!(
        column(&header, &rows, "REGISTRY_CLASSES"),
        all.registry_classes()
    );
    assert_eq!(column(&header, &rows, "PLUGIN_WARNING_CLASSES"), warnings);
    assert_eq!(
        column(&header, &rows, "MANIFEST_CLASSES"),
        all.manifest_classes()
    );
    let union: BTreeSet<String> = [
        "FIELD_CLASSES",
        "REGISTRY_CLASSES",
        "PLUGIN_WARNING_CLASSES",
        "MANIFEST_CLASSES",
    ]
    .iter()
    .flat_map(|table| column(&header, &rows, table))
    .collect();
    assert_eq!(union.len(), 110);
}

#[test]
fn the_classes_of_a_recoder_corpus_are_a_column_of_the_file() {
    let (header, rows) = rows();
    let mut expected: Vec<&str> = tables()
        .recoder_masks()
        .iter()
        .map(|(class, _)| *class)
        .collect();
    expected.extend([ARCHIVE_DBSNP_SYNONYMS, IDENTIFIER_INDEL_NOTATIONS]);
    let held = column(&header, &rows, RECODER_CORPUS_MASKS);
    assert_eq!(held, expected);
    assert_eq!(
        held,
        [
            "refseq_transcript_entries",
            "lrg_transcript_entries",
            "archive_dbsnp_synonyms",
            "identifier_indel_transcript_notations"
        ]
    );
}
