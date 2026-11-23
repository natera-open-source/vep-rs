// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class data file of this tool, `data/field_classes.tsv`, held to the class data file
//! of `tools/vep-formats-compare`, `data/classes.tsv`.
//!
//! Each tool holds the rules of the comparison's field classes and states what a class
//! reads in a data file of its own, so a class both files carry is one rule in two places.
//! For such a class both files say the same of:
//!
//! - the field names the class reads, as the list is written: the second column here
//!   against `FIELD_CLASSES.values` there, where a class without a place in
//!   `FIELD_CLASSES` reads none;
//! - the three sets a class is in or out of: the row `#transcript_reading` against the
//!   classes with a place in `TRANSCRIPT_READING_CLASSES`, `#sibling_allele` against
//!   `SIBLING_ALLELE_CLASSES` and `#record_reference` against `RECORD_REFERENCE_CLASSES`;
//! - the order: the class rows here follow the places in `FIELD_CLASSES`, and the names of
//!   a `#` row follow the places in its column.
//!
//! Not held: a class one file alone carries, as a row or as a member of a set; the places
//! as numbers, which count the classes of one file; the columns of `data/classes.tsv` that
//! `data/field_classes.tsv` has no row for; and the rules themselves, which each tool's
//! own tests hold.
//!
//! Each file is read as its tool reads it: `src/rules/fields.rs` here, and `column` of
//! `src/classes.rs` there.

use std::collections::BTreeSet;
use std::path::Path;

/// The classes both files carry.
const SHARED_CLASSES: usize = 85;

/// The `#` rows of `data/field_classes.tsv` that name a set of classes, each with the
/// column of `data/classes.tsv` whose places name the same set.
const SETS: [(&str, &str); 3] = [
    ("#transcript_reading", "TRANSCRIPT_READING_CLASSES"),
    ("#sibling_allele", "SIBLING_ALLELE_CLASSES"),
    ("#record_reference", "RECORD_REFERENCE_CLASSES"),
];

/// What one data file says of the classes.
struct Stated {
    /// Every class the file carries.
    names: BTreeSet<String>,
    /// The classes that read fields, in the file's order of them, each with its field names.
    fields: Vec<(String, Vec<String>)>,
    /// The sets of `SETS`, each in the file's order of its classes.
    sets: Vec<Vec<String>>,
}

impl Stated {
    /// The field names `class` reads; none for a class the file gives none.
    fn fields_of(&self, class: &str) -> &[String] {
        self.fields
            .iter()
            .find(|(held, _)| held == class)
            .map_or(&[], |(_, fields)| fields.as_slice())
    }
}

fn read(path: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(path);
    std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

fn owned(items: &[&str]) -> Vec<String> {
    items.iter().map(|item| item.to_string()).collect()
}

/// `data/field_classes.tsv`: a row is a name, a tab and a comma-joined list; a `#` row of
/// `SETS` lists a set's classes, any other `#` row is skipped, and every other row is a
/// class with its field names.
fn plugin_file() -> Stated {
    let text = read("data/field_classes.tsv");
    let rows: Vec<(&str, Vec<&str>)> = text
        .lines()
        .map(|line| {
            let (name, listed) = line.split_once('\t').expect("a tab in every row");
            (name, listed.split(',').collect())
        })
        .collect();
    let fields: Vec<(String, Vec<String>)> = rows
        .iter()
        .filter(|(name, _)| !name.starts_with('#'))
        .map(|(name, listed)| (name.to_string(), owned(listed)))
        .collect();
    let sets = SETS
        .iter()
        .map(|(row, _)| {
            let (_, listed) = rows
                .iter()
                .find(|(name, _)| name == row)
                .unwrap_or_else(|| panic!("data/field_classes.tsv lacks the row {row}"));
            owned(listed)
        })
        .collect();
    Stated {
        names: fields.iter().map(|(name, _)| name.clone()).collect(),
        fields,
        sets,
    }
}

/// `data/classes.tsv`: the first line names the columns and every other line is one class.
/// A table's column holds the class's place in the table, and a cell that is no whole
/// number leaves the class out of the table.
fn format_file() -> Stated {
    let text = read("../vep-formats-compare/data/classes.tsv");
    let mut lines = text
        .lines()
        .map(|line| line.split('\t').collect::<Vec<&str>>());
    let header = lines.next().expect("a header");
    let rows: Vec<Vec<&str>> = lines.collect();
    let find = |name: &str| {
        header
            .iter()
            .position(|cell| *cell == name)
            .unwrap_or_else(|| panic!("data/classes.tsv lacks the column {name}"))
    };
    let placed = |table: &str| -> Vec<&Vec<&str>> {
        let at = find(table);
        let mut found: Vec<(usize, &Vec<&str>)> = rows
            .iter()
            .filter_map(|row| Some((row.get(at)?.parse::<usize>().ok()?, row)))
            .collect();
        found.sort_by_key(|(place, row)| (*place, row[0]));
        found.into_iter().map(|(_, row)| row).collect()
    };
    let values = find("FIELD_CLASSES.values");
    let fields = placed("FIELD_CLASSES")
        .into_iter()
        .map(|row| {
            let listed: Vec<&str> = row
                .get(values)
                .copied()
                .unwrap_or("")
                .split(',')
                .filter(|item| !item.is_empty())
                .collect();
            (row[0].to_string(), owned(&listed))
        })
        .collect();
    let sets = SETS
        .iter()
        .map(|(_, column)| {
            placed(column)
                .into_iter()
                .map(|row| row[0].to_string())
                .collect()
        })
        .collect();
    Stated {
        names: rows.iter().map(|row| row[0].to_string()).collect(),
        fields,
        sets,
    }
}

/// The classes both files carry, in the order of this tool's file.
fn shared<'a>(plugin: &'a Stated, format: &Stated) -> Vec<&'a str> {
    plugin
        .fields
        .iter()
        .map(|(name, _)| name.as_str())
        .filter(|name| format.names.contains(*name))
        .collect()
}

/// `classes` without those `shared` lacks.
fn among<'a>(classes: impl Iterator<Item = &'a String>, shared: &[&str]) -> Vec<&'a str> {
    classes
        .map(String::as_str)
        .filter(|class| shared.contains(class))
        .collect()
}

#[test]
fn the_two_files_carry_the_shared_classes_counted_here() {
    let (plugin, format) = (plugin_file(), format_file());
    let shared = shared(&plugin, &format);
    assert!(!shared.is_empty());
    assert_eq!(shared.len(), SHARED_CLASSES);
}

#[test]
fn a_shared_class_reads_the_same_field_names_in_both_files() {
    let (plugin, format) = (plugin_file(), format_file());
    for class in shared(&plugin, &format) {
        assert_eq!(plugin.fields_of(class), format.fields_of(class), "{class}");
    }
}

#[test]
fn the_shared_classes_stand_in_the_same_order_in_both_files() {
    let (plugin, format) = (plugin_file(), format_file());
    let shared = shared(&plugin, &format);
    let placed = among(format.fields.iter().map(|(name, _)| name), &shared);
    assert_eq!(shared, placed);
}

#[test]
fn a_shared_class_is_in_the_same_sets_in_both_files() {
    let (plugin, format) = (plugin_file(), format_file());
    let shared = shared(&plugin, &format);
    for (at, (row, column)) in SETS.iter().enumerate() {
        assert_eq!(
            among(plugin.sets[at].iter(), &shared),
            among(format.sets[at].iter(), &shared),
            "{row} against {column}"
        );
    }
}
