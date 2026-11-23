// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class tables of the tool's data file, `data/classes.tsv`.

use vep_formats_compare::classes::{release_numbers, tables, CORRECTED_IN, NAMES};

const TABLE: &str = include_str!("../data/classes.tsv");

#[test]
fn the_data_file_names_every_class_constant_and_no_class_twice() {
    let listed: Vec<&str> = TABLE
        .lines()
        .skip(1)
        .map(|line| line.split('\t').next().unwrap_or(""))
        .collect();
    let mut once = listed.clone();
    once.sort_unstable();
    once.dedup();
    assert_eq!(once.len(), listed.len());
    for (name, value) in NAMES {
        assert!(listed.contains(value), "{name}");
    }
}

#[test]
fn every_table_of_the_data_file_numbers_its_classes_from_zero_without_a_gap() {
    let mut lines = TABLE.lines();
    let header: Vec<&str> = lines
        .next()
        .expect("the file has its column line")
        .split('\t')
        .collect();
    let rows: Vec<Vec<&str>> = lines.map(|line| line.split('\t').collect()).collect();
    for (at, name) in header.iter().enumerate().skip(1) {
        // A column named `<table>.value` or `<table>.values` holds that table's cells, and
        // `CORRECTED_IN` holds a release where the others hold a place.
        if name.contains('.') || *name == CORRECTED_IN {
            continue;
        }
        let mut places: Vec<usize> = rows
            .iter()
            .filter_map(|row| row.get(at)?.parse().ok())
            .collect();
        places.sort_unstable();
        assert!(!places.is_empty(), "{name}");
        assert_eq!(places, (0..places.len()).collect::<Vec<_>>(), "{name}");
    }
}

#[test]
fn every_class_of_a_table_that_has_a_constant_is_spelt_as_the_constant() {
    let all = tables();
    for (class, _) in all.field_classes() {
        assert!(NAMES.iter().any(|(_, value)| value == class), "{class}");
    }
    assert!(all.in_field_classes(all.field_classes()[0].0));
    assert!(!all.in_field_classes(""));
    assert_eq!(all.field_classes_of(""), None);
}

#[test]
fn the_classes_a_release_corrects_are_the_two_of_the_frequency_filter() {
    let all = tables();
    assert_eq!(
        all.corrected_classes(),
        [
            ("frequency_filter_zero_value", "116.3"),
            ("frequency_filter_allele_key_order", "116.3")
        ]
    );
    for (class, release) in all.corrected_classes() {
        assert_eq!(all.corrected_in(class), Some(*release));
        assert_eq!(release_numbers(release), Some((116, 3)));
        assert!(NAMES.iter().any(|(_, value)| value == class), "{class}");
    }
    assert_eq!(all.corrected_in("order_only"), None);
    assert_eq!(all.corrected_in(""), None);
    // The cell of every other row is empty.
    let mut lines = TABLE.lines();
    let header: Vec<&str> = lines.next().expect("a column line").split('\t').collect();
    let at = header
        .iter()
        .position(|name| *name == CORRECTED_IN)
        .expect("the column of the correcting releases");
    let filled = lines
        .filter(|line| {
            line.split('\t')
                .nth(at)
                .is_some_and(|cell| !cell.is_empty())
        })
        .count();
    assert_eq!(filled, all.corrected_classes().len());
}

#[test]
fn a_release_is_two_whole_numbers_joined_by_a_dot() {
    assert_eq!(release_numbers("116.3"), Some((116, 3)));
    assert_eq!(release_numbers("115.2"), Some((115, 2)));
    assert_eq!(release_numbers("116.10"), Some((116, 10)));
    assert!(release_numbers("116.2") < release_numbers("116.3"));
    assert!(release_numbers("116.9") < release_numbers("116.10"));
    assert!(release_numbers("115.9") < release_numbers("116.0"));
    for text in [
        "", "116", "116.", ".3", "116.3.1", "+116.3", "116.-3", "116.x", "116 .3",
    ] {
        assert_eq!(release_numbers(text), None, "{text}");
    }
}
