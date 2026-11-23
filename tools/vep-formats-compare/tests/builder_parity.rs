// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The items of the comparison that a program classifying the differences of two outputs
//! keeps copies of, held to `tests/data/builder_parity/cases.tsv`.
//!
//! The table is one line a case, three tab-separated cells: the item, its arguments as a
//! JSON list and its answer as a JSON value, in which a tuple is a list, a set a sorted list
//! and an absent answer `null`. An item is a function of the comparison by its
//! name, a constant by its name, a pattern's text as `<NAME>.pattern` and its flags as
//! `<NAME>.flags`, the interpreter's number for them. `notation_inserted_bases` takes its
//! pattern by the pattern's name, and `lone_insertion_gaps` a case of four arguments where
//! its last two are 0. `SpliceWindows.holds` takes a transcript's stable id, its strand,
//! its introns in the transcript's order as `[first, last]` pairs and the first and last
//! of a run of positions of chromosome 1, and answers the positions of the run that lie in
//! two windows of the transcript.
//!
//! A program that keeps a copy of an item asserts the same line against its copy.

use indexmap::IndexMap;
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use vep_formats_compare::entry::Entry;
use vep_formats_compare::pyjson::{self, Value};
use vep_formats_compare::rules::clipped::exonic_run;
use vep_formats_compare::rules::intron_ends as ie;
use vep_formats_compare::rules::shifted_indel::coordinate_intron_offset;
use vep_formats_compare::rules::splice_windows::SpliceWindows;
use vep_formats_compare::self_parity::reverse_complement;

/// The items of the table and the cases of each.
const ITEMS: &[(&str, usize)] = &[
    ("CLIPPED_NOTATION_RE.pattern", 1),
    ("CLIPPED_NOTATION_RE.flags", 1),
    ("coordinate_intron_offset", 14),
    ("exonic_run", 196),
    ("SpliceWindows.holds", 12),
    ("INTRON_END_COPY_RE.pattern", 1),
    ("INTRON_END_DUPLICATION_RE.pattern", 1),
    ("INTRON_END_INSERTION_RE.pattern", 1),
    ("INTRON_END_PAIRS", 1),
    ("SHORT_INTRON_BASES", 1),
    ("transcript_intron_end_gaps", 18),
    ("intron_end_insertion_form", 13),
    ("copy_of_the_intron_last_base", 13),
    ("duplication_of_the_intron_last_bases", 13),
    ("insertion_that_cannot_end_the_intron", 13),
    ("insertion_that_ends_in_an_intron_end_pair", 13),
    ("INTRON_START_INSERTION_RE.pattern", 1),
    ("coding_transcript_exons", 18),
    ("transcript_intron_start_gaps", 18),
    ("insertion_after_the_exon_last_base", 10),
    ("held_insertion_gaps", 10),
    ("reverse_complement", 4),
    ("notation_inserted_bases", 8),
    ("lone_insertion_gaps", 15),
    ("notation_duplicated_positions", 7),
];

/// The interpreter's number for the flags of a pattern compiled from a text without any:
/// `re.UNICODE`. The tool compiles each pattern from its text alone.
const NO_FLAGS: i64 = 32;

fn text(value: &Value) -> &str {
    value.as_str().expect("a text")
}

fn whole(value: &Value) -> i64 {
    match value {
        Value::Int(number) => *number,
        other => panic!("not a whole number: {other:?}"),
    }
}

fn texts(items: &[&str]) -> Value {
    Value::List(
        items
            .iter()
            .map(|item| Value::Str(item.to_string()))
            .collect(),
    )
}

fn wholes(numbers: impl IntoIterator<Item = i64>) -> Value {
    Value::List(numbers.into_iter().map(Value::Int).collect())
}

fn pairs(spans: &[(i64, i64)]) -> Value {
    Value::List(
        spans
            .iter()
            .map(|(first, last)| wholes([*first, *last]))
            .collect(),
    )
}

fn optional_text(value: Option<String>) -> Value {
    value.map_or(Value::Null, Value::Str)
}

/// The positions of `first` to `last` of chromosome 1 that lie in two windows of the
/// transcript `feature`, whose introns are `gaps` in the transcript's order.
fn held_positions(feature: &str, strand: i64, gaps: &Value, first: i64, last: i64) -> Value {
    let gaps: Vec<(i64, i64)> = gaps
        .as_list()
        .expect("the introns")
        .iter()
        .map(|gap| {
            let gap = gap.as_list().expect("an intron");
            (whole(&gap[0]), whole(&gap[1]))
        })
        .collect();
    let windows = SpliceWindows {
        root: PathBuf::new(),
        input_lines: Vec::new().into(),
        _chromosomes: BTreeSet::from(["1".to_string()]),
        _introns: RefCell::new(IndexMap::from([(
            "1".to_string(),
            IndexMap::from([(feature.to_string(), (strand, gaps))]),
        )])),
    };
    wholes((first..=last).filter(|position| {
        let e = Entry {
            key: format!("1:{position}|G|{feature}|Transcript"),
            location: Some(format!("1:{position}")),
            allele: "G".to_string(),
            feature: feature.to_string(),
            consequence: "x".to_string(),
            ..Entry::default()
        };
        windows.holds(&e).expect("the windows are read")
    }))
}

/// The tool's answer for `item` over `arguments`.
fn answer(item: &str, arguments: &[Value]) -> Value {
    let body = || text(&arguments[0]);
    match item {
        "CLIPPED_NOTATION_RE.pattern" => Value::Str(ie::CLIPPED_NOTATION_RE_PATTERN.to_string()),
        "CLIPPED_NOTATION_RE.flags" => {
            assert_eq!(
                ie::CLIPPED_NOTATION_RE.as_str(),
                ie::CLIPPED_NOTATION_RE_PATTERN
            );
            Value::Int(NO_FLAGS)
        }
        "INTRON_END_COPY_RE.pattern" => Value::Str(ie::INTRON_END_COPY_RE_PATTERN.to_string()),
        "INTRON_END_DUPLICATION_RE.pattern" => {
            Value::Str(ie::INTRON_END_DUPLICATION_RE_PATTERN.to_string())
        }
        "INTRON_END_INSERTION_RE.pattern" => {
            Value::Str(ie::INTRON_END_INSERTION_RE_PATTERN.to_string())
        }
        "INTRON_START_INSERTION_RE.pattern" => {
            Value::Str(ie::INTRON_START_INSERTION_RE_PATTERN.to_string())
        }
        "INTRON_END_PAIRS" => {
            let mut pairs = ie::INTRON_END_PAIRS.to_vec();
            pairs.sort_unstable();
            texts(&pairs)
        }
        "SHORT_INTRON_BASES" => Value::Int(ie::SHORT_INTRON_BASES),
        "coordinate_intron_offset" => Value::Str(coordinate_intron_offset(body())),
        "exonic_run" => exonic_run(body(), text(&arguments[1]))
            .expect("two coordinates")
            .map_or(Value::Null, Value::Int),
        "SpliceWindows.holds" => held_positions(
            body(),
            whole(&arguments[1]),
            &arguments[2],
            whole(&arguments[3]),
            whole(&arguments[4]),
        ),
        "transcript_intron_end_gaps" => wholes(ie::transcript_intron_end_gaps(&arguments[0])),
        "transcript_intron_start_gaps" => wholes(ie::transcript_intron_start_gaps(&arguments[0])),
        "coding_transcript_exons" => match ie::coding_transcript_exons(&arguments[0]) {
            Some((forward, exons, low, high)) => Value::List(vec![
                Value::Bool(forward),
                pairs(&exons),
                Value::Int(low),
                Value::Int(high),
            ]),
            None => Value::Null,
        },
        "intron_end_insertion_form" => optional_text(ie::intron_end_insertion_form(body())),
        "copy_of_the_intron_last_base" => Value::Bool(ie::copy_of_the_intron_last_base(body())),
        "duplication_of_the_intron_last_bases" => {
            Value::Bool(ie::duplication_of_the_intron_last_bases(body()))
        }
        "insertion_that_cannot_end_the_intron" => {
            Value::Bool(ie::insertion_that_cannot_end_the_intron(body()))
        }
        "insertion_that_ends_in_an_intron_end_pair" => {
            Value::Bool(ie::insertion_that_ends_in_an_intron_end_pair(body()))
        }
        "insertion_after_the_exon_last_base" => {
            Value::Bool(ie::insertion_after_the_exon_last_base(body()))
        }
        "held_insertion_gaps" => wholes(ie::held_insertion_gaps(
            whole(&arguments[0]),
            whole(&arguments[1]),
            text(&arguments[2]),
            text(&arguments[3]),
            text(&arguments[4]),
        )),
        "reverse_complement" => Value::Str(reverse_complement(body())),
        "notation_inserted_bases" => {
            let pattern = match body() {
                "INTRON_END_INSERTION_RE" => &*ie::INTRON_END_INSERTION_RE,
                "INTRON_START_INSERTION_RE" => &*ie::INTRON_START_INSERTION_RE,
                other => panic!("no pattern {other}"),
            };
            Value::Str(ie::notation_inserted_bases(pattern, text(&arguments[1])))
        }
        "lone_insertion_gaps" => wholes(ie::lone_insertion_gaps(
            whole(&arguments[0]),
            whole(&arguments[1]),
            text(&arguments[2]),
            text(&arguments[3]),
            arguments.get(4).map_or(0, whole),
            arguments.get(5).map_or(0, whole),
        )),
        "notation_duplicated_positions" => Value::Int(ie::notation_duplicated_positions(body())),
        other => panic!("the table names an item the test does not hold: {other}"),
    }
}

#[test]
fn every_item_answers_its_cases_as_the_table_states() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/builder_parity/cases.tsv");
    let table = std::fs::read_to_string(path).expect("the table is read");
    let mut cases: BTreeMap<&str, usize> = BTreeMap::new();
    for (at, line) in table.lines().enumerate() {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), 3, "line {}", at + 1);
        let arguments = pyjson::loads(cells[1]).expect("the arguments are JSON");
        let stated = pyjson::loads(cells[2]).expect("the answer is JSON");
        let arguments = arguments.as_list().expect("the arguments are a list");
        assert_eq!(
            answer(cells[0], arguments),
            stated,
            "line {}: {} {}",
            at + 1,
            cells[0],
            cells[1]
        );
        *cases.entry(cells[0]).or_default() += 1;
    }
    assert_eq!(cases, ITEMS.iter().copied().collect::<BTreeMap<_, _>>());
}
