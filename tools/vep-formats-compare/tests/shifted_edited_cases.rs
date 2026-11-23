// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The shifted-indel and edited-model functions held to stored answers:
//! `tests/data/shifted_edited/cases.jsonl` is a draw of calls, one JSON object a line with
//! the function (`f`), its arguments (`a`) and the answer (`w`). An entry is
//! `{"entry": [...]}` over its twelve fields in their order, the fields as pairs; an answer
//! is the returned value, `{"error": text}` for the comparator's error and
//! `{"uncaught": text}` for an exception that ends the run.

use indexmap::IndexMap;
use std::io::BufRead;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::{Error, Uncaught};
use vep_formats_compare::pyjson::{self, Value};
use vep_formats_compare::rules::edited_model as edited;
use vep_formats_compare::rules::shifted_indel as shifted;

fn text(value: &Value) -> &str {
    value.as_str().expect("a text")
}

fn optional_text(value: &Value) -> Option<String> {
    value.as_str().map(str::to_string)
}

fn whole(value: &Value) -> i64 {
    match value {
        Value::Int(number) => *number,
        other => panic!("not a whole number: {other:?}"),
    }
}

fn optional_whole(value: &Value) -> Option<i64> {
    match value {
        Value::Null => None,
        other => Some(whole(other)),
    }
}

fn optional_flag(value: &Value) -> Option<bool> {
    match value {
        Value::Null => None,
        Value::Bool(flag) => Some(*flag),
        other => panic!("not a flag: {other:?}"),
    }
}

fn entry(value: &Value) -> Entry {
    let parts = value
        .get("entry")
        .and_then(Value::as_list)
        .expect("an entry");
    Entry {
        key: text(&parts[0]).to_string(),
        record: optional_whole(&parts[1]),
        location: optional_text(&parts[2]),
        allele: text(&parts[3]).to_string(),
        feature: text(&parts[4]).to_string(),
        consequence: text(&parts[5]).to_string(),
        fields: parts[6]
            .as_list()
            .expect("the fields")
            .iter()
            .map(|pair| {
                let pair = pair.as_list().expect("a pair");
                (text(&pair[0]).to_string(), text(&pair[1]).to_string())
            })
            .collect(),
        indel: optional_flag(&parts[7]),
        insertion: optional_flag(&parts[8]),
        reference: optional_text(&parts[9]),
        line: optional_whole(&parts[10]),
        allele_location: optional_text(&parts[11]),
    }
}

fn optional_entry(value: &Value) -> Option<Entry> {
    match value {
        Value::Null => None,
        other => Some(entry(other)),
    }
}

fn tagged(tag: &str, value: Value) -> Value {
    Value::Dict(IndexMap::from([(tag.to_string(), value)]))
}

fn of_text(value: &str) -> Value {
    Value::Str(value.to_string())
}

fn of_result(result: Result<bool, Error>) -> Value {
    match result {
        Ok(flag) => Value::Bool(flag),
        Err(error) => tagged("error", Value::Str(error.to_string())),
    }
}

/// The tool's answer to the call of `function` with the arguments `a`.
fn answer(function: &str, a: &[Value]) -> Value {
    match function {
        "cdna_pair_spans" => Value::List(
            edited::cdna_pair_spans(text(&a[0]))
                .into_iter()
                .map(|span| span.map_or(Value::Null, Value::Int))
                .collect(),
        ),
        "coordinate_intron_offset" => of_text(&shifted::coordinate_intron_offset(text(&a[0]))),
        "cdna_coordinate_order" => match shifted::cdna_coordinate_order(text(&a[0])) {
            Ok((region, position, offset)) => Value::List(vec![
                Value::Int(region),
                Value::Int(position),
                Value::Int(offset),
            ]),
            Err(error) => tagged("error", Value::Str(error.to_string())),
        },
        "shifted_exon_start" => {
            shifted::shifted_exon_start(text(&a[0]), whole(&a[1])).map_or(Value::Null, Value::Int)
        }
        "shared_ends" => {
            let (front, back) = edited::shared_ends(text(&a[0]), text(&a[1]));
            Value::List(vec![Value::Int(front), Value::Int(back)])
        }
        "indel_notation_of" => shifted::indel_notation_of(&entry(&a[0]), text(&a[1]))
            .map_or(Value::Null, |kind| of_text(&kind)),
        "deletion_offsets_follow_notation" => Value::Bool(
            shifted::deletion_offsets_follow_notation(&entry(&a[0]), text(&a[1])),
        ),
        "equal_span_allele_row" => Value::Bool(shifted::equal_span_allele_row(&entry(&a[0]))),
        "shifted_indel_written_span_gate" => of_result(shifted::shifted_indel_written_span_gate(
            &entry(&a[0]),
            optional_entry(&a[1]).as_ref(),
            text(&a[2]),
            text(&a[3]),
            text(&a[4]),
        )),
        "shifted_indel_written_span_protein_gate" => {
            of_result(shifted::shifted_indel_written_span_protein_gate(
                &entry(&a[0]),
                optional_entry(&a[1]).as_ref(),
                text(&a[2]),
                text(&a[3]),
                None,
            ))
        }
        "intron_offset_abutting_insertion" => {
            Value::Bool(shifted::intron_offset_abutting_insertion(
                &entry(&a[0]),
                optional_entry(&a[1]).as_ref(),
                text(&a[2]),
                text(&a[3]),
                text(&a[4]),
            ))
        }
        "edited_model_shifted_insertion_bases" => {
            Value::Bool(edited::edited_model_shifted_insertion_bases(
                &entry(&a[0]),
                optional_entry(&a[1]).as_ref(),
                text(&a[2]),
                text(&a[3]),
            ))
        }
        "edited_model_seam_insertion_replacement" => {
            Value::Bool(edited::edited_model_seam_insertion_replacement(
                &entry(&a[0]),
                optional_entry(&a[1]).as_ref(),
                text(&a[2]),
                text(&a[3]),
            ))
        }
        "edited_model_clipped_alleles" => of_result(edited::edited_model_clipped_alleles(
            &entry(&a[0]),
            optional_entry(&a[1]).as_ref(),
            text(&a[2]),
            text(&a[3]),
            a[4] == Value::Bool(true),
        )),
        "shifted_deletion_clamped_span" => of_result(edited::shifted_deletion_clamped_span(
            &entry(&a[0]),
            text(&a[1]),
            text(&a[2]),
        )),
        "shifted_deletion_start_offsets" => of_result(edited::shifted_deletion_start_offsets(
            &entry(&a[0]),
            &entry(&a[1]),
            text(&a[2]),
            text(&a[3]),
        )),
        other => panic!("no such function: {other}"),
    }
}

/// Whether the line holds a run of 19 or more digits: a whole number past 64 bits.
fn holds_a_number_past_64_bits(line: &str) -> bool {
    let mut run = 0;
    for c in line.chars() {
        if c.is_numeric() {
            run += 1;
            if run >= 19 {
                return true;
            }
        } else {
            run = 0;
        }
    }
    false
}

fn cut(value: &str) -> String {
    let shown: String = value.chars().take(400).collect();
    if shown.len() < value.len() {
        format!("{shown}...")
    } else {
        shown
    }
}

/// How the cases of one file were answered.
#[derive(Default)]
struct Tally {
    cases: usize,
    agreed: usize,
    /// The cases answered otherwise that hold a whole number past 64 bits.
    past_64_bits: usize,
    differing: Vec<String>,
}

fn answer_all(lines: impl Iterator<Item = String>) -> Tally {
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let mut tally = Tally::default();
    let dumped = |value: &Value| pyjson::dumps(value, &pyjson::DumpOptions::default());
    for (place, line) in lines.enumerate() {
        let case = pyjson::loads(&line).expect("a case is JSON");
        let function = text(case.get("f").expect("the function"));
        let arguments = case
            .get("a")
            .and_then(Value::as_list)
            .expect("the arguments");
        let wanted = case.get("w").expect("the answer");
        tally.cases += 1;
        let got = match catch_unwind(AssertUnwindSafe(|| answer(function, arguments))) {
            Ok(value) => value,
            Err(payload) => match payload.downcast_ref::<Uncaught>() {
                Some(Uncaught(what)) => tagged("uncaught", Value::Str(what.clone())),
                None => tagged("panic", Value::Null),
            },
        };
        if &got == wanted {
            tally.agreed += 1;
        } else if holds_a_number_past_64_bits(&line) {
            tally.past_64_bits += 1;
        } else {
            tally.differing.push(format!(
                "line {} {function}: wanted {} got {} of {}",
                place + 1,
                cut(&dumped(wanted)),
                cut(&dumped(&got)),
                cut(&line)
            ));
        }
    }
    std::panic::set_hook(default_hook);
    tally
}

fn report(tally: &Tally) {
    println!("{} cases, {} agreed", tally.cases, tally.agreed);
    for line in tally.differing.iter().take(40) {
        println!("{line}");
    }
    assert!(
        tally.differing.is_empty(),
        "{} cases differ",
        tally.differing.len()
    );
}

#[test]
fn the_committed_cases_have_their_stored_answers() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/shifted_edited/cases.jsonl");
    let cases = std::fs::read_to_string(path).expect("the committed cases read");
    let tally = answer_all(cases.lines().map(str::to_string));
    assert_eq!(tally.cases, 570);
    assert_eq!(tally.past_64_bits, 0);
    report(&tally);
}

/// The whole file of cases, where `VEP_FORMATS_COMPARE_SHIFTED_EDITED_CASES` names it. A
/// whole number past 64 bits is read to the tool's own bound: a case that holds one and is
/// answered otherwise is counted apart.
#[test]
#[ignore = "reads the file VEP_FORMATS_COMPARE_SHIFTED_EDITED_CASES names"]
fn every_case_of_a_named_file_has_its_stored_answer() {
    let path = std::env::var("VEP_FORMATS_COMPARE_SHIFTED_EDITED_CASES")
        .expect("the variable names a file");
    let file = std::fs::File::open(&path).expect("the named cases open");
    let lines = std::io::BufReader::new(file)
        .lines()
        .map(|line| line.expect("a line of the cases"));
    let tally = answer_all(lines);
    println!(
        "{} cases hold a whole number past 64 bits and are answered otherwise",
        tally.past_64_bits
    );
    report(&tally);
}
