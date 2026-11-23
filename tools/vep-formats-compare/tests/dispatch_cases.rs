// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The field-class dispatch and its helpers held to stored answers:
//! `tests/data/dispatch/cases.jsonl` is a draw of calls, one JSON object a line with the
//! function, its arguments and the answer as `json.dumps` writes it (a set sorted, a tuple a
//! list, a sibling row its place among the rows; `{"error": text}` for the comparator's
//! error and `{"raise": type}` for an exception that ends the run).

use indexmap::IndexMap;
use std::io::BufRead;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::{Error, Uncaught};
use vep_formats_compare::pyjson::{self, Value};
use vep_formats_compare::rules::loftee;
use vep_formats_compare::rules::nearest_exon::nearest_exon_items;
use vep_formats_compare::rules::record_alleles as ra;
use vep_formats_compare::rules::{field_class_accepts, RuleContext, Side};

fn text(value: &Value) -> &str {
    value.as_str().expect("a text")
}

fn optional_text(value: &Value) -> Option<String> {
    value.as_str().map(str::to_string)
}

fn optional_whole(value: &Value) -> Option<i64> {
    match value {
        Value::Null => None,
        Value::Int(number) => Some(*number),
        other => panic!("not a whole number: {other:?}"),
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
    let parts = value.as_list().expect("an entry");
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

fn rows(value: &Value) -> Option<Vec<Entry>> {
    value.as_list().map(|rows| rows.iter().map(entry).collect())
}

fn texts(value: &Value) -> Vec<String> {
    value
        .as_list()
        .expect("a list")
        .iter()
        .map(|item| text(item).to_string())
        .collect()
}

fn arg<'a>(args: &'a Value, name: &str) -> &'a Value {
    args.get(name)
        .unwrap_or_else(|| panic!("no argument {name}"))
}

fn flag(value: bool) -> Value {
    Value::Bool(value)
}

fn string(value: impl Into<String>) -> Value {
    Value::Str(value.into())
}

fn optional<T>(value: Option<T>, made: impl FnOnce(T) -> Value) -> Value {
    value.map_or(Value::Null, made)
}

fn strings(values: impl IntoIterator<Item = String>) -> Value {
    Value::List(values.into_iter().map(Value::Str).collect())
}

/// The answer of the tool to one case, as the value `json.dumps` is given; `None` for a call
/// of the dispatch that names a side table, which a case does not carry.
fn answer(function: &str, args: &Value) -> Result<Option<Value>, Error> {
    let done = |value: Value| Ok(Some(value));
    if function == "field_class_accepts" {
        if !arg(args, "side").as_list().expect("a list").is_empty() {
            return Ok(None);
        }
        let e = entry(arg(args, "e"));
        let a = optional_entry(arg(args, "a"));
        let record_rows = rows(arg(args, "rows"));
        let candidate_rows = rows(arg(args, "crows"));
        let record_refs: Option<Vec<&Entry>> =
            record_rows.as_ref().map(|rows| rows.iter().collect());
        let candidate_refs: Option<Vec<&Entry>> =
            candidate_rows.as_ref().map(|rows| rows.iter().collect());
        let declared = texts(arg(args, "declared"));
        let declared: Vec<&str> = declared.iter().map(String::as_str).collect();
        let assembly = optional_text(arg(args, "assembly"));
        // The columns a run hands dbNSFP, an argument of the cases that carry one.
        let columns = args.get("columns").map(texts).unwrap_or_default();
        let side = Side::default();
        let context = RuleContext {
            record_rows: record_refs.as_deref(),
            candidate_rows: candidate_refs.as_deref(),
            declared: &declared,
            assembly: assembly.as_deref(),
            dbnsfp_columns: &columns,
            side: &side,
        };
        return done(flag(field_class_accepts(
            text(arg(args, "klass")),
            &e,
            text(arg(args, "name")),
            text(arg(args, "ev")),
            text(arg(args, "av")),
            a.as_ref(),
            &context,
        )?));
    }
    if let Some(reference) = args.get("reference") {
        let alternates = texts(arg(args, "alternates"));
        return done(match function {
            "vcf_alleles_as_spelt" => {
                let (reference, alternates) =
                    ra::vcf_alleles_as_spelt(text(reference), &alternates);
                Value::List(vec![string(reference), strings(alternates)])
            }
            "vcf_first_alternate" => string(ra::vcf_first_alternate(text(reference), &alternates)),
            other => panic!("no function {other}"),
        });
    }
    if let Some(value) = args.get("value") {
        let value = text(value);
        return done(match function {
            "nearest_exon_items" => optional(nearest_exon_items(value), |items| {
                Value::List(
                    items
                        .into_iter()
                        .map(|(exon, distance, side, length)| {
                            strings([exon, distance, side, length])
                        })
                        .collect(),
                )
            }),
            "notation_indel_length" => optional(loftee::notation_indel_length(value), Value::Int),
            "lof_info_items" => optional(loftee::lof_info_items(value), |items| {
                Value::List(
                    items
                        .into_iter()
                        .map(|(key, item)| strings([key, item]))
                        .collect(),
                )
            }),
            "lof_filter_items" => strings(loftee::lof_filter_items(value)),
            other => panic!("no function {other}"),
        });
    }
    let e = entry(arg(args, "e"));
    if let Some(a) = args.get("a") {
        let a = optional_entry(a);
        return done(match function {
            "gerp_totals" => optional(
                loftee::gerp_totals(&e, a.as_ref()),
                |(reference, candidate, items)| {
                    let items: IndexMap<String, Value> = items
                        .into_iter()
                        .map(|(key, item)| (key, Value::Str(item)))
                        .collect();
                    Value::List(vec![
                        Value::Float(reference),
                        Value::Float(candidate),
                        Value::Dict(items),
                    ])
                },
            ),
            "lof_filter_item_side" => optional(
                loftee::lof_filter_item_side(&e, a.as_ref(), text(arg(args, "item"))),
                Value::Str,
            ),
            other => panic!("no function {other}"),
        });
    }
    if function == "cdna_span_closed" {
        return done(flag(loftee::cdna_span_closed(&e)));
    }
    let held = rows(arg(args, "rows"));
    let refs: Option<Vec<&Entry>> = held.as_ref().map(|rows| rows.iter().collect());
    let given = refs.as_deref();
    let always = given.unwrap_or(&[]);
    let pair =
        |lengths: (i64, i64)| Value::List(vec![Value::Int(lengths.0), Value::Int(lengths.1)]);
    done(match function {
        "record_alleles" => strings(ra::record_alleles(&e, always)),
        "sibling_rows" => {
            let found = ra::sibling_rows_of(&e, always);
            assert_eq!(
                ra::sibling_rows(&e, always),
                found.iter().map(|row| (*row).clone()).collect::<Vec<_>>()
            );
            Value::List(
                found
                    .iter()
                    .map(|row| {
                        let place = always.iter().position(|held| std::ptr::eq(*held, *row));
                        Value::Int(place.expect("a sibling is one of the rows") as i64)
                    })
                    .collect(),
            )
        }
        "record_chromosome" => optional(ra::record_chromosome(&e, given), Value::Str),
        "record_indel_row" => flag(ra::record_indel_row(&e, always)?),
        "record_allele_lengths" => optional(ra::record_allele_lengths(&e, given), pair),
        "first_alternate" => optional(ra::first_alternate(&e, always), Value::Str),
        "alternates_in_order" => optional(ra::alternates_in_order(&e, given), strings),
        "later_alternate" => flag(ra::later_alternate(&e, given)),
        "first_allele_shift" => optional(ra::first_allele_shift(&e, given), Value::Int),
        "earlier_allele_shift" => optional(ra::earlier_allele_shift(&e, given), Value::Int),
        other => panic!("no function {other}"),
    })
}

fn cut(text: &str) -> String {
    let short: String = text.chars().take(600).collect();
    if short.len() < text.len() {
        format!("{short}...")
    } else {
        short
    }
}

/// How the cases of one file were answered.
#[derive(Default)]
struct Tally {
    cases: usize,
    agreed: usize,
    /// The cases that name a side table.
    side_tables: usize,
    /// The cases answered otherwise whose line holds a whole number past `i64`: the stored
    /// answer reads a number of any length, the tool the nearest `i64` or the comparator's
    /// error.
    past_i64: usize,
    differing: Vec<String>,
}

/// Whether a line holds a run of 19 or more digits.
fn holds_a_long_number(line: &str) -> bool {
    line.as_bytes()
        .split(|byte| !byte.is_ascii_digit())
        .any(|run| run.len() >= 19)
}

fn answer_all(lines: impl Iterator<Item = String>) -> Tally {
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let mut tally = Tally::default();
    for (place, line) in lines.enumerate() {
        let case = pyjson::loads(&line).expect("a case is JSON");
        let function = text(arg(&case, "f"));
        let args = arg(&case, "args");
        let want = text(arg(&case, "want"));
        tally.cases += 1;
        let dumped = |value: &Value| pyjson::dumps(value, &pyjson::DumpOptions::default());
        let got = match catch_unwind(AssertUnwindSafe(|| answer(function, args))) {
            Ok(Ok(Some(value))) => dumped(&value),
            Ok(Ok(None)) => {
                tally.side_tables += 1;
                continue;
            }
            Ok(Err(Error::Comparator(what))) => {
                let mut error = IndexMap::new();
                error.insert("error".to_string(), Value::Str(what));
                dumped(&Value::Dict(error))
            }
            Err(payload) => {
                match payload.downcast_ref::<Uncaught>() {
                    // A stored exception is known by its type, the tool's by its text.
                    Some(_) if want.starts_with("{\"raise\"") => want.to_string(),
                    Some(Uncaught(what)) => format!("uncaught: {what}"),
                    None => "a panic that is no exception of the stored answers".to_string(),
                }
            }
        };
        if got == want {
            tally.agreed += 1;
        } else if holds_a_long_number(&line) {
            tally.past_i64 += 1;
        } else {
            tally.differing.push(format!(
                "line {} {function}: want {} got {} of {}",
                place + 1,
                cut(want),
                cut(&got),
                cut(&dumped(args))
            ));
        }
    }
    std::panic::set_hook(default_hook);
    tally
}

fn report(tally: &Tally) {
    println!(
        "{} cases, {} agreed, {} name a side table, {} differ on a number past i64",
        tally.cases, tally.agreed, tally.side_tables, tally.past_i64
    );
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
fn the_committed_cases_are_answered_as_the_file_says() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/dispatch/cases.jsonl");
    let cases = std::fs::read_to_string(path).expect("the committed cases read");
    let tally = answer_all(cases.lines().map(str::to_string));
    assert_eq!(tally.cases, 879);
    assert_eq!((tally.side_tables, tally.past_i64), (0, 0));
    report(&tally);
}

/// The whole file of cases, where `VEP_FORMATS_COMPARE_DISPATCH_CASES` names it.
#[test]
#[ignore = "reads the file VEP_FORMATS_COMPARE_DISPATCH_CASES names"]
fn every_case_of_the_whole_file_is_answered_as_the_file_says() {
    let path =
        std::env::var("VEP_FORMATS_COMPARE_DISPATCH_CASES").expect("the variable names a file");
    let file = std::io::BufReader::new(std::fs::File::open(path).expect("the cases open"));
    let tally = answer_all(file.lines().map(|line| line.expect("a line")));
    report(&tally);
}
