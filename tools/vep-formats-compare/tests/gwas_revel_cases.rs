// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The GWAS catalog classes, the sibling-allele values and the REVEL table classes held to
//! stored answers: `tests/data/gwas_revel/cases.jsonl` is a draw of calls, one JSON object
//! a line with the function, its arguments and the answer as `json.dumps` writes it (a set
//! sorted, a tuple a list; `{"error": text}` for the comparator's error and
//! `{"raise": line}` for an exception that ends the run: the line that names it, which
//! the tool ends such a run with less the exception's name). A line
//! `{"table": name, "data": hex, "index": hex}` is a
//! table file and its index, written where the cases after it read them; `{DIR}` in an
//! answer stands for the directory they are written under.

mod compare_common;

use compare_common::scratch;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::io::BufRead;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;
use std::path::PathBuf;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::{Error, Uncaught};
use vep_formats_compare::pyjson::{self, Value};
use vep_formats_compare::rules::gwas::{
    gwas_annotated, gwas_pair_named, gwas_site, GwasCatalogPairs, GwasPair,
};
use vep_formats_compare::rules::revel::{
    class_field_values, revel_one_base_site, revel_other_reference_amino_acid_row, revel_record,
    revel_row_keys, revel_row_passes, revel_walk, sibling_allele_value, trimmed_allele_keys,
    AlleleKeys, RevelRow, RevelTable,
};
use vep_formats_compare::rules::{field_class_accepts, RuleContext, Side};

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

fn truth(value: &Value) -> bool {
    match value {
        Value::Bool(flag) => *flag,
        other => panic!("not a flag: {other:?}"),
    }
}

fn optional_flag(value: &Value) -> Option<bool> {
    match value {
        Value::Null => None,
        other => Some(truth(other)),
    }
}

fn list(value: &Value) -> &[Value] {
    value.as_list().expect("a list")
}

fn entry(value: &Value) -> Entry {
    let parts = list(value);
    Entry {
        key: text(&parts[0]).to_string(),
        record: optional_whole(&parts[1]),
        location: optional_text(&parts[2]),
        allele: text(&parts[3]).to_string(),
        feature: text(&parts[4]).to_string(),
        consequence: text(&parts[5]).to_string(),
        fields: list(&parts[6])
            .iter()
            .map(|pair| {
                let pair = list(pair);
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
    list(value)
        .iter()
        .map(|item| text(item).to_string())
        .collect()
}

fn arg<'a>(args: &'a Value, name: &str) -> &'a Value {
    args.get(name)
        .unwrap_or_else(|| panic!("no argument {name}"))
}

fn string(value: impl Into<String>) -> Value {
    Value::Str(value.into())
}

fn optional<T>(value: Option<T>, made: impl FnOnce(T) -> Value) -> Value {
    value.map_or(Value::Null, made)
}

fn strings<S: Into<String>>(values: impl IntoIterator<Item = S>) -> Value {
    Value::List(values.into_iter().map(string).collect())
}

fn object(members: Vec<(&str, Value)>) -> Value {
    Value::Dict(
        members
            .into_iter()
            .map(|(name, member)| (name.to_string(), member))
            .collect(),
    )
}

fn pairs_of(value: &Value) -> BTreeSet<GwasPair> {
    list(value)
        .iter()
        .map(|pair| {
            let pair = texts(pair);
            (
                pair[0].clone(),
                pair[1].clone(),
                pair[2].clone(),
                pair[3].clone(),
            )
        })
        .collect()
}

fn pairs_value(pairs: &BTreeSet<GwasPair>) -> Value {
    Value::List(
        pairs
            .iter()
            .map(|(chromosome, position, allele, study)| {
                strings([chromosome, position, allele, study])
            })
            .collect(),
    )
}

fn catalog_value(pairs: &GwasCatalogPairs) -> Value {
    object(vec![
        ("reference_only", pairs_value(&pairs.reference_only)),
        ("candidate_only", pairs_value(&pairs.candidate_only)),
    ])
}

fn table_row(value: &Value) -> RevelRow {
    let cells = list(value);
    RevelRow {
        position: text(&cells[0]).to_string(),
        reference: optional_text(&cells[1]),
        alternate: optional_text(&cells[2]),
        reference_amino_acid: text(&cells[3]).to_string(),
        alternate_amino_acid: text(&cells[4]).to_string(),
        score: optional_text(&cells[5]),
        transcript_ids: optional_text(&cells[6]),
    }
}

fn row_value(row: RevelRow) -> Value {
    Value::List(vec![
        string(row.position),
        optional(row.reference, Value::Str),
        optional(row.alternate, Value::Str),
        string(row.reference_amino_acid),
        string(row.alternate_amino_acid),
        optional(row.score, Value::Str),
        optional(row.transcript_ids, Value::Str),
    ])
}

fn rows_value(rows: Vec<RevelRow>) -> Value {
    Value::List(rows.into_iter().map(row_value).collect())
}

fn keys_of(value: &Value) -> AlleleKeys {
    list(value)
        .iter()
        .map(|key| {
            let key = list(key);
            (
                text(&key[0]).to_string(),
                text(&key[1]).to_string(),
                whole(&key[2]),
            )
        })
        .collect()
}

fn keys_value(keys: &AlleleKeys) -> Value {
    Value::List(
        keys.iter()
            .map(|(reference, alternate, at)| {
                Value::List(vec![
                    string(reference.as_str()),
                    string(alternate.as_str()),
                    Value::Int(*at),
                ])
            })
            .collect(),
    )
}

/// The files of the cases: the directory they are written under.
struct Planted {
    directory: PathBuf,
}

impl Planted {
    fn table_path(&self, name: &str) -> PathBuf {
        self.directory.join(name).join("revel.tsv.gz")
    }

    /// Writes the table of one line of the cases and its index.
    fn plant(&self, line: &Value) {
        let path = self.table_path(text(arg(line, "table")));
        std::fs::create_dir_all(path.parent().expect("a directory")).expect("a directory");
        if text(arg(line, "kind")) == "dir" {
            std::fs::create_dir(&path).expect("a directory");
        } else if let Some(data) = arg(line, "data").as_str() {
            std::fs::write(&path, unhex(data)).expect("the table is written");
        }
        if let Some(index) = arg(line, "index").as_str() {
            std::fs::write(format!("{}.tbi", path.display()), unhex(index))
                .expect("the index is written");
        }
    }

    /// The table a case names, read as the run that set `no_match` reads it.
    fn table(&self, args: &Value) -> Result<Option<RevelTable>, Error> {
        match arg(args, "table").as_str() {
            Some(name) => {
                RevelTable::new(&self.table_path(name), truth(arg(args, "no_match"))).map(Some)
            }
            None => Ok(None),
        }
    }

    /// `text` with the directory of the files written `{DIR}`.
    fn placeless(&self, text: &str) -> String {
        text.replace(&self.directory.display().to_string(), "{DIR}")
    }
}

fn unhex(digits: &str) -> Vec<u8> {
    (0..digits.len() / 2)
        .map(|at| u8::from_str_radix(&digits[2 * at..2 * at + 2], 16).expect("two hex digits"))
        .collect()
}

/// How a call ended that ends with the comparator's error or ends the run, as the cases
/// write it; any other panic goes on.
fn failure(planted: &Planted, ended: Result<Error, Box<dyn std::any::Any + Send>>) -> Value {
    match ended {
        Ok(Error::Comparator(what)) => object(vec![("error", string(planted.placeless(&what)))]),
        Err(payload) => match payload.downcast_ref::<Uncaught>() {
            Some(Uncaught(what)) => object(vec![("raise", string(planted.placeless(what)))]),
            None => std::panic::resume_unwind(payload),
        },
    }
}

/// A recorded value with every `{"raise": line}` in it as the tool gives it: the line of an
/// uncaught exception without the exception's name before its message.
fn without_names(value: &Value) -> Value {
    match value {
        Value::List(items) => Value::List(items.iter().map(without_names).collect()),
        Value::Dict(members) => match (members.len(), members.get("raise")) {
            (1, Some(Value::Str(line))) => {
                let message = match line.split_once(": ") {
                    Some((name, message)) if name.ends_with("Error") || name.ends_with("error") => {
                        message
                    }
                    _ => line.as_str(),
                };
                object(vec![("raise", string(message))])
            }
            _ => Value::Dict(
                members
                    .iter()
                    .map(|(name, held)| (name.clone(), without_names(held)))
                    .collect(),
            ),
        },
        other => other.clone(),
    }
}

/// A recorded answer as the tool gives it (`without_names`).
fn message_alone(want: &str) -> String {
    if !want.contains("\"raise\"") {
        return want.to_string();
    }
    match pyjson::loads(want) {
        Ok(answer) => pyjson::dumps(&without_names(&answer), &pyjson::DumpOptions::default()),
        Err(_) => want.to_string(),
    }
}

/// The answer of `call` as the cases write one: its value, or how it failed.
fn attempted(planted: &Planted, call: impl FnOnce() -> Result<Value, Error>) -> Value {
    match catch_unwind(AssertUnwindSafe(call)) {
        Ok(Ok(value)) => value,
        Ok(Err(error)) => failure(planted, Ok(error)),
        Err(payload) => failure(planted, Err(payload)),
    }
}

/// The answer of the tool to one case, as the value `json.dumps` is given.
fn answer(planted: &Planted, function: &str, args: &Value) -> Result<Value, Error> {
    let flag = |value: bool| Ok(Value::Bool(value));
    match function {
        "add_row" => {
            let mut pairs = GwasCatalogPairs::new();
            pairs.add_row(
                text(arg(args, "chromosomes")),
                text(arg(args, "positions")),
                text(arg(args, "risk_alleles")),
                text(arg(args, "study")),
                text(arg(args, "snps")),
            );
            Ok(catalog_value(&pairs))
        }
        "from_catalog" => {
            let path = planted.directory.join(text(arg(args, "name")));
            if let Some(held) = arg(args, "text").as_str() {
                std::fs::create_dir_all(path.parent().expect("a directory")).expect("a directory");
                std::fs::write(&path, held).expect("the catalog is written");
            }
            Ok(catalog_value(&GwasCatalogPairs::from_catalog(&path)?))
        }
        "trimmed_allele_keys" => Ok(keys_value(&trimmed_allele_keys(
            text(arg(args, "reference")),
            text(arg(args, "alternate")),
            whole(arg(args, "position")),
        ))),
        "revel_row_keys" => Ok(optional(
            revel_row_keys(&table_row(arg(args, "row"))),
            |keys| {
                Value::List(
                    keys.into_iter()
                        .map(|(reference, alternate, at)| strings([reference, alternate, at]))
                        .collect(),
                )
            },
        )),
        "revel_row_passes" => {
            let table = planted.table(args)?.expect("a table");
            flag(revel_row_passes(
                &table,
                &table_row(arg(args, "row")),
                text(arg(args, "feature")),
                text(arg(args, "alternate_residue")),
            ))
        }
        "revel_walk" => {
            let table = planted.table(args)?.expect("a table");
            let reached = revel_walk(
                &table,
                text(arg(args, "chromosome")),
                whole(arg(args, "start")),
                whole(arg(args, "end")),
                &keys_of(arg(args, "keys")),
                text(arg(args, "feature")),
                text(arg(args, "alternate_residue")),
                arg(args, "reference_residue").as_str(),
            )?;
            Ok(optional(reached, row_value))
        }
        "RevelTable" => {
            let table = planted.table(args)?.expect("a table");
            let read = list(arg(args, "queries"))
                .iter()
                .map(|query| {
                    let query = list(query);
                    attempted(planted, || {
                        table
                            .rows(text(&query[0]), whole(&query[1]), whole(&query[2]))
                            .map(rows_value)
                    })
                })
                .collect();
            let (reference, alternate, aaref, aaalt, score, ids) = table._columns;
            let mut columns: Vec<Value> = [reference, alternate, aaref, aaalt, score]
                .into_iter()
                .map(Value::Int)
                .collect();
            columns.push(optional(ids, Value::Int));
            Ok(object(vec![
                ("columns", Value::List(columns)),
                ("seq", Value::Int(table._sequence_column)),
                ("pos", Value::Int(table._position_column)),
                (
                    "meta",
                    optional(table._meta, |meta| Value::Int(i64::from(u32::from(meta)))),
                ),
                ("rows", Value::List(read)),
            ]))
        }
        _ => {
            let e = entry(arg(args, "e"));
            let held = args.get("rows").and_then(rows);
            let refs: Option<Vec<&Entry>> = held.as_ref().map(|rows| rows.iter().collect());
            let record_rows = refs.as_deref();
            let candidates = args.get("crows").and_then(rows);
            let candidate_refs: Option<Vec<&Entry>> =
                candidates.as_ref().map(|rows| rows.iter().collect());
            let candidate_rows = candidate_refs.as_deref();
            match function {
                "gwas_site" => Ok(optional(
                    gwas_site(&e, record_rows),
                    |(chromosome, position)| strings([chromosome, position]),
                )),
                "gwas_pair_named" => {
                    let site = texts(arg(args, "site"));
                    flag(gwas_pair_named(
                        &e,
                        text(arg(args, "studies")),
                        (&site[0], &site[1]),
                        &pairs_of(arg(args, "pairs")),
                    ))
                }
                "gwas_annotated" => flag(gwas_annotated(&e)),
                "class_field_values" => {
                    let patterns = texts(arg(args, "patterns"));
                    let patterns: Vec<&str> = patterns.iter().map(String::as_str).collect();
                    Ok(Value::List(
                        class_field_values(&e, &patterns)
                            .into_iter()
                            .map(|(name, value)| strings([name, value]))
                            .collect(),
                    ))
                }
                "sibling_allele_value" => {
                    let patterns = texts(arg(args, "patterns"));
                    let patterns: Vec<&str> = patterns.iter().map(String::as_str).collect();
                    flag(sibling_allele_value(
                        &e,
                        &patterns,
                        record_rows.expect("the record's rows"),
                        truth(arg(args, "same_amino_acids")),
                        candidate_rows,
                    ))
                }
                "revel_record" => Ok(optional(
                    revel_record(&e, record_rows)?,
                    |(chromosome, start, end, keys, place)| {
                        Value::List(vec![
                            string(chromosome),
                            Value::Int(start),
                            Value::Int(end),
                            Value::List(keys.iter().map(keys_value).collect()),
                            Value::Int(place),
                        ])
                    },
                )),
                "revel_one_base_site" => Ok(optional(
                    revel_one_base_site(&e),
                    |(chromosome, position)| {
                        Value::List(vec![string(chromosome), Value::Int(position)])
                    },
                )),
                "revel_other_reference_amino_acid_row" => {
                    let a = optional_entry(arg(args, "a"));
                    let table = planted.table(args)?;
                    flag(revel_other_reference_amino_acid_row(
                        &e,
                        a.as_ref(),
                        text(arg(args, "ev")),
                        text(arg(args, "av")),
                        record_rows,
                        table.as_ref(),
                    )?)
                }
                "field_class_accepts" => {
                    let a = optional_entry(arg(args, "a"));
                    let gwas_pairs = match arg(args, "pairs") {
                        Value::Null => None,
                        pairs => Some(GwasCatalogPairs {
                            reference_only: pairs_of(arg(pairs, "reference_only")),
                            candidate_only: pairs_of(arg(pairs, "candidate_only")),
                        }),
                    };
                    let side = Side {
                        gwas_pairs,
                        revel_table: planted.table(args)?,
                        ..Side::default()
                    };
                    let context = RuleContext {
                        record_rows,
                        candidate_rows,
                        declared: &[],
                        assembly: None,
                        dbnsfp_columns: &[],
                        side: &side,
                    };
                    flag(field_class_accepts(
                        text(arg(args, "klass")),
                        &e,
                        text(arg(args, "name")),
                        text(arg(args, "ev")),
                        text(arg(args, "av")),
                        a.as_ref(),
                        &context,
                    )?)
                }
                other => panic!("no function {other}"),
            }
        }
    }
}

fn cut(text: &str) -> String {
    let short: String = text.chars().take(500).collect();
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
    by_function: BTreeMap<String, usize>,
    differing: Vec<String>,
}

fn answer_all(directory: &Path, lines: impl Iterator<Item = String>) -> Tally {
    let planted = Planted {
        directory: directory.to_path_buf(),
    };
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let mut tally = Tally::default();
    for (place, line) in lines.enumerate() {
        let case = pyjson::loads(&line).expect("a case is JSON");
        let Some(function) = case.get("f") else {
            planted.plant(&case);
            continue;
        };
        let function = text(function);
        let args = arg(&case, "args");
        let want = message_alone(text(arg(&case, "want")));
        tally.cases += 1;
        *tally.by_function.entry(function.to_string()).or_default() += 1;
        let dumped = |value: &Value| pyjson::dumps(value, &pyjson::DumpOptions::default());
        let got = match catch_unwind(AssertUnwindSafe(|| answer(&planted, function, args))) {
            Ok(Ok(value)) => dumped(&value),
            Ok(Err(error)) => dumped(&failure(&planted, Ok(error))),
            Err(payload) => {
                if payload.downcast_ref::<Uncaught>().is_some() {
                    dumped(&failure(&planted, Err(payload)))
                } else {
                    "a panic that is no ending of the run's own".to_string()
                }
            }
        };
        if got == want {
            tally.agreed += 1;
        } else {
            tally.differing.push(format!(
                "line {} {function}: want {} got {} of {}",
                place + 1,
                cut(&want),
                cut(&got),
                cut(&dumped(args))
            ));
        }
    }
    std::panic::set_hook(default_hook);
    tally
}

fn report(tally: &Tally) {
    println!("{} cases, {} agreed", tally.cases, tally.agreed);
    for (function, count) in &tally.by_function {
        println!("{count} cases of {function}");
    }
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
fn every_stored_case_is_answered_as_its_stored_answer_has_it() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gwas_revel/cases.jsonl");
    let cases = std::fs::read_to_string(path).expect("the committed cases read");
    let directory = scratch("gwas-revel-cases");
    let tally = answer_all(&directory, cases.lines().map(str::to_string));
    report(&tally);
    assert_eq!(tally.cases, COMMITTED_CASES);
}

/// The cases of `tests/data/gwas_revel/cases.jsonl`.
const COMMITTED_CASES: usize = 440;

/// The whole file of cases, where `VEP_FORMATS_COMPARE_GWAS_REVEL_CASES` names it.
#[test]
#[ignore = "reads the file VEP_FORMATS_COMPARE_GWAS_REVEL_CASES names"]
fn every_case_of_a_named_file_is_answered_as_its_stored_answer_has_it() {
    let path =
        std::env::var("VEP_FORMATS_COMPARE_GWAS_REVEL_CASES").expect("the variable names a file");
    let file = std::io::BufReader::new(std::fs::File::open(path).expect("the cases open"));
    let directory = scratch("gwas-revel-cases-whole");
    let tally = answer_all(&directory, file.lines().map(|line| line.expect("a line")));
    report(&tally);
}
