// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the shifted-indel and edited-model classes share: planted outputs of the
//! tab, VCF and JSON formats and the readers of a comparison under one class.

#![allow(dead_code)]

use std::collections::BTreeMap;

use crate::compare_common::*;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::field_class_accepts;
use vep_formats_compare::rules::RuleContext;
use vep_formats_compare::rules::Side;

pub const CSQ_DESCRIPTION: &str = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: ";

pub fn csq_header(fields: &str) -> String {
    format!("{CSQ_DESCRIPTION}{fields}\">\n")
}

/// The header of a tab-format output whose columns after the first thirteen are `columns`.
pub fn tab_header(columns: &[&str]) -> String {
    let named: String = columns.iter().map(|column| format!("\t{column}")).collect();
    DEFAULT_HEADER.replace("\tExtra", &named)
}

/// The input record `line` as an output line: its INFO column, the last, left to follow.
pub fn record_columns(line: &str) -> &str {
    &line[..line.len() - 1]
}

pub fn text(value: &str) -> Value {
    Value::Str(value.to_string())
}

pub fn object(members: &[(&str, Value)]) -> Value {
    Value::Dict(
        members
            .iter()
            .map(|(name, member)| (name.to_string(), member.clone()))
            .collect(),
    )
}

/// `value` with each of `members` set, a member it holds keeping its place.
pub fn with_members(value: &Value, members: &[(&str, Value)]) -> Value {
    members.iter().fold(value.clone(), |held, (name, member)| {
        with_member(&held, name, member.clone())
    })
}

/// The line `json.dumps` writes for a record of one transcript consequence.
pub fn json_line(record: &[(&str, Value)], consequence: Value) -> String {
    let mut members = record.to_vec();
    members.push(("transcript_consequences", Value::List(vec![consequence])));
    pyjson::dumps(&object(&members), &pyjson::DumpOptions::default()) + "\n"
}

pub fn once(klass: &str) -> BTreeMap<String, i64> {
    count_map(&[(klass, 1)])
}

pub fn nothing() -> BTreeMap<String, i64> {
    count_map(&[])
}

/// What the report's column `column` counts as masked under `klass` once.
pub fn masked_once(klass: &str) -> Value {
    object(&[(klass, Value::Int(1))])
}

/// The report's verdict and what was applied of two entry lists under one class.
pub fn under_class(
    klass: &str,
    expected: &[Entry],
    actual: &[Entry],
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(&[klass]);
    let report = compare(expected, actual, &mut documented);
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// The verdict of two default-format texts under one class and what was applied.
pub fn default_under(
    klass: &str,
    reference: &str,
    candidate: &str,
) -> (String, BTreeMap<String, i64>) {
    let (report, documented) = compare_default(reference, candidate, &[klass]);
    (verdict(&report).to_string(), counts(&documented.applied))
}

pub fn tab_entries(text: &str) -> Vec<Entry> {
    parse_tab(text).expect("the text parses")
}

pub fn vcf_entries(text: &str, inputs: &[String]) -> Vec<Entry> {
    parse_vcf(text, inputs).expect("the text parses")
}

pub fn json_entries(text: &str, inputs: &[String]) -> Vec<Entry> {
    parse_json(text, inputs).expect("the text parses")
}

pub fn default_entry(row: &str) -> Entry {
    let mut entries = parse_default(&default_text(&[row])).expect("the row parses");
    entries.remove(0)
}

/// `value` without its member `name`.
pub fn without_member(value: &Value, name: &str) -> Value {
    let mut members = value.as_dict().expect("an object").clone();
    members.shift_remove(name);
    Value::Dict(members)
}

/// The entry of the transcript `feature` among `entries`.
pub fn transcript_entry(entries: Vec<Entry>, feature: &str) -> Entry {
    entries
        .into_iter()
        .find(|entry| entry.feature == feature)
        .expect("the transcript's entry")
}

/// The verdict of `text` against itself with `to` in the place of `from`, under `classes`.
pub fn swapped(text: &str, from: &str, to: &str, classes: &[&str]) -> String {
    verdict_of(text, &text.replace(from, to), classes)
}

/// Whether `klass` accepts the field `name` of `e` with nothing named beside the two rows.
pub fn class_accepts(
    klass: &str,
    e: &Entry,
    name: &str,
    expected: &str,
    actual: &str,
    a: Option<&Entry>,
) -> bool {
    let side = Side::default();
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    field_class_accepts(klass, e, name, expected, actual, a, &context)
        .expect("the class reads the rows")
}
