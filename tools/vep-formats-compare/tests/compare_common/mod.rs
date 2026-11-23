// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the comparison share: the committed reference outputs, planted rows and
//! the readers of a report.

#![allow(dead_code)]

use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::path::Path;
use std::path::PathBuf;

use indexmap::IndexMap;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::counter::Counter;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::entry::Fields;
use vep_formats_compare::input::input_records;
use vep_formats_compare::input::read_text;
use vep_formats_compare::input::text_lines;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::Side;

/// The records the reference annotated and its consequence rows.
pub const RECORDS: usize = 36;
pub const ROWS: i64 = 977;

pub fn golden() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116")
}

pub fn corpus() -> PathBuf {
    golden().join("GRCh38-release116")
}

/// The committed reference output of `fmt`.
pub fn expected_path(fmt: &str) -> PathBuf {
    let name = match fmt {
        "default" => "default.txt.gz",
        "tab" => "tab.txt.gz",
        "vcf" => "vcf.vcf.gz",
        _ => "json.jsonl.gz",
    };
    corpus().join("expected").join(name)
}

pub fn text_of(fmt: &str) -> String {
    read_text(&expected_path(fmt)).expect("the committed output reads")
}

pub fn inputs() -> Vec<String> {
    input_records(&corpus().join("variants.vcf")).expect("the committed input reads")
}

pub fn manifest() -> Value {
    json(&std::fs::read_to_string(corpus().join("manifest.json")).expect("the manifest reads"))
}

pub fn data_lines(text: &str) -> Vec<String> {
    text_lines(text)
        .into_iter()
        .filter(|line| !line.starts_with('#'))
        .collect()
}

pub fn header_lines(text: &str) -> Vec<String> {
    text_lines(text)
        .into_iter()
        .filter(|line| line.starts_with('#'))
        .collect()
}

/// The header lines of `text` over `rows`.
pub fn rebuild(text: &str, rows: &[String]) -> String {
    let mut lines = header_lines(text);
    lines.extend(rows.iter().cloned());
    lines.join("\n") + "\n"
}

pub fn masks(classes: &[&str]) -> BTreeSet<String> {
    classes.iter().map(|klass| klass.to_string()).collect()
}

pub fn no_documented() -> Documented {
    Documented::new()
}

/// A `Documented` of no manifest whose mask set names `classes`.
pub fn documented_with(classes: &[&str]) -> Documented {
    Documented::from_manifest(None, &masks(classes), None, None, None)
        .expect("no manifest documents nothing")
}

/// The `Documented` of `manifest` under `classes`, no run named.
pub fn from_manifest(manifest: &Value, classes: &[&str]) -> Documented {
    Documented::from_manifest(Some(manifest), &masks(classes), None, None, None)
        .expect("the manifest reads")
}

/// The `Documented` of the run `run` of `manifest` under `classes`.
pub fn from_run(manifest: &Value, run: &str, classes: &[&str]) -> Documented {
    Documented::from_manifest(Some(manifest), &masks(classes), None, Some(run), None)
        .expect("the manifest reads")
}

/// The report of two entry lists compared with nothing named beside them.
pub fn compare(expected: &[Entry], actual: &[Entry], documented: &mut Documented) -> Value {
    let side = Side::default();
    compare_entries(expected, actual, documented, &CompareOptions::new(&side))
        .expect("the entries compare")
}

/// The report of two default-format texts.
pub fn compare_texts(reference: &str, candidate: &str, documented: &mut Documented) -> Value {
    let expected = parse_default(reference).expect("the reference parses");
    let actual = parse_default(candidate).expect("the candidate parses");
    compare(&expected, &actual, documented)
}

pub fn json(text: &str) -> Value {
    pyjson::loads(text).expect("the text is JSON")
}

/// The member of `value` under the names of `path`, one object deep a name.
pub fn at<'a>(value: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(value, |value, name| {
        value
            .get(name)
            .unwrap_or_else(|| panic!("no member {name} of {path:?}"))
    })
}

pub fn verdict(report: &Value) -> &str {
    at(report, &["verdict"])
        .as_str()
        .expect("a verdict is text")
}

pub fn texts(value: &Value) -> Vec<&str> {
    value
        .as_list()
        .expect("a list")
        .iter()
        .map(|item| item.as_str().expect("a text"))
        .collect()
}

pub fn reasons(report: &Value) -> Vec<&str> {
    texts(at(report, &["reasons"]))
}

/// The whole number under the names of `path`.
pub fn number(report: &Value, path: &[&str]) -> i64 {
    match at(report, path) {
        Value::Int(value) => *value,
        other => panic!("{path:?} is {other:?}"),
    }
}

/// The counts of `counter` other than zero, as a count map equal to another holds them.
pub fn counts(counter: &Counter) -> BTreeMap<String, i64> {
    counter
        .iter()
        .filter(|(_, counted)| *counted != 0)
        .map(|(name, counted)| (name.to_string(), counted))
        .collect()
}

pub fn count_map(items: &[(&str, i64)]) -> BTreeMap<String, i64> {
    items
        .iter()
        .map(|(name, counted)| (name.to_string(), *counted))
        .collect()
}

pub fn fields(items: &[(&str, &str)]) -> Fields {
    items
        .iter()
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect()
}

/// An entry of the given key and fields.
pub fn entry(
    key: &str,
    record: Option<i64>,
    location: Option<&str>,
    allele: &str,
    feature: &str,
    consequence: &str,
    held: &[(&str, &str)],
) -> Entry {
    Entry {
        key: key.to_string(),
        record,
        location: location.map(str::to_string),
        allele: allele.to_string(),
        feature: feature.to_string(),
        consequence: consequence.to_string(),
        fields: fields(held),
        ..Entry::default()
    }
}

/// `value` with its member `name` set to `member`.
pub fn with_member(value: &Value, name: &str, member: Value) -> Value {
    let mut members: IndexMap<String, Value> = value.as_dict().expect("an object").clone();
    members.insert(name.to_string(), member);
    Value::Dict(members)
}

/// The header line of a default-format output and one row of it.
pub const DEFAULT_HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";

/// One row of a default-format output under `DEFAULT_HEADER`: a Transcript row of the gene
/// `ENSG1`.
pub fn default_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    consequence: &str,
    extra: &str,
) -> String {
    default_row_of(
        name,
        location,
        allele,
        feature,
        consequence,
        extra,
        "Transcript",
        "ENSG1",
    )
}

#[allow(clippy::too_many_arguments)]
pub fn default_row_of(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    consequence: &str,
    extra: &str,
    feature_type: &str,
    gene: &str,
) -> String {
    format!(
        "{name}\t{location}\t{allele}\t{gene}\t{feature}\t{feature_type}\t{consequence}\t-\t-\t-\t-\t-\t-\t{extra}"
    )
}

/// A default-format output of `rows` under `DEFAULT_HEADER`.
pub fn default_text<S: AsRef<str>>(rows: &[S]) -> String {
    let mut lines = vec![DEFAULT_HEADER];
    lines.extend(rows.iter().map(AsRef::as_ref));
    lines.join("\n") + "\n"
}

/// The report of two default-format texts under the built-in `classes`, and what was
/// applied.
pub fn compare_default(reference: &str, candidate: &str, classes: &[&str]) -> (Value, Documented) {
    let mut documented = documented_with(classes);
    let report = compare_texts(reference, candidate, &mut documented);
    (report, documented)
}

/// The text `json.dumps` writes for the value `text` spells.
pub fn dumped(text: &str) -> String {
    pyjson::dumps(&json(text), &pyjson::DumpOptions::default())
}

/// An empty directory of the test `name`, under the build's own directory for test files.
pub fn scratch(name: &str) -> PathBuf {
    let directory = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = std::fs::remove_dir_all(&directory);
    std::fs::create_dir_all(&directory).expect("a directory for the test's files");
    directory
}

/// Writes `text` to `name` under `directory` and returns the file's path.
pub fn planted(directory: &Path, name: &str, text: &str) -> PathBuf {
    let path = directory.join(name);
    std::fs::write(&path, text).expect("the file is written");
    path
}

/// The verdict of two default-format texts under the built-in `classes`.
pub fn verdict_of(reference: &str, candidate: &str, classes: &[&str]) -> String {
    verdict(&compare_default(reference, candidate, classes).0).to_string()
}
