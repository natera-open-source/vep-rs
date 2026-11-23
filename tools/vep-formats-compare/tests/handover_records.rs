// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A VCF or JSON comparison through the stream of record groups against the comparison of
//! the whole outputs. Each test compares planted outputs twice and holds the report, the
//! standard output and the exit code of the two runs to the same bytes. The outputs are
//! written from one table of records: per record the first five columns of its input line
//! and its consequence rows, each the CSQ fields of `PLANTED_CSQ` by name, a list written
//! with the VCF format's `&`.

mod compare_common;

use std::cell::Cell;
use std::collections::BTreeMap;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;

use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::cache::THREADS_VARIABLE;
use vep_formats_compare::classes::tables;
use vep_formats_compare::cli::parse_floors;
use vep_formats_compare::consts::CONSEQUENCE_RANK;
use vep_formats_compare::consts::CONSERVATION_KEY_SPELLING;
use vep_formats_compare::consts::MISSPELT_CONSERVATION_KEY;
use vep_formats_compare::input::read_text;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::pyfmt;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::remainder::hand_over;
use vep_formats_compare::remainder::FrontEnd;
use vep_formats_compare::remainder::STREAM;
use vep_formats_compare::remainder::WHOLE_VARIABLE;
use vep_formats_compare::shared::feature_kind;

const PLANTED_CSQ: [&str; 7] = [
    "Allele",
    "Consequence",
    "IMPACT",
    "Feature_type",
    "Feature",
    "NearestExonJB",
    "am_class",
];
const RECORD_FORMATS: [&str; 2] = ["vcf", "json"];
const BASES: [&str; 4] = ["A", "C", "G", "T"];

/// One consequence row of a record: its CSQ fields by name, and the members its JSON item
/// carries as they are.
#[derive(Clone)]
struct Row {
    fields: IndexMap<String, String>,
    members: IndexMap<String, Value>,
}

impl Row {
    fn of(fields: &[(&str, &str)]) -> Row {
        Row {
            fields: fields
                .iter()
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect(),
            members: IndexMap::new(),
        }
    }

    fn field(&self, name: &str) -> &str {
        self.fields.get(name).map_or("", String::as_str)
    }

    fn with(&self, fields: &[(&str, &str)]) -> Row {
        let mut row = self.clone();
        for (name, value) in fields {
            row.fields.insert(name.to_string(), value.to_string());
        }
        row
    }
}

/// One record of the table: the first five columns of its input line, its rows, the INFO
/// items its VCF line writes before the CSQ and the top-level members its JSON object adds.
#[derive(Clone)]
struct Record {
    columns: [String; 5],
    rows: Vec<Row>,
    info: Vec<String>,
    top: IndexMap<String, Value>,
}

impl Record {
    fn with_rows(&self, rows: Vec<Row>) -> Record {
        Record {
            rows,
            ..self.clone()
        }
    }

    fn with_top(&self, name: &str, value: Value) -> Record {
        let mut record = self.clone();
        record.top = IndexMap::from([(name.to_string(), value)]);
        record
    }
}

fn text(value: &str) -> Value {
    Value::Str(value.to_string())
}

fn dict(members: Vec<(&str, Value)>) -> Value {
    Value::Dict(
        members
            .into_iter()
            .map(|(name, value)| (name.to_string(), value))
            .collect(),
    )
}

/// Record `r` of the table with one transcript row for each of `features`.
fn planted_record_of(r: usize, features: &[String]) -> Record {
    let allele = BASES[(r + 1) % 4];
    let nearest = format!("a{r}&b{r}");
    let rows = features.iter().map(|feature| {
        Row::of(&[
            ("Allele", allele),
            ("Consequence", "intron_variant"),
            ("IMPACT", "MODIFIER"),
            ("Feature_type", "Transcript"),
            ("Feature", feature),
            ("NearestExonJB", &nearest),
        ])
    });
    Record {
        columns: [
            "21".to_string(),
            (1000 + r).to_string(),
            format!("rec{r}"),
            BASES[r % 4].to_string(),
            allele.to_string(),
        ],
        rows: rows.collect(),
        info: Vec::new(),
        top: IndexMap::new(),
    }
}

fn planted_record(r: usize) -> Record {
    planted_record_of(r, &["ENST1".to_string(), "ENST2".to_string()])
}

/// Two transcript rows for each of `records` records.
fn planted_records_of(records: usize) -> Vec<Record> {
    (0..records).map(planted_record).collect()
}

fn planted_records() -> Vec<Record> {
    planted_records_of(8)
}

/// `record` with `fields` written over those of its row `at`.
fn with_row(record: &Record, at: usize, fields: &[(&str, &str)]) -> Record {
    let mut changed = record.clone();
    changed.rows[at] = record.rows[at].with(fields);
    changed
}

/// `records` with the record at `at` replaced.
fn with_record(records: &[Record], at: usize, record: Record) -> Vec<Record> {
    let mut planted = records.to_vec();
    planted[at] = record;
    planted
}

/// `records` with `record` put before the one at `at`.
fn with_record_before(records: &[Record], at: usize, record: Record) -> Vec<Record> {
    let mut planted = records.to_vec();
    planted.insert(at, record);
    planted
}

fn input_line(record: &Record) -> String {
    format!("{}\t.\t.\t.", record.columns.join("\t"))
}

/// The input VCF of a table of records.
fn planted_input(records: &[Record]) -> String {
    records
        .iter()
        .map(|record| input_line(record) + "\n")
        .collect()
}

/// The VCF output of a table of records: one line a record, its `info` items and then its
/// rows as the CSQ entries of `fields`.
fn planted_vcf_of(records: &[Record], fields: &[&str]) -> String {
    let mut lines = vec![format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {}\">",
        fields.join("|")
    )];
    for record in records {
        let chunks: Vec<String> = record
            .rows
            .iter()
            .map(|row| {
                let values: Vec<&str> = fields.iter().map(|name| row.field(name)).collect();
                values.join("|")
            })
            .collect();
        let mut items = record.info.clone();
        if !chunks.is_empty() {
            items.push(format!("CSQ={}", chunks.join(",")));
        }
        let info = if items.is_empty() {
            ".".to_string()
        } else {
            items.join(";")
        };
        lines.push(format!("{}\t.\t.\t{info}", record.columns.join("\t")));
    }
    lines.join("\n") + "\n"
}

fn planted_vcf(records: &[Record]) -> String {
    planted_vcf_of(records, &PLANTED_CSQ)
}

/// The JSON object of one record: its `top` members, and its rows as the items of the
/// consequence list their Feature_type names, a NearestExonJB as the plugin's string, an
/// am_class inside the plugin's object and a row's `members` as they are.
fn planted_object(record: &Record) -> IndexMap<String, Value> {
    let [chromosome, position, name, reference, alternate] = &record.columns;
    let position: i64 = position.parse().expect("a position");
    let mut object: IndexMap<String, Value> = IndexMap::new();
    let mut put = |name: &str, value: Value| {
        object.insert(name.to_string(), value);
    };
    put("input", text(&input_line(record)));
    put("id", text(name));
    put("seq_region_name", text(chromosome));
    put("start", Value::Int(position));
    put("end", Value::Int(position));
    put("allele_string", text(&format!("{reference}/{alternate}")));
    put("strand", Value::Int(1));
    for (name, value) in &record.top {
        put(name, value.clone());
    }
    let rank = |term: &str| {
        let ranked = CONSEQUENCE_RANK.iter().find(|(name, _)| *name == term);
        ranked.expect("a ranked term").1
    };
    let terms = record
        .rows
        .iter()
        .flat_map(|row| row.field("Consequence").split('&'));
    let most_severe = terms.min_by_key(|term| rank(term)).unwrap_or("?");
    put("most_severe_consequence", text(most_severe));
    for row in &record.rows {
        let kind = feature_kind(row.field("Feature_type"));
        let terms = row.field("Consequence").split('&').map(text).collect();
        let mut item: IndexMap<String, Value> = IndexMap::new();
        item.insert("variant_allele".to_string(), text(row.field("Allele")));
        item.insert("consequence_terms".to_string(), Value::List(terms));
        item.insert("impact".to_string(), text(row.field("IMPACT")));
        if kind != "intergenic" {
            item.insert(format!("{kind}_id"), text(row.field("Feature")));
        }
        if !row.field("NearestExonJB").is_empty() {
            let nearest = row.field("NearestExonJB").replace('&', ",");
            item.insert("nearestexonjb".to_string(), text(&nearest));
        }
        if !row.field("am_class").is_empty() {
            let plugin = dict(vec![("am_class", text(row.field("am_class")))]);
            item.insert("alphamissense".to_string(), plugin);
        }
        for (name, value) in &row.members {
            item.insert(name.clone(), value.clone());
        }
        let list = object
            .entry(format!("{kind}_consequences"))
            .or_insert_with(|| Value::List(Vec::new()));
        if let Value::List(items) = list {
            items.push(Value::Dict(item));
        }
    }
    object
}

fn dumped_object(object: IndexMap<String, Value>) -> String {
    pyjson::dumps(&Value::Dict(object), &pyjson::DumpOptions::default())
}

/// The JSON output of a table of records, one object a line.
fn planted_json(records: &[Record]) -> String {
    records
        .iter()
        .map(|record| dumped_object(planted_object(record)) + "\n")
        .collect()
}

fn rendered(fmt: &str, records: &[Record]) -> String {
    if fmt == "vcf" {
        planted_vcf(records)
    } else {
        planted_json(records)
    }
}

fn lines_of(text: &str) -> Vec<String> {
    text.lines().map(str::to_string).collect()
}

/// What one comparison came to: the parsed report (`None` where none was written), the
/// standard output, the exit code, what the stream did with the two outputs as one line, or
/// the standard error of a run that ended in an error, and what the comparison of the whole
/// outputs wrote on the standard error.
struct Compared {
    report: Option<Value>,
    out: String,
    code: i32,
    err: String,
    whole_error: String,
}

impl Compared {
    fn report(&self) -> &Value {
        self.report.as_ref().expect("the run wrote a report")
    }

    fn verdict(&self) -> (&str, i32) {
        (self.out.as_str(), self.code)
    }
}

/// The line a handed-over VCF or JSON comparison writes about the stream.
fn kept_back(kept: usize, lines: usize, handed: usize, groups: usize) -> String {
    format!(
        "[compare_vep_formats] {kept} of {lines} reference lines were kept back by {STREAM}, {handed} of {groups} record groups handed over\n"
    )
}

/// Runs one comparison on the whole outputs and through the stream of record groups, and
/// holds the two runs' report, standard output and exit code to the same bytes.
struct BothPaths {
    directory: PathBuf,
    runs: Cell<usize>,
}

impl BothPaths {
    fn new(name: &str) -> Self {
        BothPaths {
            directory: scratch(name),
            runs: Cell::new(0),
        }
    }

    fn files(&self, expected: &Path, actual: &Path, options: &[&str], fmt: &str) -> Compared {
        self.runs.set(self.runs.get() + 1);
        let run = self.runs.get();
        let mut seen: Vec<(Option<Vec<u8>>, String, i32, String)> = Vec::new();
        for name in ["whole", "handed"] {
            let report = self.directory.join(format!("{name}.{run}.json"));
            let mut command = Command::new(env!("CARGO_BIN_EXE_vep-formats-compare"));
            command
                .args(["--format", fmt, "--expected"])
                .arg(expected)
                .arg("--actual")
                .arg(actual)
                .arg("--report")
                .arg(&report)
                .args(options)
                .env(THREADS_VARIABLE, "2");
            if name == "whole" {
                command.env(WHOLE_VARIABLE, "1");
            } else {
                command.env_remove(WHOLE_VARIABLE);
            }
            let done = command.output().expect("the tool runs");
            seen.push((
                std::fs::read(&report).ok(),
                String::from_utf8(done.stdout).expect("the standard output is text"),
                done.status.code().expect("the tool exits"),
                String::from_utf8(done.stderr).expect("the standard error is text"),
            ));
        }
        let (handed, whole) = (seen.pop().expect("two runs"), seen.pop().expect("two runs"));
        assert_eq!(
            (&handed.0, &handed.1, handed.2),
            (&whole.0, &whole.1, whole.2),
            "{}{}",
            handed.3,
            whole.3
        );
        assert!(whole.3.is_empty() || whole.2 == 2, "{}", whole.3);
        let err = if handed.2 == 2 {
            handed.3
        } else {
            assert!(handed.3.is_empty(), "{}", handed.3);
            self.note(expected, actual, options, fmt)
        };
        Compared {
            report: handed
                .0
                .map(|bytes| json(std::str::from_utf8(&bytes).expect("a report is text"))),
            out: handed.1,
            code: handed.2,
            err,
            whole_error: whole.3,
        }
    }

    /// What the stream of record groups did with the two outputs under the options of one
    /// run, as one line: how many lines were kept back and how many groups handed over.
    fn note(&self, expected: &Path, actual: &Path, options: &[&str], fmt: &str) -> String {
        let value_of = |name: &str| {
            options
                .iter()
                .position(|option| *option == name)
                .map(|at| options[at + 1])
        };
        let classes: Vec<&str> = value_of("--masks")
            .unwrap_or_default()
            .split(',')
            .map(pyfmt::strip)
            .filter(|klass| !klass.is_empty())
            .collect();
        let manifest = value_of("--manifest")
            .map(|path| json(&std::fs::read_to_string(path).expect("the manifest reads")));
        let input = value_of("--input").map(Path::new);
        let input_name = input
            .and_then(Path::file_name)
            .map(|name| name.to_string_lossy().into_owned());
        let documented = Documented::from_manifest(
            manifest.as_ref(),
            &masks(&classes),
            input_name.as_deref(),
            None,
            None,
        )
        .expect("the manifest reads");
        let spelling = CONSERVATION_KEY_SPELLING.replace("\\1", "${1}");
        let renames = classes.contains(&"conservation_short_name_spelling");
        let floors: Vec<String> = parse_floors(value_of("--floors"))
            .expect("the floors read")
            .keys()
            .map(|name| match renames {
                true => MISSPELT_CONSERVATION_KEY
                    .replace_all(name, spelling.as_str())
                    .into_owned(),
                false => name.clone(),
            })
            .collect();
        let floors: Vec<&str> = floors.iter().map(String::as_str).collect();
        let mut front_end = FrontEnd::new(STREAM);
        hand_over(
            &mut front_end,
            fmt,
            expected,
            actual,
            &documented,
            options.contains(&"--unordered"),
            value_of("--keep-only"),
            &floors,
            input,
        )
        .expect("the outputs are read");
        format!("[compare_vep_formats] {}\n", front_end.note)
    }

    /// Compares two outputs of the VCF or JSON format beside the input `variants`, written
    /// under `input_name`.
    fn texts_beside(
        &self,
        fmt: &str,
        reference: &str,
        candidate: &str,
        variants: &str,
        options: &[&str],
        input_name: &str,
    ) -> Compared {
        let directory = self.directory.join(format!("pair.{}", self.runs.get() + 1));
        std::fs::create_dir_all(&directory).expect("a directory for the pair");
        let expected = planted(&directory, &format!("expected.{fmt}"), reference);
        let actual = planted(&directory, &format!("actual.{fmt}"), candidate);
        let input = planted(&directory, input_name, variants);
        let input = input.to_str().expect("a path of text");
        let mut all = vec!["--input", input];
        all.extend_from_slice(options);
        self.files(&expected, &actual, &all, fmt)
    }

    fn texts(
        &self,
        fmt: &str,
        reference: &str,
        candidate: &str,
        variants: &str,
        options: &[&str],
    ) -> Compared {
        self.texts_beside(fmt, reference, candidate, variants, options, "variants.vcf")
    }

    /// Compares the VCF or JSON outputs of two tables of records beside the input of
    /// `variants`.
    fn records_beside(
        &self,
        fmt: &str,
        reference: &[Record],
        candidate: &[Record],
        options: &[&str],
        variants: &[Record],
    ) -> Compared {
        self.texts(
            fmt,
            &rendered(fmt, reference),
            &rendered(fmt, candidate),
            &planted_input(variants),
            options,
        )
    }

    /// `records_beside`, the input holding the reference's records.
    fn records(
        &self,
        fmt: &str,
        reference: &[Record],
        candidate: &[Record],
        options: &[&str],
    ) -> Compared {
        self.records_beside(fmt, reference, candidate, options, reference)
    }
}

fn green(fmt: &str) -> String {
    format!("{fmt}: green\n")
}

/// Eight records of two rows written alike: every line is kept back, the report counts
/// each line's record entry and its two consequence entries, and the stream counts lines.
#[test]
fn the_vcf_lines_and_json_objects_written_alike_are_kept_back() {
    let both_paths = BothPaths::new("records_kept_back_alike");
    for fmt in RECORD_FORMATS {
        let records = planted_records();
        let run = both_paths.records(fmt, &records, &records, &[]);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        for count in ["expected_entries", "actual_entries", "compared"] {
            assert_eq!(number(run.report(), &["counts", count]), 24);
        }
        let impact = if fmt == "vcf" { "IMPACT" } else { "impact" };
        for side in ["expected_populated", "actual_populated"] {
            assert_eq!(number(run.report(), &["columns", impact, side]), 16);
        }
        assert_eq!(run.err, kept_back(8, 8, 0, 8));
    }
}

/// A record whose line differs in one field, in the order of a list, or by a row one side
/// alone writes is handed over and compared as on the whole outputs; the other seven lines
/// are kept back and counted.
#[test]
fn a_vcf_line_or_json_object_that_differs_is_handed_over() {
    let both_paths = BothPaths::new("records_differing");
    for fmt in RECORD_FORMATS {
        let records = planted_records();
        let (impact, nearest) = match fmt {
            "vcf" => ("IMPACT", "NearestExonJB"),
            _ => ("impact", "nearestexonjb"),
        };
        let other = with_record(&records, 3, with_row(&records[3], 0, &[("IMPACT", "LOW")]));
        let run = both_paths.records(fmt, &records, &other, &[]);
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            format!("{fmt}: red\n  field {impact}: 1 differing entries\n")
        );
        assert_eq!(
            number(run.report(), &["columns", impact, "expected_populated"]),
            16
        );
        assert_eq!(run.err, kept_back(7, 8, 1, 8));
        let turned = with_row(&records[3], 1, &[("NearestExonJB", "b3&a3")]);
        let reordered = with_record(&records, 3, turned);
        let run = both_paths.records(fmt, &records, &reordered, &["--masks", "order_only"]);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(
            at(run.report(), &["masks", "applied"]),
            &json(r#"{"order_only": 1}"#)
        );
        let run = both_paths.records(fmt, &records, &reordered, &[]);
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            format!("{fmt}: red\n  field {nearest}: 1 differing entries\n")
        );
        let one_row = records[3].with_rows(records[3].rows[..1].to_vec());
        let fewer = with_record(&records, 3, one_row);
        let run = both_paths.records(fmt, &records, &fewer, &[]);
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            format!("{fmt}: red\n  missing in the candidate: 1 entries\n")
        );
        assert_eq!(number(run.report(), &["counts", "expected_entries"]), 24);
        assert_eq!(number(run.report(), &["counts", "compared"]), 23);
        let run = both_paths.records(fmt, &fewer, &records, &[]);
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            format!("{fmt}: red\n  unexpected in the candidate: 1 entries\n")
        );
        assert_eq!(
            at(run.report(), &["one_sided"]),
            &json(r#"{"expected_only": 0, "actual_only": 1}"#)
        );
    }
}

/// Two records that change places are an order difference at the first one's record entry,
/// its place counted in entries among every shared key, the kept-back lines' entries before
/// it included; `--unordered` reads no order; a record's two rows the other way round are
/// an order difference at the first of them.
#[test]
fn the_place_of_an_order_difference_counts_the_entries_of_vcf_and_json() {
    let both_paths = BothPaths::new("records_order_place");
    for fmt in RECORD_FORMATS {
        let records = planted_records();
        let mut swapped = records.clone();
        swapped.swap(2, 3);
        let run = both_paths.records(fmt, &records, &swapped, &[]);
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            format!("{fmt}: red\n  order differs at shared position 6\n")
        );
        assert_eq!(
            at(run.report(), &["order"]),
            &json(r#"{"first_divergence": 6, "expected": "rec2|record", "actual": "rec3|record"}"#)
        );
        assert_eq!(run.err, kept_back(6, 8, 2, 8));
        let run = both_paths.records(fmt, &records, &swapped, &["--unordered"]);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(run.err, kept_back(8, 8, 0, 8));
        let mut rows = records[5].rows.clone();
        rows.reverse();
        let turned = with_record(&records, 5, records[5].with_rows(rows));
        let run = both_paths.records(fmt, &records, &turned, &[]);
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            format!("{fmt}: red\n  order differs at shared position 16\n")
        );
        let kind = if fmt == "vcf" {
            "Transcript"
        } else {
            "transcript"
        };
        assert_eq!(
            at(run.report(), &["order"]),
            &json(&format!(
                r#"{{"first_divergence": 16, "expected": "rec5|G|ENST1|{kind}", "actual": "rec5|G|ENST2|{kind}"}}"#
            ))
        );
    }
}

fn examples(run: &Compared, side: &str) -> Value {
    at(run.report(), &[side, "", "examples"]).clone()
}

fn labels(labels: &[String]) -> Value {
    Value::List(labels.iter().map(|label| text(label)).collect())
}

/// A line no input record matches is keyed by its ordinal in the output, whatever was kept
/// back before it: unexpected where the candidate alone writes it, compared where both write
/// it at one place, and missing and unexpected where they write it at two. A second line of
/// an input record's columns is such a line; two input lines of the same columns are two
/// records, and their two lines are kept back.
#[test]
fn a_line_the_input_lacks_is_labelled_by_its_place_in_the_output() {
    let both_paths = BothPaths::new("records_labelled_by_place");
    for fmt in RECORD_FORMATS {
        let records = planted_records();
        let label = if fmt == "vcf" { "line" } else { "obj" };
        let stray = planted_record(40);
        let more = with_record_before(&records, 5, stray.clone());
        let run = both_paths.records(fmt, &records, &more, &[]);
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            format!("{fmt}: red\n  unexpected in the candidate: 3 entries\n")
        );
        assert_eq!(
            examples(&run, "unexpected"),
            labels(&[format!("{label}5|record")])
        );
        assert_eq!(run.err, kept_back(8, 8, 1, 9));
        let run = both_paths.records_beside(fmt, &more, &more, &[], &records);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(number(run.report(), &["counts", "compared"]), 27);
        assert_eq!(run.err, kept_back(8, 9, 1, 9));
        let moved = with_record_before(&records, 2, stray);
        let run = both_paths.records_beside(fmt, &more, &moved, &[], &records);
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            format!("{fmt}: red\n  missing in the candidate: 3 entries\n  unexpected in the candidate: 3 entries\n")
        );
        assert_eq!(
            examples(&run, "missing"),
            labels(&[format!("{label}5|record")])
        );
        assert_eq!(
            examples(&run, "unexpected"),
            labels(&[format!("{label}2|record")])
        );
        let twice = with_record_before(&records, 4, records[3].clone());
        let run = both_paths.records_beside(fmt, &twice, &twice, &[], &records);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(number(run.report(), &["counts", "compared"]), 27);
        assert_eq!(run.err, kept_back(7, 9, 1, 8));
        let run = both_paths.records_beside(fmt, &twice, &records, &[], &records);
        assert_eq!(run.code, 1);
        assert_eq!(
            examples(&run, "missing"),
            labels(&[format!("{label}4|record")])
        );
        let run = both_paths.records(fmt, &twice, &twice, &[]);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(number(run.report(), &["counts", "compared"]), 27);
        assert_eq!(run.err, kept_back(9, 9, 0, 8));
    }
}

/// Two lines the input lacks share their record entry's key when they stand at one place of
/// the two outputs, whatever their columns: a record written alike before the one and after
/// the other is an order difference at its record entry.
#[test]
fn the_order_of_lines_paired_by_their_place_in_the_output_is_read() {
    let both_paths = BothPaths::new("records_paired_by_place");
    for fmt in RECORD_FORMATS {
        let known = planted_record(0);
        let (first, second, third) = (planted_record(40), planted_record(41), planted_record(42));
        let reference = [known.clone(), first];
        let candidate = [second, third, known.clone()];
        let run = both_paths.records_beside(fmt, &reference, &candidate, &[], &[known]);
        assert_eq!(run.code, 1);
        assert!(run.out.contains("  order differs at shared position 0\n"));
        let label = if fmt == "vcf" { "line" } else { "obj" };
        assert_eq!(
            at(run.report(), &["order"]),
            &json(&format!(
                r#"{{"first_divergence": 0, "expected": "rec0|record", "actual": "{label}1|record"}}"#
            ))
        );
    }
}

fn first_manifest_class() -> &'static str {
    tables().manifest_classes()[0]
}

/// A key the manifest documents on a record is checked for staleness where the two sides
/// agree: the record's line is handed over though it is written alike, and the stale
/// divergence is reported; a candidate that writes the documented set and the documented
/// field value is set aside as documented.
#[test]
fn a_record_the_manifest_documents_is_handed_over_with_equal_values() {
    let both_paths = BothPaths::new("records_documented");
    let klass = first_manifest_class();
    let manifest = planted(
        &both_paths.directory,
        "manifest.json",
        &format!(
            r#"{{"divergences": [{{"location": "21:1005", "allele": "G", "feature": "ENST2", "record_indices": [5],
                 "expected_divergence": "{klass}", "vep_rs_consequence_sets": ["intron_variant,splice_region_variant"]}}],
                "field_divergences": [{{"location": "21:1006", "allele": "T", "feature": "ENST1", "record_indices": [6], "field": "IMPACT",
                 "vep": "MODIFIER", "vep_rs": "LOW", "expected_divergence": "{klass}"}}]}}"#
        ),
    );
    let options = [
        "--manifest",
        manifest.to_str().expect("a path of text"),
        "--masks",
        klass,
    ];
    for fmt in RECORD_FORMATS {
        let records = planted_records();
        let run = both_paths.records(fmt, &records, &records, &options);
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            format!("{fmt}: red\n  documented divergences not observed: 2\n")
        );
        assert_eq!(run.err, kept_back(6, 8, 2, 8));
        let consequence = [("Consequence", "splice_region_variant&intron_variant")];
        let mut documented = with_record(&records, 5, with_row(&records[5], 1, &consequence));
        documented[6] = with_row(&records[6], 0, &[("IMPACT", "LOW")]);
        let run = both_paths.records(fmt, &records, &documented, &options);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(
            number(run.report(), &["counts", "documented_divergences_skipped"]),
            1
        );
        assert_eq!(number(run.report(), &["counts", "compared"]), 23);
    }
}

/// A manifest item is counted as applied when it is read and taken out again where no entry
/// names a record its lists number. With every line kept back the item stays counted beside
/// the input its run annotates, and is taken out beside an input no run of the manifest
/// names; it is taken out as well where the input holds none of the outputs' lines, which
/// are then all handed over.
#[test]
fn the_manifest_items_are_left_uncounted_where_the_whole_outputs_leave_them() {
    let both_paths = BothPaths::new("records_unmet_items");
    let klass = first_manifest_class();
    let manifest = planted(
        &both_paths.directory,
        "manifest.json",
        &format!(
            r#"{{"divergences": [{{"location": "21:9", "allele": "A", "feature": "ENST9", "record_indices": [30],
                 "expected_divergence": "{klass}", "vep_rs_consequence_sets": ["intron_variant"]}}]}}"#
        ),
    );
    let options = [
        "--manifest",
        manifest.to_str().expect("a path of text"),
        "--masks",
        klass,
    ];
    for fmt in RECORD_FORMATS {
        let records = planted_records();
        let run = both_paths.records(fmt, &records, &records, &options);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(
            at(run.report(), &["masks", "applied"]),
            &json(&format!(r#"{{"{klass}": 1}}"#))
        );
        assert_eq!(run.err, kept_back(8, 8, 0, 8));
        let (reference, variants) = (rendered(fmt, &records), planted_input(&records));
        let run = both_paths.texts_beside(
            fmt,
            &reference,
            &reference,
            &variants,
            &options,
            "other.vcf",
        );
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(at(run.report(), &["masks", "applied"]), &json("{}"));
        assert_eq!(
            at(run.report(), &["masks", "unused"]),
            &json(&format!(r#"["{klass}"]"#))
        );
        assert_eq!(run.err, kept_back(8, 8, 0, 8));
        let strays: Vec<Record> = (20..24).map(planted_record).collect();
        let run = both_paths.records_beside(fmt, &strays, &strays, &options, &records);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(at(run.report(), &["masks", "applied"]), &json("{}"));
        assert_eq!(run.err, kept_back(0, 4, 4, 4));
    }
}

fn regulatory_row() -> Row {
    Row::of(&[
        ("Allele", "T"),
        ("Consequence", "regulatory_region_variant"),
        ("IMPACT", "MODIFIER"),
        ("Feature_type", "RegulatoryFeature"),
        ("Feature", "ENSR1"),
    ])
}

/// `records` with `row` after the rows of the record at `at`.
fn with_one_more_row(records: &[Record], at: usize, row: Row) -> Vec<Record> {
    let mut rows = records[at].rows.clone();
    rows.push(row);
    with_record(records, at, records[at].with_rows(rows))
}

/// Under `regulatory_feature_rows` the reference's regulatory row is dropped before the
/// comparison: a candidate without it agrees, and a candidate that writes the same row, the
/// record's line alike on both sides, has one unexpected entry.
#[test]
fn a_vcf_line_or_json_object_holding_a_regulatory_reference_row_is_handed_over() {
    let both_paths = BothPaths::new("records_regulatory");
    for fmt in RECORD_FORMATS {
        let records = planted_records();
        let reference = with_one_more_row(&records, 2, regulatory_row());
        let options = ["--masks", "regulatory_feature_rows"];
        let run = both_paths.records(fmt, &reference, &records, &options);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(
            at(run.report(), &["masks", "applied"]),
            &json(r#"{"regulatory_feature_rows": 1}"#)
        );
        assert_eq!(run.err, kept_back(7, 8, 1, 8));
        let run = both_paths.records(fmt, &reference, &reference, &options);
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            format!("{fmt}: red\n  unexpected in the candidate: 1 entries\n")
        );
        assert_eq!(run.err, kept_back(7, 8, 1, 8));
        let run = both_paths.records(fmt, &reference, &reference, &[]);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(run.err, kept_back(8, 8, 0, 8));
    }
}

/// Under `conservation_short_name_spelling` the reference's key is rewritten before the
/// comparison, in an INFO item of the VCF line and in a member's name of the JSON object: a
/// candidate spelling it right agrees, the rewrite counted, and a candidate that writes the
/// misspelt key too, the line alike on both sides, differs.
#[test]
fn a_vcf_line_or_json_object_holding_the_misspelt_conservation_key_is_handed_over() {
    let both_paths = BothPaths::new("records_misspelt_key");
    for fmt in RECORD_FORMATS {
        let records = planted_records();
        let scored = |key: &str| {
            let mut record = records[4].with_top(key, Value::Float(0.5));
            record.info = vec![format!("{key}=0.5")];
            with_record(&records, 4, record)
        };
        let (misspelt, spelt) = (scored("phlyoP100way"), scored("phyloP100way"));
        let options = ["--masks", "conservation_short_name_spelling"];
        let run = both_paths.records(fmt, &misspelt, &spelt, &options);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(
            at(run.report(), &["masks", "applied"]),
            &json(r#"{"conservation_short_name_spelling": 1}"#)
        );
        assert_eq!(run.err, kept_back(7, 8, 1, 8));
        let run = both_paths.records(fmt, &misspelt, &misspelt, &options);
        let differing = match fmt {
            "vcf" => "  field INFO_without_CSQ: 1 differing entries\n",
            _ => "  field phlyoP100way: 1 differing entries\n  field phyloP100way: 1 differing entries\n",
        };
        assert_eq!(run.code, 1);
        assert_eq!(run.out, format!("{fmt}: red\n{differing}"));
        assert_eq!(run.err, kept_back(7, 8, 1, 8));
        let run = both_paths.records(fmt, &misspelt, &misspelt, &[]);
        assert_eq!(run.verdict(), (green(fmt).as_str(), 0));
        assert_eq!(run.err, kept_back(8, 8, 0, 8));
    }
}

fn every_line_handed_over() -> String {
    format!("[compare_vep_formats] every line was handed over by {STREAM}: ")
}

/// `conservation_short_name_spelling` rewrites the whole reference text: a CSQ header naming
/// the misspelt key on both sides is renamed on the reference's side alone, which no line's
/// own text shows, so every line is handed over.
#[test]
fn every_vcf_line_is_handed_over_where_the_conservation_rewrite_reaches_the_csq_header() {
    let both_paths = BothPaths::new("records_rewritten_header");
    let scored = [("phlyoP100way", "0.5")];
    let records: Vec<Record> = planted_records()
        .iter()
        .map(|record| with_row(&with_row(record, 0, &scored), 1, &scored))
        .collect();
    let mut fields = PLANTED_CSQ.to_vec();
    fields.push("phlyoP100way");
    let text = planted_vcf_of(&records, &fields);
    let variants = planted_input(&records);
    let options = ["--masks", "conservation_short_name_spelling"];
    let run = both_paths.texts("vcf", &text, &text, &variants, &options);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "vcf: red\n  field phlyoP100way: 16 differing entries\n  field phyloP100way: 16 differing entries\n"
    );
    assert_eq!(
        at(run.report(), &["masks", "applied"]),
        &json(r#"{"conservation_short_name_spelling": 1}"#)
    );
    assert!(
        run.err.starts_with(&every_line_handed_over()),
        "{}",
        run.err
    );
    let run = both_paths.texts("vcf", &text, &text, &variants, &[]);
    assert_eq!(run.verdict(), ("vcf: green\n", 0));
    assert_eq!(run.err, kept_back(8, 8, 0, 8));
}

/// A second CSQ header after two data lines names the last two fields the other way round,
/// so the lines after it are read under another list than the lines before it: no line's own
/// text settles its entries, and every line is handed over where it stands.
#[test]
fn every_vcf_line_is_handed_over_under_a_header_line_after_a_data_line() {
    let both_paths = BothPaths::new("records_late_header");
    let records = planted_records_of(4);
    let lines = lines_of(&planted_vcf(&records));
    let second = lines[0].replace("NearestExonJB|am_class", "am_class|NearestExonJB");
    let mut turned = lines[..3].to_vec();
    turned.push(second);
    turned.extend_from_slice(&lines[3..]);
    let text = turned.join("\n") + "\n";
    let run = both_paths.texts("vcf", &text, &text, &planted_input(&records), &[]);
    assert_eq!(run.verdict(), ("vcf: green\n", 0));
    for column in ["NearestExonJB", "am_class"] {
        assert_eq!(
            number(run.report(), &["columns", column, "expected_populated"]),
            4
        );
    }
    assert!(
        run.err.starts_with(&every_line_handed_over()),
        "{}",
        run.err
    );
}

/// Two JSON objects give the same entries with their members in another order, a solidus
/// spelt with its escape, another engine identity, other text after the `input`'s fifth
/// column and blank lines between them: such lines are kept back. So is a line with a number
/// spelt another way that reads as the same number; a whole number spelt with a fraction is
/// another value to the comparison, and its line is handed over.
#[test]
fn a_json_object_is_read_member_by_member() {
    let both_paths = BothPaths::new("records_member_by_member");
    let records: Vec<Record> = planted_records()
        .iter()
        .map(|record| record.with_top("score", Value::Float(1.5)))
        .collect();
    let (reference, variants) = (planted_json(&records), planted_input(&records));
    let mut lines = Vec::new();
    for record in &records {
        let object = planted_object(record);
        let mut turned: IndexMap<String, Value> = object.clone().into_iter().rev().collect();
        let input = format!("{}\tGT\t0|1", input_line(record));
        turned.insert("input".to_string(), text(&input));
        turned.insert("vep_rs_version".to_string(), text("0.0.0"));
        lines.push(dumped_object(turned).replace('/', "\\/"));
    }
    let candidate = lines.join("\n\n") + "\n \n";
    let run = both_paths.texts("json", &reference, &candidate, &variants, &[]);
    assert_eq!(run.verdict(), ("json: green\n", 0));
    assert_eq!(number(run.report(), &["counts", "compared"]), 24);
    assert_eq!(run.err, kept_back(8, 8, 0, 8));
    for respelt in ["1.50", "15e-1"] {
        let candidate = reference.replacen("\"score\": 1.5,", &format!("\"score\": {respelt},"), 1);
        assert_ne!(candidate, reference);
        let run = both_paths.texts("json", &reference, &candidate, &variants, &[]);
        assert_eq!(run.verdict(), ("json: green\n", 0));
        assert_eq!(run.err, kept_back(8, 8, 0, 8));
    }
    let candidate = reference.replacen("\"strand\": 1,", "\"strand\": 1.0,", 1);
    let run = both_paths.texts("json", &reference, &candidate, &variants, &[]);
    assert_eq!(
        run.verdict(),
        ("json: red\n  field strand: 1 differing entries\n", 1)
    );
    assert_eq!(run.err, kept_back(7, 8, 1, 8));
}

const THIRTY_DIGITS: &str = "123456789012345678901234567890";

/// Record `r` of the table with one co-located variant: its frequencies by allele and
/// population, two of them below a ten-thousandth, a whole number and a whole number of
/// thirty digits; each of its rows carries a float.
fn colocated_record(r: usize) -> Record {
    let record = planted_record(r);
    let allele = record.columns[4].clone();
    let frequencies = dict(vec![(
        allele.as_str(),
        dict(vec![
            ("gnomadg", Value::Float(2.568e-5)),
            ("gnomadg_afr", Value::Float(9.391e-6)),
            ("gnomadg_nfe", Value::Float(0.25)),
            ("gnomadg_eas", Value::Float(0.5)),
        ]),
    )]);
    let colocated = dict(vec![
        ("id", text(&format!("rs{r}"))),
        ("somatic", Value::Int(1)),
        ("allele_number", Value::of_int(THIRTY_DIGITS)),
        ("frequencies", frequencies),
    ]);
    let mut record = record.with_top("colocated_variants", Value::List(vec![colocated]));
    for row in &mut record.rows {
        row.members
            .insert("cadd_raw".to_string(), Value::Float(1.25e-5));
    }
    record
}

/// Eight JSON objects with a co-located variant each, the reference spelling a frequency
/// with a padded exponent and the candidate the same value as a positional decimal or with
/// a bare exponent, in the record's own members and in a consequence item, and each side
/// spelling its own number past the range of a float: every line is kept back, and the
/// report counts what the comparison of the whole outputs counts. The nearest lines the
/// comparison reads otherwise are handed over: a frequency of another value, a whole number
/// spelt as a float, a whole number of thirty digits with another last digit, a frequency
/// inside a string. So is a line whose number past the range of a float has the other
/// sign, which the comparison reads as the same text.
#[test]
fn the_json_objects_whose_numbers_are_spelt_two_ways_are_kept_back() {
    let both_paths = BothPaths::new("records_numbers_spelt_two_ways");
    let records: Vec<Record> = (0..8).map(colocated_record).collect();
    let variants = planted_input(&records);
    let reference =
        planted_json(&records).replace("\"gnomadg_eas\": 0.5", "\"gnomadg_eas\": 1e400");
    for (spelt, times) in [
        ("2.568e-05", 8),
        ("9.391e-06", 8),
        ("1.25e-05", 16),
        ("0.25", 8),
        ("1e400", 8),
        (THIRTY_DIGITS, 8),
    ] {
        assert_eq!(reference.matches(spelt).count(), times, "{spelt}");
    }
    let respelt = reference
        .replace("2.568e-05", "0.00002568")
        .replace("9.391e-06", "9.391e-6")
        .replace("1.25e-05", "0.0000125")
        .replace("0.25", "2.5E-1")
        .replace("1e400", "17e999");
    assert_eq!(respelt.matches("e-0").count(), 0);
    for options in [
        &[][..],
        &["--floors", "Existing_variation=8,gnomADg_AF=8,IMPACT=16"][..],
        &["--unordered"][..],
    ] {
        let run = both_paths.texts("json", &reference, &respelt, &variants, options);
        assert_eq!(run.verdict(), ("json: green\n", 0), "{options:?}");
        assert_eq!(number(run.report(), &["counts", "compared"]), 24);
        for (column, populated) in [("colocated_variants", 8), ("cadd_raw", 16)] {
            for side in ["expected_populated", "actual_populated"] {
                assert_eq!(
                    number(run.report(), &["columns", column, side]),
                    populated,
                    "{column}, {options:?}"
                );
            }
        }
        assert_eq!(run.err, kept_back(8, 8, 0, 8), "{options:?}");
    }
    for (written, read_otherwise, field) in [
        ("0.00002568", "0.00002569", "colocated_variants"),
        ("\"somatic\": 1,", "\"somatic\": 1.0,", "colocated_variants"),
        (
            THIRTY_DIGITS,
            "123456789012345678901234567891",
            "colocated_variants",
        ),
        ("9.391e-6", "\"9.391e-6\"", "colocated_variants"),
        ("0.0000125", "0.0000126", "cadd_raw"),
    ] {
        let candidate = respelt.replacen(written, read_otherwise, 1);
        assert_ne!(candidate, respelt, "{written}");
        let run = both_paths.texts("json", &reference, &candidate, &variants, &[]);
        assert_eq!(
            run.verdict(),
            (
                format!("json: red\n  field {field}: 1 differing entries\n").as_str(),
                1
            ),
            "{read_otherwise}"
        );
        assert_eq!(run.err, kept_back(7, 8, 1, 8), "{read_otherwise}");
    }
    let candidate = respelt.replacen("17e999", "-17e999", 1);
    let run = both_paths.texts("json", &reference, &candidate, &variants, &[]);
    assert_eq!(run.verdict(), ("json: green\n", 0));
    assert_eq!(run.err, kept_back(7, 8, 1, 8));
}

/// A JSON object no input record matches is keyed by its line's ordinal in the output,
/// blank lines counted: the same object after two more blank lines on the candidate's side is
/// another key.
#[test]
fn blank_lines_count_in_the_place_of_a_json_object_the_input_lacks() {
    let both_paths = BothPaths::new("records_blank_lines");
    let records = planted_records();
    let mut all = records.clone();
    all.push(planted_record(40));
    let lines = lines_of(&planted_json(&all));
    let spaced = |blank: usize| {
        let mut spaced = vec![String::new(); blank];
        spaced.extend_from_slice(&lines[..4]);
        spaced.push(String::new());
        spaced.extend_from_slice(&lines[4..]);
        spaced.join("\n") + "\n"
    };
    let (reference, candidate) = (spaced(1), spaced(3));
    let variants = planted_input(&records);
    let run = both_paths.texts("json", &reference, &reference, &variants, &[]);
    assert_eq!(run.verdict(), ("json: green\n", 0));
    assert_eq!(run.err, kept_back(8, 9, 1, 9));
    let run = both_paths.texts("json", &reference, &candidate, &variants, &[]);
    assert_eq!(
        run.verdict(),
        (
            "json: red\n  missing in the candidate: 3 entries\n  unexpected in the candidate: 3 entries\n",
            1
        )
    );
    assert_eq!(
        examples(&run, "missing"),
        labels(&["obj10|record".to_string()])
    );
    assert_eq!(
        examples(&run, "unexpected"),
        labels(&["obj12|record".to_string()])
    );
}

/// Eight records of two rows, the rows of the last three carrying an `am_class`, which the
/// JSON format writes inside the plugin's object.
fn floor_records() -> Vec<Record> {
    let mut records = planted_records();
    for record in &mut records[5..] {
        let first = with_row(record, 0, &[("am_class", "likely_benign")]);
        *record = with_row(&first, 1, &[("am_class", "ambiguous")]);
    }
    records
}

/// A JSON floor is counted under the fields its column has in that format, over the
/// consequence entries of every object, the kept-back ones included: a column's own field,
/// the three feature ids, a field only a record entry carries, which counts no consequence
/// entry, and a record-level column, which has no row count.
#[test]
fn the_kept_back_entries_are_counted_under_a_json_floor() {
    let both_paths = BothPaths::new("records_json_floor");
    let records = floor_records();
    let floors = "IMPACT=16,Feature=16,NearestExonJB=17,STRAND=1,Existing_variation=1";
    let run = both_paths.records("json", &records, &records, &["--floors", floors]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "json: red\n  floor NearestExonJB=17: expected populated 16, actual populated 16\n  floor STRAND=1: expected populated 0, actual populated 0\n"
    );
    assert_eq!(
        at(run.report(), &["floors_not_applicable"]),
        &json(r#"["Existing_variation"]"#)
    );
    assert_eq!(
        at(run.report(), &["columns", "IMPACT", "floor_fields"]),
        &json(r#"["impact"]"#)
    );
    assert_eq!(
        at(run.report(), &["columns", "IMPACT", "floor_met"]),
        &Value::Bool(true)
    );
    for side in ["expected_populated", "actual_populated"] {
        assert_eq!(number(run.report(), &["columns", "Feature", side]), 16);
    }
    assert_eq!(run.err, kept_back(8, 8, 0, 8));
    let other = with_record(&records, 3, with_row(&records[3], 0, &[("IMPACT", "LOW")]));
    let floors = ["--floors", "IMPACT=16,Feature=17"];
    let run = both_paths.records("json", &records, &other, &floors);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "json: red\n  field impact: 1 differing entries\n  floor Feature=17: expected populated 16, actual populated 16\n"
    );
    assert_eq!(run.err, kept_back(7, 8, 1, 8));
}

/// A JSON floor whose own field no entry carries is looked for among the members of the
/// plugin objects and the custom sources: the lines holding either are handed over, so the
/// floor is placed and counted as on the whole outputs. A floor no entry carries anywhere is
/// unmet with no entry on either side, and one on a field holding a dot, which counts a
/// member of an object as well, is compared on the whole outputs.
#[test]
fn the_objects_a_json_floor_is_looked_for_in_are_handed_over() {
    let both_paths = BothPaths::new("records_json_floor_objects");
    let records = floor_records();
    let floors = ["--floors", "am_class=6,IMPACT=16"];
    let run = both_paths.records("json", &records, &records, &floors);
    assert_eq!(run.verdict(), ("json: green\n", 0));
    assert_eq!(
        at(run.report(), &["columns", "am_class", "floor_fields"]),
        &json(r#"["alphamissense.am_class"]"#)
    );
    for side in ["expected_populated", "actual_populated"] {
        assert_eq!(number(run.report(), &["columns", "am_class", side]), 6);
    }
    assert_eq!(
        number(run.report(), &["columns", "IMPACT", "expected_populated"]),
        16
    );
    assert_eq!(run.err, kept_back(5, 8, 3, 8));
    let run = both_paths.records("json", &records, &records, &["--floors", "am_class=7"]);
    assert_eq!(
        run.verdict(),
        (
            "json: red\n  floor am_class=7: expected populated 6, actual populated 6\n",
            1
        )
    );
    let source = json(r#"{"mysource": [{"name": "m0"}]}"#);
    let sourced = with_record(
        &records,
        0,
        records[0].with_top("custom_annotations", source),
    );
    let floors = ["--floors", "mysource_SCORE=1,NOPE=1"];
    let run = both_paths.records("json", &sourced, &sourced, &floors);
    assert_eq!(
        run.verdict(),
        (
            "json: red\n  floor NOPE=1: expected populated 0, actual populated 0\n",
            1
        )
    );
    assert_eq!(
        at(run.report(), &["floors_not_applicable"]),
        &json(r#"["mysource_SCORE"]"#)
    );
    assert_eq!(run.err, kept_back(4, 8, 4, 8));
    let floors = ["--floors", "alphamissense.am_class=6"];
    let run = both_paths.records("json", &records, &records, &floors);
    assert_eq!(run.verdict(), ("json: green\n", 0));
    assert_eq!(
        run.err,
        "[compare_vep_formats] the whole outputs were compared: a floor on a field holding a dot is outside the front end\n"
    );
}

/// Under `regulatory_feature_rows` a floor's field is looked for after the reference's
/// regulatory rows are dropped: a field only such a row carries is no place for the floor,
/// which is then found in the plugin objects, their lines handed over.
#[test]
fn a_json_floor_is_looked_for_among_the_entries_the_regulatory_class_keeps() {
    let both_paths = BothPaths::new("records_json_floor_regulatory");
    let records = floor_records();
    let mut regulatory = regulatory_row();
    regulatory
        .members
        .insert("am_class".to_string(), text("benign"));
    let reference = with_one_more_row(&records, 2, regulatory);
    let options = [
        "--masks",
        "regulatory_feature_rows",
        "--floors",
        "am_class=6",
    ];
    let run = both_paths.records("json", &reference, &records, &options);
    assert_eq!(run.verdict(), ("json: green\n", 0));
    assert_eq!(
        at(run.report(), &["columns", "am_class", "floor_fields"]),
        &json(r#"["alphamissense.am_class"]"#)
    );
    assert_eq!(
        number(run.report(), &["columns", "am_class", "expected_populated"]),
        6
    );
    assert_eq!(run.err, kept_back(4, 8, 4, 8));
    let run = both_paths.records("json", &reference, &records, &options[2..]);
    assert_eq!(run.code, 1);
    assert_eq!(
        at(run.report(), &["columns", "am_class", "floor_fields"]),
        &json(r#"["am_class"]"#)
    );
    assert_eq!(
        number(run.report(), &["columns", "am_class", "expected_populated"]),
        1
    );
    assert_eq!(
        number(run.report(), &["columns", "am_class", "actual_populated"]),
        0
    );
}

/// A VCF floor is read against the populated counts of every CSQ entry, the kept-back ones
/// included, and a floor on a field the CSQ header does not declare has no count.
#[test]
fn the_kept_back_entries_are_counted_under_a_floor_in_the_vcf_format() {
    let both_paths = BothPaths::new("records_vcf_floor");
    let records = floor_records();
    let other = with_record(&records, 3, with_row(&records[3], 0, &[("IMPACT", "LOW")]));
    let floors = ["--floors", "IMPACT=16,NearestExonJB=17,am_class=6,SYMBOL=1"];
    let run = both_paths.records("vcf", &records, &other, &floors);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "vcf: red\n  field IMPACT: 1 differing entries\n  floor NearestExonJB=17: expected populated 16, actual populated 16\n"
    );
    assert_eq!(
        at(run.report(), &["floors_not_applicable"]),
        &json(r#"["SYMBOL"]"#)
    );
    assert_eq!(
        at(run.report(), &["columns", "am_class", "floor_met"]),
        &Value::Bool(true)
    );
    assert_eq!(run.err, kept_back(7, 8, 1, 8));
}

/// A line the comparison cannot read ends both comparisons with the same error line, which
/// names the line by its place in the output or by its record, whatever stands before it.
#[test]
fn both_comparisons_end_with_the_error_of_the_whole_outputs_and_its_place_in_them() {
    let both_paths = BothPaths::new("records_error_place");
    for fmt in RECORD_FORMATS {
        let records = planted_records();
        let variants = planted_input(&records);
        let lines = lines_of(&rendered(
            fmt,
            &with_record_before(&records, 5, planted_record(40)),
        ));
        let (head, data) = lines.split_at(usize::from(fmt == "vcf"));
        // A blank line after the third line of a JSON output.
        let joined = |data: &[String]| {
            let mut all = head.to_vec();
            all.extend_from_slice(&data[..3]);
            if fmt == "json" {
                all.push(String::new());
            }
            all.extend_from_slice(&data[3..]);
            all.join("\n") + "\n"
        };
        let reference = joined(data);
        let (listed, unlisted) = (
            "\"transcript_consequences\": [",
            "\"transcript_consequences\": [1, ",
        );
        let broken: Vec<(&str, usize, String)> = if fmt == "vcf" {
            let short = data[6].rsplit_once('\t').expect("a line of columns").0;
            vec![
                ("VCF line has fewer than 8 columns: ", 6, short.to_string()),
                ("CSQ chunk width on line5: ", 5, format!("{}|more", data[5])),
                ("CSQ chunk width on rec3: ", 3, format!("{}|more", data[3])),
            ]
        } else {
            let cut = data[6][..data[6].len() - 1].to_string();
            vec![
                ("JSON line 8 does not parse: ", 6, cut),
                ("JSON line 8 is not an object", 6, format!("[{}]", data[6])),
                (
                    "obj6: a transcript_consequences element is not an object",
                    5,
                    data[5].replace(listed, unlisted),
                ),
                (
                    "rec3: a transcript_consequences element is not an object",
                    3,
                    data[3].replace(listed, unlisted),
                ),
            ]
        };
        for (error, at, line) in broken {
            let mut changed = data.to_vec();
            assert_ne!(changed[at], line);
            changed[at] = line;
            let run = both_paths.texts(fmt, &reference, &joined(&changed), &variants, &[]);
            assert!(run.report.is_none());
            assert_eq!(run.verdict(), ("", 2));
            assert_eq!(run.err, run.whole_error);
            assert!(
                run.err
                    .starts_with(&format!("ERROR: [compare_vep_formats] {error}")),
                "{}",
                run.err
            );
            assert_eq!(run.err.matches('\n').count(), 1);
        }
    }
}

fn kept_some_lines_back(err: &str) -> bool {
    err.contains(" reference lines were kept back by ")
        && !err.starts_with("[compare_vep_formats] 0 of ")
}

/// The committed reference output against itself and against a copy with one value changed,
/// beside the corpus manifest under every class a manifest documents, the regulatory and
/// order classes and floors: both paths write the same report, and the stream keeps lines
/// back. The reference release is one before 116.3, under which every such class is
/// applied.
#[test]
fn both_paths_agree_on_the_committed_vcf_and_json_outputs() {
    let both_paths = BothPaths::new("records_committed");
    let whole_classes = tables().front_end_whole_classes();
    let mut classes: Vec<&str> = tables()
        .manifest_classes()
        .iter()
        .copied()
        .filter(|klass| !whole_classes.contains(klass))
        .collect();
    classes.extend(["regulatory_feature_rows", "order_only"]);
    let classes = classes.join(",");
    let input = corpus().join("variants.vcf");
    let manifest = corpus().join("manifest.json");
    let options = [
        "--input",
        input.to_str().expect("a path of text"),
        "--manifest",
        manifest.to_str().expect("a path of text"),
        "--masks",
        &classes,
        "--floors",
        "IMPACT=100,SYMBOL=100,Feature=100,SIFT=1,Existing_variation=1,am_class=1",
        "--reference-release",
        "116.2",
    ];
    for fmt in RECORD_FORMATS {
        let committed = expected_path(fmt);
        let run = both_paths.files(&committed, &committed, &options, fmt);
        assert!(run.code == 0 || run.code == 1, "{}", run.err);
        assert_eq!(
            number(run.report(), &["counts", "expected_entries"]),
            number(run.report(), &["counts", "actual_entries"])
        );
        assert!(kept_some_lines_back(&run.err), "{}", run.err);
        let text = read_text(&committed).expect("the committed output reads");
        let (old, new, impact) = match fmt {
            "vcf" => ("|MODIFIER|", "|LOW|", "IMPACT"),
            _ => ("\"impact\":\"MODIFIER\"", "\"impact\":\"LOW\"", "impact"),
        };
        assert!(text.contains(old));
        let changed = planted(
            &both_paths.directory,
            &format!("changed.{fmt}"),
            &text.replacen(old, new, 1),
        );
        let run = both_paths.files(&committed, &changed, &options, fmt);
        assert_eq!(run.code, 1);
        assert!(run
            .out
            .contains(&format!("  field {impact}: 1 differing entries\n")));
        assert!(
            run.err.contains(" reference lines were kept back by "),
            "{}",
            run.err
        );
        let run = both_paths.files(&committed, &changed, &options[..2], fmt);
        assert!(
            run.err.contains(" reference lines were kept back by "),
            "{}",
            run.err
        );
    }
}

/// A generator of the numbers a seeded test draws.
struct Seeded(u64);

impl Seeded {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut mixed = self.0;
        mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        mixed ^ (mixed >> 31)
    }

    fn below(&mut self, bound: usize) -> usize {
        (self.next() % bound as u64) as usize
    }

    /// A number from `low` to before `high`.
    fn between(&mut self, low: usize, high: usize) -> usize {
        low + self.below(high - low)
    }

    fn chance(&mut self, chance: f64) -> bool {
        ((self.next() >> 11) as f64 / (1u64 << 53) as f64) < chance
    }

    fn pick<T: Copy>(&mut self, choices: &[T]) -> T {
        choices[self.below(choices.len())]
    }

    fn shuffle<T>(&mut self, items: &mut [T]) {
        for last in (1..items.len()).rev() {
            items.swap(last, self.below(last + 1));
        }
    }
}

/// A random table of records: some carrying an `am_class`, a regulatory row, the misspelt
/// conservation key, a custom source or a number.
fn random_table(rng: &mut Seeded) -> Vec<Record> {
    let mut table = Vec::new();
    for r in 0..rng.between(3, 11) {
        let features: Vec<String> = (0..rng.between(1, 4))
            .map(|feature| format!("ENST{feature}"))
            .collect();
        let mut record = planted_record_of(r, &features);
        if rng.chance(0.3) {
            for row in &mut record.rows {
                *row = row.with(&[("am_class", "likely_benign")]);
            }
        }
        if rng.chance(0.2) {
            let allele = record.rows[0].field("Allele").to_string();
            record.rows.push(Row::of(&[
                ("Allele", &allele),
                ("Consequence", "regulatory_region_variant"),
                ("IMPACT", "MODIFIER"),
                ("Feature_type", "RegulatoryFeature"),
                ("Feature", &format!("ENSR{r}")),
            ]));
        }
        let kind = rng.below(100);
        if kind < 15 {
            record.info = vec!["phlyoP100way=0.5".to_string()];
            record = record.with_top("phlyoP100way", Value::Float(0.5));
        } else if kind < 30 {
            let source = json(&format!(r#"{{"mysource": [{{"name": "m{r}"}}]}}"#));
            record = record.with_top("custom_annotations", source);
        } else if kind < 50 {
            record = record.with_top("score", Value::Float(1.5));
        }
        table.push(record);
    }
    table
}

/// `reference` with up to four edits: a value changed, a row dropped or added, a record
/// repeated, exchanged with another or moved, a record's rows turned round, a list in
/// another order, a record the input lacks, a record dropped.
fn random_candidate(rng: &mut Seeded, reference: &[Record]) -> Vec<Record> {
    let mut candidate = reference.to_vec();
    for _ in 0..rng.pick(&[0, 0, 0, 1, 1, 2, 3, 4]) {
        let at = rng.below(candidate.len());
        let record = candidate[at].clone();
        match rng.below(10) {
            0 => candidate[at] = with_row(&record, 0, &[("IMPACT", "LOW")]),
            1 if record.rows.len() > 1 => {
                candidate[at] = record.with_rows(record.rows[1..].to_vec());
            }
            1 => {}
            2 => {
                let mut rows = record.rows.clone();
                rows.push(record.rows[0].with(&[("Feature", "ENST9")]));
                candidate[at] = record.with_rows(rows);
            }
            3 => candidate.insert(at, record),
            4 if at + 1 < candidate.len() => candidate.swap(at, at + 1),
            4 => {}
            5 => {
                let other = rng.below(candidate.len());
                candidate.swap(at, other);
            }
            6 => {
                let mut rows = record.rows.clone();
                rows.reverse();
                candidate[at] = record.with_rows(rows);
            }
            7 if !record.rows[0].field("NearestExonJB").is_empty() => {
                let mut items: Vec<&str> =
                    record.rows[0].field("NearestExonJB").split('&').collect();
                items.reverse();
                candidate[at] = with_row(&record, 0, &[("NearestExonJB", &items.join("&"))]);
            }
            7 => {}
            8 => candidate.insert(at, planted_record(40 + rng.below(3))),
            _ if candidate.len() > 1 => {
                candidate.remove(at);
            }
            _ => {}
        }
    }
    candidate
}

/// A JSON output with some objects' members in another order, another engine identity on
/// the candidate's side, a number spelt another way and blank lines.
fn respelt_json(rng: &mut Seeded, text: &str, candidate: bool) -> String {
    let mut lines = lines_of(text);
    for line in &mut lines {
        if rng.chance(0.2) {
            let Ok(members) = json(line).into_dict() else {
                panic!("a planted line is an object");
            };
            let mut members: Vec<(String, Value)> = members.into_iter().collect();
            rng.shuffle(&mut members);
            if candidate && rng.chance(0.5) {
                members.push(("vep_rs_version".to_string(), text_value("0.0.0")));
            }
            *line = dumped_object(members.into_iter().collect());
        }
        if rng.chance(0.1) {
            let respelt = rng.pick(&["1.50", "15e-1", "1.6"]);
            *line = line.replace("\"score\": 1.5", &format!("\"score\": {respelt}"));
        }
    }
    let blank = if rng.chance(0.3) { rng.below(3) } else { 0 };
    for _ in 0..blank {
        let at = rng.below(lines.len() + 1);
        lines.insert(at, rng.pick(&["", " ", "\t"]).to_string());
    }
    lines.iter().map(|line| format!("{line}\n")).collect()
}

fn text_value(value: &str) -> Value {
    text(value)
}

/// A manifest documenting two records of `variants`: a consequence set, a field value and
/// a row the candidate alone writes.
fn random_manifest(rng: &mut Seeded, variants: &[Record], klass: &str) -> String {
    let (first, second) = (rng.below(variants.len()), rng.below(variants.len()));
    let row = &variants[first].rows[0];
    let other = variants[second].rows.last().expect("a record has a row");
    let set = rng.pick(&[
        row.field("Consequence"),
        row.field("Consequence"),
        "splice_region_variant",
    ]);
    let value = rng.pick(&["LOW", "MODIFIER", "MODIFIER"]);
    format!(
        r#"{{"divergences": [{{"location": "-", "allele": "{}", "feature": "{}", "record_indices": [{first}],
             "expected_divergence": "{klass}", "vep_rs_consequence_sets": ["{set}"]}}],
            "field_divergences": [{{"location": "-", "allele": "{}", "feature": "{}", "record_indices": [{second}],
             "field": "IMPACT", "vep": "MODIFIER", "vep_rs": "{value}", "expected_divergence": "{klass}"}}],
            "vep_rs_only_tuples": [{{"location": "-", "allele": "{}", "feature": "ENST9", "record_indices": [{second}],
             "expected_divergence": "{klass}"}}]}}"#,
        row.field("Allele"),
        row.field("Feature"),
        other.field("Allele"),
        other.field("Feature"),
        other.field("Allele"),
    )
}

/// A hundred and fifty seeded pairs of outputs written from one table of records, some
/// carrying a regulatory row, a plugin's object, a custom source, the misspelt conservation
/// key or a number. The input holds a record twice or lacks one; the candidate's table has a
/// value changed, a row dropped or added, a record dropped, repeated, exchanged with another
/// or moved, a record's rows turned round, a list in another order and a record the input
/// lacks; the reference may hold such a record too. A JSON side has blank lines, its
/// objects' members in another order and a number spelt another way. Each pair is compared
/// under a random choice of a manifest documenting records, the regulatory, conservation and
/// order classes, `--unordered` and floors: both paths write the same report, standard
/// output and exit code.
#[test]
fn both_paths_agree_on_random_vcf_and_json_pairs() {
    let both_paths = BothPaths::new("records_random_pairs");
    let klass = first_manifest_class();
    for fmt in RECORD_FORMATS {
        let mut rng = Seeded(if fmt == "vcf" { 104_729 } else { 1_299_709 });
        let floors: &[&str] = match fmt {
            "vcf" => &[
                "IMPACT",
                "NearestExonJB",
                "am_class",
                "SYMBOL",
                "phlyoP100way",
            ],
            _ => &[
                "IMPACT",
                "Feature",
                "NearestExonJB",
                "am_class",
                "mysource_SCORE",
                "NOPE",
                "STRAND",
                "Existing_variation",
                "input",
                "phlyoP100way",
                "score",
            ],
        };
        let mut verdicts: BTreeMap<String, usize> = BTreeMap::new();
        let mut kept = 0;
        for pair in 0..150 {
            let table = random_table(&mut rng);
            let (mut reference, mut variants) = (table.clone(), table);
            if rng.chance(0.2) {
                let at = rng.below(reference.len());
                reference.insert(at, reference[at].clone());
                variants.insert(at, variants[at].clone());
            }
            if rng.chance(0.2) {
                variants.remove(rng.below(variants.len()));
            }
            if rng.chance(0.1) {
                let at = rng.below(reference.len() + 1);
                reference.insert(at, planted_record(40 + rng.below(3)));
            }
            let mut candidate = random_candidate(&mut rng, &reference);
            if rng.chance(0.08) {
                let at = rng.below(reference.len() + 1);
                reference.insert(at, planted_record(40 + rng.below(3)));
            }
            let mut masks: Vec<&str> = Vec::new();
            for (name, chance) in [
                ("order_only", 0.4),
                ("regulatory_feature_rows", 0.5),
                ("conservation_short_name_spelling", 0.4),
            ] {
                if rng.chance(chance) {
                    masks.push(name);
                }
            }
            // The candidate of a class's own shape: no regulatory row, the conservation
            // key spelt right.
            if masks.contains(&"regulatory_feature_rows") && rng.chance(0.8) {
                for record in &mut candidate {
                    record
                        .rows
                        .retain(|row| row.field("Feature_type") == "Transcript");
                }
            }
            if masks.contains(&"conservation_short_name_spelling") && rng.chance(0.8) {
                for record in candidate
                    .iter_mut()
                    .filter(|record| !record.info.is_empty())
                {
                    record.info = vec!["phyloP100way=0.5".to_string()];
                    *record = record.with_top("phyloP100way", Value::Float(0.5));
                }
            }
            let mut sides = [rendered(fmt, &reference), rendered(fmt, &candidate)];
            if fmt == "json" {
                sides = [
                    respelt_json(&mut rng, &sides[0], false),
                    respelt_json(&mut rng, &sides[1], true),
                ];
            }
            let mut options: Vec<String> = Vec::new();
            let mut input_name = "variants.vcf";
            if rng.chance(0.4) {
                let manifest = planted(
                    &both_paths.directory,
                    &format!("manifest.{fmt}.{pair}.json"),
                    &random_manifest(&mut rng, &variants, klass),
                );
                options.push("--manifest".to_string());
                options.push(manifest.to_str().expect("a path of text").to_string());
                masks.push(klass);
                if rng.chance(0.15) {
                    input_name = "other.vcf";
                }
            }
            if !masks.is_empty() {
                options.push("--masks".to_string());
                options.push(masks.join(","));
            }
            if rng.chance(0.3) {
                options.push("--unordered".to_string());
            }
            let mut chosen: Vec<String> = Vec::new();
            for floor in floors {
                if rng.chance(0.25) {
                    chosen.push(format!("{floor}={}", rng.pick(&[0, 0, 0, 1, 3])));
                }
            }
            if !chosen.is_empty() {
                options.push("--floors".to_string());
                options.push(chosen.join(","));
            }
            let options: Vec<&str> = options.iter().map(String::as_str).collect();
            let variants = planted_input(&variants);
            let run =
                both_paths.texts_beside(fmt, &sides[0], &sides[1], &variants, &options, input_name);
            assert!(run.code == 0 || run.code == 1, "{}", run.err);
            let verdict = run.out.lines().next().unwrap_or_default().to_string();
            *verdicts.entry(verdict).or_default() += 1;
            kept += usize::from(kept_some_lines_back(&run.err));
        }
        let count = |verdict: &str| verdicts.get(&format!("{fmt}: {verdict}")).copied();
        assert!(
            count("green") > Some(15) && count("red") > Some(60) && kept > 100,
            "{verdicts:?}, {kept} pairs with lines kept back"
        );
    }
}

/// A JSON string spells a superscript digit, or a number longer than the interpreter
/// converts, with `\u` escapes on a line of ASCII alone whose digit runs are short. A
/// term-pair class converts the value on the whole outputs and ends the comparison with the
/// comparison's error, so the line is handed over and the handed-over comparison ends with
/// the same error; a `\u` escape of a letter is read and compared.
#[test]
fn a_json_object_whose_escape_spells_a_value_the_interpreter_does_not_convert_is_handed_over() {
    let both_paths = BothPaths::new("records_escaped_value");
    let records = planted_records();
    let (reference, variants) = (planted_json(&records), planted_input(&records));
    let long = format!("\"{}\"", format!("{}\\u0037", "7".repeat(600)).repeat(8));
    let spelt = [
        ("mature_mirna_unshifted_span", "\"\\u00b2\"".to_string()),
        ("stop_codon_unshifted_span", long),
    ];
    let impact = "\"impact\": \"MODIFIER\"";
    for (klass, value) in spelt {
        let text = reference.replacen(impact, &format!("{impact}, \"shift_length\": {value}"), 1);
        let longest = text
            .split(|c: char| !c.is_ascii_digit())
            .map(str::len)
            .max();
        assert!(text.is_ascii() && longest <= Some(640));
        let run = both_paths.texts("json", &text, &text, &variants, &["--masks", klass]);
        assert!(run.report.is_none());
        assert_eq!(run.verdict(), ("", 2));
        assert_eq!(run.err, run.whole_error);
        assert!(
            run.err
                .starts_with("ERROR: [compare_vep_formats] SHIFT_LENGTH "),
            "{}",
            run.err
        );
    }
    let escaped = "\"impact\": \"MODIF\\u0049ER\", \"shift_length\": \"\\u0037\"";
    let text = reference.replacen(impact, escaped, 1);
    let options = ["--masks", "stop_codon_unshifted_span"];
    let run = both_paths.texts("json", &text, &text, &variants, &options);
    assert_eq!(run.verdict(), ("json: green\n", 0));
    assert_eq!(run.err, kept_back(7, 8, 1, 8));
}

/// A JSON object of forty consequence items opens more than 64 brackets and nests no deeper
/// than an item's members: its line is kept back. A line is handed over for its nesting
/// alone: one whose arrays stand 65 levels deep, the object the first level, is handed over
/// and compared, and one of 64 levels is kept back.
#[test]
fn a_json_object_of_many_items_is_kept_back_and_one_nested_deep_is_handed_over() {
    let both_paths = BothPaths::new("records_nesting");
    let features: Vec<String> = (0..40).map(|feature| format!("ENST{feature}")).collect();
    let mut records = planted_records();
    records[0] = planted_record_of(0, &features);
    let run = both_paths.records("json", &records, &records, &[]);
    assert_eq!(run.verdict(), ("json: green\n", 0));
    assert_eq!(number(run.report(), &["counts", "compared"]), 62);
    let first = planted_json(&records[..1]);
    assert!(first.matches(['[', '{']).count() > 64);
    assert_eq!(run.err, kept_back(8, 8, 0, 8));
    for (levels, kept) in [(64, 8), (65, 7)] {
        let mut nested = Value::List(Vec::new());
        for _ in 0..levels - 2 {
            nested = Value::List(vec![nested]);
        }
        let deep = with_record(&records, 0, records[0].with_top("deep", nested));
        let run = both_paths.records("json", &deep, &deep, &[]);
        assert_eq!(run.verdict(), ("json: green\n", 0));
        assert_eq!(run.err, kept_back(kept, 8, 8 - kept, 8));
    }
}

/// A VCF or JSON comparison without its input is compared on the whole outputs, where the
/// parser refuses it; both paths end with that error.
#[test]
fn a_vcf_or_json_comparison_without_its_input_is_left_on_the_whole_outputs() {
    let both_paths = BothPaths::new("records_without_input");
    for fmt in RECORD_FORMATS {
        let committed = expected_path(fmt);
        let run = both_paths.files(&committed, &committed, &[], fmt);
        assert!(run.report.is_none());
        assert_eq!(run.verdict(), ("", 2));
        assert_eq!(run.err, run.whole_error);
        assert_eq!(
            run.err,
            format!("ERROR: [compare_vep_formats] --format {fmt} needs --input, the input VCF\n")
        );
        let mut front_end = FrontEnd::new(STREAM);
        let documented = no_documented();
        let kept = hand_over(
            &mut front_end,
            fmt,
            &committed,
            &committed,
            &documented,
            false,
            None,
            &[],
            None,
        )
        .expect("nothing is read");
        assert!(kept.is_none());
        assert_eq!(
            front_end.note,
            format!("the whole outputs were compared: a {fmt} comparison without --input is outside the front end")
        );
    }
    assert_eq!(
        FrontEnd::new(STREAM).note,
        "the whole outputs were compared: only a default, tab, vcf, json or parquet comparison is handed over"
    );
}
