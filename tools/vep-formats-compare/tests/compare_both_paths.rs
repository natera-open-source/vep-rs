// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison through the stream of record groups against the comparison of the whole
//! outputs. Each test compares planted outputs twice and holds the report, the standard
//! output and the exit code of the two runs to the same bytes. Every row is written here.

mod compare_common;

use std::cell::Cell;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;

use compare_common::*;
use vep_formats_compare::cache::THREADS_VARIABLE;
use vep_formats_compare::classes::tables;
use vep_formats_compare::consts::MISSPELT_CONSERVATION_KEY;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::pyfmt;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::remainder::hand_over;
use vep_formats_compare::remainder::literal_head;
use vep_formats_compare::remainder::FrontEnd;
use vep_formats_compare::remainder::STREAM;
use vep_formats_compare::remainder::WHOLE_VARIABLE;

/// One transcript row of record `record`, at a Location of the record's own.
fn planted_row(record: usize, feature: &str) -> String {
    let extra = format!("IMPACT=MODIFIER;STRAND=1;DOMAINS=db:a{record},db:b{record}");
    planted_row_with(record, feature, &extra)
}

fn planted_row_with(record: usize, feature: &str, extra: &str) -> String {
    planted_row_at(record, feature, extra, &format!("21:{}", 1000 + record))
}

fn planted_row_at(record: usize, feature: &str, extra: &str, location: &str) -> String {
    let allele = &"ACGT"[record % 4..record % 4 + 1];
    default_row(
        &format!("rec{record}"),
        location,
        allele,
        feature,
        "intron_variant",
        extra,
    )
}

/// Two transcript rows for each of eight records.
fn planted_rows() -> Vec<String> {
    (0..8)
        .flat_map(|record| ["ENST1", "ENST2"].map(|feature| planted_row(record, feature)))
        .collect()
}

/// `rows` with the row at `at` replaced.
fn with_row(rows: &[String], at: usize, row: String) -> Vec<String> {
    let mut planted = rows.to_vec();
    planted[at] = row;
    planted
}

/// `rows` with `more` put before the row at `at`.
fn with_rows_before(rows: &[String], at: usize, more: &[String]) -> Vec<String> {
    let mut planted = rows[..at].to_vec();
    planted.extend(more.iter().cloned());
    planted.extend(rows[at..].iter().cloned());
    planted
}

/// What one comparison came to: the parsed report (`None` where none was written), the
/// standard output, the exit code, and what the stream did with the two outputs as one line,
/// or the standard error of a run that ended in an error.
struct Compared {
    report: Option<Value>,
    out: String,
    code: i32,
    err: String,
}

impl Compared {
    fn report(&self) -> &Value {
        self.report.as_ref().expect("the run wrote a report")
    }
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
            self.note(expected, actual, options, fmt)
        };
        Compared {
            report: handed
                .0
                .map(|bytes| json(std::str::from_utf8(&bytes).expect("a report is text"))),
            out: handed.1,
            code: handed.2,
            err,
        }
    }

    /// What the stream of record groups did with the two outputs under the options of one
    /// run, as one line: how many rows were kept back and how many groups handed over.
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
        let documented =
            Documented::from_manifest(manifest.as_ref(), &masks(&classes), None, None, None)
                .expect("the manifest reads");
        let mut front_end = FrontEnd::new(STREAM);
        hand_over(
            &mut front_end,
            fmt,
            expected,
            actual,
            &documented,
            options.contains(&"--unordered"),
            value_of("--keep-only"),
            &[],
            None,
        )
        .expect("the outputs are read");
        format!("[compare_vep_formats] {}\n", front_end.note)
    }

    fn rows(&self, reference: &[String], candidate: &[String], options: &[&str]) -> Compared {
        let run = self.runs.get() + 1;
        let expected = planted(
            &self.directory,
            &format!("expected.{run}.txt"),
            &default_text(reference),
        );
        let actual = planted(
            &self.directory,
            &format!("actual.{run}.txt"),
            &default_text(candidate),
        );
        self.files(&expected, &actual, options, "default")
    }
}

/// Sixteen rows written alike: every row is kept back, the report counts them all, and
/// the stream says what it handed over.
#[test]
fn the_groups_written_alike_are_kept_back() {
    let both_paths = BothPaths::new("kept_back_alike");
    let rows = planted_rows();
    let run = both_paths.rows(&rows, &rows, &[]);
    assert_eq!((run.out.as_str(), run.code), ("default: green\n", 0));
    for count in ["expected_entries", "actual_entries", "compared"] {
        assert_eq!(number(run.report(), &["counts", count]), 16);
    }
    for side in ["expected_populated", "actual_populated"] {
        assert_eq!(number(run.report(), &["columns", "DOMAINS", side]), 16);
    }
    assert_eq!(
        run.err,
        format!(
            "[compare_vep_formats] 16 of 16 reference rows were kept back by {STREAM}, 0 of 8 record groups handed over\n"
        )
    );
}

/// A DOMAINS list in another order is accepted under `order_only` and a difference
/// without it; an IMPACT of another value is a difference under any class. One record group
/// of the eight is handed over.
#[test]
fn a_differing_field_is_handed_over_whether_a_class_accepts_it_or_not() {
    let both_paths = BothPaths::new("differing_field");
    let rows = planted_rows();
    let reordered = with_row(
        &rows,
        6,
        planted_row_with(3, "ENST1", "IMPACT=MODIFIER;STRAND=1;DOMAINS=db:b3,db:a3"),
    );
    let run = both_paths.rows(&rows, &reordered, &["--masks", "order_only"]);
    assert_eq!((run.out.as_str(), run.code), ("default: green\n", 0));
    assert_eq!(
        at(run.report(), &["masks", "applied"]),
        &json(r#"{"order_only": 2}"#)
    );
    assert!(run.err.contains("14 of 16 reference rows were kept back"));
    assert!(run.err.contains("1 of 8 record groups handed over"));
    let run = both_paths.rows(&rows, &reordered, &[]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "default: red\n  field DOMAINS: 1 differing entries\n  field Extra: 1 differing entries\n"
    );
    let other = with_row(
        &rows,
        6,
        planted_row_with(3, "ENST1", "IMPACT=LOW;STRAND=1;DOMAINS=db:a3,db:b3"),
    );
    let run = both_paths.rows(&rows, &other, &["--masks", "order_only"]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "default: red\n  field Extra: 1 differing entries\n  field IMPACT: 1 differing entries\n"
    );
    assert_eq!(
        number(run.report(), &["columns", "IMPACT", "expected_populated"]),
        16
    );
}

/// A row the reference alone writes is missing and a row the candidate alone writes is
/// unexpected, on a record both write and on a record one side skips.
#[test]
fn a_row_one_side_alone_writes_is_handed_over() {
    let both_paths = BothPaths::new("one_sided_row");
    let rows = planted_rows();
    let mut more = rows.clone();
    more.push(planted_row(3, "ENST3"));
    more.push(planted_row(9, "ENST1"));
    let run = both_paths.rows(&more, &rows, &[]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "default: red\n  missing in the candidate: 2 entries\n"
    );
    assert_eq!(
        at(run.report(), &["one_sided"]),
        &json(r#"{"expected_only": 2, "actual_only": 0}"#)
    );
    assert_eq!(number(run.report(), &["counts", "expected_entries"]), 18);
    assert_eq!(number(run.report(), &["counts", "compared"]), 16);
    let run = both_paths.rows(&rows, &more, &[]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "default: red\n  unexpected in the candidate: 2 entries\n"
    );
    assert_eq!(
        at(run.report(), &["one_sided"]),
        &json(r#"{"expected_only": 0, "actual_only": 2}"#)
    );
}

/// Two records whose groups change places are an order difference at the first of their
/// rows, counted among every shared key, the kept-back rows before it included; beside a row
/// one side alone writes the place counts the shared keys alone; `--unordered` reads no
/// order.
#[test]
fn the_place_of_an_order_difference_counts_every_shared_key() {
    let both_paths = BothPaths::new("order_place");
    let rows = planted_rows();
    let mut swapped = rows[..4].to_vec();
    swapped.extend_from_slice(&rows[6..8]);
    swapped.extend_from_slice(&rows[4..6]);
    swapped.extend_from_slice(&rows[8..]);
    let run = both_paths.rows(&rows, &swapped, &[]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "default: red\n  order differs at shared position 4\n"
    );
    assert_eq!(
        at(run.report(), &["order"]),
        &json(
            r#"{"first_divergence": 4, "expected": "21:1002|G|ENST1|Transcript", "actual": "21:1003|T|ENST1|Transcript"}"#
        )
    );
    assert!(run.err.contains("12 of 16 reference rows were kept back"));
    let mut reference = vec![rows[0].clone()];
    reference.extend_from_slice(&rows[2..]);
    let mut candidate = swapped.clone();
    candidate.push(planted_row(9, "ENST1"));
    let run = both_paths.rows(&reference, &candidate, &[]);
    assert_eq!(
        run.out,
        "default: red\n  unexpected in the candidate: 2 entries\n  order differs at shared position 3\n"
    );
    assert_eq!(number(run.report(), &["order", "first_divergence"]), 3);
    let run = both_paths.rows(&rows, &swapped, &["--unordered"]);
    assert_eq!((run.out.as_str(), run.code), ("default: green\n", 0));
    assert!(run.err.contains("16 of 16 reference rows were kept back"));
}

/// A record's two rows the other way round are the same rows as a set and an order
/// difference at the first of them.
#[test]
fn a_group_whose_rows_come_in_another_order_is_handed_over() {
    let both_paths = BothPaths::new("rows_in_another_order");
    let rows = planted_rows();
    let mut candidate = rows.clone();
    candidate.swap(10, 11);
    let run = both_paths.rows(&rows, &candidate, &[]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "default: red\n  order differs at shared position 10\n"
    );
}

/// A row written twice is two entries, the second under its key's counter: kept back
/// where both sides repeat it, and one missing entry where the candidate writes it once.
#[test]
fn a_repeated_key_is_kept_with_its_counter() {
    let both_paths = BothPaths::new("repeated_key");
    let rows = planted_rows();
    let twice = with_rows_before(&rows, 0, &rows[..1]);
    let run = both_paths.rows(&twice, &twice, &[]);
    assert_eq!((run.out.as_str(), run.code), ("default: green\n", 0));
    assert_eq!(number(run.report(), &["counts", "compared"]), 17);
    assert!(run.err.contains("17 of 17 reference rows were kept back"));
    let run = both_paths.rows(&twice, &rows, &[]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "default: red\n  missing in the candidate: 1 entries\n"
    );
    assert_eq!(
        texts(at(run.report(), &["missing", "intron_variant", "examples"])),
        ["21:1000|A|ENST1|Transcript#2"]
    );
}

/// A key the manifest documents is checked for staleness where the two sides agree: the
/// group is handed over though its rows are written alike, and the stale divergence is
/// reported; a candidate that writes the documented set is set aside as documented.
#[test]
fn a_key_the_manifest_documents_is_handed_over_with_equal_values() {
    let both_paths = BothPaths::new("documented_key");
    let rows = planted_rows();
    let klass = tables().manifest_classes()[0];
    let manifest = planted(
        &both_paths.directory,
        "manifest.json",
        &format!(
            r#"{{"divergences": [{{"location": "21:1005", "allele": "C", "feature": "ENST2", "expected_divergence": "{klass}",
                "vep_rs_consequence_sets": ["intron_variant,splice_region_variant"]}}],
              "field_divergences": [{{"location": "21:1006", "allele": "G", "feature": "ENST1", "field": "STRAND", "vep": "1",
                "vep_rs": "-1", "expected_divergence": "{klass}"}}]}}"#
        ),
    );
    let manifest = manifest.to_str().expect("a path of text");
    let options = ["--manifest", manifest, "--masks", klass];
    let run = both_paths.rows(&rows, &rows, &options);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "default: red\n  documented divergences not observed: 2\n"
    );
    assert!(run.err.contains("12 of 16 reference rows were kept back"));
    assert!(run.err.contains("2 of 8 record groups handed over"));
    let mut documented = rows.clone();
    documented[11] = documented[11].replace(
        "\tintron_variant\t",
        "\tsplice_region_variant,intron_variant\t",
    );
    documented[12] = documented[12].replace("STRAND=1", "STRAND=-1");
    let run = both_paths.rows(&rows, &documented, &options);
    assert_eq!((run.out.as_str(), run.code), ("default: green\n", 0));
    assert_eq!(
        number(run.report(), &["counts", "documented_divergences_skipped"]),
        1
    );
    assert_eq!(number(run.report(), &["counts", "compared"]), 15);
}

/// Under `conservation_short_name_spelling` the reference's key is rewritten before the
/// comparison: a candidate spelling it right agrees, each rewrite counted, and a candidate
/// that writes the misspelt key too, the group the same text on both sides, differs in both
/// columns.
#[test]
fn a_group_holding_the_misspelt_conservation_key_is_handed_over() {
    let both_paths = BothPaths::new("misspelt_key");
    let rows = planted_rows();
    let misspelt = with_row(
        &rows,
        8,
        planted_row_with(4, "ENST1", "IMPACT=MODIFIER;phlyoP100way=0.5"),
    );
    let spelt: Vec<String> = misspelt
        .iter()
        .map(|row| row.replace("phlyoP100way", "phyloP100way"))
        .collect();
    let options = ["--masks", "conservation_short_name_spelling"];
    let run = both_paths.rows(&misspelt, &spelt, &options);
    assert_eq!((run.out.as_str(), run.code), ("default: green\n", 0));
    assert_eq!(
        at(run.report(), &["masks", "applied"]),
        &json(r#"{"conservation_short_name_spelling": 1}"#)
    );
    assert!(run.err.contains("14 of 16 reference rows were kept back"));
    let run = both_paths.rows(&misspelt, &misspelt, &options);
    assert_eq!(run.code, 1);
    assert!(run.out.contains("field phlyoP100way: 1 differing entries"));
    assert!(run.out.contains("field phyloP100way: 1 differing entries"));
    assert!(run.err.contains("14 of 16 reference rows were kept back"));
    assert_eq!(
        literal_head(&MISSPELT_CONSERVATION_KEY).as_deref(),
        Some("phlyoP")
    );
}

/// `conservation_short_name_spelling` rewrites the whole reference text: a tab column
/// line naming the misspelt key on both sides is renamed on the reference's side alone, and a
/// misspelt key inside a reference Feature moves the row to another key. Neither shape shows
/// in a row's own text, so every row is handed over.
#[test]
fn everything_is_handed_over_where_the_conservation_rewrite_reaches_a_key() {
    let both_paths = BothPaths::new("rewrite_reaches_a_key");
    let header =
        "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tphlyoP100way";
    let mut lines = vec![header.to_string()];
    lines.extend((0..6).map(|r| {
        format!(
            "rec{r}\t21:{}\t{}\tENST1\tTranscript\tintron_variant\t0.{r}",
            1000 + r,
            &"ACGT"[r % 4..r % 4 + 1]
        )
    }));
    let text = lines.join("\n") + "\n";
    let expected = planted(&both_paths.directory, "e.tab", &text);
    let actual = planted(&both_paths.directory, "a.tab", &text);
    let options = ["--masks", "conservation_short_name_spelling"];
    let whole = format!(
        "[compare_vep_formats] every row was handed over by {STREAM}: a listed text stands in the reference's column line or a key column\n"
    );
    let run = both_paths.files(&expected, &actual, &options, "tab");
    assert_eq!(run.code, 1);
    assert_eq!(
        run.out,
        "tab: red\n  field phlyoP100way: 6 differing entries\n  field phyloP100way: 6 differing entries\n"
    );
    assert_eq!(run.err, whole);
    let run = both_paths.files(&expected, &actual, &[], "tab");
    assert_eq!((run.out.as_str(), run.code), ("tab: green\n", 0));
    assert!(run.err.contains("6 of 6 reference rows were kept back"));
    let rows = planted_rows();
    let keyed = with_rows_before(
        &rows,
        4,
        &[default_row(
            "rec2",
            "21:1002",
            "G",
            "phlyoP7way",
            "intron_variant",
            "IMPACT=MODIFIER",
        )],
    );
    let spelt: Vec<String> = keyed
        .iter()
        .map(|row| row.replace("phlyoP7way", "phyloP7way"))
        .collect();
    let run = both_paths.rows(&keyed, &spelt, &options);
    assert_eq!(
        (run.out.as_str(), run.code, run.err.as_str()),
        ("default: green\n", 0, whole.as_str())
    );
    let run = both_paths.rows(&keyed, &keyed, &options);
    assert_eq!(run.code, 1);
    assert!(run.out.contains("missing in the candidate: 1 entries"));
    assert_eq!(run.err, whole);
}

/// The text handed to the stream must be in every reference line the class rewrites:
/// a pattern with an alternation, a flag, no literal head or a quantifier on its head's
/// last character gives none or a shorter one.
#[test]
fn the_literal_head_is_one_every_match_starts_with() {
    let head =
        |pattern: &str| literal_head(&regex::Regex::new(pattern).expect("the pattern compiles"));
    assert_eq!(head(r"abc(\d+)x").as_deref(), Some("abc"));
    assert_eq!(head(r"abc?d").as_deref(), Some("ab"));
    assert_eq!(head(r"ab{2}").as_deref(), Some("a"));
    assert_eq!(head(r"abc|x"), None);
    assert_eq!(head(r"(abc)"), None);
    assert_eq!(head(r"(?i)abc"), None);
    assert_eq!(head(r"a?"), None);
}

/// A Location that returns after two other records' rows is one group: written alike it
/// is kept back with the rows between; a difference in its late row hands over the group
/// and the rows between, and an order difference inside that span is found at its place
/// among all shared keys.
#[test]
fn a_location_whose_rows_stand_apart_is_one_group() {
    let both_paths = BothPaths::new("location_apart");
    let rows = planted_rows();
    let apart = with_rows_before(&rows, 10, &[planted_row(2, "ENST9")]);
    let run = both_paths.rows(&apart, &apart, &[]);
    assert_eq!((run.out.as_str(), run.code), ("default: green\n", 0));
    assert!(run.err.contains("17 of 17 reference rows were kept back"));
    let changed = planted_row_with(2, "ENST9", "IMPACT=LOW;STRAND=1;DOMAINS=db:a2,db:b2");
    let late = with_rows_before(&rows, 10, std::slice::from_ref(&changed));
    let run = both_paths.rows(&apart, &late, &[]);
    assert_eq!(run.code, 1);
    assert!(run.out.contains("field IMPACT: 1 differing entries"));
    assert!(run.err.contains("10 of 17 reference rows were kept back"));
    assert!(run.err.contains("3 of 8 record groups handed over"));
    let mut moved = rows[..6].to_vec();
    moved.extend([
        rows[8].clone(),
        rows[9].clone(),
        rows[6].clone(),
        rows[7].clone(),
        changed,
    ]);
    moved.extend_from_slice(&rows[10..]);
    let run = both_paths.rows(&apart, &moved, &[]);
    assert_eq!(run.code, 1);
    assert!(run.out.contains("order differs at shared position 6"));
    assert_eq!(number(run.report(), &["order", "first_divergence"]), 6);
}

/// The two engines' `##` lines differ and no entry reads them: the rows are kept back.
/// Under a column line that names a column differently no row's text settles its entry, and
/// every row is handed over.
#[test]
fn note_lines_that_differ_keep_the_rows_back_and_another_column_line_hands_them_over() {
    let both_paths = BothPaths::new("note_lines");
    let rows = planted_rows();
    let text = default_text(&rows);
    let expected = planted(
        &both_paths.directory,
        "e.txt",
        &format!("## reference engine\n## IMPACT : the reference's text\n{text}"),
    );
    let actual = planted(
        &both_paths.directory,
        "a.txt",
        &format!("## candidate engine\n{text}"),
    );
    let run = both_paths.files(&expected, &actual, &[], "default");
    assert_eq!((run.out.as_str(), run.code), ("default: green\n", 0));
    assert!(run.err.contains("16 of 16 reference rows were kept back"));
    let actual = planted(
        &both_paths.directory,
        "a.txt",
        &text.replacen("\tGene\t", "\tGENE\t", 1),
    );
    let run = both_paths.files(&expected, &actual, &[], "default");
    assert_eq!(run.code, 1);
    assert!(run.out.contains("field GENE: 16 differing entries"));
    assert!(run.out.contains("field Gene: 16 differing entries"));
    assert_eq!(
        run.err,
        format!(
            "[compare_vep_formats] every row was handed over by {STREAM}: the two outputs' column lines differ\n"
        )
    );
}

fn gzipped(text: &str) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder
        .write_all(text.as_bytes())
        .expect("the text is deflated");
    encoder.finish().expect("the member is closed")
}

/// A tab comparison of gzipped outputs, one padded with zero bytes: a floor is read against
/// the populated counts of every row, the kept-back ones included, on both sides.
#[test]
fn the_kept_back_rows_count_under_a_floor_in_the_tab_format() {
    let both_paths = BothPaths::new("floor_over_kept_rows");
    let header =
        "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tIMPACT\tSYMBOL";
    let rows: Vec<String> = (0..8)
        .flat_map(|r| {
            [1, 2].map(|t| {
                format!(
                    "rec{r}\t21:{}\t{}\tENST{t}\tTranscript\tintron_variant\tMODIFIER\t{}",
                    1000 + r,
                    &"ACGT"[r % 4..r % 4 + 1],
                    if r % 2 == 1 { "-" } else { "GENE" }
                )
            })
        })
        .collect();
    let candidate = with_row(&rows, 5, rows[5].replace("MODIFIER", "LOW"));
    let reference_text = format!("## reference\n{header}\n{}\n", rows.join("\n"));
    let candidate_text = format!("{header}\n{}\n", candidate.join("\n"));
    let expected = both_paths.directory.join("e.tab.gz");
    let actual = both_paths.directory.join("a.tab.gz");
    std::fs::write(&expected, gzipped(&reference_text)).expect("the file is written");
    // Zero bytes after a gzip member are padding to both readers.
    let mut padded = gzipped(&candidate_text);
    padded.extend_from_slice(&[0; 512]);
    std::fs::write(&actual, padded).expect("the file is written");
    let run = both_paths.files(
        &expected,
        &actual,
        &["--floors", "SYMBOL=8,IMPACT=16"],
        "tab",
    );
    assert_eq!(run.code, 1);
    assert_eq!(run.out, "tab: red\n  field IMPACT: 1 differing entries\n");
    for side in ["expected_populated", "actual_populated"] {
        assert_eq!(number(run.report(), &["columns", "SYMBOL", side]), 8);
    }
    assert_eq!(
        at(run.report(), &["columns", "SYMBOL", "floor_met"]),
        &Value::Bool(true)
    );
    assert!(run.err.contains("14 of 16 reference rows were kept back"));
    let run = both_paths.files(&expected, &actual, &["--floors", "SYMBOL=9"], "tab");
    assert!(run
        .out
        .contains("  floor SYMBOL=9: expected populated 8, actual populated 8\n"));
}

/// The text of a JSON object nested `depth` deep.
fn nested_object(depth: usize) -> String {
    format!("{}1{}", r#"{"a":"#.repeat(depth), "}".repeat(depth))
}

/// A GA4GH_VRS object is put in canonical form on every row. One nested 1,000 deep ends
/// the comparison of the whole outputs with the comparator's error; the stream hands a
/// row over that holds more than 64 opening brackets, so the handed-over comparison ends
/// with the same error. 60,000 opening brackets, which are no object, end both runs alike:
/// with that error or compared as text. An object nested 64 deep is compared and its row
/// kept back. Beside a candidate that cannot be read, the nested object ends both runs with
/// exit 2 and no report.
#[test]
fn a_row_whose_object_nests_too_deep_to_read_is_handed_over() {
    let both_paths = BothPaths::new("nested_object");
    let rows = planted_rows();
    let with_first = |extra: &str| with_row(&rows, 0, planted_row_with(0, "ENST1", extra));
    let nested = with_first(&format!("GA4GH_VRS={}", nested_object(1000)));
    let run = both_paths.rows(&nested, &nested, &[]);
    assert!(run.report.is_none());
    assert_eq!((run.out.as_str(), run.code), ("", 2));
    assert_eq!(
        run.err,
        "ERROR: [compare_vep_formats] a JSON object is nested too deep to read: {\"a\":{\"a\":{\"a\":{\"a\":{\"a\"...\n"
    );
    // Beside a candidate that cannot be read, both runs end with the comparator's error.
    let reference = planted(&both_paths.directory, "nested.txt", &default_text(&nested));
    let run = both_paths.files(
        &reference,
        &both_paths.directory.join("absent.txt"),
        &[],
        "default",
    );
    assert!(run.report.is_none());
    assert_eq!((run.out.as_str(), run.code), ("", 2));
    let opened = with_first(&format!("GA4GH_VRS={}", "[".repeat(60000)));
    let run = both_paths.rows(&opened, &opened, &[]);
    assert!(
        (run.report.is_none() && run.out.is_empty() && run.code == 2)
            || (run.out == "default: green\n" && run.code == 0)
    );
    let shallow = with_first(&format!("GA4GH_VRS={}", nested_object(64)));
    let run = both_paths.rows(&shallow, &shallow, &[]);
    assert_eq!((run.out.as_str(), run.code), ("default: green\n", 0));
    assert!(run.err.contains("16 of 16 reference rows were kept back"));
    let deeper = with_first(&format!("GA4GH_VRS={}", nested_object(65)));
    let run = both_paths.rows(&deeper, &deeper, &[]);
    assert_eq!((run.out.as_str(), run.code), ("default: green\n", 0));
    assert!(run.err.contains("14 of 16 reference rows were kept back"));
}

/// The deflate stream zlib writes for `text` at its highest level, without a wrapper.
fn deflated(text: &[u8]) -> Vec<u8> {
    let mut encoder = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::best());
    encoder.write_all(text).expect("the text is deflated");
    encoder.finish().expect("the stream is closed")
}

/// One gzip member holding `text` under a header of `flags`, `fields` standing between the
/// fixed header and the deflate stream `stream`.
fn gzip_member(text: &[u8], flags: u8, fields: &[u8], stream: &[u8]) -> Vec<u8> {
    let mut checksum = flate2::Crc::new();
    checksum.update(text);
    let mut member = vec![0x1f, 0x8b, 0x08, flags, 0, 0, 0, 0, 0, 0xff];
    member.extend_from_slice(fields);
    member.extend_from_slice(stream);
    member.extend_from_slice(&checksum.sum().to_le_bytes());
    member.extend_from_slice(
        &u32::try_from(text.len())
            .expect("a short text")
            .to_le_bytes(),
    );
    member
}

/// One final deflate block of the fixed codes (RFC 1951, 3.2.6) that writes `before`,
/// copies three bytes from five bytes back and writes `after`.
fn deflate_with_a_copy(before: &[u8], after: &[u8]) -> Vec<u8> {
    let mut bits: Vec<u8> = vec![1, 1, 0];
    fn code(bits: &mut Vec<u8>, value: u32, width: u32) {
        bits.extend((0..width).map(|i| u8::from((value >> (width - 1 - i)) & 1 == 1)));
    }
    fn literals(bits: &mut Vec<u8>, text: &[u8]) {
        for byte in text {
            assert!(*byte < 144, "a literal of the eight-bit codes");
            code(bits, 0x30 + u32::from(*byte), 8);
        }
    }
    literals(&mut bits, before);
    code(&mut bits, 1, 7);
    code(&mut bits, 4, 5);
    bits.push(0);
    literals(&mut bits, after);
    code(&mut bits, 0, 7);
    bits.chunks(8)
        .map(|chunk| {
            chunk
                .iter()
                .enumerate()
                .map(|(k, bit)| bit << k)
                .sum::<u8>()
        })
        .collect()
}

/// A gzip reader reads past a member's header checksum without a check and past the flag
/// bits the format reserves, and takes an extra field, a file name and a comment: outputs
/// gzipped under such headers are compared through the stream as on the whole outputs.
#[test]
fn a_gzip_header_is_read_alike_on_both_paths() {
    let both_paths = BothPaths::new("gzip_header");
    let rows = planted_rows();
    let other = with_row(
        &rows,
        6,
        planted_row_with(3, "ENST1", "IMPACT=LOW;STRAND=1;DOMAINS=db:a3,db:b3"),
    );
    let reference = both_paths.directory.join("expected.txt.gz");
    let candidate = both_paths.directory.join("actual.txt.gz");
    let reference_text = default_text(&rows);
    let candidate_text = default_text(&other);
    std::fs::write(
        &reference,
        gzip_member(
            reference_text.as_bytes(),
            0xE2,
            b"\xff\xff",
            &deflated(reference_text.as_bytes()),
        ),
    )
    .expect("the file is written");
    let headers: [(u8, &[u8]); 3] = [
        (0x02, b"\x12\x34"),
        (0x20, b""),
        (0x1E, b"\x02\x00xyname\0note\0\0\0"),
    ];
    for (flags, fields) in headers {
        std::fs::write(
            &candidate,
            gzip_member(
                candidate_text.as_bytes(),
                flags,
                fields,
                &deflated(candidate_text.as_bytes()),
            ),
        )
        .expect("the file is written");
        let run = both_paths.files(&reference, &candidate, &[], "default");
        assert_eq!(run.code, 1);
        assert_eq!(
            run.out,
            "default: red\n  field Extra: 1 differing entries\n  field IMPACT: 1 differing entries\n"
        );
        assert!(run.err.contains("14 of 16 reference rows were kept back"));
    }
}

/// A gzip member whose deflate stream copies from before its start is refused by zlib,
/// and the comparison of the whole outputs ends in the comparator's error. The stream
/// refuses it as well, under the checksum and length of the text a decoder that takes the
/// copy writes (a zero byte for each place before the start): both runs exit 2 and write
/// no report.
#[test]
fn a_deflate_stream_that_zlib_refuses_is_refused_on_both_paths() {
    let both_paths = BothPaths::new("refused_deflate");
    let rows = planted_rows();
    let text = default_text(&rows).into_bytes();
    let mut taken = text[..3].to_vec();
    taken.extend_from_slice(b"\0\0");
    taken.extend_from_slice(&text[..1]);
    taken.extend_from_slice(&text[3..]);
    let reference = both_paths.directory.join("expected.txt");
    let candidate = both_paths.directory.join("actual.txt.gz");
    std::fs::write(&reference, &text).expect("the file is written");
    std::fs::write(
        &candidate,
        gzip_member(&taken, 0, b"", &deflate_with_a_copy(&text[..3], &text[3..])),
    )
    .expect("the file is written");
    let run = both_paths.files(&reference, &candidate, &[], "default");
    assert!(run.report.is_none());
    assert_eq!((run.out.as_str(), run.code), ("", 2));
    assert_eq!(run.err.matches('\n').count(), 1);
    assert!(run
        .err
        .starts_with("ERROR: [compare_vep_formats] cannot read "));
    // The same stream with the copy inside the text is read by both.
    let mut copied = text[..5].to_vec();
    copied.extend_from_slice(&text[..3]);
    copied.extend_from_slice(&text[5..]);
    std::fs::write(
        &candidate,
        gzip_member(
            &copied,
            0,
            b"",
            &deflate_with_a_copy(&text[..5], &text[5..]),
        ),
    )
    .expect("the file is written");
    let run = both_paths.files(&reference, &candidate, &[], "default");
    assert_eq!(run.code, 1);
    assert!(run.out.starts_with("default: red\n"));
}

/// The generator of `random.Random` of the interpreter, seeded with a whole number, so
/// that the pairs of the random test are the same on every run of a seed.
struct MersenneTwister {
    state: [u32; 624],
    index: usize,
}

impl MersenneTwister {
    fn seeded(seed: u32) -> Self {
        let mut state = [0u32; 624];
        state[0] = 19650218;
        for i in 1..624 {
            state[i] = 1812433253u32
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(u32::try_from(i).expect("a small index"));
        }
        let key = [seed];
        let (mut i, mut j) = (1, 0);
        for _ in 0..624 {
            state[i] = (state[i] ^ (state[i - 1] ^ (state[i - 1] >> 30)).wrapping_mul(1664525))
                .wrapping_add(key[j])
                .wrapping_add(u32::try_from(j).expect("a small index"));
            i += 1;
            j += 1;
            if i >= 624 {
                state[0] = state[623];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..623 {
            state[i] = (state[i] ^ (state[i - 1] ^ (state[i - 1] >> 30)).wrapping_mul(1566083941))
                .wrapping_sub(u32::try_from(i).expect("a small index"));
            i += 1;
            if i >= 624 {
                state[0] = state[623];
                i = 1;
            }
        }
        state[0] = 0x8000_0000;
        MersenneTwister { state, index: 624 }
    }

    fn word(&mut self) -> u32 {
        if self.index >= 624 {
            for k in 0..624 {
                let y = (self.state[k] & 0x8000_0000) | (self.state[(k + 1) % 624] & 0x7fff_ffff);
                let mut next = self.state[(k + 397) % 624] ^ (y >> 1);
                if y & 1 == 1 {
                    next ^= 0x9908_b0df;
                }
                self.state[k] = next;
            }
            self.index = 0;
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// `random()`: a float in [0, 1) of 53 random bits.
    fn random(&mut self) -> f64 {
        let (a, b) = (self.word() >> 5, self.word() >> 6);
        (f64::from(a) * 67108864.0 + f64::from(b)) / 9007199254740992.0
    }

    /// `randrange(n)` for 0 < n < 2**32.
    fn below(&mut self, n: usize) -> usize {
        let n = u32::try_from(n).expect("a small range");
        let bits = 32 - n.leading_zeros();
        loop {
            let drawn = self.word() >> (32 - bits);
            if drawn < n {
                return drawn as usize;
            }
        }
    }

    /// `randrange(start, stop)`.
    fn between(&mut self, start: usize, stop: usize) -> usize {
        start + self.below(stop - start)
    }
}

#[test]
fn the_generator_draws_what_the_interpreters_generator_draws() {
    let mut rng = MersenneTwister::seeded(7919);
    let drawn = [
        rng.between(3, 14),
        rng.below(8),
        rng.below(20),
        rng.below(1),
        rng.between(1, 4),
    ];
    assert_eq!(drawn, FIRST_DRAWS);
    assert_eq!(format!("{:.17}", rng.random()), FIRST_FLOAT);
}

/// What `random.Random(7919)` draws first: `randrange(3, 14)`, `randrange(8)`,
/// `randrange(20)`, `randrange(1)`, `randrange(1, 4)`, then `random()` to 17 places.
const FIRST_DRAWS: [usize; 5] = [13, 7, 5, 0, 2];
const FIRST_FLOAT: &str = "0.91928733094594528";

/// Three hundred seeded pairs of outputs, each the reference's rows with rows changed,
/// dropped, added, repeated, exchanged within a record and across records, a Location
/// returning later and a list in another order, compared with and without `order_only`
/// and `--unordered`: both paths write the same report, standard output and exit code.
#[test]
fn both_paths_agree_on_random_pairs() {
    let both_paths = BothPaths::new("random_pairs");
    let domains = regex::Regex::new(r"DOMAINS=(db:a\d+),(db:b\d+)").expect("the pattern compiles");
    let default_extra =
        |record: usize| format!("IMPACT=MODIFIER;STRAND=1;DOMAINS=db:a{record},db:b{record}");
    let mut rng = MersenneTwister::seeded(7919);
    let mut verdicts: BTreeMap<String, usize> = BTreeMap::new();
    for _ in 0..300 {
        let records = rng.between(3, 14);
        let mut reference: Vec<String> = Vec::new();
        for r in 0..records {
            let place = if rng.random() < 0.15 {
                rng.below(r + 1)
            } else {
                r
            };
            let location = format!("21:{}", 1000 + place);
            for t in 0..rng.between(1, 4) {
                reference.push(planted_row_at(
                    r,
                    &format!("ENST{t}"),
                    &default_extra(r),
                    &location,
                ));
            }
        }
        let mut candidate = reference.clone();
        for _ in 0..rng.between(0, 5) {
            let at = rng.below(candidate.len());
            let edit = rng.below(8);
            if edit == 0 {
                candidate[at] = candidate[at].replace("IMPACT=MODIFIER", "IMPACT=LOW");
            } else if edit == 1 && candidate.len() > 1 {
                candidate.remove(at);
            } else if edit == 2 {
                let record = rng.below(20);
                let feature = format!("ENST{}", rng.below(5));
                candidate.insert(at, planted_row(record, &feature));
            } else if edit == 3 {
                candidate.insert(at, candidate[at].clone());
            } else if edit == 4 && at + 1 < candidate.len() {
                candidate.swap(at, at + 1);
            } else if edit == 5 {
                let other = rng.below(candidate.len());
                candidate.swap(at, other);
            } else if edit == 6 {
                candidate[at] = domains
                    .replace_all(&candidate[at], "DOMAINS=${2},${1}")
                    .into_owned();
            } else {
                let end = (at + rng.between(1, 4)).min(candidate.len());
                let block: Vec<String> = candidate.drain(at..end).collect();
                let to = rng.below(candidate.len() + 1);
                candidate.splice(to..to, block);
            }
        }
        let mut options: Vec<&str> = Vec::new();
        if rng.random() < 0.5 {
            options.extend(["--masks", "order_only"]);
        }
        if rng.random() < 0.3 {
            options.push("--unordered");
        }
        let run = both_paths.rows(&reference, &candidate, &options);
        let first_line = run.out.split('\n').next().unwrap_or_default().to_string();
        *verdicts.entry(first_line).or_default() += 1;
    }
    let seen = |verdict: &str| verdicts.get(verdict).copied().unwrap_or_default();
    assert!(
        seen("default: green") > 30 && seen("default: red") > 150,
        "{verdicts:?}"
    );
}
