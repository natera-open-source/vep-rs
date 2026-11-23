// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `allele_number_repeated_or_reference_equal_alt` through the comparison of two
//! outputs in the default, tab, VCF and JSON formats: the allele number of a row whose
//! record writes an alternate twice, or an alternate equal to the reference allele, before
//! the row's allele, on planted rows and on the stored reference of the committed corpus
//! that holds such records.

mod compare_common;

use std::collections::BTreeMap;

use compare_common::*;
use vep_formats_compare::classes::{self, ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT};
use vep_formats_compare::entry::Entry;
use vep_formats_compare::input::input_records;
use vep_formats_compare::input::read_text;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::{parse_default, parse_tab};
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::allele_number::ALLELE_NUMBER_FIELDS;

const CLASS: &str = ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT;
const FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];
const FEATURE: &str = "ENST00000416044";

/// A made-up record: its input line, and the name and the Location its rows carry.
#[derive(Clone, Copy)]
struct Record {
    line: &'static str,
    name: &'static str,
    location: &'static str,
}

/// The record of the committed corpus `GRCh38-vcf-edge-cases` that writes one alternate
/// twice before another.
const REPEATED: Record = Record {
    line: "21\t30289183\tdup_alt\tC\tG,G,T\t.\t.\t.",
    name: "dup_alt",
    location: "21:30289183",
};

/// A record whose first alternate is its reference allele.
const REFERENCE_EQUAL: Record = Record {
    line: "21\t30289190\tref_alt\tA\tA,C\t.\t.\t.",
    name: "ref_alt",
    location: "21:30289190",
};

/// A record that writes a `*` before an alternate it writes twice.
const STAR: Record = Record {
    line: "21\t30289183\tstar_alt\tC\t*,G,G,T\t.\t.\t.",
    name: "star_alt",
    location: "21:30289183",
};

/// A record that writes each alternate once.
const PLAIN: Record = Record {
    line: "21\t30289183\tmulti_snv\tC\tG,T\t.\t.\t.",
    name: "multi_snv",
    location: "21:30289183",
};

/// An indel record that writes a deletion twice before an insertion: its alleles share
/// their first base, which Ensembl VEP's VCF parser chops, so its rows carry `-` and `AT`.
const INDEL: Record = Record {
    line: "21\t30289183\tindel_alt\tCA\tC,C,CAT\t.\t.\t.",
    name: "indel_alt",
    location: "21:30289184",
};

/// An indel record whose last alternate begins with another base: nothing is chopped, and
/// its rows carry the alleles as written.
const UNANCHORED: Record = Record {
    line: "21\t30289183\tindel_apart\tCA\tC,C,GAT\t.\t.\t.",
    name: "indel_apart",
    location: "21:30289183",
};

/// One row of a record: its allele, its consequence term, its allele number and the
/// record's written allele list as the row carries it, `None` for a row that carries none.
#[derive(Clone, Copy)]
struct Row {
    allele: &'static str,
    term: &'static str,
    number: &'static str,
    list: Option<&'static str>,
}

/// A row of `allele` with one term an allele, which both engines write.
fn row(allele: &'static str, number: &'static str, list: Option<&'static str>) -> Row {
    let term = if allele == "G" {
        "stop_lost"
    } else {
        "stop_retained_variant"
    };
    Row {
        allele,
        term,
        number,
        list,
    }
}

/// The default format: the number and the list are items of the Extra column, and a row
/// without the list carries no such item.
fn default_output(record: &Record, rows: &[Row]) -> String {
    let lines: Vec<String> = rows
        .iter()
        .map(|row| {
            let list = row
                .list
                .map_or(String::new(), |list| format!(";UPLOADED_ALLELE={list}"));
            let extra = format!("ALLELE_NUM={}{list};IMPACT=LOW;STRAND=-1", row.number);
            default_row(
                record.name,
                record.location,
                row.allele,
                FEATURE,
                row.term,
                &extra,
            )
        })
        .collect();
    default_text(&lines)
}

/// The tab format: the number and the list are columns, and a row without the list writes
/// `-` in its column.
fn tab_output(record: &Record, rows: &[Row]) -> String {
    let mut lines = vec![
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tALLELE_NUM\tUPLOADED_ALLELE\tIMPACT"
            .to_string(),
    ];
    for row in rows {
        lines.push(format!(
            "{}\t{}\t{}\tENSG1\t{FEATURE}\tTranscript\t{}\t{}\t{}\tLOW",
            record.name,
            record.location,
            row.allele,
            row.term,
            row.number,
            row.list.unwrap_or("-")
        ));
    }
    lines.join("\n") + "\n"
}

/// The VCF format: the number and the list are fields of a CSQ entry, and a row without
/// the list writes an empty field.
fn vcf_output(record: &Record, rows: &[Row]) -> String {
    let header = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|Feature_type|Feature|ALLELE_NUM|UPLOADED_ALLELE\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let entries: Vec<String> = rows
        .iter()
        .map(|row| {
            format!(
                "{}|{}|LOW|Transcript|{FEATURE}|{}|{}",
                row.allele,
                row.term,
                row.number,
                row.list.unwrap_or("")
            )
        })
        .collect();
    let fixed = record
        .line
        .strip_suffix('.')
        .expect("an input line whose INFO is a point");
    format!("{header}{fixed}CSQ={}\n", entries.join(","))
}

/// The JSON format: the number is a number under `allele_num` and the list a text under
/// `uploaded_allele`, and an object without the list has no such member.
fn json_output(record: &Record, rows: &[Row]) -> String {
    let objects: Vec<String> = rows
        .iter()
        .map(|row| {
            let list = row.list.map_or(String::new(), |list| {
                format!(r#", "uploaded_allele": "{list}""#)
            });
            format!(
                r#"{{"variant_allele": "{}", "transcript_id": "{FEATURE}", "consequence_terms": ["{}"], "impact": "LOW", "allele_num": {}{list}}}"#,
                row.allele, row.term, row.number
            )
        })
        .collect();
    let input = record.line.replace('\t', "\\t");
    dumped(&format!(
        r#"{{"input": "{input}", "id": "{}", "most_severe_consequence": "stop_lost", "transcript_consequences": [{}]}}"#,
        record.name,
        objects.join(", ")
    )) + "\n"
}

/// The entries of one engine's output of `record` in `format`.
fn entries(format: &str, record: &Record, rows: &[Row]) -> Vec<Entry> {
    let inputs = [record.line.to_string()];
    match format {
        "default" => parse_default(&default_output(record, rows)),
        "tab" => parse_tab(&tab_output(record, rows)),
        "vcf" => parse_vcf(&vcf_output(record, rows), &inputs),
        _ => parse_json(&json_output(record, rows), &inputs),
    }
    .expect("the output parses")
}

/// The report's verdict and what was applied of two entry lists under `classes`.
fn under(
    classes: &[&str],
    expected: &[Entry],
    actual: &[Entry],
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let report = compare(expected, actual, &mut documented);
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// The verdict and what was applied of the two engines' rows of `record` in `format`.
fn compared(
    classes: &[&str],
    format: &str,
    record: &Record,
    reference: &[Row],
    candidate: &[Row],
) -> (String, BTreeMap<String, i64>) {
    under(
        classes,
        &entries(format, record, reference),
        &entries(format, record, candidate),
    )
}

fn taken() -> (String, BTreeMap<String, i64>) {
    ("green".to_string(), count_map(&[(CLASS, 1)]))
}

fn left() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), BTreeMap::new())
}

/// The row of the alternate written after the repeated one is taken under the class, which
/// is credited once for the entry, and is a difference without it and under another class.
/// The first alternate's row carries one number on both sides.
#[test]
fn the_row_after_a_repeated_alternate_is_taken_under_the_class_in_every_format() {
    let list = Some("C/G/G/T");
    let reference = [row("G", "1", list), row("T", "2", list)];
    let candidate = [row("G", "1", list), row("T", "3", list)];
    for format in FORMATS {
        let outcome =
            |classes: &[&str]| compared(classes, format, &REPEATED, &reference, &candidate);
        assert_eq!(outcome(&[CLASS]), taken(), "{format}");
        assert_eq!(outcome(&[]), left(), "{format}");
        assert_eq!(outcome(&["order_only"]).0, "red", "{format}");
        assert_eq!(
            outcome(&[classes::TOTAL_LENGTH_JSON_RANGE_END]).0,
            "red",
            "{format}"
        );
    }
}

/// An alternate equal to the reference allele is dropped as a repeated one is, and a `*`
/// takes a number as any alternate does.
#[test]
fn the_row_after_an_alternate_equal_to_the_reference_allele_or_after_a_star_is_taken() {
    let equal = Some("A/A/C");
    let star = Some("C/*/G/G/T");
    for format in FORMATS {
        assert_eq!(
            compared(
                &[CLASS],
                format,
                &REFERENCE_EQUAL,
                &[row("C", "1", equal)],
                &[row("C", "2", equal)]
            ),
            taken(),
            "{format}"
        );
        let reference = [row("G", "2", star), row("T", "3", star)];
        let candidate = [row("G", "2", star), row("T", "4", star)];
        assert_eq!(
            compared(&[CLASS], format, &STAR, &reference, &candidate),
            taken(),
            "{format}"
        );
        assert_eq!(
            compared(&[], format, &STAR, &reference, &candidate),
            left(),
            "{format}"
        );
    }
}

/// A row of an indel record carries its allele without the first base the record's alleles
/// share, and is taken by that spelling of a written alternate; where the alleles share no
/// first base the row carries the written allele. A row whose allele is neither spelling
/// of an alternate stays a difference.
#[test]
fn a_row_of_an_indel_record_is_taken_by_the_allele_the_parser_spells() {
    let chopped = Some("CA/C/C/CAT");
    let written = Some("CA/C/C/GAT");
    for format in FORMATS {
        // The repeated deletion's row is numbered 1 on both sides.
        let reference = [row("-", "1", chopped), row("AT", "2", chopped)];
        let candidate = [row("-", "1", chopped), row("AT", "3", chopped)];
        assert_eq!(
            compared(&[CLASS], format, &INDEL, &reference, &candidate),
            taken(),
            "{format}"
        );
        assert_eq!(
            compared(&[], format, &INDEL, &reference, &candidate),
            left(),
            "{format}"
        );
        assert_eq!(
            compared(
                &[CLASS],
                format,
                &UNANCHORED,
                &[row("GAT", "2", written)],
                &[row("GAT", "3", written)]
            ),
            taken(),
            "{format}"
        );
        for (record, list, allele) in [
            (&INDEL, chopped, "T"),
            (&INDEL, chopped, "A"),
            (&UNANCHORED, written, "AT"),
        ] {
            assert_eq!(
                compared(
                    &[CLASS],
                    format,
                    record,
                    &[row(allele, "2", list)],
                    &[row(allele, "3", list)]
                ),
                left(),
                "{format} {list:?} {allele}"
            );
        }
    }
}

/// The class reads the list on the row itself: a pair of which one row or both carry none
/// stays a difference.
#[test]
fn a_row_without_the_allele_list_stays_a_difference_under_the_class() {
    let list = Some("C/G/G/T");
    for format in FORMATS {
        for (reference, candidate) in [(None, None), (list, None), (None, list)] {
            assert_eq!(
                compared(
                    &[CLASS],
                    format,
                    &REPEATED,
                    &[row("T", "2", reference)],
                    &[row("T", "3", candidate)]
                ),
                left(),
                "{format} {reference:?} {candidate:?}"
            );
        }
    }
}

/// Numbers that are not the place and the count the row's list gives stay a difference.
#[test]
fn numbers_that_do_not_fit_the_list_stay_a_difference_under_the_class() {
    let cases = [
        // The candidate numbered 3 where the allele's place is 2.
        (&REFERENCE_EQUAL, "A/A/C", "C", "1", "3"),
        // A reference number that is not the count, and the two sides exchanged.
        (&REPEATED, "C/G/G/T", "T", "1", "3"),
        (&REPEATED, "C/G/G/T", "T", "3", "2"),
        // The repeated alternate itself, whose place and count are one number.
        (&REPEATED, "C/G/G/T", "G", "1", "2"),
        // A list that drops no entry, and an allele the list does not hold.
        (&PLAIN, "C/G/T", "T", "2", "3"),
        (&REPEATED, "C/G/G/T", "A", "2", "3"),
    ];
    for format in FORMATS {
        for (record, list, allele, reference, candidate) in cases {
            assert_eq!(
                compared(
                    &[CLASS],
                    format,
                    record,
                    &[row(allele, reference, Some(list))],
                    &[row(allele, candidate, Some(list))]
                ),
                left(),
                "{format} {list} {allele} {reference} {candidate}"
            );
        }
        // A list the two rows write differently.
        assert_eq!(
            compared(
                &[CLASS],
                format,
                &REPEATED,
                &[row("T", "2", Some("C/G/G/T"))],
                &[row("T", "3", Some("C/G/T/T"))]
            ),
            left(),
            "{format}"
        );
    }
}

/// The class reads the two fields of the allele number and no record beside the row; the
/// --masks help names it.
#[test]
fn the_class_tables_and_the_help_carry_the_class() {
    let tables = classes::tables();
    assert_eq!(CLASS, "allele_number_repeated_or_reference_equal_alt");
    assert_eq!(
        tables.field_classes_of(CLASS),
        Some(&ALLELE_NUMBER_FIELDS[..])
    );
    assert!(!tables.in_record_reading_classes(CLASS));
    assert!(!tables.in_allele_reading_classes(CLASS));
    let help: Vec<&str> = include_str!("../data/help.txt")
        .split_whitespace()
        .collect();
    assert!(help.contains(&CLASS));
}

/// The two records of the corpus whose list is `C/G/G/T`, each with one row of the allele
/// `T` that the stored reference numbers 2.
const STORED_RECORDS: [&str; 2] = ["dup_alt", "trailing_tabs"];

/// A stored reference output with the allele number of the `T` row of each of
/// `STORED_RECORDS` written as 3, the allele's place in the record, and the count of the
/// rows rewritten.
fn with_written_places(format: &str, text: &str) -> (String, usize) {
    let mut rewritten = 0;
    let mut lines: Vec<String> = Vec::new();
    // The place of a tab column and of a CSQ field, read from the output's own header.
    let mut number_at = None;
    for line in text.lines() {
        let cells: Vec<&str> = line.split('\t').collect();
        let line = match format {
            "default" if STORED_RECORDS.contains(&cells[0]) && cells.get(2) == Some(&"T") => {
                assert!(line.contains("ALLELE_NUM=2;UPLOADED_ALLELE=C/G/G/T;"));
                rewritten += 1;
                line.replace("ALLELE_NUM=2;", "ALLELE_NUM=3;")
            }
            "tab" if line.starts_with("#Uploaded_variation") => {
                number_at = cells.iter().position(|name| *name == "ALLELE_NUM");
                line.to_string()
            }
            "tab" if STORED_RECORDS.contains(&cells[0]) && cells.get(2) == Some(&"T") => {
                let at = number_at.expect("a column of the allele number");
                assert_eq!((cells[at], cells[at + 1]), ("2", "C/G/G/T"));
                let mut cells = cells.clone();
                cells[at] = "3";
                rewritten += 1;
                cells.join("\t")
            }
            "vcf" if line.starts_with("##INFO=<ID=CSQ") => {
                let fields = line.split("Format: ").nth(1).expect("the CSQ fields");
                number_at = fields.split('|').position(|name| name == "ALLELE_NUM");
                line.to_string()
            }
            "vcf"
                if !line.starts_with('#')
                    && cells.get(2).is_some_and(|id| STORED_RECORDS.contains(id)) =>
            {
                let at = number_at.expect("a field of the allele number");
                let (before, csq) = line.split_once("CSQ=").expect("a CSQ");
                let entries: Vec<String> = csq
                    .split(',')
                    .map(|entry| {
                        let mut fields: Vec<&str> = entry.split('|').collect();
                        if fields[0] != "T" {
                            return entry.to_string();
                        }
                        assert_eq!((fields[at], fields[at + 1]), ("2", "C/G/G/T"));
                        fields[at] = "3";
                        rewritten += 1;
                        fields.join("|")
                    })
                    .collect();
                format!("{before}CSQ={}", entries.join(","))
            }
            "json" if !line.is_empty() => {
                let parsed = pyjson::loads(line).expect("a line of the stored output is JSON");
                let mut record = parsed.into_dict().expect("an object");
                let id = record["id"].as_str().expect("an id").to_string();
                if !STORED_RECORDS.contains(&id.as_str()) {
                    lines.push(line.to_string());
                    continue;
                }
                let Some(Value::List(rows)) = record.get_mut("transcript_consequences") else {
                    panic!("{id} has no transcript rows");
                };
                for row in rows {
                    let Value::Dict(row) = row else {
                        panic!("a consequence is no object");
                    };
                    if row["variant_allele"].as_str() != Some("T") {
                        continue;
                    }
                    assert_eq!(row["allele_num"], Value::Int(2));
                    assert_eq!(row["uploaded_allele"].as_str(), Some("C/G/G/T"));
                    row.insert("allele_num".to_string(), Value::Int(3));
                    rewritten += 1;
                }
                pyjson::dumps_compact(&Value::Dict(record))
            }
            _ => line.to_string(),
        };
        lines.push(line);
    }
    (lines.join("\n") + "\n", rewritten)
}

/// The stored reference of the corpus that holds the two records, in each format, against
/// itself with the two rows numbered by their place: both are taken under the class, nothing
/// else differs, and without the class the two are differences.
#[test]
fn the_stored_reference_is_green_against_its_written_places_under_the_class() {
    let corpus = golden().join("GRCh38-vcf-edge-cases");
    let inputs = input_records(&corpus.join("variants.vcf")).expect("the input reads");
    for (format, name) in [
        ("default", "default.txt.gz"),
        ("tab", "tab.txt.gz"),
        ("vcf", "vcf.vcf.gz"),
        ("json", "json.jsonl.gz"),
    ] {
        let reference_text =
            read_text(&corpus.join("expected").join(name)).expect("the stored output reads");
        let (candidate_text, rewritten) = with_written_places(format, &reference_text);
        assert_eq!(rewritten, STORED_RECORDS.len(), "{format}");
        let parsed = |text: &str| {
            match format {
                "default" => parse_default(text),
                "tab" => parse_tab(text),
                "vcf" => parse_vcf(text, &inputs),
                _ => parse_json(text, &inputs),
            }
            .expect("the output parses")
        };
        let (reference, candidate) = (parsed(&reference_text), parsed(&candidate_text));
        assert_eq!(
            under(&[CLASS], &reference, &candidate),
            ("green".to_string(), count_map(&[(CLASS, 2)])),
            "{format}"
        );
        assert_eq!(under(&[], &reference, &candidate).0, "red", "{format}");
    }
}
