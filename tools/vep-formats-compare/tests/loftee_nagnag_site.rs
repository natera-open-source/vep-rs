// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The LoFTEE class of a `NAGNAG_SITE` flag read on the first transcript's strand: two
//! made-up records whose acceptor base a forward and a reverse transcript share, read through
//! the comparison of planted rows in the four formats.

mod compare_common;
mod support;

use compare_common::*;
use std::collections::BTreeMap;
use vep_formats_compare::classes;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::Side;

const CLASS: &str = "loftee_nagnag_site_first_transcript_strand";

const ACCEPTOR: &str = "splice_acceptor_variant";

const FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

/// One transcript row of a record: its transcript, its strand, its terms and its four
/// LoFTEE values, each empty where the row carries none.
#[derive(Clone, Copy)]
struct Row<'a> {
    feature: &'a str,
    strand: &'a str,
    terms: &'a str,
    lof: &'a str,
    filter: &'a str,
    flags: &'a str,
    info: &'a str,
}

/// The forward transcript's acceptor row with the flags `flags`.
fn forward(flags: &str) -> Row<'_> {
    Row {
        feature: "FWDT1",
        strand: "1",
        terms: ACCEPTOR,
        lof: "HC",
        filter: "",
        flags,
        info: "INTRON_SIZE:200",
    }
}

/// The reverse transcript's acceptor row with the flags `flags`.
fn reverse(flags: &str) -> Row<'_> {
    Row {
        feature: "REVT1",
        strand: "-1",
        terms: ACCEPTOR,
        lof: "HC",
        filter: "",
        flags,
        info: "INTRON_SIZE:221",
    }
}

/// One input record of the made-up sequence `toy` and its transcript rows.
#[derive(Clone)]
struct Record<'a> {
    name: &'a str,
    position: &'a str,
    reference: &'a str,
    allele: &'a str,
    rows: Vec<Row<'a>>,
}

impl Record<'_> {
    /// The record's line of the input VCF.
    fn input(&self) -> String {
        format!(
            "toy\t{}\t{}\t{}\t{}\t.\t.\t.",
            self.position, self.name, self.reference, self.allele
        )
    }

    fn location(&self) -> String {
        format!("toy:{}", self.position)
    }
}

/// The record at 400, whose forward row carries `NAGNAG_SITE` on both sides, with the
/// reverse row's flags `flags`.
fn at_400(flags: &str) -> Record<'_> {
    Record {
        name: "acceptor_forward_motif",
        position: "400",
        reference: "G",
        allele: "C",
        rows: vec![forward("NAGNAG_SITE"), reverse(flags)],
    }
}

/// The record at 699, whose forward row carries `NON_CAN_SPLICE` on both sides, with the
/// reverse row's flags `flags`.
fn at_699(flags: &str) -> Record<'_> {
    Record {
        name: "acceptor_reverse_motif",
        position: "699",
        reference: "C",
        allele: "G",
        rows: vec![forward("NON_CAN_SPLICE"), reverse(flags)],
    }
}

/// The two records as the reference writes them and as the candidate does: at 400 the
/// reference's reverse row carries the item and the candidate's does not, at 699 the
/// candidate's alone carries it.
fn both_engines() -> (Vec<Record<'static>>, Vec<Record<'static>>) {
    (
        vec![at_400("NON_CAN_SPLICE,NAGNAG_SITE"), at_699("")],
        vec![at_400("NON_CAN_SPLICE"), at_699("NAGNAG_SITE")],
    )
}

/// The input lines of `records`.
fn inputs_of(records: &[Record]) -> Vec<String> {
    records.iter().map(Record::input).collect()
}

/// The LoFTEE pairs a row's Extra text carries, in the plugin's order.
fn lof_pairs(row: &Row) -> String {
    [
        ("LoF", row.lof),
        ("LoF_filter", row.filter),
        ("LoF_flags", row.flags),
        ("LoF_info", row.info),
    ]
    .iter()
    .filter(|(_, value)| !value.is_empty())
    .map(|(name, value)| format!(";{name}={value}"))
    .collect()
}

/// The default-format output of `records`, a row's strand in its Extra text where it
/// states one.
fn default_output(records: &[Record]) -> String {
    let rows: Vec<String> = records
        .iter()
        .flat_map(|record| {
            record.rows.iter().map(move |row| {
                let strand = if row.strand.is_empty() {
                    String::new()
                } else {
                    format!(";STRAND={}", row.strand)
                };
                default_row(
                    record.name,
                    &record.location(),
                    record.allele,
                    row.feature,
                    row.terms,
                    &format!("IMPACT=HIGH{strand}{}", lof_pairs(row)),
                )
            })
        })
        .collect();
    default_text(&rows)
}

/// A cell of the tab format: `-` for no value.
fn cell(value: &str) -> &str {
    if value.is_empty() {
        "-"
    } else {
        value
    }
}

/// The tab-format output of `records`.
fn tab_output(records: &[Record]) -> String {
    let mut text = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tSTRAND\tLoF\tLoF_filter\tLoF_flags\tLoF_info\n".to_string();
    for record in records {
        for row in &record.rows {
            text.push_str(&format!(
                "{}\t{}\t{}\tENSG1\t{}\tTranscript\t{}\t{}\t{}\t{}\t{}\t{}\n",
                record.name,
                record.location(),
                record.allele,
                row.feature,
                row.terms,
                cell(row.strand),
                cell(row.lof),
                cell(row.filter),
                cell(row.flags),
                cell(row.info),
            ));
        }
    }
    text
}

/// The VCF-format output of `records`: the format joins a list's items with `&`.
fn vcf_output(records: &[Record]) -> String {
    let mut text = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|STRAND|LoF|LoF_filter|LoF_flags|LoF_info\">\n".to_string();
    for record in records {
        let entries: Vec<String> = record
            .rows
            .iter()
            .map(|row| {
                format!(
                    "{}|{}|Transcript|{}|{}|{}|{}|{}|{}",
                    record.allele,
                    row.terms.replace(',', "&"),
                    row.feature,
                    row.strand,
                    row.lof,
                    row.filter.replace(',', "&"),
                    row.flags.replace(',', "&"),
                    row.info.replace(',', "&"),
                )
            })
            .collect();
        let input = record.input();
        text.push_str(&format!(
            "{}CSQ={}\n",
            &input[..input.len() - 1],
            entries.join(",")
        ));
    }
    text
}

/// The JSON-format output of `records`: a consequence object leaves out a field that
/// holds no value.
fn json_output(records: &[Record]) -> String {
    let mut text = String::new();
    for record in records {
        let consequences: Vec<String> = record
            .rows
            .iter()
            .map(|row| {
                let terms: Vec<String> = row
                    .terms
                    .split(',')
                    .map(|term| format!("\"{term}\""))
                    .collect();
                let mut members = vec![
                    format!("\"variant_allele\": \"{}\"", record.allele),
                    format!("\"transcript_id\": \"{}\"", row.feature),
                    format!("\"consequence_terms\": [{}]", terms.join(", ")),
                ];
                if !row.strand.is_empty() {
                    members.push(format!("\"strand\": {}", row.strand));
                }
                for (name, value) in [
                    ("lof", row.lof),
                    ("lof_filter", row.filter),
                    ("lof_flags", row.flags),
                    ("lof_info", row.info),
                ] {
                    if !value.is_empty() {
                        members.push(format!("\"{name}\": \"{value}\""));
                    }
                }
                format!("{{{}}}", members.join(", "))
            })
            .collect();
        text.push_str(&format!(
            "{{\"input\": \"{}\", \"id\": \"{}\", \"allele_string\": \"{}/{}\", \"strand\": 1, \"most_severe_consequence\": \"{ACCEPTOR}\", \"transcript_consequences\": [{}]}}\n",
            record.input().replace('\t', "\\t"),
            record.name,
            record.reference,
            record.allele,
            consequences.join(", ")
        ));
    }
    text
}

/// The entries of `records` in the format `format`.
fn entries(format: &str, records: &[Record]) -> Vec<Entry> {
    let inputs = inputs_of(records);
    match format {
        "default" => parse_default(&default_output(records)),
        "tab" => parse_tab(&tab_output(records)),
        "vcf" => parse_vcf(&vcf_output(records), &inputs),
        _ => parse_json(&json_output(records), &inputs),
    }
    .expect("the output parses")
}

/// What a comparison ended on: its verdict, its reasons, what it applied and, per column
/// that counts entries set aside, those counts by class.
#[derive(Debug)]
struct Outcome {
    verdict: String,
    reasons: Vec<String>,
    applied: BTreeMap<String, i64>,
    masked: BTreeMap<String, BTreeMap<String, i64>>,
}

impl Outcome {
    /// Whether the comparison is green with `count` entries counted under the class and
    /// nothing under another.
    fn takes(&self, count: i64) -> bool {
        self.verdict == "green" && self.applied == count_map(&[(CLASS, count)])
    }

    /// Whether the comparison is red and counts nothing under any class.
    fn refuses(&self) -> bool {
        self.verdict == "red" && self.applied.is_empty() && self.masked.is_empty()
    }
}

/// The columns of `report` that count entries set aside, each with its counts by class.
fn masked_columns(report: &Value) -> BTreeMap<String, BTreeMap<String, i64>> {
    at(report, &["columns"])
        .as_dict()
        .expect("the columns by name")
        .iter()
        .filter_map(|(column, held)| {
            let masked = held.get("masked")?.as_dict().expect("the counts by class");
            let counted = masked
                .iter()
                .map(|(class, count)| match count {
                    Value::Int(count) => (class.clone(), *count),
                    other => panic!("{class} counts {other:?}"),
                })
                .collect();
            Some((column.clone(), counted))
        })
        .collect()
}

/// What two entry lists compare to under the built-in `classes`, with the input's lines
/// where they are given.
fn compared(
    expected: &[Entry],
    actual: &[Entry],
    classes: &[&str],
    inputs: Option<&[String]>,
) -> Outcome {
    let mut documented = documented_with(classes);
    let side = Side::default();
    let options = CompareOptions {
        input_lines: inputs,
        ..CompareOptions::new(&side)
    };
    let report =
        compare_entries(expected, actual, &mut documented, &options).expect("the entries compare");
    Outcome {
        verdict: verdict(&report).to_string(),
        reasons: reasons(&report)
            .iter()
            .map(|reason| reason.to_string())
            .collect(),
        applied: counts(&documented.applied),
        masked: masked_columns(&report),
    }
}

/// What the reference's `reference` and the candidate's `candidate` compare to in `format`
/// under `classes`, the input's lines read.
fn compared_in(
    format: &str,
    reference: &[Record],
    candidate: &[Record],
    classes: &[&str],
) -> Outcome {
    let inputs = inputs_of(reference);
    compared(
        &entries(format, reference),
        &entries(format, candidate),
        classes,
        Some(&inputs),
    )
}

/// Whether the class takes `count` entries of the pair in each of the four formats.
fn taken(reference: &[Record], candidate: &[Record], count: i64) -> bool {
    FORMATS
        .iter()
        .all(|format| compared_in(format, reference, candidate, &[CLASS]).takes(count))
}

/// Whether the class leaves the pair a difference in each of the four formats.
fn refused(reference: &[Record], candidate: &[Record]) -> bool {
    FORMATS
        .iter()
        .all(|format| compared_in(format, reference, candidate, &[CLASS]).refuses())
}

/// The record at 400 with the rows `rows`.
fn record_of(rows: Vec<Row>) -> Vec<Record> {
    vec![Record { rows, ..at_400("") }]
}

/// The pair of the record at 400 whose reverse row is `reference` on one side and
/// `candidate` on the other, beside the forward transcript's acceptor row.
fn reverse_rows<'a>(reference: Row<'a>, candidate: Row<'a>) -> (Vec<Record<'a>>, Vec<Record<'a>>) {
    let beside = forward("NAGNAG_SITE");
    (
        record_of(vec![beside, reference]),
        record_of(vec![beside, candidate]),
    )
}

/// The flags of the reverse row of the record at 400 on the side that carries the item.
const CARRIED: &str = "NON_CAN_SPLICE,NAGNAG_SITE";

/// The flags of the same row on the side that does not.
const SPARED: &str = "NON_CAN_SPLICE";

/// At 400 the reference's reverse row carries the item and the candidate's does not, and
/// at 699 the candidate's alone carries it: each format counts the two entries under the
/// class and sets aside their `LoF_flags` alone. Without the class the same two entries
/// are the comparison's only differences, the default format reporting the pair with the
/// Extra text that holds it.
#[test]
fn the_flag_one_side_alone_carries_is_set_aside_in_every_format() {
    let (reference, candidate) = both_engines();
    let twice = count_map(&[(CLASS, 2)]);
    let columns: [(&str, &[&str]); 4] = [
        ("default", &["Extra", "LoF_flags"]),
        ("tab", &["LoF_flags"]),
        ("vcf", &["LoF_flags"]),
        ("json", &["lof_flags"]),
    ];
    for (format, differing) in columns {
        let outcome = compared_in(format, &reference, &candidate, &[CLASS]);
        assert!(outcome.takes(2), "{format}: {outcome:?}");
        assert!(outcome.reasons.is_empty(), "{format}");
        let masked: BTreeMap<String, BTreeMap<String, i64>> = differing
            .iter()
            .map(|column| (column.to_string(), twice.clone()))
            .collect();
        assert_eq!(outcome.masked, masked, "{format}");

        let unnamed = compared_in(format, &reference, &candidate, &[]);
        assert!(unnamed.refuses(), "{format}: {unnamed:?}");
        let reasons: Vec<String> = differing
            .iter()
            .map(|column| format!("field {column}: 2 differing entries"))
            .collect();
        assert_eq!(unnamed.reasons, reasons, "{format}");
    }
    // Each record is one entry: the reference's side carrying the item, then the candidate's.
    for record in 0..2 {
        assert!(
            taken(&reference[record..=record], &candidate[record..=record], 1),
            "{record}"
        );
    }
    // The row that differs stands on either strand: here the forward row, beside a reverse one.
    let lone = |flags| record_of(vec![reverse(SPARED), forward(flags)]);
    assert!(taken(&lone("NAGNAG_SITE"), &lone(""), 1));
    assert!(taken(&lone(SPARED), &lone("NAGNAG_SITE,NON_CAN_SPLICE"), 1));
}

/// A list that differs in another item, with or without the class's own, a list whose
/// other items stand in another order and a side that carries the item twice stay
/// differences; a pair whose two values are equal is no entry at all.
#[test]
fn a_flags_list_that_differs_otherwise_stays_a_difference() {
    let cases = [
        (CARRIED, "NAGNAG_SITE"),
        (CARRIED, "PHYLOCSF_WEAK"),
        ("NAGNAG_SITE", SPARED),
        (CARRIED, ""),
        (
            "NON_CAN_SPLICE,PHYLOCSF_WEAK,NAGNAG_SITE",
            "PHYLOCSF_WEAK,NON_CAN_SPLICE",
        ),
        (CARRIED, "NAGNAG_SITE,NON_CAN_SPLICE"),
        ("NAGNAG_SITE,NON_CAN_SPLICE,NAGNAG_SITE", SPARED),
        ("NAGNAG_SITE,NAGNAG_SITE", ""),
    ];
    for (reference_flags, candidate_flags) in cases {
        for (ours, theirs) in [
            (reference_flags, candidate_flags),
            (candidate_flags, reference_flags),
        ] {
            let (reference, candidate) = reverse_rows(reverse(ours), reverse(theirs));
            assert!(refused(&reference, &candidate), "{ours} against {theirs}");
        }
    }
    for flags in [CARRIED, SPARED, ""] {
        let (reference, candidate) = reverse_rows(reverse(flags), reverse(flags));
        for format in FORMATS {
            let outcome = compared_in(format, &reference, &candidate, &[CLASS]);
            assert_eq!(outcome.verdict, "green", "{format}: {flags}");
            assert!(outcome.applied.is_empty() && outcome.masked.is_empty());
        }
    }
}

/// The class reads a `splice_acceptor_variant` row alone, whatever its other terms: the
/// same pair on a `splice_donor_variant` row and on a row one side of which lacks the term
/// stays a difference. So does a pair beside a `LoF`, a `LoF_filter` or a `LoF_info`
/// that differs.
#[test]
fn a_row_of_another_term_or_beside_another_loftee_difference_stays_a_difference() {
    let carrying = reverse(CARRIED);
    let spared = reverse(SPARED);
    let (reference, candidate) = reverse_rows(carrying, spared);
    assert!(taken(&reference, &candidate, 1));
    let with_terms = |row: Row<'static>, terms: &'static str| Row { terms, ..row };
    let several = "splice_acceptor_variant,intron_variant";
    let (reference, candidate) =
        reverse_rows(with_terms(carrying, several), with_terms(spared, several));
    assert!(taken(&reference, &candidate, 1));
    let donor = "splice_donor_variant";
    let pairs = [
        (with_terms(carrying, donor), with_terms(spared, donor)),
        (carrying, with_terms(spared, "intron_variant")),
        (with_terms(carrying, "intron_variant"), spared),
        (
            carrying,
            Row {
                lof: "LC",
                ..spared
            },
        ),
        (
            carrying,
            Row {
                filter: "NON_ACCEPTOR_DISRUPTING",
                ..spared
            },
        ),
        (
            Row {
                filter: "NON_ACCEPTOR_DISRUPTING",
                lof: "LC",
                ..carrying
            },
            Row {
                filter: "NON_ACCEPTOR_DISRUPTING,5UTR_SPLICE",
                lof: "LC",
                ..spared
            },
        ),
        (
            carrying,
            Row {
                info: "INTRON_SIZE:222",
                ..spared
            },
        ),
        (carrying, Row { info: "", ..spared }),
    ];
    for (at, (reference_row, candidate_row)) in pairs.into_iter().enumerate() {
        let (reference, candidate) = reverse_rows(reference_row, candidate_row);
        assert!(refused(&reference, &candidate), "{at}");
    }
    // An equal filter and an equal item list beside the pair leave it the class's.
    let filtered = |row: Row<'static>| Row {
        filter: "NON_ACCEPTOR_DISRUPTING,5UTR_SPLICE",
        lof: "LC",
        info: "INTRON_SIZE:221,MUTANT_DONOR_MES:1",
        ..row
    };
    let (reference, candidate) = reverse_rows(filtered(carrying), filtered(spared));
    assert!(taken(&reference, &candidate, 1));
}

/// The record has an acceptor row on the other strand than the row's own: a record whose
/// other row stands on the row's strand, whose other-strand row is no acceptor row, whose
/// row or other row states no strand, or that has the one row stays a difference, and so
/// does a row whose only other-strand acceptor row is another record's.
#[test]
fn a_record_without_an_acceptor_row_on_the_other_strand_stays_a_difference() {
    let carrying = reverse(CARRIED);
    let spared = reverse(SPARED);
    let same_strand = Row {
        feature: "REVT2",
        ..reverse("NAGNAG_SITE")
    };
    let no_acceptor = Row {
        terms: "intron_variant",
        lof: "",
        info: "",
        ..forward("")
    };
    let strandless = Row {
        strand: "",
        ..forward("NAGNAG_SITE")
    };
    for (at, beside) in [same_strand, no_acceptor, strandless]
        .into_iter()
        .enumerate()
    {
        let reference = record_of(vec![beside, carrying]);
        let candidate = record_of(vec![beside, spared]);
        assert!(refused(&reference, &candidate), "{at}");
    }
    let unstated = |row: Row<'static>| Row { strand: "", ..row };
    let (reference, candidate) = reverse_rows(unstated(carrying), unstated(spared));
    assert!(refused(&reference, &candidate));
    assert!(refused(
        &record_of(vec![carrying]),
        &record_of(vec![spared])
    ));
    let elsewhere = |row: Row<'static>| {
        vec![
            Record {
                rows: vec![row],
                ..at_400("")
            },
            Record {
                rows: vec![forward("NON_CAN_SPLICE")],
                ..at_699("")
            },
        ]
    };
    assert!(refused(&elsewhere(carrying), &elsewhere(spared)));
}

/// A VCF or JSON entry names its record, so its rows are read without the input's lines.
/// A default or tab row under an identifier does not say whether its Location's rows are
/// one record's: the class takes it with the input's lines and leaves it a difference
/// without them, while a row whose Uploaded_variation spells the record's alleles is
/// taken either way.
#[test]
fn a_row_under_an_identifier_is_read_with_the_inputs_lines() {
    let (reference, candidate) = both_engines();
    for format in FORMATS {
        let outcome = compared(
            &entries(format, &reference),
            &entries(format, &candidate),
            &[CLASS],
            None,
        );
        if matches!(format, "vcf" | "json") {
            assert!(outcome.takes(2), "{format}: {outcome:?}");
        } else {
            assert!(outcome.refuses(), "{format}: {outcome:?}");
        }
    }
    let spelt = |records: &[Record<'static>]| -> Vec<Record<'static>> {
        records
            .iter()
            .cloned()
            .map(|record| Record {
                name: if record.position == "400" {
                    "toy_400_G/C"
                } else {
                    "toy_699_C/G"
                },
                ..record
            })
            .collect()
    };
    for format in ["default", "tab"] {
        let outcome = compared(
            &entries(format, &spelt(&reference)),
            &entries(format, &spelt(&candidate)),
            &[CLASS],
            None,
        );
        assert!(outcome.takes(2), "{format}: {outcome:?}");
    }
    // A line under another identifier names no record of the rows.
    let other = vec!["toy\t400\tanother\tG\tC\t.\t.\t.".to_string()];
    let outcome = compared(
        &entries("default", &reference[..1]),
        &entries("default", &candidate[..1]),
        &[CLASS],
        Some(&other),
    );
    assert!(outcome.refuses(), "{outcome:?}");
}

/// The class reads `LoF_flags` under either letter case, stands directly before the two
/// classes that read a row's term pair, which close the table, and is one of the classes
/// that read a record's rows and bring the input to a text-format comparison; the --masks
/// help names it.
#[test]
fn the_class_tables_and_the_help_carry_the_class() {
    let tables = classes::tables();
    let names: Vec<&str> = tables
        .field_classes()
        .iter()
        .map(|(class, _)| *class)
        .collect();
    assert_eq!(
        names[names.len() - 3..],
        [
            CLASS,
            "plugin_value_follows_splice_term_past_the_record_span",
            "plugin_value_follows_stop_codon_term_pair"
        ]
    );
    assert_eq!(
        tables.field_classes_of(CLASS),
        Some(&["LoF_flags", "lof_flags"][..])
    );
    assert!(tables.in_record_reading_classes(CLASS));
    assert!(tables.in_allele_reading_classes(CLASS));
    assert!(!tables.in_sibling_allele_classes(CLASS));
    assert_eq!(classes::LOFTEE_NAGNAG_SITE_FIRST_TRANSCRIPT_STRAND, CLASS);
    let help: Vec<&str> = include_str!("../data/help.txt")
        .split_whitespace()
        .collect();
    assert!(help.join(" ").contains(&format!(
        "{CLASS} (on a splice_acceptor_variant row of a record that has a splice_acceptor_variant row on the other strand, a LoF_flags whose two values hold the same items in the same order but for one NAGNAG_SITE, "
    )));
}

/// The command line names the class like any other, and reads the rows of a default-format
/// comparison under identifiers from the lines of --input.
#[test]
fn the_command_line_names_the_class_and_reads_the_input_beside_it() {
    let (reference, candidate) = both_engines();
    let directory = support::scratch("loftee-nagnag-site");
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let input = format!("{header}{}\n", inputs_of(&reference).join("\n"));
    for (name, text) in [
        ("reference.txt", default_output(&reference)),
        ("candidate.txt", default_output(&candidate)),
        ("variants.vcf", input),
    ] {
        std::fs::write(directory.join(name), text).expect("the file is written");
    }
    let exit_code = |more: &[&str]| {
        let mut args = vec![
            "--format",
            "default",
            "--expected",
            "reference.txt",
            "--actual",
            "candidate.txt",
        ];
        args.extend(more);
        support::run(&directory, &args).code
    };
    assert_eq!(exit_code(&[]), 1);
    assert_eq!(exit_code(&["--input", "variants.vcf"]), 1);
    assert_eq!(exit_code(&["--masks", CLASS]), 1);
    assert_eq!(exit_code(&["--masks", CLASS, "--input", "variants.vcf"]), 0);
}
