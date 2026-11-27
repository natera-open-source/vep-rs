// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvs_added_bases_written_outside_the_transcript` on made rows read in the default, tab,
//! VCF and JSON formats. A default or tab row carries its Location in a column; a VCF or
//! JSON entry has it from its record's position and alleles (`Entry.allele_location`). A
//! default or tab row whose first column spells its record's alleles carries them as a VCF
//! or JSON entry does; a row under an identifier carries none and is read from its Location.

mod compare_common;

use std::collections::BTreeMap;

use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::{parse_default, parse_tab};
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::{self, DumpOptions, Value};
use vep_formats_compare::rules::added_bases;
use vep_formats_compare::rules::downstream::field_of;
use vep_formats_compare::rules::record_alleles::vcf_alleles_as_spelt;
use vep_formats_compare::rules::{field_class_accepts, RuleContext, Side};

const ADDED_OUTSIDE: &str = "hgvs_added_bases_written_outside_the_transcript";
const ADDED_OUTSIDE_FEATURE: &str = "ENST00000000071";
const ADDED_OUTSIDE_FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];
/// The forms one pair of rows is read in: a default row under an identifier, a default row
/// whose first column spells its record's alleles, each also as a tab row, a VCF entry and a
/// JSON entry.
const RECORD_FORMS: [&str; 6] = [
    "identifier",
    "identifier, tab",
    "spelt",
    "spelt, tab",
    "vcf",
    "json",
];
/// The forms of `RECORD_FORMS` whose entries carry their record's alleles.
const FORMS_WITH_ALLELES: [&str; 4] = ["spelt", "spelt, tab", "vcf", "json"];
const ADDED_OUTSIDE_FIELDS: [&str; 4] = [
    "HGVSc",
    "HGVS_OFFSET",
    "HGVS_IntronStartOffset",
    "HGVS_IntronEndOffset",
];
/// The REF and ALT of the input line behind a row of the grid: one base inserted after the
/// base the line shares.
const ADDED_OUTSIDE_RECORD: (&str, &str) = ("G", "GT");
/// The Extra fields a made row may carry, in the order a tab output or a CSQ field writes
/// them.
const WRITTEN_FIELDS: [&str; 8] = [
    "IMPACT",
    "DISTANCE",
    "STRAND",
    "BAM_EDIT",
    "HGVSc",
    "HGVS_OFFSET",
    "HGVS_IntronStartOffset",
    "HGVS_IntronEndOffset",
];
const TERMS: [&str; 3] = [
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intron_variant",
];

/// The candidate's row of a one-base insertion outside the transcript: what it is, its
/// notation after the sequence name, the two intron offsets of its coordinates, its
/// HGVS_OFFSET.
type OutsideCandidate = (
    &'static str,
    Option<&'static str>,
    Option<(&'static str, &'static str)>,
    Option<&'static str>,
);
const NO_OFFSETS: Option<(&str, &str)> = Some(("0", "0"));
const ADDED_OUTSIDE_CANDIDATES: [OutsideCandidate; 13] = [
    ("dup, nothing moved", Some("n.2601dup"), NO_OFFSETS, None),
    ("dup, moved", Some("n.14dup"), NO_OFFSETS, Some("14")),
    (
        "dup, moved the other way",
        Some("n.14dup"),
        NO_OFFSETS,
        Some("-14"),
    ),
    (
        "dup into an intron",
        Some("c.20+6dup"),
        Some(("6", "6")),
        Some("26"),
    ),
    (
        "ins, moved",
        Some("c.20+8_20+9insT"),
        Some(("8", "9")),
        Some("28"),
    ),
    (
        "ins of another count",
        Some("c.20+8_20+9insTT"),
        Some(("8", "9")),
        Some("28"),
    ),
    (
        "dup of another count",
        Some("n.13_14dup"),
        NO_OFFSETS,
        Some("14"),
    ),
    (
        "fields not its notation's",
        Some("c.20+8_20+9insT"),
        Some(("9", "9")),
        Some("28"),
    ),
    (
        "offset within the distance",
        Some("n.14dup"),
        NO_OFFSETS,
        Some("3"),
    ),
    ("del", Some("n.14del"), NO_OFFSETS, Some("14")),
    ("no notation", None, NO_OFFSETS, None),
    ("no fields", Some("n.14dup"), None, Some("14")),
    (
        "offset the distance itself",
        Some("n.14dup"),
        NO_OFFSETS,
        Some("5"),
    ),
];

/// The reference's row: no field, a pair, one field alone, an HGVSc or an HGVS_OFFSET of
/// its own.
#[derive(Clone, Copy, Debug, PartialEq)]
enum OutsideReference {
    Nothing,
    Fields(&'static str, Option<&'static str>),
    Notation,
    Offset,
}
const ADDED_OUTSIDE_REFERENCES: [OutsideReference; 7] = [
    OutsideReference::Nothing,
    OutsideReference::Fields("0", Some("0")),
    OutsideReference::Fields("1", Some("0")),
    OutsideReference::Fields("9", Some("9")),
    OutsideReference::Fields("0", None),
    OutsideReference::Notation,
    OutsideReference::Offset,
];

/// One cell of the made-up grid of the class: the row's term; its strand; the candidate's
/// row; and the reference's row.
#[derive(Clone, Copy, Debug)]
struct AddedOutsideCell {
    term: &'static str,
    strand: Option<&'static str>,
    candidate: OutsideCandidate,
    reference: OutsideReference,
}

fn added_outside_grid() -> Vec<AddedOutsideCell> {
    let mut cells = Vec::new();
    for term in TERMS {
        for strand in [Some("1"), Some("-1"), None] {
            for candidate in ADDED_OUTSIDE_CANDIDATES {
                for reference in ADDED_OUTSIDE_REFERENCES {
                    cells.push(AddedOutsideCell {
                        term,
                        strand,
                        candidate,
                        reference,
                    });
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a field of a cell: the row is outside the
/// transcript, the reference prints no HGVSc and no HGVS_OFFSET and writes both plugin
/// fields or neither, the candidate prints a `dup` or `ins` of one base with the intron
/// offsets of its coordinates, or neither row writes a plugin field, beside no HGVS_OFFSET
/// on a downstream row (a `dup`) and one larger than the DISTANCE of 5 on an upstream row,
/// negative on a reverse strand; then the HGVSc, the HGVS_OFFSET and each plugin field the
/// two rows write differently.
fn added_outside_set_aside(cell: &AddedOutsideCell, field: &str) -> bool {
    let (kind, _, written, offset) = cell.candidate;
    let reference = match cell.reference {
        OutsideReference::Nothing => None,
        OutsideReference::Fields(start, Some(end)) => Some((start, end)),
        _ => return false,
    };
    if cell.term == "intron_variant" {
        return false;
    }
    if kind == "no fields" {
        // A run without the plugin writes the two fields on neither row.
        if reference.is_some() {
            return false;
        }
    } else if ![
        "dup, nothing moved",
        "dup, moved",
        "dup, moved the other way",
        "dup into an intron",
        "ins, moved",
    ]
    .contains(&kind)
    {
        return false;
    }
    if cell.term == "downstream_gene_variant" {
        if kind != "dup, nothing moved" {
            return false;
        }
    } else {
        let Some(offset) = offset else {
            return false;
        };
        if (cell.strand == Some("1") && offset.starts_with('-'))
            || (cell.strand == Some("-1") && !offset.starts_with('-'))
        {
            return false;
        }
    }
    if field == "HGVS_OFFSET" {
        return offset.is_some();
    }
    match (field, reference, written) {
        ("HGVS_IntronStartOffset", Some(held), Some(written)) => held.0 != written.0,
        ("HGVS_IntronEndOffset", Some(held), Some(written)) => held.1 != written.1,
        _ => true,
    }
}

/// The reference's and the candidate's default-format rows of a cell: a one-base
/// insertion five bases from the transcript, under a Location of the cell's own.
fn added_outside_rows(cell: &AddedOutsideCell, number: usize) -> [String; 2] {
    let (_, notation, written, offset) = cell.candidate;
    let position = 7_000_000 + 10 * number;
    let columns = format!(
        "cell{number}\t2:{position}-{}\tT\tENSG1\t{ADDED_OUTSIDE_FEATURE}\tTranscript\t{}\t-\t-\t-\t-\t-\t-",
        position + 1,
        cell.term
    );
    let extra = match cell.strand {
        Some(strand) => format!("IMPACT=MODIFIER;DISTANCE=5;STRAND={strand}"),
        None => "IMPACT=MODIFIER;DISTANCE=5".to_string(),
    };
    let mut reference_text = extra.clone();
    match cell.reference {
        OutsideReference::Notation => {
            reference_text.push_str(&format!(";HGVSc={ADDED_OUTSIDE_FEATURE}.1:n.1_2insT"));
        }
        OutsideReference::Offset => reference_text.push_str(";HGVS_OFFSET=2"),
        OutsideReference::Fields(start, end) => {
            if let Some(end) = end {
                reference_text.push_str(&format!(";HGVS_IntronEndOffset={end}"));
            }
            reference_text.push_str(&format!(";HGVS_IntronStartOffset={start}"));
        }
        OutsideReference::Nothing => {}
    }
    let mut candidate_text = extra;
    if let Some(notation) = notation {
        candidate_text.push_str(&format!(";HGVSc={ADDED_OUTSIDE_FEATURE}.1:{notation}"));
    }
    if let Some(offset) = offset {
        candidate_text.push_str(&format!(";HGVS_OFFSET={offset}"));
    }
    if let Some((start, end)) = written {
        candidate_text.push_str(&format!(
            ";HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}"
        ));
    }
    [
        format!("{columns}\t{reference_text}"),
        format!("{columns}\t{candidate_text}"),
    ]
}

/// The columns of a default-format row before its Extra text, and the fields of that text
/// in the row's order.
type SplitRow<'a> = (Vec<&'a str>, Vec<(&'a str, &'a str)>);

fn split_row(row: &str) -> SplitRow<'_> {
    let mut columns: Vec<&str> = row.split('\t').collect();
    let extra = columns.pop().expect("an Extra column");
    let fields = extra
        .split(';')
        .map(|item| item.split_once('=').expect("a name and a value"))
        .collect();
    (columns, fields)
}

fn held<'a>(fields: &[(&str, &'a str)], name: &str) -> Option<&'a str> {
    fields
        .iter()
        .find(|(field, _)| *field == name)
        .map(|(_, value)| *value)
}

fn member(name: &str, value: Value) -> (String, Value) {
    (name.to_string(), value)
}

fn text(value: &str) -> Value {
    Value::Str(value.to_string())
}

/// The reference's and the candidate's outputs of a pair of default-format rows in `fmt`,
/// and the input lines they are read against: the rows themselves, their tab form (each
/// Extra field that a row carries a column of both outputs), or the VCF lines and the JSON
/// objects of the input line that writes them, whose REF and ALT are `record` at the
/// Location's first position. A CSQ field stands for every Extra field either row carries;
/// a JSON object carries its own row's fields.
fn added_outside_outputs(
    fmt: &str,
    rows: &[String; 2],
    record: (&str, &str),
) -> (Vec<String>, Vec<String>) {
    if fmt == "default" {
        return (
            rows.iter().map(|row| default_text(&[row])).collect(),
            Vec::new(),
        );
    }
    let split: Vec<SplitRow> = rows.iter().map(|row| split_row(row)).collect();
    let names: Vec<&str> = WRITTEN_FIELDS
        .into_iter()
        .filter(|name| split.iter().any(|(_, fields)| held(fields, name).is_some()))
        .collect();
    if fmt == "tab" {
        let columns: String = names.iter().map(|name| format!("\t{name}")).collect();
        let header = DEFAULT_HEADER.replace("\tExtra", &columns);
        let outputs = split
            .iter()
            .map(|(cells, fields)| {
                let mut cells = cells.clone();
                cells.extend(names.iter().map(|name| held(fields, name).unwrap_or("-")));
                format!("{header}\n{}\n", cells.join("\t"))
            })
            .collect();
        return (outputs, Vec::new());
    }
    let first = &split[0].0;
    let (name, location, allele, feature, term) =
        (first[0], first[1], first[2], first[4], first[6]);
    let (chromosome, span) = location.split_once(':').expect("a chromosome and a span");
    let position: i64 = span
        .split('-')
        .next()
        .and_then(|start| start.parse().ok())
        .expect("a position");
    let (reference, alternates) = record;
    let line = format!("{chromosome}\t{position}\t{name}\t{reference}\t{alternates}\t.\t.\t.");
    if fmt == "vcf" {
        let header = format!(
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|{}\">\n",
            names.join("|")
        );
        let outputs = split
            .iter()
            .map(|(_, fields)| {
                let mut cells = vec![allele, term, feature, "Transcript"];
                cells.extend(names.iter().map(|name| held(fields, name).unwrap_or("")));
                format!(
                    "{header}{}CSQ={}\n",
                    line.strip_suffix('.').expect("an empty INFO column"),
                    cells.join("|")
                )
            })
            .collect();
        return (outputs, vec![line]);
    }
    let alternates: Vec<String> = alternates.split(',').map(str::to_string).collect();
    let (spelt_reference, spelt_alternates) = vcf_alleles_as_spelt(reference, &alternates);
    let mut alleles = vec![spelt_reference.clone()];
    alleles.extend(spelt_alternates);
    let outputs = split
        .iter()
        .map(|(_, fields)| {
            let mut item: IndexMap<String, Value> = IndexMap::from([
                member("variant_allele", text(allele)),
                member("transcript_id", text(feature)),
                member("consequence_terms", Value::List(vec![text(term)])),
            ]);
            for (field, value) in fields {
                let whole = value.strip_prefix('-').unwrap_or(value);
                let number = !whole.is_empty() && whole.bytes().all(|byte| byte.is_ascii_digit());
                let written = match value.parse::<i64>() {
                    Ok(number_written) if number => Value::Int(number_written),
                    _ => text(value),
                };
                item.insert(field.to_lowercase(), written);
            }
            let object = Value::Dict(IndexMap::from([
                member("input", text(&line)),
                member("id", text(name)),
                member("seq_region_name", text(chromosome)),
                member(
                    "start",
                    Value::Int(position + i64::from(spelt_reference != reference)),
                ),
                member("end", Value::Int(position + reference.len() as i64 - 1)),
                member("strand", Value::Int(1)),
                member("allele_string", text(&alleles.join("/"))),
                member("most_severe_consequence", text(term)),
                member(
                    "transcript_consequences",
                    Value::List(vec![Value::Dict(item)]),
                ),
            ]));
            pyjson::dumps(&object, &DumpOptions::default()) + "\n"
        })
        .collect();
    (outputs, vec![line])
}

/// The entries of one output of `added_outside_outputs`.
fn added_outside_parse(fmt: &str, output: &str, lines: &[String]) -> Vec<Entry> {
    match fmt {
        "vcf" => parse_vcf(output, lines),
        "json" => parse_json(output, lines),
        "default" => parse_default(output),
        _ => parse_tab(output),
    }
    .expect("the output parses")
}

/// The entry of the grid's transcript in one output read in `fmt`.
fn entry_of(fmt: &str, output: &str, lines: &[String]) -> Entry {
    added_outside_parse(fmt, output, lines)
        .into_iter()
        .find(|entry| entry.feature == ADDED_OUTSIDE_FEATURE)
        .expect("an entry of the grid's transcript")
}

/// The reference's and the candidate's entries of a pair of default-format rows read in
/// `fmt`.
fn added_outside_entries(fmt: &str, rows: &[String; 2], record: (&str, &str)) -> [Entry; 2] {
    let (outputs, lines) = added_outside_outputs(fmt, rows, record);
    [
        entry_of(fmt, &outputs[0], &lines),
        entry_of(fmt, &outputs[1], &lines),
    ]
}

/// The name Ensembl VEP gives the rows of an input line that carries no identifier: the
/// chromosome, the first position of the alleles as its entries spell them, and those
/// alleles.
fn spelt_name(rows: &[String; 2], record: (&str, &str)) -> String {
    let (columns, _) = split_row(&rows[0]);
    let (chromosome, span) = columns[1].split_once(':').expect("a chromosome and a span");
    let position: i64 = span
        .split('-')
        .next()
        .and_then(|start| start.parse().ok())
        .expect("a position");
    let (reference, alternates) = record;
    let alternates: Vec<String> = alternates.split(',').map(str::to_string).collect();
    let (spelt_reference, spelt_alternates) = vcf_alleles_as_spelt(reference, &alternates);
    let start = position + i64::from(spelt_reference != reference);
    let mut alleles = vec![spelt_reference];
    alleles.extend(spelt_alternates);
    format!("{chromosome}_{start}_{}", alleles.join("/"))
}

/// The format one of `RECORD_FORMS` is parsed in, and the reference's and the candidate's
/// outputs and input lines of a pair of default-format rows under an identifier read in that
/// form (`added_outside_outputs`), the first column of a spelt form the name of
/// `spelt_name`.
fn record_form_outputs<'a>(
    form: &'a str,
    rows: &[String; 2],
    record: (&str, &str),
) -> (&'a str, Vec<String>, Vec<String>) {
    let (spelt, fmt) = match form {
        "identifier" => (false, "default"),
        "identifier, tab" => (false, "tab"),
        "spelt" => (true, "default"),
        "spelt, tab" => (true, "tab"),
        fmt => (false, fmt),
    };
    let rows = if spelt {
        let name = spelt_name(rows, record);
        rows.clone().map(|row| {
            let (_, rest) = row.split_once('\t').expect("a first column");
            format!("{name}\t{rest}")
        })
    } else {
        rows.clone()
    };
    let (outputs, lines) = added_outside_outputs(fmt, &rows, record);
    (fmt, outputs, lines)
}

/// The reference's and the candidate's entries of a pair of default-format rows under an
/// identifier, read in one of `RECORD_FORMS`.
fn record_form_entries(form: &str, rows: &[String; 2], record: (&str, &str)) -> [Entry; 2] {
    let (fmt, outputs, lines) = record_form_outputs(form, rows, record);
    [
        entry_of(fmt, &outputs[0], &lines),
        entry_of(fmt, &outputs[1], &lines),
    ]
}

/// What the class says of the differing field `field` of the reference's entry `e` and the
/// candidate's entry `a`, asked as a comparison without side tables asks it.
fn accepts(e: &Entry, a: &Entry, field: &str) -> bool {
    let side = Side::default();
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    field_class_accepts(
        ADDED_OUTSIDE,
        e,
        field,
        &field_of(Some(e), &[field]),
        &field_of(Some(a), &[field]),
        Some(a),
        &context,
    )
    .expect("the class reads the rows")
}

type Verdicts = BTreeMap<&'static str, BTreeMap<&'static str, bool>>;

/// Per form of `RECORD_FORMS`, what the class says of each of its fields that the two
/// entries of a pair of default-format rows write differently.
fn added_outside_verdicts(rows: &[String; 2], record: (&str, &str)) -> Verdicts {
    RECORD_FORMS
        .into_iter()
        .map(|fmt| {
            let [e, a] = record_form_entries(fmt, rows, record);
            let read = ADDED_OUTSIDE_FIELDS
                .into_iter()
                .filter(|field| field_of(Some(&e), &[field]) != field_of(Some(&a), &[field]))
                .map(|field| (field, accepts(&e, &a, field)))
                .collect();
            (fmt, read)
        })
        .collect()
}

/// Each cell of the made-up grid read in the default, tab, VCF and JSON formats, a default
/// and a tab row under an identifier and under a first column that spells the record's
/// alleles: the class sets aside the same fields of the same cells in each, the ones the
/// rule's statement names, and every other cell stays a difference in each.
#[test]
fn the_grid_of_bases_written_outside_reads_alike_in_every_format() {
    let grid = added_outside_grid();
    assert_eq!(grid.len(), 819);
    let mut wrong = Vec::new();
    for (number, cell) in grid.iter().enumerate() {
        let verdicts =
            added_outside_verdicts(&added_outside_rows(cell, number), ADDED_OUTSIDE_RECORD);
        assert_eq!(verdicts.len(), RECORD_FORMS.len());
        let expected: BTreeMap<&str, bool> = verdicts["identifier"]
            .keys()
            .map(|field| (*field, added_outside_set_aside(cell, field)))
            .collect();
        for (fmt, read) in &verdicts {
            if *read != expected {
                wrong.push((*fmt, *cell));
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "{} cells, the first {:?}",
        wrong.len(),
        wrong.first()
    );

    let cell = AddedOutsideCell {
        term: "upstream_gene_variant",
        strand: Some("1"),
        candidate: ADDED_OUTSIDE_CANDIDATES[1],
        reference: OutsideReference::Nothing,
    };
    let rows = added_outside_rows(&cell, 0);
    let located: Vec<(&str, Option<String>, Option<String>)> = ADDED_OUTSIDE_FORMATS
        .into_iter()
        .map(|fmt| {
            let [e, _] = added_outside_entries(fmt, &rows, ADDED_OUTSIDE_RECORD);
            (fmt, e.location, e.allele_location)
        })
        .collect();
    let span = Some("2:7000000-7000001".to_string());
    assert_eq!(
        located,
        [
            ("default", span.clone(), None),
            ("tab", span.clone(), None),
            ("vcf", None, span.clone()),
            ("json", None, span),
        ]
    );
}

/// The upstream rows whose DISTANCE the class reads against the HGVS_OFFSET, each in every
/// form of `RECORD_FORMS` under the record that writes it. An insertion's row is
/// the class's where the 3' rule carried the base further than the DISTANCE. The row of an
/// allele that adds a base to the two bases or the one base it keeps, whose DISTANCE counts
/// from that span, is the class's from the DISTANCE on. Two bases in place of two, a smaller
/// HGVS_OFFSET, and three inserted bases carried the DISTANCE alone stay differences in each
/// format. So does a VCF entry whose line spells its allele nowhere at the line's span: its
/// record gives it no Location.
#[test]
fn an_allele_longer_than_its_span_reads_alike_in_every_format() {
    let cell = AddedOutsideCell {
        term: "upstream_gene_variant",
        strand: Some("1"),
        candidate: ADDED_OUTSIDE_CANDIDATES[1],
        reference: OutsideReference::Fields("0", Some("0")),
    };
    let inserted = added_outside_rows(&cell, 0);
    let rows_with = |changes: &[(&str, &str)]| -> [String; 2] {
        let mut rows = inserted.clone();
        for (old, new) in changes {
            rows = rows.map(|row| row.replace(old, new));
        }
        rows
    };
    let kept_two = ("\tT\t", "\tGTT\t");
    let kept_one = [("-7000001\t", "\t"), ("\tT\t", "\tGT\t")];
    let the_distance = ("HGVS_OFFSET=14", "HGVS_OFFSET=5");
    let less = ("HGVS_OFFSET=14", "HGVS_OFFSET=4");
    let one_more = ("HGVS_OFFSET=14", "HGVS_OFFSET=6");
    let three_bases = ("n.14dup", "n.12_14dup");
    let shapes = [
        (
            "an insertion carried further than the distance",
            inserted.clone(),
            ("G", "GT"),
            true,
        ),
        (
            "an insertion carried the distance",
            rows_with(&[the_distance]),
            ("G", "GT"),
            false,
        ),
        (
            "a base added to the two kept",
            rows_with(&[kept_two]),
            ("GT", "GTT,CT"),
            true,
        ),
        (
            "a base added to the two kept, carried the distance",
            rows_with(&[kept_two, the_distance]),
            ("GT", "GTT,CT"),
            true,
        ),
        (
            "a base added to the two kept, carried less",
            rows_with(&[kept_two, less]),
            ("GT", "GTT,CT"),
            false,
        ),
        (
            "a base added to the one kept, carried the distance",
            rows_with(&[kept_one[0], kept_one[1], the_distance]),
            ("G", "GT,C"),
            true,
        ),
        (
            "a base added to the one kept, carried less",
            rows_with(&[kept_one[0], kept_one[1], less]),
            ("G", "GT,C"),
            false,
        ),
        (
            "two bases in place of two",
            rows_with(&[("\tT\t", "\tGT\t")]),
            ("AC", "GT"),
            false,
        ),
        (
            "three inserted bases carried the distance",
            rows_with(&[kept_two, three_bases, the_distance]),
            ("A", "AGTT"),
            false,
        ),
        (
            "three inserted bases carried one further",
            rows_with(&[kept_two, three_bases, one_more]),
            ("A", "AGTT"),
            true,
        ),
    ];
    for (label, rows, record, set_aside) in &shapes {
        let fields: BTreeMap<&str, bool> =
            BTreeMap::from([("HGVSc", *set_aside), ("HGVS_OFFSET", *set_aside)]);
        let expected: Verdicts = RECORD_FORMS
            .into_iter()
            .map(|fmt| (fmt, fields.clone()))
            .collect();
        assert_eq!(added_outside_verdicts(rows, *record), expected, "{label}");
    }

    // Ensembl VEP shortens a lone alternate past the base it shares, so the line of the same
    // insertion written one base earlier spells the entry's allele nowhere at its span.
    let (outputs, lines) = added_outside_outputs("vcf", &inserted, ADDED_OUTSIDE_RECORD);
    let earlier =
        |written: &str| written.replace("\t7000000\tcell0\tG\tGT\t", "\t6999999\tcell0\tAG\tAGT\t");
    let lines: Vec<String> = lines.iter().map(|line| earlier(line)).collect();
    let e = entry_of("vcf", &earlier(&outputs[0]), &lines);
    let a = entry_of("vcf", &earlier(&outputs[1]), &lines);
    assert_eq!(
        (
            e.location.as_deref(),
            e.allele_location.as_deref(),
            e.allele.as_str()
        ),
        (None, None, "T")
    );
    assert!(!accepts(&e, &a, "HGVSc"));
}

/// Through the comparison itself, in every form of `RECORD_FORMS`, the class sets aside the
/// differing fields of an upstream row and of a downstream row, each row counted once, and
/// the row stays a difference without the class.
#[test]
fn a_row_of_bases_written_outside_counts_once_through_the_comparison_in_every_format() {
    let cells = [
        AddedOutsideCell {
            term: "upstream_gene_variant",
            strand: Some("-1"),
            candidate: ADDED_OUTSIDE_CANDIDATES[2],
            reference: OutsideReference::Fields("1", Some("0")),
        },
        AddedOutsideCell {
            term: "downstream_gene_variant",
            strand: Some("1"),
            candidate: ADDED_OUTSIDE_CANDIDATES[0],
            reference: OutsideReference::Nothing,
        },
    ];
    for (number, cell) in cells.iter().enumerate() {
        let rows = added_outside_rows(cell, number);
        for form in RECORD_FORMS {
            let (fmt, outputs, lines) = record_form_outputs(form, &rows, ADDED_OUTSIDE_RECORD);
            let expected = added_outside_parse(fmt, &outputs[0], &lines);
            let actual = added_outside_parse(fmt, &outputs[1], &lines);
            let with_the_class: (&[&str], &str, &[(&str, i64)]) =
                (&[ADDED_OUTSIDE], "green", &[(ADDED_OUTSIDE, 1)]);
            let without: (&[&str], &str, &[(&str, i64)]) = (&[], "red", &[]);
            for (classes, verdict_expected, applied) in [with_the_class, without] {
                let mut documented = documented_with(classes);
                let report = compare(&expected, &actual, &mut documented);
                assert_eq!(
                    (verdict(&report), counts(&documented.applied)),
                    (verdict_expected, count_map(applied)),
                    "{form} {cell:?} {classes:?}"
                );
            }
        }
    }
}

/// What a record adds to its reference: its whole allele, the surplus of its allele over
/// the bases it keeps, or nothing a `dup` describes.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Adds {
    Whole,
    Surplus,
    Nothing,
}

/// One record of the grid of records: what it is; the REF and ALT of its input line; the
/// positions its row's Location names past the first; the row's allele; and what it adds.
#[derive(Clone, Copy, Debug)]
struct OutsideRecord {
    label: &'static str,
    line: (&'static str, &'static str),
    further: i64,
    allele: &'static str,
    adds: Adds,
}

const OUTSIDE_RECORDS: [OutsideRecord; 9] = [
    OutsideRecord {
        label: "a pure insertion of five bases",
        line: ("G", "GTCTCA"),
        further: 1,
        allele: "TCTCA",
        adds: Adds::Whole,
    },
    OutsideRecord {
        label: "a pure insertion of two bases",
        line: ("G", "GTC"),
        further: 1,
        allele: "TC",
        adds: Adds::Whole,
    },
    OutsideRecord {
        label: "a pure insertion of one base",
        line: ("G", "GT"),
        further: 1,
        allele: "T",
        adds: Adds::Whole,
    },
    OutsideRecord {
        label: "an insertion kept with its shared first base",
        line: ("G", "GTCTC,C"),
        further: 0,
        allele: "GTCTC",
        adds: Adds::Surplus,
    },
    OutsideRecord {
        label: "two bases replaced by a longer sequence that keeps them",
        line: ("GT", "GTCTC,CT"),
        further: 1,
        allele: "GTCTC",
        adds: Adds::Surplus,
    },
    OutsideRecord {
        label: "two bases replaced by a longer sequence that keeps neither",
        line: ("GT", "CACAC"),
        further: 1,
        allele: "CACAC",
        adds: Adds::Nothing,
    },
    OutsideRecord {
        label: "two bases in place of two",
        line: ("AC", "GT"),
        further: 1,
        allele: "GT",
        adds: Adds::Nothing,
    },
    OutsideRecord {
        label: "a deletion spelt as a sequence, one of two bases",
        line: ("GT", "G,CT"),
        further: 1,
        allele: "G",
        adds: Adds::Nothing,
    },
    OutsideRecord {
        label: "a deletion spelt as a sequence, two of seven bases",
        line: ("GTCTCAG", "GTCTC,CTCTCAG"),
        further: 6,
        allele: "GTCTC",
        adds: Adds::Nothing,
    },
];

/// The bases a record's allele holds and its surplus over the positions of its Location.
fn held_and_surplus(record: &OutsideRecord) -> (i64, i64) {
    let held = record.allele.len() as i64;
    (held, held - (record.further + 1))
}

/// The numbers of bases a candidate `dup` counts on a record's row, in rising order: the
/// allele's length, its surplus over the positions of the Location, and the numbers one off
/// either.
fn outside_counts(record: &OutsideRecord) -> Vec<i64> {
    let (held, surplus) = held_and_surplus(record);
    let mut counts = vec![held, surplus];
    counts.extend([held - 1, held + 1, surplus - 1, surplus + 1]);
    counts.retain(|count| *count > 0);
    counts.sort_unstable();
    counts.dedup();
    counts
}

/// What the rule's statement says of a row whose entry carries its record's alleles: a pure
/// insertion's notation counts the allele's length, the notation of an allele that adds bases
/// to the ones it keeps counts those, and no count describes any other record.
fn set_aside_from_the_alleles(record: &OutsideRecord, count: i64) -> bool {
    let (held, surplus) = held_and_surplus(record);
    match record.adds {
        Adds::Whole => count == held,
        Adds::Surplus => count == surplus,
        Adds::Nothing => false,
    }
}

/// What the rule's statement says of a row under an identifier, which carries no alleles: a
/// Location of two neighbouring positions is an insertion's, whose notation counts the
/// allele's length, and an allele longer than its Location spans adds its surplus.
fn set_aside_from_the_location(record: &OutsideRecord, count: i64) -> bool {
    let (held, surplus) = held_and_surplus(record);
    (record.further == 1 && count == held) || (surplus > 0 && count == surplus)
}

/// The grid of records, a line a cell: the record; the bases the candidate's `dup` counts;
/// whether the class sets the row aside under an identifier; and whether it does in each
/// form that carries the record's alleles.
const OUTSIDE_RECORD_CELLS: [(&str, i64, bool, bool); 32] = [
    ("a pure insertion of five bases", 2, false, false),
    ("a pure insertion of five bases", 3, true, false),
    ("a pure insertion of five bases", 4, false, false),
    ("a pure insertion of five bases", 5, true, true),
    ("a pure insertion of five bases", 6, false, false),
    ("a pure insertion of two bases", 1, false, false),
    ("a pure insertion of two bases", 2, true, true),
    ("a pure insertion of two bases", 3, false, false),
    ("a pure insertion of one base", 1, true, true),
    ("a pure insertion of one base", 2, false, false),
    (
        "an insertion kept with its shared first base",
        3,
        false,
        false,
    ),
    (
        "an insertion kept with its shared first base",
        4,
        true,
        true,
    ),
    (
        "an insertion kept with its shared first base",
        5,
        false,
        false,
    ),
    (
        "an insertion kept with its shared first base",
        6,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps them",
        2,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps them",
        3,
        true,
        true,
    ),
    (
        "two bases replaced by a longer sequence that keeps them",
        4,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps them",
        5,
        true,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps them",
        6,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps neither",
        2,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps neither",
        3,
        true,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps neither",
        4,
        false,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps neither",
        5,
        true,
        false,
    ),
    (
        "two bases replaced by a longer sequence that keeps neither",
        6,
        false,
        false,
    ),
    ("two bases in place of two", 1, false, false),
    ("two bases in place of two", 2, true, false),
    ("two bases in place of two", 3, false, false),
    (
        "a deletion spelt as a sequence, one of two bases",
        1,
        true,
        false,
    ),
    (
        "a deletion spelt as a sequence, one of two bases",
        2,
        false,
        false,
    ),
    (
        "a deletion spelt as a sequence, two of seven bases",
        4,
        false,
        false,
    ),
    (
        "a deletion spelt as a sequence, two of seven bases",
        5,
        false,
        false,
    ),
    (
        "a deletion spelt as a sequence, two of seven bases",
        6,
        false,
        false,
    ),
];

/// The reference's and the candidate's default-format rows of a record under an identifier:
/// a row five bases from the transcript on which the candidate writes a `dup` of `count`
/// bases, the duplication of the transcript's last bases on a downstream row and one the 3'
/// rule carried twenty bases on an upstream row.
fn outside_record_rows(
    record: &OutsideRecord,
    term: &str,
    count: i64,
    number: usize,
) -> [String; 2] {
    let position = 8_000_000 + 100 * number as i64;
    let location = match record.further {
        0 => format!("2:{position}"),
        further => format!("2:{position}-{}", position + further),
    };
    let columns = format!(
        "cell{number}\t{location}\t{}\tENSG1\t{ADDED_OUTSIDE_FEATURE}\tTranscript\t{term}\t-\t-\t-\t-\t-\t-",
        record.allele
    );
    let extra = "IMPACT=MODIFIER;DISTANCE=5;STRAND=1";
    let (last, offset) = match term {
        "downstream_gene_variant" => (2601, ""),
        _ => (14, ";HGVS_OFFSET=20"),
    };
    let duplicated = match count {
        1 => format!("n.{last}dup"),
        count => format!("n.{}_{last}dup", last - count + 1),
    };
    [
        format!("{columns}\t{extra}"),
        format!(
            "{columns}\t{extra};HGVSc={ADDED_OUTSIDE_FEATURE}.1:{duplicated}{offset};HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0"
        ),
    ]
}

/// What the class says of the reference's entry `e` and the candidate's entry `a` of one
/// row: whether it sets aside every one of its fields the two write differently, which it
/// does of all or of none.
fn sets_the_row_aside(e: &Entry, a: &Entry) -> bool {
    let read: Vec<bool> = ADDED_OUTSIDE_FIELDS
        .into_iter()
        .filter(|field| field_of(Some(e), &[field]) != field_of(Some(a), &[field]))
        .map(|field| accepts(e, a, field))
        .collect();
    assert!(read.len() >= 2, "{e:?} {a:?}");
    assert!(read.iter().all(|taken| *taken == read[0]), "{e:?} {read:?}");
    read[0]
}

/// The grid of records, each cell as a downstream and as an upstream row in every form of
/// `RECORD_FORMS`. The forms that carry the record's alleles read every cell alike, as the
/// rule's statement reads the record: a `dup` of the whole allele on a pure insertion's row,
/// of the added bases on the row of an allele that keeps the bases it stands in for, and no
/// `dup` on the row of a replacement, a substitution or a deletion. A row under an identifier
/// is read from its Location and allele alone, and so is one that carries the bases its
/// allele stands in for (`Entry.reference`) read as the forms with alleles. The cells a row
/// under an identifier reads otherwise are the six a Location of two neighbouring positions
/// shares with an insertion or an allele longer than its span shares with added bases, and
/// in each the class sets the row aside.
#[test]
fn a_record_outside_reads_one_way_in_every_form_that_carries_its_alleles() {
    let grid: Vec<(&str, i64)> = OUTSIDE_RECORDS
        .iter()
        .flat_map(|record| {
            outside_counts(record)
                .into_iter()
                .map(|count| (record.label, count))
        })
        .collect();
    let table: Vec<(&str, i64)> = OUTSIDE_RECORD_CELLS
        .iter()
        .map(|(label, count, _, _)| (*label, *count))
        .collect();
    assert_eq!(grid, table);

    let mut number = 0;
    let mut limits = Vec::new();
    for (label, count, under_an_identifier, with_alleles) in OUTSIDE_RECORD_CELLS {
        let record = OUTSIDE_RECORDS
            .iter()
            .find(|record| record.label == label)
            .expect("a record of the grid");
        assert_eq!(
            (under_an_identifier, with_alleles),
            (
                set_aside_from_the_location(record, count),
                set_aside_from_the_alleles(record, count)
            ),
            "{label} {count}"
        );
        // The bases the allele stands in for: none on a pure insertion's row, the line's
        // REF on the row of a line whose alleles keep their first base.
        let replaced = match record.adds {
            Adds::Whole => "",
            _ => record.line.0,
        };
        for term in ["downstream_gene_variant", "upstream_gene_variant"] {
            number += 1;
            let rows = outside_record_rows(record, term, count, number);
            for form in RECORD_FORMS {
                let [e, a] = record_form_entries(form, &rows, record.line);
                let expected = if FORMS_WITH_ALLELES.contains(&form) {
                    assert!(e.indel.is_some() && e.insertion.is_some(), "{form} {label}");
                    with_alleles
                } else {
                    assert_eq!((e.indel, e.insertion, &e.reference), (None, None, &None));
                    under_an_identifier
                };
                assert_eq!(
                    sets_the_row_aside(&e, &a),
                    expected,
                    "{form} {label} {count} {term}"
                );
            }
            let [mut e, a] = record_form_entries("identifier", &rows, record.line);
            e.reference = Some(replaced.to_string());
            assert_eq!(
                sets_the_row_aside(&e, &a),
                with_alleles,
                "an identifier beside reference bases {label} {count} {term}"
            );
        }
        if under_an_identifier != with_alleles {
            assert!(under_an_identifier, "{label} {count}");
            limits.push((label, count));
        }
    }
    assert_eq!(number, 64);
    assert_eq!(
        limits,
        [
            ("a pure insertion of five bases", 3),
            ("two bases replaced by a longer sequence that keeps them", 5),
            (
                "two bases replaced by a longer sequence that keeps neither",
                3
            ),
            (
                "two bases replaced by a longer sequence that keeps neither",
                5
            ),
            ("two bases in place of two", 2),
            ("a deletion spelt as a sequence, one of two bases", 1),
        ]
    );
}

/// An upstream row whose candidate `dup` lies in an intron, where its coordinates count no
/// bases, carried as far as the row's DISTANCE. The forms that carry the record's alleles
/// read a base added to the two bases it keeps as an allele longer than its span, whose
/// HGVS_OFFSET may be the DISTANCE itself. A row under an identifier reads the same
/// Location as an insertion's, whose HGVS_OFFSET is larger, and stays a difference. A base
/// added to the one base it keeps is no insertion in any form, an inserted base is one in
/// every form, and an HGVS_OFFSET below the DISTANCE stays a difference in every form.
#[test]
fn a_dup_that_counts_no_bases_reads_the_distance_from_the_alleles_where_a_row_carries_them() {
    let cell = AddedOutsideCell {
        term: "upstream_gene_variant",
        strand: Some("1"),
        candidate: ADDED_OUTSIDE_CANDIDATES[3],
        reference: OutsideReference::Fields("0", Some("0")),
    };
    let inserted = added_outside_rows(&cell, 0);
    let rows_with = |changes: &[(&str, &str)]| -> [String; 2] {
        let mut rows = inserted.clone();
        for (old, new) in changes {
            rows = rows.map(|row| row.replace(old, new));
        }
        rows
    };
    let kept_two = ("\tT\t", "\tGTT\t");
    let kept_one = [("-7000001\t", "\t"), ("\tT\t", "\tGT\t")];
    let the_distance = ("HGVS_OFFSET=26", "HGVS_OFFSET=5");
    let less = ("HGVS_OFFSET=26", "HGVS_OFFSET=4");
    let shapes = [
        (
            "a base added to the two kept, carried the distance",
            rows_with(&[kept_two, the_distance]),
            ("GT", "GTT,CT"),
            false,
            true,
        ),
        (
            "a base added to the two kept, carried less",
            rows_with(&[kept_two, less]),
            ("GT", "GTT,CT"),
            false,
            false,
        ),
        (
            "a base added to the one kept, carried the distance",
            rows_with(&[kept_one[0], kept_one[1], the_distance]),
            ("G", "GT,C"),
            true,
            true,
        ),
        (
            "an inserted base carried the distance",
            rows_with(&[the_distance]),
            ("G", "GT"),
            false,
            false,
        ),
        (
            "an inserted base carried further",
            inserted.clone(),
            ("G", "GT"),
            true,
            true,
        ),
    ];
    for (label, rows, record, under_an_identifier, with_alleles) in &shapes {
        for form in RECORD_FORMS {
            let [e, a] = record_form_entries(form, rows, *record);
            let expected = if FORMS_WITH_ALLELES.contains(&form) {
                with_alleles
            } else {
                under_an_identifier
            };
            assert_eq!(sets_the_row_aside(&e, &a), *expected, "{form} {label}");
        }
    }
}

/// The record `chr16 15725140 AA A,AAA` on a forward-strand transcript whose
/// last exon ends inside the run: the reference's `c.*890dup` and the candidate's `c.*898dup`
/// carried 8 bases, read as VCF entries, which carry their Location as `allele_location` alone.
#[test]
fn added_bases_carried_along_a_run_reads_a_vcf_entry_at_its_records_location() {
    let line = "chr16\t15725140\tend_run\tAA\tA,AAA\t.\tPASS\t.".to_string();
    let output = |hgvsc: &str, offset: &str| -> Entry {
        let header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|HGVSc|STRAND|HGVS_OFFSET\">\n";
        let text = format!(
            "{header}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}CSQ=AA|3_prime_UTR_variant|Transcript|ENST00000911229|ENST00000911229.1:{hgvsc}|1|{offset}\n",
            &line[..line.len() - 1]
        );
        parse_vcf(&text, std::slice::from_ref(&line))
            .expect("the text parses")
            .pop()
            .expect("the CSQ entry follows the record's")
    };
    let (reference, candidate) = (output("c.*890dup", ""), output("c.*898dup", "8"));
    assert_eq!(reference.location, None);
    for field in ["HGVSc", "HGVS_OFFSET"] {
        assert!(
            added_bases::added_bases_left_at_the_span_end(&reference, Some(&candidate), field)
                .expect("the rule reads the pair"),
            "{field}"
        );
    }
}
