// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Made-up grids of rows compared in the tab format under one class: the comparison sets a
//! field aside exactly where the class's statement does.

mod compare_common;

use compare_common::*;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::pyjson::Value;

const HGVS_INTRON_OFFSET_FIELDS: [&str; 2] = ["HGVS_IntronStartOffset", "HGVS_IntronEndOffset"];

const OFFSET_GRID_TAB_COLUMNS: [&str; 8] = [
    "IMPACT",
    "DISTANCE",
    "STRAND",
    "BAM_EDIT",
    "HGVSc",
    "HGVS_OFFSET",
    "HGVS_IntronStartOffset",
    "HGVS_IntronEndOffset",
];

/// The value the pairs of an Extra text hold for `name`, the last where it is named twice.
fn value_of<'a>(held: &[(&'a str, &'a str)], name: &str) -> Option<&'a str> {
    held.iter()
        .rev()
        .find(|(key, _)| *key == name)
        .map(|(_, value)| *value)
}

/// The columns of a default-format row before its Extra text, and the pairs of that text.
type SplitRow<'a> = (Vec<&'a str>, Vec<(&'a str, &'a str)>);

/// Default-format rows of one run as the tab format writes them, one output a row: each
/// Extra field that a row carries is a column of every output.
fn offset_grid_tab(rows: &[&str]) -> Vec<String> {
    let split: Vec<SplitRow> = rows
        .iter()
        .map(|row| {
            let mut columns: Vec<&str> = row.split('\t').collect();
            let extra = columns.pop().unwrap_or_default();
            let held = extra
                .split(';')
                .map(|item| item.split_once('=').expect("a pair of the Extra text"))
                .collect();
            (columns, held)
        })
        .collect();
    let names: Vec<&str> = OFFSET_GRID_TAB_COLUMNS
        .into_iter()
        .filter(|name| split.iter().any(|(_, held)| value_of(held, name).is_some()))
        .collect();
    let named: String = names.iter().map(|name| format!("\t{name}")).collect();
    let header = DEFAULT_HEADER.replace("\tExtra", &named);
    split
        .iter()
        .map(|(columns, held)| {
            let mut cells = columns.clone();
            cells.extend(names.iter().map(|name| value_of(held, name).unwrap_or("-")));
            format!("{header}\n{}\n", cells.join("\t"))
        })
        .collect()
}

/// Whether the comparison of the two rows' tab-format outputs under `klass` sets the field
/// aside: its column counts the class and no difference.
fn set_aside_by_the_comparison(klass: &str, reference: &str, candidate: &str, field: &str) -> bool {
    let outputs = offset_grid_tab(&[reference, candidate]);
    let report = compare(
        &parse_tab(&outputs[0]).expect("the reference parses"),
        &parse_tab(&outputs[1]).expect("the candidate parses"),
        &mut documented_with(&[klass]),
    );
    let column = at(&report, &["columns", field]);
    let counted = column
        .get("masked")
        .and_then(Value::as_dict)
        .is_some_and(|masked| masked.len() == 1 && masked.get(klass) == Some(&Value::Int(1)));
    counted && column.get("differing") == Some(&Value::Int(0))
}

const ABUTTING_GRID_PAIRS: [(&str, &str); 8] = [
    ("0", "0"),
    ("1", "0"),
    ("0", "-1"),
    ("-1", "0"),
    ("0", "1"),
    ("1", "-1"),
    ("-1", "-1"),
    ("2", "0"),
];
const ABUTTING_GRID_TERMS: [&str; 3] = [
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intron_variant",
];
const ABUTTING_GRID_NOTATIONS: [Option<(&str, &str)>; 5] = [
    None,
    Some(("reference", "HGVSc=ENST00000000041.1:n.1_2insA")),
    Some(("candidate", "HGVSc=ENST00000000041.1:n.1_2insA")),
    Some(("reference", "HGVS_OFFSET=2")),
    Some(("candidate", "HGVS_OFFSET=2")),
];

/// One cell of the made-up grid of `hgvs_intron_offset_abutting_insertion`: the plugin's
/// pair, the row's term, its strand, its DISTANCE, an HGVSc or an HGVS_OFFSET on either
/// side, the candidate's two fields blank or filled, and an insertion's Location or a
/// one-base one.
#[derive(Clone, Copy)]
struct AbuttingCell {
    pair: (&'static str, &'static str),
    term: &'static str,
    strand: Option<&'static str>,
    distance: &'static str,
    notation: Option<(&'static str, &'static str)>,
    filled: bool,
    insertion: bool,
}

/// Every combination of the grid's values, the last varying fastest.
fn abutting_grid() -> Vec<AbuttingCell> {
    let mut cells = Vec::new();
    for pair in ABUTTING_GRID_PAIRS {
        for term in ABUTTING_GRID_TERMS {
            for strand in [Some("1"), Some("-1"), None] {
                for distance in ["0", "5"] {
                    for notation in ABUTTING_GRID_NOTATIONS {
                        for filled in [false, true] {
                            for insertion in [true, false] {
                                cells.push(AbuttingCell {
                                    pair,
                                    term,
                                    strand,
                                    distance,
                                    notation,
                                    filled,
                                    insertion,
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a cell: the row is an insertion's outside the
/// transcript at a DISTANCE of 0, no HGVSc and no HGVS_OFFSET stands on either side, the
/// candidate writes neither field, and the pair is 0 and 0, or 1 and 0 on a downstream row
/// whose strand is not forward, or 0 and -1 on an upstream row whose strand is not reverse.
fn abutting_grid_set_aside(cell: &AbuttingCell) -> bool {
    if !cell.insertion
        || cell.term == "intron_variant"
        || cell.distance != "0"
        || cell.notation.is_some()
        || cell.filled
    {
        return false;
    }
    match cell.pair {
        ("0", "0") => true,
        ("1", "0") => cell.term == "downstream_gene_variant" && cell.strand != Some("1"),
        ("0", "-1") => cell.term == "upstream_gene_variant" && cell.strand != Some("-1"),
        _ => false,
    }
}

/// The reference's and the candidate's default-format rows of a cell, under a Location
/// and a Feature of the cell's own.
fn abutting_grid_rows(cell: &AbuttingCell, number: usize) -> (String, String) {
    let (start, end) = cell.pair;
    let position = 1_000_000 + 10 * number;
    let location = if cell.insertion {
        format!("2:{position}-{}", position + 1)
    } else {
        format!("2:{position}")
    };
    let mut extra = format!("IMPACT=MODIFIER;DISTANCE={}", cell.distance);
    if let Some(strand) = cell.strand {
        extra.push_str(&format!(";STRAND={strand}"));
    }
    let mut reference = extra.clone();
    let mut candidate = extra;
    match cell.notation {
        Some(("reference", held)) => reference.push_str(&format!(";{held}")),
        Some((_, held)) => candidate.push_str(&format!(";{held}")),
        None => {}
    }
    reference.push_str(&format!(
        ";HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}"
    ));
    if cell.filled {
        candidate.push_str(";HGVS_IntronEndOffset=7;HGVS_IntronStartOffset=7");
    }
    let row = |extra: &str| {
        default_row(
            &format!("cell{number}"),
            &location,
            "T",
            &format!("ENST{number:011}"),
            cell.term,
            extra,
        )
    };
    (row(&reference), row(&candidate))
}

/// The format comparison run on each cell's two outputs with the class sets aside the
/// fields the rule's statement names, and no other.
#[test]
fn the_abutting_insertion_grid_is_set_aside_where_the_statement_says() {
    let klass = "hgvs_intron_offset_abutting_insertion";
    let cells = abutting_grid();
    assert_eq!(cells.len(), 2880);
    assert_eq!(
        cells
            .iter()
            .filter(|cell| abutting_grid_set_aside(cell))
            .count(),
        10
    );
    let mut wrong: Vec<(usize, &str)> = Vec::new();
    for (number, cell) in cells.iter().enumerate() {
        let (reference, candidate) = abutting_grid_rows(cell, number);
        for field in HGVS_INTRON_OFFSET_FIELDS {
            let compared = set_aside_by_the_comparison(klass, &reference, &candidate, field);
            if compared != abutting_grid_set_aside(cell) {
                wrong.push((number, field));
            }
        }
    }
    assert!(wrong.is_empty(), "{wrong:?}");
}

const CLIPPED: &str = "hgvs_clipped_reference_against_whole_alternate";

/// The made-up rows of the plugin's fields beside the class's notation: the kind of the row,
/// the reference's notation after the sequence name, and the two values the notation routine
/// stores for it, the start's and the end's.
const CLIPPED_OFFSETS_ROWS: [(&str, &str, (&str, &str)); 15] = [
    ("ins", "c.*234+1_*234insT", ("1", "0")),
    ("ins", "c.76-1_76insT", ("0", "-1")),
    ("ins", "c.76_76+1insT", ("1", "0")),
    ("ins", "c.75+50_76-50insT", ("-50", "50")),
    ("dup", "c.*234+1_*234+2dup", ("1", "2")),
    ("dup", "c.76-2_76-1dup", ("-2", "-1")),
    ("dup", "c.76+1dup", ("1", "1")),
    ("closed", "c.*234+1_*234insT", ("1", "0")),
    ("delins", "c.*234+1_*234delinsTG", ("1", "0")),
    ("deletion", "c.*234+1_*234insT", ("1", "0")),
    ("unprinted", "c.76_75insT", ("0", "0")),
    ("shorter del", "c.*232_*234del", ("0", "0")),
    ("shorter delins", "c.*233+2_*234delinsT", ("2", "0")),
    ("shorter substitution", "c.*234G>T", ("0", "0")),
    ("shorter inversion", "c.*233_*234inv", ("0", "0")),
];
const CLIPPED_OFFSETS_VALUES: [&str; 4] = ["read", "turned", "zeros", "start alone"];
const CLIPPED_OFFSETS_FEATURE: &str = "ENST00000000051";

/// One cell of the made-up grid of the plugin's fields under the class: the field asked;
/// the row; the reference's two values (the ones the routine stores, the same two the other
/// way round, 0 and 0, the start's alone); the candidate's fields (blank, the asked one
/// filled, the other one filled); the candidate's HGVSc (none, the reference's own, a
/// deletion's).
#[derive(Clone, Copy)]
struct ClippedOffsetsCell {
    field: &'static str,
    row: (&'static str, &'static str, (&'static str, &'static str)),
    values: &'static str,
    filled: Option<&'static str>,
    candidate_notation: Option<&'static str>,
}

fn clipped_offsets_grid() -> Vec<ClippedOffsetsCell> {
    let mut cells = Vec::new();
    for field in HGVS_INTRON_OFFSET_FIELDS {
        for row in CLIPPED_OFFSETS_ROWS {
            for values in CLIPPED_OFFSETS_VALUES {
                for filled in [None, Some("asked"), Some("other")] {
                    for candidate_notation in [None, Some("same"), Some("del")] {
                        cells.push(ClippedOffsetsCell {
                            field,
                            row,
                            values,
                            filled,
                            candidate_notation,
                        });
                    }
                }
            }
        }
    }
    cells
}

/// The start and end values the reference's row of a cell writes.
fn clipped_offsets_written(cell: &ClippedOffsetsCell) -> (String, String) {
    let (start, end) = cell.row.2;
    match cell.values {
        "read" => (start.to_string(), end.to_string()),
        "turned" => (end.to_string(), start.to_string()),
        "zeros" => ("0".to_string(), "0".to_string()),
        _ => {
            let end: i64 = end.parse().expect("a whole number");
            (start.to_string(), (end + 5).to_string())
        }
    }
}

/// What the rule's statement says of a cell: the row's HGVSc is one the class takes
/// against an empty candidate, the candidate prints no HGVSc and writes neither field, and
/// the reference's two values are the ones the routine stores for the printed notation.
fn clipped_offsets_set_aside(cell: &ClippedOffsetsCell) -> bool {
    let (kind, _, stored) = cell.row;
    let taken = [
        "ins",
        "dup",
        "shorter del",
        "shorter delins",
        "shorter substitution",
    ]
    .contains(&kind);
    let written = clipped_offsets_written(cell);
    taken
        && cell.filled.is_none()
        && cell.candidate_notation.is_none()
        && (written.0.as_str(), written.1.as_str()) == stored
}

/// The reference's and the candidate's default-format rows of a cell, under a Location of
/// the cell's own.
fn clipped_offsets_rows(cell: &ClippedOffsetsCell, number: usize) -> (String, String) {
    let (kind, notation, _) = cell.row;
    let (start, end) = clipped_offsets_written(cell);
    let position = 3_000_000 + 10 * number;
    let allele = if kind == "deletion" {
        "-"
    } else if kind.starts_with("shorter") {
        "TG"
    } else {
        "TGCAT"
    };
    let cdna = if kind == "closed" { "553-557" } else { "553-?" };
    let columns = format!(
        "cell{number}\t21:{position}-{}\t{allele}\tENSG1\t{CLIPPED_OFFSETS_FEATURE}\tTranscript\t3_prime_UTR_variant\t{cdna}\t-\t-\t-\t-\t-",
        position + 4
    );
    let reference = format!(
        "IMPACT=MODIFIER;STRAND=-1;HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{notation};HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}"
    );
    let mut candidate = "IMPACT=MODIFIER;STRAND=-1".to_string();
    if let Some(candidate_notation) = cell.candidate_notation {
        let written = if candidate_notation == "same" {
            notation
        } else {
            "c.*230_*233del"
        };
        candidate.push_str(&format!(";HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{written}"));
    }
    if let Some(filled) = cell.filled {
        let other = HGVS_INTRON_OFFSET_FIELDS
            .into_iter()
            .find(|f| *f != cell.field)
            .expect("the plugin writes two fields");
        let name = if filled == "asked" { cell.field } else { other };
        candidate.push_str(&format!(";{name}=99"));
    }
    (
        format!("{columns}\t{reference}"),
        format!("{columns}\t{candidate}"),
    )
}

/// The format comparison run on each cell's two outputs with the class sets aside the
/// fields of the plugin the rule's statement names, and no other.
#[test]
fn the_clipped_reference_grid_is_set_aside_where_the_statement_says() {
    let cells = clipped_offsets_grid();
    assert_eq!(cells.len(), 1080);
    let mut wrong: Vec<usize> = Vec::new();
    for (number, cell) in cells.iter().enumerate() {
        let (reference, candidate) = clipped_offsets_rows(cell, number);
        let compared = set_aside_by_the_comparison(CLIPPED, &reference, &candidate, cell.field);
        if compared != clipped_offsets_set_aside(cell) {
            wrong.push(number);
        }
    }
    assert!(wrong.is_empty(), "{wrong:?}");
}

/// The made-up rows of the plugin's fields on the candidate's side: the row's STRAND and
/// cDNA position. The allele is 20 bases over a Location of 10, so the bases it adds are 10.
const CLIPPED_CANDIDATE_SIDES: [(Option<&str>, &str); 6] = [
    (Some("-1"), "?-7"),
    (Some("1"), "7-?"),
    (Some("1"), "?-7"),
    (Some("-1"), "7-?"),
    (Some("-1"), "3-7"),
    (None, "?-7"),
];
const CLIPPED_CANDIDATE_ALLELE: &str = "ACGTACGTAGACGTACGTAG";
/// What the candidate's notation is, and the notation after the sequence name.
const CLIPPED_CANDIDATE_NOTATIONS: [(&str, &str); 8] = [
    ("an ins between cDNA 7 and 8", "n.7_8insACGTACGTAG"),
    ("an ins between cDNA 8 and 9", "n.8_9insACGTACGTAG"),
    ("an ins between cDNA 5 and 6", "n.5_6insACGTACGTAG"),
    ("an ins under coding numbering", "c.-20_-19insACGTACGTAG"),
    ("an ins of fewer bases", "n.7_8insACGTACGTA"),
    ("an ins at intron coordinates", "n.7+1_7+2insACGTACGTAG"),
    ("a del", "n.1del"),
    ("a dup", "n.7dup"),
];
const CLIPPED_CANDIDATE_FIELDS: [(&str, Option<&str>); 3] =
    [("0", Some("0")), ("0", Some("1")), ("0", None)];

/// One cell of the made-up grid of the two fields on the candidate's side alone: the field
/// asked; the row's strand and cDNA position; the candidate's notation against a reference
/// that prints none; the candidate's two fields; an HGVS_OFFSET on neither row or on one.
#[derive(Clone, Copy)]
struct ClippedCandidateCell {
    field: &'static str,
    side: (Option<&'static str>, &'static str),
    notation: (&'static str, &'static str),
    fields: (&'static str, Option<&'static str>),
    offset: Option<&'static str>,
}

fn clipped_candidate_grid() -> Vec<ClippedCandidateCell> {
    let mut cells = Vec::new();
    for field in HGVS_INTRON_OFFSET_FIELDS {
        for side in CLIPPED_CANDIDATE_SIDES {
            for notation in CLIPPED_CANDIDATE_NOTATIONS {
                for fields in CLIPPED_CANDIDATE_FIELDS {
                    for offset in [None, Some("reference"), Some("candidate")] {
                        cells.push(ClippedCandidateCell {
                            field,
                            side,
                            notation,
                            fields,
                            offset,
                        });
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a cell: the row is a reverse-strand transcript's
/// and its cDNA position is open at its start; the candidate's notation is an `ins` of
/// the ten bases the allele adds, between two exon coordinates that stand at or before the
/// position's closed end (under `n.` numbering; under `c.` their place is not read);
/// the candidate's two fields are 0 and 0, the intron parts of those coordinates; and no
/// HGVS_OFFSET stands on either row.
fn clipped_candidate_set_aside(cell: &ClippedCandidateCell) -> bool {
    let inside = [
        "an ins between cDNA 7 and 8",
        "an ins between cDNA 5 and 6",
        "an ins under coding numbering",
    ];
    cell.side == (Some("-1"), "?-7")
        && inside.contains(&cell.notation.0)
        && cell.fields == ("0", Some("0"))
        && cell.offset.is_none()
}

fn clipped_candidate_rows(cell: &ClippedCandidateCell, number: usize) -> (String, String) {
    let (strand, position) = cell.side;
    let (start, end) = cell.fields;
    let place = 5_000_000 + 20 * number;
    let columns = format!(
        "cell{number}\t21:{place}-{}\t{CLIPPED_CANDIDATE_ALLELE}\tENSG1\t{CLIPPED_OFFSETS_FEATURE}\tTranscript\tnon_coding_transcript_exon_variant\t{position}\t-\t-\t-\t-\t-",
        place + 9
    );
    let shared = match strand {
        Some(strand) => format!("IMPACT=MODIFIER;STRAND={strand}"),
        None => "IMPACT=MODIFIER".to_string(),
    };
    let mut reference = shared.clone();
    let mut candidate = format!(
        "{shared};HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{}",
        cell.notation.1
    );
    match cell.offset {
        Some("reference") => reference.push_str(";HGVS_OFFSET=2"),
        Some(_) => candidate.push_str(";HGVS_OFFSET=2"),
        None => {}
    }
    candidate.push_str(&format!(";HGVS_IntronStartOffset={start}"));
    if let Some(end) = end {
        candidate.push_str(&format!(";HGVS_IntronEndOffset={end}"));
    }
    (
        format!("{columns}\t{reference}"),
        format!("{columns}\t{candidate}"),
    )
}

/// The made-up rows of the plugin's fields on both sides: the reference's two fields beside
/// its `ins` of the outside bases before the first base, the candidate's notation and its
/// two fields. The allele is 14 bases over a Location of 7.
const CLIPPED_BOTH_REFERENCE: &str = "c.-529-1_-529insACGTACGTACGT";
const CLIPPED_BOTH_REFERENCE_FIELDS: [(&str, &str); 3] = [("0", "-1"), ("-1", "0"), ("0", "0")];
/// What the candidate writes, its notation and its two fields.
const CLIPPED_BOTH_CANDIDATES: [(&str, &str, (&str, &str)); 4] = [
    (
        "an ins within the span's reach",
        "c.-528_-527insACGTACG",
        ("0", "0"),
    ),
    (
        "the same ins, another pair of fields",
        "c.-528_-527insACGTACG",
        ("0", "1"),
    ),
    (
        "an ins past the span's reach",
        "c.-500_-499insACGTACG",
        ("0", "0"),
    ),
    ("a del", "c.-528del", ("0", "0")),
];

/// One cell of the made-up grid of the two fields on both sides: the row's strand and cDNA
/// position; the reference's two fields; the candidate's notation and fields; an
/// HGVS_OFFSET on neither row or on the reference's. The field asked is the end's.
#[derive(Clone, Copy)]
struct ClippedBothCell {
    side: (Option<&'static str>, &'static str),
    fields: (&'static str, &'static str),
    candidate: (&'static str, &'static str, (&'static str, &'static str)),
    offset: bool,
}

fn clipped_both_grid() -> Vec<ClippedBothCell> {
    let mut cells = Vec::new();
    for side in CLIPPED_CANDIDATE_SIDES {
        for fields in CLIPPED_BOTH_REFERENCE_FIELDS {
            for candidate in CLIPPED_BOTH_CANDIDATES {
                for offset in [false, true] {
                    cells.push(ClippedBothCell {
                        side,
                        fields,
                        candidate,
                        offset,
                    });
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a cell: the cDNA position has an open side, the two
/// notations are a pair the class takes, an `ins` against an `ins` within the span's
/// reach, each side's two fields are the intron parts of its own notation, and no
/// HGVS_OFFSET stands on either row.
fn clipped_both_set_aside(cell: &ClippedBothCell) -> bool {
    cell.side.1.contains('?')
        && cell.fields == ("0", "-1")
        && cell.candidate.0 == "an ins within the span's reach"
        && !cell.offset
}

fn clipped_both_rows(cell: &ClippedBothCell, number: usize) -> (String, String) {
    let (strand, position) = cell.side;
    let (start, end) = cell.fields;
    let (_, notation, (candidate_start, candidate_end)) = cell.candidate;
    let place = 7_000_000 + 20 * number;
    let columns = format!(
        "both{number}\t21:{place}-{}\tACGTACGACGTACG\tENSG1\t{CLIPPED_OFFSETS_FEATURE}\tTranscript\t5_prime_UTR_variant\t{position}\t-\t-\t-\t-\t-",
        place + 6
    );
    let shared = match strand {
        Some(strand) => format!("IMPACT=MODIFIER;STRAND={strand}"),
        None => "IMPACT=MODIFIER".to_string(),
    };
    let offset = if cell.offset { ";HGVS_OFFSET=2" } else { "" };
    let reference = format!(
        "{shared};HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{CLIPPED_BOTH_REFERENCE}{offset};HGVS_IntronStartOffset={start};HGVS_IntronEndOffset={end}"
    );
    let candidate = format!(
        "{shared};HGVSc={CLIPPED_OFFSETS_FEATURE}.1:{notation};HGVS_IntronStartOffset={candidate_start};HGVS_IntronEndOffset={candidate_end}"
    );
    (
        format!("{columns}\t{reference}"),
        format!("{columns}\t{candidate}"),
    )
}

/// Each field of the two grids that differs, with its rows, and whether the rule's
/// statement sets it aside. On the grid of both sides the end's field alone differs.
fn clipped_candidate_side_asks() -> Vec<(String, String, &'static str, bool)> {
    let mut asks = Vec::new();
    for (number, cell) in clipped_candidate_grid().iter().enumerate() {
        if cell.field == "HGVS_IntronEndOffset" && cell.fields.1.is_none() {
            // The candidate writes no end field, so none differs.
            continue;
        }
        let (reference, candidate) = clipped_candidate_rows(cell, number);
        asks.push((
            reference,
            candidate,
            cell.field,
            clipped_candidate_set_aside(cell),
        ));
    }
    for (number, cell) in clipped_both_grid().iter().enumerate() {
        if cell.fields.1 == cell.candidate.2 .1 {
            // The two end fields agree.
            continue;
        }
        let (reference, candidate) = clipped_both_rows(cell, number);
        asks.push((
            reference,
            candidate,
            "HGVS_IntronEndOffset",
            clipped_both_set_aside(cell),
        ));
    }
    asks
}

/// The fields asked of a grid that the comparison under `klass` sets aside otherwise than
/// the statement does, by their place among the asks.
fn asks_decided_otherwise(klass: &str, asks: &[(String, String, &str, bool)]) -> Vec<usize> {
    (0..asks.len())
        .filter(|at| {
            let (reference, candidate, field, expected) = &asks[*at];
            set_aside_by_the_comparison(klass, reference, candidate, field) != *expected
        })
        .collect()
}

/// The format comparison run on each cell's two outputs with the class sets aside the
/// fields the candidate alone writes, and those both sides write, where the rule's
/// statement does, and no other.
#[test]
fn the_clipped_reference_candidate_side_grids_are_set_aside_where_the_statement_says() {
    let asks = clipped_candidate_side_asks();
    assert_eq!(asks.len(), 792);
    let wrong = asks_decided_otherwise(CLIPPED, &asks);
    assert!(wrong.is_empty(), "{wrong:?}");
}

const FIRST_BASE_ALLELE: &str = "GATTCACGATTCAC";

/// The made-up row of the plugin's end offset before a forward-strand transcript's first
/// base: an allele of 14 bases over a Location of 7 whose last 2 bases are the transcript's
/// first 2; the reference's `ins` of the allele's first 12 bases from one base before the
/// first base, beside the pair 0 and -1; the candidate's duplication of the 7 bases the
/// allele adds, ending 17 bases 3' of the span's end, beside an HGVS_OFFSET of 17 and the
/// pair 0 and 0.
#[derive(Clone)]
struct FirstBaseRow {
    allele: String,
    span: usize,
    strand: Option<&'static str>,
    position: &'static str,
    edit: Option<&'static str>,
    reference: String,
    reference_fields: (&'static str, &'static str),
    reference_offset: Option<&'static str>,
    candidate: String,
    candidate_fields: (&'static str, &'static str),
    candidate_offset: Option<&'static str>,
}

impl Default for FirstBaseRow {
    fn default() -> Self {
        FirstBaseRow {
            allele: FIRST_BASE_ALLELE.to_string(),
            span: 7,
            strand: Some("1"),
            position: "?-2",
            edit: None,
            reference: format!("c.-529-1_-529ins{}", &FIRST_BASE_ALLELE[..12]),
            reference_fields: ("0", "-1"),
            reference_offset: None,
            candidate: "c.-517_-511dup".to_string(),
            candidate_fields: ("0", "0"),
            candidate_offset: Some("17"),
        }
    }
}

/// The row and the rows that depart from it in the part each names, with whether the end
/// offset is set aside.
fn first_base_rows_named() -> Vec<(&'static str, FirstBaseRow, bool)> {
    let row = FirstBaseRow::default;
    let reference = |notation: String| FirstBaseRow {
        reference: notation,
        ..row()
    };
    let candidate = |notation: &str| FirstBaseRow {
        candidate: notation.to_string(),
        ..row()
    };
    let inserted = |coordinates: &str, bases: &str| format!("{coordinates}ins{bases}");
    let allele = FIRST_BASE_ALLELE;
    vec![
        ("the row", row(), true),
        (
            "one added base that the rule does not carry",
            FirstBaseRow {
                allele: allele[..8].to_string(),
                reference: inserted("c.-529-1_-529", &allele[..6]),
                candidate: "c.-528dup".to_string(),
                candidate_offset: None,
                ..row()
            },
            true,
        ),
        (
            "an allele as long as its span",
            FirstBaseRow { span: 14, ..row() },
            false,
        ),
        (
            "a reverse-strand transcript's row",
            FirstBaseRow {
                strand: Some("-1"),
                ..row()
            },
            false,
        ),
        (
            "a row without STRAND",
            FirstBaseRow {
                strand: None,
                ..row()
            },
            false,
        ),
        (
            "a position open at its end",
            FirstBaseRow {
                position: "2-?",
                ..row()
            },
            false,
        ),
        (
            "a closed position",
            FirstBaseRow {
                position: "1-2",
                ..row()
            },
            false,
        ),
        // A span that begins in an intron and ends on an exon's seventh base is not clamped:
        // its closed end is the span's length.
        (
            "a closed end as large as the span",
            FirstBaseRow {
                position: "?-7",
                reference: inserted("c.-529-1_-529", &allele[..7]),
                candidate_offset: Some("12"),
                ..row()
            },
            false,
        ),
        (
            "a reference dup",
            FirstBaseRow {
                reference: "c.-529-1dup".to_string(),
                reference_fields: ("-1", "-1"),
                ..row()
            },
            false,
        ),
        (
            "a reference ins between two exonic bases",
            reference(inserted("c.-529_-528", &allele[..12])),
            false,
        ),
        (
            "an intron part on another coordinate",
            reference(inserted("c.-530-1_-529", &allele[..12])),
            false,
        ),
        (
            "an intron part of -2",
            reference(inserted("c.-529-2_-529", &allele[..12])),
            false,
        ),
        // The `ins` the clip leaves at a true intron's last base, nothing clamped, writes the
        // bases the allele adds and no base of the span.
        (
            "a clip-made ins at a true intron",
            FirstBaseRow {
                reference: inserted("c.88-1_88", &allele[5..12]),
                candidate: "c.100_106dup".to_string(),
                ..row()
            },
            false,
        ),
        (
            "a reference ins of the added bases alone",
            reference(inserted("c.-529-1_-529", &allele[5..12])),
            false,
        ),
        // A clamped span that reaches across a whole intron hands the clip the intron's bases too.
        (
            "a reference ins of fewer bases",
            reference(inserted("c.-529-1_-529", &allele[..11])),
            false,
        ),
        (
            "a reference ins of the allele less its first 2 bases",
            reference(inserted("c.-529-1_-529", &allele[2..])),
            false,
        ),
        (
            "a reference HGVS_OFFSET",
            FirstBaseRow {
                reference_offset: Some("2"),
                ..row()
            },
            false,
        ),
        (
            "the reference's pair -1 and -1",
            FirstBaseRow {
                reference_fields: ("-1", "-1"),
                ..row()
            },
            false,
        ),
        (
            "the reference's pair 0 and -2",
            FirstBaseRow {
                reference_fields: ("0", "-2"),
                ..row()
            },
            false,
        ),
        // A plain insertion at an exon's first base: a Location of two positions, a closed cDNA
        // position, and the pair -1 and 0, the base before the point read first.
        (
            "a plain insertion at an exon's first base",
            FirstBaseRow {
                allele: allele[..7].to_string(),
                span: 2,
                position: "900-901",
                reference: inserted("c.784-1_784", &allele[..7]),
                reference_fields: ("-1", "0"),
                candidate: "c.790_796dup".to_string(),
                candidate_offset: Some("13"),
                ..row()
            },
            false,
        ),
        (
            "a candidate dup under the other numbering",
            candidate("n.13_19dup"),
            false,
        ),
        (
            "a candidate dup of another transcript",
            candidate("ENST00000000052.1:c.-517_-511dup"),
            false,
        ),
        (
            "a candidate dup at an intron coordinate",
            candidate("c.-517_-511+1dup"),
            false,
        ),
        (
            "a candidate dup of 6 bases",
            candidate("c.-516_-511dup"),
            false,
        ),
        (
            "a candidate dup of 8 bases",
            candidate("c.-518_-511dup"),
            false,
        ),
        ("a candidate inversion", candidate("c.-517_-511inv"), false),
        (
            "a negative HGVS_OFFSET",
            FirstBaseRow {
                candidate_offset: Some("-17"),
                ..row()
            },
            false,
        ),
        (
            "the candidate's pair 0 and 1",
            FirstBaseRow {
                candidate_fields: ("0", "1"),
                ..row()
            },
            false,
        ),
        (
            "the candidate's pair 1 and 0",
            FirstBaseRow {
                candidate_fields: ("1", "0"),
                ..row()
            },
            false,
        ),
        (
            "an HGVS_OFFSET one base short",
            FirstBaseRow {
                candidate_offset: Some("16"),
                ..row()
            },
            false,
        ),
        (
            "no HGVS_OFFSET beside a duplication 17 bases on",
            FirstBaseRow {
                candidate_offset: None,
                ..row()
            },
            false,
        ),
        (
            "a duplication one base further 3'",
            candidate("c.-516_-510dup"),
            false,
        ),
        (
            "a duplication across the start codon",
            FirstBaseRow {
                candidate: "c.5_11dup".to_string(),
                candidate_offset: Some("538"),
                ..row()
            },
            false,
        ),
        (
            "a row carrying BAM_EDIT OK",
            FirstBaseRow {
                edit: Some("OK"),
                ..row()
            },
            false,
        ),
        (
            "a row carrying BAM_EDIT FAILED",
            FirstBaseRow {
                edit: Some("FAILED"),
                ..row()
            },
            false,
        ),
    ]
}

/// The reference's and the candidate's default-format rows of one made-up row, under a
/// Location of the row's own.
fn first_base_rows(row: &FirstBaseRow, number: usize) -> (String, String) {
    let place = 9_000_000 + 40 * number;
    let columns = format!(
        "first{number}\t21:{place}-{}\t{}\tENSG1\t{CLIPPED_OFFSETS_FEATURE}\tTranscript\t5_prime_UTR_variant\t{}\t-\t-\t-\t-\t-",
        place + row.span - 1,
        row.allele,
        row.position
    );
    let mut shared = "IMPACT=MODIFIER".to_string();
    if let Some(strand) = row.strand {
        shared.push_str(&format!(";STRAND={strand}"));
    }
    if let Some(edit) = row.edit {
        shared.push_str(&format!(";BAM_EDIT={edit}"));
    }
    let extra = |notation: &str, offset: Option<&str>, held: (&str, &str)| {
        let named = if notation.contains(':') {
            notation.to_string()
        } else {
            format!("{CLIPPED_OFFSETS_FEATURE}.1:{notation}")
        };
        let offset = offset.map_or(String::new(), |offset| format!(";HGVS_OFFSET={offset}"));
        format!(
            "{columns}\t{shared};HGVSc={named}{offset};HGVS_IntronStartOffset={};HGVS_IntronEndOffset={}",
            held.0, held.1
        )
    };
    (
        extra(&row.reference, row.reference_offset, row.reference_fields),
        extra(&row.candidate, row.candidate_offset, row.candidate_fields),
    )
}

/// The end offsets of the made-up first-base rows whose two sides differ, with their rows
/// and whether the rule's statement sets each aside.
fn first_base_asks() -> Vec<(String, String, &'static str, bool)> {
    let mut asks = Vec::new();
    for (number, (_, row, aside)) in first_base_rows_named().iter().enumerate() {
        if row.reference_fields.1 == row.candidate_fields.1 {
            // The two end fields agree.
            continue;
        }
        let (reference, candidate) = first_base_rows(row, number);
        asks.push((reference, candidate, "HGVS_IntronEndOffset", *aside));
    }
    asks
}

/// The format comparison run on each made-up row's two outputs with the class sets aside
/// the end offsets the rule's statement names, and no other.
#[test]
fn the_clipped_reference_first_base_rows_are_set_aside_where_the_statement_says() {
    let asks = first_base_asks();
    assert_eq!(asks.iter().filter(|ask| ask.3).count(), 2);
    assert_eq!(asks.len(), first_base_rows_named().len() - 1);
    let wrong = asks_decided_otherwise(CLIPPED, &asks);
    assert!(wrong.is_empty(), "{wrong:?}");
}

const ADDED_BASES: &str = "hgvs_added_bases_left_at_the_span_end";
const ADDED_BASES_FEATURE: &str = "ENST00000000061";
/// The made-up rows of an allele that adds two bases to the one base it keeps at a
/// reverse-strand transcript's last base: the reference's notation after the sequence name
/// and the candidate's duplication of the last two bases in the same numbering.
const ADDED_PAST_NOTATIONS: [(&str, &str); 5] = [
    ("c.*27+1_*27insAC", "c.*26_*27dup"),
    ("n.83_83+1insAC", "n.82_83dup"),
    ("c.*27+2_*27insAC", "c.*26_*27dup"),
    ("c.*27+1_*27insA", "c.*26_*27dup"),
    ("c.*27+1_*27insTT", "c.*26_*27dup"),
];
const ADDED_PAST_CANDIDATES: [&str; 8] = [
    "none",
    "duplication",
    "duplication elsewhere",
    "shorter duplication",
    "duplication with an offset",
    "duplication with other fields",
    "fields without a notation",
    "insertion",
];

/// One cell of the made-up grid of the class past a last base: whether the last added base
/// repeats the kept base; the reference's notation; the reference's pair; the row's strand;
/// and the candidate's row.
#[derive(Clone, Copy)]
struct AddedPastCell {
    repeats: bool,
    notation: (&'static str, &'static str),
    pair: (&'static str, &'static str),
    strand: Option<&'static str>,
    candidate: &'static str,
}

fn added_past_grid() -> Vec<AddedPastCell> {
    let mut cells = Vec::new();
    for repeats in [true, false] {
        for notation in ADDED_PAST_NOTATIONS {
            for pair in [("1", "0"), ("0", "0")] {
                for strand in [Some("-1"), Some("1"), None] {
                    for candidate in ADDED_PAST_CANDIDATES {
                        cells.push(AddedPastCell {
                            repeats,
                            notation,
                            pair,
                            strand,
                            candidate,
                        });
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a field of a cell: the reference's `ins` stands one
/// base past a coordinate with the pair 1 and 0 on a row that is not a forward strand's;
/// against the candidate's duplication of as many bases ending at that coordinate, where the
/// last added base repeats the kept base, the HGVSc and the start offset; against a
/// candidate that writes nothing, where it does not repeat it, the HGVSc and both offsets.
fn added_past_set_aside(cell: &AddedPastCell, field: &str) -> bool {
    let taken = ["c.*27+1_*27insAC", "n.83_83+1insAC"].contains(&cell.notation.0);
    if !taken || cell.pair != ("1", "0") || cell.strand == Some("1") {
        return false;
    }
    if cell.candidate == "none" {
        return !cell.repeats;
    }
    cell.candidate == "duplication"
        && cell.repeats
        && ["HGVSc", "HGVS_IntronStartOffset"].contains(&field)
}

/// The reference's and the candidate's default-format rows of a cell, under a Location of
/// the cell's own: one kept base, and an allele of the two added bases as the forward strand
/// reads them before it.
fn added_past_rows(cell: &AddedPastCell, number: usize) -> (String, String) {
    let (notation, duplication) = cell.notation;
    let (start, end) = cell.pair;
    let position = 5_000_000 + 10 * number;
    let allele = if cell.repeats { "GTG" } else { "GTA" };
    let columns = format!(
        "cell{number}\t2:{position}\t{allele}\tENSG1\t{ADDED_BASES_FEATURE}\tTranscript\t3_prime_UTR_variant\t900\t-\t-\t-\t-\t-"
    );
    let extra = match cell.strand {
        Some(strand) => format!("IMPACT=MODIFIER;STRAND={strand}"),
        None => "IMPACT=MODIFIER".to_string(),
    };
    let reference = format!(
        "{extra};HGVSc={ADDED_BASES_FEATURE}.1:{notation};HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}"
    );
    let named = |notation: &str| format!(";HGVSc={ADDED_BASES_FEATURE}.1:{notation}");
    let zeros = ";HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=0";
    let written = match cell.candidate {
        "none" => String::new(),
        "duplication" => named(duplication) + zeros,
        "duplication elsewhere" => {
            let elsewhere = duplication
                .replace("26_", "25_")
                .replace("27dup", "26dup")
                .replace("82_", "81_")
                .replace("83dup", "82dup");
            named(&elsewhere) + zeros
        }
        "shorter duplication" => named(&duplication.replace("*26_", "").replace("82_", "")) + zeros,
        "duplication with an offset" => named(duplication) + ";HGVS_OFFSET=-3" + zeros,
        "duplication with other fields" => {
            named(duplication) + ";HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=1"
        }
        "fields without a notation" => zeros.to_string(),
        _ => named(&duplication.replace("dup", "insAC")) + zeros,
    };
    (
        format!("{columns}\t{reference}"),
        format!("{columns}\t{extra}{written}"),
    )
}

/// The value the Extra text of a default-format row holds for `field`.
fn extra_value<'a>(row: &'a str, field: &str) -> Option<&'a str> {
    let extra = row.rsplit('\t').next().unwrap_or_default();
    let held: Vec<(&str, &str)> = extra
        .split(';')
        .filter_map(|item| item.split_once('='))
        .collect();
    value_of(&held, field)
}

/// The fields of the plugin a grid of the class asks about: those that differ between a
/// cell's two rows, a field equal on both being no difference to the comparison.
fn added_past_asks() -> Vec<(String, String, &'static str, bool)> {
    let mut asks = Vec::new();
    for (number, cell) in added_past_grid().iter().enumerate() {
        let (reference, candidate) = added_past_rows(cell, number);
        for field in HGVS_INTRON_OFFSET_FIELDS {
            if extra_value(&reference, field) != extra_value(&candidate, field) {
                asks.push((
                    reference.clone(),
                    candidate.clone(),
                    field,
                    added_past_set_aside(cell, field),
                ));
            }
        }
    }
    asks
}

/// The format comparison run on each cell's two outputs with the class sets aside the
/// fields of the plugin the rule's statement names, and no other.
#[test]
fn the_added_bases_grid_is_set_aside_where_the_statement_says() {
    let asks = added_past_asks();
    assert_eq!(asks.len(), 330);
    let wrong = asks_decided_otherwise(ADDED_BASES, &asks);
    assert!(wrong.is_empty(), "{wrong:?}");
}

const ADDED_OUTSIDE: &str = "hgvs_added_bases_written_outside_the_transcript";
const ADDED_OUTSIDE_FEATURE: &str = "ENST00000000071";

/// The candidate's row of a one-base insertion outside the transcript: what it is, its
/// notation after the sequence name, the two intron offsets of its coordinates, its
/// HGVS_OFFSET.
type OutsideCandidate = (
    &'static str,
    Option<&'static str>,
    Option<(&'static str, &'static str)>,
    Option<&'static str>,
);
const ADDED_OUTSIDE_CANDIDATES: [OutsideCandidate; 13] = [
    (
        "dup, nothing moved",
        Some("n.2601dup"),
        Some(("0", "0")),
        None,
    ),
    ("dup, moved", Some("n.14dup"), Some(("0", "0")), Some("14")),
    (
        "dup, moved the other way",
        Some("n.14dup"),
        Some(("0", "0")),
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
        Some(("0", "0")),
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
        Some(("0", "0")),
        Some("3"),
    ),
    ("del", Some("n.14del"), Some(("0", "0")), Some("14")),
    ("no notation", None, Some(("0", "0")), None),
    ("no fields", Some("n.14dup"), None, Some("14")),
    (
        "offset the distance itself",
        Some("n.14dup"),
        Some(("0", "0")),
        Some("5"),
    ),
];

/// The reference's row: no field, a pair, one field alone, an HGVSc or an HGVS_OFFSET of
/// its own.
#[derive(Clone, Copy, PartialEq)]
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
#[derive(Clone, Copy)]
struct AddedOutsideCell {
    term: &'static str,
    strand: Option<&'static str>,
    candidate: OutsideCandidate,
    reference: OutsideReference,
}

fn added_outside_grid() -> Vec<AddedOutsideCell> {
    let mut cells = Vec::new();
    for term in ABUTTING_GRID_TERMS {
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
fn added_outside_rows(cell: &AddedOutsideCell, number: usize) -> (String, String) {
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
    (
        format!("{columns}\t{reference_text}"),
        format!("{columns}\t{candidate_text}"),
    )
}

/// The plugin's fields of the grid that a cell's two rows write differently.
fn added_outside_asks() -> Vec<(String, String, &'static str, bool)> {
    let mut asks = Vec::new();
    for (number, cell) in added_outside_grid().iter().enumerate() {
        let (reference, candidate) = added_outside_rows(cell, number);
        for field in HGVS_INTRON_OFFSET_FIELDS {
            if extra_value(&reference, field) != extra_value(&candidate, field) {
                asks.push((
                    reference.clone(),
                    candidate.clone(),
                    field,
                    added_outside_set_aside(cell, field),
                ));
            }
        }
    }
    asks
}

/// The format comparison run on each cell's two outputs with the class sets aside the
/// fields of the plugin the rule's statement names, and no other.
#[test]
fn the_bases_written_outside_grid_is_set_aside_where_the_statement_says() {
    let asks = added_outside_asks();
    assert_eq!(asks.len(), 1251);
    let wrong = asks_decided_otherwise(ADDED_OUTSIDE, &asks);
    assert!(wrong.is_empty(), "{wrong:?}");
}

const SHIFT_STOP: &str = "hgvs_insertion_left_where_the_shift_stops";
const SHIFT_STOP_FEATURE: &str = "ENST00000000062";
const SHIFT_STOP_REFERENCES: [&str; 3] = ["ins", "dup", "delins"];
const SHIFT_STOP_CANDIDATES: [&str; 7] = [
    "further dup",
    "further ins",
    "ins not turned",
    "shorter dup",
    "del",
    "none",
    "fields without a notation",
];

/// An allele of 999 bases, which the reference's shift walks at most twice.
fn shift_stop_allele() -> String {
    "ACG".repeat(333)
}

/// The reverse complement of an allele of the four bases.
fn reverse_complement(allele: &str) -> String {
    allele
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// An inserted allele after `steps` steps of a forward walk: its first base moved to its
/// end at each step.
fn turned_forward(allele: &str, steps: usize) -> String {
    let turn = steps % allele.len();
    format!("{}{}", &allele[turn..], &allele[..turn])
}

/// One cell of the made-up grid of the class: the reference's notation where its walk
/// stops; its HGVS_OFFSET (the bound for the allele's length, one more, none); the
/// candidate's row; the candidate's HGVS_OFFSET (the bases its notation stands from the
/// written point, the bound, the other strand's sign, none); the strand; and the plugin's
/// fields (each row's the intron parts of its own notation, the candidate's others, the
/// reference's others, or written on neither row).
#[derive(Clone, Copy)]
struct ShiftStopCell {
    reference: &'static str,
    stopped: Option<&'static str>,
    candidate: &'static str,
    moved: Option<&'static str>,
    strand: &'static str,
    plugin: Option<&'static str>,
}

fn shift_stop_grid() -> Vec<ShiftStopCell> {
    let mut cells = Vec::new();
    for reference in SHIFT_STOP_REFERENCES {
        for stopped in [Some("the bound"), Some("one more"), None] {
            for candidate in SHIFT_STOP_CANDIDATES {
                for moved in [Some("further"), Some("the bound"), Some("other sign"), None] {
                    for strand in ["1", "-1"] {
                        for plugin in [
                            Some("own"),
                            Some("candidate other"),
                            Some("reference other"),
                            None,
                        ] {
                            cells.push(ShiftStopCell {
                                reference,
                                stopped,
                                candidate,
                                moved,
                                strand,
                                plugin,
                            });
                        }
                    }
                }
            }
        }
    }
    cells
}

/// What the rule's statement says of a field of a cell: the reference writes an `ins` or
/// a `dup` beside an HGVS_OFFSET that is the bound, and its two fields, where written, are
/// its notation's; against a candidate `dup` of as many bases or `ins` of the allele's
/// bases turned, beside a larger HGVS_OFFSET of the same sign and its own notation's fields,
/// every field that differs; against a candidate that writes nothing, every field that
/// differs, the HGVS_OFFSET where the candidate writes none or the reference's own.
fn shift_stop_set_aside(cell: &ShiftStopCell, field: &str) -> bool {
    if cell.reference == "delins"
        || cell.stopped != Some("the bound")
        || cell.plugin == Some("reference other")
    {
        return false;
    }
    if cell.candidate == "none"
        || (cell.candidate == "fields without a notation" && cell.plugin.is_none())
    {
        return cell.moved.is_none() || (cell.moved == Some("the bound") && field != "HGVS_OFFSET");
    }
    if !["further dup", "further ins"].contains(&cell.candidate) {
        return false;
    }
    cell.moved == Some("further") && cell.plugin != Some("candidate other")
}

/// The reference's and the candidate's default-format rows of a cell: an insertion of the
/// allele between two bases, described on a transcript without a coding sequence.
fn shift_stop_rows(cell: &ShiftStopCell, number: usize) -> (String, String) {
    let allele = shift_stop_allele();
    let reverse = cell.strand == "-1";
    let sign = if reverse { "-" } else { "" };
    let oriented = if reverse {
        reverse_complement(&allele)
    } else {
        allele.clone()
    };
    let position = 7_000_000 + 10 * number;
    let columns = format!(
        "cell{number}\t2:{position}-{}\t{allele}\tENSG1\t{SHIFT_STOP_FEATURE}\tTranscript\tnon_coding_transcript_variant\t-\t-\t-\t-\t-\t-",
        position + 1
    );
    let extra = format!("IMPACT=MODIFIER;STRAND={}", cell.strand);
    let fields = |pair: Option<(&str, &str)>| match (cell.plugin, pair) {
        (Some(_), Some((start, end))) => {
            format!(";HGVS_IntronEndOffset={end};HGVS_IntronStartOffset={start}")
        }
        _ => String::new(),
    };
    let (notation, pair) = match cell.reference {
        "ins" => (
            format!("n.2012+1_2012+2ins{}", turned_forward(&oriented, 2)),
            ("1", "2"),
        ),
        "dup" => ("n.1014_2012dup".to_string(), ("0", "0")),
        _ => (format!("n.2012_2013delins{oriented}"), ("0", "0")),
    };
    let stopped_text = match cell.stopped {
        Some("the bound") => format!(";HGVS_OFFSET={sign}2"),
        Some(_) => format!(";HGVS_OFFSET={sign}3"),
        None => String::new(),
    };
    let reference_pair = if cell.plugin == Some("reference other") {
        ("7", "7")
    } else {
        pair
    };
    let reference_text = format!(
        "{extra};HGVSc={SHIFT_STOP_FEATURE}.1:{notation}{stopped_text}{}",
        fields(Some(reference_pair))
    );
    let (written, mut own): (Option<String>, Option<(&str, &str)>) = match cell.candidate {
        "further dup" => (Some("n.1017_2015dup".to_string()), Some(("0", "0"))),
        "further ins" => (
            Some(format!(
                "n.2012+4_2012+5ins{}",
                turned_forward(&oriented, 5)
            )),
            Some(("4", "5")),
        ),
        "ins not turned" => (
            Some(format!(
                "n.2012+4_2012+5ins{}",
                turned_forward(&oriented, 4)
            )),
            Some(("4", "5")),
        ),
        "shorter dup" => (Some("n.1018_2015dup".to_string()), Some(("0", "0"))),
        "del" => (Some("n.1017_2015del".to_string()), Some(("0", "0"))),
        "none" => (None, None),
        _ => (None, Some(("0", "0"))),
    };
    let moved_text = match cell.moved {
        Some("further") => format!(";HGVS_OFFSET={sign}5"),
        Some("the bound") => format!(";HGVS_OFFSET={sign}2"),
        Some(_) => format!(";HGVS_OFFSET={}5", if reverse { "" } else { "-" }),
        None => String::new(),
    };
    if written.is_some() && cell.plugin == Some("candidate other") {
        own = Some(("7", "7"));
    }
    let named = written.map_or(String::new(), |written| {
        format!(";HGVSc={SHIFT_STOP_FEATURE}.1:{written}")
    });
    let candidate_text = format!("{extra}{named}{moved_text}{}", fields(own));
    (
        format!("{columns}\t{reference_text}"),
        format!("{columns}\t{candidate_text}"),
    )
}

/// The plugin's fields that differ between the two rows of each cell that writes them.
fn shift_stop_asks() -> Vec<(String, String, &'static str, bool)> {
    let cells: Vec<ShiftStopCell> = shift_stop_grid()
        .into_iter()
        .filter(|cell| cell.plugin.is_some())
        .collect();
    let mut asks = Vec::new();
    for (number, cell) in cells.iter().enumerate() {
        let (reference, candidate) = shift_stop_rows(cell, number);
        for field in HGVS_INTRON_OFFSET_FIELDS {
            if extra_value(&reference, field) != extra_value(&candidate, field) {
                asks.push((
                    reference.clone(),
                    candidate.clone(),
                    field,
                    shift_stop_set_aside(cell, field),
                ));
            }
        }
    }
    asks
}

/// The format comparison run on each cell's two outputs with the class sets aside the
/// fields of the plugin the rule's statement names, and no other.
#[test]
fn the_shift_stop_grid_is_set_aside_where_the_statement_says() {
    let asks = shift_stop_asks();
    assert_eq!(asks.len(), 2544);
    let wrong = asks_decided_otherwise(SHIFT_STOP, &asks);
    assert!(wrong.is_empty(), "{wrong:?}");
}
