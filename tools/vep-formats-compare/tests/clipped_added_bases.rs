// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The clipped-reference, added-bases and shift-stop classes on made rows.
//!
//! A file of `tests/data/clipped_added_bases` holds the questions of one test, a JSON
//! object a line, in order. An `ask` is one field of a pair of default-format rows put to a
//! class, with the answer. A `compare` is two outputs of one format compared under the
//! built-in classes it names, with the verdict, the entries each class set aside and, for
//! every column that counts a class or a difference, those counts; `cache` names the made
//! cache whose exon spans the comparison reads.

mod clipped_added_bases_common;
mod coding_end_common;
mod compare_common;

use std::collections::BTreeMap;
use std::path::Path;

use clipped_added_bases_common::{clamp_cache, CLIPPED};
use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::classes::tables;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Error;
use vep_formats_compare::input::read_text;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::{parse_default, parse_tab};
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::{self, DumpOptions, Value};
use vep_formats_compare::rules::downstream::field_of;
use vep_formats_compare::rules::partial_codon_models::PartialCodonModels;
use vep_formats_compare::rules::sequence_edges::turned_allele;
use vep_formats_compare::rules::transcript_edges::TranscriptEdges;
use vep_formats_compare::rules::{field_class_accepts, RuleContext, Side};

const ADDED_BASES: &str = "hgvs_added_bases_left_at_the_span_end";
const ADDED_OUTSIDE: &str = "hgvs_added_bases_written_outside_the_transcript";
const SHIFT_STOP: &str = "hgvs_insertion_left_where_the_shift_stops";

const SHIFTED_NOTATION_FIELDS: [&str; 6] = [
    "HGVSc",
    "hgvsc",
    "HGVS_OFFSET",
    "hgvs_offset",
    "HGVS_IntronStartOffset",
    "HGVS_IntronEndOffset",
];

/// What a class says of the differing field `field` of the reference's entry `e` and the
/// candidate's entry `a`, asked as a comparison without side tables asks it.
fn accepts(klass: &str, e: &Entry, a: &Entry, field: &str) -> bool {
    answer(klass, e, a, field).expect("the class reads the rows")
}

/// The answer of `accepts`, or the error the class ends the run with.
fn answer(klass: &str, e: &Entry, a: &Entry, field: &str) -> Result<bool, Error> {
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
        klass,
        e,
        field,
        &field_of(Some(e), &[field]),
        &field_of(Some(a), &[field]),
        Some(a),
        &context,
    )
}

fn text_of_member<'a>(record: &'a Value, name: &str) -> &'a str {
    at(record, &[name]).as_str().expect("a text")
}

fn lines_of_member(record: &Value, name: &str) -> Option<Vec<String>> {
    at(record, &[name]).as_list().map(|lines| {
        lines
            .iter()
            .map(|line| line.as_str().expect("a line").to_string())
            .collect()
    })
}

/// The classes a column counts and its differing entries, as a file of questions holds
/// them.
fn column_counts(column: &Value) -> (BTreeMap<String, i64>, i64) {
    let masked = column
        .get("masked")
        .and_then(Value::as_dict)
        .map(|masked| {
            masked
                .iter()
                .map(|(klass, count)| (klass.clone(), whole(count)))
                .collect()
        })
        .unwrap_or_default();
    (masked, column.get("differing").map_or(0, whole))
}

fn whole(value: &Value) -> i64 {
    match value {
        Value::Int(number) => *number,
        other => panic!("{other:?} is no whole number"),
    }
}

/// How many questions of each kind a file held, and how many asks a class took.
#[derive(Debug, PartialEq, Eq)]
struct Replayed {
    asks: usize,
    taken: usize,
    compares: usize,
}

/// Puts every question of the file `name` to the tool and holds each answer to the file's.
fn replay(name: &str) -> Replayed {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/clipped_added_bases")
        .join(format!("{name}.jsonl.gz"));
    let questions = read_text(&path).expect("the questions read");
    let cache = clamp_cache(&scratch(&format!("clipped-added-bases-{name}")));
    let mut replayed = Replayed {
        asks: 0,
        taken: 0,
        compares: 0,
    };
    let mut wrong: Vec<String> = Vec::new();
    for (number, line) in questions.lines().enumerate() {
        let record = json(line);
        let line_number = number + 1;
        if text_of_member(&record, "kind") == "ask" {
            let entry_of = |side: &str| {
                parse_default(&default_text(&[text_of_member(&record, side)]))
                    .expect("the row parses")
                    .remove(0)
            };
            let (klass, field) = (
                text_of_member(&record, "klass"),
                text_of_member(&record, "field"),
            );
            let expected = at(&record, &["expected"]) == &Value::Bool(true);
            replayed.asks += 1;
            replayed.taken += usize::from(expected);
            if accepts(klass, &entry_of("reference"), &entry_of("candidate"), field) != expected {
                wrong.push(format!("line {line_number}: {klass} on {field}"));
            }
            continue;
        }
        replayed.compares += 1;
        let inputs = lines_of_member(&record, "inputs").unwrap_or_default();
        let parsed = |side: &str| {
            let text = text_of_member(&record, side);
            match text_of_member(&record, "format") {
                "default" => parse_default(text),
                "tab" => parse_tab(text),
                "vcf" => parse_vcf(text, &inputs),
                _ => parse_json(text, &inputs),
            }
            .expect("the output parses")
        };
        let classes = texts(at(&record, &["classes"]));
        let mut documented = documented_with(&classes);
        let mut side = Side::default();
        if at(&record, &["cache"]).as_str().is_some() {
            let edge_inputs = lines_of_member(&record, "edge_inputs");
            side.transcript_edges = Some(
                TranscriptEdges::new(&cache, edge_inputs.as_deref()).expect("the cache reads"),
            );
        }
        let report = compare_entries(
            &parsed("expected"),
            &parsed("actual"),
            &mut documented,
            &CompareOptions::new(&side),
        )
        .expect("the outputs compare");
        if verdict(&report) != text_of_member(&record, "verdict") {
            wrong.push(format!(
                "line {line_number}: the verdict is {}",
                verdict(&report)
            ));
        }
        let applied: BTreeMap<String, i64> = at(&record, &["applied"])
            .as_dict()
            .expect("the applied classes")
            .iter()
            .map(|(klass, count)| (klass.clone(), whole(count)))
            .collect();
        if counts(&documented.applied) != applied {
            wrong.push(format!(
                "line {line_number}: applied {:?}",
                counts(&documented.applied)
            ));
        }
        let held = at(&record, &["columns"]).as_dict().expect("the columns");
        for (column, counted) in at(&report, &["columns"]).as_dict().expect("the columns") {
            let expected = held.get(column).map_or_else(Default::default, |pair| {
                let pair = pair.as_list().expect("a column's counts");
                let masked = pair[0]
                    .as_dict()
                    .map(|masked| {
                        masked
                            .iter()
                            .map(|(klass, count)| (klass.clone(), whole(count)))
                            .collect()
                    })
                    .unwrap_or_default();
                (masked, whole(&pair[1]))
            });
            if column_counts(counted) != expected {
                wrong.push(format!(
                    "line {line_number}: the column {column} counts {:?}",
                    column_counts(counted)
                ));
            }
        }
        for column in held.keys() {
            if at(&report, &["columns"]).get(column).is_none() {
                wrong.push(format!("line {line_number}: no column {column}"));
            }
        }
    }
    assert!(wrong.is_empty(), "{name}: {wrong:?}");
    replayed
}

/// On a sequence allele's row whose cDNA position has an open side, the reference's notation
/// of the row's transcript against the candidate's empty notation, or against one of the same
/// transcript and numbering whose every coordinate lies within the row's span plus one of
/// every coordinate of the reference's, is the class, counted once per row in every format.
/// An allele shorter than its span that deletes bases over the edge has no candidate
/// notation, and the reference's `del`, `delins` or substitution against it is the class. A
/// closed cDNA position, a deletion's row, another transcript or numbering, a coordinate
/// beyond the row's reach, an empty reference against a non-deletion, and an empty candidate
/// against a reference inversion or, on the row of an allele as long as its span, against
/// anything but a reference `ins` or `dup` stay differences.
#[test]
fn the_clipped_class_accepts_the_change_inside_the_transcript() {
    assert_eq!(
        replay("the_change_inside_the_transcript"),
        Replayed {
            asks: 0,
            taken: 0,
            compares: 69
        }
    );
}

/// A reference sequence masked over the row's window gives the candidate's substitution an
/// `N` reference base: it is the class as the unmasked base is, counted once per row. An `N`
/// alternate base is not a notation of the class, and the pair stays a difference without
/// the class.
#[test]
fn the_clipped_class_accepts_a_masked_reference_base() {
    assert_eq!(
        replay("a_masked_reference_base"),
        Replayed {
            asks: 0,
            taken: 0,
            compares: 3
        }
    );
}

/// On a substitution's row two `delins` notations show the clamp when each writes a run of
/// the row's allele and the reference's run holds exactly the bases that lie outside the
/// transcript's exon spans, read from the candidate cache: the class, counted once per row
/// in every format, whether the cDNA position is empty or open. A span holding the stop
/// codon is the class on an empty position and on an open one alike, as is an allele over
/// both edges. Without a cache the class takes none of these rows.
#[test]
fn the_clipped_class_accepts_notations_that_show_the_clamp() {
    assert_eq!(
        replay("notations_that_show_the_clamp"),
        Replayed {
            asks: 0,
            taken: 0,
            compares: 17
        }
    );
}

/// An allele that enters from an intron and leaves into the next one is clamped nowhere,
/// and its pair stays a difference against any candidate. So do, on the edge row: bases that
/// are no run of the row's allele; a count that is not the cache's; a candidate whose run
/// lies at the allele's other end or one base off; start coordinates on either side of the
/// stop codon; an intronic substitution; a notation of another kind or none on either side;
/// an insertion's or deletion's row; an allele in lower case or one that is a word; a row of
/// an edited model; a row whose columns leave the cDNA position out; and a closed position.
#[test]
fn the_clipped_class_refuses_notations_without_a_cut() {
    assert_eq!(
        replay("notations_without_a_cut"),
        Replayed {
            asks: 0,
            taken: 0,
            compares: 41
        }
    );
}

/// The row of two `delins` notations that show the clamp, with the HGVSc missing on either
/// side and the other side keeping its own, stays a difference and nothing is set aside.
#[test]
fn a_notation_missing_beside_the_clamp_pair_stays_a_difference() {
    assert_eq!(
        replay("a_notation_missing_beside_the_clamp_pair"),
        Replayed {
            asks: 0,
            taken: 0,
            compares: 3
        }
    );
}

/// Over the made grid the class sets aside a field of the HGVSIntronOffset plugin exactly
/// where its statement does: the reference alone writes the two fields, on a row whose HGVSc
/// the class takes against a candidate that prints none, and the two values are the intron
/// parts of that notation's coordinates as the notation routine reads them before it orders
/// the pair for printing. Every other cell stays a difference.
#[test]
fn the_clipped_class_takes_the_plugins_fields_beside_its_notation_and_nothing_wider() {
    assert_eq!(
        replay("the_plugins_fields_beside_its_notation"),
        Replayed {
            asks: 1080,
            taken: 30,
            compares: 0
        }
    );
}

/// Over the two made grids of the candidate's side the class sets aside a field of the
/// plugin exactly where its statement does: where the reference prints no notation at a
/// reverse-strand transcript's first base and the candidate prints the `ins` of the bases a
/// longer allele adds inside the span with both fields 0, and where both print an `ins` the
/// class takes as a pair, each side's fields the intron parts of its own notation. Every
/// other cell stays a difference.
#[test]
fn the_clipped_class_takes_the_fields_on_the_candidates_side_and_nothing_wider() {
    assert_eq!(
        replay("the_fields_on_the_candidates_side"),
        Replayed {
            asks: 792,
            taken: 11,
            compares: 0
        }
    );
}

/// On a row whose reference prints no notation the class takes the candidate's two fields
/// and not its HGVSc: over the grid of the candidate's side the only HGVSc the class takes is
/// a `del` against the empty reference on a row whose cDNA position has an open side, and on
/// each of the six cells whose field it sets aside the comparison counts the field under the
/// class and the HGVSc column as differing. On the grid of both sides a field is set aside
/// only beside two notations the class takes as a pair.
#[test]
fn the_clipped_class_leaves_the_hgvsc_beside_the_candidates_fields_a_difference() {
    assert_eq!(
        replay("the_hgvsc_beside_the_candidates_fields"),
        Replayed {
            asks: 869,
            taken: 95,
            compares: 6
        }
    );
}

/// Where both sides fill the plugin's two fields the class takes the differing field beside
/// two `ins` notations alone: beside a reference `dup` the candidate's `ins` and its fields
/// of 0 leave both fields differences, though the class takes the row's HGVSc pair.
#[test]
fn the_clipped_class_leaves_the_fields_beside_a_reference_dup_and_a_candidate_ins_a_difference() {
    assert_eq!(
        replay("a_reference_dup_and_a_candidate_ins"),
        Replayed {
            asks: 4,
            taken: 2,
            compares: 0
        }
    );
}

/// Before a forward-strand transcript's first base the reference's `ins` ends on the base
/// before the first exon, written as the first base's coordinate with `-1`, so the plugin's
/// end offset is -1 where the transcript has no intron. Against the candidate's duplication
/// of the bases the allele adds the class sets the end offset aside where the two rows show
/// the clamp at the first base, and the end offset alone: the start offset is not asked of
/// the shape, the HGVS_OFFSET is no field of the class, and the HGVSc is the class's only
/// where its own coordinate test takes the pair. A row that departs in one part stays a
/// difference.
#[test]
fn the_clipped_class_takes_the_end_offset_before_a_forward_strand_first_base_and_nothing_wider() {
    assert_eq!(
        replay("the_end_offset_before_a_first_base"),
        Replayed {
            asks: 41,
            taken: 3,
            compares: 1
        }
    );
}

/// The end offset's shape at a forward-strand transcript's first base takes no row of
/// another shape. Where the allele is as long as its span and the candidate prints no
/// notation, the reference's `ins` of the outside bases and the plugin's two fields beside
/// it are the class's by the shape of a reference that alone fills the fields. Where the
/// candidate prints the duplication against an `ins` between the transcript's first two
/// bases the two pairs are equal and neither field is asked: the HGVSc is the class's by its
/// coordinates and the HGVS_OFFSET stays a difference.
#[test]
fn the_clipped_class_reads_the_rows_beside_the_first_base_row_by_its_other_shapes() {
    assert_eq!(
        replay("the_rows_beside_the_first_base_row"),
        Replayed {
            asks: 24,
            taken: 22,
            compares: 1
        }
    );
}

/// The class reads the HGVSc fields and the plugin's two, and through the comparison itself
/// sets aside all three on the made row, counted once.
#[test]
fn the_clipped_class_reads_the_plugins_fields_and_the_hgvsc() {
    assert_eq!(
        tables().field_classes_of(CLIPPED),
        Some(
            &[
                "HGVSc",
                "hgvsc",
                "HGVS_IntronStartOffset",
                "HGVS_IntronEndOffset"
            ][..]
        )
    );
    assert_eq!(
        replay("the_plugins_fields_and_the_hgvsc"),
        Replayed {
            asks: 0,
            taken: 0,
            compares: 2
        }
    );
}

/// Over the made grid the class sets aside the HGVSc and the plugin's two fields exactly
/// where its statement does: past a reverse-strand transcript's last base the reference
/// writes an `ins` one base past a coordinate beside the pair 1 and 0; where the added bases
/// repeat the transcript's last bases the candidate writes their duplication and the pair 0
/// and 0, and where they repeat nothing it writes no notation and neither field. Every other
/// cell stays a difference.
#[test]
fn the_added_bases_class_takes_the_insertion_past_a_last_base_and_nothing_wider() {
    assert_eq!(
        replay("the_insertion_past_a_last_base"),
        Replayed {
            asks: 1440,
            taken: 20,
            compares: 0
        }
    );
    assert_eq!(
        tables().field_classes_of(ADDED_BASES),
        Some(&SHIFTED_NOTATION_FIELDS[..])
    );
}

/// Over the made grid the class sets aside the HGVSc and the HGVS_OFFSET exactly where its
/// statement does: the reference describes the added bases where the record's span ends,
/// with no HGVS_OFFSET; the candidate describes them the bases its HGVS_OFFSET states
/// further 3', a duplication of as many bases or an `ins` of the same bases turned. Where
/// nothing moves, the candidate's duplication of the bases that end at the reference's
/// `ins`, beside no HGVS_OFFSET, is the class on its HGVSc. Every other cell stays a
/// difference, as do an offset of the other strand's sign and an allele as long as its span.
#[test]
fn the_added_bases_class_takes_the_notation_carried_along_a_run_and_nothing_wider() {
    assert_eq!(
        replay("the_notation_carried_along_a_run"),
        Replayed {
            asks: 1991,
            taken: 37,
            compares: 0
        }
    );
}

/// Through the comparison itself the class sets aside the HGVSc and the start offset of the
/// made row past a last base, counted once, and the row stays a difference without the
/// class. A run without the plugin writes the two fields on neither row: the HGVSc alone
/// differs and is the class's past the last base, along a run and where nothing moves. One
/// row with the plugin's fields and the other without is a difference.
#[test]
fn the_added_bases_class_counts_a_row_once_through_the_comparison() {
    assert_eq!(
        replay("added_bases_through_the_comparison"),
        Replayed {
            asks: 1,
            taken: 0,
            compares: 6
        }
    );
}

/// Past a last base the class reads the kept base that stands beside the added ones. A
/// Location that ends one below its start names no position, so the allele under it holds
/// no kept base: the class ends the run with the comparator's error, which names the
/// Location. Under a Location that ends two below its start the same allele counts one
/// base more than the reference's `ins` adds, and the pair stays a difference. An allele
/// that keeps one base under a Location of one position is the class.
#[test]
fn the_span_end_class_refuses_a_location_that_ends_below_its_start() {
    let pair = |location: &str, allele: &str| {
        let entry_of = |extra: &str| {
            let row = default_row(
                "made",
                location,
                allele,
                "ENST00000000061",
                "3_prime_UTR_variant",
                extra,
            );
            parse_default(&default_text(&[row]))
                .expect("the row parses")
                .remove(0)
        };
        (
            entry_of(
                "IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000000061.1:c.*27+1_*27insAC;\
                 HGVS_IntronEndOffset=0;HGVS_IntronStartOffset=1",
            ),
            entry_of("IMPACT=MODIFIER;STRAND=-1"),
        )
    };
    let asked = |location: &str, allele: &str| {
        let (reference, candidate) = pair(location, allele);
        answer(ADDED_BASES, &reference, &candidate, "HGVSc")
    };
    assert_eq!(asked("2:5", "GTA"), Ok(true));
    assert_eq!(
        asked("2:5-4", "GT"),
        Err(Error::comparator("Location '2:5-4' spans no base"))
    );
    assert_eq!(asked("2:5-3", "GT"), Ok(false));
}

/// Over the made grid the class sets aside the HGVSc, the HGVS_OFFSET and the plugin's two
/// fields exactly where its statement does: on the row of an insertion outside the
/// transcript the reference fills no HGVSc while the candidate writes the notation the 3'
/// rule places on it. On the row of an allele that holds more bases than its span, whose
/// DISTANCE counts from that span, an HGVS_OFFSET of the DISTANCE's size is the class's and a
/// smaller one stays a difference; three inserted bases under a Location of two positions
/// are an insertion's, whose offset is larger.
#[test]
fn the_bases_written_outside_class_takes_the_candidates_notation_and_nothing_wider() {
    assert_eq!(
        replay("bases_written_outside_the_transcript"),
        Replayed {
            asks: 2738,
            taken: 151,
            compares: 0
        }
    );
    assert_eq!(
        tables().field_classes_of(ADDED_OUTSIDE),
        Some(&SHIFTED_NOTATION_FIELDS[..])
    );
}

/// An allele that reads the same on either strand is held to the row's STRAND: on the made
/// cache's forward-strand transcript the pair is the class under STRAND 1 and stays a
/// difference under STRAND -1, whatever its allele reads on the other strand.
#[test]
fn the_clipped_class_refuses_a_row_whose_strand_is_not_the_caches() {
    let cache = clamp_cache(&scratch("clipped-added-bases-strand"));
    let verdict_of_pair = |reference: &str| {
        let side = Side {
            transcript_edges: Some(TranscriptEdges::new(&cache, None).expect("the cache reads")),
            ..Side::default()
        };
        let candidate = reference.replace("c.*5-1_*6delinsACGCG", "c.*5-1_*7delinsACGC");
        let report = compare_entries(
            &parse_default(reference).expect("the reference parses"),
            &parse_default(&candidate).expect("the candidate parses"),
            &mut documented_with(&[CLIPPED]),
            &CompareOptions::new(&side),
        )
        .expect("the outputs compare");
        verdict(&report).to_string()
    };
    let same = default_text(&[default_row(
        "same",
        "21:1001-1006",
        "ACGCGT",
        "ENST00000000001",
        "3_prime_UTR_variant",
        "IMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.*5-1_*6delinsACGCG",
    )]);
    assert_eq!(verdict_of_pair(&same), "green");
    assert_eq!(
        verdict_of_pair(&same.replace("STRAND=1", "STRAND=-1")),
        "red"
    );
}

const SHIFT_STOP_FEATURE: &str = "ENST00000000062";
const SHIFT_STOP_TERM: &str = "non_coding_transcript_variant";
const UNWALKED_FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

/// An insertion of 1,001 bases, which the reference's shift does not walk going 5'.
fn unwalked() -> String {
    "ACGTA".repeat(200) + "C"
}

/// An allele of as many letters that keeps the first base of a record with a second
/// alternate.
fn unwalked_kept() -> String {
    format!("G{}", &unwalked()[..1000])
}

/// One made row of the walk that never starts: whether its allele keeps a base of its
/// record, its STRAND, the coordinates of its `dup` on the row's transcript or none for no
/// HGVSc, and its HGVS_OFFSET or none.
#[derive(Debug, Clone, Copy)]
struct UnwalkedRow {
    kept: bool,
    strand: i64,
    coordinates: Option<&'static str>,
    offset: Option<i64>,
}

fn unwalked_row(
    kept: bool,
    strand: i64,
    coordinates: Option<&'static str>,
    offset: Option<i64>,
) -> UnwalkedRow {
    UnwalkedRow {
        kept,
        strand,
        coordinates,
        offset,
    }
}

impl UnwalkedRow {
    fn allele(&self) -> String {
        if self.kept {
            unwalked_kept()
        } else {
            unwalked()
        }
    }

    fn hgvsc(&self) -> Option<String> {
        self.coordinates
            .map(|coordinates| format!("{SHIFT_STOP_FEATURE}.1:n.{coordinates}dup"))
    }

    /// The input record the row's allele is an alternate of.
    fn record(&self) -> String {
        if self.kept {
            format!("2\t7000000\tins1\tG\t{},T\t.\t.\t.", unwalked_kept())
        } else {
            format!("2\t7000000\tins1\tG\tG{}\t.\t.\t.", unwalked())
        }
    }

    /// The columns of the row in the text formats up to the Extra text.
    fn text_columns(&self) -> String {
        let named = if self.kept {
            format!("2_7000000_G/{}/T\t2:7000000", unwalked_kept())
        } else {
            "ins1\t2:7000000-7000001".to_string()
        };
        format!(
            "{named}\t{}\tENSG1\t{SHIFT_STOP_FEATURE}\tTranscript\t{SHIFT_STOP_TERM}\t-\t-\t-\t-\t-\t-",
            self.allele()
        )
    }

    /// The row's entry as the format `fmt` writes and the tool parses it.
    fn entry(&self, fmt: &str) -> Entry {
        let hgvsc = self.hgvsc();
        let entries = match fmt {
            "default" => {
                let mut extra = format!("IMPACT=MODIFIER;STRAND={}", self.strand);
                if let Some(hgvsc) = &hgvsc {
                    extra.push_str(&format!(";HGVSc={hgvsc}"));
                }
                if let Some(offset) = self.offset {
                    extra.push_str(&format!(";HGVS_OFFSET={offset}"));
                }
                parse_default(&default_text(&[format!(
                    "{}\t{extra}",
                    self.text_columns()
                )]))
            }
            "tab" => {
                let header =
                    DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tHGVSc\tHGVS_OFFSET");
                let row = format!(
                    "{}\tMODIFIER\t{}\t{}\t{}",
                    self.text_columns(),
                    self.strand,
                    hgvsc.as_deref().unwrap_or("-"),
                    self.offset
                        .map_or("-".to_string(), |offset| offset.to_string())
                );
                parse_tab(&format!("{header}\n{row}\n"))
            }
            "vcf" => {
                let header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|STRAND|HGVSc|HGVS_OFFSET\">\n";
                let record = self.record();
                let line = format!(
                    "{}CSQ={}|{SHIFT_STOP_TERM}|{SHIFT_STOP_FEATURE}|Transcript|{}|{}|{}\n",
                    &record[..record.len() - 1],
                    self.allele(),
                    self.strand,
                    hgvsc.as_deref().unwrap_or(""),
                    self.offset
                        .map_or(String::new(), |offset| offset.to_string())
                );
                parse_vcf(&format!("{header}{line}"), &[record])
            }
            _ => {
                let text = |value: &str| Value::Str(value.to_string());
                let mut item = IndexMap::from([
                    ("variant_allele".to_string(), text(&self.allele())),
                    ("transcript_id".to_string(), text(SHIFT_STOP_FEATURE)),
                    (
                        "consequence_terms".to_string(),
                        Value::List(vec![text(SHIFT_STOP_TERM)]),
                    ),
                    ("strand".to_string(), Value::Int(self.strand)),
                ]);
                if let Some(hgvsc) = &hgvsc {
                    item.insert("hgvsc".to_string(), text(hgvsc));
                }
                if let Some(offset) = self.offset {
                    item.insert("hgvs_offset".to_string(), Value::Int(offset));
                }
                let (start, allele_string) = if self.kept {
                    (7000000, format!("G/{}/T", unwalked_kept()))
                } else {
                    (7000001, format!("-/{}", unwalked()))
                };
                let record = self.record();
                let object = Value::Dict(IndexMap::from([
                    ("input".to_string(), text(&record)),
                    ("id".to_string(), text("ins1")),
                    ("seq_region_name".to_string(), text("2")),
                    ("end".to_string(), Value::Int(7000000)),
                    ("strand".to_string(), Value::Int(1)),
                    ("start".to_string(), Value::Int(start)),
                    ("allele_string".to_string(), text(&allele_string)),
                    ("most_severe_consequence".to_string(), text(SHIFT_STOP_TERM)),
                    (
                        "transcript_consequences".to_string(),
                        Value::List(vec![Value::Dict(item)]),
                    ),
                ]));
                let line = pyjson::dumps(&object, &DumpOptions::default()) + "\n";
                parse_json(&line, &[record])
            }
        };
        entries
            .expect("the row parses")
            .into_iter()
            .find(|entry| entry.feature == SHIFT_STOP_FEATURE)
            .expect("the row's entry")
    }
}

/// What the shift-stop class says of each differing field of HGVSc and HGVS_OFFSET on a
/// made pair of rows, the same pair read in the default, tab, VCF and JSON formats.
fn unwalked_verdicts(reference: UnwalkedRow, candidate: UnwalkedRow) -> Vec<Vec<bool>> {
    UNWALKED_FORMATS
        .iter()
        .map(|fmt| {
            let (e, a) = (reference.entry(fmt), candidate.entry(fmt));
            ["HGVSc", "HGVS_OFFSET"]
                .into_iter()
                .filter(|field| field_of(Some(&e), &[field]) != field_of(Some(&a), &[field]))
                .map(|field| accepts(SHIFT_STOP, &e, &a, field))
                .collect()
        })
        .collect()
}

fn in_every_format(expected: &[bool]) -> Vec<Vec<bool>> {
    vec![expected.to_vec(); UNWALKED_FORMATS.len()]
}

/// The row of an insertion of 1,001 bases on a reverse-strand transcript, which the
/// reference does not walk and writes without an HGVS_OFFSET, in the default, tab, VCF and
/// JSON formats. The first two carry the Location in a column; a VCF or JSON entry has it
/// from its record's position and alleles. In each format the candidate's duplication
/// further 3' beside a negative HGVS_OFFSET and its empty row are the class's. The same
/// allele on a forward-strand row, an allele of as many letters that keeps a base of its
/// record and a reference duplication of fewer bases than the allele holds stay differences
/// in each.
#[test]
fn the_shift_stop_class_reads_the_walk_that_never_starts_in_every_format() {
    let (stopped, further, fewer) = (Some("1014_2014"), Some("1017_2017"), Some("1016_2014"));
    let rows: [(&str, UnwalkedRow, UnwalkedRow, &[bool]); 5] = [
        (
            "further 3'",
            unwalked_row(false, -1, stopped, None),
            unwalked_row(false, -1, further, Some(-3)),
            &[true, true],
        ),
        (
            "no notation",
            unwalked_row(false, -1, stopped, None),
            unwalked_row(false, -1, None, None),
            &[true],
        ),
        (
            "forward strand",
            unwalked_row(false, 1, stopped, None),
            unwalked_row(false, 1, further, Some(3)),
            &[false, false],
        ),
        (
            "an allele that keeps a base",
            unwalked_row(true, -1, stopped, None),
            unwalked_row(true, -1, further, Some(-3)),
            &[false, false],
        ),
        (
            "a reference of fewer bases",
            unwalked_row(false, -1, fewer, None),
            unwalked_row(false, -1, None, None),
            &[false],
        ),
    ];
    for (label, reference, candidate, expected) in rows {
        assert_eq!(
            unwalked_verdicts(reference, candidate),
            in_every_format(expected),
            "{label}"
        );
    }
    let located: Vec<(Option<String>, Option<String>)> = UNWALKED_FORMATS
        .iter()
        .map(|fmt| {
            let entry = unwalked_row(false, -1, stopped, None).entry(fmt);
            (entry.location, entry.allele_location)
        })
        .collect();
    let around = Some("2:7000000-7000001".to_string());
    assert_eq!(
        located,
        [
            (around.clone(), None),
            (around.clone(), None),
            (None, around.clone()),
            (None, around)
        ]
    );
}

/// The reference maps each end of an insertion's copy to a coordinate of the transcript
/// alone. Where the 1,001 copied bases begin in one exon and end in the next, the two exonic
/// coordinates of its `dup` lie fewer than 1,001 apart: 501 around an intron of 500 bases,
/// 1,000 around an intron of one base. On the row of the walk that never starts that `dup`
/// is an insertion's own in the default, tab, VCF and JSON formats, beside a candidate that
/// writes nothing and beside a candidate `dup` of 1,001 exonic bases further 3'.
#[test]
fn the_shift_stop_class_takes_a_reference_copy_around_a_whole_intron_in_every_format() {
    let rows: [(&str, UnwalkedRow, UnwalkedRow, &[bool]); 3] = [
        (
            "no notation",
            unwalked_row(false, -1, Some("1014_1514"), None),
            unwalked_row(false, -1, None, None),
            &[true],
        ),
        (
            "further 3', past the intron",
            unwalked_row(false, -1, Some("1014_1514"), None),
            unwalked_row(false, -1, Some("1314_2314"), Some(-800)),
            &[true, true],
        ),
        (
            "an intron of one base",
            unwalked_row(false, -1, Some("1015_2014"), None),
            unwalked_row(false, -1, None, None),
            &[true],
        ),
    ];
    for (label, reference, candidate, expected) in rows {
        assert_eq!(
            unwalked_verdicts(reference, candidate),
            in_every_format(expected),
            "{label}"
        );
    }
}

/// Beside the reference `dup` whose copied bases lie around an intron, each row in the
/// default, tab, VCF and JSON formats. A `dup` of 1,001 bases of one exon is the class's. A
/// `dup` of 999 exonic bases, the count an allele that keeps the two bases of its span adds,
/// one of 1,002, more than the allele holds, and one of a single coordinate stay
/// differences. A row at the bound, the same allele walked 1,000 steps on a forward strand,
/// is the class's whatever its reference `dup` counts, against a candidate further 3' alone.
/// The row of an allele of as many letters that keeps a base of its record, and the same
/// insertion's on a forward strand without an HGVS_OFFSET, stay differences beside the `dup`
/// around an intron.
#[test]
fn the_shift_stop_class_reads_nothing_wider_beside_a_copy_around_an_intron() {
    let nothing = |kept: bool, strand: i64| unwalked_row(kept, strand, None, None);
    let rows: [(&str, UnwalkedRow, UnwalkedRow, &[bool]); 8] = [
        (
            "a copy in one exon",
            unwalked_row(false, -1, Some("1014_2014"), None),
            nothing(false, -1),
            &[true],
        ),
        (
            "two bases fewer in one exon",
            unwalked_row(false, -1, Some("1016_2014"), None),
            nothing(false, -1),
            &[false],
        ),
        (
            "one base more than the allele holds",
            unwalked_row(false, -1, Some("1014_2015"), None),
            nothing(false, -1),
            &[false],
        ),
        (
            "a single coordinate",
            unwalked_row(false, -1, Some("1014"), None),
            nothing(false, -1),
            &[false],
        ),
        (
            "at the bound",
            unwalked_row(false, 1, Some("1014_1514"), Some(1000)),
            unwalked_row(false, 1, Some("1314_2314"), Some(1800)),
            &[true, true],
        ),
        (
            "at the bound, a candidate no further",
            unwalked_row(false, 1, Some("1014_1514"), Some(1000)),
            unwalked_row(false, 1, Some("1314_2314"), Some(1000)),
            &[false],
        ),
        (
            "an allele that keeps a base",
            unwalked_row(true, -1, Some("1014_1514"), None),
            nothing(true, -1),
            &[false],
        ),
        (
            "forward strand",
            unwalked_row(false, 1, Some("1014_1514"), None),
            nothing(false, 1),
            &[false],
        ),
    ];
    for (label, reference, candidate, expected) in rows {
        assert_eq!(
            unwalked_verdicts(reference, candidate),
            in_every_format(expected),
            "{label}"
        );
    }
}

const LAST_PLACEMENT_FEATURE: &str = "ENST00000000063";
const LAST_PLACEMENT_PROTEIN: &str = "ENSP00000000063.1";
const PROTEIN_GATE: &str = "hgvsp_shifted_indel_written_span_gate";
const STOP_NAMED: &str = "hgvsp_stop_named_on_incomplete_coding_sequence";

/// An insertion of 999 bases, which the reference's shift walks two steps at most.
fn walked_twice() -> String {
    "ACG".repeat(333)
}

/// An `ins` of that insertion between the coordinates `first` and `last`, its bases turned
/// by the `steps` its HGVS_OFFSET states.
fn ins_after(steps: i64, first: &str, last: &str) -> String {
    format!(
        "{first}_{last}ins{}",
        turned_allele(&walked_twice(), steps, true)
    )
}

/// One made row of that insertion between 2:7000000 and 2:7000001, written under an
/// identifier, on a forward-strand coding transcript: its term, its CDS position or none,
/// its FLAGS or none, its HGVSc after the `c.` or none, its HGVS_OFFSET or none, and the
/// description of its HGVSp or none.
#[derive(Debug, Clone)]
struct CodingRow {
    term: &'static str,
    cds: Option<(i64, i64)>,
    flags: Option<&'static str>,
    notation: Option<String>,
    offset: Option<i64>,
    protein: Option<&'static str>,
}

/// A row of the insertion written between the coding bases 250 and 251.
fn exon_row(
    notation: Option<&str>,
    offset: Option<i64>,
    protein: Option<&'static str>,
) -> CodingRow {
    CodingRow {
        term: "frameshift_variant",
        cds: Some((250, 251)),
        flags: None,
        notation: notation.map(str::to_string),
        offset,
        protein,
    }
}

/// A row of the insertion written in an intron, 942 bases before the coding base 381.
fn intron_row(
    notation: Option<&str>,
    offset: Option<i64>,
    protein: Option<&'static str>,
) -> CodingRow {
    CodingRow {
        term: "intron_variant",
        cds: None,
        ..exon_row(notation, offset, protein)
    }
}

impl CodingRow {
    fn hgvsc(&self) -> Option<String> {
        self.notation
            .as_ref()
            .map(|notation| format!("{LAST_PLACEMENT_FEATURE}.1:c.{notation}"))
    }

    /// The row's HGVSp, its `=` escaped as the text formats and the VCF format write it.
    fn hgvsp(&self, escaped: bool) -> Option<String> {
        self.protein.map(|description| {
            let description = if escaped {
                description.replace('=', "%3D")
            } else {
                description.to_string()
            };
            format!("{LAST_PLACEMENT_PROTEIN}:p.{description}")
        })
    }

    /// The input record the row's allele is the alternate of.
    fn record(&self) -> String {
        format!("2\t7000000\tins2\tG\tG{}\t.\t.\t.", walked_twice())
    }

    /// The columns of the row in the text formats up to the Extra text.
    fn text_columns(&self) -> String {
        let (cdna, cds, protein) = match self.cds {
            Some((start, end)) => (
                format!("{}-{}", start + 100, end + 100),
                format!("{start}-{end}"),
                ((start + 2) / 3).to_string(),
            ),
            None => ("-".to_string(), "-".to_string(), "-".to_string()),
        };
        format!(
            "ins2\t2:7000000-7000001\t{}\tENSG1\t{LAST_PLACEMENT_FEATURE}\tTranscript\t{}\t{cdna}\t{cds}\t{protein}\t-\t-\t-",
            walked_twice(),
            self.term
        )
    }

    /// The row's output as the format `fmt` writes and the tool parses it.
    fn entries(&self, fmt: &str) -> Vec<Entry> {
        let hgvsc = self.hgvsc();
        let flags = self.flags.map(str::to_string);
        let offset = self.offset.map(|offset| offset.to_string());
        let entries = match fmt {
            "default" => {
                let mut extra = "IMPACT=HIGH;STRAND=1".to_string();
                for (name, value) in [
                    ("FLAGS", flags),
                    ("HGVSc", hgvsc),
                    ("HGVSp", self.hgvsp(true)),
                    ("HGVS_OFFSET", offset),
                ] {
                    if let Some(value) = value {
                        extra.push_str(&format!(";{name}={value}"));
                    }
                }
                parse_default(&default_text(&[format!(
                    "{}\t{extra}",
                    self.text_columns()
                )]))
            }
            "tab" => {
                let header = DEFAULT_HEADER.replace(
                    "\tExtra",
                    "\tIMPACT\tSTRAND\tFLAGS\tHGVSc\tHGVSp\tHGVS_OFFSET",
                );
                let cell = |value: Option<String>| value.unwrap_or_else(|| "-".to_string());
                let row = format!(
                    "{}\tHIGH\t1\t{}\t{}\t{}\t{}",
                    self.text_columns(),
                    cell(flags),
                    cell(hgvsc),
                    cell(self.hgvsp(true)),
                    cell(offset)
                );
                parse_tab(&format!("{header}\n{row}\n"))
            }
            "vcf" => {
                let header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|CDS_position|STRAND|FLAGS|HGVSc|HGVSp|HGVS_OFFSET\">\n";
                let record = self.record();
                let cell = |value: Option<String>| value.unwrap_or_default();
                let cds = self.cds.map(|(start, end)| format!("{start}-{end}"));
                let line = format!(
                    "{}CSQ={}|{}|{LAST_PLACEMENT_FEATURE}|Transcript|{}|1|{}|{}|{}|{}\n",
                    &record[..record.len() - 1],
                    walked_twice(),
                    self.term,
                    cell(cds),
                    cell(flags),
                    cell(hgvsc),
                    cell(self.hgvsp(true)),
                    cell(offset)
                );
                parse_vcf(&format!("{header}{line}"), &[record])
            }
            _ => {
                let text = |value: &str| Value::Str(value.to_string());
                let mut item = IndexMap::from([
                    ("variant_allele".to_string(), text(&walked_twice())),
                    ("transcript_id".to_string(), text(LAST_PLACEMENT_FEATURE)),
                    (
                        "consequence_terms".to_string(),
                        Value::List(vec![text(self.term)]),
                    ),
                    ("strand".to_string(), Value::Int(1)),
                ]);
                if let Some((start, end)) = self.cds {
                    item.insert("cds_start".to_string(), Value::Int(start));
                    item.insert("cds_end".to_string(), Value::Int(end));
                }
                if let Some(flags) = &flags {
                    item.insert("flags".to_string(), Value::List(vec![text(flags)]));
                }
                if let Some(hgvsc) = &hgvsc {
                    item.insert("hgvsc".to_string(), text(hgvsc));
                }
                if let Some(hgvsp) = self.hgvsp(false) {
                    item.insert("hgvsp".to_string(), text(&hgvsp));
                }
                if let Some(offset) = self.offset {
                    item.insert("hgvs_offset".to_string(), Value::Int(offset));
                }
                let record = self.record();
                let object = Value::Dict(IndexMap::from([
                    ("input".to_string(), text(&record)),
                    ("id".to_string(), text("ins2")),
                    ("seq_region_name".to_string(), text("2")),
                    ("end".to_string(), Value::Int(7000000)),
                    ("strand".to_string(), Value::Int(1)),
                    ("start".to_string(), Value::Int(7000001)),
                    (
                        "allele_string".to_string(),
                        text(&format!("-/{}", walked_twice())),
                    ),
                    ("most_severe_consequence".to_string(), text(self.term)),
                    (
                        "transcript_consequences".to_string(),
                        Value::List(vec![Value::Dict(item)]),
                    ),
                ]));
                let line = pyjson::dumps(&object, &DumpOptions::default()) + "\n";
                parse_json(&line, &[record])
            }
        };
        entries.expect("the row parses")
    }

    /// The row's entry as the format `fmt` writes and the tool parses it.
    fn entry(&self, fmt: &str) -> Entry {
        self.entries(fmt)
            .into_iter()
            .find(|entry| entry.feature == LAST_PLACEMENT_FEATURE)
            .expect("the row's entry")
    }
}

/// What the shift-stop class says of the differing HGVSp of a made pair of rows, the same
/// pair read in the default, tab, VCF and JSON formats.
fn protein_verdicts(reference: &CodingRow, candidate: &CodingRow) -> Vec<bool> {
    UNWALKED_FORMATS
        .iter()
        .map(|fmt| {
            let (e, a) = (reference.entry(fmt), candidate.entry(fmt));
            assert_ne!(
                field_of(Some(&e), &["HGVSp"]),
                field_of(Some(&a), &["HGVSp"]),
                "{fmt}"
            );
            accepts(SHIFT_STOP, &e, &a, "HGVSp")
        })
        .collect()
}

/// The report of a made pair of rows compared in the format `fmt` under `classes`, and the
/// classes applied with their counts.
fn coding_rows_compared(
    reference: &CodingRow,
    candidate: &CodingRow,
    fmt: &str,
    classes: &[&str],
) -> (Value, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let side = Side::default();
    let report = compare_entries(
        &reference.entries(fmt),
        &candidate.entries(fmt),
        &mut documented,
        &CompareOptions::new(&side),
    )
    .expect("the outputs compare");
    (report, counts(&documented.applied))
}

/// The classes the HGVSp column of a report counts and its differing entries.
fn hgvsp_column(report: &Value) -> (BTreeMap<String, i64>, i64) {
    let columns = at(report, &["columns"]).as_dict().expect("the columns");
    let (_, column) = columns
        .iter()
        .find(|(name, _)| name.eq_ignore_ascii_case("HGVSp"))
        .expect("the HGVSp column");
    column_counts(column)
}

/// The reference's row of the insertion written in the coding sequence: its walk stops two
/// steps on, where it writes an `ins` and reads the HGVSp `protein`.
fn stopped_in_the_exon(protein: &'static str) -> CodingRow {
    exon_row(Some(&ins_after(2, "252", "253")), Some(2), Some(protein))
}

/// The reference's row of the insertion written in an intron: its walk stops two steps on,
/// 939 bases before the exon, where it writes an `ins` and no HGVSp.
fn stopped_in_the_intron() -> CodingRow {
    intron_row(Some(&ins_after(2, "381-940", "381-939")), Some(2), None)
}

/// The candidate's row of the insertion written in an intron, described at its last
/// placement after the coding base 441, where it reads the HGVSp `protein`.
fn placed_in_the_exon(protein: &'static str) -> CodingRow {
    intron_row(Some("381-938_441dup"), Some(1002), Some(protein))
}

/// The reference reads an insertion's HGVSp where its bounded walk stops
/// (TranscriptVariationAllele.pm:1711) and the candidate where the HGVSc is described, at
/// the run's last placement. Where that placement lies outside the coding sequence the
/// candidate writes no HGVSp: its `dup` ends in an intron or in the 3' UTR, its `ins` stands
/// beside an intron base, or it writes no HGVSc past the transcript's last base. The
/// reference's notation against none is the class's there in the default, tab, VCF and
/// JSON formats, whatever the notation: a frameshift, an insertion, or an extension on a
/// coding sequence no flag marks incomplete.
#[test]
fn the_shift_stop_class_takes_a_reference_hgvsp_where_the_last_placement_leaves_the_coding_sequence(
) {
    let candidates: [(&str, CodingRow); 4] = [
        (
            "a dup that ends in an intron",
            exon_row(Some("255_1252+1dup"), Some(1003), None),
        ),
        (
            "a dup that ends in the 3' UTR",
            exon_row(Some("255_*1dup"), Some(1003), None),
        ),
        (
            "an ins beside an intron base",
            exon_row(Some(&ins_after(1002, "1252", "1252+1")), Some(1002), None),
        ),
        ("no HGVSc", exon_row(None, None, None)),
    ];
    for protein in [
        "Lys84AsnfsTer23",
        "Lys84_Ser85insThrArg",
        "Ter482IleextTer55",
    ] {
        for (label, candidate) in &candidates {
            assert_eq!(
                protein_verdicts(&stopped_in_the_exon(protein), candidate),
                [true; 4],
                "{protein} against {label}"
            );
        }
    }
}

/// Where the record is written outside the coding sequence and the reference's walk stops
/// there, the reference's shifted span has no translation start or end and it writes no
/// HGVSp (TranscriptVariationAllele.pm:1730-1738). The candidate reads the notation at the
/// run's last placement: where its `dup` ends on a coding base, or its `ins` stands
/// between two, a frameshift, an insertion, a duplication or a deletion-insertion against
/// none is the class's in the default, tab, VCF and JSON formats, beside a reference `ins`
/// or `dup`. A row flagged `cds_start_NF` or `cds_end_NF`, on both sides or on one alone,
/// is read as a row without a flag: the candidate writes the notation that the codons of
/// a coding sequence give, whichever end of it is not found.
#[test]
fn the_shift_stop_class_takes_a_candidate_hgvsp_where_the_last_placement_enters_the_coding_sequence(
) {
    let copied = intron_row(Some("381-1939_381-941dup"), Some(2), None);
    for protein in [
        "Asp148ValfsTer12",
        "Asp148ValfsTer?",
        "Gly147_Asp148insAlaSer",
        "Gly147dup",
        "Gly146_Gly147dup",
        "Asp148delinsGluVal",
    ] {
        for reference in [&stopped_in_the_intron(), &copied] {
            assert_eq!(
                protein_verdicts(reference, &placed_in_the_exon(protein)),
                [true; 4],
                "{protein}"
            );
        }
    }
    let inserted = intron_row(
        Some(&ins_after(1002, "441", "442")),
        Some(1002),
        Some("Asp148ValfsTer12"),
    );
    assert_eq!(
        protein_verdicts(&stopped_in_the_intron(), &inserted),
        [true; 4]
    );
    let flagged = |flags: &'static str, row: CodingRow| CodingRow {
        flags: Some(flags),
        ..row
    };
    for flags in ["cds_start_NF", "cds_end_NF"] {
        assert_eq!(
            protein_verdicts(
                &flagged(flags, stopped_in_the_intron()),
                &flagged(flags, placed_in_the_exon("Asp148ValfsTer12"))
            ),
            [true; 4],
            "{flags}"
        );
    }
    assert_eq!(
        protein_verdicts(
            &flagged("cds_end_NF", stopped_in_the_intron()),
            &placed_in_the_exon("Asp148ValfsTer?")
        ),
        [true; 4]
    );
    assert_eq!(
        protein_verdicts(
            &stopped_in_the_intron(),
            &flagged("cds_end_NF", placed_in_the_exon("Asp148ValfsTer12"))
        ),
        [true; 4]
    );
}

/// Every other pair of HGVSp values on a row at the bound stays a difference in the
/// default, tab, VCF and JSON formats: two notations; no candidate notation beside a
/// placement on coding bases; a candidate notation beside a placement in an intron or a
/// UTR; a candidate notation on a row that carries a CDS position on either side; a
/// candidate notation that names the stop alone, an extension or an unchanged residue,
/// on a row flagged `cds_end_NF` too; and a reference notation that names the
/// stop alone on a coding sequence the row shows incomplete, which is
/// `hgvsp_stop_named_on_incomplete_coding_sequence`'s. A row whose reference HGVS_OFFSET is
/// not the bound is no row of the class, whatever its HGVSp.
#[test]
fn the_shift_stop_class_leaves_every_other_hgvsp_pair_a_difference() {
    let frameshift = "Asp148ValfsTer12";
    let stopped = || stopped_in_the_exon("Lys84AsnfsTer23");
    let with_cds = |row: CodingRow| CodingRow {
        cds: Some((380, 381)),
        ..row
    };
    let flagged = |flags: &'static str, row: CodingRow| CodingRow {
        flags: Some(flags),
        ..row
    };
    let rows: Vec<(&str, CodingRow, CodingRow)> = vec![
        (
            "two notations",
            stopped(),
            exon_row(Some("255_1253dup"), Some(1003), Some("Lys418AsnfsTer23")),
        ),
        (
            "two notations, no candidate HGVSc",
            stopped(),
            exon_row(None, Some(2), Some("Lys418AsnfsTer23")),
        ),
        (
            "two notations on a row without a CDS position",
            CodingRow {
                protein: Some("Asp148GlufsTer9"),
                ..stopped_in_the_intron()
            },
            placed_in_the_exon(frameshift),
        ),
        (
            "no candidate notation beside a dup that ends on a coding base",
            stopped(),
            exon_row(Some("255_1253dup"), Some(1003), None),
        ),
        (
            "no candidate notation beside an ins between two coding bases",
            stopped(),
            exon_row(Some(&ins_after(1003, "1253", "1254")), Some(1003), None),
        ),
        (
            "a candidate notation beside a dup that ends in an intron",
            stopped_in_the_intron(),
            intron_row(Some("381-937_441+1dup"), Some(1003), Some(frameshift)),
        ),
        (
            "a candidate notation beside a dup that ends in the 3' UTR",
            stopped_in_the_intron(),
            intron_row(Some("381-937_*1dup"), Some(1003), Some(frameshift)),
        ),
        (
            "a candidate notation beside a dup that ends in the 5' UTR",
            stopped_in_the_intron(),
            intron_row(Some("-1000_-2dup"), Some(1002), Some(frameshift)),
        ),
        (
            "a candidate notation beside an ins beside an intron base",
            stopped_in_the_intron(),
            intron_row(
                Some(&ins_after(1002, "441", "441+1")),
                Some(1002),
                Some(frameshift),
            ),
        ),
        (
            "a candidate notation beside an ins before an exon's first base",
            stopped_in_the_intron(),
            intron_row(
                Some(&ins_after(941, "381-1", "381")),
                Some(941),
                Some(frameshift),
            ),
        ),
        (
            "a candidate notation and no candidate HGVSc",
            stopped_in_the_intron(),
            intron_row(None, None, Some(frameshift)),
        ),
        (
            "a CDS position on both rows",
            with_cds(stopped_in_the_intron()),
            with_cds(placed_in_the_exon(frameshift)),
        ),
        (
            "a CDS position on the reference's row",
            with_cds(stopped_in_the_intron()),
            placed_in_the_exon(frameshift),
        ),
        (
            "a CDS position on the candidate's row",
            stopped_in_the_intron(),
            with_cds(placed_in_the_exon(frameshift)),
        ),
        (
            "the stop unchanged on a row flagged cds_end_NF",
            flagged("cds_end_NF", stopped_in_the_intron()),
            flagged("cds_end_NF", placed_in_the_exon("Ter149=")),
        ),
        (
            "the stop unchanged",
            stopped_in_the_intron(),
            placed_in_the_exon("Ter149="),
        ),
        (
            "the stop extended",
            stopped_in_the_intron(),
            placed_in_the_exon("Ter149GlnextTer5"),
        ),
        (
            "the stop replaced",
            stopped_in_the_intron(),
            placed_in_the_exon("Ter149delinsGlnTer"),
        ),
        (
            "a frameshift named at the stop",
            stopped_in_the_intron(),
            placed_in_the_exon("Ter149GlnfsTer?"),
        ),
        (
            "a residue unchanged",
            stopped_in_the_intron(),
            placed_in_the_exon("Asp148="),
        ),
        (
            "a residue replaced",
            stopped_in_the_intron(),
            placed_in_the_exon("Asp148Val"),
        ),
        (
            "a deletion",
            stopped_in_the_intron(),
            placed_in_the_exon("Asp148del"),
        ),
        (
            "the stop named on an incomplete coding sequence",
            flagged("cds_end_NF", stopped_in_the_exon("Ter482IleextTer55")),
            flagged("cds_end_NF", exon_row(None, None, None)),
        ),
        (
            "a reference HGVS_OFFSET past the bound, a reference notation",
            CodingRow {
                offset: Some(3),
                ..stopped()
            },
            exon_row(Some("255_1252+1dup"), Some(1003), None),
        ),
        (
            "a reference HGVS_OFFSET past the bound, no candidate HGVSc",
            CodingRow {
                offset: Some(3),
                ..stopped()
            },
            exon_row(None, None, None),
        ),
        (
            "a reference HGVS_OFFSET past the bound, a candidate notation",
            CodingRow {
                offset: Some(3),
                ..stopped_in_the_intron()
            },
            placed_in_the_exon(frameshift),
        ),
        (
            "a candidate HGVS_OFFSET at the bound",
            stopped_in_the_intron(),
            CodingRow {
                offset: Some(2),
                ..placed_in_the_exon(frameshift)
            },
        ),
    ];
    for (label, reference, candidate) in &rows {
        assert_eq!(
            protein_verdicts(reference, candidate),
            [false; 4],
            "{label}"
        );
    }
}

/// Through the comparison itself the class sets aside the HGVSp with the HGVSc and the
/// HGVS_OFFSET of its row, the row counted once, in the default, tab, VCF and JSON
/// formats, and the row stays a difference without the class.
#[test]
fn the_shift_stop_class_counts_the_hgvsp_with_its_row_through_the_comparison() {
    let pairs = [
        (
            stopped_in_the_exon("Lys84AsnfsTer23"),
            exon_row(Some("255_1252+1dup"), Some(1003), None),
        ),
        (
            stopped_in_the_exon("Lys84AsnfsTer23"),
            exon_row(None, None, None),
        ),
        (
            stopped_in_the_intron(),
            placed_in_the_exon("Asp148ValfsTer12"),
        ),
    ];
    assert_eq!(
        tables().field_classes_of(SHIFT_STOP),
        Some(
            &[
                "HGVSc",
                "hgvsc",
                "HGVS_OFFSET",
                "hgvs_offset",
                "HGVS_IntronStartOffset",
                "HGVS_IntronEndOffset",
                "HGVSp",
                "hgvsp"
            ][..]
        )
    );
    for (reference, candidate) in &pairs {
        for fmt in UNWALKED_FORMATS {
            let (report, applied) = coding_rows_compared(reference, candidate, fmt, &[SHIFT_STOP]);
            assert_eq!(verdict(&report), "green", "{fmt}");
            assert_eq!(applied, count_map(&[(SHIFT_STOP, 1)]), "{fmt}");
            assert_eq!(
                hgvsp_column(&report),
                (count_map(&[(SHIFT_STOP, 1)]), 0),
                "{fmt}"
            );
            let (report, applied) = coding_rows_compared(reference, candidate, fmt, &[]);
            assert_eq!(verdict(&report), "red", "{fmt}");
            assert_eq!(applied, BTreeMap::new(), "{fmt}");
            assert_eq!(hgvsp_column(&report), (BTreeMap::new(), 1), "{fmt}");
        }
    }
}

/// An HGVSp another class takes is counted under that class. The protein gate stands
/// before the shift-stop class in the order of the field classes, so the candidate's
/// notation on a row both take, the reference's walk stopping between two coding bases, is
/// the gate's where both are named and the shift-stop class's where it is named alone. The
/// reference's extension at a stop that an incomplete coding sequence does not hold is
/// `hgvsp_stop_named_on_incomplete_coding_sequence`'s: under the shift-stop class alone the
/// row's HGVSp stays a difference beside the HGVSc and the HGVS_OFFSET the class takes.
#[test]
fn the_shift_stop_class_leaves_an_hgvsp_of_another_class_to_that_class() {
    let stopped = intron_row(Some(&ins_after(2, "251", "252")), Some(2), None);
    let placed = intron_row(Some("254_1252dup"), Some(1001), Some("Ala84ValfsTer?"));
    let incomplete = CodingRow {
        flags: Some("cds_end_NF"),
        ..stopped_in_the_exon("Ter482IleextTer55")
    };
    let past_the_last_base = CodingRow {
        flags: Some("cds_end_NF"),
        ..exon_row(None, None, None)
    };
    for fmt in UNWALKED_FORMATS {
        let (report, applied) =
            coding_rows_compared(&stopped, &placed, fmt, &[PROTEIN_GATE, SHIFT_STOP]);
        assert_eq!(verdict(&report), "green", "{fmt}");
        assert_eq!(
            applied,
            count_map(&[(PROTEIN_GATE, 1), (SHIFT_STOP, 1)]),
            "{fmt}"
        );
        assert_eq!(
            hgvsp_column(&report),
            (count_map(&[(PROTEIN_GATE, 1)]), 0),
            "{fmt}"
        );
        let (report, applied) = coding_rows_compared(&stopped, &placed, fmt, &[SHIFT_STOP]);
        assert_eq!(verdict(&report), "green", "{fmt}");
        assert_eq!(applied, count_map(&[(SHIFT_STOP, 1)]), "{fmt}");
        assert_eq!(
            hgvsp_column(&report),
            (count_map(&[(SHIFT_STOP, 1)]), 0),
            "{fmt}"
        );

        let (report, applied) =
            coding_rows_compared(&incomplete, &past_the_last_base, fmt, &[SHIFT_STOP]);
        assert_eq!(verdict(&report), "red", "{fmt}");
        assert_eq!(applied, count_map(&[(SHIFT_STOP, 1)]), "{fmt}");
        assert_eq!(hgvsp_column(&report), (BTreeMap::new(), 1), "{fmt}");
        let (report, applied) = coding_rows_compared(
            &incomplete,
            &past_the_last_base,
            fmt,
            &[SHIFT_STOP, STOP_NAMED],
        );
        assert_eq!(verdict(&report), "green", "{fmt}");
        assert_eq!(
            applied,
            count_map(&[(SHIFT_STOP, 1), (STOP_NAMED, 1)]),
            "{fmt}"
        );
        assert_eq!(
            hgvsp_column(&report),
            (count_map(&[(STOP_NAMED, 1)]), 0),
            "{fmt}"
        );
    }
}

/// The HGVSp column and the applied counts of the candidate's notation at the last
/// placement against none, compared in the default format beside a cache of chromosome 2
/// that holds `records`, each under the rows' transcript, its first base as far past
/// 2:6999995 as the shape's lies past `FIRST`: a shape left at `FIRST` holds the rows'
/// Location in its span; `flags` is the FLAGS of both rows.
fn last_placement_beside(
    name: &str,
    records: Vec<coding_end_common::Transcript>,
    flags: Option<&'static str>,
) -> ((BTreeMap<String, i64>, i64), BTreeMap<String, i64>) {
    let flagged = |row: CodingRow| CodingRow { flags, ..row };
    let (reference, candidate) = (
        flagged(stopped_in_the_intron()),
        flagged(placed_in_the_exon("Asp148ValfsTer12")),
    );
    let directory = scratch(&format!("shift-stop-coding-end-{name}"));
    let records: Vec<coding_end_common::Transcript> = records
        .into_iter()
        .map(|shape| coding_end_common::Transcript {
            id: LAST_PLACEMENT_FEATURE,
            first: 6_999_995 + (shape.first - coding_end_common::FIRST),
            ..shape
        })
        .collect();
    let cache = coding_end_common::planted_cache_on(&directory, "2", &records);
    let side = Side {
        partial_codon_models: Some(PartialCodonModels::new(&cache).expect("the cache")),
        ..Side::default()
    };
    let mut documented = documented_with(&[SHIFT_STOP]);
    let report = compare_entries(
        &reference.entries("default"),
        &candidate.entries("default"),
        &mut documented,
        &CompareOptions::new(&side),
    )
    .expect("the outputs compare");
    (hgvsp_column(&report), counts(&documented.applied))
}

/// The pair is the class's: one HGVSp set aside under it and none left.
fn last_placement_taken(read: ((BTreeMap<String, i64>, i64), BTreeMap<String, i64>)) -> bool {
    let ((masked, differing), applied) = read;
    masked.get(SHIFT_STOP) == Some(&1) && differing == 0 && applied.get(SHIFT_STOP) == Some(&1)
}

/// The rows' FLAGS of the two readings of a pair: none printed, and `cds_end_NF`.
const LAST_PLACEMENT_FLAGS: [(&str, Option<&str>); 2] =
    [("unflagged", None), ("flagged", Some("cds_end_NF"))];

/// Where the candidate's cache holds a record of the row's transcript that shows a stop
/// codon at its coding end, a whole one or one that polyadenylation completes, the class
/// reads the two rows alone as it does without a cache: the candidate's notation against
/// none is the class's on rows that print no flag and on rows that print `cds_end_NF`.
#[test]
fn the_shift_stop_class_reads_the_two_rows_alone_beside_a_record_with_a_closing_stop_codon() {
    use coding_end_common::{closed, shown};
    for (name, shape) in [("closed", closed()), ("completed", shown())] {
        for (rows, flags) in LAST_PLACEMENT_FLAGS {
            assert!(
                last_placement_taken(last_placement_beside(
                    &format!("{name}-rows-{rows}"),
                    vec![shape.clone()],
                    flags
                )),
                "{name}, rows {rows}"
            );
        }
    }
}

/// Where the cache holds a record of the row's transcript that shows no stop codon at its
/// coding end, a record with `cds_end_NF` or without, the class reads the two rows alone
/// as it does without a cache: the candidate's notation against none is the class's on
/// rows that print no flag and on rows that print `cds_end_NF`.
#[test]
fn the_shift_stop_class_reads_the_two_rows_alone_beside_a_record_without_a_closing_stop_codon() {
    use coding_end_common::{flagged, open, open_tail};
    for (name, shape) in [
        ("open", open()),
        ("open-tail", open_tail()),
        ("flagged", flagged()),
    ] {
        for (rows, flags) in LAST_PLACEMENT_FLAGS {
            assert!(
                last_placement_taken(last_placement_beside(
                    &format!("{name}-rows-{rows}"),
                    vec![shape.clone()],
                    flags
                )),
                "{name}, rows {rows}"
            );
        }
    }
}

/// Where the cache gives no reading of the row's transcript, it holds no record of the
/// transcript at the row's Location or a record without a coding sequence, the class
/// reads the two rows alone as it does without a cache: the candidate's notation against
/// none is the class's on rows that print no flag and on rows that print `cds_end_NF`.
#[test]
fn the_shift_stop_class_reads_the_two_rows_alone_where_the_cache_gives_no_reading() {
    use coding_end_common::{closed, Transcript, FIRST};
    let elsewhere = || Transcript {
        first: FIRST + 100_000,
        ..closed()
    };
    let without_a_coding_sequence = || Transcript {
        coding: String::new(),
        ..closed()
    };
    for (name, records) in [
        ("absent", Vec::new()),
        ("elsewhere", vec![elsewhere()]),
        ("no-coding-sequence", vec![without_a_coding_sequence()]),
    ] {
        assert!(
            last_placement_taken(last_placement_beside(name, records.clone(), None)),
            "{name}"
        );
        assert!(
            last_placement_taken(last_placement_beside(
                &format!("{name}-flagged"),
                records,
                Some("cds_end_NF")
            )),
            "{name}"
        );
    }
    // Without a cache the pair is read the same way.
    let (reference, candidate) = (
        stopped_in_the_intron(),
        placed_in_the_exon("Asp148ValfsTer12"),
    );
    let (report, counted) = coding_rows_compared(&reference, &candidate, "default", &[SHIFT_STOP]);
    assert_eq!(hgvsp_column(&report).1, 0);
    assert_eq!(counted.get(SHIFT_STOP), Some(&1));
}
