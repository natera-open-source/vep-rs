// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the UTRAnnotator classes share: planted 5' UTR rows, the plugin's two
//! fields and comparisons under named classes.

#![allow(dead_code)]

use std::collections::BTreeMap;

use crate::compare_common::*;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::entry::Fields;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::record_rows::RecordRows;
use vep_formats_compare::rules::field_class_accepts;
use vep_formats_compare::rules::utrannotator::UtrannotatorSeedSweep;
use vep_formats_compare::rules::RuleContext;
use vep_formats_compare::rules::Side;

pub const UTR_GAIN: &str = "5_prime_UTR_premature_start_codon_gain_variant";
pub const UTR_LOSS: &str = "5_prime_UTR_premature_start_codon_loss_variant";
pub const UTR_STOP_GAIN: &str = "5_prime_UTR_uORF_stop_codon_gain_variant";
pub const UTR_STOP_LOSS: &str = "5_prime_UTR_uORF_stop_codon_loss_variant";
pub const UTR_FRAMESHIFT: &str = "5_prime_UTR_uORF_frameshift_variant";

pub const PAIR_ORDER: &str = "utrannotator_annotation_pair_order";
pub const UTR_WINDOW: &str = "utrannotator_gained_start_window";
pub const UTR_PRECEDENCE: &str = "utrannotator_consequence_precedence";
pub const UTR_STRUCTURAL: &str = "utrannotator_structural_allele_text";
pub const UTR_LOST_AFTER: &str = "utrannotator_lost_start_after_insertion";
pub const UTR_REMADE: &str = "utrannotator_gained_start_remade_in_place";
pub const UTR_BOUNDARY: &str = "utrannotator_frameshift_at_stop_boundary";
pub const UTR_BEFORE_STOP: &str = "utrannotator_lost_stop_written_before_stop_codon";
pub const UTR_WITHOUT_REFERENCE: &str = "utrannotator_without_reference";

/// The plugin's three counts of the reference 5' UTR, as a row's Extra carries them.
pub const UTR_COUNTS: &str =
    "Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=1;Existing_uORFs=1";

/// The insertion 21 25734954 A ACAT on ENST00000000033 and a five-base deletion.
pub const UTR_INSERTION: (&str, &str, &str, &str) =
    ("utr_ins", "21:25734954-25734955", "CAT", "ENST00000000033");
pub const UTR_DELETION: (&str, &str, &str, &str) = ("del", "21:1000-1004", "-", "ENST2");
pub const UTR_INPUTS: [&str; 3] = [
    "21\t25734954\tutr_ins\tA\tACAT\t.\t.\t.",
    "21\t500\tmnv\tAC\tGT\t.\t.\t.",
    "21\t600\tsnv\tA\tT\t.\t.\t.",
];

/// A made-up insertion of three bases between cDNA positions 146 and 147: the base after it
/// is 5' UTR index 146.
pub const UTR_WHOLE: Row = ("whole_ins", "21:1000-1001", "CCA", "ENST1", "146-147");
pub const UTR_WHOLE_INPUTS: [&str; 3] = [
    "21\t1000\twhole_ins\tG\tGCCA\t.\t.\t.",
    "21\t500\tmnv\tAC\tGT\t.\t.\t.",
    "21\t700\tlonger\tAC\tGTTT\t.\t.\t.",
];

/// A row's Uploaded_variation, Location, Allele, Feature and cDNA_position.
pub type Row<'a> = (&'a str, &'a str, &'a str, &'a str, &'a str);

/// A row's `5UTR_consequence` and `5UTR_annotation`.
pub type Values = (String, String);

/// The classes a comparison names.
pub type Classes<'a> = &'a [&'a str];

/// What a comparison counts under each class.
pub type Applied<'a> = &'a [(&'a str, i64)];

/// The outputs of the seeded runs, one text per seed, where a comparison reads them.
pub type Seeded<'a> = Option<&'a [&'a str]>;

/// The input lines, where a comparison reads them.
pub type Inputs<'a> = Option<&'a [String]>;

pub fn values(term: &str, annotation: &str) -> Values {
    (term.to_string(), annotation.to_string())
}

/// A row on which the plugin reports nothing.
pub fn nothing() -> Values {
    values("-", "-")
}

/// A start codon lost at index 63.
pub fn utr_lost() -> Values {
    values(
        UTR_LOSS,
        "CapDistanceToStart=63:DistanceToCDS=110:KozakStrength=Moderate:type=uORF",
    )
}

/// A frameshift that names a start codon's index.
pub fn utr_shifted() -> Values {
    values(UTR_FRAMESHIFT, "CapDistanceToStart=40:type=uORF")
}

/// A frameshift of one upstream ORF.
pub fn frame_shifted() -> Values {
    values(
        UTR_FRAMESHIFT,
        "KozakStrength=Weak:alt_type=uORF:ref_type=uORF",
    )
}

/// The two fields of a start codon gained at a `CapDistanceToStart`.
pub fn gained(position: i64) -> Values {
    values(
        UTR_GAIN,
        &format!("CapDistanceToStart={position}:DistanceToCDS=110:KozakStrength=Weak:type=uORF"),
    )
}

/// The two fields of a start codon lost at a `CapDistanceToStart`.
pub fn lost_at(position: i64) -> Values {
    values(
        UTR_LOSS,
        &format!(
            "CapDistanceToStart={position}:DistanceToCDS=110:KozakStrength=Moderate:type=uORF"
        ),
    )
}

/// The four leading columns of a row with a cDNA_position.
pub fn at_cdna<'a>(row: (&'a str, &'a str, &'a str, &'a str), cdna: &'a str) -> Row<'a> {
    (row.0, row.1, row.2, row.3, cdna)
}

/// A 5' UTR row of the default format with its cDNA_position and the two UTRAnnotator
/// fields.
pub fn utr_row(row: Row, values: &Values) -> String {
    let (name, location, allele, feature, cdna) = row;
    let extra = format!(
        "IMPACT=MODIFIER;STRAND=1;5UTR_consequence={};5UTR_annotation={}",
        values.0, values.1
    );
    default_row(
        name,
        location,
        allele,
        feature,
        "5_prime_UTR_variant",
        &extra,
    )
    .replacen(
        "\t-\t-\t-\t-\t-\t-\t",
        &format!("\t{cdna}\t-\t-\t-\t-\t-\t"),
        1,
    )
}

/// `utr_row` with `counts` after its two fields, or nothing after them for an empty text.
pub fn counted_row_with(row: Row, values: &Values, counts: &str) -> String {
    let text = utr_row(row, values);
    if counts.is_empty() {
        text
    } else {
        format!("{text};{counts}")
    }
}

/// `utr_row` with the plugin's three counts of the reference 5' UTR after its two fields.
pub fn counted_row(row: Row, values: &Values) -> String {
    counted_row_with(row, values, UTR_COUNTS)
}

/// `counted_row` on a transcript of `strand` (`1` or `-1`), or with no STRAND for `None`.
pub fn stranded_row(row: Row, values: &Values, strand: Option<&str>) -> String {
    let stated = strand.map_or(String::new(), |strand| format!("STRAND={strand};"));
    counted_row(row, values).replacen("STRAND=1;", &stated, 1)
}

/// `stranded_row` on a transcript of `strand` with the HGVSc body `notation` and the
/// `HGVS_OFFSET` where one is given.
pub fn written_row(
    row: Row,
    values: &Values,
    strand: &str,
    notation: &str,
    offset: Option<i64>,
) -> String {
    let mut text = stranded_row(row, values, Some(strand));
    text.push_str(&format!(";HGVSc=ENST1.1:{notation}"));
    if let Some(offset) = offset {
        text.push_str(&format!(";HGVS_OFFSET={offset}"));
    }
    text
}

/// A default-format output of one row.
pub fn one_row(row: String) -> String {
    default_text(&[row])
}

pub fn lines(texts: &[&str]) -> Vec<String> {
    texts.iter().map(|text| text.to_string()).collect()
}

pub fn default_entries(text: &str) -> Vec<Entry> {
    parse_default(text).expect("the text parses")
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

/// The CSQ header line of a VCF output whose entries hold `fields`.
pub fn csq_header(fields: &str) -> String {
    format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {fields}\">\n"
    )
}

/// The input record `line` as an output line: its INFO column, the last, left to follow.
pub fn record_columns(line: &str) -> &str {
    &line[..line.len() - 1]
}

/// `text` as a JSON string.
pub fn quoted(text: &str) -> String {
    pyjson::dumps_str(text, true)
}

/// The seed sweep of Ensembl VEP's default-format outputs, one text per seed.
pub fn seed_sweep<S: AsRef<str>>(outputs: &[S]) -> UtrannotatorSeedSweep {
    let mut sweep = UtrannotatorSeedSweep::new();
    for output in outputs {
        sweep.add_output(&default_entries(output.as_ref()));
    }
    sweep
}

/// What a comparison reported and what it applied.
pub struct Outcome {
    pub report: Value,
    pub applied: BTreeMap<String, i64>,
}

impl Outcome {
    pub fn verdict(&self) -> &str {
        verdict(&self.report)
    }

    pub fn green(&self) -> bool {
        self.verdict() == "green"
    }

    /// Whether the comparison ended on `verdict` with exactly `applied` counted.
    pub fn is(&self, verdict: &str, applied: &[(&str, i64)]) -> bool {
        self.verdict() == verdict && self.applied == count_map(applied)
    }

    /// Whether the comparison is green with `applied` counted where it is `taken`, and red
    /// with nothing counted where it is not.
    pub fn taken_is(&self, taken: bool, applied: &[(&str, i64)]) -> bool {
        if taken {
            self.is("green", applied)
        } else {
            self.is("red", &[])
        }
    }

    /// What the report counts as set aside in `column`, by class; nothing for a column
    /// that holds no such count.
    pub fn masked(&self, column: &str) -> BTreeMap<String, i64> {
        let Some(masked) = at(&self.report, &["columns", column]).get("masked") else {
            return BTreeMap::new();
        };
        masked
            .as_dict()
            .expect("the counts by class")
            .iter()
            .map(|(class, counted)| match counted {
                Value::Int(counted) => (class.clone(), *counted),
                other => panic!("{class} counts {other:?}"),
            })
            .collect()
    }

    pub fn masks(&self, column: &str, applied: &[(&str, i64)]) -> bool {
        self.masked(column) == count_map(applied)
    }

    pub fn differing(&self, column: &str) -> i64 {
        number(&self.report, &["columns", column, "differing"])
    }
}

/// What two entry lists compare to under the built-in `classes`, with the input lines and
/// the seeded runs where they are given.
pub fn compared_entries(
    expected: &[Entry],
    actual: &[Entry],
    classes: &[&str],
    inputs: Option<&[String]>,
    sweep: Option<UtrannotatorSeedSweep>,
) -> Outcome {
    compared_in_order(expected, actual, classes, inputs, sweep, None)
}

/// `compared_entries` with the fields `order_only` compares as multisets.
pub fn compared_in_order(
    expected: &[Entry],
    actual: &[Entry],
    classes: &[&str],
    inputs: Option<&[String]>,
    sweep: Option<UtrannotatorSeedSweep>,
    order_only: Option<&Fields>,
) -> Outcome {
    let mut documented = documented_with(classes);
    let side = Side {
        utrannotator_sweep: sweep,
        ..Side::default()
    };
    let options = CompareOptions {
        input_lines: inputs,
        order_only,
        ..CompareOptions::new(&side)
    };
    let report =
        compare_entries(expected, actual, &mut documented, &options).expect("the entries compare");
    Outcome {
        report,
        applied: counts(&documented.applied),
    }
}

/// What two default-format texts compare to under `classes`.
pub fn compared(reference: &str, candidate: &str, classes: &[&str]) -> Outcome {
    compared_with(reference, candidate, classes, None)
}

/// What two default-format texts compare to under `classes` with the records of `inputs`
/// read.
pub fn compared_with(
    reference: &str,
    candidate: &str,
    classes: &[&str],
    inputs: Option<&[String]>,
) -> Outcome {
    compared_entries(
        &default_entries(reference),
        &default_entries(candidate),
        classes,
        inputs,
        None,
    )
}

/// What two default-format texts compare to under `classes` beside the seeded runs
/// `sweep`, one text per seed.
pub fn compared_swept<S: AsRef<str>>(
    reference: &str,
    candidate: &str,
    sweep: Option<&[S]>,
    classes: &[&str],
) -> Outcome {
    compared_entries(
        &default_entries(reference),
        &default_entries(candidate),
        classes,
        None,
        sweep.map(seed_sweep),
    )
}

/// Whether the class `klass` accepts the differing field `name` of `e` beside the
/// candidate's entry, the record's rows and the seeded runs where they are given.
#[allow(clippy::too_many_arguments)]
pub fn accepts(
    klass: &str,
    e: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
    a: Option<&Entry>,
    record_rows: Option<&[&Entry]>,
    sweep: Option<UtrannotatorSeedSweep>,
) -> bool {
    let side = Side {
        utrannotator_sweep: sweep,
        ..Side::default()
    };
    let context = RuleContext {
        record_rows,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    field_class_accepts(klass, e, name, expected_value, actual_value, a, &context)
        .expect("the class reads the row")
}

/// Whether `klass` accepts the two fields of the row `e`, each under the record's rows of
/// the input lines and under none: the `5UTR_consequence` and the `5utr_annotation` beside
/// the record, then beside none.
pub fn accepts_with_and_without_the_record(
    klass: &str,
    reference: &[Entry],
    inputs: &[String],
    at: usize,
    a: &Entry,
    expected_value: &str,
    actual_value: &str,
) -> [bool; 4] {
    let mut records = RecordRows::new(reference, Some(inputs));
    let e = &reference[at];
    let handed = records.of(e).expect("the record is read");
    let mut found = [false; 4];
    for (given_at, given) in [handed.as_deref(), None].into_iter().enumerate() {
        for (name_at, name) in ["5UTR_consequence", "5utr_annotation"].iter().enumerate() {
            found[given_at * 2 + name_at] = accepts(
                klass,
                e,
                name,
                expected_value,
                actual_value,
                Some(a),
                given,
                None,
            );
        }
    }
    found
}

/// The allele and the `ALT` of each entry a row is handed as its record's.
pub fn handed(records: &mut RecordRows, e: &Entry) -> Option<Vec<(String, Option<String>)>> {
    records.of(e).expect("the record is read").map(|rows| {
        rows.iter()
            .map(|row| (row.allele.clone(), row.fields.get("ALT").cloned()))
            .collect()
    })
}

/// One group of a frameshift's 5UTR_annotation, its pairs in key order: an upstream ORF
/// that starts `start` bases before the coding start, of type `kind` and length `before`
/// in the reference and of length `after` in the mutated sequence.
pub fn shifted_orf_of(start: i64, before: &str, after: &str, kind: &str, kozak: &str) -> String {
    format!(
        "Evidence=False:KozakContext={kozak}:KozakStrength=Strong:alt_type=uORF:alt_type_length={after}:ref_StartDistanceToCDS={start}:ref_type={kind}:ref_type_length={before}"
    )
}

/// `shifted_orf_of` for a uORF of the context `GCCATGG`.
pub fn shifted_orf(start: i64, before: i64, after: i64) -> String {
    shifted_orf_of(
        start,
        &before.to_string(),
        &after.to_string(),
        "uORF",
        "GCCATGG",
    )
}

/// An overlapping ORF that reads through the insertion.
pub fn orf_through() -> String {
    shifted_orf_of(120, "NA", "150", "OutOfFrame_oORF", "CCCATGC")
}

/// A uORF whose stop codon begins 42 bases before the coding start.
pub fn orf_later() -> String {
    shifted_orf_of(99, "60", "72", "uORF", "TTTATGA")
}

/// A uORF whose stop codon begins 54 bases before the coding start.
pub fn orf_at_stop() -> String {
    shifted_orf(87, 36, 45)
}

/// A second uORF whose stop codon begins 54 bases before the coding start.
pub fn orf_at_stop_too() -> String {
    shifted_orf_of(90, "39", "48", "uORF", "AAAATGG")
}

/// A made-up insertion of one base between cDNA positions 146 and 147, which its HGVSc
/// writes between the 55th and the 54th base before the coding start.
pub const UTR_BOUNDARY_ROW: Row = ("grid_ins", "21:1000-1001", "A", "ENST1", "146-147");

/// A counted 5' UTR row that reports `term` with `groups`, with the HGVSc body `notation`
/// and the `HGVS_OFFSET` where one is given.
pub fn shifted_row_of(
    groups: &[String],
    notation: Option<&str>,
    offset: Option<i64>,
    term: &str,
    row: Row,
) -> String {
    let mut text = counted_row(row, &values(term, &groups.join(",")));
    if let Some(notation) = notation {
        text.push_str(&format!(";HGVSc=ENST1.1:{notation}"));
    }
    if let Some(offset) = offset {
        text.push_str(&format!(";HGVS_OFFSET={offset}"));
    }
    text
}

/// `shifted_row_of` for a frameshift on the made-up one-base insertion, written
/// `c.-55_-54insA`.
pub fn shifted_row(groups: &[String]) -> String {
    shifted_row_of(
        groups,
        Some("c.-55_-54insA"),
        None,
        UTR_FRAMESHIFT,
        UTR_BOUNDARY_ROW,
    )
}

/// A made-up insertion of two bases between cDNA positions 146 and 147.
pub const UTR_BEFORE_STOP_ROW: Row = ("grid_ins", "21:1000-1001", "AT", "ENST1", "146-147");

/// One group of a lost stop's 5UTR_annotation, its pairs in key order: the uORF's next
/// stop codon `distance` bases before the coding start, or none.
pub fn lost_stop_group_of(distance: Option<i64>, kozak: &str, strength: &str) -> String {
    let stop = match distance {
        Some(distance) => format!("AltStop=True:AltStopDistanceToCDS={distance}"),
        None => "AltStop=False:AltStopDistanceToCDS=NA".to_string(),
    };
    format!(
        "{stop}:Evidence=False:FrameWithCDS=outOfFrame:KozakContext={kozak}:KozakStrength={strength}"
    )
}

pub fn lost_stop_group(distance: Option<i64>) -> String {
    lost_stop_group_of(distance, "GCCATGG", "Strong")
}

/// A counted 5' UTR row that reports the lost stop with `groups`.
pub fn lost_stop_row_of(groups: &[String], row: Row, counts: &str) -> String {
    counted_row_with(row, &values(UTR_STOP_LOSS, &groups.join(",")), counts)
}

pub fn lost_stop_row(groups: &[String]) -> String {
    lost_stop_row_of(groups, UTR_BEFORE_STOP_ROW, UTR_COUNTS)
}

/// The members of one annotation group as a JSON object, a whole number a number.
pub fn json_members(group: &str) -> String {
    let members: Vec<String> = group
        .split(':')
        .map(|pair| {
            let (key, value) = pair.split_once('=').expect("a pair");
            let value = if !value.is_empty() && value.bytes().all(|byte| byte.is_ascii_digit()) {
                value.to_string()
            } else {
                quoted(value)
            };
            format!("{}: {value}", quoted(key))
        })
        .collect();
    format!("{{{}}}", members.join(", "))
}

/// The groups of an annotation as the JSON format's object of numbered groups.
pub fn json_groups(groups: &[String]) -> String {
    let members: Vec<String> = groups
        .iter()
        .enumerate()
        .map(|(at, group)| format!("\"{}\": {}", at + 1, json_members(group)))
        .collect();
    format!("{{{}}}", members.join(", "))
}

/// One line of a JSON output: the record of `input` under `id` and `allele_string` with
/// one transcript consequence of the members `consequence`, the text between its braces.
pub fn json_line(input: &str, id: &str, allele_string: &str, consequence: &str) -> String {
    format!(
        "{{\"input\": {}, \"id\": {}, \"allele_string\": {}, \"most_severe_consequence\": \"5_prime_UTR_variant\", \"transcript_consequences\": [{{{consequence}}}]}}\n",
        quoted(input),
        quoted(id),
        quoted(allele_string)
    )
}

/// A directory of the test `name` that holds each text of `files` as `<name>.txt`.
pub fn planted_files(name: &str, files: &[(&str, &str)]) -> std::path::PathBuf {
    let directory = crate::support::scratch(name);
    for (file, text) in files {
        std::fs::write(directory.join(format!("{file}.txt")), text).expect("the file is written");
    }
    directory
}

/// The exit code of the tool comparing `reference.txt` and `candidate.txt` of `directory`
/// in the default format, with `more` arguments after them.
pub fn exit_code(directory: &std::path::Path, more: &[&str]) -> i32 {
    let mut args = vec![
        "--format",
        "default",
        "--expected",
        "reference.txt",
        "--actual",
        "candidate.txt",
    ];
    args.extend(more);
    crate::support::run(directory, &args).code
}

/// The plugin's classes in the order a row is put to them.
pub fn plugin_classes() -> Vec<&'static str> {
    vep_formats_compare::classes::tables()
        .field_classes()
        .iter()
        .map(|(class, _)| *class)
        .filter(|class| class.starts_with("utrannotator_"))
        .collect()
}
