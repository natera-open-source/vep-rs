// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison-mode report for one format.

#![allow(clippy::too_many_arguments)]

use crate::argparse::path_text;
use crate::classes::CONSERVATION_SHORT_NAME_SPELLING;
use crate::classes::ORDER_ONLY;
use crate::classes::REGULATORY_FEATURE_ROWS;
use crate::compare::compare_entries;
use crate::compare::CompareOptions;
use crate::consts::CONSERVATION_KEY_SPELLING;
use crate::consts::MISSPELT_CONSERVATION_KEY;
use crate::counter::Counter;
use crate::entry::Entry;
use crate::error::Error;
use crate::error::Result;
use crate::floors::json_floor_fields;
use crate::floors::json_floor_location;
use crate::floors::StrListOrStr;
use crate::input::read_text;
use crate::lists::declared_columns;
use crate::lists::order_only_fields;
use crate::manifest::Documented;
use crate::parse::json::parse_json_lines;
use crate::parse::parquet::parquet_as_tab;
use crate::parse::parquet::CanonicalText;
use crate::parse::recoder::parse_recoder;
use crate::parse::reference_lines_masked;
use crate::parse::reference_text_masked;
use crate::parse::text::parse_default;
use crate::parse::text::parse_tab;
use crate::parse::vcf::parse_vcf_lines;
use crate::parse::{parse_format, parse_tab_beside_parquet};
use crate::pyjson::Value;
use crate::regulatory::drop_regulatory_rows;
use crate::regulatory::keep_regulatory_rows;
use crate::remainder::hand_over;
use crate::remainder::hand_over_canonical;
use crate::remainder::FrontEnd;
use crate::remainder::KeptBack;
use crate::remainder::REMAINDER_LINE_UNPLACED;
use crate::rules::term_pairs::read_entry_references;
use crate::rules::Side;
use crate::stream::Text;
use crate::COMPONENT;
use indexmap::IndexMap;
use std::collections::BTreeSet;
use std::path::Path;

pub const FLOOR_UNDECLARED_NOT_APPLICABLE: &str = "a column neither output's header declares in this format (the tab header line and the CSQ Format list carry the fields the run's options name, where the default format's Extra carries every key a row has), so the floor is checked in the default format of the same case";

/// The side tables the recoder format, a comparison of no annotation rows, leaves out,
/// taken out of `side`.
fn annotation_tables(side: &mut Side) -> Side {
    Side {
        structural_variants: side.structural_variants.take(),
        sequence_edges: side.sequence_edges.take(),
        codon_tables: side.codon_tables.take(),
        downstream_sequences: side.downstream_sequences.take(),
        revel_table: side.revel_table.take(),
        splice_windows: side.splice_windows.take(),
        spliceai_table: side.spliceai_table.take(),
        ..Side::default()
    }
}

/// Puts the tables `annotation_tables` took back into `side`.
fn restore_annotation_tables(side: &mut Side, taken: Side) {
    side.structural_variants = taken.structural_variants;
    side.sequence_edges = taken.sequence_edges;
    side.codon_tables = taken.codon_tables;
    side.downstream_sequences = taken.downstream_sequences;
    side.revel_table = taken.revel_table;
    side.splice_windows = taken.splice_windows;
    side.spliceai_table = taken.spliceai_table;
}

/// The entries of a remainder's lines, read as `parse_format` reads an output of those
/// lines. `masked` names the reference side, whose text classes are applied to the lines
/// before they are parsed (`reference_text_masked`). `ordinals` gives the output ordinal
/// of each line in turn in a VCF or JSON comparison (`KeptBack::ordinals`), which needs
/// the input records as well.
fn parse_remainder(
    fmt: &str,
    mut lines: Vec<String>,
    masked: Option<&mut Documented>,
    ordinals: Option<Vec<i64>>,
    inputs: Option<&[String]>,
) -> Result<Vec<Entry>> {
    if fmt == "vcf" || fmt == "json" {
        if let Some(documented) = masked {
            reference_lines_masked(&mut lines, documented);
        }
        let Some(inputs) = inputs else {
            return Err(Error::comparator(format!(
                "--format {fmt} needs --input, the input VCF"
            )));
        };
        let ordinals = ordinals.unwrap_or_default();
        let entries = if fmt == "vcf" {
            parse_vcf_lines(&lines, &ordinals, inputs)
        } else {
            parse_json_lines(&lines, &ordinals, inputs)
        };
        // A line past the last of the remainder's line runs.
        return entries.map_err(|error| match &error {
            Error::Comparator(text) if text.ends_with("has no ordinal in its output") => {
                Error::comparator(REMAINDER_LINE_UNPLACED)
            }
            _ => error,
        });
    }
    let mut text = lines.join("\n");
    if !lines.is_empty() {
        text.push('\n');
    }
    drop(lines);
    if let Some(documented) = masked {
        text = reference_text_masked(&text, documented);
    }
    if fmt == "default" {
        parse_default(&text)
    } else {
        parse_tab(&text)
    }
}

/// `json_floor_location` with `present`, the members the entries a front end kept back
/// carry (`KeptBack::present`): a floor one of whose fields a kept-back entry carries is
/// counted under its fields.
fn json_floor_location_beside(
    column: &str,
    expected: &[Entry],
    actual: &[Entry],
    present: Option<&BTreeSet<String>>,
) -> Option<StrListOrStr> {
    let fields = json_floor_fields(column)?;
    if present.is_some_and(|present| fields.iter().any(|field| present.contains(field))) {
        return Some(StrListOrStr::StrList(fields));
    }
    json_floor_location(column, expected, actual)
}

/// The note of a Parquet comparison that ends on its candidate's failure with no line of
/// the reference held.
pub const CANDIDATE_NOT_READ: &str =
    "the comparison ended on a candidate that is not read: no line of the reference was held";

/// What the front end keeps back of a Parquet comparison, or `None` for one made on the
/// whole outputs; `front_end.note` says which.
///
/// The stream of record groups reads each output in its canonical form with the rows in
/// the output's own order (`CanonicalText`), as two tab outputs whose order is not read:
/// a record group both outputs hold line for line is never held. The whole outputs are
/// compared instead where the reference cannot be read so, where a row of either output
/// holds a line end, which the sorted canonical form keeps in its row, where a row of a
/// Parquet output opens with `#`, a header line that the sorted canonical form holds
/// among other rows than the output's own order does, and where the stream ends on an
/// error of its own. A candidate that cannot be read beside a reference that was read
/// ends the comparison with the adapter's own failure, read off the candidate alone, once
/// the reference's lines have ended it with none of their own: where no line of the
/// reference may end it, none is held and the note is `CANDIDATE_NOT_READ`.
fn hand_over_parquet(
    front_end: &mut FrontEnd,
    expected_path: &Path,
    actual_path: &Path,
    documented: &Documented,
    keep_only: Option<&str>,
) -> Result<Option<KeptBack>> {
    let expected = CanonicalText::beside(expected_path, actual_path);
    let actual = CanonicalText::of(actual_path);
    let handed = hand_over_canonical(
        front_end,
        expected.clone(),
        actual.clone(),
        documented,
        keep_only,
    );
    let mut whole = |outside: &str| {
        front_end.note =
            format!("the whole outputs were compared: {outside} is outside the front end");
        Ok(None)
    };
    if expected.unread() {
        return whole("a reference that is not read as its canonical form");
    }
    if expected.split() || actual.split() {
        return whole("a row that holds a line end");
    }
    if expected.header_row() || actual.header_row() {
        return whole("a row that opens with #");
    }
    if actual.unread() {
        let Err(failure) = parquet_as_tab(actual_path, true) else {
            return whole("a candidate that is not read as its canonical form");
        };
        // The stream words the failure of a reading with `Text::failed` of its text, and
        // ends on the candidate's where no line of the reference may end the comparison.
        let Error::Comparator(unread) = Text::failed(actual.as_ref(), "");
        if matches!(&handed, Err(Error::Comparator(error)) if error.starts_with(&unread)) {
            front_end.note = CANDIDATE_NOT_READ.to_string();
            return Err(failure);
        }
        return match handed {
            Ok(Some(mut kept)) if kept.candidate_error.is_some() => {
                kept.candidate_error = Some(failure);
                Ok(Some(kept))
            }
            _ => whole("a candidate that is not read as its canonical form"),
        };
    }
    match handed {
        Ok(kept) => Ok(kept),
        Err(_) => whole("a canonical form the stream of record groups does not read"),
    }
}

/// A count map as the report prints it: its names in byte order.
fn sorted_counts(counts: &Counter) -> Value {
    let mut items: Vec<(&str, i64)> = counts.iter().collect();
    items.sort_unstable();
    Value::Dict(
        items
            .into_iter()
            .map(|(name, counted)| (name.to_string(), Value::Int(counted)))
            .collect(),
    )
}

fn optional_text(value: Option<&str>) -> Value {
    value.map_or(Value::Null, |value| Value::Str(value.to_string()))
}

fn optional_count(value: Option<usize>) -> Value {
    value.map_or(Value::Null, |value| {
        Value::Int(i64::try_from(value).unwrap_or(i64::MAX))
    })
}

/// The comparison-mode report for one format. The classes applied are `documented.masks`
/// alone (`masks` is the same set the caller built `documented` from); the
/// structural-variant classes among them read `structural_variants`, which the recoder
/// format, a comparison of no annotation rows, leaves out, and
/// `hgvsp_intron_end_insertion` and `hgvsp_intron_start_insertion` read `intron_ends`,
/// and `maxentscan_overlapping_windows` reads `splice_windows`.
/// `unordered`
/// compares the shared entries without their order (a run whose reference rows Ensembl VEP
/// emits in hash order); a Parquet actual is compared without order in any case.
/// `hgvs_shift_sequence_edge_window` reads `sequence_edges`, which the recoder format
/// leaves out as well.
/// `organelle_peptide_paths_translation_table` reads `codon_tables`, which the recoder
/// format leaves out as well.
/// The classes of `TRANSCRIPT_READING_CLASSES` read `downstream_sequences`, which the
/// recoder format leaves out as well.
/// `keep_only` (`regulatory_feature_rows`, the one class) keeps the regulatory rows of
/// both sides and the record entries of the records that have one (`only_regulatory_rows`),
/// and is refused beside the `regulatory_feature_rows` mask, which drops them. A floor
/// is read where the format carries its column: under `conservation_short_name_spelling`
/// a floor spelt with the misspelt key is read under the renamed one; in the JSON format
/// each floor is placed by `json_floor_location`; in the tab and VCF formats a floor on a
/// column neither side's header declares (`declared_columns`) is not applicable
/// (`FLOOR_UNDECLARED_NOT_APPLICABLE`), since every row of those formats carries every
/// column its header names.
/// A default, tab, VCF or JSON comparison is handed over to the stream of record groups
/// (`hand_over`): the entries are then parsed from the remainder and `compare_entries` adds
/// what was kept back. A Parquet comparison is handed over in the canonical form of its
/// two outputs (`hand_over_parquet`), the remainder's rows in the order the whole
/// canonical form has them in. `VEP_FORMATS_COMPARE_WHOLE=1` in the environment compares
/// the whole outputs instead; the report is the same. `records_path` is the file `inputs`
/// was read from, which the stream reads beside a VCF or JSON comparison.
/// `input_path` names the input VCF under a class of `RECORD_REFERENCE_CLASSES`: the
/// default and tab entries the class can read take their record's reference bases from it
/// (`read_entry_references`).
/// A class of `documented.corrected` is applied nowhere: the report names each, with the
/// release that corrects its defect, under `masks.corrected`, a member a report of no such
/// class lacks.
///
/// `unordered` is false, `keep_only` is `None`, `input_path` is `None` and `records_path` is
/// `None` in a comparison that names none of them; the side tables are read from `side`, which
/// is left as it was given.
pub fn compare_files(
    fmt: &str,
    expected_path: &Path,
    actual_path: &Path,
    inputs: Option<&[String]>,
    documented: &mut Documented,
    masks: &BTreeSet<String>,
    floors: &IndexMap<String, i64>,
    field_pairs: &IndexMap<String, BTreeSet<(String, String)>>,
    examples: i64,
    unordered: bool,
    keep_only: Option<&str>,
    input_path: Option<&Path>,
    records_path: Option<&Path>,
    side: &mut Side,
) -> Result<Value> {
    let mut front_end = FrontEnd::stream();
    let report = compare_files_through(
        front_end.as_mut(),
        fmt,
        expected_path,
        actual_path,
        inputs,
        documented,
        masks,
        floors,
        field_pairs,
        examples,
        unordered,
        keep_only,
        input_path,
        records_path,
        side,
    );
    if let Some(restarts) = front_end.map(|front_end| front_end.restarts) {
        if restarts > 0 {
            eprintln!(
                "[{COMPONENT}] the outputs were read again {restarts} times: two record groups of one digest differed in their text"
            );
        }
    }
    report
}

/// `compare_files` with the front end named: `None` compares the whole outputs, and a
/// front end's `note` says what the comparison did with it.
pub fn compare_files_through(
    front_end: Option<&mut FrontEnd>,
    fmt: &str,
    expected_path: &Path,
    actual_path: &Path,
    inputs: Option<&[String]>,
    documented: &mut Documented,
    _masks: &BTreeSet<String>,
    floors: &IndexMap<String, i64>,
    field_pairs: &IndexMap<String, BTreeSet<(String, String)>>,
    examples: i64,
    unordered: bool,
    keep_only: Option<&str>,
    input_path: Option<&Path>,
    records_path: Option<&Path>,
    side: &mut Side,
) -> Result<Value> {
    if let Some(keep_only) = keep_only.filter(|klass| *klass != REGULATORY_FEATURE_ROWS) {
        return Err(Error::comparator(format!(
            "unknown --keep-only class: {keep_only}"
        )));
    }
    if keep_only == Some(REGULATORY_FEATURE_ROWS)
        && documented.masks.contains(REGULATORY_FEATURE_ROWS)
    {
        return Err(Error::comparator(
            "--keep-only regulatory_feature_rows conflicts with --masks regulatory_feature_rows",
        ));
    }
    // The column name as the parsed reference spells it: under
    // `conservation_short_name_spelling` the misspelt conservation key renamed.
    let spelling = CONSERVATION_KEY_SPELLING.replace("\\1", "${1}");
    let renames = documented.masks.contains(CONSERVATION_SHORT_NAME_SPELLING);
    let renamed_key = |name: &str| -> String {
        if renames {
            MISSPELT_CONSERVATION_KEY
                .replace_all(name, spelling.as_str())
                .into_owned()
        } else {
            name.to_string()
        }
    };

    let floors: IndexMap<String, i64> = floors
        .iter()
        .map(|(name, minimum)| (renamed_key(name), *minimum))
        .collect();
    let mut kept_back = match front_end {
        Some(front_end) if fmt == "parquet" => {
            hand_over_parquet(front_end, expected_path, actual_path, documented, keep_only)?
        }
        Some(front_end) => {
            let floor_columns: Vec<&str> = floors.keys().map(String::as_str).collect();
            hand_over(
                front_end,
                fmt,
                expected_path,
                actual_path,
                documented,
                unordered,
                keep_only,
                &floor_columns,
                records_path,
            )?
        }
        None => None,
    };
    let mut expected: Vec<Entry>;
    let mut actual: Vec<Entry>;
    if fmt == "recoder" {
        if keep_only.is_some() || documented.masks.contains(REGULATORY_FEATURE_ROWS) {
            return Err(Error::comparator(
                "the recoder format has no regulatory rows to keep or drop",
            ));
        }
        expected = parse_recoder(&read_text(expected_path)?, documented)?;
        actual = parse_recoder(&read_text(actual_path)?, documented)?;
    } else {
        let expected_format = if fmt == "parquet" {
            "tab-canonical"
        } else {
            fmt
        };
        match kept_back.as_mut() {
            None => {
                let mut rewrite = |text: &str| reference_text_masked(text, documented);
                expected = if fmt == "parquet" {
                    parse_tab_beside_parquet(expected_path, actual_path, Some(&mut rewrite))?
                } else {
                    parse_format(expected_format, expected_path, inputs, Some(&mut rewrite))?
                };
                actual = parse_format(fmt, actual_path, inputs, None)?;
            }
            Some(kept) => {
                let lines = std::mem::take(&mut kept.expected);
                let ordinals = kept.ordinals("expected");
                expected =
                    parse_remainder(expected_format, lines, Some(documented), ordinals, inputs)?;
                if let Some(error) = kept.candidate_error.take() {
                    return Err(error);
                }
                let lines = std::mem::take(&mut kept.actual);
                actual = parse_remainder(fmt, lines, None, kept.ordinals("actual"), inputs)?;
            }
        }
    }
    if let Some(kept) = kept_back.as_mut() {
        kept.keys = expected.iter().map(|e| e.key.clone()).collect();
    }
    if documented.masks.contains(REGULATORY_FEATURE_ROWS) {
        let written = expected.len();
        drop_regulatory_rows(&mut expected)?;
        // A class that dropped nothing is unused, as a manifest class matching no item is.
        if expected.len() < written {
            let dropped = i64::try_from(written - expected.len()).unwrap_or(i64::MAX);
            documented.applied.add(REGULATORY_FEATURE_ROWS, dropped);
        }
    }
    if keep_only == Some(REGULATORY_FEATURE_ROWS) {
        keep_regulatory_rows(&mut expected)?;
        keep_regulatory_rows(&mut actual)?;
    }
    if let Some(input_path) = input_path {
        // The entries take their record's reference bases here.
        read_entry_references(&mut expected, &mut actual, input_path)?;
    }
    let mut floor_fields: Option<IndexMap<String, Option<StrListOrStr>>> = None;
    if fmt == "json" {
        let present = kept_back.as_ref().map(|kept| &kept.present);
        floor_fields = Some(
            floors
                .keys()
                .map(|name| {
                    let location = json_floor_location_beside(name, &expected, &actual, present);
                    (name.clone(), location)
                })
                .collect(),
        );
    } else if fmt == "tab" || fmt == "vcf" {
        let mut declared: BTreeSet<String> = BTreeSet::new();
        for path in [expected_path, actual_path] {
            declared.extend(
                declared_columns(fmt, path)?
                    .iter()
                    .map(|name| renamed_key(name)),
            );
        }
        floor_fields = Some(
            floors
                .keys()
                .filter(|name| !declared.contains(*name))
                .map(|name| {
                    let note = StrListOrStr::Str(FLOOR_UNDECLARED_NOT_APPLICABLE.to_string());
                    (name.clone(), Some(note))
                })
                .collect(),
        );
    }
    let order_only = documented
        .masks
        .contains(ORDER_ONLY)
        .then(|| order_only_fields(fmt));
    let nearest_sites = side.nearest_sites.as_ref().map(|sites| sites.len());
    let json_cache_transcripts = side
        .structural_variants
        .as_ref()
        .map(|variants| variants.spans.len());
    let left_out = (fmt == "recoder").then(|| annotation_tables(side));
    let compared = compare_entries(
        &expected,
        &actual,
        documented,
        &CompareOptions {
            floors: Some(&floors),
            field_pairs: Some(field_pairs),
            check_order: fmt != "parquet" && !unordered,
            examples,
            order_only: order_only.as_ref(),
            floor_fields: floor_fields.as_ref(),
            input_lines: inputs,
            kept: kept_back.as_ref(),
            side,
        },
    );
    if let Some(taken) = left_out {
        restore_annotation_tables(side, taken);
    }
    let mut report = match compared?.into_dict() {
        Ok(report) => report,
        Err(_) => crate::uncaught!("the comparison's report is no object"),
    };
    let unused: Vec<Value> = documented
        .masks
        .iter()
        .filter(|klass| !documented.applied.contains(klass) && !documented.needed.contains(*klass))
        .map(|klass| Value::Str(klass.clone()))
        .collect();
    let mut masks: IndexMap<String, Value> = [
        (
            "requested",
            Value::List(
                documented
                    .masks
                    .iter()
                    .map(|klass| Value::Str(klass.clone()))
                    .collect(),
            ),
        ),
        ("applied", sorted_counts(&documented.applied)),
        ("unused", Value::List(unused)),
        ("available", sorted_counts(&documented.available)),
    ]
    .into_iter()
    .map(|(name, value)| (name.to_string(), value))
    .collect();
    if !documented.corrected.is_empty() {
        let corrected = documented
            .corrected
            .iter()
            .map(|(klass, release)| (klass.clone(), Value::Str(release.clone())))
            .collect();
        masks.insert("corrected".to_string(), Value::Dict(corrected));
    }
    let mut pair_names: Vec<&String> = field_pairs.keys().collect();
    pair_names.sort_unstable();
    let accepted_pairs: IndexMap<String, Value> = pair_names
        .into_iter()
        .map(|name| {
            let pairs = field_pairs[name]
                .iter()
                .map(|(reference, candidate)| {
                    Value::List(vec![
                        Value::Str(reference.clone()),
                        Value::Str(candidate.clone()),
                    ])
                })
                .collect();
            (name.clone(), Value::List(pairs))
        })
        .collect();
    let mut put = |key: &str, value: Value| {
        report.insert(key.to_string(), value);
    };
    put("mode", Value::Str("compare".to_string()));
    put("format", Value::Str(fmt.to_string()));
    put("expected", Value::Str(path_text(expected_path)));
    put("actual", Value::Str(path_text(actual_path)));
    put("masks", Value::Dict(masks));
    put("field_pairs", Value::Dict(accepted_pairs));
    put("nearest_sites", optional_count(nearest_sites));
    put(
        "json_cache_transcripts",
        optional_count(json_cache_transcripts),
    );
    put("unordered", Value::Bool(unordered));
    put("keep_only", optional_text(keep_only));
    put("manifest_run", optional_text(documented.run.as_deref()));
    put(
        "reference_sequence",
        optional_text(documented.reference_sequence.as_deref()),
    );
    Ok(Value::Dict(report))
}
