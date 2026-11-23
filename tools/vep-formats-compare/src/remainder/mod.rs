// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the front end keeps back of a comparison.
//!
//! The front end is the stream of record groups (`crate::stream`): a record group is every
//! row of one side that shares a Location (default, tab), every line of one side whose
//! first five columns are the same text (VCF) or every object of one side whose `input`
//! opens with the same five columns (JSON). `compare_entries` pairs entries by a key that
//! starts with the Location or with the label of the line's input record, groups a
//! record's rows under a label that starts with it (`record_label`) and reaches a rule only
//! through a field that differs, so a group whose lines give the same entries in the same
//! order on both sides adds its entries to the counts, its fields to the populated columns
//! and its place to the order of the shared keys, and nothing else. Three readings take an
//! equal group all the same, and each hands its groups over whole: a key the manifest
//! documents is checked for staleness on equal values too, `regulatory_feature_rows` drops
//! reference rows before the comparison, and `conservation_short_name_spelling` rewrites
//! the reference text. A VCF or JSON entry reads two things more that its own line does
//! not hold, and the front end's counts carry both: a line no input record matches is
//! labelled by its ordinal in the output (`KeptBack::ordinals`), and a JSON floor is placed
//! and counted by the members of every entry (`KeptBack::present`,
//! `KeptBack::counted_under`).
//!
//! A Parquet comparison reads the canonical form of its two outputs as two tab outputs
//! without their order (`hand_over_canonical`): its record groups are the tab format's.

use crate::classes;
use crate::classes::CONSERVATION_SHORT_NAME_SPELLING;
use crate::classes::REGULATORY_FEATURE_ROWS;
use crate::consts::JSON_CONSEQUENCE_LISTS;
use crate::consts::MISSPELT_CONSERVATION_KEY;
use crate::consts::REGULATORY_KINDS;
use crate::counter::Counter;
use crate::error::Error;
use crate::error::Result;
use crate::floors::json_floor_fields;
use crate::manifest::Documented;
use crate::pyfmt;
use crate::shared::feature_kind;
use crate::stream;
use crate::stream::Output;
use crate::stream::Text;
use indexmap::IndexMap;
use std::collections::BTreeSet;
use std::path::Path;
use std::sync::Arc;

pub const FRONT_END_FORMATS: &[&str] = &["default", "tab", "vcf", "json"];

/// The formats whose entries name an input record: the front end reads the input beside
/// them.
pub const FRONT_END_RECORD_FORMATS: &[&str] = &["vcf", "json"];

/// The Feature_type values whose rows `regulatory_feature_rows` is asked about: a reference
/// row of a value it keeps stays where it is, a row of any other value hands its group over.
pub const FRONT_END_FEATURE_TYPES: &[&str] = &["Transcript", "-", ""];

/// Set to `1`, the outputs are compared whole and no record group is kept back.
pub const WHOLE_VARIABLE: &str = "VEP_FORMATS_COMPARE_WHOLE";

/// The name the stream of record groups goes by in a `FrontEnd`'s note.
pub const STREAM: &str = "the stream of record groups";

/// The front end of a comparison by its name, and `note`, what the comparison did with
/// it.
#[derive(Debug)]
pub struct FrontEnd {
    pub executable: String,
    pub note: String,
    /// The times the outputs were read again because two record groups of one digest
    /// differed in their text.
    pub restarts: u64,
}

impl FrontEnd {
    pub fn new(executable: &str) -> Self {
        FrontEnd {
            executable: executable.to_string(),
            note: "the whole outputs were compared: only a default, tab, vcf, json or parquet comparison is handed over"
                .to_string(),
            restarts: 0,
        }
    }

    /// The front end that is the stream of record groups, or `None` where `WHOLE_VARIABLE`
    /// asks for the whole outputs.
    pub fn stream() -> Option<Self> {
        let whole = std::env::var_os(WHOLE_VARIABLE).is_some_and(|value| value == "1");
        (!whole).then(|| FrontEnd::new(STREAM))
    }
}

/// What the front end kept back of one comparison: `expected` and `actual` are the
/// remainder's lines, each side's `#` lines and the data lines handed over in file order,
/// `entries` the entries kept back on each side (a row of the default or tab format; a VCF
/// line's record entry and CSQ entries; a JSON object's record entry and consequence
/// entries), `populated` their count per field (`populated_counts`), and `runs` the
/// reference remainder's entries in file order as (entries, kept-back entries before
/// them). `keys` are the keys of those entries, set once the remainder is parsed. A VCF or
/// JSON comparison has `line_runs` as well: per side, the remainder's lines in file order
/// as (lines, lines of the output before them that the remainder lacks). A JSON comparison
/// has `present`, the members the kept-back entries carry, and `under`, per floor column
/// the fields the front end was given for it and the kept-back consequence entries
/// `populated_under` counts under them. `candidate_error` is the error reading the
/// candidate ended with beside a reference that was read: the remainder is then the whole
/// reference, and the error is the comparison's once the reference's entries are parsed.
#[derive(Debug)]
pub struct KeptBack {
    pub expected: Vec<String>,
    pub actual: Vec<String>,
    pub entries: i64,
    pub populated: Counter,
    pub runs: Vec<(i64, i64)>,
    pub keys: Vec<String>,
    pub line_runs: Option<IndexMap<String, Vec<(i64, i64)>>>,
    pub present: BTreeSet<String>,
    pub under: IndexMap<String, (Vec<String>, i64)>,
    pub candidate_error: Option<Error>,
}

impl KeptBack {
    /// The ordinal in its output of each line of the remainder of `side` (`expected` or
    /// `actual`), in order; `None` for a default or tab comparison, whose entries read no
    /// ordinal.
    pub fn ordinals(&self, side: &str) -> Option<Vec<i64>> {
        let line_runs = self.line_runs.as_ref()?;
        let Some(runs) = line_runs.get(side) else {
            crate::uncaught!("{}", pyfmt::repr_str(side))
        };
        Some(remainder_ordinals(runs))
    }

    /// The kept-back consequence entries a JSON floor on `column` counts where it is
    /// counted under `fields`: the front end's count where those are the fields it was
    /// given, none under any other (a plugin object's member, which no kept-back entry
    /// carries).
    pub fn counted_under(&self, column: &str, fields: &[String]) -> i64 {
        match self.under.get(column) {
            Some((given, count)) if given.as_slice() == fields => *count,
            _ => 0,
        }
    }

    /// The kept-back entries that precede the remainder's entry of `key` in the reference.
    pub fn before(&self, key: &str) -> i64 {
        let Some(row) = self.keys.iter().position(|held| held == key) else {
            crate::uncaught!("{} is not in list", pyfmt::repr_str(key))
        };
        let mut row = i64::try_from(row).unwrap_or(i64::MAX);
        for &(rows, before) in &self.runs {
            if row < rows {
                return before;
            }
            row -= rows;
        }
        0
    }
}

/// The error of a remainder that holds a line past the last run of its line runs.
pub const REMAINDER_LINE_UNPLACED: &str =
    "the front end's remainder holds a line its counters do not place";

/// The ordinal in its output of each line of a remainder, in order, from the runs of
/// `KeptBack::line_runs`. A line past the last run has no ordinal, and a remainder that
/// holds one is the front end's error (`REMAINDER_LINE_UNPLACED`).
pub fn remainder_ordinals(runs: &[(i64, i64)]) -> Vec<i64> {
    let mut ordinals = Vec::new();
    let mut line: i64 = 0;
    for &(lines, before) in runs {
        for _ in 0..lines {
            ordinals.push(line.saturating_add(before));
            line = line.saturating_add(1);
        }
    }
    ordinals
}

/// Whether `value` can stand on a line of a list: no tab and no line end.
pub fn on_one_line(value: &str) -> bool {
    !value.contains(['\t', '\n', '\r'])
}

/// The text every match of `pattern` starts with, read off the pattern's own text: its
/// leading letters, digits and underscores, less the last one when a quantifier follows it.
/// `None` for a pattern with an alternation or a flag, or without such a head. A flag is
/// one the pattern's text opens with, the one way a pattern of this crate carries one.
pub fn literal_head(pattern: &regex::Regex) -> Option<String> {
    let source = pattern.as_str();
    if source.contains('|') {
        return None;
    }
    let length = source
        .bytes()
        .take_while(|b| b.is_ascii_alphanumeric() || *b == b'_')
        .count();
    let mut head = &source[..length];
    let quantified = source[length..]
        .chars()
        .next()
        .is_some_and(|next| matches!(next, '*' | '+' | '?' | '{'));
    if quantified && !head.is_empty() {
        head = &head[..head.len() - 1];
    }
    (!head.is_empty()).then(|| head.to_string())
}

/// The values a list of the front end can hold, in byte order, each once: a value no line
/// of an output can hold (a tab, a line end) is left out.
fn listed<'a>(values: impl IntoIterator<Item = &'a str>) -> Vec<String> {
    let values: BTreeSet<&str> = values
        .into_iter()
        .filter(|value| on_one_line(value))
        .collect();
    values.into_iter().map(str::to_string).collect()
}

fn whole_runs(runs: &[(u64, u64)]) -> Vec<(i64, i64)> {
    runs.iter()
        .map(|(count, before)| (whole_number(*count), whole_number(*before)))
        .collect()
}

fn whole_number(value: u64) -> i64 {
    i64::try_from(value).unwrap_or(i64::MAX)
}

/// Runs the front end over one comparison and returns what it kept back, or `None` for
/// a comparison made on the whole outputs; `front_end.note` says which. The front end is
/// given the format and what marks a group the comparison must see. For a default or tab
/// comparison, the Locations of every key the manifest documents. For a VCF or JSON
/// comparison, the input `records_path` (without one the whole outputs are compared, and
/// the parser refuses them) and the ordinals of the records the manifest documents a key on
/// or names as skipped. Under `regulatory_feature_rows` the Feature_type values, in JSON
/// the consequence lists' kinds, of the rows that class keeps; under
/// `conservation_short_name_spelling` the text every misspelt key starts with
/// (`literal_head`; a pattern that shows none is compared on the whole outputs); for a
/// JSON comparison the fields each of `floors`, the floors' columns as the comparison
/// spells them, is counted under (`json_floor_fields`); and whether the order of the
/// shared keys is read. The front end hands every row over where that text stands in the
/// reference's column line or in a key column, which the rewrite would rename or move, and
/// where a Location holds a space, the separator of `record_label`. A JSON floor is
/// compared on the whole outputs where one of its fields holds a `.`, which
/// `populated_under` reads as a member of an object as well, and where a line of a list
/// cannot hold its column or a field (`on_one_line`).
///
/// `floors` is empty and `records_path` is `None` in a comparison that names neither.
#[allow(clippy::too_many_arguments)]
pub fn hand_over(
    front_end: &mut FrontEnd,
    fmt: &str,
    expected_path: &Path,
    actual_path: &Path,
    documented: &Documented,
    unordered: bool,
    keep_only: Option<&str>,
    floors: &[&str],
    records_path: Option<&Path>,
) -> Result<Option<KeptBack>> {
    hand_over_outputs(
        front_end,
        fmt,
        &Output::File(expected_path),
        &Output::File(actual_path),
        documented,
        unordered,
        keep_only,
        floors,
        records_path,
    )
}

/// `hand_over` of a Parquet comparison. `expected` and `actual` are its two outputs in
/// the canonical form with the rows in each output's own order
/// (`parse::parquet::CanonicalText`), read as two tab outputs whose order the comparison
/// does not read. What is kept back is `canonical_order`ed: the remainder's data lines
/// stand in byte order after each side's column line, the order the whole canonical form
/// has them in, so the entries parsed from the remainder are the whole comparison's
/// entries of the same groups, in its order.
pub fn hand_over_canonical(
    front_end: &mut FrontEnd,
    expected: Arc<dyn Text>,
    actual: Arc<dyn Text>,
    documented: &Documented,
    keep_only: Option<&str>,
) -> Result<Option<KeptBack>> {
    let kept = hand_over_outputs(
        front_end,
        "tab",
        &Output::Text(expected),
        &Output::Text(actual),
        documented,
        true,
        keep_only,
        &[],
        None,
    )?;
    Ok(kept.map(canonical_order))
}

/// `kept` of two canonical forms with each side's lines after its first, the column line,
/// in byte order. The runs of `KeptBack::runs` place the lines of a file's order and are
/// dropped: a comparison of canonical forms reads no order.
fn canonical_order(mut kept: KeptBack) -> KeptBack {
    for lines in [&mut kept.expected, &mut kept.actual] {
        if let Some((_, rows)) = lines.split_first_mut() {
            rows.sort_unstable();
        }
    }
    kept.runs.clear();
    kept
}

/// `hand_over` of two outputs that are each a file or a written text.
#[allow(clippy::too_many_arguments)]
fn hand_over_outputs(
    front_end: &mut FrontEnd,
    fmt: &str,
    expected: &Output,
    actual: &Output,
    documented: &Documented,
    unordered: bool,
    keep_only: Option<&str>,
    floors: &[&str],
    records_path: Option<&Path>,
) -> Result<Option<KeptBack>> {
    let mut whole: Option<String> = if !FRONT_END_FORMATS.contains(&fmt) {
        Some(format!("the {fmt} format"))
    } else if keep_only.is_some() {
        Some("--keep-only".to_string())
    } else {
        classes::tables()
            .front_end_whole_classes()
            .iter()
            .find(|k| documented.masks.contains(**k))
            .map(|k| (*k).to_string())
    };
    let by_record = FRONT_END_RECORD_FORMATS.contains(&fmt);
    if whole.is_none() && by_record && records_path.is_none() {
        whole = Some(format!("a {fmt} comparison without --input"));
    }
    let mut floor_fields: IndexMap<String, Vec<String>> = IndexMap::new();
    if whole.is_none() && fmt == "json" {
        for column in floors {
            if let Some(fields) = json_floor_fields(column) {
                floor_fields.insert((*column).to_string(), fields);
            }
        }
        let names = || {
            floor_fields
                .iter()
                .flat_map(|(column, fields)| std::iter::once(column).chain(fields))
        };
        if floor_fields
            .values()
            .flatten()
            .any(|name| name.contains('.'))
        {
            whole = Some("a floor on a field holding a dot".to_string());
        } else if !names().all(|name| on_one_line(name)) {
            whole = Some("a floor on a column no line of a list holds".to_string());
        }
    }
    let mut misspelt: Option<String> = None;
    if whole.is_none() && documented.masks.contains(CONSERVATION_SHORT_NAME_SPELLING) {
        misspelt = literal_head(&MISSPELT_CONSERVATION_KEY);
        if misspelt.is_none() {
            whole = Some(CONSERVATION_SHORT_NAME_SPELLING.to_string());
        }
    }
    if let Some(whole) = whole {
        front_end.note =
            format!("the whole outputs were compared: {whole} is outside the front end");
        return Ok(None);
    }
    let mut force_locations = Vec::new();
    let mut force_records = Vec::new();
    if by_record {
        // A set of the ordinals; one below zero names no line of an input.
        let records: BTreeSet<i64> = documented
            .by_record
            .keys()
            .map(|key| key.0)
            .chain(documented.extra_by_record.iter().map(|key| key.0))
            .chain(documented.fields_by_record.keys().map(|key| key.0))
            .chain(documented.equal_by_record.keys().map(|key| key.0))
            .chain(documented.skipped_records.iter().copied())
            .collect();
        force_records = records
            .into_iter()
            .filter_map(|ordinal| u64::try_from(ordinal).ok())
            .collect();
    } else {
        let locations = documented
            .by_location
            .keys()
            .map(|key| key.0.as_str())
            .chain(
                documented
                    .extra_by_location
                    .iter()
                    .map(|key| key.0.as_str()),
            )
            .chain(
                documented
                    .fields_by_location
                    .keys()
                    .map(|key| key.0.as_str()),
            );
        force_locations = listed(locations);
    }
    let keep_values = documented.masks.contains(REGULATORY_FEATURE_ROWS).then(|| {
        if fmt == "json" {
            let kept_kinds = JSON_CONSEQUENCE_LISTS
                .iter()
                .map(|(_, kind)| *kind)
                .filter(|kind| !REGULATORY_KINDS.contains(kind));
            return ("kind".to_string(), listed(kept_kinds));
        }
        let kept_types = FRONT_END_FEATURE_TYPES
            .iter()
            .copied()
            .filter(|t| !REGULATORY_KINDS.contains(&feature_kind(t).as_str()));
        ("Feature_type".to_string(), listed(kept_types))
    });
    let plan = stream::Plan {
        format: match fmt {
            "tab" => stream::Format::Tab,
            "vcf" => stream::Format::Vcf,
            "json" => stream::Format::Json,
            _ => stream::Format::Default,
        },
        force_locations,
        force_text: listed(misspelt.as_deref()),
        keep_values,
        unordered,
        threads: crate::cache::threads(),
        input: records_path.filter(|_| by_record).map(Path::to_path_buf),
        force_records,
        floor_fields: floor_fields
            .iter()
            .map(|(column, fields)| (column.clone(), fields.clone()))
            .collect(),
    };
    let remainder = stream::remainder_of(expected, actual, &plan)?;
    let mut kept = KeptBack {
        expected: remainder.expected,
        actual: remainder.actual,
        entries: whole_number(remainder.kept_entries),
        populated: remainder.kept_populated,
        runs: whole_runs(&remainder.kept_before),
        keys: Vec::new(),
        line_runs: None,
        present: BTreeSet::new(),
        under: IndexMap::new(),
        candidate_error: remainder.candidate_error,
    };
    if by_record {
        kept.line_runs = Some(IndexMap::from([
            (
                "expected".to_string(),
                whole_runs(&remainder.expected_line_runs),
            ),
            (
                "actual".to_string(),
                whole_runs(&remainder.actual_line_runs),
            ),
        ]));
    }
    if fmt == "json" {
        kept.present = remainder
            .kept_present
            .iter()
            .filter(|(_, count)| *count != 0)
            .map(|(name, _)| name.to_string())
            .collect();
        for (column, fields) in floor_fields {
            let count = remainder.kept_under.get(&column);
            kept.under.insert(column, (fields, count));
        }
    }
    front_end.restarts = remainder.restarts;
    let lines = remainder.expected_lines;
    let groups = format!(
        "{} of {} record groups handed over",
        remainder.groups_handed_over, remainder.groups
    );
    front_end.note = match (by_record, remainder.whole) {
        (true, Some(reason)) => format!(
            "every line was handed over by {}: {reason}",
            front_end.executable
        ),
        (true, None) => format!(
            "{} of {lines} reference lines were kept back by {}, {groups}",
            lines - remainder.expected_handed_over,
            front_end.executable,
        ),
        (false, Some(reason)) => format!(
            "every row was handed over by {}: {reason}",
            front_end.executable
        ),
        (false, None) => format!(
            "{} of {lines} reference rows were kept back by {}, {groups}",
            kept.entries, front_end.executable,
        ),
    };
    Ok(Some(kept))
}
