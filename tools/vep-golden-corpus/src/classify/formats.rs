// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The further formats of a run as `classify` reads them: the tab, VCF and JSON outputs of
//! both engines, each row under the names its format gives the fields, so that a field of a
//! documented key that differs under another name than the default format's is written
//! under that name.
//!
//! A row is found as the golden harness finds it: a tab row by its Location, Allele, Feature,
//! Feature_type and consequence set, a VCF entry and a JSON consequence object by the
//! ordinal of the record they stand on, their Allele, their Feature and their consequence
//! set. An output names its record by the first five columns of the input line, and the VCF
//! output of an input that is no VCF, which carries no input line, by its place (`Pairing`).
//! A value is written in the default format's spelling of absence and of a comma: `-`
//! for an empty VCF value and for a member a JSON object leaves out, `,` for the `&` the
//! VCF format writes in its place (OutputFactory/VCF.pm:401).

use super::fields::FieldKey;
use super::outputs::fields_agree;
use super::records::{first_five_columns, RunRecord};
use crate::constants::ORDER_ONLY_COLUMNS;
use crate::error::Failure;
use crate::json::{dumps, loads, Style, Value};
use crate::rows::{normalize_consequence_set, opener, ReadError};
use std::collections::{BTreeMap, HashMap, VecDeque};

/// The value of a field a row lacks.
pub const ABSENT: &str = "-";

/// The INFO key the VCF format writes its entries under where the run names no other.
pub const DEFAULT_INFO_FIELD: &str = "CSQ";

/// The flag that names another INFO key.
const INFO_FIELD_FLAG: &str = "--vcf_info_field";

/// The columns the golden harness compares as multisets beside `ORDER_ONLY_COLUMNS`, which
/// the default format's comparison reads: the first by its comma-separated items, the second
/// by its comma-separated groups, each a multiset of `:`-separated pairs.
const NEAREST_EXON_COLUMN: &str = "NearestExonJB";
const PAIR_GROUP_COLUMN: &str = "5UTR_annotation";

/// The column of the tab format that holds a JSON object.
const VRS_COLUMN: &str = "GA4GH_VRS";

/// The arrays of a JSON record object whose items are consequence objects, each with the
/// member that names the item's feature; an intergenic item names none.
const JSON_CONSEQUENCE_ARRAYS: [(&str, Option<&str>); 4] = [
    ("transcript_consequences", Some("transcript_id")),
    (
        "regulatory_feature_consequences",
        Some("regulatory_feature_id"),
    ),
    ("motif_feature_consequences", Some("motif_feature_id")),
    ("intergenic_consequences", None),
];

/// A row of a further format: its fields by the names the format gives them, in the byte
/// order of the names.
pub type FormatRow = BTreeMap<String, String>;

/// The key of a row whose format prints no Location: the ordinal of its record, its Allele,
/// its Feature and its consequence set.
type RecordKey = (usize, String, String, String);

/// The rows of one output of a further format.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct FormatRows {
    by_location: HashMap<FieldKey, FormatRow>,
    by_record: HashMap<RecordKey, FormatRow>,
}

impl FormatRows {
    /// The row of `key`, whose records have the ordinals `indices`: the row under the key's
    /// Location where the format prints one, else the row of the first of the records that
    /// carries the key's Allele, Feature and consequence set. A tab output that names no
    /// Feature_type column holds its rows under an empty one, and the key finds its row
    /// there: the golden harness finds a documented row without that column.
    pub fn row(&self, key: &FieldKey, indices: &[usize]) -> Option<&FormatRow> {
        if let Some(row) = self.by_location.get(key) {
            return Some(row);
        }
        if !self.by_location.is_empty() {
            let untyped = FieldKey {
                feature_type: String::new(),
                ..key.clone()
            };
            return self.by_location.get(&untyped);
        }
        indices.iter().find_map(|index| {
            let key = (
                *index,
                key.allele.clone(),
                key.feature.clone(),
                key.consequences.clone(),
            );
            self.by_record.get(&key)
        })
    }
}

/// One further format of a run: its name and the rows of the two outputs.
#[derive(Debug)]
pub struct Further {
    pub name: &'static str,
    pub vep: FormatRows,
    pub vep_rs: FormatRows,
}

/// The two outputs of one further format as the command line names them: the format, the
/// reference's file and vep-rs's.
pub type Named<'a> = (&'static str, Option<&'a str>, Option<&'a str>);

/// The refusal of a pair of which one output alone is named.
pub fn unpaired(named: &[Named<'_>]) -> Option<Failure> {
    let (name, _, _) = named
        .iter()
        .find(|(_, vep, vep_rs)| vep.is_some() != vep_rs.is_some())?;
    Some(Failure::own(format!(
        "--vep-{name} and --vep-rs-{name} name the two outputs of one format: give both or neither"
    )))
}

/// How the lines of a VCF or JSON output are paired with the records of a run.
#[derive(Debug)]
pub struct Pairing<'a> {
    /// For each record, in order, the input line whose first five columns an output writes
    /// back.
    lines: Vec<&'a str>,
    /// Whether a VCF output carries no input line, so that its lines are the records in
    /// turn.
    in_turn: bool,
}

impl<'a> Pairing<'a> {
    /// The pairing of a run with the records `records`: by the records' own lines, or, for
    /// a run over `variants.vcf` re-encoded in another format (`reencoded`, its input's
    /// lines), by those lines for a JSON output and in turn for a VCF output, as the golden
    /// harness pairs the outputs of an input that is no VCF.
    pub fn of(records: &'a [RunRecord], reencoded: Option<&'a [String]>) -> Pairing<'a> {
        match reencoded {
            Some(lines) => Pairing {
                lines: lines.iter().map(String::as_str).collect(),
                in_turn: true,
            },
            None => Pairing {
                lines: records.iter().map(|record| record.vcf.as_str()).collect(),
                in_turn: false,
            },
        }
    }
}

/// The further formats the command line names both outputs of, in the order given, read
/// under the run's `pairing` and `flags`.
pub fn load(
    named: &[Named<'_>],
    pairing: &Pairing<'_>,
    flags: &[String],
) -> Result<Vec<Further>, Failure> {
    let info_field = info_field(flags);
    let mut further = Vec::new();
    for (name, vep, vep_rs) in named {
        let (Some(vep), Some(vep_rs)) = (vep, vep_rs) else {
            continue;
        };
        let read = |path: &str| match *name {
            "tab" => tab_rows(path),
            "vcf" => vcf_rows(path, pairing, info_field),
            _ => json_rows(path, pairing),
        };
        further.push(Further {
            name,
            vep: read(vep)?,
            vep_rs: read(vep_rs)?,
        });
    }
    Ok(further)
}

/// The INFO key of the run's VCF output: the value of `--vcf_info_field` among `flags`,
/// written as one word or as two, else `CSQ`.
pub fn info_field(flags: &[String]) -> &str {
    let mut named = DEFAULT_INFO_FIELD;
    let mut words = flags.iter();
    while let Some(word) = words.next() {
        if word == INFO_FIELD_FLAG {
            named = words.next().map_or(named, String::as_str);
        } else if let Some(value) = word.strip_prefix(INFO_FIELD_FLAG) {
            named = value.strip_prefix('=').unwrap_or(named);
        }
    }
    named
}

/// Whether the two values a further format gives the field `name` agree: as the values of
/// the default format do (`fields_agree`), a column the masks compare as a multiset being
/// read under its name in any letter case, as the JSON format writes it. Under `order_only`
/// the golden harness compares two more columns of a further format as multisets: the items
/// of `NearestExonJB`, in any letter case, and the groups of `5UTR_annotation`, each group
/// by its pairs.
pub fn agree(name: &str, vep: &str, vep_rs: &str, order_only: bool) -> bool {
    let column = ORDER_ONLY_COLUMNS
        .iter()
        .find(|column| column.eq_ignore_ascii_case(name));
    if fields_agree(column.unwrap_or(&name), vep, vep_rs, order_only) {
        return true;
    }
    if !order_only || vep == ABSENT || vep_rs == ABSENT {
        return false;
    }
    let sorted = |value: &str, separator: char| {
        let mut items: Vec<String> = value.split(separator).map(str::to_string).collect();
        items.sort_unstable();
        items
    };
    if name.eq_ignore_ascii_case(NEAREST_EXON_COLUMN) {
        return sorted(vep, ',') == sorted(vep_rs, ',');
    }
    let groups = |value: &str| {
        let mut groups: Vec<Vec<String>> =
            value.split(',').map(|group| sorted(group, ':')).collect();
        groups.sort_unstable();
        groups
    };
    name == PAIR_GROUP_COLUMN && groups(vep) == groups(vep_rs)
}

/// Whether a value can stand in one cell of a line of tab-separated cells.
pub fn one_cell(value: &str) -> bool {
    !value.contains(['\t', '\n', '\r'])
}

/// A value as the default format spells an absent one.
fn dash(value: &str) -> &str {
    if value.is_empty() {
        ABSENT
    } else {
        value
    }
}

fn unreadable(path: &str, error: &ReadError) -> Failure {
    Failure::own(format!(
        "{path} is not a readable text file ({})",
        error.reason()
    ))
}

/// The input records by the columns an output writes back unchanged, each identity with the
/// ordinals of its records in order: an output line takes the next record of its identity,
/// so that identical input lines are matched in the order they stand in.
struct RecordIndex<'a> {
    unmatched: HashMap<&'a str, VecDeque<usize>>,
}

impl<'a> RecordIndex<'a> {
    fn new(pairing: &Pairing<'a>) -> RecordIndex<'a> {
        let mut unmatched: HashMap<&str, VecDeque<usize>> = HashMap::new();
        for (ordinal, line) in pairing.lines.iter().enumerate() {
            let identity = first_five_columns(line);
            unmatched.entry(identity).or_default().push_back(ordinal);
        }
        RecordIndex { unmatched }
    }

    /// The ordinal of the next record with the identity of `line`, `None` for a line of no
    /// record.
    fn record(&mut self, line: &str) -> Option<usize> {
        self.unmatched
            .get_mut(first_five_columns(line))?
            .pop_front()
    }
}

/// The rows of a tab-format output: every column under the name the `#` line gives it, the
/// first row of a key being the one read. A `GA4GH_VRS` column that holds JSON is held as
/// `compact` writes it, the golden harness comparing the object and not its text. A row
/// before the column line, or of another number of columns than that line names, is
/// refused.
pub fn tab_rows(path: &str) -> Result<FormatRows, Failure> {
    let mut lines = opener(path).map_err(|error| unreadable(path, &error))?;
    let mut names: Option<Vec<String>> = None;
    let mut rows = FormatRows::default();
    while lines.advance().map_err(|error| unreadable(path, &error))? {
        let line = lines.line();
        if let Some(header) = line.strip_prefix('#') {
            if !header.starts_with('#') {
                names = Some(header.split('\t').map(str::to_string).collect());
            }
            continue;
        }
        let Some(names) = &names else {
            return Err(Failure::own(format!(
                "{path} holds a row before the line that names its columns"
            )));
        };
        let values: Vec<&str> = line.split('\t').collect();
        if values.len() != names.len() {
            return Err(Failure::own(format!(
                "{path} holds a row of {} columns under a line that names {}",
                values.len(),
                names.len()
            )));
        }
        let mut row: FormatRow = names
            .iter()
            .cloned()
            .zip(values.iter().map(|value| value.to_string()))
            .collect();
        // Ensembl VEP writes the members of this object in an order that differs from row
        // to row.
        if let Some(object) = row.get_mut(VRS_COLUMN) {
            if let Ok(value) = loads(object) {
                *object = compact(&value);
            }
        }
        let column = |name: &str| row.get(name).map_or("", String::as_str).to_string();
        let key = FieldKey {
            location: column("Location"),
            allele: column("Allele"),
            feature: column("Feature"),
            feature_type: column("Feature_type"),
            consequences: normalize_consequence_set(&column("Consequence")),
        };
        rows.by_location.entry(key).or_insert(row);
    }
    Ok(rows)
}

/// The rows of a VCF-format output: one for each entry under `info_field` of a line that is
/// one of the run's records, its fields named by the last header line of that key, the first
/// entry of a key being the one read. A line of no record is not read. Where the lines are
/// the records in turn (`Pairing`), an output of another number of lines than the run has
/// records is refused. A line of fewer than eight columns, an entry before the header line
/// that names its fields and an entry of another number of fields than that line names are
/// refused.
pub fn vcf_rows(
    path: &str,
    pairing: &Pairing<'_>,
    info_field: &str,
) -> Result<FormatRows, Failure> {
    let mut lines = opener(path).map_err(|error| unreadable(path, &error))?;
    let header = format!("##INFO=<ID={info_field},");
    let entries_key = format!("{info_field}=");
    let mut index = RecordIndex::new(pairing);
    let mut names: Vec<String> = Vec::new();
    let mut rows = FormatRows::default();
    let mut place = 0usize;
    while lines.advance().map_err(|error| unreadable(path, &error))? {
        let line = lines.line();
        if line.starts_with(&header) {
            let format = line.split("Format: ").nth(1).unwrap_or_default();
            let format = format.trim_end_matches("\">");
            names = format.split('|').map(str::to_string).collect();
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 8 {
            return Err(Failure::own(format!(
                "{path} holds a line of fewer than eight columns"
            )));
        }
        place += 1;
        let record = if pairing.in_turn {
            Some(place - 1)
        } else {
            index.record(line)
        };
        let Some(record) = record else {
            continue;
        };
        let entries = columns[7]
            .split(';')
            .rfind(|item| item.starts_with(&entries_key));
        let Some(entries) = entries else {
            continue;
        };
        for entry in entries[entries_key.len()..].split(',') {
            let values: Vec<&str> = entry.split('|').collect();
            if values.len() != names.len() {
                return Err(Failure::own(format!(
                    "{path} holds an entry of {} fields under a header that names {} for {info_field}",
                    values.len(),
                    names.len()
                )));
            }
            let row: FormatRow = names
                .iter()
                .cloned()
                .zip(values.iter().map(|value| dash(value).replace('&', ",")))
                .collect();
            let field = |name: &str| row.get(name).map_or(ABSENT, String::as_str).to_string();
            let key = (
                record,
                field("Allele"),
                field("Feature"),
                normalize_consequence_set(&field("Consequence")),
            );
            rows.by_record.entry(key).or_insert(row);
        }
    }
    if pairing.in_turn && place != pairing.lines.len() {
        return Err(Failure::own(format!(
            "{path} holds {place} lines for the {} records of an input that is no VCF, whose lines are the records in turn",
            pairing.lines.len()
        )));
    }
    Ok(rows)
}

/// `json.dumps(value, separators=(",", ":"), ensure_ascii=False)` of a text, which is the
/// text of a JSON string as the golden harness prints one.
fn quoted(text: &str) -> String {
    dumps(&Value::Str(text.to_string()), &Style::COMPACT_UNICODE)
}

/// A float as the golden harness prints one: the shortest digits that read back as the
/// value, a whole number with `.0`, and an exponent, written with its sign, where the
/// decimal point would stand more than five places before the digits or more than sixteen
/// places into them.
pub fn float_text(value: f64) -> String {
    let scientific = format!("{value:e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let sign = if mantissa.starts_with('-') { "-" } else { "" };
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    let length = digits.len() as i32;
    // The decimal point stands behind the first `point` digits.
    let point = exponent.parse::<i32>().unwrap_or(0) + 1;
    let zeros = |count: i32| "0".repeat(count as usize);
    let body = if length <= point && point <= 16 {
        format!("{digits}{}.0", zeros(point - length))
    } else if 0 < point && point <= 16 {
        let (whole, fraction) = digits.split_at(point as usize);
        format!("{whole}.{fraction}")
    } else if -5 < point && point <= 0 {
        format!("0.{}{digits}", zeros(-point))
    } else if length == 1 {
        format!("{digits}e{:+}", point - 1)
    } else {
        let (first, rest) = digits.split_at(1);
        format!("{first}.{rest}e{:+}", point - 1)
    };
    format!("{sign}{body}")
}

/// A JSON value as the golden harness prints one inside an object or a list: no blank, the
/// members of an object in the byte order of their names, a whole number that 64 bits hold
/// as it is written and any other number as `float_text` writes it.
pub fn compact(value: &Value) -> String {
    match value {
        Value::Null => "null".to_string(),
        Value::Bool(held) => held.to_string(),
        Value::Int(number) => {
            let text = number.text();
            if text.parse::<i64>().is_ok() || text.parse::<u64>().is_ok() {
                text.to_string()
            } else {
                float_text(text.parse().unwrap_or(f64::NAN))
            }
        }
        Value::Float(number) => float_text(*number),
        Value::Str(text) => quoted(text),
        Value::List(items) => {
            let items: Vec<String> = items.iter().map(compact).collect();
            format!("[{}]", items.join(","))
        }
        Value::Object(object) => {
            let mut members: Vec<(&str, &Value)> = object.iter().collect();
            members.sort_by_key(|(name, _)| *name);
            let members: Vec<String> = members
                .into_iter()
                .map(|(name, value)| format!("{}:{}", quoted(name), compact(value)))
                .collect();
            format!("{{{}}}", members.join(","))
        }
    }
}

/// The value of a member of a consequence object as the golden harness compares it: a text
/// as it is, a list as its items joined by `,`, anything else as `compact` writes it.
pub fn member_text(value: &Value) -> String {
    match value {
        Value::Str(text) => text.clone(),
        Value::List(items) => {
            let items: Vec<String> = items.iter().map(member_text).collect();
            items.join(",")
        }
        other => compact(other),
    }
}

/// The rows of a JSON-format output: one for each item of the consequence arrays of an
/// object whose `input` is the line of one of the run's records, its members by name, the first item of a key
/// being the one read. An object of no record is not read. A line that is no JSON object,
/// and an item of a consequence array that is no object, are refused.
pub fn json_rows(path: &str, pairing: &Pairing<'_>) -> Result<FormatRows, Failure> {
    let mut lines = opener(path).map_err(|error| unreadable(path, &error))?;
    let mut index = RecordIndex::new(pairing);
    let mut rows = FormatRows::default();
    let mut place = 0usize;
    while lines.advance().map_err(|error| unreadable(path, &error))? {
        place += 1;
        let line = lines.line();
        if line.trim().is_empty() {
            continue;
        }
        let refused = |what: &str| Failure::own(format!("{path} line {place} {what}"));
        let Ok(Value::Object(object)) = loads(line) else {
            return Err(refused("is no JSON object"));
        };
        let input = object.get("input").and_then(Value::as_str);
        let Some(record) = index.record(input.unwrap_or_default()) else {
            continue;
        };
        for (array, feature_member) in JSON_CONSEQUENCE_ARRAYS {
            let items = object.get(array).and_then(Value::as_list);
            for item in items.unwrap_or_default() {
                let Some(item) = item.as_object() else {
                    return Err(refused(&format!("holds under {array} what is no object")));
                };
                let row: FormatRow = item
                    .iter()
                    .map(|(name, value)| (name.to_string(), member_text(value)))
                    .collect();
                let member = |name: &str| row.get(name).map_or("", String::as_str);
                let key = (
                    record,
                    dash(member("variant_allele")).to_string(),
                    dash(feature_member.map_or("", member)).to_string(),
                    normalize_consequence_set(member("consequence_terms")),
                );
                rows.by_record.entry(key).or_insert(row);
            }
        }
    }
    Ok(rows)
}

#[cfg(test)]
#[path = "formats_tests.rs"]
mod tests;
