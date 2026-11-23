// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The JSON end key of a position range under `--total_length`: the entries a run that
//! compares the JSON format documents for it, read from the two default-format outputs.

use super::fields::FieldKey;
use crate::error::Failure;
use crate::json::{Object, Value};
use crate::rows::iter_rows;
use std::collections::BTreeMap;

/// The class of a JSON end key under `--total_length`. Ensembl VEP writes the first
/// position of a range under both keys and neither key of a range whose first position is
/// unknown: its test of the end admits no total (OutputFactory/JSON.pm:337 at
/// release/116.3). vep-rs writes the last position with its total under the end key.
pub const TOTAL_LENGTH_JSON_RANGE_END: &str = "total_length_json_range_end";

/// The flag under which a position column carries its total.
pub const TOTAL_LENGTH_FLAG: &str = "--total_length";

/// The JSON end keys of the three position columns, in the columns' order.
pub const RANGE_END_FIELDS: [&str; 3] = ["cdna_end", "cds_end", "protein_end"];

/// Whether a run documents the class: its flags hold `--total_length` and it compares the
/// JSON format, which a run does whose entry names no `formats` or names `json` among them.
pub fn documents_range_ends(flags: &[String], run: &Object) -> bool {
    if !flags.iter().any(|flag| flag == TOTAL_LENGTH_FLAG) {
        return false;
    }
    match run.get("formats").and_then(Value::as_list) {
        None => true,
        Some(formats) => formats.iter().any(|format| format.as_str() == Some("json")),
    }
}

fn digits(text: &str) -> bool {
    !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit())
}

/// What the two engines write under the end key of a position column whose text is a range
/// with a total: Ensembl VEP's value, the first position or `-` where it is `?`, and
/// vep-rs's, the last position with the total. `None` for any other text: a single
/// position, a range without a total, a range whose last position is unknown.
pub fn range_end(column: &str) -> Option<(String, String)> {
    let (first, rest) = column.split_once('-')?;
    let (last, total) = rest.split_once('/')?;
    if !digits(last) || !digits(total) {
        return None;
    }
    let reference = match first {
        "?" => "-",
        first if digits(first) => first,
        _ => return None,
    };
    Some((reference.to_string(), rest.to_string()))
}

/// The entry among `entries` that states the field `field`, its name read in any letter
/// case, on the key and consequence set of `key`: the one a further format's reading wrote
/// from the JSON outputs themselves.
pub fn stated_end_key<'a>(
    entries: &'a mut [Object],
    key: &FieldKey,
    field: &str,
) -> Option<&'a mut Object> {
    fn text<'e>(entry: &'e Object, name: &str) -> &'e str {
        entry.get(name).and_then(Value::as_str).unwrap_or_default()
    }
    entries.iter_mut().find(|entry| {
        text(entry, "location") == key.location
            && text(entry, "allele") == key.allele
            && text(entry, "feature") == key.feature
            && text(entry, "feature_type") == key.feature_type
            && text(entry, "consequence_set") == key.consequences
            && text(entry, "field").eq_ignore_ascii_case(field)
    })
}

/// (Location, Allele, Feature, Feature_type, consequence set) -> the `cDNA_position`,
/// `CDS_position` and `Protein_position` columns of the first row that carries the key,
/// leaving out the rows whose Feature_type is in `masked`.
pub fn load_position_columns(
    path: &str,
    masked: &[&str],
) -> Result<BTreeMap<FieldKey, [String; 3]>, Failure> {
    let mut columns_of: BTreeMap<FieldKey, [String; 3]> = BTreeMap::new();
    let mut rows = iter_rows(path, "x", &[] as &[&str], None);
    while let Some(row) = rows.next_view() {
        let row = row?;
        if masked.contains(&row.feature_type) {
            continue;
        }
        let key = FieldKey {
            location: row.location.to_string(),
            allele: row.allele.to_string(),
            feature: row.feature.to_string(),
            feature_type: row.feature_type.to_string(),
            consequences: row.consequence_set.to_string(),
        };
        columns_of
            .entry(key)
            .or_insert_with(|| row.positions.map(str::to_string));
    }
    Ok(columns_of)
}

#[cfg(test)]
#[path = "total_length_tests.rs"]
mod tests;
