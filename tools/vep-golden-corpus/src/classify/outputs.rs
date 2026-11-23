// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two default-format outputs of one run as `classify` reads them: the consequence
//! sets and record names of every key, the rows that carry a tuple with their positions
//! and their coding cells, the compared fields of every row, and the masks that leave rows
//! out or compare a field as a multiset.

use super::fields::{FieldKey, FieldRow};
use super::keys::PositionRow;
use crate::constants::{ORDER_ONLY_COLUMNS, REGULATORY_FEATURE_ROWS, REGULATORY_FEATURE_TYPES};
use crate::error::Failure;
use crate::json::Object;
use crate::manifest::names_under;
use crate::pyfmt::is_py_space;
use crate::rows::{iter_rows, normalize_consequence_set, opener, parse_extra, ReadError, RowView};
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// The key of a row: Location, Allele, Feature and Feature_type.
pub type TupleKey = (String, String, String, String);

/// What the rows of one output hold under each key, the keys in their order: the
/// consequence sets, or the Uploaded_variation names.
pub type TupleSets = BTreeMap<TupleKey, BTreeSet<String>>;

/// The mask classes of a run, as the golden harness reads them: the flag set's own `masks`
/// when its entry (`target`) carries the key, else the manifest's; `target` is the manifest
/// itself for the default run.
pub fn run_masks(manifest: &Object, target: &Object) -> Vec<String> {
    if target.contains_key("masks") {
        names_under(target, "masks")
    } else {
        names_under(manifest, "masks")
    }
}

/// The `Feature_type` values a run's `masks` leave out of the reference: the regulatory
/// rows' under `regulatory_feature_rows`; `order_only` leaves out no row (`fields_agree`
/// gives it its meaning).
pub fn masked_feature_types<S: AsRef<str>>(masks: &[S]) -> &'static [&'static str] {
    if masks
        .iter()
        .any(|mask| mask.as_ref() == REGULATORY_FEATURE_ROWS)
    {
        &REGULATORY_FEATURE_TYPES
    } else {
        &[]
    }
}

/// Whether two values of the Extra field `name` agree: exactly, or, for a column of
/// `ORDER_ONLY_COLUMNS` under `order_only`, as a multiset of their comma-separated items;
/// an absent value (`-`) agrees with itself alone.
pub fn fields_agree(name: &str, vep: &str, vep_rs: &str, order_only: bool) -> bool {
    if order_only && ORDER_ONLY_COLUMNS.contains(&name) && vep != "-" && vep_rs != "-" {
        let sorted = |value: &str| {
            let mut items: Vec<String> = value.split(',').map(str::to_string).collect();
            items.sort_unstable();
            items
        };
        return sorted(vep) == sorted(vep_rs);
    }
    vep == vep_rs
}

fn tuple_key(row: &RowView<'_>) -> TupleKey {
    (
        row.location.to_string(),
        row.allele.to_string(),
        row.feature.to_string(),
        row.feature_type.to_string(),
    )
}

/// (Location, Allele, Feature, Feature_type) -> consequence-set strings, and ->
/// Uploaded_variation names, leaving out the rows whose Feature_type is in `masked`.
pub fn load_tuples(path: &str, masked: &[&str]) -> Result<(TupleSets, TupleSets), Failure> {
    let mut sets = TupleSets::new();
    let mut names = TupleSets::new();
    let mut rows = iter_rows(path, "x", &[] as &[&str], None);
    while let Some(row) = rows.next_view() {
        let row = row?;
        if masked.contains(&row.feature_type) {
            continue;
        }
        let key = tuple_key(&row);
        sets.entry(key.clone())
            .or_default()
            .insert(row.consequence_set.to_string());
        names
            .entry(key)
            .or_default()
            .insert(row.uploaded_variation.to_string());
    }
    Ok((sets, names))
}

/// (Location, Allele, Feature, Feature_type, consequence set) -> every distinct row that
/// carries the tuple, in the file's order, as its `cDNA_position` and `CDS_position` columns,
/// its `Uploaded_variation` and, where it has all twelve columns up to `Codons`, its
/// `Protein_position`, `Amino_acids` and `Codons`, leaving out the rows whose Feature_type is
/// in `masked`. A position column a row ends before is `-`.
pub fn load_positions(
    path: &str,
    masked: &[&str],
) -> Result<HashMap<FieldKey, Vec<PositionRow>>, Failure> {
    let mut rows_of: HashMap<FieldKey, Vec<PositionRow>> = HashMap::new();
    let mut rows = iter_rows(path, "x", &[] as &[&str], None);
    while let Some(row) = rows.next_view() {
        let row = row?;
        if masked.contains(&row.feature_type) {
            continue;
        }
        let positions = PositionRow {
            cells: row.cells.map(|cells| cells.map(str::to_string)),
            ..PositionRow::of_positions(row.positions[0], row.positions[1], row.uploaded_variation)
        };
        let key = FieldKey {
            location: row.location.to_string(),
            allele: row.allele.to_string(),
            feature: row.feature.to_string(),
            feature_type: row.feature_type.to_string(),
            consequences: row.consequence_set.to_string(),
        };
        let carried = rows_of.entry(key).or_default();
        if !carried.contains(&positions) {
            carried.push(positions);
        }
    }
    Ok(rows_of)
}

/// The columns of a default-format row that `--field` compares under their own names, with
/// their places: the two cells of the class `alternate_codon_past_the_coding_end`. A field
/// of `load_fields` named as one of them holds that column of the row, never an Extra value.
pub const COMPARED_COLUMNS: [(&str, usize); 2] = [("Amino_acids", 10), ("Codons", 11)];

/// (Location, Allele, Feature, Feature_type, consequence set) -> the named Extra fields (`-`
/// when absent) and the row's `cDNA_position`, `CDS_position`, `Protein_position` and
/// `Amino_acids` columns, leaving out the rows whose Feature_type is in `masked`. The first
/// row of a key is the one read. A field named `cDNA_position`, `CDS_position` or
/// `Protein_position` holds the Extra's value, and one named `Amino_acids` or `Codons` the
/// row's column.
pub fn load_fields<S: AsRef<str>>(
    path: &str,
    fields: &[S],
    masked: &[&str],
) -> Result<BTreeMap<FieldKey, FieldRow>, Failure> {
    let rows = read_field_rows(path, fields, masked, false)?;
    let first = |(key, rows): (FieldKey, Vec<FieldRow>)| Some((key, rows.into_iter().next()?));
    Ok(rows.into_iter().filter_map(first).collect())
}

/// `load_fields` with every row of a key, in the file's order: several records that share a
/// Location and an Allele write a row each on one key.
pub fn load_field_rows<S: AsRef<str>>(
    path: &str,
    fields: &[S],
    masked: &[&str],
) -> Result<BTreeMap<FieldKey, Vec<FieldRow>>, Failure> {
    read_field_rows(path, fields, masked, true)
}

/// The rows `load_fields` and `load_field_rows` read: every row of a key under `every_row`,
/// else its first alone.
fn read_field_rows<S: AsRef<str>>(
    path: &str,
    fields: &[S],
    masked: &[&str],
    every_row: bool,
) -> Result<BTreeMap<FieldKey, Vec<FieldRow>>, Failure> {
    let unreadable = |error: ReadError| {
        Failure::own(format!(
            "{path} is not a readable text file ({})",
            error.reason()
        ))
    };
    let mut lines = opener(path).map_err(unreadable)?;
    let mut out: BTreeMap<FieldKey, Vec<FieldRow>> = BTreeMap::new();
    while lines.advance().map_err(unreadable)? {
        let line = lines.line();
        if line.starts_with('#') {
            continue;
        }
        let mut columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 7 {
            continue;
        }
        columns.resize(columns.len().max(14), "-");
        if masked.contains(&columns[5]) {
            continue;
        }
        let key = FieldKey {
            location: columns[1].to_string(),
            allele: columns[2].to_string(),
            feature: columns[4].to_string(),
            feature_type: columns[5].to_string(),
            consequences: normalize_consequence_set(columns[6]),
        };
        if !every_row && out.contains_key(&key) {
            continue;
        }
        let extra = parse_extra(columns[13]);
        let mut row = FieldRow::new();
        for name in fields {
            let name = name.as_ref();
            let value = match COMPARED_COLUMNS.iter().find(|(column, _)| *column == name) {
                Some((_, place)) => columns[*place],
                None => extra.get(name).unwrap_or("-"),
            };
            row.insert(name.to_string(), value.to_string());
        }
        for (name, column) in [
            ("cDNA_position", 7),
            ("CDS_position", 8),
            ("Protein_position", 9),
            ("Amino_acids", 10),
        ] {
            row.entry(name.to_string())
                .or_insert_with(|| columns[column].to_string());
        }
        out.entry(key).or_default().push(row);
    }
    Ok(out)
}

/// What the rows of one key and consequence set state of one compared field: the place of
/// the pair of rows whose two values an entry states, and the count of the pairs beside it
/// whose two values agree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StatedPair {
    pub place: usize,
    pub rows_alike: usize,
}

/// The pair that states the field `name` on a key whose rows are `perl_rows` in Ensembl
/// VEP's output and `rust_rows` in vep-rs's, each in its file's order; `None` for a key
/// that states no difference of the field. The rows are paired in their order, the first
/// with the first, and the first pair whose two values differ (`fields_agree`) states the
/// field, with the count of the pairs that agree. Where the two outputs hold different
/// numbers of rows on the key, or two pairs that differ hold other values than each other,
/// the first row of each output is read alone: it states the field where its two values
/// differ, and no pair is counted as agreeing.
pub fn stated_pair(
    name: &str,
    perl_rows: &[FieldRow],
    rust_rows: &[FieldRow],
    order_only: bool,
) -> Option<StatedPair> {
    fn value<'r>(rows: &'r [FieldRow], place: usize, name: &str) -> &'r str {
        rows[place].get(name).map_or("-", String::as_str)
    }
    let values = |place: usize| (value(perl_rows, place, name), value(rust_rows, place, name));
    let differs = |place: usize| {
        let (vep, vep_rs) = values(place);
        !fields_agree(name, vep, vep_rs, order_only)
    };
    let first_alone = || {
        let held = !perl_rows.is_empty() && !rust_rows.is_empty();
        (held && differs(0)).then_some(StatedPair {
            place: 0,
            rows_alike: 0,
        })
    };
    if perl_rows.len() != rust_rows.len() {
        return first_alone();
    }
    let differing: Vec<usize> = (0..perl_rows.len())
        .filter(|place| differs(*place))
        .collect();
    let first = *differing.first()?;
    if differing
        .iter()
        .any(|place| values(*place) != values(first))
    {
        return first_alone();
    }
    Some(StatedPair {
        place: first,
        rows_alike: perl_rows.len() - differing.len(),
    })
}

/// `float(text)` for a text in ASCII: the blanks `float()` drops around the number, a sign,
/// digits with single underscores between them, a point and an exponent, or `inf`,
/// `infinity` and `nan` in any case. `None` where `float()` raises `ValueError`, and for a
/// digit outside ASCII, which `float()` reads.
pub fn py_float(text: &str) -> Option<f64> {
    let body = text.trim_matches(|c| is_py_space(c) && !('\u{1c}'..='\u{1f}').contains(&c));
    if !body.is_ascii() {
        return None;
    }
    // An underscore stands between two digits and nowhere else.
    let bytes = body.as_bytes();
    let mut digits = String::with_capacity(body.len());
    for (at, byte) in bytes.iter().enumerate() {
        if *byte != b'_' {
            digits.push(char::from(*byte));
            continue;
        }
        let between = at > 0
            && bytes[at - 1].is_ascii_digit()
            && bytes.get(at + 1).is_some_and(u8::is_ascii_digit);
        if !between {
            return None;
        }
    }
    let unsigned = digits.strip_prefix(['+', '-']).unwrap_or(&digits);
    let negative = digits.starts_with('-');
    let word = unsigned.to_ascii_lowercase();
    if word == "inf" || word == "infinity" {
        return Some(if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }
    if word == "nan" {
        return Some(f64::NAN);
    }
    let (mantissa, exponent) = match unsigned.find(['e', 'E']) {
        Some(at) => (&unsigned[..at], Some(&unsigned[at + 1..])),
        None => (unsigned, None),
    };
    let all_digits = |text: &str| text.bytes().all(|byte| byte.is_ascii_digit());
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    if !all_digits(whole) || !all_digits(fraction) || (whole.is_empty() && fraction.is_empty()) {
        return None;
    }
    if let Some(exponent) = exponent {
        let power = exponent.strip_prefix(['+', '-']).unwrap_or(exponent);
        if power.is_empty() || !all_digits(power) {
            return None;
        }
    }
    digits.parse().ok()
}

#[cfg(test)]
#[path = "outputs_tests.rs"]
mod tests;
