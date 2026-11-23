// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The sibling-allele classes of AlphaMissense, REVEL and SpliceVault: a plugin that matches
//! its table against every alternate allele of a record and keeps the first row that matches
//! writes one allele's values on another allele's row.

use super::entry::{blank, field_named, field_of, Entry, Row};
use super::shared::sibling_rows;
use std::collections::BTreeMap;

pub const AMINO_ACIDS_FIELDS: [&str; 1] = ["Amino_acids"];

/// The entry's filled fields among those `names` name (`field_named`), by name.
pub fn class_field_values<'a>(e: &'a Entry, names: &[&str]) -> BTreeMap<&'a str, &'a str> {
    e.fields
        .iter()
        .filter(|(name, value)| field_named(name, names) && !blank(value))
        .collect()
}

/// Whether the reference's row `e` carries, in every field `names` name, the values its row
/// of another allele of the record carries on the same feature and, under
/// `same_amino_acids`, beside the same `Amino_acids`. The row's fields are read together,
/// and the sibling counts only where the engine's row of that allele (`candidate_rows`, the
/// engine's rows of the same record) carries the values too: the named allele's own row has
/// no such sibling, and without the engine's rows no row is the shape.
pub fn sibling_allele_value(
    e: &Entry,
    names: &[&str],
    record_rows: &[Row],
    same_amino_acids: bool,
    candidate_rows: Option<&[Row]>,
) -> bool {
    let values = class_field_values(e, names);
    let amino_acids = field_of(e, &AMINO_ACIDS_FIELDS);
    let Some(candidate_rows) = candidate_rows else {
        return false;
    };
    if values.is_empty() || (same_amino_acids && blank(amino_acids)) {
        return false;
    }
    let named: Vec<Row> = sibling_rows(e, candidate_rows)
        .into_iter()
        .filter(|c| class_field_values(c, names) == values)
        .collect();
    sibling_rows(e, record_rows).iter().any(|s| {
        class_field_values(s, names) == values
            && (!same_amino_acids || field_of(s, &AMINO_ACIDS_FIELDS) == amino_acids)
            && named.iter().any(|c| c.allele == s.allele)
    })
}
