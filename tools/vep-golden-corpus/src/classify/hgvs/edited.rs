// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The transcripts of the corpus cache whose model is edited against the genome by an
//! insertion or a deletion, and the two classes of a field on such a transcript.

use super::notation::{blank, body_integers, coordinates_masked, non_zero_integer};
use crate::classify::cache::iter_cache_transcripts;
use crate::error::Failure;
use crate::json::{Object, Value};
use crate::pytext::{py_int, py_str};
use std::collections::{BTreeSet, HashMap};

/// The distances the edits of each edited transcript put between the edited and the
/// unedited cDNA numbering, by stable id.
pub type EditedModels = HashMap<String, BTreeSet<i128>>;

/// Whether `text` is an integer as an `_rna_edit` value writes a coordinate (`[+-]?[0-9]+`).
fn signed_digits(text: &str) -> bool {
    let digits = text.strip_prefix(['+', '-']).unwrap_or(text);
    !digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit())
}

/// The start and the length change of the edit one attribute states, `None` for an
/// attribute that is no `_rna_edit`, for a substitution and for a value without two leading
/// integers.
fn attribute_edit(attribute: &Object, stable_id: &str) -> Result<Option<(i128, i128)>, Failure> {
    let text = |key: &str| attribute.get(key).map(py_str).unwrap_or_default();
    if !text("code").starts_with("_rna_edit") {
        return Ok(None);
    }
    let value = text("value");
    let mut fields: Vec<&str> = value.split(' ').collect();
    while fields.last().is_some_and(|field| field.is_empty()) {
        fields.pop();
    }
    if fields.len() < 2 || !fields[..2].iter().all(|field| signed_digits(field)) {
        return Ok(None);
    }
    let coordinate = |field: &str| {
        py_int(field).map(i128::from).ok_or_else(|| {
            Failure::own(format!(
                "the cache record {stable_id} holds an _rna_edit coordinate 64 bits do not hold"
            ))
        })
    };
    let (start, end) = (coordinate(fields[0])?, coordinate(fields[1])?);
    if fields.len() == 3 {
        let bases = i128::try_from(fields[2].chars().count()).unwrap_or(i128::MAX);
        if end - start + 1 == bases {
            return Ok(None);
        }
        return Ok(Some((start, bases)));
    }
    Ok(Some((start, -1 - (end - start))))
}

/// The cache's transcripts carrying an `_rna_edit` attribute that inserts or deletes bases
/// against the genome, by stable id, each with the distances its edits put between the
/// edited and the unedited cDNA numbering: the running sums of the edits' length changes in
/// start order, as absolute values, 0 left out. The value is read as Ensembl VEP's
/// `get_misalignment_offset` reads it (TranscriptVariationAllele.pm:1605-1618 at
/// release/116): `<start> <end> <bases>` split on single spaces with trailing empty fields
/// dropped, a three-field value whose bases span exactly `end - start + 1` a substitution,
/// any other three-field value an insertion of its bases and a value of any other count of
/// fields a deletion of `end - start + 1` bases; a value without two leading integers is
/// skipped. Empty when the cache is absent.
pub fn indel_edited_models(cache_dir: &str) -> Result<EditedModels, Failure> {
    let mut edited = EditedModels::new();
    for item in iter_cache_transcripts(cache_dir) {
        let (_chromosome, transcript) = item?;
        let stable_id = transcript.get("stable_id").map(py_str).unwrap_or_default();
        let attributes: &[Value] = match transcript.get("attributes") {
            Some(Value::List(attributes)) => attributes.as_slice(),
            // The keys of an object and the characters of a text are no attributes.
            Some(Value::Object(_) | Value::Str(_)) | None => &[],
            Some(other) if !other.truthy() => &[],
            Some(_) => {
                return Err(Failure::own(format!(
                    "the cache record {stable_id} holds a number where its attributes stand"
                )))
            }
        };
        let mut edits: Vec<(i128, i128)> = Vec::new();
        for attribute in attributes {
            let Value::Object(attribute) = attribute else {
                continue;
            };
            edits.extend(attribute_edit(attribute, &stable_id)?);
        }
        if edits.is_empty() {
            continue;
        }
        edits.sort_unstable();
        let mut offsets = BTreeSet::new();
        let mut running = 0i128;
        for (_start, length_change) in edits {
            running += length_change;
            offsets.insert(running.abs());
        }
        offsets.remove(&0);
        edited.insert(stable_id, offsets);
    }
    edited.remove("");
    Ok(edited)
}

/// Whether a differing `HGVSc` pair on a transcript with the edit distances `distances` is
/// two notations of one name that differ in their cDNA coordinates alone and by one of the
/// distances: some integer of the one body lies at such a distance from some integer of the
/// other. vep-rs numbers the position against the edited model, where Ensembl VEP walks the
/// exons in their unedited numbering and corrects an exonic `NM_` or `XM_` position alone
/// by the misalignment offset (TranscriptVariationAllele.pm:1517, 1524-1525, 1608-1610).
pub fn hgvsc_an_edit_apart(distances: &BTreeSet<i128>, vep: &str, vep_rs: &str) -> bool {
    let (Some((vep_name, vep_body)), Some((vep_rs_name, vep_rs_body))) =
        (vep.split_once(':'), vep_rs.split_once(':'))
    else {
        return false;
    };
    if vep_name != vep_rs_name || coordinates_masked(vep_body) != coordinates_masked(vep_rs_body) {
        return false;
    }
    let (Some(vep_integers), Some(vep_rs_integers)) =
        (body_integers(vep_body), body_integers(vep_rs_body))
    else {
        return false;
    };
    vep_integers.iter().any(|a| {
        vep_rs_integers
            .iter()
            .any(|b| distances.contains(&(a - b).abs()))
    })
}

/// Whether a differing `HGVS_OFFSET` pair is Ensembl VEP's empty field against a non-zero
/// integer of vep-rs: on an edited `NM_` or `XM_` model `_return_3prime` copies the first
/// shift hash attached to the variant, whatever transcript or strand made it
/// (TranscriptVariationAllele.pm:164-187), a zero that prints nothing, where vep-rs shifts
/// the row on its own transcript.
pub fn offset_left_empty(vep: &str, vep_rs: &str) -> bool {
    blank(vep) && non_zero_integer(vep_rs)
}

#[cfg(test)]
#[path = "edited_tests.rs"]
mod tests;
