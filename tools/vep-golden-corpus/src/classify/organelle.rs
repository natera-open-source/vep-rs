// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The organelle field class: the transcripts the corpus cache files under a codon table
//! that reads some codon unlike the standard one, and the four peptide paths of an `HGVSp`
//! in which Ensembl VEP translates a sequence with the standard code whatever the
//! transcript's table.

use super::cache::iter_cache_transcripts;
use super::downstream::{codon_table_digits, effect_feature_cache, filed_stable_id};
use super::fields::{FieldModels, FieldQuestion, FieldRun, Verdict};
use crate::error::Failure;
use crate::genetic_codes::reads_as_the_standard_code;
use crate::json::Value;
use std::collections::HashMap;

/// The class of an `HGVSp` on a transcript under a codon table other than the standard
/// one, both engines naming the same protein, where the row's consequence set and the two
/// descriptions place the difference in one of the four paths of
/// `organelle_peptide_path_divergence`: Ensembl VEP translates a sequence of each with the
/// standard code, the alternate at TranscriptVariationAllele.pm:2333, 2492 and 2555 and the
/// reference at 2450, against a reference peptide of the transcript's code.
pub const ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE: &str =
    "organelle_peptide_paths_translation_table";

/// `dup`, `ins` or `other` for the text after `p.`; a `delins` is `other`.
pub fn peptide_description_shape(description: &str) -> &'static str {
    if description.ends_with("dup") {
        return "dup";
    }
    let inserts = description
        .match_indices("ins")
        .any(|(at, _)| !description[..at].ends_with("del"));
    if inserts {
        "ins"
    } else {
        "other"
    }
}

/// Whether a differing HGVSp pair of a row with the given consequence set comes out of one
/// of the four paths of TranscriptVariationAllele.pm that translate a sequence bare (the
/// alternate sequence in `_get_fs_peptides`, `_stop_loss_extra_AA` and `_get_del_peptides`,
/// the reference sequence in `_check_for_peptide_duplication`): `_get_fs_peptides` and
/// `_stop_loss_extra_AA` write the whole description of a frameshift (its residues, a `Ter`
/// at an immediate stop, the stop distance and the lost-start shapes);
/// `_stop_loss_extra_AA` the count of an extension after a lost stop;
/// `_check_for_peptide_duplication` the choice between a duplication and an insertion.
/// `_get_del_peptides` writes a deletion whose reference window is empty: its rows re-type
/// to an insertion or a synonymous description and carry no `ext`, no corpus row reaches
/// it, and the rule does not name it, so such a row stays unexplained. A substitution, a
/// `delins` and a plain deletion read their residues through `peptide`, which already
/// follows the transcript's table.
pub fn organelle_peptide_path_divergence(consequences: &str, vep: &str, vep_rs: &str) -> bool {
    let carries = |term: &str| consequences.split(',').any(|held| held == term);
    if carries("frameshift_variant") {
        return true;
    }
    if carries("stop_lost") {
        return vep.contains("ext") && vep_rs.contains("ext");
    }
    let shapes = [
        peptide_description_shape(vep),
        peptide_description_shape(vep_rs),
    ];
    (carries("inframe_insertion") || carries("protein_altering_variant"))
        && (shapes == ["dup", "ins"] || shapes == ["ins", "dup"])
}

/// Each transcript of a corpus cache whose `codon_table` names a genetic code under which
/// Ensembl VEP reads some codon unlike the standard code (`reads_as_the_standard_code`,
/// found from the table of the codes, 0 read as table 1), by stable id, with that table: the
/// transcripts whose peptides Ensembl VEP translates with two codes. A converted cache
/// carries the table as text, a native one as a number; a table past 64 bits is held as the
/// largest number they hold. A table that names no genetic code a table is held for files
/// no transcript, since nothing tells how Ensembl VEP reads it: its rows stay unexplained.
pub fn organelle_codon_tables(cache_dir: &str) -> Result<HashMap<String, i64>, Failure> {
    let mut tables = HashMap::new();
    for item in iter_cache_transcripts(cache_dir) {
        let (_, transcript) = item?;
        let Some(cache) = effect_feature_cache(&transcript)? else {
            continue;
        };
        let table = match cache.get("codon_table") {
            Some(Value::Int(table)) if !table.is_negative() || table.is_zero() => {
                Some(table.as_i64().unwrap_or(i64::MAX))
            }
            Some(Value::Str(text)) => codon_table_digits(text, &transcript)?,
            _ => None,
        };
        let Some(table) = table.filter(|table| *table != 0) else {
            continue;
        };
        // A code that reads as the standard one, and a table no code is held for.
        if reads_as_the_standard_code(table) != Some(false) {
            continue;
        }
        if let Some(stable_id) = filed_stable_id(&transcript)? {
            tables.insert(stable_id.to_string(), table);
        }
    }
    Ok(tables)
}

/// The transcripts of the corpus cache under a codon table that reads some codon unlike the
/// standard one.
#[derive(Default)]
pub struct Models {
    /// `organelle_codon_tables` of the corpus cache.
    pub codon_tables: Option<HashMap<String, i64>>,
}

impl Models {
    pub fn load(run: &FieldRun<'_>) -> Result<Models, Failure> {
        Ok(Models {
            codon_tables: Some(organelle_codon_tables(run.cache_dir)?),
        })
    }
}

/// `organelle_peptide_paths_translation_table`: an `HGVSp` on a transcript of the models'
/// codon tables whose two values name one protein before `:p.` and whose descriptions
/// after it are a pair of `organelle_peptide_path_divergence`.
pub fn organelle_peptide_paths_translation_table(
    question: &FieldQuestion<'_>,
    models: &FieldModels,
) -> Verdict {
    let key = question.key;
    let filed = |tables: &HashMap<String, i64>| tables.contains_key(&key.feature);
    if question.name != "HGVSp" || !models.organelle.codon_tables.as_ref().is_some_and(filed) {
        return Verdict::Pass;
    }
    let (Some((vep_name, vep_description)), Some((vep_rs_name, vep_rs_description))) = (
        question.vep.rsplit_once(":p."),
        question.vep_rs.rsplit_once(":p."),
    ) else {
        return Verdict::Pass;
    };
    if !vep_name.is_empty()
        && vep_name == vep_rs_name
        && organelle_peptide_path_divergence(&key.consequences, vep_description, vep_rs_description)
    {
        return Verdict::Class(ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE, None);
    }
    Verdict::Pass
}

#[cfg(test)]
#[path = "organelle_tests.rs"]
mod tests;
