// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of a protein notation Ensembl VEP reads from an inserted allele it left
//! unrotated: an insertion on a reverse-strand transcript that the 3' shift moved past the
//! allele's own length.

use super::notation::{blank, field, int, non_zero_integer, spells};
use crate::classify::downstream::hgvs_shifted_allele;
use crate::classify::fields::{FieldKey, FieldRow};

/// The columns and fields of a row on which two rows of one coding change agree.
pub const ONE_CODING_CHANGE_FIELDS: [&str; 4] =
    ["cDNA_position", "CDS_position", "HGVSc", "HGVS_OFFSET"];

/// The bases an inserted allele reads at the point a 3' shift of `offset` bases leaves it,
/// as the record's strand writes them: its first base moved to its end once per shifted
/// base on a forward-strand transcript, its last base moved to its start once per shifted
/// base on a reverse-strand one (`reverse`), the count taken modulo the allele's length.
pub fn shift_rotated_allele(allele: &str, offset: i128, reverse: bool) -> String {
    let mut bases: Vec<char> = allele.chars().collect();
    let held = bases.len();
    if held == 0 {
        return String::new();
    }
    let length = u128::try_from(held).unwrap_or(u128::MAX);
    let shifted = usize::try_from(offset.unsigned_abs() % length).unwrap_or(0);
    let turns = if reverse {
        (held - shifted) % held
    } else {
        shifted
    };
    bases.rotate_left(turns);
    bases.into_iter().collect()
}

/// Whether a differing `HGVSp` pair is Ensembl VEP's description read from an inserted
/// allele it left unrotated: on a reverse-strand transcript `shift_feature_seqs` turns the
/// allele by its length less the shift (TranscriptVariationAllele.pm:1337-1339), so a shift
/// longer than the allele turns nothing and the protein notation is read from bases the
/// record does not produce, where vep-rs reads the allele turned by the shift modulo its
/// length.
///
/// The row is a reverse-strand transcript's (`strand`), its Allele two or more of the bases
/// `ACGT` that the two rules leave different (`hgvs_shifted_allele`,
/// `shift_rotated_allele`), and both values name one protein before their last `:p.`. The
/// two rows state one coding change under one non-zero `HGVS_OFFSET`: the same cDNA and
/// coding positions and the same `HGVSc`, which the row carries.
///
/// The row's other fields are read from `perl_row` and `rust_row`, so a run that does not
/// compare `HGVSc` and `HGVS_OFFSET` names no row.
pub fn unrotated_shifted_insertion(
    name: &str,
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    strand: &str,
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
) -> bool {
    if name != "HGVSp" || key.feature_type != "Transcript" || strand != "-1" {
        return false;
    }
    let allele = key.allele.as_str();
    if allele.len() < 2 || !spells(allele, "ACGT") {
        return false;
    }
    let (Some((vep_name, vep_description)), Some((vep_rs_name, vep_rs_description))) =
        (vep.rsplit_once(":p."), vep_rs.rsplit_once(":p."))
    else {
        return false;
    };
    if vep_name.is_empty()
        || vep_name != vep_rs_name
        || vep_description.is_empty()
        || vep_rs_description.is_empty()
    {
        return false;
    }
    if ONE_CODING_CHANGE_FIELDS
        .iter()
        .any(|shared| field(perl_row, shared) != field(rust_row, shared))
    {
        return false;
    }
    let offset = field(perl_row, "HGVS_OFFSET");
    if blank(field(perl_row, "HGVSc")) || !non_zero_integer(offset) {
        return false;
    }
    let Some(shift) = int(offset) else {
        return false;
    };
    hgvs_shifted_allele(allele, shift, true) != shift_rotated_allele(allele, shift, true)
}

#[cfg(test)]
#[path = "unrotated_tests.rs"]
mod tests;
