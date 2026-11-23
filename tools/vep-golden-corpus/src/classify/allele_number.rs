// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The allele-number class of a record whose written allele list repeats an alternate or
//! lists an alternate equal to the reference allele: Ensembl VEP numbers a row's allele
//! over the list without those entries, vep-rs by the allele's place in the list as
//! written.

use super::fields::{FieldModels, FieldQuestion, FieldRow, Verdict};
use crate::allele_number::{
    allele_number_place_pair, ALLELE_NUMBER_FIELDS, WRITTEN_ALLELE_LIST_FIELDS,
};

/// The class of an allele number on a row of a record that writes an alternate twice, or an
/// alternate equal to the reference allele, before the row's allele:
/// VariationFeatureOverlap.pm:459-467 at release/116.3 drop such an entry from the record's
/// allele list and :477-491 number the entries left, where vep-rs numbers the allele by its
/// place among the written alternates.
pub const ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT: &str =
    "allele_number_repeated_or_reference_equal_alt";

/// The field under which a default-format row carries its record's written allele list. The
/// class reads it of both rows, whichever fields the run compares.
pub const WRITTEN_ALLELE_LIST_FIELD: &str = WRITTEN_ALLELE_LIST_FIELDS[0];

/// The record's written allele list as `row` carries it; `None` for a row without the
/// field, which is read as `-`, and for a row read without it.
fn written_list(row: &FieldRow) -> Option<&str> {
    row.get(WRITTEN_ALLELE_LIST_FIELD)
        .map(String::as_str)
        .filter(|list| *list != "-")
}

/// `allele_number_repeated_or_reference_equal_alt`: an `ALLELE_NUM`, or the `allele_num`
/// the JSON format writes, whose two values have the shape the format comparator's rule
/// states (`crate::allele_number::allele_number_place_pair`, that tool's own file). Both
/// default-format rows of the key carry the record's written allele list with one text; the
/// key's allele is a written alternate or, on an indel record, one without the first
/// character Ensembl VEP's VCF parser chops (Parser/VCF.pm:294-335); vep-rs's value is the
/// place of the first such alternate among the written alternates, and Ensembl VEP's the
/// count of the alternates up to that place that equal neither the reference allele nor an
/// alternate written before them. An allele that is a written alternate at one place and
/// the parser's spelling of an alternate at another passes: a row does not say which of the
/// two spellings its allele is. So does an allele whose alternate equals the reference
/// allele, an entry Ensembl VEP drops and writes no row for. A question without both rows
/// passes.
pub fn allele_number_repeated_or_reference_equal_alt(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let (Some(perl_row), Some(rust_row)) = (question.perl_row, question.rust_row) else {
        return Verdict::Pass;
    };
    if ALLELE_NUMBER_FIELDS.contains(&question.name)
        && allele_number_place_pair(
            written_list(perl_row),
            written_list(rust_row),
            &question.key.allele,
            question.vep,
            question.vep_rs,
        )
    {
        return Verdict::Class(ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT, None);
    }
    Verdict::Pass
}

#[cfg(test)]
#[path = "allele_number_tests.rs"]
mod tests;
