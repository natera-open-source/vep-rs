// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of a key on an edited RefSeq model whose reference span holds a whole intron
//! between two coding bases: Ensembl VEP reads the genome's reference allele, intron
//! included, as the coding reference, so the coding terms it tests on that reference
//! differ from the ones vep-rs tests on the model's own codons.

use super::fields::FieldRow;
use std::collections::BTreeSet;

/// The class of `docs/divergences/refseq-edited-model-intron-span-coding-reference.md`.
pub const REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE: &str =
    "refseq_edited_model_intron_span_coding_reference";

/// The Extra fields the class reads beside a row's `CDS_position` column.
pub const INTRON_SPAN_FIELDS: [&str; 2] = ["BAM_EDIT", "USED_REF"];

/// The terms both engines write on a row whose span holds a whole intron: the span covers
/// the intron's two splice sites and bases of the intron itself.
const WHOLE_INTRON_TERMS: [&str; 3] = [
    "splice_acceptor_variant",
    "splice_donor_variant",
    "intron_variant",
];

/// The coding terms Ensembl VEP alone writes: its reference peptide, translated through the
/// intron's bases, holds a stop the alternate peptide lacks, and its reference codons neither
/// begin nor end with the alternate codons.
const REFERENCE_ONLY_TERMS: [&str; 2] = ["stop_lost", "protein_altering_variant"];

/// The coding term vep-rs alone writes, from the model's own codons.
const CANDIDATE_ONLY_TERMS: [&str; 1] = ["inframe_deletion"];

/// The number of coding bases a `CDS_position` cell spans: `None` for a cell that is no
/// range from one number to a greater one (one position, an end written `?`, no position).
fn coding_bases(cds_position: &str) -> Option<usize> {
    let (start, end) = cds_position.split_once('-')?;
    let number = |text: &str| {
        let digits = !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit());
        digits.then(|| text.parse::<usize>().ok()).flatten()
    };
    let (start, end) = (number(start)?, number(end)?);
    (start < end).then(|| end - start + 1)
}

/// Whether the allele is bases, or `-` for a deletion's.
fn sequence_allele(allele: &str) -> bool {
    !allele.is_empty() && allele.bytes().all(|letter| b"ACGTN-".contains(&letter))
}

/// Whether a row is one of an edited model whose reference allele holds bases its coding
/// span lacks: `BAM_EDIT` reads `OK`, `CDS_position` is a range with a number at both ends,
/// and `USED_REF` spells more bases than the range holds.
fn reference_outruns_the_coding_span(row: &FieldRow) -> bool {
    let field = |name: &str| row.get(name).map_or("-", String::as_str);
    let used_ref = field("USED_REF");
    field("BAM_EDIT") == "OK"
        && used_ref.bytes().all(|byte| byte.is_ascii_alphabetic())
        && coding_bases(field("CDS_position")).is_some_and(|bases| used_ref.len() > bases)
}

/// Whether `set` holds the terms of `shape` and no other.
fn is_set(set: &BTreeSet<&str>, shape: &[&str]) -> bool {
    set.len() == shape.len() && shape.iter().all(|term| set.contains(term))
}

/// `refseq_edited_model_intron_span_coding_reference` for a key whose Ensembl VEP set
/// `perl_csq` and vep-rs set `rust_csq` differ as the class documents, read from the two
/// rows alone: `reference_row` and `candidate_row` hold each side's `CDS_position` and the
/// fields of [`INTRON_SPAN_FIELDS`].
///
/// The key is a Transcript's and its `allele` is bases, or `-` for a deletion's. Both sets
/// hold the three terms of a whole intron inside the span. Ensembl VEP's set less vep-rs's
/// is `stop_lost` and `protein_altering_variant`, and vep-rs's less Ensembl VEP's is
/// `inframe_deletion`. Both rows carry `BAM_EDIT` `OK`, the same `CDS_position` range with
/// a number at both ends and the same `USED_REF`, which spells more bases than the range
/// holds. A key that misses one of these takes no class.
pub fn classify_intron_span_pair(
    feature_type: &str,
    allele: &str,
    perl_csq: &str,
    rust_csq: &str,
    reference_row: Option<&FieldRow>,
    candidate_row: Option<&FieldRow>,
) -> Option<&'static str> {
    if feature_type != "Transcript" || !sequence_allele(allele) {
        return None;
    }
    let (reference_row, candidate_row) = (reference_row?, candidate_row?);
    let perl: BTreeSet<&str> = perl_csq.split(',').collect();
    let rust: BTreeSet<&str> = rust_csq.split(',').collect();
    let whole_intron = WHOLE_INTRON_TERMS
        .iter()
        .all(|term| perl.contains(term) && rust.contains(term));
    let perl_only: BTreeSet<&str> = perl.difference(&rust).copied().collect();
    let rust_only: BTreeSet<&str> = rust.difference(&perl).copied().collect();
    let same = |name: &str| reference_row.get(name) == candidate_row.get(name);
    (whole_intron
        && is_set(&perl_only, &REFERENCE_ONLY_TERMS)
        && is_set(&rust_only, &CANDIDATE_ONLY_TERMS)
        && reference_outruns_the_coding_span(reference_row)
        && reference_outruns_the_coding_span(candidate_row)
        && same("CDS_position")
        && same("USED_REF"))
    .then_some(REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE)
}

#[cfg(test)]
#[path = "intron_span_tests.rs"]
mod tests;
