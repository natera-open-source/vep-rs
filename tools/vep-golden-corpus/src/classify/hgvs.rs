// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HGVS field classes: an edited transcript model, the written-span gate of a shifted
//! insertion or deletion, a clamped deletion, the protein notation of an allele shifted
//! into the first codon, the protein notation of an allele shifted onto a start codon it
//! leaves standing, a reference clipped against the whole alternate allele, added bases
//! written outside the transcript, and the protein notation of an inserted allele Ensembl
//! VEP leaves unrotated.
//!
//! Each class reads the differing field with the row's other compared fields, and the two
//! classes of an edited model read the corpus cache. `notation` holds the HGVSc patterns
//! and the cDNA coordinates every class reads, `edited`, `shifted`, `clipped`, `outside`
//! and `unrotated` the classes' predicates.

pub mod clipped;
pub mod edited;
pub mod notation;
pub mod outside;
pub mod shifted;
pub mod unrotated;

use super::fields::{FieldModels, FieldQuestion, FieldRun, Verdict};
use crate::error::Failure;
pub use edited::{indel_edited_models, EditedModels};
pub use notation::{
    cdna_coordinate_index, cdna_coordinate_order, exonic_run, HGVSC_CLIPPED_BODY,
    HGVS_INTRON_OFFSET_FIELDS,
};
use std::collections::BTreeSet;

/// The class of an `HGVSc` on a transcript the corpus cache files with an `_rna_edit` that
/// inserts or deletes bases, the two notations a whole edit apart in their coordinates.
pub const REFSEQ_EDITED_MODEL_HGVSC_OFFSET: &str = "refseq_edited_model_hgvsc_offset";

/// The class of an `HGVS_OFFSET` on such a transcript that Ensembl VEP leaves empty and
/// vep-rs fills with a non-zero integer.
pub const REFSEQ_EDITED_MODEL_HGVS_OFFSET: &str = "refseq_edited_model_hgvs_offset";

/// The class of an `HGVSp` Ensembl VEP reads from an inserted allele that its rotation by
/// the transcript's strand left as written.
pub const SHIFTED_ALLELE_ROTATION_BY_STRAND: &str = "shifted_allele_rotation_by_strand";

/// The class of a notation Ensembl VEP leaves empty on a shifted insertion or deletion that
/// vep-rs describes at its shifted span.
pub const HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE: &str = "hgvs_shifted_indel_written_span_gate";

/// The class of two deletions that end together, Ensembl VEP's starting later.
pub const HGVS_SHIFTED_DELETION_CLAMPED_SPAN: &str = "hgvs_shifted_deletion_clamped_span";

/// The class of an `HGVSp` Ensembl VEP makes from the residues of an insertion or deletion
/// moved 3' into the first codon and the start codon test of the allele as written.
pub const HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST: &str =
    "hgvsp_shifted_indel_written_span_start_test";

/// The class of an `HGVSp` Ensembl VEP reads from the coding sequence edited without the
/// 5' UTR, on an insertion or deletion moved 3' onto a start codon it leaves standing.
pub const HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE: &str =
    "hgvsp_start_codon_test_past_the_shifted_allele";

/// The class of a notation Ensembl VEP reads from the whole alternate allele over the bases
/// its clamp left inside the transcript, and of the HGVSIntronOffset plugin's fields beside
/// such a notation.
pub const HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE: &str =
    "hgvs_clipped_reference_against_whole_alternate";

/// The class of vep-rs's notation of inserted bases on a row both engines annotate upstream
/// or downstream of the transcript.
pub const HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT: &str =
    "hgvs_added_bases_written_outside_the_transcript";

/// The transcripts of the corpus cache whose model is edited against the genome, each with
/// its edit distances (`indel_edited_models`). The default holds no transcript.
#[derive(Default)]
pub struct Models {
    pub edited_models: EditedModels,
}

impl Models {
    pub fn load(run: &FieldRun<'_>) -> Result<Models, Failure> {
        Ok(Models {
            edited_models: indel_edited_models(run.cache_dir)?,
        })
    }

    /// The edit distances of the transcript `feature`, `None` for a transcript without one.
    fn edit_distances(&self, feature: &str) -> Option<&BTreeSet<i128>> {
        self.edited_models
            .get(feature)
            .filter(|distances| !distances.is_empty())
    }
}

/// The class `name` where the predicate holds, a pass where it does not.
fn class_where(holds: bool, name: &'static str) -> Verdict {
    if holds {
        Verdict::Class(name, None)
    } else {
        Verdict::Pass
    }
}

/// `refseq_edited_model_hgvsc_offset`: an `HGVSc` on a transcript of the edited models
/// whose two notations differ in their cDNA coordinates alone and by one of the model's
/// edit distances.
pub fn refseq_edited_model_hgvsc_offset(
    question: &FieldQuestion<'_>,
    models: &FieldModels,
) -> Verdict {
    if question.name != "HGVSc" {
        return Verdict::Pass;
    }
    let apart = models
        .hgvs
        .edit_distances(&question.key.feature)
        .is_some_and(|distances| {
            edited::hgvsc_an_edit_apart(distances, question.vep, question.vep_rs)
        });
    class_where(apart, REFSEQ_EDITED_MODEL_HGVSC_OFFSET)
}

/// `hgvs_added_bases_written_outside_the_transcript`: an `HGVSc`, `HGVS_OFFSET` or
/// HGVSIntronOffset plugin field on the row of inserted bases both engines annotate
/// upstream or downstream of the transcript.
pub fn added_bases_written_outside_the_transcript(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let outside = outside::added_bases_written_outside_the_transcript(
        question.name,
        question.key,
        question.vep,
        question.vep_rs,
        question.strand,
        question.perl_row,
        question.rust_row,
    );
    class_where(outside, HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT)
}

/// `refseq_edited_model_hgvs_offset`: an `HGVS_OFFSET` on a transcript of the edited models
/// that Ensembl VEP leaves empty and vep-rs fills with a non-zero integer.
pub fn refseq_edited_model_hgvs_offset(
    question: &FieldQuestion<'_>,
    models: &FieldModels,
) -> Verdict {
    let left_empty = question.name == "HGVS_OFFSET"
        && models.hgvs.edit_distances(&question.key.feature).is_some()
        && edited::offset_left_empty(question.vep, question.vep_rs);
    class_where(left_empty, REFSEQ_EDITED_MODEL_HGVS_OFFSET)
}

/// `shifted_allele_rotation_by_strand`: an `HGVSp` of an inserted allele of two or more
/// bases on a reverse-strand transcript that the 3' shift moved past the allele's length
/// by no multiple of it, the two rows stating one coding change. It is asked before
/// `organelle_peptide_paths_translation_table`, which takes the `HGVSp` of every
/// frameshift on a transcript under another codon table.
pub fn unrotated_shifted_insertion(question: &FieldQuestion<'_>, _models: &FieldModels) -> Verdict {
    let unrotated = unrotated::unrotated_shifted_insertion(
        question.name,
        question.key,
        question.vep,
        question.vep_rs,
        question.strand,
        question.perl_row,
        question.rust_row,
    );
    class_where(unrotated, SHIFTED_ALLELE_ROTATION_BY_STRAND)
}

/// `hgvs_shifted_indel_written_span_gate`: an `HGVSc`, `HGVS_OFFSET` or HGVSIntronOffset
/// plugin field Ensembl VEP leaves empty on a shifted insertion or deletion vep-rs
/// describes at its shifted span.
pub fn shifted_indel_written_span_gate(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let gated = shifted::shifted_indel_written_span_gate(
        question.name,
        question.key,
        question.vep,
        question.vep_rs,
        question.perl_row,
        question.rust_row,
    );
    class_where(gated, HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE)
}

/// `hgvs_shifted_deletion_clamped_span`: an `HGVSc` whose two deletions end together and
/// Ensembl VEP's starts later, on a row Ensembl VEP gives a non-zero `HGVS_OFFSET`.
pub fn shifted_deletion_clamped_span(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let clamped = question.name == "HGVSc"
        && shifted::shifted_deletion_clamped_span(
            question.key,
            question.vep,
            question.vep_rs,
            question.perl_row,
        );
    class_where(clamped, HGVS_SHIFTED_DELETION_CLAMPED_SPAN)
}

/// `hgvsp_shifted_indel_written_span_start_test`: an `HGVSp` of an insertion or deletion
/// the 3' shift moves into the first codon, Ensembl VEP's a frameshift at residue 1 with
/// the unknown count or a deletion of residues from the first, vep-rs's the notation of
/// the allele written at its 3' place.
pub fn shifted_indel_written_span_start_test(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let tested = question.name == "HGVSp"
        && shifted::shifted_indel_written_span_start_test(
            question.key,
            question.vep,
            question.vep_rs,
            question.perl_row,
            question.rust_row,
        );
    class_where(tested, HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST)
}

/// `hgvsp_start_codon_test_past_the_shifted_allele`: an `HGVSp` of an insertion or deletion
/// the 3' shift moves onto a start codon it leaves standing, Ensembl VEP's a frameshift at
/// residue 1 with the unknown count or an insertion between a residue 0 and residue 1,
/// vep-rs's the first residue unchanged.
pub fn start_codon_test_past_the_shifted_allele(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let tested = question.name == "HGVSp"
        && shifted::start_codon_test_past_the_shifted_allele(
            question.key,
            question.vep,
            question.vep_rs,
            question.perl_row,
            question.rust_row,
        );
    class_where(tested, HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE)
}

/// `hgvs_clipped_reference_against_whole_alternate` for an `HGVSc`: a sequence allele's row
/// whose cDNA position has an open side, Ensembl VEP's notation read from the whole
/// alternate allele over the bases inside the transcript and vep-rs's the change inside.
pub fn clipped_reference_against_whole_alternate(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let clipped = question.name == "HGVSc"
        && clipped::clipped_reference_against_whole_alternate(
            question.key,
            question.vep,
            question.vep_rs,
            question.perl_row,
        );
    class_where(clipped, HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE)
}

/// `hgvs_clipped_reference_against_whole_alternate` for a field of the HGVSIntronOffset
/// plugin beside a notation the class takes.
pub fn clipped_notation_intron_offsets(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let beside = clipped::clipped_notation_intron_offsets(
        question.name,
        question.key,
        question.vep,
        question.vep_rs,
        question.perl_row,
        question.rust_row,
    );
    class_where(beside, HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE)
}

#[cfg(test)]
#[path = "hgvs_tests.rs"]
mod tests;
