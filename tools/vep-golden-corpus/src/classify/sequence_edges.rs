// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The sequence-edge class of an HGVS shift: the questions `classify` puts to the format
//! comparator, whose rule the class is, and the answers it holds.

use super::comparator::{option, Comparator, Outcome};
use super::fields::{FieldModels, FieldQuestion, FieldRun, Verdict};
use super::spliceai::Accepted;
use crate::error::Failure;
use crate::pytext::{path_join, path_str};

/// The class of an `HGVSc`, an `HGVSp` or an `HGVS_OFFSET` and, under `--shift_3prime 1`, of
/// a `SHIFT_LENGTH` or a position column on the row of an insertion or a
/// deletion less than 1,000 bases from the edge of its sequence that the transcript's strand
/// runs towards: Ensembl VEP shifts the allele along the first or the last 1,000 bases of a
/// slice it clipped at that edge (TranscriptVariationAllele.pm:427-437), which are not the
/// allele's own flank, where vep-rs shifts it as far as the sequence admits. The format
/// comparator decides it from the reference sequence `classify --reference-fasta` names
/// and the transcripts of the corpus cache (`ask`).
pub const HGVS_SHIFT_SEQUENCE_EDGE_WINDOW: &str = "hgvs_shift_sequence_edge_window";

/// The option of `classify` that names the reference sequence.
pub const OPTION: &str = "--reference-fasta";

/// The fields of a row a question carries behind the row's Feature_type and strand, in this
/// order. Each of the first four is the row's value where the run compares the field, and
/// lacking where it does not; the last three are the row's position columns.
pub const SEQUENCE_EDGE_ROW_FIELDS: [&str; 7] = [
    "HGVSc",
    "HGVSp",
    "HGVS_OFFSET",
    "SHIFT_LENGTH",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
];

/// The cell of a question for a field its row lacks.
const LACKING: &str = "-";

/// What the class reads of a run: the fields the format comparator accepts beside the
/// run's reference sequence. The default is a run that names no reference sequence.
#[derive(Default)]
pub struct Models {
    pub accepted: Accepted,
}

/// One question of the class as a line of the command's questions: the key's Location,
/// Allele and Feature, the field, Ensembl VEP's value and vep-rs's, then the key's
/// Feature_type, the row's `STRAND` and the fields of `SEQUENCE_EDGE_ROW_FIELDS` as Ensembl
/// VEP's row writes them and as vep-rs's row does, `-` for a field a row lacks. The strand
/// is Ensembl VEP's row's, and vep-rs's row states none.
pub fn question_line(question: &FieldQuestion<'_>) -> String {
    let key = question.key;
    let mut cells = vec![
        key.location.as_str(),
        key.allele.as_str(),
        key.feature.as_str(),
        question.name,
        question.vep,
        question.vep_rs,
    ];
    for (row, strand) in [
        (question.perl_row, question.strand),
        (question.rust_row, LACKING),
    ] {
        cells.extend([key.feature_type.as_str(), strand]);
        cells.extend(SEQUENCE_EDGE_ROW_FIELDS.iter().map(|name| {
            row.and_then(|row| row.get(*name))
                .map_or(LACKING, String::as_str)
        }));
    }
    cells.join("\t") + "\n"
}

/// The fields of `questions` the class accepts: the answers of one run of the format
/// comparator's `field-class` command over every question, beside the reference sequence
/// `reference` and the corpus cache. A run without a question is asked as well, since the
/// command reads the sequence's index and the cache's transcripts whatever it is asked. A
/// corpus cache without a `transcripts` directory is refused before the command starts. A
/// run of the command that leaves no answers is this run's failure: the detail of the
/// command's own error line where it ends with exit code 2 on one, as for a reference
/// sequence or a cache shard it cannot read, and else the line of `Comparator::ended`. So
/// are answers that are not one `0` or `1` a question.
pub fn ask(
    comparator: &Comparator,
    run: &FieldRun<'_>,
    reference: &str,
    questions: &[FieldQuestion<'_>],
) -> Result<Accepted, Failure> {
    let transcripts = path_join(run.cache_dir, "transcripts");
    if !std::fs::metadata(transcripts).is_ok_and(|held| held.is_dir()) {
        return Err(Failure::own(format!(
            "{OPTION}: the corpus cache {} holds no transcripts directory, which \
             {HGVS_SHIFT_SEQUENCE_EDGE_WINDOW} reads",
            run.cache_dir
        )));
    }
    let options = [
        option(OPTION, path_str(reference)),
        option("--json-cache", run.cache_dir),
    ];
    let lines: String = questions.iter().map(question_line).collect();
    let answered =
        match comparator.field_class(HGVS_SHIFT_SEQUENCE_EDGE_WINDOW, &options, &lines)? {
            Outcome::Answered(answered) => answered,
            Outcome::Ended(ended) => {
                let refused = ended.detail().filter(|_| ended.code == Some(2));
                return Err(refused.map_or_else(|| comparator.ended(&ended), Failure::own));
            }
        };
    let accepts =
        comparator.accepted(HGVS_SHIFT_SEQUENCE_EDGE_WINDOW, questions.len(), &answered)?;
    let mut accepted = Accepted::new();
    for (question, accepts) in questions.iter().zip(accepts) {
        if accepts {
            let fields = accepted.entry(question.key.clone()).or_default();
            fields.insert(question.name.to_string());
        }
    }
    Ok(accepted)
}

/// The class `hgvs_shift_sequence_edge_window` of the dispatch: the field is one the format
/// comparator accepts on the question's key (`Models::accepted`). Nothing in a run without
/// a reference sequence, whose models hold no answer.
pub fn hgvs_shift_sequence_edge_window(
    question: &FieldQuestion<'_>,
    models: &FieldModels,
) -> Verdict {
    let accepted = models.sequence_edges.accepted.get(question.key);
    if accepted.is_some_and(|fields| fields.contains(question.name)) {
        return Verdict::Class(HGVS_SHIFT_SEQUENCE_EDGE_WINDOW, None);
    }
    Verdict::Pass
}

#[cfg(test)]
#[path = "sequence_edges_tests.rs"]
mod tests;
