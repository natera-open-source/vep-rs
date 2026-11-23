// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The field classes: the question one class is asked on one field of one key, the
//! models the classes read, and the dispatch that asks the classes in their order.

use crate::json::Object;
use std::collections::HashMap;

/// The key under which a `field_divergences` entry names the class of a second cause that
/// changes the same value.
pub const SECOND_EXPECTED_DIVERGENCE: &str = "second_expected_divergence";

/// The key under which a `field_divergences` entry counts the rows of its key and
/// consequence set on which both engines print one value of the field, the rows of other
/// records that share the key with the row the entry states. An entry of a key all of whose
/// rows differ carries no such key.
pub const ROWS_ALIKE: &str = "rows_alike";

/// The class of a plugin field on a row of a multi-allelic record `--minimal` split, where
/// Ensembl VEP's output pass hands the plugin the rejoined record's span and alleles
/// (OutputFactory.pm:2300-2302, :2323, :2235) and vep-rs the allele's own.
pub const MINIMAL_PLUGIN_SEES_THE_REJOINED_RECORD: &str = "minimal_plugin_sees_the_rejoined_record";

/// The fields of the plugins that class covers whose names are fixed, and the prefixes of
/// the others.
pub const REJOINED_RECORD_PLUGIN_FIELDS: [&str; 14] = [
    "Conservation",
    "NearestExonJB",
    "TSSDistance",
    "ClinPred",
    "LOEUF",
    "PrimateAI",
    "REVEL",
    "ada_score",
    "rf_score",
    "pHaplo",
    "pTriplo",
    "Existing_InFrame_oORFs",
    "Existing_OutOfFrame_oORFs",
    "Existing_uORFs",
];
pub const REJOINED_RECORD_PLUGIN_PREFIXES: [&str; 11] = [
    "gnomAD_",
    "GWAS_",
    "AVI_",
    "OpenTargets_",
    "am_",
    "BayesDel_",
    "CADD_",
    "MaxEntScan_",
    "SpliceAI_pred",
    "SpliceVault_",
    "5UTR_",
];

/// The class of a `MOTIF_SCORE_CHANGE` on a reverse-strand motif row whose allele is two or
/// more bases over as many positions: Ensembl VEP writes the allele from the motif position
/// of the variant's genomic start, the allele's last position along a reverse-strand motif
/// (MotifFeatureVariationAllele.pm:245).
pub const MOTIF_SCORE_REVERSE_STRAND_MNV: &str = "motif_score_reverse_strand_mnv";

/// The class of a `CELL_TYPE` Ensembl VEP leaves empty under `--cell_type` where vep-rs
/// prints an activity: its comma split never matches a name holding a comma (Config.pm:760-762
/// and OutputFactory.pm:1968-1978).
pub const CELL_TYPE_COMMA_NAMES: &str = "cell_type_comma_names";

/// The class of the `Codons` Ensembl VEP cuts from a coding sequence with the 3' UTR joined
/// on (TranscriptVariationAllele.pm:2422-2433 and :877) and of the `Amino_acids` read from
/// them, where vep-rs cuts both codons from coding bases.
pub const ALTERNATE_CODON_PAST_THE_CODING_END: &str = "alternate_codon_past_the_coding_end";

/// The key of a field comparison: Location, Allele, Feature, Feature_type and consequence
/// set, ordered in that order.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct FieldKey {
    pub location: String,
    pub allele: String,
    pub feature: String,
    pub feature_type: String,
    pub consequences: String,
}

/// The compared fields of one row by name, `-` for a field the row lacks, with the row's
/// `cDNA_position`, `CDS_position`, `Protein_position` and `Amino_acids` columns.
pub type FieldRow = HashMap<String, String>;

/// What a field class is asked: one field of one key whose consequence sets agree and
/// whose two values do not.
#[derive(Clone, Copy, Debug)]
pub struct FieldQuestion<'a> {
    /// The field.
    pub name: &'a str,
    pub key: &'a FieldKey,
    /// Ensembl VEP's value.
    pub vep: &'a str,
    /// vep-rs's value.
    pub vep_rs: &'a str,
    /// The `STRAND` of Ensembl VEP's row.
    pub strand: &'a str,
    /// The run's flags.
    pub flags: &'a [String],
    pub perl_row: Option<&'a FieldRow>,
    pub rust_row: Option<&'a FieldRow>,
    /// Whether the row is one of a multi-allelic record `--minimal` split.
    pub minimal_split_record: bool,
    /// The start and end of the record the row's Location prints, for a field of the
    /// Downstream plugin.
    pub genomic_span: Option<(i64, i64)>,
}

/// What a class answers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Verdict {
    /// The class does not take the field: the next class is asked.
    Pass,
    /// The class takes it: its name, and the second cause's class where two causes change
    /// the value.
    Class(&'static str, Option<&'static str>),
    /// No class takes the field and no further class is asked.
    NoClass,
}

/// One class of the dispatch.
pub type FieldClass = fn(&FieldQuestion<'_>, &FieldModels) -> Verdict;

/// What the classes read of the corpus cache and of the run's records, and what the format
/// comparator answers under the score files and beside the reference sequence, one member a
/// class module. The default holds no model.
#[derive(Default)]
pub struct FieldModels {
    pub hgvs: super::hgvs::Models,
    pub organelle: super::organelle::Models,
    pub downstream: super::downstream::Models,
    pub protein: super::protein::Models,
    pub introns: super::introns::Models,
    pub spliceai: super::spliceai::Models,
    pub sequence_edges: super::sequence_edges::Models,
}

/// What a class module's `Models::load` is given.
#[derive(Clone, Copy)]
pub struct FieldRun<'a> {
    pub corpus: &'a str,
    /// The corpus's cache directory (`cache::corpus_cache_dir`).
    pub cache_dir: &'a str,
    pub manifest: &'a Object,
    /// The compared fields, `--field` in order.
    pub fields: &'a [String],
    /// The score files of `--spliceai-table`.
    pub spliceai_tables: &'a [String],
    /// The reference sequence of `--reference-fasta`.
    pub reference_fasta: Option<&'a str>,
}

/// The classes in the order they are asked; the first that does not pass decides.
pub const FIELD_CLASSES: [FieldClass; 27] = [
    minimal_plugin_sees_the_rejoined_record,
    motif_score_reverse_strand_mnv,
    cell_type_comma_names,
    super::hgvs::refseq_edited_model_hgvsc_offset,
    super::hgvs::added_bases_written_outside_the_transcript,
    super::hgvs::refseq_edited_model_hgvs_offset,
    super::hgvs::unrotated_shifted_insertion,
    super::organelle::organelle_peptide_paths_translation_table,
    super::hgvs::shifted_indel_written_span_gate,
    super::hgvs::shifted_deletion_clamped_span,
    super::hgvs::shifted_indel_written_span_start_test,
    super::hgvs::start_codon_test_past_the_shifted_allele,
    super::hgvs::clipped_reference_against_whole_alternate,
    super::hgvs::clipped_notation_intron_offsets,
    super::protein::stop_codon_protein_comparison,
    super::protein::stop_named_on_incomplete_coding_sequence,
    super::downstream::downstream_field_class,
    super::introns::intron_end_insertion,
    super::spliceai::spliceai_multi_entry_cell,
    super::introns::intron_start_insertion,
    super::introns::maxentscan_overlapping_windows,
    alternate_codon_past_the_coding_end,
    super::protein::frameshift_count_at_the_first_residue,
    super::allele_number::allele_number_repeated_or_reference_equal_alt,
    super::loftee::gerp_walk_without_conservation_database,
    super::protein::frameshift_with_no_changed_residue,
    super::sequence_edges::hgvs_shift_sequence_edge_window,
];

/// The place in `FIELD_CLASSES` of `spliceai_multi_entry_cell`, whose answers come from one
/// run of the format comparator over every field that reaches the class: a run with score
/// files asks the classes before this place of every differing field first.
pub const ASKED_OF_THE_COMPARATOR: usize = 18;

/// The place in `FIELD_CLASSES` of `hgvs_shift_sequence_edge_window`, the last, whose
/// answers come from one run of the format comparator as well: a run with a reference
/// sequence asks every class before it of every differing field first, so the class takes
/// a field no other class takes.
pub const SEQUENCE_EDGES_ASKED_OF_THE_COMPARATOR: usize = 26;

/// The documented class of a field divergence on a key whose consequence sets agree, with
/// the class of a second cause that changes the value, or `None` where no class takes it.
pub fn classify_field_divergence(
    question: &FieldQuestion<'_>,
    models: &FieldModels,
) -> Option<(&'static str, Option<&'static str>)> {
    for class in FIELD_CLASSES {
        match class(question, models) {
            Verdict::Pass => {}
            Verdict::Class(name, second) => return Some((name, second)),
            Verdict::NoClass => return None,
        }
    }
    None
}

/// Whether the field reaches the class at `place` of `FIELD_CLASSES`: every class asked
/// before it passes.
fn reaches(place: usize, question: &FieldQuestion<'_>, models: &FieldModels) -> bool {
    FIELD_CLASSES[..place]
        .iter()
        .all(|class| class(question, models) == Verdict::Pass)
}

/// Whether the field reaches the class at `ASKED_OF_THE_COMPARATOR`: every class asked
/// before it passes, so the field is a question for the format comparator.
pub fn reaches_the_comparator(question: &FieldQuestion<'_>, models: &FieldModels) -> bool {
    reaches(ASKED_OF_THE_COMPARATOR, question, models)
}

/// Whether the field reaches the class at `SEQUENCE_EDGES_ASKED_OF_THE_COMPARATOR`: every
/// other class passes, under the answers `models` holds for the class at
/// `ASKED_OF_THE_COMPARATOR`, so the field is a question for the format comparator.
pub fn reaches_the_sequence_edge_class(question: &FieldQuestion<'_>, models: &FieldModels) -> bool {
    reaches(SEQUENCE_EDGES_ASKED_OF_THE_COMPARATOR, question, models)
}

/// Whether `name` is a field of a plugin the class `minimal_plugin_sees_the_rejoined_record`
/// covers.
pub fn rejoined_record_plugin_field(name: &str) -> bool {
    REJOINED_RECORD_PLUGIN_FIELDS.contains(&name)
        || REJOINED_RECORD_PLUGIN_PREFIXES
            .iter()
            .any(|prefix| name.starts_with(prefix))
}

/// `minimal_plugin_sees_the_rejoined_record`: a covered plugin's field on a row of a
/// multi-allelic record `--minimal` split.
pub fn minimal_plugin_sees_the_rejoined_record(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    if question.minimal_split_record && rejoined_record_plugin_field(question.name) {
        return Verdict::Class(MINIMAL_PLUGIN_SEES_THE_REJOINED_RECORD, None);
    }
    Verdict::Pass
}

/// The count of positions a Location's span names, 1 for a span that is not two numbers in
/// ASCII digits; `None` for numbers, or a count, 128 bits do not hold.
fn span_positions(location: &str) -> Option<i128> {
    let span = location.split_once(':').map_or(location, |(_, span)| span);
    let (start, end) = span.split_once('-').unwrap_or((span, ""));
    let digits = |text: &str| !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit());
    if !digits(start) || !digits(end) {
        return Some(1);
    }
    let (start, end) = (start.parse::<i128>().ok()?, end.parse::<i128>().ok()?);
    (end - start).checked_add(1)
}

/// `motif_score_reverse_strand_mnv`: `MOTIF_SCORE_CHANGE` on a MotifFeature row of `STRAND`
/// -1 whose allele is two or more of the bases `ACGT` over as many positions.
pub fn motif_score_reverse_strand_mnv(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let key = question.key;
    if question.name != "MOTIF_SCORE_CHANGE"
        || key.feature_type != "MotifFeature"
        || question.strand != "-1"
    {
        return Verdict::Pass;
    }
    let bases = key.allele.len();
    let sequence = key.allele.bytes().all(|base| b"ACGT".contains(&base));
    if bases >= 2 && sequence && span_positions(&key.location) == i128::try_from(bases).ok() {
        return Verdict::Class(MOTIF_SCORE_REVERSE_STRAND_MNV, None);
    }
    Verdict::Pass
}

/// `cell_type_comma_names`: a `CELL_TYPE` Ensembl VEP leaves empty where vep-rs prints an
/// activity, in a run whose flags carry `--cell_type`.
pub fn cell_type_comma_names(question: &FieldQuestion<'_>, _models: &FieldModels) -> Verdict {
    let asked = |flag: &String| flag == "--cell_type" || flag.starts_with("--cell_type=");
    if question.name == "CELL_TYPE"
        && question.vep.is_empty()
        && !matches!(question.vep_rs, "" | "-")
        && question.flags.iter().any(asked)
    {
        return Verdict::Class(CELL_TYPE_COMMA_NAMES, None);
    }
    Verdict::Pass
}

/// `alternate_codon_past_the_coding_end`: a `Codons` whose two values have one of the
/// class's two shapes, or an `Amino_acids` on a row whose `Codons` have the first
/// (`crate::coding_end`). A codon window opens at residue 1 where the row's first coding
/// position is one of the first three, and the residues are read beside both rows'
/// `Codons`, so a run that compares `Amino_acids` compares `Codons` too.
pub fn alternate_codon_past_the_coding_end(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    let held = |row: Option<&FieldRow>, name: &str| -> Option<String> {
        row.and_then(|row| row.get(name)).cloned()
    };
    let taken = match question.name {
        "Codons" => {
            let position = held(question.perl_row, "CDS_position").unwrap_or_default();
            let first = position.split(['-', '/']).next().unwrap_or_default();
            let first_residue = if matches!(first, "1" | "2" | "3") {
                "1"
            } else {
                ""
            };
            crate::coding_end::codons_read_past_the_coding_end(
                first_residue,
                question.vep,
                question.vep_rs,
            )
        }
        "Amino_acids" => match (
            held(question.perl_row, "Codons"),
            held(question.rust_row, "Codons"),
        ) {
            (Some(vep), Some(vep_rs)) => crate::coding_end::residues_read_past_the_coding_end(
                &vep,
                &vep_rs,
                question.vep,
                question.vep_rs,
            ),
            _ => false,
        },
        _ => false,
    };
    if taken {
        Verdict::Class(ALTERNATE_CODON_PAST_THE_CODING_END, None)
    } else {
        Verdict::Pass
    }
}

#[cfg(test)]
#[path = "fields_tests.rs"]
mod tests;
