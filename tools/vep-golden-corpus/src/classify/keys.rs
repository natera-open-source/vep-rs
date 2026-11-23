// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The key classes: the class of a pair of consequence sets on one key, Ensembl VEP's
//! against vep-rs's. `classify` asks three readings in turn and the first that names a class
//! decides: the set rules of the tuple comparator ([`classify_pair`]), the mate-side read of
//! a point breakend ([`classify_mate_side_pair`]) and the row-aware rules of the tuple
//! comparator on the rows that carry the two sets ([`classify_stop_codon_rows`]). A pair whose
//! allele opens with its whole reference allele is put to the row-aware
//! rules also where a set rule classifies it without setting it aside, and is counted under
//! the row-aware rule that takes it ([`classify_counted_rows_after_the_span`]). A key whose
//! sets the three name under no class is asked the class of `super::intron_span`.
//!
//! The rules are the tuple comparator's own, `tools/vep-outputs-compare/src/rules/` compiled
//! into this crate as `crate::rules`, so a class is defined once. A rule's pair is named by
//! the taxonomy class of the rule's line in that tool's class file,
//! `tools/vep-outputs-compare/data/divergence_classes.tsv`.

use super::records::RunRecord;
use crate::constants::CORRECTED_BEHAVIOUR_CLASSES;
use crate::mate_context::{
    bracket_mate, mate_side_derivation, TranscriptModels, BREAKEND_MATE_LOCAL_READ,
};
use crate::rules::breakend::{strip_chr, vep_rs_pos};
use crate::rules::coding_end::PairSources;
use crate::rules::gvcf::perl_split;
use crate::rules::row::{
    may_be_written_after_the_span, pair_row, record_names_its_rows, replaced_reference,
    start_codon_stands, written_after_the_span, PairRow, RecordReadings,
};
use crate::rules::{
    classify as counting_bucket, excluded, row_pair_rules, EXCLUSION_REGISTRY,
    START_COOCCURRENCE_SWAP,
};
use std::collections::{BTreeSet, HashSet};
use std::sync::LazyLock;

/// A row that carries a tuple: its `cDNA_position`, its `CDS_position`, its
/// `Uploaded_variation`, and its `Protein_position`, `Amino_acids` and `Codons` where the row
/// has all twelve columns up to `Codons`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PositionRow {
    pub cdna_position: String,
    pub cds_position: String,
    pub uploaded_variation: String,
    pub cells: Option<[String; 3]>,
}

impl PositionRow {
    /// A row of the two positions and the first column, without its three further cells.
    pub fn of_positions(cdna_position: &str, cds_position: &str, uploaded_variation: &str) -> Self {
        PositionRow {
            cdna_position: cdna_position.to_string(),
            cds_position: cds_position.to_string(),
            uploaded_variation: uploaded_variation.to_string(),
            cells: None,
        }
    }
}

// The five classes of the stop codon tests behind `stop_retained` (Utils/VariationEffect.pm
// :1290-1324), the class of a stop codon an allele settles before its unknown bases, the
// class of a frameshift inside the codon of a stop the peptide holds, the class of a splice
// term read past the record's span and the classes of a start codon removed and of a start
// codon kept on a transcript with no 5' UTR, each the bucket of a row-aware rule of the tuple
// comparator, which defines it.
pub const STOP_RETAINED_WITHOUT_STOP_CODON_TEST: &str = "stop_retained_without_stop_codon_test";
pub const STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST: &str =
    "stop_codon_removed_without_stop_codon_test";
pub const STOP_CODON_OVERLAP_GENOMIC_SPAN: &str = "stop_codon_overlap_genomic_span";
pub const STOP_RETAINED_THREE_BASE_STOP_CODON_TEST: &str =
    "stop_retained_three_base_stop_codon_test";
pub const STOP_CODON_OVERLAP_INSERTION_RANGE: &str = "stop_codon_overlap_insertion_range";
pub const STOP_CODON_SETTLED_BEFORE_UNKNOWN_BASES: &str = "stop_codon_settled_before_unknown_bases";
pub const FRAMESHIFT_INSIDE_A_HELD_STOP_CODON: &str = "frameshift_inside_a_held_stop_codon";
pub const SPLICE_TERM_PAST_THE_RECORD_SPAN: &str = "splice_term_past_the_record_span";
pub const START_CODON_REMOVED_WITHOUT_FIVE_PRIME_UTR: &str =
    "start_codon_removed_without_five_prime_utr";
pub const START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR: &str = "start_codon_kept_without_five_prime_utr";

// The five classes of a coding end at which no stop codon or no start codon stands, or whose
// stop codon polyadenylation completes: each the bucket of a row-aware rule of the tuple
// comparator, which defines it and asks the coding ends of the row's transcript in the
// corpus's cache.
pub const STRUCTURAL_DELETION_STOP_LOST_WITHOUT_STOP_CODON: &str =
    "structural_deletion_stop_lost_without_stop_codon";
pub const STOP_CODON_COMPLETED_BY_POLYADENYLATION: &str = "stop_codon_completed_by_polyadenylation";
pub const STOP_CODON_OVERLAP_WITHOUT_STOP_CODON: &str = "stop_codon_overlap_without_stop_codon";
pub const STOP_CODON_BEHIND_LAST_RESIDUE_WITHOUT_STOP_CODON: &str =
    "stop_codon_behind_last_residue_without_stop_codon";
pub const START_CODON_OVERLAP_WITHOUT_START_CODON: &str = "start_codon_overlap_without_start_codon";

pub const STOP_LOST_TERM: &str = "stop_lost";
pub const STOP_RETAINED_TERM: &str = "stop_retained_variant";
pub const CODING_SEQUENCE_TERM: &str = "coding_sequence_variant";
pub const FRAMESHIFT_TERM: &str = "frameshift_variant";

/// What stands in vep-rs's set in the place of the reference's `stop_retained_variant` on a
/// record whose stop codon test is not made.
const UNTESTED_STOP_CODON_REPLACEMENTS: [&[&str]; 3] =
    [&[], &[CODING_SEQUENCE_TERM], &[FRAMESHIFT_TERM]];

/// What the reference writes where vep-rs writes the terms of [`RETAINED_INSERTION_TERMS`].
const FRAMESHIFT_REPLACEMENTS: [&[&str]; 2] =
    [&[FRAMESHIFT_TERM], &["stop_gained", FRAMESHIFT_TERM]];

/// What vep-rs writes for an insertion that keeps the stop codon.
const RETAINED_INSERTION_TERMS: [&[&str]; 2] = [
    &[STOP_RETAINED_TERM],
    &["inframe_insertion", STOP_RETAINED_TERM],
];

/// The registry of the class file whose rules read the two consequence sets of a pair.
pub const SET_RULES: &str = "pair";

/// The registry of the class file whose rules read the two rows of a pair.
pub const ROW_RULES: &str = "row_pair_116.2";

/// One class of the two pair registries of the tuple comparator's class file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PairClass {
    pub registry: &'static str,
    /// The name of the rule the class belongs to.
    pub bucket: &'static str,
    pub taxonomy_class: &'static str,
}

/// The tuple comparator's class file: a header line, then one tab-separated line a class.
const CLASS_FILE: &str = include_str!("../../../vep-outputs-compare/data/divergence_classes.tsv");

static PAIR_CLASSES: LazyLock<Vec<PairClass>> = LazyLock::new(|| {
    let mut lines = CLASS_FILE.lines();
    let header: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|held| *held == name)
            .unwrap_or_else(|| panic!("the class file has no column {name}"))
    };
    let columns = [
        column("registry"),
        column("bucket"),
        column("taxonomy_class"),
    ];
    lines
        .map(|line| {
            let fields: Vec<&'static str> = line.split('\t').collect();
            let [registry, bucket, taxonomy_class] =
                columns.map(|at| fields.get(at).copied().unwrap_or_default());
            PairClass {
                registry,
                bucket,
                taxonomy_class,
            }
        })
        .filter(|class| class.registry == SET_RULES || class.registry == ROW_RULES)
        .collect()
});

/// The classes of the two pair registries, in the order of the class file: the set rules,
/// then the row-aware rules, each registry in the order its rules are asked.
pub fn pair_classes() -> &'static [PairClass] {
    &PAIR_CLASSES
}

/// The taxonomy class of the rule `bucket` of `registry`.
fn taxonomy_class(registry: &str, bucket: &str) -> &'static str {
    pair_classes()
        .iter()
        .find(|class| class.registry == registry && class.bucket == bucket)
        .map(|class| class.taxonomy_class)
        .unwrap_or_else(|| panic!("the class file names no rule {bucket} of {registry}"))
}

/// The terms of a consequence set.
fn terms(csq: &str) -> BTreeSet<&str> {
    csq.split(',').collect()
}

/// Whether `set` is one of `shapes`.
fn is_one_of(set: &BTreeSet<&str>, shapes: &[&[&str]]) -> bool {
    shapes.iter().any(|shape| is_set(set, shape))
}

/// Whether `set` holds the terms of `shape` and no other.
fn is_set(set: &BTreeSet<&str>, shape: &[&str]) -> bool {
    set.len() == shape.len() && shape.iter().all(|term| set.contains(term))
}

/// The taxonomy class of the first set rule of the tuple comparator that takes Ensembl
/// VEP's set against vep-rs's, in the registry's order; the rule that sets pairs aside
/// without counting them is asked with the others.
pub fn classify_pair(perl_csq: &str, rust_csq: &str, allele: &str) -> Option<&'static str> {
    classify_pair_beside(perl_csq, rust_csq, allele, &PairSources::no_cache())
}

/// [`classify_pair`] on a key whose transcript the corpus's cache was asked for (`read`):
/// beside the cache the rule of `start_cooccurrence_swap` takes a pair only on a transcript
/// whose record shows a start codon, as the tuple comparator asks it.
pub fn classify_pair_beside(
    perl_csq: &str,
    rust_csq: &str,
    allele: &str,
    read: &PairSources<'_>,
) -> Option<&'static str> {
    EXCLUSION_REGISTRY
        .iter()
        .find(|rule| {
            rule.matches(perl_csq, rust_csq, allele)
                && (rule.bucket != START_COOCCURRENCE_SWAP || read.ask(start_codon_stands))
        })
        .map(|rule| taxonomy_class(SET_RULES, rule.bucket))
}

/// `breakend_mate_local_read` when Ensembl VEP's set is the local-coordinate read of a point
/// breakend's mate-side Transcript row and vep-rs's set is the mate-coordinate read, the two
/// differing; the derivation is `crate::mate_context`'s, over the corpus's own cache.
pub fn classify_mate_side_pair(
    location: &str,
    allele: &str,
    feature: &str,
    perl_csq: &str,
    rust_csq: &str,
    models: &TranscriptModels,
) -> Option<&'static str> {
    let (mate_chromosome, _) = bracket_mate(allele)?;
    let model = models.get(mate_chromosome, feature)?;
    let (mate_set, local_set) = mate_side_derivation(location, allele, model)?;
    if mate_set == local_set {
        return None;
    }
    (terms(perl_csq) == local_set && terms(rust_csq) == mate_set)
        .then_some(BREAKEND_MATE_LOCAL_READ)
}

/// The reference bases a row's `allele` stands in for at its `location`, read by position
/// from the alleles of a run's `records` as the tuple comparator reads them from an input
/// VCF: a record on the Location's chromosome, a `chr` prefix set aside, that writes the row
/// gives them (`replaced_reference`), the empty string for a pure insertion's row. Where two
/// or more records read a row differently, the record whose identifier is `name`, the row's
/// own first column (`Uploaded_variation`), decides; where no record carries that identifier,
/// or more than one does and they read the row differently, the row keeps none. `None` as
/// well where no record writes the row.
pub fn record_reference(
    records: &[RunRecord],
    location: &str,
    allele: &str,
    name: &str,
) -> Option<String> {
    let (chromosome, _) = location.split_once(':')?;
    let mut found = RecordReadings::default();
    for record in records {
        let fields: Vec<&str> = record.vcf.splitn(6, '\t').collect();
        if fields.len() < 5 || strip_chr(fields[0]) != strip_chr(chromosome) {
            continue;
        }
        let Some(position) = vep_rs_pos(fields[1]) else {
            continue;
        };
        let alternates = perl_split(fields[4], ',');
        let bases = replaced_reference(position, fields[3], &alternates, location, allele);
        if let Some(bases) = bases {
            found.add(bases, record_names_its_rows(fields[2], name));
        }
    }
    found.decided()
}

/// The class whose row-aware rule of the tuple comparator takes the pair of Transcript rows,
/// asked in the registry's order: the reference's row with its first column, its
/// `cDNA_position`, its `CDS_position` and, where it has all twelve columns, its
/// `Protein_position`, `Amino_acids` and `Codons`, and vep-rs's with its two positions and
/// the same three cells where it has them. The reference's row
/// carries the reference bases its allele stands in for where its first column spells the
/// alleles of a record that writes the row; else both rows carry them as the run's `records`
/// give them by position ([`record_reference`]), and neither where no record writes the row.
/// A rule is asked where its class is one a corpus documents rows under
/// (`CORRECTED_BEHAVIOUR_CLASSES`): `pachysolen_nuclear_code_read_as_standard` is none, so
/// its rule takes no pair here.
pub fn comparator_stop_codon_class(
    location: &str,
    allele: &str,
    perl_csq: &str,
    rust_csq: &str,
    reference_row: &PositionRow,
    candidate_row: &PositionRow,
    records: &[RunRecord],
) -> Option<&'static str> {
    comparator_coding_end_class(
        location,
        allele,
        perl_csq,
        rust_csq,
        reference_row,
        candidate_row,
        records,
        &PairSources::no_cache(),
    )
}

/// [`comparator_stop_codon_class`] on the rows as `classify` reads them, each with its
/// `Protein_position`, `Amino_acids` and `Codons` where it has them, in a run whose corpus
/// cache was asked for the transcript of the rows' key (`read`). A rule that asks the coding
/// ends of the transcript is asked with the reading of its cache record, and the one rule
/// that reads past a transcript's last base with the corpus's own reference sequence, as
/// the tuple comparator asks them: beside the cache a class of a coding end takes a pair
/// where vep-rs's coding terms are the ones derived for the row, and a stop codon class only
/// on a transcript whose record shows the stop codon at the allele's place.
#[allow(clippy::too_many_arguments)]
pub fn comparator_coding_end_class(
    location: &str,
    allele: &str,
    perl_csq: &str,
    rust_csq: &str,
    reference_row: &PositionRow,
    candidate_row: &PositionRow,
    records: &[RunRecord],
    read: &PairSources<'_>,
) -> Option<&'static str> {
    let transcript_row = |name: &str, csq: &str, row: &PositionRow| -> PairRow {
        let read = pair_row(&[
            name,
            location,
            allele,
            "",
            "",
            "Transcript",
            csq,
            &row.cdna_position,
            &row.cds_position,
        ]);
        match &row.cells {
            Some([protein_position, amino_acids, codons]) => {
                read.with_cells([protein_position, amino_acids, codons])
            }
            None => read,
        }
    };
    let mut reference = transcript_row(&reference_row.uploaded_variation, perl_csq, reference_row);
    let mut candidate = transcript_row("", rust_csq, candidate_row);
    if reference.reference.is_none() {
        reference.reference =
            record_reference(records, location, allele, &reference_row.uploaded_variation);
        candidate.reference = reference.reference.clone();
    }
    let (reference, candidate) = ([reference], [candidate]);
    read.ask(|at| {
        row_pair_rules()
            .into_iter()
            .map(|rule| (rule, taxonomy_class(ROW_RULES, rule.bucket)))
            .filter(|(_, class)| CORRECTED_BEHAVIOUR_CLASSES.contains(class))
            .find(|(rule, _)| rule.matches_rows_at(&reference, &candidate, at))
            .map(|(_, class)| class)
    })
}

/// The class of one reference row and one vep-rs row of a tuple whose two consequence sets
/// differ by a stop codon term, by a splice or intron term the reference alone holds or by
/// a `start_lost` or a `start_retained_variant` vep-rs alone holds, `reference_row` being
/// the reference's row, `candidate_row` vep-rs's, and `records` the run's records, whose
/// alleles a row under an identifier is read with.
///
/// A Transcript row takes the class whose rule in the tuple comparator takes the pair of
/// rows ([`comparator_stop_codon_class`]). A row that names no feature, the one row a
/// record has under `--summary` or `--most_severe`, carries the terms of all its features
/// and no position. It takes `stop_retained_without_stop_codon_test` or
/// `stop_codon_removed_without_stop_codon_test` when its sets differ in that class's way
/// (the reference's `stop_retained_variant` against nothing, `coding_sequence_variant` or
/// `frameshift_variant` for the first, against `stop_lost` for the second) and its Location
/// and Allele are in `untested_records`, the records whose Transcript rows the manifest's
/// own run documents under one of the two. It takes `stop_codon_overlap_insertion_range`
/// when its sets differ in that class's way (the reference's `coding_sequence_variant`, or
/// its `frameshift_variant` with or without `stop_gained`, against vep-rs's
/// `stop_retained_variant`, the second with or without `inframe_insertion`) and its Location
/// is in `retained_insertion_locations`, the records documented there under it. It takes
/// `stop_codon_settled_before_unknown_bases`, `frameshift_inside_a_held_stop_codon`,
/// `splice_term_past_the_record_span`, `start_codon_removed_without_five_prime_utr` and
/// `start_codon_kept_without_five_prime_utr` nowhere: the first, the third, the fourth and
/// the fifth read one transcript's CDS position, the first two its Codons and its residues,
/// the third its cDNA position as well, the fourth its cDNA position and its residues and
/// the fifth its cDNA position as well, which such a row does not carry. A row of any
/// other feature type, and a Transcript pair without either of its rows, takes no class:
/// the row rules read both rows, and the reference's positions say nothing of vep-rs's.
#[allow(clippy::too_many_arguments)]
pub fn classify_stop_codon_pair(
    location: &str,
    allele: &str,
    feature_type: &str,
    perl_csq: &str,
    rust_csq: &str,
    reference_row: Option<&PositionRow>,
    untested_records: &HashSet<(String, String)>,
    retained_insertion_locations: &HashSet<String>,
    candidate_row: Option<&PositionRow>,
    records: &[RunRecord],
) -> Option<&'static str> {
    if feature_type == "-" {
        let (perl, rust) = (terms(perl_csq), terms(rust_csq));
        let perl_only: BTreeSet<&str> = perl.difference(&rust).copied().collect();
        let rust_only: BTreeSet<&str> = rust.difference(&perl).copied().collect();
        if retained_insertion_locations.contains(location)
            && is_one_of(&rust_only, &RETAINED_INSERTION_TERMS)
            && (is_one_of(&perl_only, &FRAMESHIFT_REPLACEMENTS)
                || (is_set(&perl_only, &[CODING_SEQUENCE_TERM])
                    && is_set(&rust_only, &[STOP_RETAINED_TERM])))
        {
            return Some(STOP_CODON_OVERLAP_INSERTION_RANGE);
        }
        let record = (location.to_string(), allele.to_string());
        if !untested_records.contains(&record) || !is_set(&perl_only, &[STOP_RETAINED_TERM]) {
            return None;
        }
        if is_one_of(&rust_only, &UNTESTED_STOP_CODON_REPLACEMENTS) {
            return Some(STOP_RETAINED_WITHOUT_STOP_CODON_TEST);
        }
        if is_set(&rust_only, &[STOP_LOST_TERM]) {
            return Some(STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST);
        }
        return None;
    }
    if feature_type != "Transcript" {
        return None;
    }
    comparator_stop_codon_class(
        location,
        allele,
        perl_csq,
        rust_csq,
        reference_row?,
        candidate_row?,
        records,
    )
}

/// The class of the row-aware rule that takes a pair of Transcript rows the set rules of the
/// tuple comparator classify without setting it aside, where every reference row is that of
/// an allele that opens with its whole reference allele: the comparator
/// puts such a pair to its row-aware rules as well and counts it under the one that takes it.
/// `None` for every other pair, and where two pairs of rows are named differently.
#[allow(clippy::too_many_arguments)]
pub fn classify_counted_rows_after_the_span(
    location: &str,
    allele: &str,
    feature_type: &str,
    perl_csq: &str,
    rust_csq: &str,
    reference_rows: &[PositionRow],
    candidate_rows: &[PositionRow],
    records: &[RunRecord],
) -> Option<&'static str> {
    classify_counted_rows_after_the_span_beside(
        location,
        allele,
        feature_type,
        perl_csq,
        rust_csq,
        reference_rows,
        candidate_rows,
        records,
        &PairSources::no_cache(),
    )
}

/// [`classify_counted_rows_after_the_span`] in a run whose corpus cache was asked for the
/// transcript of the key (`read`), which the row-aware rules are asked with
/// ([`comparator_coding_end_class`]).
#[allow(clippy::too_many_arguments)]
pub fn classify_counted_rows_after_the_span_beside(
    location: &str,
    allele: &str,
    feature_type: &str,
    perl_csq: &str,
    rust_csq: &str,
    reference_rows: &[PositionRow],
    candidate_rows: &[PositionRow],
    records: &[RunRecord],
    read: &PairSources<'_>,
) -> Option<&'static str> {
    if feature_type != "Transcript"
        || excluded(perl_csq, rust_csq, allele)
        || counting_bucket(perl_csq, rust_csq, allele).is_none()
        || !may_be_written_after_the_span(location, allele)
    {
        return None;
    }
    let after_the_span = |row: &PositionRow| {
        let (cdna, cds, name) = (
            &row.cdna_position,
            &row.cds_position,
            &row.uploaded_variation,
        );
        let mut read = pair_row(&[
            name,
            location,
            allele,
            "",
            "",
            "Transcript",
            perl_csq,
            cdna,
            cds,
        ]);
        if read.reference.is_none() {
            read.reference = record_reference(records, location, allele, name);
        }
        written_after_the_span(&read, &read)
    };
    if !reference_rows.iter().all(after_the_span) {
        return None;
    }
    let mut classes = reference_rows.iter().flat_map(|reference| {
        candidate_rows.iter().map(move |candidate| {
            comparator_coding_end_class(
                location, allele, perl_csq, rust_csq, reference, candidate, records, read,
            )
        })
    });
    let first = classes.next()?;
    classes.all(|class| class == first).then_some(first)?
}

/// The class [`classify_stop_codon_pair`] names for every reference row that carries a tuple
/// against every vep-rs row that carries the other, as the tuple comparator's rule takes a
/// pair only where every such pair of rows meets it: `None` where a side has no row or two
/// pairs of rows are named differently. A row that names no feature carries no position and
/// is read by its terms alone. `records` are the run's records, whose alleles a row is read
/// with where its first column spells none ([`record_reference`]).
#[allow(clippy::too_many_arguments)]
pub fn classify_stop_codon_rows(
    location: &str,
    allele: &str,
    feature_type: &str,
    perl_csq: &str,
    rust_csq: &str,
    reference_rows: &[PositionRow],
    candidate_rows: &[PositionRow],
    untested_records: &HashSet<(String, String)>,
    retained_insertion_locations: &HashSet<String>,
    records: &[RunRecord],
) -> Option<&'static str> {
    classify_stop_codon_rows_beside(
        location,
        allele,
        feature_type,
        perl_csq,
        rust_csq,
        reference_rows,
        candidate_rows,
        untested_records,
        retained_insertion_locations,
        records,
        &PairSources::no_cache(),
    )
}

/// [`classify_stop_codon_rows`] in a run whose corpus cache was asked for the transcript of
/// the key (`read`): a pair of Transcript rows is put to the row-aware rules with their
/// five cells and that reading ([`comparator_coding_end_class`]). A row that names no feature
/// names no transcript, asks no record and is read as [`classify_stop_codon_pair`] reads
/// it: no class of a coding end takes it.
#[allow(clippy::too_many_arguments)]
pub fn classify_stop_codon_rows_beside(
    location: &str,
    allele: &str,
    feature_type: &str,
    perl_csq: &str,
    rust_csq: &str,
    reference_rows: &[PositionRow],
    candidate_rows: &[PositionRow],
    untested_records: &HashSet<(String, String)>,
    retained_insertion_locations: &HashSet<String>,
    records: &[RunRecord],
    read: &PairSources<'_>,
) -> Option<&'static str> {
    if feature_type == "Transcript" {
        let mut classes = reference_rows.iter().flat_map(|reference| {
            candidate_rows.iter().map(move |candidate| {
                comparator_coding_end_class(
                    location, allele, perl_csq, rust_csq, reference, candidate, records, read,
                )
            })
        });
        let first = classes.next()?;
        return classes.all(|class| class == first).then_some(first)?;
    }
    let class_of = |reference_row: Option<&PositionRow>,
                    candidate_row: Option<&PositionRow>,
                    records: &[RunRecord]| {
        classify_stop_codon_pair(
            location,
            allele,
            feature_type,
            perl_csq,
            rust_csq,
            reference_row,
            untested_records,
            retained_insertion_locations,
            candidate_row,
            records,
        )
    };
    if feature_type == "-" {
        return class_of(None, None, &[]);
    }
    let mut classes = reference_rows.iter().flat_map(|reference| {
        let class_of = &class_of;
        candidate_rows
            .iter()
            .map(move |candidate| class_of(Some(reference), Some(candidate), records))
    });
    let first = classes.next()?;
    classes.all(|class| class == first).then_some(first)?
}

#[cfg(test)]
#[path = "keys_tests.rs"]
mod tests;
