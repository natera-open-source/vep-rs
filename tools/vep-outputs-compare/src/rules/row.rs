// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row-aware pair rules: what they read of an output row (the row's own columns and,
//! where the first column spells the record's alleles, whether they make the row's allele an
//! insertion or a deletion, whether the record's reference allele is empty once the shared
//! first base is dropped, and the reference bases the allele stands in for) and their
//! predicates over the reference's and vep-rs's row of one pair. A row read to its twelfth
//! column carries its Protein_position, Amino_acids and Codons as well.

use super::cache_transcripts::location_span as location_bases;
use super::coding_end::{AtACodingEnd, CodingEndsRead};
use super::coding_ends::terms::{
    class_of_pair, landed_rule_may_take, Cells, CodingEndClass, LandedRule, RowFacts, WhatWasRead,
};
use regex::Regex;
use std::collections::BTreeSet;
use std::sync::LazyLock;

/// One output row as a row-aware pair rule reads it: its first column (`name`), the Location
/// and the Allele as the row writes them, its Feature_type, Consequence, cDNA_position and
/// CDS_position, and the
/// allele's shape from the record's own alleles, `None` where the first column spells none:
/// whether they make it an insertion or a deletion (`indel`), and whether every sequence
/// allele is an insertion, the reference allele empty once the shared first base is dropped
/// (`insertion`). `reference` is the reference bases the allele stands in for at the row's
/// Location, as the record's own alleles give them by position (`replaced_reference`): from
/// the first column where it spells the alleles, else from the input records where the
/// comparison is given them; empty for a pure insertion's row, `None` where no record is read.
/// `protein_position`, `amino_acids` and `codons` are the row's tenth to twelfth columns as
/// it writes them, each `None` for a row that ends before its twelfth column and for a row read
/// to its ninth column alone (`pair_row`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PairRow {
    pub name: String,
    pub location: String,
    pub allele: String,
    pub feature_type: String,
    pub consequence: String,
    pub cdna_position: String,
    pub cds_position: String,
    pub indel: Option<bool>,
    pub insertion: Option<bool>,
    pub reference: Option<String>,
    pub protein_position: Option<String>,
    pub amino_acids: Option<String>,
    pub codons: Option<String>,
}

impl PairRow {
    /// The row with its tenth to twelfth columns, Protein_position, Amino_acids and Codons.
    pub fn with_cells(self, cells: [&str; 3]) -> PairRow {
        let [protein_position, amino_acids, codons] = cells.map(|cell| Some(cell.to_string()));
        PairRow {
            protein_position,
            amino_acids,
            codons,
            ..self
        }
    }
}

/// The first nine columns of an output row, Uploaded_variation to CDS_position, in the form
/// a row-aware pair rule reads; the row carries none of its three further cells
/// (`PairRow::with_cells`).
pub fn pair_row(columns: &[&str; 9]) -> PairRow {
    let [name, location, allele, _, _, feature_type, consequence, cdna_position, cds_position] =
        *columns;
    let record_alleles = uploaded_alleles(name);
    let indel = record_alleles
        .as_ref()
        .and_then(|(reference, alternates)| allele_is_indel(reference, alternates, allele));
    let insertion = record_alleles
        .as_ref()
        .map(|(reference, alternates)| reference_allele_is_empty(reference, alternates));
    let reference = record_alleles.as_ref().and_then(|(reference, alternates)| {
        let position = uploaded_position(name)?;
        replaced_reference(position, reference, alternates, location, allele)
    });
    PairRow {
        name: name.to_string(),
        location: location.to_string(),
        allele: allele.to_string(),
        feature_type: feature_type.to_string(),
        consequence: consequence.to_string(),
        cdna_position: cdna_position.to_string(),
        cds_position: cds_position.to_string(),
        indel,
        insertion,
        reference,
        protein_position: None,
        amino_acids: None,
        codons: None,
    }
}

/// The first twelve columns of an output row, Uploaded_variation to Codons, in the form a
/// row-aware pair rule reads: the row of `pair_row` with its Protein_position, Amino_acids and
/// Codons, the tenth to twelfth columns, each as `PairRow::with_cells` sets it.
#[allow(dead_code)]
pub fn pair_row_with_codons(columns: &[&str; 12]) -> PairRow {
    let [first_nine @ .., protein_position, amino_acids, codons] = columns;
    PairRow {
        protein_position: Some(protein_position.to_string()),
        amino_acids: Some(amino_acids.to_string()),
        codons: Some(codons.to_string()),
        ..pair_row(first_nine)
    }
}

/// Ensembl VEP's Uploaded_variation for a VCF record without an ID (OutputFactory.pm:883-931):
/// the chromosome, the position and the record's own alleles joined by `/`. The classes name
/// white space as Python's `\s` does.
static UPLOADED_ALLELES: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^[^\s\x1c-\x1f]+_[0-9]+_([ACGTNacgtn*-]+(?:/[^/\s\x1c-\x1f]+)+)$").expect("regex")
});

/// The reference and alternate alleles an identifier-less Uploaded_variation spells, or
/// `None` for an identifier, which names no allele.
pub fn uploaded_alleles(uploaded: &str) -> Option<(&str, Vec<&str>)> {
    let spelt = UPLOADED_ALLELES.captures(uploaded)?.get(1)?.as_str();
    let mut parts = spelt.split('/');
    let reference = parts.next()?;
    Some((reference, parts.collect()))
}

/// The same name, the group its position.
static UPLOADED_POSITION: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^[^\s\x1c-\x1f]+_([0-9]+)_[ACGTNacgtn*-]+(?:/[^/\s\x1c-\x1f]+)+$").expect("regex")
});

/// The POS an identifier-less Uploaded_variation spells before its alleles, or `None` for an
/// identifier. A position past `i128`, or in digits other than `0` to `9`, is read as no
/// number.
pub fn uploaded_position(uploaded: &str) -> Option<i128> {
    let position = UPLOADED_POSITION.captures(uploaded)?.get(1)?;
    position.as_str().parse().ok()
}

/// The first and the last position a Location names where both are written in digits: the
/// text after the first `:`, parted at its first `-`, one position standing for both. A
/// position past `i128` is read as no number.
fn location_ends(location: &str) -> Option<(i128, i128)> {
    let span = location.split_once(':').map_or(location, |(_, span)| span);
    let (start, end) = span.split_once('-').unwrap_or((span, ""));
    let end = if end.is_empty() { start } else { end };
    let digits = |text: &str| !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit());
    if !(digits(start) && digits(end)) {
        return None;
    }
    Some((start.parse().ok()?, end.parse().ok()?))
}

/// Whether `spelt` is `head`, `written` and `tail` one after the other.
fn spells(spelt: &[char], head: &[char], written: &[char], tail: &[char]) -> bool {
    spelt.len() == head.len() + written.len() + tail.len()
        && spelt.starts_with(head)
        && spelt[head.len()..].starts_with(written)
        && spelt.ends_with(tail)
}

/// The reference bases an output row's allele stands in for, read by position from a record's
/// own alleles. The record whose reference allele `reference` starts at `position` writes
/// `allele` at `location` where one of its alternates is the reference allele with the
/// Location's bases replaced by the allele (`-` replaces them by nothing): those bases are
/// returned. Where an alternate is the reference allele with the allele inserted between the
/// Location's two neighbouring positions the row is a pure insertion's and the empty string
/// is returned. `None` where the record writes no such row, where two alternates read it
/// differently, and for an allele or a Location that is not a sequence's. The record's
/// alleles are read in upper case.
pub fn replaced_reference(
    position: i128,
    reference: &str,
    alternates: &[&str],
    location: &str,
    allele: &str,
) -> Option<String> {
    if allele != "-" && !spelt_in(allele, b"ACGTN") {
        return None;
    }
    let (start, end) = location_ends(location)?;
    let (first, last) = (start - position, end - position);
    let written: Vec<char> = match allele {
        "-" => Vec::new(),
        bases => bases.chars().collect(),
    };
    let reference: Vec<char> = reference.to_uppercase().chars().collect();
    let length = reference.len() as i128;
    let mut found: Vec<String> = Vec::new();
    let mut find = |bases: String| {
        if !found.contains(&bases) {
            found.push(bases);
        }
    };
    for alternate in alternates {
        let spelt: Vec<char> = alternate.to_uppercase().chars().collect();
        if 0 <= first && first <= last && last < length {
            let (first, last) = (first as usize, last as usize);
            if spells(
                &spelt,
                &reference[..first],
                &written,
                &reference[last + 1..],
            ) {
                find(reference[first..=last].iter().collect());
            }
        }
        if !written.is_empty() && last - first == 1 && 0 <= first && first < length {
            let after = first as usize + 1;
            if spells(&spelt, &reference[..after], &written, &reference[after..]) {
                find(String::new());
            }
        }
    }
    match found.len() {
        1 => found.pop(),
        _ => None,
    }
}

/// Whether the VCF record whose ID column is `ids` carries the identifier `name`, a row's
/// first column: the column is that text byte for byte. A record whose ID is `.` names no
/// row.
pub fn record_names_its_rows(ids: &str, name: &str) -> bool {
    ids != "." && ids == name
}

/// The readings the records of an input give one row ([`replaced_reference`]), and those
/// of the records whose rows carry the row's own first column.
#[derive(Clone, Debug, Default)]
pub struct RecordReadings {
    every: Vec<String>,
    own: Vec<String>,
}

impl RecordReadings {
    /// Adds the bases one record gives the row, `own` for a record whose rows carry the
    /// row's first column ([`record_names_its_rows`]).
    pub fn add(&mut self, bases: String, own: bool) {
        if own && !self.own.contains(&bases) {
            self.own.push(bases.clone());
        }
        if !self.every.contains(&bases) {
            self.every.push(bases);
        }
    }

    /// The row's bases: the one reading of its records where they read it alike. Where two
    /// or more records read a row differently, the record whose identifier is the row's own
    /// first column (`Uploaded_variation`) decides; where no record carries that identifier,
    /// or more than one does and they read the row differently, the row keeps none.
    pub fn decided(mut self) -> Option<String> {
        match (self.every.len(), self.own.len()) {
            (1, _) => self.every.pop(),
            (_, 1) => self.own.pop(),
            _ => None,
        }
    }
}

fn stripped_len(allele: &str) -> usize {
    allele.trim_matches('-').chars().count()
}

fn first_char(allele: &str) -> Option<char> {
    allele.chars().next()
}

fn without_first_char(allele: &str) -> &str {
    let mut chars = allele.chars();
    chars.next();
    chars.as_str()
}

/// Whether `allele` is an insertion or a deletion of the record whose reference allele is
/// `reference` and alternate alleles `alternates`, both as a row's first column spells them:
/// the alternate the allele names and the reference differ in length. A record with one
/// alternate not spelt with a `*` has that alternate's shape on every row; with several, the
/// shared first base is dropped where Ensembl VEP drops it (`Parser/VCF.pm:300-323`). `None`
/// when the allele is not a sequence, when no alternate names it, or when two alternates of
/// different shapes name it.
pub fn allele_is_indel(reference: &str, alternates: &[&str], allele: &str) -> Option<bool> {
    if allele.is_empty() || !allele.bytes().all(|b| b"ACGTN-".contains(&b)) {
        return None;
    }
    let trimmed: Vec<&str> = alternates
        .iter()
        .copied()
        .filter(|a| !a.contains('*'))
        .collect();
    if let [only] = trimmed[..] {
        return Some(stripped_len(reference) != stripped_len(only));
    }
    let reference_len = reference.chars().count();
    let holds_indel = alternates
        .iter()
        .any(|a| a.chars().count() != reference_len);
    let trim = !reference.is_empty()
        && holds_indel
        && trimmed
            .iter()
            .all(|a| first_char(a) == first_char(reference));
    let kept_reference = if trim {
        without_first_char(reference)
    } else {
        reference
    };
    let mut shapes: Vec<bool> = Vec::new();
    for alternate in alternates {
        let alt = if trim && !alternate.contains('*') {
            without_first_char(alternate)
        } else {
            alternate
        };
        if (if alt.is_empty() { "-" } else { alt }) == allele {
            let shape = stripped_len(kept_reference) != stripped_len(alt);
            if !shapes.contains(&shape) {
                shapes.push(shape);
            }
        }
    }
    match shapes[..] {
        [shape] => Some(shape),
        _ => None,
    }
}

/// Whether the reference allele of a record is empty once the record's shared first base is
/// dropped where Ensembl VEP drops it (`Parser/VCF.pm:300-323`), on the alleles as a row's
/// first column spells them. Every sequence allele of such a record is an insertion.
pub fn reference_allele_is_empty(reference: &str, alternates: &[&str]) -> bool {
    let reference_len = reference.chars().count();
    let holds_indel = alternates
        .iter()
        .any(|a| a.chars().count() != reference_len);
    let trim = !reference.is_empty()
        && holds_indel
        && alternates
            .iter()
            .filter(|a| !a.contains('*'))
            .all(|a| first_char(a) == first_char(reference));
    stripped_len(if trim {
        without_first_char(reference)
    } else {
        reference
    }) == 0
}

const CODING_SEQUENCE_TERM: &str = "coding_sequence_variant";
const STOP_RETAINED_TERM: &str = "stop_retained_variant";
const SPLICE_REGION_TERM: &str = "splice_region_variant";

/// The terms of a row's Consequence as the predicates below read them: parted at `&` and
/// `,`, the empty ones dropped, an empty Consequence reading as one empty term.
fn terms(consequence: &str) -> BTreeSet<&str> {
    let terms: BTreeSet<&str> = consequence
        .split(['&', ','])
        .filter(|term| !term.is_empty())
        .collect();
    if terms.is_empty() {
        BTreeSet::from([""])
    } else {
        terms
    }
}

/// The terms the reference's row holds, those only it holds and those only vep-rs's holds.
fn term_differences<'a>(
    e: &'a PairRow,
    a: &'a PairRow,
) -> (BTreeSet<&'a str>, BTreeSet<&'a str>, BTreeSet<&'a str>) {
    let (expected, actual) = (terms(&e.consequence), terms(&a.consequence));
    let gone = expected.difference(&actual).copied().collect();
    let come = actual.difference(&expected).copied().collect();
    (expected, gone, come)
}

fn is(terms: &BTreeSet<&str>, members: &[&str]) -> bool {
    terms.len() == members.len() && members.iter().all(|member| terms.contains(member))
}

/// Whether a Feature_type is a transcript's: read after a last `|` and up to a `#`, in any
/// letter case.
pub fn transcript_feature_type(feature_type: &str) -> bool {
    let kind = feature_type.rsplit('|').next().unwrap_or("");
    let kind = kind.split('#').next().unwrap_or("");
    kind.to_lowercase() == "transcript"
}

/// Whether the row is a Transcript row (`transcript_feature_type`).
fn transcript_row(e: &PairRow) -> bool {
    transcript_feature_type(&e.feature_type)
}

/// What the derivation of a row's coding terms (`coding_ends::terms`) reads of a pair: the
/// reference's row as the pair's own facts, `None` where its Location names no position.
fn row_facts(e: &PairRow) -> Option<RowFacts<'_>> {
    Some(RowFacts {
        transcript_row: transcript_row(e),
        allele: &e.allele,
        location: location_bases(&e.location)?,
        ensembl_vep_terms: &e.consequence,
    })
}

/// The five cells of a row that has all twelve columns up to Codons.
fn cells(row: &PairRow) -> Option<Cells<'_>> {
    Some(Cells {
        cdna_position: &row.cdna_position,
        cds_position: &row.cds_position,
        protein_position: row.protein_position.as_deref()?,
        amino_acids: row.amino_acids.as_deref()?,
        codons: row.codons.as_deref()?,
    })
}

/// What a comparison read beside the pair's key and terms: the five cells of both rows where
/// both carry them, the count of reference bases the allele stands in for, and the reference
/// sequence it was given. The count is 0 where a row's first column spells a record whose
/// reference allele is empty (`PairRow::insertion`), whatever an input record reads at the
/// same Location under the same Allele; else the length of the bases a record gave either
/// row (`PairRow::reference`).
fn what_was_read<'a>(e: &'a PairRow, a: &'a PairRow, at: &AtACodingEnd<'a>) -> WhatWasRead<'a> {
    let reference = e.reference.as_ref().or(a.reference.as_ref());
    let insertion = e.insertion == Some(true) || a.insertion == Some(true);
    WhatWasRead {
        cells: cells(e).zip(cells(a)),
        reference_length: match insertion {
            true => Some(0),
            false => reference.map(|bases| bases.chars().count() as u64),
        },
        reference_bases: at.reference_bases,
        // A row's reference bases are read from a VCF record, which is given on the forward
        // strand.
        record_on_the_reverse_strand: false,
    }
}

/// Whether the pair is the class `class`'s: the first class of a coding end, in the order
/// they are asked, that takes the two rows on the reading of the transcript's cache record
/// (`coding_ends::terms::class_of_pair`). A class takes a pair only where vep-rs's coding
/// terms are the ones derived for the row from that reading, so a comparison that was given
/// no cache, and a pair the cache answers with no one record, are no class's.
fn coding_end_class_pair(
    class: CodingEndClass,
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    let CodingEndsRead::Record(ends) = at.ends else {
        return false;
    };
    let Some(row) = row_facts(e) else {
        return false;
    };
    let read = what_was_read(e, a, at);
    class_of_pair(&row, &a.consequence, Some(ends), &read) == Some(class)
}

/// `structural_deletion_stop_lost_without_stop_codon`: on a transcript at whose coding end
/// no stop codon stands, Ensembl VEP writes `stop_lost` for a structural deletion over the
/// last three genomic positions of the coding region (Utils/VariationEffect.pm:1268-1280)
/// and vep-rs writes the same row without it (`CodingEndClass::takes`).
pub fn structural_deletion_stop_lost_without_stop_codon_pair(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    coding_end_class_pair(CodingEndClass::StructuralDeletionStopLost, e, a, at)
}

/// `stop_codon_completed_by_polyadenylation`: on a transcript whose stop codon
/// polyadenylation completes, Ensembl VEP reads the row on an incomplete terminal codon
/// (Utils/VariationEffect.pm:1588-1614) and vep-rs on the completed coding sequence
/// (`CodingEndClass::takes`).
pub fn stop_codon_completed_by_polyadenylation_pair(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    coding_end_class_pair(CodingEndClass::StopCodonCompleted, e, a, at)
}

/// `stop_codon_overlap_without_stop_codon`: on a transcript at whose coding end no stop
/// codon stands and that carries no `cds_end_NF`, Ensembl VEP writes a stop term from the
/// overlap of the allele with the last three coding positions
/// (Utils/VariationEffect.pm:1373-1381, 1400, 1408-1427) and vep-rs the terms of the row
/// without one (`CodingEndClass::takes`).
pub fn stop_codon_overlap_without_stop_codon_pair(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    coding_end_class_pair(CodingEndClass::StopCodonOverlap, e, a, at)
}

/// `stop_codon_behind_last_residue_without_stop_codon`: on a transcript at whose coding end
/// no stop codon stands, Ensembl VEP reads a stop codon an allele brings behind the last
/// residue as the transcript's own kept (Utils/VariationEffect.pm:1344, 1357) and vep-rs
/// as gained (`CodingEndClass::takes`).
pub fn stop_codon_behind_last_residue_without_stop_codon_pair(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    coding_end_class_pair(CodingEndClass::StopCodonBehindLastResidue, e, a, at)
}

/// `start_codon_overlap_without_start_codon`: on a transcript at whose coding start no start
/// codon is shown and that carries no `cds_start_NF`, Ensembl VEP writes a start term from
/// the overlap of the allele with the first three coding positions
/// (Utils/VariationEffect.pm:975-986) and vep-rs the terms of the row without one
/// (`CodingEndClass::takes`).
pub fn start_codon_overlap_without_start_codon_pair(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    coding_end_class_pair(CodingEndClass::StartCodonOverlap, e, a, at)
}

/// Whether a rule that stood before the coding ends were read may take the pair of the
/// reference's row `e` and vep-rs's row `a`, beside its own reading of the two rows. In a
/// comparison that was given no cache it stands as it stood. Beside a cache it takes a pair
/// only where the transcript's record shows the codon the rule stands on
/// (`coding_ends::terms::landed_rule_may_take`), and nothing where the cache answers the
/// pair with no one record. Where the allele lies is read by the rule's own predicate alone.
fn stands_beside_a_cache(
    rule: LandedRule,
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    let ends = match at.ends {
        CodingEndsRead::NoCache => return true,
        CodingEndsRead::NoRecord => None,
        CodingEndsRead::Record(ends) => Some(ends),
    };
    let Some(row) = row_facts(e) else {
        return false;
    };
    landed_rule_may_take(rule, &row, ends, &what_was_read(e, a, at))
}

/// Whether the set rule `start_cooccurrence_swap` stands on the transcript of a pair whose
/// two consequence sets it takes: in a comparison that was given no cache it stands as it
/// stood, and beside a cache only where the transcript's record shows a start codon at its
/// coding start (`coding_ends::terms::landed_rule_may_take`). vep-rs writes a start codon's
/// terms nowhere else, so a pair of the rule's shape on any other transcript is not the
/// rule's.
pub fn start_codon_stands(at: &AtACodingEnd<'_>) -> bool {
    let ends = match at.ends {
        CodingEndsRead::NoCache => return true,
        CodingEndsRead::NoRecord => None,
        CodingEndsRead::Record(ends) => Some(ends),
    };
    // The rule of a start codon asks the record alone: no fact of the row decides it.
    let row = RowFacts {
        transcript_row: true,
        allele: "",
        location: (0, 0),
        ensembl_vep_terms: "",
    };
    let read = WhatWasRead {
        cells: None,
        reference_length: None,
        reference_bases: None,
        record_on_the_reverse_strand: false,
    };
    landed_rule_may_take(LandedRule::StartCooccurrenceSwap, &row, ends, &read)
}

fn spelt_in(allele: &str, letters: &[u8]) -> bool {
    !allele.is_empty() && allele.bytes().all(|b| letters.contains(&b))
}

/// The number of bases a Location spans where its start and end are both written in digits:
/// the text after the first `:`, parted at its first `-`. A coordinate past `i128` is read as
/// no number.
fn location_span(location: &str) -> Option<i128> {
    let span = location.split_once(':').map_or(location, |(_, span)| span);
    let (start, end) = span.split_once('-')?;
    let digits = |text: &str| !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit());
    if !(digits(start) && digits(end)) {
        return None;
    }
    let (start, end): (i128, i128) = (start.parse().ok()?, end.parse().ok()?);
    Some((end - start).saturating_add(1))
}

/// Whether a row is an insertion's or a deletion's:
/// an allele of `-` is one; else the record's own alleles decide where the first column
/// spells them; else the Location decides, an indel's spanning another number of bases than
/// its allele of `A`, `C`, `G`, `T` and `N` has.
fn indel_row(e: &PairRow) -> bool {
    if e.allele == "-" {
        return true;
    }
    if let Some(indel) = e.indel {
        return indel;
    }
    spelt_in(&e.allele, b"ACGTN")
        && location_span(&e.location).is_some_and(|span| span != e.allele.len() as i128)
}

/// Whether a position column carries its start and its end: `?` stands for an end that maps
/// to no transcript base, one number for both, `-` for neither, each read up to a `/total`
/// suffix (OutputFactory.pm:1735-1748).
fn position_bounds(position: &str) -> (bool, bool) {
    let position = position.split_once('/').map_or(position, |(head, _)| head);
    if matches!(position, "" | "-" | "<absent>") {
        return (false, false);
    }
    match position.split_once('-') {
        Some((start, end)) => (start != "?", end != "?"),
        None => (position != "?", position != "?"),
    }
}

/// Whether the row carries the coordinates the stop codon test edits at
/// (Utils/VariationEffect.pm:1507-1508): both cDNA positions and the CDS start.
fn stop_codon_test_coordinates(e: &PairRow) -> bool {
    position_bounds(&e.cdna_position) == (true, true) && position_bounds(&e.cds_position).0
}

/// Whether the row shows that `_ins_del_stop_altered_cil` (Utils/VariationEffect.pm:1488-1539)
/// returns before it edits anything: the allele holds a letter other than `A`, `C`, `G` and
/// `T` (:1498), or the span lacks one of the test's coordinates (:1508).
fn stop_codon_test_not_made(e: &PairRow) -> bool {
    (e.allele != "-" && !spelt_in(&e.allele, b"ACGT")) || !stop_codon_test_coordinates(e)
}

/// `stop_codon_overlap_genomic_span`: a Transcript row of an insertion or deletion with both
/// cDNA positions and a CDS start on both sides, `splice_region_variant` in the reference's
/// terms, and the two sets equal but for the reference's `coding_sequence_variant` against
/// vep-rs's `stop_retained_variant` (Utils/VariationEffect.pm:1389-1431).
pub fn stop_codon_overlap_genomic_span_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || !indel_row(e) {
        return false;
    }
    let (expected, gone, come) = term_differences(e, a);
    is(&gone, &[CODING_SEQUENCE_TERM])
        && is(&come, &[STOP_RETAINED_TERM])
        && expected.contains(SPLICE_REGION_TERM)
        && !stop_codon_test_not_made(e)
        && !stop_codon_test_not_made(a)
}

/// `stop_codon_overlap_genomic_span_pair` in a comparison that may have been given the candidate's cache
/// (`stands_beside_a_cache`): beside a cache a whole stop codon stands or polyadenylation
/// completes one.
pub fn stop_codon_overlap_genomic_span_pair_at(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    stop_codon_overlap_genomic_span_pair(e, a)
        && stands_beside_a_cache(LandedRule::StopCodonOverlapGenomicSpan, e, a, at)
}

const FRAMESHIFT_TERM: &str = "frameshift_variant";

/// Whether a row is an insertion's or a deletion's as `indel_row` reads one, an allele of
/// other letters than `A`, `C`, `G`, `T` and `N` included: where the record's own alleles say
/// nothing, such an allele against a Location of another length.
fn length_changing_row(e: &PairRow) -> bool {
    if indel_row(e) {
        return true;
    }
    e.indel.is_none()
        && !e.allele.is_empty()
        && e.allele.bytes().all(|b| b.is_ascii_alphabetic())
        && location_span(&e.location).is_some_and(|span| span != e.allele.len() as i128)
}

/// `stop_retained_without_stop_codon_test`: a Transcript row of an insertion or deletion,
/// the reference's terms holding `stop_retained_variant` and vep-rs's being the same without
/// it, with nothing, `coding_sequence_variant` or `frameshift_variant` in its place, and both
/// rows showing that the stop codon test is not made (Utils/VariationEffect.pm:1318,
/// 1488-1539).
pub fn stop_retained_without_stop_codon_test_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || !length_changing_row(e) {
        return false;
    }
    let (_, gone, come) = term_differences(e, a);
    is(&gone, &[STOP_RETAINED_TERM])
        && (come.is_empty() || is(&come, &[CODING_SEQUENCE_TERM]) || is(&come, &[FRAMESHIFT_TERM]))
        && stop_codon_test_not_made(e)
        && stop_codon_test_not_made(a)
}

/// `stop_retained_without_stop_codon_test_pair` in a comparison that may have been given the candidate's cache
/// (`stands_beside_a_cache`): beside a cache a whole stop codon stands and the genome holds
/// each of its three bases.
pub fn stop_retained_without_stop_codon_test_pair_at(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    stop_retained_without_stop_codon_test_pair(e, a)
        && stands_beside_a_cache(LandedRule::StopRetainedWithoutStopCodonTest, e, a, at)
}

const STOP_GAINED_TERM: &str = "stop_gained";
/// The terms `frameshift` (Utils/VariationEffect.pm:1541-1567) reads away where
/// `stop_retained` does not hold.
const FRAMESHIFT_READ_TERMS: [&str; 3] = [
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
];

/// `stop_retained_three_base_stop_codon_test`: a Transcript row of an insertion or deletion
/// whose stop codon test is made on both sides, the reference's terms holding
/// `stop_retained_variant` and vep-rs's `coding_sequence_variant` in its place, or, where the
/// reference's hold it alone or beside one of `inframe_insertion`, `inframe_deletion` and
/// `protein_altering_variant`, `frameshift_variant` with or without `stop_gained` in the place
/// of those terms (Utils/VariationEffect.pm:1488-1539).
pub fn stop_retained_three_base_stop_codon_test_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || !length_changing_row(e) {
        return false;
    }
    let (_, gone, come) = term_differences(e, a);
    let terms_stand = if is(&come, &[CODING_SEQUENCE_TERM]) {
        is(&gone, &[STOP_RETAINED_TERM])
    } else if is(&come, &[FRAMESHIFT_TERM]) || is(&come, &[STOP_GAINED_TERM, FRAMESHIFT_TERM]) {
        gone.contains(STOP_RETAINED_TERM)
            && gone.len() <= 2
            && gone
                .iter()
                .all(|term| *term == STOP_RETAINED_TERM || FRAMESHIFT_READ_TERMS.contains(term))
    } else {
        false
    };
    terms_stand && !stop_codon_test_not_made(e) && !stop_codon_test_not_made(a)
}

/// `stop_retained_three_base_stop_codon_test_pair` in a comparison that may have been given the candidate's cache
/// (`stands_beside_a_cache`): beside a cache a whole stop codon stands.
pub fn stop_retained_three_base_stop_codon_test_pair_at(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    stop_retained_three_base_stop_codon_test_pair(e, a)
        && stands_beside_a_cache(LandedRule::StopRetainedThreeBaseStopCodonTest, e, a, at)
}

const STOP_LOST_TERM: &str = "stop_lost";

/// `stop_codon_removed_without_stop_codon_test`: a Transcript row of an insertion or
/// deletion, both rows lacking a cDNA position at one end or a CDS start, and the two sets of
/// terms equal but for the reference's `stop_retained_variant` against vep-rs's `stop_lost`
/// (Utils/VariationEffect.pm:1234-1288, 1508).
pub fn stop_codon_removed_without_stop_codon_test_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || !length_changing_row(e) {
        return false;
    }
    let (_, gone, come) = term_differences(e, a);
    is(&gone, &[STOP_RETAINED_TERM])
        && is(&come, &[STOP_LOST_TERM])
        && !stop_codon_test_coordinates(e)
        && !stop_codon_test_coordinates(a)
}

/// `stop_codon_removed_without_stop_codon_test_pair` in a comparison that may have been given the candidate's cache
/// (`stands_beside_a_cache`): beside a cache a whole stop codon stands or polyadenylation
/// completes one.
pub fn stop_codon_removed_without_stop_codon_test_pair_at(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    stop_codon_removed_without_stop_codon_test_pair(e, a)
        && stands_beside_a_cache(LandedRule::StopCodonRemovedWithoutStopCodonTest, e, a, at)
}

const INFRAME_INSERTION_TERM: &str = "inframe_insertion";

/// Whether a row's allele is a sequence of `A`, `C`, `G`, `T` and `N` and the record's own
/// alleles mark it an indel where the first column spells them, else the Location spans the
/// two bases around an insertion's point.
fn insertion_row(e: &PairRow) -> bool {
    if !spelt_in(&e.allele, b"ACGTN") {
        return false;
    }
    match e.indel {
        Some(indel) => indel,
        None => location_span(&e.location) == Some(2),
    }
}

/// Whether a row `insertion_row` takes is an insertion between two neighbouring positions.
/// Its Location is two neighbouring positions, and the first of these that applies decides:
///
/// - its reference bases (`PairRow::reference`) are empty: a pure insertion's row, whatever
///   the allele's length and however the record pads its alleles;
/// - the first column spells the alleles: the record's reference allele is empty once the
///   shared first base is dropped;
/// - the input records give the row bases that are not empty: a replacement's row, whose
///   start `_overlaps_stop_codon_cil` does not move (Utils/VariationEffect.pm:1408-1415);
/// - no record gives the row its bases: its allele has another length than two.
///
/// Under an identifier without the input records, two bases replaced by a sequence of
/// another length are therefore taken with the insertions and an insertion of two bases is
/// refused. The input records give a row its bases by Location and Allele: an insertion no
/// record places is refused where a record replaces the same two bases by the same allele.
fn inserted_between_neighbours(e: &PairRow) -> bool {
    if location_span(&e.location) != Some(2) {
        return false;
    }
    match (e.reference.as_deref(), e.insertion) {
        (Some(""), _) => true,
        (_, Some(insertion)) => insertion,
        (Some(_), None) => false,
        (None, None) => e.allele.len() != 2,
    }
}

/// `stop_codon_overlap_insertion_range`: a Transcript row of an insertion between two
/// neighbouring positions (`inserted_between_neighbours`) with both cDNA
/// positions and a CDS start on both sides, the two sets of terms equal but for the
/// reference's `coding_sequence_variant` against vep-rs's `stop_retained_variant`, or the
/// reference's `frameshift_variant` with or without `stop_gained` against vep-rs's
/// `stop_retained_variant` with or without `inframe_insertion`
/// (Utils/VariationEffect.pm:1389-1431). Beside `splice_region_variant` the first pair is
/// `stop_codon_overlap_genomic_span_pair`'s where the allele is of `A`, `C`, `G` and `T`
/// alone, and this rule's where it holds another letter, the rows that one refuses: no
/// peptide is read for such an allele (TranscriptVariationAllele.pm:701), so
/// `stop_retained` takes its second arm (Utils/VariationEffect.pm:1318) and the moved start
/// fails the overlap there too.
pub fn stop_codon_overlap_insertion_range_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || !insertion_row(e) || !inserted_between_neighbours(e) {
        return false;
    }
    let (expected, gone, come) = term_differences(e, a);
    let terms_stand = if is(&gone, &[CODING_SEQUENCE_TERM]) {
        let unknown_base = [e, a]
            .iter()
            .any(|row| row.allele != "-" && !spelt_in(&row.allele, b"ACGT"));
        is(&come, &[STOP_RETAINED_TERM]) && (!expected.contains(SPLICE_REGION_TERM) || unknown_base)
    } else if is(&gone, &[FRAMESHIFT_TERM]) || is(&gone, &[STOP_GAINED_TERM, FRAMESHIFT_TERM]) {
        is(&come, &[STOP_RETAINED_TERM]) || is(&come, &[INFRAME_INSERTION_TERM, STOP_RETAINED_TERM])
    } else {
        false
    };
    terms_stand && stop_codon_test_coordinates(e) && stop_codon_test_coordinates(a)
}

/// `stop_codon_overlap_insertion_range_pair` in a comparison that may have been given the candidate's cache
/// (`stands_beside_a_cache`): beside a cache a whole stop codon stands or polyadenylation
/// completes one.
pub fn stop_codon_overlap_insertion_range_pair_at(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    stop_codon_overlap_insertion_range_pair(e, a)
        && stands_beside_a_cache(LandedRule::StopCodonOverlapInsertionRange, e, a, at)
}

/// The stop codons of the standard genetic code.
const STANDARD_STOP_CODONS: [&[u8]; 3] = [b"TAA", b"TAG", b"TGA"];

/// The place of the stop codon the alternate half of a Codons cell settles before its
/// unknown bases, counted from 0 among the half's whole codons: the half is the text after
/// the cell's first `/`, read in upper case and cut into threes from its start, and the
/// stop codon is the first of the standard genetic code among the codons before the first
/// codon that holds a letter other than `A`, `C`, `G` and `T`. `None` for a cell without a
/// `/`, for a half that holds a character that is no letter, for a half of `A`, `C`, `G`
/// and `T` alone and for a half with no stop codon before that codon.
fn settled_stop_codon(codons: &str) -> Option<usize> {
    let (_, alternate) = codons.split_once('/')?;
    let alternate = alternate.to_ascii_uppercase();
    let bases = alternate.as_bytes();
    let known = |base: &u8| b"ACGT".contains(base);
    if !bases.iter().all(u8::is_ascii_alphabetic) || bases.iter().all(known) {
        return None;
    }
    bases
        .chunks_exact(3)
        .take_while(|codon| codon.iter().all(known))
        .position(|codon| STANDARD_STOP_CODONS.contains(&codon))
}

/// Whether the terms of the reference's row `e` and of vep-rs's row `a` differ as those of
/// a stop codon the inserted bases keep, given the reference's terms, those it alone holds
/// and those vep-rs alone holds. vep-rs alone holds `inframe_insertion`, which needs both
/// peptides (Utils/VariationEffect.pm:1113-1115) and cuts the alternate one after its first
/// stop (:1126). Where the reference holds `stop_retained_variant`, from the second arm of
/// `stop_retained` (:1318), whose edit test is not made for a letter that names no base
/// (:1498), nothing else differs. Otherwise vep-rs alone holds `stop_retained_variant` as
/// well, the reference alone `frameshift_variant` (:1564) or `coding_sequence_variant`
/// (:1616-1649), and the row is an insertion between two neighbouring positions
/// (`inserted_between_neighbours`) with the stop codon test's coordinates on both sides,
/// whose moved start fails the reference's overlap (:1408-1415).
fn kept_stop_terms(
    e: &PairRow,
    a: &PairRow,
    expected: &BTreeSet<&str>,
    gone: &BTreeSet<&str>,
    come: &BTreeSet<&str>,
) -> bool {
    if expected.contains(STOP_RETAINED_TERM) {
        return gone.is_empty() && is(come, &[INFRAME_INSERTION_TERM]);
    }
    (is(gone, &[FRAMESHIFT_TERM]) || is(gone, &[CODING_SEQUENCE_TERM]))
        && is(come, &[INFRAME_INSERTION_TERM, STOP_RETAINED_TERM])
        && insertion_row(e)
        && inserted_between_neighbours(e)
        && stop_codon_test_coordinates(e)
        && stop_codon_test_coordinates(a)
}

/// `stop_codon_settled_before_unknown_bases`: a Transcript row with one Allele and one CDS
/// position on both sides, the position carrying both ends, and one Codons cell on both
/// sides whose alternate half settles a stop codon before its unknown bases
/// (`settled_stop_codon`). The reference's Amino_acids is `-`, and the alternate half of
/// vep-rs's, the text after its first `/`, holds its first `*` at that stop codon's place.
/// The terms differ by a gained stop or by a kept one. For a gained stop, where the
/// reference's terms hold `frameshift_variant`, vep-rs's are the same terms with
/// `stop_gained` added; otherwise the reference alone writes `coding_sequence_variant` and
/// vep-rs alone `stop_gained`, with at most one of `inframe_insertion`, `inframe_deletion`
/// and `protein_altering_variant`. The reference builds the codons of an allele that holds
/// a letter other than `A`, `C`, `G` and `T` (TranscriptVariationAllele.pm:815) and reads no
/// peptide for it (:701), so `stop_gained` (Utils/VariationEffect.pm:1211-1232) is false and
/// `coding_unknown` (:1616-1649) holds where the length is whole codons; vep-rs translates
/// the allele where a stop codon stands before its first codon with such a letter. A row
/// read without its Codons is no row of the class.
///
/// The reference can write `stop_retained_variant` in the place of either term: for an
/// insertion or deletion without a peptide `stop_retained` holds where the overlap with the
/// stop codon holds, read with an insertion's start moved by the inserted length, and the
/// stop codon test returns before it reads a base (Utils/VariationEffect.pm:1318,
/// :1408-1415, :1498). On such a pair, an insertion's or a deletion's
/// (`length_changing_row`) with both rows showing that the test is not made
/// (`stop_codon_test_not_made`), the reference alone writes `stop_retained_variant` and
/// vep-rs alone `stop_gained`, with `frameshift_variant`, one of the three in-frame terms or
/// no other term.
///
/// For a kept stop the two sets differ as `kept_stop_terms` reads them;
/// `stop_codon_overlap_insertion_range_pair`, the earlier rule of the registry, takes the
/// kept stop in the place of `frameshift_variant` as well.
pub fn stop_codon_settled_before_unknown_bases_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e)
        || e.allele != a.allele
        || e.cds_position != a.cds_position
        || position_bounds(&e.cds_position) != (true, true)
        || e.codons != a.codons
        || e.amino_acids.as_deref() != Some("-")
    {
        return false;
    }
    let Some(stop) = e.codons.as_deref().and_then(settled_stop_codon) else {
        return false;
    };
    let first_stop_residue = a
        .amino_acids
        .as_deref()
        .and_then(|residues| residues.split_once('/'))
        .and_then(|(_, alternate)| alternate.chars().position(|residue| residue == '*'));
    if first_stop_residue != Some(stop) {
        return false;
    }
    let (expected, gone, come) = term_differences(e, a);
    if kept_stop_terms(e, a, &expected, &gone, &come) {
        return true;
    }
    if expected.contains(FRAMESHIFT_TERM) {
        return gone.is_empty() && is(&come, &[STOP_GAINED_TERM]);
    }
    if !come.contains(STOP_GAINED_TERM) || come.len() > 2 {
        return false;
    }
    let beside_the_stop = come.iter().find(|term| **term != STOP_GAINED_TERM);
    let inframe = beside_the_stop.is_none_or(|term| FRAMESHIFT_READ_TERMS.contains(term));
    if is(&gone, &[STOP_RETAINED_TERM]) {
        return length_changing_row(e)
            && stop_codon_test_not_made(e)
            && stop_codon_test_not_made(a)
            && (inframe || beside_the_stop == Some(&FRAMESHIFT_TERM));
    }
    is(&gone, &[CODING_SEQUENCE_TERM]) && inframe
}

/// The groups of terms Ensembl VEP writes where its `frameshift` returns 0 for a reference
/// peptide that begins with a stop (Utils/VariationEffect.pm:1554) and the length change is
/// no multiple of three: `coding_sequence_variant` alone, `inframe_deletion` alone, or
/// `inframe_insertion` with `coding_sequence_variant`.
const HELD_STOP_CODON_TERM_GROUPS: [&[&str]; 3] = [
    &[CODING_SEQUENCE_TERM],
    &[INFRAME_DELETION_TERM],
    &[INFRAME_INSERTION_TERM, CODING_SEQUENCE_TERM],
];

/// The bases in the reference half and in the alternate half of a Codons cell, the text
/// before its first `/` and the text after it: a half of `-` holds none, and any other half
/// one base a character. `None` for a cell without a `/`.
fn codon_half_bases(codons: &str) -> Option<(usize, usize)> {
    let (reference, alternate) = codons.split_once('/')?;
    let bases = |half: &str| match half {
        "-" => 0,
        _ => half.chars().count(),
    };
    Some((bases(reference), bases(alternate)))
}

/// `frameshift_inside_a_held_stop_codon`: an insertion or a deletion whose length change is
/// no multiple of three and whose first base lies inside the codon of a stop that the
/// transcript's peptide holds as a stop and reads on past (a polymorphic pseudogene's).
/// Ensembl VEP's `frameshift` returns 0 wherever the reference peptide of the allele's
/// codons begins with a stop (Utils/VariationEffect.pm:1554), which holds for the stop codon
/// no residue follows, so its row holds a group of `HELD_STOP_CODON_TERM_GROUPS` where
/// vep-rs's holds `frameshift_variant`.
///
/// The rule reads the terms, the Amino_acids and the Codons of the two rows and no other
/// column. vep-rs's terms are the reference's with `frameshift_variant` in the place of
/// exactly one of those groups, every other term standing on both rows, and the
/// reference's hold neither `stop_retained_variant` nor `stop_lost`. The two rows hold one
/// Amino_acids cell, whose half before its first `/` opens with `*` and whose half after it
/// holds `X`, and one Codons cell, whose two halves differ in bases (`codon_half_bases`) by
/// a number that is no multiple of three. A row read without its Amino_acids and Codons
/// meets none of this.
pub fn frameshift_inside_a_held_stop_codon_pair(e: &PairRow, a: &PairRow) -> bool {
    if e.amino_acids != a.amino_acids || e.codons != a.codons {
        return false;
    }
    let residues = e
        .amino_acids
        .as_deref()
        .and_then(|cell| cell.split_once('/'));
    let Some((reference, alternate)) = residues else {
        return false;
    };
    if !reference.starts_with('*') || !alternate.contains('X') {
        return false;
    }
    let bases = e.codons.as_deref().and_then(codon_half_bases);
    let Some((reference_bases, alternate_bases)) = bases else {
        return false;
    };
    if reference_bases.abs_diff(alternate_bases) % 3 == 0 {
        return false;
    }
    let (expected, gone, come) = term_differences(e, a);
    is(&come, &[FRAMESHIFT_TERM])
        && !expected.contains(STOP_RETAINED_TERM)
        && !expected.contains(STOP_LOST_TERM)
        && HELD_STOP_CODON_TERM_GROUPS
            .iter()
            .any(|group| is(&gone, group))
}

const SYNONYMOUS_TERM: &str = "synonymous_variant";
const MISSENSE_TERM: &str = "missense_variant";
/// The one codon NCBI's genetic code 26 reads unlike the standard code: alanine, where the
/// standard code gives leucine.
const CODE_26_ALANINE_CODON: &str = "CTG";

/// The reference and the alternate codon of a Codons cell, in upper case: two codons of
/// three letters around one `/`.
fn codons_of(cell: &str) -> Option<(String, String)> {
    let (reference, alternate) = cell.split_once('/')?;
    let codon = |text: &str| {
        let letters = text.len() == 3 && text.bytes().all(|b| b.is_ascii_alphabetic());
        letters.then(|| text.to_ascii_uppercase())
    };
    Some((codon(reference)?, codon(alternate)?))
}

/// `pachysolen_nuclear_code_read_as_standard`: a Transcript row on which the reference
/// writes `synonymous_variant` and vep-rs `missense_variant`, the two sets otherwise equal,
/// both rows holding one Codons cell in which exactly one codon is the alanine codon of
/// genetic code 26. The reference's Amino_acids is leucine alone, the residue the last table
/// of its codon module gives both codons (Bio/Tools/CodonTable.pm:328-338, 412-460), and
/// vep-rs's is alanine on that codon's side of a `/` and leucine on the other. Of the genetic
/// codes vep-rs holds, code 26 alone reads the codon as alanine, so vep-rs writes such a row
/// on a transcript of that code and no other. A row read without its Amino_acids and Codons
/// meets none of this.
pub fn pachysolen_nuclear_code_read_as_standard_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) || e.codons != a.codons {
        return false;
    }
    let (_, gone, come) = term_differences(e, a);
    if !is(&gone, &[SYNONYMOUS_TERM]) || !is(&come, &[MISSENSE_TERM]) {
        return false;
    }
    let Some((reference, alternate)) = e.codons.as_deref().and_then(codons_of) else {
        return false;
    };
    let residues = match (
        reference == CODE_26_ALANINE_CODON,
        alternate == CODE_26_ALANINE_CODON,
    ) {
        (true, false) => "A/L",
        (false, true) => "L/A",
        _ => return false,
    };
    e.amino_acids.as_deref() == Some("L") && a.amino_acids.as_deref() == Some(residues)
}

const START_LOST_TERM: &str = "start_lost";
const START_RETAINED_TERM: &str = "start_retained_variant";
const INFRAME_DELETION_TERM: &str = "inframe_deletion";

/// The end of a position column that opens at the first base: the cell, read up to a
/// `/total` suffix, is `1`, a `-` and a number written in digits, and the number is the end.
/// `None` for any other cell, one position alone included.
fn end_of_a_span_from_the_first_base(position: &str) -> Option<&str> {
    let head = position.split_once('/').map_or(position, |(head, _)| head);
    let (start, end) = head.split_once('-')?;
    let number = !end.is_empty() && end.bytes().all(|b| b.is_ascii_digit());
    (start == "1" && number).then_some(end)
}

/// The ends of a row's cDNA_position and of its CDS_position where each column opens at the
/// first base (`end_of_a_span_from_the_first_base`), and `None` where either does not.
fn ends_of_the_spans_from_the_first_base(row: &PairRow) -> Option<(&str, &str)> {
    Some((
        end_of_a_span_from_the_first_base(&row.cdna_position)?,
        end_of_a_span_from_the_first_base(&row.cds_position)?,
    ))
}

/// Whether the alternate half of an Amino_acids cell, the text after its first `/`, is `-`:
/// the row writes no alternate residue.
fn no_alternate_residue(amino_acids: &str) -> bool {
    amino_acids
        .split_once('/')
        .is_some_and(|(_, alternate)| alternate == "-")
}

/// `start_codon_removed_without_five_prime_utr`: a Transcript row whose cDNA_position and
/// whose CDS_position each open at 1 and end on one number on both rows
/// (`end_of_a_span_from_the_first_base`), so the coding sequence opens on the transcript's
/// first base and the removal starts there. The two sets of terms are equal but for vep-rs's
/// `start_lost`, both hold `inframe_deletion`, the reference's holds no
/// `start_retained_variant`, and neither row's Amino_acids writes an alternate residue
/// (`no_alternate_residue`). `start_lost` (Utils/VariationEffect.pm:851-904) sets its length
/// test aside for an in-frame deletion (:870), `_inv_start_altered` (:906-945) returns before
/// its edit on a transcript with no 5' UTR (:926), and the peptide test returns on an empty
/// alternate peptide (:875), so the reference writes `inframe_deletion` without `start_lost`
/// for a deletion of whole codons from the first coding base of such a transcript; vep-rs
/// reads the edited coding sequence and writes `start_lost` beside it. A row read without its
/// Amino_acids is no row of the class.
pub fn start_codon_removed_without_five_prime_utr_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) {
        return false;
    }
    let (Some(expected_ends), Some(actual_ends)) = (
        ends_of_the_spans_from_the_first_base(e),
        ends_of_the_spans_from_the_first_base(a),
    ) else {
        return false;
    };
    if expected_ends != actual_ends
        || !e.amino_acids.as_deref().is_some_and(no_alternate_residue)
        || !a.amino_acids.as_deref().is_some_and(no_alternate_residue)
    {
        return false;
    }
    let (expected, gone, come) = term_differences(e, a);
    gone.is_empty()
        && is(&come, &[START_LOST_TERM])
        && expected.contains(INFRAME_DELETION_TERM)
        && !expected.contains(START_RETAINED_TERM)
}

/// The three terms by one of which the class reads a row as an insertion's or a deletion's
/// inside the coding sequence.
const LENGTH_CHANGE_TERMS: [&str; 3] =
    [FRAMESHIFT_TERM, INFRAME_INSERTION_TERM, "inframe_deletion"];

/// The first and the last number of the span a position column writes, read up to a
/// `/total` suffix and parted at its first `-`, one number standing for both. `None` where
/// a part is no number in the digits `0` to `9` without a leading zero: an empty cell, `-`,
/// `<absent>` and an end that maps to no transcript base (`?`) among them.
fn numbered_span(position: &str) -> Option<(&str, &str)> {
    let span = position.split_once('/').map_or(position, |(span, _)| span);
    let (first, last) = span.split_once('-').unwrap_or((span, span));
    let number = |part: &str| {
        !part.is_empty() && !part.starts_with('0') && part.bytes().all(|b| b.is_ascii_digit())
    };
    (number(first) && number(last)).then_some((first, last))
}

/// `start_codon_kept_without_five_prime_utr`: a Transcript row whose cDNA_position and
/// CDS_position write one span of numbers on both rows (`numbered_span`), so the coding
/// sequence opens on the transcript's first base, and whose span opens on a base of the
/// start codon: its first number is 1, 2 or 3. vep-rs's terms are the reference's with
/// `start_retained_variant` added and without the reference's `start_lost` where it holds
/// one: the reference's hold no `start_retained_variant`, vep-rs's no `start_lost`, every
/// other term is equal, and one of `frameshift_variant`, `inframe_insertion` and
/// `inframe_deletion` is among them. `_ins_del_start_altered`
/// (Utils/VariationEffect.pm:1028-1075) reads the three bases at the start position only
/// inside its test for a 5' UTR (:1061-1065). On a transcript without one it reads every
/// deletion as an altered start (:1069), and an insertion as one where the tail of the
/// edited sequence is not the coding sequence (:1071). For such an allele that leaves a
/// start codon there `start_retained_variant` (:947-963) is not written, and `start_lost`
/// (:851-904) holds where the length or the peptide test does (:870, :878-882); vep-rs
/// makes the kept-start test at the transcript's first base and writes
/// `start_retained_variant` alone.
pub fn start_codon_kept_without_five_prime_utr_pair(e: &PairRow, a: &PairRow) -> bool {
    if !transcript_row(e) {
        return false;
    }
    let Some(span) = numbered_span(&e.cds_position) else {
        return false;
    };
    let one_span = [&e.cdna_position, &a.cdna_position, &a.cds_position]
        .iter()
        .all(|position| numbered_span(position) == Some(span));
    if !one_span || !matches!(span.0, "1" | "2" | "3") {
        return false;
    }
    let (expected, gone, come) = term_differences(e, a);
    let start_lost_goes = if expected.contains(START_LOST_TERM) {
        is(&gone, &[START_LOST_TERM])
    } else {
        gone.is_empty()
    };
    start_lost_goes
        && is(&come, &[START_RETAINED_TERM])
        && LENGTH_CHANGE_TERMS
            .iter()
            .any(|term| expected.contains(term))
}

/// A set of the seven terms a row of `intron_rows` holds, one bit a term in the order of
/// [`INTRON_ROW_TERMS`].
type TermSet = u8;

/// The five splice terms whose windows lie inside an intron: the two bases at its 5' end, the
/// two at its 3' end, its fifth base, its third to sixth base, and the bases 3 to 17 from its
/// 3' end (Utils/VariationEffect.pm:558-614); the intron term, whose flag `within_intron`
/// (Utils/VariationEffect.pm:629-634) returns and whose window is the intron itself, its end
/// bases included; and `splice_region_variant`.
const INTRON_ROW_TERMS: [&str; 7] = [
    "splice_donor_variant",
    "splice_acceptor_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "intron_variant",
    SPLICE_REGION_TERM,
];
const DONOR: TermSet = 1;
const ACCEPTOR: TermSet = 1 << 1;
const DONOR_5TH_BASE: TermSet = 1 << 2;
const DONOR_REGION: TermSet = 1 << 3;
const TRACT: TermSet = 1 << 4;
const INTRON: TermSet = 1 << 5;
const REGION: TermSet = 1 << 6;
/// The six terms `_intron_effects` sets from the differing regions and vep-rs writes only
/// where the record's span holds a base of the term's window.
const PAST_THE_SPAN_TERMS: TermSet =
    DONOR | ACCEPTOR | DONOR_5TH_BASE | DONOR_REGION | TRACT | INTRON;
/// The four beside which `splice_region` returns 0 (Utils/VariationEffect.pm:616-627).
const SPLICE_REGION_WITHHOLDING_TERMS: TermSet = DONOR | ACCEPTOR | DONOR_5TH_BASE | DONOR_REGION;
/// A gap past which no window of an intron's end reaches a span or its regions.
const FAR: i64 = 1_000_000_000;
/// The names `--terms display` prints for the seven terms, each with the terms that print it
/// (the `display_term` of each term, Utils/Constants.pm): one name for the donor and the
/// acceptor term, one for the fifth base, the donor region, the tract and
/// `splice_region_variant`, one for the intron term. No other term prints one of the three.
const DISPLAY_ROW_TERMS: [(&str, TermSet); 3] = [
    ("ESSENTIAL_SPLICE_SITE", DONOR | ACCEPTOR),
    (
        "SPLICE_SITE",
        DONOR_5TH_BASE | DONOR_REGION | TRACT | REGION,
    ),
    ("INTRONIC", INTRON),
];
/// The row-aware rules of the classes of a coding end and of the stop codon classes, whose
/// terms a row can differ by beside the six: one row can carry both differences.
const CODING_END_ROW_PAIRS: [fn(&PairRow, &PairRow, &AtACodingEnd<'_>) -> bool; 10] = [
    structural_deletion_stop_lost_without_stop_codon_pair,
    stop_codon_completed_by_polyadenylation_pair,
    stop_codon_overlap_without_stop_codon_pair,
    stop_codon_behind_last_residue_without_stop_codon_pair,
    start_codon_overlap_without_start_codon_pair,
    stop_codon_overlap_genomic_span_pair_at,
    stop_retained_without_stop_codon_test_pair_at,
    stop_retained_three_base_stop_codon_test_pair_at,
    stop_codon_removed_without_stop_codon_test_pair_at,
    stop_codon_overlap_insertion_range_pair_at,
];

/// The terms of the seven among a row's terms.
fn intron_row_terms(terms: &BTreeSet<&str>) -> TermSet {
    INTRON_ROW_TERMS
        .iter()
        .enumerate()
        .filter(|(_, term)| terms.contains(*term))
        .fold(0, |set, (bit, _)| set | 1 << bit)
}

/// The display names of the seven among a row's terms, one bit a name in the order of
/// [`DISPLAY_ROW_TERMS`].
fn display_row_names(terms: &BTreeSet<&str>) -> u8 {
    DISPLAY_ROW_TERMS
        .iter()
        .enumerate()
        .filter(|(_, (name, _))| terms.contains(name))
        .fold(0, |set, (bit, _)| set | 1 << bit)
}

/// The display names a row of the seven prints, as [`display_row_names`] reads them.
fn displayed(row: TermSet) -> u8 {
    DISPLAY_ROW_TERMS
        .iter()
        .enumerate()
        .filter(|(_, (_, members))| row & members != 0)
        .fold(0, |set, (bit, _)| set | 1 << bit)
}

/// The terms of a row outside the seven.
fn other_terms<'a>(terms: &BTreeSet<&'a str>) -> BTreeSet<&'a str> {
    terms
        .iter()
        .copied()
        .filter(|term| !INTRON_ROW_TERMS.contains(term))
        .collect()
}

/// The start and the end a position column writes: an empty cell or `-` is `-` twice; any
/// other cell is read up to a `/total` suffix and parted at its first `-`, one number standing
/// for both ends. `None` where a part is not digits, `?` or `-`.
fn position_parts(position: &str) -> Option<(&str, &str)> {
    if matches!(position, "" | "-") {
        return Some(("-", "-"));
    }
    let head = position.split_once('/').map_or(position, |(head, _)| head);
    let (start, end) = head.split_once('-').unwrap_or((head, head));
    let written = |part: &str| {
        matches!(part, "?" | "-") || (!part.is_empty() && part.bytes().all(|b| b.is_ascii_digit()))
    };
    (written(start) && written(end)).then_some((start, end))
}

/// The cDNA start and end and the CDS start and end a row carries, as written (`?` for an end
/// that maps to no transcript base, `-` for a position the row lacks); `None` where a
/// position is not of that form.
fn transcript_positions(e: &PairRow) -> Option<[&str; 4]> {
    let (cdna_start, cdna_end) = position_parts(&e.cdna_position)?;
    let (cds_start, cds_end) = position_parts(&e.cds_position)?;
    Some([cdna_start, cdna_end, cds_start, cds_end])
}

/// The number of positions a Location names, a lone position counting as one, or
/// `None` where it names no span of whole positions.
fn positions_named(location: &str) -> Option<i128> {
    let (start, end) = location_ends(location)?;
    let span = (end - start).saturating_add(1);
    (span >= 1).then_some(span)
}

/// The number of bases of the transcript's sequence from the first cDNA position of a row to
/// its second, both counted, where both are numbers. A position past `i128` is read as no
/// number.
fn sequence_width(positions: &[&str; 4]) -> Option<i128> {
    let (start, end): (i128, i128) = (positions[0].parse().ok()?, positions[1].parse().ok()?);
    Some((end - start).saturating_add(1))
}

/// Whether a row shows, without its record, a span of three or more bases that lies in the
/// transcript's sequence from end to end under an allele of more nucleotide letters than the
/// span has bases: its Location names the span, and its two cDNA positions are numbers as far
/// apart as the span's two ends, which a span that holds an intron's base does not have. A
/// Location of two neighbouring positions is refused: an insertion between them is written
/// with the same Location, allele and positions as two replaced bases, and only the record
/// tells the two apart.
fn longer_allele_inside_one_exon(e: &PairRow) -> bool {
    let Some(span) = positions_named(&e.location) else {
        return false;
    };
    if span < 3 || !spelt_in(&e.allele, b"ACGTN") || e.allele.len() as i128 <= span {
        return false;
    }
    transcript_positions(e).is_some_and(|positions| sequence_width(&positions) == Some(span))
}

/// The regions `_get_differing_regions` (VariationFeatureOverlapAllele.pm:405-460) cuts a
/// record into, as offsets from the span's first base: where both alleles hold more than one
/// base (:418) the runs of indices at which they differ, every index past the shorter allele's
/// end differing (:428-450); else one region, the span (:452-454).
fn differing_regions(reference: &[char], allele: &[char]) -> Vec<(i64, i64)> {
    let span = vec![(0, reference.len() as i64 - 1)];
    if reference.len() <= 1 || allele.len() <= 1 {
        return span;
    }
    let mut runs: Vec<(i64, i64)> = Vec::new();
    for index in 0..reference.len().max(allele.len()) {
        if index < reference.len() && index < allele.len() && reference[index] == allele[index] {
            continue;
        }
        let index = index as i64;
        match runs.last_mut() {
            Some(run) if run.1 == index - 1 => run.1 = index,
            _ => runs.push((index, index)),
        }
    }
    if runs.is_empty() {
        span
    } else {
        runs
    }
}

/// The splice and intron terms of the reference's row and of vep-rs's for a record whose span
/// is the offsets `0` to `span - 1`, on a transcript of the given strand whose one intron near
/// the span runs from offset `start` to offset `end`.
///
/// The reference's row is `_intron_effects` (BaseTranscriptVariationAllele.pm:99-224) on
/// `regions`, each placed with no bound at the span's end (:117-118). The intron is read
/// where the span overlaps it widened by three bases, and its two ends where the span
/// overlaps the eleven bases around each (BaseTranscriptVariation.pm:1027-1029); an intron of
/// thirteen bases or fewer is a frameshift intron (:1032-1035) and a region that overlaps it
/// sets nothing (BaseTranscriptVariationAllele.pm:138-141, :171-174). A region in the
/// intron's interior sets the intron term (:143-150) and one in a window of the five splice
/// terms its flag (:152-161, :176-206). The flag of `splice_region` is assigned for each
/// region in turn, whether the region lies in one of its four windows, until a region has
/// met the intron's first two or last two bases (:209-216; Utils/VariationEffect.pm:91-94),
/// so the last region assigned decides it. The tract term is written only off the exons
/// (Constants.pm:971-983), which a transcript with a frameshift intron reads twelve bases
/// wider (BaseTranscriptVariation.pm:861; `wide_exons` says another intron of the transcript
/// is one). The donor region is withheld beside the fifth base and `splice_region_variant`
/// beside the first four terms (Utils/VariationEffect.pm:597-627).
///
/// vep-rs's row keeps a term of the six only where the span itself holds a base of a window
/// of that term, the intron being the intron term's window, reads `splice_region` from any
/// region in one of its four windows, and withholds the donor region and
/// `splice_region_variant` by its own row's terms. For an allele that opens with its whole
/// reference allele (`insertion_after_the_span`) it is the row of the insertion the allele
/// writes (`insertion_after_the_span_row`).
fn intron_rows(
    regions: &[(i64, i64)],
    span: i64,
    start: i64,
    end: i64,
    forward: bool,
    wide_exons: bool,
) -> (TermSet, TermSet) {
    let last = span - 1;
    let in_list = last >= start - 3 && end + 3 >= 0;
    let at_an_end = (last >= start - 3 && start + 7 >= 0) || (last >= end - 7 && end + 3 >= 0);
    let frameshift = end - start <= 12;
    let read = || {
        let outside = move |&(from, to): &(i64, i64)| !(frameshift && to >= start && from <= end);
        regions.iter().copied().filter(outside)
    };
    let met =
        |first: i64, final_base: i64| read().any(|(from, to)| to >= first && from <= final_base);
    let held = |first: i64, final_base: i64| last >= first && final_base >= 0;
    let windows = if forward {
        [
            (DONOR, start, start + 1),
            (ACCEPTOR, end - 1, end),
            (DONOR_5TH_BASE, start + 4, start + 4),
            (DONOR_REGION, start + 2, start + 5),
            (TRACT, end - 16, end - 2),
        ]
    } else {
        [
            (DONOR, end - 1, end),
            (ACCEPTOR, start, start + 1),
            (DONOR_5TH_BASE, end - 4, end - 4),
            (DONOR_REGION, end - 5, end - 2),
            (TRACT, start + 2, start + 16),
        ]
    };
    let off_the_exons = !(start > 0
        || last > end
        || ((frameshift || wide_exons) && (start > -12 || last + 12 > end)));
    let (mut set_by_a_region, mut held_by_the_span) = (0, 0);
    for (term, first, final_base) in windows {
        let read_here = if term == TRACT {
            (in_list || at_an_end) && off_the_exons
        } else {
            at_an_end
        };
        if read_here && met(first, final_base) {
            set_by_a_region |= term;
        }
        if held(first, final_base) {
            held_by_the_span |= term;
        }
    }
    if in_list && start + 2 <= end - 2 && met(start + 2, end - 2) {
        set_by_a_region |= INTRON;
    }
    if held(start, end) {
        held_by_the_span |= INTRON;
    }
    let beside = [
        (start - 3, start - 1),
        (start + 2, start + 7),
        (end - 7, end - 2),
        (end + 1, end + 3),
    ];
    let (mut last_assigned, mut in_any, mut at_a_site) = (false, false, false);
    if at_an_end {
        for (from, to) in read() {
            let in_one = beside
                .iter()
                .any(|&(first, final_base)| to >= first && from <= final_base);
            in_any |= in_one;
            at_a_site |= (to >= start && from <= start + 1) || (to >= end - 1 && from <= end);
            if !at_a_site {
                last_assigned = in_one;
            }
        }
    }
    let row = |mut terms: TermSet, region: bool| {
        if terms & DONOR_5TH_BASE != 0 {
            terms &= !DONOR_REGION;
        }
        if region && terms & SPLICE_REGION_WITHHOLDING_TERMS == 0 {
            terms |= REGION;
        }
        terms
    };
    let candidate = if insertion_after_the_span(regions, span) {
        insertion_after_the_span_row(span, start, end, forward, wide_exons)
    } else {
        row(set_by_a_region & held_by_the_span, in_any)
    };
    (row(set_by_a_region, last_assigned), candidate)
}

/// Whether `regions` are those of an allele that opens with its whole reference allele of
/// `span` bases, two or more, and holds further bases: one run that begins at the offset
/// `span`, no base of the span differing (VariationFeatureOverlapAllele.pm:418, :428-450).
/// Such an allele changes no base of its span, and vep-rs reads it as the insertion it
/// writes, between the span's last base and the base after it.
fn insertion_after_the_span(regions: &[(i64, i64)], span: i64) -> bool {
    span >= 2 && matches!(regions, [(first, _)] if *first == span)
}

/// Whether `allele` opens with its reference allele `reference` of one base and holds
/// further bases. The reference does not compare such a pair of alleles and reads the one
/// base as changed (VariationFeatureOverlapAllele.pm:418, :452-454), so the regions do not
/// show the shape. Such an allele changes no base of its span either, and vep-rs reads it
/// as the insertion it writes, between that base and the base after it.
fn insertion_after_one_base(reference: &[char], allele: &[char]) -> bool {
    matches!(reference, [base] if allele.len() > 1 && allele[0] == *base)
}

/// Whether the rows of a pair are those of an allele that opens with its whole reference
/// allele (`insertion_after_the_span`, `insertion_after_one_base`), by the reference bases
/// either row carries; false where neither carries them.
pub fn written_after_the_span(e: &PairRow, a: &PairRow) -> bool {
    let Some(reference) = e.reference.as_ref().or(a.reference.as_ref()) else {
        return false;
    };
    let reference: Vec<char> = reference.chars().collect();
    let allele: Vec<char> = e.allele.chars().collect();
    insertion_after_one_base(&reference, &allele)
        || insertion_after_the_span(
            &differing_regions(&reference, &allele),
            reference.len() as i64,
        )
}

/// Whether a tuple's Location and allele allow an allele that opens with its whole reference
/// allele: the allele is a sequence of more letters than the Location names positions. The
/// record's alleles decide the rest (`written_after_the_span`).
pub fn may_be_written_after_the_span(location: &str, allele: &str) -> bool {
    spelt_in(allele, b"ACGTN")
        && positions_named(location).is_some_and(|positions| allele.len() as i128 > positions)
}

/// The splice and intron terms of vep-rs's row for an allele `insertion_after_the_span`
/// takes, on the intron `intron_rows` reads: those of the insertion between the span's last
/// base, the offset `span - 1`, and the base after it, read by those two neighbours as
/// `_intron_effects` reads an insertion (BaseTranscriptVariationAllele.pm:126).
///
/// The intron is read where the two neighbours overlap it widened by three bases, and its two
/// ends where they overlap the eleven bases around each (BaseTranscriptVariation.pm:1027-1029);
/// a frameshift intron that holds both sets nothing (BaseTranscriptVariationAllele.pm:138-141,
/// :171-174). A window of the donor, acceptor or donor region term that holds both neighbours
/// sets its flag (:176-206), which the one base of the fifth base's window never does. The
/// tract's window sets its flag where it holds either neighbour (:152-161), and the term is
/// written off the exons the two neighbours overlap, which a transcript with a frameshift
/// intron reads twelve bases wider. The intron term stands where the intron's interior holds
/// both neighbours, or the second is its third base, or the first its third base from the end
/// (:143-150). `splice_region` stands where one of its four windows holds both neighbours, or
/// the second neighbour is the intron's first or third base, or the first its last base or its
/// third from the end (Utils/VariationEffect.pm:91-102), and is withheld beside the donor,
/// acceptor and donor region terms.
fn insertion_after_the_span_row(
    span: i64,
    start: i64,
    end: i64,
    forward: bool,
    wide_exons: bool,
) -> TermSet {
    let (before, after) = (span - 1, span);
    let in_list = after >= start - 3 && before <= end + 3;
    let at_an_end =
        (after >= start - 3 && before <= start + 7) || (after >= end - 7 && before <= end + 3);
    let frameshift = end - start <= 12;
    let holds_both = |first: i64, final_base: i64| before >= first && after <= final_base;
    if !(in_list || at_an_end) || (frameshift && holds_both(start, end)) {
        return 0;
    }
    let (windows, tract) = if forward {
        (
            [
                (DONOR, start, start + 1),
                (ACCEPTOR, end - 1, end),
                (DONOR_REGION, start + 2, start + 5),
            ],
            (end - 16, end - 2),
        )
    } else {
        (
            [
                (DONOR, end - 1, end),
                (ACCEPTOR, start, start + 1),
                (DONOR_REGION, end - 5, end - 2),
            ],
            (start + 2, start + 16),
        )
    };
    let mut terms = 0;
    for (term, first, final_base) in windows {
        if at_an_end && holds_both(first, final_base) {
            terms |= term;
        }
    }
    let on_an_exon = before < start
        || after > end
        || ((frameshift || wide_exons) && (before < start + 12 || after + 12 > end));
    if after >= tract.0 && before <= tract.1 && !on_an_exon {
        terms |= TRACT;
    }
    if in_list && (holds_both(start + 2, end - 2) || after == start + 2 || before == end - 2) {
        terms |= INTRON;
    }
    let beside = [
        (start + 2, start + 7),
        (end - 7, end - 2),
        (start - 3, start - 1),
        (end + 1, end + 3),
    ];
    let region = beside
        .iter()
        .any(|&(first, final_base)| holds_both(first, final_base))
        || after == start
        || before == end
        || after == start + 2
        || before == end - 2;
    if at_an_end && region && terms & SPLICE_REGION_WITHHOLDING_TERMS == 0 {
        terms |= REGION;
    }
    terms
}

/// Calls `visit` with the places `(forward, start, end)` of one intron beside a span of `span`
/// bases that a row's cDNA positions allow, as offsets from the span's first base, for regions
/// that run `reach` bases past the span's end, until it returns true, and returns whether it
/// did. A numbered cDNA position is a base of the transcript's sequence and any other a base
/// of the intron; the first position is the span's first base on the forward strand and its
/// last on the reverse. Two numbers as far apart as the span put the intron after the span,
/// one to seven bases on, as far as one of its ends is still read, or before it, where an
/// intron shorter than its windows has them run on past its last base; two numbers closer
/// than the span put it inside the span. A gap wider than any window, and than the regions'
/// reach, is placed once.
fn intron_placements(
    span: i64,
    reach: i64,
    positions: &[&str; 4],
    mut visit: impl FnMut(bool, i64, i64) -> bool,
) -> bool {
    let last = span - 1;
    let before: Vec<i64> = (0..18).map(|gap| -gap).chain([-FAR]).collect();
    let gaps: BTreeSet<i64> = (0..42).chain(42.max(reach - 20)..reach + 19).collect();
    let after: Vec<i64> = gaps.iter().map(|gap| last + gap).chain([FAR]).collect();
    let numbered = |part: &str| !part.is_empty() && part.bytes().all(|b| b.is_ascii_digit());
    let (numbered_first, numbered_last) = (numbered(positions[0]), numbered(positions[1]));
    for forward in [true, false] {
        let (first_in_sequence, last_in_sequence) = match forward {
            true => (numbered_first, numbered_last),
            false => (numbered_last, numbered_first),
        };
        if first_in_sequence && last_in_sequence {
            let Some(width) = sequence_width(positions) else {
                continue;
            };
            if width == i128::from(span) {
                for start in (1..=7).map(|gap| last + gap) {
                    let far_ends = after.iter().copied().filter(|end| *end >= start);
                    let near_ends = (0..42).map(|length| start + length);
                    if far_ends
                        .chain(near_ends)
                        .any(|end| visit(forward, start, end))
                    {
                        return true;
                    }
                }
                if (-16..0).any(|start| (start..0).any(|end| visit(forward, start, end))) {
                    return true;
                }
            } else if 0 < width && width < i128::from(span) {
                let length = span - width as i64;
                if (1..=last - length).any(|start| visit(forward, start, start + length - 1)) {
                    return true;
                }
            }
        } else if first_in_sequence {
            if (1..span).any(|start| after.iter().any(|end| visit(forward, start, *end))) {
                return true;
            }
        } else if last_in_sequence {
            if (0..last).any(|end| before.iter().any(|start| visit(forward, *start, end))) {
                return true;
            }
        } else if before
            .iter()
            .any(|start| after.iter().any(|end| visit(forward, *start, *end)))
        {
            return true;
        }
    }
    false
}

/// A place of the transcript's intron, `(forward, start, end, wide_exons)` as `intron_rows`
/// reads it, at which the record's regions give the reference row's splice and intron terms
/// `expected` and the span, or the insertion written after it, gives vep-rs's `actual`, or
/// `None` where the cDNA positions `positions` allow no such place. The exons are read at
/// their own width or, for a tract term a row lacks, twelve bases wider. `reference` is the
/// reference bases the allele stands in for; the places are those of `intron_placements`, and
/// the first that fits is returned. For an allele that opens with its one reference base
/// (`insertion_after_one_base`) the reference's row is that of the base read as changed, and
/// vep-rs's the row of the insertion after it (`insertion_after_the_span_row`).
fn record_span_placement(
    positions: &[&str; 4],
    reference: &[char],
    allele: &[char],
    expected: TermSet,
    actual: TermSet,
) -> Option<(bool, i64, i64, bool)> {
    let regions = differing_regions(reference, allele);
    let span = reference.len() as i64;
    let after_one_base = insertion_after_one_base(reference, allele);
    let rows = (expected, actual);
    let mut place = None;
    intron_placements(
        span,
        allele.len() as i64 - span,
        positions,
        |forward, start, end| {
            let read = |wide_exons: bool| {
                placed_rows(
                    &regions,
                    span,
                    after_one_base,
                    (forward, start, end, wide_exons),
                )
            };
            let narrow = read(false);
            let wide_exons = if narrow == rows {
                false
            } else if ((narrow.0 & TRACT != 0 && expected & TRACT == 0)
                || (narrow.1 & TRACT != 0 && actual & TRACT == 0))
                && read(true) == rows
            {
                true
            } else {
                return false;
            };
            place = Some((forward, start, end, wide_exons));
            true
        },
    );
    place
}

/// The reference's row and vep-rs's at one place `(forward, start, end, wide_exons)` of the
/// intron: `intron_rows`, with vep-rs's row that of the insertion after the base
/// (`insertion_after_the_span_row`) for an allele that opens with its one reference base.
fn placed_rows(
    regions: &[(i64, i64)],
    span: i64,
    after_one_base: bool,
    (forward, start, end, wide_exons): (bool, i64, i64, bool),
) -> (TermSet, TermSet) {
    let (reference_row, span_row) = intron_rows(regions, span, start, end, forward, wide_exons);
    if after_one_base {
        let insertion = insertion_after_the_span_row(span, start, end, forward, wide_exons);
        (reference_row, insertion)
    } else {
        (reference_row, span_row)
    }
}

/// Whether the cDNA positions `positions` allow a place of the transcript's intron at which
/// the reference's splice and intron terms print the display names `expected` and vep-rs's
/// the names `actual` (`displayed`). The exons are read twelve bases wider as well where a
/// row at their own width holds the tract term, whose name three other terms print.
fn displayed_span_placement(
    positions: &[&str; 4],
    reference: &[char],
    allele: &[char],
    expected: u8,
    actual: u8,
) -> bool {
    let regions = differing_regions(reference, allele);
    let span = reference.len() as i64;
    let after_one_base = insertion_after_one_base(reference, allele);
    let shown = |(reference_row, candidate_row): (TermSet, TermSet)| {
        (displayed(reference_row), displayed(candidate_row)) == (expected, actual)
    };
    intron_placements(
        span,
        allele.len() as i64 - span,
        positions,
        |forward, start, end| {
            let read = |wide_exons: bool| {
                placed_rows(
                    &regions,
                    span,
                    after_one_base,
                    (forward, start, end, wide_exons),
                )
            };
            let narrow = read(false);
            shown(narrow) || ((narrow.0 | narrow.1) & TRACT != 0 && shown(read(true)))
        },
    )
}

/// The class on two rows printed under `--terms display`, or `None` where neither row holds
/// a display name of the seven (`DISPLAY_ROW_TERMS`). The rows' other terms are equal and
/// their names of the seven differ; vep-rs's row alone does not hold the name of the donor
/// and the acceptor term, as it alone holds neither term under the Sequence Ontology names;
/// the record's alleles show a longer alternate allele that is no pure insertion, as under
/// those names, and a place of one intron gives both rows' names
/// (`displayed_span_placement`). Where the allele is not written after its span, the
/// reference's row holds a name vep-rs's lacks, as it holds a term of the six alone under
/// those names. Without the record's alleles no such row is taken: three names stand for
/// seven terms, and the rows alone do not tell them apart.
fn display_vocabulary_pair(
    e: &PairRow,
    a: &PairRow,
    positions: &[&str; 4],
    expected_terms: &BTreeSet<&str>,
    actual_terms: &BTreeSet<&str>,
) -> Option<bool> {
    let expected = display_row_names(expected_terms);
    let actual = display_row_names(actual_terms);
    if expected | actual == 0 {
        return None;
    }
    let others = |terms: &BTreeSet<&str>| -> Vec<String> {
        let named = |term: &&str| DISPLAY_ROW_TERMS.iter().any(|(name, _)| name == term);
        terms
            .iter()
            .filter(|term| !named(term))
            .map(|term| term.to_string())
            .collect()
    };
    let site_name = displayed(DONOR | ACCEPTOR);
    if expected == actual
        || actual & !expected & site_name != 0
        || others(expected_terms) != others(actual_terms)
    {
        return Some(false);
    }
    let Some(reference) = e.reference.as_ref().or(a.reference.as_ref()) else {
        return Some(false);
    };
    let reference: Vec<char> = reference.chars().collect();
    let allele: Vec<char> = e.allele.chars().collect();
    let after_one_base = insertion_after_one_base(&reference, &allele);
    if (reference.len() < 2 && !after_one_base)
        || positions_named(&e.location) != Some(reference.len() as i128)
        || !spelt_in(&e.allele, b"ACGTN")
        || allele.len() <= reference.len()
    {
        return Some(false);
    }
    let after_the_span = after_one_base
        || insertion_after_the_span(
            &differing_regions(&reference, &allele),
            reference.len() as i64,
        );
    if !after_the_span && expected & !actual == 0 {
        return Some(false);
    }
    Some(displayed_span_placement(
        positions, &reference, &allele, expected, actual,
    ))
}

/// `splice_term_past_the_record_span`: a Transcript row with the same cDNA and CDS positions
/// on both sides whose two sets of terms differ by terms of the six the reference alone holds
/// and, at most, the donor region and `splice_region_variant` vep-rs alone holds, beside the
/// terms of a class of a coding end or of a stop codon class where its rule takes the rest of
/// the difference (`CODING_END_ROW_PAIRS`).
///
/// `_get_differing_regions` (VariationFeatureOverlapAllele.pm:405-460) compares the two
/// alleles index by index to the longer allele's length (:428-429), so every index past the
/// reference allele's end differs, and `_intron_effects` places each differing region on the
/// genome from the record's start with no bound at the span's end
/// (BaseTranscriptVariationAllele.pm:117-118): on a record whose alternate allele is the
/// longer one, the flags of the five splice terms and of the intron term (:143-206) are set
/// from positions the record does not replace. vep-rs writes such a term only where a base of
/// the span itself lies in a window of that term, the intron being the intron term's window.
///
/// An allele that opens with its whole reference allele of two or more bases changes no base
/// of its span: every index that differs lies past the reference allele's end (:428-429), and
/// the reference reads the further bases as a change of the bases after the span
/// (BaseTranscriptVariationAllele.pm:117-118). vep-rs reads such an allele as the insertion it
/// writes, between the span's last base and the base after it, and writes that insertion's
/// terms (`insertion_after_the_span_row`): on its row vep-rs alone may hold the tract term and
/// the intron term as well, and the reference need hold no term alone.
///
/// An allele that opens with its reference allele of one base changes no base of its span
/// either. The reference does not compare the two alleles (:418) and reads the one base as
/// changed (:452-454), so its terms are those of that base, and vep-rs's those of the
/// insertion after it (`insertion_after_one_base`): on its row the reference alone may hold
/// `splice_region_variant` as well.
///
/// With the record's alleles (`PairRow::reference`) the two rows are those of a longer
/// alternate allele that is no pure insertion, both alleles holding more than one base or
/// the allele opening with its one reference base, at a
/// place of one intron the cDNA positions allow (`record_span_placement`): no input names the
/// transcript's introns, so the place is read off the reference's row, and where two places
/// give that row and differ in what the span holds, the pair is taken at the one that gives
/// vep-rs's row. Without them the span of three or more positions lies in the transcript's
/// sequence under a longer allele (`longer_allele_inside_one_exon`), vep-rs's row holds none
/// of the six, and `splice_region_variant` alone comes, where a term that withholds it went;
/// a Location of two neighbouring positions, which an insertion's row shares, and a span that
/// holds an intron's base then stay differences, and no row is read as that of an allele that
/// opens with its whole reference allele, which the record alone shows.
///
/// Two rows printed under `--terms display` are read by the names the seven terms print
/// there, with the record's alleles alone (`display_vocabulary_pair`).
#[allow(dead_code)]
pub fn splice_term_past_the_record_span_pair(e: &PairRow, a: &PairRow) -> bool {
    splice_term_past_the_record_span_pair_at(e, a, &AtACodingEnd::NO_CACHE)
}

/// `splice_term_past_the_record_span_pair` in a comparison that may have been given the
/// candidate's cache: a row that differs by the terms of a class of a coding end or of a
/// stop codon class as well is set aside where that class's rule takes the rest of the
/// difference on the reading of the transcript's cache record `at` holds.
pub fn splice_term_past_the_record_span_pair_at(
    e: &PairRow,
    a: &PairRow,
    at: &AtACodingEnd<'_>,
) -> bool {
    if !transcript_row(e) {
        return false;
    }
    let Some(positions) = transcript_positions(e) else {
        return false;
    };
    if transcript_positions(a) != Some(positions) {
        return false;
    }
    let (expected_terms, actual_terms) = (terms(&e.consequence), terms(&a.consequence));
    if let Some(taken) = display_vocabulary_pair(e, a, &positions, &expected_terms, &actual_terms) {
        return taken;
    }
    let expected = intron_row_terms(&expected_terms);
    let actual = intron_row_terms(&actual_terms);
    let other = other_terms(&expected_terms);
    if other != other_terms(&actual_terms) {
        // The rest of the difference is a class of a coding end's or a stop codon class's:
        // its rule reads the two rows with vep-rs's splice and intron terms on both, and
        // this one the two rows with the reference's other terms on both. Those differ by
        // terms of the seven alone, the only terms of either row read below.
        let held = |bit: &usize| actual & (1 << bit) != 0;
        let splice_and_intron = (0..INTRON_ROW_TERMS.len())
            .filter(held)
            .map(|bit| INTRON_ROW_TERMS[bit]);
        let beside = PairRow {
            consequence: other
                .iter()
                .copied()
                .chain(splice_and_intron)
                .collect::<Vec<_>>()
                .join(","),
            ..e.clone()
        };
        if !CODING_END_ROW_PAIRS.iter().any(|pair| pair(&beside, a, at)) {
            return false;
        }
    }
    let (gone, come) = (expected & !actual, actual & !expected);
    let reference = e.reference.as_ref().or(a.reference.as_ref());
    let reference: Option<Vec<char>> = reference.map(|reference| reference.chars().collect());
    let allele: Vec<char> = e.allele.chars().collect();
    let after_one_base = reference
        .as_ref()
        .is_some_and(|reference| insertion_after_one_base(reference, &allele));
    let (may_come, differ) = if written_after_the_span(e, a) {
        (REGION | DONOR_REGION | TRACT | INTRON, gone | come != 0)
    } else {
        (REGION | DONOR_REGION, gone != 0)
    };
    let may_go = if after_one_base {
        PAST_THE_SPAN_TERMS | REGION
    } else {
        PAST_THE_SPAN_TERMS
    };
    if !differ || gone & !may_go != 0 || come & !may_come != 0 {
        return false;
    }
    let Some(reference) = reference else {
        return longer_allele_inside_one_exon(e)
            && longer_allele_inside_one_exon(a)
            && actual & PAST_THE_SPAN_TERMS == 0
            && (come == 0 || (come == REGION && gone & SPLICE_REGION_WITHHOLDING_TERMS != 0));
    };
    if (reference.len() < 2 && !after_one_base)
        || positions_named(&e.location) != Some(reference.len() as i128)
    {
        return false;
    }
    if !spelt_in(&e.allele, b"ACGTN") || allele.len() <= reference.len() {
        return false;
    }
    record_span_placement(&positions, &reference, &allele, expected, actual).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_row_carries_its_three_cells_after_the_positions_only_where_it_is_given_them() {
        let nine = [
            "1993671",
            "19:13208046",
            "A",
            "ENSG1",
            "ENST1",
            "Transcript",
            "missense_variant",
            "70",
            "40",
        ];
        let read = pair_row(&nine);
        assert_eq!(
            (&read.protein_position, &read.amino_acids, &read.codons),
            (&None, &None, &None)
        );
        let with_cells = read.clone().with_cells(["14", "R/K", "aGg/aAg"]);
        assert_eq!(with_cells.protein_position.as_deref(), Some("14"));
        assert_eq!(with_cells.amino_acids.as_deref(), Some("R/K"));
        assert_eq!(with_cells.codons.as_deref(), Some("aGg/aAg"));
        assert_eq!(
            PairRow {
                protein_position: None,
                amino_acids: None,
                codons: None,
                ..with_cells
            },
            read
        );
    }

    #[test]
    fn a_generated_name_with_a_position_in_other_digits_spells_no_record() {
        assert_eq!(uploaded_alleles("21_100_A/G"), Some(("A", vec!["G"])));
        assert_eq!(uploaded_position("21_100_A/G"), Some(100));
        assert_eq!(uploaded_alleles("21_١٠٠_A/G"), None);
        assert_eq!(uploaded_position("21_١٠٠_A/G"), None);
    }

    #[test]
    fn an_identifier_spells_no_alleles_and_a_generated_name_spells_the_records() {
        assert_eq!(uploaded_alleles("1993671"), None);
        assert_eq!(uploaded_alleles("rs_12_x"), None);
        assert_eq!(uploaded_alleles("21_100 _A/T"), None);
        assert_eq!(
            uploaded_alleles("21_10542414_A/T/AG"),
            Some(("A", vec!["T", "AG"]))
        );
        assert_eq!(
            uploaded_alleles("chr_un_7_100_-/<INS>"),
            Some(("-", vec!["<INS>"]))
        );
        assert_eq!(
            uploaded_alleles("21_100_A/T_5_C/G"),
            Some(("C", vec!["G"])),
            "the last position and alleles the name can end in"
        );
    }

    #[test]
    fn the_records_own_alleles_give_the_alleles_shape() {
        assert_eq!(allele_is_indel("GACCTG", &["G"], "-"), Some(true));
        assert_eq!(allele_is_indel("A", &["T"], "T"), Some(false));
        assert_eq!(allele_is_indel("TCT", &["TCTC"], "C"), Some(true));
        assert_eq!(allele_is_indel("A", &["AT", "G"], "AT"), Some(true));
        assert_eq!(allele_is_indel("A", &["AT", "G"], "G"), Some(false));
        assert_eq!(allele_is_indel("AT", &["A", "AG"], "-"), Some(true));
        assert_eq!(allele_is_indel("AT", &["A", "G"], "-"), None);
        assert_eq!(allele_is_indel("A", &["T", "G"], "C"), None);
        assert_eq!(allele_is_indel("A", &["<DEL>"], "deletion"), None);
        assert_eq!(allele_is_indel("A", &["T"], ""), None);
        assert_eq!(allele_is_indel("AT", &["A", "*"], "-"), Some(true));
    }

    /// A row of one transcript under an identifier: the Location, the allele, the terms and
    /// the two position columns.
    fn row(location: &str, allele: &str, consequence: &str, cdna: &str, cds: &str) -> PairRow {
        pair_row(&[
            "1993671",
            location,
            allele,
            "ENSG1",
            "ENST00000635895",
            "Transcript",
            consequence,
            cdna,
            cds,
        ])
    }

    #[test]
    fn a_position_column_carries_its_start_and_its_end() {
        for (position, bounds) in [
            ("1705", (true, true)),
            ("1704-1705", (true, true)),
            ("?-7057", (false, true)),
            ("1705-?", (true, false)),
            ("-", (false, false)),
            ("", (false, false)),
            ("-/7140", (false, false)),
            ("120-7100/7140", (true, true)),
            ("?-6980/7023", (false, true)),
        ] {
            assert_eq!(position_bounds(position), bounds, "{position}");
        }
    }

    #[test]
    fn a_location_and_an_allele_tell_an_insertion_or_a_deletion() {
        assert_eq!(location_span("3:36993662-36993663"), Some(2));
        assert_eq!(location_span("3:36993663"), None);
        assert_eq!(location_span("3:5-x"), None);
        assert_eq!(location_span("36993662-36993663"), Some(2));
        assert!(indel_row(&row("3:36993663", "-", "", "", "")));
        assert!(indel_row(&row("3:36993662-36993663", "A", "", "", "")));
        assert!(!indel_row(&row("3:36993663", "A", "", "", "")));
        assert!(!indel_row(&row("3:36993662-36993663", "AC", "", "", "")));
        assert!(!indel_row(&row("3:36993662-36993663", "R", "", "", "")));
        let mut named = row("3:36993663", "A", "", "", "");
        named.indel = Some(true);
        assert!(
            indel_row(&named),
            "the record's own alleles come before the Location"
        );
    }

    #[test]
    fn the_terms_of_a_consequence_are_a_set() {
        assert!(is(&terms("b,a&b,,"), &["a", "b"]));
        assert!(is(&terms(""), &[""]));
        let (expected, gone, come) = (terms("a,b"), terms("a,b"), terms("a,c"));
        assert!(expected == gone && gone != come);
    }

    #[test]
    fn the_genomic_span_class_reads_a_tested_row_at_an_exon_end() {
        let terms = |term: &str| format!("splice_region_variant,{term},NMD_transcript_variant");
        let pair = |e: &str, a: &str, cdna: &str, cds: &str| {
            stop_codon_overlap_genomic_span_pair(
                &row("3:36993663", "-", &terms(e), cdna, cds),
                &row("3:36993663", "-", &terms(a), cdna, cds),
            )
        };
        let (coding, retained) = ("coding_sequence_variant", "stop_retained_variant");
        assert!(pair(coding, retained, "1705", "1630"));
        assert!(pair(coding, retained, "1705/1800", "1630/1653"));
        assert!(!pair(retained, coding, "1705", "1630"));
        assert!(!pair(coding, retained, "1705-?", "1630-?"));
        assert!(!pair(coding, retained, "1705", "-"));
        assert!(!pair(
            coding,
            "stop_retained_variant,stop_lost",
            "1705",
            "1630"
        ));
        assert!(!stop_codon_overlap_genomic_span_pair(
            &row("3:36993663", "-", coding, "1705", "1630"),
            &row("3:36993663", "-", retained, "1705", "1630"),
        ));
        assert!(!stop_codon_overlap_genomic_span_pair(
            &row("3:36993663", "A", &terms(coding), "1705", "1630"),
            &row("3:36993663", "A", &terms(retained), "1705", "1630"),
        ));
    }

    #[test]
    fn the_untested_class_reads_rows_whose_stop_codon_test_is_not_made() {
        let pair = stop_retained_without_stop_codon_test_pair;
        let acceptor = |term: &str| format!("splice_acceptor_variant,{term}");
        let (location, retained) = ("19:13208046-13208050", "stop_retained_variant");
        for (cdna, cds) in [
            ("?-7057", "-"),
            ("-", "-"),
            ("1705-?", "1630-?"),
            ("120-7100", "-"),
        ] {
            let reference = row(location, "-", &acceptor(retained), cdna, cds);
            for term in ["coding_sequence_variant", "frameshift_variant"] {
                assert!(pair(
                    &reference,
                    &row(location, "-", &acceptor(term), cdna, cds)
                ));
            }
            let nothing = row(location, "-", "splice_acceptor_variant", cdna, cds);
            assert!(pair(&reference, &nothing), "{cdna} {cds}");
        }
        let lettered = |term: &str| {
            row(
                "7:152648662-152648663",
                "ACRT",
                term,
                "7038-7039",
                "6961-6962",
            )
        };
        assert!(pair(&lettered(retained), &lettered("frameshift_variant")));
        let tested = |term: &str| row(location, "-", &acceptor(term), "7055-7059", "6978-6982");
        assert!(!pair(&tested(retained), &tested("coding_sequence_variant")));
        let untested = |term: &str| row(location, "-", &acceptor(term), "?-7057", "-");
        assert!(!pair(&untested(retained), &untested("stop_lost")));
        assert!(!pair(
            &untested(retained),
            &untested("coding_sequence_variant,frameshift_variant")
        ));
        assert!(!pair(
            &untested("coding_sequence_variant"),
            &untested(retained)
        ));
        assert!(!pair(
            &untested(retained),
            &tested("coding_sequence_variant")
        ));
        let substitution = |term: &str| row("19:13208046", "A", term, "-", "-");
        assert!(!pair(
            &substitution(retained),
            &substitution("coding_sequence_variant")
        ));
    }

    #[test]
    fn the_three_base_class_reads_rows_whose_stop_codon_test_is_made() {
        let pair = stop_retained_three_base_stop_codon_test_pair;
        let location = "2:190970701-190970707";
        let utr = |term: &str, cdna: &str| {
            row(
                location,
                "-",
                &format!("{term},3_prime_UTR_variant"),
                cdna,
                "2249-?",
            )
        };
        let (retained, coding) = ("stop_retained_variant", "coding_sequence_variant");
        assert!(pair(&utr(retained, "2558-2564"), &utr(coding, "2558-2564")));
        let coded = |term: &str| row(location, "-", term, "961-967", "264-270");
        assert!(pair(&coded(retained), &coded("frameshift_variant")));
        assert!(pair(
            &coded("inframe_deletion,stop_retained_variant"),
            &coded("frameshift_variant")
        ));
        assert!(pair(
            &coded("protein_altering_variant,stop_retained_variant"),
            &coded("stop_gained,frameshift_variant")
        ));
        assert!(!pair(&utr(retained, "2558-?"), &utr(coding, "2558-?")));
        assert!(!pair(&utr(retained, "2558-2564"), &utr(coding, "2558-?")));
        assert!(!pair(
            &utr(coding, "2558-2564"),
            &utr(retained, "2558-2564")
        ));
        assert!(!pair(
            &coded("inframe_deletion,stop_retained_variant"),
            &coded(coding)
        ));
        assert!(!pair(
            &coded("missense_variant,stop_retained_variant"),
            &coded("frameshift_variant")
        ));
        assert!(!pair(
            &utr(retained, "2558-2564"),
            &utr("stop_lost", "2558-2564")
        ));
        let lettered = |term: &str| {
            row(
                "7:152648662-152648663",
                "ACGTACGN",
                term,
                "1213-1214",
                "1100-1101",
            )
        };
        assert!(!pair(&lettered(retained), &lettered("frameshift_variant")));
    }

    #[test]
    fn the_removed_codon_class_reads_rows_without_the_tests_coordinates() {
        let pair = stop_codon_removed_without_stop_codon_test_pair;
        let location = "19:13208046-13208050";
        let acceptor = |term: &str, cdna: &str, cds: &str| {
            row(
                location,
                "-",
                &format!("splice_acceptor_variant,{term}"),
                cdna,
                cds,
            )
        };
        let (retained, lost) = ("stop_retained_variant", "stop_lost");
        for (cdna, cds) in [
            ("?-7057", "-"),
            ("-", "-"),
            ("1705-?", "1630-?"),
            ("120-7100", "-"),
            ("21-99", "?-69"),
        ] {
            assert!(
                pair(&acceptor(retained, cdna, cds), &acceptor(lost, cdna, cds)),
                "{cdna} {cds}"
            );
        }
        for (cdna, cds) in [("7055-7059", "6978-6982"), ("7055-7059", "6978-?")] {
            assert!(
                !pair(&acceptor(retained, cdna, cds), &acceptor(lost, cdna, cds)),
                "{cdna} {cds}"
            );
        }
        assert!(!pair(
            &acceptor(retained, "?-7057", "-"),
            &acceptor("stop_lost,frameshift_variant", "?-7057", "-")
        ));
        assert!(!pair(
            &acceptor(lost, "?-7057", "-"),
            &acceptor(retained, "?-7057", "-")
        ));
        assert!(!pair(
            &acceptor(retained, "?-7057", "-"),
            &acceptor("coding_sequence_variant", "?-7057", "-")
        ));
        let substitution = |term: &str| row("19:13208046", "A", term, "-", "-");
        assert!(!pair(&substitution(retained), &substitution(lost)));
    }

    #[test]
    fn the_insertion_range_class_reads_an_insertion_that_keeps_the_stop_codon() {
        let pair = stop_codon_overlap_insertion_range_pair;
        let inserted = |allele: &str, terms: &str, cdna: &str, cds: &str| {
            row("19:13208046-13208047", allele, terms, cdna, cds)
        };
        let whole = |terms: &str| inserted("TTGA", terms, "7054-7055", "6977-6978");
        for (reference, candidate) in [
            (
                "inframe_insertion,coding_sequence_variant",
                "inframe_insertion,stop_retained_variant",
            ),
            (
                "stop_gained,frameshift_variant",
                "inframe_insertion,stop_retained_variant",
            ),
            (
                "frameshift_variant",
                "inframe_insertion,stop_retained_variant",
            ),
            ("stop_gained,frameshift_variant", "stop_retained_variant"),
        ] {
            assert!(pair(&whole(reference), &whole(candidate)), "{reference}");
            assert!(!pair(&whole(candidate), &whole(reference)), "{reference}");
            assert!(pair(
                &inserted("TNGA", reference, "7054-7055", "6978-?"),
                &inserted("TNGA", candidate, "7054-7055", "6978-?")
            ));
            assert!(!pair(
                &inserted("TTGA", reference, "7054-?", "6977-?"),
                &inserted("TTGA", candidate, "7054-?", "6977-?")
            ));
            assert!(!pair(
                &inserted("TRGA", reference, "7054-7055", "6977-6978"),
                &inserted("TRGA", candidate, "7054-7055", "6977-6978")
            ));
            let deleted =
                |terms: &str| row("19:13208046-13208049", "-", terms, "7055-7058", "6978-6981");
            assert!(!pair(&deleted(reference), &deleted(candidate)));
        }
        let region = |term: &str| whole(&format!("inframe_insertion,splice_region_variant,{term}"));
        assert!(!pair(
            &region("coding_sequence_variant"),
            &region("stop_retained_variant")
        ));
        assert!(stop_codon_overlap_genomic_span_pair(
            &region("coding_sequence_variant"),
            &region("stop_retained_variant")
        ));
        assert!(!pair(
            &whole("frameshift_variant"),
            &whole("protein_altering_variant,stop_retained_variant")
        ));
        assert!(!pair(
            &whole("stop_gained"),
            &whole("stop_retained_variant")
        ));
    }

    #[test]
    fn a_reference_allele_is_empty_once_the_shared_first_base_is_dropped() {
        for (reference, alternates) in [
            ("-", vec!["TGA"]),
            ("A", vec!["ATGA"]),
            ("-", vec!["CG", "T"]),
            ("A", vec!["AT", "AG"]),
            ("A", vec!["AT", "*"]),
        ] {
            assert!(
                reference_allele_is_empty(reference, &alternates),
                "{reference}"
            );
        }
        for (reference, alternates) in [
            ("AC", vec!["GTT"]),
            ("AC", vec!["AGT"]),
            ("AC", vec!["A", "T"]),
            ("A", vec!["AC", "T"]),
            ("A", vec!["T"]),
            ("TCT", vec!["TCTC"]),
            ("TA", vec!["TATGA"]),
            ("AC", vec!["ATGAC"]),
            ("GACCTG", vec!["G"]),
        ] {
            assert!(
                !reference_allele_is_empty(reference, &alternates),
                "{reference}"
            );
        }
        let row = |name: &str, allele: &str| {
            pair_row(&[
                name,
                "19:13208046-13208047",
                allele,
                "ENSG1",
                "ENST00000635895",
                "Transcript",
                "stop_retained_variant",
                "7054-7055",
                "6977-6978",
            ])
        };
        for (name, shape) in [
            ("19_13208045_TA/TATGA", "after the record's second base"),
            ("19_13208046_AC/ATGAC", "between the record's two bases"),
        ] {
            let padded = row(name, "TGA");
            assert_eq!(
                (padded.insertion, padded.reference.as_deref()),
                (Some(false), Some("")),
                "an insertion {shape} keeps a base of its reference allele and stands in for none"
            );
            assert!(
                inserted_between_neighbours(&padded),
                "an insertion {shape} is read by the bases it stands in for"
            );
        }
        let identified = row("1993671", "CG");
        assert_eq!((identified.insertion, identified.reference), (None, None));
    }

    #[test]
    fn the_insertion_range_class_takes_an_insertion_between_two_positions_alone() {
        // `read` is the reference bases the input records give a row that carries none.
        let taken_beside = |name: &str, location: &str, allele: &str, read: Option<&str>| {
            let row = |terms: &str| {
                let mut row = pair_row(&[
                    name,
                    location,
                    allele,
                    "ENSG1",
                    "ENST00000635895",
                    "Transcript",
                    terms,
                    "7054-7055",
                    "6977-6978",
                ]);
                if row.reference.is_none() {
                    row.reference = read.map(str::to_string);
                }
                row
            };
            stop_codon_overlap_insertion_range_pair(
                &row("coding_sequence_variant"),
                &row("stop_retained_variant"),
            )
        };
        let taken =
            |name: &str, location: &str, allele: &str| taken_beside(name, location, allele, None);
        let two = "19:13208046-13208047";
        for (name, location, allele, shape) in [
            (
                "19_13208046_AC/GTT",
                two,
                "GTT",
                "two bases replaced by three",
            ),
            (
                "19_13208046_AC/AGT",
                two,
                "AGT",
                "the same, the first base shared",
            ),
            (
                "19_13208046_AC/A/T",
                two,
                "A",
                "a deletion written as a sequence",
            ),
            (
                "19_13208046_A/AC/T",
                "19:13208046",
                "AC",
                "an insertion kept with its first base",
            ),
            (
                "19_13208046_-/TGA",
                "19:13208046",
                "TGA",
                "a spelt insertion at one position",
            ),
            (
                "19_13208046_AC/A/T",
                two,
                "GTT",
                "an allele no alternate names, the reference of two",
            ),
            (
                "1993671",
                two,
                "CG",
                "a two-base allele under an identifier",
            ),
            ("1993671", "19:13208046-13208048", "TGA", "three positions"),
        ] {
            assert!(!taken(name, location, allele), "{shape}");
        }
        for (name, allele, shape) in [
            (
                "19_13208046_-/TGA",
                "TGA",
                "an insertion, the reference allele empty",
            ),
            (
                "19_13208046_A/ATGA",
                "TGA",
                "an insertion spelt with its shared first base",
            ),
            (
                "19_13208046_-/CG/T",
                "CG",
                "a two-base insertion, the alleles spelt",
            ),
            ("1993671", "TGA", "an insertion under an identifier"),
            (
                "1993671",
                "G",
                "under an identifier, two bases replaced by one",
            ),
            (
                "1993671",
                "CGT",
                "under an identifier, two bases replaced by three",
            ),
            (
                "19_13208045_TA/TATGA",
                "TGA",
                "an insertion after the second base its record keeps",
            ),
            (
                "19_13208046_AC/ATGAC",
                "TGA",
                "an insertion between the two bases its record keeps",
            ),
            (
                "19_13208046_A/ATGA/C",
                "TGA",
                "an insertion beside an alternate of another first base",
            ),
        ] {
            assert!(taken(name, two, allele), "{shape}");
        }
        for (name, allele, read, shape) in [
            ("1993671", "CG", "", "a two-base insertion"),
            ("1993671", "TGA", "", "an insertion"),
            (
                "19_13208046_-/TGA",
                "TGA",
                "AC",
                "replaced bases beside a first column that spells an insertion",
            ),
        ] {
            assert!(
                taken_beside(name, two, allele, Some(read)),
                "the input records give {shape}"
            );
        }
        for (allele, shape) in [
            ("CGT", "two bases replaced by three"),
            ("G", "two bases replaced by one"),
            ("CG", "two bases replaced by two"),
        ] {
            assert!(
                !taken_beside("1993671", two, allele, Some("AC")),
                "under an identifier, the input records give {shape}"
            );
        }
        assert!(
            !taken_beside("1993671", "19:13208046-13208048", "TGA", Some("")),
            "three positions, whatever the bases read"
        );
    }

    #[test]
    fn the_three_base_class_takes_the_retained_term_alone_or_beside_one_read_away_term() {
        let pair = stop_retained_three_base_stop_codon_test_pair;
        let coded = |terms: &str| row("2:190970701-190970707", "-", terms, "961-967", "264-270");
        let retained = "stop_retained_variant";
        for candidate in ["frameshift_variant", "stop_gained,frameshift_variant"] {
            assert!(pair(&coded(retained), &coded(candidate)), "alone");
            for one in FRAMESHIFT_READ_TERMS {
                let beside = format!("{one},{retained}");
                assert!(pair(&coded(&beside), &coded(candidate)), "{one}");
            }
            for several in [
                "inframe_insertion,inframe_deletion",
                "inframe_insertion,protein_altering_variant",
                "inframe_deletion,protein_altering_variant",
                "inframe_insertion,inframe_deletion,protein_altering_variant",
            ] {
                let beside = format!("{several},{retained}");
                assert!(!pair(&coded(&beside), &coded(candidate)), "{several}");
            }
        }
    }

    #[test]
    fn a_row_carries_its_columns_and_the_shape_of_its_allele() {
        let named = pair_row(&[
            "19_13208045_GACCTG/G",
            "19:13208046-13208050",
            "-",
            "ENSG1",
            "ENST00000635895",
            "Transcript",
            "stop_lost",
            "?-7057",
            "-",
        ]);
        assert_eq!(
            named,
            PairRow {
                name: "19_13208045_GACCTG/G".to_string(),
                location: "19:13208046-13208050".to_string(),
                allele: "-".to_string(),
                feature_type: "Transcript".to_string(),
                consequence: "stop_lost".to_string(),
                cdna_position: "?-7057".to_string(),
                cds_position: "-".to_string(),
                indel: Some(true),
                insertion: Some(false),
                reference: Some("ACCTG".to_string()),
                protein_position: None,
                amino_acids: None,
                codons: None,
            }
        );
        let identified = pair_row(&[
            "1993671",
            "19:13208046",
            "A",
            "ENSG1",
            "ENST00000635895",
            "Transcript",
            "stop_lost",
            "-",
            "-",
        ]);
        assert_eq!(
            (identified.indel, identified.insertion, identified.reference),
            (None, None, None)
        );
    }

    #[test]
    fn a_position_column_writes_its_start_and_its_end() {
        for (position, parts) in [
            ("999-1000", Some(("999", "1000"))),
            ("999-1000/2000", Some(("999", "1000"))),
            ("1705", Some(("1705", "1705"))),
            ("1705/1800", Some(("1705", "1705"))),
            ("1000-?", Some(("1000", "?"))),
            ("?-1000", Some(("?", "1000"))),
            ("?", Some(("?", "?"))),
            ("-", Some(("-", "-"))),
            ("", Some(("-", "-"))),
            ("-/7140", None),
            ("999-", None),
            ("-1000", None),
            ("999-1000-1001", None),
            ("9x9-1000", None),
            ("<absent>", None),
        ] {
            assert_eq!(position_parts(position), parts, "{position}");
        }
        assert_eq!(
            transcript_positions(&row("7:100998-101000", "GTGTG", "", "998-1000/2000", "-")),
            Some(["998", "1000", "-", "-"])
        );
        assert_eq!(
            transcript_positions(&row("7:100998-101000", "GTGTG", "", "998-1000", "9x9")),
            None
        );
    }

    #[test]
    fn a_generated_name_spells_the_records_position_before_its_alleles() {
        assert_eq!(uploaded_position("7_100999_TA/GTGTG"), Some(100999));
        assert_eq!(uploaded_position("chr_un_7_100_-/<INS>"), Some(100));
        assert_eq!(
            uploaded_position("21_100_A/T_5_C/G"),
            Some(5),
            "the last position and alleles the name can end in"
        );
        assert_eq!(uploaded_position("1993671"), None);
        assert_eq!(
            uploaded_position(&format!("7_{}_TA/GTGTG", "9".repeat(40))),
            None,
            "a position past 128 bits"
        );
    }

    #[test]
    fn a_location_names_its_positions() {
        assert_eq!(positions_named("7:100998-101000"), Some(3));
        assert_eq!(positions_named("7:101000"), Some(1));
        assert_eq!(positions_named("100999-101000"), Some(2));
        assert_eq!(positions_named("7:101000-100998"), None);
        assert_eq!(positions_named("7:100998-x"), None);
    }

    const TWO_POSITIONS: &str = "7:100999-101000";

    fn bases(sequence: &str) -> Option<String> {
        Some(sequence.to_string())
    }

    #[test]
    fn a_replacements_allele_stands_in_for_the_bases_of_its_location() {
        let replaced = replaced_reference;
        assert_eq!(
            replaced(100999, "TA", &["GTGTG"], TWO_POSITIONS, "GTGTG"),
            bases("TA")
        );
        assert_eq!(
            replaced(100998, "CTA", &["CGTGTG"], TWO_POSITIONS, "GTGTG"),
            bases("TA"),
            "behind a shared first base"
        );
        assert_eq!(
            replaced(100999, "ta", &["gtgtg"], TWO_POSITIONS, "GTGTG"),
            bases("TA"),
            "the record's alleles in lower case"
        );
        assert_eq!(
            replaced(100999, "TA", &["C", "GTGTG"], TWO_POSITIONS, "GTGTG"),
            bases("TA"),
            "beside another alternate"
        );
        assert_eq!(
            replaced(101000, "T", &["GTG"], "7:101000", "GTG"),
            bases("T"),
            "a Location of one position"
        );
    }

    #[test]
    fn a_pure_insertions_allele_stands_in_for_no_base() {
        let replaced = replaced_reference;
        assert_eq!(
            replaced(100999, "T", &["TGTG"], TWO_POSITIONS, "GTG"),
            bases("")
        );
        assert_eq!(
            replaced(100998, "CT", &["CTGTG"], TWO_POSITIONS, "GTG"),
            bases(""),
            "after the record's last base"
        );
        assert_eq!(
            replaced(100999, "TAC", &["TGTGAC"], "7:100999-101001", "GTG"),
            None,
            "a Location of three positions"
        );
    }

    #[test]
    fn a_deletions_allele_stands_in_for_the_bases_it_removes() {
        let replaced = replaced_reference;
        assert_eq!(
            replaced(100997, "CTAC", &["C"], "7:100998-101000", "-"),
            bases("TAC")
        );
        assert_eq!(
            replaced(100998, "CTA", &["C"], TWO_POSITIONS, "-"),
            bases("TA")
        );
        assert_eq!(
            replaced(100999, "TA", &["TA"], TWO_POSITIONS, "-"),
            None,
            "nothing inserted between two neighbouring positions"
        );
    }

    #[test]
    fn a_record_that_writes_no_such_row_gives_no_bases() {
        let replaced = replaced_reference;
        assert_eq!(
            replaced(100999, "TA", &["GTGTG", "C"], TWO_POSITIONS, "GTG"),
            None,
            "no alternate writes the allele"
        );
        assert_eq!(
            replaced(100999, "TA", &["GTG", "TGTGA"], TWO_POSITIONS, "GTG"),
            None,
            "one alternate replaces the two bases and one inserts between them"
        );
        assert_eq!(
            replaced(101000, "TA", &["GTG"], TWO_POSITIONS, "GTG"),
            None,
            "a Location that starts before the record"
        );
        assert_eq!(
            replaced(100999, "TA", &["GTR"], TWO_POSITIONS, "GTR"),
            None,
            "an allele that is not a sequence"
        );
        assert_eq!(replaced(100999, "TA", &["GTG"], "7:100999-x", "GTG"), None);
    }

    fn letters(sequence: &str) -> Vec<char> {
        sequence.chars().collect()
    }

    #[test]
    fn a_record_is_cut_into_the_runs_of_indices_at_which_its_alleles_differ() {
        let regions = |reference: &str, allele: &str| {
            differing_regions(&letters(reference), &letters(allele))
        };
        assert_eq!(regions("TA", "GAGT"), [(0, 0), (2, 3)]);
        assert_eq!(regions("TA", "GTGTG"), [(0, 4)]);
        assert_eq!(regions("TA", "TA"), [(0, 1)], "no index differs");
        assert_eq!(regions("T", "GTG"), [(0, 0)], "a reference of one base");
        assert_eq!(regions("TAC", "G"), [(0, 2)], "an allele of one base");
    }

    /// Whether the class takes the reference's and vep-rs's row of one transcript, the
    /// Location, the allele and the two position columns the same on both: the reference
    /// writes the donor term and vep-rs `splice_region_variant` in its place.
    fn past_the_span_rows(location: &str, allele: &str, cdna: &str, cds: &str) -> bool {
        let row = |terms: &str| row(location, allele, terms, cdna, cds);
        splice_term_past_the_record_span_pair(
            &row("splice_donor_variant,frameshift_variant"),
            &row("frameshift_variant,splice_region_variant"),
        )
    }

    /// Whether the class takes two rows of three positions replaced by five letters in the
    /// transcript's sequence, the reference's terms and vep-rs's as given.
    fn past_the_span_terms(reference: &str, candidate: &str) -> bool {
        let row = |terms: &str| row("7:100998-101000", "GTGTG", terms, "998-1000", "898-900");
        splice_term_past_the_record_span_pair(&row(reference), &row(candidate))
    }

    #[test]
    fn the_past_the_span_class_takes_a_donor_term_gone_and_the_region_term_come() {
        assert!(past_the_span_rows(
            "7:100998-101000",
            "GTGTG",
            "998-1000",
            "898-900"
        ));
        assert!(past_the_span_rows(
            "7:100998-101000",
            "GTGTG",
            "998-1000/2000",
            "898-900/1800"
        ));
    }

    #[test]
    fn the_past_the_span_class_takes_terms_of_the_six_gone_and_none_come() {
        assert!(past_the_span_terms(
            "splice_donor_variant,frameshift_variant,splice_donor_region_variant",
            "frameshift_variant"
        ));
        assert!(past_the_span_terms(
            "splice_polypyrimidine_tract_variant,frameshift_variant",
            "frameshift_variant"
        ));
        assert!(past_the_span_terms(
            "splice_donor_variant,frameshift_variant,intron_variant",
            "frameshift_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_span_across_an_exons_end() {
        assert!(!past_the_span_rows(
            "7:100998-101000",
            "GTGTG",
            "999-?",
            "899-?"
        ));
        assert!(!past_the_span_rows(
            "7:100998-101000",
            "GTGTG",
            "?-1000",
            "?-900"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_span_inside_an_intron() {
        assert!(!past_the_span_rows("7:100998-101000", "GTGTG", "-", "-"));
    }

    #[test]
    fn the_past_the_span_class_refuses_an_allele_no_longer_than_the_span() {
        assert!(!past_the_span_rows(
            "7:100998-101000",
            "GTG",
            "998-1000",
            "898-900"
        ));
        assert!(!past_the_span_rows(
            "7:100998-101000",
            "GT",
            "998-1000",
            "898-900"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_location_of_one_position() {
        assert!(!past_the_span_rows("7:101000", "GTGTG", "1000", "900"));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_location_of_two_neighbouring_positions() {
        assert!(!past_the_span_rows(
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            "899-900"
        ));
    }

    #[test]
    fn the_first_column_tells_two_replaced_bases_from_an_insertion_between_them() {
        let spelt = |name: &str| {
            let row = |terms: &str| {
                pair_row(&[
                    name,
                    TWO_POSITIONS,
                    "GTGTG",
                    "ENSG1",
                    "ENST00000635895",
                    "Transcript",
                    terms,
                    "999-1000",
                    "899-900",
                ])
            };
            splice_term_past_the_record_span_pair(
                &row("splice_donor_variant,frameshift_variant"),
                &row("frameshift_variant,splice_region_variant"),
            )
        };
        assert!(spelt("7_100999_TA/GTGTG"), "two bases replaced by five");
        assert!(!spelt("7_100999_T/TGTGTG"), "an insertion of five");
    }

    #[test]
    fn the_past_the_span_class_refuses_a_deletion() {
        assert!(!past_the_span_rows(
            "7:100998-101000",
            "-",
            "998-1000",
            "898-900"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_candidate_that_keeps_one_of_the_six() {
        assert!(!past_the_span_terms(
            "splice_donor_variant,frameshift_variant,splice_donor_region_variant",
            "frameshift_variant,splice_donor_region_variant"
        ));
        assert!(!past_the_span_terms(
            "splice_donor_variant,frameshift_variant,intron_variant",
            "frameshift_variant,intron_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_donor_region_term_that_comes() {
        assert!(!past_the_span_terms(
            "splice_donor_variant,frameshift_variant",
            "frameshift_variant,splice_donor_region_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_the_region_term_where_the_tract_term_alone_went() {
        assert!(!past_the_span_terms(
            "splice_polypyrimidine_tract_variant,frameshift_variant",
            "frameshift_variant,splice_region_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_a_further_differing_term() {
        assert!(!past_the_span_terms(
            "splice_donor_variant,frameshift_variant",
            "stop_gained,frameshift_variant,splice_region_variant"
        ));
        assert!(!past_the_span_terms(
            "splice_donor_variant,frameshift_variant,stop_lost",
            "frameshift_variant,splice_region_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_rows_whose_transcript_positions_differ() {
        let row =
            |terms: &str, cdna: &str, cds: &str| row("7:100998-101000", "GTGTG", terms, cdna, cds);
        let reference = row(
            "splice_donor_variant,frameshift_variant",
            "998-1000",
            "898-900",
        );
        let candidate = "frameshift_variant,splice_region_variant";
        assert!(!splice_term_past_the_record_span_pair(
            &reference,
            &row(candidate, "997-999", "898-900")
        ));
        assert!(!splice_term_past_the_record_span_pair(
            &reference,
            &row(candidate, "998-1000", "897-899")
        ));
    }

    #[test]
    fn the_past_the_span_class_refuses_the_terms_the_other_way_round() {
        assert!(!past_the_span_terms(
            "frameshift_variant,splice_region_variant",
            "splice_donor_variant,frameshift_variant"
        ));
    }

    #[test]
    fn the_past_the_span_class_takes_a_row_that_differs_by_a_stop_codon_classes_terms_as_well() {
        let row = |terms: &str| row("7:100998-101000", "GTGTG", terms, "998-1000", "-");
        let reference = row("splice_donor_variant,stop_retained_variant");
        let candidate = row("stop_lost,splice_region_variant");
        assert!(splice_term_past_the_record_span_pair(
            &reference, &candidate
        ));
        assert!(
            stop_codon_removed_without_stop_codon_test_pair(
                &row("stop_retained_variant,splice_region_variant"),
                &candidate
            ),
            "the rest of the difference is that class's"
        );
        assert!(
            !splice_term_past_the_record_span_pair(
                &reference,
                &row("stop_lost,missense_variant,splice_region_variant")
            ),
            "a further term outside both classes"
        );
    }

    /// The terms of a reference row whose donor, donor-region and intron terms come from
    /// positions past the span, and of vep-rs's row with `splice_region_variant` in their
    /// place.
    const DONOR_SIDE_TERMS: &str =
        "splice_donor_variant,frameshift_variant,splice_donor_region_variant,intron_variant";
    const REGION_TERM_BACK: &str = "frameshift_variant,splice_region_variant";

    /// A row of one transcript whose allele stands in for the reference bases `reference`.
    fn record_row(
        reference: &str,
        location: &str,
        allele: &str,
        terms: &str,
        cdna: &str,
    ) -> PairRow {
        PairRow {
            reference: bases(reference),
            ..row(location, allele, terms, cdna, "-")
        }
    }

    /// Whether the class takes the reference's and vep-rs's row of one transcript on the
    /// record that replaces `reference` at `location` by `allele`, their terms as given.
    fn past_the_span_record(
        reference: &str,
        location: &str,
        allele: &str,
        cdna: &str,
        terms: (&str, &str),
    ) -> bool {
        let row = |terms: &str| record_row(reference, location, allele, terms, cdna);
        splice_term_past_the_record_span_pair(&row(terms.0), &row(terms.1))
    }

    #[test]
    fn with_its_record_the_class_takes_a_span_inside_an_exon_of_two_positions_and_of_three() {
        let terms = (DONOR_SIDE_TERMS, REGION_TERM_BACK);
        assert!(past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            terms
        ));
        assert!(past_the_span_record(
            "TAC",
            "7:100998-101000",
            "GTGTGT",
            "998-1000",
            terms
        ));
        assert!(
            past_the_span_record(
                "TA",
                TWO_POSITIONS,
                "GTGT",
                "901-902",
                (
                    "splice_acceptor_variant,frameshift_variant",
                    REGION_TERM_BACK
                )
            ),
            "the acceptor term, on the reverse strand"
        );
    }

    #[test]
    fn with_its_record_the_class_reads_the_reference_bases_of_either_row() {
        let reference = row(TWO_POSITIONS, "GTGTG", DONOR_SIDE_TERMS, "999-1000", "-");
        let candidate = record_row("TA", TWO_POSITIONS, "GTGTG", REGION_TERM_BACK, "999-1000");
        assert!(splice_term_past_the_record_span_pair(
            &reference, &candidate
        ));
    }

    #[test]
    fn with_its_record_the_class_takes_a_span_across_an_exons_end() {
        assert!(past_the_span_record(
            "TA",
            "7:101000-101001",
            "GTGTG",
            "1000-?",
            (
                "splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant",
                "splice_donor_variant,coding_sequence_variant,intron_variant"
            )
        ));
    }

    #[test]
    fn with_its_record_the_class_takes_a_span_in_an_intron() {
        assert!(past_the_span_record(
            "TA",
            "7:101003-101004",
            "GTGT",
            "-",
            (
                "splice_donor_5th_base_variant,intron_variant",
                "splice_donor_region_variant,intron_variant"
            )
        ));
    }

    #[test]
    fn with_its_record_the_class_takes_the_intron_term_gone_one_to_three_bases_before_the_intron() {
        for (location, allele, cdna, intron_start) in [
            (TWO_POSITIONS, "GTGTG", "999-1000", 2),
            ("7:100998-100999", "GTGTGTG", "998-999", 3),
            ("7:100997-100998", "GTGTGTGT", "997-998", 4),
        ] {
            assert!(
                past_the_span_record(
                    "TA",
                    location,
                    allele,
                    cdna,
                    (DONOR_SIDE_TERMS, REGION_TERM_BACK)
                ),
                "{cdna}"
            );
            let (start, end) = cdna.split_once('-').unwrap();
            let place = record_span_placement(
                &[start, end, "-", "-"],
                &letters("TA"),
                &letters(allele),
                DONOR | DONOR_REGION | INTRON,
                REGION,
            );
            assert_eq!(
                place.map(|(forward, start, ..)| (forward, start)),
                Some((true, intron_start)),
                "{cdna}"
            );
        }
    }

    #[test]
    fn with_its_record_the_class_takes_a_donor_region_term_gone_beside_a_kept_donor_term() {
        assert!(past_the_span_record(
            "TA",
            "7:101001-101002",
            "GTGT",
            "-",
            (
                "splice_donor_variant,splice_donor_region_variant,intron_variant",
                "splice_donor_variant,intron_variant"
            )
        ));
    }

    #[test]
    fn with_its_record_the_class_refuses_a_pure_insertion() {
        assert!(!past_the_span_record(
            "",
            TWO_POSITIONS,
            "GTG",
            "999-1000",
            ("splice_donor_variant,frameshift_variant", REGION_TERM_BACK)
        ));
    }

    #[test]
    fn with_its_record_the_class_refuses_an_allele_as_long_as_the_reference_bases() {
        assert!(!past_the_span_record(
            "TAC",
            "7:100998-101000",
            "GTG",
            "998-1000",
            ("splice_donor_variant,frameshift_variant", REGION_TERM_BACK)
        ));
    }

    #[test]
    fn with_its_record_the_class_refuses_one_reference_base() {
        assert!(!past_the_span_record(
            "T",
            "7:101000",
            "GTG",
            "1000",
            ("splice_donor_variant,frameshift_variant", REGION_TERM_BACK)
        ));
    }

    /// The record replaces two bases across an exon's end by four letters whose second is the
    /// reference's: the intron's first base lies in the span at an index where the alleles
    /// agree, so vep-rs keeps the donor term.
    #[test]
    fn with_its_record_the_class_refuses_a_term_gone_whose_window_holds_a_base_of_the_span() {
        let record = |candidate: &str| {
            past_the_span_record(
                "TA",
                "7:101000-101001",
                "GAGT",
                "1000-?",
                (
                    "splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant",
                    candidate,
                ),
            )
        };
        assert!(record(
            "splice_donor_variant,coding_sequence_variant,intron_variant"
        ));
        assert!(!record("coding_sequence_variant,intron_variant"));
    }

    #[test]
    fn with_its_record_the_class_refuses_rows_no_place_of_the_intron_gives() {
        assert!(!past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            (DONOR_SIDE_TERMS, "frameshift_variant")
        ));
    }

    #[test]
    fn with_its_record_the_class_refuses_a_further_differing_term() {
        assert!(!past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            (
                DONOR_SIDE_TERMS,
                "stop_gained,frameshift_variant,splice_region_variant"
            )
        ));
    }

    #[test]
    fn with_its_record_the_class_refuses_rows_whose_transcript_positions_differ() {
        let row = |terms: &str, cdna: &str| record_row("TA", TWO_POSITIONS, "GTGTG", terms, cdna);
        assert!(!splice_term_past_the_record_span_pair(
            &row(DONOR_SIDE_TERMS, "999-1000"),
            &row(REGION_TERM_BACK, "499-500")
        ));
    }

    /// The places of an intron beside a span that two cDNA positions allow, in their order.
    fn places(span: i64, cdna: (&str, &str)) -> Vec<(bool, i64, i64)> {
        let mut places = Vec::new();
        intron_placements(
            span,
            2,
            &[cdna.0, cdna.1, "-", "-"],
            |forward, start, end| {
                places.push((forward, start, end));
                false
            },
        );
        places
    }

    #[test]
    fn two_numbers_as_far_apart_as_the_span_put_the_intron_after_the_span_or_before_it() {
        let places = places(3, ("998", "1000"));
        let starts = |strand: bool, after: bool| -> BTreeSet<i64> {
            let on = |(forward, start, _): &(bool, i64, i64)| {
                *forward == strand && (*start > 0) == after
            };
            places
                .iter()
                .filter(|place| on(place))
                .map(|place| place.1)
                .collect()
        };
        for strand in [true, false] {
            assert_eq!(
                starts(strand, true),
                (3..10).collect(),
                "one to seven bases on"
            );
            assert_eq!(starts(strand, false), (-16..0).collect(), "before the span");
        }
        assert!(places.contains(&(true, 3, 3)), "an intron of one base");
        assert!(
            places.contains(&(false, 9, 10)),
            "a two-base intron seven bases on"
        );
        assert!(
            places.contains(&(true, 3, FAR)),
            "an intron longer than any window"
        );
        assert!(
            places.contains(&(true, -2, -1)),
            "a two-base intron that ends beside the span"
        );
        let before = places.iter().filter(|place| place.1 < 0);
        assert!(before
            .clone()
            .all(|(_, start, end)| start <= end && *end < 0));
        assert_eq!(before.count(), 2 * (1..=16).sum::<usize>());
        assert_eq!(places.len(), 1406);
    }

    #[test]
    fn two_numbers_closer_than_the_span_put_an_intron_of_the_missing_length_inside_it() {
        assert_eq!(
            places(5, ("998", "1000")),
            [(true, 1, 2), (true, 2, 3), (false, 1, 2), (false, 2, 3)],
            "two bases of five, the span's ends kept"
        );
        assert_eq!(
            places(3, ("999", "1000")),
            [(true, 1, 1), (false, 1, 1)],
            "the middle base of three"
        );
        assert_eq!(
            places(2, ("1000", "1000")),
            [],
            "no base inside a span of two"
        );
        assert_eq!(
            places(3, ("997", "1000")),
            [],
            "numbers further apart than the span"
        );
        assert_eq!(
            places(3, ("1000", "998")),
            [],
            "the second number below the first"
        );
    }

    /// Eight bases inside an intron that begins two bases before the span, on the reverse
    /// strand: the span's first six bases lie in a window of `splice_region` and all of them
    /// in the tract's.
    #[test]
    fn the_last_region_assigned_decides_the_references_region_term_and_any_region_vep_rs_s() {
        let rows = |regions: &[(i64, i64)]| intron_rows(regions, 8, -2, FAR, false, false);
        let two_runs = differing_regions(&letters("TACGATCA"), &letters("GACGATCAGT"));
        assert_eq!(two_runs, [(0, 0), (8, 9)]);
        let (reference, candidate) = rows(&two_runs);
        assert_eq!(
            (reference, candidate),
            (TRACT | INTRON, TRACT | INTRON | REGION),
            "the first run in a window and the last in none"
        );
        assert_eq!(candidate & TRACT, TRACT, "the span holds the tract's bases");
        assert_eq!(
            rows(&[(0, 0), (4, 9)]),
            (TRACT | INTRON | REGION, TRACT | INTRON | REGION),
            "the last run in a window"
        );
        assert_eq!(
            rows(&[(7, 9)]),
            (TRACT | INTRON, TRACT | INTRON),
            "no run in a window"
        );
        assert_eq!(
            intron_rows(&[(0, 0), (2, 2), (6, 8)], 6, 2, 101, true, false),
            (
                DONOR | DONOR_5TH_BASE | INTRON,
                DONOR | DONOR_REGION | INTRON
            ),
            "a run on the intron's first two bases between two runs in a window"
        );
    }

    /// Eight bases of an intron replaced by eleven letters alike at four of them, five runs in
    /// all: the last run lies in a window of `splice_region` and the span holds bases of the
    /// tract's window at every place that gives the reference's row.
    #[test]
    fn with_its_record_the_class_refuses_the_tract_term_gone_from_a_span_of_five_runs() {
        let (reference, allele) = ("TACGATCA", "GATGCTGAGTG");
        assert_eq!(
            differing_regions(&letters(reference), &letters(allele)),
            [(0, 0), (2, 2), (4, 4), (6, 6), (8, 10)]
        );
        assert!(!past_the_span_record(
            reference,
            "7:101007-101014",
            allele,
            "-",
            (
                "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
                "splice_region_variant,intron_variant"
            )
        ));
    }

    const FIFTH_BASE_TERMS: &str = "splice_donor_5th_base_variant,frameshift_variant";

    /// The windows of an intron of two bases run on past its last base: its fifth base is the
    /// first base past a span of two that follows it.
    #[test]
    fn with_its_record_the_class_takes_a_fifth_base_term_gone_beside_a_two_base_intron_before_the_span(
    ) {
        assert_eq!(
            intron_rows(&[(0, 4)], 2, -2, -1, true, false),
            (DONOR_5TH_BASE, DONOR_REGION)
        );
        assert!(past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            (
                FIFTH_BASE_TERMS,
                "splice_donor_region_variant,frameshift_variant"
            )
        ));
    }

    /// On the reverse strand the fifth base of an intron of two bases that begins six bases
    /// after a span of two lies three bases past the span, and the eleven bases around the
    /// intron's end hold the span's last base.
    #[test]
    fn with_its_record_the_class_takes_a_fifth_base_term_gone_beside_a_two_base_intron_six_bases_after_the_span(
    ) {
        assert_eq!(
            intron_rows(&[(0, 5)], 2, 7, 8, false, false),
            (DONOR_5TH_BASE, REGION)
        );
        assert!(past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTGT",
            "999-1000",
            (FIFTH_BASE_TERMS, REGION_TERM_BACK)
        ));
    }

    /// Six positions whose two cDNA positions are four apart hold an intron of two bases: the
    /// alleles agree on it, and its fifth base is the first base past the span.
    #[test]
    fn with_its_record_the_class_takes_a_fifth_base_term_gone_beside_an_intron_inside_the_span() {
        let record = |cdna: &str| {
            past_the_span_record(
                "TACGAT",
                "7:100999-101004",
                "GACGATGT",
                cdna,
                (
                    "splice_donor_5th_base_variant,coding_sequence_variant",
                    "splice_donor_region_variant,coding_sequence_variant",
                ),
            )
        };
        assert_eq!(
            intron_rows(&[(0, 0), (6, 7)], 6, 2, 3, true, false),
            (DONOR_5TH_BASE, DONOR_REGION)
        );
        assert!(record("999-1002"));
        assert!(!record("999-998"), "the second number below the first");
    }

    #[test]
    fn an_allele_that_opens_with_its_whole_reference_allele_writes_an_insertion_after_the_span() {
        let after_the_span = |reference: &str, allele: &str| {
            insertion_after_the_span(
                &differing_regions(&letters(reference), &letters(allele)),
                reference.len() as i64,
            )
        };
        assert!(after_the_span("TA", "TAGTG"));
        assert!(after_the_span("TA", "TAG"), "one further base");
        assert!(after_the_span("TAC", "TACGT"), "three reference bases");
        assert!(!after_the_span("T", "TGTG"), "one reference base");
        assert!(!after_the_span("TA", "TGGTG"), "the second base differs");
        assert!(!after_the_span("TA", "GAGTG"), "the first base differs");
        assert!(!after_the_span("TA", "TA"), "no further base");
        assert!(!after_the_span("TAC", "TA"), "a shorter allele");
    }

    /// A span of two bases beside one intron, under an allele that holds further bases after
    /// them: the reference's terms come from the bases after the span, the candidate's from the
    /// insertion between the span's last base and the base after it.
    #[test]
    fn an_insertion_after_the_span_is_read_by_its_two_neighbours() {
        let three_further = |start: i64, end: i64, forward: bool| {
            intron_rows(&[(2, 4)], 2, start, end, forward, false)
        };
        assert_eq!(
            three_further(-15, FAR, false),
            (INTRON, INTRON | TRACT),
            "after the tract's last base, on the reverse strand"
        );
        assert_eq!(
            three_further(-16, FAR, false),
            (INTRON, INTRON),
            "one base past the tract's last base"
        );
        assert_eq!(
            three_further(-FAR, 18, true),
            (INTRON | TRACT, INTRON | TRACT),
            "before the tract's first base, on the forward strand"
        );
        assert_eq!(
            three_further(-FAR, 3, true),
            (ACCEPTOR, INTRON | TRACT | REGION),
            "between the tract's innermost base and the acceptor bases"
        );
        assert_eq!(
            three_further(-FAR, 3, false),
            (DONOR, INTRON | REGION),
            "between the third intronic base and the donor bases, on the reverse strand"
        );
        assert_eq!(
            three_further(1, FAR, true),
            (DONOR | DONOR_REGION | INTRON, DONOR),
            "between the two donor bases"
        );
        assert_eq!(
            three_further(1, FAR, false),
            (ACCEPTOR | INTRON, ACCEPTOR),
            "between the two acceptor bases, on the reverse strand"
        );
        let two_further = |start: i64| intron_rows(&[(2, 3)], 2, start, FAR, true, false);
        assert_eq!(
            two_further(2),
            (DONOR, REGION),
            "between an exon's last base and the donor bases"
        );
        assert_eq!(
            two_further(0),
            (DONOR_REGION | INTRON, INTRON | REGION),
            "after the two donor bases"
        );
        assert_eq!(
            intron_rows(&[(0, 0)], 1, -15, FAR, false, false),
            (INTRON | TRACT, INTRON | TRACT),
            "one reference base is read as a change of that base"
        );
    }

    /// The reference's terms and the candidate's at the three places where the insertion's own
    /// neighbours give the candidate a term the bases after the span do not give the reference.
    const AFTER_THE_SPAN_PLACES: [(&str, (&str, &str)); 3] = [
        (
            "after the tract's last base, on the reverse strand",
            (
                "intron_variant",
                "splice_polypyrimidine_tract_variant,intron_variant",
            ),
        ),
        (
            "between the tract's innermost base and the acceptor bases",
            (
                "splice_acceptor_variant",
                "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
            ),
        ),
        (
            "between the third intronic base and the donor bases, on the reverse strand",
            (
                "splice_donor_variant",
                "splice_region_variant,intron_variant",
            ),
        ),
    ];

    /// Whether the class takes the two rows of an allele at the two positions of an intron
    /// whose reference bases are `reference`, their terms as given.
    fn after_the_span_record(reference: &str, allele: &str, terms: (&str, &str)) -> bool {
        past_the_span_record(reference, TWO_POSITIONS, allele, "-", terms)
    }

    #[test]
    fn with_its_record_the_class_takes_an_insertion_after_the_span_at_the_places_its_terms_differ()
    {
        for (place, terms) in AFTER_THE_SPAN_PLACES {
            assert!(after_the_span_record("TA", "TAGTG", terms), "{place}");
            assert!(
                after_the_span_record("TA", "TAG", terms),
                "{place}, one further base"
            );
        }
        assert!(
            past_the_span_record(
                "TAC",
                "7:100998-101000",
                "TACGTG",
                "-",
                AFTER_THE_SPAN_PLACES[0].1
            ),
            "three reference bases"
        );
    }

    #[test]
    fn with_its_record_the_class_takes_the_terms_an_insertion_between_two_site_bases_withholds() {
        assert!(
            past_the_span_record(
                "TA",
                TWO_POSITIONS,
                "TAGTG",
                "1000-?",
                (
                    "splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant",
                    "splice_donor_variant,coding_sequence_variant"
                )
            ),
            "between the two donor bases"
        );
        assert!(
            past_the_span_record(
                "TA",
                TWO_POSITIONS,
                "TAGTG",
                "?-1000",
                (
                    "splice_acceptor_variant,coding_sequence_variant,intron_variant",
                    "splice_acceptor_variant,coding_sequence_variant"
                )
            ),
            "between the two acceptor bases, on the reverse strand"
        );
    }

    /// The geometry of line 166 of `tests/golden/116/TAIR10-arabidopsis/variants.vcf`, in
    /// made-up bases: three reference bases, a first alternate allele of one base and a second
    /// that adds two bases to the reference allele. The second's rows stand at the two positions
    /// after the shared first base, an exon's last two bases on one transcript and the two donor
    /// bases of an intron on another.
    #[test]
    fn with_its_record_the_class_takes_both_rows_of_a_second_alternate_that_opens_with_its_reference(
    ) {
        let row = |name: &str, terms: &str, cdna: &str, cds: &str| {
            pair_row(&[
                name,
                "4:10523300-10523301",
                "TAGT",
                "AT4G19230",
                "AT4G19230.2",
                "Transcript",
                terms,
                cdna,
                cds,
            ])
        };
        let spelt = "4_10523299_CTA/C/CTAGT";
        let exon = |name: &str, terms: &str| row(name, terms, "1444-1445", "1299-1300");
        assert!(splice_term_past_the_record_span_pair(
            &exon(spelt, "splice_donor_variant,frameshift_variant"),
            &exon(spelt, "frameshift_variant,splice_region_variant")
        ));
        assert!(
            !splice_term_past_the_record_span_pair(
                &exon(spelt, "splice_donor_variant,frameshift_variant"),
                &exon(spelt, "frameshift_variant")
            ),
            "a candidate's row without the insertion's splice_region_variant"
        );
        let intron = |terms: &str| row(spelt, terms, "-", "-");
        assert!(splice_term_past_the_record_span_pair(
            &intron("splice_donor_region_variant,intron_variant"),
            &intron("intron_variant,splice_region_variant")
        ));
        assert!(
            !splice_term_past_the_record_span_pair(
                &exon("rs1", "splice_donor_variant,frameshift_variant"),
                &exon("rs1", "frameshift_variant,splice_region_variant")
            ),
            "two positions under an identifier, without the record"
        );
    }

    #[test]
    fn an_allele_that_opens_with_its_one_reference_base_writes_an_insertion_after_it() {
        let after_one_base = |reference: &str, allele: &str| {
            insertion_after_one_base(&letters(reference), &letters(allele))
        };
        assert!(after_one_base("T", "TGTG"));
        assert!(after_one_base("T", "TG"), "one further base");
        assert!(!after_one_base("T", "GTG"), "the base differs");
        assert!(!after_one_base("T", "T"), "no further base");
        assert!(!after_one_base("T", "-"), "a deletion");
        assert!(!after_one_base("TA", "TAG"), "two reference bases");
        assert!(!after_one_base("", "GTG"), "a pure insertion");
    }

    /// The reference's terms and the candidate's for an allele that opens with its one
    /// reference base, at the places where the base read as changed and the insertion after
    /// it give different terms, with the cDNA position of the base.
    const AFTER_ONE_BASE_PLACES: [(&str, &str, (&str, &str)); 10] = [
        (
            "the base before the tract's outermost base",
            "-",
            (
                "intron_variant",
                "splice_polypyrimidine_tract_variant,intron_variant",
            ),
        ),
        (
            "the intron's last base, before the exon",
            "-",
            ("splice_acceptor_variant", "splice_region_variant"),
        ),
        (
            "the exon's third base",
            "3",
            (
                "splice_region_variant,non_coding_transcript_exon_variant",
                "non_coding_transcript_exon_variant",
            ),
        ),
        (
            "the donor's second base",
            "-",
            (
                "splice_donor_variant",
                "splice_region_variant,intron_variant",
            ),
        ),
        (
            "the intron's fifth base",
            "-",
            (
                "splice_donor_5th_base_variant,intron_variant",
                "splice_donor_region_variant,intron_variant",
            ),
        ),
        (
            "the intron's sixth base",
            "-",
            (
                "splice_donor_region_variant,intron_variant",
                "splice_region_variant,intron_variant",
            ),
        ),
        (
            "the intron's eighth base",
            "-",
            ("splice_region_variant,intron_variant", "intron_variant"),
        ),
        (
            "the acceptor's base beside the tract, on the reverse strand",
            "-",
            (
                "splice_acceptor_variant",
                "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
            ),
        ),
        (
            "the eighth base from the intron's 3' end, on the reverse strand",
            "-",
            (
                "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
                "splice_polypyrimidine_tract_variant,intron_variant",
            ),
        ),
        (
            "the donor's first base, before the exon on the reverse strand",
            "-",
            ("splice_donor_variant", "splice_region_variant"),
        ),
    ];

    #[test]
    fn with_its_record_the_class_takes_an_insertion_after_one_base_at_the_places_its_terms_differ()
    {
        for (place, cdna, terms) in AFTER_ONE_BASE_PLACES {
            assert!(
                past_the_span_record("T", "7:101000", "TGTG", cdna, terms),
                "{place}"
            );
            assert!(
                past_the_span_record("T", "7:101000", "TG", cdna, terms),
                "{place}, one further base"
            );
            assert!(
                !past_the_span_record("T", "7:101000", "GTGT", cdna, terms),
                "{place}, an allele that does not open with the reference base"
            );
            let row = |terms: &str| row("7:101000", "TGTG", terms, cdna, "-");
            assert!(
                !splice_term_past_the_record_span_pair(&row(terms.0), &row(terms.1)),
                "{place}, without the record"
            );
        }
    }

    /// A record of two alternate alleles as the first column spells it: the second opens with
    /// the reference allele and adds three bases, and its row stands at the one position after
    /// the shared first base, the donor's second base.
    #[test]
    fn with_its_record_the_class_takes_the_row_of_an_alternate_that_opens_with_one_kept_base() {
        let row = |name: &str, terms: &str| {
            pair_row(&[
                name,
                "22:20996968",
                "TGTG",
                "ENSG00000099949",
                "ENST00000646124",
                "Transcript",
                terms,
                "-",
                "-",
            ])
        };
        let spelt = "22_20996967_CT/C/CTGTG";
        assert!(splice_term_past_the_record_span_pair(
            &row(spelt, "splice_donor_variant"),
            &row(spelt, "splice_region_variant,intron_variant")
        ));
        assert!(
            !splice_term_past_the_record_span_pair(
                &row(spelt, "splice_donor_variant"),
                &row(spelt, "splice_donor_region_variant,intron_variant")
            ),
            "a candidate's row with a term no place of the insertion gives beside the donor term"
        );
        assert!(
            !splice_term_past_the_record_span_pair(
                &row("rs1", "splice_donor_variant"),
                &row("rs1", "splice_region_variant,intron_variant")
            ),
            "one position under an identifier, without the record"
        );
    }

    #[test]
    fn with_its_record_the_class_refuses_a_term_an_insertion_after_one_base_does_not_explain() {
        let refused = [
            (
                "a donor term the candidate alone holds",
                "-",
                "intron_variant",
                "intron_variant,splice_donor_variant",
            ),
            (
                "another term the candidate alone holds",
                "-",
                "splice_acceptor_variant",
                "splice_region_variant,missense_variant",
            ),
            (
                "another term the reference alone holds",
                "-",
                "splice_donor_variant,non_coding_transcript_variant",
                "splice_region_variant,intron_variant",
            ),
            (
                "the intron term the insertion's place gives, missing after the donor's second base",
                "-",
                "splice_donor_variant,intron_variant",
                "splice_region_variant",
            ),
            (
                "splice_region_variant gone from a base inside the transcript's sequence, with the intron term",
                "3",
                "splice_region_variant,intron_variant",
                "intron_variant",
            ),
            (
                "the same terms on both rows",
                "-",
                "splice_region_variant,intron_variant",
                "splice_region_variant,intron_variant",
            ),
        ];
        for (shape, cdna, reference, candidate) in refused {
            assert!(
                !past_the_span_record("T", "7:101000", "TGTG", cdna, (reference, candidate)),
                "{shape}"
            );
        }
    }

    /// The rows of `--terms display`: both site terms print one name, the fifth base, the
    /// donor region, the tract and `splice_region_variant` another, the intron term a third.
    #[test]
    fn under_the_display_names_the_class_reads_the_names_the_seven_terms_print() {
        assert!(
            past_the_span_record(
                "T",
                "7:101000",
                "TT",
                "3",
                ("FRAMESHIFT_CODING,SPLICE_SITE", "FRAMESHIFT_CODING")
            ),
            "the exon's third base, one base added after it"
        );
        assert!(
            past_the_span_record(
                "T",
                "7:101000",
                "TGTG",
                "-",
                ("ESSENTIAL_SPLICE_SITE", "SPLICE_SITE,INTRONIC")
            ),
            "the donor's second base, three bases added after it"
        );
        assert!(
            past_the_span_record(
                "TA",
                TWO_POSITIONS,
                "GTGTG",
                "999-1000",
                (
                    "ESSENTIAL_SPLICE_SITE,FRAMESHIFT_CODING,SPLICE_SITE,INTRONIC",
                    "FRAMESHIFT_CODING,SPLICE_SITE"
                )
            ),
            "two replaced bases inside an exon, the donor and intron names gone"
        );
        let refused = [
            (
                "the same names on both rows, the fifth base against the donor region",
                "-",
                "SPLICE_SITE,INTRONIC",
                "SPLICE_SITE,INTRONIC",
            ),
            (
                "another term that differs",
                "3",
                "FRAMESHIFT_CODING,SPLICE_SITE",
                "STOP_GAINED",
            ),
            (
                "names no place of the intron gives",
                "-",
                "INTRONIC",
                "ESSENTIAL_SPLICE_SITE",
            ),
            (
                "a Sequence Ontology name beside a display name",
                "-",
                "ESSENTIAL_SPLICE_SITE",
                "splice_region_variant,INTRONIC",
            ),
            (
                "the name of the donor and the acceptor term on the candidate's row alone",
                "-",
                "WITHIN_NON_CODING_GENE",
                "ESSENTIAL_SPLICE_SITE,WITHIN_NON_CODING_GENE",
            ),
        ];
        for (shape, cdna, reference, candidate) in refused {
            assert!(
                !past_the_span_record("T", "7:101000", "TGTG", cdna, (reference, candidate)),
                "{shape}"
            );
        }
        assert!(
            !past_the_span_record(
                "TAC",
                "7:100998-101000",
                "GAT",
                "-",
                ("INTRONIC", "SPLICE_SITE,INTRONIC")
            ),
            "three replaced bases under an allele as long, its first and last letter differing"
        );
        assert!(
            !past_the_span_record(
                "TAC",
                "7:100998-101000",
                "GACTG",
                "-",
                ("INTRONIC", "SPLICE_SITE,INTRONIC")
            ),
            "a longer allele that does not open with its reference, no name gone from the reference's row"
        );
        let unread = |terms: &str| row("7:101000", "TT", terms, "3", "-");
        assert!(
            !splice_term_past_the_record_span_pair(
                &unread("FRAMESHIFT_CODING,SPLICE_SITE"),
                &unread("FRAMESHIFT_CODING")
            ),
            "without the record"
        );
    }

    #[test]
    fn with_its_record_the_class_refuses_an_allele_that_does_not_open_with_its_whole_reference() {
        for (place, terms) in AFTER_THE_SPAN_PLACES {
            assert!(
                !after_the_span_record("TA", "TGGTG", terms),
                "{place}, the second base differs"
            );
            assert!(
                !after_the_span_record("TA", "GAGTG", terms),
                "{place}, the first base differs"
            );
            let row = |terms: &str| row(TWO_POSITIONS, "TAGTG", terms, "-", "-");
            assert!(
                !splice_term_past_the_record_span_pair(&row(terms.0), &row(terms.1)),
                "{place}, without the record"
            );
        }
    }

    #[test]
    fn with_its_record_the_class_refuses_a_term_an_insertion_after_the_span_does_not_explain() {
        let refused = [
            (
                "a donor term the candidate alone holds",
                "intron_variant",
                "splice_polypyrimidine_tract_variant,intron_variant,splice_donor_variant",
            ),
            (
                "splice_region_variant the reference alone holds",
                "splice_region_variant,intron_variant",
                "splice_polypyrimidine_tract_variant,intron_variant",
            ),
            (
                "another term the candidate alone holds",
                "splice_acceptor_variant",
                "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant,missense_variant",
            ),
            (
                "another term the reference alone holds",
                "splice_donor_variant,non_coding_transcript_variant",
                "splice_region_variant,intron_variant",
            ),
            (
                "the intron term the insertion's place gives, missing beside the acceptor bases",
                "splice_acceptor_variant",
                "splice_region_variant,splice_polypyrimidine_tract_variant",
            ),
            (
                "the intron term the insertion's place gives, missing beside the donor bases",
                "splice_donor_variant",
                "splice_region_variant",
            ),
            (
                "the same terms on both rows",
                "intron_variant",
                "intron_variant",
            ),
        ];
        for (shape, reference, candidate) in refused {
            assert!(
                !after_the_span_record("TA", "TAGTG", (reference, candidate)),
                "{shape}"
            );
        }
    }
}

#[cfg(test)]
#[path = "row_pair_rule_rows_tests.rs"]
mod pair_rule_rows_tests;

#[cfg(test)]
#[path = "row_number_bound_tests.rs"]
mod number_bound_tests;

#[cfg(test)]
#[path = "row_coding_end_tests.rs"]
mod coding_end_tests;

#[cfg(test)]
#[path = "row_settled_stop_codon_tests.rs"]
mod settled_stop_codon_tests;

#[cfg(test)]
#[path = "row_frameshift_inside_a_held_stop_codon_tests.rs"]
mod frameshift_inside_a_held_stop_codon_tests;

#[cfg(test)]
#[path = "row_removed_start_codon_tests.rs"]
mod removed_start_codon_tests;

#[cfg(test)]
#[path = "row_kept_start_codon_tests.rs"]
mod kept_start_codon_tests;
