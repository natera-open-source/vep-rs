// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a comparison hands a row rule that reads the coding ends of the row's transcript,
//! and the row rules of the classes of pair at a coding end.
//!
//! Which pair a class takes is decided in `coding_ends::terms` (`class_of_pair`,
//! `landed_rule_may_take`); this module reads an entry into the plain values those functions
//! ask and looks the reading of the transcript's cache record up for the row.

use crate::coding_ends::terms::{
    amino_acids_on_a_completed_stop_codon, class_of_pair, codons_on_a_completed_stop_codon,
    count_alone_differs_at_the_first_residue, frameshift_counted_at_the_first_residue,
    landed_rule_may_take, lost_start_range_of_a_frameshift, position_cell,
    protein_notation_on_a_completed_stop_codon, protein_notation_where_no_start_codon_is_shown,
    Cells, CodingEndClass, LandedRule, RowFacts, WhatWasRead,
};
use crate::coding_ends::CodingEnds;
use crate::consts::{AMINO_ACIDS_FIELDS, HGVSP_FIELDS, HGVS_OFFSET_FIELDS};
use crate::entry::Entry;
use crate::error::Result;
use crate::fasta::ReferenceSequence;
use crate::fields::field_named;
use crate::rules::downstream::field_of;
use crate::rules::edited_model::CODONS_FIELDS;
use crate::rules::following::{terms_impact, terms_write_field};
use crate::rules::partial_codon_models::{location_span, PartialCodonModels};
use crate::rules::term_pairs::IMPACT_FIELDS;
use crate::self_parity::most_severe;
use crate::shared::{blank, consequence_set, feature_kind};
use std::collections::{BTreeSet, HashMap, HashSet};

/// What a comparison read of the coding ends of one row's transcript.
#[derive(Clone, Copy)]
pub enum CodingEndsRead<'a> {
    /// The comparison was given no cache: a rule that asks the coding ends takes nothing,
    /// and a rule that stood before the coding ends were read stands as it stood.
    NoCache,
    /// The comparison was given a cache, and the cache answers the row with no one record:
    /// no rule that asks the coding ends takes the pair.
    NoRecord,
    /// The reading of the cache record the row is answered with.
    Record(&'a CodingEnds),
}

/// What a row rule is asked with, beside the two entries of the pair.
#[derive(Clone, Copy)]
pub struct AtACodingEnd<'a> {
    pub ends: CodingEndsRead<'a>,
    /// The reference base at a genomic position of the row's sequence, `None` past the
    /// sequence's ends; `None` where the comparison was given no reference sequence, or
    /// one that does not hold the row's sequence.
    pub reference_bases: Option<&'a dyn Fn(i64) -> Option<u8>>,
}

impl AtACodingEnd<'static> {
    /// A pair of a comparison that was given no cache and no reference sequence.
    pub const NO_CACHE: AtACodingEnd<'static> = AtACodingEnd {
        ends: CodingEndsRead::NoCache,
        reference_bases: None,
    };
}

/// A row rule: whether the reference entry and the candidate entry of one key are a pair
/// of its class, on what the comparison read of the row's transcript.
pub type RowPair = fn(&Entry, &Entry, &AtACodingEnd<'_>) -> Result<bool>;

/// Asks `rule` with what the comparison holds for the row of the reference entry `e`: the
/// reading of the record of its transcript in the candidate's cache (`models`), looked up
/// by the record whose span holds the row's Location (`PartialCodonModels::coding_ends`),
/// and the bases of the row's sequence in `reference`. An entry that is no Transcript's
/// is answered with no record.
pub fn asked_at<T>(
    models: Option<&PartialCodonModels>,
    reference: Option<&ReferenceSequence>,
    e: &Entry,
    rule: impl FnOnce(&AtACodingEnd<'_>) -> Result<T>,
) -> Result<T> {
    let record = match models {
        Some(models) if transcript_row(e) => models.coding_ends(e)?,
        _ => None,
    };
    let ends = match (&record, models) {
        (Some(record), _) => CodingEndsRead::Record(record),
        (None, Some(_)) => CodingEndsRead::NoRecord,
        (None, None) => CodingEndsRead::NoCache,
    };
    let sequence = location_of(e)
        .and_then(|location| location.rsplit_once(':'))
        .map(|(sequence, _)| sequence);
    let held = reference
        .zip(sequence)
        .and_then(|(reference, sequence)| Some((reference, sequence, reference.length(sequence)?)));
    let base = |position: i64| {
        let (reference, sequence, length) = held?;
        if position < 1 || position > length {
            return None;
        }
        reference.bases(sequence, position, position).bytes().next()
    };
    let reference_bases: Option<&dyn Fn(i64) -> Option<u8>> = match held {
        Some(_) => Some(&base),
        None => None,
    };
    rule(&AtACodingEnd {
        ends,
        reference_bases,
    })
}

/// Whether the entry is a Transcript's: the last segment of its key, read up to a `#`, is
/// a Feature_type of the transcript list (`feature_kind`).
fn transcript_row(e: &Entry) -> bool {
    let last = e.key.rsplit('|').next().unwrap_or("");
    let feature_type = last.split('#').next().unwrap_or("");
    feature_kind(feature_type) == "transcript"
}

/// The entry's Location: a default or tab row's own, a VCF or JSON entry's the one its
/// record gives its allele (`Entry.allele_location`).
fn location_of(e: &Entry) -> Option<&str> {
    e.location
        .as_deref()
        .filter(|location| !location.is_empty())
        .or(e.allele_location.as_deref())
}

/// The five cells of a row that the derivation reads, held as the default format writes
/// them.
struct RowCells {
    cdna_position: String,
    cds_position: String,
    protein_position: String,
    amino_acids: String,
    codons: String,
}

impl RowCells {
    fn cells(&self) -> Cells<'_> {
        Cells {
            cdna_position: &self.cdna_position,
            cds_position: &self.cds_position,
            protein_position: &self.protein_position,
            amino_acids: &self.amino_acids,
            codons: &self.codons,
        }
    }
}

/// The five cells of an entry. A default or tab row and a CSQ entry carry them as the
/// columns `cDNA_position`, `CDS_position`, `Protein_position`, `Amino_acids` and `Codons`:
/// `None` where the output's header names fewer than the five (`--fields`). A JSON entry
/// carries each position as two members (`cdna_start` and `cdna_end`, and their kin),
/// read into one cell (`position_cell`), and lacks the member of a cell it leaves blank.
fn row_cells(e: &Entry) -> Option<RowCells> {
    let field = |name: &str| e.fields.get(name).map(String::as_str);
    if e.fields.contains_key("variant_allele") {
        let position = |start: &str, end: &str| position_cell(field(start), field(end));
        return Some(RowCells {
            cdna_position: position("cdna_start", "cdna_end"),
            cds_position: position("cds_start", "cds_end"),
            protein_position: position("protein_start", "protein_end"),
            amino_acids: field("amino_acids").unwrap_or("-").to_string(),
            codons: field("codons").unwrap_or("-").to_string(),
        });
    }
    Some(RowCells {
        cdna_position: field("cDNA_position")?.to_string(),
        cds_position: field("CDS_position")?.to_string(),
        protein_position: field("Protein_position")?.to_string(),
        amino_acids: field("Amino_acids")?.to_string(),
        codons: field("Codons")?.to_string(),
    })
}

/// Asks `question` with the plain values `coding_ends::terms` reads of the pair of the
/// reference entry `e` and the candidate entry `a`: the reference's row as the pair's own
/// facts, the candidate's terms, and what the comparison read beside them. The reference
/// length is that of the bases a record gave either entry (`Entry.reference`), which a
/// VCF record states on the forward strand. `None` where the reference's Location names
/// no position.
fn with_the_pair<T>(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
    question: impl FnOnce(&RowFacts<'_>, &str, &WhatWasRead<'_>) -> T,
) -> Option<T> {
    let location = location_of(e).and_then(location_span)?;
    let (theirs, ours) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    let row = RowFacts {
        transcript_row: transcript_row(e),
        allele: &e.allele,
        location,
        ensembl_vep_terms: &theirs,
    };
    let cells = row_cells(e).zip(row_cells(a));
    let reference = e.reference.as_ref().or(a.reference.as_ref());
    let read = WhatWasRead {
        cells: cells
            .as_ref()
            .map(|(theirs, ours)| (theirs.cells(), ours.cells())),
        reference_length: reference.map(|bases| bases.chars().count() as u64),
        reference_bases: at.reference_bases,
        record_on_the_reverse_strand: false,
    };
    Some(question(&row, &ours, &read))
}

/// Whether the pair is the class `class`'s: the first class of a coding end, in the order
/// they are asked, that takes the two rows on the reading of the transcript's cache record
/// (`coding_ends::terms::class_of_pair`). A class takes a pair only where vep-rs's coding
/// terms are the ones derived for the row from that reading, so a comparison that was
/// given no cache, and a row the cache answers with no one record, are no class's.
pub fn coding_end_class_pair(
    class: CodingEndClass,
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> bool {
    let CodingEndsRead::Record(ends) = at.ends else {
        return false;
    };
    with_the_pair(e, a, at, |row, ours, read| {
        class_of_pair(row, ours, Some(ends), read) == Some(class)
    })
    .unwrap_or(false)
}

/// The class of a coding end that takes the pair of the reference entry `e` and the
/// candidate entry `a` on the reading of the transcript's cache record, `None` where none
/// takes it: a comparison that was given no cache, and a row the cache answers with no one
/// record, have none.
pub fn coding_end_class(e: &Entry, a: &Entry, at: &AtACodingEnd<'_>) -> Option<CodingEndClass> {
    let CodingEndsRead::Record(ends) = at.ends else {
        return None;
    };
    with_the_pair(e, a, at, |row, ours, read| {
        class_of_pair(row, ours, Some(ends), read)
    })
    .flatten()
}

/// `structural_deletion_stop_lost_without_stop_codon`: on a transcript at whose coding end
/// no stop codon stands, or whose stop codon polyadenylation completes where the span holds
/// none of the one or two bases that begin it, Ensembl VEP writes `stop_lost` for a
/// structural deletion over the last three genomic positions of the coding region
/// (Utils/VariationEffect.pm:1268-1280) and vep-rs writes the same row without it
/// (`CodingEndClass::takes`).
pub fn structural_deletion_stop_lost_without_stop_codon_pair(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    let class = CodingEndClass::StructuralDeletionStopLost;
    Ok(coding_end_class_pair(class, e, a, at))
}

/// `stop_codon_completed_by_polyadenylation`: on a transcript whose stop codon
/// polyadenylation completes, Ensembl VEP reads the row on an incomplete terminal codon
/// (Utils/VariationEffect.pm:1588-1614) and vep-rs on the completed coding sequence
/// (`CodingEndClass::takes`).
pub fn stop_codon_completed_by_polyadenylation_pair(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    Ok(coding_end_class_pair(
        CodingEndClass::StopCodonCompleted,
        e,
        a,
        at,
    ))
}

/// `stop_codon_overlap_without_stop_codon`: on a transcript at whose coding end no stop
/// codon stands and that carries no `cds_end_NF`, Ensembl VEP writes a stop term from the
/// overlap of the allele with the last three coding positions
/// (Utils/VariationEffect.pm:1373-1381, 1400, 1408-1427) and vep-rs the terms of the row
/// without one (`CodingEndClass::takes`).
pub fn stop_codon_overlap_without_stop_codon_pair(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    Ok(coding_end_class_pair(
        CodingEndClass::StopCodonOverlap,
        e,
        a,
        at,
    ))
}

/// `stop_codon_behind_last_residue_without_stop_codon`: on a transcript at whose coding end
/// no stop codon stands, Ensembl VEP reads a stop codon an allele brings behind the last
/// residue as the transcript's own kept (Utils/VariationEffect.pm:1344, 1357) and vep-rs
/// as gained (`CodingEndClass::takes`).
pub fn stop_codon_behind_last_residue_without_stop_codon_pair(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    let class = CodingEndClass::StopCodonBehindLastResidue;
    Ok(coding_end_class_pair(class, e, a, at))
}

/// `start_codon_overlap_without_start_codon`: on a transcript at whose coding start no start
/// codon is shown and that carries no `cds_start_NF`, Ensembl VEP writes a start term from
/// the overlap of the allele with the first three coding positions
/// (Utils/VariationEffect.pm:975-986) and vep-rs the terms of the row without one
/// (`CodingEndClass::takes`).
pub fn start_codon_overlap_without_start_codon_pair(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    Ok(coding_end_class_pair(
        CodingEndClass::StartCodonOverlap,
        e,
        a,
        at,
    ))
}

/// The row rules of the five classes of a coding end, in the order they are asked
/// (`CodingEndClass::IN_ORDER`).
pub const CODING_END_PAIRS: [(CodingEndClass, RowPair); 5] = [
    (
        CodingEndClass::StructuralDeletionStopLost,
        structural_deletion_stop_lost_without_stop_codon_pair,
    ),
    (
        CodingEndClass::StopCodonCompleted,
        stop_codon_completed_by_polyadenylation_pair,
    ),
    (
        CodingEndClass::StopCodonOverlap,
        stop_codon_overlap_without_stop_codon_pair,
    ),
    (
        CodingEndClass::StopCodonBehindLastResidue,
        stop_codon_behind_last_residue_without_stop_codon_pair,
    ),
    (
        CodingEndClass::StartCodonOverlap,
        start_codon_overlap_without_start_codon_pair,
    ),
];

/// Whether a rule that stood before the coding ends were read may take the pair of the
/// reference entry `e` and the candidate entry `a`, beside its own reading of the two
/// entries. In a comparison that was given no cache it stands as it stood. Beside a cache
/// it takes a pair only where the transcript's record shows the codon the rule stands on
/// (`coding_ends::terms::landed_rule_may_take`), and nothing where the cache answers the
/// row with no one record.
pub fn stands_beside_a_cache(
    rule: LandedRule,
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> bool {
    let ends = match at.ends {
        CodingEndsRead::NoCache => return true,
        CodingEndsRead::NoRecord => None,
        CodingEndsRead::Record(ends) => Some(ends),
    };
    with_the_pair(e, a, at, |row, _, read| {
        landed_rule_may_take(rule, row, ends, read)
    })
    .unwrap_or(false)
}

/// Whether a rule that asks the record alone, and no fact of the row, stands on the
/// transcript of a pair: `start_cooccurrence_swap`, which stands where the record shows a
/// start codon at the coding start, and
/// `hgvsp_stop_named_on_incomplete_coding_sequence`, which stands but where
/// polyadenylation completes the stop codon. In a comparison that was given no cache each
/// stands as it stood (`coding_ends::terms::landed_rule_may_take`).
pub fn stands_on_the_record(rule: LandedRule, at: &AtACodingEnd<'_>) -> bool {
    let ends = match at.ends {
        CodingEndsRead::NoCache => return true,
        CodingEndsRead::NoRecord => None,
        CodingEndsRead::Record(ends) => Some(ends),
    };
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
    landed_rule_may_take(rule, &row, ends, &read)
}

/// The id of the standard genetic code, with which Ensembl VEP translates the alternate
/// of a frameshift whatever the transcript's own code.
const STANDARD_CODE: i64 = 1;

/// The description of a protein notation, after its `:p.` and with the `=` the text
/// formats escape written out, and the text before it; `None` for a value without one.
fn protein_description(value: &str) -> Option<(&str, String)> {
    let (protein, description) = value.split_once(":p.")?;
    Some((protein, description.replace("%3D", "=")))
}

/// Whether `ev` and `av` are the protein notations of a frameshift that the coding end
/// cuts short: Ensembl VEP's the deletion with its extension of unknown length, and
/// vep-rs's the frameshift named at the deletion's first residue, the one
/// `CutShortFrameshift::vep_rs_notation` gives, behind the same protein, or no notation
/// where that gives none.
fn cut_short_notations(
    cut_short: &crate::coding_ends::terms::CutShortFrameshift,
    ev: &str,
    av: &str,
) -> bool {
    let Some((protein, theirs)) = protein_description(ev) else {
        return false;
    };
    if theirs != cut_short.ensembl_vep_notation() {
        return false;
    }
    match cut_short.vep_rs_notation() {
        Some(owed) => {
            protein_description(av).is_some_and(|(named, ours)| named == protein && ours == owed)
        }
        None => blank(av),
    }
}

/// Whether the differing values `ev` and `av` of the field `name` on the reference entry
/// `e` and the candidate entry `a`, a pair the class `class` of a coding end takes, are
/// the ones the pair gives:
///
/// - the IMPACT of each side's own terms;
/// - a Downstream, NMD or LOFTEE field present exactly on the side whose terms let the
///   plugin run;
/// - under `stop_codon_overlap_without_stop_codon`, on a transcript at whose coding end
///   no stop codon stands, the two protein notations of a frameshift that the coding end
///   cuts short (`coding_ends::terms::frameshift_the_coding_end_cuts_short`): the
///   reference's deletion with an extension of unknown length against the candidate's
///   frameshift named at the first of the deleted residues;
/// - on a transcript whose stop codon polyadenylation completes, under
///   `stop_codon_completed_by_polyadenylation` and `stop_codon_overlap_without_stop_codon`,
///   the candidate's `Amino_acids` and `Codons` as they are derived on the completed coding
///   sequence, and the candidate's protein notation as it is derived there, behind the
///   protein the reference names where it writes a notation;
/// - under `start_codon_overlap_without_start_codon`, the reference's lost start at
///   residue 1 (`p.<Ref>1?`, TranscriptVariationAllele.pm:2161) against the candidate's
///   notation of the residues the cell both rows print holds, behind the same protein;
///   and for a frameshift the reference's `p.<Ref><number>?` at the residue its
///   frameshift walk stops at (:2161-2165) against the candidate's frameshift notation
///   that opens with the same residue and number, under the standard genetic code alone;
///   where that residue is not the one the row's translation span ends on and the row
///   carries no HGVS_OFFSET, the reference's `p.<Ref><number>_?<end>` instead
///   (`lost_start_range_of_a_frameshift`). A deletion from before the coding start whose
///   notation is made on the deletion moved 3' into the coding sequence, a row that
///   carries an HGVS_OFFSET and no frameshift term, is read as such a frameshift on the
///   moved span.
///
/// Any other cell is compared as a cell.
pub fn coding_end_following_cell(
    class: CodingEndClass,
    e: &Entry,
    a: &Entry,
    name: &str,
    ev: &str,
    av: &str,
    at: &AtACodingEnd<'_>,
) -> bool {
    let CodingEndsRead::Record(ends) = at.ends else {
        return false;
    };
    if ev == av {
        return false;
    }
    let (expected_set, actual_set) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    let expected_terms: Vec<&str> = expected_set.split(',').collect();
    let actual_terms: Vec<&str> = actual_set.split(',').collect();
    if field_named(name, IMPACT_FIELDS) {
        return terms_impact(&expected_terms).as_deref() == Some(ev)
            && terms_impact(&actual_terms).as_deref() == Some(av);
    }
    if let (Some(expected_written), Some(actual_written)) = (
        terms_write_field(name, &expected_terms),
        terms_write_field(name, &actual_terms),
    ) {
        return expected_written != actual_written
            && blank(ev) != expected_written
            && blank(av) != actual_written;
    }
    let completed = matches!(
        class,
        CodingEndClass::StopCodonCompleted | CodingEndClass::StopCodonOverlap
    );
    let derived =
        |cell: fn(
            &RowFacts<'_>,
            Option<&CodingEnds>,
            &WhatWasRead<'_>,
        )
            -> std::result::Result<String, crate::coding_ends::terms::NotDerived>| {
            with_the_pair(e, a, at, |row, _, read| cell(row, Some(ends), read).ok()).flatten()
        };
    if completed && field_named(name, AMINO_ACIDS_FIELDS) {
        return derived(amino_acids_on_a_completed_stop_codon).as_deref() == Some(av);
    }
    if completed && field_named(name, CODONS_FIELDS) {
        return derived(codons_on_a_completed_stop_codon).as_deref() == Some(av);
    }
    if !field_named(name, HGVSP_FIELDS) {
        return false;
    }
    if class == CodingEndClass::StopCodonOverlap {
        let cut_short = with_the_pair(e, a, at, |row, ours, read| {
            crate::coding_ends::terms::frameshift_the_coding_end_cuts_short(
                row,
                ours,
                Some(ends),
                read,
            )
            .ok()
        })
        .flatten();
        if let Some(cut_short) = cut_short {
            return cut_short_notations(&cut_short, ev, av);
        }
    }
    let Some((protein, candidate)) = protein_description(av) else {
        return false;
    };
    let reference = protein_description(ev);
    if reference
        .as_ref()
        .is_some_and(|(named, _)| *named != protein)
        || (reference.is_none() && !blank(ev))
    {
        return false;
    }
    if completed {
        return derived(protein_notation_on_a_completed_stop_codon)
            .is_some_and(|owed| owed.strip_prefix("p.") == Some(candidate.as_str()));
    }
    if class != CodingEndClass::StartCodonOverlap {
        return false;
    }
    let offsets = (
        field_of(Some(e), HGVS_OFFSET_FIELDS),
        field_of(Some(a), HGVS_OFFSET_FIELDS),
    );
    let Some(owed) = with_the_pair(e, a, at, |row, _, read| {
        let offsets = (offset_cell(&offsets.0), offset_cell(&offsets.1));
        protein_notation_where_no_start_codon_is_shown(row, Some(ends), read, offsets).ok()
    })
    .flatten() else {
        return false;
    };
    let Some(owed) = owed.strip_prefix("p.") else {
        return false;
    };
    // A notation made on an allele moved 3' is a frameshift's, read on the moved span,
    // whatever the row's terms.
    let unshifted = blank(&offsets.0);
    // Ensembl VEP translates a frameshift's alternate with the standard code and vep-rs
    // with the transcript's: under another code the two walks are not shown to stop at
    // one residue, so the notation is not read as following there.
    let another_code = ends
        .unshown_start
        .as_deref()
        .is_some_and(|kept| kept.table != STANDARD_CODE);
    if another_code && (actual_terms.contains(&"frameshift_variant") || !unshifted) {
        return false;
    }
    // Ensembl VEP's side names the residue and the number vep-rs's notation opens with,
    // then `?`: the lost start at residue 1 of a substitution and, for a frameshift, the
    // residue its frameshift walk stops at, whose alternate is set to `?` after the
    // frameshift peptides are taken (TranscriptVariationAllele.pm:2161-2165).
    let digits = owed.bytes().skip(3).take_while(u8::is_ascii_digit).count();
    let lost_start = owed
        .get(..3 + digits)
        .filter(|_| digits > 0)
        .map(|opening| format!("{opening}?"));
    let theirs = reference.map(|(_, description)| description);
    // A frameshift walk that stops at another residue than the end of the row's
    // translation span leaves Ensembl VEP's notation both numbers, which it writes as a
    // range around the `?`. The span is the one the row's cells write only where the
    // notation is not shifted 3'.
    let range = (unshifted && actual_terms.contains(&"frameshift_variant"))
        .then(|| {
            with_the_pair(e, a, at, |row, _, read| {
                lost_start_range_of_a_frameshift(owed, row, ends, read)
            })
            .flatten()
        })
        .flatten();
    owed == candidate && (theirs == lost_start || (range.is_some() && theirs == range))
}

/// An `HGVS_OFFSET` value as a cell: empty for an entry that carries none.
fn offset_cell(value: &str) -> &str {
    if blank(value) {
        ""
    } else {
        value
    }
}

/// `hgvsp_frameshift_count_at_the_first_residue` on a differing HGVSp. On a transcript
/// at whose coding start no start codon stands, flagged `cds_start_NF` or not, vep-rs
/// counts the residues of a frameshift named at residue 1 to the first stop codon of the
/// new reading frame; Ensembl VEP writes the unknown count for every frameshift named
/// at residue 1 (TranscriptVariationAllele.pm:2483). The class takes the field where
/// the two values are one frameshift named at residue 1 behind one protein, the
/// reference's with the unknown count and the candidate's with a number
/// (`count_alone_differs_at_the_first_residue`), and both are the texts derived for the
/// row from the reading of its transcript's record in the candidate's cache (`models`),
/// the row's allele, `CDS_position` and `HGVS_OFFSET` and the length of the record's
/// reference bases (`coding_ends::terms::frameshift_counted_at_the_first_residue`). The
/// cache is read for a pair of the shape alone. A comparison that was given no cache,
/// and a row the cache answers with no one record, are not the class's.
pub fn frameshift_count_at_the_first_residue(
    e: &Entry,
    a: &Entry,
    ev: &str,
    av: &str,
    models: Option<&PartialCodonModels>,
) -> Result<bool> {
    if !count_alone_differs_at_the_first_residue(ev, av) {
        return Ok(false);
    }
    let offsets = (
        field_of(Some(e), HGVS_OFFSET_FIELDS),
        field_of(Some(a), HGVS_OFFSET_FIELDS),
    );
    asked_at(models, None, e, |at| {
        let CodingEndsRead::Record(ends) = at.ends else {
            return Ok(false);
        };
        let taken = with_the_pair(e, a, at, |row, ours, read| {
            let offsets = (offset_cell(&offsets.0), offset_cell(&offsets.1));
            frameshift_counted_at_the_first_residue(row, ours, Some(ends), read, offsets)
                .is_ok_and(|counted| counted.takes(ev, av))
        });
        Ok(taken.unwrap_or(false))
    })
}

/// The most severe term among the consequence entries of each of `records` in `entries`,
/// the record entries themselves set aside, as the JSON format reads it
/// (`self_parity::most_severe`, OutputFactory/JSON.pm:387-388).
pub fn most_severe_terms(entries: &[Entry], records: &HashSet<i64>) -> HashMap<i64, String> {
    let mut terms: HashMap<i64, BTreeSet<String>> = HashMap::new();
    for entry in entries {
        let Some(record) = entry.record.filter(|record| records.contains(record)) else {
            continue;
        };
        if entry.key.ends_with("|record") {
            continue;
        }
        terms.entry(record).or_default().extend(
            entry
                .consequence
                .split(['&', ','])
                .filter(|term| !term.is_empty())
                .map(str::to_string),
        );
    }
    terms
        .into_iter()
        .map(|(record, terms)| (record, most_severe(&terms)))
        .collect()
}
