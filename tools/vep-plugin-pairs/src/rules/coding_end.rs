// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A pair of rows at a coding end of a transcript, read beside the candidate cache: the
//! class of a coding end that takes the pair's two Consequence sets, and whether a class of
//! a stop codon or a start codon that stands may take them on the transcript. The rules
//! are those of `coding_ends::terms`; this file hands them a row as this tool reads one.

use super::coding_ends::terms::{
    class_of_pair, landed_rule_may_take, Cells, CodingEndClass, LandedRule, RowFacts, WhatWasRead,
};
use super::coding_ends::CodingEnds;
use super::entry::{consequence_set, entry_kind, feature_kind, Entry};

/// The lowest and the highest genomic base of a Location, `None` for a text that names no
/// position after its sequence name.
pub fn location_bases(location: &str) -> Option<(i64, i64)> {
    let (_, position) = location.rsplit_once(':')?;
    let (first, last) = match position.split_once('-') {
        Some((first, last)) => (first.parse::<i64>().ok()?, last.parse::<i64>().ok()?),
        None => {
            let base = position.parse::<i64>().ok()?;
            (base, base)
        }
    };
    Some((first.min(last), first.max(last)))
}

/// The Location a row's allele lies at: the row's own, or for a row without one the
/// Location its record gives its allele.
pub fn row_location(e: &Entry) -> Option<&str> {
    e.location
        .as_deref()
        .filter(|location| !location.is_empty())
        .or(e.allele_location.as_deref())
}

/// The five cells of a default-format row that the derivation of a row's coding terms
/// reads; `None` for a row that lacks one of the five columns.
fn cells(e: &Entry) -> Option<Cells<'_>> {
    Some(Cells {
        cdna_position: e.fields.get("cDNA_position")?,
        cds_position: e.fields.get("CDS_position")?,
        protein_position: e.fields.get("Protein_position")?,
        amino_acids: e.fields.get("Amino_acids")?,
        codons: e.fields.get("Codons")?,
    })
}

/// Hands `ask` the pair of the reference's row `e` and the candidate's `a` as the rules of
/// a coding end read one, with the candidate's terms; `None` for a row whose Location
/// names no position.
///
/// The count of reference bases the allele stands in for is that of the reference bases
/// either row carries (`Entry.reference`); without them a deletion's is the count of the
/// positions its Location names, every one of which the allele `-` removes, and any other
/// allele's is not known.
fn asked<T>(
    e: &Entry,
    a: &Entry,
    ask: impl FnOnce(&RowFacts<'_>, &str, &WhatWasRead<'_>) -> T,
) -> Option<T> {
    let location = location_bases(row_location(e)?)?;
    let (theirs, ours) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    let reference = e
        .reference
        .borrow()
        .clone()
        .or_else(|| a.reference.borrow().clone());
    let removed = (i128::from(location.1) - i128::from(location.0)).saturating_add(1);
    let reference_length = match reference {
        Some(reference) => Some(reference.chars().count() as u64),
        None if e.allele == "-" => u64::try_from(removed).ok(),
        None => None,
    };
    let row = RowFacts {
        transcript_row: feature_kind(entry_kind(e)) == "transcript",
        allele: &e.allele,
        location,
        ensembl_vep_terms: &theirs,
    };
    let read = WhatWasRead {
        cells: cells(e).zip(cells(a)),
        reference_length,
        reference_bases: None,
        record_on_the_reverse_strand: false,
    };
    Some(ask(&row, &ours, &read))
}

/// The class of a coding end that takes the two Consequence sets of the reference's row
/// `e` and the candidate's `a` on a transcript whose coding ends are `ends`
/// (`coding_ends::terms::class_of_pair`); `None` where none takes them, and without a
/// reading of the transcript's coding ends.
pub fn coding_end_class(e: &Entry, a: &Entry, ends: Option<&CodingEnds>) -> Option<CodingEndClass> {
    asked(e, a, |row, ours, read| class_of_pair(row, ours, ends, read)).flatten()
}

/// Whether a comparison that was given the cache may let the class `rule` of a stop codon
/// or a start codon that stands take the pair of `e` and `a` on a transcript whose coding
/// ends are `ends`, beside the rule's own reading of the two rows
/// (`coding_ends::terms::landed_rule_may_take`).
pub fn landed_rule_takes_beside_a_cache(
    rule: LandedRule,
    e: &Entry,
    a: &Entry,
    ends: Option<&CodingEnds>,
) -> bool {
    asked(e, a, |row, _, read| {
        landed_rule_may_take(rule, row, ends, read)
    })
    .unwrap_or(false)
}
