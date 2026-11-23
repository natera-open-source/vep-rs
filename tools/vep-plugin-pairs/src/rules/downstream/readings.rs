// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a row says of its coding span and of its place among its record's alleles: the CDS
//! pair it prints, the genomic span of its Location, the pair a plugin reads once the HGVS
//! notation has moved the span, and the offsets the rows of the record's earlier alleles
//! carry.

use super::transcript::TranscriptModel;
use crate::pytext;
use crate::rules::entry::{field_of, uploaded_alleles, Entry, Row, ABSENT};
use crate::rules::shared::{
    digits_value, integer_of, is_non_zero_integer, record_entry, sibling_rows,
    vcf_alleles_as_spelt, HGVS_OFFSET_FIELDS,
};
use crate::rules::Res;

/// A `cds_start` and a `cds_end`, `None` for an end that is no coding base.
pub type CdsPair = (Option<i128>, Option<i128>);

/// One reading of a row the notation moved: its own CDS pair, the pair of its span moved
/// by the shift where both ends of that span are coding bases, and the record's genomic
/// end.
pub type ShiftedReading = (CdsPair, Option<(i128, i128)>, i128);

/// Whether an allele is a sequence of nucleotide letters in either case.
fn sequence_allele(allele: &str) -> bool {
    !allele.is_empty() && allele.bytes().all(|b| b"ACGTNacgtn".contains(&b))
}

/// The `cds_start` and `cds_end` an entry's coding span may spell, a gap (`?`, or a bound
/// the entry lacks) `None`, and no pair for an entry without a coding coordinate: the text
/// formats' `CDS_position` or the keys `cds_start` and `cds_end` as given and, for a
/// sequence allele at two adjacent numbers, swapped as well, the pair of an insertion,
/// whose start is one past its end and which every format spells lower first
/// (`Utils::format_coords`). Each bound is read up to the `/total` suffix `--total_length`
/// appends.
pub fn cds_pair_readings(e: &Entry) -> Vec<CdsPair> {
    let (start, end) = match e.fields.get("CDS_position") {
        Some(position) => {
            let (first, _, last) = pytext::partition(pytext::partition(position, "/").0, "-");
            (first, if last.is_empty() { first } else { last })
        }
        None => (
            pytext::partition(e.fields.get_or("cds_start", ABSENT), "/").0,
            pytext::partition(e.fields.get_or("cds_end", ABSENT), "/").0,
        ),
    };
    let (start, end) = (digits_value(start), digits_value(end));
    if start.is_none() && end.is_none() {
        return Vec::new();
    }
    let mut readings = vec![(start, end)];
    if let (Some(first), Some(last)) = (start, end) {
        if first.checked_add(1) == Some(last) && sequence_allele(&e.allele) {
            readings.push((end, start));
        }
    }
    readings
}

/// The genomic start and end a text entry's Location may spell: the span as printed and,
/// for a sequence allele at two adjacent bases, the point between them as well, an
/// insertion's span, whose start is one past its end.
pub fn record_span_readings(e: &Entry) -> Vec<(i128, i128)> {
    let location = e.location.as_deref().unwrap_or("");
    let (first, _, last) = pytext::partition(pytext::rpartition(location, ":").2, "-");
    let last = if last.is_empty() { first } else { last };
    let (Some(start), Some(end)) = (digits_value(first), digits_value(last)) else {
        return Vec::new();
    };
    let mut readings = vec![(start, end)];
    if start.checked_add(1) == Some(end) && sequence_allele(&e.allele) {
        readings.push((end, start));
    }
    readings
}

/// The readings of a row whose reference entry `e` carries the HGVS_OFFSET `shift`, each
/// the own CDS pair (one of `own_pairs`), the pair a plugin reads under --hgvs where the
/// shift reaches it (the record's span moved by `shift` has a CDS coordinate at both ends),
/// else `None`, and the record's genomic end. The record's span is a reading of the
/// Location that maps to an own pair or, for an entry without a Location, the one span an
/// own pair maps back to (`TranscriptModel::genomic_span`).
pub fn shifted_pair_readings(
    e: &Entry,
    model: &TranscriptModel,
    shift: i128,
    own_pairs: &[CdsPair],
) -> Vec<ShiftedReading> {
    let located = e
        .location
        .as_deref()
        .is_some_and(|location| !location.is_empty());
    let spans: Vec<((i128, i128), CdsPair)> = if located {
        record_span_readings(e)
            .into_iter()
            .map(|span| (span, model.cds_pair(span.0, span.1)))
            .filter(|(_, own)| own_pairs.contains(own))
            .collect()
    } else {
        own_pairs
            .iter()
            .filter_map(|own| match own {
                (Some(start), Some(end)) => Some((model.genomic_span(*start, *end)?, *own)),
                _ => None,
            })
            .collect()
    };
    spans
        .into_iter()
        .map(|((start, end), own)| {
            // An insertion's point stays one when a flank moves past the ends of `i128`; any
            // other span with an end there has no coding base at that end.
            let insertion = end.checked_add(1) == Some(start);
            let moved = match (start.checked_add(shift), end.checked_add(shift)) {
                (Some(start), Some(end)) => model.cds_pair(start, end),
                (start, end) if insertion => model.insertion_pair(start, end),
                _ => (None, None),
            };
            let shifted = match moved {
                (Some(first), Some(last)) => Some((first, last)),
                _ => None,
            };
            (own, shifted, end)
        })
        .collect()
}

/// The entry's HGVS_OFFSET where it is a non-zero integer, `None` for any other value.
pub fn non_zero_offset(e: &Entry) -> Res<Option<i128>> {
    let offset = field_of(e, &HGVS_OFFSET_FIELDS);
    if !is_non_zero_integer(offset) {
        return Ok(None);
    }
    integer_of(offset, "HGVS_OFFSET").map(Some)
}

/// The alternate alleles of `e`'s record in the record's own order, spelt as its rows spell
/// the allele, where an entry states them: a row's identifier-less `Uploaded_variation`,
/// else the `ALT` or the `allele_string` of the record entry among `record_rows`; `None`
/// where none does.
pub fn alternates_in_order(e: &Entry, record_rows: Option<&[Row]>) -> Option<Vec<String>> {
    if let Some((reference, alternates)) =
        uploaded_alleles(e.fields.get_or("Uploaded_variation", ""))
    {
        return Some(vcf_alleles_as_spelt(&reference, &alternates).1);
    }
    let record = record_entry(e, record_rows.unwrap_or(&[]))?;
    if let Some(alternates) = record.fields.get("ALT") {
        let alternates: Vec<&str> = alternates.split(',').collect();
        return Some(vcf_alleles_as_spelt(record.fields.get_or("REF", ""), &alternates).1);
    }
    let alternates: Vec<String> = record
        .fields
        .get_or("allele_string", "")
        .split('/')
        .skip(1)
        .map(str::to_string)
        .collect();
    (!alternates.is_empty()).then_some(alternates)
}

/// The alternates of the record of a second or later alternate allele's row `e`: an entry
/// states the record's alternates in order (`alternates_in_order`), one of them alone
/// spells the row's allele, and it is not the first. `None` for any other row.
fn later_alternates(e: &Entry, record_rows: Option<&[Row]>) -> Option<Vec<String>> {
    let alternates = alternates_in_order(e, record_rows)?;
    let named = alternates
        .iter()
        .filter(|allele| **allele == e.allele)
        .count();
    (named == 1 && alternates[0] != e.allele).then_some(alternates)
}

/// Whether `e` is the row of a second or later alternate allele of its record.
pub fn later_alternate(e: &Entry, record_rows: Option<&[Row]>) -> bool {
    later_alternates(e, record_rows).is_some()
}

/// The HGVS_OFFSET of the reference's row of the first alternate allele of `e`'s record on
/// `e`'s feature, for the row `e` of a later alternate (`later_alternate`); `None` for any
/// other row, without that one row among `record_rows` and where it carries no non-zero
/// offset.
pub fn first_allele_shift(e: &Entry, record_rows: Option<&[Row]>) -> Res<Option<i128>> {
    let Some(alternates) = later_alternates(e, record_rows) else {
        return Ok(None);
    };
    let siblings = sibling_rows(e, record_rows.unwrap_or(&[]));
    let mut rows = siblings.iter().filter(|r| r.allele == alternates[0]);
    match (rows.next(), rows.next()) {
        (Some(row), None) => non_zero_offset(row),
        _ => Ok(None),
    }
}

/// The non-zero HGVS_OFFSET of the reference's row, on `e`'s feature, of the last alternate
/// allele before `e`'s own in its record's order whose row carries one, for the row `e` of
/// a later alternate (`later_alternate`): the shift whose CDS pair that allele's protein
/// notation left on the TranscriptVariation the record's alleles share, where each allele's
/// own notation moves it afresh (TranscriptVariationAllele.pm:844-858). `None` for any
/// other row, where no earlier allele's row carries an offset, and where an earlier allele
/// has two rows on the feature among `record_rows`.
pub fn earlier_allele_shift(e: &Entry, record_rows: Option<&[Row]>) -> Res<Option<i128>> {
    let Some(alternates) = later_alternates(e, record_rows) else {
        return Ok(None);
    };
    let siblings = sibling_rows(e, record_rows.unwrap_or(&[]));
    let own = alternates
        .iter()
        .position(|allele| *allele == e.allele)
        .unwrap_or(0);
    for allele in alternates[..own].iter().rev() {
        let mut rows = siblings.iter().filter(|r| r.allele == *allele);
        let first = rows.next();
        if first.is_some() && rows.next().is_some() {
            return Ok(None);
        }
        if let Some(offset) = first.map(|row| non_zero_offset(row)).transpose()?.flatten() {
            return Ok(Some(offset));
        }
    }
    Ok(None)
}

/// The moves a later allele's row may print its pair under: the row's own offset and the
/// earlier allele's, each once and neither zero.
pub fn moves(shift: i128, earlier: Option<i128>) -> Vec<i128> {
    let mut moves = Vec::new();
    for by in [Some(shift), earlier].into_iter().flatten() {
        if by != 0 && !moves.contains(&by) {
            moves.push(by);
        }
    }
    moves
}

/// Whether the reference entry `e` prints the record's span moved by `by` bases
/// (`shifted_pair_readings`), the span that of a CDS pair the candidate entry `a` prints,
/// or of its pair `own` alone.
pub fn printed_pair_is_moved(
    e: &Entry,
    a: &Entry,
    model: &TranscriptModel,
    by: i128,
    own: Option<CdsPair>,
) -> bool {
    let printed = cds_pair_readings(e);
    let own_pairs = match own {
        Some(own) => vec![own],
        None => cds_pair_readings(a),
    };
    shifted_pair_readings(e, model, by, &own_pairs)
        .into_iter()
        .any(|(_, moved, _)| {
            moved.is_some_and(|(first, last)| printed.contains(&(Some(first), Some(last))))
        })
}

/// Whether the two CDS pairs of a row are those of a later allele under --hgvs: `e` is the
/// reference's row of a second or later alternate allele of its record (`later_alternate`),
/// the candidate entry `a` prints the record's own pair, and the reference prints the
/// record's span moved by its HGVS_OFFSET `shift` or by the offset of the last earlier
/// allele's row that carries one (`earlier_allele_shift`): the pair an earlier allele's
/// protein notation left on the TranscriptVariation the record's alleles share
/// (OutputFactory.pm:1746 over BaseTranscriptVariation.pm:252-267), which is the row's own
/// offset only where the two alleles move alike.
pub fn later_allele_prints_shifted_pair(
    e: &Entry,
    a: &Entry,
    model: Option<&TranscriptModel>,
    shift: i128,
    record_rows: Option<&[Row]>,
) -> Res<bool> {
    let Some(model) = model else {
        return Ok(false);
    };
    if !later_alternate(e, record_rows) {
        return Ok(false);
    }
    let earlier = earlier_allele_shift(e, record_rows)?;
    Ok(moves(shift, earlier)
        .into_iter()
        .any(|by| printed_pair_is_moved(e, a, model, by, None)))
}
