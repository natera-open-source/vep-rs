// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The coding ends of one transcript object of a JSON cache: the facts of the record that
//! `coding_ends` asks, each read as vep-rs loads the record.

use super::transcript::TranscriptModel;
use crate::pytext::{self, Json};
use crate::rules::coding_ends::{coding_ends, CodingEnds, MapperFacts, RecordFacts};
use crate::rules::shared::cache_integer;

/// An integer of a cache transcript that 64 bits hold.
fn integer(value: Option<&Json>) -> Option<i64> {
    cache_integer(value).and_then(|number| i64::try_from(number).ok())
}

/// A member that is text, `Some(None)` for one the record lacks or holds as `null`; `None`
/// for any other value, which vep-rs refuses.
fn text(member: Option<&Json>) -> Option<Option<&str>> {
    match member {
        None | Some(Json::Null) => Some(None),
        Some(Json::Str(text)) => Some(Some(text)),
        Some(_) => None,
    }
}

/// The id of the genetic code a `codon_table` member names: the standard code's for a
/// member the record lacks or holds as `null`, else the whole number the member holds as
/// a number or as text; `None` for any other value.
fn codon_table(member: Option<&Json>) -> Option<i64> {
    let digits = match member {
        None | Some(Json::Null) => return Some(1),
        Some(Json::Int(digits)) => digits.as_str(),
        Some(Json::Str(text)) => text.as_str(),
        Some(_) => return None,
    };
    if !pytext::isdecimal(digits) {
        return None;
    }
    pytext::int(digits).and_then(|table| i64::try_from(table).ok())
}

/// The record's `start` and `end`, `None` where it states either in no form an integer is
/// read from.
fn record_span(transcript: &Json) -> Option<(i64, i64)> {
    Some((
        integer(transcript.get("start"))?,
        integer(transcript.get("end"))?,
    ))
}

/// The facts of a record whose coding sequence is `coding`, read into its coding ends;
/// `None` for a record that states one of them in a form vep-rs reads otherwise or
/// refuses: `attributes` that are no list, an attribute whose `code` or `value` is no
/// text, a `ccds` or a `peptide` that is neither text nor `null`, a `codon_table` that
/// names no whole number, a `start_phase` or a coding bound of the mapper that is no small
/// whole number, a mapper beside which the record states no `TranscriptModel` of numbers
/// 64 bits hold, or `introns` that are no list of two whole numbers each. The introns are
/// the cached ones, or where the record caches none the stretches between its exons. The
/// record's 3' UTR is not handed to the reading: this comparison derives no protein
/// notation, so no transcript's sequence is kept for one.
fn read_coding_ends(transcript: &Json, cached: &Json, coding: &str) -> Option<CodingEnds> {
    let attribute_codes: Vec<&str> = match transcript.get("attributes") {
        None => Vec::new(),
        Some(Json::List(attributes)) => attributes
            .iter()
            .map(
                |attribute| match (attribute.get("value"), attribute.get("code")) {
                    (Some(Json::Str(_)), Some(Json::Str(code))) => Some(code.as_str()),
                    _ => None,
                },
            )
            .collect::<Option<_>>()?,
        Some(_) => return None,
    };
    let codon_table = codon_table(cached.get("codon_table"))?;
    let mut introns: Vec<(i64, i64)> = match cached.get("introns") {
        None | Some(Json::Null) => Vec::new(),
        Some(Json::List(introns)) => introns
            .iter()
            .map(|intron| Some((integer(intron.get("start"))?, integer(intron.get("end"))?)))
            .collect::<Option<_>>()?,
        Some(_) => return None,
    };
    let pairs: Vec<(i64, i64, i64, i64, i64)>;
    let mapper = match cached.get("mapper") {
        None | Some(Json::Null) => None,
        Some(mapper @ Json::Object(_)) => {
            let bound = |name: &str| match mapper.get(name) {
                None | Some(Json::Null) => Some(0),
                bound => integer(bound).filter(|bound| *bound >= 0),
            };
            let start_phase = match mapper.get("start_phase") {
                None => 0,
                phase => integer(phase).filter(|phase| i8::try_from(*phase).is_ok())?,
            };
            let model = TranscriptModel::from_cache(transcript, i128::from(start_phase.max(0)))?;
            let narrow = |value: i128| i64::try_from(value).ok();
            narrow(model.coding_start)?;
            narrow(model.coding_end)?;
            pairs = model
                .pairs
                .iter()
                .map(|pair| {
                    Some((
                        narrow(pair[0])?,
                        narrow(pair[1])?,
                        narrow(pair[2])?,
                        narrow(pair[3])?,
                        narrow(pair[4])?,
                    ))
                })
                .collect::<Option<_>>()?;
            let mut exons: Vec<(i64, i64)> = model
                .exons
                .iter()
                .map(|(start, end)| Some((narrow(*start)?, narrow(*end)?)))
                .collect::<Option<_>>()?;
            if introns.is_empty() {
                exons.sort_unstable();
                introns = exons
                    .windows(2)
                    .filter_map(|exons| {
                        Some((exons[0].1.checked_add(1)?, exons[1].0.checked_sub(1)?))
                    })
                    .collect();
            }
            Some(MapperFacts {
                start_phase,
                coding_start: bound("cdna_coding_start")?,
                coding_end: bound("cdna_coding_end")?,
                pairs: &pairs,
            })
        }
        Some(_) => return None,
    };
    Some(coding_ends(&RecordFacts {
        attribute_codes: &attribute_codes,
        consensus_id: text(transcript.get("ccds"))?,
        strand: integer(transcript.get("strand")).filter(|strand| matches!(strand, 1 | -1)),
        codon_table: Some(codon_table),
        coding: Some(coding),
        peptide: text(cached.get("peptide"))?,
        mapper,
        span: record_span(transcript),
        introns: &introns,
    }))
}

/// The stable id of one cache transcript object and its coding ends, read as vep-rs loads
/// the record; `None` for an object without a stable id that is text. The mapper pairs
/// are the transcript's `TranscriptModel`'s. A record without a coding sequence that is
/// text is not read, nor is one whose facts `read_coding_ends` refuses; a record that is
/// not read keeps its span.
pub fn coding_ends_of(transcript: &Json) -> Option<(&str, CodingEnds)> {
    let Some(Json::Str(stable_id)) = transcript.get("stable_id") else {
        return None;
    };
    let cached = match transcript.get("variation_effect_feature_cache") {
        Some(cached @ Json::Object(_)) => Some(cached),
        _ => None,
    };
    let coding = match cached.and_then(|cached| cached.get("translateable_seq")) {
        Some(Json::Str(coding)) => Some(coding.as_str()),
        _ => None,
    };
    let unread = |tail_length: u8| CodingEnds {
        span: record_span(transcript),
        ..CodingEnds::unread(tail_length)
    };
    let ends = match (cached, coding) {
        (Some(cached), Some(coding)) => read_coding_ends(transcript, cached, coding)
            .unwrap_or_else(|| unread((coding.chars().count() % 3) as u8)),
        _ => unread(0),
    };
    Some((stable_id, ends))
}
