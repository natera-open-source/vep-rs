// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The classes of the HGVSIntronOffset plugin's two fields, `HGVS_IntronStartOffset` and
//! `HGVS_IntronEndOffset`.
//!
//! HGVSIntronOffset.pm:62-72 ask every Transcript row for its HGVS notation and write the
//! intron parts of the notation's two coordinates, so the two fields differ wherever the two
//! engines describe a row's change by different notations. Each class is one such pair of
//! notations: its rule reads the row's HGVSc, HGVS_OFFSET, cDNA position, Location and allele
//! on both sides and decides the plugin's fields beside them, and the HGVSc and the
//! HGVS_OFFSET of the same row.
//!
//! A row is read from its own fields alone, but for `insertion_left_where_the_shift_stops`,
//! which reads the exon spans of the row's transcript where the run names a JSON cache, and
//! for `added_bases_written_outside_the_transcript`, which reads the reference bases a row
//! carries from the run's input records where it carries them.

use super::entry::{blank, consequence_set, field_named, field_of, Entry, Row, ABSENT};
use super::shared::{
    at_its_location, indel_row, insertion_row, integer_of, is_integer, is_non_zero_integer,
    is_sequence, location_ends, positive_shift, re_digit, reverse_complement, CDNA_END_FIELD,
    CDNA_START_FIELDS, DISTANCE_FIELDS, HGVSC_FIELDS, HGVS_INTRON_END_OFFSET_FIELDS,
    HGVS_INTRON_OFFSET_FIELDS, HGVS_INTRON_START_OFFSET_FIELDS, HGVS_OFFSET_FIELDS, INDEL_NOTATION,
    SHIFT_LENGTH_FIELDS, STRAND_FIELDS,
};
use super::{Refusal, Res};
use crate::pytext;
use regex::Regex;
use std::sync::LazyLock;

const BAM_EDIT_FIELDS: [&str; 2] = ["BAM_EDIT", "bam_edit"];
/// The offset of an edited model's edits, which the reference stores inside the notation's
/// routine (TranscriptVariationAllele.pm:1517, :1621) and prints on the row
/// (OutputFactory.pm:1778).
const REFSEQ_OFFSET_FIELDS: [&str; 2] = ["REFSEQ_OFFSET", "refseq_offset"];
/// The terms of a row annotated outside its transcript.
const OUTSIDE_TRANSCRIPT_TERMS: [&str; 2] = ["upstream_gene_variant", "downstream_gene_variant"];
/// The fields `hgvs_shifted_indel_written_span_gate` asks the reference to leave empty.
const SHIFTED_INDEL_GATE_FIELDS: [&str; 6] = [
    "HGVSc",
    "hgvsc",
    "HGVS_OFFSET",
    "hgvs_offset",
    "HGVS_IntronStartOffset",
    "HGVS_IntronEndOffset",
];
/// The bases of flank `perform_shift` reads on either side of an allele
/// (TranscriptVariationAllele.pm:427-437).
const SHIFT_FLANK_BASES: i128 = 1000;

/// The largest size of a number the rules read from a row: a sum of the few numbers one
/// rule adds stays inside `i128` at or below it.
const LARGEST_NUMBER: u128 = 10u128.pow(37);

const LOCATION_POSITION: &str = "Location position";
const CDNA_POSITION: &str = "cDNA position";
const HGVS_COORDINATE: &str = "HGVS coordinate";
const HGVS_OFFSET: &str = "HGVS_OFFSET";

/// One HGVS cDNA coordinate as a pattern: an optional `*`, an optional `-`, digits and an
/// optional intron offset.
fn coordinate_pattern() -> String {
    let digit = re_digit();
    format!(r"\*?-?[{digit}]+(?:[+-][{digit}]+)?")
}

/// A `del` notation after the transcript's name: the numbering, the first coordinate and the
/// second where it names two. Anchored at both ends.
static DELETION_NOTATION: LazyLock<Regex> = LazyLock::new(|| {
    let coordinate = coordinate_pattern();
    Regex::new(&format!(r"^([cn])\.({coordinate})(?:_({coordinate}))?del$")).expect("regex")
});

/// The numbering, the first coordinate and the last of a `del` notation after the
/// transcript's name, one coordinate standing for both where it names one.
fn deletion_notation(body: &str) -> Option<(&str, &str, &str)> {
    let found = DELETION_NOTATION.captures(body)?;
    let group = |index| found.get(index).map(|m| m.as_str());
    let first = group(2)?;
    Some((group(1)?, first, group(3).unwrap_or(first)))
}

/// A notation after the transcript's name whose kind the rules read: a substitution, a
/// `delins` or an `ins` with its bases, an `inv`, a `del` or a `dup`. Anchored at both ends.
static CLIPPED_NOTATION: LazyLock<Regex> = LazyLock::new(|| {
    let coordinate = coordinate_pattern();
    Regex::new(&format!(
        r"^(?P<numbering>[cn])\.(?P<start>{coordinate})(?:_(?P<end>{coordinate}))?(?:[ACGTN]>(?P<alternate_base>[ACGT])|(?P<kind>delins|ins)(?P<bases>[ACGT]+)|(?P<plain>inv|del|dup))$"
    ))
    .expect("regex")
});

/// The exon coordinate and the intron offset `hgvs_transcript` reads off a cDNA coordinate
/// (TranscriptVariationAllele.pm:1530-1531): the first integer and the one after it, the star
/// of a 3' UTR coordinate unread.
static COORDINATE_PARTS_AS_READ: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(-?[0-9]+)\+?(-?[0-9]+)?").expect("regex"));

/// A `CLIPPED_NOTATION` match.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Notation<'a> {
    /// `c` or `n`.
    numbering: &'a str,
    start: &'a str,
    end: Option<&'a str>,
    /// Whether the notation is a substitution of one base.
    substitution: bool,
    /// `delins` or `ins`.
    kind: Option<&'a str>,
    /// The bases a `delins` or an `ins` writes; empty for every other notation.
    bases: &'a str,
    /// `inv`, `del` or `dup`.
    plain: Option<&'a str>,
}

impl<'a> Notation<'a> {
    /// The second coordinate, the first standing for both where the notation names one.
    fn last(&self) -> &'a str {
        self.end.unwrap_or(self.start)
    }
}

/// The stable id of a transcript name: the text before its first `.`.
fn stable_id(name: &str) -> &str {
    name.split('.').next().unwrap_or(name)
}

/// The notation of an HGVSc that names the transcript `feature` (its stable id, with or
/// without a version); `None` for an empty value, another transcript's notation and a
/// notation of another shape.
fn clipped_notation<'a>(hgvsc: &'a str, feature: &str) -> Option<Notation<'a>> {
    if blank(hgvsc) {
        return None;
    }
    let (name, _, body) = pytext::partition(hgvsc, ":");
    if stable_id(name) != stable_id(feature) {
        return None;
    }
    let found = CLIPPED_NOTATION.captures(body)?;
    let group = |name: &str| found.name(name).map(|m| m.as_str());
    Some(Notation {
        numbering: group("numbering")?,
        start: group("start")?,
        end: group("end"),
        substitution: group("alternate_base").is_some(),
        kind: group("kind"),
        bases: group("bases").unwrap_or(""),
        plain: group("plain"),
    })
}

fn is_digit(c: char) -> bool {
    pytext::decimal_value(c).is_some()
}

/// Whether `allele` is one or more of the four bases: a sequence without `N`.
fn is_bases(allele: &str) -> bool {
    !allele.is_empty() && allele.bytes().all(|b| b"ACGT".contains(&b))
}

/// The count of characters in `text`.
fn length(text: &str) -> i128 {
    text.chars().count() as i128
}

fn not_an_integer(digits: &str, what: &str) -> Refusal {
    let shown = if digits.chars().count() <= 24 {
        digits.to_string()
    } else {
        digits.chars().take(24).collect::<String>() + "..."
    };
    Refusal(format!("{what} {} is not an integer", pytext::repr(&shown)))
}

/// `integer_of` for a number a rule adds to others: one above `LARGEST_NUMBER` in size is
/// refused as a number `integer_of` does not read is.
fn whole(digits: &str, what: &str) -> Res<i128> {
    let value = integer_of(digits, what)?;
    if value.unsigned_abs() > LARGEST_NUMBER {
        return Err(not_an_integer(digits, what));
    }
    Ok(value)
}

/// The bases a Location spans from `start` to `end`, both read as numbers.
fn span_bases(start: &str, end: &str) -> Res<i128> {
    let last = whole(end, LOCATION_POSITION)?;
    Ok(last - whole(start, LOCATION_POSITION)? + 1)
}

/// The closed end of a row's cDNA position: the part after the dash of the text formats'
/// `cDNA_position`, read up to a `/total` suffix, else `cdna_end`; empty where the row
/// writes neither.
fn cdna_closed_end(e: &Entry) -> &str {
    let position = pytext::partition(e.fields.get_or(CDNA_START_FIELDS[0], ""), "/").0;
    let high = pytext::partition(position, "-").2;
    if high.is_empty() {
        e.fields.get_or(CDNA_END_FIELD, "")
    } else {
        high
    }
}

/// The values of the plugin's two fields on a row, the start's and the end's.
fn intron_offset_pair(e: &Entry) -> [&str; 2] {
    HGVS_INTRON_OFFSET_FIELDS.map(|field| field_of(e, &[field]))
}

fn all_blank(pairs: [&[&str; 2]; 2]) -> bool {
    pairs
        .iter()
        .all(|pair| pair.iter().all(|value| blank(value)))
}

/// Whether the entry's `SHIFT_LENGTH`, where the run prints one, is positive: a row
/// `--shift_3prime 1` moved; a run without the field says nothing about the shift.
fn shift_field_consistent(e: &Entry) -> Res<bool> {
    if SHIFT_LENGTH_FIELDS
        .iter()
        .all(|field| !e.fields.contains(field))
    {
        return Ok(true);
    }
    positive_shift(e)
}

/// The kind of an HGVSc that names the entry's transcript (its stable id, with or without a
/// version) and describes an insertion or deletion: `del`, `dup` or `ins`; `None` for any
/// other notation.
fn indel_notation_of<'a>(e: &Entry, hgvsc: &'a str) -> Option<&'a str> {
    let (name, _, body) = pytext::partition(hgvsc, ":");
    if stable_id(name) != stable_id(&e.feature) {
        return None;
    }
    let found = INDEL_NOTATION.captures(body)?;
    found.get(3).map(|kind| &kind.as_str()[..3])
}

/// The intron offset a coordinate ends with: a sign that follows a digit, then digits to the
/// end of the text or to a final line feed. The leading minus of a 5' UTR position (`-3`)
/// is none.
fn intron_offset(coordinate: &str) -> Option<&str> {
    let body = coordinate.strip_suffix('\n').unwrap_or(coordinate);
    let digits = body.trim_end_matches(is_digit);
    if digits.len() == body.len() {
        return None;
    }
    let mut before = digits.chars().rev();
    if !matches!(before.next(), Some('+' | '-')) || !before.next().is_some_and(is_digit) {
        return None;
    }
    Some(&body[digits.len() - 1..])
}

/// A signed run of decimal digits as a plain integer: no `+`, no leading zero, the digits
/// of any script as `0` to `9`.
fn plain_integer(text: &str) -> String {
    let (negative, digits) = match text.strip_prefix('-') {
        Some(digits) => (true, digits),
        None => (false, text.strip_prefix('+').unwrap_or(text)),
    };
    let ascii: String = digits
        .chars()
        .filter_map(pytext::decimal_value)
        .filter_map(|value| char::from_digit(value, 10))
        .collect();
    let significant = ascii.trim_start_matches('0');
    match (significant.is_empty(), negative) {
        (true, _) => "0".to_string(),
        (false, true) => format!("-{significant}"),
        (false, false) => significant.to_string(),
    }
}

/// The intron offset of an HGVS cDNA coordinate as the HGVSIntronOffset plugin writes it:
/// `9` for `29+9`, `-944` for `30-944`, `0` for an exonic coordinate.
pub fn coordinate_intron_offset(coordinate: &str) -> String {
    intron_offset(coordinate).map_or_else(|| "0".to_string(), plain_integer)
}

/// Whether an HGVS cDNA coordinate names an intron base: it ends with an intron offset, or
/// is `*-M`, the distance before the stop codon's last base that Ensembl VEP prints without
/// that base's coordinate (TranscriptVariationAllele.pm:2854-2859).
fn intron_coordinate(coordinate: &str) -> bool {
    intron_offset(coordinate).is_some() || coordinate.starts_with("*-")
}

/// A key ordering HGVS cDNA coordinates along the transcript: the 5' UTR (`-N`, region -1),
/// the coding sequence (`N`, region 0), the 3' UTR (`*N`, region 1), each by its position,
/// then an intron offset (`+M` after the exon base, `-M` before it). `*-M` orders after
/// every coding coordinate and before `*1`. A text of another shape is the region alone.
fn cdna_coordinate_order(coordinate: &str) -> Res<(i8, i128, i128)> {
    let (region, number) = if let Some(number) = coordinate.strip_prefix('*') {
        (1, number)
    } else if let Some(number) = coordinate.strip_prefix('-') {
        (-1, number)
    } else {
        (0, coordinate)
    };
    let position = number.trim_start_matches(is_digit);
    let (position, offset) = number.split_at(number.len() - position.len());
    let signed = offset
        .strip_prefix(['+', '-'])
        .is_some_and(pytext::isdecimal);
    if !offset.is_empty() && !signed {
        return Ok((region, 0, 0));
    }
    let or_zero = |digits| if str::is_empty(digits) { "0" } else { digits };
    let position = whole(or_zero(position), HGVS_COORDINATE)?;
    let offset = whole(or_zero(offset), "HGVS coordinate offset")?;
    Ok((
        region,
        if region < 0 { -position } else { position },
        offset,
    ))
}

/// An HGVS cDNA coordinate as a position on one integer scale along the transcript: the 5'
/// UTR `-N` at `-N`, the coding sequence `N` at `N - 1` (`c.-1` and `c.1` are neighbours),
/// the 3' UTR `*N` a billion further, since its distance from the coding sequence is the
/// coding length the notation does not carry, each with its intron offset added.
fn cdna_coordinate_index(coordinate: &str) -> Res<i128> {
    let (region, position, offset) = cdna_coordinate_order(coordinate)?;
    Ok(match region {
        -1 => position + offset,
        0 => position - 1 + offset,
        _ => 1_000_000_000 + position + offset,
    })
}

/// The count of bases from the cDNA coordinate `first` to `last`, both exonic and in one
/// region of the transcript (5' UTR, coding sequence or 3' UTR); `None` where either carries
/// an intron offset or the two lie in different regions: the scale of
/// `cdna_coordinate_index` carries no distance across the stop codon or an intron.
fn exonic_run(first: &str, last: &str) -> Res<Option<i128>> {
    let (first_region, _, first_offset) = cdna_coordinate_order(first)?;
    let (last_region, _, last_offset) = cdna_coordinate_order(last)?;
    if first_region != last_region || first_offset != 0 || last_offset != 0 {
        return Ok(None);
    }
    let end = cdna_coordinate_index(last)?;
    Ok(Some(end - cdna_coordinate_index(first)? + 1))
}

/// The `cdna_coordinate_index` of a notation's one or two coordinates.
fn notation_coordinate_indexes(notation: &Notation) -> Res<Vec<i128>> {
    [Some(notation.start), notation.end]
        .into_iter()
        .flatten()
        .map(cdna_coordinate_index)
        .collect()
}

/// An inserted allele after `steps` steps of a forward walk of `perform_shift`: its first
/// base moved to its end at each step (TranscriptVariationAllele.pm:338-341).
pub fn turned_allele(allele: &str, steps: i128) -> String {
    if allele.is_empty() {
        return String::new();
    }
    let turn = steps.rem_euclid(allele.len() as i128) as usize;
    format!("{}{}", &allele[turn..], &allele[..turn])
}

/// The steps `perform_shift` walks an insertion of `inserted` bases at most: its loop runs
/// over the 1,000 bases of flank less the allele's length, inclusive, and over the whole
/// flank for a longer allele (TranscriptVariationAllele.pm:303-306). Going 5', on a
/// reverse-strand transcript, the loop counts from 1 to the flank less the allele's length
/// plus 1, which is 0 for an allele of 1,001 bases and is left at that (:304 resets a
/// negative limit alone): that allele is not walked at all.
pub fn shift_bound(inserted: i128, towards_5prime: bool) -> i128 {
    if inserted <= SHIFT_FLANK_BASES {
        SHIFT_FLANK_BASES + 1 - inserted
    } else if towards_5prime && inserted == SHIFT_FLANK_BASES + 1 {
        0
    } else {
        SHIFT_FLANK_BASES
    }
}

/// Whether a coordinate is an exon position alone: an optional `*`, an optional `-` and
/// digits.
fn exon_position(coordinate: &str) -> bool {
    let digits = coordinate.strip_prefix('*').unwrap_or(coordinate);
    pytext::isdecimal(digits.strip_prefix('-').unwrap_or(digits))
}

/// The exon position that one of two coordinates is and the other carries with the intron
/// part `part` (`+1` or `-1`); `None` where neither does.
fn position_beside<'a>(first: &'a str, last: &'a str, part: &str) -> Option<&'a str> {
    [(first, last), (last, first)]
        .into_iter()
        .find_map(|(carrying, other)| {
            let position = carrying.strip_suffix(part)?;
            (exon_position(position) && position == other).then_some(position)
        })
}

/// Whether the entry's allele is two or more bases over a Location of as many, with no
/// record alleles to read: the text an insertion and a substitution of that length both
/// print, which `indel_row` reads as a substitution.
fn equal_span_allele_row(e: &Entry) -> Res<bool> {
    let Some(location) = e.location.as_deref() else {
        return Ok(false);
    };
    if e.indel.is_some() || e.allele.len() < 2 || !is_sequence(&e.allele) {
        return Ok(false);
    }
    let (start, end) = location_ends(location);
    if !(pytext::isdigit(start) && pytext::isdigit(end)) {
        return Ok(false);
    }
    Ok(span_bases(start, end)? == length(&e.allele))
}

/// Whether the engine row's two HGVSIntronOffset fields are the intron offsets of the two
/// coordinates of its `del` notation `hgvsc`, one coordinate standing for both, in either
/// order: the plugin reads them before the coordinates are ordered for printing. A
/// coordinate before the stop codon's last base (`*-5`) prints an intron distance without
/// its exon base, so a notation carrying one says nothing here.
fn deletion_offsets_follow_notation(a: &Entry, hgvsc: &str) -> bool {
    let Some((_, first, last)) = deletion_notation(pytext::partition(hgvsc, ":").2) else {
        return false;
    };
    if first.starts_with("*-") || last.starts_with("*-") {
        return true;
    }
    let mut written = intron_offset_pair(a).map(str::to_string);
    written.sort_unstable();
    let mut read = [
        coordinate_intron_offset(first),
        coordinate_intron_offset(last),
    ];
    read.sort_unstable();
    written == read
}

/// `hgvs_shifted_indel_written_span_gate` on one differing field.
///
/// `_var2transcript_slice_coords` (TranscriptVariationAllele.pm:2708-2713) maps the span as
/// written onto the transcript slice and refuses one wholly outside it (:2729), and
/// `hgvs_transcript` refuses a clamped end the shift's offset carries past the last base
/// (:1472), so an insertion or deletion the shift carries into a transcript from wholly
/// before it, or over its last base, has no HGVSc and no HGVS_OFFSET on a row annotated
/// inside the transcript; the engine describes the shifted span.
///
/// The shape: the reference's HGVSc, HGVS_OFFSET and HGVSIntronOffset fields are all empty,
/// the engine's HGVSc is an insertion or deletion notation of the row's transcript beside a
/// non-zero HGVS_OFFSET, the differing field is one of those (an intron offset an integer,
/// and beside a `del` notation the one its coordinates give), and the row is an indel's
/// annotated inside the transcript whose SHIFT_LENGTH, where printed, is positive. An allele
/// of two or more bases at a Location of as many reads the same for an insertion and a
/// substitution; there the engine's `dup` or `ins` notation names the row an insertion's.
/// On a row carrying `BAM_EDIT` `OK` on both sides the differing field may also be the
/// REFSEQ_OFFSET, the engine's a non-zero integer: the return at :1472 comes before the
/// offset is stored (:1517, :1621), so the reference's row prints none.
pub fn shifted_indel_written_span_gate(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
) -> Res<bool> {
    if !blank(ev) || blank(av) {
        return Ok(false);
    }
    let Some(actual) = actual else {
        return Ok(false);
    };
    if !shift_field_consistent(e)? {
        return Ok(false);
    }
    if SHIFTED_INDEL_GATE_FIELDS
        .iter()
        .any(|field| !blank(field_of(e, &[*field])))
    {
        return Ok(false);
    }
    // A row annotated outside the transcript has no notation on either engine.
    if consequence_set(&e.consequence)
        .split(',')
        .any(|term| OUTSIDE_TRANSCRIPT_TERMS.contains(&term))
    {
        return Ok(false);
    }
    let hgvsc = if field_named(name, &HGVSC_FIELDS) {
        av
    } else {
        field_of(actual, &HGVSC_FIELDS)
    };
    let offset = if field_named(name, &HGVS_OFFSET_FIELDS) {
        av
    } else {
        field_of(actual, &HGVS_OFFSET_FIELDS)
    };
    let intron_field = field_named(name, &HGVS_INTRON_OFFSET_FIELDS);
    if intron_field && !is_integer(av) {
        return Ok(false);
    }
    if field_named(name, &REFSEQ_OFFSET_FIELDS)
        && !(is_non_zero_integer(av)
            && [e, actual]
                .iter()
                .all(|row| field_of(row, &BAM_EDIT_FIELDS) == "OK"))
    {
        return Ok(false);
    }
    let Some(kind) = indel_notation_of(e, hgvsc) else {
        return Ok(false);
    };
    if !is_non_zero_integer(offset) {
        return Ok(false);
    }
    if intron_field && kind == "del" && !deletion_offsets_follow_notation(actual, hgvsc) {
        return Ok(false);
    }
    Ok(indel_row(e)? || (matches!(kind, "dup" | "ins") && equal_span_allele_row(e)?))
}

/// `hgvs_intron_offset_abutting_insertion` on one differing HGVSIntronOffset field.
///
/// HGVSIntronOffset.pm:62 asks for the notation on every Transcript row, past the
/// `within_feature` test that keeps HGVSc off a row outside its transcript
/// (OutputFactory.pm:1730); an insertion directly beside the transcript's first or last base
/// has one slice coordinate inside (TranscriptVariationAllele.pm:2729) and is clamped onto
/// that base (:2732-2733), so the offsets stored at :1543-1545 belong to a notation the row
/// never prints.
///
/// The shape: the reference fills both of the plugin's fields and the engine neither, on an
/// insertion's row both engines annotate outside the transcript, its one term
/// `upstream_gene_variant` or `downstream_gene_variant` beside a DISTANCE of 0, with no
/// HGVSc and no HGVS_OFFSET on either side. The pair is 0 and 0; or 1 and 0 on a downstream
/// row of a reverse-strand transcript, where :2802 wraps the position past the last base to
/// the far exon and :2822-2826 build "last base +1"; or 0 and -1 on an upstream row of a
/// forward-strand transcript, where `_clip_alleles` moves the notation's end to the position
/// before the first base (:2231-2238). A row without STRAND admits either.
pub fn intron_offset_abutting_insertion(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
) -> bool {
    let Some(actual) = actual else {
        return false;
    };
    if !field_named(name, &HGVS_INTRON_OFFSET_FIELDS) {
        return false;
    }
    if blank(ev) || !blank(av) || field_of(e, &[name]) != ev {
        return false;
    }
    if intron_offset_pair(actual).iter().any(|value| !blank(value)) {
        return false;
    }
    let term = consequence_set(&e.consequence);
    if !OUTSIDE_TRANSCRIPT_TERMS.contains(&term.as_str())
        || consequence_set(&actual.consequence) != term
    {
        return false;
    }
    for entry in [e, actual] {
        if field_of(entry, &DISTANCE_FIELDS) != "0" {
            return false;
        }
        if HGVSC_FIELDS
            .iter()
            .chain(&HGVS_OFFSET_FIELDS)
            .any(|field| !blank(field_of(entry, &[*field])))
        {
            return false;
        }
    }
    if !insertion_row(e) {
        return false;
    }
    let pair = intron_offset_pair(e);
    if pair == ["0", "0"] {
        return true;
    }
    let strand = field_of(e, &STRAND_FIELDS);
    if term == "downstream_gene_variant" {
        return matches!(strand, "-1" | ABSENT) && pair == ["1", "0"];
    }
    matches!(strand, "1" | ABSENT) && pair == ["0", "-1"]
}

/// Whether an engine notation's coordinates agree with the row's cDNA position under a run
/// that prints `SHIFT_LENGTH`, where that position is the shifted span's. In `n.` numbering,
/// which is cDNA numbering, a plain coordinate on both sides must be equal; in `c.`
/// numbering two plain coordinates of one region span as many bases as the row's two plain
/// sides. A side the row leaves open (`?`), a coordinate with an intron offset and a `c.`
/// span across a region boundary say nothing, as does a run without `SHIFT_LENGTH`.
fn span_matches_cdna_position(e: &Entry, numbering: &str, start: &str, end: &str) -> Res<bool> {
    if SHIFT_LENGTH_FIELDS
        .iter()
        .all(|field| !e.fields.contains(field))
    {
        return Ok(true);
    }
    let (row_start, row_end) = match e.fields.get(CDNA_START_FIELDS[0]) {
        // Under --total_length the position carries the transcript's length after a slash.
        Some(position) => {
            let (row_start, separator, row_end) =
                pytext::partition(pytext::partition(position, "/").0, "-");
            (
                row_start,
                if separator.is_empty() {
                    row_start
                } else {
                    row_end
                },
            )
        }
        None => (
            e.fields.get_or(CDNA_START_FIELDS[1], ""),
            e.fields.get_or(CDNA_END_FIELD, ""),
        ),
    };
    if numbering == "n" {
        return Ok([(row_start, start), (row_end, end)].iter().all(
            |(row_position, coordinate)| {
                !(pytext::isdigit(row_position)
                    && pytext::isdigit(coordinate)
                    && row_position != coordinate)
            },
        ));
    }
    if !(pytext::isdigit(row_start) && pytext::isdigit(row_end)) {
        return Ok(true);
    }
    let (start_region, start_position, start_offset) = cdna_coordinate_order(start)?;
    let (end_region, end_position, end_offset) = cdna_coordinate_order(end)?;
    if start_offset != 0 || end_offset != 0 || start_region != end_region {
        return Ok(true);
    }
    let row_last = whole(row_end, CDNA_POSITION)?;
    Ok(end_position - start_position == row_last - whole(row_start, CDNA_POSITION)?)
}

/// `hgvs_shifted_deletion_clamped_span` on a differing HGVSc.
///
/// `_var2transcript_slice_coords` clamps a written span straddling the transcript's first
/// base onto that base (TranscriptVariationAllele.pm:2732-2733) before the offset moves it,
/// so a deletion the shift carries inside is described from the base after the shift where
/// its shifted span starts earlier.
///
/// The shape: two `del` notations of the row's transcript and one numbering, ending at the
/// same cDNA coordinate, the reference's starting later than the engine's, on a deletion's
/// row (allele `-`) whose HGVS_OFFSET is a non-zero integer on the reference side and whose
/// SHIFT_LENGTH, where printed, is positive; the engine's coordinates agree with the row's
/// cDNA position (`span_matches_cdna_position`).
pub fn shifted_deletion_clamped_span(e: &Entry, ev: &str, av: &str) -> Res<bool> {
    if blank(ev) || blank(av) || e.allele != "-" || !shift_field_consistent(e)? {
        return Ok(false);
    }
    if !is_non_zero_integer(field_of(e, &HGVS_OFFSET_FIELDS)) {
        return Ok(false);
    }
    let (expected_name, _, expected_body) = pytext::partition(ev, ":");
    let (actual_name, _, actual_body) = pytext::partition(av, ":");
    let (
        Some((reference_numbering, reference_start, reference_end)),
        Some((numbering, engine_start, engine_end)),
    ) = (
        deletion_notation(expected_body),
        deletion_notation(actual_body),
    )
    else {
        return Ok(false);
    };
    if expected_name != actual_name || indel_notation_of(e, ev).is_none() {
        return Ok(false);
    }
    if reference_numbering != numbering || reference_end != engine_end {
        return Ok(false);
    }
    let engine_order = cdna_coordinate_order(engine_start)?;
    if engine_order >= cdna_coordinate_order(reference_start)? {
        return Ok(false);
    }
    span_matches_cdna_position(e, numbering, engine_start, engine_end)
}

/// `hgvs_shifted_deletion_clamped_span` on a differing `HGVS_IntronStartOffset`: the row's
/// HGVSc pair has the class's shape (`shifted_deletion_clamped_span`) and each engine's
/// offset is the intron offset of its own notation's start coordinate.
pub fn shifted_deletion_start_offsets(e: &Entry, a: &Entry, ev: &str, av: &str) -> Res<bool> {
    let notations = [field_of(e, &HGVSC_FIELDS), field_of(a, &HGVSC_FIELDS)];
    if !shifted_deletion_clamped_span(e, notations[0], notations[1])? {
        return Ok(false);
    }
    let start_offset = |notation: &str| {
        deletion_notation(pytext::partition(notation, ":").2)
            .map(|(_, start, _)| coordinate_intron_offset(start))
    };
    Ok(start_offset(notations[0]).as_deref() == Some(ev)
        && start_offset(notations[1]).as_deref() == Some(av))
}

/// `hgvs_shifted_deletion_clamped_span` on one differing field: the start offset beside the
/// engine's row (`shifted_deletion_start_offsets`), any other field of the class as the
/// HGVSc (`shifted_deletion_clamped_span`).
pub fn shifted_deletion_clamped_span_field(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
) -> Res<bool> {
    if field_named(name, &HGVS_INTRON_START_OFFSET_FIELDS) {
        return match actual {
            Some(a) => shifted_deletion_start_offsets(e, a, ev, av),
            None => Ok(false),
        };
    }
    shifted_deletion_clamped_span(e, ev, av)
}

/// A side of a row's cDNA position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Start,
    End,
}

/// The side of the row's cDNA position that lies outside the transcript: the start for a
/// `?` before the dash of the text formats' column (`?-6`), the end for one after it
/// (`8485-?`), in the JSON format the side of `cdna_start` and `cdna_end` left out while the
/// other is present (JSON.pm:340-341); `None` for a closed position.
fn cdna_open_side(e: &Entry) -> Option<Side> {
    if let Some(position) = e.fields.get(CDNA_START_FIELDS[0]) {
        let (start, _, end) = pytext::partition(pytext::partition(position, "/").0, "-");
        if start == "?" {
            return Some(Side::Start);
        }
        return (end == "?").then_some(Side::End);
    }
    let has_start = !blank(e.fields.get_or(CDNA_START_FIELDS[1], ABSENT));
    let has_end = !blank(e.fields.get_or(CDNA_END_FIELD, ABSENT));
    if has_start == has_end {
        return None;
    }
    Some(if has_start { Side::End } else { Side::Start })
}

/// How far apart, in bases, the two HGVSc notations of an allele running over a transcript's
/// edge can lie: the bases of the allele's span, which bound the bases inside the
/// transcript, plus one base of adjacency. The span is the row's Location; without one and
/// with the start open, the closed cDNA end counts the bases inside; `None` with the end
/// open and no Location.
fn edge_row_reach(e: &Entry, side: Side) -> Res<Option<i128>> {
    if let Some(location) = e.location.as_deref() {
        let (start, end) = location_ends(location);
        if pytext::isdigit(start) {
            let bases = if pytext::isdigit(end) {
                span_bases(start, end)?
            } else {
                1
            };
            return Ok(Some(bases + 1));
        }
    }
    if side == Side::Start {
        let cdna_end = cdna_closed_end(e);
        if pytext::isdigit(cdna_end) {
            return Ok(Some(whole(cdna_end, CDNA_POSITION)? + 1));
        }
    }
    Ok(None)
}

/// Whether the reference's notation is one the class takes against an engine that prints
/// none on the row of an allele that deletes bases over the transcript's edge: the row's
/// allele holds fewer bases than its Location spans, and the notation is a `del`, a `delins`
/// or a substitution, what `hgvs_variant_notation` makes of the whole alternate against the
/// bases the clamp left inside (TranscriptVariationAllele.pm:1474-1483, 2732-2733). A row
/// without a Location and an allele as long as its span or longer stay differences.
fn deletion_over_the_edge(e: &Entry, reference: &Notation) -> Res<bool> {
    let Some(location) = e.location.as_deref() else {
        return Ok(false);
    };
    let (start, end) = location_ends(location);
    if !(pytext::isdigit(start) && pytext::isdigit(end)) {
        return Ok(false);
    }
    let spanned = span_bases(start, end)?;
    let held = length(&e.allele);
    if !(0 < held && held < spanned) {
        return Ok(false);
    }
    Ok(
        reference.plain == Some("del")
            || reference.kind == Some("delins")
            || reference.substitution,
    )
}

/// `hgvs_clipped_reference_against_whole_alternate` on a differing HGVSc.
///
/// `_var2transcript_slice_coords` clamps a written span running over a transcript's edge
/// onto the bases inside (TranscriptVariationAllele.pm:2732-2733) and `hgvs_transcript`
/// passes that span with the whole alternate allele (:1474-1483), which
/// `hgvs_variant_notation` reads against the reference bases of the span alone
/// (Utils/Sequence.pm:526); the engine cuts the alternate as the clamp cut the reference and
/// describes the change inside the transcript, none when there is none.
///
/// The shape: on a row of a sequence allele whose cDNA position has an open side
/// (`cdna_open_side`), the reference's notation names the row's transcript, or is empty
/// where the whole alternate equals the bases inside (Utils/Sequence.pm:529-531), and the
/// engine's names the same transcript and numbering, every coordinate within the row's span
/// plus one of every coordinate of the reference's (`edge_row_reach`; a row whose reach
/// cannot be read passes on the transcript and numbering alone), or is empty where the
/// change lies outside. On an edited model's row whose position is open at its end, a
/// reference coordinate written with the row's REFSEQ_OFFSET as its intron distance stands
/// for every base from the exon base it counts from to itself
/// (`slide_past_the_last_base`), and a reference substitution that names the row's whole
/// reference allele is read as the substitution of one base (`respelled_substitution`).
/// Against an empty reference the engine's is a `del`; against
/// an empty engine notation the reference's is an `ins` or a `dup`, or, on a row whose allele
/// is shorter than its span, a `del`, a `delins` or a substitution
/// (`deletion_over_the_edge`).
/// A deletion's row (allele `-`) and a closed cDNA position stay differences.
pub fn clipped_reference_against_whole_alternate(e: &Entry, ev: &str, av: &str) -> Res<bool> {
    let side = cdna_open_side(e);
    if e.allele == "-" {
        return Ok(false);
    }
    let respelled = respelled_substitution(e, ev);
    let reference = clipped_notation(ev, &e.feature).or_else(|| {
        respelled
            .as_deref()
            .and_then(|notation| clipped_notation(notation, &e.feature))
    });
    let engine = clipped_notation(av, &e.feature);
    let Some(side) = side else {
        return Ok(false);
    };
    if blank(ev) {
        return Ok(engine.is_some_and(|engine| engine.plain == Some("del")));
    }
    let Some(reference) = reference else {
        return Ok(false);
    };
    if blank(av) {
        return Ok(reference.kind == Some("ins")
            || reference.plain == Some("dup")
            || deletion_over_the_edge(e, &reference)?);
    }
    let Some(engine) = engine else {
        return Ok(false);
    };
    if pytext::partition(ev, ":").0 != pytext::partition(av, ":").0
        || engine.numbering != reference.numbering
    {
        return Ok(false);
    }
    let Some(reach) = edge_row_reach(e, side)? else {
        return Ok(true);
    };
    let engine_indexes = notation_coordinate_indexes(&engine)?;
    let slide = slide_past_the_last_base(e, side);
    let mut reference_stretches: Vec<(i128, i128)> = Vec::new();
    for coordinate in [Some(reference.start), reference.end].into_iter().flatten() {
        let index = cdna_coordinate_index(coordinate)?;
        let distance = cdna_coordinate_order(coordinate)?.2;
        let slid = if slide == Some(distance) { distance } else { 0 };
        reference_stretches.push((index - slid, index));
    }
    Ok(engine_indexes.iter().all(|c| {
        reference_stretches
            .iter()
            .all(|(first, last)| first - reach <= *c && *c <= last + reach)
    }))
}

/// The notation of a substitution whose reference allele holds two or more bases, written as
/// the substitution of the first of them: on an edited model (`BAM_EDIT` `OK`)
/// `hgvs_transcript` types the notation against the bases the clamp left inside
/// (TranscriptVariationAllele.pm:1474-1483) and then writes the row's whole reference allele
/// over it (:1510-1513), a base outside the transcript among them. `None` for any other row
/// or notation.
fn respelled_substitution(e: &Entry, hgvsc: &str) -> Option<String> {
    if field_of(e, &BAM_EDIT_FIELDS) != "OK" {
        return None;
    }
    let (head, alternate) = hgvsc.split_once('>')?;
    let position = head.trim_end_matches(['A', 'C', 'G', 'T', 'N']);
    let named = head.len() - position.len();
    if named < 2 || !position.ends_with(|last: char| last.is_ascii_digit()) {
        return None;
    }
    let first = &head[..position.len() + 1];
    Some(format!("{first}>{alternate}"))
}

/// The bases by which the reference slides a coordinate past the last base of an edited
/// model, as its row states them: its REFSEQ_OFFSET, a positive whole number, on a row
/// carrying `BAM_EDIT` `OK` whose cDNA position is open at its end. `hgvs_transcript` adds
/// the offset to the clamped span's slice positions (TranscriptVariationAllele.pm:1517,
/// :1524), and a position past the last exon is written as that exon's end with the distance
/// as an intron offset (:2798-2829), the end numbered without the model's edits. The count of
/// bases the edits add before the last base, which the engine's coordinate carries, is not
/// among the row's fields: it lies between none and the offset. `None` for any other row.
fn slide_past_the_last_base(e: &Entry, side: Side) -> Option<i128> {
    if side != Side::End || field_of(e, &BAM_EDIT_FIELDS) != "OK" {
        return None;
    }
    let offset = field_of(e, &REFSEQ_OFFSET_FIELDS);
    if offset.is_empty() || !offset.bytes().all(|digit| digit.is_ascii_digit()) {
        return None;
    }
    offset.parse().ok().filter(|slide| *slide > 0)
}

/// Whether the engine's `ins` notation writes the bases a longer allele adds inside the span
/// of a row whose cDNA position is open at its start (`?-N`): the row's allele is a sequence
/// longer than its Location, the notation inserts as many bases as it is longer, between two
/// exon coordinates next to each other, and under `n.` numbering, whose coordinates are cDNA
/// positions, at or beside a base inside: its first coordinate at most the position's closed
/// end. Under `c.` numbering the coordinates' place is not read; a row without a Location is
/// refused.
fn insertion_inside_the_span(e: &Entry, engine: &Notation) -> Res<bool> {
    let first = engine.start;
    let Some(last) = engine.end else {
        return Ok(false);
    };
    if engine.kind != Some("ins") || intron_coordinate(first) || intron_coordinate(last) {
        return Ok(false);
    }
    let after = cdna_coordinate_index(last)?;
    if after - cdna_coordinate_index(first)? != 1 {
        return Ok(false);
    }
    let Some(location) = e.location.as_deref() else {
        return Ok(false);
    };
    if !is_sequence(&e.allele) {
        return Ok(false);
    }
    let (start, end) = location_ends(location);
    if !pytext::isdigit(start) || (!end.is_empty() && !pytext::isdigit(end)) {
        return Ok(false);
    }
    let spanned = if end.is_empty() {
        1
    } else {
        span_bases(start, end)?
    };
    if length(engine.bases) != length(&e.allele) - spanned {
        return Ok(false);
    }
    if engine.numbering != "n" {
        return Ok(true);
    }
    let closed = cdna_closed_end(e);
    if !pytext::isdigit(closed) || !pytext::isdigit(first) || !pytext::isdigit(last) {
        return Ok(false);
    }
    let first = whole(first, HGVS_COORDINATE)?;
    Ok(first <= whole(closed, CDNA_POSITION)?)
}

/// Whether the row is one on which the engine prints an `ins` where the reference prints no
/// HGVSc: a sequence allele's row on a reverse-strand transcript (STRAND `-1`) whose cDNA
/// position is open at its start, the transcript's first base lying at its genomic end, the
/// engine's notation an `ins` of the bases a longer allele adds inside the span
/// (`insertion_inside_the_span`). There the clip leaves Ensembl VEP an insertion at slice
/// position 0, for which `_get_cDNA_position` returns no coordinate
/// (TranscriptVariationAllele.pm:2834).
fn insertion_without_reference_notation(e: &Entry, engine: &Notation) -> Res<bool> {
    if e.allele == "-"
        || cdna_open_side(e) != Some(Side::Start)
        || field_of(e, &STRAND_FIELDS) != "-1"
    {
        return Ok(false);
    }
    insertion_inside_the_span(e, engine)
}

/// The two integers `COORDINATE_PARTS_AS_READ` finds in a coordinate, an absent second one
/// read as 0.
fn coordinate_parts_as_read(coordinate: &str) -> Res<(i128, i128)> {
    let found = COORDINATE_PARTS_AS_READ
        .captures(coordinate)
        .ok_or_else(|| not_an_integer(coordinate, HGVS_COORDINATE))?;
    let part = |index| {
        found
            .get(index)
            .map_or(Ok(0), |digits| whole(digits.as_str(), HGVS_COORDINATE))
    };
    Ok((part(1)?, part(2)?))
}

/// The two values the HGVSIntronOffset plugin writes beside a notation, the start's and the
/// end's: the intron parts of the two coordinates as `hgvs_transcript` reads them
/// (TranscriptVariationAllele.pm:1530-1533, an absent part read as 0; stored at :1543 and
/// :1545) before it orders the pair for printing (:1535-1538). The routine holds an `ins`
/// from its 3' flank to its 5' flank, the clip leaving the start one above the end
/// (:2245-2248), and turns the two unless the end holds a star, so the start is the
/// coordinate printed second, or the one printed first where the 5' flank lies past the stop
/// codon; it holds a `dup` and every other kind in its printed order, one coordinate
/// standing for both. `None` for a pair the routine would have printed the other way round,
/// an `ins` of one coordinate, and a coordinate before the stop codon's last base (`*-5`).
fn notation_intron_offsets(notation: &Notation) -> Res<Option<[String; 2]>> {
    let (first, last) = (notation.start, notation.last());
    let insertion = notation.kind == Some("ins");
    if (insertion && notation.end.is_none()) || first.starts_with("*-") || last.starts_with("*-") {
        return Ok(None);
    }
    let first_parts = coordinate_parts_as_read(first)?;
    let last_parts = coordinate_parts_as_read(last)?;
    if insertion && last_parts > first_parts && !first.contains('*') {
        return Ok(Some([last_parts.1.to_string(), first_parts.1.to_string()]));
    }
    if first_parts > last_parts && !last.contains('*') {
        return Ok(None);
    }
    Ok(Some([first_parts.1.to_string(), last_parts.1.to_string()]))
}

/// The intron offsets of the two coordinates of a notation in their printed order, one
/// coordinate standing for both: the pair the HGVSIntronOffset plugin writes beside an
/// insertion's `ins` or `dup`, whose span `hgvs_variant_notation` holds from the 5' flank to
/// the 3' flank (Utils/Sequence.pm:577-579). `None` for a coordinate before the stop codon's
/// last base (`*-5`).
fn printed_intron_offsets(notation: &Notation) -> Option<[String; 2]> {
    let (first, last) = (notation.start, notation.last());
    if first.starts_with("*-") || last.starts_with("*-") {
        return None;
    }
    Some([
        coordinate_intron_offset(first),
        coordinate_intron_offset(last),
    ])
}

/// Whether a row's two values of the plugin's fields are the pair `read` from a notation.
fn pair_is(written: &[&str; 2], read: Option<[String; 2]>) -> bool {
    read.is_some_and(|read| written[0] == read[0] && written[1] == read[1])
}

/// Whether the two rows are the reference's `ins` from one base before a forward-strand
/// transcript's first base against the engine's duplication of the bases the allele adds,
/// read from the two rows alone in four parts, all of which hold.
///
/// The row: its allele holds more bases than its Location spans, `A` bases over `S`
/// positions; its STRAND is `1`; its cDNA position is open at its start and closed at its
/// end, `?-N`, with `N` above 0 and below `S`. The reference: an `ins` of the row's
/// transcript whose two coordinates are one coordinate `X` carrying the intron part `-1` and
/// the same `X` carrying none; it writes the allele's bases less the last `N`; it carries no
/// HGVS_OFFSET; the plugin's pair is 0 and -1. An `ins` the clip leaves has consumed every
/// reference base `hgvs_transcript` handed it (TranscriptVariationAllele.pm:2231-2242,
/// :2271), so it was handed `N` bases, fewer than the span holds, which the clamp alone does
/// (:2732-2733). The engine: a `dup` of the same transcript and numbering whose coordinates
/// carry no intron part and cover `A` - `S` bases; its HGVS_OFFSET a whole number `k`, an
/// empty one read as 0; the plugin's pair 0 and 0. The tie: from `X` to the duplication's
/// last coordinate are `N` + `k` exonic bases of one region.
///
/// A row that carries BAM_EDIT is refused: on an edited model `hgvs_transcript` respells the
/// notation's alleles after the clip (:1510-1514). A row without a Location or without
/// STRAND and a reverse-strand transcript's row are refused.
fn insertion_before_the_first_base_against_duplication(
    e: &Entry,
    actual: &Entry,
    reference: Option<&Notation>,
    engine: &Notation,
    pairs: &[[&str; 2]; 2],
) -> Res<bool> {
    let (Some(reference), Some(location)) = (reference, e.location.as_deref()) else {
        return Ok(false);
    };
    if !is_bases(&e.allele) {
        return Ok(false);
    }
    if field_of(e, &STRAND_FIELDS) != "1"
        || [e, actual]
            .iter()
            .any(|entry| !blank(field_of(entry, &BAM_EDIT_FIELDS)))
    {
        return Ok(false);
    }
    let (start, end) = location_ends(location);
    let end = if end.is_empty() { start } else { end };
    if !(pytext::isdigit(start) && pytext::isdigit(end)) {
        return Ok(false);
    }
    let spanned = span_bases(start, end)?;
    let held = length(&e.allele);
    let added = held - spanned;
    let closed = cdna_closed_end(e);
    if added < 1 || cdna_open_side(e) != Some(Side::Start) || !pytext::isdigit(closed) {
        return Ok(false);
    }
    let inside = whole(closed, CDNA_POSITION)?;
    let first = reference.start;
    let Some(last) = reference.end else {
        return Ok(false);
    };
    if !(0 < inside && inside < spanned) || reference.kind != Some("ins") {
        return Ok(false);
    }
    let Some(first_base) = position_beside(first, last, "-1") else {
        return Ok(false);
    };
    // The allele holds more bases than its span and the span more than `inside`.
    if reference.bases != &e.allele[..(held - inside) as usize] {
        return Ok(false);
    }
    if !blank(field_of(e, &HGVS_OFFSET_FIELDS)) || pairs[0] != ["0", "-1"] {
        return Ok(false);
    }
    if engine.plain != Some("dup") || engine.numbering != reference.numbering {
        return Ok(false);
    }
    let duplicated_last = engine.last();
    if exonic_run(engine.start, duplicated_last)? != Some(added) || pairs[1] != ["0", "0"] {
        return Ok(false);
    }
    let carried = field_of(actual, &HGVS_OFFSET_FIELDS);
    if !blank(carried) && !pytext::is_ascii_digits(carried) {
        return Ok(false);
    }
    let run = exonic_run(first_base, duplicated_last)?;
    let moved = if blank(carried) {
        0
    } else {
        whole(carried, HGVS_OFFSET)?
    };
    Ok(run == Some(inside + moved))
}

/// `hgvs_clipped_reference_against_whole_alternate` on one differing HGVSIntronOffset field,
/// in four shapes. The plugin writes the two fields beside whatever notation
/// `hgvs_transcript` returns, each the intron part of a coordinate as the plugin reads it
/// (`notation_intron_offsets`), and a row without a notation gets neither.
///
/// The reference alone fills both fields, on a row whose HGVSc the class takes
/// (`clipped_reference_against_whole_alternate`) against an engine that prints none. The
/// engine alone fills both, on a row whose reference prints no HGVSc where the engine prints
/// an `ins` inside the span at a reverse-strand transcript's first base
/// (`insertion_without_reference_notation`). Or both fill both, on a row whose two notations
/// the class takes as a pair, each an `ins`, and the field named differs. In each of these
/// three a side's two values are the intron parts of its own notation. Or both fill both and
/// the end's field alone is asked, on a forward-strand transcript's row whose reference
/// prints an `ins` from one base before the transcript's first base, beside the pair 0 and
/// -1, and whose engine prints the duplication of the bases the allele adds, beside the pair
/// 0 and 0 (`insertion_before_the_first_base_against_duplication`).
///
/// Outside the fourth shape a row either side's HGVS_OFFSET is written on beside the
/// engine's fields stays a difference.
pub fn clipped_notation_intron_offsets(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
) -> Res<bool> {
    let Some(actual) = actual else {
        return Ok(false);
    };
    let pairs = [intron_offset_pair(e), intron_offset_pair(actual)];
    let filled = pairs.map(|pair| pair.iter().all(|value| !blank(value)));
    let empty = pairs.map(|pair| pair.iter().all(|value| blank(value)));
    let (hgvsc, engine_hgvsc) = (field_of(e, &HGVSC_FIELDS), field_of(actual, &HGVSC_FIELDS));
    let reference = clipped_notation(hgvsc, &e.feature);
    let engine = clipped_notation(engine_hgvsc, &e.feature);
    if filled[0] && empty[1] {
        if blank(ev)
            || !blank(av)
            || field_of(e, &[name]) != ev
            || blank(hgvsc)
            || !blank(engine_hgvsc)
        {
            return Ok(false);
        }
        if !clipped_reference_against_whole_alternate(e, hgvsc, engine_hgvsc)? {
            return Ok(false);
        }
        let Some(reference) = reference else {
            return Ok(false);
        };
        return Ok(pair_is(&pairs[0], notation_intron_offsets(&reference)?));
    }
    if !filled[1] || blank(av) || field_of(actual, &[name]) != av {
        return Ok(false);
    }
    let Some(engine) = engine else {
        return Ok(false);
    };
    if engine.plain == Some("dup") {
        if !field_named(name, &HGVS_INTRON_END_OFFSET_FIELDS) || field_of(e, &[name]) != ev {
            return Ok(false);
        }
        return insertion_before_the_first_base_against_duplication(
            e,
            actual,
            reference.as_ref(),
            &engine,
            &pairs,
        );
    }
    if engine.kind != Some("ins") {
        return Ok(false);
    }
    if [e, actual]
        .iter()
        .any(|entry| !blank(field_of(entry, &HGVS_OFFSET_FIELDS)))
    {
        return Ok(false);
    }
    if !pair_is(&pairs[1], notation_intron_offsets(&engine)?) {
        return Ok(false);
    }
    if empty[0] {
        return Ok(blank(ev) && blank(hgvsc) && insertion_without_reference_notation(e, &engine)?);
    }
    if !clipped_reference_against_whole_alternate(e, hgvsc, engine_hgvsc)? {
        return Ok(false);
    }
    if !filled[0] || blank(ev) || field_of(e, &[name]) != ev || ev == av {
        return Ok(false);
    }
    let Some(reference) = reference else {
        return Ok(false);
    };
    if reference.kind != Some("ins") {
        return Ok(false);
    }
    Ok(pair_is(&pairs[0], notation_intron_offsets(&reference)?))
}

/// `hgvs_clipped_reference_against_whole_alternate` on one differing field: one of the
/// plugin's two by `clipped_notation_intron_offsets`, the HGVSc by
/// `clipped_reference_against_whole_alternate`.
pub fn clipped_reference_field(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
) -> Res<bool> {
    if field_named(name, &HGVS_INTRON_OFFSET_FIELDS) {
        return clipped_notation_intron_offsets(e, actual, name, ev, av);
    }
    clipped_reference_against_whole_alternate(e, ev, av)
}

/// `hgvs_added_bases_left_at_the_span_end` on one differing field of a row whose allele
/// holds more bases than its Location spans, the bases an allele adds to a reference it
/// keeps whole, where the reference prints an HGVSc of the row's transcript and no
/// HGVS_OFFSET: `_return_3prime` shifts a record whose class is insertion or deletion alone
/// (TranscriptVariationAllele.pm:120-123), so the added bases stay where `_clip_alleles`
/// leaves them, at the end of the record's span. Three shapes, read on the two rows
/// together.
///
/// Past a reverse-strand transcript's last base: the reference's notation is an `ins` of the
/// added bases between a coordinate and the same coordinate with an intron offset of 1
/// (:2802, :2820-2827), beside which it writes the plugin's start offset 1 and end offset 0.
/// The engine's is a duplication of as many bases ending at that coordinate, with no
/// HGVS_OFFSET and its own two intron offsets, where the last added base repeats the kept
/// base beside it: the HGVSc and the start offset are the class. Or the engine writes no
/// HGVSc, no HGVS_OFFSET and neither field, where the last added base does not repeat that
/// kept base: the HGVSc and both fields are the class.
///
/// Carried along a run: the reference's notation is an `ins` between two neighbouring exonic
/// bases, or a `dup`, and the engine's a `dup` of as many bases or an `ins` of the same
/// bases turned, a number of bases further 3' that its HGVS_OFFSET states, negative on a
/// reverse strand, the two rows writing the same intron offsets: the HGVSc and the
/// HGVS_OFFSET are the class. Where nothing moves, the engine's is the `dup` of the bases
/// that end at the reference's `ins`, beside no HGVS_OFFSET: the HGVSc is the class.
pub fn added_bases_left_at_the_span_end(
    e: &Entry,
    actual: Option<&Entry>,
    name: &str,
) -> Res<bool> {
    let (Some(actual), Some(location)) = (actual, e.location.as_deref()) else {
        return Ok(false);
    };
    if !is_bases(&e.allele) {
        return Ok(false);
    }
    let (start, end) = location_ends(location);
    let end = if end.is_empty() { start } else { end };
    if !(pytext::isdigit(start) && pytext::isdigit(end)) {
        return Ok(false);
    }
    let added = length(&e.allele) - span_bases(start, end)?;
    let Some(reference) = clipped_notation(field_of(e, &HGVSC_FIELDS), &e.feature) else {
        return Ok(false);
    };
    if added <= 0 || !blank(field_of(e, &HGVS_OFFSET_FIELDS)) {
        return Ok(false);
    }
    let engine_hgvsc = field_of(actual, &HGVSC_FIELDS);
    let engine_offset = field_of(actual, &HGVS_OFFSET_FIELDS);
    let engine = clipped_notation(engine_hgvsc, &e.feature);
    if !blank(engine_hgvsc) && engine.is_none_or(|engine| engine.numbering != reference.numbering) {
        return Ok(false);
    }
    let reference_offsets = intron_offset_pair(e);
    let engine_offsets = intron_offset_pair(actual);
    let strand = field_of(e, &STRAND_FIELDS);
    let (first, last) = (reference.start, reference.end);
    let past = last.and_then(|last| position_beside(first, last, "+1"));
    if let (Some("ins"), Some(past)) = (reference.kind, past) {
        let bases = reference.bases;
        let allele = e.allele.as_bytes();
        if length(bases) != added
            || reverse_complement(&e.allele[..allele.len().min(bases.len())]) != bases
        {
            return Ok(false);
        }
        // A run without the plugin writes the two fields on neither row.
        let without_plugin = all_blank([&reference_offsets, &engine_offsets]);
        if !(reference_offsets == ["1", "0"] || without_plugin) || !matches!(strand, "-1" | ABSENT)
        {
            return Ok(false);
        }
        let Some(kept) = allele.get(bases.len()) else {
            return Err(Refusal(format!(
                "Location {} spans no base",
                pytext::repr(location)
            )));
        };
        let repeats_the_last_base = allele[0] == *kept;
        let Some(engine) = engine else {
            if repeats_the_last_base
                || !blank(engine_offset)
                || engine_offsets.iter().any(|value| !blank(value))
            {
                return Ok(false);
            }
            return Ok(
                field_named(name, &HGVSC_FIELDS) || field_named(name, &HGVS_INTRON_OFFSET_FIELDS)
            );
        };
        if engine.plain != Some("dup")
            || engine.last() != past
            || exonic_run(engine.start, engine.last())? != Some(added)
        {
            return Ok(false);
        }
        if !repeats_the_last_base
            || !blank(engine_offset)
            || !(engine_offsets == ["0", "0"] || without_plugin)
        {
            return Ok(false);
        }
        return Ok(
            field_named(name, &HGVSC_FIELDS) || field_named(name, &HGVS_INTRON_START_OFFSET_FIELDS)
        );
    }
    let Some(engine) = engine else {
        return Ok(false);
    };
    if !(blank(engine_offset) || is_non_zero_integer(engine_offset)) {
        return Ok(false);
    }
    if reference_offsets != engine_offsets
        || blank(reference_offsets[0]) != blank(reference_offsets[1])
    {
        return Ok(false);
    }
    // The point the reference describes the added bases at, as the exonic base before it,
    // and the point the engine does. An `ins` past the stop codon prints its 3' flank first
    // (TranscriptVariationAllele.pm:1535-1538).
    let described = if reference.kind == Some("ins") {
        let Some(last) = last else {
            return Ok(false);
        };
        if length(reference.bases) != added {
            return Ok(false);
        }
        if exonic_run(first, last)? == Some(2) {
            first
        } else if exonic_run(last, first)? == Some(2) {
            last
        } else {
            return Ok(false);
        }
    } else if reference.plain == Some("dup") {
        if exonic_run(first, reference.last())? != Some(added) {
            return Ok(false);
        }
        reference.last()
    } else {
        return Ok(false);
    };
    let placed = if engine.plain == Some("dup") {
        if exonic_run(engine.start, engine.last())? != Some(added) {
            return Ok(false);
        }
        engine.last()
    } else if engine.kind == Some("ins") {
        let Some(placed_last) = engine.end else {
            return Ok(false);
        };
        if exonic_run(engine.start, placed_last)? != Some(2) || length(engine.bases) != added {
            return Ok(false);
        }
        engine.start
    } else {
        return Ok(false);
    };
    let carried = exonic_run(described, placed)?;
    if blank(engine_offset) {
        let typed =
            reference.kind == Some("ins") && engine.plain == Some("dup") && carried == Some(1);
        return Ok(typed && field_named(name, &HGVSC_FIELDS));
    }
    let Some(carried) = carried.filter(|carried| *carried >= 2) else {
        return Ok(false);
    };
    let moved = whole(engine_offset.trim_start_matches('-'), HGVS_OFFSET)?;
    if moved != carried - 1
        || (engine_offset.starts_with('-') != (strand == "-1") && strand != ABSENT)
    {
        return Ok(false);
    }
    if engine.kind == Some("ins")
        && (reference.kind != Some("ins") || engine.bases != turned_allele(reference.bases, moved))
    {
        return Ok(false);
    }
    Ok(field_named(name, &HGVSC_FIELDS) || field_named(name, &HGVS_OFFSET_FIELDS))
}

/// Whether the record's own alleles insert an entry's allele whole between two neighbouring
/// positions, a pure insertion: the bases the allele stands in for are none
/// (`Entry.reference`), or the record changes length and its reference allele is empty
/// (`Entry.indel`, `Entry.insertion`). `None` where the entry carries neither, as a row
/// under an identifier does that holds no reference bases from the lines of the input.
fn inserted_whole(e: &Entry) -> Option<bool> {
    match (e.reference.borrow().as_deref(), e.indel, e.insertion) {
        (Some(replaced), _, _) => Some(replaced.is_empty()),
        (None, Some(indel), Some(insertion)) => Some(indel && insertion),
        _ => None,
    }
}

/// Whether `allele` is the bases `replaced` with bases added before them, after them or
/// between two of them: an allele that keeps whole the bases it stands in for.
fn keeps_whole(allele: &str, replaced: &str) -> bool {
    let (allele, replaced) = (allele.as_bytes(), replaced.as_bytes());
    allele.len() > replaced.len()
        && (0..=replaced.len()).any(|kept| {
            allele.starts_with(&replaced[..kept]) && allele.ends_with(&replaced[kept..])
        })
}

/// `hgvs_added_bases_written_outside_the_transcript` on one differing field of a row both
/// engines annotate outside the transcript, its one term `upstream_gene_variant` or
/// `downstream_gene_variant`, whose allele is an insertion's or holds more bases than its
/// Location spans.
///
/// OutputFactory.pm:1730 fills the HGVSc column only on a row within the feature, so the
/// reference prints no HGVSc and no HGVS_OFFSET there and writes the HGVSIntronOffset
/// plugin's two fields or neither (HGVSIntronOffset.pm:62). The engine prints a `dup` or an
/// `ins` of the row's transcript, of as many bases as the allele holds or adds to its span
/// where its coordinates count them, and both plugin fields, the intron offsets of that
/// notation's coordinates in their printed order, or neither where the reference writes
/// neither, as in a run without the plugin. On a downstream row nothing moved: the engine's
/// is a `dup` beside no HGVS_OFFSET. On an upstream row the 3' rule carried the bases
/// inside: its HGVS_OFFSET is a non-zero integer, negative on a reverse strand, larger in
/// size than the row's DISTANCE on an insertion's row and no smaller than it on the row of
/// an allele that holds more bases than its span, whose DISTANCE counts from that span. The
/// differing HGVSc, HGVS_OFFSET and plugin fields are the class. An entry without a Location
/// is read under the one its record gives its allele (`at_its_location`) and stays a
/// difference where its record gives it none.
///
/// The bases a row adds are read from its record's own alleles where the entry carries them
/// (`inserted_whole`): a row whose first column spells the alleles, a row under an
/// identifier that holds the reference bases the lines of the input give it
/// (`Entry.reference`), and an entry a caller builds from a record. A pure insertion's
/// notation counts the whole allele. The notation of an allele longer than the bases it
/// stands in for, which it keeps whole where the entry holds them (`keeps_whole`), counts
/// the added bases, and its HGVS_OFFSET may be the DISTANCE. A deletion spelt as a sequence,
/// a substitution and a replacement that keeps neither end of its bases add none and stay
/// differences.
///
/// A row under an identifier that holds no reference bases carries no alleles and is read
/// from its Location (`insertion_row`), the limit of that form; `score` reads the lines of
/// its input for the rows of the classes `data/field_classes.tsv` lists under
/// `#record_reference`. Two neighbouring positions read as an insertion's, whose notation
/// counts the whole allele, and an allele longer than its Location spans reads as one that
/// adds its surplus. Two bases replaced by any sequence, a deletion of one of two bases
/// spelt as a sequence among them, share an insertion's Location: a `dup` of the allele's
/// length is taken on such a row, and one of the surplus on the row of a pure insertion
/// longer than two bases or of a replacement that keeps neither end. One row the alleles
/// take is left: an allele that adds to the two bases it keeps reads as an insertion's
/// where the notation's coordinates count no bases, so an HGVS_OFFSET equal to the DISTANCE
/// stays a difference.
pub fn added_bases_written_outside_the_transcript(
    e: &Row,
    actual: Option<&Entry>,
    name: &str,
    ev: &str,
    av: &str,
) -> Res<bool> {
    let e = &at_its_location(e);
    let (Some(actual), Some(location)) = (actual, e.location.as_deref()) else {
        return Ok(false);
    };
    if !is_bases(&e.allele) {
        return Ok(false);
    }
    let (span_start, span_end) = location_ends(location);
    let span_end = if span_end.is_empty() {
        span_start
    } else {
        span_end
    };
    if !(pytext::isdigit(span_start) && pytext::isdigit(span_end)) {
        return Ok(false);
    }
    let positions = span_bases(span_start, span_end)?;
    let held = length(&e.allele);
    let surplus = held - positions;
    // A row that carries none of its record's alleles is read from its Location, which a
    // replacement of two bases shares with an insertion.
    let (inserted, kept) = match inserted_whole(e) {
        Some(inserted) => (
            inserted,
            !inserted
                && e.reference
                    .borrow()
                    .as_deref()
                    .is_none_or(|replaced| keeps_whole(&e.allele, replaced)),
        ),
        None => (insertion_row(e), true),
    };
    let added = kept && surplus > 0;
    if !inserted && !added {
        return Ok(false);
    }
    let counts = |counted: i128| (inserted && counted == held) || (added && counted == surplus);
    let term = consequence_set(&e.consequence);
    if !OUTSIDE_TRANSCRIPT_TERMS.contains(&term.as_str())
        || consequence_set(&actual.consequence) != term
    {
        return Ok(false);
    }
    if HGVSC_FIELDS
        .iter()
        .chain(&HGVS_OFFSET_FIELDS)
        .any(|field| !blank(field_of(e, &[*field])))
    {
        return Ok(false);
    }
    let reference_offsets = intron_offset_pair(e);
    if blank(reference_offsets[0]) != blank(reference_offsets[1]) {
        return Ok(false);
    }
    let Some(engine) = clipped_notation(field_of(actual, &HGVSC_FIELDS), &e.feature) else {
        return Ok(false);
    };
    if !(engine.plain == Some("dup") || engine.kind == Some("ins")) {
        return Ok(false);
    }
    let (first, last) = (engine.start, engine.last());
    let counted = if engine.kind == Some("ins") {
        if engine.end.is_none() {
            return Ok(false);
        }
        Some(length(engine.bases))
    } else {
        exonic_run(first, last)?
    };
    if counted.is_some_and(|counted| !counts(counted)) {
        return Ok(false);
    }
    if first.starts_with("*-") || last.starts_with("*-") {
        return Ok(false);
    }
    let written = intron_offset_pair(actual);
    // A run without the plugin writes the two fields on neither row.
    let without_plugin = all_blank([&reference_offsets, &written]);
    if !without_plugin
        && (written[0] != coordinate_intron_offset(first)
            || written[1] != coordinate_intron_offset(last))
    {
        return Ok(false);
    }
    let offset = field_of(actual, &HGVS_OFFSET_FIELDS);
    let strand = field_of(e, &STRAND_FIELDS);
    if term == "downstream_gene_variant" {
        if !blank(offset) || engine.plain != Some("dup") {
            return Ok(false);
        }
    } else {
        let distance = field_of(e, &DISTANCE_FIELDS);
        if !is_non_zero_integer(offset) || !pytext::isdigit(distance) {
            return Ok(false);
        }
        // An insertion's DISTANCE counts from the nearer base beside its point, a longer
        // allele's from its span.
        let added_to_a_span = added && (!inserted || counted == Some(surplus));
        let moved = whole(offset.trim_start_matches('-'), HGVS_OFFSET)?;
        if moved < whole(distance, "DISTANCE")? + i128::from(!added_to_a_span) {
            return Ok(false);
        }
        if strand != ABSENT && offset.starts_with('-') != (strand == "-1") {
            return Ok(false);
        }
    }
    if field_named(name, &HGVS_INTRON_OFFSET_FIELDS) {
        return Ok(!blank(av) && field_of(actual, &[name]) == av && ev != av);
    }
    Ok(blank(ev)
        && !blank(av)
        && (field_named(name, &HGVSC_FIELDS) || field_named(name, &HGVS_OFFSET_FIELDS)))
}

/// `hgvs_insertion_left_where_the_shift_stops` on one differing field of an insertion's row
/// whose reference HGVS_OFFSET is, in size, exactly the bound of `perform_shift` for the
/// allele's length (`shift_bound`) beside an `ins` or a `dup` of the row's transcript: the
/// walk ended on its bound, not on a base that differs. Where that bound is 0, an allele of
/// 1,001 bases on a row whose STRAND is -1, the reference walked nothing and its row carries
/// no HGVS_OFFSET; that row is read only where its Location names the two bases around an
/// insertion's point (an entry without a Location under the one its record gives its
/// allele), the record's alleles, where the entry carries them, hold no reference base, and
/// the reference's notation is the `ins` of the allele's own bases or a `dup` of as many.
/// The reference maps each end of the copy to a coordinate alone
/// (TranscriptVariationAllele.pm:1524-1525), so two exonic coordinates count the copied
/// bases less those of every intron between them. Handed the exon spans of the run's JSON
/// cache (`exon_spans`), that `dup` is held to them (`copy_beside_its_point`): its two
/// coordinates are those of the first and last base of a run of as many bases as the allele
/// holds beside the point, so a count is the class's only where the introns between the two
/// bases explain it, and a row the cache does not place stays a difference, as does a `dup`
/// with one coordinate in the 3' UTR and the other before it, whose distance the notation
/// does not carry. Without the cache the class has no intron to count by: two exonic
/// coordinates of one region are taken from 2 bases up to as many as the allele holds,
/// whatever intron would explain the count, but for two fewer, the bases an allele that
/// keeps the two bases of its span adds; a `dup` of one coordinate stays a difference; and a
/// `dup` with a coordinate in an intron, or with its coordinates in two regions, is not
/// counted. Two shapes of the engine's row.
///
/// Further 3': a `dup` or an `ins` in the same numbering beside an HGVS_OFFSET of the same
/// sign, negative where the reference carries none, and a larger size. A `dup` whose
/// coordinates are exonic names as many bases as the allele holds or, handed the exon spans,
/// as many as the exons hold of a run beside its own point (`copy_beside_its_point`);
/// without the cache an engine `dup` around an intron stays a difference. An `ins` writes
/// the allele's bases as the transcript's strand reads them, turned by the steps its
/// HGVS_OFFSET states. The HGVSc, the HGVS_OFFSET and the plugin's two fields are the class.
///
/// Past the last base: no HGVSc and neither plugin field, beside no HGVS_OFFSET or the
/// reference's own. The HGVSc and the two fields are the class, and the HGVS_OFFSET where
/// the engine writes none.
///
/// Where the rows carry the plugin's fields, each row's two are the intron parts of its own
/// notation's coordinates in their printed order (`printed_intron_offsets`).
///
/// `exon_spans` is `None` where the run names no cache.
pub fn insertion_left_where_the_shift_stops(
    e: &Row,
    actual: Option<&Entry>,
    name: &str,
    exon_spans: Option<&super::exon_spans::ExonSpans>,
) -> Res<bool> {
    let Some(actual) = actual else {
        return Ok(false);
    };
    if !is_bases(&e.allele) {
        return Ok(false);
    }
    let stopped = field_of(e, &HGVS_OFFSET_FIELDS);
    let Some(reference) = clipped_notation(field_of(e, &HGVSC_FIELDS), &e.feature) else {
        return Ok(false);
    };
    if !(reference.kind == Some("ins") || reference.plain == Some("dup")) {
        return Ok(false);
    }
    let held = length(&e.allele);
    let (towards_5prime, bound);
    if blank(stopped) {
        // The one walk that never starts: its row has no HGVS_OFFSET to read a direction
        // from. An allele that keeps bases of its span carries none either and is never
        // walked (:123), so the row must be an insertion's own.
        let placed = at_its_location(e);
        let (span_start, span_end) = location_ends(placed.location.as_deref().unwrap_or(""));
        if !(pytext::isdigit(span_start) && pytext::isdigit(span_end)) {
            return Ok(false);
        }
        let after = whole(span_end, LOCATION_POSITION)?;
        if after != whole(span_start, LOCATION_POSITION)? + 1 || e.insertion == Some(false) {
            return Ok(false);
        }
        towards_5prime = field_of(e, &STRAND_FIELDS) == "-1";
        bound = shift_bound(held, towards_5prime);
        if bound != 0 {
            return Ok(false);
        }
        if reference.kind == Some("ins") {
            if reference.bases != reverse_complement(&e.allele) {
                return Ok(false);
            }
        } else if let Some(exon_spans) = exon_spans {
            if !copy_beside_its_point(exon_spans, &placed, &reference, held, 0, towards_5prime)? {
                return Ok(false);
            }
        } else {
            // Without the exon spans an allele that keeps the two bases of its span, which
            // adds two fewer than it holds, is told from an insertion by that count alone.
            let copied = exonic_run(reference.start, reference.last())?;
            if copied.is_some_and(|copied| !(2..=held).contains(&copied) || copied == held - 2) {
                return Ok(false);
            }
        }
    } else {
        if !is_non_zero_integer(stopped) {
            return Ok(false);
        }
        towards_5prime = stopped.starts_with('-');
        bound = shift_bound(held, towards_5prime);
        if whole(stopped.trim_start_matches('-'), HGVS_OFFSET)? != bound {
            return Ok(false);
        }
    }
    let reference_offsets = intron_offset_pair(e);
    let engine_offsets = intron_offset_pair(actual);
    // A run without the plugin writes the two fields on neither row.
    let without_plugin = all_blank([&reference_offsets, &engine_offsets]);
    if !without_plugin && !pair_is(&reference_offsets, printed_intron_offsets(&reference)) {
        return Ok(false);
    }
    let engine_hgvsc = field_of(actual, &HGVSC_FIELDS);
    let moved = field_of(actual, &HGVS_OFFSET_FIELDS);
    if blank(engine_hgvsc) {
        if engine_offsets.iter().any(|value| !blank(value)) || !(blank(moved) || moved == stopped) {
            return Ok(false);
        }
        if field_named(name, &HGVS_OFFSET_FIELDS) {
            return Ok(blank(moved));
        }
        return Ok(
            field_named(name, &HGVSC_FIELDS) || field_named(name, &HGVS_INTRON_OFFSET_FIELDS)
        );
    }
    let Some(engine) = clipped_notation(engine_hgvsc, &e.feature) else {
        return Ok(false);
    };
    if engine.numbering != reference.numbering {
        return Ok(false);
    }
    if !is_non_zero_integer(moved) || moved.starts_with('-') != towards_5prime {
        return Ok(false);
    }
    let further = whole(moved.trim_start_matches('-'), HGVS_OFFSET)? - bound;
    if further <= 0
        || !(without_plugin || pair_is(&engine_offsets, printed_intron_offsets(&engine)))
    {
        return Ok(false);
    }
    if engine.plain == Some("dup") {
        let duplicated = exonic_run(engine.start, engine.last())?;
        if duplicated.is_some_and(|duplicated| duplicated != held) {
            let Some(exon_spans) = exon_spans else {
                return Ok(false);
            };
            let steps = bound + further;
            let along = if towards_5prime { -steps } else { steps };
            let placed = at_its_location(e);
            if !copy_beside_its_point(exon_spans, &placed, &engine, held, along, towards_5prime)? {
                return Ok(false);
            }
        }
    } else if engine.kind == Some("ins") {
        let placed = if moved.starts_with('-') {
            reverse_complement(&e.allele)
        } else {
            e.allele.clone()
        };
        if engine.bases != turned_allele(&placed, bound + further) {
            return Ok(false);
        }
    } else {
        return Ok(false);
    }
    Ok(field_named(name, &HGVSC_FIELDS)
        || field_named(name, &HGVS_OFFSET_FIELDS)
        || field_named(name, &HGVS_INTRON_OFFSET_FIELDS))
}

/// The readings of a printed coordinate as its region, exon coordinate and intron part
/// (`cdna_coordinate_order`). A 3' UTR coordinate without an intron part has two:
/// `_get_cDNA_position` writes an intron base numbered from the stop codon's last base
/// without the `+` of its part (TranscriptVariationAllele.pm:2854-2859), so `*5` is the
/// fifth base of the 3' UTR or the fifth base of the intron after that last base. A
/// coordinate read as 0 without an intron part has none: the notation numbers no base 0,
/// and `cdna_coordinate_order` reads a text of another shape so.
fn coordinate_readings(coordinate: &str) -> Res<Vec<(i8, i128, i128)>> {
    let (region, at, part) = cdna_coordinate_order(coordinate)?;
    if at == 0 && part == 0 {
        return Ok(Vec::new());
    }
    let mut readings = vec![(region, at, part)];
    if region == 1 && part == 0 && at > 0 {
        readings.push((region, 0, at));
    }
    Ok(readings)
}

/// Whether the two coordinates of the `dup` `notation` are those of the first and last base
/// of a run of `held` consecutive bases beside its point, read from the exon spans of the
/// run's JSON cache. The notation routine makes a `dup` of the `held` bases after the point
/// or of the `held` bases before it (Utils/Sequence.pm:577-587), and `hgvs_transcript` maps
/// the first and the last of them each to a coordinate alone
/// (TranscriptVariationAllele.pm:1524-1525). `placed` is the row at its Location
/// (`at_its_location`), which names the two bases around the written point; the notation's
/// point lies `moved` bases along the chromosome from it, `moved` the HGVS_OFFSET beside the
/// notation and 0 where none is written. A run fits where both of its ends lie within the
/// transcript, each coordinate's intron part is that of its base
/// (`TranscriptExons::exon_base`: 0 for an exon base) and the two coordinates lie as far
/// apart as the exon bases they are numbered from lie along the exons: for two exonic
/// coordinates, a count of `held` less the bases of every intron between the two. The
/// cache's coding bounds are not read, so the coordinates are held to their distance and
/// their intron parts, not to the place of either, and one coordinate in the 3' UTR beside
/// one before it, whose distance the notation does not carry, fit no run. Nothing fits
/// either where the Location does not name the two bases around a point, the cache holds no
/// transcript of the row's Feature on the Location's chromosome (`ExonSpans::transcript`),
/// or that transcript's strand is not `reverse`.
fn copy_beside_its_point(
    exon_spans: &super::exon_spans::ExonSpans,
    placed: &Entry,
    notation: &Notation,
    held: i128,
    moved: i128,
    reverse: bool,
) -> Res<bool> {
    let first = coordinate_readings(notation.start)?;
    let last = coordinate_readings(notation.last())?;
    let location = placed.location.as_deref().unwrap_or("");
    let chromosome = location
        .split_once(':')
        .map_or("", |(chromosome, _)| chromosome);
    let (start, end) = location_ends(location);
    let position = |digits: &str| {
        digits
            .bytes()
            .all(|digit| digit.is_ascii_digit())
            .then(|| digits.parse::<i128>().ok())
            .flatten()
            .filter(|at| at.unsigned_abs() <= LARGEST_NUMBER)
    };
    let (Some(before), Some(after)) = (position(start), position(end)) else {
        return Ok(false);
    };
    if after != before + 1 {
        return Ok(false);
    }
    let Some(transcript) = exon_spans.transcript(chromosome, &placed.feature)? else {
        return Ok(false);
    };
    if transcript.reverse != reverse {
        return Ok(false);
    }
    let point = before + moved;
    for (low, high) in [(point - held + 1, point), (point + 1, point + held)] {
        let (five_prime, three_prime) = if reverse { (high, low) } else { (low, high) };
        let (from, from_part) = transcript.exon_base(five_prime);
        let (to, to_part) = transcript.exon_base(three_prime);
        if from == 0 || to == 0 {
            continue;
        }
        // `c.-1` and `c.1` are neighbours, and a 3' UTR coordinate counts from the stop
        // codon's last base, whose place the exon spans do not give.
        let scaled = |region: i8, at: i128| at - i128::from(region == 0);
        let fits = first.iter().any(|&(first_region, first_at, first_part)| {
            last.iter().any(|&(last_region, last_at, last_part)| {
                (first_part, last_part) == (from_part, to_part)
                    && (first_region == 1) == (last_region == 1)
                    && scaled(last_region, last_at) - scaled(first_region, first_at) == to - from
            })
        });
        if fits {
            return Ok(true);
        }
    }
    Ok(false)
}
