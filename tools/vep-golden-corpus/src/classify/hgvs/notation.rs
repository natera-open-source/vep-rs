// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HGVSc notations the HGVS field classes read: the patterns of a notation's body, the
//! place of a cDNA coordinate along its transcript, and the readings of a row the classes
//! share.
//!
//! Every pattern is read as CPython's `re` reads it, its `\d` the ten ASCII digits, and a
//! number is read within 64 bits: a notation with a longer number is no notation here.

use crate::classify::fields::FieldRow;
use crate::pytext::py_int;
use regex::Regex;
use std::sync::LazyLock;

/// An HGVSc body of any kind but a repeat count: the numbering, the one or two cDNA
/// coordinates, then a substitution (`-47A>C`; the reference base is `N` over a masked
/// reference window), a `delins` or `ins` with its bases, or a plain inversion, deletion or
/// duplication.
pub const HGVSC_CLIPPED_BODY: &str = concat!(
    r"(?P<numbering>[cn])\.(?P<start>\*?-?\d+(?:[+-]\d+)?)(?:_(?P<end>\*?-?\d+(?:[+-]\d+)?))?",
    r"(?:[ACGTN]>(?P<alternate_base>[ACGT])|(?P<kind>delins|ins)(?P<bases>[ACGT]+)|(?P<plain>inv|del|dup))",
);

/// An HGVSc body describing an insertion or deletion: the numbering, the one or two cDNA
/// coordinates and the kind (`c.520_541del`, `n.14dup`, an `ins` with its bases).
pub const HGVSC_INDEL_BODY: &str =
    r"([cn])\.(\*?-?\d+(?:[+-]\d+)?)(?:_(\*?-?\d+(?:[+-]\d+)?))?(del|dup|ins[ACGT]+)";

/// A cDNA coordinate of an HGVSc body (`-71`, `250+601`, `*12-3`); with each replaced by
/// `#` two bodies that differ in their coordinates alone compare equal.
pub const HGVSC_COORDINATE: &str = r"\*?-?\d+(?:[+-]\d+)?";

/// The signed integers of an HGVSc body: `250+601` gives 250 and 601, `*12-3` gives 12 and
/// -3.
pub const HGVSC_INTEGER: &str = r"-?\d+";

/// A non-zero `HGVS_OFFSET`, the only value Ensembl VEP prints for the field.
pub const NON_ZERO_INTEGER: &str = r"-?[1-9]\d*";

/// The exon coordinate and the intron offset `hgvs_transcript` reads off a cDNA coordinate
/// (TranscriptVariationAllele.pm:1530-1531): the first integer and the one after it.
pub const COORDINATE_PARTS_AS_READ: &str = r"(-?[0-9]+)\+?(-?[0-9]+)?";

/// A coordinate one base before an exonic one, as `_get_cDNA_position` writes the position
/// before a forward-strand transcript's first base.
pub const ONE_BEFORE_COORDINATE: &str = r"(\*?-?\d+)-1";

/// The fields the HGVSIntronOffset plugin writes beside an HGVSc, the start's and the end's.
pub const HGVS_INTRON_OFFSET_FIELDS: [&str; 2] = ["HGVS_IntronStartOffset", "HGVS_IntronEndOffset"];

/// `re.compile(pattern)` for `fullmatch`.
fn whole(pattern: &str) -> Regex {
    Regex::new(&format!(r"(?-u)\A(?:{pattern})\z")).expect("a fixed pattern compiles")
}

/// `re.compile(pattern)` for `search`, `findall` and `sub`.
fn anywhere(pattern: &str) -> Regex {
    Regex::new(&format!("(?-u){pattern}")).expect("a fixed pattern compiles")
}

static CLIPPED_BODY: LazyLock<Regex> = LazyLock::new(|| whole(HGVSC_CLIPPED_BODY));
static INDEL_BODY: LazyLock<Regex> = LazyLock::new(|| whole(HGVSC_INDEL_BODY));
static COORDINATE: LazyLock<Regex> = LazyLock::new(|| anywhere(HGVSC_COORDINATE));
static INTEGER: LazyLock<Regex> = LazyLock::new(|| anywhere(HGVSC_INTEGER));
static NON_ZERO: LazyLock<Regex> = LazyLock::new(|| whole(NON_ZERO_INTEGER));
static PARTS_AS_READ: LazyLock<Regex> = LazyLock::new(|| anywhere(COORDINATE_PARTS_AS_READ));
static ONE_BEFORE: LazyLock<Regex> = LazyLock::new(|| whole(ONE_BEFORE_COORDINATE));
static INTRON_PART: LazyLock<Regex> = LazyLock::new(|| anywhere(r"\d[+-]\d+\n?\z"));

/// `(row or {}).get(name, "")`: a field of a row, empty for a row or a field that is absent.
pub fn field<'a>(row: Option<&'a FieldRow>, name: &str) -> &'a str {
    row.and_then(|row| row.get(name)).map_or("", String::as_str)
}

/// Whether a row holds exactly `value` under `name` (`row.get(name) == value`).
pub fn holds(row: Option<&FieldRow>, name: &str, value: &str) -> bool {
    row.and_then(|row| row.get(name)).map(String::as_str) == Some(value)
}

/// The two fields of the HGVSIntronOffset plugin on a row, the start's and the end's.
pub fn intron_offset_pair(row: Option<&FieldRow>) -> (&str, &str) {
    (
        field(row, HGVS_INTRON_OFFSET_FIELDS[0]),
        field(row, HGVS_INTRON_OFFSET_FIELDS[1]),
    )
}

/// Whether a field is empty as either engine prints an empty field: no text, or `-`.
pub fn blank(value: &str) -> bool {
    matches!(value, "" | "-")
}

/// `str.isdigit()` for a text of ASCII digits: one digit or more and nothing else.
pub fn digits(text: &str) -> bool {
    !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit())
}

/// `int(text)` within 64 bits, widened for the sums the classes make.
pub fn int(text: &str) -> Option<i128> {
    py_int(text).map(i128::from)
}

/// `len(text)`: the count of code points.
pub fn length(text: &str) -> i128 {
    i128::try_from(text.chars().count()).unwrap_or(i128::MAX)
}

/// Whether `text` is one or more of the letters of `alphabet`.
pub fn spells(text: &str, alphabet: &str) -> bool {
    !text.is_empty() && text.chars().all(|letter| alphabet.contains(letter))
}

/// Whether `text` is a non-zero integer as `NON_ZERO_INTEGER` reads one.
pub fn non_zero_integer(text: &str) -> bool {
    NON_ZERO.is_match(text)
}

/// The start and the end of a Location's span as texts: what follows the first `:`, cut at
/// its first `-`; the end is empty for a Location of one position.
pub fn span_texts(location: &str) -> (&str, &str) {
    let span = location.split_once(':').map_or(location, |(_, span)| span);
    span.split_once('-').unwrap_or((span, ""))
}

/// The count of positions a Location names: one for a single position, the end less the
/// start plus one for a span. `None` for a start or an end that is not ASCII digits.
pub fn span_positions(location: &str) -> Option<i128> {
    let (start, end) = span_texts(location);
    let end = if end.is_empty() { start } else { end };
    if !digits(start) || !digits(end) {
        return None;
    }
    Some(int(end)? - int(start)? + 1)
}

/// The two ends of a `cDNA_position` column: the text before its `/`, cut at its first `-`;
/// an end that lies outside the transcript reads `?`.
pub fn position_ends(position: &str) -> (&str, &str) {
    let span = position.split_once('/').map_or(position, |(span, _)| span);
    span.split_once('-').unwrap_or((span, ""))
}

/// The name an HGVSc opens with, before its first `:`.
pub fn notation_name(hgvsc: &str) -> &str {
    hgvsc.split_once(':').map_or(hgvsc, |(name, _)| name)
}

/// The body of an HGVSc that names the row's transcript (its stable id, with or without a
/// version), `None` for another transcript's notation.
fn body_on<'a>(hgvsc: &'a str, feature: &str) -> Option<&'a str> {
    let (name, body) = hgvsc.split_once(':').unwrap_or((hgvsc, ""));
    (unversioned(name) == unversioned(feature)).then_some(body)
}

/// A stable id without its version: the text before its first `.`.
fn unversioned(id: &str) -> &str {
    id.split_once('.').map_or(id, |(stem, _)| stem)
}

/// The place of a cDNA coordinate along its transcript: the region (-1 the 5' UTR, 0 the
/// coding sequence, 1 the 3' UTR), the position in it and the intron offset.
pub type CoordinateOrder = (i8, i128, i128);

/// The start of the sign and the value of an intron offset that ends a coordinate
/// (`re.search(r"[+-]\d+$", number)`).
fn trailing_offset(number: &str) -> Option<(usize, i128)> {
    let body = number.strip_suffix('\n').unwrap_or(number);
    let run = body.len() - body.bytes().rev().take_while(u8::is_ascii_digit).count();
    if run == body.len() || run == 0 || !matches!(body.as_bytes()[run - 1], b'+' | b'-') {
        return None;
    }
    Some((run - 1, int(&body[run - 1..])?))
}

/// A key ordering HGVS cDNA coordinates along the transcript: the 5' UTR (`-N`), the coding
/// sequence (`N`), the 3' UTR (`*N`), each by its position, then an intron offset (`+M`
/// after the exon base, `-M` before it). `*-M`, an intron distance before the stop codon's
/// last base that `_get_cDNA_position` prints without that base's coordinate
/// (TranscriptVariationAllele.pm:2854-2859), orders after every coding coordinate and
/// before `*1`. `None` for a text that is no coordinate.
pub fn cdna_coordinate_order(coordinate: &str) -> Option<CoordinateOrder> {
    let (region, number) = if let Some(number) = coordinate.strip_prefix('*') {
        (1, number)
    } else if let Some(number) = coordinate.strip_prefix('-') {
        (-1, number)
    } else {
        (0, coordinate)
    };
    let (number, offset) = match trailing_offset(number) {
        Some((sign, offset)) => (&number[..sign], offset),
        None => (number, 0),
    };
    let position = if number.is_empty() { 0 } else { int(number)? };
    Some((
        region,
        if region < 0 { -position } else { position },
        offset,
    ))
}

/// A coordinate's place on one integer scale along the transcript: the 5' UTR `-N` at `-N`,
/// the coding sequence `N` at `N - 1` (`c.-1` and `c.1` are neighbours), the 3' UTR `*N` a
/// billion further, since its distance from the coding sequence is the coding length the
/// notation does not carry, each with its intron offset added.
fn order_index((region, position, offset): CoordinateOrder) -> i128 {
    match region {
        -1 => position + offset,
        0 => position - 1 + offset,
        _ => 1_000_000_000 + position + offset,
    }
}

/// `order_index` of a coordinate, `None` for a text that is no coordinate.
pub fn cdna_coordinate_index(coordinate: &str) -> Option<i128> {
    cdna_coordinate_order(coordinate).map(order_index)
}

/// The count of bases from the cDNA coordinate `first` to `last`, both exonic and in one
/// region of the transcript (5' UTR, coding sequence or 3' UTR), `None` where either
/// carries an intron offset or the two lie in different regions: the scale
/// `cdna_coordinate_index` reads carries no distance across the stop codon or an intron.
pub fn exonic_run(first: &str, last: &str) -> Option<i128> {
    let (first, last) = (cdna_coordinate_order(first)?, cdna_coordinate_order(last)?);
    if first.0 != last.0 || first.2 != 0 || last.2 != 0 {
        return None;
    }
    Some(order_index(last) - order_index(first) + 1)
}

/// The intron offset of a coordinate as the HGVSIntronOffset plugin prints one.
pub fn intron_offset_text(coordinate: &str) -> Option<String> {
    cdna_coordinate_order(coordinate).map(|order| order.2.to_string())
}

/// Whether a coordinate carries an intron offset after an exon base (`\d[+-]\d+$`).
pub fn carries_intron_part(coordinate: &str) -> bool {
    INTRON_PART.is_match(coordinate)
}

/// The coordinate one base after `coordinate`, where `coordinate` is that one with the
/// intron part `-1` (`ONE_BEFORE_COORDINATE`).
pub fn one_before(coordinate: &str) -> Option<&str> {
    let found = ONE_BEFORE.captures(coordinate)?;
    found.get(1).map(|base| base.as_str())
}

/// The body of an HGVSc with every cDNA coordinate replaced by `#`.
pub fn coordinates_masked(body: &str) -> String {
    COORDINATE.replace_all(body, "#").into_owned()
}

/// The signed integers of an HGVSc body in order, `None` where one is past 64 bits.
pub fn body_integers(body: &str) -> Option<Vec<i128>> {
    INTEGER
        .find_iter(body)
        .map(|found| int(found.as_str()))
        .collect()
}

/// An HGVSc body `HGVSC_INDEL_BODY` reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IndelBody<'a> {
    pub numbering: &'a str,
    pub start: &'a str,
    pub end: Option<&'a str>,
    /// `del`, `dup`, or `ins` with its bases.
    pub kind: &'a str,
}

impl<'a> IndelBody<'a> {
    /// The last coordinate, the first standing for both where the body names one.
    pub fn last(&self) -> &'a str {
        self.end.unwrap_or(self.start)
    }
}

/// The body `HGVSC_INDEL_BODY` reads whole, `None` for any other text.
pub fn indel_body(body: &str) -> Option<IndelBody<'_>> {
    let found = INDEL_BODY.captures(body)?;
    let group = |index: usize| found.get(index).map(|part| part.as_str());
    let notation = IndelBody {
        numbering: group(1)?,
        start: group(2)?,
        end: group(3),
        kind: group(4)?,
    };
    cdna_coordinate_order(notation.start)?;
    cdna_coordinate_order(notation.last())?;
    Some(notation)
}

/// The indel body of an HGVSc that names the row's transcript, or `None` for another
/// transcript's notation or a notation of another kind.
pub fn names_transcript<'a>(hgvsc: &'a str, feature: &str) -> Option<IndelBody<'a>> {
    indel_body(body_on(hgvsc, feature)?)
}

/// An HGVSc body `HGVSC_CLIPPED_BODY` reads, by the pattern's groups.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ClippedBody<'a> {
    pub numbering: &'a str,
    pub start: &'a str,
    pub end: Option<&'a str>,
    /// The base a substitution writes.
    pub alternate_base: Option<&'a str>,
    /// `delins` or `ins`.
    pub kind: Option<&'a str>,
    /// The bases a `delins` or an `ins` writes.
    pub bases: Option<&'a str>,
    /// `inv`, `del` or `dup`.
    pub plain: Option<&'a str>,
}

impl<'a> ClippedBody<'a> {
    /// The last coordinate, the first standing for both where the body names one.
    pub fn last(&self) -> &'a str {
        self.end.unwrap_or(self.start)
    }

    pub fn is_ins(&self) -> bool {
        self.kind == Some("ins")
    }

    pub fn is_dup(&self) -> bool {
        self.plain == Some("dup")
    }

    /// The `cdna_coordinate_index` of the body's one or two coordinates.
    pub fn coordinate_indexes(&self) -> impl Iterator<Item = i128> + 'a {
        [Some(self.start), self.end]
            .into_iter()
            .flatten()
            .filter_map(cdna_coordinate_index)
    }
}

/// The body `HGVSC_CLIPPED_BODY` reads whole, `None` for any other text.
pub fn clipped_body(body: &str) -> Option<ClippedBody<'_>> {
    let found = CLIPPED_BODY.captures(body)?;
    let group = |name: &str| found.name(name).map(|part| part.as_str());
    let notation = ClippedBody {
        numbering: group("numbering")?,
        start: group("start")?,
        end: group("end"),
        alternate_base: group("alternate_base"),
        kind: group("kind"),
        bases: group("bases"),
        plain: group("plain"),
    };
    cdna_coordinate_order(notation.start)?;
    cdna_coordinate_order(notation.last())?;
    Some(notation)
}

/// The clipped body of an HGVSc that names the row's transcript, or `None` for an empty
/// value, another transcript's notation or a notation of another shape.
pub fn clipped_notation<'a>(hgvsc: &'a str, feature: &str) -> Option<ClippedBody<'a>> {
    if blank(hgvsc) {
        return None;
    }
    clipped_body(body_on(hgvsc, feature)?)
}

/// The exon coordinate and the intron part of a coordinate as `hgvs_transcript` reads them,
/// an absent intron part read as 0.
fn parts_as_read(coordinate: &str) -> Option<(i128, i128)> {
    let found = PARTS_AS_READ.captures(coordinate)?;
    let exon = int(found.get(1)?.as_str())?;
    let intron = match found.get(2) {
        Some(part) => int(part.as_str())?,
        None => 0,
    };
    Some((exon, intron))
}

/// The two values the HGVSIntronOffset plugin writes beside a notation, the start's and
/// the end's: the intron parts of the two coordinates as `hgvs_transcript` reads them
/// (TranscriptVariationAllele.pm:1530-1533) before it orders the pair for printing
/// (:1535-1538). The routine holds the `ins` that `_clip_alleles` leaves, its start one
/// above its end (:2245-2248), from its 3' flank to its 5' flank and turns the two unless
/// the end holds a star, so the start is the coordinate printed second, or the one printed
/// first where the 5' flank lies past the stop codon; it holds a `dup` and every other kind
/// in its printed order, one coordinate standing for both. `None` for a pair the routine
/// would have printed the other way round, an `ins` of one coordinate, and a coordinate
/// before the stop codon's last base (`*-5`).
pub fn notation_intron_offsets(notation: &ClippedBody<'_>) -> Option<(String, String)> {
    let (first, last) = (notation.start, notation.last());
    let insertion = notation.is_ins();
    if (insertion && notation.end.is_none()) || first.starts_with("*-") || last.starts_with("*-") {
        return None;
    }
    let (first_parts, last_parts) = (parts_as_read(first)?, parts_as_read(last)?);
    if insertion && last_parts > first_parts && !first.contains('*') {
        return Some((last_parts.1.to_string(), first_parts.1.to_string()));
    }
    if first_parts > last_parts && !last.contains('*') {
        return None;
    }
    Some((first_parts.1.to_string(), last_parts.1.to_string()))
}

/// Whether the plugin's two values beside `notation` are exactly `pair`.
pub fn intron_offsets_are(notation: &ClippedBody<'_>, pair: (&str, &str)) -> bool {
    notation_intron_offsets(notation)
        .is_some_and(|(start, end)| (start.as_str(), end.as_str()) == pair)
}

#[cfg(test)]
#[path = "notation_tests.rs"]
mod tests;
