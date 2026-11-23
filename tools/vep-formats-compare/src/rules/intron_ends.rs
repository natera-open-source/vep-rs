// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The intron-end and intron-start insertion classes.

use crate::argparse::{path_text, repr};
use crate::cache;
use crate::consts::{
    CDNA_END_FIELD, CDNA_START_FIELDS, COMPLEMENT, HGVSC_FIELDS, HGVS_OFFSET_FIELDS,
    NON_ZERO_INTEGER_RE_PATTERN, STRAND_FIELDS,
};
use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::fields::python_type;
use crate::parse::vcf::input_line_record;
use crate::pattern::Match;
use crate::pyfmt;
use crate::pyjson::{self, Value};
use crate::report::strerror;
use crate::rules::downstream::field_of;
use crate::rules::shifted_indel::cdna_coordinate_order;
use crate::rules::{whole_number_text, whole_text};
use crate::self_parity::reverse_complement;
use crate::shared::{blank, integer_of, not_an_integer};
use indexmap::{IndexMap, IndexSet};
use std::cell::RefCell;
use std::collections::BTreeSet;
use std::path::Path;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::LazyLock;

pub const CLIPPED_NOTATION_RE_PATTERN: &str = r"(?P<numbering>[cn])\.(?P<start>\*?-?\d+(?:[+-]\d+)?)(?:_(?P<end>\*?-?\d+(?:[+-]\d+)?))?(?:[ACGTN]>(?P<alternate_base>[ACGT])|(?P<kind>delins|ins)(?P<bases>[ACGT]+)|(?P<plain>inv|del|dup))";
pub static CLIPPED_NOTATION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(CLIPPED_NOTATION_RE_PATTERN).expect("the pattern compiles"));

pub const INTRON_END_COPY_RE_PATTERN: &str = r"c\.([1-9]\d*)-1dup";
pub static INTRON_END_COPY_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(INTRON_END_COPY_RE_PATTERN).expect("the pattern compiles"));

pub const INTRON_END_DUPLICATION_RE_PATTERN: &str =
    r"c\.([1-9]\d*)([+-])([1-9]\d*)_([1-9]\d*)-1dup";
pub static INTRON_END_DUPLICATION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(INTRON_END_DUPLICATION_RE_PATTERN).expect("the pattern compiles")
});

pub const INTRON_END_INSERTION_RE_PATTERN: &str = r"c\.([1-9]\d*)-1_([1-9]\d*)ins([A-Z]+)";
pub static INTRON_END_INSERTION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(INTRON_END_INSERTION_RE_PATTERN).expect("the pattern compiles")
});

/// The two pairs of bases that end an intron, of the major kind and of the minor kind, as
/// character codes, on the transcript's strand as an HGVSc spells inserted bases.
///
/// Membership only.
pub const INTRON_END_PAIRS: &[&str] = &["AC", "AG"];

/// The longest intron at whose end the candidate keeps the notation: Ensembl VEP reads an intron
/// of this many bases or fewer with the exons beside it.
pub const SHORT_INTRON_BASES: i64 = 13;

pub const INTRON_START_INSERTION_RE_PATTERN: &str = r"c\.([1-9]\d*)_([1-9]\d*)\+1ins([A-Z]+)";
pub static INTRON_START_INSERTION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(INTRON_START_INSERTION_RE_PATTERN).expect("the pattern compiles")
});

/// The version a stable id may end in, which a row's Feature carries under --transcript_version
/// and a cache's stable id of an Ensembl transcript does not.
pub const VERSION_SUFFIX_RE_PATTERN: &str = r"\.\d+(?=\n?\z)";
pub static VERSION_SUFFIX_RE: LazyLock<fancy_regex::Regex> = LazyLock::new(|| {
    fancy_regex::Regex::new(VERSION_SUFFIX_RE_PATTERN).expect("the pattern compiles")
});

/// A cache shard's file name: the first and last position of its window.
pub const CACHE_SHARD_NAME_RE_PATTERN: &str = r"(\d+)-(\d+)\.json(?:\.gz)?";
pub static CACHE_SHARD_NAME_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(CACHE_SHARD_NAME_RE_PATTERN).expect("the pattern compiles"));

static WHOLE_CLIPPED_NOTATION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(CLIPPED_NOTATION_RE_PATTERN));
static WHOLE_INTRON_END_COPY_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(INTRON_END_COPY_RE_PATTERN));
static WHOLE_INTRON_END_DUPLICATION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(INTRON_END_DUPLICATION_RE_PATTERN));
static WHOLE_INTRON_END_INSERTION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(INTRON_END_INSERTION_RE_PATTERN));
static WHOLE_INTRON_START_INSERTION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(INTRON_START_INSERTION_RE_PATTERN));
static WHOLE_NON_ZERO_INTEGER_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(NON_ZERO_INTEGER_RE_PATTERN));
static WHOLE_CACHE_SHARD_NAME_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(CACHE_SHARD_NAME_RE_PATTERN));
static WHOLE_SEQUENCE_LETTERS_RE: LazyLock<regex::Regex> = LazyLock::new(|| whole_text("[A-Z]+"));

/// `int(text)` of a position. A number 64 bits do not hold ends the run in the words of a
/// text `int` does not read.
pub(crate) fn whole_number(text: &str) -> i64 {
    match whole_number_text(text).parse() {
        Ok(number) => number,
        Err(_) => crate::uncaught!("{}", pyfmt::IntError::Invalid(pyfmt::repr_str(text))),
    }
}

/// A position worked out from others. One that 64 bits do not hold ends the run.
fn position_of(position: i128) -> i64 {
    match i64::try_from(position) {
        Ok(position) => position,
        Err(_) => crate::uncaught!("a position of {position} is past 64 bits"),
    }
}

/// The whole number one below a positive one, both as decimal text of any length.
fn one_less(decimal: &str) -> String {
    let mut digits = decimal.as_bytes().to_vec();
    for digit in digits.iter_mut().rev() {
        if *digit > b'0' {
            *digit -= 1;
            break;
        }
        *digit = b'9';
    }
    let text = String::from_utf8(digits).expect("decimal digits");
    match text.strip_prefix('0') {
        Some(rest) if !rest.is_empty() => rest.to_string(),
        _ => text,
    }
}

/// The entry of `lines` at `ordinal` as the interpreter indexes a list: a negative ordinal
/// counts from the end, and one before the first line ends the run as its `IndexError`
/// does. `None` for an ordinal past the last line.
pub(crate) fn line_at(lines: &[String], ordinal: i64) -> Option<&String> {
    let count = i64::try_from(lines.len()).unwrap_or(i64::MAX);
    if ordinal >= count {
        return None;
    }
    let index = if ordinal < 0 {
        ordinal.checked_add(count).filter(|index| *index >= 0)
    } else {
        Some(ordinal)
    };
    match index.and_then(|index| usize::try_from(index).ok()) {
        Some(index) => lines.get(index),
        None => crate::uncaught!("list index out of range"),
    }
}

/// How far apart, in bases, the two HGVSc notations of an allele running over a transcript's
/// edge can lie: the bases of the allele's span, which bound the bases inside the transcript,
/// plus one base of adjacency on either side. The span is the row's Location; without one (the
/// VCF and JSON formats) and with the start open, the closed cDNA end counts the bases inside;
/// `None` with the end open and no Location.
pub fn edge_row_reach(e: &Entry, side: &str) -> Result<Option<i64>> {
    if let Some(location) = &e.location {
        let span = location
            .split_once(':')
            .map_or(location.as_str(), |(_, span)| span);
        let (start, _, end) = pyfmt::partition(span, "-");
        if pyfmt::is_digit(start) {
            let span_bases = if pyfmt::is_digit(end) {
                let last = integer_of(end, "Location position")?;
                let first = integer_of(start, "Location position")?;
                i128::from(last) - i128::from(first) + 1
            } else {
                1
            };
            return i64::try_from(span_bases + 1)
                .map(Some)
                .map_err(|_| not_an_integer(end, "Location position"));
        }
    }
    if side == "start" {
        let text_position = e
            .fields
            .get(CDNA_START_FIELDS[0])
            .map_or("", String::as_str);
        let mut cdna_end = pyfmt::partition(pyfmt::partition(text_position, "/").0, "-").2;
        if cdna_end.is_empty() {
            cdna_end = e.fields.get(CDNA_END_FIELD).map_or("", String::as_str);
        }
        if pyfmt::is_digit(cdna_end) {
            return integer_of(cdna_end, "cDNA position")?
                .checked_add(1)
                .map(Some)
                .ok_or_else(|| not_an_integer(cdna_end, "cDNA position"));
        }
    }
    Ok(None)
}

/// The `cdna_coordinate_index` of a `CLIPPED_NOTATION_RE` match's one or two
/// coordinates.
pub fn notation_coordinate_indexes(notation: &Match) -> Result<Vec<i64>> {
    [notation.named("start"), notation.named("end")]
        .into_iter()
        .flatten()
        .filter(|coordinate| !coordinate.is_empty())
        .map(cdna_coordinate_index)
        .collect()
}

/// An HGVS cDNA coordinate as a position on one integer scale along the transcript: the 5'
/// UTR `-N` at `-N`, the coding sequence `N` at `N - 1` (`c.-1` and `c.1` are
/// neighbours), the 3' UTR `*N` a billion further, since its distance from the coding
/// sequence is the coding length the notation does not carry, each with its intron offset added
/// (`-49-1` is one base before `-49`, `83+1` one after `83`).
pub fn cdna_coordinate_index(coordinate: &str) -> Result<i64> {
    let (region, position, offset) = cdna_coordinate_order(coordinate)?;
    let (position, offset) = (i128::from(position), i128::from(offset));
    let index = if region < 0 {
        position + offset
    } else if region == 0 {
        position - 1 + offset
    } else {
        1_000_000_000 + position + offset
    };
    i64::try_from(index).map_err(|_| not_an_integer(coordinate, "HGVS coordinate"))
}

/// The `CLIPPED_NOTATION_RE` match of an HGVSc that names the row's transcript (its
/// stable id, with or without a version), or `None` for an empty value, another
/// transcript's notation or a notation of another shape.
pub fn clipped_notation(hgvsc: &str, feature: &str) -> Option<Match> {
    if blank(hgvsc) {
        return None;
    }
    let (name, _, body) = pyfmt::partition(hgvsc, ":");
    if name.split('.').next() != feature.split('.').next() {
        return None;
    }
    WHOLE_CLIPPED_NOTATION_RE
        .captures(body)
        .map(|captures| Match::of(&WHOLE_CLIPPED_NOTATION_RE, &captures))
}

/// The copy form of `hgvsp_intron_end_insertion`: an HGVSc body duplicating the one
/// position before an exon's first coding base `N`, `c.N-1dup`.
pub fn copy_of_the_intron_last_base(body: &str) -> bool {
    WHOLE_INTRON_END_COPY_RE.is_match(body)
}

/// The duplication form of `hgvsp_intron_end_insertion`: an HGVSc body duplicating two or
/// more intron positions that end on the one before an exon's first coding base `N`, its
/// first position counted back from `N` (`c.N-K_N-1dup`, `K` at least 2) or, in the
/// intron's first half, on from the coding base before `N` (`c.M+J_N-1dup`, `M` one
/// below `N`).
pub fn duplication_of_the_intron_last_bases(body: &str) -> bool {
    let Some(m) = WHOLE_INTRON_END_DUPLICATION_RE.captures(body) else {
        return false;
    };
    let (first, sign, distance, last) = (&m[1], &m[2], &m[3], &m[4]);
    if sign == "-" {
        // A distance opens with a digit other than zero, so it is 2 or more unless it is 1.
        return first == last && whole_number_text(distance) != "1";
    }
    let first = whole_number_text(first);
    first == one_less(&whole_number_text(last))
}

/// The bases of an `ins` notation between the position before an exon's first coding base
/// and that base, `c.N-1_Nins`; `None` for any other body.
fn bases_inserted_before_the_exon(body: &str) -> Option<String> {
    let m = WHOLE_INTRON_END_INSERTION_RE.captures(body)?;
    (m[1] == m[2]).then(|| m[3].to_string())
}

/// The last two of `bases`, which are letters of one byte each; fewer where it holds fewer.
fn last_pair(bases: &str) -> &str {
    &bases[bases.len().saturating_sub(2)..]
}

/// The insertion form of `hgvsp_intron_end_insertion`: an HGVSc body inserting, between
/// the position before an exon's first coding base `N` and that base (`c.N-1_Nins`), one
/// base, or bases that do not end in one of `INTRON_END_PAIRS`.
pub fn insertion_that_cannot_end_the_intron(body: &str) -> bool {
    bases_inserted_before_the_exon(body)
        .is_some_and(|bases| bases.len() == 1 || !INTRON_END_PAIRS.contains(&last_pair(&bases)))
}

/// The pair form of `hgvsp_intron_end_insertion`: an HGVSc body inserting, between the
/// position before an exon's first coding base `N` and that base (`c.N-1_Nins`), two or
/// more bases that end in one of `INTRON_END_PAIRS`.
pub fn insertion_that_ends_in_an_intron_end_pair(body: &str) -> bool {
    bases_inserted_before_the_exon(body)
        .is_some_and(|bases| bases.len() >= 2 && INTRON_END_PAIRS.contains(&last_pair(&bases)))
}

/// The form of `hgvsp_intron_end_insertion` an HGVSc body has, each form one helper and
/// one block here, or `None` for any other notation.
pub fn intron_end_insertion_form(body: &str) -> Option<String> {
    if copy_of_the_intron_last_base(body) {
        return Some("copy".to_string());
    }

    if duplication_of_the_intron_last_bases(body) {
        return Some("duplication".to_string());
    }

    if insertion_that_cannot_end_the_intron(body) {
        return Some("insertion".to_string());
    }

    if insertion_that_ends_in_an_intron_end_pair(body) {
        return Some("pair".to_string());
    }

    None
}

/// `hgvsp_intron_end_insertion` on a differing HGVSp: the reference's value is a protein
/// notation and the candidate's is empty; the row's HGVSc is one value on both entries, a
/// notation of the row's transcript (its stable id, with or without a version) of one of the
/// four forms (`intron_end_insertion_form`); neither entry carries an HGVS_OFFSET; and the
/// record is written in a gap where an intron of the row's transcript ends before a coding exon
/// base (`IntronEnds.holds`). A record that holds the insertion with reference bases it
/// shares is placed by the bases an `ins` notation inserts, and a text row that gives such a
/// record's span with a minimised allele by those bases or by the count of positions a `dup`
/// notation duplicates (`IntronEnds.record_gaps`). A VCF or JSON entry is placed from its
/// input line by the same bases or count, however the entry spells its allele
/// (`notation_insertion_gaps`). The gaps come from the candidate's JSON cache: without
/// `intron_ends` nothing is accepted.
pub fn hgvsp_intron_end_insertion(
    e: &Entry,
    a: Option<&Entry>,
    ev: &str,
    av: &str,
    intron_ends: Option<&IntronEnds>,
) -> Result<bool> {
    let (Some(intron_ends), Some(a)) = (intron_ends, a) else {
        return Ok(false);
    };
    if !blank(av) || !ev.contains(":p.") {
        return Ok(false);
    }
    let hgvsc = field_of(Some(e), HGVSC_FIELDS);
    if blank(&hgvsc) || field_of(Some(a), HGVSC_FIELDS) != hgvsc {
        return Ok(false);
    }
    if !blank(&field_of(Some(e), HGVS_OFFSET_FIELDS))
        || !blank(&field_of(Some(a), HGVS_OFFSET_FIELDS))
    {
        return Ok(false);
    }
    let (name, _, body) = pyfmt::partition(&hgvsc, ":");
    if name.split('.').next() != e.feature.split('.').next()
        || intron_end_insertion_form(body).is_none()
    {
        return Ok(false);
    }
    let inserted = notation_inserted_bases(&INTRON_END_INSERTION_RE, body);
    intron_ends.holds(e, &inserted, notation_duplicated_positions(body))
}

/// The form of `hgvsp_intron_start_insertion`: an HGVSc body inserting bases between an
/// exon's last coding base `N` and the first position of the intron after it,
/// `c.N_N+1ins`.
pub fn insertion_after_the_exon_last_base(body: &str) -> bool {
    WHOLE_INTRON_START_INSERTION_RE
        .captures(body)
        .is_some_and(|m| m[1] == m[2])
}

/// `hgvsp_intron_start_insertion` on a differing HGVSp: the reference's value is a protein
/// notation and the candidate's is empty; the row's HGVSc is one value on both entries, a
/// notation of the row's transcript (its stable id, with or without a version) of the class's
/// form (`insertion_after_the_exon_last_base`); the two entries carry no HGVS_OFFSET, or one
/// non-zero integer, the bases the 3' rule moved the record; and the gap the record is written
/// in, moved by that offset as signed, is one where an exon of the row's transcript ends inside
/// the coding region before an intron (`IntronEnds.holds_start`). A text row that carries the
/// offset is the class as well where the gap its Location names is such a gap unmoved: a run
/// under `--shift_genomic 1` prints the Location the 3' rule carried the record to. A VCF or
/// JSON entry is placed from its input line, which is the record as written, by the bases the
/// notation inserts, however the entry spells its allele (`notation_insertion_gaps`), and is
/// held to the moved gap alone. A text row whose record holds the insertion with reference
/// bases it shares carries no offset and is placed by the bases the notation inserts, and a
/// text row that gives such a record's span with a minimised allele by those bases, turned by
/// the offset where the 3' rule moved the record (`IntronEnds.record_gaps`). A run under
/// `--shift_3prime 1` prints a moved record's allele as the rule moved it, so on a row that
/// carries the offset a minimised allele is the notation's bases as they stand or turned back
/// by the offset, its span's gaps read as named or moved, and the bases an ALT of an entry's
/// input line inserts are the notation's as they stand or turned back by the offset. The gaps
/// come from the candidate's JSON cache: without `intron_ends` nothing is accepted.
pub fn hgvsp_intron_start_insertion(
    e: &Entry,
    a: Option<&Entry>,
    ev: &str,
    av: &str,
    intron_ends: Option<&IntronEnds>,
) -> Result<bool> {
    let (Some(intron_ends), Some(a)) = (intron_ends, a) else {
        return Ok(false);
    };
    if !blank(av) || !ev.contains(":p.") {
        return Ok(false);
    }
    let hgvsc = field_of(Some(e), HGVSC_FIELDS);
    if blank(&hgvsc) || field_of(Some(a), HGVSC_FIELDS) != hgvsc {
        return Ok(false);
    }
    let offset = field_of(Some(e), HGVS_OFFSET_FIELDS);
    let candidate_offset = field_of(Some(a), HGVS_OFFSET_FIELDS);
    let mut moved = 0;
    if !(blank(&offset) && blank(&candidate_offset)) {
        if candidate_offset != offset || !WHOLE_NON_ZERO_INTEGER_RE.is_match(&offset) {
            return Ok(false);
        }
        moved = whole_number(&offset);
    }
    let (name, _, body) = pyfmt::partition(&hgvsc, ":");
    if name.split('.').next() != e.feature.split('.').next()
        || !insertion_after_the_exon_last_base(body)
    {
        return Ok(false);
    }
    let inserted = notation_inserted_bases(&INTRON_START_INSERTION_RE, body);
    intron_ends.holds_start(e, moved, &inserted)
}

/// The bases an HGVSc body inserts where it is an `ins` notation of `pattern`
/// (`INTRON_END_INSERTION_RE`, `INTRON_START_INSERTION_RE`); empty for any other body.
pub fn notation_inserted_bases(pattern: &regex::Regex, body: &str) -> String {
    let other;
    let whole_pattern: &regex::Regex = match pattern.as_str() {
        INTRON_END_INSERTION_RE_PATTERN => &WHOLE_INTRON_END_INSERTION_RE,
        INTRON_START_INSERTION_RE_PATTERN => &WHOLE_INTRON_START_INSERTION_RE,
        source => {
            other = whole_text(source);
            &other
        }
    };
    whole_pattern
        .captures(body)
        .and_then(|m| m.get(3))
        .map_or_else(String::new, |bases| bases.as_str().to_string())
}

/// The bases a row's HGVSc inserts, on the transcript's strand, read on the forward strand
/// by the row's `strand` (`1`, `-1`; any other value admits either reading).
pub fn forward_readings(inserted: &str, strand: &str) -> Vec<String> {
    match strand {
        "1" => vec![inserted.to_string()],
        "-1" => vec![reverse_complement(inserted)],
        _ => vec![inserted.to_string(), reverse_complement(inserted)],
    }
}

/// The gaps in which an insertion can lie that a text row's record holds with the reference
/// bases it shares, each the lower of the two positions either side of the gap. The record's
/// reference spans `length` positions from `first` and `allele` is its alternate allele
/// on the forward strand, those bases with the inserted ones; `inserted` are the bases the
/// row's HGVSc inserts, read on the forward strand by `strand` (`forward_readings`). The
/// allele holds as many bases more than the span as the notation inserts. For each count of
/// span bases clipped from the allele's front, the rest clipped from its back, the bases
/// between must be the inserted ones, and the gap's lower position is `first` plus that
/// count less one. A text row states no reference bases, so every clipping that leaves the
/// inserted bases names a gap.
pub fn held_insertion_gaps(
    first: i64,
    length: i64,
    allele: &str,
    inserted: &str,
    strand: &str,
) -> Vec<i64> {
    let allele: Vec<char> = allele.chars().collect();
    let width = inserted.chars().count();
    if width == 0 || length < 1 || allele.len() as i128 - i128::from(length) != width as i128 {
        return Vec::new();
    }
    let readings: Vec<Vec<char>> = forward_readings(inserted, strand)
        .iter()
        .map(|reading| reading.chars().collect())
        .collect();
    // The span is shorter than the allele here, so every count of clipped bases indexes it.
    (0..=allele.len() - width)
        .filter(|front| {
            let between = &allele[*front..front + width];
            readings.iter().any(|reading| reading.as_slice() == between)
        })
        .map(|front| position_of(i128::from(first) + front as i128 - 1))
        .collect()
}

/// Every pure insertion an alternate allele clips to against its record's reference allele,
/// which starts at `first`: the lower of the two positions either side of the gap, and the
/// bases inserted there. For each count of reference bases clipped from the alternate's front,
/// the rest clipped from its back, the clipped bases must be the reference's own and the bases
/// between are the inserted ones. An empty reference allele, an insertion written alone, names
/// the gap before `first`.
pub fn clipped_insertions(first: i64, reference: &str, alternate: &str) -> Vec<(i64, String)> {
    let reference: Vec<char> = reference.chars().collect();
    let alternate: Vec<char> = alternate.chars().collect();
    if alternate.len() <= reference.len() {
        return Vec::new();
    }
    let extra = alternate.len() - reference.len();
    (0..=reference.len())
        .filter(|front| {
            alternate[..*front] == reference[..*front]
                && alternate[front + extra..] == reference[*front..]
        })
        .map(|front| {
            (
                position_of(i128::from(first) + front as i128 - 1),
                alternate[front..front + extra].iter().collect(),
            )
        })
        .collect()
}

/// The gaps an input line's alleles give the insertion a row's HGVSc describes, each the
/// lower of the two positions either side of the gap. `reference` is the record's reference
/// allele from `first` and `alternates` its alternate alleles, on the forward strand; an
/// alternate that clips against the reference to a pure insertion names a gap with its
/// inserted bases (`clipped_insertions`), and the notation claims the gap whose bases are
/// its own: the bases an `ins` notation inserts (`inserted`), read on the forward strand
/// by `strand` (`forward_readings`), as they stand or, on a record the 3' rule moved
/// `moved` positions, turned back by the offset (`turned_back`); or, beside a `dup`
/// notation that states its count, as many bases as the positions duplicated (`duplicated`)
/// on an unmoved record. A `dup` notation that states no count names neither, and claims
/// only an insertion written alone, whose bases follow a reference allele of one base or none.
/// How an entry spells its own allele plays no part.
///
/// `duplicated` is 0 and `moved` is 0 where the caller has none.
pub fn notation_insertion_gaps(
    first: i64,
    reference: &str,
    alternates: &[String],
    inserted: &str,
    strand: &str,
    duplicated: i64,
    moved: i64,
) -> Vec<i64> {
    let wanted: Vec<String> = forward_readings(inserted, strand)
        .into_iter()
        .flat_map(|reading| {
            let turned = turned_back(&reading, moved);
            [reading, turned]
        })
        .collect();
    let reference_bases = pyfmt::len(reference) as i128;
    let mut gaps = Vec::new();
    for alternate in alternates {
        for (gap, bases) in clipped_insertions(first, reference, alternate) {
            let claimed = if !inserted.is_empty() {
                wanted.contains(&bases)
            } else if duplicated > 0 {
                moved == 0 && pyfmt::len(&bases) as i128 == i128::from(duplicated)
            } else {
                reference_bases < 2 && i128::from(gap) == i128::from(first) + reference_bases - 1
            };
            if claimed {
                gaps.push(gap);
            }
        }
    }
    gaps
}

/// The count of positions an HGVSc body duplicates where the notation states it: one for a
/// copy of the intron's last base (`c.N-1dup`), `K` for a duplication counted back from
/// the exon's first coding base (`c.N-K_N-1dup`). 0 for any other body, a duplication that
/// starts in the intron's first half among them, whose count depends on the intron's length.
pub fn notation_duplicated_positions(body: &str) -> i64 {
    if WHOLE_INTRON_END_COPY_RE.is_match(body) {
        return 1;
    }
    match WHOLE_INTRON_END_DUPLICATION_RE.captures(body) {
        // A count 64 bits do not hold is the largest they hold: no allele has as many bases.
        Some(m) if &m[2] == "-" && m[1] == m[4] => {
            whole_number_text(&m[3]).parse().unwrap_or(i64::MAX)
        }
        _ => 0,
    }
}

/// `bases` with its first `turn` bases moved to its end; itself where it holds no more.
fn turned(bases: &str, turn: usize) -> String {
    match bases.char_indices().nth(turn) {
        Some((at, _)) => format!("{}{}", &bases[at..], &bases[..at]),
        None => bases.to_string(),
    }
}

/// The gaps a text row's Location names where its Allele is the insertion alone, each the
/// lower of the two positions either side of the gap: every gap inside and beside the span of
/// `length` positions from `first`, from `first` less one to the span's last position. A
/// record that carries another allele beside the insertion keeps its reference span as the
/// row's Location, and a minimised allele spells the inserted bases alone. `allele` is the
/// insertion the row's HGVSc describes: the bases an `ins` notation inserts (`inserted`),
/// on either strand, as they stand or turned back by the bases the 3' rule moved the record
/// (`moved`), since a run prints a moved record's allele as the rule moved it or as the
/// record writes it; or, for a `dup` notation, which names no bases, as many bases as the
/// positions it duplicates (`duplicated`) on an unmoved record. None for any other allele.
///
/// `duplicated` is 0 and `moved` is 0 where the caller has none.
pub fn lone_insertion_gaps(
    first: i64,
    length: i64,
    allele: &str,
    inserted: &str,
    duplicated: i64,
    moved: i64,
) -> Vec<i64> {
    let sequence = |base: char| COMPLEMENT.iter().any(|(code, _)| *code == base as i64);
    if length < 1 || !allele.chars().all(sequence) {
        return Vec::new();
    }
    if !inserted.is_empty() {
        let turn = usize::try_from(moved.unsigned_abs() % pyfmt::len(inserted) as u64)
            .expect("a count of bases below the length of a text");
        let readings = [allele.to_string(), reverse_complement(allele)];
        let named = readings
            .iter()
            .any(|bases| bases == inserted || turned(bases, turn) == inserted);
        if !named {
            return Vec::new();
        }
    } else if duplicated < 1 || moved != 0 || pyfmt::len(allele) as i128 != i128::from(duplicated) {
        return Vec::new();
    }
    let below = i128::from(first) - 1;
    let mut gaps = Vec::new();
    // A list of more positions than a list can count, or than memory holds, ends the run as
    // the interpreter's list does.
    let Some(count) = length
        .checked_add(1)
        .and_then(|count| usize::try_from(count).ok())
    else {
        crate::uncaught!("Python int too large to convert to C ssize_t")
    };
    if gaps.try_reserve_exact(count).is_err() {
        crate::uncaught!("MemoryError");
    }
    gaps.extend((0..=i128::from(length)).map(|step| position_of(below + step)));
    gaps
}

/// An insertion's bases on the forward strand as the record writes them, from `allele`,
/// the bases as the 3' rule moved them `moved` positions along the chromosome: a run that
/// prints the moved allele turns the inserted bases as it carries them, their first base to the
/// end for each position moved up the chromosome. The allele itself for an unmoved record.
pub fn turned_back(allele: &str, moved: i64) -> String {
    let length = pyfmt::len(allele);
    if length == 0 {
        return String::new();
    }
    let forward = usize::try_from(i128::from(moved).rem_euclid(length as i128))
        .expect("a count of bases below the length of a text");
    turned(allele, length - forward)
}

/// The genomic position of the cDNA position `cdna`, counted along `exons` (in genomic
/// order) in the transcript's direction; `None` for a position outside them.
pub fn cdna_genomic_position(exons: &[(i64, i64)], forward: bool, cdna: i64) -> Option<i64> {
    let mut remaining = i128::from(cdna);
    for index in 0..exons.len() {
        let (start, end) = exons[if forward {
            index
        } else {
            exons.len() - 1 - index
        }];
        let (start, end) = (i128::from(start), i128::from(end));
        let length = end - start + 1;
        if 1 <= remaining && remaining <= length {
            // A position inside the exon, which 64 bits hold.
            return i64::try_from(if forward {
                start + remaining - 1
            } else {
                end - remaining + 1
            })
            .ok();
        }
        remaining -= length;
    }
    None
}

/// `int(value)` of a value of a cache transcript as its decimal text, at any length: `None`
/// for a value `int` refuses with a `TypeError` or a `ValueError`. An infinity ends the run
/// as its `OverflowError` does.
fn whole_number_text_of(value: &Value) -> Option<String> {
    match value {
        Value::Bool(flag) => Some(i64::from(*flag).to_string()),
        Value::Int(number) => Some(number.to_string()),
        Value::BigInt(number) => Some(number.clone()),
        Value::Float(number) => match pyfmt::int_of_f64(*number) {
            Ok(whole) => Some(whole),
            Err(pyfmt::FloatToIntError::Nan) => None,
            Err(error) => crate::uncaught!("{error}"),
        },
        Value::Str(text) => pyfmt::int_of_big(text),
        Value::Null | Value::List(_) | Value::Dict(_) => None,
    }
}

/// Whether `int(value)` of a cache transcript's strand is zero or more, the number read at any
/// length: `None` for a value `int` refuses with a `TypeError` or a `ValueError`. An infinity
/// ends the run as its `OverflowError` does.
fn not_negative(value: &Value) -> Option<bool> {
    whole_number_text_of(value).map(|whole| !whole.starts_with('-'))
}

/// What `coding_transcript_exons` gives of a coding transcript: whether it lies on the
/// forward strand, its exons in genomic order, and the lower and the higher genomic bound of
/// its coding region.
pub type CodingExons = (bool, Vec<(i64, i64)>, i64, i64);

/// Whether a cache transcript lies on the forward strand, its exons in genomic order, and
/// the lower and the higher genomic bound of its coding region, stop codon included: the
/// genomic positions of its `cdna_coding_start` and `cdna_coding_end`
/// (`cdna_genomic_position`). `None` for a transcript without a coding region and for one
/// with a coordinate that is no integer.
///
/// A coordinate 64 bits do not hold is read as one that is no integer, once every coordinate
/// has been read, so an infinite one after it still ends the run.
pub fn coding_transcript_exons(transcript: &Value) -> Option<CodingExons> {
    let Some(members) = transcript.as_dict() else {
        crate::uncaught!(
            "'{}' object has no attribute 'get'",
            python_type(transcript)
        )
    };
    let forward = match members.get("strand") {
        Some(strand) => not_negative(strand)?,
        None => true,
    };
    let mut past_64_bits = false;
    let mut whole_number_of = |value: &Value| -> Option<i64> {
        let number = whole_number_text_of(value)?.parse();
        past_64_bits |= number.is_err();
        Some(number.unwrap_or_default())
    };
    let mut exons = Vec::new();
    if let Some(listed) = members.get("exons").filter(|listed| pyjson::truthy(listed)) {
        // Anything but a list of objects has no `start` to read.
        for exon in listed.as_list()? {
            let exon = exon.as_dict()?;
            let start = whole_number_of(exon.get("start")?)?;
            exons.push((start, whole_number_of(exon.get("end")?)?));
        }
    }
    exons.sort_unstable();
    let coding_start = whole_number_of(members.get("cdna_coding_start")?)?;
    let coding_end = whole_number_of(members.get("cdna_coding_end")?)?;
    if past_64_bits {
        return None;
    }
    let first = cdna_genomic_position(&exons, forward, coding_start);
    let last = cdna_genomic_position(&exons, forward, coding_end);
    let (first, last) = (first?, last?);
    Some((forward, exons, first.min(last), first.max(last)))
}

/// The gaps of one cache transcript in which an intron of more than `SHORT_INTRON_BASES`
/// ends before an exon base inside the coding region, each the lower of the two genomic
/// positions either side of the gap: on the forward strand the intron's last base, the exon's
/// first base following it; on the reverse strand the exon's first base in the transcript's
/// direction, its highest genomic base, the intron's last base following it. Empty where
/// `coding_transcript_exons` reads no coding region.
pub fn transcript_intron_end_gaps(transcript: &Value) -> BTreeSet<i64> {
    let Some((forward, exons, low, high)) = coding_transcript_exons(transcript) else {
        return BTreeSet::new();
    };
    let mut gaps = BTreeSet::new();
    for index in 1..exons.len() {
        let (lower_end, upper_start) = (exons[index - 1].1, exons[index].0);
        if i128::from(upper_start) - i128::from(lower_end) - 1 <= i128::from(SHORT_INTRON_BASES) {
            continue;
        }
        let exon_first = if forward { upper_start } else { lower_end };
        if low <= exon_first && exon_first <= high {
            // The intron holds a base, so the exon above it starts past the lowest number.
            gaps.insert(if forward { upper_start - 1 } else { lower_end });
        }
    }
    gaps
}

/// The gaps of one cache transcript in which an exon ends inside the coding region before an
/// intron of any length, each the lower of the two genomic positions either side of the gap: on
/// the forward strand the exon's last base, the intron's first base following it; on the
/// reverse strand the intron's first base in the transcript's direction, its highest genomic
/// base, the exon's last base in that direction following it. Empty where
/// `coding_transcript_exons` reads no coding region.
pub fn transcript_intron_start_gaps(transcript: &Value) -> BTreeSet<i64> {
    let Some((forward, exons, low, high)) = coding_transcript_exons(transcript) else {
        return BTreeSet::new();
    };
    let mut gaps = BTreeSet::new();
    for index in 1..exons.len() {
        let (lower_end, upper_start) = (exons[index - 1].1, exons[index].0);
        let exon_last = if forward { lower_end } else { upper_start };
        if i128::from(upper_start) - i128::from(lower_end) > 1
            && low <= exon_last
            && exon_last <= high
        {
            gaps.insert(if forward { lower_end } else { upper_start - 1 });
        }
    }
    gaps
}

/// The chromosome and the text after it of a Location, parted at its last `:`.
fn location_span(location: &str) -> (&str, &str) {
    let (chromosome, _, span) = pyfmt::rpartition(location, ":");
    (chromosome, span)
}

/// The shards of one chromosome, each with the first and the last position of its window.
pub type ShardWindows = Vec<(i64, i64, PathBuf)>;

/// The gaps of one transcript: those at an intron's end and those at its start.
pub type IntronGaps = (BTreeSet<i64>, BTreeSet<i64>);

/// The gaps of a candidate JSON cache at the two ends of an intron, by chromosome and
/// transcript: those in which an intron ends before a coding exon base
/// (`transcript_intron_end_gaps`), which `hgvsp_intron_end_insertion` reads, and those in
/// which a coding exon ends before an intron (`transcript_intron_start_gaps`), which
/// `hgvsp_intron_start_insertion` reads. A cache files its transcripts under
/// `transcripts/<chromosome>/<first>-<last>.json` or `.json.gz`; a shard is read when a
/// row asks for a transcript no shard read so far holds, the shard whose window holds the
/// row's position first and the chromosome's others nearest first, and no shard is read
/// twice. `input_lines` are the input's data lines, VCF or Ensembl format, which place a VCF
/// or JSON entry.
#[derive(Debug)]
pub struct IntronEnds {
    pub transcripts: PathBuf,
    pub input_lines: Arc<[String]>,
    /// Per chromosome the unread shards, each with its window. A bound of a window that 64
    /// bits do not hold is the largest they hold.
    pub _unread: RefCell<IndexMap<String, ShardWindows>>,
    /// Per chromosome and stable id the gaps of every transcript read, those at an intron's
    /// end and those at its start.
    pub _gaps: RefCell<IndexMap<(String, String), IntronGaps>>,
}

impl IntronEnds {
    /// `input_lines` is `None` where the caller has no input lines.
    pub fn new(cache: &Path, input_lines: Option<&[String]>) -> Result<Self> {
        Self::sharing(cache, input_lines.unwrap_or_default().into())
    }

    /// `new` on input lines another table holds too: the lines are held once.
    pub fn sharing(cache: &Path, input_lines: Arc<[String]>) -> Result<Self> {
        let transcripts = cache::joined(cache, "transcripts");
        if !transcripts.is_dir() {
            return Err(Error::comparator(format!(
                "--json-cache: {} has no transcripts/ directory",
                path_text(cache)
            )));
        }
        Ok(IntronEnds {
            transcripts,
            input_lines,
            _unread: RefCell::new(IndexMap::new()),
            _gaps: RefCell::new(IndexMap::new()),
        })
    }

    /// Whether one of the gaps the entry's insertion can be written in (`record_gaps`, by
    /// `inserted`, the bases the row's HGVSc inserts, or `duplicated`, the count of
    /// positions it duplicates) is a gap of the entry's transcript in which an intron ends
    /// before a coding exon base.
    ///
    /// `inserted` is empty and `duplicated` is 0 where the caller has none.
    pub fn holds(&self, e: &Entry, inserted: &str, duplicated: i64) -> Result<bool> {
        for (chromosome, position) in self.record_gaps(e, inserted, duplicated, 0) {
            let held =
                self.with_transcript_gaps(&chromosome, &e.feature, i128::from(position), |gaps| {
                    gaps.is_some_and(|(ends, _)| ends.contains(&position))
                })?;
            if held {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Whether one of the gaps the entry's insertion can be written in (`record_gaps`, by
    /// `inserted`, the bases the row's HGVSc inserts), moved `moved` positions along the
    /// chromosome, is a gap of the entry's transcript in which a coding exon ends before an
    /// intron. A text entry's gap is read as named as well: a run that prints the moved
    /// location names the gap the insertion is left in.
    ///
    /// `moved` is 0 and `inserted` is empty where the caller has none.
    pub fn holds_start(&self, e: &Entry, moved: i64, inserted: &str) -> Result<bool> {
        let gaps = self.record_gaps(e, inserted, 0, moved);
        let moves: &[i64] = if e.location.is_none() || moved == 0 {
            &[moved]
        } else {
            &[moved, 0]
        };
        for (chromosome, position) in gaps {
            for step in moves {
                let place = i128::from(position) + i128::from(*step);
                let held = self.with_transcript_gaps(&chromosome, &e.feature, place, |gaps| {
                    gaps.is_some_and(|(_, starts)| {
                        i64::try_from(place).is_ok_and(|place| starts.contains(&place))
                    })
                })?;
                if held {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    /// The chromosome and lower position of each gap the entry's insertion can be written
    /// in. A VCF or JSON entry is placed from its input line by the insertion its HGVSc
    /// describes (`notation_insertion_gaps`), however the entry spells its allele: an ALT of
    /// the line, the ALT less a base every allele of the record shares, a minimised allele or
    /// the allele as the 3' rule moved it. A text entry is placed by its Location: the gap of
    /// an insertion written alone first (`written_gap`); a record that carries another
    /// allele beside the insertion keeps the reference base they share and gives its reference
    /// span, one position or more, where an unmoved allele that holds the span's bases with
    /// the inserted ones names the gaps its clippings leave (`held_insertion_gaps`), as
    /// written or reverse-complemented, since the row spells it on the record's strand, and an
    /// allele that is the insertion alone, a minimised allele, names every gap of the span
    /// (`lone_insertion_gaps`), a span of two adjacent positions among them, which is an
    /// insertion written alone or a reference allele of two bases. `inserted` are the bases
    /// the row's HGVSc inserts, `duplicated` the count of positions it duplicates and
    /// `moved` the positions the 3' rule moved the record, by which a moved allele is read;
    /// the gaps are the record's own.
    ///
    /// `inserted` is empty, `duplicated` is 0 and `moved` is 0 where the caller has none.
    pub fn record_gaps(
        &self,
        e: &Entry,
        inserted: &str,
        duplicated: i64,
        moved: i64,
    ) -> Vec<(String, i64)> {
        let strand = field_of(Some(e), STRAND_FIELDS);
        let each_once = |chromosome: &str, positions: Vec<i64>| -> Vec<(String, i64)> {
            positions
                .into_iter()
                .collect::<IndexSet<i64>>()
                .into_iter()
                .map(|position| (chromosome.to_string(), position))
                .collect()
        };
        let Some(location) = &e.location else {
            let Some((chromosome, first, reference, alternates)) = self.input_record(e) else {
                return Vec::new();
            };
            let positions = notation_insertion_gaps(
                first,
                &reference,
                &alternates,
                inserted,
                &strand,
                duplicated,
                moved,
            );
            return each_once(&chromosome, positions);
        };
        let mut gaps: Vec<(String, i64)> = self.written_gap(e).into_iter().collect();
        let allele = pyfmt::upper(&e.allele);
        let (chromosome, span) = location_span(location);
        let (start, _, end) = pyfmt::partition(span, "-");
        let end = if end.is_empty() { start } else { end };
        if !WHOLE_SEQUENCE_LETTERS_RE.is_match(&allele)
            || !(!chromosome.is_empty() && pyfmt::is_digit(start) && pyfmt::is_digit(end))
        {
            return gaps;
        }
        let first = whole_number(start);
        // A span of 2**63 positions is one position short here; a list of either length is
        // more than a list can count.
        let length = i64::try_from(i128::from(whole_number(end)) - i128::from(first) + 1)
            .unwrap_or(i64::MAX);
        let mut positions =
            lone_insertion_gaps(first, length, &allele, inserted, duplicated, moved);
        if moved == 0 {
            for bases in [allele.clone(), reverse_complement(&allele)] {
                positions.extend(held_insertion_gaps(
                    first, length, &bases, inserted, &strand,
                ));
            }
        }
        gaps.extend(each_once(chromosome, positions));
        gaps
    }

    /// The chromosome, first position, reference allele and alternate alleles of a VCF or
    /// JSON entry's input line (`input_line_record`): the line of the entry's record, or the
    /// Ensembl-format line a VCF entry's own line is written from (`Entry.line`). `None`
    /// for an entry without one and for a line of neither format.
    pub fn input_record(&self, e: &Entry) -> Option<(String, i64, String, Vec<String>)> {
        let ordinal = e.record.or(e.line)?;
        input_line_record(line_at(&self.input_lines, ordinal)?)
    }

    /// The chromosome and the lower of the two positions either side of the gap a text
    /// entry's record is written in: the first position of its Location where that names two
    /// adjacent positions. `None` for any other Location and for a VCF or JSON entry, which
    /// is placed from its input line.
    pub fn written_gap(&self, e: &Entry) -> Option<(String, i64)> {
        let (chromosome, span) = location_span(e.location.as_ref()?);
        let (start, _, end) = pyfmt::partition(span, "-");
        if chromosome.is_empty() || !pyfmt::is_digit(start) || !pyfmt::is_digit(end) {
            return None;
        }
        let last = whole_number(end);
        let first = whole_number(start);
        (i128::from(last) == i128::from(first) + 1).then(|| (chromosome.to_string(), first))
    }

    /// The gaps of the transcript `stable_id` on `chromosome` in which an intron ends
    /// before a coding exon base (`transcript_gaps`).
    pub fn gaps(&self, chromosome: &str, stable_id: &str, position: i64) -> Result<BTreeSet<i64>> {
        Ok(self.transcript_gaps(chromosome, stable_id, position)?.0)
    }

    /// The gaps of the transcript `stable_id` on `chromosome` in which a coding exon ends
    /// before an intron (`transcript_gaps`).
    pub fn start_gaps(
        &self,
        chromosome: &str,
        stable_id: &str,
        position: i64,
    ) -> Result<BTreeSet<i64>> {
        Ok(self.transcript_gaps(chromosome, stable_id, position)?.1)
    }

    /// The gaps of the transcript `stable_id` on `chromosome`, those at an intron's end
    /// and those at an intron's start, none for a transcript the cache does not hold there:
    /// the transcript the cache files under that id or, where it files none, under the id
    /// without a trailing version (a row under `--transcript_version`). The chromosome's
    /// directory is the one of that name or, failing that, of the name without a leading
    /// `chr`; its unread shards are read nearest to `position` first until one holds the
    /// transcript.
    pub fn transcript_gaps(
        &self,
        chromosome: &str,
        stable_id: &str,
        position: i64,
    ) -> Result<IntronGaps> {
        self.with_transcript_gaps(chromosome, stable_id, i128::from(position), |gaps| {
            gaps.cloned().unwrap_or_default()
        })
    }

    /// `transcript_gaps` handed to `read` where they are kept, `None` for a transcript the
    /// cache does not hold; `position` may lie past 64 bits.
    fn with_transcript_gaps<T>(
        &self,
        chromosome: &str,
        stable_id: &str,
        position: i128,
        read: impl FnOnce(Option<&IntronGaps>) -> T,
    ) -> Result<T> {
        let named = cache::joined(&self.transcripts, chromosome).is_dir();
        let directory = if named {
            chromosome
        } else {
            pyfmt::removeprefix(chromosome, "chr")
        };
        if !self._unread.borrow().contains_key(directory) {
            let shards = self.shards(directory);
            self._unread
                .borrow_mut()
                .insert(directory.to_string(), shards);
        }
        let mut key = (directory.to_string(), stable_id.to_string());
        while !self._gaps.borrow().contains_key(&key) {
            let nearest = {
                let mut unread = self._unread.borrow_mut();
                let unread = &mut unread[directory];
                // The first of the shards equally near.
                let Some(index) = (0..unread.len()).min_by_key(|index| {
                    let (first, last, _) = &unread[*index];
                    (i128::from(*first) - position)
                        .max(position - i128::from(*last))
                        .max(0)
                }) else {
                    break;
                };
                unread.remove(index).2
            };
            self.read(directory, &nearest)?;
        }
        let gaps = self._gaps.borrow();
        if !gaps.contains_key(&key) {
            key.1 = VERSION_SUFFIX_RE.replace_all(stable_id, "").into_owned();
        }
        Ok(read(gaps.get(&key)))
    }

    /// The shards of the chromosome directory `directory` with their windows; none for a
    /// chromosome the cache does not hold.
    pub fn shards(&self, directory: &str) -> ShardWindows {
        let path = cache::joined(&self.transcripts, directory);
        if directory.is_empty() || !path.is_dir() {
            return Vec::new();
        }
        if let Err(error) = std::fs::read_dir(&path) {
            crate::uncaught!(
                "[Errno {}] {}: {}",
                error.raw_os_error().unwrap_or(0),
                strerror(&error),
                repr(&path_text(&path))
            );
        }
        // Every name of the pattern ends in one of the two suffixes the listing keeps.
        let mut listed = cache::shards(&self.transcripts, Some(directory)).unwrap_or_default();
        listed.sort_by(|a, b| a.file_name().cmp(&b.file_name()));
        let bound = |digits: &str| pyfmt::int_of(digits).unwrap_or(i64::MAX);
        listed
            .into_iter()
            .filter_map(|shard| {
                let name = shard.file_name()?.to_str()?;
                let m = WHOLE_CACHE_SHARD_NAME_RE.captures(name)?;
                let window = (bound(&m[1]), bound(&m[2]));
                Some((window.0, window.1, shard))
            })
            .collect()
    }

    /// Reads one shard: the gaps of each of its transcripts, kept under `chromosome` and
    /// the transcript's stable id.
    pub fn read(&self, chromosome: &str, shard: &Path) -> Result<()> {
        for transcript in cache::shard_transcripts(shard)? {
            let stable_id = match transcript.get("stable_id") {
                Some(stable_id) if pyjson::truthy(stable_id) => pyjson::py_str(stable_id),
                _ => continue,
            };
            let ends = transcript_intron_end_gaps(&transcript);
            let starts = transcript_intron_start_gaps(&transcript);
            let mut gaps = self._gaps.borrow_mut();
            let kept = gaps.entry((chromosome.to_string(), stable_id)).or_default();
            kept.0.extend(ends);
            kept.1.extend(starts);
        }
        Ok(())
    }
}
