// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the rules of more than one plugin read: a record's alleles and rows, a row's shape
//! as an insertion or a deletion, its positions on the transcript, and the files of a cache.

use super::entry::{
    allele_is_indel, drops_first_base, entry_kind, field_of, uploaded_alleles, Entry, Row, ABSENT,
};
use super::record::record_label;
use super::Refusal;
use crate::pytext::{self, Json, JsonError};
use flate2::read::MultiGzDecoder;
use regex::Regex;
use std::borrow::Cow;
use std::collections::BTreeSet;
use std::fs;
use std::io::{self, Read};
use std::path::Path;
use std::rc::Rc;
use std::sync::LazyLock;

pub const STRAND_FIELDS: [&str; 2] = ["STRAND", "strand"];
pub const HGVSC_FIELDS: [&str; 2] = ["HGVSc", "hgvsc"];
pub const HGVS_OFFSET_FIELDS: [&str; 2] = ["HGVS_OFFSET", "hgvs_offset"];
pub const SHIFT_LENGTH_FIELDS: [&str; 2] = ["SHIFT_LENGTH", "shift_length"];
pub const CDNA_START_FIELDS: [&str; 2] = ["cDNA_position", "cdna_start"];
pub const CDNA_END_FIELD: &str = "cdna_end";
pub const DISTANCE_FIELDS: [&str; 1] = ["DISTANCE"];
pub const HGVS_INTRON_START_OFFSET_FIELDS: [&str; 1] = ["HGVS_IntronStartOffset"];
pub const HGVS_INTRON_END_OFFSET_FIELDS: [&str; 1] = ["HGVS_IntronEndOffset"];
pub const HGVS_INTRON_OFFSET_FIELDS: [&str; 2] = ["HGVS_IntronStartOffset", "HGVS_IntronEndOffset"];

pub const CODING_SEQUENCE_TERM: &str = "coding_sequence_variant";
pub const FRAMESHIFT_TERM: &str = "frameshift_variant";
pub const STOP_GAINED_TERM: &str = "stop_gained";
pub const STOP_RETAINED_TERM: &str = "stop_retained_variant";
pub const SPLICE_ACCEPTOR_TERM: &str = "splice_acceptor_variant";
pub const SPLICE_SITE_TERMS: [&str; 2] = [SPLICE_ACCEPTOR_TERM, "splice_donor_variant"];

/// The default format's spelling of absence for an empty value.
pub fn dash(s: &str) -> &str {
    if s.is_empty() {
        "-"
    } else {
        s
    }
}

/// The integer a string of digits writes. `str.isdigit` takes digits `int` does not read (a
/// superscript) and `int` refuses a number of more than 4,300 digits: such a text is
/// refused, `what` naming it, and so is a number outside `i128`.
pub fn integer_of(digits: &str, what: &str) -> Result<i128, Refusal> {
    pytext::int(digits).ok_or_else(|| {
        let shown = if digits.chars().count() <= 24 {
            digits.to_string()
        } else {
            digits.chars().take(24).collect::<String>() + "..."
        };
        Refusal(format!("{what} {} is not an integer", pytext::repr(&shown)))
    })
}

/// The body of a regular-expression class that matches what `\d` matches in a Python `str`
/// pattern: the decimal digits of every script.
pub fn re_digit() -> &'static str {
    static CLASS: LazyLock<String> = LazyLock::new(|| {
        (0..=0x10ffffu32)
            .filter_map(char::from_u32)
            .filter(|c| pytext::decimal_value(*c) == Some(0))
            .map(|zero| format!(r"\x{{{:x}}}-\x{{{:x}}}", zero as u32, zero as u32 + 9))
            .collect()
    });
    &CLASS
}

/// Whether `text` is an integer as the rules write one: decimal digits of any script after
/// an optional `-` (`-?\d+` matched whole).
pub fn is_integer(text: &str) -> bool {
    pytext::isdecimal(text.strip_prefix('-').unwrap_or(text))
}

/// Whether `text` is an integer whose first digit is `1` to `9` (`-?[1-9]\d*` matched whole).
pub fn is_non_zero_integer(text: &str) -> bool {
    let digits = text.strip_prefix('-').unwrap_or(text);
    matches!(digits.as_bytes().first(), Some(b'1'..=b'9')) && pytext::isdecimal(digits)
}

/// The coding or non-coding notation of an insertion, a duplication or a deletion after the
/// transcript's name: its first position, its second where it names two, and `dup`, `del`
/// or `ins` with the inserted bases. A position may open with `*` and close with an intron
/// offset. The pattern is anchored at both ends.
pub static INDEL_NOTATION: LazyLock<Regex> = LazyLock::new(|| {
    let digit = re_digit();
    Regex::new(&format!(
        r"^[cn]\.(\*?-?[{digit}]+(?:[+-][{digit}]+)?)(?:_(\*?-?[{digit}]+(?:[+-][{digit}]+)?))?(dup|del|ins[ACGT]+)$"
    ))
    .expect("regex")
});

pub fn is_sequence(allele: &str) -> bool {
    !allele.is_empty() && allele.bytes().all(|b| b"ACGTN".contains(&b))
}

/// The start and the end a Location names after its first `:`, parted at the first `-`;
/// the end is empty where the Location names one position.
pub fn location_ends(location: &str) -> (&str, &str) {
    let span = location.split_once(':').map_or(location, |(_, span)| span);
    let (start, _, end) = pytext::partition(span, "-");
    (start, end)
}

/// Whether an entry is an insertion's or a deletion's: an allele of `-` is one; else the
/// record's own alleles decide where the parser read them (`Entry.indel`); else the Location
/// decides, an indel's spanning another number of bases than the allele has (an insertion as
/// long as its span, or one at a one-base Location, reads as a substitution here, so a row
/// under an identifier may miss the class).
pub fn indel_row(e: &Entry) -> Result<bool, Refusal> {
    if e.allele == "-" {
        return Ok(true);
    }
    if let Some(indel) = e.indel {
        return Ok(indel);
    }
    let Some(location) = e.location.as_deref() else {
        return Ok(false);
    };
    if !is_sequence(&e.allele) {
        return Ok(false);
    }
    let (start, end) = location_ends(location);
    if !(pytext::is_ascii_digits(start) && pytext::is_ascii_digits(end)) {
        return Ok(false);
    }
    let between = integer_of(end, "Location position")? - integer_of(start, "Location position")?;
    Ok(between.saturating_add(1) != e.allele.len() as i128)
}

/// Whether an entry's allele is a sequence and the record's own alleles mark it an indel
/// where the parser read them (`Entry.indel`, true of a deletion spelt as a sequence too),
/// else the Location spans the two bases around an insertion's point. A position outside
/// `i128` is no position.
pub fn insertion_row(e: &Entry) -> bool {
    if !is_sequence(&e.allele) {
        return false;
    }
    if let Some(indel) = e.indel {
        return indel;
    }
    let Some(location) = e.location.as_deref() else {
        return false;
    };
    let (start, end) = location_ends(location);
    if !(pytext::is_ascii_digits(start) && pytext::is_ascii_digits(end)) {
        return false;
    }
    match (pytext::int(start), pytext::int(end)) {
        (Some(start), Some(end)) => start.checked_add(1) == Some(end),
        _ => false,
    }
}

/// The number of positions an entry's Location names, or `None` where it names no span of
/// whole positions; a number past the largest `i128` is read as that one.
pub fn location_span(e: &Entry) -> Result<Option<i128>, Refusal> {
    let Some(location) = e.location.as_deref() else {
        return Ok(None);
    };
    let (start, end) = location_ends(location);
    let end = if end.is_empty() { start } else { end };
    if !(pytext::is_ascii_digits(start) && pytext::is_ascii_digits(end)) {
        return Ok(None);
    }
    let between = integer_of(end, "Location position")? - integer_of(start, "Location position")?;
    let span = between.saturating_add(1);
    Ok((span >= 1).then_some(span))
}

/// The Location both engines print for a span: `chrom:start`, or `chrom:lo-hi` with an
/// insertion's `end = start - 1` written first.
pub fn sequence_location(chrom: &str, start: i128, end: i128) -> String {
    let (lo, hi) = (start.min(end), start.max(end));
    if lo == hi {
        format!("{chrom}:{lo}")
    } else {
        format!("{chrom}:{lo}-{hi}")
    }
}

/// The rows of `record_rows` that are no record entry and stand under `e`'s label.
fn record_members<'a>(e: &Entry, record_rows: &'a [Row]) -> impl Iterator<Item = &'a Row> {
    let label = record_label(e);
    record_rows
        .iter()
        .filter(move |r| entry_kind(r) != "record" && record_label(r) == label)
}

/// The record entry of `e`'s record among `record_rows`: the first under `e`'s label.
pub fn record_entry<'a>(e: &Entry, record_rows: &'a [Row]) -> Option<&'a Row> {
    let label = record_label(e);
    record_rows
        .iter()
        .find(|r| entry_kind(r) == "record" && record_label(r) == label)
}

/// The alleles the rows of `e`'s record carry, its own among them; an entry of
/// `record_rows` under another record's label (`record_label`) is not read.
pub fn record_alleles(e: &Entry, record_rows: &[Row]) -> BTreeSet<String> {
    let mut alleles: BTreeSet<String> = record_members(e, record_rows)
        .map(|r| r.allele.clone())
        .collect();
    alleles.insert(e.allele.clone());
    alleles
}

/// The rows of `e`'s record on its feature under another allele; an entry of `record_rows`
/// under another record's label (`record_label`) is none of them.
pub fn sibling_rows(e: &Entry, record_rows: &[Row]) -> Vec<Row> {
    let (kind, label) = (entry_kind(e), record_label(e));
    record_rows
        .iter()
        .filter(|r| {
            r.allele != e.allele
                && r.feature == e.feature
                && entry_kind(r) == kind
                && record_label(r) == label
        })
        .cloned()
        .collect()
}

/// The chromosome of `e`'s record as its output spells it: the Location's, else the `CHROM`
/// or the `seq_region_name` of the record entry among `record_rows` under `e`'s own label;
/// `None` where none is carried.
pub fn record_chromosome(e: &Entry, record_rows: Option<&[Row]>) -> Option<String> {
    if let Some(location) = e.location.as_deref() {
        let chromosome = pytext::rpartition(location, ":").0;
        return (!chromosome.is_empty()).then(|| chromosome.to_string());
    }
    let record = record_entry(e, record_rows.unwrap_or(&[]))?;
    ["CHROM", "seq_region_name"]
        .iter()
        .find_map(|name| record.fields.get(name))
        .map(str::to_string)
}

/// `indel_row` read from the record's own alleles where the parser did not read them and
/// `record_rows` carries them: the `REF` and `ALT` of the record entry built from the input
/// line of a row under an identifier. The Location alone reads an insertion written at a
/// one-base Location as a substitution.
pub fn record_indel_row(e: &Entry, record_rows: &[Row]) -> Result<bool, Refusal> {
    if e.allele != "-" && e.indel.is_none() {
        let label = record_label(e);
        for r in record_rows {
            if entry_kind(r) != "record" || record_label(r) != label {
                continue;
            }
            let (Some(reference), Some(alternates)) = (r.fields.get("REF"), r.fields.get("ALT"))
            else {
                continue;
            };
            let alternates: Vec<&str> = alternates.split(',').collect();
            if let Some(shape) = allele_is_indel(reference, &alternates, &e.allele, false) {
                return Ok(shape);
            }
        }
    }
    indel_row(e)
}

/// The integer of a text `str.isdigit()` takes; `None` for any other text, and for digits
/// `int()` does not read or a number outside `i128`.
pub fn digits_value(text: &str) -> Option<i128> {
    if pytext::isdigit(text) {
        pytext::int(text)
    } else {
        None
    }
}

fn length(text: &str) -> i128 {
    text.chars().count() as i128
}

/// The lengths of the reference allele and of `e`'s allele as Ensembl VEP hands them to a
/// plugin, an empty allele 0: read from the record's own alleles (an identifier-less
/// `Uploaded_variation`, else the record entry among `record_rows`: its `REF` and `ALT`, or
/// its `allele_string`), else from the row's text where it states them. A deletion's
/// Location spans the deleted bases and a Location of one base or of three or more is the
/// reference allele's; a two-base Location is an insertion's two flanking bases or a
/// two-base reference allele, so `None` there, as for an allele that is no sequence and one
/// that no alternate, or more than one, spells. A Location's span past the largest `i128`
/// is read as that one.
pub fn record_allele_lengths(e: &Entry, record_rows: Option<&[Row]>) -> Option<(i128, i128)> {
    if !is_sequence(&e.allele) && e.allele != "-" {
        return None;
    }
    let mut alleles: Option<(String, Vec<String>)> = None;
    if let Some((reference, alternates)) =
        uploaded_alleles(e.fields.get_or("Uploaded_variation", ""))
    {
        alleles = Some(vcf_alleles_as_spelt(&reference, &alternates));
    } else if let Some(record) = record_entry(e, record_rows.unwrap_or(&[])) {
        if let (Some(reference), Some(alternates)) =
            (record.fields.get("REF"), record.fields.get("ALT"))
        {
            let alternates: Vec<&str> = alternates.split(',').collect();
            alleles = Some(vcf_alleles_as_spelt(reference, &alternates));
        } else if let Some(allele_string) = record.fields.get("allele_string") {
            let mut parts = allele_string.split('/').map(str::to_string);
            let reference = parts.next().unwrap_or_default();
            alleles = Some((reference, parts.collect()));
        }
    }
    if let Some((reference, alternates)) = alleles {
        let named = alternates.iter().filter(|a| **a == e.allele).count();
        return (named == 1).then(|| {
            (
                length(reference.trim_matches('-')),
                length(e.allele.trim_matches('-')),
            )
        });
    }
    let (start, end) = location_ends(e.location.as_deref()?);
    if !(pytext::isdigit(start) && (pytext::isdigit(end) || end.is_empty())) {
        return None;
    }
    let first = pytext::int(start)?;
    let last = if end.is_empty() {
        first
    } else {
        pytext::int(end)?
    };
    let span = (last - first).saturating_add(1);
    if e.allele == "-" {
        return Some((span, 0));
    }
    (span != 2).then(|| (span, length(&e.allele)))
}

/// The start and the end a row's cDNA position writes: the two parts of the text formats'
/// `cDNA_position`, read up to a `/total` suffix, or the `cdna_start` and `cdna_end` keys.
fn cdna_ends(e: &Entry) -> (&str, &str) {
    match e.fields.get("cDNA_position") {
        Some(position) => {
            let (start, _, end) = pytext::partition(pytext::partition(position, "/").0, "-");
            (start, end)
        }
        None => (
            pytext::partition(e.fields.get_or("cdna_start", ABSENT), "/").0,
            pytext::partition(e.fields.get_or(CDNA_END_FIELD, ABSENT), "/").0,
        ),
    }
}

/// The row's first and last cDNA positions, the lower first: the `a-b` or the single
/// position of `cDNA_position`, read up to a `/total` suffix, or `cdna_start` and
/// `cdna_end`. `None` where either is not a number.
pub fn cdna_span(e: &Entry) -> Option<(i128, i128)> {
    let (start, mut end) = cdna_ends(e);
    if e.fields.contains("cDNA_position") && end.is_empty() {
        end = start;
    }
    let (start, end) = (digits_value(start)?, digits_value(end)?);
    Some((start.min(end), start.max(end)))
}

/// Whether the row's cDNA position names two positions, both numbers. A row that begins or
/// ends outside an exon prints `?` or nothing there.
pub fn cdna_span_closed(e: &Entry) -> bool {
    let (start, end) = cdna_ends(e);
    pytext::isdigit(start) && pytext::isdigit(end) && start != end
}

/// The reverse complement of a sequence allele; any other allele unchanged.
pub fn reverse_complement(allele: &str) -> String {
    if allele.is_empty() || !allele.bytes().all(|b| b"ACGTNacgtn-".contains(&b)) {
        return allele.to_string();
    }
    allele
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'a' => 't',
            b'c' => 'g',
            b'g' => 'c',
            b't' => 'a',
            other => char::from(other),
        })
        .collect()
}

/// The entry's sequence allele as its transcript reads it: as written where the row's
/// `STRAND` is `1`, reverse-complemented where it is `-1`. `None` for an allele that is no
/// sequence and for a row that states no strand.
pub fn allele_on_transcript_strand(e: &Entry) -> Option<String> {
    if !is_sequence(&e.allele) {
        return None;
    }
    match field_of(e, &STRAND_FIELDS) {
        "1" => Some(e.allele.clone()),
        "-1" => Some(reverse_complement(&e.allele)),
        _ => None,
    }
}

/// A VCF line's REF and ALTs as its entries spell them: Parser/VCF.pm:300-323 drop the first
/// base of a record that holds an ALT of another length than its REF, a `*` among them, and
/// whose alleles not spelt with `*` all share it, an empty remainder written `-`.
pub fn vcf_alleles_as_spelt<S: AsRef<str>>(
    reference: &str,
    alternates: &[S],
) -> (String, Vec<String>) {
    let as_given = |a: &S| a.as_ref().to_string();
    if !drops_first_base(reference, alternates, false) {
        return (
            reference.to_string(),
            alternates.iter().map(as_given).collect(),
        );
    }
    let rest = |allele: &str| {
        let mut chars = allele.chars();
        chars.next();
        dash(chars.as_str()).to_string()
    };
    let spelt = alternates
        .iter()
        .map(|a| {
            if a.as_ref().contains('*') {
                as_given(a)
            } else {
                rest(a.as_ref())
            }
        })
        .collect();
    (rest(reference), spelt)
}

/// The first ALT of a VCF line as its entries spell the allele (`vcf_alleles_as_spelt`);
/// `None` for a line without an alternate.
pub fn vcf_first_alternate<S: AsRef<str>>(reference: &str, alternates: &[S]) -> Option<String> {
    vcf_alleles_as_spelt(reference, alternates)
        .1
        .into_iter()
        .next()
}

/// The first alternate allele of `e`'s record where an entry states the record's alleles in
/// order, spelt as its rows spell the allele: a row's identifier-less `Uploaded_variation`,
/// else the `ALT` or the `allele_string` of the record entry among `record_rows`; `None`
/// where none does.
pub fn first_alternate(e: &Entry, record_rows: &[Row]) -> Option<String> {
    if let Some((reference, alternates)) =
        uploaded_alleles(e.fields.get_or("Uploaded_variation", ""))
    {
        return vcf_first_alternate(&reference, &alternates);
    }
    let record = record_entry(e, record_rows)?;
    if let Some(alternates) = record.fields.get("ALT") {
        let alternates: Vec<&str> = alternates.split(',').collect();
        return vcf_first_alternate(record.fields.get_or("REF", ""), &alternates);
    }
    record
        .fields
        .get_or("allele_string", "")
        .split('/')
        .nth(1)
        .map(str::to_string)
}

/// The alternate alleles of `e`'s record in the record's own order, spelt as its rows spell
/// the allele, where an entry states them (`first_alternate` names the entries); `None`
/// where none does.
pub fn alternates_in_order(e: &Entry, record_rows: &[Row]) -> Option<Vec<String>> {
    if let Some((reference, alternates)) =
        uploaded_alleles(e.fields.get_or("Uploaded_variation", ""))
    {
        return Some(vcf_alleles_as_spelt(&reference, &alternates).1);
    }
    let record = record_entry(e, record_rows)?;
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

/// The number of bases an HGVSc inserts, duplicates or deletes: the `ins` bases, one for a
/// single position, else the span between two positions without an intron offset or a `*`
/// (the numbering has no position 0); `None` for any other notation, and for a position or
/// a length outside `i128`.
pub fn notation_indel_length(hgvsc: &str) -> Option<i128> {
    let found = INDEL_NOTATION.captures(pytext::partition(hgvsc, ":").2)?;
    let kind = found.get(3)?.as_str();
    if let Some(bases) = kind.strip_prefix("ins") {
        return Some(length(bases));
    }
    let Some(end) = found.get(2) else {
        return Some(1);
    };
    let (start, end) = (found.get(1)?.as_str(), end.as_str());
    if !(is_integer(start) && is_integer(end)) {
        return None;
    }
    let (first, last) = (pytext::int(start)?, pytext::int(end)?);
    if first < 0 && 0 < last {
        last.checked_sub(first)
    } else if last < first {
        // The one is added first: the difference alone lies one below `i128` where the
        // length is the smallest `i128`.
        last.checked_add(1)?.checked_sub(first)
    } else {
        last.checked_sub(first)?.checked_add(1)
    }
}

/// Whether the entry's `SHIFT_LENGTH` is a positive integer: a row `--shift_3prime 1` moved,
/// printed under `--shift_length`. Digits `int` does not read are refused (`integer_of`).
pub fn positive_shift(e: &Entry) -> Result<bool, Refusal> {
    let value = e.fields.first_of(&SHIFT_LENGTH_FIELDS, ABSENT);
    Ok(pytext::isdigit(value) && integer_of(value, "SHIFT_LENGTH")? > 0)
}

/// The items of two list-valued field texts of one format, split at the format's list
/// separator: `&` where the reference's text carries one, as a VCF value does, `,`
/// otherwise; no items for a text carrying both, which neither separator splits into the
/// items its engine wrote.
pub fn list_parts<'a>(expected: &'a str, actual: &'a str) -> (Vec<&'a str>, Vec<&'a str>) {
    if [expected, actual]
        .iter()
        .any(|text| text.contains(',') && text.contains('&'))
    {
        return (Vec::new(), Vec::new());
    }
    let separator = if expected.contains('&') { '&' } else { ',' };
    (
        expected.split(separator).collect(),
        actual.split(separator).collect(),
    )
}

/// An entry without a Location under the Location its record gives its allele
/// (`Entry.allele_location`), the entry a rule over a row's Location reads; any other entry
/// as it is.
pub fn at_its_location(e: &Row) -> Row {
    match (&e.location, &e.allele_location) {
        (None, Some(location)) => Rc::new(Entry {
            location: Some(location.clone()),
            ..Entry::clone(e)
        }),
        _ => Rc::clone(e),
    }
}

/// A Feature without a version: the text less a final `.` and the decimal digits after it,
/// read before a final line feed where the text ends in one (`re.sub(r"\.\d+$", "", ...)`).
pub fn without_version(feature: &str) -> Cow<'_, str> {
    let (body, ending) = match feature.strip_suffix('\n') {
        Some(body) => (body, "\n"),
        None => (feature, ""),
    };
    let stem = body.trim_end_matches(|c| pytext::decimal_value(c).is_some());
    match stem.strip_suffix('.') {
        Some(name) if stem.len() < body.len() => {
            if ending.is_empty() {
                Cow::Borrowed(name)
            } else {
                Cow::Owned(format!("{name}{ending}"))
            }
        }
        _ => Cow::Borrowed(feature),
    }
}

/// The line of the `UnicodeDecodeError` CPython's UTF-8 decoder raises at the first
/// ill-formed sequence of `bytes`, which starts at `at`.
pub fn utf8_error(bytes: &[u8], at: usize) -> String {
    let rest = &bytes[at..];
    let lead = rest[0];
    let continuation = |offset: usize| rest.get(offset).is_some_and(|b| b & 0xc0 == 0x80);
    let second_refused = |lowest: u8, low_lead: u8, high_lead: u8| {
        !continuation(1)
            || if rest[1] < lowest {
                lead == low_lead
            } else {
                lead == high_lead
            }
    };
    // The bytes of the refused part, zero for a sequence the input's end cuts short.
    let (refused, reason) = match lead {
        0xc2..=0xdf if rest.len() < 2 => (0, "unexpected end of data"),
        0xc2..=0xdf => (1, "invalid continuation byte"),
        0xe0..=0xef if rest.len() < 2 => (0, "unexpected end of data"),
        0xe0..=0xef if second_refused(0xa0, 0xe0, 0xed) => (1, "invalid continuation byte"),
        0xe0..=0xef if rest.len() < 3 => (0, "unexpected end of data"),
        0xe0..=0xef => (2, "invalid continuation byte"),
        0xf0..=0xf4 if rest.len() < 2 => (0, "unexpected end of data"),
        0xf0..=0xf4 if second_refused(0x90, 0xf0, 0xf4) => (1, "invalid continuation byte"),
        0xf0..=0xf4 if rest.len() < 3 => (0, "unexpected end of data"),
        0xf0..=0xf4 if !continuation(2) => (2, "invalid continuation byte"),
        0xf0..=0xf4 if rest.len() < 4 => (0, "unexpected end of data"),
        0xf0..=0xf4 => (3, "invalid continuation byte"),
        _ => (1, "invalid start byte"),
    };
    let end = if refused == 0 {
        bytes.len()
    } else {
        at + refused
    };
    if end == at + 1 {
        format!("'utf-8' codec can't decode byte 0x{lead:02x} in position {at}: {reason}")
    } else {
        format!(
            "'utf-8' codec can't decode bytes in position {at}-{}: {reason}",
            end - 1
        )
    }
}

/// `repr()` of a `bytes` value: `b'...'`, in double quotes where the bytes hold a `'` and no
/// `"`, a byte outside printable ASCII as its escape.
pub fn bytes_repr(bytes: &[u8]) -> String {
    let quote = if bytes.contains(&b'\'') && !bytes.contains(&b'"') {
        '"'
    } else {
        '\''
    };
    let shown: String = bytes
        .iter()
        .map(|&b| match b {
            b'\\' => "\\\\".to_string(),
            b'\t' => "\\t".to_string(),
            b'\n' => "\\n".to_string(),
            b'\r' => "\\r".to_string(),
            _ if char::from(b) == quote => format!("\\{quote}"),
            0x20..=0x7e => char::from(b).to_string(),
            _ => format!("\\x{b:02x}"),
        })
        .collect();
    format!("b{quote}{shown}{quote}")
}

/// The reason of a file error as Python's `strerror` words it.
pub fn strerror(error: &io::Error) -> String {
    let shown = error.to_string();
    match shown.rfind(" (os error ") {
        Some(at) if error.raw_os_error().is_some() => shown[..at].to_string(),
        _ => shown,
    }
}

/// The text of `path`, gunzipped when its name ends in `.gz`, decoded as UTF-8, every
/// `\r\n` and lone `\r` read as `\n`. A file that cannot be read, or holds other bytes than
/// UTF-8, is refused as `cannot read <path>: <reason>`.
pub fn read_text(path: &Path) -> Result<String, Refusal> {
    let refuse = |reason: String| Refusal(format!("cannot read {}: {reason}", path.display()));
    let mut bytes = fs::read(path).map_err(|error| refuse(strerror(&error)))?;
    if path.as_os_str().as_encoded_bytes().ends_with(b".gz") && !bytes.is_empty() {
        if !bytes.starts_with(&[0x1f, 0x8b]) {
            return Err(refuse(format!(
                "Not a gzipped file ({})",
                bytes_repr(&bytes[..bytes.len().min(2)])
            )));
        }
        let mut plain = Vec::new();
        MultiGzDecoder::new(bytes.as_slice())
            .read_to_end(&mut plain)
            .map_err(|error| {
                refuse(match error.kind() {
                    io::ErrorKind::UnexpectedEof => {
                        "Compressed file ended before the end-of-stream marker was reached"
                            .to_string()
                    }
                    _ => error.to_string(),
                })
            })?;
        bytes = plain;
    }
    let text = String::from_utf8(bytes).map_err(|error| {
        let at = error.utf8_error().valid_up_to();
        refuse(utf8_error(error.as_bytes(), at))
    })?;
    if text.contains('\r') {
        return Ok(text.replace("\r\n", "\n").replace('\r', "\n"));
    }
    Ok(text)
}

/// The JSON value of a cache shard (`read_text`), refused as `<shard> is not JSON` with the
/// reason where the text is none.
pub fn read_shard(shard: &Path) -> Result<Json, Refusal> {
    pytext::json_loads(&read_text(shard)?).map_err(|error| {
        let reason = match error {
            JsonError::Invalid => "the text is no JSON value",
            JsonError::TooDeep => "the value is nested too deep",
        };
        Refusal(format!("{} is not JSON: {reason}", shard.display()))
    })
}

/// The transcript objects of a parsed shard: a bare array, or a range-keyed map of
/// arrays; a `null` in either is a placeholder and is dropped.
pub fn shard_transcripts(shard: &Path, parsed: Json) -> Result<Vec<Json>, Refusal> {
    let refuse = |what: &str| Refusal(format!("{}{what}", shard.display()));
    let items = match parsed {
        Json::Object(members) => {
            let mut items = Vec::new();
            if !members
                .iter()
                .all(|(_, wrapped)| matches!(wrapped, Json::List(_)))
            {
                return Err(refuse(
                    ": a range-keyed map value is not a transcript array",
                ));
            }
            for (_, wrapped) in members {
                if let Json::List(wrapped) = wrapped {
                    items.extend(wrapped);
                }
            }
            items
        }
        Json::List(items) => items,
        _ => {
            return Err(refuse(
                " is neither a transcript array nor a range-keyed map",
            ))
        }
    };
    let transcripts: Vec<Json> = items
        .into_iter()
        .filter(|item| *item != Json::Null)
        .collect();
    if !transcripts
        .iter()
        .all(|transcript| matches!(transcript, Json::Object(_)))
    {
        return Err(refuse(": a transcript is not an object"));
    }
    Ok(transcripts)
}

/// An integer of a cache transcript, a number in a native cache and the dump's string in a
/// converted one; `None` for any other value, and for one outside `i128`.
pub fn cache_integer(value: Option<&Json>) -> Option<i128> {
    match value? {
        Json::Int(digits) => pytext::int(digits),
        Json::Str(text) if is_integer(text) => pytext::int(text),
        _ => None,
    }
}
