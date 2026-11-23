// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the comparison shares with the tuple comparison (`tools/vep-outputs-compare`).

mod lines;

use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::pyfmt;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;
use std::sync::LazyLock;

/// Ensembl VEP's Uploaded_variation for a VCF record without an ID (OutputFactory.pm:883-931):
/// the chromosome, the position and the record's own alleles joined by `/`.
pub const UPLOADED_ALLELES_RE_PATTERN: &str =
    r"^[^\s\x1c-\x1f]+_[0-9]+_([ACGTNacgtn*-]+(?:/[^/\s\x1c-\x1f]+)+)(?:\n?\z)";
pub static UPLOADED_ALLELES_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(UPLOADED_ALLELES_RE_PATTERN).expect("the pattern compiles"));

pub const UPLOADED_POSITION_RE_PATTERN: &str =
    r"^[^\s\x1c-\x1f]+_([0-9]+)_[ACGTNacgtn*-]+(?:/[^/\s\x1c-\x1f]+)+(?:\n?\z)";
pub static UPLOADED_POSITION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(UPLOADED_POSITION_RE_PATTERN).expect("the pattern compiles")
});

/// The field value of an entry that lacks the field.
pub const ABSENT: &str = "<absent>";

pub const CODING_SEQUENCE_TERM: &str = "coding_sequence_variant";

pub const STOP_RETAINED_TERM: &str = "stop_retained_variant";

pub const SPLICE_REGION_TERM: &str = "splice_region_variant";

pub const FRAMESHIFT_TERM: &str = "frameshift_variant";

pub const STOP_GAINED_TERM: &str = "stop_gained";

/// The terms `frameshift` (Utils/VariationEffect.pm:1541-1567) reads away where `stop_retained`
/// does not hold: the in-frame pair through `_get_codon_alleles` (:1107, :1170) and
/// `protein_altering_variant` (:390).
pub const FRAMESHIFT_READ_TERMS: &[&str] = &[
    "inframe_deletion",
    "inframe_insertion",
    "protein_altering_variant",
];

pub const STOP_LOST_TERM: &str = "stop_lost";

pub const INFRAME_INSERTION_TERM: &str = "inframe_insertion";

/// The four of them beside which `splice_region` returns 0 (Utils/VariationEffect.pm:616-627).
pub const SPLICE_REGION_WITHHOLDING_TERMS: &[&str] = &[
    "splice_acceptor_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_donor_variant",
];

/// `within_intron` (Utils/VariationEffect.pm:629-634) returns the flag `intronic`; the
/// intron itself, its end bases included, is that term's window.
pub const INTRON_TERM: &str = "intron_variant";

/// The six terms `_intron_effects` sets from the differing regions and the candidate
/// writes only where the record's span holds a base of the term's window.
pub const PAST_THE_SPAN_TERMS: &[&str] = &[
    "intron_variant",
    "splice_acceptor_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_donor_variant",
    "splice_polypyrimidine_tract_variant",
];

pub const SPLICE_DONOR_5TH_BASE_TERM: &str = "splice_donor_5th_base_variant";

pub const SPLICE_DONOR_REGION_TERM: &str = "splice_donor_region_variant";

/// A gap past which no window of an intron's end reaches a span or its regions.
pub const _FAR: i64 = 1000000000;

/// The row-aware rules of the stop codon classes, whose terms a row can differ by beside the
/// six: one row can carry both differences.
pub const STOP_CODON_ROW_PAIRS: &[fn(&Entry, &Entry) -> Result<bool>] = &[
    stop_codon_overlap_genomic_span_pair,
    stop_retained_without_stop_codon_test_pair,
    stop_retained_three_base_stop_codon_test_pair,
    stop_codon_removed_without_stop_codon_test_pair,
    stop_codon_overlap_insertion_range_pair,
];

/// The largest coordinate vep-rs holds; a POS past it fails its line in every input format,
/// where the reference numifies the text and annotates the record.
pub const VEP_RS_COORDINATE_CEILING: i64 = 4611686018427387903;

/// The POS spellings vep-rs reads besides a decimal integer: what Rust's `f64` parser accepts,
/// ASCII digits with an optional fraction and exponent and an optional `+`, blanks around.
pub const _VEP_RS_POS_NUMBER_PATTERN: &str = r"\+?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?";
pub static _VEP_RS_POS_NUMBER: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(_VEP_RS_POS_NUMBER_PATTERN).expect("the pattern compiles"));

/// `_VEP_RS_POS_NUMBER` held to the whole of a text.
static VEP_RS_POS_WHOLE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(&format!(r"\A(?:{_VEP_RS_POS_NUMBER_PATTERN})\z"))
        .expect("the pattern compiles")
});

/// The seven terms a row of `intron_rows` holds: `PAST_THE_SPAN_TERMS` and
/// `SPLICE_REGION_TERM`, one bit of a `TermBits` each in this order.
const INTRON_ROW_TERMS: [&str; 7] = [
    "splice_donor_variant",
    "splice_acceptor_variant",
    SPLICE_DONOR_5TH_BASE_TERM,
    SPLICE_DONOR_REGION_TERM,
    SPLICE_POLYPYRIMIDINE_TRACT_TERM,
    INTRON_TERM,
    SPLICE_REGION_TERM,
];
const SPLICE_POLYPYRIMIDINE_TRACT_TERM: &str = "splice_polypyrimidine_tract_variant";
type TermBits = u8;
const DONOR: TermBits = 1;
const ACCEPTOR: TermBits = 1 << 1;
const DONOR_5TH_BASE: TermBits = 1 << 2;
const DONOR_REGION: TermBits = 1 << 3;
const TRACT: TermBits = 1 << 4;
const INTRON: TermBits = 1 << 5;
const REGION: TermBits = 1 << 6;
const WITHHOLDING: TermBits = DONOR | ACCEPTOR | DONOR_5TH_BASE | DONOR_REGION;
/// The names `--terms display` prints for the seven terms, each with the terms that print it
/// (the `display_term` of each term, Utils/Constants.pm): one name for the donor and the
/// acceptor term, one for the fifth base, the donor region, the tract and
/// `splice_region_variant`, one for the intron term. No other term prints one of the three.
const DISPLAY_ROW_TERMS: [(&str, TermBits); 3] = [
    ("ESSENTIAL_SPLICE_SITE", DONOR | ACCEPTOR),
    (
        "SPLICE_SITE",
        DONOR_5TH_BASE | DONOR_REGION | TRACT | REGION,
    ),
    ("INTRONIC", INTRON),
];

/// Whether `text` is one or more of the ASCII `letters` and nothing else.
fn spelt_in(text: &str, letters: &str) -> bool {
    !text.is_empty() && text.bytes().all(|b| letters.as_bytes().contains(&b))
}

/// The start and the end a Location writes after its first `:`, parted at the first `-`; the
/// end is empty where the Location names one position.
fn location_bounds(location: &str) -> (&str, &str) {
    let coords = location
        .split_once(':')
        .map_or(location, |(_, coords)| coords);
    coords.split_once('-').unwrap_or((coords, ""))
}

/// The number a run of the digits `0` to `9` writes plus `addend`, as decimal text without
/// leading zeros, for a run of any length `int` converts; a longer run ends the run in the
/// words of `int`'s `ValueError`.
fn plus(digits: &str, addend: usize) -> String {
    if digits.len() > pyfmt::INT_MAX_STR_DIGITS {
        crate::uncaught!("{}", pyfmt::IntError::TooManyDigits(digits.len()));
    }
    let mut carry = addend as u128;
    let mut sum: Vec<u8> = Vec::with_capacity(digits.len() + 1);
    for digit in digits.bytes().rev() {
        carry += u128::from(digit - b'0');
        sum.push(b'0' + (carry % 10) as u8);
        carry /= 10;
    }
    while carry > 0 {
        sum.push(b'0' + (carry % 10) as u8);
        carry /= 10;
    }
    while sum.len() > 1 && sum.last() == Some(&b'0') {
        sum.pop();
    }
    sum.reverse();
    String::from_utf8(sum).expect("decimal digits")
}

/// Whether a Location names two neighbouring positions, each written in the digits `0` to
/// `9`.
fn neighbouring_positions(location: &str) -> bool {
    let (start, end) = location_bounds(location);
    _ten_digits(start) && _ten_digits(end) && {
        let end = plus(end, 0);
        end == plus(start, 1)
    }
}

/// The number of positions an entry's Location names (`location_span`), held past `i64`.
fn positions_named(e: &Entry) -> Result<Option<i128>> {
    let Some(location) = e.location.as_deref() else {
        return Ok(None);
    };
    let (start, end) = location_bounds(location);
    let end = if end.is_empty() { start } else { end };
    if !(_ten_digits(start) && _ten_digits(end)) {
        return Ok(None);
    }
    let end = i128::from(integer_of(end, "Location position")?);
    let start = i128::from(integer_of(start, "Location position")?);
    let span = end - start + 1;
    Ok((span >= 1).then_some(span))
}

/// Whether Ensembl VEP drops the first base of a record's alleles (`Parser/VCF.pm` lines
/// 300-323): they are not already spelt without it (`chopped`), the reference allele holds a
/// base, an alternate has another length than it, and every alternate not spelt with a `*`
/// starts with its first base.
fn drops_shared_first_base(reference: &str, alternates: &[&str], chopped: bool) -> bool {
    let length = pyfmt::len(reference);
    !chopped
        && !reference.is_empty()
        && alternates.iter().any(|a| pyfmt::len(a) != length)
        && alternates
            .iter()
            .filter(|a| !a.contains('*'))
            .all(|a| pyfmt::prefix(a, 1) == pyfmt::prefix(reference, 1))
}

/// An allele without its first base.
fn without_first(allele: &str) -> &str {
    pyfmt::slice(allele, Some(1), None)
}

/// Whether `spelt` is `head`, `written` and `tail` one after the other.
fn spells(spelt: &[char], head: &[char], written: &[char], tail: &[char]) -> bool {
    spelt.len() == head.len() + written.len() + tail.len()
        && spelt.starts_with(head)
        && spelt[head.len()..].starts_with(written)
        && spelt.ends_with(tail)
}

/// The comparator's error for digits `int` does not read or that lie past `i64`.
pub(crate) fn not_an_integer(digits: &str, what: &str) -> Error {
    let shown = if pyfmt::len(digits) <= 24 {
        digits.to_string()
    } else {
        format!("{}...", pyfmt::prefix(digits, 24))
    };
    Error::comparator(format!(
        "{what} {} is not an integer",
        pyfmt::repr_str(&shown)
    ))
}

/// The terms of a Consequence as the pair rules read them: `consequence_set` parted at `,`,
/// so a Consequence without a term reads as the one empty term.
fn term_set(consequence: &str) -> BTreeSet<&str> {
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

/// The terms the reference's entry holds, those it alone holds and those the candidate's
/// alone holds.
fn term_differences<'a>(
    e: &'a Entry,
    a: &'a Entry,
) -> (BTreeSet<&'a str>, BTreeSet<&'a str>, BTreeSet<&'a str>) {
    let (expected, actual) = (term_set(&e.consequence), term_set(&a.consequence));
    let gone = expected.difference(&actual).copied().collect();
    let come = actual.difference(&expected).copied().collect();
    (expected, gone, come)
}

/// Whether `terms` holds `members` and nothing else.
fn is(terms: &BTreeSet<&str>, members: &[&str]) -> bool {
    terms.len() == members.len() && members.iter().all(|member| terms.contains(member))
}

/// Whether an entry is a Transcript's: the Feature_type its key ends in after the last `|`,
/// read up to a `#`, goes to the transcript list (`feature_kind`).
fn transcript_entry(e: &Entry) -> bool {
    let feature_type = e.key.rsplit('|').next().unwrap_or("");
    let feature_type = feature_type.split('#').next().unwrap_or("");
    feature_kind(feature_type) == "transcript"
}

/// The terms of the seven among `terms`.
fn intron_row_bits_of(terms: &BTreeSet<&str>) -> TermBits {
    INTRON_ROW_TERMS
        .iter()
        .enumerate()
        .filter(|(_, term)| terms.contains(*term))
        .fold(0, |bits, (bit, _)| bits | 1 << bit)
}

/// The display names of the seven among `terms`, one bit a name in the order of
/// `DISPLAY_ROW_TERMS`.
fn display_row_names(terms: &BTreeSet<&str>) -> u8 {
    DISPLAY_ROW_TERMS
        .iter()
        .enumerate()
        .filter(|(_, (name, _))| terms.contains(name))
        .fold(0, |names, (bit, _)| names | 1 << bit)
}

/// The display names a row of the seven prints, as `display_row_names` reads them.
fn displayed(row: TermBits) -> u8 {
    DISPLAY_ROW_TERMS
        .iter()
        .enumerate()
        .filter(|(_, (_, members))| row & members != 0)
        .fold(0, |names, (bit, _)| names | 1 << bit)
}

/// The terms `bits` holds.
fn intron_row_terms(bits: TermBits) -> BTreeSet<String> {
    INTRON_ROW_TERMS
        .iter()
        .enumerate()
        .filter(|(bit, _)| bits & (1 << bit) != 0)
        .map(|(_, term)| term.to_string())
        .collect()
}

/// `intron_rows` with each row as its bits.
fn intron_row_bits(
    regions: &[(i64, i64)],
    span: i64,
    start: i64,
    end: i64,
    forward: bool,
    wide_exons: bool,
) -> (TermBits, TermBits) {
    let last = span - 1;
    let in_list = last >= start - 3 && end + 3 >= 0;
    let at_an_end = (last >= start - 3 && start + 7 >= 0) || (last >= end - 7 && end + 3 >= 0);
    let frameshift = end - start <= 12;
    let read = || {
        let outside = move |&(from, to): &(i64, i64)| !(frameshift && to >= start && from <= end);
        regions.iter().copied().filter(outside)
    };
    let met =
        |first: i64, last_base: i64| read().any(|(from, to)| to >= first && from <= last_base);
    let held = |first: i64, last_base: i64| last >= first && last_base >= 0;
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
    for (term, first, last_base) in windows {
        let read_here = if term == TRACT {
            (in_list || at_an_end) && off_the_exons
        } else {
            at_an_end
        };
        if read_here && met(first, last_base) {
            set_by_a_region |= term;
        }
        if held(first, last_base) {
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
                .any(|&(first, last_base)| to >= first && from <= last_base);
            in_any |= in_one;
            at_a_site |= (to >= start && from <= start + 1) || (to >= end - 1 && from <= end);
            if !at_a_site {
                last_assigned = in_one;
            }
        }
    }
    let row = |mut terms: TermBits, region: bool| {
        if terms & DONOR_5TH_BASE != 0 {
            terms &= !DONOR_REGION;
        }
        if region && terms & WITHHOLDING == 0 {
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
/// Such an allele changes no base of its span, and the candidate reads it as the insertion it
/// writes, between the span's last base and the base after it.
fn insertion_after_the_span(regions: &[(i64, i64)], span: i64) -> bool {
    span >= 2 && matches!(regions, [(first, _)] if *first == span)
}

/// Whether `allele` opens with its reference allele `reference` of one base and holds
/// further bases. The reference does not compare such a pair of alleles and reads the one
/// base as changed (VariationFeatureOverlapAllele.pm:418, :452-454), so the regions do not
/// show the shape. Such an allele changes no base of its span either, and the candidate
/// reads it as the insertion it writes, between that base and the base after it.
fn insertion_after_one_base(reference: &str, allele: &str) -> bool {
    let mut bases = reference.chars();
    let (Some(base), None) = (bases.next(), bases.next()) else {
        return false;
    };
    let mut letters = allele.chars();
    letters.next() == Some(base) && letters.next().is_some()
}

/// The splice and intron terms of the candidate's row for an allele
/// `insertion_after_the_span` takes, on the intron `intron_rows` reads: those of the
/// insertion between the span's last base, the offset `span - 1`, and the base after it, read
/// by those two neighbours as `_intron_effects` reads an insertion
/// (BaseTranscriptVariationAllele.pm:126).
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
) -> TermBits {
    let (before, after) = (span - 1, span);
    let in_list = after >= start - 3 && before <= end + 3;
    let at_an_end =
        (after >= start - 3 && before <= start + 7) || (after >= end - 7 && before <= end + 3);
    let frameshift = end - start <= 12;
    let holds_both = |first: i64, last_base: i64| before >= first && after <= last_base;
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
    for (term, first, last_base) in windows {
        if at_an_end && holds_both(first, last_base) {
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
        .any(|&(first, last_base)| holds_both(first, last_base))
        || after == start
        || before == end
        || after == start + 2
        || before == end - 2;
    if at_an_end && region && terms & WITHHOLDING == 0 {
        terms |= REGION;
    }
    terms
}

/// Calls `visit` with the places of `intron_placements`, in its order, until it returns
/// true, and returns whether it did.
fn visit_intron_placements(
    span: i64,
    reach: i64,
    positions: (&str, &str, &str, &str),
    mut visit: impl FnMut(bool, i64, i64) -> bool,
) -> Result<bool> {
    let last = span - 1;
    let before: Vec<i64> = (0..18).map(|gap| -gap).chain([-_FAR]).collect();
    let gaps: BTreeSet<i64> = (0..42).chain(42.max(reach - 20)..reach + 19).collect();
    let after: Vec<i64> = gaps.iter().map(|gap| last + gap).chain([_FAR]).collect();
    let (numbered_first, numbered_last) = (_ten_digits(positions.0), _ten_digits(positions.1));
    let width = if numbered_first && numbered_last {
        let end = i128::from(integer_of(positions.1, "cDNA position")?);
        let start = i128::from(integer_of(positions.0, "cDNA position")?);
        Some(end - start + 1)
    } else {
        None
    };
    for forward in [true, false] {
        let (first_in_sequence, last_in_sequence) = if forward {
            (numbered_first, numbered_last)
        } else {
            (numbered_last, numbered_first)
        };
        if let Some(width) = width {
            if width == i128::from(span) {
                for start in (1..8).map(|gap| last + gap) {
                    let far_ends = after.iter().copied().filter(|end| *end >= start);
                    let near_ends = (0..42).map(|length| start + length);
                    if far_ends
                        .chain(near_ends)
                        .any(|end| visit(forward, start, end))
                    {
                        return Ok(true);
                    }
                }
                if (-16..0).any(|start| (start..0).any(|end| visit(forward, start, end))) {
                    return Ok(true);
                }
            } else if 0 < width && width < i128::from(span) {
                let length = span - width as i64;
                if (1..last - length + 1).any(|start| visit(forward, start, start + length - 1)) {
                    return Ok(true);
                }
            }
        } else if first_in_sequence {
            if (1..span).any(|start| after.iter().any(|end| visit(forward, start, *end))) {
                return Ok(true);
            }
        } else if last_in_sequence {
            if (0..last).any(|end| before.iter().any(|start| visit(forward, *start, end))) {
                return Ok(true);
            }
        } else if before
            .iter()
            .any(|start| after.iter().any(|end| visit(forward, *start, *end)))
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether `text` is written in the digits `0` to `9` alone, the digits a variant
/// file writes a coordinate in; `str.isdigit` takes the digits of other scripts as well.
pub fn _ten_digits(text: &str) -> bool {
    !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit())
}

/// The default format's spelling of absence for an empty VCF or JSON value.
pub fn dash(s: &str) -> String {
    if s.is_empty() { "-" } else { s }.to_string()
}

/// Whether `allele` is an insertion or a deletion of the record whose reference allele
/// is `reference` and alternate alleles `alternates`: the alternate the allele names and
/// the reference differ in length, a difference no trimming of shared bases changes. A record
/// with one alternate not spelt with a `*` has that alternate's shape on every row, since
/// Ensembl VEP shortens a lone allele past the shared base (`minimise_alleles`: `TCT` to
/// `TCTC` is written `C`) and the allele then spells neither the REF nor the ALT. With
/// several, Ensembl VEP drops the record's shared leading base only when every alternate not
/// spelt with a `*` shares it (`Parser/VCF.pm` lines 306-321 exempt `<DEL:*>` as they
/// exempt `*`), so `A` to `AT,G` keeps `AT` as the insertion's allele; a VCF line's
/// REF and ALT are read the same way, and the parts of a JSON `allele_string`, already so
/// spelt (`chopped`), as they are, since dropping a base again from `T/TT/-` leaves no
/// allele spelt `TT`. `None` when the allele is not a sequence (a symbolic or a breakend
/// allele), when no alternate names it, or when two alternates of different shapes name it.
pub fn allele_is_indel(
    reference: &str,
    alternates: &[&str],
    allele: &str,
    chopped: bool,
) -> Option<bool> {
    if !spelt_in(allele, "ACGTN-") {
        return None;
    }
    let bases = |text: &str| pyfmt::len(pyfmt::strip_chars(text, "-"));
    let trimmed: Vec<&str> = alternates
        .iter()
        .copied()
        .filter(|a| !a.contains('*'))
        .collect();
    if let [only] = trimmed[..] {
        return Some(bases(reference) != bases(only));
    }
    let trim = drops_shared_first_base(reference, alternates, chopped);
    let reference = if trim {
        without_first(reference)
    } else {
        reference
    };
    let mut shapes = BTreeSet::new();
    for alternate in alternates {
        let alt = if trim && !alternate.contains('*') {
            without_first(alternate)
        } else {
            alternate
        };
        if (if alt.is_empty() { "-" } else { alt }) == allele {
            shapes.insert(bases(reference) != bases(alt));
        }
    }
    if shapes.len() == 1 {
        shapes.pop_first()
    } else {
        None
    }
}

/// Whether the reference allele of a record is empty once the record's shared first base is
/// dropped, `reference` and `alternates` read as `allele_is_indel` reads them
/// (`Parser/VCF.pm` lines 300-323 drop the base where every alternate not spelt with a `*`
/// shares it and one ALT has another length than the REF). Every sequence allele of such a
/// record is an insertion, which Ensembl VEP stores with a start one above its end. A
/// reference allele that keeps a base past the shared first one is not empty here, a lone
/// alternate Ensembl VEP shortens further (`trim_sequences`) included.
pub fn reference_allele_is_empty(reference: &str, alternates: &[&str], chopped: bool) -> bool {
    let kept = if drops_shared_first_base(reference, alternates, chopped) {
        without_first(reference)
    } else {
        reference
    };
    pyfmt::strip_chars(kept, "-").is_empty()
}

/// The reference and alternate alleles an identifier-less `Uploaded_variation` spells
/// (`21_10542414_A/T/AG`), or `None` for an identifier, which names no allele.
pub fn uploaded_alleles(uploaded: &str) -> Option<(String, Vec<String>)> {
    let spelt = UPLOADED_ALLELES_RE.captures(uploaded)?.get(1)?.as_str();
    let mut alleles = spelt.split('/').map(str::to_string);
    let reference = alleles.next()?;
    Some((reference, alleles.collect()))
}

/// The POS an identifier-less `Uploaded_variation` spells before its alleles
/// (`21_10542414_A/T/AG`), or `None` for an identifier.
pub fn uploaded_position(uploaded: &str) -> Result<Option<i64>> {
    match UPLOADED_POSITION_RE
        .captures(uploaded)
        .and_then(|found| found.get(1))
    {
        Some(position) => integer_of(position.as_str(), "Uploaded_variation position").map(Some),
        None => Ok(None),
    }
}

/// The reference bases an output row's allele stands in for, read by position from a
/// record's own alleles. The record whose reference allele `reference` starts at
/// `position` writes `allele` at `location` where one of its alternates is the
/// reference allele with the Location's bases replaced by the allele (`-` replaces them by
/// nothing): those bases are returned. Where an alternate is the reference allele with the
/// allele inserted between the Location's two neighbouring positions the row is a pure
/// insertion's and the empty string is returned. `None` where the record writes no such
/// row, where two alternates read it differently, and for an allele or a Location that is
/// not a sequence's.
pub fn replaced_reference(
    position: i64,
    reference: &str,
    alternates: &[&str],
    location: Option<&str>,
    allele: &str,
) -> Result<Option<String>> {
    let Some(location) = location else {
        return Ok(None);
    };
    if allele != "-" && !spelt_in(allele, "ACGTN") {
        return Ok(None);
    }
    let (start, end) = location_bounds(location);
    let end = if end.is_empty() { start } else { end };
    if !(_ten_digits(start) && _ten_digits(end)) {
        return Ok(None);
    }
    let position = i128::from(position);
    let first = i128::from(integer_of(start, "Location position")?) - position;
    let last = i128::from(integer_of(end, "Location position")?) - position;
    let written: Vec<char> = if allele == "-" {
        Vec::new()
    } else {
        allele.chars().collect()
    };
    let reference: Vec<char> = pyfmt::upper(reference).chars().collect();
    let length = reference.len() as i128;
    let mut found: BTreeSet<String> = BTreeSet::new();
    for alternate in alternates {
        let spelt: Vec<char> = pyfmt::upper(alternate).chars().collect();
        if 0 <= first && first <= last && last < length {
            let (first, last) = (first as usize, last as usize);
            if spells(
                &spelt,
                &reference[..first],
                &written,
                &reference[last + 1..],
            ) {
                found.insert(reference[first..=last].iter().collect());
            }
        }
        if !written.is_empty() && last == first + 1 && 0 <= first && first < length {
            let after = first as usize + 1;
            if spells(&spelt, &reference[..after], &written, &reference[after..]) {
                found.insert(String::new());
            }
        }
    }
    Ok(if found.len() == 1 {
        found.pop_first()
    } else {
        None
    })
}

/// The count of leading letters two alleles share.
fn shared_prefix(one: &[char], other: &[char]) -> usize {
    one.iter()
        .zip(other)
        .take_while(|(one, other)| one == other)
        .count()
}

/// The count of trailing letters two alleles share.
fn shared_suffix(one: &[char], other: &[char]) -> usize {
    one.iter()
        .rev()
        .zip(other.iter().rev())
        .take_while(|(one, other)| one == other)
        .count()
}

/// `bases` turned by `places`: the bases from the place `places` on, then the bases before
/// it, the place counted around the run's length in either direction.
fn turned(bases: &[char], places: i128) -> Vec<char> {
    if bases.is_empty() {
        return Vec::new();
    }
    let from = places.rem_euclid(bases.len() as i128) as usize;
    bases[from..]
        .iter()
        .chain(&bases[..from])
        .copied()
        .collect()
}

/// The reference bases an output row's allele stands in for where the row is a record's
/// insertion or deletion that the 3' rule moved. A run under `--shift_genomic 1` writes
/// such a row's Location as the record's span moved by the count of bases the rule moved
/// the allele, on along the chromosome for a transcript of the forward strand and back for
/// one of the reverse strand (OutputFactory.pm:1711-1715; an intergenic row's on along it,
/// :2001-2006), and the moved allele as the row's allele (:1349-1352). That Location is no
/// span of the record's alleles, and `replaced_reference` reads none for it.
///
/// The record whose reference allele `reference` starts at `position` writes `allele` at
/// `location` so where `allele` is `-` and one of its alternates is the reference allele
/// less a run of as many bases as the Location names. A deletion is moved only past
/// bases that repeat it, so the bases at the Location are the run turned by the count of
/// places it lies from the run as written, whichever way: those bases are returned. Where
/// the Location names two neighbouring positions, an alternate is the reference allele
/// with bases added before one of its bases or after its last, and `allele` is the added
/// bases turned by the count of places moved, the row is a pure insertion's and the empty
/// string is returned.
///
/// The record's own bases bound the move. A Location before the reference allele is read
/// only where the alternate lacks the run, or adds the bases, at the allele's first place,
/// and one past the allele only where it does so at the allele's last place: at any other
/// place a base of the record that does not repeat the run or the added bases stands in
/// the way. `None` for such a Location; for a Location on the reference allele's own
/// positions, every position of a deletion's and the first of an insertion's, where
/// `replaced_reference` reads what the record's alleles place and no row is a moved one;
/// where no alternate writes the row so or two write it differently; and for an allele or
/// a Location that is not a sequence's.
pub fn moved_reference(
    position: i64,
    reference: &str,
    alternates: &[&str],
    location: Option<&str>,
    allele: &str,
) -> Result<Option<String>> {
    let Some(location) = location else {
        return Ok(None);
    };
    if allele != "-" && !spelt_in(allele, "ACGTN") {
        return Ok(None);
    }
    let (start, end) = location_bounds(location);
    let end = if end.is_empty() { start } else { end };
    if !(_ten_digits(start) && _ten_digits(end)) {
        return Ok(None);
    }
    let position = i128::from(position);
    let first = i128::from(integer_of(start, "Location position")?) - position;
    let last = i128::from(integer_of(end, "Location position")?) - position;
    let written: Vec<char> = if allele == "-" {
        Vec::new()
    } else {
        allele.chars().collect()
    };
    let reference: Vec<char> = pyfmt::upper(reference).chars().collect();
    // A Location on the reference allele's own positions, every position of a deletion's
    // and the first of an insertion's, is `replaced_reference`'s to read.
    let on_the_allele = if written.is_empty() { last } else { first };
    if 0 <= first && on_the_allele < reference.len() as i128 {
        return Ok(None);
    }
    // Any other Location starts before the reference allele or reaches past it.
    let before_the_allele = first < 0;
    let mut found: BTreeSet<String> = BTreeSet::new();
    for alternate in alternates {
        let spelt: Vec<char> = pyfmt::upper(alternate).chars().collect();
        let before = shared_prefix(&reference, &spelt);
        let after = shared_suffix(&reference, &spelt);
        if written.is_empty() {
            let run = last - first + 1;
            if run < 1 || reference.len() as i128 - spelt.len() as i128 != run {
                continue;
            }
            // The places of the reference allele at which the alternate lacks the run,
            // the last of them the place at which the run closes the allele.
            let lowest = spelt.len().saturating_sub(after);
            let highest = before.min(spelt.len());
            let unbounded = if before_the_allele {
                lowest == 0
            } else {
                highest == spelt.len()
            };
            if lowest > highest || !unbounded {
                continue;
            }
            let run = &reference[lowest..lowest + run as usize];
            found.insert(turned(run, first - lowest as i128).into_iter().collect());
        } else if last == first + 1 && spelt.len() == reference.len() + written.len() {
            // The places of the reference allele before which the alternate adds the
            // bases, the last of them the place after the allele's last base: the
            // Location's first position is the base before the added bases.
            let lowest = reference.len().saturating_sub(after);
            let highest = before.min(reference.len());
            let unbounded = if before_the_allele {
                lowest == 0
            } else {
                highest == reference.len()
            };
            if lowest > highest || !unbounded {
                continue;
            }
            let added = &spelt[lowest..lowest + written.len()];
            if turned(added, first - (lowest as i128 - 1)) == written {
                found.insert(String::new());
            }
        }
    }
    Ok(if found.len() == 1 {
        found.pop_first()
    } else {
        None
    })
}

/// Whether a field value is empty in any format's spelling: `-`, the empty string or
/// absent from the entry.
pub fn blank(value: &str) -> bool {
    matches!(value, "" | "-" | ABSENT)
}

/// The integer a string of digits writes. `str.isdigit` takes digits `int` does not
/// read (a superscript) and `int` refuses a number longer than the interpreter converts:
/// either is a malformed input and ends the run with the comparator's error. So does a
/// number `i64` does not hold; `integer_text_of` reads one of any length.
pub fn integer_of(digits: &str, what: &str) -> Result<i64> {
    pyfmt::int_of(digits).ok_or_else(|| not_an_integer(digits, what))
}

/// The integer a string of digits writes, as its decimal text: `integer_of` for a number of
/// any length `int` converts.
pub fn integer_text_of(digits: &str, what: &str) -> Result<String> {
    pyfmt::int_of_big(digits).ok_or_else(|| not_an_integer(digits, what))
}

/// Consequence terms as a sorted, de-duplicated, comma-joined set.
pub fn consequence_set(s: &str) -> String {
    let terms: BTreeSet<&str> = s
        .split(['&', ','])
        .filter(|term| !term.is_empty())
        .collect();
    terms.into_iter().collect::<Vec<_>>().join(",")
}

/// `JSON.pm:348-349`: the consequence list a feature type's entries go to
/// (`Transcript` to `transcript`, `RegulatoryFeature` to `regulatory_feature`, an
/// undefined type to `intergenic`); a list name is returned unchanged.
pub fn feature_kind(feature_type: &str) -> String {
    let kind = pyfmt::lower(if matches!(feature_type, "" | "-") {
        "intergenic"
    } else {
        feature_type
    });
    if kind == "intergenic" || kind.contains("_feature") {
        kind
    } else {
        kind.replace("feature", "_feature")
    }
}

/// Whether an entry is an insertion's or a deletion's: an allele of `-` is one; else the
/// record's own alleles decide where the parser read them (`Entry.indel`: the VCF and JSON
/// entries, and a text row whose Uploaded_variation spells the alleles); else the Location
/// decides, an indel's spanning more bases than the allele has (an insertion's two flanking
/// bases against a one-base allele, a deletion's span against the kept bases; an insertion as
/// long as its span, or an insertion at a one-base Location, reads as a substitution here, so
/// a row under an identifier may miss the class).
pub fn indel_row(e: &Entry) -> Result<bool> {
    if e.allele == "-" {
        return Ok(true);
    }
    if let Some(indel) = e.indel {
        return Ok(indel);
    }
    let Some(location) = e.location.as_deref() else {
        return Ok(false);
    };
    if !spelt_in(&e.allele, "ACGTN") {
        return Ok(false);
    }
    let (start, end) = location_bounds(location);
    if !(_ten_digits(start) && _ten_digits(end)) {
        return Ok(false);
    }
    let end = i128::from(integer_of(end, "Location position")?);
    let start = i128::from(integer_of(start, "Location position")?);
    Ok(end - start + 1 != e.allele.len() as i128)
}

/// Whether the entry carries the start and the end of a transcript position: the two bounds
/// of the text formats' `column` (`?` for an end that maps to no transcript base, one number
/// standing for both, `-` for neither, each read up to the `/total` suffix `--total_length`
/// appends, OutputFactory.pm:1735-1748), or the JSON `start_key` and `end_key`, which an
/// unmapped end leaves out.
pub fn position_bounds(e: &Entry, column: &str, start_key: &str, end_key: &str) -> (bool, bool) {
    let Some(position) = e.fields.get(column) else {
        let carried = |key: &str| !blank(e.fields.get(key).map_or(ABSENT, String::as_str));
        return (carried(start_key), carried(end_key));
    };
    let position = position
        .split_once('/')
        .map_or(position.as_str(), |(head, _)| head);
    if blank(position) {
        return (false, false);
    }
    match position.split_once('-') {
        Some((start, end)) => (start != "?", end != "?"),
        None => (position != "?", position != "?"),
    }
}

/// Whether the entry carries the coordinates the stop codon test edits at
/// (Utils/VariationEffect.pm:1507-1508, `return 0 unless $cdna_start && $cdna_end &&
/// $cds_start`): both cDNA positions and the CDS start. A CDS end alone, the `?-N` of a
/// span that begins in the 5' UTR, is not one of them.
pub fn stop_codon_test_coordinates(e: &Entry) -> bool {
    position_bounds(e, "cDNA_position", "cdna_start", "cdna_end") == (true, true)
        && position_bounds(e, "CDS_position", "cds_start", "cds_end").0
}

/// Whether the entry shows that `_ins_del_stop_altered_cil` (Utils/VariationEffect.pm
/// :1488-1539) returns before it edits anything: the allele holds a letter other than `A`,
/// `C`, `G` and `T` (:1498), the span lacks a cDNA position at one end, or it has no CDS
/// start (:1508, `stop_codon_test_coordinates`).
pub fn stop_codon_test_not_made(e: &Entry) -> bool {
    (e.allele != "-" && !spelt_in(&e.allele, "ACGT")) || !stop_codon_test_coordinates(e)
}

/// Whether the reference entry `e` and the candidate entry `a` of one key are a row on
/// which Ensembl VEP's genomic stop codon overlap misses a stop base that ends an exon: a
/// Transcript entry of an insertion or deletion with both cDNA positions and a CDS start
/// on both sides, `splice_region_variant` in both sets of terms, and the two sets equal but
/// for the reference's `coding_sequence_variant` against the candidate's
/// `stop_retained_variant`.
pub fn stop_codon_overlap_genomic_span_pair(e: &Entry, a: &Entry) -> Result<bool> {
    if !transcript_entry(e) || !indel_row(e)? {
        return Ok(false);
    }
    let (expected, gone, come) = term_differences(e, a);
    Ok(is(&gone, &[CODING_SEQUENCE_TERM])
        && is(&come, &[STOP_RETAINED_TERM])
        && expected.contains(SPLICE_REGION_TERM)
        && !stop_codon_test_not_made(e)
        && !stop_codon_test_not_made(a))
}

/// Whether an entry is an insertion's or a deletion's as `indel_row` reads one, an allele
/// of other nucleotide letters than `A`, `C`, `G`, `T` and `N` included: where the
/// record's own alleles say nothing, such an allele against a Location of another length.
pub fn length_changing_row(e: &Entry) -> Result<bool> {
    if indel_row(e)? {
        return Ok(true);
    }
    let Some(location) = e.location.as_deref() else {
        return Ok(false);
    };
    let letters = !e.allele.is_empty() && e.allele.bytes().all(|b| b.is_ascii_alphabetic());
    if e.indel.is_some() || !letters {
        return Ok(false);
    }
    let (start, end) = location_bounds(location);
    Ok(_ten_digits(start) && _ten_digits(end) && {
        let end = plus(end, 1);
        end != plus(start, e.allele.len())
    })
}

/// Whether the reference entry `e` and the candidate entry `a` of one key are a row on
/// which Ensembl VEP writes `stop_retained_variant` for a stop codon test it does not make:
/// a Transcript entry of an insertion or deletion, the reference's terms holding
/// `stop_retained_variant` and the candidate's being the same without it, with nothing,
/// `coding_sequence_variant` or `frameshift_variant` in its place, and both entries
/// showing that the test is not made.
pub fn stop_retained_without_stop_codon_test_pair(e: &Entry, a: &Entry) -> Result<bool> {
    if !transcript_entry(e) || !length_changing_row(e)? {
        return Ok(false);
    }
    let (_, gone, come) = term_differences(e, a);
    Ok(is(&gone, &[STOP_RETAINED_TERM])
        && (come.is_empty() || is(&come, &[CODING_SEQUENCE_TERM]) || is(&come, &[FRAMESHIFT_TERM]))
        && stop_codon_test_not_made(e)
        && stop_codon_test_not_made(a))
}

/// Whether the reference entry `e` and the candidate entry `a` of one key are a row on
/// which Ensembl VEP writes `stop_retained_variant` from the three bases at the old stop
/// position alone: a Transcript entry of an insertion or deletion with both cDNA positions, a
/// CDS start and an allele of bases on both sides, the reference's terms holding
/// `stop_retained_variant` and the candidate's `coding_sequence_variant` in its place, or,
/// where the reference's hold it alone or beside one of `inframe_insertion`,
/// `inframe_deletion` and `protein_altering_variant`, `frameshift_variant` with or
/// without `stop_gained` in the place of those terms.
pub fn stop_retained_three_base_stop_codon_test_pair(e: &Entry, a: &Entry) -> Result<bool> {
    if !transcript_entry(e) || !length_changing_row(e)? {
        return Ok(false);
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
    Ok(terms_stand && !stop_codon_test_not_made(e) && !stop_codon_test_not_made(a))
}

/// Whether the reference entry `e` and the candidate entry `a` of one key are a row on
/// which Ensembl VEP writes `stop_retained_variant` for a stop codon the record removes with
/// no test made: a Transcript entry of an insertion or deletion, both entries lacking a cDNA
/// position at one end or a CDS start, and the two sets of terms equal but for the
/// reference's `stop_retained_variant` against the candidate's `stop_lost`.
pub fn stop_codon_removed_without_stop_codon_test_pair(e: &Entry, a: &Entry) -> Result<bool> {
    if !transcript_entry(e) || !length_changing_row(e)? {
        return Ok(false);
    }
    let (_, gone, come) = term_differences(e, a);
    Ok(is(&gone, &[STOP_RETAINED_TERM])
        && is(&come, &[STOP_LOST_TERM])
        && !stop_codon_test_coordinates(e)
        && !stop_codon_test_coordinates(a))
}

/// Whether an entry's allele is a sequence (`-`, a symbolic or a breakend allele is
/// refused) and the record's own alleles mark it an indel where the parser read them
/// (`Entry.indel`, true of a deletion spelt as a sequence too, the `A` row of `AT` to
/// `A,G`, whose `del` notation `clamped_insertion_position` and `shifted_insertion_end`
/// both refuse), else the Location spans the two bases around an insertion's point.
///
/// A VCF or JSON entry and a default or tab row whose first column spells the alleles
/// answer from the record's alleles, alike in each form. A row under an identifier carries
/// none and answers from its Location, the limit of that form: two bases replaced by any
/// sequence, a deletion of one of two bases spelt as a sequence among them, share an
/// insertion's Location and are taken, and an insertion kept with the base it shares with
/// its reference, whose Location is that one position, is not. A rule that counts the bases
/// a row adds reads the alleles itself where an entry carries them
/// (`added_bases_written_outside_the_transcript`).
pub fn insertion_row(e: &Entry) -> bool {
    if !spelt_in(&e.allele, "ACGTN") {
        return false;
    }
    match e.indel {
        Some(indel) => indel,
        None => e.location.as_deref().is_some_and(neighbouring_positions),
    }
}

/// Whether an entry `insertion_row` takes is an insertion between two neighbouring
/// positions, the allele Ensembl VEP stores with a start one above its end and the only one
/// whose start `_overlaps_stop_codon_cil` moves (Utils/VariationEffect.pm:1408-1415). Its
/// Location, where the format writes one, is two neighbouring positions, and the first of
/// these that applies decides:
///
/// - its reference bases (`Entry.reference`) are empty: a pure insertion's entry, whatever
///   the allele's length and however the record pads its alleles;
/// - the parser read the record's own alleles: its reference allele is empty once the
///   shared first base is dropped (`Entry.insertion`), which refuses a deletion spelt as a
///   sequence, a replacement by a sequence of another length, an insertion kept with its
///   shared first base and an allele none of them names beside a reference allele that
///   keeps a base;
/// - the records of the input give the entry bases that are not empty: a replacement's
///   entry;
/// - no record gives the entry its bases: it has a Location and its allele has another
///   length than two, which refuses the two-base substitution of the same Location.
///
/// Under an identifier without the records of the input, two bases replaced by a sequence
/// of another length have an insertion's Location and allele and are taken with the
/// insertions, and an insertion of two bases is refused. The records give an entry its
/// bases by Location and allele: an insertion no record places is refused where a record
/// replaces the same two bases by the same allele.
pub fn inserted_between_neighbours(e: &Entry) -> bool {
    if e.location
        .as_deref()
        .is_some_and(|location| !neighbouring_positions(location))
    {
        return false;
    }
    match (e.reference.as_deref(), e.insertion) {
        (Some(""), _) => true,
        (_, Some(insertion)) => insertion,
        (Some(_), None) => false,
        (None, None) => e.location.is_some() && pyfmt::len(&e.allele) != 2,
    }
}

/// Whether the reference entry `e` and the candidate entry `a` of one key are a row on
/// which Ensembl VEP's genomic stop codon overlap misses an insertion that keeps the stop
/// codon: a Transcript entry of an insertion between two neighbouring positions
/// (`inserted_between_neighbours`) with both cDNA positions and a CDS start on
/// both sides, the two sets of terms equal but for the reference's
/// `coding_sequence_variant` against the candidate's `stop_retained_variant`, or the
/// reference's `frameshift_variant` with or without `stop_gained` against the candidate's
/// `stop_retained_variant` with or without `inframe_insertion`. Beside
/// `splice_region_variant` the first pair is `stop_codon_overlap_genomic_span_pair`'s where
/// the allele is of `A`, `C`, `G` and `T` alone, and this predicate's where it holds
/// another letter, the rows that one refuses: `peptide` returns nothing for such an allele
/// (TranscriptVariationAllele.pm:701), so `stop_retained` takes its second arm
/// (Utils/VariationEffect.pm:1318) and the moved start fails the overlap there too.
pub fn stop_codon_overlap_insertion_range_pair(e: &Entry, a: &Entry) -> Result<bool> {
    if !transcript_entry(e) || !insertion_row(e) || !inserted_between_neighbours(e) {
        return Ok(false);
    }
    let (expected, gone, come) = term_differences(e, a);
    let terms_stand = if is(&gone, &[CODING_SEQUENCE_TERM]) {
        let unknown_base = [e, a]
            .iter()
            .any(|entry| entry.allele != "-" && !spelt_in(&entry.allele, "ACGT"));
        is(&come, &[STOP_RETAINED_TERM]) && (!expected.contains(SPLICE_REGION_TERM) || unknown_base)
    } else if is(&gone, &[FRAMESHIFT_TERM]) || is(&gone, &[STOP_GAINED_TERM, FRAMESHIFT_TERM]) {
        is(&come, &[STOP_RETAINED_TERM]) || is(&come, &[INFRAME_INSERTION_TERM, STOP_RETAINED_TERM])
    } else {
        false
    };
    Ok(terms_stand && stop_codon_test_coordinates(e) && stop_codon_test_coordinates(a))
}

/// The cDNA start and end and the CDS start and end an entry carries, as written (`?`
/// for an end that maps to no transcript base, `-` for a position the entry lacks), or
/// `None` where a position is not of that form. A text format's column is read up to the
/// `/total` suffix `--total_length` appends, one number standing for both ends; the JSON
/// keys are read as they are.
pub fn transcript_positions(e: &Entry) -> Option<(String, String, String, String)> {
    let mut written: Vec<String> = Vec::with_capacity(4);
    for (column, start_key, end_key) in [
        ("cDNA_position", "cdna_start", "cdna_end"),
        ("CDS_position", "cds_start", "cds_end"),
    ] {
        let pair = match e.fields.get(column).map(String::as_str) {
            None => {
                [start_key, end_key].map(|key| dash(e.fields.get(key).map_or("-", String::as_str)))
            }
            Some("" | "-") => ["-".to_string(), "-".to_string()],
            Some(position) => {
                let head = position.split_once('/').map_or(position, |(head, _)| head);
                let (start, end) = head.split_once('-').unwrap_or((head, head));
                [start.to_string(), end.to_string()]
            }
        };
        if !pair
            .iter()
            .all(|part| matches!(part.as_str(), "?" | "-") || _ten_digits(part))
        {
            return None;
        }
        written.extend(pair);
    }
    let [cdna_start, cdna_end, cds_start, cds_end]: [String; 4] = written.try_into().ok()?;
    Some((cdna_start, cdna_end, cds_start, cds_end))
}

/// The number of positions an entry's Location names, or `None` where it names no span
/// of whole positions.
pub fn location_span(e: &Entry) -> Result<Option<i64>> {
    Ok(positions_named(e)?.map(|span| i64::try_from(span).unwrap_or(i64::MAX)))
}

/// Whether an entry shows, without its record, a span of three or more bases that lies in
/// the transcript's sequence from end to end under an allele of more nucleotide letters than
/// the span has bases: its Location names the span, and its two cDNA positions are numbers
/// as far apart as the span's two ends, which a span that holds an intron's base does not
/// have. A Location of two neighbouring positions is refused: an insertion between them is
/// written with the same Location, allele and positions as two replaced bases, and only
/// the record tells the two apart.
pub fn longer_allele_inside_one_exon(e: &Entry) -> Result<bool> {
    let Some(span) = positions_named(e)? else {
        return Ok(false);
    };
    if span < 3 || !spelt_in(&e.allele, "ACGTN") || e.allele.len() as i128 <= span {
        return Ok(false);
    }
    let Some((cdna_start, cdna_end, _, _)) = transcript_positions(e) else {
        return Ok(false);
    };
    if !(_ten_digits(&cdna_start) && _ten_digits(&cdna_end)) {
        return Ok(false);
    }
    let end = i128::from(integer_of(&cdna_end, "cDNA position")?);
    let start = i128::from(integer_of(&cdna_start, "cDNA position")?);
    Ok(end - start + 1 == span)
}

/// The regions `_get_differing_regions` (VariationFeatureOverlapAllele.pm:405-460)
/// cuts a record into, as offsets from the span's first base: where both alleles hold more
/// than one base (:418) the runs of indices at which they differ, every index past the
/// shorter allele's end differing (:428-450); else one region, the span (:452-454).
pub fn differing_regions(reference: &str, allele: &str) -> Vec<(i64, i64)> {
    let reference: Vec<char> = reference.chars().collect();
    let allele: Vec<char> = allele.chars().collect();
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

/// The splice and intron terms of the reference's row and of the candidate's for a record
/// whose span is the offsets `0` to `span - 1`, on a transcript of the given strand whose
/// one intron near the span runs from offset `start` to offset `end`.
///
/// The reference's row is `_intron_effects` (BaseTranscriptVariationAllele.pm:99-224) on
/// `regions`, each placed with no bound at the span's end (:117-118). The intron is read
/// where the span overlaps it widened by three bases, and its two ends where the span
/// overlaps the eleven bases around each (BaseTranscriptVariation.pm:1027-1029); an intron
/// of thirteen bases or fewer is a frameshift intron (:1032-1035) and a region that
/// overlaps it sets nothing (BaseTranscriptVariationAllele.pm:138-141, :171-174). A region
/// in the intron's interior sets the intron term (:143-150) and one in a window of the five
/// splice terms its flag (:152-161, :176-206). The flag of `splice_region` is assigned for
/// each region in turn, whether the region lies in one of its four windows, until a region has
/// met the intron's first two or last two bases (:209-216; Utils/VariationEffect.pm:91-94),
/// so the last region assigned decides it. The tract term is written only
/// off the exons (Constants.pm:971-983), which a transcript with a frameshift intron reads
/// twelve bases wider (BaseTranscriptVariation.pm:861; `wide_exons` says another intron
/// of the transcript is one). The donor region is withheld beside the fifth base and
/// `splice_region_variant` beside the first four terms (Utils/VariationEffect.pm:597-627).
///
/// The candidate's row keeps a term of the six only where the span itself holds a base of
/// a window of that term, the intron being the intron term's window, reads
/// `splice_region` from any region in one of its four windows, and withholds the donor
/// region and `splice_region_variant` by its own row's terms. For an allele that opens with
/// its whole reference allele (`insertion_after_the_span`) it is the row of the insertion the
/// allele writes (`insertion_after_the_span_row`).
pub fn intron_rows(
    regions: &[(i64, i64)],
    span: i64,
    start: i64,
    end: i64,
    forward: bool,
    wide_exons: bool,
) -> (BTreeSet<String>, BTreeSet<String>) {
    let (reference, candidate) = intron_row_bits(regions, span, start, end, forward, wide_exons);
    (intron_row_terms(reference), intron_row_terms(candidate))
}

/// The places `(forward, start, end)` of one intron beside a span of `span` bases that
/// an entry's cDNA positions allow, as offsets from the span's first base, for regions that
/// run `reach` bases past the span's end. A numbered cDNA position is a base of the
/// transcript's sequence and any other a base of the intron; the first position is the
/// span's first base on the forward strand and its last on the reverse. Two numbers as far
/// apart as the span put the intron after the span, one to seven bases on, as far as one of
/// its ends is still read, or before it, where an intron shorter than its windows has
/// them run on past its last base; two numbers closer than the span put it inside the
/// span. A gap wider than any window, and than the regions' reach, is placed once.
pub fn intron_placements(
    span: i64,
    reach: i64,
    positions: (&str, &str, &str, &str),
) -> Result<Vec<(bool, i64, i64)>> {
    let mut places = Vec::new();
    visit_intron_placements(span, reach, positions, |forward, start, end| {
        places.push((forward, start, end));
        false
    })?;
    Ok(places)
}

/// A place of the transcript's intron, `(forward, start, end, wide_exons)` as
/// `intron_rows` reads it, at which the record's regions give the reference entry's
/// splice and intron terms and the span, or the insertion written after it, gives the
/// candidate entry's, or `None` where the entry's cDNA positions allow no such place. The
/// exons are read at their own width or, for a tract term a row lacks, twelve bases wider.
/// `reference` is the reference bases the allele stands in for; the places are those of
/// `intron_placements`. For an allele that opens with its one reference base
/// (`insertion_after_one_base`) the reference's row is that of the base read as changed, and
/// the candidate's the row of the insertion after it (`insertion_after_the_span_row`).
pub fn record_span_placement(
    e: &Entry,
    a: &Entry,
    reference: &str,
) -> Result<Option<(bool, i64, i64, bool)>> {
    let Some(positions) = transcript_positions(e) else {
        return Ok(None);
    };
    let rows = (
        intron_row_bits_of(&term_set(&e.consequence)),
        intron_row_bits_of(&term_set(&a.consequence)),
    );
    let regions = differing_regions(reference, &e.allele);
    let span = pyfmt::len(reference) as i64;
    let reach = pyfmt::len(&e.allele) as i64 - span;
    let after_one_base = insertion_after_one_base(reference, &e.allele);
    let positions = (
        positions.0.as_str(),
        positions.1.as_str(),
        positions.2.as_str(),
        positions.3.as_str(),
    );
    let mut place = None;
    visit_intron_placements(span, reach, positions, |forward, start, end| {
        let read = |wide_exons: bool| {
            placed_row_bits(
                &regions,
                span,
                after_one_base,
                (forward, start, end, wide_exons),
            )
        };
        let narrow = read(false);
        let wide_exons = if narrow == rows {
            false
        } else if ((narrow.0 & TRACT != 0 && rows.0 & TRACT == 0)
            || (narrow.1 & TRACT != 0 && rows.1 & TRACT == 0))
            && read(true) == rows
        {
            true
        } else {
            return false;
        };
        place = Some((forward, start, end, wide_exons));
        true
    })?;
    Ok(place)
}

/// The reference's row and the candidate's at one place `(forward, start, end, wide_exons)`
/// of the intron: `intron_row_bits`, with the candidate's row that of the insertion after the
/// base (`insertion_after_the_span_row`) for an allele that opens with its one reference base.
fn placed_row_bits(
    regions: &[(i64, i64)],
    span: i64,
    after_one_base: bool,
    (forward, start, end, wide_exons): (bool, i64, i64, bool),
) -> (TermBits, TermBits) {
    let (reference_row, span_row) = intron_row_bits(regions, span, start, end, forward, wide_exons);
    if after_one_base {
        let insertion = insertion_after_the_span_row(span, start, end, forward, wide_exons);
        (reference_row, insertion)
    } else {
        (reference_row, span_row)
    }
}

/// Whether the cDNA positions `positions` allow a place of the transcript's intron at which
/// the reference's splice and intron terms print the display names `expected` and the
/// candidate's the names `actual` (`displayed`). The exons are read twelve bases wider as
/// well where a row at their own width holds the tract term, whose name three other terms
/// print.
fn displayed_span_placement(
    positions: (&str, &str, &str, &str),
    reference: &str,
    allele: &str,
    expected: u8,
    actual: u8,
) -> Result<bool> {
    let regions = differing_regions(reference, allele);
    let span = pyfmt::len(reference) as i64;
    let reach = pyfmt::len(allele) as i64 - span;
    let after_one_base = insertion_after_one_base(reference, allele);
    let shown = |(reference_row, candidate_row): (TermBits, TermBits)| {
        (displayed(reference_row), displayed(candidate_row)) == (expected, actual)
    };
    visit_intron_placements(span, reach, positions, |forward, start, end| {
        let read = |wide_exons: bool| {
            placed_row_bits(
                &regions,
                span,
                after_one_base,
                (forward, start, end, wide_exons),
            )
        };
        let narrow = read(false);
        shown(narrow) || ((narrow.0 | narrow.1) & TRACT != 0 && shown(read(true)))
    })
}

/// The class on two entries printed under `--terms display`, or `None` where neither holds
/// a display name of the seven (`DISPLAY_ROW_TERMS`). The entries' other terms are equal and
/// their names of the seven differ; the candidate's entry alone does not hold the name of the
/// donor and the acceptor term, so its site term at an intron of one base stays a difference
/// here; the record's alleles show a longer alternate allele that is no pure insertion, as
/// under those names, and a place of one intron gives both entries' names
/// (`displayed_span_placement`). Where the allele is not written after its span, the
/// reference's entry holds a name the candidate's lacks, as it holds a term of the six alone
/// under those names. Without the record's alleles no such entry is taken: three names stand
/// for seven terms, and the rows alone do not tell them apart.
fn display_vocabulary_pair(
    e: &Entry,
    a: &Entry,
    positions: (&str, &str, &str, &str),
    expected_terms: &BTreeSet<&str>,
    actual_terms: &BTreeSet<&str>,
) -> Result<Option<bool>> {
    let expected = display_row_names(expected_terms);
    let actual = display_row_names(actual_terms);
    if expected | actual == 0 {
        return Ok(None);
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
        return Ok(Some(false));
    }
    let Some(reference) = e.reference.as_deref().or(a.reference.as_deref()) else {
        return Ok(Some(false));
    };
    let length = pyfmt::len(reference) as i128;
    let after_one_base = insertion_after_one_base(reference, &e.allele);
    if positions_named(e)? != Some(length)
        || (length < 2 && !after_one_base)
        || !spelt_in(&e.allele, "ACGTN")
        || e.allele.len() as i128 <= length
    {
        return Ok(Some(false));
    }
    let after_the_span = after_one_base
        || insertion_after_the_span(&differing_regions(reference, &e.allele), length as i64);
    if !after_the_span && expected & !actual == 0 {
        return Ok(Some(false));
    }
    Ok(Some(displayed_span_placement(
        positions, reference, &e.allele, expected, actual,
    )?))
}

/// Whether the reference entry `e` and the candidate entry `a` of one key are a row on
/// which Ensembl VEP writes one or more of the five intron-window splice terms or the intron
/// term from positions past the end of the record's span and the candidate does not: a
/// Transcript entry with the same cDNA and CDS positions on both sides whose two sets of
/// terms differ by terms of the six the reference alone holds and, at most, the donor
/// region and `splice_region_variant` the candidate alone holds, beside the terms of a
/// stop codon class where its rule takes the rest of the difference.
///
/// An allele that opens with its whole reference allele of two or more bases changes no base
/// of its span: every index that differs lies past the reference allele's end
/// (VariationFeatureOverlapAllele.pm:428-429), and the reference reads the further bases as a
/// change of the bases after the span (BaseTranscriptVariationAllele.pm:117-118). The
/// candidate reads such an allele as the insertion it writes, between the span's last base
/// and the base after it, and writes that insertion's terms (`insertion_after_the_span_row`):
/// on its row the candidate alone may hold the tract term and the intron term as well, and
/// the reference need hold no term alone.
///
/// An allele that opens with its reference allele of one base changes no base of its span
/// either. The reference does not compare the two alleles (:418) and reads the one base as
/// changed (:452-454), so its terms are those of that base, and the candidate's those of the
/// insertion after it (`insertion_after_one_base`): on its row the reference alone may hold
/// `splice_region_variant` as well. Where the reference's row holds none of the seven and is
/// no row beside the transcript, the candidate's may alone hold the donor or the acceptor
/// term: the base is an intron of one base, and both engines write that term for the same
/// insertion written as a record of one alternate allele.
///
/// With the record's alleles the two rows are those of a longer alternate allele that is no
/// pure insertion at a place of the intron the cDNA positions allow
/// (`record_span_placement`); without them the span lies in the transcript's sequence under
/// a longer allele (`longer_allele_inside_one_exon`), the candidate holds none of the six,
/// and `splice_region_variant` alone comes, where a term that withholds it went; no row is
/// then read as that of an allele that opens with its whole reference allele, which the
/// record alone shows.
///
/// Two entries printed under `--terms display` are read by the names the seven terms print
/// there, with the record's alleles alone (`display_vocabulary_pair`).
pub fn splice_term_past_the_record_span_pair(e: &Entry, a: &Entry) -> Result<bool> {
    splice_term_past_the_record_span_pair_beside(e, a, &|beside: &Entry, a: &Entry| {
        for pair in STOP_CODON_ROW_PAIRS {
            if pair(beside, a)? {
                return Ok(true);
            }
        }
        Ok(false)
    })
}

/// `splice_term_past_the_record_span_pair` with the rules that may take the rest of a
/// difference given: `takes_the_rest` is asked with the reference's entry rewritten to the
/// candidate's splice and intron terms beside its own other terms, and the candidate's
/// entry.
pub fn splice_term_past_the_record_span_pair_beside(
    e: &Entry,
    a: &Entry,
    takes_the_rest: &dyn Fn(&Entry, &Entry) -> Result<bool>,
) -> Result<bool> {
    if !transcript_entry(e) {
        return Ok(false);
    }
    let Some(positions) = transcript_positions(e) else {
        return Ok(false);
    };
    if transcript_positions(a).as_ref() != Some(&positions) {
        return Ok(false);
    }
    let expected_terms = term_set(&e.consequence);
    let mut actual_terms = term_set(&a.consequence);
    let cdna_and_cds = (
        positions.0.as_str(),
        positions.1.as_str(),
        positions.2.as_str(),
        positions.3.as_str(),
    );
    if let Some(taken) =
        display_vocabulary_pair(e, a, cdna_and_cds, &expected_terms, &actual_terms)?
    {
        return Ok(taken);
    }
    let read = |term: &&str| INTRON_ROW_TERMS.contains(term);
    let unread = |term: &&str| !INTRON_ROW_TERMS.contains(term);
    let expected_others: BTreeSet<&str> = expected_terms.iter().copied().filter(unread).collect();
    let actual_others: BTreeSet<&str> = actual_terms.iter().copied().filter(unread).collect();
    let rewritten;
    let mut a = a;
    if expected_others != actual_others {
        // The rest of the difference is a stop codon class's: its rule reads the two rows
        // with the candidate's splice and intron terms on both, and this one the two rows
        // with the reference's other terms on both.
        let joined = |terms: &BTreeSet<&str>| terms.iter().copied().collect::<Vec<_>>().join(",");
        let both: BTreeSet<&str> = actual_terms
            .iter()
            .copied()
            .filter(read)
            .chain(expected_others.iter().copied())
            .collect();
        let beside = Entry {
            consequence: joined(&both),
            ..e.clone()
        };
        if !takes_the_rest(&beside, a)? {
            return Ok(false);
        }
        rewritten = Entry {
            consequence: joined(&both),
            ..a.clone()
        };
        actual_terms = both;
        a = &rewritten;
    }
    let gone: BTreeSet<&str> = expected_terms.difference(&actual_terms).copied().collect();
    let come: BTreeSet<&str> = actual_terms.difference(&expected_terms).copied().collect();
    let reference = e.reference.as_deref().or(a.reference.as_deref());
    let after_one_base =
        reference.is_some_and(|reference| insertion_after_one_base(reference, &e.allele));
    let after_the_span = after_one_base
        || reference.is_some_and(|reference| {
            insertion_after_the_span(
                &differing_regions(reference, &e.allele),
                pyfmt::len(reference) as i64,
            )
        });
    let may_come: &[&str] = if after_the_span {
        &[
            SPLICE_REGION_TERM,
            SPLICE_DONOR_REGION_TERM,
            SPLICE_POLYPYRIMIDINE_TRACT_TERM,
            INTRON_TERM,
        ]
    } else {
        &[SPLICE_REGION_TERM, SPLICE_DONOR_REGION_TERM]
    };
    let differ = !gone.is_empty() || (after_the_span && !come.is_empty());
    let may_go = |term: &&str| {
        PAST_THE_SPAN_TERMS.contains(term) || (after_one_base && *term == SPLICE_REGION_TERM)
    };
    // An intron of one base that is the kept base: the window of its first two bases holds it
    // and the exon base after it, so the insertion draws a site term
    // (BaseTranscriptVariationAllele.pm:176-178) where the base, read as changed inside a
    // frameshift intron, draws none of the seven (:171-174). No other place gives such rows.
    let beside_the_transcript = ["upstream_gene_variant", "downstream_gene_variant"];
    let site = after_one_base
        && !expected_terms.iter().any(read)
        && !expected_others
            .iter()
            .any(|term| beside_the_transcript.contains(term));
    let may_come = |term: &&str| {
        may_come.contains(term)
            || (site && ["splice_donor_variant", "splice_acceptor_variant"].contains(term))
    };
    if !differ || !gone.iter().all(may_go) || !come.iter().all(may_come) {
        return Ok(false);
    }
    let Some(reference) = reference else {
        if !(longer_allele_inside_one_exon(e)? && longer_allele_inside_one_exon(a)?) {
            return Ok(false);
        }
        if actual_terms
            .iter()
            .any(|term| PAST_THE_SPAN_TERMS.contains(term))
        {
            return Ok(false);
        }
        return Ok(come.is_empty()
            || (is(&come, &[SPLICE_REGION_TERM])
                && gone
                    .iter()
                    .any(|term| SPLICE_REGION_WITHHOLDING_TERMS.contains(term))));
    };
    let length = pyfmt::len(reference) as i128;
    if positions_named(e)? != Some(length) || (length < 2 && !after_one_base) {
        return Ok(false);
    }
    if !spelt_in(&e.allele, "ACGTN") || e.allele.len() as i128 <= length {
        return Ok(false);
    }
    Ok(record_span_placement(e, a, reference)?.is_some())
}

/// The POS column as vep-rs's parsers read it: an unpadded decimal integer up to
/// `VEP_RS_COORDINATE_CEILING` exactly, or another spelling of an exact non-negative integer
/// below the ceiling that Rust's `f64` parser accepts, blanks around it ignored (`3e7`,
/// `30425650.0`, `30425650.`, `+30425650`, ` 30425650 `), read as the double it names
/// (` 72057594037927937 ` is 72057594037927936); None where vep-rs fails the line (`.`, `abc`,
/// `-5`, `0x1D04A1A`, `30_425_650`, fullwidth digits, `inf`, `1e400`, `30425650.5`).
pub fn vep_rs_pos(text: &str) -> Option<i64> {
    let stripped = text.trim_matches([' ', '\t', '\n', '\r', '\x0c']);
    if !VEP_RS_POS_WHOLE.is_match(stripped) {
        return None;
    }
    let digits = text.strip_prefix('+').unwrap_or(text);
    if _ten_digits(digits) {
        let written = plus(digits, 0);
        let value: Option<i64> = written.parse().ok();
        return value.filter(|value| *value <= VEP_RS_COORDINATE_CEILING);
    }
    let number: f64 = stripped.parse().ok()?;
    (number.is_finite() && number.fract() == 0.0 && number < 2f64.powi(62)).then_some(number as i64)
}

/// Give each of `rows` that carries no reference bases (`PairRow.reference`) the
/// bases its allele stands in for, read by position from the records of the input VCF
/// `path` (plain or gzipped): a record is asked about a row whose Location starts on one
/// of the record's own positions on its chromosome, a `chr` prefix set aside
/// (`replaced_reference`). Where two or more records read a row differently, the record
/// whose identifier is the row's own first column (`Uploaded_variation`) decides; where no
/// record carries that identifier, or more than one does and they read the row differently,
/// the row keeps none. A row that no record reads by position is read from the record on
/// its chromosome whose identifier is the row's own first column, as that record's
/// insertion or deletion moved by the 3' rule, whose Location lies off the record's
/// alleles (`moved_reference`); it keeps none where two such records read it differently.
/// The file is not opened when no row asks.
pub fn read_record_references(path: &Path, rows: &mut [Entry]) -> Result<()> {
    let mut rows: Vec<&mut Entry> = rows.iter_mut().collect();
    read_record_references_of(path, &mut rows)
}

/// The reference's own normalisation of a breakend's chromosome name at the slice test
/// (AnnotationType/Transcript.pm:134): the `chr` prefix removed and nothing else, so
/// `M` stays `M` and is not the cache's `MT`.
pub fn _strip_chr(chrom: &str) -> String {
    if chrom.len() >= 3 && chrom.is_char_boundary(3) && chrom[..3].eq_ignore_ascii_case("chr") {
        chrom[3..].to_string()
    } else {
        chrom.to_string()
    }
}

/// The fields Perl's list-context `split` keeps: the trailing empty fields dropped and the
/// middle ones kept, so `G,,T,,` is `G`, the empty field, `T`.
pub fn perl_fields(fields: &[String]) -> Vec<String> {
    let kept = fields
        .iter()
        .rposition(|field| !field.is_empty())
        .map_or(0, |last| last + 1);
    fields[..kept].to_vec()
}

/// Perl's list-context `split` of `text` on the separator `sep`, as the reference
/// splits the ALT column (BaseVCF4.pm:489): `<DEL>,` is the one field `<DEL>` and the empty
/// string has no field.
pub fn perl_split(text: &str, sep: &str) -> Vec<String> {
    if sep.is_empty() {
        crate::uncaught!("empty separator");
    }
    let fields: Vec<String> = text.split(sep).map(str::to_string).collect();
    perl_fields(&fields)
}

/// The Location both engines print for a span: `chrom:start`, or `chrom:lo-hi` with an
/// insertion's `end = start - 1` written first (`21:30425650-30425651`).
pub fn sequence_location(chrom: &str, start: i64, end: i64) -> String {
    let (lo, hi) = (start.min(end), start.max(end));
    if lo == hi {
        format!("{chrom}:{lo}")
    } else {
        format!("{chrom}:{lo}-{hi}")
    }
}

/// The Location with each coordinate written as its number: Ensembl VEP prints an
/// INFO/END as given, so `21:30425651-030425949` and `21:30425651-+30425949` name the
/// same span as `21:30425651-30425949`. A coordinate that is not an integer (`30425949.5`)
/// stays as written.
pub fn canonical_location(location: &str) -> String {
    let Some((chrom, coords)) = location.split_once(':') else {
        return location.to_string();
    };
    let coords: Vec<String> = coords.split('-').map(_canonical_coordinate).collect();
    format!("{chrom}:{}", coords.join("-"))
}

pub fn _canonical_coordinate(text: &str) -> String {
    let digits = text.strip_prefix('+').unwrap_or(text);
    if _ten_digits(digits) {
        plus(digits, 0)
    } else {
        text.to_string()
    }
}

/// The first column of a default or tab row: the record's identifier, or its position and
/// alleles where the record carries none.
pub const UPLOADED_VARIATION: &str = "Uploaded_variation";

/// `read_record_references` over entries that stand in more than one list. An input that
/// cannot be opened or read is the comparator's error `cannot read <path>: <why>`.
pub fn read_record_references_of(path: &Path, rows: &mut [&mut Entry]) -> Result<()> {
    let mut asking: HashMap<String, BTreeMap<i64, Vec<usize>>> = HashMap::new();
    // The asking rows of each chromosome by their first column, which names a row's record
    // where the record carries an identifier.
    let mut asking_by_name: HashMap<String, HashMap<Vec<u8>, Vec<usize>>> = HashMap::new();
    for (index, row) in rows.iter().enumerate() {
        let Some(location) = row.location.as_deref() else {
            continue;
        };
        if row.reference.is_some() || !location.contains(':') {
            continue;
        }
        let canonical = canonical_location(location);
        let (chrom, coords) = canonical.split_once(':').unwrap_or((&canonical, ""));
        let start = coords.split_once('-').map_or(coords, |(start, _)| start);
        if _ten_digits(start) {
            let chrom = _strip_chr(chrom);
            if let Some(name) = row.fields.get(UPLOADED_VARIATION) {
                asking_by_name
                    .entry(chrom.clone())
                    .or_default()
                    .entry(name.as_bytes().to_vec())
                    .or_default()
                    .push(index);
            }
            asking
                .entry(chrom)
                .or_default()
                .entry(integer_of(start, "Location position")?)
                .or_default()
                .push(index);
        }
    }
    if asking.is_empty() {
        return Ok(());
    }
    let mut read: Vec<BTreeSet<String>> = vec![BTreeSet::new(); rows.len()];
    // The readings of the records whose identifier is the row's own first column.
    let mut named: Vec<BTreeSet<String>> = vec![BTreeSet::new(); rows.len()];
    // The readings of a row as an insertion or deletion of its named record, moved.
    let mut moved: Vec<BTreeSet<String>> = vec![BTreeSet::new(); rows.len()];
    lines::for_each_line(path, |line| {
        if line.first() == Some(&b'#') {
            return Ok(());
        }
        let mut fields = line.splitn(6, |byte| *byte == b'\t');
        let (Some(chrom), Some(position), Some(identifier), Some(reference), Some(alternates)) = (
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
        ) else {
            return Ok(());
        };
        let chrom = _strip_chr(&pyfmt::decode_utf8_replace(chrom));
        let Some(starts) = asking.get(&chrom) else {
            return Ok(());
        };
        let Some(position) = vep_rs_pos(&pyfmt::decode_utf8_replace(position)) else {
            return Ok(());
        };
        let reference = pyfmt::decode_utf8_replace(reference);
        let alternates = perl_split(&pyfmt::decode_utf8_replace(alternates), ",");
        let alternates: Vec<&str> = alternates.iter().map(String::as_str).collect();
        // An insertion's Location starts on the base before the inserted bases, the
        // record's last where the shared bases are dropped.
        let last = position + pyfmt::len(&reference).max(1) as i64 - 1;
        for group in starts.range(position..=last).map(|(_, group)| group) {
            for &index in group {
                let row = &rows[index];
                let found = replaced_reference(
                    position,
                    &reference,
                    &alternates,
                    row.location.as_deref(),
                    &row.allele,
                )?;
                let Some(found) = found else {
                    continue;
                };
                // A record without an identifier is named by its alleles, which the row's
                // parser reads.
                let names_the_row = identifier != &b"."[..]
                    && row
                        .fields
                        .get(UPLOADED_VARIATION)
                        .is_some_and(|name| name.as_bytes() == identifier);
                if names_the_row {
                    named[index].insert(found.clone());
                }
                read[index].insert(found);
            }
        }
        // A row the 3' rule moved lies off its record's alleles: the record is the one
        // whose identifier the row's first column is.
        let named_here = asking_by_name
            .get(&chrom)
            .filter(|_| identifier != &b"."[..])
            .and_then(|names| names.get(identifier));
        for &index in named_here.into_iter().flatten() {
            let row = &rows[index];
            let found = moved_reference(
                position,
                &reference,
                &alternates,
                row.location.as_deref(),
                &row.allele,
            )?;
            if let Some(found) = found {
                moved[index].insert(found);
            }
        }
        Ok(())
    })?;
    for (((row, found), named), moved) in rows.iter_mut().zip(read).zip(named).zip(moved) {
        let decided = match (found.len(), named.len()) {
            (1, _) => found.into_iter().next(),
            (0, _) if moved.len() == 1 => moved.into_iter().next(),
            (0, _) => None,
            (_, 1) => named.into_iter().next(),
            _ => None,
        };
        if decided.is_some() {
            row.reference = decided;
        }
    }
    Ok(())
}

pub const _NON_REF_BLOCK: &str = "<NON_REF>";

pub fn _is_non_ref_block(alt: &str) -> bool {
    pyfmt::upper(alt) == _NON_REF_BLOCK
}

/// Whether an ALT names a variant: not the `<NON_REF>` block, `*` or `.`.
pub fn _names_an_allele(alt: &str) -> bool {
    !_is_non_ref_block(alt) && !matches!(alt, "*" | ".")
}

/// Whether two letters are one letter in upper case.
fn same_letter(a: char, b: char) -> bool {
    if a.is_ascii() && b.is_ascii() {
        return a.eq_ignore_ascii_case(&b);
    }
    let upper = |letter: char| pyfmt::upper(letter.encode_utf8(&mut [0; 4]));
    upper(a) == upper(b)
}

/// A bi-allelic indel minimised as vep-rs minimises it (`trim_alleles`): the shared
/// prefix is trimmed while both alleles have bases, then the shared suffix, and an emptied
/// allele reads `-`. Returns the two alleles and the bases trimmed in front, by which the
/// start moves.
pub fn _minimised(r#ref: &str, alt: &str) -> (String, String, i64) {
    let reference: Vec<char> = r#ref.chars().collect();
    let alt: Vec<char> = alt.chars().collect();
    let mut shared = 0;
    while shared < reference.len().min(alt.len()) && same_letter(reference[shared], alt[shared]) {
        shared += 1;
    }
    let (mut reference, mut alt) = (&reference[shared..], &alt[shared..]);
    while let (Some((r, kept_reference)), Some((a, kept_alt))) =
        (reference.split_last(), alt.split_last())
    {
        if !same_letter(*r, *a) {
            break;
        }
        (reference, alt) = (kept_reference, kept_alt);
    }
    let or_dash = |bases: &[char]| {
        if bases.is_empty() {
            "-".to_string()
        } else {
            bases.iter().collect()
        }
    };
    (or_dash(reference), or_dash(alt), shared as i64)
}

/// The (allele, start, end) vep-rs writes for each real allele of a gVCF variant site,
/// the sequence path of its parser: `<NON_REF>`, `*` and `.` name no allele; when any ALT's
/// length differs from REF's, a `*` of one base counted (Parser/VCF.pm:296-298), the first
/// base REF shares with every ALT is trimmed and the start moves one base, an emptied allele
/// reading `-`, an ALT carrying a `*` (`<*>`) kept whole and, beside other ALTs, left out of
/// the comparison (:306-316 with several ALTs, :328 with one; a `.` counts and never shares
/// the base); a single remaining allele whose length still differs from REF's is minimised,
/// its shared prefix then its shared suffix trimmed and the start moved by the prefix; an
/// insertion's end is its start minus one; the allele is written in upper case (`a g,<non_ref>`
/// is `G`, as Parser.pm:554 writes it). The minimisation needs the one ALT beside the block to
/// name the allele: a `*` or `.` beside the indel keeps the anchor-trimmed alleles, as Ensembl
/// VEP's `minimise_alleles` leaves a string of three parts as it is (Parser.pm:978-984). A
/// trailing empty field of the ALT list is no field (`perl_fields`).
///
/// A start or an end past `i64` is written as the nearest number `i64` holds.
pub fn gvcf_site_rows<S: AsRef<str>>(pos: i64, r#ref: &str, alts: &[S]) -> Vec<(String, i64, i64)> {
    let kept = alts
        .iter()
        .rposition(|alt| !alt.as_ref().is_empty())
        .map_or(0, |last| last + 1);
    let named: Vec<&str> = alts[..kept]
        .iter()
        .map(AsRef::as_ref)
        .filter(|alt| !_is_non_ref_block(alt))
        .collect();
    let real: Vec<&str> = named
        .iter()
        .copied()
        .filter(|alt| _names_an_allele(alt))
        .collect();
    if real.is_empty() {
        return Vec::new();
    }
    let reference_length = pyfmt::len(r#ref);
    let is_indel = named.iter().any(|alt| pyfmt::len(alt) != reference_length);
    // With one ALT its own first base decides (:328): a lone `<*>` filtered out would leave
    // nothing to compare, and `all` of nothing anchors. With several, the alleles carrying a
    // `*` are left out of the comparison (:306).
    let first_base = r#ref.chars().next();
    let anchored = first_base.is_some_and(|base| {
        named
            .iter()
            .filter(|alt| named.len() == 1 || !alt.contains('*'))
            .all(|alt| alt.chars().next().is_some_and(|c| same_letter(c, base)))
    });
    let prefix = i64::from(is_indel && anchored);
    let without_prefix = |text: &str| dash(pyfmt::slice(text, Some(prefix), None));
    let trimmed_ref = without_prefix(r#ref);
    let within = |number: i128| number.clamp(i128::from(i64::MIN), i128::from(i64::MAX)) as i64;
    let mut rows = Vec::with_capacity(real.len());
    for alt in real {
        // An allele carrying a `*` (`<*>`) is kept as written (:315).
        let mut allele = if alt.contains('*') {
            alt.to_string()
        } else {
            without_prefix(alt)
        };
        let mut start = i128::from(pos) + i128::from(prefix);
        let mut ref_part = trimmed_ref.clone();
        if named.len() == 1
            && ref_part != "-"
            && allele != "-"
            && pyfmt::len(&ref_part) != pyfmt::len(&allele)
        {
            let shared;
            (ref_part, allele, shared) = _minimised(&ref_part, &allele);
            start += i128::from(shared);
        }
        let end = if ref_part == "-" {
            start - 1
        } else {
            start + pyfmt::len(&ref_part) as i128 - 1
        };
        rows.push((pyfmt::upper(&allele), within(start), within(end)));
    }
    rows
}

pub const _START_TERMS: &[&str] = &["start_lost", "start_retained_variant"];

/// `set(csq.split(","))`: a plain split, so an empty text holds the one empty term.
fn split_terms(csq: &str) -> BTreeSet<&str> {
    csq.split(',').collect()
}

/// Whether a tuple's allele column names a structural allele.
///
/// VEP writes a sequence allele as its bases or `-` and a symbolic structural
/// allele as its Sequence Ontology class (`deletion`, `duplication`,
/// `tandem_repeat` and so on), so anything outside `-` and the nucleotide
/// alphabet is structural.
pub fn is_structural_allele(allele: &str) -> bool {
    allele != "-" && !pyfmt::upper(allele).chars().all(|c| "ACGTN".contains(c))
}

/// Whether the reference writes both `start_lost` and `start_retained_variant` and the
/// candidate the same consequence set without the member the allele kind does not support:
/// `start_retained_variant` on a structural allele (Utils/VariationEffect.pm:886-896, :961,
/// :1037), `start_lost` on a sequence one (:879-882, :906-944, :1028-1075).
pub fn _is_start_cooccurrence_swap(perl_csq: &str, rust_csq: &str, allele: &str) -> bool {
    let (pset, rset) = (split_terms(perl_csq), split_terms(rust_csq));
    if !_START_TERMS.iter().all(|term| pset.contains(term)) {
        return false;
    }
    if rset.difference(&pset).next().is_some() {
        return false;
    }
    let dropped = if is_structural_allele(allele) {
        "start_retained_variant"
    } else {
        "start_lost"
    };
    pset.difference(&rset).eq([&dropped])
}

/// The one shape `_intron_effects` (BaseTranscriptVariationAllele.pm:99-224) produces by its
/// assignment at line 215 of a value `_intron_overlap` returns as a defined 0
/// (Utils/VariationEffect.pm:107): the reference lacks no term of the candidate's set, and
/// the candidate adds `splice_region_variant`, with `intron_variant` at most beside it.
pub fn _is_covered_splice_region_swap(perl_csq: &str, rust_csq: &str) -> bool {
    let (pset, rset) = (split_terms(perl_csq), split_terms(rust_csq));
    if pset.difference(&rset).next().is_some() {
        return false;
    }
    let diff_r: Vec<&str> = rset.difference(&pset).copied().collect();
    diff_r.contains(&"splice_region_variant")
        && diff_r
            .iter()
            .all(|term| matches!(*term, "splice_region_variant" | "intron_variant"))
}

/// The recogniser of an `ExclusionRule`: over the two consequence sets, or over them and
/// the tuple's allele.
#[derive(Debug, Clone, Copy)]
pub enum ExclusionPredicate {
    Sets(fn(&str, &str) -> bool),
    SetsAndAllele(fn(&str, &str, &str) -> bool),
}

/// One divergence shape of the tuple comparison: its recogniser, its reporting name,
/// whether it assigns a pair to that name, whether the adjusted figure removes the pair, and
/// the class it belongs to.
#[derive(Debug, Clone, Copy)]
pub struct ExclusionRule {
    pub bucket: &'static str,
    pub predicate: ExclusionPredicate,
    pub classifies: bool,
    pub excludes: bool,
    pub taxonomy_class: &'static str,
}

impl ExclusionRule {
    /// Whether the rule accepts the reference's and the candidate's consequence sets of one
    /// tuple whose allele is `allele`.
    pub fn matches(&self, perl_csq: &str, rust_csq: &str, allele: &str) -> bool {
        match self.predicate {
            ExclusionPredicate::Sets(predicate) => predicate(perl_csq, rust_csq),
            ExclusionPredicate::SetsAndAllele(predicate) => predicate(perl_csq, rust_csq, allele),
        }
    }
}

/// The rules whose pairs the tuple comparison's adjusted figure removes, in its order.
pub const EXCLUDING_RULES: &[ExclusionRule] = &[
    ExclusionRule {
        bucket: "start_cooccurrence_swap",
        predicate: ExclusionPredicate::SetsAndAllele(_is_start_cooccurrence_swap),
        classifies: true,
        excludes: true,
        taxonomy_class: "start_cooccurrence_swap",
    },
    ExclusionRule {
        bucket: "covered_splice_region_swap",
        predicate: ExclusionPredicate::Sets(_is_covered_splice_region_swap),
        classifies: false,
        excludes: true,
        taxonomy_class: "splice_family_swap",
    },
];
