// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Span and term pairs and the stop-codon pairs.

use crate::coding_ends::terms::{count_alone_differs_at_the_first_residue, LandedRule};
use crate::compare::CONSEQUENCE_FIELDS;
use crate::consts::AMINO_ACIDS_FIELDS;
use crate::consts::BAM_EDIT_FIELDS;
use crate::consts::CDNA_END_FIELD;
use crate::consts::CDNA_POSITION_FIELDS;
use crate::consts::CDNA_START_FIELDS;
use crate::consts::FREQS_FIELDS;
use crate::consts::HGVSP_FIELDS;
use crate::consts::INDEL_NOTATION_RE_PATTERN;
use crate::consts::INTRON_OFFSET_RE;
use crate::consts::SHIFT_LENGTH_FIELDS;
use crate::consts::STRAND_FIELDS;
use crate::consts::USED_REF_FIELDS;
use crate::entry::Entry;
use crate::entry::Fields;
use crate::error::Result;
use crate::fields::as_number;
use crate::fields::field_named;
use crate::genetic_codes::{MODULE_CODES, STANDARD_CODE};
use crate::pyfmt;
use crate::rules::coding_end_pairs::{
    stands_beside_a_cache, start_codon_overlap_without_start_codon_pair,
    stop_codon_behind_last_residue_without_stop_codon_pair,
    stop_codon_completed_by_polyadenylation_pair, stop_codon_overlap_without_stop_codon_pair,
    structural_deletion_stop_lost_without_stop_codon_pair, AtACodingEnd, RowPair,
};
use crate::rules::downstream::positive_shift;
use crate::rules::edited_model::CODONS_FIELDS;
use crate::rules::edited_model::GIVEN_REF_FIELDS;
use crate::rules::frameshift_at_a_held_stop::frameshift_inside_a_held_stop_codon_pair;
use crate::rules::kept_start::start_codon_kept_without_five_prime_utr_pair;
use crate::rules::removed_start::start_codon_removed_without_five_prime_utr_pair;
use crate::rules::settled_stop::stop_codon_settled_before_unknown_bases_pair;
use crate::rules::shifted_indel::cdna_coordinate_order;
use crate::rules::{first_field, whole_number_text, whole_text};
use crate::self_parity::reverse_complement;
use crate::shared;
use crate::shared::blank;
use crate::shared::consequence_set;
use crate::shared::feature_kind;
use crate::shared::inserted_between_neighbours;
use crate::shared::insertion_row;
use crate::shared::location_span;
use crate::shared::read_record_references_of;
use crate::shared::ABSENT;
use crate::shared::CODING_SEQUENCE_TERM;
use crate::shared::STOP_LOST_TERM;
use crate::shared::STOP_RETAINED_TERM;
use std::borrow::Cow;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::path::Path;
use std::sync::LazyLock;

pub use crate::rules::following::hgvsp_frameshift_form;
pub use crate::rules::following::term_pair_following_cell;
pub use crate::rules::following::terms_impact;
pub use crate::rules::following::terms_write_field;
pub use crate::rules::following::DOWNSTREAM_FIELDS;
pub use crate::rules::following::TERM_IMPACT;
pub use crate::rules::following::TRUNCATING_TERMS;
pub use crate::rules::following::TRUNCATING_TERM_FIELDS;

pub const CLAMPED_NOTATION_RE_PATTERN: &str =
    r"([cn])\.(-?[1-9]\d*)(?:[ACGTN-]>([ACGT])|delins([ACGT]+))";
pub static CLAMPED_NOTATION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(CLAMPED_NOTATION_RE_PATTERN).expect("the pattern compiles"));

pub const INSERTION_NOTATION_RE_PATTERN: &str = r"([cn])\.(-?[1-9]\d*)_(-?[1-9]\d*)ins([ACGT]+)";
pub static INSERTION_NOTATION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(INSERTION_NOTATION_RE_PATTERN).expect("the pattern compiles")
});

pub const SHIFTED_DUPLICATION_RE_PATTERN: &str = r"([cn])\.(?:(-?[1-9]\d*)_)?(-?[1-9]\d*)dup";
pub static SHIFTED_DUPLICATION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(SHIFTED_DUPLICATION_RE_PATTERN).expect("the pattern compiles")
});

static CLAMPED_NOTATION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(CLAMPED_NOTATION_RE_PATTERN));
static INSERTION_NOTATION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(INSERTION_NOTATION_RE_PATTERN));
static SHIFTED_DUPLICATION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(SHIFTED_DUPLICATION_RE_PATTERN));
static INDEL_NOTATION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(INDEL_NOTATION_RE_PATTERN));

pub const ACCEPTOR_OFFSETS: &[&str] = &["-1", "-2"];

pub const MATURE_MIRNA_TERM: &str = "mature_miRNA_variant";

pub const NON_CODING_EXON_TERM: &str = "non_coding_transcript_exon_variant";

pub const IMPACT_FIELDS: &[&str] = &["IMPACT", "impact"];

pub const INTERGENIC_TERM: &str = "intergenic_variant";

pub const MISSENSE_TERM: &str = "missense_variant";

/// Membership only.
pub const DELETED_SPAN_REGION_TERMS: &[&str] = &[
    "3_prime_UTR_variant",
    "5_prime_UTR_variant",
    "coding_sequence_variant",
];

pub const THREE_PRIME_UTR_TERM: &str = "3_prime_UTR_variant";

/// The terms both engines write on a row whose span holds a whole intron. Membership only.
pub const INTRON_SPAN_TERMS: &[&str] = &[
    "splice_acceptor_variant",
    "splice_donor_variant",
    "intron_variant",
];

pub const PROTEIN_ALTERING_TERM: &str = "protein_altering_variant";

pub const INFRAME_DELETION_TERM: &str = "inframe_deletion";

/// The stop-lost extension with an unknown count, which TranscriptVariationAllele.pm:1916 and
/// :1926 write after the notation of a deletion on a row that passes `stop_lost`.
pub const STOP_LOST_EXTENSION: &str = "extTer?";

/// The row rules beside whose pair a field of a plugin a term gates that one engine alone
/// writes is `plugin_value_follows_stop_codon_term_pair`'s: three rules of a stop codon that
/// stands, the rule of a stop codon settled before unknown bases, the rule of a frameshift
/// inside a held stop codon, and the five classes of a coding end, on whose rows the engines
/// read different coding terms.
pub const STOP_CODON_PLUGIN_VALUE_PAIRS: &[RowPair] = &[
    stop_retained_without_stop_codon_test_pair_at,
    stop_codon_overlap_insertion_range_pair_at,
    stop_retained_three_base_stop_codon_test_pair_at,
    stop_codon_settled_before_unknown_bases_pair_at,
    frameshift_inside_a_held_stop_codon_pair_at,
    structural_deletion_stop_lost_without_stop_codon_pair,
    stop_codon_completed_by_polyadenylation_pair,
    stop_codon_overlap_without_stop_codon_pair,
    stop_codon_behind_last_residue_without_stop_codon_pair,
    start_codon_overlap_without_start_codon_pair,
];

/// The built-in row classes that accept a transcript row's Consequence, each with the predicate
/// deciding whether a reference entry and the candidate's of one key have its shape, on what
/// the comparison read of the row's transcript (`AtACodingEnd`); a row is the first class's
/// whose predicate holds. The five classes of a coding end stand in the order they are asked
/// (`CodingEndClass::IN_ORDER`), before the classes of a stop codon that stands.
pub const TERM_PAIR_CLASSES: &[(&str, RowPair)] = &[
    (
        "refseq_edited_model_deleted_span_term",
        deleted_span_term_pair_at,
    ),
    ("refseq_edited_model_tail_term", tail_term_pair_at),
    (
        "mature_mirna_unshifted_span",
        mature_mirna_unshifted_span_pair_at,
    ),
    (
        "stop_codon_unshifted_span",
        stop_codon_unshifted_span_pair_at,
    ),
    (
        "structural_deletion_stop_lost_without_stop_codon",
        structural_deletion_stop_lost_without_stop_codon_pair,
    ),
    (
        "stop_codon_completed_by_polyadenylation",
        stop_codon_completed_by_polyadenylation_pair,
    ),
    (
        "stop_codon_overlap_without_stop_codon",
        stop_codon_overlap_without_stop_codon_pair,
    ),
    (
        "stop_codon_behind_last_residue_without_stop_codon",
        stop_codon_behind_last_residue_without_stop_codon_pair,
    ),
    (
        "start_codon_overlap_without_start_codon",
        start_codon_overlap_without_start_codon_pair,
    ),
    (
        "stop_retained_without_stop_codon_test",
        stop_retained_without_stop_codon_test_pair_at,
    ),
    (
        "stop_codon_overlap_genomic_span",
        stop_codon_overlap_genomic_span_pair_at,
    ),
    (
        "stop_retained_three_base_stop_codon_test",
        stop_retained_three_base_stop_codon_test_pair_at,
    ),
    (
        "stop_codon_removed_without_stop_codon_test",
        stop_codon_removed_without_stop_codon_test_pair_at,
    ),
    (
        "stop_codon_overlap_insertion_range",
        stop_codon_overlap_insertion_range_pair_at,
    ),
    (
        "refseq_edited_model_equal_span_insertion_term",
        equal_span_insertion_term_pair_at,
    ),
    (
        "cnv_tr_symbolic_loss_named_a_gain",
        symbolic_loss_named_a_gain_pair_at,
    ),
    (
        "splice_term_past_the_record_span",
        splice_term_past_the_record_span_pair_at,
    ),
    (
        "stop_codon_settled_before_unknown_bases",
        stop_codon_settled_before_unknown_bases_pair_at,
    ),
    (
        "start_codon_removed_without_five_prime_utr",
        start_codon_removed_without_five_prime_utr_pair_at,
    ),
    (
        "frameshift_inside_a_held_stop_codon",
        frameshift_inside_a_held_stop_codon_pair_at,
    ),
    (
        "start_codon_kept_without_five_prime_utr",
        start_codon_kept_without_five_prime_utr_pair_at,
    ),
];

/// `deleted_span_term_pair`, which reads nothing of the row's transcript.
fn deleted_span_term_pair_at(e: &Entry, a: &Entry, _: &AtACodingEnd<'_>) -> Result<bool> {
    deleted_span_term_pair(e, a)
}

/// `tail_term_pair`, which reads nothing of the row's transcript.
fn tail_term_pair_at(e: &Entry, a: &Entry, _: &AtACodingEnd<'_>) -> Result<bool> {
    tail_term_pair(e, a)
}

/// `equal_span_insertion_term_pair`, which reads nothing of the row's transcript.
fn equal_span_insertion_term_pair_at(e: &Entry, a: &Entry, _: &AtACodingEnd<'_>) -> Result<bool> {
    equal_span_insertion_term_pair(e, a)
}

/// `symbolic_loss_named_a_gain_pair`, which reads nothing of the row's transcript.
fn symbolic_loss_named_a_gain_pair_at(e: &Entry, a: &Entry, _: &AtACodingEnd<'_>) -> Result<bool> {
    symbolic_loss_named_a_gain_pair(e, a)
}

/// `stop_codon_settled_before_unknown_bases_pair`, which reads nothing of the row's
/// transcript.
fn stop_codon_settled_before_unknown_bases_pair_at(
    e: &Entry,
    a: &Entry,
    _: &AtACodingEnd<'_>,
) -> Result<bool> {
    stop_codon_settled_before_unknown_bases_pair(e, a)
}

/// `start_codon_removed_without_five_prime_utr_pair`, which reads nothing of the row's
/// transcript.
fn start_codon_removed_without_five_prime_utr_pair_at(
    e: &Entry,
    a: &Entry,
    _: &AtACodingEnd<'_>,
) -> Result<bool> {
    start_codon_removed_without_five_prime_utr_pair(e, a)
}

/// `frameshift_inside_a_held_stop_codon_pair`, which reads nothing of the row's transcript.
fn frameshift_inside_a_held_stop_codon_pair_at(
    e: &Entry,
    a: &Entry,
    _: &AtACodingEnd<'_>,
) -> Result<bool> {
    frameshift_inside_a_held_stop_codon_pair(e, a)
}

/// `start_codon_kept_without_five_prime_utr_pair`, which reads nothing of the row's
/// transcript.
fn start_codon_kept_without_five_prime_utr_pair_at(
    e: &Entry,
    a: &Entry,
    _: &AtACodingEnd<'_>,
) -> Result<bool> {
    start_codon_kept_without_five_prime_utr_pair(e, a)
}

/// `mature_mirna_unshifted_span_pair`, which reads nothing of the row's transcript.
fn mature_mirna_unshifted_span_pair_at(e: &Entry, a: &Entry, _: &AtACodingEnd<'_>) -> Result<bool> {
    mature_mirna_unshifted_span_pair(e, a)
}

/// `stop_codon_unshifted_span_pair` in a comparison that may have been given the
/// candidate's cache (`stands_beside_a_cache`): beside a cache a stop codon stands.
pub fn stop_codon_unshifted_span_pair_at(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    Ok(stop_codon_unshifted_span_pair(e, a)?
        && stands_beside_a_cache(LandedRule::StopCodonUnshiftedSpan, e, a, at))
}

/// `stop_retained_without_stop_codon_test_pair` in a comparison that may have been given
/// the candidate's cache (`stands_beside_a_cache`): beside a cache a stop codon stands.
pub fn stop_retained_without_stop_codon_test_pair_at(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    Ok(stop_retained_without_stop_codon_test_pair(e, a)?
        && stands_beside_a_cache(LandedRule::StopRetainedWithoutStopCodonTest, e, a, at))
}

/// `stop_codon_overlap_genomic_span_pair` in a comparison that may have been given the
/// candidate's cache (`stands_beside_a_cache`): beside a cache a stop codon stands.
pub fn stop_codon_overlap_genomic_span_pair_at(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    Ok(stop_codon_overlap_genomic_span_pair(e, a)?
        && stands_beside_a_cache(LandedRule::StopCodonOverlapGenomicSpan, e, a, at))
}

/// `stop_retained_three_base_stop_codon_test_pair` in a comparison that may have been given
/// the candidate's cache (`stands_beside_a_cache`): beside a cache the stop codon is
/// whole, since on one that polyadenylation completes the engines read two triplets.
pub fn stop_retained_three_base_stop_codon_test_pair_at(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    Ok(stop_retained_three_base_stop_codon_test_pair(e, a)?
        && stands_beside_a_cache(LandedRule::StopRetainedThreeBaseStopCodonTest, e, a, at))
}

/// `stop_codon_removed_without_stop_codon_test_pair` in a comparison that may have been
/// given the candidate's cache (`stands_beside_a_cache`): beside a cache a stop codon
/// stands.
pub fn stop_codon_removed_without_stop_codon_test_pair_at(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    Ok(stop_codon_removed_without_stop_codon_test_pair(e, a)?
        && stands_beside_a_cache(LandedRule::StopCodonRemovedWithoutStopCodonTest, e, a, at))
}

/// `stop_codon_overlap_insertion_range_pair` in a comparison that may have been given the
/// candidate's cache (`stands_beside_a_cache`): beside a cache a stop codon stands.
pub fn stop_codon_overlap_insertion_range_pair_at(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    Ok(stop_codon_overlap_insertion_range_pair(e, a)?
        && stands_beside_a_cache(LandedRule::StopCodonOverlapInsertionRange, e, a, at))
}

/// `int(end) - int(start)` of two whole numbers written in digits alone without a leading
/// zero, at any length; a difference past `i128` is the nearest bound.
fn difference(end: &str, start: &str) -> i128 {
    if let (Ok(end), Ok(start)) = (end.parse::<i128>(), start.parse::<i128>()) {
        return end - start;
    }
    let negative = (end.len(), end) < (start.len(), start);
    let (larger, smaller) = if negative { (start, end) } else { (end, start) };
    let mut lower = smaller.bytes().rev();
    let mut borrow = 0;
    let mut digits: Vec<u8> = Vec::with_capacity(larger.len());
    for upper in larger.bytes().rev() {
        let taken = lower.next().map_or(0, |digit| digit - b'0') + borrow;
        let upper = upper - b'0';
        borrow = u8::from(upper < taken);
        digits.push(upper + 10 * borrow - taken);
    }
    let magnitude = digits
        .iter()
        .rev()
        .skip_while(|digit| **digit == 0)
        .try_fold(0i128, |number, digit| {
            number.checked_mul(10)?.checked_add(i128::from(*digit))
        });
    match (magnitude, negative) {
        (Some(magnitude), false) => magnitude,
        (Some(magnitude), true) => -magnitude,
        (None, false) => i128::MAX,
        (None, true) => i128::MIN,
    }
}

/// Whether the entry is a Transcript's: the last segment of its key, read up to a `#`,
/// is a Feature_type of the transcript list (`feature_kind`).
pub(crate) fn transcript_row(e: &Entry) -> bool {
    let last = e.key.rsplit('|').next().unwrap_or("");
    let feature_type = last.split('#').next().unwrap_or("");
    feature_kind(feature_type) == "transcript"
}

/// Whether the allele is bases, or `-` for a deletion's.
fn sequence_allele(allele: &str) -> bool {
    !allele.is_empty() && allele.bytes().all(|letter| b"ACGTN-".contains(&letter))
}

/// Whether the entry carries `BAM_EDIT` `OK`.
pub(crate) fn edited_model_row(e: &Entry) -> bool {
    BAM_EDIT_FIELDS
        .iter()
        .any(|name| e.fields.get(*name).is_some_and(|value| value == "OK"))
}

/// Whether the entry carries the cDNA position field `name`.
fn positioned(e: &Entry, name: &str) -> bool {
    !blank(e.fields.get(name).map_or(ABSENT, String::as_str))
}

/// Whether a candidate notation's coordinates agree with the row's cDNA position under a
/// run that prints `SHIFT_LENGTH`, where that position is the shifted span's. In `n.`
/// numbering, which is cDNA numbering, a plain coordinate on both sides must be equal; in
/// `c.` numbering two plain coordinates of one region (the 5' UTR, the coding sequence or
/// the 3' UTR) span as many bases as the row's two plain sides. A side the row leaves open
/// (`?`), a coordinate with an intron offset and a `c.` span across a region boundary say
/// nothing, as does a run without `SHIFT_LENGTH`, which positions the row at the written
/// span. The text formats write both sides in one `cDNA_position` (`23-43`, `?-31`,
/// `23` for one base, `23-43/1500` under --total_length); the JSON format writes
/// `cdna_start` and `cdna_end` apart and leaves an open side out (JSON.pm:340-341).
pub fn span_matches_cdna_position(
    e: &Entry,
    numbering: &str,
    start: &str,
    end: &str,
) -> Result<bool> {
    if SHIFT_LENGTH_FIELDS
        .iter()
        .all(|name| !e.fields.contains_key(*name))
    {
        return Ok(true);
    }
    let (text_position_field, json_start_field) = (CDNA_START_FIELDS[0], CDNA_START_FIELDS[1]);
    let (row_start, row_end) = match e.fields.get(text_position_field) {
        Some(position) => {
            // Under --total_length the position carries the transcript's length after a slash.
            let written = pyfmt::partition(position, "/").0;
            let (row_start, separator, row_end) = pyfmt::partition(written, "-");
            let row_end = if separator.is_empty() {
                row_start
            } else {
                row_end
            };
            (row_start, row_end)
        }
        None => (
            e.fields.get(json_start_field).map_or("", String::as_str),
            e.fields.get(CDNA_END_FIELD).map_or("", String::as_str),
        ),
    };
    if numbering == "n" {
        for (row_position, coordinate) in [(row_start, start), (row_end, end)] {
            if pyfmt::is_digit(row_position)
                && pyfmt::is_digit(coordinate)
                && row_position != coordinate
            {
                return Ok(false);
            }
        }
        return Ok(true);
    }
    if !(pyfmt::is_digit(row_start) && pyfmt::is_digit(row_end)) {
        return Ok(true);
    }
    let (start_region, start_position, start_offset) = cdna_coordinate_order(start)?;
    let (end_region, end_position, end_offset) = cdna_coordinate_order(end)?;
    if start_offset != 0 || end_offset != 0 || start_region != end_region {
        return Ok(true);
    }
    let row_end = whole_number_text(row_end);
    let row_start = whole_number_text(row_start);
    Ok(i128::from(end_position) - i128::from(start_position) == difference(&row_end, &row_start))
}

/// The intron offsets of the coordinates that place an insertion's or deletion's shifted
/// span, where an HGVSc body places it at an acceptor site's two bases or in the exon,
/// touching no intron base away from the splice site: every coordinate is exonic or carries
/// the intron offset -1 or -2, a duplication's first coordinate excepted, since the
/// duplicated bases run back into the intron while the insertion sits after the last. An
/// empty tuple is a notation whose coordinates are all exonic; `None` a coordinate at
/// offset -3 or beyond, or on the donor side, or any other notation.
pub fn acceptor_site_offsets(body: &str) -> Option<Vec<String>> {
    let m = INDEL_NOTATION_WHOLE.captures(body)?;
    let (first, kind) = (&m[1], &m[3]);
    let coordinates = match m.get(2).map(|last| last.as_str()) {
        None => vec![first],
        Some(last) if kind == "dup" => vec![last],
        Some(last) => vec![first, last],
    };
    let mut offsets = Vec::new();
    for coordinate in coordinates {
        if let Some(offset) = INTRON_OFFSET_RE.find(coordinate).ok().flatten() {
            if !ACCEPTOR_OFFSETS.contains(&offset.as_str()) {
                return None;
            }
            offsets.push(offset.as_str().to_string());
        }
    }
    Some(offsets)
}

/// The cDNA coordinate one base from `position`, `None` where it is past `i64`.
fn stepped(position: i64, step: i64) -> Option<i64> {
    let moved = position.checked_add(step)?;
    if moved == 0 {
        moved.checked_add(step)
    } else {
        Some(moved)
    }
}

/// The cDNA coordinate one base (`step` 1 or -1) from `position`: the `c.` numbering
/// has no 0, so -1 is followed by 1. A coordinate past `i64` is the nearest one it holds.
pub fn cdna_step(position: i64, step: i64) -> i64 {
    stepped(position, step).unwrap_or(if step < 0 { i64::MIN } else { i64::MAX })
}

/// The numbering, the position P, the inserted bases and the count of bases stripped from
/// them of the notation Ensembl VEP writes for an insertion clamped to one base of the
/// transcript slice (`hgvs_shifted_insertion_clamped_end`): `P X>Y` and `P delins`
/// carry the bases whole, `P-1_P ins` all but the last; `None` for any other notation,
/// and for a position past `i64`. `X` is `-` where a model that carries an edit has the
/// row's reference allele, an insertion's `-`, written over the base the notation was typed
/// with (TranscriptVariationAllele.pm:1508-1514).
pub fn clamped_insertion_position(body: &str) -> Option<(String, i64, String, i64)> {
    if let Some(m) = CLAMPED_NOTATION_WHOLE.captures(body) {
        let position = whole_number_text(&m[2]);
        let bases = m.get(3).or(m.get(4))?.as_str();
        return Some((
            m[1].to_string(),
            position.parse().ok()?,
            bases.to_string(),
            0,
        ));
    }
    let m = INSERTION_NOTATION_WHOLE.captures(body)?;
    let last = whole_number_text(&m[3]);
    let first = whole_number_text(&m[2]);
    let (last, first): (i64, i64) = (last.parse().ok()?, first.parse().ok()?);
    (stepped(last, -1) == Some(first)).then(|| (m[1].to_string(), last, m[4].to_string(), 1))
}

/// The numbering, the position P after the notation of an insertion at its shifted span,
/// its inserted bases and their count: a duplication ending at P-1, whose notation carries
/// no bases (`None`) and whose span is their count (a span crossing from c.-1 to c.1 holds
/// one base fewer than its bounds say, since the `c.` numbering has no 0), or an insertion
/// between P-1 and P; `None` for any other notation, and for a position past `i64`.
pub fn shifted_insertion_end(body: &str) -> Option<(String, i64, Option<String>, i64)> {
    if let Some(m) = SHIFTED_DUPLICATION_WHOLE.captures(body) {
        let end = whole_number_text(&m[3]);
        let start = m.get(2).map(|start| whole_number_text(start.as_str()));
        let end: i64 = end.parse().ok()?;
        let start: i64 = match start {
            Some(start) => start.parse().ok()?,
            None => end,
        };
        let mut duplicated_count = i128::from(end) - i128::from(start) + 1;
        if start < 0 && 0 < end {
            duplicated_count -= 1;
        }
        return Some((
            m[1].to_string(),
            stepped(end, 1)?,
            None,
            i64::try_from(duplicated_count).ok()?,
        ));
    }
    let m = INSERTION_NOTATION_WHOLE.captures(body)?;
    let last = whole_number_text(&m[3]);
    let first = whole_number_text(&m[2]);
    let (last, first): (i64, i64) = (last.parse().ok()?, first.parse().ok()?);
    (stepped(last, -1) == Some(first)).then(|| {
        (
            m[1].to_string(),
            last,
            Some(m[4].to_string()),
            m[4].len() as i64,
        )
    })
}

/// Whether the reference entry `e` and the candidate entry `a` of one key have the
/// shape of `refseq_edited_model_deleted_span_term`: a Transcript entry of a sequence allele
/// (bases, or `-` for a deletion's), `BAM_EDIT` `OK` on both, the reference's terms
/// `intergenic_variant` alone and the candidate's exactly one of
/// `DELETED_SPAN_REGION_TERMS`, no cDNA position on either and one IMPACT on both.
pub fn deleted_span_term_pair(e: &Entry, a: &Entry) -> Result<bool> {
    if !transcript_row(e) || !sequence_allele(&e.allele) {
        return Ok(false);
    }
    if consequence_set(&e.consequence) != INTERGENIC_TERM
        || !DELETED_SPAN_REGION_TERMS.contains(&consequence_set(&a.consequence).as_str())
    {
        return Ok(false);
    }
    for entry in [e, a] {
        if !edited_model_row(entry)
            || CDNA_POSITION_FIELDS
                .iter()
                .any(|name| positioned(entry, name))
        {
            return Ok(false);
        }
    }
    Ok(first_field(e, IMPACT_FIELDS) == first_field(a, IMPACT_FIELDS))
}

/// Whether the reference entry `e` and the candidate entry `a` of one key have the
/// shape of `refseq_edited_model_tail_term`: a Transcript entry of a sequence allele,
/// `BAM_EDIT` `OK` on both, the reference's terms `intergenic_variant` alone and the
/// candidate's `3_prime_UTR_variant` alone, a cDNA position on both and one IMPACT on both.
pub fn tail_term_pair(e: &Entry, a: &Entry) -> Result<bool> {
    if !transcript_row(e) || !sequence_allele(&e.allele) {
        return Ok(false);
    }
    if consequence_set(&e.consequence) != INTERGENIC_TERM
        || consequence_set(&a.consequence) != THREE_PRIME_UTR_TERM
    {
        return Ok(false);
    }
    for entry in [e, a] {
        if !edited_model_row(entry)
            || CDNA_POSITION_FIELDS
                .iter()
                .all(|name| !positioned(entry, name))
        {
            return Ok(false);
        }
    }
    Ok(first_field(e, IMPACT_FIELDS) == first_field(a, IMPACT_FIELDS))
}

/// Whether the reference entry `e` and the candidate entry `a` of one key have the
/// shape of `refseq_edited_model_equal_span_insertion_term`: a Transcript entry of an
/// insertion between two neighbouring positions (`inserted_between_neighbours`),
/// `BAM_EDIT` `OK` on both, the reference's terms `intergenic_variant` alone and the
/// candidate's `missense_variant` alone, and one CDS span (`cds_bounds`), the same on both
/// sides, that holds as many bases as the allele. The IMPACT is no part of the shape: it
/// follows the two terms (`term_pair_following_cell`).
pub fn equal_span_insertion_term_pair(e: &Entry, a: &Entry) -> Result<bool> {
    if !transcript_row(e) || !insertion_row(e) || !inserted_between_neighbours(e) {
        return Ok(false);
    }
    if consequence_set(&e.consequence) != INTERGENIC_TERM
        || consequence_set(&a.consequence) != MISSENSE_TERM
    {
        return Ok(false);
    }
    if !edited_model_row(e) || !edited_model_row(a) {
        return Ok(false);
    }
    let span = cds_bounds(e);
    if span.is_none() || span != cds_bounds(a) {
        return Ok(false);
    }
    Ok(cds_span_length(e).is_some_and(|bases| i128::from(bases) == pyfmt::len(&e.allele) as i128))
}

/// The first and last coding base of the entry's span: the `CDS_position` of the text
/// formats (one position or a range) or the JSON `cds_start` and `cds_end`, each bound read
/// up to the `/total` suffix `--total_length` appends; `None` where a bound is missing, is
/// not a number (`?`) or is past `i64`.
pub fn cds_bounds(e: &Entry) -> Option<(i64, i64)> {
    let (start, end) = match e.fields.get("CDS_position") {
        Some(position) => {
            let (start, _, end) = pyfmt::partition(pyfmt::partition(position, "/").0, "-");
            (start, if end.is_empty() { start } else { end })
        }
        None => {
            let bound = |name: &str| {
                pyfmt::partition(e.fields.get(name).map_or(ABSENT, String::as_str), "/").0
            };
            (bound("cds_start"), bound("cds_end"))
        }
    };
    if !pyfmt::is_digit(start) || !pyfmt::is_digit(end) {
        return None;
    }
    let start = whole_number_text(start);
    let end = whole_number_text(end);
    Some((start.parse().ok()?, end.parse().ok()?))
}

/// The number of coding bases the entry's span holds (`cds_bounds`), `None` where a
/// bound is missing or is not a number.
pub fn cds_span_length(e: &Entry) -> Option<i64> {
    let (start, end) = cds_bounds(e)?;
    i64::try_from(i128::from(end) - i128::from(start) + 1).ok()
}

/// Whether the Codons `expected_value` of the reference entry `e` and `actual_value`
/// of the candidate entry `a` have the shape of `refseq_edited_model_deleted_span_codons`:
/// a Transcript entry of a sequence allele carrying `BAM_EDIT` `OK` on both whose USED_REF
/// has more bases than its CDS span, the two values' alternate codons equal, their reference
/// codons of one length agreeing at every lower-case base and differing in an upper-case one,
/// the reference's upper-case bases one run that ends its codons and spells, on the
/// transcript's strand, the leading bases of the allele the splice writes, which USED_REF
/// holds or, where the row gives one of USED_REF's length, GIVEN_REF. The Codons are written
/// (OutputFactory.pm:1743) before the protein notation of a deletion that the 3' shift moves
/// (OutputFactory.pm:1769) puts the shifted bases over the reference allele
/// (TranscriptVariationAllele.pm:1785-1786), which USED_REF then prints
/// (OutputFactory.pm:1785), so on such a row GIVEN_REF alone holds the spliced allele.
pub fn deleted_span_codons_pair(
    e: &Entry,
    a: &Entry,
    expected_value: &str,
    actual_value: &str,
) -> bool {
    if !transcript_row(e) {
        return false;
    }
    if !sequence_allele(&e.allele) || blank(expected_value) || blank(actual_value) {
        return false;
    }
    if !edited_model_row(e) || !edited_model_row(a) {
        return false;
    }
    let used_ref = first_field(e, USED_REF_FIELDS);
    let span = cds_span_length(e);
    let Some(span) = span.filter(|_| !blank(used_ref)) else {
        return false;
    };
    if pyfmt::len(used_ref) as i128 <= i128::from(span) {
        return false;
    }
    let (expected_reference_codons, _, expected_alternate_codons) =
        pyfmt::partition(expected_value, "/");
    let (actual_reference_codons, _, actual_alternate_codons) = pyfmt::partition(actual_value, "/");
    if expected_alternate_codons != actual_alternate_codons {
        return false;
    }
    if pyfmt::len(expected_reference_codons) != pyfmt::len(actual_reference_codons)
        || expected_reference_codons == actual_reference_codons
    {
        return false;
    }
    let lower = |letter: char| pyfmt::is_lower(letter.encode_utf8(&mut [0; 4]));
    let upper = |letter: char| pyfmt::is_upper(letter.encode_utf8(&mut [0; 4]));
    if expected_reference_codons
        .chars()
        .zip(actual_reference_codons.chars())
        .any(|(x, y)| (lower(x) || lower(y)) && x != y)
    {
        return false;
    }
    // The splice writes the allele from the variant's first base to the window's end, so the
    // reference's upper-case bases are one run that ends the window and spells the allele's
    // leading bases on the transcript's strand.
    let Some((upper_start, _)) = expected_reference_codons
        .char_indices()
        .find(|(_, letter)| upper(*letter))
    else {
        return false;
    };
    let reverse_strand = first_field(e, STRAND_FIELDS) == "-1";
    // A lower-case base inside the run fails the prefix test, since the allele is upper-cased.
    let leads_with_the_run = |allele: &str| {
        let allele = if reverse_strand {
            Cow::Owned(reverse_complement(allele))
        } else {
            Cow::Borrowed(allele)
        };
        pyfmt::upper(&allele).starts_with(&expected_reference_codons[upper_start..])
    };
    if leads_with_the_run(used_ref) {
        return true;
    }
    // USED_REF holds a shifted deletion's bases where the protein notation follows the
    // Codons (TranscriptVariationAllele.pm:1785-1786); GIVEN_REF keeps the allele as the
    // record gives it, which is the spliced one.
    let given_ref = first_field(e, GIVEN_REF_FIELDS);
    !blank(given_ref)
        && pyfmt::len(given_ref) == pyfmt::len(used_ref)
        && leads_with_the_run(given_ref)
}

/// Whether the reference entry `e` and the candidate entry `a` of one key have the row
/// shape of `refseq_edited_model_intron_span_coding_reference`: a Transcript entry of a
/// sequence allele carrying `BAM_EDIT` `OK` on both, the coding span (`cds_bounds`) one
/// range of two or more positions on both, the USED_REF one text of letters on both with
/// more of them than the range has positions, and each of `INTRON_SPAN_TERMS` among the
/// terms of both.
pub fn intron_span_row(e: &Entry, a: &Entry) -> bool {
    if !transcript_row(e) || !sequence_allele(&e.allele) {
        return false;
    }
    if !edited_model_row(e) || !edited_model_row(a) {
        return false;
    }
    let Some((start, end)) = cds_bounds(e).filter(|bounds| cds_bounds(a) == Some(*bounds)) else {
        return false;
    };
    let used_ref = first_field(e, USED_REF_FIELDS);
    if start >= end
        || used_ref != first_field(a, USED_REF_FIELDS)
        || !used_ref.bytes().all(|letter| letter.is_ascii_alphabetic())
    {
        return false;
    }
    if used_ref.len() as i128 <= i128::from(end) - i128::from(start) + 1 {
        return false;
    }
    [e, a].iter().all(|entry| {
        let terms = consequence_set(&entry.consequence);
        INTRON_SPAN_TERMS
            .iter()
            .all(|term| terms.split(',').any(|held| held == *term))
    })
}

/// Whether two Consequence values differ as the term pair of
/// `refseq_edited_model_intron_span_coding_reference`: the reference's terms less the
/// candidate's are `stop_lost` and `protein_altering_variant`, and the candidate's less the
/// reference's are `inframe_deletion` alone.
fn intron_span_term_pair(expected_value: &str, actual_value: &str) -> bool {
    let (expected_set, actual_set) = (
        consequence_set(expected_value),
        consequence_set(actual_value),
    );
    let expected_terms: BTreeSet<&str> = expected_set.split(',').collect();
    let actual_terms: BTreeSet<&str> = actual_set.split(',').collect();
    let reference_only: Vec<&str> = expected_terms.difference(&actual_terms).copied().collect();
    let candidate_only: Vec<&str> = actual_terms.difference(&expected_terms).copied().collect();
    // The terms one side alone holds come in their sorted order.
    reference_only == [PROTEIN_ALTERING_TERM, STOP_LOST_TERM]
        && candidate_only == [INFRAME_DELETION_TERM]
}

/// The residues of `codons`, whole codons of the letters `A`, `C`, `G` and `T` in either
/// case, under the standard genetic code, a stop as `*`; `None` for bases that are no whole
/// number of codons and for a codon with any other letter.
fn standard_translation(codons: &str) -> Option<String> {
    let code = MODULE_CODES.get(&STANDARD_CODE)?;
    let bases = codons.to_ascii_uppercase();
    if !bases.len().is_multiple_of(3) {
        return None;
    }
    bases
        .as_bytes()
        .chunks_exact(3)
        .map(|codon| code.amino_acid(codon).map(char::from))
        .collect()
}

/// Whether the Amino_acids `expected_value` of the reference entry `e` and `actual_value` of
/// the candidate entry `a` are the residues of a row whose two Codons are a
/// `deleted_span_codons_pair`: the alternate halves equal, the reference halves of one
/// length and different, each the translation of its own row's reference codons
/// (`standard_translation`), which TranscriptVariationAllele.pm:703 and :724-736 make of the
/// reference allele's codons.
fn deleted_span_residues_pair(
    e: &Entry,
    a: &Entry,
    expected_value: &str,
    actual_value: &str,
) -> bool {
    let expected_codons = first_field(e, CODONS_FIELDS);
    let actual_codons = first_field(a, CODONS_FIELDS);
    if !deleted_span_codons_pair(e, a, expected_codons, actual_codons) {
        return false;
    }
    let (expected_reference, _, expected_alternate) = pyfmt::partition(expected_value, "/");
    let (actual_reference, _, actual_alternate) = pyfmt::partition(actual_value, "/");
    if expected_alternate != actual_alternate
        || pyfmt::len(expected_reference) != pyfmt::len(actual_reference)
        || expected_reference == actual_reference
    {
        return false;
    }
    let translates_to = |codons: &str, residues: &str| {
        standard_translation(pyfmt::partition(codons, "/").0).as_deref() == Some(residues)
    };
    translates_to(expected_codons, expected_reference)
        && translates_to(actual_codons, actual_reference)
}

/// Whether the two codon readings of Ensembl VEP's `inframe_deletion`
/// (Utils/VariationEffect.pm:1173-1185) accept a Codons value, its bases compared in one
/// case and a `-` half read as no base: the alternate codons are fewer than the reference's,
/// and either the reference codons begin or end with the alternate's (:1178), or trimming
/// the bases the two share at the start and then at the end (Utils/Sequence.pm:965) leaves
/// none of the alternate's and whole codons of the reference's (:1181, :1185).
fn inframe_deletion_reads_codons(codons: &str) -> bool {
    let (reference, _, alternate) = pyfmt::partition(codons, "/");
    let bases = |half: &str| {
        if half == "-" {
            Vec::new()
        } else {
            half.to_ascii_uppercase().into_bytes()
        }
    };
    let (reference, alternate) = (bases(reference), bases(alternate));
    if alternate.len() >= reference.len() {
        return false;
    }
    if reference.starts_with(&alternate) || reference.ends_with(&alternate) {
        return true;
    }
    let (mut reference, mut alternate) = (reference.as_slice(), alternate.as_slice());
    while let (Some((first, rest)), Some((other, other_rest))) =
        (reference.split_first(), alternate.split_first())
    {
        if first != other {
            break;
        }
        (reference, alternate) = (rest, other_rest);
    }
    while let (Some((last, rest)), Some((other, other_rest))) =
        (reference.split_last(), alternate.split_last())
    {
        if last != other {
            break;
        }
        (reference, alternate) = (rest, other_rest);
    }
    alternate.is_empty() && reference.len().is_multiple_of(3)
}

/// Whether the Consequence `expected_value` of the reference entry `e` and `actual_value`
/// of the candidate entry `a` differ as the term pair of
/// `refseq_edited_model_deleted_span_codons`: the reference's terms less the candidate's are
/// `protein_altering_variant` alone and the candidate's less the reference's
/// `inframe_deletion` alone, the row's two Codons are a `deleted_span_codons_pair`, and
/// `inframe_deletion_reads_codons` refuses the reference's Codons and accepts the
/// candidate's, so `protein_altering_variant` (Utils/VariationEffect.pm:375-392) stands
/// where the in-frame test read the allele's bases.
pub fn deleted_span_inframe_term_pair(
    e: &Entry,
    a: &Entry,
    expected_value: &str,
    actual_value: &str,
) -> bool {
    let (expected_set, actual_set) = (
        consequence_set(expected_value),
        consequence_set(actual_value),
    );
    let expected_terms: BTreeSet<&str> = expected_set.split(',').collect();
    let actual_terms: BTreeSet<&str> = actual_set.split(',').collect();
    if !expected_terms
        .difference(&actual_terms)
        .eq([&PROTEIN_ALTERING_TERM])
        || !actual_terms
            .difference(&expected_terms)
            .eq([&INFRAME_DELETION_TERM])
    {
        return false;
    }
    let expected_codons = first_field(e, CODONS_FIELDS);
    let actual_codons = first_field(a, CODONS_FIELDS);
    deleted_span_codons_pair(e, a, expected_codons, actual_codons)
        && !inframe_deletion_reads_codons(expected_codons)
        && inframe_deletion_reads_codons(actual_codons)
}

/// Whether the differing field `name` of the reference entry `e` and the candidate entry
/// `a` of one key, with the values `expected_value` and `actual_value`, is a cell
/// `refseq_edited_model_deleted_span_codons` takes. Both values are written, and the field
/// is one of three:
///
/// - the Codons `deleted_span_codons_pair` accepts;
/// - the Amino_acids of a row whose two Codons are such a pair, each reference half the
///   translation of its own row's reference codons (`deleted_span_residues_pair`);
/// - the Consequence of a row whose two Codons are such a pair and whose two sets of terms
///   are the class's term pair (`deleted_span_inframe_term_pair`).
pub fn deleted_span_coding_reference_cell(
    e: &Entry,
    a: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
) -> bool {
    if blank(expected_value) || blank(actual_value) {
        return false;
    }
    if field_named(name, CODONS_FIELDS) {
        return deleted_span_codons_pair(e, a, expected_value, actual_value);
    }
    if field_named(name, AMINO_ACIDS_FIELDS) {
        return deleted_span_residues_pair(e, a, expected_value, actual_value);
    }
    if field_named(name, CONSEQUENCE_FIELDS) {
        return deleted_span_inframe_term_pair(e, a, expected_value, actual_value);
    }
    false
}

/// Whether the differing field `name` of the reference entry `e` and the candidate entry
/// `a` of one key, with the values `expected_value` and `actual_value`, is a cell
/// `refseq_edited_model_intron_span_coding_reference` takes. The row has the class's shape
/// (`intron_span_row`), both values are written, and the field is one of four:
///
/// - the Consequence whose two values are the class's term pair (`intron_span_term_pair`);
/// - the Codons `deleted_span_codons_pair` accepts: the alternate halves equal, the
///   reference halves of one length, the reference's the model's bases before the row's
///   first coding base and then, on the transcript's strand, the leading bases of USED_REF
///   or of a GIVEN_REF of its length;
/// - the Amino_acids of a row whose two Codons are such a pair: the alternate halves equal,
///   the reference halves of one length and different, each the translation of its own
///   row's reference codons (`deleted_span_residues_pair`);
/// - the HGVSp of a row whose two sets of terms are the class's term pair: the reference's
///   value the candidate's followed by `STOP_LOST_EXTENSION`.
pub fn intron_span_coding_reference_cell(
    e: &Entry,
    a: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
) -> bool {
    if blank(expected_value) || blank(actual_value) || !intron_span_row(e, a) {
        return false;
    }
    if field_named(name, CONSEQUENCE_FIELDS) {
        return intron_span_term_pair(expected_value, actual_value);
    }
    if field_named(name, CODONS_FIELDS) {
        return deleted_span_codons_pair(e, a, expected_value, actual_value);
    }
    if field_named(name, AMINO_ACIDS_FIELDS) {
        return deleted_span_residues_pair(e, a, expected_value, actual_value);
    }
    if field_named(name, HGVSP_FIELDS) {
        return intron_span_term_pair(&e.consequence, &a.consequence)
            && expected_value.strip_suffix(STOP_LOST_EXTENSION) == Some(actual_value);
    }
    false
}

/// Whether the reference entry `e` and the candidate entry `a` of one key have the
/// shape of a class whose one test Ensembl VEP reads at the unshifted span
/// (`mature_mirna_unshifted_span`, `stop_codon_unshifted_span`): a Transcript entry with a
/// positive `SHIFT_LENGTH`, which only an insertion or deletion is given, whose two
/// Consequence sets are equal but for `reference_term` on the reference's side against
/// `candidate_term` on the candidate's, or, under `either_way`, the two terms swapped. A
/// set differing by any other term, or by these two beside another, stays a difference.
pub fn unshifted_span_term_pair(
    e: &Entry,
    a: &Entry,
    reference_term: &str,
    candidate_term: &str,
    either_way: bool,
) -> Result<bool> {
    if !transcript_row(e) || !positive_shift(e)? {
        return Ok(false);
    }
    let (expected_set, actual_set) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    // A Consequence without a term is the one empty term.
    let expected_terms: BTreeSet<&str> = expected_set.split(',').collect();
    let actual_terms: BTreeSet<&str> = actual_set.split(',').collect();
    let reference_only: Vec<&str> = expected_terms.difference(&actual_terms).copied().collect();
    let candidate_only: Vec<&str> = actual_terms.difference(&expected_terms).copied().collect();
    if reference_only == [reference_term] && candidate_only == [candidate_term] {
        return Ok(true);
    }
    Ok(either_way && reference_only == [candidate_term] && candidate_only == [reference_term])
}

/// `mature_mirna_unshifted_span`: `mature_miRNA_variant` against
/// `non_coding_transcript_exon_variant` in either direction.
pub fn mature_mirna_unshifted_span_pair(e: &Entry, a: &Entry) -> Result<bool> {
    unshifted_span_term_pair(e, a, MATURE_MIRNA_TERM, NON_CODING_EXON_TERM, true)
}

/// `stop_codon_unshifted_span`: the reference's `coding_sequence_variant` against the
/// candidate's `stop_retained_variant`.
pub fn stop_codon_unshifted_span_pair(e: &Entry, a: &Entry) -> Result<bool> {
    unshifted_span_term_pair(e, a, CODING_SEQUENCE_TERM, STOP_RETAINED_TERM, false)
}

/// `shared::stop_codon_test_coordinates`, the test `tools/vep-outputs-compare` makes as well.
pub fn stop_codon_test_coordinates(e: &Entry) -> bool {
    shared::stop_codon_test_coordinates(e)
}

/// `shared::stop_codon_test_not_made`, the test `tools/vep-outputs-compare` makes as well.
pub fn stop_codon_test_not_made(e: &Entry) -> bool {
    shared::stop_codon_test_not_made(e)
}

/// `shared::stop_retained_without_stop_codon_test_pair`, which decides the class as
/// `tools/vep-outputs-compare` decides it.
pub fn stop_retained_without_stop_codon_test_pair(e: &Entry, a: &Entry) -> Result<bool> {
    shared::stop_retained_without_stop_codon_test_pair(e, a)
}

/// `shared::stop_codon_overlap_genomic_span_pair`, which decides the class as
/// `tools/vep-outputs-compare` decides it.
pub fn stop_codon_overlap_genomic_span_pair(e: &Entry, a: &Entry) -> Result<bool> {
    shared::stop_codon_overlap_genomic_span_pair(e, a)
}

/// `shared::stop_codon_overlap_insertion_range_pair`, which decides the class as
/// `tools/vep-outputs-compare` decides it.
pub fn stop_codon_overlap_insertion_range_pair(e: &Entry, a: &Entry) -> Result<bool> {
    shared::stop_codon_overlap_insertion_range_pair(e, a)
}

/// `shared::stop_retained_three_base_stop_codon_test_pair`, which decides the class as
/// `tools/vep-outputs-compare` decides it.
pub fn stop_retained_three_base_stop_codon_test_pair(e: &Entry, a: &Entry) -> Result<bool> {
    shared::stop_retained_three_base_stop_codon_test_pair(e, a)
}

/// `shared::stop_codon_removed_without_stop_codon_test_pair`, which decides the class as
/// `tools/vep-outputs-compare` decides it.
pub fn stop_codon_removed_without_stop_codon_test_pair(e: &Entry, a: &Entry) -> Result<bool> {
    shared::stop_codon_removed_without_stop_codon_test_pair(e, a)
}

/// The Allele both engines write on the row of a `<CNV:TR>` allele read symbolically.
pub const TANDEM_REPEAT_ALLELE: &str = "tandem_repeat";

pub const FEATURE_ELONGATION_TERM: &str = "feature_elongation";

pub const FEATURE_TRUNCATION_TERM: &str = "feature_truncation";

/// Whether `actual` is `expected` without its one `feature_elongation`, the terms left in
/// their order, and holds no `feature_truncation`.
fn names_no_length_term(expected: &[&str], actual: &[&str]) -> bool {
    let left: Vec<&str> = expected
        .iter()
        .copied()
        .filter(|term| *term != FEATURE_ELONGATION_TERM)
        .collect();
    left.len() + 1 == expected.len() && left == actual && !actual.contains(&FEATURE_TRUNCATION_TERM)
}

/// `cnv_tr_symbolic_loss_named_a_gain`: a Transcript row whose Allele is `tandem_repeat` on
/// both sides and whose two lists of terms, each in the order its side writes it, differ in
/// the reference's one `feature_elongation` alone. Either the two lists are of one length
/// and the candidate writes `feature_truncation` in that place, both terms carrying the
/// impact HIGH and neighbouring ranks (Utils/Constants.pm:803-834), so the one stands where
/// the other stood. Or the candidate writes no term in that place
/// (`names_no_length_term`), and its IMPACT, where the row carries one, is the one its own
/// terms give (`terms_impact`).
pub fn symbolic_loss_named_a_gain_pair(e: &Entry, a: &Entry) -> Result<bool> {
    if !transcript_row(e) || e.allele != TANDEM_REPEAT_ALLELE || a.allele != TANDEM_REPEAT_ALLELE {
        return Ok(false);
    }
    // The text formats join a row's terms with `,` and the VCF format with `&`.
    let expected: Vec<&str> = e.consequence.split(['&', ',']).collect();
    let actual: Vec<&str> = a.consequence.split(['&', ',']).collect();
    if names_no_length_term(&expected, &actual) {
        let impact = first_field(a, IMPACT_FIELDS);
        return Ok(blank(impact) || terms_impact(&actual).as_deref() == Some(impact));
    }
    if expected.len() != actual.len() || actual.contains(&FEATURE_ELONGATION_TERM) {
        return Ok(false);
    }
    let mut named_a_gain = 0;
    for (theirs, ours) in expected.into_iter().zip(actual) {
        if theirs == ours {
            continue;
        }
        if theirs != FEATURE_ELONGATION_TERM || ours != FEATURE_TRUNCATION_TERM {
            return Ok(false);
        }
        named_a_gain += 1;
    }
    Ok(named_a_gain == 1)
}

/// `at_its_location` without a copy of an entry it returns as it is.
fn located(e: &Entry) -> Cow<'_, Entry> {
    match (&e.location, &e.allele_location) {
        (None, Some(location)) => Cow::Owned(Entry {
            location: Some(location.clone()),
            ..e.clone()
        }),
        _ => Cow::Borrowed(e),
    }
}

/// A VCF or JSON entry with the Location its record gives its allele
/// (`Entry.allele_location`) as its `location`, the entry a predicate over a default
/// row's Location reads; any other entry as it is.
pub fn at_its_location(e: &Entry) -> Entry {
    located(e).into_owned()
}

/// `shared::splice_term_past_the_record_span_pair`, which decides the class as
/// `tools/vep-outputs-compare` decides it, on the two entries at their Locations
/// (`at_its_location`): a VCF or JSON entry whose record gives its allele no Location is refused.
pub fn splice_term_past_the_record_span_pair(e: &Entry, a: &Entry) -> Result<bool> {
    shared::splice_term_past_the_record_span_pair(&located(e), &located(a))
}

/// The row rules of the classes of a coding end and of the stop codon classes, whose terms
/// a row can differ by beside the six of `splice_term_past_the_record_span`: one row can
/// carry both differences.
pub const CODING_END_ROW_PAIRS: &[RowPair] = &[
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

/// `splice_term_past_the_record_span_pair` in a comparison that may have been given the
/// candidate's cache: a row that differs by the terms of a class of a coding end or of a
/// stop codon class as well is set aside where that class's rule (`CODING_END_ROW_PAIRS`)
/// takes the rest of the difference on what `at` holds of the row's transcript.
pub fn splice_term_past_the_record_span_pair_at(
    e: &Entry,
    a: &Entry,
    at: &AtACodingEnd<'_>,
) -> Result<bool> {
    shared::splice_term_past_the_record_span_pair_beside(
        &located(e),
        &located(a),
        &|beside: &Entry, a: &Entry| {
            for pair in CODING_END_ROW_PAIRS {
                if pair(beside, a, at)? {
                    return Ok(true);
                }
            }
            Ok(false)
        },
    )
}

/// The entries of `entries` whose place `asks` marks.
fn asking_of(entries: &mut [Entry], asks: Vec<bool>) -> impl Iterator<Item = &mut Entry> {
    entries
        .iter_mut()
        .zip(asks)
        .filter_map(|(entry, asks)| asks.then_some(entry))
}

/// Gives the default and tab entries of a key whose two Consequence texts differ, and those
/// of a key whose two protein notations are one frameshift named at residue 1 with the
/// unknown count and with a number (`count_alone_differs_at_the_first_residue`), where
/// they carry no reference bases (`Entry.reference`), the bases their allele stands in
/// for, read by position from the records of the input VCF `path` as the tuple comparator
/// reads them (`read_record_references` of `tools/vep-outputs-compare`), and for an entry of
/// an insertion or deletion that the 3' rule moved off its record's alleles from the
/// record whose identifier is the entry's own first column (`moved_reference`). An entry no
/// record states keeps `None`; where two or more records read an entry differently, the
/// record whose
/// identifier is the entry's own first column decides, and the entry keeps `None` where no
/// such record reads it or such records read it differently. An entry whose Location names
/// no span of whole positions (`location_span`), which no record gives, asks for none.
pub fn read_entry_references(
    expected: &mut [Entry],
    actual: &mut [Entry],
    path: &Path,
) -> Result<()> {
    // Of two candidate entries under one key the last is the key's.
    let actual_by_key: HashMap<&str, usize> = actual
        .iter()
        .enumerate()
        .map(|(index, a)| (a.key.as_str(), index))
        .collect();
    let mut expected_asks = vec![false; expected.len()];
    let mut actual_asks = vec![false; actual.len()];
    for (index, e) in expected.iter().enumerate() {
        let Some(&paired) = actual_by_key.get(e.key.as_str()) else {
            continue;
        };
        let a = &actual[paired];
        let counted_at_the_first_residue = || {
            count_alone_differs_at_the_first_residue(
                first_field(e, HGVSP_FIELDS),
                first_field(a, HGVSP_FIELDS),
            )
        };
        if e.consequence != a.consequence || counted_at_the_first_residue() {
            expected_asks[index] = location_span(e)?.is_some();
            // A candidate entry two reference entries pair with asks once.
            actual_asks[paired] = location_span(a)?.is_some() || actual_asks[paired];
        }
    }
    drop(actual_by_key);
    let mut asking: Vec<&mut Entry> = asking_of(expected, expected_asks)
        .chain(asking_of(actual, actual_asks))
        .collect();
    read_record_references_of(path, &mut asking)
}

/// The `population:allele` keys of a `FREQS` value with their frequencies, or `None`
/// when an item is not a `population:allele:value` triple; the items are joined with `,`
/// in the default and tab formats and `&` in the VCF (`OutputFactory.pm:1152-1159`).
pub fn freqs_items(value: &str) -> Option<Fields> {
    let mut items = Fields::new();
    for item in value.split([',', '&']) {
        let (population, _, rest) = pyfmt::partition(item, ":");
        let (allele, separator, frequency) = pyfmt::partition(rest, ":");
        if separator.is_empty() {
            return None;
        }
        items.insert(format!("{population}:{allele}"), frequency.to_string());
    }
    Some(items)
}

/// Whether two `FREQS` values carry the same items and differ only where the reference
/// reads `NA` and the candidate `0` (`frequency_filter_zero_value`).
pub fn freqs_differ_by_zero_read_as_na(expected_value: &str, actual_value: &str) -> bool {
    if blank(expected_value) || blank(actual_value) {
        return false;
    }
    let (Some(expected_items), Some(actual_items)) =
        (freqs_items(expected_value), freqs_items(actual_value))
    else {
        return false;
    };
    if expected_items.len() != actual_items.len()
        || !expected_items
            .keys()
            .all(|key| actual_items.contains_key(key))
    {
        return false;
    }
    let mut changed = expected_items
        .iter()
        .filter(|(key, value)| actual_items.get(*key) != Some(*value))
        .peekable();
    changed.peek().is_some()
        && changed.all(|(key, value)| {
            value == "NA" && actual_items.get(key).and_then(|value| as_number(value)) == Some(0.0)
        })
}

/// Whether the entry's `FREQS` carries at least one item and every item's frequency is
/// 0: the rows of a record `--freq_filter include` keeps for a zero frequency.
pub fn freqs_all_zero(e: &Entry) -> bool {
    let value = first_field(e, FREQS_FIELDS);
    if blank(value) {
        return false;
    }
    freqs_items(value).is_some_and(|items| {
        !items.is_empty() && items.values().all(|value| as_number(value) == Some(0.0))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Uncaught;

    fn entry(fields: &[(&str, &str)]) -> Entry {
        Entry {
            key: "21:100|A|ENST1|Transcript".to_string(),
            location: Some("21:100".to_string()),
            allele: "A".to_string(),
            feature: "ENST1".to_string(),
            consequence: "intron_variant".to_string(),
            fields: fields
                .iter()
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect(),
            ..Entry::default()
        }
    }

    /// What a run that ends on an input no reading of the comparison expects says.
    fn uncaught<T>(run: impl FnOnce() -> T + std::panic::UnwindSafe) -> String {
        let payload = std::panic::catch_unwind(run).err().expect("the run ends");
        let Uncaught(what) = payload.downcast_ref::<Uncaught>().expect("an exception");
        what.clone()
    }

    #[test]
    fn a_span_in_cdna_numbering_is_held_to_the_rows_position_under_a_printed_shift() {
        let matches = |fields: &[(&str, &str)], start: &str, end: &str| {
            span_matches_cdna_position(&entry(fields), "n", start, end).expect("no coordinate")
        };
        let shifted = |position: &'static str| [("SHIFT_LENGTH", "3"), ("cDNA_position", position)];
        assert!(matches(&[("cDNA_position", "23-43")], "24", "44"));
        assert!(matches(&shifted("23-43"), "23", "43"));
        assert!(!matches(&shifted("23-43"), "24", "43"));
        assert!(matches(&shifted("?-43"), "24", "43"));
        assert!(!matches(&shifted("23-43/1500"), "23", "44"));
        // The two sides are compared as they are written, and one position is both sides.
        assert!(!matches(&shifted("048"), "48", "48"));
        assert!(!matches(&shifted("48"), "48", "49"));
        // A position that ends in its separator has an empty end.
        assert!(matches(&shifted("23-"), "23", "44"));
        assert!(!matches(&shifted("\u{b2}-5"), "3", "5"));
        let json = [
            ("shift_length", "3"),
            ("cdna_start", "23"),
            ("cdna_end", "43"),
        ];
        assert!(!matches(&json, "23", "44"));
        assert!(matches(&json[..2], "23", "44"));
        assert!(!matches(
            &[("shift_length", "0"), ("cdna_end", "43")],
            "23",
            "44"
        ));
    }

    #[test]
    fn a_difference_of_whole_numbers_is_exact_at_any_length() {
        assert_eq!(difference("43", "23"), 20);
        assert_eq!(difference("23", "43"), -20);
        assert_eq!(difference("0", "0"), 0);
        let wide = format!("1{}", "0".repeat(41));
        let wider = format!("1{}20", "0".repeat(39));
        assert_eq!(difference(&wider, &wide), 20);
        assert_eq!(difference(&wide, &wider), -20);
        assert_eq!(difference(&wide, &wide), 0);
        assert_eq!(difference(&wide, "5"), i128::MAX);
        assert_eq!(difference("5", &wide), i128::MIN);
        let past = "170141183460469231731687303715884105728";
        assert_eq!(difference(past, "1"), i128::MAX);
        assert_eq!(difference("1", past), -i128::MAX);
    }

    #[test]
    fn an_acceptor_sites_offsets_are_those_of_the_coordinates_that_place_the_span() {
        let offsets = |body: &str| acceptor_site_offsets(body);
        let some = |offsets: &[&str]| Some(offsets.iter().map(|o| o.to_string()).collect());
        assert_eq!(offsets("c.5-1_5insA"), some(&["-1"]));
        assert_eq!(offsets("c.5-2_5-1insACG"), some(&["-2", "-1"]));
        // A duplication is placed by its last coordinate alone.
        assert_eq!(offsets("c.5-3_5-2dup"), some(&["-2"]));
        assert_eq!(offsets("c.5-4_5-1dup"), some(&["-1"]));
        assert_eq!(offsets("c.5-1dup"), some(&["-1"]));
        assert_eq!(offsets("c.5-3dup"), None);
        assert_eq!(offsets("c.5+1dup"), None);
        assert_eq!(offsets("c.5_6del"), some(&[]));
        assert_eq!(offsets("c.5del"), some(&[]));
        assert_eq!(offsets("c.5_6insA"), some(&[]));
        assert_eq!(offsets("c.5-2_5del"), some(&["-2"]));
        assert_eq!(offsets("c.-5-1_-5del"), some(&["-1"]));
        assert_eq!(offsets("c.*5-2_*5del"), some(&["-2"]));
        assert_eq!(offsets("c.5-3_5del"), None);
        assert_eq!(offsets("c.5-1_6+1del"), None);
        // An offset is one of the two as they are written: another script's digit is none.
        assert_eq!(offsets("c.5-\u{661}del"), None);
        for other in ["c.5delins", "c.5>T", "n.5-1del\n", ""] {
            assert_eq!(offsets(other), None, "{other:?}");
        }
    }

    #[test]
    fn the_coding_numbering_steps_over_its_missing_zero() {
        for (position, step, stepped) in [
            (-2, 1, -1),
            (-1, 1, 1),
            (1, -1, -1),
            (1, 1, 2),
            (5, -1, 4),
            (0, 1, 1),
            (0, -1, -1),
            (i64::MAX, 1, i64::MAX),
            (i64::MIN, -1, i64::MIN),
        ] {
            assert_eq!(cdna_step(position, step), stepped, "{position} {step}");
        }
    }

    #[test]
    fn a_clamped_insertion_and_a_shifted_one_are_read_from_their_notations() {
        let clamped = |body: &str| clamped_insertion_position(body);
        let shifted = |body: &str| shifted_insertion_end(body);
        let at = |numbering: &str, position: i64, bases: &str, stripped: i64| {
            Some((numbering.to_string(), position, bases.to_string(), stripped))
        };
        let end = |numbering: &str, position: i64, bases: Option<&str>, count: i64| {
            Some((
                numbering.to_string(),
                position,
                bases.map(str::to_string),
                count,
            ))
        };
        assert_eq!(clamped("c.5A>T"), at("c", 5, "T", 0));
        assert_eq!(clamped("n.-5delinsACG"), at("n", -5, "ACG", 0));
        assert_eq!(clamped("c.5N>A"), at("c", 5, "A", 0));
        assert_eq!(clamped("c.5_6insA"), at("c", 6, "A", 1));
        assert_eq!(shifted("c.5_6insA"), end("c", 6, Some("A"), 1));
        assert_eq!(clamped("n.1_2insTTGA"), at("n", 2, "TTGA", 1));
        assert_eq!(shifted("n.1_2insTTGA"), end("n", 2, Some("TTGA"), 4));
        // -1 and 1 are neighbours, and a duplication across them holds one base fewer.
        assert_eq!(clamped("c.-1_1insA"), at("c", 1, "A", 1));
        assert_eq!(clamped("c.-5_-4insA"), at("c", -4, "A", 1));
        assert_eq!(shifted("c.-1_1dup"), end("c", 2, None, 2));
        assert_eq!(shifted("c.-2_3dup"), end("c", 4, None, 5));
        assert_eq!(shifted("c.4_6dup"), end("c", 7, None, 3));
        assert_eq!(shifted("c.6dup"), end("c", 7, None, 1));
        // A digit of another script is a digit after the first.
        assert_eq!(clamped("c.1\u{665}A>T"), at("c", 15, "T", 0));
        for other in [
            "c.5_7insA",
            "c.5del",
            "c.5_6insA\n",
            "c.05A>T",
            "c.5A>N",
            "g.5A>T",
        ] {
            assert_eq!(clamped(other), None, "{other:?}");
            assert_eq!(shifted(other), None, "{other:?}");
        }
        assert_eq!(shifted("c.5A>T"), None);
        assert_eq!(clamped("c.4_6dup"), None);
        let last = i64::MAX;
        assert_eq!(
            clamped(&format!("c.{}_{last}insA", last - 1)),
            at("c", last, "A", 1)
        );
        // A position `i64` does not hold is no notation of either.
        assert_eq!(clamped("c.9223372036854775808A>T"), None);
        assert_eq!(shifted(&format!("c.{last}dup")), None);
        let long = "7".repeat(4301);
        for body in [format!("c.{long}A>T"), format!("c.5_{long}insA")] {
            let ending = uncaught(|| clamped_insertion_position(&body));
            assert!(
                ending.starts_with("Exceeds the limit (4300 digits)"),
                "{ending}"
            );
        }
        let ending = uncaught(|| shifted_insertion_end(&format!("c.{long}dup")));
        assert!(
            ending.starts_with("Exceeds the limit (4300 digits)"),
            "{ending}"
        );
    }

    #[test]
    fn the_coding_bounds_are_read_from_one_position_or_from_two_keys() {
        let bounds = |fields: &[(&str, &str)]| {
            let e = entry(fields);
            (cds_bounds(&e), cds_span_length(&e))
        };
        let span = (Some((3371, 3378)), Some(8));
        assert_eq!(bounds(&[("CDS_position", "3371-3378")]), span);
        assert_eq!(bounds(&[("CDS_position", "3371-3378/4551")]), span);
        assert_eq!(
            bounds(&[("cds_start", "3371/4551"), ("cds_end", "3378/4551")]),
            span
        );
        let one = (Some((3371, 3371)), Some(1));
        assert_eq!(bounds(&[("CDS_position", "3371")]), one);
        // A position that ends in its separator ends where it starts.
        assert_eq!(bounds(&[("CDS_position", "3371-")]), one);
        assert_eq!(bounds(&[("CDS_position", "5-3")]), (Some((5, 3)), Some(-1)));
        assert_eq!(
            bounds(&[("CDS_position", "\u{663}-\u{665}")]),
            (Some((3, 5)), Some(3))
        );
        for none in [
            &[("CDS_position", "?-3378")][..],
            &[("CDS_position", "-")],
            &[("cds_start", "3371")],
            &[],
            // The position column is read where the row has one, whatever the keys say.
            &[("CDS_position", "-"), ("cds_start", "1"), ("cds_end", "2")],
            // A bound `i64` does not hold is no bound.
            &[("CDS_position", "9223372036854775808")],
        ] {
            assert_eq!(bounds(none), (None, None), "{none:?}");
        }
        let superscript = entry(&[("CDS_position", "\u{b2}")]);
        assert_eq!(
            uncaught(|| cds_bounds(&superscript)),
            "invalid literal for int() with base 10: '\u{b2}'"
        );
    }

    #[test]
    fn a_vcf_or_json_entry_is_read_at_the_location_of_its_allele() {
        let placed = |location: Option<&str>, allele_location: Option<&str>| {
            let e = Entry {
                location: location.map(str::to_string),
                allele_location: allele_location.map(str::to_string),
                ..entry(&[("IMPACT", "LOW")])
            };
            let at = at_its_location(&e);
            assert_eq!(
                Entry {
                    location: e.location.clone(),
                    ..at.clone()
                },
                e
            );
            at.location
        };
        assert_eq!(placed(None, Some("7:5-6")).as_deref(), Some("7:5-6"));
        assert_eq!(placed(Some("21:100"), None).as_deref(), Some("21:100"));
        assert_eq!(
            placed(Some("21:100"), Some("7:9")).as_deref(),
            Some("21:100")
        );
        assert_eq!(placed(None, None), None);
    }

    #[test]
    fn the_items_of_a_freqs_value_are_keyed_by_population_and_allele() {
        let items = |value: &str| {
            freqs_items(value).map(|items| items.into_iter().collect::<Vec<(String, String)>>())
        };
        let pairs = |pairs: &[(&str, &str)]| {
            Some(
                pairs
                    .iter()
                    .map(|(key, value)| (key.to_string(), value.to_string()))
                    .collect::<Vec<_>>(),
            )
        };
        assert_eq!(items("GNOMADE_AFR:C:NA"), pairs(&[("GNOMADE_AFR:C", "NA")]));
        assert_eq!(
            items("GNOMADE_AFR:C:NA,AF:C:0.2"),
            pairs(&[("GNOMADE_AFR:C", "NA"), ("AF:C", "0.2")])
        );
        assert_eq!(items("a:b:c&d:e:f"), pairs(&[("a:b", "c"), ("d:e", "f")]));
        // What follows the second colon is the frequency, and a key written twice has its
        // last frequency.
        assert_eq!(items("a:b:c:d"), pairs(&[("a:b", "c:d")]));
        assert_eq!(items("a:b:1,a:b:2"), pairs(&[("a:b", "2")]));
        assert_eq!(items("a:b:"), pairs(&[("a:b", "")]));
        assert_eq!(items(":b:c"), pairs(&[(":b", "c")]));
        assert_eq!(items("a:b"), None);
        assert_eq!(items(""), None);
    }

    #[test]
    fn two_freqs_values_differ_by_a_zero_the_reference_reads_as_na() {
        for (expected, actual, taken) in [
            ("p:C:NA", "p:C:0", true),
            ("p:C:NA", "p:C:0.0", true),
            ("p:C:NA", "p:C:-0", true),
            ("p:C:NA", "p:C:\u{660}", true),
            ("p:C:NA", "p:C:0.01", false),
            ("p:C:NA,q:C:0.2", "p:C:0,q:C:0.2", true),
            ("p:C:NA,q:C:0.2", "q:C:0.2,p:C:0", true),
            ("p:C:NA,q:C:0.2", "p:C:0,q:C:0.3", false),
            ("p:C:0", "p:C:NA", false),
            ("p:C:NA", "p:C:NA", false),
            ("p:C:NA", "p:C:0&q:C:0", false),
            ("-", "p:C:0", false),
            ("p:C:NA", "p:C", false),
        ] {
            assert_eq!(
                freqs_differ_by_zero_read_as_na(expected, actual),
                taken,
                "{expected} {actual}"
            );
        }
    }

    #[test]
    fn a_row_kept_for_a_zero_frequency_has_every_item_at_zero() {
        for (fields, zero) in [
            (&[("FREQS", "p:C:0")][..], true),
            (&[("FREQS", "p:C:0,q:C:0.0")], true),
            (&[("freqs", "p:C:0")], true),
            (&[("FREQS", "p:C:0,q:C:0.1")], false),
            (&[("FREQS", "p:C:NA")], false),
            // The first of the two names the row carries is the one read.
            (&[("FREQS", "-"), ("freqs", "p:C:0")], false),
            (&[], false),
            (&[("FREQS", "p:C")], false),
            (&[("FREQS", "")], false),
        ] {
            assert_eq!(freqs_all_zero(&entry(fields)), zero, "{fields:?}");
        }
    }

    const INTRON_SPAN_REFERENCE: &str =
        "splice_acceptor_variant,splice_donor_variant,stop_lost,protein_altering_variant,intron_variant";
    const INTRON_SPAN_CANDIDATE: &str =
        "splice_acceptor_variant,splice_donor_variant,inframe_deletion,intron_variant";

    /// The two entries of a made-up row of `refseq_edited_model_intron_span_coding_reference`,
    /// the reference's first: a deletion of the coding bases 5 to 7 of a forward-strand model
    /// whose intron of six bases follows coding base 5. The codon window is the coding bases 4
    /// to 9, and the nine bases of the allele are coding base 5, the intron and coding bases 6
    /// and 7.
    fn intron_span_entries() -> (Entry, Entry) {
        let row = |terms: &str, codons: &str, amino_acids: &str, hgvsp: &str| Entry {
            key: "7:1005-1013|-|NM_990005.1|Transcript".to_string(),
            location: Some("7:1005-1013".to_string()),
            allele: "-".to_string(),
            feature: "NM_990005.1".to_string(),
            consequence: terms.to_string(),
            ..entry(&[
                ("Consequence", terms),
                ("CDS_position", "5-7"),
                ("Amino_acids", amino_acids),
                ("Codons", codons),
                ("STRAND", "1"),
                ("USED_REF", "AGTAAAGTT"),
                ("BAM_EDIT", "OK"),
                ("HGVSp", hgvsp),
            ])
        };
        (
            row(
                INTRON_SPAN_REFERENCE,
                "cAGTAA/cca",
                "Q*/P",
                "NP_990005.1:p.His2_Ser3delextTer?",
            ),
            row(
                INTRON_SPAN_CANDIDATE,
                "cATTCA/cca",
                "HS/P",
                "NP_990005.1:p.His2_Ser3del",
            ),
        )
    }

    /// `entry` with its field `name` set to `value`, or without the field for `None`.
    fn with_field(entry: &Entry, name: &str, value: Option<&str>) -> Entry {
        let mut changed = entry.clone();
        match value {
            Some(value) => {
                changed.fields.insert(name.to_string(), value.to_string());
            }
            None => {
                changed.fields.shift_remove(name);
            }
        }
        changed
    }

    /// `entry` with the terms `terms`, as its own and in its Consequence field.
    fn with_terms(entry: &Entry, terms: &str) -> Entry {
        Entry {
            consequence: terms.to_string(),
            ..with_field(entry, "Consequence", Some(terms))
        }
    }

    /// Whether the class takes the field `name` of the two entries, each with its own value.
    fn intron_span_cell(e: &Entry, a: &Entry, name: &str) -> bool {
        let value = |entry: &Entry| {
            entry
                .fields
                .get(name)
                .map_or("-", String::as_str)
                .to_string()
        };
        intron_span_coding_reference_cell(e, a, name, &value(e), &value(a))
    }

    const INTRON_SPAN_CELLS: [&str; 4] = ["Consequence", "Codons", "Amino_acids", "HGVSp"];

    #[test]
    fn a_row_of_the_intron_span_class_is_an_edited_models_with_a_longer_allele_and_the_intron_terms(
    ) {
        let (e, a) = intron_span_entries();
        assert!(intron_span_row(&e, &a));
        let both = |name: &str, value: Option<&str>| {
            intron_span_row(&with_field(&e, name, value), &with_field(&a, name, value))
        };
        let reference_alone =
            |name: &str, value: Option<&str>| intron_span_row(&with_field(&e, name, value), &a);
        let candidate_alone =
            |name: &str, value: Option<&str>| intron_span_row(&e, &with_field(&a, name, value));
        assert!(!reference_alone("BAM_EDIT", None));
        assert!(!candidate_alone("BAM_EDIT", None));
        assert!(!both("BAM_EDIT", Some("FAILED")));
        // The coding span is one range of numbers on both sides, read past a /total suffix.
        assert!(both("CDS_position", Some("5-7/30")));
        assert!(!candidate_alone("CDS_position", Some("5-8")));
        for no_range in ["?-7", "5-?", "5", "7-5", "-"] {
            assert!(!both("CDS_position", Some(no_range)), "{no_range}");
        }
        // The allele is one text of letters on both sides, longer than the range.
        assert!(both("USED_REF", Some("AGTA")));
        assert!(!candidate_alone("USED_REF", Some("AGTAAAGTA")));
        for no_longer_allele in [Some("AGT"), Some("-"), Some("AGT.AAGTT"), None] {
            assert!(!both("USED_REF", no_longer_allele), "{no_longer_allele:?}");
        }
        // Each side holds the three terms of an intron inside the span.
        for term in INTRON_SPAN_TERMS {
            let without = |terms: &str| {
                let kept: Vec<&str> = terms.split(',').filter(|held| held != term).collect();
                kept.join(",")
            };
            let reference = with_terms(&e, &without(INTRON_SPAN_REFERENCE));
            let candidate = with_terms(&a, &without(INTRON_SPAN_CANDIDATE));
            assert!(!intron_span_row(&reference, &a), "{term}");
            assert!(!intron_span_row(&e, &candidate), "{term}");
            assert!(!intron_span_row(&reference, &candidate), "{term}");
        }
        let of_kind = |entry: &Entry, kind: &str| Entry {
            key: entry.key.replace("|Transcript", kind),
            ..entry.clone()
        };
        assert!(!intron_span_row(
            &of_kind(&e, "|RegulatoryFeature"),
            &of_kind(&a, "|RegulatoryFeature")
        ));
        let structural = |entry: &Entry| Entry {
            allele: "deletion".to_string(),
            ..entry.clone()
        };
        assert!(!intron_span_row(&structural(&e), &structural(&a)));
        // The JSON format writes the span under two keys and every key in lower case.
        let json = |entry: &Entry, end: &str| {
            let mut fields = Fields::new();
            for (name, value) in [
                ("cds_start", "5"),
                ("cds_end", end),
                ("bam_edit", "OK"),
                ("used_ref", "AGTAAAGTT"),
            ] {
                fields.insert(name.to_string(), value.to_string());
            }
            Entry {
                fields,
                ..entry.clone()
            }
        };
        assert!(intron_span_row(&json(&e, "7"), &json(&a, "7")));
        assert!(!intron_span_row(&json(&e, "5"), &json(&a, "5")));
        assert!(!intron_span_row(&json(&e, "7"), &json(&a, "8")));
    }

    #[test]
    fn the_intron_span_class_takes_the_four_cells_of_its_row_and_no_other() {
        let (e, a) = intron_span_entries();
        for name in INTRON_SPAN_CELLS {
            assert!(intron_span_cell(&e, &a, name), "{name}");
            // The two engines swapped are another difference.
            assert!(!intron_span_cell(&a, &e, name), "{name}");
            // A value one side leaves out is no pair.
            let left_out = with_field(&a, name, None);
            assert!(!intron_span_cell(&e, &left_out, name), "{name}");
        }
        for (name, text_name) in [
            ("consequence_terms", "Consequence"),
            ("codons", "Codons"),
            ("amino_acids", "Amino_acids"),
            ("hgvsp", "HGVSp"),
        ] {
            let (expected, actual) = (&e.fields[text_name], &a.fields[text_name]);
            assert!(
                intron_span_coding_reference_cell(&e, &a, name, expected, actual),
                "{name}"
            );
        }
        for (name, expected, actual) in [
            ("IMPACT", "HIGH", "MODERATE"),
            ("CDS_position", "5-7", "5-8"),
            ("HGVSc", "NM_990005.1:c.5_7del", "NM_990005.1:c.5_8del"),
            ("USED_REF", "AGTAAAGTT", "AGTAAAGTA"),
            ("Extra", "IMPACT=HIGH", "IMPACT=MODERATE"),
        ] {
            assert!(
                !intron_span_coding_reference_cell(&e, &a, name, expected, actual),
                "{name}"
            );
        }
        // On a row without the shape each of the four cells stays a difference.
        let (unedited_e, unedited_a) = (
            with_field(&e, "BAM_EDIT", None),
            with_field(&a, "BAM_EDIT", None),
        );
        let (no_intron_e, no_intron_a) = (
            with_terms(
                &e,
                "splice_donor_variant,stop_lost,protein_altering_variant",
            ),
            with_terms(&a, "splice_donor_variant,inframe_deletion"),
        );
        for name in INTRON_SPAN_CELLS {
            assert!(!intron_span_cell(&unedited_e, &unedited_a, name), "{name}");
            assert!(
                !intron_span_cell(&no_intron_e, &no_intron_a, name),
                "{name}"
            );
        }
    }

    #[test]
    fn the_intron_span_class_takes_one_pair_of_term_sets_and_the_notation_beside_it() {
        let (e, a) = intron_span_entries();
        let taken = |reference: &str, candidate: &str| {
            let (e, a) = (with_terms(&e, reference), with_terms(&a, candidate));
            (
                intron_span_cell(&e, &a, "Consequence"),
                intron_span_cell(&e, &a, "HGVSp"),
            )
        };
        let fifth_base = |terms: &str| format!("{terms},splice_donor_5th_base_variant");
        // A term both sides hold beside the pair, the order of the terms and the VCF format's
        // separator change nothing.
        assert_eq!(
            taken(
                &fifth_base(INTRON_SPAN_REFERENCE),
                &fifth_base(INTRON_SPAN_CANDIDATE)
            ),
            (true, true)
        );
        assert_eq!(
            taken(
                "intron_variant,stop_lost,splice_donor_variant,protein_altering_variant,splice_acceptor_variant",
                INTRON_SPAN_CANDIDATE
            ),
            (true, true)
        );
        assert_eq!(
            taken(
                &INTRON_SPAN_REFERENCE.replace(',', "&"),
                &INTRON_SPAN_CANDIDATE.replace(',', "&")
            ),
            (true, true)
        );
        let frameshift = INTRON_SPAN_CANDIDATE.replace("inframe_deletion", "frameshift_variant");
        for (reference, candidate) in [
            // One of the reference's two terms alone.
            (
                INTRON_SPAN_REFERENCE.replace("stop_lost,", ""),
                INTRON_SPAN_CANDIDATE.to_string(),
            ),
            (
                INTRON_SPAN_REFERENCE.replace("protein_altering_variant,", ""),
                INTRON_SPAN_CANDIDATE.to_string(),
            ),
            // A further term on one side alone.
            (
                fifth_base(INTRON_SPAN_REFERENCE),
                INTRON_SPAN_CANDIDATE.to_string(),
            ),
            (
                INTRON_SPAN_REFERENCE.to_string(),
                fifth_base(INTRON_SPAN_CANDIDATE),
            ),
            // Another term in the place of the candidate's.
            (INTRON_SPAN_REFERENCE.to_string(), frameshift),
            // One set on both sides.
            (
                INTRON_SPAN_REFERENCE.to_string(),
                INTRON_SPAN_REFERENCE.to_string(),
            ),
            (
                INTRON_SPAN_CANDIDATE.to_string(),
                INTRON_SPAN_CANDIDATE.to_string(),
            ),
        ] {
            assert_eq!(
                taken(&reference, &candidate),
                (false, false),
                "{reference} {candidate}"
            );
        }
        // The notation is the candidate's followed by the extension, and nothing else.
        let notation = |reference: &str, candidate: &str| {
            intron_span_coding_reference_cell(&e, &a, "HGVSp", reference, candidate)
        };
        let deletion = "NP_990005.1:p.His2_Ser3del";
        assert!(notation(&format!("{deletion}extTer?"), deletion));
        for (reference, candidate) in [
            (format!("{deletion}extTer4"), deletion),
            (format!("{deletion}insProextTer?"), deletion),
            (format!("{deletion}extTer?"), "NP_990005.1:p.His2del"),
            (deletion.to_string(), deletion),
            ("-extTer?".to_string(), "-"),
            ("extTer?".to_string(), ""),
        ] {
            assert!(!notation(&reference, candidate), "{reference} {candidate}");
        }
    }

    #[test]
    fn the_intron_span_class_holds_each_sides_amino_acids_to_its_own_reference_codons() {
        let (e, a) = intron_span_entries();
        let codons = |reference: &str, candidate: &str| {
            let e = with_field(&e, "Codons", Some(reference));
            let a = with_field(&a, "Codons", Some(candidate));
            (
                intron_span_cell(&e, &a, "Codons"),
                intron_span_cell(&e, &a, "Amino_acids"),
            )
        };
        assert_eq!(codons("cAGTAA/cca", "cATTCA/cca"), (true, true));
        for (reference, candidate) in [
            // Another alternate half, another length, and one value on both sides.
            ("cAGTAA/cca", "cATTCA/ccc"),
            ("cAGTAA/cca", "cATTCAGGG/cca"),
            ("cAGTAA/cca", "cAGTAA/cca"),
            // Upper-case bases that are not the allele's leading bases.
            ("cAGTAC/cca", "cATTCA/cca"),
            // The model's base before the span differs between the two sides.
            ("gAGTAA/cca", "cATTCA/cca"),
        ] {
            assert_eq!(
                codons(reference, candidate),
                (false, false),
                "{reference} {candidate}"
            );
        }
        let amino_acids = |reference: &str, candidate: &str| {
            intron_span_coding_reference_cell(&e, &a, "Amino_acids", reference, candidate)
        };
        assert!(amino_acids("Q*/P", "HS/P"));
        for (reference, candidate) in [
            // Another alternate half.
            ("Q*/P", "HS/L"),
            // Residues that are not the translation of the side's own reference codons.
            ("QW/P", "HS/P"),
            ("Q*/P", "HT/P"),
            // Reference halves of two lengths, and one reference half on both sides.
            ("Q*/P", "H/P"),
            ("Q*/P", "Q*/P"),
            ("-", "HS/P"),
        ] {
            assert!(
                !amino_acids(reference, candidate),
                "{reference} {candidate}"
            );
        }
        // A codon the standard code does not read ties no residue to it.
        let unread = |entry: &Entry, codons: &str| {
            with_field(
                &with_field(entry, "USED_REF", Some("AGTANAGTT")),
                "Codons",
                Some(codons),
            )
        };
        let (unread_e, unread_a) = (unread(&e, "cAGTAN/cca"), unread(&a, "cATTCA/cca"));
        assert!(intron_span_cell(&unread_e, &unread_a, "Codons"));
        assert!(!intron_span_coding_reference_cell(
            &unread_e,
            &unread_a,
            "Amino_acids",
            "QX/P",
            "HS/P"
        ));
        // On a reverse-strand model the allele is read as its reverse complement.
        let reverse = |entry: &Entry, allele: &str| {
            with_field(
                &with_field(entry, "STRAND", Some("-1")),
                "USED_REF",
                Some(allele),
            )
        };
        let (reverse_e, reverse_a) = (reverse(&e, "AACTTTACT"), reverse(&a, "AACTTTACT"));
        let (forward_e, forward_a) = (reverse(&e, "AGTAAAGTT"), reverse(&a, "AGTAAAGTT"));
        for name in ["Codons", "Amino_acids"] {
            assert!(intron_span_cell(&reverse_e, &reverse_a, name), "{name}");
            assert!(!intron_span_cell(&forward_e, &forward_a, name), "{name}");
        }
        // The terms and the notation read neither the codons nor the strand.
        for name in ["Consequence", "HGVSp"] {
            assert!(intron_span_cell(&forward_e, &forward_a, name), "{name}");
        }
    }

    #[test]
    fn whole_codons_are_translated_under_the_standard_code() {
        assert_eq!(standard_translation("cAGTAA").as_deref(), Some("Q*"));
        assert_eq!(standard_translation("cattca").as_deref(), Some("HS"));
        assert_eq!(standard_translation("").as_deref(), Some(""));
        for unread in ["cAGTA", "cAGTAN", "c-GTAA", "cAGTA\u{e9}"] {
            assert_eq!(standard_translation(unread), None, "{unread}");
        }
    }

    #[test]
    fn the_in_frame_test_reads_a_codons_value_by_its_two_readings() {
        // The reference codons begin or end with the alternate's, a `-` half holding no
        // base, or the two trim to whole codons of the reference and nothing of the
        // alternate.
        for accepted in [
            "GGCCCACGT/-",
            "AAAGGG/aaa",
            "cCACGT/cgt",
            "cCGCCGTCC/ccc",
            "aaaGGGccc/aaaccc",
        ] {
            assert!(inframe_deletion_reads_codons(accepted), "{accepted}");
        }
        // Trimming leaves a base of the alternate's, or bases of the reference's that are no
        // whole codons, and an alternate no shorter than the reference is no deletion.
        for refused in [
            "cCACCGCCG/ccc",
            "cCCCCA/cgt",
            "aCGt/at",
            "AAA/aaa",
            "AAA/aaaa",
            "-/-",
        ] {
            assert!(!inframe_deletion_reads_codons(refused), "{refused}");
        }
    }

    /// The two entries of a row of `refseq_edited_model_deleted_span_codons`, the reference's
    /// first: a deletion of twelve bases over the six a forward-strand model lacks, six
    /// coding bases inside the three codons of the coding bases 136 to 144.
    fn deleted_span_entries() -> (Entry, Entry) {
        let row = |terms: &str, codons: &str, amino_acids: &str| Entry {
            key: "12:122064784-122064795|-|NM_032790.4|Transcript".to_string(),
            location: Some("12:122064784-122064795".to_string()),
            allele: "-".to_string(),
            feature: "NM_032790.4".to_string(),
            consequence: terms.to_string(),
            ..entry(&[
                ("Consequence", terms),
                ("CDS_position", "137-142"),
                ("Amino_acids", amino_acids),
                ("Codons", codons),
                ("STRAND", "1"),
                ("USED_REF", "CACCGCCGCCGT"),
                ("BAM_EDIT", "OK"),
            ])
        };
        (
            row("protein_altering_variant", "cCACCGCCG/ccc", "PPP/P"),
            row("inframe_deletion", "cCGCCGTCC/ccc", "PPS/P"),
        )
    }

    /// Whether the class takes the field `name` of the two entries, each with its own value.
    fn deleted_span_cell(e: &Entry, a: &Entry, name: &str) -> bool {
        let value = |entry: &Entry| {
            entry
                .fields
                .get(name)
                .map_or("-", String::as_str)
                .to_string()
        };
        deleted_span_coding_reference_cell(e, a, name, &value(e), &value(a))
    }

    #[test]
    fn the_deleted_span_class_takes_the_codons_the_residues_and_the_term_pair_of_its_row() {
        let (e, a) = deleted_span_entries();
        for name in ["Consequence", "Codons", "Amino_acids"] {
            assert!(deleted_span_cell(&e, &a, name), "{name}");
            assert!(!deleted_span_cell(&a, &e, name), "{name}");
        }
        let noted = |entry: &Entry, value: &str| with_field(entry, "HGVSp", Some(value));
        assert!(!deleted_span_cell(
            &noted(&e, "NP_1.1:p.Pro47_Ser48delextTer?"),
            &noted(&a, "NP_1.1:p.Pro47_Ser48del"),
            "HGVSp"
        ));
        // The terms the two sides share stand beside the pair; any other difference of the
        // two sets, and the pair the other way round, is not the class's.
        assert!(deleted_span_cell(
            &with_terms(&e, "protein_altering_variant,splice_region_variant"),
            &with_terms(&a, "inframe_deletion,splice_region_variant"),
            "Consequence"
        ));
        for (reference, candidate) in [
            (
                "protein_altering_variant,splice_region_variant",
                "inframe_deletion",
            ),
            (
                "protein_altering_variant",
                "inframe_deletion,splice_region_variant",
            ),
            ("missense_variant", "inframe_deletion"),
            ("protein_altering_variant", "frameshift_variant"),
            ("inframe_deletion", "protein_altering_variant"),
        ] {
            assert!(
                !deleted_span_cell(
                    &with_terms(&e, reference),
                    &with_terms(&a, candidate),
                    "Consequence"
                ),
                "{reference} {candidate}"
            );
        }
        // Reference codons the in-frame test accepts explain no such term: a deletion of
        // whole codons leaves no alternate codon, which every reference begins with.
        let whole_codons = |entry: &Entry, codons: &str, amino_acids: &str| {
            let spanned = with_field(entry, "CDS_position", Some("3373-3381"));
            let allele = with_field(&spanned, "USED_REF", Some("GGCCCCCCAGGCCCACGT"));
            let residues = with_field(&allele, "Amino_acids", Some(amino_acids));
            with_field(&residues, "Codons", Some(codons))
        };
        let (whole_e, whole_a) = (
            whole_codons(&e, "GGCCCCCCA/-", "GPP/-"),
            whole_codons(&a, "GGCCCACGT/-", "GPR/-"),
        );
        assert!(!deleted_span_cell(&whole_e, &whole_a, "Consequence"));
        for name in ["Codons", "Amino_acids"] {
            assert!(deleted_span_cell(&whole_e, &whole_a, name), "{name}");
        }
        // Beside Codons that are no pair of the class it takes neither other cell.
        let agreeing = with_field(&a, "Codons", Some("cCACCGCCG/ccc"));
        for name in ["Consequence", "Codons", "Amino_acids"] {
            assert!(!deleted_span_cell(&e, &agreeing, name), "{name}");
        }
        // Each side's reference residues are the translation of its own reference codons,
        // and the alternate halves are equal.
        for (reference, candidate) in [("PPL/P", "PPS/P"), ("PPP/P", "PPT/P"), ("PPP/P", "PPS/S")] {
            assert!(
                !deleted_span_cell(
                    &with_field(&e, "Amino_acids", Some(reference)),
                    &with_field(&a, "Amino_acids", Some(candidate)),
                    "Amino_acids"
                ),
                "{reference} {candidate}"
            );
        }
    }
}
