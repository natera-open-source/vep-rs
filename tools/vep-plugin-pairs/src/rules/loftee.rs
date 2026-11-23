// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The LoFTEE classes: a `LoF_filter` item one engine alone writes with the row's `LoF`
//! call, a `GERP_DIST` total read through a summary level of the conservation file, and the
//! plugin's fields on a row whose GERP walk asks a conservation database the run lacks.

use super::entry::{as_number, blank, consequence_set, field_named, field_of, Entry, Row};
use super::shared::{first_alternate, record_indel_row, sibling_rows};
use super::Res;

pub const LOF_FIELDS: [&str; 1] = ["LoF"];
pub const LOF_FILTER_FIELDS: [&str; 1] = ["LoF_filter"];
pub const LOF_INFO_FIELDS: [&str; 1] = ["LoF_info"];
/// The fields the plugin writes on a row (LoF.pm:333-340).
pub const LOFTEE_FIELDS: [&str; 4] = ["LoF", "LoF_info", "LoF_filter", "LoF_flags"];
pub const GERP_DIST_ITEM: &str = "GERP_DIST";
pub const ANC_ALLELE_FILTER: &str = "ANC_ALLELE";
pub const LOF_SUBSTITUTION_TERMS: [&str; 3] = [
    "stop_gained",
    "splice_acceptor_variant",
    "splice_donor_variant",
];
pub const GC_TO_GT_DONOR_FILTER: &str = "GC_TO_GT_DONOR";
pub const LOF_DONOR_TERM: &str = "splice_donor_variant";
pub const LOF_TRUNCATING_TERMS: [&str; 2] = ["stop_gained", "frameshift_variant"];
pub const END_TRUNC_FILTER: &str = "END_TRUNC";
/// LoFTEE's default `gerp_end_trunc_cutoff`.
pub const GERP_END_TRUNC_CUTOFF: f64 = -58.0;
/// The interval, in bases, from which a total is read through a summary level.
pub const GERP_EXACT_BELOW_BASES: f64 = 80.0;
pub const LAST_EXON_REACH: f64 = 50.0;

/// The side whose `LoF_filter` carries an item.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    Reference,
    Candidate,
}

/// The `KEY:VALUE` items of a `LoF_info` value in order, joined by `,` in the default, tab
/// and JSON formats and by `&` in the VCF format; `None` for an empty value.
pub fn lof_info_items(value: &str) -> Option<Vec<(&str, &str)>> {
    if blank(value) {
        return None;
    }
    Some(
        value
            .split([',', '&'])
            .map(|part| part.split_once(':').unwrap_or((part, "")))
            .collect(),
    )
}

/// The value of the last of `items` under `key`, the one a mapping built from them holds.
fn item_of<'a>(items: &[(&'a str, &'a str)], key: &str) -> Option<&'a str> {
    items
        .iter()
        .rev()
        .find_map(|(name, value)| (*name == key).then_some(*value))
}

/// The two `GERP_DIST` totals of a row and the reference's `LoF_info` items.
#[derive(Debug, Clone, PartialEq)]
pub struct GerpTotals<'a> {
    pub reference: f64,
    pub candidate: f64,
    pub items: Vec<(&'a str, &'a str)>,
}

/// The reference's and the engine's `GERP_DIST` on a row whose two `LoF_info` values carry
/// the same items in the same order and differ in that item alone, two different numbers,
/// with the reference's items; `None` for any other row.
pub fn gerp_totals<'a>(e: &'a Entry, a: Option<&Entry>) -> Option<GerpTotals<'a>> {
    let items = lof_info_items(field_of(e, &LOF_INFO_FIELDS))?;
    let actual = lof_info_items(field_of(a?, &LOF_INFO_FIELDS))?;
    if !items
        .iter()
        .map(|(key, _)| key)
        .eq(actual.iter().map(|(key, _)| key))
    {
        return None;
    }
    let differing = items
        .iter()
        .zip(&actual)
        .filter(|((_, value), (_, other))| value != other)
        .map(|((key, _), _)| *key);
    if !differing.eq(["GERP_DIST"]) {
        return None;
    }
    let reference = as_number(item_of(&items, "GERP_DIST")?)?;
    let candidate = as_number(item_of(&actual, "GERP_DIST")?)?;
    (reference != candidate).then_some(GerpTotals {
        reference,
        candidate,
        items,
    })
}

/// The items of a `LoF_filter` value, joined by `,` in the default, tab and JSON formats and
/// by `&` in the VCF format; none for an empty value.
pub fn lof_filter_items(value: &str) -> Vec<&str> {
    if blank(value) {
        Vec::new()
    } else {
        value.split([',', '&']).collect()
    }
}

/// The side whose `LoF_filter` carries `item` on a row whose two `LoF_filter` values differ
/// by that item alone and whose `LoF` reads `LC` beside the item and, on the other side,
/// `LC` beside another item and `HC` beside none; `None` for any other row.
pub fn lof_filter_item_side(e: &Entry, a: Option<&Entry>, item: &str) -> Option<Side> {
    let a = a?;
    let expected_items = lof_filter_items(field_of(e, &LOF_FILTER_FIELDS));
    let actual_items = lof_filter_items(field_of(a, &LOF_FILTER_FIELDS));
    if expected_items.contains(&item) == actual_items.contains(&item) {
        return None;
    }
    let (expected_call, actual_call) = (field_of(e, &LOF_FIELDS), field_of(a, &LOF_FIELDS));
    let (side, carrying_items, carrying_call, other_items, other_call) =
        if actual_items.contains(&item) {
            (
                Side::Candidate,
                actual_items,
                actual_call,
                expected_items,
                expected_call,
            )
        } else {
            (
                Side::Reference,
                expected_items,
                expected_call,
                actual_items,
                actual_call,
            )
        };
    let carried = carrying_items.iter().filter(|i| **i == item).count();
    if carried != 1
        || !carrying_items
            .iter()
            .filter(|i| **i != item)
            .eq(other_items.iter())
    {
        return None;
    }
    let other_wanted = if other_items.is_empty() { "HC" } else { "LC" };
    (carrying_call == "LC" && other_call == other_wanted).then_some(side)
}

/// Whether one of `terms` is among the row's consequence terms.
fn carries_term(e: &Entry, terms: &[&str]) -> bool {
    consequence_set(&e.consequence)
        .split(',')
        .any(|term| terms.contains(&term))
}

/// LoF.pm:531 returns before the ancestral comparison when the record's allele string holds
/// `-` or is not three characters long, so no allele of a multi-allelic record is filtered
/// `ANC_ALLELE`. The shape: on a one-base substitution's stop-gained or splice-site row of a
/// record that has a row of another allele on the same feature, a `LoF_filter` the engine
/// writes as the reference's plus `ANC_ALLELE`, with the row's `LoF`. Nothing is accepted
/// without the record's rows.
pub fn anc_allele_multiallelic_record(
    e: &Entry,
    a: Option<&Entry>,
    record_rows: Option<&[Row]>,
) -> Res<bool> {
    let Some(record_rows) = record_rows else {
        return Ok(false);
    };
    if lof_filter_item_side(e, a, ANC_ALLELE_FILTER) != Some(Side::Candidate)
        || !carries_term(e, &LOF_SUBSTITUTION_TERMS)
    {
        return Ok(false);
    }
    if !matches!(e.allele.as_str(), "A" | "C" | "G" | "T") || record_indel_row(e, record_rows)? {
        return Ok(false);
    }
    Ok(!sibling_rows(e, record_rows).is_empty())
}

/// LoF.pm:363-367 test the reference and the first alternate allele of the record for every
/// allele's row, so a later allele of a multi-allelic donor record is filtered
/// `GC_TO_GT_DONOR`, or spared, by its sibling's base. The shape: on a
/// `splice_donor_variant` row of a record that has a row of another allele on the same
/// feature, and that is not the row of the record's first alternate, a `LoF_filter` whose
/// two values differ by `GC_TO_GT_DONOR` alone, either side carrying it. Where no entry
/// states the record's alleles in order (`first_alternate`) nothing is accepted.
pub fn first_alt_donor_filter(e: &Entry, a: Option<&Entry>, record_rows: Option<&[Row]>) -> bool {
    let Some(record_rows) = record_rows else {
        return false;
    };
    if lof_filter_item_side(e, a, GC_TO_GT_DONOR_FILTER).is_none()
        || !carries_term(e, &[LOF_DONOR_TERM])
        || sibling_rows(e, record_rows).is_empty()
    {
        return false;
    }
    first_alternate(e, record_rows).is_some_and(|first| first != e.allele)
}

/// gerp_dist.pl:70-79 total each exon interval's GERP score through one summary bin of the
/// bigWig, which the library serves from a zoom level once the interval is 80 bases, so
/// `GERP_DIST` approximates the per-base sum, and LoF.pm:254 pushes `END_TRUNC` by that
/// total. The shape: on a stop-gained or frameshift row whose `BP_DIST` is at least 79 (an
/// interval of N bases walks N - 1), a `LoF_info` whose two values differ in the `GERP_DIST`
/// item alone, both numbers; and, where the row's `DIST_FROM_LAST_EXON` is at most 50 and
/// the two totals straddle the cutoff, a `LoF_filter` or `LoF` differing by `END_TRUNC`
/// alone on the side of the lower total.
pub fn gerp_zoom_level_summary(e: &Entry, a: Option<&Entry>, name: &str) -> bool {
    let Some(totals) = gerp_totals(e, a) else {
        return false;
    };
    if !carries_term(e, &LOF_TRUNCATING_TERMS) {
        return false;
    }
    let number = |key: &str| as_number(item_of(&totals.items, key).unwrap_or(""));
    if number("BP_DIST").is_none_or(|walked| walked + 1.0 < GERP_EXACT_BELOW_BASES) {
        return false;
    }
    if field_named(name, &LOF_INFO_FIELDS) {
        return true;
    }
    let Some(side) = lof_filter_item_side(e, a, END_TRUNC_FILTER) else {
        return false;
    };
    if number("DIST_FROM_LAST_EXON").is_none_or(|reach| reach > LAST_EXON_REACH) {
        return false;
    }
    let (lower, higher) = match side {
        Side::Reference => (totals.reference, totals.candidate),
        Side::Candidate => (totals.candidate, totals.reference),
    };
    lower <= GERP_END_TRUNC_CUTOFF && GERP_END_TRUNC_CUTOFF < higher
}

/// On the plugin's `master` branch gerp_dist.pl:59-60 hand an exon that holds neither the
/// variant nor the stop codon to `get_exon_gerp`, which reads the exon's total from the
/// conservation database (gerp_dist.pl:129). LoF.pm:119-120 leave the word `false` in that
/// database's place in a run without `conservation_file`, so the walk that LoF.pm:258
/// starts dies at that read and Ensembl VEP writes the row without the plugin's fields.
/// The shape: on a row whose terms hold `stop_gained` or `frameshift_variant`, the
/// engine's as well where its entry states terms, each of `LOFTEE_FIELDS` the engine alone
/// writes, where the reference's row carries none of the four and the engine's carries
/// `LoF` beside a `LoF_info` holding a `GERP_DIST` item, or beside no `LoF_info`, as a run
/// under `skip_lof_info` writes the row (LoF.pm:243 and LoF.pm:339).
pub fn gerp_walk_without_conservation_database(
    e: &Entry,
    a: Option<&Entry>,
    ev: &str,
    av: &str,
) -> bool {
    let Some(a) = a else {
        return false;
    };
    if !blank(ev) || blank(av) {
        return false;
    }
    let truncating = |row: &Entry| carries_term(row, &LOF_TRUNCATING_TERMS);
    if !truncating(e) || !(a.consequence.is_empty() || truncating(a)) {
        return false;
    }
    if LOFTEE_FIELDS
        .iter()
        .any(|field| !blank(field_of(e, &[field])))
        || blank(field_of(a, &LOF_FIELDS))
    {
        return false;
    }
    lof_info_items(field_of(a, &LOF_INFO_FIELDS))
        .is_none_or(|items| item_of(&items, GERP_DIST_ITEM).is_some())
}
