// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The NearestExonJB classes: a value the plugin gives or withholds by the records before a
//! row, a value read at the 3'-shifted span, one item folded from exons without a stable id,
//! and items printed in the order of a hash's keys.

use super::entry::{blank, consequence_set, entry_kind, feature_kind, field_of, Entry};
use super::shared::{
    cdna_span_closed, list_parts, positive_shift, re_digit, DISTANCE_FIELDS, SPLICE_SITE_TERMS,
};
use super::Res;
use crate::pytext;
use regex::Regex;
use std::collections::BTreeSet;
use std::sync::LazyLock;

/// The bases outside an exon within which BaseTranscriptVariation.pm:858-865 count a record
/// as overlapping it once a frameshift intron was seen.
pub const FRAMESHIFT_INTRON_STRETCH: i128 = 12;

/// One item of a `NearestExonJB` value, the exon id, distance, side (`start_end` on a tie)
/// and length: `|` joins the parts and `,` the items in the default, tab and JSON formats,
/// `+` and `&` in the VCF format.
fn item_pattern() -> String {
    let digit = re_digit();
    format!(r"([^|+&,]*)[|+]([{digit}]+)[|+](start|end|start_end)[|+]([{digit}]+)")
}

static NEAREST_EXON_ITEM: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(&item_pattern()).expect("regex"));

static NEAREST_EXON_LIST: LazyLock<Regex> = LazyLock::new(|| {
    let item = item_pattern();
    Regex::new(&format!(r"^{item}(?:[,&]{item})*$")).expect("regex")
});

/// One item of a `NearestExonJB` value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NearestExonItem<'a> {
    pub exon_id: &'a str,
    pub distance: &'a str,
    pub side: &'a str,
    pub length: &'a str,
}

/// The items of a `NearestExonJB` value in any format's spelling; `None` for a value that
/// is not a list of such items.
pub fn nearest_exon_items(value: &str) -> Option<Vec<NearestExonItem<'_>>> {
    if !NEAREST_EXON_LIST.is_match(value) {
        return None;
    }
    let part = |found: &regex::Captures<'_>, group: usize| {
        found.get(group).map_or("", |m| &value[m.range()])
    };
    Some(
        NEAREST_EXON_ITEM
            .captures_iter(value)
            .map(|found| NearestExonItem {
                exon_id: part(&found, 1),
                distance: part(&found, 2),
                side: part(&found, 3),
                length: part(&found, 4),
            })
            .collect(),
    )
}

/// NearestExonJB.pm:184 reads the overlapped exons, whose 12-base stretch applies only once
/// a record inside the transcript has set its frameshift-intron flag, so a record within 12
/// bases outside the first or last exon gets a value or none by the records before it. The
/// shape: a value the engine alone fills on a row whose `DISTANCE` is at most 12.
pub fn stretch_priming(e: &Entry, ev: &str, av: &str) -> bool {
    let distance = field_of(e, &DISTANCE_FIELDS);
    blank(ev)
        && !blank(av)
        && pytext::isdigit(distance)
        && pytext::int(distance).is_some_and(|distance| distance <= FRAMESHIFT_INTRON_STRETCH)
}

/// Under `--shift_3prime 1` the consequence tests fill the overlapped-exons memo at the
/// shifted span and NearestExonJB.pm:184 reads that memo while :225-226 measure from the
/// record's span. The shape: on a transcript row with a positive `SHIFT_LENGTH`, a value one
/// side fills and the other leaves empty, or whose two values name different exons.
pub fn shifted_span(e: &Entry, ev: &str, av: &str) -> Res<bool> {
    if !positive_shift(e)? || feature_kind(entry_kind(e)) != "transcript" {
        return Ok(false);
    }
    if blank(ev) || blank(av) {
        return Ok(blank(ev) != blank(av));
    }
    let (Some(expected_items), Some(actual_items)) =
        (nearest_exon_items(ev), nearest_exon_items(av))
    else {
        return Ok(false);
    };
    let exons = |items: &[NearestExonItem]| -> BTreeSet<String> {
        items.iter().map(|item| item.exon_id.to_string()).collect()
    };
    Ok(exons(&expected_items) != exons(&actual_items))
}

/// NearestExonJB.pm:220-237 file each exon's distance, side and length under its stable id,
/// so the exons of a row that carry none share the empty key and one item is printed from
/// their parts. The shape: both values hold items of an empty exon id, the reference's one
/// item whose distance and side are those of one of the engine's items, on a row the engine
/// gives two or more items, or one whose Consequence carries both splice-site terms and
/// whose cDNA position names two positions.
pub fn id_less_exon_fold(e: &Entry, ev: &str, av: &str) -> bool {
    let (Some(expected_items), Some(actual_items)) =
        (nearest_exon_items(ev), nearest_exon_items(av))
    else {
        return false;
    };
    let [folded] = expected_items[..] else {
        return false;
    };
    if !folded.exon_id.is_empty() || actual_items.iter().any(|item| !item.exon_id.is_empty()) {
        return false;
    }
    if !actual_items
        .iter()
        .any(|item| (item.distance, item.side) == (folded.distance, folded.side))
    {
        return false;
    }
    let terms = consequence_set(&e.consequence);
    let terms: Vec<&str> = terms.split(',').collect();
    actual_items.len() > 1
        || (SPLICE_SITE_TERMS.iter().all(|term| terms.contains(term)) && cdna_span_closed(e))
}

/// Whether two list-valued field texts carry the same items the same number of times in
/// different orders, so two or more of them: NearestExonJB.pm:267-271 build the value by
/// walking the keys of a hash.
pub fn reordered_items(expected: &str, actual: &str) -> bool {
    let (mut expected_items, mut actual_items) = list_parts(expected, actual);
    if expected_items == actual_items {
        return false;
    }
    expected_items.sort_unstable();
    actual_items.sort_unstable();
    expected_items == actual_items
}
