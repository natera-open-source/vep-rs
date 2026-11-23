// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The LOFTEE items and filters.

use crate::consts::{
    CDNA_END_FIELD, INDEL_NOTATION_RE_PATTERN, INTEGER_RE_PATTERN, LOF_INFO_FIELDS, STRAND_FIELDS,
};
use crate::entry::Entry;
use crate::entry::Fields;
use crate::fields::{as_number, entry_kind, record_label};
use crate::pyfmt;
use crate::rules::downstream::field_of;
use crate::rules::{nearest_i64, whole_number, whole_text, SPLICE_ACCEPTOR_TERM};
use crate::shared::{blank, ABSENT};
use std::sync::LazyLock;

pub const LOF_FIELDS: &[&str] = &["LoF"];

pub const LOF_FILTER_FIELDS: &[&str] = &["LoF_filter"];

pub const LOF_FLAGS_FIELDS: &[&str] = &["LoF_flags"];

/// The `LoF_info` item that holds the GERP total of the walk to the stop codon.
pub const GERP_DIST_ITEM: &str = "GERP_DIST";

pub const NAGNAG_SITE_FLAG: &str = "NAGNAG_SITE";

/// `INDEL_NOTATION_RE` held to the whole text, as `fullmatch` holds it.
static INDEL_NOTATION_WHOLE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(INDEL_NOTATION_RE_PATTERN));

/// `INTEGER_RE` held to the whole text, as `fullmatch` holds it.
static INTEGER_WHOLE_RE: LazyLock<regex::Regex> = LazyLock::new(|| whole_text(INTEGER_RE_PATTERN));

/// The parts of a list value between its `,` and `&`.
fn list_parts(value: &str) -> impl Iterator<Item = &str> {
    value.split([',', '&'])
}

/// Whether the row's cDNA position names two positions, both numbers: the `a-b` of the
/// text formats' `cDNA_position`, read up to a `/total` suffix, or the JSON `cdna_start`
/// and `cdna_end`. A row that begins or ends outside an exon prints `?` or nothing there.
pub fn cdna_span_closed(e: &Entry) -> bool {
    fn before_total(text: &str) -> &str {
        pyfmt::partition(text, "/").0
    }
    let field = |name: &str| e.fields.get(name).map_or(ABSENT, String::as_str);
    let (start, end) = match e.fields.get("cDNA_position") {
        Some(position) => {
            let (start, _, end) = pyfmt::partition(before_total(position), "-");
            (start, end)
        }
        None => (
            before_total(field("cdna_start")),
            before_total(field(CDNA_END_FIELD)),
        ),
    };
    pyfmt::is_digit(start) && pyfmt::is_digit(end) && start != end
}

/// The number of bases an HGVSc inserts, duplicates or deletes: the `ins` bases, one for
/// a single position, else the span between two positions without an intron offset or a
/// `*` (the numbering has no position 0); `None` for any other notation.
///
/// A span past `i64` is the nearest `i64`.
pub fn notation_indel_length(hgvsc: &str) -> Option<i64> {
    let found = INDEL_NOTATION_WHOLE_RE.captures(pyfmt::partition(hgvsc, ":").2)?;
    let kind = &found[3];
    if let Some(inserted) = kind.strip_prefix("ins") {
        return Some(pyfmt::len(inserted) as i64);
    }
    let Some(end) = found.get(2) else {
        return Some(1);
    };
    let (start, end) = (&found[1], end.as_str());
    if !(INTEGER_WHOLE_RE.is_match(start) && INTEGER_WHOLE_RE.is_match(end)) {
        return None;
    }
    let (first, last) = (whole_number(start), whole_number(end));
    let mut length = last.saturating_sub(first).saturating_add(1);
    if first < 0 && 0 < last {
        length -= 1;
    }
    Some(nearest_i64(length))
}

/// The `KEY:VALUE` items of a `LoF_info` value in order, joined by `,` in the default,
/// tab and JSON formats and by `&` in the VCF format; `None` for an empty value.
pub fn lof_info_items(value: &str) -> Option<Vec<(String, String)>> {
    if blank(value) {
        return None;
    }
    Some(
        list_parts(value)
            .map(|part| {
                let (key, _, item) = pyfmt::partition(part, ":");
                (key.to_string(), item.to_string())
            })
            .collect(),
    )
}

/// The reference's and the candidate's `GERP_DIST` on a row whose two `LoF_info` values
/// carry the same items in the same order and differ in that item alone, two different
/// numbers, with the row's other items by key; `None` for any other row.
pub fn gerp_totals(e: &Entry, a: Option<&Entry>) -> Option<(f64, f64, Fields)> {
    let a = a?;
    let expected = lof_info_items(&field_of(Some(e), LOF_INFO_FIELDS))?;
    let actual = lof_info_items(&field_of(Some(a), LOF_INFO_FIELDS))?;
    if !expected
        .iter()
        .map(|(key, _)| key)
        .eq(actual.iter().map(|(key, _)| key))
    {
        return None;
    }
    let mut differing = expected
        .iter()
        .zip(&actual)
        .filter(|((_, value), (_, other))| value != other)
        .map(|((key, _), _)| key.as_str());
    if !(differing.next() == Some("GERP_DIST") && differing.next().is_none()) {
        return None;
    }
    let items: Fields = expected.into_iter().collect();
    let candidate_items: Fields = actual.into_iter().collect();
    let reference_total = as_number(&items["GERP_DIST"])?;
    let candidate_total = as_number(&candidate_items["GERP_DIST"])?;
    if reference_total == candidate_total {
        return None;
    }
    Some((reference_total, candidate_total, items))
}

/// Whether the reference's row carries none of the plugin's four fields and the candidate's
/// row carries a `LoF` beside a `LoF_info` that holds a `GERP_DIST` item, or beside no
/// `LoF_info`, as a run under `skip_lof_info` writes the row (LoF.pm:243 and :339): the row
/// of a GERP walk the reference's plugin died in (gerp_dist.pl:129) and the candidate
/// finished.
pub fn gerp_walk_the_candidate_alone_finished(e: &Entry, a: &Entry) -> bool {
    let carries = |row: &Entry, fields: &[&str]| !blank(&field_of(Some(row), fields));
    let reference_carries_none = [
        LOF_FIELDS,
        LOF_FILTER_FIELDS,
        LOF_FLAGS_FIELDS,
        LOF_INFO_FIELDS,
    ]
    .iter()
    .all(|fields| !carries(e, fields));
    let walked_or_unwritten = lof_info_items(&field_of(Some(a), LOF_INFO_FIELDS))
        .is_none_or(|items| items.iter().any(|(key, _)| key == GERP_DIST_ITEM));
    reference_carries_none && carries(a, LOF_FIELDS) && walked_or_unwritten
}

/// The items of a `LoF_filter` value, joined by `,` in the default, tab and JSON
/// formats and by `&` in the VCF format; none for an empty value.
pub fn lof_filter_items(value: &str) -> Vec<String> {
    if blank(value) {
        return Vec::new();
    }
    list_parts(value).map(str::to_string).collect()
}

/// The side, `reference` or `candidate`, whose `LoF_filter` carries `item` on a row
/// whose two `LoF_filter` values differ by that item alone and whose `LoF` reads `LC`
/// beside the item and, on the other side, `LC` beside another item and `HC` beside
/// none; `None` for any other row.
pub fn lof_filter_item_side(e: &Entry, a: Option<&Entry>, item: &str) -> Option<String> {
    let a = a?;
    let expected_items = lof_filter_items(&field_of(Some(e), LOF_FILTER_FIELDS));
    let actual_items = lof_filter_items(&field_of(Some(a), LOF_FILTER_FIELDS));
    let carried = |items: &[String]| items.iter().any(|held| held == item);
    if carried(&expected_items) == carried(&actual_items) {
        return None;
    }
    let expected_call = field_of(Some(e), LOF_FIELDS);
    let actual_call = field_of(Some(a), LOF_FIELDS);
    let (side, carrying_items, carrying_call, other_items, other_call) = if carried(&actual_items) {
        (
            "candidate",
            &actual_items,
            &actual_call,
            &expected_items,
            &expected_call,
        )
    } else {
        (
            "reference",
            &expected_items,
            &expected_call,
            &actual_items,
            &actual_call,
        )
    };
    let once = carrying_items.iter().filter(|held| *held == item).count() == 1;
    if !(once
        && carrying_items
            .iter()
            .filter(|held| *held != item)
            .eq(other_items))
    {
        return None;
    }
    let other_wanted = if other_items.is_empty() { "HC" } else { "LC" };
    if carrying_call != "LC" || other_call != other_wanted {
        return None;
    }
    Some(side.to_string())
}

/// Whether the two `LoF_flags` values of a row hold the same items in the same order but
/// for one `item` that one side alone carries, beside a `LoF`, a `LoF_filter` and a
/// `LoF_info` that are each equal on both sides. A `LoF_flags` value is read as
/// `lof_filter_items` reads a list, and an empty value equals an absent one.
pub fn lof_flags_differ_by(e: &Entry, a: &Entry, item: &str) -> bool {
    let expected_items = lof_filter_items(&field_of(Some(e), LOF_FLAGS_FIELDS));
    let actual_items = lof_filter_items(&field_of(Some(a), LOF_FLAGS_FIELDS));
    let carried = |items: &[String]| items.iter().filter(|held| *held == item).count();
    let (carrying_items, other_items) = match (carried(&expected_items), carried(&actual_items)) {
        (1, 0) => (&expected_items, &actual_items),
        (0, 1) => (&actual_items, &expected_items),
        _ => return false,
    };
    if !carrying_items
        .iter()
        .filter(|held| *held != item)
        .eq(other_items)
    {
        return false;
    }
    [LOF_FIELDS, LOF_FILTER_FIELDS, LOF_INFO_FIELDS]
        .iter()
        .all(|fields| {
            let expected = field_of(Some(e), fields);
            let actual = field_of(Some(a), fields);
            expected == actual || (blank(&expected) && blank(&actual))
        })
}

/// The strand a row states in its `STRAND` field, `1` or `-1`; `None` for a row that
/// states neither.
fn stated_strand(e: &Entry) -> Option<&str> {
    STRAND_FIELDS
        .iter()
        .find_map(|name| e.fields.get(*name))
        .map(String::as_str)
        .filter(|strand| matches!(*strand, "1" | "-1"))
}

/// Whether the terms of a row include `splice_acceptor_variant`.
fn acceptor_row(e: &Entry) -> bool {
    e.consequence
        .split(['&', ','])
        .any(|term| term == SPLICE_ACCEPTOR_TERM)
}

/// Whether the reference's entry `e` and the candidate's entry `a` of one row differ as
/// `loftee_nagnag_site_first_transcript_strand` reads it: the terms of both include
/// `splice_acceptor_variant`; their `LoF_flags` differ by `NAGNAG_SITE` alone beside an
/// equal `LoF`, `LoF_filter` and `LoF_info` (`lof_flags_differ_by`); and `record_rows`,
/// the reference's rows of the row's record, hold a `splice_acceptor_variant` row of the
/// row's kind on the other strand than its own. A row that states no strand has none; an
/// entry of `record_rows` under another record's label (`record_label`) is not read.
pub fn nagnag_site_first_transcript_strand_row(
    e: &Entry,
    a: &Entry,
    record_rows: &[&Entry],
) -> bool {
    if !(acceptor_row(e) && acceptor_row(a) && lof_flags_differ_by(e, a, NAGNAG_SITE_FLAG)) {
        return false;
    }
    let Some(own) = stated_strand(e) else {
        return false;
    };
    let (kind, label) = (entry_kind(e), record_label(e));
    record_rows.iter().any(|r| {
        entry_kind(r) == kind
            && record_label(r) == label
            && acceptor_row(r)
            && stated_strand(r).is_some_and(|strand| strand != own)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(pairs: &[(&str, &str)]) -> Entry {
        Entry {
            key: "21:100|A|ENST1|Transcript".to_string(),
            location: Some("21:100".to_string()),
            allele: "A".to_string(),
            feature: "ENST1".to_string(),
            consequence: "stop_gained".to_string(),
            fields: pairs
                .iter()
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect(),
            ..Entry::default()
        }
    }

    fn pairs(items: &[(&str, &str)]) -> Vec<(String, String)> {
        items
            .iter()
            .map(|(key, item)| (key.to_string(), item.to_string()))
            .collect()
    }

    #[test]
    fn a_cdna_position_is_closed_where_it_names_two_numbers() {
        let closed: Vec<bool> = [
            &[("cDNA_position", "1266-1297/2000")][..],
            &[("cdna_start", "1266"), ("cdna_end", "1297")],
            &[("cDNA_position", "?-1297")],
            &[("cDNA_position", "1266")],
            &[("cdna_end", "1297")],
            &[],
        ]
        .iter()
        .map(|fields| cdna_span_closed(&row(fields)))
        .collect();
        assert_eq!(closed, [true, true, false, false, false, false]);
    }

    #[test]
    fn a_notation_states_the_bases_it_inserts_duplicates_or_deletes() {
        let lengths: Vec<Option<i64>> = [
            "c.271_272dup",
            "c.5dup",
            "c.100_102del",
            "c.7_8insACG",
            "c.-2_3del",
            "c.10+1_10+2del",
            "c.9_*2del",
            "c.5A>G",
            "c.7_5del",
            "c.-5_-3del",
        ]
        .iter()
        .map(|notation| notation_indel_length(&format!("ENST1:{notation}")))
        .collect();
        assert_eq!(
            lengths,
            [
                Some(2),
                Some(1),
                Some(3),
                Some(3),
                Some(5),
                None,
                None,
                None,
                Some(-1),
                Some(3)
            ]
        );
        // The notation is the text after the sequence name, to its end.
        assert_eq!(notation_indel_length("c.5del"), None);
        assert_eq!(notation_indel_length("ENST1:c.5del\n"), None);
    }

    #[test]
    fn the_items_of_a_value_are_read_at_either_separator() {
        let two = pairs(&[("GERP_DIST", "-60.5"), ("BP_DIST", "79")]);
        assert_eq!(
            lof_info_items("GERP_DIST:-60.5,BP_DIST:79"),
            Some(two.clone())
        );
        assert_eq!(lof_info_items("GERP_DIST:-60.5&BP_DIST:79"), Some(two));
        assert_eq!(
            lof_info_items("A:b:c,D"),
            Some(pairs(&[("A", "b:c"), ("D", "")]))
        );
        assert_eq!(lof_info_items("-"), None);
        assert_eq!(lof_info_items(""), None);
        assert_eq!(
            lof_filter_items("END_TRUNC&ANC_ALLELE"),
            ["END_TRUNC", "ANC_ALLELE"]
        );
        assert_eq!(lof_filter_items("A:b:c,D"), ["A:b:c", "D"]);
        assert!(lof_filter_items("-").is_empty());
        assert!(lof_filter_items("").is_empty());
    }

    #[test]
    fn the_two_totals_are_read_where_they_alone_differ() {
        let whole = "GERP_DIST:-60.5,BP_DIST:79,DIST_FROM_LAST_EXON:40";
        let totals = |expected: &str, actual: &str| {
            gerp_totals(
                &row(&[("LoF_info", expected)]),
                Some(&row(&[("LoF_info", actual)])),
            )
            .map(|(reference, candidate, items)| (reference, candidate, items.len()))
        };
        assert_eq!(
            totals(whole, &whole.replace("-60.5", "-55.5")),
            Some((-60.5, -55.5, 3))
        );
        assert_eq!(totals(whole, whole), None);
        assert_eq!(
            totals(whole, &whole.replace("79", "80").replace("-60.5", "-1")),
            None
        );
        assert_eq!(totals(whole, &whole.replace("-60.5", "x")), None);
        // Two spellings of one number are no two totals.
        assert_eq!(totals(whole, &whole.replace("-60.5", "-60.50")), None);
        assert_eq!(
            totals(whole, "BP_DIST:79,GERP_DIST:-55.5,DIST_FROM_LAST_EXON:40"),
            None
        );
        // An item written twice is read at its last value.
        assert_eq!(
            totals(
                "GERP_DIST:-1,BP_DIST:79,GERP_DIST:-60.5",
                "GERP_DIST:-1,BP_DIST:79,GERP_DIST:-2"
            ),
            Some((-60.5, -2.0, 2))
        );
        assert!(gerp_totals(&row(&[("LoF_info", whole)]), None).is_none());
        assert!(gerp_totals(&row(&[]), Some(&row(&[("LoF_info", whole)]))).is_none());
        let lower = gerp_totals(
            &row(&[("lof_info", whole)]),
            Some(&row(&[("lof_info", &whole.replace("-60.5", "-58"))])),
        );
        assert_eq!(
            lower.map(|(reference, candidate, _)| (reference, candidate)),
            Some((-60.5, -58.0))
        );
    }

    #[test]
    fn a_filter_item_names_the_side_that_alone_carries_it() {
        let side = |expected: (&str, &str), actual: (&str, &str)| {
            lof_filter_item_side(
                &row(&[("LoF_filter", expected.0), ("LoF", expected.1)]),
                Some(&row(&[("LoF_filter", actual.0), ("LoF", actual.1)])),
                "ANC_ALLELE",
            )
        };
        let candidate = Some("candidate".to_string());
        assert_eq!(side(("", "HC"), ("ANC_ALLELE", "LC")), candidate);
        assert_eq!(
            side(("END_TRUNC", "LC"), ("END_TRUNC,ANC_ALLELE", "LC")),
            candidate
        );
        assert_eq!(
            side(("END_TRUNC", "LC"), ("ANC_ALLELE,END_TRUNC", "LC")),
            candidate
        );
        assert_eq!(
            side(("END_TRUNC", "LC"), ("ANC_ALLELE&END_TRUNC", "LC")),
            candidate
        );
        assert_eq!(
            side(("ANC_ALLELE", "LC"), ("-", "HC")),
            Some("reference".to_string())
        );
        // The other side reads HC beside no item and LC beside another.
        assert_eq!(side(("ANC_ALLELE", "LC"), ("-", "LC")), None);
        assert_eq!(side(("", "LC"), ("ANC_ALLELE", "LC")), None);
        assert_eq!(side(("", "HC"), ("ANC_ALLELE", "HC")), None);
        assert_eq!(side(("", "HC"), ("ANC_ALLELE,ANC_ALLELE", "LC")), None);
        assert_eq!(side(("END_TRUNC", "LC"), ("ANC_ALLELE", "LC")), None);
        assert_eq!(side(("ANC_ALLELE", "LC"), ("ANC_ALLELE", "LC")), None);
        assert_eq!(
            lof_filter_item_side(
                &row(&[("LoF_filter", "ANC_ALLELE"), ("LoF", "LC")]),
                None,
                "ANC_ALLELE"
            ),
            None
        );
        assert_eq!(
            lof_filter_item_side(
                &row(&[("LoF", "HC")]),
                Some(&row(&[("LoF_filter", "ANC_ALLELE"), ("LoF", "LC")])),
                "ANC_ALLELE"
            ),
            candidate
        );
    }

    #[test]
    fn a_walk_the_candidate_alone_finished_leaves_the_reference_without_the_four_fields() {
        let info = "GERP_DIST:705,BP_DIST:702,PERCENTILE:0.15,DIST_FROM_LAST_EXON:266";
        let finished = |expected: &[(&str, &str)], actual: &[(&str, &str)]| {
            gerp_walk_the_candidate_alone_finished(&row(expected), &row(actual))
        };
        let written = [("LoF", "HC"), ("LoF_info", info)];
        assert!(finished(&[], &written));
        assert!(finished(&[("IMPACT", "HIGH")], &written));
        // Each format's spelling of an empty value, and the JSON format's keys.
        let empty = [
            ("LoF", "-"),
            ("LoF_filter", ""),
            ("LoF_flags", "-"),
            ("LoF_info", ""),
        ];
        assert!(finished(&empty, &written));
        assert!(finished(&[], &[("lof", "LC"), ("lof_info", info)]));
        let joined_as_the_vcf_format_joins = info.replace(',', "&");
        assert!(finished(
            &[],
            &[("LoF", "HC"), ("LoF_info", &joined_as_the_vcf_format_joins)]
        ));
        // One field of the four on the reference's row is a row the plugin finished.
        for held in [
            ("LoF", "HC"),
            ("LoF_filter", "END_TRUNC"),
            ("LoF_flags", "SINGLE_EXON"),
            ("LoF_info", info),
            ("lof_flags", "SINGLE_EXON"),
        ] {
            assert!(!finished(&[held], &written), "{}", held.0);
        }
        // A run under `skip_lof_info` writes the call and its filter or flag without the items.
        assert!(finished(&[], &[("LoF", "HC")]));
        assert!(finished(&[], &[("LoF", "HC"), ("LoF_info", "-")]));
        assert!(finished(
            &empty,
            &[("LoF", "LC"), ("LoF_filter", "END_TRUNC"), ("LoF_info", "")]
        ));
        assert!(finished(
            &[],
            &[("lof", "HC"), ("lof_flags", "NO_EXON_NUMBER")]
        ));
        // A candidate without a call, or with items and no total among them.
        assert!(!finished(&[], &[("LoF_info", info)]));
        assert!(!finished(&[], &[("LoF", "-"), ("LoF_info", info)]));
        assert!(!finished(&[], &[("LoF_filter", "END_TRUNC")]));
        assert!(!finished(
            &[],
            &[("LoF", "HC"), ("LoF_info", "BP_DIST:702,PERCENTILE:0.15")]
        ));
        assert!(!finished(
            &[],
            &[("LoF", "HC"), ("LoF_info", "XGERP_DIST:705,BP_DIST:702")]
        ));
        assert!(!finished(&[], &[]));
    }

    #[test]
    fn two_flag_lists_differ_by_an_item_one_side_alone_carries_once() {
        let differs = |expected: &[(&str, &str)], actual: &[(&str, &str)]| {
            lof_flags_differ_by(&row(expected), &row(actual), NAGNAG_SITE_FLAG)
        };
        let flags = |value: &'static str| [("LoF_flags", value)];
        assert!(differs(
            &flags("NON_CAN_SPLICE,NAGNAG_SITE"),
            &flags("NON_CAN_SPLICE")
        ));
        assert!(differs(
            &flags("NON_CAN_SPLICE"),
            &flags("NAGNAG_SITE&NON_CAN_SPLICE")
        ));
        assert!(differs(&[], &flags("NAGNAG_SITE")));
        assert!(differs(&flags("NAGNAG_SITE"), &flags("-")));
        assert!(differs(&[("lof_flags", "NAGNAG_SITE")], &flags("")));
        // The item on both sides, on neither, twice on one, and other items that differ or
        // stand in another order.
        assert!(!differs(&flags("NAGNAG_SITE"), &flags("NAGNAG_SITE")));
        assert!(!differs(&flags("NON_CAN_SPLICE"), &[]));
        assert!(!differs(&flags("NAGNAG_SITE,NAGNAG_SITE"), &[]));
        assert!(!differs(
            &flags("NAGNAG_SITE,NON_CAN_SPLICE"),
            &flags("PHYLOCSF_WEAK")
        ));
        assert!(!differs(
            &flags("NON_CAN_SPLICE,NAGNAG_SITE,PHYLOCSF_WEAK"),
            &flags("PHYLOCSF_WEAK,NON_CAN_SPLICE")
        ));
        // The row's other three fields are equal, an empty value equal to an absent one.
        let carrying = [
            ("LoF_flags", "NAGNAG_SITE"),
            ("LoF", "HC"),
            ("LoF_filter", "-"),
            ("LoF_info", "INTRON_SIZE:221"),
        ];
        assert!(differs(
            &carrying,
            &[("LoF", "HC"), ("LoF_info", "INTRON_SIZE:221")]
        ));
        assert!(!differs(
            &carrying,
            &[("LoF", "LC"), ("LoF_info", "INTRON_SIZE:221")]
        ));
        assert!(!differs(
            &carrying,
            &[
                ("LoF", "HC"),
                ("LoF_filter", "5UTR_SPLICE"),
                ("LoF_info", "INTRON_SIZE:221")
            ]
        ));
        assert!(!differs(&carrying, &[("LoF", "HC")]));
    }

    /// A row of the record at `21:100` on `feature`, with its strand where it states one.
    fn stranded(feature: &str, strand: &str, consequence: &str, flags: &str) -> Entry {
        let mut entry = row(&[("LoF_flags", flags)]);
        entry.key = format!("21:100|A|{feature}|Transcript");
        entry.feature = feature.to_string();
        entry.consequence = consequence.to_string();
        if !strand.is_empty() {
            entry
                .fields
                .insert("STRAND".to_string(), strand.to_string());
        }
        entry
    }

    #[test]
    fn the_strand_class_reads_an_acceptor_row_on_the_other_strand_of_the_record() {
        let acceptor = "splice_acceptor_variant";
        let e = stranded("ENST2", "-1", acceptor, "NON_CAN_SPLICE,NAGNAG_SITE");
        let a = stranded("ENST2", "-1", acceptor, "NON_CAN_SPLICE");
        let taken = |e: &Entry, a: &Entry, beside: &Entry| {
            nagnag_site_first_transcript_strand_row(e, a, &[beside, e])
        };
        let forward = stranded("ENST1", "1", "intron_variant&splice_acceptor_variant", "");
        assert!(taken(&e, &a, &forward));
        assert!(taken(&a, &e, &forward));
        assert!(!nagnag_site_first_transcript_strand_row(&e, &a, &[&e]));
        assert!(!nagnag_site_first_transcript_strand_row(&e, &a, &[]));
        // The other row on the row's own strand, with no strand, with another term, of
        // another kind and of another record.
        assert!(!taken(&e, &a, &stranded("ENST1", "-1", acceptor, "")));
        assert!(!taken(&e, &a, &stranded("ENST1", "", acceptor, "")));
        assert!(!taken(&e, &a, &stranded("ENST1", "+", acceptor, "")));
        assert!(!taken(
            &e,
            &a,
            &stranded("ENST1", "1", "splice_donor_variant", "")
        ));
        let mut regulatory = forward.clone();
        regulatory.key = "21:100|A|ENSR1|RegulatoryFeature".to_string();
        assert!(!taken(&e, &a, &regulatory));
        let mut elsewhere = forward.clone();
        elsewhere.key = "21:200|A|ENST1|Transcript".to_string();
        assert!(!taken(&e, &a, &elsewhere));
        // The row's own terms and strand.
        let donor = |entry: &Entry| Entry {
            consequence: "splice_donor_variant".to_string(),
            ..entry.clone()
        };
        assert!(!taken(&donor(&e), &a, &forward));
        assert!(!taken(&e, &donor(&a), &forward));
        let unstated = stranded("ENST2", "", acceptor, "NON_CAN_SPLICE,NAGNAG_SITE");
        assert!(!taken(&unstated, &a, &forward));
        // The JSON format writes the strand under a lower-case key.
        let mut lower = forward.clone();
        let strand = lower.fields.shift_remove("STRAND").expect("a strand");
        lower.fields.insert("strand".to_string(), strand);
        assert!(taken(&e, &a, &lower));
    }
}
