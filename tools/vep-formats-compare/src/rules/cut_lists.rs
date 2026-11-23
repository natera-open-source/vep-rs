// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The co-located lists cut at the cache cap.

use crate::entry::Entry;
use crate::entry::Fields;
use crate::fields::extra_pairs;
use crate::lists::field_objects;
use crate::parse::json::json_compact;
use crate::parse::json::json_scalar;
use crate::parse::json::scalar_frames;
use crate::pyfmt;
use crate::pyjson;
use crate::pyjson::Value;
use crate::shared::blank;
use indexmap::IndexMap;
use std::borrow::Cow;

pub const CACHE_LIST_CAP: i64 = 1024;

/// The separators of each list column's stored value, the first of them the one that joins the
/// values of a row's known variants.
pub const CUT_LIST_SEPARATORS: &[(&str, &[&[u8]])] =
    &[("PUBMED", &[b","]), ("VAR_SYNONYMS", &[b"--", b"::", b","])];

pub const COLOCATED_VARIANTS_FIELD: &str = "colocated_variants";

/// Whether a member of the reference's entry and the candidate's have a list's cut shape.
pub type CutShape = fn(&Value, &Value) -> bool;

/// The members of a JSON `colocated_variants` entry that hold a stored list, each with the
/// predicate of its cut shape.
pub const JSON_CUT_LISTS: &[(&str, CutShape)] = &[
    ("pubmed", crate::rules::cut_lists::json_cut_pubmed),
    ("var_synonyms", crate::rules::cut_lists::json_cut_synonyms),
];

/// `CACHE_LIST_CAP` as a count of bytes.
const CAP: usize = CACHE_LIST_CAP as usize;

/// The frame of `cut_list_reading` on the interpreter's call stack where `compare_entries`
/// reads a differing field: above the module, `main`, `run`, `compare_files` and
/// `compare_entries`.
const READING_FRAME: usize = 6;

/// The frame of `json_cut_lists_pair` under `cut_list_reading`, and of the predicate of a
/// cut shape under it.
const JSON_PAIR_FRAME: usize = READING_FRAME + 1;
const CUT_SHAPE_FRAME: usize = JSON_PAIR_FRAME + 1;

/// The frame of the first `row_cut_lists_pair` under `cut_list_reading` and its `cut_pair`.
const ROW_PAIR_FRAME: usize = READING_FRAME + 2;

/// What the interpreter says where a call needs a frame past its limit.
const RECURSION: &str = "maximum recursion depth exceeded";

fn separators_of(column: &str) -> Option<&'static [&'static [u8]]> {
    CUT_LIST_SEPARATORS
        .iter()
        .find(|(name, _)| *name == column)
        .map(|(_, separators)| *separators)
}

/// `CUT_LIST_SEPARATORS[column]`, which has no entry for any other column.
fn column_separators(column: &str) -> &'static [&'static [u8]] {
    match separators_of(column) {
        Some(separators) => separators,
        None => crate::uncaught!("{}", pyfmt::repr_str(column)),
    }
}

fn find(text: &[u8], pattern: &[u8]) -> Option<usize> {
    memchr::memmem::find(text, pattern)
}

fn rfind(text: &[u8], pattern: &[u8]) -> Option<usize> {
    memchr::memmem::rfind(text, pattern)
}

/// Whether `stored` is a value `cut_list_whole_items` reads: one of `CACHE_LIST_CAP`
/// bytes that ends in no separator.
fn cut_inside_an_item(stored: &[u8], separators: &[&[u8]]) -> bool {
    stored.len() == CAP
        && !separators
            .iter()
            .any(|separator| stored.ends_with(separator))
}

/// The whole items of a value `cut_inside_an_item` passes.
fn whole_items<'a>(stored: &'a [u8], column: &str, separators: &[&[u8]]) -> &'a [u8] {
    let last = separators
        .iter()
        .filter_map(|separator| rfind(stored, separator))
        .max()
        .unwrap_or(0);
    let mut whole = &stored[..last];
    if column == "VAR_SYNONYMS" {
        if stored[whole.len()..].starts_with(b"::") {
            whole = &whole[..rfind(whole, b"--").unwrap_or(0)];
        }
        let lists_an_id = split(whole, b"--").any(|group| {
            find(group, b"::").is_some_and(|at| group[at + 2..].iter().any(|byte| *byte != b','))
        });
        if !lists_an_id {
            return b"";
        }
    }
    whole
}

/// `text.split(separator)` of bytes.
fn split<'a>(text: &'a [u8], separator: &'a [u8]) -> impl Iterator<Item = &'a [u8]> {
    let mut rest = Some(text);
    std::iter::from_fn(move || {
        let text = rest?;
        match find(text, separator) {
            Some(at) => {
                rest = Some(&text[at + separator.len()..]);
                Some(&text[..at])
            }
            None => {
                rest = None;
                Some(text)
            }
        }
    })
}

/// The whole items of a stored `column` value the cache cut inside its last item: the
/// bytes before its last separator, or before the `--` of the last group when that
/// separator is the group's `::` (a `var_synonyms` source whose only id is cut goes with
/// it), none of them when they hold no id (a `pubmed` value without a `,`, a
/// `var_synonyms` value without an id after a `::`). `None` for a value not cut so:
/// one of another length than `CACHE_LIST_CAP` bytes, or one ending in a separator.
pub fn cut_list_whole_items(stored: &[u8], column: &str) -> Option<Vec<u8>> {
    let separators = column_separators(column);
    cut_inside_an_item(stored, separators).then(|| whole_items(stored, column, separators).to_vec())
}

/// A default, tab or VCF list value as the bytes of the stored values it joins: none for
/// an empty or absent value, `;` for the `%3B` the default and VCF formats write
/// (VEP_output.pm:118, VCF.pm:402) and, in the VCF format, `,` for `&` (VCF.pm:401).
pub fn row_list_bytes(value: &str, vcf: bool) -> Vec<u8> {
    if blank(value) {
        return Vec::new();
    }
    let text = value.replace("%3B", ";");
    if vcf {
        text.replace('&', ",").into_bytes()
    } else {
        text.into_bytes()
    }
}

/// Whether `candidate` is `reference`, a row's `column` value (`row_list_bytes`),
/// without the cut last item of one or more of the stored values it joins
/// (`colocated_list_cut_at_cache_cap`). A run of `CACHE_LIST_CAP` bytes that begins at the
/// value's start or after a joining separator, at `start` or later, and ends at the
/// value's end or before a joining separator is read as `cut_list_whole_items` reads it,
/// a run left with no id going with one joining separator; the bytes up to the run's last
/// whole item are the candidate's, and the rest is the candidate's as it stands or with a
/// later run read the same way.
pub fn row_cut_lists_pair(reference: &[u8], candidate: &[u8], column: &str, start: i64) -> bool {
    let start = usize::try_from(start).unwrap_or(0);
    cut_lists_pair_at(reference, candidate, column, start, ROW_PAIR_FRAME)
}

/// `row_cut_lists_pair` called at the frame `frame` of the interpreter's call stack. The
/// rest of a value behind each run read is read a frame deeper a run, and a run's last
/// separator is read two frames above the frame that reads the run: the comparison ends
/// where a run is read with fewer frames left than that.
fn cut_lists_pair_at(
    reference: &[u8],
    candidate: &[u8],
    column: &str,
    start: usize,
    frame: usize,
) -> bool {
    let separators = column_separators(column);
    let join = separators[0];
    let size = reference.len();
    // A run begins at the value's start and after every joining separator, one that
    // overlaps the separator before it included; no run before `start` is read.
    let from = start.saturating_sub(join.len()).min(size);
    let begins: Vec<usize> = std::iter::once(0)
        .chain(
            memchr::memchr_iter(join[0], &reference[from..])
                .map(|at| from + at)
                .filter(|at| reference[*at..].starts_with(join))
                .map(|at| at + join.len()),
        )
        .collect();
    for begin in begins {
        let end = begin + CAP;
        if begin < start || !(end == size || (end < size && reference[end..].starts_with(join))) {
            continue;
        }
        let run = &reference[begin..end];
        if !cut_inside_an_item(run, separators) {
            continue;
        }
        if frame + 2 > pyjson::RECURSION_LIMIT {
            crate::uncaught!("{RECURSION}");
        }
        // What is read of the value is a prefix of it: the bytes before the run and the
        // run's whole items, or the bytes before the separator that joins an emptied run.
        let whole = whole_items(run, column, separators).len();
        let (read, rest) = if whole > 0 {
            (begin + whole, end)
        } else if end < size {
            (begin, end + join.len())
        } else {
            (begin.saturating_sub(join.len()), size)
        };
        let (read, rest) = (&reference[..read], &reference[rest..]);
        let Some(left) = candidate.strip_prefix(read) else {
            continue;
        };
        if left == rest {
            return true;
        }
        let later = read.len() + if whole > 0 { join.len() } else { 0 };
        if cut_lists_pair_at(&[read, rest].concat(), candidate, column, later, frame + 1) {
            return true;
        }
    }
    false
}

/// `json_compact(value)` called by a function at the frame `frame`: the run ends on a
/// value nested deeper than the frames the interpreter has left.
fn compact_at(value: &Value, frame: usize) -> String {
    if frame + pyjson::compact_frames(value) > pyjson::RECURSION_LIMIT {
        crate::uncaught!("{RECURSION}");
    }
    json_compact(value)
}

/// `json_scalar(value)` called by a function at the frame `frame`.
fn scalar_at(value: &Value, frame: usize) -> String {
    if frame + scalar_frames(value) > pyjson::RECURSION_LIMIT {
        crate::uncaught!("{RECURSION}");
    }
    json_scalar(value)
}

/// Whether what a cut leaves of the reference's list is the candidate's `kept`: no member
/// where nothing is left, else the same compact text.
fn left_is_kept(left: Option<Value>, kept: &Value) -> bool {
    let left = left.unwrap_or(Value::Null);
    if matches!(left, Value::Null) != matches!(kept, Value::Null) {
        return false;
    }
    let left = compact_at(&left, CUT_SHAPE_FRAME);
    left == compact_at(kept, CUT_SHAPE_FRAME)
}

/// Whether a JSON `pubmed` array `ids` of the reference and the candidate's `kept`
/// (`None` for no member) have the shape of `colocated_list_cut_at_cache_cap`: the candidate
/// lists the reference's ids but the last, and has no member when none is left; the
/// reference's ids joined with `,`, the stored value, are cut inside the last
/// (`cut_list_whole_items`).
pub fn json_cut_pubmed(ids: &Value, kept: &Value) -> bool {
    let Some((_, left)) = ids.as_list().and_then(|ids| ids.split_last()) else {
        return false;
    };
    // A list beside no member of the candidate's may be nested deeper than a copy walks:
    // it is refused before anything of it is copied.
    if left.is_empty() != matches!(kept, Value::Null) {
        return false;
    }
    let left = (!left.is_empty()).then(|| Value::List(left.to_vec()));
    if !left_is_kept(left, kept) {
        return false;
    }
    cut_list_whole_items(scalar_at(ids, CUT_SHAPE_FRAME).as_bytes(), "PUBMED").is_some()
}

/// Whether a JSON `var_synonyms` object `groups` of the reference and the
/// candidate's `kept` (`None` for no member) have the shape of
/// `colocated_list_cut_at_cache_cap`. One source of the reference is the stored value's last
/// group: the candidate lists its ids but the last, or lacks the source when the reference
/// lists one id under it or none (a name with no `::`), and has no member when no id is
/// left under any source. Every other source lists an id, so the groups, `source::id,id`,
/// joined with `--` are the stored value, cut inside its last item
/// (`cut_list_whole_items`). The object's keys are not in the stored order, so any source
/// the candidate's object fits this way is taken as the last.
///
/// The reference prints a piece of digits that closes on a point as the number of those
/// digits (Utils.pm:241 adds 0 to each id its number test admits), a byte shorter than the
/// stored piece. Where the last id of the source taken as the last is a `Value::Int`, a
/// number with no decimal part within a signed 64 bits, the groups with a point after them are read
/// as the stored value too. The object does not tell that value from one a byte shorter
/// that ends in the digits, or from one of the same length that ends in a `,` after them.
pub fn json_cut_synonyms(groups: &Value, kept: &Value) -> bool {
    let Some(groups) = groups.as_dict() else {
        return false;
    };
    if !groups.values().all(|listed| listed.as_list().is_some()) {
        return false;
    }
    let ids_of = |listed: &Value| -> usize { listed.as_list().map_or(0, <[Value]>::len) };
    for (last, ids) in groups {
        let others = || groups.iter().filter(|(source, _)| *source != last);
        if !others().all(|(_, listed)| ids_of(listed) > 0) {
            continue;
        }
        // As in `json_cut_pubmed`: nothing of an object beside no member is copied.
        let nothing_left = others().next().is_none() && ids_of(ids) < 2;
        if nothing_left != matches!(kept, Value::Null) {
            continue;
        }
        let mut left: IndexMap<String, Value> = others()
            .map(|(source, listed)| (source.clone(), listed.clone()))
            .collect();
        if let Some((_, before)) = ids.as_list().and_then(|ids| ids.split_last()) {
            if !before.is_empty() {
                left.insert(last.clone(), Value::List(before.to_vec()));
            }
        }
        // Every other source lists an id, so an object left with a member lists one.
        let left = (!left.is_empty()).then_some(Value::Dict(left));
        if !left_is_kept(left, kept) {
            continue;
        }
        let mut texts: Vec<String> = others()
            .map(|(name, listed)| format!("{name}::{}", scalar_at(listed, CUT_SHAPE_FRAME)))
            .collect();
        texts.push(if ids_of(ids) > 0 {
            format!("{last}::{}", scalar_at(ids, CUT_SHAPE_FRAME))
        } else {
            last.clone()
        });
        let mut stored = texts.join("--");
        if cut_list_whole_items(stored.as_bytes(), "VAR_SYNONYMS").is_some() {
            return true;
        }
        let cut_piece = ids.as_list().and_then(|ids| ids.last());
        if matches!(cut_piece, Some(Value::Int(_))) {
            stored.push('.');
            if cut_list_whole_items(stored.as_bytes(), "VAR_SYNONYMS").is_some() {
                return true;
            }
        }
    }
    false
}

/// Whether two JSON `colocated_variants` values have the shape of
/// `colocated_list_cut_at_cache_cap`: the same known variants in the same order, every member
/// of each equal but the `pubmed` (`json_cut_pubmed`) or the `var_synonyms`
/// (`json_cut_synonyms`) of at least one.
///
/// The members of a variant are walked in the reference's order and then the candidate's
/// own. The answer is the same in any order but where one member of a variant is refused
/// and another is nested past the interpreter's frame limit: the member walked first
/// decides whether the pair is refused or the run ends.
pub fn json_cut_lists_pair(expected: &str, actual: &str) -> bool {
    let expected_variants = field_objects(expected);
    let actual_variants = field_objects(actual);
    if expected_variants.len() != actual_variants.len() {
        return false;
    }
    let mut cut = false;
    for (reference, candidate) in expected_variants.iter().zip(&actual_variants) {
        let (Some(reference), Some(candidate)) = (reference.as_dict(), candidate.as_dict()) else {
            continue;
        };
        let keys = reference
            .keys()
            .chain(candidate.keys().filter(|key| !reference.contains_key(*key)));
        for key in keys {
            let (ours, theirs) = (reference.get(key), candidate.get(key));
            if let (Some(ours), Some(theirs)) = (ours, theirs) {
                let ours = compact_at(ours, JSON_PAIR_FRAME);
                if ours == compact_at(theirs, JSON_PAIR_FRAME) {
                    continue;
                }
            }
            let cut_shape = JSON_CUT_LISTS
                .iter()
                .find(|(member, _)| *member == key.as_str())
                .map(|(_, cut_shape)| cut_shape);
            let Some(cut_shape) = cut_shape else {
                return false;
            };
            if !cut_shape(ours.unwrap_or(&Value::Null), theirs.unwrap_or(&Value::Null)) {
                return false;
            }
            cut = true;
        }
    }
    cut
}

/// The reference's differing field `name` of the entry `e` as the candidate reads the
/// cut lists in it (`colocated_list_cut_at_cache_cap`), `None` for a field holding none:
/// the candidate's value `av` for a `PUBMED` or `VAR_SYNONYMS` value
/// (`row_cut_lists_pair`; an entry without a Location is a VCF format's) and for a JSON
/// `colocated_variants` value (`json_cut_lists_pair`) of the class's shape; for the
/// default format's `Extra` text, the reference's `ev` with each such pair carrying the
/// candidate's value, or left out where the candidate writes none, the text's other pairs
/// as they stand.
pub fn cut_list_reading(e: &Entry, name: &str, ev: &str, av: &str) -> Option<String> {
    if name == COLOCATED_VARIANTS_FIELD {
        return json_cut_lists_pair(ev, av).then(|| av.to_string());
    }
    let vcf = e.location.is_none();
    let cut_pair = |column: &str, reference: &str, candidate: &str| -> bool {
        separators_of(column).is_some()
            && row_cut_lists_pair(
                &row_list_bytes(reference, vcf),
                &row_list_bytes(candidate, vcf),
                column,
                0,
            )
    };
    if name != "Extra" {
        return cut_pair(name, ev, av).then(|| av.to_string());
    }
    // The candidate's pairs are read once a pair of the reference names a list column.
    let mut written: Option<Fields> = None;
    let mut pairs: Vec<Cow<'_, str>> = Vec::new();
    let mut read = false;
    for pair in ev.split(';') {
        let (key, _, value) = pyfmt::partition(pair, "=");
        if separators_of(key).is_some() {
            let candidate = written.get_or_insert_with(|| extra_pairs(av)).get(key);
            if cut_pair(key, value, candidate.map_or("-", String::as_str)) {
                read = true;
                if let Some(candidate) = candidate {
                    pairs.push(Cow::Owned(format!("{key}={candidate}")));
                }
                continue;
            }
        }
        pairs.push(Cow::Borrowed(pair));
    }
    if !read {
        return None;
    }
    let joined = pairs.join(";");
    Some(if joined.is_empty() {
        "-".to_string()
    } else {
        joined
    })
}

/// The `mask_examples` line of a field `colocated_list_cut_at_cache_cap` accepts: the two
/// values' lengths, since a value at the cap runs to a thousand characters of ids.
pub fn cut_list_example(e: &Entry, name: &str, ev: &str, av: &str) -> String {
    let written = if blank(av) { 0 } else { pyfmt::len(av) };
    format!(
        "{}: {name} expected {} characters, got {written}",
        e.key,
        pyfmt::len(ev)
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Uncaught;

    /// What a run says that ends on an input no reading of the comparison expects.
    fn uncaught<T>(run: impl FnOnce() -> T + std::panic::UnwindSafe) -> String {
        let payload = std::panic::catch_unwind(run).err().expect("the run ends");
        let Uncaught(what) = payload.downcast_ref::<Uncaught>().expect("an exception");
        what.clone()
    }

    /// A value of `CACHE_LIST_CAP` bytes that ends in `tail`.
    fn capped(tail: &str) -> Vec<u8> {
        let mut stored = "9".repeat(CAP - tail.len()).into_bytes();
        stored.extend_from_slice(tail.as_bytes());
        stored
    }

    fn row() -> Entry {
        Entry {
            key: "21:500|G|ENST1|Transcript".to_string(),
            location: Some("21:500".to_string()),
            ..Entry::default()
        }
    }

    #[test]
    fn the_whole_items_of_a_stored_value_end_before_its_cut_item() {
        let whole = |tail: &str, column: &str| {
            cut_list_whole_items(&capped(tail), column)
                .map(|whole| String::from_utf8(whole).unwrap())
        };
        let kept = |tail: &str, kept: &str| Some(format!("{}{kept}", "9".repeat(CAP - tail.len())));
        assert_eq!(whole(",77,7", "PUBMED"), kept(",77,7", ",77"));
        assert_eq!(whole("", "PUBMED"), Some(String::new()));
        assert_eq!(whole(",7,", "PUBMED"), None);
        assert_eq!(cut_list_whole_items(&capped("")[1..], "PUBMED"), None);
        assert_eq!(whole("--A::1,2", "PUBMED"), kept("--A::1,2", "--A::1"));
        // A value whose items hold no id after a `::` keeps none of them.
        assert_eq!(whole("--A::,2", "VAR_SYNONYMS"), Some(String::new()));
        let stored = |tail: &str| format!("S::{}{tail}", "9".repeat(CAP - 3 - tail.len()));
        let synonyms = |tail: &str| {
            cut_list_whole_items(stored(tail).as_bytes(), "VAR_SYNONYMS")
                .map(|whole| String::from_utf8(whole).unwrap())
        };
        let body =
            |tail: &str, kept: &str| Some(format!("S::{}{kept}", "9".repeat(CAP - 3 - tail.len())));
        assert_eq!(synonyms(",R1,R"), body(",R1,R", ",R1"));
        // A group whose only id is cut goes with its source name.
        assert_eq!(synonyms("--A::1"), body("--A::1", ""));
        assert_eq!(synonyms("--A::1,2"), body("--A::1,2", "--A::1"));
        assert_eq!(synonyms("--A"), body("--A", ""));
        assert_eq!(synonyms("--A:"), body("--A:", ""));
        assert_eq!(synonyms("--A::"), None);
        assert_eq!(synonyms("--"), None);
        assert_eq!(
            uncaught(|| cut_list_whole_items(&capped(""), "OTHER")),
            "'OTHER'"
        );
    }

    #[test]
    fn a_row_value_is_read_as_the_bytes_of_the_stored_values_it_joins() {
        assert_eq!(row_list_bytes("-", false), b"");
        assert_eq!(row_list_bytes("", true), b"");
        assert_eq!(row_list_bytes(crate::shared::ABSENT, true), b"");
        assert_eq!(row_list_bytes("a%3Bb&c", false), b"a;b&c");
        assert_eq!(row_list_bytes("a%3Bb&c", true), b"a;b,c");
        assert_eq!(row_list_bytes("a%3bb", true), b"a%3bb");
    }

    #[test]
    fn an_extra_text_carries_the_candidates_value_of_each_cut_list() {
        let cut = format!("{},7", "9".repeat(CAP - 2));
        let kept = "9".repeat(CAP - 2);
        let e = row();
        assert_eq!(
            cut_list_reading(
                &e,
                "Extra",
                &format!("A=1;PUBMED={cut};;B"),
                &format!("PUBMED={kept};A=2")
            ),
            Some(format!("A=1;PUBMED={kept};;B"))
        );
        // A list left with no id is a pair the candidate does not write.
        let alone = "9".repeat(CAP);
        assert_eq!(
            cut_list_reading(&e, "Extra", &format!("A=1;PUBMED={alone}"), "A=2"),
            Some("A=1".to_string())
        );
        assert_eq!(
            cut_list_reading(&e, "Extra", &format!("PUBMED={alone}"), "-"),
            Some("-".to_string())
        );
        assert_eq!(
            cut_list_reading(&e, "Extra", &format!("PUBMED={cut}"), "PUBMED=9"),
            None
        );
        assert_eq!(cut_list_reading(&e, "Extra", "A=1", "A=2"), None);
        assert_eq!(
            cut_list_example(&e, "PUBMED", "9\u{e9}9", "-"),
            "21:500|G|ENST1|Transcript: PUBMED expected 3 characters, got 0"
        );
    }

    #[test]
    fn a_value_of_more_runs_than_the_interpreter_has_frames_for_ends_the_run() {
        let runs = |count: usize, join: &str| vec!["9".repeat(CAP); count].join(join);
        let reading = |name: &'static str, ev: String, av: &'static str| {
            move || cut_list_reading(&row(), name, &ev, av)
        };
        assert_eq!(
            reading("PUBMED", runs(991, ","), "-")(),
            Some("-".to_string())
        );
        for (column, join) in [("PUBMED", ","), ("VAR_SYNONYMS", "--")] {
            assert_eq!(
                uncaught(reading(column, runs(992, join), "-")),
                "maximum recursion depth exceeded"
            );
        }
        // The deepest call reads a value that is no run of the class.
        let short = format!("{},{}", runs(991, ","), "9".repeat(CAP - 1));
        assert_eq!(reading("PUBMED", short, "x")(), None);
    }

    #[test]
    fn a_json_member_nested_deeper_than_the_interpreter_walks_ends_the_run() {
        let nested = |levels: usize, leaf: &str| {
            format!(
                "{{\"id\":\"rs1\",\"pubmed\":{}{leaf}{}}}",
                "[".repeat(levels),
                "]".repeat(levels)
            )
        };
        assert!(!json_cut_lists_pair(&nested(496, "1"), &nested(496, "2")));
        assert_eq!(
            uncaught(|| json_cut_lists_pair(&nested(497, "1"), &nested(497, "2"))),
            "maximum recursion depth exceeded"
        );
    }

    #[test]
    fn a_member_the_candidate_lacks_is_refused_however_deep_it_is_nested() {
        // A member of the reference alone is refused without a walk of it, at a nesting
        // a recursive walk would not end on.
        let nested = format!("{}1{}", "[".repeat(80_000), "]".repeat(80_000));
        let lacking = "{\"id\":\"rs1\"}";
        let pubmed = format!("{{\"id\":\"rs1\",\"pubmed\":[{nested},2]}}");
        assert!(!json_cut_lists_pair(&pubmed, lacking));
        let synonyms =
            format!("{{\"id\":\"rs1\",\"var_synonyms\":{{\"A\":[{nested}],\"B\":[\"x\"]}}}}");
        assert!(!json_cut_lists_pair(&synonyms, lacking));
    }
}
