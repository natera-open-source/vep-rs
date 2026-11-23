// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The JSON parser and its value forms.

use crate::consts::JSON_CONSEQUENCE_LISTS;
use crate::entry::{Entry, Fields};
use crate::error::{Error, Result};
use crate::input::{dedup_keys, line_slices};
use crate::parse::text::COMPARISON_CALLERS;
use crate::parse::vcf::{line_key, vcf_line_key, InputIndex};
use crate::pyfmt;
use crate::pyjson::{self, Value};
use crate::self_parity::record_location_under;
use crate::shared::{
    allele_is_indel, is_structural_allele, reference_allele_is_empty, replaced_reference,
};
use std::panic::{self, AssertUnwindSafe};

/// The engine identity vep-rs writes in every JSON object and Ensembl VEP never writes, set
/// aside like `input` rather than compared; its header forms (`## vep-rs <version>`,
/// `##vep-rs=<version>`) fall with the other `##` lines the text parsers skip.
pub const JSON_IDENTITY_KEY: &str = "vep_rs_version";

/// The call depth `parse_json` runs at in the self-parity mode, one less than in a
/// comparison (`COMPARISON_CALLERS`). A value ends the run where the calls that read it
/// (`scalar_frames`) would pass the interpreter's recursion limit
/// (`pyjson::RECURSION_LIMIT`) less the depth the parser runs at: a text inside lists at
/// 497 lists in a comparison and at 498 in the self-parity mode.
pub const SELF_PARITY_CALLERS: usize = COMPARISON_CALLERS - 1;

/// A float as a JSON text, shortest round trip, the exponent without sign padding.
pub fn _json_float(f: f64) -> String {
    pyjson::json_float(f)
}

/// A JSON value as compact text with object keys sorted.
pub fn json_compact(v: &Value) -> String {
    pyjson::dumps_compact(v)
}

/// A JSON value as one comparable string: a string itself, an array its items joined
/// with `,`, anything else its compact JSON text.
pub fn json_scalar(v: &Value) -> String {
    pyjson::json_scalar(v)
}

/// The call depth `json_scalar(value)` reaches, counted as CPython counts the frames of
/// the recursive reading, its own first frame included: one for a string, two
/// a list around its items (the function and the generator that joins them), and its own
/// above those of `json_compact` (`pyjson::compact_frames`) for any other value.
pub(crate) fn scalar_frames(value: &Value) -> usize {
    let mut most = 0;
    let mut pending: Vec<(&Value, usize)> = vec![(value, 0)];
    while let Some((value, lists)) = pending.pop() {
        let own = match value {
            Value::Str(_) => 1,
            Value::List(items) => {
                pending.extend(items.iter().map(|item| (item, lists + 1)));
                2
            }
            other => 1 + pyjson::compact_frames(other),
        };
        most = most.max(2 * lists + own);
    }
    most
}

/// `json_scalar(value)` called by a function `frames` frames deep, where `frames` is
/// given: a value nested deeper than the frames the interpreter has left ends the run.
fn scalar_under(value: &Value, frames: Option<usize>) -> String {
    if let Some(frames) = frames {
        if frames + scalar_frames(value) > pyjson::RECURSION_LIMIT {
            crate::uncaught!("maximum recursion depth exceeded");
        }
    }
    pyjson::json_scalar(value)
}

/// The object of the line `line`, the `n`th of its output counted from 0.
///
/// A line the interpreter reads into a string that holds half a surrogate pair is refused
/// with the comparator's error.
fn json_object(line: &str, n: i64) -> Result<Value> {
    let number = n.saturating_add(1);
    match pyjson::loads(line) {
        Ok(value) if value.as_dict().is_some() => Ok(value),
        Ok(_) => Err(Error::comparator(format!(
            "JSON line {number} is not an object"
        ))),
        Err(error) if error.is_recursion_error() => crate::uncaught!("{error}"),
        Err(error) => Err(Error::comparator(format!(
            "JSON line {number} does not parse: {error}"
        ))),
    }
}

/// The (line ordinal, object) of every non-empty line of a JSON Lines output.
pub fn json_objects(text: &str) -> Result<Vec<(i64, Value)>> {
    let mut out = Vec::new();
    for (n, line) in (0i64..).zip(line_slices(text)) {
        if pyfmt::strip(line).is_empty() {
            continue;
        }
        out.push((n, json_object(line, n)?));
    }
    Ok(out)
}

/// The input record ordinal of a JSON object and its label. An object names its record
/// by the first five columns of its `input`, as a VCF line does: Ensembl VEP stores a
/// structural variant's line with an empty or `0` QUAL, FILTER or INFO written as `.`
/// (`Parser/VCF.pm:464-466`), so the whole `input` may differ from the line as read.
pub fn json_record_label(index: &mut InputIndex, n: i64, obj: &Value) -> (Option<i64>, String) {
    let key = obj
        .get("input")
        .and_then(Value::as_str)
        .map_or("", line_key);
    let record = index.record(key);
    let label = match record {
        Some(record) => format!("rec{record}"),
        None => format!("obj{n}"),
    };
    (record, label)
}

/// A whole number that is no `bool`; one past `i64` as the nearest `i64`.
fn whole_number(value: Option<&Value>) -> Option<i64> {
    match value? {
        Value::Int(number) => Some(*number),
        Value::BigInt(decimal) if decimal.starts_with('-') => Some(i64::MIN),
        Value::BigInt(_) => Some(i64::MAX),
        _ => None,
    }
}

/// A JSON object's `start` and the Location the default format writes for its rows
/// (`record_location`); `None` where `start` or `end` is not a whole number.
///
/// A `start` past `i64` is given as the nearest `i64`: the Location spells every digit of
/// such a bound, and `replaced_reference` reads the Location's bounds before the `start`.
pub fn json_object_span(obj: &Value) -> Option<(i64, String)> {
    object_span(obj, None)
}

/// The calls an object's Location is written above `parse_json`: its nested values are
/// read that much nearer the interpreter's recursion limit than the parser's own.
const LOCATION_FRAMES: usize = 2;

/// `json_object_span(obj)` called by `parse_json` where that stands `frames` frames deep
/// on the interpreter's call stack: the Location is written `LOCATION_FRAMES` higher, so a
/// `seq_region_name` that `parse_json` itself still reads may end the run there.
fn object_span(obj: &Value, frames: Option<usize>) -> Option<(i64, String)> {
    let start = whole_number(obj.get("start"))?;
    whole_number(obj.get("end"))?;
    let callers = frames.map(|frames| frames + LOCATION_FRAMES - 1);
    Some((start, record_location_under(obj, callers)))
}

/// The JSON format: one record entry per object (the top-level scalars other than the
/// engine identity, keyed by the input record the `input` line names) and one entry per
/// element of each consequence list, keyed by the record, `variant_allele`, the feature id
/// and the list's feature type. An object whose `input` matches no record is keyed by its
/// output ordinal. A consequence entry carries the reference bases its allele stands in for
/// and the Location of the default format (`Entry.reference`, `Entry.allele_location`)
/// where the object's `start`, `end` and `allele_string` give the allele at that
/// span. A structural allele's entry stands in for no bases and carries the Location of
/// the object's own `seq_region_name`, `start` and `end`.
pub fn parse_json(text: &str, inputs: &[String]) -> Result<Vec<Entry>> {
    json_entries(line_slices(text), None, inputs, COMPARISON_CALLERS)
}

/// `parse_json` called with `callers` frames of the interpreter's call stack below it
/// (`COMPARISON_CALLERS`, `SELF_PARITY_CALLERS`): the frames the interpreter has left
/// decide which nested value ends the run.
pub fn parse_json_under(text: &str, inputs: &[String], callers: usize) -> Result<Vec<Entry>> {
    json_entries(line_slices(text), None, inputs, callers)
}

/// `parse_json` of some lines of an output, as `text_lines` gives them and in file order,
/// where `compare_files` parses an output. `ordinals` holds one number a line of `lines`,
/// a blank line included: the line's ordinal among all the lines of the whole output,
/// which labels an object that matches no input record and numbers a line that is no
/// object.
pub fn parse_json_lines(
    lines: &[String],
    ordinals: &[i64],
    inputs: &[String],
) -> Result<Vec<Entry>> {
    json_entries(
        lines.iter().map(String::as_str),
        Some(ordinals),
        inputs,
        COMPARISON_CALLERS,
    )
}

fn json_entries<'a>(
    lines: impl Iterator<Item = &'a str>,
    ordinals: Option<&[i64]>,
    inputs: &[String],
    callers: usize,
) -> Result<Vec<Entry>> {
    // `parse_json` stands one frame above its callers.
    let frames = callers + 1;
    let mut index = InputIndex::new(inputs, &mut vcf_line_key);
    let mut entries: Vec<Entry> = Vec::new();
    // A whole output numbers its own lines.
    let mut counted = 0i64..;
    let mut given = ordinals.map(|ordinals| ordinals.iter().copied());
    let mut objects = lines.filter_map(|line| {
        let n = match &mut given {
            None => counted.next(),
            Some(given) => given.next(),
        };
        let Some(n) = n else {
            return Some(Err(Error::comparator(
                "a line has no ordinal in its output",
            )));
        };
        if pyfmt::strip(line).is_empty() {
            return None;
        }
        Some(json_object(line, n).map(|obj| (n, line, obj)))
    });
    while let Some(read) = objects.next() {
        let (n, line, obj) = read?;
        let frames = may_reach_the_limit(line, &obj, frames + LOCATION_FRAMES).then_some(frames);
        let built = panic::catch_unwind(AssertUnwindSafe(|| {
            object_entries(&mut index, n, &obj, frames, &mut entries)
        }));
        if !matches!(built, Ok(Ok(()))) {
            // The errors come in the order of reading every line into an object before
            // an entry is built: a later line that is no object ends the run before this
            // object does.
            for later in objects.by_ref() {
                later?;
            }
        }
        built.unwrap_or_else(|ended| panic::resume_unwind(ended))?;
    }
    dedup_keys(&mut entries);
    Ok(entries)
}

/// Whether a value of `obj`, the object of `line`, may be nested deep enough for
/// `json_scalar`, called `frames` frames deep, to reach the interpreter's recursion limit:
/// a value under `n` lists and objects stacks at most `2 * n + 5` frames. The line opens
/// every list and object of its values, and its brackets are counted before the object is
/// walked.
fn may_reach_the_limit(line: &str, obj: &Value, frames: usize) -> bool {
    let reaches = |nested: usize| frames + 2 * nested + 5 > pyjson::RECURSION_LIMIT;
    reaches(memchr::memchr2_iter(b'[', b'{', line.as_bytes()).count())
        && reaches(pyjson::depth(obj))
}

/// Appends the record entry of `obj` and the entry of every item of its consequence
/// lists. `frames` is the depth of `parse_json` on the interpreter's call stack where a
/// value of the object may reach the recursion limit there.
fn object_entries(
    index: &mut InputIndex,
    n: i64,
    obj: &Value,
    frames: Option<usize>,
    entries: &mut Vec<Entry>,
) -> Result<()> {
    let (record, rec_label) = json_record_label(index, n, obj);
    let set_aside = |name: &str| {
        name == "input"
            || name == JSON_IDENTITY_KEY
            || JSON_CONSEQUENCE_LISTS.iter().any(|(list, _)| *list == name)
    };
    let mut fields = Fields::new();
    for (name, value) in obj.as_dict().into_iter().flatten() {
        if !set_aside(name) {
            fields.insert(name.clone(), scalar_under(value, frames));
        }
    }
    entries.push(Entry {
        key: format!("{rec_label}|record"),
        record,
        fields,
        ..Entry::default()
    });
    let allele_string = obj
        .get("allele_string")
        .map(|spelt| scalar_under(spelt, frames))
        .unwrap_or_default();
    let mut spelt = allele_string.split('/');
    let reference_allele = spelt.next().unwrap_or_default();
    let alternate_alleles: Vec<&str> = spelt.collect();
    let span = object_span(obj, frames);
    // The allele string is spelt without the base before `start`, on which an insertion's
    // Location begins: `N`, which stands for any base, takes that base's place in the
    // alleles `replaced_reference` reads, an allele of `-` being empty.
    let anchored = |spelt: &str| format!("N{}", if spelt == "-" { "" } else { spelt });
    let anchored_reference = anchored(reference_allele);
    let anchored_alternates: Vec<String> = alternate_alleles
        .iter()
        .map(|spelt| anchored(spelt))
        .collect();
    let anchored_alternates: Vec<&str> = anchored_alternates.iter().map(String::as_str).collect();
    for (list_name, kind) in JSON_CONSEQUENCE_LISTS {
        let Some(items) = obj.get(list_name).and_then(Value::as_list) else {
            continue;
        };
        let feature_id = format!("{kind}_id");
        for item in items {
            let Some(members) = item.as_dict() else {
                return Err(Error::comparator(format!(
                    "{rec_label}: a {list_name} element is not an object"
                )));
            };
            let mut cfields = Fields::with_capacity(members.len());
            for (name, value) in members {
                cfields.insert(name.clone(), scalar_under(value, frames));
            }
            let field = |name: &str| cfields.get(name).cloned().unwrap_or_default();
            let allele = field("variant_allele");
            let feature = if *kind == "intergenic" {
                String::new()
            } else {
                field(&feature_id)
            };
            let reference = match &span {
                Some((start, location)) => replaced_reference(
                    start.saturating_sub(1),
                    &anchored_reference,
                    &anchored_alternates,
                    Some(location),
                    &allele,
                )?,
                None => None,
            };
            let allele_location = match (&span, &reference) {
                (Some((_, location)), Some(_)) => Some(location.clone()),
                (Some((_, location)), None) if is_structural_allele(&allele) => {
                    Some(location.clone())
                }
                _ => None,
            };
            entries.push(Entry {
                key: format!("{rec_label}|{allele}|{feature}|{kind}"),
                record,
                location: None,
                consequence: field("consequence_terms"),
                indel: allele_is_indel(reference_allele, &alternate_alleles, &allele, true),
                insertion: Some(reference_allele_is_empty(
                    reference_allele,
                    &alternate_alleles,
                    true,
                )),
                reference,
                line: None,
                allele_location,
                allele,
                feature,
                fields: cfields,
            });
        }
    }
    Ok(())
}
