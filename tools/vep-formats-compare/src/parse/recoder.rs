// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The recoder parser.

use crate::classes::{self, ARCHIVE_DBSNP_SYNONYMS};
use crate::entry::{Entry, Fields};
use crate::error::{Error, Result};
use crate::input::dedup_keys;
use crate::manifest::Documented;
use crate::parse::json::json_compact;
use crate::pyfmt;
use crate::pyjson::{self, Value};
use std::collections::BTreeSet;
use std::sync::LazyLock;

pub const RECODER_TRANSCRIPT_KEYS: &[&str] = &["hgvsc", "hgvsp"];

/// The keys whose values name a transcript, protein or region by prefix, the ones the classes read.
pub const RECODER_NOTATION_KEYS: &[&str] = &["hgvsc", "hgvsp", "hgvsg", "spdi", "vcf_string"];

pub const ARCHIVE_DBSNP_SOURCE: &str = "Archive dbSNP: ";

pub const IDENTIFIER_INDEL_NOTATIONS: &str = "identifier_indel_transcript_notations";

pub const IDENTIFIER_INDEL_KEYS: &[&str] = &["hgvsc", "hgvsp", "mane_select"];

pub const INDEL_NOTATION_PATTERN: &str = r"del|dup|ins";
pub static INDEL_NOTATION: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(INDEL_NOTATION_PATTERN).expect("the pattern compiles"));

/// The call depth `recoder_value` runs at where a comparison parses a recoder output. A
/// value ends the run where the calls that write it as text (`pyjson::compact_frames`)
/// would pass the interpreter's recursion limit (`pyjson::RECURSION_LIMIT`) less this
/// depth: objects around a whole number are written 495 deep and end the run at 496.
const RECODER_VALUE_CALLERS: usize = 7;

/// `json_compact` of a value `recoder_value` reads. A value nested past the frames the
/// interpreter has left ends the run and is no comparator's error.
fn compact(value: &Value) -> String {
    if matches!(value, Value::List(_) | Value::Dict(_))
        && RECODER_VALUE_CALLERS + pyjson::compact_frames(value) > pyjson::RECURSION_LIMIT
    {
        crate::uncaught!("maximum recursion depth exceeded");
    }
    json_compact(value)
}

/// Whether an applied `RECODER_MASKS` class drops `text`, counted under the first class of
/// the table that does.
fn masked_text(text: &str, documented: &mut Documented) -> bool {
    for (class, prefixes) in classes::tables().recoder_masks() {
        if documented.masks.contains(*class)
            && prefixes.iter().any(|prefix| text.starts_with(prefix))
        {
            documented.applied.add(class, 1);
            return true;
        }
    }
    false
}

/// One allele key's value in comparable form: a list as the sorted set of its items'
/// canonical texts, a `var_synonyms` entry with its identifiers sorted and a `mane_select`
/// triple as its canonical JSON; a scalar as its canonical text. An item an applied class
/// masks is dropped and counted under the class: under a `RECODER_MASKS` class, a string
/// under one of `RECODER_NOTATION_KEYS` starting with one of its prefixes, or a `mane_select`
/// triple whose `hgvsc` or `hgvsp` does; under `archive_dbsnp_synonyms`, a `var_synonyms`
/// entry of the `Archive dbSNP` source. A list the classes emptied is `None`, the key
/// compared as absent, since the other side writes no key for an empty list.
pub fn recoder_value(key: &str, value: &Value, documented: &mut Documented) -> Option<String> {
    let Some(list) = value.as_list() else {
        return Some(compact(value));
    };
    let mut items: BTreeSet<String> = BTreeSet::new();
    for item in list {
        if RECODER_NOTATION_KEYS.contains(&key)
            && item
                .as_str()
                .is_some_and(|text| masked_text(text, documented))
        {
            continue;
        }
        // A triple is dropped at its first masked notation, and counted once.
        if key == "mane_select"
            && RECODER_TRANSCRIPT_KEYS.iter().any(|name| {
                item.get(name)
                    .and_then(Value::as_str)
                    .is_some_and(|text| masked_text(text, documented))
            })
        {
            continue;
        }
        let synonym = match item.as_str() {
            Some(text) if key == "var_synonyms" => text.split_once(": ").map(|cut| (text, cut)),
            _ => None,
        };
        let Some((text, (source, identifiers))) = synonym else {
            items.insert(compact(item));
            continue;
        };
        if documented.masks.contains(ARCHIVE_DBSNP_SYNONYMS)
            && text.starts_with(ARCHIVE_DBSNP_SOURCE)
        {
            documented.applied.add(ARCHIVE_DBSNP_SYNONYMS, 1);
            continue;
        }
        let mut identifiers: Vec<&str> = identifiers.split(',').collect();
        identifiers.sort_unstable();
        let sorted = Value::Str(format!("{source}: {}", identifiers.join(",")));
        items.insert(compact(&sorted));
    }
    if !list.is_empty() && items.is_empty() {
        return None;
    }
    let items: Vec<String> = items.into_iter().collect();
    Some(format!("[{}]", items.join(",")))
}

/// Whether a list holds a text whose third and fourth `separator`-joined fields, the
/// deleted and inserted sequences of an SPDI or the reference and alternate of a VCF form,
/// differ in length; a value that is not a list holds none.
fn unbalanced(values: Option<&Value>, separator: char) -> bool {
    let Some(values) = values.and_then(Value::as_list) else {
        return false;
    };
    values.iter().filter_map(Value::as_str).any(|text| {
        let parts: Vec<&str> = text.split(separator).collect();
        parts.len() == 4 && pyfmt::len(parts[2]) != pyfmt::len(parts[3])
    })
}

/// Whether an allele object is an insertion or deletion of an identifier input: the input
/// has no notation separator and no space, and the allele is `-`, or an `hgvsg` names a
/// deletion, duplication or insertion, or an `spdi` deletes and inserts sequences of
/// different lengths, or a `vcf_string` has a reference and alternate of different lengths.
pub fn recoder_identifier_indel(label: &str, allele: &str, fields: &Value) -> bool {
    if label.contains(':') || label.chars().any(pyfmt::char_is_space) {
        return false;
    }
    if allele == "-" {
        return true;
    }
    let names_an_indel = |notation: &Value| {
        notation
            .as_str()
            .is_some_and(|h| INDEL_NOTATION.is_match(h))
    };
    if fields
        .get("hgvsg")
        .and_then(Value::as_list)
        .is_some_and(|hgvsg| hgvsg.iter().any(names_an_indel))
    {
        return true;
    }
    unbalanced(fields.get("spdi"), ':') || unbalanced(fields.get("vcf_string"), '-')
}

fn entry(key: String, allele: &str, fields: Fields) -> Entry {
    Entry {
        key,
        allele: allele.to_string(),
        fields,
        ..Entry::default()
    }
}

/// The fields read through `recoder_value`, a key it returns `None` for left out.
fn comparable_fields<'a>(
    fields: impl Iterator<Item = (&'a str, &'a Value)>,
    documented: &mut Documented,
) -> Fields {
    fields
        .filter_map(|(key, value)| Some((key.to_string(), recoder_value(key, value, documented)?)))
        .collect()
}

/// The recoder format: one entry per allele of each object of the JSON array, in allele
/// order, keyed by the object's `input` and the allele, plus one `|warnings` entry for an
/// object's warnings; a `null` object is one entry keyed `null`, an object of warnings alone
/// one keyed `warnings`. The keys `identifier_indel_transcript_notations` drops are left out
/// first; every other value is read through `recoder_value`.
pub fn parse_recoder(text: &str, documented: &mut Documented) -> Result<Vec<Entry>> {
    let objects = match pyjson::loads(text) {
        Ok(objects) => objects,
        Err(error) if error.is_value_error() => {
            return Err(Error::comparator(format!(
                "the recoder output does not parse: {error}"
            )))
        }
        Err(error) if error.is_recursion_error() => crate::uncaught!("{error}"),
        Err(_) => crate::uncaught!("the recoder output holds half a surrogate pair in a string"),
    };
    let Some(objects) = objects.as_list() else {
        return Err(Error::comparator("the recoder output is not a JSON array"));
    };
    let mut entries: Vec<Entry> = Vec::new();
    for (n, obj) in objects.iter().enumerate() {
        let n = n + 1;
        if matches!(obj, Value::Null) {
            entries.push(entry("null".to_string(), "", Fields::new()));
            continue;
        }
        let Some(obj) = obj.as_dict() else {
            return Err(Error::comparator(format!(
                "recoder object {n} is neither an object nor null"
            )));
        };
        let mut alleles: Vec<(&str, &Value)> = obj
            .iter()
            .filter(|(key, _)| *key != "warnings")
            .map(|(key, value)| (key.as_str(), value))
            .collect();
        let inputs: BTreeSet<&str> = alleles
            .iter()
            .filter_map(|(_, value)| value.get("input").and_then(Value::as_str))
            .collect();
        if inputs.len() > 1 {
            let inputs: Vec<String> = inputs.iter().map(|input| pyfmt::repr_str(input)).collect();
            return Err(Error::comparator(format!(
                "recoder object {n} names more than one input: [{}]",
                inputs.join(", ")
            )));
        }
        let label = inputs.first().copied().unwrap_or("warnings");
        // An object's allele keys come out of a Perl hash, in no order.
        alleles.sort_unstable_by_key(|(allele, _)| *allele);
        for (allele, fields) in alleles {
            let Some(members) = fields.as_dict() else {
                return Err(Error::comparator(format!(
                    "recoder object {n}: allele {} is not an object",
                    pyfmt::repr_str(allele)
                )));
            };
            let drops = documented.masks.contains(IDENTIFIER_INDEL_NOTATIONS)
                && recoder_identifier_indel(label, allele, fields);
            let dropped = |key: &str| drops && IDENTIFIER_INDEL_KEYS.contains(&key);
            if drops {
                let count = members.keys().filter(|key| dropped(key)).count();
                // Counted even where the allele holds none of the keys.
                documented
                    .applied
                    .add(IDENTIFIER_INDEL_NOTATIONS, count as i64);
            }
            let kept = members
                .iter()
                .filter(|(key, _)| !dropped(key))
                .map(|(key, value)| (key.as_str(), value));
            let fields = comparable_fields(kept, documented);
            entries.push(entry(format!("{label}|{allele}"), allele, fields));
        }
        if let Some(warnings) = obj.get("warnings") {
            let fields = comparable_fields(std::iter::once(("warnings", warnings)), documented);
            entries.push(entry(format!("{label}|warnings"), "", fields));
        }
    }
    dedup_keys(&mut entries);
    Ok(entries)
}
