// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Populated counts, Extra pairs, the organelle codon tables, field names and labels.

use crate::argparse::path_text;
use crate::cache;
use crate::counter::Counter;
use crate::entry::Entry;
use crate::entry::Fields;
use crate::error::{Error, Result};
use crate::genetic_codes::{reads_as_the_standard_code, CodonTable};
use crate::lists::same_items;
use crate::pyfmt;
use crate::pyjson::{self, Value};
use crate::rules::intron_ends::VERSION_SUFFIX_RE;
use crate::shared::ABSENT;
use indexmap::IndexMap;
use std::cell::RefCell;
use std::collections::BTreeSet;
use std::path::Path;
use std::path::PathBuf;

/// Per field, the entries carrying a value other than `-`.
pub fn populated_counts(entries: &[Entry]) -> Counter {
    let mut counts = Counter::new();
    for e in entries {
        for (name, value) in &e.fields {
            // `dash(value) != "-"`: an empty value is the default format's `-`.
            if !value.is_empty() && value != "-" {
                counts.add(name, 1);
            }
        }
    }
    counts
}

/// The `key=value` pairs of a default-format `Extra` value.
pub fn extra_pairs(extra: &str) -> Fields {
    let mut pairs = Fields::new();
    if extra == "-" || extra == ABSENT {
        return pairs;
    }
    for pair in extra.split(';') {
        let (key, value) = pair.split_once('=').unwrap_or((pair, ""));
        pairs.insert(key.to_string(), value.to_string());
    }
    pairs
}

/// The transcripts a candidate JSON cache files under a codon table that reads some codon
/// unlike the standard one, by stable id: a transcript's
/// `variation_effect_feature_cache.codon_table`, a number in a native cache and the dump's
/// string in a converted one, under which Ensembl VEP reads some codon unlike the standard
/// code (`reads_as_the_standard_code`, the reading the class's reference peptide is
/// translated under). A transcript under a code no table is held for is kept as the cache
/// spells its table. A chromosome's shards are read on the first lookup of an entry whose
/// Location names its directory, every shard on the first lookup of an entry without one (a
/// VCF or JSON entry).
#[derive(Debug)]
pub struct OrganelleCodonTables {
    pub root: PathBuf,
    pub _tables: RefCell<IndexMap<Option<String>, IndexMap<String, CodonTable>>>,
}

impl OrganelleCodonTables {
    pub fn new(cache: &Path) -> Result<Self> {
        let root = cache::joined(cache, "transcripts");
        if !root.is_dir() {
            return Err(Error::comparator(format!(
                "--json-cache: {} has no transcripts/ directory",
                path_text(cache)
            )));
        }
        Ok(OrganelleCodonTables {
            root,
            _tables: RefCell::new(IndexMap::new()),
        })
    }

    pub fn _read(&self, chromosome: Option<&str>) -> Result<IndexMap<String, CodonTable>> {
        let shards = cache::shards(&self.root, chromosome)?;
        let read = cache::map_shards(&shards, cache::threads(), |_, transcripts| {
            Ok(transcripts
                .iter()
                .filter_map(codon_table)
                .collect::<Vec<_>>())
        })?;
        let mut tables = IndexMap::new();
        for (stable_id, table) in read.into_iter().flatten() {
            tables.insert(stable_id, table);
        }
        Ok(tables)
    }

    /// The codon table of the entry's transcript, or `None` for a transcript the cache
    /// files under a table reading as the standard one or does not hold; a Feature carrying
    /// a version the cache's stable id lacks (`--transcript_version`) is looked up without
    /// it. A transcript under a code no table is held for ends the run by name.
    pub fn table(&self, e: &Entry) -> Result<Option<i64>> {
        let mut chromosome = e
            .location
            .as_deref()
            .filter(|location| !location.is_empty())
            .map(|location| location.split(':').next().unwrap_or("").to_string());
        if chromosome
            .as_deref()
            .is_some_and(|name| !cache::joined(&self.root, name).is_dir())
        {
            chromosome = None;
        }
        if !self._tables.borrow().contains_key(&chromosome) {
            let read = self._read(chromosome.as_deref())?;
            self._tables.borrow_mut().insert(chromosome.clone(), read);
        }
        let tables = self._tables.borrow();
        let tables = &tables[&chromosome];
        let filed = match tables.get(&e.feature) {
            Some(table) => Some(table),
            None => tables.get(VERSION_SUFFIX_RE.replace_all(&e.feature, "").as_ref()),
        };
        filed
            .map(|table| table.held(&e.feature).map_err(Error::comparator))
            .transpose()
    }
}

/// The stable id and codon table of a transcript the cache files under a table that reads
/// some codon unlike the standard one, or under a code no table is held for: a `codon_table`
/// that is absent, null or 0 names the standard code, and any other value that is not the
/// number of a held table is such a code, kept as the cache spells it.
fn codon_table(transcript: &Value) -> Option<(String, CodonTable)> {
    let table = match transcript
        .get("variation_effect_feature_cache")
        .filter(|cached| pyjson::truthy(cached))
    {
        None => None,
        Some(cached) => match cached.as_dict() {
            Some(members) => members
                .get("codon_table")
                .filter(|table| !matches!(table, Value::Null))
                .map(pyjson::py_str),
            None => crate::uncaught!("'{}' object has no attribute 'get'", python_type(cached)),
        },
    };
    let stable_id = transcript
        .get("stable_id")
        .filter(|stable_id| pyjson::truthy(stable_id))?;
    let text = table?;
    let number = match pyfmt::is_decimal(&text) {
        true => match pyfmt::int_checked(&text) {
            Ok(digits) => digits.parse::<i64>().ok(),
            Err(error) => crate::uncaught!("{error}"),
        },
        false => None,
    };
    let table = match number.map(|table| if table == 0 { 1 } else { table }) {
        Some(table) => match reads_as_the_standard_code(table) {
            Some(true) => return None,
            Some(false) => CodonTable::Held(table),
            None => CodonTable::Lacking(text),
        },
        None => CodonTable::Lacking(text),
    };
    match stable_id {
        Value::Str(stable_id) => Some((stable_id.clone(), table)),
        Value::List(_) | Value::Dict(_) => crate::uncaught!(
            "cannot use '{0}' as a dict key (unhashable type: '{0}')",
            python_type(stable_id)
        ),
        // A number is a key no Feature names.
        _ => None,
    }
}

/// The member `name` of a value read as an object: any other value ends the run with
/// `'<type>' object has no attribute 'get'`.
pub(crate) fn member<'a>(value: &'a Value, name: &str) -> Option<&'a Value> {
    match value.as_dict() {
        Some(members) => members.get(name),
        None => crate::uncaught!("'{}' object has no attribute 'get'", python_type(value)),
    }
}

/// The interpreter's name of the type `json.loads` gives a value.
pub(crate) fn python_type(value: &Value) -> &'static str {
    match value {
        Value::Null => "NoneType",
        Value::Bool(_) => "bool",
        Value::Int(_) | Value::BigInt(_) => "int",
        Value::Float(_) => "float",
        Value::Str(_) => "str",
        Value::List(_) => "list",
        Value::Dict(_) => "dict",
    }
}

/// Whether `actual` is `expected` with at least one differing `key=value` pair that
/// `accept(key, expected value, actual value)` allows substituted, the keys and their
/// order kept; a difference in the order of the pairs is not accepted.
pub fn extra_accepted(
    expected: &str,
    actual: &str,
    accept: &mut dyn FnMut(&str, &str, &str) -> bool,
) -> bool {
    let actual_pairs = extra_pairs(actual);
    let mut rewritten = String::with_capacity(actual.len());
    let mut substituted = false;
    for (index, pair) in expected.split(';').enumerate() {
        if index > 0 {
            rewritten.push(';');
        }
        let other = pair.split_once('=').and_then(|(key, value)| {
            let other = actual_pairs.get(key)?;
            (value != other && accept(key, value, other)).then_some((key, other))
        });
        match other {
            Some((key, other)) => {
                rewritten.push_str(key);
                rewritten.push('=');
                rewritten.push_str(other);
                substituted = true;
            }
            None => rewritten.push_str(pair),
        }
    }
    substituted && rewritten == actual
}

/// `extra_accepted` for the accepted field pairs.
pub fn extra_pairs_accepted(
    expected: &str,
    actual: &str,
    field_pairs: &IndexMap<String, BTreeSet<(String, String)>>,
) -> bool {
    extra_accepted(expected, actual, &mut |key, value, other| {
        field_pairs
            .get(key)
            .is_some_and(|pairs| pairs.contains(&(value.to_string(), other.to_string())))
    })
}

/// `extra_accepted` for the set-valued columns, their items agreeing as a multiset.
pub fn extra_order_only_accepted(expected: &str, actual: &str, order_only: &Fields) -> bool {
    extra_accepted(expected, actual, &mut |key, value, other| {
        order_only
            .get(key)
            .is_some_and(|separator| same_items(value, other, separator))
    })
}

/// Whether the two `Extra` texts differ only in pairs `accept(key, expected value,
/// actual value)` allows: a value that differs, a pair one side alone carries (its other
/// value read as `-`, since Ensembl VEP and the candidate both leave an undefined key out of
/// the text), the keys both sides carry in the same order. Unlike `extra_accepted` it admits
/// a pair the reference lacks, the shape of a field the candidate alone fills.
pub fn extra_class_accepted(
    expected: &str,
    actual: &str,
    accept: &mut dyn FnMut(&str, &str, &str) -> bool,
) -> bool {
    let expected_pairs = extra_pairs(expected);
    let actual_pairs = extra_pairs(actual);
    let keys: BTreeSet<&str> = expected_pairs
        .keys()
        .chain(actual_pairs.keys())
        .map(String::as_str)
        .collect();
    let mut changed = false;
    for key in keys {
        let value = expected_pairs.get(key).map_or("-", String::as_str);
        let other = actual_pairs.get(key).map_or("-", String::as_str);
        if value == other {
            continue;
        }
        if !accept(key, value, other) {
            return false;
        }
        changed = true;
    }
    changed
        && expected_pairs
            .keys()
            .filter(|key| actual_pairs.contains_key(*key))
            .eq(actual_pairs
                .keys()
                .filter(|key| expected_pairs.contains_key(*key)))
}

/// The value as a finite number, or `None` when it is not one.
pub fn as_number(value: &str) -> Option<f64> {
    pyfmt::float_of(value).filter(|number| number.is_finite())
}

/// Whether a field name is one of `patterns`, each a full name or a prefix ending in
/// `_`, the case of both set aside (the JSON format lower-cases every key).
pub fn field_named(name: &str, patterns: &[&str]) -> bool {
    if name.is_ascii() && patterns.iter().all(|pattern| pattern.is_ascii()) {
        let name = name.as_bytes();
        return patterns.iter().any(|pattern| {
            let pattern = pattern.as_bytes();
            if pattern.ends_with(b"_") {
                name.len() >= pattern.len() && name[..pattern.len()].eq_ignore_ascii_case(pattern)
            } else {
                name.eq_ignore_ascii_case(pattern)
            }
        });
    }
    let lowered = pyfmt::lower(name);
    patterns.iter().any(|pattern| {
        let pattern_lowered = pyfmt::lower(pattern);
        if pattern.ends_with('_') {
            lowered.starts_with(&pattern_lowered)
        } else {
            lowered == pattern_lowered
        }
    })
}

/// The input record an entry belongs to: the record ordinal of a VCF or JSON entry, and
/// in the text formats the Location with the `Uploaded_variation`, which tells two records
/// at one Location apart (two records sharing both are read as one, which only narrows the
/// rules that group a record's rows).
pub fn record_label(e: &Entry) -> String {
    let label = e.key.split('|').next().unwrap_or("");
    if e.location.is_none() {
        return label.to_string();
    }
    match e.fields.get("Uploaded_variation").map(String::as_str) {
        None | Some("") | Some("-") => label.to_string(),
        Some(uploaded) => format!("{label} {uploaded}"),
    }
}

/// The last segment of an entry's key without its duplicate counter: the Feature_type
/// (default, tab, VCF), the JSON list's kind, or `record`.
pub fn entry_kind(e: &Entry) -> String {
    let last = e.key.rsplit('|').next().unwrap_or("");
    last.split('#').next().unwrap_or("").to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(key: &str, location: Option<&str>, fields: &[(&str, &str)]) -> Entry {
        Entry {
            key: key.to_string(),
            location: location.map(str::to_string),
            fields: fields
                .iter()
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect(),
            ..Entry::default()
        }
    }

    #[test]
    fn an_extra_pair_one_side_alone_carries_is_read_against_a_dash() {
        let mut accept_any = |_: &str, _: &str, _: &str| true;
        assert!(extra_class_accepted(
            "IMPACT=LOW;STRAND=1",
            "IMPACT=LOW;STRAND=1;am_class=benign",
            &mut accept_any
        ));
        assert!(extra_class_accepted(
            "IMPACT=LOW;STRAND=1;GWAS_study=X",
            "IMPACT=LOW;STRAND=1",
            &mut accept_any
        ));
        assert!(extra_class_accepted(
            "IMPACT=LOW;Conservation=-5.2",
            "IMPACT=LOW;Conservation=-2.1",
            &mut accept_any
        ));
        assert!(!extra_class_accepted(
            "IMPACT=LOW;STRAND=1",
            "IMPACT=LOW;STRAND=1",
            &mut accept_any
        ));
        assert!(!extra_class_accepted(
            "IMPACT=LOW;STRAND=1",
            "STRAND=1;IMPACT=LOW",
            &mut accept_any
        ));
        let mut seen = Vec::new();
        extra_class_accepted(
            "IMPACT=LOW",
            "IMPACT=LOW;am_class=benign",
            &mut |key, value, other| {
                seen.push((key.to_string(), value.to_string(), other.to_string()));
                true
            },
        );
        assert_eq!(
            seen,
            [(
                "am_class".to_string(),
                "-".to_string(),
                "benign".to_string()
            )]
        );
        assert!(!extra_class_accepted(
            "IMPACT=LOW;STRAND=1",
            "IMPACT=HIGH;STRAND=1",
            &mut |key, _, _| key != "IMPACT"
        ));
    }

    #[test]
    fn extra_pairs_keep_the_first_place_and_the_last_value_of_a_repeated_key() {
        let pairs = extra_pairs("A=1;B;A=2;C=x=y;");
        assert_eq!(
            pairs
                .iter()
                .map(|(key, value)| (key.as_str(), value.as_str()))
                .collect::<Vec<_>>(),
            [("A", "2"), ("B", ""), ("C", "x=y"), ("", "")]
        );
        assert!(extra_pairs("-").is_empty());
        assert!(extra_pairs(ABSENT).is_empty());
        assert_eq!(extra_pairs("").len(), 1);
    }

    #[test]
    fn a_substituted_extra_pair_keeps_the_keys_and_their_order() {
        let mut accept_any = |_: &str, _: &str, _: &str| true;
        assert!(extra_accepted("A=1;B=2", "A=1;B=3", &mut accept_any));
        assert!(!extra_accepted("A=1;B=2", "B=3;A=1", &mut accept_any));
        assert!(!extra_accepted("A=1;B=2", "A=1;B=2", &mut accept_any));
        assert!(!extra_accepted("A=1;B=2", "A=1;B=3;C=4", &mut accept_any));
        assert!(!extra_accepted("A=1;B=2", "A=1;B=3", &mut |_, _, _| false));
        let mut pairs = IndexMap::new();
        pairs.insert(
            "B".to_string(),
            BTreeSet::from([("2".to_string(), "3".to_string())]),
        );
        assert!(extra_pairs_accepted("A=1;B=2", "A=1;B=3", &pairs));
        assert!(!extra_pairs_accepted("A=1;B=2", "A=1;B=4", &pairs));
        assert!(!extra_pairs_accepted("A=1;B=2", "A=5;B=3", &pairs));
    }

    #[test]
    fn a_populated_field_is_neither_empty_nor_a_dash() {
        let entries = [
            entry("k1", None, &[("B", "x"), ("A", "-"), ("C", "")]),
            entry("k2", None, &[("A", "y"), ("B", "z")]),
        ];
        let counts = populated_counts(&entries);
        assert_eq!(counts.iter().collect::<Vec<_>>(), [("B", 2), ("A", 1)]);
    }

    #[test]
    fn a_field_name_matches_a_full_name_or_a_prefix_whatever_the_case() {
        assert!(field_named("HGVSp", &["hgvsp"]));
        assert!(field_named("cadd_phred", &["CADD_"]));
        assert!(field_named("CADD_", &["cadd_"]));
        assert!(!field_named("cad", &["cadd_"]));
        assert!(!field_named("hgvsp_x", &["hgvsp"]));
        assert!(!field_named("x", &[]));
    }

    #[test]
    fn a_label_is_the_record_ordinal_or_the_location_with_the_uploaded_name() {
        assert_eq!(record_label(&entry("7|A|ENST1|Transcript", None, &[])), "7");
        let located = |fields: &[(&str, &str)]| {
            record_label(&entry("21:5|A|ENST1|Transcript", Some("21:5"), fields))
        };
        assert_eq!(located(&[]), "21:5");
        assert_eq!(located(&[("Uploaded_variation", "-")]), "21:5");
        assert_eq!(located(&[("Uploaded_variation", "")]), "21:5");
        assert_eq!(located(&[("Uploaded_variation", "rs1")]), "21:5 rs1");
        assert_eq!(
            entry_kind(&entry("21:5|A|ENST1|Transcript#2", None, &[])),
            "Transcript"
        );
        assert_eq!(entry_kind(&entry("7|record", None, &[])), "record");
        assert_eq!(entry_kind(&entry("plain", None, &[])), "plain");
    }

    #[test]
    fn a_number_is_a_finite_float_as_the_interpreter_reads_one() {
        assert_eq!(as_number("0.1261"), Some(0.1261));
        assert_eq!(as_number(" 1e3 "), Some(1000.0));
        assert_eq!(as_number("1_0"), Some(10.0));
        assert_eq!(as_number("-0"), Some(-0.0));
        for text in [
            "",
            "-",
            "abc",
            "nan",
            "inf",
            "-Infinity",
            "1e400",
            "1,5",
            "0x10",
        ] {
            assert_eq!(as_number(text), None, "{text}");
        }
    }

    #[test]
    fn a_codon_table_is_read_for_the_chromosome_of_the_first_row_that_asks() {
        use crate::cache::testing::{plant, scratch};
        let root = scratch("codon-tables");
        plant(
            &root,
            "transcripts/MT/1-5.json",
            r#"[{"stable_id": "ENST1", "variation_effect_feature_cache": {"codon_table": 2}},
                {"stable_id": "ENST2", "variation_effect_feature_cache": {"codon_table": "2"}},
                {"stable_id": "ENST3", "variation_effect_feature_cache": {"codon_table": 1}},
                {"stable_id": "ENST4", "variation_effect_feature_cache": {"codon_table": "11"}},
                {"stable_id": "ENST5", "variation_effect_feature_cache": {"codon_table": 2.0}},
                {"stable_id": "ENST6"},
                {"stable_id": "ENST8", "variation_effect_feature_cache": {"codon_table": 24}},
                {"stable_id": "ENST9", "variation_effect_feature_cache": {"codon_table": "26"}},
                {"stable_id": "ENST10", "variation_effect_feature_cache": {"codon_table": 0}},
                {"stable_id": "ENST11", "variation_effect_feature_cache": {"codon_table": null}},
                {"stable_id": "ENST12", "variation_effect_feature_cache": {"codon_table": 7}},
                {"stable_id": "ENST13", "variation_effect_feature_cache": {"codon_table": "13"}},
                null]"#,
        );
        plant(
            &root,
            "transcripts/21/1-5.json",
            r#"{"1-5": [{"stable_id": "ENST7", "variation_effect_feature_cache": {"codon_table": 4}}]}"#,
        );
        plant(&root, "transcripts/22/1-5.json", "[");
        let tables = OrganelleCodonTables::new(&root).unwrap();
        let table = |location: Option<&str>, feature: &str| {
            let mut row = entry("k", location, &[]);
            row.feature = feature.to_string();
            tables.table(&row).map_err(|refused| refused.to_string())
        };
        assert_eq!(table(Some("MT:7"), "ENST1"), Ok(Some(2)));
        assert_eq!(table(Some("MT:7"), "ENST2"), Ok(Some(2)));
        assert_eq!(table(Some("MT:7"), "ENST1.9"), Ok(Some(2)));
        assert_eq!(table(Some("MT:7"), "ENST13"), Ok(Some(13)));
        // The standard code, the codes that read every codon as it does (11 and 24), the
        // code Ensembl VEP reads under its last table (26), no table, table 0 and a null.
        for feature in [
            "ENST3", "ENST4", "ENST6", "ENST7", "ENST8", "ENST9", "ENST10", "ENST11", "ENST1.",
        ] {
            assert_eq!(table(Some("MT:7"), feature), Ok(None), "{feature}");
        }
        // A value that names no held table is refused by the transcript and the value.
        for (feature, spelt) in [("ENST5", "2.0"), ("ENST12", "7"), ("ENST12.3", "7")] {
            let refused = table(Some("MT:7"), feature).unwrap_err();
            assert!(
                refused.starts_with(&format!(
                    "transcript {feature}: the cache's codon_table {spelt} names no genetic \
                     code the comparison holds a table for (it holds 1, 2, "
                )),
                "{refused}"
            );
        }
        assert_eq!(table(Some("21:5-9"), "ENST7"), Ok(Some(4)));
        // The malformed shard of another chromosome is read by a row without a Location, or
        // one whose Location names no directory of the cache.
        for location in [None, Some(""), Some("nowhere:1")] {
            let refused = table(location, "ENST7").unwrap_err();
            assert!(
                refused.contains("transcripts/22/1-5.json is not JSON: "),
                "{refused}"
            );
        }
        assert_eq!(table(Some("21:5"), "ENST1"), Ok(None));
        let refused = OrganelleCodonTables::new(&root.join("absent"))
            .unwrap_err()
            .to_string();
        assert!(
            refused.starts_with("--json-cache: ")
                && refused.ends_with("absent has no transcripts/ directory")
        );
        std::fs::remove_dir_all(root).unwrap();
    }
}
