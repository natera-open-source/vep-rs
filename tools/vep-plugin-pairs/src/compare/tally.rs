// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Compares the plugin values of a file pair, a shard at a time: per plugin, the keys
//! either output annotates in sorted order, counted as matching, differing or present on
//! one side only. A shard's tally holds its counts, the fields it met and its first
//! mismatches; the shards' tallies are joined into the pair's result.

use super::fields::Run;
use super::shards::{KeyParts, Record, Shard};
use super::text::py_float;
use crate::rules::utrannotator::reordered_pair_groups;
use crate::rules::ORDER_ONLY_PAIR_GROUP_COLUMN;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::sync::Arc;

/// How a key that is not a match differs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Difference {
    Mismatch,
    PerlOnly,
    RustOnly,
}

impl Difference {
    pub fn name(self) -> &'static str {
        match self {
            Difference::Mismatch => "mismatch",
            Difference::PerlOnly => "perl_only",
            Difference::RustOnly => "rust_only",
        }
    }
}

/// One key that is not a match, with each side's fields in its row's order.
pub struct Mismatch {
    pub key: [String; 3],
    pub difference: Difference,
    pub perl: Vec<(String, String)>,
    pub rust: Vec<(String, String)>,
}

/// One plugin's counts over a file pair.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Stats {
    pub total_keys: u64,
    pub matched: u64,
    pub mismatch: u64,
    pub perl_only: u64,
    pub rust_only: u64,
    pub both_empty: u64,
    /// `None` where neither output annotates a key: unmeasured, never a perfect score.
    pub concordance_pct: Option<f64>,
    pub annotated_keys: u64,
    pub perl_annotated: u64,
    pub rust_annotated: u64,
}

/// One plugin's result over a file pair.
pub struct PluginResult {
    /// The fields either output carries, or the entries looked for when neither carries one.
    pub fields: Vec<String>,
    /// The compared fields that are entries of another requested plugin as well.
    pub shared_fields: Vec<String>,
    pub perl_annotated_records: u64,
    pub rust_annotated_records: u64,
    pub stats: Stats,
    /// The first mismatches in key order, as many as the run keeps.
    pub mismatches_sample: Vec<Mismatch>,
    pub total_mismatches: u64,
}

/// One side's fields under one key: their numbers in the run and their values.
type Data<'a> = Vec<(u32, &'a str)>;

/// Whether two sides of one key agree: every field of either side has the same text on
/// both, an absent field read as empty, or both texts are numbers less than the tolerance
/// apart. The two texts of the pair-group field agree as well where they hold the same
/// groups of the same pairs in another order, as the scorer's class
/// `utrannotator_annotation_pair_order` reads them: UTRAnnotator.pm:282-283 walk the groups
/// and each group's pairs as hashes.
fn agree(perl: &Data, rust: &Data, terms: &Terms) -> bool {
    fn value_in<'a>(data: &Data<'a>, field: u32) -> Option<&'a str> {
        data.iter()
            .find(|(number, _)| *number == field)
            .map(|(_, value)| *value)
    }
    let same = |field: u32, perl_value: &str, rust_value: &str| {
        perl_value == rust_value
            || matches!(
                (py_float(perl_value), py_float(rust_value)),
                (Some(a), Some(b)) if (a - b).abs() < terms.tolerance
            )
            || (&*terms.names[field as usize] == ORDER_ONLY_PAIR_GROUP_COLUMN
                && reordered_pair_groups(perl_value, rust_value))
    };
    perl.iter()
        .all(|(field, value)| same(*field, value, value_in(rust, *field).unwrap_or("")))
        && rust
            .iter()
            .filter(|(field, _)| value_in(perl, *field).is_none())
            .all(|(field, value)| same(*field, "", value))
}

/// One plugin's tally over one shard of a file pair.
#[derive(Default)]
pub struct Counting {
    stats: Stats,
    perl_records: u64,
    rust_records: u64,
    total_mismatches: u64,
    /// The shard's first mismatches in key order, as many as the run keeps.
    sample: Vec<Mismatch>,
    /// Per field number, whether either side carries the field.
    found: Vec<bool>,
}

/// What a comparison is made under.
pub struct Terms<'a> {
    /// The run's field names by their numbers.
    pub names: &'a [Arc<str>],
    pub tolerance: f64,
    pub max_mismatches: usize,
}

impl Counting {
    fn take(&mut self, key: KeyParts, perl: Option<&Data>, rust: Option<&Data>, terms: &Terms) {
        self.stats.total_keys += 1;
        if self.found.is_empty() {
            self.found.resize(terms.names.len(), false);
        }
        for (number, _) in perl.into_iter().chain(rust).flatten() {
            self.found[*number as usize] = true;
        }
        let difference = match (perl, rust) {
            (Some(_), None) => {
                self.perl_records += 1;
                self.stats.perl_only += 1;
                Difference::PerlOnly
            }
            (None, Some(_)) => {
                self.rust_records += 1;
                self.stats.rust_only += 1;
                Difference::RustOnly
            }
            (Some(perl), Some(rust)) => {
                self.perl_records += 1;
                self.rust_records += 1;
                if agree(perl, rust, terms) {
                    self.stats.matched += 1;
                    return;
                }
                self.stats.mismatch += 1;
                Difference::Mismatch
            }
            (None, None) => return,
        };
        self.total_mismatches += 1;
        if self.sample.len() < terms.max_mismatches {
            let owned = |data: Option<&Data>| -> Vec<(String, String)> {
                data.into_iter()
                    .flatten()
                    .map(|(number, value)| {
                        (terms.names[*number as usize].to_string(), value.to_string())
                    })
                    .collect()
            };
            self.sample.push(Mismatch {
                key: [key.0.to_string(), key.1.to_string(), key.2.to_string()],
                difference,
                perl: owned(perl),
                rust: owned(rust),
            });
        }
    }
}

/// The tally of each distinct plugin of the run over one shard of a file pair: the two
/// sides' records walked together in the order of their keys.
pub fn compare_shard(perl: &Shard, rust: &Shard, run: &Run, terms: &Terms) -> Vec<Counting> {
    let mut counting: Vec<Counting> = run.distinct.iter().map(|_| Counting::default()).collect();
    let (mut perl_records, mut rust_records) =
        (perl.records().peekable(), rust.records().peekable());
    // The fields of the record in hand, a side each; the lists are filled anew per record.
    let (mut perl_fields, mut rust_fields): (Data, Data) = (Vec::new(), Vec::new());
    loop {
        let order = match (perl_records.peek(), rust_records.peek()) {
            (Some(perl_record), Some(rust_record)) => perl_record.cmp_place(rust_record),
            (Some(_), None) => Ordering::Less,
            (None, Some(_)) => Ordering::Greater,
            (None, None) => break,
        };
        let perl_record = (order != Ordering::Greater)
            .then(|| perl_records.next())
            .flatten();
        let rust_record = (order != Ordering::Less)
            .then(|| rust_records.next())
            .flatten();
        let Some(record) = perl_record.or(rust_record) else {
            break;
        };
        perl_fields.clear();
        rust_fields.clear();
        perl_fields.extend(perl_record.iter().flat_map(Record::fields));
        rust_fields.extend(rust_record.iter().flat_map(Record::fields));
        counting[record.plugin as usize].take(
            record.key,
            perl_record.map(|_| &perl_fields),
            rust_record.map(|_| &rust_fields),
            terms,
        );
    }
    counting
}

/// The result of each distinct plugin of the run over one file pair, in the run's order,
/// from the tallies of the pair's shards. A key lies in one shard, so the pair's first
/// mismatches are the first, in key order, of the shards' first mismatches.
pub fn join(shards: Vec<Vec<Counting>>, run: &Run, terms: &Terms) -> Vec<PluginResult> {
    let mut joined: Vec<Counting> = run.distinct.iter().map(|_| Counting::default()).collect();
    for shard in shards {
        for (total, counted) in joined.iter_mut().zip(shard) {
            total.stats.total_keys += counted.stats.total_keys;
            total.stats.matched += counted.stats.matched;
            total.stats.mismatch += counted.stats.mismatch;
            total.stats.perl_only += counted.stats.perl_only;
            total.stats.rust_only += counted.stats.rust_only;
            total.perl_records += counted.perl_records;
            total.rust_records += counted.rust_records;
            total.total_mismatches += counted.total_mismatches;
            total.sample.extend(counted.sample);
            if total.found.len() < counted.found.len() {
                total.found.resize(counted.found.len(), false);
            }
            for (met, here) in total.found.iter_mut().zip(counted.found) {
                *met |= here;
            }
        }
    }
    fn found_names<'a>(counted: &Counting, names: &'a [Arc<str>]) -> BTreeSet<&'a str> {
        counted
            .found
            .iter()
            .zip(names)
            .filter(|(met, _)| **met)
            .map(|(_, name)| &**name)
            .collect()
    }

    // A compared field that is an entry of two requested plugins is marked on both, the one
    // that owns its records and the one that does not.
    let contested: BTreeSet<&str> = joined
        .iter()
        .flat_map(|counted| found_names(counted, terms.names))
        .filter(|field| {
            run.distinct
                .iter()
                .filter(|plugin| plugin.writes.claims(field))
                .count()
                > 1
        })
        .collect();
    run.distinct
        .iter()
        .zip(joined)
        .map(|(plugin, mut counted)| {
            let found = found_names(&counted, terms.names);
            let mut stats = std::mem::take(&mut counted.stats);
            let annotated = stats.matched + stats.mismatch + stats.perl_only + stats.rust_only;
            stats.concordance_pct =
                (annotated > 0).then(|| stats.matched as f64 / annotated as f64 * 100.0);
            stats.annotated_keys = annotated;
            stats.perl_annotated = stats.matched + stats.mismatch + stats.perl_only;
            stats.rust_annotated = stats.matched + stats.mismatch + stats.rust_only;
            let mut fields: Vec<String> = plugin
                .registered
                .iter()
                .filter(|entry| found.contains(entry.as_str()))
                .cloned()
                .collect();
            fields.extend(
                found
                    .iter()
                    .filter(|field| !plugin.registered.iter().any(|entry| entry == *field))
                    .map(|field| field.to_string()),
            );
            if fields.is_empty() {
                fields = plugin.registered.clone();
            }
            let mut sample = std::mem::take(&mut counted.sample);
            sample.sort_by(|a, b| a.key.cmp(&b.key));
            sample.truncate(terms.max_mismatches);
            PluginResult {
                fields,
                shared_fields: contested
                    .iter()
                    .filter(|field| plugin.writes.claims(field))
                    .map(|field| field.to_string())
                    .collect(),
                perl_annotated_records: counted.perl_records,
                rust_annotated_records: counted.rust_records,
                stats,
                mismatches_sample: sample,
                total_mismatches: counted.total_mismatches,
            }
        })
        .collect()
}
