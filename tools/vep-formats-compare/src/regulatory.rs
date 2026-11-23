// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The regulatory-row filters.

use crate::compare::key_kind;
use crate::compare::key_record;
use crate::compare::record_key;
use crate::consts::CONSEQUENCE_RANK;
use crate::consts::REGULATORY_KINDS;
use crate::entry::Entry;
use crate::error::Error;
use crate::error::Result;
use crate::manifest::repr_text_list;
use crate::self_parity::most_severe;
use crate::shared::feature_kind;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;

/// A RegulatoryFeature or MotifFeature consequence entry in any format: the key's last
/// segment is the Feature_type (default, tab, VCF) or the JSON list's kind.
pub fn is_regulatory_row(e: &Entry) -> bool {
    REGULATORY_KINDS.contains(&feature_kind(key_kind(&e.key)).as_str())
}

/// Per entry of `entries`, whether `filter_consequence_rows` keeps it; the record entries
/// it keeps of a record that lost a row carry the recomputed `most_severe_consequence`.
fn kept_rows(
    entries: &mut [Entry],
    keep: &mut dyn FnMut(&Entry) -> bool,
    keep_bare_records: bool,
) -> Result<Vec<bool>> {
    let mut terms_by_record: HashMap<String, BTreeSet<String>> = HashMap::new();
    let mut lost_rows: HashSet<String> = HashSet::new();
    for e in entries.iter() {
        if record_key(&e.key) {
            continue;
        }
        let record = key_record(&e.key);
        if keep(e) {
            terms_by_record
                .entry(record.to_string())
                .or_default()
                .extend(
                    e.consequence
                        .split(['&', ','])
                        .filter(|term| !term.is_empty())
                        .map(str::to_string),
                );
        } else {
            lost_rows.insert(record.to_string());
        }
    }
    let mut kept: Vec<bool> = Vec::with_capacity(entries.len());
    for e in entries.iter_mut() {
        if !record_key(&e.key) {
            kept.push(keep(e));
            continue;
        }
        let record = key_record(&e.key);
        if !keep_bare_records && !terms_by_record.contains_key(record) {
            kept.push(false);
            continue;
        }
        if lost_rows.contains(record) && e.fields.contains_key("most_severe_consequence") {
            let none = BTreeSet::new();
            let terms = terms_by_record.get(record).unwrap_or(&none);
            let unranked: Vec<&str> = terms
                .iter()
                .map(String::as_str)
                .filter(|term| !CONSEQUENCE_RANK.iter().any(|(ranked, _)| ranked == term))
                .collect();
            if !unranked.is_empty() {
                return Err(Error::comparator(format!(
                    "{}: no rank for the consequence terms {}",
                    e.key,
                    repr_text_list(&unranked)
                )));
            }
            e.fields
                .insert("most_severe_consequence".to_string(), most_severe(terms));
        }
        kept.push(true);
    }
    Ok(kept)
}

/// The entries whose consequence rows `keep` accepts, in their order, with the VCF/JSON
/// record entries: every one when `keep_bare_records`, else those of a record that kept a
/// row. A record that lost a row has its `most_severe_consequence` recomputed from the
/// rows it kept, as the candidate computes it from the rows it writes. Mutates the kept
/// record entries; a term `CONSEQUENCE_RANK` lacks is the comparator's error.
pub fn filter_consequence_rows(
    entries: &mut [Entry],
    keep: &mut dyn FnMut(&Entry) -> bool,
    keep_bare_records: bool,
) -> Result<Vec<Entry>> {
    let kept = kept_rows(entries, keep, keep_bare_records)?;
    Ok(entries
        .iter()
        .zip(kept)
        .filter(|(_, kept)| *kept)
        .map(|(e, _)| e.clone())
        .collect())
}

/// `filter_consequence_rows` on a list that keeps the result: the entries it does not keep
/// leave the list.
fn retain_consequence_rows(
    entries: &mut Vec<Entry>,
    keep: &mut dyn FnMut(&Entry) -> bool,
    keep_bare_records: bool,
) -> Result<()> {
    let mut kept = kept_rows(entries, keep, keep_bare_records)?.into_iter();
    entries.retain(|_| kept.next().unwrap_or(false));
    Ok(())
}

/// The reference entries without their regulatory rows (`regulatory_feature_rows`): a
/// record the reference annotates against a regulatory feature and no transcript compares
/// as its intergenic row, every record entry stays (`filter_consequence_rows`).
pub fn without_regulatory_rows(expected: &mut [Entry]) -> Result<Vec<Entry>> {
    filter_consequence_rows(expected, &mut |e| !is_regulatory_row(e), true)
}

/// The regulatory rows of either side and the record entries of the records that have
/// one (`--keep-only regulatory_feature_rows`; `filter_consequence_rows`).
pub fn only_regulatory_rows(entries: &mut [Entry]) -> Result<Vec<Entry>> {
    filter_consequence_rows(entries, &mut is_regulatory_row, false)
}

/// `without_regulatory_rows` on a list that keeps the result.
pub(crate) fn drop_regulatory_rows(expected: &mut Vec<Entry>) -> Result<()> {
    retain_consequence_rows(expected, &mut |e| !is_regulatory_row(e), true)
}

/// `only_regulatory_rows` on a list that keeps the result.
pub(crate) fn keep_regulatory_rows(entries: &mut Vec<Entry>) -> Result<()> {
    retain_consequence_rows(entries, &mut is_regulatory_row, false)
}
