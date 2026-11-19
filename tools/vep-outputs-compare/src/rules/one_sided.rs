// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The 116.2 one-sided pass over the discordant rows no pair rule removed.

use super::breakend::{
    bracket_mate, dot_breakend_points, is_breakend_allele, is_own_breakend_allele, location_span,
    non_integer_end, reference_mate_chromosome, strip_chr, INTERGENIC_FEATURE,
};
use super::scope::{Group, ReferenceRecordScope};
use super::table::{DIVERGENCE_CLASSES_116_2, PAIR_RULE_COUNT};
use crate::input::CACHE_CONTIG_NAMES;
use std::collections::{HashMap, HashSet};

pub const REFERENCE_SKIPPED_RECORD_BUCKET: &str = "reference_skipped_record";
pub const BND_OWN_ROWS_LOST_BUCKET: &str = "bnd_own_chromosome_rows_lost";
pub const BND_SYNONYM_MATE_ROWS_LOST_BUCKET: &str = "bnd_synonym_mate_rows_lost";
pub const FASTA_NAMED_SLICE_INTERGENIC_BUCKET: &str = "fasta_named_slice_intergenic";
pub const BND_BARE_REFERENCE_ALLELE_ROW_BUCKET: &str = "bnd_bare_reference_allele_row";
pub const BND_MATE_32_BIT_BUFFER_BUCKET: &str = "bnd_mate_32bit_buffer_intergenic";
pub const GVCF_VARIANT_SITE_DROPPED_BUCKET: &str = "gvcf_variant_site_dropped";
pub const DOT_BREAKEND_PAIR_BUCKET: &str = "dot_breakend_without_svtype_pair";
pub const NON_INTEGER_END_BUCKET: &str = "non_integer_end_annotated";
pub const NON_NUMERIC_POS_NUMIFIED_BUCKET: &str = "non_numeric_pos_numified";

/// One row of `discordant.tsv`.
pub struct Row {
    pub file_name: String,
    pub source: String,
    pub location: String,
    pub allele: String,
    pub feature: String,
    pub feature_type: String,
    pub consequence_set: String,
    /// Every column as read, for the open-rows file.
    pub fields: Vec<String>,
}

impl Row {
    fn group(&self) -> Group {
        (self.file_name.clone(), self.location.clone())
    }
}

/// A bucket's statistics: rows recognised, rows set aside, distinct (file, Location) pairs
/// reached, and the set-aside rows per side.
#[derive(Default, Clone)]
pub struct OneSidedStats {
    pub count: u64,
    pub excluded: u64,
    pub records: u64,
    pub perl: u64,
    pub rust: u64,
}

struct Tally<'a> {
    stats: Vec<(&'static str, OneSidedStats)>,
    records: HashMap<&'static str, HashSet<Group>>,
    masked: &'a mut HashSet<usize>,
}

/// The statistics of a registered bucket.
fn bucket_entry<'s>(
    stats: &'s mut [(&'static str, OneSidedStats)],
    bucket: &str,
) -> &'s mut OneSidedStats {
    &mut stats
        .iter_mut()
        .find(|(b, _)| *b == bucket)
        .expect("registered bucket")
        .1
}

impl Tally<'_> {
    fn hit(&mut self, bucket: &'static str, idx: usize, row: &Row, exclude: bool) {
        let entry = bucket_entry(&mut self.stats, bucket);
        entry.count += 1;
        self.records.entry(bucket).or_default().insert(row.group());
        if exclude {
            entry.excluded += 1;
            if row.source == "missing_in_rust" {
                entry.perl += 1;
            } else {
                entry.rust += 1;
            }
            self.masked.insert(idx);
        }
    }
}

/// Marks the set-aside rows in `masked`: vep-rs's rows of a skipped record and of a dropped
/// gVCF variant site, the reference's rows at a non-integer coordinate and at a POS numified
/// below zero, and both rows of a dot-breakend pair; counts without masking the two lost-row
/// breakend shapes, the FASTA-named slice shape and the reference's rows at a record whose
/// coordinate overflows the buffer. Returns the per-bucket statistics in registry order.
pub fn apply_116_2_one_sided_rules(
    rows: &[Row],
    masked: &mut HashSet<usize>,
    scope: &ReferenceRecordScope,
) -> Vec<(&'static str, OneSidedStats)> {
    let stats = DIVERGENCE_CLASSES_116_2[PAIR_RULE_COUNT..]
        .iter()
        .map(|class| (class.bucket, OneSidedStats::default()))
        .collect();
    let mut tally = Tally {
        stats,
        records: HashMap::new(),
        masked,
    };
    let bare_rows: u64 = scope.bare_reference_rows.values().sum();
    let bare = bucket_entry(&mut tally.stats, BND_BARE_REFERENCE_ALLELE_ROW_BUCKET);
    bare.count = bare_rows;
    bare.excluded = bare_rows;
    tally
        .records
        .entry(BND_BARE_REFERENCE_ALLELE_ROW_BUCKET)
        .or_default()
        .extend(scope.bare_reference_rows.keys().cloned());

    dot_breakend_and_coordinate_rows(rows, &mut tally, scope);
    let (perl_keys, rust_keys) = side_keys(rows, tally.masked);
    record_rules(rows, &mut tally, scope);
    fasta_named_slice(rows, &mut tally, &perl_keys, &rust_keys);

    let Tally {
        mut stats, records, ..
    } = tally;
    for (bucket, entry) in &mut stats {
        entry.records = records.get(bucket).map_or(0, |set| set.len() as u64);
    }
    stats
}

/// The dot-breakend pair, the reference's rows whose Location carries a non-integer
/// coordinate, and its rows at a Location vep-rs wrote nothing for whose start the reference
/// numified below zero.
fn dot_breakend_and_coordinate_rows(rows: &[Row], tally: &mut Tally, scope: &ReferenceRecordScope) {
    let mut point_rows: HashMap<(String, String, String, String), Vec<usize>> = HashMap::new();
    for (idx, row) in rows.iter().enumerate() {
        if tally.masked.contains(&idx)
            || row.source != "extra_in_rust"
            || !is_breakend_allele(&row.allele)
        {
            continue;
        }
        if let Some((_, start, end)) = location_span(&row.location) {
            if start == end {
                point_rows
                    .entry((
                        row.file_name.clone(),
                        row.location.clone(),
                        row.allele.clone(),
                        row.feature.clone(),
                    ))
                    .or_default()
                    .push(idx);
            }
        }
    }
    for (idx, row) in rows.iter().enumerate() {
        if tally.masked.contains(&idx) || row.source != "missing_in_rust" {
            continue;
        }
        if non_integer_end(&row.location) {
            tally.hit(NON_INTEGER_END_BUCKET, idx, row, true);
            continue;
        }
        let span = location_span(&row.location);
        if let Some((_, start, _)) = &span {
            if scope.reference_only.contains(&row.group()) && start.is_negative() {
                tally.hit(NON_NUMERIC_POS_NUMIFIED_BUCKET, idx, row, true);
                continue;
            }
        }
        if !is_breakend_allele(&row.allele) {
            continue;
        }
        let Some(span) = span else { continue };
        if span.1 >= span.2 {
            continue;
        }
        for point in dot_breakend_points(&span) {
            let key = (
                row.file_name.clone(),
                point,
                row.allele.clone(),
                row.feature.clone(),
            );
            if let Some(partners) = point_rows.remove(&key) {
                if partners.is_empty() {
                    continue;
                }
                tally.hit(DOT_BREAKEND_PAIR_BUCKET, idx, row, true);
                for rust_idx in partners {
                    tally.hit(DOT_BREAKEND_PAIR_BUCKET, rust_idx, &rows[rust_idx], true);
                }
                break;
            }
        }
    }
}

type SideKeys<'a> = HashMap<(&'a str, &'a str, &'a str), Vec<usize>>;

/// The one-sided rows of each side per (file, Location, allele), read before the record rules
/// mask anything further.
fn side_keys<'a>(rows: &'a [Row], masked: &HashSet<usize>) -> (SideKeys<'a>, SideKeys<'a>) {
    let mut perl_keys: SideKeys = HashMap::new();
    let mut rust_keys: SideKeys = HashMap::new();
    for (idx, row) in rows.iter().enumerate() {
        if masked.contains(&idx) {
            continue;
        }
        let key = (
            row.file_name.as_str(),
            row.location.as_str(),
            row.allele.as_str(),
        );
        match row.source.as_str() {
            "missing_in_rust" => perl_keys.entry(key).or_default().push(idx),
            "extra_in_rust" => rust_keys.entry(key).or_default().push(idx),
            _ => {}
        }
    }
    (perl_keys, rust_keys)
}

/// The record-level rules: the reference's rows at a buffer-overflowing coordinate, vep-rs's
/// rows at a dropped gVCF site or a skipped record, and the breakend rows the reference lost.
fn record_rules(rows: &[Row], tally: &mut Tally, scope: &ReferenceRecordScope) {
    let mut own_rows_lost: HashSet<&Group> = HashSet::new();
    for (group, breakends) in &scope.breakends {
        let own = strip_chr(group.1.split_once(':').map_or(group.1.as_str(), |(c, _)| c));
        if breakends
            .perl_alleles
            .iter()
            .any(|a| bracket_mate(a).is_some())
            && !breakends
                .perl_alleles
                .iter()
                .any(|a| is_own_breakend_allele(a))
            && breakends.mates.iter().any(|m| m != own)
        {
            own_rows_lost.insert(group);
        }
    }
    for (idx, row) in rows.iter().enumerate() {
        if tally.masked.contains(&idx) {
            continue;
        }
        let group = row.group();
        if row.source == "missing_in_rust" {
            if scope.mate_overflow.contains(&group) {
                tally.hit(BND_MATE_32_BIT_BUFFER_BUCKET, idx, row, false);
            }
            continue;
        }
        if row.source != "extra_in_rust" {
            continue;
        }
        if scope.gvcf_absent.contains(&group) {
            tally.hit(GVCF_VARIANT_SITE_DROPPED_BUCKET, idx, row, true);
            continue;
        }
        if scope.absent.contains(&group) {
            tally.hit(REFERENCE_SKIPPED_RECORD_BUCKET, idx, row, true);
            continue;
        }
        let Some(breakend) = scope.breakends.get(&group) else {
            continue;
        };
        if row.feature == INTERGENIC_FEATURE {
            continue;
        }
        if is_own_breakend_allele(&row.allele) {
            if own_rows_lost.contains(&group) {
                tally.hit(BND_OWN_ROWS_LOST_BUCKET, idx, row, false);
            }
            continue;
        }
        if let Some(mate) = reference_mate_chromosome(&row.allele) {
            if !breakend.perl_alleles.is_empty()
                && !breakend.perl_alleles.contains(&row.allele)
                && !CACHE_CONTIG_NAMES.contains(&mate)
            {
                tally.hit(BND_SYNONYM_MATE_ROWS_LOST_BUCKET, idx, row, false);
            }
        }
    }
}

/// The reference's only row at a key is its intergenic row while vep-rs wrote transcript rows
/// there and no intergenic row: a slice named after a FASTA whose names differ from the
/// cache's. Counted on the reference's row.
fn fasta_named_slice(rows: &[Row], tally: &mut Tally, perl_keys: &SideKeys, rust_keys: &SideKeys) {
    let mut keys: Vec<_> = perl_keys.iter().collect();
    keys.sort_by_key(|(_, idx)| idx[0]);
    for (key, perl_idx) in keys {
        let Some(rust_idx) = rust_keys.get(key) else {
            continue;
        };
        if rust_idx.is_empty() || perl_idx.len() != 1 {
            continue;
        }
        let perl_row = &rows[perl_idx[0]];
        if perl_row.feature == INTERGENIC_FEATURE
            && rust_idx
                .iter()
                .all(|i| rows[*i].feature != INTERGENIC_FEATURE)
        {
            tally.hit(
                FASTA_NAMED_SLICE_INTERGENIC_BUCKET,
                perl_idx[0],
                perl_row,
                false,
            );
        }
    }
}
