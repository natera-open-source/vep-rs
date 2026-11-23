// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `merged_same_dbid_collapse` over the discordant rows: vep-rs's rows of a
//! transcript the reference dropped under `--merged` because a transcript of the other set
//! with its `dbID` stood earlier in the record's buffer
//! (`AnnotationType/Transcript.pm:283-287`), and the reference's `intergenic_variant` row at
//! a record the pass left with no transcript. The rows alone cannot tell such a row from one
//! vep-rs writes too many, so the class is applied only where the comparison's collector
//! holds the transcripts a replay of the reference's buffers found dropped.

use super::breakend::INTERGENIC_FEATURE;
use super::one_sided::{OneSidedStats, Row};
use super::scope::{Group, ReferenceRecordScope};
use crate::input::for_each_row_at;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io;
use std::path::Path;

/// vep-rs's rows at one (file, Location, allele): those of the class, and the others.
#[derive(Default)]
struct RustRows {
    of_the_class: u64,
    others: u64,
}

/// The Feature of every row the reference wrote at the group of each of `named`, indices of
/// `rows` with their groups, read from its outputs under `reference_dir`. A row both engines
/// wrote is not among the discordant rows, so each output that holds such a group is read
/// once more.
pub(super) fn features_written_at(
    rows: &[Row],
    named: &[(usize, Group)],
    reference_dir: &Path,
) -> io::Result<HashMap<Group, HashSet<String>>> {
    let mut asked: HashMap<&str, HashSet<&[u8]>> = HashMap::new();
    for (idx, _) in named {
        let row = &rows[*idx];
        asked
            .entry(row.file_name.as_str())
            .or_default()
            .insert(row.location.as_bytes());
    }
    let mut written: HashMap<Group, HashSet<String>> = HashMap::new();
    let text = |cell: &[u8]| String::from_utf8_lossy(cell).into_owned();
    for (file_name, locations) in asked {
        let path = reference_dir.join(file_name);
        let with_path = |e: io::Error| io::Error::new(e.kind(), format!("{}: {e}", path.display()));
        let reference = File::open(&path).map_err(with_path)?;
        for_each_row_at(reference, &locations, |key, _| {
            let mut cells = key.split(|b| *b == b'\t');
            if let (Some(location), Some(_), Some(feature)) =
                (cells.next(), cells.next(), cells.next())
            {
                written
                    .entry((file_name.to_string(), text(location)))
                    .or_default()
                    .insert(text(feature));
            }
        })
        .map_err(with_path)?;
    }
    Ok(written)
}

/// Marks the rows of the class in `masked` and returns their statistics, or `None` for a
/// comparison that names no converted cache.
///
/// A row vep-rs alone writes is of the class when it is a transcript's, the replay found
/// that transcript dropped at the row's Location, and the reference wrote no row of the
/// transcript there under any allele or consequence: none among `rows` and, where
/// `reference_dir` names the reference's outputs, none that vep-rs wrote too, which those
/// outputs are read for. The reference's intergenic row is of the class when vep-rs wrote
/// rows at its Location and allele and each of them is.
pub fn apply_shared_dbid_rule(
    rows: &[Row],
    masked: &mut HashSet<usize>,
    scope: &ReferenceRecordScope,
    reference_dir: Option<&Path>,
) -> io::Result<Option<OneSidedStats>> {
    let Some(drops) = scope.shared_dbid_drops.as_ref() else {
        return Ok(None);
    };
    let mut stats = OneSidedStats::default();
    let mut records: HashSet<Group> = HashSet::new();
    let reference_features: HashSet<(&str, &str, &str)> = rows
        .iter()
        .filter(|row| row.source == "missing_in_rust")
        .map(|row| {
            (
                row.file_name.as_str(),
                row.location.as_str(),
                row.feature.as_str(),
            )
        })
        .collect();
    let mut rust_rows: HashMap<(&str, &str, &str), RustRows> = HashMap::new();
    let mut named: Vec<(usize, Group)> = Vec::new();
    for (idx, row) in rows.iter().enumerate() {
        if row.source != "extra_in_rust" {
            continue;
        }
        let group: Group = (row.file_name.clone(), row.location.clone());
        let at_key = rust_rows
            .entry((
                row.file_name.as_str(),
                row.location.as_str(),
                row.allele.as_str(),
            ))
            .or_default();
        let dropped_there = !masked.contains(&idx)
            && row.feature_type == "Transcript"
            && !scope.absent.contains(&group)
            && !scope.gvcf_absent.contains(&group)
            && drops
                .get(&group)
                .is_some_and(|stable_ids| stable_ids.contains(&row.feature))
            && !reference_features.contains(&(
                row.file_name.as_str(),
                row.location.as_str(),
                row.feature.as_str(),
            ));
        if dropped_there {
            named.push((idx, group));
        } else {
            at_key.others += 1;
        }
    }
    let written = match reference_dir {
        Some(reference_dir) if !named.is_empty() => {
            features_written_at(rows, &named, reference_dir)?
        }
        _ => HashMap::new(),
    };
    for (idx, group) in named {
        let row = &rows[idx];
        let at_key = rust_rows
            .entry((
                row.file_name.as_str(),
                row.location.as_str(),
                row.allele.as_str(),
            ))
            .or_default();
        if written
            .get(&group)
            .is_some_and(|stable_ids| stable_ids.contains(&row.feature))
        {
            at_key.others += 1;
            continue;
        }
        at_key.of_the_class += 1;
        stats.count += 1;
        stats.excluded += 1;
        stats.rust += 1;
        masked.insert(idx);
        records.insert(group);
    }
    for (idx, row) in rows.iter().enumerate() {
        if masked.contains(&idx)
            || row.source != "missing_in_rust"
            || row.feature != INTERGENIC_FEATURE
            || row.consequence_set != "intergenic_variant"
        {
            continue;
        }
        let key = (
            row.file_name.as_str(),
            row.location.as_str(),
            row.allele.as_str(),
        );
        if rust_rows
            .get(&key)
            .is_some_and(|at_key| at_key.of_the_class > 0 && at_key.others == 0)
        {
            stats.count += 1;
            stats.excluded += 1;
            stats.perl += 1;
            masked.insert(idx);
            records.insert((row.file_name.clone(), row.location.clone()));
        }
    }
    stats.records = records.len() as u64;
    Ok(Some(stats))
}

#[cfg(test)]
#[path = "shared_dbid_tests.rs"]
mod tests;
