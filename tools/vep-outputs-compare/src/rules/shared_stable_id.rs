// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `merged_same_stable_id_collapse` over the discordant rows: vep-rs's rows of a
//! transcript the reference dropped under `--refseq` or `--merged` because another
//! transcript of the record's buffer carries its stable id
//! (`AnnotationType/Transcript.pm:314-343`), and the reference's `intergenic_variant` row at
//! a record the pass left with no transcript. The rows alone cannot tell such a row from one
//! vep-rs writes too many, so the class is applied only where the comparison's collector
//! holds the transcripts a replay of the reference's buffers found dropped.

use super::breakend::INTERGENIC_FEATURE;
use super::one_sided::{OneSidedStats, Row};
use super::scope::{Group, ReferenceRecordScope};
use super::shared_dbid::features_written_at;
use std::collections::{HashMap, HashSet};
use std::io;
use std::path::Path;

/// vep-rs's rows at one (file, Location, allele): those of the class, and the others. A row
/// another rule took of a transcript the `dbID` pass dropped there is neither.
#[derive(Default)]
struct RustRows {
    of_the_class: u64,
    others: u64,
}

/// Marks the rows of the class in `masked` and returns their statistics, or `None` for a
/// comparison that names no converted cache.
///
/// The rows vep-rs alone writes of one transcript id at one Location and allele are of the
/// class when all of these hold: they are a transcript's; the replay found transcripts of
/// the id dropped within reach of the record at that Location, at least as many as the
/// rows; the reference wrote no row of the id there that vep-rs lacks; and where the
/// reference wrote a row of the id there, which `reference_dir` names its outputs to be read
/// for, the replay found the transcript the pass kept within reach as well. The reference's
/// intergenic row is of the class when vep-rs wrote rows at its Location and allele, one of
/// them is of the class and each other is of a transcript the `dbID` pass dropped there.
pub fn apply_shared_stable_id_rule(
    rows: &[Row],
    masked: &mut HashSet<usize>,
    scope: &ReferenceRecordScope,
    reference_dir: Option<&Path>,
) -> io::Result<Option<OneSidedStats>> {
    let Some(drops) = scope.shared_stable_id_drops.as_ref() else {
        return Ok(None);
    };
    let dropped_for_a_dbid = |group: &Group, feature: &String| {
        let drops = scope.shared_dbid_drops.as_ref();
        drops.is_some_and(|d| d.get(group).is_some_and(|ids| ids.contains(feature)))
    };
    let mut stats = OneSidedStats::default();
    let mut records: HashSet<Group> = HashSet::new();
    let reference_alone: HashSet<(&str, &str, &str)> = rows
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
        let of_a_transcript = row.feature_type == "Transcript";
        if masked.contains(&idx) {
            if !(of_a_transcript && dropped_for_a_dbid(&group, &row.feature)) {
                at_key.others += 1;
            }
            continue;
        }
        let dropped_there = of_a_transcript
            && !scope.absent.contains(&group)
            && !scope.gvcf_absent.contains(&group)
            && drops
                .get(&group)
                .is_some_and(|stable_ids| stable_ids.contains_key(&row.feature))
            && !reference_alone.contains(&(
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
    // The rows of one id at one (file, Location, allele) stand or fall together.
    let mut of_one_id: HashMap<(&str, &str, &str, &str), Vec<usize>> = HashMap::new();
    for (idx, _) in &named {
        let row = &rows[*idx];
        of_one_id
            .entry((
                row.file_name.as_str(),
                row.location.as_str(),
                row.allele.as_str(),
                row.feature.as_str(),
            ))
            .or_default()
            .push(*idx);
    }
    for ((file_name, location, allele, feature), indices) in of_one_id {
        let group: Group = (file_name.to_string(), location.to_string());
        let seen = drops
            .get(&group)
            .and_then(|stable_ids| stable_ids.get(feature));
        let reference_wrote = written
            .get(&group)
            .is_some_and(|stable_ids| stable_ids.contains(feature));
        let of_the_class = seen.is_some_and(|seen| {
            indices.len() as u64 <= seen.dropped && (seen.kept_in_reach || !reference_wrote)
        });
        let at_key = rust_rows.entry((file_name, location, allele)).or_default();
        if !of_the_class {
            at_key.others += indices.len() as u64;
            continue;
        }
        for idx in indices {
            at_key.of_the_class += 1;
            stats.count += 1;
            stats.excluded += 1;
            stats.rust += 1;
            masked.insert(idx);
        }
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
#[path = "shared_stable_id_tests.rs"]
mod tests;
