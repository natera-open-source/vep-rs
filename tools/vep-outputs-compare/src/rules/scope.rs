// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the one-sided passes need from the sorted merge, collected in the one read of both
//! key streams: the 115.2 collector of per-record transcript sets and the 116.2 collector
//! of per-Location record facts.

use super::breakend::{
    batch_dependent_annotation_scope, bracket_mate, is_own_breakend_allele, location_span,
    BND_MATE_32_BIT_LIMIT, BND_MATE_32_BIT_LIMIT_DIGITS, TANDEM_REPEAT_ALLELE,
};
use super::pair::is_structural_allele;
use std::collections::{HashMap, HashSet};

/// A (file name, Location) pair.
pub type Group = (String, String);

/// The first three columns of a key line.
fn key_columns(key_line: &str) -> (&str, &str, &str) {
    let mut it = key_line.splitn(4, '\t');
    let loc = it.next().unwrap_or("");
    let allele = it.next().unwrap_or("");
    let feature = it.next().unwrap_or("");
    (loc, allele, feature)
}

/// The features each engine named on every record in `batch_dependent_annotation_scope`,
/// keyed on (Location, allele); the 115.2 collector.
#[derive(Default)]
pub struct ScopeFeatures {
    pub records: HashMap<(String, String), (HashSet<String>, HashSet<String>)>,
}

impl ScopeFeatures {
    /// Adds the feature of a key line to its record's set for one side.
    pub fn record(&mut self, key_line: &str, rust: bool) {
        let (loc, allele, feature) = key_columns(key_line);
        if batch_dependent_annotation_scope(loc, allele) {
            let entry = self
                .records
                .entry((loc.to_string(), allele.to_string()))
                .or_default();
            let side = if rust { &mut entry.1 } else { &mut entry.0 };
            side.insert(feature.to_string());
        }
    }

    /// Adds the records another collector gathered.
    pub fn absorb(&mut self, part: ScopeFeatures) {
        for (record, (perl, rust)) in part.records {
            let entry = self.records.entry(record).or_default();
            entry.0.extend(perl);
            entry.1.extend(rust);
        }
    }
}

/// The alleles each side wrote at a breakend Location and the mate chromosomes its bracket
/// alleles name.
#[derive(Default)]
pub struct Breakends {
    pub perl_alleles: HashSet<String>,
    pub rust_alleles: HashSet<String>,
    pub mates: HashSet<String>,
}

/// What the reference's stable id pass did to the transcripts of one stable id, as one
/// record sees it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SharedStableId {
    /// The transcripts of the id that the pass dropped and that lie in reach of the record.
    pub dropped: u64,
    /// Whether the transcript of the id that the pass kept lies in reach of the record.
    pub kept_in_reach: bool,
}

/// The 116.2 collector. `absent`: the groups at which vep-rs wrote a structural row and the
/// reference nothing, less those `record_reference_named` takes out because a row there
/// carries the name of a record the reference wrote a row for: a record the reference
/// dropped as `vep_skip`; `gvcf_absent`: a group at a gVCF variant site that the reference
/// wrote nothing at; `breakends`: the breakend alleles of every group either
/// engine wrote one at; `bare_reference_rows`: the reference-allele rows of a bare
/// breakpoint the unique key stream dropped; `mate_overflow`: the groups whose reference row
/// carries a coordinate at or past `BND_MATE_32_BIT_LIMIT`; `reference_only`: the groups the
/// reference wrote and vep-rs nothing; `gvcf_sites`: the groups of vep-rs rows at a gVCF
/// variant site, handed in before the merge; `shared_dbid_drops`: for a comparison that
/// names the converted cache, the groups of the records in reach of a transcript the
/// reference dropped for a transcript of the other set with its `dbID`, each with the stable
/// ids dropped there, and `None` for a comparison that names no converted cache;
/// `shared_stable_id_drops`: for such a comparison, the groups of the records in reach of a
/// transcript the reference dropped for another transcript of its stable id, each with what
/// the record sees of the transcripts of each such id, and `None` likewise;
/// `shared_synonyms`: for a comparison that names the reference's cache, where that cache's
/// synonyms file gives a name to two or more sequences the cache holds, those names, and
/// `None` for any other comparison.
#[derive(Default)]
pub struct ReferenceRecordScope {
    pub absent: HashSet<Group>,
    pub gvcf_absent: HashSet<Group>,
    pub breakends: HashMap<Group, Breakends>,
    pub bare_reference_rows: HashMap<Group, u64>,
    pub mate_overflow: HashSet<Group>,
    pub gvcf_sites: HashSet<Group>,
    pub reference_only: HashSet<Group>,
    pub shared_dbid_drops: Option<HashMap<Group, HashSet<String>>>,
    pub shared_stable_id_drops: Option<HashMap<Group, HashMap<String, SharedStableId>>>,
    pub shared_synonyms: Option<HashSet<String>>,
    group: Option<Group>,
    reference_wrote: bool,
    rust_wrote: bool,
    rust_wrote_structural: bool,
}

impl ReferenceRecordScope {
    /// Observes one key line with the side or sides that wrote it.
    pub fn observe(&mut self, file_name: &str, key_line: &str, perl: bool, rust: bool) {
        let (loc, allele, _) = key_columns(key_line);
        let same_group = self
            .group
            .as_ref()
            .is_some_and(|(f, l)| f == file_name && l == loc);
        if !same_group {
            self.close();
            self.group = Some((file_name.to_string(), loc.to_string()));
        }
        let group = self.group.as_ref().expect("open group");
        if perl {
            self.reference_wrote = true;
            let coords = loc.split_once(':').map_or("", |(_, c)| c);
            if coords.len() >= BND_MATE_32_BIT_LIMIT_DIGITS
                && coords
                    .split('-')
                    .map(str::len)
                    .max()
                    .is_some_and(|n| n >= BND_MATE_32_BIT_LIMIT_DIGITS)
            {
                if let Some((_, start, end)) = location_span(loc) {
                    if start.max(end).at_least(BND_MATE_32_BIT_LIMIT) {
                        self.mate_overflow.insert(group.clone());
                    }
                }
            }
        }
        if rust {
            self.rust_wrote = true;
            if self.gvcf_sites.contains(group)
                || (is_structural_allele(allele) && allele != TANDEM_REPEAT_ALLELE)
            {
                self.rust_wrote_structural = true;
            }
        }
        let mate = bracket_mate(allele);
        if mate.is_some() || is_own_breakend_allele(allele) {
            let entry = self.breakends.entry(group.clone()).or_default();
            if perl {
                entry.perl_alleles.insert(allele.to_string());
            }
            if rust {
                entry.rust_alleles.insert(allele.to_string());
            }
            if let Some((chrom, position)) = mate {
                entry.mates.insert(chrom.to_string());
                if perl && position.at_least(BND_MATE_32_BIT_LIMIT) {
                    self.mate_overflow.insert(group.clone());
                }
            }
        }
    }

    /// Closes the open Location group; the merge calls it when its file ends.
    pub fn close(&mut self) {
        if let Some(group) = self.group.take() {
            if self.rust_wrote_structural && !self.reference_wrote {
                if self.gvcf_sites.contains(&group) {
                    self.gvcf_absent.insert(group.clone());
                } else {
                    self.absent.insert(group.clone());
                }
            }
            if self.reference_wrote && !self.rust_wrote {
                self.reference_only.insert(group);
            }
        }
        self.reference_wrote = false;
        self.rust_wrote = false;
        self.rust_wrote_structural = false;
    }

    /// Records one file's reference-allele rows of bare breakpoints, per Location.
    pub fn record_bare_reference_rows(
        &mut self,
        file_name: &str,
        rows_by_location: &[(String, u64)],
    ) {
        for (loc, rows) in rows_by_location {
            if *rows > 0 {
                self.bare_reference_rows
                    .insert((file_name.to_string(), loc.clone()), *rows);
            }
        }
    }

    /// Records one file's gVCF variant-site Locations.
    pub fn record_gvcf_sites<'a>(
        &mut self,
        file_name: &str,
        locations: impl IntoIterator<Item = &'a String>,
    ) {
        for loc in locations {
            self.gvcf_sites.insert((file_name.to_string(), loc.clone()));
        }
    }

    /// Records one file's transcripts dropped for a transcript of the other set with the
    /// same `dbID`, per Location.
    pub fn record_shared_dbid_drops(
        &mut self,
        file_name: &str,
        dropped: HashMap<String, HashSet<String>>,
    ) {
        let drops = self.shared_dbid_drops.get_or_insert_with(HashMap::new);
        for (location, stable_ids) in dropped {
            drops
                .entry((file_name.to_string(), location))
                .or_default()
                .extend(stable_ids);
        }
    }

    /// Records one file's transcripts dropped for another transcript of the same stable
    /// id, per Location. Two files' Locations stand apart, and what two records of one
    /// Location see of an id is kept where it holds for both.
    pub fn record_shared_stable_id_drops(
        &mut self,
        file_name: &str,
        dropped: HashMap<String, HashMap<String, SharedStableId>>,
    ) {
        let drops = self.shared_stable_id_drops.get_or_insert_with(HashMap::new);
        for (location, stable_ids) in dropped {
            let at_group = drops.entry((file_name.to_string(), location)).or_default();
            for (stable_id, seen) in stable_ids {
                at_group
                    .entry(stable_id)
                    .and_modify(|held| {
                        held.dropped = held.dropped.min(seen.dropped);
                        held.kept_in_reach &= seen.kept_in_reach;
                    })
                    .or_insert(seen);
            }
        }
    }

    /// Records the names the cache's synonyms file gives to two or more sequences the cache
    /// holds.
    pub fn record_shared_synonyms(&mut self, names: HashSet<String>) {
        self.shared_synonyms = Some(names);
    }

    /// One file's Locations in `absent`.
    pub fn absent_locations(&self, file_name: &str) -> HashSet<String> {
        self.absent
            .iter()
            .filter(|(file, _)| file == file_name)
            .map(|(_, location)| location.clone())
            .collect()
    }

    /// Takes out of `absent` one file's Locations at which a vep-rs row carries the name of
    /// a record the reference wrote a row for.
    pub fn record_reference_named(&mut self, file_name: &str, locations: &HashSet<String>) {
        self.absent
            .retain(|(file, location)| file != file_name || !locations.contains(location));
    }

    /// Adds what another collector, closed, observed of the key stream: the facts of the
    /// Location groups it read.
    pub fn absorb(&mut self, part: ReferenceRecordScope) {
        self.absent.extend(part.absent);
        self.gvcf_absent.extend(part.gvcf_absent);
        self.mate_overflow.extend(part.mate_overflow);
        self.reference_only.extend(part.reference_only);
        for (group, breakends) in part.breakends {
            let entry = self.breakends.entry(group).or_default();
            entry.perl_alleles.extend(breakends.perl_alleles);
            entry.rust_alleles.extend(breakends.rust_alleles);
            entry.mates.extend(breakends.mates);
        }
    }
}
