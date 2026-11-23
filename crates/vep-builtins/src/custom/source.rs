// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The source behind one `--custom` specification and the set a run holds.
//!
//! A source is opened once, before the input is read, and annotates each batch
//! in two phases as a built-in plugin does ([`crate::BuiltinPlugin`]): a
//! [`prefetch`](CustomSource::prefetch) over the shared batch, which the set
//! runs for every source in parallel, then a sequential
//! [`annotate`](CustomSource::annotate) that writes into the variants.
//!
//! The format-independent rules of `AnnotationSource/File.pm` live here: the
//! overlap predicates of `_record_overlaps_VF` (`:534-575`,
//! `base_overlap`), the record list and the statistics one input record
//! gathers over its matching records (`annotate_InputBuffer` `:303-355` and
//! `annotate_VariationFeature` `:389-443`, `MatchedRecords`), and the walk
//! over a batch's input records (`TabixSource`). Ensembl VEP annotates one
//! feature per input record, every allele at once, so a batch is walked record
//! group by record group and every allele of a group takes the group's list.
//!
//! A bigWig opens as `BigWigSource`, a `gencode_promoter` GFF as
//! `GffOverlaySource`, and a GFF in any other mode as a source that reads no
//! record and leaves each variant as it found it, since Ensembl VEP annotates
//! nothing from one (`GFF.pm:155-165`).

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::path::Path;
use std::sync::Arc;

use indexmap::IndexMap;
use rayon::prelude::*;
use tracing::warn;
use vep_core::custom_annotation::{
    CustomAnnotations, CustomFieldValue, CustomRecord, CustomStats, CustomValue,
};
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::InputVariant;
use vep_core::warning_hold::Held;

use super::bigwig::BigWigSource;
use super::gff_overlay::GffOverlaySource;
use super::grammar::{
    perl_number, CustomFormat, CustomSpec, GffType, MatchType, RecordLimit, SummaryStat,
};
use super::CustomRunOptions;
use crate::pool_pass::{collect_in_order, map_in_order, raise_again};
use crate::tabix::{self, BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord};
use crate::traits::{PluginError, PrefetchData};

/// Why a `--custom` source cannot serve a run.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum CustomError {
    /// The file cannot be opened: it is not there, has no index beside it, or
    /// is not of its format.
    Open { file: String, reason: String },
    /// A read of the file failed while a batch was annotated.
    Read { file: String, reason: String },
}

impl fmt::Display for CustomError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CustomError::Open { file, reason } => {
                write!(f, "cannot open the --custom source {file}: {reason}")
            }
            CustomError::Read { file, reason } => {
                write!(f, "cannot read the --custom source {file}: {reason}")
            }
        }
    }
}

impl std::error::Error for CustomError {}

/// One `--custom` source: its specification and the two annotation phases.
pub trait CustomSource: Send + Sync {
    /// The specification the source was opened from.
    fn spec(&self) -> &CustomSpec;

    /// The chromosome names the source's file lists, as the source reaches
    /// them: a tabix index's sequence names (`File.pm:369-372`), a bigWig's
    /// chromosomes less a leading `chr` (`BigWig.pm:147-149`). Ensembl VEP
    /// admits a record on any chromosome an annotation source lists
    /// (`BaseRunner.pm:444-460`).
    fn chromosomes(&self) -> Vec<String>;

    /// I/O phase over a shared batch: the records the batch can reach, kept
    /// for [`annotate`](Self::annotate). Runs in parallel with the other
    /// sources' prefetches.
    fn prefetch(&self, variants: &[InputVariant]) -> Result<PrefetchData, CustomError>;

    /// Annotation phase: writes the prefetched records' values into the
    /// variants. Runs source by source.
    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), CustomError>;

    /// The `Description` of each `##INFO` key of a VCF source whose file name
    /// ends in `.vcf`, `.vcf.gz` or `.vcf.bgz` (`OutputFactory/BaseTab.pm:86-101`
    /// read them for the header through that name test alone); `None` for
    /// every other source.
    fn field_descriptions(&self) -> Option<&IndexMap<String, String>> {
        None
    }
}

/// A source whose file is open and from which no record is annotated: a GFF
/// in a mode Ensembl VEP does not know, which it opens and then passes over
/// (`GFF.pm:155-165`), and a transcript model
/// ([`CustomSpec::is_transcript_source`]), whose transcripts the run's
/// transcript index reads. The batch passes through unannotated.
struct UnreadSource {
    spec: CustomSpec,
    sequence_names: HashSet<String>,
}

impl CustomSource for UnreadSource {
    fn spec(&self) -> &CustomSpec {
        &self.spec
    }

    fn chromosomes(&self) -> Vec<String> {
        self.sequence_names.iter().cloned().collect()
    }

    fn prefetch(&self, _variants: &[InputVariant]) -> Result<PrefetchData, CustomError> {
        Ok(PrefetchData::new(()))
    }

    fn annotate(
        &self,
        _variants: &mut [InputVariant],
        _data: PrefetchData,
    ) -> Result<(), CustomError> {
        Ok(())
    }
}

/// The tabix column layout of a format, as `tabix -p` presets it: the
/// chromosome, start and end columns, and whether the start is 0-based.
fn tabix_config(spec: &CustomSpec) -> TabixConfig {
    let (start_col, end_col, ref_col, alt_col, zero_based) = match spec.format {
        CustomFormat::Vcf => (1, Some(7), Some(3), Some(4), false),
        CustomFormat::Bed => (1, Some(2), None, None, true),
        CustomFormat::Gff | CustomFormat::Gtf => (3, Some(4), None, None, false),
        CustomFormat::BigWig => unreachable!("a bigWig is not a tabix file"),
    };
    TabixConfig {
        file_path: spec.file.clone().into(),
        chr_col: 0,
        start_col,
        end_col,
        ref_col,
        alt_col,
        zero_based,
    }
}

/// The engine's own words for a tabix failure, without the plugin framing
/// its error type adds.
pub(crate) fn tabix_reason(error: PluginError) -> String {
    match error {
        PluginError::Init(reason) | PluginError::Run(reason) => reason,
        PluginError::Io(error) => error.to_string(),
    }
}

/// Whether the sources pass over `variant`'s record: one Ensembl VEP skips
/// (`File.pm:318`, `vep_skip`) or one past `--max_sv_size`.
fn skipped(variant: &InputVariant) -> bool {
    variant.vep_skip || variant.oversize_sv
}

/// The input record one source annotates at a time: Ensembl VEP's variation
/// feature, which the batch holds as the adjacent alleles of one record
/// ([`InputVariant::same_record`]).
pub(crate) struct RecordGroup<'a> {
    /// The feature's start, `$vf->{start}`: one past its end for an insertion.
    pub start: u64,
    /// The feature's end, `$vf->{end}`.
    pub end: u64,
    /// The allele whose class, strand and structural kind are the feature's.
    pub lead: &'a InputVariant,
    /// The alternate alleles an `exact` VCF record is matched against
    /// (`$vf->alt_alleles`, `VCF.pm:409`): every alternate of the record, or
    /// the first alone for a `--minimal` split record, whose later split
    /// features lose their custom annotations in the rejoin
    /// (`OutputFactory.pm:2291-2349` copy the feature variations alone).
    pub alleles: Vec<&'a InputVariant>,
}

impl<'a> RecordGroup<'a> {
    /// The feature of `group`, or `None` for a record the sources pass over
    /// (`skipped`). The reference row of `--process_ref_homs` shares the
    /// record's feature and names no alternate allele.
    pub fn of(group: &'a [InputVariant]) -> Option<Self> {
        if group.iter().any(skipped) {
            return None;
        }
        let lead = group
            .iter()
            .find(|v| !v.reference_allele_row)
            .unwrap_or(group.first()?);
        let alleles = if lead.original_start.is_some() {
            vec![lead]
        } else {
            group.iter().filter(|v| !v.reference_allele_row).collect()
        };
        Some(RecordGroup {
            start: lead.start,
            end: lead.end,
            lead,
            alleles,
        })
    }
}

/// A record's start or end as Perl reads the cell ensembl-io hands Ensembl
/// VEP as written (`GXF.pm:163-166, 187-190`, `Bed.pm:159-162, 197-200`,
/// `BaseVCF4.pm:244-247, 307-309`): a whole number, or a number Perl holds as
/// a float, one with a fraction, one past the 64-bit range or a whole number
/// from 2^53 that a point spells without an exponent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Coordinate {
    /// A whole number inside the 64-bit range.
    Whole(i64),
    /// A finite number Perl holds as a float.
    Float(f64),
}

impl From<i64> for Coordinate {
    fn from(whole: i64) -> Self {
        Coordinate::Whole(whole)
    }
}

impl Coordinate {
    /// The number Perl reads from `cell` ([`perl_number`]): `1e2` is 100,
    /// `2000e-1` is 200, `100.5` and `1.25e1` keep their fraction, and a
    /// whole number is exact while it fits 64 bits, a float past them and,
    /// where a point spells it without an exponent (`100000000000000000.0`),
    /// from 2^53. `None` for a cell Perl reads as an infinity or as no number
    /// at all (`inf`, `nan`, `1e999`), which no coordinate holds.
    pub(crate) fn read(cell: &str) -> Option<Self> {
        let number = perl_number(cell);
        if !number.is_finite() {
            return None;
        }
        let text = cell.trim_start();
        let bytes = text.as_bytes();
        let sign = usize::from(matches!(bytes.first(), Some(b'+' | b'-')));
        let digits = bytes[sign..]
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count();
        let after = &bytes[sign + digits..];
        let opens_exponent = |rest: &[u8]| {
            matches!(
                rest,
                [b'e' | b'E', b'0'..=b'9', ..] | [b'e' | b'E', b'+' | b'-', b'0'..=b'9', ..]
            )
        };
        let exponent = digits > 0 && opens_exponent(after);
        if after.first() != Some(&b'.') && !exponent {
            // Perl reads the run of digits alone, exactly while it fits.
            return Some(match text[..sign + digits].parse() {
                Ok(whole) => Coordinate::Whole(whole),
                Err(_) if digits == 0 => Coordinate::Whole(0),
                Err(_) => Coordinate::Float(number),
            });
        }
        // Perl adds to a whole number an exponent spells as to an integer
        // (`-1e17` plus one is -99999999999999999) and to one a point alone
        // spells in the float, exact below 2^53 (`-100000000000000000.0`
        // plus one is itself).
        const EXACT: f64 = 9_007_199_254_740_992.0;
        const RANGE: std::ops::Range<f64> =
            -9_223_372_036_854_775_808.0..9_223_372_036_854_775_808.0;
        // Where no exponent follows the digits a point does, and the exponent
        // follows the point's digits (`1.5e2`).
        let exponent = exponent || {
            let fraction = after[1..]
                .iter()
                .take_while(|byte| byte.is_ascii_digit())
                .count();
            opens_exponent(&after[1 + fraction..])
        };
        let whole = number.fract() == 0.0
            && (number.abs() < EXACT || (exponent && RANGE.contains(&number)));
        Some(if whole {
            Coordinate::Whole(number as i64)
        } else {
            Coordinate::Float(number)
        })
    }

    /// The whole number Ensembl VEP names a record by and hands its compiled
    /// `overlap` (`%i`, `File.pm:489-495`; `:569`): a fraction is dropped
    /// toward zero, and a number past the 64-bit range is held at the end of
    /// the range, where Perl prints one past the upper end wrapped below zero.
    pub fn whole(self) -> i64 {
        match self {
            Coordinate::Whole(whole) => whole,
            // `as` drops the fraction toward zero and holds at the range's ends.
            Coordinate::Float(number) => number as i64,
        }
    }

    /// The number itself, as Perl's arithmetic reads it once a float is in
    /// the sum.
    pub fn value(self) -> f64 {
        match self {
            Coordinate::Whole(whole) => whole as f64,
            Coordinate::Float(number) => number,
        }
    }

    /// The number plus `count`, whole while the sum fits 64 bits, as Perl
    /// adds (`Bed.pm:161`, `BaseVCF4.pm:263, 270, 305, 314`).
    pub fn plus(self, count: i64) -> Self {
        match self {
            Coordinate::Whole(whole) => match whole.checked_add(count) {
                Some(sum) => Coordinate::Whole(sum),
                None => Coordinate::Float(whole as f64 + count as f64),
            },
            Coordinate::Float(number) => Coordinate::Float(number + count as f64),
        }
    }

    /// Whether the number lies past `position`, as Perl compares the two
    /// (`File.pm:335`).
    pub fn is_past(self, position: i64) -> bool {
        match self {
            Coordinate::Whole(whole) => whole > position,
            Coordinate::Float(number) => number > position as f64,
        }
    }
}

/// The predicates of `_record_overlaps_VF` on a record spanning
/// `ref_start..=ref_end` against the feature `vf_start..=vf_end`
/// (`File.pm:534-575`; the caller has applied the format's own shift to the
/// record's ends, `:524-532`). `None` when the record does not match, else
/// the overlap percentage as `PC` prints it: under `overlap`, `within` and
/// `surrounding`, the share of the record the feature covers, or of the
/// feature under `reciprocal`, to three decimals (`:552-568`); under `exact`,
/// which asks for the same start and the same end, `100` (`:571-574`). A
/// record with a float for an end is tested as Ensembl VEP tests one
/// (`float_overlap`).
pub(crate) fn base_overlap(
    spec: &CustomSpec,
    ref_start: impl Into<Coordinate>,
    ref_end: impl Into<Coordinate>,
    vf_start: i64,
    vf_end: i64,
) -> Option<String> {
    let (ref_start, ref_end) = match (ref_start.into(), ref_end.into()) {
        (Coordinate::Whole(start), Coordinate::Whole(end)) => (start, end),
        (start, end) => return float_overlap(spec, start, end, vf_start, vf_end),
    };
    if spec.match_type == MatchType::Exact {
        return (ref_start == vf_start && ref_end == vf_end).then(|| "100".to_string());
    }
    // An insertion's feature starts one past its end.
    let (vs, ve) = if vf_start > vf_end {
        (vf_end, vf_start)
    } else {
        (vf_start, vf_end)
    };
    let length = ve - vs + 1;
    if spec.match_type == MatchType::Within && (vs > ref_start || ve < ref_end) {
        return None;
    }
    if spec.match_type == MatchType::Surrounding && (vs < ref_start || ve > ref_end) {
        return None;
    }
    if let Some(distance) = spec.distance {
        // Written as Ensembl VEP tests it (`abs(..) <= $distance`), so a distance
        // that is not a number admits nothing.
        let within = |a: i64, b: i64| distance.admits(a.abs_diff(b));
        if !(within(vs, ref_start) && within(ve, ref_end)) {
            return None;
        }
    }
    // A start a minus sign leads is read as written, so a record's ends can
    // lie near the lowest integer: the count of shared bases and the record's
    // span saturate there.
    let shared = 1i64
        .saturating_add(ve.min(ref_end))
        .saturating_sub(vs.max(ref_start));
    let share_of = |span: i64| {
        if span != 0 {
            100.0 * shared as f64 / span as f64
        } else {
            100.0
        }
    };
    let percentage = if spec.reciprocal {
        share_of(length)
    } else {
        share_of(ref_end.saturating_sub(ref_start).saturating_add(1))
    };
    if percentage < spec.overlap_cutoff {
        return None;
    }
    // `overlap()` of `Utils/VariationEffect.pm`: the two spans share a base.
    (ref_end >= vs && ref_start <= ve).then(|| format!("{percentage:.3}"))
}

/// [`base_overlap`] for a record one of whose ends Perl holds as a float:
/// `within`, `surrounding`, `distance` and the percentage read the number as
/// written (`File.pm:541-566`), so a record from 100.5 does not surround a
/// feature at 100 and covers half of it; the last test, the compiled
/// `overlap` (`:569`), reads each end's whole number
/// ([`Coordinate::whole`]).
fn float_overlap(
    spec: &CustomSpec,
    ref_start: Coordinate,
    ref_end: Coordinate,
    vf_start: i64,
    vf_end: i64,
) -> Option<String> {
    let (start, end) = (ref_start.value(), ref_end.value());
    if spec.match_type == MatchType::Exact {
        return (start == vf_start as f64 && end == vf_end as f64).then(|| "100".to_string());
    }
    let (vs, ve) = (vf_start.min(vf_end), vf_start.max(vf_end));
    let (feature_start, feature_end) = (vs as f64, ve as f64);
    if spec.match_type == MatchType::Within && (feature_start > start || feature_end < end) {
        return None;
    }
    if spec.match_type == MatchType::Surrounding && (feature_start < start || feature_end > end) {
        return None;
    }
    if let Some(distance) = spec.distance {
        // Perl subtracts two whole numbers as integers and tests the whole
        // difference; only a float end makes the difference a float.
        let within = |feature: i64, bound: Coordinate| match bound {
            Coordinate::Whole(whole) => distance.admits(feature.abs_diff(whole)),
            Coordinate::Float(number) => distance.admits_float((feature as f64 - number).abs()),
        };
        if !(within(vs, ref_start) && within(ve, ref_end)) {
            return None;
        }
    }
    let shared = 1.0 + feature_end.min(end) - feature_start.max(start);
    let share_of = |span: f64| {
        if span != 0.0 {
            100.0 * shared / span
        } else {
            100.0
        }
    };
    let percentage = if spec.reciprocal {
        share_of(feature_end - feature_start + 1.0)
    } else {
        share_of(end - start + 1.0)
    };
    if percentage < spec.overlap_cutoff {
        return None;
    }
    (ref_end.whole() >= vs && ref_start.whole() <= ve).then(|| format!("{percentage:.3}"))
}

/// The records and statistics one source gathers for one input record while
/// its candidate records are walked in file order (`File.pm:303-443`).
pub(crate) struct MatchedRecords<'a> {
    spec: &'a CustomSpec,
    short_name: &'a Arc<str>,
    /// The entries the sources run before this one left on the input record:
    /// Ensembl VEP keeps one record list per short name, which every source
    /// of that name reads and extends (`File.pm:406, 437`).
    earlier: &'a [CustomAnnotations],
    /// Records that matched, the ones past the cap included (`:339`).
    record_count: u64,
    records: Vec<CustomRecord>,
    /// The `COSV` records whose namesake an earlier source listed under the
    /// short name: lent, never listed (`:405-415`).
    lent: Vec<CustomRecord>,
    stats: Option<StatsAccumulator>,
}

/// What one source's walk of one input record leaves on the record
/// ([`MatchedRecords::finish`]).
pub(crate) struct Listed {
    /// The source's entry: `None` when no record was listed and no statistic
    /// was asked for.
    pub entry: Option<CustomAnnotations>,
    /// The `COSV` records that lend their fields to the namesake an earlier
    /// source listed under the short name (`File.pm:409-413`).
    pub lent: Vec<CustomRecord>,
}

impl Listed {
    /// Leaves the walk's result on one allele of the record: each lent record
    /// gives its namesake the fields it lacks, then the entry is appended.
    fn leave_on(&self, short_name: &str, entries: &mut Vec<CustomAnnotations>) {
        for record in &self.lent {
            let namesake = entries
                .iter_mut()
                .filter(|entry| &*entry.short_name == short_name)
                .flat_map(|entry| entry.records.iter_mut())
                .find(|listed| listed.name == record.name);
            if let Some(namesake) = namesake {
                lend_fields(namesake, record.fields.clone());
            }
        }
        if let Some(entry) = &self.entry {
            entries.push(entry.clone());
        }
    }
}

/// Gives `namesake` each of `fields` it lacks (`File.pm:409-413`).
fn lend_fields(
    namesake: &mut CustomRecord,
    fields: impl IntoIterator<Item = (String, CustomFieldValue)>,
) {
    for (field, value) in fields {
        namesake.fields.entry(field).or_insert(value);
    }
}

/// The running statistics of `annotate_VariationFeature` (`File.pm:423-434`).
#[derive(Default)]
struct StatsAccumulator {
    /// The least score as a number and as the record carried it.
    min: Option<(f64, CustomValue)>,
    max: Option<(f64, CustomValue)>,
    /// Kept while `sum` or `mean` is asked for.
    sum: Option<f64>,
}

impl<'a> MatchedRecords<'a> {
    /// `earlier` are the entries already on the input record, the lists of
    /// the sources run before this one.
    pub fn new(
        spec: &'a CustomSpec,
        short_name: &'a Arc<str>,
        earlier: &'a [CustomAnnotations],
    ) -> Self {
        MatchedRecords {
            spec,
            short_name,
            earlier,
            record_count: 0,
            records: Vec::new(),
            lent: Vec::new(),
            stats: spec
                .summary_stats
                .as_ref()
                .map(|_| StatsAccumulator::default()),
        }
    }

    /// Whether a record's score is read (`$get_scores`, `File.pm:399`).
    pub fn scores(&self) -> bool {
        self.stats.is_some()
    }

    /// Whether the walk ends before the next candidate: past the cap, nothing
    /// more is listed, unless a statistic counts every record (`File.pm:337`).
    pub fn stops_before_next(&self) -> bool {
        self.stats.is_none()
            && matches!(self.spec.num_records, RecordLimit::Count(cap) if self.record_count > cap)
    }

    fn asks(&self, stat: SummaryStat) -> bool {
        self.spec
            .summary_stats
            .as_ref()
            .is_some_and(|stats| stats.contains(&stat))
    }

    /// Whether a source run before this one listed `name` under the short name
    /// (`File.pm:406` greps the one list the short name owns).
    fn listed_earlier(&self, name: &str) -> bool {
        self.earlier
            .iter()
            .filter(|entry| entry.short_name == *self.short_name)
            .flat_map(|entry| entry.records.iter())
            .any(|record| record.name == name)
    }

    /// One matching record's output records and its overlap percentage
    /// (`annotate_VariationFeature`, `File.pm:404-442`). A `COSV` record named
    /// like one already listed under the short name, by this source or an
    /// earlier one, lends it the fields it lacks and is listed no second time
    /// (`:405-415`); otherwise the first record takes `PC` when `fields` names
    /// it (`:418-420`), its score enters the statistics (`:423-434`), and the
    /// records are listed while the cap allows, one `...` at the cap
    /// (`:436-440`). Every matching record counts (`:339`).
    pub fn add(&mut self, mut records: Vec<CustomRecord>, percentage: Option<String>) {
        let Some(first) = records.first_mut() else {
            return;
        };
        let cosmic = first.name.starts_with("COSV");
        let listed_namesake = if cosmic {
            self.records.iter_mut().find(|r| r.name == first.name)
        } else {
            None
        };
        if let Some(listed) = listed_namesake {
            lend_fields(listed, first.fields.drain(..));
        } else if cosmic && self.listed_earlier(&first.name) {
            self.lent.push(records.swap_remove(0));
        } else {
            if self
                .spec
                .fields
                .as_ref()
                .is_some_and(|fields| fields.iter().any(|f| f == "PC"))
            {
                first.fields.insert(
                    "PC".to_string(),
                    CustomFieldValue::Text(percentage.unwrap_or_default()),
                );
            }
            if let Some(score) = first.score.clone() {
                self.add_score(score);
            }
            match self.spec.num_records {
                RecordLimit::All => self.records.append(&mut records),
                RecordLimit::Count(cap) if self.record_count < cap => {
                    self.records.append(&mut records)
                }
                RecordLimit::Count(cap) if self.record_count == cap => self.records.push(
                    CustomRecord::new("...".to_string(), None, IndexMap::new(), None),
                ),
                RecordLimit::Count(_) => {}
            }
        }
        self.record_count += 1;
    }

    /// `File.pm:425-432`: a defined score moves the least and greatest, which
    /// keep the score as the record carried it, and adds to the sum.
    fn add_score(&mut self, score: CustomValue) {
        let value = match &score {
            CustomValue::Text(text) => perl_number(text),
            CustomValue::Number(number) => *number,
            CustomValue::Undefined => return,
        };
        let (asks_min, asks_max) = (self.asks(SummaryStat::Min), self.asks(SummaryStat::Max));
        let asks_sum = self.asks(SummaryStat::Sum) || self.asks(SummaryStat::Mean);
        let Some(stats) = &mut self.stats else {
            return;
        };
        if asks_min && stats.min.as_ref().is_none_or(|(least, _)| value < *least) {
            stats.min = Some((value, score.clone()));
        }
        if asks_max
            && stats
                .max
                .as_ref()
                .is_none_or(|(greatest, _)| value > *greatest)
        {
            stats.max = Some((value, score));
        }
        if asks_sum {
            stats.sum = Some(stats.sum.unwrap_or(0.0) + value);
        }
    }

    /// The source's entry for the input record, `None` when no record was
    /// listed and no statistic was asked for, with the records it lends. The
    /// statistics close as `annotate_InputBuffer` closes them
    /// (`File.pm:344-352`): `count` is the number of matching records, `mean`
    /// the sum over that count where a sum exists, and the sum stays only when
    /// asked for.
    pub fn finish(self) -> Listed {
        let stats = self.stats.as_ref().map(|acc| {
            let mut values: IndexMap<String, CustomValue> = IndexMap::new();
            if let Some((_, least)) = &acc.min {
                values.insert("min".to_string(), least.clone());
            }
            if let Some((_, greatest)) = &acc.max {
                values.insert("max".to_string(), greatest.clone());
            }
            if let (true, Some(sum)) = (self.asks(SummaryStat::Mean), acc.sum) {
                values.insert(
                    "mean".to_string(),
                    CustomValue::Number(sum / self.record_count as f64),
                );
            }
            if self.asks(SummaryStat::Count) {
                values.insert(
                    "count".to_string(),
                    CustomValue::Number(self.record_count as f64),
                );
            }
            if let (true, Some(sum)) = (self.asks(SummaryStat::Sum), acc.sum) {
                values.insert("sum".to_string(), CustomValue::Number(sum));
            }
            // Ensembl VEP builds the hash through `count`, the `mean` test and
            // the deletion of an unasked `sum`; `sum` asked for alone, with no
            // score to add, leaves it undefined.
            if values.is_empty()
                && self.asks(SummaryStat::Sum)
                && !self.asks(SummaryStat::Count)
                && !self.asks(SummaryStat::Mean)
            {
                CustomStats::Undefined
            } else {
                CustomStats::Values(values)
            }
        });
        let entry = (!self.records.is_empty() || stats.is_some()).then(|| {
            CustomAnnotations::new(
                Arc::clone(self.short_name),
                self.records,
                self.record_count > 0,
                stats,
            )
        });
        Listed {
            entry,
            lent: self.lent,
        }
    }
}

/// What one candidate record is to the input record walked
/// (`annotate_InputBuffer`, `File.pm:333-341`).
pub(crate) enum Verdict {
    /// The record's start, as ensembl-io's `get_start` returns it, lies past
    /// the feature's end plus one: the walk ends here (`File.pm:333-335`).
    PastEnd,
    /// The record does not match the feature (`File.pm:396`).
    NoMatch,
    /// The output records the record yields (`_create_records`) and the
    /// overlap percentage, `None` for an `exact` VCF record (`VCF.pm:420`
    /// returns the matches alone).
    Match(Vec<CustomRecord>, Option<String>),
}

/// How one format reads a record of its file against an input record
/// (`_record_overlaps_VF` and `_create_records`, which each format class
/// defines over `File.pm`'s).
pub(crate) trait FormatRules: Send + Sync {
    /// The record's verdict against the feature; `scores` asks for the
    /// record's score (`File.pm:399`).
    fn match_record(&self, record: &TabixRecord, group: &RecordGroup<'_>, scores: bool) -> Verdict;
}

/// A source over a tabix-indexed file: the engine, the file's sequence names
/// and the synonyms that map an input chromosome onto them.
pub(crate) struct TabixSource {
    spec: CustomSpec,
    short_name: Arc<str>,
    annotator: TabixAnnotator,
    sequence_names: HashSet<String>,
    synonyms: Arc<ChromosomeSynonyms>,
}

/// Input records per match chunk: a batch is matched chunk by chunk on the
/// run's pool.
const RECORDS_PER_CHUNK: usize = 256;

/// The index's linear window: tabix resolves a position to a 16 kb window
/// whose first record's offset is where a read of that position begins, so
/// a read of any window inflates its 16 kb window from the start. A batch's
/// windows are read one 16 kb window at a time, each window once, in
/// parallel on the run's pool.
const LINEAR_WINDOW_SHIFT: u32 = 14;

/// One batch's records, read window by window.
struct Prefetched {
    /// The variant index at which each match chunk starts, then the batch's
    /// length.
    bounds: Vec<usize>,
    /// One query result per (chromosome, linear window) pair of the batch: the
    /// records of the windows sharing the pair, indexed by the file's spelling
    /// of the chromosome.
    reads: HashMap<(usize, u64), BatchQueryResult>,
    /// The file's name for each input chromosome of the batch, and its index
    /// in `reads`.
    names: HashMap<String, (String, usize)>,
}

/// The variant indices at which the batch's match chunks start, then the
/// batch's length: a chunk holds [`RECORDS_PER_CHUNK`] whole input records.
fn chunk_bounds(variants: &[InputVariant]) -> Vec<usize> {
    let mut bounds = vec![0];
    let mut records_in_chunk = 0;
    let mut index = 0;
    for group in variants.chunk_by(InputVariant::same_record) {
        if records_in_chunk == RECORDS_PER_CHUNK {
            bounds.push(index);
            records_in_chunk = 0;
        }
        index += group.len();
        records_in_chunk += 1;
    }
    bounds.push(variants.len());
    bounds
}

/// The read window of an allele's record: Ensembl VEP seeks each feature's
/// `[start - 1, end + 1]` (`File.pm:321`).
fn window(variant: &InputVariant) -> (u64, u64) {
    (
        variant.start.saturating_sub(1),
        variant.end.saturating_add(1),
    )
}

/// `number` as a coordinate, held at the upper end of the 64-bit range.
pub(super) fn whole_number(number: u64) -> Coordinate {
    Coordinate::Whole(i64::try_from(number).unwrap_or(i64::MAX))
}

/// The number of the cell at `column` of `record`, which the reader read to
/// the leading integer `leading`: that integer for a cell that is a plain
/// whole number, else the number Perl reads from the cell
/// ([`Coordinate::read`]), the leading integer again where Perl reads an
/// infinity or no number. A cell that is its number's decimal digits and
/// nothing else has their count for its length, so a plain cell is told by
/// its length alone and is not read a second time.
fn cell_number(record: &TabixRecord, column: usize, leading: u64) -> Coordinate {
    let decimal_digits = leading.checked_ilog10().map_or(1, |log| log as usize + 1);
    match record.columns.get(column) {
        Some(cell) if cell.len() != decimal_digits => {
            Coordinate::read(cell).unwrap_or_else(|| whole_number(leading))
        }
        _ => whole_number(leading),
    }
}

impl TabixSource {
    /// Opens the file of `spec` on the run's `assembly`.
    ///
    /// The file's name for a record's chromosome is Ensembl VEP's
    /// (`source_chr_name`); a name the file lacks is then tried under
    /// its other spellings by the reader. On GRCh37, where `MT` and `chrMT`
    /// name one mitochondrial sequence and `chrM` and `M` another, the reader
    /// tries such a name inside its own pair alone, so a contig of the other
    /// pair is read only under a name `source_chr_name` gives the record: a
    /// synonym the file lists, or `MT` for a record written bare `M`
    /// (`Parser.pm:794-796`).
    pub fn open(
        spec: CustomSpec,
        synonyms: Arc<ChromosomeSynonyms>,
        assembly: Option<vep_core::assembly::Assembly>,
    ) -> Result<Self, CustomError> {
        let open_error = |reason: String| CustomError::Open {
            file: spec.file.clone(),
            reason,
        };
        let mut annotator = TabixAnnotator::open(tabix_config(&spec))
            .map_err(|e| open_error(tabix_reason(e)))?
            .reading_rows_as_tabix_indexes_them()
            .on_assembly(assembly);
        annotator.frame_mut().keep_to_own_pair();
        let sequence_names = index_sequence_names(&spec)?;
        Ok(TabixSource {
            short_name: spec.short_name.as_str().into(),
            spec,
            annotator,
            sequence_names,
            synonyms,
        })
    }

    pub fn spec(&self) -> &CustomSpec {
        &self.spec
    }

    /// The start of `record` as the formats name and match it: the number
    /// Perl reads from the start cell, which the parser of ensembl-io hands
    /// Ensembl VEP as written (`GXF.pm:163-166`, `BaseVCF4.pm:244-247`), one
    /// more for a zero-based start (`Bed.pm:159-162`). Ensembl VEP walks to
    /// the record, names it and matches it by that number
    /// (`File.pm:333-335, 489-495, 523`), while the index returned the row by
    /// the leading integer `tabix` read ([`TabixConfig::record_start`]), the
    /// first base for a start a minus sign leads (htslib 1.9 `tbx.c:99`).
    pub fn start_of(&self, record: &TabixRecord) -> Coordinate {
        let config = self.annotator.config();
        let shift = u64::from(config.zero_based);
        cell_number(record, config.start_col, record.start - shift)
            .plus(i64::from(config.zero_based))
    }

    /// The end of `record` as the formats name and match it: the number Perl
    /// reads from the generic end cell (`GXF.pm:187-190`, `Bed.pm:197-200`),
    /// where the index returned the row by the cell's leading integer
    /// ([`TabixConfig::record_end`]).
    pub fn end_of(&self, record: &TabixRecord) -> Coordinate {
        match self.annotator.config().end_col {
            Some(end_col) => cell_number(record, end_col, record.end),
            None => whole_number(record.end),
        }
    }

    /// `self` reading its file through the reader `configure` returns for the
    /// one it opened, for a format that reads a cell otherwise than as text.
    pub fn with_reader(mut self, configure: impl FnOnce(TabixAnnotator) -> TabixAnnotator) -> Self {
        self.annotator = configure(self.annotator);
        self
    }

    /// The index's sequence names ([`CustomSource::chromosomes`]).
    pub fn chromosomes(&self) -> Vec<String> {
        self.sequence_names.iter().cloned().collect()
    }

    /// The file's name for the input chromosome `chr` (`BaseVEP.pm:689-738`).
    fn source_chr_name(&self, chr: &str) -> String {
        self.synonyms
            .source_chr_name(chr, |name| self.sequence_names.contains(name))
    }

    /// The records every input record of the batch can reach: Ensembl VEP
    /// seeks each feature's window in turn (`File.pm:321`); the same windows,
    /// grouped by chromosome and by the index's 16 kb linear window and merged
    /// within each group, are read in parallel on the current pool, every
    /// linear window inflated once, and looked up per record in
    /// [`Self::annotate`].
    pub fn prefetch(&self, variants: &[InputVariant]) -> Result<PrefetchData, CustomError> {
        let mut names: HashMap<String, (String, usize)> = HashMap::new();
        let mut by_window: IndexMap<(usize, u64), Vec<(u64, u64)>> = IndexMap::new();
        for variant in variants.iter().filter(|v| !skipped(v)) {
            let chr_index = match names.get(&variant.original_chr) {
                Some((_, chr_index)) => *chr_index,
                None => {
                    let chr_index = names.len();
                    let file_name = self.source_chr_name(&variant.original_chr);
                    names.insert(variant.original_chr.clone(), (file_name, chr_index));
                    chr_index
                }
            };
            let (start, end) = window(variant);
            by_window
                .entry((chr_index, start >> LINEAR_WINDOW_SHIFT))
                .or_default()
                .push((start, end));
        }
        let file_names: HashMap<usize, &str> = names
            .values()
            .map(|(name, index)| (*index, name.as_str()))
            .collect();
        let reads = by_window
            .iter()
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|(&(chr_index, linear_window), windows)| {
                let chr = file_names[&chr_index];
                let regions: Vec<(&str, u64, u64)> = windows
                    .iter()
                    .map(|&(start, end)| (chr, start, end))
                    .collect();
                let records =
                    self.annotator
                        .query_batch(&regions)
                        .map_err(|e| CustomError::Read {
                            file: self.spec.file.clone(),
                            reason: tabix_reason(e),
                        })?;
                Ok(((chr_index, linear_window), records))
            })
            .collect::<Result<HashMap<_, _>, _>>()?;
        Ok(PrefetchData::new(Prefetched {
            bounds: chunk_bounds(variants),
            reads,
            names,
        }))
    }

    /// Walks the batch input record by input record (`annotate_InputBuffer`,
    /// `File.pm:316-353`), the chunks of [`Self::prefetch`] in parallel on the
    /// current pool: each record's list and statistics (`record_annotation`)
    /// go to every allele of the record, as do the fields its `COSV` records
    /// lend to the namesakes earlier sources listed under the short name.
    /// The warnings the chunks raise join the hold of the calling thread in
    /// record order (`vep_core::warning_hold`).
    pub fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
        rules: &dyn FormatRules,
    ) -> Result<(), CustomError> {
        let data: Prefetched = data.downcast().map_err(|e| CustomError::Read {
            file: self.spec.file.clone(),
            reason: tabix_reason(e),
        })?;
        let mut slices: Vec<&mut [InputVariant]> = Vec::with_capacity(data.bounds.len());
        let mut rest = variants;
        for bound in data.bounds.windows(2) {
            let (chunk, tail) = rest.split_at_mut(bound[1] - bound[0]);
            slices.push(chunk);
            rest = tail;
        }
        map_in_order(slices, |chunk| {
            for group in chunk.chunk_by_mut(InputVariant::same_record) {
                if let Some(listed) = self.record_annotation(group, &data, rules) {
                    for variant in group.iter_mut() {
                        listed.leave_on(&self.short_name, &mut variant.custom_annotations);
                    }
                }
            }
        });
        Ok(())
    }

    /// What the source leaves on the input record `group`: the prefetched
    /// records reaching the feature's `[start - 1, end + 1]` are taken in
    /// file order while their start lies at or before the feature's end plus
    /// one and the cap allows, each read by `rules`; the entry closes as
    /// [`MatchedRecords::finish`] does. `None` for a record the sources pass
    /// over.
    fn record_annotation(
        &self,
        group: &[InputVariant],
        data: &Prefetched,
        rules: &dyn FormatRules,
    ) -> Option<Listed> {
        let feature = RecordGroup::of(group)?;
        let (chr, chr_index) = data.names.get(&feature.lead.original_chr)?;
        let (start, end) = window(feature.lead);
        let read = data
            .reads
            .get(&(*chr_index, start >> LINEAR_WINDOW_SHIFT))?;
        let candidates = read.lookup(chr, start, end);
        let mut matched = MatchedRecords::new(
            &self.spec,
            &self.short_name,
            &feature.lead.custom_annotations,
        );
        for record in candidates {
            if matched.stops_before_next() {
                break;
            }
            match rules.match_record(record, &feature, matched.scores()) {
                Verdict::PastEnd => break,
                Verdict::NoMatch => {}
                Verdict::Match(records, percentage) => matched.add(records, percentage),
            }
        }
        Some(matched.finish())
    }
}

/// The sequence names of the tabix index beside the file of `spec`, the
/// source's valid chromosomes (`File.pm:369-372`).
fn index_sequence_names(spec: &CustomSpec) -> Result<HashSet<String>, CustomError> {
    let index = tabix::find_index_path(Path::new(&spec.file))
        .and_then(|index_path| tabix::load_index(&index_path))
        .map_err(|e| CustomError::Open {
            file: spec.file.clone(),
            reason: tabix_reason(e),
        })?;
    Ok(index
        .header()
        .map(|header| {
            header
                .reference_sequence_names()
                .iter()
                .map(|name| name.to_string())
                .collect()
        })
        .unwrap_or_default())
}

/// Opens the file of `spec` in its format and returns the source reading it;
/// the warning its `num_records` earned is written once the file is open, so
/// a run that refuses the source, or opens none, leaves no warnings file.
fn open_source(
    spec: CustomSpec,
    synonyms: &Arc<ChromosomeSynonyms>,
    options: CustomRunOptions,
) -> Result<Box<dyn CustomSource>, CustomError> {
    let open_error = |reason: String| CustomError::Open {
        file: spec.file.clone(),
        reason,
    };
    let source: Box<dyn CustomSource> = match spec.format {
        CustomFormat::Vcf => Box::new(super::vcf::VcfSource::open(
            spec,
            Arc::clone(synonyms),
            options,
        )?),
        CustomFormat::Bed => Box::new(super::bed::BedSource::open(
            spec,
            Arc::clone(synonyms),
            options.assembly,
        )?),
        CustomFormat::Gff if spec.gff_type == GffType::GencodePromoter => Box::new(
            GffOverlaySource::open(spec, Arc::clone(synonyms), options.assembly)?,
        ),
        CustomFormat::Gff | CustomFormat::Gtf => {
            TabixAnnotator::open(tabix_config(&spec)).map_err(|e| open_error(tabix_reason(e)))?;
            let sequence_names = index_sequence_names(&spec)?;
            Box::new(UnreadSource {
                spec,
                sequence_names,
            })
        }
        CustomFormat::BigWig => Box::new(BigWigSource::open(spec, Arc::clone(synonyms))?),
    };
    if let Some(warning) = &source.spec().num_records_warning {
        warn!("{warning}");
    }
    Ok(source)
}

/// One batch's prefetch results, one per source in the order given
/// ([`CustomSources::prefetch`]), with the warnings the prefetches raised, in
/// that order, which [`CustomSources::annotate`] raises when the batch is
/// annotated (`vep_core::warning_hold`).
pub struct CustomPrefetch(Vec<Result<PrefetchData, CustomError>>, Vec<Held>);

/// The `--custom` sources of one run, in the order given.
pub struct CustomSources {
    sources: Vec<Box<dyn CustomSource>>,
}

impl fmt::Debug for CustomSources {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list()
            .entries(self.sources.iter().map(|source| source.spec()))
            .finish()
    }
}

impl CustomSources {
    /// Opens every specification's file; the first that cannot be opened ends
    /// the run before any input is read. `synonyms` is the run's chromosome
    /// synonym table, through which each source maps an input's chromosome
    /// name to its own (`BaseVEP.pm:689-738`); `options` are the run-wide
    /// flags a source reads.
    pub fn open(
        specs: Vec<CustomSpec>,
        synonyms: Arc<ChromosomeSynonyms>,
        options: CustomRunOptions,
    ) -> Result<Self, CustomError> {
        let sources = specs
            .into_iter()
            .map(|spec| open_source(spec, &synonyms, options))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(CustomSources { sources })
    }

    pub fn len(&self) -> usize {
        self.sources.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sources.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &dyn CustomSource> {
        self.sources.iter().map(AsRef::as_ref)
    }

    /// Annotates a batch from every source: the prefetches in parallel on the
    /// current rayon pool, then each source's annotate phase in order.
    pub fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), CustomError> {
        self.annotate(variants, self.prefetch(variants))
    }

    /// Every source's prefetch over a batch, in parallel on the current rayon
    /// pool. The sources read the records' coordinates and the parse's marks
    /// alone, so a batch may be prefetched before it is annotated, and the
    /// result is positional: it is applied to the batch it was read for, in
    /// the same order and length. A warning a prefetch raises is kept with
    /// the result, so a batch read ahead and never annotated leaves none.
    pub fn prefetch(&self, variants: &[InputVariant]) -> CustomPrefetch {
        let (read, raised) = collect_in_order(&self.sources, |source| source.prefetch(variants));
        CustomPrefetch(read, raised)
    }

    /// Each source's annotate phase in order over its prefetched records; the
    /// first source, in that order, whose prefetch or annotate failed is the
    /// error. The warnings `prefetch` carries join the hold of the calling
    /// thread ahead of those the annotate phases raise.
    pub fn annotate(
        &self,
        variants: &mut [InputVariant],
        prefetch: CustomPrefetch,
    ) -> Result<(), CustomError> {
        let CustomPrefetch(read, raised_ahead) = prefetch;
        debug_assert_eq!(read.len(), self.sources.len());
        raise_again(raised_ahead);
        for (source, data) in self.sources.iter().zip(read) {
            source.annotate(variants, data?)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::custom::grammar::parse_custom;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec, Warnings};

    fn spec(argument: &str) -> CustomSpec {
        parse_custom(argument).unwrap()
    }

    fn open(specs: Vec<CustomSpec>) -> Result<CustomSources, CustomError> {
        CustomSources::open(specs, Arc::default(), CustomRunOptions::default())
    }

    fn record(name: &str, score: Option<CustomValue>) -> CustomRecord {
        CustomRecord::new(name.to_string(), None, IndexMap::new(), score)
    }

    /// A file that is not there, or has no index beside it, is refused with the
    /// file named; Ensembl VEP opens the file for its chromosome names while
    /// the input parser is built, before the output opens (`BaseRunner.pm:453`,
    /// `File.pm:371`, `File/VCF.pm:154`, ensembl-io `TabixParser.pm:52-63`).
    #[test]
    fn open_refuses_a_missing_or_unindexed_file() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("absent.vcf.gz");
        let spec = parse_custom(&format!("file={},format=vcf", missing.display())).unwrap();
        let error = open(vec![spec]).unwrap_err();
        assert!(
            matches!(&error, CustomError::Open { file, .. } if *file == missing.display().to_string()),
            "{error}"
        );
        assert!(
            error.to_string().starts_with(&format!(
                "cannot open the --custom source {}: ",
                missing.display()
            )),
            "{error}"
        );

        let unindexed = dir.path().join("plain.bed.gz");
        std::fs::write(&unindexed, b"").unwrap();
        let spec = parse_custom(&format!("{},regions,bed", unindexed.display())).unwrap();
        let error = open(vec![spec]).unwrap_err();
        assert!(error.to_string().contains("index"), "{error}");

        let not_bigwig = dir.path().join("plain.bw");
        std::fs::write(&not_bigwig, b"not a bigWig").unwrap();
        let spec = parse_custom(&format!("{},cons,bigwig,exact", not_bigwig.display())).unwrap();
        assert!(matches!(
            open(vec![spec]).unwrap_err(),
            CustomError::Open { .. }
        ));
    }

    /// A bgzipped, indexed source opens by its content whatever its name: one
    /// named `sites_plain.vcf` and one named `regions.bed.GZ` annotate the
    /// variant their record overlaps, and a plain-text file with an index
    /// beside it is refused by its bytes. Ensembl VEP 116.3 tests the name
    /// alone before opening any tabix source (ensembl-io `TabixParser.pm:52`,
    /// `/\.b?gz$/`), so it reads neither of the two and ends the run, the
    /// documented class `tabix_source_named_without_gz`.
    #[test]
    fn a_bgzipped_source_opens_under_a_name_without_gz() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = write_tabix_fixture(
            dir.path(),
            "sites_plain.vcf",
            &[
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "21\t100\tmulti_1\tA\tG\t.\tPASS\tAF=0.003",
            ],
            IndexSpec::standard(),
        );
        let bed = write_tabix_fixture(
            dir.path(),
            "regions.bed.GZ",
            &["21\t90\t110\texact_1\t7"],
            IndexSpec::standard(),
        );
        let sources = open(vec![
            spec(&format!(
                "file={},short_name=sites,format=vcf,type=exact,fields=AF",
                vcf.display()
            )),
            spec(&format!("{},regions,bed,overlap", bed.display())),
        ])
        .unwrap();
        let mut batch = vec![InputVariant::new(
            "21".into(),
            100,
            100,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        sources.run_batch(&mut batch).unwrap();
        let listed: Vec<(&str, Vec<&str>)> = batch[0]
            .custom_annotations
            .iter()
            .map(|a| {
                (
                    a.short_name.as_ref(),
                    a.records.iter().map(|r| r.name.as_str()).collect(),
                )
            })
            .collect();
        assert_eq!(
            listed,
            [("sites", vec!["multi_1"]), ("regions", vec!["exact_1"])]
        );
        assert_eq!(
            batch[0].custom_annotations[0].records[0].fields.get("AF"),
            Some(&CustomFieldValue::Text("0.003".to_string()))
        );

        // The index is read before the data file, so a valid one sits beside
        // the plain text and the refusal is the data file's, not the index's.
        let lines = [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "21\t100\tplain_1\tA\tG\t.\tPASS\t.",
        ];
        let plain = write_tabix_fixture(dir.path(), "plain.vcf", &lines, IndexSpec::standard());
        std::fs::write(&plain, lines.join("\n") + "\n").unwrap();
        let error = open(vec![spec(&format!("file={},format=vcf", plain.display()))]).unwrap_err();
        let CustomError::Open { file, reason } = &error else {
            panic!("{error}");
        };
        assert_eq!(*file, plain.display().to_string());
        assert!(!reason.contains("index"), "{error}");
    }

    /// The warning a fractional or negative `num_records` earns (Ensembl VEP
    /// 116.3 compares the text as a number, `File.pm:337`, `:436-439`) is
    /// written when the source opens, once per source, and not for a source
    /// that cannot be opened, so a refused run leaves no warnings file.
    #[test]
    fn a_num_records_warning_is_written_once_the_source_opens() {
        let dir = tempfile::tempdir().unwrap();
        let bed = write_tabix_fixture(
            dir.path(),
            "regions.bed.gz",
            &["21\t90\t110\tregion_1\t7"],
            IndexSpec::standard().zero_based(),
        );
        let missing = dir.path().join("absent.bed.gz");
        let warnings = Arc::new(Warnings::default());
        let (opened, refused) = tracing::subscriber::with_default(Arc::clone(&warnings), || {
            let opened = open(vec![
                spec(&format!(
                    "file={},format=bed,num_records=2.5",
                    bed.display()
                )),
                spec(&format!("file={},format=bed,num_records=3", bed.display())),
            ]);
            let refused = open(vec![spec(&format!(
                "file={},format=bed,num_records=-1",
                missing.display()
            ))]);
            (opened, refused)
        });
        assert_eq!(opened.unwrap().len(), 2);
        assert!(matches!(refused.unwrap_err(), CustomError::Open { .. }));
        assert_eq!(
            *warnings.0.lock().unwrap(),
            [format!(
                "num_records=2.5 for {} is read as 2: the limit is a whole number of at least 0",
                bed.display()
            )]
        );
    }

    /// Indexed VCF, BED and GFF files open under their short names and annotate
    /// in the order given: the VCF, BED and promoter sources annotate the
    /// variant their record overlaps, while a batch passes unchanged through a
    /// GFF source in a mode Ensembl VEP does not know (`GFF.pm:155-165`
    /// annotate nothing from one).
    #[test]
    fn opened_sources_annotate_in_order_and_an_unknown_gff_mode_leaves_the_batch_as_found() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = write_tabix_fixture(
            dir.path(),
            "sites.vcf.gz",
            &[
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "21\t100\trs1\tA\tG\t.\tPASS\tAF=0.5",
            ],
            IndexSpec::standard(),
        );
        let bed = write_tabix_fixture(
            dir.path(),
            "regions.bed.gz",
            &["21\t90\t110\tregion_1\t7"],
            IndexSpec::standard(),
        );
        let gff = write_tabix_fixture(
            dir.path(),
            "promoters.gff.gz",
            &["21\t.\tpromoter_window\t90\t110\t.\t+\t.\tID=PW1"],
            IndexSpec::standard().with_span_cols(1, 4, 5),
        );
        let specs = vec![
            parse_custom(&format!(
                "file={},format=vcf,type=exact,fields=AF",
                vcf.display()
            ))
            .unwrap(),
            parse_custom(&format!("{},regions,bed", bed.display())).unwrap(),
            parse_custom(&format!(
                "file={},format=gff,gff_type=gencode_promoter",
                gff.display()
            ))
            .unwrap(),
        ];
        let sources = open(specs).unwrap();
        assert_eq!(sources.len(), 3);
        let names: Vec<&str> = sources
            .iter()
            .map(|s| s.spec().short_name.as_str())
            .collect();
        assert_eq!(names, ["sites.vcf.gz", "regions", "GENCODE_PROMOTER"]);

        let before = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let mut batch = vec![before.clone()];
        sources.run_batch(&mut batch).unwrap();
        let listed: Vec<(&str, Vec<&str>)> = batch[0]
            .custom_annotations
            .iter()
            .map(|a| {
                (
                    a.short_name.as_ref(),
                    a.records.iter().map(|r| r.name.as_str()).collect(),
                )
            })
            .collect();
        assert_eq!(
            listed,
            [
                ("sites.vcf.gz", vec!["rs1"]),
                ("regions", vec!["region_1"]),
                ("GENCODE_PROMOTER", vec!["21:90-110"]),
            ]
        );
        // The prefetch read before the batch is annotated and applied after
        // leaves the batch as `run_batch` does.
        let mut ahead = vec![before.clone()];
        let prefetch = sources.prefetch(&ahead);
        assert_eq!(prefetch.0.len(), 3);
        assert!(ahead[0].custom_annotations.is_empty());
        sources.annotate(&mut ahead, prefetch).unwrap();
        assert_eq!(format!("{ahead:?}"), format!("{batch:?}"));
        let mut unannotated = before.clone();
        open(vec![parse_custom(&format!(
            "file={},format=gff,gff_type=enhancer_window",
            gff.display()
        ))
        .unwrap()])
        .unwrap()
        .run_batch(std::slice::from_mut(&mut unannotated))
        .unwrap();
        assert_eq!(format!("{unannotated:?}"), format!("{before:?}"));
        assert!(open(Vec::new()).unwrap().is_empty());
    }

    /// A source that raises one warning in its prefetch, which runs the longer
    /// the earlier the source was given, and one in its annotate, about the
    /// batch's first record: each as the layer that writes a run's warnings
    /// gives it to the hold of the thread raising it.
    struct WarningSource {
        spec: CustomSpec,
        name: &'static str,
        prefetch_ms: u64,
    }

    impl CustomSource for WarningSource {
        fn spec(&self) -> &CustomSpec {
            &self.spec
        }

        fn chromosomes(&self) -> Vec<String> {
            Vec::new()
        }

        fn prefetch(&self, _variants: &[InputVariant]) -> Result<PrefetchData, CustomError> {
            std::thread::sleep(std::time::Duration::from_millis(self.prefetch_ms));
            let _ = vep_core::warning_hold::offer(0, format!("{} prefetch", self.name));
            Ok(PrefetchData::new(()))
        }

        fn annotate(
            &self,
            variants: &mut [InputVariant],
            _data: PrefetchData,
        ) -> Result<(), CustomError> {
            let _ = vep_core::warning_hold::offer(
                variants[0].input_record,
                format!("{} annotate", self.name),
            );
            Ok(())
        }
    }

    /// The warnings the sources raise over a batch join the hold of the thread
    /// that runs the batch in a fixed order, whatever thread ran which
    /// prefetch and in whatever order they ended: the prefetches' in the order
    /// the sources were given, then each source's annotate in that order,
    /// under the record it names. A prefetch read ahead keeps its warnings
    /// from the hold of the thread that reads ahead, and the annotate phase
    /// that applies it raises them ahead of its own; one that is dropped
    /// leaves none.
    #[test]
    fn source_warnings_join_the_callers_hold_in_the_order_given() {
        use vep_core::warning_hold::holding;
        const NAMES: [&str; 3] = ["first", "second", "third"];
        let sources = CustomSources {
            sources: NAMES
                .iter()
                .enumerate()
                .map(|(index, name)| {
                    Box::new(WarningSource {
                        spec: spec("file=made_up.bed.gz,short_name=made_up,format=bed"),
                        name,
                        prefetch_ms: 20 * (NAMES.len() - index) as u64,
                    }) as Box<dyn CustomSource>
                })
                .collect(),
        };
        let held = |input_record: u64, text: String| Held { input_record, text };
        let prefetches: Vec<Held> = NAMES
            .iter()
            .map(|name| held(0, format!("{name} prefetch")))
            .collect();
        let annotates: Vec<Held> = NAMES
            .iter()
            .map(|name| held(41, format!("{name} annotate")))
            .collect();
        let all = [prefetches, annotates].concat();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(NAMES.len())
            .build()
            .unwrap();
        let mut batch = vec![InputVariant::new(
            "21".into(),
            100,
            100,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        batch[0].input_record = 41;

        let (ran, raised) = pool.install(|| holding(|| sources.run_batch(&mut batch)));
        ran.unwrap();
        assert_eq!(raised, all);

        let (ahead, raised) = pool.install(|| holding(|| sources.prefetch(&batch)));
        assert!(raised.is_empty(), "{raised:?}");
        let (ran, raised) = pool.install(|| holding(|| sources.annotate(&mut batch, ahead)));
        ran.unwrap();
        assert_eq!(raised, all);

        let (dropped, raised) = pool.install(|| holding(|| sources.prefetch(&batch)));
        drop(dropped);
        assert!(raised.is_empty(), "{raised:?}");
    }

    /// Rules that match nothing and raise one warning for every record
    /// walked, about the input record walking it, as the layer that writes a
    /// run's warnings gives it to the hold of the thread raising it. A walk
    /// takes long enough for the chunks of a batch to run on several threads.
    struct WarningRules;

    impl FormatRules for WarningRules {
        fn match_record(
            &self,
            _record: &TabixRecord,
            group: &RecordGroup<'_>,
            _scores: bool,
        ) -> Verdict {
            std::thread::sleep(std::time::Duration::from_micros(200));
            let record = group.lead.input_record;
            let _ = vep_core::warning_hold::offer(record, format!("record {record}"));
            Verdict::NoMatch
        }
    }

    /// A batch of four chunks, matched on a pool of four threads, raises the
    /// warnings of its records into the hold of the thread that runs the
    /// batch in record order, each under its record, whatever thread matched
    /// which chunk.
    #[test]
    fn a_batch_matched_in_chunks_raises_its_warnings_in_record_order() {
        use vep_core::warning_hold::holding;
        let dir = tempfile::tempdir().unwrap();
        let bed = write_tabix_fixture(
            dir.path(),
            "whole.bed.gz",
            &["21\t0\t5000\twhole"],
            IndexSpec::standard(),
        );
        let source = TabixSource::open(
            spec(&format!("file={},format=bed", bed.display())),
            Arc::default(),
            None,
        )
        .unwrap();
        let records = 3 * RECORDS_PER_CHUNK as u64 + 10;
        let mut batch: Vec<InputVariant> = (1..=records)
            .map(|record| {
                let mut variant =
                    InputVariant::new("21".into(), record, record, b"A".to_vec(), b"G".to_vec());
                variant.input_record = record;
                variant
            })
            .collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let (matched, raised) = pool.install(|| {
            holding(|| {
                let data = source.prefetch(&batch)?;
                source.annotate(&mut batch, data, &WarningRules)
            })
        });
        matched.unwrap();
        let expected: Vec<Held> = (1..=records)
            .map(|record| Held {
                input_record: record,
                text: format!("record {record}"),
            })
            .collect();
        assert_eq!(raised, expected);
    }

    /// `File.pm:534-569` on a record spanning 100-200: `overlap` admits a
    /// feature sharing a base; `within` (`:541`) a feature holding the whole
    /// record; `surrounding` (`:544`) a feature the record holds; `distance`
    /// (`:546-549`) a feature whose two ends each lie within the distance of
    /// the record's; `overlap_cutoff` (`:566`) a feature covering at least
    /// that share of the record (51 of 101 bases is 50.495 percent).
    #[test]
    fn overlap_within_surrounding_distance_and_cutoff_follow_file_pm() {
        let features = [(150, 250), (120, 180), (90, 210)];
        let matched = |argument: &str| -> Vec<bool> {
            let spec = spec(argument);
            features
                .iter()
                .map(|&(start, end)| base_overlap(&spec, 100, 200, start, end).is_some())
                .collect()
        };
        assert_eq!(matched("file=r.bed.gz,format=bed"), [true, true, true]);
        assert_eq!(
            matched("file=r.bed.gz,format=bed,type=within"),
            [false, false, true]
        );
        assert_eq!(
            matched("file=r.bed.gz,format=bed,type=surrounding"),
            [false, true, false]
        );
        let near = spec("file=r.bed.gz,format=bed,distance=10");
        assert!(base_overlap(&near, 100, 200, 95, 205).is_some());
        assert!(base_overlap(&near, 100, 200, 95, 215).is_none());
        assert!(base_overlap(&near, 100, 200, 85, 205).is_none());
        let half = spec("file=r.bed.gz,format=bed,overlap_cutoff=50");
        assert_eq!(
            base_overlap(&half, 100, 200, 150, 250).as_deref(),
            Some("50.495")
        );
        assert!(base_overlap(&half, 100, 200, 151, 250).is_none());
        assert!(base_overlap(&spec("file=r.bed.gz,format=bed"), 100, 200, 201, 300).is_none());
        // An insertion's feature, start one past its end, spans its two
        // flanking bases (`:536-538`).
        assert_eq!(
            base_overlap(&spec("file=r.bed.gz,format=bed"), 100, 200, 201, 200).as_deref(),
            Some("0.990")
        );
        let exact = spec("file=r.bed.gz,format=bed,type=exact");
        assert_eq!(
            base_overlap(&exact, 100, 200, 100, 200).as_deref(),
            Some("100")
        );
        assert!(base_overlap(&exact, 100, 200, 100, 201).is_none());
    }

    /// `File.pm:546-549` on a record that starts near the lowest integer: the
    /// distance is tested on the whole difference between the two starts, as
    /// Perl's subtraction widens to hold it, so a feature at 5000 lies
    /// 2^63 + 4999 bases from a start of -(2^63 - 1): past a `distance` of
    /// 2^63, within one of 2^63 + 8192.
    #[test]
    fn distance_is_tested_on_a_difference_past_the_integer_range() {
        let matched = |distance: &str| {
            let spec = spec(&format!("file=r.bed.gz,format=bed,distance={distance}"));
            base_overlap(&spec, -i64::MAX, 10_000, 5_000, 5_000)
        };
        assert_eq!(matched("9223372036854775808"), None);
        assert_eq!(matched("9223372036854784000").as_deref(), Some("0.000"));
    }

    /// `File.pm:547-548` tests the difference of two whole ends against the
    /// distance as Perl reads the option's text. Digits alone and an exponent
    /// form with a whole value are an integer, and the test is exact: a start
    /// 2^53 + 1 bases from the feature's lies past a distance of 2^53. A point
    /// without an exponent or trailing text makes the distance a float, which
    /// the difference rounds to. A record with a float for one end has its
    /// whole end tested the same exact way, and its float end in floats,
    /// where an integer distance of 2^53 + 3 rounds to the 2^53 + 4 it admits.
    #[test]
    fn a_whole_distance_is_compared_exactly_past_two_to_the_fifty_third() {
        const TWO_53: i64 = 1 << 53;
        let record = TWO_53 + 2;
        let matched = |distance: &str, start: Coordinate, end: Coordinate, feature_start: i64| {
            let spec = spec(&format!("file=r.bed.gz,format=bed,distance={distance}"));
            base_overlap(&spec, start, end, feature_start, record).is_some()
        };
        let whole = |distance: &str| matched(distance, record.into(), record.into(), 1);
        assert!(!whole(&TWO_53.to_string()));
        assert!(whole(&(TWO_53 + 1).to_string()));
        assert!(!whole(&format!("+{TWO_53}")));
        assert!(!whole("9.007199254740992e15"));
        assert!(whole("9007199254740992.0"));
        assert!(whole("9007199254740992abc"));

        let float_end = Coordinate::Float(record as f64);
        let mixed = |distance: &str| matched(distance, record.into(), float_end, 1);
        assert!(!mixed(&TWO_53.to_string()));
        assert!(mixed(&(TWO_53 + 1).to_string()));

        let float_start = Coordinate::Float((TWO_53 + 4) as f64);
        let far_end = Coordinate::Whole(TWO_53 + 4);
        let floats = |distance: &str| {
            let spec = spec(&format!("file=r.bed.gz,format=bed,distance={distance}"));
            base_overlap(&spec, float_start, far_end, 0, TWO_53 + 4).is_some()
        };
        assert!(floats(&(TWO_53 + 3).to_string()));
        assert!(!floats(&(TWO_53 + 2).to_string()));
    }

    /// A record whose two ends both lie near the lowest integer, as a VCF
    /// insertion's do when its `POS` does (`File.pm:526-528` ends it at its
    /// start plus its reference's length less one), shares no base with a
    /// feature at 500 and matches nothing, under a cutoff or without one.
    #[test]
    fn a_record_ending_near_the_lowest_integer_shares_no_base_with_a_feature() {
        let near_lowest = -i64::MAX + 1;
        for options in ["", ",overlap_cutoff=50", ",reciprocal=1"] {
            let spec = spec(&format!("file=r.vcf.gz,format=vcf{options}"));
            assert_eq!(
                base_overlap(&spec, near_lowest, near_lowest, 500, 500),
                None,
                "{options:?}"
            );
        }
    }

    /// `File.pm:552-568`: the percentage is the share of the record the
    /// feature covers, or under `reciprocal` the share of the feature the
    /// record covers; the two agree for spans of one length (51 of 101 bases)
    /// and differ for a 50-base record inside a 101-base feature.
    #[test]
    fn reciprocal_zero_reports_the_record_share_and_one_the_variant_share() {
        let plain = spec("file=r.bed.gz,format=bed");
        let reciprocal = spec("file=r.bed.gz,format=bed,reciprocal=1");
        assert_eq!(
            base_overlap(&plain, 100, 200, 150, 250).as_deref(),
            Some("50.495")
        );
        assert_eq!(
            base_overlap(&reciprocal, 100, 200, 150, 250).as_deref(),
            Some("50.495")
        );
        assert_eq!(
            base_overlap(&plain, 150, 199, 150, 250).as_deref(),
            Some("100.000")
        );
        assert_eq!(
            base_overlap(&reciprocal, 150, 199, 150, 250).as_deref(),
            Some("49.505")
        );
    }

    /// `File.pm:436-440`: records are listed while fewer than `num_records`
    /// have matched, the record at the cap is listed as `...`, and the walk
    /// ends after it (`:337`); `num_records=0` lists the `...` alone. With a
    /// statistic asked for the walk never ends early (`:337`) and `count`
    /// counts every matching record.
    #[test]
    fn num_records_cap_appends_an_ellipsis_and_zero_lists_only_it() {
        let names = |argument: &str, matches: usize| -> (Vec<String>, u64, bool) {
            let spec = spec(argument);
            let short_name: Arc<str> = "r".into();
            let mut matched = MatchedRecords::new(&spec, &short_name, &[]);
            let mut stopped = false;
            for i in 1..=matches {
                if matched.stops_before_next() {
                    stopped = true;
                    break;
                }
                matched.add(
                    vec![record(
                        &format!("r{i}"),
                        Some(CustomValue::Text("1".into())),
                    )],
                    Some("100.000".into()),
                );
            }
            let count = matched.record_count;
            let listed = matched
                .finish()
                .entry
                .map(|a| a.records.into_iter().map(|r| r.name).collect())
                .unwrap_or_default();
            (listed, count, stopped)
        };
        assert_eq!(
            names("file=r.bed.gz,format=bed,num_records=2", 4),
            (vec!["r1".into(), "r2".into(), "...".into()], 3, true)
        );
        assert_eq!(
            names("file=r.bed.gz,format=bed,num_records=0", 3),
            (vec!["...".into()], 1, true)
        );
        assert_eq!(
            names(
                "file=r.bed.gz,format=bed,num_records=2,summary_stats=count",
                4
            ),
            (vec!["r1".into(), "r2".into(), "...".into()], 4, false)
        );
        assert_eq!(
            names("file=r.bed.gz,format=bed,num_records=all", 3),
            (vec!["r1".into(), "r2".into(), "r3".into()], 3, false)
        );
        let spec = spec("file=r.bed.gz,format=bed");
        let short_name: Arc<str> = "r".into();
        assert!(MatchedRecords::new(&spec, &short_name, &[])
            .finish()
            .entry
            .is_none());
    }

    /// `File.pm:405-415` against the list an earlier source left under the
    /// short name: the `COSV` namesake is lent, not listed, takes no `PC` and
    /// enters no statistic (`:417-440`), and counts (`:339`, `:442`); a list
    /// under another short name is not searched. [`Listed::leave_on`] lends
    /// the namesake the fields it lacks on the entries it is given and then
    /// appends the entry.
    #[test]
    fn a_cosv_namesake_an_earlier_source_listed_is_lent_and_counted() {
        let spec = spec("file=r.vcf.gz,format=vcf,fields=PC,summary_stats=count%max");
        let short_name: Arc<str> = "COSMIC".into();
        let other: Arc<str> = "OTHER".into();
        let listed_by = |name: &Arc<str>, fields: &[(&str, &str)]| {
            CustomAnnotations::new(
                Arc::clone(name),
                vec![{
                    let mut custom_record = record("COSV1", None);
                    custom_record.fields = fields
                        .iter()
                        .map(|(k, v)| (k.to_string(), CustomFieldValue::Text(v.to_string())))
                        .collect();
                    custom_record
                }],
                true,
                None,
            )
        };
        let score = |text: &str| Some(CustomValue::Text(text.to_string()));
        let mut lent_fields = IndexMap::new();
        lent_fields.insert("GENE".to_string(), CustomFieldValue::Text("Y".into()));
        lent_fields.insert("AA".to_string(), CustomFieldValue::Text("p.1".into()));

        let earlier = [
            listed_by(&other, &[]),
            listed_by(&short_name, &[("GENE", "X")]),
        ];
        let mut matched = MatchedRecords::new(&spec, &short_name, &earlier);
        matched.add(
            vec![{
                let mut custom_record = record("COSV1", score("5"));
                custom_record.fields = lent_fields.clone();
                custom_record
            }],
            Some("100.000".into()),
        );
        matched.add(vec![record("COSV2", score("3"))], Some("100.000".into()));
        let listed = matched.finish();
        assert_eq!(listed.lent.len(), 1);
        assert_eq!(listed.lent[0].name, "COSV1");
        assert_eq!(listed.lent[0].fields, lent_fields);
        let entry = listed.entry.as_ref().unwrap();
        assert_eq!(entry.records.len(), 1);
        assert_eq!(entry.records[0].name, "COSV2");
        assert_eq!(
            entry.stats,
            Some(CustomStats::Values(
                [
                    ("max".to_string(), CustomValue::Text("3".into())),
                    ("count".to_string(), CustomValue::Number(2.0)),
                ]
                .into_iter()
                .collect()
            ))
        );

        let mut entries = earlier.to_vec();
        listed.leave_on(&short_name, &mut entries);
        assert_eq!(entries.len(), 3);
        assert_eq!(entries[0], earlier[0]);
        let fields: Vec<(&str, &CustomFieldValue)> = entries[1].records[0]
            .fields
            .iter()
            .map(|(k, v)| (k.as_str(), v))
            .collect();
        assert_eq!(
            fields,
            [
                ("GENE", &CustomFieldValue::Text("X".into())),
                ("AA", &CustomFieldValue::Text("p.1".into())),
            ]
        );
        assert_eq!(entries[2], *entry);

        let earlier = [listed_by(&other, &[])];
        let mut matched = MatchedRecords::new(&spec, &short_name, &earlier);
        matched.add(vec![record("COSV1", None)], Some("100.000".into()));
        let listed = matched.finish();
        assert!(listed.lent.is_empty());
        assert_eq!(listed.entry.unwrap().records[0].name, "COSV1");

        // A source whose only match is lent lists nothing and still matched;
        // one that matched nothing leaves its statistics unmatched.
        let earlier = [listed_by(&short_name, &[])];
        let mut matched = MatchedRecords::new(&spec, &short_name, &earlier);
        matched.add(vec![record("COSV1", score("5"))], Some("100.000".into()));
        let entry = matched.finish().entry.unwrap();
        assert!(entry.records.is_empty() && entry.matched, "{entry:?}");
        assert_eq!(
            entry.stats,
            Some(CustomStats::Values(
                [("count".to_string(), CustomValue::Number(1.0))]
                    .into_iter()
                    .collect()
            ))
        );
        let entry = MatchedRecords::new(&spec, &short_name, &[])
            .finish()
            .entry
            .unwrap();
        assert!(entry.records.is_empty() && !entry.matched, "{entry:?}");
    }

    /// `File.pm:426-432, 344-352` over the scores 0 and 5: a least score of 0
    /// is kept (`!defined || <`), the greatest is 5, the mean is the sum over
    /// the count of matching records, and `min` and `max` are the scores as
    /// the records carry them. A record without a score counts but adds no
    /// value; `sum` asked for alone with no value is Ensembl VEP's undefined
    /// hash, and `min` alone an empty one.
    #[test]
    fn summary_stats_keep_a_zero_min_and_max() {
        let stats = |argument: &str, scores: &[Option<CustomValue>]| -> Option<CustomStats> {
            let spec = spec(argument);
            let short_name: Arc<str> = "r".into();
            let mut matched = MatchedRecords::new(&spec, &short_name, &[]);
            for (i, score) in scores.iter().enumerate() {
                matched.add(vec![record(&format!("r{i}"), score.clone())], None);
            }
            matched.finish().entry.unwrap().stats
        };
        let text = |s: &str| Some(CustomValue::Text(s.to_string()));
        let all = stats(
            "file=r.bed.gz,format=bed,summary_stats=min%max%mean%count%sum",
            &[text("0"), text("5")],
        );
        let expected: IndexMap<String, CustomValue> = [
            ("min", CustomValue::Text("0".into())),
            ("max", CustomValue::Text("5".into())),
            ("mean", CustomValue::Number(2.5)),
            ("count", CustomValue::Number(2.0)),
            ("sum", CustomValue::Number(5.0)),
        ]
        .into_iter()
        .map(|(k, v)| (k.to_string(), v))
        .collect();
        assert_eq!(all, Some(CustomStats::Values(expected)));
        assert_eq!(
            stats(
                "file=r.bed.gz,format=bed,summary_stats=count%mean",
                &[text("3"), None, Some(CustomValue::Undefined)]
            ),
            Some(CustomStats::Values(
                [
                    ("mean".to_string(), CustomValue::Number(1.0)),
                    ("count".to_string(), CustomValue::Number(3.0)),
                ]
                .into_iter()
                .collect()
            ))
        );
        assert_eq!(
            stats("file=r.bed.gz,format=bed,summary_stats=sum", &[None]),
            Some(CustomStats::Undefined)
        );
        assert_eq!(
            stats("file=r.bed.gz,format=bed,summary_stats=min", &[None]),
            Some(CustomStats::Values(IndexMap::new()))
        );
        assert_eq!(
            stats(
                "file=r.bed.gz,format=bed,summary_stats=max",
                &[text("1e1"), text("9")]
            ),
            Some(CustomStats::Values(
                [("max".to_string(), CustomValue::Text("1e1".into()))]
                    .into_iter()
                    .collect()
            ))
        );
    }

    /// [`TabixSource::prefetch`] reads a batch by chromosome and 16 kb linear
    /// window and [`TabixSource::annotate`] matches it in chunks of
    /// [`RECORDS_PER_CHUNK`] input records: every input record lists the
    /// records a walk of its own `[start - 1, end + 1]` lists, in file order,
    /// whichever window a record it reaches starts in and whichever chunk
    /// holds it, and every allele of a record carries the record's list.
    #[test]
    fn windows_and_chunks_list_what_a_walk_of_each_record_lists() {
        let dir = tempfile::tempdir().unwrap();
        let boundary = 1u64 << LINEAR_WINDOW_SHIFT;
        let dense = RECORDS_PER_CHUNK as u64 + 4;
        // Three-base records `r<p>` over `p..=p + 2`, then `s1` across the
        // window boundary, `s3` ending just before it and `s2` just past it.
        let mut lines: Vec<String> = (1..=dense)
            .map(|p| format!("21\t{}\t{}\tr{p}", p - 1, p + 2))
            .collect();
        lines.push(format!("21\t{}\t{}\ts1", boundary - 101, boundary + 100));
        lines.push(format!("21\t{}\t{}\ts3", boundary - 6, boundary - 1));
        lines.push(format!("21\t{}\t{}\ts2", boundary + 9, boundary + 20));
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        let bed = write_tabix_fixture(dir.path(), "walk.bed.gz", &lines, IndexSpec::standard());
        let sources = open(vec![spec(&format!("file={},format=bed", bed.display()))]).unwrap();

        let mut positions: Vec<u64> = (1..=dense).collect();
        positions.extend([boundary - 1, boundary, boundary + 15, boundary + 200]);
        let mut batch: Vec<InputVariant> = Vec::new();
        for (record, &position) in positions.iter().enumerate() {
            for alt in [b"G", b"T"] {
                let mut variant =
                    InputVariant::new("21".into(), position, position, b"A".to_vec(), alt.to_vec());
                variant.input_record = record as u64 + 1;
                batch.push(variant);
            }
        }
        sources.run_batch(&mut batch).unwrap();

        let expected = |position: u64| -> Vec<String> {
            if position <= dense {
                (position.saturating_sub(2).max(1)..=position)
                    .map(|q| format!("r{q}"))
                    .collect()
            } else if position == boundary - 1 {
                vec!["s1".into(), "s3".into()]
            } else if position == boundary {
                vec!["s1".into()]
            } else if position == boundary + 15 {
                vec!["s1".into(), "s2".into()]
            } else {
                Vec::new()
            }
        };
        for alleles in batch.chunks(2) {
            assert_eq!(alleles[0].custom_annotations, alleles[1].custom_annotations);
            let listed: Vec<String> = alleles[0]
                .custom_annotations
                .iter()
                .flat_map(|a| a.records.iter().map(|r| r.name.clone()))
                .collect();
            assert_eq!(
                listed,
                expected(alleles[0].start),
                "at {}",
                alleles[0].start
            );
        }
    }

    /// A row that starts below the first base, its start cell `0` or led by a
    /// minus sign, is listed among the rows written at base 1 in the file's
    /// order, as Ensembl VEP 116.3 lists it: `tabix` files such a row at the
    /// first base (htslib 1.9 `tbx.c:99`) and Ensembl VEP walks the rows as
    /// the index hands them out (`File.pm:333-345`).
    #[test]
    fn a_row_below_the_first_base_is_listed_in_the_files_order() {
        let dir = tempfile::tempdir().unwrap();
        let listed_from = |name: &str, below: &str| -> Vec<String> {
            let vcf = write_tabix_fixture(
                dir.path(),
                name,
                &[
                    "##fileformat=VCFv4.2",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                    "21\t1\tat_one\tA\tG\t.\t.\t.",
                    below,
                    "21\t1\tat_one_again\tA\tG\t.\t.\t.",
                ],
                IndexSpec::standard(),
            );
            let sources = open(vec![spec(&format!(
                "file={},format=vcf,type=overlap,coords=0",
                vcf.display()
            ))])
            .unwrap();
            let mut batch = vec![InputVariant::new(
                "21".into(),
                1,
                1,
                b"A".to_vec(),
                b"G".to_vec(),
            )];
            sources.run_batch(&mut batch).unwrap();
            batch[0]
                .custom_annotations
                .iter()
                .flat_map(|entry| entry.records.iter().map(|record| record.name.clone()))
                .collect()
        };
        assert_eq!(
            listed_from("minus.vcf.gz", "21\t-5\tbelow\tAAAAAAAAAA\tA\t.\t.\t."),
            ["at_one", "below", "at_one_again"]
        );
        assert_eq!(
            listed_from("zero.vcf.gz", "21\t0\tzero\tAA\tA\t.\t.\t."),
            ["at_one", "zero", "at_one_again"]
        );
    }

    /// `Coordinate::read` reads a cell as the pinned image's Perl reads the
    /// string: an exponent or a `+` gives the number it spells, a fraction
    /// stays, text after the number is ignored, a hexadecimal prefix reads 0
    /// and a zero-padded number in decimal; a run of digits, or a whole number
    /// an exponent spells, is exact while it fits 64 bits and a float past
    /// them, and a whole number a point alone spells is a float from 2^53. A
    /// cell Perl reads as an infinity or as no number at all has no
    /// coordinate.
    #[test]
    fn a_cell_reads_as_the_number_perl_reads_from_it() {
        use Coordinate::{Float, Whole};
        for (cell, number) in [
            ("100", Whole(100)),
            ("+100", Whole(100)),
            (" 100", Whole(100)),
            ("100x", Whole(100)),
            ("1e2", Whole(100)),
            ("+1e2", Whole(100)),
            ("-1e2", Whole(-100)),
            ("1E2", Whole(100)),
            ("1e+2", Whole(100)),
            ("1.e2", Whole(100)),
            ("1.5e2", Whole(150)),
            ("2000e-1", Whole(200)),
            ("1e18", Whole(1_000_000_000_000_000_000)),
            ("-1e17", Whole(-100_000_000_000_000_000)),
            ("1.5e17", Whole(150_000_000_000_000_000)),
            ("100000000000000000.5e0", Whole(100_000_000_000_000_000)),
            ("9007199254740991.0", Whole(9_007_199_254_740_991)),
            ("5.", Whole(5)),
            ("1e", Whole(1)),
            ("1e+", Whole(1)),
            ("0x1A", Whole(0)),
            ("010", Whole(10)),
            ("-0", Whole(0)),
            ("xyz", Whole(0)),
            ("", Whole(0)),
            ("1.25e1", Float(12.5)),
            ("15e-1", Float(1.5)),
            ("1e-2", Float(0.01)),
            ("44107419.5", Float(44107419.5)),
            ("-5.5", Float(-5.5)),
            (".5", Float(0.5)),
            ("9007199254740992.0", Float(9_007_199_254_740_992.0)),
            ("-100000000000000000.0", Float(-1e17)),
            ("9223372036854775807", Whole(i64::MAX)),
            ("-9223372036854775807", Whole(-i64::MAX)),
            ("-9223372036854775808", Whole(i64::MIN)),
            ("9223372036854775808", Float(9_223_372_036_854_775_808.0)),
            ("-9223372036854775809", Float(-9_223_372_036_854_775_808.0)),
            ("18446744073709551616", Float(18_446_744_073_709_551_616.0)),
            ("1e19", Float(1e19)),
            ("-1e30", Float(-1e30)),
        ] {
            assert_eq!(Coordinate::read(cell), Some(number), "{cell:?}");
        }
        for cell in ["inf", "-inf", "Infinity", "nan", "NaN", "1e999", "-1e999"] {
            assert_eq!(Coordinate::read(cell), None, "{cell:?}");
        }
    }

    /// The whole number of a coordinate is what `%i` prints for it on the
    /// pinned image (`File.pm:489-495`): a fraction dropped toward zero, a
    /// number below the 64-bit range held at its lower end. Past the upper end
    /// Perl prints a number wrapped below zero (`9223372036854775808` as
    /// `-9223372036854775808`, `1e19` as `-8446744073709551616`), where the
    /// whole number here holds at the upper end. One more is added as Perl
    /// adds it (`Bed.pm:161`): exactly inside the range, in the float past it
    /// and, from 2^53, to a whole number a point alone spells
    /// (`-100000000000000000.0` stays itself).
    #[test]
    fn a_whole_number_drops_the_fraction_and_holds_at_the_ends_of_the_range() {
        let whole = |cell: &str| Coordinate::read(cell).unwrap().whole();
        for (cell, number) in [
            ("1.25e1", 12),
            ("15e-1", 1),
            ("1e-2", 0),
            ("-1e-2", 0),
            ("44107419.5", 44107419),
            ("-5.5", -5),
            ("-9223372036854775809", i64::MIN),
            ("-1e19", i64::MIN),
            ("-1e30", i64::MIN),
            ("9223372036854775808", i64::MAX),
            ("18446744073709551615", i64::MAX),
            ("1e19", i64::MAX),
            ("1e30", i64::MAX),
        ] {
            assert_eq!(whole(cell), number, "{cell:?}");
        }
        let one_more = |cell: &str| Coordinate::read(cell).unwrap().plus(1).whole();
        for (cell, number) in [
            ("-5.5", -4),
            ("-0.5", 0),
            ("-1e17", -99_999_999_999_999_999),
            ("1.5e17", 150_000_000_000_000_001),
            ("100000000000000000.5e0", 100_000_000_000_000_001),
            ("-100000000000000000.0", -100_000_000_000_000_000),
            ("9007199254740992.0", 9_007_199_254_740_992),
            ("9007199254740993.0", 9_007_199_254_740_992),
            ("44107419.5", 44107420),
            ("-9223372036854775808", -i64::MAX),
            ("-9223372036854775809", i64::MIN),
            ("9223372036854775807", i64::MAX),
        ] {
            assert_eq!(one_more(cell), number, "{cell:?}");
        }
        assert!(Coordinate::Float(100.5).is_past(100));
        assert!(!Coordinate::Float(100.5).is_past(101));
        assert!(!Coordinate::Whole(100).is_past(100));
    }

    /// `File.pm:541-566` on a record one of whose ends has a fraction, each
    /// row as Ensembl VEP 116.3 lists it on the pinned image: a record from
    /// 100.5 to 102 does not surround a feature at 100 or one from 100 to
    /// 102, and a feature at 100 covers 20 percent of it and one at 101 40, so
    /// a cutoff of 35 keeps the second alone; a record from 100 to 102.5 does
    /// not lie within a feature from 100 to 102; neither lies at a distance of
    /// 0 from that feature. The compiled `overlap` reads whole numbers
    /// (`:569`), so both records share a base with a feature at 100 or 102.
    #[test]
    fn a_record_with_a_fraction_is_matched_on_the_number_as_written() {
        let from_fraction = (Coordinate::Float(100.5), Coordinate::Whole(102));
        let to_fraction = (Coordinate::Whole(100), Coordinate::Float(102.5));
        let overlap = |options: &str, record: (Coordinate, Coordinate), start: i64, end: i64| {
            let spec = spec(&format!("file=r.bed.gz,format=bed{options}"));
            base_overlap(&spec, record.0, record.1, start, end)
        };
        assert_eq!(
            overlap("", from_fraction, 100, 100).as_deref(),
            Some("20.000")
        );
        assert_eq!(
            overlap("", from_fraction, 101, 101).as_deref(),
            Some("40.000")
        );
        assert_eq!(
            overlap("", to_fraction, 102, 102).as_deref(),
            Some("28.571")
        );
        assert_eq!(overlap("", to_fraction, 103, 103), None);
        assert_eq!(overlap(",type=surrounding", from_fraction, 100, 100), None);
        assert_eq!(overlap(",type=surrounding", from_fraction, 100, 102), None);
        assert!(overlap(",type=surrounding", from_fraction, 101, 102).is_some());
        assert!(overlap(",type=surrounding", to_fraction, 100, 102).is_some());
        assert_eq!(overlap(",type=within", to_fraction, 100, 102), None);
        assert!(overlap(",type=within", to_fraction, 2, 150).is_some());
        assert!(overlap(",type=within", from_fraction, 100, 102).is_some());
        assert_eq!(overlap(",distance=0", from_fraction, 100, 102), None);
        assert_eq!(overlap(",distance=0", to_fraction, 100, 102), None);
        assert!(overlap(",distance=0.5", from_fraction, 100, 102).is_some());
        assert_eq!(overlap(",overlap_cutoff=35", from_fraction, 100, 100), None);
        assert!(overlap(",overlap_cutoff=35", from_fraction, 101, 101).is_some());
        let whole = (Coordinate::Whole(100), Coordinate::Whole(102));
        assert_eq!(overlap(",overlap_cutoff=35", whole, 101, 101), None);
        assert_eq!(overlap(",type=exact", from_fraction, 100, 102), None);
        // The percentages Ensembl VEP prints for a BED record from -5.5 to 10
        // and for one from 44107419.5, its start one more (`Bed.pm:161`).
        let below_zero = (Coordinate::Float(-4.5), Coordinate::Whole(10));
        assert_eq!(overlap("", below_zero, 1, 1).as_deref(), Some("6.452"));
        assert_eq!(overlap("", below_zero, 2, 150).as_deref(), Some("58.065"));
        let far = (Coordinate::Float(44107420.5), Coordinate::Whole(44107430));
        assert_eq!(
            overlap("", far, 44107424, 44107424).as_deref(),
            Some("9.524")
        );
    }

    /// The runs of `tests/data/custom_coordinate_corners.tsv`, one row each.
    const CORNERS: &str = include_str!("../../tests/data/custom_coordinate_corners.tsv");

    /// The input records every run of the corner table reads, with the name
    /// the table gives each: single bases of sequence 20 and of 21, and three
    /// deletions of sequence 20 as structural records.
    fn corner_probes() -> Vec<(String, InputVariant)> {
        let base = |chr: &str, position: u64| {
            InputVariant::new(chr.into(), position, position, b"A".to_vec(), b"G".to_vec())
        };
        let deletion = |start: u64, end: u64| {
            let mut variant =
                InputVariant::new("20".into(), start, end, b"N".to_vec(), b"<DEL>".to_vec());
            variant.is_structural = true;
            (format!("del {start}-{end}"), variant)
        };
        let on_20 = |position: u64| (position.to_string(), base("20", position));
        let mut probes = vec![on_20(1), deletion(2, 150)];
        probes.extend([2, 3, 12, 13, 55, 99].map(on_20));
        probes.push(deletion(100, 102));
        probes.extend([100, 101, 102, 103].map(on_20));
        probes.push(deletion(120, 300));
        probes.extend(
            [
                125, 149, 150, 151, 199, 200, 201, 250, 251, 300, 301, 405, 1999, 2000, 2001, 3000,
                3001, 5000, 5001, 44107424, 44107431,
            ]
            .map(on_20),
        );
        probes.extend(
            [1, 100, 125, 42395401, 44107424]
                .map(|position| (format!("21:{position}"), base("21", position))),
        );
        for (record, (_, variant)) in probes.iter_mut().enumerate() {
            variant.input_record = record as u64;
        }
        probes
    }

    /// One line of a corner file from its `sequence:start:end:name` row: a
    /// promoter window, a BED record, or a VCF row that is a deletion with an
    /// `END` or, for the end `SEQ250`, 250 reference bases.
    fn corner_line(route: &str, row: &str) -> String {
        let cells: Vec<&str> = row.split(':').collect();
        let (seq, start, end, name) = (cells[0], cells[1], cells[2], cells[3]);
        match route {
            "gff" => format!(
                "{seq}\t.\tpromoter_window\t{start}\t{end}\t.\t+\t.\tID={name}_{seq};gene_id=G{name}"
            ),
            "bed" => format!("{seq}\t{start}\t{end}\t{name}_{seq}"),
            _ if end == "SEQ250" => format!(
                "{seq}\t{start}\t{name}_{seq}\t{}\tC\t.\t.\tAF=0.1",
                "A".repeat(250)
            ),
            _ => format!("{seq}\t{start}\t{name}_{seq}\tN\t<DEL>\t.\t.\tSVTYPE=DEL;{end}"),
        }
    }

    /// What each probe lists, as the corner table writes a listing: the
    /// records' names and, for a source that prints one, their percentages,
    /// then the probes that list them, in the order the values first appear.
    fn corner_listing(route: &str, probes: &[(String, InputVariant)]) -> String {
        let mut groups: Vec<(String, Vec<&str>)> = Vec::new();
        for (name, variant) in probes {
            let records = variant
                .custom_annotations
                .first()
                .map_or(&[][..], |annotations| &annotations.records[..]);
            let names: Vec<&str> = records.iter().map(|record| record.name.as_str()).collect();
            let percentages: Vec<&str> = records
                .iter()
                .filter_map(|record| match record.fields.get("PC") {
                    Some(CustomFieldValue::Text(percentage)) if route != "gff" => {
                        Some(percentage.as_str())
                    }
                    _ => None,
                })
                .collect();
            let mut value = if names.is_empty() {
                "-".to_string()
            } else {
                names.join(",")
            };
            if !percentages.is_empty() {
                value = format!("{value} [{}]", percentages.join(","));
            }
            match groups.iter_mut().find(|(listed, _)| *listed == value) {
                Some((_, at)) => at.push(name),
                None => groups.push((value, vec![name])),
            }
        }
        groups
            .iter()
            .map(|(value, at)| format!("{value} at {}", at.join(",")))
            .collect::<Vec<_>>()
            .join(" | ")
    }

    /// The cell of a record's start or end column is read once more only
    /// when it is no plain whole number: a plain cell keeps the reader's own
    /// number, held at the upper end of the range from 2^63 up, where it
    /// would wrap below zero; any other cell is the number Perl reads from
    /// it, the reader's leading integer again for one Perl reads as an
    /// infinity; a row short of the column keeps the reader's number.
    #[test]
    fn a_plain_cell_keeps_the_readers_number_and_any_other_is_read_as_perl_reads_it() {
        use Coordinate::{Float, Whole};
        let number = |cell: &str, leading: u64| {
            let record = TabixRecord {
                chr: "20".to_string(),
                start: 1,
                end: 1,
                columns: vec![String::new(), cell.to_string()],
                header: None,
            };
            cell_number(&record, 1, leading)
        };
        for (cell, leading, read) in [
            ("100", 100, Whole(100)),
            ("0", 0, Whole(0)),
            ("9223372036854775807", i64::MAX as u64, Whole(i64::MAX)),
            ("9223372036854775808", 1 << 63, Whole(i64::MAX)),
            ("18446744073709551615", u64::MAX, Whole(i64::MAX)),
            ("007", 7, Whole(7)),
            ("+100", 100, Whole(100)),
            ("100x", 100, Whole(100)),
            ("1e2", 1, Whole(100)),
            ("2000e-1", 2000, Whole(200)),
            ("100.5", 100, Float(100.5)),
            ("-5", 0, Whole(-5)),
            (
                "-9223372036854775809",
                0,
                Float(-9_223_372_036_854_775_808.0),
            ),
            ("1e999", 1, Whole(1)),
        ] {
            assert_eq!(number(cell, leading), read, "{cell:?}");
        }
        let short = TabixRecord {
            chr: "20".to_string(),
            start: 100,
            end: 100,
            columns: vec![String::new(), "100".to_string()],
            header: None,
        };
        assert_eq!(cell_number(&short, 2, 100), Whole(100));
    }

    /// A record that starts below the 32-bit range is listed by the number
    /// written, here a deletion row whose `POS` is -9223372036854775807: for
    /// every record from the first base to its `END`, under its start of
    /// -9223372036854775806. Ensembl VEP 116.3 lists it from base 2 alone:
    /// the compiled `overlap` its image swaps in
    /// (`Utils/VariationEffect.pm:73-74`) takes each end as a 32-bit `int`
    /// (ensembl-xs 2.3.2 `Bio-EnsEMBL.xs:468-472`) and so reads the start as
    /// its low 32 bits, 2, where the Perl function of the same module
    /// compares the numbers.
    #[test]
    fn a_start_below_the_32_bit_range_is_listed_by_the_number_written() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = write_tabix_fixture(
            dir.path(),
            "low.vcf.gz",
            &[
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "20\t-9223372036854775807\tlow\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=100",
            ],
            IndexSpec::standard(),
        );
        let sources = open(vec![spec(&format!(
            "file={},format=vcf,type=overlap,coords=1",
            vcf.display()
        ))])
        .unwrap();
        let mut batch: Vec<InputVariant> = [1, 2, 3, 100, 101]
            .into_iter()
            .enumerate()
            .map(|(record, position)| {
                let mut variant = InputVariant::new(
                    "20".into(),
                    position,
                    position,
                    b"A".to_vec(),
                    b"G".to_vec(),
                );
                variant.input_record = record as u64;
                variant
            })
            .collect();
        sources.run_batch(&mut batch).unwrap();
        let listed: Vec<Vec<&str>> = batch
            .iter()
            .map(|variant| {
                variant
                    .custom_annotations
                    .iter()
                    .flat_map(|annotations| &annotations.records)
                    .map(|record| record.name.as_str())
                    .collect()
            })
            .collect();
        let low = vec!["20:-9223372036854775806-100"];
        assert_eq!(listed, [low.clone(), low.clone(), low.clone(), low, vec![]]);
    }

    /// Every run of the corner table lists what the table says. Each row is a
    /// `--custom` file of two sequences, one of whose lines spells its start
    /// or end as no plain whole number, read by 40 input records; its listing
    /// is the one Ensembl VEP 116.3 wrote on the pinned image over the same
    /// lines indexed by the image's `tabix`, for every row the table marks
    /// `matched`: an exponent, a fraction, a `+`, digits below the 64-bit
    /// range, in the start or the end of a GFF window, a BED record or a VCF
    /// row, under `overlap`, `within`, `surrounding`, `distance`, a cutoff and
    /// `reciprocal`. A row with another mark carries Ensembl VEP's listing
    /// beside the one held here: an end the index build fails, after which
    /// Ensembl VEP loses the window that follows for a record that listed an
    /// earlier one; an end past the 32 bits Ensembl VEP's compiled `overlap`
    /// reads; a start past the 64-bit range, which its index files at the
    /// sequence's first base; a VCF `END` the index reads below its `POS`;
    /// and a start Perl reads as an infinity.
    #[test]
    fn every_corner_run_lists_what_its_row_says() {
        let dir = tempfile::tempdir().unwrap();
        let mut differing = Vec::new();
        let mut matched = 0;
        for (number, row) in CORNERS.lines().enumerate().skip(1) {
            let cells: Vec<&str> = row.split('\t').collect();
            let (route, rows, options, relation, listing) =
                (cells[0], cells[3], cells[4], cells[5], cells[6]);
            let mut lines: Vec<String> = if route == "vcf" {
                vec![
                    "##fileformat=VCFv4.2".to_string(),
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
                ]
            } else {
                Vec::new()
            };
            lines.extend(rows.split(';').map(|row| corner_line(route, row)));
            let index = if route == "gff" {
                IndexSpec::standard().with_span_cols(1, 4, 5)
            } else {
                IndexSpec::standard()
            };
            let file = write_tabix_fixture(
                dir.path(),
                &format!("row{number}.{route}.gz"),
                &lines,
                index,
            );
            let sources = open(vec![spec(&format!("file={},{options}", file.display()))]).unwrap();
            let (names, mut batch): (Vec<String>, Vec<InputVariant>) =
                corner_probes().into_iter().unzip();
            sources.run_batch(&mut batch).unwrap();
            let probes: Vec<(String, InputVariant)> = names.into_iter().zip(batch).collect();
            let listed = corner_listing(route, &probes);
            if listed != listing {
                differing.push(format!(
                    "row {number} ({route} {} {:?}, {options}):\n  listed {listed}\n  table  {listing}",
                    cells[1], cells[2]
                ));
            }
            matched += usize::from(relation == "matched");
            // A row that is not matched carries the listing it differs from.
            assert_eq!(
                relation == "matched",
                cells.get(7).is_none_or(|ensembl| ensembl.is_empty())
            );
        }
        assert!(differing.is_empty(), "{}", differing.join("\n"));
        assert_eq!(matched, 78);
    }
}
