// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A bigWig `--custom` source, read through `bigtools` with Ensembl VEP's rules
//! for the format.
//!
//! A record is an interval of the file as the query for a variant returns it:
//! clipped to the variant's window, the bases from one before the variant's
//! start to one past its end (`AnnotationSource/File.pm:321`, kent
//! `bwgQuery.c`), its 0-based start read as a 1-based one
//! (`File.pm:529-532`). The record matches by `type` on those bounds
//! (`File.pm:534-575`), is named by its value or, under `coords`, by the
//! file's chromosome name with its 1-based start and end (the form of
//! `File/BigWig.pm:391-393`), enters the summary statistics (`File.pm:325-352,
//! 422-434`) and is listed up to `num_records` then as `...`
//! (`File.pm:436-440`). Ensembl VEP's per-record path names a record with the
//! chromosome name less its `chr` and the 0-based start (`BigWig.pm:416-428`
//! over ensembl-io `Parser/BigWig.pm:68, 107-110`), a mixed convention its
//! span path does not share and this source does not reproduce.
//!
//! A span of more than 65,536 bases under `type=overlap` with statistics of
//! `min` and `max` alone takes the span rule (`BigWig.pm:175-207, 239-249,
//! 264-296`): the extrema over the whole span, the records of the span's first
//! 32,768 bases, then `...` where the span continues (`BigWig.pm:369-401`).
//! Ensembl VEP annotates every shorter span twice on that path
//! (`BigWig.pm:195-204`); this source annotates each variant once.

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use bigtools::utils::file::reopen::ReopenableFile;
use bigtools::{BBIReadError, BigWigRead, Value};
use indexmap::IndexMap;
use vep_core::custom_annotation::{
    CustomAnnotations, CustomFieldValue, CustomRecord, CustomStats, CustomValue,
};
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::InputVariant;
use vep_core::variation::perl_number_string;

use super::grammar::{CustomSpec, MatchType, RecordLimit, SummaryStat};
use super::source::{CustomError, CustomSource};
use crate::bigwig::open_bigwig;
use crate::traits::PrefetchData;

/// The widest span the per-record walk covers; a longer one takes the span
/// rule (`BigWig.pm:86-88`: twice the exact end).
const SPAN_RULE_MIN_SPAN: u64 = 65_536;

/// The bases at the start of a long span whose records are listed
/// (`BigWig.pm:87`).
const SPAN_RECORD_WINDOW: u32 = 32_768;

/// The record Ensembl VEP lists in place of the ones past `num_records`.
const ELLIPSIS: &str = "...";

/// Two variants' windows this many bases apart or closer are read in one query:
/// a bigWig data block holds about a thousand items, so the bases between two
/// windows within one block cost less than a second walk of the file's index.
const CLUSTER_GAP: u32 = 1_000;

/// A chromosome of the file: its own name and its length.
#[derive(Debug, Clone)]
struct Chromosome {
    name: String,
    length: u32,
}

/// Which statistics the source computes (`summary_stats`).
#[derive(Debug, Clone, Copy, Default)]
struct Statistics {
    min: bool,
    max: bool,
    mean: bool,
    sum: bool,
    count: bool,
}

impl Statistics {
    fn from_spec(stats: &[SummaryStat]) -> Self {
        let mut wanted = Statistics::default();
        for stat in stats {
            match stat {
                SummaryStat::Min => wanted.min = true,
                SummaryStat::Max => wanted.max = true,
                SummaryStat::Mean => wanted.mean = true,
                SummaryStat::Sum => wanted.sum = true,
                SummaryStat::Count => wanted.count = true,
            }
        }
        wanted
    }

    /// Whether the statistics are a non-empty subset of `min` and `max`
    /// (`BigWig.pm:241`).
    fn extrema_only(self) -> bool {
        (self.min || self.max) && !(self.mean || self.sum || self.count)
    }
}

/// The running statistics of one variant (`File.pm:423-433`), closed into the
/// output map once its records are walked (`File.pm:344-351`).
#[derive(Debug, Default)]
struct Accumulator {
    min: Option<f64>,
    max: Option<f64>,
    sum: Option<f64>,
    count: u64,
}

impl Accumulator {
    fn add(&mut self, value: f64, wanted: Statistics) {
        if wanted.min && self.min.is_none_or(|min| value < min) {
            self.min = Some(value);
        }
        if wanted.max && self.max.is_none_or(|max| value > max) {
            self.max = Some(value);
        }
        if wanted.sum || wanted.mean {
            self.sum = Some(self.sum.unwrap_or(0.0) + value);
        }
    }

    /// The statistics in the order `summary_stats` named them: `mean` is the
    /// sum over the matched records divided by their count, and `sum` prints
    /// only when asked for. With no value added, Perl's closing (`File.pm:344-352`)
    /// creates the hash through `count`, the `mean` test or the deletion of an
    /// unasked `sum`, so `sum` asked for without `count` or `mean` leaves it
    /// undefined.
    fn close(self, order: &[SummaryStat]) -> CustomStats {
        let mut stats = IndexMap::new();
        for stat in order {
            let value = match stat {
                SummaryStat::Min => self.min,
                SummaryStat::Max => self.max,
                SummaryStat::Mean => self.sum.map(|sum| sum / self.count as f64),
                SummaryStat::Sum => self.sum,
                SummaryStat::Count => Some(self.count as f64),
            };
            if let Some(value) = value {
                stats.insert(stat.as_str().to_string(), CustomValue::Number(value));
            }
        }
        let creates_the_hash =
            |stat: &SummaryStat| matches!(stat, SummaryStat::Count | SummaryStat::Mean);
        if stats.is_empty()
            && order.contains(&SummaryStat::Sum)
            && !order.iter().any(creates_the_hash)
        {
            return CustomStats::Undefined;
        }
        CustomStats::Values(stats)
    }
}

/// How the source matches and lists records, read once from the specification.
#[derive(Debug, Clone)]
struct Rules {
    match_type: MatchType,
    report_coords: bool,
    overlap_cutoff: f64,
    distance: Option<f64>,
    reciprocal: bool,
    limit: RecordLimit,
    /// `summary_stats` as given, `None` when none were asked for.
    stats_order: Option<Vec<SummaryStat>>,
    stats: Statistics,
    /// `fields` names `PC`, so each record carries its overlap percentage.
    percentage_field: bool,
    /// Whether a long span takes the span rule (`BigWig.pm:239-249`).
    span_rule: bool,
}

impl Rules {
    fn from_spec(spec: &CustomSpec) -> Self {
        let stats_order = spec.summary_stats.clone();
        let stats = stats_order
            .as_deref()
            .map(Statistics::from_spec)
            .unwrap_or_default();
        let percentage_field = spec
            .fields
            .as_ref()
            .is_some_and(|fields| fields.iter().any(|field| field == "PC"));
        let span_rule = stats_order.is_some()
            && stats.extrema_only()
            && spec.match_type == MatchType::Overlap
            && spec.overlap_cutoff == 0.0
            && spec.distance.is_none()
            && !spec.reciprocal
            && !spec.same_type
            && !percentage_field;
        Rules {
            match_type: spec.match_type,
            report_coords: spec.report_coords,
            overlap_cutoff: spec.overlap_cutoff,
            distance: spec.distance,
            reciprocal: spec.reciprocal,
            limit: spec.num_records,
            stats_order,
            stats,
            percentage_field,
            span_rule,
        }
    }

    /// Whether the record whose window-clipped bounds are `ref_start..=ref_end`
    /// (1-based) matches the variant on `start..=end` (the variant's own
    /// start and end, an insertion's start one past its end), and the overlap
    /// percentage the match reports (`File.pm:534-575`).
    fn matches(&self, ref_start: u64, ref_end: u64, start: u64, end: u64) -> Option<String> {
        if self.match_type == MatchType::Exact {
            return (ref_start == start && ref_end == end).then(|| "100".to_string());
        }
        let (vs, ve) = (start.min(end), start.max(end));
        match self.match_type {
            MatchType::Within if vs > ref_start || ve < ref_end => return None,
            MatchType::Surrounding if vs < ref_start || ve > ref_end => return None,
            _ => {}
        }
        if let Some(distance) = self.distance {
            if vs.abs_diff(ref_start) as f64 > distance || ve.abs_diff(ref_end) as f64 > distance {
                return None;
            }
        }
        let overlap_start = vs.max(ref_start);
        let overlap_end = ve.min(ref_end);
        let covered = (1 + overlap_end as i64 - overlap_start as i64) as f64;
        let measured_length = if self.reciprocal {
            ve as i64 - vs as i64 + 1
        } else {
            ref_end as i64 - ref_start as i64 + 1
        };
        let percentage = if measured_length != 0 {
            100.0 * covered / measured_length as f64
        } else {
            100.0
        };
        if percentage < self.overlap_cutoff {
            return None;
        }
        (ref_end >= vs && ref_start <= ve).then(|| format!("{percentage:.3}"))
    }
}

/// One `--custom` bigWig source.
pub struct BigWigSource {
    spec: CustomSpec,
    short_name: Arc<str>,
    rules: Rules,
    reader: Mutex<BigWigRead<ReopenableFile>>,
    /// The file's chromosomes by the name a variant reaches them under: the
    /// file's own name with a leading `chr` removed (ensembl-io
    /// `BigFileParser.pm:104-109`).
    chromosomes: HashMap<String, Chromosome>,
    synonyms: Arc<ChromosomeSynonyms>,
}

impl std::fmt::Debug for BigWigSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BigWigSource")
            .field("spec", &self.spec)
            .finish_non_exhaustive()
    }
}

/// The annotations of a batch, by the index of the variant they belong to.
type Prefetched = Vec<(usize, CustomAnnotations)>;

/// A variant's index in its batch with its query window ([`window`]).
type Window = (usize, u32, u32);

impl BigWigSource {
    /// Opens the file of `spec`; a file that is not there or is not a bigWig
    /// is refused with the file named.
    pub fn open(spec: CustomSpec, synonyms: Arc<ChromosomeSynonyms>) -> Result<Self, CustomError> {
        let reader = open_bigwig(&spec.file).map_err(|e| CustomError::Open {
            file: spec.file.clone(),
            reason: e.to_string(),
        })?;
        let chromosomes = reader
            .chroms()
            .iter()
            .map(|chromosome| {
                let reached_as = chromosome
                    .name
                    .strip_prefix("chr")
                    .unwrap_or(&chromosome.name);
                (
                    reached_as.to_string(),
                    Chromosome {
                        name: chromosome.name.clone(),
                        length: chromosome.length,
                    },
                )
            })
            .collect();
        Ok(BigWigSource {
            short_name: Arc::from(spec.short_name.as_str()),
            rules: Rules::from_spec(&spec),
            reader: Mutex::new(reader),
            chromosomes,
            synonyms,
            spec,
        })
    }

    /// The file's chromosome an input chromosome name reaches, through the
    /// run's synonyms and the `chr` rule (`BaseVEP.pm:689-738`); `None` when
    /// the file has no such chromosome.
    fn chromosome(&self, chr: &str) -> Option<&Chromosome> {
        let name = self
            .synonyms
            .source_chr_name(chr, |name| self.chromosomes.contains_key(name));
        self.chromosomes.get(&name)
    }

    fn read_error(&self, error: BBIReadError) -> CustomError {
        CustomError::Read {
            file: self.spec.file.clone(),
            reason: error.to_string(),
        }
    }

    /// The intervals of `chromosome` the file holds over the 0-based half-open
    /// `start..end`, clipped to it, in position order.
    fn intervals(
        &self,
        reader: &mut BigWigRead<ReopenableFile>,
        chromosome: &Chromosome,
        start: u32,
        end: u32,
    ) -> Result<Vec<Value>, CustomError> {
        if end <= start {
            return Ok(Vec::new());
        }
        reader
            .get_interval(&chromosome.name, start, end)
            .map_err(|e| self.read_error(e))?
            .map(|value| value.map_err(|e| self.read_error(e)))
            .collect()
    }

    /// The name a record prints: its value as Ensembl VEP stringifies the
    /// number, or under `coords` the file's chromosome name with the record's
    /// 1-based start and end (`BigWig.pm:391-393`, on both paths).
    fn record_name(&self, chromosome: &Chromosome, interval: Value) -> String {
        if self.rules.report_coords {
            format!(
                "{}:{}-{}",
                chromosome.name,
                interval.start + 1,
                interval.end
            )
        } else {
            perl_number_string(f64::from(interval.value))
        }
    }

    /// The annotation of one variant whose span the per-record walk covers
    /// (`File.pm:316-352`), from the intervals of its window
    /// ([`window`]) in position order, each clipped to it: `None` when no
    /// record matches and no statistics were asked for, since `File.pm:344-352`
    /// closes the statistics of every variant the walk visits.
    fn annotate_window(
        &self,
        chromosome: &Chromosome,
        variant: &InputVariant,
        intervals: impl Iterator<Item = Value>,
    ) -> Option<CustomAnnotations> {
        let rules = &self.rules;
        let mut records = Vec::new();
        let mut accumulator = Accumulator::default();
        for interval in intervals {
            let ref_start = u64::from(interval.start) + 1;
            let ref_end = u64::from(interval.end);
            let Some(percentage) = rules.matches(ref_start, ref_end, variant.start, variant.end)
            else {
                continue;
            };
            let value = f64::from(interval.value);
            if rules.stats_order.is_some() {
                accumulator.add(value, rules.stats);
            }
            let listed = accumulator.count;
            accumulator.count += 1;
            match rules.limit {
                RecordLimit::Count(limit) if listed == limit => {
                    records.push(ellipsis());
                    if rules.stats_order.is_none() {
                        break;
                    }
                }
                RecordLimit::Count(limit) if listed > limit => {}
                _ => {
                    let mut fields = IndexMap::new();
                    if rules.percentage_field {
                        fields.insert("PC".to_string(), CustomFieldValue::Text(percentage));
                    }
                    records.push(CustomRecord {
                        name: self.record_name(chromosome, interval),
                        allele: None,
                        fields,
                        score: rules
                            .stats_order
                            .is_some()
                            .then_some(CustomValue::Number(value)),
                    });
                }
            }
        }
        let stats = rules
            .stats_order
            .as_deref()
            .map(|order| accumulator.close(order));
        if records.is_empty() && stats.is_none() {
            return None;
        }
        Some(self.annotations(records, stats))
    }

    /// The annotations of the variants of one chromosome that the per-record
    /// walk covers, read in one query per cluster of windows within
    /// [`CLUSTER_GAP`] of each other, each variant's records being the cluster's
    /// intervals that reach into its window, clipped to it: the intervals a
    /// query over the window alone returns.
    fn annotate_clusters(
        &self,
        reader: &mut BigWigRead<ReopenableFile>,
        chromosome: &Chromosome,
        windows: &mut [Window],
        variants: &[InputVariant],
        annotated: &mut Prefetched,
    ) -> Result<(), CustomError> {
        windows.sort_unstable_by_key(|&(_, start, _)| start);
        let mut first = 0;
        while first < windows.len() {
            let (_, cluster_start, mut cluster_end) = windows[first];
            let mut past = first + 1;
            while past < windows.len() && windows[past].1 <= cluster_end.saturating_add(CLUSTER_GAP)
            {
                cluster_end = cluster_end.max(windows[past].2);
                past += 1;
            }
            let intervals = self.intervals(reader, chromosome, cluster_start, cluster_end)?;
            for &(index, start, end) in &windows[first..past] {
                let from = intervals.partition_point(|interval| interval.end <= start);
                let in_window = intervals[from..]
                    .iter()
                    .take_while(|interval| interval.start < end)
                    .map(|interval| Value {
                        start: interval.start.max(start),
                        end: interval.end.min(end),
                        value: interval.value,
                    });
                if let Some(annotation) =
                    self.annotate_window(chromosome, &variants[index], in_window)
                {
                    annotated.push((index, annotation));
                }
            }
            first = past;
        }
        Ok(())
    }

    /// The annotation of one variant under the span rule
    /// (`BigWig.pm:264-296`): the extrema over the whole span, the records of
    /// its first window and `...` where the span continues; `None` when the
    /// span holds no value.
    fn annotate_long_span(
        &self,
        reader: &mut BigWigRead<ReopenableFile>,
        chromosome: &Chromosome,
        variant: &InputVariant,
    ) -> Result<Option<CustomAnnotations>, CustomError> {
        let (vs, ve) = (
            variant.start.min(variant.end),
            variant.start.max(variant.end),
        );
        let start = u32::try_from(vs.saturating_sub(1)).unwrap_or(u32::MAX);
        let end = u32::try_from(ve).unwrap_or(u32::MAX).min(chromosome.length);
        if end <= start {
            return Ok(None);
        }
        let mut accumulator = Accumulator::default();
        for interval in self.intervals(reader, chromosome, start, end)? {
            accumulator.add(f64::from(interval.value), self.rules.stats);
        }
        if accumulator.min.is_none() && accumulator.max.is_none() {
            return Ok(None);
        }
        let order = self.rules.stats_order.as_deref().unwrap_or_default();
        let stats = Some(accumulator.close(order));
        let mut records = Vec::new();
        let window_end = start.saturating_add(SPAN_RECORD_WINDOW).min(end);
        match self.rules.limit {
            RecordLimit::Count(0) => records.push(ellipsis()),
            limit => {
                for interval in self.intervals(reader, chromosome, start, window_end)? {
                    if matches!(limit, RecordLimit::Count(limit) if records.len() as u64 >= limit) {
                        records.push(ellipsis());
                        return Ok(Some(self.annotations(records, stats)));
                    }
                    records.push(CustomRecord {
                        name: self.record_name(chromosome, interval),
                        allele: None,
                        fields: IndexMap::new(),
                        score: None,
                    });
                }
                if end > window_end {
                    records.push(ellipsis());
                }
            }
        }
        Ok(Some(self.annotations(records, stats)))
    }

    fn annotations(
        &self,
        records: Vec<CustomRecord>,
        stats: Option<CustomStats>,
    ) -> CustomAnnotations {
        CustomAnnotations {
            short_name: Arc::clone(&self.short_name),
            matched: !records.is_empty(),
            records,
            stats,
        }
    }

    /// Whether `variant` takes the span rule: the source allows it and the
    /// variant's span exceeds [`SPAN_RULE_MIN_SPAN`] (`BigWig.pm:186`).
    fn takes_span_rule(&self, variant: &InputVariant) -> bool {
        self.rules.span_rule && variant.start.abs_diff(variant.end) + 1 > SPAN_RULE_MIN_SPAN
    }
}

/// The query window of a variant, 0-based half-open: from one base before its
/// start to one past its end (`File.pm:321`); an insertion's start is one past
/// its end, so its window is the one base after the insertion point.
fn window(variant: &InputVariant) -> (u32, u32) {
    (
        u32::try_from(variant.start.saturating_sub(1)).unwrap_or(u32::MAX),
        u32::try_from(variant.end.saturating_add(1)).unwrap_or(u32::MAX),
    )
}

fn ellipsis() -> CustomRecord {
    CustomRecord {
        name: ELLIPSIS.to_string(),
        allele: None,
        fields: IndexMap::new(),
        score: None,
    }
}

impl CustomSource for BigWigSource {
    fn spec(&self) -> &CustomSpec {
        &self.spec
    }

    /// The file's chromosomes less a leading `chr`, the names the parser's
    /// cache keys them by (`BigWig.pm:147-149`, `BigFileParser.pm:104-109`).
    fn chromosomes(&self) -> Vec<String> {
        self.chromosomes.keys().cloned().collect()
    }

    /// Reads and matches the batch's records: one query per cluster of nearby
    /// windows on a chromosome, and one per variant under the span rule. A
    /// variant Ensembl VEP skips (`vep_skip`) is left alone (`File.pm:318`);
    /// one on a chromosome the file lacks keeps the statistics entry the
    /// closing writes for every other variant (`:344-352`), unless it takes the
    /// span rule, whose path returns before any entry (`BigWig.pm:270`).
    fn prefetch(&self, variants: &[InputVariant]) -> Result<PrefetchData, CustomError> {
        let mut reader = self
            .reader
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let mut annotated: Prefetched = Vec::new();
        let mut windows_by_chromosome: HashMap<&str, (&Chromosome, Vec<Window>)> = HashMap::new();
        for (index, variant) in variants.iter().enumerate() {
            if variant.vep_skip {
                continue;
            }
            let chr = if variant.original_chr.is_empty() {
                variant.chr.as_str()
            } else {
                variant.original_chr.as_str()
            };
            let Some(chromosome) = self.chromosome(chr) else {
                if let (false, Some(order)) = (
                    self.takes_span_rule(variant),
                    self.rules.stats_order.as_deref(),
                ) {
                    let stats = Accumulator::default().close(order);
                    annotated.push((index, self.annotations(Vec::new(), Some(stats))));
                }
                continue;
            };
            if self.takes_span_rule(variant) {
                if let Some(annotation) =
                    self.annotate_long_span(&mut reader, chromosome, variant)?
                {
                    annotated.push((index, annotation));
                }
                continue;
            }
            let (start, end) = window(variant);
            windows_by_chromosome
                .entry(chromosome.name.as_str())
                .or_insert_with(|| (chromosome, Vec::new()))
                .1
                .push((index, start, end));
        }
        for (chromosome, mut windows) in windows_by_chromosome.into_values() {
            self.annotate_clusters(
                &mut reader,
                chromosome,
                &mut windows,
                variants,
                &mut annotated,
            )?;
        }
        annotated.sort_unstable_by_key(|(index, _)| *index);
        Ok(PrefetchData::new(annotated))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), CustomError> {
        let annotated: Prefetched = data.downcast().map_err(|e| CustomError::Read {
            file: self.spec.file.clone(),
            reason: e.to_string(),
        })?;
        for (index, annotation) in annotated {
            variants[index].custom_annotations.push(annotation);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::custom::grammar::parse_custom;
    use std::path::Path;

    /// A committed fixture under `tests/data`, each written by the `bigtools`
    /// command line's `bedgraphtobigwig` over a bedGraph and a `chrom.sizes` of
    /// `chr21` (46,709,983 bases) and `CM000684.2` (50,818,468 bases).
    /// `cons.bw` holds, on `chr21`, the one-base intervals at 0-based 99, 100,
    /// 101, 102 and 105 with the values 10, 20, 30, 0 and -1.5, the interval
    /// 199-210 with the value 7.25 and the interval 300-302 with the value 0.1,
    /// and on `CM000684.2` the one-base interval at 49 with the value 5.
    /// `sparse.bw` holds, on `chr21`, one-base intervals at 1,000,010,
    /// 1,000,500, 1,050,000, 1,060,000 and 1,090,000 with the values 3, 1,
    /// -4, 2 and 9.
    fn fixture(name: &str) -> String {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
            .display()
            .to_string()
    }

    fn opened(argument: &str, synonyms: ChromosomeSynonyms) -> BigWigSource {
        BigWigSource::open(parse_custom(argument).unwrap(), Arc::new(synonyms)).unwrap()
    }

    fn source_with_stats(file: &str, match_type: &str, stats: &str) -> BigWigSource {
        opened(
            &format!(
                "file={file},short_name=cons,format=bigwig,type={match_type},summary_stats={stats}"
            ),
            ChromosomeSynonyms::new(),
        )
    }

    fn variant(chr: &str, start: u64, end: u64) -> InputVariant {
        let mut variant = InputVariant::new(chr.into(), start, end, b"A".to_vec(), b"G".to_vec());
        variant.original_chr = chr.to_string();
        variant
    }

    /// The annotations the source leaves on each variant of `batch`, by index.
    fn annotate(source: &BigWigSource, batch: &mut [InputVariant]) -> Vec<Vec<CustomAnnotations>> {
        let data = source.prefetch(batch).unwrap();
        source.annotate(batch, data).unwrap();
        batch
            .iter()
            .map(|variant| variant.custom_annotations.clone())
            .collect()
    }

    fn names(annotation: &CustomAnnotations) -> Vec<&str> {
        annotation
            .records
            .iter()
            .map(|record| record.name.as_str())
            .collect()
    }

    /// The statistics by name, each a number.
    fn stats(annotation: &CustomAnnotations) -> IndexMap<String, f64> {
        let Some(CustomStats::Values(values)) = &annotation.stats else {
            panic!("statistics were asked for");
        };
        values
            .iter()
            .map(|(name, value)| match value {
                CustomValue::Number(number) => (name.clone(), *number),
                other => panic!("a bigWig statistic is a number, not {other:?}"),
            })
            .collect()
    }

    fn number(value: f64) -> Option<CustomValue> {
        Some(CustomValue::Number(value))
    }

    fn text(value: &str) -> CustomFieldValue {
        CustomFieldValue::Text(value.to_string())
    }

    /// Under `exact` a record matches when its window-clipped bounds, the
    /// 0-based start read as 1-based, equal the variant's (`File.pm:529-532,
    /// 571-575`): the one-base interval at the variant lists its value as the
    /// record's name and score (`t/AnnotationSource_File_BigWig.t:104-106`,
    /// `{ name => 10, score => 10 }` with statistics on), and an interval
    /// spanning the variant lists nothing, since the query window reaches one
    /// base past the variant's end and the clipped record ends there.
    #[test]
    fn exact_lists_the_interval_value_at_the_variant() {
        let cons = fixture("cons.bw");
        let source = opened(
            &format!("{cons},cons,bigwig,exact"),
            ChromosomeSynonyms::new(),
        );
        let mut batch = vec![variant("21", 100, 100), variant("21", 205, 205)];
        let annotated = annotate(&source, &mut batch);
        assert_eq!(names(&annotated[0][0]), ["10"]);
        assert_eq!(annotated[0][0].records[0].score, None);
        assert_eq!(annotated[0][0].stats, None);
        assert!(annotated[1].is_empty());

        let with_stats = source_with_stats(&cons, "exact", "min%max%mean%sum%count");
        let mut batch = vec![variant("21", 100, 100)];
        let annotated = annotate(&with_stats, &mut batch);
        let annotation = &annotated[0][0];
        assert_eq!(&*annotation.short_name, "cons");
        assert_eq!(names(annotation), ["10"]);
        assert_eq!(annotation.records[0].score, number(10.0));
        assert_eq!(
            stats(annotation),
            IndexMap::from([
                ("min".to_string(), 10.0),
                ("max".to_string(), 10.0),
                ("mean".to_string(), 10.0),
                ("sum".to_string(), 10.0),
                ("count".to_string(), 1.0),
            ])
        );
    }

    /// `File.pm:344-352` closes the statistics for every variant the walk
    /// visits once `summary_stats` is given, a variant no interval matches
    /// included, one on a chromosome the file lacks among them, so such a
    /// variant keeps an entry with no record: `count` 0 when asked for; nothing
    /// when `min`, `max` or `mean` alone were asked for, whose hash the closing
    /// creates empty; undefined when `sum` was asked for without `count` or
    /// `mean`, since no step creates the hash. Without statistics the variant
    /// keeps no entry, nor does a variant that takes the span rule on a
    /// chromosome the file lacks, since that path returns before any entry
    /// (`BigWig.pm:270`).
    #[test]
    fn a_variant_no_interval_matches_keeps_its_statistics_entry() {
        let cons = fixture("cons.bw");
        let unmatched = |stats: &str| {
            let source = source_with_stats(&cons, "overlap", stats);
            let mut batch = vec![variant("21", 150, 150), variant("9", 150, 150)];
            let mut annotated = annotate(&source, &mut batch);
            assert_eq!(annotated[0], annotated[1], "{stats}");
            annotated.remove(0)
        };
        let entry = unmatched("count").remove(0);
        assert!(entry.records.is_empty() && !entry.matched);
        assert_eq!(stats(&entry), IndexMap::from([("count".to_string(), 0.0)]));
        for empty in ["min", "max%min", "mean", "sum%mean"] {
            let entry = unmatched(empty).remove(0);
            assert!(entry.records.is_empty() && !entry.matched);
            assert!(stats(&entry).is_empty(), "{empty}");
        }
        for undefined in ["sum", "sum%min%max"] {
            let entry = unmatched(undefined).remove(0);
            assert_eq!(entry.stats, Some(CustomStats::Undefined), "{undefined}");
        }
        let without_stats = opened(
            &format!("{cons},cons,bigwig,overlap"),
            ChromosomeSynonyms::new(),
        );
        let mut batch = vec![variant("21", 150, 150), variant("9", 150, 150)];
        assert!(annotate(&without_stats, &mut batch)
            .iter()
            .all(Vec::is_empty));
        let span_rule = source_with_stats(&cons, "overlap", "min%max");
        let mut batch = vec![variant("9", 1, 100_000)];
        assert!(annotate(&span_rule, &mut batch)[0].is_empty());
    }

    /// Under `overlap` every interval sharing a base with the variant is one
    /// record, named by its value (`BigWig.pm:427`, `t/AnnotationSource_File_BigWig.t:79-81`),
    /// and the statistics run over the matched values (`File.pm:422-434,
    /// 344-352`): the intervals 10, 20 and 30 over a three-base variant give
    /// min 10, max 30, mean 20, sum 60 and count 3; a value of 0 is a minimum
    /// (`File.pm:427` tests definedness, not truth). A value the file stores
    /// as the single-precision 0.1 prints as Ensembl VEP stringifies the
    /// widened number, `0.100000001490116`.
    #[test]
    fn overlap_lists_one_record_per_interval_with_stats() {
        let cons = fixture("cons.bw");
        let source = source_with_stats(&cons, "overlap", "min%max%mean%sum%count");
        let mut batch = vec![
            variant("21", 100, 102),
            variant("21", 100, 103),
            variant("21", 301, 301),
        ];
        let annotated = annotate(&source, &mut batch);
        let three = &annotated[0][0];
        assert_eq!(names(three), ["10", "20", "30"]);
        assert_eq!(
            stats(three),
            IndexMap::from([
                ("min".to_string(), 10.0),
                ("max".to_string(), 30.0),
                ("mean".to_string(), 20.0),
                ("sum".to_string(), 60.0),
                ("count".to_string(), 3.0),
            ])
        );
        let with_zero = &annotated[1][0];
        assert_eq!(names(with_zero), ["10", "20", "30", "0"]);
        assert_eq!(stats(with_zero)["min"], 0.0);
        assert_eq!(stats(with_zero)["count"], 4.0);
        assert_eq!(names(&annotated[2][0]), ["0.100000001490116"]);
        assert_eq!(annotated[2][0].records[0].score, number(f64::from(0.1f32)));
    }

    /// `coords` names a record by the file's chromosome name and the record's
    /// 1-based start and end (`BigWig.pm:391-392`; Ensembl VEP's per-record
    /// path prints `21:99-100` here, the `chr`-stripped name with the 0-based
    /// start, `BigWig.pm:420-427`), and `fields=PC` adds the share of the
    /// record the variant covers to three decimals, 100 for an `exact` match
    /// (`File.pm:418-420, 554-574`); `num_records` caps the list with `...`
    /// (`File.pm:436-440`) while the statistics run over every match.
    #[test]
    fn coords_percentage_and_record_cap_follow_file_pm() {
        let cons = fixture("cons.bw");
        let source = opened(
            &format!(
                "file={cons},short_name=cons,format=bigwig,type=overlap,coords=1,fields=PC,num_records=2,summary_stats=count"
            ),
            ChromosomeSynonyms::new(),
        );
        let mut batch = vec![variant("21", 100, 103), variant("21", 205, 205)];
        let annotated = annotate(&source, &mut batch);
        let capped = &annotated[0][0];
        assert_eq!(names(capped), ["chr21:100-100", "chr21:101-101", "..."]);
        assert_eq!(capped.records[0].fields["PC"], text("100.000"));
        assert_eq!(stats(capped)["count"], 4.0);
        let partial = &annotated[1][0];
        assert_eq!(names(partial), ["chr21:205-206"]);
        assert_eq!(partial.records[0].fields["PC"], text("50.000"));

        let exact = opened(
            &format!("file={cons},short_name=cons,format=bigwig,type=exact,fields=PC"),
            ChromosomeSynonyms::new(),
        );
        let mut batch = vec![variant("21", 100, 100)];
        let annotated = annotate(&exact, &mut batch);
        assert_eq!(annotated[0][0].records[0].fields["PC"], text("100"));
    }

    /// A span over 65,536 bases under `overlap` with `min` and `max` alone
    /// takes the span rule (`BigWig.pm:186, 239-249, 264-296`): the extrema
    /// are exact over the whole span, and the record list holds the intervals
    /// of the span's first 32,768 bases then `...` (`BigWig.pm:369-401`). The
    /// values -4 and 9 lie past that window and reach the statistics only.
    #[test]
    fn large_span_min_max_lists_records_from_the_first_window_then_ellipsis() {
        let sparse = fixture("sparse.bw");
        let source = source_with_stats(&sparse, "overlap", "min%max");
        let mut batch = vec![variant("21", 1_000_001, 1_100_000)];
        let annotated = annotate(&source, &mut batch);
        let annotation = &annotated[0][0];
        assert_eq!(names(annotation), ["3", "1", "..."]);
        assert_eq!(
            stats(annotation),
            IndexMap::from([("min".to_string(), -4.0), ("max".to_string(), 9.0)])
        );
        for record in &annotation.records {
            assert_eq!(record.score, None);
        }

        let only_ellipsis = opened(
            &format!(
                "file={sparse},short_name=cons,format=bigwig,type=overlap,summary_stats=max,num_records=0"
            ),
            ChromosomeSynonyms::new(),
        );
        let mut batch = vec![
            variant("21", 1_000_001, 1_100_000),
            variant("21", 2_000_001, 2_100_000),
        ];
        let annotated = annotate(&only_ellipsis, &mut batch);
        assert_eq!(names(&annotated[0][0]), ["..."]);
        assert_eq!(
            stats(&annotated[0][0]),
            IndexMap::from([("max".to_string(), 9.0)])
        );
        assert!(annotated[1].is_empty());
    }

    /// On the `min` and `max` path a span of at most 65,536 bases is walked
    /// record by record once: one record per interval. Ensembl VEP runs the
    /// walk twice there (`BigWig.pm:195-204`) and lists each record twice.
    #[test]
    fn small_spans_are_annotated_once_on_the_min_max_path() {
        let cons = fixture("cons.bw");
        let source = source_with_stats(&cons, "overlap", "min%max");
        let mut batch = vec![variant("21", 100, 102)];
        let annotated = annotate(&source, &mut batch);
        let annotation = &annotated[0][0];
        assert_eq!(names(annotation), ["10", "20", "30"]);
        assert_eq!(
            stats(annotation),
            IndexMap::from([("min".to_string(), 10.0), ("max".to_string(), 30.0)])
        );
    }

    /// A batch is read in one query per cluster of windows, and every variant
    /// of it is annotated as a batch of itself alone annotates it: the same
    /// records, in the same order, whatever the batch's order and however the
    /// windows share, overlap or straddle a cluster's intervals. Spans under
    /// the span rule read on their own in the same batch.
    #[test]
    fn clustered_reads_equal_one_query_per_variant() {
        let cons = fixture("cons.bw");
        for argument in [
            format!("{cons},cons,bigwig,exact"),
            format!("file={cons},short_name=cons,format=bigwig,type=overlap,coords=1,fields=PC,summary_stats=min%max%mean%sum%count"),
            format!("file={cons},short_name=cons,format=bigwig,type=within,num_records=1"),
            format!("file={cons},short_name=cons,format=bigwig,type=overlap,summary_stats=min%max"),
        ] {
            let mut synonyms = ChromosomeSynonyms::new();
            synonyms.add_synonym("22", "CM000684.2");
            let source = opened(&argument, synonyms);
            let mut batch: Vec<InputVariant> = (95..=320)
                .rev()
                .map(|position| variant("21", position, position))
                .chain((98..=104).map(|start| variant("21", start, start + 2)))
                .chain([
                    variant("21", 101, 100),
                    variant("21", 1, 100_000),
                    variant("22", 50, 50),
                    variant("22", 49, 51),
                    variant("X", 100, 100),
                ])
                .collect();
            let mut skipped = variant("21", 100, 100);
            skipped.vep_skip = true;
            batch.push(skipped);
            let alone: Vec<Vec<CustomAnnotations>> = batch
                .iter()
                .map(|variant| annotate(&source, &mut [variant.clone()]).remove(0))
                .collect();
            let together = annotate(&source, &mut batch);
            assert_eq!(together, alone, "{argument}");
            assert!(together.iter().any(|annotations| !annotations.is_empty()));
        }
    }

    /// The file's chromosomes are reached by their names with a leading `chr`
    /// removed (ensembl-io `BigFileParser.pm:104-109`), through the run's
    /// synonyms and the `chr` rule (`BaseVEP.pm:689-738`): `21` and `chr21`
    /// reach the file's `chr21`; `22` reaches `CM000684.2` only through a
    /// synonym; a chromosome the file lacks annotates nothing.
    #[test]
    fn bigwig_chromosomes_map_through_synonyms_and_chr_prefix() {
        let cons = fixture("cons.bw");
        let mut synonyms = ChromosomeSynonyms::new();
        synonyms.add_synonym("22", "CM000684.2");
        let source = opened(&format!("{cons},cons,bigwig,exact"), synonyms);
        let mut batch = vec![
            variant("21", 100, 100),
            variant("chr21", 100, 100),
            variant("22", 50, 50),
            variant("CM000684.2", 50, 50),
            variant("X", 100, 100),
        ];
        let annotated = annotate(&source, &mut batch);
        assert_eq!(names(&annotated[0][0]), ["10"]);
        assert_eq!(names(&annotated[1][0]), ["10"]);
        assert_eq!(names(&annotated[2][0]), ["5"]);
        assert_eq!(names(&annotated[3][0]), ["5"]);
        assert!(annotated[4].is_empty());

        let unmapped = opened(
            &format!("{cons},cons,bigwig,exact"),
            ChromosomeSynonyms::new(),
        );
        let mut batch = vec![variant("22", 50, 50)];
        assert!(annotate(&unmapped, &mut batch)[0].is_empty());
    }
}
