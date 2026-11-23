// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Conservation plugin: the GERP conservation score over a variant's span, read
//! from a bigWig of Ensembl Compara's per-base scores (`Conservation.pm`,
//! VEP_plugins release/116, in its bigWig mode).
//!
//! The plugin takes the bigWig as its first parameter (or `file=`) and the word
//! `MAX` as an optional last one of the positional list (`:82-89`: the method is
//! read from the last parameter alone, `MAX` for the maximum, anything else the
//! average; the keyed reading of [`PluginParams`] has no flag token, so `file=`
//! goes without `MAX`). Another positional parameter after the bigWig is not
//! read, as `:91` and `:192` read the first alone; it is logged. The plugin
//! writes one field, `Conservation`, at the record level, so every row of a
//! record carries it, intergenic rows included (`feature_types
//! ['Feature','Intergenic']`, `:133-135`). The header names the method
//! (`:137-144`).
//!
//! The span read is the record's own bases, or for an insertion the two bases
//! that flank it (`:205-211`). Every base the file covers within the span
//! contributes the value of its interval once (`:216-236`: an interval of `n`
//! bases contributes its value `n` times; a base without data contributes
//! nothing); a span the file covers nowhere gets no field (`:213`). One covered
//! base prints its value; more print their average or, under `MAX`, their
//! maximum, each as `%.3f` (`:239-254`). Chromosomes are matched as ensembl-io's
//! bigWig parser matches them: the record's name less a leading `chr` in any
//! case (`:202-203`) against the file's names less a leading `chr`
//! (`BigFileParser.pm:104-109`). Structural variants get no value: the plugin
//! keeps the default `variant_feature_types` of `['VariationFeature']`
//! (BaseVepPlugin.pm:100), so Ensembl VEP never hands it one.
//!
//! The other two modes of the Ensembl VEP plugin are refused by name: a first
//! parameter `database` (`:91-93`, the Compara database) and a first parameter
//! that is not a local file (`:155-190`, a species set or FTP address looked up
//! on Ensembl's FTP site at run time).

use std::collections::HashMap;
use std::path::Path;
use std::sync::Mutex;

use bigtools::utils::file::reopen::ReopenableFile;
use bigtools::{BigWigRead, Value};
use indexmap::IndexMap;
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::variant::InputVariant;

use crate::bigwig::open_bigwig;
use crate::params::{warn_unread_parameters, Grammar, PluginParams};
use crate::perl_format::fixed;
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};

const NAME: &str = "Conservation";

/// The output field (`:142`).
const FIELD: &str = "Conservation";

/// The parameter that selects the maximum over the average (`:83`).
const MAX: &str = "MAX";

/// The first parameter of the refused database mode (`:91`).
const DATABASE: &str = "database";

/// How several covered bases are reduced to one score (`:81-89`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum Method {
    #[default]
    Average,
    Max,
}

impl Method {
    /// The word of the header text (`:140`).
    fn header_word(self) -> &'static str {
        match self {
            Method::Average => "average",
            Method::Max => "maximum",
        }
    }
}

/// Conservation plugin.
pub struct ConservationPlugin {
    method: Method,
    /// The opened bigWig, after `init`.
    source: Option<Source>,
}

/// An opened bigWig with its chromosome names.
struct Source {
    file: String,
    /// `BigWigRead::get_interval` takes `&mut self`; the plugin runs its batches
    /// sequentially, so the lock is never contended.
    reader: Mutex<BigWigRead<ReopenableFile>>,
    /// The file's chromosome names by the name a record reaches them under:
    /// the file's own name with a leading `chr` removed, a later chromosome
    /// of the file winning a shared key, as the hash of
    /// `BigFileParser.pm:104-109` has it.
    chromosomes: HashMap<String, String>,
}

impl Source {
    fn open(file: &str) -> Result<Self, PluginError> {
        let reader = open_bigwig(file).map_err(|error| {
            PluginError::Init(format!("{NAME}: {file} is not a readable bigWig: {error}"))
        })?;
        let chromosomes = reader
            .chroms()
            .iter()
            .map(|chromosome| {
                let name = chromosome.name.as_str();
                (
                    name.strip_prefix("chr").unwrap_or(name).to_string(),
                    name.to_string(),
                )
            })
            .collect();
        Ok(Self {
            file: file.to_string(),
            reader: Mutex::new(reader),
            chromosomes,
        })
    }

    /// The file's chromosome a record's name reaches: the name less a leading
    /// `chr` in any case (`:202-203`), looked up among the file's `chr`-less
    /// names.
    fn chromosome(&self, chr: &str) -> Option<&str> {
        let key = match chr.get(..3) {
            Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &chr[3..],
            _ => chr,
        };
        self.chromosomes.get(key).map(String::as_str)
    }

    /// The scores of a batch, by the index of the variant each belongs to: one
    /// query per cluster of nearby windows on a chromosome, each window then
    /// reading its own bases from the cluster's intervals, so the data block
    /// a run of close records lies in is inflated once, not once per record.
    fn scores(&self, variants: &[InputVariant], method: Method) -> Result<Scores, PluginError> {
        let mut windows_by_chromosome: HashMap<&str, Vec<Window>> = HashMap::new();
        for (index, variant) in variants.iter().enumerate() {
            if variant.variant_class.is_structural() {
                continue;
            }
            let Some(chromosome) = self.chromosome(&variant.original_chr) else {
                continue;
            };
            if let Some((start, end)) = window(variant) {
                windows_by_chromosome
                    .entry(chromosome)
                    .or_default()
                    .push((index, start, end));
            }
        }
        let mut reader = self
            .reader
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let mut scores = Vec::new();
        for (chromosome, mut windows) in windows_by_chromosome {
            windows.sort_unstable_by_key(|&(_, start, _)| start);
            let mut first = 0;
            while first < windows.len() {
                let (_, cluster_start, mut cluster_end) = windows[first];
                let mut past = first + 1;
                while past < windows.len()
                    && windows[past].1 <= cluster_end.saturating_add(CLUSTER_GAP)
                {
                    cluster_end = cluster_end.max(windows[past].2);
                    past += 1;
                }
                let intervals: Vec<Value> = reader
                    .get_interval(chromosome, cluster_start, cluster_end)
                    .map_err(|error| self.read_error(error))?
                    .collect::<Result<_, _>>()
                    .map_err(|error| self.read_error(error))?;
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
                    if let Some(score) = score(in_window, method) {
                        scores.push((index, score));
                    }
                }
                first = past;
            }
        }
        scores.sort_unstable_by_key(|(index, _)| *index);
        Ok(scores)
    }

    fn read_error(&self, error: bigtools::BBIReadError) -> PluginError {
        PluginError::Run(format!("{NAME}: reading {}: {error}", self.file))
    }
}

/// Two variants' windows this many bases apart or closer are read in one
/// query: a bigWig data block holds about a thousand items and is inflated
/// once per query, so the bases between two windows within one block cost
/// less than a second walk of the file's index and a second inflation.
const CLUSTER_GAP: u32 = 1_000;

/// A variant's index in its batch with its 0-based half-open query window.
type Window = (usize, u32, u32);

/// A batch's scores, each with the index of the variant it belongs to.
type Scores = Vec<(usize, String)>;

/// The 0-based half-open window a variant reads (`:205-211`): its own bases,
/// or for an insertion (start one past end) the two bases it falls between;
/// `None` for the empty window of a record whose start and end are both 0. An
/// allele of a `--minimal` split record reads its own trimmed span (the
/// documented class `minimal_plugin_sees_the_rejoined_record`).
fn window(variant: &InputVariant) -> Option<(u32, u32)> {
    let (start, end) = (variant.start, variant.end);
    let (first, last) = if start > end {
        (end, start)
    } else {
        (start, end)
    };
    let start = first.saturating_sub(1) as u32;
    let end = last as u32;
    (end > start).then_some((start, end))
}

/// The score over the covered bases of one window (`:216-254`), from the
/// intervals clipped to it in position order: `None` where the file covers no
/// base of it. One value per covered base, summed in position order and
/// compared one by one, so the double arithmetic is the plugin's own. One
/// covered base is not summed: its stored value is the score, so a stored zero
/// keeps its sign.
fn score(intervals: impl Iterator<Item = Value>, method: Method) -> Option<String> {
    let mut covered = 0u64;
    let mut total = 0.0f64;
    let mut maximum = f64::NEG_INFINITY;
    let mut latest = 0.0f64;
    for interval in intervals {
        let value = f64::from(interval.value);
        for _ in interval.start..interval.end {
            total += value;
        }
        if value > maximum {
            maximum = value;
        }
        if interval.end > interval.start {
            latest = value;
        }
        covered += u64::from(interval.end - interval.start);
    }
    if covered == 0 {
        return None;
    }
    let score = match method {
        Method::Max => maximum,
        // One covered base prints its value as stored (`:252-254`); the sum
        // from zero of `:246-248` would drop the sign of a stored zero.
        Method::Average if covered == 1 => latest,
        Method::Average => total / covered as f64,
    };
    Some(fixed(score, 3))
}

impl Default for ConservationPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl ConservationPlugin {
    /// Create an uninitialized plugin.
    pub fn new() -> Self {
        Self {
            method: Method::default(),
            source: None,
        }
    }
}

impl crate::traits::sealed::Sealed for ConservationPlugin {}

impl BuiltinPlugin for ConservationPlugin {
    fn name(&self) -> &str {
        NAME
    }

    /// The Ensembl VEP module declares `['Feature', 'Intergenic']`
    /// (`Conservation.pm:133-134`): every row, the intergenic row included.
    fn feature_types(&self) -> &[RowKind] {
        RowKind::ALL
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            FIELD.to_string(),
            format!(
                "The {} conservation score for this site",
                self.method.header_word()
            ),
        )]
    }

    /// The bigWig is the first parameter (`:192`) or `file=`, and `MAX` closes
    /// the list (`:82-89`); another positional parameter after the bigWig is
    /// not read, as `:91` and `:192` read the first alone. The database mode
    /// and a first parameter that is not a local file are refused.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        if params.first().is_some_and(|first| first == DATABASE) {
            return Err(PluginError::Init(format!(
                "{NAME}: the database mode (the Ensembl Compara database) is not supported; \
                 give the path of a GERP bigWig"
            )));
        }
        let max_before_the_last_token = params
            .split_last()
            .is_some_and(|(_, before_last)| before_last.iter().any(|param| param == MAX));
        let params = PluginParams::parse_with(
            NAME,
            Grammar::PositionalOrKeyValue { flags: &[MAX] },
            &["file"],
            params,
        )?;
        // Every `MAX` but a closing one sits where `:91` and `:192` read the
        // file or where `:82-89` forget it.
        if max_before_the_last_token {
            return Err(PluginError::Init(format!(
                "{NAME}: {MAX} must be the last parameter"
            )));
        }
        let max_closes_the_list = params.trailing_flag() == Some(MAX);
        let (file, unread) = match (params.get("file"), params.positional()) {
            (Some(file), _) => (file, &[][..]),
            (None, []) => {
                return Err(PluginError::Init(format!(
                    "{NAME}: no bigWig file given; use {NAME},<file.bw> or {NAME},<file.bw>,{MAX}"
                )));
            }
            (None, [file, unread @ ..]) => (file.as_str(), unread),
        };
        if !Path::new(file).is_file() {
            return Err(PluginError::Init(format!(
                "{NAME}: bigWig file {file} not found; give the path of a local bigWig \
                 (a species set or FTP address is not looked up)"
            )));
        }
        warn_unread_parameters(NAME, unread, "the bigWig and a closing MAX");
        self.method = if max_closes_the_list {
            Method::Max
        } else {
            Method::Average
        };
        self.source = Some(Source::open(file)?);
        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(source) = &self.source else {
            return Ok(None);
        };
        Ok(Some(PrefetchData::new(
            source.scores(variants, self.method)?,
        )))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        for (index, score) in data.downcast::<Scores>()? {
            variants[index].set_plugin_value(FIELD.to_string(), score, false, self.feature_types());
        }
        Ok(())
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        if let Some(data) = self.prefetch(variants)? {
            self.annotate(variants, data)?;
        }
        Ok(())
    }

    fn run(
        &self,
        _consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // Record-level plugin: the value is written in `annotate`.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::variant::VariantClass;

    /// `tests/data/cons.bw`: on `chr21` the one-base intervals at 1-based 100,
    /// 101, 102, 103 and 106 with the values 10, 20, 30, 0 and -1.5, the
    /// interval 200-210 with the value 7.25 and 301-302 with 0.1; on
    /// `CM000684.2` the base 50 with the value 5.
    fn fixture() -> String {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/cons.bw")
            .display()
            .to_string()
    }

    /// `tests/data/cons_negative_zero.bw`, made-up scores on the made-up
    /// chromosome `made_up`: the one-base intervals at 1-based 10, 20 and 21
    /// with the values -0, -0 and -3.5, each zero stored with its sign bit set.
    fn negative_zero_fixture() -> String {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/cons_negative_zero.bw")
            .display()
            .to_string()
    }

    /// One-base intervals at consecutive bases, one per value.
    fn one_base_intervals(values: &[f32]) -> impl Iterator<Item = Value> + '_ {
        values.iter().zip(0u32..).map(|(&value, start)| Value {
            start,
            end: start + 1,
            value,
        })
    }

    fn plugin(params: &[&str]) -> ConservationPlugin {
        let params: Vec<String> = params.iter().map(|p| p.to_string()).collect();
        let mut plugin = ConservationPlugin::new();
        plugin.init(&params).unwrap();
        plugin
    }

    fn refusal(params: &[&str]) -> String {
        let params: Vec<String> = params.iter().map(|p| p.to_string()).collect();
        ConservationPlugin::new()
            .init(&params)
            .unwrap_err()
            .to_string()
    }

    fn variant(chr: &str, start: u64, end: u64) -> InputVariant {
        InputVariant::new(chr.to_string(), start, end, b"A".to_vec(), b"G".to_vec())
    }

    fn score(plugin: &ConservationPlugin, mut variant: InputVariant) -> Option<String> {
        plugin
            .run_batch(std::slice::from_mut(&mut variant))
            .unwrap();
        variant.plugin_data.get(FIELD).cloned()
    }

    #[test]
    fn average_over_a_span() {
        let plugin = plugin(&[&fixture()]);
        // One base prints its value; several their average.
        assert_eq!(
            score(&plugin, variant("21", 106, 106)),
            Some("-1.500".into())
        );
        assert_eq!(
            score(&plugin, variant("21", 100, 102)),
            Some("20.000".into())
        );
        assert_eq!(
            score(&plugin, variant("21", 100, 103)),
            Some("15.000".into())
        );
        // A base without data contributes nothing: 103 (0) and 106 (-1.5) alone.
        assert_eq!(
            score(&plugin, variant("21", 103, 106)),
            Some("-0.750".into())
        );
        // An interval of several bases counts once per base it covers within
        // the span: two bases of the eleven-base 7.25 interval.
        assert_eq!(
            score(&plugin, variant("21", 199, 201)),
            Some("7.250".into())
        );
        // 102 (30), 103 (0) and 106 (-1.5) across two bases without data.
        assert_eq!(
            score(&plugin, variant("21", 102, 106)),
            Some("9.500".into())
        );
        // A span the file covers nowhere, and a chromosome it lacks, get nothing.
        assert_eq!(score(&plugin, variant("21", 400, 402)), None);
        assert_eq!(score(&plugin, variant("22", 100, 100)), None);
        // The record's `chr` spelling, in any case, reaches the file's `chr21`;
        // a name without the prefix reaches the file's name as written.
        assert_eq!(
            score(&plugin, variant("chr21", 101, 101)),
            Some("20.000".into())
        );
        assert_eq!(
            score(&plugin, variant("CHR21", 101, 101)),
            Some("20.000".into())
        );
        assert_eq!(
            score(&plugin, variant("CM000684.2", 50, 50)),
            Some("5.000".into())
        );
        // A structural variant is not handed to the plugin.
        let mut deletion = variant("21", 100, 103);
        deletion.variant_class = VariantClass::StructuralDeletion;
        assert_eq!(score(&plugin, deletion), None);
        assert_eq!(
            plugin.header_info(),
            [(
                FIELD.to_string(),
                "The average conservation score for this site".to_string()
            )]
        );
    }

    #[test]
    fn insertion_averages_the_two_flanks() {
        let plugin = plugin(&[&fixture()]);
        // Between 100 (10) and 101 (20).
        assert_eq!(
            score(&plugin, variant("21", 101, 100)),
            Some("15.000".into())
        );
        // Between 103 (0) and 104 (no data): the one covered flank alone.
        assert_eq!(
            score(&plugin, variant("21", 104, 103)),
            Some("0.000".into())
        );
        // Between two bases without data.
        assert_eq!(score(&plugin, variant("21", 105, 104)), None);
    }

    /// An allele of a `--minimal` split record averages its own trimmed bases,
    /// the ones its row's alleles describe, as the same allele does when it is
    /// the record's only alternate; the record's span, the one its `Location`
    /// prints, is not read. An insertion allele reads the two flanks of its own
    /// insertion point.
    #[test]
    fn a_split_allele_reads_its_own_span() {
        let average = plugin(&[&fixture()]);
        // `-` of `XYY/X/Z` at 100: own bases 101 (20) and 102 (30); the record
        // spans 100 (10) to 102.
        let mut deletion = variant("21", 101, 102);
        deletion.original_start = Some(100);
        deletion.original_end = Some(102);
        assert_eq!(window(&deletion), Some((100, 102)));
        assert_eq!(score(&average, deletion), Some("25.000".into()));
        // The insertion allele of `X/XY/Z` at 100: its own flanks 100 (10) and
        // 101 (20); the record's one base is 100.
        let mut insertion = variant("21", 101, 100);
        insertion.original_start = Some(100);
        insertion.original_end = Some(100);
        assert_eq!(window(&insertion), Some((99, 101)));
        assert_eq!(score(&average, insertion.clone()), Some("15.000".into()));
        assert_eq!(
            score(&average, variant("21", 101, 102)),
            Some("25.000".into())
        );
        let maximum = plugin(&[&fixture(), "MAX"]);
        assert_eq!(score(&maximum, insertion), Some("20.000".into()));
    }

    #[test]
    fn max_option_reports_the_maximum() {
        let plugin = plugin(&[&fixture(), "MAX"]);
        assert_eq!(
            score(&plugin, variant("21", 100, 103)),
            Some("30.000".into())
        );
        assert_eq!(
            score(&plugin, variant("21", 103, 106)),
            Some("0.000".into())
        );
        assert_eq!(
            score(&plugin, variant("21", 106, 106)),
            Some("-1.500".into())
        );
        assert_eq!(
            score(&plugin, variant("21", 101, 100)),
            Some("20.000".into())
        );
        assert_eq!(
            plugin.header_info(),
            [(
                FIELD.to_string(),
                "The maximum conservation score for this site".to_string()
            )]
        );
        // The method is read from the last parameter alone.
        assert_eq!(
            refusal(&["MAX", &fixture()]),
            "plugin initialization failed: Conservation: MAX must be the last parameter"
        );
    }

    #[test]
    fn database_mode_is_refused() {
        let database = refusal(&["database", "GERP_CONSERVATION_SCORE", "mammals"]);
        assert!(
            database.contains("Conservation: the database mode"),
            "{database}"
        );
        assert!(refusal(&["database"]).contains("Conservation: the database mode"));
        // The remote lookup of a species set is a file that is not there.
        assert!(refusal(&["mammals"]).contains("Conservation: bigWig file mammals not found"));
        assert!(refusal(&[]).contains("Conservation: no bigWig file given"));
        assert!(refusal(&["path=x"]).contains("Conservation: unknown parameter 'path'"));
        // A file that is not a bigWig is refused with its name.
        let dir = tempfile::tempdir().unwrap();
        let text = dir.path().join("scores.bw");
        std::fs::write(&text, b"not a bigWig\n").unwrap();
        let text = text.display().to_string();
        assert!(
            refusal(&[&text]).contains(&format!("Conservation: {text} is not a readable bigWig"))
        );
        // The keyed spelling of the file.
        let keyed = plugin(&[&format!("file={}", fixture())]);
        assert_eq!(
            score(&keyed, variant("21", 100, 100)),
            Some("10.000".into())
        );
    }

    /// A positional token after the bigWig is not read (`Conservation.pm:91`
    /// reads the first parameter alone), and `MAX` still closes the list after
    /// it; `MAX` followed by a token is the stated refusal, where the Ensembl
    /// VEP plugin drops the requested maximum without a word (`:82-89`).
    #[test]
    fn a_token_after_the_bigwig_is_not_read() {
        let average = plugin(&[&fixture(), "extra"]);
        assert_eq!(
            score(&average, variant("21", 100, 103)),
            Some("15.000".into())
        );
        let maximum = plugin(&[&fixture(), "extra", "MAX"]);
        assert_eq!(
            score(&maximum, variant("21", 100, 103)),
            Some("30.000".into())
        );
        assert_eq!(
            refusal(&[&fixture(), "MAX", "extra"]),
            "plugin initialization failed: Conservation: MAX must be the last parameter"
        );
        // A second `MAX` closing the list does not excuse the first: Ensembl
        // VEP's plugin reads that first token as the file (`:192`).
        assert_eq!(
            refusal(&["MAX", &fixture(), "MAX"]),
            "plugin initialization failed: Conservation: MAX must be the last parameter"
        );
    }

    /// One covered base prints its stored value (`Conservation.pm:252-254`), so
    /// a stored zero whose sign bit is set prints with its sign.
    #[test]
    fn one_covered_base_prints_a_stored_negative_zero_with_its_sign() {
        assert_eq!(
            super::score(one_base_intervals(&[-0.0]), Method::Average),
            Some("-0.000".into())
        );
    }

    /// Several covered bases are summed from zero and divided by their count
    /// (`Conservation.pm:246-249`), so a negative zero among them adds nothing.
    #[test]
    fn a_negative_zero_beside_another_value_is_averaged() {
        assert_eq!(
            super::score(one_base_intervals(&[-0.0, -3.5]), Method::Average),
            Some("-1.750".into())
        );
        assert_eq!(
            super::score(one_base_intervals(&[-3.5, -0.0]), Method::Average),
            Some("-1.750".into())
        );
    }

    /// The sum of several covered bases starts at zero (`Conservation.pm:246`),
    /// not at the first value, so two negative zeros average to a plain zero,
    /// as two one-base intervals and as one interval of two bases.
    #[test]
    fn two_negative_zeros_average_to_a_plain_zero() {
        assert_eq!(
            super::score(one_base_intervals(&[-0.0, -0.0]), Method::Average),
            Some("0.000".into())
        );
        let two_bases = Value {
            start: 0,
            end: 2,
            value: -0.0,
        };
        assert_eq!(
            super::score(std::iter::once(two_bases), Method::Average),
            Some("0.000".into())
        );
    }

    /// Under `MAX` the maximum of one covered base is that base's stored value.
    #[test]
    fn the_maximum_of_one_negative_zero_keeps_its_sign() {
        assert_eq!(
            super::score(one_base_intervals(&[-0.0]), Method::Max),
            Some("-0.000".into())
        );
    }

    /// Read from a bigWig, a stored negative zero prints with its sign where it
    /// is the one covered base and enters the average where a second base is
    /// covered beside it.
    #[test]
    fn a_negative_zero_read_from_the_file_keeps_its_sign_on_one_base() {
        let plugin = plugin(&[&negative_zero_fixture()]);
        assert_eq!(
            score(&plugin, variant("made_up", 10, 10)),
            Some("-0.000".into())
        );
        assert_eq!(
            score(&plugin, variant("made_up", 20, 21)),
            Some("-1.750".into())
        );
    }
}
