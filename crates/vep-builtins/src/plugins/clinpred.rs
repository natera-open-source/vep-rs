// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! ClinPred plugin: pathogenicity scores for nonsynonymous single-base
//! substitutions, read from the tabix-indexed score table prepared as
//! `ClinPred.pm:36-44` prepares it (space to tab, the first line commented
//! and `Chr` lowered to `chr`, `tabix -s 1 -b 2 -e 2`).
//!
//! # Parameters
//!
//! The score file is `file=<path>`, or the first positional parameter when no
//! parameter is `key=value` (`ClinPred.pm:79-86`). The file must exist; Ensembl
//! VEP reads the run's assembly and does nothing with it (`:90`), so no build
//! marker is looked for here either.
//!
//! # Data file
//!
//! Five positional columns, `chr start ref alt score` (`:127`); the recipe
//! comments the header line (`:41`), which the index skips. The reference
//! column is parsed and not compared
//! (`:116-120` match the alternate allele and the coordinates alone), so a row
//! whose reference disagrees with the record's still scores, as it does there.
//! A matched row lacking its fifth cell scores as the empty value and a present
//! empty fifth cell as the empty string, the two values `:127`'s five-variable
//! split hands the writers.
//!
//! # Scoring
//!
//! A transcript row scores when its consequence terms include `missense_variant`,
//! `stop_lost`, `stop_gained` or `start_lost` (`:69`, `:106`) and the alternate
//! allele is one base (`:110`); the plugin declares `['Transcript']` (`:96`), so
//! the score is a per-row value. The record's alternate allele is compared on the
//! forward strand with the row's `alt` cell as written (`:109`, `:117` change neither).
//!
//! # Output fields
//!
//! - `ClinPred`: the score of the matching row (`:101`, `:133`).
//!
//! `ClinPred.pm` is the plugin's Ensembl VEP counterpart in Ensembl/VEP_plugins
//! release/116.

use std::path::Path;

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{open_best_store_with, AnnotationStore};
use crate::params::PluginParams;
use crate::tabix::{
    allele_match, table_contig, BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord,
};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::assembly::Assembly;
use vep_core::consequence::{Consequence, TranscriptConsequence};
use vep_core::variant::InputVariant;

const NAME: &str = "ClinPred";

/// The output field.
const FIELD: &str = "ClinPred";

/// The consequence terms a row scores under (`ClinPred.pm:69`).
const GATE_TERMS: [Consequence; 4] = [
    Consequence::MissenseVariant,
    Consequence::StopLost,
    Consequence::StopGained,
    Consequence::StartLost,
];

/// 0-based positional columns (`ClinPred.pm:127`).
const CHR_COL: usize = 0;
const START_COL: usize = 1;
const REF_COL: usize = 2;
const ALT_COL: usize = 3;
const SCORE_COL: usize = 4;

/// ClinPred score plugin.
pub struct ClinPredPlugin {
    /// Annotation store (tabix or binary), opened by `init`.
    store: Option<Box<dyn AnnotationStore>>,
    /// The run's assembly, which names a mitochondrial record's contig to the
    /// table ([`table_contig`]).
    assembly: Option<Assembly>,
}

impl Default for ClinPredPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl ClinPredPlugin {
    /// Create an uninitialized ClinPred plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            assembly: None,
        }
    }

    /// Whether a row is one the plugin scores: a transcript row (`ClinPred.pm:96`)
    /// whose consequence terms include a gate term (`:106`).
    fn scores_row(&self, consequence: &TranscriptConsequence) -> bool {
        self.annotates_row(consequence)
            && consequence
                .consequences
                .iter()
                .any(|term| GATE_TERMS.contains(term))
    }

    /// The forward-strand alternate allele of a record the file can score: one
    /// base of `ACGT` (`ClinPred.pm:110`) at a single position (`:118-119`).
    fn scorable_alt(variant: &InputVariant) -> Option<u8> {
        if variant.start != variant.end {
            return None;
        }
        let (_, alt) = variant.forward_strand_alleles();
        match alt.as_ref() {
            [base @ (b'A' | b'C' | b'G' | b'T')] => Some(*base),
            _ => None,
        }
    }

    /// The score cell of the first of `records` at `position` carrying the
    /// alternate allele `alt` (`ClinPred.pm:116-122`), the row's position its
    /// Start cell as `$_->{start} == $vf_start` (`:118`) numifies it
    /// ([`allele_match::perl_position`]: `30396091.0` and `30396091x` are
    /// 30396091, `30396091.5` is no position): `None` when no row matches,
    /// `Some(None)` for a matched row lacking its fifth cell, which `:127`'s
    /// five-variable split leaves undefined, and `Some(Some(cell))` for the cell
    /// as the line holds it, the empty string for a trailing empty cell, which
    /// that split keeps.
    fn matching_score<'a>(
        records: impl IntoIterator<Item = &'a TabixRecord>,
        position: u64,
        alt: u8,
    ) -> Option<Option<String>> {
        records
            .into_iter()
            .filter(|record| allele_match::perl_position(record, START_COL) == Some(position))
            .find(|record| {
                record
                    .columns
                    .get(ALT_COL)
                    .is_some_and(|file_alt| file_alt.as_bytes() == [alt])
            })
            .map(|record| record.columns.get(SCORE_COL).cloned())
    }
}

impl crate::traits::sealed::Sealed for ClinPredPlugin {}

impl BuiltinPlugin for ClinPredPlugin {
    fn name(&self) -> &str {
        NAME
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            FIELD.to_string(),
            "Prediction tool to identify disease-relevant nonsynonymous single-nucleotide variants"
                .to_string(),
        )]
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn set_reference_fasta(&mut self, fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>) {
        if let Some(store) = self.store.as_mut() {
            store.set_reference_fasta(fasta);
        }
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(NAME, &["file"], params)?;
        let file = params
            .get("file")
            .or_else(|| params.positional().first().map(String::as_str))
            .filter(|file| !file.is_empty())
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: no file specified (--plugin {NAME},file=<path>)"
                ))
            })?;
        let path = Path::new(file);
        if !path.is_file() {
            return Err(PluginError::Init(format!("{NAME}: file {file} not found")));
        }
        let config = TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: CHR_COL,
            start_col: START_COL,
            end_col: None,
            ref_col: Some(REF_COL),
            alt_col: Some(ALT_COL),
            zero_based: false,
        };
        self.store = Some(open_best_store_with(
            path,
            config,
            self.assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )?);
        debug!(file, "ClinPred plugin initialized");
        Ok(())
    }

    fn prefetch_reads_consequences(&self) -> bool {
        true
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|v| {
                self.annotates_variant(v)
                    && Self::scorable_alt(v).is_some()
                    && v.transcript_consequences
                        .iter()
                        .any(|tc| self.scores_row(tc))
            })
            .map(|v| (table_contig(v, self.assembly), v.start, v.start))
            .collect();
        if regions.is_empty() {
            return Ok(Some(PrefetchData::new(BatchQueryResult::empty())));
        }
        Ok(Some(PrefetchData::new(store.query_batch(&regions)?)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch: BatchQueryResult = data.downcast()?;
        if batch.is_empty() {
            return Ok(());
        }
        for variant in variants.iter_mut().filter(|v| self.annotates_variant(v)) {
            let Some(alt) = Self::scorable_alt(variant) else {
                continue;
            };
            let position = variant.start;
            let Some(cell) = Self::matching_score(
                batch.lookup(table_contig(variant, self.assembly), position, position),
                position,
                alt,
            ) else {
                continue;
            };
            // A matched row lacking its score cell scores as Ensembl VEP's
            // undefined value, an unmarked empty value here; a present empty
            // cell is its empty string, which the row names (`ClinPred.pm:127`).
            let is_empty_string = cell.as_ref().is_some_and(String::is_empty);
            let score = cell.unwrap_or_default();
            for tc in variant.transcript_consequences.iter_mut() {
                if self.scores_row(tc) {
                    tc.plugin_data.insert(FIELD.to_string(), score.clone());
                    // A later plugin's value replaces an earlier one under the
                    // same key (OutputFactory.pm:2239-2240), its mark too.
                    tc.plugin_empty_strings.retain(|key| key != FIELD);
                    if is_empty_string {
                        tc.plugin_empty_strings.push(FIELD.to_string());
                    }
                }
            }
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
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let (fields, _) = self.run_naming_empty_strings(consequence, variant)?;
        Ok(fields)
    }

    /// The row's score from one query of its position, with the field named
    /// where the matched row's score cell is present and empty, the empty
    /// string `annotate` names in the row's `plugin_empty_strings`; a matched
    /// row lacking the cell scores as the unnamed empty value
    /// (`ClinPred.pm:127`).
    fn run_naming_empty_strings(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<(IndexMap<String, String>, Vec<String>), PluginError> {
        let mut fields = IndexMap::new();
        let mut empty_strings = Vec::new();
        let Some(store) = &self.store else {
            return Ok((fields, empty_strings));
        };
        if !self.scores_row(consequence) {
            return Ok((fields, empty_strings));
        }
        let Some(alt) = Self::scorable_alt(variant) else {
            return Ok((fields, empty_strings));
        };
        let records = store.query(
            table_contig(variant, self.assembly),
            variant.start,
            variant.start,
        )?;
        if let Some(cell) = Self::matching_score(&records, variant.start, alt) {
            if cell.as_ref().is_some_and(String::is_empty) {
                empty_strings.push(FIELD.to_string());
            }
            fields.insert(FIELD.to_string(), cell.unwrap_or_default());
        }
        Ok((fields, empty_strings))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use smallvec::smallvec;
    use vep_core::consequence::FeatureType;

    /// The heading line of this module's made-up score files.
    const HEADER: &str = "#chr\tStart\tRef\tAlt\tClinPred_Score";

    fn fixture(dir: &Path, name: &str, header: &str, rows: &[&str]) -> String {
        let mut lines = vec![header];
        lines.extend(rows);
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
            .display()
            .to_string()
    }

    fn loaded(params: &[String]) -> ClinPredPlugin {
        let mut plugin = ClinPredPlugin::new();
        plugin.init(params).expect("fixture loads");
        plugin
    }

    fn row_of(term: Consequence) -> TranscriptConsequence {
        let mut consequence = TranscriptConsequence::default();
        consequence.consequences = smallvec![term];
        consequence
    }

    /// A single-base record at `pos` with one transcript row per term.
    fn snv(pos: u64, reference: &[u8], alt: &[u8], terms: &[Consequence]) -> InputVariant {
        let mut variant =
            InputVariant::new("21".into(), pos, pos, reference.to_vec(), alt.to_vec());
        variant
            .transcript_consequences
            .extend(terms.iter().map(|term| row_of(*term)));
        variant
    }

    /// The `ClinPred` value of every row of a record, in row order.
    fn scores(plugin: &ClinPredPlugin, variant: InputVariant) -> Vec<Option<String>> {
        let mut variants = vec![variant];
        plugin.run_batch(&mut variants).unwrap();
        assert!(
            variants[0].plugin_data.is_empty(),
            "ClinPred writes per row"
        );
        variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.get(FIELD).cloned())
            .collect()
    }

    fn init_message(result: Result<(), PluginError>) -> String {
        match result {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected PluginError::Init, got {other:?}"),
        }
    }

    /// `ClinPred.pm:79-86`: `file=` when any parameter is keyed, else the first
    /// positional parameter.
    #[test]
    fn file_is_keyed_or_the_first_positional_parameter() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "synthetic_fixture.tsv.gz",
            HEADER,
            &["21\t5030000\tA\tG\t0.987654"],
        );
        let missense = || snv(5030000, b"A", b"G", &[Consequence::MissenseVariant]);
        let keyed = loaded(&[format!("file={path}")]);
        assert_eq!(scores(&keyed, missense()), [Some("0.987654".to_string())]);
        let positional = loaded(std::slice::from_ref(&path));
        assert_eq!(
            scores(&positional, missense()),
            [Some("0.987654".to_string())]
        );
        let message =
            init_message(ClinPredPlugin::new().init(&[format!("file={path}"), "x".into()]));
        assert_eq!(
            message,
            format!(
                "ClinPred: parameter 'x' is not key=value while 'file={path}' is; \
                 give every parameter as key=value, or every one by position"
            )
        );
    }

    /// `ClinPred.pm:88` hands the file to `add_file`, whose `check_file` dies on
    /// no file and on a path that does not exist.
    #[test]
    fn a_missing_file_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        assert_eq!(
            init_message(ClinPredPlugin::new().init(&[])),
            "ClinPred: no file specified (--plugin ClinPred,file=<path>)"
        );
        assert_eq!(
            init_message(ClinPredPlugin::new().init(&["file=".to_string()])),
            "ClinPred: no file specified (--plugin ClinPred,file=<path>)"
        );
        let missing = format!("{}/no_such_synthetic.tsv.gz", dir.path().display());
        assert_eq!(
            init_message(ClinPredPlugin::new().init(&[format!("file={missing}")])),
            format!("ClinPred: file {missing} not found")
        );
        assert_eq!(
            init_message(ClinPredPlugin::new().init(std::slice::from_ref(&missing))),
            format!("ClinPred: file {missing} not found")
        );
    }

    /// `ClinPred.pm:106`: a row scores when any of its terms is one of the four;
    /// a row of other terms, and a regulatory row, get nothing, and the score is
    /// never written on the record.
    #[test]
    fn only_the_four_gate_terms_are_scored() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "gate_synthetic.tsv.gz",
            HEADER,
            &["21\t5030000\tA\tG\t0.5"],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let mut variant = snv(
            5030000,
            b"A",
            b"G",
            &[
                Consequence::MissenseVariant,
                Consequence::StopLost,
                Consequence::StopGained,
                Consequence::StartLost,
                Consequence::SynonymousVariant,
                Consequence::IntronVariant,
            ],
        );
        variant.transcript_consequences.push({
            let mut consequence = TranscriptConsequence::default();
            consequence.consequences = smallvec![
                Consequence::SpliceRegionVariant,
                Consequence::MissenseVariant
            ];
            consequence
        });
        variant.transcript_consequences.push({
            let mut consequence = TranscriptConsequence::default();
            consequence.consequences = smallvec![Consequence::MissenseVariant];
            consequence.feature_type = FeatureType::RegulatoryFeature;
            consequence
        });
        let hit = Some("0.5".to_string());
        assert_eq!(
            scores(&plugin, variant),
            [
                hit.clone(),
                hit.clone(),
                hit.clone(),
                hit.clone(),
                None,
                None,
                hit,
                None
            ]
        );
    }

    /// `ClinPred.pm:110` scores an allele of one base alone and `:118-119`
    /// compare both coordinates with a row's one position.
    #[test]
    fn only_single_base_alleles_are_scored() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "shapes_synthetic.tsv.gz",
            HEADER,
            &[
                "21\t5030000\tA\tG\t0.5",
                "21\t5030001\tT\t-\t0.6",
                "21\t5030002\t-\tC\t0.7",
                "21\t5030003\tAG\tCT\t0.8",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let gate = [Consequence::MissenseVariant];
        assert_eq!(
            scores(&plugin, snv(5030000, b"A", b"G", &gate)),
            [Some("0.5".to_string())]
        );
        let mut deletion =
            InputVariant::new("21".into(), 5030001, 5030001, b"T".to_vec(), b"-".to_vec());
        deletion
            .transcript_consequences
            .push(row_of(Consequence::StopLost));
        assert_eq!(scores(&plugin, deletion), [None]);
        let mut insertion =
            InputVariant::new("21".into(), 5030003, 5030002, b"-".to_vec(), b"C".to_vec());
        insertion
            .transcript_consequences
            .push(row_of(Consequence::StopGained));
        assert_eq!(scores(&plugin, insertion), [None]);
        let mut substitution = InputVariant::new(
            "21".into(),
            5030003,
            5030004,
            b"AG".to_vec(),
            b"CT".to_vec(),
        );
        substitution
            .transcript_consequences
            .push(row_of(Consequence::MissenseVariant));
        assert_eq!(scores(&plugin, substitution), [None]);
    }

    /// `ClinPred.pm:127` splits a row into five columns by position and never
    /// reads the header, so the names on the header line do not matter, and the
    /// reference column is parsed and not compared (`:116-120`).
    #[test]
    fn columns_are_positional_chr_start_ref_alt_score() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "names_synthetic.tsv.gz",
            "#Chromosome\tPosition\tReference\tAlternate\tScore",
            &[
                "21\t5030000\tA\tC\t0.111111",
                "21\t5030000\tA\tG\t0.222222",
                "21\t5030000\tA\tT\t0.333333",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let gate = [Consequence::MissenseVariant];
        assert_eq!(
            scores(&plugin, snv(5030000, b"A", b"T", &gate)),
            [Some("0.333333".to_string())]
        );
        assert_eq!(
            scores(&plugin, snv(5030000, b"C", b"G", &gate)),
            [Some("0.222222".to_string())],
            "the reference column is not compared"
        );
        assert_eq!(scores(&plugin, snv(5030000, b"A", b"A", &gate)), [None]);
        assert_eq!(scores(&plugin, snv(5030001, b"A", b"G", &gate)), [None]);
    }

    /// The commented first line of the prepared file (`ClinPred.pm:41`) is a
    /// header the index skips: it yields no record and no score, and its names
    /// are read as the file's header.
    #[test]
    fn header_line_is_skipped() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "header_synthetic.tsv.gz",
            HEADER,
            &["21\t5030000\tA\tG\t0.5"],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let store = plugin.store.as_ref().unwrap();
        assert_eq!(
            store.header(),
            Some(
                ["chr", "Start", "Ref", "Alt", "ClinPred_Score"]
                    .map(String::from)
                    .as_slice()
            )
        );
        let records = store.query("21", 1, 6_000_000).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].columns[ALT_COL], "G");
        let gate = [Consequence::MissenseVariant];
        assert_eq!(
            scores(&plugin, snv(5030000, b"A", b"G", &gate)),
            [Some("0.5".to_string())]
        );
    }

    /// Off GRCh37 the record's chromosome is the cache's name, `MT`, and the
    /// store tries the four spellings of the mitochondrial contig
    /// (`chromosome_query_candidates`), so a table naming it `MT`, `chrMT`, `M`
    /// or `chrM` scores the record. Ensembl VEP 116.3's base class tries the
    /// record's spelling and its `chr` variant alone
    /// (`BaseVepTabixPlugin.pm:693-708`), the divergence
    /// `plugin_contig_mitochondrial_spelling` states.
    #[test]
    fn a_mitochondrial_record_reads_a_file_under_any_of_the_four_contig_spellings() {
        let dir = tempfile::tempdir().unwrap();
        for contig in ["MT", "chrMT", "M", "chrM"] {
            let path = fixture(
                dir.path(),
                &format!("clinpred_{contig}_synthetic.tsv.gz"),
                HEADER,
                &[&format!("{contig}\t3460\tG\tA\t0.5")],
            );
            let plugin = loaded(&[format!("file={path}")]);
            let mut record =
                InputVariant::new("MT".into(), 3460, 3460, b"G".to_vec(), b"A".to_vec());
            record
                .transcript_consequences
                .push(row_of(Consequence::MissenseVariant));
            assert_eq!(
                scores(&plugin, record),
                [Some("0.5".to_string())],
                "table contig {contig}"
            );
        }
        let path = fixture(
            dir.path(),
            "clinpred_chrM_against_21_synthetic.tsv.gz",
            HEADER,
            &["chrM\t3460\tG\tA\t0.5"],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let gate = [Consequence::MissenseVariant];
        assert_eq!(scores(&plugin, snv(3460, b"G", b"A", &gate)), [None]);
    }

    /// `ClinPred.pm:127` splits the line into five variables, so Ensembl VEP keeps a
    /// trailing empty cell as the empty string and leaves a missing cell
    /// undefined, and `:116-122` returns the matched row's result either way.
    /// Both print as `ClinPred=` in the default format; the tab format prints
    /// an empty cell against `-` and JSON `""` against no key (`Tab.pm:128`,
    /// `JSON.pm:327`), so the row names the field of a present empty cell in
    /// `plugin_empty_strings` and not that of a missing one.
    #[test]
    fn a_row_lacking_its_score_cell_is_the_empty_value_and_an_empty_cell_the_empty_string() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "cells_synthetic.tsv.gz",
            HEADER,
            &[
                "21\t5030000\tA\tG",
                "21\t5030001\tA\tG\t",
                "21\t5030002\tA\tG\t\tx",
                "21\t5030003\tA\tG\t0.5",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let gate = [Consequence::MissenseVariant];
        let empty = Some(String::new());
        let marks = |pos: u64| {
            let mut variants = vec![snv(pos, b"A", b"G", &gate)];
            plugin.run_batch(&mut variants).unwrap();
            let row = &variants[0].transcript_consequences[0];
            (
                row.plugin_data.get(FIELD).cloned(),
                row.plugin_empty_strings.clone(),
            )
        };
        assert_eq!(marks(5030000), (empty.clone(), vec![]), "four cells");
        assert_eq!(
            marks(5030001),
            (empty.clone(), vec![FIELD.to_string()]),
            "five cells, the fifth empty"
        );
        assert_eq!(
            marks(5030002),
            (empty, vec![FIELD.to_string()]),
            "six cells, the fifth empty"
        );
        assert_eq!(marks(5030003), (Some("0.5".to_string()), vec![]));
    }

    /// The per-row route reads a row's score cell as the batch route does.
    /// `run` hands the score back and `run_naming_empty_strings` names the
    /// field of a present empty cell beside it, so the trait's `run_batch`,
    /// which writes each row from that call, leaves the value and the mark
    /// the plugin's own batch route leaves: the unmarked empty value for a
    /// missing cell and the named empty string for an empty one.
    #[test]
    fn the_per_row_route_names_an_empty_score_cell_as_the_batch_route_does() {
        /// The plugin through its per-row methods alone, each row written by
        /// the trait's `run_batch`.
        struct PerRow<'a>(&'a ClinPredPlugin);

        impl crate::traits::sealed::Sealed for PerRow<'_> {}

        impl BuiltinPlugin for PerRow<'_> {
            fn name(&self) -> &str {
                NAME
            }
            fn header_info(&self) -> Vec<(String, String)> {
                self.0.header_info()
            }
            fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
                Ok(())
            }
            fn run(
                &self,
                consequence: &TranscriptConsequence,
                variant: &InputVariant,
            ) -> Result<IndexMap<String, String>, PluginError> {
                self.0.run(consequence, variant)
            }
            fn run_naming_empty_strings(
                &self,
                consequence: &TranscriptConsequence,
                variant: &InputVariant,
            ) -> Result<(IndexMap<String, String>, Vec<String>), PluginError> {
                self.0.run_naming_empty_strings(consequence, variant)
            }
        }

        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "per_row_cells_synthetic.tsv.gz",
            HEADER,
            &[
                "21\t5030000\tA\tG",
                "21\t5030001\tA\tG\t",
                "21\t5030002\tA\tG\t\tx",
                "21\t5030003\tA\tG\t0.5",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let gate = [Consequence::MissenseVariant];
        let left_by = |route: &dyn BuiltinPlugin, pos: u64| {
            let mut variants = vec![snv(pos, b"A", b"G", &gate)];
            route.run_batch(&mut variants).unwrap();
            let row = &variants[0].transcript_consequences[0];
            (
                row.plugin_data.get(FIELD).cloned(),
                row.plugin_empty_strings.clone(),
            )
        };
        let empty = Some(String::new());
        let named = vec![FIELD.to_string()];
        for (pos, cells, expected) in [
            (5030000, "four cells", (empty.clone(), vec![])),
            (
                5030001,
                "five cells, the fifth empty",
                (empty.clone(), named.clone()),
            ),
            (5030002, "six cells, the fifth empty", (empty, named)),
            (5030003, "a score", (Some("0.5".to_string()), vec![])),
        ] {
            assert_eq!(left_by(&PerRow(&plugin), pos), expected, "per row: {cells}");
            assert_eq!(left_by(&plugin, pos), expected, "batch: {cells}");
            let record = snv(pos, b"A", b"G", &gate);
            let fields = plugin
                .run(&record.transcript_consequences[0], &record)
                .unwrap();
            assert_eq!(fields.get(FIELD).cloned(), expected.0, "run: {cells}");
        }
    }

    /// `ClinPred.pm:127` splits the line htslib hands over, the line end alone
    /// removed, with `split("\t", $line)`, so the score cell keeps the
    /// whitespace before the line end and a whitespace-only cell is the score:
    /// `0.91 ` prints as `ClinPred=0.91 ` and ` ` as `ClinPred= `.
    #[test]
    fn a_trailing_space_on_the_last_cell_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "trailing_synthetic.tsv.gz",
            HEADER,
            &["21\t5030000\tA\tG\t0.91 ", "21\t5030001\tA\tG\t "],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let gate = [Consequence::MissenseVariant];
        assert_eq!(
            scores(&plugin, snv(5030000, b"A", b"G", &gate)),
            [Some("0.91 ".to_string())]
        );
        assert_eq!(
            scores(&plugin, snv(5030001, b"A", b"G", &gate)),
            [Some(" ".to_string())]
        );
    }

    /// Ensembl VEP 116.3 never hands ClinPred a structural variant: the module
    /// keeps the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2229). A `<DEL>` row gets no score,
    /// and nor does a scored substitution carrying the structural mark, which
    /// pins the gate rather than the allele test; the plain one is scored.
    #[test]
    fn a_structural_variant_row_gets_no_score() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "synthetic_fixture.tsv.gz",
            HEADER,
            &["21\t5030000\tA\tG\t0.987654"],
        );
        let plugin = loaded(&[format!("file={path}")]);
        assert!(!plugin.annotates_structural_variants());
        let gate = [Consequence::MissenseVariant];
        let mut deletion = snv(5030000, b"N", b"<DEL>", &gate);
        deletion.end = 5030100;
        deletion.is_structural = true;
        assert_eq!(scores(&plugin, deletion), [None]);
        let mut flagged = snv(5030000, b"A", b"G", &gate);
        flagged.is_structural = true;
        assert_eq!(scores(&plugin, flagged), [None]);
        assert_eq!(
            scores(&plugin, snv(5030000, b"A", b"G", &gate)),
            [Some("0.987654".to_string())]
        );
    }

    /// A row's Start cell is read as `$_->{start} == $vf_start` (ClinPred.pm:118)
    /// numifies it: `5030000.0`, `5030000x` and `5.03e6` are 5030000 and the
    /// row scores the variant there; `5030000.5`, `abc` and `0` are no
    /// position. The tabix index places a cell at its leading integer, so the
    /// first three file shapes reach the match and the exponent and the rest
    /// are read from made-up records built here.
    #[test]
    fn a_position_cell_is_matched_as_perl_numifies_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "cells_synthetic.tsv.gz",
            HEADER,
            &[
                "21\t5030000.0\tA\tG\t0.5",
                "21\t5030000.5\tA\tT\t0.6",
                "21\t5030000x\tA\tC\t0.7",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let gate = [Consequence::MissenseVariant];
        assert_eq!(
            scores(&plugin, snv(5030000, b"A", b"G", &gate)),
            [Some("0.5".to_string())]
        );
        assert_eq!(scores(&plugin, snv(5030000, b"A", b"T", &gate)), [None]);
        assert_eq!(
            scores(&plugin, snv(5030000, b"A", b"C", &gate)),
            [Some("0.7".to_string())]
        );

        let record = |start: &str| TabixRecord {
            chr: "21".to_string(),
            start: 5030000,
            end: 5030000,
            columns: ["21", start, "A", "G", "0.9"]
                .iter()
                .map(|cell| cell.to_string())
                .collect(),
            header: None,
        };
        let scored =
            |start: &str| ClinPredPlugin::matching_score([&record(start)], 5030000, b'G').is_some();
        assert!(scored("5.03e6"));
        assert!(scored(" 5030000"));
        assert!(!scored("abc"));
        assert!(!scored("0"));
        assert!(!scored("5030000.5"));
    }

    /// The row scores of a missense record spelt `record` at the site of a
    /// one-row table naming the mitochondrial contig `contig`, opened on
    /// `assembly`.
    fn mitochondrial_scores(
        dir: &Path,
        assembly: Option<vep_core::assembly::Assembly>,
        record: &str,
        contig: &str,
    ) -> Vec<Option<String>> {
        let path = fixture(
            dir,
            &format!("ClinPred_{contig}_synthetic.tsv.gz"),
            HEADER,
            &[&format!("{contig}\t3460\tG\tA\t0.5")],
        );
        let mut plugin = ClinPredPlugin::new();
        plugin.set_assembly(assembly);
        plugin
            .init(&[format!("file={path}")])
            .expect("fixture loads");
        let mut record = InputVariant::new(record.into(), 3460, 3460, b"G".to_vec(), b"A".to_vec());
        record
            .transcript_consequences
            .push(row_of(Consequence::MissenseVariant));
        scores(&plugin, record)
    }

    /// On GRCh37 `MT` and `chrMT` name the assembly's sequence and `chrM` and
    /// `M` another one. Without a reference sequence the table's rows cannot
    /// be compared with it, so a mitochondrial record reads the table under
    /// its own spelling and its `chr` variant alone, the name `ClinPred.pm:120`
    /// hands `BaseVepTabixPlugin.pm:693-708`.
    #[test]
    fn on_grch37_a_mitochondrial_record_reads_the_table_within_its_pair_of_names() {
        let dir = tempfile::tempdir().unwrap();
        let grch37 = Some(vep_core::assembly::Assembly::Grch37);
        for (record, read) in [
            ("MT", ["MT", "chrMT"]),
            ("chrMT", ["chrMT", "MT"]),
            ("chrM", ["chrM", "M"]),
            ("M", ["M", "chrM"]),
        ] {
            for contig in ["M", "MT", "chrM", "chrMT"] {
                let expected = read.contains(&contig).then(|| "0.5".to_string());
                assert_eq!(
                    mitochondrial_scores(dir.path(), grch37, record, contig),
                    [expected],
                    "record {record}, table {contig}"
                );
            }
        }
    }

    /// On GRCh38 the four names are one sequence: every spelling of the record
    /// reads every spelling of the table, as a run without an assembly does.
    #[test]
    fn off_grch37_a_mitochondrial_record_reads_every_spelling_of_the_table() {
        let dir = tempfile::tempdir().unwrap();
        for assembly in [Some(vep_core::assembly::Assembly::Grch38), None] {
            for record in ["MT", "chrMT", "chrM", "M"] {
                for contig in ["M", "MT", "chrM", "chrMT"] {
                    assert_eq!(
                        mitochondrial_scores(dir.path(), assembly, record, contig),
                        [Some("0.5".to_string())],
                        "record {record}, table {contig}, {assembly:?}"
                    );
                }
            }
        }
    }

    /// On GRCh37, with the run's reference sequence, a table naming its
    /// mitochondrial contig under the other pair of names than the record's
    /// is read when its rows state the reference bases, and not when they
    /// state the bases one position on: there the row at the record's site
    /// carries the record's alleles, the reference repeating its base, and
    /// its score is the neighbouring base's.
    #[test]
    fn on_grch37_a_table_under_the_other_pair_is_read_only_when_its_rows_state_the_reference() {
        use crate::test_fixtures::{
            made_up_alternate, made_up_reference, made_up_rows, repeated_base,
        };
        let dir = tempfile::tempdir().unwrap();
        let (reference, bases) = made_up_reference(dir.path(), 1600);
        let site = repeated_base(&bases, 100);
        let base = bases[site - 1];
        for (frame, offset, expected) in [
            ("stating_the_reference", 0, Some(format!("0.{site}"))),
            ("one_position_on", 1, None),
        ] {
            let rows = made_up_rows(
                "chrM",
                &bases,
                10..1590,
                offset,
                |contig, position, reference, alt| {
                    format!("{contig}\t{position}\t{reference}\t{alt}\t0.{position}")
                },
            );
            let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
            let path = fixture(
                dir.path(),
                &format!("ClinPred_{frame}_synthetic.tsv.gz"),
                HEADER,
                &rows,
            );
            let mut plugin = ClinPredPlugin::new();
            plugin.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
            plugin
                .init(&[format!("file={path}")])
                .expect("fixture loads");
            plugin.set_reference_fasta(Some(std::sync::Arc::clone(&reference)));
            let record_spelt = |spelling: &str| {
                let mut record = InputVariant::new(
                    spelling.into(),
                    site as u64,
                    site as u64,
                    vec![base],
                    vec![made_up_alternate(base)],
                );
                record
                    .transcript_consequences
                    .push(row_of(Consequence::MissenseVariant));
                record
            };
            for record in ["MT", "chrMT"] {
                assert_eq!(
                    scores(&plugin, record_spelt(record)),
                    std::slice::from_ref(&expected),
                    "{frame}, record {record}"
                );
            }
            assert_eq!(
                scores(&plugin, record_spelt("chrM")),
                [Some(format!("0.{site}"))],
                "{frame}: a record of the table's own pair"
            );
        }
    }
}
