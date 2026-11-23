// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! LoFtool plugin: gene-level loss-of-function intolerance percentiles, as
//! Ensembl VEP's LoFtool plugin (LoFtool.pm) writes them.
//!
//! **Usage:** `--plugin LoFtool[,<scores file>]`
//!
//! Without a file, or with a file parameter Perl reads as false (`0` or the
//! empty string), the plugin reads `LoFtool_scores.txt` from the first plugins
//! directory that holds one, as the Ensembl VEP plugin reads the copy bundled
//! beside its module (LoFtool.pm:74-78). The file has two whitespace-separated
//! columns, a gene symbol and its percentile, under a header whose score cell
//! reads `LoFtool_percentile` (`:87-88`); the bundled file's lines end in CRLF.
//!
//! Each transcript row whose gene symbol has a score, matched
//! case-insensitively (`:89`, `:118`), carries it under `LoFtool` as Perl's
//! `%g` writes the number (`:89`): the file's `6.89E-05` prints `6.89e-05`, its
//! `1.38E-04` prints `0.000138`.
//!
//! Per-transcript plugin: returns annotation from [`crate::BuiltinPlugin::run`] via
//! the default [`crate::BuiltinPlugin::run_batch`] dispatcher.

use indexmap::IndexMap;
use tracing::debug;

use crate::perl_format::sprintf_g;
use crate::plugins::gene_scores::{bundled_file, MissingCell, ScoreTable};
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::consequence::TranscriptConsequence;
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

/// The output field (LoFtool.pm:107, 118).
const OUTPUT_FIELD: &str = "LoFtool";

/// The scores file read from the plugins directory when none is given
/// (LoFtool.pm:77).
const BUNDLED_FILE: &str = "LoFtool_scores.txt";

/// LoFtool gene intolerance score plugin.
#[derive(Default)]
pub struct LoFtoolPlugin {
    /// Lower-cased gene symbol to the percentile as it prints.
    scores: ScoreTable,
    /// The run's plugins directories, in search order.
    plugin_dirs: Vec<String>,
}

impl LoFtoolPlugin {
    /// Create a new uninitialized LoFtool plugin.
    pub fn new() -> Self {
        Self::default()
    }

    /// The row's score: the one filed under its gene symbol (LoFtool.pm:115-118).
    fn score_of(&self, consequence: &TranscriptConsequence) -> Option<&str> {
        self.scores.get(consequence.gene_symbol()?)
    }
}

impl crate::traits::sealed::Sealed for LoFtoolPlugin {}

impl BuiltinPlugin for LoFtoolPlugin {
    fn name(&self) -> &str {
        "LoFtool"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            OUTPUT_FIELD.to_string(),
            "LoFtool score for gene".to_string(),
        )]
    }

    fn set_plugin_dirs(&mut self, dirs: &[String]) {
        self.plugin_dirs = dirs.to_vec();
    }

    /// The first parameter is the scores file (LoFtool.pm:72); without one, or
    /// with one Perl reads as false (`0` or the empty string, `if(!$file)` at
    /// `:74`), the bundled file is read from the plugins directory.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let path = match params.first().filter(|file| perl_true(file)) {
            Some(file) => std::path::PathBuf::from(file),
            None => bundled_file(self.name(), "scores", &self.plugin_dirs, BUNDLED_FILE)?,
        };
        self.scores = ScoreTable::load(
            self.name(),
            "scores",
            &path,
            "LoFtool_percentile",
            MissingCell::Skipped,
            sprintf_g,
        )?;
        debug!(
            genes = self.scores.len(),
            path = %path.display(),
            "LoFtool plugin loaded scores"
        );
        Ok(())
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut result = IndexMap::new();
        if let Some(score) = self.score_of(consequence) {
            result.insert(OUTPUT_FIELD.to_string(), score.to_string());
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use vep_core::transcript::Transcript;

    fn test_variant() -> InputVariant {
        InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec())
    }

    fn row(symbol: &str) -> TranscriptConsequence {
        let mut consequence = TranscriptConsequence::default();
        consequence.transcript = Some(Arc::new({
            let mut transcript = Transcript::default();
            transcript.gene_symbol = Some(symbol.into());
            transcript
        }));
        consequence
    }

    /// A plugin loaded from `text` written to a file in a fresh directory, the
    /// directory kept alive with it.
    fn loaded(text: &[u8]) -> (LoFtoolPlugin, tempfile::TempDir) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("LoFtool_scores.txt");
        std::fs::write(&path, text).unwrap();
        let mut plugin = LoFtoolPlugin::new();
        plugin.init(&[path.display().to_string()]).unwrap();
        (plugin, dir)
    }

    fn score(plugin: &LoFtoolPlugin, symbol: &str) -> Option<String> {
        plugin
            .run(&row(symbol), &test_variant())
            .unwrap()
            .get(OUTPUT_FIELD)
            .cloned()
    }

    #[test]
    fn plugin_name_is_correct() {
        assert_eq!(LoFtoolPlugin::new().name(), "LoFtool");
    }

    /// The field is `LoFtool`, described as LoFtool.pm:107 describes it.
    #[test]
    fn header_info_field_is_loftool() {
        assert_eq!(
            LoFtoolPlugin::new().header_info(),
            [("LoFtool".to_string(), "LoFtool score for gene".to_string())]
        );
    }

    /// The scores are keyed by `lc($gene)` and looked up by `lc($symbol)`
    /// (LoFtool.pm:89, 118), so a symbol matches in any case.
    #[test]
    fn symbol_lookup_is_case_insensitive() {
        let (plugin, _dir) = loaded(b"Gene\tLoFtool_percentile\r\nSod1\t0.85\r\n");
        assert_eq!(score(&plugin, "SOD1").as_deref(), Some("0.85"));
        assert_eq!(score(&plugin, "sod1").as_deref(), Some("0.85"));
        assert_eq!(score(&plugin, "Sod1").as_deref(), Some("0.85"));
        assert_eq!(score(&plugin, "SOD2"), None);
    }

    /// A score prints as `sprintf("%g")` writes it (LoFtool.pm:89): the two
    /// spellings the bundled file carries.
    #[test]
    fn score_formats_like_perl_percent_g() {
        let (plugin, _dir) = loaded(
            b"Gene\tLoFtool_percentile\r\nNIPBL\t6.89E-05\r\nSCN1A\t1.38E-04\r\nRB1\t1.00E+00\r\n",
        );
        assert_eq!(score(&plugin, "NIPBL").as_deref(), Some("6.89e-05"));
        assert_eq!(score(&plugin, "SCN1A").as_deref(), Some("0.000138"));
        assert_eq!(score(&plugin, "RB1").as_deref(), Some("1"));
    }

    /// A cell is read as `sprintf("%g", $score)` (LoFtool.pm:89) numifies it:
    /// `1.5abc` scores `1.5`, a cell with no leading number (`NA`) scores
    /// nothing, and the later of two lines for one gene wins, so a later `NA`
    /// takes an earlier score away as the Ensembl VEP hash's `0` does under the
    /// truth test at `:118`.
    #[test]
    fn a_cell_is_read_as_perl_numifies_it() {
        let (plugin, _dir) = loaded(
            b"Gene\tLoFtool_percentile\r\nAPP\t1.5abc\r\nSOD1\tNA\r\nCXADR\t1.38E-04\r\n\
              CXADR\tNA\r\nSCAF4\tNA\r\nSCAF4\t0.85\r\n",
        );
        assert_eq!(score(&plugin, "APP").as_deref(), Some("1.5"));
        assert_eq!(score(&plugin, "SOD1"), None);
        assert_eq!(score(&plugin, "CXADR"), None);
        assert_eq!(score(&plugin, "SCAF4").as_deref(), Some("0.85"));
    }

    /// Without a file parameter the plugin reads `LoFtool_scores.txt` from the
    /// plugins directory (LoFtool.pm:74-78); a directory without it is refused
    /// naming the path looked for (`:80`).
    #[test]
    fn default_scores_file_resolves_from_plugins_dir() {
        let dir = tempfile::tempdir().unwrap();
        let dirs = [dir.path().display().to_string()];
        let mut plugin = LoFtoolPlugin::new();
        plugin.set_plugin_dirs(&dirs);
        let error = plugin.init(&[]).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "plugin initialization failed: LoFtool: scores file {} not found",
                dir.path().join("LoFtool_scores.txt").display()
            )
        );
        std::fs::write(
            dir.path().join("LoFtool_scores.txt"),
            b"Gene\tLoFtool_percentile\r\nSOD1\t0.85\r\n",
        )
        .unwrap();
        let mut plugin = LoFtoolPlugin::new();
        plugin.set_plugin_dirs(&dirs);
        plugin.init(&[]).unwrap();
        assert_eq!(score(&plugin, "SOD1").as_deref(), Some("0.85"));
    }

    /// The file parameter is read as Perl's `if(!$file)` reads it
    /// (LoFtool.pm:74): `0` and the empty string are no file, so the bundled
    /// file is read, while `00` and `0.0` are paths and refused as such.
    #[test]
    fn a_false_file_parameter_reads_the_bundled_file() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(
            dir.path().join("LoFtool_scores.txt"),
            b"Gene\tLoFtool_percentile\r\nSOD1\t0.85\r\n",
        )
        .unwrap();
        let dirs = [dir.path().display().to_string()];
        for absent in ["0", ""] {
            let mut plugin = LoFtoolPlugin::new();
            plugin.set_plugin_dirs(&dirs);
            plugin.init(&[absent.to_string()]).unwrap();
            assert_eq!(
                score(&plugin, "SOD1").as_deref(),
                Some("0.85"),
                "{absent:?}"
            );
        }
        for path in ["00", "0.0"] {
            let mut plugin = LoFtoolPlugin::new();
            plugin.set_plugin_dirs(&dirs);
            assert_eq!(
                plugin.init(&[path.to_string()]).unwrap_err().to_string(),
                format!("plugin initialization failed: LoFtool: scores file {path} not found")
            );
        }
    }

    /// Without a plugins directory the default file is refused by its bare name,
    /// never read from the working directory, and a named file that does not
    /// exist is refused by its path.
    #[test]
    fn missing_scores_file_is_refused_by_path() {
        let mut plugin = LoFtoolPlugin::new();
        assert_eq!(
            plugin.init(&[]).unwrap_err().to_string(),
            "plugin initialization failed: LoFtool: scores file LoFtool_scores.txt not found"
        );
        assert_eq!(
            plugin
                .init(&["/no/such/scores.txt".to_string()])
                .unwrap_err()
                .to_string(),
            "plugin initialization failed: LoFtool: scores file /no/such/scores.txt not found"
        );
    }

    /// A scores file whose third line is one byte over the cap is refused as
    /// the plugin opens, naming the file and the line, and the scores read
    /// before it are not kept.
    #[test]
    fn scores_line_past_the_cap_is_refused() {
        let cap = crate::tabix::MAX_LINE_BYTES;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("LoFtool_scores.txt");
        let mut text = b"Gene\tLoFtool_percentile\nSOD1\t0.85\n".to_vec();
        text.extend(std::iter::repeat_n(b'x', cap + 1));
        text.push(b'\n');
        std::fs::write(&path, text).unwrap();
        let mut plugin = LoFtoolPlugin::new();
        let error = plugin.init(&[path.display().to_string()]).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "plugin initialization failed: {}: line 3 is longer than the {cap}-byte limit \
                 for one line",
                path.display()
            )
        );
        assert!(plugin
            .run(&row("SOD1"), &test_variant())
            .unwrap()
            .is_empty());
    }

    /// A row without a gene symbol, or with one the file lacks, carries nothing
    /// (LoFtool.pm:116, 118).
    #[test]
    fn rows_without_a_scored_symbol_carry_nothing() {
        let (plugin, _dir) = loaded(b"Gene\tLoFtool_percentile\r\nSOD1\t0.85\r\n");
        let variant = test_variant();
        assert!(plugin
            .run(&TranscriptConsequence::default(), &variant)
            .unwrap()
            .is_empty());
        assert!(plugin
            .run(&row("UNKNOWN_GENE"), &variant)
            .unwrap()
            .is_empty());
    }

    /// `run_batch` writes the same entry `run` returns, on transcript rows alone.
    #[test]
    fn run_batch_writes_the_field_on_transcript_rows() {
        let (plugin, _dir) = loaded(b"Gene\tLoFtool_percentile\r\nSOD1\t0.85\r\n");
        let mut variants = vec![test_variant()];
        variants[0]
            .transcript_consequences
            .extend([row("SOD1"), row("SOD2"), {
                let mut built = row("SOD1");
                built.feature_type = vep_core::consequence::FeatureType::RegulatoryFeature;
                built
            }]);
        plugin.run_batch(&mut variants).unwrap();
        let written: Vec<Option<&String>> = variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.get(OUTPUT_FIELD))
            .collect();
        assert_eq!(written, [Some(&"0.85".to_string()), None, None]);
    }

    /// Ensembl VEP 116.3 never hands LoFtool a structural variant: the module
    /// keeps the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2229). The transcript rows of a
    /// `<DEL>` over a scored gene carry no field, while the sequence variant's
    /// row of the same gene does.
    #[test]
    fn a_structural_variant_row_gets_no_field() {
        let (plugin, _dir) = loaded(b"Gene\tLoFtool_percentile\r\nSOD1\t0.85\r\n");
        assert!(!plugin.annotates_structural_variants());
        let mut deletion =
            InputVariant::new("1".into(), 100, 5000, b"N".to_vec(), b"<DEL>".to_vec());
        deletion.is_structural = true;
        deletion.transcript_consequences.push(row("SOD1"));
        let mut sequence = test_variant();
        sequence.transcript_consequences.push(row("SOD1"));
        let mut variants = vec![deletion, sequence];
        plugin.run_batch(&mut variants).unwrap();
        assert!(variants[0].transcript_consequences[0]
            .plugin_data
            .is_empty());
        assert_eq!(
            variants[1].transcript_consequences[0]
                .plugin_data
                .get(OUTPUT_FIELD)
                .map(String::as_str),
            Some("0.85")
        );
    }
}
