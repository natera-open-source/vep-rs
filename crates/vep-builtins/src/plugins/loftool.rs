// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! LoFtool plugin: flat-file gene-level loss-of-function intolerance scores.
//!
//! Loads a simple TSV file mapping gene symbols to LoFtool percentile scores
//! into a [`HashMap`] at init time. At annotation time, looks up the gene
//! symbol from each transcript consequence.
//!
//! **Usage:** `--plugin LoFtool,/path/to/LoFtool_scores.txt`
//!
//! The scores file is a two-column TSV: `Gene\tLoFtool_percentile`; the path is
//! required.
//!
//! Output field: `LoFtool_percentile`
//!
//! Per-transcript plugin: returns annotation from [`run`] via the default
//! [`run_batch`] dispatcher.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use indexmap::IndexMap;
use tracing::debug;

use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// Output field name.
const OUTPUT_FIELD: &str = "LoFtool_percentile";

/// LoFtool gene intolerance score plugin.
pub struct LoFtoolPlugin {
    /// Gene symbol -> LoFtool percentile score.
    scores: HashMap<String, f64>,
    /// Path to the scores file (for diagnostics).
    file_path: Option<PathBuf>,
}

impl LoFtoolPlugin {
    /// Create a new uninitialized LoFtool plugin.
    pub fn new() -> Self {
        Self {
            scores: HashMap::new(),
            file_path: None,
        }
    }

    /// Load scores from a two-column TSV file (Gene\tScore).
    fn load_scores(path: &str) -> Result<HashMap<String, f64>, PluginError> {
        let file = std::fs::File::open(path)
            .map_err(|e| PluginError::Init(format!("LoFtool: failed to open {path}: {e}")))?;
        let reader = BufReader::new(file);
        let mut scores = HashMap::new();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result.map_err(|e| {
                PluginError::Init(format!(
                    "LoFtool: failed to read line {}: {e}",
                    line_num + 1
                ))
            })?;

            if line.starts_with('#') || line.starts_with("Gene") || line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 2 {
                continue;
            }

            let gene = fields[0].trim();
            if gene.is_empty() {
                continue;
            }

            if let Ok(score) = fields[1].trim().parse::<f64>() {
                scores.insert(gene.to_string(), score);
            }
        }

        Ok(scores)
    }
}

impl Default for LoFtoolPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl BuiltinPlugin for LoFtoolPlugin {
    fn name(&self) -> &str {
        "LoFtool"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            OUTPUT_FIELD.to_string(),
            "LoFtool gene loss-of-function intolerance percentile".to_string(),
        )]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let file_path = params.first().ok_or_else(|| {
            PluginError::Init("LoFtool requires a path to the scores file".into())
        })?;

        self.scores = Self::load_scores(file_path)?;
        self.file_path = Some(PathBuf::from(file_path));

        debug!(
            genes = self.scores.len(),
            path = file_path.as_str(),
            "LoFtool plugin loaded scores"
        );

        if self.scores.is_empty() {
            return Err(PluginError::Init(
                "LoFtool scores file is empty or has no valid entries".into(),
            ));
        }

        Ok(())
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut result = IndexMap::new();

        if let Some(gene_symbol) = &consequence.gene_symbol {
            if let Some(&score) = self.scores.get(&**gene_symbol) {
                result.insert(OUTPUT_FIELD.to_string(), format!("{score}"));
            }
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_variant() -> InputVariant {
        InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec())
    }

    #[test]
    fn plugin_name_is_correct() {
        let plugin = LoFtoolPlugin::new();
        assert_eq!(plugin.name(), "LoFtool");
    }

    #[test]
    fn header_info_has_expected_field() {
        let plugin = LoFtoolPlugin::new();
        let headers = plugin.header_info();
        assert_eq!(headers.len(), 1);
        assert_eq!(headers[0].0, "LoFtool_percentile");
    }

    #[test]
    fn init_fails_without_params() {
        let mut plugin = LoFtoolPlugin::new();
        let result = plugin.init(&[]);
        assert!(result.is_err());
    }

    #[test]
    fn run_returns_empty_for_no_gene_symbol() {
        let plugin = LoFtoolPlugin {
            scores: {
                let mut m = HashMap::new();
                m.insert("SOD1".to_string(), 0.85);
                m
            },
            file_path: None,
        };

        let consequence = TranscriptConsequence::default();
        let variant = test_variant();
        let result = plugin.run(&consequence, &variant).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn run_returns_score_for_matching_gene() {
        let plugin = LoFtoolPlugin {
            scores: {
                let mut m = HashMap::new();
                m.insert("SOD1".to_string(), 0.85);
                m
            },
            file_path: None,
        };

        let consequence = TranscriptConsequence {
            gene_symbol: Some("SOD1".into()),
            ..Default::default()
        };
        let variant = test_variant();

        let result = plugin.run(&consequence, &variant).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result.get("LoFtool_percentile").unwrap(), "0.85");
    }

    #[test]
    fn run_returns_empty_for_unknown_gene() {
        let plugin = LoFtoolPlugin {
            scores: {
                let mut m = HashMap::new();
                m.insert("SOD1".to_string(), 0.85);
                m
            },
            file_path: None,
        };

        let consequence = TranscriptConsequence {
            gene_symbol: Some("UNKNOWN_GENE".into()),
            ..Default::default()
        };
        let variant = test_variant();

        let result = plugin.run(&consequence, &variant).unwrap();
        assert!(result.is_empty());
    }
}
