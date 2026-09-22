// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! pLI plugin: flat-file gene-level probability of loss-of-function intolerance.
//!
//! Loads a TSV file of gene constraint scores (pLI values from gnomAD/ExAC) into
//! a [`HashMap`] at init time. At annotation time, looks up the gene symbol from
//! each transcript consequence to retrieve the pLI score.
//!
//! **Usage:** `--plugin pLI,/path/to/pLI_values.txt`
//!
//! The scores file is a tab-separated file with at least a gene column and a pLI
//! column. The plugin auto-detects column positions from the header.
//!
//! Output field: `pLI_gene_value`
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
const OUTPUT_FIELD: &str = "pLI_gene_value";

/// Common column header names for the gene symbol in pLI files.
const GENE_COLUMN_NAMES: &[&str] = &["gene", "Gene", "gene_symbol", "GENE", "#gene"];

/// Common column header names for the pLI score.
const PLI_COLUMN_NAMES: &[&str] = &["pLI", "pli", "pLI_score"];

/// pLI gene constraint score plugin.
pub struct PliPlugin {
    /// Gene symbol -> pLI score.
    scores: HashMap<String, f64>,
    /// Path to the scores file (for diagnostics).
    file_path: Option<PathBuf>,
}

impl PliPlugin {
    /// Create a new uninitialized pLI plugin.
    pub fn new() -> Self {
        Self {
            scores: HashMap::new(),
            file_path: None,
        }
    }

    /// Find the column index for a given set of candidate header names.
    fn find_column(headers: &[&str], candidates: &[&str]) -> Option<usize> {
        for candidate in candidates {
            if let Some(idx) = headers.iter().position(|h| h == candidate) {
                return Some(idx);
            }
        }
        None
    }

    /// Load pLI scores from a TSV file, auto-detecting gene and pLI columns
    /// from the header line.
    fn load_scores(path: &str) -> Result<HashMap<String, f64>, PluginError> {
        let file = std::fs::File::open(path)
            .map_err(|e| PluginError::Init(format!("pLI: failed to open {path}: {e}")))?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let header_line = lines
            .next()
            .ok_or_else(|| PluginError::Init("pLI: scores file is empty".into()))?
            .map_err(|e| PluginError::Init(format!("pLI: failed to read header: {e}")))?;

        let header_line = header_line.trim_start_matches('#');
        let headers: Vec<&str> = header_line.split('\t').collect();

        let gene_col = Self::find_column(&headers, GENE_COLUMN_NAMES).ok_or_else(|| {
            PluginError::Init(format!(
                "pLI: could not find gene column in header. Expected one of: {:?}. \
                 Got: {:?}",
                GENE_COLUMN_NAMES, headers
            ))
        })?;

        let pli_col = Self::find_column(&headers, PLI_COLUMN_NAMES).ok_or_else(|| {
            PluginError::Init(format!(
                "pLI: could not find pLI column in header. Expected one of: {:?}. \
                 Got: {:?}",
                PLI_COLUMN_NAMES, headers
            ))
        })?;

        let mut scores = HashMap::new();

        for (line_num, line_result) in lines.enumerate() {
            let line = line_result.map_err(|e| {
                PluginError::Init(format!("pLI: failed to read line {}: {e}", line_num + 2))
            })?;

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            let gene = match fields.get(gene_col) {
                Some(g) => g.trim(),
                None => continue,
            };
            if gene.is_empty() {
                continue;
            }

            let pli_str = match fields.get(pli_col) {
                Some(s) => s.trim(),
                None => continue,
            };

            if let Ok(score) = pli_str.parse::<f64>() {
                scores.insert(gene.to_string(), score);
            }
        }

        Ok(scores)
    }
}

impl Default for PliPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl BuiltinPlugin for PliPlugin {
    fn name(&self) -> &str {
        "pLI"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            OUTPUT_FIELD.to_string(),
            "pLI (probability of loss-of-function intolerance) gene constraint score".to_string(),
        )]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let file_path = params
            .first()
            .ok_or_else(|| PluginError::Init("pLI requires a path to the scores file".into()))?;

        self.scores = Self::load_scores(file_path)?;
        self.file_path = Some(PathBuf::from(file_path));

        debug!(
            genes = self.scores.len(),
            path = file_path.as_str(),
            "pLI plugin loaded scores"
        );

        if self.scores.is_empty() {
            return Err(PluginError::Init(
                "pLI scores file is empty or has no valid entries".into(),
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
        let plugin = PliPlugin::new();
        assert_eq!(plugin.name(), "pLI");
    }

    #[test]
    fn header_info_has_expected_field() {
        let plugin = PliPlugin::new();
        let headers = plugin.header_info();
        assert_eq!(headers.len(), 1);
        assert_eq!(headers[0].0, "pLI_gene_value");
    }

    #[test]
    fn init_fails_without_params() {
        let mut plugin = PliPlugin::new();
        let result = plugin.init(&[]);
        assert!(result.is_err());
    }

    #[test]
    fn find_column_matches_candidates() {
        let headers = vec!["chr", "gene", "pLI", "other"];
        assert_eq!(PliPlugin::find_column(&headers, GENE_COLUMN_NAMES), Some(1));
        assert_eq!(PliPlugin::find_column(&headers, PLI_COLUMN_NAMES), Some(2));
    }

    #[test]
    fn find_column_returns_none_for_missing() {
        let headers = vec!["foo", "bar"];
        assert_eq!(PliPlugin::find_column(&headers, GENE_COLUMN_NAMES), None);
    }

    #[test]
    fn run_returns_empty_for_no_gene_symbol() {
        let plugin = PliPlugin {
            scores: {
                let mut m = HashMap::new();
                m.insert("BRCA1".to_string(), 1.0);
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
        let plugin = PliPlugin {
            scores: {
                let mut m = HashMap::new();
                m.insert("BRCA1".to_string(), 1.0);
                m
            },
            file_path: None,
        };

        let consequence = TranscriptConsequence {
            gene_symbol: Some("BRCA1".into()),
            ..Default::default()
        };
        let variant = test_variant();

        let result = plugin.run(&consequence, &variant).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result.get("pLI_gene_value").unwrap(), "1");
    }

    #[test]
    fn run_returns_empty_for_unknown_gene() {
        let plugin = PliPlugin {
            scores: {
                let mut m = HashMap::new();
                m.insert("BRCA1".to_string(), 1.0);
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
