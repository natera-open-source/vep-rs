// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! dbNSFP plugin: tabix-based functional prediction scores for nonsynonymous SNVs.
//!
//! dbNSFP aggregates 100+ prediction/conservation scores from many tools. This
//! plugin lets the user select which columns to extract via parameters, making it
//! the most configurable tabix plugin in the suite.
//!
//! **Usage:** `--plugin dbNSFP,/path/to/dbNSFP.gz,field1,field2,...`
//!
//! Two-phase plugin: [`BuiltinPlugin::prefetch`] batches the tabix query for a
//! buffer and [`BuiltinPlugin::annotate`] writes per-transcript results, matching
//! each consequence on chromosome, position, ref, alt, and optionally transcript
//! (Ensembl_transcriptid).

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{open_best_store, AnnotationStore};
use crate::tabix::{BatchQueryResult, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// Standard column layout in dbNSFP files.
const CHR_COL: usize = 0;
const POS_COL: usize = 1;
const REF_COL: usize = 2;
const ALT_COL: usize = 3;

/// dbNSFP multi-column annotation plugin.
pub struct DbNsfpPlugin {
    /// Annotation store, initialized in `init()`.
    store: Option<Box<dyn AnnotationStore>>,
    /// User-selected field names to extract from dbNSFP columns.
    selected_fields: Vec<String>,
    /// Column header names resolved at init (from the tabix file header).
    column_headers: Vec<String>,
}

impl DbNsfpPlugin {
    /// Create a new uninitialized dbNSFP plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            selected_fields: Vec::new(),
            column_headers: Vec::new(),
        }
    }

    /// Look up the 0-based column index for a field name.
    fn field_index(&self, name: &str) -> Option<usize> {
        self.column_headers.iter().position(|h| h == name)
    }

    /// Check whether a dbNSFP record matches the variant's alleles.
    fn allele_matches(&self, columns: &[String], variant: &InputVariant) -> bool {
        let file_ref = match columns.get(REF_COL) {
            Some(r) => r.as_str(),
            None => return false,
        };
        let file_alt = match columns.get(ALT_COL) {
            Some(a) => a.as_str(),
            None => return false,
        };

        let var_ref = String::from_utf8_lossy(&variant.ref_allele);
        let var_alt = String::from_utf8_lossy(variant.alt_allele());

        file_ref.eq_ignore_ascii_case(&var_ref) && file_alt.eq_ignore_ascii_case(&var_alt)
    }

    /// For a matched record, find the index within the semicolon-delimited
    /// `Ensembl_transcriptid` column that matches the consequence's transcript,
    /// so multi-valued fields yield the matching sub-value.
    fn transcript_index(&self, columns: &[String], transcript_id: &str) -> Option<usize> {
        let enst_col = self.field_index("Ensembl_transcriptid")?;
        let enst_value = columns.get(enst_col)?;
        enst_value
            .split(';')
            .position(|tid| tid == transcript_id || tid.starts_with(transcript_id))
    }

    /// Extract a (possibly multi-valued) column value. If `tx_index` is Some,
    /// pick the value at that index from semicolon-delimited values; otherwise
    /// return the full value.
    fn extract_field_value(
        &self,
        columns: &[String],
        field: &str,
        tx_index: Option<usize>,
    ) -> Option<String> {
        let col_idx = self.field_index(field)?;
        let raw = columns.get(col_idx)?;
        if raw.is_empty() || raw == "." {
            return None;
        }

        if let Some(idx) = tx_index {
            let parts: Vec<&str> = raw.split(';').collect();
            let val = parts.get(idx).copied().unwrap_or(parts.last().copied()?);
            if val.is_empty() || val == "." {
                None
            } else {
                Some(val.to_string())
            }
        } else {
            Some(raw.clone())
        }
    }
}

impl Default for DbNsfpPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl BuiltinPlugin for DbNsfpPlugin {
    fn name(&self) -> &str {
        "dbNSFP"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        self.selected_fields
            .iter()
            .map(|f| {
                (
                    format!("dbNSFP_{f}"),
                    format!("dbNSFP annotation field: {f}"),
                )
            })
            .collect()
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        if params.is_empty() {
            return Err(PluginError::Init(
                "dbNSFP requires at least a file path argument".into(),
            ));
        }

        let file_path = &params[0];

        self.selected_fields = params[1..].to_vec();
        if self.selected_fields.is_empty() {
            return Err(PluginError::Init(
                "dbNSFP requires at least one field name (e.g., SIFT_score,Polyphen2_HDIV_score)"
                    .into(),
            ));
        }

        let config = TabixConfig {
            file_path: file_path.into(),
            chr_col: CHR_COL,
            start_col: POS_COL,
            end_col: None,
            ref_col: Some(REF_COL),
            alt_col: Some(ALT_COL),
            zero_based: false,
        };

        let store = open_best_store(std::path::Path::new(file_path.as_str()), config)?;

        self.column_headers = store.header().map(|h| h.to_vec()).unwrap_or_default();

        if !self.column_headers.is_empty() {
            for field in &self.selected_fields {
                if self.field_index(field).is_none() {
                    return Err(PluginError::Init(format!(
                        "dbNSFP field '{field}' not found in file header. \
                         Available columns: {}",
                        self.column_headers.join(", ")
                    )));
                }
            }
        }

        debug!(
            fields = ?self.selected_fields,
            columns = self.column_headers.len(),
            "dbNSFP plugin initialized with {}",
            file_path
        );

        self.store = Some(store);
        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        if variants.is_empty() {
            return Ok(None);
        }

        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .map(|v| (v.chr.as_str(), v.start, v.end))
            .collect();
        let batch = store.query_batch(&regions)?;
        Ok(Some(PrefetchData::new(batch)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch: BatchQueryResult = data.downcast()?;

        for variant in variants.iter_mut() {
            let records = batch.lookup(&variant.chr, variant.start, variant.end);
            if records.is_empty() {
                continue;
            }

            for tc_idx in 0..variant.transcript_consequences.len() {
                let transcript_id = variant.transcript_consequences[tc_idx]
                    .transcript_id
                    .clone();
                let mut result = IndexMap::new();

                for record in &records {
                    if !self.allele_matches(&record.columns, variant) {
                        continue;
                    }

                    let tx_index = self.transcript_index(&record.columns, &transcript_id);
                    for field in &self.selected_fields {
                        if let Some(value) =
                            self.extract_field_value(&record.columns, field, tx_index)
                        {
                            result.insert(format!("dbNSFP_{field}"), value);
                        }
                    }

                    if !result.is_empty() {
                        break;
                    }
                }

                if !result.is_empty() {
                    variant.transcript_consequences[tc_idx]
                        .plugin_data
                        .extend(result);
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
        let store = self
            .store
            .as_ref()
            .ok_or_else(|| PluginError::Run("dbNSFP plugin not initialized".into()))?;

        let records = store.query(&variant.chr, variant.start, variant.end)?;

        let mut result = IndexMap::new();

        for record in &records {
            if !self.allele_matches(&record.columns, variant) {
                continue;
            }

            let tx_index = self.transcript_index(&record.columns, &consequence.transcript_id);

            for field in &self.selected_fields {
                if let Some(value) = self.extract_field_value(&record.columns, field, tx_index) {
                    result.insert(format!("dbNSFP_{field}"), value);
                }
            }

            if !result.is_empty() {
                break;
            }
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plugin_name_is_correct() {
        let plugin = DbNsfpPlugin::new();
        assert_eq!(plugin.name(), "dbNSFP");
    }

    #[test]
    fn header_info_prefixes_fields() {
        let mut plugin = DbNsfpPlugin::new();
        plugin.selected_fields = vec!["SIFT_score".into(), "Polyphen2_HDIV_score".into()];
        let headers = plugin.header_info();
        assert_eq!(headers.len(), 2);
        assert_eq!(headers[0].0, "dbNSFP_SIFT_score");
        assert_eq!(headers[1].0, "dbNSFP_Polyphen2_HDIV_score");
    }

    #[test]
    fn init_fails_without_params() {
        let mut plugin = DbNsfpPlugin::new();
        assert!(plugin.init(&[]).is_err());
    }

    #[test]
    fn init_fails_without_fields() {
        let mut plugin = DbNsfpPlugin::new();
        let result = plugin.init(&["/path/to/file.gz".to_string()]);
        assert!(result.is_err());
    }

    #[test]
    fn field_index_lookup() {
        let mut plugin = DbNsfpPlugin::new();
        plugin.column_headers = vec![
            "chr".into(),
            "pos".into(),
            "ref".into(),
            "alt".into(),
            "SIFT_score".into(),
        ];
        assert_eq!(plugin.field_index("SIFT_score"), Some(4));
        assert_eq!(plugin.field_index("nonexistent"), None);
    }

    #[test]
    fn extract_field_value_single() {
        let mut plugin = DbNsfpPlugin::new();
        plugin.column_headers = vec!["chr".into(), "score".into()];
        let columns = vec!["1".into(), "0.95".into()];
        assert_eq!(
            plugin.extract_field_value(&columns, "score", None),
            Some("0.95".into())
        );
    }

    #[test]
    fn extract_field_value_multi_with_index() {
        let mut plugin = DbNsfpPlugin::new();
        plugin.column_headers = vec!["chr".into(), "score".into()];
        let columns = vec!["1".into(), "0.1;0.5;0.9".into()];
        assert_eq!(
            plugin.extract_field_value(&columns, "score", Some(1)),
            Some("0.5".into())
        );
    }

    #[test]
    fn extract_field_value_dot_is_none() {
        let mut plugin = DbNsfpPlugin::new();
        plugin.column_headers = vec!["chr".into(), "score".into()];
        let columns = vec!["1".into(), ".".into()];
        assert_eq!(plugin.extract_field_value(&columns, "score", None), None);
    }
}
