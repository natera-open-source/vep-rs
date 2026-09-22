// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! gnomADc plugin: gnomAD coverage statistics.
//!
//! Simplest tabix plugin - position-only lookup, no allele matching.
//! Results are stored at the variant level (plugin_data on InputVariant).

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{open_best_store, AnnotationStore};
use crate::tabix::{BatchQueryResult, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// Output field prefix for gnomADc annotations.
const FIELD_PREFIX: &str = "gnomADc";

/// gnomAD coverage statistics plugin.
///
/// Annotates variants with coverage statistics from gnomAD coverage files.
/// This is the simplest tabix plugin: it queries by position only (no allele
/// matching) and stores results at the variant level.
pub struct GnomadcPlugin {
    /// Annotation store, initialized after `init()` is called.
    store: Option<Box<dyn AnnotationStore>>,
    /// Column names from the tabix file header to output.
    output_fields: Vec<String>,
}

impl Default for GnomadcPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl GnomadcPlugin {
    /// Create an uninitialized gnomADc plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            output_fields: Vec::new(),
        }
    }
}

impl BuiltinPlugin for GnomadcPlugin {
    fn name(&self) -> &str {
        "gnomADc"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        self.output_fields
            .iter()
            .map(|field| {
                let key = format!("{FIELD_PREFIX}_{field}");
                let desc = format!("gnomAD coverage statistic: {field}");
                (key, desc)
            })
            .collect()
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let file_path = params
            .first()
            .ok_or_else(|| PluginError::Init("gnomADc requires a file path parameter".into()))?;

        let config = TabixConfig {
            file_path: file_path.into(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        };

        let store = open_best_store(std::path::Path::new(file_path.as_str()), config)?;

        if let Some(header) = store.header() {
            self.output_fields = header
                .iter()
                .skip(2) // skip chrom + pos columns
                .cloned()
                .collect();
        }

        debug!(
            file = file_path,
            fields = ?self.output_fields,
            "gnomADc plugin initialized"
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

        let batch_result = store.query_batch(&regions)?;
        Ok(Some(PrefetchData::new(batch_result)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch_result: BatchQueryResult = data.downcast()?;

        if batch_result.is_empty() {
            return Ok(());
        }

        for variant in variants.iter_mut() {
            let records = batch_result.lookup(&variant.chr, variant.start, variant.end);

            if records.is_empty() {
                continue;
            }

            // Coverage is per-position, so the first record is the only one.
            let record = records[0];

            for field in &self.output_fields {
                if let Some(value) = record.get(field) {
                    if !value.is_empty() && value != "." {
                        let key = format!("{FIELD_PREFIX}_{field}");
                        variant.plugin_data.insert(key, value.to_string());
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
        _consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // Variant-level plugin: annotation happens in `run_batch`.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_plugin_is_uninitialized() {
        let plugin = GnomadcPlugin::new();
        assert_eq!(plugin.name(), "gnomADc");
        assert!(plugin.store.is_none());
        assert!(plugin.output_fields.is_empty());
    }

    #[test]
    fn init_requires_file_path() {
        let mut plugin = GnomadcPlugin::new();
        let result = plugin.init(&[]);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("requires a file path"));
    }

    #[test]
    fn header_info_reflects_output_fields() {
        let mut plugin = GnomadcPlugin::new();
        plugin.output_fields = vec!["mean".into(), "median".into()];
        let headers = plugin.header_info();
        assert_eq!(headers.len(), 2);
        assert_eq!(headers[0].0, "gnomADc_mean");
        assert_eq!(headers[1].0, "gnomADc_median");
    }

    #[test]
    fn run_batch_skips_when_uninitialized() {
        let plugin = GnomadcPlugin::new();
        let mut variants = vec![];
        assert!(plugin.run_batch(&mut variants).is_ok());
    }
}
