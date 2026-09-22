// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! SpliceAI plugin: splice site prediction scores.
//!
//! Queries tabix-indexed SpliceAI VCF files (separate SNV and indel files)
//! to annotate transcript consequences with delta scores and positions.
//! Matches by position + allele + gene symbol.
//!
//! # Parameters
//!
//! - `snv=<path>`: path to the tabix-indexed SNV SpliceAI VCF
//! - `indel=<path>`: path to the tabix-indexed indel SpliceAI VCF
//! - `cutoff=<float>`: minimum max delta score to report (optional)
//!
//! # Output fields
//!
//! - `SpliceAI_pred_DS_AG`: delta score for acceptor gain
//! - `SpliceAI_pred_DS_AL`: delta score for acceptor loss
//! - `SpliceAI_pred_DS_DG`: delta score for donor gain
//! - `SpliceAI_pred_DS_DL`: delta score for donor loss
//! - `SpliceAI_pred_DP_AG`: delta position for acceptor gain
//! - `SpliceAI_pred_DP_AL`: delta position for acceptor loss
//! - `SpliceAI_pred_DP_DG`: delta position for donor gain
//! - `SpliceAI_pred_DP_DL`: delta position for donor loss
//! - `SpliceAI_pred_SYMBOL`: gene symbol

use indexmap::IndexMap;
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::{InputVariant, VariantClass};

use crate::annotation_store::{open_best_store, AnnotationStore};
use crate::tabix::{allele_match, BatchQueryResult, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};

/// Parsed SpliceAI scores from a single INFO field entry.
///
/// The SpliceAI VCF INFO field contains pipe-delimited values:
/// `ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`
///
/// Delta score (DS) fields are stored as raw strings to preserve the
/// original formatting in the output, while still supporting numeric
/// comparison for cutoff filtering.
struct SpliceAiScores {
    /// The alternate allele from the SpliceAI annotation.
    ///
    /// Parsed for completeness but not used in output (allele matching is
    /// done at the VCF record level via `allele_match::matches_allele`).
    #[allow(dead_code)]
    allele: String,
    /// Gene symbol.
    symbol: String,
    /// Delta score strings (AG, AL, DG, DL).
    ds_ag: String,
    ds_al: String,
    ds_dg: String,
    ds_dl: String,
    /// Delta position strings (AG, AL, DG, DL).
    dp_ag: String,
    dp_al: String,
    dp_dg: String,
    dp_dl: String,
}

impl SpliceAiScores {
    /// Parse from a pipe-delimited SpliceAI value.
    ///
    /// Returns `None` if the value has fewer than 10 fields or if the
    /// delta score fields are not valid floating-point numbers.
    fn parse(value: &str) -> Option<Self> {
        let parts: Vec<&str> = value.split('|').collect();
        if parts.len() < 10 {
            return None;
        }

        for part in &parts[2..6] {
            part.parse::<f64>().ok()?;
        }

        Some(Self {
            allele: parts[0].to_string(),
            symbol: parts[1].to_string(),
            ds_ag: parts[2].to_string(),
            ds_al: parts[3].to_string(),
            ds_dg: parts[4].to_string(),
            ds_dl: parts[5].to_string(),
            dp_ag: parts[6].to_string(),
            dp_al: parts[7].to_string(),
            dp_dg: parts[8].to_string(),
            dp_dl: parts[9].to_string(),
        })
    }

    /// Maximum delta score across all four splice categories.
    ///
    /// Returns 0.0 if no score can be parsed (should not happen after
    /// validation in [`parse`]).
    fn max_ds(&self) -> f64 {
        [&self.ds_ag, &self.ds_al, &self.ds_dg, &self.ds_dl]
            .iter()
            .filter_map(|s| s.parse::<f64>().ok())
            .fold(0.0_f64, f64::max)
    }

    /// Convert to output key-value pairs for the VEP Extra column.
    fn to_output(&self) -> IndexMap<String, String> {
        let mut out = IndexMap::with_capacity(9);
        out.insert("SpliceAI_pred_DS_AG".to_string(), self.ds_ag.clone());
        out.insert("SpliceAI_pred_DS_AL".to_string(), self.ds_al.clone());
        out.insert("SpliceAI_pred_DS_DG".to_string(), self.ds_dg.clone());
        out.insert("SpliceAI_pred_DS_DL".to_string(), self.ds_dl.clone());
        out.insert("SpliceAI_pred_DP_AG".to_string(), self.dp_ag.clone());
        out.insert("SpliceAI_pred_DP_AL".to_string(), self.dp_al.clone());
        out.insert("SpliceAI_pred_DP_DG".to_string(), self.dp_dg.clone());
        out.insert("SpliceAI_pred_DP_DL".to_string(), self.dp_dl.clone());
        out.insert("SpliceAI_pred_SYMBOL".to_string(), self.symbol.clone());
        out
    }
}

/// VCF REF column index (CHROM=0, POS=1, ID=2, REF=3).
const VCF_REF_COL: usize = 3;
/// VCF ALT column index (CHROM=0, POS=1, ID=2, REF=3, ALT=4).
const VCF_ALT_COL: usize = 4;

/// SpliceAI splice prediction plugin.
///
/// SpliceAI predicts splicing effects for SNVs and indels using a deep
/// neural network. Predictions are distributed as two tabix-indexed VCF
/// files (one for SNVs, one for indels). This plugin queries the
/// appropriate file based on the variant class, matches by allele and
/// gene symbol, and optionally applies a delta score cutoff.
pub struct SpliceAiPlugin {
    /// Annotation store for the SNV SpliceAI VCF file.
    snv_store: Option<Box<dyn AnnotationStore>>,
    /// Annotation store for the indel SpliceAI VCF file.
    indel_store: Option<Box<dyn AnnotationStore>>,
    /// Minimum max delta score to include results.
    cutoff: Option<f64>,
}

impl Default for SpliceAiPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl SpliceAiPlugin {
    /// Create an uninitialized SpliceAI plugin.
    ///
    /// Call [`init`](BuiltinPlugin::init) with SNV/indel file paths before use.
    pub fn new() -> Self {
        Self {
            snv_store: None,
            indel_store: None,
            cutoff: None,
        }
    }

    /// Build a [`TabixConfig`] for a VCF-format SpliceAI file.
    ///
    /// VCF columns: CHROM(0), POS(1), ID(2), REF(3), ALT(4), QUAL(5),
    /// FILTER(6), INFO(7).
    fn vcf_config(path: &str) -> TabixConfig {
        TabixConfig {
            file_path: path.into(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: Some(3),
            alt_col: Some(4),
            zero_based: false,
        }
    }

    /// Select the correct annotation store based on variant class.
    ///
    /// SNVs use the SNV file; insertions, deletions, indels, and other
    /// variant types use the indel file.
    fn store_for(&self, variant: &InputVariant) -> Option<&dyn AnnotationStore> {
        match variant.variant_class {
            VariantClass::Snv => self.snv_store.as_deref(),
            _ => self.indel_store.as_deref(),
        }
    }

    /// Parse all SpliceAI score entries from a VCF INFO field.
    ///
    /// The INFO field may contain multiple annotations separated by commas
    /// within the `SpliceAI=` key (one per overlapping gene/allele).
    fn parse_info(info: &str) -> Vec<SpliceAiScores> {
        let mut results = Vec::new();
        for field in info.split(';') {
            if let Some(value) = field.strip_prefix("SpliceAI=") {
                for entry in value.split(',') {
                    if let Some(scores) = SpliceAiScores::parse(entry) {
                        results.push(scores);
                    }
                }
            }
        }
        results
    }

    fn query_bounds(variant: &InputVariant) -> (u64, u64) {
        let query_start = variant.start.saturating_sub(1).max(1);
        let query_end = variant.end.max(variant.start);
        (query_start, query_end)
    }

    fn record_info(record: &crate::tabix::TabixRecord) -> Option<&str> {
        record
            .get("INFO")
            .or_else(|| record.columns.get(7).map(|s| s.as_str()))
    }
}

impl BuiltinPlugin for SpliceAiPlugin {
    fn name(&self) -> &str {
        "SpliceAI"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            (
                "SpliceAI_pred_DS_AG".into(),
                "SpliceAI delta score for acceptor gain".into(),
            ),
            (
                "SpliceAI_pred_DS_AL".into(),
                "SpliceAI delta score for acceptor loss".into(),
            ),
            (
                "SpliceAI_pred_DS_DG".into(),
                "SpliceAI delta score for donor gain".into(),
            ),
            (
                "SpliceAI_pred_DS_DL".into(),
                "SpliceAI delta score for donor loss".into(),
            ),
            (
                "SpliceAI_pred_DP_AG".into(),
                "SpliceAI delta position for acceptor gain".into(),
            ),
            (
                "SpliceAI_pred_DP_AL".into(),
                "SpliceAI delta position for acceptor loss".into(),
            ),
            (
                "SpliceAI_pred_DP_DG".into(),
                "SpliceAI delta position for donor gain".into(),
            ),
            (
                "SpliceAI_pred_DP_DL".into(),
                "SpliceAI delta position for donor loss".into(),
            ),
            ("SpliceAI_pred_SYMBOL".into(), "SpliceAI gene symbol".into()),
        ]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let mut snv_path = None;
        let mut indel_path = None;

        for param in params {
            if let Some(path) = param.strip_prefix("snv=") {
                snv_path = Some(path.to_string());
            } else if let Some(path) = param.strip_prefix("indel=") {
                indel_path = Some(path.to_string());
            } else if let Some(val) = param.strip_prefix("cutoff=") {
                self.cutoff = Some(val.parse::<f64>().map_err(|e| {
                    PluginError::Init(format!("invalid cutoff value '{val}': {e}"))
                })?);
            }
        }

        if snv_path.is_none() && indel_path.is_none() {
            return Err(PluginError::Init(
                "SpliceAI plugin requires at least one of 'snv' or 'indel' parameters".into(),
            ));
        }

        if let Some(path) = snv_path {
            self.snv_store = Some(open_best_store(
                std::path::Path::new(&path),
                Self::vcf_config(&path),
            )?);
        }
        if let Some(path) = indel_path {
            self.indel_store = Some(open_best_store(
                std::path::Path::new(&path),
                Self::vcf_config(&path),
            )?);
        }

        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        if variants.is_empty() {
            return Ok(None);
        }
        if self.snv_store.is_none() && self.indel_store.is_none() {
            return Ok(None);
        }

        let mut snv_regions = Vec::new();
        let mut indel_regions = Vec::new();

        for variant in variants.iter() {
            let (query_start, query_end) = Self::query_bounds(variant);
            match variant.variant_class {
                VariantClass::Snv => {
                    snv_regions.push((variant.chr.as_str(), query_start, query_end))
                }
                _ => indel_regions.push((variant.chr.as_str(), query_start, query_end)),
            }
        }

        let snv_batch = match &self.snv_store {
            Some(s) => Some(s.query_batch(&snv_regions)?),
            None => None,
        };
        let indel_batch = match &self.indel_store {
            Some(s) => Some(s.query_batch(&indel_regions)?),
            None => None,
        };

        Ok(Some(PrefetchData::new((snv_batch, indel_batch))))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let (snv_batch, indel_batch): (Option<BatchQueryResult>, Option<BatchQueryResult>) =
            data.downcast()?;

        for variant in variants.iter_mut() {
            let batch = match variant.variant_class {
                VariantClass::Snv => match (&self.snv_store, &snv_batch) {
                    (Some(_), Some(b)) => b,
                    _ => continue,
                },
                _ => match (&self.indel_store, &indel_batch) {
                    (Some(_), Some(b)) => b,
                    _ => continue,
                },
            };

            let (query_start, query_end) = Self::query_bounds(variant);
            let records = batch.lookup(&variant.chr, query_start, query_end);
            if records.is_empty() {
                continue;
            }

            for tc_idx in 0..variant.transcript_consequences.len() {
                let gene_symbol = variant.transcript_consequences[tc_idx]
                    .gene_symbol
                    .as_deref();
                let mut output = None;

                'records: for record in &records {
                    if !allele_match::matches_allele(record, variant, VCF_REF_COL, VCF_ALT_COL) {
                        continue;
                    }

                    let Some(info) = Self::record_info(record) else {
                        continue;
                    };
                    let entries = Self::parse_info(info);
                    for entry in &entries {
                        if let Some(symbol) = gene_symbol {
                            if !entry.symbol.eq_ignore_ascii_case(symbol) {
                                continue;
                            }
                        }
                        if let Some(cutoff) = self.cutoff {
                            if entry.max_ds() < cutoff {
                                continue;
                            }
                        }
                        output = Some(entry.to_output());
                        break 'records;
                    }
                }

                if let Some(values) = output {
                    for (k, v) in values {
                        variant.transcript_consequences[tc_idx]
                            .plugin_data
                            .insert(k, v);
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
        let store = match self.store_for(variant) {
            Some(s) => s,
            None => return Ok(IndexMap::new()),
        };

        // A VCF indel position is one less than VEP's (VCF includes the anchor
        // base), so the range covers both SNV and indel positions.
        let (query_start, query_end) = Self::query_bounds(variant);

        let records = store.query(&variant.chr, query_start, query_end)?;
        let gene_symbol = consequence.gene_symbol.as_deref();

        for record in &records {
            if !allele_match::matches_allele(record, variant, VCF_REF_COL, VCF_ALT_COL) {
                continue;
            }

            let info = match Self::record_info(record) {
                Some(info) => info,
                None => continue,
            };

            let entries = Self::parse_info(info);

            for entry in &entries {
                if let Some(symbol) = gene_symbol {
                    if !entry.symbol.eq_ignore_ascii_case(symbol) {
                        continue;
                    }
                }

                if let Some(cutoff) = self.cutoff {
                    if entry.max_ds() < cutoff {
                        continue;
                    }
                }

                return Ok(entry.to_output());
            }
        }

        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_plugin_is_uninitialized() {
        let plugin = SpliceAiPlugin::new();
        assert_eq!(plugin.name(), "SpliceAI");
        assert!(plugin.snv_store.is_none());
        assert!(plugin.indel_store.is_none());
        assert!(plugin.cutoff.is_none());
    }

    #[test]
    fn header_info_returns_all_fields() {
        let plugin = SpliceAiPlugin::new();
        let headers = plugin.header_info();
        assert_eq!(headers.len(), 9);

        let field_names: Vec<&str> = headers.iter().map(|(k, _)| k.as_str()).collect();
        assert!(field_names.contains(&"SpliceAI_pred_DS_AG"));
        assert!(field_names.contains(&"SpliceAI_pred_DS_AL"));
        assert!(field_names.contains(&"SpliceAI_pred_DS_DG"));
        assert!(field_names.contains(&"SpliceAI_pred_DS_DL"));
        assert!(field_names.contains(&"SpliceAI_pred_DP_AG"));
        assert!(field_names.contains(&"SpliceAI_pred_DP_AL"));
        assert!(field_names.contains(&"SpliceAI_pred_DP_DG"));
        assert!(field_names.contains(&"SpliceAI_pred_DP_DL"));
        assert!(field_names.contains(&"SpliceAI_pred_SYMBOL"));
    }

    #[test]
    fn run_returns_empty_when_not_initialized() {
        let plugin = SpliceAiPlugin::new();
        let consequence = TranscriptConsequence::default();
        let variant = InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let result = plugin.run(&consequence, &variant).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn parse_info_basic() {
        let info = "SpliceAI=T|SOD1|0.00|0.01|0.00|0.00|2|37|-22|37";
        let entries = SpliceAiPlugin::parse_info(info);
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].allele, "T");
        assert_eq!(entries[0].symbol, "SOD1");
        assert_eq!(entries[0].ds_ag, "0.00");
        assert_eq!(entries[0].ds_al, "0.01");
    }

    #[test]
    fn parse_info_multiple_entries() {
        let info =
            "SpliceAI=T|GENE1|0.10|0.20|0.30|0.40|1|2|3|4,T|GENE2|0.50|0.60|0.70|0.80|5|6|7|8";
        let entries = SpliceAiPlugin::parse_info(info);
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].symbol, "GENE1");
        assert_eq!(entries[1].symbol, "GENE2");
    }

    #[test]
    fn parse_info_with_other_fields() {
        let info = "AF=0.01;SpliceAI=G|TP53|0.05|0.00|0.00|0.10|10|-5|3|20;DB";
        let entries = SpliceAiPlugin::parse_info(info);
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].allele, "G");
        assert_eq!(entries[0].symbol, "TP53");
        assert_eq!(entries[0].ds_dl, "0.10");
    }

    #[test]
    fn parse_info_invalid() {
        assert!(SpliceAiPlugin::parse_info("").is_empty());
        assert!(SpliceAiPlugin::parse_info("AF=0.01").is_empty());
        assert!(SpliceAiPlugin::parse_info("SpliceAI=T|GENE").is_empty());
    }

    #[test]
    fn max_ds_calculation() {
        let scores = SpliceAiScores::parse("T|SOD1|0.10|0.50|0.20|0.30|1|2|3|4").unwrap();
        assert!((scores.max_ds() - 0.50).abs() < f64::EPSILON);
    }

    #[test]
    fn max_ds_all_zero() {
        let scores = SpliceAiScores::parse("T|SOD1|0.00|0.00|0.00|0.00|1|2|3|4").unwrap();
        assert!((scores.max_ds() - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn to_output_has_all_fields() {
        let scores = SpliceAiScores::parse("T|SOD1|0.10|0.20|0.30|0.40|1|2|3|4").unwrap();
        let output = scores.to_output();
        assert_eq!(output.len(), 9);
        assert_eq!(output["SpliceAI_pred_DS_AG"], "0.10");
        assert_eq!(output["SpliceAI_pred_DS_AL"], "0.20");
        assert_eq!(output["SpliceAI_pred_DS_DG"], "0.30");
        assert_eq!(output["SpliceAI_pred_DS_DL"], "0.40");
        assert_eq!(output["SpliceAI_pred_DP_AG"], "1");
        assert_eq!(output["SpliceAI_pred_DP_AL"], "2");
        assert_eq!(output["SpliceAI_pred_DP_DG"], "3");
        assert_eq!(output["SpliceAI_pred_DP_DL"], "4");
        assert_eq!(output["SpliceAI_pred_SYMBOL"], "SOD1");
    }

    #[test]
    fn store_selection_by_variant_class() {
        let plugin = SpliceAiPlugin::new();

        let snv = InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert!(plugin.store_for(&snv).is_none());

        let del = InputVariant::new("1".into(), 100, 102, b"ACG".to_vec(), b"-".to_vec());
        assert!(plugin.store_for(&del).is_none());
    }
}
