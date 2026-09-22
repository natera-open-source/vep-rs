// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Variation types from the VEP variation cache (_var.gz files).

use std::collections::HashMap;

/// A variation from the cache (one record from _var.gz files).
///
/// These are the known variants stored in the VEP variation cache,
/// used for co-location lookups during annotation.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct CachedVariation {
    /// Variation name / identifier (e.g., "rs699").
    pub variation_name: String,
    /// Failed status: 0 = not failed.
    pub failed: u8,
    /// Somatic flag: 0 = germline, 1 = somatic.
    pub somatic: u8,
    /// 1-based start position.
    pub start: u64,
    /// 1-based end position (inclusive).
    pub end: u64,
    /// Allele string (e.g., "A/G", "C/A/T").
    pub allele_string: String,
    /// Strand: 1 or -1.
    pub strand: i8,
    /// Minor allele (e.g., "G").
    pub minor_allele: Option<String>,
    /// Minor allele frequency.
    pub minor_allele_freq: Option<f64>,
    /// Clinical significance terms (e.g., "benign,likely_benign").
    pub clin_sig: Option<String>,
    /// Phenotype or disease association flag: 0 or 1.
    pub phenotype_or_disease: u8,
    /// PubMed IDs (comma-separated).
    pub pubmed: Option<String>,
    /// Population allele frequencies.
    /// Key = population name, value = "allele:freq,allele:freq" encoded string.
    pub frequencies: HashMap<String, String>,
    /// Variant synonyms (encoded string from cache).
    pub var_synonyms: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cached_variation_default() {
        let v = CachedVariation::default();
        assert_eq!(v.variation_name, "");
        assert_eq!(v.failed, 0);
        assert_eq!(v.somatic, 0);
        assert_eq!(v.start, 0);
        assert_eq!(v.end, 0);
        assert_eq!(v.allele_string, "");
        assert_eq!(v.strand, 0);
        assert!(v.minor_allele.is_none());
        assert!(v.minor_allele_freq.is_none());
        assert!(v.clin_sig.is_none());
        assert_eq!(v.phenotype_or_disease, 0);
        assert!(v.pubmed.is_none());
        assert!(v.frequencies.is_empty());
        assert!(v.var_synonyms.is_none());
    }

    #[test]
    fn test_cached_variation_construction() {
        let v = CachedVariation {
            variation_name: "rs699".to_string(),
            start: 230710048,
            end: 230710048,
            allele_string: "A/G".to_string(),
            strand: 1,
            minor_allele: Some("G".to_string()),
            minor_allele_freq: Some(0.2345),
            clin_sig: Some("benign".to_string()),
            ..Default::default()
        };
        assert_eq!(v.variation_name, "rs699");
        assert_eq!(v.minor_allele_freq, Some(0.2345));
    }
}
