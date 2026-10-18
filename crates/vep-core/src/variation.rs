// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Variation records of the Ensembl VEP variation cache.

use std::collections::HashMap;

/// One known variant read from the variation cache, the hash Ensembl VEP's
/// `BaseCacheVariation.pm::parse_variation` builds at release/116.2.
///
/// The named fields are the cache's fixed columns; every other column of the
/// cache's `variation_cols` is a population frequency and lands in
/// [`frequencies`](Self::frequencies) under its own name, so the record carries
/// whatever populations the cache declares.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CachedVariation {
    /// Variation name / identifier (e.g., "rs699").
    pub variation_name: String,
    /// Failed status: 0 = not failed.
    pub failed: u8,
    /// Somatic flag: 0 = germline, 1 = somatic.
    pub somatic: u8,
    /// 1-based start position.
    pub start: u64,
    /// 1-based end position (inclusive); an insertion has `end = start - 1`.
    pub end: u64,
    /// Allele string (e.g., "A/G", "C/A/T"); a null allele has no `/`.
    pub allele_string: String,
    /// Strand: 1 or -1.
    pub strand: i8,
    /// Clinical significance terms (e.g., "benign,likely_benign").
    pub clin_sig: Option<String>,
    /// Phenotype or disease association flag: 0 or 1.
    pub phenotype_or_disease: u8,
    /// The ClinVar allele and somatic fields, present on the few records that
    /// carry any of them; boxed so the records without them pay one pointer.
    pub clinvar: Option<Box<ClinVarDetail>>,
    /// PubMed IDs (comma-separated).
    pub pubmed: Option<String>,
    /// Variant synonyms (encoded string from cache).
    pub var_synonyms: Option<String>,
    /// Population allele frequencies.
    /// Key = population column name, value = "allele:freq,allele:freq" encoded string.
    pub frequencies: HashMap<String, String>,
}

impl CachedVariation {
    /// The `(allele, term)` pairs of `clin_sig_allele`, in the cache's order.
    ///
    /// An entry is `allele:term`; the allele is the text before its first `:`
    /// and the term the text between the first and the second (`OutputFactory.pm:1106-1108`
    /// at release/116.2), so an entry with no `:` has an empty term.
    pub fn clin_sig_allele_entries(&self) -> Vec<(&str, &str)> {
        let Some(entries) = self
            .clinvar
            .as_deref()
            .and_then(|clinvar| clinvar.clin_sig_allele.as_deref())
        else {
            return Vec::new();
        };
        entries
            .split(';')
            .map(|entry| {
                let (allele, rest) = entry.split_once(':').unwrap_or((entry, ""));
                let term = rest.split(':').next().unwrap_or("");
                (allele, term)
            })
            .collect()
    }
}

/// The ClinVar fields of a variation record beyond `clin_sig`.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ClinVarDetail {
    /// Significance per risk allele, `allele:term` pairs joined by `;`.
    pub clin_sig_allele: Option<String>,
    /// The reference allele behind `clin_sig_allele`.
    pub clin_sig_ref_allele: Option<String>,
    /// Somatic classification, `phenotype(classification)` entries joined by `;`.
    pub clinical_impact: Option<String>,
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
        assert!(v.clin_sig.is_none());
        assert_eq!(v.phenotype_or_disease, 0);
        assert!(v.clinvar.is_none());
        assert!(v.pubmed.is_none());
        assert!(v.frequencies.is_empty());
        assert!(v.var_synonyms.is_none());
        assert!(v.clin_sig_allele_entries().is_empty());
    }

    /// The `rs142566111` fixture is that record's `clin_sig_allele` in the
    /// release-116 GRCh38 cache.
    #[test]
    fn clin_sig_allele_entries_split_on_semicolon_then_colon() {
        let v = CachedVariation {
            variation_name: "rs142566111".to_string(),
            allele_string: "T/A/G".to_string(),
            clinvar: Some(Box::new(ClinVarDetail {
                clin_sig_allele: Some("G:likely_benign;A:likely_benign".to_string()),
                ..Default::default()
            })),
            ..Default::default()
        };
        assert_eq!(
            v.clin_sig_allele_entries(),
            vec![("G", "likely_benign"), ("A", "likely_benign")]
        );
        let odd = CachedVariation {
            clinvar: Some(Box::new(ClinVarDetail {
                clin_sig_allele: Some(":pathogenic;C:a:b;G".to_string()),
                ..Default::default()
            })),
            ..Default::default()
        };
        assert_eq!(
            odd.clin_sig_allele_entries(),
            vec![("", "pathogenic"), ("C", "a"), ("G", "")]
        );
    }
}
