// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VEP cache metadata from info.txt files.

use std::collections::BTreeMap;

/// Metadata from a VEP cache info.txt file.
///
/// Each VEP cache directory contains an `info.txt` file describing the cache
/// contents, assembly, version, and available annotations. A JSON cache carries
/// the same object as `info.json`; every field is optional there.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct CacheInfo {
    /// Species name (e.g., "homo_sapiens").
    pub species: String,
    /// Genome assembly (e.g., "GRCh38"); `None` when the key is absent, and
    /// `Some("")` for a present but empty value, which Ensembl VEP's reader keeps
    /// and compares (`CacheDir.pm:466, 475, 407` at release/116.2). An absent key
    /// stays absent when the object is written.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assembly: Option<String>,
    /// Cache version number (e.g., 115).
    pub cache_version: Option<u32>,
    /// SIFT predictions availability: "b" = both prediction + score, "p" = prediction only, "s" = score only.
    pub sift: Option<String>,
    /// PolyPhen predictions availability: "b" = both prediction + score.
    pub polyphen: Option<String>,
    /// Whether regulatory features are included.
    pub regulatory: bool,
    /// Available cell types for regulatory features.
    pub cell_types: Vec<String>,
    /// Column names for the variation cache TSV files.
    pub variation_cols: Vec<String>,
    /// Serializer type: "sereal" or None for Storable.
    pub serialiser_type: Option<String>,
    /// Variation file type: "tabix" for tabix-indexed, None for flat files.
    pub var_type: Option<String>,
    /// Build identifier.
    pub build: Option<String>,
    /// Source database versions (e.g., "ClinVar" -> "202306"), in the byte
    /// order the output headers print them.
    pub source_versions: BTreeMap<String, String>,
    /// The cache's chromosome directories in the order the converter and the
    /// builder list them: autosomes by number, X, Y, MT, then every other
    /// contig by name. Empty when the cache's metadata predates the key.
    pub contig_order: Vec<String>,
}

impl CacheInfo {
    /// Whether the variation files are one indexed `all_vars.gz` per chromosome
    /// (`var_type tabix`, `CacheDir.pm:184` at release/116.2) rather than one
    /// `<start>-<end>_var.gz` per region.
    pub fn variation_is_tabix(&self) -> bool {
        self.var_type.as_deref() == Some("tabix")
    }

    /// The column map of the variation files, when the cache declares one.
    pub fn variation_columns(&self) -> Option<VariationColumns> {
        (!self.variation_cols.is_empty())
            .then(|| VariationColumns::new(self.variation_cols.clone()))
    }
}

/// The columns Ensembl VEP's variation dumper writes before the population
/// frequencies (`Database/Variation.pm:70-80` and
/// `Pipeline/DumpVEP/Dumper/Variation.pm:181-188, 195` at release/116.2); every
/// other name in `variation_cols` is a population frequency column.
pub const FIXED_VARIATION_COLUMNS: &[&str] = &[
    "chr",
    "variation_name",
    "failed",
    "somatic",
    "start",
    "end",
    "allele_string",
    "strand",
    "clin_sig",
    "phenotype_or_disease",
    "clin_sig_allele",
    "clin_sig_ref_allele",
    "clinical_impact",
    "pubmed",
    "var_synonyms",
];

/// The column map of a variation cache: the names of its `variation_cols` line
/// in file order, which is the order of the fields on every record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariationColumns {
    names: Vec<String>,
}

impl VariationColumns {
    pub fn new(names: Vec<String>) -> Self {
        Self { names }
    }

    /// The column names in file order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// The 0-based field index of a column.
    pub fn position(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }

    /// The population frequency columns: every name outside
    /// [`FIXED_VARIATION_COLUMNS`], in file order.
    pub fn frequency_columns(&self) -> impl Iterator<Item = &str> {
        self.names
            .iter()
            .map(String::as_str)
            .filter(|name| !FIXED_VARIATION_COLUMNS.contains(name))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The `variation_cols` line of the release-116 GRCh38 cache's `info.txt`:
    /// the fixed columns, then the 1000 Genomes, gnomAD exome and gnomAD genome
    /// populations.
    const RELEASE_116_GRCH38_COLS: &str = "chr,variation_name,failed,somatic,start,end,allele_string,strand,clin_sig,phenotype_or_disease,clin_sig_allele,clin_sig_ref_allele,clinical_impact,pubmed,var_synonyms,AF,AFR,AMR,EAS,EUR,SAS,gnomADe,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_NFE,gnomADe_REMAINING,gnomADe_SAS,gnomADe_MID,gnomADg,gnomADg_AFR,gnomADg_AMI,gnomADg_AMR,gnomADg_ASJ,gnomADg_EAS,gnomADg_FIN,gnomADg_MID,gnomADg_NFE,gnomADg_REMAINING,gnomADg_SAS";

    #[test]
    fn variation_columns_split_fixed_from_frequency_names() {
        let info = CacheInfo {
            variation_cols: RELEASE_116_GRCH38_COLS
                .split(',')
                .map(str::to_string)
                .collect(),
            var_type: Some("tabix".to_string()),
            ..Default::default()
        };
        assert!(info.variation_is_tabix());
        let columns = info.variation_columns().unwrap();
        assert_eq!(columns.names().len(), 42);
        assert_eq!(columns.position("chr"), Some(0));
        assert_eq!(columns.position("start"), Some(4));
        assert_eq!(columns.position("gnomADg_AMI"), Some(33));
        assert_eq!(columns.position("ASN"), None);
        let frequencies: Vec<&str> = columns.frequency_columns().collect();
        assert_eq!(frequencies.len(), 27);
        assert_eq!(frequencies[0], "AF");
        assert_eq!(frequencies[6], "gnomADe");
        assert_eq!(frequencies[26], "gnomADg_SAS");
        assert!(CacheInfo::default().variation_columns().is_none());
        assert!(!CacheInfo::default().variation_is_tabix());
    }

    #[test]
    fn test_cache_info_default() {
        let info = CacheInfo::default();
        assert_eq!(info.species, "");
        assert_eq!(info.assembly, None);
        assert!(info.cache_version.is_none());
        assert!(info.sift.is_none());
        assert!(info.polyphen.is_none());
        assert!(!info.regulatory);
        assert!(info.cell_types.is_empty());
        assert!(info.variation_cols.is_empty());
        assert!(info.serialiser_type.is_none());
        assert!(info.var_type.is_none());
        assert!(info.build.is_none());
        assert!(info.source_versions.is_empty());
        assert!(info.contig_order.is_empty());
    }

    #[test]
    fn test_cache_info_construction() {
        let mut source_versions = BTreeMap::new();
        source_versions.insert("ensembl-variation".to_string(), "115.0".to_string());

        let info = CacheInfo {
            species: "homo_sapiens".to_string(),
            assembly: Some("GRCh38".to_string()),
            cache_version: Some(115),
            sift: Some("b".to_string()),
            polyphen: Some("b".to_string()),
            regulatory: true,
            cell_types: vec!["A549".to_string(), "DND-41".to_string()],
            variation_cols: vec![
                "variation_name".to_string(),
                "failed".to_string(),
                "somatic".to_string(),
            ],
            serialiser_type: Some("sereal".to_string()),
            var_type: Some("tabix".to_string()),
            build: Some("all".to_string()),
            source_versions,
            contig_order: vec!["1".to_string(), "X".to_string(), "MT".to_string()],
        };
        assert_eq!(info.species, "homo_sapiens");
        assert_eq!(info.assembly.as_deref(), Some("GRCh38"));
        assert_eq!(info.cache_version, Some(115));
        assert!(info.regulatory);
        assert_eq!(info.cell_types.len(), 2);
        assert_eq!(
            info.source_versions.get("ensembl-variation"),
            Some(&"115.0".to_string())
        );
        assert_eq!(info.contig_order, ["1", "X", "MT"]);
    }
}
