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
    /// Genome assembly (e.g., "GRCh38").
    pub assembly: String,
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cache_info_default() {
        let info = CacheInfo::default();
        assert_eq!(info.species, "");
        assert_eq!(info.assembly, "");
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
    }

    #[test]
    fn test_cache_info_construction() {
        let mut source_versions = BTreeMap::new();
        source_versions.insert("ensembl-variation".to_string(), "115.0".to_string());

        let info = CacheInfo {
            species: "homo_sapiens".to_string(),
            assembly: "GRCh38".to_string(),
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
        };
        assert_eq!(info.species, "homo_sapiens");
        assert_eq!(info.assembly, "GRCh38");
        assert_eq!(info.cache_version, Some(115));
        assert!(info.regulatory);
        assert_eq!(info.cell_types.len(), 2);
        assert_eq!(
            info.source_versions.get("ensembl-variation"),
            Some(&"115.0".to_string())
        );
    }
}
