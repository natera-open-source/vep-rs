// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VEP cache metadata from info.txt files.

use std::collections::BTreeMap;

/// The regulatory shard shape this build reads: per regulatory feature the
/// stable id, span, feature type, database id and activity codes; per motif
/// feature the stable id, span, strand, database id, score, stranded sequence
/// and the stable id of its binding matrix; the matrices once in
/// `regulatory/binding_matrices.json`.
pub const REGULATORY_SCHEMA: u32 = 1;

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
    /// The epigenomes of the regulatory build, as the converter lists them: the
    /// sorted union of every regulatory feature's activity-map keys (epigenome
    /// short names with spaces as `_`, `Database/RegFeat.pm:196-203` at
    /// release/116.2). A regulatory feature's `activity` string holds one code
    /// per name in this order. Ensembl VEP's own `cell_types` line is a
    /// comma-joined string the names cannot be split back from (some contain
    /// commas), so a cache converted without the list has none.
    pub cell_types: Vec<String>,
    /// The shape of the regulatory shards and of
    /// `regulatory/binding_matrices.json`, written by the converter beside
    /// `regulatory`; the current shape is [`REGULATORY_SCHEMA`], and a cache
    /// whose regulatory tree declares none holds an unversioned shape. An
    /// absent key stays absent when the object is written.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub regulatory_schema: Option<u32>,
    /// Column names for the variation cache TSV files.
    pub variation_cols: Vec<String>,
    /// Serializer type: "sereal" or None for Storable.
    pub serialiser_type: Option<String>,
    /// Variation file type: "tabix" for tabix-indexed, None for flat files.
    pub var_type: Option<String>,
    /// Build identifier.
    pub build: Option<String>,
    /// The alignment file the transcript models were edited against, the `bam`
    /// line of a RefSeq or merged cache's `info.txt` (`Dumper/Otherfeatures.pm:80`
    /// at release/116.2); a cache with the key switches `--use_transcript_ref`
    /// on (`CacheDir.pm:412-416`). An absent key stays absent when the object
    /// is written.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bam: Option<String>,
    /// The transcript set the cache holds, `ensembl`, `refseq` or `merged`, from
    /// the species directory the cache was converted from: a key of this cache
    /// layout with no counterpart in Ensembl VEP's `info.txt`, `None` when the
    /// cache's metadata lacks the key.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_type: Option<String>,
    /// Source database versions (e.g., "ClinVar" -> "202306"), in the byte
    /// order the output headers print them.
    pub source_versions: BTreeMap<String, String>,
    /// The cache's chromosome directories in the order the converter and the
    /// builder list them: autosomes by number, X, Y, MT, then every other
    /// contig by name. Empty when the cache's metadata lacks the key.
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

    /// Two keys of an Ensembl VEP cache's `info.txt`, carried into a converted
    /// cache's `info.json`, change nothing here: `cache_region_size`, the span
    /// of one Storable region file (`CacheDir.pm:149` reads it, default
    /// `Config.pm:296`'s 1000000), where a JSON shard's span is its own; and
    /// `build`, which `CacheDir.pm:454-483` stores and no module reads. The
    /// object parses to the same metadata with and without `cache_region_size`,
    /// and `build` is kept as the metadata it is.
    #[test]
    fn region_size_and_build_keys_are_accepted() {
        let with: CacheInfo = serde_json::from_str(
            r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116,"cache_region_size":1000000,"build":"all"}"#,
        )
        .unwrap();
        let without: CacheInfo = serde_json::from_str(
            r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116,"build":"all"}"#,
        )
        .unwrap();
        assert_eq!(format!("{with:?}"), format!("{without:?}"));
        assert_eq!(with.build.as_deref(), Some("all"));
        assert_eq!(with.cache_version, Some(116));
        assert!(!serde_json::to_string(&with)
            .unwrap()
            .contains("cache_region_size"));
    }

    /// The `bam` line of a RefSeq cache's `info.txt` and the converter's
    /// `source_type` key are read and written as the strings they are; a cache
    /// without either key (every Ensembl cache, every converted cache whose
    /// metadata lacks the key) parses to `None` and writes neither key back.
    #[test]
    fn bam_line_round_trips_through_info_json() {
        let edited: CacheInfo = serde_json::from_str(
            r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":84,
                "bam":"interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam","source_type":"refseq"}"#,
        )
        .unwrap();
        assert_eq!(
            edited.bam.as_deref(),
            Some("interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam")
        );
        assert_eq!(edited.source_type.as_deref(), Some("refseq"));
        let written = serde_json::to_string(&edited).unwrap();
        let again: CacheInfo = serde_json::from_str(&written).unwrap();
        assert_eq!(again.bam, edited.bam);
        assert_eq!(again.source_type, edited.source_type);

        let plain: CacheInfo = serde_json::from_str(
            r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116}"#,
        )
        .unwrap();
        assert!(plain.bam.is_none());
        assert!(plain.source_type.is_none());
        let written = serde_json::to_string(&plain).unwrap();
        assert!(!written.contains("\"bam\""), "{written}");
        assert!(!written.contains("source_type"), "{written}");
    }

    /// `cell_types` is a JSON array of the epigenome names, so a name holding a
    /// comma (77 of the 112 GRCh38 names) survives a write and a read intact and
    /// in place, where Ensembl VEP's comma-joined `cell_types` line cannot be
    /// split back into them (CacheDir.pm:174); an `info.json` without the key
    /// reads as an empty list.
    #[test]
    fn cell_types_round_trip() {
        let names = [
            "A549_(m,_58_y)",
            "GM12878_(f)",
            "common_myeloid_progenitor,_CD34-positive_(27-36_y)",
            "astrocyte",
        ];
        let info = CacheInfo {
            regulatory: true,
            regulatory_schema: Some(REGULATORY_SCHEMA),
            cell_types: names.map(str::to_string).to_vec(),
            ..Default::default()
        };
        let written = serde_json::to_string(&info).unwrap();
        assert!(
            written.contains(
                r#""cell_types":["A549_(m,_58_y)","GM12878_(f)","common_myeloid_progenitor,_CD34-positive_(27-36_y)","astrocyte"]"#
            ),
            "{written}"
        );
        let again: CacheInfo = serde_json::from_str(&written).unwrap();
        assert_eq!(again.cell_types, names);
        assert_eq!(again.cell_types.join(",").split(',').count(), 6);

        let without: CacheInfo =
            serde_json::from_str(r#"{"species":"homo_sapiens","regulatory":true}"#).unwrap();
        assert!(without.cell_types.is_empty());
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
        assert!(info.regulatory_schema.is_none());
        assert!(info.variation_cols.is_empty());
        assert!(info.serialiser_type.is_none());
        assert!(info.var_type.is_none());
        assert!(info.build.is_none());
        assert!(info.bam.is_none());
        assert!(info.source_type.is_none());
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
            regulatory_schema: Some(REGULATORY_SCHEMA),
            variation_cols: vec![
                "variation_name".to_string(),
                "failed".to_string(),
                "somatic".to_string(),
            ],
            serialiser_type: Some("sereal".to_string()),
            var_type: Some("tabix".to_string()),
            build: Some("all".to_string()),
            bam: None,
            source_type: Some("ensembl".to_string()),
            source_versions,
            contig_order: vec!["1".to_string(), "X".to_string(), "MT".to_string()],
        };
        assert_eq!(info.species, "homo_sapiens");
        assert_eq!(info.assembly.as_deref(), Some("GRCh38"));
        assert_eq!(info.cache_version, Some(115));
        assert!(info.regulatory);
        assert_eq!(info.cell_types.len(), 2);
        assert_eq!(info.regulatory_schema, Some(1));
        assert_eq!(
            info.source_versions.get("ensembl-variation"),
            Some(&"115.0".to_string())
        );
        assert_eq!(info.contig_order, ["1", "X", "MT"]);
    }
}
