// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Annotation store trait: unified interface for tabix and binary backends.
//!
//! [`AnnotationStore`] abstracts over the data access pattern shared by all
//! tabix-based plugins: query a set of genomic regions and retrieve annotation
//! records indexed by position. Both [`TabixAnnotator`](crate::tabix::TabixAnnotator)
//! and [`BinaryAnnotator`](crate::binary_store::BinaryAnnotator) implement this
//! trait, so plugins can be backend-agnostic.
//!
//! Format selection happens at plugin init time: a `.vpd` sibling is selected
//! only when its provenance verifies the source, settings, data and index.
//! Otherwise the tabix backend is the fallback.

use std::path::Path;

use crate::tabix::{BatchQueryResult, TabixRecord};
use crate::PluginError;

/// Unified interface for querying annotation data files.
///
/// Implementors provide region-based queries that return [`TabixRecord`]s
/// indexed by position. The [`query_batch`](AnnotationStore::query_batch)
/// method is the primary hot path used by plugin prefetch.
pub trait AnnotationStore: Send + Sync {
    /// Batch query: accepts a slice of `(chr, start, end)` regions, merges
    /// overlapping regions before querying, and returns all matching records
    /// indexed by position for fast per-variant lookup.
    fn query_batch(&self, regions: &[(&str, u64, u64)]) -> Result<BatchQueryResult, PluginError>;

    /// Single-region query convenience method.
    fn query(&self, chr: &str, start: u64, end: u64) -> Result<Vec<TabixRecord>, PluginError>;

    /// Column header names from the data file, if available.
    fn header(&self) -> Option<&[String]>;
}

/// Try to find a pre-converted binary store (`.vpd`) for a tabix file path.
///
/// Given a path like `/path/to/cadd_snvs.tsv.gz`, checks if
/// `/path/to/cadd_snvs.tsv.vpd` exists. Returns the path if found.
pub fn find_binary_store_path(tabix_path: &Path) -> Option<std::path::PathBuf> {
    let stem = tabix_path.to_str()?;
    let base = stem
        .strip_suffix(".gz")
        .or_else(|| stem.strip_suffix(".bgz"))
        .unwrap_or(stem);
    let vpd_path = std::path::PathBuf::from(format!("{base}.vpd"));
    if vpd_path.exists() {
        Some(vpd_path)
    } else {
        None
    }
}

/// Open the best available annotation store for a given tabix file path.
///
/// If a pre-converted `.vpd` file has matching content and conversion
/// provenance, opens a [`crate::BinaryAnnotator`]. Otherwise, warns and falls
/// back to [`crate::TabixAnnotator`].
///
/// Plugins call this during `init()` to select a verified binary sibling
/// without changing their configuration.
pub fn open_best_store(
    tabix_path: &Path,
    config: crate::tabix::TabixConfig,
) -> Result<Box<dyn AnnotationStore>, PluginError> {
    if config.file_path != tabix_path {
        return Err(PluginError::Init(format!(
            "annotation source {} differs from configured path {}",
            tabix_path.display(),
            config.file_path.display()
        )));
    }
    if let Some(vpd_path) = find_binary_store_path(tabix_path) {
        let verified =
            crate::provenance::Provenance::load(&vpd_path, &config).and_then(|identity| {
                let store = crate::binary_store::BinaryAnnotator::open_checked(
                    &vpd_path,
                    |data, index| identity.verify(tabix_path, data, index),
                )?;
                Ok((store, identity))
            });
        match verified {
            Ok((store, identity)) => {
                tracing::info!(
                    backend = "binary",
                    source = %tabix_path.display(),
                    vpd = %vpd_path.display(),
                    source_sha256 = %identity.source_sha256,
                    data_sha256 = %identity.data_sha256,
                    index_sha256 = %identity.index_sha256,
                    "using verified binary annotation store"
                );
                return Ok(Box::new(store));
            }
            Err(e) => {
                tracing::warn!(
                    source = %tabix_path.display(),
                    vpd = %vpd_path.display(),
                    error = %e,
                    "binary annotation identity is unverifiable, falling back to requested tabix source"
                );
            }
        }
    }

    let annotator = crate::tabix::TabixAnnotator::open(config)?;
    tracing::info!(
        backend = "tabix",
        source = %tabix_path.display(),
        identity = "requested source; no binary substitution",
        "using tabix annotation store"
    );
    Ok(Box::new(annotator))
}

impl AnnotationStore for crate::tabix::TabixAnnotator {
    fn query_batch(&self, regions: &[(&str, u64, u64)]) -> Result<BatchQueryResult, PluginError> {
        self.query_batch(regions)
    }

    fn query(&self, chr: &str, start: u64, end: u64) -> Result<Vec<TabixRecord>, PluginError> {
        self.query(chr, start, end)
    }

    fn header(&self) -> Option<&[String]> {
        self.header()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_store::BinaryStoreWriter;
    use crate::provenance::{manifest_path, ConversionSettings, Publication};
    use crate::tabix::TabixConfig;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use std::fs;
    use std::path::PathBuf;

    fn fixture(dir: &Path, name: &str, score: &str, position: u64) -> TabixConfig {
        let file_path = write_tabix_fixture(
            dir,
            name,
            &["#chrom\tpos\tscore", &format!("21\t{position}\t{score}")],
            IndexSpec::standard(),
        );
        TabixConfig {
            file_path,
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        }
    }

    fn convert(config: &TabixConfig) -> PathBuf {
        let output = config.file_path.with_extension("vpd");
        let publication =
            Publication::begin(&config.file_path, &output, ConversionSettings::from(config))
                .unwrap();
        let digest = publication
            .source_reader()
            .unwrap()
            .digest_to_end()
            .unwrap();
        let tabix = crate::TabixAnnotator::open(config.clone()).unwrap();
        let mut records = tabix.query("21", 1, 1000).unwrap();
        for record in &mut records {
            record.columns[config.chr_col] = record.chr.clone();
        }
        BinaryStoreWriter::new(publication.data_path())
            .write(tabix.header().unwrap(), &[("21".into(), records)])
            .unwrap();
        publication
            .publish(digest, ConversionSettings::from(config))
            .unwrap();
        output
    }

    fn assert_requested_source(config: &TabixConfig, score: &str) {
        let store = open_best_store(&config.file_path, config.clone()).unwrap();
        let records = store.query("21", 200, 200).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].get("score"), Some(score));
        // Tabix moves the chromosome field into `chr`; binary keeps its column.
        assert!(records[0].columns[config.chr_col].is_empty());
    }

    #[test]
    fn find_binary_store_path_returns_none_for_missing() {
        let path = PathBuf::from("/nonexistent/data/cadd_snvs.tsv.gz");
        assert!(find_binary_store_path(&path).is_none());
    }

    #[test]
    fn verified_sibling_works_without_the_tabix_index() {
        let dir = tempfile::tempdir().unwrap();
        let config = fixture(dir.path(), "data.tsv.gz", "SRC_A", 200);
        convert(&config);
        fs::remove_file(dir.path().join("data.tsv.gz.tbi")).unwrap();
        let store = open_best_store(&config.file_path, config.clone()).unwrap();
        let records = store.query("21", 200, 200).unwrap();
        assert_eq!(records[0].get("score"), Some("SRC_A"));
        assert_eq!(records[0].columns[0], "21");
    }

    #[test]
    fn source_change_with_same_length_and_timestamp_falls_back() {
        let dir = tempfile::tempdir().unwrap();
        let config = fixture(dir.path(), "data.tsv.gz", "SRC_A", 200);
        convert(&config);
        let old = fs::metadata(&config.file_path).unwrap();
        fixture(dir.path(), "data.tsv.gz", "SRC_B", 200);
        let file = fs::File::options()
            .write(true)
            .open(&config.file_path)
            .unwrap();
        file.set_times(fs::FileTimes::new().set_modified(old.modified().unwrap()))
            .unwrap();
        assert_eq!(fs::metadata(&config.file_path).unwrap().len(), old.len());
        assert_eq!(
            fs::metadata(&config.file_path).unwrap().modified().unwrap(),
            old.modified().unwrap()
        );
        assert_requested_source(&config, "SRC_B");
    }

    #[test]
    fn swapping_data_index_or_both_falls_back() {
        for extensions in [vec!["vpd"], vec!["vpdi"], vec!["vpd", "vpdi"]] {
            let dir = tempfile::tempdir().unwrap();
            let config = fixture(dir.path(), "data.tsv.gz", "SRC_A", 200);
            let output = convert(&config);
            let other = fixture(dir.path(), "other.tsv.gz", "FOREIGN_SCORE", 300);
            let other_output = convert(&other);
            for extension in extensions {
                fs::copy(
                    other_output.with_extension(extension),
                    output.with_extension(extension),
                )
                .unwrap();
            }
            assert_requested_source(&config, "SRC_A");
        }
    }

    #[test]
    fn missing_legacy_or_corrupt_metadata_falls_back() {
        let dir = tempfile::tempdir().unwrap();
        let config = fixture(dir.path(), "data.tsv.gz", "SRC_A", 200);
        let output = convert(&config);
        let path = manifest_path(&output);
        let good = fs::read(&path).unwrap();
        fs::remove_file(&path).unwrap();
        assert_requested_source(&config, "SRC_A");
        for corrupt in [b"{".to_vec(), b"{}".to_vec(), vec![b' '; 16385]] {
            fs::write(&path, corrupt).unwrap();
            assert_requested_source(&config, "SRC_A");
        }
        for (key, value) in [
            ("version", serde_json::json!(999)),
            ("source_sha256", serde_json::json!("invalid")),
            ("data_sha256", serde_json::json!("0".repeat(64))),
            ("index_sha256", serde_json::json!("0".repeat(64))),
        ] {
            let mut manifest: serde_json::Value = serde_json::from_slice(&good).unwrap();
            manifest[key] = value;
            fs::write(&path, serde_json::to_vec(&manifest).unwrap()).unwrap();
            assert_requested_source(&config, "SRC_A");
        }
    }

    #[test]
    fn every_coordinate_setting_mismatch_falls_back() {
        let dir = tempfile::tempdir().unwrap();
        let config = fixture(dir.path(), "data.tsv.gz", "SRC_A", 200);
        let output = convert(&config);
        let path = manifest_path(&output);
        let good = fs::read(&path).unwrap();
        for (key, value) in [
            ("chr_col", serde_json::json!(2)),
            ("start_col", serde_json::json!(2)),
            ("end_col", serde_json::json!(1)),
            ("zero_based", serde_json::json!(true)),
        ] {
            let mut manifest: serde_json::Value = serde_json::from_slice(&good).unwrap();
            manifest["settings"][key] = value;
            fs::write(&path, serde_json::to_vec(&manifest).unwrap()).unwrap();
            assert_requested_source(&config, "SRC_A");
        }
    }

    #[test]
    fn a_different_configured_source_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let config = fixture(dir.path(), "data.tsv.gz", "SRC_A", 200);
        let error = open_best_store(&dir.path().join("other.tsv.gz"), config)
            .err()
            .unwrap();
        assert!(error.to_string().contains("differs from configured path"));
    }

    #[test]
    fn missing_source_cannot_be_replaced_by_a_verified_sibling() {
        let dir = tempfile::tempdir().unwrap();
        let config = fixture(dir.path(), "data.tsv.gz", "SRC_A", 200);
        convert(&config);
        fs::remove_file(&config.file_path).unwrap();
        assert!(open_best_store(&config.file_path, config.clone()).is_err());
    }

    #[test]
    fn verification_uses_the_same_bytes_as_the_opened_store() {
        let dir = tempfile::tempdir().unwrap();
        let config = fixture(dir.path(), "data.tsv.gz", "SRC_A", 200);
        let output = convert(&config);
        let other = fixture(dir.path(), "other.tsv.gz", "FOREIGN_SCORE", 300);
        let other_output = convert(&other);
        let identity = crate::provenance::Provenance::load(&output, &config).unwrap();
        let store = crate::BinaryAnnotator::open_checked(&output, |data, index| {
            fs::rename(&other_output, &output).unwrap();
            fs::rename(
                other_output.with_extension("vpdi"),
                output.with_extension("vpdi"),
            )
            .unwrap();
            identity.verify(&config.file_path, data, index)
        })
        .unwrap();
        assert_eq!(
            store.query("21", 200, 200).unwrap()[0].get("score"),
            Some("SRC_A")
        );
        assert!(store.query("21", 300, 300).unwrap().is_empty());
    }
}
