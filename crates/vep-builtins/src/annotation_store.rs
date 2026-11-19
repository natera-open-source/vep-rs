// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Annotation store trait: unified interface for tabix and binary backends.
//!
//! [`AnnotationStore`] abstracts over the data access pattern shared by all
//! tabix-based plugins: query a set of genomic regions and retrieve annotation
//! records indexed by position. Both [`TabixAnnotator`](crate::tabix::TabixAnnotator)
//! and [`BinaryAnnotator`] implement this trait, so plugins can be
//! backend-agnostic.
//!
//! The backend is chosen at plugin init time, by [`open_backend`] alone: a
//! `.vpd` store beside the table is read when [`BinaryAnnotator::open`]
//! accepts it for the table's columns and content, and the table is read
//! through tabix otherwise. Either choice is logged.

use std::path::Path;

use vep_core::assembly::Assembly;

use crate::binary_store::BinaryAnnotator;
use crate::tabix::{BatchQueryResult, TabixAnnotator, TabixRecord};
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

/// Ensembl VEP's check on a plugin's data table ([`crate::tabix::check_plugin_table`])
/// for a plugin that reads the table's header or index before it opens its
/// store, as Ensembl VEP's modules run `add_file` before reading anything: a
/// `.vpd` store that is accepted replaces tabix and needs no `.tbi`, so the
/// check is skipped when a store file is there, and [`open_backend`] makes it
/// when it refuses the store and opens the table.
pub fn check_table_or_store(tabix_path: &Path) -> Result<(), PluginError> {
    if find_binary_store_path(tabix_path).is_none() {
        crate::tabix::check_plugin_table(tabix_path)?;
    }
    Ok(())
}

/// What a plugin's table is read through, as [`open_backend`] chose it.
pub enum Backend {
    /// The `.vpd` store beside the table, accepted for the table's columns
    /// and content.
    Store(BinaryAnnotator),
    /// The table itself, read through tabix.
    Table(TabixAnnotator),
}

/// The one choice between a table and the `.vpd` store beside it.
///
/// A store that [`BinaryAnnotator::open`] accepts is read, and an `INFO` line
/// names the table, the store and the SHA256 of each. A store it refuses is
/// left with a warning carrying the refusal, and the table is read: opened as
/// a plugin's table ([`TabixAnnotator::open_plugin_table`]: `<file>.tbi` must
/// exist) with `configure` applied to its reader, under an `INFO` line naming
/// the table. The answer is therefore the table's whichever is read. Either
/// is opened on `assembly`, the run's, which decides the spellings a
/// mitochondrial contig is asked under ([`crate::tabix::plugin_table_candidates`]).
pub fn open_backend(
    tabix_path: &Path,
    config: crate::tabix::TabixConfig,
    assembly: Option<Assembly>,
    configure: impl FnOnce(TabixAnnotator) -> TabixAnnotator,
) -> Result<Backend, PluginError> {
    if let Some(vpd_path) = find_binary_store_path(tabix_path) {
        match BinaryAnnotator::open(&vpd_path, &config) {
            Ok(store) => {
                tracing::info!(
                    table = %tabix_path.display(),
                    store = %vpd_path.display(),
                    table_sha256 = %store.source().sha256_hex(),
                    store_sha256 = %store.data().sha256_hex(),
                    "plugin data read from the binary store"
                );
                return Ok(Backend::Store(store.on_assembly(assembly)));
            }
            Err(e) => {
                tracing::warn!(
                    store = %vpd_path.display(),
                    table = %tabix_path.display(),
                    reason = %e,
                    "binary store refused, the table is read instead"
                );
            }
        }
    }

    let annotator = configure(TabixAnnotator::open_plugin_table(config)?.on_assembly(assembly));
    tracing::info!(table = %tabix_path.display(), "plugin data read from the table");
    Ok(Backend::Table(annotator))
}

/// Open the best available annotation store for a given tabix file path.
///
/// [`open_backend`]'s choice behind one interface: the `.vpd` store when it
/// is accepted, else the table, reading each record line with the whitespace
/// before its line end trimmed, where a store keeps it;
/// [`open_best_store_with`] configures that reader.
pub fn open_best_store(
    tabix_path: &Path,
    config: crate::tabix::TabixConfig,
) -> Result<Box<dyn AnnotationStore>, PluginError> {
    open_best_store_with(tabix_path, config, None, |annotator| annotator)
}

/// [`open_best_store`] on the run's `assembly`, with `configure` applied to
/// the tabix reader when the file is read through tabix (a `.vpd` store is
/// read as it was built).
pub fn open_best_store_with(
    tabix_path: &Path,
    config: crate::tabix::TabixConfig,
    assembly: Option<Assembly>,
    configure: impl FnOnce(TabixAnnotator) -> TabixAnnotator,
) -> Result<Box<dyn AnnotationStore>, PluginError> {
    Ok(
        match open_backend(tabix_path, config, assembly, configure)? {
            Backend::Store(store) => Box::new(store),
            Backend::Table(table) => Box::new(table),
        },
    )
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
    use std::path::PathBuf;

    #[test]
    fn find_binary_store_path_returns_none_for_missing() {
        let path = PathBuf::from("/nonexistent/data/cadd_snvs_synthetic.tsv.gz");
        assert!(find_binary_store_path(&path).is_none());
    }

    fn config(data: &Path) -> crate::tabix::TabixConfig {
        crate::tabix::TabixConfig {
            file_path: data.to_path_buf(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        }
    }

    /// A plugin table without `<file>.tbi` beside it is refused as Ensembl VEP's
    /// tabix plugins refuse it, whatever other index htslib would find.
    #[test]
    fn tabix_fallback_needs_the_tables_tbi() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = ["#chrom\tpos\tname", "21\t99\tA", "21\t100\tB"];
        let data = write_tabix_fixture(dir.path(), "scores.tsv.gz", &lines, IndexSpec::standard());
        assert!(open_best_store(&data, config(&data)).is_ok());

        std::fs::rename(
            dir.path().join("scores.tsv.gz.tbi"),
            dir.path().join("scores.tsv.tbi"),
        )
        .unwrap();
        let refusal = open_best_store(&data, config(&data))
            .err()
            .map_or_else(|| "the table opened".to_string(), |e| e.to_string());
        assert!(refusal.contains("scores.tsv.gz.tbi not found"), "{refusal}");
    }

    /// A `.vpd` store beside the table is the documented exception: it is read
    /// in place of tabix, so the table needs no `.tbi` at all. The store holds
    /// the table's own row and keeps the chromosome in its cells, where the
    /// tabix reader moves it out, so the first cell tells the store answered.
    #[test]
    fn a_vpd_store_beside_the_table_needs_no_tbi() {
        use crate::binary_store::{BinaryStoreWriter, ContentIdentity};
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        use std::sync::Arc;

        let dir = tempfile::tempdir().unwrap();
        let header = ["chrom", "pos", "name"].map(String::from).to_vec();
        let lines = ["#chrom\tpos\tname", "21\t100\tB"];
        let data = write_tabix_fixture(dir.path(), "scores.tsv.gz", &lines, IndexSpec::standard());
        std::fs::remove_file(dir.path().join("scores.tsv.gz.tbi")).unwrap();
        let record = TabixRecord {
            chr: "21".to_string(),
            start: 100,
            end: 100,
            columns: ["21", "100", "B"].map(String::from).to_vec(),
            header: Some(Arc::new(header.clone())),
        };
        BinaryStoreWriter::new(dir.path().join("scores.tsv.vpd"))
            .write(
                &header,
                &[("21".to_string(), vec![record])],
                &config(&data),
                &ContentIdentity::of_file(&data).unwrap(),
            )
            .unwrap();

        let store = open_best_store(&data, config(&data)).unwrap();
        assert_eq!(
            store.query("21", 100, 100).unwrap()[0].columns,
            ["21", "100", "B"]
        );
    }
}
