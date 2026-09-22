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
//! Format selection happens at plugin init time: if a `.vpd` file exists
//! alongside the original `.tsv.gz`, the binary backend is used; otherwise
//! the tabix backend is the fallback.

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
/// If a pre-converted `.vpd` file exists, opens a [`BinaryAnnotator`].
/// Otherwise, falls back to [`TabixAnnotator`].
///
/// This function is used by plugins in their `init()` to transparently
/// select the fastest available backend.
pub fn open_best_store(
    tabix_path: &Path,
    config: crate::tabix::TabixConfig,
) -> Result<Box<dyn AnnotationStore>, PluginError> {
    if let Some(vpd_path) = find_binary_store_path(tabix_path) {
        match crate::binary_store::BinaryAnnotator::open(&vpd_path) {
            Ok(store) => {
                tracing::info!(
                    vpd = %vpd_path.display(),
                    "using binary annotation store"
                );
                return Ok(Box::new(store));
            }
            Err(e) => {
                tracing::warn!(
                    vpd = %vpd_path.display(),
                    error = %e,
                    "failed to open binary store, falling back to tabix"
                );
            }
        }
    }

    let annotator = crate::tabix::TabixAnnotator::open(config)?;
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
    use std::path::PathBuf;

    #[test]
    fn find_binary_store_path_returns_none_for_missing() {
        let path = PathBuf::from("/nonexistent/data/cadd_snvs.tsv.gz");
        assert!(find_binary_store_path(&path).is_none());
    }
}
