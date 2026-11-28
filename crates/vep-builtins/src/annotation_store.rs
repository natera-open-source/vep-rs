// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Annotation store trait: unified interface for tabix and binary backends.
//!
//! [`AnnotationStore`] abstracts over the data access pattern shared by all
//! tabix-based plugins: query a set of genomic regions and retrieve annotation
//! records indexed by position. Both [`TabixAnnotator`]
//! and [`BinaryAnnotator`] implement this trait, so plugins can be
//! backend-agnostic.
//!
//! The backend is chosen at plugin init time, by `open_backend` alone: a
//! `.vpd` store beside the table is read when [`BinaryAnnotator::open`]
//! accepts it for the table's columns and content, and the table is read
//! through tabix otherwise. Either choice is logged. A table that
//! `open_backend` or `check_table_or_store` cannot read beside a store
//! that is refused ends the plugin's load, and the refusal either returns
//! names the store and why the store is refused. Where that names the
//! conversion that rebuilds the store and the table has no tabix index, the
//! refusal ends with the flag the conversion needs on such a table.

use std::path::Path;
use std::sync::Arc;

use vep_core::assembly::Assembly;

use crate::binary_store::BinaryAnnotator;
use crate::tabix::{BatchQueryResult, MitochondrialFrame, TabixAnnotator, TabixRecord};
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

    /// The inputs of the store's frame check ([`MitochondrialFrame`]), for a
    /// store that makes one; none by default.
    #[doc(hidden)]
    fn mitochondrial_frame_mut(&mut self) -> Option<&mut MitochondrialFrame> {
        None
    }

    /// Hands the store the run's reference sequence, which its frame check
    /// compares a table's mitochondrial rows with on GRCh37
    /// ([`MitochondrialFrame`]). Without one a mitochondrial contig is read
    /// there under its own pair of names alone.
    #[doc(hidden)]
    fn set_reference_fasta(&mut self, fasta: Option<Arc<vep_fasta::IndexedFasta>>) {
        if let Some(frame) = self.mitochondrial_frame_mut() {
            frame.set_reference(fasta);
        }
    }
}

/// Try to find a pre-converted binary store (`.vpd`) for a tabix file path.
///
/// Given a path like `/path/to/cadd_snvs.tsv.gz`, checks if
/// `/path/to/cadd_snvs.tsv.vpd` exists. Returns the path if found.
pub(crate) fn find_binary_store_path(tabix_path: &Path) -> Option<std::path::PathBuf> {
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
/// store, as Ensembl VEP's modules run `add_file` before reading anything. A
/// `.vpd` store that is accepted replaces tabix and needs no `.tbi`, so the
/// check is skipped when the store's two files are there with an index of the
/// current version (`BinaryAnnotator::has_current_index`, which reads no
/// data). Any other store is refused when it is opened, so the table is
/// checked here, before the plugin reads its header, and a table refused
/// here is refused naming that store and why it is refused. A store refused
/// for its columns or its content, which only [`BinaryAnnotator::open`]
/// reads, passes here, and [`open_backend`] makes the check when it opens
/// the table.
pub(crate) fn check_table_or_store(tabix_path: &Path) -> Result<(), PluginError> {
    let Some(vpd_path) = find_binary_store_path(tabix_path) else {
        return crate::tabix::check_plugin_table(tabix_path);
    };
    match BinaryAnnotator::current_index(&vpd_path) {
        Ok(()) => Ok(()),
        Err(refusal) => crate::tabix::check_plugin_table(tabix_path)
            .map_err(|error| beside_refused_store(error, tabix_path, &vpd_path, &refusal)),
    }
}

/// What the refusal of a store ends with beside a table that has no tabix
/// index, where the refusal names the conversion that rebuilds the store:
/// how that conversion is run on such a table. The converter checks a store
/// through the table's index and ends without one unless it is given this
/// flag.
const NO_INDEX_CONVERSION: &str = "a data file with no index is converted with --no-verify";

/// `error`, the refusal of `table`, which cannot be read, naming also the
/// store at `vpd_path` and `refusal`, why that store is refused: the table is
/// all that is left to read beside a refused store, so the line that ends the
/// plugin's load says what is wrong with the store as well.
///
/// Where `refusal` names the conversion that rebuilds the store and `table`
/// stands with no index at any path the converter's check reads
/// ([`crate::tabix::find_index_path`]), the line ends with
/// [`NO_INDEX_CONVERSION`], so the conversion it names runs on the table as
/// it stands. A table whose index stands under another name than
/// `<file>.tbi` converts without the flag, and a table that is missing is
/// not converted: neither line gets the ending.
fn beside_refused_store(
    error: PluginError,
    table: &Path,
    vpd_path: &Path,
    refusal: &PluginError,
) -> PluginError {
    let named = error.with_note(&format!(
        "the store {} beside the table is refused: {}",
        vpd_path.display(),
        refusal.detail()
    ));
    let converts_with_the_flag_alone = crate::binary_store::names_the_conversion(refusal)
        && table.exists()
        && crate::tabix::find_index_path(table).is_err();
    if converts_with_the_flag_alone {
        named.with_note(NO_INDEX_CONVERSION)
    } else {
        named
    }
}

/// The query configuration of a table that has no tabix index, from the
/// `.vpd` store beside it, for a plugin that takes its coordinate columns
/// from the index before it opens its store. The store's index records the
/// columns it was built on, and [`open_backend`] accepts the store for those
/// columns only where it is bound to the table's content, so the store stands
/// in for the index as it does for a plugin that names its own columns.
/// `no_index`, the index lookup's refusal, is returned as it is where no
/// store lies beside the table, and names the store and why it is refused
/// where the store's index does not give its columns.
pub(crate) fn config_from_store_beside(
    table: &Path,
    no_index: PluginError,
) -> Result<crate::tabix::TabixConfig, PluginError> {
    let Some(vpd_path) = find_binary_store_path(table) else {
        return Err(no_index);
    };
    BinaryAnnotator::index_named_columns(&vpd_path, table)
        .map_err(|refusal| beside_refused_store(no_index, table, &vpd_path, &refusal))
}

/// What a plugin's table is read through, as [`open_backend`] chose it.
#[non_exhaustive]
pub(crate) enum Backend {
    /// The `.vpd` store beside the table, accepted for the table's columns
    /// and content.
    Store(BinaryAnnotator),
    /// The table itself, read through tabix.
    Table(TabixAnnotator),
}

/// The one choice between a table and the `.vpd` store beside it.
///
/// A store that [`BinaryAnnotator::open`] accepts is read, and an `INFO` line
/// names the table, the store and the SHA256 of each. A store it refuses
/// leaves the table, opened as a plugin's table
/// ([`TabixAnnotator::open_plugin_table`]: `<file>.tbi` must exist). A table
/// that does not open ends the plugin's load on its own refusal, which then
/// names the store and why it is refused. A table that opens is read with
/// `configure` applied to its reader: a warning carries the store's refusal
/// in its own text, without the words of a failed plugin load, since the
/// plugin loads, and an `INFO` line names
/// the table. The answer is therefore the table's whichever is read. Either
/// is opened on `assembly`, the run's, which decides the spellings a
/// mitochondrial contig is asked under ([`crate::tabix::plugin_table_candidates`]),
/// and either makes the same frame check of the same rows on GRCh37
/// ([`MitochondrialFrame`]).
pub(crate) fn open_backend(
    tabix_path: &Path,
    config: crate::tabix::TabixConfig,
    assembly: Option<Assembly>,
    configure: impl FnOnce(TabixAnnotator) -> TabixAnnotator,
) -> Result<Backend, PluginError> {
    let mut refused_store = None;
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
            Err(refusal) => refused_store = Some((vpd_path, refusal)),
        }
    }

    let table = match TabixAnnotator::open_plugin_table(config) {
        Ok(table) => table,
        Err(error) => {
            return Err(match &refused_store {
                Some((vpd_path, refusal)) => {
                    beside_refused_store(error, tabix_path, vpd_path, refusal)
                }
                None => error,
            })
        }
    };
    if let Some((vpd_path, refusal)) = &refused_store {
        tracing::warn!(
            store = %vpd_path.display(),
            table = %tabix_path.display(),
            reason = %refusal.detail(),
            "binary store refused, the table is read instead"
        );
    }
    let annotator = configure(table.on_assembly(assembly));
    tracing::info!(table = %tabix_path.display(), "plugin data read from the table");
    Ok(Backend::Table(annotator))
}

/// Open the best available annotation store for a given tabix file path.
///
/// `open_backend`'s choice behind one interface: the `.vpd` store when it
/// is accepted, else the table, reading each record line with the whitespace
/// before its line end trimmed, where a store keeps it;
/// `open_best_store_with` configures that reader.
pub fn open_best_store(
    tabix_path: &Path,
    config: crate::tabix::TabixConfig,
) -> Result<Box<dyn AnnotationStore>, PluginError> {
    open_best_store_with(tabix_path, config, None, |annotator| annotator)
}

/// [`open_best_store`] on the run's `assembly`, with `configure` applied to
/// the tabix reader when the file is read through tabix (a `.vpd` store is
/// read as it was built).
pub(crate) fn open_best_store_with(
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

    fn mitochondrial_frame_mut(&mut self) -> Option<&mut MitochondrialFrame> {
        Some(self.frame_mut())
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

    /// A table without `<table>.tbi` beside a store that is refused ends in
    /// Ensembl VEP's `.tbi` refusal, which then names the store and why it is
    /// refused. A store without its index file, or with an index of another
    /// version, earlier or later, is refused whatever it holds, and
    /// `check_table_or_store` makes the table's check itself, before a plugin
    /// reads the table's header. A store with both files and a current index
    /// is refused only by its content, here a table replaced under its name:
    /// the check passes, and `open_best_store` returns the refusal when it
    /// opens the table. Either function returns the same text for one store.
    ///
    /// A refusal that names the conversion ends with the flag the conversion
    /// needs on a table with no index. One that names no conversion does
    /// not, nor does one beside a table whose index stands on the file's
    /// stem, where the conversion runs without the flag, or beside a table
    /// that is gone.
    #[test]
    fn a_refused_store_beside_a_table_without_its_tbi_ends_in_the_index_refusal() {
        use crate::binary_store::{BinaryStoreWriter, ContentIdentity};
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        use std::sync::Arc;

        let dir = tempfile::tempdir().unwrap();
        let header = ["chrom", "pos", "name"].map(String::from).to_vec();
        let table_without_tbi = |row: &str| {
            let lines = ["#chrom\tpos\tname", row];
            let data =
                write_tabix_fixture(dir.path(), "scores.tsv.gz", &lines, IndexSpec::standard());
            std::fs::remove_file(dir.path().join("scores.tsv.gz.tbi")).unwrap();
            data
        };
        let data = table_without_tbi("21\t100\tB");
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
        let checked = || {
            check_table_or_store(&data)
                .err()
                .map_or_else(|| "the check passed".to_string(), |e| e.to_string())
        };
        let opened = || {
            open_best_store(&data, config(&data))
                .err()
                .map_or_else(|| "a backend opened".to_string(), |e| e.to_string())
        };
        assert_eq!(checked(), "the check passed");
        assert_eq!(opened(), "a backend opened");

        let store = dir.path().join("scores.tsv.vpd");
        let named = |store_refusal: &str| {
            format!(
                "plugin initialization failed: tabix index {} not found; a plugin data file \
                 needs its .tbi index beside it, so index it with tabix without -C; the store {} \
                 beside the table is refused: {store_refusal}",
                dir.path().join("scores.tsv.gz.tbi").display(),
                store.display()
            )
        };

        let vpdi = dir.path().join("scores.tsv.vpdi");
        let index = std::fs::read(&vpdi).unwrap();
        for other in [4u32, 7] {
            let mut patched = index.clone();
            patched[4..8].copy_from_slice(&other.to_le_bytes());
            std::fs::write(&vpdi, patched).unwrap();
            let refusal = named(&format!(
                "unsupported index version: {other} (expected 6); convert the table again with \
                 vep-cache-converter convert-plugin-data; a data file with no index is converted \
                 with --no-verify"
            ));
            assert_eq!(checked(), refusal);
            assert_eq!(opened(), refusal);
        }
        // The same store beside a table indexed on its stem: the conversion's
        // check reads that index, so the line ends on the conversion.
        let on_the_stem = dir.path().join("scores.tsv.tbi");
        write_tabix_fixture(
            dir.path(),
            "scores.tsv.gz",
            &["#chrom\tpos\tname", "21\t100\tB"],
            IndexSpec::standard(),
        );
        std::fs::rename(dir.path().join("scores.tsv.gz.tbi"), &on_the_stem).unwrap();
        let refusal = named(
            "unsupported index version: 7 (expected 6); convert the table again with \
             vep-cache-converter convert-plugin-data",
        );
        assert_eq!(checked(), refusal);
        assert_eq!(opened(), refusal);
        std::fs::remove_file(&on_the_stem).unwrap();

        std::fs::remove_file(&vpdi).unwrap();
        let refusal = named(&format!("binary index file not found: {}", vpdi.display()));
        assert_eq!(checked(), refusal);
        assert_eq!(opened(), refusal);

        std::fs::write(&vpdi, index).unwrap();
        // A replaced table of another length, which every open refuses.
        table_without_tbi("21\t100\tCC");
        assert_eq!(checked(), "the check passed");
        let differs = crate::binary_store::BinaryAnnotator::open(&store, &config(&data))
            .err()
            .map_or_else(|| "the store opened".to_string(), |e| e.detail());
        assert!(
            differs.starts_with("store source differs")
                && differs.ends_with("vep-cache-converter convert-plugin-data"),
            "{differs}"
        );
        assert_eq!(
            opened(),
            named(&format!(
                "{differs}; a data file with no index is converted with --no-verify"
            ))
        );

        // A table that is gone, beside the store of another index version:
        // nothing is left to convert, so the line ends on the conversion.
        let mut patched = std::fs::read(&vpdi).unwrap();
        patched[4..8].copy_from_slice(&7u32.to_le_bytes());
        std::fs::write(&vpdi, patched).unwrap();
        std::fs::remove_file(&data).unwrap();
        assert_eq!(
            checked(),
            format!(
                "plugin initialization failed: data file {} not found; the store {} beside the \
                 table is refused: unsupported index version: 7 (expected 6); convert the table \
                 again with vep-cache-converter convert-plugin-data",
                data.display(),
                store.display()
            )
        );
    }

    /// A plugin that takes its coordinate columns from its table's index
    /// takes them from the store where the table has none:
    /// `indexed_table_config` gives the columns the store was built on, an
    /// end column and a zero-based start among them, and `open_best_store`
    /// reads the store with them. With no store beside the table the index
    /// lookup's refusal stands alone, and beside a store whose index is cut
    /// short inside its columns it names the store and why it is refused.
    #[test]
    fn a_table_with_no_index_takes_its_columns_from_the_store_beside_it() {
        use crate::binary_store::{BinaryStoreWriter, ContentIdentity};
        use crate::tabix::indexed_table_config;
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        use std::sync::Arc;

        let dir = tempfile::tempdir().unwrap();
        let lines = ["#chrom\tstart\tend\tname", "21\t99\t120\tB"];
        let spans = IndexSpec::standard().with_span_cols(1, 2, 3).zero_based();
        let data = write_tabix_fixture(dir.path(), "spans.tsv.gz", &lines, spans);
        std::fs::remove_file(dir.path().join("spans.tsv.gz.tbi")).unwrap();
        let columns = || indexed_table_config(&data, |index| index.display().to_string());
        let refusal = || match columns() {
            Err(error) => error.to_string(),
            Ok(_) => "the table gave its columns".to_string(),
        };
        let no_index = format!(
            "plugin initialization failed: tabix index not found for {table}. Expected \
             {table}.csi or {table}.tbi",
            table = data.display()
        );
        assert_eq!(refusal(), no_index);

        let header = ["chrom", "start", "end", "name"].map(String::from).to_vec();
        let record = TabixRecord {
            chr: "21".to_string(),
            start: 100,
            end: 120,
            columns: ["21", "99", "120", "B"].map(String::from).to_vec(),
            header: Some(Arc::new(header.clone())),
        };
        let built_on = crate::tabix::TabixConfig {
            end_col: Some(2),
            zero_based: true,
            ..config(&data)
        };
        let store = dir.path().join("spans.tsv.vpd");
        BinaryStoreWriter::new(store.clone())
            .write(
                &header,
                &[("21".to_string(), vec![record])],
                &built_on,
                &ContentIdentity::of_file(&data).unwrap(),
            )
            .unwrap();

        let Ok(read) = columns() else {
            panic!("{}", refusal())
        };
        assert_eq!(
            (
                read.file_path.as_path(),
                read.chr_col,
                read.start_col,
                read.end_col,
                read.ref_col,
                read.alt_col,
                read.zero_based
            ),
            (data.as_path(), 0, 1, Some(2), None, None, true)
        );
        let answers = open_best_store(&data, read).unwrap();
        assert_eq!(
            answers.query("21", 110, 110).unwrap()[0].columns,
            ["21", "99", "120", "B"]
        );

        let vpdi = dir.path().join("spans.tsv.vpdi");
        let index = std::fs::read(&vpdi).unwrap();
        std::fs::write(&vpdi, &index[..4 + 4 + 2]).unwrap();
        assert_eq!(
            refusal(),
            format!(
                "{no_index}; the store {} beside the table is refused: failed to read binary \
                 store: failed to fill whole buffer",
                store.display()
            )
        );
    }
}
