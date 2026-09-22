// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! dbscSNV plugin: tabix-based splice site variant predictions.
//!
//! Queries a tabix-indexed dbscSNV file and annotates transcripts with
//! ADA and RF splice-site prediction scores. These two ensemble methods predict
//! whether SNVs at splice sites affect splicing.
//!
//! **Usage:** `--plugin dbscSNV,/path/to/dbscSNV.txt.gz`
//!
//! Output fields:
//! - `dbscSNV_ADA_SCORE`: AdaBoost splice-site prediction score
//! - `dbscSNV_RF_SCORE`: Random Forest splice-site prediction score
//!
//! Two-phase plugin: [`BuiltinPlugin::prefetch`] batches the tabix query for a
//! buffer and [`BuiltinPlugin::annotate`] writes per-transcript results.

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{open_best_store, AnnotationStore};
use crate::tabix::{allele_match, BatchQueryResult, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// Standard column layout in dbscSNV files.
const CHR_COL: usize = 0;
const POS_COL: usize = 1;
const REF_COL: usize = 2;
const ALT_COL: usize = 3;

/// Field names for output annotation.
const ADA_FIELD: &str = "dbscSNV_ADA_SCORE";
const RF_FIELD: &str = "dbscSNV_RF_SCORE";

/// Column names in the dbscSNV file header to extract scores from.
const ADA_HEADER_NAME: &str = "ada_score";
const RF_HEADER_NAME: &str = "rf_score";

/// The GRCh38 file names its coordinate columns `hg38_chr` / `hg38_pos`; the
/// GRCh37 file uses plain `chr` / `pos`. That difference is the only in-file
/// signal of which build a dbscSNV copy holds.
const GRCH38_CHR_HEADER_NAME: &str = "hg38_chr";

/// dbscSNV splice-site prediction plugin.
pub struct DbscSnvPlugin {
    /// Annotation store, initialized in `init()`.
    store: Option<Box<dyn AnnotationStore>>,
    /// Assembly from `--assembly`, injected before `init()`.
    assembly: Option<vep_core::assembly::Assembly>,
}

impl DbscSnvPlugin {
    /// Create a new uninitialized dbscSNV plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            assembly: None,
        }
    }

    /// Which assembly a dbscSNV file holds, from its coordinate column names.
    ///
    /// The two upstream files are pre-split per assembly and are otherwise
    /// identical in shape, so pointing a run at the wrong one produces a silent
    /// zero: the coordinates simply do not exist in the other build.
    fn assembly_from_header(header: &[String]) -> Option<vep_core::assembly::Assembly> {
        let first = header.first()?.trim().trim_start_matches('#');
        if first.eq_ignore_ascii_case(GRCH38_CHR_HEADER_NAME) {
            Some(vep_core::assembly::Assembly::Grch38)
        } else if first.eq_ignore_ascii_case("chr") {
            Some(vep_core::assembly::Assembly::Grch37)
        } else {
            None
        }
    }
}

impl Default for DbscSnvPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl BuiltinPlugin for DbscSnvPlugin {
    fn name(&self) -> &str {
        "dbscSNV"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            (
                ADA_FIELD.to_string(),
                "dbscSNV AdaBoost splice-site prediction score".to_string(),
            ),
            (
                RF_FIELD.to_string(),
                "dbscSNV Random Forest splice-site prediction score".to_string(),
            ),
        ]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let file_path = params
            .first()
            .ok_or_else(|| PluginError::Init("dbscSNV requires a file path argument".into()))?;

        let config = TabixConfig {
            file_path: file_path.into(),
            chr_col: CHR_COL,
            start_col: POS_COL,
            end_col: None,
            ref_col: Some(REF_COL),
            alt_col: Some(ALT_COL),
            zero_based: false,
        };

        let store = open_best_store(std::path::Path::new(file_path.as_str()), config)?;

        if let (Some(requested), Some(header)) = (self.assembly, store.header()) {
            if let Some(found) = Self::assembly_from_header(header) {
                if found != requested {
                    return Err(PluginError::Init(format!(
                        "dbscSNV: --assembly is {requested} but '{file_path}' is a {found} file \
                         (its first column is '{}'). dbscSNV publishes one file per assembly; \
                         point --plugin dbscSNV at the {requested} file.",
                        header.first().map(String::as_str).unwrap_or("?"),
                    )));
                }
            }
        }

        self.store = Some(store);
        debug!("dbscSNV plugin initialized with {}", file_path);
        Ok(())
    }

    fn set_assembly(&mut self, assembly: Option<vep_core::assembly::Assembly>) {
        self.assembly = assembly;
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        if variants.is_empty() {
            return Ok(None);
        }

        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .map(|v| (v.chr.as_str(), v.start, v.end))
            .collect();
        let batch = store.query_batch(&regions)?;
        Ok(Some(PrefetchData::new(batch)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch: BatchQueryResult = data.downcast()?;

        for variant in variants.iter_mut() {
            let records = batch.lookup(&variant.chr, variant.start, variant.end);
            if records.is_empty() {
                continue;
            }

            let mut variant_result = IndexMap::new();
            for record in &records {
                if !allele_match::matches_allele(record, variant, REF_COL, ALT_COL) {
                    continue;
                }

                if let Some(ada) = record.get(ADA_HEADER_NAME) {
                    if !ada.is_empty() && ada != "." {
                        variant_result.insert(ADA_FIELD.to_string(), ada.to_string());
                    }
                }

                if let Some(rf) = record.get(RF_HEADER_NAME) {
                    if !rf.is_empty() && rf != "." {
                        variant_result.insert(RF_FIELD.to_string(), rf.to_string());
                    }
                }

                if !variant_result.is_empty() {
                    break;
                }
            }

            if variant_result.is_empty() {
                continue;
            }

            for tc in &mut variant.transcript_consequences {
                for (k, v) in &variant_result {
                    tc.plugin_data.insert(k.clone(), v.clone());
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
        _consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let store = self
            .store
            .as_ref()
            .ok_or_else(|| PluginError::Run("dbscSNV plugin not initialized".into()))?;

        let records = store.query(&variant.chr, variant.start, variant.end)?;

        for record in &records {
            if !allele_match::matches_allele(record, variant, REF_COL, ALT_COL) {
                continue;
            }

            let mut result = IndexMap::new();

            if let Some(ada) = record.get(ADA_HEADER_NAME) {
                if !ada.is_empty() && ada != "." {
                    result.insert(ADA_FIELD.to_string(), ada.to_string());
                }
            }

            if let Some(rf) = record.get(RF_HEADER_NAME) {
                if !rf.is_empty() && rf != "." {
                    result.insert(RF_FIELD.to_string(), rf.to_string());
                }
            }

            if !result.is_empty() {
                return Ok(result);
            }
        }

        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plugin_name_is_correct() {
        let plugin = DbscSnvPlugin::new();
        assert_eq!(plugin.name(), "dbscSNV");
    }

    #[test]
    fn header_info_has_expected_fields() {
        let plugin = DbscSnvPlugin::new();
        let headers = plugin.header_info();
        assert_eq!(headers.len(), 2);
        assert_eq!(headers[0].0, "dbscSNV_ADA_SCORE");
        assert_eq!(headers[1].0, "dbscSNV_RF_SCORE");
    }

    #[test]
    fn init_fails_without_file_path() {
        let mut plugin = DbscSnvPlugin::new();
        let result = plugin.init(&[]);
        assert!(result.is_err());
    }

    /// dbscSNV1.1 ships its header with **no** leading `#` on both assemblies
    /// (GRCh37 `chr pos ref alt ada_score rf_score`, GRCh38 `hg38_chr hg38_pos
    /// ...`) and uses CRLF line endings. A `#`-only header check leaves `header`
    /// `None`, so `TabixRecord::get("ada_score")` returns `None` for every record
    /// and the plugin annotates **nothing while exiting 0**. The fixture is a real
    /// bgzip+tabix file shaped like the upstream one.
    #[test]
    fn annotates_from_headerless_crlf_file_as_published() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "chr\tpos\tref\talt\tada_score\trf_score\r",
            "21\t10907039\tG\tA\t0.00106\t0.028\r",
            "21\t10907040\tC\tT\t0.91234\t0.874\r",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "dbscsnv_grch37.txt.gz",
            &lines,
            IndexSpec::standard().bare_header(),
        );

        let mut plugin = DbscSnvPlugin::new();
        plugin
            .init(&[path.to_string_lossy().to_string()])
            .expect("init should succeed against the published file layout");

        let mut variants = vec![InputVariant::new(
            "21".into(),
            10907039,
            10907039,
            b"G".to_vec(),
            b"A".to_vec(),
        )];
        variants[0]
            .transcript_consequences
            .push(TranscriptConsequence::default());

        plugin.run_batch(&mut variants).expect("run_batch");

        let annotated = &variants[0].transcript_consequences[0].plugin_data;
        assert_eq!(
            annotated.get(ADA_FIELD).map(String::as_str),
            Some("0.00106"),
            "ADA score must be read from a headerless file; got {annotated:?}"
        );
        // The trailing \r must not survive onto the final column's value.
        assert_eq!(
            annotated.get(RF_FIELD).map(String::as_str),
            Some("0.028"),
            "RF score must be CR-stripped; got {annotated:?}"
        );
    }

    /// The GRCh38 file uses different header names (`hg38_chr`, `hg38_pos`) for
    /// the same first two columns, and the same code path reads it.
    #[test]
    fn annotates_from_grch38_headerless_file() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "hg38_chr\thg38_pos\tref\talt\tada_score\trf_score\r",
            "21\t20998616\tA\tC\t0.00386\t0.026\r",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "dbscsnv_grch38.txt.gz",
            &lines,
            IndexSpec::standard().bare_header(),
        );

        let mut plugin = DbscSnvPlugin::new();
        plugin.init(&[path.to_string_lossy().to_string()]).unwrap();

        let mut variants = vec![InputVariant::new(
            "21".into(),
            20998616,
            20998616,
            b"A".to_vec(),
            b"C".to_vec(),
        )];
        variants[0]
            .transcript_consequences
            .push(TranscriptConsequence::default());

        plugin.run_batch(&mut variants).unwrap();

        let annotated = &variants[0].transcript_consequences[0].plugin_data;
        assert_eq!(
            annotated.get(ADA_FIELD).map(String::as_str),
            Some("0.00386")
        );
        assert_eq!(annotated.get(RF_FIELD).map(String::as_str), Some("0.026"));
    }

    /// A GRCh38 run pointed at the GRCh37 file (or the reverse) must be refused:
    /// without the guard both directions annotate nothing at exit 0, which is
    /// indistinguishable from an input with no splice-region variants.
    #[test]
    fn rejects_the_other_assemblys_file() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();

        let g37 = write_tabix_fixture(
            dir.path(),
            "g37.txt.gz",
            &[
                "chr\tpos\tref\talt\tada_score\trf_score",
                "21\t10907039\tG\tA\t0.00106\t0.028",
            ],
            IndexSpec::standard().bare_header(),
        );
        let g38 = write_tabix_fixture(
            dir.path(),
            "g38.txt.gz",
            &[
                "hg38_chr\thg38_pos\tref\talt\tada_score\trf_score",
                "21\t10413724\tG\tA\t0.00002\t0.004",
            ],
            IndexSpec::standard().bare_header(),
        );

        let mut plugin = DbscSnvPlugin::new();
        plugin.set_assembly(Some(vep_core::assembly::Assembly::Grch38));
        let err = plugin
            .init(&[g37.to_string_lossy().to_string()])
            .expect_err("a GRCh37 file must be refused for a GRCh38 run");
        assert!(
            err.to_string().contains("GRCh37 file"),
            "error must name the file's actual assembly, got: {err}"
        );

        let mut plugin = DbscSnvPlugin::new();
        plugin.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        let err = plugin
            .init(&[g38.to_string_lossy().to_string()])
            .expect_err("a GRCh38 file must be refused for a GRCh37 run");
        assert!(
            err.to_string().contains("GRCh38 file"),
            "error must name the file's actual assembly, got: {err}"
        );

        let mut plugin = DbscSnvPlugin::new();
        plugin.set_assembly(Some(vep_core::assembly::Assembly::Grch37));
        plugin.init(&[g37.to_string_lossy().to_string()]).unwrap();
        let mut plugin = DbscSnvPlugin::new();
        plugin.set_assembly(Some(vep_core::assembly::Assembly::Grch38));
        plugin.init(&[g38.to_string_lossy().to_string()]).unwrap();

        let mut plugin = DbscSnvPlugin::new();
        plugin.init(&[g37.to_string_lossy().to_string()]).unwrap();
    }
}
