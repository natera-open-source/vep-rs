// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! AlphaMissense plugin: tabix-based pathogenicity scores for missense variants.
//!
//! Queries a tabix-indexed AlphaMissense file and annotates variants with
//! `am_class` (likely_pathogenic/likely_benign/ambiguous) and `am_pathogenicity` (0.0-1.0).
//! Uses allele matching to ensure the correct alt allele is matched.
//!
//! This is a variant-level plugin that writes directly to `variant.plugin_data`
//! via [`run_batch`], similar to CADD.

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{open_best_store, AnnotationStore};
use crate::tabix::{allele_match, BatchQueryResult, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// Column indices in the AlphaMissense TSV.
const CHR_COL: usize = 0;
const POS_COL: usize = 1;
const REF_COL: usize = 2;
const ALT_COL: usize = 3;
/// Every data row names its own build here, in UCSC spelling (`hg19` / `hg38`).
const GENOME_COL: usize = 4;
const PATHOGENICITY_COL: usize = 8;
const CLASS_COL: usize = 9;

/// AlphaMissense pathogenicity annotation plugin.
pub struct AlphaMissensePlugin {
    /// Annotation store, initialized in `init()`.
    store: Option<Box<dyn AnnotationStore>>,
    /// Target assembly from `--assembly`, injected before `init()`.
    assembly: Option<vep_core::assembly::Assembly>,
}

impl AlphaMissensePlugin {
    /// Create a new uninitialized AlphaMissense plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            assembly: None,
        }
    }
}

impl Default for AlphaMissensePlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl BuiltinPlugin for AlphaMissensePlugin {
    fn name(&self) -> &str {
        "AlphaMissense"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            (
                "am_class".to_string(),
                "AlphaMissense classification (likely_pathogenic, likely_benign, ambiguous)"
                    .to_string(),
            ),
            (
                "am_pathogenicity".to_string(),
                "AlphaMissense pathogenicity score (0.0-1.0)".to_string(),
            ),
        ]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let mut file_path: Option<&str> = None;
        for param in params {
            if let Some(path) = param.strip_prefix("file=") {
                file_path = Some(path);
                break;
            }
            if !param.contains('=') && file_path.is_none() {
                file_path = Some(param.as_str());
            }
        }
        let file_path = file_path.ok_or_else(|| {
            PluginError::Init("AlphaMissense requires a file path argument".into())
        })?;

        let config = TabixConfig {
            file_path: file_path.into(),
            chr_col: CHR_COL,
            start_col: POS_COL,
            end_col: None,
            ref_col: Some(REF_COL),
            alt_col: Some(ALT_COL),
            zero_based: false,
        };

        // The two upstream files are pre-split per build, identical in shape and
        // both UCSC-contig, so a wrong-build file resolves contigs and yields zero
        // hits at exit 0. The `genome` column is the only in-file build signal.
        if let Some(requested) = self.assembly {
            let declared =
                crate::tabix::first_data_row_field(std::path::Path::new(file_path), GENOME_COL)?;
            if let Some(declared) = declared {
                let declared = declared.trim();
                if !declared.is_empty()
                    && !declared.eq_ignore_ascii_case(requested.ucsc_name())
                    && !declared.eq_ignore_ascii_case(requested.as_str())
                {
                    return Err(PluginError::Init(format!(
                        "AlphaMissense: --assembly is {requested} but '{file_path}' declares \
                         genome '{declared}'. AlphaMissense publishes one file per build; \
                         point --plugin AlphaMissense at the {requested} file."
                    )));
                }
            }
        }

        self.store = Some(open_best_store(std::path::Path::new(file_path), config)?);
        debug!("AlphaMissense plugin initialized with {}", file_path);
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

        let batch_result: BatchQueryResult = store.query_batch(&regions)?;
        Ok(Some(PrefetchData::new(batch_result)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch_result: BatchQueryResult = data.downcast()?;

        if batch_result.is_empty() {
            return Ok(());
        }

        for variant in variants.iter_mut() {
            let records = batch_result.lookup(&variant.chr, variant.start, variant.end);

            for record in records {
                if !allele_match::matches_allele(record, variant, REF_COL, ALT_COL) {
                    continue;
                }

                if let Some(pathogenicity) = record.columns.get(PATHOGENICITY_COL) {
                    if !pathogenicity.is_empty() && pathogenicity != "." {
                        variant
                            .plugin_data
                            .insert("am_pathogenicity".to_string(), pathogenicity.clone());
                    }
                }

                if let Some(class) = record.columns.get(CLASS_COL) {
                    if !class.is_empty() && class != "." {
                        variant
                            .plugin_data
                            .insert("am_class".to_string(), class.clone());
                    }
                }

                break;
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
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // Variant-level plugin: annotation happens in `run_batch`.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plugin_name_is_correct() {
        let plugin = AlphaMissensePlugin::new();
        assert_eq!(plugin.name(), "AlphaMissense");
    }

    #[test]
    fn header_info_has_expected_fields() {
        let plugin = AlphaMissensePlugin::new();
        let headers = plugin.header_info();
        assert_eq!(headers.len(), 2);
        assert_eq!(headers[0].0, "am_class");
        assert_eq!(headers[1].0, "am_pathogenicity");
    }

    #[test]
    fn init_fails_without_file_path() {
        let mut plugin = AlphaMissensePlugin::new();
        let result = plugin.init(&[]);
        assert!(result.is_err());
    }

    #[test]
    fn init_accepts_file_prefix() {
        let mut plugin = AlphaMissensePlugin::new();
        let result = plugin.init(&["file=/path/to/alpha.tsv.gz".to_string()]);
        assert!(result.is_err());
    }

    /// Both AlphaMissense builds use UCSC contigs and identical column layouts,
    /// so a wrong-build file resolves contigs and then annotates nothing at exit
    /// 0. The `genome` column is the only in-file signal, so the guard must read
    /// it and refuse, and must still accept the matching build.
    #[test]
    fn refuses_the_other_builds_file_and_accepts_its_own() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        use vep_core::assembly::Assembly;

        let dir = tempfile::tempdir().unwrap();
        let hg38_rows = [
            "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class",
            "chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST1\tA1B\t0.9123\tlikely_pathogenic",
        ];
        let hg19_rows = [
            "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class",
            "chr21\t9907250\tC\tA\thg19\tQ9NQ38\tENST1\tA1B\t0.9123\tlikely_pathogenic",
        ];

        let hg38 = write_tabix_fixture(
            dir.path(),
            "am_hg38.tsv.gz",
            &hg38_rows,
            IndexSpec::standard(),
        );
        let hg19 = write_tabix_fixture(
            dir.path(),
            "am_hg19.tsv.gz",
            &hg19_rows,
            IndexSpec::standard(),
        );

        let mut plugin = AlphaMissensePlugin::new();
        plugin.set_assembly(Some(Assembly::Grch37));
        let err = plugin
            .init(&[hg38.to_string_lossy().to_string()])
            .expect_err("a GRCh37 run must refuse an hg38 AlphaMissense file");
        assert!(format!("{err}").contains("hg38"), "got: {err}");

        let mut plugin = AlphaMissensePlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        let err = plugin
            .init(&[hg19.to_string_lossy().to_string()])
            .expect_err("a GRCh38 run must refuse an hg19 AlphaMissense file");
        assert!(format!("{err}").contains("hg19"), "got: {err}");

        let mut plugin = AlphaMissensePlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        plugin
            .init(&[hg38.to_string_lossy().to_string()])
            .expect("a GRCh38 run must accept an hg38 AlphaMissense file");

        let mut plugin = AlphaMissensePlugin::new();
        plugin
            .init(&[hg19.to_string_lossy().to_string()])
            .expect("no --assembly must not block loading");
    }
}
