// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! CADD plugin: Combined Annotation Dependent Depletion scores.
//!
//! Tabix plugin with allele matching. Annotates variants with CADD
//! pathogenicity scores (PHRED-scaled and raw) from tabix-indexed files.
//! Supports separate SNV and indel score files.
//! Results are stored at the variant level (plugin_data on InputVariant).
//!
//! `CADD.pm` is the plugin's Perl counterpart in Ensembl/VEP_plugins release/115.

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{open_best_store, AnnotationStore};
use crate::tabix::allele_match::matches_allele_with_position;
use crate::tabix::{BatchQueryResult, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::assembly::Assembly;
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// A configured annotation source with column indices resolved from its header.
struct CaddSource {
    store: Box<dyn AnnotationStore>,
    /// 0-based column index for the reference allele.
    ref_col: usize,
    /// 0-based column index for the alternate allele.
    alt_col: usize,
    /// 0-based column index for the raw CADD score.
    raw_score_col: usize,
    /// 0-based column index for the PHRED-scaled CADD score.
    phred_col: usize,
}

/// CADD pathogenicity score plugin.
///
/// Annotates variants with CADD PHRED and raw scores by querying one or
/// more tabix-indexed score files. Allele matching is performed to ensure
/// the annotation corresponds to the correct ref/alt pair.
pub struct CaddPlugin {
    /// Tabix sources (SNV, indels, or generic), initialized after `init()`.
    sources: Vec<CaddSource>,
    /// Whether to force annotation even if existing annotations are present.
    force_annotate: bool,
    /// Target assembly from `--assembly`, injected before `init()`; a data file
    /// whose own provenance line declares the other build is rejected.
    assembly: Option<Assembly>,
}

impl Default for CaddPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl CaddPlugin {
    /// Create an uninitialized CADD plugin.
    pub fn new() -> Self {
        Self {
            sources: Vec::new(),
            force_annotate: false,
            assembly: None,
        }
    }

    /// Reject a data file whose own assembly marker contradicts `--assembly`.
    ///
    /// CADD's first line is a provenance comment naming the build
    /// (`##CADD GRCh38-v1.7 (c) University of Washington...`). Perl VEP's
    /// `CADD.pm` greps for the requested assembly and dies on a mismatch; without
    /// the guard a GRCh37 file in a GRCh38 run yields nothing or wrong-genome
    /// scores. A file with no marker is accepted: absence of evidence is not a
    /// mismatch.
    fn check_assembly(path: &str, requested: Option<Assembly>) -> Result<(), PluginError> {
        let Some(requested) = requested else {
            return Ok(());
        };
        let found = crate::tabix::detect_assembly_marker(std::path::Path::new(path))?;
        match found {
            Some(found) if found != requested => Err(PluginError::Init(format!(
                "CADD: --assembly is {requested} but the data file '{path}' declares {found} \
                 in its header. Supply the {requested} CADD file."
            ))),
            _ => Ok(()),
        }
    }

    /// Open an annotation source and resolve CADD column indices from its header.
    fn open_source(path: &str) -> Result<CaddSource, PluginError> {
        let config = TabixConfig {
            file_path: path.into(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: Some(2),
            alt_col: Some(3),
            zero_based: false,
        };

        let store = open_best_store(std::path::Path::new(path), config)?;

        // Column names and positions vary across CADD releases, so indices are
        // resolved from the header.
        let header = store
            .header()
            .ok_or_else(|| PluginError::Init(format!("CADD file {path} has no header")))?;

        let ref_col = find_column(header, &["Ref", "REF", "ref"])
            .ok_or_else(|| PluginError::Init(format!("CADD file {path}: no Ref column found")))?;
        let alt_col = find_column(header, &["Alt", "ALT", "alt", "Alt_allele"])
            .ok_or_else(|| PluginError::Init(format!("CADD file {path}: no Alt column found")))?;
        let raw_score_col = find_column(header, &["RawScore", "CADD_RAW", "raw", "RawScore-v1.7"])
            .ok_or_else(|| {
                PluginError::Init(format!("CADD file {path}: no RawScore column found"))
            })?;
        let phred_col = find_column(header, &["PHRED", "CADD_PHRED", "phred", "PHRED-v1.7"])
            .ok_or_else(|| PluginError::Init(format!("CADD file {path}: no PHRED column found")))?;

        debug!(
            file = path,
            ref_col, alt_col, raw_score_col, phred_col, "opened CADD source"
        );

        Ok(CaddSource {
            store,
            ref_col,
            alt_col,
            raw_score_col,
            phred_col,
        })
    }
}

/// Find a column index by trying multiple possible header names.
fn find_column(header: &[String], candidates: &[&str]) -> Option<usize> {
    for candidate in candidates {
        if let Some(idx) = header.iter().position(|h| h == candidate) {
            return Some(idx);
        }
    }
    None
}

impl BuiltinPlugin for CaddPlugin {
    fn name(&self) -> &str {
        "CADD"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            (
                "CADD_PHRED".to_string(),
                "CADD PHRED-scaled pathogenicity score".to_string(),
            ),
            (
                "CADD_RAW".to_string(),
                "CADD raw pathogenicity score".to_string(),
            ),
        ]
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        if params.is_empty() {
            return Err(PluginError::Init(
                "CADD requires at least one file path parameter".into(),
            ));
        }

        let mut snv_paths: Vec<String> = Vec::new();
        let mut indel_paths: Vec<String> = Vec::new();
        let mut generic_paths: Vec<String> = Vec::new();

        for param in params {
            if let Some(path) = param.strip_prefix("snv=") {
                snv_paths.push(path.to_string());
            } else if let Some(path) = param.strip_prefix("indels=") {
                indel_paths.push(path.to_string());
            } else if let Some(path) = param.strip_prefix("file=") {
                generic_paths.push(path.to_string());
            } else if let Some(val) = param.strip_prefix("force_annotate=") {
                self.force_annotate = val == "1" || val == "true";
            } else {
                generic_paths.push(param.clone());
            }
        }

        for path in snv_paths
            .iter()
            .chain(indel_paths.iter())
            .chain(generic_paths.iter())
        {
            Self::check_assembly(path, self.assembly)?;
            let source = Self::open_source(path)?;
            self.sources.push(source);
        }

        if self.sources.is_empty() {
            return Err(PluginError::Init(
                "CADD: no valid file paths provided".into(),
            ));
        }

        debug!(
            sources = self.sources.len(),
            force_annotate = self.force_annotate,
            "CADD plugin initialized"
        );

        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        if self.sources.is_empty() || variants.is_empty() {
            return Ok(None);
        }

        let mut batch_results = Vec::with_capacity(self.sources.len());
        for source in &self.sources {
            let regions: Vec<(&str, u64, u64)> = variants
                .iter()
                .map(|v| {
                    let query_start = v.start.saturating_sub(2);
                    (v.chr.as_str(), query_start, v.end)
                })
                .collect();

            batch_results.push(source.store.query_batch(&regions)?);
        }

        Ok(Some(PrefetchData::new(batch_results)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch_results: Vec<BatchQueryResult> = data.downcast()?;

        for (source, batch_result) in self.sources.iter().zip(batch_results.iter()) {
            if batch_result.is_empty() {
                continue;
            }

            for variant in variants.iter_mut() {
                if !self.force_annotate && variant.plugin_data.contains_key("CADD_PHRED") {
                    continue;
                }

                let query_start = variant.start.saturating_sub(2);
                let records = batch_result.lookup(&variant.chr, query_start, variant.end);

                for record in &records {
                    if matches_allele_with_position(record, variant, source.ref_col, source.alt_col)
                    {
                        if let Some(phred) = record.columns.get(source.phred_col) {
                            if !phred.is_empty() && phred != "." {
                                variant
                                    .plugin_data
                                    .insert("CADD_PHRED".to_string(), phred.clone());
                            }
                        }

                        if let Some(raw) = record.columns.get(source.raw_score_col) {
                            if !raw.is_empty() && raw != "." {
                                variant
                                    .plugin_data
                                    .insert("CADD_RAW".to_string(), raw.clone());
                            }
                        }

                        break;
                    }
                }
            }
        }

        Ok(())
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        // Combined prefetch and annotate for callers that skip `prefetch`.
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
    fn new_plugin_is_uninitialized() {
        let plugin = CaddPlugin::new();
        assert_eq!(plugin.name(), "CADD");
        assert!(plugin.sources.is_empty());
        assert!(!plugin.force_annotate);
    }

    #[test]
    fn init_requires_params() {
        let mut plugin = CaddPlugin::new();
        let result = plugin.init(&[]);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("requires at least one"));
    }

    #[test]
    fn header_info_returns_both_fields() {
        let plugin = CaddPlugin::new();
        let headers = plugin.header_info();
        assert_eq!(headers.len(), 2);
        let names: Vec<&str> = headers.iter().map(|(k, _)| k.as_str()).collect();
        assert!(names.contains(&"CADD_PHRED"));
        assert!(names.contains(&"CADD_RAW"));
    }

    #[test]
    fn run_batch_skips_when_no_sources() {
        let plugin = CaddPlugin::new();
        let mut variants = vec![];
        assert!(plugin.run_batch(&mut variants).is_ok());
    }

    #[test]
    fn find_column_resolves_candidates() {
        let header: Vec<String> = vec![
            "Chrom".into(),
            "Pos".into(),
            "Ref".into(),
            "Alt".into(),
            "RawScore".into(),
            "PHRED".into(),
        ];
        assert_eq!(find_column(&header, &["Ref", "REF"]), Some(2));
        assert_eq!(find_column(&header, &["PHRED", "phred"]), Some(5));
        assert_eq!(find_column(&header, &["missing"]), None);
    }

    #[test]
    fn parse_named_params() {
        let mut plugin = CaddPlugin::new();
        // The files do not exist, so `init` fails after the parameters are parsed.
        let result = plugin.init(&[
            "snv=/path/to/snv.tsv.gz".into(),
            "indels=/path/to/indels.tsv.gz".into(),
            "force_annotate=1".into(),
        ]);
        assert!(result.is_err());
        assert!(plugin.force_annotate);
    }

    #[test]
    fn parse_file_param() {
        let mut plugin = CaddPlugin::new();
        let result = plugin.init(&["file=/path/to/snv.tsv.gz".into()]);
        assert!(result.is_err());
    }

    /// Build a CADD fixture with the upstream two-line preamble: a provenance
    /// comment naming the build, then the column-name row.
    fn cadd_fixture(dir: &std::path::Path, name: &str, assembly_token: &str) -> std::path::PathBuf {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        let provenance =
            format!("##CADD {assembly_token}-v1.7 (c) University of Washington 2013-2020");
        let lines = vec![
            provenance.as_str(),
            "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED",
            "21\t5030000\tA\tC\t0.754016\t7.867",
            "21\t5030000\tA\tG\t0.780435\t8.106",
            "21\t5030001\tC\tA\t0.625243\t6.694",
        ];
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
    }

    /// A GRCh37 CADD file in a GRCh38 run must be refused: Perl VEP's CADD.pm
    /// greps the header for the assembly and dies.
    #[test]
    fn rejects_a_data_file_from_the_other_assembly() {
        let dir = tempfile::tempdir().unwrap();
        let g37 = cadd_fixture(dir.path(), "cadd_grch37.tsv.gz", "GRCh37");

        let mut plugin = CaddPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        let err = plugin
            .init(&[format!("snv={}", g37.display())])
            .expect_err("a GRCh37 file must be rejected in a GRCh38 run");

        let msg = err.to_string();
        assert!(
            msg.contains("GRCh38") && msg.contains("GRCh37"),
            "error must name both assemblies: {msg}"
        );
    }

    /// The matching case must load and annotate, so the guard is not simply
    /// rejecting everything.
    #[test]
    fn accepts_and_annotates_a_matching_assembly_file() {
        let dir = tempfile::tempdir().unwrap();
        let g38 = cadd_fixture(dir.path(), "cadd_grch38.tsv.gz", "GRCh38");

        let mut plugin = CaddPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        plugin
            .init(&[format!("snv={}", g38.display())])
            .expect("a matching-assembly file must load");

        let mut variants = vec![InputVariant::new(
            "21".into(),
            5030000,
            5030000,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        plugin.run_batch(&mut variants).expect("run_batch");

        // The A>G allele at 5030000 must match, not the A>C on the preceding
        // line, which also exercises the multi-record read path.
        assert_eq!(
            variants[0]
                .plugin_data
                .get("CADD_PHRED")
                .map(String::as_str),
            Some("8.106"),
            "expected the A>G score; got {:?}",
            variants[0].plugin_data
        );
        assert_eq!(
            variants[0].plugin_data.get("CADD_RAW").map(String::as_str),
            Some("0.780435")
        );
    }

    /// A file with no assembly marker is accepted: absence of evidence is not a
    /// mismatch, and several upstream files carry no provenance line.
    #[test]
    fn accepts_a_file_with_no_assembly_marker() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED",
            "21\t5030000\tA\tG\t0.780435\t8.106",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "cadd_nomarker.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );

        let mut plugin = CaddPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        plugin
            .init(&[format!("snv={}", path.display())])
            .expect("an unmarked file must not be rejected");
    }

    /// Without `--assembly` there is nothing to contradict, so any file loads.
    #[test]
    fn no_requested_assembly_skips_the_guard() {
        let dir = tempfile::tempdir().unwrap();
        let g37 = cadd_fixture(dir.path(), "cadd_grch37b.tsv.gz", "GRCh37");

        let mut plugin = CaddPlugin::new();
        plugin
            .init(&[format!("snv={}", g37.display())])
            .expect("no --assembly means no mismatch to detect");
    }
}
