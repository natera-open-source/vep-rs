// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! REVEL plugin: missense variant pathogenicity scores.
//!
//! Queries a tabix-indexed REVEL TSV file to annotate transcript consequences
//! with REVEL pathogenicity scores. Only missense variants are annotated.
//! Matches by position + allele, optionally by Ensembl transcript ID.
//!
//! # Parameters
//!
//! - `file=<path>` (or first positional): path to the tabix-indexed REVEL TSV
//! - `no_match=1`: report score even when transcript ID doesn't match
//!
//! # Assembly handling
//!
//! A single REVEL download serves **both** assemblies: its columns are
//! `chr, hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL,
//! Ensembl_transcriptid`, so GRCh37 coordinates live in column 2 and GRCh38
//! coordinates in column 3. The upstream preparation instructions index each
//! assembly on its own column (`tabix -s 1 -b 2 -e 2` for GRCh37,
//! `-s 1 -b 3 -e 3` for GRCh38 after re-sorting).
//!
//! This plugin therefore resolves the position column **by header name** for the
//! requested assembly (`hg19_pos` / `grch38_pos`), exactly as Perl VEP's
//! `REVEL.pm` does, and errors when the required column is absent. Reading a
//! fixed column instead would stamp every GRCh38 record with its hg19 position;
//! because the engine filters records on the re-parsed position, that discards
//! every record and annotates nothing while still exiting 0.
//!
//! # Output fields
//!
//! - `REVEL_score`: REVEL pathogenicity score (0 to 1, higher = more pathogenic)
//!
//! `REVEL.pm` is the plugin's Perl counterpart in Ensembl/VEP_plugins release/115.

use indexmap::IndexMap;
use vep_core::assembly::Assembly;
use vep_core::consequence::{Consequence, TranscriptConsequence};
use vep_core::variant::InputVariant;

use crate::annotation_store::{open_best_store, AnnotationStore};
use crate::tabix::{allele_match, BatchQueryResult, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};

/// Header name of the GRCh37 position column in a REVEL TSV.
const HG19_POS_COLUMN: &str = "hg19_pos";

/// Header name of the GRCh38 position column in a REVEL TSV.
const GRCH38_POS_COLUMN: &str = "grch38_pos";

/// REVEL pathogenicity score plugin for missense variants.
///
/// REVEL (Rare Exome Variant Ensemble Learner) integrates scores from 13
/// individual pathogenicity tools into a single meta-predictor. This plugin
/// reads a tabix-indexed TSV file containing pre-computed REVEL scores and
/// annotates each missense transcript consequence with its score.
pub struct RevelPlugin {
    /// Annotation store (tabix or binary).
    store: Option<Box<dyn AnnotationStore>>,
    /// If true, report scores even when transcript ID doesn't match.
    no_match: bool,
    /// 0-based column index for the reference allele.
    ref_col: usize,
    /// 0-based column index for the alternate allele.
    alt_col: usize,
    /// Target assembly from `--assembly`, injected before `init()`. Selects
    /// which of the file's two position columns tabix indexed.
    assembly: Option<Assembly>,
}

impl Default for RevelPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl RevelPlugin {
    /// Create an uninitialized REVEL plugin.
    ///
    /// Call [`init`](BuiltinPlugin::init) with file path before use.
    pub fn new() -> Self {
        Self {
            store: None,
            no_match: false,
            // Defaults for a `chr, pos, ref, alt` layout; `init()` overrides them
            // from the header.
            ref_col: 2,
            alt_col: 3,
            assembly: None,
        }
    }

    /// Header name of the position column for the requested assembly.
    fn position_column_for(assembly: Assembly) -> &'static str {
        match assembly {
            Assembly::Grch37 => HG19_POS_COLUMN,
            Assembly::Grch38 => GRCH38_POS_COLUMN,
        }
    }

    /// Resolve the 0-based position column index for `assembly` from `header`.
    ///
    /// Mirrors Perl `REVEL.pm`, which dies when the assembly's column is absent
    /// rather than silently reading the other assembly's coordinates.
    fn resolve_position_column(
        header: &[String],
        assembly: Assembly,
    ) -> Result<usize, PluginError> {
        let wanted = Self::position_column_for(assembly);
        header
            .iter()
            .position(|h| h.trim().eq_ignore_ascii_case(wanted))
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "REVEL: assembly is {assembly} but the data file has no '{wanted}' column \
                     (header: {}). REVEL ships one file carrying both assemblies; index it for \
                     this assembly and confirm the header retains its column names.",
                    header.join(",")
                ))
            })
    }

    /// Check whether a transcript consequence includes missense_variant.
    fn is_missense(consequence: &TranscriptConsequence) -> bool {
        consequence
            .consequences
            .contains(&Consequence::MissenseVariant)
    }

    fn transcript_matches(&self, record: &crate::tabix::TabixRecord, tc_id: &str) -> bool {
        if self.no_match {
            return true;
        }

        let Some(file_transcript) = record.get("Ensembl_transcriptid") else {
            return true;
        };

        let tc_base = tc_id.split('.').next().unwrap_or(tc_id);
        file_transcript.split(';').any(|t| {
            let t_trimmed = t.trim();
            let t_base = t_trimmed.split('.').next().unwrap_or(t_trimmed);
            t_base == tc_base
        })
    }

    fn extract_score(record: &crate::tabix::TabixRecord) -> Option<&str> {
        record
            .get("REVEL")
            .or_else(|| record.get("REVEL_score"))
            .filter(|score| !score.is_empty() && *score != ".")
    }
}

impl BuiltinPlugin for RevelPlugin {
    fn name(&self) -> &str {
        "REVEL"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            "REVEL_score".to_string(),
            "REVEL pathogenicity score for missense variants".to_string(),
        )]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let mut file_path = None;

        for param in params {
            if let Some(path) = param.strip_prefix("file=") {
                file_path = Some(path.to_string());
            } else if param == "no_match=1" || param == "no_match" {
                self.no_match = true;
            } else if file_path.is_none() {
                file_path = Some(param.to_string());
            }
        }

        let path = file_path
            .ok_or_else(|| PluginError::Init("REVEL plugin requires 'file' parameter".into()))?;

        let file_path_buf: std::path::PathBuf = path.into();

        // The position column is fixed in `TabixConfig` at construction and the
        // engine re-parses each record's position from it, so it is resolved
        // before the store opens. The REVEL header is `#`-prefixed, so the
        // `start_col` of 1 passed to the header peek cannot mis-detect it.
        let peeked = crate::tabix::parse_header(&file_path_buf, 1)?;
        let start_col = match (self.assembly, peeked.as_deref()) {
            (Some(assembly), Some(header)) => Self::resolve_position_column(header, assembly)?,
            // Without `--assembly` the column-2 (hg19) default applies and is
            // logged: a GRCh38-indexed file annotates nothing in that configuration.
            (None, _) => {
                tracing::warn!(
                    "REVEL: no --assembly given; assuming the GRCh37 position column \
                     ('{HG19_POS_COLUMN}', column 2). Pass --assembly GRCh38 when using a \
                     GRCh38-indexed REVEL file, or it will annotate nothing."
                );
                1
            }
            // No parseable header, so the column cannot be verified: fall back to
            // the index the upstream preparation instructions give the assembly.
            (Some(assembly), None) => {
                let fallback = match assembly {
                    Assembly::Grch37 => 1,
                    Assembly::Grch38 => 2,
                };
                tracing::warn!(
                    "REVEL: data file has no parseable header; assuming the documented \
                     {assembly} position column (column {}). Verify the file was indexed \
                     for {assembly}.",
                    fallback + 1
                );
                fallback
            }
        };

        // The two assembly copies of REVEL are byte-identical; only the column
        // recorded in the `.tbi` distinguishes them. An index built on the other
        // column stamps every record with the wrong position, and the engine's
        // position filter then discards them all: a silent zero at exit 0.
        if let (Some(assembly), Some(indexed)) = (
            self.assembly,
            crate::tabix::indexed_start_column(&file_path_buf)?,
        ) {
            if indexed != start_col {
                let indexed_name = peeked
                    .as_deref()
                    .and_then(|h| h.get(indexed))
                    .map(|s| s.trim().to_string())
                    .unwrap_or_else(|| format!("column {}", indexed + 1));
                return Err(PluginError::Init(format!(
                    "REVEL: assembly is {assembly}, which needs the '{}' column (column {}), \
                     but '{}' is indexed on '{indexed_name}' (column {}). REVEL ships ONE file \
                     carrying both assemblies' coordinates, so the two prepared copies are \
                     byte-identical and differ only in their tabix index. Point --plugin REVEL \
                     at the copy indexed for {assembly}.",
                    Self::position_column_for(assembly),
                    start_col + 1,
                    file_path_buf.display(),
                    indexed + 1,
                )));
            }
        }

        let config = TabixConfig {
            file_path: file_path_buf.clone(),
            chr_col: 0,
            start_col,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        };

        let store = open_best_store(&file_path_buf, config)?;

        if let Some(header) = store.header() {
            if let Some(idx) = header.iter().position(|h| h.eq_ignore_ascii_case("ref")) {
                self.ref_col = idx;
            }
            if let Some(idx) = header.iter().position(|h| h.eq_ignore_ascii_case("alt")) {
                self.alt_col = idx;
            }
        }

        self.store = Some(store);
        Ok(())
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
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

            for tc_idx in 0..variant.transcript_consequences.len() {
                if !Self::is_missense(&variant.transcript_consequences[tc_idx]) {
                    continue;
                }

                let transcript_id = variant.transcript_consequences[tc_idx]
                    .transcript_id
                    .clone();
                let mut score: Option<String> = None;
                for record in &records {
                    if !allele_match::matches_allele(record, variant, self.ref_col, self.alt_col) {
                        continue;
                    }
                    if !self.transcript_matches(record, &transcript_id) {
                        continue;
                    }
                    if let Some(value) = Self::extract_score(record) {
                        score = Some(value.to_string());
                        break;
                    }
                }

                if let Some(value) = score {
                    variant.transcript_consequences[tc_idx]
                        .plugin_data
                        .insert("REVEL_score".to_string(), value);
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
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        if !Self::is_missense(consequence) {
            return Ok(IndexMap::new());
        }

        let store = match &self.store {
            Some(s) => s,
            None => return Ok(IndexMap::new()),
        };

        let records = store.query(&variant.chr, variant.start, variant.end)?;

        for record in &records {
            if !allele_match::matches_allele(record, variant, self.ref_col, self.alt_col) {
                continue;
            }

            if !self.transcript_matches(record, &consequence.transcript_id) {
                continue;
            }

            if let Some(score) = Self::extract_score(record) {
                let mut result = IndexMap::new();
                result.insert("REVEL_score".to_string(), score.to_string());
                return Ok(result);
            }
        }

        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;

    #[test]
    fn new_plugin_is_uninitialized() {
        let plugin = RevelPlugin::new();
        assert_eq!(plugin.name(), "REVEL");
        assert!(plugin.store.is_none());
        assert!(!plugin.no_match);
    }

    #[test]
    fn header_info_returns_revel_score() {
        let plugin = RevelPlugin::new();
        let headers = plugin.header_info();
        assert_eq!(headers.len(), 1);
        assert_eq!(headers[0].0, "REVEL_score");
    }

    #[test]
    fn run_returns_empty_when_not_initialized() {
        let plugin = RevelPlugin::new();
        let consequence = TranscriptConsequence::default();
        let variant = InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let result = plugin.run(&consequence, &variant).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn run_skips_non_missense() {
        let mut plugin = RevelPlugin::new();
        plugin.store = None; // would fail if it tried to query

        let consequence = TranscriptConsequence {
            consequences: smallvec![Consequence::SynonymousVariant],
            ..Default::default()
        };

        let variant = InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let result = plugin.run(&consequence, &variant).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn is_missense_detects_correctly() {
        let mut tc = TranscriptConsequence {
            consequences: smallvec![Consequence::SynonymousVariant],
            ..Default::default()
        };
        assert!(!RevelPlugin::is_missense(&tc));

        tc.consequences = smallvec![Consequence::MissenseVariant];
        assert!(RevelPlugin::is_missense(&tc));

        tc.consequences = smallvec![
            Consequence::SpliceRegionVariant,
            Consequence::MissenseVariant,
        ];
        assert!(RevelPlugin::is_missense(&tc));
    }

    /// A REVEL row as published: one file, two coordinate columns.
    /// `chr, hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL, Ensembl_transcriptid`
    const REVEL_HEADER: &str =
        "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid";

    fn missense_variant(chr: &str, pos: u64, r: &[u8], a: &[u8]) -> InputVariant {
        let mut v = InputVariant::new(chr.into(), pos, pos, r.to_vec(), a.to_vec());
        let tc = TranscriptConsequence {
            consequences: smallvec![Consequence::MissenseVariant],
            ..Default::default()
        };
        v.transcript_consequences.push(tc);
        v
    }

    /// On GRCh38 the file is indexed on `grch38_pos` (column 3). The engine
    /// re-parses each record's position from the configured column and filters on
    /// it, so a hardcoded hg19 column stamps every record with its hg19 coordinate
    /// and discards it, annotating nothing while exiting 0. The fixture is indexed
    /// on column 3, mirroring the upstream GRCh38 preparation
    /// (`tabix -s 1 -b 3 -e 3`).
    #[test]
    fn grch38_reads_the_grch38_position_column() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        // hg19 pos 100200, GRCh38 pos 500600, deliberately far apart so
        // reading the wrong column cannot accidentally match.
        let lines = [
            REVEL_HEADER,
            "21\t100200\t500600\tA\tG\tM\tV\t0.742\tENST00000400000",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "revel_grch38.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(3),
        );

        let mut plugin = RevelPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        plugin
            .init(&[format!("file={}", path.display()), "no_match=1".to_string()])
            .expect("init");

        let mut variants = vec![missense_variant("21", 500600, b"A", b"G")];
        plugin.run_batch(&mut variants).expect("run_batch");

        assert_eq!(
            variants[0].transcript_consequences[0]
                .plugin_data
                .get("REVEL_score")
                .map(String::as_str),
            Some("0.742"),
            "GRCh38 run must read grch38_pos (column 3); got {:?}",
            variants[0].transcript_consequences[0].plugin_data
        );
    }

    /// The GRCh37 path reads `hg19_pos` (column 2).
    #[test]
    fn grch37_reads_the_hg19_position_column() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = [
            REVEL_HEADER,
            "21\t100200\t500600\tA\tG\tM\tV\t0.742\tENST00000400000",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "revel_grch37.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );

        let mut plugin = RevelPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch37));
        plugin
            .init(&[format!("file={}", path.display()), "no_match=1".to_string()])
            .expect("init");

        let mut variants = vec![missense_variant("21", 100200, b"A", b"G")];
        plugin.run_batch(&mut variants).expect("run_batch");

        assert_eq!(
            variants[0].transcript_consequences[0]
                .plugin_data
                .get("REVEL_score")
                .map(String::as_str),
            Some("0.742"),
            "GRCh37 run must read hg19_pos (column 2)"
        );
    }

    /// Perl's `REVEL.pm` dies when the requested assembly's column is missing
    /// rather than silently reading the other assembly's coordinates. Match that.
    #[test]
    fn missing_assembly_column_is_a_hard_error() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        // A GRCh37-only seven-column REVEL file with no `grch38_pos`.
        let lines = [
            "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL",
            "21\t100200\tA\tG\tM\tV\t0.742",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "revel_hg19_only.tsv.gz",
            &lines,
            IndexSpec::standard(),
        );

        let mut plugin = RevelPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        let err = plugin
            .init(&[format!("file={}", path.display())])
            .expect_err("GRCh38 against an hg19-only file must fail loudly");

        let msg = err.to_string();
        assert!(
            msg.contains("grch38_pos"),
            "error must name the missing column: {msg}"
        );
    }

    #[test]
    fn position_column_names_map_to_assemblies() {
        assert_eq!(
            RevelPlugin::position_column_for(Assembly::Grch37),
            "hg19_pos"
        );
        assert_eq!(
            RevelPlugin::position_column_for(Assembly::Grch38),
            "grch38_pos"
        );
    }

    /// REVEL ships one table carrying both `hg19_pos` (col 2) and `grch38_pos`
    /// (col 3), so the two prepared copies are byte-identical and differ only in
    /// which column their `.tbi` indexes. Pointing a GRCh38 run at the
    /// GRCh37-indexed copy passes every header check and annotates nothing at exit
    /// 0; the index's recorded start column is the only signal that catches it.
    #[test]
    fn rejects_a_copy_indexed_for_the_other_assembly() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};

        let dir = tempfile::tempdir().unwrap();
        let lines = [
            "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid",
            "21\t9907250\t9068417\tC\tA\tQ\tH\t0.193\tENST00000400754",
        ];

        let g37 = write_tabix_fixture(
            dir.path(),
            "revel_g37.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(2),
        );
        let g38 = write_tabix_fixture(
            dir.path(),
            "revel_g38.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(3),
        );

        let mut plugin = RevelPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        let err = plugin
            .init(&[g37.to_string_lossy().to_string()])
            .expect_err("a GRCh37-indexed copy must be refused for a GRCh38 run");
        let msg = err.to_string();
        assert!(
            msg.contains("grch38_pos") && msg.contains("hg19_pos"),
            "error must name both the needed and the indexed column, got: {msg}"
        );

        let mut plugin = RevelPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch37));
        let err = plugin
            .init(&[g38.to_string_lossy().to_string()])
            .expect_err("a GRCh38-indexed copy must be refused for a GRCh37 run");
        assert!(
            err.to_string().contains("grch38_pos"),
            "error must name the indexed column, got: {err}"
        );

        let mut plugin = RevelPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch37));
        plugin.init(&[g37.to_string_lossy().to_string()]).unwrap();
        let mut plugin = RevelPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        plugin.init(&[g38.to_string_lossy().to_string()]).unwrap();
    }
}
