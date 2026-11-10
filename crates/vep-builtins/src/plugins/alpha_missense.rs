// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! AlphaMissense plugin: tabix-based pathogenicity scores for missense variants.
//!
//! Queries a tabix-indexed AlphaMissense file and annotates variants with
//! `am_class` (likely_pathogenic/likely_benign/ambiguous) and `am_pathogenicity` (0.0-1.0).
//! Matches genomic alleles and each missense transcript's reference and alternate
//! amino acids, writing the first matching data row's scores on that transcript.

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{open_best_store, AnnotationStore};
use crate::tabix::{allele_match, BatchQueryResult, TabixConfig};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::consequence::{Consequence, TranscriptConsequence};
use vep_core::variant::InputVariant;

/// Column indices in the AlphaMissense TSV.
const CHR_COL: usize = 0;
const POS_COL: usize = 1;
const REF_COL: usize = 2;
const ALT_COL: usize = 3;
/// Every data row names its own build here, in UCSC spelling (`hg19` / `hg38`).
const GENOME_COL: usize = 4;
const PROTEIN_VARIANT_COL: usize = 7;
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
            let records: Vec<_> = batch_result
                .lookup(&variant.chr, variant.start, variant.end)
                .into_iter()
                .filter(|record| allele_match::matches_allele(record, variant, REF_COL, ALT_COL))
                .collect();

            for tc in &mut variant.transcript_consequences {
                // AlphaMissense.pm:242 requires missense among the transcript's terms.
                if !crate::traits::annotates_row(tc)
                    || !tc.consequences.contains(&Consequence::MissenseVariant)
                {
                    continue;
                }
                // AlphaMissense.pm:224-235 compares whole, transcript-oriented peptides
                // with the first and last protein-variant letters, ignoring their positions.
                let Some([aa_ref, b'/', aa_alt]) = tc.amino_acids.as_deref().map(str::as_bytes)
                else {
                    continue;
                };
                for record in &records {
                    let Some(protein_variant) = record.columns.get(PROTEIN_VARIANT_COL) else {
                        continue;
                    };
                    if protein_variant.as_bytes().first() != Some(aa_ref)
                        || protein_variant.as_bytes().last() != Some(aa_alt)
                    {
                        continue;
                    }

                    for (key, column) in [
                        ("am_pathogenicity", PATHOGENICITY_COL),
                        ("am_class", CLASS_COL),
                    ] {
                        if let Some(value) = record.columns.get(column) {
                            if !value.is_empty() && value != "." {
                                tc.plugin_data.insert(key.into(), value.clone());
                            }
                        }
                    }

                    break;
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
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // Variant-level plugin: annotation happens in `run_batch`.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use smallvec::smallvec;
    use vep_core::consequence::Consequence;
    use vep_core::coordinate::Strand;

    const HEADER: &str = "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class";

    fn missense(transcript_id: &str, amino_acids: &str) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript_id: transcript_id.into(),
            consequences: smallvec![Consequence::MissenseVariant],
            amino_acids: Some(amino_acids.into()),
            ..TranscriptConsequence::default()
        }
    }

    fn prediction(tc: &TranscriptConsequence) -> (Option<&str>, Option<&str>) {
        (
            tc.plugin_data.get("am_pathogenicity").map(String::as_str),
            tc.plugin_data.get("am_class").map(String::as_str),
        )
    }

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

    /// AlphaMissense.pm:195-197 limits annotation to transcripts, and
    /// AlphaMissense.pm:242,253 requires missense with matching peptides.
    #[test]
    fn writes_on_transcript_rows_alone() {
        use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
        use vep_core::consequence::FeatureType;

        let dir = tempfile::tempdir().unwrap();
        let rows = [
            "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class",
            "chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST1\tA1D\t0.9123\tlikely_pathogenic",
        ];
        let path = write_tabix_fixture(dir.path(), "am.tsv.gz", &rows, IndexSpec::standard());

        let mut plugin = AlphaMissensePlugin::new();
        plugin.init(&[path.to_string_lossy().to_string()]).unwrap();

        let mut variants = vec![InputVariant::new(
            "21".into(),
            9068417,
            9068417,
            b"C".to_vec(),
            b"A".to_vec(),
        )];
        for feature_type in [
            FeatureType::Transcript,
            FeatureType::RegulatoryFeature,
            FeatureType::MotifFeature,
            FeatureType::Transcript,
        ] {
            variants[0]
                .transcript_consequences
                .push(TranscriptConsequence {
                    feature_type,
                    ..missense("ENST1", "A/D")
                });
        }

        plugin.run_batch(&mut variants).unwrap();

        let classes: Vec<Option<&str>> = variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.get("am_class").map(String::as_str))
            .collect();
        assert_eq!(
            classes,
            [
                Some("likely_pathogenic"),
                None,
                None,
                Some("likely_pathogenic")
            ]
        );
        assert_eq!(
            variants[0].transcript_consequences[0]
                .plugin_data
                .get("am_pathogenicity")
                .map(String::as_str),
            Some("0.9123")
        );
        assert!(variants[0].plugin_data.is_empty());
    }

    /// AlphaMissense.pm:242-287 selects the first allele-and-peptide match
    /// separately for each missense transcript; AlphaMissense.pm:232-233
    /// leaves protein positions and transcript identifiers unrestricted by default.
    #[test]
    fn selects_matching_records_per_transcript_on_both_input_strands() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            HEADER,
            "chr21\t100\tA\tT\thg38\tP1\tENST1\tK11R\t0.01\tlikely_benign",
            "chr21\t100\tT\tG\thg38\tP1\tENST1\tK11R\t0.02\tlikely_benign",
            "chr21\t100\tA\tG\thg38\tP1\tENST1\tK11R\t0.99\tlikely_pathogenic",
            "chr21\t100\tA\tG\thg38\tP2\tENST2\tN99S\t0.45\tambiguous",
            "chr21\t100\tA\tG\thg38\tP3\tENST3\tF18S\t0.12\tlikely_benign",
            "chr21\t100\tA\tG\thg38\tP4\tENST4\tK77R\t0.22\tlikely_benign",
        ];
        let path = write_tabix_fixture(dir.path(), "am.tsv.gz", &rows, IndexSpec::standard());
        let mut plugin = AlphaMissensePlugin::new();
        plugin.init(&[path.to_string_lossy().into_owned()]).unwrap();

        for input_strand in [Strand::Forward, Strand::Reverse] {
            let (reference, alternate) = match input_strand {
                Strand::Reverse => (b"T", b"C"),
                _ => (b"A", b"G"),
            };
            let mut variant = InputVariant::new(
                "21".into(),
                100,
                100,
                reference.to_vec(),
                alternate.to_vec(),
            );
            variant.strand = input_strand;
            let mut first = missense("ENST_DIFFERENT", "K/R");
            first.protein_position = Some("700".into());
            first.strand = 1;
            first.consequences.push(Consequence::SpliceRegionVariant);
            let mut reverse = missense("ENST_REVERSE", "F/S");
            reverse.strand = -1;
            variant.transcript_consequences = vec![
                first,
                missense("ENST_OTHER", "N/S"),
                reverse,
                missense("ENST_BAD_REF", "A/R"),
                missense("ENST_BAD_ALT", "K/W"),
                TranscriptConsequence {
                    consequences: smallvec![Consequence::IntronVariant],
                    ..missense("ENST_INTRON", "K/R")
                },
                TranscriptConsequence {
                    consequences: smallvec![Consequence::SynonymousVariant],
                    ..missense("ENST_SYNONYMOUS", "K/R")
                },
                TranscriptConsequence {
                    amino_acids: None,
                    ..missense("ENST_NO_PEPTIDE", "K/R")
                },
            ];
            for absent_or_multiple in ["", ".", "-", "K", "/", "K/", "/R", "KR/R", "K/RS"] {
                variant
                    .transcript_consequences
                    .push(missense("ENST_INCOMPATIBLE", absent_or_multiple));
            }
            plugin
                .run_batch(std::slice::from_mut(&mut variant))
                .unwrap();
            let predictions: Vec<_> = variant
                .transcript_consequences
                .iter()
                .map(prediction)
                .collect();
            assert_eq!(
                &predictions[..3],
                &[
                    (Some("0.99"), Some("likely_pathogenic")),
                    (Some("0.45"), Some("ambiguous")),
                    (Some("0.12"), Some("likely_benign")),
                ],
                "{input_strand:?}"
            );
            assert!(
                predictions[3..].iter().all(|value| *value == (None, None)),
                "{input_strand:?}: {predictions:?}"
            );
            assert!(variant.plugin_data.is_empty());
        }
    }

    /// AlphaMissense.pm:224-235 compares both whole peptides with the first
    /// and last protein-variant letters; AlphaMissense.pm:280-282 returns the
    /// first matching row even when its selected values are absent.
    #[test]
    fn skips_absent_protein_changes_and_keeps_values_from_one_matching_row() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            HEADER,
            "chr21\t100\tA\tG\thg38\tP1\tENST1",
            "chr21\t100\tA\tG\thg38\tP1\tENST1\t.\t0.01\tlikely_benign",
            "chr21\t100\tA\tG\thg38\tP1\tENST1\t\t0.02\tlikely_benign",
            "chr21\t100\tA\tG\thg38\tP1\tENST1\tK11R\t.\t",
            "chr21\t100\tA\tG\thg38\tP1\tENST1\tK11R\t0.99\tlikely_pathogenic",
            "chr21\t100\tA\tG\thg38\tP2\tENST2\tN99S\t0.45\t.",
            "chr21\t100\tA\tG\thg38\tP3\tENST3\tF18S\t\tlikely_benign",
            "chr21\t100\tA\tG\thg38\tP4\tENST4\tD42G",
            "chr21\t101\tA\tG\thg38\tP1\tENST1\tK11R\t0.88\tlikely_pathogenic",
        ];
        let path = write_tabix_fixture(dir.path(), "am.tsv.gz", &rows, IndexSpec::standard());
        let mut plugin = AlphaMissensePlugin::new();
        plugin.init(&[path.to_string_lossy().into_owned()]).unwrap();
        let mut variants = [100, 101, 102].map(|position| {
            let mut variant = InputVariant::new(
                "21".into(),
                position,
                position,
                b"A".to_vec(),
                b"G".to_vec(),
            );
            variant.transcript_consequences = ["K/R", "N/S", "F/S", "D/G"]
                .map(|change| missense("ENST_OTHER", change))
                .into();
            variant
        });
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(
            variants[0]
                .transcript_consequences
                .iter()
                .map(prediction)
                .collect::<Vec<_>>(),
            [
                (None, None),
                (Some("0.45"), None),
                (None, Some("likely_benign")),
                (None, None)
            ]
        );
        assert_eq!(
            prediction(&variants[1].transcript_consequences[0]),
            (Some("0.88"), Some("likely_pathogenic"))
        );
        for tc in variants[1].transcript_consequences[1..]
            .iter()
            .chain(&variants[2].transcript_consequences)
        {
            assert_eq!(prediction(tc), (None, None));
        }
    }
}
