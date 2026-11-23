// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! AVI plugin: AlphaGenome Variant Impact scores for single-base substitutions.
//!
//! Reads a bgzipped, tabix-indexed six-column TSV (chromosome, position,
//! reference allele, alternate allele, raw score, PHRED-scaled score) with no
//! header the plugin needs: the columns are read by position (`AVI.pm:161`). A variant that
//! is one base replaced by one base takes the row matching its position and
//! alleles; every other variant gets nothing. The values are written at the
//! variant level, so they print on every row of the variant, as a plugin whose
//! Ensembl VEP `feature_types` are `Feature` and `Intergenic` prints them.
//!
//! # Parameters
//!
//! - `file=<path>` (or the path as the first positional parameter): the score file
//! - `raw=1`: also write the raw score
//! - `links=1`: also write a link to the variant in the AlphaGenome Atlas
//!
//! A flag is on when its value is `1`, `true` or `yes` in any letter case.
//!
//! # Output fields
//!
//! - `AVI_PHRED`: the PHRED-scaled quantile (0 to about 80, `AVI.pm:53-58`)
//! - `AVI_RAW`: the raw score (about -1.3 to +4.0, `AVI.pm:53-58`), under `raw=1`
//! - `AVI_LINK`: the Atlas URL of the matching row, under `links=1`
//!
//! `AVI.pm` is the plugin's Ensembl VEP counterpart in Ensembl/VEP_plugins release/116.

use std::collections::{BTreeMap, HashMap};

use indexmap::IndexMap;
use rayon::prelude::*;
use vep_core::assembly::Assembly;
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

use crate::annotation_store::{open_backend, AnnotationStore, Backend};
use crate::binary_store::BinaryAnnotator;
use crate::params::PluginParams;
use crate::tabix::allele_match::matched_variant_alleles;
use crate::tabix::{table_contig, BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};

/// The parameters `AVI.pm` reads from `--plugin AVI,...`.
const PARAMETERS: &[&str] = &["file", "raw", "links"];

const PHRED_FIELD: &str = "AVI_PHRED";
const RAW_FIELD: &str = "AVI_RAW";
const LINK_FIELD: &str = "AVI_LINK";

/// The header descriptions vep-rs writes, byte for byte the ones Ensembl VEP writes.
const PHRED_DESCRIPTION: &str = "PHRED-scaled quantile (0 to approximately 80). A phred of 20 \
     means the variant is among the top 1% most deleterious substitutions genome-wide. Higher \
     => more deleterious.";
const RAW_DESCRIPTION: &str = "Raw score (approximately -1.3 to +4.0). Higher => more deleterious.";
const LINK_DESCRIPTION: &str = "Link to the variant in the Google DeepMind AlphaGenome Atlas";

const ATLAS_URL: &str = "https://deepmind.google.com/science/alphagenome/atlas";

/// 0-based columns of the six-column file (`AVI.pm:161`): chromosome, position, reference
/// allele, alternate allele, raw score, PHRED score.
const POS_COL: usize = 1;
const REF_COL: usize = 2;
const ALT_COL: usize = 3;
const RAW_COL: usize = 4;
const PHRED_COL: usize = 5;

/// The score table as `init` opened it.
enum ScoreTable {
    /// The bgzipped table, read window by window: a batch's substitutions are
    /// single-base windows, and where a table holds rows between two
    /// substitutions, a read of the merged span between them would hand each
    /// of those rows back.
    Tabix(TabixAnnotator),
    /// A `.vpd` store beside the table, read by region.
    Binary(BinaryAnnotator),
}

impl ScoreTable {
    /// [`open_backend`]'s choice between the table and the `.vpd` store beside
    /// it, kept under its own types so the tabix table keeps its window read.
    fn open(
        tabix_path: &std::path::Path,
        config: TabixConfig,
        assembly: Option<Assembly>,
    ) -> Result<Self, PluginError> {
        let backend = open_backend(
            tabix_path,
            config,
            assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )?;
        Ok(match backend {
            Backend::Store(store) => Self::Binary(store),
            Backend::Table(table) => Self::Tabix(table),
        })
    }

    /// The rows at the positions of `substitutions`, indexed by chromosome and
    /// start for `annotate`'s lookup.
    fn rows_at<'a>(
        &self,
        substitutions: impl Iterator<Item = &'a InputVariant>,
        assembly: Option<Assembly>,
    ) -> Result<BatchQueryResult, PluginError> {
        match self {
            Self::Binary(store) => {
                let regions: Vec<(&str, u64, u64)> = substitutions
                    .map(|v| {
                        let contig = table_contig(v, assembly);
                        (contig, v.start, v.end)
                    })
                    .collect();
                store.query_batch(&regions)
            }
            Self::Tabix(tabix) => {
                let mut windows_by_chr: BTreeMap<&str, Vec<(u64, u64)>> = BTreeMap::new();
                for variant in substitutions {
                    windows_by_chr
                        .entry(table_contig(variant, assembly))
                        .or_default()
                        .push((variant.start, variant.start));
                }
                let mut records_by_pos = HashMap::new();
                for (chr, mut windows) in windows_by_chr {
                    windows.sort_unstable();
                    windows.dedup();
                    // The pipeline reads a batch ahead on one task per plugin, so
                    // the batch's blocks are inflated side by side on the current
                    // pool; the rows come back in file order.
                    let mut by_start: BTreeMap<u64, Vec<TabixRecord>> = BTreeMap::new();
                    for record in tabix.query_windows_parallel(chr, &windows)? {
                        by_start.entry(record.start).or_default().push(record);
                    }
                    if !by_start.is_empty() {
                        records_by_pos.insert(chr.to_string(), by_start);
                    }
                }
                Ok(BatchQueryResult::from_records(records_by_pos))
            }
        }
    }
}

/// AlphaGenome Variant Impact score plugin.
pub struct AviPlugin {
    /// The score table, opened by `init`.
    table: Option<ScoreTable>,
    /// Write `AVI_RAW` beside `AVI_PHRED` (`raw=1`).
    include_raw: bool,
    /// Write `AVI_LINK` beside `AVI_PHRED` (`links=1`).
    include_links: bool,
    /// The run's assembly (the cache's own, else the value of `--assembly`): the
    /// table is opened on it, the rows are asked for under it, and `init` warns
    /// of a GRCh37 run.
    assembly: Option<Assembly>,
}

impl Default for AviPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl AviPlugin {
    /// Create an uninitialized AVI plugin; `init` opens the file.
    pub fn new() -> Self {
        Self {
            table: None,
            include_raw: false,
            include_links: false,
            assembly: None,
        }
    }

    /// `AVI.pm` `_enabled`: a flag is on for `1`, `true` or `yes` in any case.
    fn option_enabled(value: Option<&str>) -> bool {
        value.is_some_and(|v| {
            ["1", "true", "yes"]
                .iter()
                .any(|on| v.eq_ignore_ascii_case(on))
        })
    }

    /// `AVI.pm` scores a variant whose start equals its end and whose reference
    /// and alternate alleles are each one of `A`, `C`, `G`, `T` (`run`, the
    /// three early returns); Ensembl VEP upper-cases the alleles before that
    /// test, so the case is ignored here.
    fn is_single_base_substitution(variant: &InputVariant) -> bool {
        let is_base = |allele: &[u8]| {
            allele.len() == 1 && matches!(allele[0].to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
        };
        variant.start == variant.end
            && is_base(&variant.ref_allele)
            && is_base(variant.alt_allele())
    }

    /// The row of `records` whose alleles match `variant` at its position
    /// (`get_matched_variant_alleles`, `AVI.pm` `run`), with its score columns
    /// present (`AVI.pm:163` drops a line without a sixth column
    /// before any match).
    fn matching_record<'a>(
        records: &[&'a TabixRecord],
        variant: &InputVariant,
    ) -> Option<&'a TabixRecord> {
        records
            .iter()
            .copied()
            .filter(|record| record.columns.len() > PHRED_COL)
            .find(|record| matched_variant_alleles(record, variant, POS_COL, REF_COL, ALT_COL))
    }

    /// Write the fields of `record` on `variant`: the PHRED score, the raw
    /// score under `raw=1` and the Atlas link under `links=1`. The link names the
    /// row's own chromosome, position and alleles, with `chr` prefixed when the
    /// file's name lacks it, as `AVI.pm` `parse_data` builds it. A score cell
    /// the row has and leaves empty is Ensembl VEP's empty string, which the
    /// seven-way split of `AVI.pm:161` keeps, so the record names its field.
    fn write_fields(&self, record: &TabixRecord, variant: &mut InputVariant) {
        self.write_cell(variant, PHRED_FIELD, &record.columns[PHRED_COL]);
        if self.include_raw {
            self.write_cell(variant, RAW_FIELD, &record.columns[RAW_COL]);
        }
        if self.include_links {
            let chr = &record.chr;
            let has_prefix = chr
                .get(..3)
                .is_some_and(|prefix| prefix.eq_ignore_ascii_case("chr"));
            let link_chr = if has_prefix {
                chr.clone()
            } else {
                format!("chr{chr}")
            };
            variant.set_plugin_value(
                LINK_FIELD.to_string(),
                format!(
                    "{ATLAS_URL}?q={link_chr}:{}:{}>{}&m=variant",
                    record.start, record.columns[REF_COL], record.columns[ALT_COL]
                ),
                false,
                self.feature_types(),
            );
        }
    }

    /// Write `cell` under `field` on `variant`, naming the field in the record's
    /// `plugin_empty_strings` when the cell is empty. A later plugin's value
    /// replaces an earlier one under the same key (OutputFactory.pm:2239-2240),
    /// its mark too.
    fn write_cell(&self, variant: &mut InputVariant, field: &str, cell: &str) {
        variant.set_plugin_value(
            field.to_string(),
            cell.to_string(),
            true,
            self.feature_types(),
        );
    }
}

impl crate::traits::sealed::Sealed for AviPlugin {}

impl BuiltinPlugin for AviPlugin {
    fn name(&self) -> &str {
        "AVI"
    }

    /// The Ensembl VEP module declares `['Feature', 'Intergenic']`: every row,
    /// the intergenic row included.
    fn feature_types(&self) -> &[RowKind] {
        RowKind::ALL
    }

    fn header_info(&self) -> Vec<(String, String)> {
        let mut fields = vec![(PHRED_FIELD.to_string(), PHRED_DESCRIPTION.to_string())];
        if self.include_raw {
            fields.push((RAW_FIELD.to_string(), RAW_DESCRIPTION.to_string()));
        }
        if self.include_links {
            fields.push((LINK_FIELD.to_string(), LINK_DESCRIPTION.to_string()));
        }
        fields
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn set_reference_fasta(&mut self, fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>) {
        match self.table.as_mut() {
            Some(ScoreTable::Tabix(table)) => table.frame_mut().set_reference(fasta),
            Some(ScoreTable::Binary(store)) => store.frame_mut().set_reference(fasta),
            None => {}
        }
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let parsed = PluginParams::parse(self.name(), PARAMETERS, params)?;
        // `$params->{file} || $positional_file`, then `unless $file` (AVI.pm:96-98):
        // a value Perl reads as false, `0` or the empty string, is no file.
        let path = parsed
            .get("file")
            .or_else(|| parsed.positional().first().map(String::as_str))
            .filter(|path| perl_true(path))
            .ok_or_else(|| {
                PluginError::Init(
                    "AVI: no data file given; pass file=<path> or the path as the first \
                     parameter"
                        .to_string(),
                )
            })?;
        self.include_raw = Self::option_enabled(parsed.get("raw"));
        self.include_links = Self::option_enabled(parsed.get("links"));

        let file_path = std::path::PathBuf::from(path);
        if !file_path.is_file() {
            return Err(PluginError::Init(format!(
                "AVI: data file {} not found",
                file_path.display()
            )));
        }
        if self.assembly == Some(Assembly::Grch37) {
            tracing::warn!(
                "AVI: no assembly is checked for this plugin; a GRCh37 run matches the file's rows \
                 by position and alleles as they stand, so the file must be on GRCh37 coordinates"
            );
        }

        let config = TabixConfig {
            file_path: file_path.clone(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: Some(REF_COL),
            alt_col: Some(ALT_COL),
            zero_based: false,
        };
        self.table = Some(ScoreTable::open(&file_path, config, self.assembly)?);
        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(table) = &self.table else {
            return Ok(None);
        };
        // A variant the gate refuses is never looked up, as `AVI.pm` returns
        // before `get_data` for it.
        let substitutions = variants
            .iter()
            .filter(|v| self.annotates_variant(v) && Self::is_single_base_substitution(v));
        Ok(Some(PrefetchData::new(
            table.rows_at(substitutions, self.assembly)?,
        )))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batch: BatchQueryResult = data.downcast()?;
        if batch.is_empty() {
            return Ok(());
        }
        // Each variant reads the shared batch and writes its own fields, so
        // the variants are annotated side by side on the current pool.
        variants.par_iter_mut().for_each(|variant| {
            if !self.annotates_variant(variant) || !Self::is_single_base_substitution(variant) {
                return;
            }
            let records = batch.lookup(
                table_contig(variant, self.assembly),
                variant.start,
                variant.end,
            );
            if let Some(record) = Self::matching_record(&records, variant) {
                self.write_fields(record, variant);
            }
        });
        // Dropping the batch's records costs as much as parsing them, so one
        // pool job frees them off the batch's path; freed on several threads
        // at once, they cost more cycles than on one.
        rayon::spawn(move || drop(batch));
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
        // Variant-level plugin: annotation happens in `annotate`.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use vep_core::coordinate::Strand;

    /// Made-up rows: six columns, no header, `chr`-prefixed names: two alleles at 5030192 and
    /// one row each at the positions of the non-substitution records below.
    const ROWS: &[&str] = &[
        "chr21\t5030192\tA\tG\t1.2345\t12.34",
        "chr21\t5030192\tA\tT\t-0.4321\t3.21",
        "chr21\t5030602\tA\tC\t0.5\t8.0",
        "chr21\t5030625\tG\tA\t0.6\t9.0",
        "chr21\t5032970\tT\tA\t0.7\t10.0",
        "chr21\t5131441\tT\tC\t2.75\t30.5",
    ];

    fn fixture(dir: &std::path::Path) -> std::path::PathBuf {
        write_tabix_fixture(
            dir,
            "avi_scores.synthetic.tsv.gz",
            ROWS,
            IndexSpec::standard(),
        )
    }

    fn plugin(params: &[String]) -> AviPlugin {
        let mut plugin = AviPlugin::new();
        plugin.init(params).expect("init");
        plugin
    }

    fn variant(pos: u64, reference: &[u8], alt: &[u8]) -> InputVariant {
        InputVariant::new("21".into(), pos, pos, reference.to_vec(), alt.to_vec())
    }

    fn fields(variant: &InputVariant) -> Vec<(&str, &str)> {
        variant
            .plugin_data
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect()
    }

    /// `--plugin AVI,file=<path>`: the PHRED score alone, from the row whose
    /// alleles match; a substitution the file lacks gets nothing.
    #[test]
    fn phred_only_by_default() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(&[format!("file={}", fixture(dir.path()).display())]);
        assert_eq!(
            plugin.header_info(),
            [(PHRED_FIELD.to_string(), PHRED_DESCRIPTION.to_string())]
        );

        let mut variants = vec![
            variant(5030192, b"A", b"G"),
            variant(5030192, b"A", b"T"),
            variant(5030192, b"A", b"C"),
            variant(5131441, b"T", b"C"),
            variant(6000000, b"A", b"G"),
        ];
        plugin.run_batch(&mut variants).expect("run_batch");

        assert_eq!(fields(&variants[0]), [("AVI_PHRED", "12.34")]);
        assert_eq!(fields(&variants[1]), [("AVI_PHRED", "3.21")]);
        assert!(fields(&variants[2]).is_empty(), "A>C has no row");
        assert_eq!(fields(&variants[3]), [("AVI_PHRED", "30.5")]);
        assert!(fields(&variants[4]).is_empty(), "no row at 6000000");
    }

    /// `AVI.pm:102` hands the file to the tabix base class, which needs
    /// `<file>.tbi` beside it: an index named on the stem is a refusal
    /// naming the missing `.tbi`, though htslib would find it.
    #[test]
    fn table_without_its_tbi_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(dir.path());
        std::fs::rename(
            dir.path().join("avi_scores.synthetic.tsv.gz.tbi"),
            dir.path().join("avi_scores.synthetic.tsv.tbi"),
        )
        .unwrap();
        let error = AviPlugin::new()
            .init(&[format!("file={}", path.display())])
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("avi_scores.synthetic.tsv.gz.tbi not found"),
            "{error}"
        );
    }

    /// `raw=1` adds the raw score and `links=1` the Atlas URL built from the
    /// row's own chromosome, position and alleles; a flag is on for `1`, `true`
    /// or `yes` in any case and off otherwise; the file may be positional.
    #[test]
    fn raw_and_links_follow_the_options() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(dir.path());
        let file = format!("file={}", path.display());

        let both = plugin(&[file.clone(), "raw=1".to_string(), "links=1".to_string()]);
        assert_eq!(
            both.header_info(),
            [
                (PHRED_FIELD.to_string(), PHRED_DESCRIPTION.to_string()),
                (RAW_FIELD.to_string(), RAW_DESCRIPTION.to_string()),
                (LINK_FIELD.to_string(), LINK_DESCRIPTION.to_string()),
            ]
        );
        let mut variants = vec![variant(5030192, b"A", b"T")];
        both.run_batch(&mut variants).unwrap();
        assert_eq!(
            fields(&variants[0]),
            [
                ("AVI_PHRED", "3.21"),
                ("AVI_RAW", "-0.4321"),
                (
                    "AVI_LINK",
                    "https://deepmind.google.com/science/alphagenome/atlas\
                     ?q=chr21:5030192:A>T&m=variant"
                ),
            ]
        );

        let raw_only = plugin(&[file.clone(), "raw=TRUE".to_string(), "links=0".to_string()]);
        let mut variants = vec![variant(5030192, b"A", b"G")];
        raw_only.run_batch(&mut variants).unwrap();
        assert_eq!(
            fields(&variants[0]),
            [("AVI_PHRED", "12.34"), ("AVI_RAW", "1.2345")]
        );

        let links_only = plugin(&[file, "raw=no".to_string(), "links=yes".to_string()]);
        assert!(links_only.include_links && !links_only.include_raw);

        let positional = plugin(&[path.display().to_string()]);
        assert!(!positional.include_raw && !positional.include_links);
        let mut variants = vec![variant(5131441, b"T", b"C")];
        positional.run_batch(&mut variants).unwrap();
        assert_eq!(fields(&variants[0]), [("AVI_PHRED", "30.5")]);

        let bare = write_tabix_fixture(
            dir.path(),
            "avi_bare.synthetic.tsv.gz",
            &["21\t5131441\tT\tC\t2.75\t30.5"],
            IndexSpec::standard(),
        );
        let bare_links = plugin(&[format!("file={}", bare.display()), "links=1".to_string()]);
        let mut variants = vec![variant(5131441, b"T", b"C")];
        bare_links.run_batch(&mut variants).unwrap();
        assert_eq!(
            variants[0].plugin_data.get("AVI_LINK").map(String::as_str),
            Some(
                "https://deepmind.google.com/science/alphagenome/atlas\
                 ?q=chr21:5131441:T>C&m=variant"
            )
        );
    }

    /// An insertion, a deletion, a two-base substitution and a substitution of
    /// `N` take no score even where the file has a row at their position.
    #[test]
    fn non_snv_rows_get_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(&[format!("file={}", fixture(dir.path()).display())]);

        let insertion =
            InputVariant::new("21".into(), 5030603, 5030602, b"-".to_vec(), b"A".to_vec());
        let deletion = InputVariant::new(
            "21".into(),
            5030626,
            5030628,
            b"CAT".to_vec(),
            b"-".to_vec(),
        );
        let two_bases = InputVariant::new(
            "21".into(),
            5032970,
            5032971,
            b"TC".to_vec(),
            b"AG".to_vec(),
        );
        let n_base = variant(5131441, b"N", b"C");

        let mut variants = vec![insertion, deletion, two_bases, n_base];
        plugin.run_batch(&mut variants).expect("run_batch");
        for variant in &variants {
            assert!(
                variant.plugin_data.is_empty(),
                "{} {}..{} got {:?}",
                variant.allele_string,
                variant.start,
                variant.end,
                variant.plugin_data
            );
        }
    }

    /// Ensembl VEP 116.3 scores a record whose start equals its end before any
    /// row is matched (`AVI.pm:134`), so a substitution written with an
    /// unchanged leading base (`CA/CG` at 5030191) takes no score although
    /// `get_matched_variant_alleles` (ensembl-variation Utils/Sequence.pm:1199-1208)
    /// would trim it to the `A G` row at 5030192, which the single-base `A/G`
    /// there takes.
    #[test]
    fn a_substitution_with_an_unchanged_leading_base_is_not_scored() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(&[format!("file={}", fixture(dir.path()).display())]);
        let leading = InputVariant::new(
            "21".into(),
            5030191,
            5030192,
            b"CA".to_vec(),
            b"CG".to_vec(),
        );
        let single = variant(5030192, b"A", b"G");
        let mut variants = vec![leading, single];
        plugin.run_batch(&mut variants).expect("run_batch");
        assert!(variants[0].plugin_data.is_empty());
        assert_eq!(fields(&variants[1]), [("AVI_PHRED", "12.34")]);
    }

    /// `AVI.pm:161` splits the line htslib hands over, the line end alone
    /// removed, with `split /\t/, $line, 7`, so the PHRED cell keeps the
    /// whitespace before the line end and a whitespace-only cell is a defined
    /// score (`:163`): `12.34 ` prints as `AVI_PHRED=12.34 ` and ` ` as
    /// `AVI_PHRED= `.
    #[test]
    fn a_trailing_space_on_the_last_cell_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "avi_trailing.synthetic.tsv.gz",
            &[
                "chr21\t5030192\tA\tG\t1.2345\t12.34 ",
                "chr21\t5030193\tC\tT\t0.5\t ",
            ],
            IndexSpec::standard(),
        );
        let plugin = plugin(&[format!("file={}", path.display())]);
        let mut variants = vec![variant(5030192, b"A", b"G"), variant(5030193, b"C", b"T")];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(fields(&variants[0]), [("AVI_PHRED", "12.34 ")]);
        assert_eq!(fields(&variants[1]), [("AVI_PHRED", " ")]);
    }

    /// `AVI.pm:161` splits the line into six variables with a limit of seven,
    /// so Ensembl VEP keeps a trailing empty cell as the empty string, and `:163`
    /// drops a row lacking its sixth cell before any match. An empty PHRED or
    /// raw cell prints as `AVI_PHRED=` in the default format, an empty tab cell
    /// against `-` and JSON `""` against no key (`Tab.pm:128`, `JSON.pm:327`),
    /// so the record names the field of a present empty cell in
    /// `plugin_empty_strings`; a row lacking the cell scores nothing, as the
    /// next matching row, if any, is read instead.
    #[test]
    fn a_row_lacking_its_score_cell_is_the_empty_value_and_an_empty_cell_the_empty_string() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "avi_cells.synthetic.tsv.gz",
            &[
                "chr21\t5030000\tA\tG\t1.2345",
                "chr21\t5030001\tA\tG\t1.2345\t",
                "chr21\t5030002\tA\tG\t1.2345\t\tx",
                "chr21\t5030003\tA\tG\t\t12.34",
                "chr21\t5030004\tA\tG\t1.2345\t12.34",
            ],
            IndexSpec::standard(),
        );
        let plugin = plugin(&[format!("file={}", path.display()), "raw=1".to_string()]);
        let marks = |pos: u64| {
            let mut variants = vec![variant(pos, b"A", b"G")];
            plugin.run_batch(&mut variants).unwrap();
            (
                fields(&variants[0])
                    .into_iter()
                    .map(|(k, v)| (k.to_string(), v.to_string()))
                    .collect::<Vec<_>>(),
                variants[0].plugin_empty_strings.clone(),
            )
        };
        let pair = |k: &str, v: &str| (k.to_string(), v.to_string());
        assert_eq!(marks(5030000), (vec![], vec![]), "five cells");
        assert_eq!(
            marks(5030001),
            (
                vec![pair("AVI_PHRED", ""), pair("AVI_RAW", "1.2345")],
                vec!["AVI_PHRED".to_string()]
            ),
            "six cells, the sixth empty"
        );
        assert_eq!(
            marks(5030002),
            (
                vec![pair("AVI_PHRED", ""), pair("AVI_RAW", "1.2345")],
                vec!["AVI_PHRED".to_string()]
            ),
            "seven cells, the sixth empty"
        );
        assert_eq!(
            marks(5030003),
            (
                vec![pair("AVI_PHRED", "12.34"), pair("AVI_RAW", "")],
                vec!["AVI_RAW".to_string()]
            ),
            "six cells, the fifth empty"
        );
        assert_eq!(
            marks(5030004),
            (
                vec![pair("AVI_PHRED", "12.34"), pair("AVI_RAW", "1.2345")],
                vec![]
            )
        );
    }

    /// The header texts are the ones Ensembl VEP writes for the three fields,
    /// byte for byte.
    #[test]
    fn header_descriptions_match_the_perl_text() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(&[
            format!("file={}", fixture(dir.path()).display()),
            "raw=1".to_string(),
            "links=1".to_string(),
        ]);
        let header: std::collections::BTreeMap<String, String> =
            plugin.header_info().into_iter().collect();
        assert_eq!(
            header.get("AVI_PHRED").map(String::as_str),
            Some(
                "PHRED-scaled quantile (0 to approximately 80). A phred of 20 means the variant \
                 is among the top 1% most deleterious substitutions genome-wide. Higher => more \
                 deleterious."
            )
        );
        assert_eq!(
            header.get("AVI_RAW").map(String::as_str),
            Some("Raw score (approximately -1.3 to +4.0). Higher => more deleterious.")
        );
        assert_eq!(
            header.get("AVI_LINK").map(String::as_str),
            Some("Link to the variant in the Google DeepMind AlphaGenome Atlas")
        );
    }

    /// A record given on the reverse strand is matched through its forward-strand
    /// alleles (`get_matched_variant_alleles` reverse-complements them), so `T/C -1`
    /// at 5030192 takes the `A G` row.
    #[test]
    fn a_reverse_strand_record_matches_through_its_forward_alleles() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(&[format!("file={}", fixture(dir.path()).display())]);
        let mut reverse = variant(5030192, b"T", b"C");
        reverse.strand = Strand::Reverse;
        let mut variants = vec![reverse];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(fields(&variants[0]), [("AVI_PHRED", "12.34")]);
    }

    /// A `.vpd` store beside the table is read instead of it, by region, and
    /// answers with the same fields.
    #[test]
    fn a_vpd_store_beside_the_table_is_read_by_region() {
        use crate::binary_store::{BinaryStoreWriter, ContentIdentity};
        use std::sync::Arc;

        let dir = tempfile::tempdir().unwrap();
        let path = fixture(dir.path());
        let header: Arc<Vec<String>> = Arc::new(
            ["chr", "pos", "ref", "alt", "raw", "phred"]
                .iter()
                .map(|c| c.to_string())
                .collect(),
        );
        let rows: Vec<TabixRecord> = ROWS
            .iter()
            .map(|line| {
                let columns: Vec<String> = line.split('\t').map(str::to_string).collect();
                let start: u64 = columns[1].parse().unwrap();
                TabixRecord {
                    chr: columns[0].clone(),
                    start,
                    end: start,
                    columns,
                    header: Some(Arc::clone(&header)),
                }
            })
            .collect();
        let vpd = path.with_extension("vpd");
        BinaryStoreWriter::new(vpd)
            .write(
                &header,
                &[("chr21".to_string(), rows)],
                &TabixConfig::default(),
                &ContentIdentity::of_file(&path).unwrap(),
            )
            .expect("write the binary store");

        let plugin = plugin(&[format!("file={}", path.display()), "raw=1".to_string()]);
        assert!(matches!(plugin.table, Some(ScoreTable::Binary(_))));
        let mut variants = vec![variant(5030192, b"A", b"T"), variant(5030192, b"A", b"C")];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(
            fields(&variants[0]),
            [("AVI_PHRED", "3.21"), ("AVI_RAW", "-0.4321")]
        );
        assert!(fields(&variants[1]).is_empty());
    }

    /// No file, an empty `file=`, a missing file and an unknown key are refused
    /// with a message opening `AVI: `.
    #[test]
    fn init_refusals_name_the_plugin() {
        let message = |params: &[&str]| {
            let params: Vec<String> = params.iter().map(|p| p.to_string()).collect();
            match AviPlugin::new().init(&params) {
                Err(PluginError::Init(message)) => message,
                other => panic!("expected an init error for {params:?}, got {other:?}"),
            }
        };
        assert!(message(&[]).starts_with("AVI: no data file given"));
        assert!(message(&["file="]).starts_with("AVI: no data file given"));
        // `$params->{file} || $positional_file` then `unless $file` (AVI.pm:96-98):
        // `0` is no file in either place.
        assert!(message(&["file=0"]).starts_with("AVI: no data file given"));
        assert!(message(&["0"]).starts_with("AVI: no data file given"));
        assert_eq!(
            message(&["file=/no/such/avi.synthetic.tsv.gz"]),
            "AVI: data file /no/such/avi.synthetic.tsv.gz not found"
        );
        assert!(message(&["file=/no/such/avi.synthetic.tsv.gz", "link=1"])
            .starts_with("AVI: unknown parameter 'link'"));
    }

    /// `init` on a GRCh37 run over the fixture warns once: no assembly is
    /// checked, and the file's rows are matched by position and alleles as
    /// they stand.
    #[test]
    fn a_grch37_run_is_warned_that_no_assembly_is_checked() {
        use crate::test_fixtures::Warnings;
        use std::sync::Arc;

        let dir = tempfile::tempdir().unwrap();
        let path = fixture(dir.path());
        let warnings = Arc::new(Warnings::default());
        let mut plugin = AviPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch37));
        tracing::subscriber::with_default(Arc::clone(&warnings), || {
            plugin
                .init(&[format!("file={}", path.display())])
                .expect("init");
        });
        assert_eq!(
            *warnings.0.lock().unwrap(),
            ["AVI: no assembly is checked for this plugin; a GRCh37 run matches the file's rows \
              by position and alleles as they stand, so the file must be on GRCh37 coordinates"]
        );
    }

    /// Ensembl VEP 116.3 never hands AVI a structural variant: the module keeps
    /// the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2229). A `<DEL>` gets no field, and
    /// nor does a scored substitution carrying the structural mark, which pins
    /// the gate rather than the allele test; the plain substitution is scored.
    #[test]
    fn a_structural_variant_gets_no_field() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(&[format!("file={}", fixture(dir.path()).display())]);
        assert!(!plugin.annotates_structural_variants());
        let mut deletion = InputVariant::new(
            "21".into(),
            5030192,
            5030300,
            b"N".to_vec(),
            b"<DEL>".to_vec(),
        );
        deletion.is_structural = true;
        let mut flagged = variant(5030192, b"A", b"G");
        flagged.is_structural = true;
        let mut variants = vec![deletion, flagged, variant(5030192, b"A", b"G")];
        plugin.run_batch(&mut variants).unwrap();
        assert!(fields(&variants[0]).is_empty());
        assert!(fields(&variants[1]).is_empty());
        assert_eq!(fields(&variants[2]), [("AVI_PHRED", "12.34")]);
    }

    /// A row's position cell reaches `get_matched_variant_alleles` as text and
    /// is keyed as text (`"$ref\_$alt\_$pos"`, ensembl-variation
    /// Utils/Sequence.pm:1232), since a one-base pair loses no leading base:
    /// `5030192` scores the variant at 5030192, `5030192.0` and `5030192.5` do
    /// not, and `5030193x` does not score the variant at 5030193. The tabix
    /// index places each at its leading integer, so every row is returned and
    /// the text alone decides.
    #[test]
    fn a_position_cell_is_matched_as_the_text_the_key_holds() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            "chr21\t5030192\tA\tG\t1.2\t12.0",
            "chr21\t5030192.0\tA\tT\t1.3\t13.0",
            "chr21\t5030192.5\tA\tC\t1.4\t14.0",
            "chr21\t5030193x\tA\tG\t1.5\t15.0",
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "cells.synthetic.tsv.gz",
            &rows,
            IndexSpec::standard(),
        );
        let plugin = plugin(&[format!("file={}", path.display())]);
        let mut variants = vec![
            variant(5030192, b"A", b"G"),
            variant(5030192, b"A", b"T"),
            variant(5030192, b"A", b"C"),
            variant(5030193, b"A", b"G"),
        ];
        plugin.run_batch(&mut variants).expect("run_batch");
        assert_eq!(fields(&variants[0]), [("AVI_PHRED", "12.0")]);
        assert!(
            fields(&variants[1]).is_empty(),
            "5030192.0 is not the text 5030192"
        );
        assert!(
            fields(&variants[2]).is_empty(),
            "5030192.5 is not the text 5030192"
        );
        assert!(
            fields(&variants[3]).is_empty(),
            "5030193x is not the text 5030193"
        );
    }
}
