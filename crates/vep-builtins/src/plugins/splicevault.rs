// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! SpliceVault plugin: the most common mis-splicing events around a splice
//! site predicted to be lost, read from the tabix-indexed SpliceVault table
//! (`SpliceVault_data_GRCh38.tsv.gz` or `SpliceVault_data_hg19.tsv.gz`,
//! the names of `SpliceVault.pm:75-80`, indexed on its position column).
//!
//! # Parameters
//!
//! Every parameter is `key=value`. Every key but `transcript_match` names a
//! data file (`SpliceVault.pm:119-122`), `file=<path>` by convention, and a
//! list naming none is refused (`:124-126`); each file is opened in the order
//! written, where Perl's hash order is unspecified, and a record's rows are
//! read file by file in that order (BaseVepTabixPlugin.pm:591-610). A path
//! Perl reads as false (`file=`, `file=0`) is refused before it is looked for
//! (`check_file`, BaseVepTabixPlugin.pm:269). `transcript_match=1` (the
//! default, `SpliceVault.pm:116-117`) annotates a transcript row from the rows
//! naming its stable id alone, and `transcript_match=0` from the first row
//! matching the allele whatever transcript it names (`:161-163`). The value is
//! read as Perl's truth test `if ($self->{transcript_match})` (`:161`) reads
//! it: `0` and the empty value turn the filter off, every other value (`2`,
//! `yes`, `0.0`) leaves it on. Under `--refseq` or `--merged` with
//! `transcript_match` on, the plugin warns on standard error that the match
//! may fail (`:129-131`).
//!
//! # Data file
//!
//! Twelve positional columns, `chrom pos ref alt feature site_type site_pos
//! spliceai_delta out_of_frame_events top_events sample_count max_depth`
//! (`:201`); the reader takes the alleles as a VCF writes them (an insertion
//! carries its anchor base), and a `#` header line is skipped as a comment. A
//! row matches a
//! record when its alleles and the record's forward-strand alleles, both
//! trimmed of their shared bases, are one change at one start
//! (`get_matched_variant_alleles`, `:165-178`, through
//! [`allele_match::matched_variant_alleles`]); the rows of the window from two
//! bases before the record's start to its end are read (`:158`). The
//! twelve-scalar split of `:201` keeps a trailing empty cell as the defined
//! empty string and leaves a cell past the row's end undefined, so a row of
//! fewer than four cells has no reference or alternate cell, and the match dies
//! on it (`trim_sequences`, ensembl-variation Utils/Sequence.pm:968, reached
//! from `:1230`): under `transcript_match=0` every row of the window is put to
//! the match in file order until one matches, so such a row before the record's
//! first match fails the record, which keeps none of the seven fields on any of
//! its transcript rows under the warning `Plugin 'SpliceVault' went wrong`,
//! while a row after the match is never read; under `transcript_match` the
//! feature test of `:161-163` runs first and skips it, since a row too short to
//! hold an allele cell holds no feature cell either. `top_events`
//! is rewritten with every space as `_` and every `;` as `:` (`:203-204`) and
//! split on `|` into one event per element (`:216`); the field carries that
//! list with its `|` kept, so an event's text stays whole whatever it holds,
//! and the writers join it as `convert_arrayref` joins a list (Utils.pm:185-197):
//! with `,` in the default and tab formats (`OutputFactory/VEP_output.pm:117`,
//! `Tab.pm:128`) and `&` in the VCF (`OutputFactory/VCF.pm:391`). In the JSON
//! format the seven fields nest under `splicevault`, with the events keyed by
//! rank (`:180-191`); the JSON writer builds that object from the fields
//! written here.
//!
//! # Output fields
//!
//! The seven `SpliceVault_*` fields of `:142-148`, per transcript row (the
//! plugin declares `['Transcript']`, `:136-138`): `SpliceVault_top_events`,
//! `SpliceVault_site_sample_count`, `SpliceVault_site_max_depth`,
//! `SpliceVault_out_of_frame_events`, `SpliceVault_site_pos`,
//! `SpliceVault_site_type` and `SpliceVault_SpliceAI_delta`.
//!
//! Ensembl VEP compares a row's alleles with every alternate allele of the
//! record (`:156`), so a row naming a sibling allele annotates every allele's
//! rows there; here a row annotates the allele it names.
//!
//! `SpliceVault.pm` is the plugin's Ensembl VEP counterpart in
//! Ensembl/VEP_plugins release/116.

use std::path::Path;

use tracing::debug;

use crate::annotation_store::{open_best_store_with, AnnotationStore};
use crate::params::{Grammar, PluginParams};
use crate::tabix::{
    allele_match, table_contig, BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord,
};
use crate::traits::{warn_on_stderr, BuiltinPlugin, PluginError, PrefetchData};
use vep_core::assembly::Assembly;
use vep_core::consequence::TranscriptConsequence;
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

const NAME: &str = "SpliceVault";

/// 0-based positional columns (`SpliceVault.pm:201`).
const CHR_COL: usize = 0;
const START_COL: usize = 1;
const REF_COL: usize = 2;
const ALT_COL: usize = 3;
const FEATURE_COL: usize = 4;
const SITE_TYPE_COL: usize = 5;
const SITE_POS_COL: usize = 6;
const SPLICEAI_DELTA_COL: usize = 7;
const OUT_OF_FRAME_COL: usize = 8;
const TOP_EVENTS_COL: usize = 9;
const SAMPLE_COUNT_COL: usize = 10;
const MAX_DEPTH_COL: usize = 11;

/// The output fields with the column each is read from, in the order
/// `SpliceVault.pm:211-219` lists them.
const FIELDS: [(&str, usize); 7] = [
    ("SpliceVault_site_sample_count", SAMPLE_COUNT_COL),
    ("SpliceVault_site_max_depth", MAX_DEPTH_COL),
    ("SpliceVault_SpliceAI_delta", SPLICEAI_DELTA_COL),
    ("SpliceVault_out_of_frame_events", OUT_OF_FRAME_COL),
    ("SpliceVault_top_events", TOP_EVENTS_COL),
    ("SpliceVault_site_pos", SITE_POS_COL),
    ("SpliceVault_site_type", SITE_TYPE_COL),
];

/// The bases before a record's start the window opens with (`SpliceVault.pm:158`).
const WINDOW_BEFORE_START: u64 = 2;

/// One data file (`SpliceVault.pm:122`) with its annotation store (tabix or
/// binary).
struct Source {
    /// The file as given, named by the error of a row the match dies on.
    file: String,
    store: Box<dyn AnnotationStore>,
}

/// A row of a record's window with the file it came from and its allele match
/// ([`SpliceVaultPlugin::allele_matches`]).
struct Row<'a> {
    file: &'a str,
    record: &'a TabixRecord,
    matched: Option<bool>,
}

/// SpliceVault plugin.
pub struct SpliceVaultPlugin {
    /// The data files in parameter order, opened by `init`.
    sources: Vec<Source>,
    /// Whether a row annotates the transcript it names alone (`SpliceVault.pm:161-163`):
    /// `transcript_match` as Perl's truth test reads it, on when absent (`:117`).
    transcript_match: bool,
    /// Whether the run reads the RefSeq or merged transcript set.
    refseq_or_merged: bool,
    /// The run's assembly, which names a mitochondrial record's contig to the
    /// table ([`table_contig`]).
    assembly: Option<Assembly>,
}

impl Default for SpliceVaultPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl SpliceVaultPlugin {
    /// Create an uninitialized SpliceVault plugin.
    pub fn new() -> Self {
        Self {
            sources: Vec::new(),
            transcript_match: true,
            refseq_or_merged: false,
            assembly: None,
        }
    }

    /// Opens one data file as `add_file` does (BaseVepTabixPlugin.pm:247-273):
    /// a path Perl reads as false is refused as no file, then a missing one by
    /// name.
    fn open_source(
        key: &str,
        file: &str,
        assembly: Option<Assembly>,
    ) -> Result<Source, PluginError> {
        if !perl_true(file) {
            return Err(PluginError::Init(format!(
                "{NAME}: {key}= names no file; use {key}=<path>"
            )));
        }
        let path = Path::new(file);
        if !path.is_file() {
            return Err(PluginError::Init(format!("{NAME}: file {file} not found")));
        }
        let config = TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: CHR_COL,
            start_col: START_COL,
            end_col: None,
            ref_col: Some(REF_COL),
            alt_col: Some(ALT_COL),
            zero_based: false,
        };
        // A row's last cell is written as `split /\t/` reads it off the line
        // (`SpliceVault.pm:201`), its trailing whitespace included.
        let store = open_best_store_with(
            path,
            config,
            assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )?;
        Ok(Source {
            file: file.to_string(),
            store,
        })
    }

    /// The window `SpliceVault.pm:158` reads for a record: two bases before its
    /// start to its end, which for an insertion (start = end + 1) is the base
    /// before its anchor base and the anchor base.
    fn window(variant: &InputVariant) -> (u64, u64) {
        (
            variant.start.saturating_sub(WINDOW_BEFORE_START),
            variant.end,
        )
    }

    /// Each row of `records`, file by file and in file order, with whether its
    /// alleles and position match the record's, both sides trimmed of their
    /// shared bases (`SpliceVault.pm:165-178`), or `None` for a row the match
    /// dies on: `get_matched_variant_alleles` trims the row's reference and
    /// alternate cells (ensembl-variation Utils/Sequence.pm:1230) and
    /// `trim_sequences` throws when either is undefined (`:968`). The
    /// twelve-scalar split of `:201` keeps a trailing empty cell as the defined
    /// empty string, which trims to `-` and compares, so a cell is undefined
    /// only when the row ends before it: a row of fewer than three cells has no
    /// reference, one of fewer than four no alternate. Compared once per
    /// record, since the match reads nothing of the transcript.
    fn allele_matches<'a>(
        records: &[(&'a str, Vec<&'a TabixRecord>)],
        variant: &InputVariant,
    ) -> Vec<Row<'a>> {
        records
            .iter()
            .flat_map(|(file, rows)| rows.iter().map(move |&record| (*file, record)))
            .map(|(file, record)| {
                let defined =
                    record.columns.get(REF_COL).is_some() && record.columns.get(ALT_COL).is_some();
                // A Perl-false position cell returns the match before an undefined
                // allele dies (Sequence.pm:1173 before :1230).
                let matched = if allele_match::perl_false_position(record, START_COL) {
                    Some(false)
                } else {
                    defined.then(|| {
                        allele_match::matched_variant_alleles(
                            record, variant, START_COL, REF_COL, ALT_COL,
                        )
                    })
                };
                Row {
                    file,
                    record,
                    matched,
                }
            })
            .collect()
    }

    /// The row annotating `consequence`, the rows read in file order as the loop
    /// of `SpliceVault.pm:160-196` reads them: under `transcript_match` a row
    /// whose feature column is not the transcript's stable id is skipped before
    /// the match (`:161-163`; a row too short to hold an allele cell holds no
    /// feature cell either, so it is skipped here and never reaches the match);
    /// the first row matching the record's alleles is returned (`:178-193`) and
    /// the rows after it are not read. A row the match dies on
    /// ([`Self::allele_matches`]) reached before that match is the error, as
    /// Ensembl VEP's `run` dies there for the transcript row.
    fn matching_row<'a>(
        &self,
        rows: &[Row<'a>],
        consequence: &TranscriptConsequence,
    ) -> Result<Option<&'a TabixRecord>, PluginError> {
        for row in rows {
            if self.transcript_match
                && !row
                    .record
                    .columns
                    .get(FEATURE_COL)
                    .is_some_and(|feature| **feature == *consequence.transcript_id())
            {
                continue;
            }
            match row.matched {
                None => {
                    return Err(allele_match::undefined_allele_row(
                        NAME, row.file, row.record,
                    ))
                }
                Some(true) => return Ok(Some(row.record)),
                Some(false) => {}
            }
        }
        Ok(None)
    }

    /// Writes the seven output fields of a row (`SpliceVault.pm:206-220`) on the
    /// transcript row `consequence`, naming in its `plugin_empty_strings` the
    /// fields whose empty value is Ensembl VEP's empty string. The
    /// twelve-variable split of `:201` keeps an empty cell as the empty string
    /// and leaves a cell the row lacks undefined; both read as empty here, the
    /// former named. The events column is a list either way (`:216`), which
    /// `Tab.pm:128` prints as its join, an empty cell when no event is left, so
    /// its field is named whenever its list is empty, the cell present or not. A
    /// later plugin's value replaces an earlier one under the same key
    /// (OutputFactory.pm:2239-2240), its mark too.
    fn write_fields(record: &TabixRecord, consequence: &mut TranscriptConsequence) {
        consequence
            .plugin_empty_strings
            .retain(|key| !FIELDS.iter().any(|&(field, _)| field == key));
        for &(field, column) in &FIELDS {
            let cell = record.columns.get(column).map(String::as_str);
            let (value, named) = if column == TOP_EVENTS_COL {
                let events = top_events_list(cell.unwrap_or(""));
                let no_event = events.is_empty();
                (events, no_event)
            } else {
                (cell.unwrap_or("").to_string(), cell == Some(""))
            };
            if named {
                consequence.plugin_empty_strings.push(field.to_string());
            }
            consequence.plugin_data.insert(field.to_string(), value);
        }
    }

    /// Writes the matching row's fields on every transcript row of `variant`, or
    /// none of them when a transcript row's read reaches a row the match dies on
    /// ([`Self::matching_row`]), the error returned: under `transcript_match=0`
    /// every transcript row of the record reads the same rows and dies alike, so
    /// the first annotated row's read decides before anything is written, and
    /// under `transcript_match` none reaches the die.
    fn annotate_variant(
        &self,
        variant: &mut InputVariant,
        records: &[(&str, Vec<&TabixRecord>)],
    ) -> Result<(), PluginError> {
        let rows = Self::allele_matches(records, variant);
        for consequence in variant.transcript_consequences.iter_mut() {
            if !self.annotates_row(consequence) {
                continue;
            }
            if let Some(record) = self.matching_row(&rows, consequence)? {
                Self::write_fields(record, consequence);
            }
        }
        Ok(())
    }
}

/// The `top_events` column as the list `SpliceVault.pm:216` makes of it: every
/// space `_` and every `;` `:` (`:203-204`), the events kept `|`-separated
/// (Perl's `split` dropping the trailing empty ones, which are the trailing `|`
/// of the cell). The `|` is the column's own separator, so no event holds one
/// and the list reads back whatever an event's text carries; the writers join it
/// with their format's separator.
fn top_events_list(cell: &str) -> String {
    cell.trim_end_matches('|')
        .chars()
        .map(|c| match c {
            ' ' => '_',
            ';' => ':',
            other => other,
        })
        .collect()
}

impl crate::traits::sealed::Sealed for SpliceVaultPlugin {}

impl BuiltinPlugin for SpliceVaultPlugin {
    fn name(&self) -> &str {
        NAME
    }

    fn header_info(&self) -> Vec<(String, String)> {
        [
            (
                "SpliceVault_top_events",
                "List of SpliceVault top events with the following colon-separated fields per event: rank:type:transcript_impact:percent_of_supporting_samples:frame",
            ),
            (
                "SpliceVault_site_sample_count",
                "Number of SpliceVault annotated splicing samples",
            ),
            (
                "SpliceVault_site_max_depth",
                "Maximum number of reads seen in any one sample representing annotated splicing in GTEx",
            ),
            (
                "SpliceVault_out_of_frame_events",
                "Fraction of SpliceVault top events that cause a frameshift",
            ),
            (
                "SpliceVault_site_pos",
                "Genomic position of splice-site predicted to be lost by SpliceAI",
            ),
            (
                "SpliceVault_site_type",
                "Type of splice-site predicted to be lost by SpliceAI",
            ),
            (
                "SpliceVault_SpliceAI_delta",
                "SpliceVault-provided SpliceAI delta score",
            ),
        ]
        .into_iter()
        .map(|(field, description)| (field.to_string(), description.to_string()))
        .collect()
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn set_reference_fasta(&mut self, fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>) {
        for source in &mut self.sources {
            source.store.set_reference_fasta(fasta.clone());
        }
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        // `params_to_hash` reads the list once any parameter carries `=`
        // (BaseVepPlugin.pm:123-128); a list without one is the empty hash,
        // which `SpliceVault.pm:124-126` refuses as no file.
        let params =
            PluginParams::parse_any_key(NAME, Grammar::KeyValueOnly { file_key: "file" }, params)?;
        self.transcript_match = params.get("transcript_match").is_none_or(perl_true);
        // Every key but `transcript_match` names a data file (`:119-122`).
        let files: Vec<&(String, String)> = params
            .keyed()
            .iter()
            .filter(|(key, _)| key != "transcript_match")
            .collect();
        if files.is_empty() {
            return Err(PluginError::Init(format!(
                "{NAME}: no file specified; use file=<path>"
            )));
        }
        for (key, file) in files {
            self.sources
                .push(Self::open_source(key, file, self.assembly)?);
        }
        // After the file checks, as `SpliceVault.pm:127-131` orders them.
        if self.refseq_or_merged && self.transcript_match {
            warn_on_stderr(
                "SpliceVault transcript_match may fail with --refseq/--merged; consider \
                 transcript_match=0",
            );
        }
        debug!(
            files = ?self.sources.iter().map(|source| &source.file).collect::<Vec<_>>(),
            transcript_match = self.transcript_match,
            "SpliceVault plugin initialized"
        );
        Ok(())
    }

    fn set_transcript_set(&mut self, refseq: bool, merged: bool) {
        self.refseq_or_merged = refseq || merged;
    }

    fn prefetch_reads_consequences(&self) -> bool {
        true
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        if self.sources.is_empty() {
            return Ok(None);
        }
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|v| {
                self.annotates_variant(v)
                    && v.transcript_consequences
                        .iter()
                        .any(|tc| self.annotates_row(tc))
            })
            .map(|v| {
                let (start, end) = Self::window(v);
                (table_contig(v, self.assembly), start, end)
            })
            .collect();
        if regions.is_empty() {
            return Ok(Some(PrefetchData::new(Vec::<BatchQueryResult>::new())));
        }
        // One batch per file, in file order (BaseVepTabixPlugin.pm:591-610).
        let batches = self
            .sources
            .iter()
            .map(|source| source.store.query_batch(&regions))
            .collect::<Result<Vec<BatchQueryResult>, PluginError>>()?;
        Ok(Some(PrefetchData::new(batches)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batches: Vec<BatchQueryResult> = data.downcast()?;
        if batches.iter().all(BatchQueryResult::is_empty) {
            return Ok(());
        }
        let mut failure = None;
        for variant in variants.iter_mut().filter(|v| self.annotates_variant(v)) {
            let (start, end) = Self::window(variant);
            let records: Vec<(&str, Vec<&TabixRecord>)> = self
                .sources
                .iter()
                .zip(&batches)
                .map(|(source, batch)| {
                    (
                        source.file.as_str(),
                        batch.lookup(table_contig(variant, self.assembly), start, end),
                    )
                })
                .collect();
            if records.iter().all(|(_, rows)| rows.is_empty()) {
                continue;
            }
            if let Err(error) = self.annotate_variant(variant, &records) {
                failure.get_or_insert(error);
            }
        }
        match failure {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        if let Some(data) = self.prefetch(variants)? {
            self.annotate(variants, data)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use std::sync::Arc;
    use vep_core::consequence::FeatureType;
    use vep_core::transcript::Transcript;

    /// The made-up fixture's header line: the twelve columns in the order the plugin
    /// splits a row into (SpliceVault.pm:201).
    const HEADER: &str = "#chrom\tpos\tref\talt\ttranscript_id\tsplicevault_site_type\tsplicevault_site_loc\tspliceai_delta\tsplicevault_out_of_frame_events\tsplicevault_top4_events\tsplicevault_site_sample_count\tsplicevault_site_max_depth";

    /// One made-up row of those twelve columns.
    fn row(
        pos: u64,
        reference: &str,
        alt: &str,
        transcript: &str,
        top_events: &str,
        delta: &str,
    ) -> String {
        [
            "chr21",
            &pos.to_string(),
            reference,
            alt,
            transcript,
            "Donor_loss",
            "chr21:5030031",
            delta,
            "Frameshift:2/2",
            top_events,
            "16054",
            "17",
        ]
        .join("\t")
    }

    fn fixture(dir: &Path, rows: &[String]) -> String {
        let mut lines = vec![HEADER];
        lines.extend(rows.iter().map(String::as_str));
        write_tabix_fixture(
            dir,
            "synthetic_table_of_events.tsv.gz",
            &lines,
            IndexSpec::standard(),
        )
        .display()
        .to_string()
    }

    fn loaded(params: &[String]) -> SpliceVaultPlugin {
        let mut plugin = SpliceVaultPlugin::new();
        plugin.init(params).expect("fixture loads");
        plugin
    }

    fn transcript_row(id: &str) -> TranscriptConsequence {
        let mut consequence = TranscriptConsequence::default();
        consequence.transcript = Some(Arc::new({
            let mut transcript = Transcript::default();
            transcript.stable_id = id.into();
            transcript
        }));
        consequence
    }

    /// A record on chromosome 21 with one transcript row per id.
    fn variant(
        start: u64,
        end: u64,
        reference: &[u8],
        alt: &[u8],
        transcripts: &[&str],
    ) -> InputVariant {
        let mut variant =
            InputVariant::new("21".into(), start, end, reference.to_vec(), alt.to_vec());
        variant
            .transcript_consequences
            .extend(transcripts.iter().map(|id| transcript_row(id)));
        variant
    }

    /// The value of `field` on every row of a record, in row order.
    fn values(
        plugin: &SpliceVaultPlugin,
        variant: InputVariant,
        field: &str,
    ) -> Vec<Option<String>> {
        let mut variants = vec![variant];
        plugin.run_batch(&mut variants).unwrap();
        assert!(
            variants[0].plugin_data.is_empty(),
            "SpliceVault writes per transcript row"
        );
        variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.get(field).cloned())
            .collect()
    }

    fn init_message(result: Result<(), PluginError>) -> String {
        match result {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected PluginError::Init, got {other:?}"),
        }
    }

    /// A made-up top-events cell of two events.
    const EVENTS: &str = "Top1:CD;-35;31%;Frameshift|Top2:CD;+48;0.006%;Frameshift";

    /// `SpliceVault.pm:119-122` reads every key but `transcript_match` as a data
    /// file: `data=` annotates as `file=` does, and two file keys read both
    /// files, a record's rows taken file by file in the order written
    /// (BaseVepTabixPlugin.pm:591-610), the error of a row the match dies on
    /// naming its file. A list naming no file is refused (`SpliceVault.pm:124-126`),
    /// as a path Perl reads as false is (`check_file`, BaseVepTabixPlugin.pm:269), a
    /// missing file by name (BaseVepTabixPlugin.pm:272), and a bare path beside a
    /// `key=value` as `params_to_hash` refuses it (BaseVepPlugin.pm:127).
    #[test]
    fn any_key_but_transcript_match_names_a_data_file() {
        let dir = tempfile::tempdir().unwrap();
        let first = write_tabix_fixture(
            dir.path(),
            "first_synthetic.tsv.gz",
            &[
                HEADER,
                &row(5030000, "G", "T", "ENST00000300305", EVENTS, "0.36"),
            ],
            IndexSpec::standard(),
        )
        .display()
        .to_string();
        let second = write_tabix_fixture(
            dir.path(),
            "second_synthetic.tsv.gz",
            &[
                HEADER,
                &row(5030000, "G", "T", "ENST00000300305", EVENTS, "0.52"),
                &row(5030100, "C", "A", "ENST00000300305", EVENTS, "0.77"),
            ],
            IndexSpec::standard(),
        )
        .display()
        .to_string();
        let delta = |plugin: &SpliceVaultPlugin, start: u64, reference: &[u8], alt: &[u8]| {
            values(
                plugin,
                variant(start, start, reference, alt, &["ENST00000300305"]),
                "SpliceVault_SpliceAI_delta",
            )
        };

        let data = loaded(&[format!("data={first}")]);
        assert_eq!(data.sources.len(), 1);
        assert_eq!(
            delta(&data, 5030000, b"G", b"T"),
            [Some("0.36".to_string())]
        );
        let two = loaded(&[format!("file={first}"), format!("file2={second}")]);
        assert_eq!(two.sources.len(), 2);
        assert_eq!(delta(&two, 5030000, b"G", b"T"), [Some("0.36".to_string())]);
        assert_eq!(delta(&two, 5030100, b"C", b"A"), [Some("0.77".to_string())]);
        let reversed = loaded(&[format!("file={second}"), format!("file2={first}")]);
        assert_eq!(
            delta(&reversed, 5030000, b"G", b"T"),
            [Some("0.52".to_string())]
        );

        // A made-up row that ends after its reference cell.
        let short = TabixRecord {
            chr: "21".to_string(),
            start: 5030000,
            end: 5030000,
            columns: vec!["chr21".into(), "5030000".into(), "G".into()],
            header: None,
        };
        let snv = variant(5030000, 5030000, b"G", b"T", &["ENST00000300305"]);
        let plugin = loaded(&[format!("file={first}"), "transcript_match=0".into()]);
        let rows =
            SpliceVaultPlugin::allele_matches(&[("second_synthetic.tsv.gz", vec![&short])], &snv);
        assert_eq!(
            plugin
                .matching_row(&rows, &snv.transcript_consequences[0])
                .unwrap_err()
                .to_string(),
            "plugin execution failed: SpliceVault: the row of second_synthetic.tsv.gz at \
             21:5030000 has no reference or alternate allele cell"
        );

        let message = |params: &[String]| init_message(SpliceVaultPlugin::new().init(params));
        for params in [
            vec![first.clone()],
            vec![],
            vec!["transcript_match=0".to_string()],
        ] {
            assert_eq!(
                message(&params),
                "SpliceVault: no file specified; use file=<path>",
                "{params:?}"
            );
        }
        for value in ["", "0"] {
            assert_eq!(
                message(&[format!("file={value}")]),
                "SpliceVault: file= names no file; use file=<path>",
                "file={value}"
            );
        }
        assert_eq!(
            message(&[format!("file={first}"), "data=0".to_string()]),
            "SpliceVault: data= names no file; use data=<path>"
        );
        let missing = dir
            .path()
            .join("absent_synthetic.tsv.gz")
            .display()
            .to_string();
        assert_eq!(
            message(&[format!("file={missing}")]),
            format!("SpliceVault: file {missing} not found")
        );
        assert_eq!(
            message(&[format!("file={first}"), first.clone()]),
            format!(
                "SpliceVault: parameter '{first}' is not key=value while 'file={first}' is; \
                 give every parameter as key=value"
            )
        );
        let on = loaded(&[format!("file={first}"), "transcript_match=1".into()]);
        assert!(on.transcript_match);
        let off = loaded(&[format!("file={first}"), "transcript_match=0".into()]);
        assert!(!off.transcript_match);
    }

    /// `transcript_match` is read as Perl's truth test `if ($self->{transcript_match})`
    /// (`SpliceVault.pm:161`) reads it, `1` when absent (`:117`): `2`, `yes` and
    /// `0.0` keep the feature-column filter on as `1` does, while `0` and the
    /// empty value turn it off, so every transcript row takes the first allele
    /// match in file order.
    #[test]
    fn transcript_match_is_read_by_perl_truth() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[
                row(5030000, "G", "T", "ENST00000000001", EVENTS, "0.11"),
                row(5030000, "G", "T", "ENST00000300305", EVENTS, "0.36"),
            ],
        );
        let delta = |value: &str| {
            let plugin = loaded(&[format!("file={path}"), format!("transcript_match={value}")]);
            let record = variant(5030000, 5030000, b"G", b"T", &["ENST00000300305"]);
            let deltas = values(&plugin, record, "SpliceVault_SpliceAI_delta");
            (plugin.transcript_match, deltas)
        };
        for on in ["1", "2", "yes", "0.0", "01", "false"] {
            assert_eq!(
                delta(on),
                (true, vec![Some("0.36".to_string())]),
                "transcript_match={on}"
            );
        }
        for off in ["0", ""] {
            assert_eq!(
                delta(off),
                (false, vec![Some("0.11".to_string())]),
                "transcript_match={off}"
            );
        }
    }

    /// `SpliceVault.pm:129` warns under `--refseq` or `--merged` alike, and the
    /// plugin loads under either set as under the Ensembl set; the line itself
    /// is read off the binary's standard error by the command-line tests.
    #[test]
    fn either_transcript_set_flag_is_the_warning_condition() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[row(5030000, "G", "T", "ENST00000300305", EVENTS, "0.36")],
        );
        for (refseq, merged, other_set) in [
            (true, false, true),
            (false, true, true),
            (true, true, true),
            (false, false, false),
        ] {
            let mut plugin = SpliceVaultPlugin::new();
            plugin.set_transcript_set(refseq, merged);
            assert_eq!(plugin.refseq_or_merged, other_set, "{refseq} {merged}");
            plugin.init(&[format!("file={path}")]).expect("loads");
            assert!(plugin.transcript_match);
            assert_eq!(plugin.sources.len(), 1);
        }
    }

    /// `SpliceVault.pm:161-163`: with `transcript_match` on, a transcript row is
    /// annotated from the row naming its stable id in the feature column, a row
    /// naming another id is skipped, and a transcript no row names gets nothing.
    #[test]
    fn transcript_match_default_on_filters_on_the_feature_column() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[
                row(5030000, "G", "T", "ENST00000000001", EVENTS, "0.11"),
                row(5030000, "G", "T", "ENST00000300305", EVENTS, "0.36"),
                row(5030000, "G", "T", "ENST00000346798", EVENTS, "0.47"),
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let record = variant(
            5030000,
            5030000,
            b"G",
            b"T",
            &["ENST00000346798", "ENST00000300305", "ENST00000999999"],
        );
        assert_eq!(
            values(&plugin, record, "SpliceVault_SpliceAI_delta"),
            [Some("0.47".to_string()), Some("0.36".to_string()), None]
        );
        let other_allele = variant(5030000, 5030000, b"G", b"A", &["ENST00000300305"]);
        assert_eq!(
            values(&plugin, other_allele, "SpliceVault_SpliceAI_delta"),
            [None]
        );
    }

    /// With `transcript_match=0` the feature column is not read: every
    /// transcript row takes the first row, in file order, whose alleles and
    /// position match the record (`SpliceVault.pm:165-178`), an insertion's
    /// anchor base trimmed as `get_matched_variant_alleles` trims it.
    #[test]
    fn transcript_match_0_takes_the_first_allele_match() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[
                row(5030000, "G", "A", "ENST00000300305", EVENTS, "0.22"),
                row(5030000, "G", "T", "ENST00000000001", EVENTS, "0.11"),
                row(5030000, "G", "T", "ENST00000300305", EVENTS, "0.36"),
                row(5030000, "G", "GT", "ENST00000300305", EVENTS, "0.45"),
            ],
        );
        let plugin = loaded(&[format!("file={path}"), "transcript_match=0".into()]);
        let snv = variant(
            5030000,
            5030000,
            b"G",
            b"T",
            &["ENST00000346798", "ENST00000300305"],
        );
        assert_eq!(
            values(&plugin, snv, "SpliceVault_SpliceAI_delta"),
            [Some("0.11".to_string()), Some("0.11".to_string())]
        );
        let insertion = variant(5030001, 5030000, b"-", b"T", &["ENST00000346798"]);
        assert_eq!(
            values(&plugin, insertion, "SpliceVault_SpliceAI_delta"),
            [Some("0.45".to_string())]
        );
        let with_match = loaded(&[format!("file={path}")]);
        let insertion = variant(5030001, 5030000, b"-", b"T", &["ENST00000346798"]);
        assert_eq!(
            values(&with_match, insertion, "SpliceVault_SpliceAI_delta"),
            [None]
        );
    }

    /// Ensembl VEP 116.3 trims the record's alleles as well as the row's inside
    /// `get_matched_variant_alleles` (ensembl-variation Utils/Sequence.pm:1199-1208
    /// and :1223-1232), so a substitution written with an unchanged leading
    /// base (`GG/GA` at 5029999) or trailing base (`GT/AT` at 5030000) takes
    /// the `G A` row at 5030000; a substitution of another base takes nothing.
    #[test]
    fn a_substitution_with_an_unchanged_flanking_base_matches_its_single_base_row() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[row(5030000, "G", "A", "ENST00000300305", EVENTS, "0.22")],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let leading = variant(5029999, 5030000, b"GG", b"GA", &["ENST00000300305"]);
        assert_eq!(
            values(&plugin, leading, "SpliceVault_SpliceAI_delta"),
            [Some("0.22".to_string())]
        );
        let trailing = variant(5030000, 5030001, b"GT", b"AT", &["ENST00000300305"]);
        assert_eq!(
            values(&plugin, trailing, "SpliceVault_SpliceAI_delta"),
            [Some("0.22".to_string())]
        );
        let other_base = variant(5029999, 5030000, b"GG", b"GT", &["ENST00000300305"]);
        assert_eq!(
            values(&plugin, other_base, "SpliceVault_SpliceAI_delta"),
            [None]
        );
    }

    /// `SpliceVault.pm:203-204` rewrite every space as `_` and every `;` as `:`
    /// before the events are split on `|` (`:216`); the field carries that list
    /// `|`-separated, every other character of an event kept, a `,` included,
    /// which the text formats' `,` join (Utils.pm:185-197) would otherwise make
    /// an event boundary.
    #[test]
    fn top_events_rewrite_spaces_and_semicolons() {
        // The cells below are made up.
        assert_eq!(
            top_events_list(
                "Top1:CD;+284;15 840%;Frameshift|Top2:CD;+332;10 028%;Frameshift|Top3:ES;5-6;31%;inFrame"
            ),
            "Top1:CD:+284:15_840%:Frameshift|Top2:CD:+332:10_028%:Frameshift|Top3:ES:5-6:31%:inFrame"
        );
        assert_eq!(
            top_events_list("Top1:CA;-25;1%;inFrame"),
            "Top1:CA:-25:1%:inFrame"
        );
        assert_eq!(
            top_events_list("Top1:CD;-330;73 314%;inFrame|Top2:CA;+112;1,5%;inFrame"),
            "Top1:CD:-330:73_314%:inFrame|Top2:CA:+112:1,5%:inFrame",
            "a `,` inside an event is text"
        );
        assert_eq!(top_events_list(""), "");
        assert_eq!(
            top_events_list("a||b"),
            "a||b",
            "an interior empty event stays"
        );
        assert_eq!(top_events_list("|a"), "|a", "a leading empty event stays");
        assert_eq!(
            top_events_list("a||"),
            "a",
            "trailing empty events are dropped"
        );
        assert_eq!(top_events_list("||"), "");
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[row(
                5030000,
                "G",
                "T",
                "ENST00000300305",
                "Top1:CD;+284;15 840%;Frameshift|Top2:ES;5-6;1,5%;inFrame",
                "0.36",
            )],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let record = variant(5030000, 5030000, b"G", b"T", &["ENST00000300305"]);
        assert_eq!(
            values(&plugin, record, "SpliceVault_top_events"),
            [Some(
                "Top1:CD:+284:15_840%:Frameshift|Top2:ES:5-6:1,5%:inFrame".to_string()
            )]
        );
    }

    /// `SpliceVault.pm:201` splits a row into twelve columns by position and
    /// `:211-219` map seven of them to the output fields; the header's names are
    /// never read, and a deletion matches through its trimmed alleles at the
    /// position after its anchor base.
    #[test]
    fn twelve_columns_map_to_seven_fields() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[[
                "chr21",
                "5030010",
                "GAC",
                "G",
                "ENST00000300305",
                "Acceptor_loss",
                "chr21:5030041",
                "0.23",
                "Frameshift:1/2",
                EVENTS,
                "102996",
                "65",
            ]
            .join("\t")],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let deletion = variant(5030011, 5030012, b"AC", b"-", &["ENST00000300305"]);
        let mut variants = vec![deletion];
        plugin.run_batch(&mut variants).unwrap();
        let fields: Vec<(&str, &str)> = variants[0].transcript_consequences[0]
            .plugin_data
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
        assert_eq!(
            fields,
            [
                ("SpliceVault_site_sample_count", "102996"),
                ("SpliceVault_site_max_depth", "65"),
                ("SpliceVault_SpliceAI_delta", "0.23"),
                ("SpliceVault_out_of_frame_events", "Frameshift:1/2"),
                (
                    "SpliceVault_top_events",
                    "Top1:CD:-35:31%:Frameshift|Top2:CD:+48:0.006%:Frameshift"
                ),
                ("SpliceVault_site_pos", "chr21:5030041"),
                ("SpliceVault_site_type", "Acceptor_loss"),
            ]
        );
        let mut written: Vec<String> = fields.iter().map(|(k, _)| k.to_string()).collect();
        written.sort();
        let mut declared: Vec<String> = plugin.header_info().into_iter().map(|(f, _)| f).collect();
        declared.sort();
        assert_eq!(written, declared);
        let mut regulatory = variant(5030011, 5030012, b"AC", b"-", &["ENST00000300305"]);
        regulatory.transcript_consequences[0].feature_type = FeatureType::RegulatoryFeature;
        assert_eq!(values(&plugin, regulatory, "SpliceVault_site_pos"), [None]);
    }

    /// Ensembl VEP 116.3 never hands SpliceVault a structural variant: the
    /// module keeps the default `variant_feature_types` of
    /// `['VariationFeature']` (BaseVepPlugin.pm:100), which the dispatch checks
    /// before running a plugin on a line (OutputFactory.pm:2229). A `<DEL>`
    /// row gets no field, and nor does a matching substitution carrying the
    /// structural mark, which pins the gate rather than the allele match; the
    /// plain substitution is annotated.
    #[test]
    fn a_structural_variant_row_gets_no_field() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[row(5030000, "G", "T", "ENST00000300305", EVENTS, "0.36")],
        );
        let plugin = loaded(&[format!("file={path}")]);
        assert!(!plugin.annotates_structural_variants());
        let field = "SpliceVault_SpliceAI_delta";
        let mut deletion = variant(5030000, 5030100, b"N", b"<DEL>", &["ENST00000300305"]);
        deletion.is_structural = true;
        assert_eq!(values(&plugin, deletion, field), [None]);
        let mut flagged = variant(5030000, 5030000, b"G", b"T", &["ENST00000300305"]);
        flagged.is_structural = true;
        assert_eq!(values(&plugin, flagged, field), [None]);
        let plain = variant(5030000, 5030000, b"G", b"T", &["ENST00000300305"]);
        assert_eq!(values(&plugin, plain, field), [Some("0.36".to_string())]);
    }

    /// `SpliceVault.pm:201` splits the line into twelve variables, so Ensembl VEP keeps
    /// a trailing empty cell as the empty string and leaves a missing cell
    /// undefined, and `:206-220` write both. Both print as
    /// `SpliceVault_site_max_depth=` in the default format; the tab format
    /// prints an empty cell against `-` and JSON `""` against `null`
    /// (`Tab.pm:128`; the nested hash is encoded as it stands), so the row names
    /// the field of a present empty cell in `plugin_empty_strings`, an interior
    /// one included, and not that of a missing one. The events column is a list
    /// either way (`:216`), which `Tab.pm:128` prints as an empty cell when it
    /// is empty, so its field is named when no event is left, the cell present
    /// or not.
    #[test]
    fn a_row_lacking_its_last_cell_is_the_empty_value_and_an_empty_cell_the_empty_string() {
        let dir = tempfile::tempdir().unwrap();
        let ending = |alt: &str, tail: &str| {
            let full = row(5030000, "G", alt, "ENST00000300305", EVENTS, "0.36");
            let (head, _) = full.rsplit_once('\t').expect("twelve columns");
            format!("{head}{tail}")
        };
        let delta_empty = row(5030000, "G", "GA", "ENST00000300305", "", "");
        let nine_cells = row(5030000, "G", "GT", "ENST00000300305", EVENTS, "0.36")
            .split('\t')
            .take(9)
            .collect::<Vec<_>>()
            .join("\t");
        let path = fixture(
            dir.path(),
            &[
                ending("A", ""),
                ending("C", "\t"),
                ending("T", "\t\tx"),
                delta_empty,
                nine_cells,
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let marks = |alt: &[u8]| {
            let mut variants = vec![variant(5030000, 5030000, b"G", alt, &["ENST00000300305"])];
            plugin.run_batch(&mut variants).unwrap();
            let tc = &variants[0].transcript_consequences[0];
            (
                tc.plugin_data.get("SpliceVault_site_max_depth").cloned(),
                tc.plugin_empty_strings.clone(),
            )
        };
        let empty = Some(String::new());
        let depth = vec!["SpliceVault_site_max_depth".to_string()];
        assert_eq!(marks(b"A"), (empty.clone(), vec![]), "eleven cells");
        assert_eq!(
            marks(b"C"),
            (empty.clone(), depth.clone()),
            "twelve cells, the twelfth empty"
        );
        assert_eq!(
            marks(b"T"),
            (empty.clone(), depth),
            "thirteen cells, the twelfth empty"
        );
        let mut variants = vec![variant(5030001, 5030000, b"-", b"A", &["ENST00000300305"])];
        plugin.run_batch(&mut variants).unwrap();
        let tc = &variants[0].transcript_consequences[0];
        assert_eq!(
            tc.plugin_data.get("SpliceVault_SpliceAI_delta"),
            Some(&String::new()),
            "an interior empty cell is written empty"
        );
        assert_eq!(
            tc.plugin_data.get("SpliceVault_top_events"),
            Some(&String::new()),
            "an empty events cell is the empty list"
        );
        assert_eq!(
            tc.plugin_empty_strings,
            vec![
                "SpliceVault_SpliceAI_delta".to_string(),
                "SpliceVault_top_events".to_string()
            ],
            "the interior empty cell and the empty events list are named"
        );
        let mut variants = vec![variant(5030001, 5030000, b"-", b"T", &["ENST00000300305"])];
        plugin.run_batch(&mut variants).unwrap();
        let tc = &variants[0].transcript_consequences[0];
        assert_eq!(
            tc.plugin_data.get("SpliceVault_site_max_depth"),
            Some(&String::new()),
            "a missing cell is written empty"
        );
        assert_eq!(
            tc.plugin_empty_strings,
            vec!["SpliceVault_top_events".to_string()],
            "a missing events cell is the empty list and is named; the other missing cells are not"
        );
    }

    /// `SpliceVault.pm:201` splits the line htslib hands over, its line end
    /// alone removed, on `\t`, and `:213` copies the last cell into
    /// `SpliceVault_site_max_depth` as it stands: `93 ` is written with its
    /// trailing space, and a cell holding a single space as that space.
    #[test]
    fn a_trailing_space_on_the_last_cell_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let ending = |alt: &str, max_depth: &str| {
            let full = row(5030000, "G", alt, "ENST00000300305", EVENTS, "0.36");
            let (head, _) = full.rsplit_once('\t').expect("twelve columns");
            format!("{head}\t{max_depth}")
        };
        let path = fixture(dir.path(), &[ending("A", "93 "), ending("T", " ")]);
        let plugin = loaded(&[format!("file={path}")]);
        let field = "SpliceVault_site_max_depth";
        let spaced = variant(5030000, 5030000, b"G", b"A", &["ENST00000300305"]);
        assert_eq!(values(&plugin, spaced, field), [Some("93 ".to_string())]);
        let blank = variant(5030000, 5030000, b"G", b"T", &["ENST00000300305"]);
        assert_eq!(values(&plugin, blank, field), [Some(" ".to_string())]);
    }

    /// Per record of a batch, per transcript row: the number of the plugin's
    /// fields written and the fields named in `plugin_empty_strings`.
    type RowFields = Vec<Vec<(usize, Vec<String>)>>;

    /// The [`RowFields`] of the batch `variants` after it runs, with the batch's
    /// result beside them.
    fn fields_of(
        plugin: &SpliceVaultPlugin,
        variants: &mut [InputVariant],
    ) -> (Result<(), PluginError>, RowFields) {
        let result = plugin.run_batch(variants);
        let fields = variants
            .iter()
            .map(|variant| {
                variant
                    .transcript_consequences
                    .iter()
                    .map(|tc| (tc.plugin_data.len(), tc.plugin_empty_strings.clone()))
                    .collect()
            })
            .collect();
        (result, fields)
    }

    /// With `transcript_match=0` every row of the window is put to the match in
    /// file order (`SpliceVault.pm:160-178`), and `get_matched_variant_alleles`
    /// dies on a row whose reference or alternate cell is undefined
    /// (ensembl-variation Utils/Sequence.pm:968). The twelve-scalar split of
    /// `:201` keeps a trailing empty cell, so a cell is undefined only when the
    /// row ends before it: a two-cell row has no reference, a three-cell row no
    /// alternate. Such a row reached before the record's first allele match, at
    /// its position or two bases before it (`:158`), fails the record: none of
    /// the seven fields on any of its transcript rows, the batch's error naming
    /// the first such row, while a record at a clean position in the same batch
    /// keeps its fields.
    #[test]
    fn transcript_match_0_fails_a_record_on_a_short_row_before_its_match() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[
                "chr21\t5029998\tG".to_string(),
                row(5030000, "G", "T", "ENST00000300305", EVENTS, "0.36"),
                "chr21\t5030010".to_string(),
                row(5030010, "G", "T", "ENST00000300305", EVENTS, "0.47"),
                row(5030020, "G", "T", "ENST00000300305", EVENTS, "0.58"),
            ],
        );
        let plugin = loaded(&[format!("file={path}"), "transcript_match=0".into()]);
        let mut variants = vec![
            variant(
                5030000,
                5030000,
                b"G",
                b"T",
                &["ENST00000300305", "ENST00000346798"],
            ),
            variant(5030010, 5030010, b"G", b"T", &["ENST00000300305"]),
            variant(5030020, 5030020, b"G", b"T", &["ENST00000300305"]),
        ];
        let (result, fields) = fields_of(&plugin, &mut variants);
        match result {
            Err(PluginError::Run(message)) => {
                assert!(
                    message.contains(":5029998 has no reference or alternate allele cell"),
                    "{message}"
                );
                assert!(message.contains(&path), "{message}");
            }
            other => panic!("expected PluginError::Run, got {other:?}"),
        }
        assert_eq!(
            fields,
            [
                vec![(0, vec![]), (0, vec![])],
                vec![(0, vec![])],
                vec![(7, vec![])],
            ],
            "the two failed records keep nothing, the clean one its seven fields"
        );
        assert_eq!(
            variants[2].transcript_consequences[0]
                .plugin_data
                .get("SpliceVault_SpliceAI_delta"),
            Some(&"0.58".to_string())
        );
    }

    /// The loop returns at the first allele match (`SpliceVault.pm:178-193`),
    /// so a short row after it is never read and fails nothing; and the
    /// twelve-scalar split of `:201` keeps an empty cell, the row's last or
    /// before a filled one, as the defined empty string, which the match trims
    /// to `-` and compares with no die (Utils/Sequence.pm:968 dies on an
    /// undefined cell alone), so a row with an empty alternate cell before the
    /// matching row matches no substitution and fails nothing.
    #[test]
    fn transcript_match_0_reads_a_short_row_after_the_match_and_an_empty_alternate_cell_as_no_failure(
    ) {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[
                "chr21\t5030000\tG\t".to_string(),
                "chr21\t5030000\tG\t\tENST00000300305".to_string(),
                row(5030000, "G", "T", "ENST00000300305", EVENTS, "0.36"),
                "chr21\t5030000\tG".to_string(),
                "chr21\t5030000".to_string(),
                "chr21\t5030001".to_string(),
            ],
        );
        let plugin = loaded(&[format!("file={path}"), "transcript_match=0".into()]);
        let record = variant(
            5030000,
            5030000,
            b"G",
            b"T",
            &["ENST00000300305", "ENST00000346798"],
        );
        assert_eq!(
            values(&plugin, record, "SpliceVault_SpliceAI_delta"),
            [Some("0.36".to_string()), Some("0.36".to_string())]
        );
        let deletion = variant(5030000, 5030000, b"G", b"-", &["ENST00000300305"]);
        let mut variants = vec![deletion];
        let (result, fields) = fields_of(&plugin, &mut variants);
        result
            .expect("the empty alternate cell is the deletion's match, read before any short row");
        assert_eq!(
            fields,
            [vec![(7, vec!["SpliceVault_top_events".to_string()])]],
            "the four-cell row's missing events cell is the empty list and is named"
        );
        assert_eq!(
            variants[0].transcript_consequences[0]
                .plugin_data
                .get("SpliceVault_SpliceAI_delta"),
            Some(&"".to_string()),
            "the four-cell row's cells past its alternate are undefined and written empty"
        );
    }

    /// Under `transcript_match` the feature test runs before the match
    /// (`SpliceVault.pm:161-163`), and a row of fewer than five cells has no
    /// feature cell to equal the transcript's stable id, so a two-, three- or
    /// four-cell row before the matching row is skipped and the match's die
    /// (Utils/Sequence.pm:968) is never reached: the transcript's own row
    /// annotates it, a transcript no row names gets nothing, and the batch
    /// succeeds either way.
    #[test]
    fn transcript_match_on_skips_a_short_row_before_the_match() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            &[
                "chr21\t5029998\tG".to_string(),
                "chr21\t5030000".to_string(),
                "chr21\t5030000\tG".to_string(),
                "chr21\t5030000\tG\t".to_string(),
                row(5030000, "G", "T", "ENST00000300305", EVENTS, "0.36"),
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let record = variant(
            5030000,
            5030000,
            b"G",
            b"T",
            &["ENST00000300305", "ENST00000999999"],
        );
        assert_eq!(
            values(&plugin, record, "SpliceVault_SpliceAI_delta"),
            [Some("0.36".to_string()), None]
        );
    }

    /// A row's `pos` cell is keyed by `get_matched_variant_alleles`
    /// (`SpliceVault.pm:165-178`) as text when the alleles lose no leading base
    /// (ensembl-variation Utils/Sequence.pm:1232), and as its number plus the
    /// bases moved once a shared leading base leaves them (`:971, :1013`): for a
    /// substitution `5030030` annotates and `5030030.0` and `5030030.5` do not,
    /// nor does `5030031x` at 5030031; for the deletion `GGC>G` at `5030040.0`
    /// and `5030060x` the moved start is the variant's, while `5030050.5` moves
    /// to 5030051.5 and matches nothing.
    #[test]
    fn a_position_cell_is_matched_as_the_key_holds_it() {
        let dir = tempfile::tempdir().unwrap();
        let cell_row = |pos: &str, reference: &str, alt: &str, delta: &str| {
            row(1, reference, alt, "ENST00000400001", EVENTS, delta).replacen(
                "\t1\t",
                &format!("\t{pos}\t"),
                1,
            )
        };
        let path = fixture(
            dir.path(),
            &[
                cell_row("5030030", "G", "A", "0.11"),
                cell_row("5030030.0", "G", "C", "0.12"),
                cell_row("5030030.5", "G", "T", "0.13"),
                cell_row("5030031x", "G", "A", "0.14"),
                cell_row("5030040.0", "GGC", "G", "0.21"),
                cell_row("5030050.5", "GGC", "G", "0.22"),
                cell_row("5030060x", "GGC", "G", "0.23"),
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let delta = |start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            values(
                &plugin,
                variant(start, end, reference, alt, &["ENST00000400001"]),
                "SpliceVault_SpliceAI_delta",
            )
        };
        assert_eq!(
            delta(5030030, 5030030, b"G", b"A"),
            [Some("0.11".to_string())]
        );
        assert_eq!(delta(5030030, 5030030, b"G", b"C"), [None]);
        assert_eq!(delta(5030030, 5030030, b"G", b"T"), [None]);
        assert_eq!(delta(5030031, 5030031, b"G", b"A"), [None]);
        assert_eq!(
            delta(5030041, 5030042, b"GC", b"-"),
            [Some("0.21".to_string())]
        );
        assert_eq!(delta(5030051, 5030052, b"GC", b"-"), [None]);
        assert_eq!(
            delta(5030061, 5030062, b"GC", b"-"),
            [Some("0.23".to_string())]
        );

        // A made-up row short of its allele cells fails the variant, unless its position
        // cell is Perl-false, which returns the match first (Sequence.pm:1173).
        let short = |pos: &str| TabixRecord {
            chr: "21".to_string(),
            start: 5030030,
            end: 5030030,
            columns: vec!["chr21".into(), pos.into(), "G".into()],
            header: None,
        };
        let snv = variant(5030030, 5030030, b"G", b"A", &["ENST00000400001"]);
        let matched = |pos: &str| {
            let record = short(pos);
            SpliceVaultPlugin::allele_matches(&[("short_synthetic.tsv.gz", vec![&record])], &snv)[0]
                .matched
        };
        assert_eq!(matched("5030030"), None);
        assert_eq!(matched("0"), Some(false));
        assert_eq!(matched(""), Some(false));
    }
}
