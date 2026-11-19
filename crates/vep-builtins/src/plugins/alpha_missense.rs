// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! AlphaMissense plugin: the pre-computed AlphaMissense pathogenicity
//! predictions of missense rows, read from the tabix-indexed prediction table
//! and written as Ensembl VEP's `AlphaMissense.pm` writes them.
//!
//! # Parameters
//!
//! Every parameter is `key=value` (`AlphaMissense.pm:175`); a file given by
//! position is refused as the file given without its key.
//!
//! - `file=<path>`: the prediction table (mandatory, `:179-181`).
//! - `cols=<name>[:<name>...]` or `cols=all`: the columns written, by their
//!   output names or the file's own, each prefixed `am_` unless it already is
//!   (`:146-164`); a name the header lacks is refused with the valid names
//!   listed. Default `am_pathogenicity:am_class` (`:185`).
//! - `transcript_match=<value>`: write a prediction only when its transcript is
//!   the row's (`:256-261`), read as Ensembl VEP reads a Perl value (every value
//!   but the empty string and `0` sets it). Default off (`:176`).
//!
//! # Data file
//!
//! The header's last line names the columns; the first four identify the
//! variant (`CHROM POS REF ALT`) and the rest are the prediction's columns
//! (`:117-137`): `genome uniprot_id transcript_id protein_variant
//! am_pathogenicity am_class` in this module's fixtures. The reader reads the
//! build in the `genome` column, and a file whose rows name the other build is
//! refused under `--assembly`, since the two builds' files share one layout and
//! the wrong one annotates nothing at exit 0.
//!
//! # Matching
//!
//! The records are read from two bases before the variant's start to its end
//! (`:248`), the window clamped at position 1. A row is annotated when it is a
//! missense row (`:242`) and a record in that window carries the variant's
//! alleles at its position once both are trimmed of the bases they share
//! (`:263-275`, `get_matched_variant_alleles`) and a `protein_variant` whose
//! first and last letters are the row's reference and alternate peptides
//! (`:221-236`); under `transcript_match` the record's `transcript_id`, its
//! version suffix removed, must be the row's transcript as well (`:256-261`).
//! The first matching record's selected columns are written as the file spells
//! them. `AlphaMissense.pm:259` removes a version suffix with `s/\.[0-9]//g`,
//! which turns a two-digit version (`ENST00000000233.10`) into another id
//! (`ENST0000000023310`) and so never matches it; vep-rs removes the whole
//! suffix, so such a record matches its transcript.
//!
//! # Output fields
//!
//! The selected columns (`:199-219`), each described as `column from <file>`
//! with the five published columns carrying their own text. In the JSON format
//! the fields nest under `alphamissense` (`:281-282`), which the JSON writer
//! does for every `am_`-prefixed plugin field.
//!
//! `AlphaMissense.pm` is the plugin's Ensembl VEP counterpart in
//! Ensembl/VEP_plugins release/116.

use std::path::Path;

use tracing::debug;
use vep_core::assembly::Assembly;
use vep_core::consequence::{Consequence, TranscriptConsequence};
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

use crate::annotation_store::{check_table_or_store, open_best_store_with, AnnotationStore};
use crate::params::{Grammar, PluginParams};
use crate::tabix::{
    allele_match, parse_header, split_field_count, BatchQueryResult, TabixAnnotator, TabixConfig,
    TabixRecord,
};
use crate::traits::{warn_on_stderr, BuiltinPlugin, PluginError, PrefetchData};

const NAME: &str = "AlphaMissense";

/// The columns identifying the variant, before the prediction's columns
/// (`AlphaMissense.pm:132`, `:290`).
const VARIANT_COLUMNS: usize = 4;
const CHR_COL: usize = 0;
const POS_COL: usize = 1;
const REF_COL: usize = 2;
const ALT_COL: usize = 3;

/// The prefix every output name carries (`AlphaMissense.pm:135`, `:152`).
const FIELD_PREFIX: &str = "am_";

/// The bases before a variant's start its rows are read from
/// (`AlphaMissense.pm:248`, `$vf->{start} - 2`).
const WINDOW_BEFORE: u64 = 2;

/// The columns written when `cols` is not given (`AlphaMissense.pm:185`).
const DEFAULT_COLS: &str = "am_pathogenicity:am_class";

const PROTEIN_VARIANT: &str = "am_protein_variant";
const TRANSCRIPT_ID: &str = "am_transcript_id";
const GENOME: &str = "am_genome";

/// AlphaMissense pathogenicity annotation plugin.
pub struct AlphaMissensePlugin {
    /// Annotation store, initialized in `init()`.
    store: Option<Box<dyn AnnotationStore>>,
    /// Target assembly from `--assembly` or the cache, injected before `init()`.
    assembly: Option<Assembly>,
    /// The file, as given, named in every field description.
    file: String,
    /// The output name of every prediction column, in file order from the fifth
    /// column (`AlphaMissense.pm:117-137`).
    colnames: Vec<String>,
    /// The selected columns in the order given: output name and 0-based index
    /// in a record's columns.
    selected: Vec<(String, usize)>,
    /// The 0-based column of the amino-acid change and of the transcript id.
    protein_variant_col: Option<usize>,
    transcript_col: Option<usize>,
    /// Write a prediction only when its transcript is the row's.
    transcript_match: bool,
    /// Whether the run reads the RefSeq or merged transcript set.
    other_transcript_set: bool,
}

impl AlphaMissensePlugin {
    /// Create a new uninitialized AlphaMissense plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            assembly: None,
            file: String::new(),
            colnames: Vec::new(),
            selected: DEFAULT_COLS
                .split(':')
                .enumerate()
                .map(|(index, name)| (name.to_string(), VARIANT_COLUMNS + index))
                .collect(),
            protein_variant_col: None,
            transcript_col: None,
            transcript_match: false,
            other_transcript_set: false,
        }
    }

    /// Whether a transcript consequence includes missense_variant
    /// (`AlphaMissense.pm:242`).
    fn is_missense(consequence: &TranscriptConsequence) -> bool {
        consequence
            .consequences
            .contains(&Consequence::MissenseVariant)
    }

    /// A column name with the `am_` prefix, added unless already there
    /// (`AlphaMissense.pm:135`, `:152`).
    fn prefixed(name: &str) -> String {
        if name.starts_with(FIELD_PREFIX) {
            name.to_string()
        } else {
            format!("{FIELD_PREFIX}{name}")
        }
    }

    /// A stable id without its `.<version>` suffix.
    fn unversioned(id: &str) -> &str {
        id.split('.').next().unwrap_or(id)
    }

    /// The 0-based column, in a record, of the prediction column whose output
    /// name is `name`; `None` when the header lacks it.
    fn column(&self, name: &str) -> Option<usize> {
        self.colnames
            .iter()
            .position(|colname| colname == name)
            .map(|index| VARIANT_COLUMNS + index)
    }

    /// The columns `cols` selects (`AlphaMissense.pm:139-165`): every column
    /// under `all`, else the `:`-separated names, prefixed, each of which the
    /// header must carry, a name listed twice selected once (`:216` and `:279`
    /// key the header and the row by a hash over the list). Trailing empty names
    /// are dropped before the prefix is added, as Perl's `split` drops them
    /// (`:149`), so a literal `am_` is a name the header must carry.
    fn select_columns(&self, cols: &str) -> Result<Vec<(String, usize)>, PluginError> {
        if cols == "all" {
            return Ok(self
                .colnames
                .iter()
                .enumerate()
                .map(|(index, name)| (name.clone(), VARIANT_COLUMNS + index))
                .collect());
        }
        let mut names: Vec<&str> = cols.split(':').collect();
        while names.last() == Some(&"") {
            names.pop();
        }
        let names: Vec<String> = names.into_iter().map(Self::prefixed).collect();
        let mut selected: Vec<(String, usize)> = Vec::with_capacity(names.len());
        let mut invalid = Vec::new();
        for name in names {
            match self.column(&name) {
                Some(index) if selected.iter().any(|(_, column)| *column == index) => {}
                Some(index) => selected.push((name, index)),
                None => invalid.push(name),
            }
        }
        if !invalid.is_empty() {
            return Err(PluginError::Init(format!(
                "{NAME}: the following columns were not found in file header: {}; valid columns \
                 are: {}",
                invalid.join(", "),
                self.colnames.join(", ")
            )));
        }
        Ok(selected)
    }

    /// The description of an output column (`AlphaMissense.pm:204-213`).
    fn description(&self, name: &str) -> String {
        let suffix = format!("column from {}", self.file);
        let text = match name {
            "am_pathogenicity" => {
                "Continuous AlphaMissense score between 0 and 1 which can be interpreted as the \
                 predicted probability of the variant being pathogenic"
            }
            "am_class" => {
                "The AlphaMissense thresholds are: 'Likely benign' if score < 0.34, 'Likely \
                 pathogenic' if score > 0.564, 'ambiguous' otherwise -- see \
                 doi.org/10.1126/science.adg7492 for details"
            }
            PROTEIN_VARIANT => "Amino acid change used in AlphaMissense prediction",
            "am_uniprot_id" => "Protein isoform used in AlphaMissense prediction",
            TRANSCRIPT_ID => "Transcript sequence in AlphaMissense prediction",
            _ => return suffix,
        };
        format!("{text}; {suffix}")
    }

    /// Whether `record`'s amino-acid change is the row's (`AlphaMissense.pm:221-236`):
    /// the first and last letters of its `protein_variant` are the row's reference
    /// and alternate peptides. A file without that column matches nothing.
    fn amino_acid_changes_match(
        &self,
        record: &TabixRecord,
        reference_peptide: &str,
        alternate_peptide: &str,
    ) -> bool {
        let Some(protein_variant) = self
            .protein_variant_col
            .and_then(|column| record.columns.get(column))
        else {
            return false;
        };
        let mut letters = protein_variant.chars();
        let (Some(first), Some(last)) = (letters.next(), protein_variant.chars().next_back())
        else {
            return false;
        };
        let mut reference = reference_peptide.chars();
        let mut alternate = alternate_peptide.chars();
        reference.next() == Some(first)
            && reference.next().is_none()
            && alternate.next() == Some(last)
            && alternate.next().is_none()
    }

    /// Whether `record`'s transcript is the row's, when `transcript_match` asks
    /// (`AlphaMissense.pm:256-261`), the version suffix set aside on both.
    fn transcript_matches(&self, record: &TabixRecord, transcript_id: &str) -> bool {
        if !self.transcript_match {
            return true;
        }
        self.transcript_col
            .and_then(|column| record.columns.get(column))
            .is_some_and(|id| Self::unversioned(id) == Self::unversioned(transcript_id))
    }

    /// The first of `records` (the records at the variant carrying its alleles,
    /// in file order) whose amino-acid change and, when matched, transcript are
    /// the row's (`AlphaMissense.pm:250-284`).
    fn matching_record<'a>(
        &self,
        records: &[&'a TabixRecord],
        row: &TranscriptConsequence,
    ) -> Option<&'a TabixRecord> {
        let (reference_peptide, alternate_peptide) = row.peptides()?;
        records.iter().copied().find(|record| {
            self.amino_acid_changes_match(record, reference_peptide, alternate_peptide)
                && self.transcript_matches(record, row.transcript_id())
        })
    }

    /// The region a variant's rows are read from (`AlphaMissense.pm:248`,
    /// `get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})`): two bases before
    /// its start to its end, so a multi-base row written against an earlier
    /// anchor base is read and the match trims it to the variant's position.
    /// The start is clamped at position 1, the rows a record at the first
    /// bases of a chromosome has (the documented class
    /// `tabix_window_below_position_one`).
    fn window(variant: &InputVariant) -> (u64, u64) {
        (
            variant.start.saturating_sub(WINDOW_BEFORE).max(1),
            variant.end,
        )
    }

    /// The records of `batch` in `variant`'s window carrying its alleles, in
    /// file order (`get_matched_variant_alleles`, `AlphaMissense.pm:263-275`:
    /// a row is the variant's when its alleles, trimmed of the bases they
    /// share, are the variant's own at the variant's position), each row first
    /// adjusted as `parse_data` adjusts an indel-shaped one (`:293-300`,
    /// [`allele_match::vcf_adjusted_row`]): the shared first base dropped and
    /// the position cell stepped as a string.
    fn matching_records<'a>(
        &self,
        batch: &'a BatchQueryResult,
        variant: &InputVariant,
    ) -> Vec<&'a TabixRecord> {
        let (start, end) = Self::window(variant);
        batch
            .lookup(&variant.chr, start, end)
            .into_iter()
            .filter(|record| {
                allele_match::matched_variant_alleles(
                    &allele_match::vcf_adjusted_row(record, POS_COL, REF_COL, ALT_COL),
                    variant,
                    POS_COL,
                    REF_COL,
                    ALT_COL,
                )
            })
            .collect()
    }

    /// Refuse a file whose rows name the other build. The two upstream files are
    /// pre-split per build, identical in shape and both UCSC-contig, so a
    /// wrong-build file resolves contigs and yields zero hits at exit 0; the
    /// `genome` column is the only in-file build signal.
    fn check_build(&self, path: &Path) -> Result<(), PluginError> {
        let (Some(requested), Some(genome_col)) = (self.assembly, self.column(GENOME)) else {
            return Ok(());
        };
        let Some(declared) = crate::tabix::first_data_row_field(path, genome_col)? else {
            return Ok(());
        };
        let declared = declared.trim();
        if !declared.is_empty()
            && !declared.eq_ignore_ascii_case(requested.ucsc_name())
            && !declared.eq_ignore_ascii_case(requested.as_str())
        {
            return Err(PluginError::Init(format!(
                "{NAME}: --assembly is {requested} but '{}' declares genome '{declared}'. \
                 AlphaMissense publishes one file per build; point --plugin {NAME} at the \
                 {requested} file.",
                self.file
            )));
        }
        Ok(())
    }
}

impl Default for AlphaMissensePlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl BuiltinPlugin for AlphaMissensePlugin {
    fn name(&self) -> &str {
        NAME
    }

    fn header_info(&self) -> Vec<(String, String)> {
        self.selected
            .iter()
            .map(|(name, _)| (name.clone(), self.description(name)))
            .collect()
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse_with(
            NAME,
            Grammar::KeyValueOnly { file_key: "file" },
            &["file", "cols", "transcript_match"],
            params,
        )?;
        let file = params
            .get("file")
            .filter(|file| !file.is_empty())
            .ok_or_else(|| {
                PluginError::Init(format!("{NAME}: no file specified; use file=<path>"))
            })?;
        self.transcript_match = params.get("transcript_match").is_some_and(perl_true);
        self.file = file.to_string();
        let path = Path::new(file);
        // `add_file` runs before the header is read (AlphaMissense.pm:182).
        check_table_or_store(path)?;

        // The last header line names the columns; the prediction's columns follow
        // the four that identify the variant (`AlphaMissense.pm:117-137`).
        self.colnames = parse_header(path, POS_COL)?
            .unwrap_or_default()
            .iter()
            .skip(VARIANT_COLUMNS)
            .map(|name| Self::prefixed(name))
            .collect();
        self.protein_variant_col = self.column(PROTEIN_VARIANT);
        self.transcript_col = self.column(TRANSCRIPT_ID);
        let cols = params
            .get("cols")
            .filter(|cols| perl_true(cols))
            .unwrap_or(DEFAULT_COLS);
        self.selected = self.select_columns(cols)?;
        self.check_build(path)?;

        let config = TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: CHR_COL,
            start_col: POS_COL,
            end_col: None,
            ref_col: Some(REF_COL),
            alt_col: Some(ALT_COL),
            zero_based: false,
        };
        // Opened on no assembly, so a mitochondrial record reads the table under
        // every spelling on GRCh37 as on GRCh38: the published GRCh37 table spells
        // the contig `chrM` and states each row at the assembly's `MT` position.
        self.store = Some(open_best_store_with(
            path,
            config,
            None,
            TabixAnnotator::keeping_trailing_whitespace,
        )?);
        // After the index check, as `AlphaMissense.pm:182-190` orders them.
        if self.other_transcript_set && self.transcript_match {
            warn_on_stderr(
                "AlphaMissense transcript_match may fail with --refseq/--merged; consider \
                 transcript_match=0",
            );
        }
        debug!(
            file,
            cols,
            transcript_match = self.transcript_match,
            "AlphaMissense plugin initialized"
        );
        Ok(())
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn set_transcript_set(&mut self, refseq: bool, merged: bool) {
        self.other_transcript_set = refseq || merged;
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        if variants.is_empty() {
            return Ok(None);
        }
        // Every sequence variant's window is queried, missense row or not: the
        // windows of a batch merge into a few contiguous regions the file is
        // read through once, where the missense windows alone form many small
        // regions, each a seek and a block inflated on its own. A structural
        // variant's alleles are class words no single-base row matches, so its
        // span, which can cover megabases of rows, is left out.
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|v| !v.is_structural)
            .map(|v| {
                let (start, end) = Self::window(v);
                (v.chr.as_str(), start, end)
            })
            .collect();
        if regions.is_empty() {
            return Ok(Some(PrefetchData::new(BatchQueryResult::empty())));
        }
        Ok(Some(PrefetchData::new(store.query_batch(&regions)?)))
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
        for variant in variants.iter_mut().filter(|v| !v.is_structural) {
            // The records come first: most variants of a genome-wide input lie
            // outside the file's coding positions and are done at the lookup,
            // before any of their rows is read.
            let records = self.matching_records(&batch, variant);
            if records.is_empty() {
                continue;
            }
            for row in &mut variant.transcript_consequences {
                if !self.annotates_row(row) || !Self::is_missense(row) {
                    continue;
                }
                let Some(record) = self.matching_record(&records, row) else {
                    continue;
                };
                // Every selected column is a key of the row's result
                // (AlphaMissense.pm:279): the `split /\t/` of `:290` drops the
                // trailing empty cells, so a cell past the row's last filled
                // one is the undefined value, the empty value here, and an
                // empty cell before a filled one the empty string, named in
                // the row's `plugin_empty_strings`.
                let present = split_field_count(&record.columns);
                for (name, column) in &self.selected {
                    let cell = record.columns.get(*column).filter(|_| *column < present);
                    row.plugin_data
                        .insert(name.clone(), cell.cloned().unwrap_or_default());
                    row.plugin_empty_strings.retain(|key| key != name);
                    if cell.is_some_and(String::is_empty) {
                        row.plugin_empty_strings.push(name.clone());
                    }
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use smallvec::smallvec;
    use std::collections::{BTreeMap, HashMap};
    use std::sync::Arc;
    use vep_core::consequence::FeatureType;
    use vep_core::transcript::Transcript;

    const HEADER: &str = "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class";

    /// A made-up fixture: one made-up comment line, the header line of the ten
    /// columns, then `rows`.
    fn fixture(dir: &Path, name: &str, rows: &[&str]) -> String {
        let mut lines = vec!["# A made-up comment line of the fixture", HEADER];
        lines.extend_from_slice(rows);
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
            .display()
            .to_string()
    }

    fn loaded<S: AsRef<str>>(params: &[S]) -> AlphaMissensePlugin {
        let mut plugin = AlphaMissensePlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        let params: Vec<String> = params.iter().map(|p| p.as_ref().to_string()).collect();
        plugin.init(&params).expect("fixture loads");
        plugin
    }

    fn init_message<S: AsRef<str>>(assembly: Option<Assembly>, params: &[S]) -> String {
        let mut plugin = AlphaMissensePlugin::new();
        plugin.set_assembly(assembly);
        let params: Vec<String> = params.iter().map(|p| p.as_ref().to_string()).collect();
        match plugin.init(&params) {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected PluginError::Init, got {other:?}"),
        }
    }

    /// A row of `transcript_id` with the consequence and amino-acid change given.
    fn row(
        transcript_id: &str,
        consequence: Consequence,
        amino_acids: &str,
    ) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: transcript_id.into(),
                ..Transcript::default()
            })),
            consequences: smallvec![consequence],
            amino_acids: Some(amino_acids.to_string()),
            ..Default::default()
        }
    }

    fn snv(pos: u64, reference: &[u8], alt: &[u8]) -> InputVariant {
        InputVariant::new("21".into(), pos, pos, reference.to_vec(), alt.to_vec())
    }

    /// The plugin fields of every row of `variant` after a batch run, as
    /// `key=value` lists in the plugin's order.
    fn fields(plugin: &AlphaMissensePlugin, variant: InputVariant) -> Vec<Vec<String>> {
        let mut variants = vec![variant];
        plugin.run_batch(&mut variants).expect("run_batch");
        variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| {
                tc.plugin_data
                    .iter()
                    .map(|(key, value)| format!("{key}={value}"))
                    .collect()
            })
            .collect()
    }

    /// `add_file` checks for `<file>.tbi` before the header is read
    /// (AlphaMissense.pm:182): a table lacking both is refused naming the index.
    #[test]
    fn table_without_its_tbi_is_refused_before_its_header_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "am.tsv.gz",
            &["#CHROM\tPOS", "21\t100"],
            IndexSpec::standard(),
        );
        std::fs::remove_file(dir.path().join("am.tsv.gz.tbi")).unwrap();
        let message = init_message(
            Some(Assembly::Grch38),
            &[format!("file={}", path.display())],
        );
        assert!(message.contains("am.tsv.gz.tbi not found"), "{message}");
    }

    #[test]
    fn plugin_name_is_correct() {
        assert_eq!(AlphaMissensePlugin::new().name(), "AlphaMissense");
    }

    /// `cols=all` writes every column after the four variant columns, each
    /// under its `am_`-prefixed name in file order (`AlphaMissense.pm:132-135`,
    /// `:146-147`), and the header describes each with its file
    /// (`:199-219`).
    #[test]
    fn cols_all_emits_every_am_column() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &["chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA2B\t0.9123\tlikely_pathogenic"],
        );
        let plugin = loaded(&[format!("file={path}"), "cols=all".to_string()]);
        assert_eq!(
            plugin.header_info(),
            vec![
                ("am_genome".to_string(), format!("column from {path}")),
                (
                    "am_uniprot_id".to_string(),
                    format!("Protein isoform used in AlphaMissense prediction; column from {path}")
                ),
                (
                    "am_transcript_id".to_string(),
                    format!("Transcript sequence in AlphaMissense prediction; column from {path}")
                ),
                (
                    "am_protein_variant".to_string(),
                    format!(
                        "Amino acid change used in AlphaMissense prediction; column from {path}"
                    )
                ),
                (
                    "am_pathogenicity".to_string(),
                    format!(
                        "Continuous AlphaMissense score between 0 and 1 which can be interpreted \
                         as the predicted probability of the variant being pathogenic; column \
                         from {path}"
                    )
                ),
                (
                    "am_class".to_string(),
                    format!(
                        "The AlphaMissense thresholds are: 'Likely benign' if score < 0.34, \
                         'Likely pathogenic' if score > 0.564, 'ambiguous' otherwise -- see \
                         doi.org/10.1126/science.adg7492 for details; column from {path}"
                    )
                ),
            ]
        );

        let mut v = snv(9068417, b"C", b"A");
        v.transcript_consequences
            .push(row("ENST00000400001", Consequence::MissenseVariant, "A/B"));
        assert_eq!(
            fields(&plugin, v),
            [[
                "am_genome=hg38",
                "am_uniprot_id=Q9NQ38",
                "am_transcript_id=ENST00000400001.4",
                "am_protein_variant=A2B",
                "am_pathogenicity=0.9123",
                "am_class=likely_pathogenic",
            ]]
        );

        // The default selection, and a selection by the file's own names, which
        // take the prefix (`:152`); a trailing `:` adds no column.
        assert_eq!(
            loaded(&[format!("file={path}")])
                .header_info()
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            ["am_pathogenicity", "am_class"]
        );
        assert_eq!(
            loaded(&[format!("file={path}"), "cols=class:uniprot_id:".to_string()])
                .header_info()
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            ["am_class", "am_uniprot_id"]
        );
    }

    /// `AlphaMissense.pm:216` and `:279` key the header and the row's values by
    /// a hash over `cols`, so a column named twice is one key: one header line
    /// in its first position and one cell.
    #[test]
    fn a_column_named_twice_in_cols_is_declared_and_written_once() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &["chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA2B\t0.9123\tlikely_pathogenic"],
        );
        let plugin = loaded(&[
            format!("file={path}"),
            "cols=am_class:am_class:pathogenicity".to_string(),
        ]);
        assert_eq!(
            plugin
                .header_info()
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            ["am_class", "am_pathogenicity"]
        );
        let mut v = snv(9068417, b"C", b"A");
        v.transcript_consequences
            .push(row("ENST00000400001", Consequence::MissenseVariant, "A/B"));
        assert_eq!(
            fields(&plugin, v),
            [["am_class=likely_pathogenic", "am_pathogenicity=0.9123"]]
        );
    }

    /// A column the header lacks is refused with the valid names listed
    /// (`AlphaMissense.pm:155-163`).
    #[test]
    fn cols_unknown_column_is_an_init_error_listing_valid_ones() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &["chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA2B\t0.9123\tlikely_pathogenic"],
        );
        assert_eq!(
            init_message(
                Some(Assembly::Grch38),
                &[
                    format!("file={path}"),
                    "cols=am_class:score:am_pathogenicity:pathogenic".to_string()
                ],
            ),
            "AlphaMissense: the following columns were not found in file header: am_score, \
             am_pathogenic; valid columns are: am_genome, am_uniprot_id, am_transcript_id, \
             am_protein_variant, am_pathogenicity, am_class"
        );
    }

    /// `AlphaMissense.pm:149` splits `cols` before `:152` adds the prefix, so
    /// only a trailing empty name is dropped: a literal `am_` is looked up in
    /// the header and refused by name, in its place among the invalid columns.
    #[test]
    fn a_literal_bare_prefix_in_cols_is_refused_as_a_column_the_header_lacks() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &["chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA2B\t0.9123\tlikely_pathogenic"],
        );
        for cols in [
            "cols=am_class:am_",
            "cols=am_class:am_:",
            "cols=am_class::am_pathogenicity",
        ] {
            assert_eq!(
                init_message(
                    Some(Assembly::Grch38),
                    &[format!("file={path}"), cols.to_string()],
                ),
                "AlphaMissense: the following columns were not found in file header: am_; valid \
                 columns are: am_genome, am_uniprot_id, am_transcript_id, am_protein_variant, \
                 am_pathogenicity, am_class",
                "{cols}"
            );
        }
        assert_eq!(
            loaded(&[format!("file={path}"), "cols=am_class::".to_string()])
                .header_info()
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            ["am_class"],
            "every trailing empty name is dropped",
        );
    }

    /// Under `transcript_match=1` a record is written on the row of its own
    /// transcript alone (`AlphaMissense.pm:256-261`). The version suffix is
    /// removed whole, so `ENST00000000233.10` is the transcript
    /// `ENST00000000233`; `AlphaMissense.pm:259`'s `s/\.[0-9]//g` leaves
    /// `ENST0000000023310` and never matches a two-digit version.
    #[test]
    fn transcript_match_strips_full_version() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &[
                "chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST00000000233.10\tA2B\t0.9123\tlikely_pathogenic",
                "chr21\t9068417\tC\tA\thg38\tQ9NQ39\tENST00000000234.2\tA2B\t0.1234\tlikely_benign",
            ],
        );
        let rows = || {
            let mut v = snv(9068417, b"C", b"A");
            v.transcript_consequences.push(row(
                "ENST00000000233",
                Consequence::MissenseVariant,
                "A/B",
            ));
            v.transcript_consequences.push(row(
                "ENST00000000234",
                Consequence::MissenseVariant,
                "A/B",
            ));
            v.transcript_consequences.push(row(
                "ENST00000000235",
                Consequence::MissenseVariant,
                "A/B",
            ));
            v
        };
        let matched = loaded(&[format!("file={path}"), "transcript_match=1".to_string()]);
        assert_eq!(
            fields(&matched, rows()),
            [
                vec!["am_pathogenicity=0.9123", "am_class=likely_pathogenic"],
                vec!["am_pathogenicity=0.1234", "am_class=likely_benign"],
                vec![],
            ]
        );
        // Without it every row with the amino-acid change takes the first
        // record (`:250-253`); `transcript_match=0` is the default (`:176`).
        for params in [
            vec![format!("file={path}")],
            vec![format!("file={path}"), "transcript_match=0".to_string()],
        ] {
            let first = vec!["am_pathogenicity=0.9123", "am_class=likely_pathogenic"];
            assert_eq!(
                fields(&loaded(&params), rows()),
                [first.clone(), first.clone(), first]
            );
        }
    }

    /// A record is written on a missense row whose reference and alternate
    /// peptides are the first and last letters of its `protein_variant`
    /// (`AlphaMissense.pm:224-235`, `:242`): a synonymous row, a missense row of
    /// another change and a regulatory row take nothing, and the alleles must
    /// be the record's (`:263-275`).
    #[test]
    fn only_missense_rows_with_matching_peptide_change_are_annotated() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &[
                "chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tQ12H\t0.9123\tlikely_pathogenic",
                "chr21\t9068417\tC\tG\thg38\tQ9NQ38\tENST00000400001.4\tQ12P\t0.4000\tambiguous",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);

        let mut v = snv(9068417, b"C", b"A");
        v.transcript_consequences.extend([
            row("ENST00000400001", Consequence::MissenseVariant, "Q/H"),
            row("ENST00000400002", Consequence::MissenseVariant, "Q/H"),
            row("ENST00000400003", Consequence::MissenseVariant, "Q/L"),
            row("ENST00000400004", Consequence::SynonymousVariant, "Q"),
            row("ENST00000400005", Consequence::MissenseVariant, "QH/HH"),
            TranscriptConsequence {
                feature_type: FeatureType::RegulatoryFeature,
                consequences: smallvec![Consequence::RegulatoryRegionVariant],
                ..TranscriptConsequence::default()
            },
        ]);
        let scored = vec!["am_pathogenicity=0.9123", "am_class=likely_pathogenic"];
        assert_eq!(
            fields(&plugin, v),
            [scored.clone(), scored, vec![], vec![], vec![], vec![]]
        );

        // The other alternate allele takes its own record.
        let mut v = snv(9068417, b"C", b"G");
        v.transcript_consequences
            .push(row("ENST00000400001", Consequence::MissenseVariant, "Q/P"));
        assert_eq!(
            fields(&plugin, v),
            [["am_pathogenicity=0.4000", "am_class=ambiguous"]]
        );

        // An allele the file lacks takes nothing, whatever its peptides.
        let mut v = snv(9068417, b"C", b"T");
        v.transcript_consequences
            .push(row("ENST00000400001", Consequence::MissenseVariant, "Q/H"));
        assert_eq!(fields(&plugin, v), [Vec::<String>::new()]);
    }

    /// A row without `missense_variant` takes nothing (`AlphaMissense.pm:242`),
    /// even from a record carrying its alleles and both letters of its
    /// amino-acid change: a start-lost row of `M/V`, which Ensembl VEP never
    /// calls missense (VariationEffect.pm:1090), beside the `M54V` record that
    /// the missense `M/V` row of another transcript takes.
    #[test]
    fn a_start_lost_row_takes_nothing_from_a_record_of_its_amino_acid_change() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &["chr21\t9068417\tA\tG\thg38\tQ9NQ38\tENST00000400001.4\tM54V\t0.9123\tlikely_pathogenic"],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let mut v = snv(9068417, b"A", b"G");
        v.transcript_consequences.extend([
            row("ENST00000400001", Consequence::MissenseVariant, "M/V"),
            row("ENST00000400002", Consequence::StartLost, "M/V"),
        ]);
        assert_eq!(
            fields(&plugin, v),
            [
                vec!["am_pathogenicity=0.9123", "am_class=likely_pathogenic"],
                vec![],
            ]
        );
    }

    /// A missense row whose reference peptide is not the first letter of the
    /// record's `protein_variant` takes nothing (`AlphaMissense.pm:234`), its
    /// alternate peptide being the record's last letter: an `A/R` row beside
    /// the `K11R` record that the `K/R` row takes.
    #[test]
    fn a_missense_row_of_another_reference_amino_acid_takes_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &["chr21\t9068417\tA\tG\thg38\tQ9NQ38\tENST00000400001.4\tK11R\t0.9123\tlikely_pathogenic"],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let mut v = snv(9068417, b"A", b"G");
        v.transcript_consequences.extend([
            row("ENST00000400001", Consequence::MissenseVariant, "K/R"),
            row("ENST00000400002", Consequence::MissenseVariant, "A/R"),
        ]);
        assert_eq!(
            fields(&plugin, v),
            [
                vec!["am_pathogenicity=0.9123", "am_class=likely_pathogenic"],
                vec![],
            ]
        );
    }

    /// Each allele of a multi-allelic record is matched against the table alone
    /// (`AlphaMissense.pm:263-275` on this allele, where `:245` hands Ensembl
    /// VEP every alternate allele of the record): two alleles of one codon
    /// that give one amino acid each take their own row, whichever allele the
    /// batch runs first and whichever row the table writes first.
    #[test]
    fn each_allele_of_a_multi_allelic_record_takes_its_own_row() {
        let dir = tempfile::tempdir().unwrap();
        let c_row =
            "chr21\t13610388\tA\tC\thg38\tQ86YR6\tENST00000299443.5\tM54L\t0.1079\tlikely_benign";
        let t_row = "chr21\t13610388\tA\tT\thg38\tQ86YR6\tENST00000299443.5\tM54L\t0.9000\tlikely_pathogenic";
        let allele = |alt: &str| {
            let mut v = snv(13610388, b"A", alt.as_bytes());
            v.transcript_consequences.push(row(
                "ENST00000299443",
                Consequence::MissenseVariant,
                "M/L",
            ));
            v
        };
        let own_prediction = |alt: &str| {
            if alt == "C" {
                ["am_pathogenicity=0.1079", "am_class=likely_benign"]
            } else {
                ["am_pathogenicity=0.9000", "am_class=likely_pathogenic"]
            }
        };
        for (name, rows) in [
            ("c_first.tsv.gz", [c_row, t_row]),
            ("t_first.tsv.gz", [t_row, c_row]),
        ] {
            let plugin = loaded(&[format!("file={}", fixture(dir.path(), name, &rows))]);
            for order in [["C", "T"], ["T", "C"]] {
                let mut variants = order.map(allele);
                plugin.run_batch(&mut variants).expect("run_batch");
                for (variant, alt) in variants.iter().zip(order) {
                    let written: Vec<String> = variant.transcript_consequences[0]
                        .plugin_data
                        .iter()
                        .map(|(key, value)| format!("{key}={value}"))
                        .collect();
                    assert_eq!(written, own_prediction(alt), "table {name}, allele {alt}");
                }
            }
        }
    }

    /// A structural variant enters neither the batch's regions nor the annotate
    /// loop: a `<DEL>` alone, its span covering the record, yields an empty
    /// prefetch and takes no field even on a missense row; beside a missense
    /// SNV inside its span, the deletion's end lying past the index's last
    /// addressable position (2^29 for a `.tbi`), the SNV keeps its fields
    /// and the deletion, its span left out by the flag, takes none, nor does a
    /// copy of the SNV carrying the flag.
    #[test]
    fn structural_variants_are_left_out_of_the_sweep() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &["chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA2B\t0.9123\tlikely_pathogenic"],
        );
        let plugin = loaded(&[format!("file={path}")]);
        // A deletion as the VCF parser builds one, carrying the row the record
        // would annotate on an SNV.
        let deletion = |end: u64| {
            let mut sv =
                InputVariant::new("21".into(), 9068000, end, b"N".to_vec(), b"<DEL>".to_vec());
            sv.is_structural = true;
            sv.transcript_consequences.push(row(
                "ENST00000400001",
                Consequence::MissenseVariant,
                "A/B",
            ));
            sv
        };

        let batch: BatchQueryResult = plugin
            .prefetch(&[deletion(9069000)])
            .expect("prefetch")
            .expect("the store is open")
            .downcast()
            .expect("a batch result");
        assert!(batch.is_empty());
        assert_eq!(fields(&plugin, deletion(9069000)), [Vec::<String>::new()]);

        let mut v = snv(9068417, b"C", b"A");
        v.transcript_consequences
            .push(row("ENST00000400001", Consequence::MissenseVariant, "A/B"));
        // The same SNV carrying the flag pins the loop's filter: its lookup
        // finds the record and its alleles match it, so only the flag leaves it
        // out.
        let mut flagged = v.clone();
        flagged.is_structural = true;
        let mut variants = vec![deletion(600_000_000), v, flagged];
        plugin.run_batch(&mut variants).expect("run_batch");
        let written: Vec<Vec<String>> = variants
            .iter()
            .map(|v| {
                v.transcript_consequences[0]
                    .plugin_data
                    .iter()
                    .map(|(key, value)| format!("{key}={value}"))
                    .collect()
            })
            .collect();
        assert_eq!(
            written,
            [
                vec![],
                vec![
                    "am_pathogenicity=0.9123".to_string(),
                    "am_class=likely_pathogenic".to_string(),
                ],
                vec![],
            ]
        );
    }

    /// `AlphaMissense.pm:290` splits the line htslib hands over, the line end
    /// alone removed, with `split /\t/`, and `:304` copies the cells into the
    /// prediction's columns as they are, so the `am_class` cell keeps the
    /// whitespace before the line end and a whitespace-only cell is the class:
    /// `likely_pathogenic ` prints as `am_class=likely_pathogenic ` and ` ` as
    /// `am_class= `.
    #[test]
    fn a_trailing_space_on_the_last_cell_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &[
                "chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA2B\t0.9123\tlikely_pathogenic ",
                "chr21\t9068418\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA3B\t0.1234\t ",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let mut v = snv(9068417, b"C", b"A");
        v.transcript_consequences
            .push(row("ENST00000400001", Consequence::MissenseVariant, "A/B"));
        assert_eq!(
            fields(&plugin, v),
            [["am_pathogenicity=0.9123", "am_class=likely_pathogenic "]]
        );
        let mut v = snv(9068418, b"C", b"A");
        v.transcript_consequences
            .push(row("ENST00000400001", Consequence::MissenseVariant, "A/B"));
        assert_eq!(
            fields(&plugin, v),
            [["am_pathogenicity=0.1234", "am_class= "]]
        );
    }

    /// `AlphaMissense.pm:175-181` reads `key=value` pairs alone and dies without
    /// `file`, so a file given by position is the file given without its key.
    #[test]
    fn positional_file_is_refused() {
        assert_eq!(
            init_message(Some(Assembly::Grch38), &["/data/AlphaMissense_hg38.tsv.gz"]),
            "AlphaMissense: no file specified; use file=<path>"
        );
        assert_eq!(
            init_message::<&str>(Some(Assembly::Grch38), &[]),
            "AlphaMissense: no file specified; use file=<path>"
        );
        assert_eq!(
            init_message(Some(Assembly::Grch38), &["cols=all"]),
            "AlphaMissense: no file specified; use file=<path>"
        );
        assert_eq!(
            init_message(Some(Assembly::Grch38), &["file=a.tsv.gz", "no_match=1"]),
            "AlphaMissense: unknown parameter 'no_match'; AlphaMissense accepts file, cols, \
             transcript_match"
        );
    }

    /// Both AlphaMissense builds use UCSC contigs and identical column layouts,
    /// so a wrong-build file resolves contigs and then annotates nothing at exit
    /// 0. The `genome` column is the only in-file signal, so the guard must read
    /// it and refuse, and must still accept the matching build.
    #[test]
    fn refuses_the_other_builds_file_and_accepts_its_own() {
        let dir = tempfile::tempdir().unwrap();
        let hg38 = fixture(
            dir.path(),
            "am_hg38.tsv.gz",
            &["chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST1\tA1B\t0.9123\tlikely_pathogenic"],
        );
        let hg19 = fixture(
            dir.path(),
            "am_hg19.tsv.gz",
            &["chr21\t9907250\tC\tA\thg19\tQ9NQ38\tENST1\tA1B\t0.9123\tlikely_pathogenic"],
        );

        let message = init_message(Some(Assembly::Grch37), &[format!("file={hg38}")]);
        assert!(message.contains("hg38"), "got: {message}");
        let message = init_message(Some(Assembly::Grch38), &[format!("file={hg19}")]);
        assert!(message.contains("hg19"), "got: {message}");

        loaded(&[format!("file={hg38}")]);
        let mut plugin = AlphaMissensePlugin::new();
        plugin
            .init(&[format!("file={hg19}")])
            .expect("no --assembly must not block loading");
    }

    /// Every selected column is a key of the matched row's result
    /// (AlphaMissense.pm:279, `%res = map { $_ => $res{$_} } @{ $self->{cols} }`):
    /// a row short of its `am_class` cell writes the key with the undefined
    /// value, the empty value here (the `split /\t/` of `:290` drops the
    /// trailing empty cells), and a row whose `am_pathogenicity` cell is empty
    /// before a filled `am_class` writes the empty string, named in the row's
    /// `plugin_empty_strings`.
    #[test]
    fn a_matched_row_short_of_a_selected_cell_writes_the_key_with_the_empty_value() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &[
                "chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA2B\t0.9123",
                "chr21\t9068420\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA3B\t\tlikely_benign",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let mut short = snv(9068417, b"C", b"A");
        short.transcript_consequences.push(row(
            "ENST00000400001",
            Consequence::MissenseVariant,
            "A/B",
        ));
        let mut empty = snv(9068420, b"C", b"A");
        empty.transcript_consequences.push(row(
            "ENST00000400001",
            Consequence::MissenseVariant,
            "A/B",
        ));
        let mut variants = vec![short, empty];
        plugin.run_batch(&mut variants).unwrap();
        let written = |variant: &InputVariant| -> Vec<String> {
            variant.transcript_consequences[0]
                .plugin_data
                .iter()
                .map(|(key, value)| format!("{key}={value}"))
                .collect()
        };
        assert_eq!(
            written(&variants[0]),
            ["am_pathogenicity=0.9123", "am_class="]
        );
        assert!(variants[0].transcript_consequences[0]
            .plugin_empty_strings
            .is_empty());
        assert_eq!(
            written(&variants[1]),
            ["am_pathogenicity=", "am_class=likely_benign"]
        );
        assert_eq!(
            variants[1].transcript_consequences[0].plugin_empty_strings,
            ["am_pathogenicity"]
        );
    }

    /// A row's POS cell is keyed as text by `get_matched_variant_alleles`
    /// (`AlphaMissense.pm:263-275`; ensembl-variation Utils/Sequence.pm:1232)
    /// when the alleles lose no leading base: `9068417` annotates the variant
    /// at 9068417, `9068417.0` and `9068417.5` do not, and `9068418x` does not
    /// annotate the variant at 9068418. The tabix index places each at its
    /// leading integer, so every row is returned and the text alone decides.
    #[test]
    fn a_position_cell_is_matched_as_the_text_the_key_holds() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "cells.tsv.gz",
            &[
                "chr21\t9068417\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA2B\t0.91\tlikely_pathogenic",
                "chr21\t9068417.0\tC\tG\thg38\tQ9NQ38\tENST00000400001.4\tA2C\t0.92\tlikely_pathogenic",
                "chr21\t9068417.5\tC\tT\thg38\tQ9NQ38\tENST00000400001.4\tA2D\t0.93\tlikely_pathogenic",
                "chr21\t9068418x\tC\tA\thg38\tQ9NQ38\tENST00000400001.4\tA2E\t0.94\tlikely_pathogenic",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let classes = |pos: u64, alt: &[u8], change: &str| {
            let mut variant = snv(pos, b"C", alt);
            variant.transcript_consequences.push(row(
                "ENST00000400001",
                Consequence::MissenseVariant,
                change,
            ));
            fields(&plugin, variant)
        };
        assert_eq!(
            classes(9068417, b"A", "A/B"),
            [["am_pathogenicity=0.91", "am_class=likely_pathogenic"]]
        );
        assert_eq!(classes(9068417, b"G", "A/C"), [Vec::<String>::new()]);
        assert_eq!(classes(9068417, b"T", "A/D"), [Vec::<String>::new()]);
        assert_eq!(classes(9068418, b"A", "A/E"), [Vec::<String>::new()]);
        // An indel-shaped row steps its cell as `parse_data` does (`:297`,
        // `$start++` on a string): `9068430` and `9068430.0` reach the deletion
        // at 9068431, `09068430` steps to the text `09068431` and `abc` to
        // `abd`, which no start equals, and `0` steps to `1`.
        // The row sits where tabix indexes its cell, at the leading integer.
        let indel_row = |cell: &str, indexed: u64| TabixRecord {
            chr: "chr21".to_string(),
            start: indexed,
            end: indexed + 1,
            columns: vec!["chr21".into(), cell.into(), "GA".into(), "G".into()],
            header: None,
        };
        let matched = |cell: &str, indexed: u64, start: u64| {
            let batch = BatchQueryResult::from_records(HashMap::from([(
                "21".to_string(),
                BTreeMap::from([(indexed, vec![indel_row(cell, indexed)])]),
            )]));
            let deletion =
                InputVariant::new("21".into(), start, start, b"A".to_vec(), b"-".to_vec());
            !plugin.matching_records(&batch, &deletion).is_empty()
        };
        assert!(matched("9068430", 9068430, 9068431));
        assert!(matched("9068430.0", 9068430, 9068431));
        assert!(!matched("09068430", 9068430, 9068431));
        assert!(!matched("abc", 0, 1));
        assert!(matched("0", 0, 1));
    }

    /// The rows are read from two bases before the variant's start to its end
    /// (`AlphaMissense.pm:248`), so a multi-base row written one base before
    /// the variant, `GG/GA` at 13610378, is read and trimmed of its shared
    /// anchor base to `G/A` at 13610379 (`get_matched_variant_alleles`,
    /// Utils/Sequence.pm:1199-1240), the key the variant's `G/A` holds: the A
    /// allele takes its prediction, the T allele nothing from that row and its
    /// own single-base row. A single-base row one base before a variant with
    /// the same alleles keys its own position and matches nothing.
    #[test]
    fn a_multi_base_row_one_base_before_the_variant_is_trimmed_to_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &[
                "chr21\t13610378\tGG\tGA\thg38\tQ86YR6\tENST00000299443.5\tD51N\t0.0926\tlikely_benign",
                "chr21\t13610379\tG\tT\thg38\tQ86YR6\tENST00000299443.5\tD51Y\t0.1494\tlikely_benign",
                "chr21\t13610390\tG\tA\thg38\tQ86YR6\tENST00000299443.5\tD55N\t0.5000\tambiguous",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let allele = |alt: &[u8], change: &str| {
            let mut v = snv(13610379, b"G", alt);
            v.transcript_consequences.push(row(
                "ENST00000299443",
                Consequence::MissenseVariant,
                change,
            ));
            v
        };
        assert_eq!(
            fields(&plugin, allele(b"A", "D/N")),
            [["am_pathogenicity=0.0926", "am_class=likely_benign"]]
        );
        assert_eq!(
            fields(&plugin, allele(b"T", "D/Y")),
            [["am_pathogenicity=0.1494", "am_class=likely_benign"]]
        );
        assert_eq!(fields(&plugin, allele(b"T", "D/N")), [Vec::<String>::new()]);
        // The single-base row at 13610390 lies in the window of the variant at
        // 13610391 and spells its alleles, and the variant takes nothing: the
        // row keys 13610390, the variant 13610391.
        let mut v = snv(13610391, b"G", b"A");
        v.transcript_consequences
            .push(row("ENST00000299443", Consequence::MissenseVariant, "D/N"));
        assert_eq!(fields(&plugin, v), [Vec::<String>::new()]);
    }

    /// A variant at the first positions of a chromosome reads the rows from
    /// position 1, the window clamped where `$vf->{start} - 2` reaches 0 or
    /// below (`AlphaMissense.pm:248`; the documented class
    /// `tabix_window_below_position_one`): a row at 1 annotates the variant at
    /// 1, and a two-base row at 1 trimmed to position 2 annotates the variant at
    /// 2.
    #[test]
    fn the_window_is_clamped_at_position_one() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "am.tsv.gz",
            &[
                "chr21\t1\tC\tG\thg38\tQ9NQ38\tENST00000400001.4\tA1B\t0.5000\tambiguous",
                "chr21\t1\tCC\tCG\thg38\tQ9NQ38\tENST00000400001.4\tA2B\t0.6000\tlikely_pathogenic",
            ],
        );
        let plugin = loaded(&[format!("file={path}")]);
        let at = |pos: u64| {
            let mut v = snv(pos, b"C", b"G");
            v.transcript_consequences.push(row(
                "ENST00000400001",
                Consequence::MissenseVariant,
                "A/B",
            ));
            v
        };
        assert_eq!(
            fields(&plugin, at(1)),
            [["am_pathogenicity=0.5000", "am_class=ambiguous"]]
        );
        assert_eq!(
            fields(&plugin, at(2)),
            [["am_pathogenicity=0.6000", "am_class=likely_pathogenic"]]
        );
        let batch: BatchQueryResult = plugin
            .prefetch(&[at(1), at(2)])
            .expect("prefetch")
            .expect("the store is open")
            .downcast()
            .expect("a batch result");
        assert_eq!(batch.lookup("21", 1, 2).len(), 2);
    }
}
