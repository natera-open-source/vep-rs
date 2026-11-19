// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! REVEL plugin: the Rare Exome Variant Ensemble Learner missense score, read
//! from the tabix-indexed REVEL table and written as Ensembl VEP's `REVEL.pm`
//! writes it.
//!
//! # Parameters
//!
//! `--plugin REVEL,<file>[,<no_match>]` by position or
//! `--plugin REVEL,file=<path>[,no_match=<value>]` by key (`REVEL.pm:105-110`).
//! `no_match` is read as Ensembl VEP reads a Perl value: every value but the
//! empty string and `0` sets it. The run's assembly (`--assembly`, else the
//! cache's) is required (`:132-133`).
//!
//! # Data file
//!
//! One download serves both assemblies. The header is the last of the file's
//! leading `#` lines, its first `#` removed and its names split on whitespace
//! (`:116-125`); a table whose header line lacks its `#`, or whose `#` line
//! names nothing, is refused, as a table with no header line is. The header's
//! column count fixes the layout (`:127-130`, `:209-241`): 7 columns are `chr,
//! hg19_pos, ref, alt, aaref, aaalt, REVEL`; 8 and 9 columns are `chr,
//! hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL, Ensembl_transcriptid`,
//! the ninth read from each row that has it (`:210`), whether or not the header
//! names it. Any other count is refused. The header must name the run's
//! assembly's position column (`hg19_pos` under GRCh37, `grch38_pos` under
//! GRCh38; `:135-142`), and that column must be the one the tabix index was
//! built on: the two prepared copies of the table carry the same header and
//! pass the same checks, and reading the other column would discard every
//! record at exit 0.
//!
//! # Matching
//!
//! A row is annotated when it is a missense row (`:162`) and a record at the
//! variant's position carries its alleles (`:176-188`), an `aaref` equal to
//! the row's reference peptide and an `aaalt` equal to the row's alternate
//! peptide (`:193`, `:197`). A record scores one substitution, from its
//! `aaref` to its `aaalt`: where the transcript model reads another residue at
//! the codon than the record's `aaref`, the record scores a substitution the
//! row does not report and is passed over. `REVEL.pm` reads the `aaref` cell
//! (`:210`, `:229`) and tests nothing with it, so it writes such a record's
//! score. Both cells are compared as written: the published table spells its
//! residues with the twenty one-letter codes alone, as the row does. When the
//! record names transcripts and `no_match` is unset, one of its `;`-separated
//! ids must be the row's transcript as well (`:191-195`); a record without
//! transcript ids, and every record under `no_match`, matches on the two
//! peptides alone. A file id is
//! compared as written with the transcript's stable id, which carries no
//! version, as the published table's ids do not. The first matching record's
//! score is written as the file spells it; a matching record whose score cell
//! lies at or past its last filled cell, the undefined value `split /\t/`
//! leaves (`:207`, `:223`), writes the empty value.
//!
//! The records of the variant's span are put to the match in file order, each
//! before the next is read (`:170-188`), and `get_matched_variant_alleles`
//! dies on a record whose reference or alternate cell is undefined
//! (`trim_sequences`, ensembl-variation Utils/Sequence.pm:968, from `:1230`):
//! a missense row whose walk meets such a record before its first full match
//! takes no field and the plugin reports the failure, while the variant's
//! other rows keep theirs, as under Ensembl VEP, where each output line's
//! `run` dies on its own. A record after a row's first full match is never
//! read for that row.
//!
//! # Output fields
//!
//! - `REVEL`: Rare Exome Variant Ensemble Learner (`:156`).
//!
//! `REVEL.pm` is the plugin's Ensembl VEP counterpart in Ensembl/VEP_plugins
//! release/116.

use std::path::Path;

use vep_core::assembly::Assembly;
use vep_core::consequence::{Consequence, TranscriptConsequence};
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

use crate::annotation_store::{check_table_or_store, open_best_store_with, AnnotationStore};
use crate::params::PluginParams;
use crate::tabix::{
    allele_match, last_preamble_line, resolve_indexed_position_column, split_field_count,
    table_contig, BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord,
};
use crate::traits::{warn_on_stderr, BuiltinPlugin, PluginError, PrefetchData};

const NAME: &str = "REVEL";

/// The output field (`REVEL.pm:156`).
const FIELD: &str = "REVEL";

/// Header name of the GRCh37 position column (`REVEL.pm:135`).
const HG19_POS_COLUMN: &str = "hg19_pos";

/// Header name of the GRCh38 position column (`REVEL.pm:136`).
const GRCH38_POS_COLUMN: &str = "grch38_pos";

/// The 0-based columns of one of the table's layouts, fixed by the header's
/// column count as `REVEL.pm:209-241` fixes them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Layout {
    reference: usize,
    alternate: usize,
    /// The `aaref` column, which `REVEL.pm:210` and `:229` read into `$refaa`
    /// and no line tests.
    reference_peptide: usize,
    alternate_peptide: usize,
    score: usize,
    /// The column of a row's transcript ids, read from a row that has it
    /// (`:210` lists nine variables whatever the header's count).
    transcript_ids: Option<usize>,
}

impl Layout {
    /// The layout of a header with `columns` columns; `None` for a count the
    /// plugin refuses (`REVEL.pm:127-129`).
    fn for_column_count(columns: usize) -> Option<Self> {
        match columns {
            7 => Some(Self {
                reference: 2,
                alternate: 3,
                reference_peptide: 4,
                alternate_peptide: 5,
                score: 6,
                transcript_ids: None,
            }),
            8 | 9 => Some(Self {
                reference: 3,
                alternate: 4,
                reference_peptide: 5,
                alternate_peptide: 6,
                score: 7,
                transcript_ids: Some(8),
            }),
            _ => None,
        }
    }
}

/// A record of the variant's span as the loop of `REVEL.pm:170-199` meets it,
/// in file order: one carrying the variant's alleles, whose transcript and
/// peptide cells the loop then tests, or one with no reference or alternate
/// allele cell, which `get_matched_variant_alleles` dies on before any test
/// (`trim_sequences`, ensembl-variation Utils/Sequence.pm:968, from `:1230`).
/// A record with other alleles is passed over and left out.
enum Candidate<'a> {
    Alleles(&'a TabixRecord),
    Undefined(&'a TabixRecord),
}

/// A matched record's score cell as `REVEL.pm:223` copies it into the result:
/// the cell's text, and whether an empty text is the defined empty string
/// `split /\t/` leaves for an empty cell before a filled one, rather than the
/// undefined value of a cell at or past the row's last filled one; the tab and
/// JSON writers print the two apart (`Tab.pm:128`, `OutputFactory/JSON.pm:327`).
struct Score<'a> {
    text: &'a str,
    empty_string: bool,
}

/// REVEL score plugin for missense rows.
pub struct RevelPlugin {
    /// Annotation store (tabix or binary), opened by `init`.
    store: Option<Box<dyn AnnotationStore>>,
    /// The table's path as given, named by a failure on one of its rows.
    file: String,
    /// Write the score without matching the record's transcripts to the row's.
    no_match: bool,
    /// The table's columns, from the header's column count.
    layout: Layout,
    /// Column of the assembly's position, the one the index is built on and the
    /// one the match reads (`$revel_start`, `REVEL.pm:172`).
    pos_col: usize,
    /// Target assembly from `--assembly` or the cache, injected before `init()`;
    /// selects which of the table's position columns the index was built on.
    assembly: Option<Assembly>,
    /// Whether the run reads the RefSeq or merged transcript set.
    other_transcript_set: bool,
}

impl Default for RevelPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl RevelPlugin {
    /// Create an uninitialized REVEL plugin.
    ///
    /// Call [`init`](BuiltinPlugin::init) with the file path before use.
    pub fn new() -> Self {
        Self {
            store: None,
            file: String::new(),
            no_match: false,
            layout: Layout::for_column_count(9).expect("nine columns is a layout"),
            pos_col: 0,
            assembly: None,
            other_transcript_set: false,
        }
    }

    /// Header name of the position column for the requested assembly.
    fn position_column_for(assembly: Assembly) -> &'static str {
        match assembly {
            Assembly::Grch37 => HG19_POS_COLUMN,
            Assembly::Grch38 => GRCH38_POS_COLUMN,
        }
    }

    /// Whether a transcript consequence includes missense_variant (`REVEL.pm:162`).
    fn is_missense(consequence: &TranscriptConsequence) -> bool {
        consequence
            .consequences
            .contains(&Consequence::MissenseVariant)
    }

    /// Whether `record` names the row's transcript, when it has to
    /// (`REVEL.pm:191-197`): a record whose transcript column is empty or
    /// absent, and every record under `no_match`, needs no transcript match.
    fn transcript_matches(&self, record: &TabixRecord, transcript_id: &str) -> bool {
        if self.no_match {
            return true;
        }
        let Some(file_transcripts) = self
            .layout
            .transcript_ids
            .and_then(|column| record.columns.get(column))
            .filter(|ids| perl_true(ids))
        else {
            return true;
        };
        // Each id is compared as `split /;/` leaves it with the transcript's
        // stable id (`:166`, `:192-193`), whitespace and version included.
        file_transcripts.split(';').any(|id| id == transcript_id)
    }

    /// The score of the first of `records` whose reference and alternate
    /// peptides are `reference_peptide` and `peptide` and whose transcripts,
    /// when matched, name `transcript_id`, the records walked in file order as
    /// `REVEL.pm:170-201` walks the span's records for one row: `Err` when the
    /// walk meets a record with no allele cell before that match, the record
    /// `get_matched_variant_alleles` dies on (Utils/Sequence.pm:968), and
    /// `Ok(None)` when it ends without a match. The reference peptide is the
    /// one test `REVEL.pm:193` and `:197` leave out: a record whose `aaref`
    /// cell is another residue than the row's, or is empty, scores another
    /// substitution than the row's and is passed over like a record of another
    /// alternate peptide.
    /// The matching record's result is returned whatever its score cell holds
    /// (`:193`, `:197`, `:223`): a cell at or past the row's last filled one is
    /// the undefined value, written as the unmarked empty value.
    fn score_for<'a>(
        &self,
        records: &[Candidate<'a>],
        transcript_id: &str,
        reference_peptide: &str,
        peptide: &str,
    ) -> Result<Option<Score<'a>>, PluginError> {
        for candidate in records {
            let record = match candidate {
                Candidate::Undefined(record) => {
                    return Err(allele_match::undefined_allele_row(NAME, &self.file, record));
                }
                Candidate::Alleles(record) => record,
            };
            let cell_is = |column: usize, residue: &str| {
                record
                    .columns
                    .get(column)
                    .is_some_and(|cell| cell == residue)
            };
            if !cell_is(self.layout.reference_peptide, reference_peptide)
                || !cell_is(self.layout.alternate_peptide, peptide)
                || !self.transcript_matches(record, transcript_id)
            {
                continue;
            }
            let present = split_field_count(&record.columns);
            let text = record
                .columns
                .get(self.layout.score)
                .filter(|_| self.layout.score < present)
                .map_or("", String::as_str);
            return Ok(Some(Score {
                text,
                empty_string: text.is_empty() && self.layout.score < present,
            }));
        }
        Ok(None)
    }

    /// The records of `batch` at `variant` that the loop acts on, in file order
    /// (`REVEL.pm:168-188`): each record of the span is put to
    /// `get_matched_variant_alleles`, which dies on one with no reference or
    /// alternate allele cell and passes over one with other alleles
    /// ([`allele_match::matched_or_undefined`]).
    fn candidate_records<'a>(
        &self,
        batch: &'a BatchQueryResult,
        variant: &InputVariant,
    ) -> Vec<Candidate<'a>> {
        batch
            .lookup(
                table_contig(variant, self.assembly),
                variant.start,
                variant.end,
            )
            .into_iter()
            .filter_map(|record| {
                match allele_match::matched_or_undefined(
                    record,
                    variant,
                    self.pos_col,
                    self.layout.reference,
                    self.layout.alternate,
                ) {
                    None => Some(Candidate::Undefined(record)),
                    Some(true) => Some(Candidate::Alleles(record)),
                    Some(false) => None,
                }
            })
            .collect()
    }
}

impl BuiltinPlugin for RevelPlugin {
    fn name(&self) -> &str {
        NAME
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            FIELD.to_string(),
            "Rare Exome Variant Ensemble Learner".to_string(),
        )]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(NAME, &["file", "no_match"], params)?;
        let file = params
            .get("file")
            .or_else(|| params.positional().first().map(String::as_str))
            .filter(|file| !file.is_empty())
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: no file specified; use file=<path> or give the file first"
                ))
            })?;
        self.no_match = params
            .get("no_match")
            .or_else(|| params.positional().get(1).map(String::as_str))
            .is_some_and(perl_true);
        self.file = file.to_string();
        let path = Path::new(file);
        // `add_file` runs before the header is read (REVEL.pm:115-116).
        check_table_or_store(path)?;

        // The position column is fixed in `TabixConfig` at construction and the
        // engine re-parses each record's position from it, so it is resolved
        // before the store opens. The header is the last of the file's leading
        // `#` lines, the lines `tabix -h` prints, with its first `#` removed and
        // its names split as Perl's bare `split` splits a byte string, on ASCII
        // whitespace (`REVEL.pm:116-125`), so a header line without its `#` is
        // no header, nor is a `#` line naming nothing, and a second `#` is a
        // name.
        let header: Vec<String> = last_preamble_line(path)?
            .as_deref()
            .and_then(|line| line.strip_prefix('#'))
            .unwrap_or_default()
            .split_ascii_whitespace()
            .map(str::to_string)
            .collect();
        if header.is_empty() {
            return Err(PluginError::Init(format!(
                "{NAME}: could not read headers from {file}"
            )));
        }
        self.layout = Layout::for_column_count(header.len()).ok_or_else(|| {
            PluginError::Init(format!(
                "{NAME}: column count must be 8 or 9 for REVEL files with GRCh38 positions \
                 or 7 for REVEL files with GRCh37 positions only; {file} has {} columns",
                header.len()
            ))
        })?;

        let Some(assembly) = self.assembly else {
            return Err(PluginError::Init(format!(
                "{NAME}: specify assembly using --assembly GRCh37 or --assembly GRCh38"
            )));
        };
        let wanted = Self::position_column_for(assembly);
        let start_col = resolve_indexed_position_column(path, &header, wanted).map_err(
            |error| match error {
                PluginError::Init(message) => PluginError::Init(format!(
                    "{NAME}: assembly is {assembly}, which reads the '{wanted}' column, but {message}"
                )),
                other => other,
            },
        )?;
        self.pos_col = start_col;

        let config = TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: 0,
            start_col,
            end_col: None,
            ref_col: None,
            alt_col: None,
            zero_based: false,
        };
        // A row's last cell is read as `split /\t/` reads it off the line
        // (`REVEL.pm:207`), its trailing whitespace included: the score of an
        // eight-column file is written with it, and the transcript cell of a
        // nine-column file holding whitespace alone is a cell that matches no id.
        self.store = Some(open_best_store_with(
            path,
            config,
            self.assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )?);
        // After the index check, as `REVEL.pm:115-146` orders them.
        if self.other_transcript_set && !self.no_match {
            warn_on_stderr(
                "REVEL transcript matching may fail with --refseq/--merged; consider setting \
                 no_match=1",
            );
        }
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
        // Every sequence variant's span is queried, missense row or not: the
        // spans of a batch merge into a few contiguous regions the file is read
        // through once, where the missense spans alone form many small regions,
        // each a seek and a block inflated on its own. A structural variant's
        // alleles are class words no single-base row matches, so its span,
        // which can cover megabases of rows, is left out.
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|v| !v.is_structural)
            .map(|v| (table_contig(v, self.assembly), v.start, v.end))
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
        let mut failure = None;
        for variant in variants.iter_mut().filter(|v| !v.is_structural) {
            // The records come first: most variants of a genome-wide input lie
            // outside the file's coding positions and are done at the lookup,
            // before any of their rows is read.
            let records = self.candidate_records(&batch, variant);
            if records.is_empty() {
                continue;
            }
            for tc in &mut variant.transcript_consequences {
                if !Self::is_missense(tc) {
                    continue;
                }
                let Some((reference_peptide, peptide)) = tc.peptides() else {
                    continue;
                };
                match self.score_for(&records, tc.transcript_id(), reference_peptide, peptide) {
                    Ok(Some(score)) => {
                        tc.plugin_data
                            .insert(FIELD.to_string(), score.text.to_string());
                        // A later plugin's value replaces an earlier one under
                        // the same key (OutputFactory.pm:2230-2231), its mark
                        // too.
                        tc.plugin_empty_strings.retain(|key| key != FIELD);
                        if score.empty_string {
                            tc.plugin_empty_strings.push(FIELD.to_string());
                        }
                    }
                    Ok(None) => {}
                    Err(error) => {
                        failure.get_or_insert(error);
                    }
                }
            }
        }
        match failure {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }

    /// Each output line's `run` walks the span's records and dies on its own
    /// (`REVEL.pm:170-201`, each line's call evaluated apart by Ensembl VEP's
    /// plugin runner), so a line matched before the undefined row keeps
    /// `REVEL`: a failure costs the rows it failed on alone.
    fn annotate_fails_per_row(&self) -> bool {
        true
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
    use std::sync::Arc;
    use vep_core::transcript::Transcript;

    /// A REVEL row as published: one file, two coordinate columns.
    /// `chr, hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL, Ensembl_transcriptid`
    const REVEL_HEADER: &str =
        "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid";

    /// A GRCh38-indexed fixture (`tabix -s 1 -b 3 -e 3`) of `lines`.
    fn grch38_fixture(dir: &Path, name: &str, lines: &[&str]) -> String {
        write_tabix_fixture(dir, name, lines, IndexSpec::standard().with_start_col(3))
            .display()
            .to_string()
    }

    fn loaded<S: AsRef<str>>(assembly: Assembly, params: &[S]) -> RevelPlugin {
        let mut plugin = RevelPlugin::new();
        plugin.set_assembly(Some(assembly));
        let params: Vec<String> = params.iter().map(|p| p.as_ref().to_string()).collect();
        plugin.init(&params).expect("fixture loads");
        plugin
    }

    fn init_message<S: AsRef<str>>(assembly: Option<Assembly>, params: &[S]) -> String {
        let mut plugin = RevelPlugin::new();
        plugin.set_assembly(assembly);
        let params: Vec<String> = params.iter().map(|p| p.as_ref().to_string()).collect();
        match plugin.init(&params) {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected PluginError::Init, got {other:?}"),
        }
    }

    /// A missense row of `transcript_id` whose amino-acid change is `amino_acids`.
    fn missense_row(transcript_id: &str, amino_acids: &str) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: transcript_id.into(),
                ..Transcript::default()
            })),
            consequences: smallvec![Consequence::MissenseVariant],
            amino_acids: Some(amino_acids.to_string()),
            ..Default::default()
        }
    }

    fn variant(chr: &str, pos: u64, r: &[u8], a: &[u8]) -> InputVariant {
        InputVariant::new(chr.into(), pos, pos, r.to_vec(), a.to_vec())
    }

    /// The `REVEL` field of every row of `variant` after a batch run.
    fn scores(plugin: &RevelPlugin, variant: InputVariant) -> Vec<Option<String>> {
        let mut variants = vec![variant];
        plugin.run_batch(&mut variants).expect("run_batch");
        variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.get(FIELD).cloned())
            .collect()
    }

    /// `add_file` checks for `<file>.tbi` before the header is read
    /// (REVEL.pm:115-116): a table lacking both is refused naming the index.
    #[test]
    fn table_without_its_tbi_is_refused_before_its_header_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "revel.tsv.gz",
            &["#chr\thg19_pos", "21\t100"],
            IndexSpec::standard(),
        );
        std::fs::remove_file(dir.path().join("revel.tsv.gz.tbi")).unwrap();
        let message = init_message(
            Some(Assembly::Grch38),
            &[format!("file={}", path.display())],
        );
        assert!(message.contains("revel.tsv.gz.tbi not found"), "{message}");
    }

    #[test]
    fn new_plugin_is_uninitialized() {
        let plugin = RevelPlugin::new();
        assert_eq!(plugin.name(), "REVEL");
        assert!(plugin.store.is_none());
        assert!(!plugin.no_match);
    }

    /// `REVEL.pm:156` names the one field `REVEL`.
    #[test]
    fn header_info_field_is_revel() {
        let headers = RevelPlugin::new().header_info();
        assert_eq!(
            headers,
            vec![(
                "REVEL".to_string(),
                "Rare Exome Variant Ensemble Learner".to_string()
            )]
        );
    }

    /// `REVEL.pm:193` and `:197` write a record's score only when its `aaalt`
    /// equals the row's alternate peptide: two records of one allele with
    /// different `aaalt` score the rows whose peptides they name, and a
    /// synonymous row is never scored (`:162`).
    #[test]
    fn run_requires_alternate_amino_acid_match() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001",
                "21\t100200\t500600\tA\tG\tQ\tL\t0.311\tENST00000400002",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&format!("file={path}")]);

        let mut v = variant("21", 500600, b"A", b"G");
        v.transcript_consequences
            .push(missense_row("ENST00000400001", "Q/H"));
        v.transcript_consequences
            .push(missense_row("ENST00000400002", "Q/L"));
        v.transcript_consequences
            .push(missense_row("ENST00000400001", "Q/L"));
        v.transcript_consequences.push(TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: "ENST00000400001".into(),
                ..Transcript::default()
            })),
            consequences: smallvec![Consequence::SynonymousVariant],
            amino_acids: Some("Q".to_string()),
            ..Default::default()
        });
        assert_eq!(
            scores(&plugin, v),
            [
                Some("0.742".to_string()),
                Some("0.311".to_string()),
                None,
                None
            ]
        );
    }

    /// A record scores the substitution from its `aaref` to its `aaalt`
    /// (`REVEL.pm:210` and `:229` read both cells; `:193` and `:197` test the
    /// second alone): a record naming the row's alternate peptide beside
    /// another reference peptide than the row's scores nothing, in the
    /// seven-column layout and in the eight- and nine-column one, on either
    /// assembly's position column, with the transcript match and under
    /// `no_match`, while a record of the row's own pair scores.
    #[test]
    fn run_requires_reference_amino_acid_match() {
        let dir = tempfile::tempdir().unwrap();
        let nine = [
            REVEL_HEADER,
            "21\t100200\t500600\tA\tG\tR\tH\t0.742\tENST00000400001",
            "21\t100300\t500700\tA\tG\tQ\tH\t0.311\tENST00000400001",
        ];
        let eight = [
            "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL",
            "21\t100200\t500600\tA\tG\tR\tH\t0.742",
            "21\t100300\t500700\tA\tG\tQ\tH\t0.311",
        ];
        let seven = [
            "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL",
            "21\t100200\tA\tG\tR\tH\t0.742",
            "21\t100300\tA\tG\tQ\tH\t0.311",
        ];
        let grch37 = |name: &str, lines: &[&str]| {
            write_tabix_fixture(dir.path(), name, lines, IndexSpec::standard())
                .display()
                .to_string()
        };
        let tables = [
            (
                Assembly::Grch38,
                grch38_fixture(dir.path(), "nine38.tsv.gz", &nine),
                500600,
            ),
            (
                Assembly::Grch38,
                grch38_fixture(dir.path(), "eight38.tsv.gz", &eight),
                500600,
            ),
            (Assembly::Grch37, grch37("nine37.tsv.gz", &nine), 100200),
            (Assembly::Grch37, grch37("eight37.tsv.gz", &eight), 100200),
            (Assembly::Grch37, grch37("seven37.tsv.gz", &seven), 100200),
        ];
        let row = |pos: u64| {
            let mut v = variant("21", pos, b"A", b"G");
            v.transcript_consequences
                .push(missense_row("ENST00000400001", "Q/H"));
            v
        };
        for (assembly, file, other_reference) in &tables {
            for params in [vec![file.as_str()], vec![file.as_str(), "1"]] {
                let plugin = loaded(*assembly, &params);
                assert_eq!(
                    scores(&plugin, row(*other_reference)),
                    [None],
                    "{file}, no_match {}: the record names another reference peptide",
                    plugin.no_match
                );
                assert_eq!(
                    scores(&plugin, row(other_reference + 100)),
                    [Some("0.311".to_string())],
                    "{file}, no_match {}: the record names the row's own pair",
                    plugin.no_match
                );
            }
        }
    }

    /// The walk goes on in file order past a record of another reference
    /// peptide, as it does past one of another alternate peptide
    /// (`REVEL.pm:170-199`): the row takes the score of the first record of
    /// its own pair, under `no_match` whatever transcripts that record names
    /// and with the transcript match the first one naming its transcript.
    #[test]
    fn a_later_record_of_the_rows_own_pair_scores_past_one_of_another_reference_peptide() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\tR\tH\t0.742\tENST00000400009",
                "21\t100200\t500600\tA\tG\tQ\tH\t0.311\tENST00000400008",
                "21\t100200\t500600\tA\tG\tQ\tH\t0.512\tENST00000400007",
                "21\t100300\t500700\tA\tG\tR\tH\t0.742\tENST00000400001",
                "21\t100300\t500700\tA\tG\tQ\tH\t0.311\tENST00000400001",
            ],
        );
        let row = |pos: u64, transcript: &str| {
            let mut v = variant("21", pos, b"A", b"G");
            v.transcript_consequences
                .push(missense_row(transcript, "Q/H"));
            v
        };
        let unmatched = loaded(Assembly::Grch38, &[path.as_str(), "1"]);
        assert_eq!(
            scores(&unmatched, row(500600, "ENST00000400001")),
            [Some("0.311".to_string())]
        );
        assert_eq!(
            scores(&unmatched, row(500700, "ENST00000400002")),
            [Some("0.311".to_string())]
        );
        let matched = loaded(Assembly::Grch38, &[&path]);
        assert_eq!(scores(&matched, row(500600, "ENST00000400001")), [None]);
        assert_eq!(scores(&matched, row(500600, "ENST00000400009")), [None]);
        assert_eq!(
            scores(&matched, row(500600, "ENST00000400007")),
            [Some("0.512".to_string())]
        );
        assert_eq!(
            scores(&matched, row(500700, "ENST00000400001")),
            [Some("0.311".to_string())]
        );
    }

    /// Among records of the row's own pair the first in file order scores
    /// (`REVEL.pm:193`, `:197`), with the transcript match and under
    /// `no_match`.
    #[test]
    fn the_first_record_of_the_rows_own_pair_scores() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001;ENST00000400002",
                "21\t100200\t500600\tA\tG\tQ\tH\t0.311\tENST00000400001;ENST00000400003",
            ],
        );
        let row = |transcript: &str| {
            let mut v = variant("21", 500600, b"A", b"G");
            v.transcript_consequences
                .push(missense_row(transcript, "Q/H"));
            v
        };
        let matched = loaded(Assembly::Grch38, &[&path]);
        assert_eq!(
            scores(&matched, row("ENST00000400001")),
            [Some("0.742".to_string())]
        );
        assert_eq!(
            scores(&matched, row("ENST00000400003")),
            [Some("0.311".to_string())]
        );
        let unmatched = loaded(Assembly::Grch38, &[path.as_str(), "1"]);
        assert_eq!(
            scores(&unmatched, row("ENST00000400003")),
            [Some("0.742".to_string())]
        );
    }

    /// Each alternate of a record is put to its own records: the second
    /// alternate's row takes no score from its record of another reference
    /// peptide, beside the first alternate's row, whose record names its pair.
    #[test]
    fn a_second_alternates_record_of_another_reference_peptide_scores_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001",
                "21\t100200\t500600\tA\tT\tR\tL\t0.311\tENST00000400001",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&path]);
        let mut first = variant("21", 500600, b"A", b"G");
        first
            .transcript_consequences
            .push(missense_row("ENST00000400001", "Q/H"));
        let mut second = variant("21", 500600, b"A", b"T");
        second
            .transcript_consequences
            .push(missense_row("ENST00000400001", "Q/L"));
        let mut variants = vec![first, second];
        plugin.run_batch(&mut variants).expect("run_batch");
        assert_eq!(fields(&variants[0]), [Some("0.742")]);
        assert_eq!(fields(&variants[1]), [None]);
    }

    /// The `aaref` cell is compared with the row's reference peptide as the
    /// two are written, the test `REVEL.pm:193` makes of the `aaalt` cell: the
    /// published table spells its residues with the twenty one-letter codes
    /// alone, so a selenocysteine letter, a stop's sign, an unknown residue's
    /// letter and a lower-case letter each name another cell than the row's,
    /// and an empty `aaref` cell names no residue.
    #[test]
    fn the_reference_amino_acid_cell_is_compared_as_written() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\t*\tH\t0.742\tENST00000400001",
                "21\t100300\t500700\tA\tG\tX\tH\t0.742\tENST00000400001",
                "21\t100400\t500800\tA\tG\tC\tH\t0.742\tENST00000400001",
                "21\t100500\t500900\tA\tG\tq\tH\t0.742\tENST00000400001",
                "21\t100600\t501000\tA\tG\t\tH\t0.742\tENST00000400001",
                "21\t100700\t501100\tA\tG\tU\tH\t0.311\tENST00000400001",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&path]);
        let score = |pos: u64, change: &str| {
            let mut v = variant("21", pos, b"A", b"G");
            v.transcript_consequences
                .push(missense_row("ENST00000400001", change));
            scores(&plugin, v)
        };
        assert_eq!(score(500600, "U/H"), [None]);
        assert_eq!(score(500700, "*/H"), [None]);
        assert_eq!(score(500800, "X/H"), [None]);
        assert_eq!(score(500900, "Q/H"), [None]);
        assert_eq!(score(501000, "Q/H"), [None]);
        assert_eq!(score(501100, "U/H"), [Some("0.311".to_string())]);
    }

    /// `REVEL.pm:209-210` read a row's ninth cell as its transcript ids under an
    /// eight-name header as under a nine-name one, so a row that carries one is
    /// matched on it (`:191-195`) and a row without one on the peptides alone;
    /// a seven-name header reads no ninth cell (`:227-228`).
    #[test]
    fn a_rows_ninth_cell_is_its_transcript_ids_whatever_the_headers_count() {
        let dir = tempfile::tempdir().unwrap();
        let eight_names = "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL";
        let path = grch38_fixture(
            dir.path(),
            "revel8h9.tsv.gz",
            &[
                eight_names,
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001",
                "21\t100200\t500601\tA\tG\tQ\tH\t0.311",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&format!("file={path}")]);
        assert_eq!(plugin.layout.transcript_ids, Some(8));

        let mut named = variant("21", 500600, b"A", b"G");
        named
            .transcript_consequences
            .push(missense_row("ENST00000400001", "Q/H"));
        named
            .transcript_consequences
            .push(missense_row("ENST00000400002", "Q/H"));
        assert_eq!(
            scores(&plugin, named),
            [Some("0.742".to_string()), None],
            "the ninth cell names the first row's transcript alone"
        );

        let mut unnamed = variant("21", 500601, b"A", b"G");
        unnamed
            .transcript_consequences
            .push(missense_row("ENST00000400002", "Q/H"));
        assert_eq!(
            scores(&plugin, unnamed),
            [Some("0.311".to_string())],
            "a row with no ninth cell matches on the peptides alone"
        );

        let seven = write_tabix_fixture(
            dir.path(),
            "revel7.tsv.gz",
            &[
                "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL",
                "21\t100200\tA\tG\tQ\tH\t0.742\tENST00000400001",
            ],
            IndexSpec::standard(),
        )
        .display()
        .to_string();
        let plugin = loaded(Assembly::Grch37, &[&format!("file={seven}")]);
        assert_eq!(plugin.layout.transcript_ids, None);
        let mut other = variant("21", 100200, b"A", b"G");
        other
            .transcript_consequences
            .push(missense_row("ENST00000400002", "Q/H"));
        assert_eq!(scores(&plugin, other), [Some("0.742".to_string())]);
    }

    /// Under `no_match` (`REVEL.pm:197`) the transcript ids are not read, so a
    /// record of another transcript scores the row when its peptides match;
    /// the second positional parameter sets it (`:107`), as `no_match=1` does
    /// (`:110`), and `0` or an empty value leaves it unset.
    #[test]
    fn positional_one_sets_no_match() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001;ENST00000400009",
            ],
        );
        let row = || {
            let mut v = variant("21", 500600, b"A", b"G");
            v.transcript_consequences
                .push(missense_row("ENST00000400002", "Q/H"));
            v
        };

        assert_eq!(scores(&loaded(Assembly::Grch38, &[&path]), row()), [None]);
        assert_eq!(
            scores(&loaded(Assembly::Grch38, &[path.as_str(), "1"]), row()),
            [Some("0.742".to_string())]
        );
        assert_eq!(
            scores(
                &loaded(
                    Assembly::Grch38,
                    &[format!("file={path}").as_str(), "no_match=1"]
                ),
                row()
            ),
            [Some("0.742".to_string())]
        );
        assert_eq!(
            scores(
                &loaded(
                    Assembly::Grch38,
                    &[format!("file={path}").as_str(), "no_match=0"]
                ),
                row()
            ),
            [None]
        );
        assert!(!loaded(Assembly::Grch38, &[path.as_str(), ""]).no_match);

        // A `;`-separated id list matches on any member (`:192`), each compared
        // as written: a versioned id is not the unversioned stable id (`:193`).
        let row_of = |id: &str| {
            let mut v = variant("21", 500600, b"A", b"G");
            v.transcript_consequences.push(missense_row(id, "Q/H"));
            v
        };
        let plugin = loaded(Assembly::Grch38, &[&path]);
        assert_eq!(
            scores(&plugin, row_of("ENST00000400009")),
            [Some("0.742".to_string())]
        );
        assert_eq!(scores(&plugin, row_of("ENST00000400009.3")), [None]);
    }

    /// A structural variant enters neither the batch's regions nor the annotate
    /// loop: a `<DEL>` alone, its span covering the record, yields an empty
    /// prefetch and takes no field even on a missense row; beside a missense
    /// SNV inside its span, the deletion's end lying past the index's last
    /// addressable position (2^29 for a `.tbi`), the SNV keeps its score
    /// and the deletion, its span left out by the flag, takes none, nor does a
    /// copy of the SNV carrying the flag.
    #[test]
    fn structural_variants_are_left_out_of_the_sweep() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&format!("file={path}")]);
        // A deletion as the VCF parser builds one, carrying the row the record
        // would score on an SNV.
        let deletion = |end: u64| {
            let mut sv =
                InputVariant::new("21".into(), 500000, end, b"N".to_vec(), b"<DEL>".to_vec());
            sv.is_structural = true;
            sv.transcript_consequences
                .push(missense_row("ENST00000400001", "Q/H"));
            sv
        };

        let batch: BatchQueryResult = plugin
            .prefetch(&[deletion(501000)])
            .expect("prefetch")
            .expect("the store is open")
            .downcast()
            .expect("a batch result");
        assert!(batch.is_empty());
        assert_eq!(scores(&plugin, deletion(501000)), [None]);

        let mut snv = variant("21", 500600, b"A", b"G");
        snv.transcript_consequences
            .push(missense_row("ENST00000400001", "Q/H"));
        // The same SNV carrying the flag pins the loop's filter: its lookup
        // finds the record and its alleles match it, so only the flag leaves it
        // out.
        let mut flagged = snv.clone();
        flagged.is_structural = true;
        let mut variants = vec![deletion(600_000_000), snv, flagged];
        plugin.run_batch(&mut variants).expect("run_batch");
        let fields: Vec<Option<&str>> = variants
            .iter()
            .map(|v| {
                v.transcript_consequences[0]
                    .plugin_data
                    .get(FIELD)
                    .map(String::as_str)
            })
            .collect();
        assert_eq!(fields, [None, Some("0.742"), None]);
    }

    /// `REVEL.pm:127-129` dies on a header of any count but 7, 8 or 9; a
    /// seven-column file is the GRCh37-only layout and an eight-column file the
    /// two-assembly layout whose rows carry no ninth cell, both matched on the
    /// peptides alone.
    #[test]
    fn column_count_outside_7_8_9_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let six = grch38_fixture(
            dir.path(),
            "six.tsv.gz",
            &[
                "#chr\thg19_pos\tgrch38_pos\tref\talt\tREVEL",
                "21\t100200\t500600\tA\tG\t0.742",
            ],
        );
        let message = init_message(Some(Assembly::Grch38), &[&six]);
        assert_eq!(
            message,
            format!(
                "REVEL: column count must be 8 or 9 for REVEL files with GRCh38 positions or 7 \
                 for REVEL files with GRCh37 positions only; {six} has 6 columns"
            )
        );

        let seven = write_tabix_fixture(
            dir.path(),
            "seven.tsv.gz",
            &[
                "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL",
                "21\t100200\tA\tG\tQ\tH\t0.742",
            ],
            IndexSpec::standard(),
        )
        .display()
        .to_string();
        let mut v = variant("21", 100200, b"A", b"G");
        v.transcript_consequences
            .push(missense_row("ENST00000400002", "Q/H"));
        assert_eq!(
            scores(&loaded(Assembly::Grch37, &[&seven]), v),
            [Some("0.742".to_string())]
        );

        let eight = grch38_fixture(
            dir.path(),
            "eight.tsv.gz",
            &[
                "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL",
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742",
            ],
        );
        let mut v = variant("21", 500600, b"A", b"G");
        v.transcript_consequences
            .push(missense_row("ENST00000400002", "Q/H"));
        assert_eq!(
            scores(&loaded(Assembly::Grch38, &[&eight]), v),
            [Some("0.742".to_string())]
        );
    }

    /// The GRCh37 path reads `hg19_pos` (column 2); without an assembly the
    /// plugin refuses to guess (`REVEL.pm:133`).
    #[test]
    fn grch37_reads_the_hg19_position_column() {
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
        )
        .display()
        .to_string();

        let mut v = variant("21", 100200, b"A", b"G");
        v.transcript_consequences
            .push(missense_row("ENST00000400000", "M/V"));
        assert_eq!(
            scores(&loaded(Assembly::Grch37, &[&path]), v),
            [Some("0.742".to_string())]
        );
        assert_eq!(
            init_message(None, &[&path]),
            "REVEL: specify assembly using --assembly GRCh37 or --assembly GRCh38"
        );
    }

    /// `REVEL.pm:138-142` dies when the requested assembly's column is missing
    /// rather than silently reading the other assembly's coordinates.
    #[test]
    fn missing_assembly_column_is_a_hard_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "revel_hg19_only.tsv.gz",
            &[
                "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL",
                "21\t100200\tA\tG\tM\tV\t0.742",
            ],
            IndexSpec::standard(),
        );
        let message = init_message(Some(Assembly::Grch38), &[&path.display().to_string()]);
        assert!(
            message.starts_with("REVEL: assembly is GRCh38, which reads the 'grch38_pos' column")
                && message.contains("has no 'grch38_pos' column"),
            "{message}"
        );
    }

    /// REVEL ships one table carrying both `hg19_pos` (col 2) and `grch38_pos`
    /// (col 3), so the two prepared copies share the header and differ in
    /// which column their `.tbi` indexes. Pointing a GRCh38 run at the
    /// GRCh37-indexed copy passes every header check and annotates nothing at exit
    /// 0; the index's recorded start column is the only signal that catches it.
    #[test]
    fn rejects_a_copy_indexed_for_the_other_assembly() {
        let dir = tempfile::tempdir().unwrap();
        let lines = [
            REVEL_HEADER,
            "21\t9907250\t9068417\tC\tA\tQ\tH\t0.193\tENST00000400754",
        ];
        let g37 = write_tabix_fixture(
            dir.path(),
            "revel_g37.tsv.gz",
            &lines,
            IndexSpec::standard().with_start_col(2),
        )
        .display()
        .to_string();
        let g38 = grch38_fixture(dir.path(), "revel_g38.tsv.gz", &lines);

        let message = init_message(Some(Assembly::Grch38), &[&g37]);
        assert!(
            message.contains("grch38_pos") && message.contains("indexed on 'hg19_pos'"),
            "error must name both the needed and the indexed column, got: {message}"
        );
        let message = init_message(Some(Assembly::Grch37), &[&g38]);
        assert!(
            message.contains("indexed on 'grch38_pos'"),
            "error must name the indexed column, got: {message}"
        );
        loaded(Assembly::Grch37, &[&g37]);
        loaded(Assembly::Grch38, &[&g38]);
    }

    /// `REVEL.pm:116-125` reads the header through `tabix -h`, which prints the
    /// file's leading `#` lines alone, strips one `#`, splits the names on
    /// whitespace and dies `Could not read headers` when none is left: a table
    /// whose header line lacks its `#` (one indexed with `tabix -S 1`) is
    /// refused with the file named, as a table with no header line or a `#`
    /// line naming nothing is, while the same rows under the `#` header score,
    /// a tab between the `#` and the first name being no empty column. A `##`
    /// line alone keeps its second `#` as its one name and is refused on its
    /// column count (`:127-129`).
    #[test]
    fn bare_header_table_is_refused_as_ensembl_vep_refuses_it() {
        let dir = tempfile::tempdir().unwrap();
        let row = "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001";
        let bare_header = REVEL_HEADER
            .strip_prefix('#')
            .expect("REVEL_HEADER starts with #");
        let bare = write_tabix_fixture(
            dir.path(),
            "bare.tsv.gz",
            &[bare_header, row],
            IndexSpec::standard().with_start_col(3).bare_header(),
        )
        .display()
        .to_string();
        assert_eq!(
            init_message(Some(Assembly::Grch38), &[&format!("file={bare}")]),
            format!("REVEL: could not read headers from {bare}")
        );

        let none = grch38_fixture(dir.path(), "none.tsv.gz", &[row]);
        assert_eq!(
            init_message(Some(Assembly::Grch38), &[&format!("file={none}")]),
            format!("REVEL: could not read headers from {none}")
        );

        let nameless = grch38_fixture(dir.path(), "nameless.tsv.gz", &["# \t", row]);
        assert_eq!(
            init_message(Some(Assembly::Grch38), &[&format!("file={nameless}")]),
            format!("REVEL: could not read headers from {nameless}")
        );

        let double = grch38_fixture(dir.path(), "double.tsv.gz", &["##", row]);
        assert_eq!(
            init_message(Some(Assembly::Grch38), &[&format!("file={double}")]),
            format!(
                "REVEL: column count must be 8 or 9 for REVEL files with GRCh38 positions or 7 \
                 for REVEL files with GRCh37 positions only; {double} has 1 columns"
            )
        );

        let hashed = grch38_fixture(dir.path(), "hashed.tsv.gz", &[REVEL_HEADER, row]);
        let padded = grch38_fixture(
            dir.path(),
            "padded.tsv.gz",
            &[&format!("#\t{bare_header}"), row],
        );
        for file in [hashed, padded] {
            let plugin = loaded(Assembly::Grch38, &[&format!("file={file}")]);
            let mut v = variant("21", 500600, b"A", b"G");
            v.transcript_consequences
                .push(missense_row("ENST00000400001", "Q/H"));
            assert_eq!(scores(&plugin, v), [Some("0.742".to_string())]);
        }
    }

    /// The plugin's grammar is positional or keyed, never mixed, and `file` is
    /// mandatory under both readings.
    #[test]
    fn parameters_follow_the_positional_or_keyed_grammar() {
        assert_eq!(
            init_message::<&str>(Some(Assembly::Grch38), &[]),
            "REVEL: no file specified; use file=<path> or give the file first"
        );
        assert_eq!(
            init_message(Some(Assembly::Grch38), &["file=a.tsv.gz", "1"]),
            "REVEL: parameter '1' is not key=value while 'file=a.tsv.gz' is; give every \
             parameter as key=value, or every one by position"
        );
        assert_eq!(
            init_message(
                Some(Assembly::Grch38),
                &["file=a.tsv.gz", "transcript_match=1"]
            ),
            "REVEL: unknown parameter 'transcript_match'; REVEL accepts file, no_match"
        );
    }

    /// `REVEL.pm:207` splits the line htslib hands over, its line end alone
    /// removed, on `\t`, and `:223` copies the last cell of an eight-column
    /// row into `REVEL` as it stands: `0.176 ` is written with its trailing
    /// space, and a cell holding a single space as that space. In a
    /// nine-column row the last cell is the transcript list (`:192-193`), each
    /// id compared as `split /;/` leaves it: one holding a single space is true
    /// and names no transcript, and `ENST00000400001 ` is not the row's id, so
    /// neither row scores.
    #[test]
    fn a_trailing_space_on_the_last_cell_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let eight = grch38_fixture(
            dir.path(),
            "eight.tsv.gz",
            &[
                "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL",
                "21\t100200\t500600\tA\tG\tQ\tH\t0.176 ",
                "21\t100300\t500700\tA\tG\tQ\tH\t ",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&eight]);
        let scored = |pos: u64| {
            let mut v = variant("21", pos, b"A", b"G");
            v.transcript_consequences
                .push(missense_row("ENST00000400001", "Q/H"));
            v
        };
        assert_eq!(
            scores(&plugin, scored(500600)),
            [Some("0.176 ".to_string())]
        );
        assert_eq!(scores(&plugin, scored(500700)), [Some(" ".to_string())]);

        let nine = grch38_fixture(
            dir.path(),
            "nine.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\tQ\tH\t0.176\tENST00000400001 ",
                "21\t100300\t500700\tA\tG\tQ\tH\t0.020\t ",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&nine]);
        assert_eq!(scores(&plugin, scored(500600)), [None]);
        assert_eq!(scores(&plugin, scored(500700)), [None]);
    }

    /// The text of the `Run` error a batch run of `variants` fails with.
    fn failed_batch(plugin: &RevelPlugin, variants: &mut [InputVariant]) -> String {
        match plugin.run_batch(variants) {
            Err(PluginError::Run(message)) => message,
            other => panic!("expected PluginError::Run, got {other:?}"),
        }
    }

    /// The `REVEL` field of every row of `variant`.
    fn fields(variant: &InputVariant) -> Vec<Option<&str>> {
        variant
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.get(FIELD).map(String::as_str))
            .collect()
    }

    /// `REVEL.pm:170-188` puts every record of the variant's span to
    /// `get_matched_variant_alleles` in file order, and `trim_sequences` dies
    /// on a reference or alternate cell `split /\t/` left undefined, a cell at
    /// or past the row's last filled one (Utils/Sequence.pm:968, from `:1230`),
    /// before any later record is read: a missense row whose walk meets such a
    /// record before its first full match takes no field, and the batch
    /// reports the failure after annotating every other row. Each row's `run`
    /// dies on its own, so a sibling transcript's row matched before the
    /// record keeps its score while the row whose match lies after it takes
    /// none, and a variant at a clean position in the same batch keeps its
    /// fields.
    #[test]
    fn a_row_with_no_allele_cell_before_the_match_fails_the_rows_that_reach_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA",
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001",
                "21\t100300\t500700\tA\tG\tQ\tH\t0.311\tENST00000400002",
                "21\t100400\t500800\tA\tG\tQ\tH\t0.500\tENST00000400003",
                "21\t100400\t500800\tA",
                "21\t100400\t500800\tA\tG\tQ\tH\t0.600\tENST00000400004",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&format!("file={path}")]);
        let mut dying = variant("21", 500600, b"A", b"G");
        dying
            .transcript_consequences
            .push(missense_row("ENST00000400001", "Q/H"));
        let mut clean = variant("21", 500700, b"A", b"G");
        clean
            .transcript_consequences
            .push(missense_row("ENST00000400002", "Q/H"));
        let mut split = variant("21", 500800, b"A", b"G");
        split
            .transcript_consequences
            .push(missense_row("ENST00000400003", "Q/H"));
        split
            .transcript_consequences
            .push(missense_row("ENST00000400004", "Q/H"));
        let mut variants = vec![dying, clean, split];

        let message = failed_batch(&plugin, &mut variants);
        assert!(
            message.contains("has no reference or alternate allele cell")
                && message.contains("21:500600"),
            "{message}"
        );
        assert_eq!(fields(&variants[0]), [None]);
        assert_eq!(fields(&variants[1]), [Some("0.311")]);
        assert_eq!(fields(&variants[2]), [Some("0.500"), None]);
    }

    /// `REVEL.pm:193` and `:197` return at the first full match, so a record
    /// after it is never put to the match: a row with no allele cell after the
    /// matching record dies nowhere when every missense row of the variant
    /// matches before it.
    #[test]
    fn a_row_with_no_allele_cell_after_the_match_is_never_reached() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001;ENST00000400002",
                "21\t100200\t500600\tA",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&path]);
        let mut v = variant("21", 500600, b"A", b"G");
        v.transcript_consequences
            .push(missense_row("ENST00000400001", "Q/H"));
        v.transcript_consequences
            .push(missense_row("ENST00000400002", "Q/H"));
        assert_eq!(
            scores(&plugin, v),
            [Some("0.742".to_string()), Some("0.742".to_string())]
        );
    }

    /// An empty cell before a filled one is the defined empty string to
    /// `split /\t/`, which `trim_sequences` compares (Utils/Sequence.pm:968
    /// dies on the undefined value alone): a record whose `alt` cell is empty
    /// matches no allele and dies nowhere, and the record after it scores.
    #[test]
    fn an_empty_alt_cell_before_a_filled_cell_matches_nothing_and_dies_nowhere() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\t\tQ\tH\t0.500\tENST00000400001",
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&path]);
        let mut v = variant("21", 500600, b"A", b"G");
        v.transcript_consequences
            .push(missense_row("ENST00000400001", "Q/H"));
        assert_eq!(scores(&plugin, v), [Some("0.742".to_string())]);
    }

    /// `REVEL.pm:162` returns before the loop for a row that is not missense,
    /// so a variant with no missense row never puts a record to the match: a
    /// record with no allele cell at its position costs it nothing.
    #[test]
    fn a_row_with_no_allele_cell_is_harmless_to_a_variant_with_no_missense_row() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[REVEL_HEADER, "21\t100200\t500600\tA"],
        );
        let plugin = loaded(Assembly::Grch38, &[&path]);
        let mut v = variant("21", 500600, b"A", b"G");
        v.transcript_consequences.push(TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: "ENST00000400001".into(),
                ..Transcript::default()
            })),
            consequences: smallvec![Consequence::SynonymousVariant],
            amino_acids: Some("Q".to_string()),
            ..Default::default()
        });
        assert_eq!(scores(&plugin, v), [None]);
    }

    /// `REVEL.pm:197` returns the matching record's result whatever its score
    /// cell holds (`:223`): a record of seven cells matches on its alleles and
    /// peptides alone (its transcript cell undefined, `:191`) and writes the
    /// undefined score as the empty value, which the writers print as Ensembl
    /// VEP prints an undefined value, and the full record after it is never
    /// read; a score cell that is empty and last is undefined the same way. A
    /// score cell empty before a filled transcript cell is the defined empty
    /// string, which the row names as a string (`Tab.pm:128` and
    /// `OutputFactory/JSON.pm:327` print the two apart).
    #[test]
    fn a_matching_row_without_a_score_cell_writes_the_empty_value() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "revel.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\tQ\tH",
                "21\t100200\t500600\tA\tG\tQ\tH\t0.900\tENST00000400001",
                "21\t100300\t500700\tA\tG\tQ\tH\t",
                "21\t100400\t500800\tA\tG\tQ\tH\t\tENST00000400001",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[&path]);
        let scored = |pos: u64| {
            let mut v = variant("21", pos, b"A", b"G");
            v.transcript_consequences
                .push(missense_row("ENST00000400001", "Q/H"));
            v
        };
        for pos in [500600, 500700] {
            let mut variants = vec![scored(pos)];
            plugin.run_batch(&mut variants).expect("run_batch");
            assert_eq!(fields(&variants[0]), [Some("")]);
            assert!(variants[0].transcript_consequences[0]
                .plugin_empty_strings
                .is_empty());
        }
        let mut variants = vec![scored(500800)];
        plugin.run_batch(&mut variants).expect("run_batch");
        assert_eq!(fields(&variants[0]), [Some("")]);
        assert_eq!(
            variants[0].transcript_consequences[0].plugin_empty_strings,
            [FIELD.to_string()]
        );
    }

    /// A row's position cell for the run's assembly reaches
    /// `get_matched_variant_alleles` as text (`$revel_start`, REVEL.pm:172) and
    /// is keyed as text (ensembl-variation Utils/Sequence.pm:1232), a one-base
    /// pair losing no leading base: `500600` scores the variant at 500600,
    /// `500600.0` and `500600.5` do not, and `500601x` does not score the
    /// variant at 500601. The tabix index places each at its leading integer,
    /// so every row is returned and the text decides.
    #[test]
    fn a_position_cell_is_matched_as_the_text_the_key_holds() {
        let dir = tempfile::tempdir().unwrap();
        let path = grch38_fixture(
            dir.path(),
            "cells.tsv.gz",
            &[
                REVEL_HEADER,
                "21\t100200\t500600\tA\tG\tQ\tH\t0.742\tENST00000400001",
                "21\t100200\t500600.0\tA\tT\tQ\tL\t0.311\tENST00000400001",
                "21\t100200\t500600.5\tA\tC\tQ\tP\t0.512\tENST00000400001",
                "21\t100201\t500601x\tA\tG\tQ\tH\t0.613\tENST00000400001",
            ],
        );
        let plugin = loaded(Assembly::Grch38, &[format!("file={path}")]);
        let score = |pos: u64, alt: &[u8], change: &str| {
            let mut v = variant("21", pos, b"A", alt);
            v.transcript_consequences
                .push(missense_row("ENST00000400001", change));
            scores(&plugin, v)
        };
        assert_eq!(score(500600, b"G", "Q/H"), [Some("0.742".to_string())]);
        assert_eq!(score(500600, b"T", "Q/L"), [None]);
        assert_eq!(score(500600, b"C", "Q/P"), [None]);
        assert_eq!(score(500601, b"G", "Q/H"), [None]);
    }
}
