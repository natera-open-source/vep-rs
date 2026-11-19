// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! TranscriptAnnotator plugin: the columns of a tabix-indexed table keyed by
//! variant and transcript, appended to the transcript rows that match.
//!
//! The table's first five columns are the chromosome, position, reference
//! allele, alternate allele and transcript id; every later column is an output
//! field named `<prefix>_<column>`, the prefix the file's basename without its
//! extensions unless `prefix=` names another or `prefix=0` drops it
//! (TranscriptAnnotator.pm:158-186). The column names come from the last header
//! line (`:95-124`); `cols=` selects a colon-separated subset of them (`:126-156`),
//! a value the truth test `if ( $param_hash->{cols} )` (`:132`) reads as false
//! (`0` or empty) selecting every column as an absent `cols=` does (`:153`).
//! A row annotates a transcript row when its transcript id is the transcript's
//! stable id or one of its RefSeq cross-references (`:210-211`) and its alleles
//! match the variant's (`:224-236`); the first such row of the region
//! `start-2..end` wins (`:218-237`). Every row of the region before that one
//! is put to the allele match, which dies on a row whose reference or
//! alternate cell is undefined (`:245-247`; ensembl-variation
//! Utils/Sequence.pm:968): a transcript row whose scan reaches such a row
//! keeps none of the plugin's fields, while one matched before it keeps its
//! own, as each output line's `run` dies on its own.

use std::borrow::Cow;
use std::path::Path;

use indexmap::IndexMap;
use tracing::{debug, warn};

use crate::annotation_store::{open_best_store_with, AnnotationStore};
use crate::params::PluginParams;
use crate::tabix::{
    allele_match, split_field_count, table_contig, BatchQueryResult, TabixAnnotator, TabixConfig,
    TabixRecord,
};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::assembly::Assembly;
use vep_core::consequence::TranscriptConsequence;
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

const NAME: &str = "TranscriptAnnotator";

/// The fixed leading columns: chromosome, position, reference and alternate
/// allele, transcript id (`TranscriptAnnotator.pm:247`).
const POS_COL: usize = 1;
const REF_COL: usize = 2;
const ALT_COL: usize = 3;
const TRANSCRIPT_COL: usize = 4;
const FIRST_VALUE_COL: usize = 5;

/// What the first five column names are expected to start with, checked
/// case-insensitively and warned about when no column does (`:113-119`).
const EXPECTED_NAME_PREFIXES: [&str; 5] = ["chr", "pos", "ref", "alt", "tr"];

/// Bases read before the variant's start, as the plugin widens its query
/// (`:218`) so a row written with the base before an indel is found.
const WINDOW_BEFORE: u64 = 2;

/// The characters `trim=1` removes: Ensembl VEP's `_trim_whitespaces` runs
/// `\s` on the file's bytes (`:91`), which matches the ASCII whitespace alone,
/// so a Unicode space such as U+00A0 stays part of the cell.
const TRIMMED: [char; 6] = [' ', '\t', '\n', '\x0B', '\x0C', '\r'];

/// One appended column of the file: its output key and whether `cols=` kept it.
#[derive(Debug, Clone)]
struct Column {
    key: String,
    selected: bool,
}

/// TranscriptAnnotator plugin.
pub struct TranscriptAnnotatorPlugin {
    store: Option<Box<dyn AnnotationStore>>,
    /// The data file as given, the text every header description repeats.
    file: String,
    /// Whether every header name and cell is trimmed of surrounding whitespace:
    /// `trim=` is read as the truth test `if $self->{trim}` reads it (`:110`,
    /// `:246`), so only `0` and the empty value turn it off, and it is on when
    /// absent (`:164`).
    trim: bool,
    /// The appended columns in header order, from the sixth column on.
    columns: Vec<Column>,
    /// The run's assembly, which names a mitochondrial record's contig to the
    /// table ([`table_contig`]).
    assembly: Option<Assembly>,
}

impl Default for TranscriptAnnotatorPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl TranscriptAnnotatorPlugin {
    /// Create an uninitialized plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            file: String::new(),
            trim: true,
            columns: Vec::new(),
            assembly: None,
        }
    }

    /// `text` with surrounding whitespace removed under `trim=1`, as given
    /// under `trim=0` (`:89-93`, `:110`, `:246`).
    fn cell<'a>(&self, text: &'a str) -> &'a str {
        if self.trim {
            text.trim_matches(TRIMMED)
        } else {
            text
        }
    }

    /// `record` with every cell trimmed when `trim=1` and any cell carries
    /// surrounding whitespace; the record itself otherwise.
    fn trimmed<'a>(&self, record: &'a TabixRecord) -> Cow<'a, TabixRecord> {
        if !self.trim
            || record
                .columns
                .iter()
                .all(|c| c.trim_matches(TRIMMED).len() == c.len())
        {
            return Cow::Borrowed(record);
        }
        let mut trimmed = record.clone();
        for cell in &mut trimmed.columns {
            *cell = cell.trim_matches(TRIMMED).to_string();
        }
        Cow::Owned(trimmed)
    }

    /// The region read for a variant: from two bases before its span to its
    /// end (`:218`), the span taken in either orientation since an insertion's
    /// end precedes its start.
    fn window(variant: &InputVariant) -> (u64, u64) {
        let low = variant.start.min(variant.end);
        let high = variant.start.max(variant.end);
        (low.saturating_sub(WINDOW_BEFORE).max(1), high)
    }

    /// The transcript ids a row may name for `consequence`: the stable id and
    /// the comma-separated RefSeq cross-references, none when they are `-`
    /// (`:210-211`).
    fn transcript_ids(consequence: &TranscriptConsequence) -> Vec<&str> {
        let mut ids = vec![consequence.transcript_id()];
        if let Some(refseq) = consequence.refseq().filter(|r| *r != "-") {
            ids.extend(refseq.split(',').filter(|id| !id.is_empty()));
        }
        ids
    }

    /// The selected columns' values of the first row of `records` that names
    /// one of `ids` and matches the variant's alleles, both sides trimmed of
    /// their shared bases (`:221-238`); `Ok(None)` when no row does. Ensembl
    /// VEP puts every row to `get_matched_variant_alleles` before it reads the
    /// match and the transcript test together (`:222-237`), and the match dies
    /// on a row whose reference or alternate cell `split /\t/` left undefined
    /// (`:245-247`; ensembl-variation Utils/Sequence.pm:1230, :968), whatever
    /// the row's transcript cell: such a row before the first matching row is
    /// the transcript row's failure, and one after it is never reached. The die
    /// is confined to the rows Ensembl VEP reads, two bases before the variant's
    /// start to its end (`:218`): the rows of the wider window an insertion is
    /// read with here are matched but never fail a row.
    fn matching_values(
        &self,
        records: &[&TabixRecord],
        variant: &InputVariant,
        ids: &[&str],
    ) -> Result<Option<RowValues>, PluginError> {
        let read_by_ensembl_vep = variant.start.saturating_sub(WINDOW_BEFORE).max(1)..=variant.end;
        for record in records {
            // The cells are counted before trimming: a whitespace-only trailing
            // cell is a field `split /\t/` keeps and `trim=1` then empties
            // (`:245-246`), compared as `-`.
            let present = split_field_count(&record.columns);
            let trimmed = self.trimmed(record);
            // A Perl-false position cell, read after the trim, returns the match
            // before an undefined allele dies (Sequence.pm:1173 before :1230).
            if read_by_ensembl_vep.contains(&record.start)
                && !allele_match::perl_false_position(&trimmed, POS_COL)
                && allele_match::defined_alleles(record, REF_COL, ALT_COL).is_none()
            {
                return Err(allele_match::undefined_allele_row(NAME, &self.file, record));
            }
            let record = trimmed;
            let same_transcript = record
                .columns
                .get(TRANSCRIPT_COL)
                .is_some_and(|id| ids.contains(&id.as_str()));
            if !same_transcript
                || !allele_match::matched_variant_alleles(
                    &record, variant, POS_COL, REF_COL, ALT_COL,
                )
            {
                continue;
            }
            // A column past the row's fields is undefined (`:250`), which
            // Ensembl VEP writes as an empty value; a field the row has and
            // leaves empty, or trims to nothing, is the empty string.
            let mut values = Vec::new();
            let mut empty_strings = Vec::new();
            for (i, column) in self.columns.iter().enumerate() {
                if !column.selected {
                    continue;
                }
                let index = FIRST_VALUE_COL + i;
                let value = record.columns.get(index).map_or("", String::as_str);
                if value.is_empty() && index < present {
                    empty_strings.push(column.key.clone());
                }
                values.push((column.key.clone(), value.to_string()));
            }
            return Ok(Some(RowValues {
                values,
                empty_strings,
            }));
        }
        Ok(None)
    }
}

/// The selected columns' values of a matched row, and the keys among them
/// whose value is Ensembl VEP's empty string rather than its undefined value.
struct RowValues {
    values: Vec<(String, String)>,
    empty_strings: Vec<String>,
}

/// The file's basename without its extensions, the default prefix (`:83-87`):
/// the basename cut at the first `.` whose remainder holds no other `.`, or
/// holds one only as the `.gz` of a `.<ext>.gz` tail, so `scores.tsv.gz` and
/// `scores.gz` both give `scores` and `calls.vcf.bgz` gives `calls.vcf`.
fn basename_without_extensions(file: &str) -> String {
    let base = Path::new(file)
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    for (dot, _) in base.match_indices('.') {
        let rest = &base[dot + 1..];
        let plain = !rest.contains('.');
        let gz_tail = rest
            .strip_suffix(".gz")
            .is_some_and(|extension| !extension.contains('.'));
        if plain || gz_tail {
            return base[..dot].to_string();
        }
    }
    base
}

impl BuiltinPlugin for TranscriptAnnotatorPlugin {
    fn name(&self) -> &str {
        NAME
    }

    fn header_info(&self) -> Vec<(String, String)> {
        self.columns
            .iter()
            .filter(|column| column.selected)
            .map(|column| (column.key.clone(), format!("column from {}", self.file)))
            .collect()
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(NAME, &["file", "cols", "prefix", "trim"], params)?;
        if !params.positional().is_empty() {
            return Err(PluginError::Init(format!(
                "{NAME}: parameters are key=value: file=path/to/file.txt.gz[,cols=a:b][,prefix=text][,trim=0|1]"
            )));
        }
        self.trim = params.get("trim").is_none_or(perl_true);
        let Some(file) = params.get("file").filter(|f| !f.is_empty()) else {
            return Err(PluginError::Init(format!(
                "{NAME}: no file specified; give file=path/to/file.txt.gz"
            )));
        };
        if !Path::new(file).is_file() {
            return Err(PluginError::Init(format!(
                "{NAME}: data file {file} not found"
            )));
        }
        self.file = file.to_string();

        let config = TabixConfig {
            file_path: file.into(),
            chr_col: 0,
            start_col: POS_COL,
            end_col: None,
            ref_col: Some(REF_COL),
            alt_col: Some(ALT_COL),
            zero_based: false,
        };
        // The rows keep their trailing cells, which Ensembl VEP's `split /\t/`
        // counts when they hold whitespace (`:245-246`).
        let store = open_best_store_with(
            Path::new(file),
            config,
            self.assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )?;

        let Some(header) = store.header() else {
            return Err(PluginError::Init(format!(
                "{NAME}: {file} has no header line naming its columns"
            )));
        };
        let header = &header[..split_field_count(header)];
        let names: Vec<&str> = header.iter().map(|name| self.cell(name)).collect();
        for (i, prefix) in EXPECTED_NAME_PREFIXES.iter().enumerate() {
            if !names
                .iter()
                .any(|name| name.to_ascii_lowercase().starts_with(prefix))
            {
                warn!(
                    "{NAME}: Name of column {i} ('{}') should start with '{prefix}'",
                    names.get(i).copied().unwrap_or("")
                );
            }
        }
        let appended: Vec<&str> = names.get(FIRST_VALUE_COL..).unwrap_or(&[]).to_vec();

        let selected: Vec<&str> = match params.get("cols").filter(|cols| perl_true(cols)) {
            None => appended.clone(),
            Some(cols) => {
                let (valid, invalid): (Vec<&str>, Vec<&str>) =
                    cols.split(':').partition(|name| appended.contains(name));
                let valid_names = format!("Valid column names: {}", appended.join(", "));
                if valid.is_empty() {
                    return Err(PluginError::Init(format!(
                        "{NAME}: no valid columns selected. {valid_names}"
                    )));
                }
                if !invalid.is_empty() {
                    warn!(
                        "{NAME}: columns {} not found in {file}. {valid_names}",
                        invalid.join(", ")
                    );
                }
                valid
            }
        };

        let prefix = match params.get("prefix") {
            Some("0") => String::new(),
            Some(prefix) if !prefix.is_empty() => format!("{prefix}_"),
            _ => format!("{}_", basename_without_extensions(file)),
        };
        self.columns = appended
            .iter()
            .map(|name| Column {
                key: format!("{prefix}{name}"),
                selected: selected.contains(name),
            })
            .collect();

        debug!(
            file,
            columns = self.columns.len(),
            selected = selected.len(),
            "TranscriptAnnotator plugin initialized"
        );
        self.store = Some(store);
        Ok(())
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
            .filter(|v| self.annotates_variant(v))
            .map(|v| {
                let (low, high) = Self::window(v);
                (table_contig(v, self.assembly), low, high)
            })
            .collect();
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
        for variant in variants.iter_mut().filter(|v| self.annotates_variant(v)) {
            let (low, high) = Self::window(variant);
            let records = batch.lookup(table_contig(variant, self.assembly), low, high);
            if records.is_empty() {
                continue;
            }
            let mut matched: Vec<(usize, RowValues)> = Vec::new();
            for (i, consequence) in variant.transcript_consequences.iter().enumerate() {
                if !self.annotates_row(consequence) {
                    continue;
                }
                let ids = Self::transcript_ids(consequence);
                match self.matching_values(&records, variant, &ids) {
                    Ok(Some(row)) => matched.push((i, row)),
                    Ok(None) => {}
                    Err(error) => {
                        failure.get_or_insert(error);
                    }
                }
            }
            for (i, row) in matched {
                let consequence = &mut variant.transcript_consequences[i];
                // A later plugin's value replaces an earlier one under the same
                // key (OutputFactory.pm:2230-2231), its empty-string mark included.
                consequence
                    .plugin_empty_strings
                    .retain(|key| !row.values.iter().any(|(written, _)| written == key));
                consequence.plugin_data.extend(row.values);
                consequence.plugin_empty_strings.extend(row.empty_strings);
            }
        }
        match failure {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }

    /// Ensembl VEP runs the module once per output line and catches each
    /// line's failure on its own (`run`, `:205-240`; OutputFactory.pm:2225-2245),
    /// so a transcript row whose matching row precedes the row the match dies
    /// on keeps the fields while the rows whose scan reaches it lose them.
    fn annotate_fails_per_row(&self) -> bool {
        true
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
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use std::path::PathBuf;
    use std::sync::Arc;
    use vep_core::consequence::FeatureType;
    use vep_core::transcript::Transcript;

    const HEADER: &str = "#Chrom\tPos\tRef\tAlt\tTranscript\tSIFT_score\tSIFT_pred\tComment";

    /// The table of the module documentation's shape: two SNVs on Ensembl
    /// transcripts, a RefSeq-keyed row, a row written with the base before a
    /// deletion, and cells with surrounding whitespace.
    fn fixture(dir: &Path) -> PathBuf {
        let rows = [
            "## an annotation table",
            HEADER,
            "21\t436154\tA\tG\tNM_001347882.2\t0.03\tDeleterious\tBad",
            "21\t1887470\tCT\tC\tENST00000421485\t0.5\tTolerated\tDeleted",
            "21\t1887471\tC\tT\tENST00000421485\t0.86\t Tolerated \tGood",
            "21\t1887471\tC\tT\tENST00000999999\t0.11\tDeleterious\tOther",
            "21\t1887471\tC\tA\tENST00000421485\t0.42\tTolerated\tSecond allele",
        ];
        write_tabix_fixture(dir, "scores.tsv.gz", &rows, IndexSpec::standard())
    }

    /// A table whose header line ends in a tab and whose rows carry the
    /// empty-cell shapes: two interior empty cells, a trailing zero-length
    /// cell, a whitespace-only interior cell, a whitespace-only trailing
    /// cell, and a cell holding one no-break space (U+00A0).
    fn empty_cells_fixture(dir: &Path) -> PathBuf {
        let header = format!("{HEADER}\t");
        let rows = [
            header.as_str(),
            "21\t3000001\tC\tT\tENST00000421485\t\t\tEmptyMiddle",
            "21\t3000002\tC\tT\tENST00000421485\t0.44\tTolerated\t",
            "21\t3000003\tC\tT\tENST00000421485\t \tTolerated\tSpaceMiddle",
            "21\t3000004\tC\tT\tENST00000421485\t0.86\tTolerated\t  ",
            "21\t3000005\tC\tT\tENST00000421485\t0.03\t\u{a0}\tNbsp",
        ];
        write_tabix_fixture(dir, "cells.tsv.gz", &rows, IndexSpec::standard())
    }

    /// A table whose positions each hold a valid row and a row short of some
    /// cell, ten bases apart so no window reads two positions: a three-cell row
    /// before the valid row, a row of empty trailing cells after the valid row,
    /// an empty alternate cell before a filled transcript cell, a four-cell row
    /// with alleles but no transcript cell, a whitespace-only trailing
    /// alternate cell, and a three-cell row between the matching rows of two
    /// transcripts.
    fn short_rows_fixture(dir: &Path) -> PathBuf {
        let rows = [
            HEADER,
            "21\t5000001\tC",
            "21\t5000001\tC\tT\tENST00000421485\t0.5\tTolerated\tAfterShort",
            "21\t5000011\tC\tT\tENST00000421485\t0.6\tTolerated\tClean",
            "21\t5000021\tC\tT\tENST00000421485\t0.7\tTolerated\tBeforeShort",
            "21\t5000021\tC\t\t\t",
            "21\t5000031\tC\t\tENST00000421485\t0.8\tTolerated\tEmptyAlt",
            "21\t5000031\tC\tT\tENST00000421485\t0.9\tTolerated\tAfterEmptyAlt",
            "21\t5000041\tC\tT",
            "21\t5000041\tC\tT\tENST00000421485\t1.0\tTolerated\tAfterNoTranscript",
            "21\t5000051\tC\t ",
            "21\t5000051\tC\tT\tENST00000421485\t0.2\tTolerated\tAfterSpaceAlt",
            "21\t5000061\tC\tT\tENST00000421485\t0.3\tTolerated\tFirstTranscript",
            "21\t5000061\tC",
            "21\t5000061\tC\tT\tENST00000999999\t0.4\tTolerated\tSecondTranscript",
        ];
        write_tabix_fixture(dir, "short.tsv.gz", &rows, IndexSpec::standard())
    }

    fn init_table(file: &Path, extra: &[&str]) -> Result<TranscriptAnnotatorPlugin, PluginError> {
        let mut params = vec![format!("file={}", file.display())];
        params.extend(extra.iter().map(|p| p.to_string()));
        let mut plugin = TranscriptAnnotatorPlugin::new();
        plugin.init(&params)?;
        Ok(plugin)
    }

    fn init_with(dir: &Path, extra: &[&str]) -> Result<TranscriptAnnotatorPlugin, PluginError> {
        init_table(&fixture(dir), extra)
    }

    /// The `Init` message `result` failed with; any other outcome fails the test.
    fn expect_init_error<T>(result: Result<T, PluginError>) -> String {
        match result {
            Err(PluginError::Init(message)) => message,
            Ok(_) => panic!("expected PluginError::Init"),
            Err(other) => panic!("expected PluginError::Init, got {other:?}"),
        }
    }

    /// The `Init` message of initializing a fresh plugin with `params`.
    fn init_error(params: &[&str]) -> String {
        let params: Vec<String> = params.iter().map(|p| p.to_string()).collect();
        expect_init_error(TranscriptAnnotatorPlugin::new().init(&params))
    }

    /// The `Init` message of initializing the fixture with `extra` parameters.
    fn init_message(dir: &Path, extra: &[&str]) -> String {
        expect_init_error(init_with(dir, extra))
    }

    fn transcript_row(id: &str, refseq: Option<&str>) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: id.into(),
                refseq: refseq.map(str::to_string),
                ..Transcript::default()
            })),
            feature_type: FeatureType::Transcript,
            ..TranscriptConsequence::default()
        }
    }

    fn snv(pos: u64, alt: &str) -> InputVariant {
        InputVariant::new(
            "21".into(),
            pos,
            pos,
            b"C".to_vec(),
            alt.as_bytes().to_vec(),
        )
    }

    fn keys(plugin: &TranscriptAnnotatorPlugin) -> Vec<String> {
        plugin.header_info().into_iter().map(|(k, _)| k).collect()
    }

    fn strings(texts: &[&str]) -> Vec<String> {
        texts.iter().map(|text| text.to_string()).collect()
    }

    #[test]
    fn file_is_mandatory() {
        assert_eq!(
            init_error(&[]),
            "TranscriptAnnotator: no file specified; give file=path/to/file.txt.gz"
        );
        let without_file = init_error(&["cols=SIFT_score"]);
        assert!(
            without_file.starts_with("TranscriptAnnotator: no file specified"),
            "{without_file}"
        );
        assert_eq!(
            init_error(&["file=/no/such/table.tsv.gz"]),
            "TranscriptAnnotator: data file /no/such/table.tsv.gz not found"
        );
        let positional = init_error(&["/some/table.tsv.gz"]);
        assert!(
            positional.starts_with("TranscriptAnnotator: parameters are key=value"),
            "{positional}"
        );
    }

    /// The column names are the last `#` line's (`:102-105`), so the preamble
    /// line above it names nothing, and the five leading columns are not
    /// output fields (`:122`).
    #[test]
    fn columns_come_from_the_last_header_line() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &[]).unwrap();
        assert_eq!(
            keys(&plugin),
            ["scores_SIFT_score", "scores_SIFT_pred", "scores_Comment"]
        );
        let description = &plugin.header_info()[0].1;
        assert_eq!(
            description,
            &format!("column from {}", dir.path().join("scores.tsv.gz").display())
        );
    }

    #[test]
    fn default_prefix_is_the_basename_without_extensions() {
        assert_eq!(
            basename_without_extensions("/a/b/transcript_annotator.tsv.gz"),
            "transcript_annotator"
        );
        assert_eq!(basename_without_extensions("file.txt.gz"), "file");
        assert_eq!(basename_without_extensions("file.gz"), "file");
        assert_eq!(basename_without_extensions("file.bgz"), "file");
        assert_eq!(basename_without_extensions("calls.vcf.bgz"), "calls.vcf");
        assert_eq!(basename_without_extensions("a.b.c"), "a.b");
        assert_eq!(basename_without_extensions("x.tar.gz"), "x");
        assert_eq!(basename_without_extensions("plain"), "plain");
        assert_eq!(basename_without_extensions("dir.d/plain"), "plain");

        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["prefix=SIFT"]).unwrap();
        assert_eq!(
            keys(&plugin),
            ["SIFT_SIFT_score", "SIFT_SIFT_pred", "SIFT_Comment"]
        );
        // An empty prefix falls back to the basename (`:178`).
        let plugin = init_with(dir.path(), &["prefix="]).unwrap();
        assert_eq!(keys(&plugin)[0], "scores_SIFT_score");
    }

    #[test]
    fn prefix_0_disables_the_prefix() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["prefix=0"]).unwrap();
        assert_eq!(keys(&plugin), ["SIFT_score", "SIFT_pred", "Comment"]);
    }

    /// The selected columns keep the header's order whatever order `cols=`
    /// names them in, and a name the header lacks is warned about and skipped
    /// (`:135-151`).
    #[test]
    fn cols_subset_keeps_header_order_and_unknown_names_warn() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &["cols=Comment:Missing:SIFT_score"]).unwrap();
        assert_eq!(keys(&plugin), ["scores_SIFT_score", "scores_Comment"]);
        let selected: Vec<bool> = plugin.columns.iter().map(|c| c.selected).collect();
        assert_eq!(selected, [true, false, true]);
    }

    #[test]
    fn no_valid_cols_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let message = init_message(dir.path(), &["cols=Missing:Absent"]);
        assert_eq!(
            message,
            "TranscriptAnnotator: no valid columns selected. Valid column names: SIFT_score, SIFT_pred, Comment"
        );
    }

    /// `cols` is gated by the truth test `if ( $param_hash->{cols} )` (`:132`),
    /// so `cols=0` and `cols=` select every column as an absent `cols` does
    /// (`:153`), while a true value naming no column is refused (`:146`).
    #[test]
    fn cols_0_selects_every_column() {
        let dir = tempfile::tempdir().unwrap();
        let all = ["scores_SIFT_score", "scores_SIFT_pred", "scores_Comment"];
        assert_eq!(keys(&init_with(dir.path(), &["cols=0"]).unwrap()), all);
        assert_eq!(keys(&init_with(dir.path(), &["cols="]).unwrap()), all);
        // `0:Comment` is a true list: `0` is an unknown name, warned about.
        assert_eq!(
            keys(&init_with(dir.path(), &["cols=0:Comment"]).unwrap()),
            ["scores_Comment"]
        );
        let message = init_message(dir.path(), &["cols=foo"]);
        assert!(
            message.starts_with("TranscriptAnnotator: no valid columns selected."),
            "{message}"
        );
    }

    /// A row keyed on a RefSeq accession annotates the Ensembl transcript that
    /// carries the accession among its cross-references (`:210-211`).
    #[test]
    fn refseq_xref_ids_match_as_transcript_ids() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &[]).unwrap();
        let mut variants = vec![InputVariant::new(
            "21".into(),
            436154,
            436154,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        variants[0].transcript_consequences.extend([
            transcript_row("ENST00000300305", Some("NM_001347882.2,NM_000001.1")),
            transcript_row("ENST00000300306", Some("-")),
            transcript_row("ENST00000300307", None),
        ]);
        plugin.run_batch(&mut variants).unwrap();
        let scores: Vec<Option<&str>> = variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.get("scores_SIFT_score").map(String::as_str))
            .collect();
        assert_eq!(scores, [Some("0.03"), None, None]);
        assert_eq!(
            variants[0].transcript_consequences[0]
                .plugin_data
                .get("scores_Comment")
                .map(String::as_str),
            Some("Bad")
        );
        assert!(variants[0].plugin_data.is_empty());
    }

    /// Rows are matched through `allele_match::matched_variant_alleles`: the
    /// row with this allele at this position on this transcript wins over a
    /// row of another allele or another transcript, and a deletion written with
    /// its leading base matches the variant's trimmed form.
    #[test]
    fn allele_match_uses_the_shared_matcher() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &[]).unwrap();
        let mut variants = vec![
            snv(1887471, "T"),
            snv(1887471, "A"),
            snv(1887471, "G"),
            InputVariant::new("21".into(), 1887471, 1887471, b"T".to_vec(), b"-".to_vec()),
        ];
        for variant in &mut variants {
            variant
                .transcript_consequences
                .push(transcript_row("ENST00000421485", None));
            variant
                .transcript_consequences
                .push(transcript_row("ENST00000999999", None));
        }
        plugin.run_batch(&mut variants).unwrap();
        let comments: Vec<Vec<Option<&str>>> = variants
            .iter()
            .map(|v| {
                v.transcript_consequences
                    .iter()
                    .map(|tc| tc.plugin_data.get("scores_Comment").map(String::as_str))
                    .collect()
            })
            .collect();
        assert_eq!(
            comments,
            [
                vec![Some("Good"), Some("Other")],
                vec![Some("Second allele"), None],
                vec![None, None],
                vec![Some("Deleted"), None],
            ]
        );
    }

    /// Ensembl VEP 116.2 trims the record's alleles as well as the row's inside
    /// `get_matched_variant_alleles` (ensembl-variation Utils/Sequence.pm:1199-1208
    /// and :1223-1232), so a substitution written with an unchanged leading
    /// base (`GC/GT` at 1887470) or trailing base (`CA/TA` at 1887471) takes
    /// the rows of `C T` at 1887471; a substitution of another base takes none.
    #[test]
    fn a_substitution_with_an_unchanged_flanking_base_matches_its_single_base_row() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &[]).unwrap();
        let two_bases = |start: u64, reference: &[u8], alt: &[u8]| {
            InputVariant::new(
                "21".into(),
                start,
                start + 1,
                reference.to_vec(),
                alt.to_vec(),
            )
        };
        let mut variants = vec![
            two_bases(1887470, b"GC", b"GT"),
            two_bases(1887471, b"CA", b"TA"),
            two_bases(1887470, b"GC", b"GG"),
        ];
        for variant in &mut variants {
            variant
                .transcript_consequences
                .push(transcript_row("ENST00000421485", None));
            variant
                .transcript_consequences
                .push(transcript_row("ENST00000999999", None));
        }
        plugin.run_batch(&mut variants).unwrap();
        let comments: Vec<Vec<Option<&str>>> = variants
            .iter()
            .map(|v| {
                v.transcript_consequences
                    .iter()
                    .map(|tc| tc.plugin_data.get("scores_Comment").map(String::as_str))
                    .collect()
            })
            .collect();
        assert_eq!(
            comments,
            [
                vec![Some("Good"), Some("Other")],
                vec![Some("Good"), Some("Other")],
                vec![None, None],
            ]
        );
    }

    /// `trim=1` (the default) trims header names and cells, so ` Tolerated `
    /// is written as `Tolerated`; `trim=0` keeps the cell as the file has it
    /// (`:110`, `:246`).
    #[test]
    fn trim_0_keeps_whitespace() {
        let dir = tempfile::tempdir().unwrap();
        let run = |extra: &[&str]| {
            let plugin = init_with(dir.path(), extra).unwrap();
            let mut variants = vec![snv(1887471, "T")];
            variants[0]
                .transcript_consequences
                .push(transcript_row("ENST00000421485", None));
            plugin.run_batch(&mut variants).unwrap();
            variants[0].transcript_consequences[0]
                .plugin_data
                .get("scores_SIFT_pred")
                .cloned()
        };
        assert_eq!(run(&[]).as_deref(), Some("Tolerated"));
        assert_eq!(run(&["trim=1"]).as_deref(), Some("Tolerated"));
        assert_eq!(run(&["trim=0"]).as_deref(), Some(" Tolerated "));
    }

    /// `trim` is read as the truth test `if $self->{trim}` reads it (`:110`,
    /// `:246`), the value defaulting to 1 (`:164`): `2`, `yes` and `0.0` trim
    /// as `trim=1` does, and the empty value leaves the cell as `trim=0` does.
    #[test]
    fn trim_is_read_as_perl_truth() {
        let dir = tempfile::tempdir().unwrap();
        let run = |extra: &[&str]| {
            let plugin = init_with(dir.path(), extra).unwrap();
            let mut variants = vec![snv(1887471, "T")];
            variants[0]
                .transcript_consequences
                .push(transcript_row("ENST00000421485", None));
            plugin.run_batch(&mut variants).unwrap();
            variants[0].transcript_consequences[0]
                .plugin_data
                .get("scores_SIFT_pred")
                .cloned()
        };
        for on in ["trim=2", "trim=yes", "trim=0.0", "trim=01"] {
            assert_eq!(run(&[on]).as_deref(), Some("Tolerated"), "{on}");
        }
        assert_eq!(run(&["trim="]).as_deref(), Some(" Tolerated "));
    }

    /// A header line ending in a tab names no column past `Comment`: Ensembl
    /// VEP's `split /\t/` drops the trailing empty name (`:109`).
    #[test]
    fn a_trailing_empty_header_name_is_not_a_column() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_table(&empty_cells_fixture(dir.path()), &[]).unwrap();
        assert_eq!(
            keys(&plugin),
            ["cells_SIFT_score", "cells_SIFT_pred", "cells_Comment"]
        );
    }

    /// Ensembl VEP 116.2 splits a row on tabs, dropping its trailing empty
    /// fields, then trims each field under `trim=1` (`:245-246`), and assigns
    /// the fields to the column names (`:250`): a cell the row has and leaves
    /// empty, or trims to nothing, is the defined empty string, written as an
    /// empty tab cell and `""` in JSON, where a column past the row's fields is
    /// undefined, `-` and no key. The row names its empty-string fields in
    /// `plugin_empty_strings`; under `trim=0` a whitespace-only cell is its
    /// own text and no field is an empty string.
    #[test]
    fn an_empty_cell_the_row_has_is_a_string_and_a_missing_one_undefined() {
        let dir = tempfile::tempdir().unwrap();
        let file = empty_cells_fixture(dir.path());
        let run = |extra: &[&str]| {
            let plugin = init_table(&file, extra).unwrap();
            let mut variants: Vec<InputVariant> =
                (3000001..=3000004).map(|pos| snv(pos, "T")).collect();
            for variant in &mut variants {
                variant
                    .transcript_consequences
                    .push(transcript_row("ENST00000421485", None));
            }
            plugin.run_batch(&mut variants).unwrap();
            variants
                .into_iter()
                .map(|v| {
                    let tc = &v.transcript_consequences[0];
                    let values: Vec<String> = tc.plugin_data.values().cloned().collect();
                    (values, tc.plugin_empty_strings.clone())
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(
            run(&[]),
            [
                (
                    strings(&["", "", "EmptyMiddle"]),
                    strings(&["cells_SIFT_score", "cells_SIFT_pred"])
                ),
                (strings(&["0.44", "Tolerated", ""]), Vec::new()),
                (
                    strings(&["", "Tolerated", "SpaceMiddle"]),
                    strings(&["cells_SIFT_score"])
                ),
                (
                    strings(&["0.86", "Tolerated", ""]),
                    strings(&["cells_Comment"])
                ),
            ]
        );
        assert_eq!(
            run(&["trim=0"]),
            [
                (
                    strings(&["", "", "EmptyMiddle"]),
                    strings(&["cells_SIFT_score", "cells_SIFT_pred"])
                ),
                (strings(&["0.44", "Tolerated", ""]), Vec::new()),
                (strings(&[" ", "Tolerated", "SpaceMiddle"]), Vec::new()),
                (strings(&["0.86", "Tolerated", "  "]), Vec::new()),
            ]
        );
        // `cols=` keeps the rule to the selected columns.
        let selected = run(&["cols=Comment"]);
        assert_eq!(selected[0], (strings(&["EmptyMiddle"]), Vec::new()));
        assert_eq!(selected[3], (strings(&[""]), strings(&["cells_Comment"])));
    }

    /// Ensembl VEP's `_trim_whitespaces` runs `\s` on the file's bytes (`:91`),
    /// which the ASCII whitespace alone matches, so a cell holding one no-break
    /// space (U+00A0) is written as that character under `trim=1`, not emptied.
    #[test]
    fn a_no_break_space_is_not_whitespace_to_trim() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_table(&empty_cells_fixture(dir.path()), &[]).unwrap();
        let mut variants = vec![snv(3000005, "T")];
        variants[0]
            .transcript_consequences
            .push(transcript_row("ENST00000421485", None));
        plugin.run_batch(&mut variants).unwrap();
        let tc = &variants[0].transcript_consequences[0];
        assert_eq!(tc.plugin_data["cells_SIFT_pred"], "\u{a0}");
        assert!(tc.plugin_empty_strings.is_empty());
    }

    /// Two instances under one prefix write the same keys, and the later
    /// value replaces the earlier one (OutputFactory.pm:2230-2231): a row whose
    /// cell is the empty string in the first table and absent from the second
    /// ends undefined, and the reverse order ends as the empty string.
    #[test]
    fn a_later_instance_replaces_an_earlier_empty_string_mark() {
        let dir = tempfile::tempdir().unwrap();
        let empty_pred = write_tabix_fixture(
            dir.path(),
            "first.tsv.gz",
            &[HEADER, "21\t4000001\tC\tT\tENST00000421485\t0.1\t\tGood"],
            IndexSpec::standard(),
        );
        let short_row = write_tabix_fixture(
            dir.path(),
            "second.tsv.gz",
            &[HEADER, "21\t4000001\tC\tT\tENST00000421485\t0.2"],
            IndexSpec::standard(),
        );
        let run = |order: [&PathBuf; 2]| {
            let mut variants = vec![snv(4000001, "T")];
            variants[0]
                .transcript_consequences
                .push(transcript_row("ENST00000421485", None));
            for file in order {
                let plugin = init_table(file, &["prefix=X"]).unwrap();
                plugin.run_batch(&mut variants).unwrap();
            }
            let tc = &variants[0].transcript_consequences[0];
            (
                tc.plugin_data["X_SIFT_pred"].clone(),
                tc.plugin_empty_strings.clone(),
            )
        };
        assert_eq!(run([&empty_pred, &short_row]), (String::new(), Vec::new()));
        assert_eq!(
            run([&short_row, &empty_pred]),
            (String::new(), strings(&["X_SIFT_pred"]))
        );
    }

    /// A regulatory or motif row is not a transcript row, and a transcript the
    /// table does not name gets nothing.
    #[test]
    fn only_transcript_rows_of_a_named_transcript_are_annotated() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &[]).unwrap();
        let mut variants = vec![snv(1887471, "T")];
        variants[0].transcript_consequences.extend([
            transcript_row("ENST00000421485", None),
            TranscriptConsequence {
                feature_id: Some(Arc::from("ENST00000421485")),
                feature_type: FeatureType::RegulatoryFeature,
                ..TranscriptConsequence::default()
            },
            transcript_row("ENST00000000001", None),
        ]);
        plugin.run_batch(&mut variants).unwrap();
        let written: Vec<bool> = variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| !tc.plugin_data.is_empty())
            .collect();
        assert_eq!(written, [true, false, false]);
    }

    /// Ensembl VEP 116.2 never hands TranscriptAnnotator a structural variant:
    /// the module keeps the default `variant_feature_types` of
    /// `['VariationFeature']` (BaseVepPlugin.pm:100), which the dispatch checks
    /// before running a plugin on a line (OutputFactory.pm:2220). A `<DEL>`
    /// row of a named transcript gets no field, and nor does a matching
    /// substitution carrying the structural mark, which pins the gate rather
    /// than the allele match; the plain substitution is annotated.
    #[test]
    fn a_structural_variant_row_gets_no_field() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_with(dir.path(), &[]).unwrap();
        assert!(!plugin.annotates_structural_variants());
        let mut deletion = snv(1887471, "<DEL>");
        deletion.ref_allele = b"N".to_vec();
        deletion.end = 1887600;
        deletion.is_structural = true;
        let mut flagged = snv(1887471, "T");
        flagged.is_structural = true;
        let mut variants = vec![deletion, flagged, snv(1887471, "T")];
        for variant in &mut variants {
            variant
                .transcript_consequences
                .push(transcript_row("ENST00000421485", None));
        }
        plugin.run_batch(&mut variants).unwrap();
        let written: Vec<bool> = variants
            .iter()
            .map(|v| !v.transcript_consequences[0].plugin_data.is_empty())
            .collect();
        assert_eq!(written, [false, false, true]);
    }

    /// The substitutions of the short-rows table at `positions`, each on the
    /// transcript its valid row names.
    fn short_rows_variants(positions: &[u64]) -> Vec<InputVariant> {
        positions
            .iter()
            .map(|&pos| {
                let mut variant = snv(pos, "T");
                variant
                    .transcript_consequences
                    .push(transcript_row("ENST00000421485", None));
                variant
            })
            .collect()
    }

    /// The `Comment` field of the first transcript row of each variant.
    fn comments(variants: &[InputVariant]) -> Vec<Option<&str>> {
        variants
            .iter()
            .map(|v| {
                v.transcript_consequences[0]
                    .plugin_data
                    .get("short_Comment")
                    .map(String::as_str)
            })
            .collect()
    }

    /// The `Run` message `result` failed with; any other outcome fails the test.
    fn expect_run_error(result: Result<(), PluginError>) -> String {
        match result {
            Err(PluginError::Run(message)) => message,
            Ok(()) => panic!("expected PluginError::Run"),
            Err(other) => panic!("expected PluginError::Run, got {other:?}"),
        }
    }

    /// Ensembl VEP puts every row of the window to `get_matched_variant_alleles`
    /// before it reads the match and the transcript test together
    /// (TranscriptAnnotator.pm:222-237), and the match dies on a row whose
    /// reference or alternate cell `split /\t/` left undefined (`:245-247`;
    /// ensembl-variation Utils/Sequence.pm:1230, :968), whatever the row's
    /// transcript cell: a three-cell row before the matching row of its
    /// position fails the transcript row, which keeps none of the plugin's
    /// fields, while a variant at a position of valid rows alone keeps its own.
    #[test]
    fn a_row_with_no_allele_cell_before_the_matching_row_fails_the_transcript_row() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_table(&short_rows_fixture(dir.path()), &[]).unwrap();
        let mut variants = short_rows_variants(&[5000001, 5000011]);
        let message = expect_run_error(plugin.run_batch(&mut variants));
        assert!(
            message.contains("at 21:5000001 has no reference or alternate allele cell"),
            "{message}"
        );
        assert_eq!(comments(&variants), [None, Some("Clean")]);
        assert!(variants[0].transcript_consequences[0]
            .plugin_data
            .is_empty());
    }

    /// Ensembl VEP runs the module once per output line and catches each
    /// line's failure on its own (`run`, `:205-240`; OutputFactory.pm:2225-2245),
    /// so a transcript whose matching row precedes the short row of its
    /// position keeps its fields while a transcript whose scan reaches the
    /// short row loses them, whether the table never names it or names it in
    /// a row after the short one; the batch still fails with the short row.
    #[test]
    fn a_transcript_row_matched_before_the_short_row_keeps_its_fields_and_a_later_scan_loses_them()
    {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_table(&short_rows_fixture(dir.path()), &[]).unwrap();
        let mut variants = short_rows_variants(&[5000021, 5000061]);
        for variant in &mut variants {
            variant
                .transcript_consequences
                .push(transcript_row("ENST00000999999", None));
        }
        let message = expect_run_error(plugin.run_batch(&mut variants));
        assert!(
            message.contains("at 21:5000021 has no reference or alternate allele cell"),
            "{message}"
        );
        let comments: Vec<Vec<Option<&str>>> = variants
            .iter()
            .map(|v| {
                v.transcript_consequences
                    .iter()
                    .map(|tc| tc.plugin_data.get("short_Comment").map(String::as_str))
                    .collect()
            })
            .collect();
        assert_eq!(
            comments,
            [
                vec![Some("BeforeShort"), None],
                vec![Some("FirstTranscript"), None],
            ]
        );
        assert!(variants[1].transcript_consequences[1]
            .plugin_data
            .is_empty());
    }

    /// The loop returns at the first row that names the transcript and matches
    /// the alleles (`:237`), so a row of empty trailing cells after that row is
    /// never put to the match and the variant keeps the row's values.
    #[test]
    fn a_row_with_no_allele_cell_after_the_matching_row_does_not_fail_the_variant() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_table(&short_rows_fixture(dir.path()), &[]).unwrap();
        let mut variants = short_rows_variants(&[5000021]);
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(comments(&variants), [Some("BeforeShort")]);
    }

    /// A cell the row has is defined, so the match reads it instead of dying:
    /// an empty alternate cell before a filled cell and a whitespace-only
    /// trailing alternate cell, which `split /\t/` keeps and `trim=1` empties
    /// (`:245-246`), both compare as `-` (Utils/Sequence.pm:1033-1034) and
    /// match no substitution, and a four-cell row has its alleles but no
    /// transcript cell, which fails the transcript test alone; the valid row
    /// after each annotates the variant.
    #[test]
    fn a_defined_empty_alt_or_a_missing_transcript_cell_is_a_row_that_does_not_match() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = init_table(&short_rows_fixture(dir.path()), &[]).unwrap();
        let mut variants = short_rows_variants(&[5000031, 5000041, 5000051]);
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(
            comments(&variants),
            [
                Some("AfterEmptyAlt"),
                Some("AfterNoTranscript"),
                Some("AfterSpaceAlt")
            ]
        );
    }

    /// Ensembl VEP reads the rows from two bases before the variant's start to
    /// its end (TranscriptAnnotator.pm:218), two bases for an insertion, whose
    /// end precedes its start; the row written at the insertion's start, which
    /// vep-rs alone reads and matches, can therefore never die there: a
    /// three-cell row at that position fails nothing, while one inside the two
    /// bases Ensembl VEP reads fails the transcript row.
    #[test]
    fn a_row_with_no_allele_cell_fails_only_inside_the_window_ensembl_vep_reads() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            HEADER,
            "21\t5000098\tC",
            "21\t5000101\tC",
            "21\t5000101\t-\tA\tENST00000421485\t0.5\tTolerated\tAtStart",
            "21\t5000200\tC",
            "21\t5000201\t-\tA\tENST00000421485\t0.6\tTolerated\tInsideWindow",
        ];
        // Named as the other short-row table so the fields keep the `short_` prefix.
        let table = write_tabix_fixture(dir.path(), "short.tsv.gz", &rows, IndexSpec::standard());
        let plugin = init_table(&table, &[]).unwrap();
        let insertion = |start: u64| {
            let mut variant =
                InputVariant::new("21".into(), start, start - 1, b"-".to_vec(), b"A".to_vec());
            variant
                .transcript_consequences
                .push(transcript_row("ENST00000421485", None));
            variant
        };

        // The three-cell rows at 5000098 and 5000101 lie outside Ensembl VEP's
        // window 5000099..=5000100 for the insertion starting at 5000101: no
        // failure, and the row at the start is matched as documented.
        let mut variants = vec![insertion(5000101)];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(comments(&variants), [Some("AtStart")]);

        // The three-cell row at 5000200 lies inside 5000199..=5000200 for the
        // insertion starting at 5000201 and is reached before its matching row.
        let mut variants = vec![insertion(5000201)];
        let message = expect_run_error(plugin.run_batch(&mut variants));
        assert!(message.contains("21:5000200"), "{message}");
        assert_eq!(comments(&variants), [None]);
    }

    /// A row's Pos cell is keyed by `get_matched_variant_alleles`
    /// (`TranscriptAnnotator.pm:222-237`) as text when the alleles lose no
    /// leading base (ensembl-variation Utils/Sequence.pm:1232), and as its
    /// number plus the bases moved once a shared leading base leaves them
    /// (`:971, :1013`): for a substitution `1887471` annotates and `1887471.0`
    /// and `1887471.5` do not, nor does `1887472x` at 1887472; for the deletion
    /// `CT>C` at `1887480.0` and `1887500x` the moved start is the variant's,
    /// while `1887490.5` moves to 1887491.5 and matches nothing.
    #[test]
    fn a_position_cell_is_matched_as_the_key_holds_it() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            HEADER,
            "21\t1887471\tC\tT\tENST00000421485\t0.86\tTolerated\tText",
            "21\t1887471.0\tC\tA\tENST00000421485\t0.11\tTolerated\tPoint zero",
            "21\t1887471.5\tC\tG\tENST00000421485\t0.12\tTolerated\tPoint five",
            "21\t1887472x\tC\tT\tENST00000421485\t0.13\tTolerated\tTrailing",
            "21\t1887480.0\tCT\tC\tENST00000421485\t0.21\tTolerated\tMoved",
            "21\t1887490.5\tCT\tC\tENST00000421485\t0.22\tTolerated\tMoved half",
            "21\t1887500x\tCT\tC\tENST00000421485\t0.23\tTolerated\tMoved trailing",
        ];
        let file = write_tabix_fixture(dir.path(), "cells.tsv.gz", &rows, IndexSpec::standard());
        let plugin = init_table(&file, &[]).unwrap();
        let sift = |mut variant: InputVariant| {
            variant
                .transcript_consequences
                .push(transcript_row("ENST00000421485", None));
            let mut variants = vec![variant];
            plugin.run_batch(&mut variants).unwrap();
            variants[0].transcript_consequences[0]
                .plugin_data
                .get("cells_SIFT_score")
                .cloned()
        };
        assert_eq!(sift(snv(1887471, "T")), Some("0.86".to_string()));
        assert_eq!(sift(snv(1887471, "A")), None);
        assert_eq!(sift(snv(1887471, "G")), None);
        assert_eq!(sift(snv(1887472, "T")), None);
        let deletion =
            |start: u64| InputVariant::new("21".into(), start, start, b"T".to_vec(), b"-".to_vec());
        assert_eq!(sift(deletion(1887481)), Some("0.21".to_string()));
        assert_eq!(sift(deletion(1887491)), None);
        assert_eq!(sift(deletion(1887501)), Some("0.23".to_string()));

        // A row short of its allele cells inside the read window fails the
        // variant, unless its position cell is Perl-false, which returns the
        // match first (Sequence.pm:1173).
        let short = |pos: &str| TabixRecord {
            chr: "21".to_string(),
            start: 1887471,
            end: 1887471,
            columns: vec!["21".into(), pos.into(), "C".into()],
            header: None,
        };
        let matching = |pos: &str| {
            plugin.matching_values(&[&short(pos)], &snv(1887471, "T"), &["ENST00000421485"])
        };
        assert!(matches!(matching("1887471"), Err(PluginError::Run(_))));
        assert!(matches!(matching("0"), Ok(None)));
        assert!(matches!(matching(" 0 "), Ok(None)));
    }
}
