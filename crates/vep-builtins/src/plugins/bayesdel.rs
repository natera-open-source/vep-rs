// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! BayesDel plugin: a deleteriousness meta-score read from the tabix-indexed
//! score table prepared as `BayesDel.pm:54-62` prepares it (the per-chromosome
//! files concatenated under their `#` header lines, sorted, `tabix -s 1 -b 2 -e 2`).
//!
//! # Parameters
//!
//! The score file is `file=<path>` (`BayesDel.pm:96-99`) and must exist. Ensembl
//! VEP reads its header with `tabix -H` and dies when the file is not named or
//! the header is empty (`:103-111`), so the file is required here as well; the
//! header's names are not checked, as they are not checked there.
//!
//! # Data file
//!
//! Five positional columns, `chr pos ref alt score` (`:151`). A row matches a
//! record when its position is the record's start (the smaller of start and end,
//! `:129-131`) and its reference and alternate alleles both equal the record's
//! (`:140-144`), in the representation Ensembl VEP stores: the shared leading
//! base of a VCF indel removed and `-` for the empty allele, so an insertion is
//! matched at its anchor base and a deletion at its first deleted base. A
//! matched row lacking its fifth cell scores as the empty value and a present
//! empty fifth cell as the empty string, the two values `:151`'s five-variable
//! split hands the writers.
//!
//! The table's contig is read under the record's name with or without `chr`
//! (`:144` hands the name to `BaseVepTabixPlugin.pm:693-708`). On GRCh37 a
//! mitochondrial record reads a contig of its own spelling's pair alone, `MT`
//! and `chrMT` or `chrM` and `M`, which name two sequences there
//! ([`table_contig`]).
//!
//! # Scoring
//!
//! The record's alleles are compared on the forward strand, the strand of the
//! file's alleles: `BayesDel.pm:135-138` reverse-complements the alleles of a
//! record given on the reverse strand (through a `reverse_comp` the module does
//! not import, so that path dies there; a VCF record is never on that strand).
//! The plugin declares `['Feature', 'Intergenic']` (`:117`), so the score is a
//! record-level value printed on every row of the allele, intergenic rows included.
//!
//! # Output fields
//!
//! - `BayesDel`: the score of the matching row (`:121`, `:158`).
//!
//! `BayesDel.pm` is the plugin's Ensembl VEP counterpart in Ensembl/VEP_plugins
//! release/116.

use std::path::Path;

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{open_best_store_with, AnnotationStore};
use crate::params::PluginParams;
use crate::tabix::{
    allele_match, last_preamble_line, table_contig, BatchQueryResult, TabixAnnotator, TabixConfig,
};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::assembly::Assembly;
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::variant::InputVariant;

const NAME: &str = "BayesDel";

/// The output field.
const FIELD: &str = "BayesDel";

/// 0-based positional columns (`BayesDel.pm:151`).
const CHR_COL: usize = 0;
const POS_COL: usize = 1;
const REF_COL: usize = 2;
const ALT_COL: usize = 3;
const SCORE_COL: usize = 4;

/// BayesDel meta-score plugin.
pub struct BayesDelPlugin {
    /// Annotation store (tabix or binary), opened by `init`.
    store: Option<Box<dyn AnnotationStore>>,
    /// The run's assembly, which names a mitochondrial record's contig to the
    /// table ([`table_contig`]).
    assembly: Option<Assembly>,
}

impl Default for BayesDelPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl BayesDelPlugin {
    /// Create an uninitialized BayesDel plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            assembly: None,
        }
    }

    /// The coordinates `BayesDel.pm:129-131` queries: start and end in order,
    /// so an insertion (start = end + 1) is queried from its anchor base.
    fn span(variant: &InputVariant) -> (u64, u64) {
        if variant.start <= variant.end {
            (variant.start, variant.end)
        } else {
            (variant.end, variant.start)
        }
    }

    /// The score cell of the first row at the record's start carrying both of
    /// its forward-strand alleles (`BayesDel.pm:140-146`): `None` when no row
    /// matches, `Some(None)` for a matched row lacking its fifth cell, which
    /// `:151`'s five-variable split leaves undefined, and `Some(Some(cell))` for
    /// the cell as the line holds it, the empty string for a trailing empty
    /// cell, which that split keeps. The position cell is compared as `:141`'s
    /// `==` numifies it ([`allele_match::perl_position`]): a cell with a numeric
    /// prefix (`100x`, `100.0`) matches the start its number equals, which is
    /// where tabix indexed the row, and one that is no number reads 0.
    fn score_for(
        &self,
        batch: &BatchQueryResult,
        variant: &InputVariant,
    ) -> Option<Option<String>> {
        let (start, end) = Self::span(variant);
        let (reference, alt) = variant.forward_strand_alleles();
        batch
            .lookup(table_contig(variant, self.assembly), start, end)
            .into_iter()
            .filter(|record| allele_match::perl_position(record, POS_COL) == Some(start))
            .find(|record| {
                record
                    .columns
                    .get(REF_COL)
                    .is_some_and(|file_ref| file_ref.as_bytes() == reference.as_ref())
                    && record
                        .columns
                        .get(ALT_COL)
                        .is_some_and(|file_alt| file_alt.as_bytes() == alt.as_ref())
            })
            .map(|record| record.columns.get(SCORE_COL).cloned())
    }
}

impl BuiltinPlugin for BayesDelPlugin {
    fn name(&self) -> &str {
        NAME
    }

    /// The Ensembl VEP module declares `['Feature', 'Intergenic']`
    /// (`BayesDel.pm:116-117`): every row, the intergenic row included.
    fn feature_types(&self) -> &[RowKind] {
        &RowKind::ALL
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(FIELD.to_string(), "Deleteriousness meta-score".to_string())]
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(NAME, &["file"], params)?;
        let file = params
            .get("file")
            .filter(|file| !file.is_empty())
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: the scores file is given as file=<path> (--plugin {NAME},file=<path>)"
                ))
            })?;
        let path = Path::new(file);
        if !path.is_file() {
            return Err(PluginError::Init(format!("{NAME}: file {file} not found")));
        }
        let config = TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: CHR_COL,
            start_col: POS_COL,
            end_col: None,
            ref_col: Some(REF_COL),
            alt_col: Some(ALT_COL),
            zero_based: false,
        };
        // The score cell keeps its trailing whitespace, which Ensembl VEP's
        // `split /\t/` leaves in the score it prints (`BayesDel.pm:151, 158`).
        let store = open_best_store_with(
            path,
            config,
            self.assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )?;
        // `BayesDel.pm:103-111`: the header is the last line `tabix -H` prints,
        // with its first `#` removed and split on whitespace (`split` on a byte
        // string: ASCII whitespace), and must hold at least one token; the
        // token is never read, so a second `#` is one.
        let has_header = last_preamble_line(path)?
            .as_deref()
            .and_then(|line| line.strip_prefix('#'))
            .is_some_and(|names| names.split_ascii_whitespace().next().is_some());
        if !has_header {
            return Err(PluginError::Init(format!(
                "{NAME}: could not read headers from {file}"
            )));
        }
        self.store = Some(store);
        debug!(file, "BayesDel plugin initialized");
        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        if variants.is_empty() {
            return Ok(Some(PrefetchData::new(BatchQueryResult::empty())));
        }
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|v| self.annotates_variant(v))
            .map(|v| {
                let (start, end) = Self::span(v);
                (table_contig(v, self.assembly), start, end)
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
        for variant in variants.iter_mut().filter(|v| self.annotates_variant(v)) {
            let Some(cell) = self.score_for(&batch, variant) else {
                continue;
            };
            // A matched row lacking its score cell scores as Ensembl VEP's
            // undefined value, an unmarked empty value here; a present empty
            // cell is its empty string, which the record names (`BayesDel.pm:151`).
            let is_empty_string = cell.as_ref().is_some_and(String::is_empty);
            variant
                .plugin_data
                .insert(FIELD.to_string(), cell.unwrap_or_default());
            // A later plugin's value replaces an earlier one under the same key
            // (OutputFactory.pm:2230-2231), its mark too.
            variant.plugin_empty_strings.retain(|key| key != FIELD);
            if is_empty_string {
                variant.plugin_empty_strings.push(FIELD.to_string());
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
        // Record-level plugin: annotation happens in `annotate`.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use vep_core::consequence::FeatureType;

    /// The made-up score files' heading line; no name in it is checked (`BayesDel.pm:106-107`).
    const HEADER: &str = "#Chr\tPos\tRef\tAlt\tBayesDel_addAF";

    fn fixture(dir: &Path, name: &str, header: Option<&str>, rows: &[&str]) -> String {
        let mut lines = Vec::new();
        lines.extend(header);
        lines.extend(rows);
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
            .display()
            .to_string()
    }

    fn loaded(path: &str) -> BayesDelPlugin {
        let mut plugin = BayesDelPlugin::new();
        plugin
            .init(&[format!("file={path}")])
            .expect("fixture loads");
        plugin
    }

    fn variant(start: u64, end: u64, reference: &[u8], alt: &[u8]) -> InputVariant {
        InputVariant::new("21".into(), start, end, reference.to_vec(), alt.to_vec())
    }

    fn score(plugin: &BayesDelPlugin, variant: InputVariant) -> Option<String> {
        let mut variants = vec![variant];
        plugin.run_batch(&mut variants).unwrap();
        variants[0].plugin_data.get(FIELD).cloned()
    }

    fn init_message(result: Result<(), PluginError>) -> String {
        match result {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected PluginError::Init, got {other:?}"),
        }
    }

    /// `BayesDel.pm:96-111`: the file is `file=<path>`, it must exist, and a
    /// positional parameter names no file.
    #[test]
    fn file_is_mandatory_and_must_exist() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "a_synthetic_bayesdel_score_fixture.txt.gz",
            Some(HEADER),
            &["21\t5030000\tA\tG\t0.31337"],
        );
        let no_file =
            "BayesDel: the scores file is given as file=<path> (--plugin BayesDel,file=<path>)";
        assert_eq!(init_message(BayesDelPlugin::new().init(&[])), no_file);
        assert_eq!(
            init_message(BayesDelPlugin::new().init(std::slice::from_ref(&path))),
            no_file,
            "a positional parameter is not the file"
        );
        assert_eq!(
            init_message(BayesDelPlugin::new().init(&["file=".to_string()])),
            no_file
        );
        let missing = format!("{}/no_such.synthetic.txt.gz", dir.path().display());
        assert_eq!(
            init_message(BayesDelPlugin::new().init(&[format!("file={missing}")])),
            format!("BayesDel: file {missing} not found")
        );
        assert_eq!(
            init_message(BayesDelPlugin::new().init(&[format!("file={path}"), "cols=x".into()])),
            "BayesDel: unknown parameter 'cols'; BayesDel accepts file"
        );
        assert_eq!(
            score(&loaded(&path), variant(5030000, 5030000, b"A", b"G")),
            Some("0.31337".to_string())
        );
    }

    /// `BayesDel.pm:103-111`: the header is the last `#` line split on
    /// whitespace; a file whose `#` line holds no name, or that has no `#`
    /// line, is refused.
    #[test]
    fn an_empty_header_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let rows = ["21\t5030000\tA\tG\t0.31337"];
        let blank = fixture(dir.path(), "blank.synthetic.txt.gz", Some("#"), &rows);
        assert_eq!(
            init_message(BayesDelPlugin::new().init(&[format!("file={blank}")])),
            format!("BayesDel: could not read headers from {blank}")
        );
        let spaces = fixture(
            dir.path(),
            "spaces.synthetic.txt.gz",
            Some("#   \t  "),
            &rows,
        );
        assert_eq!(
            init_message(BayesDelPlugin::new().init(&[format!("file={spaces}")])),
            format!("BayesDel: could not read headers from {spaces}")
        );
        let none = fixture(dir.path(), "none.synthetic.txt.gz", None, &rows);
        assert_eq!(
            init_message(BayesDelPlugin::new().init(&[format!("file={none}")])),
            format!("BayesDel: could not read headers from {none}")
        );
        // The names themselves are not checked (`:106-107`).
        let other_names = fixture(dir.path(), "names.synthetic.txt.gz", Some("#a\tb"), &rows);
        assert_eq!(
            score(&loaded(&other_names), variant(5030000, 5030000, b"A", b"G")),
            Some("0.31337".to_string())
        );
    }

    /// Ensembl VEP 116.2 removes one `#` from the last header line and tests
    /// that a token is left (`BayesDel.pm:106`, `:111`), so a line of hash
    /// marks alone is a header: its second `#` is the token.
    #[test]
    fn a_line_of_hash_marks_alone_is_a_header() {
        let dir = tempfile::tempdir().unwrap();
        let rows = ["21\t5030000\tA\tG\t0.31337"];
        for (name, header) in [
            ("two.synthetic.txt.gz", "##"),
            ("three.synthetic.txt.gz", "###"),
            ("padded.synthetic.txt.gz", "##  "),
        ] {
            let path = fixture(dir.path(), name, Some(header), &rows);
            assert_eq!(
                score(&loaded(&path), variant(5030000, 5030000, b"A", b"G")),
                Some("0.31337".to_string()),
                "header {header:?} loads"
            );
        }
        // The last header line decides (`:104-108`): a preamble whose last line
        // is a lone `#` is refused behind any number of named lines.
        let last_blank = fixture(
            dir.path(),
            "last_blank.synthetic.txt.gz",
            Some(HEADER),
            &["#", rows[0]],
        );
        assert_eq!(
            init_message(BayesDelPlugin::new().init(&[format!("file={last_blank}")])),
            format!("BayesDel: could not read headers from {last_blank}")
        );
    }

    /// `BayesDel.pm:151` splits a row into five columns by position; the header's
    /// names play no part in it.
    #[test]
    fn columns_are_positional_chr_pos_ref_alt_score() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "columns.synthetic.txt.gz",
            Some("#Score\tAlt\tRef\tPos\tChr"),
            &[
                "21\t5030000\tA\tC\t-0.11111",
                "21\t5030000\tA\tG\t0.22222",
                "21\t5030000\tA\tT\t0.33333",
            ],
        );
        let plugin = loaded(&path);
        assert_eq!(
            score(&plugin, variant(5030000, 5030000, b"A", b"C")),
            Some("-0.11111".to_string())
        );
        assert_eq!(
            score(&plugin, variant(5030000, 5030000, b"A", b"T")),
            Some("0.33333".to_string())
        );
        assert_eq!(score(&plugin, variant(5030001, 5030001, b"A", b"G")), None);
    }

    /// `BayesDel.pm:140-144`: the row's reference and alternate alleles must
    /// both equal the record's; a record given on the reverse strand is compared
    /// through its forward-strand alleles.
    #[test]
    fn ref_and_alt_must_both_match() {
        use vep_core::coordinate::Strand;
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "alleles.synthetic.txt.gz",
            Some(HEADER),
            &["21\t5030000\tA\tG\t0.5", "21\t5030000\tC\tG\t0.6"],
        );
        let plugin = loaded(&path);
        assert_eq!(
            score(&plugin, variant(5030000, 5030000, b"A", b"G")),
            Some("0.5".to_string())
        );
        assert_eq!(
            score(&plugin, variant(5030000, 5030000, b"C", b"G")),
            Some("0.6".to_string())
        );
        assert_eq!(score(&plugin, variant(5030000, 5030000, b"T", b"G")), None);
        assert_eq!(score(&plugin, variant(5030000, 5030000, b"A", b"T")), None);
        let mut reverse = variant(5030000, 5030000, b"T", b"C");
        reverse.strand = Strand::Reverse;
        assert_eq!(score(&plugin, reverse), Some("0.5".to_string()));
    }

    /// `BayesDel.pm:117` declares `['Feature', 'Intergenic']`: a record with no
    /// transcript is scored, and the score is the record's, printed on every row.
    #[test]
    fn intergenic_variants_are_scored() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "intergenic.synthetic.txt.gz",
            Some(HEADER),
            &["21\t19268555\tC\tT\t-0.98765"],
        );
        let plugin = loaded(&path);
        let intergenic = variant(19268555, 19268555, b"C", b"T");
        assert!(intergenic.transcript_consequences.is_empty());
        let mut genic = variant(19268555, 19268555, b"C", b"T");
        genic.transcript_consequences.extend([
            TranscriptConsequence {
                feature_type: FeatureType::Transcript,
                ..TranscriptConsequence::default()
            },
            TranscriptConsequence {
                feature_type: FeatureType::MotifFeature,
                ..TranscriptConsequence::default()
            },
        ]);
        let mut variants = vec![intergenic, genic];
        plugin.run_batch(&mut variants).unwrap();
        for variant in &variants {
            assert_eq!(
                variant.plugin_data.get(FIELD).map(String::as_str),
                Some("-0.98765")
            );
            assert!(variant
                .transcript_consequences
                .iter()
                .all(|tc| tc.plugin_data.is_empty()));
        }
    }

    /// An indel is matched in the representation Ensembl VEP stores and
    /// `BayesDel.pm:140-144` compares: a deletion at its first deleted base with
    /// `-` as the alternate allele, an insertion at its anchor base (the smaller
    /// coordinate, `:129-131`) with `-` as the reference; a row in the VCF
    /// spelling, anchor base included, matches neither.
    #[test]
    fn indel_alleles_match_only_the_file_representation() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "indels.synthetic.txt.gz",
            Some(HEADER),
            &[
                "21\t5030000\tAT\tA\t0.1",
                "21\t5030001\tT\t-\t0.2",
                "21\t5030100\tA\tAG\t0.3",
                "21\t5030100\t-\tG\t0.4",
                "21\t5030200\tAG\tCT\t0.5",
            ],
        );
        let plugin = loaded(&path);
        // VCF `21 5030000 AT A`: start 5030001, end 5030001, T/-.
        assert_eq!(
            score(&plugin, variant(5030001, 5030001, b"T", b"-")),
            Some("0.2".to_string())
        );
        // VCF `21 5030100 A AG`: start 5030101, end 5030100, -/G.
        assert_eq!(
            score(&plugin, variant(5030101, 5030100, b"-", b"G")),
            Some("0.4".to_string())
        );
        // A two-base substitution matches its whole alleles.
        assert_eq!(
            score(&plugin, variant(5030200, 5030201, b"AG", b"CT")),
            Some("0.5".to_string())
        );
        // A file holding the VCF spellings alone, anchor base included, scores
        // neither record: the record never carries the anchor base.
        let vcf_spelled = fixture(
            dir.path(),
            "vcf_spelled.synthetic.txt.gz",
            Some(HEADER),
            &["21\t5030000\tAT\tA\t0.1", "21\t5030100\tA\tAG\t0.3"],
        );
        let plugin = loaded(&vcf_spelled);
        assert_eq!(score(&plugin, variant(5030001, 5030001, b"T", b"-")), None);
        assert_eq!(score(&plugin, variant(5030101, 5030100, b"-", b"G")), None);
    }

    /// `BayesDel.pm:151` splits the line htslib hands over, the line end alone
    /// removed, and `:158` prints the score cell as split: a trailing space
    /// stays in the score and a whitespace-only cell is a score of one space.
    #[test]
    fn a_trailing_space_on_the_score_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "trailing_whitespace.synthetic.txt.gz",
            Some(HEADER),
            &["21\t30396091\tC\tT\t-1.18855 ", "21\t30396092\tC\tA\t "],
        );
        let plugin = loaded(&path);
        assert_eq!(
            score(&plugin, variant(30396091, 30396091, b"C", b"T")),
            Some("-1.18855 ".to_string())
        );
        assert_eq!(
            score(&plugin, variant(30396092, 30396092, b"C", b"A")),
            Some(" ".to_string())
        );
    }

    /// `BayesDel.pm:151` splits the line into five variables, so Perl keeps a
    /// trailing empty cell as the empty string and leaves a missing cell
    /// undefined, and `:146` returns the matched row's result either way. Both
    /// print as `BayesDel=` in the default format; the tab format prints an
    /// empty cell against `-` and JSON `""` against no key (`Tab.pm:128`,
    /// `JSON.pm:327`), so the record names the field of a present empty cell in
    /// `plugin_empty_strings` and not that of a missing one.
    #[test]
    fn a_row_lacking_its_score_cell_is_the_empty_value_and_an_empty_cell_the_empty_string() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "cells.synthetic.txt.gz",
            Some(HEADER),
            &[
                "21\t5030000\tA\tG",
                "21\t5030001\tA\tG\t",
                "21\t5030002\tA\tG\t\tx",
                "21\t5030003\tA\tG\t0.5",
            ],
        );
        let plugin = loaded(&path);
        let empty = Some(String::new());
        let marks = |pos: u64| {
            let mut variants = vec![variant(pos, pos, b"A", b"G")];
            plugin.run_batch(&mut variants).unwrap();
            (
                variants[0].plugin_data.get(FIELD).cloned(),
                variants[0].plugin_empty_strings.clone(),
            )
        };
        assert_eq!(marks(5030000), (empty.clone(), vec![]), "four cells");
        assert_eq!(
            marks(5030001),
            (empty.clone(), vec![FIELD.to_string()]),
            "five cells, the fifth empty"
        );
        assert_eq!(
            marks(5030002),
            (empty, vec![FIELD.to_string()]),
            "six cells, the fifth empty"
        );
        assert_eq!(marks(5030003), (Some("0.5".to_string()), vec![]));
    }

    /// The position cell is compared as `BayesDel.pm:141`'s `==` numifies it:
    /// tabix indexes `5030000x` and `5030001.0` at the integer they open with
    /// (`strtoll` from htslib 1.10), and Perl reads them as 5030000 and
    /// 5030001, so their rows match the variants there; `5030002.5` is indexed
    /// at 5030002 and reads as 5030002.5, so its row matches nothing.
    #[test]
    fn a_position_cell_is_matched_as_perl_numifies_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "positions.synthetic.txt.gz",
            Some(HEADER),
            &[
                "21\t5030000x\tA\tG\t0.1",
                "21\t5030001.0\tA\tG\t0.2",
                "21\t5030002.5\tA\tG\t0.3",
            ],
        );
        let plugin = loaded(&path);
        let score = |pos: u64| {
            let mut variants = vec![variant(pos, pos, b"A", b"G")];
            plugin.run_batch(&mut variants).unwrap();
            variants[0].plugin_data.get(FIELD).cloned()
        };
        assert_eq!(score(5030000), Some("0.1".to_string()));
        assert_eq!(score(5030001), Some("0.2".to_string()));
        assert_eq!(score(5030002), None);
    }

    /// Ensembl VEP 116.2 never hands BayesDel a structural variant: the module
    /// keeps the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2220). A `<DEL>` gets no score, and
    /// nor does a scored substitution carrying the structural mark, which pins
    /// the gate rather than the allele match; the plain substitution is scored.
    #[test]
    fn a_structural_variant_gets_no_score() {
        let dir = tempfile::tempdir().unwrap();
        let path = fixture(
            dir.path(),
            "scores.synthetic.txt.gz",
            Some(HEADER),
            &["21\t5030000\tA\tG\t0.31337"],
        );
        let plugin = loaded(&path);
        assert!(!plugin.annotates_structural_variants());
        let mut deletion = variant(5030000, 5030100, b"N", b"<DEL>");
        deletion.is_structural = true;
        assert_eq!(score(&plugin, deletion), None);
        let mut flagged = variant(5030000, 5030000, b"A", b"G");
        flagged.is_structural = true;
        assert_eq!(score(&plugin, flagged), None);
        assert_eq!(
            score(&plugin, variant(5030000, 5030000, b"A", b"G")),
            Some("0.31337".to_string())
        );
    }

    /// A one-row table naming the mitochondrial contig `contig`, opened on
    /// `assembly`, and a record at the row's site spelt `record`.
    fn mitochondrial_score(
        dir: &Path,
        assembly: Option<vep_core::assembly::Assembly>,
        record: &str,
        contig: &str,
    ) -> Option<String> {
        let path = fixture(
            dir,
            &format!("scores_{contig}.synthetic.txt.gz"),
            Some(HEADER),
            &[&format!("{contig}\t3460\tG\tA\t0.31337")],
        );
        let mut plugin = BayesDelPlugin::new();
        plugin.set_assembly(assembly);
        plugin
            .init(&[format!("file={path}")])
            .expect("fixture loads");
        let record = InputVariant::new(record.into(), 3460, 3460, b"G".to_vec(), b"A".to_vec());
        score(&plugin, record)
    }

    /// On GRCh37 `MT` and `chrMT` name the assembly's mitochondrial sequence
    /// and `chrM` and `M` another one, so a record reads a table's
    /// mitochondrial rows under its own spelling and its `chr` variant alone:
    /// the name `BayesDel.pm:144` hands `BaseVepTabixPlugin.pm:693-708`.
    #[test]
    fn on_grch37_a_mitochondrial_record_reads_the_table_within_its_pair_of_names() {
        let dir = tempfile::tempdir().unwrap();
        let grch37 = Some(vep_core::assembly::Assembly::Grch37);
        for (record, read) in [
            ("MT", ["MT", "chrMT"]),
            ("chrMT", ["chrMT", "MT"]),
            ("chrM", ["chrM", "M"]),
            ("M", ["M", "chrM"]),
        ] {
            for contig in ["M", "MT", "chrM", "chrMT"] {
                let expected = read.contains(&contig).then(|| "0.31337".to_string());
                assert_eq!(
                    mitochondrial_score(dir.path(), grch37, record, contig),
                    expected,
                    "record {record}, table {contig}"
                );
            }
        }
    }

    /// On GRCh38 the four names are one sequence, and a run without an assembly
    /// reads the same way: every spelling of the record reads every spelling of
    /// the table.
    #[test]
    fn off_grch37_a_mitochondrial_record_reads_every_spelling_of_the_table() {
        let dir = tempfile::tempdir().unwrap();
        for assembly in [Some(vep_core::assembly::Assembly::Grch38), None] {
            for record in ["MT", "chrMT", "chrM", "M"] {
                for contig in ["M", "MT", "chrM", "chrMT"] {
                    assert_eq!(
                        mitochondrial_score(dir.path(), assembly, record, contig),
                        Some("0.31337".to_string()),
                        "record {record}, table {contig}, {assembly:?}"
                    );
                }
            }
        }
    }
}
