// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! PrimateAI plugin: deep-learning pathogenicity scores for single-base
//! substitutions, read from the tabix-indexed score file prepared as
//! `PrimateAI.pm:65-69` prepares it (the twelfth line, the column names,
//! prefixed with `#`; sorted; `tabix -s 1 -b 2 -e 2`).
//!
//! # Parameters
//!
//! The score file is the first positional parameter (`--plugin PrimateAI,<file>`,
//! `PrimateAI.pm:87`); `file=<path>` is accepted as well. The file must exist. A
//! path carrying `GRCh37` or `GRCh38` must agree with the run's assembly
//! (`:94-103`); a path carrying neither is logged with the run's assembly where
//! Ensembl VEP prints a line on standard output (`:104`).
//!
//! # Data file
//!
//! The header is the `#chr` line (`:114`); the `pos`, `alt` and `primateDL_score`
//! columns are located by name and each is required (`:123-138`). Records are
//! matched on the position alone (`:174-178`): the reference column is not read,
//! as `PrimateAI.pm` does not read it. Every row of the position is read before
//! the match, as the base module parses each line tabix returns before `run`
//! sees one (BaseVepTabixPlugin.pm:605-609): a row with no `alt` or score cell
//! leaves every allele at the position unscored (`:198`), and among intact rows
//! the first in file order carrying the allele supplies the score (`:174`).
//!
//! # Scoring
//!
//! Only a single-base alternate allele is scored (`:161`), on the forward strand
//! (`:169-171` reverse-complements the allele of a record given on the reverse
//! strand). A record with such an allele over more than one base (an insertion,
//! or several bases replaced by one) is queried over its span (`:164-166`), so a
//! short row inside the span unscores it as one at a position does, and matches
//! nothing, as the row's one position never equals both of its coordinates
//! (`:175-176`). The plugin declares `['Feature', 'Intergenic']` (`:144`), so the
//! score is a record-level value printed on every row of the allele, intergenic
//! rows included.
//!
//! # Output fields
//!
//! - `PrimateAI`: the `primateDL_score` of the matching row (`:151`, `:206`).
//!
//! `PrimateAI.pm` is the plugin's Ensembl VEP counterpart in Ensembl/VEP_plugins
//! release/116.

use std::path::Path;

use indexmap::IndexMap;
use tracing::{debug, info};

use crate::annotation_store::{check_table_or_store, open_best_store_with, AnnotationStore};
use crate::params::PluginParams;
use crate::tabix::{
    allele_match, header_lines, split_field_count, table_contig, BatchQueryResult, TabixAnnotator,
    TabixConfig, TabixRecord,
};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::assembly::Assembly;
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::variant::InputVariant;

const NAME: &str = "PrimateAI";

/// The output field.
const FIELD: &str = "PrimateAI";

/// PrimateAI score plugin.
pub struct PrimateAiPlugin {
    /// Annotation store (tabix or binary), opened by `init`.
    store: Option<Box<dyn AnnotationStore>>,
    /// 0-based index of the `pos` column.
    pos_col: usize,
    /// 0-based index of the `alt` column.
    alt_col: usize,
    /// 0-based index of the `primateDL_score` column.
    score_col: usize,
    /// The data file, for error messages.
    file: String,
    /// The run's assembly (the cache's own, else the value of `--assembly`),
    /// injected before `init()`.
    assembly: Option<Assembly>,
}

impl Default for PrimateAiPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl PrimateAiPlugin {
    /// Create an uninitialized PrimateAI plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            pos_col: 0,
            alt_col: 0,
            score_col: 0,
            file: String::new(),
            assembly: None,
        }
    }

    /// The assembly a path names: `GRCh37` is tested before `GRCh38`, as
    /// `PrimateAI.pm:95-99` tests them, on the whole path as given.
    fn path_assembly_tag(file: &str) -> Option<Assembly> {
        if file.contains("GRCh37") {
            Some(Assembly::Grch37)
        } else if file.contains("GRCh38") {
            Some(Assembly::Grch38)
        } else {
            None
        }
    }

    /// Refuse a score file whose path names the other assembly (`PrimateAI.pm:94-103`).
    /// A path naming neither is accepted and logged; a path naming an assembly is
    /// accepted without a log line when the run gave no assembly to compare with.
    fn check_assembly_tag(&self, file: &str) -> Result<(), PluginError> {
        match (Self::path_assembly_tag(file), self.assembly) {
            (Some(tag), Some(run)) if tag != run => Err(PluginError::Init(format!(
                "{NAME}: the path of the scores file {file} names {tag}, but the run's \
                 assembly is {run}; check the file and the --assembly option"
            ))),
            (Some(_), _) => Ok(()),
            (None, run) => {
                info!(
                    file,
                    assembly = run.map(|a| a.as_str()).unwrap_or("unset"),
                    "the PrimateAI scores file names no assembly; it is used as given"
                );
                Ok(())
            }
        }
    }

    /// The forward-strand alternate allele of a record the module queries the
    /// file for: one base of `ACGT` (`PrimateAI.pm:161`), whatever the record's
    /// span, which `:174-178` compare with the row's one position afterwards.
    fn scorable_alt(variant: &InputVariant) -> Option<u8> {
        let (_, alt) = variant.forward_strand_alleles();
        match alt.as_ref() {
            [base @ (b'A' | b'C' | b'G' | b'T')] => Some(*base),
            _ => None,
        }
    }

    /// The span the module queries for a record (`PrimateAI.pm:164-166`): its
    /// start and end in order, so an insertion is queried from its anchor base,
    /// and a substitution of several bases for one over every base it replaces.
    fn span(variant: &InputVariant) -> (u64, u64) {
        if variant.start <= variant.end {
            (variant.start, variant.end)
        } else {
            (variant.end, variant.start)
        }
    }

    /// The score of the first row of the span whose position cell, as `:175-176`'s
    /// `==` numifies it ([`allele_match::perl_position`]), is both `start` and
    /// `end` and whose `alt` is `alt`, in file order (`my ($res) = grep`,
    /// `PrimateAI.pm:174-178`), once every row of the span has been read: the
    /// base module parses each line tabix returns before `run` sees one
    /// (BaseVepTabixPlugin.pm:605-609), so a row with no `alt` or score cell
    /// dies for the record (`:198`) whatever the other rows hold, the ones
    /// before it included, and a span of more than one base, which no row's one
    /// position equals, is read and dies the same way while matching nothing.
    fn first_score<'a>(
        &self,
        records: impl IntoIterator<Item = &'a TabixRecord>,
        (start, end): (u64, u64),
        alt: u8,
    ) -> Result<Option<&'a String>, PluginError> {
        let mut score = None;
        for record in records {
            // A cell past the row's last filled one is the undefined value
            // `split /\t/` leaves (`:188`, `:198`); a whitespace-only one is a
            // defined score.
            let present = split_field_count(&record.columns);
            let cell = |column: usize| record.columns.get(column).filter(|_| column < present);
            let (Some(file_alt), Some(file_score)) = (cell(self.alt_col), cell(self.score_col))
            else {
                return Err(PluginError::Run(format!(
                    "{NAME}: one or more necessary values (pos, alt or primateDL_score) \
                     are missing from a row of {}",
                    self.file
                )));
            };
            if score.is_none()
                && file_alt.as_bytes() == [alt]
                && allele_match::perl_position(record, self.pos_col) == Some(start)
                && start == end
            {
                score = Some(file_score);
            }
        }
        Ok(score)
    }
}

impl crate::traits::sealed::Sealed for PrimateAiPlugin {}

impl BuiltinPlugin for PrimateAiPlugin {
    fn name(&self) -> &str {
        NAME
    }

    /// The Ensembl VEP module declares `['Feature', 'Intergenic']`
    /// (`PrimateAI.pm:143-144`): every row, the intergenic row included.
    fn feature_types(&self) -> &[RowKind] {
        RowKind::ALL
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            FIELD.to_string(),
            "PrimateAI score for variants".to_string(),
        )]
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn set_reference_fasta(&mut self, fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>) {
        if let Some(store) = self.store.as_mut() {
            store.set_reference_fasta(fasta);
        }
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(NAME, &["file"], params)?;
        let file = params
            .get("file")
            .or_else(|| params.positional().first().map(String::as_str))
            .filter(|file| !file.is_empty())
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: the scores file is the first parameter (--plugin {NAME},<file>)"
                ))
            })?;
        let path = Path::new(file);
        if !path.is_file() {
            return Err(PluginError::Init(format!(
                "{NAME}: scores file {file} not found"
            )));
        }
        // `add_file` runs before the header is read (PrimateAI.pm:92, 111).
        check_table_or_store(path)?;
        self.check_assembly_tag(file)?;

        // `PrimateAI.pm:113-130` keeps the last header line opening `#chr`,
        // whatever `#` lines sit before or after it, and `:118` strips its one
        // `#`; the column names it needs are found by name, wherever the file
        // puts them, and each is required, in the order `:136` checks them.
        let preamble = header_lines(path)?;
        let column_line = preamble
            .iter()
            .rev()
            .find(|line| line.starts_with("#chr"))
            .ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: could not read the #chr header line from {file}"
                ))
            })?;
        let header: Vec<&str> = column_line[1..].split('\t').collect();
        let column = |name: &str| {
            header.iter().position(|h| *h == name).ok_or_else(|| {
                PluginError::Init(format!(
                    "{NAME}: the required {name} header is missing from {file}"
                ))
            })
        };
        let pos_col = column("pos")?;
        let alt_col = column("alt")?;
        let score_col = column("primateDL_score")?;

        let config = TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: 0,
            start_col: pos_col,
            end_col: None,
            ref_col: None,
            alt_col: Some(alt_col),
            zero_based: false,
        };
        // A row keeps the whitespace before its line end, which Ensembl VEP's
        // `split /\t/` leaves on the score cell (`parse_data`, `:188`, `:196`),
        // so a whitespace-only score is a defined score, not a missing one.
        let mut store = open_best_store_with(
            path,
            config,
            self.assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )?;
        if let Some(frame) = store.mitochondrial_frame_mut() {
            frame.set_allele_column(header.iter().position(|name| *name == "ref"));
        }
        self.store = Some(store);
        self.pos_col = pos_col;
        self.alt_col = alt_col;
        self.score_col = score_col;
        self.file = file.to_string();
        debug!(
            file,
            pos_col, alt_col, score_col, "PrimateAI plugin initialized"
        );
        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|v| self.annotates_variant(v) && Self::scorable_alt(v).is_some())
            .map(|v| {
                let (start, end) = Self::span(v);
                (table_contig(v, self.assembly), start, end)
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
        let mut failure = None;
        for variant in variants.iter_mut().filter(|v| self.annotates_variant(v)) {
            let Some(alt) = Self::scorable_alt(variant) else {
                continue;
            };
            let (start, end) = Self::span(variant);
            match self.first_score(
                batch.lookup(table_contig(variant, self.assembly), start, end),
                (start, end),
                alt,
            ) {
                // An empty score `first_score` returns is a cell with a filled
                // cell after it, Ensembl VEP's empty string (`PrimateAI.pm:196`,
                // `:204-207`), which the record names.
                Ok(Some(score)) => {
                    variant.set_plugin_value(
                        FIELD.to_string(),
                        score.clone(),
                        true,
                        self.feature_types(),
                    );
                }
                Ok(None) => {}
                Err(error) => {
                    failure.get_or_insert(error);
                }
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
    use vep_core::coordinate::Strand;

    /// The heading line of this module's made-up score files.
    const HEADER: &str = "#chr\tpos\tref\talt\trefAA\taltAA\tstrand_1pos_0neg\
        \ttrinucleotide_context\tUCSC_gene\tExAC_coverage\tprimateDL_score";

    /// A made-up score file in the shape the plugin's recipe leaves
    /// (PrimateAI.pm:65, the twelfth line prefixed with `#`): eleven comment
    /// lines, the `#chr` header on line twelve, then `rows`.
    fn fixture(dir: &Path, name: &str, header: &str, rows: &[impl AsRef<str>]) -> String {
        let preamble: Vec<String> = (1..=11).map(|i| format!("# preamble line {i}")).collect();
        let mut lines: Vec<&str> = preamble.iter().map(String::as_str).collect();
        lines.push(header);
        lines.extend(rows.iter().map(|row| row.as_ref()));
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
            .display()
            .to_string()
    }

    /// A made-up row of `HEADER`.
    fn row(chr: &str, pos: u64, reference: &str, alt: &str, score: &str) -> String {
        format!("{chr}\t{pos}\t{reference}\t{alt}\tM\tL\t1\tAAT\tNM_000001.1\t50.00\t{score}")
    }

    fn loaded(path: &str) -> PrimateAiPlugin {
        let mut plugin = PrimateAiPlugin::new();
        plugin.init(&[path.to_string()]).expect("fixture loads");
        plugin
    }

    fn snv(pos: u64, reference: &[u8], alt: &[u8]) -> InputVariant {
        InputVariant::new("21".into(), pos, pos, reference.to_vec(), alt.to_vec())
    }

    fn score(plugin: &PrimateAiPlugin, variant: InputVariant) -> Option<String> {
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

    /// `add_file` checks for `<file>.tbi` before the header is read
    /// (PrimateAI.pm:92, 111): a table lacking both is refused naming the index.
    #[test]
    fn table_without_its_tbi_is_refused_before_its_header_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "scores_synthetic.tsv.bgz",
            &["#chr\tpos", "21\t100"],
            IndexSpec::standard(),
        );
        std::fs::remove_file(dir.path().join("scores_synthetic.tsv.bgz.tbi")).unwrap();
        let message = init_message(PrimateAiPlugin::new().init(&[path.display().to_string()]));
        assert!(
            message.contains("scores_synthetic.tsv.bgz.tbi not found"),
            "{message}"
        );
    }

    /// `PrimateAI.pm:87-88`: the file is `params->[0]` and must exist.
    #[test]
    fn file_is_the_first_positional_parameter_and_must_exist() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [row("chr21", 5030000, "A", "G", "0.8043")];
        let path = fixture(dir.path(), "scores_synthetic.tsv.bgz", HEADER, &rows);

        let message = init_message(PrimateAiPlugin::new().init(&[]));
        assert_eq!(
            message,
            "PrimateAI: the scores file is the first parameter (--plugin PrimateAI,<file>)"
        );
        let missing = format!("{}/no_such_synthetic.tsv.bgz", dir.path().display());
        let message = init_message(PrimateAiPlugin::new().init(std::slice::from_ref(&missing)));
        assert_eq!(
            message,
            format!("PrimateAI: scores file {missing} not found")
        );
        let message = init_message(
            PrimateAiPlugin::new().init(&[format!("file={path}"), "cols=x".to_string()]),
        );
        assert_eq!(
            message,
            "PrimateAI: unknown parameter 'cols'; PrimateAI accepts file"
        );

        let positional = loaded(&path);
        assert_eq!(
            score(&positional, snv(5030000, b"A", b"G")),
            Some("0.8043".to_string())
        );
        let mut keyed = PrimateAiPlugin::new();
        keyed.init(&[format!("file={path}")]).unwrap();
        assert_eq!(
            score(&keyed, snv(5030000, b"A", b"G")),
            Some("0.8043".to_string())
        );
    }

    /// `PrimateAI.pm:94-103`: a path naming an assembly must name the run's; a
    /// path naming none is used as given (`:104`), as is any path when the run
    /// gave no assembly to compare with.
    #[test]
    fn filename_assembly_tag_must_match_the_run() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [row("chr21", 5030000, "A", "G", "0.8043")];
        let g37 = fixture(
            dir.path(),
            "synthetic_score_fixture_GRCh37_sorted.tsv.bgz",
            HEADER,
            &rows,
        );
        let untagged = fixture(dir.path(), "synthetic_score_fixture.tsv.bgz", HEADER, &rows);

        let mut plugin = PrimateAiPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        let message = init_message(plugin.init(std::slice::from_ref(&g37)));
        assert_eq!(
            message,
            format!(
                "PrimateAI: the path of the scores file {g37} names GRCh37, but the run's \
                 assembly is GRCh38; check the file and the --assembly option"
            )
        );

        let mut plugin = PrimateAiPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch37));
        plugin
            .init(std::slice::from_ref(&g37))
            .expect("a GRCh37 file loads in a GRCh37 run");

        let mut plugin = PrimateAiPlugin::new();
        plugin.set_assembly(Some(Assembly::Grch38));
        plugin
            .init(&[untagged])
            .expect("a path naming no assembly loads");

        PrimateAiPlugin::new()
            .init(&[g37])
            .expect("without an assembly there is nothing to contradict");
    }

    /// `PrimateAI.pm:123-127`: `pos`, `alt` and `primateDL_score` are read by
    /// name, so a file that orders its columns differently scores the same.
    #[test]
    fn header_columns_are_located_by_name() {
        let dir = tempfile::tempdir().unwrap();
        let header = "#chr\tpos\tprimateDL_score\tUCSC_gene\talt\tref";
        let rows = [
            "chr21\t5030000\t0.1111\tNM_1.1\tC\tA",
            "chr21\t5030000\t0.2222\tNM_1.1\tG\tA",
            "chr21\t5030000\t0.3333\tNM_1.1\tT\tA",
        ];
        let path = fixture(dir.path(), "reordered_synthetic.tsv.bgz", header, &rows);
        let plugin = loaded(&path);
        assert_eq!(plugin.alt_col, 4);
        assert_eq!(plugin.score_col, 2);
        assert_eq!(
            score(&plugin, snv(5030000, b"A", b"G")),
            Some("0.2222".to_string())
        );
        assert_eq!(
            score(&plugin, snv(5030000, b"A", b"T")),
            Some("0.3333".to_string())
        );
        assert_eq!(score(&plugin, snv(5030001, b"A", b"G")), None);
    }

    /// `PrimateAI.pm:134-138`: no `#chr` line, or a required name missing from
    /// it, refuses the file.
    #[test]
    fn a_missing_required_header_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let rows = ["chr21\t5030000\tA\tG\t0.8043"];
        let no_score = fixture(
            dir.path(),
            "no_score_synthetic.tsv.bgz",
            "#chr\tpos\tref\talt\tscore",
            &rows,
        );
        let message = init_message(PrimateAiPlugin::new().init(std::slice::from_ref(&no_score)));
        assert_eq!(
            message,
            format!("PrimateAI: the required primateDL_score header is missing from {no_score}")
        );

        let no_pos = fixture(
            dir.path(),
            "no_pos_synthetic.tsv.bgz",
            "#chr\tposition\tref\talt\tprimateDL_score",
            &rows,
        );
        let message = init_message(PrimateAiPlugin::new().init(std::slice::from_ref(&no_pos)));
        assert_eq!(
            message,
            format!("PrimateAI: the required pos header is missing from {no_pos}")
        );

        let no_chr_line = fixture(
            dir.path(),
            "no_chr_line_synthetic.tsv.bgz",
            "#CHROM\tpos\tref\talt\tprimateDL_score",
            &rows,
        );
        let message = init_message(PrimateAiPlugin::new().init(std::slice::from_ref(&no_chr_line)));
        assert_eq!(
            message,
            format!("PrimateAI: could not read the #chr header line from {no_chr_line}")
        );

        // `#chrom` opens with `#chr` and so passes `PrimateAI.pm:114`'s test.
        let chrom = fixture(
            dir.path(),
            "chrom_synthetic.tsv.bgz",
            "#chrom\tpos\tref\talt\tprimateDL_score",
            &rows,
        );
        assert_eq!(
            score(&loaded(&chrom), snv(5030000, b"A", b"G")),
            Some("0.8043".to_string())
        );
    }

    /// `PrimateAI.pm:114` keeps the header lines opening `#chr` alone, so a
    /// comment after the column line is skipped, and a `##chr` line, whose
    /// second character is not `c`, is no header.
    #[test]
    fn the_header_is_the_last_line_opening_hash_chr_wherever_the_comments_sit() {
        let dir = tempfile::tempdir().unwrap();
        let rows = ["chr21\t5030000\tA\tG\t0.8043"];
        let header = "#chr\tpos\tref\talt\tprimateDL_score";
        let lines = [
            "# preamble",
            header,
            "# a comment after the column line",
            rows[0],
        ];
        let comment_after = write_tabix_fixture(
            dir.path(),
            "comment_after_synthetic.tsv.bgz",
            &lines,
            IndexSpec::standard(),
        )
        .display()
        .to_string();
        assert_eq!(
            score(&loaded(&comment_after), snv(5030000, b"A", b"G")),
            Some("0.8043".to_string())
        );

        let lines = [
            "# preamble",
            "##chr\tpos\tref\talt\tprimateDL_score",
            rows[0],
        ];
        let double_hash = write_tabix_fixture(
            dir.path(),
            "double_hash_synthetic.tsv.bgz",
            &lines,
            IndexSpec::standard(),
        )
        .display()
        .to_string();
        let message = init_message(PrimateAiPlugin::new().init(std::slice::from_ref(&double_hash)));
        assert_eq!(
            message,
            format!("PrimateAI: could not read the #chr header line from {double_hash}")
        );
    }

    /// `PrimateAI.pm:161` scores an allele of one base alone, and `:174-178`
    /// match a row's one position against both coordinates, so a deletion, an
    /// insertion and a two-base substitution score nothing whatever the file holds.
    #[test]
    fn only_single_base_alleles_are_scored() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            row("chr21", 5030000, "A", "G", "0.8043"),
            row("chr21", 5030001, "T", "-", "0.9000"),
            row("chr21", 5030002, "-", "C", "0.9100"),
            row("chr21", 5030003, "AG", "CT", "0.9200"),
        ];
        let plugin = loaded(&fixture(
            dir.path(),
            "shapes_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));

        assert_eq!(
            score(&plugin, snv(5030000, b"A", b"G")),
            Some("0.8043".to_string())
        );
        let deletion =
            InputVariant::new("21".into(), 5030001, 5030001, b"T".to_vec(), b"-".to_vec());
        assert_eq!(score(&plugin, deletion), None);
        let insertion =
            InputVariant::new("21".into(), 5030003, 5030002, b"-".to_vec(), b"C".to_vec());
        assert_eq!(score(&plugin, insertion), None);
        let substitution = InputVariant::new(
            "21".into(),
            5030003,
            5030004,
            b"AG".to_vec(),
            b"CT".to_vec(),
        );
        assert_eq!(score(&plugin, substitution), None);
        assert_eq!(score(&plugin, snv(5030000, b"A", b"N")), None);
    }

    /// A record with a single-base allele over more than one base is queried
    /// over its span (`PrimateAI.pm:164-166`, the start and end in order) and
    /// matches no row, whose one position never equals both coordinates
    /// (`:175-176`): a row short of its score cell inside the span dies for the
    /// record (`:198`), as it does at a one-base position, where a span over
    /// intact rows alone is left unscored without a failure. A substitution of
    /// one base for one at the span's start is scored from its own row.
    #[test]
    fn a_single_base_allele_over_a_span_reads_the_span_and_matches_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            row("chr21", 5030000, "C", "A", "0.8043"),
            "chr21\t5030001\tT\tG\tL\tM\t1\tCTG\tNM_000001.1".to_string(),
            row("chr21", 5030005, "A", "C", "0.7000"),
            row("chr21", 5030006, "A", "C", "0.7100"),
        ];
        let plugin = loaded(&fixture(
            dir.path(),
            "span_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));

        // `CT/A` from 5030000 to 5030001: the short row at 5030001 is in the span.
        let two_for_one =
            InputVariant::new("21".into(), 5030000, 5030001, b"CT".to_vec(), b"A".to_vec());
        let mut variants = vec![two_for_one];
        let error = plugin.run_batch(&mut variants).unwrap_err().to_string();
        assert!(
            error
                .contains("one or more necessary values (pos, alt or primateDL_score) are missing"),
            "{error}"
        );
        assert_eq!(variants[0].plugin_data.get(FIELD), None);

        // An insertion of `C` between 5030005 and 5030006: both rows intact, no
        // row's position equals both coordinates, nothing scored, no failure.
        let insertion =
            InputVariant::new("21".into(), 5030006, 5030005, b"-".to_vec(), b"C".to_vec());
        assert_eq!(score(&plugin, insertion), None);
        assert_eq!(
            score(&plugin, snv(5030000, b"C", b"A")),
            Some("0.8043".to_string())
        );
    }

    /// The position cell is compared as `PrimateAI.pm:175-176`'s `==` numifies
    /// it: tabix indexes `5030010x` and `5030011.5` at their leading integers,
    /// and the module reads them as 5030010, which scores the variant there,
    /// and 5030011.5, which equals no position, so that row scores nothing.
    #[test]
    fn a_position_cell_is_matched_as_perl_numifies_it() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            "chr21\t5030010x\tA\tG\tM\tL\t1\tAAT\tNM_000001.1\t50.00\t0.6100".to_string(),
            "chr21\t5030011.5\tA\tG\tM\tL\t1\tAAT\tNM_000001.1\t50.00\t0.6200".to_string(),
        ];
        let plugin = loaded(&fixture(
            dir.path(),
            "positions_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));
        assert_eq!(
            score(&plugin, snv(5030010, b"A", b"G")),
            Some("0.6100".to_string())
        );
        assert_eq!(score(&plugin, snv(5030011, b"A", b"G")), None);
    }

    /// Ensembl VEP splits the line htslib hands it, the line end alone removed,
    /// on tabs (`parse_data`, PrimateAI.pm:188) and prints the score cell as
    /// read (`:196`, `:206`), so the whitespace before the line end stays on the
    /// score, and a whitespace-only score cell is a defined score rather than
    /// the missing value `:198` refuses the file for.
    #[test]
    fn a_trailing_space_on_the_score_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            row("chr21", 5030000, "A", "G", "0.8043 "),
            row("chr21", 5030001, "T", "C", " "),
        ];
        let plugin = loaded(&fixture(
            dir.path(),
            "trailing_space_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));
        assert_eq!(
            score(&plugin, snv(5030000, b"A", b"G")),
            Some("0.8043 ".to_string())
        );
        assert_eq!(
            score(&plugin, snv(5030001, b"T", b"C")),
            Some(" ".to_string())
        );
    }

    /// A row whose score cell is empty and last has no score cell to Ensembl
    /// VEP: `split /\t/` drops the trailing empty field and `:198` dies on the
    /// undefined value, which the plugin runner turns into the row's warning.
    #[test]
    fn an_empty_trailing_score_cell_is_the_missing_value() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [row("chr21", 5030000, "A", "G", "")];
        let plugin = loaded(&fixture(
            dir.path(),
            "empty_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));
        let mut variants = vec![snv(5030000, b"A", b"G")];
        match plugin.run_batch(&mut variants) {
            Err(PluginError::Run(message)) => assert!(
                message.contains("necessary values (pos, alt or primateDL_score) are missing"),
                "{message}"
            ),
            other => panic!("expected PluginError::Run, got {other:?}"),
        }
        assert!(variants[0].plugin_data.is_empty());
    }

    /// Ensembl VEP parses every row tabix returns for the position before the
    /// match is made (BaseVepTabixPlugin.pm:605-609, then PrimateAI.pm:174), so
    /// a row with no score cell AFTER the matching row dies for the whole
    /// position (`:198`): no allele at that position is scored, the ones whose
    /// rows precede it included, and a position with intact rows is unaffected.
    #[test]
    fn a_missing_score_cell_after_the_matching_row_unscores_the_position() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            row("chr21", 30282928, "G", "A", "0.4488"),
            row("chr21", 30282928, "G", "C", "0.6365"),
            row("chr21", 30282928, "G", "T", "0.1673"),
            row("chr21", 30282928, "G", "G", ""),
            row("chr21", 30396091, "C", "T", "0.9618"),
        ];
        let plugin = loaded(&fixture(
            dir.path(),
            "empty_after_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));
        let mut variants = vec![
            snv(30282928, b"G", b"T"),
            snv(30282928, b"G", b"A"),
            snv(30396091, b"C", b"T"),
        ];
        match plugin.run_batch(&mut variants) {
            Err(PluginError::Run(message)) => assert!(
                message.contains("necessary values (pos, alt or primateDL_score) are missing"),
                "{message}"
            ),
            other => panic!("expected PluginError::Run, got {other:?}"),
        }
        assert!(
            variants[0].plugin_data.is_empty(),
            "{:?}",
            variants[0].plugin_data
        );
        assert!(
            variants[1].plugin_data.is_empty(),
            "{:?}",
            variants[1].plugin_data
        );
        assert_eq!(
            variants[2].plugin_data.get(FIELD).map(String::as_str),
            Some("0.9618")
        );
    }

    /// The same die for a row short of the score column (`split /\t/` leaves
    /// the cell undefined, `:198`), after the matching row.
    #[test]
    fn a_short_row_after_the_matching_row_unscores_the_position() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            row("chr21", 30282928, "G", "A", "0.4488"),
            row("chr21", 30282928, "G", "T", "0.1673"),
            "chr21\t30282928\tG\tC".to_string(),
        ];
        let plugin = loaded(&fixture(
            dir.path(),
            "short_after_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));
        let mut variants = vec![snv(30282928, b"G", b"T")];
        assert!(matches!(
            plugin.run_batch(&mut variants),
            Err(PluginError::Run(_))
        ));
        assert!(
            variants[0].plugin_data.is_empty(),
            "{:?}",
            variants[0].plugin_data
        );
    }

    /// `my ($res) = grep {...}` (`:174-178`) takes the first row in file order
    /// whose `alt` is the allele: a second row for the same allele with another
    /// score changes nothing, and nor does a later row for another allele.
    #[test]
    fn the_first_row_for_the_allele_supplies_the_score() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            row("chr21", 30282928, "G", "A", "0.4488"),
            row("chr21", 30282928, "G", "T", "0.1673"),
            row("chr21", 30282928, "G", "T", "0.9999"),
            row("chr21", 30282928, "G", "C", "0.6365"),
        ];
        let plugin = loaded(&fixture(
            dir.path(),
            "two_rows_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));
        assert_eq!(
            score(&plugin, snv(30282928, b"G", b"T")),
            Some("0.1673".to_string())
        );
        assert_eq!(
            score(&plugin, snv(30282928, b"G", b"A")),
            Some("0.4488".to_string())
        );
    }

    /// A score cell left empty before a filled cell is Ensembl VEP's empty
    /// string: `split /\t/` keeps it (`PrimateAI.pm:188`), `:198` passes it as
    /// defined and `:204-207` return it, so the record names the field in its
    /// `plugin_empty_strings`. A scored row leaves no mark.
    #[test]
    fn an_empty_score_cell_before_a_filled_cell_is_the_empty_string() {
        let dir = tempfile::tempdir().unwrap();
        let header = format!("{HEADER}\tnote");
        let rows = [
            format!("{}\tmade_up", row("chr21", 5030000, "A", "G", "")),
            format!("{}\tmade_up", row("chr21", 5030010, "A", "G", "0.25")),
        ];
        let plugin = loaded(&fixture(
            dir.path(),
            "empty_before_filled_synthetic.tsv.bgz",
            &header,
            &rows,
        ));
        let mut variants = vec![snv(5030000, b"A", b"G"), snv(5030010, b"A", b"G")];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(variants[0].plugin_data[FIELD], "");
        assert_eq!(variants[0].plugin_empty_strings, [FIELD]);
        assert_eq!(variants[1].plugin_data[FIELD], "0.25");
        assert!(variants[1].plugin_empty_strings.is_empty());
    }

    /// An empty cell before a filled one is a defined empty string to
    /// `split /\t/`, so a row whose `alt` cell is empty dies nowhere (`:198`)
    /// and matches no allele: the other rows of the position score as usual.
    #[test]
    fn an_empty_alt_cell_before_the_score_is_a_row_that_matches_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            row("chr21", 30282928, "G", "A", "0.4488"),
            row("chr21", 30282928, "G", "T", "0.1673"),
            row("chr21", 30282928, "G", "", "0.5000"),
        ];
        let plugin = loaded(&fixture(
            dir.path(),
            "empty_alt_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));
        assert_eq!(
            score(&plugin, snv(30282928, b"G", b"T")),
            Some("0.1673".to_string())
        );
        assert_eq!(
            score(&plugin, snv(30282928, b"G", b"A")),
            Some("0.4488".to_string())
        );
    }

    /// `PrimateAI.pm:169-171`: the allele of a record given on the reverse strand
    /// is reverse-complemented before the file, whose alleles are on the forward
    /// strand, is read; `T/C` on the reverse strand takes the `G` row, not the `C` row.
    #[test]
    fn negative_strand_alleles_are_complemented() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            row("chr21", 43774726, "A", "C", "0.1000"),
            row("chr21", 43774726, "A", "G", "0.2500"),
        ];
        let plugin = loaded(&fixture(
            dir.path(),
            "strand_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));

        let mut reverse = snv(43774726, b"T", b"C");
        reverse.strand = Strand::Reverse;
        assert_eq!(score(&plugin, reverse), Some("0.2500".to_string()));
        assert_eq!(
            score(&plugin, snv(43774726, b"T", b"C")),
            Some("0.1000".to_string())
        );
    }

    /// `PrimateAI.pm:144` declares `['Feature', 'Intergenic']`: a record with no
    /// transcript is scored, and the score is the record's, printed on every row.
    #[test]
    fn intergenic_variants_are_scored() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [row("chr21", 19268555, "C", "T", "0.4321")];
        let plugin = loaded(&fixture(
            dir.path(),
            "intergenic_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));

        let intergenic = snv(19268555, b"C", b"T");
        assert!(intergenic.transcript_consequences.is_empty());
        let mut genic = snv(19268555, b"C", b"T");
        genic.transcript_consequences.extend([
            {
                let mut consequence = TranscriptConsequence::default();
                consequence.feature_type = FeatureType::Transcript;
                consequence
            },
            {
                let mut consequence = TranscriptConsequence::default();
                consequence.feature_type = FeatureType::RegulatoryFeature;
                consequence
            },
        ]);
        let mut variants = vec![intergenic, genic];
        plugin.run_batch(&mut variants).unwrap();
        for variant in &variants {
            assert_eq!(
                variant.plugin_data.get(FIELD).map(String::as_str),
                Some("0.4321")
            );
            assert!(variant
                .transcript_consequences
                .iter()
                .all(|tc| tc.plugin_data.is_empty()));
        }
    }

    /// Ensembl VEP 116.3 never hands PrimateAI a structural variant: the
    /// module keeps the default `variant_feature_types` of
    /// `['VariationFeature']` (BaseVepPlugin.pm:100), which the dispatch checks
    /// before running a plugin on a line (OutputFactory.pm:2229). A `<DEL>`
    /// gets no score, and nor does a scored substitution carrying the
    /// structural mark, which pins the gate rather than the allele test; the
    /// plain substitution is scored.
    #[test]
    fn a_structural_variant_gets_no_score() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [row("chr21", 5030000, "A", "G", "0.8043")];
        let plugin = loaded(&fixture(
            dir.path(),
            "scores_synthetic.tsv.bgz",
            HEADER,
            &rows,
        ));
        assert!(!plugin.annotates_structural_variants());
        let mut deletion = snv(5030000, b"N", b"<DEL>");
        deletion.end = 5030100;
        deletion.is_structural = true;
        assert_eq!(score(&plugin, deletion), None);
        let mut flagged = snv(5030000, b"A", b"G");
        flagged.is_structural = true;
        assert_eq!(score(&plugin, flagged), None);
        assert_eq!(
            score(&plugin, snv(5030000, b"A", b"G")),
            Some("0.8043".to_string())
        );
    }
}
