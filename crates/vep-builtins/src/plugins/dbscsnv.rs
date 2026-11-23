// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! dbscSNV plugin: the AdaBoost and random-forest splice-site scores of a
//! single-base substitution from a tabix-indexed dbscSNV file, on its transcript
//! rows, as Ensembl VEP's `dbscSNV.pm` writes them.
//!
//! **Usage:** `--plugin dbscSNV,/path/to/dbscSNV1.1_GRCh38.txt.gz[,GRCh37|GRCh38]`
//!
//! | Parameter | Meaning |
//! | --- | --- |
//! | first positional | the bgzip-compressed, tabix-indexed dbscSNV file |
//! | second positional | the assembly the file was prepared for, in place of the run's (dbscSNV.pm:98-100, 214) |
//!
//! Output fields `ada_score` and `rf_score` (dbscSNV.pm:120-123, 137-140), read
//! by header name; a `.` value is not written (dbscSNV.pm:181). A record matches
//! on position and the alternate allele on the forward strand (dbscSNV.pm:152-168);
//! only a single-base substitution whose allele is `A`, `C`, `G` or `T` is
//! annotated (dbscSNV.pm:148, 155), and the first matching record is the
//! variant's (dbscSNV.pm:172-173), its scores printed on every transcript row.
//!
//! The plugin's documentation gives dbscSNV 1.1 coordinates for both assemblies
//! (dbscSNV.pm:46-47) and a header line that opens with `c` in place of `#`
//! (dbscSNV.pm:61-62). The reader takes a header without a leading `#` as the
//! column names and reads a line ending in CRLF with the CR dropped, as the
//! plugin drops it (dbscSNV.pm:191-193). The documentation prepares the file per
//! assembly: for GRCh37 indexed on columns 1 and 2 (`tabix -s 1 -b 2 -e 2 -c c`),
//! for GRCh38 sorted and indexed on columns 5 and 6
//! (`tabix -s 5 -b 6 -e 6 -c c`) (dbscSNV.pm:53-59). The position column is
//! `pos` under GRCh37 and `hg38_pos` under GRCh38 (dbscSNV.pm:209-231), resolved
//! by name and checked against the column the index was built on, since the two
//! preparations share the header and differ in their index. A file that keeps
//! only the columns the plugin reads (`chr pos ref alt ada_score rf_score`, or
//! `hg38_chr hg38_pos ...` for GRCh38) is read the same way. Ensembl VEP's
//! module refuses a header without `Ensembl_gene` (dbscSNV.pm:111-113), a
//! column no line of it reads (the test of dbscSNV.pm:169-170 is commented
//! out); the reader does not require it.
//!
//! Two-phase plugin: [`BuiltinPlugin::prefetch`] batches the tabix query for a
//! buffer and [`BuiltinPlugin::annotate`] writes the per-variant scores onto
//! the transcript rows.

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::{check_table_or_store, open_best_store_with, AnnotationStore};
use crate::params::{warn_unread_parameters, PluginParams};
use crate::tabix::{
    allele_match, indexed_position, parse_header, split_field_count, table_contig,
    BatchQueryResult, IndexedPosition, TabixAnnotator, TabixConfig, TabixRecord,
};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::assembly::Assembly;
use vep_core::consequence::TranscriptConsequence;
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

const NAME: &str = "dbscSNV";
/// The output fields and the header columns they are read from (dbscSNV.pm:120-123).
const SCORE_COLUMNS: [&str; 2] = ["ada_score", "rf_score"];

/// dbscSNV splice-site prediction plugin.
pub struct DbscSnvPlugin {
    /// Annotation store, initialized in `init()`.
    store: Option<Box<dyn AnnotationStore>>,
    /// Column of the assembly's position (dbscSNV.pm:215-220), the one the
    /// index is built on.
    pos_col: usize,
    /// Column of the alternate allele.
    alt_col: usize,
    /// Columns of `ada_score` and `rf_score`, where the header carries them.
    score_cols: [Option<usize>; 2],
    /// The run's assembly (the cache's own, else the value of `--assembly`),
    /// injected before `init()`.
    assembly: Option<Assembly>,
}

impl DbscSnvPlugin {
    /// Create a new uninitialized dbscSNV plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            pos_col: 0,
            alt_col: 0,
            score_cols: [None, None],
            assembly: None,
        }
    }

    /// The header name of the position column for `assembly` (dbscSNV.pm:215-220),
    /// `None` for an assembly the file carries no position column for.
    fn position_column_for(assembly: Assembly) -> Option<&'static str> {
        match assembly {
            Assembly::Grch37 => Some("pos"),
            Assembly::Grch38 => Some("hg38_pos"),
            _ => None,
        }
    }

    /// The header name of the sequence column for `assembly`, `None` for an
    /// assembly the file carries no sequence column for.
    fn sequence_column_for(assembly: Assembly) -> Option<&'static str> {
        match assembly {
            Assembly::Grch37 => Some("chr"),
            Assembly::Grch38 => Some("hg38_chr"),
            _ => None,
        }
    }

    /// The scores of the first record at the variant's position with its
    /// alternate allele (dbscSNV.pm:163-185); empty for a variant outside the
    /// plugin's shape or without a record. The position is the record's cell as
    /// `$tmp_data->{$pos_column} == $vf->{start}` (`:166`) numifies it
    /// ([`allele_match::perl_position`]): `20204701.0` and `20204701x` are
    /// 20204701, `20204701.5` is no position.
    fn scores(&self, variant: &InputVariant, records: &[&TabixRecord]) -> Scores {
        let mut scores = Scores::default();
        if variant.start != variant.end {
            return scores;
        }
        let (_, alt) = variant.forward_strand_alleles();
        if !matches!(alt.as_ref(), b"A" | b"C" | b"G" | b"T") {
            return scores;
        }
        let Some(record) = records.iter().find(|record| {
            allele_match::perl_position(record, self.pos_col) == Some(variant.start)
                && record.columns.get(self.alt_col).map(String::as_bytes) == Some(alt.as_ref())
        }) else {
            return scores;
        };
        // Every selected column the header names is written unless its value
        // is `.` (dbscSNV.pm:179-184): a cell past the row's last filled one
        // is the undefined value, which Ensembl VEP prints as `<name>=` in the
        // default format and `-` in the tab format; a cell the row has and
        // leaves empty is its empty string (:196 keys every cell `split /\t/`
        // reads), an empty tab cell and `""` in JSON, named on the row.
        let present = split_field_count(&record.columns);
        for (name, col) in SCORE_COLUMNS.iter().zip(self.score_cols) {
            let Some(col) = col else {
                continue;
            };
            let value = record.columns.get(col).map_or("", String::as_str);
            if value == "." {
                continue;
            }
            if value.is_empty() && col < present {
                scores.empty_strings.push(name.to_string());
            }
            scores.values.insert(name.to_string(), value.to_string());
        }
        scores
    }
}

/// A record's selected scores, and the keys among them whose value is Ensembl
/// VEP's empty string rather than its undefined value.
#[derive(Default)]
struct Scores {
    values: IndexMap<String, String>,
    empty_strings: Vec<String>,
}

impl Default for DbscSnvPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Debug for DbscSnvPlugin {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DbscSnvPlugin")
            .field("pos_col", &self.pos_col)
            .field("alt_col", &self.alt_col)
            .field("score_cols", &self.score_cols)
            .field("assembly", &self.assembly)
            .finish_non_exhaustive()
    }
}

impl crate::traits::sealed::Sealed for DbscSnvPlugin {}

impl BuiltinPlugin for DbscSnvPlugin {
    fn name(&self) -> &str {
        NAME
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            ("ada_score".to_string(), "dbscSNV ADA score".to_string()),
            ("rf_score".to_string(), "dbscSNV RF score".to_string()),
        ]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse(NAME, &[], params)?;
        let mut positional = params.positional().iter();
        let Some(file) = positional.next() else {
            return Err(PluginError::Init(format!(
                "{NAME}: no file given; the first parameter is the tabix-indexed dbscSNV file"
            )));
        };
        let path = std::path::PathBuf::from(file);
        // `if(my $assembly = $self->params->[1])` (dbscSNV.pm:98): an empty or
        // `0` second token is no assembly; a later token is never read there.
        let assembly = match positional.next().filter(|given| perl_true(given)) {
            Some(given) => match given.as_str() {
                "GRCh37" => Some(Assembly::Grch37),
                "GRCh38" => Some(Assembly::Grch38),
                other => {
                    return Err(PluginError::Init(format!(
                        "{NAME}: assembly \"{other}\" is not compatible with this plugin"
                    )))
                }
            },
            None => self.assembly,
        };
        // An assembly the file carries no position column for is refused as a
        // run without an assembly is.
        let position =
            assembly.and_then(|assembly| Some((assembly, Self::position_column_for(assembly)?)));
        let Some((assembly, position_name)) = position else {
            return Err(PluginError::Init(format!(
                "{NAME}: could not establish which position column to use based on assembly; \
                 try setting assembly manually with --assembly"
            )));
        };
        warn_unread_parameters(NAME, positional.as_slice(), "the file and an assembly");

        // `add_file` runs before the header is read (dbscSNV.pm:96, 103).
        check_table_or_store(&path)?;
        // The documented header has no `#` (dbscSNV.pm:61-62), so the header is the
        // first row whose second cell leads with no integer, `pos` in the fixtures'
        // header on both assemblies.
        let header = parse_header(&path, 1)?.ok_or_else(|| {
            PluginError::Init(format!("{NAME}: could not read headers from {file}"))
        })?;
        let column = |name: &str| header.iter().position(|h| h == name);
        // dbscSNV.pm:111-113 requires `alt` and `Ensembl_gene`. Of the two the
        // match reads `alt` alone (dbscSNV.pm:166-168; the gene test of :169-170
        // is commented out), so a header without `Ensembl_gene` is read.
        let Some(alt_col) = column("alt") else {
            return Err(PluginError::Init(format!(
                "{NAME}: could not find required column alt in {file}"
            )));
        };
        let Some(start_col) = column(position_name) else {
            return Err(PluginError::Init(format!(
                "{NAME}: could not find {position_name} column in {file}"
            )));
        };
        let mut chr_col = Self::sequence_column_for(assembly)
            .and_then(column)
            .unwrap_or(0);
        // The two preparations share the header and differ in their index, so
        // the index must name the assembly's position column: one built on the
        // other column would place every record at the other assembly's
        // position and match nothing.
        match indexed_position(&path, start_col)? {
            IndexedPosition::Unchecked => {}
            IndexedPosition::Agrees { sequence_column } => chr_col = sequence_column,
            IndexedPosition::Disagrees {
                indexed_column: indexed,
            } => {
                let indexed_name = header.get(indexed).map_or_else(
                    || format!("column {}", indexed + 1),
                    |name| format!("'{name}'"),
                );
                return Err(PluginError::Init(format!(
                    "{NAME}: a {assembly} run reads '{position_name}' (column {}), but {file} is indexed on \
                     {indexed_name} (column {}); prepare and index the file for {assembly} as the plugin's \
                     documentation describes, or point --plugin {NAME} at the copy indexed for it",
                    start_col + 1,
                    indexed + 1,
                )));
            }
        }

        let config = TabixConfig {
            file_path: path.clone(),
            chr_col,
            start_col,
            end_col: None,
            ref_col: column("ref"),
            alt_col: Some(alt_col),
            zero_based: false,
        };
        let store = open_best_store_with(
            &path,
            config,
            self.assembly,
            TabixAnnotator::keeping_trailing_whitespace,
        )?;
        self.pos_col = start_col;
        self.alt_col = alt_col;
        self.score_cols = [column(SCORE_COLUMNS[0]), column(SCORE_COLUMNS[1])];
        self.store = Some(store);
        debug!("{NAME} plugin initialized with {file} for {assembly}, position column '{position_name}'");
        Ok(())
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn set_reference_fasta(&mut self, fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>) {
        if let Some(store) = self.store.as_mut() {
            store.set_reference_fasta(fasta);
        }
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        // Only a single-base substitution is annotated, so only those are read.
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|v| self.annotates_variant(v) && v.start == v.end)
            .map(|v| (table_contig(v, self.assembly), v.start, v.end))
            .collect();
        if regions.is_empty() {
            return Ok(Some(PrefetchData::new(BatchQueryResult::empty())));
        }
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
            if !self.annotates_variant(variant) || variant.start != variant.end {
                continue;
            }
            let records = batch.lookup(
                table_contig(variant, self.assembly),
                variant.start,
                variant.end,
            );
            if records.is_empty() {
                continue;
            }
            let scores = self.scores(variant, &records);
            if scores.values.is_empty() {
                continue;
            }
            for tc in &mut variant.transcript_consequences {
                if self.annotates_row(tc) {
                    // A later plugin's value replaces an earlier one under the
                    // same key (OutputFactory.pm:2239-2240), its mark too.
                    tc.plugin_empty_strings
                        .retain(|key| !scores.values.contains_key(key));
                    tc.plugin_data
                        .extend(scores.values.iter().map(|(k, v)| (k.clone(), v.clone())));
                    tc.plugin_empty_strings
                        .extend(scores.empty_strings.iter().cloned());
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
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let store = self
            .store
            .as_ref()
            .ok_or_else(|| PluginError::Run(format!("{NAME} plugin not initialized")))?;
        if variant.start != variant.end {
            return Ok(IndexMap::new());
        }
        let records = store.query(
            table_contig(variant, self.assembly),
            variant.start,
            variant.end,
        )?;
        let records: Vec<&TabixRecord> = records.iter().collect();
        Ok(self.scores(variant, &records).values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use vep_core::consequence::FeatureType;

    /// The made-up fixture's header, without a `#` and with a CR before its line end; the
    /// fixture's rows carry both assemblies' coordinates.
    const HEADER: &str = "chr\tpos\tref\talt\thg38_chr\thg38_pos\tRefSeq?\tEnsembl?\tRefSeq_region\tRefSeq_gene\tRefSeq_functional_consequence\tRefSeq_id_c.change_p.change\tEnsembl_region\tEnsembl_gene\tEnsembl_functional_consequence\tEnsembl_id_c.change_p.change\tada_score\trf_score\r";

    /// A made-up row of `HEADER`.
    fn row(
        chr: &str,
        pos: u64,
        hg38_pos: u64,
        reference: &str,
        alt: &str,
        ada: &str,
        rf: &str,
    ) -> String {
        format!("{chr}\t{pos}\t{reference}\t{alt}\t{chr}\t{hg38_pos}\tY\tY\tintron\tGENE1\tintron\tNM_1:c.1A>G\tintron\tENSG1\tintron\tENST1:c.1A>G\t{ada}\t{rf}\r")
    }

    /// Rows at hg19 position 100 (hg38 300) and hg19 200 (hg38 400).
    fn rows() -> Vec<String> {
        vec![
            row("21", 100, 300, "G", "A", "0.00106", "0.028"),
            row("21", 100, 300, "G", "T", "0.91234", "."),
            row("21", 200, 400, "C", "T", ".", "0.874"),
        ]
    }

    /// The GRCh37 preparation the plugin documents (dbscSNV.pm:58-59): indexed on
    /// columns 1 and 2, the fixture's `chr` and `pos`.
    fn grch37_file(dir: &std::path::Path) -> String {
        let rows = rows();
        let mut lines: Vec<&str> = vec![HEADER];
        lines.extend(rows.iter().map(String::as_str));
        write_tabix_fixture(
            dir,
            "synthetic_scores_GRCh37.gz",
            &lines,
            IndexSpec::standard().bare_header(),
        )
        .to_string_lossy()
        .to_string()
    }

    /// The GRCh38 preparation the plugin documents (dbscSNV.pm:54-55): sorted and
    /// indexed on columns 5 and 6, the fixture's `hg38_chr` and `hg38_pos`.
    fn grch38_file(dir: &std::path::Path) -> String {
        let rows = rows();
        let mut lines: Vec<&str> = vec![HEADER];
        lines.extend(rows.iter().map(String::as_str));
        write_tabix_fixture(
            dir,
            "synthetic_scores_GRCh38.gz",
            &lines,
            IndexSpec::standard().with_span_cols(5, 6, 6).bare_header(),
        )
        .to_string_lossy()
        .to_string()
    }

    fn load(assembly: Option<Assembly>, params: &[&str]) -> Result<DbscSnvPlugin, PluginError> {
        let mut plugin = DbscSnvPlugin::new();
        plugin.set_assembly(assembly);
        let params: Vec<String> = params.iter().map(|p| p.to_string()).collect();
        plugin.init(&params)?;
        Ok(plugin)
    }

    fn variant(pos: u64, reference: &str, alt: &str) -> InputVariant {
        let mut variant = InputVariant::new(
            "21".into(),
            pos,
            pos,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        );
        variant
            .transcript_consequences
            .push(TranscriptConsequence::default());
        variant
    }

    fn scores_of(plugin: &DbscSnvPlugin, variant: &mut InputVariant) -> Vec<(String, String)> {
        plugin.run_batch(std::slice::from_mut(variant)).unwrap();
        variant.transcript_consequences[0]
            .plugin_data
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect()
    }

    /// `add_file` checks for `<file>.tbi` before the header is read
    /// (dbscSNV.pm:96, 103): a table lacking both is refused naming the index.
    #[test]
    fn table_without_its_tbi_is_refused_before_its_header_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_tabix_fixture(
            dir.path(),
            "synthetic_scores_GRCh38.gz",
            &["chr\tpos", "21\t100"],
            IndexSpec::standard().bare_header(),
        );
        std::fs::remove_file(dir.path().join("synthetic_scores_GRCh38.gz.tbi")).unwrap();
        let message = load(Some(Assembly::Grch38), &[&path.display().to_string()])
            .err()
            .map_or_else(|| "the table opened".to_string(), |e| e.to_string());
        assert!(
            message.contains("synthetic_scores_GRCh38.gz.tbi not found"),
            "{message}"
        );
    }

    /// The output fields are `ada_score` and `rf_score` with Ensembl VEP's
    /// descriptions (dbscSNV.pm:137-140).
    #[test]
    fn header_fields_are_ada_score_and_rf_score() {
        assert_eq!(
            DbscSnvPlugin::new().header_info(),
            vec![
                ("ada_score".to_string(), "dbscSNV ADA score".to_string()),
                ("rf_score".to_string(), "dbscSNV RF score".to_string()),
            ]
        );
    }

    /// A GRCh38 run reads `hg38_pos` (dbscSNV.pm:218-219), column 6 of the fixture,
    /// so a record is found at its hg38 position; a `.` score is not written
    /// (dbscSNV.pm:181); the fixture's CRLF line ends and bare header are read.
    #[test]
    fn grch38_file_resolves_hg38_pos_column_by_name() {
        let dir = tempfile::tempdir().unwrap();
        let file = grch38_file(dir.path());
        let plugin = load(Some(Assembly::Grch38), &[&file]).unwrap();
        assert_eq!(
            scores_of(&plugin, &mut variant(300, "G", "A")),
            [
                ("ada_score".to_string(), "0.00106".to_string()),
                ("rf_score".to_string(), "0.028".to_string())
            ]
        );
        assert_eq!(
            scores_of(&plugin, &mut variant(300, "G", "T")),
            [("ada_score".to_string(), "0.91234".to_string())]
        );
        assert_eq!(
            scores_of(&plugin, &mut variant(400, "C", "T")),
            [("rf_score".to_string(), "0.874".to_string())]
        );
        assert!(scores_of(&plugin, &mut variant(100, "G", "A")).is_empty());

        let file37 = grch37_file(dir.path());
        let plugin37 = load(Some(Assembly::Grch37), &[&file37]).unwrap();
        assert_eq!(
            scores_of(&plugin37, &mut variant(100, "G", "A")),
            [
                ("ada_score".to_string(), "0.00106".to_string()),
                ("rf_score".to_string(), "0.028".to_string())
            ]
        );
        assert!(scores_of(&plugin37, &mut variant(300, "G", "A")).is_empty());
    }

    /// A file without `alt`, or without the assembly's position column, is refused
    /// (dbscSNV.pm:111-118); so is the file indexed for the other assembly, naming both
    /// columns, and a run with no assembly at all (dbscSNV.pm:226).
    #[test]
    fn missing_required_column_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let no_alt = write_tabix_fixture(
            dir.path(),
            "no_alt_synthetic.txt.gz",
            &["chr\tpos\tref\tada_score\trf_score", "21\t100\tG\t0.5\t0.5"],
            IndexSpec::standard().bare_header(),
        );
        let err = load(Some(Assembly::Grch37), &[&no_alt.to_string_lossy()]).unwrap_err();
        assert!(
            err.to_string()
                .contains("could not find required column alt"),
            "{err}"
        );

        let six_columns_grch37 = write_tabix_fixture(
            dir.path(),
            "six37_synthetic.txt.gz",
            &[
                "chr\tpos\tref\talt\tada_score\trf_score\r",
                "21\t100\tG\tA\t0.00106\t0.028\r",
            ],
            IndexSpec::standard().bare_header(),
        );
        let err = load(
            Some(Assembly::Grch38),
            &[&six_columns_grch37.to_string_lossy()],
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("could not find hg38_pos column"),
            "{err}"
        );
        let plugin = load(
            Some(Assembly::Grch37),
            &[&six_columns_grch37.to_string_lossy()],
        )
        .unwrap();
        assert_eq!(scores_of(&plugin, &mut variant(100, "G", "A")).len(), 2);

        let file38 = grch38_file(dir.path());
        let err = load(Some(Assembly::Grch37), &[&file38]).unwrap_err();
        let message = err.to_string();
        assert!(
            message.contains("reads 'pos' (column 2)")
                && message.contains("indexed on 'hg38_pos' (column 6)"),
            "{message}"
        );
        let file37 = grch37_file(dir.path());
        let err = load(Some(Assembly::Grch38), &[&file37]).unwrap_err();
        assert!(
            err.to_string().contains("indexed on 'pos' (column 2)"),
            "{err}"
        );

        let err = load(None, &[&file37]).unwrap_err();
        assert!(
            err.to_string()
                .contains("could not establish which position column"),
            "{err}"
        );
    }

    /// dbscSNV.pm:111-113 dies on a header without `Ensembl_gene`, and the one
    /// test that reads the cell (dbscSNV.pm:169-170) is commented out, so a
    /// record matches on its position and alternate allele alone
    /// (dbscSNV.pm:166-168). A file of the six columns the module reads is
    /// loaded under either assembly's column names, and each of its rows scores
    /// as the same row of the eighteen-column fixture does.
    #[test]
    fn a_header_without_the_gene_column_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let whole = load(Some(Assembly::Grch37), &[&grch37_file(dir.path())]).unwrap();
        let six37 = write_tabix_fixture(
            dir.path(),
            "six_columns_synthetic_GRCh37.txt.gz",
            &[
                "chr\tpos\tref\talt\tada_score\trf_score\r",
                "21\t100\tG\tA\t0.00106\t0.028\r",
                "21\t100\tG\tT\t0.91234\t.\r",
                "21\t200\tC\tT\t.\t0.874\r",
            ],
            IndexSpec::standard().bare_header(),
        );
        let plugin = load(Some(Assembly::Grch37), &[&six37.to_string_lossy()]).unwrap();
        let expected = [
            (
                100,
                "G",
                "A",
                vec![("ada_score", "0.00106"), ("rf_score", "0.028")],
            ),
            (100, "G", "T", vec![("ada_score", "0.91234")]),
            (200, "C", "T", vec![("rf_score", "0.874")]),
        ];
        for (pos, reference, alt, scores) in expected {
            let scores: Vec<(String, String)> = scores
                .into_iter()
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect();
            assert_eq!(
                scores_of(&plugin, &mut variant(pos, reference, alt)),
                scores,
                "{pos} {alt}"
            );
            assert_eq!(
                scores_of(&whole, &mut variant(pos, reference, alt)),
                scores,
                "{pos} {alt} in the eighteen-column fixture"
            );
        }

        let six38 = write_tabix_fixture(
            dir.path(),
            "six_columns_synthetic_GRCh38.txt.gz",
            &[
                "hg38_chr\thg38_pos\tref\talt\tada_score\trf_score\r",
                "21\t300\tG\tA\t0.00106\t0.028\r",
            ],
            IndexSpec::standard().bare_header(),
        );
        let plugin = load(Some(Assembly::Grch38), &[&six38.to_string_lossy()]).unwrap();
        assert_eq!(
            scores_of(&plugin, &mut variant(300, "G", "A")),
            [
                ("ada_score".to_string(), "0.00106".to_string()),
                ("rf_score".to_string(), "0.028".to_string())
            ]
        );
    }

    /// The second positional parameter names the assembly the file was prepared for
    /// and outranks the run's (dbscSNV.pm:98-100, 214); a name other than `GRCh37` or
    /// `GRCh38` is refused (dbscSNV.pm:222), and a `key=value` parameter is not read
    /// by position.
    #[test]
    fn positional_assembly_override_wins() {
        let dir = tempfile::tempdir().unwrap();
        let file38 = grch38_file(dir.path());
        let plugin = load(Some(Assembly::Grch37), &[&file38, "GRCh38"]).unwrap();
        assert_eq!(scores_of(&plugin, &mut variant(300, "G", "A")).len(), 2);
        let plugin = load(None, &[&file38, "GRCh38"]).unwrap();
        assert_eq!(scores_of(&plugin, &mut variant(300, "G", "A")).len(), 2);
        let err = load(Some(Assembly::Grch38), &[&file38, "hg38"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "plugin initialization failed: dbscSNV: assembly \"hg38\" is not compatible with this plugin"
        );
        // A second token that is empty (as `<file>,,X` leaves it) or `0` is no
        // assembly: the Ensembl VEP `if(my $assembly = $self->params->[1])` reads
        // both as false and keeps the run's, and the third token is never read.
        let plugin = load(Some(Assembly::Grch38), &[&file38, "", "X"]).unwrap();
        assert_eq!(scores_of(&plugin, &mut variant(300, "G", "A")).len(), 2);
        let plugin = load(Some(Assembly::Grch38), &[&file38, "0"]).unwrap();
        assert_eq!(scores_of(&plugin, &mut variant(300, "G", "A")).len(), 2);
        let err = load(Some(Assembly::Grch38), &[&format!("file={file38}")]).unwrap_err();
        assert!(
            err.to_string().contains("unknown parameter 'file'"),
            "{err}"
        );
    }

    /// Every assembly `Assembly::ALL` lists has its position and sequence
    /// columns, so none is refused as an assembly the file carries no columns
    /// for.
    #[test]
    fn every_listed_assembly_has_its_columns() {
        for &assembly in Assembly::ALL {
            assert!(
                DbscSnvPlugin::position_column_for(assembly).is_some(),
                "{assembly}"
            );
            assert!(
                DbscSnvPlugin::sequence_column_for(assembly).is_some(),
                "{assembly}"
            );
        }
    }

    /// Only a single-base substitution whose allele is a base is annotated
    /// (dbscSNV.pm:148, 155), on the forward strand (dbscSNV.pm:153), and its scores
    /// land on every transcript row and on no regulatory or motif row (`feature_types`
    /// of `['Transcript']`).
    #[test]
    fn single_base_substitutions_on_transcript_rows_alone() {
        let dir = tempfile::tempdir().unwrap();
        let file = grch37_file(dir.path());
        let plugin = load(Some(Assembly::Grch37), &[&file]).unwrap();
        let mut deletion = InputVariant::new("21".into(), 100, 101, b"GA".to_vec(), b"G".to_vec());
        deletion
            .transcript_consequences
            .push(TranscriptConsequence::default());
        assert!(scores_of(&plugin, &mut deletion).is_empty());
        let mut reverse = variant(100, "C", "T");
        reverse.strand = vep_core::coordinate::Strand::Reverse;
        assert_eq!(scores_of(&plugin, &mut reverse).len(), 2);

        let mut variant = variant(100, "G", "A");
        variant.transcript_consequences.clear();
        for feature_type in [
            FeatureType::RegulatoryFeature,
            FeatureType::Transcript,
            FeatureType::MotifFeature,
            FeatureType::Transcript,
        ] {
            variant.transcript_consequences.push({
                let mut consequence = TranscriptConsequence::default();
                consequence.feature_type = feature_type;
                consequence
            });
        }
        plugin
            .run_batch(std::slice::from_mut(&mut variant))
            .unwrap();
        let scores: Vec<Option<&str>> = variant
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.get("ada_score").map(String::as_str))
            .collect();
        assert_eq!(scores, [None, Some("0.00106"), None, Some("0.00106")]);
        assert!(variant.plugin_data.is_empty());
    }

    /// `dbscSNV.pm:191-193` takes the line htslib hands over, the line end alone
    /// removed, drops a `\r` and splits it with `split /\t/`, so the `rf_score`
    /// cell keeps the whitespace before the line end, and a whitespace-only cell
    /// is not `.` (`:181`): `0.028 ` prints as `rf_score=0.028 ` and ` ` as
    /// `rf_score= `.
    #[test]
    fn a_trailing_space_on_the_last_cell_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            row("21", 100, 300, "G", "A", "0.00106", "0.028 "),
            row("21", 100, 300, "G", "T", "0.91234", " "),
        ];
        let mut lines: Vec<&str> = vec![HEADER];
        lines.extend(rows.iter().map(String::as_str));
        let file = write_tabix_fixture(
            dir.path(),
            "synthetic_scores_trailing_GRCh37.gz",
            &lines,
            IndexSpec::standard().bare_header(),
        );
        let plugin = load(Some(Assembly::Grch37), &[&file.to_string_lossy()]).unwrap();
        assert_eq!(
            scores_of(&plugin, &mut variant(100, "G", "A")),
            [
                ("ada_score".to_string(), "0.00106".to_string()),
                ("rf_score".to_string(), "0.028 ".to_string())
            ]
        );
        assert_eq!(
            scores_of(&plugin, &mut variant(100, "G", "T")),
            [
                ("ada_score".to_string(), "0.91234".to_string()),
                ("rf_score".to_string(), " ".to_string())
            ]
        );
    }

    /// `dbscSNV.pm:196` keys every header name, a cell the row lacks included,
    /// and `:181` drops a `.` alone, so an empty or absent `rf_score` is written
    /// as the empty value, `rf_score=` in the default format; only `.` is left
    /// out. A cell past the row's last filled one is Ensembl VEP's undefined
    /// value (`-` in the tab format), where a cell the row has and leaves empty
    /// is its empty string (an empty tab cell, `""` in JSON), which the row
    /// names in `plugin_empty_strings`.
    #[test]
    fn an_empty_or_absent_score_cell_is_the_empty_value() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            row("21", 100, 300, "G", "A", "0.00106", ""),
            row("21", 100, 300, "G", "T", "0.91234", ".").replace("\t.\r", "\r"),
            row("21", 200, 400, "C", "T", ".", "."),
            row("21", 250, 450, "A", "G", "", "0.874"),
        ];
        let mut lines: Vec<&str> = vec![HEADER];
        lines.extend(rows.iter().map(String::as_str));
        let file = write_tabix_fixture(
            dir.path(),
            "synthetic_scores_empty_GRCh37.gz",
            &lines,
            IndexSpec::standard().bare_header(),
        );
        let plugin = load(Some(Assembly::Grch37), &[&file.to_string_lossy()]).unwrap();
        for alt in ["A", "T"] {
            let mut absent = variant(100, "G", alt);
            assert_eq!(
                scores_of(&plugin, &mut absent),
                [
                    (
                        "ada_score".to_string(),
                        if alt == "A" { "0.00106" } else { "0.91234" }.to_string()
                    ),
                    ("rf_score".to_string(), String::new())
                ],
                "{alt}"
            );
            assert!(
                absent.transcript_consequences[0]
                    .plugin_empty_strings
                    .is_empty(),
                "{alt}"
            );
        }
        assert!(scores_of(&plugin, &mut variant(200, "C", "T")).is_empty());

        let mut interior = variant(250, "A", "G");
        assert_eq!(
            scores_of(&plugin, &mut interior),
            [
                ("ada_score".to_string(), String::new()),
                ("rf_score".to_string(), "0.874".to_string())
            ]
        );
        assert_eq!(
            interior.transcript_consequences[0].plugin_empty_strings,
            ["ada_score".to_string()]
        );
    }

    /// Ensembl VEP 116.3 never hands dbscSNV a structural variant: the module
    /// keeps the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2229). A one-base `<DEL>` row, whose
    /// start equals its end as a substitution's does, gets no score, and nor
    /// does a scored substitution carrying the structural mark, which pins the
    /// gate rather than the span test; the plain substitution is scored.
    #[test]
    fn a_structural_variant_row_gets_no_score() {
        let dir = tempfile::tempdir().unwrap();
        let file = grch37_file(dir.path());
        let plugin = load(Some(Assembly::Grch37), &[&file]).unwrap();
        assert!(!plugin.annotates_structural_variants());
        let mut deletion = variant(100, "N", "<DEL>");
        deletion.is_structural = true;
        assert!(scores_of(&plugin, &mut deletion).is_empty());
        let mut flagged = variant(100, "G", "A");
        flagged.is_structural = true;
        assert!(scores_of(&plugin, &mut flagged).is_empty());
        assert_eq!(scores_of(&plugin, &mut variant(100, "G", "A")).len(), 2);
    }

    /// A row's position cell is read as `$tmp_data->{$pos_column} == $vf->{start}`
    /// (dbscSNV.pm:166) numifies it: `300.0`, `300x` and `3e2` are 300 and the
    /// row scores the variant at 300; `300.5`, `abc` and `0` are no position.
    /// The tabix index places a cell at its leading integer, so the `300.0`,
    /// `300.5` and `300x` shapes reach the match from a file and the exponent
    /// and the rest are read from records built here.
    #[test]
    fn a_position_cell_is_matched_as_perl_numifies_it() {
        let dir = tempfile::tempdir().unwrap();
        let cell_row = |hg38_pos: &str, alt: &str, ada: &str| {
            format!("21\t100\tG\t{alt}\t21\t{hg38_pos}\tY\tY\tintron\tGENE1\tintron\tNM_1:c.1A>G\tintron\tENSG1\tintron\tENST1:c.1A>G\t{ada}\t0.5\r")
        };
        let rows = [
            cell_row("300.0", "A", "0.1"),
            cell_row("300.5", "T", "0.2"),
            cell_row("300x", "C", "0.3"),
        ];
        let mut lines: Vec<&str> = vec![HEADER];
        lines.extend(rows.iter().map(String::as_str));
        let file = write_tabix_fixture(
            dir.path(),
            "synthetic_scores_GRCh38.gz",
            &lines,
            IndexSpec::standard().with_span_cols(5, 6, 6).bare_header(),
        )
        .to_string_lossy()
        .to_string();
        let plugin = load(Some(Assembly::Grch38), &[&file]).unwrap();
        assert_eq!(
            scores_of(&plugin, &mut variant(300, "G", "A")),
            [
                ("ada_score".to_string(), "0.1".to_string()),
                ("rf_score".to_string(), "0.5".to_string())
            ]
        );
        assert!(scores_of(&plugin, &mut variant(300, "G", "T")).is_empty());
        assert_eq!(scores_of(&plugin, &mut variant(300, "G", "C"))[0].1, "0.3");

        let record = |hg38_pos: &str| TabixRecord {
            chr: "21".to_string(),
            start: 300,
            end: 300,
            columns: cell_row(hg38_pos, "A", "0.9")
                .trim_end_matches('\r')
                .split('\t')
                .map(str::to_string)
                .collect(),
            header: None,
        };
        let scored = |hg38_pos: &str| {
            !plugin
                .scores(&variant(300, "G", "A"), &[&record(hg38_pos)])
                .values
                .is_empty()
        };
        assert!(scored("3e2"));
        assert!(scored(" 300"));
        assert!(!scored("abc"));
        assert!(!scored("0"));
        assert!(!scored("300.5"));
    }
}
