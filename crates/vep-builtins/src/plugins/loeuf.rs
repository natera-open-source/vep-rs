// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! LOEUF plugin: gnomAD's loss-of-function observed/expected upper bound
//! fraction, matched to the row's transcript or taken as the gene's minimum.
//!
//! The port of Ensembl VEP's `LOEUF.pm` (VEP_plugins release/116). The data
//! file is gnomAD's constraint table, tabix-indexed on its gene coordinate
//! columns (`:56-59` for the v2.1.1 table, `:65-66` for the v4 table), so a
//! region query returns every gene whose span overlaps the variant and the row
//! then picks the record whose transcript or gene id matches. The two gnomAD
//! layouts are told apart by the names at fixed header indices (`:84-93`,
//! `:122-143`); the coordinate columns are read from the tabix index itself,
//! the columns `tabix` queries on.
//!
//! **Usage:** `--plugin LOEUF,file=/path/to/loeuf_dataset.tsv.gz,match_by=gene`
//!
//! Output field: `LOEUF`, printed on transcript rows of sequence variants
//! (`variant_feature_types` keeps Ensembl VEP's default of `VariationFeature`),
//! never on an `upstream_gene_variant` or `downstream_gene_variant` row (`:169`).
//! `match_by=transcript` prints the matched record's value as written, `NA`
//! included (`:177-182`); `match_by=gene` prints the smallest numeric value over
//! the gene's records (`:184-205`).

use std::path::Path;

use indexmap::IndexMap;
use tracing::debug;

use crate::annotation_store::AnnotationStore;
use crate::params::PluginParams;
use crate::tabix::{
    check_plugin_table, find_index_path, load_index, split_field_count, table_contig,
    BatchQueryResult, TabixAnnotator, TabixConfig, TabixRecord,
};
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use noodles::csi::binning_index::index::header::format::CoordinateSystem;
use vep_core::assembly::Assembly;
use vep_core::consequence::{Consequence, TranscriptConsequence};
use vep_core::variant::InputVariant;

/// Output field name.
const OUTPUT_FIELD: &str = "LOEUF";

/// Where a gnomAD constraint table keeps the transcript id, the LOEUF value and
/// the gene id: 0-based column indices and the header names expected there
/// (`LOEUF.pm:84-93`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Layout {
    /// The gnomAD release the layout belongs to.
    pub version: &'static str,
    pub transcript: usize,
    pub loeuf: usize,
    pub gene_id: usize,
    /// The header names at `transcript`, `loeuf` and `gene_id`, in that order.
    pub names: [&'static str; 3],
}

/// gnomAD v2.1.1 `lof_metrics.by_transcript` table. The `by_gene` table holds
/// `oe_lof_upper` and `gene_id` one column earlier, so its header matches no
/// layout.
pub const V2_LAYOUT: Layout = Layout {
    version: "v2",
    transcript: 1,
    loeuf: 30,
    gene_id: 64,
    names: ["transcript", "oe_lof_upper", "gene_id"],
};

/// gnomAD v4 `constraint_metrics` table.
pub const V4_LAYOUT: Layout = Layout {
    version: "v4",
    transcript: 2,
    loeuf: 93,
    gene_id: 1,
    names: ["transcript", "lof.oe_ci.upper", "gene_id"],
};

const LAYOUTS: [Layout; 2] = [V2_LAYOUT, V4_LAYOUT];

impl Layout {
    /// The layout whose three names sit at its indices of `header`, v2 before
    /// v4 (`LOEUF.pm:126-141`).
    pub fn detect(header: &[String]) -> Option<Layout> {
        LAYOUTS
            .into_iter()
            .find(|layout| layout.missing(header).is_empty())
    }

    /// The column names a header is missing for this layout, as the plugin's
    /// refusal lists them.
    fn missing(&self, header: &[String]) -> Vec<&'static str> {
        [self.transcript, self.loeuf, self.gene_id]
            .iter()
            .zip(self.names)
            .filter(|(&index, name)| header.get(index).is_none_or(|cell| cell != name))
            .map(|(_, name)| name)
            .collect()
    }
}

/// Which identifier of the row a record must match (`LOEUF.pm:177, :184`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchBy {
    Gene,
    Transcript,
}

/// LOEUF gene-constraint plugin.
pub struct LoeufPlugin {
    /// Annotation store over the constraint table, set by `init`.
    store: Option<Box<dyn AnnotationStore>>,
    layout: Option<Layout>,
    match_by: Option<MatchBy>,
    /// The run's assembly, which names a mitochondrial record's contig to the
    /// table ([`table_contig`]).
    assembly: Option<Assembly>,
}

impl Default for LoeufPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl LoeufPlugin {
    /// Create an uninitialized LOEUF plugin.
    pub fn new() -> Self {
        Self {
            store: None,
            layout: None,
            match_by: None,
            assembly: None,
        }
    }

    /// The column configuration the file's tabix index was built with
    /// (`tabix -s`, `-b`, `-e`): the columns Ensembl VEP's `tabix` call
    /// queries on, whatever the table's own header says.
    fn indexed_config(path: &Path) -> Result<TabixConfig, PluginError> {
        let index_path = find_index_path(path)?;
        let index = load_index(&index_path)?;
        let header = index.header().ok_or_else(|| {
            PluginError::Init(format!(
                "LOEUF: index {} carries no tabix header, so its coordinate columns are \
                 unknown; index the file with tabix",
                index_path.display()
            ))
        })?;
        Ok(TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: header.reference_sequence_name_index(),
            start_col: header.start_position_index(),
            end_col: header.end_position_index(),
            ref_col: None,
            alt_col: None,
            zero_based: header.format().coordinate_system() == CoordinateSystem::Bed,
        })
    }

    /// The value the row prints for its record set: with `match_by=transcript`
    /// the first record naming the row's transcript, its LOEUF cell as written
    /// (`LOEUF.pm:177-182`), the empty string for an empty cell before a filled
    /// one (`true` beside it), which the array slice of `:218` keeps where the
    /// `split /\t/` of `:215` drops a trailing empty cell into the undefined
    /// value; with `match_by=gene` the smallest numeric LOEUF over the records
    /// naming the row's gene, as that record wrote it, ties to the earlier
    /// record (`:184-205`). `None` when nothing matches.
    fn value_for(
        &self,
        records: &[&TabixRecord],
        tc: &TranscriptConsequence,
    ) -> Option<(String, bool)> {
        let layout = self.layout?;
        match self.match_by? {
            MatchBy::Transcript => records
                .iter()
                .find(|record| cell(record, layout.transcript) == Some(tc.transcript_id()))
                .map(|record| {
                    let present = split_field_count(&record.columns);
                    let value = cell(record, layout.loeuf).filter(|_| layout.loeuf < present);
                    (value.unwrap_or("").to_string(), value == Some(""))
                }),
            MatchBy::Gene => {
                let mut minimum: Option<(f64, &str)> = None;
                for record in records {
                    if cell(record, layout.gene_id) != Some(tc.gene_id()) {
                        continue;
                    }
                    let Some(text) = cell(record, layout.loeuf) else {
                        continue;
                    };
                    let Ok(number) = text.parse::<f64>() else {
                        continue;
                    };
                    if minimum.is_none_or(|(least, _)| number < least) {
                        minimum = Some((number, text));
                    }
                }
                minimum.map(|(_, text)| (text.to_string(), false))
            }
        }
    }
}

/// The record's cell at a 0-based column, `None` past its last tab.
fn cell(record: &TabixRecord, index: usize) -> Option<&str> {
    record.columns.get(index).map(String::as_str)
}

/// Whether the row is an `upstream_gene_variant` or `downstream_gene_variant`
/// row, which gets no LOEUF (`LOEUF.pm:169`).
fn is_flanking_row(tc: &TranscriptConsequence) -> bool {
    tc.consequences.iter().any(|term| {
        matches!(
            term,
            Consequence::UpstreamGeneVariant | Consequence::DownstreamGeneVariant
        )
    })
}

impl BuiltinPlugin for LoeufPlugin {
    fn name(&self) -> &str {
        "LOEUF"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            OUTPUT_FIELD.to_string(),
            "Loss-of-function observed/expected upper bound fraction".to_string(),
        )]
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse("LOEUF", &["file", "match_by"], params)?;
        let file = params
            .get("file")
            .filter(|file| !file.is_empty())
            .ok_or_else(|| {
                PluginError::Init("LOEUF: file not provided; give file=<path>".into())
            })?;
        let path = Path::new(file);
        if !path.exists() {
            return Err(PluginError::Init(format!("LOEUF: file not found: {file}")));
        }
        let match_by = match params.get("match_by") {
            Some("gene") => MatchBy::Gene,
            Some("transcript") => MatchBy::Transcript,
            Some(other) => {
                return Err(PluginError::Init(format!(
                    "LOEUF: match_by must be gene or transcript, not '{other}'"
                )))
            }
            None => {
                return Err(PluginError::Init(
                    "LOEUF: argument 'match_by' is undefined; give match_by=gene or \
                     match_by=transcript"
                        .into(),
                ))
            }
        };

        // `add_file` runs before the index is read (LOEUF.pm:110).
        check_plugin_table(path)?;
        // Not `open_best_store`: the table's coordinate columns are read from
        // its tabix index (the `-s 114 -b 115 -e 116` of LOEUF.pm:66), which a
        // `.vpd` store beside the table does not carry.
        let store: Box<dyn AnnotationStore> = Box::new(
            TabixAnnotator::open_plugin_table(Self::indexed_config(path)?)?
                .on_assembly(self.assembly),
        );
        let header = store.header().ok_or_else(|| {
            PluginError::Init(format!("LOEUF: could not read headers from {file}"))
        })?;
        let layout = Layout::detect(header).ok_or_else(|| {
            let missing: Vec<String> = LAYOUTS
                .iter()
                .map(|layout| {
                    format!(
                        "{} layout: {}",
                        layout.version,
                        layout.missing(header).join(", ")
                    )
                })
                .collect();
            PluginError::Init(format!(
                "LOEUF: missing columns in {file} ({})",
                missing.join("; ")
            ))
        })?;

        debug!(
            file,
            layout = layout.version,
            match_by = ?match_by,
            "LOEUF plugin initialized"
        );
        self.store = Some(store);
        self.layout = Some(layout);
        self.match_by = Some(match_by);
        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        let Some(store) = &self.store else {
            return Ok(None);
        };
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|v| !v.is_structural)
            .map(|v| (table_contig(v, self.assembly), v.start, v.end))
            .collect();
        if regions.is_empty() {
            return Ok(None);
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
            let records = batch.lookup(
                table_contig(variant, self.assembly),
                variant.start,
                variant.end,
            );
            if records.is_empty() {
                continue;
            }
            for tc in &mut variant.transcript_consequences {
                if !self.annotates_row(tc) || is_flanking_row(tc) {
                    continue;
                }
                if let Some((value, empty_string)) = self.value_for(&records, tc) {
                    tc.plugin_data.insert(OUTPUT_FIELD.to_string(), value);
                    tc.plugin_empty_strings.retain(|key| key != OUTPUT_FIELD);
                    if empty_string {
                        tc.plugin_empty_strings.push(OUTPUT_FIELD.to_string());
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

    fn run(
        &self,
        _consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        // The rows are written in `annotate`, from one region query per batch.
        Ok(IndexMap::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_tabix_fixture, IndexSpec};
    use std::path::PathBuf;
    use std::sync::Arc;
    use vep_core::transcript::Transcript;

    /// A table row of `width` tab-separated cells, `.` except the given
    /// 0-based cells.
    fn row(width: usize, cells: &[(usize, &str)]) -> String {
        let mut fields = vec!["."; width];
        for &(index, value) in cells {
            fields[index] = value;
        }
        fields.join("\t")
    }

    /// The v2.1.1 table (78 columns, indexed `-s 76 -b 77 -e 78`) with SOD1's two
    /// transcripts as the GRCh37 table carries them (`oe_lof_upper` 9.7800e-01
    /// for ENST00000270142, 1.0000e+00 for ENST00000389995; gene 21:33031935-33041244)
    /// and a BACH1 transcript whose value is `NA`, plus a gene on another
    /// chromosome.
    fn v2_fixture(dir: &Path) -> PathBuf {
        let header = row(
            78,
            &[
                (0, "#gene"),
                (1, "transcript"),
                (30, "oe_lof_upper"),
                (64, "gene_id"),
                (75, "chromosome"),
                (76, "start_position"),
                (77, "end_position"),
            ],
        );
        let sod1 = |transcript: &str, loeuf: &str| {
            row(
                78,
                &[
                    (0, "SOD1"),
                    (1, transcript),
                    (30, loeuf),
                    (64, "ENSG00000142168"),
                    (75, "21"),
                    (76, "33031935"),
                    (77, "33041244"),
                ],
            )
        };
        let bach1 = row(
            78,
            &[
                (0, "BACH1"),
                (1, "ENST00000546469"),
                (30, "NA"),
                (64, "ENSG00000156273"),
                (75, "21"),
                (76, "30566392"),
                (77, "31003071"),
            ],
        );
        let other = row(
            78,
            &[
                (0, "OR4F5"),
                (1, "ENST00000335137"),
                (30, "1.8400e+00"),
                (64, "ENSG00000186092"),
                (75, "1"),
                (76, "69091"),
                (77, "70008"),
            ],
        );
        let lines = [
            header.as_str(),
            other.as_str(),
            bach1.as_str(),
            &sod1("ENST00000270142", "9.7800e-01"),
            &sod1("ENST00000389995", "1.0000e+00"),
        ];
        write_tabix_fixture(
            dir,
            "loeuf_dataset.tsv.gz",
            &lines,
            IndexSpec::standard().with_span_cols(76, 77, 78),
        )
    }

    /// The v4 table (116 columns, indexed `-s 114 -b 115 -e 116`) with SOD1's
    /// rows as the GRCh38 table carries them (`lof.oe_ci.upper` 9.5524e-01 for
    /// ENST00000389995, 9.0772e-01 for ENST00000270142 and for the RefSeq row
    /// NM_000454.5 under gene id 6647; gene 21:31659666-31668931).
    fn v4_fixture(dir: &Path) -> PathBuf {
        let header = row(
            116,
            &[
                (0, "#gene"),
                (1, "gene_id"),
                (2, "transcript"),
                (93, "lof.oe_ci.upper"),
                (113, "chrom"),
                (114, "start"),
                (115, "end"),
            ],
        );
        let sod1 = |gene_id: &str, transcript: &str, loeuf: &str, start: &str| {
            row(
                116,
                &[
                    (0, "SOD1"),
                    (1, gene_id),
                    (2, transcript),
                    (93, loeuf),
                    (113, "21"),
                    (114, start),
                    (115, "31668931"),
                ],
            )
        };
        let lines = [
            header.as_str(),
            &sod1(
                "ENSG00000142168",
                "ENST00000389995",
                "9.5524e-01",
                "31659666",
            ),
            &sod1("6647", "NM_000454.5", "9.0772e-01", "31659693"),
            &sod1(
                "ENSG00000142168",
                "ENST00000270142",
                "9.0772e-01",
                "31659693",
            ),
        ];
        write_tabix_fixture(
            dir,
            "loeuf_dataset_grch38.tsv.gz",
            &lines,
            IndexSpec::standard().with_span_cols(114, 115, 116),
        )
    }

    fn plugin(file: &Path, match_by: &str) -> LoeufPlugin {
        let mut plugin = LoeufPlugin::new();
        plugin
            .init(&[
                format!("file={}", file.display()),
                format!("match_by={match_by}"),
            ])
            .unwrap();
        plugin
    }

    fn consequence(transcript: &str, gene: &str, terms: &[Consequence]) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: transcript.into(),
                gene_stable_id: gene.into(),
                ..Transcript::default()
            })),
            consequences: terms.iter().copied().collect(),
            ..Default::default()
        }
    }

    /// A chromosome 21 SNV with one row per `(transcript, gene, terms)`.
    fn variant(pos: u64, rows: &[(&str, &str, &[Consequence])]) -> InputVariant {
        let mut v = InputVariant::new("21".into(), pos, pos, b"A".to_vec(), b"G".to_vec());
        for (transcript, gene, terms) in rows {
            v.transcript_consequences
                .push(consequence(transcript, gene, terms));
        }
        v
    }

    fn loeuf_of(variants: &[InputVariant], variant: usize, row: usize) -> Option<&str> {
        variants[variant].transcript_consequences[row]
            .plugin_data
            .get(OUTPUT_FIELD)
            .map(String::as_str)
    }

    const MISSENSE: &[Consequence] = &[Consequence::MissenseVariant];
    const SOD1: &str = "ENSG00000142168";

    #[test]
    fn plugin_name_and_header() {
        let plugin = LoeufPlugin::new();
        assert_eq!(plugin.name(), "LOEUF");
        assert_eq!(
            plugin.header_info(),
            [(
                "LOEUF".to_string(),
                "Loss-of-function observed/expected upper bound fraction".to_string()
            )]
        );
    }

    /// `LOEUF.pm:84-93`: the v2 table names `transcript`, `oe_lof_upper` and `gene_id` at
    /// indices 1, 30 and 64, the v4 table `transcript`, `lof.oe_ci.upper` and `gene_id`
    /// at 2, 93 and 1; the header is matched by index, never by searching for the names.
    #[test]
    fn v2_and_v4_header_layouts_are_detected_by_column_index() {
        let dir = tempfile::tempdir().unwrap();
        let v2 = plugin(&v2_fixture(dir.path()), "gene");
        assert_eq!(v2.layout, Some(V2_LAYOUT));
        let v4 = plugin(&v4_fixture(dir.path()), "gene");
        assert_eq!(v4.layout, Some(V4_LAYOUT));

        // The v2 names at the v4 indices, `transcript` and `gene_id` swapped, are no
        // layout: the refusal names the columns each layout misses (`:142-143`).
        let shuffled = row(
            116,
            &[
                (0, "#gene"),
                (1, "transcript"),
                (2, "gene_id"),
                (93, "oe_lof_upper"),
                (113, "chrom"),
                (114, "start"),
                (115, "end"),
            ],
        );
        let data = row(
            116,
            &[
                (0, "SOD1"),
                (113, "21"),
                (114, "31659666"),
                (115, "31668931"),
            ],
        );
        let file = write_tabix_fixture(
            dir.path(),
            "shuffled.tsv.gz",
            &[shuffled.as_str(), data.as_str()],
            IndexSpec::standard().with_span_cols(114, 115, 116),
        );
        let error = LoeufPlugin::new()
            .init(&[format!("file={}", file.display()), "match_by=gene".into()])
            .unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "plugin initialization failed: LOEUF: missing columns in {} \
                 (v2 layout: oe_lof_upper, gene_id; v4 layout: transcript, lof.oe_ci.upper, gene_id)",
                file.display()
            )
        );
    }

    /// `LOEUF.pm:110` hands the file to the tabix base class, which needs
    /// `<file>.tbi` beside it: an index named on the stem is a refusal
    /// naming the missing `.tbi`, though htslib would find it.
    #[test]
    fn table_without_its_tbi_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let file = v2_fixture(dir.path());
        let tbi = PathBuf::from(format!("{}.tbi", file.display()));
        std::fs::rename(&tbi, file.with_extension("tbi")).unwrap();
        let error = LoeufPlugin::new()
            .init(&[format!("file={}", file.display()), "match_by=gene".into()])
            .unwrap_err()
            .to_string();
        assert!(
            error.contains(&format!("{} not found", tbi.display())),
            "{error}"
        );
    }

    /// `LOEUF.pm:184-205`: under `match_by=gene` every record of the span whose
    /// gene id is the row's is read and the smallest numeric value is printed as
    /// that record wrote it; the gene's rows all print it, whichever transcript
    /// they are on, a non-numeric `NA` record is passed over, and the RefSeq
    /// record under another gene id does not count.
    #[test]
    fn match_by_gene_reports_the_minimum_loeuf_across_transcripts() {
        let dir = tempfile::tempdir().unwrap();
        let v2 = plugin(&v2_fixture(dir.path()), "gene");
        let mut variants = vec![
            variant(
                33035000,
                &[
                    ("ENST00000389995", SOD1, MISSENSE),
                    ("ENST00000270142", SOD1, MISSENSE),
                ],
            ),
            variant(
                30600000,
                &[("ENST00000546469", "ENSG00000156273", MISSENSE)],
            ),
        ];
        v2.run_batch(&mut variants).unwrap();
        assert_eq!(loeuf_of(&variants, 0, 0), Some("9.7800e-01"));
        assert_eq!(loeuf_of(&variants, 0, 1), Some("9.7800e-01"));
        assert_eq!(
            loeuf_of(&variants, 1, 0),
            None,
            "a gene whose only record is NA"
        );

        let v4 = plugin(&v4_fixture(dir.path()), "gene");
        let mut variants = vec![variant(
            31660000,
            &[
                ("ENST00000389995", SOD1, MISSENSE),
                ("ENST00000270142", SOD1, MISSENSE),
            ],
        )];
        v4.run_batch(&mut variants).unwrap();
        assert_eq!(loeuf_of(&variants, 0, 0), Some("9.0772e-01"));
        assert_eq!(loeuf_of(&variants, 0, 1), Some("9.0772e-01"));
    }

    /// `LOEUF.pm:177-182`: under `match_by=transcript` the row prints the value
    /// of the record whose transcript id equals its stable id, as written, so a
    /// transcript whose record reads `NA` prints `NA`, and a transcript the
    /// table lacks prints nothing.
    #[test]
    fn match_by_transcript_matches_the_stable_id() {
        let dir = tempfile::tempdir().unwrap();
        let v4 = plugin(&v4_fixture(dir.path()), "transcript");
        let mut variants = vec![variant(
            31660000,
            &[
                ("ENST00000389995", SOD1, MISSENSE),
                ("ENST00000270142", SOD1, MISSENSE),
                ("ENST00000999999", SOD1, MISSENSE),
            ],
        )];
        v4.run_batch(&mut variants).unwrap();
        assert_eq!(loeuf_of(&variants, 0, 0), Some("9.5524e-01"));
        assert_eq!(loeuf_of(&variants, 0, 1), Some("9.0772e-01"));
        assert_eq!(loeuf_of(&variants, 0, 2), None);

        let v2 = plugin(&v2_fixture(dir.path()), "transcript");
        let mut variants = vec![
            variant(33035000, &[("ENST00000389995", SOD1, MISSENSE)]),
            variant(
                30600000,
                &[("ENST00000546469", "ENSG00000156273", MISSENSE)],
            ),
        ];
        v2.run_batch(&mut variants).unwrap();
        assert_eq!(loeuf_of(&variants, 0, 0), Some("1.0000e+00"));
        assert_eq!(loeuf_of(&variants, 1, 0), Some("NA"));
    }

    /// Under `match_by=transcript` the matched record's LOEUF cell is written as
    /// it stands (`LOEUF.pm:177-182`, `:218`): an empty cell before a filled one
    /// is the empty string, named in the row's `plugin_empty_strings`, so the
    /// tab format writes an empty cell and the JSON format `""`; under
    /// `match_by=gene` such a record is not a number and takes part in no
    /// minimum (`:192`).
    #[test]
    fn an_empty_loeuf_cell_before_a_filled_one_is_the_empty_string() {
        let dir = tempfile::tempdir().unwrap();
        let header = row(
            116,
            &[
                (0, "#gene"),
                (1, "gene_id"),
                (2, "transcript"),
                (93, "lof.oe_ci.upper"),
                (113, "chrom"),
                (114, "start"),
                (115, "end"),
            ],
        );
        let sod1 = |transcript: &str, loeuf: &str| {
            row(
                116,
                &[
                    (0, "SOD1"),
                    (1, SOD1),
                    (2, transcript),
                    (93, loeuf),
                    (113, "21"),
                    (114, "31659666"),
                    (115, "31668931"),
                ],
            )
        };
        let lines = [
            header.as_str(),
            &sod1("ENST00000389995", ""),
            &sod1("ENST00000270142", "9.0772e-01"),
        ];
        let path = write_tabix_fixture(
            dir.path(),
            "loeuf_empty_cell.tsv.gz",
            &lines,
            IndexSpec::standard().with_span_cols(114, 115, 116),
        );
        let by_transcript = plugin(&path, "transcript");
        let mut variants = vec![variant(
            31660000,
            &[
                ("ENST00000389995", SOD1, MISSENSE),
                ("ENST00000270142", SOD1, MISSENSE),
            ],
        )];
        by_transcript.run_batch(&mut variants).unwrap();
        assert_eq!(loeuf_of(&variants, 0, 0), Some(""));
        assert_eq!(
            variants[0].transcript_consequences[0].plugin_empty_strings,
            [OUTPUT_FIELD]
        );
        assert_eq!(loeuf_of(&variants, 0, 1), Some("9.0772e-01"));
        assert!(variants[0].transcript_consequences[1]
            .plugin_empty_strings
            .is_empty());

        let by_gene = plugin(&path, "gene");
        let mut variants = vec![variant(31660000, &[("ENST00000389995", SOD1, MISSENSE)])];
        by_gene.run_batch(&mut variants).unwrap();
        assert_eq!(loeuf_of(&variants, 0, 0), Some("9.0772e-01"));
        assert!(variants[0].transcript_consequences[0]
            .plugin_empty_strings
            .is_empty());
    }

    /// `LOEUF.pm:169`: a row carrying `upstream_gene_variant` or
    /// `downstream_gene_variant` gets no value, under either `match_by`, even
    /// when the variant lies inside the table's gene span (the table's span is
    /// the gene's, which can be wider than the transcript's). A structural
    /// variant is never handed to the plugin (`BaseVepPlugin.pm` keeps
    /// `variant_feature_types` at `VariationFeature`).
    #[test]
    fn upstream_and_downstream_rows_get_no_loeuf() {
        let dir = tempfile::tempdir().unwrap();
        let file = v4_fixture(dir.path());
        for match_by in ["gene", "transcript"] {
            let plugin = plugin(&file, match_by);
            let mut variants = vec![variant(
                31660000,
                &[
                    ("ENST00000389995", SOD1, &[Consequence::UpstreamGeneVariant]),
                    (
                        "ENST00000270142",
                        SOD1,
                        &[Consequence::DownstreamGeneVariant],
                    ),
                    (
                        "ENST00000270142",
                        SOD1,
                        &[Consequence::IntronVariant, Consequence::UpstreamGeneVariant],
                    ),
                    ("ENST00000270142", SOD1, &[Consequence::IntronVariant]),
                ],
            )];
            plugin.run_batch(&mut variants).unwrap();
            assert_eq!(loeuf_of(&variants, 0, 0), None, "{match_by}");
            assert_eq!(loeuf_of(&variants, 0, 1), None, "{match_by}");
            assert_eq!(loeuf_of(&variants, 0, 2), None, "{match_by}");
            assert_eq!(loeuf_of(&variants, 0, 3), Some("9.0772e-01"), "{match_by}");
        }

        let plugin = plugin(&file, "gene");
        let mut deletion = variant(31659000, &[("ENST00000270142", SOD1, MISSENSE)]);
        deletion.end = 31670000;
        deletion.is_structural = true;
        let mut variants = vec![deletion];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(loeuf_of(&variants, 0, 0), None);
    }

    /// `LOEUF.pm:146-153` makes `match_by` mandatory and `:207-209` silently
    /// annotates nothing under any other value; vep-rs refuses both at
    /// initialization, as it refuses a missing or absent file (`:109`).
    #[test]
    fn match_by_outside_gene_or_transcript_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let file = v4_fixture(dir.path());
        let init = |params: &[String]| LoeufPlugin::new().init(params).unwrap_err().to_string();
        let file_param = format!("file={}", file.display());
        assert_eq!(
            init(&[file_param.clone(), "match_by=exon".into()]),
            "plugin initialization failed: LOEUF: match_by must be gene or transcript, not 'exon'"
        );
        assert_eq!(
            init(std::slice::from_ref(&file_param)),
            "plugin initialization failed: LOEUF: argument 'match_by' is undefined; give \
             match_by=gene or match_by=transcript"
        );
        assert_eq!(
            init(&["match_by=gene".into()]),
            "plugin initialization failed: LOEUF: file not provided; give file=<path>"
        );
        assert_eq!(
            init(&[file.display().to_string(), "gene".into()]),
            "plugin initialization failed: LOEUF: file not provided; give file=<path>"
        );
        assert_eq!(
            init(&["file=/no/such/loeuf.tsv.gz".into(), "match_by=gene".into()]),
            "plugin initialization failed: LOEUF: file not found: /no/such/loeuf.tsv.gz"
        );
        assert_eq!(
            init(&[file_param, "match_by=gene".into(), "cover=1".into()]),
            "plugin initialization failed: LOEUF: unknown parameter 'cover'; LOEUF accepts \
             file, match_by"
        );
    }

    #[test]
    fn uninitialized_plugin_annotates_nothing() {
        let plugin = LoeufPlugin::new();
        let mut variants = vec![variant(31660000, &[("ENST00000270142", SOD1, MISSENSE)])];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(loeuf_of(&variants, 0, 0), None);
        let tc = consequence("ENST00000270142", SOD1, &[Consequence::MissenseVariant]);
        assert!(plugin.run(&tc, &variants[0]).unwrap().is_empty());
    }
}
