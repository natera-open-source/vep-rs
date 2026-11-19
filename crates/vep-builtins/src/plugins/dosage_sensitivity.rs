// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! DosageSensitivity plugin: the probabilities that a gene is haploinsufficient
//! (`pHaplo`) and triplosensitive (`pTriplo`) from Collins et al. 2022.
//!
//! The port of Ensembl VEP's `DosageSensitivity.pm` (VEP_plugins release/116).
//! The whole three-column table (`gene`, `pHaplo`, `pTriplo`) is read into a
//! symbol-keyed map at initialization (`:127-146`), and each transcript row
//! prints the two values of its gene symbol (`:169, :180`), or of the symbol
//! fourteen GRCh38 gene stable ids map to where the annotation and the table
//! name the gene differently (`:74-89, :170-173`).
//!
//! **Usage:** `--plugin DosageSensitivity,file=/path/to/dosage_sensitivity_scores.tsv.gz[,cover=1]`
//!
//! Output fields: `pHaplo`, `pTriplo`, on transcript rows of sequence and
//! structural variants (`variant_feature_types` names both, `:113-115`;
//! `feature_types` is `Feature`, `:117-119`, and only a transcript row carries
//! the gene the lookup reads). Under `cover=1` a row prints them only when the
//! variant spans the whole transcript (`:148-164, :176-178`); `cover` is on for
//! the string `1` alone, as `$param_hash->{cover} eq '1'` (`:101`) reads it.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use indexmap::IndexMap;
use tracing::debug;

use crate::params::PluginParams;
use crate::tabix::BoundedLines;
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// Output field names (`DosageSensitivity.pm:103-106`).
const P_HAPLO: &str = "pHaplo";
const P_TRIPLO: &str = "pTriplo";

/// The gene stable ids whose symbol in the annotation differs from the clone-based
/// symbol the table carries (`DosageSensitivity.pm:74-89`).
const SYMBOL_REMAP: [(&str, &str); 14] = [
    ("ENSG00000235478", "AC006946.15"),
    ("ENSG00000218672", "AC008060.7"),
    ("ENSG00000177340", "AC024940.1"),
    ("ENSG00000222007", "AC064874.1"),
    ("ENSG00000228919", "AC097381.1"),
    ("ENSG00000268218", "AC137932.1"),
    ("ENSG00000267270", "AC139100.2"),
    ("ENSG00000245317", "CTC-241N9.1"),
    ("ENSG00000264278", "RP11-162A12.2"),
    ("ENSG00000272297", "RP11-215A19.2"),
    ("ENSG00000271949", "RP11-302M6.4"),
    ("ENSG00000204398", "RP11-65D24.2"),
    ("ENSG00000267127", "RP11-795F19.5"),
    ("ENSG00000259658", "RP11-89K11.1"),
];

/// A gene's `pHaplo` and `pTriplo` cells as its table line carries them:
/// `None` for a cell the line lacks, Ensembl VEP's undefined value, and
/// `Some("")` for an empty one, its empty string.
type Cells = (Option<String>, Option<String>);

/// DosageSensitivity gene-score plugin.
pub struct DosageSensitivityPlugin {
    /// Gene symbol to its `(pHaplo, pTriplo)` cells as the table writes them.
    scores: HashMap<String, Cells>,
    /// `cover=1`: a row prints only when the variant spans its transcript. On
    /// for the string `1` alone, as `$param_hash->{cover} eq '1' ? 1 : 0`
    /// (`:101`) reads it; every other value, and no `cover`, is off.
    cover: bool,
}

impl Default for DosageSensitivityPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl DosageSensitivityPlugin {
    /// Create an uninitialized DosageSensitivity plugin.
    pub fn new() -> Self {
        Self {
            scores: HashMap::new(),
            cover: false,
        }
    }

    /// Read the table, gzip-compressed or plain, as `_process_DS_file` does
    /// (`DosageSensitivity.pm:127-146`): every line splits on tabs into gene,
    /// `pHaplo` and `pTriplo`, a line whose gene opens with `#` is skipped, and a
    /// gene listed twice keeps its last line. The three-variable assignment of
    /// `:138` keeps a gene whose line is short, each missing cell undefined, and
    /// an empty cell, trailing or not, as the empty string; the `chomp` of
    /// `:136` leaves a `\r` before the line end on the last cell. The lines are
    /// split as the bytes the module reads; a byte that is not UTF-8 in a kept
    /// cell is written as U+FFFD where the module prints the byte.
    fn load_scores(path: &str) -> Result<HashMap<String, Cells>, PluginError> {
        let file = std::fs::File::open(path).map_err(|e| {
            PluginError::Init(format!("DosageSensitivity: failed to open {path}: {e}"))
        })?;
        let mut reader = BufReader::new(file);
        let is_gzip = reader
            .fill_buf()
            .map_err(|e| {
                PluginError::Init(format!("DosageSensitivity: failed to read {path}: {e}"))
            })?
            .starts_with(&[0x1f, 0x8b]);
        let reader: Box<dyn BufRead> = if is_gzip {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(reader)))
        } else {
            Box::new(reader)
        };
        let mut scores = HashMap::new();
        let text = |cell: &[u8]| String::from_utf8_lossy(cell).into_owned();
        let mut lines = BoundedLines::chomped(reader, Path::new(path));
        while let Some(line) = lines.next_line() {
            let line = line?;
            let mut cells = line.split(|byte| *byte == b'\t');
            let gene = cells.next().unwrap_or_default();
            if gene.starts_with(b"#") {
                continue;
            }
            let p_haplo = cells.next().map(text);
            let p_triplo = cells.next().map(text);
            scores.insert(text(gene), (p_haplo, p_triplo));
        }
        Ok(scores)
    }

    /// The table's symbol for the row's gene: the row's own symbol when the
    /// table lists it, else the clone-based symbol its gene id maps to, else
    /// the row's symbol as it stands (`DosageSensitivity.pm:169-174`).
    fn table_symbol<'a>(&self, tc: &'a TranscriptConsequence) -> Option<&'a str> {
        let symbol = tc.gene_symbol();
        if symbol.is_some_and(|symbol| self.scores.contains_key(symbol)) {
            return symbol;
        }
        SYMBOL_REMAP
            .iter()
            .find(|(gene_id, _)| *gene_id == tc.gene_id())
            .map(|(_, clone_symbol)| *clone_symbol)
            .or(symbol)
    }

    /// The cells the row prints (`DosageSensitivity.pm:169-180`): its gene's
    /// table line, or none when the gene has no symbol or no line, or when
    /// `cover=1` and the variant does not span the transcript.
    fn gene_cells(&self, tc: &TranscriptConsequence, variant: &InputVariant) -> Option<&Cells> {
        let symbol = self.table_symbol(tc)?;
        if self.cover && !covers_transcript(variant, tc) {
            return None;
        }
        self.scores.get(symbol)
    }
}

/// The two fields in output order, each with its cell.
fn fields_of(cells: &Cells) -> [(&'static str, Option<&str>); 2] {
    [
        (P_HAPLO, cells.0.as_deref()),
        (P_TRIPLO, cells.1.as_deref()),
    ]
}

/// Whether the variant's span covers the row's transcript whole
/// (`DosageSensitivity.pm:163`: `vf.start <= feature.start && vf.end >= feature.end`,
/// with the record's span as written, so an insertion, whose start lies past
/// its end, never covers).
fn covers_transcript(variant: &InputVariant, tc: &TranscriptConsequence) -> bool {
    tc.feature_end >= tc.feature_start
        && tc.feature_end != 0
        && variant.start <= tc.feature_start
        && variant.end >= tc.feature_end
}

impl BuiltinPlugin for DosageSensitivityPlugin {
    fn name(&self) -> &str {
        "DosageSensitivity"
    }

    /// The Ensembl VEP module declares `['VariationFeature',
    /// 'StructuralVariationFeature']` (DosageSensitivity.pm:113-115).
    fn annotates_structural_variants(&self) -> bool {
        true
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            (
                P_HAPLO.to_string(),
                "Probability of haploinsufficiency (deletion intolerance) of the affected gene"
                    .to_string(),
            ),
            (
                P_TRIPLO.to_string(),
                "Probability of triplosensitivity (duplication intolerance) of the affected gene"
                    .to_string(),
            ),
        ]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse("DosageSensitivity", &["file", "cover"], params)?;
        let file = params.get("file").ok_or_else(|| {
            PluginError::Init(
                "DosageSensitivity: file must be specified, such as \
                 file=/FULL_PATH_TO/dosage_sensitivity_scores.tsv.gz"
                    .into(),
            )
        })?;
        self.cover = params.get("cover") == Some("1");
        self.scores = Self::load_scores(file)?;
        debug!(
            genes = self.scores.len(),
            path = file,
            cover = self.cover,
            "DosageSensitivity plugin loaded scores"
        );
        Ok(())
    }

    /// [`run`](Self::run) on every transcript row, naming in the row's
    /// `plugin_empty_strings` the field of an empty cell, whose empty value is
    /// Ensembl VEP's empty string where a missing cell's is its undefined value.
    /// A later plugin's value replaces an earlier one under the same key
    /// (OutputFactory.pm:2230-2231), its mark too.
    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        for variant in variants.iter_mut().filter(|v| self.annotates_variant(v)) {
            for idx in 0..variant.transcript_consequences.len() {
                let tc = &variant.transcript_consequences[idx];
                if !self.annotates_row(tc) {
                    continue;
                }
                let Some(cells) = self.gene_cells(tc, variant) else {
                    continue;
                };
                let tc = &mut variant.transcript_consequences[idx];
                tc.plugin_empty_strings
                    .retain(|key| key != P_HAPLO && key != P_TRIPLO);
                for (field, cell) in fields_of(cells) {
                    if cell == Some("") {
                        tc.plugin_empty_strings.push(field.to_string());
                    }
                    tc.plugin_data
                        .insert(field.to_string(), cell.unwrap_or_default().to_string());
                }
            }
        }
        Ok(())
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut result = IndexMap::new();
        if let Some(cells) = self.gene_cells(consequence, variant) {
            for (field, cell) in fields_of(cells) {
                result.insert(field.to_string(), cell.unwrap_or_default().to_string());
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;
    use std::sync::Arc;
    use vep_core::transcript::Transcript;

    /// Rows of the Collins et al. table as published (`#gene`, `pHaplo`,
    /// `pTriplo`), with the clone-named gene the stable id ENSG00000235478 maps to.
    const TABLE: &str = "#gene\tpHaplo\tpTriplo\n\
        CACNA1C\t0.99898184581082\t1\n\
        ZNF462\t1\t0.987995995573708\n\
        SOD1\t0.0587254800331713\t0.0117470702227125\n\
        AC006946.15\t0.979788404622403\t0.999999999973713\n";

    fn write_table(dir: &Path, name: &str, gzip: bool) -> PathBuf {
        let path = dir.join(name);
        if gzip {
            let mut writer = flate2::write::GzEncoder::new(
                std::fs::File::create(&path).unwrap(),
                flate2::Compression::default(),
            );
            writer.write_all(TABLE.as_bytes()).unwrap();
            writer.finish().unwrap();
        } else {
            std::fs::write(&path, TABLE).unwrap();
        }
        path
    }

    fn plugin(dir: &Path, cover: Option<&str>) -> DosageSensitivityPlugin {
        let file = write_table(
            dir,
            "Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz",
            true,
        );
        let mut params = vec![format!("file={}", file.display())];
        if let Some(cover) = cover {
            params.push(format!("cover={cover}"));
        }
        let mut plugin = DosageSensitivityPlugin::new();
        plugin.init(&params).unwrap();
        plugin
    }

    /// A plugin loaded from the plain-text table `lines`.
    fn plugin_of_table(dir: &Path, lines: &str) -> DosageSensitivityPlugin {
        let path = dir.join("scores.tsv");
        std::fs::write(&path, lines).unwrap();
        let mut plugin = DosageSensitivityPlugin::new();
        plugin.init(&[format!("file={}", path.display())]).unwrap();
        plugin
    }

    /// The two fields and the `plugin_empty_strings` marks `run_batch` leaves
    /// on a transcript row of `symbol`.
    fn batch_fields(
        plugin: &DosageSensitivityPlugin,
        symbol: &str,
    ) -> (Vec<(String, String)>, Vec<String>) {
        let mut variants = vec![snv(150)];
        variants[0]
            .transcript_consequences
            .push(row(Some(symbol), "ENSG00000000001", 100, 200));
        plugin.run_batch(&mut variants).unwrap();
        let tc = &variants[0].transcript_consequences[0];
        (
            tc.plugin_data
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect(),
            tc.plugin_empty_strings.clone(),
        )
    }

    fn fields(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect()
    }

    /// A transcript row of `gene` at `symbol`, spanning `start..=end`.
    fn row(symbol: Option<&str>, gene_id: &str, start: u64, end: u64) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: "ENST00000270142".into(),
                gene_stable_id: gene_id.into(),
                gene_symbol: symbol.map(Into::into),
                ..Transcript::default()
            })),
            feature_start: start,
            feature_end: end,
            ..Default::default()
        }
    }

    fn snv(pos: u64) -> InputVariant {
        InputVariant::new("21".into(), pos, pos, b"A".to_vec(), b"G".to_vec())
    }

    fn scores_of(result: &IndexMap<String, String>) -> Option<(&str, &str)> {
        Some((
            result.get(P_HAPLO)?.as_str(),
            result.get(P_TRIPLO)?.as_str(),
        ))
    }

    const SOD1_SCORES: (&str, &str) = ("0.0587254800331713", "0.0117470702227125");

    /// `DosageSensitivity.pm:98` dies without `file=`; a positional path is no
    /// `file=` either, since `params_to_hash` reads a bare token as a key with no
    /// value.
    #[test]
    fn file_is_mandatory() {
        let dir = tempfile::tempdir().unwrap();
        let file = write_table(dir.path(), "scores.tsv.gz", true);
        let message = "plugin initialization failed: DosageSensitivity: file must be \
                       specified, such as file=/FULL_PATH_TO/dosage_sensitivity_scores.tsv.gz";
        for params in [
            vec![],
            vec![file.display().to_string()],
            vec!["cover=1".to_string()],
        ] {
            let error = DosageSensitivityPlugin::new().init(&params).unwrap_err();
            assert_eq!(error.to_string(), message, "{params:?}");
        }
        let error = DosageSensitivityPlugin::new()
            .init(&["file=/no/such/scores.tsv.gz".to_string()])
            .unwrap_err()
            .to_string();
        assert!(
            error.starts_with(
                "plugin initialization failed: DosageSensitivity: failed to open /no/such/scores.tsv.gz: "
            ),
            "{error}"
        );
    }

    /// `DosageSensitivity.pm:135-143`: the header, whose gene cell opens with
    /// `#`, is skipped and every other line is keyed by its first cell; the
    /// gzip-compressed file and the plain text file load alike
    /// (`get_compressed_filehandle`).
    #[test]
    fn hash_lines_are_skipped_and_symbols_keyed() {
        let dir = tempfile::tempdir().unwrap();
        for (name, gzip) in [("scores.tsv.gz", true), ("scores.tsv", false)] {
            let file = write_table(dir.path(), name, gzip);
            let mut plugin = DosageSensitivityPlugin::new();
            plugin.init(&[format!("file={}", file.display())]).unwrap();
            assert_eq!(plugin.scores.len(), 4, "{name}");
            assert!(!plugin.scores.contains_key("#gene"), "{name}");
            let cells = |p_haplo: &str, p_triplo: &str| {
                Some((Some(p_haplo.to_string()), Some(p_triplo.to_string())))
            };
            assert_eq!(
                plugin.scores.get("CACNA1C").cloned(),
                cells("0.99898184581082", "1"),
                "{name}"
            );
            assert_eq!(
                plugin.scores.get("ZNF462").cloned(),
                cells("1", "0.987995995573708"),
                "{name}"
            );
        }
    }

    /// A byte that is not UTF-8 (Latin-1 `é`, 0xE9) is a byte of its cell, as
    /// `split /\t/` reads the bytes of a line (`DosageSensitivity.pm:135-138`):
    /// in a `#` line and in a fourth cell the table loads unchanged, in the gene
    /// cell it names another gene, and in a score cell the module prints the
    /// byte, written here as U+FFFD.
    #[test]
    fn a_byte_that_is_not_utf8_is_a_byte_of_its_cell() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.tsv");
        std::fs::write(
            &path,
            b"#gene\tpHaplo\tpTriplo\n#caf\xE9\t1\t1\nSOD1\t0.0587254800331713\t0.0117470702227125\tcaf\xE9\n\
              AP\xE9P\t0.9\t0.8\nZNF462\t1\xE9\t0.987995995573708\n",
        )
        .unwrap();
        let mut plugin = DosageSensitivityPlugin::new();
        plugin.init(&[format!("file={}", path.display())]).unwrap();
        assert_eq!(plugin.scores.len(), 3);
        let sod1 = row(Some("SOD1"), "ENSG00000142168", 100, 200);
        assert_eq!(
            scores_of(&plugin.run(&sod1, &snv(150)).unwrap()),
            Some(SOD1_SCORES)
        );
        let app = row(Some("APP"), "ENSG00000142192", 100, 200);
        assert!(plugin.run(&app, &snv(150)).unwrap().is_empty());
        let znf462 = row(Some("ZNF462"), "ENSG00000148143", 100, 200);
        assert_eq!(
            scores_of(&plugin.run(&znf462, &snv(150)).unwrap()),
            Some(("1\u{FFFD}", "0.987995995573708"))
        );
    }

    /// `DosageSensitivity.pm:138-142` assigns the split line to three
    /// variables and stores both keys for the gene whatever the line holds, so a
    /// line of one or two cells keeps its gene with each missing score
    /// undefined: `run` (`:180`) returns both keys, the default format writes
    /// `pTriplo=` (`VEP_output.pm:117-119`), the tab format `-` and the JSON
    /// format no key (`Tab.pm:128`, `OutputFactory/JSON.pm:327`). The row
    /// carries both fields, the missing one empty and unmarked; a four-cell
    /// line keeps its first three cells, as the three-variable assignment does.
    #[test]
    fn a_line_lacking_a_score_cell_keeps_the_gene_with_the_empty_value() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin_of_table(
            dir.path(),
            "#gene\tpHaplo\tpTriplo\nSOD1\t0.5\nBACH1\nKRTAP20-4\t0.1\t0.2\textra\n",
        );
        assert_eq!(plugin.scores.len(), 3);
        assert_eq!(
            batch_fields(&plugin, "SOD1"),
            (fields(&[("pHaplo", "0.5"), ("pTriplo", "")]), vec![]),
            "two cells"
        );
        assert_eq!(
            batch_fields(&plugin, "BACH1"),
            (fields(&[("pHaplo", ""), ("pTriplo", "")]), vec![]),
            "one cell"
        );
        assert_eq!(
            batch_fields(&plugin, "KRTAP20-4"),
            (fields(&[("pHaplo", "0.1"), ("pTriplo", "0.2")]), vec![]),
            "four cells"
        );
        assert_eq!(
            batch_fields(&plugin, "APP"),
            (vec![], vec![]),
            "not in the table"
        );
        let short = row(Some("SOD1"), "ENSG00000142168", 100, 200);
        assert_eq!(
            scores_of(&plugin.run(&short, &snv(150)).unwrap()),
            Some(("0.5", ""))
        );
    }

    /// `DosageSensitivity.pm:138`: the three-variable list assignment supplies
    /// `split` a limit, so an empty cell, interior or trailing, is the empty
    /// string, which the tab format writes as an empty cell and the JSON
    /// format as `""` (`Tab.pm:128`, `OutputFactory/JSON.pm:327`); the row
    /// names that field in `plugin_empty_strings` and leaves a missing cell's
    /// field unnamed. A later value under the same key replaces the mark
    /// (OutputFactory.pm:2230-2231).
    #[test]
    fn an_empty_cell_is_the_empty_string() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin_of_table(
            dir.path(),
            "#gene\tpHaplo\tpTriplo\nAPP\t\t0.3\nKRTAP20-4\t0.1\t\nCECR6\t\t\n",
        );
        assert_eq!(
            batch_fields(&plugin, "APP"),
            (
                fields(&[("pHaplo", ""), ("pTriplo", "0.3")]),
                vec!["pHaplo".to_string()]
            ),
            "an empty first cell"
        );
        assert_eq!(
            batch_fields(&plugin, "KRTAP20-4"),
            (
                fields(&[("pHaplo", "0.1"), ("pTriplo", "")]),
                vec!["pTriplo".to_string()]
            ),
            "an empty trailing cell"
        );
        assert_eq!(
            batch_fields(&plugin, "CECR6"),
            (
                fields(&[("pHaplo", ""), ("pTriplo", "")]),
                vec!["pHaplo".to_string(), "pTriplo".to_string()]
            ),
            "both cells empty"
        );
        let mut variants = vec![snv(150)];
        let mut marked = row(Some("KRTAP20-4"), "ENSG00000206105", 100, 200);
        marked.plugin_empty_strings = vec!["pHaplo".to_string(), "Other".to_string()];
        variants[0].transcript_consequences.push(marked);
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(
            variants[0].transcript_consequences[0].plugin_empty_strings,
            vec!["Other".to_string(), "pTriplo".to_string()],
            "an earlier mark under the plugin's keys is replaced"
        );
    }

    /// `DosageSensitivity.pm:136` chomps the `\n` alone, so a table whose lines
    /// end in `\r\n` keeps the `\r` on the last cell of each line (`:138`): a
    /// score cell reads `0.4\r`, and an empty trailing cell reads `\r`, a
    /// non-empty string that is not marked.
    #[test]
    fn a_carriage_return_stays_on_the_last_cell_as_chomp_leaves_it() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin_of_table(
            dir.path(),
            "#gene\tpHaplo\tpTriplo\r\nSOD1\t0.3\t0.4\r\nAPP\t0.1\t\r\nBACH1\t0.2\r\n",
        );
        assert_eq!(plugin.scores.len(), 3);
        assert_eq!(
            batch_fields(&plugin, "SOD1"),
            (fields(&[("pHaplo", "0.3"), ("pTriplo", "0.4\r")]), vec![]),
            "a full line"
        );
        assert_eq!(
            batch_fields(&plugin, "APP"),
            (fields(&[("pHaplo", "0.1"), ("pTriplo", "\r")]), vec![]),
            "an empty trailing cell before the line end"
        );
        assert_eq!(
            batch_fields(&plugin, "BACH1"),
            (fields(&[("pHaplo", "0.2\r"), ("pTriplo", "")]), vec![]),
            "a two-cell line"
        );
    }

    /// `DosageSensitivity.pm:74-89, :170-173`: a row whose symbol the table
    /// lacks is looked up by its gene stable id in the fourteen-entry map, each
    /// id resolving to the table's clone-based symbol; a row whose own symbol
    /// the table carries keeps it, and a row with neither prints nothing.
    #[test]
    fn stable_id_remap_covers_the_fourteen_clone_symbols() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(dir.path(), None);
        assert_eq!(SYMBOL_REMAP.len(), 14);
        for (gene_id, clone_symbol) in SYMBOL_REMAP {
            let tc = row(Some("AC006946.1"), gene_id, 100, 200);
            assert_eq!(plugin.table_symbol(&tc), Some(clone_symbol), "{gene_id}");
            let tc = row(None, gene_id, 100, 200);
            assert_eq!(plugin.table_symbol(&tc), Some(clone_symbol), "{gene_id}");
        }
        let remapped = row(Some("AC006946.1"), "ENSG00000235478", 100, 200);
        assert_eq!(
            scores_of(&plugin.run(&remapped, &snv(150)).unwrap()),
            Some(("0.979788404622403", "0.999999999973713"))
        );
        let own = row(Some("SOD1"), "ENSG00000142168", 100, 200);
        assert_eq!(
            scores_of(&plugin.run(&own, &snv(150)).unwrap()),
            Some(SOD1_SCORES)
        );
        let unknown = row(Some("TTN"), "ENSG00000155657", 100, 200);
        assert!(plugin.run(&unknown, &snv(150)).unwrap().is_empty());
        let nameless = row(None, "ENSG00000155657", 100, 200);
        assert!(plugin.run(&nameless, &snv(150)).unwrap().is_empty());
    }

    /// `DosageSensitivity.pm:159-163, :176-178`: under `cover=1` the row prints
    /// only when the variant's start is at or before the transcript's start and
    /// its end at or after the transcript's end; a variant inside the
    /// transcript, or reaching one end alone, prints nothing. `cover=0` and no
    /// `cover` print regardless.
    #[test]
    fn cover_requires_the_allele_to_span_the_transcript() {
        let dir = tempfile::tempdir().unwrap();
        let tc = row(Some("SOD1"), "ENSG00000142168", 31659666, 31668931);
        let inside = snv(31660000);
        let spanning = InputVariant::new(
            "21".into(),
            31659000,
            31670000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let mut exact = spanning.clone();
        exact.start = 31659666;
        exact.end = 31668931;
        let mut short_at_end = spanning.clone();
        short_at_end.end = 31668930;
        let mut insertion = snv(31659666);
        insertion.end = 31659665;

        let cover = plugin(dir.path(), Some("1"));
        assert!(cover.run(&tc, &inside).unwrap().is_empty());
        assert!(cover.run(&tc, &short_at_end).unwrap().is_empty());
        assert!(cover.run(&tc, &insertion).unwrap().is_empty());
        assert_eq!(
            scores_of(&cover.run(&tc, &spanning).unwrap()),
            Some(SOD1_SCORES)
        );
        assert_eq!(
            scores_of(&cover.run(&tc, &exact).unwrap()),
            Some(SOD1_SCORES)
        );

        for cover in [None, Some("0")] {
            let plain = plugin(dir.path(), cover);
            assert_eq!(
                scores_of(&plain.run(&tc, &inside).unwrap()),
                Some(SOD1_SCORES)
            );
            assert_eq!(
                scores_of(&plain.run(&tc, &spanning).unwrap()),
                Some(SOD1_SCORES)
            );
        }
    }

    /// `DosageSensitivity.pm:113-115, :152-154`: a structural variant is handed
    /// to the plugin (its `variant_feature_types` lists
    /// `StructuralVariationFeature`, which the dispatch gate of
    /// `OutputFactory.pm:2220` reads) and its feature's span is the one the
    /// cover test reads, so a deletion spanning the transcript prints the
    /// scores under `cover=1` and a duplication ending inside it does not.
    #[test]
    fn cover_applies_to_a_structural_variant_record() {
        let dir = tempfile::tempdir().unwrap();
        let cover = plugin(dir.path(), Some("1"));
        assert!(cover.annotates_structural_variants());
        let tc = row(Some("SOD1"), "ENSG00000142168", 31659666, 31668931);

        let mut deletion = InputVariant::new(
            "21".into(),
            31650001,
            31680000,
            b"N".to_vec(),
            b"<DEL>".to_vec(),
        );
        deletion.is_structural = true;
        let mut variants = vec![deletion];
        variants[0].transcript_consequences.push(tc.clone());
        cover.run_batch(&mut variants).unwrap();
        assert_eq!(
            scores_of(&variants[0].transcript_consequences[0].plugin_data),
            Some(SOD1_SCORES)
        );

        let mut duplication = InputVariant::new(
            "21".into(),
            31650001,
            31660000,
            b"N".to_vec(),
            b"<DUP>".to_vec(),
        );
        duplication.is_structural = true;
        let mut variants = vec![duplication];
        variants[0].transcript_consequences.push(tc);
        cover.run_batch(&mut variants).unwrap();
        assert!(variants[0].transcript_consequences[0]
            .plugin_data
            .is_empty());
    }

    /// `cover` is on for the string `1` alone and off under every other value,
    /// as `$param_hash->{cover} eq '1' ? 1 : 0` (`DosageSensitivity.pm:101`)
    /// reads it: `2`, `01`, `1.0`, `yes`, `true` and the empty value print the
    /// scores of a variant inside the transcript as `cover=0` does.
    #[test]
    fn cover_is_on_for_the_string_1_alone() {
        let dir = tempfile::tempdir().unwrap();
        let tc = row(Some("SOD1"), "ENSG00000142168", 31659666, 31668931);
        let inside = snv(31660000);
        for value in ["0", "2", "01", "1.0", "yes", "true", ""] {
            let off = plugin(dir.path(), Some(value));
            assert!(!off.cover, "cover={value}");
            assert_eq!(
                scores_of(&off.run(&tc, &inside).unwrap()),
                Some(SOD1_SCORES),
                "cover={value}"
            );
        }
        let on = plugin(dir.path(), Some("1"));
        assert!(on.cover);
        assert!(on.run(&tc, &inside).unwrap().is_empty());
    }

    /// An unknown key and a positional token beside a `key=value` are refused by
    /// the shared parameter grammar.
    #[test]
    fn unknown_key_and_mixed_grammar_are_init_errors() {
        let dir = tempfile::tempdir().unwrap();
        let file = write_table(dir.path(), "scores.tsv.gz", true);
        let file_param = format!("file={}", file.display());
        let init = |params: &[String]| {
            DosageSensitivityPlugin::new()
                .init(params)
                .unwrap_err()
                .to_string()
        };
        assert_eq!(
            init(&[file_param.clone(), "covers=1".into()]),
            "plugin initialization failed: DosageSensitivity: unknown parameter 'covers'; \
             DosageSensitivity accepts file, cover"
        );
        assert_eq!(
            init(&[file_param, "1".into()]),
            format!(
                "plugin initialization failed: DosageSensitivity: parameter '1' is not key=value \
                 while 'file={}' is; give every parameter as key=value, or every one by position",
                file.display()
            )
        );
    }

    /// `DosageSensitivity.pm:103-106, :180`: the row carries both keys, named
    /// `pHaplo` and `pTriplo`, with the table's values as written; the header
    /// describes each.
    #[test]
    fn both_scores_land_as_extra_keys() {
        let dir = tempfile::tempdir().unwrap();
        let plugin = plugin(dir.path(), None);
        assert_eq!(plugin.name(), "DosageSensitivity");
        assert_eq!(
            plugin.header_info(),
            [
                (
                    "pHaplo".to_string(),
                    "Probability of haploinsufficiency (deletion intolerance) of the affected gene"
                        .to_string()
                ),
                (
                    "pTriplo".to_string(),
                    "Probability of triplosensitivity (duplication intolerance) of the affected gene"
                        .to_string()
                ),
            ]
        );
        let mut variants = vec![snv(31660000)];
        variants[0]
            .transcript_consequences
            .push(row(Some("CACNA1C"), "ENSG00000151067", 100, 200));
        plugin.run_batch(&mut variants).unwrap();
        let extra: Vec<(&str, &str)> = variants[0].transcript_consequences[0]
            .plugin_data
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
        assert_eq!(extra, [("pHaplo", "0.99898184581082"), ("pTriplo", "1")]);
        assert!(variants[0].plugin_data.is_empty());
    }
}
