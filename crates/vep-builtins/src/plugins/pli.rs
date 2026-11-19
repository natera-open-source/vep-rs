// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! pLI plugin: the probability that a gene is loss-of-function intolerant, as
//! Ensembl VEP's pLI plugin (pLI.pm) writes it.
//!
//! **Usage:** `--plugin pLI[,<values file>[,<mode>]][,sigfig]`
//!
//! The first parameter is the values file; without one, or with one Perl reads
//! as false (`0` or the empty string), the plugin reads `pLI_values.txt` from
//! the first plugins directory that holds one, as the Ensembl VEP plugin reads
//! the copy bundled beside its module (pLI.pm:115-119).
//! The second is the mode, `gene` (the default) or `transcript` (`:109`):
//! `gene` files the scores under the file's gene symbols and writes
//! `pLI_gene_value` on each row whose symbol has one (`:161-172`, `:218-222`);
//! `transcript` files them under transcript stable ids and writes
//! `pLI_transcript_value` on each row whose stable id has one (`:174-185`,
//! `:224-228`). `sigfig` anywhere in the list prints a score with three
//! significant figures (`%.3g`) instead of two decimals (`%.2f`) (`:105-110`).
//! A positional parameter after the mode is not read, as `:108-109` read the
//! first two alone; it is logged. A mode other than the two is refused: Ensembl
//! VEP's plugin accepts it and then returns no hash from any row (`:218-229`),
//! which its output factory reports row by row (OutputFactory.pm:2228-2236).
//!
//! The file has two whitespace-separated columns under a header naming its key
//! column, `gene` or `transcript`: the last line containing either word in
//! lower case (`:125-131`); a header without the mode's column is refused
//! (`:162`, `:175`). Keys match case-insensitively (`:145`,
//! `:221`, `:227`). In transcript mode under `--refseq` or `--merged` the plugin
//! warns on standard error that stable ids may not match (`:187-189`).
//!
//! Per-transcript plugin: returns annotation from [`crate::BuiltinPlugin::run`] via
//! the default [`crate::BuiltinPlugin::run_batch`] dispatcher.

use std::io::BufReader;
use std::path::{Path, PathBuf};

use bstr::Finder;
use indexmap::IndexMap;
use tracing::debug;

use crate::params::{warn_unread_parameters, Grammar, PluginParams};
use crate::perl_format::{fixed, sig3};
use crate::plugins::gene_scores::{bundled_file, perl_space, MissingCell, ScoreTable};
use crate::tabix::BoundedLines;
use crate::traits::{warn_on_stderr, BuiltinPlugin, PluginError};
use vep_core::consequence::TranscriptConsequence;
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

/// The values file read from the plugins directory when none is given
/// (pLI.pm:118).
const BUNDLED_FILE: &str = "pLI_values.txt";

/// The flag selecting three significant figures (pLI.pm:105).
const SIGFIG_FLAG: &str = "sigfig";

/// What the scores are filed under (pLI.pm:91-98, 109).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum Mode {
    #[default]
    Gene,
    Transcript,
}

impl Mode {
    /// The header column the mode needs (pLI.pm:162, 175).
    fn column(self) -> &'static str {
        match self {
            Self::Gene => "gene",
            Self::Transcript => "transcript",
        }
    }

    /// The output field (pLI.pm:93, 96).
    fn field(self) -> &'static str {
        match self {
            Self::Gene => "pLI_gene_value",
            Self::Transcript => "pLI_transcript_value",
        }
    }

    /// The field's header description (pLI.pm:163, 176).
    fn description(self) -> &'static str {
        match self {
            Self::Gene => "pLI value by gene",
            Self::Transcript => "pLI value by transcript",
        }
    }

    /// The refusal of a header without the mode's column (pLI.pm:162, 175).
    fn missing_column(self) -> &'static str {
        match self {
            Self::Gene => "pLI: File does not have a gene column",
            Self::Transcript => "pLI: Could not find transcript in the headers",
        }
    }
}

/// `sprintf("%.2f")`, the default score format (pLI.pm:110).
fn two_decimals(value: f64) -> String {
    fixed(value, 2)
}

/// pLI gene constraint score plugin.
#[derive(Default)]
pub struct PliPlugin {
    /// Lower-cased key to the score as it prints.
    scores: ScoreTable,
    mode: Mode,
    /// The run's plugins directories, in search order.
    plugin_dirs: Vec<String>,
    /// Whether the run reads the RefSeq or the merged transcript set.
    refseq_or_merged: bool,
}

impl PliPlugin {
    /// Create a new uninitialized pLI plugin.
    pub fn new() -> Self {
        Self::default()
    }

    /// The header's columns: those of the last line containing `gene` or
    /// `transcript` in lower case, as pLI.pm:125-131's `next unless
    /// m/gene|transcript/` loop assigns every matching line over the one
    /// before. The file's existence is checked first (`:121`), and a file
    /// with no such line is refused (`:136`). The lines are matched and split
    /// as the bytes the module reads, so a byte that is not UTF-8 anywhere in
    /// the file changes nothing but the text of the column it sits in.
    fn header_columns(path: &Path) -> Result<Vec<String>, PluginError> {
        if !path.exists() {
            return Err(PluginError::Init(format!(
                "pLI: values file {} not found",
                path.display()
            )));
        }
        let file = std::fs::File::open(path).map_err(|e| {
            PluginError::Init(format!("pLI: failed to open {}: {e}", path.display()))
        })?;
        let mut header = None;
        // The two searchers are built once: a one-shot search builds its own
        // per call, which the file's thousands of lines would repeat.
        let (gene, transcript) = (Finder::new("gene"), Finder::new("transcript"));
        let mut lines = BoundedLines::new(BufReader::new(file), path);
        while let Some(line) = lines.next_line() {
            let line = line?;
            if gene.find(line).is_some() || transcript.find(line).is_some() {
                header = Some(line.to_vec());
            }
        }
        let header = header.ok_or_else(|| {
            PluginError::Init(format!(
                "pLI: Could not read headers from {}",
                path.display()
            ))
        })?;
        Ok(header
            .split(|byte| perl_space(*byte))
            .filter(|column| !column.is_empty())
            .map(|column| String::from_utf8_lossy(column).into_owned())
            .collect())
    }

    /// The row's score: the one filed under its gene symbol (pLI.pm:219-221) or
    /// its transcript stable id (`:225-227`).
    fn score_of(&self, consequence: &TranscriptConsequence) -> Option<&str> {
        match self.mode {
            Mode::Gene => self.scores.get(consequence.gene_symbol()?),
            Mode::Transcript => self.scores.get(consequence.transcript_id()),
        }
    }
}

impl BuiltinPlugin for PliPlugin {
    fn name(&self) -> &str {
        "pLI"
    }

    fn header_info(&self) -> Vec<(String, String)> {
        vec![(
            self.mode.field().to_string(),
            self.mode.description().to_string(),
        )]
    }

    fn set_plugin_dirs(&mut self, dirs: &[String]) {
        self.plugin_dirs = dirs.to_vec();
    }

    fn set_transcript_set(&mut self, refseq: bool, merged: bool) {
        self.refseq_or_merged = refseq || merged;
    }

    /// The values file, the mode and the `sigfig` flag (pLI.pm:104-110), read
    /// by position with the flag anywhere in the list; a later positional
    /// parameter is not read, as `:108-109` read the first two alone.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let parsed = PluginParams::parse_with(
            self.name(),
            Grammar::PositionalOrKeyValue {
                flags: &[SIGFIG_FLAG],
            },
            &[],
            params,
        )?;
        let positional = parsed.positional();
        // `if(!$file)` (pLI.pm:115): `0` and the empty string are no file.
        let path = match positional.first().filter(|file| perl_true(file)) {
            Some(file) => PathBuf::from(file),
            None => bundled_file(self.name(), "values", &self.plugin_dirs, BUNDLED_FILE)?,
        };
        self.mode = match positional.get(1).map(String::as_str) {
            None | Some("gene") => Mode::Gene,
            Some("transcript") => Mode::Transcript,
            Some(other) => {
                return Err(PluginError::Init(format!(
                    "pLI: mode '{other}' is not gene or transcript"
                )))
            }
        };
        let format: fn(f64) -> String = if parsed.flag(SIGFIG_FLAG) {
            sig3
        } else {
            two_decimals
        };

        let columns = Self::header_columns(&path)?;
        if !columns.iter().any(|column| column == self.mode.column()) {
            return Err(PluginError::Init(self.mode.missing_column().to_string()));
        }
        if positional.len() > 2 {
            warn_unread_parameters(
                self.name(),
                &positional[2..],
                "a values file, a mode (gene or transcript) and the flag sigfig",
            );
        }
        if self.mode == Mode::Transcript && self.refseq_or_merged {
            warn_on_stderr(
                "pLI transcript mode may fail with --refseq/--merged; consider gene mode",
            );
        }
        self.scores = ScoreTable::load(
            self.name(),
            "values",
            &path,
            "pLI",
            MissingCell::NoScore,
            format,
        )?;
        debug!(
            keys = self.scores.len(),
            mode = ?self.mode,
            path = %path.display(),
            "pLI plugin loaded scores"
        );
        Ok(())
    }

    fn run(
        &self,
        consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let mut result = IndexMap::new();
        if let Some(score) = self.score_of(consequence) {
            result.insert(self.mode.field().to_string(), score.to_string());
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use vep_core::transcript::Transcript;

    fn test_variant() -> InputVariant {
        InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec())
    }

    fn gene_row(symbol: &str) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                gene_symbol: Some(symbol.into()),
                ..Transcript::default()
            })),
            ..Default::default()
        }
    }

    fn transcript_row(stable_id: &str) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript: Some(Arc::new(Transcript {
                stable_id: stable_id.into(),
                ..Transcript::default()
            })),
            ..Default::default()
        }
    }

    fn params(tokens: &[&str]) -> Vec<String> {
        tokens.iter().map(|t| t.to_string()).collect()
    }

    /// A plugin loaded from `text` written to `pLI_values.txt` in a fresh
    /// directory, with `extra` after the file parameter; the directory is kept
    /// alive with it.
    fn loaded(text: &[u8], extra: &[&str]) -> (PliPlugin, tempfile::TempDir) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pLI_values.txt");
        std::fs::write(&path, text).unwrap();
        let mut plugin = PliPlugin::new();
        let mut params = vec![path.display().to_string()];
        params.extend(extra.iter().map(|t| t.to_string()));
        plugin.init(&params).unwrap();
        (plugin, dir)
    }

    fn value(plugin: &PliPlugin, row: &TranscriptConsequence) -> Option<String> {
        plugin
            .run(row, &test_variant())
            .unwrap()
            .get(plugin.mode.field())
            .cloned()
    }

    fn init_error(text: &[u8], extra: &[&str]) -> String {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pLI_values.txt");
        std::fs::write(&path, text).unwrap();
        let mut params = vec![path.display().to_string()];
        params.extend(extra.iter().map(|t| t.to_string()));
        PliPlugin::new().init(&params).unwrap_err().to_string()
    }

    const GENE_FILE: &[u8] = b"gene\tpLI\nAGRN\t0.17335234048116\nNOC2L\t1.33038194561114e-19\n";

    #[test]
    fn plugin_name_is_correct() {
        assert_eq!(PliPlugin::new().name(), "pLI");
    }

    /// Gene mode, the default, writes `pLI_gene_value` under `%.2f` (pLI.pm:109-110,
    /// 163): two values of the bundled file.
    #[test]
    fn gene_mode_formats_two_decimals() {
        let (plugin, _dir) = loaded(GENE_FILE, &[]);
        assert_eq!(
            plugin.header_info(),
            [(
                "pLI_gene_value".to_string(),
                "pLI value by gene".to_string()
            )]
        );
        assert_eq!(value(&plugin, &gene_row("AGRN")).as_deref(), Some("0.17"));
        assert_eq!(value(&plugin, &gene_row("NOC2L")).as_deref(), Some("0.00"));
        assert_eq!(value(&plugin, &gene_row("SOD1")), None);
        assert_eq!(value(&plugin, &TranscriptConsequence::default()), None);
        let (explicit, _dir) = loaded(GENE_FILE, &["gene"]);
        assert_eq!(value(&explicit, &gene_row("AGRN")).as_deref(), Some("0.17"));
    }

    /// `sigfig` anywhere in the list selects `%.3g` (pLI.pm:105-106, 110).
    #[test]
    fn sigfig_formats_three_significant_figures() {
        let (plugin, _dir) = loaded(GENE_FILE, &["sigfig"]);
        assert_eq!(value(&plugin, &gene_row("AGRN")).as_deref(), Some("0.173"));
        assert_eq!(
            value(&plugin, &gene_row("NOC2L")).as_deref(),
            Some("1.33e-19")
        );
        let (before_mode, _dir) = loaded(GENE_FILE, &["sigfig", "gene"]);
        assert_eq!(
            value(&before_mode, &gene_row("agrn")).as_deref(),
            Some("0.173")
        );
    }

    /// Transcript mode files the scores under `lc($transcript)` and looks a row
    /// up by `lc(stable_id)` (pLI.pm:182, 225-227), under `pLI_transcript_value`.
    #[test]
    fn transcript_mode_matches_stable_id_case_insensitively() {
        let (plugin, _dir) = loaded(
            b"transcript\tpLI\nENST00000379370\t0.17335234048116\nenst00000327044\t1.33038194561114e-19\n",
            &["transcript"],
        );
        assert_eq!(
            plugin.header_info(),
            [(
                "pLI_transcript_value".to_string(),
                "pLI value by transcript".to_string()
            )]
        );
        assert_eq!(
            value(&plugin, &transcript_row("ENST00000379370")).as_deref(),
            Some("0.17")
        );
        assert_eq!(
            value(&plugin, &transcript_row("enst00000379370")).as_deref(),
            Some("0.17")
        );
        assert_eq!(
            value(&plugin, &transcript_row("ENST00000327044")).as_deref(),
            Some("0.00")
        );
        assert_eq!(value(&plugin, &transcript_row("ENST00000000001")), None);
        let (sigfig, _dir) = loaded(
            b"transcript\tpLI\nENST00000379370\t0.17335234048116\n",
            &["transcript", "sigfig"],
        );
        assert_eq!(
            value(&sigfig, &transcript_row("ENST00000379370")).as_deref(),
            Some("0.173")
        );
    }

    /// Each mode needs its column in the header (pLI.pm:162, 175), and a file
    /// whose no line names either column has no header (`:136`).
    #[test]
    fn transcript_mode_without_transcript_column_is_an_error() {
        assert_eq!(
            init_error(GENE_FILE, &["transcript"]),
            "plugin initialization failed: pLI: Could not find transcript in the headers"
        );
        assert_eq!(
            init_error(b"transcript\tpLI\nENST00000379370\t0.5\n", &[]),
            "plugin initialization failed: pLI: File does not have a gene column"
        );
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pLI_values.txt");
        std::fs::write(&path, b"symbol\tpLI\nAGRN\t0.5\n").unwrap();
        assert_eq!(
            PliPlugin::new()
                .init(&[path.display().to_string()])
                .unwrap_err()
                .to_string(),
            format!(
                "plugin initialization failed: pLI: Could not read headers from {}",
                path.display()
            )
        );
    }

    /// The header is the last line containing `gene` or `transcript`, as
    /// pLI.pm:125-131's loop keeps the last line its `m/gene|transcript/`
    /// matches: a comment line before the table is read over, so the file
    /// loads, and a comment line after the table naming `transcripts` is the
    /// header, so the file is refused for its mode's column. The match is on
    /// the lower-case text alone (`:126`), so a `Gene` header is no header.
    #[test]
    fn the_header_is_the_last_line_naming_gene_or_transcript() {
        let (plugin, _dir) = loaded(
            b"# gene-level pLI scores\ngene\tpLI\nSOD1\t0.442188723698817\n",
            &[],
        );
        assert_eq!(value(&plugin, &gene_row("SOD1")).as_deref(), Some("0.44"));
        assert_eq!(
            init_error(
                b"gene\tpLI\nSOD1\t0.442188723698817\n# transcripts not listed\n",
                &[]
            ),
            "plugin initialization failed: pLI: File does not have a gene column"
        );
        assert_eq!(
            init_error(
                b"transcript\tpLI\nENST00000270142\t0.442188723698817\n# genes not listed\n",
                &["transcript"]
            ),
            "plugin initialization failed: pLI: Could not find transcript in the headers"
        );
        let (plugin, _dir) = loaded(
            b"# transcripts not listed\ngene\tpLI\nSOD1\t0.442188723698817\n# gene scores end\n",
            &[],
        );
        assert_eq!(value(&plugin, &gene_row("SOD1")).as_deref(), Some("0.44"));
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pLI_values.txt");
        std::fs::write(&path, b"Gene\tpLI\nSOD1\t0.442188723698817\n").unwrap();
        assert_eq!(
            PliPlugin::new()
                .init(&[path.display().to_string()])
                .unwrap_err()
                .to_string(),
            format!(
                "plugin initialization failed: pLI: Could not read headers from {}",
                path.display()
            )
        );
    }

    /// A byte that is not UTF-8 (Latin-1 `é`, 0xE9) is a byte of the line it
    /// sits in, as `<IN>` and `<$fh>` (pLI.pm:125, 141) read it: in a comment
    /// line, in a third cell `my ($gene, $score) = split` (`:143`) never reads
    /// and in a header column the mode does not need, the file loads and SOD1
    /// scores 0.44; in the key cell the key is another gene, so SOD1 has no
    /// score; after the score it is the text `sprintf` numifies past (`:145`),
    /// so SOD1 scores 0.44. The header is still refused when the byte sits in
    /// the mode's column name.
    #[test]
    fn a_byte_that_is_not_utf8_is_a_byte_of_its_line() {
        let shapes: [(&[u8], Option<&str>); 5] = [
            (
                b"# caf\xE9 scores\ngene\tpLI\nSOD1\t0.442188723698817\n",
                Some("0.44"),
            ),
            (
                b"gene\tpLI\nSOD1\t0.442188723698817\tcaf\xE9\n",
                Some("0.44"),
            ),
            (
                b"gene\tpLI\tcaf\xE9\nSOD1\t0.442188723698817\n",
                Some("0.44"),
            ),
            (b"gene\tpLI\nSOD\xE91\t0.442188723698817\n", None),
            (b"gene\tpLI\nSOD1\t0.4421\xE9\n", Some("0.44")),
        ];
        for (text, score) in shapes {
            let (plugin, _dir) = loaded(text, &[]);
            assert_eq!(
                value(&plugin, &gene_row("SOD1")).as_deref(),
                score,
                "{}",
                String::from_utf8_lossy(text)
            );
        }
        assert_eq!(
            init_error(b"gene\xE9\tpLI\nSOD1\t0.442188723698817\n", &[]),
            "plugin initialization failed: pLI: File does not have a gene column"
        );
    }

    /// The columns are whitespace-separated (`split`, pLI.pm:129, 143), so the
    /// space-separated extract the plugin's documentation produces loads as the
    /// tab-separated bundled file does, and a later line spelling a symbol in
    /// another case replaces the earlier (`:145`).
    #[test]
    fn space_separated_file_loads() {
        let (plugin, _dir) = loaded(
            b"gene pLI\nAGRN 0.17335234048116\nC2ORF15 0.013119592808789\nC2orf15 0.543252452898348\n",
            &[],
        );
        assert_eq!(value(&plugin, &gene_row("AGRN")).as_deref(), Some("0.17"));
        assert_eq!(
            value(&plugin, &gene_row("C2orf15")).as_deref(),
            Some("0.54")
        );
        assert_eq!(
            value(&plugin, &gene_row("C2ORF15")).as_deref(),
            Some("0.54")
        );
    }

    /// Without a file parameter the plugin reads `pLI_values.txt` from the
    /// plugins directory (pLI.pm:115-119), the mode and the flag still read
    /// (`:104-109`); a directory without it is refused naming the path (`:121`).
    #[test]
    fn default_values_file_resolves_from_plugins_dir() {
        let dir = tempfile::tempdir().unwrap();
        let dirs = [dir.path().display().to_string()];
        let mut plugin = PliPlugin::new();
        plugin.set_plugin_dirs(&dirs);
        assert_eq!(
            plugin.init(&[]).unwrap_err().to_string(),
            format!(
                "plugin initialization failed: pLI: values file {} not found",
                dir.path().join("pLI_values.txt").display()
            )
        );
        std::fs::write(dir.path().join("pLI_values.txt"), GENE_FILE).unwrap();
        let mut plugin = PliPlugin::new();
        plugin.set_plugin_dirs(&dirs);
        plugin.init(&params(&["sigfig"])).unwrap();
        assert_eq!(value(&plugin, &gene_row("AGRN")).as_deref(), Some("0.173"));
        let mut plugin = PliPlugin::new();
        assert_eq!(
            plugin.init(&[]).unwrap_err().to_string(),
            "plugin initialization failed: pLI: values file pLI_values.txt not found"
        );
    }

    /// The file parameter is read as Perl's `if(!$file)` reads it (pLI.pm:115):
    /// `0` and the empty string are no file, so the bundled file is read with
    /// the mode and the flag after them, while `0.0` is a path and refused.
    #[test]
    fn a_false_file_parameter_reads_the_bundled_file() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("pLI_values.txt"), GENE_FILE).unwrap();
        let dirs = [dir.path().display().to_string()];
        for (tokens, expected) in [
            (&["0"][..], "0.17"),
            (&[""][..], "0.17"),
            (&["0", "gene"][..], "0.17"),
            (&["0", "sigfig"][..], "0.173"),
            (&["", "sigfig"][..], "0.173"),
        ] {
            let mut plugin = PliPlugin::new();
            plugin.set_plugin_dirs(&dirs);
            plugin.init(&params(tokens)).unwrap();
            assert_eq!(
                value(&plugin, &gene_row("AGRN")).as_deref(),
                Some(expected),
                "{tokens:?}"
            );
        }
        let mut plugin = PliPlugin::new();
        plugin.set_plugin_dirs(&dirs);
        assert_eq!(
            plugin.init(&params(&["0.0"])).unwrap_err().to_string(),
            "plugin initialization failed: pLI: values file 0.0 not found"
        );
    }

    /// A mode other than `gene` or `transcript`, the empty mode of `<file>,,sigfig`
    /// among them, and a `key=value` parameter are refused (pLI.pm:104-109 reads
    /// by position; an empty second token is a mode of `''` there, under which
    /// no row gets the field).
    #[test]
    fn unknown_mode_and_keyed_parameters_are_refused() {
        assert_eq!(
            init_error(GENE_FILE, &["genes"]),
            "plugin initialization failed: pLI: mode 'genes' is not gene or transcript"
        );
        assert_eq!(
            init_error(GENE_FILE, &["", "sigfig"]),
            "plugin initialization failed: pLI: mode '' is not gene or transcript"
        );
        assert_eq!(
            PliPlugin::new()
                .init(&params(&["file=pLI_values.txt"]))
                .unwrap_err()
                .to_string(),
            "plugin initialization failed: pLI: unknown parameter 'file'; pLI reads its \
             parameters by position"
        );
    }

    /// Ensembl VEP's plugin reads the first two positional parameters and the
    /// flag alone (pLI.pm:105-109), so a later token is not read: `sigfigs` is
    /// no flag and leaves the two-decimal format, and `extra` before `sigfig`
    /// leaves the flag in force.
    #[test]
    fn surplus_positional_parameters_are_ignored() {
        let (misspelt, _dir) = loaded(GENE_FILE, &["gene", "sigfigs"]);
        assert_eq!(value(&misspelt, &gene_row("AGRN")).as_deref(), Some("0.17"));
        let (surplus, _dir) = loaded(GENE_FILE, &["gene", "extra", "sigfig"]);
        assert_eq!(value(&surplus, &gene_row("AGRN")).as_deref(), Some("0.173"));
        let (transcript, _dir) = loaded(
            b"transcript\tpLI\nENST00000379370\t0.17335234048116\n",
            &["transcript", "extra"],
        );
        assert_eq!(
            value(&transcript, &transcript_row("ENST00000379370")).as_deref(),
            Some("0.17")
        );
    }

    /// An empty values file has no header line and is refused as the plugin
    /// opens, naming the file (pLI.pm:136).
    #[test]
    fn empty_scores_file_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pLI_values.txt");
        std::fs::write(&path, b"").unwrap();
        let error = PliPlugin::new()
            .init(&[path.display().to_string()])
            .unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "plugin initialization failed: pLI: Could not read headers from {}",
                path.display()
            )
        );
    }

    /// A scores file with CRLF line ends loads as one with LF line ends: the
    /// line end is stripped whole, so the last header cell is `pLI`, not
    /// `pLI\r`, and the score of the last column reads. A header split by a
    /// vertical tab reads too, as pLI.pm:129's `split` breaks on one.
    #[test]
    fn crlf_scores_file_loads() {
        let (plugin, _dir) = loaded(b"gene\tpLI\r\nSOD1\t0.95\r\n", &[]);
        assert_eq!(value(&plugin, &gene_row("SOD1")).as_deref(), Some("0.95"));
        let (plugin, _dir) = loaded(b"gene\x0BpLI\nSOD1\x0B0.95\n", &[]);
        assert_eq!(value(&plugin, &gene_row("SOD1")).as_deref(), Some("0.95"));
    }

    /// A scores file whose third line is one byte over the cap is refused as
    /// the plugin opens, naming the file and the line, and the scores read
    /// before it are not kept.
    #[test]
    fn scores_line_past_the_cap_is_refused() {
        let cap = crate::tabix::MAX_LINE_BYTES;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pLI_values.txt");
        let mut text = b"gene\tpLI\nSOD1\t0.95\n".to_vec();
        text.extend(std::iter::repeat_n(b'x', cap + 1));
        text.push(b'\n');
        std::fs::write(&path, text).unwrap();
        let mut plugin = PliPlugin::new();
        let error = plugin.init(&[path.display().to_string()]).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "plugin initialization failed: {}: line 3 is longer than the {cap}-byte limit \
                 for one line",
                path.display()
            )
        );
        assert!(plugin
            .run(&gene_row("SOD1"), &test_variant())
            .unwrap()
            .is_empty());
    }

    /// `run_batch` writes the mode's field as `run` returns it, on transcript
    /// rows alone.
    #[test]
    fn run_batch_writes_the_field_on_transcript_rows() {
        let (plugin, _dir) = loaded(GENE_FILE, &[]);
        let mut variants = vec![test_variant()];
        variants[0].transcript_consequences.extend([
            gene_row("AGRN"),
            gene_row("SOD1"),
            TranscriptConsequence {
                feature_type: vep_core::consequence::FeatureType::RegulatoryFeature,
                ..gene_row("AGRN")
            },
        ]);
        plugin.run_batch(&mut variants).unwrap();
        let written: Vec<Option<&String>> = variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.get("pLI_gene_value"))
            .collect();
        assert_eq!(written, [Some(&"0.17".to_string()), None, None]);
    }

    /// Ensembl VEP 116.2 never hands pLI a structural variant: the module keeps
    /// the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2220). The transcript rows of a
    /// `<DEL>` over a scored gene carry no field, while the sequence variant's
    /// row of the same gene does.
    #[test]
    fn a_structural_variant_row_gets_no_field() {
        let (plugin, _dir) = loaded(GENE_FILE, &[]);
        assert!(!plugin.annotates_structural_variants());
        let mut deletion =
            InputVariant::new("1".into(), 100, 5000, b"N".to_vec(), b"<DEL>".to_vec());
        deletion.is_structural = true;
        deletion.transcript_consequences.push(gene_row("AGRN"));
        let mut sequence = test_variant();
        sequence.transcript_consequences.push(gene_row("AGRN"));
        let mut variants = vec![deletion, sequence];
        plugin.run_batch(&mut variants).unwrap();
        assert!(variants[0].transcript_consequences[0]
            .plugin_data
            .is_empty());
        assert_eq!(
            variants[1].transcript_consequences[0]
                .plugin_data
                .get("pLI_gene_value")
                .map(String::as_str),
            Some("0.17")
        );
    }
}
