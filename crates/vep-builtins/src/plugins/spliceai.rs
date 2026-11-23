// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! SpliceAI plugin: splice-site prediction scores read from the tabix-indexed
//! SpliceAI score files (one for single-nucleotide variants, one for indels),
//! matched by position, allele and gene symbol.
//!
//! # Parameters
//!
//! Every parameter is `key=value` (`SpliceAI.pm:160-163` reads `params_to_hash`
//! alone): `snv=<path>` and `indel=<path>` name the score files, of which at
//! least one is given (Ensembl VEP's plugin dies without both,
//! `:162-163`);
//! `cutoff=<value>` adds the
//! `SpliceAI_cutoff` flag (`:165-171`, `:222-224`, `:295-303`), `PASS` when the
//! largest delta score reaches the cutoff and `FAIL` otherwise, and drops no
//! record: the cutoff is the value's number as `$cutoff < 0 || $cutoff > 1`
//! (`:167`) numifies it (`0.5x` is 0.5, `abc` is 0), refused outside `0..=1`
//! (`:167-168`), and the empty value sets no cutoff, as it fails
//! `if($self->{cutoff})` (`:222`), while the literal `0` sets a cutoff of 0,
//! the documented divergence from Ensembl VEP; `split_output=<anything>` prints
//! the nine fields (`:177-179`).
//!
//! # Output shape
//!
//! Ensembl VEP's plugin shapes its fields by the run's output format
//! (`:173-175`, `:208-220`, `:275-292`, `:316`), read here through
//! [`BuiltinPlugin::set_output_format`]:
//!
//! - the default, tab and Parquet formats print one field, `SpliceAI_pred`,
//!   the record's value after its allele: `SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|
//!   DP_AL|DP_DG|DP_DL`;
//! - the VCF format, and any format under `split_output`, prints the nine
//!   `SpliceAI_pred_<part>` fields;
//! - the JSON format nests the nine parts and the flag under `spliceai`, built
//!   by the JSON writer from the nine fields written here.
//!
//! A record's entries reach the transcript rows whose gene symbol equals an
//! entry's `SYMBOL` exactly (`:313-317`); a row without a symbol takes nothing.
//! A row carries every entry of the record that names its gene, in the
//! record's order: the one field holds them as the record spells them from the
//! first of them on, each later one with its allele part, joined by a comma;
//! each of the nine fields holds its values over them, joined by a comma; and
//! the flag is `PASS` where any of them reaches the cutoff. When several
//! records naming one gene match the variant the last one wins (`:305`). Both
//! files are read for every variant (`BaseVepTabixPlugin.pm:591`).
//! The header descriptions carry the tool version read from the first file's
//! `##INFO` line (`:193-206`).
//!
//! `SpliceAI.pm` is the plugin's Ensembl VEP counterpart in
//! Ensembl/VEP_plugins release/116.

use std::borrow::Cow;
use std::io::BufReader;
use std::path::Path;

use indexmap::IndexMap;
use noodles::bgzf;
use vep_core::assembly::Assembly;
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

use crate::annotation_store::{open_best_store_with, AnnotationStore};
use crate::params::{perl_numeric, Grammar, PluginParams};
use crate::tabix::{
    allele_match, lossy_text, table_contig, BatchQueryResult, BoundedLines, TabixAnnotator,
    TabixConfig, TabixRecord,
};
use crate::traits::{BuiltinPlugin, OutputFormat, PluginError, PrefetchData};

const NAME: &str = "SpliceAI";
const PARAMETERS: &[&str] = &["snv", "indel", "cutoff", "split_output"];

/// The one field of the default and tab formats (`SpliceAI.pm:219, 291`).
const JOINED_FIELD: &str = "SpliceAI_pred";
/// The flag field `cutoff=` adds (`SpliceAI.pm:223, 302`).
const CUTOFF_FIELD: &str = "SpliceAI_cutoff";
/// The nine fields of the VCF format and of `split_output`, in the order of the
/// parts of a record's value after its allele (`SpliceAI.pm:279-287`).
const SPLIT_FIELDS: [&str; 9] = [
    "SpliceAI_pred_SYMBOL",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL",
    "SpliceAI_pred_DP_DG",
    "SpliceAI_pred_DP_DL",
];
/// The parts of a value after its allele that are delta scores: after the
/// symbol, before the four delta positions.
const DELTA_SCORE_PARTS: std::ops::Range<usize> = 1..5;

/// VCF POS column index (CHROM=0, POS=1).
const VCF_POS_COL: usize = 1;
/// VCF REF column index (CHROM=0, POS=1, ID=2, REF=3).
const VCF_REF_COL: usize = 3;
/// VCF ALT column index (CHROM=0, POS=1, ID=2, REF=3, ALT=4).
const VCF_ALT_COL: usize = 4;
/// VCF INFO column index.
const VCF_INFO_COL: usize = 7;

/// One entry of a record's INFO column as `SpliceAiPlugin::score_values` reads
/// it: the entry as the record spells it, its allele part included, the gene
/// symbol it names, and the value after the allele, `SYMBOL|DS_AG|...|DP_DL`.
struct ScoreEntry<'a> {
    spelled: &'a str,
    symbol: &'a str,
    scores: Cow<'a, str>,
}

/// The entries of one matching record that name one gene symbol, in the
/// record's order, and the cutoff flag when a cutoff is set: `PASS` where any
/// of the entries reaches the cutoff by its own largest delta score, `FAIL`
/// where none does.
struct Hit<'a> {
    symbol: &'a str,
    entries: Vec<ScoreEntry<'a>>,
    flag: Option<&'static str>,
}

/// SpliceAI splice prediction plugin.
pub struct SpliceAiPlugin {
    /// The score files, the SNV file first, as Ensembl VEP's plugin adds them.
    stores: Vec<Box<dyn AnnotationStore>>,
    /// The cutoff the flag compares the largest delta score with: the Perl
    /// number of the `cutoff=` text (`$cutoff < 0`, `SpliceAI.pm:167`); `None`
    /// without the parameter or with the empty value, which fails
    /// `if($self->{cutoff})` (`:222, 295, 335`).
    cutoff: Option<f64>,
    /// Whether `split_output` was given.
    split_output: bool,
    output_format: OutputFormat,
    /// The tool version the first file's `##INFO` line names, such as `v1.3`.
    version: Option<String>,
    /// The run's assembly, which names a mitochondrial record's contig to the
    /// table ([`table_contig`]).
    assembly: Option<Assembly>,
}

impl Default for SpliceAiPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl SpliceAiPlugin {
    /// Create an uninitialized SpliceAI plugin.
    ///
    /// Call [`init`](BuiltinPlugin::init) with the score file paths before use.
    pub fn new() -> Self {
        Self {
            stores: Vec::new(),
            cutoff: None,
            split_output: false,
            output_format: OutputFormat::default(),
            version: None,
            assembly: None,
        }
    }

    /// Build a [`TabixConfig`] for a VCF-format SpliceAI file.
    ///
    /// VCF columns: CHROM(0), POS(1), ID(2), REF(3), ALT(4), QUAL(5),
    /// FILTER(6), INFO(7).
    fn vcf_config(path: &Path) -> TabixConfig {
        TabixConfig {
            file_path: path.to_path_buf(),
            chr_col: 0,
            start_col: 1,
            end_col: None,
            ref_col: Some(VCF_REF_COL),
            alt_col: Some(VCF_ALT_COL),
            zero_based: false,
        }
    }

    /// Whether the run prints the nine fields: the VCF format, the JSON format
    /// (whose writer nests them) and `split_output` (`SpliceAI.pm:208, 275`).
    fn splits(&self) -> bool {
        self.split_output || matches!(self.output_format, OutputFormat::Vcf | OutputFormat::Json)
    }

    /// The `SpliceAI (<version>) ` opening of every description; without a
    /// version token in the file, `SpliceAI ` alone.
    fn description_prefix(&self) -> String {
        match &self.version {
            Some(version) => format!("{NAME} ({version}) "),
            None => format!("{NAME} "),
        }
    }

    /// The region read for a variant: one base before its start, so an indel's
    /// VCF record, anchored on the base before the change, is found, to its end.
    fn query_bounds(variant: &InputVariant) -> (u64, u64) {
        let query_start = variant.start.saturating_sub(1).max(1);
        let query_end = variant.end.max(variant.start);
        (query_start, query_end)
    }

    /// The `SpliceAI=` entries of a record's INFO column, each the value after
    /// its allele as `parse_data` keeps it (`SpliceAI.pm:326-349`): the entry
    /// split on `|` into the allele, the symbol and the rest, the value the
    /// symbol and the rest joined by `|`, so a truncated entry is kept with the
    /// parts it has (`G|SYM|0.5|0.6` is `SYM|0.5|0.6`), an entry of two parts
    /// is its symbol and a trailing `|` (`:331` joins an undefined rest), and
    /// the allele part is never read, the match being on the REF and ALT
    /// columns (`:272`). An entry of one part has no symbol, which no row's
    /// gene symbol equals (`:314`), so it is left out. One INFO may carry
    /// several entries joined by `,`, and each is read on its own.
    fn score_values(record: &TabixRecord) -> impl Iterator<Item = ScoreEntry<'_>> {
        let info = record
            .get("INFO")
            .or_else(|| record.columns.get(VCF_INFO_COL).map(String::as_str))
            .unwrap_or("");
        info.split(';')
            .filter_map(|field| field.strip_prefix("SpliceAI="))
            .flat_map(|value| value.split(','))
            .filter_map(|entry| {
                let (_allele, rest) = entry.split_once('|')?;
                let (symbol, scores) = match rest.split_once('|') {
                    Some((symbol, _)) => (symbol, Cow::Borrowed(rest)),
                    None => (rest, Cow::Owned(format!("{rest}|"))),
                };
                Some(ScoreEntry {
                    spelled: entry,
                    symbol,
                    scores,
                })
            })
    }

    /// The largest delta score of a value (`SpliceAI.pm:336-338`); a part that
    /// is not a number counts as zero, as Perl's `max` reads it.
    fn max_delta_score(scores: &str) -> f64 {
        scores
            .split('|')
            .skip(DELTA_SCORE_PARTS.start)
            .take(DELTA_SCORE_PARTS.len())
            .filter_map(|part| part.parse::<f64>().ok())
            .fold(0.0, f64::max)
    }

    /// Whether a record is the variant's (`SpliceAI.pm:259-272`): the record's
    /// position and alleles equal the variant's in the file's anchored form,
    /// which Ensembl VEP builds through `to_VCF_record` for an allele string
    /// carrying `-`, or equal the variant's start and alleles as written, the
    /// comparison Ensembl VEP makes for every other allele string, so the
    /// insertion or deletion allele of a multi-allelic record whose alleles
    /// share no first base (`C/G/CA`, kept unchopped) takes its record.
    fn is_variant_record(record: &TabixRecord, variant: &InputVariant) -> bool {
        allele_match::matches_allele_with_position(
            record,
            variant,
            VCF_POS_COL,
            VCF_REF_COL,
            VCF_ALT_COL,
        ) || Self::written_alleles_match(record, variant)
    }

    /// Whether the record's position, REF and ALT equal the variant's start and
    /// forward-strand alleles as written, case-insensitively. The position is
    /// the POS cell as `$start == $data_value->{start}` (`:272`) numifies it
    /// ([`allele_match::perl_position`]): `21418500.0` and `21418500x` are
    /// 21418500, `21418500.5` is no position.
    fn written_alleles_match(record: &TabixRecord, variant: &InputVariant) -> bool {
        if allele_match::perl_position(record, VCF_POS_COL) != Some(variant.start) {
            return false;
        }
        let (Some(file_ref), Some(file_alt)) = (
            record.columns.get(VCF_REF_COL),
            record.columns.get(VCF_ALT_COL),
        ) else {
            return false;
        };
        let (var_ref, var_alt) = variant.forward_strand_alleles();
        file_ref.as_bytes().eq_ignore_ascii_case(&var_ref)
            && file_alt.as_bytes().eq_ignore_ascii_case(&var_alt)
    }

    /// The entries of the matching records of `records`, one hit per gene
    /// symbol: every entry of a record that names the symbol, in the record's
    /// order, the last record that names a gene taking the place of an earlier
    /// one's entries (`SpliceAI.pm:248-307`, the hash of `:305` keeping one
    /// record a gene). A record matches when it is the variant's
    /// ([`Self::is_variant_record`]).
    fn hits<'a>(&self, variant: &InputVariant, records: &[&'a TabixRecord]) -> Vec<Hit<'a>> {
        let mut hits: Vec<Hit> = Vec::new();
        for record in records {
            if !Self::is_variant_record(record, variant) {
                continue;
            }
            let mut of_record: Vec<Hit> = Vec::new();
            for entry in Self::score_values(record) {
                match of_record.iter_mut().find(|hit| hit.symbol == entry.symbol) {
                    Some(hit) => hit.entries.push(entry),
                    None => of_record.push(Hit {
                        symbol: entry.symbol,
                        entries: vec![entry],
                        flag: None,
                    }),
                }
            }
            for mut hit in of_record {
                hit.flag = self.cutoff.map(|cutoff| {
                    let reaches =
                        |entry: &ScoreEntry| Self::max_delta_score(&entry.scores) >= cutoff;
                    if hit.entries.iter().any(reaches) {
                        "PASS"
                    } else {
                        "FAIL"
                    }
                });
                match hits.iter_mut().find(|earlier| earlier.symbol == hit.symbol) {
                    Some(earlier) => *earlier = hit,
                    None => hits.push(hit),
                }
            }
        }
        hits
    }

    /// A hit's fields in the run's shape, with the names of the fields whose
    /// empty value is the empty string. The nine fields of the split shape are
    /// an entry's `split /\|/` (`SpliceAI.pm:276-287`), which drops the
    /// trailing empty parts, read into nine variables: a part the entry lacks is
    /// the undefined value, written as the empty value, and an empty part before
    /// a filled one the empty string. Over several entries each of the nine
    /// fields is its values joined by a comma, in the entries' order. The one
    /// field of the other shape is a lone entry's value, and over several
    /// entries the first without its allele part and each later one as the
    /// record spells it, joined by a comma, the text `:328-332` keep of a
    /// record that names the gene alone.
    fn fields_of(&self, hit: &Hit<'_>) -> (IndexMap<String, String>, Vec<String>) {
        let mut fields = IndexMap::new();
        let mut empty_strings = Vec::new();
        if self.splits() {
            fields.reserve(SPLIT_FIELDS.len() + usize::from(hit.flag.is_some()));
            let parts_of: Vec<Vec<&str>> = hit
                .entries
                .iter()
                .map(|entry| {
                    let mut parts: Vec<&str> = entry.scores.split('|').collect();
                    while parts.last() == Some(&"") {
                        parts.pop();
                    }
                    parts
                })
                .collect();
            for (index, name) in SPLIT_FIELDS.iter().enumerate() {
                let values: Vec<&str> = parts_of
                    .iter()
                    .map(|parts| parts.get(index).copied().unwrap_or(""))
                    .collect();
                fields.insert((*name).to_string(), values.join(","));
                if let [parts] = parts_of.as_slice() {
                    if parts.get(index) == Some(&"") {
                        empty_strings.push((*name).to_string());
                    }
                }
            }
        } else {
            let value = match hit.entries.as_slice() {
                [entry] => entry.scores.to_string(),
                entries => entries
                    .iter()
                    .enumerate()
                    .map(|(place, entry)| match place {
                        0 => entry.spelled.split_once('|').map_or("", |(_, rest)| rest),
                        _ => entry.spelled,
                    })
                    .collect::<Vec<_>>()
                    .join(","),
            };
            fields.insert(JOINED_FIELD.to_string(), value);
        }
        if let Some(flag) = hit.flag {
            fields.insert(CUTOFF_FIELD.to_string(), flag.to_string());
        }
        (fields, empty_strings)
    }

    /// Writes the hit of each transcript row's gene symbol, every entry of the
    /// row's record that names it, on that row (`SpliceAI.pm:313-317`).
    fn annotate_variant(&self, variant: &mut InputVariant, records: &[&TabixRecord]) {
        let hits = self.hits(variant, records);
        if hits.is_empty() {
            return;
        }
        for consequence in &mut variant.transcript_consequences {
            if !self.annotates_row(consequence) {
                continue;
            }
            let Some(symbol) = consequence.gene_symbol() else {
                continue;
            };
            if let Some(hit) = hits.iter().find(|hit| hit.symbol == symbol) {
                let (fields, empty_strings) = self.fields_of(hit);
                consequence
                    .plugin_empty_strings
                    .retain(|key| !SPLIT_FIELDS.contains(&key.as_str()));
                consequence.plugin_empty_strings.extend(empty_strings);
                consequence.plugin_data.extend(fields);
            }
        }
    }

    /// The records of every store in the variant's region, the SNV file's first.
    fn query_stores(&self, variant: &InputVariant) -> Result<Vec<TabixRecord>, PluginError> {
        let (start, end) = Self::query_bounds(variant);
        let mut records = Vec::new();
        for store in &self.stores {
            records.extend(store.query(table_contig(variant, self.assembly), start, end)?);
        }
        Ok(records)
    }

    /// The tool version named by the last `##INFO` header line carrying a
    /// `SpliceAIv<digits and dots>` token (`SpliceAI.pm:196-205`), with the
    /// `SpliceAI` prefix removed: `v1.3` from `SpliceAIv1.3`. `None` when no header
    /// line carries one.
    fn read_version(path: &Path) -> Result<Option<String>, PluginError> {
        let file = std::fs::File::open(path).map_err(|e| {
            PluginError::Init(format!("{NAME}: failed to open {}: {e}", path.display()))
        })?;
        let reader = bgzf::io::Reader::new(BufReader::new(file));
        let mut version = None;
        for line in BoundedLines::new(reader, path) {
            let line = line?;
            if !line.starts_with(b"#") {
                break;
            }
            if !line.starts_with(b"##INFO") {
                continue;
            }
            let line = lossy_text(line);
            if let Some(token) = line.split("SpliceAIv").nth(1) {
                let digits: String = token
                    .chars()
                    .take_while(|c| c.is_ascii_digit() || *c == '.')
                    .collect();
                if !digits.is_empty() {
                    version = Some(format!("v{digits}"));
                }
            }
        }
        Ok(version)
    }
}

impl crate::traits::sealed::Sealed for SpliceAiPlugin {}

impl BuiltinPlugin for SpliceAiPlugin {
    fn name(&self) -> &str {
        NAME
    }

    /// The descriptions of `SpliceAI.pm:208-224`, the run's shape deciding
    /// between the one field and the nine.
    fn header_info(&self) -> Vec<(String, String)> {
        let prefix = self.description_prefix();
        let mut header: Vec<(String, String)> = if self.splits() {
            [
                ("SpliceAI_pred_SYMBOL", "gene symbol"),
                ("SpliceAI_pred_DS_AG", "delta score (0–1) for acceptor gain: max inc. in predicted acceptor probability within ±50 bp. Guide cutoffs: <0.01 none; ≥0.2 high-recall; ≥0.5 recommended; ≥0.8 high-precision."),
                ("SpliceAI_pred_DS_AL", "delta score (0–1) for acceptor loss: max dec. in predicted acceptor probability within ±50 bp. Guide cutoffs: <0.01 none; ≥0.2 high-recall; ≥0.5 recommended; ≥0.8 high-precision."),
                ("SpliceAI_pred_DS_DG", "delta score (0–1) for donor gain: max inc. in predicted donor probability within ±50 bp. Guide cutoffs: <0.01 none; ≥0.2 high-recall; ≥0.5 recommended; ≥0.8 high-precision."),
                ("SpliceAI_pred_DS_DL", "delta score (0–1) for donor loss: max dec. in predicted donor probability within ±50 bp. Guide cutoffs: <0.01 none; ≥0.2 high-recall; ≥0.5 recommended; ≥0.8 high-precision."),
                ("SpliceAI_pred_DP_AG", "delta position (bp) for acceptor gain: position (relative to variant) of max AG change (+ downstream, − upstream)."),
                ("SpliceAI_pred_DP_AL", "delta position (bp) for acceptor loss: position (relative to variant) of max AL change (+ downstream, − upstream)."),
                ("SpliceAI_pred_DP_DG", "delta position (bp) for donor gain: position (relative to variant) of max DG change (+ downstream, − upstream)."),
                ("SpliceAI_pred_DP_DL", "delta position (bp) for donor loss: position (relative to variant) of max DL change (+ downstream, − upstream)."),
            ]
            .iter()
            .map(|(field, text)| ((*field).to_string(), format!("{prefix}{text}")))
            .collect()
        } else {
            vec![(
                JOINED_FIELD.to_string(),
                format!(
                    "{prefix}predicted effect on splicing. These include delta scores (DS) and \
                     delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain \
                     (DG), and donor loss (DL). Format: SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|\
                     DP_AL|DP_DG|DP_DL"
                ),
            )]
        };
        if self.cutoff.is_some() {
            header.push((
                CUTOFF_FIELD.to_string(),
                "Flag if delta score pass the cutoff (PASS) or if it does not (FAIL)".to_string(),
            ));
        }
        header
    }

    fn set_output_format(&mut self, format: OutputFormat) {
        self.output_format = format;
    }

    fn set_assembly(&mut self, assembly: Option<Assembly>) {
        self.assembly = assembly;
    }

    fn set_reference_fasta(&mut self, fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>) {
        for store in &mut self.stores {
            store.set_reference_fasta(fasta.clone());
        }
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse_with(
            NAME,
            Grammar::KeyValueOnly { file_key: "snv" },
            PARAMETERS,
            params,
        )?;
        let paths: Vec<&str> = ["snv", "indel"]
            .iter()
            .filter_map(|key| params.get(key))
            .filter(|path| !path.is_empty())
            .collect();
        if paths.is_empty() {
            return Err(PluginError::Init(format!(
                "{NAME}: no score file given; use snv=<path> and indel=<path>"
            )));
        }
        if let Some(value) = params.get("cutoff") {
            // `$cutoff < 0 || $cutoff > 1` (SpliceAI.pm:167) numifies the text
            // and dies outside `0..=1`, both comparisons false for NaN, which
            // is kept as a cutoff every score fails; the kept text then meets
            // `if($self->{cutoff})` (:222, :295, :335), which the empty value
            // fails. The literal `0` sets a cutoff of 0 here, the documented
            // divergence from Ensembl VEP, whose truth test drops it.
            let number = perl_numeric(NAME, "cutoff", value);
            if !number.is_nan() && !(0.0..=1.0).contains(&number) {
                return Err(PluginError::Init(format!(
                    "{NAME}: cutoff must be between 0 and 1, got '{value}'"
                )));
            }
            self.cutoff = (!value.is_empty()).then_some(number);
        }
        self.split_output = params.get("split_output").is_some();
        for path in &paths {
            let path = Path::new(path);
            if !path.is_file() {
                return Err(PluginError::Init(format!(
                    "{NAME}: score file {} not found",
                    path.display()
                )));
            }
            // The INFO cell keeps its trailing whitespace, which Ensembl VEP's
            // `split /\t/` leaves in the value it prints (`SpliceAI.pm:326, 331`).
            self.stores.push(open_best_store_with(
                path,
                Self::vcf_config(path),
                self.assembly,
                TabixAnnotator::keeping_trailing_whitespace,
            )?);
        }
        self.version = Self::read_version(Path::new(paths[0]))?;
        Ok(())
    }

    fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        if variants.is_empty() || self.stores.is_empty() {
            return Ok(None);
        }
        let regions: Vec<(&str, u64, u64)> = variants
            .iter()
            .filter(|variant| self.annotates_variant(variant))
            .map(|variant| {
                let (start, end) = Self::query_bounds(variant);
                (table_contig(variant, self.assembly), start, end)
            })
            .collect();
        let batches = self
            .stores
            .iter()
            .map(|store| store.query_batch(&regions))
            .collect::<Result<Vec<BatchQueryResult>, PluginError>>()?;
        Ok(Some(PrefetchData::new(batches)))
    }

    fn annotate(
        &self,
        variants: &mut [InputVariant],
        data: PrefetchData,
    ) -> Result<(), PluginError> {
        let batches: Vec<BatchQueryResult> = data.downcast()?;
        for variant in variants
            .iter_mut()
            .filter(|variant| self.annotates_variant(variant))
        {
            let (start, end) = Self::query_bounds(variant);
            let records: Vec<&TabixRecord> = batches
                .iter()
                .flat_map(|batch| batch.lookup(table_contig(variant, self.assembly), start, end))
                .collect();
            if !records.is_empty() {
                self.annotate_variant(variant, &records);
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
        consequence: &TranscriptConsequence,
        variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        let Some(symbol) = consequence.gene_symbol() else {
            return Ok(IndexMap::new());
        };
        if self.stores.is_empty() || !self.annotates_row(consequence) {
            return Ok(IndexMap::new());
        }
        let records = self.query_stores(variant)?;
        let records: Vec<&TabixRecord> = records.iter().collect();
        // The fields alone: the row's empty-string marks are the batch path's.
        Ok(self
            .hits(variant, &records)
            .iter()
            .find(|hit| hit.symbol == symbol)
            .map(|hit| self.fields_of(hit).0)
            .unwrap_or_default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_fixtures::{write_bgzf_bytes, write_tabix_fixture, IndexSpec};
    use std::path::PathBuf;
    use std::sync::Arc;
    use vep_core::consequence::FeatureType;
    use vep_core::transcript::Transcript;

    const INFO_HEADER: &str = "##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAIv1.3 made-up scores of this module's tests\">";
    const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    /// A made-up score file of `rows` under `header`, bgzipped and indexed.
    fn score_file(dir: &Path, name: &str, header: &[&str], rows: &[&str]) -> PathBuf {
        let mut lines: Vec<&str> = header.to_vec();
        lines.extend_from_slice(rows);
        write_tabix_fixture(dir, name, &lines, IndexSpec::standard())
    }

    fn snv_file(dir: &Path) -> PathBuf {
        score_file(
            dir,
            "spliceai_snv.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
            &[
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|SOD1|0.10|0.20|0.30|0.42|1|2|-3|4",
                "21\t100\t.\tA\tT\t.\t.\tSpliceAI=T|SOD1|0.90|0.00|0.00|0.00|5|6|7|8",
                "21\t101\t.\tG\tA\t.\t.\tSpliceAI=A|SOD1|0.00|0.00|0.00|0.00|0|0|0|0",
                "21\t200\t.\tC\tT\t.\t.\tSpliceAI=T|GENE1|0.50|0.00|0.00|0.00|1|1|1|1",
                "21\t200\t.\tC\tT\t.\t.\tSpliceAI=T|GENE2|0.00|0.60|0.00|0.00|2|2|2|2",
                "21\t200\t.\tC\tT\t.\t.\tSpliceAI=T|GENE2|0.00|0.70|0.00|0.00|3|3|3|3",
            ],
        )
    }

    fn loaded(format: OutputFormat, params: &[String]) -> SpliceAiPlugin {
        let mut plugin = SpliceAiPlugin::new();
        plugin.set_output_format(format);
        plugin.init(params).unwrap();
        plugin
    }

    /// A transcript row whose gene symbol is `symbol`; `None` leaves the row
    /// without one.
    fn gene_row(symbol: Option<&str>) -> TranscriptConsequence {
        let mut consequence = TranscriptConsequence::default();
        consequence.transcript = Some(Arc::new({
            let mut transcript = Transcript::default();
            transcript.gene_symbol = symbol.map(Into::into);
            transcript
        }));
        consequence.feature_type = FeatureType::Transcript;
        consequence
    }

    fn snv(pos: u64, reference: &[u8], alt: &[u8], symbols: &[Option<&str>]) -> InputVariant {
        let mut variant =
            InputVariant::new("21".into(), pos, pos, reference.to_vec(), alt.to_vec());
        for symbol in symbols {
            variant.transcript_consequences.push(gene_row(*symbol));
        }
        variant
    }

    fn fields_of(variant: &InputVariant, row: usize) -> Vec<(&str, &str)> {
        variant.transcript_consequences[row]
            .plugin_data
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect()
    }

    #[test]
    fn new_plugin_is_uninitialized() {
        let plugin = SpliceAiPlugin::new();
        assert_eq!(plugin.name(), "SpliceAI");
        assert!(plugin.stores.is_empty());
        assert!(plugin.cutoff.is_none());
        assert!(!plugin.split_output);
        assert!(plugin.version.is_none());
    }

    /// The version is read from the `##INFO` lines as the bytes tabix prints
    /// them (SpliceAI.pm:196-205): a byte that is not UTF-8 (Latin-1 `é`,
    /// 0xE9) in the description leaves the `SpliceAIv1.3` token readable.
    #[test]
    fn a_header_that_is_not_utf8_still_names_the_version() {
        let dir = tempfile::tempdir().unwrap();
        let file = write_bgzf_bytes(
            dir.path(),
            "spliceai_snv.synthetic.vcf.gz",
            b"##fileformat=VCFv4.2\n##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAIv1.3 \xE9valuation\">\n#CHROM\tPOS\n",
        );
        assert_eq!(
            SpliceAiPlugin::read_version(&file).unwrap().as_deref(),
            Some("v1.3")
        );
    }

    #[test]
    fn run_returns_empty_when_not_initialized() {
        let plugin = SpliceAiPlugin::new();
        let consequence = gene_row(Some("SOD1"));
        let variant = InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert!(plugin.run(&consequence, &variant).unwrap().is_empty());
    }

    /// The default and tab formats print one field, the record's value after
    /// its allele, under the one header of `SpliceAI.pm:219` with the file's
    /// version (`:193-206`); the key=value grammar is `params_to_hash` alone.
    #[test]
    fn default_output_is_one_pipe_joined_field() {
        let dir = tempfile::tempdir().unwrap();
        let path = snv_file(dir.path());
        for format in [OutputFormat::Vep, OutputFormat::Tab, OutputFormat::Parquet] {
            let plugin = loaded(format, &[format!("snv={}", path.display())]);
            assert_eq!(
                plugin.header_info(),
                vec![(
                    "SpliceAI_pred".to_string(),
                    "SpliceAI (v1.3) predicted effect on splicing. These include delta scores (DS) \
                     and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor \
                     gain (DG), and donor loss (DL). Format: \
                     SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"
                        .to_string()
                )]
            );
            let mut variants = vec![snv(100, b"A", b"G", &[Some("SOD1")])];
            plugin.run_batch(&mut variants).unwrap();
            assert_eq!(
                fields_of(&variants[0], 0),
                [("SpliceAI_pred", "SOD1|0.10|0.20|0.30|0.42|1|2|-3|4")]
            );
        }
        let positional = SpliceAiPlugin::new().init(&[path.display().to_string()]);
        assert_eq!(
            positional.unwrap_err().to_string(),
            "plugin initialization failed: SpliceAI: no file specified; use snv=<path>"
        );
    }

    /// The VCF format and `split_output` print the nine fields under the nine
    /// headers of `SpliceAI.pm:209-217`; the JSON format writes the same nine
    /// fields for its writer to nest (`:275-288, 316`).
    #[test]
    fn split_output_and_vcf_emit_nine_fields() {
        let dir = tempfile::tempdir().unwrap();
        let path = snv_file(dir.path());
        let snv_param = format!("snv={}", path.display());
        let nine = [
            ("SpliceAI_pred_SYMBOL", "SOD1"),
            ("SpliceAI_pred_DS_AG", "0.10"),
            ("SpliceAI_pred_DS_AL", "0.20"),
            ("SpliceAI_pred_DS_DG", "0.30"),
            ("SpliceAI_pred_DS_DL", "0.42"),
            ("SpliceAI_pred_DP_AG", "1"),
            ("SpliceAI_pred_DP_AL", "2"),
            ("SpliceAI_pred_DP_DG", "-3"),
            ("SpliceAI_pred_DP_DL", "4"),
        ];
        for (format, params) in [
            (OutputFormat::Vcf, vec![snv_param.clone()]),
            (OutputFormat::Json, vec![snv_param.clone()]),
            (
                OutputFormat::Vep,
                vec![snv_param.clone(), "split_output=1".to_string()],
            ),
            (
                OutputFormat::Tab,
                vec![snv_param.clone(), "split_output=1".to_string()],
            ),
        ] {
            let plugin = loaded(format, &params);
            let header = plugin.header_info();
            let names: Vec<&str> = header.iter().map(|(name, _)| name.as_str()).collect();
            assert_eq!(names, SPLIT_FIELDS, "{format:?}");
            assert_eq!(header[0].1, "SpliceAI (v1.3) gene symbol");
            assert_eq!(
                header[1].1,
                "SpliceAI (v1.3) delta score (0–1) for acceptor gain: max inc. in predicted \
                 acceptor probability within ±50 bp. Guide cutoffs: <0.01 none; ≥0.2 \
                 high-recall; ≥0.5 recommended; ≥0.8 high-precision."
            );
            assert_eq!(
                header[8].1,
                "SpliceAI (v1.3) delta position (bp) for donor loss: position (relative to \
                 variant) of max DL change (+ downstream, − upstream)."
            );
            let mut variants = vec![snv(100, b"A", b"G", &[Some("SOD1")])];
            plugin.run_batch(&mut variants).unwrap();
            assert_eq!(fields_of(&variants[0], 0), nine, "{format:?}");
        }
    }

    /// `cutoff=` adds the `SpliceAI_cutoff` flag beside the scores and drops no
    /// record (`SpliceAI.pm:222-224, 295-303, 334-339`): the largest delta score
    /// `0.42` under `cutoff=0.5` is `FAIL`, `0.90` is `PASS`; the file's second
    /// matching record of a gene replaces the first (`:305`).
    #[test]
    fn cutoff_adds_pass_fail_flag_and_keeps_records() {
        let dir = tempfile::tempdir().unwrap();
        let path = snv_file(dir.path());
        let plugin = loaded(
            OutputFormat::Vep,
            &[format!("snv={}", path.display()), "cutoff=0.5".to_string()],
        );
        assert_eq!(
            plugin.header_info()[1],
            (
                "SpliceAI_cutoff".to_string(),
                "Flag if delta score pass the cutoff (PASS) or if it does not (FAIL)".to_string()
            )
        );
        let mut variants = vec![
            snv(100, b"A", b"G", &[Some("SOD1")]),
            snv(100, b"A", b"T", &[Some("SOD1")]),
            snv(200, b"C", b"T", &[Some("GENE1"), Some("GENE2")]),
        ];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(
            fields_of(&variants[0], 0),
            [
                ("SpliceAI_pred", "SOD1|0.10|0.20|0.30|0.42|1|2|-3|4"),
                ("SpliceAI_cutoff", "FAIL"),
            ]
        );
        assert_eq!(
            fields_of(&variants[1], 0),
            [
                ("SpliceAI_pred", "SOD1|0.90|0.00|0.00|0.00|5|6|7|8"),
                ("SpliceAI_cutoff", "PASS"),
            ]
        );
        assert_eq!(
            fields_of(&variants[2], 0),
            [
                ("SpliceAI_pred", "GENE1|0.50|0.00|0.00|0.00|1|1|1|1"),
                ("SpliceAI_cutoff", "PASS"),
            ]
        );
        assert_eq!(
            fields_of(&variants[2], 1),
            [
                ("SpliceAI_pred", "GENE2|0.00|0.70|0.00|0.00|3|3|3|3"),
                ("SpliceAI_cutoff", "PASS"),
            ]
        );
        // The flag sits inside the nine-field shape too, and a cutoff of zero
        // flags every record PASS instead of switching the flag off.
        let split = loaded(
            OutputFormat::Vcf,
            &[format!("snv={}", path.display()), "cutoff=0".to_string()],
        );
        let mut variants = vec![snv(101, b"G", b"A", &[Some("SOD1")])];
        split.run_batch(&mut variants).unwrap();
        let fields = fields_of(&variants[0], 0);
        assert_eq!(fields.len(), 10);
        assert_eq!(fields[9], ("SpliceAI_cutoff", "PASS"));
    }

    /// `cutoff=` sets no cutoff, the empty text passing the range check
    /// `$cutoff < 0 || $cutoff > 1` (`SpliceAI.pm:167`) and failing every later
    /// `if($self->{cutoff})` (`:222, 295, 335`); any other text is a cutoff at
    /// its number, as `$data_value->{info} >= $self->{cutoff}` (`:296`) reads
    /// it: `0.0`, `00` and `abc` are 0, `0.5x` and `0.5 ` are 0.5, `1abc` is 1,
    /// and `nan` fails every score, as no number reaches NaN. The literal `0` is
    /// a cutoff of 0 too, the documented divergence from Ensembl VEP, whose
    /// truth test drops it.
    #[test]
    fn empty_cutoff_is_none_and_any_other_text_is_its_perl_number() {
        let dir = tempfile::tempdir().unwrap();
        let path = snv_file(dir.path());
        let snv_param = format!("snv={}", path.display());
        // The flags of the records whose largest delta scores are 0.42, 0.90,
        // 0.00 and 0.50; `None` when the text sets no cutoff.
        let cases: [(&str, Option<[&str; 4]>); 9] = [
            ("", None),
            ("0", Some(["PASS", "PASS", "PASS", "PASS"])),
            ("0.0", Some(["PASS", "PASS", "PASS", "PASS"])),
            ("00", Some(["PASS", "PASS", "PASS", "PASS"])),
            ("abc", Some(["PASS", "PASS", "PASS", "PASS"])),
            ("0.5x", Some(["FAIL", "PASS", "FAIL", "PASS"])),
            ("0.5 ", Some(["FAIL", "PASS", "FAIL", "PASS"])),
            ("1abc", Some(["FAIL", "FAIL", "FAIL", "FAIL"])),
            ("nan", Some(["FAIL", "FAIL", "FAIL", "FAIL"])),
        ];
        for (value, flags) in cases {
            let plugin = loaded(
                OutputFormat::Vep,
                &[snv_param.clone(), format!("cutoff={value}")],
            );
            let header = plugin.header_info();
            assert_eq!(
                header.iter().any(|(name, _)| name == CUTOFF_FIELD),
                flags.is_some(),
                "{value:?}"
            );
            let mut variants = vec![
                snv(100, b"A", b"G", &[Some("SOD1")]),
                snv(100, b"A", b"T", &[Some("SOD1")]),
                snv(101, b"G", b"A", &[Some("SOD1")]),
                snv(200, b"C", b"T", &[Some("GENE1")]),
            ];
            plugin.run_batch(&mut variants).unwrap();
            for (index, variant) in variants.iter().enumerate() {
                let fields = fields_of(variant, 0);
                assert_eq!(fields[0].0, JOINED_FIELD, "{value:?} record {index}");
                assert_eq!(
                    fields.len(),
                    1 + usize::from(flags.is_some()),
                    "{value:?} record {index}"
                );
                let flag = fields
                    .iter()
                    .find(|(name, _)| *name == CUTOFF_FIELD)
                    .map(|(_, flag)| *flag);
                assert_eq!(
                    flag,
                    flags.map(|flags| flags[index]),
                    "{value:?} record {index}"
                );
            }
        }
    }

    /// `SpliceAI.pm:167-168` dies on a cutoff whose number is outside `0..=1`,
    /// the text numified by `$cutoff < 0 || $cutoff > 1`: `2abc` is 2.
    #[test]
    fn cutoff_outside_unit_interval_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = snv_file(dir.path());
        for value in ["1.5", "-0.1", "2abc"] {
            let error = SpliceAiPlugin::new()
                .init(&[format!("snv={}", path.display()), format!("cutoff={value}")])
                .unwrap_err();
            assert_eq!(
                error.to_string(),
                format!(
                    "plugin initialization failed: SpliceAI: cutoff must be between 0 and 1, got '{value}'"
                )
            );
        }
        let unknown = SpliceAiPlugin::new()
            .init(&[format!("snv={}", path.display()), "cut=0.5".to_string()])
            .unwrap_err();
        assert_eq!(
            unknown.to_string(),
            "plugin initialization failed: SpliceAI: unknown parameter 'cut'; SpliceAI accepts snv, indel, cutoff, split_output"
        );
        let none = SpliceAiPlugin::new().init(&[]).unwrap_err();
        assert_eq!(
            none.to_string(),
            "plugin initialization failed: SpliceAI: no score file given; use snv=<path> and indel=<path>"
        );
    }

    /// A record's scores reach the rows whose symbol equals the record's exactly
    /// (`SpliceAI.pm:313-317`, a hash lookup): a row of another gene, a row
    /// differing in case and a row without a symbol take nothing; a record at
    /// another position with the same alleles is not a match (`:272`).
    #[test]
    fn gene_symbol_match_is_exact_and_required() {
        let dir = tempfile::tempdir().unwrap();
        let path = snv_file(dir.path());
        let plugin = loaded(OutputFormat::Vep, &[format!("snv={}", path.display())]);
        let mut variants = vec![
            snv(
                100,
                b"A",
                b"G",
                &[Some("SOD1"), Some("sod1"), Some("APP"), None],
            ),
            snv(99, b"A", b"G", &[Some("SOD1")]),
        ];
        plugin.run_batch(&mut variants).unwrap();
        let symbols: Vec<Option<&str>> = (0..4)
            .map(|row| {
                variants[0].transcript_consequences[row]
                    .plugin_data
                    .get("SpliceAI_pred")
                    .map(String::as_str)
            })
            .collect();
        assert_eq!(
            symbols,
            [Some("SOD1|0.10|0.20|0.30|0.42|1|2|-3|4"), None, None, None]
        );
        assert!(variants[1].transcript_consequences[0]
            .plugin_data
            .is_empty());
        assert!(variants[0].plugin_data.is_empty());
    }

    /// Both files are read for every variant (`BaseVepTabixPlugin.pm:591`), an
    /// indel's VCF record matching from its anchor base, and the header version
    /// comes from the first file, `v1.3.1` keeping every dot of the token
    /// (`SpliceAI.pm:195-202`); a file whose header names no version writes the
    /// descriptions without the parenthetical.
    #[test]
    fn indel_file_is_read_for_every_variant_and_names_the_version() {
        let dir = tempfile::tempdir().unwrap();
        let snv_path = score_file(
            dir.path(),
            "snv.synthetic.vcf.gz",
            &[
                "##fileformat=VCFv4.2",
                &INFO_HEADER.replace("SpliceAIv1.3 ", "SpliceAIv1.3.1 "),
                COLUMNS,
            ],
            &["21\t300\t.\tA\tG\t.\t.\tSpliceAI=G|APP|0.10|0.10|0.10|0.10|1|1|1|1"],
        );
        let indel_path = score_file(
            dir.path(),
            "indel.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", COLUMNS],
            &[
                "21\t100\t.\tAC\tA\t.\t.\tSpliceAI=A|SOD1|0.11|0.12|0.13|0.14|-1|-2|-3|-4",
                "21\t100\t.\tA\tAT\t.\t.\tSpliceAI=AT|SOD1|0.21|0.22|0.23|0.24|5|6|7|8",
            ],
        );
        let plugin = loaded(
            OutputFormat::Vep,
            &[
                format!("indel={}", indel_path.display()),
                format!("snv={}", snv_path.display()),
            ],
        );
        assert!(plugin.header_info()[0]
            .1
            .starts_with("SpliceAI (v1.3.1) predicted"));
        let mut deletion = InputVariant::new("21".into(), 101, 101, b"C".to_vec(), b"-".to_vec());
        let mut insertion = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"T".to_vec());
        for variant in [&mut deletion, &mut insertion] {
            variant.transcript_consequences.push(gene_row(Some("SOD1")));
        }
        let mut variants = vec![deletion, insertion, snv(300, b"A", b"G", &[Some("APP")])];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(
            fields_of(&variants[0], 0),
            [("SpliceAI_pred", "SOD1|0.11|0.12|0.13|0.14|-1|-2|-3|-4")]
        );
        assert_eq!(
            fields_of(&variants[1], 0),
            [("SpliceAI_pred", "SOD1|0.21|0.22|0.23|0.24|5|6|7|8")]
        );
        assert_eq!(
            fields_of(&variants[2], 0),
            [("SpliceAI_pred", "APP|0.10|0.10|0.10|0.10|1|1|1|1")]
        );
        let unversioned = loaded(
            OutputFormat::Vcf,
            &[format!("indel={}", indel_path.display())],
        );
        assert_eq!(unversioned.header_info()[0].1, "SpliceAI gene symbol");
    }

    /// A multi-allelic record whose alleles share no first base keeps its
    /// insertion or deletion allele as written (`C/G/CA` at one site, `CT/C/GT`
    /// at another), and Ensembl VEP's plugin compares that allele with the
    /// file's record verbatim (`SpliceAI.pm:272`): the record matches as it
    /// does in the anchored form of a bi-allelic record; another insertion at
    /// the site, and the record's SNV allele, which no record of this file
    /// names, take nothing.
    #[test]
    fn multi_allelic_indel_allele_matches_as_written() {
        let dir = tempfile::tempdir().unwrap();
        let indel_path = score_file(
            dir.path(),
            "indel.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
            &[
                "21\t100\t.\tC\tCA\t.\t.\tSpliceAI=CA|SOD1|0.21|0.22|0.23|0.24|5|6|7|8",
                "21\t200\t.\tCT\tC\t.\t.\tSpliceAI=C|SOD1|0.11|0.12|0.13|0.14|-1|-2|-3|-4",
            ],
        );
        let plugin = loaded(
            OutputFormat::Vep,
            &[format!("indel={}", indel_path.display())],
        );
        let mut variants = vec![
            InputVariant::new("21".into(), 100, 100, b"C".to_vec(), b"CA".to_vec()),
            InputVariant::new("21".into(), 100, 100, b"C".to_vec(), b"CG".to_vec()),
            InputVariant::new("21".into(), 100, 100, b"C".to_vec(), b"G".to_vec()),
            InputVariant::new("21".into(), 200, 201, b"CT".to_vec(), b"C".to_vec()),
            InputVariant::new("21".into(), 200, 201, b"CT".to_vec(), b"GT".to_vec()),
        ];
        for variant in &mut variants {
            variant.transcript_consequences.push(gene_row(Some("SOD1")));
        }
        plugin.run_batch(&mut variants).unwrap();
        let written: Vec<Option<&str>> = variants
            .iter()
            .map(|variant| {
                variant.transcript_consequences[0]
                    .plugin_data
                    .get("SpliceAI_pred")
                    .map(String::as_str)
            })
            .collect();
        assert_eq!(
            written,
            [
                Some("SOD1|0.21|0.22|0.23|0.24|5|6|7|8"),
                None,
                None,
                Some("SOD1|0.11|0.12|0.13|0.14|-1|-2|-3|-4"),
                None,
            ]
        );
    }

    /// `run` answers one row as `run_batch` writes it, and `SpliceAI.pm:184-186`
    /// declares `['Transcript']`: a regulatory or motif row takes nothing.
    #[test]
    fn writes_on_transcript_rows_alone() {
        let dir = tempfile::tempdir().unwrap();
        let path = snv_file(dir.path());
        let plugin = loaded(OutputFormat::Vep, &[format!("snv={}", path.display())]);
        let mut variants = vec![snv(100, b"A", b"G", &[Some("SOD1")])];
        for feature_type in [FeatureType::RegulatoryFeature, FeatureType::MotifFeature] {
            variants[0].transcript_consequences.push({
                let mut consequence = gene_row(Some("SOD1"));
                consequence.feature_type = feature_type;
                consequence
            });
        }
        let through_run = plugin
            .run(&variants[0].transcript_consequences[0], &variants[0])
            .unwrap();
        plugin.run_batch(&mut variants).unwrap();
        let written: Vec<bool> = variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.contains_key("SpliceAI_pred"))
            .collect();
        assert_eq!(written, [true, false, false]);
        assert_eq!(
            through_run,
            variants[0].transcript_consequences[0].plugin_data
        );
        assert!(plugin
            .run(&variants[0].transcript_consequences[1], &variants[0])
            .unwrap()
            .is_empty());
    }

    /// Ensembl VEP 116.3 never hands SpliceAI a structural variant: the module
    /// keeps the default `variant_feature_types` of `['VariationFeature']`
    /// (BaseVepPlugin.pm:100), which the dispatch checks before running a
    /// plugin on a line (OutputFactory.pm:2229). A `<DEL>` row gets no field,
    /// and nor does a scored substitution carrying the structural mark, which
    /// pins the gate rather than the allele match; the plain one is scored.
    #[test]
    fn a_structural_variant_row_gets_no_field() {
        let dir = tempfile::tempdir().unwrap();
        let path = snv_file(dir.path());
        let plugin = loaded(OutputFormat::Vep, &[format!("snv={}", path.display())]);
        assert!(!plugin.annotates_structural_variants());
        let mut deletion = snv(100, b"N", b"<DEL>", &[Some("SOD1")]);
        deletion.end = 150;
        deletion.is_structural = true;
        let mut flagged = snv(100, b"A", b"G", &[Some("SOD1")]);
        flagged.is_structural = true;
        let mut variants = vec![deletion, flagged, snv(100, b"A", b"G", &[Some("SOD1")])];
        plugin.run_batch(&mut variants).unwrap();
        assert!(fields_of(&variants[0], 0).is_empty());
        assert!(fields_of(&variants[1], 0).is_empty());
        assert_eq!(
            fields_of(&variants[2], 0),
            [("SpliceAI_pred", "SOD1|0.10|0.20|0.30|0.42|1|2|-3|4")]
        );
    }

    /// `SpliceAI.pm:326` splits the line htslib hands over, the line end alone
    /// removed, and `:331` prints the value after the allele as split: a
    /// trailing space on the INFO cell stays in `SpliceAI_pred`, and in the
    /// last of the nine fields of the split shape (`:276, 287`). A
    /// whitespace-only INFO cell names no gene and writes nothing.
    #[test]
    fn a_trailing_space_on_the_info_cell_is_printed_as_ensembl_vep_prints_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = score_file(
            dir.path(),
            "trailing.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
            &[
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|SOD1|0.10|0.20|0.30|0.42|1|2|-3|4 ",
                "21\t200\t.\tC\tT\t.\t.\t ",
            ],
        );
        let joined = loaded(OutputFormat::Vep, &[format!("snv={}", path.display())]);
        let mut variants = vec![
            snv(100, b"A", b"G", &[Some("SOD1")]),
            snv(200, b"C", b"T", &[Some("SOD1")]),
        ];
        joined.run_batch(&mut variants).unwrap();
        assert_eq!(
            fields_of(&variants[0], 0),
            [("SpliceAI_pred", "SOD1|0.10|0.20|0.30|0.42|1|2|-3|4 ")]
        );
        assert!(fields_of(&variants[1], 0).is_empty());

        let split = loaded(OutputFormat::Vcf, &[format!("snv={}", path.display())]);
        let mut variants = vec![snv(100, b"A", b"G", &[Some("SOD1")])];
        split.run_batch(&mut variants).unwrap();
        let fields = fields_of(&variants[0], 0);
        assert_eq!(fields.len(), 9);
        assert_eq!(fields[8], ("SpliceAI_pred_DP_DL", "4 "));
    }

    /// A `Number=.` INFO cell joining two genes' entries with a comma scores
    /// each gene from its own entry, where `SpliceAI.pm:328-332` reads the
    /// whole cell as one entry and Ensembl VEP 116.3 prints
    /// `GENE1|0.10|0.20|0.30|0.40|1|2|3|4,G|GENE2|0.50|0.60|0.70|0.80|5|6|7|8`
    /// on the GENE1 rows, `DP_DL` as `4,G`, and nothing on the GENE2 rows.
    #[test]
    fn two_gene_entries_in_one_record_score_each_gene() {
        let dir = tempfile::tempdir().unwrap();
        let path = score_file(
            dir.path(),
            "two_genes.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
            &["21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|0.10|0.20|0.30|0.40|1|2|3|4,G|GENE2|0.50|0.60|0.70|0.80|5|6|7|8"],
        );
        let snv_param = format!("snv={}", path.display());
        let rows = &[Some("GENE1"), Some("GENE2"), Some("GENE3")];

        let joined = loaded(OutputFormat::Vep, std::slice::from_ref(&snv_param));
        let mut variants = vec![snv(100, b"A", b"G", rows)];
        joined.run_batch(&mut variants).unwrap();
        assert_eq!(
            fields_of(&variants[0], 0),
            [("SpliceAI_pred", "GENE1|0.10|0.20|0.30|0.40|1|2|3|4")]
        );
        assert_eq!(
            fields_of(&variants[0], 1),
            [("SpliceAI_pred", "GENE2|0.50|0.60|0.70|0.80|5|6|7|8")]
        );
        assert!(fields_of(&variants[0], 2).is_empty());

        let split = loaded(OutputFormat::Vcf, std::slice::from_ref(&snv_param));
        let mut variants = vec![snv(100, b"A", b"G", rows)];
        split.run_batch(&mut variants).unwrap();
        assert_eq!(fields_of(&variants[0], 0)[8], ("SpliceAI_pred_DP_DL", "4"));
        assert_eq!(fields_of(&variants[0], 1)[8], ("SpliceAI_pred_DP_DL", "8"));
        assert!(fields_of(&variants[0], 2).is_empty());

        let flagged = loaded(OutputFormat::Vep, &[snv_param, "cutoff=0.5".to_string()]);
        let mut variants = vec![snv(100, b"A", b"G", rows)];
        flagged.run_batch(&mut variants).unwrap();
        assert_eq!(fields_of(&variants[0], 0)[1], ("SpliceAI_cutoff", "FAIL"));
        assert_eq!(fields_of(&variants[0], 1)[1], ("SpliceAI_cutoff", "PASS"));
        assert!(fields_of(&variants[0], 2).is_empty());
    }

    /// Three entries of one made-up record, two naming one gene and one
    /// another, in each of their three orders: a row carries every entry that
    /// names its gene, in the record's order, and no entry of the other gene,
    /// where `SpliceAI.pm:328-332` files the whole cell under the first
    /// entry's gene.
    #[test]
    fn every_entry_of_a_gene_stands_on_its_rows_in_the_record_s_order() {
        let first = "G|GENE1|0.10|0.20|0.30|0.40|1|2|3|4";
        let second = "G|GENE1|0.11|0.21|0.31|0.41|-1|-2|-3|-4";
        let other = "G|GENE2|0.50|0.60|0.70|0.80|5|6|7|8";
        let rows = &[Some("GENE1"), Some("GENE2"), Some("GENE3")];
        for order in [
            [first, second, other],
            [first, other, second],
            [other, first, second],
        ] {
            let dir = tempfile::tempdir().unwrap();
            let line = format!("21\t100\t.\tA\tG\t.\t.\tSpliceAI={}", order.join(","));
            let path = score_file(
                dir.path(),
                "three_entries.synthetic.vcf.gz",
                &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
                &[&line],
            );
            let plugin = loaded(OutputFormat::Vep, &[format!("snv={}", path.display())]);
            let mut variants = vec![snv(100, b"A", b"G", rows)];
            plugin.run_batch(&mut variants).unwrap();
            assert_eq!(
                fields_of(&variants[0], 0),
                [(
                    "SpliceAI_pred",
                    "GENE1|0.10|0.20|0.30|0.40|1|2|3|4,G|GENE1|0.11|0.21|0.31|0.41|-1|-2|-3|-4"
                )],
                "{order:?}"
            );
            assert_eq!(
                fields_of(&variants[0], 1),
                [("SpliceAI_pred", "GENE2|0.50|0.60|0.70|0.80|5|6|7|8")],
                "{order:?}"
            );
            assert!(fields_of(&variants[0], 2).is_empty(), "{order:?}");
        }
    }

    /// A made-up record that names one gene alone in two and in four entries:
    /// the one field of the default and tab formats is the record's value
    /// behind its first allele, each later entry with its allele part, the
    /// value `SpliceAI.pm:328-332` keeps and `:291` prints.
    #[test]
    fn a_later_entry_of_a_gene_keeps_its_allele_part_in_the_one_field_form() {
        let dir = tempfile::tempdir().unwrap();
        let path = score_file(
            dir.path(),
            "one_gene.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
            &[
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|SOD1|0.10|0.20|0.30|0.40|1|2|3|4,G|SOD1|0.11|0.21|0.31|0.41|-1|-2|-3|-4",
                "21\t200\t.\tC\tT\t.\t.\tSpliceAI=T|SOD1|0.01|0.02|0.03|0.04|1|2|3|4,T|SOD1|0.01|0.02|0.03|0.04|1|2|3|4,T|SOD1|0.05|0.06|0.07|0.08|5|6|7|8,T|SOD1|0.09|0.10|0.11|0.12|9|10|11|12",
            ],
        );
        for format in [OutputFormat::Vep, OutputFormat::Tab] {
            let plugin = loaded(format, &[format!("snv={}", path.display())]);
            let mut variants = vec![
                snv(100, b"A", b"G", &[Some("SOD1")]),
                snv(200, b"C", b"T", &[Some("SOD1")]),
            ];
            plugin.run_batch(&mut variants).unwrap();
            assert_eq!(
                fields_of(&variants[0], 0),
                [(
                    "SpliceAI_pred",
                    "SOD1|0.10|0.20|0.30|0.40|1|2|3|4,G|SOD1|0.11|0.21|0.31|0.41|-1|-2|-3|-4"
                )],
                "{format:?}"
            );
            assert_eq!(
                fields_of(&variants[1], 0),
                [(
                    "SpliceAI_pred",
                    "SOD1|0.01|0.02|0.03|0.04|1|2|3|4,T|SOD1|0.01|0.02|0.03|0.04|1|2|3|4,T|SOD1|0.05|0.06|0.07|0.08|5|6|7|8,T|SOD1|0.09|0.10|0.11|0.12|9|10|11|12"
                )],
                "{format:?}"
            );
        }
    }

    /// In the nine-field shape each member carries its values over the
    /// entries of a made-up record that name the row's gene, in the record's
    /// order, joined by a comma, where the ninth variable of `SpliceAI.pm:287`
    /// holds `4,G` and the later entries are dropped; a gene named by one
    /// entry keeps its nine values.
    #[test]
    fn the_nine_members_join_their_values_over_a_gene_s_entries() {
        let dir = tempfile::tempdir().unwrap();
        let path = score_file(
            dir.path(),
            "members.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
            &["21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|0.10|0.20|0.30|0.40|1|2|3|4,G|GENE2|0.50|0.60|0.70|0.80|5|6|7|8,G|GENE1|0.11|0.21|0.31|0.41|-1|-2|-3|-4"],
        );
        let snv_param = format!("snv={}", path.display());
        let rows = &[Some("GENE1"), Some("GENE2")];
        let joined = [
            ("SpliceAI_pred_SYMBOL", "GENE1,GENE1"),
            ("SpliceAI_pred_DS_AG", "0.10,0.11"),
            ("SpliceAI_pred_DS_AL", "0.20,0.21"),
            ("SpliceAI_pred_DS_DG", "0.30,0.31"),
            ("SpliceAI_pred_DS_DL", "0.40,0.41"),
            ("SpliceAI_pred_DP_AG", "1,-1"),
            ("SpliceAI_pred_DP_AL", "2,-2"),
            ("SpliceAI_pred_DP_DG", "3,-3"),
            ("SpliceAI_pred_DP_DL", "4,-4"),
        ];
        let single = [
            ("SpliceAI_pred_SYMBOL", "GENE2"),
            ("SpliceAI_pred_DS_AG", "0.50"),
            ("SpliceAI_pred_DS_AL", "0.60"),
            ("SpliceAI_pred_DS_DG", "0.70"),
            ("SpliceAI_pred_DS_DL", "0.80"),
            ("SpliceAI_pred_DP_AG", "5"),
            ("SpliceAI_pred_DP_AL", "6"),
            ("SpliceAI_pred_DP_DG", "7"),
            ("SpliceAI_pred_DP_DL", "8"),
        ];
        for (format, split_output) in [
            (OutputFormat::Vcf, false),
            (OutputFormat::Json, false),
            (OutputFormat::Vep, true),
            (OutputFormat::Tab, true),
        ] {
            let mut params = vec![snv_param.clone()];
            if split_output {
                params.push("split_output=1".to_string());
            }
            let plugin = loaded(format, &params);
            let mut variants = vec![snv(100, b"A", b"G", rows)];
            plugin.run_batch(&mut variants).unwrap();
            assert_eq!(fields_of(&variants[0], 0), joined, "{format:?}");
            assert_eq!(fields_of(&variants[0], 1), single, "{format:?}");
            assert!(variants[0].transcript_consequences[0]
                .plugin_empty_strings
                .is_empty());
        }
    }

    /// Under a cutoff a row carries one verdict, in the one-field and the
    /// nine-field shape: `PASS` where any entry of a made-up record that names
    /// the row's gene reaches the cutoff by its own largest delta score, also
    /// a later entry alone, and `FAIL` where none does, an entry of another
    /// gene counting for nothing. `SpliceAI.pm:336-338` reads the first four
    /// score parts of the whole cell, the first entry's.
    #[test]
    fn the_cutoff_flag_is_one_verdict_over_a_gene_s_entries() {
        let dir = tempfile::tempdir().unwrap();
        let path = score_file(
            dir.path(),
            "verdicts.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
            &[
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|SOD1|0.10|0.20|0.30|0.40|1|2|3|4,G|SOD1|0.10|0.90|0.30|0.40|1|2|3|4",
                "21\t200\t.\tC\tT\t.\t.\tSpliceAI=T|SOD1|0.10|0.20|0.30|0.40|1|2|3|4,T|SOD1|0.11|0.21|0.31|0.41|1|2|3|4",
                "21\t300\t.\tG\tA\t.\t.\tSpliceAI=A|SOD1|0.70|0.20|0.30|0.40|1|2|3|4,A|SOD1|0.11|0.21|0.31|0.41|1|2|3|4",
                "21\t400\t.\tT\tC\t.\t.\tSpliceAI=C|SOD1|0.10|0.20|0.30|0.40|1|2|3|4,C|GENE2|0.95|0.21|0.31|0.41|1|2|3|4,C|SOD1|0.11|0.21|0.31|0.41|1|2|3|4",
            ],
        );
        let snv_param = format!("snv={}", path.display());
        let verdicts = ["PASS", "FAIL", "PASS", "FAIL"];
        for format in [OutputFormat::Vep, OutputFormat::Vcf, OutputFormat::Json] {
            let plugin = loaded(format, &[snv_param.clone(), "cutoff=0.5".to_string()]);
            let mut variants = vec![
                snv(100, b"A", b"G", &[Some("SOD1")]),
                snv(200, b"C", b"T", &[Some("SOD1")]),
                snv(300, b"G", b"A", &[Some("SOD1")]),
                snv(400, b"T", b"C", &[Some("SOD1"), Some("GENE2")]),
            ];
            plugin.run_batch(&mut variants).unwrap();
            for (variant, verdict) in variants.iter().zip(verdicts) {
                let fields = fields_of(variant, 0);
                assert_eq!(
                    fields.last(),
                    Some(&("SpliceAI_cutoff", verdict)),
                    "{format:?}"
                );
                let flags = fields
                    .iter()
                    .filter(|(name, _)| *name == "SpliceAI_cutoff")
                    .count();
                assert_eq!(flags, 1, "{format:?}");
            }
            assert_eq!(
                fields_of(&variants[3], 1).last(),
                Some(&("SpliceAI_cutoff", "PASS")),
                "{format:?}"
            );
        }
    }

    /// Where several made-up record lines stand at one allele, the last line
    /// that names a gene is its rows' record (`SpliceAI.pm:305`): its entries
    /// of the gene take the place of an earlier line's, and a gene the last
    /// line leaves out keeps the entries of the last line that names it.
    #[test]
    fn a_later_record_line_of_a_gene_replaces_the_earlier_line_s_entries() {
        let dir = tempfile::tempdir().unwrap();
        let path = score_file(
            dir.path(),
            "lines.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
            &[
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|0.10|0.20|0.30|0.40|1|2|3|4,G|GENE1|0.11|0.21|0.31|0.41|-1|-2|-3|-4",
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|0.12|0.22|0.32|0.42|9|9|9|9",
                "21\t200\t.\tC\tT\t.\t.\tSpliceAI=T|GENE1|0.10|0.20|0.30|0.40|1|2|3|4,T|GENE2|0.50|0.60|0.70|0.80|5|6|7|8",
                "21\t200\t.\tC\tT\t.\t.\tSpliceAI=T|GENE1|0.13|0.23|0.33|0.43|6|6|6|6,T|GENE1|0.14|0.24|0.34|0.44|7|7|7|7",
            ],
        );
        let plugin = loaded(OutputFormat::Vep, &[format!("snv={}", path.display())]);
        let rows = &[Some("GENE1"), Some("GENE2")];
        let mut variants = vec![snv(100, b"A", b"G", rows), snv(200, b"C", b"T", rows)];
        plugin.run_batch(&mut variants).unwrap();
        assert_eq!(
            fields_of(&variants[0], 0),
            [("SpliceAI_pred", "GENE1|0.12|0.22|0.32|0.42|9|9|9|9")]
        );
        assert!(fields_of(&variants[0], 1).is_empty());
        assert_eq!(
            fields_of(&variants[1], 0),
            [(
                "SpliceAI_pred",
                "GENE1|0.13|0.23|0.33|0.43|6|6|6|6,T|GENE1|0.14|0.24|0.34|0.44|7|7|7|7"
            )]
        );
        assert_eq!(
            fields_of(&variants[1], 1),
            [("SpliceAI_pred", "GENE2|0.50|0.60|0.70|0.80|5|6|7|8")]
        );
    }

    /// A truncated entry is kept as `parse_data` keeps it (SpliceAI.pm:326-349)
    /// and the last entry of a gene wins (`:305`), so a four-part entry after a
    /// whole one leaves the gene's rows its three parts: one joined field in
    /// the default format, and in the split shape the nine fields with the
    /// missing parts as the undefined value and an empty part before a filled
    /// one as the empty string, named in the row's `plugin_empty_strings`; a
    /// two-part entry is its symbol and a trailing `|`; a one-part entry reaches
    /// no row; an entry with an empty allele part is kept, the match being on
    /// the REF and ALT columns (`:272`).
    #[test]
    fn a_truncated_entry_is_kept_with_the_parts_it_has() {
        let dir = tempfile::tempdir().unwrap();
        let path = score_file(
            dir.path(),
            "spliceai_snv.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
            &[
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|SOD1|0.10|0.20|0.30|0.42|1|2|-3|4",
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|SOD1|0.50|0.60",
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=G|GENE2",
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=G",
                "21\t100\t.\tA\tG\t.\t.\tSpliceAI=|GENE3|0.98||0.49|0.83|4|-33|-28|",
            ],
        );
        let snv_param = format!("snv={}", path.display());
        let symbols = [Some("SOD1"), Some("GENE2"), Some("GENE3")];

        let joined = loaded(OutputFormat::Vep, std::slice::from_ref(&snv_param));
        let mut variants = vec![snv(100, b"A", b"G", &symbols)];
        joined.run_batch(&mut variants).unwrap();
        assert_eq!(
            fields_of(&variants[0], 0),
            [("SpliceAI_pred", "SOD1|0.50|0.60")]
        );
        assert_eq!(fields_of(&variants[0], 1), [("SpliceAI_pred", "GENE2|")]);
        assert_eq!(
            fields_of(&variants[0], 2),
            [("SpliceAI_pred", "GENE3|0.98||0.49|0.83|4|-33|-28|")]
        );

        let split = loaded(OutputFormat::Vcf, &[snv_param]);
        let mut variants = vec![snv(100, b"A", b"G", &symbols)];
        split.run_batch(&mut variants).unwrap();
        assert_eq!(
            fields_of(&variants[0], 0),
            [
                ("SpliceAI_pred_SYMBOL", "SOD1"),
                ("SpliceAI_pred_DS_AG", "0.50"),
                ("SpliceAI_pred_DS_AL", "0.60"),
                ("SpliceAI_pred_DS_DG", ""),
                ("SpliceAI_pred_DS_DL", ""),
                ("SpliceAI_pred_DP_AG", ""),
                ("SpliceAI_pred_DP_AL", ""),
                ("SpliceAI_pred_DP_DG", ""),
                ("SpliceAI_pred_DP_DL", ""),
            ]
        );
        assert!(variants[0].transcript_consequences[0]
            .plugin_empty_strings
            .is_empty());
        assert_eq!(
            fields_of(&variants[0], 1)[0],
            ("SpliceAI_pred_SYMBOL", "GENE2")
        );
        assert_eq!(fields_of(&variants[0], 1)[1], ("SpliceAI_pred_DS_AG", ""));
        assert_eq!(fields_of(&variants[0], 1).len(), 9);
        let gene3 = fields_of(&variants[0], 2);
        assert_eq!(gene3[2], ("SpliceAI_pred_DS_AL", ""));
        assert_eq!(gene3[7], ("SpliceAI_pred_DP_DG", "-28"));
        assert_eq!(gene3[8], ("SpliceAI_pred_DP_DL", ""));
        assert_eq!(
            variants[0].transcript_consequences[2].plugin_empty_strings,
            ["SpliceAI_pred_DS_AL"]
        );
    }

    /// Every record and INFO cell of this test is made up.
    #[test]
    fn score_values_read_every_entry_of_a_well_formed_info() {
        let record = |info: &str| TabixRecord {
            chr: "21".into(),
            start: 1,
            end: 1,
            columns: vec![
                "21".into(),
                "1".into(),
                ".".into(),
                "A".into(),
                "G".into(),
                ".".into(),
                ".".into(),
                info.into(),
            ],
            header: None,
        };
        let values = |info: &str| -> Vec<String> {
            SpliceAiPlugin::score_values(&record(info))
                .map(|entry| entry.scores.into_owned())
                .collect()
        };
        assert_eq!(
            values(
                "SpliceAI=G|GENE1|0.10|0.20|0.30|0.40|1|2|3|4,G|GENE2|0.50|0.60|0.70|0.80|5|6|7|8"
            ),
            [
                "GENE1|0.10|0.20|0.30|0.40|1|2|3|4",
                "GENE2|0.50|0.60|0.70|0.80|5|6|7|8"
            ]
        );
        assert_eq!(
            values("AF=0.01;SpliceAI=G|TP53|0.05|0.00|0.00|0.10|10|-5|3|20;DB"),
            ["TP53|0.05|0.00|0.00|0.10|10|-5|3|20"]
        );
        assert!(values("").is_empty());
        assert!(values("AF=0.01").is_empty());
        // A truncated entry is kept with the parts it has (SpliceAI.pm:326-349):
        // two parts are the symbol and a trailing `|`; one part has no symbol.
        assert_eq!(values("SpliceAI=T|GENE"), ["GENE|"]);
        assert_eq!(values("SpliceAI=T|GENE|0.5|0.6"), ["GENE|0.5|0.6"]);
        assert_eq!(
            values("SpliceAI=|GENE|0.5|0.6|0.7|0.8|1|2|3|4"),
            ["GENE|0.5|0.6|0.7|0.8|1|2|3|4"]
        );
        assert!(values("SpliceAI=T").is_empty());
        let symbols = |info: &str| -> Vec<String> {
            SpliceAiPlugin::score_values(&record(info))
                .map(|entry| entry.symbol.to_string())
                .collect()
        };
        assert_eq!(symbols("SpliceAI=T|GENE,T|GENE2|0.5"), ["GENE", "GENE2"]);
        assert_eq!(
            SpliceAiPlugin::max_delta_score("SOD1|0.10|0.50|0.20|0.30|1|2|3|4"),
            0.5
        );
        assert_eq!(SpliceAiPlugin::max_delta_score("SOD1|.|.|.|.|1|2|3|4"), 0.0);
    }

    /// A record's POS cell is read as `$start == $data_value->{start}`
    /// (SpliceAI.pm:272) numifies it: `100.0`, `100x` and `1e2` are 100 and the
    /// record scores the variant there, for a substitution and for a deletion
    /// read through its anchored form alike; `100.5`, `abc` and `0` are no
    /// position. The tabix index places a cell at its leading integer, so the
    /// file carries the first shapes and the exponent and the rest are read
    /// from records built here.
    #[test]
    fn a_position_cell_is_matched_as_perl_numifies_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = score_file(
            dir.path(),
            "cells.synthetic.vcf.gz",
            &["##fileformat=VCFv4.2", INFO_HEADER, COLUMNS],
            &[
                "21\t100.0\t.\tA\tG\t.\t.\tSpliceAI=G|SOD1|0.10|0.20|0.30|0.40|1|2|3|4",
                "21\t100.5\t.\tA\tT\t.\t.\tSpliceAI=T|SOD1|0.50|0.60|0.70|0.80|5|6|7|8",
                "21\t100x\t.\tA\tC\t.\t.\tSpliceAI=C|SOD1|0.11|0.21|0.31|0.41|1|1|1|1",
                "21\t200.0\t.\tCT\tC\t.\t.\tSpliceAI=C|SOD1|0.90|0.00|0.00|0.00|9|9|9|9",
                "21\t300.5\t.\tCT\tC\t.\t.\tSpliceAI=C|SOD1|0.91|0.00|0.00|0.00|9|9|9|9",
            ],
        );
        let plugin = loaded(OutputFormat::Vep, &[format!("snv={}", path.display())]);
        let scored = |variant: InputVariant| {
            let mut variants = vec![variant];
            plugin.run_batch(&mut variants).unwrap();
            fields_of(&variants[0], 0)
                .first()
                .map(|(_, value)| value.to_string())
        };
        let symbol = &[Some("SOD1")];
        assert_eq!(
            scored(snv(100, b"A", b"G", symbol)).as_deref(),
            Some("SOD1|0.10|0.20|0.30|0.40|1|2|3|4")
        );
        assert_eq!(scored(snv(100, b"A", b"T", symbol)), None);
        assert_eq!(
            scored(snv(100, b"A", b"C", symbol)).as_deref(),
            Some("SOD1|0.11|0.21|0.31|0.41|1|1|1|1")
        );
        let deletion = |start: u64| {
            let mut variant =
                InputVariant::new("21".into(), start, start, b"T".to_vec(), b"-".to_vec());
            variant.transcript_consequences.push(gene_row(Some("SOD1")));
            variant
        };
        assert_eq!(
            scored(deletion(201)).as_deref(),
            Some("SOD1|0.90|0.00|0.00|0.00|9|9|9|9")
        );
        assert_eq!(scored(deletion(301)), None);

        let record = |pos: &str| TabixRecord {
            chr: "21".into(),
            start: 100,
            end: 100,
            columns: [
                "21",
                pos,
                ".",
                "A",
                "G",
                ".",
                ".",
                "SpliceAI=G|SOD1|0.1|0.2|0.3|0.4|1|2|3|4",
            ]
            .iter()
            .map(|cell| cell.to_string())
            .collect(),
            header: None,
        };
        let is_variant_record = |pos: &str| {
            SpliceAiPlugin::is_variant_record(&record(pos), &snv(100, b"A", b"G", symbol))
        };
        assert!(is_variant_record("1e2"));
        assert!(is_variant_record(" 100"));
        assert!(!is_variant_record("abc"));
        assert!(!is_variant_record("0"));
        assert!(!is_variant_record("100.5"));
    }
}
