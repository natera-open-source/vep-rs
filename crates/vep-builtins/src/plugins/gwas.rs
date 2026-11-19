// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! GWAS plugin: the NHGRI-EBI GWAS Catalog's association at a variant's site, as
//! Ensembl VEP's GWAS plugin (`GWAS.pm` at VEP_plugins release/116) reads it from
//! the catalog's curated associations download or from a harmonised
//! summary-statistics file.
//!
//! Every parameter is `key=value` (GWAS.pm:94-104): `file=<path>`, mandatory;
//! `type=curated` (the default) or `type=sstate`; `verbose=1` names every
//! catalog row the plugin skips; `db=1`, which reads variant locations from the
//! Ensembl database, is refused.
//!
//! Output fields, written on every row of a record whose site the file names
//! with a matching allele (GWAS.pm:150-152, 170-199): `GWAS_accessions`,
//! `GWAS_associated_gene`, `GWAS_beta_coef`, `GWAS_odds_ratio`, `GWAS_p_value`,
//! `GWAS_pmid`, `GWAS_risk_allele` and `GWAS_study`, each a string that may be
//! empty.
//!
//! Ensembl VEP's plugin writes each file once into a processed, sorted and
//! indexed copy under its data directory and reads that copy; this plugin
//! reads the file itself at load and keeps its records in memory, ordered as
//! the processed copy orders them, so the record it reports at a site with
//! several matching records is the one Ensembl VEP reports.

use std::cmp::Ordering;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use indexmap::IndexMap;

use crate::params::{Grammar, PluginParams};
use crate::tabix::allele_match::{dashed, minimal_forms};
use crate::tabix::{lossy_text, BoundedLines};
use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::codon::reverse_complement;
use vep_core::consequence::{RowKind, TranscriptConsequence};
use vep_core::coordinate::Strand;
use vep_core::sv_routing::perl_true;
use vep_core::variant::InputVariant;

/// The output fields in name order, the order the headers list them.
const FIELDS: [&str; 8] = [
    "GWAS_accessions",
    "GWAS_associated_gene",
    "GWAS_beta_coef",
    "GWAS_odds_ratio",
    "GWAS_p_value",
    "GWAS_pmid",
    "GWAS_risk_allele",
    "GWAS_study",
];

/// The index of `GWAS_risk_allele` in [`FIELDS`], the field a record is
/// matched to a variant by.
const RISK_ALLELE: usize = 6;

/// The field descriptions (GWAS.pm:158-165), in the order of [`FIELDS`].
const DESCRIPTIONS: [&str; 8] = [
    "URI with trait mapped to a ontology term",
    "Gene(s) reported by author in the NHGRI-EBI GWAS Catalog",
    "Beta co-efficient reported for the variant",
    "Odds ratio reported for the variant",
    "P-value reported for the variant",
    "Pubmed identifier of the paper published",
    "Allele associated with the variant that is strongly associated with the trait",
    "GWAS study id",
];

/// The parameters the plugin reads (GWAS.pm:96-104).
const PARAMETERS: [&str; 4] = ["file", "type", "verbose", "db"];

/// The columns of the catalog download the plugin reads (GWAS.pm:257-269),
/// as the upper-cased header names them (`:382`). A header may lack any
/// column outside [`REQUIRED_CATALOG_COLUMNS`]: the plugin reads the columns
/// the header names (`:388`) and an absent one as an empty value
/// (`:396-404`), so the catalog's v1.0 associations download, which has no
/// `STUDY ACCESSION` or `MAPPED_TRAIT_URI`, loads with `GWAS_study` and
/// `GWAS_accessions` empty.
const CATALOG_COLUMNS: [&str; 12] = [
    "PUBMEDID",
    "STUDY ACCESSION",
    "DISEASE/TRAIT",
    "REPORTED GENE(S)",
    "STRONGEST SNP-RISK ALLELE",
    "SNPS",
    "P-VALUE",
    "OR OR BETA",
    "95% CI (TEXT)",
    "MAPPED_TRAIT_URI",
    "CHR_ID",
    "CHR_POS",
];

/// The columns a row's site and identity come from (GWAS.pm:233-237, 272,
/// 330, 395): a header lacking one gives Ensembl VEP's plugin no record to
/// write, and is refused.
const REQUIRED_CATALOG_COLUMNS: [&str; 5] = [
    "DISEASE/TRAIT",
    "STRONGEST SNP-RISK ALLELE",
    "SNPS",
    "CHR_ID",
    "CHR_POS",
];

/// The columns of a harmonised summary-statistics file the plugin reads
/// (GWAS.pm:441), `hm_chrom` and `hm_pos` the two it cannot do without.
const SUMMARY_COLUMNS: [&str; 6] = [
    "hm_chrom",
    "hm_pos",
    "p_value",
    "hm_beta",
    "hm_odds_ratio",
    "hm_effect_allele",
];

/// One association at one site. For a curated record, `text` is the line of
/// Ensembl VEP's processed file after its four coordinate columns: the gene,
/// risk allele, p-value, study, PubMed id, trait URIs, beta coefficient and
/// odds ratio, tab-separated (GWAS.pm:406-416). For a summary-statistics record
/// it is the file's own line, read through the header's column map
/// (`:494-501`).
#[derive(Debug)]
struct Record {
    pos: u64,
    text: Box<str>,
}

/// The columns of a harmonised summary-statistics file, by index, and the
/// study the file name carries (GWAS.pm:139-142, 489-492).
#[derive(Debug)]
struct SummaryStatistics {
    p_value: Option<usize>,
    beta: Option<usize>,
    odds_ratio: Option<usize>,
    effect_allele: Option<usize>,
    pmid: String,
    study: String,
    accession: String,
}

#[derive(Debug)]
enum Source {
    Curated,
    SummaryStatistics(SummaryStatistics),
}

/// The loaded file: its records per chromosome, each list ordered as the
/// processed file orders them ([`collate`]).
#[derive(Debug)]
struct Table {
    source: Source,
    by_chr: HashMap<String, Vec<Record>>,
}

/// GWAS annotation plugin.
pub struct GwasPlugin {
    table: Option<Table>,
}

impl GwasPlugin {
    pub fn new() -> Self {
        Self { table: None }
    }

    /// Open the file, decompressing gzip detected by its magic bytes (`1f 8b`),
    /// so a compressed file saved without a `.gz` suffix loads as well. Every
    /// member of the file is read, as `zcat` reads a `bgzip` output
    /// (GWAS.pm:367), which writes one member per 64 KB.
    fn open_maybe_gzipped(path: &str) -> Result<Box<dyn BufRead>, PluginError> {
        let file = std::fs::File::open(path)
            .map_err(|e| PluginError::Init(format!("GWAS: failed to open {path}: {e}")))?;
        let mut reader = BufReader::new(file);
        // A file too short for the magic is not gzip, and a read failure
        // surfaces from the header read that follows.
        let gzipped = reader
            .fill_buf()
            .is_ok_and(|head| head.starts_with(&[0x1f, 0x8b]));
        if gzipped {
            let decoder = flate2::read::MultiGzDecoder::new(reader);
            Ok(Box::new(BufReader::new(decoder)))
        } else {
            Ok(Box::new(reader))
        }
    }

    /// The chromosome a file names, as the engine names it: without a `chr`
    /// prefix and with `M` as `MT`; `None` for an empty field.
    fn normalize_chr(raw: &str) -> Option<String> {
        let s = raw.trim();
        let without_chr = match s.get(..3) {
            Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &s[3..],
            _ => s,
        };
        if without_chr.is_empty() {
            None
        } else if without_chr.eq_ignore_ascii_case("m") {
            Some("MT".to_string())
        } else {
            Some(without_chr.to_string())
        }
    }

    /// `REPORTED GENE(S)` as the plugin keeps it (GWAS.pm:260, 274-277): empty
    /// when the column carries a `?`, else without whitespace, with an en dash
    /// as `-`, without any character outside ASCII, and empty for `-` or `NR`.
    fn normalize_gene(raw_gene: &str) -> String {
        if raw_gene.contains('?') {
            return String::new();
        }
        let mut gene = raw_gene
            .replace(char::is_whitespace, "")
            .replace('\u{2013}', "-");
        gene.retain(|c| c.is_ascii());
        if gene == "-" || gene == "NR" {
            return String::new();
        }
        gene
    }

    /// `OR or BETA` and `95% CI (TEXT)` as the beta coefficient and the odds
    /// ratio (GWAS.pm:288-319). The ratio is the column's closing decimal,
    /// `(\d+)?(\.\d+)$`, with `0` ahead of a bare fraction and `0.00` read as
    /// `0`; a column without one carries neither field. The unit is the
    /// confidence-interval text after a leading `[...]`, without parentheses and
    /// with `µ` as `micro`: a unit naming an increase or a decrease makes the
    /// ratio a beta coefficient, `<ratio> <unit>`, any other an odds ratio,
    /// which a `0` empties (`0 || ""`, `:404`).
    fn parse_ratio(raw_ratio: &str, raw_ratio_info: &str) -> (String, String) {
        let Some(mut ratio) = Self::closing_decimal(raw_ratio) else {
            return (String::new(), String::new());
        };
        if ratio == "0.00" {
            ratio = "0".to_string();
        }
        match Self::ratio_unit(raw_ratio_info) {
            Some(unit) if unit.contains("decrease") || unit.contains("increase") => {
                (format!("{ratio} {unit}"), String::new())
            }
            _ if perl_true(&ratio) => (String::new(), ratio),
            _ => (String::new(), String::new()),
        }
    }

    /// The match of `(\d+)?(\.\d+)$` joined, with `0` ahead of a bare fraction.
    fn closing_decimal(raw: &str) -> Option<String> {
        let before_fraction = raw.trim_end_matches(|c: char| c.is_ascii_digit());
        if before_fraction.len() == raw.len() {
            return None;
        }
        let before_dot = before_fraction.strip_suffix('.')?;
        let before_whole = before_dot.trim_end_matches(|c: char| c.is_ascii_digit());
        let whole = if before_whole.len() == before_dot.len() {
            "0"
        } else {
            &before_dot[before_whole.len()..]
        };
        Some(format!("{whole}{}", &raw[before_dot.len()..]))
    }

    /// The unit of `^\s*(\[.+\])?\s*(.+)$` (GWAS.pm:302-306): the text after
    /// the last `]` that closes a bracketed interval and leaves a character,
    /// else the whole text after its leading whitespace; `None` for a blank
    /// column. Whitespace alone after the `]` is the unit's last character,
    /// since `\s*` gives one back for `(.+)` to match.
    fn ratio_unit(raw: &str) -> Option<String> {
        let text = raw.trim_start_matches(perl_space);
        let mut unit = text;
        if text.starts_with('[') {
            for (close, _) in text.rmatch_indices(']') {
                if close < 2 {
                    break;
                }
                let rest = &text[close + 1..];
                if rest.is_empty() {
                    continue;
                }
                let after_space = rest.trim_start_matches(perl_space);
                unit = if after_space.is_empty() {
                    let last = rest.chars().next_back().map_or(0, char::len_utf8);
                    &rest[rest.len() - last..]
                } else {
                    after_space
                };
                break;
            }
        }
        if unit.is_empty() {
            return None;
        }
        Some(unit.replace(['(', ')'], "").replace('\u{b5}', "micro"))
    }

    /// Every `rs<digits>` of the `SNPS` column, in order and with repeats
    /// (GWAS.pm:324-326).
    fn extract_rs_ids(snps: &str) -> Vec<&str> {
        let mut ids = Vec::new();
        let bytes = snps.as_bytes();
        let mut i = 0usize;
        while i + 2 < bytes.len() {
            if bytes[i] == b'r' && bytes[i + 1] == b's' && bytes[i + 2].is_ascii_digit() {
                let mut j = i + 2;
                while j < bytes.len() && bytes[j].is_ascii_digit() {
                    j += 1;
                }
                ids.push(&snps[i..j]);
                i = j;
            } else {
                i += 1;
            }
        }
        ids
    }

    /// The risk allele of `rs_id` in `STRONGEST SNP-RISK ALLELE` (GWAS.pm:343-348):
    /// the second `-`-separated field of the last `;`, `,` or `x`-separated
    /// token naming the id, without surrounding whitespace, which the plugin
    /// keeps from a token before an ` x ` and so never matches. `None` when no
    /// token names the id or the last one has no `-`.
    fn risk_allele_for_rs_id(rs_risk_allele: &str, rs_id: &str) -> Option<String> {
        let token = rs_risk_allele
            .split([';', ',', 'x'])
            .rfind(|token| token.contains(rs_id))?;
        let allele = token.split('-').nth(1)?.trim_matches(perl_space);
        perl_true(allele).then(|| allele.to_string())
    }

    /// The locations of a catalog row (GWAS.pm:230-252): `CHR_ID` and `CHR_POS`
    /// split on `;`, `,` or `x` with their whitespace removed, so a row for a
    /// haplotype or an interaction names several.
    fn locations(chr_id: &str, chr_pos: &str) -> Option<Vec<(String, String)>> {
        let split = |field: &str| -> Vec<String> {
            field
                .split([';', ',', 'x'])
                .map(|part| part.replace(perl_space, ""))
                .collect()
        };
        let chrs = split(chr_id);
        let starts = split(chr_pos);
        (chrs.len() == starts.len()).then(|| chrs.into_iter().zip(starts).collect())
    }

    /// The records of a curated catalog (GWAS.pm:254-356, 373-422), each the
    /// association of one SNP at one location. The row's SNPs and locations
    /// pair in order: Ensembl VEP's plugin gives every location of the row to
    /// every SNP of the row (`:338-340`), which at a haplotype's sites reports
    /// another SNP's risk allele; a row whose counts differ puts every SNP at
    /// its one location, and is skipped when it names several.
    fn load_curated(
        path: &str,
        verbose: bool,
    ) -> Result<HashMap<String, Vec<Record>>, PluginError> {
        let reader = Self::open_maybe_gzipped(path)?;
        let mut lines = BoundedLines::new(reader, Path::new(path));
        let header_line = lines.next().transpose()?.map(lossy_text).ok_or_else(|| {
            PluginError::Init(format!(
                "GWAS: {path} is empty; a catalog opens with its header"
            ))
        })?;
        if !is_catalog_header(&header_line) {
            return Err(PluginError::Init(format!(
                "GWAS: {path} does not open with the catalog header line, whose first column is \
                 DATE ADDED TO CATALOG"
            )));
        }
        let header: Vec<String> = header_line
            .split('\t')
            .map(|name| name.to_ascii_uppercase())
            .collect();
        let column = |name: &str| header.iter().rposition(|h| h == name);
        let missing: Vec<&str> = REQUIRED_CATALOG_COLUMNS
            .iter()
            .copied()
            .filter(|name| column(name).is_none())
            .collect();
        if !missing.is_empty() {
            return Err(PluginError::Init(format!(
                "GWAS: the header of {path} lacks the column(s) {}",
                missing.join(", ")
            )));
        }
        let [pmid_col, study_col, trait_col, gene_col, risk_col, snps_col, pvalue_col, ratio_col, ratio_info_col, uri_col, chr_col, pos_col] =
            CATALOG_COLUMNS.map(column);

        let mut by_chr: HashMap<String, Vec<Record>> = HashMap::new();
        let mut line_number = 1usize;
        let mut skipped = Skipped::default();
        for line in lines {
            line_number += 1;
            let line = lossy_text(line?);
            let cols: Vec<&str> = line.split('\t').collect();
            let field =
                |index: Option<usize>| index.and_then(|i| cols.get(i)).copied().unwrap_or("");
            if field(trait_col).is_empty() {
                if verbose {
                    tracing::warn!(
                        "GWAS: 'DISEASE/TRAIT' entry is empty for '{}'",
                        field(snps_col)
                    );
                }
                continue;
            }
            let rs_ids = Self::extract_rs_ids(field(snps_col));
            if rs_ids.is_empty() {
                if verbose {
                    tracing::warn!(
                        "GWAS: Could not parse any rsIds from string '{}'",
                        field(snps_col)
                    );
                }
                continue;
            }
            let rs_risk_allele = if field(risk_col).contains('?') {
                ""
            } else {
                field(risk_col)
            };
            let Some(locations) = Self::locations(field(chr_col), field(pos_col)) else {
                skipped.note(
                    line_number,
                    path,
                    "its chromosomes and positions differ in count",
                );
                continue;
            };
            let gene = Self::normalize_gene(field(gene_col));
            let (beta, odds) = Self::parse_ratio(field(ratio_col), field(ratio_info_col));
            let accessions = field(uri_col).trim_end_matches(',');
            let payload = |risk_allele: &str| -> Box<str> {
                [
                    gene.as_str(),
                    risk_allele,
                    field(pvalue_col),
                    field(study_col),
                    field(pmid_col),
                    accessions,
                    beta.as_str(),
                    odds.as_str(),
                ]
                .map(|value| if perl_true(value) { value } else { "" })
                .join("\t")
                .into_boxed_str()
            };
            let pairs: Vec<(&str, &(String, String))> = if locations.len() == rs_ids.len() {
                rs_ids.iter().copied().zip(locations.iter()).collect()
            } else if locations.len() == 1 {
                rs_ids.iter().map(|id| (*id, &locations[0])).collect()
            } else {
                skipped.note(
                    line_number,
                    path,
                    "its SNPs and its locations differ in count",
                );
                continue;
            };
            for (rs_id, (chr, start)) in pairs {
                let Some(risk_allele) = Self::risk_allele_for_rs_id(rs_risk_allele, rs_id) else {
                    continue;
                };
                let Some(chr) = Self::normalize_chr(chr) else {
                    continue;
                };
                let Some(pos) = start.parse::<u64>().ok().filter(|pos| *pos > 0) else {
                    skipped.note(line_number, path, "a position is not a number");
                    continue;
                };
                by_chr.entry(chr).or_default().push(Record {
                    pos,
                    text: payload(&risk_allele),
                });
            }
        }
        skipped.summarize(path);
        Ok(by_chr)
    }

    /// The records of a harmonised summary-statistics file (GWAS.pm:434-480,
    /// 489-513): every line after the header whose `hm_chrom` and `hm_pos` are
    /// not `NA`. The header's last column is read like the others, where
    /// `:437-446` matches each name against the `head -n 1` line with its line
    /// end still on and so never maps that column.
    fn load_summary_statistics(
        path: &str,
    ) -> Result<(SummaryStatistics, HashMap<String, Vec<Record>>), PluginError> {
        let (pmid, study, accession) = Self::summary_statistics_name(path)?;
        let reader = Self::open_maybe_gzipped(path)?;
        let mut lines = BoundedLines::new(reader, Path::new(path));
        let header_line = lines.next().transpose()?.map(lossy_text).ok_or_else(|| {
            PluginError::Init(format!(
                "GWAS: {path} is empty; a summary-statistics file opens with its header"
            ))
        })?;
        let header: Vec<&str> = header_line
            .strip_prefix('#')
            .unwrap_or(&header_line)
            .split('\t')
            .collect();
        let column = |name: &str| header.iter().position(|h| *h == name);
        let [chrom, pos, p_value, beta, odds_ratio, effect_allele] = SUMMARY_COLUMNS.map(column);
        let (Some(chrom_col), Some(pos_col)) = (chrom, pos) else {
            return Err(PluginError::Init(format!(
                "GWAS: the header of {path} lacks hm_chrom or hm_pos, the harmonised \
                 coordinates a summary-statistics file is read by"
            )));
        };
        let mut by_chr: HashMap<String, Vec<Record>> = HashMap::new();
        let mut line_number = 1usize;
        let mut skipped = Skipped::default();
        for line in lines {
            line_number += 1;
            let line = lossy_text(line?);
            let cols: Vec<&str> = line.split('\t').collect();
            let field = |i: usize| cols.get(i).copied().unwrap_or("");
            if field(chrom_col) == "NA" || field(pos_col) == "NA" {
                continue;
            }
            let Some(chr) = Self::normalize_chr(field(chrom_col)) else {
                continue;
            };
            let Some(pos) = field(pos_col).parse::<u64>().ok().filter(|pos| *pos > 0) else {
                skipped.note(line_number, path, "hm_pos is not a number");
                continue;
            };
            by_chr.entry(chr).or_default().push(Record {
                pos,
                text: line.into_boxed_str(),
            });
        }
        skipped.summarize(path);
        let columns = SummaryStatistics {
            p_value,
            beta,
            odds_ratio,
            effect_allele,
            pmid,
            study,
            accession,
        };
        Ok((columns, by_chr))
    }

    /// The PubMed id, study accession and trait accession a summary-statistics
    /// file name carries, `<pmid>-<study>-<accession>.<suffix>`
    /// (GWAS.pm:139-142; the catalog names its files so).
    fn summary_statistics_name(path: &str) -> Result<(String, String, String), PluginError> {
        let name = Path::new(path)
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_default();
        Self::split_summary_statistics_name(&name).ok_or_else(|| {
            PluginError::Init(format!(
                "GWAS: the summary-statistics file name '{name}' is not \
                 <pmid>-<study>-<accession>.<suffix>, the catalog's form the study is read from"
            ))
        })
    }

    /// The three parts of `<pmid>-<study>-<accession>.<suffix>`; `None` when
    /// `name` lacks one of them or the suffix.
    fn split_summary_statistics_name(name: &str) -> Option<(String, String, String)> {
        let (pmid, rest) = name.split_once('-')?;
        let (study, rest) = rest.split_once('-')?;
        // `(.+?)\.(.+)`: the accession is at least one character, the suffix too.
        let first = rest.chars().next()?.len_utf8();
        let dot = rest[first..].find('.')? + first;
        let accession = &rest[..dot];
        if pmid.is_empty() || study.is_empty() || rest.len() == dot + 1 {
            return None;
        }
        Some((pmid.to_string(), study.to_string(), accession.to_string()))
    }

    /// The output values of `record` in the order of [`FIELDS`].
    fn values<'a>(source: &'a Source, record: &'a Record) -> [&'a str; 8] {
        match source {
            Source::Curated => {
                let mut columns = record.text.split('\t');
                let mut column = || columns.next().unwrap_or("");
                let gene = column();
                let risk_allele = column();
                let p_value = column();
                let study = column();
                let pmid = column();
                let accessions = column();
                let beta = column();
                let odds_ratio = column();
                [
                    accessions,
                    gene,
                    beta,
                    odds_ratio,
                    p_value,
                    pmid,
                    risk_allele,
                    study,
                ]
            }
            Source::SummaryStatistics(file) => {
                let columns: Vec<&str> = record.text.split('\t').collect();
                let column = |index: Option<usize>| {
                    index
                        .and_then(|i| columns.get(i).copied())
                        .filter(|value| perl_true(value))
                        .unwrap_or("")
                };
                [
                    file.accession.as_str(),
                    "",
                    column(file.beta),
                    column(file.odds_ratio),
                    column(file.p_value),
                    file.pmid.as_str(),
                    column(file.effect_allele),
                    file.study.as_str(),
                ]
            }
        }
    }

    /// The record Ensembl VEP's plugin reports for `variant` (GWAS.pm:170-199):
    /// the first, in the processed file's order, of the records within two
    /// bases ahead of the variant's start and up to its end whose risk allele
    /// matches an alternate allele of the variant under
    /// `get_matched_variant_alleles`, with the variant's own reference as the
    /// file record's. A record given on the reverse strand compares its
    /// reverse-complemented alleles against the file's (`Utils/Sequence.pm:1182,
    /// 1196`), the reference it lends the file record as given. An allele of a
    /// `--minimal` split record is matched over its own trimmed span (the
    /// documented class `minimal_plugin_sees_the_rejoined_record`).
    fn matching_record<'a>(table: &'a Table, variant: &InputVariant) -> Option<&'a Record> {
        let records = table.by_chr.get(&variant.chr)?;
        let window_start = variant.start.saturating_sub(2);
        let first = records.partition_point(|record| record.pos < window_start);
        let candidates = records[first..]
            .iter()
            .take_while(|record| record.pos <= variant.end);
        // The record's alleles are built once a candidate exists: most records
        // of a whole-chromosome input sit at no catalog site.
        let mut alleles = None;
        for record in candidates {
            let alleles = alleles.get_or_insert_with(|| RecordAlleles::of(variant));
            let risk_allele = Self::values(&table.source, record)[RISK_ALLELE];
            if alleles.matches(risk_allele, record.pos) {
                return Some(record);
            }
        }
        None
    }
}

/// The alleles of the record a row belongs to, as Ensembl VEP's plugin hands
/// them to `get_matched_variant_alleles` (GWAS.pm:182-188): the record's
/// reference and every alternate on the strand the record was given, and
/// their minimal forms at the record's start on the forward strand. An allele
/// of a `--minimal` split record is its own trimmed reference and alternate
/// (the documented class `minimal_plugin_sees_the_rejoined_record`).
struct RecordAlleles {
    reference: Vec<u8>,
    forms: Vec<(Vec<u8>, Vec<u8>, i64)>,
}

impl RecordAlleles {
    fn of(variant: &InputVariant) -> Self {
        let split_allele = variant.original_start.is_some();
        let (reference, alts): (Vec<u8>, Vec<Vec<u8>>) = match &variant.record_allele_string_multi {
            Some(record) if !split_allele => {
                let mut alleles = record.split('/').map(|a| a.as_bytes().to_vec());
                let reference = alleles.next().unwrap_or_default();
                (reference, alleles.collect())
            }
            _ => (
                variant.ref_allele.clone(),
                vec![variant.alt_allele().to_vec()],
            ),
        };
        let forward = |allele: &[u8]| {
            if variant.strand == Strand::Reverse {
                reverse_complement(allele)
            } else {
                allele.to_vec()
            }
        };
        let forward_reference = forward(&reference);
        let forms = alts
            .iter()
            .flat_map(|alt| minimal_forms(&forward_reference, &forward(alt), variant.start as i64))
            .collect();
        Self { reference, forms }
    }

    /// Whether a file record at `pos` with `risk_allele` matches: the plugin
    /// lends the record the variant's reference (GWAS.pm:180, 190) and
    /// compares the minimal forms of that pair with the record's.
    fn matches(&self, risk_allele: &str, pos: u64) -> bool {
        minimal_forms(&self.reference, &dashed(risk_allele.as_bytes()), pos as i64)
            .iter()
            .any(|form| self.forms.contains(form))
    }
}

/// The rows a load skipped for a reason Ensembl VEP's plugin has no path for:
/// the first is named with its line and reason, the total is reported after.
#[derive(Default)]
struct Skipped {
    count: u64,
}

impl Skipped {
    fn note(&mut self, line_number: usize, path: &str, reason: &str) {
        self.count += 1;
        if self.count == 1 {
            tracing::warn!("GWAS: skipping line {line_number} of {path}: {reason}");
        }
    }

    fn summarize(&self, path: &str) {
        if self.count > 1 {
            tracing::warn!("GWAS: skipped {} rows of {path}", self.count);
        }
    }
}

/// Whether `line` is the catalog's header, `^DATE\s+ADDED\s+TO\s+CATALOG`
/// (GWAS.pm:381).
fn is_catalog_header(line: &str) -> bool {
    let mut words = line.split(perl_space).filter(|word| !word.is_empty());
    line.starts_with("DATE")
        && ["DATE", "ADDED", "TO"]
            .iter()
            .all(|expected| words.next() == Some(expected))
        && words.next().is_some_and(|word| word.starts_with("CATALOG"))
}

/// Perl's `\s` on a byte string.
fn perl_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r' | '\x0b' | '\x0c')
}

/// Orders two lines as GNU `sort` orders them in the English UTF-8 locale of
/// the Ensembl VEP image, where the plugin sorts its processed file
/// (GWAS.pm:426, 464) and so fixes which of several records at one site it
/// reports: letters and digits compared first, without case and with every
/// other character ignored; then the case of each letter, lower before upper;
/// then every character, one that is not an ASCII letter or digit before one
/// that is, those among themselves by code point. The order the Ensembl VEP
/// reference is written in; a byte order differs on about one site in eight
/// of the catalog.
fn collate(a: &str, b: &str) -> Ordering {
    let alnum = |s: &str| {
        s.chars()
            .filter(char::is_ascii_alphanumeric)
            .map(|c| c.to_ascii_lowercase())
            .collect::<Vec<char>>()
    };
    let case = |s: &str| {
        s.chars()
            .filter(char::is_ascii_alphabetic)
            .map(|c| c.is_ascii_uppercase())
            .collect::<Vec<bool>>()
    };
    let shape = |s: &str| {
        s.chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() {
                    u32::MAX
                } else {
                    c as u32
                }
            })
            .collect::<Vec<u32>>()
    };
    alnum(a)
        .cmp(&alnum(b))
        .then_with(|| case(a).cmp(&case(b)))
        .then_with(|| shape(a).cmp(&shape(b)))
}

impl Default for GwasPlugin {
    fn default() -> Self {
        Self::new()
    }
}

impl BuiltinPlugin for GwasPlugin {
    fn name(&self) -> &str {
        "GWAS"
    }

    /// The Ensembl VEP module declares `['Feature', 'Intergenic']`: every row,
    /// the intergenic row included.
    fn feature_types(&self) -> &[RowKind] {
        &RowKind::ALL
    }

    fn header_info(&self) -> Vec<(String, String)> {
        FIELDS
            .iter()
            .zip(DESCRIPTIONS)
            .map(|(field, description)| (field.to_string(), description.to_string()))
            .collect()
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let params = PluginParams::parse_with(
            "GWAS",
            Grammar::KeyValueOnly { file_key: "file" },
            &PARAMETERS,
            params,
        )?;
        let Some(file) = params.get("file").filter(|file| !file.is_empty()) else {
            return Err(PluginError::Init(
                "GWAS: no file specified; use file=<path>".to_string(),
            ));
        };
        let file_type = params
            .get("type")
            .filter(|t| perl_true(t))
            .unwrap_or("curated");
        if !matches!(file_type, "curated" | "sstate") {
            return Err(PluginError::Init(format!(
                "GWAS: provided type ({file_type}) is not recognized; the types are curated and \
                 sstate"
            )));
        }
        if params.get("db").is_some_and(perl_true) {
            return Err(PluginError::Init(
                "GWAS: db=1 reads variant locations from the Ensembl database, which vep-rs \
                 does not connect to; omit db or pass db=0"
                    .to_string(),
            ));
        }
        let verbose = params.get("verbose").is_some_and(perl_true);
        let (source, mut by_chr) = if file_type == "curated" {
            (Source::Curated, Self::load_curated(file, verbose)?)
        } else {
            let (columns, by_chr) = Self::load_summary_statistics(file)?;
            (Source::SummaryStatistics(columns), by_chr)
        };
        for records in by_chr.values_mut() {
            records.sort_by(|a, b| a.pos.cmp(&b.pos).then_with(|| collate(&a.text, &b.text)));
        }
        self.table = Some(Table { source, by_chr });
        Ok(())
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        let Some(table) = &self.table else {
            return Ok(());
        };
        for variant in variants.iter_mut() {
            // Ensembl VEP hands the plugin sequence variants alone
            // (BaseVepPlugin.pm `variant_feature_types`).
            if variant.is_structural || variant.variant_class.is_structural() {
                continue;
            }
            let Some(record) = Self::matching_record(table, variant) else {
                continue;
            };
            for (field, value) in FIELDS.iter().zip(Self::values(&table.source, record)) {
                variant
                    .plugin_data
                    .insert(field.to_string(), value.to_string());
            }
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

    fn params(tokens: &[&str]) -> Vec<String> {
        tokens.iter().map(|t| t.to_string()).collect()
    }

    fn init_error(params_list: &[&str]) -> String {
        match GwasPlugin::new().init(&params(params_list)) {
            Err(PluginError::Init(message)) => message,
            other => panic!("expected an init error, got {other:?}"),
        }
    }

    /// A made-up minimal catalog: a leading `DATE ADDED TO CATALOG`, a column
    /// the parser does not read, then the columns it reads by name, the ones
    /// `GWAS.pm:257-269` reads.
    const HEADER: &str = "DATE ADDED TO CATALOG\tPUBMEDID\tDISEASE/TRAIT\tCHR_ID\tCHR_POS\tREPORTED GENE(S)\tSTRONGEST SNP-RISK ALLELE\tSNPS\tP-VALUE\tOR or BETA\t95% CI (TEXT)\tSTUDY ACCESSION\tMAPPED_TRAIT_URI";
    const ROW: &str = "2020-01-01\t12345678\tType 2 diabetes\t21\t33000100\tAPP\trs123-A\trs123\t1e-8\t1.35\t\tGCST000001\thttp://www.ebi.ac.uk/efo/EFO_0001360";

    fn catalog(rows: &[&str]) -> String {
        let mut text = HEADER.to_string();
        for row in rows {
            text.push('\n');
            text.push_str(row);
        }
        text.push('\n');
        text
    }

    /// `text` without the columns the header names in `names`, compared
    /// without case as the plugin reads the header.
    fn without_columns(text: &str, names: &[&str]) -> String {
        let header: Vec<&str> = text.lines().next().unwrap().split('\t').collect();
        let kept: Vec<usize> = (0..header.len())
            .filter(|&i| !names.iter().any(|n| n.eq_ignore_ascii_case(header[i])))
            .collect();
        assert_eq!(kept.len(), header.len() - names.len(), "{names:?}");
        text.lines()
            .map(|line| {
                let cols: Vec<&str> = line.split('\t').collect();
                kept.iter().map(|&i| cols[i]).collect::<Vec<_>>().join("\t")
            })
            .map(|line| format!("{line}\n"))
            .collect()
    }

    fn write(dir: &std::path::Path, name: &str, text: &str, gzip: bool) -> String {
        use std::io::Write;
        let path = dir.join(name);
        if gzip {
            let f = std::fs::File::create(&path).unwrap();
            let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
            enc.write_all(text.as_bytes()).unwrap();
            enc.finish().unwrap();
        } else {
            std::fs::write(&path, text).unwrap();
        }
        path.to_str().unwrap().to_string()
    }

    fn loaded(params_list: &[&str]) -> GwasPlugin {
        let mut plugin = GwasPlugin::new();
        plugin.init(&params(params_list)).unwrap();
        plugin
    }

    fn variant(chr: &str, pos: u64, reference: &str, alt: &str) -> InputVariant {
        InputVariant::new(
            chr.into(),
            pos,
            pos + reference.len() as u64 - 1,
            reference.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        )
    }

    fn annotate(plugin: &GwasPlugin, variant: InputVariant) -> IndexMap<String, String> {
        let mut variants = vec![variant];
        plugin.run_batch(&mut variants).unwrap();
        variants.remove(0).plugin_data
    }

    fn value<'a>(data: &'a IndexMap<String, String>, field: &str) -> Option<&'a str> {
        data.get(field).map(String::as_str)
    }

    /// The position and risk allele of each record of `chr`, in the table's order.
    fn risk_alleles<'a>(table: &'a Table, chr: &str) -> Vec<(u64, &'a str)> {
        table.by_chr[chr]
            .iter()
            .map(|record| {
                (
                    record.pos,
                    GwasPlugin::values(&table.source, record)[RISK_ALLELE],
                )
            })
            .collect()
    }

    #[test]
    fn plugin_name_and_headers() {
        let plugin = GwasPlugin::new();
        assert_eq!(plugin.name(), "GWAS");
        let headers = plugin.header_info();
        assert_eq!(headers.len(), 8);
        assert_eq!(
            headers[0],
            (
                "GWAS_accessions".to_string(),
                "URI with trait mapped to a ontology term".to_string()
            )
        );
        assert_eq!(
            headers[6].1,
            "Allele associated with the variant that is strongly associated with the trait"
        );
    }

    /// `params_to_hash` alone reads the parameters (GWAS.pm:94-96): a file given
    /// by position is a file given without its key.
    #[test]
    fn positional_file_is_refused() {
        assert_eq!(
            init_error(&["/data/gwas_catalog.tsv"]),
            "GWAS: no file specified; use file=<path>"
        );
        assert_eq!(init_error(&[]), "GWAS: no file specified; use file=<path>");
        assert_eq!(
            init_error(&["type=curated"]),
            "GWAS: no file specified; use file=<path>"
        );
    }

    /// `type` takes `curated` or `sstate` and nothing else (GWAS.pm:99-101).
    #[test]
    fn unknown_type_is_an_init_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = write(dir.path(), "gwas_catalog.tsv", &catalog(&[ROW]), false);
        assert_eq!(
            init_error(&[&format!("file={path}"), "type=summary"]),
            "GWAS: provided type (summary) is not recognized; the types are curated and sstate"
        );
        assert_eq!(
            init_error(&[&format!("file={path}"), "kind=curated"]),
            "GWAS: unknown parameter 'kind'; GWAS accepts file, type, verbose, db"
        );
        loaded(&[&format!("file={path}"), "type=curated", "verbose=0", "db=0"]);
    }

    /// `db=1` fetches each SNP's locations from the Ensembl database
    /// (GWAS.pm:104, 201-228), which an offline engine has no connection to.
    #[test]
    fn db_option_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let path = write(dir.path(), "gwas_catalog.tsv", &catalog(&[ROW]), false);
        assert_eq!(
            init_error(&[&format!("file={path}"), "db=1"]),
            "GWAS: db=1 reads variant locations from the Ensembl database, which vep-rs does \
             not connect to; omit db or pass db=0"
        );
    }

    /// The risk allele matches an alternate allele of the record as written on
    /// the forward strand (GWAS.pm:182-194); the complement of the risk allele
    /// is another allele.
    #[test]
    fn risk_allele_complement_does_not_match() {
        let dir = tempfile::tempdir().unwrap();
        let path = write(dir.path(), "gwas_catalog.tsv", &catalog(&[ROW]), false);
        let plugin = loaded(&[&format!("file={path}")]);
        let matched = annotate(&plugin, variant("21", 33000100, "G", "A"));
        assert_eq!(value(&matched, "GWAS_risk_allele"), Some("A"));
        assert_eq!(value(&matched, "GWAS_study"), Some("GCST000001"));
        assert!(annotate(&plugin, variant("21", 33000100, "G", "T")).is_empty());
        assert!(annotate(&plugin, variant("21", 33000101, "G", "A")).is_empty());
    }

    /// The rows read for a record lie from two bases before its start to its
    /// end (GWAS.pm:170-199), a window clamped at position 1 for a record at 1
    /// or 2: a record at 1 reports the row at 1 alone and a record at 2 the
    /// row at 2 alone, while a record at 3 reports neither, since both rows
    /// lie in its window and a row matches at its own position
    /// (`Utils/Sequence.pm:1190-1241`).
    #[test]
    fn records_at_the_first_two_positions_of_a_chromosome_match_their_rows() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            ROW.replace("\t33000100\t", "\t1\t")
                .replace("\trs123-A\trs123\t", "\trs2822999-G\trs2822999\t")
                .replace("GCST000001", "GCST004988"),
            ROW.replace("\t33000100\t", "\t2\t")
                .replace("\trs123-A\trs123\t", "\trs2823265-G\trs2823265\t")
                .replace("GCST000001", "GCST90002367"),
        ];
        let path = write(
            dir.path(),
            "gwas_catalog.tsv",
            &catalog(&rows.iter().map(String::as_str).collect::<Vec<_>>()),
            false,
        );
        let plugin = loaded(&[&format!("file={path}")]);
        assert_eq!(
            risk_alleles(plugin.table.as_ref().unwrap(), "21"),
            [(1, "G"), (2, "G")]
        );

        let at_one = annotate(&plugin, variant("21", 1, "C", "G"));
        assert_eq!(value(&at_one, "GWAS_study"), Some("GCST004988"));
        assert_eq!(value(&at_one, "GWAS_risk_allele"), Some("G"));
        let at_two = annotate(&plugin, variant("21", 2, "C", "G"));
        assert_eq!(value(&at_two, "GWAS_study"), Some("GCST90002367"));
        assert_eq!(value(&at_two, "GWAS_risk_allele"), Some("G"));
        assert!(
            annotate(&plugin, variant("21", 3, "C", "G")).is_empty(),
            "the rows at 1 and 2 are in the window of a record at 3"
        );
    }

    /// A record of several alternate alleles matches through any of them, as
    /// the plugin reads every alternate of the record (GWAS.pm:185). A record
    /// given on the reverse strand compares its complemented alleles while it
    /// lends the file its reference as given (`:180, 184`), so its single-base
    /// reference never agrees with its complement and no risk allele matches.
    #[test]
    fn every_alternate_of_the_record_is_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = write(dir.path(), "gwas_catalog.tsv", &catalog(&[ROW]), false);
        let plugin = loaded(&[&format!("file={path}")]);
        let mut row = variant("21", 33000100, "G", "C");
        row.record_allele_string_multi = Some("G/C/A".into());
        assert_eq!(
            value(&annotate(&plugin, row), "GWAS_risk_allele"),
            Some("A")
        );
        for (reference, alt) in [("C", "T"), ("G", "T"), ("G", "A")] {
            let mut reverse = variant("21", 33000100, reference, alt);
            reverse.strand = Strand::Reverse;
            assert!(annotate(&plugin, reverse).is_empty(), "{reference}/{alt}");
        }
    }

    /// An allele of a `--minimal` split record is matched by its own trimmed
    /// alleles over its own span, as the same allele is when it is the record's
    /// only alternate: the insertion allele of `G/A/GT` at the catalog's
    /// `G>A` site carries no association, which its sibling `A` carries. The
    /// record's whole allele list and span would match it through `A`.
    #[test]
    fn a_split_allele_is_matched_by_its_own_alleles() {
        let dir = tempfile::tempdir().unwrap();
        let path = write(dir.path(), "gwas_catalog.tsv", &catalog(&[ROW]), false);
        let plugin = loaded(&[&format!("file={path}")]);
        let mut insertion = InputVariant::new(
            "21".into(),
            33000101,
            33000100,
            b"-".to_vec(),
            b"T".to_vec(),
        );
        insertion.original_start = Some(33000100);
        insertion.original_end = Some(33000100);
        insertion.record_allele_string_multi = Some("G/A/GT".into());
        assert!(annotate(&plugin, insertion).is_empty());
        let mut substitution = variant("21", 33000100, "G", "A");
        substitution.original_start = Some(33000100);
        substitution.original_end = Some(33000100);
        substitution.record_allele_string_multi = Some("G/A/GT".into());
        assert_eq!(
            value(&annotate(&plugin, substitution), "GWAS_risk_allele"),
            Some("A")
        );
    }

    /// A haplotype or interaction row names several SNPs, chromosomes and
    /// positions (GWAS.pm:233-235); each SNP is recorded at its own location
    /// with its own risk allele, and a row of several SNPs at one location puts
    /// each there.
    #[test]
    fn multi_location_row_yields_one_record_per_location() {
        let dir = tempfile::tempdir().unwrap();
        let haplotype = ROW
            .replace("\t21\t33000100\t", "\t21;21\t33000100;33000200\t")
            .replace("\trs123-A\trs123\t", "\trs123-A; rs456-C\trs123; rs456\t");
        let interaction = ROW
            .replace("\t21\t33000100\t", "\t11 x 21\t92975544 x 33000300\t")
            .replace(
                "\trs123-A\trs123\t",
                "\trs10830963-G x rs789-T\trs10830963 x rs789\t",
            )
            .replace("GCST000001", "GCST000002");
        let one_location = ROW
            .replace("\t33000100\t", "\t33000400\t")
            .replace("\trs123-A\trs123\t", "\trs1-G; rs2-T\trs1; rs2\t")
            .replace("GCST000001", "GCST000003");
        let path = write(
            dir.path(),
            "gwas_catalog.tsv",
            &catalog(&[&haplotype, &interaction, &one_location]),
            false,
        );
        let plugin = loaded(&[&format!("file={path}")]);
        let table = plugin.table.as_ref().unwrap();
        assert_eq!(
            risk_alleles(table, "21"),
            [
                (33000100, "A"),
                (33000200, "C"),
                (33000300, "T"),
                (33000400, "G"),
                (33000400, "T")
            ]
        );
        assert_eq!(risk_alleles(table, "11"), [(92975544, "G")]);
        assert!(annotate(&plugin, variant("21", 33000100, "G", "C")).is_empty());
        assert_eq!(
            value(
                &annotate(&plugin, variant("21", 33000200, "G", "C")),
                "GWAS_study"
            ),
            Some("GCST000001")
        );
        assert_eq!(
            value(
                &annotate(&plugin, variant("21", 33000400, "A", "T")),
                "GWAS_study"
            ),
            Some("GCST000003")
        );
    }

    /// A byte that is not UTF-8 (Latin-1 `é`, 0xE9) is a byte of its cell, as
    /// `<$input_FH>` and `split /\t/` (GWAS.pm:375-377, 491-494) read the bytes
    /// of a line: in a column the plugin does not read, in another row, and in
    /// the gene, whose non-ASCII bytes the plugin strips (`:276`), the catalog
    /// annotates as without it; in a chromosome it names another chromosome,
    /// so the row annotates nothing; in a printed cell the module prints the
    /// byte, written here as U+FFFD. A summary-statistics file reads alike.
    #[test]
    fn a_byte_that_is_not_utf8_is_a_byte_of_its_cell() {
        let dir = tempfile::tempdir().unwrap();
        let catalog_path = dir.path().join("gwas_catalog.tsv");
        let row = |date: &[u8], chr: &[u8], gene: &[u8], study: &[u8]| -> Vec<u8> {
            let mut row = date.to_vec();
            row.extend_from_slice(b"\t12345678\tType 2 diabetes\t");
            row.extend_from_slice(chr);
            row.extend_from_slice(b"\t33000100\t");
            row.extend_from_slice(gene);
            row.extend_from_slice(b"\trs123-A\trs123\t1e-8\t1.35\t\t");
            row.extend_from_slice(study);
            row.extend_from_slice(b"\thttp://www.ebi.ac.uk/efo/EFO_0001360\n");
            row
        };
        let found = Some(("APP", "GCST000001"));
        let shapes = [
            (row(b"caf\xE9", b"21", b"APP", b"GCST000001"), found),
            (
                [
                    row(b"2020-01-01", b"22", b"caf\xE9", b"GCST000002"),
                    row(b"2020-01-01", b"21", b"APP", b"GCST000001"),
                ]
                .concat(),
                found,
            ),
            (row(b"2020-01-01", b"21", b"AP\xE9P", b"GCST000001"), found),
            (row(b"2020-01-01", b"2\xE91", b"APP", b"GCST000001"), None),
            (
                row(b"2020-01-01", b"21", b"APP", b"GCST00000\xE91"),
                Some(("APP", "GCST00000\u{FFFD}1")),
            ),
        ];
        for (rows, expected) in shapes {
            let mut text = HEADER.as_bytes().to_vec();
            text.push(b'\n');
            text.extend_from_slice(&rows);
            std::fs::write(&catalog_path, &text).unwrap();
            let plugin = loaded(&[&format!("file={}", catalog_path.display())]);
            let data = annotate(&plugin, variant("21", 33000100, "G", "A"));
            assert_eq!(
                value(&data, "GWAS_associated_gene").zip(value(&data, "GWAS_study")),
                expected,
                "{}",
                String::from_utf8_lossy(&rows)
            );
        }
        let sstate_path = dir.path().join("17463246-GCST000028-EFO_0001360.h.tsv");
        std::fs::write(
            &sstate_path,
            b"hm_variant_id\thm_rsid\thm_chrom\thm_pos\thm_other_allele\thm_effect_allele\thm_beta\thm_odds_ratio\tp_value\n\
              21_33000100_G_A\trs123\t21\t33000100\tG\xE9\tA\tNA\t1.2\xE9\t8.088e-05\n",
        )
        .unwrap();
        let plugin = loaded(&["type=sstate", &format!("file={}", sstate_path.display())]);
        let data = annotate(&plugin, variant("21", 33000100, "G", "A"));
        assert_eq!(value(&data, "GWAS_odds_ratio"), Some("1.2\u{FFFD}"));
        assert_eq!(value(&data, "GWAS_risk_allele"), Some("A"));
    }

    /// A harmonised summary-statistics file is read by its header's `hm_`
    /// columns (GWAS.pm:441, 494-509), rows whose `hm_chrom` or `hm_pos` is `NA`
    /// left out (`:460`), and its PubMed id, study and trait accession come from
    /// the file name (`:139-142`); the plugin reads it compressed or not.
    #[test]
    fn sstate_file_parses_harmonised_columns_and_filename() {
        let dir = tempfile::tempdir().unwrap();
        let text = "hm_variant_id\thm_rsid\thm_chrom\thm_pos\thm_other_allele\thm_effect_allele\thm_beta\thm_odds_ratio\tp_value\n\
                    NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0.5\n\
                    21_33000100_G_A\trs123\t21\t33000100\tG\tA\tNA\t1.2\t8.088e-05\n\
                    21_33000200_C_T\trs456\t21\t33000200\tC\tT\t0\tNA\t0.0001\n";
        for (name, gzip) in [
            ("17463246-GCST000028-EFO_0001360.h.tsv.gz", true),
            ("17463246-GCST000028-EFO_0001360.h.tsv", false),
        ] {
            let path = write(dir.path(), name, text, gzip);
            let plugin = loaded(&["type=sstate", &format!("file={path}")]);
            let matched = annotate(&plugin, variant("21", 33000100, "G", "A"));
            assert_eq!(
                matched
                    .iter()
                    .map(|(k, v)| (k.as_str(), v.as_str()))
                    .collect::<Vec<_>>(),
                [
                    ("GWAS_accessions", "EFO_0001360"),
                    ("GWAS_associated_gene", ""),
                    ("GWAS_beta_coef", "NA"),
                    ("GWAS_odds_ratio", "1.2"),
                    ("GWAS_p_value", "8.088e-05"),
                    ("GWAS_pmid", "17463246"),
                    ("GWAS_risk_allele", "A"),
                    ("GWAS_study", "GCST000028"),
                ]
            );
            // `hm_beta || ""` reads a `0` as no value.
            let zero_beta = annotate(&plugin, variant("21", 33000200, "C", "T"));
            assert_eq!(value(&zero_beta, "GWAS_beta_coef"), Some(""));
            assert_eq!(value(&zero_beta, "GWAS_odds_ratio"), Some("NA"));
            assert!(annotate(&plugin, variant("21", 33000100, "G", "T")).is_empty());
        }
        for name in ["summary.tsv", "1-GCST1-", "1-GCST1-EFO_1"] {
            let unnamed = write(dir.path(), name, text, false);
            assert_eq!(
                init_error(&["type=sstate", &format!("file={unnamed}")]),
                format!(
                    "GWAS: the summary-statistics file name '{name}' is not \
                     <pmid>-<study>-<accession>.<suffix>, the catalog's form the study is read from"
                )
            );
        }
        assert_eq!(
            GwasPlugin::split_summary_statistics_name("1-GCST1-\u{e9}.tsv"),
            Some(("1".to_string(), "GCST1".to_string(), "\u{e9}".to_string()))
        );
        let no_coordinates = write(
            dir.path(),
            "1-GCST1-EFO_1.tsv",
            "variant_id\tp_value\nrs1\t0.5\n",
            false,
        );
        assert_eq!(
            init_error(&["type=sstate", &format!("file={no_coordinates}")]),
            format!(
                "GWAS: the header of {no_coordinates} lacks hm_chrom or hm_pos, the harmonised \
                 coordinates a summary-statistics file is read by"
            )
        );
    }

    /// A summary-statistics row with an empty cell before `hm_pos` is read by
    /// its tab-separated columns like any other, so a file whose rows open with
    /// an empty `hm_variant_id` loads and matches its records whatever their
    /// positions' digit counts. Ensembl VEP's plugin filters and sorts the file
    /// counting blank-separated fields (GWAS.pm:460, :464), so the empty cell
    /// shifts its sort onto the position as text, `14482906` before `9053097`,
    /// and the tabix build of :471 refuses the copy, dropping the plugin.
    #[test]
    fn sstate_rows_with_an_empty_leading_cell_are_read_by_their_tab_columns() {
        let dir = tempfile::tempdir().unwrap();
        let text = "hm_variant_id\thm_rsid\thm_chrom\thm_pos\thm_other_allele\thm_effect_allele\thm_beta\thm_odds_ratio\tp_value\tci_upper\n\
                    \trs552964654\t21\t9053097\tG\tC\tNA\tNA\t0.5\tNA\n\
                    \trs2822690\t21\t14482906\tC\tT\tNA\tNA\t0.641\tNA\n";
        let path = write(
            dir.path(),
            "17463246-GCST000028-EFO_0001360.h.tsv",
            text,
            false,
        );
        let plugin = loaded(&["type=sstate", &format!("file={path}")]);
        let seven_digits = annotate(&plugin, variant("21", 9053097, "G", "C"));
        assert_eq!(value(&seven_digits, "GWAS_p_value"), Some("0.5"));
        assert_eq!(value(&seven_digits, "GWAS_risk_allele"), Some("C"));
        assert_eq!(value(&seven_digits, "GWAS_study"), Some("GCST000028"));
        let eight_digits = annotate(&plugin, variant("21", 14482906, "C", "T"));
        assert_eq!(value(&eight_digits, "GWAS_p_value"), Some("0.641"));
        assert_eq!(value(&eight_digits, "GWAS_risk_allele"), Some("T"));
        // The filled form of the same rows reads the same.
        let filled_dir = tempfile::tempdir().unwrap();
        let filled = write(
            filled_dir.path(),
            "17463246-GCST000028-EFO_0001360.h.tsv",
            &text.replace("\n\trs", "\n21_x\trs"),
            false,
        );
        let filled_plugin = loaded(&["type=sstate", &format!("file={filled}")]);
        assert_eq!(
            annotate(&filled_plugin, variant("21", 9053097, "G", "C")),
            seven_digits
        );
    }

    /// Of several records at one site, the one reported is the first in the
    /// order of Ensembl VEP's processed file, GNU `sort` in the English UTF-8
    /// locale: letters and digits without case first, punctuation ignored, so
    /// `APP-X`, first in byte order, sorts after `APPA` and `APPC`, and a
    /// lower-case letter before its upper-case twin.
    #[test]
    fn first_record_in_the_processed_files_order_is_reported() {
        assert_eq!(collate("a\tb", "ab\ta"), Ordering::Less);
        assert_eq!(collate("ab", "Ab"), Ordering::Less);
        assert_eq!(collate("Ab", "AB"), Ordering::Less);
        assert_eq!(collate("a_b\ty", "a-b\tz"), Ordering::Less);
        assert_eq!(collate("a-b", "a_b"), Ordering::Less);
        assert_eq!(collate("ab", "ab-"), Ordering::Less);
        assert_eq!(collate("a~b", "ab"), Ordering::Less);
        assert_eq!(collate("x01", "x1"), Ordering::Less);
        assert_eq!(collate("http://x", "http:/y"), Ordering::Less);
        assert_eq!(collate("B\tA", "ab\ta"), Ordering::Greater);
        assert_eq!(collate("same", "same"), Ordering::Equal);

        let dir = tempfile::tempdir().unwrap();
        let rows = [
            ROW.replace("\tAPP\t", "\tAPPC\t")
                .replace("GCST000001", "GCST000011"),
            ROW.replace("\tAPP\t", "\tAPP-X\t")
                .replace("GCST000001", "GCST000012"),
            ROW.replace("\tAPP\t", "\tAPPA\t")
                .replace("GCST000001", "GCST000013"),
        ];
        let path = write(
            dir.path(),
            "gwas_catalog.tsv",
            &catalog(&rows.iter().map(String::as_str).collect::<Vec<_>>()),
            false,
        );
        let plugin = loaded(&[&format!("file={path}")]);
        let matched = annotate(&plugin, variant("21", 33000100, "G", "A"));
        assert_eq!(value(&matched, "GWAS_associated_gene"), Some("APPA"));
        assert_eq!(value(&matched, "GWAS_study"), Some("GCST000013"));
    }

    /// The catalog's fields reach the output as the plugin writes them
    /// (GWAS.pm:257-319, 396-404): the trait URIs as written, a beta
    /// coefficient as `<ratio> <unit>`, an odds ratio of `0.00` as no value,
    /// and every empty field present and empty.
    #[test]
    fn fields_follow_the_plugins_reading_of_the_catalog() {
        assert_eq!(
            GwasPlugin::parse_ratio("0.237", "[0.11-0.36] unit decrease"),
            ("0.237 unit decrease".to_string(), String::new())
        );
        assert_eq!(
            GwasPlugin::parse_ratio("1.35", "[1.2-1.5]"),
            (String::new(), "1.35".to_string())
        );
        assert_eq!(
            GwasPlugin::parse_ratio(".35", ""),
            (String::new(), "0.35".to_string())
        );
        assert_eq!(
            GwasPlugin::parse_ratio("0.00", "[0.0-0.0] unit increase"),
            ("0 unit increase".to_string(), String::new())
        );
        assert_eq!(
            GwasPlugin::parse_ratio("0.00", "[0.0-0.0]"),
            (String::new(), String::new())
        );
        assert_eq!(
            GwasPlugin::parse_ratio("12", "odds ratio"),
            (String::new(), String::new())
        );
        assert_eq!(
            GwasPlugin::parse_ratio("0.5", "[0.3-0.7] (\u{b5}mol/l increase)"),
            ("0.5 micromol/l increase".to_string(), String::new())
        );
        assert_eq!(
            GwasPlugin::parse_ratio("1.5", " [1.2-1.9] "),
            (String::new(), "1.5".to_string())
        );
        // `(.+)` keeps one trailing space as the unit, so the ratio is an odds
        // ratio; without the space the bracket text itself is the unit.
        assert_eq!(
            GwasPlugin::parse_ratio("0.5", "[1 increase] "),
            (String::new(), "0.5".to_string())
        );
        assert_eq!(
            GwasPlugin::parse_ratio("0.5", "[1 increase]"),
            ("0.5 [1 increase]".to_string(), String::new())
        );
        assert_eq!(GwasPlugin::normalize_gene("NPPA - NPPB"), "NPPA-NPPB");
        assert_eq!(
            GwasPlugin::normalize_gene("APP \u{2013} CYYR1"),
            "APP-CYYR1"
        );
        assert_eq!(GwasPlugin::normalize_gene("NR"), "");
        assert_eq!(GwasPlugin::normalize_gene("-"), "");
        assert_eq!(GwasPlugin::normalize_gene("APP?"), "");
        assert_eq!(
            GwasPlugin::extract_rs_ids("rs123 x rs456; RS789"),
            ["rs123", "rs456"]
        );
        assert_eq!(
            GwasPlugin::risk_allele_for_rs_id("rs123-A; rs456-C", "rs456"),
            Some("C".to_string())
        );
        assert_eq!(
            GwasPlugin::risk_allele_for_rs_id("rs123-G x rs456-T", "rs123"),
            Some("G".to_string())
        );
        assert_eq!(
            GwasPlugin::risk_allele_for_rs_id("rs123-A; rs123", "rs123"),
            None
        );
        assert_eq!(GwasPlugin::risk_allele_for_rs_id("rs456-A", "rs123"), None);

        let dir = tempfile::tempdir().unwrap();
        let row = ROW
            .replace("\t1.35\t\t", "\t0.237\t[0.11-0.36] unit decrease\t")
            .replace(
                "http://www.ebi.ac.uk/efo/EFO_0001360",
                "http://purl.obolibrary.org/obo/MONDO_0004975, http://www.ebi.ac.uk/efo/EFO_0009268,",
            )
            .replace("\tAPP\t", "\t\t");
        let path = write(dir.path(), "gwas_catalog.tsv", &catalog(&[&row]), false);
        let plugin = loaded(&[&format!("file={path}")]);
        let matched = annotate(&plugin, variant("21", 33000100, "G", "A"));
        assert_eq!(
            matched.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect::<Vec<_>>(),
            [
                (
                    "GWAS_accessions",
                    "http://purl.obolibrary.org/obo/MONDO_0004975, http://www.ebi.ac.uk/efo/EFO_0009268"
                ),
                ("GWAS_associated_gene", ""),
                ("GWAS_beta_coef", "0.237 unit decrease"),
                ("GWAS_odds_ratio", ""),
                ("GWAS_p_value", "1e-8"),
                ("GWAS_pmid", "12345678"),
                ("GWAS_risk_allele", "A"),
                ("GWAS_study", "GCST000001"),
            ]
        );
    }

    /// Rows the plugin skips (GWAS.pm:272, 330, 395): an empty `DISEASE/TRAIT`,
    /// no `rs` id, a `?` anywhere in the risk-allele column; and a row whose
    /// chromosome field is empty names no site.
    #[test]
    fn rows_the_plugin_skips_load_no_record() {
        let dir = tempfile::tempdir().unwrap();
        let rows = [
            ROW.replace("\tType 2 diabetes\t", "\t\t"),
            ROW.replace("\trs123-A\trs123\t", "\tchr21:33000100-A\tchr21:33000100\t"),
            ROW.replace("\trs123-A\t", "\trs123-?\t"),
            ROW.replace("\t21\t33000100\t", "\t\t\t"),
            ROW.replace("\t33000100\t", "\t33000101\t"),
        ];
        let path = write(
            dir.path(),
            "gwas_catalog.tsv",
            &catalog(&rows.iter().map(String::as_str).collect::<Vec<_>>()),
            false,
        );
        let plugin = loaded(&[&format!("file={path}")]);
        let table = plugin.table.as_ref().unwrap();
        assert_eq!(table.by_chr.len(), 1);
        assert_eq!(table.by_chr["21"].len(), 1);
        assert_eq!(table.by_chr["21"][0].pos, 33000101);
    }

    /// The catalog opens with its header (GWAS.pm:381-385) and carries every
    /// column a row's site and identity come from.
    #[test]
    fn catalog_without_the_header_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let path = write(
            dir.path(),
            "gwas_catalog.tsv",
            &format!("# catalog\n{}", catalog(&[ROW])),
            false,
        );
        assert_eq!(
            init_error(&[&format!("file={path}")]),
            format!(
                "GWAS: {path} does not open with the catalog header line, whose first column \
                 is DATE ADDED TO CATALOG"
            )
        );
        let path = write(
            dir.path(),
            "narrow.tsv",
            &catalog(&[ROW]).replace("\tCHR_POS\t", "\tPOSITION\t"),
            false,
        );
        assert_eq!(
            init_error(&[&format!("file={path}")]),
            format!("GWAS: the header of {path} lacks the column(s) CHR_POS")
        );
        for name in REQUIRED_CATALOG_COLUMNS {
            let path = write(
                dir.path(),
                "required.tsv",
                &without_columns(&catalog(&[ROW]), &[name]),
                false,
            );
            assert_eq!(
                init_error(&[&format!("file={path}")]),
                format!("GWAS: the header of {path} lacks the column(s) {name}")
            );
        }
        let path = write(dir.path(), "empty.tsv", "", false);
        assert_eq!(
            init_error(&[&format!("file={path}")]),
            format!("GWAS: {path} is empty; a catalog opens with its header")
        );
    }

    /// A header lacking a column outside the five a row's site and identity
    /// come from loads, that field empty on every record: Ensembl VEP's plugin
    /// maps the columns the header names (GWAS.pm:388) and reads an absent one
    /// as `""` (`:396-404`), the odds ratio taking a ratio whose unit column
    /// is absent (`:302`, `:315`). The catalog's v1.0 associations download,
    /// without `STUDY ACCESSION` and `MAPPED_TRAIT_URI`, so annotates with
    /// `GWAS_study` and `GWAS_accessions` empty.
    #[test]
    fn header_without_a_value_column_loads_that_field_empty() {
        let dir = tempfile::tempdir().unwrap();
        let fields = |catalog_text: &str| -> Vec<(String, String)> {
            let path = write(dir.path(), "gwas_catalog.tsv", catalog_text, false);
            let plugin = loaded(&[&format!("file={path}")]);
            annotate(&plugin, variant("21", 33000100, "G", "A"))
                .into_iter()
                .collect()
        };
        let pairs = |values: [&str; 8]| -> Vec<(String, String)> {
            FIELDS
                .iter()
                .zip(values)
                .map(|(field, value)| (field.to_string(), value.to_string()))
                .collect()
        };
        let v1_0_catalog =
            without_columns(&catalog(&[ROW]), &["STUDY ACCESSION", "MAPPED_TRAIT_URI"]);
        assert_eq!(
            fields(&v1_0_catalog),
            pairs(["", "APP", "", "1.35", "1e-8", "12345678", "A", ""])
        );
        let beta_row = ROW.replace("\t1.35\t\t", "\t0.237\t[0.11-0.36] unit decrease\t");
        for (column, expected) in [
            (
                "PUBMEDID",
                [
                    "http://www.ebi.ac.uk/efo/EFO_0001360",
                    "APP",
                    "0.237 unit decrease",
                    "",
                    "1e-8",
                    "",
                    "A",
                    "GCST000001",
                ],
            ),
            (
                "REPORTED GENE(S)",
                [
                    "http://www.ebi.ac.uk/efo/EFO_0001360",
                    "",
                    "0.237 unit decrease",
                    "",
                    "1e-8",
                    "12345678",
                    "A",
                    "GCST000001",
                ],
            ),
            (
                "P-VALUE",
                [
                    "http://www.ebi.ac.uk/efo/EFO_0001360",
                    "APP",
                    "0.237 unit decrease",
                    "",
                    "",
                    "12345678",
                    "A",
                    "GCST000001",
                ],
            ),
            (
                "OR OR BETA",
                [
                    "http://www.ebi.ac.uk/efo/EFO_0001360",
                    "APP",
                    "",
                    "",
                    "1e-8",
                    "12345678",
                    "A",
                    "GCST000001",
                ],
            ),
            (
                "95% CI (TEXT)",
                [
                    "http://www.ebi.ac.uk/efo/EFO_0001360",
                    "APP",
                    "",
                    "0.237",
                    "1e-8",
                    "12345678",
                    "A",
                    "GCST000001",
                ],
            ),
            (
                "MAPPED_TRAIT_URI",
                [
                    "",
                    "APP",
                    "0.237 unit decrease",
                    "",
                    "1e-8",
                    "12345678",
                    "A",
                    "GCST000001",
                ],
            ),
            (
                "STUDY ACCESSION",
                [
                    "http://www.ebi.ac.uk/efo/EFO_0001360",
                    "APP",
                    "0.237 unit decrease",
                    "",
                    "1e-8",
                    "12345678",
                    "A",
                    "",
                ],
            ),
        ] {
            assert_eq!(
                fields(&without_columns(&catalog(&[&beta_row]), &[column])),
                pairs(expected),
                "{column}"
            );
        }
        let value_columns: Vec<&str> = CATALOG_COLUMNS
            .iter()
            .copied()
            .filter(|name| !REQUIRED_CATALOG_COLUMNS.contains(name))
            .collect();
        assert_eq!(
            fields(&without_columns(&catalog(&[&beta_row]), &value_columns)),
            pairs(["", "", "", "", "", "", "A", ""])
        );
    }

    /// The catalog is published compressed; it loads through gzip detected by
    /// content, so a compressed file without a `.gz` suffix loads too, and a
    /// file of several gzip members, as `bgzip` writes one, loads whole.
    #[test]
    fn loads_a_gzipped_catalog() {
        let dir = tempfile::tempdir().unwrap();
        for name in ["gwas_catalog.tsv.gz", "gwas_catalog_no_suffix.tsv"] {
            let path = write(dir.path(), name, &catalog(&[ROW]), true);
            let plugin = loaded(&[&format!("file={path}")]);
            assert_eq!(plugin.table.as_ref().unwrap().by_chr["21"].len(), 1);
        }
        let second_row = ROW.replace("\t33000100\t", "\t33000200\t");
        let first_member =
            std::fs::read(write(dir.path(), "a.gz", &catalog(&[ROW]), true)).unwrap();
        let second_member =
            std::fs::read(write(dir.path(), "b.gz", &format!("{second_row}\n"), true)).unwrap();
        let path = dir.path().join("two_members.tsv.gz");
        std::fs::write(&path, [first_member, second_member].concat()).unwrap();
        let plugin = loaded(&[&format!("file={}", path.display())]);
        assert_eq!(
            risk_alleles(plugin.table.as_ref().unwrap(), "21"),
            [(33000100, "A"), (33000200, "A")]
        );
    }

    /// A catalog whose second row is one byte over the cap is refused as the
    /// plugin opens, naming the file and the line (3, after the header and the
    /// first row).
    #[test]
    fn catalog_line_past_the_cap_is_refused() {
        let cap = crate::tabix::MAX_LINE_BYTES;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("gwas_catalog.tsv");
        let mut text = catalog(&[ROW]).into_bytes();
        text.extend(std::iter::repeat_n(b'x', cap + 1));
        text.push(b'\n');
        std::fs::write(&path, text).unwrap();
        assert_eq!(
            init_error(&[&format!("file={}", path.display())]),
            format!(
                "{}: line 3 is longer than the {cap}-byte limit for one line",
                path.display()
            )
        );
    }
}
