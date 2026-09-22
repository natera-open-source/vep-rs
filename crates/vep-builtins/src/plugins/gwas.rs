// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! GWAS plugin: Perl-compatible GWAS catalog semantics.
//!
//! Parses curated GWAS TSV rows and emits plugin fields aligned to Perl VEP's
//! GWAS plugin behavior.
//!
//! Output fields:
//! - `GWAS_accessions`
//! - `GWAS_associated_gene`
//! - `GWAS_beta_coef`
//! - `GWAS_odds_ratio`
//! - `GWAS_p_value`
//! - `GWAS_pmid`
//! - `GWAS_risk_allele`
//! - `GWAS_study`
//!
//! `GWAS.pm` is the plugin's Perl counterpart in Ensembl/VEP_plugins release/115.

use std::collections::{BTreeSet, HashMap};
use std::io::{BufRead, BufReader};

use indexmap::IndexMap;

use crate::traits::{BuiltinPlugin, PluginError};
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

const ACCESSIONS_FIELD: &str = "GWAS_accessions";
const ASSOCIATED_GENE_FIELD: &str = "GWAS_associated_gene";
const BETA_FIELD: &str = "GWAS_beta_coef";
const ODDS_RATIO_FIELD: &str = "GWAS_odds_ratio";
const P_VALUE_FIELD: &str = "GWAS_p_value";
const PMID_FIELD: &str = "GWAS_pmid";
const RISK_ALLELE_FIELD: &str = "GWAS_risk_allele";
const STUDY_FIELD: &str = "GWAS_study";

#[derive(Clone, Debug, Default)]
struct GwasRecord {
    accessions: String,
    associated_gene: String,
    beta_coef: String,
    odds_ratio: String,
    p_value: String,
    pmid: String,
    risk_allele: String,
    study: String,
}

/// GWAS annotation plugin.
pub struct GwasPlugin {
    by_pos: HashMap<(String, u64), Vec<GwasRecord>>,
}

impl GwasPlugin {
    pub fn new() -> Self {
        Self {
            by_pos: HashMap::new(),
        }
    }

    fn find_column(headers: &[String], candidates: &[&str]) -> Option<usize> {
        for c in candidates {
            if let Some(idx) = headers.iter().position(|h| h.eq_ignore_ascii_case(c)) {
                return Some(idx);
            }
        }
        None
    }

    fn normalize_chr(raw: &str) -> String {
        let s = raw.trim();
        if s.is_empty() {
            return String::new();
        }
        let without_chr = if s.len() >= 3 && s[..3].eq_ignore_ascii_case("chr") {
            &s[3..]
        } else {
            s
        };
        if without_chr.eq_ignore_ascii_case("m") {
            "MT".to_string()
        } else {
            without_chr.to_string()
        }
    }

    fn clean_field(value: Option<&str>) -> String {
        value
            .map(|v| v.trim())
            .filter(|v| !v.is_empty() && *v != ".")
            .unwrap_or("")
            .to_string()
    }

    fn extract_risk_allele_symbol(raw: &str) -> Option<String> {
        let cleaned = raw.trim();
        if cleaned.is_empty() {
            return None;
        }
        let allele = cleaned
            .rsplit_once('-')
            .map(|(_, a)| a)
            .unwrap_or(cleaned)
            .trim();
        if allele.is_empty() {
            return None;
        }
        if allele == "?" {
            return None;
        }
        let allele_upper = allele.to_ascii_uppercase();
        if matches!(allele_upper.as_str(), "A" | "C" | "G" | "T") {
            Some(allele_upper)
        } else {
            None
        }
    }

    fn extract_rs_ids(snps: &str) -> Vec<String> {
        let mut ids = BTreeSet::new();
        let bytes = snps.as_bytes();
        let mut i = 0usize;
        while i + 2 < bytes.len() {
            if (bytes[i] == b'r' || bytes[i] == b'R')
                && (bytes[i + 1] == b's' || bytes[i + 1] == b'S')
                && bytes[i + 2].is_ascii_digit()
            {
                let mut j = i + 2;
                while j < bytes.len() && bytes[j].is_ascii_digit() {
                    j += 1;
                }
                ids.insert(snps[i..j].to_ascii_lowercase());
                i = j;
            } else {
                i += 1;
            }
        }
        ids.into_iter().collect()
    }

    fn risk_allele_for_rs_id(rs_risk_allele: &str, rs_id: &str) -> Option<String> {
        let rs_id_lower = rs_id.to_ascii_lowercase();
        for token in rs_risk_allele.split([';', ',', 'x']) {
            let tok = token.trim();
            if tok.is_empty() {
                continue;
            }
            if !tok.to_ascii_lowercase().contains(&rs_id_lower) {
                continue;
            }
            if let Some(allele) = Self::extract_risk_allele_symbol(tok) {
                return Some(allele);
            }
        }
        None
    }

    fn normalize_gene(raw_gene: &str) -> String {
        if raw_gene.contains('?') {
            return String::new();
        }
        let mut gene = raw_gene.replace(char::is_whitespace, "");
        gene = gene.replace('–', "-");
        gene.retain(|c| c.is_ascii());
        if gene == "-" || gene.eq_ignore_ascii_case("NR") {
            return String::new();
        }
        gene
    }

    fn parse_ratio(raw_ratio: &str, raw_ratio_info: &str) -> (String, String) {
        let ratio = raw_ratio.trim();
        if ratio.is_empty() {
            return (String::new(), String::new());
        }

        let mut normalized_ratio = None;
        if let Some((pre, post)) = ratio.split_once('.') {
            if post.chars().all(|c| c.is_ascii_digit())
                && (pre.is_empty() || pre.chars().all(|c| c.is_ascii_digit()))
            {
                let mut value = if pre.is_empty() {
                    format!("0.{post}")
                } else {
                    format!("{pre}.{post}")
                };
                if value == "0.00" {
                    value = "0".to_string();
                }
                normalized_ratio = Some(value);
            }
        }

        let Some(ratio_num) = normalized_ratio else {
            return (String::new(), String::new());
        };

        let mut unit = raw_ratio_info.trim().to_string();
        if unit.starts_with('[') && unit.contains(']') {
            if let Some(idx) = unit.find(']') {
                unit = unit[idx + 1..].trim().to_string();
            }
        }
        unit = unit.replace(['(', ')'], "");
        unit = unit.replace('µ', "micro");

        if !unit.is_empty() && (unit.contains("decrease") || unit.contains("increase")) {
            (format!("{ratio_num} {unit}"), String::new())
        } else {
            (String::new(), ratio_num)
        }
    }

    fn complement_base(allele: &str) -> Option<String> {
        if allele.len() != 1 {
            return None;
        }
        let b = allele.as_bytes()[0].to_ascii_uppercase() as char;
        let comp = match b {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => return None,
        };
        Some(comp.to_string())
    }

    fn allele_matches_variant(risk_allele: &str, variant_alt: &str) -> bool {
        if risk_allele.eq_ignore_ascii_case(variant_alt) {
            return true;
        }
        if let Some(comp) = Self::complement_base(risk_allele) {
            return comp.eq_ignore_ascii_case(variant_alt);
        }
        false
    }

    /// Open the catalog, transparently decompressing gzip.
    ///
    /// The GWAS Catalog is distributed as a `.zip`/`.gz`, and Perl's `GWAS.pm`
    /// pipes the file through `zcat` when its name ends in `gz`. Reading a
    /// compressed file as plain text here would yield zero parseable rows, so
    /// detect gzip by magic bytes (`1f 8b`) rather than by extension: that also
    /// covers a gzipped file saved without a `.gz` suffix.
    fn open_maybe_gzipped(path: &str) -> Result<Box<dyn BufRead>, PluginError> {
        let mut magic = [0u8; 2];
        {
            use std::io::Read;
            let mut probe = std::fs::File::open(path)
                .map_err(|e| PluginError::Init(format!("GWAS: failed to open {path}: {e}")))?;
            // Under 2 bytes cannot be gzip; the header check below reports the error.
            let _ = probe.read(&mut magic);
        }

        let file = std::fs::File::open(path)
            .map_err(|e| PluginError::Init(format!("GWAS: failed to open {path}: {e}")))?;

        if magic == [0x1f, 0x8b] {
            Ok(Box::new(BufReader::new(flate2::read::GzDecoder::new(file))))
        } else {
            Ok(Box::new(BufReader::new(file)))
        }
    }

    fn load_records(path: &str) -> Result<HashMap<(String, u64), Vec<GwasRecord>>, PluginError> {
        let reader = Self::open_maybe_gzipped(path)?;
        let mut lines = reader.lines();

        let mut header_line = None;
        for line in lines.by_ref() {
            let line =
                line.map_err(|e| PluginError::Init(format!("GWAS: failed to read header: {e}")))?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            header_line = Some(line);
            break;
        }

        let header_line = header_line.ok_or_else(|| {
            PluginError::Init("GWAS: file appears empty (no header line)".to_string())
        })?;

        let headers: Vec<String> = header_line
            .split('\t')
            .map(|s| s.trim().to_string())
            .collect();

        let chr_col = Self::find_column(&headers, &["CHR_ID", "CHR", "CHROM", "CHROMOSOME"])
            .ok_or_else(|| PluginError::Init("GWAS: missing chromosome column".to_string()))?;
        let pos_col = Self::find_column(&headers, &["CHR_POS", "POS", "POSITION"])
            .ok_or_else(|| PluginError::Init("GWAS: missing position column".to_string()))?;

        let trait_uri_col = Self::find_column(&headers, &["MAPPED_TRAIT_URI", "MAPPED TRAIT URI"]);
        let study_accession_col = Self::find_column(
            &headers,
            &["STUDY ACCESSION", "STUDY_ACCESSION", "ACCESSION"],
        );
        let disease_trait_col = Self::find_column(&headers, &["DISEASE/TRAIT", "DISEASE_TRAIT"]);
        let associated_gene_col = Self::find_column(
            &headers,
            &[
                "REPORTED GENE(S)",
                "REPORTED_GENE_S",
                "REPORTED_GENE",
                "MAPPED_GENE",
            ],
        );
        let ratio_col = Self::find_column(&headers, &["OR OR BETA", "OR_OR_BETA", "OR or BETA"]);
        let ratio_info_col = Self::find_column(&headers, &["95% CI (TEXT)", "95_CI_TEXT"]);
        let pvalue_col = Self::find_column(&headers, &["P-VALUE", "P_VALUE", "PVALUE"]);
        let pmid_col = Self::find_column(&headers, &["PUBMEDID", "PMID"]);
        let risk_col = Self::find_column(
            &headers,
            &[
                "STRONGEST SNP-RISK ALLELE",
                "STRONGEST_SNP_RISK_ALLELE",
                "RISK_ALLELE",
            ],
        );
        let snps_col = Self::find_column(&headers, &["SNPS", "SNP", "RSID"]);

        let mut by_pos: HashMap<(String, u64), Vec<GwasRecord>> = HashMap::new();
        for line in lines {
            let line =
                line.map_err(|e| PluginError::Init(format!("GWAS: failed to read row: {e}")))?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            let chr = cols
                .get(chr_col)
                .map(|s| Self::normalize_chr(s))
                .unwrap_or_default();
            let pos = cols
                .get(pos_col)
                .and_then(|p| p.trim().parse::<u64>().ok())
                .unwrap_or(0);
            if chr.is_empty() || pos == 0 {
                continue;
            }

            let disease_trait =
                Self::clean_field(disease_trait_col.and_then(|i| cols.get(i).copied()));
            if disease_trait.is_empty() {
                continue;
            }

            let accessions = trait_uri_col
                .and_then(|i| cols.get(i).copied())
                .map(|raw| {
                    raw.split(',')
                        .map(|s| s.trim())
                        .filter(|s| !s.is_empty())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default();
            let associated_gene = associated_gene_col
                .and_then(|i| cols.get(i).copied())
                .map(Self::normalize_gene)
                .unwrap_or_default();
            let p_value = Self::clean_field(pvalue_col.and_then(|i| cols.get(i).copied()));
            let pmid = Self::clean_field(pmid_col.and_then(|i| cols.get(i).copied()));
            let study = Self::clean_field(study_accession_col.and_then(|i| cols.get(i).copied()));
            let ratio = Self::clean_field(ratio_col.and_then(|i| cols.get(i).copied()));
            let ratio_info = Self::clean_field(ratio_info_col.and_then(|i| cols.get(i).copied()));
            let (beta_coef, odds_ratio) = Self::parse_ratio(&ratio, &ratio_info);
            let snps = Self::clean_field(snps_col.and_then(|i| cols.get(i).copied()));
            let rs_risk_allele = Self::clean_field(risk_col.and_then(|i| cols.get(i).copied()));
            let rs_ids = Self::extract_rs_ids(&snps);
            if rs_ids.is_empty() {
                continue;
            }

            for rs_id in rs_ids {
                let Some(risk_allele) = Self::risk_allele_for_rs_id(&rs_risk_allele, &rs_id) else {
                    continue;
                };
                let record = GwasRecord {
                    accessions: accessions.clone(),
                    associated_gene: associated_gene.clone(),
                    beta_coef: beta_coef.clone(),
                    odds_ratio: odds_ratio.clone(),
                    p_value: p_value.clone(),
                    pmid: pmid.clone(),
                    risk_allele,
                    study: study.clone(),
                };
                by_pos.entry((chr.clone(), pos)).or_default().push(record);
            }
        }

        Ok(by_pos)
    }
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

    fn header_info(&self) -> Vec<(String, String)> {
        vec![
            (
                ACCESSIONS_FIELD.to_string(),
                "GWAS catalog accession(s)".to_string(),
            ),
            (
                ASSOCIATED_GENE_FIELD.to_string(),
                "GWAS associated gene(s)".to_string(),
            ),
            (
                BETA_FIELD.to_string(),
                "GWAS beta coefficient(s)".to_string(),
            ),
            (
                ODDS_RATIO_FIELD.to_string(),
                "GWAS odds ratio(s)".to_string(),
            ),
            (P_VALUE_FIELD.to_string(), "GWAS p-value(s)".to_string()),
            (PMID_FIELD.to_string(), "GWAS PubMed id(s)".to_string()),
            (
                RISK_ALLELE_FIELD.to_string(),
                "GWAS risk allele annotation(s)".to_string(),
            ),
            (
                STUDY_FIELD.to_string(),
                "GWAS study description(s)".to_string(),
            ),
        ]
    }

    fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let mut file_path = None;
        for param in params {
            if let Some(path) = param.strip_prefix("file=") {
                file_path = Some(path.to_string());
            } else if file_path.is_none() {
                file_path = Some(param.to_string());
            }
        }

        let file_path = file_path
            .ok_or_else(|| PluginError::Init("GWAS requires a file path argument".to_string()))?;
        self.by_pos = Self::load_records(&file_path)?;
        if self.by_pos.is_empty() {
            return Err(PluginError::Init(
                "GWAS data file loaded but produced no usable records".to_string(),
            ));
        }
        Ok(())
    }

    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        for variant in variants.iter_mut() {
            let key = (variant.chr.clone(), variant.start);
            let Some(records) = self.by_pos.get(&key) else {
                continue;
            };

            let alt = variant.display_allele();
            let mut first_match: Option<&GwasRecord> = None;
            for record in records {
                if record.risk_allele.is_empty() {
                    continue;
                }
                if !Self::allele_matches_variant(&record.risk_allele, &alt) {
                    continue;
                }
                first_match = Some(record);
                break;
            }

            let Some(record) = first_match else {
                continue;
            };

            if !record.accessions.is_empty() {
                variant
                    .plugin_data
                    .insert(ACCESSIONS_FIELD.to_string(), record.accessions.clone());
            }
            if !record.associated_gene.is_empty() {
                variant.plugin_data.insert(
                    ASSOCIATED_GENE_FIELD.to_string(),
                    record.associated_gene.clone(),
                );
            }
            if !record.beta_coef.is_empty() {
                variant
                    .plugin_data
                    .insert(BETA_FIELD.to_string(), record.beta_coef.clone());
            }
            if !record.odds_ratio.is_empty() {
                variant
                    .plugin_data
                    .insert(ODDS_RATIO_FIELD.to_string(), record.odds_ratio.clone());
            }
            if !record.p_value.is_empty() {
                variant
                    .plugin_data
                    .insert(P_VALUE_FIELD.to_string(), record.p_value.clone());
            }
            if !record.pmid.is_empty() {
                variant
                    .plugin_data
                    .insert(PMID_FIELD.to_string(), record.pmid.clone());
            }
            if !record.risk_allele.is_empty() {
                variant
                    .plugin_data
                    .insert(RISK_ALLELE_FIELD.to_string(), record.risk_allele.clone());
            }
            if !record.study.is_empty() {
                variant
                    .plugin_data
                    .insert(STUDY_FIELD.to_string(), record.study.clone());
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

    #[test]
    fn plugin_name_and_headers() {
        let plugin = GwasPlugin::new();
        assert_eq!(plugin.name(), "GWAS");
        let headers = plugin.header_info();
        assert!(headers.iter().any(|(k, _)| k == "GWAS_accessions"));
        assert!(headers.iter().any(|(k, _)| k == "GWAS_pmid"));
    }

    #[test]
    fn risk_allele_extraction() {
        assert_eq!(
            GwasPlugin::extract_risk_allele_symbol("rs123-A"),
            Some("A".to_string())
        );
        assert_eq!(
            GwasPlugin::extract_risk_allele_symbol("C"),
            Some("C".to_string())
        );
        assert_eq!(GwasPlugin::extract_risk_allele_symbol("rs123-?"), None);
        assert_eq!(GwasPlugin::extract_risk_allele_symbol("rs123-DEL"), None);
    }

    #[test]
    fn gene_normalization_perl_style() {
        assert_eq!(GwasPlugin::normalize_gene("NPPA - NPPB"), "NPPA-NPPB");
        assert_eq!(GwasPlugin::normalize_gene("NR"), "");
        assert_eq!(GwasPlugin::normalize_gene("?"), "");
    }

    #[test]
    fn ratio_parsing_matches_expected_semantics() {
        let (beta, odds) = GwasPlugin::parse_ratio("0.237", "[95% CI] unit decrease");
        assert_eq!(beta, "0.237 unit decrease");
        assert_eq!(odds, "");

        let (beta, odds) = GwasPlugin::parse_ratio("8.0", "odds ratio");
        assert_eq!(beta, "");
        assert_eq!(odds, "8.0");
    }

    #[test]
    fn extract_rs_ids_parses_multiple() {
        let ids = GwasPlugin::extract_rs_ids("foo rs123 bar;rs456");
        assert_eq!(ids, vec!["rs123".to_string(), "rs456".to_string()]);
    }

    #[test]
    fn allele_match_supports_complement() {
        assert!(GwasPlugin::allele_matches_variant("A", "T"));
        assert!(!GwasPlugin::allele_matches_variant("A", "C"));
    }

    /// A minimal catalog with the columns the parser resolves by name. The real
    /// NHGRI-EBI file detects its header via `DATE ADDED TO CATALOG` and carries
    /// ~40 columns; these are the ones this plugin reads.
    const CATALOG: &str = "DATE ADDED TO CATALOG\tPUBMEDID\tDISEASE/TRAIT\tCHR_ID\tCHR_POS\tREPORTED GENE(S)\tSTRONGEST SNP-RISK ALLELE\tSNPS\tP-VALUE\tOR or BETA\t95% CI (TEXT)\tSTUDY ACCESSION\tMAPPED_TRAIT_URI
2020-01-01\t12345678\tType 2 diabetes\t21\t33000100\tAPP\trs123-A\trs123\t1e-8\t1.35\t\tGCST000001\thttp://www.ebi.ac.uk/efo/EFO_0001360";

    fn write_catalog(dir: &std::path::Path, name: &str, gzip: bool) -> std::path::PathBuf {
        use std::io::Write;
        let path = dir.join(name);
        if gzip {
            let f = std::fs::File::create(&path).unwrap();
            let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
            enc.write_all(CATALOG.as_bytes()).unwrap();
            enc.finish().unwrap();
        } else {
            std::fs::write(&path, CATALOG).unwrap();
        }
        path
    }

    /// The GWAS Catalog is published compressed and Perl's GWAS.pm pipes it
    /// through `zcat`; a bare `File::open` on a gzipped catalog yields zero
    /// parseable rows and annotates nothing.
    #[test]
    fn loads_a_gzipped_catalog() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_catalog(dir.path(), "gwas_catalog.tsv.gz", true);

        let records =
            GwasPlugin::load_records(path.to_str().unwrap()).expect("a gzipped catalog must load");
        assert!(
            records.contains_key(&("21".to_string(), 33000100)),
            "expected the chr21 record to be indexed; got {} keys",
            records.len()
        );
    }

    /// A plain (uncompressed) catalog loads through the same path.
    #[test]
    fn loads_a_plain_catalog() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_catalog(dir.path(), "gwas_catalog.tsv", false);

        let records =
            GwasPlugin::load_records(path.to_str().unwrap()).expect("a plain catalog must load");
        assert!(records.contains_key(&("21".to_string(), 33000100)));
    }

    /// Gzip is detected by magic bytes, so a compressed catalog saved without a
    /// `.gz` suffix still loads rather than silently parsing as garbage.
    #[test]
    fn detects_gzip_without_a_gz_extension() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_catalog(dir.path(), "gwas_catalog_no_suffix.tsv", true);

        let records = GwasPlugin::load_records(path.to_str().unwrap())
            .expect("gzip must be detected by content, not extension");
        assert!(records.contains_key(&("21".to_string(), 33000100)));
    }
}
