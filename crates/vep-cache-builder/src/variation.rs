// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Parse Ensembl variation VCF files and build variation JSON cache.
//!
//! Reads bgzipped VCF files from Ensembl (germline or somatic), extracts
//! variation records, and writes them as JSON files in the directory structure
//! expected by `vep-cli`'s `json_cache.rs::load_all_variations()`:
//!
//! ```text
//! {output_dir}/variations/{chr}/{region_start}-{region_end}.json
//! ```

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;
use tracing::{debug, info};

/// Region size for binning variations (1 Mb), matching the transcript builder.
const REGION_SIZE: u64 = 1_000_000;

/// First 1-based position of the 1 Mb region containing `pos`.
///
/// Extracted so the binning tests exercise this function rather than re-inlining the
/// same arithmetic on literals: a test that recomputes the formula it is checking
/// cannot fail when the formula changes.
fn region_start(pos: u64) -> u64 {
    ((pos - 1) / REGION_SIZE) * REGION_SIZE + 1
}

/// VCF INFO field name to JSON frequency key mapping.
///
/// The Ensembl VCF stores population frequencies as INFO fields (e.g., `AFR_AF`),
/// while the JSON cache stores them as top-level keys (e.g., `AFR`).
const FREQ_FIELD_MAPPING: &[(&str, &str)] = &[
    ("AF", "AF"),
    ("AFR_AF", "AFR"),
    ("AMR_AF", "AMR"),
    ("EAS_AF", "EAS"),
    ("EUR_AF", "EUR"),
    ("SAS_AF", "SAS"),
    ("AA_AF", "AA"),
    ("EA_AF", "EA"),
    ("gnomADe_AF", "gnomADe"),
    ("gnomADe_AFR_AF", "gnomADe_AFR"),
    ("gnomADe_AMR_AF", "gnomADe_AMR"),
    ("gnomADe_ASJ_AF", "gnomADe_ASJ"),
    ("gnomADe_EAS_AF", "gnomADe_EAS"),
    ("gnomADe_FIN_AF", "gnomADe_FIN"),
    ("gnomADe_MID_AF", "gnomADe_MID"),
    ("gnomADe_NFE_AF", "gnomADe_NFE"),
    ("gnomADe_REMAINING_AF", "gnomADe_REMAINING"),
    ("gnomADe_SAS_AF", "gnomADe_SAS"),
    ("gnomADg_AF", "gnomADg"),
    ("gnomADg_AFR_AF", "gnomADg_AFR"),
    ("gnomADg_AMR_AF", "gnomADg_AMR"),
    ("gnomADg_ASJ_AF", "gnomADg_ASJ"),
    ("gnomADg_EAS_AF", "gnomADg_EAS"),
    ("gnomADg_FIN_AF", "gnomADg_FIN"),
    ("gnomADg_MID_AF", "gnomADg_MID"),
    ("gnomADg_NFE_AF", "gnomADg_NFE"),
    ("gnomADg_REMAINING_AF", "gnomADg_REMAINING"),
    ("gnomADg_SAS_AF", "gnomADg_SAS"),
];

/// A variation record parsed from an Ensembl VCF, ready for JSON serialization.
///
/// Field names and types match `json_cache.rs::JsonVariation` so that
/// `load_all_variations()` can deserialize these without changes.
///
/// Frequency fields are serialized as top-level keys (e.g., `"AF"`, `"gnomADe"`)
/// matching the Perl cache JSON format that `JsonVariation` expects.
#[derive(Debug, Clone, Serialize)]
#[allow(non_snake_case)] // Field names match Perl JSON cache format (e.g., "gnomADe_AFR")
pub struct VariationRecord {
    pub variation_name: String,
    pub start: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end: Option<u64>,
    pub allele_string: String,
    pub strand: i8,
    pub failed: u8,
    pub somatic: u8,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub minor_allele: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub minor_allele_freq: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub clin_sig: Option<String>,
    pub phenotype_or_disease: u8,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pubmed: Option<String>,

    #[serde(rename = "AF", skip_serializing_if = "Option::is_none")]
    pub af: Option<String>,
    #[serde(rename = "AFR", skip_serializing_if = "Option::is_none")]
    pub afr: Option<String>,
    #[serde(rename = "AMR", skip_serializing_if = "Option::is_none")]
    pub amr: Option<String>,
    #[serde(rename = "EAS", skip_serializing_if = "Option::is_none")]
    pub eas: Option<String>,
    #[serde(rename = "EUR", skip_serializing_if = "Option::is_none")]
    pub eur: Option<String>,
    #[serde(rename = "SAS", skip_serializing_if = "Option::is_none")]
    pub sas: Option<String>,
    #[serde(rename = "AA", skip_serializing_if = "Option::is_none")]
    pub aa: Option<String>,
    #[serde(rename = "EA", skip_serializing_if = "Option::is_none")]
    pub ea: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_AFR: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_AMR: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_ASJ: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_EAS: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_FIN: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_MID: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_NFE: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_REMAINING: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADe_SAS: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_AFR: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_AMR: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_ASJ: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_EAS: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_FIN: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_MID: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_NFE: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_REMAINING: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gnomADg_SAS: Option<String>,
}

impl Default for VariationRecord {
    fn default() -> Self {
        Self {
            variation_name: String::new(),
            start: 0,
            end: None,
            allele_string: String::new(),
            strand: 1,
            failed: 0,
            somatic: 0,
            minor_allele: None,
            minor_allele_freq: None,
            clin_sig: None,
            phenotype_or_disease: 0,
            pubmed: None,
            af: None,
            afr: None,
            amr: None,
            eas: None,
            eur: None,
            sas: None,
            aa: None,
            ea: None,
            gnomADe: None,
            gnomADe_AFR: None,
            gnomADe_AMR: None,
            gnomADe_ASJ: None,
            gnomADe_EAS: None,
            gnomADe_FIN: None,
            gnomADe_MID: None,
            gnomADe_NFE: None,
            gnomADe_REMAINING: None,
            gnomADe_SAS: None,
            gnomADg: None,
            gnomADg_AFR: None,
            gnomADg_AMR: None,
            gnomADg_ASJ: None,
            gnomADg_EAS: None,
            gnomADg_FIN: None,
            gnomADg_MID: None,
            gnomADg_NFE: None,
            gnomADg_REMAINING: None,
            gnomADg_SAS: None,
        }
    }
}

impl VariationRecord {
    /// Set a frequency field by its JSON key name (e.g., "AF", "gnomADe_AFR").
    fn set_freq(&mut self, key: &str, value: String) {
        match key {
            "AF" => self.af = Some(value),
            "AFR" => self.afr = Some(value),
            "AMR" => self.amr = Some(value),
            "EAS" => self.eas = Some(value),
            "EUR" => self.eur = Some(value),
            "SAS" => self.sas = Some(value),
            "AA" => self.aa = Some(value),
            "EA" => self.ea = Some(value),
            "gnomADe" => self.gnomADe = Some(value),
            "gnomADe_AFR" => self.gnomADe_AFR = Some(value),
            "gnomADe_AMR" => self.gnomADe_AMR = Some(value),
            "gnomADe_ASJ" => self.gnomADe_ASJ = Some(value),
            "gnomADe_EAS" => self.gnomADe_EAS = Some(value),
            "gnomADe_FIN" => self.gnomADe_FIN = Some(value),
            "gnomADe_MID" => self.gnomADe_MID = Some(value),
            "gnomADe_NFE" => self.gnomADe_NFE = Some(value),
            "gnomADe_REMAINING" => self.gnomADe_REMAINING = Some(value),
            "gnomADe_SAS" => self.gnomADe_SAS = Some(value),
            "gnomADg" => self.gnomADg = Some(value),
            "gnomADg_AFR" => self.gnomADg_AFR = Some(value),
            "gnomADg_AMR" => self.gnomADg_AMR = Some(value),
            "gnomADg_ASJ" => self.gnomADg_ASJ = Some(value),
            "gnomADg_EAS" => self.gnomADg_EAS = Some(value),
            "gnomADg_FIN" => self.gnomADg_FIN = Some(value),
            "gnomADg_MID" => self.gnomADg_MID = Some(value),
            "gnomADg_NFE" => self.gnomADg_NFE = Some(value),
            "gnomADg_REMAINING" => self.gnomADg_REMAINING = Some(value),
            "gnomADg_SAS" => self.gnomADg_SAS = Some(value),
            _ => {}
        }
    }
}

/// Statistics from a variation cache build.
pub struct VariationBuildStats {
    pub variation_count: usize,
    pub region_file_count: usize,
}

/// Parse an Ensembl variation VCF file into variation records.
///
/// Supports both plain-text and bgzipped (.vcf.gz) VCF files. Records where
/// FILTER is not PASS, ".", or empty are skipped. The `somatic` flag is applied
/// to all records from the file (Ensembl distributes germline and somatic VCFs
/// separately).
pub fn parse_variation_vcf(path: &Path, somatic: bool) -> Result<Vec<VariationRecord>> {
    use noodles::vcf;

    let file = std::fs::File::open(path)
        .with_context(|| format!("Failed to open VCF: {}", path.display()))?;

    let reader: Box<dyn BufRead + Send> = if path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "bgz")
    {
        Box::new(std::io::BufReader::new(noodles::bgzf::io::Reader::new(
            file,
        )))
    } else {
        Box::new(std::io::BufReader::new(file))
    };

    let mut vcf_reader = vcf::io::Reader::new(reader);
    let header = vcf_reader
        .read_header()
        .with_context(|| format!("Failed to read VCF header: {}", path.display()))?;

    let mut records = Vec::new();
    let mut record = vcf::Record::default();
    let mut skipped = 0u64;

    loop {
        match vcf_reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => match parse_vcf_record(&record, &header, somatic) {
                Ok(Some(var)) => records.push(var),
                Ok(None) => skipped += 1,
                Err(e) => {
                    debug!("Skipping malformed VCF record: {}", e);
                    skipped += 1;
                }
            },
            Err(e) => {
                debug!("Error reading VCF record, skipping: {}", e);
                skipped += 1;
            }
        }
    }

    info!(
        "Parsed {} variations from {} ({} skipped)",
        records.len(),
        path.display(),
        skipped
    );

    Ok(records)
}

/// Parse a single VCF record into a `VariationRecord`, or `None` if it should be skipped.
fn parse_vcf_record(
    record: &noodles::vcf::Record,
    header: &noodles::vcf::Header,
    somatic: bool,
) -> Result<Option<VariationRecord>> {
    let filters = record.filters();
    let filter_str: &str = filters.as_ref();
    if !filter_str.is_empty() && filter_str != "." && filter_str != "PASS" {
        return Ok(None);
    }

    let pos: u64 = record
        .variant_start()
        .context("Missing POS")?
        .context("Invalid POS")?
        .get() as u64;

    let ids = record.ids();
    let ids_str: &str = ids.as_ref();
    let variation_name = if ids_str.is_empty() || ids_str == "." {
        String::new()
    } else {
        ids_str.split(';').next().unwrap_or("").to_string()
    };

    let ref_allele = record.reference_bases().to_string();

    use noodles::vcf::variant::record::AlternateBases;
    let alts: Vec<String> = record
        .alternate_bases()
        .iter()
        .filter_map(|r| r.ok())
        .filter(|alt| *alt != "." && *alt != "*")
        .map(|s| s.to_string())
        .collect();

    if alts.is_empty() {
        return Ok(None);
    }

    let allele_string = if alts.len() == 1 {
        format!("{}/{}", ref_allele, alts[0])
    } else {
        format!("{}/{}", ref_allele, alts.join("/"))
    };

    let end = if ref_allele.len() <= 1 {
        pos
    } else {
        pos + ref_allele.len() as u64 - 1
    };

    let info = record.info();
    let info_raw: &str = info.as_ref();

    let minor_allele = get_info_string_value(record, header, "MA", info_raw);
    let minor_allele_freq =
        get_info_string_value(record, header, "MAF", info_raw).and_then(|s| s.parse::<f64>().ok());

    let clin_sig = get_info_string_value(record, header, "CLIN_SIG", info_raw);
    let phenotype_or_disease = if clin_sig.is_some() { 1 } else { 0 };

    let pubmed = get_info_string_value(record, header, "PUBMED", info_raw);

    let mut var = VariationRecord {
        variation_name,
        start: pos,
        end: Some(end),
        allele_string,
        strand: 1, // Ensembl VCFs are always forward-strand
        failed: 0,
        somatic: if somatic { 1 } else { 0 },
        minor_allele,
        minor_allele_freq,
        clin_sig,
        phenotype_or_disease,
        pubmed,
        ..Default::default()
    };

    for &(vcf_field, json_key) in FREQ_FIELD_MAPPING {
        if let Some(raw_value) = get_info_string_value(record, header, vcf_field, info_raw) {
            let freq_str = format_freq_value(&raw_value, &alts);
            if !freq_str.is_empty() {
                var.set_freq(json_key, freq_str);
            }
        }
    }

    Ok(Some(var))
}

/// Format a VCF frequency INFO value into the Perl cache "ALT:freq" format.
///
/// For single-allelic sites, this produces "ALT:0.123".
/// For multi-allelic sites with comma-separated frequencies, it produces
/// "ALT1:freq1,ALT2:freq2".
fn format_freq_value(raw_value: &str, alts: &[String]) -> String {
    let freqs: Vec<&str> = raw_value.split(',').collect();

    let mut parts = Vec::new();
    for (i, freq_str) in freqs.iter().enumerate() {
        let freq_str = freq_str.trim();
        if freq_str.is_empty() || freq_str == "." {
            continue;
        }
        let allele = if i < alts.len() {
            &alts[i]
        } else if alts.len() == 1 {
            &alts[0]
        } else {
            continue;
        };
        parts.push(format!("{}:{}", allele, freq_str));
    }

    parts.join(",")
}

/// Try to extract a string INFO field value. Uses the noodles typed API first,
/// then falls back to raw string parsing.
fn get_info_string_value(
    record: &noodles::vcf::Record,
    header: &noodles::vcf::Header,
    key: &str,
    info_raw: &str,
) -> Option<String> {
    use noodles::vcf::variant::record::info::field::Value;

    if let Some(Ok(Some(value))) = record.info().get(header, key) {
        return match value {
            Value::String(s) => {
                let s = s.to_string();
                if s.is_empty() || s == "." {
                    None
                } else {
                    Some(s)
                }
            }
            Value::Integer(n) => Some(n.to_string()),
            Value::Float(f) => Some(f.to_string()),
            Value::Array(arr) => {
                use noodles::vcf::variant::record::info::field::value::Array;
                match arr {
                    Array::Integer(iter) => {
                        let vals: Vec<String> = iter
                            .iter()
                            .filter_map(|r| r.ok())
                            .filter_map(|v| v.map(|n| n.to_string()))
                            .collect();
                        if vals.is_empty() {
                            None
                        } else {
                            Some(vals.join(","))
                        }
                    }
                    Array::Float(iter) => {
                        let vals: Vec<String> = iter
                            .iter()
                            .filter_map(|r| r.ok())
                            .filter_map(|v| v.map(|f| f.to_string()))
                            .collect();
                        if vals.is_empty() {
                            None
                        } else {
                            Some(vals.join(","))
                        }
                    }
                    Array::String(iter) => {
                        let vals: Vec<String> = iter
                            .iter()
                            .filter_map(|r| r.ok())
                            .filter_map(|v| v.map(|s| s.to_string()))
                            .collect();
                        if vals.is_empty() {
                            None
                        } else {
                            Some(vals.join(","))
                        }
                    }
                    Array::Character(iter) => {
                        let vals: Vec<String> = iter
                            .iter()
                            .filter_map(|r| r.ok())
                            .filter_map(|v| v.map(|c| c.to_string()))
                            .collect();
                        if vals.is_empty() {
                            None
                        } else {
                            Some(vals.join(","))
                        }
                    }
                }
            }
            Value::Flag => Some("1".to_string()),
            Value::Character(c) => Some(c.to_string()),
        };
    }

    get_raw_info_value(info_raw, key).map(|s| s.to_string())
}

/// Parse a single INFO field value from the raw semicolon-delimited INFO string.
fn get_raw_info_value<'a>(info_raw: &'a str, key: &str) -> Option<&'a str> {
    if info_raw == "." || info_raw.is_empty() {
        return None;
    }
    for field in info_raw.split(';') {
        if let Some((k, v)) = field.split_once('=') {
            if k == key && !v.is_empty() && v != "." {
                return Some(v);
            }
        }
    }
    None
}

/// Build variation JSON cache files from parsed variation records.
///
/// Groups variations into 1Mb regions and writes one JSON file per region,
/// matching the directory structure expected by `json_cache.rs::load_all_variations()`:
/// `{output_dir}/variations/{chr}/{region_start}-{region_end}.json`
pub fn build_variation_cache(
    records: &[VariationRecord],
    output_dir: &Path,
    chr: &str,
) -> Result<VariationBuildStats> {
    let mut by_region: HashMap<u64, Vec<&VariationRecord>> = HashMap::new();
    for rec in records {
        by_region
            .entry(region_start(rec.start))
            .or_default()
            .push(rec);
    }

    let variations_dir = output_dir.join("variations").join(chr);
    std::fs::create_dir_all(&variations_dir)
        .with_context(|| format!("Failed to create {}", variations_dir.display()))?;

    let mut region_file_count = 0usize;

    for (region_start, region_records) in &by_region {
        let region_end = region_start + REGION_SIZE - 1;
        let filename = format!("{}-{}.json", region_start, region_end);
        let filepath = variations_dir.join(&filename);

        let json =
            serde_json::to_string(region_records).context("Failed to serialize variations")?;

        std::fs::write(&filepath, json)
            .with_context(|| format!("Failed to write {}", filepath.display()))?;

        region_file_count += 1;
    }

    info!(
        "chr {}: {} variations in {} regions",
        chr,
        records.len(),
        region_file_count
    );

    Ok(VariationBuildStats {
        variation_count: records.len(),
        region_file_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal VCF as bytes (plain text, not bgzipped) for testing.
    fn make_vcf_bytes(info_headers: &str, records: &str) -> Vec<u8> {
        format!(
            "##fileformat=VCFv4.1\n\
             ##INFO=<ID=MA,Number=1,Type=String,Description=\"Minor allele\">\n\
             ##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n\
             ##INFO=<ID=CLIN_SIG,Number=.,Type=String,Description=\"Clinical significance\">\n\
             ##INFO=<ID=PUBMED,Number=.,Type=String,Description=\"PubMed IDs\">\n\
             ##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
             ##INFO=<ID=AFR_AF,Number=A,Type=Float,Description=\"AFR allele frequency\">\n\
             ##INFO=<ID=EUR_AF,Number=A,Type=Float,Description=\"EUR allele frequency\">\n\
             ##INFO=<ID=gnomADe_AF,Number=A,Type=Float,Description=\"gnomAD exome AF\">\n\
             {info_headers}\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             {records}"
        )
        .into_bytes()
    }

    /// Write a VCF to a temp file and parse it.
    fn parse_vcf_str(info_headers: &str, records: &str, somatic: bool) -> Vec<VariationRecord> {
        let dir = tempfile::tempdir().unwrap();
        let vcf_path = dir.path().join("test.vcf");
        let data = make_vcf_bytes(info_headers, records);
        std::fs::write(&vcf_path, &data).unwrap();
        parse_variation_vcf(&vcf_path, somatic).unwrap()
    }

    #[test]
    fn test_parse_snv() {
        let records = parse_vcf_str(
            "",
            "21\t25000100\trs699\tA\tG\t.\tPASS\tMA=G;MAF=0.2345\n",
            false,
        );
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.variation_name, "rs699");
        assert_eq!(r.start, 25000100);
        assert_eq!(r.end, Some(25000100));
        assert_eq!(r.allele_string, "A/G");
        assert_eq!(r.strand, 1);
        assert_eq!(r.failed, 0);
        assert_eq!(r.somatic, 0);
        assert_eq!(r.minor_allele, Some("G".to_string()));
        assert_eq!(r.minor_allele_freq, Some(0.2345));
    }

    #[test]
    fn test_parse_deletion() {
        let records = parse_vcf_str("", "21\t100\trs1\tATG\tA\t.\t.\t.\n", false);
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.start, 100);
        assert_eq!(r.end, Some(102)); // start + len("ATG") - 1
        assert_eq!(r.allele_string, "ATG/A");
    }

    #[test]
    fn test_parse_multi_allelic() {
        let records = parse_vcf_str("", "21\t200\trs2\tA\tG,T\t.\tPASS\tAF=0.3,0.1\n", false);
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.allele_string, "A/G/T");
        assert_eq!(r.af, Some("G:0.3,T:0.1".to_string()));
    }

    #[test]
    fn test_somatic_flag() {
        let records = parse_vcf_str("", "21\t100\tCOSV1\tA\tG\t.\t.\t.\n", true);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].somatic, 1);
    }

    #[test]
    fn test_filter_skips_non_pass() {
        let records = parse_vcf_str(
            "",
            "21\t100\trs1\tA\tG\t.\tPASS\t.\n\
             21\t200\trs2\tA\tG\t.\tLowQual\t.\n\
             21\t300\trs3\tA\tG\t.\t.\t.\n",
            false,
        );
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].variation_name, "rs1");
        assert_eq!(records[1].variation_name, "rs3");
    }

    #[test]
    fn test_clin_sig_and_pubmed() {
        let records = parse_vcf_str(
            "",
            "21\t100\trs1\tA\tG\t.\tPASS\tCLIN_SIG=benign,likely_benign;PUBMED=12345,67890\n",
            false,
        );
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.clin_sig, Some("benign,likely_benign".to_string()));
        assert_eq!(r.pubmed, Some("12345,67890".to_string()));
        assert_eq!(r.phenotype_or_disease, 1);
    }

    #[test]
    fn test_no_clin_sig_phenotype_zero() {
        let records = parse_vcf_str("", "21\t100\trs1\tA\tG\t.\tPASS\t.\n", false);
        assert_eq!(records[0].phenotype_or_disease, 0);
    }

    #[test]
    fn test_frequency_formatting_single_allele() {
        let result = format_freq_value("0.123", &["G".to_string()]);
        assert_eq!(result, "G:0.123");
    }

    #[test]
    fn test_frequency_formatting_multi_allele() {
        let result = format_freq_value("0.3,0.1", &["G".to_string(), "T".to_string()]);
        assert_eq!(result, "G:0.3,T:0.1");
    }

    #[test]
    fn test_frequency_formatting_missing_dot() {
        let result = format_freq_value(".", &["G".to_string()]);
        assert_eq!(result, "");
    }

    #[test]
    fn test_population_frequencies() {
        let records = parse_vcf_str(
            "",
            "21\t100\trs1\tA\tG\t.\tPASS\tAF=0.5;AFR_AF=0.3;EUR_AF=0.6;gnomADe_AF=0.45\n",
            false,
        );
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.af, Some("G:0.5".to_string()));
        assert_eq!(r.afr, Some("G:0.3".to_string()));
        assert_eq!(r.eur, Some("G:0.6".to_string()));
        assert_eq!(r.gnomADe, Some("G:0.45".to_string()));
    }

    #[test]
    fn test_build_variation_cache() {
        let dir = tempfile::tempdir().unwrap();
        let records = vec![
            VariationRecord {
                variation_name: "rs1".to_string(),
                start: 500_000,
                end: Some(500_000),
                allele_string: "A/G".to_string(),
                strand: 1,
                ..Default::default()
            },
            VariationRecord {
                variation_name: "rs2".to_string(),
                start: 1_500_000,
                end: Some(1_500_000),
                allele_string: "C/T".to_string(),
                strand: 1,
                ..Default::default()
            },
            VariationRecord {
                variation_name: "rs3".to_string(),
                start: 999_999,
                end: Some(999_999),
                allele_string: "G/A".to_string(),
                strand: 1,
                ..Default::default()
            },
        ];

        let stats = build_variation_cache(&records, dir.path(), "21").unwrap();
        assert_eq!(stats.variation_count, 3);
        assert_eq!(stats.region_file_count, 2); // 1-1000000 and 1000001-2000000

        let region1 = dir.path().join("variations/21/1-1000000.json");
        let region2 = dir.path().join("variations/21/1000001-2000000.json");
        assert!(region1.exists());
        assert!(region2.exists());

        let data1: Vec<serde_json::Value> =
            serde_json::from_str(&std::fs::read_to_string(&region1).unwrap()).unwrap();
        assert_eq!(data1.len(), 2);

        let data2: Vec<serde_json::Value> =
            serde_json::from_str(&std::fs::read_to_string(&region2).unwrap()).unwrap();
        assert_eq!(data2.len(), 1);
        assert_eq!(data2[0]["variation_name"], "rs2");
    }

    #[test]
    fn test_region_binning_boundary() {
        // Position 1 -> region 1-1000000
        assert_eq!(region_start(1), 1);
        // Position 1000000 -> region 1-1000000
        assert_eq!(region_start(1_000_000), 1);
        // Position 1000001 -> region 1000001-2000000
        assert_eq!(region_start(1_000_001), 1_000_001);
    }

    #[test]
    fn test_json_roundtrip_compatible_with_loader() {
        // The serialized JSON must deserialize through the serde structure that
        // json_cache.rs uses.
        let rec = VariationRecord {
            variation_name: "rs699".to_string(),
            start: 230710048,
            end: Some(230710048),
            allele_string: "A/G".to_string(),
            strand: 1,
            failed: 0,
            somatic: 0,
            minor_allele: Some("G".to_string()),
            minor_allele_freq: Some(0.2345),
            clin_sig: Some("benign".to_string()),
            phenotype_or_disease: 1,
            pubmed: Some("12345".to_string()),
            af: Some("G:0.5".to_string()),
            afr: Some("G:0.3".to_string()),
            gnomADe: Some("G:0.45".to_string()),
            ..Default::default()
        };

        let json = serde_json::to_string(&rec).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();

        assert_eq!(parsed["variation_name"], "rs699");
        assert_eq!(parsed["start"], 230710048);
        assert_eq!(parsed["end"], 230710048);
        assert_eq!(parsed["allele_string"], "A/G");
        assert_eq!(parsed["strand"], 1);
        assert_eq!(parsed["failed"], 0);
        assert_eq!(parsed["somatic"], 0);
        assert_eq!(parsed["minor_allele"], "G");
        assert_eq!(parsed["minor_allele_freq"], 0.2345);
        assert_eq!(parsed["clin_sig"], "benign");
        assert_eq!(parsed["phenotype_or_disease"], 1);
        assert_eq!(parsed["pubmed"], "12345");
        assert_eq!(parsed["AF"], "G:0.5");
        assert_eq!(parsed["AFR"], "G:0.3");
        assert_eq!(parsed["gnomADe"], "G:0.45");

        assert!(parsed.get("gnomADe_AFR").is_none());
        assert!(parsed.get("gnomADg").is_none());
        assert!(parsed.get("var_synonyms").is_none());
    }

    #[test]
    fn test_empty_variation_name() {
        let records = parse_vcf_str("", "21\t100\t.\tA\tG\t.\tPASS\t.\n", false);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].variation_name, "");
    }

    #[test]
    fn test_multiple_ids_takes_first() {
        let records = parse_vcf_str("", "21\t100\trs1;rs2;rs3\tA\tG\t.\tPASS\t.\n", false);
        assert_eq!(records[0].variation_name, "rs1");
    }

    #[test]
    fn test_monomorphic_site_skipped() {
        let records = parse_vcf_str("", "21\t100\trs1\tA\t.\t.\tPASS\t.\n", false);
        assert!(records.is_empty());
    }

    #[test]
    fn test_raw_info_value() {
        assert_eq!(
            get_raw_info_value("MA=G;MAF=0.5;CLIN_SIG=benign", "MA"),
            Some("G")
        );
        assert_eq!(
            get_raw_info_value("MA=G;MAF=0.5;CLIN_SIG=benign", "MAF"),
            Some("0.5")
        );
        assert_eq!(
            get_raw_info_value("MA=G;MAF=0.5;CLIN_SIG=benign", "PUBMED"),
            None
        );
        assert_eq!(get_raw_info_value(".", "MA"), None);
        assert_eq!(get_raw_info_value("", "MA"), None);
    }

    #[test]
    fn test_build_empty_records() {
        let dir = tempfile::tempdir().unwrap();
        let stats = build_variation_cache(&[], dir.path(), "21").unwrap();
        assert_eq!(stats.variation_count, 0);
        assert_eq!(stats.region_file_count, 0);
    }
}
