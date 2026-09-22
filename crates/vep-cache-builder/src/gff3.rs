// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Parse Ensembl GFF3 files into gene -> transcript -> exon/CDS hierarchy.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use tracing::debug;

/// A parsed gene record with its child transcripts.
#[derive(Debug)]
#[allow(dead_code)]
pub struct GeneRecord {
    pub gene_id: String,
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub biotype: String,
    pub name: Option<String>,
    pub description: Option<String>,
    pub hgnc_id: Option<String>,
    pub transcripts: Vec<TranscriptRecord>,
}

/// A parsed transcript record with its child exons and CDS regions.
#[derive(Debug)]
pub struct TranscriptRecord {
    pub transcript_id: String,
    pub gene_id: String,
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub biotype: String,
    pub source: String,
    pub version: Option<u32>,
    pub is_canonical: bool,
    pub mane_select: Option<String>,
    pub mane_plus_clinical: Option<String>,
    pub tsl: Option<u8>,
    pub appris: Option<String>,
    /// Whether the GFF3 `tag` attribute lists `basic` (GENCODE basic set).
    pub gencode_basic: bool,
    pub exons: Vec<ExonRecord>,
    pub cds_regions: Vec<CdsRecord>,
}

/// A parsed exon record.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ExonRecord {
    pub exon_id: Option<String>,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub rank: Option<u32>,
}

/// A parsed CDS record.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct CdsRecord {
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub phase: i8,
}

/// Parse a GFF3 file (optionally gzipped) into a gene map.
///
/// Returns `HashMap<gene_id, GeneRecord>`.
pub fn parse_gff3(path: &Path) -> Result<HashMap<String, GeneRecord>> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;

    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut genes: HashMap<String, GeneRecord> = HashMap::new();
    let mut transcript_map: HashMap<String, TranscriptRecord> = HashMap::new();
    let mut exons_by_transcript: HashMap<String, Vec<ExonRecord>> = HashMap::new();
    let mut cds_by_transcript: HashMap<String, Vec<CdsRecord>> = HashMap::new();
    let mut transcript_to_gene: HashMap<String, String> = HashMap::new();

    for line in reader.lines() {
        let line = line.context("Failed to read GFF3 line")?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let chr = normalize_chr(fields[0]);
        let feature_type = fields[2];
        let start: u64 = fields[3].parse().unwrap_or(0);
        let end: u64 = fields[4].parse().unwrap_or(0);
        let strand: i8 = match fields[6] {
            "-" => -1,
            _ => 1,
        };
        let phase: i8 = fields[7].parse().unwrap_or(-1);
        let attrs = parse_attributes(fields[8]);

        match feature_type {
            "gene" | "pseudogene" => {
                if let Some(raw_id) = attrs.get("ID") {
                    let gene_id = strip_prefix(raw_id, "gene:");
                    let biotype = attrs.get("biotype").cloned().unwrap_or_default();
                    let name = attrs.get("Name").cloned();
                    let description = attrs.get("description").cloned();
                    let hgnc_id = extract_hgnc_id(&attrs);

                    genes.insert(
                        gene_id.clone(),
                        GeneRecord {
                            gene_id: gene_id.clone(),
                            chr: chr.to_string(),
                            start,
                            end,
                            strand,
                            biotype,
                            name,
                            description,
                            hgnc_id,
                            transcripts: Vec::new(),
                        },
                    );
                }
            }
            "mRNA"
            | "transcript"
            | "lnc_RNA"
            | "miRNA"
            | "ncRNA"
            | "rRNA"
            | "snRNA"
            | "snoRNA"
            | "scRNA"
            | "tRNA"
            | "pseudogenic_transcript"
            | "unconfirmed_transcript"
            | "mRNA_TE_gene"
            | "V_gene_segment"
            | "D_gene_segment"
            | "J_gene_segment"
            | "C_gene_segment"
            | "NMD_transcript_variant"
            | "pre_miRNA" => {
                if let Some(raw_id) = attrs.get("ID") {
                    let transcript_id = strip_prefix(raw_id, "transcript:");
                    let parent = attrs
                        .get("Parent")
                        .map(|p| strip_prefix(p, "gene:"))
                        .unwrap_or_default();
                    let biotype = attrs.get("biotype").cloned().unwrap_or_default();
                    let source = fields[1].to_string();
                    let version = attrs.get("version").and_then(|v| v.parse::<u32>().ok());
                    let tag = attrs.get("tag").cloned().unwrap_or_default();
                    let is_canonical = tag.contains("Ensembl_canonical");
                    let mane_select = if tag.contains("MANE_Select") {
                        extract_mane_id(&attrs, "MANE_Select")
                    } else {
                        None
                    };
                    let mane_plus_clinical = if tag.contains("MANE_Plus_Clinical") {
                        extract_mane_id(&attrs, "MANE_Plus_Clinical")
                    } else {
                        None
                    };
                    let tsl = extract_tsl(&tag);
                    let appris = extract_appris(&tag);
                    let gencode_basic = tag.split(',').any(|t| t == "basic");

                    transcript_to_gene.insert(transcript_id.clone(), parent.clone());
                    transcript_map.insert(
                        transcript_id.clone(),
                        TranscriptRecord {
                            transcript_id: transcript_id.clone(),
                            gene_id: parent,
                            chr: chr.to_string(),
                            start,
                            end,
                            strand,
                            biotype,
                            source,
                            version,
                            is_canonical,
                            mane_select,
                            mane_plus_clinical,
                            tsl,
                            appris,
                            gencode_basic,
                            exons: Vec::new(),
                            cds_regions: Vec::new(),
                        },
                    );
                }
            }
            "exon" => {
                let parent_transcripts = attrs
                    .get("Parent")
                    .map(|p| {
                        p.split(',')
                            .map(|s| strip_prefix(s, "transcript:"))
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default();
                let exon_id = attrs.get("Name").cloned();
                let rank = attrs.get("rank").and_then(|r| r.parse::<u32>().ok());

                let exon = ExonRecord {
                    exon_id,
                    start,
                    end,
                    strand,
                    rank,
                };

                for parent in parent_transcripts {
                    exons_by_transcript
                        .entry(parent)
                        .or_default()
                        .push(exon.clone());
                }
            }
            "CDS" => {
                let parent_transcripts = attrs
                    .get("Parent")
                    .map(|p| {
                        p.split(',')
                            .map(|s| strip_prefix(s, "transcript:"))
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default();

                let cds = CdsRecord {
                    start,
                    end,
                    strand,
                    phase,
                };

                for parent in parent_transcripts {
                    cds_by_transcript
                        .entry(parent)
                        .or_default()
                        .push(cds.clone());
                }
            }
            _ => {
                debug!("Skipping feature type: {}", feature_type);
            }
        }
    }

    for (tid, mut transcript) in transcript_map {
        if let Some(mut exons) = exons_by_transcript.remove(&tid) {
            exons.sort_by(|a, b| match (a.rank, b.rank) {
                (Some(ra), Some(rb)) => ra.cmp(&rb),
                _ => a.start.cmp(&b.start),
            });
            transcript.exons = exons;
        }
        if let Some(mut cds) = cds_by_transcript.remove(&tid) {
            cds.sort_by_key(|c| c.start);
            transcript.cds_regions = cds;
        }

        if let Some(gene) = genes.get_mut(&transcript.gene_id) {
            gene.transcripts.push(transcript);
        }
    }

    Ok(genes)
}

/// Normalize chromosome names: remove "chr" prefix, "chrM" -> "MT".
fn normalize_chr(chr: &str) -> &str {
    let stripped = chr.strip_prefix("chr").unwrap_or(chr);
    if stripped == "M" {
        "MT"
    } else {
        stripped
    }
}

/// Strip a known prefix (e.g., "gene:", "transcript:") from an ID.
fn strip_prefix(s: &str, prefix: &str) -> String {
    s.strip_prefix(prefix).unwrap_or(s).to_string()
}

/// Parse GFF3 attribute column (key=value pairs separated by semicolons).
fn parse_attributes(attrs_str: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for pair in attrs_str.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        if let Some((key, value)) = pair.split_once('=') {
            let decoded = value
                .replace("%20", " ")
                .replace("%2C", ",")
                .replace("%3B", ";")
                .replace("%3D", "=")
                .replace("%26", "&")
                .replace("%25", "%");
            map.insert(key.to_string(), decoded);
        }
    }
    map
}

/// Extract HGNC ID from GFF3 gene attributes.
///
/// Ensembl GFF3 gene records include `Dbxref=GeneID:123,HGNC:HGNC:12345`.
/// The numeric HGNC ID (e.g., "HGNC:12345") is extracted from this attribute.
fn extract_hgnc_id(attrs: &HashMap<String, String>) -> Option<String> {
    let dbxref = attrs.get("Dbxref")?;
    for entry in dbxref.split(',') {
        // Format: "HGNC:HGNC:12345"; the first "HGNC:" is the Dbxref source
        // prefix, the second is part of the identifier.
        if let Some(rest) = entry.strip_prefix("HGNC:") {
            return Some(rest.to_string());
        }
    }
    None
}

/// Extract Transcript Support Level from GFF3 tag attribute.
///
/// Ensembl GFF3 transcript tags include values like `tsl:1`, `tsl:2`, etc.
/// The TSL is a numeric value 1-5 (or NA).
fn extract_tsl(tag: &str) -> Option<u8> {
    for part in tag.split(',') {
        let part = part.trim();
        if let Some(rest) = part.strip_prefix("tsl:") {
            if let Ok(level) = rest.parse::<u8>() {
                return Some(level);
            }
        }
    }
    None
}

/// Extract APPRIS annotation from GFF3 tag attribute.
///
/// Ensembl GFF3 transcript tags include values like `appris_principal_1`,
/// `appris_alternative_1`, etc., returned in VEP short form (e.g., "P1", "A1").
fn extract_appris(tag: &str) -> Option<String> {
    for part in tag.split(',') {
        let part = part.trim();
        if part.starts_with("appris_") {
            if let Some(rest) = part.strip_prefix("appris_principal_") {
                return Some(format!("P{rest}"));
            } else if let Some(rest) = part.strip_prefix("appris_alternative_") {
                return Some(format!("A{rest}"));
            } else if let Some(rest) = part.strip_prefix("appris_candidate_longest_") {
                return Some(format!("CL{rest}"));
            } else if let Some(rest) = part.strip_prefix("appris_candidate_") {
                return Some(format!("C{rest}"));
            }
            return Some(part.to_string());
        }
    }
    None
}

/// The value stored for a MANE-tagged transcript: its GFF3 `Name` attribute.
///
/// Ensembl GFF3 marks MANE membership in `tag` (`MANE_Select`,
/// `MANE_Plus_Clinical`) and carries the RefSeq identifier (e.g.
/// `NM_001005484.2`) elsewhere, so `Name` is what this builder records for
/// both MANE columns.
fn extract_mane_id(attrs: &HashMap<String, String>, _mane_type: &str) -> Option<String> {
    attrs.get("Name").cloned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_chr() {
        assert_eq!(normalize_chr("chr1"), "1");
        assert_eq!(normalize_chr("chrX"), "X");
        assert_eq!(normalize_chr("chrM"), "MT");
        assert_eq!(normalize_chr("1"), "1");
        assert_eq!(normalize_chr("MT"), "MT");
    }

    #[test]
    fn test_strip_prefix() {
        assert_eq!(strip_prefix("gene:ENSG00000123", "gene:"), "ENSG00000123");
        assert_eq!(
            strip_prefix("transcript:ENST00000456", "transcript:"),
            "ENST00000456"
        );
        assert_eq!(strip_prefix("ENSG00000123", "gene:"), "ENSG00000123");
    }

    #[test]
    fn test_parse_attributes() {
        let attrs = parse_attributes("ID=gene:ENSG00000123;Name=BRCA2;biotype=protein_coding");
        assert_eq!(attrs.get("ID").unwrap(), "gene:ENSG00000123");
        assert_eq!(attrs.get("Name").unwrap(), "BRCA2");
        assert_eq!(attrs.get("biotype").unwrap(), "protein_coding");
    }

    #[test]
    fn test_parse_attributes_url_decode() {
        let attrs = parse_attributes("description=Gene%20name%2C%20variant");
        assert_eq!(attrs.get("description").unwrap(), "Gene name, variant");
    }

    #[test]
    fn test_extract_hgnc_id() {
        let attrs = parse_attributes("ID=gene:ENSG00000139618;Dbxref=GeneID:675,HGNC:HGNC:1100");
        assert_eq!(extract_hgnc_id(&attrs), Some("HGNC:1100".to_string()));
    }

    #[test]
    fn test_extract_hgnc_id_no_hgnc() {
        let attrs = parse_attributes("ID=gene:ENSG00000139618;Dbxref=GeneID:675");
        assert_eq!(extract_hgnc_id(&attrs), None);
    }

    #[test]
    fn test_extract_hgnc_id_no_dbxref() {
        let attrs = parse_attributes("ID=gene:ENSG00000139618;Name=BRCA2");
        assert_eq!(extract_hgnc_id(&attrs), None);
    }

    #[test]
    fn test_extract_tsl() {
        assert_eq!(extract_tsl("basic,Ensembl_canonical,tsl:1"), Some(1));
        assert_eq!(extract_tsl("basic,tsl:5"), Some(5));
        assert_eq!(extract_tsl("basic,Ensembl_canonical"), None);
        assert_eq!(extract_tsl("tsl:NA"), None);
    }

    #[test]
    fn test_extract_appris() {
        assert_eq!(
            extract_appris("basic,appris_principal_1,Ensembl_canonical"),
            Some("P1".to_string())
        );
        assert_eq!(
            extract_appris("basic,appris_alternative_2"),
            Some("A2".to_string())
        );
        assert_eq!(
            extract_appris("basic,appris_candidate_longest_1"),
            Some("CL1".to_string())
        );
        assert_eq!(
            extract_appris("basic,appris_candidate_3"),
            Some("C3".to_string())
        );
        assert_eq!(extract_appris("basic,Ensembl_canonical"), None);
    }
}
