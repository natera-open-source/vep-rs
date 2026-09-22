// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Build the JSON cache directory structure from parsed GFF3/FASTA data.
//!
//! Output:
//! ```text
//! {output_dir}/
//!   transcripts/
//!     1/
//!       1-1000000.json
//!       1000001-2000000.json
//!     ...
//! ```

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;
use tracing::info;

use crate::gff3::GeneRecord;
use crate::mapper;

/// Region size for binning transcripts (1 Mb).
const REGION_SIZE: u64 = 1_000_000;

/// First 1-based position of the 1 Mb region containing `pos`.
///
/// Extracted so the binning test exercises this function rather than re-inlining the
/// same arithmetic on literals: a test that recomputes the formula it is checking
/// cannot fail when the formula changes.
fn region_start(pos: u64) -> u64 {
    ((pos - 1) / REGION_SIZE) * REGION_SIZE + 1
}

/// Statistics from the cache build.
pub struct BuildStats {
    pub transcript_count: usize,
    pub chromosome_count: usize,
    pub region_file_count: usize,
}

// These must serialize to field names matching what json_cache.rs expects.

#[derive(Serialize)]
struct CacheTranscript {
    stable_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    version: Option<u32>,
    #[serde(rename = "dbID", skip_serializing_if = "Option::is_none")]
    db_id: Option<u64>,
    gene_stable_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    chr: Option<String>,
    start: u64,
    end: u64,
    strand: i8,
    biotype: String,
    source: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    description: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    gene_symbol: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    gene_symbol_source: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none", rename = "gene_hgnc_id")]
    hgnc_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    gene_phenotype: Option<u8>,
    #[serde(rename = "is_canonical")]
    canonical: u8,
    #[serde(skip_serializing_if = "Option::is_none")]
    mane_select: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    mane_plus_clinical: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    tsl: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    appris: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    protein: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    exons: Vec<CacheExon>,
    #[serde(skip_serializing_if = "Option::is_none")]
    cdna_coding_start: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    cdna_coding_end: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    translation: Option<CacheTranslation>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    attributes: Vec<CacheAttribute>,
    variation_effect_feature_cache: CacheVEFC,
}

#[derive(Serialize, Clone)]
struct CacheExon {
    #[serde(skip_serializing_if = "Option::is_none")]
    stable_id: Option<String>,
    start: u64,
    end: u64,
    rank: u32,
    phase: i8,
    end_phase: i8,
    strand: i8,
}

#[derive(Serialize)]
struct CacheTranslation {
    stable_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    version: Option<u32>,
    #[serde(rename = "dbID", skip_serializing_if = "Option::is_none")]
    db_id: Option<u64>,
    start: u64,
    end: u64,
}

#[derive(Serialize)]
struct CacheAttribute {
    code: String,
    value: String,
}

#[derive(Serialize)]
struct CacheVEFC {
    codon_table: u8,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    introns: Vec<CacheIntron>,
    #[serde(skip_serializing_if = "Option::is_none")]
    mapper: Option<CacheMapper>,
    #[serde(skip_serializing_if = "Option::is_none")]
    peptide: Option<String>,
    /// Pre-computed CDS nucleotide sequence (concatenated coding exons from reference FASTA).
    /// Required for codon translation, hence missense/synonymous/stop consequence assignment.
    #[serde(skip_serializing_if = "Option::is_none")]
    translateable_seq: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    sorted_exons: Vec<CacheExon>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    protein_features: Vec<CacheProteinFeature>,
    #[serde(skip_serializing_if = "Option::is_none")]
    protein_function_predictions: Option<CachePredictions>,
}

/// Serialized SIFT/PolyPhen prediction matrices for JSON cache output.
#[derive(Serialize)]
struct CachePredictions {
    #[serde(skip_serializing_if = "Option::is_none")]
    sift: Option<CachePredictionMatrix>,
    #[serde(skip_serializing_if = "Option::is_none")]
    polyphen_humdiv: Option<CachePredictionMatrix>,
    #[serde(skip_serializing_if = "Option::is_none")]
    polyphen_humvar: Option<CachePredictionMatrix>,
}

/// A single prediction matrix in JSON-serializable form.
#[derive(Serialize)]
struct CachePredictionMatrix {
    analysis: String,
    peptide_length: usize,
    /// Base64-encoded gzip-compressed binary matrix.
    predictions_data: String,
}

#[derive(Serialize)]
struct CacheIntron {
    start: u64,
    end: u64,
    rank: u32,
    strand: i8,
}

#[derive(Serialize)]
struct CacheMapper {
    start_phase: i8,
    #[serde(skip_serializing_if = "Option::is_none")]
    cdna_coding_start: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    cdna_coding_end: Option<u64>,
    pair_count: usize,
    pairs: Vec<CacheMapperPair>,
}

#[derive(Serialize)]
struct CacheMapperPair {
    from_start: u64,
    from_end: u64,
    from_id: &'static str,
    to_start: u64,
    to_end: u64,
    to_id: &'static str,
    ori: i8,
}

#[derive(Serialize)]
struct CacheProteinFeature {
    start: u64,
    end: u64,
    hseqname: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    analysis: Option<String>,
}

/// Build the full JSON cache directory from parsed GFF3 genes and protein sequences.
pub fn build_cache(
    genes: &HashMap<String, GeneRecord>,
    protein_sequences: &HashMap<String, String>,
    output_dir: &Path,
    prediction_data: Option<&HashMap<String, crate::predictions::PeptidePredictions>>,
    genome_fasta: Option<&vep_fasta::IndexedFasta>,
    gtf_tags: Option<&crate::gtf::GtfTags>,
) -> Result<BuildStats> {
    let mut by_chr: HashMap<String, Vec<CacheTranscript>> = HashMap::new();
    let mut transcript_count = 0usize;

    for gene in genes.values() {
        for tr in &gene.transcripts {
            let cache_tr = build_transcript(
                gene,
                tr,
                protein_sequences,
                prediction_data,
                genome_fasta,
                gtf_tags,
            );
            by_chr.entry(tr.chr.clone()).or_default().push(cache_tr);
            transcript_count += 1;
        }
    }

    let transcripts_dir = output_dir.join("transcripts");
    let mut region_file_count = 0usize;

    for (chr, transcripts) in &by_chr {
        let mut by_region: HashMap<u64, Vec<&CacheTranscript>> = HashMap::new();
        for tr in transcripts {
            by_region
                .entry(region_start(tr.start))
                .or_default()
                .push(tr);
        }

        let chr_dir = transcripts_dir.join(chr);
        std::fs::create_dir_all(&chr_dir)
            .with_context(|| format!("Failed to create {}", chr_dir.display()))?;

        for (region_start, region_transcripts) in &by_region {
            let region_end = region_start + REGION_SIZE - 1;
            let filename = format!("{}-{}.json", region_start, region_end);
            let filepath = chr_dir.join(&filename);

            let json = serde_json::to_string_pretty(region_transcripts)
                .context("Failed to serialize transcripts")?;

            std::fs::write(&filepath, json)
                .with_context(|| format!("Failed to write {}", filepath.display()))?;

            region_file_count += 1;
        }

        info!(
            "chr {}: {} transcripts in {} regions",
            chr,
            transcripts.len(),
            by_region.len()
        );
    }

    Ok(BuildStats {
        transcript_count,
        chromosome_count: by_chr.len(),
        region_file_count,
    })
}

/// Generate `info.json` in the output directory.
///
/// This file describes the cache contents and is read by `vep-cli` to
/// determine available annotations and cache version.
pub fn generate_info_json(
    output_dir: &Path,
    species: &str,
    assembly: &str,
    cache_version: u32,
    include_predictions: bool,
    source_versions: HashMap<String, String>,
) -> Result<()> {
    let info = CacheInfoJson {
        species: species.to_string(),
        assembly: assembly.to_string(),
        cache_version,
        sift: if include_predictions {
            Some("b".to_string())
        } else {
            None
        },
        polyphen: if include_predictions {
            Some("b".to_string())
        } else {
            None
        },
        regulatory: false,
        variation_cols: vec![
            "variation_name".to_string(),
            "failed".to_string(),
            "somatic".to_string(),
            "start".to_string(),
            "end".to_string(),
            "allele_string".to_string(),
            "strand".to_string(),
            "minor_allele".to_string(),
            "minor_allele_freq".to_string(),
            "clin_sig".to_string(),
            "phenotype_or_disease".to_string(),
            "pubmed".to_string(),
        ],
        source_versions,
    };

    let filepath = output_dir.join("info.json");
    let json = serde_json::to_string_pretty(&info).context("Failed to serialize info.json")?;
    std::fs::write(&filepath, json)
        .with_context(|| format!("Failed to write {}", filepath.display()))?;

    info!("Wrote {}", filepath.display());
    Ok(())
}

/// JSON structure for the cache info.json file.
#[derive(Serialize)]
struct CacheInfoJson {
    species: String,
    assembly: String,
    cache_version: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    sift: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    polyphen: Option<String>,
    regulatory: bool,
    variation_cols: Vec<String>,
    source_versions: HashMap<String, String>,
}

/// Build a single CacheTranscript from parsed GFF3 data.
fn build_transcript(
    gene: &GeneRecord,
    tr: &crate::gff3::TranscriptRecord,
    protein_sequences: &HashMap<String, String>,
    prediction_data: Option<&HashMap<String, crate::predictions::PeptidePredictions>>,
    genome_fasta: Option<&vep_fasta::IndexedFasta>,
    gtf_tags: Option<&crate::gtf::GtfTags>,
) -> CacheTranscript {
    let strand = tr.strand;

    let mapper_result = mapper::build_mapper(&tr.exons, &tr.cds_regions, strand);

    let exons = build_exons(&tr.exons, &tr.cds_regions, strand);

    let sorted_exons = {
        let mut se = exons.clone();
        se.sort_by_key(|e| e.rank);
        se
    };

    let introns = build_introns(&exons, strand);

    let peptide = protein_sequences.get(&tr.transcript_id).cloned();

    let protein_id = if !tr.cds_regions.is_empty() {
        derive_protein_id(&tr.transcript_id)
    } else {
        None
    };

    let translation = if !tr.cds_regions.is_empty() {
        let cdna_start = mapper_result.cdna_coding_start.unwrap_or(1);
        let cdna_end = mapper_result.cdna_coding_end.unwrap_or(1);
        Some(CacheTranslation {
            stable_id: protein_id.clone().unwrap_or_default(),
            version: None,
            db_id: None,
            start: cdna_start,
            end: cdna_end,
        })
    } else {
        None
    };

    let cache_mapper = if !mapper_result.pairs.is_empty() {
        Some(CacheMapper {
            start_phase: mapper_result.start_phase,
            cdna_coding_start: mapper_result.cdna_coding_start,
            cdna_coding_end: mapper_result.cdna_coding_end,
            pair_count: mapper_result.pairs.len(),
            pairs: mapper_result
                .pairs
                .iter()
                .map(|p| CacheMapperPair {
                    from_start: p.from_start,
                    from_end: p.from_end,
                    from_id: "cdna",
                    to_start: p.to_start,
                    to_end: p.to_end,
                    to_id: "genome",
                    ori: p.ori,
                })
                .collect(),
        })
    } else {
        None
    };

    let gene_symbol_source = gene.name.as_ref().map(|_| {
        // Ensembl GFF3 doesn't explicitly provide symbol source;
        // most human gene symbols come from HGNC.
        "HGNC".to_string()
    });

    // Attribute codes as VEP's cache dumper keeps them: `gencode_basic` from
    // the GFF3 `basic` tag, and the GTF's `cds_start_NF` / `cds_end_NF` /
    // `gencode_primary` when a GTF is given (the GFF3 lacks them).
    let mut attributes = Vec::new();
    let mut codes: Vec<&str> = Vec::new();
    if tr.gencode_basic {
        codes.push("gencode_basic");
    }
    if let Some(gtf) = gtf_tags {
        for code in gtf.kept_attribute_codes(&tr.transcript_id) {
            if !codes.contains(&code) {
                codes.push(code);
            }
        }
    }
    for code in codes {
        attributes.push(CacheAttribute {
            code: code.to_string(),
            value: "1".to_string(),
        });
    }
    if let Some(tsl) = tr.tsl {
        attributes.push(CacheAttribute {
            code: "TSL".to_string(),
            value: format!("tsl{tsl}"),
        });
    }
    if let Some(ref appris) = tr.appris {
        attributes.push(CacheAttribute {
            code: "appris".to_string(),
            value: appris.clone(),
        });
    }

    let translateable_seq = genome_fasta
        .and_then(|fasta| compute_translateable_seq(&tr.chr, &tr.cds_regions, strand, fasta));

    CacheTranscript {
        stable_id: tr.transcript_id.clone(),
        version: tr.version,
        db_id: None,
        gene_stable_id: gene.gene_id.clone(),
        chr: None, // Not emitted; derived from directory structure.
        start: tr.start,
        end: tr.end,
        strand,
        biotype: tr.biotype.clone(),
        source: tr.source.clone(),
        description: gene.description.clone(),
        gene_symbol: gene.name.clone(),
        gene_symbol_source,
        hgnc_id: gene.hgnc_id.clone(),
        gene_phenotype: Some(0),
        canonical: if tr.is_canonical { 1 } else { 0 },
        mane_select: tr.mane_select.clone(),
        mane_plus_clinical: tr.mane_plus_clinical.clone(),
        tsl: tr.tsl,
        appris: tr.appris.clone(),
        protein: protein_id,
        exons,
        cdna_coding_start: mapper_result.cdna_coding_start,
        cdna_coding_end: mapper_result.cdna_coding_end,
        translation,
        attributes,
        variation_effect_feature_cache: CacheVEFC {
            codon_table: 1,
            introns,
            mapper: cache_mapper,
            peptide,
            translateable_seq,
            sorted_exons,
            // Protein features (InterPro domains, Pfam) are absent from Ensembl GFF3
            // and would need a separate source (Ensembl BioMart/REST, InterPro GFF3,
            // or the Perl VEP cache's domain data), so this stays empty and DOMAINS
            // is not annotated.
            protein_features: Vec::new(),
            protein_function_predictions: prediction_data
                .and_then(|pd| pd.get(&tr.transcript_id))
                .map(convert_predictions),
        },
    }
}

/// Compute the translateable (CDS) nucleotide sequence from the reference genome FASTA.
///
/// This concatenates the coding exon sequences in transcript order (5'→3'),
/// reverse-complementing for reverse-strand genes. The result is the same
/// `translateable_seq` field that Perl VEP embeds in its Storable cache.
///
/// Without this field, the VEP runtime cannot determine codons and falls back
/// to the generic `coding_sequence_variant` instead of specific consequences
/// like `missense_variant`, `synonymous_variant`, or `stop_gained`.
fn compute_translateable_seq(
    chr: &str,
    cds_regions: &[crate::gff3::CdsRecord],
    strand: i8,
    fasta: &vep_fasta::IndexedFasta,
) -> Option<String> {
    if cds_regions.is_empty() {
        return None;
    }

    let mut sorted: Vec<&crate::gff3::CdsRecord> = cds_regions.iter().collect();
    sorted.sort_by_key(|c| c.start);

    // For reverse strand, process in reverse genomic order (3'→5' genomically = 5'→3' in mRNA).
    if strand < 0 {
        sorted.reverse();
    }

    let mut seq = String::new();
    for cds in &sorted {
        let bases = fasta.sequence(chr, cds.start, cds.end)?;
        if strand < 0 {
            let rc = vep_core::codon::reverse_complement(&bases);
            seq.push_str(&String::from_utf8_lossy(&rc));
        } else {
            seq.push_str(&String::from_utf8_lossy(&bases));
        }
    }

    if seq.is_empty() {
        None
    } else {
        Some(seq)
    }
}

/// Convert fetched prediction data into the JSON-serializable cache format.
fn convert_predictions(preds: &crate::predictions::PeptidePredictions) -> CachePredictions {
    fn convert_blob(blob: &crate::predictions::PredictionBlob) -> CachePredictionMatrix {
        CachePredictionMatrix {
            analysis: blob.analysis.clone(),
            peptide_length: blob.peptide_length,
            predictions_data: blob.base64_data.clone(),
        }
    }

    CachePredictions {
        sift: preds.sift.as_ref().map(convert_blob),
        polyphen_humdiv: preds.polyphen_humdiv.as_ref().map(convert_blob),
        polyphen_humvar: preds.polyphen_humvar.as_ref().map(convert_blob),
    }
}

/// Build exon list with phases computed from CDS overlap.
fn build_exons(
    gff3_exons: &[crate::gff3::ExonRecord],
    cds_regions: &[crate::gff3::CdsRecord],
    strand: i8,
) -> Vec<CacheExon> {
    let mut sorted: Vec<&crate::gff3::ExonRecord> = gff3_exons.iter().collect();
    sorted.sort_by(|a, b| match (a.rank, b.rank) {
        (Some(ra), Some(rb)) => ra.cmp(&rb),
        _ => a.start.cmp(&b.start),
    });

    sorted
        .iter()
        .enumerate()
        .map(|(i, exon)| {
            let rank = exon.rank.unwrap_or((i + 1) as u32);
            let (phase, end_phase) = compute_exon_phases(exon, cds_regions, strand);

            CacheExon {
                stable_id: exon.exon_id.clone(),
                start: exon.start,
                end: exon.end,
                rank,
                phase,
                end_phase,
                strand,
            }
        })
        .collect()
}

/// Compute phase and end_phase for an exon based on CDS overlap.
fn compute_exon_phases(
    exon: &crate::gff3::ExonRecord,
    cds_regions: &[crate::gff3::CdsRecord],
    _strand: i8,
) -> (i8, i8) {
    let overlapping: Vec<&crate::gff3::CdsRecord> = cds_regions
        .iter()
        .filter(|cds| cds.start <= exon.end && cds.end >= exon.start)
        .collect();

    if overlapping.is_empty() {
        return (-1, -1);
    }

    // Phase at the start of the exon comes from the CDS phase attribute.
    let phase = overlapping[0].phase;

    // End phase: (phase + CDS_length_in_this_exon) % 3.
    let cds_start = overlapping
        .iter()
        .map(|c| c.start.max(exon.start))
        .min()
        .unwrap();
    let cds_end = overlapping
        .iter()
        .map(|c| c.end.min(exon.end))
        .max()
        .unwrap();
    let cds_len_in_exon = cds_end - cds_start + 1;

    let effective_phase = if phase < 0 { 0u64 } else { phase as u64 };
    let end_phase = ((effective_phase + cds_len_in_exon) % 3) as i8;

    (phase, end_phase)
}

/// Build introns from exon list.
fn build_introns(exons: &[CacheExon], strand: i8) -> Vec<CacheIntron> {
    if exons.len() < 2 {
        return Vec::new();
    }

    let mut sorted: Vec<&CacheExon> = exons.iter().collect();
    sorted.sort_by_key(|e| e.start);

    sorted
        .windows(2)
        .enumerate()
        .map(|(i, pair)| CacheIntron {
            start: pair[0].end + 1,
            end: pair[1].start - 1,
            rank: (i + 1) as u32,
            strand,
        })
        .collect()
}

/// Derive a protein (ENSP) ID from a transcript (ENST) ID.
///
/// Convention: replace the ENST prefix with ENSP. Ensembl assigns ENSP accessions
/// independently of ENST accessions, so this derived id is a placeholder.
fn derive_protein_id(transcript_id: &str) -> Option<String> {
    if transcript_id.starts_with("ENST") {
        Some(transcript_id.replacen("ENST", "ENSP", 1))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_derive_protein_id() {
        assert_eq!(
            derive_protein_id("ENST00000309812"),
            Some("ENSP00000309812".to_string())
        );
        assert_eq!(derive_protein_id("NM_001234"), None);
    }

    #[test]
    fn test_region_binning() {
        // Position 1 -> region 1-1000000
        assert_eq!(region_start(1), 1);
        // Position 1000000 -> region 1-1000000
        assert_eq!(region_start(1_000_000), 1);
        // Position 1000001 -> region 1000001-2000000
        assert_eq!(region_start(1_000_001), 1_000_001);
    }

    #[test]
    fn test_generate_info_json() {
        let tmpdir = std::env::temp_dir().join("vep_cache_builder_test_info_json");
        let _ = std::fs::remove_dir_all(&tmpdir);
        std::fs::create_dir_all(&tmpdir).unwrap();

        let versions = HashMap::from([("assembly".to_string(), "GRCh37.p13".to_string())]);
        generate_info_json(&tmpdir, "homo_sapiens", "GRCh37", 115, false, versions).unwrap();

        let info_path = tmpdir.join("info.json");
        assert!(info_path.exists());

        let content = std::fs::read_to_string(&info_path).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();
        assert_eq!(parsed["species"], "homo_sapiens");
        assert_eq!(parsed["assembly"], "GRCh37");
        assert_eq!(parsed["cache_version"], 115);
        assert_eq!(parsed["regulatory"], false);
        assert!(parsed["variation_cols"].is_array());
        assert_eq!(parsed["variation_cols"][0], "variation_name");
        assert_eq!(parsed["source_versions"]["assembly"], "GRCh37.p13");

        let _ = std::fs::remove_dir_all(&tmpdir);
    }
}
