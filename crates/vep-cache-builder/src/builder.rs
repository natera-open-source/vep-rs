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

use anyhow::{bail, Context, Result};
use serde::Serialize;
use tracing::info;
use vep_core::coordinate::is_cache_path_component;

use crate::core_dumps::{CoreDumps, TranslationFacts};
use crate::gff3::{CdsRecord, ExonRecord, GeneRecord};
use crate::mapper;

/// Region size for binning transcripts (1 Mb).
const REGION_SIZE: u64 = 1_000_000;

/// First 1-based position of the 1 Mb region containing `pos`.
fn region_start(pos: u64) -> u64 {
    ((pos - 1) / REGION_SIZE) * REGION_SIZE + 1
}

/// Statistics from the cache build.
pub struct BuildStats {
    pub transcript_count: usize,
    pub chromosome_count: usize,
    pub region_file_count: usize,
    /// The chromosome directories written, in the order `info.json` lists them.
    pub chromosomes: Vec<String>,
    /// Transcripts written with a SIFT matrix.
    pub sift_matrices: usize,
    /// Transcripts written with a PolyPhen matrix, HumVar or HumDiv.
    pub polyphen_matrices: usize,
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
    gene_version: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    chr: Option<String>,
    start: u64,
    end: u64,
    strand: i8,
    biotype: String,
    source: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    gene_symbol: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    gene_symbol_source: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none", rename = "gene_hgnc_id")]
    hgnc_id: Option<String>,
    /// The gene's significant phenotype feature count.
    #[serde(skip_serializing_if = "Option::is_none")]
    gene_phenotype: Option<u32>,
    /// `1` on the gene's canonical transcript; the Ensembl cache carries no
    /// key on the others.
    #[serde(rename = "is_canonical", skip_serializing_if = "Option::is_none")]
    canonical: Option<u8>,
    /// The first CCDS xref's display label, or the dumper's `-` sentinel.
    ccds: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    protein: Option<String>,
    /// Comma-joined xref labels, or the dumper's `-` sentinel; absent when
    /// no core dumps were read.
    #[serde(skip_serializing_if = "Option::is_none")]
    refseq: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    swissprot: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    trembl: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    uniparc: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    uniprot_isoform: Option<String>,
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
    /// The coding slice of the spliced cDNA, led by one `N` per base of the
    /// translation start exon's phase; the sequence behind codon-level consequences.
    #[serde(skip_serializing_if = "Option::is_none")]
    translateable_seq: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    five_prime_utr: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    three_prime_utr: Option<String>,
    /// The exons in genomic order.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    sorted_exons: Vec<CacheExon>,
    /// The translation's sequence edits; always written, empty for most translations.
    seq_edits: Vec<CacheSeqEdit>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    protein_features: Vec<CacheProteinFeature>,
    #[serde(skip_serializing_if = "Option::is_none")]
    protein_function_predictions: Option<CachePredictions>,
}

/// A translation sequence edit (`Bio::EnsEMBL::SeqEdit`): `start` and `end` are
/// protein positions, `alt_seq` the residues written over them and `length_diff`
/// the growth they cause.
#[derive(Serialize, Clone, Debug, PartialEq, Eq)]
pub struct CacheSeqEdit {
    pub code: String,
    pub name: String,
    pub start: u64,
    pub end: u64,
    pub alt_seq: String,
    pub length_diff: i64,
}

impl CacheSeqEdit {
    /// A seq edit from a `translation_attrib` row of one of the edit codes
    /// `Translation::get_all_SeqEdits` reads (core `Translation.pm:1111`), whose value
    /// is `start end alt_seq` split on whitespace (core `SeqEdit.pm:144`); the name
    /// is the code's `attrib_type.name`. `None` when the value lacks two integer
    /// coordinates.
    pub fn from_attrib(code: &str, value: &str) -> Option<Self> {
        let name = match code {
            "initial_met" => "Initial methionine",
            "_selenocysteine" => "Selenocysteine",
            "amino_acid_sub" => "Amino acid substitution",
            "_stop_codon_rt" => "Stop Codon Readthrough",
            _ => "",
        };
        let mut parts = value.split_whitespace();
        let start: u64 = parts.next()?.parse().ok()?;
        let end: u64 = parts.next()?.parse().ok()?;
        let alt_seq = parts.next().unwrap_or("").to_string();
        Some(Self {
            code: code.to_string(),
            name: name.to_string(),
            start,
            end,
            length_diff: alt_seq.len() as i64 - (end as i64 - start as i64 + 1),
            alt_seq,
        })
    }
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

/// A single prediction matrix in JSON-serializable form: the keys of the
/// Ensembl matrix object, with the gzip blob base64-encoded as
/// `predictions_data`.
#[derive(Serialize)]
struct CachePredictionMatrix {
    analysis: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    sub_analysis: Option<String>,
    translation_md5: String,
    predictions_data: String,
}

#[derive(Serialize)]
struct CacheIntron {
    start: u64,
    end: u64,
    rank: u32,
}

#[derive(Serialize)]
struct CacheMapper {
    start_phase: i8,
    #[serde(skip_serializing_if = "Option::is_none")]
    cdna_coding_start: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    cdna_coding_end: Option<u64>,
    pairs: Vec<CacheMapperPair>,
}

#[derive(Serialize)]
struct CacheMapperPair {
    from_start: u64,
    from_end: u64,
    to_start: u64,
    to_end: u64,
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
    core_dumps: Option<&CoreDumps>,
) -> Result<BuildStats> {
    let mut by_chr: HashMap<String, Vec<CacheTranscript>> = HashMap::new();
    let mut transcript_count = 0usize;
    let mut readthrough_dropped = 0usize;

    for gene in genes.values() {
        for tr in &gene.transcripts {
            // The dumper drops a readthrough transcript by its attribute; the
            // GRCh37 release-87 GFF3 has no tag for it, the core dumps do.
            if core_dumps
                .and_then(|d| d.transcript(&tr.transcript_id))
                .is_some_and(|f| f.readthrough)
            {
                readthrough_dropped += 1;
                continue;
            }
            let cache_tr = build_transcript(
                gene,
                tr,
                protein_sequences,
                prediction_data,
                genome_fasta,
                gtf_tags,
                core_dumps,
            );
            by_chr.entry(tr.chr.clone()).or_default().push(cache_tr);
            transcript_count += 1;
        }
    }
    if readthrough_dropped > 0 {
        info!(
            "Dropped {} readthrough transcripts named by the core dumps",
            readthrough_dropped
        );
    }

    // Position order: the dumper's region fetch returns genes that way, the
    // same-named HGNC copy below depends on it, and the shards then come out
    // the same on every build.
    for transcripts in by_chr.values_mut() {
        transcripts
            .sort_by(|a, b| (a.start, a.end, &a.stable_id).cmp(&(b.start, b.end, &b.stable_id)));
        copy_hgnc_ids_by_symbol(transcripts);
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
        // Refused at the GFF3 parse as well; the function that makes the
        // directory does not rely on its caller to have done so.
        if !is_cache_path_component(chr) || chr_dir.parent() != Some(transcripts_dir.as_path()) {
            bail!("invalid chromosome name '{chr}'");
        }
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

    let mut sift_matrices = 0usize;
    let mut polyphen_matrices = 0usize;
    for tr in by_chr.values().flatten() {
        let Some(predictions) = &tr
            .variation_effect_feature_cache
            .protein_function_predictions
        else {
            continue;
        };
        if predictions.sift.is_some() {
            sift_matrices += 1;
        }
        if predictions.polyphen_humvar.is_some() || predictions.polyphen_humdiv.is_some() {
            polyphen_matrices += 1;
        }
    }
    Ok(BuildStats {
        transcript_count,
        chromosome_count: by_chr.len(),
        region_file_count,
        chromosomes: crate::info::contig_order(by_chr.keys().cloned()),
        sift_matrices,
        polyphen_matrices,
    })
}

/// Ensembl VEP's dumper writes each 1 Mb region through `merge_features`
/// (AnnotationType/Transcript.pm), which gives every transcript whose gene
/// symbol matches that of a transcript carrying an HGNC id, among the
/// transcripts overlapping the region, that HGNC id; a gene whose display
/// xref is EntrezGene so prints the HGNC id of a same-named HGNC gene nearby.
/// The map is built from the transcripts overlapping each region, in the
/// order given (a later transcript's id wins, as in `merge_features`), and
/// applied to the transcripts the region's shard holds.
fn copy_hgnc_ids_by_symbol(transcripts: &mut [CacheTranscript]) {
    let mut by_region: HashMap<u64, HashMap<String, String>> = HashMap::new();
    for tr in transcripts.iter() {
        let (Some(symbol), Some(hgnc_id)) = (&tr.gene_symbol, &tr.hgnc_id) else {
            continue;
        };
        let mut region = region_start(tr.start);
        while region <= tr.end {
            by_region
                .entry(region)
                .or_default()
                .insert(symbol.clone(), hgnc_id.clone());
            region += REGION_SIZE;
        }
    }
    for tr in transcripts.iter_mut() {
        let Some(symbol) = &tr.gene_symbol else {
            continue;
        };
        if let Some(hgnc_id) = by_region
            .get(&region_start(tr.start))
            .and_then(|m| m.get(symbol))
        {
            tr.hgnc_id = Some(hgnc_id.clone());
        }
    }
}

/// Build a single CacheTranscript from parsed GFF3 data.
fn build_transcript(
    gene: &GeneRecord,
    tr: &crate::gff3::TranscriptRecord,
    protein_sequences: &HashMap<String, String>,
    prediction_data: Option<&HashMap<String, crate::predictions::PeptidePredictions>>,
    genome_fasta: Option<&vep_fasta::IndexedFasta>,
    gtf_tags: Option<&crate::gtf::GtfTags>,
    core_dumps: Option<&CoreDumps>,
) -> CacheTranscript {
    let strand = tr.strand;
    let facts = core_dumps.and_then(|d| d.transcript(&tr.transcript_id));
    let gene_facts = core_dumps.and_then(|d| d.gene(&gene.gene_id));
    let translation_facts = facts.and_then(|f| f.translation.as_ref());

    let ordered = transcript_order(&tr.exons, strand);
    let phases = exon_phases(&ordered, &tr.cds_regions, strand);
    // The mapper's start phase is exon 1's phase (core TranscriptMapper.pm:110).
    let first_exon_phase = phases.first().map_or(-1, |p| p.0);
    let mapper_result = mapper::build_mapper(&tr.exons, &tr.cds_regions, strand, first_exon_phase);

    let exons = build_exons(&ordered, &phases);

    // Database/Transcript.pm:350 sorts the exons by genomic start.
    let sorted_exons = {
        let mut se = exons.clone();
        se.sort_by_key(|e| (e.start, e.end));
        se
    };

    let introns = build_introns(&exons);

    let peptide = protein_sequences.get(&tr.transcript_id).cloned();

    let coding = match (
        genome_fasta,
        mapper_result.cdna_coding_start,
        mapper_result.cdna_coding_end,
    ) {
        (Some(fasta), Some(cdna_start), Some(cdna_end)) => {
            // The padding follows the translation start exon's phase (core
            // Transcript.pm:917-919), which is exon 1's only when exon 1 codes.
            let span = cds_span(&tr.cds_regions);
            let start_exon_phase = ordered
                .iter()
                .zip(&phases)
                .find(|(e, _)| span.is_some_and(|(s, end)| e.start <= end && e.end >= s))
                .map_or(-1, |(_, p)| p.0);
            spliced_seq(&tr.chr, &ordered, strand, fasta)
                .map(|spliced| coding_sequences(&spliced, cdna_start, cdna_end, start_exon_phase))
        }
        _ => None,
    };
    // A transcript without a translation has the empty translateable sequence
    // (core Transcript.pm:906-908), whether or not a genome is given.
    let (translateable_seq, five_prime_utr, three_prime_utr) = match coding {
        Some((t, five, three)) => (Some(t), five, three),
        None if tr.cds_regions.is_empty() => (Some(String::new()), None, None),
        None => (None, None, None),
    };

    // The translation's stable identifier and version are the CDS line's
    // `protein_id` and `version`; every CDS segment of a transcript carries
    // the same pair.
    let protein_id = tr
        .cds_regions
        .iter()
        .find_map(|c| c.protein_id.clone())
        .or_else(|| translation_facts.map(|t| t.stable_id.clone()));
    let protein_version = tr
        .cds_regions
        .iter()
        .find_map(|c| c.protein_version)
        .or_else(|| translation_facts.and_then(|t| t.version));

    let translation = translation_bounds(&ordered, &tr.cds_regions, strand).map(|(start, end)| {
        CacheTranslation {
            stable_id: protein_id.clone().unwrap_or_default(),
            version: protein_version,
            db_id: None,
            start,
            end,
        }
    });

    let cache_mapper = if !mapper_result.pairs.is_empty() {
        Some(CacheMapper {
            start_phase: mapper_result.start_phase,
            cdna_coding_start: mapper_result.cdna_coding_start,
            cdna_coding_end: mapper_result.cdna_coding_end,
            pairs: mapper_result
                .pairs
                .iter()
                .map(|p| CacheMapperPair {
                    from_start: p.from_start,
                    from_end: p.from_end,
                    to_start: p.to_start,
                    to_end: p.to_end,
                    ori: p.ori,
                })
                .collect(),
        })
    } else {
        None
    };

    // Attribute codes and values as Ensembl VEP's cache dumper keeps them
    // (Pipeline/DumpVEP/Dumper/Core.pm `clean_transcript`): the core dumps'
    // kept transcript attributes verbatim when read, else what the GFF3 and
    // the GTF state.
    let attributes: Vec<CacheAttribute> = match facts {
        Some(f) => f
            .attributes
            .iter()
            .map(|a| CacheAttribute {
                code: a.code.clone(),
                value: a.value.clone(),
            })
            .collect(),
        None => gff3_attributes(tr, gtf_tags),
    };

    // Xref lists as `prefetch_transcript_ids` and `prefetch_translation_ids`
    // store them: comma-joined labels, `-` when there are none. Without the
    // core dumps the transcript's are unknown; `ccds` is on the GFF3.
    let xref_list = |labels: &[String]| -> String {
        if labels.is_empty() {
            "-".to_string()
        } else {
            labels.join(",")
        }
    };
    let ccds = facts
        .and_then(|f| f.ccds.clone())
        .or_else(|| tr.ccds.clone())
        .unwrap_or_else(|| "-".to_string());
    let refseq = facts.map(|f| xref_list(&f.refseq));
    let translation_xrefs = |pick: fn(&TranslationFacts) -> &[String]| -> Option<String> {
        if tr.cds_regions.is_empty() || facts.is_none() {
            return None;
        }
        match translation_facts {
            Some(t) => Some(xref_list(pick(t))),
            None => Some("-".to_string()),
        }
    };

    let protein_features = translation_facts
        .map(|t| {
            t.protein_features
                .iter()
                .map(|pf| CacheProteinFeature {
                    start: pf.start,
                    end: pf.end,
                    hseqname: pf.hseqname.clone(),
                    analysis: pf.display_label.clone(),
                })
                .collect()
        })
        .unwrap_or_default();

    let (gene_symbol, gene_symbol_source, hgnc_id, gene_phenotype) = match gene_facts {
        Some(g) => (
            g.symbol.clone(),
            g.symbol_source.clone(),
            g.hgnc_id.clone(),
            g.phenotype_count,
        ),
        None => (gene.name.clone(), None, None, None),
    };

    let seq_edits: Vec<CacheSeqEdit> = translation_facts
        .map(|t| {
            t.seq_edits
                .iter()
                .filter_map(|e| CacheSeqEdit::from_attrib(&e.code, &e.value))
                .collect()
        })
        .unwrap_or_default();

    CacheTranscript {
        stable_id: tr.transcript_id.clone(),
        version: tr.version,
        db_id: None,
        gene_stable_id: gene.gene_id.clone(),
        gene_version: gene.version,
        chr: None, // Not emitted; derived from directory structure.
        start: tr.start,
        end: tr.end,
        strand,
        biotype: tr.biotype.clone(),
        source: tr.source.clone(),
        gene_symbol,
        gene_symbol_source,
        hgnc_id,
        gene_phenotype,
        canonical: (tr.is_canonical || facts.is_some_and(|f| f.is_canonical)).then_some(1),
        ccds,
        protein: protein_id,
        refseq,
        swissprot: translation_xrefs(|t| &t.swissprot),
        trembl: translation_xrefs(|t| &t.trembl),
        uniparc: translation_xrefs(|t| &t.uniparc),
        uniprot_isoform: translation_xrefs(|t| &t.uniprot_isoform),
        exons,
        cdna_coding_start: mapper_result.cdna_coding_start,
        cdna_coding_end: mapper_result.cdna_coding_end,
        translation,
        attributes,
        variation_effect_feature_cache: CacheVEFC {
            codon_table: codon_table(&tr.chr, core_dumps),
            introns,
            mapper: cache_mapper,
            peptide,
            translateable_seq,
            five_prime_utr,
            three_prime_utr,
            sorted_exons,
            seq_edits,
            protein_features,
            protein_function_predictions: prediction_data
                .and_then(|pd| pd.get(&tr.transcript_id))
                .map(convert_predictions),
        },
    }
}

/// The transcript attributes a build without the core dumps can state, each
/// code once: the GFF3's GENCODE set labels, the GTF's `cds_start_NF`,
/// `cds_end_NF` and `gencode_primary` when a GTF is given, the `TSL` value
/// with its `tsl` prefix, and a MANE membership whose accession the GFF3 does
/// not carry (the attribute present, its value empty).
fn gff3_attributes(
    tr: &crate::gff3::TranscriptRecord,
    gtf_tags: Option<&crate::gtf::GtfTags>,
) -> Vec<CacheAttribute> {
    let mut attributes: Vec<CacheAttribute> = Vec::new();
    let mut push_attribute = |code: &str, value: &str| {
        if !attributes.iter().any(|a| a.code == code) {
            attributes.push(CacheAttribute {
                code: code.to_string(),
                value: value.to_string(),
            });
        }
    };
    if tr.gencode_basic {
        push_attribute("gencode_basic", "GENCODE basic");
    }
    if tr.gencode_primary {
        push_attribute("gencode_primary", "GENCODE primary");
    }
    if let Some(gtf) = gtf_tags {
        for code in gtf.kept_attribute_codes(&tr.transcript_id) {
            let value = match code {
                "gencode_basic" => "GENCODE basic",
                "gencode_primary" => "GENCODE primary",
                _ => "1",
            };
            push_attribute(code, value);
        }
    }
    if let Some(tsl) = &tr.tsl {
        push_attribute("TSL", &format!("tsl{tsl}"));
    }
    if tr.mane_select {
        push_attribute("MANE_Select", "");
    }
    if tr.mane_plus_clinical {
        push_attribute("MANE_Plus_Clinical", "");
    }
    attributes
}

/// A transcript's exons in transcript order, 5' to 3': by `rank` when every exon
/// carries one, else by genomic start, descending on the reverse strand.
fn transcript_order(exons: &[ExonRecord], strand: i8) -> Vec<&ExonRecord> {
    let mut ordered: Vec<&ExonRecord> = exons.iter().collect();
    if ordered.iter().all(|e| e.rank.is_some()) {
        ordered.sort_by_key(|e| e.rank);
    } else if strand < 0 {
        ordered.sort_by_key(|e| std::cmp::Reverse(e.start));
    } else {
        ordered.sort_by_key(|e| e.start);
    }
    ordered
}

/// The genomic span of the coding sequence, `None` for a non-coding transcript.
fn cds_span(cds_regions: &[CdsRecord]) -> Option<(u64, u64)> {
    let start = cds_regions.iter().map(|c| c.start).min()?;
    let end = cds_regions.iter().map(|c| c.end).max()?;
    Some((start, end))
}

/// The Ensembl phase and end phase of every exon of `ordered`, in transcript order.
///
/// Each is the GFF3 exon line's `ensembl_phase` and `ensembl_end_phase` pair, the
/// exon's own `phase` and `end_phase` (core `Exon.pm:1651-1659`). An exon line
/// without the pair takes the derived one: -1 at an exon end that lies in
/// untranslated sequence; at a coding 5' end the frame of the exon's CDS segment,
/// the GFF3 phase column with 1 and 2 swapped (the column counts the bases before
/// the next full codon, the Ensembl phase the bases of the first codon that sit in
/// the previous exon; core `Transcript.pm:1977` applies the same swap the other
/// way); at a coding 3' end `(phase + coding length) % 3` with -1 counted as 0.
fn exon_phases(ordered: &[&ExonRecord], cds_regions: &[CdsRecord], strand: i8) -> Vec<(i8, i8)> {
    let span = cds_span(cds_regions);
    ordered
        .iter()
        .map(|exon| {
            if let (Some(phase), Some(end_phase)) = (exon.ensembl_phase, exon.ensembl_end_phase) {
                return (phase, end_phase);
            }
            let Some((cds_start, cds_end)) = span else {
                return (-1, -1);
            };
            let coding_start = exon.start.max(cds_start);
            let coding_end = exon.end.min(cds_end);
            if coding_start > coding_end {
                return (-1, -1);
            }
            let (five_prime_coding, three_prime_coding) = if strand < 0 {
                (cds_end >= exon.end, cds_start <= exon.start)
            } else {
                (cds_start <= exon.start, cds_end >= exon.end)
            };
            let phase = if five_prime_coding {
                let gff3_phase = cds_regions
                    .iter()
                    .find(|c| c.start <= exon.end && c.end >= exon.start)
                    .map_or(0, |c| c.phase);
                match gff3_phase {
                    1 => 2,
                    2 => 1,
                    p => p.max(0),
                }
            } else {
                -1
            };
            let end_phase = if three_prime_coding {
                ((phase.max(0) as u64 + (coding_end - coding_start + 1)) % 3) as i8
            } else {
                -1
            };
            (phase, end_phase)
        })
        .collect()
}

/// The spliced cDNA of the exons in transcript order, read from the genome FASTA
/// and reverse-complemented on the reverse strand (core `Transcript.pm` `spliced_seq`
/// before its sequence edits); `None` when an exon's sequence is missing.
fn spliced_seq(
    chr: &str,
    ordered: &[&ExonRecord],
    strand: i8,
    fasta: &vep_fasta::IndexedFasta,
) -> Option<String> {
    let mut seq = String::new();
    for exon in ordered {
        let bases = fasta.sequence(chr, exon.start, exon.end)?;
        if strand < 0 {
            let rc = vep_core::codon::reverse_complement(&bases);
            seq.push_str(&String::from_utf8_lossy(&rc));
        } else {
            seq.push_str(&String::from_utf8_lossy(&bases));
        }
    }
    (!seq.is_empty()).then_some(seq)
}

/// The translateable sequence and the two UTRs of a coding transcript, cut from
/// its spliced cDNA.
///
/// `translateable_seq` is the slice from `cdna_coding_start` to `cdna_coding_end`
/// led by one `N` per base of the translation start exon's phase (core
/// `Transcript.pm:903-926`), `five_prime_utr` the bases before it (`:1746-1762`) and
/// `three_prime_utr` the bases after it (`:1781-1797`), each `None` when empty.
fn coding_sequences(
    spliced: &str,
    cdna_coding_start: u64,
    cdna_coding_end: u64,
    start_exon_phase: i8,
) -> (String, Option<String>, Option<String>) {
    let coding_start = (cdna_coding_start as usize)
        .saturating_sub(1)
        .min(spliced.len());
    let coding_end = (cdna_coding_end as usize)
        .min(spliced.len())
        .max(coding_start);
    let padding = "N".repeat(start_exon_phase.max(0) as usize);
    let translateable = format!("{padding}{}", &spliced[coding_start..coding_end]);
    let five = &spliced[..coding_start];
    let three = &spliced[coding_end..];
    (
        translateable,
        (!five.is_empty()).then(|| five.to_string()),
        (!three.is_empty()).then(|| three.to_string()),
    )
}

/// The translation's `start` and `end`: the first coding base's position within the
/// translation start exon and the last coding base's within the end exon, each
/// 1-based in transcript orientation (core `Translation.pm:190-226`); `None` for a
/// non-coding transcript.
fn translation_bounds(
    ordered: &[&ExonRecord],
    cds_regions: &[CdsRecord],
    strand: i8,
) -> Option<(u64, u64)> {
    let (cds_start, cds_end) = cds_span(cds_regions)?;
    let overlaps = |e: &&&ExonRecord| e.start <= cds_end && e.end >= cds_start;
    let start_exon = ordered.iter().find(overlaps)?;
    let end_exon = ordered.iter().rev().find(overlaps)?;
    Some(if strand < 0 {
        (start_exon.end - cds_end + 1, end_exon.end - cds_start + 1)
    } else {
        (
            cds_start - start_exon.start + 1,
            cds_end - end_exon.start + 1,
        )
    })
}

/// The codon table of the transcript's sequence region: its `codon_table` attribute
/// when the core dumps carry one for the region, else 2 on `MT` and 1 elsewhere.
///
/// Ensembl VEP reads the slice attribute and defaults to 1
/// (`Database/Transcript.pm:359-363`); on both human assemblies the mitochondrion
/// is the only region carrying the attribute, with value 2.
fn codon_table(chr: &str, core_dumps: Option<&CoreDumps>) -> u8 {
    core_dumps
        .and_then(|d| d.codon_table(chr))
        .unwrap_or(if chr == "MT" { 2 } else { 1 })
}

/// Convert fetched prediction data into the JSON-serializable cache format.
fn convert_predictions(preds: &crate::predictions::PeptidePredictions) -> CachePredictions {
    fn convert_blob(blob: &crate::predictions::PredictionBlob) -> CachePredictionMatrix {
        CachePredictionMatrix {
            analysis: blob.analysis.clone(),
            sub_analysis: blob.sub_analysis.clone(),
            translation_md5: blob.translation_md5.clone(),
            predictions_data: blob.base64_data.clone(),
        }
    }

    CachePredictions {
        sift: preds.sift.as_ref().map(convert_blob),
        polyphen_humdiv: preds.polyphen_humdiv.as_ref().map(convert_blob),
        polyphen_humvar: preds.polyphen_humvar.as_ref().map(convert_blob),
    }
}

/// The cached exon list in transcript order, each exon with its Ensembl phases
/// and its `rank`: the line's when every exon of the transcript carries one,
/// else its position in the list, the order [`transcript_order`] fell back to.
fn build_exons(ordered: &[&ExonRecord], phases: &[(i8, i8)]) -> Vec<CacheExon> {
    let ranks_from_file = ordered.iter().all(|e| e.rank.is_some());
    ordered
        .iter()
        .zip(phases)
        .enumerate()
        .map(|(i, (exon, &(phase, end_phase)))| CacheExon {
            stable_id: exon.exon_id.clone(),
            start: exon.start,
            end: exon.end,
            rank: match exon.rank {
                Some(rank) if ranks_from_file => rank,
                _ => (i + 1) as u32,
            },
            phase,
            end_phase,
        })
        .collect()
}

/// The introns between consecutive exons of the transcript-ordered list, in that
/// order (core `Transcript.pm:1613-1619`), so the first intron follows exon 1
/// whatever the strand.
fn build_introns(ordered: &[CacheExon]) -> Vec<CacheIntron> {
    ordered
        .windows(2)
        .enumerate()
        .map(|(i, pair)| {
            let (upstream, downstream) = if pair[0].start <= pair[1].start {
                (&pair[0], &pair[1])
            } else {
                (&pair[1], &pair[0])
            };
            CacheIntron {
                start: upstream.end + 1,
                end: downstream.start - 1,
                rank: (i + 1) as u32,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exon(start: u64, end: u64, strand: i8, rank: u32) -> ExonRecord {
        ExonRecord {
            exon_id: Some(format!("E{rank}")),
            start,
            end,
            strand,
            rank: Some(rank),
            ensembl_phase: None,
            ensembl_end_phase: None,
        }
    }

    fn cds(start: u64, end: u64, strand: i8, phase: i8) -> CdsRecord {
        CdsRecord {
            start,
            end,
            strand,
            phase,
            protein_id: None,
            protein_version: None,
        }
    }

    /// A one-contig genome on `21`: positions 1 to 60, one 60-base line.
    fn genome(dir: &Path) -> vep_fasta::IndexedFasta {
        let fa = dir.join("ref.fa");
        let seq = "GGGATGCCCAAATTTGGGTAGCCCTTTAAAGGGCCCATGAAACCCGGGTTTTAGCCCAAA";
        assert_eq!(seq.len(), 60);
        std::fs::write(&fa, format!(">21\n{seq}\n")).unwrap();
        std::fs::write(dir.join("ref.fa.fai"), "21\t60\t4\t60\t61\n").unwrap();
        vep_fasta::IndexedFasta::from_path(&fa).unwrap()
    }

    /// An exon's phases are the GFF3 line's `ensembl_phase` and `ensembl_end_phase`
    /// (core `Exon.pm:1651-1659`) whatever the CDS says: the pair 0/0 on an exon that
    /// opens in the 5' UTR is kept where the derivation would read -1, a disagreement
    /// the GRCh38 GFF3 carries. Without the pair, the phases follow Ensembl's
    /// convention: an exon starting in the 5' UTR has phase -1 and end phase
    /// `coding length % 3`; a fully coding exon's phase is its CDS segment's GFF3 phase
    /// with 1 and 2 swapped; an exon ending in the 3' UTR has end phase -1; a
    /// non-coding exon is -1/-1. The same transcript laid out on the reverse strand
    /// reads the same phases exon by exon.
    #[test]
    fn exon_phases_from_ensembl_phase_attributes() {
        let mut attributed = exon(1, 10, 1, 1);
        attributed.ensembl_phase = Some(0);
        attributed.ensembl_end_phase = Some(0);
        let ordered = [&attributed];
        assert_eq!(exon_phases(&ordered, &[cds(6, 10, 1, 0)], 1), vec![(0, 0)]);

        // Forward: exon 1 (1-10) with the CDS from 6; exon 2 (21-30) coding
        // throughout with GFF3 phase 1 (two leading bases of a codon); exon 3
        // (41-50) coding to 44 then UTR; exon 4 (55-60) untranslated.
        let exons = [
            exon(1, 10, 1, 1),
            exon(21, 30, 1, 2),
            exon(41, 50, 1, 3),
            exon(55, 60, 1, 4),
        ];
        let cds_regions = [cds(6, 10, 1, 0), cds(21, 30, 1, 1), cds(41, 44, 1, 0)];
        let ordered = transcript_order(&exons, 1);
        assert_eq!(
            exon_phases(&ordered, &cds_regions, 1),
            vec![(-1, 2), (2, 0), (0, -1), (-1, -1)]
        );

        // Reverse: the same layout with exon 1 at the high coordinates.
        let exons = [
            exon(51, 60, -1, 1),
            exon(31, 40, -1, 2),
            exon(11, 20, -1, 3),
            exon(1, 6, -1, 4),
        ];
        let cds_regions = [cds(51, 55, -1, 0), cds(31, 40, -1, 1), cds(17, 20, -1, 0)];
        let ordered = transcript_order(&exons, -1);
        assert_eq!(
            exon_phases(&ordered, &cds_regions, -1),
            vec![(-1, 2), (2, 0), (0, -1), (-1, -1)]
        );

        // A non-coding transcript is -1 throughout.
        let ordered = transcript_order(&exons, -1);
        assert_eq!(exon_phases(&ordered, &[], -1), vec![(-1, -1); 4]);
    }

    /// `sorted_exons` is the exon list by genomic start (`Database/Transcript.pm:350`),
    /// so a reverse-strand transcript's exon 3 comes first; the `exons` list keeps
    /// transcript order and the introns follow exon 1 first (`Transcript.pm:1613-1619`).
    #[test]
    fn sorted_exons_in_genomic_order() {
        let exons = [exon(41, 50, -1, 1), exon(21, 30, -1, 2), exon(1, 10, -1, 3)];
        let ordered = transcript_order(&exons, -1);
        let phases = exon_phases(&ordered, &[], -1);
        let cache_exons = build_exons(&ordered, &phases);
        assert_eq!(
            cache_exons.iter().map(|e| e.rank).collect::<Vec<_>>(),
            vec![1, 2, 3]
        );
        assert_eq!(
            cache_exons.iter().map(|e| e.start).collect::<Vec<_>>(),
            vec![41, 21, 1]
        );

        let mut sorted = cache_exons.clone();
        sorted.sort_by_key(|e| (e.start, e.end));
        assert_eq!(
            sorted.iter().map(|e| (e.start, e.rank)).collect::<Vec<_>>(),
            vec![(1, 3), (21, 2), (41, 1)]
        );

        let introns = build_introns(&cache_exons);
        assert_eq!(
            introns
                .iter()
                .map(|i| (i.rank, i.start, i.end))
                .collect::<Vec<_>>(),
            vec![(1, 31, 40), (2, 11, 20)]
        );
    }

    /// Once any exon of a transcript lacks a `rank`, the list is in start order,
    /// descending on the reverse strand, and every exon is numbered by its place
    /// in it: the file's ranks on the other exons would otherwise sit among
    /// numbers they were never ordered against.
    #[test]
    fn mixed_rank_exons_are_renumbered_in_transcript_order() {
        let mut exons = [exon(41, 50, -1, 7), exon(21, 30, -1, 2), exon(1, 10, -1, 9)];
        exons[1].rank = None;
        let ordered = transcript_order(&exons, -1);
        let phases = exon_phases(&ordered, &[], -1);
        let cache_exons = build_exons(&ordered, &phases);
        assert_eq!(
            cache_exons
                .iter()
                .map(|e| (e.start, e.rank))
                .collect::<Vec<_>>(),
            vec![(41, 1), (21, 2), (1, 3)]
        );
    }

    /// `translateable_seq` is led by one `N` per base of the translation start exon's
    /// phase (core `Transcript.pm:917-919`): a `cds_start_NF` transcript whose exon 1
    /// opens mid-codon with Ensembl phase 1 (GFF3 phase 2) gets `N` plus the coding
    /// bases, not the `NN` that padding by the raw GFF3 phase column would give; the
    /// mapper's `start_phase` is that same exon-1 phase.
    #[test]
    fn translateable_seq_n_padded_on_cds_start_nf() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = genome(dir.path());
        // Exon 1 at 4-12, coding from its first base: GFF3 phase 2 is Ensembl phase 1.
        let exons = [exon(4, 12, 1, 1), exon(21, 30, 1, 2)];
        let cds_regions = [cds(4, 12, 1, 2), cds(21, 26, 1, 0)];
        let ordered = transcript_order(&exons, 1);
        let phases = exon_phases(&ordered, &cds_regions, 1);
        assert_eq!(phases[0].0, 1);

        let mapper_result = mapper::build_mapper(&exons, &cds_regions, 1, phases[0].0);
        assert_eq!(mapper_result.start_phase, 1);
        assert_eq!(mapper_result.cdna_coding_start, Some(1));
        assert_eq!(mapper_result.cdna_coding_end, Some(15));

        let spliced = spliced_seq("21", &ordered, 1, &fasta).unwrap();
        assert_eq!(spliced, "ATGCCCAAAGCCCTTTAAA");
        let (translateable, five, three) = coding_sequences(&spliced, 1, 15, phases[0].0);
        assert_eq!(translateable, "NATGCCCAAAGCCCTT");
        assert_eq!(five, None);
        assert_eq!(three, Some("TAAA".to_string()));
    }

    /// `five_prime_utr` and `three_prime_utr` are the spliced cDNA before
    /// `cdna_coding_start` and after `cdna_coding_end` (core `Transcript.pm:1746-1797`),
    /// read from the genome and reverse-complemented on the reverse strand, and the
    /// translateable sequence between them carries no padding when the start exon's
    /// phase is -1 or 0. A transcript whose CDS fills its exons has neither UTR.
    #[test]
    fn utr_sequences_written_from_fasta() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = genome(dir.path());

        // Forward: exon 1 at 1-12 with the CDS from 4, exon 2 at 21-30 with the CDS to 26.
        let exons = [exon(1, 12, 1, 1), exon(21, 30, 1, 2)];
        let cds_regions = [cds(4, 12, 1, 0), cds(21, 26, 1, 0)];
        let ordered = transcript_order(&exons, 1);
        let phases = exon_phases(&ordered, &cds_regions, 1);
        let mapper_result = mapper::build_mapper(&exons, &cds_regions, 1, phases[0].0);
        let spliced = spliced_seq("21", &ordered, 1, &fasta).unwrap();
        let (translateable, five, three) = coding_sequences(
            &spliced,
            mapper_result.cdna_coding_start.unwrap(),
            mapper_result.cdna_coding_end.unwrap(),
            phases[0].0,
        );
        assert_eq!(five, Some("GGG".to_string()));
        assert_eq!(translateable, "ATGCCCAAAGCCCTT");
        assert_eq!(three, Some("TAAA".to_string()));

        // Reverse: exon 1 at 51-60 and exon 2 at 31-40 on the minus strand; the CDS
        // runs from 58 down to 35, so the 5' UTR is the reverse complement of 59-60.
        let exons = [exon(51, 60, -1, 1), exon(31, 40, -1, 2)];
        let cds_regions = [cds(51, 58, -1, 0), cds(35, 40, -1, 0)];
        let ordered = transcript_order(&exons, -1);
        let phases = exon_phases(&ordered, &cds_regions, -1);
        let mapper_result = mapper::build_mapper(&exons, &cds_regions, -1, phases[0].0);
        let spliced = spliced_seq("21", &ordered, -1, &fasta).unwrap();
        assert_eq!(spliced, "TTTGGGCTAATCATGGGCCC");
        let (translateable, five, three) = coding_sequences(
            &spliced,
            mapper_result.cdna_coding_start.unwrap(),
            mapper_result.cdna_coding_end.unwrap(),
            phases[0].0,
        );
        assert_eq!(five, Some("TT".to_string()));
        assert_eq!(translateable, "TGGGCTAATCATGG");
        assert_eq!(three, Some("GCCC".to_string()));

        // A CDS covering every exon base leaves both UTRs absent.
        let exons = [exon(4, 12, 1, 1)];
        let ordered = transcript_order(&exons, 1);
        let spliced = spliced_seq("21", &ordered, 1, &fasta).unwrap();
        let (translateable, five, three) = coding_sequences(&spliced, 1, 9, 0);
        assert_eq!(translateable, "ATGCCCAAA");
        assert_eq!((five, three), (None, None));
    }

    /// The translation's `start` and `end` are positions within the start and end
    /// exons (core `Translation.pm:190-226`), not cDNA coordinates: on the forward
    /// strand the CDS opening at the fourth base of exon 1 gives `start` 4; on the
    /// reverse strand the CDS opening three bases below exon 1's end gives 3.
    #[test]
    fn translation_bounds_are_positions_within_the_exons() {
        let exons = [exon(1, 12, 1, 1), exon(21, 30, 1, 2)];
        let cds_regions = [cds(4, 12, 1, 0), cds(21, 26, 1, 0)];
        let ordered = transcript_order(&exons, 1);
        assert_eq!(translation_bounds(&ordered, &cds_regions, 1), Some((4, 6)));

        let exons = [exon(51, 60, -1, 1), exon(31, 40, -1, 2)];
        let cds_regions = [cds(51, 58, -1, 0), cds(35, 40, -1, 0)];
        let ordered = transcript_order(&exons, -1);
        assert_eq!(translation_bounds(&ordered, &cds_regions, -1), Some((3, 6)));

        assert_eq!(translation_bounds(&ordered, &[], -1), None);
    }

    fn transcript(
        chr: &str,
        strand: i8,
        exons: Vec<ExonRecord>,
        cds_regions: Vec<CdsRecord>,
    ) -> crate::gff3::TranscriptRecord {
        crate::gff3::TranscriptRecord {
            transcript_id: "ENST1".to_string(),
            gene_id: "ENSG1".to_string(),
            chr: chr.to_string(),
            start: exons.iter().map(|e| e.start).min().unwrap(),
            end: exons.iter().map(|e| e.end).max().unwrap(),
            strand,
            biotype: if cds_regions.is_empty() {
                "lncRNA"
            } else {
                "protein_coding"
            }
            .to_string(),
            source: "ensembl".to_string(),
            version: Some(1),
            is_canonical: true,
            mane_select: false,
            mane_plus_clinical: false,
            gencode_basic: true,
            gencode_primary: false,
            tsl: None,
            ccds: None,
            exons,
            cds_regions,
        }
    }

    fn gene(tr: &crate::gff3::TranscriptRecord) -> GeneRecord {
        GeneRecord {
            gene_id: tr.gene_id.clone(),
            chr: tr.chr.clone(),
            start: tr.start,
            end: tr.end,
            strand: tr.strand,
            biotype: tr.biotype.clone(),
            name: Some("GENE1".to_string()),
            version: Some(3),
            transcripts: Vec::new(),
        }
    }

    /// A cached record carries the keys Ensembl VEP's dumper writes and none of the
    /// keys it never writes: no exon or intron `strand`, no mapper `pair_count`,
    /// `from_id` or `to_id`; `seq_edits` is always an array; a non-coding transcript's
    /// `translateable_seq` is the empty string (core `Transcript.pm:906-908`), with no
    /// UTR and no translation; `codon_table` on MT is 2.
    #[test]
    fn cached_record_carries_the_dumper_keys_only() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = genome(dir.path());
        let coding = transcript(
            "21",
            1,
            vec![exon(1, 12, 1, 1), exon(21, 30, 1, 2)],
            vec![cds(4, 12, 1, 0), cds(21, 26, 1, 0)],
        );
        let record = build_transcript(
            &gene(&coding),
            &coding,
            &HashMap::new(),
            None,
            Some(&fasta),
            None,
            None,
        );
        let json = serde_json::to_value(&record).unwrap();
        let vefc = &json["variation_effect_feature_cache"];
        assert_eq!(vefc["codon_table"], 1);
        assert_eq!(vefc["translateable_seq"], "ATGCCCAAAGCCCTT");
        assert_eq!(vefc["five_prime_utr"], "GGG");
        assert_eq!(vefc["three_prime_utr"], "TAAA");
        assert_eq!(vefc["seq_edits"], serde_json::json!([]));
        assert_eq!(json["translation"]["start"], 4);
        assert_eq!(json["translation"]["end"], 6);
        assert_eq!(vefc["mapper"]["start_phase"], -1);
        for exon in json["exons"].as_array().unwrap() {
            assert!(exon.get("strand").is_none(), "{exon}");
        }
        for intron in vefc["introns"].as_array().unwrap() {
            assert!(intron.get("strand").is_none(), "{intron}");
        }
        assert!(vefc["mapper"].get("pair_count").is_none());
        for pair in vefc["mapper"]["pairs"].as_array().unwrap() {
            assert!(
                pair.get("from_id").is_none() && pair.get("to_id").is_none(),
                "{pair}"
            );
        }
        assert!(json.get("description").is_none());

        let non_coding = transcript(
            "MT",
            -1,
            vec![exon(31, 40, -1, 1), exon(11, 20, -1, 2)],
            vec![],
        );
        let record = build_transcript(
            &gene(&non_coding),
            &non_coding,
            &HashMap::new(),
            None,
            Some(&fasta),
            None,
            None,
        );
        let json = serde_json::to_value(&record).unwrap();
        let vefc = &json["variation_effect_feature_cache"];
        assert_eq!(vefc["codon_table"], 2);
        assert_eq!(vefc["translateable_seq"], "");
        assert!(vefc.get("five_prime_utr").is_none() && vefc.get("three_prime_utr").is_none());
        assert!(json.get("translation").is_none());
        assert!(vefc.get("mapper").is_some());
    }

    /// The codon table is the sequence region's attribute when the dumps carry one,
    /// else 2 on `MT` and 1 elsewhere (`Database/Transcript.pm:359-363`).
    #[test]
    fn codon_table_2_on_mt() {
        assert_eq!(codon_table("MT", None), 2);
        assert_eq!(codon_table("21", None), 1);
    }

    /// A seq edit reads a `translation_attrib` value as `start end alt_seq`
    /// (`SeqEdit.pm:144`): `_selenocysteine` at residue 7, `initial_met` over residue 1
    /// and an `amino_acid_sub`; a value without two coordinates is no edit.
    #[test]
    fn seq_edit_from_translation_attrib_value() {
        let sec = CacheSeqEdit::from_attrib("_selenocysteine", "7 7 U").unwrap();
        assert_eq!(
            sec,
            CacheSeqEdit {
                code: "_selenocysteine".to_string(),
                name: "Selenocysteine".to_string(),
                start: 7,
                end: 7,
                alt_seq: "U".to_string(),
                length_diff: 0,
            }
        );
        let met = CacheSeqEdit::from_attrib("initial_met", "1 1 M").unwrap();
        assert_eq!(
            (met.start, met.end, met.alt_seq.as_str(), met.length_diff),
            (1, 1, "M", 0)
        );
        let json = serde_json::to_value(&met).unwrap();
        assert_eq!(json["code"], "initial_met");
        assert_eq!(json["alt_seq"], "M");
        assert_eq!(json["length_diff"], 0);
        assert!(CacheSeqEdit::from_attrib("amino_acid_sub", "x").is_none());
    }

    /// A same-named HGNC gene's id is lent once, at build time, from the
    /// transcripts overlapping the region the transcript starts in, so a transcript
    /// prints one `HGNC_ID` whatever the input. Ensembl VEP's `merge_features`
    /// lends it per 1 Mb region file (Dumper/Core.pm:109; AnnotationType/Transcript.pm:290,
    /// 305) and again over every buffer's features (AnnotationSource.pm:139-143),
    /// so its value follows the region file read and the buffer's other records;
    /// `docs/divergences/hgnc-id-copy.md` states the difference.
    #[test]
    fn hgnc_id_is_copied_once_from_the_start_region() {
        fn named(
            id: &str,
            start: u64,
            end: u64,
            symbol: &str,
            hgnc_id: Option<&str>,
        ) -> CacheTranscript {
            let mut tr = transcript("3", 1, vec![exon(start, end, 1, 1)], vec![]);
            tr.transcript_id = id.to_string();
            let mut record =
                build_transcript(&gene(&tr), &tr, &HashMap::new(), None, None, None, None);
            record.gene_symbol = Some(symbol.to_string());
            record.hgnc_id = hgnc_id.map(str::to_string);
            record
        }
        let mut transcripts = vec![
            // The HGNC-named gene, inside the region 185000001-186000000.
            named(
                "HGNC_NAMED",
                185_960_039,
                185_979_837,
                "NMRAL2P",
                Some("HGNC:52332"),
            ),
            // The same symbol, starting in that region and ending in the next.
            named("BOUNDARY", 185_959_973, 186_012_561, "NMRAL2P", None),
            // The same symbol in the next region alone.
            named("NEXT_REGION", 186_005_000, 186_020_000, "NMRAL2P", None),
            // The same symbol, far outside both regions.
            named("FAR", 273_000_000, 273_001_000, "NMRAL2P", None),
            // Another symbol in the HGNC gene's region.
            named("OTHER", 185_970_000, 185_971_000, "OTHER", None),
            // An HGNC-named gene spanning two regions lends its id in both.
            named("WIDE", 190_990_000, 191_010_000, "WIDE", Some("HGNC:9")),
            named("BEFORE_WIDE", 190_500_000, 190_600_000, "WIDE", None),
            named("AFTER_WIDE", 191_500_000, 191_600_000, "WIDE", None),
            named("PAST_WIDE", 192_000_001, 192_000_100, "WIDE", None),
            // Two lenders of one symbol in one region: the later one's id is copied
            // onto every transcript of the symbol, the earlier lender included.
            named("TWO_FIRST", 200_100_000, 200_100_500, "TWO", Some("HGNC:1")),
            named("TWO_LATER", 200_200_000, 200_200_500, "TWO", Some("HGNC:2")),
            named("TWO_BORROWER", 200_300_000, 200_300_500, "TWO", None),
        ];
        copy_hgnc_ids_by_symbol(&mut transcripts);
        let ids: Vec<(&str, Option<&str>)> = transcripts
            .iter()
            .map(|t| (t.stable_id.as_str(), t.hgnc_id.as_deref()))
            .collect();
        assert_eq!(
            ids,
            [
                ("HGNC_NAMED", Some("HGNC:52332")),
                ("BOUNDARY", Some("HGNC:52332")),
                ("NEXT_REGION", None),
                ("FAR", None),
                ("OTHER", None),
                ("WIDE", Some("HGNC:9")),
                ("BEFORE_WIDE", Some("HGNC:9")),
                ("AFTER_WIDE", Some("HGNC:9")),
                ("PAST_WIDE", None),
                ("TWO_FIRST", Some("HGNC:2")),
                ("TWO_LATER", Some("HGNC:2")),
                ("TWO_BORROWER", Some("HGNC:2")),
            ]
        );
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
}
