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

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

use anyhow::{bail, Context, Result};
use serde::Serialize;
use tracing::{info, warn};
use vep_core::codon::{is_start_codon, is_ter_codon, known_table, translate_codon_with_table};
use vep_core::coordinate::is_cache_path_component;

use crate::core_dumps::{CoreDumps, TranslationFacts};
use crate::fasta::{record_of, Peptide};
use crate::gff3::{located, CdsRecord, ExonRecord, GeneRecord};
use crate::mapper;

/// Region size for binning transcripts (1 Mb).
pub(crate) const REGION_SIZE: u64 = 1_000_000;

/// First 1-based position of the 1 Mb region containing `pos`.
pub(crate) fn region_start(pos: u64) -> u64 {
    ((pos - 1) / REGION_SIZE) * REGION_SIZE + 1
}

/// Statistics from the cache build.
pub struct BuildStats {
    pub transcript_count: usize,
    pub chromosome_count: usize,
    pub region_file_count: usize,
    /// The chromosome directories written, in `info::contig_order`. `info.json`
    /// lists them in `info::karyotype_order`: those the core dumps rank in a
    /// karyotype first, by rank.
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
    /// The gene's significant phenotype feature count, 0 when the core dumps
    /// came without the phenotype tables, as the Ensembl VEP dumper stores it
    /// (`Database/Transcript.pm` lines 202, 243-248 and 284); absent when no
    /// core dumps were read.
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
    /// What the build makes of the transcript's protein FASTA record, read
    /// for the build's own counts and never written.
    #[serde(skip)]
    record_use: RecordUse,
}

/// What a build makes of the protein FASTA record that names a transcript.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RecordUse {
    /// No record names the transcript, or no protein FASTA is read.
    Absent,
    /// The record's peptide is the transcript's.
    Taken,
    /// The transcript has a coding sequence and the record's length is not
    /// one it translates to ([`fits_coding_sequence`]): the transcript is
    /// built as one with no record.
    SetAside,
}

/// Whether a protein FASTA record of `residues` can be the peptide of a
/// coding sequence of `coding_bases`, the translateable sequence's length:
/// its whole codons, or one fewer. The sequence carries its stop codon where
/// the model's CDS lines do (an Ensembl GFF3's CDS line, a GTF's `stop_codon`
/// line) and a record never carries the stop's residue, so a complete
/// transcript's record is one fewer; one whose sequence ends without a stop
/// codon has a record of exactly its whole codons.
fn fits_coding_sequence(residues: usize, coding_bases: u64) -> bool {
    let codons = coding_bases / 3;
    let residues = residues as u64;
    residues == codons || residues + 1 == codons
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

/// What a build reads beside its transcript model: the inputs a build may
/// lack, each `None` when the build was given no such file.
#[derive(Clone, Copy)]
pub struct BuildSources<'a> {
    /// The protein FASTA's records, each under its header's `transcript:`
    /// token whole ([`record_of`] finds a transcript's): the peptide as
    /// Ensembl VEP's cache carries it and the header's protein identifier and
    /// version; `None` when no protein FASTA was read and each peptide is
    /// translated from the genome.
    pub protein_sequences: Option<&'a HashMap<String, Peptide>>,
    /// The prediction matrices of each protein FASTA record, under the
    /// record's key in `protein_sequences`.
    pub prediction_data: Option<&'a HashMap<String, crate::predictions::PeptidePredictions>>,
    pub genome_fasta: Option<&'a vep_fasta::IndexedFasta>,
    pub gtf_tags: Option<&'a crate::gtf::GtfTags>,
    pub core_dumps: Option<&'a CoreDumps>,
}

/// Whether a build reading `protein_sequences` and `core_dumps` is the Ensembl
/// cache's content rather than what Ensembl VEP's `--gff` or `--gtf` source
/// annotates: the source builds its transcripts from the model and a genome
/// FASTA alone (`BaseGXF.pm` `new`, `_add_translation`), so a build reading
/// the release's protein FASTA or core table dumps, inputs the source takes
/// no part of, is the cache's.
pub fn is_cache_content(
    protein_sequences: Option<&HashMap<String, Peptide>>,
    core_dumps: Option<&CoreDumps>,
) -> bool {
    protein_sequences.is_some() || core_dumps.is_some()
}

/// The transcripts of a model that Ensembl VEP's cache holds none of, counted
/// by reason: those of biotype `artifact`, and those of another biotype
/// carrying the `readthrough_tra` attribute, on their line or in the core
/// dumps (`Database/Transcript.pm` `get_features_by_regions_uncached`, the
/// dumper's source). Ensembl VEP's `--gff` and `--gtf` sources build both
/// (`BaseGXF.pm` `_create_transcript`), so a build of the cache's content
/// ([`is_cache_content`]) alone sets them aside. The GRCh37 release-87 GFF3
/// carries no readthrough tag; its dumps do.
#[derive(Default)]
pub struct SetAside {
    pub artifact: usize,
    pub readthrough: usize,
}

impl SetAside {
    /// The count over every transcript of `genes`.
    pub fn of(genes: &HashMap<String, GeneRecord>, core_dumps: Option<&CoreDumps>) -> Self {
        let mut set_aside = Self::default();
        for tr in genes.values().flat_map(|gene| &gene.transcripts) {
            set_aside.counts(tr, core_dumps);
        }
        set_aside
    }

    /// Whether the dumper leaves `tr` out, counting it under its reason when
    /// it does.
    fn counts(
        &mut self,
        tr: &crate::gff3::TranscriptRecord,
        core_dumps: Option<&CoreDumps>,
    ) -> bool {
        if tr.biotype == "artifact" {
            self.artifact += 1;
        } else if tr.readthrough
            || core_dumps
                .and_then(|d| d.transcript(&tr.transcript_id))
                .is_some_and(|f| f.readthrough)
        {
            self.readthrough += 1;
        } else {
            return false;
        }
        true
    }

    /// Every transcript counted, under either reason.
    pub fn total(&self) -> usize {
        self.artifact + self.readthrough
    }
}

/// The version of a translation whose first CDS line in transcript order
/// carries `cds`: that `version` when the line has one; else, with the core
/// table dumps read and holding the translation, theirs, present or absent;
/// else the version the protein FASTA header's identifier carries; else none,
/// as Ensembl VEP's cache holds a translation whose identifier the release
/// writes unversioned. The 1 Ensembl VEP's `--gff` and `--gtf` sources give a
/// line without one is set at the parse ([`crate::gff3::Reading::GffSource`]).
fn translation_version(
    cds: Option<u32>,
    dumps: Option<&TranslationFacts>,
    header: Option<u32>,
) -> Option<u32> {
    cds.or(match dumps {
        Some(facts) => facts.version,
        None => header,
    })
}

/// Build the full JSON cache directory from parsed GFF3 genes and `sources`;
/// `codon_tables` is [`codon_tables`]'s map over the same genes. `model` is
/// the file the genes were read from: a transcript whose coding sequence does
/// not sit in its exons ([`Unplaced`]) fails the build before anything is
/// written, naming that file, the CDS record's line and the transcript.
pub fn build_cache(
    genes: &HashMap<String, GeneRecord>,
    model: &Path,
    output_dir: &Path,
    sources: &BuildSources,
    codon_tables: &HashMap<String, u8>,
) -> Result<BuildStats> {
    // Name order: the chromosomes are logged, and names that tie in
    // `contig_order` listed, in one order on every run.
    let mut by_chr: BTreeMap<String, Vec<CacheTranscript>> = BTreeMap::new();
    let mut transcript_count = 0usize;
    let mut set_aside = SetAside::default();
    let mut translated = 0usize;
    let mut without_peptide_record = 0usize;
    let mut records_set_aside = 0usize;

    let cache_content = is_cache_content(sources.protein_sequences, sources.core_dumps);
    // Key order: of several transcripts the build refuses, the same one is
    // named on every run.
    let mut gene_ids: Vec<&String> = genes.keys().collect();
    gene_ids.sort_unstable();
    for gene in gene_ids.into_iter().map(|id| &genes[id]) {
        for tr in &gene.transcripts {
            if cache_content && set_aside.counts(tr, sources.core_dumps) {
                continue;
            }
            let codon_table = codon_tables
                .get(&tr.chr)
                .copied()
                .expect("codon_tables settles every sequence region of genes");
            let cache_tr = match build_transcript(gene, tr, sources, codon_table) {
                Ok(cache_tr) => cache_tr,
                Err(unplaced) => bail!("{}", unplaced.refusal(model, tr)),
            };
            if cache_tr.translation.is_some() {
                translated += 1;
                if sources.protein_sequences.is_some() && cache_tr.record_use == RecordUse::Absent {
                    without_peptide_record += 1;
                }
            }
            if cache_tr.record_use == RecordUse::SetAside {
                records_set_aside += 1;
            }
            by_chr.entry(tr.chr.clone()).or_default().push(cache_tr);
            transcript_count += 1;
        }
    }
    if set_aside.total() > 0 {
        info!(
            "Dropped {} artifact or readthrough transcripts, as Ensembl VEP's dumper does",
            set_aside.total()
        );
    }
    // A protein FASTA of another release or assembly, or one whose headers
    // name no transcript, leaves the cache's content to the other sources, as
    // does a record that is not its transcript's peptide.
    if without_peptide_record > 0 || records_set_aside > 0 {
        warn!(
            "{without_peptide_record} of {translated} translated transcripts have no record in \
             the protein FASTA, and {records_set_aside} records were set aside for their length, \
             neither the coding sequence's whole codons nor one fewer: their peptides are \
             translated from the genome when one is given, else absent, and their versions come \
             from the model or the core dumps, else none"
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

/// Build a single CacheTranscript from parsed GFF3 data; `codon_table` is the
/// transcript's sequence region's, as [`codon_table`] settles it. A transcript
/// whose coding sequence does not sit in its exons ([`Unplaced`]) builds no
/// record.
fn build_transcript(
    gene: &GeneRecord,
    tr: &crate::gff3::TranscriptRecord,
    sources: &BuildSources,
    codon_table: u8,
) -> Result<CacheTranscript, Unplaced> {
    if let Some(outside) = translation_end_outside_exons(&tr.exons, &tr.cds_regions, tr.strand) {
        return Err(outside);
    }
    let BuildSources {
        protein_sequences,
        prediction_data,
        genome_fasta,
        gtf_tags,
        core_dumps,
    } = *sources;
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

    // The padding follows the translation start exon's phase (core
    // Transcript.pm:917-919), which is exon 1's only when exon 1 codes.
    let span = cds_span(&tr.cds_regions);
    let start_exon_phase = ordered
        .iter()
        .zip(&phases)
        .find(|(e, _)| span.is_some_and(|(s, end)| e.start <= end && e.end >= s))
        .map_or(-1, |(_, p)| p.0);
    let coding = match (
        genome_fasta,
        mapper_result.cdna_coding_start,
        mapper_result.cdna_coding_end,
    ) {
        (Some(fasta), Some(cdna_start), Some(cdna_end)) => {
            spliced_seq(&tr.chr, &ordered, strand, fasta)
                .map(|spliced| coding_sequences(&spliced, cdna_start, cdna_end, start_exon_phase))
        }
        _ => None,
    };
    // The translateable sequence's length, known with or without a genome:
    // the cDNA coding span led by the padding; none without a coding span.
    let coding_bases = mapper_result
        .cdna_coding_start
        .zip(mapper_result.cdna_coding_end)
        .map(|(cdna_start, cdna_end)| {
            (cdna_end + 1).saturating_sub(cdna_start) + start_exon_phase.max(0) as u64
        });
    // A transcript without a translation has the empty translateable sequence
    // (core Transcript.pm:906-908), whether or not a genome is given.
    let (translateable_seq, five_prime_utr, three_prime_utr) = match coding {
        Some((t, five, three)) => (Some(t), five, three),
        None if tr.cds_regions.is_empty() => (Some(String::new()), None, None),
        None => (None, None, None),
    };

    // The protein FASTA's peptide carries the translation's sequence edits;
    // without one, or with one whose length the coding sequence cannot
    // translate to, the peptide is the translateable sequence's own
    // translation, as Ensembl VEP reads it off a GFF3 or GTF
    // (`BaseGXF.pm` `_add_translation`, `$translation->seq`).
    let named =
        protein_sequences.and_then(|peptides| record_of(peptides, &tr.transcript_id, tr.version));
    let (record_key, fasta_record) = named
        .filter(|(_, record)| {
            coding_bases.is_none_or(|bases| fits_coding_sequence(record.sequence.len(), bases))
        })
        .unzip();
    let record_use = match (named, fasta_record) {
        (None, _) => RecordUse::Absent,
        (Some(_), Some(_)) => RecordUse::Taken,
        (Some(_), None) => RecordUse::SetAside,
    };
    let peptide = fasta_record.map(|p| p.sequence.clone()).or_else(|| {
        translateable_seq
            .as_deref()
            .filter(|s| !s.is_empty())
            .and_then(|s| translate(s, codon_table))
    });

    // The translation's stable identifier and version start from the
    // `protein_id` and `version` of its first CDS line in transcript order,
    // the segment of lowest start on the forward strand and of highest start
    // on the reverse (`BaseGXF.pm` `_add_translation`, the first of the
    // ordered CDS records); the other sources fill in what the line lacks.
    let first_cds = if strand > 0 {
        tr.cds_regions.iter().min_by_key(|c| c.start)
    } else {
        tr.cds_regions.iter().max_by_key(|c| c.start)
    };
    let protein_id = first_cds
        .and_then(|c| c.protein_id.clone())
        .or_else(|| translation_facts.map(|t| t.stable_id.clone()));
    let protein_version = translation_version(
        first_cds.and_then(|c| c.protein_version),
        translation_facts,
        fasta_record.and_then(|p| p.version_beside(protein_id.as_deref())),
    );

    let translation = translation_bounds(&ordered, &tr.cds_regions, strand)?.map(|(start, end)| {
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
            Some(g.phenotype_count),
        ),
        None => (tr.symbol.clone(), None, None, None),
    };

    let seq_edits: Vec<CacheSeqEdit> = translation_facts
        .map(|t| {
            t.seq_edits
                .iter()
                .filter_map(|e| CacheSeqEdit::from_attrib(&e.code, &e.value))
                .collect()
        })
        .unwrap_or_default();

    Ok(CacheTranscript {
        stable_id: tr.transcript_id.clone(),
        version: tr.version,
        db_id: None,
        gene_stable_id: tr.gene_id.clone(),
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
            codon_table,
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
                .zip(record_key)
                .and_then(|(predictions, key)| predictions.get(key))
                .map(convert_predictions),
        },
        record_use,
    })
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
        for code in gtf.kept_attribute_codes(&tr.transcript_id, tr.version) {
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
                // The segment at the exon's 5' end: a stop codon written as
                // its own segment sits 3' of the exon's CDS segment.
                let mut overlapping = cds_regions
                    .iter()
                    .filter(|c| c.start <= exon.end && c.end >= exon.start);
                let gff3_phase = if strand < 0 {
                    overlapping.next_back()
                } else {
                    overlapping.next()
                }
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

/// A coding sequence that does not sit in its transcript's exons, one variant
/// for each shape a build refuses.
#[derive(Debug, PartialEq, Eq)]
enum Unplaced {
    /// An end of the translation that lies in no exon of the transcript
    /// ([`translation_end_outside_exons`]): `point` is the coding sequence's
    /// 5' base when `starts` holds and its 3' base otherwise, and `line` and
    /// `cds` are those of the CDS record that sets it, the lowest line when
    /// several do.
    End {
        starts: bool,
        point: u64,
        line: usize,
        cds: (u64, u64),
    },
    /// A coding sequence that opens outside the first exon it overlaps in
    /// transcript order ([`translation_bounds`]), as one does whose records
    /// each lie inside an exon while the exons' ranks run against the strand:
    /// `cds` is the whole coding span and `line` the line of the record that
    /// opens it, the lowest when several open it at one base.
    Opening {
        line: usize,
        cds: (u64, u64),
        exon: (u64, u64),
    },
    /// An end of the translation whose offset lies outside the exon it is
    /// counted in ([`translation_bounds`]), as the 3' end's does where the
    /// exons' ranks run against position there or an exon nests in the last
    /// coding one and comes after it in transcript order: `offset` is the
    /// position the translation would carry for the end, `exon` the exon it
    /// is counted in, and `line` and `cds` are those of the CDS record that
    /// sets the end, the lowest line when several do.
    Offset {
        starts: bool,
        offset: u64,
        exon: (u64, u64),
        line: usize,
        cds: (u64, u64),
    },
}

impl Unplaced {
    /// The refusal `tr`, read from the file `model`, draws: the file and the
    /// CDS record's line, the transcript, and either the record's span with
    /// the end of the translation it sets, in no exon or at an offset outside
    /// the exon it is counted in, or the coding sequence's span with the
    /// exon's.
    fn refusal(&self, model: &Path, tr: &crate::gff3::TranscriptRecord) -> String {
        let id = &tr.transcript_id;
        match *self {
            Unplaced::End {
                starts,
                point,
                line,
                cds: (start, end),
            } => format!(
                "{}: CDS {start}-{end} of transcript {id} {} the translation at {point}, in no \
                 exon of the transcript",
                located(model, line),
                if starts { "starts" } else { "ends" },
            ),
            Unplaced::Opening {
                line,
                cds: (start, end),
                exon: (exon_start, exon_end),
            } => format!(
                "{}: the coding sequence of transcript {id}, {start}-{end}, {} \
                 {exon_start}-{exon_end}, the first exon in transcript order that it overlaps",
                located(model, line),
                if tr.strand < 0 {
                    "ends after"
                } else {
                    "starts before"
                },
            ),
            Unplaced::Offset {
                starts,
                offset,
                exon: (exon_start, exon_end),
                line,
                cds: (start, end),
            } => format!(
                "{}: CDS {start}-{end} of transcript {id} {} the translation at base {offset} of \
                 the {}-base exon {exon_start}-{exon_end}, the {} in transcript order that the \
                 coding sequence overlaps",
                located(model, line),
                if starts { "starts" } else { "ends" },
                exon_end - exon_start + 1,
                if starts { "first" } else { "last" },
            ),
        }
    }
}

/// The end of the translation that lies in no exon of `exons`, its start
/// looked at before its end: the coding sequence's 5' base and its 3' base
/// each lie inside an exon of the transcript. A CDS record's other end, and a
/// record that sets neither base, are not held to the exons.
fn translation_end_outside_exons(
    exons: &[ExonRecord],
    cds_regions: &[CdsRecord],
    strand: i8,
) -> Option<Unplaced> {
    let (lowest, highest) = cds_span(cds_regions)?;
    [true, false].into_iter().find_map(|starts| {
        // The translation starts at the span's lowest base on the forward
        // strand and at its highest on the reverse.
        let at_lowest = starts == (strand >= 0);
        let point = if at_lowest { lowest } else { highest };
        if exons.iter().any(|e| e.start <= point && point <= e.end) {
            return None;
        }
        let sets_it = |c: &&CdsRecord| point == if at_lowest { c.start } else { c.end };
        let record = cds_regions
            .iter()
            .filter(sets_it)
            .min_by_key(|c| c.line)
            .expect("the span's ends are those of CDS records");
        Some(Unplaced::End {
            starts,
            point,
            line: record.line,
            cds: (record.start, record.end),
        })
    })
}

/// The translation's `start` and `end`: the first coding base's position within the
/// translation start exon and the last coding base's within the end exon, each
/// 1-based in transcript orientation (core `Translation.pm:190-226`); `None` for a
/// non-coding transcript. The start exon is the first in transcript order that the
/// coding sequence overlaps, and a coding sequence that opens outside it, before
/// its start on the forward strand or after its end on the reverse, has no
/// position within it and is refused ([`Unplaced::Opening`]). The end exon is
/// the last in transcript order that the coding sequence overlaps, and neither
/// position is returned outside the exon it is counted in: one that would be,
/// as the end is when that exon does not hold the last coding base, is refused
/// ([`Unplaced::Offset`]).
fn translation_bounds(
    ordered: &[&ExonRecord],
    cds_regions: &[CdsRecord],
    strand: i8,
) -> Result<Option<(u64, u64)>, Unplaced> {
    let Some((cds_start, cds_end)) = cds_span(cds_regions) else {
        return Ok(None);
    };
    let overlaps = |e: &&&ExonRecord| e.start <= cds_end && e.end >= cds_start;
    let (Some(start_exon), Some(end_exon)) = (
        ordered.iter().find(overlaps),
        ordered.iter().rev().find(overlaps),
    ) else {
        return Ok(None);
    };
    let (before_start, end) = if strand < 0 {
        (
            start_exon.end.checked_sub(cds_end),
            end_exon.end - cds_start + 1,
        )
    } else {
        (
            cds_start.checked_sub(start_exon.start),
            cds_end - end_exon.start + 1,
        )
    };
    // The CDS record that sets an end of the coding sequence, its 5' end when
    // `starts` holds: the lowest line of those that share the base.
    let sets = |starts: bool| {
        cds_regions
            .iter()
            .filter(|c| {
                if starts == (strand < 0) {
                    c.end == cds_end
                } else {
                    c.start == cds_start
                }
            })
            .min_by_key(|c| c.line)
            .expect("the span's ends are those of CDS records")
    };
    let Some(before_start) = before_start else {
        return Err(Unplaced::Opening {
            cds: (cds_start, cds_end),
            exon: (start_exon.start, start_exon.end),
            line: sets(true).line,
        });
    };
    let start = before_start + 1;
    for (starts, offset, exon) in [(true, start, start_exon), (false, end, end_exon)] {
        let length = exon.end - exon.start + 1;
        if !(1..=length).contains(&offset) {
            let record = sets(starts);
            return Err(Unplaced::Offset {
                starts,
                offset,
                exon: (exon.start, exon.end),
                line: record.line,
                cds: (record.start, record.end),
            });
        }
    }
    Ok(Some((start, end)))
}

/// The peptide `Transcript::translate` reads off a translateable sequence under
/// `codon_table` (core `Transcript.pm:2092-2176`): the sequence cut to whole
/// codons; each codon's residue, an ambiguous codon resolving when every
/// expansion of its IUPAC codes agrees, to `B` for Asp or Asn and `Z` for Glu
/// or Gln, else to `X` (`Bio::Tools::CodonTable::translate`); the terminal
/// codon's residue dropped when it is a stop codon of the table
/// ([`is_ter_codon`]); and the first residue `M` when the first codon is a
/// start codon of the table ([`is_start_codon`]). `None` for a sequence
/// shorter than one codon.
fn translate(seq: &str, codon_table: u8) -> Option<String> {
    let bytes = seq.as_bytes();
    let whole = &bytes[..bytes.len() - bytes.len() % 3];
    if whole.is_empty() {
        return None;
    }
    let mut peptide: Vec<u8> = whole
        .as_chunks::<3>()
        .0
        .iter()
        .map(|codon| translate_codon(codon, codon_table))
        .collect();
    if is_ter_codon(&whole[whole.len() - 3..], codon_table) {
        peptide.pop();
    }
    if peptide.first().is_some_and(|&aa| aa != b'M') && is_start_codon(&whole[..3], codon_table) {
        peptide[0] = b'M';
    }
    Some(String::from_utf8(peptide).expect("residues are ASCII"))
}

fn is_unambiguous(codon: &[u8]) -> bool {
    codon
        .iter()
        .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

/// The bases an IUPAC code stands for; `U` is `T`.
fn iupac_bases(code: u8) -> &'static [u8] {
    match code.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => b"",
    }
}

/// One codon's residue as `Bio::Tools::CodonTable::translate` gives it.
fn translate_codon(codon: &[u8], codon_table: u8) -> u8 {
    if is_unambiguous(codon) {
        return translate_codon_with_table(codon, codon_table);
    }
    let mut residues: Vec<u8> = Vec::new();
    for &a in iupac_bases(codon[0]) {
        for &b in iupac_bases(codon[1]) {
            for &c in iupac_bases(codon[2]) {
                let aa = translate_codon_with_table(&[a, b, c], codon_table);
                if !residues.contains(&aa) {
                    residues.push(aa);
                }
            }
        }
    }
    residues.sort_unstable();
    match residues.as_slice() {
        [aa] => *aa,
        b"DN" => b'B',
        b"EQ" => b'Z',
        _ => b'X',
    }
}

/// Names an organelle genome goes by in the assemblies Ensembl serves and in
/// GFF3 or GTF files from other sources, compared in any letter case; a region
/// so named translates under a table of its own, which nothing but the core
/// dumps or `--codon-table` can state.
const ORGANELLE_REGION_NAMES: [&str; 12] = [
    "MT",
    "chrM",
    "chrMT",
    "Mito",
    "MtDNA",
    "mitochondrion_genome",
    "Mt",
    "Pt",
    "ChrC",
    "ChrM",
    "chloroplast",
    "plastid",
];

fn is_organelle_name(chr: &str) -> bool {
    ORGANELLE_REGION_NAMES
        .iter()
        .any(|name| name.eq_ignore_ascii_case(chr))
}

/// The `--codon-table` values as a map of sequence region name to table:
/// each value `<region>=<table>`, the region named as the GFF3 or GTF spells
/// it (`chrM` and `chrMT` read as `MT`, as the model readers read them), the
/// table one of the genetic codes `Bio::Tools::CodonTable` ships
/// ([`vep_core::codon::known_table`]). A region named twice with two tables,
/// a table the module lacks, an empty region or a value without `=` is refused.
pub fn codon_table_flag(values: &[String]) -> Result<HashMap<String, u8>> {
    let mut tables: HashMap<String, u8> = HashMap::new();
    for value in values {
        let Some((region, table_text)) = value.split_once('=') else {
            bail!("--codon-table {value} is not <region>=<table>");
        };
        let region = crate::gff3::normalize_chr(region.trim());
        if region.is_empty() {
            bail!("--codon-table {value} names no sequence region");
        }
        let table: u8 = match table_text.trim().parse::<u8>() {
            Ok(table) if known_table(table) => table,
            _ => bail!(
                "--codon-table {value}: {table_text} is not a genetic code Ensembl VEP translates with"
            ),
        };
        if let Some(&earlier) = tables.get(region) {
            if earlier != table {
                bail!(
                    "--codon-table names sequence region {region} twice, as {earlier} and {table}"
                );
            }
        }
        tables.insert(region.to_string(), table);
    }
    Ok(tables)
}

/// The codon table of every sequence region a transcript of `genes` sits on,
/// by region name, each settled by [`codon_table`]; a `--codon-table` region
/// no transcript sits on is logged. Resolved before any shard is built, so a
/// region whose table is unknown refuses the build before it starts.
pub fn codon_tables(
    genes: &HashMap<String, GeneRecord>,
    core_dumps: Option<&CoreDumps>,
    named: &HashMap<String, u8>,
) -> Result<HashMap<String, u8>> {
    // Name order, so two unknown regions refuse on the same one every run.
    let regions: BTreeSet<&str> = genes
        .values()
        .flat_map(|g| &g.transcripts)
        .map(|tr| tr.chr.as_str())
        .collect();
    let mut tables: HashMap<String, u8> = HashMap::new();
    for chr in regions {
        tables.insert(chr.to_string(), codon_table(chr, core_dumps, named)?);
    }
    let unused: BTreeSet<&str> = named
        .keys()
        .map(String::as_str)
        .filter(|r| !tables.contains_key(*r))
        .collect();
    for region in unused {
        warn!("--codon-table names {region}, a sequence region no transcript sits on");
    }
    Ok(tables)
}

/// The codon table of one sequence region: its `codon_table` attribute when
/// the core dumps carry one for the region, else the `--codon-table` value
/// `named` for it, else 1. Ensembl VEP reads the slice attribute and defaults
/// to 1 (`Database/Transcript.pm:359-364`); on both human assemblies the
/// mitochondrion is the only region carrying the attribute, with value 2.
///
/// Refused: a region the dumps and the flag give two different tables, and a
/// region named like an organelle genome ([`ORGANELLE_REGION_NAMES`]) that
/// neither names, since the default of 1 would translate its transcripts under
/// the nuclear code.
fn codon_table(
    chr: &str,
    core_dumps: Option<&CoreDumps>,
    named: &HashMap<String, u8>,
) -> Result<u8> {
    let from_dumps = core_dumps.and_then(|d| d.codon_table(chr));
    let from_flag = named.get(chr).copied();
    match (from_dumps, from_flag) {
        (Some(dumps), Some(flag)) if dumps != flag => bail!(
            "codon table of sequence region {chr} is {dumps} in the core dumps and {flag} in --codon-table"
        ),
        (Some(table), _) | (None, Some(table)) => Ok(table),
        (None, None) if is_organelle_name(chr) => {
            if core_dumps.is_some() {
                bail!(
                    "codon table of sequence region {chr} unknown: the core dumps carry no codon_table attribute for it; pass --codon-table {chr}=<n>"
                );
            }
            bail!("codon table of sequence region {chr} unknown: pass --core-dumps or --codon-table {chr}=<n>")
        }
        (None, None) => Ok(1),
    }
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
            line: 1,
        }
    }

    /// A protein FASTA map of one record under `id`, the header's identifier
    /// written without a version.
    fn peptides(id: &str, sequence: &str) -> HashMap<String, Peptide> {
        HashMap::from([(
            id.to_string(),
            Peptide {
                sequence: sequence.to_string(),
                identifier: "ENSP1".to_string(),
                version: None,
            },
        )])
    }

    /// A build's sources with the peptides and genome given and nothing else;
    /// an empty peptide map stands for a build without a protein FASTA.
    fn sources<'a>(
        proteins: &'a HashMap<String, Peptide>,
        genome: Option<&'a vep_fasta::IndexedFasta>,
    ) -> BuildSources<'a> {
        BuildSources {
            protein_sequences: (!proteins.is_empty()).then_some(proteins),
            prediction_data: None,
            genome_fasta: genome,
            gtf_tags: None,
            core_dumps: None,
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
    /// reverse strand the CDS opening at the third base from exon 1's end gives 3.
    #[test]
    fn translation_bounds_are_positions_within_the_exons() {
        let exons = [exon(1, 12, 1, 1), exon(21, 30, 1, 2)];
        let cds_regions = [cds(4, 12, 1, 0), cds(21, 26, 1, 0)];
        let ordered = transcript_order(&exons, 1);
        assert_eq!(
            translation_bounds(&ordered, &cds_regions, 1),
            Ok(Some((4, 6)))
        );

        let exons = [exon(51, 60, -1, 1), exon(31, 40, -1, 2)];
        let cds_regions = [cds(51, 58, -1, 0), cds(35, 40, -1, 0)];
        let ordered = transcript_order(&exons, -1);
        assert_eq!(
            translation_bounds(&ordered, &cds_regions, -1),
            Ok(Some((3, 6)))
        );

        assert_eq!(translation_bounds(&ordered, &[], -1), Ok(None));
    }

    /// A coding sequence opens within the first exon in transcript order that
    /// it overlaps, at that exon's first base at the earliest. One that opens
    /// before the exon's start on the forward strand, or after its end on the
    /// reverse, is refused with both spans and the line of the CDS record that
    /// opens it, the lowest of two lines sharing the opening base, whether the
    /// exon is the transcript's first or a later one.
    #[test]
    fn a_coding_sequence_that_opens_outside_its_start_exon_is_refused() {
        let at_line = |line, record: CdsRecord| CdsRecord { line, ..record };

        let exons = [exon(101, 200, 1, 1), exon(301, 400, 1, 2)];
        let ordered = transcript_order(&exons, 1);
        assert_eq!(
            translation_bounds(&ordered, &[cds(101, 150, 1, 0)], 1),
            Ok(Some((1, 50)))
        );
        let before_exon_one = [
            at_line(7, cds(100, 150, 1, 0)),
            at_line(8, cds(301, 350, 1, 0)),
        ];
        assert_eq!(
            translation_bounds(&ordered, &before_exon_one, 1),
            Err(Unplaced::Opening {
                cds: (100, 350),
                exon: (101, 200),
                line: 7,
            })
        );
        let in_the_intron = [
            at_line(9, cds(250, 350, 1, 0)),
            at_line(5, cds(250, 320, 1, 0)),
        ];
        assert_eq!(
            translation_bounds(&ordered, &in_the_intron, 1),
            Err(Unplaced::Opening {
                cds: (250, 350),
                exon: (301, 400),
                line: 5,
            })
        );

        let exons = [exon(301, 400, -1, 1), exon(101, 200, -1, 2)];
        let ordered = transcript_order(&exons, -1);
        assert_eq!(
            translation_bounds(&ordered, &[cds(150, 200, -1, 0)], -1),
            Ok(Some((1, 51)))
        );
        assert_eq!(
            translation_bounds(&ordered, &[at_line(4, cds(150, 250, -1, 0))], -1),
            Err(Unplaced::Opening {
                cds: (150, 250),
                exon: (101, 200),
                line: 4,
            })
        );
        let after_exon_one = [
            at_line(3, cds(350, 401, -1, 0)),
            at_line(6, cds(150, 200, -1, 0)),
        ];
        assert_eq!(
            translation_bounds(&ordered, &after_exon_one, -1),
            Err(Unplaced::Opening {
                cds: (150, 401),
                exon: (301, 400),
                line: 3,
            })
        );
    }

    /// Neither end of the translation is counted outside its exon. Where the
    /// exons' ranks run against position at the translation's end, or an
    /// exon nests in the last coding one and is ranked after it, the last
    /// exon in transcript order that the coding sequence overlaps does not
    /// hold its last base: the bounds are refused with the offset, that exon
    /// and the CDS record that sets the end, on either strand. The same
    /// exons ranked by position give bounds. The start's refusal reads as
    /// the end's does.
    #[test]
    fn a_translation_offset_outside_its_exon_is_refused() {
        let at_line = |line, record: CdsRecord| CdsRecord { line, ..record };
        let past_its_exon = |offset, exon, line, cds| {
            Err(Unplaced::Offset {
                starts: false,
                offset,
                exon,
                line,
                cds,
            })
        };

        let coding = |strand| {
            [
                at_line(7, cds(1100, 1200, strand, 0)),
                at_line(8, cds(1400, 1500, strand, 0)),
                at_line(9, cds(1700, 1800, strand, 0)),
            ]
        };
        let ranked = |strand, ranks: [u32; 3]| {
            [
                exon(1000, 1200, strand, ranks[0]),
                exon(1400, 1500, strand, ranks[1]),
                exon(1700, 2000, strand, ranks[2]),
            ]
        };
        let exons = ranked(1, [1, 3, 2]);
        assert_eq!(
            translation_bounds(&transcript_order(&exons, 1), &coding(1), 1),
            past_its_exon(401, (1400, 1500), 9, (1700, 1800))
        );
        let exons = ranked(1, [1, 2, 3]);
        assert_eq!(
            translation_bounds(&transcript_order(&exons, 1), &coding(1), 1),
            Ok(Some((101, 101)))
        );
        let exons = ranked(-1, [2, 3, 1]);
        assert_eq!(
            translation_bounds(&transcript_order(&exons, -1), &coding(-1), -1),
            past_its_exon(401, (1400, 1500), 7, (1100, 1200))
        );
        let exons = ranked(-1, [3, 2, 1]);
        assert_eq!(
            translation_bounds(&transcript_order(&exons, -1), &coding(-1), -1),
            Ok(Some((201, 101)))
        );

        let across = |strand| [at_line(6, cds(1100, 1450, strand, 0))];
        let nested = [exon(1000, 1500, 1, 1), exon(1300, 1350, 1, 2)];
        assert_eq!(
            translation_bounds(&transcript_order(&nested, 1), &across(1), 1),
            past_its_exon(151, (1300, 1350), 6, (1100, 1450))
        );
        let nested = [exon(1000, 1500, -1, 1), exon(1300, 1350, -1, 2)];
        assert_eq!(
            translation_bounds(&transcript_order(&nested, -1), &across(-1), -1),
            past_its_exon(251, (1300, 1350), 6, (1100, 1450))
        );

        let tr = transcript("21", 1, vec![exon(11, 15, 1, 1)], Vec::new());
        let start_past_its_exon = Unplaced::Offset {
            starts: true,
            offset: 7,
            exon: (11, 15),
            line: 3,
            cds: (17, 20),
        };
        assert_eq!(
            start_past_its_exon.refusal(Path::new("model.gff3"), &tr),
            "model.gff3 line 3: CDS 17-20 of transcript ENST1 starts the translation at base 7 \
             of the 5-base exon 11-15, the first in transcript order that the coding sequence \
             overlaps"
        );
    }

    /// The translation's start and end each lie inside an exon of the
    /// transcript: the coding sequence's lowest and highest bases on the
    /// forward strand, its highest and lowest on the reverse. An end in no
    /// exon is named with the CDS record that sets it, the lowest line when
    /// two set it, and the start when both ends lie in none. What lies
    /// between the two ends is not looked at: a record in the intron, a
    /// record across it, a record's other end.
    #[test]
    fn an_end_of_the_translation_in_no_exon_is_named() {
        let at_line = |line, record: CdsRecord| CdsRecord { line, ..record };
        let end = |starts, point, line, cds| {
            Some(Unplaced::End {
                starts,
                point,
                line,
                cds,
            })
        };
        for strand in [1, -1] {
            let forward = strand > 0;
            let exons = [exon(101, 200, strand, 1), exon(301, 400, strand, 2)];
            let named =
                |records: &[CdsRecord]| translation_end_outside_exons(&exons, records, strand);
            assert_eq!(named(&[]), None);
            let sound = [cds(150, 200, strand, 0), cds(301, 350, strand, 0)];
            assert_eq!(named(&sound), None);

            let lowest_before_its_exon = [
                at_line(4, cds(100, 200, strand, 0)),
                at_line(5, cds(301, 350, strand, 0)),
            ];
            assert_eq!(
                named(&lowest_before_its_exon),
                end(forward, 100, 4, (100, 200))
            );
            let highest_past_its_exon = [
                at_line(4, cds(150, 200, strand, 0)),
                at_line(5, cds(301, 401, strand, 0)),
            ];
            assert_eq!(
                named(&highest_past_its_exon),
                end(!forward, 401, 5, (301, 401))
            );
            let in_no_exon = [at_line(7, cds(210, 290, strand, 0))];
            let start = if forward { 210 } else { 290 };
            assert_eq!(named(&in_no_exon), end(true, start, 7, (210, 290)));
            let two_at_one_base = [
                at_line(9, cds(250, 350, strand, 0)),
                at_line(6, cds(250, 320, strand, 0)),
            ];
            assert_eq!(named(&two_at_one_base), end(forward, 250, 6, (250, 320)));

            let one_more_in_the_intron = [
                cds(150, 200, strand, 0),
                cds(220, 280, strand, 0),
                cds(301, 350, strand, 0),
            ];
            assert_eq!(named(&one_more_in_the_intron), None);
            assert_eq!(named(&[cds(150, 350, strand, 0)]), None);
            let other_end_in_the_intron = [cds(150, 250, strand, 0), cds(301, 350, strand, 0)];
            assert_eq!(named(&other_end_in_the_intron), None);
        }
        assert_eq!(
            translation_end_outside_exons(&[], &[at_line(3, cds(10, 20, 1, 0))], 1),
            end(true, 10, 3, (10, 20))
        );
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
            symbol: Some("GENE1".to_string()),
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
            readthrough: false,
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
    /// UTR and no translation; `codon_table` on `MT` is the 2 that
    /// `--codon-table MT=2` names.
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
        let no_flag = HashMap::new();
        let no_proteins = HashMap::new();
        let record = build_transcript(
            &gene(&coding),
            &coding,
            &sources(&no_proteins, Some(&fasta)),
            codon_table("21", None, &no_flag).unwrap(),
        )
        .unwrap();
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
        let mt_is_two = codon_table_flag(&["MT=2".to_string()]).unwrap();
        let record = build_transcript(
            &gene(&non_coding),
            &non_coding,
            &sources(&no_proteins, Some(&fasta)),
            codon_table("MT", None, &mt_is_two).unwrap(),
        )
        .unwrap();
        let json = serde_json::to_value(&record).unwrap();
        let vefc = &json["variation_effect_feature_cache"];
        assert_eq!(vefc["codon_table"], 2);
        assert_eq!(vefc["translateable_seq"], "");
        assert!(vefc.get("five_prime_utr").is_none() && vefc.get("three_prime_utr").is_none());
        assert!(json.get("translation").is_none());
        assert!(vefc.get("mapper").is_some());
    }

    /// Without a protein FASTA the peptide is the translation of
    /// `translateable_seq` under the region's codon table, as Ensembl VEP sets it
    /// from a GFF3 or GTF (`BaseGXF.pm:764-768`, `$translation->seq`): `MPKAL`
    /// for `ATGCCCAAAGCCCTTTAA`, the terminal `TAA` dropped (core
    /// `Transcript.pm:2092-2176`). A protein FASTA holding the transcript wins,
    /// as its residues carry the sequence edits; a coding transcript with no
    /// genome has no peptide.
    #[test]
    fn peptide_is_translated_when_no_protein_fasta_is_given() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = genome(dir.path());
        let coding = transcript(
            "21",
            1,
            vec![exon(1, 12, 1, 1), exon(21, 30, 1, 2)],
            vec![cds(4, 12, 1, 0), cds(21, 29, 1, 0)],
        );
        let peptide = |proteins: &HashMap<String, Peptide>, genome: Option<&_>| {
            let record =
                build_transcript(&gene(&coding), &coding, &sources(proteins, genome), 1).unwrap();
            let vefc = &record.variation_effect_feature_cache;
            (vefc.translateable_seq.clone(), vefc.peptide.clone())
        };
        assert_eq!(
            peptide(&HashMap::new(), Some(&fasta)),
            (
                Some("ATGCCCAAAGCCCTTTAA".to_string()),
                Some("MPKAL".to_string())
            )
        );
        let from_fasta = peptides("ENST1", "MPKALU");
        assert_eq!(
            peptide(&from_fasta, Some(&fasta)).1,
            Some("MPKALU".to_string())
        );
        assert_eq!(peptide(&HashMap::new(), None), (None, None));
    }

    /// A transcript takes the protein FASTA record held under its identifier
    /// and version (`ENST1.1` for `ENST1` at version 1), else under its
    /// identifier alone: the peptide and, under the same key, the record's
    /// prediction matrices. The record of another version (`ENST1.2`) gives
    /// it neither.
    #[test]
    fn a_record_under_the_transcripts_version_gives_its_peptide_and_predictions() {
        use crate::predictions::{PeptidePredictions, PredictionBlob};
        let coding = transcript(
            "21",
            1,
            vec![exon(1, 12, 1, 1), exon(21, 30, 1, 2)],
            vec![cds(4, 12, 1, 0), cds(21, 29, 1, 0)],
        );
        let built = |key: &str| {
            let proteins = peptides(key, "MPKALU");
            let predictions = HashMap::from([(
                key.to_string(),
                PeptidePredictions {
                    sift: Some(PredictionBlob {
                        analysis: "sift".to_string(),
                        sub_analysis: None,
                        translation_md5: "md5".to_string(),
                        base64_data: "matrix".to_string(),
                    }),
                    ..PeptidePredictions::default()
                },
            )]);
            let with_predictions = BuildSources {
                prediction_data: Some(&predictions),
                ..sources(&proteins, None)
            };
            let record = build_transcript(&gene(&coding), &coding, &with_predictions, 1).unwrap();
            let vefc = record.variation_effect_feature_cache;
            (vefc.peptide, vefc.protein_function_predictions.is_some())
        };
        assert_eq!(built("ENST1.1"), (Some("MPKALU".to_string()), true));
        assert_eq!(built("ENST1"), (Some("MPKALU".to_string()), true));
        assert_eq!(built("ENST1.2"), (None, false));
    }

    /// A record fits a coding sequence of as many whole codons as it has
    /// residues, or of one more, the stop codon; bases short of a codon are
    /// not counted. The record of 385 residues fits a sequence of 1,158
    /// bases and not one of 1,321.
    #[test]
    fn a_record_fits_its_whole_codons_or_one_fewer() {
        assert!(fits_coding_sequence(5, 18));
        assert!(fits_coding_sequence(6, 18));
        assert!(fits_coding_sequence(6, 20));
        assert!(!fits_coding_sequence(4, 18));
        assert!(!fits_coding_sequence(7, 18));
        assert!(fits_coding_sequence(385, 1158));
        assert!(!fits_coding_sequence(385, 1321));
    }

    /// The length rule holds a record to a coding sequence: a transcript
    /// with none keeps the record that names it, and a coding transcript
    /// whose record fits takes it and whose record does not is built as one
    /// with no record.
    #[test]
    fn the_length_rule_reads_a_record_against_a_coding_sequence() {
        let used = |cds: Vec<CdsRecord>, residues: &str| {
            let tr = transcript("21", 1, vec![exon(1, 12, 1, 1), exon(21, 30, 1, 2)], cds);
            let proteins = peptides("ENST1", residues);
            let record = build_transcript(&gene(&tr), &tr, &sources(&proteins, None), 1).unwrap();
            (
                record.record_use,
                record.variation_effect_feature_cache.peptide,
            )
        };
        let six_codons = || vec![cds(4, 12, 1, 0), cds(21, 29, 1, 0)];
        assert_eq!(
            used(Vec::new(), "MPKALU"),
            (RecordUse::Taken, Some("MPKALU".to_string()))
        );
        assert_eq!(
            used(six_codons(), "MPKAL"),
            (RecordUse::Taken, Some("MPKAL".to_string()))
        );
        assert_eq!(used(six_codons(), "MPK"), (RecordUse::SetAside, None));
    }

    /// A record set aside for its length gives the transcript neither its
    /// peptide nor the prediction matrices held under its key; a record that
    /// fits, held under the same key, gives both.
    #[test]
    fn a_record_set_aside_gives_no_predictions() {
        use crate::predictions::{PeptidePredictions, PredictionBlob};
        let coding = transcript(
            "21",
            1,
            vec![exon(1, 12, 1, 1), exon(21, 30, 1, 2)],
            vec![cds(4, 12, 1, 0), cds(21, 29, 1, 0)],
        );
        let built = |residues: &str| {
            let proteins = peptides("ENST1", residues);
            let predictions = HashMap::from([(
                "ENST1".to_string(),
                PeptidePredictions {
                    sift: Some(PredictionBlob {
                        analysis: "sift".to_string(),
                        sub_analysis: None,
                        translation_md5: "md5".to_string(),
                        base64_data: "matrix".to_string(),
                    }),
                    ..PeptidePredictions::default()
                },
            )]);
            let with_predictions = BuildSources {
                prediction_data: Some(&predictions),
                ..sources(&proteins, None)
            };
            let record = build_transcript(&gene(&coding), &coding, &with_predictions, 1).unwrap();
            let vefc = record.variation_effect_feature_cache;
            (
                record.record_use,
                vefc.peptide,
                vefc.protein_function_predictions.is_some(),
            )
        };
        assert_eq!(built("MPK"), (RecordUse::SetAside, None, false));
        assert_eq!(
            built("MPKAL"),
            (RecordUse::Taken, Some("MPKAL".to_string()), true)
        );
    }

    /// `Transcript::translate` cuts the sequence to whole codons, drops the
    /// terminal stop codon's residue and writes `M` for a first codon that is a
    /// start codon of the table (`CTG` on table 1, `ATT` on table 2); table 2
    /// reads `TGA` as `W` and `ATA` as `M`; an ambiguous codon resolves when
    /// every expansion agrees (`GCN` is `A`), to `B` for Asp or Asn, `Z` for
    /// Glu or Gln, else `X`, and is a start or stop codon when every expansion
    /// is one (`YTG` and `TAR` under table 1; `TRA` is not under table 2, where
    /// `TGA` is Trp) (`Bio::Tools::CodonTable`).
    #[test]
    fn translation_follows_transcript_translate() {
        assert_eq!(translate("ATGCCCAAAGCCCTT", 1).as_deref(), Some("MPKAL"));
        assert_eq!(translate("ATGCCCAAAGCCCTTTAA", 1).as_deref(), Some("MPKAL"));
        assert_eq!(translate("ATGCC", 1).as_deref(), Some("M"));
        assert_eq!(translate("AT", 1), None);
        assert_eq!(translate("CTGAAA", 1).as_deref(), Some("MK"));
        assert_eq!(translate("TTGAAA", 1).as_deref(), Some("MK"));
        assert_eq!(translate("GTGAAA", 1).as_deref(), Some("VK"));
        assert_eq!(translate("NTGAAA", 1).as_deref(), Some("XK"));
        assert_eq!(translate("NATGCNRATSAATGA", 1).as_deref(), Some("XABZ"));
        assert_eq!(translate("TGAAAACCC", 1).as_deref(), Some("*KP"));
        assert_eq!(translate("TGAAAACCC", 2).as_deref(), Some("WKP"));
        assert_eq!(translate("ATTAAAAGA", 2).as_deref(), Some("MK"));
        assert_eq!(translate("ATATAA", 2).as_deref(), Some("M"));
        assert_eq!(translate("GTGAAA", 2).as_deref(), Some("MK"));
        assert_eq!(translate("TGA", 1).as_deref(), Some(""));
        assert_eq!(translate("atgcccaaa", 1).as_deref(), Some("MPK"));
        assert_eq!(translate("YTGAAATAR", 1).as_deref(), Some("MK"));
        assert_eq!(translate("NTGAAATRA", 2).as_deref(), Some("XKX"));
    }

    /// The translation's identifier and version are its first CDS line's in
    /// transcript order (`BaseGXF.pm` `_add_translation`, the first of the
    /// ordered CDS records): the segment of lowest start on the forward
    /// strand, of highest start on the reverse, so CDS lines that disagree
    /// yield the 5'-most line's pair on either strand.
    #[test]
    fn translation_identifiers_come_from_the_five_prime_cds_line() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = genome(dir.path());
        let versioned = |start, end, strand, version| CdsRecord {
            protein_id: Some(format!("ENSP{version}")),
            protein_version: Some(version),
            ..cds(start, end, strand, 0)
        };
        let pair_on = |strand: i8| {
            let (first, second) = if strand > 0 { (1, 2) } else { (2, 1) };
            let tr = transcript(
                "21",
                strand,
                vec![exon(1, 12, strand, first), exon(21, 30, strand, second)],
                vec![versioned(4, 12, strand, 1), versioned(21, 29, strand, 2)],
            );
            let record =
                build_transcript(&gene(&tr), &tr, &sources(&HashMap::new(), Some(&fasta)), 1)
                    .unwrap();
            let translation = record.translation.unwrap();
            (translation.stable_id, translation.version)
        };
        assert_eq!(pair_on(1), ("ENSP1".to_string(), Some(1)));
        assert_eq!(pair_on(-1), ("ENSP2".to_string(), Some(2)));
    }

    /// One case per arm of the version's precedence: the CDS line's `version`
    /// over the dumps' and the header's; the dumps' translation version,
    /// present or absent, over the header's; the header's when no dumps hold
    /// the translation; none when no source carries one. On a build from a
    /// protein FASTA, the header's version reaches the written translation
    /// when the CDS line carries none and stays out when it does.
    #[test]
    fn translation_version_comes_from_the_first_source_that_carries_one() {
        let facts = |version| TranslationFacts {
            stable_id: "ENSP1".to_string(),
            version,
            ..TranslationFacts::default()
        };
        assert_eq!(
            translation_version(Some(2), Some(&facts(Some(3))), Some(4)),
            Some(2)
        );
        assert_eq!(
            translation_version(None, Some(&facts(Some(3))), Some(4)),
            Some(3)
        );
        assert_eq!(translation_version(None, Some(&facts(None)), Some(4)), None);
        assert_eq!(translation_version(None, None, Some(4)), Some(4));
        assert_eq!(translation_version(None, None, None), None);

        let dir = tempfile::tempdir().unwrap();
        let fasta = genome(dir.path());
        let from_fasta = HashMap::from([(
            "ENST1".to_string(),
            Peptide {
                sequence: "MPKAL".to_string(),
                identifier: "ENSP1.4".to_string(),
                version: Some(4),
            },
        )]);
        let version_with_cds = |cds_version| {
            let tr = transcript(
                "21",
                1,
                vec![exon(1, 12, 1, 1), exon(21, 30, 1, 2)],
                vec![CdsRecord {
                    protein_id: Some("ENSP1".to_string()),
                    protein_version: cds_version,
                    ..cds(4, 29, 1, 0)
                }],
            );
            build_transcript(&gene(&tr), &tr, &sources(&from_fasta, Some(&fasta)), 1)
                .unwrap()
                .translation
                .unwrap()
                .version
        };
        assert_eq!(version_with_cds(None), Some(4));
        assert_eq!(version_with_cds(Some(1)), Some(1));
    }

    /// The release 116 GRCh38 core dump slice of the integration tests, whose
    /// `seq_region_attrib` gives `MT` codon table 2 and `21` none.
    fn chr21_dumps() -> CoreDumps {
        CoreDumps::load(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21"))
            .unwrap()
    }

    /// The codon table is the sequence region's `codon_table` attribute when the
    /// dumps carry one (`Database/Transcript.pm:359-364`: `MT` reads 2 with no
    /// flag), else the `--codon-table` value (`mitochondrion_genome=5` with no
    /// dumps reads 5; the flag's `chrM` is the readers' `MT`), else 1 (`21` with
    /// neither). The region's name never decides: the dumps' 2 beside a flag's
    /// `MT=3` is refused, as is a flag whose table BioPerl lacks or whose
    /// value names no region.
    #[test]
    fn codon_table_is_the_regions_attribute_or_the_flag_never_the_name() {
        let dumps = chr21_dumps();
        let no_flag = HashMap::new();
        assert_eq!(codon_table("MT", Some(&dumps), &no_flag).unwrap(), 2);
        assert_eq!(codon_table("21", Some(&dumps), &no_flag).unwrap(), 1);
        assert_eq!(codon_table("21", None, &no_flag).unwrap(), 1);

        let fly = codon_table_flag(&["mitochondrion_genome=5".to_string()]).unwrap();
        assert_eq!(codon_table("mitochondrion_genome", None, &fly).unwrap(), 5);
        assert_eq!(codon_table("2L", None, &fly).unwrap(), 1);

        let joined = codon_table_flag(&["Pt=11".to_string(), "chrM=2".to_string()]).unwrap();
        assert_eq!(
            joined,
            HashMap::from([("Pt".to_string(), 11), ("MT".to_string(), 2)])
        );
        assert_eq!(codon_table("MT", Some(&dumps), &joined).unwrap(), 2);

        let mt_is_three = codon_table_flag(&["MT=3".to_string()]).unwrap();
        assert_eq!(
            codon_table("MT", Some(&dumps), &mt_is_three)
                .unwrap_err()
                .to_string(),
            "codon table of sequence region MT is 2 in the core dumps and 3 in --codon-table"
        );

        for refused in ["MT=7", "MT=25", "MT=x", "MT", "=2", "MT=2=3"] {
            assert!(
                codon_table_flag(&[refused.to_string()]).is_err(),
                "{refused} parses"
            );
        }
        assert!(codon_table_flag(&["MT=2".to_string(), "MT=3".to_string()]).is_err());
        assert_eq!(
            codon_table_flag(&["MT=2".to_string(), "MT=2".to_string()]).unwrap(),
            HashMap::from([("MT".to_string(), 2)])
        );
    }

    /// A region named like an organelle genome, in any letter case, whose table
    /// neither the dumps nor the flag states refuses the build with the line
    /// naming both remedies; with dumps that carry no attribute for it, the
    /// line names the flag alone. A nuclear name reads 1, and `codon_tables`
    /// over a model whose one region is `Mito` refuses with the same line.
    #[test]
    fn an_organelle_named_region_without_a_table_refuses() {
        let no_flag = HashMap::new();
        for name in ORGANELLE_REGION_NAMES {
            for spelled in [name.to_string(), name.to_uppercase(), name.to_lowercase()] {
                assert_eq!(
                    codon_table(&spelled, None, &no_flag).unwrap_err().to_string(),
                    format!(
                        "codon table of sequence region {spelled} unknown: pass --core-dumps or --codon-table {spelled}=<n>"
                    )
                );
            }
        }
        let dumps = chr21_dumps();
        assert_eq!(
            codon_table("Mito", Some(&dumps), &no_flag).unwrap_err().to_string(),
            "codon table of sequence region Mito unknown: the core dumps carry no codon_table attribute for it; pass --codon-table Mito=<n>"
        );
        assert_eq!(codon_table("X", None, &no_flag).unwrap(), 1);
        assert_eq!(codon_table("KI270728.1", None, &no_flag).unwrap(), 1);

        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/mito_codon_table");
        let genes =
            crate::gff3::parse_gff3(&fixture.join("mito.gff3"), crate::gff3::Reading::GffSource)
                .unwrap();
        assert_eq!(
            codon_tables(&genes, None, &no_flag).unwrap_err().to_string(),
            "codon table of sequence region Mito unknown: pass --core-dumps or --codon-table Mito=<n>"
        );
    }

    /// A synthetic `Mito` GFF3 built under `--codon-table Mito=3` writes the
    /// record with `codon_table` 3 and the CDS line's bases as its
    /// `translateable_seq`. The peptide is the protein FASTA's `MTWMK`, and
    /// without that FASTA the translation of the CDS under table 3 reads the
    /// same (`CTT` Thr, `TGA` Trp and `ATA` Met, where table 1 reads `ML*IK`).
    #[test]
    fn a_build_from_a_synthetic_gff3_under_table_three_writes_three() {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/mito_codon_table");
        let genes =
            crate::gff3::parse_gff3(&fixture.join("mito.gff3"), crate::gff3::Reading::GffSource)
                .unwrap();
        let genome = vep_fasta::IndexedFasta::from_path(fixture.join("mito.fa")).unwrap();
        let peptides =
            crate::fasta::parse_protein_fasta(&fixture.join("mito_peptides.fa")).unwrap();
        let named = codon_table_flag(&["Mito=3".to_string()]).unwrap();
        let tables = codon_tables(&genes, None, &named).unwrap();
        assert_eq!(tables, HashMap::from([("Mito".to_string(), 3)]));

        let cds = "ATGCTTTGAATAAAATAA";
        assert_eq!(translate(cds, 3).as_deref(), Some("MTWMK"));
        assert_eq!(translate(cds, 1).as_deref(), Some("ML*IK"));

        let build = |proteins: &HashMap<String, Peptide>| {
            let out = tempfile::tempdir().unwrap();
            let stats = build_cache(
                &genes,
                &fixture.join("mito.gff3"),
                out.path(),
                &sources(proteins, Some(&genome)),
                &tables,
            )
            .unwrap();
            assert_eq!(
                (stats.transcript_count, stats.chromosomes),
                (1, vec!["Mito".to_string()])
            );
            let shard = std::fs::read_to_string(out.path().join("transcripts/Mito/1-1000000.json"))
                .unwrap();
            let records: serde_json::Value = serde_json::from_str(&shard).unwrap();
            records[0]["variation_effect_feature_cache"].clone()
        };
        let vefc = build(&peptides);
        assert_eq!(vefc["codon_table"], 3);
        assert_eq!(vefc["translateable_seq"], cds);
        assert_eq!(vefc["peptide"], "MTWMK");
        assert_eq!(vefc["five_prime_utr"], "GGGGGCCCCC");
        assert_eq!(vefc["three_prime_utr"], "TT");
        let vefc = build(&HashMap::new());
        assert_eq!(vefc["codon_table"], 3);
        assert_eq!(vefc["peptide"], "MTWMK");
    }

    /// A build from the model and the genome alone writes every transcript
    /// of the tree, an `artifact` one and a readthrough one included, with
    /// the biotype the line carries, as Ensembl VEP's `--gff` and `--gtf`
    /// sources annotate them (`BaseGXF.pm` `_create_transcript`). A build
    /// reading a protein FASTA, or the core dumps alone, writes neither, nor a
    /// transcript the dumps mark `readthrough_tra`, as Ensembl VEP's cache
    /// carries none of them (`Database/Transcript.pm`
    /// `get_features_by_regions_uncached`).
    #[test]
    fn artifact_and_readthrough_transcripts_stay_on_a_genome_only_build() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = genome(dir.path());
        let named = |id: &str| {
            let mut tr = transcript("21", 1, vec![exon(1, 12, 1, 1), exon(21, 30, 1, 2)], vec![]);
            tr.transcript_id = id.to_string();
            tr
        };
        let mut artifact = named("ENST1");
        artifact.biotype = "artifact".to_string();
        let mut readthrough = named("ENST2");
        readthrough.readthrough = true;
        // Marked `readthrough_tra` by the fixture dumps, untagged on its line.
        let dump_marked = named("ENST00000612267");
        let plain = named("ENST3");
        let mut gene_record = gene(&plain);
        gene_record.transcripts = vec![artifact, readthrough, dump_marked, plain];
        let genes = HashMap::from([("ENSG1".to_string(), gene_record)]);
        let tables = HashMap::from([("21".to_string(), 1u8)]);

        let build = |sources: BuildSources| {
            let out = tempfile::tempdir().unwrap();
            let stats = build_cache(
                &genes,
                Path::new("model.gff3"),
                out.path(),
                &sources,
                &tables,
            )
            .unwrap();
            let shard =
                std::fs::read_to_string(out.path().join("transcripts/21/1-1000000.json")).unwrap();
            let records: Vec<serde_json::Value> = serde_json::from_str(&shard).unwrap();
            let written: Vec<(String, String)> = records
                .iter()
                .map(|r| {
                    (
                        r["stable_id"].as_str().unwrap().to_string(),
                        r["biotype"].as_str().unwrap().to_string(),
                    )
                })
                .collect();
            (stats.transcript_count, written)
        };
        let pair = |id: &str, biotype: &str| (id.to_string(), biotype.to_string());
        let no_fasta = HashMap::new();
        assert_eq!(
            build(sources(&no_fasta, Some(&fasta))),
            (
                4,
                vec![
                    pair("ENST00000612267", "lncRNA"),
                    pair("ENST1", "artifact"),
                    pair("ENST2", "lncRNA"),
                    pair("ENST3", "lncRNA"),
                ]
            )
        );
        let from_fasta = peptides("ENST3", "MPKAL");
        assert_eq!(
            build(sources(&from_fasta, Some(&fasta))),
            (
                2,
                vec![pair("ENST00000612267", "lncRNA"), pair("ENST3", "lncRNA")]
            )
        );
        let dumps = chr21_dumps();
        assert!(dumps.transcript("ENST00000612267").unwrap().readthrough);
        let dumps_alone = BuildSources {
            core_dumps: Some(&dumps),
            ..sources(&no_fasta, Some(&fasta))
        };
        assert_eq!(build(dumps_alone), (1, vec![pair("ENST3", "lncRNA")]));
    }

    /// Chromosome names that read as one number tie in `contig_order`, which
    /// keeps tied names in the order it is given them: a build lists them in
    /// name order, the same on every run.
    #[test]
    fn chromosome_names_reading_as_one_number_are_listed_in_name_order() {
        let names = ["1", "001", "01"];
        let genes: HashMap<String, GeneRecord> = names
            .iter()
            .map(|name| {
                let tr = transcript(name, 1, vec![exon(1, 12, 1, 1)], vec![]);
                let mut gene_record = gene(&tr);
                gene_record.transcripts = vec![tr];
                (format!("ENSG{name}"), gene_record)
            })
            .collect();
        let tables = names.iter().map(|name| (name.to_string(), 1u8)).collect();
        let out = tempfile::tempdir().unwrap();
        let stats = build_cache(
            &genes,
            Path::new("model.gff3"),
            out.path(),
            &sources(&HashMap::new(), None),
            &tables,
        )
        .unwrap();
        assert_eq!(stats.chromosomes, ["001", "01", "1"]);
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
                build_transcript(&gene(&tr), &tr, &sources(&HashMap::new(), None), 1).unwrap();
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
