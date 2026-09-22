// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Transcript and related types for cached annotation data.
//!
//! Defines the [`Transcript`] type and its components (exons, introns, translation,
//! coordinate mapper, VEFC cache) as loaded from the VEP JSON or binary cache.
//! These types are the primary input to the consequence calculator in `vep-effects`.

use std::sync::{Arc, OnceLock};

use crate::coordinate::Strand;

/// An Ensembl transcript loaded from the VEP cache.
///
/// Represents a single transcript isoform with its full annotation: genomic
/// coordinates, exon/intron structure, coding region boundaries, translation
/// product, and the pre-computed VEFC cache (translateable sequence, coordinate
/// mapper, protein features). This is the central data structure for consequence
/// calculation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Transcript {
    /// Ensembl stable ID (e.g., "ENST00000366667").
    pub stable_id: Arc<str>,
    /// Stable ID version.
    pub version: Option<u32>,
    /// Database ID.
    pub db_id: Option<u64>,
    /// Gene stable ID (e.g., "ENSG00000142168").
    pub gene_stable_id: Arc<str>,
    /// Chromosome.
    pub chr: String,
    /// 1-based start position.
    pub start: u64,
    /// 1-based end position (inclusive).
    pub end: u64,
    /// Strand.
    pub strand: Strand,
    /// Biotype (e.g., "protein_coding", "lncRNA").
    pub biotype: Arc<str>,
    /// Source (Ensembl or RefSeq).
    pub source: String,
    /// Human-readable description.
    pub description: Option<String>,
    /// Gene symbol (e.g., "SOD1").
    pub gene_symbol: Option<Arc<str>>,
    /// Gene symbol source (e.g., "HGNC").
    pub gene_symbol_source: Option<Arc<str>>,
    /// HGNC ID.
    pub hgnc_id: Option<Arc<str>>,
    /// Gene phenotype association count (number of phenotypes associated with
    /// this gene; values can exceed 255 for highly-studied genes like TP53,
    /// hence u32 not u8).
    pub gene_phenotype: Option<u32>,
    /// Whether this is the canonical transcript for the gene.
    pub canonical: bool,
    /// MANE Select transcript ID.
    pub mane_select: Option<String>,
    /// MANE Plus Clinical transcript ID.
    pub mane_plus_clinical: Option<String>,
    /// Transcript Support Level (1-5, NA).
    pub tsl: Option<u8>,
    /// APPRIS annotation.
    pub appris: Option<String>,
    /// CCDS ID.
    pub ccds: Option<String>,
    /// Protein (ENSP) ID.
    pub protein_id: Option<Arc<str>>,
    /// RefSeq transcript ID.
    pub refseq: Option<String>,
    /// UniProtKB/Swiss-Prot accession(s).
    pub swissprot: Option<String>,
    /// UniProtKB/TrEMBL accession(s).
    pub trembl: Option<String>,
    /// UniParc accession(s).
    pub uniparc: Option<String>,
    /// Exons.
    pub exons: Vec<Exon>,
    /// Introns from the cache's intron list, else the gaps between exons, filled at load time.
    pub introns: Vec<Intron>,
    /// cDNA coding start position (1-based offset within cDNA).
    pub cdna_coding_start: Option<u64>,
    /// cDNA coding end position (1-based offset within cDNA).
    pub cdna_coding_end: Option<u64>,
    /// Coding region start (genomic position).
    pub coding_region_start: Option<u64>,
    /// Coding region end (genomic position).
    pub coding_region_end: Option<u64>,
    /// Translation start (genomic position).
    pub translation_start: Option<u64>,
    /// Translation end (genomic position).
    pub translation_end: Option<u64>,
    /// Translation object (protein product).
    pub translation: Option<Translation>,
    /// cDNA sequence (populated from cache).
    pub cdna_sequence: Option<Vec<u8>>,
    /// Protein sequence (populated from cache).
    pub protein_sequence: Option<Vec<u8>>,
    /// Transcript flags (e.g., "cds_start_NF", "cds_end_NF").
    pub flags: Arc<[String]>,
    /// GENCODE primary flag.
    pub gencode_primary: bool,
    /// Transcript attributes (code/value pairs from the Perl cache).
    pub attributes: Vec<Attribute>,
    /// Cached data for variant effect calculation (_variation_effect_feature_cache).
    pub vefc: Option<TranscriptVEFC>,
    /// Facts derived from the fields above, filled on first use and read by
    /// every later row against this transcript ([`Transcript::facts`],
    /// [`Transcript::feature_ordinals`]).
    #[serde(skip)]
    pub derived: TranscriptDerived,
}

/// The lazily filled cells behind [`Transcript::facts`] and
/// [`Transcript::feature_ordinals`]. A fresh or cloned transcript starts empty,
/// so a transcript edited after construction is read as edited.
#[derive(Debug, Default)]
pub struct TranscriptDerived {
    facts: OnceLock<TranscriptFacts>,
    ordinals: OnceLock<Box<FeatureOrdinals>>,
}

impl Clone for TranscriptDerived {
    fn clone(&self) -> Self {
        Self::default()
    }
}

/// Per-transcript facts that every row against the transcript reads; the
/// counterparts of what Perl VEP stores on the transcript's
/// `_variation_effect_feature_cache` once and reads per allele.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TranscriptFacts {
    /// An intron of the transcript list (or, with no intron list, a gap between
    /// consecutive exons) spans at most 12 bp: Perl's `_has_frameshift_intron`,
    /// which stretches every exon by 12 bases in the exon-overlap test
    /// (`BaseTranscriptVariation::_overlapped_exons`).
    pub has_frameshift_intron: bool,
    /// The same test read from the VEFC intron list, falling back to the
    /// transcript list. Equal to `has_frameshift_intron` for every cache-loaded
    /// transcript, whose two lists are one; a hand-built transcript may set them
    /// apart, and the consequence path reads this one.
    pub vefc_has_frameshift_intron: bool,
    /// The `cds_start_NF` flag is present.
    pub cds_start_nf: bool,
    /// The `cds_end_NF` flag is present.
    pub cds_end_nf: bool,
    /// The transcript has a translation and a mapper with a coding range.
    pub has_coding_model: bool,
    /// Phase of the 5' exon when positive, else 0: the number of padding bases
    /// the cached translateable sequence starts with.
    pub start_phase_offset: i8,
}

impl TranscriptFacts {
    fn compute(tx: &Transcript) -> Self {
        let has_frameshift_intron = if tx.introns.is_empty() {
            let mut exons: Vec<(i64, i64)> = tx
                .exons
                .iter()
                .map(|e| (e.start as i64, e.end as i64))
                .collect();
            exons.sort_unstable();
            exons
                .windows(2)
                .any(|w| ((w[1].0 - 1) - (w[0].1 + 1)).abs() <= 12)
        } else {
            tx.introns
                .iter()
                .any(|i| (i.end as i64 - i.start as i64).abs() <= 12)
        };
        let vefc_introns = tx
            .vefc
            .as_ref()
            .filter(|v| !v.introns.is_empty())
            .map(|v| v.introns.as_slice())
            .unwrap_or(tx.introns.as_slice());
        let vefc_has_frameshift_intron = vefc_introns
            .iter()
            .any(|i| i.end.saturating_sub(i.start) <= 12);
        let has_coding_model = tx.translation.is_some()
            && tx
                .vefc
                .as_ref()
                .and_then(|v| v.mapper.as_ref())
                .is_some_and(|m| {
                    m.cdna_coding_start > 0 && m.cdna_coding_end >= m.cdna_coding_start
                });
        let first_exon = match tx.strand {
            Strand::Forward => tx.exons.iter().min_by_key(|e| e.start),
            Strand::Reverse => tx.exons.iter().max_by_key(|e| e.end),
        };
        let start_phase_offset = match first_exon {
            Some(e) if e.phase > 0 => e.phase,
            _ => 0,
        };
        TranscriptFacts {
            has_frameshift_intron,
            vefc_has_frameshift_intron,
            cds_start_nf: tx.flags.iter().any(|f| f == "cds_start_NF"),
            cds_end_nf: tx.flags.iter().any(|f| f == "cds_end_NF"),
            has_coding_model,
            start_phase_offset,
        }
    }
}

/// The transcript's exons and introns as `(start, end)` in transcript order
/// (5' to 3'), the order VEP numbers them in for the `EXON` and `INTRON`
/// columns.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureOrdinals {
    pub exons: Box<[(i64, i64)]>,
    /// From the transcript's intron list when present, else the gaps between
    /// consecutive exons.
    pub introns: Box<[(i64, i64)]>,
}

impl FeatureOrdinals {
    fn compute(tx: &Transcript) -> Self {
        let reverse = tx.strand == Strand::Reverse;
        let mut exons: Vec<(i64, i64)> = tx
            .exons
            .iter()
            .map(|e| (e.start as i64, e.end as i64))
            .collect();
        exons.sort_unstable();
        if reverse {
            exons.reverse();
        }
        let introns: Vec<(i64, i64)> = if tx.introns.is_empty() {
            exons
                .windows(2)
                .map(|w| {
                    let (a, b) = (w[0], w[1]);
                    if a.0 < b.0 {
                        (a.1 + 1, b.0 - 1)
                    } else {
                        (b.1 + 1, a.0 - 1)
                    }
                })
                .collect()
        } else {
            let mut introns: Vec<(i64, i64)> = tx
                .introns
                .iter()
                .map(|i| (i.start as i64, i.end as i64))
                .collect();
            introns.sort_unstable();
            if reverse {
                introns.reverse();
            }
            introns
        };
        FeatureOrdinals {
            exons: exons.into_boxed_slice(),
            introns: introns.into_boxed_slice(),
        }
    }
}

/// An exon within a transcript.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Exon {
    /// Exon stable ID (e.g., "ENSE00001472862").
    pub stable_id: Option<String>,
    /// 1-based start position.
    pub start: u64,
    /// 1-based end position (inclusive).
    pub end: u64,
    /// Exon rank (1-based).
    pub rank: u32,
    /// Phase at the start of the exon (0, 1, 2, or -1 for non-coding).
    pub phase: i8,
    /// Phase at the end of the exon.
    pub end_phase: i8,
}

/// The protein product of a transcript.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Translation {
    /// Ensembl protein stable ID (e.g., "ENSP00000355627").
    pub stable_id: String,
    /// Stable ID version.
    pub version: Option<u32>,
    /// Database ID.
    pub db_id: Option<u64>,
    /// Translation start position within the start exon (1-based).
    pub start: u64,
    /// Translation end position within the end exon (1-based).
    pub end: u64,
    /// Index into the parent transcript's exon array for the start exon.
    pub start_exon_index: usize,
    /// Index into the parent transcript's exon array for the end exon.
    pub end_exon_index: usize,
    /// Protein sequence from the Perl cache.
    pub seq: Option<String>,
}

/// A code/value attribute pair from the Perl cache.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Attribute {
    pub code: String,
    pub value: String,
}

/// Pre-computed cache for variant effect calculation (`_variation_effect_feature_cache`).
///
/// Contains the CDS (translateable) sequence, UTR sequences, coordinate mapper
/// pairs, and protein features needed for consequence calculation. Populated
/// during cache loading from the Perl Storable or JSON cache.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct TranscriptVEFC {
    /// Codon table number (usually 1 for the standard genetic code).
    pub codon_table: u8,
    /// 5' UTR sequence.
    pub five_prime_utr: Option<String>,
    /// 3' UTR sequence.
    pub three_prime_utr: Option<String>,
    /// CDS (translateable) sequence.
    pub translateable_seq: Option<String>,
    /// Protein (peptide) sequence.
    pub peptide: Option<String>,
    /// Introns from the VEFC (may duplicate transcript-level introns).
    pub introns: Vec<Intron>,
    /// Exons sorted by rank for the VEFC.
    pub sorted_exons: Vec<Exon>,
    /// Coordinate mapper between cDNA and genomic coordinates.
    pub mapper: Option<TranscriptMapper>,
    /// Protein domain features.
    pub protein_features: Vec<ProteinFeature>,
    /// SIFT/PolyPhen prediction matrices.
    pub protein_function_predictions: Option<ProteinFunctionPredictions>,
    /// Sequence edits (selenocysteine insertions, RNA edits, etc.).
    pub seq_edits: Vec<SeqEdit>,
}

/// Coordinate mapper between cDNA and genomic coordinates.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TranscriptMapper {
    /// Phase at the start of the CDS.
    pub start_phase: i8,
    /// cDNA position of the coding start (1-based).
    pub cdna_coding_start: u64,
    /// cDNA position of the coding end (1-based).
    pub cdna_coding_end: u64,
    /// Exon-level coordinate mapping pairs.
    pub exon_coord_mapper: ExonCoordMapper,
}

/// Collection of coordinate mapping pairs for exon-level cDNA <-> genomic mapping.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ExonCoordMapper {
    /// Number of mapping pairs.
    pub pair_count: usize,
    /// The individual mapping pairs, stored in cDNA (transcript 5'->3') order:
    /// `pairs[i].from_start` ascends with `i`. On the reverse strand the genomic
    /// `to_start` therefore *descends* with `i`, so genomic lookups cannot binary
    /// search `pairs` directly; use `by_genomic_start` below.
    pub pairs: Vec<MapperPair>,
    /// Permutation of `pairs` indices sorted ascending by `pairs[i].to_start`
    /// (genomic). Enables O(log n) genomic->cDNA lookup regardless of strand.
    /// Derived at load time, never serialized; an empty/length-mismatched vec
    /// signals callers to fall back to the linear scan.
    #[serde(skip)]
    pub by_genomic_start: Vec<u32>,
}

impl ExonCoordMapper {
    /// Build an `ExonCoordMapper` from its pairs, computing the genomic-sorted
    /// index. This is the single construction funnel so the index can never be
    /// silently forgotten by a caller.
    #[must_use]
    pub fn new(pairs: Vec<MapperPair>) -> Self {
        let by_genomic_start = Self::build_genomic_index(&pairs);
        Self {
            pair_count: pairs.len(),
            pairs,
            by_genomic_start,
        }
    }

    /// Sort `0..pairs.len()` ascending by `pairs[i].to_start` (genomic start).
    /// Exon genomic intervals are disjoint, so the resulting order is a total
    /// order with no ties that matter for containment lookups.
    #[must_use]
    pub fn build_genomic_index(pairs: &[MapperPair]) -> Vec<u32> {
        let mut idx: Vec<u32> = (0..pairs.len() as u32).collect();
        idx.sort_unstable_by_key(|&i| pairs[i as usize].to_start);
        idx
    }

    /// Rebuild the genomic index from the current `pairs`. Call after any path
    /// that constructs/deserializes `pairs` without going through [`new`].
    pub fn rebuild_genomic_index(&mut self) {
        self.by_genomic_start = Self::build_genomic_index(&self.pairs);
    }
}

/// A single coordinate mapping pair (cDNA <-> genomic).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MapperPair {
    /// cDNA start (1-based).
    pub from_start: u64,
    /// cDNA end (1-based).
    pub from_end: u64,
    /// Genomic start (1-based).
    pub to_start: u64,
    /// Genomic end (1-based).
    pub to_end: u64,
    /// Orientation: 1 or -1.
    pub ori: i8,
}

/// A protein domain feature overlapping the translation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ProteinFeature {
    /// Start position in the protein (1-based).
    pub start: u64,
    /// End position in the protein (1-based).
    pub end: u64,
    /// Domain ID (e.g., "cd01667", "PF00080").
    pub hseqname: String,
    /// Analysis name (e.g., "Pfam", "PROSITE_profiles").
    pub analysis: Option<String>,
}

/// Pre-computed SIFT/PolyPhen prediction matrices for a transcript.
///
/// These matrices are sourced from Ensembl's variation database (public MySQL at
/// `ensembldb.ensembl.org`, keyed by MD5 of the protein sequence) and embedded in
/// the JSON cache by `vep-cache-builder --include-predictions`.
///
/// At runtime, predictions are decoded by [`crate::prediction::decode_prediction`].
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ProteinFunctionPredictions {
    pub sift: Option<PredictionMatrix>,
    pub polyphen_humdiv: Option<PredictionMatrix>,
    pub polyphen_humvar: Option<PredictionMatrix>,
}

/// A SIFT or PolyPhen prediction matrix for one transcript.
///
/// Stores a binary matrix encoding predictions for every possible amino acid
/// substitution at every protein position. The format is defined by
/// `ProteinFunctionPredictionMatrix.pm` in ensembl-variation release/115:
///
/// - Gzip-compressed blob (decompressed at load time or lazily on first lookup)
/// - Header: 3 bytes `"VEP"`
/// - Per position: 20 amino acids × 2 bytes = 40 bytes
/// - Each cell: `u16` little-endian; top 2 bits (15-14) = category index, bits 13-10
///   unused, bottom 10 = score × 1000
/// - `0xFFFF` = no prediction
/// - Amino acid order: A C D E F G H I K L M N P Q R S T V W Y
///
/// See [`crate::prediction`] for the decoder implementation and format documentation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PredictionMatrix {
    /// Analysis name: "sift", "polyphen_humdiv", or "polyphen_humvar".
    pub analysis: String,
    /// Optional sub-analysis identifier.
    pub sub_analysis: Option<String>,
    /// Length of the protein (number of residues).
    pub peptide_length: usize,
    /// Decompressed prediction matrix bytes (header + position data).
    /// Use [`crate::prediction::decode_prediction`] to look up values.
    pub predictions_data: Vec<u8>,
}

impl PredictionMatrix {
    /// Look up the prediction for a given protein position and amino acid.
    ///
    /// Returns `None` if the position/AA is out of range or no prediction exists.
    pub fn lookup(
        &self,
        position: usize,
        amino_acid: u8,
        analysis: crate::prediction::AnalysisType,
    ) -> Option<crate::prediction::Prediction> {
        crate::prediction::decode_prediction(&self.predictions_data, position, amino_acid, analysis)
    }
}

/// A sequence edit (e.g., selenocysteine insertion, RNA editing).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SeqEdit {
    /// Start position (1-based, in cDNA or protein coordinates according to context).
    pub start: u64,
    /// End position (1-based).
    pub end: u64,
    /// Replacement sequence.
    pub alt_seq: String,
}

/// An intron within a transcript.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Intron {
    /// 1-based start position (first intronic base).
    pub start: u64,
    /// 1-based end position (last intronic base).
    pub end: u64,
    /// Intron rank (1-based).
    pub rank: u32,
}

impl Transcript {
    /// The per-transcript facts, computed on the first call and shared by every
    /// later one.
    pub fn facts(&self) -> &TranscriptFacts {
        self.derived
            .facts
            .get_or_init(|| TranscriptFacts::compute(self))
    }

    /// The exon and intron ordinal lists, computed on the first call and
    /// shared by every later one.
    pub fn feature_ordinals(&self) -> &FeatureOrdinals {
        self.derived
            .ordinals
            .get_or_init(|| Box::new(FeatureOrdinals::compute(self)))
    }

    /// Length of the transcript on the genome.
    pub fn genomic_length(&self) -> u64 {
        if self.end >= self.start {
            self.end - self.start + 1
        } else {
            0
        }
    }

    /// Total number of exons.
    pub fn exon_count(&self) -> usize {
        self.exons.len()
    }

    /// Total number of introns.
    pub fn intron_count(&self) -> usize {
        if self.exons.len() > 1 {
            self.exons.len() - 1
        } else {
            0
        }
    }

    /// Whether this transcript is protein-coding.
    pub fn is_protein_coding(&self) -> bool {
        &*self.biotype == "protein_coding"
    }

    /// Whether this transcript has a CDS (coding sequence).
    /// This includes protein_coding and NMD transcripts, which have CDS/UTR
    /// regions even though they are targeted for degradation. Perl VEP's SV
    /// annotation engine treats NMD transcripts as coding (they have CDS exons,
    /// 5'UTR, 3'UTR) and assigns coding_sequence_variant, stop_lost, etc.
    pub fn has_cds(&self) -> bool {
        self.coding_region_start.is_some()
    }

    /// Whether this transcript is annotated as an NMD transcript.
    pub fn is_nmd_transcript(&self) -> bool {
        self.biotype.eq_ignore_ascii_case("nonsense_mediated_decay")
    }

    /// Whether this transcript is non-coding for consequence-context purposes.
    ///
    /// NMD is a separate context term (`NMD_transcript_variant`), mirroring Perl
    /// VEP output, so NMD transcripts are excluded here.
    pub fn is_non_coding_context_transcript(&self) -> bool {
        !self.is_protein_coding() && !self.is_nmd_transcript()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mapper_pair_construction() {
        let pair = MapperPair {
            from_start: 1,
            from_end: 100,
            to_start: 25000001,
            to_end: 25000100,
            ori: 1,
        };
        assert_eq!(pair.from_start, 1);
        assert_eq!(pair.from_end, 100);
        assert_eq!(pair.to_start, 25000001);
        assert_eq!(pair.to_end, 25000100);
        assert_eq!(pair.ori, 1);
    }

    #[test]
    fn test_transcript_mapper_construction() {
        let mapper = TranscriptMapper {
            start_phase: 0,
            cdna_coding_start: 155,
            cdna_coding_end: 622,
            exon_coord_mapper: ExonCoordMapper::new(vec![
                MapperPair {
                    from_start: 1,
                    from_end: 200,
                    to_start: 25000001,
                    to_end: 25000200,
                    ori: 1,
                },
                MapperPair {
                    from_start: 201,
                    from_end: 500,
                    to_start: 25001001,
                    to_end: 25001300,
                    ori: 1,
                },
            ]),
        };
        assert_eq!(mapper.start_phase, 0);
        assert_eq!(mapper.cdna_coding_start, 155);
        assert_eq!(mapper.cdna_coding_end, 622);
        assert_eq!(mapper.exon_coord_mapper.pair_count, 2);
        assert_eq!(mapper.exon_coord_mapper.pairs.len(), 2);
    }

    #[test]
    fn test_transcript_vefc_default() {
        let vefc = TranscriptVEFC::default();
        assert_eq!(vefc.codon_table, 0);
        assert!(vefc.five_prime_utr.is_none());
        assert!(vefc.three_prime_utr.is_none());
        assert!(vefc.translateable_seq.is_none());
        assert!(vefc.peptide.is_none());
        assert!(vefc.introns.is_empty());
        assert!(vefc.sorted_exons.is_empty());
        assert!(vefc.mapper.is_none());
        assert!(vefc.protein_features.is_empty());
        assert!(vefc.protein_function_predictions.is_none());
        assert!(vefc.seq_edits.is_empty());
    }

    #[test]
    fn test_protein_function_predictions_default() {
        let preds = ProteinFunctionPredictions::default();
        assert!(preds.sift.is_none());
        assert!(preds.polyphen_humdiv.is_none());
        assert!(preds.polyphen_humvar.is_none());
    }

    #[test]
    fn test_translation_construction() {
        let t = Translation {
            stable_id: "ENSP00000355627".to_string(),
            version: Some(2),
            db_id: Some(12345),
            start: 1,
            end: 155,
            start_exon_index: 0,
            end_exon_index: 4,
            seq: Some("MKAVILF...".to_string()),
        };
        assert_eq!(t.stable_id, "ENSP00000355627");
        assert_eq!(t.version, Some(2));
        assert_eq!(t.start_exon_index, 0);
        assert_eq!(t.end_exon_index, 4);
    }

    #[test]
    fn test_exon_stable_id() {
        let exon = Exon {
            stable_id: Some("ENSE00001472862".to_string()),
            start: 25000001,
            end: 25000200,
            rank: 1,
            phase: 0,
            end_phase: 2,
        };
        assert_eq!(exon.stable_id.as_deref(), Some("ENSE00001472862"));
    }

    fn exon(start: u64, end: u64, phase: i8) -> Exon {
        Exon {
            stable_id: None,
            start,
            end,
            rank: 0,
            phase,
            end_phase: -1,
        }
    }

    fn intron(start: u64, end: u64) -> Intron {
        Intron {
            start,
            end,
            rank: 0,
        }
    }

    /// A bare transcript with the given strand, exons and introns; every
    /// other field is empty.
    fn transcript(strand: Strand, exons: Vec<Exon>, introns: Vec<Intron>) -> Transcript {
        Transcript {
            stable_id: "ENST00000000001".into(),
            version: None,
            db_id: None,
            gene_stable_id: "ENSG00000000001".into(),
            chr: "1".into(),
            start: exons.iter().map(|e| e.start).min().unwrap_or(1),
            end: exons.iter().map(|e| e.end).max().unwrap_or(1),
            strand,
            biotype: "protein_coding".into(),
            source: "Ensembl".into(),
            description: None,
            gene_symbol: None,
            gene_symbol_source: None,
            hgnc_id: None,
            gene_phenotype: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            refseq: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            exons,
            introns,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            translation_start: None,
            translation_end: None,
            translation: None,
            cdna_sequence: None,
            protein_sequence: None,
            flags: Arc::from([]),
            gencode_primary: false,
            attributes: vec![],
            vefc: None,
            derived: Default::default(),
        }
    }

    #[test]
    fn facts_frameshift_intron_reads_the_intron_list() {
        let short = transcript(
            Strand::Forward,
            vec![exon(100, 200, -1), exon(213, 300, -1)],
            vec![intron(201, 212)],
        );
        assert!(short.facts().has_frameshift_intron);
        assert!(short.facts().vefc_has_frameshift_intron);

        let long = transcript(
            Strand::Forward,
            vec![exon(100, 200, -1), exon(215, 300, -1)],
            vec![intron(201, 214)],
        );
        assert!(!long.facts().has_frameshift_intron);
        assert!(!long.facts().vefc_has_frameshift_intron);
    }

    #[test]
    fn facts_frameshift_intron_falls_back_to_exon_gaps_without_an_intron_list() {
        let tx = transcript(
            Strand::Reverse,
            vec![exon(213, 300, -1), exon(100, 200, -1)],
            vec![],
        );
        assert!(tx.facts().has_frameshift_intron);
        assert!(!tx.facts().vefc_has_frameshift_intron);
    }

    #[test]
    fn facts_vefc_intron_list_is_read_before_the_transcript_list() {
        let mut tx = transcript(
            Strand::Forward,
            vec![exon(100, 200, -1), exon(215, 300, -1)],
            vec![intron(201, 214)],
        );
        tx.vefc = Some(TranscriptVEFC {
            introns: vec![intron(201, 212)],
            ..Default::default()
        });
        assert!(!tx.facts().has_frameshift_intron);
        assert!(tx.facts().vefc_has_frameshift_intron);
    }

    #[test]
    fn facts_read_flags_and_the_coding_model() {
        let mut tx = transcript(Strand::Forward, vec![exon(100, 200, 0)], vec![]);
        let facts = *tx.facts();
        assert!(!facts.cds_start_nf);
        assert!(!facts.cds_end_nf);
        assert!(!facts.has_coding_model);

        tx.flags = Arc::from(["cds_end_NF".to_string()]);
        tx.translation = Some(Translation {
            stable_id: "ENSP00000000001".into(),
            version: None,
            db_id: None,
            start: 1,
            end: 10,
            start_exon_index: 0,
            end_exon_index: 0,
            seq: None,
        });
        tx.vefc = Some(TranscriptVEFC {
            mapper: Some(TranscriptMapper {
                start_phase: 0,
                cdna_coding_start: 10,
                cdna_coding_end: 60,
                exon_coord_mapper: ExonCoordMapper::new(vec![]),
            }),
            ..Default::default()
        });
        tx.derived = Default::default();
        let facts = *tx.facts();
        assert!(!facts.cds_start_nf);
        assert!(facts.cds_end_nf);
        assert!(facts.has_coding_model);
    }

    #[test]
    fn facts_start_phase_offset_is_the_first_exon_by_strand() {
        let forward = transcript(
            Strand::Forward,
            vec![exon(100, 200, 2), exon(300, 400, 0)],
            vec![intron(201, 299)],
        );
        assert_eq!(forward.facts().start_phase_offset, 2);
        let reverse = transcript(
            Strand::Reverse,
            vec![exon(100, 200, 2), exon(300, 400, 1)],
            vec![intron(201, 299)],
        );
        assert_eq!(reverse.facts().start_phase_offset, 1);
        let negative = transcript(Strand::Forward, vec![exon(100, 200, -1)], vec![]);
        assert_eq!(negative.facts().start_phase_offset, 0);
    }

    #[test]
    fn feature_ordinals_follow_transcript_order() {
        let tx = transcript(
            Strand::Reverse,
            vec![exon(100, 200, -1), exon(500, 600, -1), exon(300, 400, -1)],
            vec![intron(201, 299), intron(401, 499)],
        );
        let ordinals = tx.feature_ordinals();
        assert_eq!(&*ordinals.exons, &[(500, 600), (300, 400), (100, 200)]);
        assert_eq!(&*ordinals.introns, &[(401, 499), (201, 299)]);

        let derived = transcript(
            Strand::Forward,
            vec![exon(300, 400, -1), exon(100, 200, -1)],
            vec![],
        );
        assert_eq!(&*derived.feature_ordinals().introns, &[(201, 299)]);
    }

    /// The empty cells every loaded transcript carries, touched or not.
    #[test]
    fn derived_cells_cost_at_most_32_bytes_per_transcript() {
        assert!(std::mem::size_of::<TranscriptDerived>() <= 32);
    }

    #[test]
    fn derived_cache_is_per_value_and_empty_after_clone() {
        let tx = transcript(
            Strand::Forward,
            vec![exon(100, 200, -1), exon(215, 300, -1)],
            vec![intron(201, 214)],
        );
        assert!(!tx.facts().has_frameshift_intron);
        let mut copy = tx.clone();
        copy.introns = vec![intron(201, 212)];
        assert!(copy.facts().has_frameshift_intron);
        assert!(!tx.facts().has_frameshift_intron);
    }
}
