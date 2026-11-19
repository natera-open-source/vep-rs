// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Transcript and related types for cached annotation data.
//!
//! Defines the [`Transcript`] type and its components (exons, introns, translation,
//! coordinate mapper, VEFC cache) as loaded from the VEP JSON or binary cache.
//! These types are the primary input to the consequence calculator in `vep-effects`.

use std::borrow::Cow;
use std::sync::{Arc, OnceLock};

use crate::coordinate::Strand;
use crate::sv_routing::perl_split;

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
    /// Gene stable ID version.
    pub gene_version: Option<u32>,
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
    /// Gene symbol (e.g., "SOD1").
    pub gene_symbol: Option<Arc<str>>,
    /// Gene symbol source (e.g., "HGNC").
    pub gene_symbol_source: Option<Arc<str>>,
    /// The primary identifier of the EntrezGene entry a RefSeq gene's symbol is
    /// taken from, Ensembl VEP's `_gene_symbol_id` key
    /// (`AnnotationSource/Database/Transcript.pm:563` at release/116.2); `None`
    /// on an Ensembl cache's transcripts and on a RefSeq gene without exactly
    /// one EntrezGene entry. Read by `--transcript_filter` alone.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_symbol_id: Option<Arc<str>>,
    /// HGNC ID.
    pub hgnc_id: Option<Arc<str>>,
    /// Gene phenotype association count (number of phenotypes associated with
    /// this gene; values can exceed 255 for highly-studied genes like TP53,
    /// hence u32 not u8).
    pub gene_phenotype: Option<u32>,
    /// Whether this is the canonical transcript for the gene.
    pub canonical: bool,
    /// The `MANE_Select` attribute's value (a RefSeq accession), `Some` when the
    /// attribute is present, empty when it carries no value.
    pub mane_select: Option<String>,
    /// The `MANE_Plus_Clinical` attribute's value, `Some` on presence as above.
    pub mane_plus_clinical: Option<String>,
    /// Transcript Support Level: the digit run after `tsl` in a `TSL`
    /// attribute's value, kept as the text it was matched as ([`tsl_level`]), since
    /// Ensembl VEP prints that text (`1`; `007` from a GFF3 that wrote it) and
    /// reads it as a number only in the JSON format and the pick orders. `None`
    /// for `tslNA`, for the one digit `0`, and without the attribute.
    pub tsl: Option<String>,
    /// APPRIS annotation as the cache spells it (`principal1`, `alternativeM`).
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
    /// UniProt isoform accession(s).
    pub uniprot_isoform: Option<String>,
    /// Exons.
    pub exons: Vec<Exon>,
    /// The exons' stable ids in rank order, loaded only for a run whose plugins
    /// print them: no output or consequence reads an exon's id, and parsing
    /// them is a cost of every load.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exon_ids: Option<Arc<[Arc<str>]>>,
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
    /// A `gencode_basic` attribute is present.
    pub gencode_basic: bool,
    /// A `gencode_primary` attribute is present.
    pub gencode_primary: bool,
    /// Transcript attributes (code/value pairs from the Perl cache).
    pub attributes: Vec<Attribute>,
    /// The cache the transcript came from in a merged cache: the `_source_cache`
    /// key stamped on the transcripts of a region file both caches hold
    /// (`Pipeline/DumpVEP/MergeVEP.pm:103, 110` at release/116.2); `None` on a
    /// single-source cache and on a region file one cache alone held.
    pub source_cache: Option<SourceCache>,
    /// Whether the model was corrected against the cache's alignment file: the
    /// `_bam_edit_status` key `apply_edits` sets
    /// (`AnnotationType/Transcript.pm:601, 621` at release/116.2); `None` on a
    /// transcript the alignment never touched and on every cache without one.
    pub bam_edit_status: Option<BamEditStatus>,
    /// Cached data for variant effect calculation (_variation_effect_feature_cache).
    pub vefc: Option<TranscriptVEFC>,
    /// Facts derived from the fields above, filled on first use and read by
    /// every later row against this transcript ([`Transcript::facts`],
    /// [`Transcript::feature_ordinals`], [`Transcript::spliced_seq`]).
    #[serde(skip)]
    pub derived: TranscriptDerived,
    /// Position of the transcript's Ensembl VEP sort key, the stable id followed
    /// by `_` and the dbID when the record carries one
    /// (`BaseVariationFeature.pm:310-321`), in the byte order of the
    /// chromosome's keys, counted from 1 and assigned by the loader once the
    /// chromosome is complete: two transcripts of one chromosome compare by rank
    /// as their keys compare, and equal keys share a rank. Never serialized; 0
    /// until assigned.
    #[serde(skip)]
    pub id_rank: u32,
}

/// The cache a transcript of a merged cache came from: Ensembl VEP's
/// `_source_cache` key, `Ensembl` or `RefSeq` (`Pipeline/DumpVEP/MergeVEP.pm`
/// at release/116.2), spelled as the cache spells it; or the short name of the
/// custom annotation source a transcript was read from, which the `SOURCE`
/// column prints alike (OutputFactory.pm:1583, Constants.pm:98 and :101).
/// Serialized as that text: `"Ensembl"`, `"RefSeq"`, or the custom name.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SourceCache {
    Ensembl,
    RefSeq,
    Custom(Arc<str>),
}

impl SourceCache {
    /// The `SOURCE` column's text.
    pub fn column_text(&self) -> &str {
        match self {
            SourceCache::Ensembl => "Ensembl",
            SourceCache::RefSeq => "RefSeq",
            SourceCache::Custom(name) => name,
        }
    }
}

impl serde::Serialize for SourceCache {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(self.column_text())
    }
}

impl<'de> serde::Deserialize<'de> for SourceCache {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let text: String = serde::Deserialize::deserialize(deserializer)?;
        Ok(match text.as_str() {
            "Ensembl" => SourceCache::Ensembl,
            "RefSeq" => SourceCache::RefSeq,
            name => SourceCache::Custom(Arc::from(name)),
        })
    }
}

/// The outcome of editing a transcript model against the cache's alignment
/// file: Ensembl VEP's `_bam_edit_status`, `ok` when the edited sequence equals
/// the aligned one and the model carries the edits, `failed` when it did not and
/// the model was dumped unedited (`AnnotationType/Transcript.pm:599-627` at
/// release/116.2), spelled as the cache spells it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum BamEditStatus {
    Ok,
    Failed,
}

impl BamEditStatus {
    /// The `BAM_EDIT` column's text, the status upper-cased
    /// (`uc($tr->{_bam_edit_status})`, OutputFactory.pm:1551-1553 at release/116.2).
    pub fn column_text(self) -> &'static str {
        match self {
            BamEditStatus::Ok => "OK",
            BamEditStatus::Failed => "FAILED",
        }
    }
}

/// The lazily filled cells behind [`Transcript::facts`],
/// [`Transcript::feature_ordinals`], [`Transcript::spliced_seq`] and
/// [`Transcript::refseq_alignment`]. A fresh or cloned transcript starts empty,
/// so a transcript edited after construction is read as edited.
#[derive(Debug, Default)]
pub struct TranscriptDerived {
    facts: OnceLock<TranscriptFacts>,
    ordinals: OnceLock<Box<FeatureOrdinals>>,
    spliced_seq: OnceLock<Option<Box<str>>>,
    refseq_alignment: OnceLock<Option<Arc<RefSeqAlignment>>>,
}

/// What a RefSeq model's attributes record of its alignment to the genome,
/// read once per transcript: the match codes the `REFSEQ_MATCH` column lists
/// and the `_rna_edit` entries, which Ensembl VEP reads per allele from the
/// attribute list (`TranscriptVariationAllele::hgvs_transcript` and
/// `::codon`, ensembl-variation release/116). A transcript carrying neither
/// kind of attribute has none ([`Transcript::refseq_alignment`] is `None`).
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct RefSeqAlignment {
    /// The codes of the attributes matching `^rseq`, in attribute order
    /// (`grep {$_->code =~ /^rseq/}` over the attribute list, OutputFactory.pm:1547 at
    /// release/116.2); `enst_refseq_compare` and `_rna_edit` are not among them.
    pub match_codes: Arc<[Box<str>]>,
    /// The `_rna_edit` attributes that change the cDNA's length, in order of
    /// their start, as [`get_misalignment_offset`](Self::misalignment_offset)
    /// reads them.
    pub indel_edits: Box<[IndelEdit]>,
    /// The bases of the first `_rna_edit` attribute in attribute order, its
    /// `<start> <end> ` prefix removed (`s/\d+\s\d+\s//`,
    /// TranscriptVariationAllele.pm:832-833): empty for a deletion; `None` on a
    /// transcript without an edit.
    pub first_edit_bases: Option<Box<str>>,
}

/// One `_rna_edit` attribute that inserts or deletes bases against the genome,
/// read as `get_misalignment_offset` reads the attribute's `<start> <end>
/// <bases>` value (TranscriptVariationAllele.pm:1605-1618): split on single
/// spaces with trailing empty fields dropped ([`perl_split`]), so a deletion's
/// value `50 52 ` has two fields and an insertion's `10 9 AAA` three.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct IndelEdit {
    /// The edit's start in the unedited cDNA, the first field.
    pub start: i64,
    /// The bases the edit adds to the cDNA numbering downstream of it: the
    /// length of the third field for an edit with three fields
    /// (`length($split_val[2])`, :1617), `-1 - (end - start)` for one with two
    /// (:1614).
    pub length_change: i64,
}

impl RefSeqAlignment {
    /// The transcript's alignment facts, or `None` when its attributes carry
    /// no `rseq` code and no `_rna_edit`.
    fn read(attributes: &[Attribute]) -> Option<Self> {
        let match_codes: Arc<[Box<str>]> = attributes
            .iter()
            .filter(|a| a.code.starts_with("rseq"))
            .map(|a| Box::from(&*a.code))
            .collect();
        let edits: Vec<&Attribute> = attributes
            .iter()
            .filter(|a| a.code.starts_with("_rna_edit"))
            .collect();
        if match_codes.is_empty() && edits.is_empty() {
            return None;
        }
        let first_edit_bases = edits
            .first()
            .map(|edit| Box::from(strip_edit_span(&edit.value)));
        let mut indel_edits: Vec<IndelEdit> = edits.into_iter().filter_map(indel_edit).collect();
        indel_edits.sort_by_key(|edit| edit.start);
        Some(RefSeqAlignment {
            match_codes,
            indel_edits: indel_edits.into_boxed_slice(),
            first_edit_bases,
        })
    }

    /// `get_misalignment_offset` (TranscriptVariationAllele.pm:1594-1624): the
    /// sum of the length changes of the indel edits whose start is at or before
    /// `variant_cdna_start`, the allele's cDNA start, 0 standing for an
    /// undefined one (:1608-1610) so that an allele outside the exons counts no
    /// edit.
    pub fn misalignment_offset(&self, variant_cdna_start: Option<i64>) -> i64 {
        let location = variant_cdna_start.unwrap_or(0);
        self.indel_edits
            .iter()
            .filter(|edit| location >= edit.start)
            .map(|edit| edit.length_change)
            .sum()
    }
}

/// The attribute value past its `<digits><space><digits><space>` prefix, the
/// bases of the edit, which `s/\d+\s\d+\s//` leaves (TranscriptVariationAllele.pm:833);
/// a value not of that shape is kept whole, as the substitution then matches
/// nothing.
fn strip_edit_span(value: &str) -> &str {
    let is_digits = |field: &str| !field.is_empty() && field.bytes().all(|b| b.is_ascii_digit());
    let Some((start, rest)) = value.split_once(|c: char| c.is_ascii_whitespace()) else {
        return value;
    };
    let Some((end, bases)) = rest.split_once(|c: char| c.is_ascii_whitespace()) else {
        return value;
    };
    if is_digits(start) && is_digits(end) {
        bases
    } else {
        value
    }
}

/// The indel an `_rna_edit` attribute describes, or `None` for a substitution,
/// whose three fields replace as many bases as they span (`:1606`). A value
/// without two leading integers, which no cache writes, is skipped here;
/// Ensembl VEP reads a missing or non-numeric field as 0 and applies the same
/// arithmetic, counting an empty value as a deletion of one (`:1614`).
fn indel_edit(attribute: &Attribute) -> Option<IndelEdit> {
    let fields: Vec<&str> = perl_split(&attribute.value, &[' ']).collect();
    let start: i64 = fields.first()?.parse().ok()?;
    let end: i64 = fields.get(1)?.parse().ok()?;
    let length_change = if fields.len() == 3 {
        let bases = fields[2].len() as i64;
        if end - start + 1 == bases {
            return None;
        }
        bases
    } else {
        -1 - (end - start)
    };
    Some(IndelEdit {
        start,
        length_change,
    })
}

impl Clone for TranscriptDerived {
    fn clone(&self) -> Self {
        Self::default()
    }
}

/// Per-transcript facts that every row against the transcript reads; the
/// counterparts of what Ensembl VEP stores on the transcript's
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
    /// The biotype is `miRNA`, the one biotype whose `miRNA` attributes name
    /// mature products.
    pub is_mirna: bool,
    /// The biotype is `protein_coding`, the one biotype whose rows take
    /// `coding_transcript_variant` and never the non-coding context term.
    pub is_protein_coding: bool,
    /// The biotype is `nonsense_mediated_decay` in any ASCII case, the one
    /// biotype whose rows take `NMD_transcript_variant`.
    pub is_nmd: bool,
    /// Neither `is_protein_coding` nor `is_nmd`: every other biotype, whose
    /// rows take `non_coding_transcript_variant` as the context term.
    pub is_non_coding_context: bool,
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
            is_mirna: &*tx.biotype == "miRNA",
            is_protein_coding: tx.is_protein_coding(),
            is_nmd: tx.is_nmd_transcript(),
            is_non_coding_context: tx.is_non_coding_context_transcript(),
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

/// An exon within a transcript: its span, rank and phases. The cache's exon
/// stable id is not loaded onto the exon, since no output column or
/// consequence reads it; a run whose plugins print the ids loads them onto
/// [`Transcript::exon_ids`].
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Exon {
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

/// The support level a `TSL` attribute's value carries, as Ensembl VEP reads
/// it for the `TSL` field (`OutputFactory.pm:1610-1611` at release/116.2) and
/// for the pick orders (`:752-753`; `VariationFeatureAdaptor.pm:2182-2183` at
/// ensembl-variation release/116): `m/tsl(\d+)/` then `$1 if $1`. The match is
/// the first `tsl` a digit follows and the run of digits after it, kept as
/// text, so `tsl1 (assigned to previous version 8)` is `1`, `tslNA then tsl4`
/// is `4` and `tsl007` is `007`; `None` when no `tsl` is followed by a digit
/// (`tslNA`) and when the run is the one digit `0`, the only digit run Perl's
/// truth test drops (`00` is true).
pub fn tsl_level(value: &str) -> Option<&str> {
    let digits = value.match_indices("tsl").find_map(|(at, tsl)| {
        let rest = &value[at + tsl.len()..];
        let digit_run_len = rest
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(rest.len());
        (digit_run_len > 0).then_some(&rest[..digit_run_len])
    })?;
    (digits != "0").then_some(digits)
}

/// A code/value attribute pair from the Perl cache.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Attribute {
    /// The attribute code (`TSL`, `gencode_basic`, `_rna_edit`): one of a few
    /// dozen strings a cache repeats on most of its transcripts, so the loader
    /// shares one allocation per distinct code.
    pub code: Arc<str>,
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
    /// The spliced cDNA with the model's edits applied, which `apply_edits`
    /// stores on an `ok` BAM-edited transcript alone
    /// (`AnnotationType/Transcript.pm:617` at release/116.2); every other
    /// transcript's cDNA is assembled from the pieces above or read from the
    /// FASTA.
    pub spliced_seq: Option<String>,
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
    /// that constructs/deserializes `pairs` without going through [`Self::new`].
    pub fn rebuild_genomic_index(&mut self) {
        self.by_genomic_start = Self::build_genomic_index(&self.pairs);
    }

    /// The genomic pieces of the cDNA range `start..=end`, as
    /// `Bio::EnsEMBL::Mapper::map_coordinates` lists them for the forward
    /// strand (Mapper.pm, `TranscriptMapper::cdna2genomic`): for each pair the
    /// range overlaps, in cDNA order, a [`MappedPiece::Gap`] over the cDNA
    /// bases before the pair that no pair covers, then a
    /// [`MappedPiece::Coordinate`] over the bases the pair places on the
    /// genome, on the pair's strand; a gap after the last pair when the range
    /// runs past it; one gap over the whole range when no pair overlaps it. A
    /// range of the shape `end + 1 .. end` is an insertion point, mapped as
    /// `map_insert` maps it: the two bases around it are mapped and the piece
    /// is turned around (`start = end + 1` on the genome too); when the point
    /// falls between two pairs each piece is moved to its side of the point
    /// and a gap among them is dropped. A range whose start lies more than one
    /// past its end is one gap: the mapper's interval arithmetic has no
    /// reading for it.
    pub fn cdna2genomic(&self, start: i64, end: i64) -> Vec<MappedPiece> {
        if start == end + 1 {
            return self.map_insert(start, end);
        }
        if start > end {
            return vec![MappedPiece::Gap { start, end }];
        }
        let mut pieces = Vec::new();
        let mut cursor = start;
        let mut last_used: Option<&MapperPair> = None;
        let overlapping = self
            .pairs
            .iter()
            .filter(|pair| pair.from_start as i64 <= end && pair.from_end as i64 >= start);
        for pair in overlapping {
            let (from_start, from_end) = (pair.from_start as i64, pair.from_end as i64);
            let (to_start, to_end) = (pair.to_start as i64, pair.to_end as i64);
            if cursor < from_start {
                pieces.push(MappedPiece::Gap {
                    start: cursor,
                    end: from_start - 1,
                });
                cursor = from_start;
            }
            let forward = pair.ori == 1;
            let enveloped = end > from_end;
            let (target_start, target_end) = if forward {
                let target_end = if enveloped {
                    to_end
                } else {
                    to_start + (end - from_start)
                };
                (to_start + (cursor - from_start), target_end)
            } else {
                let target_start = if enveloped {
                    to_start
                } else {
                    to_end - (end - from_start)
                };
                (target_start, to_end - (cursor - from_start))
            };
            pieces.push(MappedPiece::Coordinate {
                start: target_start,
                end: target_end,
                strand: if forward {
                    Strand::Forward
                } else {
                    Strand::Reverse
                },
            });
            last_used = Some(pair);
            cursor = from_end + 1;
        }
        match last_used {
            None => pieces.push(MappedPiece::Gap { start, end }),
            Some(pair) if (pair.from_end as i64) < end => pieces.push(MappedPiece::Gap {
                start: pair.from_end as i64 + 1,
                end,
            }),
            Some(_) => {}
        }
        pieces
    }

    /// `Mapper::map_insert` for the forward strand: the two bases around the
    /// insertion point `start = end + 1` are mapped, a single piece is turned
    /// around, and two pieces are each moved off the point onto their own side,
    /// a gap among them dropped. Any other count of pieces is returned as
    /// mapped, where the Perl mapper throws.
    fn map_insert(&self, start: i64, end: i64) -> Vec<MappedPiece> {
        let mut pieces = self.cdna2genomic(end, start);
        match pieces.as_mut_slice() {
            [MappedPiece::Coordinate { start, end, .. }] | [MappedPiece::Gap { start, end }] => {
                std::mem::swap(start, end);
                pieces
            }
            [first, second] => {
                let mut moved = Vec::with_capacity(2);
                if let MappedPiece::Coordinate { start, end, strand } = first {
                    match strand {
                        Strand::Reverse => *end -= 1,
                        Strand::Forward => *start += 1,
                    }
                    moved.push(*first);
                }
                if let MappedPiece::Coordinate { start, end, strand } = second {
                    match strand {
                        Strand::Reverse => *start += 1,
                        Strand::Forward => *end -= 1,
                    }
                    moved.push(*second);
                }
                moved
            }
            _ => pieces,
        }
    }
}

/// One piece of a cDNA range mapped onto the genome by
/// [`ExonCoordMapper::cdna2genomic`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MappedPiece {
    /// `Bio::EnsEMBL::Mapper::Coordinate`: the genomic bases `start..=end` on
    /// `strand`; an insertion point has `start = end + 1`.
    Coordinate {
        start: i64,
        end: i64,
        strand: Strand,
    },
    /// `Bio::EnsEMBL::Mapper::Gap`: the cDNA bases `start..=end` no pair of the
    /// mapper covers.
    Gap { start: i64, end: i64 },
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
/// - Gzip-compressed blob, held compressed and inflated on the first lookup
/// - Header: 3 bytes `"VEP"`
/// - Per position: 20 amino acids × 2 bytes = 40 bytes
/// - Each cell: `u16` little-endian; top 2 bits (15-14) = category index, bits 13-10
///   unused, bottom 10 = score × 1000
/// - `0xFFFF` = no prediction
/// - Amino acid order: A C D E F G H I K L M N P Q R S T V W Y
///
/// A run without `--sift` or `--polyphen` never looks a matrix up, so every
/// stream stays as the cache stores it; a run that does pays the inflate once
/// per matrix it looks up, on the thread of the first lookup.
///
/// Serialized as `analysis`, `sub_analysis`, `peptide_length` and the inflated
/// bytes under `predictions_data`, whichever way the matrix was built.
///
/// See [`crate::prediction`] for the decoder implementation and format documentation.
#[derive(Debug, Clone)]
pub struct PredictionMatrix {
    /// Analysis name: "sift", "polyphen_humdiv", or "polyphen_humvar".
    pub analysis: String,
    /// Optional sub-analysis identifier.
    pub sub_analysis: Option<String>,
    /// Length of the protein (number of residues).
    pub peptide_length: usize,
    /// The gzip stream of the matrix as the cache stores it; empty for a matrix
    /// built from inflated bytes.
    compressed: Box<[u8]>,
    /// The inflated matrix (header + position data), filled on the first
    /// lookup; `None` inside when the stream did not inflate.
    inflated: OnceLock<Option<Box<[u8]>>>,
}

/// The serialized shape of a [`PredictionMatrix`].
#[derive(serde::Serialize, serde::Deserialize)]
struct PredictionMatrixRepr<'a> {
    analysis: Cow<'a, str>,
    sub_analysis: Option<Cow<'a, str>>,
    peptide_length: usize,
    predictions_data: Cow<'a, [u8]>,
}

impl serde::Serialize for PredictionMatrix {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        PredictionMatrixRepr {
            analysis: Cow::Borrowed(self.analysis.as_str()),
            sub_analysis: self.sub_analysis.as_deref().map(Cow::Borrowed),
            peptide_length: self.peptide_length,
            predictions_data: Cow::Borrowed(self.predictions_data().unwrap_or_default()),
        }
        .serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for PredictionMatrix {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let repr = PredictionMatrixRepr::deserialize(deserializer)?;
        Ok(Self::from_inflated(
            repr.analysis.into_owned(),
            repr.sub_analysis.map(Cow::into_owned),
            repr.peptide_length,
            repr.predictions_data.into_owned(),
        ))
    }
}

impl PredictionMatrix {
    /// A matrix from the gzip stream the cache stores, inflated on the first
    /// lookup.
    pub fn from_compressed(
        analysis: String,
        sub_analysis: Option<String>,
        peptide_length: usize,
        compressed: Vec<u8>,
    ) -> Self {
        Self {
            analysis,
            sub_analysis,
            peptide_length,
            compressed: compressed.into_boxed_slice(),
            inflated: OnceLock::new(),
        }
    }

    /// A matrix from its inflated bytes (header + position data).
    pub fn from_inflated(
        analysis: String,
        sub_analysis: Option<String>,
        peptide_length: usize,
        inflated: Vec<u8>,
    ) -> Self {
        Self {
            analysis,
            sub_analysis,
            peptide_length,
            compressed: Box::default(),
            inflated: OnceLock::from(Some(inflated.into_boxed_slice())),
        }
    }

    /// The inflated matrix bytes (header + position data), inflating the
    /// stored stream on the first call; `None` when the stream does not
    /// inflate, which [`Self::lookup`] reads as no prediction anywhere. Use
    /// [`crate::prediction::decode_prediction`] to look up values.
    pub fn predictions_data(&self) -> Option<&[u8]> {
        self.inflated
            .get_or_init(|| {
                crate::prediction::gunzip_matrix(&self.compressed)
                    .ok()
                    .map(Vec::into_boxed_slice)
            })
            .as_deref()
    }

    /// Whether the inflate has run, on a lookup or at construction from inflated
    /// bytes; a stream that failed to inflate counts as run.
    pub fn is_inflated(&self) -> bool {
        self.inflated.get().is_some()
    }

    /// Look up the prediction for a given protein position and amino acid.
    ///
    /// Returns `None` if the position/AA is out of range or no prediction exists.
    pub fn lookup(
        &self,
        position: usize,
        amino_acid: u8,
        analysis: crate::prediction::AnalysisType,
    ) -> Option<crate::prediction::Prediction> {
        crate::prediction::decode_prediction(
            self.predictions_data()?,
            position,
            amino_acid,
            analysis,
        )
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

/// The empty transcript: empty identifiers, chromosome, biotype and source,
/// every optional value absent, the forward strand, a span of 0 to 0, no
/// exons, introns, flags or attributes, and rank 0.
impl Default for Transcript {
    fn default() -> Self {
        Self {
            stable_id: Arc::from(""),
            version: None,
            db_id: None,
            gene_stable_id: Arc::from(""),
            gene_version: None,
            chr: String::new(),
            start: 0,
            end: 0,
            strand: Strand::Forward,
            biotype: Arc::from(""),
            source: String::new(),
            gene_symbol: None,
            gene_symbol_source: None,
            gene_symbol_id: None,
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
            uniprot_isoform: None,
            exons: Vec::new(),
            exon_ids: None,
            introns: Vec::new(),
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
            gencode_basic: false,
            gencode_primary: false,
            attributes: Vec::new(),
            source_cache: None,
            bam_edit_status: None,
            vefc: None,
            derived: TranscriptDerived::default(),
            id_rank: 0,
        }
    }
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

    /// The spliced cDNA the cache holds for the transcript: the edited cDNA a
    /// BAM-edited `ok` model stores whole ([`TranscriptVEFC::spliced_seq`], read
    /// before any assembly as `TranscriptVariation::_get_ref_allele` reads it,
    /// TranscriptVariation.pm:462 at release/116), else a coding transcript's
    /// cDNA assembled on the first call from the cached 5' UTR, translateable
    /// sequence and 3' UTR and shared by every later one; `None` for any other
    /// transcript without a coding model or whose cached pieces do not fit
    /// together.
    ///
    /// The cache holds the three pieces Ensembl's `Transcript::translateable_seq`
    /// and the two UTR accessors cut from one spliced sequence, so their
    /// concatenation is that sequence once the `N` padding the translateable
    /// sequence starts with (one `N` per base of the start exon's phase) is
    /// removed. The padding's length is the excess of the cached sequence over
    /// the coding span, so no exon rule is needed to find it; a cache whose 5'
    /// UTR does not end at the coding start, or whose padding is not 0 to 2
    /// bases, yields `None`. An unedited non-coding transcript's cDNA is not in
    /// the cache.
    pub fn spliced_seq(&self) -> Option<&str> {
        if let Some(edited) = self.vefc.as_ref().and_then(|v| v.spliced_seq.as_deref()) {
            return Some(edited);
        }
        self.derived
            .spliced_seq
            .get_or_init(|| self.assemble_spliced_seq().map(String::into_boxed_str))
            .as_deref()
    }

    /// What the transcript's attributes record of a RefSeq model's alignment to
    /// the genome ([`RefSeqAlignment`]), read on the first call and shared by
    /// every later one; `None` for a transcript whose attributes carry no
    /// `rseq` code and no `_rna_edit`, every unedited Ensembl model included.
    pub fn refseq_alignment(&self) -> Option<&Arc<RefSeqAlignment>> {
        self.derived
            .refseq_alignment
            .get_or_init(|| RefSeqAlignment::read(&self.attributes).map(Arc::new))
            .as_ref()
    }

    /// The `rseq` match codes of [`Transcript::refseq_alignment`], empty without
    /// one.
    pub fn rseq_codes(&self) -> &[Box<str>] {
        self.refseq_alignment()
            .map_or(&[], |alignment| &alignment.match_codes)
    }

    /// The indel edits of [`Transcript::refseq_alignment`], empty without one.
    pub fn indel_edits(&self) -> &[IndelEdit] {
        self.refseq_alignment()
            .map_or(&[], |alignment| &alignment.indel_edits)
    }

    fn assemble_spliced_seq(&self) -> Option<String> {
        if !self.facts().has_coding_model {
            return None;
        }
        let vefc = self.vefc.as_ref()?;
        let mapper = vefc.mapper.as_ref()?;
        let coding = vefc.translateable_seq.as_deref()?;
        let five = vefc.five_prime_utr.as_deref().unwrap_or("");
        let three = vefc.three_prime_utr.as_deref().unwrap_or("");
        let coding_len =
            usize::try_from(mapper.cdna_coding_end - mapper.cdna_coding_start + 1).ok()?;
        let padding = coding.len().checked_sub(coding_len)?;
        if padding > 2 || five.len() as u64 != mapper.cdna_coding_start - 1 {
            return None;
        }
        let mut seq = String::with_capacity(five.len() + coding_len + three.len());
        seq.push_str(five);
        seq.push_str(&coding[padding..]);
        seq.push_str(three);
        Some(seq)
    }

    /// Length of the cached translateable sequence of a transcript with a
    /// translation, padding included
    /// (`length($tr->{_variation_effect_feature_cache}->{translateable_seq})`,
    /// OutputFactory.pm:743-748 at release/116.2); `None` without a translation.
    pub fn cds_length(&self) -> Option<u32> {
        self.translation.is_some().then(|| {
            self.vefc
                .as_ref()
                .and_then(|v| v.translateable_seq.as_deref())
                .map_or(0, |s| s.len() as u32)
        })
    }

    /// The sum of the exon lengths (`Bio::EnsEMBL::Transcript::length`).
    pub fn spliced_length(&self) -> u32 {
        self.exons
            .iter()
            .map(|e| (e.end.saturating_sub(e.start) + 1) as u32)
            .sum()
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
    /// regions even though they are targeted for degradation. Ensembl VEP's SV
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
    /// NMD is a separate context term (`NMD_transcript_variant`), mirroring
    /// Ensembl VEP output, so NMD transcripts are excluded here.
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

    /// A SIFT matrix built from its gzip stream inflates on the first lookup
    /// and once only; a stream that is not gzip reads as no prediction
    /// anywhere; a matrix built from inflated bytes answers without a stream;
    /// both constructions serialize to the inflated bytes and read back whole.
    #[test]
    fn prediction_matrix_inflates_on_first_lookup() {
        use crate::prediction::AnalysisType;
        use std::io::Write as _;

        let mut raw = b"VEP".to_vec();
        raw.extend(std::iter::repeat_n(0xFFu8, 2 * 20 * 2));
        // Residue 2, amino acid C (column 1): category 1 (deleterious), score 0.01.
        let offset = 3 + (20 + 1) * 2;
        raw[offset..offset + 2].copy_from_slice(&((1u16 << 14) | 10).to_le_bytes());
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(&raw).unwrap();
        let compressed = gz.finish().unwrap();

        let lazy = PredictionMatrix::from_compressed("sift".to_string(), None, 2, compressed);
        assert!(!lazy.is_inflated());
        let hit = lazy.lookup(2, b'C', AnalysisType::Sift).unwrap();
        assert_eq!((hit.label(), hit.score), ("deleterious", 0.01));
        assert!(lazy.is_inflated());
        let first = lazy.predictions_data().unwrap().as_ptr();
        assert!(lazy.lookup(1, b'A', AnalysisType::Sift).is_none());
        assert_eq!(lazy.predictions_data().unwrap().as_ptr(), first);
        assert_eq!(lazy.predictions_data(), Some(raw.as_slice()));

        let broken =
            PredictionMatrix::from_compressed("sift".to_string(), None, 2, b"VEP".to_vec());
        assert!(broken.lookup(2, b'C', AnalysisType::Sift).is_none());
        assert!(broken.is_inflated() && broken.predictions_data().is_none());

        let eager = PredictionMatrix::from_inflated("sift".to_string(), None, 2, raw.clone());
        assert!(eager.is_inflated());
        assert_eq!(eager.predictions_data(), Some(raw.as_slice()));
        assert!(eager.lookup(2, b'C', AnalysisType::Sift).is_some());

        let lazy_json = serde_json::to_value(&lazy).unwrap();
        assert_eq!(lazy_json, serde_json::to_value(&eager).unwrap());
        assert_eq!(
            lazy_json["predictions_data"].as_array().unwrap().len(),
            raw.len()
        );
        let read_back: PredictionMatrix = serde_json::from_value(lazy_json).unwrap();
        assert_eq!(read_back.predictions_data(), Some(raw.as_slice()));
        assert_eq!(
            (read_back.analysis.as_str(), read_back.peptide_length),
            ("sift", 2)
        );
    }

    fn exon(start: u64, end: u64, phase: i8) -> Exon {
        Exon {
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
            gene_version: None,
            chr: "1".into(),
            start: exons.iter().map(|e| e.start).min().unwrap_or(1),
            end: exons.iter().map(|e| e.end).max().unwrap_or(1),
            strand,
            biotype: "protein_coding".into(),
            source: "Ensembl".into(),
            gene_symbol: None,
            gene_symbol_source: None,
            gene_symbol_id: None,
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
            uniprot_isoform: None,
            exons,
            exon_ids: None,
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
            gencode_basic: false,
            gencode_primary: false,
            attributes: vec![],
            vefc: None,
            source_cache: None,
            bam_edit_status: None,
            derived: Default::default(),
            id_rank: 0,
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

    /// The empty cells every loaded transcript carries, touched or not. Each
    /// `OnceLock` holds a `Once` beside its value: where `Once` is a 4-byte
    /// futex word (Linux) the four cells take 16 + 16 + 24 + 16 = 72 bytes, and
    /// where it is an 8-byte queue pointer (macOS) the 10-byte facts cell pads
    /// to 24, for 80.
    #[test]
    fn derived_cells_cost_at_most_72_bytes_per_transcript() {
        let bound = if std::mem::size_of::<std::sync::Once>() == 4 {
            72
        } else {
            80
        };
        assert!(std::mem::size_of::<TranscriptDerived>() <= bound);
    }

    /// A coding transcript with a translation and a mapper over `cds`, the
    /// cached translateable sequence padded with `padding` leading `N`s as
    /// Ensembl pads a start exon of that phase.
    fn coding_transcript(five: &str, cds: &str, three: &str, padding: usize) -> Transcript {
        let mut tx = transcript(
            Strand::Forward,
            vec![exon(
                100,
                100 + (five.len() + cds.len() + three.len()) as u64 - 1,
                0,
            )],
            vec![],
        );
        tx.translation = Some(Translation {
            stable_id: "ENSP00000000001".into(),
            version: None,
            db_id: None,
            start: 1,
            end: cds.len() as u64,
            start_exon_index: 0,
            end_exon_index: 0,
            seq: None,
        });
        let cdna_coding_start = five.len() as u64 + 1;
        tx.vefc = Some(TranscriptVEFC {
            five_prime_utr: (!five.is_empty()).then(|| five.to_string()),
            three_prime_utr: (!three.is_empty()).then(|| three.to_string()),
            translateable_seq: Some(format!("{}{cds}", "N".repeat(padding))),
            mapper: Some(TranscriptMapper {
                start_phase: padding as i8,
                cdna_coding_start,
                cdna_coding_end: cdna_coding_start + cds.len() as u64 - 1,
                exon_coord_mapper: ExonCoordMapper::new(vec![]),
            }),
            ..Default::default()
        });
        tx
    }

    /// The lengths the pick family compares (OutputFactory.pm:743-748): the
    /// cached translateable sequence's, padding included, for a translated
    /// transcript, and the exon sum for every transcript.
    #[test]
    fn cds_length_counts_the_padded_cached_sequence_and_spliced_length_the_exons() {
        let padded = coding_transcript("", "TGCCCTAA", "TTT", 2);
        assert_eq!(padded.cds_length(), Some(10));
        assert_eq!(padded.spliced_length(), 11);
        let mut non_coding = coding_transcript("GGAAC", "ATGCCCTAA", "TTTAAACC", 0);
        non_coding.translation = None;
        assert_eq!(non_coding.cds_length(), None);
        assert_eq!(non_coding.spliced_length(), 22);
    }

    /// The cDNA is the 5' UTR, the translateable sequence without its
    /// start-phase padding, and the 3' UTR, in that order: Ensembl's
    /// `translateable_seq` prepends one `N` per base of the start exon's phase
    /// to the coding span it cuts from `spliced_seq`.
    #[test]
    fn spliced_seq_joins_the_utrs_around_the_unpadded_coding_sequence() {
        let tx = coding_transcript("GGAAC", "ATGCCCTAA", "TTTAAACC", 0);
        assert_eq!(tx.spliced_seq(), Some("GGAACATGCCCTAATTTAAACC"));

        let padded = coding_transcript("", "TGCCCTAA", "TTT", 2);
        assert_eq!(
            padded.vefc.as_ref().unwrap().translateable_seq.as_deref(),
            Some("NNTGCCCTAA")
        );
        assert_eq!(padded.spliced_seq(), Some("TGCCCTAATTT"));

        let no_utrs = coding_transcript("", "ATGTAA", "", 0);
        assert_eq!(no_utrs.spliced_seq(), Some("ATGTAA"));
    }

    /// A transcript without a coding model has no cached cDNA, and cached
    /// pieces that do not fit together (a 5' UTR not ending at the coding
    /// start, more than two bases of padding) yield none rather than a
    /// sequence with shifted coordinates.
    #[test]
    fn spliced_seq_is_none_without_a_coding_model_or_with_pieces_that_do_not_fit() {
        let non_coding = transcript(Strand::Forward, vec![exon(100, 200, -1)], vec![]);
        assert_eq!(non_coding.spliced_seq(), None);

        let mut utr_too_short = coding_transcript("GGAAC", "ATGCCCTAA", "TTT", 0);
        utr_too_short.vefc.as_mut().unwrap().five_prime_utr = Some("GGA".into());
        assert_eq!(utr_too_short.spliced_seq(), None);

        let mut over_padded = coding_transcript("", "ATGCCCTAA", "TTT", 0);
        over_padded.vefc.as_mut().unwrap().translateable_seq = Some("NNNATGCCCTAA".into());
        assert_eq!(over_padded.spliced_seq(), None);
    }

    /// The edited cDNA a BAM-edited `ok` model stores whole is the spliced
    /// sequence, read before any assembly: on a non-coding transcript, whose
    /// pieces assemble nothing, and over the pieces of a coding one, leaving the
    /// assembly cell empty. A transcript without it assembles its cDNA from the
    /// pieces into the cell.
    #[test]
    fn spliced_seq_reads_the_cached_edited_cdna_before_assembling_one() {
        let mut non_coding = transcript(Strand::Forward, vec![exon(100, 111, -1)], vec![]);
        non_coding.vefc = Some(TranscriptVEFC {
            spliced_seq: Some("ACGTACGTAAAC".into()),
            ..Default::default()
        });
        assert_eq!(non_coding.spliced_seq(), Some("ACGTACGTAAAC"));
        assert!(non_coding.derived.spliced_seq.get().is_none());

        let mut coding = coding_transcript("GGAAC", "ATGCCCTAA", "TTTAAACC", 0);
        coding.vefc.as_mut().unwrap().spliced_seq = Some("GGAACATGCCCTAATTTAAACG".into());
        assert_eq!(coding.spliced_seq(), Some("GGAACATGCCCTAATTTAAACG"));
        assert!(coding.derived.spliced_seq.get().is_none());

        let assembled = coding_transcript("GGAAC", "ATGCCCTAA", "TTTAAACC", 0);
        assert_eq!(assembled.spliced_seq(), Some("GGAACATGCCCTAATTTAAACC"));
        assert_eq!(
            assembled
                .derived
                .spliced_seq
                .get()
                .and_then(|s| s.as_deref()),
            Some("GGAACATGCCCTAATTTAAACC")
        );
    }

    /// The two transcript keys serialize as the cache spells them, `Ensembl` or
    /// `RefSeq` for `_source_cache` and `ok` or `failed` for `_bam_edit_status`,
    /// and the `BAM_EDIT` column prints the status upper-cased.
    #[test]
    fn source_cache_and_bam_edit_status_spell_as_the_cache_does() {
        assert_eq!(
            serde_json::to_string(&SourceCache::Ensembl).unwrap(),
            r#""Ensembl""#
        );
        assert_eq!(
            serde_json::to_string(&SourceCache::RefSeq).unwrap(),
            r#""RefSeq""#
        );
        assert_eq!(
            serde_json::to_string(&BamEditStatus::Ok).unwrap(),
            r#""ok""#
        );
        assert_eq!(
            serde_json::to_string(&BamEditStatus::Failed).unwrap(),
            r#""failed""#
        );
        assert_eq!(
            serde_json::from_str::<SourceCache>(r#""RefSeq""#).unwrap(),
            SourceCache::RefSeq
        );
        assert_eq!(
            serde_json::from_str::<BamEditStatus>(r#""failed""#).unwrap(),
            BamEditStatus::Failed
        );
        assert!(serde_json::from_str::<BamEditStatus>(r#""OK""#).is_err());
        assert_eq!(BamEditStatus::Ok.column_text(), "OK");
        assert_eq!(BamEditStatus::Failed.column_text(), "FAILED");
    }

    /// The alignment cell reads the `rseq` codes in attribute order, leaving
    /// `enst_refseq_compare` and `_rna_edit` out (OutputFactory.pm:1547), and
    /// the `_rna_edit` values as `get_misalignment_offset` splits them
    /// (TranscriptVariationAllele.pm:1605-1618): `10 9 AAA` inserts three bases,
    /// `50 52 ` deletes three, `70 70 T` replaces one base with one and counts
    /// for nothing. The offset sums the edits starting at or before the allele's
    /// cDNA start, none for an undefined start (:1608-1610); the first edit's
    /// bases are kept for the mismatch rule (:832-833); a transcript with
    /// neither kind of attribute has no alignment and empty lists.
    #[test]
    fn refseq_alignment_reads_the_match_codes_and_the_indel_edits() {
        let mut tx = transcript(Strand::Forward, vec![exon(100, 200, -1)], vec![]);
        assert!(tx.refseq_alignment().is_none());
        assert!(tx.rseq_codes().is_empty() && tx.indel_edits().is_empty());

        let attribute = |code: &str, value: &str| Attribute {
            code: code.into(),
            value: value.to_string(),
        };
        tx.attributes = vec![
            attribute("_rna_edit", "50 52 "),
            attribute("rseq_5p_mismatch", ""),
            attribute("enst_refseq_compare", "1"),
            attribute("_rna_edit", "10 9 AAA"),
            attribute("rseq_cds_mismatch", ""),
            attribute("_rna_edit", "70 70 T"),
        ];
        tx.derived = Default::default();
        let alignment = tx.refseq_alignment().unwrap();
        assert_eq!(
            alignment
                .match_codes
                .iter()
                .map(|c| &**c)
                .collect::<Vec<_>>(),
            ["rseq_5p_mismatch", "rseq_cds_mismatch"]
        );
        assert_eq!(
            &*alignment.indel_edits,
            &[
                IndelEdit {
                    start: 10,
                    length_change: 3
                },
                IndelEdit {
                    start: 50,
                    length_change: -3
                },
            ]
        );
        assert_eq!(alignment.first_edit_bases.as_deref(), Some(""));
        assert_eq!(alignment.misalignment_offset(Some(30)), 3);
        assert_eq!(alignment.misalignment_offset(Some(60)), 0);
        assert_eq!(alignment.misalignment_offset(Some(5)), 0);
        assert_eq!(alignment.misalignment_offset(Some(10)), 3);
        assert_eq!(alignment.misalignment_offset(None), 0);
        assert_eq!(tx.rseq_codes().len(), 2);
        assert_eq!(tx.indel_edits().len(), 2);

        let mut insertion_first = transcript(Strand::Forward, vec![exon(100, 200, -1)], vec![]);
        insertion_first.attributes = vec![attribute("_rna_edit", "10 9 AAA")];
        let alignment = insertion_first.refseq_alignment().unwrap();
        assert!(alignment.match_codes.is_empty());
        assert_eq!(alignment.first_edit_bases.as_deref(), Some("AAA"));

        let mut codes_only = transcript(Strand::Forward, vec![exon(100, 200, -1)], vec![]);
        codes_only.attributes = vec![attribute("rseq_mrna_match", "")];
        let alignment = codes_only.refseq_alignment().unwrap();
        assert_eq!(alignment.first_edit_bases, None);
        assert!(alignment.indel_edits.is_empty());
        assert_eq!(alignment.misalignment_offset(Some(500)), 0);
    }

    /// A three-exon mapper on `strand`: cDNA 1-100 on the first exon, 101-300
    /// on the second, 301-350 on the third, the exons at 1001-1100, 2001-2200
    /// and 3001-3050 read in genomic order on the forward strand and from the
    /// top on the reverse.
    fn three_exon_mapper(strand: Strand) -> ExonCoordMapper {
        let ori = strand.as_i8();
        let genomic = [(1001, 1100), (2001, 2200), (3001, 3050)];
        let mut pairs = Vec::new();
        let mut cdna_start = 1;
        let order: Vec<(u64, u64)> = match strand {
            Strand::Forward => genomic.to_vec(),
            Strand::Reverse => genomic.iter().rev().copied().collect(),
        };
        for (to_start, to_end) in order {
            let len = to_end - to_start + 1;
            pairs.push(MapperPair {
                from_start: cdna_start,
                from_end: cdna_start + len - 1,
                to_start,
                to_end,
                ori,
            });
            cdna_start += len;
        }
        ExonCoordMapper::new(pairs)
    }

    /// The genomic base a cDNA position names, read from the pairs directly.
    fn genomic_of(mapper: &ExonCoordMapper, cdna: u64) -> u64 {
        let pair = mapper
            .pairs
            .iter()
            .find(|p| p.from_start <= cdna && cdna <= p.from_end)
            .unwrap();
        match pair.ori {
            1 => pair.to_start + (cdna - pair.from_start),
            _ => pair.to_end - (cdna - pair.from_start),
        }
    }

    /// Every cDNA base of a transcript on either strand maps to the genomic
    /// base its pair places it on and to nothing else, and a range within one
    /// exon maps to one coordinate over the same bases, so the mapping inverts
    /// the pairs the cache stores (`TranscriptMapper::cdna2genomic` over the
    /// `exon_coord_mapper` pairs `_load_mapper` adds, one per exon).
    #[test]
    fn cdna2genomic_roundtrip() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let mapper = three_exon_mapper(strand);
            for cdna in 1..=350u64 {
                let genomic = genomic_of(&mapper, cdna);
                assert_eq!(
                    mapper.cdna2genomic(cdna as i64, cdna as i64),
                    [MappedPiece::Coordinate {
                        start: genomic as i64,
                        end: genomic as i64,
                        strand
                    }],
                    "{strand:?} cDNA {cdna}"
                );
            }
            let (lo, hi) = (genomic_of(&mapper, 110), genomic_of(&mapper, 120));
            assert_eq!(
                mapper.cdna2genomic(110, 120),
                [MappedPiece::Coordinate {
                    start: lo.min(hi) as i64,
                    end: lo.max(hi) as i64,
                    strand
                }]
            );
        }
    }

    /// A range over an exon boundary is one coordinate per exon in cDNA order,
    /// so on the reverse strand the genomic coordinates descend; a range that
    /// starts before the first cDNA base opens with a gap, one that runs past
    /// the last base closes with one, and one wholly outside the pairs is a
    /// single gap, as `map_coordinates` reports them.
    #[test]
    fn cdna2genomic_lists_gaps_and_exon_pieces_in_cdna_order() {
        let forward = three_exon_mapper(Strand::Forward);
        assert_eq!(
            forward.cdna2genomic(99, 102),
            [
                MappedPiece::Coordinate {
                    start: 1099,
                    end: 1100,
                    strand: Strand::Forward
                },
                MappedPiece::Coordinate {
                    start: 2001,
                    end: 2002,
                    strand: Strand::Forward
                },
            ]
        );
        let reverse = three_exon_mapper(Strand::Reverse);
        assert_eq!(
            reverse.cdna2genomic(49, 52),
            [
                MappedPiece::Coordinate {
                    start: 3001,
                    end: 3002,
                    strand: Strand::Reverse
                },
                MappedPiece::Coordinate {
                    start: 2199,
                    end: 2200,
                    strand: Strand::Reverse
                },
            ]
        );
        assert_eq!(
            forward.cdna2genomic(-2, 3),
            [
                MappedPiece::Gap { start: -2, end: 0 },
                MappedPiece::Coordinate {
                    start: 1001,
                    end: 1003,
                    strand: Strand::Forward
                },
            ]
        );
        assert_eq!(
            forward.cdna2genomic(349, 355),
            [
                MappedPiece::Coordinate {
                    start: 3049,
                    end: 3050,
                    strand: Strand::Forward
                },
                MappedPiece::Gap {
                    start: 351,
                    end: 355
                },
            ]
        );
        assert_eq!(
            forward.cdna2genomic(400, 410),
            [MappedPiece::Gap {
                start: 400,
                end: 410
            }]
        );
        assert_eq!(
            forward.cdna2genomic(-5, -1),
            [MappedPiece::Gap { start: -5, end: -1 }]
        );
        assert_eq!(
            forward.cdna2genomic(10, 5),
            [MappedPiece::Gap { start: 10, end: 5 }]
        );
    }

    /// An insertion point (`start = end + 1`) inside an exon maps to the two
    /// genomic bases around it turned around, on either strand; one on an exon
    /// boundary maps to a piece on each side, each moved off the point
    /// (`map_insert`).
    #[test]
    fn cdna2genomic_maps_an_insertion_point_as_map_insert_does() {
        let forward = three_exon_mapper(Strand::Forward);
        assert_eq!(
            forward.cdna2genomic(11, 10),
            [MappedPiece::Coordinate {
                start: 1011,
                end: 1010,
                strand: Strand::Forward
            }]
        );
        let reverse = three_exon_mapper(Strand::Reverse);
        assert_eq!(
            reverse.cdna2genomic(11, 10),
            [MappedPiece::Coordinate {
                start: 3041,
                end: 3040,
                strand: Strand::Reverse
            }]
        );
        assert_eq!(
            forward.cdna2genomic(101, 100),
            [
                MappedPiece::Coordinate {
                    start: 1101,
                    end: 1100,
                    strand: Strand::Forward
                },
                MappedPiece::Coordinate {
                    start: 2001,
                    end: 2000,
                    strand: Strand::Forward
                },
            ]
        );
        assert_eq!(
            reverse.cdna2genomic(51, 50),
            [
                MappedPiece::Coordinate {
                    start: 3001,
                    end: 3000,
                    strand: Strand::Reverse
                },
                MappedPiece::Coordinate {
                    start: 2201,
                    end: 2200,
                    strand: Strand::Reverse
                },
            ]
        );
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
