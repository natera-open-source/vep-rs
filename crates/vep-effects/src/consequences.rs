// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Consequence calculator: assigns SO consequence terms to variant-transcript overlaps.
//!
//! Perl citations name modules of ensembl-variation release/116
//! (`Bio/EnsEMBL/Variation/...`); the `VEP::...` modules are in ensembl-vep release/116.2.
//!
//! Under `--use_transcript_ref` every row reads the reference allele
//! [`transcript_reference_allele`] gives. Not ported: the edited transcript
//! sequences of a BAM-edited cache, which switch the flag on by themselves
//! (`VEP::CacheDir`, the `bam` entry of the cache's info file) and supply the
//! cDNA the reference allele is cut from (`VEP::AnnotationType::Transcript`,
//! `apply_edits`); the `--bam` flag that builds such a cache is refused, and
//! the `invalid_alleles` warning of an edited RefSeq transcript is not
//! written.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::sync::Arc;

use smallvec::{smallvec, SmallVec};
use vep_core::consequence::{Consequence, ConsequenceList, Impact, TranscriptConsequence};
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

use crate::coding::{
    classify_inframe_indel_by_peptide, compute_codon_window_peptide_alleles,
    compute_peptide_alleles, get_codon_change, InframeIndelKind,
};
use crate::coding::{frameshift_stop_gained_full_cds, frameshift_stop_gained_in_codon_window};
use crate::mapper::{map_genomic_to_transcript, CdsSpanBounds, TranscriptPosition};

/// Configuration for consequence calculation.
pub struct EffectsConfig {
    /// Maximum upstream distance in base pairs (default 5000).
    pub upstream_distance: u64,
    /// Maximum downstream distance in base pairs (default 5000).
    pub downstream_distance: u64,
    /// Optional reference FASTA for sequence-aware behaviors (e.g. 3' shifting indels).
    pub reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    /// `--shift_3prime 1`: move each insertion or deletion to its most 3' position
    /// on the transcript, within [`perform_shift_bound`], before its terms,
    /// positions, codons and amino acids are computed; needs `reference_fasta`.
    ///
    /// Off by default because it changes consequence semantics.
    pub enable_indel_3prime_shift: bool,
    /// The flags that shape the HGVSc and HGVSp notations: `--shift_hgvs`,
    /// `--hgvsp_use_prediction`, `--remove_hgvsp_version` and, for the
    /// predicates `hgvs_protein` reads, `--shift_3prime`, which an embedder
    /// sets as it sets `enable_indel_3prime_shift`.
    pub hgvs: crate::hgvs::HgvsOptions,
    /// Attach [`LofteeContext`](vep_core::consequence::LofteeContext) to each
    /// transcript consequence at construction time.
    ///
    /// Off by default. Set only when a consequence-filtering plugin (LoFTEE)
    /// is active, so the exon/intron model and genomic CDS bounds travel to the
    /// plugin phase. When off, `TranscriptConsequence.loftee_ctx` stays `None`
    /// and output is byte-identical to a non-plugin run.
    pub populate_loftee_context: bool,
    /// Compute HGVSc / HGVSp notations for each transcript consequence.
    ///
    /// Off by default. On the default `-o vep` path HGVS is emitted only under
    /// `--hgvs` (Ensembl VEP parity: Perl emits HGVS solely under that flag), so
    /// computing it unconditionally spends significant CPU building `String`s
    /// that are then discarded. When off, `TranscriptConsequence.hgvsc`,
    /// `.hgvsp` and `.hgvs_offset` stay `None`.
    ///
    /// Callers enable this when the consumer prints the fields: `vep-cli` turns
    /// it on under `--hgvs` with a reference FASTA present, the condition under
    /// which Ensembl VEP fills `HGVSc`/`HGVSp` in every output format. Embedders
    /// that always populate HGVS pass `true`.
    pub compute_hgvs: bool,
    /// Compute the `EXON` and `INTRON` ordinals for each transcript consequence.
    ///
    /// Off by default. Every output format prints them only under `--numbers`
    /// (which `--vcf` and `--everything` imply), so callers turn this on
    /// whenever that flag is set, and whenever a plugin runs: LoFTEE reads
    /// `TranscriptConsequence.exon` / `.intron` directly and a dylib plugin
    /// receives the whole consequence as JSON. When off, both fields stay
    /// `None`.
    pub compute_exon_intron_numbers: bool,
    /// Append the translation's version to the protein identifier
    /// (`ENSP00000286808.3`), as `--protein_version` does.
    pub protein_version: bool,
    /// Compute `InputVariant::hgvsg` for every sequence variant, as `--hgvsg`
    /// does; needs `reference_fasta`. `None` leaves the field `None`.
    pub hgvsg: Option<crate::hgvsg::HgvsgOptions>,
    /// Compute `InputVariant::spdi`, and the GA4GH VRS allele the options ask
    /// for, for every sequence variant, as `--spdi` and `--ga4gh_vrs` do;
    /// `reference_fasta` names the chromosome of a variant overlapping no
    /// transcript when it is set. `None` leaves both fields `None`.
    pub spdi: Option<crate::spdi::SpdiOptions>,
    /// List the protein features the variant's residues overlap on each coding
    /// consequence, as `--domains` does; off, `TranscriptConsequence.domains`
    /// stays empty.
    pub domains: bool,
    /// `--use_transcript_ref`: each row is analysed against the reference allele
    /// [`transcript_reference_allele`] reads from the transcript, and carries it
    /// as `TranscriptConsequence.used_ref`; a non-coding transcript's cDNA is
    /// read from `reference_fasta`.
    pub use_transcript_ref: bool,
    /// `--process_ref_homs`: the record's reference row is written, and its
    /// codon is read before any alternate allele's, which keeps the codon path
    /// from rotating `used_ref` with a 3' shift (see `used_ref_rotation`).
    pub process_ref_homs: bool,
}

impl Default for EffectsConfig {
    fn default() -> Self {
        Self {
            upstream_distance: 5000,
            downstream_distance: 5000,
            reference_fasta: None,
            enable_indel_3prime_shift: false,
            hgvs: crate::hgvs::HgvsOptions::default(),
            populate_loftee_context: false,
            compute_hgvs: false,
            compute_exon_intron_numbers: false,
            protein_version: false,
            hgvsg: None,
            spdi: None,
            domains: false,
            use_transcript_ref: false,
            process_ref_homs: false,
        }
    }
}

/// Classify the stop-codon effect from full-CDS peptide alleles.
///
/// Returns `StopRetainedVariant` when both ref and alt contain a stop,
/// `StopLost` when only ref contains a stop, or `None` when the ref
/// has no stop (unresolvable).
fn classify_stop_from_peptides(ref_pep: &[u8], alt_pep: &[u8]) -> Option<Consequence> {
    if ref_pep.contains(&b'*') && alt_pep.contains(&b'*') {
        Some(Consequence::StopRetainedVariant)
    } else if ref_pep.contains(&b'*') && !alt_pep.contains(&b'*') {
        Some(Consequence::StopLost)
    } else {
        None
    }
}

/// The one output field the term analysis itself supplies; the display columns
/// come from [`crate::display`].
#[derive(Default)]
struct PositionFields {
    distance: Option<u64>,
}

fn insertion_analysis_cds_pos(cds_pos: u64, coding_bounds: Option<&CdsSpanBounds>) -> u64 {
    coding_bounds
        .filter(|bounds| bounds.cds_start > bounds.cds_end)
        .map(|bounds| bounds.cds_start)
        .unwrap_or(cds_pos)
}

/// Reverse-strand-aware CDS analysis position for the start/stop-codon peptide
/// blocks.
///
/// `*cds_pos` is the per-endpoint CDS index from `apply_position`'s start
/// genomic anchor (`populate_fields=true`). On a reverse-strand deletion that
/// anchor maps to the higher (3'-most) CDS index of the span, e.g. `cds_pos=4`
/// for a span covering CDS 2-4, or `cds_pos=1416` for a span covering CDS
/// 1415-1416. Feeding that anchor to `compute_peptide_alleles` either excises
/// the wrong CDS bases or returns `None` (`idx + ref_len > cds.len()`), so the
/// start_lost / stop_lost peptide path silently undercalls on reverse-strand
/// deletions.
///
/// Perl's `start_lost` / `stop_lost` predicates work in cDNA/CDS space, which is
/// strand-normalized, so they always anchor at the 5'-most CDS index of the
/// span. Here, for deletions/substitutions (`cds_start <= cds_end`)
/// the 5'-most index is `cds_start.min(cds_end)`; for insertions
/// (`cds_start > cds_end`) the insertion anchor (`cds_start`, the
/// higher value) is preserved. When bounds are unavailable, fall back to the
/// raw per-endpoint `cds_pos` (forward-strand behavior, where the start anchor
/// already is the 5'-most index).
fn analysis_cds_pos_5prime(cds_pos: u64, coding_bounds: Option<&CdsSpanBounds>) -> u64 {
    match coding_bounds {
        Some(bounds) if bounds.cds_start > bounds.cds_end => bounds.cds_start, // insertion anchor
        Some(bounds) => bounds.cds_start.min(bounds.cds_end), // deletion/sub: 5'-most CDS index
        None => cds_pos,
    }
}

fn should_use_insertion_anchor(raw_insertion: bool, coding_bounds: Option<&CdsSpanBounds>) -> bool {
    raw_insertion || coding_bounds.is_some_and(|bounds| bounds.cds_start > bounds.cds_end)
}

fn push_unique(consequences: &mut ConsequenceList, consequence: Consequence) {
    if !consequences.contains(&consequence) {
        consequences.push(consequence);
    }
}

fn remove_consequence(consequences: &mut ConsequenceList, consequence: Consequence) {
    consequences.retain(|c| *c != consequence);
}

/// Calculate all consequences for a variant against a transcript.
///
/// Maps the variant's start position to transcript coordinates and assigns
/// the appropriate SO consequence terms based on where the variant falls
/// relative to the transcript structure (UTR, coding, intron, etc.).
pub fn calculate_consequences(
    variant: &InputVariant,
    transcript: &Transcript,
    config: &EffectsConfig,
) -> Option<TranscriptConsequence> {
    // Under `--shift_3prime 1` the allele is analysed where the shift leaves it:
    // every predicate reads the moved span (`get_all_OverlapConsequences` moves
    // the feature's start and end for each one, BaseVariationFeatureOverlapAllele.pm:269-282,
    // and `_bvfo_preds` adds the offset, :454-463) and the rotated allele
    // (`shift_feature_seqs`, TranscriptVariationAllele.pm:1330, which rotates the
    // allele as the input wrote it, on the input's strand), as do the cDNA,
    // CDS and translation coordinates (`clear_shifting_variables`, :502-528) and
    // the codons (:808-858). The input span keeps DISTANCE
    // (BaseTranscriptVariation.pm:587-601) and the HGVS notations, which shift on
    // their own; the exon and intron numbers read both spans. An allele given on
    // the reverse strand is analysed on the forward strand, as `feature_seq`
    // reads it against each transcript; the HGVS notations take it as given,
    // because `hgvs_transcript` compares the input's strand with the
    // transcript's (TranscriptVariationAllele.pm:1463-1466).
    let placed: Cow<InputVariant> = match variant.strand {
        vep_core::coordinate::Strand::Forward => Cow::Borrowed(variant),
        vep_core::coordinate::Strand::Reverse => Cow::Owned(variant.on_forward_strand()),
    };
    let placed = placed.as_ref();
    // Under `--use_transcript_ref` the reference allele's object is built from
    // the transcript's sequence (`-use_feature_ref`, AnnotationType/Transcript.pm:161;
    // VariationFeatureOverlap.pm:418-425), so every predicate that compares the
    // alleles reads that allele, the reference row is that allele, and the HGVS
    // shift moves the deleted bases it names (`_genomic_shift` reads the
    // reference allele's object, TranscriptVariationAllele.pm:439-449). `USED_REF`
    // prints it relative to the variant (`variation_feature_seq`,
    // OutputFactory.pm:1776), so an allele read on the reverse strand is given
    // back on it.
    let used_ref = if config.use_transcript_ref
        && !placed.is_structural
        && !placed.variant_class.is_structural()
    {
        Some(
            transcript_reference_allele(placed, transcript, config.reference_fasta.as_deref())
                .unwrap_or_else(|| placed.ref_allele.clone()),
        )
    } else {
        None
    };
    let on_input_strand = |reference: &[u8]| match variant.strand {
        vep_core::coordinate::Strand::Reverse if vep_core::codon::seq_is_dna(reference) => {
            vep_core::codon::reverse_complement(reference)
        }
        _ => reference.to_vec(),
    };
    let (placed, input): (Cow<InputVariant>, Cow<InputVariant>) = match &used_ref {
        Some(reference) if *reference != placed.ref_allele => (
            Cow::Owned(placed.against_reference_allele(reference)),
            Cow::Owned(variant.against_reference_allele(&on_input_strand(reference))),
        ),
        _ => (Cow::Borrowed(placed), Cow::Borrowed(variant)),
    };
    let (placed, input) = (placed.as_ref(), input.as_ref());
    let shifted = consequence_shift(placed, variant.strand, transcript, config);
    let shifted_allele = shifted.and_then(|span| {
        if placed.reference_allele_row {
            Some(shifted_reference_row(
                placed,
                span,
                config.reference_fasta.as_deref(),
            ))
        } else {
            crate::hgvs::shifted_variant(
                placed,
                variant.strand,
                transcript,
                span,
                config.reference_fasta.as_deref(),
            )
        }
    });
    let analysed = shifted_allele.as_ref().unwrap_or(placed);
    let mut tc = calculate_terms(analysed, input, transcript, config)?;
    let mut fields = crate::display::display_fields(analysed, transcript);
    let translation_span = fields.translation_span;
    if config.domains {
        tc.domains = domains_overlapping_protein_position(transcript, translation_span);
    }
    fields.record_span_positions(variant, transcript);
    fields.write_to(&mut tc);
    if config.compute_exon_intron_numbers {
        (tc.exon, tc.intron) = crate::display::exon_intron_numbers(analysed, placed, transcript);
    }
    // A reference row prints `SHIFT_LENGTH` 0 whatever span it was read at: the
    // reference allele's object carries no shift of its own (OutputFactory.pm:1690).
    tc.shift_length = match shifted {
        Some((start, _)) if !placed.reference_allele_row => start.abs_diff(placed.start),
        _ => 0,
    };
    tc.shifted_allele = shifted_display_allele(
        variant,
        transcript.strand,
        translation_span,
        tc.shift_length,
    );
    if config.protein_version {
        append_translation_version(&mut tc, transcript);
    }
    tc.mirna = mirna_structure_terms(analysed, transcript);
    // `hgvs_protein` takes the shifted bases on the forward strand and gives
    // them back on the input's (`create_shift_hash`, TranscriptVariationAllele.pm:373);
    // `shift_feature_seqs` rotates the allele as the input wrote it (:1334-1343).
    tc.used_ref = used_ref.map(|mut reference| {
        let as_given = match used_ref_rotation(&tc, translation_span, config.process_ref_homs) {
            Some(rotation @ UsedRefRotation::Protein(_)) => {
                rotate_along_transcript(&mut reference, rotation, transcript.strand);
                on_input_strand(&reference)
            }
            Some(rotation @ UsedRefRotation::Codon(_)) => {
                let mut as_given = on_input_strand(&reference);
                rotate_along_transcript(&mut as_given, rotation, transcript.strand);
                as_given
            }
            None => on_input_strand(&reference),
        };
        Arc::from(String::from_utf8_lossy(&as_given).as_ref())
    });
    Some(tc)
}

/// The path that rotates a shifted deletion's `USED_REF`, with the shift in
/// bases. Ensembl VEP rotates the reference allele's sequence along with a
/// shifted alternate allele only through the codon and protein paths, so only
/// where the shifted span has both translation coordinates: `Protein` when
/// `hgvs_protein` writes the deleted bases `perform_shift` moved to into the
/// reference allele (TranscriptVariationAllele.pm:1771-1786) after the
/// translation test at the shifted offset (:1730-1737), which a protein
/// notation attests; otherwise `Codon` under `--shift_3prime 1`, where
/// `display_codon_allele_string` hands the alternate allele's shift to the
/// reference allele (:674) and `codon` rotates its sequence
/// (`shift_feature_seqs`, :1330-1345) once `translation_start` and
/// `translation_end` are defined (:813-820). A non-coding row, a row in a UTR
/// and a deletion whose shifted span leaves the coding sequence keep the bases
/// as the input wrote them. `codon` rotates only on the call that memoises the
/// reference allele's codon (:802-804), and under `--process_ref_homs` the
/// reference row reads that codon, with no shift, before any alternate
/// allele's row hands it one, so the codon path then rotates nothing
/// (`reference_row_first`); the protein path writes the sequence outright.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum UsedRefRotation {
    Codon(u64),
    Protein(u64),
}

fn used_ref_rotation(
    tc: &TranscriptConsequence,
    translation_span: Option<(i64, i64)>,
    reference_row_first: bool,
) -> Option<UsedRefRotation> {
    let protein_shift = match (&tc.hgvsp, tc.hgvs_offset) {
        (Some(_), Some(offset)) if offset != 0 => Some(offset.unsigned_abs()),
        _ => None,
    };
    if let Some(shift) = protein_shift {
        return Some(UsedRefRotation::Protein(shift));
    }
    if tc.shift_length > 0 && translation_span.is_some() && !reference_row_first {
        return Some(UsedRefRotation::Codon(tc.shift_length));
    }
    None
}

/// Rotates `allele` as the path moves it. `perform_shift`
/// (TranscriptVariationAllele.pm:291-350) moves the deleted bases one place per
/// step towards the transcript's 3' end, the base that falls off the end coming
/// back at the start, so a `Protein` rotation wraps on either strand;
/// `shift_feature_seqs` (:1330-1345) moves the bases `shift` places on a
/// forward-strand transcript and `length - shift` places on a reverse-strand
/// one (:1338-1339), so a `Codon` rotation on the reverse strand moves nothing
/// when the shift exceeds the allele's length.
fn rotate_along_transcript(
    allele: &mut [u8],
    rotation: UsedRefRotation,
    strand: vep_core::coordinate::Strand,
) {
    let length = allele.len() as u64;
    if length < 2 {
        return;
    }
    let shift = match rotation {
        UsedRefRotation::Protein(shift) => shift,
        UsedRefRotation::Codon(shift) => {
            if strand == vep_core::coordinate::Strand::Reverse && shift > length {
                return;
            }
            shift
        }
    };
    let steps = (shift % length) as usize;
    match strand {
        vep_core::coordinate::Strand::Forward => allele.rotate_left(steps),
        vep_core::coordinate::Strand::Reverse => allele.rotate_right(steps),
    }
}

/// [`TranscriptConsequence::shifted_allele`]: the inserted bases as the record
/// wrote them, rotated as `shift_feature_seqs` rotates them when `codon` runs
/// on a shifted allele (TranscriptVariationAllele.pm:820, 1330-1347), which it
/// does on a row whose shifted span has both translation coordinates: one
/// place per base of the shift along a forward-strand transcript, the allele's
/// length less the shift along a reverse-strand one (:1338-1339), where a shift
/// of the allele's length or more leaves it as written. `None` wherever the
/// `Allele` is the allele as given: a row nothing moved, a row with no
/// translation coordinate, a reference row, any allele but an insertion of DNA,
/// and a rotation by a multiple of the allele's length.
fn shifted_display_allele(
    variant: &InputVariant,
    transcript_strand: vep_core::coordinate::Strand,
    translation_span: Option<(i64, i64)>,
    shift: u64,
) -> Option<Arc<str>> {
    if shift == 0
        || translation_span.is_none()
        || variant.reference_allele_row
        || variant.variant_class != vep_core::variant::VariantClass::Insertion
    {
        return None;
    }
    let mut written = variant.alt_allele().to_ascii_uppercase();
    if written == b"-" || !vep_core::codon::seq_is_dna(&written) {
        return None;
    }
    let length = written.len() as u64;
    let steps = match transcript_strand {
        vep_core::coordinate::Strand::Forward => shift,
        vep_core::coordinate::Strand::Reverse => length.checked_sub(shift)?,
    } % length;
    if steps == 0 {
        return None;
    }
    written.rotate_left(steps as usize);
    Some(Arc::from(String::from_utf8_lossy(&written).as_ref()))
}

/// The reference allele a row reads under `--use_transcript_ref`
/// (`TranscriptVariation::_get_ref_allele`, TranscriptVariation.pm:441-468 at
/// release/116, reached with `-no_ref_check` from AnnotationType/Transcript.pm:159):
/// `-` for an insertion; the transcript's spliced cDNA over the variant's cDNA
/// span when the span maps onto one stretch of it (at most two mapped pieces,
/// the first and last both exonic), on the forward strand, an empty stretch
/// read as `-` (VariationFeatureOverlap.pm:425); and `None` for any other
/// span, because the genomic slice branch is skipped under `-no_ref_check`
/// (`:368`) and `_raw_allele_hashes` then keeps the input's reference allele
/// (`:418-424`). A coding transcript's cDNA is the cache's
/// ([`Transcript::spliced_seq`]); any other transcript's is read from `fasta`
/// over the exons the span covers, as Ensembl's `Transcript::spliced_seq`
/// reads every exon, and is `None` without one. `variant` is the allele on the
/// forward strand.
pub fn transcript_reference_allele(
    variant: &InputVariant,
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<Vec<u8>> {
    use crate::display::Segment;
    if variant.start > variant.end {
        return Some(b"-".to_vec());
    }
    if variant.end < transcript.start || variant.start > transcript.end {
        return None;
    }
    let mapper = transcript.vefc.as_ref()?.mapper.as_ref()?;
    let segments = crate::display::genomic2cdna(
        mapper,
        variant.start as i64,
        variant.end as i64,
        transcript.strand.as_i8(),
    );
    if segments.len() > 2 {
        return None;
    }
    let (Some(Segment::Coord { start, .. }), Some(Segment::Coord { end, .. })) =
        (segments.first(), segments.last())
    else {
        return None;
    };
    let (cdna_start, cdna_end) = (u64::try_from(*start).ok()?, u64::try_from(*end).ok()?);
    let mut reference = match transcript.spliced_seq() {
        // Perl's `substr` reads to the end of a shorter sequence.
        Some(cdna) => {
            let from = (cdna_start.saturating_sub(1) as usize).min(cdna.len());
            let to = (cdna_end as usize).clamp(from, cdna.len());
            let bases = cdna.as_bytes()[from..to].to_vec();
            match transcript.strand {
                vep_core::coordinate::Strand::Reverse => {
                    vep_core::codon::reverse_complement(&bases)
                }
                vep_core::coordinate::Strand::Forward => bases,
            }
        }
        None => cdna_bases_from_fasta(
            fasta?,
            &transcript.chr,
            &mapper.exon_coord_mapper,
            cdna_start,
            cdna_end,
        )?,
    };
    if reference.is_empty() {
        reference = b"-".to_vec();
    }
    Some(reference)
}

/// The bases of the cDNA span `cdna_start..=cdna_end` read from the FASTA on
/// `chr` over the exon pieces `mapper` projects it to, on the forward strand;
/// `None` when a piece cannot be read. The buffer grows with the bases the
/// FASTA returns, which its records bound, never with the span's cDNA ends,
/// which a mapper from a malformed cache can set any distance apart.
fn cdna_bases_from_fasta(
    fasta: &vep_fasta::IndexedFasta,
    chr: &str,
    mapper: &vep_core::transcript::ExonCoordMapper,
    cdna_start: u64,
    cdna_end: u64,
) -> Option<Vec<u8>> {
    use vep_core::transcript::MappedPiece;
    let mut pieces: Vec<(u64, u64)> = mapper
        .cdna2genomic(cdna_start as i64, cdna_end as i64)
        .into_iter()
        .filter_map(|piece| match piece {
            MappedPiece::Coordinate { start, end, .. } => Some((start as u64, end as u64)),
            MappedPiece::Gap { .. } => None,
        })
        .collect();
    // Genomic order is forward-strand order whichever strand the transcript is on.
    pieces.sort_unstable();
    let mut bases = Vec::new();
    for (genomic_start, genomic_end) in pieces {
        bases.extend(fasta.sequence(chr, genomic_start, genomic_end)?);
    }
    Some(bases)
}

/// The `miRNA` column's terms (OutputFactory.pm:1624-1664 at release/116.2).
///
/// On a transcript of biotype `miRNA` whose first `ncRNA` attribute reads
/// `<start>:<end> <structure>` (three fields split on whitespace or `:`), the
/// dot-bracket structure is expanded (`.7(16` is seven `.` then sixteen `(`)
/// and every cDNA position from the variant's cDNA start to its end (an
/// insertion's two flanking bases, put in order) that falls inside the
/// structure contributes its element's term: `miRNA_stem` for a bracket,
/// `miRNA_loop` for a dot. The terms come sorted, one per distinct structure
/// character, so a span over both an opening and a closing bracket lists
/// `miRNA_stem` twice. Nothing is produced when the attribute lacks a
/// structure, when either cDNA end is undefined (the variant starts or ends in
/// an intron or outside the transcript) or when the variant lies outside
/// `<start>..<end>`.
///
/// Ensembl VEP's index test reads one position past the structure's end as an
/// undefined element, which sorts first; it is kept as an empty term, which
/// beside another term the text formats print as an empty list element and
/// the JSON format as `null`. A list of that term alone joins to an empty
/// value, which the writers drop as they drop every empty value, where
/// Ensembl VEP prints an empty `miRNA=`, an empty tab cell and `[null]`; no
/// Ensembl structure is shorter than its span, so the case is not reached on
/// its caches.
fn mirna_structure_terms(variant: &InputVariant, transcript: &Transcript) -> Vec<String> {
    if !transcript.facts().is_mirna {
        return Vec::new();
    }
    let Some(attribute) = transcript.attributes.iter().find(|a| a.code == "ncRNA") else {
        return Vec::new();
    };
    let fields = ncrna_fields(&attribute.value);
    let (Some(start), Some(end), Some(structure)) = (
        fields.first().and_then(|f| f.parse::<i64>().ok()),
        fields.get(1).and_then(|f| f.parse::<i64>().ok()),
        fields.get(2).copied(),
    ) else {
        return Vec::new();
    };
    if start == 0 || end == 0 || !structure.contains(['(', '.', ')']) {
        return Vec::new();
    }
    let Some((cdna_start, cdna_end)) = cdna_span(variant, transcript) else {
        return Vec::new();
    };
    let overlaps_structure = end >= cdna_start && start <= cdna_end;
    if cdna_start == 0 || cdna_end == 0 || !overlaps_structure {
        return Vec::new();
    }
    let (cdna_start, cdna_end) = (cdna_start.min(cdna_end), cdna_start.max(cdna_end));

    let mut elements: Vec<char> = Vec::with_capacity(structure.len());
    let mut chars = structure.chars().peekable();
    while let Some(c) = chars.next() {
        if !matches!(c, '(' | '.' | ')') {
            continue;
        }
        let mut digits = String::new();
        while let Some(d) = chars.next_if(|d| d.is_ascii_digit()) {
            digits.push(d);
        }
        let count = digits.parse::<usize>().unwrap_or(1).max(1);
        elements.extend(std::iter::repeat_n(c, count));
    }

    // `sort map {$map{$_}} keys %chars` (OutputFactory.pm:1662): the set is of
    // structure characters, not of terms, so `(` and `)` each yield a `miRNA_stem`.
    let mut overlapped: BTreeSet<Option<char>> = BTreeSet::new();
    for position in cdna_start..=cdna_end {
        let offset = position - start;
        if offset < 0 || offset > elements.len() as i64 {
            continue;
        }
        overlapped.insert(elements.get(offset as usize).copied());
    }
    let mut terms: Vec<String> = overlapped
        .into_iter()
        .map(|element| match element {
            Some('.') => "miRNA_loop".to_string(),
            Some(_) => "miRNA_stem".to_string(),
            None => String::new(),
        })
        .collect();
    terms.sort();
    terms
}

/// The first three fields of an `ncRNA` attribute value, split as Perl's
/// `split /\s+|\:/` splits it: a colon or a run of whitespace ends a field, so
/// a colon beside whitespace leaves an empty field between them.
fn ncrna_fields(value: &str) -> Vec<&str> {
    let mut fields = Vec::with_capacity(3);
    let mut field_start = 0;
    let bytes = value.as_bytes();
    let mut i = 0;
    while i < bytes.len() && fields.len() < 3 {
        let separator_len = if bytes[i] == b':' {
            1
        } else if bytes[i].is_ascii_whitespace() {
            bytes[i..]
                .iter()
                .take_while(|b| b.is_ascii_whitespace())
                .count()
        } else {
            i += 1;
            continue;
        };
        fields.push(&value[field_start..i]);
        i += separator_len;
        field_start = i;
    }
    if fields.len() < 3 && field_start < value.len() {
        fields.push(&value[field_start..]);
    }
    fields
}

/// The variant's cDNA start and end on the transcript as `BaseTranscriptVariation`'s
/// `cdna_start` and `cdna_end` give them: the first mapped segment's start and
/// the last's end, `None` when either segment is a gap. An insertion maps its
/// two flanking bases, so its start exceeds its end.
fn cdna_span(variant: &InputVariant, transcript: &Transcript) -> Option<(i64, i64)> {
    use crate::display::Segment;
    let mapper = transcript.vefc.as_ref()?.mapper.as_ref()?;
    let segments = crate::display::genomic2cdna(
        mapper,
        variant.start as i64,
        variant.end as i64,
        transcript.strand.as_i8(),
    );
    let Segment::Coord { start, .. } = segments.first()? else {
        return None;
    };
    let Segment::Coord { end, .. } = segments.last()? else {
        return None;
    };
    Some((*start, *end))
}

/// The span an insertion or deletion is annotated at under `--shift_3prime 1`:
/// `TranscriptVariation::new` shifts every allele unless `-no_shift`
/// (TranscriptVariation.pm:132-135), and `_return_3prime` moves nothing for a
/// structural or symbolic allele (TranscriptVariationAllele.pm:123, 129). The
/// reference allele is never shifted itself (:116), but the reference row of
/// `--process_ref_homs` is read at the span the record's last alternate allele
/// shifts to ([`InputVariant::reference_row_shift_allele`]), the coordinates
/// `reset_shifted_positions` leaves on the transcript variation the alleles
/// share (OutputFactory.pm:306-323, TranscriptVariationAllele.pm:502-528). That
/// alternate is kept as the record wrote it, on `input_strand`, while `variant`
/// is read on the forward strand, so an alternate given on the reverse strand
/// is complemented before it is shifted. `None` when the run does not shift, no
/// FASTA is loaded, or the allele stays put.
fn consequence_shift(
    variant: &InputVariant,
    input_strand: vep_core::coordinate::Strand,
    transcript: &Transcript,
    config: &EffectsConfig,
) -> Option<(u64, u64)> {
    if !config.enable_indel_3prime_shift
        || variant.is_structural
        || variant.variant_class.is_structural()
    {
        return None;
    }
    let fasta = config.reference_fasta.as_deref()?;
    if variant.reference_allele_row {
        let alternate = variant.reference_row_shift_allele.as_deref()?.as_bytes();
        let alternate = match input_strand {
            vep_core::coordinate::Strand::Reverse if vep_core::codon::seq_is_dna(alternate) => {
                vep_core::codon::reverse_complement(alternate)
            }
            _ => alternate.to_vec(),
        };
        let mut shifting = InputVariant::new(
            variant.chr.clone(),
            variant.start,
            variant.end,
            variant.ref_allele.clone(),
            alternate,
        );
        shifting.strand = variant.strand;
        shifting.original_start = variant.original_start;
        shifting
            .record_allele_string_multi
            .clone_from(&variant.record_allele_string_multi);
        return shift_indel_3prime_coords(&shifting, transcript, fasta);
    }
    shift_indel_3prime_coords(variant, transcript, fasta)
}

/// The reference row of `--process_ref_homs` at `span`, the shifted span of the
/// record's alternate allele: the row's alleles are the reference bases there,
/// which the transcript's translateable sequence spells for its codon
/// (TranscriptVariationAllele.pm:868, 877); an insertion's reference row keeps
/// its empty allele. Without a FASTA the bases as given stand.
fn shifted_reference_row(
    placed: &InputVariant,
    span: (u64, u64),
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> InputVariant {
    let (start, end) = span;
    let mut row = placed.clone();
    (row.start, row.end) = (start, end);
    if placed.variant_class != vep_core::variant::VariantClass::Insertion && start <= end {
        if let Some(bases) = fasta.and_then(|f| f.sequence(&placed.chr, start, end)) {
            let bases = bases.to_ascii_uppercase();
            row.ref_allele.clone_from(&bases);
            row.alt_alleles = vec![bases];
        }
    }
    row
}

/// `--protein_version` (OutputFactory.pm `BaseTranscriptVariationAllele_to_output_hash`):
/// the protein identifier gains `.<translation version>` unless the translation
/// has no version or the identifier already ends in `.<digits>`.
fn append_translation_version(tc: &mut TranscriptConsequence, transcript: &Transcript) {
    let Some(version) = transcript.translation.as_ref().and_then(|t| t.version) else {
        return;
    };
    let Some(id) = tc.protein_id.as_deref() else {
        return;
    };
    let versioned = id
        .rsplit_once('.')
        .is_some_and(|(_, tail)| !tail.is_empty() && tail.bytes().all(|b| b.is_ascii_digit()));
    if !versioned {
        tc.protein_id = Some(Arc::from(format!("{id}.{version}")));
    }
}

/// The `DOMAINS` pairs of a coding consequence (`OutputFactory.pm:1494-1512`
/// at release/116.2; `BaseTranscriptVariation.pm` `get_overlapping_ProteinFeatures`,
/// release/116): every protein feature of the transcript whose residue span
/// overlaps `translation_span`, in the cache's order, each as its analysis
/// label and its identifier with every whitespace character, `;` and `=`
/// replaced by `_`. Nothing when either residue end is undefined, as for a
/// variant that leaves the coding sequence. The overlap is Ensembl's
/// `feature.end >= start && feature.start <= end` on the raw span, so an
/// insertion between two residues (start one past end) lists only the features
/// covering both.
fn domains_overlapping_protein_position(
    transcript: &Transcript,
    translation_span: Option<(i64, i64)>,
) -> Vec<(String, String)> {
    let Some((start, end)) = translation_span else {
        return Vec::new();
    };
    let Some(vefc) = transcript.vefc.as_ref() else {
        return Vec::new();
    };
    let printable = |text: &str| -> String {
        // Perl's `\s` matches the vertical tab, which `is_ascii_whitespace` does not.
        text.chars()
            .map(|c| match c {
                ' ' | '\t' | '\n' | '\r' | '\u{0b}' | '\u{0c}' | ';' | '=' => '_',
                other => other,
            })
            .collect()
    };
    vefc.protein_features
        .iter()
        .filter(|feature| feature.end as i64 >= start && feature.start as i64 <= end)
        .map(|feature| {
            (
                printable(feature.analysis.as_deref().unwrap_or_default()),
                printable(&feature.hseqname),
            )
        })
        .collect()
}

/// The terms Ensembl VEP's `_intron_effects` decides (VariationEffect.pm
/// `within_intron`, the splice-site predicates and `splice_region`).
const INTRON_EFFECT_TERMS: [Consequence; 7] = [
    Consequence::IntronVariant,
    Consequence::SpliceAcceptorVariant,
    Consequence::SpliceDonorVariant,
    Consequence::SpliceDonor5thBaseVariant,
    Consequence::SpliceDonorRegionVariant,
    Consequence::SplicePolypyrimidineTractVariant,
    Consequence::SpliceRegionVariant,
];

/// Consequence terms and the fields the analysis derives alongside them.
/// `variant` is the allele as analysed, `input` the allele as the input gave
/// it, span and strand; the spans differ only under `--shift_3prime 1`. The
/// transcript window, DISTANCE and the HGVS notations read `input`.
fn calculate_terms(
    variant: &InputVariant,
    input: &InputVariant,
    transcript: &Transcript,
    config: &EffectsConfig,
) -> Option<TranscriptConsequence> {
    if variant.is_structural || variant.variant_class.is_structural() {
        return crate::sv::calculate_sv_consequences(
            variant,
            transcript,
            config.upstream_distance,
            config.downstream_distance,
        );
    }

    let insertion = matches!(
        variant.variant_class,
        vep_core::variant::VariantClass::Insertion
    ) && variant.start == variant.end + 1;

    // Perl keeps a transcript for a variant only when its non-normalising
    // `overlap(vf.start, vf.end, tr.start - dist, tr.end + dist)` holds
    // (`VEP::AnnotationType::Transcript::annotate_InputBuffer`,
    // `VEP::InputBuffer::get_overlapping_vfs`). The left-hand test reads
    // `vf.end >= tr.start - dist`, so an insertion (start = end + 1) anchored
    // exactly `dist` bases before the transcript is one base outside the window;
    // the right-hand test reads `vf.start <= tr.end + dist`, which the anchor
    // mapping below already applies.
    if insertion {
        let left_dist = match transcript.strand {
            vep_core::coordinate::Strand::Forward => config.upstream_distance,
            vep_core::coordinate::Strand::Reverse => config.downstream_distance,
        };
        if input.end < transcript.start.saturating_sub(left_dist) {
            return None;
        }
    }

    // Perl's `complete_overlap` pre-predicate
    // (`BaseVariationFeatureOverlapAllele::_bvf_preds`): `vf.start <= tr.start
    // && vf.end >= tr.end` on the raw coordinates, which an insertion
    // (start = end + 1) cannot satisfy.
    let complete_overlap =
        !insertion && variant.start <= transcript.start && variant.end >= transcript.end;

    // `transcript_ablation` is a candidate only under `include => {complete_overlap
    // => 1, deletion => 1}` (Utils/Constants.pm), and it is the sole term because
    // it is tier 1 (`get_all_OverlapConsequences` stops before tier 3 once a tier
    // 1 or 2 term fires). A same-length substitution or a longer ALT covering the
    // transcript is `snp` or `insertion` there, never `deletion`, so it takes the
    // tier-3 predicates below.
    if complete_overlap && perl_deletion_pre_predicate(variant) {
        return Some(TranscriptConsequence {
            consequences: smallvec![Consequence::TranscriptAblation],
            impact: Impact::HIGH,
            ..TranscriptConsequence::for_transcript(transcript)
        });
    }

    // The window test reads the input span, as the fetch does; a 3' shift moves
    // the allele towards the transcript's downstream edge, so it can carry the
    // analysed span up to `shift` bases past that edge, where its predicates
    // still read `downstream`: the analysed span is mapped with the downstream
    // window widened by that much.
    let shift = variant.start.abs_diff(input.start);
    if shift > 0 && !complete_overlap {
        let within_window = |pos| {
            map_genomic_to_transcript(
                pos,
                transcript,
                config.upstream_distance,
                config.downstream_distance,
            )
            .is_some()
        };
        let input_within_window = within_window(input.start) || within_window(input.end);
        if !input_within_window {
            return None;
        }
    }
    let upstream_window = config.upstream_distance;
    let downstream_window = config.downstream_distance.saturating_add(shift);
    let start_pos = map_genomic_to_transcript(
        variant.start,
        transcript,
        upstream_window,
        downstream_window,
    );
    let end_pos = if variant.end != variant.start {
        map_genomic_to_transcript(variant.end, transcript, upstream_window, downstream_window)
    } else {
        None
    };

    // `upstream_gene_variant` and `downstream_gene_variant` are candidates only
    // for `within_feature => 0` (Utils/Constants.pm), and a span covering the
    // transcript is within it, so its endpoints in the flanks contribute nothing:
    // every term of such a span comes from the whole-span passes below, whichever
    // side of the neighbourhood bound the endpoints fall on.
    let (start_pos, end_pos) = if complete_overlap {
        (
            start_pos.filter(|p| !is_flanking_position(p)),
            end_pos.filter(|p| !is_flanking_position(p)),
        )
    } else {
        (start_pos, end_pos)
    };

    let mut consequences: ConsequenceList = SmallVec::new();
    let mut fields = PositionFields::default();
    {
        let mut ctx = ApplyPositionContext {
            variant,
            transcript,
            consequences: &mut consequences,
            fields: &mut fields,
            reference_fasta: config.reference_fasta.as_deref(),
        };
        match (start_pos.as_ref(), end_pos.as_ref()) {
            (None, None) => {
                if !complete_overlap {
                    return None;
                }
            }
            (Some(start_p), None) => {
                apply_position(start_p, &mut ctx, None, true);
            }
            (None, Some(end_p)) => {
                // Insertions are anchored at `start`; if that anchor is out of range,
                // ignore a lone `end` mapping.
                if insertion {
                    return None;
                }
                apply_position(end_p, &mut ctx, None, true);
            }
            (Some(start_p), Some(end_p)) => {
                let start_is_flank = is_flanking_position(start_p);
                let end_is_flank = is_flanking_position(end_p);
                let start_is_intron = matches!(start_p, TranscriptPosition::Intron { .. });
                let end_is_intron = matches!(end_p, TranscriptPosition::Intron { .. });

                if insertion && start_is_flank && end_is_flank {
                    // Pure-flank insertions use the start anchor distance.
                    apply_position(start_p, &mut ctx, Some(end_p), true);
                } else if insertion && (start_is_intron ^ end_is_intron) {
                    // The non-intronic side is the primary anchor, keeping the
                    // coding/exonic consequences while the paired intron context
                    // feeds the boundary checks.
                    let keep_start = !start_is_intron;
                    if keep_start {
                        apply_position(start_p, &mut ctx, Some(end_p), true);
                    } else {
                        apply_position(end_p, &mut ctx, Some(start_p), true);
                    }
                } else if insertion && start_is_intron && end_is_intron {
                    // For purely intronic insertions, Perl behavior is effectively
                    // single-anchor; using both sides overcalls mixed splice labels
                    // (e.g. donor_5th_base + donor_region).
                    apply_position(start_p, &mut ctx, Some(end_p), true);
                } else if start_is_flank ^ end_is_flank {
                    // Perl parity around transcript edges:
                    // - insertions straddling body/flank are treated as flank-only
                    // - non-insertions straddling body/flank are treated as body-only
                    let keep_start = if insertion {
                        start_is_flank
                    } else {
                        !start_is_flank
                    };
                    if keep_start {
                        apply_position(start_p, &mut ctx, Some(end_p), true);
                    } else {
                        apply_position(end_p, &mut ctx, Some(start_p), true);
                    }
                } else {
                    apply_position(start_p, &mut ctx, Some(end_p), true);
                    apply_position(end_p, &mut ctx, Some(start_p), false);
                }
            }
        }
    }

    // Perl's coding predicates, ported whole for every non-SNV allele and for
    // an SNV inside the start codon (`coding::perl_coding_terms`): the
    // per-endpoint passes above keep their display fields, and the coding terms
    // they emitted are replaced by the set Perl assigns.
    if let Some(terms) = crate::coding::perl_coding_terms(
        variant,
        transcript,
        variant.start,
        variant.end,
        config.reference_fasta.as_deref(),
    ) {
        consequences.retain(|c| !crate::coding::is_perl_coding_term(*c));
        for term in terms {
            push_unique(&mut consequences, term);
        }
    }

    // Perl reports DISTANCE as the smallest of the four gaps between either end
    // of the variant and either end of the transcript (BaseTranscriptVariation
    // distance_to_transcript), not the gap from one anchor position.
    if consequences.iter().any(|c| {
        matches!(
            c,
            Consequence::UpstreamGeneVariant | Consequence::DownstreamGeneVariant
        )
    }) {
        fields.distance = Some(distance_to_transcript(input, transcript));
    }

    // Perl checks splice sites against all overlapping introns
    // (_overlapped_introns / _overlapped_introns_boundary), while the per-endpoint
    // apply_position calls above cover only the introns containing the endpoints;
    // push_unique keeps the supplemental pass duplicate-free.
    if !insertion {
        add_splice_for_overlapping_introns(&mut consequences, variant, transcript);
        // Perl uses full-span overlap for UTR terms (_before_coding /
        // _after_coding), so a deletion from coding into UTR gets the UTR term.
        add_utr_for_overlapping_span(&mut consequences, variant, transcript);

        // Perl uses full-span exon overlap for non_coding_transcript_exon_variant
        // (_overlapped_exons); with both endpoints intronic the endpoint mapping
        // never enters the NonCodingExon branch, and the term added here makes
        // add_transcript_context_consequences skip the generic one.
        if !transcript.facts().has_coding_model {
            let lo = variant.start.min(variant.end);
            let hi = variant.start.max(variant.end);
            if allele_overlaps_any_exon(transcript, lo, hi) {
                let term = if variant_overlaps_mature_mirna(variant, transcript) {
                    Consequence::MatureMirnaVariant
                } else {
                    Consequence::NonCodingTranscriptExonVariant
                };
                push_unique(&mut consequences, term);
            }
        }
    } else {
        // Perl's _intron_effects iterates all overlapping introns for insertions
        // too; the per-endpoint apply_position calls cover only the endpoint introns.
        add_splice_for_overlapping_introns_insertion(&mut consequences, variant, transcript);
    }

    // Perl's `VariationEffect::splice_donor_region_variant` and `splice_region`
    // predicates read the transcript-wide `_intron_effects` flags: a 5th-base
    // hit in any intron suppresses the donor region, and a donor, acceptor,
    // 5th-base or donor-region hit in any intron suppresses splice_region,
    // whichever intron or exonic window raised it.
    if consequences.contains(&Consequence::SpliceDonor5thBaseVariant) {
        remove_consequence(&mut consequences, Consequence::SpliceDonorRegionVariant);
    }
    if consequences.iter().any(|c| {
        matches!(
            c,
            Consequence::SpliceDonorVariant
                | Consequence::SpliceAcceptorVariant
                | Consequence::SpliceDonor5thBaseVariant
                | Consequence::SpliceDonorRegionVariant
        )
    }) {
        remove_consequence(&mut consequences, Consequence::SpliceRegionVariant);
    }

    // Ensembl VEP never pairs mature_miRNA_variant with
    // non_coding_transcript_exon_variant; an insertion straddling the miRNA
    // subfeature boundary gets one from each apply_position, and the miRNA term
    // takes precedence.
    if consequences.contains(&Consequence::MatureMirnaVariant) {
        remove_consequence(
            &mut consequences,
            Consequence::NonCodingTranscriptExonVariant,
        );
    }

    // The predicates `VariationEffect` answers from `complete_overlap_feature`
    // alone: `coding_unknown` returns 0 before its peptide tests, so no
    // `coding_sequence_variant`; `non_coding_exon_variant` returns 0, so
    // `within_non_coding_gene` yields the transcript-context term instead; and
    // `coding_transcript_variant` is `!coding_unknown && complete_overlap_feature
    // && within_coding_gene`, included for biotype `protein_coding` and true
    // when the transcript has a translation. The splice, intron and UTR
    // predicates read the whole-span windows above unchanged.
    if complete_overlap {
        remove_consequence(&mut consequences, Consequence::CodingSequenceVariant);
        remove_consequence(
            &mut consequences,
            Consequence::NonCodingTranscriptExonVariant,
        );
        if transcript.facts().is_protein_coding && transcript.translation.is_some() {
            push_unique(&mut consequences, Consequence::CodingTranscriptVariant);
        }
    }

    // Perl `within_feature`: `overlap(vf.start, vf.end, tr.start, tr.end)` on the
    // raw coordinates, so an insertion butting against a transcript edge
    // (end = tr.start - 1 or start = tr.end + 1) is outside.
    let within_feature = variant.end >= transcript.start && variant.start <= transcript.end;
    add_transcript_context_consequences(transcript, within_feature, &mut consequences);

    // Ensembl VEP reads the intron and splice-site predicates over the regions
    // where the allele differs from the reference (BaseTranscriptVariationAllele.pm:117-218
    // over `_get_differing_regions`, VariationFeatureOverlapAllele.pm:405-460): an
    // allele of two or more bases equal to the reference allele, the reference
    // row's or an alternate allele equal to the transcript's sequence under
    // `--use_transcript_ref`, differs nowhere, so its row carries no intron or
    // splice term, while a one-base or empty allele keeps the one region its
    // span is. Such a row left with no term takes the default consequence on
    // its transcript (BaseVariationFeatureOverlapAllele.pm:285): an insertion's
    // empty reference allele in a coding exon at a codon boundary earns no
    // coding term, since its reference peptide is empty.
    if variant.reference_allele_row || variant.alt_allele() == variant.ref_allele.as_slice() {
        if variant.ref_allele.len() > 1 {
            consequences.retain(|c| !INTRON_EFFECT_TERMS.contains(c));
        }
        if consequences.is_empty() {
            consequences.push(Consequence::IntergenicVariant);
        }
    }

    if consequences.is_empty() {
        return None;
    }

    // Ensembl VEP outputs consequence terms in severity order (most severe first).
    consequences.sort_by_key(|c| c.rank());

    let impact = consequences
        .iter()
        .map(|c| c.impact())
        .min()
        .unwrap_or(Impact::MODIFIER);

    // HGVS strings are consumed solely by the output layer, so they are built
    // only when `compute_hgvs` is set.
    // The HGVS notations are built only `within_feature` (OutputFactory.pm:1721),
    // the predicate read at the shifted span, so an allele that shifts off the
    // transcript has none even where its input span overlapped it.
    let hgvs = if config.compute_hgvs && within_feature {
        crate::hgvs::generate_hgvs(
            input,
            transcript,
            config.reference_fasta.as_deref(),
            config.hgvs,
        )
    } else {
        crate::hgvs::HgvsNotation::default()
    };

    Some(TranscriptConsequence {
        consequences,
        impact,
        distance: fields.distance,
        hgvsc: hgvs.hgvsc,
        hgvsp: hgvs.hgvsp,
        hgvs_offset: hgvs.offset,
        ..TranscriptConsequence::for_transcript(transcript)
    })
}

fn is_flanking_position(position: &TranscriptPosition) -> bool {
    matches!(
        position,
        TranscriptPosition::Upstream { .. } | TranscriptPosition::Downstream { .. }
    )
}

/// Shared context for position-based consequence application.
struct ApplyPositionContext<'a> {
    variant: &'a InputVariant,
    transcript: &'a Transcript,
    consequences: &'a mut ConsequenceList,
    fields: &'a mut PositionFields,
    reference_fasta: Option<&'a vep_fasta::IndexedFasta>,
}

/// Perl's `deletion` pre-predicate for a sequence variant
/// (`BaseVariationFeatureOverlapAllele::_pre_consequence_predicates`): the
/// reference span `vf.end - vf.start + 1` (`_bvf_preds` `ref_length`) is longer
/// than the ALT sequence, `-` counting as empty. `snp` is the equal-length case
/// and `insertion` the longer ALT. A symbolic `<DEL>` never reaches this test:
/// it is a `StructuralVariationFeature`, whose `deletion` flag comes from its
/// class term in `_bvf_preds`, and here it takes the structural path.
fn perl_deletion_pre_predicate(variant: &InputVariant) -> bool {
    let ref_length = (variant.end as i64 - variant.start as i64 + 1).max(0);
    let alt = variant.alt_allele();
    let alt_length = if alt == b"-" { 0 } else { alt.len() as i64 };
    ref_length > alt_length
}

/// The smallest absolute gap between {variant start, variant end} and
/// {transcript start, transcript end}, which is what VEP prints as DISTANCE for
/// an up/downstream row. An adjacent base gives 1.
fn distance_to_transcript(variant: &InputVariant, transcript: &Transcript) -> u64 {
    let ends = [variant.start, variant.end];
    let bounds = [transcript.start, transcript.end];
    ends.iter()
        .flat_map(|v| bounds.iter().map(move |t| v.abs_diff(*t)))
        .min()
        .unwrap_or(0)
}

fn apply_position(
    position: &TranscriptPosition,
    ctx: &mut ApplyPositionContext<'_>,
    paired_position: Option<&TranscriptPosition>,
    populate_fields: bool,
) {
    let variant = ctx.variant;
    let transcript = ctx.transcript;
    let consequences = &mut *ctx.consequences;
    let fields = &mut *ctx.fields;
    let insertion = matches!(
        variant.variant_class,
        vep_core::variant::VariantClass::Insertion
    ) && variant.start == variant.end + 1;

    match position {
        TranscriptPosition::Upstream { distance: d } => {
            let mut out_distance = *d;
            if insertion
                && paired_position.is_some()
                && !matches!(
                    paired_position,
                    Some(
                        TranscriptPosition::Upstream { .. } | TranscriptPosition::Downstream { .. }
                    )
                )
            {
                out_distance = 0;
            }
            push_unique(consequences, Consequence::UpstreamGeneVariant);
            if populate_fields {
                fields.distance = Some(out_distance);
            }
        }
        TranscriptPosition::Downstream { distance: d } => {
            let mut out_distance = *d;
            if insertion
                && paired_position.is_some()
                && !matches!(
                    paired_position,
                    Some(
                        TranscriptPosition::Upstream { .. } | TranscriptPosition::Downstream { .. }
                    )
                )
            {
                out_distance = 0;
            }
            push_unique(consequences, Consequence::DownstreamGeneVariant);
            if populate_fields {
                fields.distance = Some(out_distance);
            }
        }
        TranscriptPosition::FivePrimeUtr { .. } => {
            push_unique(consequences, Consequence::FivePrimeUtrVariant);
            if insertion {
                if matches!(
                    paired_position,
                    Some(
                        TranscriptPosition::Intron {
                            dist_to_donor: 0,
                            ..
                        } | TranscriptPosition::Intron {
                            dist_to_acceptor: 0,
                            ..
                        }
                    )
                ) {
                    push_unique(consequences, Consequence::SpliceRegionVariant);
                }
                // Perl adds exonic splice_region only through _intron_effects() for an
                // insertion at an exon/intron boundary; with an intronic paired
                // position the Intron branch already evaluates the splice windows.
                let paired_is_intron =
                    matches!(paired_position, Some(TranscriptPosition::Intron { .. }));
                if !paired_is_intron
                    && !paired_intron_has_essential_splice(variant, transcript, paired_position)
                {
                    add_exonic_splice_region(consequences, variant, transcript);
                }
            } else if !paired_intron_has_essential_splice(variant, transcript, paired_position) {
                add_exonic_splice_region(consequences, variant, transcript);
            }
        }
        TranscriptPosition::ThreePrimeUtr { .. } => {
            push_unique(consequences, Consequence::ThreePrimeUtrVariant);
            if insertion {
                if matches!(
                    paired_position,
                    Some(
                        TranscriptPosition::Intron {
                            dist_to_donor: 0,
                            ..
                        } | TranscriptPosition::Intron {
                            dist_to_acceptor: 0,
                            ..
                        }
                    )
                ) {
                    push_unique(consequences, Consequence::SpliceRegionVariant);
                }
                // Perl: exonic splice_region for insertions comes from _intron_effects()
                // only. Skip when paired with an intron (intron branch handles it).
                let paired_is_intron =
                    matches!(paired_position, Some(TranscriptPosition::Intron { .. }));
                if !paired_is_intron
                    && !paired_intron_has_essential_splice(variant, transcript, paired_position)
                {
                    add_exonic_splice_region(consequences, variant, transcript);
                }
            } else if !paired_intron_has_essential_splice(variant, transcript, paired_position) {
                add_exonic_splice_region(consequences, variant, transcript);
            }
        }
        TranscriptPosition::Coding { cds_pos, .. } => {
            if insertion {
                // An insertion between the coding region's edge base and the base
                // outside it is a UTR variant in Perl: `_before_coding` and
                // `_after_coding` special-case an insertion whose start is the coding
                // region start or whose end is the coding region end, `within_cdna`
                // holds through the coding flank, and `_bvfo_preds` sets `exon` and
                // `utr` on the sorted flank pair. When the outer flank is intronic
                // this is the only apply_position call the allele gets.
                push_perl_utr_terms(consequences, variant, transcript);
            }

            // Ensembl VEP gives a variant spanning an intron boundary or extending
            // beyond the CDS (UTR, flanks) generic coding_sequence_variant instead
            // of codon analysis: one endpoint coding and the other intronic; both
            // exonic with an intervening intron; coding into UTR; coding past the
            // transcript.
            let spans_non_coding = !insertion
                && (matches!(
                    paired_position,
                    Some(
                        TranscriptPosition::Intron { .. }
                            | TranscriptPosition::ThreePrimeUtr { .. }
                            | TranscriptPosition::FivePrimeUtr { .. }
                            | TranscriptPosition::Downstream { .. }
                            | TranscriptPosition::Upstream { .. }
                    )
                ) || variant_span_overlaps_intron(variant, transcript));

            if spans_non_coding {
                // Perl's stop_lost/start_lost for these alleles go through
                // `_ins_del_stop_altered` / `_ins_del_start_altered`, which require
                // `VariationFeatureOverlapAllele::seq_is_unambiguous_dna` (true for
                // "-"; regex /^[ACGT-]+$/i), a length change (increase_length OR
                // decrease_length), and `_overlaps_stop_codon` /
                // `_overlaps_start_codon`, which return 0 for `cds_end_NF` /
                // `cds_start_NF`.
                let also_spans_intron =
                    matches!(paired_position, Some(TranscriptPosition::Intron { .. }))
                        || variant_span_overlaps_intron(variant, transcript);

                let cds_len = transcript
                    .vefc
                    .as_ref()
                    .and_then(|v| v.translateable_seq.as_ref())
                    .map(|s| s.len() as u64)
                    .unwrap_or(0);

                // For an MNV (same ref/alt lengths) the only Perl stop_lost path is
                // via defined peptide alleles, and `_get_peptide_alleles` is empty
                // past the CDS boundary because `codon()` requires `tv_tr_start &&
                // tv_tr_end` and `tv_tr_end` is undef.
                let ref_len = if variant.ref_allele.as_slice() == b"-" {
                    0usize
                } else {
                    variant.ref_allele.len()
                };
                let alt_slice = variant.alt_allele();
                let alt_len = if alt_slice == b"-" {
                    0usize
                } else {
                    alt_slice.len()
                };
                let is_length_changing = ref_len != alt_len;

                let cds_start_nf = transcript.facts().cds_start_nf;
                let cds_end_nf = transcript.facts().cds_end_nf;

                // A deletion from CDS into the 3' UTR covers the stop codon by
                // construction; Perl still requires a length change and excludes
                // `cds_end_NF` via `_overlaps_stop_codon`.
                //
                // `Downstream` is not a stop_lost shape. A reference span past the
                // transcript's 3' boundary leaves Perl's `cdna_end` undef (trailing
                // `Gap`), so `VariationEffect::_overlaps_stop_codon` returns 0, which
                // gates `_ins_del_stop_altered`, the only `stop_lost` path for these
                // alleles (`_get_peptide_alleles` is undef, so `stop_lost` delegates
                // to it); Perl emits `3_prime_UTR_variant,coding_sequence_variant`
                // via `coding_unknown`. The mapper here snaps the out-of-transcript
                // endpoint to the nearest exon boundary, so `cdna_hi` alone would
                // open the gate. Symmetric with the `extends_upstream_of_transcript`
                // exclusion on `is_start_lost`.
                //
                // The structural-variant arm of Perl's `stop_lost` (its
                // `TranscriptStructuralVariationAllele` branch) does emit stop_lost
                // for a Downstream-extending deletion with no cdna_end requirement,
                // but only symbolic/structural classes reach it: the dispatch in
                // `calculate_consequences` tests `variant.is_structural` (symbolic
                // ALTs) and `VariantClass::is_structural()`, and
                // `VariantClass::Deletion` is in neither, so an explicit-sequence
                // long deletion arrives here on Perl's `coding_unknown` path.
                let is_stop_lost = !also_spans_intron
                    && is_length_changing
                    && !cds_end_nf
                    && cds_len > 0
                    && *cds_pos <= cds_len
                    && matches!(
                        paired_position,
                        Some(TranscriptPosition::ThreePrimeUtr { .. })
                    );

                // start_lost, symmetric on the start codon side; excluded when
                // `cds_start_NF`.
                //
                // A reference span extending upstream of the transcript (forward:
                // `span_start < transcript.start`; reverse: `span_end >
                // transcript.end`) never gets `start_lost` from Perl: its mapper
                // yields a leading `Gap`, so `cds_start` and `cdna_start` are undef
                // and `_ins_del_start_altered` is bypassed, leaving
                // `5_prime_UTR_variant,coding_sequence_variant`. The mapper here
                // snaps the endpoint to `cdna_lo = 1`, so the shape is excluded from
                // `is_start_lost` and takes the `else` branch below.
                let extends_upstream_of_transcript = !insertion
                    && match transcript.strand {
                        vep_core::coordinate::Strand::Forward => variant.start < transcript.start,
                        vep_core::coordinate::Strand::Reverse => variant.end > transcript.end,
                    };

                // For the upstream-extension shape Perl emits
                // `5_prime_UTR_variant,coding_sequence_variant`. The 5'UTR term
                // normally comes from `add_utr_for_overlapping_span`, which
                // short-circuits when `genomic_coding_bounds` or
                // `transcript.translation` is None or the span fails
                // `overlaps_any_exon` (common in Storable-derived transcripts), so
                // it is also emitted here when the paired endpoint resolves to
                // 5'UTR or Upstream; `push_unique` keeps that idempotent.
                if extends_upstream_of_transcript
                    && matches!(
                        paired_position,
                        Some(
                            TranscriptPosition::FivePrimeUtr { .. }
                                | TranscriptPosition::Upstream { .. }
                        )
                    )
                {
                    push_unique(consequences, Consequence::FivePrimeUtrVariant);
                }

                let near_start_codon = cds_len > 0 && *cds_pos <= cds_len;
                let is_start_lost = !also_spans_intron
                    && is_length_changing
                    && !cds_start_nf
                    && near_start_codon
                    && !extends_upstream_of_transcript
                    && matches!(
                        paired_position,
                        Some(
                            TranscriptPosition::FivePrimeUtr { .. }
                                | TranscriptPosition::Upstream { .. }
                        )
                    );

                if is_stop_lost {
                    // Perl's stop_lost fires only when alt_pep lacks '*' while ref_pep
                    // has '*'. `_get_peptide_alleles` returns undef for a deletion
                    // extending past the CDS end, which routes Perl into
                    // `_ins_del_stop_altered`: it concatenates CDS + 3'UTR, applies
                    // the edit and re-translates the codon at the original stop, so
                    // shifted UTR bases spelling another stop yield
                    // `stop_retained_variant`. `compute_full_peptide_alleles_clamped()`
                    // would truncate the deletion at the CDS end and report StopLost,
                    // hiding that shift, so the past-CDS-end case takes
                    // `replicate_ins_del_stop_altered()`; without FASTA or cached UTR
                    // the clamped path remains.
                    let extends_past_cds_end = crate::mapper::map_genomic_span_to_cds_bounds(
                        variant.start,
                        variant.end,
                        transcript,
                    )
                    .is_some_and(|bounds| {
                        bounds.cds_start.max(bounds.cds_end) >= cds_len
                            && matches!(
                                paired_position,
                                Some(
                                    TranscriptPosition::ThreePrimeUtr { .. }
                                        | TranscriptPosition::Downstream { .. }
                                )
                            )
                    });

                    let utr_aware = || -> Option<Consequence> {
                        // The helper prefers the cached `vefc.three_prime_utr` and uses
                        // FASTA only when the cache lacks the field.
                        let fasta_opt = ctx.reference_fasta;
                        let bounds = crate::mapper::map_genomic_span_to_cds_bounds(
                            variant.start,
                            variant.end,
                            transcript,
                        )?;
                        let cds_lo = bounds.cds_start.min(bounds.cds_end);
                        if cds_lo == 0 || cds_lo > cds_len {
                            return None;
                        }
                        // Perl's `_ins_del_stop_altered` anchors the splice at
                        // `$cds_start - 1` but uses the full cDNA length
                        // `($cdna_end - $cdna_start) + 1`, not clamped to the CDS end:
                        // a 4bp deletion with 3bp in the stop codon and 1bp in the
                        // 3'UTR still splices 4 bytes out of CDS+UTR, so the re-read
                        // codon pulls in a UTR base.
                        let (cdna_lo, cdna_hi) = crate::mapper::map_genomic_span_to_cdna_bounds(
                            variant.start,
                            variant.end,
                            transcript,
                        )?;
                        let variant_cdna_len = (cdna_hi - cdna_lo + 1) as usize;
                        if variant_cdna_len == 0 {
                            return None;
                        }
                        crate::coding::replicate_ins_del_stop_altered(
                            variant,
                            transcript,
                            cds_lo,
                            variant_cdna_len,
                            fasta_opt,
                        )
                        .map(|altered| {
                            if altered {
                                Consequence::StopLost
                            } else {
                                Consequence::StopRetainedVariant
                            }
                        })
                    };

                    let consequence = if extends_past_cds_end {
                        // The clamped peptide would produce a spurious StopLost here.
                        utr_aware()
                            .or_else(|| {
                                // No FASTA and no cached UTR: the clamped peptide path.
                                crate::coding::compute_full_peptide_alleles_clamped(
                                    variant, transcript, *cds_pos,
                                )
                                .and_then(|(ref_pep, alt_pep)| {
                                    classify_stop_from_peptides(&ref_pep, &alt_pep)
                                })
                            })
                            .unwrap_or(Consequence::StopLost)
                    } else {
                        crate::coding::compute_full_peptide_alleles(variant, transcript, *cds_pos)
                            .and_then(|(ref_pep, alt_pep)| {
                                classify_stop_from_peptides(&ref_pep, &alt_pep)
                            })
                            .or_else(|| {
                                crate::coding::compute_full_peptide_alleles_clamped(
                                    variant, transcript, *cds_pos,
                                )
                                .and_then(|(ref_pep, alt_pep)| {
                                    classify_stop_from_peptides(&ref_pep, &alt_pep)
                                })
                            })
                            .or_else(utr_aware)
                            .unwrap_or(Consequence::StopLost)
                    };
                    push_unique(consequences, consequence);
                } else if is_start_lost {
                    // Perl's `VariationEffect::_ins_del_start_altered`
                    // checks DNA-level ATG integrity: when the triplet at
                    // `length($utr->seq)` is still 'ATG' after the edit and the
                    // 5'UTR bases are unchanged, Perl returns 0 and suppresses
                    // `start_lost`; the widened gate alone would fire it for every
                    // UTR-spanning length-changing deletion. The helper prefers the
                    // cached `vefc.five_prime_utr` and uses FASTA only as fallback.
                    // Upstream-extending deletions never reach here (see
                    // `extends_upstream_of_transcript` above).
                    let utr_aware_start = || -> Option<bool> {
                        let fasta_opt = ctx.reference_fasta;
                        let (cdna_lo, cdna_hi) = crate::mapper::map_genomic_span_to_cdna_bounds(
                            variant.start,
                            variant.end,
                            transcript,
                        )?;
                        let variant_cdna_len = (cdna_hi - cdna_lo + 1) as usize;
                        if variant_cdna_len == 0 {
                            return None;
                        }
                        crate::coding::replicate_ins_del_start_altered(
                            variant,
                            transcript,
                            cdna_lo,
                            variant_cdna_len,
                            fasta_opt,
                        )
                    };

                    match utr_aware_start() {
                        Some(true) => {
                            push_unique(consequences, Consequence::StartLost);
                        }
                        Some(false) => {
                            // Start preserved: Perl emits no start_lost; the 5'UTR term
                            // comes from the paired_position UTR emission.
                        }
                        None => {
                            // No UTR sequence available: the widened gate stands.
                            push_unique(consequences, Consequence::StartLost);
                        }
                    }
                } else {
                    // `frameshift_variant` / `inframe_insertion` / `inframe_deletion`
                    // for the `[Coordinate, Gap, Coordinate]` shape (a deletion
                    // spanning introns with both endpoints in exonic CDS): Perl's
                    // `VariationEffect::frameshift` fires when `defined cds_start &&
                    // defined cds_end`, which `BaseTranscriptVariation::cds_start` /
                    // `cds_end` take from the first and last `cds_coords` entries; a
                    // Gap at either end short-circuits to 0. `coding_sequence_variant`
                    // (`coding_unknown`) fires only when frameshift / inframe_* /
                    // protein_altering / start_* / stop_* all return 0, so the
                    // specific term suppresses it. The gate reads the projection's
                    // segment shape rather than
                    // `map_genomic_span_to_cds_bounds`, which clamps intronic
                    // endpoints and cannot tell them from coding ones.
                    let mut emitted_specific_term = false;
                    if let Some(proj) = crate::mapper::map_genomic_span_to_cds_projection(
                        variant.start,
                        variant.end,
                        transcript,
                    ) {
                        // Perl gate 1: `defined cds_start && defined cds_end`, both
                        // endpoints `EndpointClass::InCds`.
                        let perl_gate_1 = proj.both_endpoints_in_cds();

                        // Perl gate 2: `return 0 if partial_codon`.
                        let perl_gate_2 = !proj.has_partial_terminal_codon;

                        // Perl gate 3: `_overlaps_stop_codon` / `_overlaps_start_codon`
                        // return 0 for cds_end_NF / cds_start_NF, and the frameshift
                        // chain reaches one of them via `stop_retained`.
                        let perl_gate_3 = !proj.cds_start_nf && !proj.cds_end_nf;

                        // Perl gate 4: the is_stop_lost / is_start_lost branches above
                        // already handled the exon-to-UTR and exon-to-flank cases.
                        let perl_gate_4 = !is_stop_lost && !is_start_lost;

                        if perl_gate_1 && perl_gate_2 && perl_gate_3 && perl_gate_4 {
                            if let Some((cs, ce)) = proj.coordinate_bounds() {
                                let var_len = ce.saturating_sub(cs) + 1;
                                // Perl: allele_len = bvfoa->seq_length, 0 for a pure deletion.
                                let alt_bytes = variant.alt_allele();
                                let allele_len: u64 = if alt_bytes == b"-" {
                                    0
                                } else {
                                    alt_bytes.len() as u64
                                };
                                let diff = var_len.abs_diff(allele_len);
                                if diff % 3 != 0 {
                                    // Perl `frameshift`: abs(allele_len - var_len) % 3
                                    push_unique(consequences, Consequence::FrameshiftVariant);
                                    emitted_specific_term = true;
                                } else if allele_len != var_len {
                                    if allele_len > var_len {
                                        push_unique(consequences, Consequence::InframeInsertion);
                                    } else {
                                        push_unique(consequences, Consequence::InframeDeletion);
                                    }
                                    emitted_specific_term = true;
                                }
                                // allele_len == var_len: Perl emits only
                                // coding_sequence_variant for a same-length change in the
                                // spans_non_coding branch.
                            }
                        }

                        // Perl's `VariationEffect::partial_codon` is both a suppressor
                        // (`perl_gate_2` above) and a positive emitter of
                        // `incomplete_terminal_codon_variant`, including for spans that
                        // leave the CDS (`CDS_position = N-?`). The `segments.first()`
                        // guard ports its early return (`return 0 unless defined
                        // $bvfo->translation_start`): a leading Gap suppresses the
                        // term. `segments.first()` is 5'->3' transcript order, which a
                        // genomic-ordered endpoint would invert on reverse strand.
                        if proj.has_partial_terminal_codon
                            && matches!(
                                proj.segments.first(),
                                Some(crate::mapper::MapperSegment::Coordinate { .. })
                            )
                        {
                            push_unique(consequences, Consequence::IncompleteTerminalCodonVariant);
                        }
                    }

                    // `VariationEffect::coding_unknown` fires only when no specific
                    // coding term has been emitted.
                    if !emitted_specific_term {
                        push_unique(consequences, Consequence::CodingSequenceVariant);
                    }
                }
            } else {
                let coding_bounds = crate::mapper::map_genomic_span_to_cds_bounds(
                    variant.start,
                    variant.end,
                    transcript,
                );
                let use_insertion_anchor =
                    should_use_insertion_anchor(insertion, coding_bounds.as_ref());
                let analysis_cds_pos = if use_insertion_anchor {
                    insertion_analysis_cds_pos(*cds_pos, coding_bounds.as_ref())
                } else {
                    *cds_pos
                };

                // Perl emits `incomplete_terminal_codon_variant` + `coding_sequence_variant`
                // for a variant in the last partial codon of a CDS whose length is not
                // divisible by 3.
                let in_incomplete_terminal =
                    is_in_incomplete_terminal_codon(transcript, analysis_cds_pos);

                // Perl evaluates the coding predicates once per variant-transcript
                // pair over the full span, so the secondary endpoint of a multi-base
                // edit (which would apply the same edit at another cds_pos) is skipped.
                {
                    let ref_a = &variant.ref_allele;
                    let alt_a = variant.alt_allele();
                    let is_mnv = ref_a.len() > 1
                        && alt_a.len() > 1
                        && ref_a.len() == alt_a.len()
                        && ref_a != b"-"
                        && alt_a != b"-";
                    let is_length_changing_indel = ref_a.len() != alt_a.len();
                    if (is_mnv || is_length_changing_indel || insertion) && !populate_fields {
                        return;
                    }
                }

                if let Some(cc) = get_codon_change(variant, transcript, analysis_cds_pos) {
                    // Perl's frameshift predicate uses the CDS-mapped length
                    // (cds_end - cds_start + 1), which differs from the raw ref length
                    // at exon/intron boundaries and UTR spans and changes the mod-3
                    // result; raw lengths apply only when no bounds exist.
                    let is_frameshift = match coding_bounds.as_ref() {
                        Some(bounds) => crate::coding::is_frameshift_cds_aware(variant, bounds),
                        None => cc.is_frameshift,
                    };

                    if in_incomplete_terminal {
                        push_unique(consequences, Consequence::IncompleteTerminalCodonVariant);
                        push_unique(consequences, Consequence::CodingSequenceVariant);
                    } else if is_frameshift {
                        let frameshift_bounds = coding_bounds;

                        // Perl's `frameshift` returns 0 if `stop_retained` is true, so a
                        // stop-preserving variant at the stop codon is inframe_insertion +
                        // stop_retained_variant. The codon window extends by the alt
                        // allele length, so a large insertion near the stop includes it;
                        // the full-CDS check applies only when bounds cannot be computed.
                        let stop_retained = if let Some(b) = frameshift_bounds.as_ref() {
                            is_stop_retained_insertion(variant, transcript, b, true)
                        } else {
                            is_stop_retained_insertion_full_cds(variant, transcript, *cds_pos)
                        };

                        if stop_retained {
                            push_unique(consequences, Consequence::InframeInsertion);
                            push_unique(consequences, Consequence::StopRetainedVariant);
                        } else {
                            push_unique(consequences, Consequence::FrameshiftVariant);
                            // Perl's frameshift stop_gained reads the local codon-window
                            // peptides, where the stop can sit in a later codon of the
                            // window (`AGC` -> `ACTTAGC` => `T*X`), so a single-codon AA
                            // check is insufficient.
                            if let Some(bounds) = frameshift_bounds.as_ref() {
                                if frameshift_stop_gained_in_codon_window(
                                    variant, transcript, bounds,
                                ) {
                                    push_unique(consequences, Consequence::StopGained);
                                }
                            } else {
                                // Without codon-window bounds (span at a CDS/UTR boundary),
                                // the full-CDS translation with the same scan limit.
                                if frameshift_stop_gained_full_cds(
                                    variant,
                                    transcript,
                                    analysis_cds_pos,
                                ) {
                                    push_unique(consequences, Consequence::StopGained);
                                }
                            }
                        }
                    } else {
                        // Perl's codon windows use translation_start/end from the full
                        // variant span, not a single cds_pos.
                        assign_coding_consequence(
                            consequences,
                            cc.ref_amino_acid,
                            cc.alt_amino_acid,
                            analysis_cds_pos,
                            variant,
                            transcript,
                            coding_bounds.as_ref(),
                        );

                        // Perl's stop predicates are peptide-based and fire for inframe
                        // indels too. `stop_retained` -> `ref_eq_alt_sequence` reads
                        // codon-window peptides (`_get_peptide_alleles`), not the full CDS,
                        // whose terminal `*` would make every distant insertion
                        // stop_retained: the stop must appear in the local window.
                        if populate_fields
                            && consequences.contains(&Consequence::InframeInsertion)
                            && !consequences.contains(&Consequence::StopRetainedVariant)
                        {
                            // Codon-window check first; near the stop codon the full-CDS
                            // check follows (Perl's condition 2 in ref_eq_alt_sequence, full
                            // peptide substitution).
                            let stop_retained_inframe = if let Some(b) = coding_bounds.as_ref() {
                                if is_stop_retained_insertion(variant, transcript, b, false) {
                                    true
                                } else {
                                    // Perl's _overlaps_stop_codon: the variant cDNA overlaps
                                    // the last 3 CDS positions.
                                    let near_stop = transcript
                                        .vefc
                                        .as_ref()
                                        .and_then(|v| v.translateable_seq.as_ref())
                                        .map(|s| {
                                            let cds_len = s.len() as u64;
                                            cds_len >= 3 && *cds_pos >= cds_len - 2
                                        })
                                        .unwrap_or(false);
                                    near_stop
                                        && is_stop_retained_insertion_full_cds(
                                            variant, transcript, *cds_pos,
                                        )
                                }
                            } else {
                                is_stop_retained_insertion_full_cds(
                                    variant,
                                    transcript,
                                    analysis_cds_pos,
                                )
                            };
                            if stop_retained_inframe {
                                push_unique(consequences, Consequence::StopRetainedVariant);
                            }
                        }

                        // stop_gained from codon-window peptides: the full CDS always
                        // reaches the transcript's terminal stop after an insertion.
                        if populate_fields
                            && consequences.contains(&Consequence::InframeInsertion)
                            && !consequences.contains(&Consequence::StopGained)
                            && !consequences.contains(&Consequence::StopRetainedVariant)
                        {
                            if let Some(b) = coding_bounds.as_ref() {
                                if let Some((ref_pep, alt_pep)) =
                                    compute_codon_window_peptide_alleles(variant, transcript, b)
                                {
                                    if alt_pep.contains(&b'*') && !ref_pep.contains(&b'*') {
                                        push_unique(consequences, Consequence::StopGained);
                                    }
                                }
                            }
                        }

                        // stop_gained/stop_lost from codon-window peptides: the full
                        // translation always includes the terminal stop, which flank
                        // trimming can expose asymmetrically for a distant deletion.
                        if populate_fields
                            && consequences.contains(&Consequence::InframeDeletion)
                            && !consequences.contains(&Consequence::StopGained)
                            && !consequences.contains(&Consequence::StopLost)
                            && !consequences.contains(&Consequence::StopRetainedVariant)
                        {
                            if let Some(b) = coding_bounds.as_ref() {
                                if let Some((ref_pep, alt_pep)) =
                                    compute_codon_window_peptide_alleles(variant, transcript, b)
                                {
                                    let ref_has_stop = ref_pep.contains(&b'*');
                                    let alt_has_stop = alt_pep.contains(&b'*');
                                    if alt_has_stop && !ref_has_stop {
                                        push_unique(consequences, Consequence::StopGained);
                                    } else if ref_has_stop && !alt_has_stop {
                                        push_unique(consequences, Consequence::StopLost);
                                    }
                                }
                            }
                        }

                        // No StartLost here: the start-codon block above covers
                        // InframeDeletion and runs the alt-CDS ATG scanner before
                        // deciding, so a push here would short-circuit its
                        // StartRetainedVariant emission.
                        //
                        // Perl's stop predicates fire for any coding variant, so a
                        // ProteinAltering alt peptide with a premature stop still gets
                        // stop_gained, scoped by the codon window.
                        if populate_fields
                            && consequences.contains(&Consequence::ProteinAlteringVariant)
                            && !consequences.contains(&Consequence::StopGained)
                            && !consequences.contains(&Consequence::StopLost)
                            && !consequences.contains(&Consequence::StopRetainedVariant)
                        {
                            if let Some(b) = coding_bounds.as_ref() {
                                if let Some((ref_pep, alt_pep)) =
                                    compute_codon_window_peptide_alleles(variant, transcript, b)
                                {
                                    let ref_has_stop = ref_pep.contains(&b'*');
                                    let alt_has_stop = alt_pep.contains(&b'*');
                                    if alt_has_stop && !ref_has_stop {
                                        push_unique(consequences, Consequence::StopGained);
                                    } else if ref_has_stop && !alt_has_stop {
                                        push_unique(consequences, Consequence::StopLost);
                                    }
                                }
                            }
                        }

                        // An MNV spanning codons can change a stop across a codon boundary
                        // (GC->AA producing MQ->I*), which the single-codon AA comparison
                        // misses; indels are excluded because a frame change always hits
                        // a stop.
                        let ref_allele = &variant.ref_allele;
                        let alt_allele_bytes = variant.alt_allele();
                        let is_mnv = ref_allele.len() > 1
                            && alt_allele_bytes.len() > 1
                            && ref_allele.len() == alt_allele_bytes.len()
                            && ref_allele != b"-"
                            && alt_allele_bytes != b"-";
                        if populate_fields
                            && is_mnv
                            && !consequences.contains(&Consequence::StopGained)
                            && !consequences.contains(&Consequence::StopLost)
                            && !consequences.contains(&Consequence::StopRetainedVariant)
                        {
                            if let Some(bounds) = crate::mapper::map_genomic_span_to_cds_bounds(
                                variant.start,
                                variant.end,
                                transcript,
                            ) {
                                if let Some((ref_pep, alt_pep)) =
                                    compute_codon_window_peptide_alleles(
                                        variant, transcript, &bounds,
                                    )
                                {
                                    let ref_has_stop = ref_pep.contains(&b'*');
                                    let alt_has_stop = alt_pep.contains(&b'*');
                                    if alt_has_stop && !ref_has_stop {
                                        push_unique(consequences, Consequence::StopGained);
                                    } else if ref_has_stop && !alt_has_stop {
                                        push_unique(consequences, Consequence::StopLost);
                                    }
                                }
                            }
                        }

                        // Perl suppresses missense/synonymous when stop_gained or stop_lost
                        // is true; both anchors of a multi-base allele are evaluated here,
                        // so the invariant is enforced after stop detection.
                        if consequences.contains(&Consequence::StopGained)
                            || consequences.contains(&Consequence::StopLost)
                        {
                            remove_consequence(consequences, Consequence::MissenseVariant);
                            remove_consequence(consequences, Consequence::SynonymousVariant);
                        }
                    }
                } else if in_incomplete_terminal {
                    push_unique(consequences, Consequence::IncompleteTerminalCodonVariant);
                    push_unique(consequences, Consequence::CodingSequenceVariant);
                } else {
                    push_unique(consequences, Consequence::CodingSequenceVariant);
                }
            }

            // Perl co-emits frameshift_variant + stop_lost for a terminal frameshift
            // but drops the generic inframe labels once stop_lost fires.
            if !consequences.contains(&Consequence::StopLost)
                && !consequences.contains(&Consequence::StopGained)
                && !consequences.contains(&Consequence::StopRetainedVariant)
                && (consequences.contains(&Consequence::FrameshiftVariant)
                    || consequences.contains(&Consequence::InframeDeletion)
                    || consequences.contains(&Consequence::InframeInsertion))
            {
                if let Some(cds_len) = transcript
                    .vefc
                    .as_ref()
                    .and_then(|v| v.translateable_seq.as_ref())
                    .map(|s| s.len() as u64)
                {
                    if cds_len >= 3 {
                        let stop_codon_start = cds_len - 2; // 1-based
                        let bounds = crate::mapper::map_genomic_span_to_cds_bounds(
                            variant.start,
                            variant.end,
                            transcript,
                        );
                        let use_insertion_anchor =
                            should_use_insertion_anchor(insertion, bounds.as_ref());
                        // The 5'-most CDS index of the span, not the per-endpoint start
                        // anchor (`*cds_pos`), which on reverse strand is the 3'-most
                        // index and makes `compute_peptide_alleles` return None.
                        let stop_loss_cds_pos = if use_insertion_anchor {
                            insertion_analysis_cds_pos(*cds_pos, bounds.as_ref())
                        } else {
                            analysis_cds_pos_5prime(*cds_pos, bounds.as_ref())
                        };
                        // The stop-loss test reads both the start position and the CDS
                        // span end, so a deletion starting before the stop codon but
                        // extending through it counts.
                        let span_overlaps_stop = stop_loss_cds_pos >= stop_codon_start
                            || bounds
                                .as_ref()
                                .is_some_and(|b| b.cds_end >= stop_codon_start);
                        if span_overlaps_stop {
                            if let Some((ref_pep, alt_pep)) =
                                compute_peptide_alleles(variant, transcript, stop_loss_cds_pos)
                            {
                                let ref_has_stop = ref_pep.contains(&b'*');
                                let alt_has_stop = alt_pep.contains(&b'*');
                                if ref_has_stop && !alt_has_stop {
                                    push_unique(consequences, Consequence::StopLost);
                                    remove_consequence(consequences, Consequence::InframeDeletion);
                                    remove_consequence(consequences, Consequence::InframeInsertion);
                                }
                            }
                        }
                    }
                }
            }

            // Perl fires start_lost when the first codon does not translate to M,
            // and start_retained when the original ATG triplet is still present in
            // the altered codon window. The gate uses the 5'-most CDS index of the
            // span: on a reverse-strand deletion `*cds_pos` is the 3'-most index (4
            // for a CDS 2-4 span) and the `<= 3` gate would miss it.
            let start_bounds = if populate_fields {
                crate::mapper::map_genomic_span_to_cds_bounds(
                    variant.start,
                    variant.end,
                    transcript,
                )
            } else {
                None
            };
            let start_5prime_cds_pos = analysis_cds_pos_5prime(*cds_pos, start_bounds.as_ref());
            if populate_fields
                && start_5prime_cds_pos <= 3
                && !consequences.contains(&Consequence::StartLost)
                && !consequences.contains(&Consequence::StartRetainedVariant)
                && (consequences.contains(&Consequence::FrameshiftVariant)
                    || consequences.contains(&Consequence::InframeInsertion)
                    || consequences.contains(&Consequence::InframeDeletion)
                    || consequences.contains(&Consequence::ProteinAlteringVariant))
            {
                let ref_is_met = transcript
                    .vefc
                    .as_ref()
                    .and_then(|v| v.translateable_seq.as_ref())
                    .map(|s| s.len() >= 3 && s.as_bytes()[..3].eq_ignore_ascii_case(b"ATG"))
                    .unwrap_or(false);

                if ref_is_met {
                    // The alt-CDS ATG scan runs first: Perl's start_retained fires when
                    // an in-frame ATG remains at the canonical start regardless of
                    // whether the first-codon peptide translates to M (a 1bp deletion
                    // of the G in ATG when codon 2 begins with G leaves ATG at position
                    // 0), and a StartLost emitted before the scan would suppress it.
                    // Only codon-aligned position 0 counts; a downstream ATG is a
                    // coincidental Met residue, not the preserved start codon.
                    let alt_has_codon_zero_atg = transcript
                        .vefc
                        .as_ref()
                        .and_then(|v| v.translateable_seq.as_ref())
                        .map(|cds_seq| {
                            let cds = cds_seq.as_bytes();
                            let ref_allele = &variant.ref_allele;
                            let alt_allele = variant.alt_allele();
                            let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
                            let ref_len = if ref_is_dash { 0 } else { ref_allele.len() };
                            let alt_seq = crate::coding::feature_seq(
                                alt_allele,
                                variant.strand,
                                transcript.strand,
                            );
                            // `idx` derives from the 5'-most CDS index, as the gate above
                            // does: the per-endpoint `*cds_pos` is the 3'-most index on
                            // reverse strand and would excise the wrong bases, preserving
                            // the ATG prefix and falsely reading start_retained.
                            let idx = (start_5prime_cds_pos - 1) as usize;
                            if idx + ref_len > cds.len() {
                                return false;
                            }
                            let mut alt_cds = Vec::with_capacity(
                                cds.len().saturating_sub(ref_len) + alt_seq.len(),
                            );
                            alt_cds.extend_from_slice(&cds[..idx]);
                            alt_cds.extend_from_slice(&alt_seq);
                            alt_cds.extend_from_slice(&cds[idx + ref_len..]);
                            if alt_cds.len() < 3 {
                                return false;
                            }
                            let alt_upper: [u8; 3] = [
                                alt_cds[0].to_ascii_uppercase(),
                                alt_cds[1].to_ascii_uppercase(),
                                alt_cds[2].to_ascii_uppercase(),
                            ];
                            &alt_upper == b"ATG"
                        })
                        .unwrap_or(false);

                    if alt_has_codon_zero_atg {
                        push_unique(consequences, Consequence::StartRetainedVariant);
                    } else {
                        // StartLost when the first-codon peptide does not translate to
                        // M and no codon-aligned ATG remains;
                        // `compute_codon_window_peptide_alleles` reads the strand-normalized
                        // `start_bounds`, not the per-endpoint `*cds_pos`.
                        if let Some(bounds) = start_bounds.as_ref() {
                            if let Some((ref_pep, alt_pep)) =
                                crate::coding::compute_codon_window_peptide_alleles(
                                    variant, transcript, bounds,
                                )
                            {
                                let ref_starts_m = ref_pep.first() == Some(&b'M');
                                let alt_starts_m = alt_pep.first() == Some(&b'M');
                                if ref_starts_m && !alt_starts_m {
                                    push_unique(consequences, Consequence::StartLost);
                                }
                            }
                        }
                    }
                }
            }

            if insertion {
                if matches!(
                    paired_position,
                    Some(
                        TranscriptPosition::Intron {
                            dist_to_donor: 0,
                            ..
                        } | TranscriptPosition::Intron {
                            dist_to_acceptor: 0,
                            ..
                        }
                    )
                ) {
                    push_unique(consequences, Consequence::SpliceRegionVariant);
                }
                // Perl: exonic splice_region for insertions comes from _intron_effects()
                // only. Skip when paired with an intron (intron branch handles it).
                let paired_is_intron =
                    matches!(paired_position, Some(TranscriptPosition::Intron { .. }));
                if !paired_is_intron
                    && !paired_intron_has_essential_splice(variant, transcript, paired_position)
                {
                    add_exonic_splice_region(consequences, variant, transcript);
                }
            } else if !paired_intron_has_essential_splice(variant, transcript, paired_position) {
                add_exonic_splice_region(consequences, variant, transcript);
            }
        }
        TranscriptPosition::Intron {
            intron_number,
            dist_to_donor,
            dist_to_acceptor,
            ..
        } => {
            let intron_bounds = get_intron_bounds(transcript, *intron_number);

            if let Some((intron_start, intron_end)) = intron_bounds {
                if is_frameshift_intron(intron_start, intron_end) {
                    // Perl's _intron_effects sets within_frameshift_intron for a
                    // frameshift intron (<= 12 bp) and skips every intron and splice
                    // term. On a coding transcript `within_cds` then accepts the span
                    // only where it overlaps the coding region
                    // (coding_sequence_variant; a paired exonic endpoint adds the
                    // specific terms in its own apply_position call). A frameshift
                    // intron on the UTR side of the CDS gets the UTR term instead:
                    // `_bvfo_preds` sets `utr` for a span clear of the coding region,
                    // the 12 bp exon stretch satisfies `exon`, and `within_cdna`
                    // accepts the intron through `within_transcript`.
                    if transcript.facts().has_coding_model {
                        if allele_overlaps_coding_bounds(variant, transcript) {
                            push_unique(consequences, Consequence::CodingSequenceVariant);
                        } else {
                            push_perl_utr_terms(consequences, variant, transcript);
                        }
                    } else {
                        // Non-coding transcript: `VariationEffect::non_coding_exon_variant`
                        // needs the raw span to overlap an exon, which a variant inside
                        // the intron does not, so `within_non_coding_gene` yields
                        // non_coding_transcript_variant.
                        let lo = variant.start.min(variant.end);
                        let hi = variant.start.max(variant.end);
                        if !insertion && allele_overlaps_any_exon(transcript, lo, hi) {
                            push_unique(consequences, Consequence::NonCodingTranscriptExonVariant);
                        } else if transcript.facts().is_non_coding_context {
                            push_unique(consequences, Consequence::NonCodingTranscriptVariant);
                        }
                    }
                    return;
                }
            }

            let lo = variant.start.min(variant.end);
            let hi = variant.start.max(variant.end);
            let overlaps_exon = !insertion && allele_overlaps_any_exon(transcript, lo, hi);
            let overlaps_coding =
                overlaps_exon && allele_overlaps_coding_region(transcript, lo, hi);
            if overlaps_coding {
                push_unique(consequences, Consequence::CodingSequenceVariant);
            }

            let spans_exon_boundary = !insertion
                && matches!(
                    paired_position,
                    Some(
                        TranscriptPosition::Coding { .. }
                            | TranscriptPosition::FivePrimeUtr { .. }
                            | TranscriptPosition::ThreePrimeUtr { .. }
                            | TranscriptPosition::NonCodingExon { .. }
                    )
                );
            let paired_intron_same = matches!(
                paired_position,
                Some(TranscriptPosition::Intron {
                    intron_number: paired_intron_number,
                    ..
                }) if *paired_intron_number == *intron_number
            );
            if !insertion && (spans_exon_boundary || paired_intron_same) {
                if let Some((intron_start, intron_end)) = intron_bounds {
                    // Perl computes intronic splice labels by overlap across the full
                    // allele span; evaluating both endpoints independently over-calls
                    // mixed labels (donor_5th_base + donor_region).
                    if paired_intron_same && !populate_fields {
                        return;
                    }
                    // Perl's `include` gate: splice_polypyrimidine_tract_variant requires
                    // `exon => 0`, so suppress polypyr when the variant also overlaps an exon.
                    // Perl also extends the exon overlap window by 12bp when the transcript
                    // has any frameshift intron (`BaseTranscriptVariation::_overlapped_exons`).
                    let lo = variant.start.min(variant.end);
                    let hi = variant.start.max(variant.end);
                    let frameshift_stretch = !spans_exon_boundary
                        && !overlaps_exon
                        && transcript.facts().vefc_has_frameshift_intron
                        && allele_overlaps_any_exon_stretched(transcript, lo, hi);
                    let suppress_polypyr =
                        spans_exon_boundary || overlaps_exon || frameshift_stretch;
                    add_intronic_non_insertion_splice_consequences(
                        consequences,
                        transcript,
                        variant,
                        intron_start,
                        intron_end,
                        suppress_polypyr,
                    );
                    return;
                }
            }

            // For intronic non-insertions Perl uses overlap-based intron effects;
            // that path applies whenever intron bounds exist.
            if !insertion {
                if let Some((intron_start, intron_end)) = intron_bounds {
                    let lo2 = variant.start.min(variant.end);
                    let hi2 = variant.start.max(variant.end);
                    let fs_stretch = !overlaps_exon
                        && transcript.facts().vefc_has_frameshift_intron
                        && allele_overlaps_any_exon_stretched(transcript, lo2, hi2);
                    add_intronic_non_insertion_splice_consequences(
                        consequences,
                        transcript,
                        variant,
                        intron_start,
                        intron_end,
                        overlaps_exon || fs_stretch,
                    );
                    return;
                }
            }

            if insertion
                && matches!(
                    paired_position,
                    Some(TranscriptPosition::Intron {
                        intron_number: paired_intron_number,
                        ..
                    }) if *paired_intron_number == *intron_number
                )
            {
                if let Some((intron_start, intron_end)) = intron_bounds {
                    add_intronic_insertion_splice_consequences(
                        consequences,
                        transcript,
                        variant,
                        intron_start,
                        intron_end,
                    );
                    return;
                }
            }

            // Ensembl VEP treats insertions exactly at exon/intron boundaries as
            // non-essential for donor/acceptor calls (splice_region-only).
            let boundary_insertion_against_exon = insertion
                && paired_position.is_some()
                && !matches!(paired_position, Some(TranscriptPosition::Intron { .. }));
            let suppress_donor = boundary_insertion_against_exon && *dist_to_donor == 0;
            let suppress_acceptor = boundary_insertion_against_exon && *dist_to_acceptor == 0;
            if boundary_insertion_against_exon && (*dist_to_donor == 0 || *dist_to_acceptor == 0) {
                // Insertions exactly at exon/intron joins are non-essential splice
                // events but still counted as splice_region in Ensembl VEP.
                push_unique(consequences, Consequence::SpliceRegionVariant);
            }

            add_splice_consequences(
                consequences,
                *dist_to_donor,
                *dist_to_acceptor,
                suppress_donor,
                suppress_acceptor,
            );
        }
        TranscriptPosition::NonCodingExon { .. } => {
            let insertion_at_exon_intron_boundary = insertion
                && matches!(
                    paired_position,
                    Some(
                        TranscriptPosition::Intron {
                            dist_to_donor: 0,
                            ..
                        } | TranscriptPosition::Intron {
                            dist_to_acceptor: 0,
                            ..
                        }
                    )
                );

            if variant_overlaps_mature_mirna(variant, transcript) {
                push_unique(consequences, Consequence::MatureMirnaVariant);
            } else if !insertion_at_exon_intron_boundary {
                push_unique(consequences, Consequence::NonCodingTranscriptExonVariant);
            }

            if insertion {
                if matches!(
                    paired_position,
                    Some(
                        TranscriptPosition::Intron {
                            dist_to_donor: 0,
                            ..
                        } | TranscriptPosition::Intron {
                            dist_to_acceptor: 0,
                            ..
                        }
                    )
                ) {
                    push_unique(consequences, Consequence::SpliceRegionVariant);
                }
                // Perl: exonic splice_region for insertions comes from _intron_effects()
                // only. Skip when paired with an intron (intron branch handles it).
                let paired_is_intron =
                    matches!(paired_position, Some(TranscriptPosition::Intron { .. }));
                if !paired_is_intron
                    && !paired_intron_has_essential_splice(variant, transcript, paired_position)
                {
                    add_exonic_splice_region(consequences, variant, transcript);
                }
            } else if !paired_intron_has_essential_splice(variant, transcript, paired_position) {
                add_exonic_splice_region(consequences, variant, transcript);
            }
        }
    }
}

fn add_transcript_context_consequences(
    transcript: &Transcript,
    within_feature: bool,
    consequences: &mut ConsequenceList,
) {
    // The context predicates (`within_nmd_transcript`, `within_non_coding_gene`)
    // are gated on `within_feature => 1`; a variant lying wholly inside a
    // frameshift intron of a non-coding transcript reaches here with no other
    // term and still gets its context term.
    if !within_feature {
        return;
    }

    // Ensembl VEP does not generally pair transcript-context labels with
    // `non_coding_transcript_exon_variant` rows, so keep those exon rows clean.
    if consequences.contains(&Consequence::NonCodingTranscriptExonVariant) {
        return;
    }

    // Ensembl VEP also keeps `mature_miRNA_variant` rows clean (no additional
    // transcript-context labels).
    if consequences.contains(&Consequence::MatureMirnaVariant) {
        return;
    }

    let facts = transcript.facts();
    if facts.is_nmd {
        push_unique(consequences, Consequence::NmdTranscriptVariant);
    } else if facts.is_non_coding_context && !facts.has_coding_model {
        push_unique(consequences, Consequence::NonCodingTranscriptVariant);
    }
}

fn get_intron_bounds(transcript: &Transcript, intron_number: u32) -> Option<(u64, u64)> {
    let introns = transcript
        .vefc
        .as_ref()
        .filter(|v| !v.introns.is_empty())
        .map(|v| v.introns.as_slice())
        .unwrap_or(transcript.introns.as_slice());
    introns
        .iter()
        .find(|i| i.rank == intron_number)
        .map(|i| (i.start, i.end))
}

fn allele_overlaps_any_exon(transcript: &Transcript, lo: u64, hi: u64) -> bool {
    let exons = transcript
        .vefc
        .as_ref()
        .map(|v| v.sorted_exons.as_slice())
        .unwrap_or(transcript.exons.as_slice());
    exons.iter().any(|e| hi >= e.start && lo <= e.end)
}

/// Check if the variant overlaps any exon with a 12bp stretch (Perl's frameshift-intron
/// exon extension). Suppresses polypyrimidine for transcripts with frameshift introns.
fn allele_overlaps_any_exon_stretched(transcript: &Transcript, lo: u64, hi: u64) -> bool {
    let exons = transcript
        .vefc
        .as_ref()
        .map(|v| v.sorted_exons.as_slice())
        .unwrap_or(transcript.exons.as_slice());
    exons
        .iter()
        .any(|e| hi >= e.start.saturating_sub(12) && lo <= e.end.saturating_add(12))
}

/// Check if the variant's full genomic span overlaps any intron of the transcript.
/// Used for large deletions that span from one exon through introns to another exon.
fn variant_span_overlaps_intron(variant: &InputVariant, transcript: &Transcript) -> bool {
    let lo = variant.start.min(variant.end);
    let hi = variant.start.max(variant.end);
    let introns = transcript
        .vefc
        .as_ref()
        .filter(|v| !v.introns.is_empty())
        .map(|v| v.introns.as_slice())
        .unwrap_or(transcript.introns.as_slice());
    introns.iter().any(|i| hi >= i.start && lo <= i.end)
}

fn allele_overlaps_coding_region(transcript: &Transcript, lo: u64, hi: u64) -> bool {
    let (Some(cs), Some(ce)) = (transcript.coding_region_start, transcript.coding_region_end)
    else {
        return false;
    };
    let start = cs.min(ce);
    let end = cs.max(ce);
    hi >= start && lo <= end
}

/// Perl's `overlap` (`Utils/VariationEffect.pm`) on the allele's own coordinates,
/// which for an insertion run `start == end + 1`: an insertion overlaps a region
/// only when both flanks lie inside it, so one sitting on the region's edge does not.
fn perl_overlap(f1_start: u64, f1_end: u64, f2_start: u64, f2_end: u64) -> bool {
    f1_end >= f2_start && f1_start <= f2_end
}

/// Perl's `within_cds` frameshift-intron arm and the `coding` pre-predicate of
/// `_bvfo_preds`: the allele overlaps the genomic coding region.
fn allele_overlaps_coding_bounds(variant: &InputVariant, transcript: &Transcript) -> bool {
    match genomic_coding_bounds(transcript) {
        Some((cds_start, cds_end)) => perl_overlap(variant.start, variant.end, cds_start, cds_end),
        None => false,
    }
}

/// The positional half of Perl's `within_5_prime_utr` and `within_3_prime_utr`
/// (`Utils/VariationEffect.pm`). `_before_coding` is an overlap with
/// `[transcript start, coding_region_start - 1]`, or an insertion whose start is
/// the coding region start; `_after_coding` is an overlap with
/// `[coding_region_end + 1, transcript end]`, or an insertion whose end is the
/// coding region end. Both read genomic coordinates, and the strand names the
/// UTR. Pushes whichever terms hold; the caller has settled `within_cdna`.
fn push_perl_utr_terms(
    consequences: &mut ConsequenceList,
    variant: &InputVariant,
    transcript: &Transcript,
) {
    let Some((cds_start, cds_end)) = genomic_coding_bounds(transcript) else {
        return;
    };
    let (s, e) = (variant.start, variant.end);
    let insertion = s == e + 1;
    let before_coding = (insertion && s == cds_start)
        || perl_overlap(s, e, transcript.start, cds_start.saturating_sub(1));
    let after_coding = (insertion && e == cds_end)
        || perl_overlap(s, e, cds_end.saturating_add(1), transcript.end);
    let (five_prime, three_prime) = match transcript.strand {
        vep_core::coordinate::Strand::Forward => (before_coding, after_coding),
        vep_core::coordinate::Strand::Reverse => (after_coding, before_coding),
    };
    if five_prime {
        push_unique(consequences, Consequence::FivePrimeUtrVariant);
    }
    if three_prime {
        push_unique(consequences, Consequence::ThreePrimeUtrVariant);
    }
}

fn paired_intron_has_essential_splice(
    variant: &InputVariant,
    transcript: &Transcript,
    paired_position: Option<&TranscriptPosition>,
) -> bool {
    let intron_number = match paired_position {
        Some(TranscriptPosition::Intron { intron_number, .. }) => *intron_number,
        _ => return false,
    };
    let (intron_start, intron_end) = match get_intron_bounds(transcript, intron_number) {
        Some(v) => v,
        None => return false,
    };
    let lo = variant.start.min(variant.end);
    let hi = variant.start.max(variant.end);
    let (donor_s, donor_e, acceptor_s, acceptor_e) = match transcript.strand {
        vep_core::coordinate::Strand::Forward => (
            intron_start,
            intron_start.saturating_add(1),
            intron_end.saturating_sub(1),
            intron_end,
        ),
        vep_core::coordinate::Strand::Reverse => (
            intron_end.saturating_sub(1),
            intron_end,
            intron_start,
            intron_start.saturating_add(1),
        ),
    };
    let overlaps =
        |a_start: u64, a_end: u64, b_start: u64, b_end: u64| a_start <= b_end && a_end >= b_start;
    overlaps(lo, hi, donor_s, donor_e) || overlaps(lo, hi, acceptor_s, acceptor_e)
}

/// Assign the specific coding consequence based on amino acid change.
///
/// Uses `compute_peptide_alleles` for multi-base variants (MNVs, complex indels)
/// so that stop_gained/stop_lost detection works across codon boundaries,
/// matching Ensembl VEP's peptide-level `_get_peptide_alleles` approach.
fn assign_coding_consequence(
    consequences: &mut ConsequenceList,
    ref_aa: u8,
    alt_aa: u8,
    cds_pos: u64,
    variant: &InputVariant,
    transcript: &Transcript,
    cds_bounds: Option<&CdsSpanBounds>,
) {
    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();

    // Perl's `TranscriptVariationAllele::peptide` returns undef when
    // `VariationFeatureOverlapAllele::seq_is_unambiguous_dna` is false (allele
    // contains non-ACGT characters like N, R, Y), and `_get_peptide_alleles` needs
    // the reference allele's peptide as well as the alternate's
    // (Utils/VariationEffect.pm:772-790), so a reference base of `N` leaves every
    // peptide predicate false: Ensembl VEP 116.2 writes `coding_sequence_variant`
    // for `N G` where `A G` at the same base is `missense_variant`.
    let ambiguous = |allele: &[u8]| {
        !allele.is_empty() && allele != b"-" && !crate::coding::is_unambiguous_dna(allele)
    };
    if ambiguous(alt_allele) || ambiguous(ref_allele) {
        push_unique(consequences, Consequence::CodingSequenceVariant);
        return;
    }

    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();

    // Pure insertions and complex indels with a net length change take
    // classify_inframe_indel_by_peptide's peptide-level containment checks,
    // matching the containment gate in Perl's
    // `VariationEffect::protein_altering_variant`.
    if ref_is_dash && !alt_is_dash {
        let kind = classify_inframe_indel_by_peptide(variant, transcript, cds_pos, cds_bounds);
        match kind {
            Some(InframeIndelKind::LeadingStopGained) => {
                push_unique(consequences, Consequence::StopGained);
            }
            Some(InframeIndelKind::ProteinAltering) => {
                push_unique(consequences, Consequence::ProteinAlteringVariant);
            }
            Some(InframeIndelKind::Insertion) => {
                push_unique(consequences, Consequence::InframeInsertion);
            }
            Some(InframeIndelKind::Deletion) => {
                push_unique(consequences, Consequence::InframeInsertion);
            }
            None => {
                // Perl falls back to coding_sequence_variant when the peptide
                // comparison fails (no translateable_seq, invalid CDS bounds).
                push_unique(consequences, Consequence::CodingSequenceVariant);
            }
        }
        return;
    }
    if !ref_is_dash && alt_is_dash {
        // Perl classifies a pure deletion as protein_altering_variant when the alt
        // peptide does not contain the ref as a sub/super-string.
        let kind = classify_inframe_indel_by_peptide(variant, transcript, cds_pos, cds_bounds);
        match kind {
            Some(InframeIndelKind::LeadingStopGained) => {
                push_unique(consequences, Consequence::StopGained);
            }
            Some(InframeIndelKind::ProteinAltering) => {
                push_unique(consequences, Consequence::ProteinAlteringVariant);
            }
            Some(InframeIndelKind::Deletion) | Some(InframeIndelKind::Insertion) => {
                push_unique(consequences, Consequence::InframeDeletion);
            }
            None => {
                // Peptide comparison unavailable: inframe_deletion, the structural
                // classification.
                push_unique(consequences, Consequence::InframeDeletion);
            }
        }
        return;
    }
    if ref_allele.len() != alt_allele.len() {
        let kind = classify_inframe_indel_by_peptide(variant, transcript, cds_pos, cds_bounds);
        match kind {
            Some(InframeIndelKind::Insertion) => {
                push_unique(consequences, Consequence::InframeInsertion);
            }
            Some(InframeIndelKind::Deletion) => {
                push_unique(consequences, Consequence::InframeDeletion);
            }
            Some(InframeIndelKind::LeadingStopGained) => {
                push_unique(consequences, Consequence::StopGained);
            }
            Some(InframeIndelKind::ProteinAltering) => {
                push_unique(consequences, Consequence::ProteinAlteringVariant);
            }
            None => {
                // Perl falls back to coding_sequence_variant when translateable_seq
                // is unavailable or CDS bounds fail.
                push_unique(consequences, Consequence::CodingSequenceVariant);
            }
        }
        return;
    }

    // Perl's coding_unknown fires when either peptide contains 'X' (an
    // ambiguous codon, or an `X` SeqEdit on the reference residue). Its
    // exclusion list names neither stop_gained nor missense_variant, and
    // `missense_variant` is plain `ref_pep ne alt_pep`, so an edited `X`
    // reference co-emits one of them with coding_sequence_variant.
    if ref_aa == b'X' || alt_aa == b'X' {
        push_unique(consequences, Consequence::CodingSequenceVariant);
        if alt_aa == b'*' && ref_aa != b'*' {
            push_unique(consequences, Consequence::StopGained);
        } else if ref_aa != alt_aa && ref_aa != b'*' && alt_aa != b'*' {
            push_unique(consequences, Consequence::MissenseVariant);
        }
        return;
    }

    let ref_is_stop = ref_aa == b'*';
    let alt_is_stop = alt_aa == b'*';

    if ref_is_stop && alt_is_stop {
        push_unique(consequences, Consequence::StopRetainedVariant);
    } else if alt_is_stop {
        push_unique(consequences, Consequence::StopGained);
    } else if ref_is_stop {
        push_unique(consequences, Consequence::StopLost);
    } else if cds_pos <= 3 && is_start_codon_affected(ref_aa, cds_pos) {
        if ref_aa == alt_aa {
            push_unique(consequences, Consequence::StartRetainedVariant);
        } else {
            push_unique(consequences, Consequence::StartLost);
        }
    } else if ref_aa == alt_aa {
        push_unique(consequences, Consequence::SynonymousVariant);
    } else {
        push_unique(consequences, Consequence::MissenseVariant);
    }
}

/// Check if a variant affects the start codon.
///
/// Returns true only when the CDS position is within the first codon (positions 1-3)
/// and the reference amino acid is Met, confirming a canonical ATG start. This avoids
/// misclassifying variants in non-ATG start codons or transcripts with `cds_start_NF`.
fn is_start_codon_affected(ref_aa: u8, cds_pos: u64) -> bool {
    (1..=3).contains(&cds_pos) && ref_aa == b'M'
}

/// Check if a variant is an insertion/indel at/near the stop codon that preserves the stop.
///
/// Perl's `stop_retained` predicate suppresses `frameshift` when an insertion
/// at or just before the stop codon does not destroy the downstream stop
/// (e.g., `taa/taAGTGa` → `*/*VX`).
///
/// This uses the **codon-window peptide** (local codon context around the variant),
/// not the full-CDS peptide, whose terminal `*` would fire for an insertion anywhere
/// in the gene.
///
/// This also handles complex indels (ref and alt both have sequence) where the
/// alt is longer than the ref (net insertion), which VEP classifies as Indel
/// rather than Insertion.
pub(crate) fn is_stop_retained_insertion(
    variant: &InputVariant,
    transcript: &Transcript,
    bounds: &CdsSpanBounds,
    for_frameshift: bool,
) -> bool {
    use vep_core::variant::VariantClass;

    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();
    let is_net_insertion = matches!(variant.variant_class, VariantClass::Insertion)
        || ref_allele == b"-"
        || ref_allele.is_empty()
        || (alt_allele.len() > ref_allele.len() && !ref_allele.is_empty());
    if !is_net_insertion {
        return false;
    }

    // Perl's `TranscriptVariationAllele::peptide` reads the local codon window;
    // the stop must appear there to count as retained.
    let Some((ref_pep, alt_pep)) =
        compute_codon_window_peptide_alleles(variant, transcript, bounds)
    else {
        return false;
    };

    // Perl guard in `VariationEffect::ref_eq_alt_sequence`:
    //   return 0 if $ref_pep eq "X" && $alt_pep eq "X";  # incomplete coding terminal
    // Both peptides a single 'X' means an incomplete coding terminal that cannot
    // be translated either way, and Perl emits no stop_retained. The wider codon
    // window here produces 'X' at boundary codons more readily, so the guard fires
    // only when both resolve to exactly one 'X'.
    if ref_pep == b"X" && alt_pep == b"X" {
        return false;
    }

    // Perl's `stop_retained` → `ref_eq_alt_sequence` has three conditions (any true → retained):
    //
    // 1. ref_pep == first char of alt_pep AND alt_pep contains '*'
    //    (large insertion preserves the ref AA and introduces a stop in the inserted sequence)
    //
    // 2. The full translation is preserved after substitution with overflow < 3 chars
    //    (checked below against the full CDS translation, near the CDS end)
    //
    // 3. Both ref and alt contain '*' at the same peptide position
    //    (stop codon is at the same offset in both peptides)
    //
    // Condition 1: Perl's `$ref_pep eq substr($alt_pep, 0, 1)` requires the entire
    // ref_pep to be one character equal to the first character of alt_pep, so it
    // fires only when ref_pep is exactly 1 AA. CdsSpanBounds can produce a wider
    // window than Perl's, whose alt peptide may include shifted-frame stops Perl
    // never sees, so the scan is capped to ref AAs + inserted AAs, with no +1 for
    // the boundary codon: Perl's window does not extend into it.
    if ref_pep.len() == 1 && !alt_pep.is_empty() && ref_pep[0] == alt_pep[0] {
        let ref_len_nt = if ref_allele == b"-" || ref_allele.is_empty() {
            0usize
        } else {
            ref_allele.len()
        };
        let net_insertion_nt = alt_allele.len().saturating_sub(ref_len_nt);
        // For a frameshift (net insertion not divisible by 3) the partial boundary
        // codon translates to 'X' in Perl's `TranscriptVariationAllele::peptide`,
        // hiding shifted-frame stops, so floor division excludes it; an inframe
        // indel scans the full window.
        let insertion_aa_count = if for_frameshift {
            net_insertion_nt / 3 // floor: exclude partial boundary codon
        } else {
            net_insertion_nt.div_ceil(3)
        };
        let scan_limit = ref_pep.len() + insertion_aa_count;
        let scan_end = scan_limit.min(alt_pep.len());
        if alt_pep[..scan_end].contains(&b'*') {
            return true;
        }
    }

    // Condition 3: both peptides contain '*' (a small insertion near the stop
    // codon, where both windows include the stop).
    if ref_pep.contains(&b'*') && alt_pep.contains(&b'*') {
        // Perl checks same position: index(ref, '*') == index(alt, '*')
        let ref_stop_pos = ref_pep.iter().position(|&b| b == b'*');
        let alt_stop_pos = alt_pep.iter().position(|&b| b == b'*');
        if ref_stop_pos == alt_stop_pos {
            return true;
        }
    }

    // Condition 2: full-peptide substitution check.
    //
    // Perl's ref_eq_alt_sequence condition 2:
    //   $ref_seq = $bvfo->_peptide;  # full protein (from cache, without terminal *)
    //   $mut_seq = $ref_seq;
    //   substr($mut_seq, $tl_start-1, $tl_end - $tl_start + 1) = $alt_pep;
    //   $mut_substring = substr($mut_seq, 0, length($ref_seq));
    //   $final_stop = substr($mut_seq, length($ref_seq));
    //   return 1 if ($ref_seq eq $mut_substring && length($final_stop) < 3);
    //
    // Perl's cached peptide (_peptide) does not include the terminal `*`, so the
    // substitution check compares against a protein ending at the last AA: for an
    // insertion in the penultimate codon the first ref_len chars of mut_seq still
    // match because the stop was never there. Hence the CDS is translated, the
    // terminal `*` stripped, and the codon-window alt_pep substituted at
    // translation_start. Only the last ~15 CDS positions are checked: an insertion
    // further upstream displaces AAs and cannot satisfy the equality.
    //
    // Perl's stop_retained checks `return 0 if stop_lost(@_)` first, and stop_lost
    // reads the codon-window peptides (`ref_pep =~ /\*/ && alt_pep !~ /\*/`), so an
    // insertion in the stop codon destroys it and skips condition 2.
    let ref_has_stop = ref_pep.contains(&b'*');
    let alt_has_stop = alt_pep.contains(&b'*');
    let is_stop_lost = ref_has_stop && !alt_has_stop;
    if !alt_pep.is_empty() && !is_stop_lost {
        if let Some(vefc) = transcript.vefc.as_ref() {
            if let Some(tseq) = vefc.translateable_seq.as_ref() {
                let cds = tseq.as_bytes();
                let cds_len = cds.len();
                let near_stop =
                    cds_len >= 3 && bounds.cds_start as usize >= cds_len.saturating_sub(15);
                if near_stop {
                    let tl_start = bounds.translation_start as usize;
                    let tl_end = bounds.translation_end as usize;
                    if tl_start > 0 && tl_end > 0 {
                        let mut ref_full_pep = crate::coding::translate_to_peptide(cds);
                        if ref_full_pep.last() == Some(&b'*') {
                            ref_full_pep.pop();
                        }
                        let ref_len = ref_full_pep.len();

                        // Perl: substr($mut_seq, $tl_start-1, $tl_end - $tl_start + 1) = $alt_pep
                        let sub_start = tl_start - 1; // 0-based index
                        let sub_len = if tl_end >= tl_start {
                            tl_end - tl_start + 1
                        } else {
                            0 // boundary insertion: tl_start > tl_end
                        };
                        if sub_start <= ref_len {
                            let sub_end = (sub_start + sub_len).min(ref_len);
                            let mut mut_seq = Vec::with_capacity(ref_len + alt_pep.len());
                            mut_seq.extend_from_slice(&ref_full_pep[..sub_start]);
                            mut_seq.extend_from_slice(&alt_pep);
                            if sub_end < ref_len {
                                mut_seq.extend_from_slice(&ref_full_pep[sub_end..]);
                            }

                            if mut_seq.len() >= ref_len {
                                let matches_ref = mut_seq[..ref_len] == ref_full_pep[..];
                                let overflow = mut_seq.len() - ref_len;
                                if matches_ref && overflow < 3 {
                                    return true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    false
}

/// Full-CDS variant of `is_stop_retained_insertion`, for the **inframe** supplement path.
///
/// For inframe insertions (length divisible by 3), Perl's stop_retained predicate fires
/// whenever both the ref and alt full-CDS peptides contain a stop codon. This is correct
/// because an inframe insertion shifts the stop downstream but preserves it.
///
/// This is distinct from the codon-window version used in the frameshift branch, where
/// the narrower window correctly prevents false stop_retained on frameshifts far from
/// the stop codon.
fn is_stop_retained_insertion_full_cds(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
) -> bool {
    use crate::coding::compute_full_peptide_alleles;
    use vep_core::variant::VariantClass;

    // Perl skips stop_retained analysis for an incomplete CDS, whose truncated
    // sequence gives an unreliable stop position.
    let facts = transcript.facts();
    if facts.cds_start_nf || facts.cds_end_nf {
        return false;
    }

    let has_translateable_seq = transcript
        .vefc
        .as_ref()
        .and_then(|v| v.translateable_seq.as_ref())
        .map(|s| !s.is_empty())
        .unwrap_or(false);
    if !has_translateable_seq {
        return false;
    }

    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();
    let is_net_insertion = matches!(variant.variant_class, VariantClass::Insertion)
        || ref_allele == b"-"
        || ref_allele.is_empty()
        || (alt_allele.len() > ref_allele.len() && !ref_allele.is_empty());
    if !is_net_insertion {
        return false;
    }

    let Some((ref_pep, alt_pep)) = compute_full_peptide_alleles(variant, transcript, cds_pos)
    else {
        return false;
    };

    ref_pep.contains(&b'*') && alt_pep.contains(&b'*')
}

/// Check if a CDS position falls within an incomplete terminal codon.
///
/// An incomplete terminal codon exists when the translateable sequence length
/// is not divisible by 3. Variants in the last 1-2 bases get this annotation
/// instead of a standard coding consequence.
pub(crate) fn is_in_incomplete_terminal_codon(transcript: &Transcript, cds_pos: u64) -> bool {
    let cds_len = transcript
        .vefc
        .as_ref()
        .and_then(|v| v.translateable_seq.as_ref())
        .map(|s| s.len() as u64)
        .unwrap_or(0);
    if cds_len == 0 || cds_len.is_multiple_of(3) {
        return false;
    }
    let last_complete_codon_end = (cds_len / 3) * 3;
    cds_pos > last_complete_codon_end
}

/// Add splice site consequences for intronic variants.
fn add_splice_consequences(
    consequences: &mut ConsequenceList,
    dist_to_donor: u64,
    dist_to_acceptor: u64,
    suppress_donor: bool,
    suppress_acceptor: bool,
) {
    let donor_splice_site = dist_to_donor <= 1 && !suppress_donor;
    let acceptor_splice_site = dist_to_acceptor <= 1 && !suppress_acceptor;
    let donor_5th_base = dist_to_donor == 4;

    // Ensembl VEP donor-region window: intron positions +3..+6.
    // With 0-based intronic distance from donor, that's 2..=5.
    // `splice_donor_5th_base_variant` takes precedence over donor-region.
    let donor_region = (2..=5).contains(&dist_to_donor) && !donor_5th_base;

    // Ensembl VEP polypyrimidine tract window: `$intron_end - 16 .. $intron_end - 2`,
    // offsets from the last intron base. With that base at distance 0, that's 2..=16.
    let polypyrimidine = (2..=16).contains(&dist_to_acceptor);

    if donor_splice_site {
        push_unique(consequences, Consequence::SpliceDonorVariant);
    }
    if acceptor_splice_site {
        push_unique(consequences, Consequence::SpliceAcceptorVariant);
    }

    if donor_5th_base {
        push_unique(consequences, Consequence::SpliceDonor5thBaseVariant);
    }
    if donor_region {
        push_unique(consequences, Consequence::SpliceDonorRegionVariant);
    }
    if polypyrimidine {
        push_unique(consequences, Consequence::SplicePolypyrimidineTractVariant);
    }

    // Ensembl VEP precedence: splice_region is suppressed by donor/acceptor,
    // donor_region, and donor_5th_base classifications.
    let donor_side_splice_region = (2..=7).contains(&dist_to_donor);
    let acceptor_side_splice_region = (2..=7).contains(&dist_to_acceptor);
    if (donor_side_splice_region || acceptor_side_splice_region)
        && !donor_splice_site
        && !acceptor_splice_site
        && !donor_region
        && !donor_5th_base
    {
        push_unique(consequences, Consequence::SpliceRegionVariant);
    }
}

fn overlap_unsorted(v_start: u64, v_end: u64, f_start: u64, f_end: u64) -> bool {
    v_end >= f_start && v_start <= f_end
}

fn overlap_sorted(v_start: u64, v_end: u64, f_start: u64, f_end: u64) -> bool {
    let lo = v_start.min(v_end);
    let hi = v_start.max(v_end);
    hi >= f_start && lo <= f_end
}

fn intron_overlap_for_insertion(
    v_start: u64,
    v_end: u64,
    intron_start: u64,
    intron_end: u64,
) -> bool {
    overlap_unsorted(
        v_start,
        v_end,
        intron_start.saturating_add(2),
        intron_start.saturating_add(7),
    ) || overlap_unsorted(
        v_start,
        v_end,
        intron_end.saturating_sub(7),
        intron_end.saturating_sub(2),
    ) || overlap_unsorted(
        v_start,
        v_end,
        intron_start.saturating_sub(3),
        intron_start.saturating_sub(1),
    ) || overlap_unsorted(
        v_start,
        v_end,
        intron_end.saturating_add(1),
        intron_end.saturating_add(3),
    ) || v_start == intron_start
        || v_end == intron_end
        || v_start == intron_start.saturating_add(2)
        || v_end == intron_end.saturating_sub(2)
}

/// Splice terms for one intron of a pure insertion, following Perl's
/// `_intron_effects`: the single differing region is the inverted anchor pair
/// (`start = end + 1`), so the non-normalising `overlap` holds only when both
/// anchors sit inside a window, and the `$insertion` equality clauses admit the
/// exact boundary positions. As for non-insertions, the intron must be in the
/// memoised intron list, the boundary flags need the boundary list, and a
/// frameshift intron (span <= 12) the anchors overlap is skipped (the `next` in
/// `BaseTranscriptVariationAllele::_intron_effects`).
fn add_intronic_insertion_splice_consequences(
    consequences: &mut ConsequenceList,
    transcript: &Transcript,
    variant: &InputVariant,
    intron_start: u64,
    intron_end: u64,
) {
    let (v_start, v_end) = (variant.start, variant.end);
    let (in_intron_list, in_boundary_list) = perl_intron_list_membership(
        v_start.min(v_end),
        v_start.max(v_end),
        intron_start,
        intron_end,
    );
    if !in_intron_list && !in_boundary_list {
        return;
    }
    if is_frameshift_intron(intron_start, intron_end)
        && overlap_unsorted(v_start, v_end, intron_start, intron_end)
    {
        return;
    }

    let boundary_flag =
        |fs: u64, fe: u64| -> bool { in_boundary_list && overlap_unsorted(v_start, v_end, fs, fe) };
    let start_splice_site = boundary_flag(intron_start, intron_start.saturating_add(1));
    let end_splice_site = boundary_flag(intron_end.saturating_sub(1), intron_end);
    let fifth_base_splice_site = boundary_flag(
        intron_start.saturating_add(4),
        intron_start.saturating_add(4),
    );
    let donor_region_splice_site = boundary_flag(
        intron_start.saturating_add(2),
        intron_start.saturating_add(5),
    );
    // The intron loop tests the polypyrimidine windows on the sorted anchors,
    // so either anchor inside the window counts.
    let polypyrimidine_splice_site = overlap_sorted(
        v_start,
        v_end,
        intron_end.saturating_sub(16),
        intron_end.saturating_sub(2),
    );
    let fifth_base_splice_site_reverse =
        boundary_flag(intron_end.saturating_sub(4), intron_end.saturating_sub(4));
    let donor_region_splice_site_reverse =
        boundary_flag(intron_end.saturating_sub(5), intron_end.saturating_sub(2));
    let polypyrimidine_splice_site_reverse = overlap_sorted(
        v_start,
        v_end,
        intron_start.saturating_add(2),
        intron_start.saturating_add(16),
    );

    let (donor_splice_site, acceptor_splice_site, donor_5th_base, donor_region, polypyrimidine) =
        match transcript.strand {
            vep_core::coordinate::Strand::Forward => (
                start_splice_site,
                end_splice_site,
                fifth_base_splice_site,
                donor_region_splice_site,
                polypyrimidine_splice_site,
            ),
            vep_core::coordinate::Strand::Reverse => (
                end_splice_site,
                start_splice_site,
                fifth_base_splice_site_reverse,
                donor_region_splice_site_reverse,
                polypyrimidine_splice_site_reverse,
            ),
        };

    if donor_splice_site {
        push_unique(consequences, Consequence::SpliceDonorVariant);
    }
    if acceptor_splice_site {
        push_unique(consequences, Consequence::SpliceAcceptorVariant);
    }

    if donor_5th_base {
        push_unique(consequences, Consequence::SpliceDonor5thBaseVariant);
    }
    let donor_region = donor_region && !donor_5th_base;
    if donor_region {
        push_unique(consequences, Consequence::SpliceDonorRegionVariant);
    }
    // Perl's `include` gate for polypyrimidine is `exon => 0`: PPT is filtered
    // when the pre-predicate `exon=1`, and `_overlapped_exons` stretches exon
    // bounds by 12bp on a transcript with any frameshift intron (span <= 12bp),
    // so a variant near but outside an exon can still register `exon=1`
    // (`BaseTranscriptVariation::_overlapped_exons`,
    // `BaseVariationFeatureOverlapAllele::_skip_oc`).
    let lo = variant.start.min(variant.end);
    let hi = variant.start.max(variant.end);
    let overlaps_exon = if transcript.facts().vefc_has_frameshift_intron {
        allele_overlaps_any_exon_stretched(transcript, lo, hi)
    } else {
        allele_overlaps_any_exon(transcript, lo, hi)
    };
    if polypyrimidine && !overlaps_exon {
        push_unique(consequences, Consequence::SplicePolypyrimidineTractVariant);
    }

    // Perl `within_intron` tests the core intron (intron_start+2..intron_end-2)
    // plus the two anchor equalities; for an intron shorter than 4bp the
    // equalities still hold.
    let core_start = intron_start.saturating_add(2);
    let core_end = intron_end.saturating_sub(2);
    let within_core = overlap_unsorted(v_start, v_end, core_start, core_end)
        || v_start == core_start
        || v_end == core_end;
    if in_intron_list && within_core {
        push_unique(consequences, Consequence::IntronVariant);
    }

    let splice_region =
        in_boundary_list && intron_overlap_for_insertion(v_start, v_end, intron_start, intron_end);
    if splice_region
        && !donor_splice_site
        && !acceptor_splice_site
        && !donor_region
        && !donor_5th_base
    {
        push_unique(consequences, Consequence::SpliceRegionVariant);
    }
}

/// For non-insertion variants, add splice consequences for every intron Perl's
/// memoised intron list holds (`perl_intron_list_membership`), not only the
/// introns containing the variant endpoints. Perl reaches `_intron_effects`
/// only when `within_feature` holds for the raw variant span
/// (`BaseVariationFeatureOverlapAllele::_bvfo_preds`), and the inserted bases
/// of a delins play no part in which introns are fetched.
///
/// This is called after the per-endpoint `apply_position` calls; `push_unique`
/// keeps the re-evaluation of the endpoint introns duplicate-free.
fn add_splice_for_overlapping_introns(
    consequences: &mut ConsequenceList,
    variant: &InputVariant,
    transcript: &Transcript,
) {
    let vefc = match transcript.vefc.as_ref() {
        Some(v) => v,
        None => return,
    };
    let introns = if !vefc.introns.is_empty() {
        &vefc.introns
    } else {
        &transcript.introns
    };

    let vf_lo = variant.start.min(variant.end);
    let vf_hi = variant.start.max(variant.end);
    if vf_hi < transcript.start || vf_lo > transcript.end {
        return;
    }

    let overlaps_exon = allele_overlaps_any_exon(transcript, vf_lo, vf_hi);
    // Perl extends exon overlap by 12bp when the transcript has any frameshift intron.
    let suppress_polypyr = overlaps_exon
        || (transcript.facts().vefc_has_frameshift_intron
            && allele_overlaps_any_exon_stretched(transcript, vf_lo, vf_hi));

    for intron in introns {
        let (in_intron_list, in_boundary_list) =
            perl_intron_list_membership(vf_lo, vf_hi, intron.start, intron.end);
        if !in_intron_list && !in_boundary_list {
            continue;
        }
        add_intronic_non_insertion_splice_consequences(
            consequences,
            transcript,
            variant,
            intron.start,
            intron.end,
            suppress_polypyr,
        );
    }
}

/// Supplemental splice check for insertions, over every intron Perl's memoised
/// intron list holds for the anchor pair (`perl_intron_list_membership`), after
/// `within_feature` on the raw anchors (both inside the transcript). The
/// per-endpoint `apply_position` calls cover only the anchor introns.
fn add_splice_for_overlapping_introns_insertion(
    consequences: &mut ConsequenceList,
    variant: &InputVariant,
    transcript: &Transcript,
) {
    let vefc = match transcript.vefc.as_ref() {
        Some(v) => v,
        None => return,
    };
    let introns = if !vefc.introns.is_empty() {
        &vefc.introns
    } else {
        &transcript.introns
    };

    let (v_start, v_end) = (variant.start, variant.end);
    if !overlap_unsorted(v_start, v_end, transcript.start, transcript.end) {
        return;
    }
    let (vf_lo, vf_hi) = (v_start.min(v_end), v_start.max(v_end));

    for intron in introns {
        let (in_intron_list, in_boundary_list) =
            perl_intron_list_membership(vf_lo, vf_hi, intron.start, intron.end);
        if !in_intron_list && !in_boundary_list {
            continue;
        }
        add_intronic_insertion_splice_consequences(
            consequences,
            transcript,
            variant,
            intron.start,
            intron.end,
        );
    }
}

/// For non-insertion variants, add `5_prime_UTR_variant` / `3_prime_UTR_variant`
/// exactly as Perl's `VariationEffect::within_5_prime_utr` /
/// `within_3_prime_utr` do.
///
/// Perl tests `overlap(bvf_start, bvf_end, region_start, region_end)` as
/// `bvf_end >= region_start && bvf_start <= region_end` with no check that the
/// region is non-empty. `_before_coding` uses the region
/// `transcript.start .. coding_region_start - 1` and `_after_coding` uses
/// `coding_region_end + 1 .. transcript.end`; when the CDS reaches the transcript
/// boundary the region is inverted and the test still holds for a span that
/// crosses that boundary point, and fails for one that merely touches it. The
/// `cds_start_NF` / `cds_end_NF` flags play no part.
fn add_utr_for_overlapping_span(
    consequences: &mut ConsequenceList,
    variant: &InputVariant,
    transcript: &Transcript,
) {
    let (crs, cre) = match genomic_coding_bounds(transcript) {
        Some(bounds) => bounds,
        None => return,
    };
    if transcript.translation.is_none() {
        return;
    }

    let lo = variant.start.min(variant.end);
    let hi = variant.start.max(variant.end);

    // An SNV is already classified by its `apply_position` call; only a
    // multi-base variant can span from coding into UTR with neither endpoint
    // landing in the UTR.
    if lo == hi {
        return;
    }

    // `within_cdna`: any exon overlap, not specifically a UTR exon, so a deletion
    // from a CDS exon across an intron into UTR still gets the term.
    let exons = transcript
        .vefc
        .as_ref()
        .map(|v| v.sorted_exons.as_slice())
        .unwrap_or(transcript.exons.as_slice());
    if !exons.iter().any(|e| hi >= e.start && lo <= e.end) {
        return;
    }

    // `overlap(lo, hi, transcript.start, crs - 1)` and
    // `overlap(lo, hi, cre + 1, transcript.end)`.
    let before_coding = hi >= transcript.start && lo < crs;
    let after_coding = hi > cre && lo <= transcript.end;
    let (five, three) = match transcript.strand {
        vep_core::coordinate::Strand::Forward => (before_coding, after_coding),
        vep_core::coordinate::Strand::Reverse => (after_coding, before_coding),
    };
    if five {
        push_unique(consequences, Consequence::FivePrimeUtrVariant);
    }
    if three {
        push_unique(consequences, Consequence::ThreePrimeUtrVariant);
    }
}

/// Genomic coding bounds, preferring the transcript's own fields and falling back to the
/// cDNA-coding offsets in the cache's coordinate mapper.
///
/// `pub(crate)` so the SV predicates in `sv/mod.rs` share it: Perl derives the
/// bounds through its mapper when a cache entry carries the CDS only as cDNA
/// offsets, and reading `coding_region_start` / `_end` alone would withhold
/// `coding_sequence_variant` and both UTR terms on such a transcript.
pub(crate) fn genomic_coding_bounds(transcript: &Transcript) -> Option<(u64, u64)> {
    match (transcript.coding_region_start, transcript.coding_region_end) {
        (Some(start), Some(end)) if start > 0 && end > 0 => return Some((start, end)),
        _ => {}
    }

    let mapper = transcript.vefc.as_ref()?.mapper.as_ref()?;
    if mapper.cdna_coding_start == 0 || mapper.cdna_coding_end < mapper.cdna_coding_start {
        return None;
    }

    let pairs = &mapper.exon_coord_mapper.pairs;
    let start = cdna_to_genomic_pos(mapper.cdna_coding_start, pairs)?;
    let end = cdna_to_genomic_pos(mapper.cdna_coding_end, pairs)?;
    Some((start.min(end), start.max(end)))
}

fn cdna_to_genomic_pos(cdna_pos: u64, pairs: &[vep_core::transcript::MapperPair]) -> Option<u64> {
    for pair in pairs {
        if cdna_pos >= pair.from_start && cdna_pos <= pair.from_end {
            let offset = cdna_pos - pair.from_start;
            let genomic_pos = if pair.ori == 1 {
                pair.to_start + offset
            } else {
                pair.to_end - offset
            };
            return Some(genomic_pos);
        }
    }
    None
}

/// Perl's `_frameshift` flag on an intron: `abs(end - start) <= 12`
/// (`BaseTranscriptVariation::_create_intron_trees`). `_intron_effects` skips
/// such an intron, in both its loops, for every differing region that overlaps
/// it, so the region gets none of that intron's splice or intron terms.
fn is_frameshift_intron(intron_start: u64, intron_end: u64) -> bool {
    intron_end.saturating_sub(intron_start) <= 12
}

/// Perl's `BaseTranscriptVariation::_overlapped_introns` /
/// `_overlapped_introns_boundary` are memoised on the first call, which
/// `_bvfo_preds` makes with the raw variant span, so every later per-region check
/// in `_intron_effects` runs over the introns that span fetched. The interval
/// trees hold `[intron_start-3, intron_end+3]` and the two boundary windows
/// `[intron_start-3, intron_start+7]`, `[intron_end-7, intron_end+3]` (1-based
/// closed, after the half-open insert/fetch adjustment of Set::IntervalTree
/// 0.12, the version Ensembl VEP release 115.2 runs). Returns (in intron list,
/// in boundary list); for an intron shorter than 5bp the boundary windows reach
/// past the intron window, so neither list contains the other.
fn perl_intron_list_membership(
    vf_lo: u64,
    vf_hi: u64,
    intron_start: u64,
    intron_end: u64,
) -> (bool, bool) {
    let overlap_perl = |rs: u64, re: u64, fs: u64, fe: u64| -> bool { re >= fs && rs <= fe };
    let in_intron_list = overlap_perl(
        vf_lo,
        vf_hi,
        intron_start.saturating_sub(3),
        intron_end.saturating_add(3),
    );
    let in_boundary_list = overlap_perl(
        vf_lo,
        vf_hi,
        intron_start.saturating_sub(3),
        intron_start.saturating_add(7),
    ) || overlap_perl(
        vf_lo,
        vf_hi,
        intron_end.saturating_sub(7),
        intron_end.saturating_add(3),
    );
    (in_intron_list, in_boundary_list)
}

/// Splice terms for one intron of a non-insertion variant, following Perl's
/// `BaseTranscriptVariationAllele::_intron_effects`: every window is
/// tested per differing region (`_get_differing_regions`, which for a variant
/// with either allele of 1 base or less is the raw reference span and otherwise
/// the XOR of the alleles, as long as the longer allele), with the
/// non-normalising `overlap`. The intron-loop flags (intron_variant,
/// polypyrimidine) need the intron in the memoised intron list; the
/// boundary-loop flags (donor, acceptor, 5th base, donor region, splice_region,
/// polypyrimidine again) need it in the boundary list.
fn add_intronic_non_insertion_splice_consequences(
    consequences: &mut ConsequenceList,
    transcript: &Transcript,
    variant: &InputVariant,
    intron_start: u64,
    intron_end: u64,
    suppress_polypyr: bool,
) {
    let vf_lo = variant.start.min(variant.end);
    let vf_hi = variant.start.max(variant.end);
    let (in_intron_list, in_boundary_list) =
        perl_intron_list_membership(vf_lo, vf_hi, intron_start, intron_end);
    if !in_intron_list && !in_boundary_list {
        return;
    }

    let overlap_perl = |rs: u64, re: u64, fs: u64, fe: u64| -> bool { re >= fs && rs <= fe };
    // A region overlapping a frameshift intron (span <= 12) skips that intron
    // (the `next` in `_intron_effects`); regions beside it are still tested, so
    // an exonic base adjacent to a 1bp intron can sit in its donor or acceptor
    // window.
    let frameshift_intron = is_frameshift_intron(intron_start, intron_end);
    let mut regions = get_differing_regions_perl_unclamped(variant);
    if frameshift_intron {
        regions.retain(|&mut (rs, re)| !overlap_perl(rs, re, intron_start, intron_end));
    }
    let any_region_overlap =
        |fs: u64, fe: u64| -> bool { regions.iter().any(|&(rs, re)| overlap_perl(rs, re, fs, fe)) };

    let core_start = intron_start.saturating_add(2);
    let core_end = intron_end.saturating_sub(2);
    // The intron loop of `_intron_effects` alone sets `intronic`.
    let within_core_intron = in_intron_list && any_region_overlap(core_start, core_end);
    // Both loops test the polypyrimidine windows.
    let polypyrimidine_splice_site =
        any_region_overlap(intron_end.saturating_sub(16), intron_end.saturating_sub(2));
    let polypyrimidine_splice_site_reverse = any_region_overlap(
        intron_start.saturating_add(2),
        intron_start.saturating_add(16),
    );

    let boundary_flag =
        |fs: u64, fe: u64| -> bool { in_boundary_list && any_region_overlap(fs, fe) };
    let start_splice_site = boundary_flag(intron_start, intron_start.saturating_add(1));
    let end_splice_site = boundary_flag(intron_end.saturating_sub(1), intron_end);
    let fifth_base_splice_site = boundary_flag(
        intron_start.saturating_add(4),
        intron_start.saturating_add(4),
    );
    let donor_region_splice_site = boundary_flag(
        intron_start.saturating_add(2),
        intron_start.saturating_add(5),
    );
    let fifth_base_splice_site_reverse =
        boundary_flag(intron_end.saturating_sub(4), intron_end.saturating_sub(4));
    let donor_region_splice_site_reverse =
        boundary_flag(intron_end.saturating_sub(5), intron_end.saturating_sub(2));
    // `VariationEffect::_intron_overlap`: the two exonic 3bp windows and the two
    // intronic +2..+7 / -7..-2 windows.
    let splice_region =
        boundary_flag(
            intron_start.saturating_sub(3),
            intron_start.saturating_sub(1),
        ) || boundary_flag(intron_end.saturating_add(1), intron_end.saturating_add(3))
            || boundary_flag(
                intron_start.saturating_add(2),
                intron_start.saturating_add(7),
            )
            || boundary_flag(intron_end.saturating_sub(7), intron_end.saturating_sub(2));

    let (donor_splice_site, acceptor_splice_site, donor_5th_base, donor_region, polypyrimidine) =
        match transcript.strand {
            vep_core::coordinate::Strand::Forward => (
                start_splice_site,
                end_splice_site,
                fifth_base_splice_site,
                donor_region_splice_site,
                polypyrimidine_splice_site,
            ),
            vep_core::coordinate::Strand::Reverse => (
                end_splice_site,
                start_splice_site,
                fifth_base_splice_site_reverse,
                donor_region_splice_site_reverse,
                polypyrimidine_splice_site_reverse,
            ),
        };

    if donor_splice_site {
        push_unique(consequences, Consequence::SpliceDonorVariant);
    }
    if acceptor_splice_site {
        push_unique(consequences, Consequence::SpliceAcceptorVariant);
    }
    if donor_5th_base {
        push_unique(consequences, Consequence::SpliceDonor5thBaseVariant);
    }
    let donor_region = donor_region && !donor_5th_base;
    if donor_region {
        push_unique(consequences, Consequence::SpliceDonorRegionVariant);
    }
    // Perl's `include` gate: splice_polypyrimidine_tract_variant requires
    // `exon => 0` and `intron => 1`, so exon overlap suppresses it.
    if polypyrimidine && !suppress_polypyr {
        push_unique(consequences, Consequence::SplicePolypyrimidineTractVariant);
    }
    if within_core_intron {
        push_unique(consequences, Consequence::IntronVariant);
    }
    if splice_region {
        if !donor_splice_site && !acceptor_splice_site && !donor_region && !donor_5th_base {
            push_unique(consequences, Consequence::SpliceRegionVariant);
        } else {
            // Perl's suppression rules also apply to a splice_region the exonic side added.
            remove_consequence(consequences, Consequence::SpliceRegionVariant);
        }
    }
}

/// The bases of reference `perform_shift` (TranscriptVariationAllele.pm:291)
/// reads on either side of the allele: `_genomic_shift` (:427) and the RefSeq
/// edited-sequence shift in `_return_3prime` (:217) take 1,000 of flank.
const SHIFT_SEARCH_AREA: usize = 1000;

/// How far `perform_shift` can move an insertion or deletion of `len` bases: its
/// loop runs to `length(post_seq) - len` inclusive going 3' and to
/// `length(pre_seq) - len + 1` from 1 going 5' (:303-306), so a pattern that fits
/// the flank moves at most `1001 - len` positions; a longer one gets a bound of
/// `length(post_seq)` (:304), 1,000 moves, except that going 5' `len == 1001`
/// leaves the loop bound at zero and nothing moves.
pub fn perform_shift_bound(len: usize, reverse: bool) -> u64 {
    match len {
        0..=SHIFT_SEARCH_AREA => (SHIFT_SEARCH_AREA + 1 - len) as u64,
        _ if reverse && len == SHIFT_SEARCH_AREA + 1 => 0,
        _ => SHIFT_SEARCH_AREA as u64,
    }
}

/// The span of `variant` after the 3' shift along `transcript`'s strand
/// (`_genomic_shift`, TranscriptVariationAllele.pm:411), or `None` when the
/// variant is not an insertion or deletion of plain bases or nothing moves.
pub fn shift_indel_3prime_coords(
    variant: &InputVariant,
    transcript: &Transcript,
    fasta: &vep_fasta::IndexedFasta,
) -> Option<(u64, u64)> {
    shift_indel_3prime(variant, transcript.strand, fasta)
}

/// Whether the record's alleles shift at all: `_return_3prime` and
/// `_genomic_shift` read the class of the whole record
/// (`$vf->var_class`, TranscriptVariationAllele.pm:119-123, 420), the
/// `SO_variation_class` of its allele string (Utils/Sequence.pm:335-398), so a
/// multi-allelic record is an `insertion` only when every allele inserts
/// (`-/A/AA`) and a `deletion` only when every allele deletes, and one that
/// mixes classes (`TT/-/T`) shifts none of its alleles. A record with one
/// allele takes that allele's class, as does an allele of a `--minimal` split,
/// which is a feature of its own, with its own allele string, until the rejoin
/// (InputBuffer.pm:638-718) and is annotated before it (OutputFactory.pm:2289,
/// 2303).
fn record_class_shifts(variant: &InputVariant) -> bool {
    if variant.original_start.is_some() {
        return true;
    }
    let Some(record) = variant.record_allele_string_multi.as_deref() else {
        return true;
    };
    let mut alleles = record.split('/');
    let Some(reference) = alleles.next() else {
        return false;
    };
    let bases = |allele: &str| !allele.is_empty() && allele.bytes().all(|b| b.is_ascii_uppercase());
    if reference == "-" {
        alleles.all(bases)
    } else if bases(reference) {
        alleles.all(|allele| allele == "-")
    } else {
        false
    }
}

/// The bases `perform_shift` walks along the reference, on the forward strand:
/// the inserted allele of an insertion, the deleted allele of a deletion
/// (`_genomic_shift`, TranscriptVariationAllele.pm:420, 445-459, which
/// reverse-complements an allele given on the reverse strand); `None` for every
/// other class and for an allele of a record whose class shifts nothing.
pub fn shift_motif(variant: &InputVariant) -> Option<Vec<u8>> {
    use vep_core::coordinate::Strand;
    use vep_core::variant::VariantClass;

    if !record_class_shifts(variant) {
        return None;
    }
    let bases = match variant.variant_class {
        VariantClass::Insertion => variant.alt_allele(),
        VariantClass::Deletion => variant.ref_allele.as_slice(),
        _ => return None,
    };
    if bases == b"-" || bases.is_empty() {
        return None;
    }
    let motif = bases.to_ascii_uppercase();
    Some(match variant.strand {
        Strand::Forward => motif,
        Strand::Reverse => vep_core::codon::reverse_complement(&motif),
    })
}

/// [`shift_indel_3prime_coords`] towards the 3' end of `strand`: the span an
/// insertion or deletion reaches when every base it passes matches the pattern
/// cyclically, within [`perform_shift_bound`].
pub fn shift_indel_3prime(
    variant: &InputVariant,
    strand: vep_core::coordinate::Strand,
    fasta: &vep_fasta::IndexedFasta,
) -> Option<(u64, u64)> {
    use vep_core::coordinate::Strand;
    use vep_core::variant::VariantClass;

    let (mut start, mut end) = (variant.start, variant.end);
    let mut motif = shift_motif(variant)?;
    let is_insertion = variant.variant_class == VariantClass::Insertion;

    // VEP shifts indels in the transcript 3' direction (genomic right for + strand,
    // genomic left for - strand transcripts).
    let shift_right = strand == Strand::Forward;
    let max_shift = perform_shift_bound(motif.len(), !shift_right);
    let mut shifted = 0u64;

    while shifted < max_shift {
        if shift_right {
            // Right shift by 1:
            // - insertion: compare reference base at `start` (base after insertion site)
            // - deletion: compare reference base at `end + 1` (base after deleted segment)
            let check_pos = if is_insertion {
                start
            } else {
                end.saturating_add(1)
            };

            let Some(ref_base) = fasta.base(&variant.chr, check_pos) else {
                break;
            };
            if ref_base != motif[0] {
                break;
            }

            start = start.saturating_add(1);
            end = end.saturating_add(1);
            motif.rotate_left(1);
            shifted += 1;
        } else {
            // Left shift by 1:
            // - insertion: compare reference base at `end` (base before insertion site)
            // - deletion: compare reference base at `start - 1` (base before deleted segment)
            let check_pos = if is_insertion {
                end
            } else {
                start.saturating_sub(1)
            };
            if check_pos == 0 {
                break;
            }

            let Some(ref_base) = fasta.base(&variant.chr, check_pos) else {
                break;
            };
            if ref_base != *motif.last().unwrap() {
                break;
            }

            start = start.saturating_sub(1);
            end = end.saturating_sub(1);
            motif.rotate_right(1);
            shifted += 1;
        }
    }

    if shifted == 0 {
        None
    } else {
        Some((start, end))
    }
}

/// Add splice_region_variant for exonic positions near splice junctions.
///
/// Mirrors Perl's `VariationEffect::_intron_overlap`, which checks the
/// narrowed variant span (via `_get_differing_regions`) against four windows. Two
/// are the 3bp exonic windows adjacent to each intron:
///   - donor-side:    [intron_start - 3, intron_start - 1]
///   - acceptor-side: [intron_end + 1,   intron_end + 3]
///
/// The other two are nominally intronic:
///   - [intron_start + 2, intron_start + 7]
///   - [intron_end - 7,   intron_end - 2]
///
/// Perl does not clamp those two to the intron, so for an intron shorter than 8bp
/// they spill past the far boundary into the flanking exon and a purely exonic
/// variant 4-7bp away satisfies the predicate. `spill_windows` below reproduces
/// that, which is why this function owns them rather than the intronic path: the
/// variants in question never overlap the intron at all, so no intronic caller is
/// ever reached for them. The test is on the narrowed genomic anchor, not a cDNA
/// distance, which for a large insertion can sit near a junction while the
/// anchor lies outside the 3bp window.
///
/// All of these windows are evaluated inside the boundary loop of
/// `BaseTranscriptVariationAllele::_intron_effects`, which `next`s past a
/// frameshift intron (`is_frameshift_intron`) for every differing region that
/// overlaps it, so a deletion covering a 1bp intron gets no splice term from that
/// intron (the exonic 3bp windows and the spill included) even though it covers
/// them; a region beside the intron is still tested against it.
fn add_exonic_splice_region(
    consequences: &mut ConsequenceList,
    variant: &InputVariant,
    transcript: &Transcript,
) {
    let vefc = match transcript.vefc.as_ref() {
        Some(v) => v,
        None => return,
    };

    let introns = &vefc.introns;
    if introns.is_empty() {
        return;
    }

    // Every differing region is visited, as Perl's _get_differing_regions XORs
    // ref/alt character by character and groups consecutive differences; a complex indel
    // with matching bases inside the pair (ACCCCA to CTTCC shares CC at positions
    // 3-4) yields multiple non-contiguous regions.
    let regions = get_differing_regions(variant);
    // The same regions as Perl maps them, `$vf_start + $region->{e}` with no
    // clamp to the reference span: the frameshift skip reads these, because a
    // net-insertion delins whose region runs past the reference span overlaps
    // an intron the clamped region stops short of, and Perl skips that intron.
    // The window tests keep the clamped regions.
    let perl_regions = get_differing_regions_perl_unclamped(variant);

    // Perl's overlap does not normalise, so start > end is no overlap.
    let overlap_perl = |vs: u64, ve: u64, fs: u64, fe: u64| -> bool { ve >= fs && vs <= fe };

    let skipped_for = |intron: &vep_core::transcript::Intron, (rs, re): (u64, u64)| -> bool {
        is_frameshift_intron(intron.start, intron.end)
            && overlap_perl(rs, re, intron.start, intron.end)
    };

    let is_insertion = matches!(
        variant.variant_class,
        vep_core::variant::VariantClass::Insertion
    );

    // Perl's two nominally-intronic windows, restricted to the part that spills
    // out of a short intron into the flanking exon.
    //
    // For an intron of length L, `[start+2, start+7]` extends `7 - L + 1` bases
    // past `end` when L < 8, and `[end-7, end-2]` symmetrically extends before
    // `start`. Only that spill is returned, so a variant inside the intron is
    // unaffected: the intronic path already owns those.
    //
    // Insertions are excluded: their anchor pair straddles the boundary by
    // construction (`ref = "-"`, `end = start - 1`), and Perl reaches them
    // through `_intron_overlap`'s dedicated `$insertion` clause with its four
    // exact-equality tests rather than through these windows.
    let spill_windows = |intron: &vep_core::transcript::Intron| -> [Option<(u64, u64)>; 2] {
        if is_insertion || intron.end < intron.start {
            return [None, None];
        }
        let intron_len = intron.end - intron.start + 1;
        if intron_len > 7 {
            return [None, None];
        }
        // Variant must lie wholly outside the intron; anything overlapping it is
        // the intronic path's business.
        let donor_spill = {
            let lo = intron.end + 1;
            let hi = intron.start + 7;
            (hi >= lo).then_some((lo, hi))
        };
        let acceptor_spill = {
            let lo = intron.end.saturating_sub(7);
            let hi = intron.start - 1;
            (lo <= hi && intron.end >= 7).then_some((lo, hi))
        };
        [donor_spill, acceptor_spill]
    };

    // A variant lies in a splice region if any of its differing regions overlaps any
    // intron's exonic splice window: a disjunction over (region, intron) pairs.
    // VEP's `BaseTranscriptVariationAllele::_intron_effects` assigns the
    // flag inside its nested loops instead of OR-ing it, so a later region that
    // does not qualify overwrites an earlier one that did; that overwrite is a VEP
    // defect and is not reproduced here.
    let mut splice_region_hit = false;

    for (&(v_start, v_end), &perl_region) in regions.iter().zip(&perl_regions) {
        let mut region_hit = false;
        for intron in introns {
            if skipped_for(intron, perl_region) {
                continue;
            }

            let donor_exonic_start = intron.start.saturating_sub(3);
            let donor_exonic_end = intron.start.saturating_sub(1);
            let donor_exonic_hit = donor_exonic_end > 0
                && overlap_perl(v_start, v_end, donor_exonic_start, donor_exonic_end);

            let acceptor_exonic_start = intron.end + 1;
            let acceptor_exonic_end = intron.end + 3;
            let acceptor_exonic_hit =
                overlap_perl(v_start, v_end, acceptor_exonic_start, acceptor_exonic_end);

            // Short-intron spill only applies when the variant misses the intron
            // itself; otherwise the intronic path classifies it.
            let outside_intron = !overlap_perl(v_start, v_end, intron.start, intron.end);
            let spill_hit = outside_intron
                && spill_windows(intron)
                    .into_iter()
                    .flatten()
                    .any(|(fs, fe)| overlap_perl(v_start, v_end, fs, fe));

            if donor_exonic_hit || acceptor_exonic_hit || spill_hit {
                region_hit = true;
            }
        }

        splice_region_hit = splice_region_hit || region_hit;
    }

    if !splice_region_hit {
        return;
    }

    // Short-intron donor terms take the first hitting region; the first hit is
    // what selects the sub-term.
    let (v_start_first, v_end_first) = exonic_splice_differing_region(variant);

    for intron in introns {
        // This span stands in for every region at once, so the intron is skipped
        // here when any region skips it.
        if perl_regions.iter().any(|&r| skipped_for(intron, r)) {
            continue;
        }

        let donor_exonic_start = intron.start.saturating_sub(3);
        let donor_exonic_end = intron.start.saturating_sub(1);
        let donor_hit = donor_exonic_end > 0
            && overlap_perl(
                v_start_first,
                v_end_first,
                donor_exonic_start,
                donor_exonic_end,
            );

        let acceptor_exonic_start = intron.end + 1;
        let acceptor_exonic_end = intron.end + 3;
        let acceptor_hit = overlap_perl(
            v_start_first,
            v_end_first,
            acceptor_exonic_start,
            acceptor_exonic_end,
        );

        let outside_intron = !overlap_perl(v_start_first, v_end_first, intron.start, intron.end);
        let spill_hit = outside_intron
            && spill_windows(intron)
                .into_iter()
                .flatten()
                .any(|(fs, fe)| overlap_perl(v_start_first, v_end_first, fs, fe));

        if !(donor_hit || acceptor_hit || spill_hit) {
            continue;
        }

        if !add_short_intron_far_boundary_donor_terms(
            consequences,
            transcript,
            v_start_first,
            v_end_first,
            intron.start,
            intron.end,
        ) {
            push_unique(consequences, Consequence::SpliceRegionVariant);
        }
    }
}

fn add_short_intron_far_boundary_donor_terms(
    consequences: &mut ConsequenceList,
    transcript: &Transcript,
    v_start: u64,
    v_end: u64,
    intron_start: u64,
    intron_end: u64,
) -> bool {
    let intron_len = intron_end.saturating_sub(intron_start).saturating_add(1);
    if intron_len > 5 {
        return false;
    }

    let overlap_perl = |fs: u64, fe: u64| -> bool { v_end >= fs && v_start <= fe };

    let (donor_5th_start, donor_5th_end, donor_region_start, donor_region_end) =
        match transcript.strand {
            vep_core::coordinate::Strand::Forward => (
                intron_start.saturating_add(4),
                intron_start.saturating_add(4),
                intron_start.saturating_add(2),
                intron_start.saturating_add(5),
            ),
            vep_core::coordinate::Strand::Reverse => (
                intron_end.saturating_sub(4),
                intron_end.saturating_sub(4),
                intron_end.saturating_sub(5),
                intron_end.saturating_sub(2),
            ),
        };

    // Only treat these as exonic donor sub-terms when the donor window crosses
    // the intron boundary into the opposite exon. Otherwise the normal intronic
    // path should own the classification.
    let donor_5th_crosses_exon = donor_5th_start < intron_start || donor_5th_end > intron_end;
    let donor_region_crosses_exon =
        donor_region_start < intron_start || donor_region_end > intron_end;

    if donor_5th_crosses_exon && overlap_perl(donor_5th_start, donor_5th_end) {
        push_unique(consequences, Consequence::SpliceDonor5thBaseVariant);
        return true;
    }

    if donor_region_crosses_exon && overlap_perl(donor_region_start, donor_region_end) {
        push_unique(consequences, Consequence::SpliceDonorRegionVariant);
        return true;
    }

    false
}

/// The differing regions of one allele pair as genomic `(start, end)` spans.
/// The XOR of two alleles whose longer one has `n` bases groups into at most
/// `n.div_ceil(2)` runs, so four inline slots hold every pair of up to eight
/// bases and a longer delins spills to the heap.
type DifferingRegions = SmallVec<[(u64, u64); 4]>;

/// Compute all differing regions between ref and alt alleles, matching Perl's
/// `VariationFeatureOverlapAllele::_get_differing_regions`.
///
/// Perl XORs ref and alt character-by-character, finds positions where they differ,
/// and groups consecutive differing positions into contiguous regions. This produces
/// multiple regions for complex indels where interior characters match.
///
/// For SNPs or when one allele is 1 base or less, returns a single region covering
/// the full ref span (matching Perl's `else` branch).
///
/// Each region is returned as absolute genomic coordinates (start, end). For pure
/// insertions, returns a single inverted-coordinate region.
fn get_differing_regions(variant: &InputVariant) -> DifferingRegions {
    let lo = variant.start.min(variant.end);
    let hi = variant.start.max(variant.end);

    let ref_allele = variant.ref_allele.as_slice();
    let alt_allele = variant.alt_allele();
    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();

    if ref_is_dash {
        return smallvec![(hi, lo)];
    }
    if alt_is_dash {
        return smallvec![(lo, hi)];
    }

    let ref_len = ref_allele.len();
    let alt_len = alt_allele.len();

    // Perl's _get_differing_regions only does multi-region splitting when both
    // alleles are > 1 base. Otherwise, returns a single region.
    if ref_len <= 1 || alt_len <= 1 {
        return smallvec![(lo, hi)];
    }

    // Perl: `$m = $al1 ^ $al2;` then finds non-zero positions. XOR pads the
    // shorter string with null bytes, so positions beyond it always differ.
    let max_len = ref_len.max(alt_len);
    let mut diff_positions = Vec::new();
    for i in 0..max_len {
        let r = if i < ref_len { ref_allele[i] } else { 0u8 };
        let a = if i < alt_len { alt_allele[i] } else { 0u8 };
        if r ^ a != 0 {
            diff_positions.push(i);
        }
    }

    if diff_positions.is_empty() {
        // Identical alleles.
        return smallvec![(lo, hi)];
    }

    let mut regions = DifferingRegions::new();
    let mut region_start = diff_positions[0];
    let mut prev = diff_positions[0];

    for &pos in &diff_positions[1..] {
        if pos == prev + 1 {
            prev = pos;
        } else {
            // Perl: r_end = vf_start + region.e, offsets relative to the ref allele
            // and clamped to ref_len - 1.
            regions.push((
                lo.saturating_add(region_start as u64),
                lo.saturating_add(prev.min(ref_len - 1) as u64),
            ));
            region_start = pos;
            prev = pos;
        }
    }
    regions.push((
        lo.saturating_add(region_start as u64),
        lo.saturating_add(prev.min(ref_len - 1) as u64),
    ));

    regions
}

/// Compute differing regions without clamping the end offset to `ref_len - 1`.
///
/// This matches Perl's `_get_differing_regions` exactly: the XOR byte string is
/// `max(ref_len, alt_len)` bytes long (shorter string padded with null bytes),
/// so positions beyond the ref span correspond to genomic positions past
/// `ref_end`. Perl's `_intron_effects` maps these positions directly to
/// `vf_start + region->{e}` without clamping, which extends the region into
/// the downstream genomic window, exactly the behaviour needed for
/// splice-site overlap checks on net-insertion delins at exon boundaries.
///
/// Use this variant when the caller needs Perl-identical per-region coordinates
/// for intronic splice sub-term detection. The clamping version
/// (`get_differing_regions`) is the one `add_exonic_splice_region` uses, which
/// relies on inverted regions to signal "no overlap" against the exonic splice
/// windows.
fn get_differing_regions_perl_unclamped(variant: &InputVariant) -> DifferingRegions {
    let lo = variant.start.min(variant.end);
    let hi = variant.start.max(variant.end);

    let ref_allele = variant.ref_allele.as_slice();
    let alt_allele = variant.alt_allele();
    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();

    if ref_is_dash {
        // Pure insertion: inverted coordinates (Perl emits {s: 0, e: -1}).
        return smallvec![(hi, lo)];
    }
    if alt_is_dash {
        return smallvec![(lo, hi)];
    }

    let ref_len = ref_allele.len();
    let alt_len = alt_allele.len();

    // Perl's multi-region splitting only runs when both alleles are > 1 base.
    if ref_len <= 1 || alt_len <= 1 {
        return smallvec![(lo, hi)];
    }

    let max_len = ref_len.max(alt_len);
    let mut diff_positions = Vec::new();
    for i in 0..max_len {
        let r = if i < ref_len { ref_allele[i] } else { 0u8 };
        let a = if i < alt_len { alt_allele[i] } else { 0u8 };
        if r ^ a != 0 {
            diff_positions.push(i);
        }
    }

    if diff_positions.is_empty() {
        return smallvec![(lo, hi)];
    }

    let mut regions = DifferingRegions::new();
    let mut region_start = diff_positions[0];
    let mut prev = diff_positions[0];

    for &pos in &diff_positions[1..] {
        if pos == prev + 1 {
            prev = pos;
        } else {
            // Perl: $r_end = $vf_start + $region->{e}; no clamp.
            regions.push((
                lo.saturating_add(region_start as u64),
                lo.saturating_add(prev as u64),
            ));
            region_start = pos;
            prev = pos;
        }
    }
    regions.push((
        lo.saturating_add(region_start as u64),
        lo.saturating_add(prev as u64),
    ));

    regions
}

/// Compute differing-region coordinates for exonic splice region checks,
/// preserving Perl's inversion semantics for insertions.
///
/// This function returns inverted coordinates (start > end) for pure insertions
/// and for complex indels whose entire ref allele is consumed by shared
/// prefix/suffix. This matches Perl's `_get_differing_regions` output.
///
/// note: This returns a single span (prefix/suffix trimmed). For complex indels,
/// `get_differing_regions` should be used instead to get all non-contiguous segments.
fn exonic_splice_differing_region(variant: &InputVariant) -> (u64, u64) {
    let lo = variant.start.min(variant.end);
    let hi = variant.start.max(variant.end);

    let ref_allele = variant.ref_allele.as_slice();
    let alt_allele = variant.alt_allele();
    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();

    if ref_is_dash {
        // Pure insertion: Perl's _get_differing_regions returns {s: 0, e: -1},
        // which gives r_start = vf_start, r_end = vf_start - 1 (inverted).
        // Use (hi, lo), which is (max, min): inverted so overlap_perl returns false
        // for exonic boundary positions.
        return (hi, lo);
    }
    if alt_is_dash {
        return (lo, hi);
    }

    let min_len = ref_allele.len().min(alt_allele.len());

    let mut prefix = 0usize;
    while prefix < min_len && ref_allele[prefix].eq_ignore_ascii_case(&alt_allele[prefix]) {
        prefix += 1;
    }

    let mut suffix = 0usize;
    while suffix < min_len.saturating_sub(prefix)
        && ref_allele[ref_allele.len() - 1 - suffix]
            .eq_ignore_ascii_case(&alt_allele[alt_allele.len() - 1 - suffix])
    {
        suffix += 1;
    }

    // Not normalised: when prefix + suffix >= ref_allele.len() the result is
    // inverted (start > end), as in Perl.
    (
        lo.saturating_add(prefix as u64),
        hi.saturating_sub(suffix as u64),
    )
}

/// Perl `within_mature_miRNA`: the transcript is biotype `miRNA` and the raw
/// variant span overlaps, by the non-normalising `overlap()`, the genomic
/// projection (`cdna2genomic`) of a `miRNA` attribute's `a-b` cDNA range.
/// Raw coordinates matter: an insertion just past the last mature base has
/// `start > c.end`, so it is outside, and only the whole span decides, so a
/// deletion whose endpoints both lie outside the mature range but cover it is
/// inside.
fn variant_overlaps_mature_mirna(variant: &InputVariant, transcript: &Transcript) -> bool {
    if !transcript.facts().is_mirna {
        return false;
    }
    let Some(pairs) = transcript
        .vefc
        .as_ref()
        .and_then(|v| v.mapper.as_ref())
        .map(|m| m.exon_coord_mapper.pairs.as_slice())
    else {
        return false;
    };
    transcript
        .attributes
        .iter()
        .filter(|a| a.code == "miRNA")
        .filter_map(|a| mirna_attribute_cdna_range(&a.value))
        .any(|(cdna_lo, cdna_hi)| {
            pairs.iter().any(|pair| {
                let lo = cdna_lo.max(pair.from_start);
                let hi = cdna_hi.min(pair.from_end);
                if lo > hi {
                    return false;
                }
                let (g_lo, g_hi) = if pair.ori == 1 {
                    (
                        pair.to_start + (lo - pair.from_start),
                        pair.to_start + (hi - pair.from_start),
                    )
                } else {
                    (
                        pair.to_end - (hi - pair.from_start),
                        pair.to_end - (lo - pair.from_start),
                    )
                };
                variant.end >= g_lo && variant.start <= g_hi
            })
        })
}

/// The first `(\d+)-(\d+)` in a `miRNA` attribute value, which is all Perl's
/// regex captures of it.
fn mirna_attribute_cdna_range(value: &str) -> Option<(u64, u64)> {
    let bytes = value.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        while i < bytes.len() && !bytes[i].is_ascii_digit() {
            i += 1;
        }
        let left_start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if left_start == i || i >= bytes.len() || bytes[i] != b'-' {
            continue;
        }
        let right_start = i + 1;
        i = right_start;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if right_start == i {
            continue;
        }
        let left = value[left_start..right_start - 1].parse::<u64>().ok()?;
        let right = value[right_start..i].parse::<u64>().ok()?;
        return Some((left.min(right), left.max(right)));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::make_test_transcript;
    use vep_core::consequence::FeatureType;

    /// A one-exon 85-base miRNA transcript at 1000-1084 carrying the `ncRNA`
    /// structure of ENST00000612014 in the Ensembl VEP 116 GRCh38 cache,
    /// `1:85 .7(16.(4.(3.2(4.10)4.2)3.).)14.6)5`: cDNA positions 1-7 are loop,
    /// 8-23 stem, 24 loop, 25-28 stem and so on. On the forward strand cDNA
    /// position `c` is genomic `999 + c`; on the reverse strand it is `1085 - c`.
    fn mirna_transcript(strand: vep_core::coordinate::Strand, structure: &str) -> Transcript {
        use vep_core::transcript::{
            Attribute, Exon, ExonCoordMapper, MapperPair, TranscriptMapper,
        };
        let mut tx = make_test_transcript();
        tx.stable_id = "ENST00000612014".into();
        tx.biotype = "miRNA".into();
        tx.strand = strand;
        tx.start = 1000;
        tx.end = 1084;
        tx.exons = vec![Exon {
            stable_id: None,
            start: 1000,
            end: 1084,
            rank: 1,
            phase: -1,
            end_phase: -1,
        }];
        tx.introns = Vec::new();
        tx.translation = None;
        tx.protein_id = None;
        tx.cdna_coding_start = None;
        tx.cdna_coding_end = None;
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation_start = None;
        tx.translation_end = None;
        tx.attributes = vec![Attribute {
            code: "ncRNA".to_string(),
            value: structure.to_string(),
        }];
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.translateable_seq = None;
        vefc.introns = Vec::new();
        vefc.sorted_exons = tx.exons.clone();
        vefc.mapper = Some(TranscriptMapper {
            start_phase: -1,
            cdna_coding_start: 0,
            cdna_coding_end: 0,
            exon_coord_mapper: ExonCoordMapper::new(vec![MapperPair {
                from_start: 1,
                from_end: 85,
                to_start: 1000,
                to_end: 1084,
                ori: strand.as_i8(),
            }]),
        });
        tx
    }

    const HAIRPIN: &str = "1:85 .7(16.(4.(3.2(4.10)4.2)3.).)14.6)5";

    fn snv(pos: u64) -> InputVariant {
        InputVariant::new("21".into(), pos, pos, b"A".to_vec(), b"G".to_vec())
    }

    fn terms(variant: &InputVariant, tx: &Transcript) -> Vec<String> {
        mirna_structure_terms(variant, tx)
    }

    /// `miRNA` reads the structure element under each cDNA position of the
    /// variant (OutputFactory.pm:1624-1664 at release/116.2): a base in a run
    /// of dots is `miRNA_loop`, a base in a bracket run `miRNA_stem`, a
    /// deletion across the boundary both, sorted (Ensembl VEP 116.2 writes
    /// `miRNA=miRNA_loop,miRNA_stem` for the deletion `21 8249508 TGCG T` on
    /// MIR6724-2 and `miRNA=miRNA_stem` for the insertion `21 8249510 C CTTA`,
    /// whose two flanking bases are both stem). A deletion from one arm across
    /// the terminal loop into the other arm (cDNA 35-52, `(4.10)4`) lists
    /// `miRNA_stem` once per bracket kind, as the Perl `sort map` over the
    /// distinct characters does. The same cDNA positions decide on the reverse
    /// strand; a variant outside the transcript, a transcript of another
    /// biotype, and an attribute without a structure give nothing.
    #[test]
    fn mirna_stem_and_loop_from_structure() {
        use vep_core::coordinate::Strand;
        let forward = mirna_transcript(Strand::Forward, HAIRPIN);
        assert_eq!(terms(&snv(1003), &forward), ["miRNA_loop"]);
        assert_eq!(terms(&snv(1010), &forward), ["miRNA_stem"]);
        let deletion = InputVariant::new("21".into(), 1006, 1008, b"ACG".to_vec(), b"-".to_vec());
        assert_eq!(terms(&deletion, &forward), ["miRNA_loop", "miRNA_stem"]);
        let insertion = InputVariant::new("21".into(), 1010, 1009, b"-".to_vec(), b"TTA".to_vec());
        assert_eq!(terms(&insertion, &forward), ["miRNA_stem"]);
        let arm_to_arm = InputVariant::new(
            "21".into(),
            1034,
            1051,
            b"ACGTACGTACGTACGTAC".to_vec(),
            b"-".to_vec(),
        );
        assert_eq!(
            terms(&arm_to_arm, &forward),
            ["miRNA_loop", "miRNA_stem", "miRNA_stem"]
        );
        assert!(terms(&snv(990), &forward).is_empty());
        assert!(terms(&snv(1090), &forward).is_empty());

        let reverse = mirna_transcript(Strand::Reverse, HAIRPIN);
        assert_eq!(terms(&snv(1084), &reverse), ["miRNA_loop"]);
        assert_eq!(terms(&snv(1075), &reverse), ["miRNA_stem"]);

        let mut snorna = mirna_transcript(Strand::Forward, HAIRPIN);
        snorna.biotype = "snoRNA".into();
        assert!(terms(&snv(1003), &snorna).is_empty());
        let bare = mirna_transcript(Strand::Forward, "1:85");
        assert!(terms(&snv(1003), &bare).is_empty());
        let no_structure = mirna_transcript(Strand::Forward, "1:85 x");
        assert!(terms(&snv(1003), &no_structure).is_empty());

        let config = EffectsConfig::default();
        let tc = calculate_consequences(&deletion, &forward, &config).unwrap();
        assert_eq!(tc.mirna, ["miRNA_loop", "miRNA_stem"]);
        assert!(
            calculate_consequences(&snv(1010), &make_test_transcript(), &config)
                .is_none_or(|tc| tc.mirna.is_empty())
        );
    }

    /// Ensembl VEP's index test skips an offset below zero or above the
    /// structure's length, so the element one past a structure shorter than
    /// its `<start>..<end>` span is read as undefined, which the list keeps as
    /// an empty first term; a variant outside the span contributes nothing.
    #[test]
    fn mirna_structure_read_past_its_end_keeps_the_undefined_term() {
        use vep_core::coordinate::Strand;
        let short = mirna_transcript(Strand::Forward, "1:85 .7(16.(4.(3.2(4.10)4.2)3.).)14.6)4");
        assert_eq!(terms(&snv(1084), &short), [""]);
        let span = InputVariant::new("21".into(), 1083, 1084, b"AC".to_vec(), b"-".to_vec());
        assert_eq!(terms(&span, &short), ["", "miRNA_stem"]);
        let offset = mirna_transcript(Strand::Forward, "10:85 .7(16.(4.(3.2(4.10)4.2)3.).)14.6)5");
        assert!(terms(&snv(1003), &offset).is_empty());
        assert_eq!(terms(&snv(1012), &offset), ["miRNA_loop"]);
    }

    /// Under `--domains` a coding row lists the protein features whose residue
    /// span overlaps the variant's, in the cache's order, as `label:identifier`
    /// with whitespace, `;` and `=` written `_` (OutputFactory.pm:1494-1512 at
    /// release/116.2 over `get_overlapping_ProteinFeatures`, whose test is
    /// Ensembl's `overlap` on the raw residue span). On the test transcript
    /// residue n is genomic 25_000_050 + 3(n - 1) to + 2: the residue-3 SNV
    /// overlaps the Pfam and the AFDB features and not the residue-5 Gene3D one;
    /// the deletion of residues 4 and 5 overlaps all three; the insertion
    /// between residues 4 and 5 maps to the span 5-4, which Ensembl VEP
    /// overlaps only with a feature covering both residues, so the Gene3D and
    /// AFDB features drop out; a 5' UTR variant is not coding and lists none;
    /// without the flag the list stays empty.
    #[test]
    fn domains_overlapping_protein_position() {
        use vep_core::transcript::ProteinFeature;
        let mut tx = make_test_transcript();
        tx.vefc.as_mut().unwrap().protein_features = vec![
            ProteinFeature {
                start: 1,
                end: 10,
                hseqname: "PF00080".to_string(),
                analysis: Some("Pfam".to_string()),
            },
            ProteinFeature {
                start: 5,
                end: 5,
                hseqname: "3.30.160.20".to_string(),
                analysis: Some("Gene3D".to_string()),
            },
            ProteinFeature {
                start: 2,
                end: 4,
                hseqname: "AF-P60331-F1".to_string(),
                analysis: Some("AFDB-ENSP mappings".to_string()),
            },
            ProteinFeature {
                start: 20,
                end: 30,
                hseqname: "seg;a=b".to_string(),
                analysis: Some("Low complexity (Seg)".to_string()),
            },
        ];
        let config = EffectsConfig {
            domains: true,
            ..EffectsConfig::default()
        };
        let pairs = |variant: &InputVariant, config: &EffectsConfig| {
            calculate_consequences(variant, &tx, config)
                .unwrap()
                .domains
                .iter()
                .map(|(label, id)| format!("{label}:{id}"))
                .collect::<Vec<_>>()
        };
        let snv_residue_3 = InputVariant::new(
            "21".into(),
            25_000_057,
            25_000_057,
            b"G".to_vec(),
            b"T".to_vec(),
        );
        assert_eq!(
            pairs(&snv_residue_3, &config),
            ["Pfam:PF00080", "AFDB-ENSP_mappings:AF-P60331-F1"]
        );
        let deletion_residues_4_5 = InputVariant::new(
            "21".into(),
            25_000_059,
            25_000_064,
            b"AAATTC".to_vec(),
            b"-".to_vec(),
        );
        assert_eq!(
            pairs(&deletion_residues_4_5, &config),
            [
                "Pfam:PF00080",
                "Gene3D:3.30.160.20",
                "AFDB-ENSP_mappings:AF-P60331-F1"
            ]
        );
        let insertion_between_4_and_5 = InputVariant::new(
            "21".into(),
            25_000_062,
            25_000_061,
            b"-".to_vec(),
            b"GGG".to_vec(),
        );
        assert_eq!(pairs(&insertion_between_4_and_5, &config), ["Pfam:PF00080"]);
        let snv_residue_25 = InputVariant::new(
            "21".into(),
            25_000_123,
            25_000_123,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        assert_eq!(
            pairs(&snv_residue_25, &config),
            ["Low_complexity_(Seg):seg_a_b"]
        );
        let utr = InputVariant::new(
            "21".into(),
            25_000_010,
            25_000_010,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        assert!(pairs(&utr, &config).is_empty());
        assert!(pairs(&snv_residue_3, &EffectsConfig::default()).is_empty());
    }

    /// `split /\s+|\:/`: a colon or a whitespace run ends a field, a colon
    /// beside whitespace leaves an empty field, and only the first three fields
    /// are read.
    #[test]
    fn ncrna_fields_split_on_colon_or_whitespace_run() {
        assert_eq!(ncrna_fields("1:85 .7(16"), ["1", "85", ".7(16"]);
        assert_eq!(ncrna_fields("1:85  .7(16 more"), ["1", "85", ".7(16"]);
        assert_eq!(ncrna_fields("1: 85 x"), ["1", "", "85"]);
        assert_eq!(ncrna_fields(" 1:85 x"), ["", "1", "85"]);
        assert_eq!(ncrna_fields("1:85"), ["1", "85"]);
        assert_eq!(ncrna_fields("1:85 "), ["1", "85"]);
    }

    #[test]
    fn test_upstream_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            24_997_000,
            24_997_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::UpstreamGeneVariant));
        assert_eq!(tc.distance, Some(3000));
    }

    #[test]
    fn test_downstream_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_008_000,
            25_008_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::DownstreamGeneVariant));
        assert_eq!(tc.distance, Some(2000));
    }

    #[test]
    fn distance_is_the_nearest_gap_between_variant_and_transcript_ends() {
        // Transcript 25,000,000-25,006,000 on the forward strand. A 2-bp deletion
        // ending one base short of the transcript start is 1 bp away by its end and
        // 2 bp away by its start; VEP prints 1. A deletion downstream measures from
        // its start, the end being farther.
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let upstream_del = InputVariant::new(
            "21".into(),
            24_999_998,
            24_999_999,
            b"AC".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&upstream_del, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::UpstreamGeneVariant));
        assert_eq!(tc.distance, Some(1));

        let downstream_del = InputVariant::new(
            "21".into(),
            25_006_010,
            25_006_024,
            b"ACGTACGTACGTACG".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&downstream_del, &tx, &config).unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::DownstreamGeneVariant));
        assert_eq!(tc.distance, Some(10));
    }

    #[test]
    fn test_five_prime_utr_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // cDNA 11 -> 5' UTR (cdna_coding_start = 51)
        let variant = InputVariant::new(
            "21".into(),
            25_000_010,
            25_000_010,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::FivePrimeUtrVariant));
    }

    #[test]
    fn test_three_prime_utr_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // cDNA 901 -> 3' UTR (cdna_coding_end = 900)
        let variant = InputVariant::new(
            "21".into(),
            25_004_300,
            25_004_300,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::ThreePrimeUtrVariant));
    }

    #[test]
    fn test_intron_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_001_000,
            25_001_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_spanning_exon_intron_deletion_adds_intronic_consequences() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_295,
            25_000_305,
            b"ACGTACGTACG".to_vec(),
            b"-".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(tc.consequences.iter().any(|c| {
            matches!(
                c,
                Consequence::SpliceDonorVariant
                    | Consequence::SpliceAcceptorVariant
                    | Consequence::SpliceDonor5thBaseVariant
                    | Consequence::SpliceDonorRegionVariant
                    | Consequence::SpliceRegionVariant
            )
        }));
    }

    #[test]
    fn test_splice_donor_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // First base of intron 1 (dist_to_donor = 0)
        // Intron 1 starts at 25_000_300
        let variant = InputVariant::new(
            "21".into(),
            25_000_300,
            25_000_300,
            b"G".to_vec(),
            b"T".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::SpliceDonorVariant));
        assert!(!tc.consequences.contains(&Consequence::IntronVariant));
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_splice_acceptor_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // Last base of intron 1 (dist_to_acceptor = 0)
        // Intron 1 ends at 25_001_999
        let variant = InputVariant::new(
            "21".into(),
            25_001_999,
            25_001_999,
            b"G".to_vec(),
            b"T".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::SpliceAcceptorVariant));
        assert!(!tc.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_missense_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // CDS pos 5 = codon 2, pos 1 = GCT -> GAT = Ala -> Asp
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::MissenseVariant));
        assert_eq!(tc.amino_acids.as_deref(), Some("A/D"));
        assert_eq!(tc.impact, Impact::MODERATE);
    }

    #[test]
    fn test_ambiguous_alt_allele_n_gets_coding_sequence_variant() {
        // Perl's `TranscriptVariationAllele::peptide` returns undef for a non-ACGT
        // allele (`seq_is_unambiguous_dna` fails), so only coding_unknown fires:
        // coding_sequence_variant, not missense.
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"N".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config)
            .expect("ambiguous ALT=N should still produce a consequence");
        assert!(
            tc.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "Ambiguous ALT=N should get coding_sequence_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences.contains(&Consequence::MissenseVariant),
            "Ambiguous ALT=N should not get missense_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_synonymous_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // CDS pos 6 = codon 2, pos 2 = GCT -> GCC = Ala -> Ala
        let variant = InputVariant::new(
            "21".into(),
            25_000_055,
            25_000_055,
            b"T".to_vec(),
            b"C".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::SynonymousVariant));
        assert_eq!(tc.impact, Impact::LOW);
    }

    #[test]
    fn test_stop_gained() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // CDS pos 7 = codon 3, pos 0 = GGA -> TGA = Gly -> Stop
        let variant = InputVariant::new(
            "21".into(),
            25_000_056,
            25_000_056,
            b"G".to_vec(),
            b"T".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::StopGained));
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_frameshift_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // 2bp insertion at CDS pos 3 -> frameshift
        let variant = InputVariant::new(
            "21".into(),
            25_000_052,
            25_000_052,
            b"-".to_vec(),
            b"AA".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::FrameshiftVariant));
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_complex_inframe_indel_calls_protein_altering_variant() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // Transcript-matching complex inframe indel where neither allele is '-'.
        // Codon 3 in the test CDS is GGA. Replacing it with 9bp `TTTGGAGCT`
        // yields a net +6bp inframe event whose altered peptide does not start
        // or end with the ref peptide, so Perl falls through to
        // protein_altering_variant.
        let variant = InputVariant::new(
            "21".into(),
            25_000_056,
            25_000_058,
            b"GGA".to_vec(),
            b"TTTGGAGCT".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::ProteinAlteringVariant));
    }

    #[test]
    fn test_no_consequence_too_far() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            24_000_000,
            24_000_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_none());
    }

    #[test]
    fn test_non_coding_transcript_variant_added_for_non_coding_intron() {
        let mut tx = make_test_transcript();
        tx.biotype = "processed_transcript".into();
        tx.translation = None;
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(mapper) = vefc.mapper.as_mut() {
                mapper.cdna_coding_start = 0;
                mapper.cdna_coding_end = 0;
            }
        }

        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_001_000,
            25_001_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptVariant));
        assert!(!tc.consequences.contains(&Consequence::NmdTranscriptVariant));
    }

    #[test]
    fn test_non_coding_transcript_variant_not_added_for_non_coding_exon() {
        let mut tx = make_test_transcript();
        tx.biotype = "processed_transcript".into();
        tx.translation = None;
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(mapper) = vefc.mapper.as_mut() {
                mapper.cdna_coding_start = 0;
                mapper.cdna_coding_end = 0;
            }
        }

        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_010,
            25_000_010,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptVariant));
    }

    #[test]
    fn test_nmd_transcript_variant_added_for_nmd_intron() {
        let mut tx = make_test_transcript();
        tx.biotype = "nonsense_mediated_decay".into();
        tx.translation = None;
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(mapper) = vefc.mapper.as_mut() {
                mapper.cdna_coding_start = 0;
                mapper.cdna_coding_end = 0;
            }
        }

        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_001_000,
            25_001_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(tc.consequences.contains(&Consequence::NmdTranscriptVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptVariant));
    }

    #[test]
    fn test_non_coding_context_not_added_for_upstream_only_hit() {
        let mut tx = make_test_transcript();
        tx.biotype = "processed_transcript".into();
        tx.translation = None;
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            24_997_000,
            24_997_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::UpstreamGeneVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptVariant));
        assert!(!tc.consequences.contains(&Consequence::NmdTranscriptVariant));
    }

    #[test]
    fn test_non_coding_exon_near_junction_gets_splice_region() {
        let mut tx = make_test_transcript();
        tx.biotype = "processed_transcript".into();
        tx.translation = None;
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(mapper) = vefc.mapper.as_mut() {
                mapper.cdna_coding_start = 0;
                mapper.cdna_coding_end = 0;
            }
        }
        let config = EffectsConfig::default();

        let variant = InputVariant::new(
            "21".into(),
            25_000_299,
            25_000_299,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant));
        assert!(tc.consequences.contains(&Consequence::SpliceRegionVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptVariant));
    }

    #[test]
    fn test_splice_donor_region_upper_bound_is_five() {
        let mut csq: ConsequenceList = smallvec![Consequence::IntronVariant];
        add_splice_consequences(&mut csq, 5, 100, false, false);
        assert!(csq.contains(&Consequence::SpliceDonorRegionVariant));
        assert!(!csq.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn test_splice_donor_region_not_called_at_six() {
        let mut csq: ConsequenceList = smallvec![Consequence::IntronVariant];
        add_splice_consequences(&mut csq, 6, 100, false, false);
        assert!(!csq.contains(&Consequence::SpliceDonorRegionVariant));
        assert!(csq.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn test_splice_donor_5th_base_suppresses_donor_region_and_splice_region() {
        let mut csq: ConsequenceList = smallvec![Consequence::IntronVariant];
        add_splice_consequences(&mut csq, 4, 100, false, false);
        assert!(csq.contains(&Consequence::SpliceDonor5thBaseVariant));
        assert!(!csq.contains(&Consequence::SpliceDonorRegionVariant));
        assert!(!csq.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn test_polypyrimidine_window_extends_to_sixteen() {
        let mut csq: ConsequenceList = smallvec![Consequence::IntronVariant];
        add_splice_consequences(&mut csq, 100, 16, false, false);
        assert!(csq.contains(&Consequence::SplicePolypyrimidineTractVariant));
    }

    #[test]
    fn test_polypyrimidine_not_called_at_seventeen() {
        let mut csq: ConsequenceList = smallvec![Consequence::IntronVariant];
        add_splice_consequences(&mut csq, 100, 17, false, false);
        assert!(!csq.contains(&Consequence::SplicePolypyrimidineTractVariant));
    }

    #[test]
    fn test_insertion_at_exon_intron_boundary_is_not_donor() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_300,
            25_000_299,
            b"-".to_vec(),
            b"A".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::SpliceRegionVariant));
        assert!(!tc.consequences.contains(&Consequence::SpliceDonorVariant));
    }

    #[test]
    fn test_insertion_three_bases_into_exon_does_not_get_splice_region() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // Exon 2 starts at 25_002_000. Insert between bases 3 and 4 of the exon
        // (between 25_002_002 and 25_002_003). Perl represents pure insertions
        // with an inverted differing region, so the raw overlap check does not
        // count this as overlapping the acceptor-side exonic 3bp window.
        let variant = InputVariant::new(
            "21".into(),
            25_002_003,
            25_002_002,
            b"-".to_vec(),
            b"A".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::FrameshiftVariant));
        assert!(!tc.consequences.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn test_short_forward_intron_can_emit_donor_region_from_far_exon() {
        let mut tx = make_test_transcript();
        if let Some(vefc) = tx.vefc.as_mut() {
            vefc.introns[0].start = 25_000_312;
            vefc.introns[0].end = 25_000_313;
        }

        let variant = InputVariant::new(
            "21".into(),
            25_000_314,
            25_000_314,
            b"A".to_vec(),
            b"G".to_vec(),
        );

        let mut csq: ConsequenceList = SmallVec::new();
        add_exonic_splice_region(&mut csq, &variant, &tx);

        assert!(csq.contains(&Consequence::SpliceDonorRegionVariant));
        assert!(!csq.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn test_short_reverse_intron_can_emit_donor_region_from_far_exon() {
        let mut tx = make_test_transcript();
        tx.strand = vep_core::coordinate::Strand::Reverse;
        if let Some(vefc) = tx.vefc.as_mut() {
            vefc.introns[0].start = 25_000_312;
            vefc.introns[0].end = 25_000_313;
        }

        let variant = InputVariant::new(
            "21".into(),
            25_000_310,
            25_000_310,
            b"A".to_vec(),
            b"G".to_vec(),
        );

        let mut csq: ConsequenceList = SmallVec::new();
        add_exonic_splice_region(&mut csq, &variant, &tx);

        assert!(csq.contains(&Consequence::SpliceDonorRegionVariant));
        assert!(!csq.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn test_get_differing_regions_complex_indel_multi_segment() {
        // ACCCCA → CTTCC: positions 0-2 differ, 3-4 match, 5 differs (net deletion)
        // Perl's XOR produces two regions: (0,2) and (5,5).
        let variant = InputVariant::new(
            "10".into(),
            13112454,
            13112459,
            b"ACCCCA".to_vec(),
            b"CTTCC".to_vec(),
        );
        let regions = get_differing_regions(&variant);
        assert_eq!(regions.len(), 2, "should produce 2 differing regions");
        assert_eq!(regions[0], (13112454, 13112456)); // positions 0-2
        assert_eq!(regions[1], (13112459, 13112459)); // position 5
    }

    #[test]
    fn test_get_differing_regions_snp_single_region() {
        let variant = InputVariant::new("10".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let regions = get_differing_regions(&variant);
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0], (100, 100));
    }

    #[test]
    fn test_get_differing_regions_pure_insertion() {
        let variant = InputVariant::new("10".into(), 101, 100, b"-".to_vec(), b"ACG".to_vec());
        let regions = get_differing_regions(&variant);
        assert_eq!(regions.len(), 1);
        // Pure insertion: inverted coordinates (hi, lo) = (101, 100)
        assert_eq!(regions[0], (101, 100));
    }

    /// ACG to ATGTT: the XOR differs at offsets 1, 3 and 4, so the second region
    /// starts past the 3-base reference span. The clamped form ends it at
    /// `ref_len - 1` (an inverted span, no overlap), the Perl form at offset 4.
    #[test]
    fn test_get_differing_regions_clamps_only_the_exonic_form() {
        let variant = InputVariant::new("10".into(), 500, 502, b"ACG".to_vec(), b"ATGTT".to_vec());
        let clamped = get_differing_regions(&variant);
        assert_eq!(clamped.as_slice(), &[(501, 501), (503, 502)]);
        let unclamped = get_differing_regions_perl_unclamped(&variant);
        assert_eq!(unclamped.as_slice(), &[(501, 501), (503, 504)]);
    }

    /// Ten alternating differences give five regions, one past the inline
    /// capacity, in XOR order.
    #[test]
    fn test_get_differing_regions_past_the_inline_capacity() {
        let variant = InputVariant::new(
            "10".into(),
            200,
            209,
            b"ACACACACAC".to_vec(),
            b"GCGCGCGCGC".to_vec(),
        );
        let expected = [(200, 200), (202, 202), (204, 204), (206, 206), (208, 208)];
        let regions = get_differing_regions(&variant);
        assert_eq!(regions.as_slice(), &expected);
        assert!(regions.spilled());
        let unclamped = get_differing_regions_perl_unclamped(&variant);
        assert_eq!(unclamped.as_slice(), &expected);
    }

    /// The biotype facts agree with the string predicates they stand for on every
    /// biotype of the golden caches (`tests/golden/116/*/json_cache`) and on the
    /// case variants only the NMD test admits.
    #[test]
    fn facts_biotype_fields_agree_with_the_string_predicates() {
        const GOLDEN_BIOTYPES: [&str; 32] = [
            "protein_coding",
            "lncRNA",
            "nonsense_mediated_decay",
            "retained_intron",
            "processed_transcript",
            "protein_coding_CDS_not_defined",
            "antisense",
            "lincRNA",
            "processed_pseudogene",
            "miRNA",
            "snoRNA",
            "misc_RNA",
            "snRNA",
            "Mt_tRNA",
            "unprocessed_pseudogene",
            "sense_intronic",
            "transcribed_processed_pseudogene",
            "TEC",
            "transcribed_unprocessed_pseudogene",
            "sense_overlapping",
            "pseudogene",
            "rRNA",
            "non_stop_decay",
            "rRNA_pseudogene",
            "Mt_rRNA",
            "transcribed_unitary_pseudogene",
            "ribozyme",
            "polymorphic_pseudogene",
            "LRG_gene",
            "IG_J_gene",
            "IG_C_gene",
            "IG_V_gene",
        ];
        const CASE_VARIANTS: [&str; 6] = [
            "Nonsense_Mediated_Decay",
            "NONSENSE_MEDIATED_DECAY",
            "Protein_Coding",
            "MIRNA",
            "mirna",
            "",
        ];
        for biotype in GOLDEN_BIOTYPES.iter().chain(CASE_VARIANTS.iter()) {
            let mut tx = make_test_transcript();
            tx.biotype = (*biotype).into();
            let facts = *tx.facts();
            assert_eq!(facts.is_mirna, &*tx.biotype == "miRNA", "{biotype}");
            assert_eq!(facts.is_protein_coding, tx.is_protein_coding(), "{biotype}");
            assert_eq!(facts.is_nmd, tx.is_nmd_transcript(), "{biotype}");
            assert_eq!(
                facts.is_non_coding_context,
                tx.is_non_coding_context_transcript(),
                "{biotype}"
            );
        }
    }

    #[test]
    fn test_exonic_splice_region_accumulates_across_regions() {
        // ACCCCA to CTTCC at 13112454-13112459 with intron 3 ending at 13112452: the
        // acceptor exonic window is 13112453-13112455, the first differing region
        // (13112454-13112456) overlaps it and the last (13112459) does not. VEP's
        // `_intron_effects` assigns the flag per region, so its last region
        // overwrites the first and it drops the term; the disjunction here emits it.
        let mut tx = make_test_transcript();
        if let Some(vefc) = tx.vefc.as_mut() {
            // Set intron to end at 13112452 (acceptor at 13112453-13112455)
            vefc.introns[0].start = 13110477;
            vefc.introns[0].end = 13112452;
        }

        let variant = InputVariant::new(
            "10".into(),
            13112454,
            13112459,
            b"ACCCCA".to_vec(),
            b"CTTCC".to_vec(),
        );

        let mut csq: ConsequenceList = SmallVec::new();
        add_exonic_splice_region(&mut csq, &variant, &tx);

        assert!(
            csq.contains(&Consequence::SpliceRegionVariant),
            "a multi-region indel whose first region overlaps the exonic splice window lies in a splice region; a later non-overlapping region must not erase that"
        );
    }

    #[test]
    fn test_frameshift_insertion_can_create_immediate_stop_gained() {
        // Perl output `tat/tAat` produces `Y/*` and adds `stop_gained` alongside
        // `frameshift_variant`.
        let mut tx = make_test_transcript();
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(seq) = vefc.translateable_seq.as_mut() {
                // Make the first codon TAT (Tyr). Then insert 'A' after the first base
                // (CDS pos 2) to get TAA (Stop) as the first altered codon.
                seq.replace_range(0..3, "TAT");
            }
        }

        let config = EffectsConfig::default();
        // CDS starts at cDNA pos 51, which is genomic 25_000_050 in the test transcript.
        // Insert between CDS bases 1 and 2 (between genomic 25_000_050 and 25_000_051).
        let variant = InputVariant::new(
            "21".into(),
            25_000_051,
            25_000_050,
            b"-".to_vec(),
            b"A".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::FrameshiftVariant));
        assert!(tc.consequences.contains(&Consequence::StopGained));
    }

    #[test]
    fn test_frameshift_insertion_stop_in_second_codon_is_stop_gained() {
        // An insertion can create a stop in the *second* codon of Perl's local
        // window (AGC -> A CTTA GC => ACT TAG C => T*X), and Perl reports stop_gained.
        let mut tx = make_test_transcript();
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(seq) = vefc.translateable_seq.as_mut() {
                seq.replace_range(0..3, "AGC"); // Ser (S)
            }
        }

        let config = EffectsConfig::default();
        // Insert CTTA between CDS bases 1 and 2 (between genomic 25_000_050 and 25_000_051).
        let variant = InputVariant::new(
            "21".into(),
            25_000_051,
            25_000_050,
            b"-".to_vec(),
            b"CTTA".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::FrameshiftVariant));
        assert!(tc.consequences.contains(&Consequence::StopGained));
        assert_eq!(
            tc.consequences
                .iter()
                .filter(|c| **c == Consequence::StopGained)
                .count(),
            1
        );
    }

    #[test]
    fn test_frameshift_deletion_partial_codon_does_not_overcall_stop_gained() {
        // A 1-base deletion shortens the codon window to a partial codon (a VX
        // peptide); Perl does not call stop_gained because the partial codon
        // translates to X, not *.
        let mut tx = make_test_transcript();
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(seq) = vefc.translateable_seq.as_mut() {
                // First codon GAG (E). Delete 1 base at CDS pos 2 -> "G_G..." = partial.
                seq.replace_range(0..3, "GAG");
            }
        }

        let config = EffectsConfig::default();
        // Delete the 2nd CDS base (genomic 25_000_051).
        let variant = InputVariant::new(
            "21".into(),
            25_000_051,
            25_000_051,
            b"A".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::FrameshiftVariant));
        assert!(
            !tc.consequences.contains(&Consequence::StopGained),
            "1-base deletion producing partial codon should not call stop_gained"
        );
    }

    #[test]
    fn test_reverse_strand_deletion_does_not_overcall_stop_gained() {
        // A reverse-strand 1-base deletion: Perl reports frameshift_variant only
        // (gag/gTag -> E/VX), the shortened alt CDS giving a partial-codon X.
        let mut tx = make_test_transcript();
        tx.strand = vep_core::coordinate::Strand::Reverse;
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(seq) = vefc.translateable_seq.as_mut() {
                // First codon GAG (E = Glu).
                seq.replace_range(0..3, "GAG");
            }
        }

        let config = EffectsConfig::default();
        // Delete 1 base at CDS position 2 on a reverse-strand transcript.
        // Genomic coordinates: the variant spans 2 genomic positions.
        // On reverse strand, genomic start maps to higher CDS pos.
        let variant = InputVariant::new(
            "21".into(),
            25_000_051,
            25_000_051,
            b"A".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::FrameshiftVariant));
        assert!(
            !tc.consequences.contains(&Consequence::StopGained),
            "reverse-strand 1-base deletion should not overcall stop_gained"
        );
    }

    #[test]
    fn test_mnv_within_single_codon_does_not_false_positive_stop_gained() {
        // An MNV changing two bases within one codon without creating a stop guards
        // the secondary stop check against applying the edit at the wrong CDS offset
        // on the second anchor.
        let mut tx = make_test_transcript();
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(seq) = vefc.translateable_seq.as_mut() {
                // CDS begins with ATG (M) then GAC (D).
                // Variant changes AC->GT at CDS pos 5 (bases 5-6), turning GAC->GGT (G),
                // but should not introduce any stop codon. Codon 2 keeps the edit off
                // the start codon, where Perl's `_inv_start_altered` would call any
                // non-ATG result start_lost.
                seq.replace_range(0..6, "ATGGAC");
            }
        }

        let config = EffectsConfig::default();
        // Change CDS bases 5-6 (genomic 25_000_054..25_000_055).
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_055,
            b"AC".to_vec(),
            b"GT".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::MissenseVariant));
        assert!(!tc.consequences.contains(&Consequence::StopGained));
        assert!(!tc.consequences.contains(&Consequence::StopLost));
    }

    #[test]
    fn test_cross_codon_mnv_stop_gained_suppresses_missense() {
        // Cross-codon MNV: ATG CAG (MQ) -> ATT TAG (I*), so Perl reports stop_gained
        // (and does not co-emit missense_variant).
        let mut tx = make_test_transcript();
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(seq) = vefc.translateable_seq.as_mut() {
                seq.replace_range(0..6, "ATGCAG");
            }
        }

        let config = EffectsConfig::default();
        // Change CDS bases 3-4 (genomic 25_000_052..25_000_053): GC->TT.
        let variant = InputVariant::new(
            "21".into(),
            25_000_052,
            25_000_053,
            b"GC".to_vec(),
            b"TT".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::StopGained));
        assert!(!tc.consequences.contains(&Consequence::MissenseVariant));
        assert!(!tc.consequences.contains(&Consequence::SynonymousVariant));
    }

    #[test]
    fn test_inframe_deletion_can_create_stop_gained() {
        // Inframe deletion that creates a stop codon at the first affected codon.
        // Ref begins with TCC CAA (S Q). Delete 3 bases starting at CDS pos 2 (CCC),
        // leaving TAA ... which is Stop.
        let mut tx = make_test_transcript();
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(seq) = vefc.translateable_seq.as_mut() {
                seq.replace_range(0..6, "TCCCAA");
            }
        }

        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_051,
            25_000_053,
            b"CCC".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::InframeDeletion));
        assert!(tc.consequences.contains(&Consequence::StopGained));
    }

    #[test]
    fn test_unminimised_multi_allelic_deletion_is_inframe_deletion() {
        // VEP does not minimise the alleles of a multi-allelic record, so the allele
        // `TGGA/T` (CDS 6-9, keeping the T and deleting codon 3) reaches the coding
        // analysis with its shared base attached. VEP's `inframe_deletion` compares
        // codon windows (`GCTGGA` against `GCT`, a prefix match), so the term is
        // inframe_deletion, as it is for the minimised `GGA/-`.
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let unminimised = InputVariant::new(
            "21".into(),
            25_000_055,
            25_000_058,
            b"TGGA".to_vec(),
            b"T".to_vec(),
        );
        let minimised = InputVariant::new(
            "21".into(),
            25_000_056,
            25_000_058,
            b"GGA".to_vec(),
            b"-".to_vec(),
        );
        let a = calculate_consequences(&unminimised, &tx, &config).unwrap();
        let b = calculate_consequences(&minimised, &tx, &config).unwrap();
        assert!(a.consequences.contains(&Consequence::InframeDeletion));
        assert!(!a.consequences.contains(&Consequence::CodingSequenceVariant));
        assert_eq!(a.consequences, b.consequences);
    }

    #[test]
    fn test_spanning_deletion_suppresses_polypyr_and_adds_coding_sequence_variant() {
        // Deletion whose endpoints map to different introns, but the allele span overlaps
        // exon 2. Perl suppresses splice_polypyrimidine_tract_variant when exon overlap
        // exists, and emits coding_sequence_variant when the span overlaps the coding region.
        let tx = make_test_transcript();
        let config = EffectsConfig::default();

        // Delete from intron 1 start through intron 2 end, spanning exon 2 entirely.
        let variant = InputVariant::new(
            "21".into(),
            25_000_300,
            25_003_999,
            b"AAAA".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();

        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(tc.consequences.contains(&Consequence::SpliceDonorVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::SpliceAcceptorVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::SplicePolypyrimidineTractVariant));
    }

    #[test]
    fn test_insertion_at_transcript_end_is_downstream_only() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_006_001,
            25_006_000,
            b"-".to_vec(),
            b"A".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::DownstreamGeneVariant));
        assert!(!tc.consequences.contains(&Consequence::ThreePrimeUtrVariant));
        assert_eq!(tc.distance, Some(0));
    }

    #[test]
    fn test_non_coding_insertion_at_exon_intron_boundary_uses_context_variant() {
        let mut tx = make_test_transcript();
        tx.biotype = "processed_transcript".into();
        tx.translation = None;
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(mapper) = vefc.mapper.as_mut() {
                mapper.cdna_coding_start = 0;
                mapper.cdna_coding_end = 0;
            }
        }
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_300,
            25_000_299,
            b"-".to_vec(),
            b"A".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::SpliceRegionVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant));
    }

    #[test]
    fn test_intronic_insertion_does_not_combine_donor_5th_and_donor_region() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // Insertion fully within intron 1: one side maps to donor+5, the other to donor+4.
        // Perl-style single-anchor behavior should not combine both donor labels.
        let variant = InputVariant::new(
            "21".into(),
            25_000_305,
            25_000_304,
            b"-".to_vec(),
            b"A".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::SpliceDonorRegionVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::SpliceDonor5thBaseVariant));
        assert!(!tc.consequences.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn test_spanning_deletion_calls_donor_and_donor_5th_base() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();

        // Delete across the exon1/intron1 join, covering intron positions +0..+4.
        // Ensembl VEP emits both splice_donor_variant and splice_donor_5th_base_variant.
        let variant = InputVariant::new(
            "21".into(),
            25_000_299,
            25_000_304,
            b"AAAAAA".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::SpliceDonorVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::SpliceDonor5thBaseVariant));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::SpliceDonorRegionVariant));
        assert!(!tc.consequences.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn test_intronic_deletion_donor_5th_base_suppresses_donor_region() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();

        // Delete intron positions +3..+4; endpoint-union logic would incorrectly combine
        // donor_region (+3) with donor_5th_base (+4). Overlap-based logic should suppress
        // donor_region when donor_5th_base is present.
        let variant = InputVariant::new(
            "21".into(),
            25_000_303,
            25_000_304,
            b"AA".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::SpliceDonor5thBaseVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::SpliceDonorRegionVariant));
        assert!(!tc.consequences.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn test_intronic_deletion_donor_region_suppresses_splice_region() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();

        // Delete intron positions +2..+3: within donor_region and splice_region windows.
        // Perl precedence keeps splice_region suppressed when donor_region is present.
        let variant = InputVariant::new(
            "21".into(),
            25_000_302,
            25_000_303,
            b"AA".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::SpliceDonorRegionVariant));
        assert!(!tc.consequences.contains(&Consequence::SpliceRegionVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::SpliceDonor5thBaseVariant));
        assert!(!tc.consequences.contains(&Consequence::SpliceDonorVariant));
    }

    #[test]
    fn test_coding_model_non_protein_biotype_skips_non_coding_context_term() {
        let mut tx = make_test_transcript();
        tx.biotype = "non_stop_decay".into();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::MissenseVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptVariant));
    }

    #[test]
    fn test_short_noncoding_intron_produces_non_coding_transcript_variant() {
        // A deletion inside a frameshift intron (<= 12 bp) of a non-coding
        // transcript overlaps no exon, so Perl's `non_coding_exon_variant` is
        // false and `within_non_coding_gene` yields non_coding_transcript_variant
        // (both in `VariationEffect`); the intron flags are skipped for a region
        // overlapping the frameshift intron, so no intron_variant either.
        let mut tx = make_test_transcript();
        tx.biotype = "processed_transcript".into();
        tx.translation = None;
        if let Some(vefc) = tx.vefc.as_mut() {
            if let Some(mapper) = vefc.mapper.as_mut() {
                mapper.cdna_coding_start = 0;
                mapper.cdna_coding_end = 0;
            }
            if let Some(intron) = vefc.introns.get_mut(0) {
                intron.start = 25_000_300;
                intron.end = 25_000_304;
            }
        }
        if let Some(intron) = tx.introns.get_mut(0) {
            intron.start = 25_000_300;
            intron.end = 25_000_304;
        }

        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_301,
            25_000_302,
            b"AA".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        crate::test_helpers::assert_consequence_set_eq(&tc, &["non_coding_transcript_variant"]);
    }

    #[test]
    fn test_transcript_metadata_populated() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_005_000,
            25_005_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = calculate_consequences(&variant, &tx, &config);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(&*tc.transcript_id, "ENST00000000001");
        assert_eq!(&*tc.gene_id, "ENSG00000000001");
        assert_eq!(tc.gene_symbol.as_deref(), Some("TEST1"));
        assert_eq!(tc.biotype.as_deref(), Some("protein_coding"));
        assert_eq!(tc.strand, 1);
        assert_eq!(tc.feature_type, FeatureType::Transcript);
    }

    /// Ensembl VEP 116.2 on the GRCh38-release116 corpus prints
    /// `ENSP00000286808.3` for ENST00000286808 (translation version 3) under
    /// `--protein --protein_version`, and `ENSP00000286808` without the flag.
    #[test]
    fn protein_version_appends_the_translation_version() {
        let mut tx = make_test_transcript();
        tx.protein_id = Some("ENSP00000286808".into());
        tx.translation.as_mut().unwrap().version = Some(3);
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let off = calculate_consequences(&variant, &tx, &EffectsConfig::default()).unwrap();
        assert_eq!(off.protein_id.as_deref(), Some("ENSP00000286808"));
        let on = EffectsConfig {
            protein_version: true,
            ..Default::default()
        };
        let tc = calculate_consequences(&variant, &tx, &on).unwrap();
        assert_eq!(tc.protein_id.as_deref(), Some("ENSP00000286808.3"));
        assert_eq!(off.consequences, tc.consequences);
    }

    /// OutputFactory.pm appends the version only when the identifier does not
    /// already end in `.<digits>` and the translation has one.
    #[test]
    fn protein_version_leaves_a_versioned_or_unversioned_identifier_alone() {
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let on = EffectsConfig {
            protein_version: true,
            ..Default::default()
        };
        let mut tx = make_test_transcript();
        tx.protein_id = Some("NP_001234.2".into());
        tx.translation.as_mut().unwrap().version = Some(3);
        let tc = calculate_consequences(&variant, &tx, &on).unwrap();
        assert_eq!(tc.protein_id.as_deref(), Some("NP_001234.2"));

        let mut tx = make_test_transcript();
        tx.translation.as_mut().unwrap().version = None;
        let tc = calculate_consequences(&variant, &tx, &on).unwrap();
        assert_eq!(tc.protein_id.as_deref(), Some("ENSP00000000001"));
    }

    #[test]
    fn test_hgvs_suppressed_when_compute_hgvs_off() {
        let tx = make_test_transcript();
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );

        let off = EffectsConfig::default();
        assert!(!off.compute_hgvs, "compute_hgvs must default to false");
        let tc_off = calculate_consequences(&variant, &tx, &off).unwrap();
        assert!(tc_off.hgvsc.is_none(), "HGVSc must be None when gated off");
        assert!(tc_off.hgvsp.is_none(), "HGVSp must be None when gated off");

        let on = EffectsConfig {
            compute_hgvs: true,
            ..Default::default()
        };
        let tc_on = calculate_consequences(&variant, &tx, &on).unwrap();
        assert!(tc_on.hgvsc.is_some(), "HGVSc must be Some when gated on");
        assert!(tc_on.hgvsp.is_some(), "HGVSp must be Some when gated on");

        // The gate must change only the HGVS fields, never the consequence call.
        assert_eq!(
            tc_off.consequences, tc_on.consequences,
            "compute_hgvs must not affect the consequence set"
        );
    }

    /// The shift reads the record's class (TranscriptVariationAllele.pm:119-123,
    /// 420; Utils/Sequence.pm:335-398): every allele of `-/A/AA` or `TT/-/-`
    /// walks, no allele of `TT/-/T` or `GT/G/AT` does, a lone allele walks by its
    /// own class, and an allele given on the reverse strand walks as its
    /// forward-strand bases (:459).
    #[test]
    fn shift_motif_follows_the_record_class_and_the_forward_strand() {
        let allele = |record: Option<&str>, reference: &str, alt: &str| {
            let (start, end) = if reference == "-" {
                (101, 100)
            } else {
                (100, 99 + reference.len() as u64)
            };
            let mut v = InputVariant::new(
                "1".into(),
                start,
                end,
                reference.as_bytes().to_vec(),
                alt.as_bytes().to_vec(),
            );
            v.record_allele_string_multi = record.map(Arc::from);
            v
        };
        assert_eq!(shift_motif(&allele(None, "-", "A")), Some(b"A".to_vec()));
        assert_eq!(shift_motif(&allele(None, "TT", "-")), Some(b"TT".to_vec()));
        assert_eq!(shift_motif(&allele(None, "T", "A")), None);
        assert_eq!(
            shift_motif(&allele(Some("-/A/AA"), "-", "AA")),
            Some(b"AA".to_vec())
        );
        assert_eq!(
            shift_motif(&allele(Some("TT/-/-"), "TT", "-")),
            Some(b"TT".to_vec())
        );
        assert_eq!(shift_motif(&allele(Some("TT/-/T"), "TT", "-")), None);
        assert_eq!(shift_motif(&allele(Some("TT/-/T"), "T", "-")), None);
        assert_eq!(shift_motif(&allele(Some("GT/G/AT"), "T", "-")), None);
        let mut reverse = allele(None, "-", "TG");
        reverse.strand = vep_core::coordinate::Strand::Reverse;
        assert_eq!(shift_motif(&reverse), Some(b"CA".to_vec()));
    }

    /// An allele given on the reverse strand analyses as its forward-strand form
    /// (`feature_seq`, VariationFeatureOverlapAllele.pm:253-256): the same terms,
    /// codons, amino acids and HGVS notations, shifted or not.
    #[test]
    fn reverse_strand_input_analyses_as_its_forward_form() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        let tx = make_compact_test_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(163, b'A')]);
        let config = EffectsConfig {
            reference_fasta: Some(Arc::new(fasta)),
            enable_indel_3prime_shift: true,
            compute_hgvs: true,
            ..Default::default()
        };
        let pairs = [
            ((160, 160, "T", "A"), (160, 160, "A", "T")),
            ((160, 160, "T", "-"), (160, 160, "A", "-")),
            ((161, 160, "-", "TG"), (161, 160, "-", "CA")),
        ];
        for ((s, e, r, a), (rs, re, rr, ra)) in pairs {
            let forward = InputVariant::new(
                "1".into(),
                s,
                e,
                r.as_bytes().to_vec(),
                a.as_bytes().to_vec(),
            );
            let mut reverse = InputVariant::new(
                "1".into(),
                rs,
                re,
                rr.as_bytes().to_vec(),
                ra.as_bytes().to_vec(),
            );
            reverse.strand = vep_core::coordinate::Strand::Reverse;
            let f = calculate_consequences(&forward, &tx, &config).unwrap();
            let g = calculate_consequences(&reverse, &tx, &config).unwrap();
            assert_eq!(f.consequences, g.consequences, "{forward:?}");
            assert_eq!(
                (f.codons, f.amino_acids),
                (g.codons, g.amino_acids),
                "{forward:?}"
            );
            assert!(f.hgvsc.is_some());
            assert_eq!(
                (f.hgvsc, f.hgvsp, f.hgvs_offset),
                (g.hgvsc, g.hgvsp, g.hgvs_offset)
            );
            assert_eq!(f.shift_length, g.shift_length);
        }
    }

    /// `hgvs_transcript` compares the input's strand with the transcript's
    /// (TranscriptVariationAllele.pm:1463-1466), so under `--shift_hgvs 0` an
    /// inserted `C` given on the forward strand against a reverse-strand transcript
    /// duplicates the base before it, `c.60dup`, while the same insertion given on
    /// the reverse strand, `-/G -1`, keeps the default order, `c.61dup`; shifted,
    /// both are `c.62dup` with the same protein notation.
    #[test]
    fn shift_hgvs_0_reads_the_lookup_order_from_the_input_strand() {
        use crate::concordance_tests::{
            compact_on_chr1, hgvsc_fasta_fixture, make_descending_reverse_strand_transcript,
        };
        let tx = compact_on_chr1(make_descending_reverse_strand_transcript());
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(5989, b'C'), (5990, b'C'), (5991, b'C')]);
        let fasta = Arc::new(fasta);
        let config = |shift_hgvs| EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            compute_hgvs: true,
            hgvs: crate::hgvs::HgvsOptions {
                shift_hgvs,
                ..Default::default()
            },
            ..Default::default()
        };
        let forward = InputVariant::new("1".into(), 5991, 5990, b"-".to_vec(), b"C".to_vec());
        let mut reverse = InputVariant::new("1".into(), 5991, 5990, b"-".to_vec(), b"G".to_vec());
        reverse.strand = vep_core::coordinate::Strand::Reverse;
        let hgvs = |variant: &InputVariant, shift_hgvs: bool| {
            let tc = calculate_consequences(variant, &tx, &config(shift_hgvs)).unwrap();
            (tc.hgvsc, tc.hgvsp, tc.hgvs_offset)
        };
        let (f_shifted, r_shifted) = (hgvs(&forward, true), hgvs(&reverse, true));
        assert_eq!(f_shifted.0.as_deref(), Some("ENST00000000001.1:c.62dup"));
        assert_eq!(f_shifted, r_shifted);
        assert!(f_shifted.1.is_some());
        assert_eq!(
            hgvs(&forward, false).0.as_deref(),
            Some("ENST00000000001.1:c.60dup")
        );
        assert_eq!(
            hgvs(&reverse, false).0.as_deref(),
            Some("ENST00000000001.1:c.61dup")
        );
    }

    /// `perform_shift` (TranscriptVariationAllele.pm:291) walks the 1,000-base
    /// flank `_genomic_shift` reads (:427-437): on a reference that is one long run
    /// a 1-base deletion moves 1,000 bases and a 2-base one 999 (:303, the loop
    /// bound `length(post_seq) - indel_length` inclusive from 0), in either
    /// direction; an allele longer than the flank moves 1,000 going 3' (:304) and a
    /// 1,001-base allele going 5' does not move at all (:303, a loop from 1 to 0).
    #[test]
    fn shift_3prime_bound_matches_perform_shift() {
        use crate::concordance_tests::{
            compact_on_chr1, hgvsc_fasta_fixture, make_compact_test_transcript,
            make_descending_reverse_strand_transcript,
        };
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let forward = make_compact_test_transcript();
        let reverse = compact_on_chr1(make_descending_reverse_strand_transcript());
        let deletion = |start: u64, len: u64| {
            InputVariant::new(
                "1".into(),
                start,
                start + len - 1,
                vec![b'T'; len as usize],
                b"-".to_vec(),
            )
        };

        assert_eq!(perform_shift_bound(1, false), 1000);
        assert_eq!(perform_shift_bound(2, true), 999);
        assert_eq!(perform_shift_bound(1000, false), 1);
        assert_eq!(perform_shift_bound(1001, false), 1000);
        assert_eq!(perform_shift_bound(1001, true), 0);
        assert_eq!(perform_shift_bound(1500, true), 1000);

        assert_eq!(
            shift_indel_3prime_coords(&deletion(300, 1), &forward, &fasta),
            Some((1300, 1300))
        );
        assert_eq!(
            shift_indel_3prime_coords(&deletion(300, 2), &forward, &fasta),
            Some((1299, 1300))
        );
        let insertion = InputVariant::new("1".into(), 300, 299, b"-".to_vec(), b"T".to_vec());
        assert_eq!(
            shift_indel_3prime_coords(&insertion, &forward, &fasta),
            Some((1300, 1299))
        );
        assert_eq!(
            shift_indel_3prime_coords(&deletion(5000, 1), &reverse, &fasta),
            Some((4000, 4000))
        );
        assert_eq!(
            shift_indel_3prime_coords(&deletion(5000, 2), &reverse, &fasta),
            Some((4001, 4002))
        );
        assert_eq!(
            shift_indel_3prime_coords(&deletion(2000, 1001), &forward, &fasta),
            Some((3000, 4000))
        );
        assert_eq!(
            shift_indel_3prime_coords(&deletion(3000, 1001), &reverse, &fasta),
            None
        );
    }

    /// Under `--shift_3prime 1` the display columns describe the shifted allele
    /// and `shift_length` counts the move (OutputFactory.pm:1690-1695, 1725-1741);
    /// a substitution moves nothing and reports 0, and without the flag an indel
    /// keeps the input span and reports 0 too.
    #[test]
    fn shift_3prime_moves_the_display_columns_and_counts_the_shift() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        let tx = make_compact_test_transcript();
        // CDS 1 is genomic 150: deleting the `T` at 160 (CDS 11) before `TT` and an
        // `A` at 163 moves two bases to CDS 13.
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(163, b'A')]);
        let fasta = Arc::new(fasta);
        let deletion = InputVariant::new("1".into(), 160, 160, b"T".to_vec(), b"-".to_vec());
        let snv = InputVariant::new("1".into(), 160, 160, b"T".to_vec(), b"A".to_vec());

        let unshifted = EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            ..Default::default()
        };
        let tc = calculate_consequences(&deletion, &tx, &unshifted).unwrap();
        assert_eq!(tc.shift_length, 0);
        assert_eq!(tc.cds_position.as_deref(), Some("11"));

        let shifted = EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            enable_indel_3prime_shift: true,
            ..Default::default()
        };
        let tc = calculate_consequences(&deletion, &tx, &shifted).unwrap();
        assert_eq!(tc.shift_length, 2);
        assert_eq!(tc.cds_position.as_deref(), Some("13"));
        let tc = calculate_consequences(&snv, &tx, &shifted).unwrap();
        assert_eq!(tc.shift_length, 0);
        assert_eq!(tc.cds_position.as_deref(), Some("11"));
    }

    /// The reference row of `--process_ref_homs` is read at the span the record's
    /// alternate shifts to and counts no shift of its own (OutputFactory.pm:306-323,
    /// 1690): the `T` deletion at 160 moves two bases to CDS 13 and its reference
    /// row reads CDS 13 as a `synonymous_variant`; an inserted `T` before `TTA`
    /// moves two bases to CDS 13-14, as does the reference row of the same
    /// insertion written on the reverse strand as `-/A`, whose alternate is
    /// complemented before the shift. Without the flag the reference row keeps the
    /// record's span.
    #[test]
    fn a_reference_row_reads_the_alternates_shifted_span_with_shift_length_0() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        let tx = make_compact_test_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(163, b'A')]);
        let fasta = Arc::new(fasta);
        let shifted = EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            enable_indel_3prime_shift: true,
            ..Default::default()
        };
        let deletion = InputVariant::new("1".into(), 160, 160, b"T".to_vec(), b"-".to_vec());
        let row = deletion.reference_allele_row();
        assert_eq!(row.reference_row_shift_allele.as_deref(), Some("-"));
        let tc = calculate_consequences(&row, &tx, &shifted).unwrap();
        assert_eq!(tc.shift_length, 0);
        assert_eq!(tc.cds_position.as_deref(), Some("13"));
        assert_eq!(tc.consequences.as_slice(), [Consequence::SynonymousVariant]);

        let forward = InputVariant::new("1".into(), 161, 160, b"-".to_vec(), b"T".to_vec());
        let mut reverse = InputVariant::new("1".into(), 161, 160, b"-".to_vec(), b"A".to_vec());
        reverse.strand = vep_core::coordinate::Strand::Reverse;
        for insertion in [&forward, &reverse] {
            let tc = calculate_consequences(insertion, &tx, &shifted).unwrap();
            assert_eq!(tc.shift_length, 2, "{insertion:?}");
            assert_eq!(tc.cds_position.as_deref(), Some("13-14"), "{insertion:?}");
            let tc =
                calculate_consequences(&insertion.reference_allele_row(), &tx, &shifted).unwrap();
            assert_eq!(tc.shift_length, 0, "{insertion:?}");
            assert_eq!(tc.cds_position.as_deref(), Some("13-14"), "{insertion:?}");
            assert_eq!(
                tc.consequences.as_slice(),
                [Consequence::SynonymousVariant],
                "{insertion:?}"
            );
        }

        let unshifted = EffectsConfig {
            reference_fasta: Some(fasta),
            ..Default::default()
        };
        let tc = calculate_consequences(&row, &tx, &unshifted).unwrap();
        assert_eq!(tc.shift_length, 0);
        assert_eq!(tc.cds_position.as_deref(), Some("11"));
    }

    /// The transcript window is read on the input span: a deletion at the last
    /// base of the downstream window whose shift carries it past the window keeps
    /// its row, `downstream_gene_variant` at the input's distance, while one
    /// placed a base beyond the window has none, shifted or not.
    #[test]
    fn shift_3prime_keeps_a_row_the_input_span_earned_past_the_window() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        let tx = make_compact_test_transcript();
        assert_eq!(tx.end, 6100);
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(6603, b'A'), (6604, b'A')]);
        let fasta = Arc::new(fasta);
        let config = |shift_3prime| EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            enable_indel_3prime_shift: shift_3prime,
            downstream_distance: 500,
            ..Default::default()
        };
        let at_edge = InputVariant::new("1".into(), 6600, 6600, b"T".to_vec(), b"-".to_vec());
        for (shift_3prime, shift_length) in [(false, 0), (true, 2)] {
            let tc = calculate_consequences(&at_edge, &tx, &config(shift_3prime)).unwrap();
            assert_eq!(
                tc.consequences.as_slice(),
                [Consequence::DownstreamGeneVariant]
            );
            assert_eq!(tc.distance, Some(500));
            assert_eq!(tc.shift_length, shift_length);
        }
        let past_edge = InputVariant::new("1".into(), 6601, 6601, b"T".to_vec(), b"-".to_vec());
        assert!(calculate_consequences(&past_edge, &tx, &config(false)).is_none());
        assert!(calculate_consequences(&past_edge, &tx, &config(true)).is_none());

        // Upstream, the shift moves the allele towards the transcript, here the
        // whole 1,000-base bound along the run of `T`: a deletion a base beyond
        // the upstream window shifts into the transcript and still has no row,
        // because the fetch never produced one; at the edge it keeps its row,
        // described where the shift leaves it.
        assert_eq!(tx.start, 100);
        let upstream_config = |shift_3prime| EffectsConfig {
            upstream_distance: 50,
            ..config(shift_3prime)
        };
        let before_edge = InputVariant::new("1".into(), 49, 49, b"T".to_vec(), b"-".to_vec());
        assert!(calculate_consequences(&before_edge, &tx, &upstream_config(false)).is_none());
        assert!(calculate_consequences(&before_edge, &tx, &upstream_config(true)).is_none());
        let at_upstream_edge = InputVariant::new("1".into(), 50, 50, b"T".to_vec(), b"-".to_vec());
        let tc = calculate_consequences(&at_upstream_edge, &tx, &upstream_config(false)).unwrap();
        assert_eq!(
            tc.consequences.as_slice(),
            [Consequence::UpstreamGeneVariant]
        );
        assert_eq!(tc.distance, Some(50));
        let tc = calculate_consequences(&at_upstream_edge, &tx, &upstream_config(true)).unwrap();
        assert_eq!(tc.shift_length, 1000);
        assert!(!tc.consequences.contains(&Consequence::UpstreamGeneVariant));
    }

    #[test]
    fn exon_intron_numbers_follow_the_config_flag() {
        let tx = make_test_transcript();
        let exonic = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let intronic = InputVariant::new(
            "21".into(),
            25_001_000,
            25_001_000,
            b"C".to_vec(),
            b"A".to_vec(),
        );

        let off = EffectsConfig::default();
        assert!(!off.compute_exon_intron_numbers);
        let exonic_off = calculate_consequences(&exonic, &tx, &off).unwrap();
        let intronic_off = calculate_consequences(&intronic, &tx, &off).unwrap();
        assert!(exonic_off.exon.is_none() && exonic_off.intron.is_none());
        assert!(intronic_off.exon.is_none() && intronic_off.intron.is_none());

        let on = EffectsConfig {
            compute_exon_intron_numbers: true,
            ..Default::default()
        };
        let exonic_on = calculate_consequences(&exonic, &tx, &on).unwrap();
        let intronic_on = calculate_consequences(&intronic, &tx, &on).unwrap();
        assert_eq!(exonic_on.exon.as_deref(), Some("1/3"));
        assert_eq!(exonic_on.intron, None);
        assert_eq!(intronic_on.exon, None);
        assert_eq!(intronic_on.intron.as_deref(), Some("1/2"));

        // The flag changes only the two ordinals.
        assert_eq!(exonic_off.consequences, exonic_on.consequences);
        assert_eq!(exonic_off.cdna_position, exonic_on.cdna_position);
        assert_eq!(exonic_off.codons, exonic_on.codons);
        assert_eq!(intronic_off.consequences, intronic_on.consequences);
    }

    #[test]
    fn test_hgvsc_populated_for_coding_snv() {
        let tx = make_test_transcript();
        let config = EffectsConfig {
            compute_hgvs: true,
            ..Default::default()
        };
        // Missense SNV: CDS pos 5 (GCT -> GAT = Ala -> Asp)
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(
            tc.hgvsc.is_some(),
            "HGVSc should be populated for coding SNV"
        );
        let hgvsc = tc.hgvsc.unwrap();
        assert!(
            hgvsc.starts_with("ENST00000000001.1:c."),
            "HGVSc should start with transcript prefix, got: {hgvsc}"
        );
        assert!(
            hgvsc.contains(">"),
            "SNV HGVSc should contain '>', got: {hgvsc}"
        );
    }

    #[test]
    fn test_hgvsp_populated_for_missense() {
        let tx = make_test_transcript();
        let config = EffectsConfig {
            compute_hgvs: true,
            ..Default::default()
        };
        // Missense SNV: CDS pos 5 (GCT -> GAT = Ala -> Asp)
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::MissenseVariant),
            "Should be missense, got: {:?}",
            tc.consequences
        );
        assert!(tc.hgvsp.is_some(), "HGVSp should be populated for missense");
        let hgvsp = tc.hgvsp.unwrap();
        assert!(
            hgvsp.starts_with("ENSP00000000001.1:p."),
            "HGVSp should start with protein prefix, got: {hgvsp}"
        );
        assert!(
            hgvsp.contains("Ala2Asp"),
            "HGVSp should contain Ala2Asp, got: {hgvsp}"
        );
    }

    #[test]
    fn test_hgvsp_not_populated_for_intronic() {
        let tx = make_test_transcript();
        let config = EffectsConfig {
            compute_hgvs: true,
            ..Default::default()
        };
        let variant = InputVariant::new(
            "21".into(),
            25_001_000,
            25_001_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(
            tc.hgvsp.is_none(),
            "HGVSp should not be populated for intronic variants"
        );
        assert!(
            tc.hgvsc.is_some(),
            "HGVSc should be populated even for intronic variants"
        );
    }

    #[test]
    fn test_transcript_ablation_deletion_encompasses_transcript() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            24_999_000, // before transcript start
            25_007_000, // after transcript end
            b"AAAA".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::TranscriptAblation],
            "Deletion encompassing entire transcript should return transcript_ablation only"
        );
        assert_eq!(tc.impact, Impact::HIGH);
        assert!(tc.cdna_position.is_none());
        assert!(tc.cds_position.is_none());
        assert!(tc.protein_position.is_none());
    }

    #[test]
    fn test_transcript_ablation_exact_transcript_boundaries() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_000, // exact transcript start
            25_006_000, // exact transcript end
            b"AAAA".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::TranscriptAblation],
            "Deletion matching exact transcript boundaries should return transcript_ablation"
        );
    }

    #[test]
    fn test_no_transcript_ablation_for_partial_deletion() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            24_999_000,
            25_003_000, // inside transcript
            b"AAAA".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(
            !tc.consequences.contains(&Consequence::TranscriptAblation),
            "Partial deletion should not produce transcript_ablation"
        );
    }

    #[test]
    fn test_no_transcript_ablation_for_insertion() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_006_001,
            25_006_000,
            b"-".to_vec(),
            b"AAAA".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(
            !tc.consequences.contains(&Consequence::TranscriptAblation),
            "Insertions should never produce transcript_ablation"
        );
    }

    #[test]
    fn test_consequences_sorted_by_so_rank() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_300,
            25_003_999,
            b"AAAA".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        // Verify consequences are sorted by SO rank (most severe first).
        let ranks: Vec<u32> = tc.consequences.iter().map(|c| c.rank()).collect();
        let mut sorted_ranks = ranks.clone();
        sorted_ranks.sort();
        assert_eq!(
            ranks, sorted_ranks,
            "Consequences should be sorted by SO rank (most severe first)"
        );
    }

    #[test]
    fn test_large_spanning_deletion_from_cds_into_utr_is_stop_lost() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        // A deletion from coding exon 1 into the 3' UTR of exon 3: Perl's
        // `cds_end` is undef (the span leaves the CDS), so the peptide alleles are
        // empty and `stop_lost` delegates to `_ins_del_stop_altered`, which
        // splices the cDNA span out of CDS + 3'UTR and re-reads the codon at the
        // original stop position. Without cached UTR bases that codon is unknown
        // and not a stop, so `stop_lost` fires and `coding_unknown` is suppressed.
        let variant = InputVariant::new(
            "21".into(),
            25_000_100,
            25_004_500,
            b"AAAA".to_vec(),
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::StopLost),
            "Deletion from the CDS into the 3' UTR should emit stop_lost. got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "Perl's coding_unknown is 0 when stop_lost fires. got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_insertion_analysis_cds_pos_uses_perl_insertion_start() {
        let bounds = CdsSpanBounds {
            cds_start: 748,
            cds_end: 747,
            translation_start: 249,
            translation_end: 250,
        };

        assert_eq!(insertion_analysis_cds_pos(747, Some(&bounds)), 748);
        assert_eq!(insertion_analysis_cds_pos(747, None), 747);
    }

    #[test]
    fn test_should_use_insertion_anchor_for_anchored_insertion_bounds() {
        let bounds = CdsSpanBounds {
            cds_start: 748,
            cds_end: 747,
            translation_start: 250,
            translation_end: 249,
        };
        let same_codon_bounds = CdsSpanBounds {
            cds_start: 747,
            cds_end: 747,
            translation_start: 249,
            translation_end: 249,
        };

        assert!(should_use_insertion_anchor(false, Some(&bounds)));
        assert!(!should_use_insertion_anchor(
            false,
            Some(&same_codon_bounds)
        ));
        assert!(should_use_insertion_anchor(true, Some(&same_codon_bounds)));
    }

    #[test]
    fn test_boundary_inframe_insertion_uses_inserted_field_display() {
        let tx = make_test_transcript();
        let config = EffectsConfig::default();
        let variant = InputVariant::new(
            "21".into(),
            25_000_059,
            25_000_058,
            b"-".to_vec(),
            b"ATT".to_vec(),
        );

        let tc = calculate_consequences(&variant, &tx, &config).unwrap();
        assert!(tc.consequences.contains(&Consequence::InframeInsertion));
        assert_eq!(tc.amino_acids.as_deref(), Some("-/I"));
        assert_eq!(tc.codons.as_deref(), Some("-/ATT"));
    }

    /// The compact transcript with its UTRs cached, so its cDNA is the cache's:
    /// 50 `C`, the coding sequence (`ATGGCTGGAAAATTCGAT`, then `GCT`), 1,701 `T`.
    fn compact_transcript_with_cached_cdna() -> Transcript {
        let mut tx = crate::concordance_tests::make_compact_test_transcript();
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.five_prime_utr = Some("C".repeat(50));
        vefc.three_prime_utr = Some("T".repeat(1701));
        assert_eq!(tx.spliced_seq().map(|s| &s[49..56]), Some("CATGGCT"));
        tx
    }

    /// Under `--use_transcript_ref` the reference allele is cut from the
    /// transcript's spliced cDNA when the variant maps onto one stretch of it
    /// (TranscriptVariation.pm:441-468 at release/116): a coding transcript's
    /// from the cache, any other transcript's from the FASTA. A variant in an
    /// intron or across an exon boundary keeps the input's allele, because
    /// Ensembl VEP builds every TranscriptVariation with `-no_ref_check`
    /// (AnnotationType/Transcript.pm:159), which skips the genomic branch
    /// (VariationFeatureOverlap.pm:368) and leaves the input's allele as the
    /// reference (`:418-424`). The row's terms follow the allele it is read
    /// against, as Ensembl VEP's codons do (`_get_alternate_cds` writes the
    /// alternate allele over the transcript's own bases), and `USED_REF` is
    /// given on the input's strand (`variation_feature_seq`).
    #[test]
    fn used_ref_comes_from_spliced_cdna_inside_exons_and_the_input_elsewhere() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        // Exon 1 is 100-399, CDS 1 is genomic 150; the FASTA is `T` throughout.
        let mismatch = InputVariant::new("1".into(), 152, 152, b"A".to_vec(), b"C".to_vec());
        let no_cdna = make_compact_test_transcript();
        assert_eq!(no_cdna.spliced_seq(), None);
        assert_eq!(
            transcript_reference_allele(&mismatch, &no_cdna, Some(&fasta)),
            Some(b"T".to_vec())
        );
        assert_eq!(transcript_reference_allele(&mismatch, &no_cdna, None), None);

        let tx = compact_transcript_with_cached_cdna();
        assert_eq!(
            transcript_reference_allele(&mismatch, &tx, None),
            Some(b"G".to_vec())
        );
        let intronic = InputVariant::new("1".into(), 1000, 1000, b"A".to_vec(), b"C".to_vec());
        assert_eq!(
            transcript_reference_allele(&intronic, &tx, Some(&fasta)),
            None
        );
        let into_intron = InputVariant::new("1".into(), 398, 402, b"AAAAA".to_vec(), b"-".to_vec());
        assert_eq!(
            transcript_reference_allele(&into_intron, &tx, Some(&fasta)),
            None
        );
        let over_intron = InputVariant::new("1".into(), 398, 2101, vec![b'A'; 1704], b"-".to_vec());
        assert_eq!(
            transcript_reference_allele(&over_intron, &tx, Some(&fasta)),
            None
        );
        let upstream = InputVariant::new("1".into(), 50, 50, b"A".to_vec(), b"C".to_vec());
        assert_eq!(
            transcript_reference_allele(&upstream, &tx, Some(&fasta)),
            None
        );

        let on = EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            use_transcript_ref: true,
            ..Default::default()
        };
        let off = EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            ..Default::default()
        };
        let tc = calculate_consequences(&mismatch, &tx, &on).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("G"));
        assert_eq!(
            calculate_consequences(&mismatch, &tx, &off)
                .unwrap()
                .used_ref,
            None
        );
        let tc = calculate_consequences(&intronic, &tx, &on).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("A"));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));

        // Given `GA` over the transcript's `GC` at 396-397, the last coding bases
        // before the exon's splice region, the alternate `GC` differs from the
        // input at 397 (a splice-region base) and from the transcript nowhere,
        // so the row keeps its coding term and loses the splice-region one.
        let mnv = InputVariant::new("1".into(), 396, 397, b"GA".to_vec(), b"GC".to_vec());
        let tc = calculate_consequences(&mnv, &tx, &on).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("GC"));
        assert_eq!(tc.consequences.as_slice(), [Consequence::SynonymousVariant]);
        let off_terms = calculate_consequences(&mnv, &tx, &off)
            .unwrap()
            .consequences;
        assert!(
            off_terms.contains(&Consequence::SpliceRegionVariant),
            "{off_terms:?}"
        );

        // The reference row is the transcript's allele.
        let reference_row = mismatch.reference_allele_row();
        assert_eq!(reference_row.alt_allele(), b"A");
        let reread = reference_row.against_reference_allele(b"G");
        assert_eq!(
            (reread.ref_allele.as_slice(), reread.alt_allele()),
            (&b"G"[..], &b"G"[..])
        );
        assert_eq!(reread.allele_string, "G/G");
        let tc = calculate_consequences(&reference_row, &tx, &on).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("G"));
        assert!(
            tc.consequences.contains(&Consequence::StartRetainedVariant),
            "{:?}",
            tc.consequences
        );

        // An allele given on the reverse strand is read on the forward strand
        // and its reference given back on the reverse strand.
        let mut reverse = InputVariant::new("1".into(), 152, 152, b"T".to_vec(), b"G".to_vec());
        reverse.strand = vep_core::coordinate::Strand::Reverse;
        let tc = calculate_consequences(&reverse, &tx, &on).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("C"));
    }

    /// The bases read from the FASTA for a cDNA span are sized by what the
    /// FASTA returns, never by the span's cDNA ends: a mapper from a malformed
    /// cache can number two abutting exon pieces with a jump of nearly 2^63
    /// between them, so a variant across the junction has a cDNA span no
    /// buffer could hold, while its eleven bases read from the FASTA as they
    /// are.
    #[test]
    fn used_ref_from_the_fasta_is_not_sized_by_a_jumping_cdna_numbering() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        use vep_core::transcript::{ExonCoordMapper, MapperPair};
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let mut tx = make_compact_test_transcript();
        assert_eq!(tx.spliced_seq(), None);
        let mapper = tx.vefc.as_mut().unwrap().mapper.as_mut().unwrap();
        let mut pairs = mapper.exon_coord_mapper.pairs.clone();
        // The compact transcript's exon 1 is genomic 100..=399, cDNA 1..=300;
        // it is split in two abutting pieces, and the second piece and the
        // later exons take the jumped numbering.
        let jump = i64::MAX as u64 - 10_000;
        for pair in &mut pairs[1..] {
            pair.from_start += jump - 151;
            pair.from_end += jump - 151;
        }
        let first = pairs.remove(0);
        let split_first = [
            MapperPair {
                from_start: 1,
                from_end: 150,
                to_start: 100,
                to_end: 249,
                ..first.clone()
            },
            MapperPair {
                from_start: jump,
                from_end: jump + 149,
                to_start: 250,
                to_end: 399,
                ..first
            },
        ];
        let pairs: Vec<MapperPair> = split_first.into_iter().chain(pairs).collect();
        mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
        let across = InputVariant::new("1".into(), 245, 255, vec![b'A'; 11], b"-".to_vec());
        assert_eq!(
            transcript_reference_allele(&across, &tx, Some(&fasta)),
            Some(vec![b'T'; 11])
        );
    }

    /// An insertion's reference allele is `-` (TranscriptVariation.pm:451) and a
    /// stretch of `N` in the transcript's sequence is kept as the reference
    /// allele, since the cDNA branch returns it without the `N` test of the
    /// genomic branch (VariationFeatureOverlap.pm:392, skipped under
    /// `-no_ref_check`). A reverse-strand transcript's cDNA is read back onto
    /// the forward strand (`:465`).
    #[test]
    fn used_ref_is_dash_for_an_insertion_and_keeps_an_all_n_exon_reference() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        use crate::test_helpers::{make_single_exon_transcript, SINGLE_EXON_UTR5};
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(152, b'N'), (153, b'N'), (154, b'N')]);
        let fasta = Arc::new(fasta);
        let tx = make_compact_test_transcript();
        let insertion = InputVariant::new("1".into(), 153, 152, b"-".to_vec(), b"GGG".to_vec());
        assert_eq!(
            transcript_reference_allele(&insertion, &tx, Some(&fasta)),
            Some(b"-".to_vec())
        );
        let in_ns = InputVariant::new("1".into(), 152, 154, b"TTT".to_vec(), b"-".to_vec());
        assert_eq!(
            transcript_reference_allele(&in_ns, &tx, Some(&fasta)),
            Some(b"NNN".to_vec())
        );
        let config = EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            use_transcript_ref: true,
            ..Default::default()
        };
        let tc = calculate_consequences(&insertion, &tx, &config).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("-"));
        let tc = calculate_consequences(&in_ns, &tx, &config).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("NNN"));

        let mut reverse =
            make_single_exon_transcript(vep_core::coordinate::Strand::Reverse, "ATGGCTTAA", "GGGG");
        reverse.vefc.as_mut().unwrap().five_prime_utr = Some("C".repeat(SINGLE_EXON_UTR5 as usize));
        assert!(reverse.spliced_seq().is_some());
        // The transcript's first base, a `C` of its 5' UTR, sits at the span's end.
        let first = InputVariant::new(
            "21".into(),
            reverse.end,
            reverse.end,
            b"A".to_vec(),
            b"T".to_vec(),
        );
        assert_eq!(
            transcript_reference_allele(&first, &reverse, None),
            Some(b"G".to_vec())
        );
        // The last base of the 3' UTR, a `G`, sits at the span's start.
        let last = InputVariant::new(
            "21".into(),
            reverse.start,
            reverse.start,
            b"A".to_vec(),
            b"T".to_vec(),
        );
        assert_eq!(
            transcript_reference_allele(&last, &reverse, None),
            Some(b"C".to_vec())
        );
        // The start codon, cDNA 31-33, read back onto the forward strand.
        let start_codon_end = reverse.end - SINGLE_EXON_UTR5;
        let codon = InputVariant::new(
            "21".into(),
            start_codon_end - 2,
            start_codon_end,
            b"AAA".to_vec(),
            b"-".to_vec(),
        );
        assert_eq!(
            transcript_reference_allele(&codon, &reverse, None),
            Some(b"CAT".to_vec())
        );
    }

    /// A non-coding transcript has no cached cDNA, so its reference allele is
    /// read from the FASTA over the exon piece the span maps to, as Ensembl's
    /// `Transcript::spliced_seq` reads every exon from the genome: on a
    /// reverse-strand transcript the bases come back on the forward strand (the
    /// strand the variant is analysed on), and a span reaching into an intron
    /// maps to a gap, so the input's allele stands.
    #[test]
    fn used_ref_of_a_non_coding_reverse_transcript_is_read_from_the_fasta() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        use vep_core::transcript::{Exon, ExonCoordMapper, MapperPair, TranscriptMapper};
        // Two exons at 1201-1300 (cDNA 1-100, cDNA 1 at genomic 1300) and 1001-1100
        // (cDNA 101-200); the FASTA is `T` but for 1300 `A`, 1299 `C`, 1100 `G`.
        let mut tx = make_compact_test_transcript();
        tx.strand = vep_core::coordinate::Strand::Reverse;
        tx.biotype = "lncRNA".into();
        tx.start = 1001;
        tx.end = 1300;
        tx.translation = None;
        tx.protein_id = None;
        tx.cdna_coding_start = None;
        tx.cdna_coding_end = None;
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation_start = None;
        tx.translation_end = None;
        let exon = |start, end, rank| Exon {
            stable_id: None,
            start,
            end,
            rank,
            phase: -1,
            end_phase: -1,
        };
        tx.exons = vec![exon(1201, 1300, 1), exon(1001, 1100, 2)];
        tx.introns = Vec::new();
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.translateable_seq = None;
        vefc.introns = Vec::new();
        vefc.sorted_exons = tx.exons.clone();
        let pair = |from_start, from_end, to_start, to_end| MapperPair {
            from_start,
            from_end,
            to_start,
            to_end,
            ori: -1,
        };
        vefc.mapper = Some(TranscriptMapper {
            start_phase: -1,
            cdna_coding_start: 0,
            cdna_coding_end: 0,
            exon_coord_mapper: ExonCoordMapper::new(vec![
                pair(1, 100, 1201, 1300),
                pair(101, 200, 1001, 1100),
            ]),
        });
        assert_eq!(tx.spliced_seq(), None);
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(1300, b'A'), (1299, b'C'), (1100, b'G')]);

        let first_base = InputVariant::new("1".into(), 1300, 1300, b"G".to_vec(), b"T".to_vec());
        assert_eq!(
            transcript_reference_allele(&first_base, &tx, Some(&fasta)),
            Some(b"A".to_vec())
        );
        let three = InputVariant::new("1".into(), 1298, 1300, b"GGG".to_vec(), b"-".to_vec());
        assert_eq!(
            transcript_reference_allele(&three, &tx, Some(&fasta)),
            Some(b"TCA".to_vec())
        );
        let second_exon = InputVariant::new("1".into(), 1100, 1100, b"A".to_vec(), b"C".to_vec());
        assert_eq!(
            transcript_reference_allele(&second_exon, &tx, Some(&fasta)),
            Some(b"G".to_vec())
        );
        let into_intron = InputVariant::new("1".into(), 1100, 1101, b"AA".to_vec(), b"-".to_vec());
        assert_eq!(
            transcript_reference_allele(&into_intron, &tx, Some(&fasta)),
            None
        );
        let intronic = InputVariant::new("1".into(), 1150, 1150, b"A".to_vec(), b"C".to_vec());
        assert_eq!(
            transcript_reference_allele(&intronic, &tx, Some(&fasta)),
            None
        );
        assert_eq!(transcript_reference_allele(&first_base, &tx, None), None);

        let config = EffectsConfig {
            reference_fasta: Some(Arc::new(fasta)),
            use_transcript_ref: true,
            ..Default::default()
        };
        let tc = calculate_consequences(&three, &tx, &config).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("TCA"));
        assert!(tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant));
    }

    /// `USED_REF` of a deletion follows its 3' shift where Ensembl VEP rotates
    /// the reference allele's sequence along with the alternate allele's, through
    /// the codon path under `--shift_3prime 1` (TranscriptVariationAllele.pm:674,
    /// 813-820, 1330-1346) and the protein-notation path under `--hgvs`
    /// (:1730-1737, 1771-1786), so only on a row whose shifted span has both
    /// protein coordinates: a coding deletion of `TA` moved one base prints
    /// `AT`, measured on Ensembl VEP 116.2 (`CA` deleted at chromosome 21
    /// 34716447 and moved one base on ENST00000349499 prints `USED_REF=AC` under
    /// either flag, `CA` under neither), while a deletion in the 3' UTR and a
    /// run without either flag print the bases as the input wrote them. Under
    /// `--process_ref_homs` the codon path rotates nothing (the reference row
    /// memoises the codon first, :802-804) and the protein path still does.
    #[test]
    fn used_ref_of_a_shifted_coding_deletion_is_rotated_by_the_shift() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        // The FASTA is `T` throughout but for an `A` at 163 and 4500 and a `G` at
        // 165 and 4502, so deleting `TA` at 162-163 (CDS 13-14) and at 4499-4500
        // (the 3' UTR, cDNA 1000-1001) each move one base, onto `AT`.
        let (_dir, fasta) =
            hgvsc_fasta_fixture(&[(163, b'A'), (165, b'G'), (4500, b'A'), (4502, b'G')]);
        let fasta = Arc::new(fasta);
        let tx = crate::concordance_tests::make_compact_test_transcript();
        let coding = InputVariant::new("1".into(), 162, 163, b"TA".to_vec(), b"-".to_vec());
        let utr = InputVariant::new("1".into(), 4499, 4500, b"TA".to_vec(), b"-".to_vec());
        let config = |shift_3prime: bool, hgvs: bool| EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            enable_indel_3prime_shift: shift_3prime,
            compute_hgvs: hgvs,
            use_transcript_ref: true,
            ..Default::default()
        };
        let tc = calculate_consequences(&coding, &tx, &config(true, false)).unwrap();
        assert_eq!(tc.shift_length, 1);
        assert!(tc.protein_position.is_some());
        assert_eq!(tc.used_ref.as_deref(), Some("AT"));
        let tc = calculate_consequences(&utr, &tx, &config(true, false)).unwrap();
        assert_eq!(tc.shift_length, 1);
        assert!(tc.consequences.contains(&Consequence::ThreePrimeUtrVariant));
        assert_eq!(tc.used_ref.as_deref(), Some("TA"));

        let tc = calculate_consequences(&coding, &tx, &config(false, true)).unwrap();
        assert_eq!(tc.shift_length, 0);
        assert_eq!(tc.hgvs_offset, Some(1));
        assert!(tc.hgvsp.is_some());
        assert_eq!(tc.used_ref.as_deref(), Some("AT"));
        let tc = calculate_consequences(&utr, &tx, &config(false, true)).unwrap();
        assert_eq!(tc.hgvs_offset, Some(1));
        assert_eq!(tc.used_ref.as_deref(), Some("TA"));

        let tc = calculate_consequences(&coding, &tx, &config(false, false)).unwrap();
        assert_eq!(tc.shift_length, 0);
        assert_eq!(tc.used_ref.as_deref(), Some("TA"));

        // Under `--process_ref_homs` the reference row reads the codon before any
        // alternate allele hands it a shift, so only the protein path rotates.
        let with_reference_row = |shift_3prime: bool, hgvs: bool| EffectsConfig {
            process_ref_homs: true,
            ..config(shift_3prime, hgvs)
        };
        let tc = calculate_consequences(&coding, &tx, &with_reference_row(true, false)).unwrap();
        assert_eq!(tc.shift_length, 1);
        assert_eq!(tc.used_ref.as_deref(), Some("TA"));
        let tc = calculate_consequences(&coding, &tx, &with_reference_row(false, true)).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("AT"));
        let tc = calculate_consequences(&coding, &tx, &with_reference_row(true, true)).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("AT"));
    }

    /// The `Allele` a coding row prints for a shifted insertion is the written
    /// allele rotated as `shift_feature_seqs` rotates it: one place per base
    /// along a forward-strand transcript, the length less the shift along a
    /// reverse-strand one, so `AAG` moved one base on the reverse-strand APP
    /// prints `GAA` as Ensembl VEP 116.2 does, an allele given on the reverse
    /// strand rotates as written (`CTT` to `TCT`), and a shift of the allele's
    /// length or more on the reverse strand leaves it. A row with no translation
    /// coordinate, an unmoved row, a reference row and a deletion keep the
    /// allele as given.
    #[test]
    fn a_shifted_coding_insertion_prints_the_rotated_allele() {
        use vep_core::coordinate::Strand;
        let insertion = |alt: &str, strand: Strand| {
            let mut v = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), alt.into());
            v.strand = strand;
            v
        };
        let rotated = |v: &InputVariant, strand, span, shift| {
            shifted_display_allele(v, strand, span, shift).map(|a| a.to_string())
        };
        let span = Some((242, 242));
        let forward_given = insertion("AAG", Strand::Forward);
        assert_eq!(
            rotated(&forward_given, Strand::Reverse, span, 1).as_deref(),
            Some("GAA")
        );
        assert_eq!(
            rotated(&forward_given, Strand::Forward, span, 1).as_deref(),
            Some("AGA")
        );
        assert_eq!(
            rotated(&forward_given, Strand::Forward, span, 4).as_deref(),
            Some("AGA")
        );
        assert_eq!(rotated(&forward_given, Strand::Reverse, span, 3), None);
        assert_eq!(rotated(&forward_given, Strand::Reverse, span, 4), None);
        assert_eq!(
            rotated(&insertion("CTT", Strand::Reverse), Strand::Reverse, span, 1).as_deref(),
            Some("TCT")
        );
        assert_eq!(rotated(&forward_given, Strand::Reverse, None, 1), None);
        assert_eq!(rotated(&forward_given, Strand::Reverse, span, 0), None);
        let mut reference_row = insertion("-", Strand::Forward);
        reference_row.reference_allele_row = true;
        assert_eq!(rotated(&reference_row, Strand::Reverse, span, 1), None);
        let deletion = InputVariant::new("21".into(), 101, 103, b"TCT".to_vec(), b"-".to_vec());
        assert_eq!(rotated(&deletion, Strand::Reverse, span, 1), None);
    }

    /// Both paths move the bases towards the transcript's 3' end, so on the
    /// reverse strand a shift of one moves the last base to the front. The
    /// protein path wraps a shift longer than the allele on either strand; the
    /// codon path wraps it on the forward strand and moves nothing on the
    /// reverse strand, where `shift_feature_seqs` moves `length - shift` places:
    /// Ensembl VEP 116.2 prints `USED_REF=GA` for `GA` deleted at chromosome 21
    /// 30348449 and moved three bases on ENST00000334160 under `--shift_3prime 1`
    /// alone, and `AG` once `--hgvs` is added.
    #[test]
    fn used_ref_rotation_follows_the_transcript_strand_and_the_path() {
        use vep_core::coordinate::Strand;
        let rotated = |allele: &[u8], rotation, strand| {
            let mut allele = allele.to_vec();
            rotate_along_transcript(&mut allele, rotation, strand);
            allele
        };
        for rotation in [UsedRefRotation::Codon, UsedRefRotation::Protein] {
            assert_eq!(rotated(b"CA", rotation(1), Strand::Forward), b"AC");
            assert_eq!(rotated(b"AC", rotation(1), Strand::Reverse), b"CA");
            assert_eq!(rotated(b"GCA", rotation(2), Strand::Reverse), b"CAG");
            assert_eq!(rotated(b"GATT", rotation(6), Strand::Forward), b"TTGA");
            assert_eq!(rotated(b"GA", rotation(2), Strand::Reverse), b"GA");
            assert_eq!(rotated(b"A", rotation(3), Strand::Forward), b"A");
            assert_eq!(rotated(b"-", rotation(3), Strand::Reverse), b"-");
        }
        assert_eq!(
            rotated(b"GA", UsedRefRotation::Protein(3), Strand::Reverse),
            b"AG"
        );
        assert_eq!(
            rotated(b"GA", UsedRefRotation::Codon(3), Strand::Reverse),
            b"GA"
        );
        assert_eq!(
            rotated(b"CAC", UsedRefRotation::Protein(4), Strand::Reverse),
            b"CCA"
        );
        assert_eq!(
            rotated(b"CAC", UsedRefRotation::Codon(4), Strand::Reverse),
            b"CAC"
        );

        let mut tc = TranscriptConsequence {
            shift_length: 3,
            ..Default::default()
        };
        assert_eq!(
            used_ref_rotation(&tc, Some((2, 3)), false),
            Some(UsedRefRotation::Codon(3))
        );
        assert_eq!(used_ref_rotation(&tc, Some((2, 3)), true), None);
        assert_eq!(used_ref_rotation(&tc, None, false), None);
        tc.hgvs_offset = Some(-3);
        assert_eq!(
            used_ref_rotation(&tc, Some((2, 3)), false),
            Some(UsedRefRotation::Codon(3))
        );
        tc.hgvsp = Some("ENSP00000000001.1:p.Leu2fs".into());
        assert_eq!(
            used_ref_rotation(&tc, Some((2, 3)), false),
            Some(UsedRefRotation::Protein(3))
        );
        assert_eq!(
            used_ref_rotation(&tc, Some((2, 3)), true),
            Some(UsedRefRotation::Protein(3))
        );
        tc.shift_length = 0;
        assert_eq!(
            used_ref_rotation(&tc, None, false),
            Some(UsedRefRotation::Protein(3))
        );
        tc.hgvs_offset = None;
        assert_eq!(used_ref_rotation(&tc, Some((2, 3)), false), None);
    }

    /// A single-exon reverse-strand coding transcript on chromosome 1 at
    /// 1001-1046: the 3' UTR `GGGG` at 1001-1004, the coding sequence
    /// `ATGCTCTCTTAA` at 1005-1016 (its start codon at 1014-1016) and the 5' UTR
    /// at 1017-1046, with the FASTA spelling the same bases on the forward
    /// strand (`CCCC TTAAGAGAGCAT`, then `T`).
    fn reverse_coding_transcript_on_chr1(
    ) -> (tempfile::TempDir, vep_fasta::IndexedFasta, Transcript) {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        use crate::test_helpers::make_single_exon_transcript;
        let mut tx = make_single_exon_transcript(
            vep_core::coordinate::Strand::Reverse,
            "ATGCTCTCTTAA",
            "GGGG",
        );
        let delta = tx.start - 1001;
        tx.chr = "1".into();
        tx.start -= delta;
        tx.end -= delta;
        for exon in &mut tx.exons {
            exon.start -= delta;
            exon.end -= delta;
        }
        for coord in [
            &mut tx.coding_region_start,
            &mut tx.coding_region_end,
            &mut tx.translation_start,
            &mut tx.translation_end,
        ] {
            *coord = coord.map(|c| c - delta);
        }
        let vefc = tx.vefc.as_mut().unwrap();
        for exon in &mut vefc.sorted_exons {
            exon.start -= delta;
            exon.end -= delta;
        }
        let mapper = &mut vefc.mapper.as_mut().unwrap().exon_coord_mapper;
        for pair in &mut mapper.pairs {
            pair.to_start -= delta;
            pair.to_end -= delta;
        }
        mapper.rebuild_genomic_index();
        assert_eq!((tx.start, tx.end), (1001, 1046));
        assert_eq!(
            (tx.coding_region_start, tx.coding_region_end),
            (Some(1005), Some(1016))
        );
        let mut edits: Vec<(usize, u8)> = (1001..=1004).map(|p| (p, b'C')).collect();
        edits.extend(
            b"TTAAGAGAGCAT"
                .iter()
                .enumerate()
                .map(|(i, &base)| (1005 + i, base)),
        );
        let (dir, fasta) = hgvsc_fasta_fixture(&edits);
        (dir, fasta, tx)
    }

    /// On a reverse-strand transcript the 3' shift runs towards lower
    /// coordinates and the rotation the other way: `AGA` shifted one base
    /// prints `AAG`. A shift longer than the deletion rotates only through the
    /// protein path: `GA` moved three bases prints `GA` under `--shift_3prime 1`
    /// alone and `AG` under `--hgvs`.
    #[test]
    fn used_ref_of_a_reverse_strand_deletion_rotates_the_other_way() {
        let (_dir, fasta, tx) = reverse_coding_transcript_on_chr1();
        assert_eq!(tx.spliced_seq(), None);
        let fasta = Arc::new(fasta);
        let config = |shift_3prime: bool, hgvs: bool| EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            enable_indel_3prime_shift: shift_3prime,
            compute_hgvs: hgvs,
            use_transcript_ref: true,
            ..Default::default()
        };
        // `AGA` at 1008-1010 moves one base, onto `AAG` at 1007-1009.
        let one = InputVariant::new("1".into(), 1008, 1010, b"AGA".to_vec(), b"-".to_vec());
        let tc = calculate_consequences(&one, &tx, &config(true, false)).unwrap();
        assert_eq!(tc.shift_length, 1);
        assert!(tc.protein_position.is_some());
        assert_eq!(tc.used_ref.as_deref(), Some("AAG"));
        let tc = calculate_consequences(&one, &tx, &config(false, true)).unwrap();
        assert_eq!(tc.hgvs_offset, Some(-1));
        assert!(tc.hgvsp.is_some());
        assert_eq!(tc.used_ref.as_deref(), Some("AAG"));
        let tc = calculate_consequences(&one, &tx, &config(false, false)).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("AGA"));

        // `GA` at 1011-1012 moves three bases, onto `AG` at 1008-1009.
        let three = InputVariant::new("1".into(), 1011, 1012, b"GA".to_vec(), b"-".to_vec());
        let tc = calculate_consequences(&three, &tx, &config(true, false)).unwrap();
        assert_eq!(tc.shift_length, 3);
        assert!(tc.protein_position.is_some());
        assert_eq!(tc.used_ref.as_deref(), Some("GA"));
        let tc = calculate_consequences(&three, &tx, &config(false, true)).unwrap();
        assert_eq!(tc.hgvs_offset, Some(-3));
        assert!(tc.hgvsp.is_some());
        assert_eq!(tc.used_ref.as_deref(), Some("AG"));
        let tc = calculate_consequences(&three, &tx, &config(true, true)).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("AG"));
        let misspelt = InputVariant::new("1".into(), 1011, 1012, b"TA".to_vec(), b"-".to_vec());
        let tc = calculate_consequences(&misspelt, &tx, &config(true, false)).unwrap();
        assert_eq!(tc.shift_length, 3);
        assert_eq!(tc.used_ref.as_deref(), Some("GA"));

        // Given on the reverse strand, `AGA` is `TCT`: the codon path rotates the
        // allele as written, the protein path hands back the forward-strand
        // rotation.
        let mut given_on_reverse = one.clone();
        given_on_reverse.strand = vep_core::coordinate::Strand::Reverse;
        given_on_reverse.ref_allele = b"TCT".to_vec();
        given_on_reverse.allele_string = "TCT/-".into();
        let tc = calculate_consequences(&given_on_reverse, &tx, &config(true, false)).unwrap();
        assert_eq!(tc.shift_length, 1);
        assert_eq!(tc.used_ref.as_deref(), Some("TTC"));
        let tc = calculate_consequences(&given_on_reverse, &tx, &config(false, true)).unwrap();
        assert_eq!(tc.hgvs_offset, Some(-1));
        assert_eq!(tc.used_ref.as_deref(), Some("CTT"));
        let tc = calculate_consequences(&given_on_reverse, &tx, &config(false, false)).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("TCT"));
    }

    /// The HGVS shift of a deletion moves the bases the reference allele's
    /// object names (`_genomic_shift`, TranscriptVariationAllele.pm:439-449), so
    /// under `--use_transcript_ref` a deletion whose input spells the deleted
    /// bases wrongly shifts along the transcript's bases: Ensembl VEP 116.2
    /// writes `c.2027_2028del` with `HGVS_OFFSET=2` for `21 34716445 TGC T` on
    /// ENST00000349499, whose sequence there is `ACAC`, where the input's `GC`
    /// shifts nowhere without the flag.
    #[test]
    fn hgvs_shifts_the_bases_the_transcript_reference_names() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        // CDS 11-16 (genomic 160-165) read `ACACAC` from the FASTA, the transcript's
        // cDNA here, since the compact transcript caches no UTR.
        let (_dir, fasta) = hgvsc_fasta_fixture(&[
            (160, b'A'),
            (161, b'C'),
            (162, b'A'),
            (163, b'C'),
            (164, b'A'),
            (165, b'C'),
        ]);
        let fasta = Arc::new(fasta);
        let tx = make_compact_test_transcript();
        let misspelt = InputVariant::new("1".into(), 160, 161, b"GC".to_vec(), b"-".to_vec());
        let on = EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            compute_hgvs: true,
            use_transcript_ref: true,
            ..Default::default()
        };
        let tc = calculate_consequences(&misspelt, &tx, &on).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("AC"));
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:c.15_16del"));
        assert_eq!(tc.hgvs_offset, Some(4));
        let off = EffectsConfig {
            use_transcript_ref: false,
            ..on
        };
        let tc = calculate_consequences(&misspelt, &tx, &off).unwrap();
        assert_eq!(tc.used_ref, None);
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:c.11_12del"));
        assert_eq!(tc.hgvs_offset, None);
    }

    /// `hgvs_protein` returns nothing only for the reference allele and a
    /// non-ACGT allele (TranscriptVariationAllele.pm:1684, 1687 at release/116),
    /// so an alternate allele equal to the transcript's reference, which the
    /// flag reads as a synonymous change, is named `p.<Aa><pos>=`, while the
    /// reference row stays without a protein notation.
    #[test]
    fn hgvsp_names_an_alternate_equal_to_the_transcript_reference() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let tx = compact_transcript_with_cached_cdna();
        let on = EffectsConfig {
            reference_fasta: Some(Arc::new(fasta)),
            compute_hgvs: true,
            use_transcript_ref: true,
            ..Default::default()
        };
        // Codon 83 of the coding sequence (genomic 396-398) reads `GCT`.
        assert_eq!(tx.spliced_seq().map(|s| &s[296..299]), Some("GCT"));
        let mnv = InputVariant::new("1".into(), 396, 397, b"GA".to_vec(), b"GC".to_vec());
        let tc = calculate_consequences(&mnv, &tx, &on).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("GC"));
        assert_eq!(tc.consequences.as_slice(), [Consequence::SynonymousVariant]);
        assert_eq!(tc.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ala83="));
        let tc = calculate_consequences(&mnv.reference_allele_row(), &tx, &on).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("GC"));
        assert_eq!(tc.hgvsp, None);
        // Without the flag the codons are the transcript's as well, so the
        // notation is the same while the input's `GA` keeps the splice-region term.
        let off = EffectsConfig {
            use_transcript_ref: false,
            ..on
        };
        let tc = calculate_consequences(&mnv, &tx, &off).unwrap();
        assert_eq!(tc.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ala83="));
        assert!(tc.consequences.contains(&Consequence::SpliceRegionVariant));
    }

    /// The module documentation names what the `--use_transcript_ref` port
    /// leaves out: the edited transcript sequences of a BAM-edited cache, which
    /// switch the flag on by themselves, and the warning of an edited RefSeq
    /// transcript.
    #[test]
    fn documented_subset_names_the_bam_edited_cache_auto_enable() {
        let module_docs: String = include_str!("consequences.rs")
            .lines()
            .take_while(|l| l.starts_with("//") || l.trim().is_empty())
            .filter(|l| l.starts_with("//!"))
            .map(|l| l.trim_start_matches("//!").trim())
            .collect::<Vec<_>>()
            .join(" ");
        for phrase in [
            "BAM-edited cache",
            "switch the flag on by themselves",
            "`VEP::CacheDir`",
            "`apply_edits`",
            "`--bam`",
            "`invalid_alleles`",
        ] {
            assert!(module_docs.contains(phrase), "{phrase}: {module_docs}");
        }
    }
}
