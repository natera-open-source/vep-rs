// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Consequence calculator: assigns SO consequence terms to variant-transcript overlaps.
//!
//! Perl citations name modules of ensembl-variation release/116
//! (`Bio/EnsEMBL/Variation/...`); the `VEP::...` modules are in ensembl-vep release/116.2.
//!
//! Under `--use_transcript_ref` every row reads the reference allele
//! [`transcript_reference_allele`] gives; a BAM-edited cache switches the flag
//! on by itself (`VEP::CacheDir`, the `bam` entry of the cache's info file), and
//! an edited transcript's reference allele is cut from the edited cDNA the cache
//! stores (`VEP::AnnotationType::Transcript`, `apply_edits`), coding or not. A
//! row whose alternate allele equals the bases of the transcript's first
//! `_rna_edit` on the codon path is marked for the `Transcript-assembly
//! mismatch` warning ([`transcript_assembly_mismatch`]), which the caller
//! writes. Not ported: the `--bam` flag that edits the models at run time,
//! which is refused.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::sync::Arc;

use smallvec::{smallvec, SmallVec};
use vep_core::consequence::{
    Consequence, ConsequenceList, Impact, ProteinIdDisplay, TranscriptConsequence,
};
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

use crate::coding::{
    classify_inframe_indel_by_peptide, compute_codon_window_peptide_alleles,
    compute_peptide_alleles, get_codon_change, InframeIndelKind,
};
use crate::coding::{frameshift_stop_gained_full_cds, frameshift_stop_gained_in_codon_window};
use crate::mapper::{map_genomic_to_transcript, partial_codon, CdsSpanBounds, TranscriptPosition};

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
    /// `--hgvsp_use_prediction`, `--remove_hgvsp_version` and, for the span the
    /// HGVSc gate reads and the predicates `hgvs_protein` reads,
    /// `--shift_3prime`, which an embedder sets as it sets
    /// `enable_indel_3prime_shift`.
    pub hgvs: crate::hgvs::HgvsOptions,
    /// The transcript consequences that carry a
    /// [`LofteeContext`](vep_core::consequence::LofteeContext) from construction,
    /// so the exon/intron model and genomic CDS bounds travel to the plugin
    /// phase, where the transcript is out of scope.
    ///
    /// [`LofteeContextRows::None`] by default: `TranscriptConsequence.loftee_ctx`
    /// stays `None` and output is byte-identical to a non-plugin run. The
    /// context clones the transcript's exon and intron lists per row, so a
    /// plugin that reads one consequence class asks for those rows alone.
    pub loftee_context: LofteeContextRows,
    /// The transcript consequences that carry the transcript's coding, 3' UTR
    /// and protein sequences
    /// ([`TranscriptSequences`](vep_core::consequence::TranscriptSequences))
    /// from construction, for a plugin that translates a changed coding
    /// sequence in the plugin phase. [`LofteeContextRows::None`] by default:
    /// `TranscriptConsequence.transcript_sequences` stays `None`. The sequences
    /// are copied per row, so a plugin asks for the rows of its one term.
    pub context_sequences: LofteeContextRows,
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
    /// as `TranscriptConsequence.used_ref`; an unedited non-coding transcript's
    /// cDNA is read from `reference_fasta`.
    pub use_transcript_ref: bool,
    /// `--process_ref_homs`: the record's reference row is written, and its
    /// codon is read before any alternate allele's, which keeps the codon path
    /// from rotating `used_ref` with a 3' shift (see `used_ref_rotation`).
    pub process_ref_homs: bool,
}

/// Which transcript consequences carry a
/// [`LofteeContext`](vep_core::consequence::LofteeContext) ([`EffectsConfig::loftee_context`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LofteeContextRows {
    /// No row.
    #[default]
    None,
    /// Every row: LoFTEE reads every protein-coding row.
    Every,
    /// The rows whose consequence terms include this one: UTRAnnotator reads
    /// the `5_prime_UTR_variant` rows alone.
    WithTerm(Consequence),
}

impl LofteeContextRows {
    /// Whether the row `consequences` carries the context.
    pub fn wants(self, consequences: &[Consequence]) -> bool {
        match self {
            Self::None => false,
            Self::Every => true,
            Self::WithTerm(term) => consequences.contains(&term),
        }
    }
}

impl Default for EffectsConfig {
    fn default() -> Self {
        Self {
            upstream_distance: 5000,
            downstream_distance: 5000,
            reference_fasta: None,
            enable_indel_3prime_shift: false,
            hgvs: crate::hgvs::HgvsOptions::default(),
            loftee_context: LofteeContextRows::None,
            context_sequences: LofteeContextRows::None,
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
/// relative to the transcript structure (UTR, coding, intron, etc.). The
/// consequence holds a shared handle on `transcript` and reads the transcript's
/// attributes through it.
///
/// On a model the genome is the reference of
/// ([`crate::coding::genome_is_reference`]) an allele that keeps the bases on
/// both sides of bases the model holds and the genome lacks adds its other
/// bases between two neighbours of the transcript and replaces nothing: the
/// row takes the term, the codons, the residues, the protein position and the
/// HGVS notations, with the offset of their 3' shift, of the same bases
/// written as an insertion at that point
/// ([`insertion_between_kept_bases`]), and its cDNA and CDS position cells
/// keep the record's span, the model's bases inside it. Ensembl VEP reads the
/// allele in place of the whole span, the model's bases with it
/// (BaseTranscriptVariation.pm:263-264; TranscriptVariationAllele.pm:2409-2410,
/// 2421).
pub fn calculate_consequences(
    variant: &InputVariant,
    transcript: &Arc<Transcript>,
    config: &EffectsConfig,
) -> Option<TranscriptConsequence> {
    let mut row = consequences_of_the_record(variant, transcript, config)?;
    let as_an_insertion = insertion_between_kept_bases(variant, transcript)
        .and_then(|insertion| consequences_of_the_record(&insertion, transcript, config));
    if let Some(read) = as_an_insertion {
        row.consequences = read.consequences;
        row.impact = read.impact;
        row.protein_position = read.protein_position;
        row.protein_start = read.protein_start;
        row.amino_acids = read.amino_acids;
        row.codons = read.codons;
        row.domains = read.domains;
        row.hgvsc = read.hgvsc;
        row.hgvsp = read.hgvsp;
        row.hgvs_offset = read.hgvs_offset;
        row.hgvs_intron_offsets = read.hgvs_intron_offsets;
    }
    Some(row)
}

/// The insertion a record states on a model the genome is the reference of
/// ([`crate::coding::genome_is_reference`]) when its allele keeps the bases on
/// both sides of bases the model holds and the genome lacks: the record's
/// span lies on two pieces of the mapper that the genome joins end to end and
/// the model numbers more than one apart, and the allele starts with the
/// reference bases before the join, ends with those after it and holds other
/// bases between. The insertion is those bases at the join, on the strand the
/// record gives its alleles on. `None` for any other record.
fn insertion_between_kept_bases(
    variant: &InputVariant,
    transcript: &Transcript,
) -> Option<InputVariant> {
    use crate::display::Segment;
    use vep_core::coordinate::Strand;
    if !crate::coding::genome_is_reference(transcript)
        || variant.reference_allele_row
        || variant.is_structural
        || variant.start >= variant.end
    {
        return None;
    }
    let mapper = transcript.vefc.as_ref()?.mapper.as_ref()?;
    let segments = crate::display::genomic2cdna(
        mapper,
        variant.start as i64,
        variant.end as i64,
        transcript.strand.as_i8(),
    );
    let [Segment::Coord {
        start: first,
        end: first_end,
        ..
    }, Segment::Coord { start: second, .. }] = segments.as_slice()
    else {
        return None;
    };
    if *second <= first_end + 1 {
        return None;
    }
    let (reference, alt) = variant.forward_strand_alleles();
    let length = reference.len();
    if length as u64 != variant.end - variant.start + 1 || alt.len() <= length {
        return None;
    }
    // The pieces come in transcript order: the first holds the record's bases
    // 5' of the join on the transcript's strand.
    let before_the_join = usize::try_from(first_end - first + 1)
        .ok()
        .filter(|bases| *bases < length)?;
    let lead = match transcript.strand {
        Strand::Forward => before_the_join,
        Strand::Reverse => length - before_the_join,
    };
    if !alt.starts_with(&reference[..lead]) || !alt.ends_with(&reference[lead..]) {
        return None;
    }
    let added = &alt[lead..alt.len() - (length - lead)];
    let added = match variant.strand {
        Strand::Reverse if vep_core::codon::seq_is_dna(added) => {
            vep_core::codon::reverse_complement(added)
        }
        _ => added.to_vec(),
    };
    let mut insertion = variant.clone();
    insertion.start = variant.start + lead as u64;
    insertion.end = insertion.start - 1;
    insertion.ref_allele = b"-".to_vec();
    insertion.allele_string = format!("-/{}", String::from_utf8_lossy(&added));
    insertion.alt_alleles = vec![added];
    insertion.variant_class = vep_core::variant::VariantClass::Insertion;
    // The insertion is a record of its own: the 3' shift reads its class,
    // where the record's other alleles would hold it still
    // (`record_class_shifts`).
    insertion.record_allele_string_multi = None;
    Some(insertion)
}

/// The row of `variant` against `transcript`, read on the record as it is
/// written; [`calculate_consequences`] states what it is given.
fn consequences_of_the_record(
    variant: &InputVariant,
    transcript: &Arc<Transcript>,
    config: &EffectsConfig,
) -> Option<TranscriptConsequence> {
    // Under `--shift_3prime 1` the allele is analysed where the shift leaves it:
    // every predicate reads the moved span (`get_all_OverlapConsequences` moves
    // the feature's start and end for each one, BaseVariationFeatureOverlapAllele.pm:269-282,
    // and `_bvfo_preds` adds the offset, :454-463) and the rotated allele
    // (`shift_feature_seqs`, TranscriptVariationAllele.pm:1330, applied here as
    // `coding::rotate_shifted_allele` reads it on the transcript's strand), as do the cDNA,
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
    tc.transcript_assembly_mismatch =
        transcript_assembly_mismatch(analysed, transcript, translation_span);
    // `hgvs_protein` takes the shifted bases on the forward strand and gives
    // them back on the input's (`create_shift_hash`, TranscriptVariationAllele.pm:373);
    // `shift_feature_seqs` rotates the allele as the input wrote it (:1334-1343),
    // here as read on the transcript's strand.
    tc.used_ref = used_ref.map(|mut reference| {
        let reference_rewritten =
            !crate::coding::annotated_in_coding_intron(transcript, analysed.start, analysed.end)
                || intronic_allele_reads_as_coding(transcript, analysed);
        let as_given = match used_ref_rotation(
            &tc,
            translation_span,
            config.process_ref_homs,
            reference_rewritten,
        ) {
            Some(UsedRefRotation::Protein(shift)) => {
                crate::coding::rotate_shifted_allele(
                    &mut reference,
                    shift,
                    vep_core::coordinate::Strand::Forward,
                    transcript.strand,
                );
                on_input_strand(&reference)
            }
            Some(UsedRefRotation::Codon(shift)) => {
                let mut as_given = on_input_strand(&reference);
                crate::coding::rotate_shifted_allele(
                    &mut as_given,
                    shift,
                    variant.strand,
                    transcript.strand,
                );
                as_given
            }
            None => on_input_strand(&reference),
        };
        Arc::from(String::from_utf8_lossy(&as_given).as_ref())
    });
    Some(tc)
}

/// `invalid_alleles` (`TranscriptVariationAllele::codon`,
/// TranscriptVariationAllele.pm:813-837 at release/116): on a transcript with an
/// `_rna_edit` attribute, an alternate allele's row whose span has both
/// translation coordinates and whose allele is DNA (`:815`, the one path that
/// reaches the test) compares the allele, read on the transcript's strand with
/// `-` as the empty sequence (`feature_seq`, `:821-822`), with the bases of the
/// transcript's first `_rna_edit` (`:832-833`); equal, the row is one Ensembl
/// VEP warns about as `Transcript-assembly mismatch in <Uploaded_variation>`
/// (OutputFactory.pm:1717-1718), unless `--use_given_ref`. `variant` is the
/// allele as analysed, on the forward strand. A deletion's empty allele equals
/// the empty bases of a first edit that deletes, as the rule reads it.
pub fn transcript_assembly_mismatch(
    variant: &InputVariant,
    transcript: &Transcript,
    translation_span: Option<(i64, i64)>,
) -> bool {
    if variant.reference_allele_row || translation_span.is_none() {
        return false;
    }
    let Some(first_edit_bases) = transcript
        .refseq_alignment()
        .and_then(|alignment| alignment.first_edit_bases.as_deref())
    else {
        return false;
    };
    let alt = variant.alt_allele();
    if !vep_core::codon::seq_is_dna(alt) {
        return false;
    }
    let alt = alt.to_ascii_uppercase();
    let feature_seq = crate::coding::feature_seq(
        &alt,
        vep_core::coordinate::Strand::Forward,
        transcript.strand,
    );
    *feature_seq == *first_edit_bases.as_bytes()
}

/// Whether Ensembl VEP's `coding` pre-consequence predicate holds for an allele
/// annotated in an intron inside the coding region, so that `hgvs_protein`
/// passes its test (TranscriptVariationAllele.pm:1731) and rewrites the
/// reference allele: the transcript has a frameshift intron, `_overlapped_exons`
/// reads its exons 12 bases wider (BaseTranscriptVariation.pm:861-866) and the
/// allele lies within that stretch of one
/// (BaseVariationFeatureOverlapAllele.pm:511-557).
pub(crate) fn intronic_allele_reads_as_coding(
    transcript: &Transcript,
    allele: &InputVariant,
) -> bool {
    let (lo, hi) = (allele.start.min(allele.end), allele.start.max(allele.end));
    transcript.facts().vefc_has_frameshift_intron
        && allele_overlaps_any_exon_stretched(transcript, lo, hi)
}

/// The path that rotates a shifted deletion's `USED_REF`, with the shift in
/// bases. Ensembl VEP rotates the reference allele's sequence along with a
/// shifted alternate allele only through the codon and protein paths, so only
/// where the shifted span has both translation coordinates: `Protein` when
/// `hgvs_protein` writes the deleted bases `perform_shift` moved to into the
/// reference allele (TranscriptVariationAllele.pm:1771-1786) after the
/// translation test at the shifted offset (:1730-1737), which a protein
/// notation attests on a row whose allele as annotated passes `coding`
/// (`reference_rewritten`); otherwise `Codon` under `--shift_3prime 1`, where
/// `display_codon_allele_string` hands the alternate allele's shift to the
/// reference allele (:674) and `codon` rotates its sequence
/// (`shift_feature_seqs`, :1330-1345) once `translation_start` and
/// `translation_end` are defined (:813-820). Either path moves the bases
/// towards the transcript's 3' end, `perform_shift` one place per step with
/// the base that falls off the end coming back at the start (:291-350), and
/// vep-rs rotates both by [`crate::coding::rotate_shifted_allele`]; the
/// protein path rotates the bases on the forward strand before they are given
/// back on the input's, the codon path rotates them as the input wrote them. A
/// non-coding row, a row in a UTR and a deletion whose shifted span leaves the
/// coding sequence keep the bases as the input wrote them. `codon` rotates only
/// on the call that memoises the reference allele's codon (:802-804), and under
/// `--process_ref_homs` the reference row reads that codon, with no shift,
/// before any alternate allele's row hands it one, so the codon path then
/// rotates nothing (`reference_row_first`); the protein path writes the
/// sequence outright. An offset that counts a move of the HGVSc alone
/// (`hgvs_offset_of_hgvsc_alone`) is no shift of the record, and rotates
/// nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum UsedRefRotation {
    Codon(u64),
    Protein(u64),
}

fn used_ref_rotation(
    tc: &TranscriptConsequence,
    translation_span: Option<(i64, i64)>,
    reference_row_first: bool,
    reference_rewritten: bool,
) -> Option<UsedRefRotation> {
    let protein_shift = match (&tc.hgvsp, tc.hgvs_offset) {
        (Some(_), Some(offset))
            if offset != 0 && !tc.hgvs_offset_of_hgvsc_alone && reference_rewritten =>
        {
            Some(offset.unsigned_abs())
        }
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

/// [`TranscriptConsequence::shifted_allele`]: the inserted bases as the record
/// wrote them, rotated as `shift_feature_seqs` rotates them when `codon` runs
/// on a shifted allele (TranscriptVariationAllele.pm:820, 1330-1345), which it
/// does on a row whose shifted span has both translation coordinates, read on
/// the transcript's strand ([`crate::coding::rotate_shifted_allele`]). `None`
/// wherever the `Allele` is the allele as given: a row nothing moved, a row with
/// no translation coordinate, a reference row, any allele but an insertion of
/// DNA, and a rotation by a multiple of the allele's length.
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
    let steps = crate::coding::rotate_shifted_allele(
        &mut written,
        shift,
        variant.strand,
        transcript_strand,
    );
    if steps == 0 {
        return None;
    }
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
/// (`:418-424`). A BAM-edited model's cDNA and a coding transcript's are the
/// cache's ([`Transcript::spliced_seq`]); any other transcript's is read from
/// `fasta` over the exons the span covers, as Ensembl's `Transcript::spliced_seq`
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
    let Some(attribute) = transcript.attributes.iter().find(|a| &*a.code == "ncRNA") else {
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
/// the protein identifier the row prints gains `.<translation version>` unless
/// the row prints none, the translation has no version or the identifier
/// already ends in `.<digits>`.
fn append_translation_version(tc: &mut TranscriptConsequence, transcript: &Transcript) {
    let Some(version) = transcript.translation.as_ref().and_then(|t| t.version) else {
        return;
    };
    let versioned = match tc.protein_id() {
        None => return,
        Some(id) => id
            .rsplit_once('.')
            .is_some_and(|(_, tail)| !tail.is_empty() && tail.bytes().all(|b| b.is_ascii_digit())),
    };
    if !versioned {
        tc.protein_id_display = ProteinIdDisplay::Versioned(version);
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
    transcript: &Arc<Transcript>,
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
        // A span covering the transcript is `within_feature`, so Ensembl VEP
        // asks `hgvs_transcript` for this row like any other
        // (OutputFactory.pm:1721, 1750-1754): `_var2transcript_slice_coords`
        // refuses only a span with both ends outside on one side and clamps
        // each end onto the transcript (TranscriptVariationAllele.pm:2729-2733),
        // so a deletion covering it at a zero 3' shift is written end to end.
        // `hgvs_protein` has no peptide position for such a span (:1730-1738).
        let hgvs = if config.compute_hgvs {
            crate::hgvs::generate_hgvs(
                input,
                transcript,
                config.reference_fasta.as_deref(),
                config.hgvs,
            )
        } else {
            crate::hgvs::HgvsNotation::default()
        };
        return Some(TranscriptConsequence {
            consequences: smallvec![Consequence::TranscriptAblation],
            impact: Impact::HIGH,
            hgvsc: hgvs.hgvsc,
            hgvsp: hgvs.hgvsp,
            hgvs_offset: hgvs.offset,
            hgvs_offset_of_hgvsc_alone: hgvs.offset_of_hgvsc_alone,
            hgvs_intron_offsets: hgvs_intron_offsets(hgvs.intron_offsets),
            refseq_offset: hgvs.refseq_offset,
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
    let mut splice_evaluations = SpliceEvaluations::new();
    {
        let mut ctx = ApplyPositionContext {
            variant,
            transcript,
            consequences: &mut consequences,
            fields: &mut fields,
            reference_fasta: config.reference_fasta.as_deref(),
            splice_evaluations: &mut splice_evaluations,
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
    // the evaluations made so far skip the introns those calls read.
    if !insertion {
        add_splice_for_overlapping_introns(
            &mut consequences,
            variant,
            transcript,
            &mut splice_evaluations,
        );
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
        add_splice_for_overlapping_introns_insertion(
            &mut consequences,
            variant,
            transcript,
            &mut splice_evaluations,
        );
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
    // transcript has none even where its input span overlapped it. A row
    // outside the transcript whose allele adds bases is asked all the same:
    // it has an HGVSc where the 3' rule places those bases on the transcript.
    let hgvs = if config.compute_hgvs && (within_feature || crate::hgvs::adds_bases(input)) {
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
        hgvs_offset_of_hgvsc_alone: hgvs.offset_of_hgvsc_alone,
        hgvs_intron_offsets: hgvs_intron_offsets(hgvs.intron_offsets),
        refseq_offset: hgvs.refseq_offset,
        ..TranscriptConsequence::for_transcript(transcript)
    })
}

/// The HGVSc intron offsets in the row's pair. An intron offset is a distance
/// inside one transcript, so both fit the narrower pair; a pair that did not
/// would be a defect.
fn hgvs_intron_offsets(offsets: Option<(i64, i64)>) -> Option<(i32, i32)> {
    offsets.and_then(|(start, end)| {
        let start = i32::try_from(start).ok()?;
        let end = i32::try_from(end).ok()?;
        Some((start, end))
    })
}

fn is_flanking_position(position: &TranscriptPosition) -> bool {
    matches!(
        position,
        TranscriptPosition::Upstream { .. } | TranscriptPosition::Downstream { .. }
    )
}

/// The intron splice evaluations made for one allele, each keyed by the
/// intron's bounds and the polypyrimidine suppression it ran under (`false`
/// for an insertion, whose evaluation reads the exon overlap itself). An
/// evaluation whose key is already present returns at once: its window flags
/// are a function of the key and the allele, every term it adds is pushed
/// unique, and the one term it removes, `splice_region_variant`, is removed
/// again by the transcript-wide suppression once every window has been read.
/// Two slots hold the endpoint evaluations without a heap allocation.
type SpliceEvaluations = SmallVec<[(u64, u64, bool); 2]>;

/// Shared context for position-based consequence application.
struct ApplyPositionContext<'a> {
    variant: &'a InputVariant,
    transcript: &'a Transcript,
    consequences: &'a mut ConsequenceList,
    fields: &'a mut PositionFields,
    reference_fasta: Option<&'a vep_fasta::IndexedFasta>,
    splice_evaluations: &'a mut SpliceEvaluations,
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
    let splice_evaluations = &mut *ctx.splice_evaluations;
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
                        // leave the CDS (`CDS_position = N-?`); its early return
                        // (`return 0 unless defined $bvfo->translation_start`, a span
                        // whose transcript-5' base is outside the CDS) is the
                        // predicate's own.
                        if proj.has_partial_terminal_codon {
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
                let translation_start = (analysis_cds_pos - 1) / 3 + 1;
                let in_incomplete_terminal = partial_codon(transcript, Some(translation_start));

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
            let lo = variant.start.min(variant.end);
            let hi = variant.start.max(variant.end);
            // Every endpoint unmapped inside an exon's genomic span: the allele lies
            // in bases the model's alignment deleted, which carry no transcript
            // coordinate, so it takes the term of the region whose span holds it,
            // insertion or not, and reads no intron. An allele reaching one mapped
            // base takes its terms through that endpoint's arm instead.
            if position.is_unmapped()
                && paired_position.is_none_or(TranscriptPosition::is_unmapped)
                && allele_overlaps_any_exon(transcript, lo, hi)
            {
                push_unmapped_region_terms(consequences, variant, transcript);
                return;
            }

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
                        if !insertion && allele_overlaps_any_exon(transcript, lo, hi) {
                            push_unique(consequences, Consequence::NonCodingTranscriptExonVariant);
                        } else if transcript.facts().is_non_coding_context {
                            push_unique(consequences, Consequence::NonCodingTranscriptVariant);
                        }
                    }
                    return;
                }
            }

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
                        (intron_start, intron_end),
                        suppress_polypyr,
                        splice_evaluations,
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
                        (intron_start, intron_end),
                        overlaps_exon || fs_stretch,
                        splice_evaluations,
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
                        (intron_start, intron_end),
                        splice_evaluations,
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

/// The term of the transcript region whose genomic span holds an allele with no
/// transcript coordinate, every base of it one the model's alignment deleted:
/// on a coding model, `coding_sequence_variant` where the span overlaps the
/// coding region (Perl's `coding` pre-predicate) and the UTR terms of
/// [`push_perl_utr_terms`], so a span across the coding boundary takes both, as
/// a mapped span does; on any other model the exon term (a mature miRNA product
/// lies on mapped bases, so none holds the span). The bases have no cDNA, CDS or
/// protein position, so no field is written.
fn push_unmapped_region_terms(
    consequences: &mut ConsequenceList,
    variant: &InputVariant,
    transcript: &Transcript,
) {
    if transcript.facts().has_coding_model {
        if allele_overlaps_coding_bounds(variant, transcript) {
            push_unique(consequences, Consequence::CodingSequenceVariant);
        }
        push_perl_utr_terms(consequences, variant, transcript);
    } else {
        push_unique(consequences, Consequence::NonCodingTranscriptExonVariant);
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
/// `BaseTranscriptVariationAllele::_intron_effects`). Evaluated once per
/// intron per allele through `evaluated`.
fn add_intronic_insertion_splice_consequences(
    consequences: &mut ConsequenceList,
    transcript: &Transcript,
    variant: &InputVariant,
    (intron_start, intron_end): (u64, u64),
    evaluated: &mut SpliceEvaluations,
) {
    let key = (intron_start, intron_end, false);
    if evaluated.contains(&key) {
        return;
    }
    evaluated.push(key);
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
/// This is called after the per-endpoint `apply_position` calls; `evaluated`
/// skips the endpoint introns they evaluated.
fn add_splice_for_overlapping_introns(
    consequences: &mut ConsequenceList,
    variant: &InputVariant,
    transcript: &Transcript,
    evaluated: &mut SpliceEvaluations,
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
            (intron.start, intron.end),
            suppress_polypyr,
            evaluated,
        );
    }
}

/// Supplemental splice check for insertions, over every intron Perl's memoised
/// intron list holds for the anchor pair (`perl_intron_list_membership`), after
/// `within_feature` on the raw anchors (both inside the transcript). The
/// per-endpoint `apply_position` calls cover only the anchor introns, which
/// `evaluated` skips here.
fn add_splice_for_overlapping_introns_insertion(
    consequences: &mut ConsequenceList,
    variant: &InputVariant,
    transcript: &Transcript,
    evaluated: &mut SpliceEvaluations,
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
            (intron.start, intron.end),
            evaluated,
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

/// The windows of the donor, acceptor, donor 5th base, donor region,
/// polypyrimidine tract and intron terms on one intron, in that order, for a
/// transcript on `strand`: the five splice windows as `_intron_effects` places
/// them (BaseTranscriptVariationAllele.pm:155-204) and the predicates pick them
/// by strand (VariationEffect.pm:560-614), and for the intron term the whole
/// intron, its two end bases included.
fn strand_splice_windows(
    strand: vep_core::coordinate::Strand,
    intron_start: u64,
    intron_end: u64,
) -> [(u64, u64); 6] {
    let first_two = (intron_start, intron_start.saturating_add(1));
    let last_two = (intron_end.saturating_sub(1), intron_end);
    let whole = (intron_start, intron_end);
    match strand {
        vep_core::coordinate::Strand::Forward => [
            first_two,
            last_two,
            (
                intron_start.saturating_add(4),
                intron_start.saturating_add(4),
            ),
            (
                intron_start.saturating_add(2),
                intron_start.saturating_add(5),
            ),
            (intron_end.saturating_sub(16), intron_end.saturating_sub(2)),
            whole,
        ],
        vep_core::coordinate::Strand::Reverse => [
            last_two,
            first_two,
            (intron_end.saturating_sub(4), intron_end.saturating_sub(4)),
            (intron_end.saturating_sub(5), intron_end.saturating_sub(2)),
            (
                intron_start.saturating_add(2),
                intron_start.saturating_add(16),
            ),
            whole,
        ],
    }
}

/// For each of the six terms of `strand_splice_windows`, whether the span
/// `vf_lo..=vf_hi` holds a base of a window of that term on `intron` or on any
/// other intron of the transcript.
fn span_holds_splice_window_base(
    transcript: &Transcript,
    intron: (u64, u64),
    vf_lo: u64,
    vf_hi: u64,
) -> [bool; 6] {
    let introns = match transcript.vefc.as_ref() {
        Some(vefc) if !vefc.introns.is_empty() => vefc.introns.as_slice(),
        _ => transcript.introns.as_slice(),
    };
    let mut held = [false; 6];
    let spans = std::iter::once(intron).chain(introns.iter().map(|i| (i.start, i.end)));
    for (intron_start, intron_end) in spans {
        let windows = strand_splice_windows(transcript.strand, intron_start, intron_end);
        for (slot, (window_start, window_end)) in held.iter_mut().zip(windows) {
            *slot |= vf_hi >= window_start && vf_lo <= window_end;
        }
    }
    held
}

/// Splice terms for one intron of a non-insertion variant, following Perl's
/// `BaseTranscriptVariationAllele::_intron_effects`: every window is
/// tested per differing region (`_get_differing_regions`, which for a variant
/// with either allele of 1 base or less is the raw reference span and otherwise
/// the XOR of the alleles, as long as the longer allele), with the
/// non-normalising `overlap`. The intron-loop flags (intron_variant,
/// polypyrimidine) need the intron in the memoised intron list; the
/// boundary-loop flags (donor, acceptor, 5th base, donor region, splice_region,
/// polypyrimidine again) need it in the boundary list. Evaluated once per
/// intron and suppression per allele through `evaluated`.
///
/// One difference from `_intron_effects`: a donor, acceptor, 5th-base,
/// donor-region, polypyrimidine or intron term that only positions past the
/// variant's span reach is written only when the span holds a base of a window
/// of that term in the transcript (`span_holds_splice_window_base`); the
/// intron term's window is the whole intron.
fn add_intronic_non_insertion_splice_consequences(
    consequences: &mut ConsequenceList,
    transcript: &Transcript,
    variant: &InputVariant,
    (intron_start, intron_end): (u64, u64),
    suppress_polypyr: bool,
    evaluated: &mut SpliceEvaluations,
) {
    let key = (intron_start, intron_end, suppress_polypyr);
    if evaluated.contains(&key) {
        return;
    }
    evaluated.push(key);
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
    let mut within_core_intron = in_intron_list && any_region_overlap(core_start, core_end);
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

    let (
        mut donor_splice_site,
        mut acceptor_splice_site,
        mut donor_5th_base,
        mut donor_region,
        mut polypyrimidine,
    ) = match transcript.strand {
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

    // A region of an alternate allele longer than the reference allele runs
    // past the variant's span, onto bases the variant does not replace. Each
    // of these six terms then stands only where the span itself holds a base
    // of a window of that term on an intron of the transcript, which a region
    // inside the span always does; Ensembl VEP writes the term from the
    // positions past the span alone.
    let any_of_the_six = donor_splice_site
        || acceptor_splice_site
        || donor_5th_base
        || donor_region
        || polypyrimidine
        || within_core_intron;
    if any_of_the_six && regions.iter().any(|&(_, region_end)| region_end > vf_hi) {
        let [donor_held, acceptor_held, fifth_base_held, donor_region_held, tract_held, intron_held] =
            span_holds_splice_window_base(transcript, (intron_start, intron_end), vf_lo, vf_hi);
        donor_splice_site &= donor_held;
        acceptor_splice_site &= acceptor_held;
        donor_5th_base &= fifth_base_held;
        donor_region &= donor_region_held;
        polypyrimidine &= tract_held;
        within_core_intron &= intron_held;
    }

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
    shift_indel_3prime_within(variant, strand, fasta, true)
}

/// The span an insertion reaches when every base it passes matches its bases
/// cyclically, however far that is: its most 3' placement on `transcript`'s
/// strand, which the HGVS 3' rule describes it at. `perform_shift` reads 1,000
/// bases of flank (TranscriptVariationAllele.pm:427-437) and stops where they
/// end (:303-306), short of that placement in a longer run. `None` for any
/// other class and for an insertion nothing moves.
pub fn shift_insertion_to_its_last_placement(
    variant: &InputVariant,
    transcript: &Transcript,
    fasta: &vep_fasta::IndexedFasta,
) -> Option<(u64, u64)> {
    if variant.variant_class != vep_core::variant::VariantClass::Insertion {
        return None;
    }
    shift_indel_3prime_within(variant, transcript.strand, fasta, false)
}

/// [`shift_indel_3prime`], within [`perform_shift_bound`] when `bounded` and
/// to the end of the run otherwise.
fn shift_indel_3prime_within(
    variant: &InputVariant,
    strand: vep_core::coordinate::Strand,
    fasta: &vep_fasta::IndexedFasta,
    bounded: bool,
) -> Option<(u64, u64)> {
    use vep_core::coordinate::Strand;
    use vep_core::variant::VariantClass;

    let (mut start, mut end) = (variant.start, variant.end);
    let mut motif = shift_motif(variant)?;
    let is_insertion = variant.variant_class == VariantClass::Insertion;

    // VEP shifts indels in the transcript 3' direction (genomic right for + strand,
    // genomic left for - strand transcripts).
    let shift_right = strand == Strand::Forward;
    let max_shift = if bounded {
        perform_shift_bound(motif.len(), !shift_right)
    } else {
        u64::MAX
    };
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
/// `vf_start + region->{e}` without clamping, which extends the region onto
/// the bases after the span, bases the variant does not replace.
///
/// Use this variant when the caller needs Perl-identical per-region
/// coordinates. `add_intronic_non_insertion_splice_consequences` reads the
/// intronic windows with it and then holds the donor, acceptor, 5th-base,
/// donor-region, polypyrimidine and intron terms to the span;
/// `add_exonic_splice_region` reads its frameshift-intron skip with it. The
/// clamping version (`get_differing_regions`) is the one
/// `add_exonic_splice_region` tests its windows with, which relies on inverted
/// regions to signal "no overlap" against the exonic splice windows.
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
        .filter(|a| &*a.code == "miRNA")
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
    use vep_core::coordinate::Strand;

    /// [`super::calculate_consequences`] on a transcript the test owns: the row
    /// shares a copy of it, as the index shares its handle.
    fn calculate_consequences(
        variant: &InputVariant,
        transcript: &Transcript,
        config: &EffectsConfig,
    ) -> Option<TranscriptConsequence> {
        super::calculate_consequences(variant, &Arc::new(transcript.clone()), config)
    }

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
            code: "ncRNA".into(),
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
        assert_eq!(tc.transcript_id(), "ENST00000000001");
        assert_eq!(tc.gene_id(), "ENSG00000000001");
        assert_eq!(tc.gene_symbol(), Some("TEST1"));
        assert_eq!(tc.biotype(), Some("protein_coding"));
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
        assert_eq!(off.protein_id().as_deref(), Some("ENSP00000286808"));
        let on = EffectsConfig {
            protein_version: true,
            ..Default::default()
        };
        let tc = calculate_consequences(&variant, &tx, &on).unwrap();
        assert_eq!(tc.protein_id().as_deref(), Some("ENSP00000286808.3"));
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
        assert_eq!(tc.protein_id().as_deref(), Some("NP_001234.2"));

        let mut tx = make_test_transcript();
        tx.translation.as_mut().unwrap().version = None;
        let tc = calculate_consequences(&variant, &tx, &on).unwrap();
        assert_eq!(tc.protein_id().as_deref(), Some("ENSP00000000001"));
    }

    /// A breakend row against a transcript without a coding region prints no
    /// protein id, and `--protein_version` has nothing to append to.
    #[test]
    fn protein_version_leaves_a_row_without_a_protein_id_alone() {
        let mut breakend = InputVariant::new(
            "21".into(),
            25_001_000,
            25_001_000,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        breakend.variant_class = vep_core::variant::VariantClass::Translocation;
        breakend.is_structural = true;
        breakend.mate_id = Some("bnd_mate".into());
        let on = EffectsConfig {
            protein_version: true,
            ..Default::default()
        };
        let mut tx = make_test_transcript();
        tx.translation.as_mut().unwrap().version = Some(3);
        tx.coding_region_start = None;
        let tc = calculate_consequences(&breakend, &tx, &on).unwrap();
        assert_eq!(tc.protein_id(), None);
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

    /// A `T` inserted on the base before a non-coding transcript's first base,
    /// with n.1-14 `T` and n.15 `A` (GRCh38 21:25928745-25928775 around
    /// ENST00000617755), is `upstream_gene_variant` under `--hgvs` alone, its
    /// row carrying the notation and offset the 3' rule gives the inserted
    /// base and no HGVSp; under `--shift_3prime 1` it is annotated 14 bases in,
    /// exonic, with the same notation and offset.
    #[test]
    fn an_insertion_shifted_into_exon_1_carries_hgvsc_and_the_offset() {
        use crate::concordance_tests::{
            bases_from, hgvsc_fasta_fixture, make_compact_test_transcript, without_coding_model,
        };
        let tx = without_coding_model(make_compact_test_transcript());
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(91, "TTTGATGCCTTTTTTTTTTTTTTAACTGGAG"));
        let fasta = Arc::new(fasta);
        let config = |shift_3prime| EffectsConfig {
            reference_fasta: Some(fasta.clone()),
            enable_indel_3prime_shift: shift_3prime,
            compute_hgvs: true,
            hgvs: crate::hgvs::HgvsOptions {
                shift_3prime,
                ..Default::default()
            },
            ..Default::default()
        };
        let insertion = InputVariant::new("1".into(), 100, 99, b"-".to_vec(), b"T".to_vec());

        let tc = calculate_consequences(&insertion, &tx, &config(false)).unwrap();
        assert_eq!(
            tc.consequences.as_slice(),
            [Consequence::UpstreamGeneVariant]
        );
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:n.14dup"));
        assert_eq!(tc.hgvs_offset, Some(14));
        assert_eq!(tc.hgvsp, None);
        assert_eq!(tc.hgvs_intron_offsets, Some((0, 0)));

        let tc = calculate_consequences(&insertion, &tx, &config(true)).unwrap();
        assert_eq!(
            tc.consequences.as_slice(),
            [Consequence::NonCodingTranscriptExonVariant]
        );
        assert_eq!(tc.shift_length, 14);
        assert_eq!(tc.cdna_position.as_deref(), Some("14-15"));
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:n.14dup"));
        assert_eq!(tc.hgvs_offset, Some(14));
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

    /// The row carries the intron offsets of its HGVSc coordinates whenever it
    /// carries the notation, and neither without `compute_hgvs`.
    #[test]
    fn hgvs_intron_offsets_travel_with_the_hgvsc_on_the_row() {
        let tx = make_test_transcript();
        // Three bases into intron 1, c.250+3.
        let variant = InputVariant::new(
            "21".into(),
            25_000_302,
            25_000_302,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let with_hgvs = EffectsConfig {
            compute_hgvs: true,
            ..Default::default()
        };
        let tc = calculate_consequences(&variant, &tx, &with_hgvs).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:c.250+3A>G"));
        assert_eq!(tc.hgvs_intron_offsets, Some((3, 3)));

        let tc = calculate_consequences(&variant, &tx, &EffectsConfig::default()).unwrap();
        assert_eq!(tc.hgvsc, None);
        assert_eq!(tc.hgvs_intron_offsets, None);
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

    /// Ensembl VEP 116.2 writes HGVSc on the `transcript_ablation` row of a
    /// deletion covering the whole transcript: the row is `within_feature`
    /// (OutputFactory.pm:1721), `_var2transcript_slice_coords` clamps both ends
    /// onto the transcript (TranscriptVariationAllele.pm:2729-2733) and, at a
    /// zero 3' shift, `hgvs_transcript` describes the deletion from the first
    /// base to the last (:1472-1478, :1524-1525); `hgvs_protein` finds no
    /// peptide position for it (:1730-1738), so the row has no HGVSp and no
    /// `HGVS_OFFSET` (OutputFactory.pm:1755). Without `compute_hgvs` the row
    /// carries nothing.
    #[test]
    fn transcript_ablation_row_carries_the_hgvsc_clamped_to_the_whole_transcript() {
        let tx = make_test_transcript();
        // Ten bases before the transcript to ten bases past it.
        let variant = InputVariant::new(
            "21".into(),
            24_999_990,
            25_006_010,
            vec![b'A'; 6021],
            b"-".to_vec(),
        );
        let with_hgvs = EffectsConfig {
            compute_hgvs: true,
            ..Default::default()
        };
        let tc = calculate_consequences(&variant, &tx, &with_hgvs).unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::TranscriptAblation]
        );
        assert_eq!(
            tc.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.-50_*1701del")
        );
        assert_eq!(tc.hgvsp, None);
        assert_eq!(tc.hgvs_offset, None);
        assert_eq!(tc.hgvs_intron_offsets, Some((0, 0)));

        let tc = calculate_consequences(&variant, &tx, &EffectsConfig::default()).unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::TranscriptAblation]
        );
        assert_eq!(tc.hgvsc, None);
        assert_eq!(tc.hgvs_intron_offsets, None);
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

    /// A non-coding reverse-strand transcript with two exons, 1201-1300 (cDNA
    /// 1-100, cDNA 1 at genomic 1300) and 1001-1100 (cDNA 101-200), no
    /// translation and no cached sequence.
    fn non_coding_reverse_two_exon_transcript() -> Transcript {
        use crate::concordance_tests::make_compact_test_transcript;
        use vep_core::transcript::{Exon, ExonCoordMapper, MapperPair, TranscriptMapper};
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
        tx
    }

    /// A non-coding transcript has no cached cDNA, so its reference allele is
    /// read from the FASTA over the exon piece the span maps to, as Ensembl's
    /// `Transcript::spliced_seq` reads every exon from the genome: on a
    /// reverse-strand transcript the bases come back on the forward strand (the
    /// strand the variant is analysed on), and a span reaching into an intron
    /// maps to a gap, so the input's allele stands.
    #[test]
    fn used_ref_of_a_non_coding_reverse_transcript_is_read_from_the_fasta() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        // The FASTA is `T` but for 1300 `A`, 1299 `C`, 1100 `G`.
        let tx = non_coding_reverse_two_exon_transcript();
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

    /// A BAM-edited `ok` transcript stores its edited cDNA whole, and
    /// `TranscriptVariation::_get_ref_allele` reads that before anything else
    /// (`$tr->{_variation_effect_feature_cache}->{spliced_seq} ||=
    /// $tr->spliced_seq`, TranscriptVariation.pm:462), so a non-coding edited
    /// model's reference allele is the cached base where the FASTA disagrees,
    /// and needs no FASTA at all; the row carries the edit status.
    #[test]
    fn used_ref_of_a_bam_edited_non_coding_transcript_is_the_cached_cdna() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        use vep_core::transcript::BamEditStatus;
        // The transcript of the FASTA case above, with the edited cDNA stored:
        // its first three bases (genomic 1300 down to 1298) read `GCC` where
        // the FASTA reads `A`, `C`, `T`.
        let mut tx = non_coding_reverse_two_exon_transcript();
        tx.stable_id = "NR_000001.1".into();
        tx.bam_edit_status = Some(BamEditStatus::Ok);
        tx.vefc.as_mut().unwrap().spliced_seq = Some(format!("GCC{}", "T".repeat(197)));
        assert!(!tx.facts().has_coding_model);
        assert_eq!(tx.spliced_seq().map(|s| &s[..3]), Some("GCC"));
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(1300, b'A'), (1299, b'C'), (1100, b'G')]);

        // cDNA 1-3 is `GCC` on the transcript strand, `GGC` on the forward strand.
        let three = InputVariant::new("1".into(), 1298, 1300, b"TCA".to_vec(), b"-".to_vec());
        assert_eq!(
            transcript_reference_allele(&three, &tx, Some(&fasta)),
            Some(b"GGC".to_vec())
        );
        assert_eq!(
            transcript_reference_allele(&three, &tx, None),
            Some(b"GGC".to_vec())
        );
        let first_base = InputVariant::new("1".into(), 1300, 1300, b"A".to_vec(), b"T".to_vec());
        assert_eq!(
            transcript_reference_allele(&first_base, &tx, None),
            Some(b"C".to_vec())
        );

        let config = EffectsConfig {
            reference_fasta: None,
            use_transcript_ref: true,
            ..Default::default()
        };
        let tc = calculate_consequences(&three, &tx, &config).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("GGC"));
        assert_eq!(tc.bam_edit_status, Some(BamEditStatus::Ok));
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
    /// allele rotated as `shift_feature_seqs` rotates it, read on the
    /// transcript's strand: one place per base of the shift, so `AAG` moved one
    /// base on the reverse-strand APP prints `GAA` as Ensembl VEP 116.2 does,
    /// the same insertion given on the reverse strand as `CTT` prints `TTC`
    /// (Ensembl VEP 116.2 turns it by the transcript's strand alone to `TCT`,
    /// TranscriptVariationAllele.pm:1337-1338, a defect), a shift of the
    /// allele's length or of a multiple leaves it, and a longer shift rotates by
    /// the remainder (Ensembl VEP 116.2 leaves a reverse-strand allele as
    /// written there, `:1338-1339`, a defect). A row with no translation
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
        assert_eq!(
            rotated(&forward_given, Strand::Reverse, span, 4).as_deref(),
            Some("GAA")
        );
        assert_eq!(rotated(&forward_given, Strand::Reverse, span, 6), None);
        assert_eq!(
            rotated(&insertion("CTT", Strand::Reverse), Strand::Reverse, span, 1).as_deref(),
            Some("TTC")
        );
        assert_eq!(rotated(&forward_given, Strand::Reverse, None, 1), None);
        assert_eq!(rotated(&forward_given, Strand::Reverse, span, 0), None);
        let mut reference_row = insertion("-", Strand::Forward);
        reference_row.reference_allele_row = true;
        assert_eq!(rotated(&reference_row, Strand::Reverse, span, 1), None);
        let deletion = InputVariant::new("21".into(), 101, 103, b"TCT".to_vec(), b"-".to_vec());
        assert_eq!(rotated(&deletion, Strand::Reverse, span, 1), None);
    }

    /// The path a shifted deletion's `USED_REF` rotates through: the protein
    /// path on a row with a protein notation and a non-zero HGVS offset whose
    /// allele as annotated passes `coding`, the
    /// codon path on a row the consequence shift moved whose span has
    /// translation coordinates, unless the reference row read the codon first,
    /// and neither otherwise. Both paths rotate the deleted bases the same way,
    /// so `GA` deleted at chromosome 21 30348449 and moved three bases on
    /// ENST00000334160 (reverse strand) prints `USED_REF=AG` under
    /// `--shift_3prime 1` alone as under `--hgvs`, the bases the shifted span
    /// holds; Ensembl VEP 116.2 prints `GA` under `--shift_3prime 1` alone,
    /// where `shift_feature_seqs` counts `seq_length - shift_length` places
    /// (TranscriptVariationAllele.pm:1338-1339), a defect.
    #[test]
    fn used_ref_rotation_follows_the_path_and_both_paths_rotate_alike() {
        use vep_core::coordinate::Strand;
        let mut deleted = b"GA".to_vec();
        crate::coding::rotate_shifted_allele(&mut deleted, 3, Strand::Forward, Strand::Reverse);
        assert_eq!(deleted, b"AG");

        let mut tc = TranscriptConsequence {
            shift_length: 3,
            ..Default::default()
        };
        assert_eq!(
            used_ref_rotation(&tc, Some((2, 3)), false, true),
            Some(UsedRefRotation::Codon(3))
        );
        assert_eq!(used_ref_rotation(&tc, Some((2, 3)), true, true), None);
        assert_eq!(used_ref_rotation(&tc, None, false, true), None);
        tc.hgvs_offset = Some(-3);
        assert_eq!(
            used_ref_rotation(&tc, Some((2, 3)), false, true),
            Some(UsedRefRotation::Codon(3))
        );
        tc.hgvsp = Some("ENSP00000000001.1:p.Leu2fs".into());
        assert_eq!(
            used_ref_rotation(&tc, Some((2, 3)), false, true),
            Some(UsedRefRotation::Protein(3))
        );
        assert_eq!(
            used_ref_rotation(&tc, Some((2, 3)), true, true),
            Some(UsedRefRotation::Protein(3))
        );
        tc.shift_length = 0;
        assert_eq!(
            used_ref_rotation(&tc, None, false, true),
            Some(UsedRefRotation::Protein(3))
        );
        // A notation on a row whose reference allele Ensembl VEP does not
        // rewrite: an allele annotated in an intron.
        assert_eq!(used_ref_rotation(&tc, None, false, false), None);
        tc.hgvs_offset = None;
        assert_eq!(used_ref_rotation(&tc, Some((2, 3)), false, true), None);
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
    /// prints `AAG`. A shift longer than the deletion rotates by the remainder
    /// through either path: `GA` moved three bases prints `AG`, the bases of the
    /// shifted span, under `--shift_3prime 1` alone as under `--hgvs` (Ensembl
    /// VEP 116.2 prints `GA` under `--shift_3prime 1` alone, where
    /// `shift_feature_seqs` counts `seq_length - shift_length` places,
    /// TranscriptVariationAllele.pm:1338-1339, a defect).
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
        assert_eq!(tc.used_ref.as_deref(), Some("AG"));
        let tc = calculate_consequences(&three, &tx, &config(false, true)).unwrap();
        assert_eq!(tc.hgvs_offset, Some(-3));
        assert!(tc.hgvsp.is_some());
        assert_eq!(tc.used_ref.as_deref(), Some("AG"));
        let tc = calculate_consequences(&three, &tx, &config(true, true)).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("AG"));
        let misspelt = InputVariant::new("1".into(), 1011, 1012, b"TA".to_vec(), b"-".to_vec());
        let tc = calculate_consequences(&misspelt, &tx, &config(true, false)).unwrap();
        assert_eq!(tc.shift_length, 3);
        assert_eq!(tc.used_ref.as_deref(), Some("AG"));

        // Given on the reverse strand, `AGA` is `TCT`: the codon path rotates the
        // allele as written, read on the transcript's strand, and the protein
        // path hands back the forward-strand rotation, so both print the shifted
        // span's `AAG` on the record's strand, `CTT` (Ensembl VEP's codon path
        // turns the written `TCT` by the transcript's strand alone,
        // `seq_length - shift_length` places, TranscriptVariationAllele.pm:1337-1338,
        // which gives `TTC`, a defect).
        let mut given_on_reverse = one.clone();
        given_on_reverse.strand = vep_core::coordinate::Strand::Reverse;
        given_on_reverse.ref_allele = b"TCT".to_vec();
        given_on_reverse.allele_string = "TCT/-".into();
        let tc = calculate_consequences(&given_on_reverse, &tx, &config(true, false)).unwrap();
        assert_eq!(tc.shift_length, 1);
        assert_eq!(tc.used_ref.as_deref(), Some("CTT"));
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
    /// covers and leaves out: a BAM-edited cache switches the flag on and
    /// supplies the edited cDNA, an edited RefSeq transcript's mismatch row is
    /// marked for the warning; run-time editing (`--bam`) is not ported.
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
            "BAM-edited cache switches the flag on by itself",
            "`VEP::CacheDir`",
            "cut from the edited cDNA the cache stores",
            "`apply_edits`",
            "`--bam`",
            "`Transcript-assembly mismatch` warning",
        ] {
            assert!(module_docs.contains(phrase), "{phrase}: {module_docs}");
        }
    }

    /// One mapper pair in cDNA order on the forward strand.
    fn pair(
        from_start: u64,
        from_end: u64,
        to_start: u64,
        to_end: u64,
    ) -> vep_core::transcript::MapperPair {
        vep_core::transcript::MapperPair {
            from_start,
            from_end,
            to_start,
            to_end,
            ori: 1,
        }
    }

    /// The compact coding transcript as a RefSeq model named `stable_id` whose
    /// alignment inserted `AAA` after unedited cDNA 9 (`_rna_edit` `10 9 AAA`),
    /// as the cache holds such a model: the mapper folds the edit (cDNA 1-9 and
    /// 13-303 on exon 1, 304-603 and 604-2604 on the other two), the coding
    /// start moves from 51 to 54, and the cached 5' UTR carries the three bases
    /// (the cache loader's `check_mapper_bounds` reads the attribute the same
    /// way).
    fn insertion_edited_model(stable_id: &str) -> Transcript {
        use vep_core::transcript::{Attribute, ExonCoordMapper};
        let mut tx = compact_transcript_with_cached_cdna();
        tx.stable_id = stable_id.into();
        tx.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: "10 9 AAA".to_string(),
        });
        tx.attributes.push(Attribute {
            code: "rseq_mrna_match".into(),
            value: String::new(),
        });
        tx.cdna_coding_start = Some(54);
        tx.cdna_coding_end = Some(903);
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.five_prime_utr = Some(format!("{}AAA{}", "C".repeat(9), "C".repeat(41)));
        let mapper = vefc.mapper.as_mut().unwrap();
        mapper.cdna_coding_start = 54;
        mapper.cdna_coding_end = 903;
        mapper.exon_coord_mapper = ExonCoordMapper::new(vec![
            pair(1, 9, 100, 108),
            pair(13, 303, 109, 399),
            pair(304, 603, 2100, 2399),
            pair(604, 2604, 4100, 6100),
        ]);
        tx.derived = Default::default();
        assert_eq!(tx.spliced_seq().map(|s| &s[9..12]), Some("AAA"));
        assert_eq!(tx.spliced_seq().map(|s| &s[53..56]), Some("ATG"));
        tx
    }

    /// The compact coding transcript as an `NM_` model with the three edits of
    /// `get_misalignment_offset`'s reading: `10 9 AAA` inserts three bases,
    /// `50 52 ` deletes three, `70 70 T` replaces one base with one. Its 5' UTR
    /// is 100 bases long so every edit sits upstream of the coding start, which
    /// the edits leave at 101; the mapper folds the first two (cDNA 1-9, 13-52
    /// and 53-300 on exon 1, the deleted bases at 149-151 mapping nowhere).
    fn three_edit_model() -> Transcript {
        use vep_core::transcript::{Attribute, ExonCoordMapper};
        let mut tx = compact_transcript_with_cached_cdna();
        tx.stable_id = "NM_000003.1".into();
        for value in ["10 9 AAA", "50 52 ", "70 70 T"] {
            tx.attributes.push(Attribute {
                code: "_rna_edit".into(),
                value: value.to_string(),
            });
        }
        tx.cdna_coding_start = Some(101);
        tx.cdna_coding_end = Some(950);
        tx.coding_region_start = Some(200);
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.five_prime_utr = Some(format!(
            "{}AAA{}{}",
            "C".repeat(9),
            "C".repeat(40),
            "C".repeat(48)
        ));
        let mapper = vefc.mapper.as_mut().unwrap();
        mapper.cdna_coding_start = 101;
        mapper.cdna_coding_end = 950;
        mapper.exon_coord_mapper = ExonCoordMapper::new(vec![
            pair(1, 9, 100, 108),
            pair(13, 52, 109, 148),
            pair(53, 300, 152, 399),
            pair(301, 600, 2100, 2399),
            pair(601, 2601, 4100, 6100),
        ]);
        tx.derived = Default::default();
        tx
    }

    /// An effects configuration computing HGVS against `fasta` without the 3'
    /// shift, which an all-`T` reference would carry to the transcript's end.
    fn hgvs_config(
        fasta: &Arc<vep_fasta::IndexedFasta>,
        use_transcript_ref: bool,
    ) -> EffectsConfig {
        EffectsConfig {
            reference_fasta: Some(Arc::clone(fasta)),
            compute_hgvs: true,
            use_transcript_ref,
            hgvs: crate::hgvs::HgvsOptions {
                shift_hgvs: false,
                ..Default::default()
            },
            ..Default::default()
        }
    }

    /// `get_misalignment_offset` (TranscriptVariationAllele.pm:1594-1624) sums
    /// the length changes of the indel edits whose start is at or before the
    /// allele's cDNA start (`:1608-1610`), the insertion `10 9 AAA` counting +3,
    /// the deletion `50 52 ` -3 and the substitution `70 70 T` nothing (`:1606`):
    /// an allele at cDNA 30 carries +3 as `REFSEQ_OFFSET`, one at cDNA 60 sums
    /// to 0 and carries no column (OutputFactory.pm:1769), one at cDNA 5 sums
    /// nothing. The offset is computed for an `NM_` or `XM_` model alone
    /// (`:1517`) and only where `hgvs_transcript` runs, under `--hgvs`. A
    /// genomic base the deletion edit removed from the model maps to no cDNA
    /// position, so its notation is described from the nearest modelled base
    /// as an intronic base is, `c.-49+2`; Ensembl VEP numbers it `c.-50`, as if
    /// the model held it.
    #[test]
    fn refseq_offset_sums_the_indel_edits_upstream_of_the_variant() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let tx = three_edit_model();
        let config = hgvs_config(&fasta, false);
        let at = |pos: u64| InputVariant::new("1".into(), pos, pos, b"T".to_vec(), b"A".to_vec());
        // Genomic 126 is folded cDNA 30, 159 is 60, 104 is 5; 150 is deleted.
        let at_30 = calculate_consequences(&at(126), &tx, &config).unwrap();
        assert_eq!(at_30.refseq_offset, Some(3));
        assert_eq!(at_30.hgvsc.as_deref(), Some("NM_000003.1:c.-71T>A"));
        let at_60 = calculate_consequences(&at(159), &tx, &config).unwrap();
        assert_eq!(at_60.refseq_offset, None);
        assert_eq!(at_60.hgvsc.as_deref(), Some("NM_000003.1:c.-41T>A"));
        let at_5 = calculate_consequences(&at(104), &tx, &config).unwrap();
        assert_eq!(at_5.refseq_offset, None);
        assert_eq!(at_5.hgvsc.as_deref(), Some("NM_000003.1:c.-96T>A"));
        assert_eq!(
            crate::hgvs::generate_hgvsc(&at(150), &tx, Some(&fasta)).as_deref(),
            Some("NM_000003.1:c.-49+2T>A")
        );
        assert_eq!(
            at_30.refseq_match().map(|codes| codes.len()),
            None,
            "no rseq code on this model"
        );
        let without_hgvs =
            calculate_consequences(&at(126), &tx, &EffectsConfig::default()).unwrap();
        assert_eq!(without_hgvs.refseq_offset, None);
        let mut ensembl = three_edit_model();
        ensembl.stable_id = "ENST00000000003".into();
        ensembl.derived = Default::default();
        assert_eq!(
            calculate_consequences(&at(126), &ensembl, &config)
                .unwrap()
                .refseq_offset,
            None
        );
    }

    /// HGVSc on a model with an indel edit is numbered against the edited model
    /// through the cached mapper pairs, which fold the edit: a 5' UTR base at
    /// folded cDNA 33 is `c.-21` on the `NM_` and `XM_` models, and an intronic
    /// allele's distance counts from the edited exon end, `c.250+601`, the
    /// reading of the unedited transcript too since the insertion sits before
    /// the coding start. Ensembl VEP walks the exons unedited and moves an
    /// `NM_` or `XM_` position by the misalignment offset
    /// (TranscriptVariationAllele.pm:1517, 1524-1525), which an intronic allele,
    /// having no cDNA start, never gets (`:1609`), so it reads `c.247+601`
    /// there; a deletion of the first exon's last two bases, which the offset
    /// pushes past the unedited exon end into `c.247+2_247+3del`, reads
    /// `c.249_250del` here. `REFSEQ_OFFSET` keeps Ensembl VEP's value.
    #[test]
    fn hgvsc_of_an_edited_model_is_numbered_against_the_edited_model() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, false);
        let utr = InputVariant::new("1".into(), 129, 129, b"T".to_vec(), b"A".to_vec());
        let intronic = InputVariant::new("1".into(), 1000, 1000, b"T".to_vec(), b"A".to_vec());
        let exon_end = InputVariant::new("1".into(), 398, 399, b"TT".to_vec(), b"-".to_vec());

        let nm = insertion_edited_model("NM_000001.1");
        let tc = calculate_consequences(&utr, &nm, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000001.1:c.-21T>A"));
        assert_eq!(tc.refseq_offset, Some(3));
        assert_eq!(
            tc.refseq_match().map(|codes| codes.to_vec()),
            Some(vec![Box::from("rseq_mrna_match")])
        );
        let tc = calculate_consequences(&intronic, &nm, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000001.1:c.250+601T>A"));
        assert_eq!(tc.refseq_offset, None);
        let tc = calculate_consequences(&exon_end, &nm, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000001.1:c.249_250del"));
        assert_eq!(tc.refseq_offset, Some(3));

        let xm = insertion_edited_model("XM_000001.1");
        assert_eq!(
            calculate_consequences(&utr, &xm, &config)
                .unwrap()
                .hgvsc
                .as_deref(),
            Some("XM_000001.1:c.-21T>A")
        );

        let unedited = compact_transcript_with_cached_cdna();
        assert_eq!(
            calculate_consequences(&intronic, &unedited, &config)
                .unwrap()
                .hgvsc
                .as_deref(),
            Some("ENST00000000001.1:c.250+601T>A")
        );
        assert_eq!(
            calculate_consequences(&exon_end, &unedited, &config)
                .unwrap()
                .hgvsc
                .as_deref(),
            Some("ENST00000000001.1:c.249_250del")
        );
    }

    /// An Ensembl model with an indel `_rna_edit` is numbered against the
    /// edited model like a RefSeq one: the 5' UTR base at folded cDNA 33 reads
    /// `c.-21`, where Ensembl VEP, which computes the misalignment offset for an
    /// `NM_` or `XM_` model alone (TranscriptVariationAllele.pm:1517), reads the
    /// unedited `c.-24`. `REFSEQ_OFFSET` stays the RefSeq models' column.
    #[test]
    fn ensembl_model_with_an_indel_edit_is_numbered_against_the_edited_model() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, false);
        let enst = insertion_edited_model("ENST00000000001");
        let utr = InputVariant::new("1".into(), 129, 129, b"T".to_vec(), b"A".to_vec());
        let tc = calculate_consequences(&utr, &enst, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:c.-21T>A"));
        assert_eq!(tc.refseq_offset, None);
        let intronic = InputVariant::new("1".into(), 1000, 1000, b"T".to_vec(), b"A".to_vec());
        let tc = calculate_consequences(&intronic, &enst, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:c.250+601T>A"));
    }

    /// The edited numbering holds on the reverse strand, where the mapper folds
    /// an insertion into the genomic end of the first exon: on the non-coding
    /// two-exon model edited by `21 20 AAA` (cDNA 1-20 at 1281-1300, 24-103 at
    /// 1201-1280, 104-203 at 1001-1100), the exonic base at 1250 is `n.54` and
    /// the intronic base at 1180, 21 past the first exon's genomic start, is
    /// `n.103+21`, where the unedited model reads `n.51` and `n.100+21`, the
    /// strings Ensembl VEP writes for the edited model too, since an `NR_`
    /// model gets no misalignment offset (TranscriptVariationAllele.pm:1517).
    /// The alleles are spelled on the transcript strand.
    #[test]
    fn reverse_strand_edited_model_is_numbered_against_the_edited_model() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        use vep_core::transcript::{Attribute, ExonCoordMapper, MapperPair};
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, false);
        let exonic = InputVariant::new("1".into(), 1250, 1250, b"T".to_vec(), b"A".to_vec());
        let intronic = InputVariant::new("1".into(), 1180, 1180, b"T".to_vec(), b"A".to_vec());

        let unedited = non_coding_reverse_two_exon_transcript();
        let tc = calculate_consequences(&exonic, &unedited, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:n.51A>T"));
        let tc = calculate_consequences(&intronic, &unedited, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:n.100+21A>T"));

        let mut edited = non_coding_reverse_two_exon_transcript();
        edited.stable_id = "NR_000001.1".into();
        edited.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: "21 20 AAA".to_string(),
        });
        let pair = |from_start, from_end, to_start, to_end| MapperPair {
            from_start,
            from_end,
            to_start,
            to_end,
            ori: -1,
        };
        let mapper = edited.vefc.as_mut().unwrap().mapper.as_mut().unwrap();
        mapper.exon_coord_mapper = ExonCoordMapper::new(vec![
            pair(1, 20, 1281, 1300),
            pair(24, 103, 1201, 1280),
            pair(104, 203, 1001, 1100),
        ]);
        edited.derived = Default::default();
        let tc = calculate_consequences(&exonic, &edited, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("NR_000001.1:n.54A>T"));
        assert_eq!(tc.refseq_offset, None);
        let tc = calculate_consequences(&intronic, &edited, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("NR_000001.1:n.103+21A>T"));
        assert_eq!(tc.refseq_offset, None);
    }

    /// A single-base substitution in a coding exon takes the CDS coordinate the
    /// edited mapper gives (TranscriptVariationAllele.pm:1519-1522), and every
    /// other allele goes through the edited-model walk of `_get_cDNA_position`,
    /// which agrees with it inside an exon: on the model edited after cDNA 9, a
    /// substitution at genomic 159 is `c.10` and a deletion of 159-160 is
    /// `c.10_11del`, both carrying the offset 3.
    #[test]
    fn coding_exon_snp_hgvsc_takes_the_cds_position_unshifted() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, false);
        let nm = insertion_edited_model("NM_000001.1");
        let snp = InputVariant::new("1".into(), 159, 159, b"T".to_vec(), b"A".to_vec());
        let tc = calculate_consequences(&snp, &nm, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000001.1:c.10T>A"));
        assert_eq!(tc.cds_position.as_deref(), Some("10"));
        assert_eq!(tc.refseq_offset, Some(3));
        let deletion = InputVariant::new("1".into(), 159, 160, b"TT".to_vec(), b"-".to_vec());
        let tc = calculate_consequences(&deletion, &nm, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000001.1:c.10_11del"));
        assert_eq!(tc.cds_position.as_deref(), Some("10-11"));
        assert_eq!(tc.refseq_offset, Some(3));
    }

    /// On a transcript with an `_rna_edit` the notation spells the alleles
    /// relative to the transcript (`$ref->feature_seq`, `$self->feature_seq`,
    /// TranscriptVariationAllele.pm:1510-1514), the reference being the allele
    /// the row was analysed against, so under `--use_transcript_ref` a
    /// substitution at CDS 3 reads the transcript's `G` where the genome has
    /// `T`, and the alleles are clipped as on any other model, so a two-base
    /// allele sharing its first base with the reference prints the one base
    /// that changes. The unedited transcript keeps the reference the genomic
    /// slice gives.
    #[test]
    fn hgvsc_ref_and_alt_follow_the_transcript_sequence_when_edits_exist() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, true);
        let nm = insertion_edited_model("NM_000001.1");
        // Genomic 152 is folded cDNA 56, CDS 3, the `G` of `ATG`.
        let snp = InputVariant::new("1".into(), 152, 152, b"T".to_vec(), b"C".to_vec());
        let tc = calculate_consequences(&snp, &nm, &config).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("G"));
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000001.1:c.3G>C"));
        let mnp = InputVariant::new("1".into(), 151, 152, b"TT".to_vec(), b"TC".to_vec());
        let tc = calculate_consequences(&mnp, &nm, &config).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("TG"));
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000001.1:c.3G>C"));

        let unedited = compact_transcript_with_cached_cdna();
        let tc = calculate_consequences(&snp, &unedited, &config).unwrap();
        assert_eq!(tc.used_ref.as_deref(), Some("G"));
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:c.3T>C"));
        let tc = calculate_consequences(&mnp, &unedited, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("ENST00000000001.1:c.3T>C"));
    }

    /// `invalid_alleles` (TranscriptVariationAllele.pm:826-837) marks an
    /// alternate allele equal to the bases of the transcript's first
    /// `_rna_edit`, read on the codon path alone: an insertion of `AAA` in the
    /// coding sequence of the model edited by `10 9 AAA` is marked, an
    /// insertion of `AAC`, the same insertion in the 5' UTR (no translation
    /// coordinate) and the reference row are not, and a transcript without an
    /// edit marks nothing. The deletion of a model whose first edit deletes
    /// matches its empty bases, as the rule reads `-` as the empty sequence
    /// (`:822`).
    #[test]
    fn transcript_assembly_mismatch_follows_the_first_rna_edit() {
        use vep_core::transcript::Attribute;
        let nm = insertion_edited_model("NM_000001.1");
        let config = EffectsConfig::default();
        let coding = InputVariant::new("1".into(), 160, 159, b"-".to_vec(), b"AAA".to_vec());
        let tc = calculate_consequences(&coding, &nm, &config).unwrap();
        assert!(tc.transcript_assembly_mismatch);
        let other = InputVariant::new("1".into(), 160, 159, b"-".to_vec(), b"AAC".to_vec());
        assert!(
            !calculate_consequences(&other, &nm, &config)
                .unwrap()
                .transcript_assembly_mismatch
        );
        let utr = InputVariant::new("1".into(), 120, 119, b"-".to_vec(), b"AAA".to_vec());
        assert!(
            !calculate_consequences(&utr, &nm, &config)
                .unwrap()
                .transcript_assembly_mismatch
        );
        let reference_row = coding.reference_allele_row();
        assert!(
            !calculate_consequences(&reference_row, &nm, &config)
                .unwrap()
                .transcript_assembly_mismatch
        );
        let unedited = compact_transcript_with_cached_cdna();
        assert!(
            !calculate_consequences(&coding, &unedited, &config)
                .unwrap()
                .transcript_assembly_mismatch
        );

        let mut deleting = compact_transcript_with_cached_cdna();
        deleting.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: "901 903 ".to_string(),
        });
        deleting.derived = Default::default();
        let deletion = InputVariant::new("1".into(), 159, 160, b"TT".to_vec(), b"-".to_vec());
        assert!(
            calculate_consequences(&deletion, &deleting, &config)
                .unwrap()
                .transcript_assembly_mismatch
        );
        assert!(
            !calculate_consequences(&coding, &deleting, &config)
                .unwrap()
                .transcript_assembly_mismatch
        );
    }

    /// The compact coding transcript as a RefSeq model whose alignment inserted
    /// `ACC` after unedited cDNA 58 (`_rna_edit` `59 58 ACC`), inside the coding
    /// sequence, as the cache holds such a model: the mapper folds the edit
    /// (cDNA 1-58 and 62-303 on exon 1), so genomic 157 is CDS 8, genomic 158
    /// is CDS 12, and CDS 9-11, the bases the genome lacks, lie between them.
    /// The coding sequence reads `ATG GCT GGA CCA AAA`.
    fn coding_seam_model() -> Transcript {
        let tx = coding_model_holding(8, "ACC");
        assert_eq!(
            tx.spliced_seq().map(|s| &s[50..65]),
            Some("ATGGCTGGACCAAAA")
        );
        tx
    }

    /// [`coding_seam_model`] with `held` as the bases the alignment inserted
    /// after unedited cDNA 58: genomic 157 is CDS 8, the held bases follow from
    /// CDS 9, and genomic 158 is the CDS position after the last of them.
    fn coding_seam_model_holding(held: &str) -> Transcript {
        use vep_core::transcript::{Attribute, ExonCoordMapper};
        let n = held.len() as u64;
        let mut tx = compact_transcript_with_cached_cdna();
        tx.stable_id = "NM_000004.1".into();
        tx.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: format!("59 58 {held}"),
        });
        tx.cdna_coding_end = Some(900 + n);
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.translateable_seq.as_mut().unwrap().insert_str(8, held);
        let mapper = vefc.mapper.as_mut().unwrap();
        mapper.cdna_coding_end = 900 + n;
        mapper.exon_coord_mapper = ExonCoordMapper::new(vec![
            pair(1, 58, 100, 157),
            pair(59 + n, 300 + n, 158, 399),
            pair(301 + n, 600 + n, 2100, 2399),
            pair(601 + n, 2601 + n, 4100, 6100),
        ]);
        tx.derived = Default::default();
        assert_eq!(
            tx.spliced_seq().map(|s| &s[58..58 + held.len()]),
            Some(held)
        );
        tx
    }

    /// The compact coding transcript as a RefSeq model whose alignment
    /// inserted `bases` after CDS `after_cds`, in exon 1, as the cache holds
    /// such a model: the mapper folds the edit, so the genomic base of CDS
    /// `after_cds` and the next genomic base lie on either side of `bases`,
    /// which the genome lacks.
    fn coding_model_holding(after_cds: usize, bases: &str) -> Transcript {
        use vep_core::transcript::{Attribute, ExonCoordMapper};
        let held = bases.len() as u64;
        let cdna = after_cds as u64 + 50;
        let mut tx = compact_transcript_with_cached_cdna();
        tx.stable_id = "NM_000004.1".into();
        tx.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: format!("{} {cdna} {bases}", cdna + 1),
        });
        tx.cdna_coding_end = Some(900 + held);
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.translateable_seq
            .as_mut()
            .unwrap()
            .insert_str(after_cds, bases);
        let mapper = vefc.mapper.as_mut().unwrap();
        mapper.cdna_coding_end = 900 + held;
        mapper.exon_coord_mapper = ExonCoordMapper::new(vec![
            pair(1, cdna, 100, 99 + cdna),
            pair(cdna + held + 1, 300 + held, 100 + cdna, 399),
            pair(301 + held, 600 + held, 2100, 2399),
            pair(601 + held, 2601 + held, 4100, 6100),
        ]);
        tx.derived = Default::default();
        tx
    }

    /// [`coding_model_holding`] as a model whose alignment also dropped the
    /// genomic base `dropped`, further along exon 1 (an `_rna_edit` with no
    /// bases): the mapper splits the exon's pair once more around it, numbers
    /// every later base one lower, and the exon list keeps the one exon.
    fn coding_model_holding_and_dropping(
        after_cds: usize,
        bases: &str,
        dropped: u64,
    ) -> Transcript {
        use vep_core::transcript::{Attribute, ExonCoordMapper};
        let held = bases.len() as u64;
        let cdna = after_cds as u64 + 50;
        let dropped_cdna = dropped - 99 + held;
        let mut tx = coding_model_holding(after_cds, bases);
        tx.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: format!("{0} {0} ", dropped - 99),
        });
        tx.cdna_coding_end = Some(899 + held);
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.translateable_seq
            .as_mut()
            .unwrap()
            .remove(dropped_cdna as usize - 51);
        let mapper = vefc.mapper.as_mut().unwrap();
        mapper.cdna_coding_end = 899 + held;
        mapper.exon_coord_mapper = ExonCoordMapper::new(vec![
            pair(1, cdna, 100, 99 + cdna),
            pair(cdna + held + 1, dropped_cdna - 1, 100 + cdna, dropped - 1),
            pair(dropped_cdna, 299 + held, dropped + 1, 399),
            pair(300 + held, 599 + held, 2100, 2399),
            pair(600 + held, 2600 + held, 4100, 6100),
        ]);
        tx.derived = Default::default();
        tx
    }

    /// An insertion between the two genomic bases around the bases an edited
    /// model holds and the genome lacks maps onto those bases: its CDS span
    /// is theirs (`Bio::EnsEMBL::Mapper::map_insert` shrinks the two mapped
    /// bases towards each other), and the alternate coding sequence carries
    /// the allele in their place. An allele as long as that span changes the
    /// residues and keeps the length, a `missense_variant`. Ensembl VEP reads
    /// the length class from the genomic span, where every insertion is
    /// longer than its reference, and `missense_variant` excludes
    /// `increase_length` (Utils/Constants.pm), so the row earns no term there
    /// and takes the default `intergenic_variant`
    /// (BaseVariationFeatureOverlapAllele.pm:285). The terms both engines
    /// give on the same span stand: the model's own bases are
    /// `synonymous_variant`, a stop codon `stop_gained`, an allele of another
    /// length `frameshift_variant`, and the same three bases one position
    /// further on, off the edit, an `inframe_insertion`.
    #[test]
    fn insertion_as_long_as_the_edited_model_bases_it_stands_on_is_a_missense_variant() {
        let tx = coding_seam_model();
        let config = EffectsConfig::default();
        let on_the_seam =
            |alt: &[u8]| InputVariant::new("1".into(), 158, 157, b"-".to_vec(), alt.to_vec());
        let row = calculate_consequences(&on_the_seam(b"TTT"), &tx, &config)
            .expect("an insertion inside a coding exon has a row");
        assert_eq!(row.consequences.as_slice(), [Consequence::MissenseVariant]);
        assert_eq!(row.impact, Impact::MODERATE);
        assert_eq!(row.cds_position.as_deref(), Some("9-11"));
        assert_eq!(row.protein_position.as_deref(), Some("3-4"));
        assert_eq!(row.amino_acids.as_deref(), Some("GP/GL"));
        assert_eq!(row.codons.as_deref(), Some("ggacca/ggTTTa"));

        let terms = |variant: &InputVariant| {
            calculate_consequences(variant, &tx, &config)
                .unwrap()
                .consequences
                .to_vec()
        };
        assert_eq!(
            terms(&on_the_seam(b"ACC")),
            [Consequence::SynonymousVariant]
        );
        assert_eq!(terms(&on_the_seam(b"ATG")), [Consequence::StopGained]);
        assert_eq!(terms(&on_the_seam(b"T")), [Consequence::FrameshiftVariant]);
        let past_the_seam = InputVariant::new("1".into(), 159, 158, b"-".to_vec(), b"TTT".to_vec());
        assert_eq!(terms(&past_the_seam), [Consequence::InframeInsertion]);
    }

    /// An insertion standing on the bases an edited model holds and the genome
    /// lacks that brings the model's own bases with one run of three left out
    /// leaves the model three bases short, in frame: an `inframe_deletion`,
    /// whether the record is minimised to the inserted bases (CDS 9-14 of the
    /// model holding six) or keeps both genomic neighbours (CDS 8-15). Ensembl
    /// VEP asks `inframe_deletion` only of an allele shorter than its genomic
    /// span (`deletion`, Utils/Constants.pm) and withholds
    /// `protein_altering_variant` where that predicate holds
    /// (VariationEffect.pm:388), so its row for the spelling with both
    /// neighbours takes the default `intergenic_variant`
    /// (BaseVariationFeatureOverlapAllele.pm:285); on the minimised spelling
    /// it reads the reference codons with the span cut out
    /// (TranscriptVariationAllele.pm:872-873), where the predicate fails, and
    /// writes `protein_altering_variant`. Three bases that are not the model's
    /// own change a residue as well and stay `protein_altering_variant`, and
    /// an allele short of the model's bases by no multiple of three stays a
    /// `frameshift_variant`.
    #[test]
    fn insertion_three_bases_short_of_the_edited_model_bases_is_an_inframe_deletion() {
        let tx = coding_seam_model_holding("ACCGTC");
        let config = EffectsConfig::default();
        let minimised =
            |alt: &[u8]| InputVariant::new("1".into(), 158, 157, b"-".to_vec(), alt.to_vec());
        let with_both_neighbours =
            |alt: &[u8]| InputVariant::new("1".into(), 157, 158, b"GA".to_vec(), alt.to_vec());
        let lengths = |cell: Option<&str>| -> Vec<usize> {
            cell.map(|c| c.split('/').map(str::len).collect())
                .unwrap_or_default()
        };
        // Each spelling by position and kind: the variant, its CDS span, and
        // what it leaves out of the six held bases.
        let spellings = [
            (minimised(b"ACC"), "9-14", "minimised at 157|158, tail"),
            (minimised(b"GTC"), "9-14", "minimised at 157|158, head"),
            (
                with_both_neighbours(b"GACCA"),
                "8-15",
                "both neighbours kept at 157-158, tail",
            ),
            (
                with_both_neighbours(b"GGTCA"),
                "8-15",
                "both neighbours kept at 157-158, head",
            ),
        ];
        for (variant, cds_position, spelling) in &spellings {
            let allele_len = variant.alt_allele().len();
            let row = calculate_consequences(variant, &tx, &config).unwrap_or_else(|| {
                panic!("{spelling} left out, an allele of {allele_len} bases: no row is written")
            });
            assert_eq!(
                row.consequences.as_slice(),
                [Consequence::InframeDeletion],
                "{spelling} left out, an allele of {allele_len} bases"
            );
            assert_eq!(row.impact, Impact::MODERATE, "{spelling} left out");
            assert_eq!(
                row.cds_position.as_deref(),
                Some(*cds_position),
                "{spelling} left out"
            );
            assert_eq!(
                row.protein_position.as_deref(),
                Some("3-5"),
                "{spelling} left out"
            );
            assert_eq!(
                lengths(row.amino_acids.as_deref()),
                [3, 2],
                "{spelling} left out: residues before and after, counted"
            );
            assert_eq!(
                lengths(row.codons.as_deref()),
                [9, 6],
                "{spelling} left out: codon bases before and after, counted"
            );
        }

        let terms = |variant: &InputVariant| {
            calculate_consequences(variant, &tx, &config)
                .map(|row| row.consequences.to_vec())
                .unwrap_or_default()
        };
        assert_eq!(
            terms(&minimised(b"TTT")),
            [Consequence::ProteinAlteringVariant],
            "minimised at 157|158, three bases that are not the model's"
        );
        assert_eq!(
            terms(&minimised(b"ACCG")),
            [Consequence::FrameshiftVariant],
            "minimised at 157|158, four of the six held bases"
        );
        assert_eq!(
            terms(&minimised(b"ACCGTC")),
            [Consequence::SynonymousVariant],
            "minimised at 157|158, the six held bases themselves"
        );
    }

    /// The HGVSc of an insertion standing on the bases an edited model holds
    /// and the genome lacks names the replacement of those bases, the span the
    /// row's `CDS_position` gives: on the model whose CDS 9-11 is `ACC`
    /// between genomic 157 and 158, an allele sharing no end with them is a
    /// `delins` of the three, and one sharing an end is clipped to what
    /// changes, as any allele replacing reference bases is: a deletion of the
    /// bases it leaves out, a substitution of the one base it changes, an
    /// insertion after the bases it keeps. The bases repeated are a
    /// duplication, their reverse complement an inversion, and the bases
    /// themselves, which change nothing on the model, have no notation.
    /// `ins` between the two bases around them would read `c.8_12ins`,
    /// positions that are not neighbours; Ensembl VEP slides the two by the
    /// misalignment offset (TranscriptVariationAllele.pm:1524-1525). An
    /// allele the genome shows as a duplication of the base after it
    /// replaces the model's bases like any other, and an insertion one base
    /// on, off the edit, is an `ins` between neighbours.
    #[test]
    fn hgvsc_of_an_insertion_on_edited_model_bases_names_their_replacement() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, false);
        let tx = coding_seam_model();
        let hgvsc = |start: u64, alt: &[u8]| {
            let variant =
                InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec());
            calculate_consequences(&variant, &tx, &config)
                .expect("an insertion inside a coding exon has a row")
                .hgvsc
        };
        let on_the_seam = |alt: &[u8]| hgvsc(158, alt);
        assert_eq!(
            on_the_seam(b"G").as_deref(),
            Some("NM_000004.1:c.9_11delinsG")
        );
        assert_eq!(
            on_the_seam(b"GGG").as_deref(),
            Some("NM_000004.1:c.9_11delinsGGG")
        );
        assert_eq!(on_the_seam(b"A").as_deref(), Some("NM_000004.1:c.10_11del"));
        assert_eq!(on_the_seam(b"GCC").as_deref(), Some("NM_000004.1:c.9A>G"));
        assert_eq!(
            on_the_seam(b"ACCG").as_deref(),
            Some("NM_000004.1:c.11_12insG")
        );
        assert_eq!(
            on_the_seam(b"ACCACC").as_deref(),
            Some("NM_000004.1:c.9_11dup")
        );
        assert_eq!(
            on_the_seam(b"GGT").as_deref(),
            Some("NM_000004.1:c.9_11inv")
        );
        assert_eq!(on_the_seam(b"ACC"), None);
        assert_eq!(
            on_the_seam(b"T").as_deref(),
            Some("NM_000004.1:c.9_11delinsT")
        );
        assert_eq!(hgvsc(159, b"G").as_deref(), Some("NM_000004.1:c.12_13insG"));
    }

    /// The replacement is numbered like any exonic position: before the start
    /// codon on the model whose 5' UTR holds `AAA` at cDNA 10-12, 44 to 42
    /// bases before a coding start at 54; and on the reverse strand, where the
    /// non-coding model edited by `21 20 AAA` holds cDNA 21-23 between genomic
    /// 1281 and 1280 and the allele is spelled on the transcript's strand. The
    /// model's bases are read from the cDNA the cache holds, so a model
    /// without one has no notation there.
    #[test]
    fn hgvsc_of_an_insertion_on_edited_model_bases_outside_the_coding_sequence() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, false);
        let insertion = |start: u64, alt: &[u8]| {
            InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec())
        };

        let utr = insertion_edited_model("NM_000001.1");
        let tc = calculate_consequences(&insertion(109, b"G"), &utr, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000001.1:c.-44_-42delinsG"));
        assert_eq!(tc.cdna_position.as_deref(), Some("10-12"));

        let mut reverse = reverse_non_coding_seam_model();
        reverse.vefc.as_mut().unwrap().spliced_seq = None;
        reverse.derived = Default::default();
        let without_cdna = calculate_consequences(&insertion(1281, b"G"), &reverse, &config)
            .expect("an insertion inside an exon has a row");
        assert_eq!(without_cdna.hgvsc, None);

        let reverse = reverse_non_coding_seam_model();
        let tc = calculate_consequences(&insertion(1281, b"C"), &reverse, &config).unwrap();
        assert_eq!(tc.hgvsc.as_deref(), Some("NR_000001.1:n.21_23delinsG"));
        assert_eq!(tc.cdna_position.as_deref(), Some("21-23"));
    }

    /// The reverse-strand non-coding two-exon model as an `NR_` model edited by
    /// `21 20 AAA`, as the cache holds such a model: the mapper folds the edit
    /// (cDNA 1-20 at 1281-1300, 24-103 at 1201-1280, 104-203 at 1001-1100), so
    /// cDNA 21-23 lie between genomic 1281 and 1280, and the cached cDNA reads
    /// 20 `A`, then `C` to its end.
    fn reverse_non_coding_seam_model() -> Transcript {
        use vep_core::transcript::{Attribute, ExonCoordMapper, MapperPair};
        let mut tx = non_coding_reverse_two_exon_transcript();
        tx.stable_id = "NR_000001.1".into();
        tx.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: "21 20 AAA".to_string(),
        });
        let pair = |from_start, from_end, to_start, to_end| MapperPair {
            from_start,
            from_end,
            to_start,
            to_end,
            ori: -1,
        };
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.mapper.as_mut().unwrap().exon_coord_mapper = ExonCoordMapper::new(vec![
            pair(1, 20, 1281, 1300),
            pair(24, 103, 1201, 1280),
            pair(104, 203, 1001, 1100),
        ]);
        vefc.spliced_seq = Some(format!("{}{}", "A".repeat(20), "C".repeat(183)));
        tx.derived = Default::default();
        tx
    }

    /// What the clip leaves of a replacement of model bases is placed by the 3'
    /// rule on the model's own sequence only while the 3' shift is on: under
    /// `--shift_hgvs 0` every notation stays where it is, among bases the
    /// genome lacks as anywhere else. With the shift, the duplication of `ACC`
    /// before an `A` is that of `CCA` one base on, an `A` inserted before the
    /// four `A` of CDS 12-15 duplicates the last of them, a `C` inserted after
    /// the `AAA` of the 5' UTR duplicates the last `C` before the start codon,
    /// and a `C` deleted from a run of `C` that fills the rest of the reverse
    /// model is the exon's last base, 103, though the next exon goes on with
    /// `C`: the intron between the two begins `A`. Without it each stays where
    /// the clip leaves it: the duplication of CDS 9-11, an insertion after CDS
    /// 11, an insertion after the `AAA`, the deletion of cDNA 23. The genome,
    /// `T` throughout, moves none of these alleles under either setting.
    #[test]
    fn replaced_model_bases_are_placed_3prime_only_while_the_shift_is_on() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let unshifted = hgvs_config(&fasta, false);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        assert!(shifted.hgvs.shift_hgvs && !unshifted.hgvs.shift_hgvs);
        let hgvsc = |tx: &Transcript, start: u64, alt: &[u8], config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an insertion inside an exon has a row")
                .hgvsc
        };
        let coding = coding_seam_model();
        let utr = insertion_edited_model("NM_000001.1");
        let reverse = reverse_non_coding_seam_model();
        let cases: [(&Transcript, u64, &[u8], &str, &str); 4] = [
            (
                &coding,
                158,
                b"ACCACC",
                "NM_000004.1:c.10_12dup",
                "NM_000004.1:c.9_11dup",
            ),
            (
                &coding,
                158,
                b"ACCA",
                "NM_000004.1:c.15dup",
                "NM_000004.1:c.11_12insA",
            ),
            (
                &utr,
                109,
                b"AAAC",
                "NM_000001.1:c.-1dup",
                "NM_000001.1:c.-42_-41insC",
            ),
            (
                &reverse,
                1281,
                b"GG",
                "NR_000001.1:n.103del",
                "NR_000001.1:n.23del",
            ),
        ];
        for (tx, start, alt, with_the_shift, without_it) in cases {
            assert_eq!(
                hgvsc(tx, start, alt, &shifted).as_deref(),
                Some(with_the_shift)
            );
            assert_eq!(
                hgvsc(tx, start, alt, &unshifted).as_deref(),
                Some(without_it)
            );
        }
    }

    /// The compact coding transcript as a RefSeq model whose alignment
    /// inserted `copies` of `GCA` after unedited cDNA 58, in a coding sequence
    /// that reads `ATG CAG CA` before them and `A AAA TTC` after: the mapper
    /// folds the edit, so genomic 152-157, CDS 3-8, repeat `GCA` twice, the
    /// model's own bases go on with the repeat from CDS 9, and genomic 158, an
    /// `A`, ends it on the genome.
    fn repeat_seam_model(copies: usize) -> Transcript {
        use vep_core::transcript::{Attribute, ExonCoordMapper};
        let added = 3 * copies as u64;
        let mut tx = compact_transcript_with_cached_cdna();
        tx.stable_id = "NM_000005.1".into();
        tx.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: format!("59 58 {}", "GCA".repeat(copies)),
        });
        tx.cdna_coding_end = Some(900 + added);
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.translateable_seq
            .as_mut()
            .unwrap()
            .replace_range(..8, &format!("ATGCAGCA{}", "GCA".repeat(copies)));
        let mapper = vefc.mapper.as_mut().unwrap();
        mapper.cdna_coding_end = 900 + added;
        mapper.exon_coord_mapper = ExonCoordMapper::new(vec![
            pair(1, 58, 100, 157),
            pair(59 + added, 300 + added, 158, 399),
            pair(301 + added, 600 + added, 2100, 2399),
            pair(601 + added, 2601 + added, 4100, 6100),
        ]);
        tx.derived = Default::default();
        assert_eq!(tx.spliced_seq().map(|s| &s[50..58]), Some("ATGCAGCA"));
        tx
    }

    /// An insertion the genome shows as a duplication, standing on the bases
    /// an edited model holds and the genome lacks, is described on the model
    /// as any other insertion there is: the allele in place of the model's
    /// bases, the reading the row's HGVSp takes. `GCA` written after CDS 5 is
    /// carried by the 3' shift through the genome's second copy to the bases
    /// the model adds after CDS 8. Where those are one more `GCA` the allele
    /// is the model's own sequence and has no notation; two copies of the
    /// allele duplicate the model's CDS 9-11; and one copy standing on a
    /// model that adds two deletes CDS 12-14. Ensembl VEP types the notation
    /// on the genomic slice (TranscriptVariationAllele.pm:1474-1478), where
    /// the allele repeats the bases before the point
    /// (Utils/Sequence.pm:577-585), and numbers those bases on the model: a
    /// duplication of CDS 6-8 or 3-8, a sequence longer than the one its
    /// HGVSp describes by the model's bases. The term and the position
    /// columns read the allele where the record writes it. An allele written
    /// on the model's bases is read the same way under either shift setting,
    /// and one the shift leaves off them keeps the genome's duplication.
    #[test]
    fn a_duplication_the_genome_shows_on_edited_model_bases_is_read_on_the_model() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let genome = [
            bases_from(150, "ATGCAGCA"),
            bases_from(158, "AAAATTCGAT"),
            bases_from(1278, "GGG"),
        ]
        .concat();
        let (_dir, fasta) = hgvsc_fasta_fixture(&genome);
        let fasta = Arc::new(fasta);
        let unshifted = hgvs_config(&fasta, false);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        let row = |tx: &Transcript, start: u64, alt: &[u8], config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an insertion inside an exon has a row")
        };
        let one_copy = repeat_seam_model(1);
        let two_copies = repeat_seam_model(2);

        let moved = row(&one_copy, 155, b"GCA", &shifted);
        assert_eq!(moved.hgvs_offset, Some(3));
        assert_eq!(moved.hgvsc, None);
        assert!(moved.hgvsp.as_deref().is_some_and(|p| p.ends_with('=')));
        assert_eq!(
            moved.consequences.as_slice(),
            [Consequence::InframeInsertion]
        );
        assert_eq!(moved.cds_position.as_deref(), Some("5-6"));
        assert_eq!(
            row(&one_copy, 155, b"GCAGCA", &shifted).hgvsc.as_deref(),
            Some("NM_000005.1:c.9_11dup")
        );
        assert_eq!(
            row(&two_copies, 155, b"GCA", &shifted).hgvsc.as_deref(),
            Some("NM_000005.1:c.12_14del")
        );

        for config in [&shifted, &unshifted] {
            assert_eq!(row(&one_copy, 158, b"GCA", config).hgvsc, None);
        }
        assert_eq!(
            row(&one_copy, 155, b"GCA", &unshifted).hgvsc.as_deref(),
            Some("NM_000005.1:c.6_8dup")
        );

        // On the reverse strand the model's cDNA 21-23, `C` in the cached
        // cDNA, follow the genome's run of `A`, cDNA 1-20: the shift carries
        // an `A` inserted at genomic 1291 to them, where the genome shows the
        // duplication of cDNA 20.
        let reverse = reverse_non_coding_seam_model();
        assert_eq!(
            row(&reverse, 1291, b"T", &shifted).hgvsc.as_deref(),
            Some("NR_000001.1:n.21_23delinsA")
        );
    }

    /// An allele that keeps the base it shares with its reference, as a
    /// record with a second alternate leaves it, spans that base and not the
    /// point after it: on the model whose CDS 9-11, `ACC`, lie between
    /// genomic 157 and 158, `T` to `TGGG` at genomic 157 is read at CDS 8,
    /// and its three inserted bases stand between CDS 8 and CDS 9, the
    /// model's own base next to it, two neighbouring positions with the
    /// model's bases untouched. Inserting the model's own `ACC` there
    /// duplicates them: CDS 9-11 where the record writes them and, by the 3'
    /// rule, CDS 10-12 before the `A` of CDS 12. An allele that repeats its
    /// shared base stays the duplication of that base. Ensembl VEP numbers
    /// the two genomic bases around the point from the unedited exons
    /// (TranscriptVariationAllele.pm:1524-1525) and respells the allele with
    /// its shared base (`:1513`). On the reverse strand the shared base,
    /// cDNA 24, follows the model's cDNA 21-23, and the inserted bases stand
    /// before it; three more `C` there duplicate the model's own, or the
    /// last three of the exon's run by the 3' rule.
    #[test]
    fn an_allele_that_keeps_a_shared_base_is_inserted_beside_edited_model_bases() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let unshifted = hgvs_config(&fasta, false);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        let coding = coding_seam_model();
        let reverse = reverse_non_coding_seam_model();
        let keeping = |tx: &Transcript, at: u64, alt: &[u8], config: &EffectsConfig| {
            let variant = InputVariant::new("1".into(), at, at, b"T".to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an allele inside an exon has a row")
        };
        for config in [&unshifted, &shifted] {
            let row = keeping(&coding, 157, b"TGGG", config);
            assert_eq!(row.cds_position.as_deref(), Some("8"));
            assert_eq!(row.hgvsc.as_deref(), Some("NM_000004.1:c.8_9insGGG"));
            assert_eq!(
                keeping(&coding, 157, b"TT", config).hgvsc.as_deref(),
                Some("NM_000004.1:c.8dup")
            );
            assert_eq!(
                keeping(&reverse, 1280, b"TA", config).hgvsc.as_deref(),
                Some("NR_000001.1:n.23_24insT")
            );
        }
        let own_bases: [(&Transcript, u64, &[u8], &str, &str); 2] = [
            (
                &coding,
                157,
                b"TACC",
                "NM_000004.1:c.10_12dup",
                "NM_000004.1:c.9_11dup",
            ),
            (
                &reverse,
                1280,
                b"TGGG",
                "NR_000001.1:n.101_103dup",
                "NR_000001.1:n.21_23dup",
            ),
        ];
        for (tx, at, alt, with_the_shift, without_it) in own_bases {
            assert_eq!(
                keeping(tx, at, alt, &shifted).hgvsc.as_deref(),
                Some(with_the_shift)
            );
            assert_eq!(
                keeping(tx, at, alt, &unshifted).hgvsc.as_deref(),
                Some(without_it)
            );
        }
    }

    /// On a model with an `_rna_edit` the notation spells the alleles relative
    /// to the transcript, and it spells what the clip leaves of them, as on
    /// any other model: of `T` to `TGGG` the three inserted bases, of `TT` to
    /// `TGGG` the `GGG` that replaces the second base, of `TTT` to `TGT` the
    /// one base that changes. Ensembl VEP clips the alleles
    /// (TranscriptVariationAllele.pm:1492) and then writes both over from
    /// `feature_seq` (`:1510-1513`), so its notation on such a model carries
    /// the bases the allele shares with its reference: `ins` of four bases
    /// where three are inserted, a substitution that names three bases a
    /// side. An allele with no shared base is spelled the same either way.
    #[test]
    fn hgvsc_on_an_edited_model_spells_the_alleles_the_clip_leaves() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, false);
        let hgvsc = |tx: &Transcript, start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            let variant =
                InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, &config)
                .expect("an allele inside a coding exon has a row")
                .hgvsc
        };
        // Genomic 160 is CDS 11 on both models: the edit lies in the 5' UTR.
        let edited = insertion_edited_model("NM_000001.1");
        let on_edited = |start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            hgvsc(&edited, start, end, reference, alt)
        };
        assert_eq!(
            on_edited(160, 160, b"T", b"TGGG").as_deref(),
            Some("NM_000001.1:c.11_12insGGG")
        );
        assert_eq!(
            on_edited(160, 161, b"TT", b"TGGG").as_deref(),
            Some("NM_000001.1:c.12delinsGGG")
        );
        assert_eq!(
            on_edited(160, 162, b"TTT", b"TGT").as_deref(),
            Some("NM_000001.1:c.12T>G")
        );
        assert_eq!(
            on_edited(160, 160, b"T", b"G").as_deref(),
            Some("NM_000001.1:c.11T>G")
        );
        assert_eq!(
            on_edited(161, 160, b"-", b"GGG").as_deref(),
            Some("NM_000001.1:c.11_12insGGG")
        );
        let unedited = compact_transcript_with_cached_cdna();
        assert_eq!(
            hgvsc(&unedited, 160, 160, b"T", b"TGGG").as_deref(),
            Some("ENST00000000001.1:c.11_12insGGG")
        );
    }

    /// The respelled alleles lose the bases the notation's span lost and no
    /// other. Genomic 150-152 is `TTT` where the model's CDS 1-3 is `ATG`. Of
    /// `TT` to `GG` at CDS 2-3 the genome shares no base, so the span keeps
    /// both positions and the `G` the allele shares with the model alone
    /// stays in the notation: two bases replace two, one base replaces two
    /// for `TT` to `G`, and two replace the model's `G` for `T` to `GC`.
    /// `TT` to `TC` shares its first base with the genome and with the model,
    /// and the notation names the one base that changes.
    #[test]
    fn hgvsc_on_an_edited_model_drops_only_the_bases_the_span_lost() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, true);
        let edited = insertion_edited_model("NM_000001.1");
        let hgvsc = |start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            let variant =
                InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            calculate_consequences(&variant, &edited, &config)
                .expect("an allele inside a coding exon has a row")
                .hgvsc
        };
        assert_eq!(
            hgvsc(151, 152, b"TT", b"GG").as_deref(),
            Some("NM_000001.1:c.2_3delinsGG")
        );
        assert_eq!(
            hgvsc(151, 152, b"TT", b"G").as_deref(),
            Some("NM_000001.1:c.2_3delinsG")
        );
        assert_eq!(
            hgvsc(152, 152, b"T", b"GC").as_deref(),
            Some("NM_000001.1:c.3delinsGC")
        );
        assert_eq!(
            hgvsc(150, 152, b"TTT", b"ACG").as_deref(),
            Some("NM_000001.1:c.1_3delinsACG")
        );
        assert_eq!(
            hgvsc(151, 152, b"TT", b"TC").as_deref(),
            Some("NM_000001.1:c.3G>C")
        );
    }

    /// A base is clipped from a substitution's respelled pair only where the
    /// model itself has it. Genomic 150-158 is all `T` where the model's CDS
    /// 1-4 is `ATGG` and its CDS 8-9 `GA`. Of `TT` to `TC` over CDS 1-2 the
    /// clip takes the
    /// first base on the genome's alleles, which the model spells `A`: the
    /// record's two bases stand in place of the model's two, the span the
    /// row's own position cells cover. The same holds at the other end (`TT`
    /// to `CT` over CDS 2-3). A three-base allele keeps every base from the
    /// one the model spells differently (`TTT` to `TTC` over CDS 1-3) and
    /// loses an end base the model has (`TTT` to `TCT` over CDS 2-4, whose
    /// first base is `T` on both). Two bases that are the reverse complement
    /// of the model's own are an inversion, as any allele standing in place
    /// of reference bases is typed (`TT` to `TC` over CDS 8-9, the model's
    /// `GA`). Ensembl VEP writes each as a substitution at the one position
    /// the clip leaves, the model's whole allele against the record's
    /// (TranscriptVariationAllele.pm:1492, 1510-1513; Utils/Sequence.pm:656).
    #[test]
    fn a_substitution_keeps_a_clipped_base_the_edited_model_spells_differently() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, true);
        let edited = insertion_edited_model("NM_000001.1");
        let row = |start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            let variant =
                InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            calculate_consequences(&variant, &edited, &config)
                .expect("an allele inside a coding exon has a row")
        };
        let first_kept = row(150, 151, b"TT", b"TC");
        assert_eq!(first_kept.used_ref.as_deref(), Some("AT"));
        assert_eq!(first_kept.cds_position.as_deref(), Some("1-2"));
        assert_eq!(
            first_kept.hgvsc.as_deref(),
            Some("NM_000001.1:c.1_2delinsTC")
        );
        assert_eq!(
            row(151, 152, b"TT", b"CT").hgvsc.as_deref(),
            Some("NM_000001.1:c.2_3delinsCT")
        );
        assert_eq!(
            row(150, 152, b"TTT", b"TTC").hgvsc.as_deref(),
            Some("NM_000001.1:c.1_3delinsTTC")
        );
        assert_eq!(
            row(151, 153, b"TTT", b"TCT").hgvsc.as_deref(),
            Some("NM_000001.1:c.3_4delinsCT")
        );
        assert_eq!(
            row(157, 158, b"TT", b"TC").hgvsc.as_deref(),
            Some("NM_000001.1:c.8_9inv")
        );
    }

    /// A substitution beside the bases an edited model holds and the genome
    /// lacks spans them. The genome carries the model's own bases around the
    /// seam, `TGG` at genomic 155-157 (CDS 6-8) and `AAA` at 158-160 (CDS
    /// 12-14), and the model's reference allele over genomic 157-158 is CDS
    /// 8-12, `GACCA`. Of `GA` to `GC` the clip takes the first base, which
    /// the model has, and the record's one base stands in place of the model's
    /// four from CDS 9, the first base the genome lacks, to CDS 12, the base
    /// the record changes. A three-base allele that ends at the seam reads
    /// the same. The reverse model's cDNA 21-23 lie between genomic 1281 and
    /// 1280: `GT` to `CT` keeps genomic 1281, cDNA 20, and replaces cDNA
    /// 21-24, and `GT` to `GA` keeps genomic 1280, cDNA 24, and replaces cDNA
    /// 20-23. An allele that runs across the seam stops the clip there,
    /// though the model's first own base is the record's next one (`GGT` to
    /// `CGT` over cDNA 20-25); read against the record's own reference the
    /// pair holds none of the model's bases, the clip runs across the seam
    /// and the notation names the one base the record changes. A base beside
    /// the seam that the model spells differently keeps the whole span.
    /// Ensembl VEP writes each as a substitution at the one position the clip
    /// leaves, the model's whole reference allele against the record's
    /// (TranscriptVariation.pm:463; TranscriptVariationAllele.pm:1492,
    /// 1510-1513).
    #[test]
    fn a_substitution_beside_edited_model_bases_spans_them() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let mut genome = bases_from(155, "TGGAAA");
        genome.extend(bases_from(1279, "GG"));
        let (_dir, fasta) = hgvsc_fasta_fixture(&genome);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, true);
        let coding = coding_seam_model();
        let reverse = reverse_non_coding_seam_model();
        let row = |tx: &Transcript, start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            let variant =
                InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, &config)
                .expect("an allele inside an exon has a row")
        };
        let after_the_seam = row(&coding, 157, 158, b"GA", b"GC");
        assert_eq!(after_the_seam.used_ref.as_deref(), Some("GACCA"));
        assert_eq!(after_the_seam.cds_position.as_deref(), Some("8-12"));
        assert_eq!(
            after_the_seam.hgvsc.as_deref(),
            Some("NM_000004.1:c.9_12delinsC")
        );
        assert_eq!(
            row(&coding, 156, 158, b"GGA", b"GGC").hgvsc.as_deref(),
            Some("NM_000004.1:c.9_12delinsC")
        );
        let on_the_reverse_strand = row(&reverse, 1280, 1281, b"GT", b"CT");
        assert_eq!(on_the_reverse_strand.used_ref.as_deref(), Some("GGGGT"));
        assert_eq!(
            on_the_reverse_strand.hgvsc.as_deref(),
            Some("NR_000001.1:n.21_24delinsG")
        );
        assert_eq!(
            row(&reverse, 1280, 1281, b"GT", b"GA").hgvsc.as_deref(),
            Some("NR_000001.1:n.20_23delinsT")
        );
        assert_eq!(
            row(&reverse, 1279, 1281, b"GGT", b"CGT").hgvsc.as_deref(),
            Some("NR_000001.1:n.21_25delinsCG")
        );
        let against_the_record = hgvs_config(&fasta, false);
        let variant = InputVariant::new("1".into(), 1279, 1281, b"GGT".to_vec(), b"CGT".to_vec());
        assert_eq!(
            calculate_consequences(&variant, &reverse, &against_the_record)
                .expect("an allele inside an exon has a row")
                .hgvsc
                .as_deref(),
            Some("NR_000001.1:n.25C>G")
        );
        let (_dir, all_t) = hgvsc_fasta_fixture(&[]);
        let config = hgvs_config(&Arc::new(all_t), true);
        let variant = InputVariant::new("1".into(), 157, 158, b"TT".to_vec(), b"TC".to_vec());
        assert_eq!(
            calculate_consequences(&variant, &coding, &config)
                .expect("an allele inside a coding exon has a row")
                .hgvsc
                .as_deref(),
            Some("NM_000004.1:c.8_12delinsTC")
        );
    }

    /// A `delins` the clip leaves beside the bases an edited model holds and
    /// the genome lacks spans them, as the substitution there does. On the
    /// coding model, whose reference allele over genomic 157-158 is CDS 8-12,
    /// `GA` to `GCT` keeps its first base, which the model has, and the two
    /// bases the clip leaves stand in place of the model's four from CDS 9,
    /// the first base the genome lacks, to CDS 12, the last base the record
    /// changes; a longer allele and a three-base reference read the same
    /// way. The reverse model's cDNA 21-23 lie between genomic 1281 and
    /// 1280: `GT` to `CAT` keeps genomic 1281, cDNA 20, and replaces cDNA
    /// 21-24, and `GT` to `GAC` keeps genomic 1280, cDNA 24, and replaces
    /// cDNA 20-23. An allele that keeps two bases across the seam stops the
    /// clip there (`GGT` to `CCGT` over cDNA 20-25). A `del` the clip
    /// leaves beside the seam takes the model's bases with it: `GGA` to
    /// `GG` over genomic 156-158 deletes CDS 9-12, and `GGA` to `GA`, which
    /// keeps a base on each side, deletes CDS 8-11. Read against the
    /// record's own reference the pair holds none of the model's bases, and
    /// the notation names the span the clip leaves. Ensembl VEP writes each
    /// over the span the clip leaves, a `delins` as the record's whole
    /// alternate allele (TranscriptVariation.pm:463;
    /// TranscriptVariationAllele.pm:1492, 1513; Utils/Sequence.pm:659-667).
    #[test]
    fn a_delins_or_a_del_beside_edited_model_bases_spans_them() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let mut genome = bases_from(155, "TGGAAA");
        genome.extend(bases_from(1279, "GG"));
        let (_dir, fasta) = hgvsc_fasta_fixture(&genome);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, true);
        let coding = coding_seam_model();
        let reverse = reverse_non_coding_seam_model();
        let row_under = |tx: &Transcript,
                         start: u64,
                         end: u64,
                         reference: &[u8],
                         alt: &[u8],
                         config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an allele inside an exon has a row")
        };
        let row = |tx: &Transcript, start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            row_under(tx, start, end, reference, alt, &config)
        };
        let after_the_seam = row(&coding, 157, 158, b"GA", b"GCT");
        assert_eq!(after_the_seam.used_ref.as_deref(), Some("GACCA"));
        assert_eq!(after_the_seam.cds_position.as_deref(), Some("8-12"));
        assert_eq!(
            after_the_seam.hgvsc.as_deref(),
            Some("NM_000004.1:c.9_12delinsCT")
        );
        assert_eq!(
            row(&coding, 157, 158, b"GA", b"GCTT").hgvsc.as_deref(),
            Some("NM_000004.1:c.9_12delinsCTT")
        );
        assert_eq!(
            row(&coding, 157, 159, b"GAA", b"GC").hgvsc.as_deref(),
            Some("NM_000004.1:c.9_13delinsC")
        );
        assert_eq!(
            row(&coding, 157, 159, b"GAA", b"GCT").hgvsc.as_deref(),
            Some("NM_000004.1:c.9_13delinsCT")
        );
        let on_the_reverse_strand = row(&reverse, 1280, 1281, b"GT", b"CAT");
        assert_eq!(on_the_reverse_strand.used_ref.as_deref(), Some("GGGGT"));
        assert_eq!(
            on_the_reverse_strand.cdna_position.as_deref(),
            Some("20-24")
        );
        assert_eq!(
            on_the_reverse_strand.hgvsc.as_deref(),
            Some("NR_000001.1:n.21_24delinsTG")
        );
        assert_eq!(
            row(&reverse, 1280, 1281, b"GT", b"GAC").hgvsc.as_deref(),
            Some("NR_000001.1:n.20_23delinsGT")
        );
        assert_eq!(
            row(&reverse, 1279, 1281, b"GGT", b"CCGT").hgvsc.as_deref(),
            Some("NR_000001.1:n.21_25delinsCGG")
        );
        assert_eq!(
            row(&coding, 156, 158, b"GGA", b"GG").hgvsc.as_deref(),
            Some("NM_000004.1:c.9_12del")
        );
        assert_eq!(
            row(&coding, 156, 158, b"GGA", b"GA").hgvsc.as_deref(),
            Some("NM_000004.1:c.8_11del")
        );
        let against_the_record = hgvs_config(&fasta, false);
        assert_eq!(
            row_under(&reverse, 1280, 1281, b"GT", b"CAT", &against_the_record)
                .hgvsc
                .as_deref(),
            Some("NR_000001.1:n.24delinsTG")
        );
    }

    /// The model's bases are replaced only where the row's span holds them:
    /// the minimised insertion, whose span is the model's CDS 9-11, and an
    /// allele that keeps a base on either side of the point, `TT` to `TGT`
    /// over genomic 157-158, whose span runs from CDS 8 to CDS 12. Both name
    /// the replacement of CDS 9-11 by the inserted bases.
    #[test]
    fn inserted_bases_replace_edited_model_bases_only_inside_the_row_span() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let config = hgvs_config(&fasta, false);
        let tx = coding_seam_model();
        let row = |start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            let variant =
                InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            calculate_consequences(&variant, &tx, &config)
                .expect("an allele inside a coding exon has a row")
        };
        let minimised = row(158, 157, b"-", b"GGG");
        assert_eq!(minimised.cds_position.as_deref(), Some("9-11"));
        assert_eq!(
            minimised.hgvsc.as_deref(),
            Some("NM_000004.1:c.9_11delinsGGG")
        );
        assert_eq!(row(158, 157, b"-", b"ACC").hgvsc, None);
        let across = row(157, 158, b"TT", b"TGT");
        assert_eq!(across.cds_position.as_deref(), Some("8-12"));
        assert_eq!(across.hgvsc.as_deref(), Some("NM_000004.1:c.9_11delinsG"));
    }

    /// Under `--use_transcript_ref` the row is read against the model's own
    /// bases, `G` at CDS 8 and `A` at CDS 12 where the genome has `T`, so the
    /// alleles respelled relative to the transcript share no base and stay
    /// whole. The bases standing beside the model's CDS 9-11 are the ones the
    /// clip leaves of the allele all the same: `T` to `TGGG` at genomic 157
    /// inserts `GGG` after CDS 8, and `TT` to `TGT` over genomic 157-158
    /// replaces CDS 9-11 by `G`, the notations of the same rows read against
    /// the record's own reference allele.
    #[test]
    fn bases_beside_edited_model_bases_are_the_ones_the_clip_leaves_under_either_reference() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let fasta = Arc::new(fasta);
        let tx = coding_seam_model();
        for use_transcript_ref in [false, true] {
            let config = hgvs_config(&fasta, use_transcript_ref);
            let row = |start: u64, end: u64, reference: &[u8], alt: &[u8]| {
                let variant =
                    InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
                calculate_consequences(&variant, &tx, &config)
                    .expect("an allele inside a coding exon has a row")
            };
            let keeping = row(157, 157, b"T", b"TGGG");
            assert_eq!(
                keeping.used_ref.as_deref(),
                use_transcript_ref.then_some("G"),
                "{use_transcript_ref}"
            );
            assert_eq!(
                keeping.hgvsc.as_deref(),
                Some("NM_000004.1:c.8_9insGGG"),
                "{use_transcript_ref}"
            );
            let across = row(157, 158, b"TT", b"TGT");
            assert_eq!(
                across.used_ref.as_deref(),
                use_transcript_ref.then_some("GACCA"),
                "{use_transcript_ref}"
            );
            assert_eq!(
                across.hgvsc.as_deref(),
                Some("NM_000004.1:c.9_11delinsG"),
                "{use_transcript_ref}"
            );
        }
    }

    /// The compact coding transcript as a RefSeq model that came with an
    /// attribute replacing unedited cDNA 59 by two of itself (`_rna_edit`
    /// `59 59 AA`), the shape of a coding sequence that joins over one
    /// position read twice: the mapper folds the edit (cDNA 1-58 and 60-301
    /// on exon 1), so genomic 157 is CDS 8, genomic 158 is CDS 10, and CDS
    /// 9, the base the genome lacks, lies between them. The coding sequence
    /// reads `ATG GCT GGA AAA ATT` where the genome's exons read
    /// `ATG GCT GGA AAA TTC`. `status` is the outcome of the cache's
    /// alignment on the model.
    fn twice_read_base_model(status: vep_core::transcript::BamEditStatus) -> Transcript {
        use vep_core::transcript::{Attribute, ExonCoordMapper};
        let mut tx = compact_transcript_with_cached_cdna();
        tx.stable_id = "NM_000006.1".into();
        tx.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: "59 59 AA".to_string(),
        });
        tx.bam_edit_status = Some(status);
        tx.cdna_coding_end = Some(901);
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.translateable_seq.as_mut().unwrap().insert(8, 'A');
        let mapper = vefc.mapper.as_mut().unwrap();
        mapper.cdna_coding_end = 901;
        mapper.exon_coord_mapper = ExonCoordMapper::new(vec![
            pair(1, 58, 100, 157),
            pair(60, 301, 158, 399),
            pair(302, 601, 2100, 2399),
            pair(602, 2602, 4100, 6100),
        ]);
        tx.derived = Default::default();
        assert_eq!(
            tx.spliced_seq().map(|s| &s[50..65]),
            Some("ATGGCTGGAAAAATT")
        );
        tx
    }

    /// On a model whose correction against the cache's alignment failed the
    /// genome is the reference: the model's sequence with the attribute it
    /// came with is not the aligned one (AnnotationType/Transcript.pm:579,
    /// 619-625), so the base its mapper holds between genomic 157 and 158 is
    /// no base of the transcript, and an insertion written there adds a base
    /// to the coding sequence. The row stands at the point after CDS 8 and
    /// is a `frameshift_variant`, the term, the HGVSc and the HGVSp of the
    /// same sequence written one base on, and an inserted base that is not
    /// the run's own is an `ins` between CDS 8 and CDS 9. Ensembl VEP gives
    /// the insertion the model's base as its span (`map_insert`,
    /// Mapper.pm:522-542), reads the allele in place of that base and so no
    /// length change (VariationEffect.pm:1556-1564): `synonymous_variant`
    /// beside the frameshift its own HGVSp states. A model the alignment
    /// corrected holds the base, and its rows are read on it; no other row
    /// of the failed model differs from that model's.
    #[test]
    fn an_insertion_on_a_base_a_failed_edit_model_holds_is_read_on_the_genome() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        use vep_core::transcript::BamEditStatus;
        let genome = [bases_from(150, "ATGGCTGG"), bases_from(158, "AAAATTCGAT")].concat();
        let (_dir, fasta) = hgvsc_fasta_fixture(&genome);
        let fasta = Arc::new(fasta);
        let unshifted = hgvs_config(&fasta, false);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        let failed = twice_read_base_model(BamEditStatus::Failed);
        let corrected = twice_read_base_model(BamEditStatus::Ok);
        let insertion = |tx: &Transcript, start: u64, alt: &[u8], config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an insertion inside a coding exon has a row")
        };

        // The run's own base, written on the point and one base on.
        let on_the_point = insertion(&failed, 158, b"A", &shifted);
        let one_base_on = insertion(&failed, 159, b"A", &shifted);
        assert_eq!(
            on_the_point.consequences.as_slice(),
            [Consequence::FrameshiftVariant]
        );
        assert_eq!(on_the_point.impact, Impact::HIGH);
        assert_eq!(on_the_point.cdna_position.as_deref(), Some("58-59"));
        assert_eq!(on_the_point.cds_position.as_deref(), Some("8-9"));
        assert_eq!(on_the_point.protein_position.as_deref(), Some("3"));
        assert_eq!(on_the_point.amino_acids.as_deref(), Some("G/GX"));
        assert_eq!(on_the_point.codons.as_deref(), Some("gga/ggAa"));
        assert_eq!(on_the_point.hgvs_offset, Some(4));
        assert_eq!(on_the_point.hgvsc.as_deref(), Some("NM_000006.1:c.13dup"));
        assert!(on_the_point
            .hgvsp
            .as_deref()
            .is_some_and(|p| p.contains("fs")));
        assert_eq!(one_base_on.consequences, on_the_point.consequences);
        assert_eq!(one_base_on.cds_position.as_deref(), Some("10-11"));
        assert_eq!(one_base_on.hgvsc, on_the_point.hgvsc);
        assert_eq!(one_base_on.hgvsp, on_the_point.hgvsp);
        assert_eq!(
            insertion(&failed, 158, b"A", &unshifted).hgvsc.as_deref(),
            Some("NM_000006.1:c.10dup")
        );

        // Another base, and the base before the point, which the genome
        // shows as a duplication.
        for config in [&shifted, &unshifted] {
            let other = insertion(&failed, 158, b"C", config);
            assert_eq!(
                other.consequences.as_slice(),
                [Consequence::FrameshiftVariant]
            );
            assert_eq!(other.cds_position.as_deref(), Some("8-9"));
            assert_eq!(other.hgvsc.as_deref(), Some("NM_000006.1:c.8_9insC"));
            assert!(other.hgvsp.as_deref().is_some_and(|p| p.contains("fs")));
            let repeated = insertion(&failed, 158, b"G", config);
            assert_eq!(
                repeated.consequences.as_slice(),
                [Consequence::FrameshiftVariant]
            );
            assert_eq!(repeated.hgvsc.as_deref(), Some("NM_000006.1:c.8dup"));
        }

        // The corrected model holds the base: the allele stands in its place.
        let own = insertion(&corrected, 158, b"A", &shifted);
        assert_eq!(
            own.consequences.as_slice(),
            [Consequence::SynonymousVariant]
        );
        assert_eq!(own.cds_position.as_deref(), Some("9"));
        let other = insertion(&corrected, 158, b"C", &shifted);
        assert_eq!(
            other.consequences.as_slice(),
            [Consequence::SynonymousVariant]
        );
        assert_eq!(other.cds_position.as_deref(), Some("9"));
        assert_eq!(other.hgvsc.as_deref(), Some("NM_000006.1:c.9A>C"));
        assert!(other.hgvsp.as_deref().is_some_and(|p| p.ends_with('=')));

        // Off the point the two models write the same rows.
        let off_the_point: [(u64, u64, &[u8], &[u8]); 5] = [
            (159, 158, b"-", b"A"),
            (157, 156, b"-", b"C"),
            (157, 157, b"G", b"C"),
            (158, 158, b"A", b"-"),
            (160, 160, b"A", b"C"),
        ];
        for (start, end, reference, alt) in off_the_point {
            let variant =
                InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            let rows = [&failed, &corrected]
                .map(|tx| calculate_consequences(&variant, tx, &shifted).unwrap());
            let cells = |row: &TranscriptConsequence| {
                (
                    row.consequences.clone(),
                    row.cdna_position.clone(),
                    row.cds_position.clone(),
                    row.protein_position.clone(),
                    row.amino_acids.clone(),
                    row.codons.clone(),
                    row.hgvsc.clone(),
                    row.hgvsp.clone(),
                )
            };
            assert_eq!(cells(&rows[0]), cells(&rows[1]), "{start}-{end}");
        }
    }

    /// On a model the genome is the reference of, an allele that keeps a
    /// base it shares with its reference is one insertion at the point after
    /// CDS 8, whichever side the kept base stands on, and its HGVSc is the
    /// one the same bases take written as a minimised insertion there: the
    /// 3' rule and the duplication test read the genome's `G AAAA TT` around
    /// the point, not the model's `G AAAAA TT`. Inserted `AC` moves over one
    /// genomic `A`, CDS 10, where the model's own base would stop it at CDS
    /// 9; `AG` then copies CDS 8 and CDS 10; a kept base after the point
    /// leaves `C` at the point, not after the model's base. A model the
    /// alignment corrected holds the base, and its notations are read on it.
    #[test]
    fn an_allele_that_keeps_a_shared_base_on_a_failed_edit_model_is_placed_on_the_genome() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        use vep_core::transcript::BamEditStatus;
        let genome = [bases_from(150, "ATGGCTGG"), bases_from(158, "AAAATTCGAT")].concat();
        let (_dir, fasta) = hgvsc_fasta_fixture(&genome);
        let fasta = Arc::new(fasta);
        let unshifted = hgvs_config(&fasta, false);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        let failed = twice_read_base_model(BamEditStatus::Failed);
        let corrected = twice_read_base_model(BamEditStatus::Ok);
        let hgvsc = |tx: &Transcript,
                     start: u64,
                     end: u64,
                     reference: &[u8],
                     alt: &[u8],
                     config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an allele inside a coding exon has a row")
                .hgvsc
        };
        // The inserted bases, with their notation under the 3' shift and
        // without it.
        let inserted: [(&[u8], &str, &str); 5] = [
            (b"A", "c.13dup", "c.10dup"),
            (b"C", "c.8_9insC", "c.8_9insC"),
            (b"G", "c.8dup", "c.8dup"),
            (b"AC", "c.10_11insCA", "c.8_9insAC"),
            (b"AG", "c.8_10dup", "c.8_9insAG"),
        ];
        for (bases, with_the_shift, without) in inserted {
            let kept_before = [b"G".as_slice(), bases].concat();
            let kept_after = [bases, b"A".as_slice()].concat();
            for (config, notation) in [(&shifted, with_the_shift), (&unshifted, without)] {
                let expected = Some(format!("NM_000006.1:{notation}"));
                let label = String::from_utf8_lossy(bases);
                assert_eq!(
                    hgvsc(&failed, 158, 157, b"-", bases, config),
                    expected,
                    "minimised {label}"
                );
                assert_eq!(
                    hgvsc(&failed, 157, 157, b"G", &kept_before, config),
                    expected,
                    "kept before {label}"
                );
                // Bases that start with the base after the point share it
                // at their front, and the clip leaves them after it.
                if bases[0] != b'A' {
                    assert_eq!(
                        hgvsc(&failed, 158, 158, b"A", &kept_after, config),
                        expected,
                        "kept after {label}"
                    );
                }
            }
        }

        // The corrected model holds the base: the 3' rule stops on it.
        assert_eq!(
            hgvsc(&corrected, 157, 157, b"G", b"GAC", &shifted).as_deref(),
            Some("NM_000006.1:c.9_10insCA")
        );
        assert_eq!(
            hgvsc(&corrected, 158, 158, b"A", b"CA", &shifted).as_deref(),
            Some("NM_000006.1:c.9_10insC")
        );
    }

    /// An allele that keeps the bases on both sides of the base a failed-edit
    /// model holds adds its other bases between two neighbours of the
    /// transcript: `GA` to `GCA` over genomic 157-158, CDS 8 and CDS 10 with
    /// the model's own CDS 9 between them, is the insertion of `C` written at
    /// that point, and takes its term, codons, residues, protein position,
    /// HGVSc, HGVSp and shift offset, as one, two and three other bases and a
    /// longer record that keeps two bases before the point do; a second
    /// alternate of another class does not hold it still. The cDNA and CDS position
    /// cells keep the record's span of three positions, the cells the row of
    /// the same record carries on a model the alignment corrected, and
    /// `USED_REF` keeps the model's three bases. Ensembl VEP reads the allele
    /// in place of all three positions, the model's base with them
    /// (BaseTranscriptVariation.pm:263-264;
    /// TranscriptVariationAllele.pm:2409-2410, 2421): a single base between
    /// the two is a change of no length, `synonymous_variant` here. The
    /// corrected model holds the base, and its row keeps that reading.
    #[test]
    fn bases_added_between_kept_bases_on_a_failed_edit_model_are_an_insertion() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        use vep_core::transcript::BamEditStatus;
        let failed = twice_read_base_model(BamEditStatus::Failed);
        let corrected = twice_read_base_model(BamEditStatus::Ok);
        // The genome carries the model's bases without its own CDS 9.
        let cdna = failed
            .spliced_seq()
            .expect("the model caches its cDNA")
            .to_string();
        let genome = [
            bases_from(150, &cdna[50..58]),
            bases_from(158, &cdna[59..69]),
        ]
        .concat();
        let (_dir, fasta) = hgvsc_fasta_fixture(&genome);
        let fasta = Arc::new(fasta);
        let unshifted = hgvs_config(&fasta, true);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, true)
        };
        let row = |tx: &Transcript,
                   start: u64,
                   end: u64,
                   reference: &[u8],
                   alt: &[u8],
                   config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an allele inside a coding exon has a row")
        };
        let read = |row: &TranscriptConsequence| {
            (
                row.consequences.clone(),
                row.impact,
                row.protein_position.clone(),
                row.amino_acids.clone(),
                row.codons.clone(),
                row.hgvsc.clone(),
                row.hgvsp.clone(),
                row.hgvs_offset,
            )
        };
        let span = |row: &TranscriptConsequence| {
            (
                row.cdna_position.clone(),
                row.cds_position.clone(),
                row.used_ref.clone(),
            )
        };

        let added: [(&[u8], &[u8]); 5] = [
            (b"C", b"GCA"),
            (b"A", b"GAA"),
            (b"G", b"GGA"),
            (b"CC", b"GCCA"),
            (b"CCC", b"GCCCA"),
        ];
        for config in [&shifted, &unshifted] {
            for (bases, keeping_both) in added {
                let kept = row(&failed, 157, 158, b"GA", keeping_both, config);
                let insertion = row(&failed, 158, 157, b"-", bases, config);
                assert_eq!(read(&kept), read(&insertion), "{}", bases.len());
                assert_eq!(kept.cdna_position.as_deref(), Some("58-60"));
                assert_eq!(kept.cds_position.as_deref(), Some("8-10"));
                assert_eq!(kept.used_ref.as_deref(), Some("GAA"));
                let on_the_corrected_model = row(&corrected, 157, 158, b"GA", keeping_both, config);
                assert_eq!(span(&kept), span(&on_the_corrected_model));
            }
        }
        let one_other_base = row(&failed, 157, 158, b"GA", b"GCA", &shifted);
        assert_eq!(
            one_other_base.consequences.as_slice(),
            [Consequence::FrameshiftVariant]
        );
        assert_eq!(
            one_other_base.hgvsc.as_deref(),
            Some("NM_000006.1:c.8_9insC")
        );
        // The run's own base moves 3' as the insertion does, though the record
        // carries a second alternate of another class.
        let mut beside_another_allele =
            InputVariant::new("1".into(), 157, 158, b"GA".to_vec(), b"GAA".to_vec());
        beside_another_allele.record_allele_string_multi = Some(Arc::from("GA/GAA/CC"));
        let moved = calculate_consequences(&beside_another_allele, &failed, &shifted)
            .expect("an allele inside a coding exon has a row");
        assert_eq!(moved.hgvsc.as_deref(), Some("NM_000006.1:c.13dup"));
        assert_eq!(moved.hgvs_offset, Some(4));
        assert_eq!(
            row(&failed, 157, 158, b"GA", b"GCCCA", &shifted)
                .consequences
                .as_slice(),
            [Consequence::InframeInsertion]
        );
        let two_kept_before = row(&failed, 156, 158, b"GGA", b"GGCA", &shifted);
        assert_eq!(read(&two_kept_before), read(&one_other_base));
        assert_eq!(two_kept_before.cds_position.as_deref(), Some("7-10"));

        let on_the_corrected_model = row(&corrected, 157, 158, b"GA", b"GCA", &shifted);
        assert_eq!(
            on_the_corrected_model.consequences.as_slice(),
            [Consequence::SynonymousVariant]
        );
        assert_eq!(on_the_corrected_model.cds_position.as_deref(), Some("8-10"));
    }

    /// The reference that carries the compact transcript's first 138 coding
    /// bases from genomic 150, `ATGGCTGGAAAATTCGAT` and then `GCT` repeated,
    /// with no base for the ones an edited model holds besides: under
    /// [`coding_seam_model`] CDS 1-8 lie at genomic 150-157 and CDS 12-21,
    /// `AAAATTCGAT`, from genomic 158.
    fn coding_seam_genome() -> (tempfile::TempDir, vep_fasta::IndexedFasta) {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let coding = format!("ATGGCTGGAAAATTCGAT{}", "GCT".repeat(40));
        hgvsc_fasta_fixture(&bases_from(150, &coding))
    }

    /// An insertion written on the bases an edited model holds and the genome
    /// lacks is described there even where the genome's 3' shift moves it. The
    /// genome has no such bases, so an allele that starts with the base after
    /// them shifts past them, and the notation Ensembl VEP types at the
    /// shifted point (TranscriptVariationAllele.pm:1446, 1474-1478) adds the
    /// allele after the model's bases and keeps them all, where the row's own
    /// `CDS_position` is those bases. On the model whose CDS 9-11 is `ACC`
    /// between genomic 157 and 158, before the `AAAA` of CDS 12-15:
    ///
    /// - `ATT` shifts one base and reads `c.12_13insTTA` at the shifted point,
    ///   an inserted residue. Where it is written it stands in place of `ACC`:
    ///   the shared `A` is clipped, `TT` replaces `CC`, and `CCA` becomes
    ///   `TTA`, `Pro4Leu`, the row's `missense_variant`.
    /// - `AA` shifts four bases and duplicates CDS 14-15 on the genome. In
    ///   place of `ACC` it leaves the `A` and replaces `CC` by `A`; the coding
    ///   sequence is one base shorter and codon 4 reads `AAA`.
    /// - `A` shifts four bases and duplicates CDS 15 on the genome. In place of
    ///   `ACC` it deletes `CC`, and codon 4 reads `AAA` again.
    ///
    /// Neither notation of such a row moved, so the row carries no
    /// `HGVS_OFFSET`. The bases `TTA` written one place on, off the model's
    /// bases, are an `ins` between neighbours, and an insertion the shift
    /// moves there keeps its shifted notation and its offset. Under
    /// `--shift_3prime 1` the row itself is annotated at the shifted point and
    /// keeps the notation of that point; so does the row of a model whose edit
    /// failed, which is read against the genome.
    #[test]
    fn insertion_the_shift_moves_off_edited_model_bases_is_described_where_it_is_written() {
        use vep_core::transcript::BamEditStatus;
        let (_dir, fasta) = coding_seam_genome();
        let fasta = Arc::new(fasta);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        let row = |tx: &Transcript, start: u64, alt: &[u8], config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an insertion inside a coding exon has a row")
        };
        let tx = coding_seam_model();
        let on_the_seam: [(&[u8], &str, &str); 3] = [
            (b"ATT", "c.10_11delinsTT", "p.Pro4Leu"),
            (b"AA", "c.10_11delinsA", "p.Pro4LysfsTer?"),
            (b"A", "c.10_11del", "p.Pro4LysfsTer?"),
        ];
        for (alt, hgvsc, hgvsp) in on_the_seam {
            let tc = row(&tx, 158, alt, &shifted);
            assert_eq!(tc.cds_position.as_deref(), Some("9-11"));
            assert_eq!(
                tc.hgvsc.as_deref(),
                Some(format!("NM_000004.1:{hgvsc}").as_str())
            );
            assert_eq!(
                tc.hgvsp.as_deref(),
                Some(format!("ENSP00000000001.1:{hgvsp}").as_str())
            );
            assert_eq!(tc.hgvs_offset, None);
        }
        assert_eq!(
            row(&tx, 158, b"ATT", &shifted).consequences.as_slice(),
            [Consequence::MissenseVariant]
        );

        let one_place_on = row(&tx, 159, b"TTA", &shifted);
        assert_eq!(
            one_place_on.hgvsc.as_deref(),
            Some("NM_000004.1:c.12_13insTTA")
        );
        assert_eq!(one_place_on.hgvs_offset, None);
        let moved_off_the_seam = row(&tx, 162, b"TC", &shifted);
        assert_eq!(
            moved_off_the_seam.hgvsc.as_deref(),
            Some("NM_000004.1:c.16_17insCT")
        );
        assert_eq!(moved_off_the_seam.hgvs_offset, Some(1));

        let annotated_where_shifted = EffectsConfig {
            enable_indel_3prime_shift: true,
            hgvs: crate::hgvs::HgvsOptions {
                shift_3prime: true,
                ..Default::default()
            },
            ..hgvs_config(&fasta, false)
        };
        let tc = row(&tx, 158, b"ATT", &annotated_where_shifted);
        assert_eq!(tc.cds_position.as_deref(), Some("12-13"));
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000004.1:c.12_13insTTA"));

        let mut failed = coding_seam_model();
        failed.bam_edit_status = Some(BamEditStatus::Failed);
        let tc = row(&failed, 158, b"ATT", &shifted);
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000004.1:c.12_13insTTA"));
        assert_eq!(tc.hgvs_offset, Some(1));
    }

    /// An insertion the genome's 3' shift moves off the bases an edited model
    /// holds and the genome lacks has its HGVSp read on those bases where
    /// they open the coding sequence. `map_insert` ends the insertion's first
    /// half on the base before them (Mapper.pm:522-542), a base of the 5'
    /// UTR, so the row's `CDS_position` is `?-9` and its own span has no
    /// translation start (TranscriptMapper.pm:441-443), although every base
    /// it stands on is coding.
    ///
    /// The model holds `ATGCCACGT` as CDS 1-9 before the genome's `ATGG` of
    /// CDS 10-13, and the allele is those nine bases twice. The genome's
    /// shift moves it three bases, past them, where Ensembl VEP reads both
    /// notations (TranscriptVariationAllele.pm:1441, 1711): the whole allele
    /// added after the model's bases, 18 bases and 6 residues more than the
    /// model has. In place of the model's nine bases the allele adds nine,
    /// the duplication of CDS 1-9 that the 3' rule places as `c.4_12dup`,
    /// and three residues, the duplication of `Pro2` to `Met4`. Neither
    /// notation moved, so the row carries no `HGVS_OFFSET`.
    ///
    /// The same record under `--shift_hgvs 0`, and an insertion on the same
    /// bases that the shift leaves in place, are read at the row's own span
    /// and have no HGVSp (1730-1738).
    #[test]
    fn insertion_the_shift_moves_off_bases_that_open_the_coding_sequence_has_its_hgvsp_read_there()
    {
        let (_dir, fasta) = coding_seam_genome();
        let fasta = Arc::new(fasta);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        let tx = coding_model_holding(0, "ATGCCACGT");
        let row = |alt: &[u8], config: &EffectsConfig| {
            let variant = InputVariant::new("1".into(), 150, 149, b"-".to_vec(), alt.to_vec());
            calculate_consequences(&variant, &tx, config)
                .expect("an insertion on the first coding bases has a row")
        };

        let twice = b"ATGCCACGTATGCCACGT";
        let tc = row(twice, &shifted);
        assert_eq!(tc.cds_position.as_deref(), Some("?-9"));
        assert_eq!(tc.hgvsc.as_deref(), Some("NM_000004.1:c.4_12dup"));
        assert_eq!(
            tc.hgvsp.as_deref(),
            Some("ENSP00000000001.1:p.Pro2_Met4dup")
        );
        assert_eq!(tc.hgvs_offset, None);

        let unshifted = row(twice, &hgvs_config(&fasta, false));
        assert_eq!(unshifted.hgvsc.as_deref(), Some("NM_000004.1:c.1_9dup"));
        assert_eq!(unshifted.hgvsp, None);

        let left_in_place = row(b"CCA", &shifted);
        assert_eq!(left_in_place.cds_position.as_deref(), Some("?-9"));
        assert_eq!(
            left_in_place.hgvsc.as_deref(),
            Some("NM_000004.1:c.1_9delinsCCA")
        );
        assert_eq!(left_in_place.hgvsp, None);
        assert_eq!(left_in_place.hgvs_offset, None);
    }

    /// An insertion beside the bases an edited model holds and the genome
    /// lacks, which the genome shows as a duplication of bases on both sides
    /// of them, is numbered on the model. `hgvs_variant_notation` types the
    /// duplication from the genomic bases before the point
    /// (Utils/Sequence.pm:577-585); numbered from the first copied base to
    /// the last, the range holds the model's own bases too, more positions
    /// than bases copied. On the model a duplication copies the bases beside
    /// the point, as many as the allele.
    ///
    /// The two forward-strand models carry their bases where `NM_015120.4`
    /// and `NM_002111.8` do, CDS 75-77 and CDS 111-116, each ending in the
    /// base before them, `C`, with `T` and `G` after. `CT` inserted one base
    /// past them copies, on the genome, the `C` before them and the `T`
    /// after; on the model it copies their last base and that `T`:
    /// `c.77_78dup` and `c.116_117dup`, two positions for two bases, the
    /// ranges Ensembl VEP 116.2 writes on those models
    /// (TranscriptVariationAllele.pm:1524-1525 moves both ends by one
    /// offset). Written two bases earlier, the same allele is carried there
    /// by the shift and keeps its offset. Where the model's bases end in `A`
    /// the allele copies nothing on the model and is an `ins` at the point.
    /// A duplication of two bases past them is untouched.
    ///
    /// Under `--shift_hgvs 0` the genome reads the bases after the point
    /// first: `CTG` before CDS 74 copies CDS 74 and the two genomic bases
    /// after the model's, which the model does not hold side by side, and
    /// then the three bases before the point, CDS 71-73, which it does.
    ///
    /// On the reverse strand the copied bases are cDNA 20 and 24, and the
    /// model shows the duplication only where cDNA 23 is the base of cDNA 20.
    /// A model without a cached cDNA cannot show it and has no notation. A
    /// model whose edit failed is read against the genome and keeps the
    /// genome's range.
    ///
    /// Written on the model's bases, the same forms follow the rule of an
    /// insertion the shift moves off them: on `ACC` at CDS 9-11, `AG` and
    /// `AGG` shift one base, where the genome duplicates CDS 8 and 12 or CDS
    /// 7, 8 and 12, and are described in place of `ACC`.
    #[test]
    fn duplication_across_edited_model_bases_is_numbered_on_the_model() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        use vep_core::transcript::BamEditStatus;
        let (_dir, fasta) = coding_seam_genome();
        let fasta = Arc::new(fasta);
        let unshifted = hgvs_config(&fasta, false);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        let row = |tx: &Transcript, start: u64, alt: &[u8], config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config).expect("an insertion in an exon has a row")
        };
        let hgvsc = |tx: &Transcript, start: u64, alt: &[u8]| row(tx, start, alt, &shifted).hgvsc;

        let three_held = coding_model_holding(74, "AAC");
        let six_held = coding_model_holding(110, "AAAAAC");
        assert_eq!(
            hgvsc(&three_held, 225, b"CT").as_deref(),
            Some("NM_000004.1:c.77_78dup")
        );
        assert_eq!(
            hgvsc(&six_held, 261, b"CT").as_deref(),
            Some("NM_000004.1:c.116_117dup")
        );
        let carried = row(&three_held, 223, b"CT", &shifted);
        assert_eq!(carried.hgvsc.as_deref(), Some("NM_000004.1:c.77_78dup"));
        assert_eq!(carried.hgvs_offset, Some(2));
        assert_eq!(
            hgvsc(&coding_model_holding(74, "AAA"), 225, b"CT").as_deref(),
            Some("NM_000004.1:c.78_79insCT")
        );
        assert_eq!(
            hgvsc(&three_held, 226, b"TG").as_deref(),
            Some("NM_000004.1:c.78_79dup")
        );
        assert_eq!(
            row(&three_held, 223, b"CTG", &unshifted).hgvsc.as_deref(),
            Some("NM_000004.1:c.71_73dup")
        );

        let (_dir, reverse_fasta) = hgvsc_fasta_fixture(&bases_from(1201, &"G".repeat(80)));
        let reverse_shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&Arc::new(reverse_fasta), false)
        };
        let on_reverse = |tx: &Transcript| row(tx, 1280, b"GT", &reverse_shifted).hgvsc;
        let mut reverse = reverse_non_coding_seam_model();
        assert_eq!(
            on_reverse(&reverse).as_deref(),
            Some("NR_000001.1:n.24_25insAC")
        );
        reverse.vefc.as_mut().unwrap().spliced_seq =
            Some(format!("{}CCA{}", "A".repeat(20), "C".repeat(180)));
        reverse.derived = Default::default();
        assert_eq!(
            on_reverse(&reverse).as_deref(),
            Some("NR_000001.1:n.23_24dup")
        );
        reverse.vefc.as_mut().unwrap().spliced_seq = None;
        reverse.derived = Default::default();
        assert_eq!(on_reverse(&reverse), None);

        let mut failed = coding_model_holding(74, "AAC");
        failed.bam_edit_status = Some(BamEditStatus::Failed);
        assert_eq!(
            hgvsc(&failed, 225, b"CT").as_deref(),
            Some("NM_000004.1:c.74_78dup")
        );

        let on_the_bases = coding_seam_model();
        assert_eq!(
            hgvsc(&on_the_bases, 158, b"AG").as_deref(),
            Some("NM_000004.1:c.10_11delinsG")
        );
        assert_eq!(
            hgvsc(&on_the_bases, 158, b"AGG").as_deref(),
            Some("NM_000004.1:c.10_11delinsGG")
        );
    }

    /// A record that keeps its reference allele whole, with bases on both
    /// sides of the bases an edited model holds and the genome lacks, spans
    /// those bases on the model: the row's `CDS_position` runs from the first
    /// reference base to the last and its reference is cut from the model
    /// over that span (TranscriptVariation.pm:463). An alternate that is the
    /// genomic reference twice is typed a duplication on the genome
    /// (Utils/Sequence.pm:602-606) and skips the clip
    /// (TranscriptVariationAllele.pm:1492); numbered through the model its
    /// range is the whole span, the model's own bases named as copied. On
    /// the model the alternate stands in place of the span.
    ///
    /// The forward-strand model carries `AGC` where `NM_015120.4` holds its
    /// bases, CDS 75-77, between the genomic `C` of CDS 74 and the `T` of
    /// CDS 78:
    ///
    /// - `CT` to `CTCT` stands in place of `CAGCT`, CDS 74-78: the shared
    ///   ends go and `T` replaces `AG`, one base fewer, the row's frameshift.
    /// - `CTG` to `CTGCTG` stands in place of CDS 74-79, and `GCT` to
    ///   `GCTGCT` of CDS 73-78: six bases for six that differ at CDS 75
    ///   alone, a change of one base that leaves codon 25 `Ala`.
    ///
    /// `TC` inserted between the two genomic bases is the first change in
    /// its minimal spelling and has its notation, with the shift on or off.
    /// A reference that ends before the model's bases, or starts after them,
    /// is a duplication on the model as on the genome. A model whose edit
    /// failed is read against the genome and keeps the genome's range.
    ///
    /// On the reverse strand the reference `GT` is cDNA 20 and 24 around the
    /// model's `CCC`: `ACAC` stands in place of `ACCCC`, and `A` replaces
    /// `CC`. A model without a cached cDNA cannot show its bases and has no
    /// notation.
    #[test]
    fn duplication_spelt_across_edited_model_bases_replaces_the_model_span() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        use vep_core::transcript::BamEditStatus;
        let (_dir, fasta) = coding_seam_genome();
        let fasta = Arc::new(fasta);
        let unshifted = hgvs_config(&fasta, false);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        let row = |tx: &Transcript, start: u64, reference: &[u8], config: &EffectsConfig| {
            let end = start + reference.len() as u64 - 1;
            let variant = InputVariant::new(
                "1".into(),
                start,
                end,
                reference.to_vec(),
                reference.repeat(2),
            );
            calculate_consequences(&variant, tx, config).expect("an allele in an exon has a row")
        };
        let inserted = |tx: &Transcript, start: u64, alt: &[u8], config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an insertion in an exon has a row")
                .hgvsc
        };

        let held = coding_model_holding(74, "AGC");
        for config in [&shifted, &unshifted] {
            let two_bases = row(&held, 223, b"CT", config);
            assert_eq!(two_bases.cds_position.as_deref(), Some("74-78"));
            assert_eq!(
                two_bases.consequences.as_slice(),
                [Consequence::FrameshiftVariant]
            );
            assert_eq!(
                two_bases.hgvsc.as_deref(),
                Some("NM_000004.1:c.75_76delinsT")
            );
            assert_eq!(
                two_bases.hgvsp.as_deref(),
                Some("ENSP00000000001.1:p.Ala26LeufsTer?")
            );
            assert_eq!(two_bases.hgvs_offset, None);
            assert_eq!(
                inserted(&held, 224, b"TC", config).as_deref(),
                Some("NM_000004.1:c.75_76delinsT")
            );
        }

        let one_base: [(u64, &[u8], &str); 2] = [(223, b"CTG", "74-79"), (222, b"GCT", "73-78")];
        for (start, reference, cds_position) in one_base {
            let three_bases = row(&held, start, reference, &shifted);
            assert_eq!(three_bases.cds_position.as_deref(), Some(cds_position));
            assert_eq!(
                three_bases.consequences.as_slice(),
                [Consequence::SynonymousVariant]
            );
            assert_eq!(three_bases.hgvsc.as_deref(), Some("NM_000004.1:c.75A>T"));
        }

        assert_eq!(
            row(&held, 221, b"TGC", &shifted).hgvsc.as_deref(),
            Some("NM_000004.1:c.72_74dup")
        );
        assert_eq!(
            row(&held, 224, b"TG", &shifted).hgvsc.as_deref(),
            Some("NM_000004.1:c.78_79dup")
        );

        let mut failed = coding_model_holding(74, "AGC");
        failed.bam_edit_status = Some(BamEditStatus::Failed);
        assert_eq!(
            row(&failed, 223, b"CT", &shifted).hgvsc.as_deref(),
            Some("NM_000004.1:c.74_78dup")
        );

        let (_dir, reverse_fasta) = hgvsc_fasta_fixture(&bases_from(1201, &"G".repeat(80)));
        let reverse_shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&Arc::new(reverse_fasta), false)
        };
        let mut reverse = reverse_non_coding_seam_model();
        let across = row(&reverse, 1280, b"GT", &reverse_shifted);
        assert_eq!(across.cdna_position.as_deref(), Some("20-24"));
        assert_eq!(across.hgvsc.as_deref(), Some("NR_000001.1:n.22_23delinsA"));
        assert_eq!(
            inserted(&reverse, 1281, b"TG", &reverse_shifted).as_deref(),
            Some("NR_000001.1:n.22_23delinsA")
        );
        reverse.vefc.as_mut().unwrap().spliced_seq = None;
        reverse.derived = Default::default();
        assert_eq!(row(&reverse, 1280, b"GT", &reverse_shifted).hgvsc, None);
    }

    /// The span a duplication spelt across an edited model's bases replaces
    /// is bounded by the exon, not by the mapper pairs the genome joins end
    /// to end. An alignment that also drops a genome base splits the exon's
    /// pairs once more, so a reference kept whole across the model's bases
    /// and the dropped base has its two ends in two runs of pairs and in one
    /// exon. The model's span between them holds the bases the model adds
    /// and lacks the base it drops, and the alternate, the genomic reference
    /// twice, stands in place of it.
    ///
    /// The forward-strand model carries `AGC` at CDS 75-77 and drops genomic
    /// 227, the `T` after CDS 80. The six bases of genomic 223-228 twice
    /// stand in place of `CAGCTGCG`, CDS 74-81: the first and the last base
    /// are shared, and ten bases replace CDS 75-80. Where the model carries
    /// one base, `A`, its span `CATGCG` is as long as the reference and the
    /// alternate twice as long as both, so only the bases show that the
    /// alternate is no copy of the model's span: ten bases replace CDS
    /// 75-78.
    ///
    /// A reference whose last base is the dropped one, and a reference that
    /// runs through the intron into exon 2, keep the duplication typed on
    /// the genome.
    ///
    /// On the reverse strand the model carries `CCC` at cDNA 21-23 and drops
    /// genomic 1278: the reference of genomic 1277-1281 is cDNA 20-26,
    /// `ACCCCCC`, and `ACCCCACCCC` in its place leaves `ACC` inserted after
    /// cDNA 24.
    #[test]
    fn duplication_spelt_across_held_and_dropped_bases_replaces_the_span_in_one_exon() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        use vep_core::transcript::{ExonCoordMapper, MapperPair};
        let (_dir, fasta) = coding_seam_genome();
        let fasta = Arc::new(fasta);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        let row = |tx: &Transcript, start: u64, reference: &[u8], config: &EffectsConfig| {
            let end = start + reference.len() as u64 - 1;
            let variant = InputVariant::new(
                "1".into(),
                start,
                end,
                reference.to_vec(),
                reference.repeat(2),
            );
            calculate_consequences(&variant, tx, config).expect("an allele in an exon has a row")
        };

        let three_held = coding_model_holding_and_dropping(74, "AGC", 227);
        assert_eq!(
            three_held.spliced_seq().map(|s| &s[123..133]),
            Some("CAGCTGCGCT")
        );
        let across = row(&three_held, 223, b"CTGCTG", &shifted);
        assert_eq!(across.cds_position.as_deref(), Some("74-81"));
        assert_eq!(
            across.hgvsc.as_deref(),
            Some("NM_000004.1:c.75_80delinsTGCTGCTGCT")
        );
        assert_eq!(across.hgvs_offset, None);

        let one_held = coding_model_holding_and_dropping(74, "A", 227);
        let as_long = row(&one_held, 223, b"CTGCTG", &shifted);
        assert_eq!(as_long.cds_position.as_deref(), Some("74-79"));
        assert_eq!(
            as_long.hgvsc.as_deref(),
            Some("NM_000004.1:c.75_78delinsTGCTGCTGCT")
        );

        let typed_on_the_genome = |hgvsc: Option<String>| hgvsc.is_some_and(|h| h.ends_with("dup"));
        assert!(typed_on_the_genome(
            row(&three_held, 223, b"CTGCT", &shifted).hgvsc
        ));
        let into_exon_2 = fasta
            .sequence("1", 223, 2100)
            .expect("the fixture holds both exons");
        assert!(typed_on_the_genome(
            row(&three_held, 223, &into_exon_2, &shifted).hgvsc
        ));

        let (_dir, reverse_fasta) = hgvsc_fasta_fixture(&bases_from(1201, &"G".repeat(80)));
        let reverse_shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&Arc::new(reverse_fasta), false)
        };
        let reverse_pair = |from_start, from_end, to_start, to_end| MapperPair {
            ori: -1,
            ..pair(from_start, from_end, to_start, to_end)
        };
        let mut reverse = reverse_non_coding_seam_model();
        let vefc = reverse.vefc.as_mut().unwrap();
        vefc.mapper.as_mut().unwrap().exon_coord_mapper = ExonCoordMapper::new(vec![
            reverse_pair(1, 20, 1281, 1300),
            reverse_pair(24, 25, 1279, 1280),
            reverse_pair(26, 102, 1201, 1277),
            reverse_pair(103, 202, 1001, 1100),
        ]);
        vefc.spliced_seq = Some(format!("{}{}", "A".repeat(20), "C".repeat(182)));
        reverse.derived = Default::default();
        let on_reverse = row(&reverse, 1277, b"GGGGT", &reverse_shifted);
        assert_eq!(on_reverse.cdna_position.as_deref(), Some("20-26"));
        assert_eq!(
            on_reverse.hgvsc.as_deref(),
            Some("NR_000001.1:n.24_25insACC")
        );
    }

    /// One change spelt two ways reads one notation where the model also
    /// drops a genome base further along the exon: the reference kept whole
    /// and doubled, and the insertion alone. The 3' rule places what the
    /// clip leaves of either within the exon, past the dropped base, where
    /// the model's own sequence runs on unbroken.
    ///
    /// The forward-strand model carries one `A` after CDS 9, the first of
    /// the four `A` at genomic 158-161, and drops genomic 160, so CDS 9-12
    /// read `AAAA` before the `T` of CDS 13. `AA` added between genomic 158
    /// and 159 leaves one `A` inserted, the duplication of CDS 12, the last
    /// `A` of the run and the base after the dropped one. Where the model
    /// carries `AAA` there, the same two bases stand in place of three: one
    /// `A` is deleted, CDS 14, the last of that run. Four `A` added before
    /// genomic 158, off the model's base, copy genomic 158-161 on the genome
    /// and CDS 9-12 on the model, the range the doubled reference of those
    /// four bases names.
    ///
    /// On the reverse strand the model carries `CCC` at cDNA 21-23, drops
    /// genomic 1278 and reads `C` to the end of its first exon, cDNA 102:
    /// `CC` in place of the three deletes one `C`, the exon's last. A model
    /// whose exon list holds the bases on either side of its own in two
    /// exons keeps the mapper pairs the genome joins end to end as the
    /// bound, and the duplication stops before the dropped base.
    #[test]
    fn insertion_on_held_bases_is_placed_within_the_exon_as_its_doubled_reference_is() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        use vep_core::transcript::{ExonCoordMapper, MapperPair};
        let (_dir, fasta) = coding_seam_genome();
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&Arc::new(fasta), false)
        };
        let doubled = |tx: &Transcript, start: u64, reference: &[u8], config: &EffectsConfig| {
            let end = start + reference.len() as u64 - 1;
            let variant = InputVariant::new(
                "1".into(),
                start,
                end,
                reference.to_vec(),
                reference.repeat(2),
            );
            calculate_consequences(&variant, tx, config)
                .expect("an allele in an exon has a row")
                .hgvsc
        };
        let inserted = |tx: &Transcript, start: u64, alt: &[u8], config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an insertion in an exon has a row")
                .hgvsc
        };

        let one_held = coding_model_holding_and_dropping(9, "A", 160);
        let three_held = coding_model_holding_and_dropping(9, "AAA", 160);
        assert_eq!(one_held.spliced_seq().map(|s| &s[57..63]), Some("GAAAAT"));
        assert_eq!(
            three_held.spliced_seq().map(|s| &s[57..65]),
            Some("GAAAAAAT")
        );
        // The reference, a run of one base, is the allele its doubling adds.
        let spellings: [(&Transcript, u64, &[u8], u64, &str); 3] = [
            (&one_held, 158, b"AA", 159, "c.12dup"),
            (&three_held, 158, b"AA", 159, "c.14del"),
            (&one_held, 158, b"AAAA", 158, "c.9_12dup"),
        ];
        for (tx, whole_at, reference, inserted_at, notation) in spellings {
            let notation = format!("NM_000004.1:{notation}");
            assert_eq!(
                doubled(tx, whole_at, reference, &shifted).as_deref(),
                Some(notation.as_str())
            );
            assert_eq!(
                inserted(tx, inserted_at, reference, &shifted).as_deref(),
                Some(notation.as_str())
            );
        }

        let mut two_exons = coding_model_holding_and_dropping(9, "A", 160);
        let first_exon = two_exons.exons[0].clone();
        assert_eq!((first_exon.start, first_exon.end), (100, 399));
        two_exons.exons[0].end = 158;
        two_exons.exons.insert(
            1,
            vep_core::transcript::Exon {
                start: 159,
                ..first_exon
            },
        );
        assert_eq!(
            inserted(&two_exons, 159, b"AA", &shifted).as_deref(),
            Some("NM_000004.1:c.11dup")
        );

        let (_dir, reverse_fasta) = hgvsc_fasta_fixture(&bases_from(1201, &"G".repeat(100)));
        let reverse_shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&Arc::new(reverse_fasta), false)
        };
        let reverse_pair = |from_start, from_end, to_start, to_end| MapperPair {
            ori: -1,
            ..pair(from_start, from_end, to_start, to_end)
        };
        let mut reverse = reverse_non_coding_seam_model();
        let vefc = reverse.vefc.as_mut().unwrap();
        vefc.mapper.as_mut().unwrap().exon_coord_mapper = ExonCoordMapper::new(vec![
            reverse_pair(1, 20, 1281, 1300),
            reverse_pair(24, 25, 1279, 1280),
            reverse_pair(26, 102, 1201, 1277),
            reverse_pair(103, 202, 1001, 1100),
        ]);
        vefc.spliced_seq = Some("C".repeat(202));
        reverse.derived = Default::default();
        assert_eq!(
            doubled(&reverse, 1280, b"GG", &reverse_shifted).as_deref(),
            Some("NR_000001.1:n.102del")
        );
        assert_eq!(
            inserted(&reverse, 1281, b"GG", &reverse_shifted).as_deref(),
            Some("NR_000001.1:n.102del")
        );
    }

    /// The HGVSc of an allele that keeps a base on a model the genome is the
    /// reference of carries the offset of its move. Placed by the 3' rule
    /// along the genome's bases, it is the notation the same bases take
    /// written as a minimised insertion, a record that shifts and prints
    /// its offset: `A` after the kept `G` of CDS 8 moves four bases to the
    /// last `A` of the genome's run, `AC` and `AG` one base. `C` stays at
    /// the point, `G` copies the kept base where it stands, a base kept
    /// after the point leaves `C` there, and an `A` written before a kept
    /// `A` duplicates that base as both engines type it: none has an
    /// offset. The record itself does not shift, so the fields read from it
    /// stay where it is written: of three kept bases the row's reference
    /// allele is the three as given, which a shift of one would turn. On
    /// the reverse strand the offset is negative: a `C` after the last `A`
    /// of the reverse model moves over the exon's 80 `C` towards lower
    /// coordinates. Without the 3' shift nothing moves, and a model the
    /// alignment corrected, whose bases are placed along its own, prints no
    /// offset. Ensembl VEP leaves the bases at the point, an `ins` with no
    /// offset (TranscriptVariationAllele.pm:1492; OutputFactory.pm:1750-1758).
    #[test]
    fn a_notation_moved_beside_a_kept_base_on_a_failed_edit_model_carries_its_offset() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        use vep_core::transcript::BamEditStatus;
        let genome = [bases_from(150, "ATGGCTGG"), bases_from(158, "AAAATTCGAT")].concat();
        let (_dir, fasta) = hgvsc_fasta_fixture(&genome);
        let fasta = Arc::new(fasta);
        let unshifted = hgvs_config(&fasta, true);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, true)
        };
        let failed = twice_read_base_model(BamEditStatus::Failed);
        let corrected = twice_read_base_model(BamEditStatus::Ok);
        let row = |tx: &Transcript,
                   start: u64,
                   end: u64,
                   reference: &[u8],
                   alt: &[u8],
                   config: &EffectsConfig| {
            let variant =
                InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            calculate_consequences(&variant, tx, config)
                .expect("an allele inside an exon has a row")
        };
        // The inserted bases, their notation and the offset of its move.
        let moved: [(&[u8], &str, Option<i64>); 5] = [
            (b"A", "c.13dup", Some(4)),
            (b"AC", "c.10_11insCA", Some(1)),
            (b"AG", "c.8_10dup", Some(1)),
            (b"C", "c.8_9insC", None),
            (b"G", "c.8dup", None),
        ];
        for (bases, notation, offset) in moved {
            let label = String::from_utf8_lossy(bases);
            let minimised = row(&failed, 158, 157, b"-", bases, &shifted);
            let kept_before = [b"G".as_slice(), bases].concat();
            let kept = row(&failed, 157, 157, b"G", &kept_before, &shifted);
            let expected = Some(format!("NM_000006.1:{notation}"));
            assert_eq!(kept.hgvsc, expected, "kept before {label}");
            assert_eq!(minimised.hgvsc, expected, "minimised {label}");
            assert_eq!(kept.hgvs_offset, offset, "kept before {label}");
            assert_eq!(minimised.hgvs_offset, offset, "minimised {label}");
            let without_the_shift = row(&failed, 157, 157, b"G", &kept_before, &unshifted);
            assert_eq!(without_the_shift.hgvs_offset, None, "unshifted {label}");
            assert_eq!(kept.used_ref, without_the_shift.used_ref, "{label}");
            assert_eq!(kept.cds_position, without_the_shift.cds_position, "{label}");
        }
        let kept_after = row(&failed, 158, 158, b"A", b"CA", &shifted);
        assert_eq!(kept_after.hgvsc.as_deref(), Some("NM_000006.1:c.8_9insC"));
        assert_eq!(kept_after.hgvs_offset, None);
        let typed_on_the_kept_base = row(&failed, 158, 158, b"A", b"AA", &shifted);
        assert_eq!(
            typed_on_the_kept_base.hgvsc.as_deref(),
            Some("NM_000006.1:c.10dup")
        );
        assert_eq!(typed_on_the_kept_base.hgvs_offset, None);

        let three_kept = row(&failed, 155, 157, b"TGG", b"TGGAC", &shifted);
        assert_eq!(
            three_kept.hgvsc.as_deref(),
            Some("NM_000006.1:c.10_11insCA")
        );
        assert_eq!(three_kept.hgvs_offset, Some(1));
        assert!(three_kept.hgvsp.is_some());
        assert_eq!(three_kept.used_ref.as_deref(), Some("TGG"));
        assert_eq!(three_kept.cds_position.as_deref(), Some("6-8"));

        let mut reverse = reverse_non_coding_seam_model();
        reverse.bam_edit_status = Some(BamEditStatus::Failed);
        reverse.derived = Default::default();
        let on_the_reverse_strand = row(&reverse, 1281, 1281, b"T", b"GT", &shifted);
        assert_eq!(
            on_the_reverse_strand.hgvsc.as_deref(),
            Some("NR_000001.1:n.103dup")
        );
        assert_eq!(on_the_reverse_strand.hgvs_offset, Some(-80));

        let on_the_corrected_model = row(&corrected, 157, 157, b"G", b"GAC", &shifted);
        assert_eq!(
            on_the_corrected_model.hgvsc.as_deref(),
            Some("NM_000006.1:c.9_10insCA")
        );
        assert_eq!(on_the_corrected_model.hgvs_offset, None);
    }

    /// A non-coding `NR_` model of two exons on `strand` with 100 genomic bases
    /// between them, as the cache holds a model edited by `21 20 GCA`: the first
    /// exon is cDNA 1-20, the model's own `GCA` at cDNA 21-23 and its last
    /// seven bases, `TGATTCG`, at cDNA 24-30; the second exon, which begins
    /// `AC`, follows the bases `own_after` that the model alone holds directly
    /// after cDNA 30. On the forward strand the first exon is genomic 1001-1027
    /// with the model's bases between 1020 and 1021, the intron 1028-1127 and
    /// the second exon 1128-1227; on the reverse strand the first exon is
    /// 1227-1201 with the model's bases between 1208 and 1207, the intron
    /// 1200-1101 and the second exon 1100-1001.
    fn seam_before_an_intron_model(strand: Strand, own_after: &str) -> Transcript {
        use vep_core::transcript::{Attribute, Exon, ExonCoordMapper, MapperPair};
        let mut tx = non_coding_reverse_two_exon_transcript();
        tx.stable_id = "NR_000002.1".into();
        tx.strand = strand;
        tx.start = 1001;
        tx.end = 1227;
        tx.attributes.push(Attribute {
            code: "_rna_edit".into(),
            value: "21 20 GCA".to_string(),
        });
        let own = own_after.len() as u64;
        let (first_exon, second_exon, pairs) = if strand == Strand::Forward {
            (
                (1001, 1027),
                (1128, 1227),
                [
                    (1, 20, 1001, 1020),
                    (24, 30, 1021, 1027),
                    (31 + own, 130 + own, 1128, 1227),
                ],
            )
        } else {
            (
                (1201, 1227),
                (1001, 1100),
                [
                    (1, 20, 1208, 1227),
                    (24, 30, 1201, 1207),
                    (31 + own, 130 + own, 1001, 1100),
                ],
            )
        };
        let exon = |(start, end): (u64, u64), rank| Exon {
            start,
            end,
            rank,
            phase: -1,
            end_phase: -1,
        };
        tx.exons = vec![exon(first_exon, 1), exon(second_exon, 2)];
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.sorted_exons = tx.exons.clone();
        vefc.mapper.as_mut().unwrap().exon_coord_mapper = ExonCoordMapper::new(
            pairs
                .into_iter()
                .map(|(from_start, from_end, to_start, to_end)| MapperPair {
                    from_start,
                    from_end,
                    to_start,
                    to_end,
                    ori: strand.as_i8(),
                })
                .collect(),
        );
        vefc.spliced_seq = Some(format!(
            "{}GCA{}{own_after}AC{}",
            "C".repeat(20),
            "TGATTCG",
            "G".repeat(98)
        ));
        tx.derived = Default::default();
        tx
    }

    /// The edits that make the intron of [`seam_before_an_intron_model`] begin
    /// `intron` on the transcript's strand.
    fn seam_model_intron_begins(strand: Strand, intron: &str) -> Vec<(usize, u8)> {
        use crate::concordance_tests::bases_from;
        if strand == Strand::Forward {
            return bases_from(1028, intron);
        }
        let forward = crate::coding::reverse_complement(intron.as_bytes());
        bases_from(
            1201 - intron.len(),
            std::str::from_utf8(&forward).expect("ASCII bases"),
        )
    }

    /// A `T` reference with `edits` and two effects configurations computing
    /// HGVS against it, the first with the 3' shift and the second without.
    fn hgvs_configs_with_and_without_the_shift(
        edits: &[(usize, u8)],
    ) -> (tempfile::TempDir, EffectsConfig, EffectsConfig) {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let (dir, fasta) = hgvsc_fasta_fixture(edits);
        let fasta = Arc::new(fasta);
        let shifted = EffectsConfig {
            hgvs: crate::hgvs::HgvsOptions::default(),
            ..hgvs_config(&fasta, false)
        };
        (dir, shifted, hgvs_config(&fasta, false))
    }

    /// How a record spells bases added at a point: minimised, or keeping the
    /// genomic base before the point, the one after it, or both, as an allele
    /// of a record with a second alternate does.
    #[derive(Debug, Clone, Copy)]
    enum Spelling {
        Minimised,
        KeepsBefore,
        KeepsAfter,
        KeepsBoth,
    }

    const SPELLINGS: [Spelling; 4] = [
        Spelling::Minimised,
        Spelling::KeepsBefore,
        Spelling::KeepsAfter,
        Spelling::KeepsBoth,
    ];

    /// `inserted`, on the transcript's strand, standing on the bases
    /// [`seam_before_an_intron_model`] holds on `strand`, in one of its four
    /// spellings; the genomic bases beside the point are `T`.
    fn on_the_seam_models_bases(
        strand: Strand,
        inserted: &str,
        spelling: Spelling,
    ) -> InputVariant {
        let (before, allele) = if strand == Strand::Forward {
            (1020, inserted.as_bytes().to_vec())
        } else {
            (1207, crate::coding::reverse_complement(inserted.as_bytes()))
        };
        let (keeps_before, keeps_after) = match spelling {
            Spelling::Minimised => {
                return InputVariant::new("1".into(), before + 1, before, b"-".to_vec(), allele);
            }
            Spelling::KeepsBefore => (true, false),
            Spelling::KeepsAfter => (false, true),
            Spelling::KeepsBoth => (true, true),
        };
        let kept = |keeps: bool| {
            if keeps {
                b"T".as_slice()
            } else {
                b"".as_slice()
            }
        };
        let reference = [kept(keeps_before), kept(keeps_after)].concat();
        let alt = [kept(keeps_before), &allele, kept(keeps_after)].concat();
        let start = before + u64::from(!keeps_before);
        let end = before + u64::from(keeps_after);
        InputVariant::new("1".into(), start, end, reference, alt)
    }

    /// An HGVSc with its inserted bases replaced by their count.
    fn hgvsc_shape(hgvsc: Option<String>) -> Option<String> {
        hgvsc.map(|text| match text.rsplit_once("ins") {
            Some((head, bases)) => format!("{head}ins{}", bases.len()),
            None => text,
        })
    }

    /// The 3' rule carries what the clip leaves of a replacement of model bases
    /// past the exon's last base onto the genome's intron, where the model
    /// numbers the next exon's first base directly after it: HGVS excepts from
    /// the rule a change around an exon/exon junction, never one around an
    /// exon/intron junction. The model's first exon ends with its own `GCA` and
    /// then `TGATTCG`, before an intron that begins `ACG`. The allele adding
    /// the seven and `A` after the model's three copies cDNA 24 to the
    /// intron's first base; the one adding the seven and `AC` copies cDNA 24
    /// to the intron's second base; the one adding the seven and `AG`, which
    /// copies nothing, is nine bases inserted after the intron's first base.
    /// Each holds on either strand and in each of the four spellings of the
    /// record: one change takes one notation. The walk ends on the exon's last
    /// base, the change a copy ending there, where the intron begins `C`
    /// though the next exon begins `A`. Under `--shift_hgvs 0` nothing moves.
    #[test]
    fn a_copy_ending_on_an_exons_last_base_is_carried_into_the_intron_on_an_edited_model() {
        let row = |tx: &Transcript, variant: &InputVariant, config: &EffectsConfig| {
            calculate_consequences(variant, tx, config).expect("an allele inside an exon has a row")
        };
        let one = format!("GCA{}A", "TGATTCG");
        let two = format!("GCA{}AC", "TGATTCG");
        let uncopied = format!("GCA{}AG", "TGATTCG");

        for strand in [Strand::Forward, Strand::Reverse] {
            let model = seam_before_an_intron_model(strand, "");
            let (_ending_dir, ending, _) =
                hgvs_configs_with_and_without_the_shift(&seam_model_intron_begins(strand, "CCGG"));
            let (_continuing_dir, continuing, unshifted) =
                hgvs_configs_with_and_without_the_shift(&seam_model_intron_begins(strand, "ACGG"));
            for spelling in SPELLINGS {
                let stays: [(&EffectsConfig, &str, &str); 3] = [
                    (&ending, &one, "NR_000002.1:n.23_30dup"),
                    (&unshifted, &one, "NR_000002.1:n.23_24ins8"),
                    (&unshifted, &two, "NR_000002.1:n.23_24ins9"),
                ];
                for (case, (config, inserted, stayed)) in stays.into_iter().enumerate() {
                    let variant = on_the_seam_models_bases(strand, inserted, spelling);
                    assert_eq!(
                        hgvsc_shape(row(&model, &variant, config).hgvsc).as_deref(),
                        Some(stayed),
                        "{strand:?} strand, {spelling:?}, case {case} of those that stay"
                    );
                }
                let carried = [
                    (&one, "NR_000002.1:n.24_30+1dup", (0, 1)),
                    (&two, "NR_000002.1:n.24_30+2dup", (0, 2)),
                    (&uncopied, "NR_000002.1:n.30+1_30+2ins9", (1, 2)),
                ];
                for (inserted, into_the_intron, offsets) in carried {
                    let variant = on_the_seam_models_bases(strand, inserted, spelling);
                    let tc = row(&model, &variant, &continuing);
                    assert_eq!(
                        hgvsc_shape(tc.hgvsc).as_deref(),
                        Some(into_the_intron),
                        "{strand:?} strand, {spelling:?}, {} bases after the model's",
                        inserted.len() - 3
                    );
                    assert_eq!(tc.hgvs_intron_offsets, Some(offsets));
                }
            }
        }
    }

    /// The walk goes on past the exon's last base for an insertion alone: a
    /// deletion or a duplication of the model's own bases ends on that base
    /// though its run goes on after it. On the reverse model whose cDNA is `C`
    /// from 21 on, before an intron that begins `CCA`, a `C` deleted from the
    /// run is the exon's last base, 103, and the model's three bases
    /// duplicated are cDNA 101 to 103; so are they on the model that numbers
    /// five bases of its own, `C` as well, between that exon and the next.
    #[test]
    fn a_deletion_or_a_duplication_of_model_bases_ends_on_the_exons_last_base() {
        use vep_core::transcript::{ExonCoordMapper, MapperPair};
        let (_dir, shifted, _) =
            hgvs_configs_with_and_without_the_shift(&[(1199, b'G'), (1200, b'G')]);
        let before_the_intron = reverse_non_coding_seam_model();
        let mut before_its_own_bases = reverse_non_coding_seam_model();
        let vefc = before_its_own_bases.vefc.as_mut().unwrap();
        vefc.mapper.as_mut().unwrap().exon_coord_mapper = ExonCoordMapper::new(
            [
                (1, 20, 1281, 1300),
                (24, 103, 1201, 1280),
                (109, 208, 1001, 1100),
            ]
            .into_iter()
            .map(|(from_start, from_end, to_start, to_end)| MapperPair {
                from_start,
                from_end,
                to_start,
                to_end,
                ori: -1,
            })
            .collect(),
        );
        vefc.spliced_seq = Some(format!("{}{}", "A".repeat(20), "C".repeat(188)));
        before_its_own_bases.derived = Default::default();

        for (model, tx) in [
            ("an intron that continues the run", &before_the_intron),
            ("bases of the model's own", &before_its_own_bases),
        ] {
            let hgvsc = |alt: &[u8]| {
                let variant =
                    InputVariant::new("1".into(), 1281, 1280, b"-".to_vec(), alt.to_vec());
                calculate_consequences(&variant, tx, &shifted)
                    .expect("an insertion inside an exon has a row")
                    .hgvsc
            };
            assert_eq!(
                hgvsc(b"GG").as_deref(),
                Some("NR_000001.1:n.103del"),
                "a deletion before {model}"
            );
            assert_eq!(
                hgvsc(b"GGGGGG").as_deref(),
                Some("NR_000001.1:n.101_103dup"),
                "a duplication before {model}"
            );
        }
    }

    /// Where the model numbers bases of its own directly after the exon's last
    /// base, the 3' rule reads those bases and not the genome's intron: the
    /// notation is written against the model's sequence. On the model that
    /// holds `ACGTT` there, the allele adding the exon's last seven bases and
    /// `A` after the model's three copies cDNA 24 to 31, the first base after
    /// the exon's last; the one adding the seven and `AC` copies cDNA 24 to
    /// 32; the one adding the seven and `AG`, which copies nothing, is nine
    /// bases inserted between cDNA 31 and 32. Each holds on either strand and
    /// in each of the four spellings of the record, whether the intron begins
    /// `ACG` or `CCG`: one change takes one notation. Under
    /// `--shift_hgvs 0` nothing moves.
    #[test]
    fn a_copy_ending_on_an_exons_last_base_is_carried_onto_the_bases_the_model_holds_after_it() {
        let written = |tx: &Transcript, variant: &InputVariant, config: &EffectsConfig| {
            let tc = calculate_consequences(variant, tx, config)
                .expect("an allele inside an exon has a row");
            (hgvsc_shape(tc.hgvsc), tc.hgvs_intron_offsets)
        };
        let one = format!("GCA{}A", "TGATTCG");
        let two = format!("GCA{}AC", "TGATTCG");
        let uncopied = format!("GCA{}AG", "TGATTCG");
        let cases = [
            (&one, "NR_000002.1:n.24_31dup", "NR_000002.1:n.23_24ins8"),
            (&two, "NR_000002.1:n.24_32dup", "NR_000002.1:n.23_24ins9"),
            (
                &uncopied,
                "NR_000002.1:n.31_32ins9",
                "NR_000002.1:n.23_24ins9",
            ),
        ];

        for strand in [Strand::Forward, Strand::Reverse] {
            let own_after = seam_before_an_intron_model(strand, "ACGTT");
            for intron in ["ACGG", "CCGG"] {
                let (_dir, shifted, unshifted) = hgvs_configs_with_and_without_the_shift(
                    &seam_model_intron_begins(strand, intron),
                );
                for spelling in SPELLINGS {
                    for (inserted, onto_the_models_bases, unmoved) in cases {
                        let variant = on_the_seam_models_bases(strand, inserted, spelling);
                        assert_eq!(
                            written(&own_after, &variant, &shifted),
                            (Some(onto_the_models_bases.to_string()), Some((0, 0))),
                            "{strand:?} strand, {spelling:?}, intron begins {intron}, {} bases \
                             after the model's",
                            inserted.len() - 3
                        );
                        assert_eq!(
                            written(&own_after, &variant, &unshifted).0.as_deref(),
                            Some(unmoved)
                        );
                    }
                }
            }
        }
    }

    /// An insertion the walk leaves directly after an exon's last base lies
    /// between that base and the intron's first position, where the model
    /// numbers the next exon's first base directly after the exon's last: the
    /// next cDNA number is the next exon's first base, the intron's length
    /// further on the genome. On the model whose exon ends `TGATTCG` at cDNA
    /// 30, the allele adding the seven and `AC` after the model's three is
    /// nine bases inserted there where the intron begins `C`, though the next
    /// exon begins `AC`, and the one adding the seven and `G` is eight,
    /// whichever way the intron begins, on either strand and in each of the
    /// four spellings of the record. Where the model holds `ACGTT` of its own
    /// directly after cDNA 30, the eight are inserted between cDNA 30 and 31,
    /// the first of them.
    #[test]
    fn an_insertion_left_after_an_exons_last_base_is_written_before_the_introns_first_position() {
        let written = |tx: &Transcript, variant: &InputVariant, config: &EffectsConfig| {
            let tc = calculate_consequences(variant, tx, config)
                .expect("an allele inside an exon has a row");
            (hgvsc_shape(tc.hgvsc), tc.hgvs_intron_offsets)
        };
        let next_exon_begins = format!("GCA{}AC", "TGATTCG");
        let other = format!("GCA{}G", "TGATTCG");

        for strand in [Strand::Forward, Strand::Reverse] {
            let model = seam_before_an_intron_model(strand, "");
            let own_after = seam_before_an_intron_model(strand, "ACGTT");
            let (_ending_dir, ending, _) =
                hgvs_configs_with_and_without_the_shift(&seam_model_intron_begins(strand, "CCGG"));
            let (_continuing_dir, continuing, _) =
                hgvs_configs_with_and_without_the_shift(&seam_model_intron_begins(strand, "ACGG"));
            let cases = [
                (
                    &model,
                    &ending,
                    &next_exon_begins,
                    "NR_000002.1:n.30_30+1ins9",
                    (0, 1),
                ),
                (&model, &ending, &other, "NR_000002.1:n.30_30+1ins8", (0, 1)),
                (
                    &model,
                    &continuing,
                    &other,
                    "NR_000002.1:n.30_30+1ins8",
                    (0, 1),
                ),
                (
                    &own_after,
                    &ending,
                    &other,
                    "NR_000002.1:n.30_31ins8",
                    (0, 0),
                ),
                (
                    &own_after,
                    &continuing,
                    &other,
                    "NR_000002.1:n.30_31ins8",
                    (0, 0),
                ),
            ];
            for spelling in SPELLINGS {
                for (case, (tx, config, inserted, after_the_exon, offsets)) in
                    cases.into_iter().enumerate()
                {
                    let variant = on_the_seam_models_bases(strand, inserted, spelling);
                    assert_eq!(
                        written(tx, &variant, config),
                        (Some(after_the_exon.to_string()), Some(offsets)),
                        "{strand:?} strand, {spelling:?}, case {case}"
                    );
                }
            }
        }
    }

    /// A four-exon transcript without a translation on `strand`: exons
    /// 10_000-10_100, 10_201-10_300, 10_401-10_500 and 10_601-10_700, the three
    /// 100-base introns between them, exons and introns numbered in transcript
    /// order.
    fn four_exon_non_coding_transcript(strand: vep_core::coordinate::Strand) -> Transcript {
        non_coding_transcript(
            strand,
            &[
                (10_000, 10_100),
                (10_201, 10_300),
                (10_401, 10_500),
                (10_601, 10_700),
            ],
        )
    }

    /// A transcript without a translation on `strand` over `exons`, given in
    /// genomic order, with the introns between them; exons and introns are
    /// numbered in transcript order.
    fn non_coding_transcript(
        strand: vep_core::coordinate::Strand,
        exons: &[(u64, u64)],
    ) -> Transcript {
        use vep_core::coordinate::Strand;
        use vep_core::transcript::{
            Exon, ExonCoordMapper, Intron, MapperPair, TranscriptMapper, TranscriptVEFC,
        };
        let (tx_start, tx_end) = (exons[0].0, exons[exons.len() - 1].1);
        let mut exon_spans = exons.to_vec();
        let mut intron_spans: Vec<(u64, u64)> = exons
            .windows(2)
            .map(|pair| (pair[0].1 + 1, pair[1].0 - 1))
            .collect();
        if strand == Strand::Reverse {
            exon_spans.reverse();
            intron_spans.reverse();
        }
        let mut exons = Vec::new();
        let mut pairs = Vec::new();
        let mut cdna = 1u64;
        for (i, &(start, end)) in exon_spans.iter().enumerate() {
            exons.push(Exon {
                start,
                end,
                rank: i as u32 + 1,
                phase: -1,
                end_phase: -1,
            });
            pairs.push(MapperPair {
                from_start: cdna,
                from_end: cdna + (end - start),
                to_start: start,
                to_end: end,
                ori: strand.as_i8(),
            });
            cdna += end - start + 1;
        }
        let introns: Vec<Intron> = intron_spans
            .iter()
            .enumerate()
            .map(|(i, &(start, end))| Intron {
                start,
                end,
                rank: i as u32 + 1,
            })
            .collect();
        let mut tx = make_test_transcript();
        tx.strand = strand;
        tx.start = tx_start;
        tx.end = tx_end;
        tx.biotype = "lncRNA".into();
        tx.exons = exons.clone();
        tx.introns = introns.clone();
        tx.cdna_coding_start = None;
        tx.cdna_coding_end = None;
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation_start = None;
        tx.translation_end = None;
        tx.translation = None;
        tx.protein_id = None;
        tx.vefc = Some(TranscriptVEFC {
            codon_table: 1,
            introns,
            sorted_exons: exons,
            mapper: Some(TranscriptMapper {
                start_phase: 0,
                cdna_coding_start: 0,
                cdna_coding_end: 0,
                exon_coord_mapper: ExonCoordMapper::new(pairs),
            }),
            ..Default::default()
        });
        tx.derived = Default::default();
        tx
    }

    /// A base named by its offset from one junction: `In(k)` is the k-th
    /// intronic base counted into the intron from the junction (0 is the first
    /// or last intronic base), `Ex(j)` the j-th exonic base counted away from
    /// it (1 is the exonic base beside the junction).
    #[derive(Debug, Clone, Copy)]
    enum At {
        In(i64),
        Ex(i64),
    }

    /// An allele placed by junction offsets: an SNV at a base, a deletion of the
    /// bases between two offsets inclusive, or an insertion between two
    /// adjacent bases.
    #[derive(Debug, Clone, Copy)]
    enum Shape {
        Snv(At),
        Del(At, At),
        Ins(At, At),
    }

    /// The sorted term names of the one row `variant` writes on `tx`.
    fn sorted_terms(
        tx: &Transcript,
        variant: &InputVariant,
        config: &EffectsConfig,
    ) -> Vec<String> {
        let tc = calculate_consequences(variant, tx, config).expect("a row");
        let mut got: Vec<String> = tc.consequences.iter().map(|c| c.to_string()).collect();
        got.sort();
        got
    }

    /// Every splice window of every junction on both strands, for an SNV, a
    /// deletion and an insertion, against the term sets the Perl predicates
    /// assign: the windows of `_intron_effects`
    /// (BaseTranscriptVariationAllele.pm:99-224: donor and acceptor at the two
    /// terminal intronic bases, the 5th base, the donor region at intronic bases
    /// 3 to 6, the polypyrimidine tract at the last 17 to 3, the intron body
    /// from base 3 to the third-last) and of `_intron_overlap`
    /// (VariationEffect.pm:87-108: intronic bases 3 to 8 and exonic bases 1 to
    /// 3, plus the four equalities an insertion satisfies at a junction and two
    /// bases inside it), read through the predicates of VariationEffect.pm:560-633
    /// (donor and acceptor by strand, the 5th base over the donor region, all
    /// five over `splice_region`) and the `include` gates of Utils/Constants.pm
    /// (the polypyrimidine tract needs no exon overlap). The context term is
    /// `non_coding_transcript_exon_variant` on an exon overlap, else
    /// `non_coding_transcript_variant` (VariationEffect.pm:495-520). The
    /// insertion rows read Perl's inverted differing region, so a window holds an
    /// insertion only when both flanks lie inside it, and the polypyrimidine
    /// window either flank. Each shape's endpoint intron is evaluated by the
    /// endpoint pass and met again by the whole-span pass, which the evaluation
    /// set skips.
    #[test]
    fn one_splice_pass_matches_two_on_every_intron_shape() {
        use vep_core::coordinate::Strand;
        use At::{Ex, In};
        use Shape::{Del, Ins, Snv};
        const N: &str = "non_coding_transcript_variant";
        const X: &str = "non_coding_transcript_exon_variant";
        const SR: &str = "splice_region_variant";
        const PPT: &str = "splice_polypyrimidine_tract_variant";
        const IV: &str = "intron_variant";
        const DONOR: &str = "splice_donor_variant";
        const ACCEPTOR: &str = "splice_acceptor_variant";
        const FIFTH: &str = "splice_donor_5th_base_variant";
        const REGION: &str = "splice_donor_region_variant";

        const DONOR_CASES: &[(Shape, &[&str])] = &[
            (Snv(Ex(4)), &[X]),
            (Snv(Ex(3)), &[SR, X]),
            (Snv(Ex(2)), &[SR, X]),
            (Snv(Ex(1)), &[SR, X]),
            (Snv(In(0)), &[DONOR, N]),
            (Snv(In(1)), &[DONOR, N]),
            (Snv(In(2)), &[REGION, IV, N]),
            (Snv(In(3)), &[REGION, IV, N]),
            (Snv(In(4)), &[FIFTH, IV, N]),
            (Snv(In(5)), &[REGION, IV, N]),
            (Snv(In(6)), &[SR, IV, N]),
            (Snv(In(7)), &[SR, IV, N]),
            (Snv(In(8)), &[IV, N]),
            (Snv(In(50)), &[IV, N]),
            (Del(Ex(2), In(1)), &[DONOR, X]),
            (Del(In(3), In(6)), &[FIFTH, IV, N]),
            (Del(In(6), In(9)), &[SR, IV, N]),
            (Ins(Ex(1), In(0)), &[SR, N]),
            (Ins(In(0), In(1)), &[DONOR, N]),
            (Ins(In(1), In(2)), &[SR, IV, N]),
            (Ins(In(2), In(3)), &[REGION, IV, N]),
            (Ins(In(3), In(4)), &[REGION, IV, N]),
            (Ins(In(4), In(5)), &[REGION, IV, N]),
            (Ins(In(5), In(6)), &[SR, IV, N]),
            (Ins(In(6), In(7)), &[SR, IV, N]),
            (Ins(In(7), In(8)), &[IV, N]),
        ];
        const ACCEPTOR_CASES: &[(Shape, &[&str])] = &[
            (Snv(Ex(4)), &[X]),
            (Snv(Ex(3)), &[SR, X]),
            (Snv(Ex(2)), &[SR, X]),
            (Snv(Ex(1)), &[SR, X]),
            (Snv(In(0)), &[ACCEPTOR, N]),
            (Snv(In(1)), &[ACCEPTOR, N]),
            (Snv(In(2)), &[SR, PPT, IV, N]),
            (Snv(In(5)), &[SR, PPT, IV, N]),
            (Snv(In(7)), &[SR, PPT, IV, N]),
            (Snv(In(8)), &[PPT, IV, N]),
            (Snv(In(12)), &[PPT, IV, N]),
            (Snv(In(16)), &[PPT, IV, N]),
            (Snv(In(17)), &[IV, N]),
            (Del(In(3), Ex(1)), &[ACCEPTOR, IV, X]),
            (Del(In(6), In(3)), &[SR, PPT, IV, N]),
            (Del(In(12), In(9)), &[PPT, IV, N]),
            (Ins(In(0), Ex(1)), &[SR, N]),
            (Ins(In(1), In(0)), &[ACCEPTOR, N]),
            (Ins(In(2), In(1)), &[SR, PPT, IV, N]),
            (Ins(In(3), In(2)), &[SR, PPT, IV, N]),
            (Ins(In(7), In(6)), &[SR, PPT, IV, N]),
            (Ins(In(8), In(7)), &[PPT, IV, N]),
            (Ins(In(17), In(16)), &[PPT, IV, N]),
            (Ins(In(18), In(17)), &[IV, N]),
        ];

        let config = EffectsConfig::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let tx = four_exon_non_coding_transcript(strand);
            for (s, e) in [(10_101i64, 10_200i64), (10_301, 10_400), (10_501, 10_600)] {
                // The junction base and the direction into the intron from it.
                let (donor, acceptor) = match strand {
                    Strand::Forward => ((s, 1i64), (e, -1i64)),
                    Strand::Reverse => ((e, -1i64), (s, 1i64)),
                };
                for ((base, into), cases) in [(donor, DONOR_CASES), (acceptor, ACCEPTOR_CASES)] {
                    let at = |a: At| -> u64 {
                        match a {
                            In(k) => (base + into * k) as u64,
                            Ex(j) => (base - into * j) as u64,
                        }
                    };
                    for &(shape, expected) in cases {
                        let variant = match shape {
                            Snv(p) => InputVariant::new(
                                "21".into(),
                                at(p),
                                at(p),
                                b"A".to_vec(),
                                b"C".to_vec(),
                            ),
                            Del(a, b) => {
                                let (lo, hi) = (at(a).min(at(b)), at(a).max(at(b)));
                                InputVariant::new(
                                    "21".into(),
                                    lo,
                                    hi,
                                    vec![b'A'; (hi - lo + 1) as usize],
                                    b"-".to_vec(),
                                )
                            }
                            Ins(a, b) => {
                                let (lo, hi) = (at(a).min(at(b)), at(a).max(at(b)));
                                assert_eq!(hi, lo + 1, "{shape:?}");
                                InputVariant::new(
                                    "21".into(),
                                    hi,
                                    lo,
                                    b"-".to_vec(),
                                    b"ACG".to_vec(),
                                )
                            }
                        };
                        let label = format!("{strand:?} intron {s}-{e} {shape:?}");
                        let mut want: Vec<&str> = expected.to_vec();
                        want.sort_unstable();
                        assert_eq!(sorted_terms(&tx, &variant, &config), want, "{label}");
                    }
                }
            }
        }
    }

    /// A second evaluation of an intron and suppression already made adds no
    /// term and no entry; the other suppression of the same intron is its own
    /// evaluation.
    #[test]
    fn splice_evaluations_run_each_intron_and_suppression_once() {
        use vep_core::coordinate::Strand;
        let tx = four_exon_non_coding_transcript(Strand::Forward);
        // Intron 1 is 10_101-10_200: a deletion of its bases 7 to 8 sits in the
        // intronic splice_region window.
        let variant = InputVariant::new("21".into(), 10_107, 10_108, b"AA".to_vec(), b"-".to_vec());
        let intron = (10_101, 10_200);
        let mut evaluated = SpliceEvaluations::new();
        let mut csq: ConsequenceList = SmallVec::new();
        add_intronic_non_insertion_splice_consequences(
            &mut csq,
            &tx,
            &variant,
            intron,
            false,
            &mut evaluated,
        );
        let once = csq.clone();
        assert_eq!(
            once.as_slice(),
            &[Consequence::IntronVariant, Consequence::SpliceRegionVariant]
        );
        add_intronic_non_insertion_splice_consequences(
            &mut csq,
            &tx,
            &variant,
            intron,
            false,
            &mut evaluated,
        );
        assert_eq!(csq, once);
        assert_eq!(evaluated.as_slice(), &[(10_101, 10_200, false)]);

        add_intronic_non_insertion_splice_consequences(
            &mut csq,
            &tx,
            &variant,
            intron,
            true,
            &mut evaluated,
        );
        assert_eq!(csq, once);
        assert_eq!(
            evaluated.as_slice(),
            &[(10_101, 10_200, false), (10_101, 10_200, true)]
        );

        let insertion =
            InputVariant::new("21".into(), 10_102, 10_101, b"-".to_vec(), b"ACG".to_vec());
        let mut evaluated = SpliceEvaluations::new();
        let mut csq: ConsequenceList = SmallVec::new();
        add_intronic_insertion_splice_consequences(
            &mut csq,
            &tx,
            &insertion,
            intron,
            &mut evaluated,
        );
        assert_eq!(csq.as_slice(), &[Consequence::SpliceDonorVariant]);
        add_intronic_insertion_splice_consequences(
            &mut csq,
            &tx,
            &insertion,
            intron,
            &mut evaluated,
        );
        assert_eq!(csq.as_slice(), &[Consequence::SpliceDonorVariant]);
        assert_eq!(evaluated.as_slice(), &[(10_101, 10_200, false)]);
    }

    /// The `splice_*` terms, sorted, of the one row `variant` writes on `tx`.
    fn splice_terms(tx: &Transcript, variant: &InputVariant) -> Vec<String> {
        sorted_terms(tx, variant, &EffectsConfig::default())
            .into_iter()
            .filter(|term| term.starts_with("splice_"))
            .collect()
    }

    /// A record replacing the `ref_len` bases from `position` by `alt_len`
    /// bases: the alternate allele holds the reference's base at each index of
    /// `same` and another base at every other index.
    fn replacement(position: u64, ref_len: usize, alt_len: usize, same: &[usize]) -> InputVariant {
        let alt = (0..alt_len)
            .map(|i| if same.contains(&i) { b'C' } else { b'G' })
            .collect();
        InputVariant::new(
            "21".into(),
            position,
            position + ref_len as u64 - 1,
            vec![b'C'; ref_len],
            alt,
        )
    }

    /// An insertion of three bases between `before` and the base after it.
    fn insertion_after(before: u64) -> InputVariant {
        InputVariant::new(
            "21".into(),
            before + 1,
            before,
            b"-".to_vec(),
            b"GGG".to_vec(),
        )
    }

    /// One record on the four-exon transcript of its strand: the position, the
    /// two alleles' lengths, the indices at which the alleles hold the same
    /// base, and the `splice_*` terms of its row.
    type SpliceCase = (
        vep_core::coordinate::Strand,
        u64,
        usize,
        usize,
        &'static [usize],
        &'static [&'static str],
    );

    /// Every case's row holds exactly its `splice_*` terms; a failure names
    /// each case that differs.
    fn assert_splice_terms(cases: &[SpliceCase]) {
        assert_case_terms(cases, splice_terms);
    }

    /// Every case's row holds exactly its terms, the splice terms and every
    /// other; a failure names each case that differs.
    fn assert_whole_rows(cases: &[SpliceCase]) {
        assert_case_terms(cases, |tx, variant| {
            sorted_terms(tx, variant, &EffectsConfig::default())
        });
    }

    fn assert_case_terms(
        cases: &[SpliceCase],
        terms: impl Fn(&Transcript, &InputVariant) -> Vec<String>,
    ) {
        let wrong: Vec<String> = cases
            .iter()
            .filter_map(|&(strand, position, ref_len, alt_len, same, expected)| {
                let tx = four_exon_non_coding_transcript(strand);
                let got = terms(&tx, &replacement(position, ref_len, alt_len, same));
                (got != expected).then(|| {
                    format!(
                        "{strand:?} strand, {ref_len} bases from {position} replaced by \
                         {alt_len}: wrote {got:?}, expected {expected:?}"
                    )
                })
            })
            .collect();
        assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    }

    const DONOR_TERM: &str = "splice_donor_variant";
    const ACCEPTOR_TERM: &str = "splice_acceptor_variant";
    const FIFTH_BASE_TERM: &str = "splice_donor_5th_base_variant";
    const DONOR_REGION_TERM: &str = "splice_donor_region_variant";
    const TRACT_TERM: &str = "splice_polypyrimidine_tract_variant";
    const REGION_TERM: &str = "splice_region_variant";
    const INTRON_TERM: &str = "intron_variant";
    const EXON_TERM: &str = "non_coding_transcript_exon_variant";
    const TRANSCRIPT_TERM: &str = "non_coding_transcript_variant";

    /// An alternate allele longer than the reference allele is compared index
    /// by index to its own length, so its indices past the reference allele
    /// stand on the bases after the record's span
    /// (VariationFeatureOverlapAllele.pm:428-429, placed with no clamp at
    /// BaseTranscriptVariationAllele.pm:117-118). A donor term whose two bases
    /// only those positions reach is not written: the span replaces neither
    /// base. Intron 2 is 10_301-10_400; its donor bases are 10_301-10_302 on
    /// the forward strand and 10_399-10_400 on the reverse strand. What the
    /// span itself holds stands: the exonic or intronic `splice_region` window,
    /// or the donor region.
    #[test]
    fn a_donor_term_resting_on_bases_past_the_span_is_not_written() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        assert_splice_terms(&[
            // The span ends on the exon's last base, flush with the donor bases.
            (Forward, 10_299, 2, 5, &[], &[REGION_TERM]),
            // The span ends two bases short of the donor bases.
            (Forward, 10_297, 2, 7, &[], &[REGION_TERM]),
            // The span is intronic bases 4 and 3, in the donor region, flush
            // with the donor bases.
            (Reverse, 10_397, 2, 4, &[], &[DONOR_REGION_TERM]),
            // The span ends four bases short of the donor bases.
            (Reverse, 10_393, 2, 8, &[], &[REGION_TERM]),
        ]);
    }

    /// The acceptor bases of intron 2 are 10_399-10_400 on the forward strand
    /// and 10_301-10_302 on the reverse strand. An acceptor term that rests
    /// only on positions past the span is not written, whether the span ends
    /// beside the acceptor bases, a few bases before them, or at the other end
    /// of a 24-base intron.
    #[test]
    fn an_acceptor_term_resting_on_bases_past_the_span_is_not_written() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        assert_splice_terms(&[
            // The span ends two bases short of the acceptor bases, in the tract.
            (Forward, 10_394, 3, 7, &[], &[TRACT_TERM, REGION_TERM]),
            // The span ends on the exon's base beside the acceptor bases.
            (Reverse, 10_299, 2, 4, &[], &[REGION_TERM]),
            // The span ends two bases short of the acceptor bases.
            (Reverse, 10_297, 2, 6, &[], &[REGION_TERM]),
        ]);

        // A 24-base intron, 20_301-20_324: the span is its bases 3 to 5 and the
        // alternate allele's 24 bases run on over the acceptor bases
        // 20_323-20_324 and the tract. The span holds the 5th base alone.
        let short_intron = non_coding_transcript(Forward, &[(20_000, 20_300), (20_325, 20_700)]);
        assert_eq!(
            splice_terms(&short_intron, &replacement(20_303, 3, 24, &[])),
            [FIFTH_BASE_TERM]
        );
    }

    /// The 5th base of intron 2 is 10_305 on the forward strand and 10_396 on
    /// the reverse strand. A 5th-base term that rests only on positions past
    /// the span is not written; the donor region is written in its place where
    /// the span holds one of its bases (intronic bases 3 to 6).
    #[test]
    fn a_donor_fifth_base_term_resting_on_bases_past_the_span_is_not_written() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        assert_splice_terms(&[
            // The span is intronic bases 3 and 4, flush with the 5th base.
            (Forward, 10_303, 2, 4, &[], &[DONOR_REGION_TERM]),
            // The span is the exon's last two bases.
            (Forward, 10_299, 2, 8, &[], &[REGION_TERM]),
            // The span is intronic bases 7 and 6, flush with the 5th base.
            (Reverse, 10_394, 2, 4, &[], &[DONOR_REGION_TERM]),
            // The span is intronic bases 8 and 7.
            (Reverse, 10_393, 2, 5, &[], &[REGION_TERM]),
        ]);
    }

    /// The donor region of intron 2 is 10_303-10_306 on the forward strand and
    /// 10_395-10_398 on the reverse strand. A donor-region term that rests only
    /// on positions past the span is not written.
    #[test]
    fn a_donor_region_term_resting_on_bases_past_the_span_is_not_written() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        assert_splice_terms(&[
            // The span is the two donor bases, flush with the donor region.
            (Forward, 10_301, 2, 4, &[], &[DONOR_TERM]),
            // The span is the exon's last two bases.
            (Forward, 10_299, 2, 6, &[], &[REGION_TERM]),
            // The span is intronic bases 9 and 8; the alternate allele reaches
            // the donor region's outermost base.
            (Reverse, 10_392, 2, 4, &[], &[REGION_TERM]),
            // The same span with one base more, reaching the 5th base too.
            (Reverse, 10_392, 2, 5, &[], &[REGION_TERM]),
        ]);
    }

    /// The polypyrimidine tract of intron 2 is 10_384-10_398 on the forward
    /// strand and 10_303-10_317 on the reverse strand. A tract term that rests
    /// only on positions past the span is not written.
    #[test]
    fn a_polypyrimidine_tract_term_resting_on_bases_past_the_span_is_not_written() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        assert_splice_terms(&[
            // The span ends eleven bases short of the tract.
            (Forward, 10_370, 3, 20, &[], &[]),
            // The span is the two acceptor bases, flush with the tract.
            (Reverse, 10_301, 2, 5, &[], &[ACCEPTOR_TERM]),
        ]);
    }

    /// A term stands beside positions past the span when the span itself holds
    /// a base of a window of that term, even at an index where the two alleles
    /// hold the same base, as Ensembl VEP writes it
    /// (BaseTranscriptVariationAllele.pm:176-204).
    #[test]
    fn a_term_whose_window_the_span_holds_stands_beside_bases_past_the_span() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        assert_splice_terms(&[
            // The span holds the first donor base at an index of equal bases;
            // the alternate allele's third base stands on the second.
            (Forward, 10_300, 2, 3, &[1], &[DONOR_TERM]),
            // The span holds the first donor base on the reverse strand the
            // same way, and three bases of the donor region.
            (
                Reverse,
                10_392,
                8,
                10,
                &[2, 3, 4, 7],
                &[DONOR_REGION_TERM, DONOR_TERM],
            ),
            // The span holds the tract's outermost base at an index of equal
            // bases, after the two acceptor bases.
            (Reverse, 10_301, 3, 4, &[1, 2], &[ACCEPTOR_TERM, TRACT_TERM]),
            (Reverse, 10_301, 3, 26, &[2], &[ACCEPTOR_TERM, TRACT_TERM]),
            // The span holds the second acceptor base at an index of equal bases.
            (Reverse, 10_300, 2, 12, &[1], &[ACCEPTOR_TERM]),
        ]);

        // With one base more the alternate allele also reaches the donor
        // region, which the span does not hold: the donor term stands.
        let tx = four_exon_non_coding_transcript(Forward);
        let terms = splice_terms(&tx, &replacement(10_300, 2, 4, &[1]));
        assert!(terms.iter().any(|term| term == DONOR_TERM), "{terms:?}");
    }

    /// A record of 204 bases from intron 2 into intron 3 holds the acceptor
    /// bases of intron 2 (10_399-10_400) at indices of equal bases, and its
    /// alternate allele's bases past the span stand on the acceptor bases of
    /// intron 3 (10_599-10_600). The span holds a base of an acceptor window
    /// of the transcript, so the acceptor term stands.
    #[test]
    fn a_long_record_over_two_introns_keeps_its_acceptor_term() {
        use vep_core::coordinate::Strand::Forward;
        let tx = four_exon_non_coding_transcript(Forward);
        let same: Vec<usize> = (1..203).filter(|&i| i != 100).collect();
        assert_eq!(
            splice_terms(&tx, &replacement(10_390, 204, 300, &same)),
            [ACCEPTOR_TERM]
        );
    }

    /// An insertion replaces no base and is read by its two neighbours
    /// (VariationEffect.pm:97-100, BaseTranscriptVariationAllele.pm:152-159):
    /// with one neighbour inside a window it keeps the `splice_region` and
    /// polypyrimidine tract terms Ensembl VEP writes.
    #[test]
    fn an_insertion_with_one_neighbour_in_a_window_keeps_its_splice_terms() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        let forward = four_exon_non_coding_transcript(Forward);
        let reverse = four_exon_non_coding_transcript(Reverse);
        // Between the intron's second and third base at its lower end.
        assert_eq!(
            splice_terms(&forward, &insertion_after(10_302)),
            [REGION_TERM]
        );
        assert_eq!(
            splice_terms(&reverse, &insertion_after(10_302)),
            [TRACT_TERM, REGION_TERM]
        );
        // Between the intron's last base and the exon after it.
        assert_eq!(
            splice_terms(&forward, &insertion_after(10_400)),
            [REGION_TERM]
        );
        assert_eq!(
            splice_terms(&reverse, &insertion_after(10_400)),
            [REGION_TERM]
        );
        // Between the tract's outermost base and the base outside it.
        assert_eq!(
            splice_terms(&forward, &insertion_after(10_383)),
            [TRACT_TERM]
        );
        assert_eq!(
            splice_terms(&reverse, &insertion_after(10_317)),
            [TRACT_TERM]
        );
        // Between a two-base intron's upper base and the exon after it.
        let exons = [(30_000, 30_300), (30_303, 30_700)];
        for strand in [Forward, Reverse] {
            let tx = non_coding_transcript(strand, &exons);
            assert_eq!(
                splice_terms(&tx, &insertion_after(30_302)),
                [REGION_TERM],
                "{strand:?}"
            );
        }
    }

    /// A record whose alternate allele is as long as its reference allele, or
    /// shorter, places no index past its span and keeps its terms.
    #[test]
    fn a_record_without_a_longer_alternate_allele_keeps_its_splice_terms() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        assert_splice_terms(&[
            (Forward, 10_299, 3, 3, &[], &[DONOR_TERM]),
            (Forward, 10_303, 2, 2, &[], &[DONOR_REGION_TERM]),
            (Reverse, 10_300, 3, 3, &[], &[ACCEPTOR_TERM]),
            (Forward, 10_299, 5, 2, &[], &[DONOR_REGION_TERM, DONOR_TERM]),
            (Reverse, 10_302, 6, 2, &[], &[ACCEPTOR_TERM, TRACT_TERM]),
        ]);
    }

    /// The intron term is read from the same regions as the splice terms
    /// (BaseTranscriptVariationAllele.pm:143-149, on intronic bases 3 to the
    /// third-last), so an alternate allele's bases past the span set it too. A
    /// span that lies wholly in the exon beside intron 2 (10_301-10_400), one
    /// to three bases short of it, takes out no intronic base, and the term is
    /// not written. The exon term and the `splice_region` the span's own exonic
    /// bases give stand.
    #[test]
    fn an_intron_term_resting_on_bases_past_a_span_in_the_exon_is_not_written() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        let cases: Vec<SpliceCase> = [Forward, Reverse]
            .into_iter()
            .flat_map(|strand| -> [SpliceCase; 3] {
                [
                    // The span ends on the exon's base beside the intron.
                    (strand, 10_299, 2, 6, &[], &[EXON_TERM, REGION_TERM]),
                    // The span ends one base before it.
                    (strand, 10_298, 2, 6, &[], &[EXON_TERM, REGION_TERM]),
                    // The span ends two bases before it.
                    (strand, 10_297, 2, 7, &[], &[EXON_TERM, REGION_TERM]),
                ]
            })
            .collect();
        assert_whole_rows(&cases);
    }

    /// The intron term stands where the span itself takes out a base of an
    /// intron of the transcript, the intron's two end bases included, as
    /// Ensembl VEP writes it.
    #[test]
    fn an_intron_term_stands_where_the_span_holds_a_base_of_an_intron() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        assert_whole_rows(&[
            // The exon's base beside intron 2 and the intron's first base.
            (
                Forward,
                10_300,
                2,
                5,
                &[],
                &[INTRON_TERM, EXON_TERM, DONOR_TERM],
            ),
            (
                Reverse,
                10_300,
                2,
                5,
                &[],
                &[INTRON_TERM, EXON_TERM, ACCEPTOR_TERM],
            ),
            // The exon's base beside the intron and the intron's two end bases.
            (
                Forward,
                10_300,
                3,
                6,
                &[],
                &[INTRON_TERM, EXON_TERM, DONOR_TERM],
            ),
            (
                Reverse,
                10_300,
                3,
                6,
                &[],
                &[INTRON_TERM, EXON_TERM, ACCEPTOR_TERM],
            ),
            // The intron's two end bases alone.
            (
                Forward,
                10_301,
                2,
                5,
                &[],
                &[INTRON_TERM, TRANSCRIPT_TERM, DONOR_TERM],
            ),
            (
                Reverse,
                10_301,
                2,
                5,
                &[],
                &[INTRON_TERM, TRANSCRIPT_TERM, ACCEPTOR_TERM],
            ),
            // The last base of intron 2 and the whole exon after it: the bases
            // past the span stand on intron 3, and the span holds a base of
            // intron 2.
            (
                Forward,
                10_400,
                101,
                104,
                &[],
                &[INTRON_TERM, EXON_TERM, ACCEPTOR_TERM],
            ),
        ]);
    }

    /// An insertion and a record of equal allele lengths place no index past
    /// their span and keep the intron term.
    #[test]
    fn an_insertion_and_an_equal_length_record_keep_their_intron_term() {
        use vep_core::coordinate::Strand::{Forward, Reverse};
        let config = EffectsConfig::default();
        let forward = four_exon_non_coding_transcript(Forward);
        let reverse = four_exon_non_coding_transcript(Reverse);
        assert_eq!(
            sorted_terms(&forward, &insertion_after(10_302), &config),
            [INTRON_TERM, TRANSCRIPT_TERM, REGION_TERM]
        );
        assert_eq!(
            sorted_terms(&reverse, &insertion_after(10_302), &config),
            [INTRON_TERM, TRANSCRIPT_TERM, TRACT_TERM, REGION_TERM]
        );
        assert_whole_rows(&[
            // Two intronic bases, 3 and 4, replaced by two.
            (
                Forward,
                10_303,
                2,
                2,
                &[],
                &[INTRON_TERM, TRANSCRIPT_TERM, DONOR_REGION_TERM],
            ),
            // Five bases from the exon's end into the intron replaced by five.
            (
                Forward,
                10_299,
                5,
                5,
                &[],
                &[INTRON_TERM, EXON_TERM, DONOR_REGION_TERM, DONOR_TERM],
            ),
        ]);
    }

    /// A reference whose bases 10_107 to 10_110 are `A` inside a run of `C`,
    /// as a one-sequence FASTA with its index.
    fn homopolymer_fasta() -> (tempfile::TempDir, vep_fasta::IndexedFasta) {
        let dir = tempfile::tempdir().expect("tempdir");
        let fa = dir.path().join("ref.fa");
        let mut seq = vec![b'C'; 10_800];
        for base in &mut seq[10_106..10_110] {
            *base = b'A';
        }
        let header = ">21\n";
        let mut contents = header.as_bytes().to_vec();
        contents.extend_from_slice(&seq);
        contents.push(b'\n');
        std::fs::write(&fa, &contents).expect("write FASTA");
        std::fs::write(
            dir.path().join("ref.fa.fai"),
            format!(
                "21\t{}\t{}\t{}\t{}\n",
                seq.len(),
                header.len(),
                seq.len(),
                seq.len() + 1
            ),
        )
        .expect("write FASTA index");
        let fasta = vep_fasta::IndexedFasta::from_path(&fa).expect("load FASTA");
        (dir, fasta)
    }

    /// With the 3' shift on, the splice windows read the shifted span while the
    /// endpoint passes anchor on the input span. A deletion of two `A`s in a
    /// four-`A` run at intronic bases 7 to 8 of intron 1 shifts to bases 9 to
    /// 10 on the forward strand, out of the splice_region window (intronic 3
    /// to 8); on the reverse strand a deletion at bases 9 to 10 from the
    /// acceptor shifts to 7 to 8, into the splice_region window of the
    /// polypyrimidine tract. Without the shift each keeps the terms of its
    /// input span.
    #[test]
    fn shifted_deletion_reads_the_splice_windows_of_its_shifted_span() {
        use vep_core::coordinate::Strand;
        let (_dir, fasta) = homopolymer_fasta();
        let shifted = EffectsConfig {
            reference_fasta: Some(Arc::new(fasta)),
            enable_indel_3prime_shift: true,
            ..Default::default()
        };
        let unshifted = EffectsConfig::default();

        let forward = four_exon_non_coding_transcript(Strand::Forward);
        let at_7_8 = InputVariant::new("21".into(), 10_107, 10_108, b"AA".to_vec(), b"-".to_vec());
        assert_eq!(
            sorted_terms(&forward, &at_7_8, &unshifted),
            [
                "intron_variant",
                "non_coding_transcript_variant",
                "splice_region_variant"
            ]
        );
        assert_eq!(
            sorted_terms(&forward, &at_7_8, &shifted),
            ["intron_variant", "non_coding_transcript_variant"]
        );

        let reverse = four_exon_non_coding_transcript(Strand::Reverse);
        let at_9_10 = InputVariant::new("21".into(), 10_109, 10_110, b"AA".to_vec(), b"-".to_vec());
        assert_eq!(
            sorted_terms(&reverse, &at_9_10, &unshifted),
            [
                "intron_variant",
                "non_coding_transcript_variant",
                "splice_polypyrimidine_tract_variant"
            ]
        );
        assert_eq!(
            sorted_terms(&reverse, &at_9_10, &shifted),
            [
                "intron_variant",
                "non_coding_transcript_variant",
                "splice_polypyrimidine_tract_variant",
                "splice_region_variant"
            ]
        );
    }

    /// A single-exon coding model whose alignment deleted three genomic
    /// stretches, as a RefSeq model's deletion edits leave them in the cache:
    /// two bases inside the 5' UTR (between cDNA 10 and 11), fifteen inside the
    /// coding sequence (between cDNA 60 and 61, codons 10 and 11) and two inside
    /// the 3' UTR (between cDNA 150 and 151). The one exon spans the whole model,
    /// 188 genomic bases for 169 of cDNA (a 30-base 5' UTR, a 99-base coding
    /// sequence, a 40-base 3' UTR), and the mapper pairs the deleted bases with
    /// nothing. On the reverse strand the model's first base is the span's last.
    fn edited_single_exon_model(strand: Strand) -> Transcript {
        use vep_core::transcript::{BamEditStatus, Exon, ExonCoordMapper, MapperPair};
        const START: u64 = 40_000_001;
        const END: u64 = 40_000_188;
        let cds = format!("ATG{}TAA", "GCT".repeat(31));
        let mut tx =
            crate::test_helpers::make_single_exon_transcript(strand, &cds, &"A".repeat(40));
        let pair = |from: (u64, u64), to: (u64, u64)| MapperPair {
            from_start: from.0,
            from_end: from.1,
            to_start: to.0,
            to_end: to.1,
            ori: strand.as_i8(),
        };
        let (pairs, coding) = match strand {
            Strand::Forward => (
                vec![
                    pair((1, 10), (START, 40_000_010)),
                    pair((11, 60), (40_000_013, 40_000_062)),
                    pair((61, 150), (40_000_078, 40_000_167)),
                    pair((151, 169), (40_000_170, END)),
                ],
                (40_000_033, 40_000_146),
            ),
            Strand::Reverse => (
                vec![
                    pair((151, 169), (START, 40_000_019)),
                    pair((61, 150), (40_000_022, 40_000_111)),
                    pair((11, 60), (40_000_127, 40_000_176)),
                    pair((1, 10), (40_000_179, END)),
                ],
                (40_000_043, 40_000_156),
            ),
        };
        tx.end = END;
        tx.exons = vec![Exon {
            start: START,
            end: END,
            rank: 1,
            phase: -1,
            end_phase: -1,
        }];
        tx.coding_region_start = Some(coding.0);
        tx.coding_region_end = Some(coding.1);
        tx.translation_start = Some(coding.0);
        tx.translation_end = Some(coding.1);
        tx.bam_edit_status = Some(BamEditStatus::Ok);
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.sorted_exons = tx.exons.clone();
        vefc.mapper.as_mut().unwrap().exon_coord_mapper = ExonCoordMapper::new(pairs);
        tx
    }

    /// The deleted stretches of [`edited_single_exon_model`] as genomic spans:
    /// the one in the 5' UTR, the one in the coding sequence, the one in the
    /// 3' UTR.
    fn deleted_spans(strand: Strand) -> [(u64, u64); 3] {
        match strand {
            Strand::Forward => [
                (40_000_011, 40_000_012),
                (40_000_063, 40_000_077),
                (40_000_168, 40_000_169),
            ],
            Strand::Reverse => [
                (40_000_177, 40_000_178),
                (40_000_112, 40_000_126),
                (40_000_020, 40_000_021),
            ],
        }
    }

    /// A substitution at the first base of `span`, an insertion between its
    /// first two bases and a deletion of all of it.
    fn alleles_inside(span: (u64, u64)) -> [(&'static str, InputVariant); 3] {
        let (lo, hi) = span;
        let bases = (hi - lo + 1) as usize;
        [
            (
                "substitution",
                InputVariant::new("21".into(), lo, lo, b"A".to_vec(), b"G".to_vec()),
            ),
            (
                "insertion",
                InputVariant::new("21".into(), lo + 1, lo, b"-".to_vec(), b"T".to_vec()),
            ),
            (
                "deletion",
                InputVariant::new("21".into(), lo, hi, vec![b'A'; bases], b"-".to_vec()),
            ),
        ]
    }

    /// On a model whose alignment deleted bases, a substitution, an insertion
    /// and a deletion wholly inside the deleted bases take the term of the
    /// transcript region whose genomic span holds them, with no cDNA, CDS or
    /// protein position, on both strands: the bases carry no transcript
    /// coordinate, so the model gives no closer reading. Ensembl VEP 116.2
    /// writes `intergenic_variant` on such a row, the default it assigns when
    /// no predicate fires on a mapper gap.
    #[test]
    fn alleles_inside_deleted_edit_bases_take_the_region_term() {
        use crate::test_helpers::assert_consequence_set_eq;
        let config = EffectsConfig::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let tx = edited_single_exon_model(strand);
            let [utr5, cds, utr3] = deleted_spans(strand);
            let regions = [
                (utr5, "5_prime_UTR_variant"),
                (cds, "coding_sequence_variant"),
                (utr3, "3_prime_UTR_variant"),
            ];
            for (span, term) in regions {
                for (kind, variant) in alleles_inside(span) {
                    let at = format!("{strand:?} {kind} at {}-{}", span.0, span.1);
                    let tc = calculate_consequences(&variant, &tx, &config)
                        .unwrap_or_else(|| panic!("{at}: no row"));
                    assert_consequence_set_eq(&tc, &[term]);
                    assert_eq!(tc.impact, Impact::MODIFIER, "{at}");
                    let fields = crate::display::display_fields(&variant, &tx);
                    assert_eq!(
                        (
                            fields.cdna_position,
                            fields.cds_position,
                            fields.protein_position
                        ),
                        (None, None, None),
                        "{at}: the deleted bases have no transcript coordinate"
                    );
                }
            }
        }
    }

    /// An allele reaching one mapped base keeps its mapped endpoint's reading,
    /// as Ensembl VEP gives it: a deletion from the last mapped coding base into
    /// the deleted bases is `coding_sequence_variant` with its CDS position open
    /// at the deleted end, and an insertion between the last mapped base and the
    /// first deleted one is read at the mapped codon.
    #[test]
    fn alleles_reaching_a_mapped_base_keep_the_mapped_reading() {
        use crate::test_helpers::assert_consequence_set_eq;
        let config = EffectsConfig::default();
        let tx = edited_single_exon_model(Strand::Forward);
        let into_deleted = InputVariant::new(
            "21".into(),
            40_000_062,
            40_000_070,
            vec![b'A'; 9],
            b"-".to_vec(),
        );
        let tc = calculate_consequences(&into_deleted, &tx, &config).unwrap();
        assert_consequence_set_eq(&tc, &["coding_sequence_variant"]);
        let fields = crate::display::display_fields(&into_deleted, &tx);
        assert_eq!(fields.cdna_position.as_deref(), Some("60-?"));
        assert_eq!(fields.cds_position.as_deref(), Some("30-?"));
        assert_eq!(fields.protein_position.as_deref(), Some("10-?"));

        let after_mapped = InputVariant::new(
            "21".into(),
            40_000_063,
            40_000_062,
            b"-".to_vec(),
            b"T".to_vec(),
        );
        let tc = calculate_consequences(&after_mapped, &tx, &config).unwrap();
        assert_consequence_set_eq(&tc, &["frameshift_variant"]);
    }

    /// On a model with no coding sequence the deleted bases lie in a non-coding
    /// exon, whose term every allele kind inside them takes.
    #[test]
    fn alleles_inside_deleted_edit_bases_of_a_non_coding_model_take_the_exon_term() {
        use crate::test_helpers::assert_consequence_set_eq;
        let config = EffectsConfig::default();
        let mut tx = edited_single_exon_model(Strand::Forward);
        tx.biotype = "lncRNA".into();
        tx.translation = None;
        tx.protein_id = None;
        tx.cdna_coding_start = None;
        tx.cdna_coding_end = None;
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation_start = None;
        tx.translation_end = None;
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.translateable_seq = None;
        vefc.three_prime_utr = None;
        let mapper = vefc.mapper.as_mut().unwrap();
        mapper.cdna_coding_start = 0;
        mapper.cdna_coding_end = 0;
        for span in deleted_spans(Strand::Forward) {
            for (kind, variant) in alleles_inside(span) {
                let tc = calculate_consequences(&variant, &tx, &config)
                    .unwrap_or_else(|| panic!("{kind} at {}-{}: no row", span.0, span.1));
                assert_consequence_set_eq(&tc, &["non_coding_transcript_exon_variant"]);
            }
        }
    }
}
