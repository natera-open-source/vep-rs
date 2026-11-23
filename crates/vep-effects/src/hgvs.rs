// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! HGVS notation generation for variant-transcript consequences.
//!
//! Produces HGVSc (coding DNA) and HGVSp (protein) notation strings following
//! the Human Genome Variation Society nomenclature conventions, matching
//! Ensembl VEP's output.
//!
//! Ensembl VEP source citations name modules of ensembl-variation release/116
//! (`Bio/EnsEMBL/Variation/...`); `OutputFactory.pm`, `BaseRunner.pm` and
//! `Pipeline/DumpVEP/Dumper/Core.pm` are in ensembl-vep release/116.3
//! (`Bio/EnsEMBL/VEP/...`), and `Exon` is in ensembl core release/116.
//!
//! The port describes the transcripts the cache holds, every notation naming
//! the cached identifiers; on a model whose alignment to the genome inserted or
//! deleted bases (an `_rna_edit` attribute), HGVSc is numbered against the
//! edited model: the cached mapper pairs fold the edits, so an exonic position
//! takes its edited coordinate and an intronic distance counts from the edited
//! exon end. Ensembl VEP's `hgvs_transcript` walks the exons in their unedited
//! numbering, since the cache dumper dropped the transcript's edit-aware mapper,
//! and corrects an exonic position on an `NM_` or `XM_` model alone by the
//! misalignment offset (TranscriptVariationAllele.pm 1507-1526, 1594-1624);
//! vep-rs keeps that offset as the `REFSEQ_OFFSET` column and documents the
//! numbering as an intended divergence. Two paths of `hgvs_transcript` and
//! `hgvs_protein` stay out: the second shift along a RefSeq transcript's own
//! edited sequence (the tail of `_return_3prime`, 211-256) and the GenBank
//! protein name of such a transcript (1750-1756). Away from a transcript's last
//! bases and its edited bases Ensembl VEP returns before that shift, copying
//! a shift it already has (182-187). The port shifts an insertion or a deletion
//! on an edited transcript 3' along the genome (`hgvs_shift`), and along the
//! model's own sequence only on or beside bases the model holds and the genome
//! lacks (`place_3prime_on_model`). An LRG reference name carries no version,
//! as there (1500, 1758-1760). An edited transcript's HGVSc spells a
//! shifted insertion's bases as the shift rotates them, the rotation
//! `shift_feature_seqs` (1330-1345) gives `feature_seq` on a coding row under
//! `--shift_3prime 1` (`codon`, 813-820); on every other row `hgvs_transcript`
//! copies the unrotated `feature_seq` over the shifted notation's alternate
//! allele (1513), a defect vep-rs does not reproduce and documents as an
//! intended divergence.
//!
//! Reference: <https://varnomen.hgvs.org/>

use std::borrow::Cow;

use vep_core::coordinate::Strand;
use vep_core::transcript::{BamEditStatus, MapperPair, Transcript};
use vep_core::variant::{InputVariant, VariantClass};

/// Convert a single-letter amino acid code to its three-letter abbreviation.
pub(crate) fn amino_acid_three_letter(one_letter: u8) -> &'static str {
    match one_letter {
        b'A' => "Ala",
        b'R' => "Arg",
        b'N' => "Asn",
        b'D' => "Asp",
        b'C' => "Cys",
        b'E' => "Glu",
        b'Q' => "Gln",
        b'G' => "Gly",
        b'H' => "His",
        b'I' => "Ile",
        b'L' => "Leu",
        b'K' => "Lys",
        b'M' => "Met",
        b'F' => "Phe",
        b'P' => "Pro",
        b'S' => "Ser",
        b'T' => "Thr",
        b'W' => "Trp",
        b'Y' => "Tyr",
        b'V' => "Val",
        b'*' => "Ter",
        b'X' => "Xaa",
        b'U' => "Sec",
        b'B' => "Asx",
        b'Z' => "Glx",
        b'J' => "Xle",
        b'O' => "Pyl",
        _ => "Xaa",
    }
}

/// The 3'-shifted span of an insertion or deletion, Ensembl VEP's `_return_3prime(1)`
/// (TranscriptVariationAllele.pm 109): along the genome in the transcript's 3'
/// direction, whatever `--shift_3prime` says, when a reference FASTA is present
/// and the shift moves the variant. `None` under `--shift_hgvs 0`, which
/// `hgvs_transcript` (1441) and `hgvs_protein` (1711) read as no shift at all.
fn hgvs_shift(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
    shift_hgvs: bool,
) -> Option<(u64, u64)> {
    if !shift_hgvs {
        return None;
    }
    let fasta = reference_fasta?;
    crate::consequences::shift_indel_3prime_coords(variant, transcript, fasta)
}

/// The span both notations of a variant are described at, its HGVSc and its
/// HGVSp: its most 3' placement on the transcript's strand, which the HGVS
/// 3' rule assigns the change to. For an insertion whose run of placements
/// is longer than `perform_shift` reads (`shifted` stands where that roll
/// stops), the placement the run ends at, or the first one wholly past the
/// transcript's last base where the run goes on beyond it, which has no
/// notation as every later one has none; `shifted` for every other variant,
/// and under `--shift_3prime 1`, where the row is annotated at the span
/// `perform_shift` reaches. Ensembl VEP reads both notations at `shifted`
/// (TranscriptVariationAllele.pm 1441-1446, 1711-1720).
fn hgvsc_shift(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
    shifted: Option<(u64, u64)>,
    options: HgvsOptions,
) -> Option<(u64, u64)> {
    if !options.shift_hgvs
        || options.shift_3prime
        || variant.variant_class != VariantClass::Insertion
    {
        return shifted;
    }
    let Some(inserted) = crate::consequences::shift_motif(variant) else {
        return shifted;
    };
    let towards_lower = transcript.strand == Strand::Reverse;
    let bound = crate::consequences::perform_shift_bound(inserted.len(), towards_lower);
    let moved = shifted.map_or(0, |(start, _)| start.abs_diff(variant.start));
    if moved < bound {
        return shifted;
    }
    crate::consequences::shift_insertion_to_its_last_placement(
        variant,
        transcript,
        reference_fasta?,
    )
    .or(shifted)
}

/// The allele as the shift leaves it for the codons and `hgvs_protein`: the span
/// `shifted`, the reference bases read there for a deletion, the inserted bases
/// as [`shift_feature_seq`] rotates them for an insertion. `variant` is the
/// allele on the forward strand and `input_strand` the strand the input record
/// wrote it on. `None` for any other class.
pub(crate) fn shifted_variant(
    variant: &InputVariant,
    input_strand: Strand,
    transcript: &Transcript,
    shifted: (u64, u64),
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<InputVariant> {
    let (start, end) = shifted;
    if !matches!(
        variant.variant_class,
        VariantClass::Deletion | VariantClass::Insertion
    ) {
        return None;
    }
    let ref_allele = normalized_hgvs_ref_allele(variant, Some(shifted), reference_fasta)
        .unwrap_or_else(|| b"-".to_vec());
    let alt_allele = shifted_feature_alt(variant, input_strand, transcript.strand, Some(shifted))
        .unwrap_or_else(|| b"-".to_vec());
    let mut moved = InputVariant::new(variant.chr.clone(), start, end, ref_allele, alt_allele);
    // The alleles as the record writes them, which an allele of unknown
    // bases is held to wherever it is read
    // (`PerlCodingEval::reading_unknown_bases`): the row's own two, or every
    // allele of a record of several alternate alleles.
    moved.original_allele_string =
        Some(variant.original_allele_string.clone().unwrap_or_else(|| {
            let alternates: Vec<&[u8]> = if variant.alt_alleles.len() > 1 {
                variant.alt_alleles.iter().map(Vec::as_slice).collect()
            } else {
                vec![variant.alt_allele()]
            };
            std::iter::once(variant.ref_allele.as_slice())
                .chain(alternates)
                .map(String::from_utf8_lossy)
                .collect::<Vec<_>>()
                .join("/")
        }));
    Some(moved)
}

/// The inserted allele of a shifted insertion as the codons and `hgvs_protein`
/// read it, on the forward strand: [`shift_feature_seq`] of the forward-strand
/// allele that `shift_motif` gives, moved `shifted` from the record's position.
/// `None` for any other class, an unmoved insertion or an allele that shifts
/// nothing.
fn shifted_feature_alt(
    variant: &InputVariant,
    input_strand: Strand,
    transcript_strand: Strand,
    shifted: Option<(u64, u64)>,
) -> Option<Vec<u8>> {
    let (motif, shift) = moved_insertion_motif(variant, shifted)?;
    Some(shift_feature_seq(
        &motif,
        input_strand,
        transcript_strand,
        shift,
    ))
}

/// The forward-strand bases of an insertion the shift moved
/// ([`crate::consequences::shift_motif`], whatever strand the input gave them
/// on) and the bases it moved them by. `None` for any other class, an unmoved
/// insertion or an allele of a record whose class shifts nothing.
fn moved_insertion_motif(
    variant: &InputVariant,
    shifted: Option<(u64, u64)>,
) -> Option<(Vec<u8>, usize)> {
    if variant.variant_class != VariantClass::Insertion {
        return None;
    }
    let shift = shifted?.0.abs_diff(variant.start) as usize;
    if shift == 0 {
        return None;
    }
    let motif = crate::consequences::shift_motif(variant)?;
    Some((motif, shift))
}

/// `shift_feature_seqs` (TranscriptVariationAllele.pm 1330-1345) on an inserted
/// allele, given and returned on the forward strand: the allele as the record
/// writes it, on `input_strand`, rotates by `shift` as read on the transcript's
/// strand ([`crate::coding::rotate_shifted_allele`]), so a record written on
/// either strand gives the same forward-strand bases. `hgvs_protein`
/// (1777-1778) and `codon` (820) read the rotated allele back through
/// `feature_seq`, which complements it onto the transcript's strand only when
/// the two strands differ.
fn shift_feature_seq(
    forward_allele: &[u8],
    input_strand: Strand,
    transcript_strand: Strand,
    shift: usize,
) -> Vec<u8> {
    // Reverse-complementing is its own inverse, so one step moves an allele
    // either way between the forward strand and the input's.
    let between_strands = |allele: &[u8]| match input_strand {
        Strand::Forward => allele.to_ascii_uppercase(),
        Strand::Reverse => crate::coding::reverse_complement(&allele.to_ascii_uppercase()),
    };
    let mut written = between_strands(forward_allele);
    crate::coding::rotate_shifted_allele(
        &mut written,
        shift as u64,
        input_strand,
        transcript_strand,
    );
    between_strands(&written)
}

/// The transcript's genomic span read in transcript orientation (Ensembl VEP's
/// transcript feature slice): position 1 is the transcript's first base on its
/// own strand. Bases come from the reference FASTA, complemented on the reverse
/// strand; without a FASTA nothing can be read.
struct TranscriptSlice<'a> {
    chr: &'a str,
    tr_start: i64,
    tr_end: i64,
    reverse: bool,
    fasta: Option<&'a vep_fasta::IndexedFasta>,
}

impl TranscriptSlice<'_> {
    fn len(&self) -> i64 {
        self.tr_end - self.tr_start + 1
    }

    /// Genomic position of 1-based slice position `pos`.
    fn genomic(&self, pos: i64) -> i64 {
        if self.reverse {
            self.tr_end - pos + 1
        } else {
            self.tr_start + pos - 1
        }
    }

    /// 1-based slice position of genomic position `g`, the inverse of `genomic`.
    fn position(&self, g: i64) -> i64 {
        if self.reverse {
            self.tr_end - g + 1
        } else {
            g - self.tr_start + 1
        }
    }

    /// Ensembl VEP `substr($slice->seq, $start - 1, $len)`, or `None` when the FASTA is
    /// absent or the span leaves the slice.
    fn substr(&self, start: i64, len: usize) -> Option<Vec<u8>> {
        let fasta = self.fasta?;
        if len == 0 {
            return Some(Vec::new());
        }
        let end = start + len as i64 - 1;
        if start < 1 || end > self.len() {
            return None;
        }
        let (g_lo, g_hi) = if self.reverse {
            (self.genomic(end), self.genomic(start))
        } else {
            (self.genomic(start), self.genomic(end))
        };
        let mut seq = fasta.sequence(self.chr, g_lo as u64, g_hi as u64)?;
        seq.make_ascii_uppercase();
        if self.reverse {
            seq = crate::coding::reverse_complement(&seq);
        }
        Some(seq)
    }

    /// The base at slice position `pos` on the transcript's strand, read along
    /// the chromosome past either end of the transcript too; `None` without a
    /// FASTA or off the chromosome.
    fn base_at(&self, pos: i64) -> Option<u8> {
        let g = u64::try_from(self.genomic(pos)).ok()?;
        let base = self.fasta?.base(self.chr, g)?;
        Some(if self.reverse {
            vep_core::codon::complement_base(base)
        } else {
            base
        })
    }
}

/// Ensembl VEP's HGVS variant type (`hgvs_variant_notation`, Utils/Sequence.pm 493).
#[derive(Debug, Clone, PartialEq, Eq)]
enum HgvscKind {
    Del,
    Sub,
    Inv,
    Delins,
    Dup,
    Ins,
    /// `[n]`: the alternate allele is the reference repeated `n` times, `n > 2`.
    Multiple(usize),
}

/// Ensembl VEP's `$hgvs_notation` hash for a transcript-level description.
struct HgvscNotation {
    start: i64,
    end: i64,
    ref_seq: Vec<u8>,
    alt_seq: Vec<u8>,
    kind: HgvscKind,
}

/// Which side of an insertion `hgvs_variant_notation` reads first for a
/// duplication (Utils/Sequence.pm 507-508, 572-573): the bases after the site
/// (`lookup_order` 1, the default) or the bases before it (`-1`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DupLookup {
    AfterFirst,
    BeforeFirst,
}

/// `hgvs_variant_notation` (Utils/Sequence.pm 493): the type and displayed span
/// of `alt_seq` replacing slice positions `ref_start..=ref_end`. `None` where
/// Ensembl VEP returns `undef` (the alleles are equal) or the reference bases cannot be
/// read.
fn hgvs_variant_notation(
    slice: &TranscriptSlice<'_>,
    fallback_ref: &[u8],
    alt_seq: Vec<u8>,
    ref_start: i64,
    ref_end: i64,
    lookup: DupLookup,
) -> Option<HgvscNotation> {
    let ref_length = (ref_end - ref_start + 1).max(0) as usize;
    let ref_seq = match slice.substr(ref_start, ref_length) {
        Some(s) => s,
        None if slice.fasta.is_none() => fallback_ref.to_vec(),
        None => return None,
    };
    if ref_seq == alt_seq {
        return None;
    }
    let alt_length = alt_seq.len();
    let mut n = HgvscNotation {
        start: ref_start,
        end: ref_end,
        ref_seq,
        alt_seq,
        kind: HgvscKind::Del,
    };
    if ref_length == 0 && alt_length > 0 {
        let dup_after = |n: &mut HgvscNotation| {
            let after = slice.substr(ref_end + 1, alt_length);
            if after.as_deref() != Some(n.alt_seq.as_slice()) {
                return false;
            }
            n.end = n.start + alt_length as i64 - 1;
            n.kind = HgvscKind::Dup;
            true
        };
        let dup_before = |n: &mut HgvscNotation| {
            let before = slice.substr(ref_end - alt_length as i64 + 1, alt_length);
            if before.as_deref() != Some(n.alt_seq.as_slice()) {
                return false;
            }
            n.start = n.end - alt_length as i64 + 1;
            n.kind = HgvscKind::Dup;
            true
        };
        let is_dup = match lookup {
            DupLookup::AfterFirst => dup_after(&mut n) || dup_before(&mut n),
            DupLookup::BeforeFirst => dup_before(&mut n) || dup_after(&mut n),
        };
        if is_dup {
            return Some(n);
        }
        (n.start, n.end) = (ref_end, ref_start);
        n.kind = HgvscKind::Ins;
        return Some(n);
    }
    n.kind = replacement_kind(&n.ref_seq, ref_length, &n.alt_seq);
    Some(n)
}

/// The type `hgvs_variant_notation` gives `alt_seq` standing in place of the
/// `ref_length` reference bases `ref_seq` (Utils/Sequence.pm 493, every branch
/// but the insertion's): a deletion for an empty allele; for an allele of the
/// reference's length a substitution of one base, an inversion of the
/// reference's reverse complement, a `delins` otherwise; a duplication or a
/// repeat count for the reference repeated; a `delins` for any other allele.
fn replacement_kind(ref_seq: &[u8], ref_length: usize, alt_seq: &[u8]) -> HgvscKind {
    let alt_length = alt_seq.len();
    if alt_length == 0 {
        return HgvscKind::Del;
    }
    if ref_length == alt_length {
        return if ref_length == 1 {
            HgvscKind::Sub
        } else if alt_seq == crate::coding::reverse_complement(ref_seq) {
            HgvscKind::Inv
        } else {
            HgvscKind::Delins
        };
    }
    if alt_length.is_multiple_of(ref_length) {
        let multiple = alt_length / ref_length;
        if alt_seq == ref_seq.repeat(multiple) {
            return if multiple == 2 {
                HgvscKind::Dup
            } else {
                HgvscKind::Multiple(multiple)
            };
        }
    }
    HgvscKind::Delins
}

/// Bases an edited model holds and the genome lacks, or a reference allele's
/// span on the model that holds such bases, in cDNA positions: the bases or
/// the span, and the exon they lie in, every mapper pair inside their exon
/// of the exon list ([`exon_cdna_bounds`]), a pair past a genome base the
/// model drops included.
/// `own_last` is the last of the bases
/// the model numbers of its own between that exon's last base and the next
/// exon's first, the exon's last base where it numbers none. `intron_after`
/// is the intron between that exon and the next, its first slice position and
/// its length, where the model numbers the next exon's first base directly
/// after the exon's last; a model that holds bases of its own after that base
/// has none.
struct ModelBases {
    first: i64,
    last: i64,
    exon_first: i64,
    exon_last: i64,
    own_last: i64,
    intron_after: Option<(i64, i64)>,
}

/// Where a position of [`place_3prime_on_model`]'s walk lies: on the model's
/// cDNA, or in the intron after the exon, at a slice position.
enum WalkPlace {
    Cdna(i64),
    Intron(i64),
}

impl ModelBases {
    /// The walk counts on from the exon's last base through `intron_after`
    /// and then through the next exon: a position past that base is an intron
    /// base while the intron lasts and a cDNA position after it.
    fn place(&self, pos: i64) -> WalkPlace {
        match self.intron_after {
            Some((_, length)) if pos > self.exon_last + length => WalkPlace::Cdna(pos - length),
            Some((first, _)) if pos > self.exon_last => {
                WalkPlace::Intron(first + pos - self.exon_last - 1)
            }
            _ => WalkPlace::Cdna(pos),
        }
    }

    /// Whether `n`, in the walk's positions, ends on the exon's last base or
    /// past it: an insertion by the base it follows, any other notation by
    /// its last position.
    fn reaches_the_exons_last_base(&self, n: &HgvscNotation) -> bool {
        let last = if n.kind == HgvscKind::Ins {
            n.start
        } else {
            n.end
        };
        last >= self.exon_last
    }

    /// These bases for a walk that ends on the exon's last base.
    fn within_the_exon(&self) -> ModelBases {
        ModelBases {
            own_last: self.exon_last,
            intron_after: None,
            ..*self
        }
    }
}

/// The first and last cDNA positions of the exon of the model's exon list
/// that holds genomic `start..=end`, read from the mapper pairs inside the
/// exon: an edit that adds or drops bases splits the exon's pair around it.
/// `None` where no exon of the list holds both positions.
fn exon_cdna_bounds(
    transcript: &Transcript,
    exons: &[ExonSpan],
    start: i64,
    end: i64,
) -> Option<(i64, i64)> {
    let exon = transcript
        .exons
        .iter()
        .find(|exon| exon.start as i64 <= start && end <= exon.end as i64)?;
    let pairs = || {
        exons
            .iter()
            .filter(|pair| exon.start as i64 <= pair.start && pair.end <= exon.end as i64)
    };
    Some((
        pairs().map(|pair| pair.cdna_start).min()?,
        pairs().map(|pair| pair.cdna_end).max()?,
    ))
}

/// The exon holding slice position `pos`, as its index in `exons`, and the
/// position's cDNA number; `None` for a position in no exon.
fn exonic_cdna_position(
    slice: &TranscriptSlice<'_>,
    exons: &[ExonSpan],
    pos: i64,
) -> Option<(usize, i64)> {
    let g = slice.genomic(pos);
    let index = exons.iter().position(|e| e.start <= g && g <= e.end)?;
    let exon = &exons[index];
    let cdna = if slice.reverse {
        exon.cdna_start + (exon.end - g)
    } else {
        exon.cdna_start + (g - exon.start)
    };
    Some((index, cdna))
}

/// The model's bases between slice position `pos` and the position after it,
/// when both are exonic and the model numbers them more than one apart: the
/// mapper pairs fold such bases in between the two. Their exon is the exon
/// of the exon list that holds both positions, so a base the 3' rule places
/// may lie past a genome base the model drops; where the list holds the two
/// in no one exon, it is every mapper pair the genome joins end to end with
/// theirs.
fn model_bases_between(
    transcript: &Transcript,
    slice: &TranscriptSlice<'_>,
    exons: &[ExonSpan],
    pos: i64,
) -> Option<ModelBases> {
    let ((index_before, before), (index_after, after)) = (
        exonic_cdna_position(slice, exons, pos)?,
        exonic_cdna_position(slice, exons, pos + 1)?,
    );
    if after <= before + 1 {
        return None;
    }
    let mut lowest = index_before.min(index_after);
    let mut highest = index_before.max(index_after);
    while lowest > 0 && exons[lowest - 1].end + 1 == exons[lowest].start {
        lowest -= 1;
    }
    while highest + 1 < exons.len() && exons[highest].end + 1 == exons[highest + 1].start {
        highest += 1;
    }
    let joined_pairs = || {
        let joined = &exons[lowest..=highest];
        Some((
            joined.iter().map(|e| e.cdna_start).min()?,
            joined.iter().map(|e| e.cdna_end).max()?,
        ))
    };
    let (start, end) = (slice.genomic(pos), slice.genomic(pos + 1));
    let (exon_first, exon_last) =
        exon_cdna_bounds(transcript, exons, start.min(end), start.max(end))
            .or_else(joined_pairs)?;
    // The next exon on the transcript, the genomic base after this exon's
    // last on the transcript's strand, and the bases between the two exons.
    let next_exon = if slice.reverse {
        lowest.checked_sub(1).map(|next| {
            let next = &exons[next];
            let first = exons[lowest].start - 1;
            (next, first, first - next.end)
        })
    } else {
        exons.get(highest + 1).map(|next| {
            let first = exons[highest].end + 1;
            (next, first, next.start - first)
        })
    };
    let own_last = next_exon.map_or(exon_last, |(next, _, _)| exon_last.max(next.cdna_start - 1));
    let intron_after = next_exon
        .filter(|(_, _, length)| own_last == exon_last && *length > 0)
        .map(|(_, first, length)| (slice.position(first), length));
    Some(ModelBases {
        first: before + 1,
        last: after - 1,
        exon_first,
        exon_last,
        own_last,
        intron_after,
    })
}

/// The model's bases between two neighbours among slice positions
/// `first..=last`. The model holds bases the genome lacks only where two
/// mapper pairs meet end to end on the genome.
fn model_bases_within(
    transcript: &Transcript,
    slice: &TranscriptSlice<'_>,
    exons: &[ExonSpan],
    first: i64,
    last: i64,
) -> Option<ModelBases> {
    exons
        .windows(2)
        .filter(|pair| pair[0].end + 1 == pair[1].start)
        .find_map(|pair| {
            let before = slice.position(if slice.reverse {
                pair[1].start
            } else {
                pair[0].end
            });
            (first..last)
                .contains(&before)
                .then(|| model_bases_between(transcript, slice, exons, before))?
        })
}

/// The span on the model of a reference allele at slice positions
/// `first..=last` whose bases lie on both sides of bases the model holds and
/// the genome lacks ([`model_bases_within`]): from the allele's first base to
/// its last, both held by the model in one exon of its exon list, the span
/// the mapper gives the row as its cDNA and CDS positions. The model's own
/// bases lie in it, and a genome base the model drops between the two ends
/// does not. `None` for an allele that spans no such bases, or whose first
/// or last base the model drops or holds in another exon.
fn model_span_across_model_bases(
    transcript: &Transcript,
    slice: &TranscriptSlice<'_>,
    exons: &[ExonSpan],
    first: i64,
    last: i64,
) -> Option<ModelBases> {
    model_bases_within(transcript, slice, exons, first, last)?;
    let (_, span_first) = exonic_cdna_position(slice, exons, first)?;
    let (_, span_last) = exonic_cdna_position(slice, exons, last)?;
    let (start, end) = (slice.genomic(first), slice.genomic(last));
    let (exon_first, exon_last) =
        exon_cdna_bounds(transcript, exons, start.min(end), start.max(end))?;
    Some(ModelBases {
        first: span_first,
        last: span_last,
        exon_first,
        exon_last,
        own_last: exon_last,
        intron_after: None,
    })
}

/// The notation of `alt_seq` standing in place of the model's `bases`, typed
/// and clipped as an allele replacing reference bases is
/// ([`replacement_kind`], [`hgvsc_clip_alleles`]) and then, while the 3'
/// shift is on (`shift_hgvs`), placed by the 3' rule on the model's own
/// sequence ([`place_3prime_on_model`]), its span in the walk's positions
/// ([`ModelBases::place`]), cDNA positions up to the exon's last base. Under
/// `--shift_hgvs 0` what the clip leaves stays where it is, as every other
/// notation does. `None` when the allele is those bases, which changes
/// nothing on the model, or when the cache holds no cDNA to read them from.
fn replaced_model_bases_notation(
    transcript: &Transcript,
    slice: &TranscriptSlice<'_>,
    bases: &ModelBases,
    alt_seq: Vec<u8>,
    shift_hgvs: bool,
) -> Option<HgvscNotation> {
    let cdna = transcript.spliced_seq()?.as_bytes();
    let ref_seq = cdna
        .get(usize::try_from(bases.first - 1).ok()?..usize::try_from(bases.last).ok()?)?
        .to_ascii_uppercase();
    if ref_seq == alt_seq {
        return None;
    }
    let kind = replacement_kind(&ref_seq, ref_seq.len(), &alt_seq);
    let mut n = HgvscNotation {
        start: bases.first,
        end: bases.last,
        ref_seq,
        alt_seq,
        kind,
    };
    if n.kind != HgvscKind::Dup {
        hgvsc_clip_alleles(&mut n);
    }
    if shift_hgvs {
        place_3prime_on_model(cdna, slice, bases, &mut n);
    }
    if n.start > n.end {
        (n.start, n.end) = (n.end, n.start);
    }
    Some(n)
}

/// Intron bases [`WalkedBases`] reads from the FASTA at a time.
const INTRON_STRETCH: i64 = 64;

/// The bases [`place_3prime_on_model`] walks over: the model's cDNA to the
/// last base it numbers before the next exon ([`ModelBases::own_last`]),
/// then the genome's bases of the intron after the exon
/// ([`ModelBases::intron_after`]), read a stretch at a time.
struct WalkedBases<'a> {
    cdna: &'a [u8],
    slice: &'a TranscriptSlice<'a>,
    bases: &'a ModelBases,
    intron: Vec<u8>,
}

impl WalkedBases<'_> {
    /// The base at walk position `pos`, `None` past the last base the walk
    /// reads or where the FASTA gives no base.
    fn base(&mut self, pos: i64) -> Option<u8> {
        match self.bases.place(pos) {
            WalkPlace::Cdna(cdna) if cdna <= self.bases.own_last => usize::try_from(cdna - 1)
                .ok()
                .and_then(|index| self.cdna.get(index))
                .map(u8::to_ascii_uppercase),
            WalkPlace::Cdna(_) => None,
            WalkPlace::Intron(slice_pos) => {
                let (first, length) = self.bases.intron_after?;
                let index = usize::try_from(slice_pos - first).ok()?;
                while self.intron.len() <= index {
                    let read = self.intron.len() as i64;
                    let stretch = usize::try_from((length - read).min(INTRON_STRETCH)).ok()?;
                    self.intron
                        .extend(self.slice.substr(first + read, stretch)?);
                }
                self.intron.get(index).copied()
            }
        }
    }
}

/// HGVS's 3' rule for what the clip leaves of a replacement of model bases,
/// read on the model's cDNA within the exon: a deletion or a duplication
/// moves on while the base after it is its first base; an insertion moves on
/// while the base after it is the allele's first base, the allele turning
/// with each step, and is then a duplication when it copies the bases before
/// it. An insertion alone goes on past the exon's last base, on the model's
/// own sequence: through the bases the model numbers of its own directly
/// after that base and, where it numbers the next exon's first base directly
/// after the exon's last, onto the genome's intron. HGVS excepts from the
/// rule a change around an exon/exon junction, never one around an
/// exon/intron junction. It ends at the latest on the last base the model
/// numbers of its own after the exon or on the intron's last base, and on
/// the exon's last base where no exon follows. The genome cannot place the
/// model's bases: they are not in it.
fn place_3prime_on_model(
    cdna: &[u8],
    slice: &TranscriptSlice<'_>,
    bases: &ModelBases,
    n: &mut HgvscNotation,
) {
    let mut walked = WalkedBases {
        cdna,
        slice,
        bases,
        intron: Vec::new(),
    };
    match n.kind {
        HgvscKind::Del | HgvscKind::Dup => {
            while n.end < bases.exon_last
                && walked
                    .base(n.end + 1)
                    .is_some_and(|after| Some(after) == walked.base(n.start))
            {
                n.start += 1;
                n.end += 1;
                n.ref_seq.rotate_left(1);
                if n.kind == HgvscKind::Dup {
                    n.alt_seq.rotate_left(1);
                }
            }
        }
        HgvscKind::Ins => {
            // The clip leaves an insertion with `start` one past `end`, the
            // base it follows.
            let mut follows = n.end;
            while walked
                .base(follows + 1)
                .is_some_and(|after| Some(&after) == n.alt_seq.first())
            {
                n.alt_seq.rotate_left(1);
                follows += 1;
            }
            let copied_from = follows - n.alt_seq.len() as i64 + 1;
            let copies_the_bases_before = copied_from >= bases.exon_first
                && (copied_from..=follows)
                    .map(|pos| walked.base(pos))
                    .eq(n.alt_seq.iter().map(|b| Some(*b)));
            if copies_the_bases_before {
                n.kind = HgvscKind::Dup;
                (n.start, n.end) = (copied_from, follows);
            } else {
                (n.start, n.end) = (follows + 1, follows);
            }
        }
        _ => {}
    }
}

/// The notation of `alt_seq` inserted after cDNA position `point`, beside the
/// model's `bases`: an `ins` between `point` and the position after it, or
/// the duplication the model's cDNA shows there. While the 3' shift is on
/// (`shift_hgvs`) the 3' rule places it within the exon
/// ([`place_3prime_on_model`], [`ModelBases::within_the_exon`]); under
/// `--shift_hgvs 0` it stays at `point`, a duplication when it repeats the
/// bases on the side `lookup` reads first, then the other, within the exon.
/// A model with no cached cDNA shows no duplication.
fn inserted_beside_model_bases_notation(
    transcript: &Transcript,
    slice: &TranscriptSlice<'_>,
    bases: &ModelBases,
    point: i64,
    alt_seq: Vec<u8>,
    shift_hgvs: bool,
    lookup: DupLookup,
) -> HgvscNotation {
    let cdna = transcript.spliced_seq().map(str::as_bytes);
    inserted_after_notation(cdna, slice, bases, point, alt_seq, shift_hgvs, lookup)
}

/// The notation of `alt_seq` inserted where the model's `bases` stand, on a
/// model the genome is the reference of (`coding::genome_is_reference`): the
/// bases are no bases of the transcript, so the point is the one after the
/// base 5' of them, and the 3' rule and the duplication test read the
/// model's cDNA without them, as [`inserted_beside_model_bases_notation`]
/// reads it with them. An insertion left at the point is an `ins` after the
/// base 5' of it; a position past the point keeps the number the mapper
/// gives its genomic base. Beside the notation stand the bases the 3' rule
/// moves it from the point, each a base of the genome, 0 under
/// `--shift_hgvs 0`.
fn inserted_at_model_bases_on_the_genome(
    transcript: &Transcript,
    slice: &TranscriptSlice<'_>,
    bases: &ModelBases,
    alt_seq: Vec<u8>,
    shift_hgvs: bool,
    lookup: DupLookup,
) -> (HgvscNotation, i64) {
    let point = bases.first - 1;
    let held = bases.last - point;
    let genome = transcript.spliced_seq().and_then(|cdna| {
        let cdna = cdna.as_bytes();
        let before = cdna.get(..usize::try_from(point).ok()?)?;
        let after = cdna.get(usize::try_from(bases.last).ok()?..)?;
        Some([before, after].concat())
    });
    let exon = ModelBases {
        first: bases.first,
        last: point,
        exon_first: bases.exon_first,
        exon_last: bases.exon_last - held,
        own_last: bases.exon_last - held,
        intron_after: None,
    };
    let mut n = inserted_after_notation(
        genome.as_deref(),
        slice,
        &exon,
        point,
        alt_seq,
        shift_hgvs,
        lookup,
    );
    // The base the placed bases follow: the last one a duplication copies,
    // the first of the two an `ins` stands between.
    let follows = if n.kind == HgvscKind::Dup {
        n.end
    } else {
        n.start
    };
    let moved = if shift_hgvs { follows - point } else { 0 };
    let left_at_the_point = n.kind == HgvscKind::Ins && n.start == point;
    if !left_at_the_point {
        for position in [&mut n.start, &mut n.end] {
            if *position > point {
                *position += held;
            }
        }
    }
    (n, moved)
}

/// The notation of `alt_seq` inserted after position `point` of `cdna`, the
/// exon's first and last positions read from `bases`: an `ins`, or the
/// duplication `cdna` shows there, placed as
/// [`inserted_beside_model_bases_notation`] states.
fn inserted_after_notation(
    cdna: Option<&[u8]>,
    slice: &TranscriptSlice<'_>,
    bases: &ModelBases,
    point: i64,
    alt_seq: Vec<u8>,
    shift_hgvs: bool,
    lookup: DupLookup,
) -> HgvscNotation {
    let length = alt_seq.len() as i64;
    let mut n = HgvscNotation {
        start: point + 1,
        end: point,
        ref_seq: Vec::new(),
        alt_seq,
        kind: HgvscKind::Ins,
    };
    if let Some(cdna) = cdna {
        if shift_hgvs {
            place_3prime_on_model(cdna, slice, &bases.within_the_exon(), &mut n);
        } else {
            let repeats = |first: i64, last: i64| {
                first >= bases.exon_first
                    && last <= bases.exon_last
                    && usize::try_from(first - 1)
                        .ok()
                        .and_then(|from| cdna.get(from..last as usize))
                        .is_some_and(|model| model.eq_ignore_ascii_case(&n.alt_seq))
            };
            let after = (point + 1, point + length);
            let before = (point - length + 1, point);
            let sides = match lookup {
                DupLookup::AfterFirst => [after, before],
                DupLookup::BeforeFirst => [before, after],
            };
            if let Some(&(first, last)) = sides.iter().find(|(first, last)| repeats(*first, *last))
            {
                n.kind = HgvscKind::Dup;
                (n.start, n.end) = (first, last);
            }
        }
    }
    if n.start > n.end {
        (n.start, n.end) = (n.end, n.start);
    }
    n
}

/// A duplication typed on the genomic slice at an insertion point, read on
/// the model when its copied bases lie on both sides of bases the model
/// holds and the genome lacks. Numbered against the edited model from the
/// first copied base to the last, the range would name the model's own bases
/// between them as well, more positions than the allele has bases; Ensembl
/// VEP numbers both ends from the unedited exons, where the copied bases are
/// neighbours (TranscriptVariationAllele.pm 1524-1525). On the model a
/// duplication copies the bases beside the point, as many as the allele:
/// where the cached cDNA shows them equal to the allele, within the exon, `n`
/// becomes that range in cDNA positions and the result is `Some(true)`. Where
/// it does not, the other side of the point is read, if `lookup` has it still
/// to read (Utils/Sequence.pm 577-585), on the genome unless its bases are
/// parted the same way, and an allele that copies neither side is an `ins` at
/// the point. `n` is then left in slice positions, `Some(false)`, as it is,
/// untouched, when the model holds the copied bases side by side. `None` when
/// the cache holds no cDNA to read the model's bases from.
///
/// `before_the_point` is the slice position the allele is inserted after and
/// `typed_alt` the allele as it was typed, where `n` holds another spelling.
fn duplication_across_model_bases(
    transcript: &Transcript,
    slice: &TranscriptSlice<'_>,
    exons: &[ExonSpan],
    before_the_point: i64,
    lookup: DupLookup,
    typed_alt: Option<&[u8]>,
    n: &mut HgvscNotation,
) -> Option<bool> {
    let parted_by =
        |first: i64, last: i64| model_bases_within(transcript, slice, exons, first, last);
    if parted_by(n.start, n.end).is_none() {
        return Some(false);
    }
    let cdna = transcript.spliced_seq()?.as_bytes();
    let inserted = typed_alt.map_or_else(|| n.alt_seq.clone(), <[u8]>::to_vec);
    let length = inserted.len() as i64;
    let typed_before = n.end == before_the_point;
    // The side the genome typed, then the other one where the lookup had not
    // reached it.
    let sides: &[bool] = match (lookup, typed_before) {
        (DupLookup::AfterFirst, false) => &[false, true],
        (DupLookup::BeforeFirst, true) => &[true, false],
        (_, typed) => &[typed],
    };
    for (index, &before) in sides.iter().enumerate() {
        let (first, last) = if before {
            (before_the_point - length + 1, before_the_point)
        } else {
            (before_the_point + 1, before_the_point + length)
        };
        let copied_on_the_genome = index == 0
            || slice.substr(first, inserted.len()).as_deref() == Some(inserted.as_slice());
        if !copied_on_the_genome {
            continue;
        }
        let Some(bases) = parted_by(first, last) else {
            (n.start, n.end) = (first, last);
            return Some(false);
        };
        let Some((_, beside)) =
            exonic_cdna_position(slice, exons, if before { last } else { first })
        else {
            continue;
        };
        let (model_first, model_last) = if before {
            (beside - length + 1, beside)
        } else {
            (beside, beside + length - 1)
        };
        let copied_on_the_model = model_first >= bases.exon_first
            && model_last <= bases.exon_last
            && usize::try_from(model_first - 1)
                .ok()
                .and_then(|from| cdna.get(from..model_last as usize))
                .is_some_and(|model| model.eq_ignore_ascii_case(&inserted));
        if copied_on_the_model {
            (n.start, n.end) = (model_first, model_last);
            return Some(true);
        }
    }
    (n.start, n.end) = (before_the_point, before_the_point + 1);
    n.kind = HgvscKind::Ins;
    Some(false)
}

/// An exonic cDNA position as `_get_cDNA_position` prints it
/// (TranscriptVariationAllele.pm 2755): counted from the start codon on a
/// coding transcript, with no position 0, and from the stop codon behind a `*`
/// past it.
fn exonic_hgvs_coordinate(coding: Option<(i64, i64)>, cdna: i64) -> String {
    match coding {
        Some((_, stop_codon)) if cdna > stop_codon => format!("*{}", cdna - stop_codon),
        Some((start_codon, _)) => (cdna + i64::from(cdna >= start_codon) - start_codon).to_string(),
        None => cdna.to_string(),
    }
}

/// `_clip_alleles` (TranscriptVariationAllele.pm 2188) as `hgvs_transcript`
/// calls it, before `numbering` is set: the shared residues are trimmed from the
/// front and the back, moving `start` and `end`, and the type is re-read from
/// what remains; the `=` and `dup` re-sets, which need `numbering`, never fire.
fn hgvsc_clip_alleles(n: &mut HgvscNotation) {
    let mut ref_seq: &[u8] = &n.ref_seq;
    let mut alt_seq: &[u8] = &n.alt_seq;
    while let (Some(r), Some(a)) = (ref_seq.first(), alt_seq.first()) {
        if r != a {
            break;
        }
        n.start += 1;
        ref_seq = &ref_seq[1..];
        alt_seq = &alt_seq[1..];
    }
    while let (Some(r), Some(a)) = (ref_seq.last(), alt_seq.last()) {
        if r != a {
            break;
        }
        ref_seq = &ref_seq[..ref_seq.len() - 1];
        alt_seq = &alt_seq[..alt_seq.len() - 1];
        n.end -= 1;
    }
    let kind = if ref_seq != b"-" && ref_seq.len() == 1 && alt_seq.len() == 1 && ref_seq != alt_seq
    {
        Some(HgvscKind::Sub)
    } else if ref_seq.is_empty() && !alt_seq.is_empty() {
        Some(HgvscKind::Ins)
    } else if !ref_seq.is_empty() && alt_seq.is_empty() {
        Some(HgvscKind::Del)
    } else {
        None
    };
    n.ref_seq = ref_seq.to_vec();
    n.alt_seq = alt_seq.to_vec();
    if let Some(kind) = kind {
        n.kind = kind;
    }
}

/// Drops from a pair of alleles respelled relative to the transcript the
/// `front` leading and `back` trailing bases the notation's span lost to
/// `_clip_alleles` (TranscriptVariationAllele.pm 2188), where the pair
/// shares them at both ends; the bases the clamp onto the transcript cut are
/// gone before it. A pair that does not share them stays whole, as
/// `hgvs_transcript` writes it (1510-1513), and a base the pair shares
/// inside the span stays too: the notation's positions still name it.
fn clip_respelled_alleles(ref_seq: &mut Vec<u8>, alt_seq: &mut Vec<u8>, front: usize, back: usize) {
    let lost = front + back;
    if lost == 0 || ref_seq.len() < lost || alt_seq.len() < lost {
        return;
    }
    let shared = ref_seq[..front] == alt_seq[..front]
        && ref_seq[ref_seq.len() - back..] == alt_seq[alt_seq.len() - back..];
    if shared {
        for allele in [ref_seq, alt_seq] {
            allele.truncate(allele.len() - back);
            allele.drain(..front);
        }
    }
}

/// Drops from the pair of alleles of a substitution, respelled relative to
/// the transcript, the `front` leading and `back` trailing bases the clip
/// takes from the notation's span, each end from the outside in for as long
/// as the pair shares the base: a base is clipped only where the transcript
/// itself has it. Returns how many of the `front` and of the `back` bases
/// the pair keeps, from the first one the transcript spells differently.
fn clip_substitution_alleles(
    ref_seq: &mut Vec<u8>,
    alt_seq: &mut Vec<u8>,
    front: usize,
    back: usize,
) -> (usize, usize) {
    if ref_seq.len() < front + back || alt_seq.len() < front + back {
        return (front, back);
    }
    let dropped_front = ref_seq
        .iter()
        .zip(alt_seq.iter())
        .take(front)
        .take_while(|(r, a)| r == a)
        .count();
    let dropped_back = ref_seq
        .iter()
        .rev()
        .zip(alt_seq.iter().rev())
        .take(back)
        .take_while(|(r, a)| r == a)
        .count();
    for allele in [ref_seq, alt_seq] {
        allele.truncate(allele.len() - dropped_back);
        allele.drain(..dropped_front);
    }
    (front - dropped_front, back - dropped_back)
}

/// The span on the model, in cDNA positions, of a pair of alleles respelled
/// relative to the transcript whose reference holds bases the model has and
/// the genome lacks beside the span the clip leaves, with the leading and
/// trailing bases the pair loses to it. The clip took `lost.0` bases before
/// slice positions `clipped` and `lost.1` after them, of the record's span
/// `typed`; on the model it stops at the first of the model's bases from
/// either end, and the span runs from the first base the pair keeps to the
/// last, the model's bases between a base the clip takes and the bases it
/// leaves inside it. `None` when the reference holds no such bases, when
/// they all lie inside the span the clip leaves, which then names them
/// already, when the pair does not share the bases the clip takes, or when
/// an end of the record's span is not exonic.
#[allow(clippy::too_many_arguments)]
fn span_beside_model_bases(
    transcript: &Transcript,
    slice: &TranscriptSlice<'_>,
    exons: &[ExonSpan],
    typed: (i64, i64),
    clipped: (i64, i64),
    lost: (usize, usize),
    ref_seq: &[u8],
    alt_seq: &[u8],
) -> Option<((i64, i64), (usize, usize))> {
    if ref_seq.len() as i64 <= typed.1 - typed.0 + 1 {
        return None;
    }
    let joined = |pos: i64| model_bases_between(transcript, slice, exons, pos).is_none();
    let front = (0..lost.0 as i64)
        .take_while(|i| *i == 0 || joined(typed.0 + i - 1))
        .count();
    let back = (0..lost.1 as i64)
        .take_while(|i| *i == 0 || joined(typed.1 - i))
        .count();
    let kept = ref_seq.len().checked_sub(front + back)?;
    let shared = alt_seq.len() >= front + back
        && ref_seq[..front] == alt_seq[..front]
        && ref_seq[ref_seq.len() - back..] == alt_seq[alt_seq.len() - back..];
    let cdna = |pos: i64| exonic_cdna_position(slice, exons, pos).map(|(_, cdna)| cdna);
    let span = (cdna(typed.0)? + front as i64, cdna(typed.1)? - back as i64);
    let names_them_already = cdna(clipped.0).zip(cdna(clipped.1)) == Some(span);
    (shared && !names_them_already && span.1 - span.0 + 1 == kept as i64)
        .then_some((span, (front, back)))
}

/// The alternate allele as read inside the transcript, without the bases the
/// clamp cut from the reference allele, both on the transcript's strand:
/// `front_cut` bases before the first base and `back_cut` after the last. A
/// pair of equal length pairs its bases by position, so the alternate's outside
/// bases go whatever they are. A pair of unequal length has no pairing: its
/// change is aligned 3' within the span, the shared leading bases first and the
/// shared trailing bases after them, as `_clip_alleles`
/// (TranscriptVariationAllele.pm 2204-2242) and the trimming of a minimised
/// record align it. A change lying wholly outside the transcript, a base
/// deleted from or inserted into a run across the transcript's end, changes no
/// base of the transcript, and the reference's bases inside are returned.
/// Otherwise the change's alternate bases replace the bases inside the
/// transcript the change covers, the rest of the bases inside staying the
/// reference's: a change lying wholly inside leaves the alternate without the
/// outside bases, a deletion running over an edge deletes the bases inside it,
/// and a change of other bases straddling an edge has its alternate bases read
/// inside, where Ensembl VEP reads them, without the shared bases outside.
fn alt_cut_as_reference(
    alt: Vec<u8>,
    reference: &[u8],
    front_cut: usize,
    back_cut: usize,
) -> Vec<u8> {
    let cut = front_cut + back_cut;
    if cut == 0 || reference.len() < cut {
        return alt;
    }
    let inside = front_cut..reference.len() - back_cut;
    if alt.len() == reference.len() {
        return alt[inside].to_vec();
    }
    let (lead, tail) = shared_ends_aligned_3prime(&alt, reference);
    let (change_start, change_end) = (lead, reference.len() - tail);
    let past_the_last_base = back_cut > 0 && change_start >= inside.end;
    let before_the_first_base = front_cut > 0 && change_end <= inside.start;
    if past_the_last_base || before_the_first_base {
        return reference[inside].to_vec();
    }
    let kept_before = &reference[inside.start..change_start.max(inside.start)];
    let kept_after = &reference[change_end.min(inside.end)..inside.end];
    let changed = &alt[lead..alt.len() - tail];
    [kept_before, changed, kept_after].concat()
}

/// The bases a pair of alleles of unequal length shares at its two ends, the
/// leading ones counted first and the trailing ones from what is left, which
/// places the change between them most 3' within the reference's span.
fn shared_ends_aligned_3prime(alt: &[u8], reference: &[u8]) -> (usize, usize) {
    let lead = alt
        .iter()
        .zip(reference)
        .take_while(|(a, r)| a == r)
        .count();
    let tail = alt[lead..]
        .iter()
        .rev()
        .zip(reference[lead..].iter().rev())
        .take_while(|(a, r)| a == r)
        .count();
    (lead, tail)
}

/// Whether an allele that keeps bases it shares with its reference, written
/// over a transcript's first or last base, deletes bases that run over that
/// base and can be written nowhere wholly inside the transcript. `reference`
/// holds the span's bases from slice position `raw_start`, of which the clamp
/// cuts `front_cut` before the first base and `back_cut` past the last. The
/// deleted window is aligned 3' within the span, then carried each way along
/// the chromosome while the base it leaves equals the base it takes in; the
/// window nearest the transcript's inside is the most 3' one at the first
/// base and the most 5' one at the last. Where that window still runs over
/// the base, the deletion extends beyond the transcript however it is
/// written, and the HGVS nomenclature gives it to the genomic description
/// alone (DNA deletion, "deletions extending beyond the transcribed
/// region"): the transcript has no description of it. A deletion with a
/// window wholly inside the transcript keeps its row.
///
/// `base_at` reads the base at a slice position of a transcript of `tr_len`
/// bases, along the chromosome past either end of it too, and each walk ends
/// where the answer is settled. The window is carried towards the 3' end only
/// where the clamp cut bases before the first base, until its first base is
/// the transcript's first, and towards the 5' end only where the clamp cut
/// bases past the last base, until its last base is the transcript's last: a
/// window carried that far does not run over that base, and none carried
/// farther does. So each walk makes at most as many moves as the clamp cut
/// bases at that end, two bases read a move, however long the run is.
fn kept_allele_deletes_over_a_transcript_end(
    tr_len: i64,
    base_at: impl Fn(i64) -> Option<u8>,
    reference: &[u8],
    alt: &[u8],
    raw_start: i64,
    front_cut: usize,
    back_cut: usize,
) -> bool {
    if front_cut + back_cut == 0 || alt.is_empty() || alt.len() >= reference.len() {
        return false;
    }
    let (lead, tail) = shared_ends_aligned_3prime(alt, reference);
    if lead + tail != alt.len() {
        return false;
    }
    let raw_end = raw_start + reference.len() as i64 - 1;
    let base = |pos: i64| {
        if (raw_start..=raw_end).contains(&pos) {
            Some(reference[(pos - raw_start) as usize])
        } else {
            base_at(pos)
        }
    };
    let (first, last) = (raw_start + lead as i64, raw_end - tail as i64);
    let to_the_first_base = if front_cut > 0 { 1 - first } else { 0 };
    let mut towards_3prime = 0;
    while towards_3prime < to_the_first_base
        && base(first + towards_3prime)
            .is_some_and(|leaving| base(last + towards_3prime + 1) == Some(leaving))
    {
        towards_3prime += 1;
    }
    let to_the_last_base = if back_cut > 0 { last - tr_len } else { 0 };
    let mut towards_5prime = 0;
    while towards_5prime < to_the_last_base
        && base(last - towards_5prime)
            .is_some_and(|leaving| base(first - towards_5prime - 1) == Some(leaving))
    {
        towards_5prime += 1;
    }
    let innermost_at_the_first_base = (first + towards_3prime, last + towards_3prime);
    let innermost_at_the_last_base = (first - towards_5prime, last - towards_5prime);
    let over_the_first_base =
        front_cut > 0 && innermost_at_the_first_base.0 < 1 && innermost_at_the_first_base.1 >= 1;
    let over_the_last_base = back_cut > 0
        && innermost_at_the_last_base.0 <= tr_len
        && innermost_at_the_last_base.1 > tr_len;
    over_the_first_base || over_the_last_base
}

/// The description the 3' rule gives the bases an allele adds, where the run
/// of points at which they give the same sequence reaches the transcript's
/// first or last base, with the bases the rule carries them from `point`.
/// `inserted` stands on the transcript's strand directly after slice position
/// `point` of a transcript of `tr_len` bases; the run is read through
/// `base_at`, the base at a slice position along the chromosome, past either
/// end of the transcript.
///
/// Each walk ends where the answer is settled. A run that passes the point
/// directly after the last base has no description here wherever it ends, so
/// the walk towards the 3' end reads no base beyond the one directly past the
/// last base. Whether the run starts at or before the first base, and whether
/// a whole copy of the added bases lies between its first point and its last,
/// read the same for every first point as many bases before the first base as
/// the allele adds or farther, so the walk towards the 5' end reads no more
/// bases before the first base than the allele adds.
///
/// A run whose most 3' point lies directly past the last base, the added
/// bases repeating the transcript's last bases, is a duplication of those
/// bases (HGVS nomenclature, DNA duplication: a change that can be described
/// as a duplication is described as one). A run that starts at or before the
/// first base and ends between two bases of the transcript is described at
/// its most 3' point (general recommendations, the 3' rule): a duplication
/// where the copied bases all lie in the transcript, an insertion where the
/// added bases copy nothing. `None` for every other run, among them one
/// whose copied bases start before the first base and one that runs past
/// the last base.
fn added_bases_at_a_transcript_end(
    tr_len: i64,
    base_at: impl Fn(i64) -> Option<u8>,
    point: i64,
    inserted: &[u8],
) -> Option<(HgvscNotation, i64)> {
    let len = inserted.len() as i64;
    if len == 0 {
        return None;
    }
    let mut carried = 0;
    while point + carried <= tr_len
        && base_at(point + carried + 1) == Some(inserted[(carried % len) as usize])
    {
        carried += 1;
    }
    let mut back = 0;
    while point - back > -len
        && base_at(point - back) == Some(inserted[(len - 1 - back % len) as usize])
    {
        back += 1;
    }
    let (first_point, last_point) = (point - back, point + carried);
    let mut alt_seq = inserted.to_vec();
    alt_seq.rotate_left((carried % len) as usize);
    let (start, end, kind) = if last_point == tr_len {
        if len > tr_len || tr_len - first_point < len {
            return None;
        }
        (tr_len - len + 1, tr_len, HgvscKind::Dup)
    } else if first_point <= 0 && 0 < last_point && last_point < tr_len {
        if last_point >= len {
            (last_point - len + 1, last_point, HgvscKind::Dup)
        } else if last_point - first_point < len {
            (last_point, last_point + 1, HgvscKind::Ins)
        } else {
            return None;
        }
    } else {
        return None;
    };
    let notation = HgvscNotation {
        start,
        end,
        ref_seq: Vec::new(),
        alt_seq,
        kind,
    };
    Some((notation, carried))
}

/// The bases an allele adds to a reference it keeps whole, on the
/// transcript's strand, with the slice position they stand directly after:
/// the most 3' point within the reference's span, which starts at slice
/// position `raw_start`. `None` for an allele that changes or drops a base of
/// its reference.
fn bases_added_to_a_kept_reference<'a>(
    reference: &[u8],
    alt: &'a [u8],
    raw_start: i64,
) -> Option<(i64, &'a [u8])> {
    if reference.is_empty() || alt.len() <= reference.len() {
        return None;
    }
    let (lead, tail) = shared_ends_aligned_3prime(alt, reference);
    (lead + tail == reference.len())
        .then(|| (raw_start + lead as i64 - 1, &alt[lead..alt.len() - tail]))
}

/// One exon of the transcript as `_get_cDNA_position` walks them: genomic span
/// and cDNA span, from the transcript mapper's pairs.
struct ExonSpan {
    start: i64,
    end: i64,
    cdna_start: i64,
    cdna_end: i64,
}

impl ExonSpan {
    /// The mapper's exon pairs sorted by genomic start, Ensembl VEP's `_sorted_exons`.
    /// On a model with an `_rna_edit` that inserts or deletes bases the cached
    /// pairs fold the edit (an exon split around it, the pairs downstream
    /// renumbered), so the walk gives the edited model's numbering, where
    /// `_exon_cdna_start` and `_exon_cdna_end` (TranscriptVariationAllele.pm
    /// 2881-2905) give the unedited one: `Exon::cdna_start` maps through a
    /// transcript mapper built at run time, and the cache dumper dropped the
    /// transcript's `edits_enabled` and `transcript_mapper`
    /// (`Pipeline/DumpVEP/Dumper/Core.pm` 155-173, ensembl-vep release/116.3).
    fn in_genomic_order(pairs: &[MapperPair]) -> Vec<ExonSpan> {
        let mut exons: Vec<ExonSpan> = pairs
            .iter()
            .map(|p| ExonSpan {
                start: p.to_start as i64,
                end: p.to_end as i64,
                cdna_start: p.from_start as i64,
                cdna_end: p.from_end as i64,
            })
            .collect();
        exons.sort_by_key(|e| e.start);
        exons
    }
}

/// The numbering of an HGVSc on a model the genome is the reference of
/// (`coding::genome_is_reference`): each exon's bases counted along the genome
/// alone, and the coding bounds at the genomic bases the cached ones map to.
/// The `_rna_edit` such a model keeps marks a base that translation reads
/// twice or passes over. It adds no base to the transcript and takes none
/// out, so a position past it is counted as the transcript's own sequence
/// holds it, where the cached mapper pairs fold the edit in. Ensembl VEP
/// numbers a one-base coding substitution there at the row's CDS start
/// (TranscriptVariationAllele.pm 1519-1522) and any other notation along the
/// exons plus the misalignment offset (1524-1525, 1594-1624), each of which
/// counts the edit.
struct GenomeCount {
    exons: Vec<ExonSpan>,
    coding: Option<(i64, i64)>,
}

impl GenomeCount {
    fn of(
        transcript: &Transcript,
        cached: &[ExonSpan],
        coding: Option<(i64, i64)>,
        reverse: bool,
    ) -> Self {
        let mut spans: Vec<(i64, i64)> = transcript
            .exons
            .iter()
            .map(|exon| (exon.start as i64, exon.end as i64))
            .collect();
        spans.sort_by_key(|(start, _)| *start);
        if reverse {
            spans.reverse();
        }
        let mut counted = 0;
        let mut exons: Vec<ExonSpan> = spans
            .into_iter()
            .map(|(start, end)| {
                let span = ExonSpan {
                    start,
                    end,
                    cdna_start: counted + 1,
                    cdna_end: counted + (end - start + 1),
                };
                counted = span.cdna_end;
                span
            })
            .collect();
        exons.sort_by_key(|e| e.start);
        let mut count = GenomeCount {
            exons,
            coding: None,
        };
        count.coding = coding.map(|(start, end)| {
            (
                count.number(cached, start, reverse),
                count.number(cached, end, reverse),
            )
        });
        count
    }

    /// The count of the genomic base the cached cDNA position `cdna` maps to.
    /// A position the cached pairs give no genomic base, the one an edit
    /// adds, stands for the base of the next cached position, read once more.
    fn number(&self, cached: &[ExonSpan], cdna: i64, reverse: bool) -> i64 {
        let holding = cached
            .iter()
            .find(|e| e.cdna_start <= cdna && cdna <= e.cdna_end);
        let next = || {
            cached
                .iter()
                .filter(|e| e.cdna_start > cdna)
                .min_by_key(|e| e.cdna_start)
        };
        let Some((pair, cached_cdna)) = holding
            .map(|pair| (pair, cdna))
            .or_else(|| next().map(|pair| (pair, pair.cdna_start)))
        else {
            return cdna;
        };
        let genomic = if reverse {
            pair.end - (cached_cdna - pair.cdna_start)
        } else {
            pair.start + (cached_cdna - pair.cdna_start)
        };
        self.exons
            .iter()
            .find(|e| e.start <= genomic && genomic <= e.end)
            .map_or(cdna, |e| {
                if reverse {
                    e.cdna_start + (e.end - genomic)
                } else {
                    e.cdna_start + (genomic - e.start)
                }
            })
    }
}

/// `_get_cDNA_position` (TranscriptVariationAllele.pm 2755): the HGVS cDNA
/// coordinate of 1-based slice position `pos`, exonic as a plain coordinate,
/// intronic as the nearest exon boundary with a `+` or `-` distance (the
/// upstream exon on a tie), then made relative to the start codon (`-` before it)
/// and the stop codon (`*` after it) on a coding transcript. A position outside
/// the transcript has no coordinate: Ensembl VEP reads the exons in genomic
/// order, so a position before the first exon in that order is intronic to it
/// (2800-2827) and prints an offset of one off the transcript's first base on a
/// forward strand (`-49-1`) or its last base on a reverse strand (`83+1`), while
/// a position past the last exon prints nothing (2834)
/// (`docs/divergences/hgvs-clipped-reference-against-whole-alternate.md`).
fn perl_cdna_position(
    slice: &TranscriptSlice<'_>,
    exons: &[ExonSpan],
    coding: Option<(i64, i64)>,
    pos: i64,
) -> Option<String> {
    let g = slice.genomic(pos);
    let mut coord: Option<i64> = None;
    let mut offset: Option<(u8, i64)> = None;
    for (i, exon) in exons.iter().enumerate() {
        if g > exon.end {
            continue;
        }
        if g >= exon.start {
            coord = Some(if slice.reverse {
                exon.cdna_start + (exon.end - g)
            } else {
                exon.cdna_start + (g - exon.start)
            });
            break;
        }
        let prev = exons.get(i.checked_sub(1)?)?;
        let updist = (g - prev.end).abs();
        let downdist = (exon.start - g).abs();
        if updist < downdist || (updist == downdist && !slice.reverse) {
            coord = Some(if slice.reverse {
                prev.cdna_start
            } else {
                prev.cdna_end
            });
            offset = Some((if slice.reverse { b'-' } else { b'+' }, updist));
        } else {
            coord = Some(if slice.reverse {
                exon.cdna_end
            } else {
                exon.cdna_start
            });
            offset = Some((if slice.reverse { b'+' } else { b'-' }, downdist));
        }
        break;
    }
    let mut coord = coord?;
    let mut prefix = "";
    let mut offset_text = offset.map(|(sign, d)| format!("{}{d}", sign as char));
    if let Some((start_codon, stop_codon)) = coding {
        if coord > stop_codon {
            coord -= stop_codon;
            prefix = "*";
        } else if coord == stop_codon && offset.is_some() {
            // Ensembl VEP drops the coordinate and the `+` here: the stop codon's last base
            // plus an intronic distance prints as `*` and the bare distance.
            prefix = "*";
            offset_text = offset_text.map(|t| t.replace('+', ""));
            return Some(format!("{prefix}{}", offset_text.unwrap_or_default()));
        }
        if prefix.is_empty() {
            coord += i64::from(coord >= start_codon);
            coord -= start_codon;
        }
    }
    Some(format!(
        "{prefix}{coord}{}",
        offset_text.unwrap_or_default()
    ))
}

/// The `(exon coordinate, intron offset)` pair Ensembl VEP reads back from an HGVS
/// coordinate string with `m/(\-?[0-9]+)\+?(\-?[0-9]+)?/`: the first signed
/// integer (a `*` prefix is skipped) and the signed integer after it, if any.
fn hgvs_coordinate_parts(text: &str) -> (i64, i64) {
    let bytes = text.as_bytes();
    let mut i = 0;
    while i < bytes.len() && !(bytes[i].is_ascii_digit() || bytes[i] == b'-') {
        i += 1;
    }
    let read_int = |from: usize| -> (Option<i64>, usize) {
        let mut j = from;
        if j < bytes.len() && bytes[j] == b'-' {
            j += 1;
        }
        let digits_from = j;
        while j < bytes.len() && bytes[j].is_ascii_digit() {
            j += 1;
        }
        if j == digits_from {
            return (None, from);
        }
        (text[from..j].parse().ok(), j)
    };
    let (exon, mut next) = read_int(i);
    let Some(exon) = exon else {
        return (0, 0);
    };
    if next < bytes.len() && bytes[next] == b'+' {
        next += 1;
    }
    let (offset, _) = read_int(next);
    (exon, offset.unwrap_or(0))
}

/// `format_hgvs_string` (Utils/Sequence.pm 635).
fn format_hgvs_string(
    ref_name: &str,
    numbering: char,
    start: &str,
    end: &str,
    n: &HgvscNotation,
) -> String {
    let coordinates = if start == end {
        start.to_string()
    } else {
        format!("{start}_{end}")
    };
    let alt = String::from_utf8_lossy(&n.alt_seq);
    let substitution = || format!("{start}{}>{alt}", String::from_utf8_lossy(&n.ref_seq));
    let body = match &n.kind {
        HgvscKind::Sub => substitution(),
        HgvscKind::Inv if n.ref_seq.len() == 1 => substitution(),
        HgvscKind::Del => format!("{coordinates}del"),
        HgvscKind::Inv => format!("{coordinates}inv"),
        HgvscKind::Dup => format!("{coordinates}dup"),
        HgvscKind::Delins => format!("{coordinates}delins{alt}"),
        HgvscKind::Ins => format!("{coordinates}ins{alt}"),
        HgvscKind::Multiple(m) => format!("{coordinates}[{m}]"),
    };
    format!("{ref_name}:{numbering}.{body}")
}

/// The reference-sequence name of an HGVS string, `hgvs_transcript`
/// (TranscriptVariationAllele.pm 1499) and `hgvs_protein` (1759): the stable id
/// with its version appended unless there is none, the id already ends in
/// `.<digits>`, or the id is an LRG identifier.
fn versioned_reference_name(stable_id: &str, version: Option<u32>) -> String {
    let already_versioned = stable_id
        .rsplit_once('.')
        .is_some_and(|(_, tail)| !tail.is_empty() && tail.bytes().all(|b| b.is_ascii_digit()));
    match version {
        Some(version) if !already_versioned && !stable_id.contains("LRG") => {
            format!("{stable_id}.{version}")
        }
        _ => stable_id.to_string(),
    }
}

/// The flags that shape the HGVS notations of a run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct HgvsOptions {
    /// `--shift_hgvs`: describe an insertion or deletion in HGVSc and HGVSp at its
    /// most 3' position on the transcript (on by default); off, the notations
    /// describe the allele where the input placed it and a reverse-strand
    /// duplication is looked up 5' first (TranscriptVariationAllele.pm 1465).
    pub shift_hgvs: bool,
    /// `--hgvsp_use_prediction`: the protein description after `p.` is wrapped in
    /// parentheses, `ENSP00000000001.1:p.(Ala2Asp)` (`_get_hgvs_protein_format`,
    /// TranscriptVariationAllele.pm 1908 and 2040).
    pub hgvsp_use_prediction: bool,
    /// `--remove_hgvsp_version`: the HGVSp reference name is the translation's
    /// stable id without its version appended (TranscriptVariationAllele.pm
    /// 1759-1760); HGVSc keeps its version.
    pub remove_hgvsp_version: bool,
    /// `--shift_3prime 1`: the consequences were computed on the shifted allele,
    /// so the row has a notation when the shifted span overlaps the transcript
    /// (OutputFactory.pm's `within_feature`, 1730), and `hgvs_transcript` clears
    /// the predicate cache once the shift moved it (TranscriptVariationAllele.pm
    /// 1469), so `hgvs_protein` reads `stop_lost`, `start_lost` and the
    /// frameshift guards on the shifted allele too; off, both read the allele as
    /// the input placed it, and `coding::perl_hgvs_protein` the shifted one where
    /// the input placed it in an intron or before the coding start.
    pub shift_3prime: bool,
}

impl Default for HgvsOptions {
    fn default() -> Self {
        Self {
            shift_hgvs: true,
            hgvsp_use_prediction: false,
            remove_hgvsp_version: false,
            shift_3prime: false,
        }
    }
}

/// The HGVS output of one variant-transcript pair under `--hgvs`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct HgvsNotation {
    /// HGVSc, `generate_hgvsc`.
    pub hgvsc: Option<String>,
    /// HGVSp, `generate_hgvsp`.
    pub hgvsp: Option<String>,
    /// VEP's `HGVS_OFFSET`: the bases the insertion or deletion was shifted 3'
    /// along the transcript, negative on the reverse strand, present only when
    /// the shift moved it and at least one notation was produced.
    pub offset: Option<i64>,
    /// Whether `offset` counts a move of the HGVSc alone: the record does not
    /// shift, and the 3' rule places the bases its allele adds beside a base
    /// it keeps, on a model the genome is the reference of
    /// (`coding::genome_is_reference`), along the genome's bases, or carries
    /// them along a run that reaches a transcript's end
    /// (`added_bases_at_a_transcript_end`). The notation is the one the same
    /// bases take written as a minimised insertion, and it carries that
    /// record's offset.
    pub offset_of_hgvsc_alone: bool,
    /// VEP's `REFSEQ_OFFSET`: the misalignment offset `hgvs_transcript` computed
    /// for an `NM_` or `XM_` model with an `_rna_edit` (`refseq_misalignment_offset`,
    /// TranscriptVariationAllele.pm 1516-1517, 1621), when it is not 0
    /// (OutputFactory.pm:1778).
    pub refseq_offset: Option<i64>,
    /// The intron offsets of the HGVSc start and end coordinates
    /// (`_hgvs_intron_start_offset`, `_hgvs_intron_end_offset`,
    /// TranscriptVariationAllele.pm 1530-1545): the signed distance after each
    /// position's exon coordinate as read before the two positions are ordered,
    /// 0 for an exonic position; present whenever `hgvsc` is.
    pub intron_offsets: Option<(i64, i64)>,
}

/// What `hgvs_transcript` leaves on the allele: the notation, the
/// misalignment offset it computed on the way (`refseq_misalignment_offset`),
/// which stays set when the position lookup after it yields no notation
/// (TranscriptVariationAllele.pm 1517, 1527), and the intron offsets of the
/// notation's two coordinates. `read_where_written` is set for an insertion
/// the 3' shift moved that is described where the record wrote it all the
/// same, on the bases the model holds there and the genome lacks.
/// `moved_on_the_genome` is vep-rs's own: the
/// bases the 3' rule moves the notation of an unshifted record
/// (`inserted_at_model_bases_on_the_genome`,
/// `added_bases_at_a_transcript_end`), 0 for every other notation.
#[derive(Debug, Default)]
struct Hgvsc {
    notation: Option<String>,
    refseq_offset: Option<i64>,
    intron_offsets: Option<(i64, i64)>,
    read_where_written: bool,
    moved_on_the_genome: i64,
}

/// Whether the row is annotated outside the transcript: `within_feature`
/// (BaseVariationFeatureOverlapAllele.pm 460-481) read on the span the row was
/// annotated at, the shifted span under `--shift_3prime 1` and the span as
/// written otherwise, so an insertion butting against either end is outside.
fn annotated_outside(
    variant: &InputVariant,
    transcript: &Transcript,
    shifted: Option<(u64, u64)>,
    options: HgvsOptions,
) -> bool {
    let (annotated_start, annotated_end) = match shifted {
        Some(span) if options.shift_3prime => span,
        _ => (variant.start, variant.end),
    };
    annotated_end < transcript.start || annotated_start > transcript.end
}

/// Whether a row outside the transcript can carry an HGVSc: its allele adds
/// bases, which the 3' rule may place on the transcript.
pub(crate) fn adds_bases(variant: &InputVariant) -> bool {
    let alt = variant.alt_allele();
    let reference = variant.ref_allele.as_slice();
    alt != b"-" && (reference == b"-" || alt.len() > reference.len())
}

/// HGVSc, HGVSp and the shift offset for a variant-transcript pair, the indel
/// shift computed once for both notations; under `--shift_hgvs 0` nothing shifts
/// and no offset is reported. `variant` is the allele as the input gave it, on
/// its own strand: HGVSc reads that strand against the transcript's, HGVSp the
/// allele on the forward strand. An insertion that HGVSc reads where the
/// record wrote it, on bases the model holds and the genome lacks, has its
/// HGVSp read there too and reports no offset: neither notation moved. Where
/// those bases open the coding sequence the HGVSp reads their own coding
/// span (`coding::hgvs_protein_where_written`).
///
/// Both notations are read at one span, the one `hgvsc_shift` gives. An
/// insertion whose run of placements is longer than `perform_shift` walks
/// has its HGVSp read at the run's last placement, where its HGVSc is
/// described: none where that placement lies in an intron, in a UTR or past
/// the transcript's last base, and the notation `coding::perl_hgvs_protein`
/// gives the insertion there where it lies between two coding bases.
/// Ensembl VEP's `hgvs_protein` reads the span its own walk reaches
/// (TranscriptVariationAllele.pm 1711, 1715-1721) and has no notation where
/// that span has no translation start or end (1730-1738), so the walk's
/// bound, a count that is no property of the sequence, decides whether its
/// row carries one.
///
/// A record that does not shift, whose HGVSc the 3' rule moves along the
/// genome's bases all the same (`inserted_at_model_bases_on_the_genome`,
/// `added_bases_at_a_transcript_end`),
/// carries the offset of that move: the notation is the one the minimised
/// insertion of the same bases takes there, whose own shift both engines
/// print. Ensembl VEP leaves such a record's bases at the point and prints no
/// offset.
pub fn generate_hgvs(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
    options: HgvsOptions,
) -> HgvsNotation {
    let shifted = hgvs_shift(variant, transcript, reference_fasta, options.shift_hgvs);
    let described_at = hgvsc_shift(variant, transcript, reference_fasta, shifted, options);
    let hgvsc = hgvsc_at(variant, transcript, reference_fasta, described_at, options);
    // An insertion the HGVSc reads where its record writes it is read there
    // by both notations.
    let described_at = described_at.filter(|_| !hgvsc.read_where_written);
    // A row outside the transcript has no protein notation, whatever the
    // 3' rule makes of its HGVSc.
    let hgvsp = if annotated_outside(variant, transcript, described_at, options) {
        None
    } else {
        hgvsp_at(
            variant,
            transcript,
            reference_fasta,
            described_at,
            hgvsc.read_where_written,
            options,
        )
    };
    let along_the_transcript = |length: i64| {
        if transcript.strand == Strand::Reverse {
            -length
        } else {
            length
        }
    };
    // The offset is the bases the 3' rule moves the span both notations are
    // read at, on a row that carries either of them.
    let offset = described_at
        .filter(|_| hgvsc.notation.is_some() || hgvsp.is_some())
        .map(|(start, _)| along_the_transcript(start.abs_diff(variant.start) as i64));
    let offset_of_hgvsc_alone = offset.is_none() && hgvsc.moved_on_the_genome != 0;
    HgvsNotation {
        hgvsc: hgvsc.notation,
        hgvsp,
        offset: offset.or_else(|| {
            offset_of_hgvsc_alone.then(|| along_the_transcript(hgvsc.moved_on_the_genome))
        }),
        offset_of_hgvsc_alone,
        refseq_offset: hgvsc.refseq_offset,
        intron_offsets: hgvsc.intron_offsets,
    }
}

/// Generate HGVSc notation for a variant-transcript pair: Ensembl VEP's
/// `hgvs_transcript` (TranscriptVariationAllele.pm 1375) for an Ensembl
/// transcript.
///
/// Returns `None` when the alternate allele carries a character outside
/// `ACGT-` or equals the reference, when the variant does not overlap the
/// transcript's genomic span (so an upstream or downstream variant has no HGVSc),
/// when the 3'-shifted span leaves it, or when the alleles agree at the shifted
/// position. A variant whose span, after the shift, runs over either edge of the
/// transcript is described over the bases inside it, each end clamped onto the
/// transcript (`_var2transcript_slice_coords`, 2732-2733), and its alternate
/// allele loses the bases the clamp cut (`alt_cut_as_reference`): a
/// substitution running over an edge is described by the bases that change
/// inside, and one whose bases inside are unchanged has no notation, where
/// Ensembl VEP reads the whole alternate allele against the bases inside
/// (1474-1483). An insertion or
/// deletion is described at its most 3' position on the transcript strand when
/// a FASTA is present; the type and the displayed span follow
/// `hgvs_variant_notation`,
/// shared bases are clipped, an exonic single-base substitution inside the CDS
/// takes its CDS coordinate and every other position goes through
/// `_get_cDNA_position`. Without a FASTA the reference bases come from the
/// variant and no duplication is detected.
///
/// On a transcript with an `_rna_edit` attribute the notation spells the
/// reference and alternate alleles relative to the transcript
/// (`$ref->feature_seq`, `$self->feature_seq`, TranscriptVariationAllele.pm
/// 1510-1514), the reference being the allele the row was analysed against and
/// a shifted insertion's alternate allele rotated along the shift, as
/// `shift_feature_seqs` (1330-1345) leaves it on a coding row under
/// `--shift_3prime 1` (`codon`, 813-820); on every other row Ensembl VEP spells
/// the unrotated allele at the shifted position (1513), a defect vep-rs does not
/// reproduce. Every position is numbered
/// against the edited model through the cached mapper pairs. Ensembl VEP
/// numbers the exons unedited and moves an `NM_` or
/// `XM_` model's position by the misalignment offset of the indel edits
/// upstream of the allele (1516-1517, 1524-1525), an offset an intronic allele
/// never gets (1608-1610), so its intronic distance counts from the unedited
/// exon end; vep-rs counts it from the edited one. The CDS coordinate of a
/// coding-exon substitution is read from the edited mapper on both engines
/// (1519-1522). An insertion written on bases such a model holds and the
/// genome lacks is described on them, where it is written even when the
/// genome's shift moves it past them, and a duplication the genome shows
/// across such bases is numbered from the model's own. A model the genome is
/// the reference of (`coding::genome_is_reference`) is numbered along the
/// genome's exons instead, that substitution with it (`GenomeCount`).
pub fn generate_hgvsc(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<String> {
    let options = HgvsOptions::default();
    let shifted = hgvs_shift(variant, transcript, reference_fasta, options.shift_hgvs);
    let described_at = hgvsc_shift(variant, transcript, reference_fasta, shifted, options);
    hgvsc_at(variant, transcript, reference_fasta, described_at, options).notation
}

/// `generate_hgvsc` with the 3' shift already computed. Under `--shift_hgvs 0`
/// (`options.shift_hgvs` false) an allele given on the strand opposite the
/// transcript's looks a duplication up 5' first (`$lookup_order = -1` when
/// `$vf->strand` differs from the transcript's, TranscriptVariationAllele.pm
/// 1462-1466); under `--shift_3prime 1` (`options.shift_3prime`) the row was
/// annotated at `shifted`, which decides whether it gets a notation at all.
fn hgvsc_at(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
    shifted: Option<(u64, u64)>,
    options: HgvsOptions,
) -> Hgvsc {
    let mut hgvsc = Hgvsc::default();
    let (notation, moved_on_the_genome) = hgvsc_notation(
        variant,
        transcript,
        reference_fasta,
        shifted,
        options,
        &mut hgvsc,
    )
    .unzip();
    hgvsc.notation = notation;
    hgvsc.refseq_offset = hgvsc.refseq_offset.filter(|offset| *offset != 0);
    hgvsc.moved_on_the_genome = moved_on_the_genome.unwrap_or(0);
    hgvsc
}

/// The notation of [`hgvsc_at`], `None` where `hgvs_transcript` returns
/// `undef`, with the bases the 3' rule moves it where the record itself
/// does not shift (`inserted_at_model_bases_on_the_genome`), 0 anywhere
/// else; what else the row carries is left in `found`.
/// `found.refseq_offset` is set once the whole alternate allele is known to
/// differ from the bases of the clamped span, before the notation is typed,
/// so it stands on a row whose cut alleles or position lookup yield no
/// notation;
/// `found.intron_offsets` is set once both cDNA positions are known, which is
/// exactly when a notation follows; `found.read_where_written` is set as soon
/// as the shift is set aside, whether or not a notation follows.
fn hgvsc_notation(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
    shifted: Option<(u64, u64)>,
    options: HgvsOptions,
    found: &mut Hgvsc,
) -> Option<(String, i64)> {
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let pairs = &mapper.exon_coord_mapper.pairs;
    if pairs.is_empty() {
        return None;
    }

    let alt_allele = variant.alt_allele();
    if alt_allele
        .iter()
        .any(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'-'))
    {
        return None;
    }
    if alt_allele == variant.ref_allele.as_slice() {
        return None;
    }

    let reverse = transcript.strand == Strand::Reverse;
    let has_an_edit = transcript
        .refseq_alignment()
        .is_some_and(|alignment| alignment.first_edit_bases.is_some());

    // `feature_seq` (VariationFeatureOverlapAllele.pm 246-262): the allele as
    // given, reverse-complemented when the input's strand is not the
    // transcript's, its `-` kept.
    let feature_seq = |mut allele: Vec<u8>, given_on: Strand| {
        allele.make_ascii_uppercase();
        if given_on != transcript.strand && allele != b"-" {
            allele = crate::coding::reverse_complement(&allele);
        }
        allele
    };
    // `hgvs_transcript` reverse-complements an allele given on the strand
    // opposite the transcript's (TranscriptVariationAllele.pm 1463-1464) and
    // `hgvs_variant_notation` strips its gap characters before reading it.
    let onto_transcript = |allele: Vec<u8>, given_on: Strand| {
        let mut allele = feature_seq(allele, given_on);
        allele.retain(|&b| b != b'-');
        allele
    };

    // OutputFactory.pm asks for HGVS only `within_feature` (1730), the overlap of
    // the transcript with the span the row was annotated at: under
    // `--shift_3prime 1` the shifted span, which `_bvfo_preds` moves before the
    // test (BaseVariationFeatureOverlapAllele.pm 460-481), otherwise the span as
    // written. A row annotated outside the transcript has a notation only for
    // added bases that the 3' rule places on the transcript
    // (`added_bases_at_a_transcript_end`): the description of a change does
    // not turn on which of its placements the record writes.
    let outside = annotated_outside(variant, transcript, shifted, options);
    if outside && (!options.shift_hgvs || options.shift_3prime || has_an_edit) {
        return None;
    }

    let slice = TranscriptSlice {
        chr: &variant.chr,
        tr_start: transcript.start as i64,
        tr_end: transcript.end as i64,
        reverse,
        fasta: reference_fasta,
    };
    let (vf_start, vf_end) = (variant.start as i64, variant.end as i64);
    let exons = ExonSpan::in_genomic_order(pairs);

    // An insertion written between two exonic bases the model numbers more
    // than one apart stands on the bases the model holds between them and
    // the genome lacks: the row's cDNA, CDS and protein positions, its codons
    // and its term are read there. The 3' shift runs along the genome
    // (`_genomic_shift`, TranscriptVariationAllele.pm 422-465), which has no
    // such bases, so it can carry the insertion past them, and Ensembl VEP
    // types and places the notation at the shifted point (1446, 1474-1478):
    // the allele added after the model's bases with every one of them kept,
    // a sequence longer by those bases than the one the row's own span
    // states. Such a row is described where its record is written, the shift
    // set aside. Under `--shift_3prime 1` the row itself is annotated at the
    // shifted point, and a model whose edit failed is read against the
    // genome; both keep the shifted notation.
    let written_on_model_bases = shifted
        .filter(|&(shifted_start, _)| {
            shifted_start != variant.start
                && variant.variant_class == VariantClass::Insertion
                && !options.shift_3prime
                && transcript.bam_edit_status != Some(BamEditStatus::Failed)
        })
        .and_then(|_| {
            let before_the_point = slice.position(vf_start).min(slice.position(vf_end));
            model_bases_between(transcript, &slice, &exons, before_the_point)
        });
    found.read_where_written = written_on_model_bases.is_some();
    let shifted = shifted.filter(|_| !found.read_where_written);

    // A shifted insertion is described by its rotated allele, on the forward
    // strand whatever strand the input gave it on; every other allele as given.
    let (described_alt, described_alt_strand) =
        match normalized_hgvs_alt_allele(variant, transcript, shifted) {
            Some(rotated) => (rotated, Strand::Forward),
            None => (alt_allele.to_vec(), variant.strand),
        };
    let alt_seq = onto_transcript(described_alt.clone(), described_alt_strand);
    let fallback_ref = onto_transcript(variant.ref_allele.clone(), variant.strand);

    // The notation describes the span the 3' shift moved the variant to, the
    // span as written when nothing moved. Ensembl VEP reads the written span
    // instead: `_var2transcript_slice_coords` maps it
    // (TranscriptVariationAllele.pm 2708-2713), refuses one wholly outside the
    // transcript (2729) and clamps its ends onto the transcript (2732-2733), then
    // `hgvs_transcript` refuses a clamped end the offset carries past the last
    // base (1472) and adds the offset to both ends (1477-1478). So an insertion
    // or deletion the shift carries into a transcript from wholly before it, or
    // over its last base, gets no notation, and a deletion straddling the first
    // base is clamped onto that base before the offset moves it, two deleted
    // bases printing as one (`n.14del` for the span `n.13_14del`). Read at the
    // shifted span, each row takes the notation both engines print for the same
    // allele written there.
    let (hgvs_start, hgvs_end) = shifted.unwrap_or((variant.start, variant.end));
    let (hgvs_start, hgvs_end) = (hgvs_start as i64, hgvs_end as i64);
    // On the reverse strand the ends swap.
    let (raw_start, raw_end) = if reverse {
        (slice.position(hgvs_end), slice.position(hgvs_start))
    } else {
        (slice.position(hgvs_start), slice.position(hgvs_end))
    };
    // An empty slice has no base to describe, which Ensembl VEP never checks.
    let tr_len = slice.len();
    if tr_len < 1 {
        return None;
    }
    // The bases an allele adds to a reference it keeps whole are placed by the
    // 3' rule like any inserted bases, though its record does not shift.
    // Ensembl VEP leaves them where the clip does, most 3' within the
    // record's span, and types them an insertion there (2204-2284), with a
    // coordinate one past the last base where the span ends on it and one
    // before the first base where the span starts before it (2802,
    // 2820-2827).
    let base_at = |pos: i64| slice.base_at(pos);
    let carried_to_an_end = (options.shift_hgvs && shifted.is_none() && !has_an_edit)
        .then(|| bases_added_to_a_kept_reference(&fallback_ref, &alt_seq, raw_start))
        .flatten()
        .and_then(|(point, inserted)| {
            added_bases_at_a_transcript_end(tr_len, base_at, point, inserted)
        });
    // Of the rows outside, one whose added bases stand, as the record's own
    // shift leaves them, at the point the rule gives them.
    if outside && carried_to_an_end.is_none() {
        let placed_by_the_shift = fallback_ref.is_empty()
            && added_bases_at_a_transcript_end(tr_len, base_at, raw_end, &alt_seq)
                .is_some_and(|(_, carried)| carried == 0);
        if !placed_by_the_shift {
            return None;
        }
    }
    // A span wholly before or wholly past the transcript has no notation (2729):
    // of the spans the gate above admits, one the shift carried wholly past the
    // last base under `--hgvs` alone.
    if carried_to_an_end.is_none()
        && ((raw_start < 1 && raw_end < 1) || (raw_start > tr_len && raw_end > tr_len))
    {
        return None;
    }
    // An end outside the slice is clamped onto its edge (2732-2733), so a
    // deletion running over an edge is described over the bases inside the
    // transcript. A start can only run over the first base and an end over the
    // last, so an insertion's span, the base after the site then the site, stays
    // whole at either edge; Ensembl VEP clamps the written span (1, 0) of an
    // insertion anchored against the first base onto that base and, adding the
    // offset to both ends, writes a substitution of the base the shift moved it
    // to (`n.15A>T` for a `T` carried into a run of `T`) where the sequence
    // shows the duplication `n.14dup`.
    let (slice_start, slice_end) = (raw_start.max(1), raw_end.min(tr_len));
    let front_cut = (slice_start - raw_start) as usize;
    let back_cut = (raw_end - slice_end) as usize;
    // The misalignment offset of an `NM_` or `XM_` model, the length change
    // of the indel edits at or before the allele's cDNA start
    // (`get_misalignment_offset`, TranscriptVariationAllele.pm 1517,
    // 1594-1624). Ensembl VEP computes it once the whole alternate allele,
    // read against the bases of the clamped span (1474-1483), gives a
    // notation, so it stands on a row whose bases inside the transcript are
    // unchanged while a base outside changes, where the alternate cut as the
    // clamp cut the reference types none below. An alternate that equals the
    // bases of the clamped span has no offset.
    if let Some(alignment) = transcript.refseq_alignment().filter(|_| has_an_edit) {
        if transcript.stable_id.starts_with("NM_") || transcript.stable_id.starts_with("XM_") {
            let inside =
                match slice.substr(slice_start, (slice_end - slice_start + 1).max(0) as usize) {
                    Some(inside) => Some(inside),
                    None if slice.fasta.is_none() => fallback_ref
                        .get(front_cut..fallback_ref.len().saturating_sub(back_cut))
                        .map(<[u8]>::to_vec),
                    None => None,
                };
            if inside.is_some_and(|inside| inside != alt_seq) {
                let cdna_start = crate::display::genomic2cdna(
                    mapper,
                    vf_start,
                    vf_end,
                    transcript.strand.as_i8(),
                )
                .first()
                .and_then(|segment| match segment {
                    crate::display::Segment::Coord { start, .. } => Some(*start),
                    crate::display::Segment::Gap { .. } => None,
                });
                found.refseq_offset = Some(alignment.misalignment_offset(cdna_start));
            }
        }
    }
    // A deletion that runs over the transcript's first or last base and has
    // no window wholly inside it has no description on the transcript.
    // Ensembl VEP clamps the span and describes the bases inside
    // (`docs/divergences/hgvs-clipped-reference-against-whole-alternate.md`).
    if options.shift_hgvs
        && kept_allele_deletes_over_a_transcript_end(
            tr_len,
            base_at,
            &fallback_ref,
            &alt_seq,
            raw_start,
            front_cut,
            back_cut,
        )
    {
        return None;
    }
    // `hgvs_transcript` passes the clamped span with the whole alternate allele
    // (1474-1483), and `hgvs_variant_notation` reads the reference from the span
    // alone (Utils/Sequence.pm 526), so Ensembl VEP describes an allele running
    // over an edge as every alternate base replacing the bases inside, a `delins`
    // of a substitution's whole allele, or, when the bases inside are unchanged,
    // an `ins` of the outside bases at a coordinate before the first exon, a
    // defect vep-rs does not reproduce
    // (`docs/divergences/hgvs-clipped-reference-against-whole-alternate.md`).
    let lookup = if variant.strand != transcript.strand && !options.shift_hgvs {
        DupLookup::BeforeFirst
    } else {
        DupLookup::AfterFirst
    };
    let at_a_transcript_end = carried_to_an_end.is_some();
    let (mut n, (typed_start, typed_end), carried_by_the_rule) = match carried_to_an_end {
        Some((n, carried)) => {
            let typed = (n.start, n.end);
            (n, typed, carried)
        }
        None => {
            let alt_seq = alt_cut_as_reference(alt_seq, &fallback_ref, front_cut, back_cut);
            // Without a FASTA the reference bases come from the variant: the
            // bases the clamp cut off are dropped from that fallback too.
            let fallback_ref = fallback_ref
                .get(front_cut..fallback_ref.len().saturating_sub(back_cut))
                .unwrap_or(&fallback_ref);
            let mut n = hgvs_variant_notation(
                &slice,
                fallback_ref,
                alt_seq,
                slice_start,
                slice_end,
                lookup,
            )?;
            let typed = (n.start, n.end);
            if n.kind != HgvscKind::Dup {
                hgvsc_clip_alleles(&mut n);
            }
            (n, typed, 0)
        }
    };
    // The bases the notation's span lost at each end, to the clamp onto the
    // transcript and to the clip.
    let lost_front = front_cut + (n.start - typed_start) as usize;
    let lost_back = back_cut + (typed_end - n.end) as usize;

    let ref_name = versioned_reference_name(&transcript.stable_id, transcript.version);

    let coding = (transcript.translation.is_some() && mapper.cdna_coding_start > 0).then_some((
        mapper.cdna_coding_start as i64,
        mapper.cdna_coding_end as i64,
    ));

    // An `_rna_edit` on the transcript (TranscriptVariationAllele.pm 1507-1517):
    // the alleles are respelled relative to the transcript, the reference as the
    // row was analysed against it, a shifted insertion's alternate
    // allele as the shift rotates it, the bases `shift_feature_seqs` (1330-1345)
    // leaves in `feature_seq` on a coding row under `--shift_3prime 1` (`codon`,
    // 813-820). On every other row `hgvs_transcript` copies the unrotated
    // `feature_seq` over the shifted notation's alternate allele (1513), spelling
    // a sequence the variant never produces; vep-rs writes the rotated allele on
    // every row (`docs/divergences/refseq-edited-model-shifted-insertion-bases.md`).
    // The respelled alleles are clipped as the notation's were: `hgvs_transcript`
    // writes them over the clipped pair whole (1492, 1510-1513), so an allele
    // that shares a base with its reference keeps that base in Ensembl VEP's
    // notation, an `ins` or a `delins` one base too long. They lose the bases
    // the notation's span lost and no other (`clip_respelled_alleles`): the
    // clip reads the genome's bases, which the transcript's differ from at a
    // substitution edit. An
    // `NM_` or `XM_` model's misalignment offset, the length change of the indel
    // edits at or before the allele's cDNA start, is kept for the `REFSEQ_OFFSET`
    // column alone and is not added to the positions below as `hgvs_transcript`
    // adds it (1524-1525): the cached mapper pairs already fold the edits.
    let mut clipped_alt = None;
    let mut numbered_on_model = false;
    if has_an_edit {
        n.ref_seq = feature_seq(variant.ref_allele.clone(), variant.strand);
        clipped_alt = Some(std::mem::replace(
            &mut n.alt_seq,
            feature_seq(described_alt, described_alt_strand),
        ));
        // A base the clamp cut from the span is no base of the transcript,
        // whether or not the alternate changes it: the respelled pair loses
        // those bases as the notation's own alleles do
        // (`alt_cut_as_reference`), which leaves the clip below the bases
        // `_clip_alleles` took. Ensembl VEP writes the whole pair over the
        // clamped span (1510-1513).
        if front_cut + back_cut > 0 && n.ref_seq.len() >= front_cut + back_cut {
            let whole_alt = std::mem::take(&mut n.alt_seq);
            n.alt_seq = alt_cut_as_reference(whole_alt, &n.ref_seq, front_cut, back_cut);
            n.ref_seq.truncate(n.ref_seq.len() - back_cut);
            n.ref_seq.drain(..front_cut);
        }
        match n.kind {
            HgvscKind::Dup => {}
            // A substitution names one base a side (Utils/Sequence.pm 656).
            // A base the clip takes that the model spells differently stays
            // in the respelled pair, and the span keeps its position: the
            // record's bases stand in place of the model's own from that
            // base to the one the record changes. Ensembl VEP writes the
            // whole pair at the one position the clip leaves (1492,
            // 1510-1513).
            //
            // The model's reference allele holds the bases the model has
            // between two bases of the record and the genome lacks
            // (TranscriptVariation.pm 463). They are no bases of the
            // record, so the clip stops at the first of them from either
            // end, and the span of a reference that holds them is numbered
            // on the model, from the first base the pair keeps to the last:
            // it starts at the first base the genome lacks where the clip
            // takes the base before them. Ensembl VEP writes them into a
            // substitution at the one position the clip leaves.
            HgvscKind::Sub if front_cut + back_cut == 0 => {
                let reference_holds_model_bases = n.ref_seq.len() > n.alt_seq.len();
                let joined = |pos: i64| {
                    !reference_holds_model_bases
                        || model_bases_between(transcript, &slice, &exons, pos).is_none()
                };
                let front = (0..lost_front as i64)
                    .take_while(|i| *i == 0 || joined(typed_start + i - 1))
                    .count();
                let back = (0..lost_back as i64)
                    .take_while(|i| *i == 0 || joined(typed_end - i))
                    .count();
                let (kept_front, kept_back) =
                    clip_substitution_alleles(&mut n.ref_seq, &mut n.alt_seq, front, back);
                let dropped = ((front - kept_front) as i64, (back - kept_back) as i64);
                n.start = typed_start + dropped.0;
                n.end = typed_end - dropped.1;
                let cdna = |pos: i64| exonic_cdna_position(&slice, &exons, pos);
                let on_model = cdna(typed_start)
                    .zip(cdna(typed_end))
                    .map(|((_, first), (_, last))| (first + dropped.0, last - dropped.1))
                    .filter(|(first, last)| {
                        n.ref_seq.len() > n.alt_seq.len()
                            && last - first + 1 == n.ref_seq.len() as i64
                    });
                if let Some(span) = on_model {
                    (n.start, n.end) = span;
                    numbered_on_model = true;
                }
                let keeps_a_lost_base = dropped != (lost_front as i64, lost_back as i64);
                if keeps_a_lost_base || numbered_on_model {
                    n.kind = replacement_kind(&n.ref_seq, n.ref_seq.len(), &n.alt_seq);
                }
            }
            // A `delins` or a `del` the clip leaves beside the bases the
            // model holds and the genome lacks takes them too: the
            // record's bases stand in place of the model's own over the
            // row's span, so the span runs on the model from the first
            // base the pair keeps to the last (`span_beside_model_bases`),
            // with the alternate bases the clip leaves. Ensembl VEP keeps
            // the span the clip leaves, a `delins` of the whole alternate
            // allele (1492, 1513; Utils/Sequence.pm 666-667) or a `del` of
            // that span alone (659-663).
            HgvscKind::Delins | HgvscKind::Del if front_cut + back_cut == 0 => {
                let beside = span_beside_model_bases(
                    transcript,
                    &slice,
                    &exons,
                    (typed_start, typed_end),
                    (n.start, n.end),
                    (lost_front, lost_back),
                    &n.ref_seq,
                    &n.alt_seq,
                );
                let (front, back) = beside.map_or((lost_front, lost_back), |(_, taken)| taken);
                clip_respelled_alleles(&mut n.ref_seq, &mut n.alt_seq, front, back);
                if let Some((span, _)) = beside {
                    (n.start, n.end) = span;
                    numbered_on_model = true;
                    n.kind = replacement_kind(&n.ref_seq, n.ref_seq.len(), &n.alt_seq);
                }
            }
            _ => clip_respelled_alleles(
                &mut n.ref_seq,
                &mut n.alt_seq,
                lost_front - front_cut,
                lost_back - back_cut,
            ),
        }
    }

    // An insertion between two exonic bases the model numbers more than one
    // apart stands on the bases the model holds between them and the genome
    // lacks: the mapper gives the row those bases as its cDNA and CDS span
    // (`coding::perl_map_insert`), and the notation names their replacement by
    // the allele, where `ins` between the two would name positions that are
    // not neighbours. Ensembl VEP slides the insertion's two positions by the
    // misalignment offset (TranscriptVariationAllele.pm 1524-1525) and names
    // neither the bases nor their positions. An insertion that repeats the
    // genomic bases beside that point stands on the model's bases all the
    // same: typed on the genomic slice it is a duplication of those genomic
    // bases (Utils/Sequence.pm 577-585), a sequence longer by the model's
    // bases than the one `hgvs_protein` reads there (1732-1733).
    //
    // An allele that keeps a base it shares with its reference spans that
    // base, and the clip leaves it an `ins` of its other bases beside it
    // (1492). Where the model's bases lie on the far side of that point the
    // row is read at the kept base, with the model's bases untouched: the
    // inserted bases stand between the kept base and the model's base next
    // to it, two neighbouring positions of the model. Ensembl VEP numbers
    // the point's two genomic bases from the unedited exons (1524-1525). A
    // span that holds a base on either
    // side of the point holds the model's bases too, which the inserted
    // bases replace.
    //
    // Where the genome is the model's reference
    // (`coding::genome_is_reference`) the bases the model holds at the point
    // are no bases of the transcript: a duplication the genome shows there
    // keeps the genome's notation, and any other insertion at the point is
    // an `ins` after the base 5' of it, the point the row's span gives. The
    // bases an allele that keeps a shared base leaves are inserted at that
    // point too, placed by the 3' rule on the model's sequence without the
    // model's bases (`inserted_at_model_bases_on_the_genome`).
    //
    // An insertion the shift moved is
    // read on those bases whatever the genome types it where it is written,
    // an `ins` or a duplication of the bases after the point.
    //
    // A reference allele kept whole with bases on both sides of the model's
    // spans them on the model, and the row's reference is cut from the model
    // over that span (TranscriptVariation.pm 463). An alternate that is the
    // genomic reference twice is typed a duplication on the genome
    // (Utils/Sequence.pm 602-606), which Ensembl VEP numbers along the
    // unedited exon (TranscriptVariationAllele.pm 1524-1525); on the model
    // the alternate stands in place of the span, the model's own bases among
    // them. A model whose edit failed is read against the genome and keeps
    // the duplication.
    let described_at_a_point = slice_start > slice_end;
    let on_the_genome = crate::coding::genome_is_reference(transcript);
    let at_a_point_on_the_genome = described_at_a_point && on_the_genome;
    let mut moved_on_the_genome = carried_by_the_rule;
    let mut in_walk_positions = false;
    let on_model_bases = written_on_model_bases.or_else(|| match n.kind {
        _ if at_a_transcript_end => None,
        HgvscKind::Dup if at_a_point_on_the_genome => None,
        HgvscKind::Ins | HgvscKind::Dup if described_at_a_point => {
            model_bases_between(transcript, &slice, &exons, slice_end)
        }
        HgvscKind::Ins => model_bases_between(transcript, &slice, &exons, n.end),
        HgvscKind::Dup
            if slice_start <= slice_end
                && transcript.bam_edit_status != Some(BamEditStatus::Failed) =>
        {
            model_span_across_model_bases(transcript, &slice, &exons, slice_start, slice_end)
        }
        _ => None,
    });
    let numbered_on_model = if let Some(bases) = &on_model_bases {
        // On an edited model the bases standing at the point are the ones
        // the clip leaves: a respelled pair that does not share the bases
        // the span lost stays whole.
        let inserted = clipped_alt.unwrap_or_else(|| std::mem::take(&mut n.alt_seq));
        let kept_before = !described_at_a_point && n.end >= slice_start;
        let kept_after = !described_at_a_point && n.start <= slice_end;
        // An allele that keeps one base is the replacement its minimised
        // spelling is where the 3' rule carries that replacement to the
        // exon's last base or past it, so that one change has one notation
        // there however the record spells it. The rule is read with the
        // shift on: the reading of an allele does not turn on
        // `--shift_hgvs`.
        let replaces_to_the_exons_last_base = kept_before != kept_after
            && !on_the_genome
            && replaced_model_bases_notation(transcript, &slice, bases, inserted.clone(), true)
                .is_some_and(|placed| bases.reaches_the_exons_last_base(&placed));
        let beside = |point: i64, inserted: Vec<u8>| {
            inserted_beside_model_bases_notation(
                transcript,
                &slice,
                bases,
                point,
                inserted,
                options.shift_hgvs,
                lookup,
            )
        };
        n = match (kept_before, kept_after) {
            _ if at_a_point_on_the_genome => HgvscNotation {
                start: bases.first - 1,
                end: bases.first,
                ref_seq: Vec::new(),
                alt_seq: inserted,
                kind: HgvscKind::Ins,
            },
            (true, false) | (false, true) if on_the_genome => {
                let (placed, moved) = inserted_at_model_bases_on_the_genome(
                    transcript,
                    &slice,
                    bases,
                    inserted,
                    options.shift_hgvs,
                    lookup,
                );
                moved_on_the_genome = moved;
                placed
            }
            (true, false) if !replaces_to_the_exons_last_base => beside(bases.first - 1, inserted),
            (false, true) if !replaces_to_the_exons_last_base => beside(bases.last, inserted),
            _ => {
                in_walk_positions = true;
                replaced_model_bases_notation(
                    transcript,
                    &slice,
                    bases,
                    inserted,
                    options.shift_hgvs,
                )?
            }
        };
        true
    } else if n.kind == HgvscKind::Dup
        && slice_start > slice_end
        && transcript.bam_edit_status != Some(BamEditStatus::Failed)
    {
        // An insertion beside the model's bases that the genome shows as a
        // duplication of bases on both sides of them is a duplication on the
        // model only of the bases the model holds beside the point. A model
        // whose edit failed is read against the genome and keeps its range.
        duplication_across_model_bases(
            transcript,
            &slice,
            &exons,
            slice_end,
            lookup,
            clipped_alt.as_deref(),
            &mut n,
        )?
    } else {
        numbered_on_model
    };

    let same_pos = n.start == n.end;
    let is_snp = variant.ref_allele.len() == 1
        && alt_allele.len() == 1
        && variant.ref_allele != b"-"
        && alt_allele != b"-";
    let exonic = exons
        .iter()
        .any(|e| vf_start.min(vf_end) <= e.end && vf_start.max(vf_end) >= e.start);
    // A model the genome is the reference of is numbered along the genome's
    // exons, the transcript's own bases; every other model, and one whose
    // record holds no exon list, on the cached pairs.
    let genome_count = (on_the_genome && !transcript.exons.is_empty())
        .then(|| GenomeCount::of(transcript, &exons, coding, reverse));
    let (numbered_exons, numbered_coding) = match &genome_count {
        Some(count) => (count.exons.as_slice(), count.coding),
        None => (exons.as_slice(), coding),
    };
    let numbered = |cdna: i64| {
        genome_count
            .as_ref()
            .map_or(cdna, |count| count.number(&exons, cdna, reverse))
    };
    // The CDS start counts what the cached pairs count: on a model numbered
    // along the genome it moves by what the count differs by at the row's
    // cDNA start, less what it differs by at the coding start.
    let cds_snp_position = if is_snp && exonic && coding.is_some() {
        crate::coding::perl_span(transcript, variant.start, variant.end)
            .and_then(|s| {
                s.cds_start
                    .zip(s.cds_end)
                    .map(|(cds, _)| (cds, s.cdna_start))
            })
            .map(
                |(cds_start, cdna_start)| match (cdna_start, coding, numbered_coding) {
                    (Some(cdna), Some((cached_start, _)), Some((counted_start, _))) => {
                        cds_start - (cdna - numbered(cdna)) + (cached_start - counted_start)
                    }
                    _ => cds_start,
                },
            )
    } else {
        None
    };
    // The notation of replaced model bases counts in the walk's positions
    // (`ModelBases::place`); every other notation on or beside the model's
    // bases counts in cDNA positions.
    let (start, end) = match (&on_model_bases, cds_snp_position) {
        (Some(bases), _) if in_walk_positions => {
            let coordinate = |pos: i64| match bases.place(pos) {
                WalkPlace::Cdna(cdna) => {
                    Some(exonic_hgvs_coordinate(numbered_coding, numbered(cdna)))
                }
                WalkPlace::Intron(pos) => {
                    perl_cdna_position(&slice, numbered_exons, numbered_coding, pos)
                }
            };
            (coordinate(n.start)?, coordinate(n.end)?)
        }
        _ if numbered_on_model => (
            exonic_hgvs_coordinate(numbered_coding, numbered(n.start)),
            exonic_hgvs_coordinate(numbered_coding, numbered(n.end)),
        ),
        (_, Some(cds_start)) => (cds_start.to_string(), cds_start.to_string()),
        (_, None) => {
            let start = perl_cdna_position(&slice, numbered_exons, numbered_coding, n.start)?;
            let end = if same_pos {
                start.clone()
            } else {
                perl_cdna_position(&slice, numbered_exons, numbered_coding, n.end)?
            };
            (start, end)
        }
    };

    let (exon_start, intron_start) = hgvs_coordinate_parts(&start);
    let (exon_end, intron_end) = if same_pos {
        (exon_start, intron_start)
    } else {
        hgvs_coordinate_parts(&end)
    };
    found.intron_offsets = Some((intron_start, intron_end));
    let (start, end) = if (exon_start > exon_end
        || (exon_start == exon_end && intron_start > intron_end))
        && !end.contains('*')
    {
        (end, start)
    } else {
        (start, end)
    };

    let numbering = if coding.is_some() { 'c' } else { 'n' };
    let notation = format_hgvs_string(&ref_name, numbering, &start, &end, &n);
    Some((notation, moved_on_the_genome))
}

fn normalized_hgvs_ref_allele(
    variant: &InputVariant,
    shifted: Option<(u64, u64)>,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<Vec<u8>> {
    if variant.variant_class != VariantClass::Deletion {
        return None;
    }

    let (shifted_start, shifted_end) = shifted?;
    if shifted_start == variant.start && shifted_end == variant.end {
        return None;
    }

    let mut seq = reference_fasta?.sequence(&variant.chr, shifted_start, shifted_end)?;
    for base in &mut seq {
        *base = base.to_ascii_uppercase();
    }
    Some(seq)
}

fn normalized_hgvs_alt_allele(
    variant: &InputVariant,
    transcript: &Transcript,
    shifted: Option<(u64, u64)>,
) -> Option<Vec<u8>> {
    let (motif, shift_steps) = moved_insertion_motif(variant, shifted)?;
    Some(rotate_sequence_for_transcript_3prime(
        &motif,
        transcript.strand,
        shift_steps,
    ))
}

fn rotate_sequence_for_transcript_3prime(
    seq: &[u8],
    strand: Strand,
    shift_steps: usize,
) -> Vec<u8> {
    let mut rotated: Vec<u8> = seq.iter().map(|base| base.to_ascii_uppercase()).collect();
    if rotated.is_empty() {
        return rotated;
    }

    let normalized_steps = shift_steps % rotated.len();
    if normalized_steps == 0 {
        return rotated;
    }

    match strand {
        Strand::Forward => rotated.rotate_left(normalized_steps),
        Strand::Reverse => rotated.rotate_right(normalized_steps),
    }
    rotated
}

/// `hgvs_protein` (`coding::perl_hgvs_protein`) behind the `ENSP...:p.` prefix.
///
/// Returns `None` when the transcript has no coding model or protein id, the
/// alternate allele is no DNA (`vep_core::codon::seq_is_dna`), the row is the
/// reference allele's, or `coding::perl_hgvs_protein` has no notation for the
/// span. An insertion or deletion is first shifted to its most 3' position on
/// the transcript strand
/// (Ensembl VEP's `_return_3prime(1)`, applied for HGVS whatever `--shift_3prime`
/// says) when a reference FASTA is available; the shifted alleles are the
/// deleted reference bases at the shifted position and the rotated insertion.
/// An insertion whose run of placements is longer than that walk is read at
/// the run's last placement, the span its HGVSc is described at
/// (`hgvsc_shift`), where `_return_3prime(1)` stops short of it
/// (TranscriptVariationAllele.pm 303-306).
pub fn generate_hgvsp(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<String> {
    generate_hgvs(variant, transcript, reference_fasta, HgvsOptions::default()).hgvsp
}

/// `generate_hgvsp` read at `described_at`, the span the HGVSc of the same
/// pair is described at (`hgvsc_shift`); `None` reads the allele where the
/// input placed it. `hgvs_protein` reads the
/// allele as `feature_seq`, on the transcript's strand
/// (VariationFeatureOverlapAllele.pm 253-256), so an allele given on the reverse
/// strand is taken on the forward strand. The reference name carries the
/// translation's version unless `remove_hgvsp_version` (TranscriptVariationAllele.pm
/// 1759-1760), and `hgvsp_use_prediction` wraps the description after `p.` in
/// parentheses (1908, 2040); the `=` of a synonymous description is escaped by the
/// output layer afterwards, so it prints as `p.(Ala2%3D)`. `read_where_written`
/// is `Hgvsc::read_where_written` of the same pair, whose `described_at` is
/// then `None`.
fn hgvsp_at(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
    described_at: Option<(u64, u64)>,
    read_where_written: bool,
    options: HgvsOptions,
) -> Option<String> {
    let input_strand = variant.strand;
    let forward: Cow<InputVariant> = match input_strand {
        Strand::Forward => Cow::Borrowed(variant),
        Strand::Reverse => Cow::Owned(variant.on_forward_strand()),
    };
    let variant = forward.as_ref();
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let has_coding_model = transcript.translation.is_some()
        && mapper.cdna_coding_start > 0
        && mapper.cdna_coding_end >= mapper.cdna_coding_start;
    if !has_coding_model {
        return None;
    }

    // An alternate allele that is no DNA has no notation; one that holds
    // unknown bases has the notation `coding::hgvs_protein` reads where its
    // known bases settle a stop codon before the first of them.
    let alt_allele = variant.alt_allele();
    if !alt_allele.is_empty() && !vep_core::codon::seq_is_dna(alt_allele) {
        return None;
    }
    // `hgvs_protein` returns nothing for the reference allele's object
    // (TranscriptVariationAllele.pm:1687); an alternate allele equal to the
    // allele the row is read against takes the synonymous notation.
    if variant.reference_allele_row {
        return None;
    }

    let protein_id = transcript.protein_id.as_ref()?;
    let version = if options.remove_hgvsp_version {
        None
    } else {
        transcript.translation.as_ref().and_then(|t| t.version)
    };
    let ref_name = versioned_reference_name(protein_id, version);

    let shifted_variant = described_at
        .and_then(|span| shifted_variant(variant, input_strand, transcript, span, reference_fasta));
    let notation = if read_where_written {
        crate::coding::hgvs_protein_where_written(variant, transcript, reference_fasta)
    } else {
        crate::coding::perl_hgvs_protein(
            variant,
            shifted_variant.as_ref().map(|v| (v, v.start, v.end)),
            transcript,
            reference_fasta,
            options.shift_3prime,
        )
    }?;
    Some(if options.hgvsp_use_prediction {
        format!("{ref_name}:p.({notation})")
    } else {
        format!("{ref_name}:p.{notation}")
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::{make_single_exon_transcript, make_test_transcript};

    const SHIFTED: HgvsOptions = HgvsOptions {
        shift_hgvs: true,
        hgvsp_use_prediction: false,
        remove_hgvsp_version: false,
        shift_3prime: false,
    };
    const UNSHIFTED: HgvsOptions = HgvsOptions {
        shift_hgvs: false,
        ..SHIFTED
    };
    /// `--shift_3prime 1`: the consequences read the shifted span.
    const SHIFT_3PRIME: HgvsOptions = HgvsOptions {
        shift_3prime: true,
        ..SHIFTED
    };

    /// The count along the genome's exons of a model the genome is the
    /// reference of, for the two shapes of edit such a model keeps. A forward
    /// model of exons 101-200 and 301-500 whose cached pairs hold one position
    /// the genome lacks, cDNA 151: a position before it keeps its number, the
    /// extra position stands for the base after it, and a position past it,
    /// the coding end with it, counts one lower. A reverse model of the same
    /// exons whose cached pairs pass over the genomic base 450: a position
    /// past it, the coding end with it, counts one higher.
    #[test]
    fn a_model_the_genome_is_the_reference_of_is_counted_along_its_exons() {
        use vep_core::transcript::{Exon, MapperPair};
        let exon = |start, end, rank| Exon {
            start,
            end,
            rank,
            phase: -1,
            end_phase: -1,
        };
        let pair = |from_start, from_end, to_start, to_end, ori| MapperPair {
            from_start,
            from_end,
            to_start,
            to_end,
            ori,
        };
        let mut model = make_test_transcript();
        model.exons = vec![exon(101, 200, 1), exon(301, 500, 2)];

        let cached = ExonSpan::in_genomic_order(&[
            pair(1, 100, 101, 200, 1),
            pair(101, 150, 301, 350, 1),
            pair(152, 301, 351, 500, 1),
        ]);
        let count = GenomeCount::of(&model, &cached, Some((21, 261)), false);
        let numbers = |cdna: &[i64], reverse: bool, count: &GenomeCount, cached: &[ExonSpan]| {
            cdna.iter()
                .map(|&position| count.number(cached, position, reverse))
                .collect::<Vec<_>>()
        };
        assert_eq!(
            numbers(&[1, 100, 101, 150, 151, 152, 301], false, &count, &cached),
            [1, 100, 101, 150, 151, 151, 300]
        );
        assert_eq!(count.coding, Some((21, 260)));
        let counted: Vec<_> = count
            .exons
            .iter()
            .map(|e| (e.start, e.end, e.cdna_start, e.cdna_end))
            .collect();
        assert_eq!(counted, [(101, 200, 1, 100), (301, 500, 101, 300)]);

        model.exons = vec![exon(301, 500, 1), exon(101, 200, 2)];
        let cached = ExonSpan::in_genomic_order(&[
            pair(1, 50, 451, 500, -1),
            pair(51, 199, 301, 449, -1),
            pair(200, 299, 101, 200, -1),
        ]);
        let count = GenomeCount::of(&model, &cached, Some((11, 250)), true);
        assert_eq!(
            numbers(&[1, 50, 51, 199, 200, 299], true, &count, &cached),
            [1, 50, 52, 200, 201, 300]
        );
        assert_eq!(count.coding, Some((11, 251)));
        let counted: Vec<_> = count
            .exons
            .iter()
            .map(|e| (e.start, e.end, e.cdna_start, e.cdna_end))
            .collect();
        assert_eq!(counted, [(101, 200, 201, 300), (301, 500, 1, 200)]);
    }

    /// A model the genome is the reference of whose record holds no exon list
    /// has no exons to count along and is numbered on the cached pairs: a
    /// deletion in the 5' UTR and a substitution in the 3' UTR keep the
    /// notations the model takes with its exon list.
    #[test]
    fn a_failed_edit_model_without_an_exon_list_is_numbered_on_the_cached_pairs() {
        let mut model = make_test_transcript();
        model.bam_edit_status = Some(BamEditStatus::Failed);
        let records = [
            InputVariant::new(
                "21".into(),
                25_000_010,
                25_000_010,
                b"A".to_vec(),
                b"-".to_vec(),
            ),
            InputVariant::new(
                "21".into(),
                25_005_000,
                25_005_000,
                b"A".to_vec(),
                b"C".to_vec(),
            ),
        ];
        let notations = |model: &Transcript| {
            records
                .each_ref()
                .map(|record| generate_hgvsc(record, model, None))
        };
        let with_the_list = notations(&model);
        assert_eq!(
            with_the_list,
            [
                Some("ENST00000000001.1:c.-40del".to_string()),
                Some("ENST00000000001.1:c.*701A>C".to_string()),
            ]
        );
        model.exons.clear();
        model.derived = Default::default();
        assert_eq!(notations(&model), with_the_list);
    }

    /// `shift_feature_seqs` on the inserted allele `CAG` (forward strand), as
    /// it is left for the codons and HGVSp: read on the transcript's strand the
    /// allele turns one base per step of the shift, so along a forward
    /// transcript the forward-written allele turns left and along a reverse
    /// transcript the other way, a shift of the allele's length or of a
    /// multiple leaves it as written, a shift of four turns a three-base allele
    /// as a shift of one does, and the reverse-strand spelling `CTG` of the same
    /// insertion gives the same forward-strand bases. Ensembl VEP 116.3 counts
    /// the turns from the transcript's strand alone, `seq_length -
    /// shift_length` of them on a reverse-strand transcript
    /// (TranscriptVariationAllele.pm:1337-1339), so it turns nothing once the
    /// shift exceeds the allele (`21 46302557 46302556 -/CAG 1` on
    /// ENST00000291691 shifts 4: the turned `GCA`, `TGC` on the transcript,
    /// spells the repeat's own `p.Leu248dup` where the unturned `CAG` spells
    /// `p.Gln249delinsProGlu`, a protein the variant never produces) and turns
    /// a reverse-written allele the wrong way (`-/CTT -1` on APP,
    /// ENST00000346798, reverse, gives `p.Ala242_Asp243insLeu` where its twin
    /// `-/AAG 1` gives `insSer`).
    #[test]
    fn shift_feature_seq_turns_the_written_allele_by_the_transcript_strand_rule() {
        use Strand::{Forward as F, Reverse as R};
        let turned = |input, transcript, shift| {
            String::from_utf8(shift_feature_seq(b"CAG", input, transcript, shift)).unwrap()
        };
        for (input, transcript, shift, expected) in [
            (F, F, 1, "AGC"),
            (F, F, 2, "GCA"),
            (F, F, 4, "AGC"),
            (F, R, 1, "GCA"),
            (F, R, 3, "CAG"),
            (F, R, 4, "GCA"),
            (F, R, 6, "CAG"),
            (R, F, 1, "AGC"),
            (R, F, 4, "AGC"),
            (R, R, 1, "GCA"),
            (R, R, 2, "AGC"),
            (R, R, 4, "GCA"),
        ] {
            assert_eq!(
                turned(input, transcript, shift),
                expected,
                "input {input:?}, transcript {transcript:?}, shift {shift}"
            );
        }
        assert_eq!(shift_feature_seq(b"", F, F, 2), b"");
    }

    /// The shifted allele of an insertion follows [`shift_feature_seq`] (`CAG`
    /// moved four bases turns by the remainder, one base, read on the
    /// transcript's strand: the other way along a reverse-strand transcript,
    /// and along a forward-strand one when the record wrote it on the reverse
    /// strand); a deletion's is `-`, an unmoved insertion's `None`.
    #[test]
    fn shifted_feature_alt_reads_the_input_strand() {
        let insertion = InputVariant::new("21".into(), 101, 100, b"-".to_vec(), b"CAG".to_vec());
        let moved = Some((105u64, 104u64));
        assert_eq!(
            shifted_feature_alt(&insertion, Strand::Forward, Strand::Reverse, moved),
            Some(b"GCA".to_vec())
        );
        assert_eq!(
            shifted_feature_alt(&insertion, Strand::Reverse, Strand::Forward, moved),
            Some(b"AGC".to_vec())
        );
        assert_eq!(
            shifted_feature_alt(
                &insertion,
                Strand::Forward,
                Strand::Forward,
                Some((101, 100))
            ),
            None
        );
        let deletion = InputVariant::new("21".into(), 101, 103, b"CAG".to_vec(), b"-".to_vec());
        assert_eq!(
            shifted_feature_alt(
                &deletion,
                Strand::Forward,
                Strand::Forward,
                Some((105, 107))
            ),
            None
        );
    }

    #[test]
    fn test_amino_acid_three_letter() {
        assert_eq!(amino_acid_three_letter(b'A'), "Ala");
        assert_eq!(amino_acid_three_letter(b'M'), "Met");
        assert_eq!(amino_acid_three_letter(b'*'), "Ter");
        assert_eq!(amino_acid_three_letter(b'X'), "Xaa");
        assert_eq!(amino_acid_three_letter(b'B'), "Asx");
        assert_eq!(amino_acid_three_letter(b'Z'), "Glx");
        assert_eq!(amino_acid_three_letter(b'J'), "Xle");
        assert_eq!(amino_acid_three_letter(b'O'), "Pyl");
    }

    /// The fixture transcript as `_get_cDNA_position` walks it: a forward-strand
    /// slice with no FASTA and the mapper's exon pairs.
    fn fixture_slice_and_exons(tx: &Transcript) -> (TranscriptSlice<'static>, Vec<ExonSpan>) {
        let mapper = tx.vefc.as_ref().unwrap().mapper.as_ref().unwrap();
        let exons = ExonSpan::in_genomic_order(&mapper.exon_coord_mapper.pairs);
        let slice = TranscriptSlice {
            chr: "21",
            tr_start: tx.start as i64,
            tr_end: tx.end as i64,
            reverse: false,
            fasta: None,
        };
        (slice, exons)
    }

    #[test]
    fn test_perl_cdna_position_coding() {
        let tx = make_test_transcript();
        let (slice, exons) = fixture_slice_and_exons(&tx);
        let coding = Some((51, 900));
        // cDNA 51 is CDS 1; cDNA 60 is CDS 10 (slice position = cDNA in exon 1).
        assert_eq!(
            perl_cdna_position(&slice, &exons, coding, 51).as_deref(),
            Some("1")
        );
        assert_eq!(
            perl_cdna_position(&slice, &exons, coding, 60).as_deref(),
            Some("10")
        );
    }

    #[test]
    fn test_perl_cdna_position_five_prime_utr() {
        let tx = make_test_transcript();
        let (slice, exons) = fixture_slice_and_exons(&tx);
        let coding = Some((51, 900));
        // The base before the start codon is c.-1; the transcript's first base c.-50.
        assert_eq!(
            perl_cdna_position(&slice, &exons, coding, 50).as_deref(),
            Some("-1")
        );
        assert_eq!(
            perl_cdna_position(&slice, &exons, coding, 1).as_deref(),
            Some("-50")
        );
    }

    #[test]
    fn test_perl_cdna_position_three_prime_utr() {
        let tx = make_test_transcript();
        let (slice, exons) = fixture_slice_and_exons(&tx);
        let coding = Some((51, 900));
        // cDNA 901 (slice 25_004_300 - 25_000_000 + 1 = 4_301) is the first base
        // after the stop codon, c.*1; ten bases on is c.*10.
        assert_eq!(
            perl_cdna_position(&slice, &exons, coding, 4_301).as_deref(),
            Some("*1")
        );
        assert_eq!(
            perl_cdna_position(&slice, &exons, coding, 4_310).as_deref(),
            Some("*10")
        );
    }

    #[test]
    fn test_perl_cdna_position_intronic_offsets() {
        let tx = make_test_transcript();
        let (slice, exons) = fixture_slice_and_exons(&tx);
        let coding = Some((51, 900));
        // Intron 1 spans genomic 25_000_300..25_001_999 between cDNA 300 and 301:
        // two bases into it is c.250+2, two bases before exon 2 is c.251-2, and
        // the midpoint (850 bases from either exon) goes to the upstream exon.
        assert_eq!(
            perl_cdna_position(&slice, &exons, coding, 302).as_deref(),
            Some("250+2")
        );
        assert_eq!(
            perl_cdna_position(&slice, &exons, coding, 1_999).as_deref(),
            Some("251-2")
        );
        assert_eq!(
            perl_cdna_position(&slice, &exons, coding, 1_150).as_deref(),
            Some("250+850")
        );
        assert_eq!(
            perl_cdna_position(&slice, &exons, coding, 1_151).as_deref(),
            Some("251-850")
        );
    }

    #[test]
    fn test_perl_cdna_position_non_coding() {
        let tx = make_test_transcript();
        let (slice, exons) = fixture_slice_and_exons(&tx);
        // Without a CDS the coordinate is the raw cDNA position.
        assert_eq!(
            perl_cdna_position(&slice, &exons, None, 42).as_deref(),
            Some("42")
        );
        assert_eq!(
            perl_cdna_position(&slice, &exons, None, 302).as_deref(),
            Some("300+2")
        );
    }

    #[test]
    fn test_hgvs_coordinate_parts() {
        assert_eq!(hgvs_coordinate_parts("250"), (250, 0));
        assert_eq!(hgvs_coordinate_parts("-50"), (-50, 0));
        assert_eq!(hgvs_coordinate_parts("250+2"), (250, 2));
        assert_eq!(hgvs_coordinate_parts("251-2"), (251, -2));
        assert_eq!(hgvs_coordinate_parts("*10"), (10, 0));
        assert_eq!(hgvs_coordinate_parts("*3"), (3, 0));
        assert_eq!(hgvs_coordinate_parts("-12+7"), (-12, 7));
    }

    #[test]
    fn test_hgvsc_snv_coding() {
        let tx = make_test_transcript();
        // CDS pos 5 = cDNA 55, genomic 25_000_054 (GCT -> GAT = Ala->Asp)
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        assert_eq!(
            generate_hgvsc(&variant, &tx, None).as_deref(),
            Some("ENST00000000001.1:c.5C>A")
        );
    }

    #[test]
    fn test_hgvsc_snv_five_prime_utr() {
        let tx = make_test_transcript();
        // cDNA 10 = 5' UTR, genomic 25_000_009
        let variant = InputVariant::new(
            "21".into(),
            25_000_009,
            25_000_009,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        assert_eq!(
            generate_hgvsc(&variant, &tx, None).as_deref(),
            Some("ENST00000000001.1:c.-41A>G")
        );
    }

    #[test]
    fn test_hgvsc_snv_three_prime_utr() {
        let tx = make_test_transcript();
        // cDNA 901 = first base of 3' UTR
        // genomic = 25_004_000 + (901 - 601) = 25_004_300
        let variant = InputVariant::new(
            "21".into(),
            25_004_300,
            25_004_300,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        assert_eq!(
            generate_hgvsc(&variant, &tx, None).as_deref(),
            Some("ENST00000000001.1:c.*1A>G")
        );
    }

    #[test]
    fn test_hgvsc_intronic() {
        let tx = make_test_transcript();
        // Intron 1: genomic 25_000_300..25_001_999
        // Position near the donor end (close to intron start for forward strand)
        let variant = InputVariant::new(
            "21".into(),
            25_000_302,
            25_000_302,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        assert_eq!(
            generate_hgvsc(&variant, &tx, None).as_deref(),
            Some("ENST00000000001.1:c.250+3A>G")
        );
    }

    #[test]
    fn test_hgvsc_deletion() {
        let tx = make_test_transcript();
        // 3-base deletion at CDS pos 4-6 (genomic 25_000_053..25_000_055)
        let variant = InputVariant::new(
            "21".into(),
            25_000_053,
            25_000_055,
            b"GCT".to_vec(),
            b"-".to_vec(),
        );
        assert_eq!(
            generate_hgvsc(&variant, &tx, None).as_deref(),
            Some("ENST00000000001.1:c.4_6del")
        );
    }

    #[test]
    fn test_hgvsc_non_coding() {
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.translation = None;
        let variant = InputVariant::new(
            "21".into(),
            25_000_050,
            25_000_050,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        assert_eq!(
            generate_hgvsc(&variant, &tx, None).as_deref(),
            Some("ENST00000000001.1:n.51A>G")
        );
    }

    /// A deletion of `start..=end` with a reference allele of the right length,
    /// as the FASTA-less path reads it.
    fn deletion(start: u64, end: u64) -> InputVariant {
        InputVariant::new(
            "21".into(),
            start,
            end,
            vec![b'A'; (end - start + 1) as usize],
            b"-".to_vec(),
        )
    }

    /// `_var2transcript_slice_coords` (TranscriptVariationAllele.pm 2729-2733): a
    /// deletion that begins ten bases before the transcript's first base has its
    /// start clamped onto that base, so it is described from c.-50 (slice
    /// position 1) to c.-41 (slice position 10). Rows of the same shape in
    /// `tests/golden/116/GRCh37-hgvs`: `21:34914375-34916723 -` against
    /// ENST00000290239, `ENST00000290239.6:c.-49_77+1248del`, and
    /// `21:9937919-9937971 -` against ENST00000459307, `ENST00000459307.1:n.48_83del`.
    #[test]
    fn test_hgvsc_deletion_starting_before_transcript_is_clamped_to_first_base() {
        let tx = make_test_transcript();
        assert_eq!(
            generate_hgvsc(&deletion(24_999_990, 25_000_009), &tx, None).as_deref(),
            Some("ENST00000000001.1:c.-50_-41del")
        );
    }

    /// The same rule at the 3' end: a deletion running ten bases past the
    /// transcript's last base (cDNA 2601, c.*1701) has its end clamped onto it,
    /// and starts at cDNA 2592 (c.*1692). Row of the same shape in
    /// `tests/golden/116/GRCh37-hgvs`: `21:34961913-34962053 -` against
    /// ENST00000445393, `ENST00000445393.1:c.*44_*174del`.
    #[test]
    fn test_hgvsc_deletion_running_past_transcript_end_is_clamped_to_last_base() {
        let tx = make_test_transcript();
        assert_eq!(
            generate_hgvsc(&deletion(25_005_991, 25_006_010), &tx, None).as_deref(),
            Some("ENST00000000001.1:c.*1692_*1701del")
        );
    }

    /// Both ends before the transcript, or both past it, give no HGVSc: the
    /// variant is not `within_feature`, and `_var2transcript_slice_coords`
    /// (TranscriptVariationAllele.pm 2729) returns `undef` for the same shape.
    #[test]
    fn test_hgvsc_deletion_wholly_outside_transcript_has_no_notation() {
        let tx = make_test_transcript();
        assert_eq!(
            generate_hgvsc(&deletion(24_999_980, 24_999_999), &tx, None),
            None
        );
        assert_eq!(
            generate_hgvsc(&deletion(25_006_001, 25_006_020), &tx, None),
            None
        );
    }

    /// Both ends outside on opposite sides is not "wholly outside": the whole
    /// transcript is described, c.-50 to c.*1701, when no 3' shift applies.
    #[test]
    fn test_hgvsc_deletion_covering_whole_transcript_is_described_end_to_end() {
        let tx = make_test_transcript();
        assert_eq!(
            generate_hgvsc(&deletion(24_999_990, 25_006_010), &tx, None).as_deref(),
            Some("ENST00000000001.1:c.-50_*1701del")
        );
    }

    /// `hgvs_transcript` keeps the intron offset of each cDNA coordinate as it
    /// read them (`_hgvs_intron_start_offset`, `_hgvs_intron_end_offset`,
    /// TranscriptVariationAllele.pm 1530-1545), 0 for an exonic coordinate, and
    /// keeps neither when it returns `undef`. Exon 1 ends at c.250 (genomic
    /// 25_000_299) and exon 2 opens at c.251 (25_002_000).
    #[test]
    fn intron_offsets_follow_the_two_hgvsc_coordinates() {
        let tx = make_test_transcript();
        let hgvs = |variant: &InputVariant| generate_hgvs(variant, &tx, None, SHIFTED);

        let intronic = InputVariant::new(
            "21".into(),
            25_000_302,
            25_000_302,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let notation = hgvs(&intronic);
        assert_eq!(
            notation.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.250+3A>G")
        );
        assert_eq!(notation.intron_offsets, Some((3, 3)));

        let coding = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let notation = hgvs(&coding);
        assert_eq!(notation.hgvsc.as_deref(), Some("ENST00000000001.1:c.5C>A"));
        assert_eq!(notation.intron_offsets, Some((0, 0)));

        let into_intron = hgvs(&deletion(25_000_298, 25_000_302));
        assert_eq!(
            into_intron.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.249_250+3del")
        );
        assert_eq!(into_intron.intron_offsets, Some((0, 3)));

        let out_of_intron = hgvs(&deletion(25_001_998, 25_002_001));
        assert_eq!(
            out_of_intron.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.251-2_252del")
        );
        assert_eq!(out_of_intron.intron_offsets, Some((-2, 0)));

        let outside = hgvs(&deletion(24_999_980, 24_999_999));
        assert_eq!(outside.hgvsc, None);
        assert_eq!(outside.intron_offsets, None);
    }

    /// An insertion butting against either transcript edge (VEP's start = end + 1
    /// convention, so one of its ends lies inside) is not `within_feature`
    /// (OutputFactory.pm asks for HGVS only then), so the clamp never sees it and
    /// it has no HGVSc.
    #[test]
    fn test_hgvsc_insertion_adjacent_to_transcript_edge_has_no_notation() {
        let tx = make_test_transcript();
        let before = InputVariant::new(
            "21".into(),
            25_000_000,
            24_999_999,
            b"-".to_vec(),
            b"GG".to_vec(),
        );
        let after = InputVariant::new(
            "21".into(),
            25_006_001,
            25_006_000,
            b"-".to_vec(),
            b"GG".to_vec(),
        );
        assert_eq!(generate_hgvsc(&before, &tx, None), None);
        assert_eq!(generate_hgvsc(&after, &tx, None), None);
    }

    /// A substitution given with its unchanged flanks, running over a
    /// transcript's edge with every base inside unchanged, changes no base of the
    /// transcript and has no HGVSc. Ensembl VEP 116.3 reads the whole alternate
    /// allele against the bases inside (TranscriptVariationAllele.pm 1474-1483,
    /// Utils/Sequence.pm 526), strips the shared bases and writes the outside
    /// bases as an insertion at the coordinate `_get_cDNA_position` gives a
    /// position before the first exon in genomic order, intronic to that exon
    /// (2800-2827):
    /// `ENST00000290239.6:c.-49-1_-49insCGGCATGCTG` for
    /// `21 34915340 AGGCATGCTGGGAGCC CGGCATGCTGGGAGCC`, a change ten bases
    /// before the transcript's first base, and
    /// `ENST00000459307.1:n.83_83+1insCAAACATTCAGACCCTCT` for a change past
    /// the last base of a reverse-strand transcript, the genomic-low edge of
    /// either. The forward fixture's first base is c.-50 and the reverse
    /// fixture's span starts at 40,000,001, its first genomic base the
    /// transcript's last; a one-base change inside, c.-48, shows the notation
    /// path live.
    #[test]
    fn test_hgvsc_allele_changing_only_bases_before_the_first_genomic_base_has_no_notation() {
        let tx = make_test_transcript();
        let forward = InputVariant::new(
            "21".into(),
            24_999_990,
            25_000_005,
            b"AAAAAAAAAAAAAAAA".to_vec(),
            b"CAAAAAAAAAAAAAAA".to_vec(),
        );
        assert_eq!(generate_hgvsc(&forward, &tx, None), None);

        let reverse_tx = make_single_exon_transcript(Strand::Reverse, "ATGAAATAA", "AAAAA");
        let reverse = InputVariant::new(
            "21".into(),
            39_999_996,
            40_000_003,
            b"AAAAAAAA".to_vec(),
            b"CAAAAAAA".to_vec(),
        );
        assert_eq!(generate_hgvsc(&reverse, &reverse_tx, None), None);

        let inside = InputVariant::new(
            "21".into(),
            25_000_002,
            25_000_002,
            b"A".to_vec(),
            b"C".to_vec(),
        );
        assert_eq!(
            generate_hgvsc(&inside, &tx, None).as_deref(),
            Some("ENST00000000001.1:c.-48A>C")
        );
    }

    /// A substitution running over a transcript's edge is described by the
    /// bases that change inside the transcript: the alternate allele loses the
    /// bases the clamp cut from the reference, position by position, and
    /// `hgvs_variant_notation` and the clip read the pair inside. Ensembl VEP
    /// 116.3 reads the whole alternate allele against the clamped reference
    /// and writes a `delins` of every alternate base over the bases inside,
    /// `ENST00000290239.6:c.-49_-47delinsTGCTGGGC` for
    /// `21 34915345 TGCTGGGAGCC TGCTGGGCGCC`, whose one changed base is c.-47,
    /// the notation both engines write for `21 34915352 A C`. The forward
    /// fixture's first base is c.-50 and its last c.*1701; the reverse
    /// fixture's last base is genomic 40,000,001, its 3' UTR c.*1 to c.*5. A
    /// pair of two changed bases keeps the `delins` of the two, and a pair that
    /// is its own reverse complement is an inversion, as `hgvs_variant_notation`
    /// (Utils/Sequence.pm 547-553) classes the minimised allele.
    #[test]
    fn test_hgvsc_substitution_running_over_an_edge_describes_the_bases_that_change_inside() {
        let tx = make_test_transcript();
        let over_first_base = |alt: &[u8]| {
            let v = InputVariant::new(
                "21".into(),
                24_999_995,
                25_000_005,
                b"AAAAAAAAAAA".to_vec(),
                alt.to_vec(),
            );
            generate_hgvsc(&v, &tx, None)
        };
        assert_eq!(
            over_first_base(b"AAAAAAACAAA").as_deref(),
            Some("ENST00000000001.1:c.-48A>C")
        );
        assert_eq!(
            over_first_base(b"AAAAAAACGAA").as_deref(),
            Some("ENST00000000001.1:c.-48_-47delinsCG")
        );

        let two_outside_two_inside = |alt: &[u8]| {
            let v = InputVariant::new(
                "21".into(),
                24_999_998,
                25_000_001,
                b"AAAA".to_vec(),
                alt.to_vec(),
            );
            generate_hgvsc(&v, &tx, None)
        };
        assert_eq!(
            two_outside_two_inside(b"CCGC").as_deref(),
            Some("ENST00000000001.1:c.-50_-49delinsGC")
        );
        assert_eq!(
            two_outside_two_inside(b"CCTT").as_deref(),
            Some("ENST00000000001.1:c.-50_-49inv")
        );

        let over_last_base = InputVariant::new(
            "21".into(),
            25_005_995,
            25_006_005,
            b"AAAAAAAAAAA".to_vec(),
            b"AACAAAAAAAA".to_vec(),
        );
        assert_eq!(
            generate_hgvsc(&over_last_base, &tx, None).as_deref(),
            Some("ENST00000000001.1:c.*1698A>C")
        );

        let reverse_tx = make_single_exon_transcript(Strand::Reverse, "ATGAAATAA", "AAAAA");
        let over_reverse_last_base = InputVariant::new(
            "21".into(),
            39_999_996,
            40_000_003,
            b"AAAAAAAA".to_vec(),
            b"AAAAAACA".to_vec(),
        );
        assert_eq!(
            generate_hgvsc(&over_reverse_last_base, &reverse_tx, None).as_deref(),
            Some("ENST00000000001.1:c.*4T>G")
        );
    }

    /// An allele of unequal length straddling a transcript's first base, given
    /// with an unchanged outside base, loses that base from its alternate allele
    /// where the alternate begins with it, so the change inside is described as
    /// the minimised record's would be: an inserted `C` after c.-50 is
    /// `c.-50_-49insC`, a deleted c.-49 is `c.-49del`. Ensembl VEP 116.3 keeps
    /// the outside base in the alternate and writes `insAC` for the insertion
    /// and the substitution `c.-49C>A` for the deletion (the pair `AC` against
    /// `AA` inside). A change straddling the first base has its alternate bases
    /// read inside without the shared bases outside. A deletion that runs
    /// over the first base and has no window wholly inside the transcript has
    /// no notation, the HGVS nomenclature giving it to the genomic description
    /// alone: `CTGA`/`C`, where Ensembl VEP reads the unchanged `C` against
    /// the `A` inside as `c.-50A>C`, and three `A` deleted from a run of five
    /// with two inside, where Ensembl VEP reads `AA` against `AA` and writes
    /// nothing either. A change whose bases lie wholly outside has no notation, where
    /// Ensembl VEP writes the outside bases as an insertion before exon 1.
    #[test]
    fn test_hgvsc_unequal_allele_straddling_the_first_base_keeps_the_change_inside() {
        let tx = make_test_transcript();
        let hgvsc = |start, end, reference: &[u8], alt: &[u8]| {
            let v = InputVariant::new("21".into(), start, end, reference.to_vec(), alt.to_vec());
            generate_hgvsc(&v, &tx, None)
        };
        assert_eq!(
            hgvsc(24_999_999, 25_000_000, b"AA", b"AAC").as_deref(),
            Some("ENST00000000001.1:c.-50_-49insC")
        );
        assert_eq!(
            hgvsc(24_999_999, 25_000_001, b"AAC", b"AA").as_deref(),
            Some("ENST00000000001.1:c.-49del")
        );
        assert_eq!(
            hgvsc(24_999_998, 25_000_001, b"AAAA", b"CC").as_deref(),
            Some("ENST00000000001.1:c.-50_-49delinsCC")
        );
        assert_eq!(hgvsc(24_999_997, 25_000_000, b"AAAA", b"CAAAA"), None);
        // An alternate shorter than the bases the clamp cut: a change the 3'
        // alignment places outside has no notation (Ensembl VEP duplicates the
        // base inside, `c.-49dup` for the `GG` of `21 34915346 GCTGG GGG,GA`).
        assert_eq!(hgvsc(24_999_997, 25_000_000, b"CTGA", b"AA"), None);
        assert_eq!(hgvsc(24_999_996, 25_000_000, b"ACTGA", b"AAA"), None);
        assert_eq!(hgvsc(25_006_000, 25_006_003, b"ATCA", b"AA"), None);
        // A deletion over the first base with no window wholly inside has no
        // notation; any other change straddling it has the bases inside
        // replaced by the change's alternate bases, never by a shared base.
        assert_eq!(hgvsc(24_999_997, 25_000_000, b"CTGA", b"C"), None);
        assert_eq!(hgvsc(24_999_997, 25_000_001, b"AAAAA", b"AA"), None);
        assert_eq!(
            hgvsc(24_999_997, 25_000_000, b"CTGA", b"CGC").as_deref(),
            Some("ENST00000000001.1:c.-50delinsGC")
        );
        assert_eq!(
            hgvsc(24_999_997, 25_000_000, b"CTGA", b"CTC").as_deref(),
            Some("ENST00000000001.1:c.-50A>C")
        );
    }

    /// An allele of unequal length straddling a transcript's last base is
    /// aligned 3' within its span before the outside bases are cut: a base
    /// deleted from or inserted into a run that continues past the last base
    /// lies outside the transcript and has no notation, as both engines write
    /// none for the minimised record (`21 27087234 TCC TC` is
    /// `downstream_gene_variant` on ENST00000312957, whose last base 27087235
    /// is a `C` followed by a `C`, and Ensembl VEP 116.3 writes no HGVSc for
    /// the alleles of `21 27087234 TCC TC,TCCC` either). A deletion reaching
    /// past the last base deletes the bases inside it, `c.*1701del` for two
    /// `C` deleted from a run of three, and one whose deleted base is inside,
    /// the `T` of `TA` or the `C` of `TCA`, is described there. The forward
    /// fixture's last base is c.*1701 at 25,006,000; the reverse fixture's is
    /// c.*5 at 40,000,001, its alleles read as their reverse complements.
    #[test]
    fn test_hgvsc_unequal_allele_straddling_the_last_base_aligns_the_change_3prime() {
        let tx = make_test_transcript();
        let hgvsc = |start, end, reference: &[u8], alt: &[u8]| {
            let v = InputVariant::new("21".into(), start, end, reference.to_vec(), alt.to_vec());
            generate_hgvsc(&v, &tx, None)
        };
        assert_eq!(hgvsc(25_006_000, 25_006_001, b"AA", b"A"), None);
        assert_eq!(hgvsc(25_006_000, 25_006_001, b"AA", b"AAA"), None);
        assert_eq!(hgvsc(25_005_999, 25_006_001, b"CCC", b"CC"), None);
        assert_eq!(
            hgvsc(25_005_999, 25_006_001, b"CCC", b"C").as_deref(),
            Some("ENST00000000001.1:c.*1701del")
        );
        assert_eq!(
            hgvsc(25_006_000, 25_006_001, b"TA", b"A").as_deref(),
            Some("ENST00000000001.1:c.*1701del")
        );
        assert_eq!(
            hgvsc(25_005_999, 25_006_001, b"TCA", b"TA").as_deref(),
            Some("ENST00000000001.1:c.*1701del")
        );

        let reverse_tx = make_single_exon_transcript(Strand::Reverse, "ATGAAATAA", "AAAAA");
        let reverse_hgvsc = |reference: &[u8], alt: &[u8]| {
            let v = InputVariant::new(
                "21".into(),
                40_000_000,
                40_000_001,
                reference.to_vec(),
                alt.to_vec(),
            );
            generate_hgvsc(&v, &reverse_tx, None)
        };
        assert_eq!(reverse_hgvsc(b"GG", b"G"), None);
        assert_eq!(
            reverse_hgvsc(b"GT", b"G").as_deref(),
            Some("ENST00000000001.1:c.*5del")
        );
    }

    /// A deletion written with bases it keeps that runs over a transcript's
    /// first or last base and has no window wholly inside the transcript has
    /// no HGVSc, no offset and no intron offsets: the HGVS nomenclature gives a deletion that extends
    /// beyond the transcribed region to the genomic description alone (DNA
    /// deletion, "deletions extending beyond the transcribed region"). The
    /// compact forward fixture runs from genomic 100, n.1, to 6100, n.2601;
    /// the reverse fixture has n.1 at 6100 and n.2601 at 100. `CTGAC`/`CC` at
    /// 97-101 deletes 98-100, the two bases before 100 and 100 itself, and
    /// 98 differs from 101 and 100 from 97, so the window has no other
    /// placement; `AGTC`/`AC` at 6099-6102 deletes 6100-6101, and 6100
    /// differs from 6102 and 6101 from 6099. Each record
    /// runs over the first base of one fixture and the last base of the
    /// other. Ensembl VEP 116.3 clamps the span onto the transcript
    /// (TranscriptVariationAllele.pm 2732-2733) and reads the whole alternate
    /// allele against the bases inside (1474-1483). The allele is no
    /// insertion or deletion record, so `--shift_3prime 1` reads it the same
    /// way; under `--shift_hgvs 0` no 3' rule is applied and the bases inside
    /// are described as deleted. A deletion with a second window wholly past
    /// the last base has none either: `AGG`/`A` at 6099-6101 in a run of `G`
    /// from 6100 to 6102.
    #[test]
    fn a_kept_allele_deleting_over_a_transcript_end_with_no_window_inside_has_no_notation() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let mut edits = bases_from(94, "GCACTGACG");
        edits.extend(bases_from(6097, "GCAGTCA"));
        let (_dir, fasta) = hgvsc_fasta_fixture(&edits);
        let low = InputVariant::new("1".into(), 97, 101, b"CTGAC".to_vec(), b"CC".to_vec());
        let high = InputVariant::new("1".into(), 6099, 6102, b"AGTC".to_vec(), b"AC".to_vec());
        let forward = compact_non_coding_transcript();
        let reverse = compact_non_coding_reverse_transcript();

        for (variant, tx, unshifted) in [
            (&low, &forward, "ENST00000000001.1:n.1del"),
            (&high, &forward, "ENST00000000001.1:n.2601del"),
            (&high, &reverse, "ENST00000000001.1:n.1del"),
            (&low, &reverse, "ENST00000000001.1:n.2601del"),
        ] {
            for options in [SHIFTED, SHIFT_3PRIME] {
                let notation = generate_hgvs(variant, tx, Some(&fasta), options);
                assert_eq!(notation.hgvsc, None, "{unshifted} {options:?}");
                assert_eq!(notation.offset, None, "{unshifted} {options:?}");
                assert_eq!(notation.intron_offsets, None, "{unshifted} {options:?}");
            }
            let notation = generate_hgvs(variant, tx, Some(&fasta), UNSHIFTED);
            assert_eq!(notation.hgvsc.as_deref(), Some(unshifted));
            assert_eq!(notation.intron_offsets, Some((0, 0)), "{unshifted}");
        }

        let (_dir, run) = hgvsc_fasta_fixture(&bases_from(6098, "CAGGGT"));
        let past = InputVariant::new("1".into(), 6099, 6101, b"AGG".to_vec(), b"A".to_vec());
        assert_eq!(
            generate_hgvs(&past, &forward, Some(&run), SHIFTED).hgvsc,
            None
        );
    }

    /// The deletions beside it keep their notation. One whose deleted bases
    /// lie wholly inside is described there, `n.1_2del` for `GACG`/`GG` at
    /// 99-102; one whose deleted bases lie wholly outside has none. One whose
    /// window runs over the first base within its span and moves wholly
    /// inside along the chromosome, `CAAA`/`CA` at 97-100 in a run of `A` from
    /// 98 to 101, keeps the deletion of the base inside its span, `n.1del`;
    /// one with a window over the last base and another wholly inside,
    /// `GCGT`/`GT` at 6099-6102 with the same base at 6099 and 6101, keeps
    /// the deletion of the base inside its span, `n.2601del`. A deletion
    /// written without kept bases is clamped onto the transcript as Ensembl
    /// VEP clamps it: `n.1_2del` for 98-101 deleted.
    #[test]
    fn a_deletion_beside_a_transcript_end_keeps_its_notation() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_transcript();
        let hgvsc = |fasta: &vep_fasta::IndexedFasta, start, end, reference: &[u8], alt: &[u8]| {
            let v = InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            generate_hgvsc(&v, &tx, Some(fasta))
        };

        let mut edits = bases_from(94, "GCACTGACG");
        edits.extend(bases_from(6097, "GCAGTCA"));
        let (_dir, fasta) = hgvsc_fasta_fixture(&edits);
        assert_eq!(
            hgvsc(&fasta, 99, 102, b"GACG", b"GG").as_deref(),
            Some("ENST00000000001.1:n.1_2del")
        );
        assert_eq!(hgvsc(&fasta, 96, 100, b"ACTGA", b"AGA"), None);
        assert_eq!(
            hgvsc(&fasta, 98, 101, b"TGAC", b"-").as_deref(),
            Some("ENST00000000001.1:n.1_2del")
        );

        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(96, "GCAAAAG"));
        assert_eq!(
            hgvsc(&fasta, 97, 100, b"CAAA", b"CA").as_deref(),
            Some("ENST00000000001.1:n.1del")
        );

        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(6098, "AGCGT"));
        assert_eq!(
            hgvsc(&fasta, 6099, 6102, b"GCGT", b"GT").as_deref(),
            Some("ENST00000000001.1:n.2601del")
        );
    }

    /// The bases an allele adds to a reference it keeps whole, where their
    /// most 3' placement lies directly past a transcript's last base and they
    /// repeat the transcript's last bases, are a duplication of those bases
    /// with intron offsets 0 and 0: the HGVS nomenclature describes a change
    /// that can be a duplication as one (DNA duplication, notes). The compact
    /// forward fixture ends at genomic 6100, n.2601, with 6098-6101 laid as
    /// `GTCA`: `GTC`/`GTCC` at 6098-6100 adds the last base again,
    /// `GTC`/`GTCTC` the last two, and `TCA`/`TCCA` at 6099-6101, a span over
    /// the last base, the last base again. The reverse fixture ends at genomic
    /// 100, where 99-103 are laid as `AGTCA`: `GTC`/`GGTC` and `GTC`/`GTGTC`
    /// at 100-102 add its last base and its last two. Nothing moves, so no
    /// offset is reported. Ensembl VEP 116.3 does not shift an allele of such
    /// a record (TranscriptVariationAllele.pm 120-123) and its clip leaves an
    /// insertion (2204-2284) past the last base, which has a coordinate only
    /// on a reverse-strand transcript, one base off the last (2802,
    /// 2820-2827): an insertion there, nothing on the forward strand.
    #[test]
    fn bases_added_at_the_last_base_that_repeat_it_are_a_duplication_of_the_last_bases() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let mut edits = bases_from(6096, "CAGTCA");
        edits.extend(bases_from(99, "AGTCA"));
        let (_dir, fasta) = hgvsc_fasta_fixture(&edits);
        let forward = compact_non_coding_transcript();
        let reverse = compact_non_coding_reverse_transcript();
        let kept = |start, end, reference: &[u8], alt: &[u8]| {
            InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec())
        };

        for (variant, tx, notation) in [
            (kept(6098, 6100, b"GTC", b"GTCC"), &forward, "n.2601dup"),
            (
                kept(6098, 6100, b"GTC", b"GTCTC"),
                &forward,
                "n.2600_2601dup",
            ),
            (kept(6099, 6101, b"TCA", b"TCCA"), &forward, "n.2601dup"),
            (kept(100, 102, b"GTC", b"GGTC"), &reverse, "n.2601dup"),
            (kept(100, 102, b"GTC", b"GTGTC"), &reverse, "n.2600_2601dup"),
        ] {
            let expected = format!("ENST00000000001.1:{notation}");
            for options in [SHIFTED, SHIFT_3PRIME] {
                let hgvs = generate_hgvs(&variant, tx, Some(&fasta), options);
                assert_eq!(
                    hgvs.hgvsc.as_deref(),
                    Some(expected.as_str()),
                    "{options:?}"
                );
                assert_eq!(hgvs.offset, None, "{notation} {options:?}");
                assert_eq!(hgvs.intron_offsets, Some((0, 0)), "{notation} {options:?}");
            }
            let unshifted = generate_hgvs(&variant, tx, Some(&fasta), UNSHIFTED);
            assert_eq!(unshifted.hgvsc, None, "{notation}");
        }
    }

    /// The same duplication where the record's span ends before the last
    /// base and a run carries the added bases to it: the offset counts the
    /// bases the 3' rule moves them, negative on the reverse strand. Forward,
    /// 6096-6101 laid as `GCCCCA`: `GC`/`GCC` at 6096-6097 adds a `C` that the
    /// run carries three bases to the last base, `n.2601dup`, and `CC`/`CCCC`
    /// at 6097-6098, the span twice, moves two bases to `n.2600_2601dup`.
    /// Reverse, 99-105 laid as `AGGGGCT`: `CT`/`GCT` at 104-105 adds a `C` of
    /// the transcript's strand that moves four bases. Ensembl VEP 116.3
    /// leaves each where its record's span ends, `n.2598_2599insC` and
    /// `n.2598_2599dup` on the forward fixture and `n.2597_2598insC` on the
    /// reverse one, as the notations read under `--shift_hgvs 0`.
    #[test]
    fn bases_a_run_carries_to_the_last_base_are_a_duplication_of_the_last_bases() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let mut edits = bases_from(6095, "AGCCCCA");
        edits.extend(bases_from(99, "AGGGGCT"));
        let (_dir, fasta) = hgvsc_fasta_fixture(&edits);
        let forward = compact_non_coding_transcript();
        let reverse = compact_non_coding_reverse_transcript();
        let kept = |start, end, reference: &[u8], alt: &[u8]| {
            InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec())
        };

        for (variant, tx, notation, offset, where_written) in [
            (
                kept(6096, 6097, b"GC", b"GCC"),
                &forward,
                "n.2601dup",
                3,
                "n.2598_2599insC",
            ),
            (
                kept(6097, 6098, b"CC", b"CCCC"),
                &forward,
                "n.2600_2601dup",
                2,
                "n.2598_2599dup",
            ),
            (
                kept(104, 105, b"CT", b"GCT"),
                &reverse,
                "n.2601dup",
                -4,
                "n.2597_2598insC",
            ),
        ] {
            let expected = format!("ENST00000000001.1:{notation}");
            for options in [SHIFTED, SHIFT_3PRIME] {
                let hgvs = generate_hgvs(&variant, tx, Some(&fasta), options);
                assert_eq!(
                    hgvs.hgvsc.as_deref(),
                    Some(expected.as_str()),
                    "{options:?}"
                );
                assert_eq!(hgvs.offset, Some(offset), "{notation} {options:?}");
                assert_eq!(hgvs.intron_offsets, Some((0, 0)), "{notation} {options:?}");
            }
            let unshifted = generate_hgvs(&variant, tx, Some(&fasta), UNSHIFTED);
            let expected = format!("ENST00000000001.1:{where_written}");
            assert_eq!(unshifted.hgvsc.as_deref(), Some(expected.as_str()));
            assert_eq!(unshifted.offset, None, "{where_written}");
        }
    }

    /// Added bases beside the last base that are no such duplication keep
    /// their row. On the layout `GTCA` at 6098-6101: a base that repeats
    /// nothing, `GTC`/`GTCG`, has one placement alone, directly past the last
    /// base, and no notation; two bases of which only the second repeats the
    /// last base, `GTC`/`GTCGC`, have no whole copy before them and no
    /// notation. Where the base past the last base continues the run, `GTCC`
    /// at 6098-6101, the placements of `GTC`/`GTCC` go on past the last base
    /// and it has no notation. An allele whose placements all lie inside the
    /// transcript, `GC`/`GCC` at 5000-5001 before a second `C`, keeps the
    /// insertion where its span ends, `n.1502_1503insC`, as Ensembl VEP 116.3
    /// writes it.
    #[test]
    fn added_bases_beside_the_last_base_that_are_no_duplication_keep_their_row() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_transcript();
        let hgvs = |fasta: &vep_fasta::IndexedFasta, start, end, reference: &[u8], alt: &[u8]| {
            let v = InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            generate_hgvs(&v, &tx, Some(fasta), SHIFTED)
        };

        let mut edits = bases_from(6096, "CAGTCA");
        edits.extend(bases_from(5000, "GCCA"));
        let (_dir, fasta) = hgvsc_fasta_fixture(&edits);
        assert_eq!(hgvs(&fasta, 6098, 6100, b"GTC", b"GTCG").hgvsc, None);
        assert_eq!(hgvs(&fasta, 6098, 6100, b"GTC", b"GTCGC").hgvsc, None);
        let inside = hgvs(&fasta, 5000, 5001, b"GC", b"GCC");
        assert_eq!(
            inside.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.1502_1503insC")
        );
        assert_eq!(inside.offset, None);

        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(6096, "CAGTCCA"));
        assert_eq!(hgvs(&fasta, 6098, 6100, b"GTC", b"GTCC").hgvsc, None);
    }

    /// The bases an allele adds to a reference it keeps whole, where a run of
    /// placements that starts at or before a transcript's first base carries
    /// them to a point between two of its bases, are described at that most
    /// 3' point (HGVS nomenclature, general recommendations, the 3' rule): a
    /// duplication where the copied bases all lie in the transcript, an
    /// insertion where the added bases copy nothing. The offset counts the
    /// bases the rule moves them. Forward, n.1 at genomic 100: on `CA` four
    /// times at 97-104, `CACA`/`CACACA` at 99-102 adds a `CA` that the run
    /// carries two bases to `n.4_5dup`; on a run of `A` at 97-102,
    /// `AA`/`AAA` at 100-101, a span inside, moves one base to `n.3dup`; on
    /// `GTACGTA` at 97-103, `AC`/`ACGTCC` at 99-100 adds four bases whose
    /// placements run from the point before n.1 to the point after n.3 and
    /// copy nothing, `n.3_4insCCGT`. Reverse, n.1 at genomic 6100: the first
    /// layout's reverse complement at 6095-6104 gives `n.4_5dup` with offset
    /// -2. Ensembl VEP 116.3 leaves each allele where its record's span ends
    /// (TranscriptVariationAllele.pm 120-123, 2204-2284), as the notations
    /// read under `--shift_hgvs 0`.
    #[test]
    fn bases_a_run_carries_across_the_first_base_are_described_at_their_most_3prime_point() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let forward = compact_non_coding_transcript();
        let reverse = compact_non_coding_reverse_transcript();
        let kept = |start, end, reference: &[u8], alt: &[u8]| {
            InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec())
        };
        let mut repeats = bases_from(95, "GGCACACACAG");
        repeats.extend(bases_from(6095, "CTGTGTGTGC"));

        for (edits, variant, tx, notation, offset, where_written) in [
            (
                repeats.clone(),
                kept(99, 102, b"CACA", b"CACACA"),
                &forward,
                "n.4_5dup",
                2,
                "n.3_4insCA",
            ),
            (
                bases_from(95, "GCAAAAAAC"),
                kept(100, 101, b"AA", b"AAA"),
                &forward,
                "n.3dup",
                1,
                "n.2_3insA",
            ),
            (
                bases_from(97, "GTACGTA"),
                kept(99, 100, b"AC", b"ACGTCC"),
                &forward,
                "n.3_4insCCGT",
                2,
                "n.1_2insGTCC",
            ),
            (
                repeats.clone(),
                kept(6098, 6101, b"TGTG", b"TGTGTG"),
                &reverse,
                "n.4_5dup",
                -2,
                "n.3_4insCA",
            ),
        ] {
            let (_dir, fasta) = hgvsc_fasta_fixture(&edits);
            let expected = format!("ENST00000000001.1:{notation}");
            for options in [SHIFTED, SHIFT_3PRIME] {
                let hgvs = generate_hgvs(&variant, tx, Some(&fasta), options);
                assert_eq!(
                    hgvs.hgvsc.as_deref(),
                    Some(expected.as_str()),
                    "{options:?}"
                );
                assert_eq!(hgvs.offset, Some(offset), "{notation} {options:?}");
                assert_eq!(hgvs.intron_offsets, Some((0, 0)), "{notation} {options:?}");
            }
            let unshifted = generate_hgvs(&variant, tx, Some(&fasta), UNSHIFTED);
            let expected = format!("ENST00000000001.1:{where_written}");
            assert_eq!(unshifted.hgvsc.as_deref(), Some(expected.as_str()));
            assert_eq!(unshifted.offset, None, "{where_written}");
        }
    }

    /// Added bases beside the first base that the rule does not decide keep
    /// their row. On `GCATG` at 97-101, `AT`/`ATCAT` at 99-100 adds three
    /// bases that copy the base before n.1, the base before that and n.1: the
    /// copy starts before the transcript, and the insertion stays where the
    /// span ends, `n.1_2insCAT`. On `GTACGTA` at 97-103, `ACGT`/`ACGTCCGT` at
    /// 99-102 writes its four added bases at their most 3' point already and
    /// keeps `n.3_4insCCGT` with no offset.
    #[test]
    fn added_bases_beside_the_first_base_that_the_rule_does_not_move_keep_their_row() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_transcript();
        let hgvs = |fasta: &vep_fasta::IndexedFasta, start, end, reference: &[u8], alt: &[u8]| {
            let v = InputVariant::new("1".into(), start, end, reference.to_vec(), alt.to_vec());
            generate_hgvs(&v, &tx, Some(fasta), SHIFTED)
        };

        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(97, "GCATG"));
        let copy_from_outside = hgvs(&fasta, 99, 100, b"AT", b"ATCAT");
        assert_eq!(
            copy_from_outside.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.1_2insCAT")
        );
        assert_eq!(copy_from_outside.offset, None);

        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(97, "GTACGTA"));
        let at_its_point = hgvs(&fasta, 99, 102, b"ACGT", b"ACGTCCGT");
        assert_eq!(
            at_its_point.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.3_4insCCGT")
        );
        assert_eq!(at_its_point.offset, None);
        assert_eq!(at_its_point.intron_offsets, Some((0, 0)));
    }

    /// An insertion is described at its most 3' placement however far the run
    /// of its placements goes: the HGVS 3' rule assigns the change to the most
    /// 3' position possible. `perform_shift` reads 1,000 bases of flank
    /// (TranscriptVariationAllele.pm 427-437) and stops where they end
    /// (303-306), so Ensembl VEP 116.3 describes an insertion whose run is
    /// longer where that roll stops. The compact fixture's reference is `T`
    /// throughout; with an `A` at 299 and at 2500 a `T` inserted after 300,
    /// n.201, has placements up to the point after 2499, the hundredth base of
    /// the intron after n.600: `n.600+100dup` with offset 2199 and the intron
    /// offsets 100 and 100, where the roll of 1,000 bases stops in the intron
    /// before n.301 and reads `n.301-799dup` with offset 1000. A run within
    /// the flank and a deletion are outside this case, and under `--shift_3prime 1`
    /// the row is annotated where `perform_shift` stops and described there.
    /// The bases an allele adds to a reference it keeps are carried as far:
    /// with an `A` at 6101, past the last base, the `T` that `TT`/`TTT` at
    /// 4000-4001 adds moves 2,099 bases and duplicates the last base, n.2601.
    #[test]
    fn an_insertion_is_described_at_its_last_placement_past_the_flank_perform_shift_reads() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let tx = compact_non_coding_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(299, b'A'), (2500, b'A')]);
        let insertion = InputVariant::new("1".into(), 301, 300, b"-".to_vec(), b"T".to_vec());
        let hgvs = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFTED);
        assert_eq!(
            hgvs.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.600+100dup")
        );
        assert_eq!(hgvs.offset, Some(2199));
        assert_eq!(hgvs.intron_offsets, Some((100, 100)));
        assert_eq!(
            generate_hgvsc(&insertion, &tx, Some(&fasta)).as_deref(),
            Some("ENST00000000001.1:n.600+100dup")
        );
        let annotated_shifted = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFT_3PRIME);
        assert_eq!(
            annotated_shifted.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.301-799dup")
        );
        assert_eq!(annotated_shifted.offset, Some(1000));

        let (_dir, short_run) = hgvsc_fasta_fixture(&[(299, b'A'), (350, b'A')]);
        let within = generate_hgvs(&insertion, &tx, Some(&short_run), SHIFTED);
        assert_eq!(within.hgvsc.as_deref(), Some("ENST00000000001.1:n.250dup"));
        assert_eq!(within.offset, Some(49));

        let deletion = InputVariant::new("1".into(), 301, 301, b"T".to_vec(), b"-".to_vec());
        let deleted = generate_hgvs(&deletion, &tx, Some(&fasta), SHIFTED);
        assert_eq!(
            deleted.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.301-799del")
        );
        assert_eq!(deleted.offset, Some(1000));

        let (_dir, to_the_last_base) = hgvsc_fasta_fixture(&[(6101, b'A')]);
        let kept = InputVariant::new("1".into(), 4000, 4001, b"TT".to_vec(), b"TTT".to_vec());
        let carried = generate_hgvs(&kept, &tx, Some(&to_the_last_base), SHIFTED);
        assert_eq!(
            carried.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.2601dup")
        );
        assert_eq!(carried.offset, Some(2099));
        assert_eq!(carried.intron_offsets, Some((0, 0)));
    }

    /// A made-up reference for a walk along a transcript: the run of `CAG`
    /// over the slice positions `run`, a `T` at every other position, and
    /// the positions a walk read from it.
    struct RunSource {
        run: std::ops::RangeInclusive<i64>,
        read: std::cell::RefCell<Vec<i64>>,
    }

    impl RunSource {
        fn over(run: std::ops::RangeInclusive<i64>) -> Self {
            Self {
                run,
                read: Default::default(),
            }
        }

        /// The base the run holds at `pos`, whether or not it reaches there.
        fn run_base(pos: i64) -> u8 {
            b"CAG"[pos.rem_euclid(3) as usize]
        }

        /// The bases the run holds at `first..=last`.
        fn run_bases(first: i64, last: i64) -> Vec<u8> {
            (first..=last).map(Self::run_base).collect()
        }

        fn base_at(&self, pos: i64) -> Option<u8> {
            self.read.borrow_mut().push(pos);
            Some(if self.run.contains(&pos) {
                Self::run_base(pos)
            } else {
                b'T'
            })
        }

        /// The count of the bases read, with the lowest and the highest
        /// position among them; `None` when nothing was read.
        fn reads(&self) -> Option<(usize, i64, i64)> {
            let read = self.read.borrow();
            Some((read.len(), *read.iter().min()?, *read.iter().max()?))
        }
    }

    /// The window a kept allele deletes is carried along a run only as far
    /// as the transcript end it runs over, where the answer is settled, with
    /// at most two bases read for each base the clamp cut. The made-up
    /// reference is one run of `CAG`, 100,000 bases each way around a
    /// transcript of 600 bases. Eight bases of the run with three of them
    /// deleted, written at -6 to 1 over the first base, have a window wholly
    /// inside the transcript two moves on, at 1 to 3; written at 599 to 606
    /// over the last base, they have one six moves back, at 598 to 600: both
    /// keep their row. Where the run stops at slice position 2, or starts at
    /// 599, the window cannot be carried inside, and the deletion runs over
    /// that end however it is written.
    #[test]
    fn a_deleted_window_is_carried_no_farther_than_the_transcript_end_it_runs_over() {
        const TR_LEN: i64 = 600;
        let runs_over_an_end = |source: &RunSource, raw_start: i64, cut: (usize, usize)| {
            let reference = RunSource::run_bases(raw_start, raw_start + 7);
            kept_allele_deletes_over_a_transcript_end(
                TR_LEN,
                |pos| source.base_at(pos),
                &reference,
                &reference[..5],
                raw_start,
                cut.0,
                cut.1,
            )
        };
        let (at_the_first_base, at_the_last_base) = ((-6, (7, 0)), (599, (0, 6)));
        let cases = [
            (-100_000..=100_000, at_the_first_base, false),
            (-100_000..=100_000, at_the_last_base, false),
            (-100_000..=2, at_the_first_base, true),
            (599..=100_000, at_the_last_base, true),
        ]
        .map(|(run, (raw_start, cut), over_the_end)| {
            let source = RunSource::over(run);
            assert_eq!(
                runs_over_an_end(&source, raw_start, cut),
                over_the_end,
                "{:?} at {raw_start}",
                source.run
            );
            (source, raw_start, cut)
        });
        // A walk takes in one base beyond the span a move: after it towards the
        // 3' end, before it towards the 5' end.
        for (source, raw_start, (front_cut, back_cut)) in cases {
            let (count, lowest, highest) = source.reads().expect("the walk reads the run");
            let farthest_before = raw_start - back_cut as i64;
            let farthest_after = raw_start + 7 + front_cut as i64;
            assert!(
                count <= 2 * (front_cut + back_cut)
                    && farthest_before <= lowest
                    && highest <= farthest_after,
                "{count} bases read from {lowest} to {highest} at {raw_start}"
            );
        }
    }

    /// The walks that carry added bases along a run end where the
    /// description is settled: the one towards the 3' end at the base
    /// directly past the transcript's last base, the one towards the 5' end
    /// as many bases before its first base as the allele adds. The made-up
    /// reference is one run of `CAG` from 100,000 bases before a transcript
    /// of 600 bases. Three bases added after slice position 300 that
    /// continue the run have no description where the run goes on 100,000
    /// bases past the last base; where it stops at the last base they
    /// duplicate the last three bases, carried 300 bases; where it stops at
    /// slice position 450 they duplicate 448 to 450, carried 150. A source
    /// that holds no base gives no description at any point.
    #[test]
    fn added_bases_on_a_long_run_are_read_no_farther_than_the_transcript_ends() {
        const TR_LEN: i64 = 600;
        let inserted = RunSource::run_bases(301, 303);
        let described = |source: &RunSource| {
            added_bases_at_a_transcript_end(TR_LEN, |pos| source.base_at(pos), 300, &inserted)
                .map(|(n, carried)| (n.start, n.end, n.kind, n.alt_seq, carried))
        };
        let duplication =
            |start, end, carried| Some((start, end, HgvscKind::Dup, inserted.clone(), carried));
        let cases = [
            (-100_000..=100_000, None),
            (-100_000..=TR_LEN, duplication(598, 600, 300)),
            (-100_000..=450, duplication(448, 450, 150)),
        ]
        .map(|(run, description)| {
            let source = RunSource::over(run);
            assert_eq!(described(&source), description, "{:?}", source.run);
            source
        });
        for source in cases {
            let (count, lowest, highest) = source.reads().expect("the walk reads the run");
            let added = inserted.len() as i64;
            assert!(
                count as i64 <= TR_LEN + added + 1 && -added < lowest && highest <= TR_LEN + 1,
                "{count} bases read from {lowest} to {highest} on {:?}",
                source.run
            );
        }
        for point in [0, 300, TR_LEN] {
            assert!(added_bases_at_a_transcript_end(TR_LEN, |_| None, point, &inserted).is_none());
        }
    }

    /// Inserted bases written directly past a transcript's last base that
    /// repeat its last bases are a duplication of those bases, though the row
    /// is annotated outside the transcript: the description of a change does
    /// not turn on which of its placements the record writes. On the layout
    /// `GTCA` at 6098-6101 a `C` or a `TC` written after 6100, the forward
    /// fixture's n.2601, cannot move and reads `n.2601dup` or
    /// `n.2600_2601dup`; on `AGT` at 99-101 a `G` written before 100, the
    /// reverse fixture's n.2601, is the transcript strand's `C` after its last
    /// base. The pair of intron offsets is 0 and 0 and nothing moved, so no
    /// offset is reported, and a row outside has no HGVSp. Ensembl VEP 116.3
    /// asks for no notation on a row outside (OutputFactory.pm 1730). Under
    /// `--shift_3prime 1` and under `--shift_hgvs 0` the row has no notation.
    #[test]
    fn bases_written_past_the_last_base_that_repeat_it_are_a_duplication_of_the_last_bases() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let mut edits = bases_from(6096, "CAGTCA");
        edits.extend(bases_from(99, "AGT"));
        let (_dir, fasta) = hgvsc_fasta_fixture(&edits);
        let forward = compact_non_coding_transcript();
        let reverse = compact_non_coding_reverse_transcript();
        let insertion = |start, alt: &[u8]| {
            InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec())
        };

        for (variant, tx, notation) in [
            (insertion(6101, b"C"), &forward, "n.2601dup"),
            (insertion(6101, b"TC"), &forward, "n.2600_2601dup"),
            (insertion(100, b"G"), &reverse, "n.2601dup"),
        ] {
            let expected = format!("ENST00000000001.1:{notation}");
            let hgvs = generate_hgvs(&variant, tx, Some(&fasta), SHIFTED);
            assert_eq!(hgvs.hgvsc.as_deref(), Some(expected.as_str()));
            assert_eq!(hgvs.offset, None, "{notation}");
            assert_eq!(hgvs.intron_offsets, Some((0, 0)), "{notation}");
            assert_eq!(hgvs.hgvsp, None, "{notation}");
            for options in [SHIFT_3PRIME, UNSHIFTED] {
                let hgvs = generate_hgvs(&variant, tx, Some(&fasta), options);
                assert_eq!(hgvs.hgvsc, None, "{notation} {options:?}");
            }
        }
    }

    /// Inserted bases written before a transcript's first base that a run
    /// carries to a point between two of its bases and that copy nothing are
    /// an insertion at that point. On `GTACGTA` at 97-103, `CGTC` written
    /// before 100, n.1, gives the same sequence up to the point after n.3,
    /// where it reads `CCGT`: `n.3_4insCCGT` with offset 3.
    #[test]
    fn bases_written_before_the_first_base_that_a_run_carries_inside_are_inserted_there() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(97, "GTACGTA"));
        let insertion = InputVariant::new("1".into(), 100, 99, b"-".to_vec(), b"CGTC".to_vec());
        let hgvs = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFTED);
        assert_eq!(
            hgvs.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.3_4insCCGT")
        );
        assert_eq!(hgvs.offset, Some(3));
        assert_eq!(hgvs.intron_offsets, Some((0, 0)));
        assert_eq!(hgvs.hgvsp, None);
        assert_eq!(
            generate_hgvs(&insertion, &tx, Some(&fasta), UNSHIFTED).hgvsc,
            None
        );
    }

    /// Inserted bases written beside a transcript's end that the rule does not
    /// place on it keep a row with no notation. Past the last base, on `GTCA`
    /// at 6098-6101: a `G`, which repeats nothing; `GC`, of which only the
    /// second base repeats the last base; and, where 6101 continues the run,
    /// a `C` whose placements go on past the last base. Before the first
    /// base: a `G` on `GTACGTA` at 97-103, which cannot move inside, and
    /// `TCA` on `GCATG` at 97-101, which moves one base inside as a copy of
    /// the two bases before n.1 and n.1 itself, a copy that starts before
    /// the transcript.
    #[test]
    fn bases_written_beside_a_transcript_end_that_the_rule_does_not_place_have_no_notation() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_transcript();
        let hgvsc = |edits: Vec<(usize, u8)>, start: u64, alt: &[u8]| {
            let (_dir, fasta) = hgvsc_fasta_fixture(&edits);
            let v = InputVariant::new("1".into(), start, start - 1, b"-".to_vec(), alt.to_vec());
            generate_hgvs(&v, &tx, Some(&fasta), SHIFTED).hgvsc
        };
        assert_eq!(hgvsc(bases_from(6096, "CAGTCA"), 6101, b"G"), None);
        assert_eq!(hgvsc(bases_from(6096, "CAGTCA"), 6101, b"GC"), None);
        assert_eq!(hgvsc(bases_from(6096, "CAGTCCA"), 6101, b"C"), None);
        assert_eq!(hgvsc(bases_from(97, "GTACGTA"), 100, b"G"), None);
        assert_eq!(hgvsc(bases_from(97, "GCATG"), 100, b"TCA"), None);
    }

    /// With a FASTA the reference bases inside come from it and the alternate
    /// is cut the same way: a `C` at genomic 102 of the compact fixture, c.-48,
    /// given inside eleven bases from 95, is `c.-48T>C`, and the same eleven
    /// bases changed at 95 alone have no notation.
    #[test]
    fn test_hgvsc_substitution_over_the_first_base_reads_the_bases_inside_from_the_fasta() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        let tx = make_compact_test_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&[]);
        let hgvsc = |alt: &[u8]| {
            let v = InputVariant::new("1".into(), 95, 105, b"TTTTTTTTTTT".to_vec(), alt.to_vec());
            generate_hgvsc(&v, &tx, Some(&fasta))
        };
        assert_eq!(
            hgvsc(b"TTTTTTTCTTT").as_deref(),
            Some("ENST00000000001.1:c.-48T>C")
        );
        assert_eq!(hgvsc(b"CTTTTTTTTTT"), None);
    }

    /// `--shift_hgvs 0` sets the adaptor flag to 0 (BaseRunner.pm:498-500), so
    /// `hgvs_transcript` skips `_return_3prime(1)` (TranscriptVariationAllele.pm
    /// 1441) and describes the allele where the input placed it with no offset
    /// (1444, 1465; 1711 for the protein). The same fixtures shift under the default.
    ///
    /// Forward strand, CDS position = genomic - 149: an inserted `G` before
    /// genomic `GG` at 211-212 is `c.63dup` with offset 2 when shifted; unshifted,
    /// `hgvs_variant_notation` reads the base after the site first (Utils/Sequence.pm
    /// 573-578), the `G` at 211, so the duplication is `c.62dup`. Deleting `GG` at
    /// 210-211 before another `G` is `c.62_63del` with offset 1 shifted and
    /// `c.61_62del` unshifted. HGVSp follows the same allele, so unshifted it equals
    /// the notation computed with no FASTA at all.
    #[test]
    fn shift_hgvs_0_leaves_an_indel_in_a_repeat_unshifted() {
        use crate::concordance_tests::{hgvsc_fasta_fixture, make_compact_test_transcript};
        let tx = make_compact_test_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(210, b'G'), (211, b'G'), (212, b'G')]);

        let insertion = InputVariant::new("1".into(), 211, 210, b"-".to_vec(), b"G".to_vec());
        let shifted = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFTED);
        assert_eq!(shifted.hgvsc.as_deref(), Some("ENST00000000001.1:c.63dup"));
        assert_eq!(shifted.offset, Some(2));
        let unshifted = generate_hgvs(&insertion, &tx, Some(&fasta), UNSHIFTED);
        assert_eq!(
            unshifted.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.62dup")
        );
        assert_eq!(unshifted.offset, None);

        let deletion = InputVariant::new("1".into(), 210, 211, b"GG".to_vec(), b"-".to_vec());
        let shifted = generate_hgvs(&deletion, &tx, Some(&fasta), SHIFTED);
        assert_eq!(
            shifted.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.62_63del")
        );
        assert_eq!(shifted.offset, Some(1));
        let unshifted = generate_hgvs(&deletion, &tx, Some(&fasta), UNSHIFTED);
        assert_eq!(
            unshifted.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.61_62del")
        );
        assert_eq!(unshifted.offset, None);
        let without_fasta = generate_hgvs(&deletion, &tx, None, SHIFTED);
        assert!(unshifted.hgvsp.is_some());
        assert_eq!(unshifted.hgvsp, without_fasta.hgvsp);
        assert_ne!(unshifted.hgvsp, shifted.hgvsp);
    }

    /// An allele given on the strand opposite the transcript's turns the
    /// `--shift_hgvs 0` duplication lookup round (`$lookup_order = -1` when
    /// `$vf->strand` differs from the transcript's, TranscriptVariationAllele.pm
    /// 1462-1466): the bases before the site in transcript orientation are tried
    /// first (Utils/Sequence.pm 573, 579, 584). An inserted `C` between genomic
    /// 5990 and 5991 inside `CCC` at 5989-5991, read as `G` inside `GGG` on the
    /// reverse-strand transcript (CDS 60-62, CDS position = 6051 - genomic),
    /// shifts two bases to `c.62dup` with offset -2; unshifted it duplicates the
    /// `G` before it, `c.60dup`, where the default order names the `G` after it,
    /// `c.61dup`. The same insertion given on the reverse strand, `-/G -1`, is on
    /// the transcript's own strand and keeps the default order: `c.61dup`
    /// unshifted, `c.62dup` shifted.
    #[test]
    fn shift_hgvs_0_on_the_reverse_strand_looks_a_duplication_up_5prime_first() {
        use crate::concordance_tests::{
            compact_on_chr1, hgvsc_fasta_fixture, make_descending_reverse_strand_transcript,
        };
        let tx = compact_on_chr1(make_descending_reverse_strand_transcript());
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(5989, b'C'), (5990, b'C'), (5991, b'C')]);
        let insertion = InputVariant::new("1".into(), 5991, 5990, b"-".to_vec(), b"C".to_vec());

        let shifted = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFTED);
        assert_eq!(shifted.hgvsc.as_deref(), Some("ENST00000000001.1:c.62dup"));
        assert_eq!(shifted.offset, Some(-2));
        let unshifted = generate_hgvs(&insertion, &tx, Some(&fasta), UNSHIFTED);
        assert_eq!(
            unshifted.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.60dup")
        );
        assert_eq!(unshifted.offset, None);

        let mut given_on_reverse =
            InputVariant::new("1".into(), 5991, 5990, b"-".to_vec(), b"G".to_vec());
        given_on_reverse.strand = Strand::Reverse;
        let shifted = generate_hgvs(&given_on_reverse, &tx, Some(&fasta), SHIFTED);
        assert_eq!(shifted.hgvsc.as_deref(), Some("ENST00000000001.1:c.62dup"));
        assert_eq!(shifted.offset, Some(-2));
        let unshifted = generate_hgvs(&given_on_reverse, &tx, Some(&fasta), UNSHIFTED);
        assert_eq!(
            unshifted.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.61dup")
        );
        assert_eq!(unshifted.offset, None);
    }

    /// The compact forward-strand fixture as a non-coding transcript: n.1 is
    /// its first base, genomic 100, and n.2601 its last, genomic 6100.
    fn compact_non_coding_transcript() -> Transcript {
        use crate::concordance_tests::{make_compact_test_transcript, without_coding_model};
        without_coding_model(make_compact_test_transcript())
    }

    /// The compact reverse-strand fixture as a non-coding transcript: n.1 is
    /// its first base, genomic 6100, and n.2601 its last, genomic 100.
    fn compact_non_coding_reverse_transcript() -> Transcript {
        use crate::concordance_tests::{
            compact_on_chr1, make_descending_reverse_strand_transcript, without_coding_model,
        };
        compact_on_chr1(without_coding_model(
            make_descending_reverse_strand_transcript(),
        ))
    }

    /// An insertion anchored on the base before a forward-strand transcript's
    /// first base is `upstream_gene_variant` as written; a run of the inserted
    /// base carries it into exon 1, and the row takes the notation and offset
    /// of the exonic insertion the 3' rule makes of it, whether it is
    /// annotated there, under `--shift_3prime 1`, or where the record writes
    /// it. The fixture copies GRCh38 21:25928745-25928775,
    /// `TTTGATGCCTTTTTTTTTTTTTTAACTGGAG`, with a non-coding transcript starting
    /// at its tenth base, as ENST00000617755 starts at 21:25928754: n.1-14 are
    /// `T`, n.15 is `A`. A `T` inserted before n.1 moves 14 bases and
    /// duplicates n.14, `TT` duplicates n.13_14, each with offset 14. Under
    /// `--shift_3prime 1` Ensembl
    /// VEP writes these rows `n.15A>T` and `n.15delinsTT`: it clamps the
    /// insertion's slice span (1, 0) onto the first base before adding the
    /// offset (TranscriptVariationAllele.pm 2733, 1477-1478), while for the same
    /// alleles written at the shifted position both engines print `n.14dup` and
    /// `n.13_14dup`.
    #[test]
    fn an_insertion_the_shift_carries_into_exon_1_is_the_duplication_the_sequence_shows() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_transcript();
        assert_eq!(tx.start, 100);
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(91, "TTTGATGCCTTTTTTTTTTTTTTAACTGGAG"));

        for (alt, notation) in [
            ("T", "ENST00000000001.1:n.14dup"),
            ("TT", "ENST00000000001.1:n.13_14dup"),
        ] {
            let insertion = InputVariant::new("1".into(), 100, 99, b"-".to_vec(), alt.into());
            let annotated_shifted = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFT_3PRIME);
            assert_eq!(annotated_shifted.hgvsc.as_deref(), Some(notation), "{alt}");
            assert_eq!(annotated_shifted.offset, Some(14), "{alt}");
            assert_eq!(annotated_shifted.hgvsp, None, "{alt}");

            let annotated_as_written = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFTED);
            assert_eq!(
                annotated_as_written.hgvsc.as_deref(),
                Some(notation),
                "{alt}"
            );
            assert_eq!(annotated_as_written.offset, Some(14), "{alt}");
            assert_eq!(annotated_as_written.hgvsp, None, "{alt}");
        }
    }

    /// The same shape along a longer run: GRCh38 21:16420331-16420365 is
    /// `AATTTCTTTCTTTTTTTTTTTTTTTTTAAATCAGG` and ENST00000698005 starts at its
    /// eleventh base, 21:16420341, so n.1-17 are `T` and n.18 is `A`. A `T`
    /// inserted before n.1 moves 17 bases to `n.17dup` with offset 17, the
    /// notation both engines print for `T/TT` written at 21:16420357.
    #[test]
    fn an_insertion_carried_along_a_longer_run_duplicates_the_runs_last_base() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_transcript();
        let (_dir, fasta) =
            hgvsc_fasta_fixture(&bases_from(90, "AATTTCTTTCTTTTTTTTTTTTTTTTTAAATCAGG"));

        let insertion = InputVariant::new("1".into(), 100, 99, b"-".to_vec(), b"T".to_vec());
        let annotated_shifted = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFT_3PRIME);
        assert_eq!(
            annotated_shifted.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.17dup")
        );
        assert_eq!(annotated_shifted.offset, Some(17));
        assert_eq!(
            generate_hgvs(&insertion, &tx, Some(&fasta), SHIFTED),
            annotated_shifted
        );
    }

    /// The same shape on a reverse-strand transcript, whose first base is its
    /// genomic end: the 31 bases of the first fixture reverse-complemented and
    /// laid so that n.1 sits at the transcript's end, the inserted `T` of the
    /// transcript strand written as the forward-strand `A` between that base
    /// and the next. The shift carries it 14 bases towards lower coordinates
    /// into exon 1: `n.14dup`, offset -14 on the reverse strand, on the row
    /// annotated there and on the row annotated outside, where it is written.
    #[test]
    fn an_insertion_the_shift_carries_into_a_reverse_strand_exon_1_is_described_as_forward() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_reverse_transcript();
        assert_eq!(tx.end, 6100);
        // Reverse complement of `TTTGATGCCTTTTTTTTTTTTTTAACTGGAG`, its 22nd base
        // the transcript's first.
        let (_dir, fasta) =
            hgvsc_fasta_fixture(&bases_from(6079, "CTCCAGTTAAAAAAAAAAAAAAGGCATCAAA"));
        assert_eq!(fasta.base("1", 6100), Some(b'A'));
        assert_eq!(fasta.base("1", 6101), Some(b'G'));

        let insertion = InputVariant::new("1".into(), 6101, 6100, b"-".to_vec(), b"A".to_vec());
        let annotated_shifted = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFT_3PRIME);
        assert_eq!(
            annotated_shifted.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.14dup")
        );
        assert_eq!(annotated_shifted.offset, Some(-14));
        let annotated_as_written = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFTED);
        assert_eq!(annotated_as_written, annotated_shifted);
    }

    /// The first fixture with the run of `T` starting five bases before the
    /// transcript: `TTTGTTTTTTTTTTTTTTTTTTTAACTGGAG` at 91, so 94 is `G`, 95-113
    /// are `T` and n.1-14 are `T`, n.15 `A`, as before. An anchor anywhere along
    /// the run before the first base carries an inserted `T` into exon 1.
    const RUN_FROM_FIVE_BEFORE: &str = "TTTGTTTTTTTTTTTTTTTTTTTAACTGGAG";

    /// The reverse complement of [`RUN_FROM_FIVE_BEFORE`] laid from 6079 for the
    /// reverse-strand transcript ending at 6100: 6086 is `T`, 6087-6105 are `A`,
    /// 6106 `C`, so n.1-14 are `T` and n.15 is `A` on the transcript strand and
    /// the run reaches five bases before the first base.
    const RUN_FROM_FIVE_BEFORE_REVERSE: &str = "CTCCAGTTAAAAAAAAAAAAAAAAAAACAAA";

    /// An insertion anchored two or more bases before a forward-strand
    /// transcript's first base has, as written, a slice span wholly before the
    /// transcript, for which `_var2transcript_slice_coords` returns nothing
    /// (TranscriptVariationAllele.pm 2729) in Ensembl VEP, though the shift
    /// carries it into exon 1 as it carries the one-base anchor. Read at the
    /// shifted span, every anchor along the run gives the duplication the
    /// sequence shows, `n.14dup`, with the offset of its own shift: 14 from the
    /// base before n.1, 15 from two bases before, 19 from the run's first base.
    #[test]
    fn an_insertion_anchored_further_before_exon_1_takes_the_notation_of_its_shifted_span() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(91, RUN_FROM_FIVE_BEFORE));
        assert_eq!(fasta.base("1", 94), Some(b'G'));
        assert_eq!(fasta.base("1", 95), Some(b'T'));
        assert_eq!(fasta.base("1", 113), Some(b'T'));
        assert_eq!(fasta.base("1", 114), Some(b'A'));

        for (anchor, offset) in [(100, 14), (99, 15), (95, 19)] {
            let insertion =
                InputVariant::new("1".into(), anchor, anchor - 1, b"-".to_vec(), b"T".to_vec());
            let annotated_shifted = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFT_3PRIME);
            assert_eq!(
                annotated_shifted.hgvsc.as_deref(),
                Some("ENST00000000001.1:n.14dup"),
                "anchor {anchor}"
            );
            assert_eq!(annotated_shifted.offset, Some(offset), "anchor {anchor}");
            let annotated_as_written = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFTED);
            assert_eq!(annotated_as_written, annotated_shifted, "anchor {anchor}");
        }

        let two_inserted = InputVariant::new("1".into(), 99, 98, b"-".to_vec(), b"TT".to_vec());
        let hgvs = generate_hgvs(&two_inserted, &tx, Some(&fasta), SHIFT_3PRIME);
        assert_eq!(hgvs.hgvsc.as_deref(), Some("ENST00000000001.1:n.13_14dup"));
        assert_eq!(hgvs.offset, Some(15));
    }

    /// The same anchors on the reverse-strand transcript, the inserted `T` of
    /// the transcript strand written as the forward-strand `A` after the anchor
    /// base: `n.14dup` with offsets -14, -15 and -19 towards lower coordinates.
    #[test]
    fn an_insertion_anchored_further_before_a_reverse_strand_exon_1_takes_its_shifted_span() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_reverse_transcript();
        assert_eq!(tx.end, 6100);
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(6079, RUN_FROM_FIVE_BEFORE_REVERSE));
        assert_eq!(fasta.base("1", 6086), Some(b'T'));
        assert_eq!(fasta.base("1", 6087), Some(b'A'));
        assert_eq!(fasta.base("1", 6105), Some(b'A'));
        assert_eq!(fasta.base("1", 6106), Some(b'C'));

        for (anchor, offset) in [(6101, -14), (6102, -15), (6106, -19)] {
            let insertion =
                InputVariant::new("1".into(), anchor, anchor - 1, b"-".to_vec(), b"A".to_vec());
            let annotated_shifted = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFT_3PRIME);
            assert_eq!(
                annotated_shifted.hgvsc.as_deref(),
                Some("ENST00000000001.1:n.14dup"),
                "anchor {anchor}"
            );
            assert_eq!(annotated_shifted.offset, Some(offset), "anchor {anchor}");
            let annotated_as_written = generate_hgvs(&insertion, &tx, Some(&fasta), SHIFTED);
            assert_eq!(annotated_as_written, annotated_shifted, "anchor {anchor}");
        }
    }

    /// A deletion wholly before the first base that the shift carries into
    /// exon 1 is the same shape: nothing as written (2729), the deleted bases at
    /// the shifted span once read there. One `T` of the run deleted at 99 or at
    /// 95 is the run's last base, `n.14del`, moved 14 or 18; `TT` deleted at
    /// 98-99 is `n.13_14del`, moved 14.
    #[test]
    fn a_deletion_the_shift_carries_into_exon_1_is_described_at_its_shifted_span() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(91, RUN_FROM_FIVE_BEFORE));

        for (start, end, deleted, notation, offset) in [
            (99, 99, "T", "ENST00000000001.1:n.14del", 14),
            (95, 95, "T", "ENST00000000001.1:n.14del", 18),
            (98, 99, "TT", "ENST00000000001.1:n.13_14del", 14),
        ] {
            let deletion = InputVariant::new("1".into(), start, end, deleted.into(), b"-".to_vec());
            let annotated_shifted = generate_hgvs(&deletion, &tx, Some(&fasta), SHIFT_3PRIME);
            assert_eq!(
                annotated_shifted.hgvsc.as_deref(),
                Some(notation),
                "{deleted} at {start}"
            );
            assert_eq!(
                annotated_shifted.offset,
                Some(offset),
                "{deleted} at {start}"
            );
            let annotated_as_written = generate_hgvs(&deletion, &tx, Some(&fasta), SHIFTED);
            assert_eq!(annotated_as_written.hgvsc, None, "{deleted} at {start}");
            assert_eq!(annotated_as_written.offset, None, "{deleted} at {start}");
        }
    }

    /// The same deletions on the reverse strand, written as the forward-strand
    /// `A` bases after the transcript's genomic end: `n.14del` moved -14 from
    /// 6101 and -18 from 6105, `n.13_14del` moved -14 from 6101-6102.
    #[test]
    fn a_deletion_the_shift_carries_into_a_reverse_strand_exon_1_is_described_at_its_shifted_span()
    {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_reverse_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(6079, RUN_FROM_FIVE_BEFORE_REVERSE));

        for (start, end, deleted, notation, offset) in [
            (6101, 6101, "A", "ENST00000000001.1:n.14del", -14),
            (6105, 6105, "A", "ENST00000000001.1:n.14del", -18),
            (6101, 6102, "AA", "ENST00000000001.1:n.13_14del", -14),
        ] {
            let deletion = InputVariant::new("1".into(), start, end, deleted.into(), b"-".to_vec());
            let annotated_shifted = generate_hgvs(&deletion, &tx, Some(&fasta), SHIFT_3PRIME);
            assert_eq!(
                annotated_shifted.hgvsc.as_deref(),
                Some(notation),
                "{deleted} at {start}"
            );
            assert_eq!(
                annotated_shifted.offset,
                Some(offset),
                "{deleted} at {start}"
            );
            let annotated_as_written = generate_hgvs(&deletion, &tx, Some(&fasta), SHIFTED);
            assert_eq!(annotated_as_written.hgvsc, None, "{deleted} at {start}");
            assert_eq!(annotated_as_written.offset, None, "{deleted} at {start}");
        }
    }

    /// A deletion straddling the first base overlaps the transcript as written,
    /// so it is described under `--hgvs` alone as well. Ensembl VEP clamps the
    /// written span onto the first base before adding the offset
    /// (TranscriptVariationAllele.pm 2732-2733, 1477-1478): `TT` deleted at
    /// 99-100 is clamped to the one base n.1 and moved 13 to `n.14del`, where the
    /// shifted span covers two bases, `n.13_14del`; `TTT` at 98-100 is
    /// `n.12_14del`, not `n.14del`. Along a run ending at 101 the shifted span
    /// of `TTT` at 98-100 still straddles the first base, 99-101, and the clamp
    /// describes the two bases inside, `n.1_2del` moved 1, where Ensembl VEP
    /// writes `n.2del`. A deletion the shift leaves where it was is clamped the
    /// same way on both engines: `TTTT` at 98-101, the run's last four bases, is
    /// `n.1_2del` with no offset.
    #[test]
    fn a_deletion_straddling_the_first_base_is_described_over_the_bases_the_shift_moved_it_to() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(91, RUN_FROM_FIVE_BEFORE));

        for (start, end, deleted, notation) in [
            (99, 100, "TT", "ENST00000000001.1:n.13_14del"),
            (98, 100, "TTT", "ENST00000000001.1:n.12_14del"),
        ] {
            let deletion = InputVariant::new("1".into(), start, end, deleted.into(), b"-".to_vec());
            for options in [SHIFT_3PRIME, SHIFTED] {
                let hgvs = generate_hgvs(&deletion, &tx, Some(&fasta), options);
                assert_eq!(
                    hgvs.hgvsc.as_deref(),
                    Some(notation),
                    "{deleted} at {start}"
                );
                assert_eq!(hgvs.offset, Some(13), "{deleted} at {start}");
            }
        }

        // `TTTGTTTTTTTAACTGGAG` at 91: the run is 95-101, 102 is `A`.
        let (_dir, short_run) = hgvsc_fasta_fixture(&bases_from(91, "TTTGTTTTTTTAACTGGAG"));
        assert_eq!(short_run.base("1", 101), Some(b'T'));
        assert_eq!(short_run.base("1", 102), Some(b'A'));
        let still_straddling =
            InputVariant::new("1".into(), 98, 100, b"TTT".to_vec(), b"-".to_vec());
        let hgvs = generate_hgvs(&still_straddling, &tx, Some(&short_run), SHIFT_3PRIME);
        assert_eq!(hgvs.hgvsc.as_deref(), Some("ENST00000000001.1:n.1_2del"));
        assert_eq!(hgvs.offset, Some(1));

        let unmoved = InputVariant::new("1".into(), 98, 101, b"TTTT".to_vec(), b"-".to_vec());
        let hgvs = generate_hgvs(&unmoved, &tx, Some(&short_run), SHIFT_3PRIME);
        assert_eq!(hgvs.hgvsc.as_deref(), Some("ENST00000000001.1:n.1_2del"));
        assert_eq!(hgvs.offset, None);
    }

    /// The straddling deletions on the reverse strand, written as forward-strand
    /// `A` bases over the transcript's genomic end 6100: `AA` at 6100-6101 is
    /// `n.13_14del` and `AAA` at 6100-6102 `n.12_14del`, each moved -13; along a
    /// run ending at 6099 (`CTCCAGTTAAAAAAACAAA` at 6091, 6099-6105 `A`) `AAA` at
    /// 6100-6102 moves -1 to 6099-6101 and is `n.1_2del`.
    #[test]
    fn a_deletion_straddling_a_reverse_strand_first_base_is_described_over_its_shifted_span() {
        use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture};
        let tx = compact_non_coding_reverse_transcript();
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(6079, RUN_FROM_FIVE_BEFORE_REVERSE));

        for (start, end, deleted, notation) in [
            (6100, 6101, "AA", "ENST00000000001.1:n.13_14del"),
            (6100, 6102, "AAA", "ENST00000000001.1:n.12_14del"),
        ] {
            let deletion = InputVariant::new("1".into(), start, end, deleted.into(), b"-".to_vec());
            for options in [SHIFT_3PRIME, SHIFTED] {
                let hgvs = generate_hgvs(&deletion, &tx, Some(&fasta), options);
                assert_eq!(
                    hgvs.hgvsc.as_deref(),
                    Some(notation),
                    "{deleted} at {start}"
                );
                assert_eq!(hgvs.offset, Some(-13), "{deleted} at {start}");
            }
        }

        let (_dir, short_run) = hgvsc_fasta_fixture(&bases_from(6091, "CTCCAGTTAAAAAAACAAA"));
        assert_eq!(short_run.base("1", 6098), Some(b'T'));
        assert_eq!(short_run.base("1", 6099), Some(b'A'));
        assert_eq!(short_run.base("1", 6105), Some(b'A'));
        assert_eq!(short_run.base("1", 6106), Some(b'C'));
        let still_straddling =
            InputVariant::new("1".into(), 6100, 6102, b"AAA".to_vec(), b"-".to_vec());
        let hgvs = generate_hgvs(&still_straddling, &tx, Some(&short_run), SHIFT_3PRIME);
        assert_eq!(hgvs.hgvsc.as_deref(), Some("ENST00000000001.1:n.1_2del"));
        assert_eq!(hgvs.offset, Some(-1));
    }

    /// A deletion inside the last exon that the shift carries over the last base
    /// is clamped onto that base like one written there, which both engines
    /// describe over the bases inside; Ensembl VEP refuses the shifted span
    /// instead (TranscriptVariationAllele.pm 1472). On the reference of `T` with
    /// an `A` at 6103, `TTT` deleted at 6097-6099 moves 3 to 6100-6102 and is the
    /// last base alone, `n.2601del`; `TTTT` at 6096-6099 moves 3 to 6099-6102,
    /// `n.2600_2601del`. `TT` at 6098-6099 moves 3 to 6101-6102, wholly past the
    /// transcript, and has no notation (2729).
    #[test]
    fn a_deletion_the_shift_carries_over_the_last_base_is_described_over_the_bases_inside() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let tx = compact_non_coding_transcript();
        assert_eq!(tx.end, 6100);
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(6103, b'A')]);

        for (start, end, deleted, notation) in [
            (6097, 6099, "TTT", "ENST00000000001.1:n.2601del"),
            (6096, 6099, "TTTT", "ENST00000000001.1:n.2600_2601del"),
        ] {
            let deletion = InputVariant::new("1".into(), start, end, deleted.into(), b"-".to_vec());
            for options in [SHIFT_3PRIME, SHIFTED] {
                let hgvs = generate_hgvs(&deletion, &tx, Some(&fasta), options);
                assert_eq!(
                    hgvs.hgvsc.as_deref(),
                    Some(notation),
                    "{deleted} at {start}"
                );
                assert_eq!(hgvs.offset, Some(3), "{deleted} at {start}");
            }
        }

        let carried_past = InputVariant::new("1".into(), 6098, 6099, b"TT".to_vec(), b"-".to_vec());
        for options in [SHIFT_3PRIME, SHIFTED] {
            let hgvs = generate_hgvs(&carried_past, &tx, Some(&fasta), options);
            assert_eq!(hgvs.hgvsc, None);
            assert_eq!(hgvs.offset, None);
        }
    }

    /// The same on the reverse strand, whose last base n.2601 is the genomic
    /// start 100: with an `A` at 97 the run of `T` ends at 98, so `TTT` at
    /// 101-103 moves -3 to 98-100, `n.2601del`, `TTTT` at 101-104 to 98-101,
    /// `n.2600_2601del`, and `TT` at 101-102 to 98-99, wholly past.
    #[test]
    fn a_deletion_the_shift_carries_over_a_reverse_strand_last_base_is_described_inside() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let tx = compact_non_coding_reverse_transcript();
        assert_eq!(tx.start, 100);
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(97, b'A')]);

        for (start, end, deleted, notation) in [
            (101, 103, "TTT", "ENST00000000001.1:n.2601del"),
            (101, 104, "TTTT", "ENST00000000001.1:n.2600_2601del"),
        ] {
            let deletion = InputVariant::new("1".into(), start, end, deleted.into(), b"-".to_vec());
            for options in [SHIFT_3PRIME, SHIFTED] {
                let hgvs = generate_hgvs(&deletion, &tx, Some(&fasta), options);
                assert_eq!(
                    hgvs.hgvsc.as_deref(),
                    Some(notation),
                    "{deleted} at {start}"
                );
                assert_eq!(hgvs.offset, Some(-3), "{deleted} at {start}");
            }
        }

        let carried_past = InputVariant::new("1".into(), 101, 102, b"TT".to_vec(), b"-".to_vec());
        for options in [SHIFT_3PRIME, SHIFTED] {
            let hgvs = generate_hgvs(&carried_past, &tx, Some(&fasta), options);
            assert_eq!(hgvs.hgvsc, None);
            assert_eq!(hgvs.offset, None);
        }
    }

    /// The other edge: an insertion inside the last exon that the shift carries
    /// past the transcript's last base is annotated outside it under
    /// `--shift_3prime 1` (`downstream_gene_variant`), so it has no notation
    /// there, while under `--hgvs` alone the row is exonic and the notation
    /// duplicates the last base, `n.2601dup` with the offset. The same on the
    /// reverse strand, whose last base is the transcript's genomic start.
    #[test]
    fn an_insertion_the_shift_carries_past_the_last_base_has_no_notation_at_the_annotated_span() {
        use crate::concordance_tests::hgvsc_fasta_fixture;
        let forward = compact_non_coding_transcript();
        assert_eq!(forward.end, 6100);
        let (_dir, fasta) = hgvsc_fasta_fixture(&[(6101, b'A'), (99, b'A')]);
        let insertion = InputVariant::new("1".into(), 6095, 6094, b"-".to_vec(), b"T".to_vec());
        let annotated_shifted = generate_hgvs(&insertion, &forward, Some(&fasta), SHIFT_3PRIME);
        assert_eq!(annotated_shifted.hgvsc, None);
        assert_eq!(annotated_shifted.offset, None);
        let annotated_as_written = generate_hgvs(&insertion, &forward, Some(&fasta), SHIFTED);
        assert_eq!(
            annotated_as_written.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.2601dup")
        );
        assert_eq!(annotated_as_written.offset, Some(6));

        let reverse = compact_non_coding_reverse_transcript();
        assert_eq!(reverse.start, 100);
        let insertion = InputVariant::new("1".into(), 106, 105, b"-".to_vec(), b"T".to_vec());
        let annotated_shifted = generate_hgvs(&insertion, &reverse, Some(&fasta), SHIFT_3PRIME);
        assert_eq!(annotated_shifted.hgvsc, None);
        assert_eq!(annotated_shifted.offset, None);
        let annotated_as_written = generate_hgvs(&insertion, &reverse, Some(&fasta), SHIFTED);
        assert_eq!(
            annotated_as_written.hgvsc.as_deref(),
            Some("ENST00000000001.1:n.2601dup")
        );
        assert_eq!(annotated_as_written.offset, Some(-6));
    }

    /// `transcript` as a RefSeq model named `stable_id` whose alignment replaced
    /// one base with one (`_rna_edit` `70 70 T`), the shape of NM_006933.7's
    /// three edits: the mapper is unchanged and `refseq_alignment` carries the
    /// edit's bases.
    fn substitution_edited_model(mut transcript: Transcript, stable_id: &str) -> Transcript {
        transcript.stable_id = stable_id.into();
        transcript.attributes.push(vep_core::transcript::Attribute {
            code: "_rna_edit".into(),
            value: "70 70 T".to_string(),
        });
        transcript.derived = Default::default();
        assert!(transcript
            .refseq_alignment()
            .is_some_and(|alignment| alignment.first_edit_bases.is_some()));
        transcript
    }

    /// On a model with an `_rna_edit`, `hgvs_transcript` respells the alleles
    /// from `feature_seq` (TranscriptVariationAllele.pm 1512-1513). For an
    /// insertion the shift moved, `feature_seq` holds the rotated allele on a
    /// coding row under `--shift_3prime 1` (`codon`, 813-820, has run
    /// `shift_feature_seqs`, 1338-1342) and the allele as written on every other
    /// row, at the shifted coordinates either way. The fixture lays
    /// NM_006933.7's bases around
    /// c.203, `CTGGCA` (GRCh37 21:35467699-35467704), on c.202-c.207 of a
    /// substitution-edited model: `GA` inserted after c.203 shifts one base and
    /// both runs write `c.204_205insAG`, the notation Ensembl VEP 116.3 writes
    /// for `21 35467700 T TGA` on that model under `--refseq --hgvs
    /// --shift_3prime 1`; its run without the flag writes `c.204_205insGA`, the
    /// sequence `CTG GA GC`, which the variant (`CTG AG GC`) does not produce.
    /// The unedited model reads `insAG` on both engines in both runs, and under
    /// `--shift_hgvs 0` the allele is written where the input placed it.
    #[test]
    fn edited_model_hgvsc_spells_a_shifted_insertion_from_the_rotated_allele() {
        use crate::concordance_tests::{
            bases_from, hgvsc_fasta_fixture, make_compact_test_transcript,
        };
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(351, "CTGGCA"));
        let insertion = InputVariant::new("1".into(), 353, 352, b"-".to_vec(), b"GA".to_vec());
        let edited = substitution_edited_model(make_compact_test_transcript(), "NM_000001.1");
        for options in [SHIFTED, SHIFT_3PRIME] {
            let notation = generate_hgvs(&insertion, &edited, Some(&fasta), options);
            assert_eq!(
                notation.hgvsc.as_deref(),
                Some("NM_000001.1:c.204_205insAG"),
                "{options:?}"
            );
            assert_eq!(notation.offset, Some(1), "{options:?}");
        }
        let unedited = make_compact_test_transcript();
        assert_eq!(
            generate_hgvs(&insertion, &unedited, Some(&fasta), SHIFTED)
                .hgvsc
                .as_deref(),
            Some("ENST00000000001.1:c.204_205insAG")
        );
        assert_eq!(
            generate_hgvs(&insertion, &edited, Some(&fasta), UNSHIFTED)
                .hgvsc
                .as_deref(),
            Some("NM_000001.1:c.203_204insGA")
        );
    }

    /// An alternate allele equal to the edited base of a substitution-edited
    /// model gets no HGVSc: with the transcript's reference in force the row's
    /// reference allele is the edited base, and an allele equal to its
    /// reference is no change. The compact model's edit `70 70 T` puts `T` at
    /// c.20 (genomic 169) over a genome `G`: `T/T` writes nothing, `T/A` writes
    /// `c.20T>A`, and the given reference `G/T` writes `c.20G>T`. Ensembl VEP
    /// 116.3 compares the alternate with the genome's base, finds `G` and `T`
    /// different, then respells both alleles from `feature_seq`
    /// (TranscriptVariationAllele.pm 1474-1483, 1508-1514) and writes `c.20T>T`,
    /// the documented class `edited_model_no_change_hgvsc`.
    #[test]
    fn edited_model_hgvsc_is_absent_for_an_alternate_equal_to_the_edited_base() {
        use crate::concordance_tests::{
            bases_from, hgvsc_fasta_fixture, make_compact_test_transcript,
        };
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(169, "G"));
        let edited = substitution_edited_model(make_compact_test_transcript(), "NM_000001.1");
        let hgvsc = |reference: &[u8], alt: &[u8]| {
            let variant = InputVariant::new("1".into(), 169, 169, reference.to_vec(), alt.to_vec());
            generate_hgvs(&variant, &edited, Some(&fasta), SHIFTED).hgvsc
        };
        assert_eq!(hgvsc(b"T", b"T"), None);
        assert_eq!(hgvsc(b"T", b"A").as_deref(), Some("NM_000001.1:c.20T>A"));
        assert_eq!(hgvsc(b"G", b"T").as_deref(), Some("NM_000001.1:c.20G>T"));
    }

    /// On a model with an `_rna_edit`, an allele whose span runs over the
    /// transcript's first or last base is described by the bases that change
    /// inside the transcript, as on a model without one: the pair respelled
    /// relative to the transcript loses the bases the clamp cut from the
    /// span, whatever they are. The forward fixture's first base is c.-50 and
    /// its last c.*1701. `AA` to `CG` over the last base and the base after
    /// it changes c.*1701 from `A` to `C`; `AA` to `GC` over the base before
    /// the first base and the first base changes c.-50 from `A` to `C`; `AAA`
    /// to `CGC` over the last two bases and the base after them puts `CG` in
    /// place of c.*1700 and c.*1701; `AAA` to `ACG` there keeps c.*1700 and
    /// changes c.*1701 from `A` to `C`; `AA` to `G` over the last base and
    /// the base after it puts `G` in place of c.*1701; and `AA` to `CA`,
    /// whose outside base is unchanged, is the substitution of c.*1701 alone.
    /// The reverse fixture's last base, c.*1701, is its genomic start, where
    /// the forward strand's `AA` to `CG` over the base before it and that
    /// base changes the transcript's `T` to `C`. Ensembl VEP 116.3 reads the
    /// whole alternate allele against the bases inside
    /// (TranscriptVariationAllele.pm 1474-1483) and writes both alleles over
    /// whole (1508-1514): `c.*1701delinsCG` for the first record.
    #[test]
    fn edited_model_hgvsc_of_an_allele_over_an_edge_names_the_bases_inside() {
        use crate::concordance_tests::make_descending_reverse_strand_transcript;
        let hgvsc = |model: &Transcript, start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            let variant =
                InputVariant::new("21".into(), start, end, reference.to_vec(), alt.to_vec());
            generate_hgvsc(&variant, model, None).unwrap_or_default()
        };
        let edited = substitution_edited_model(make_test_transcript(), "NM_000001.1");
        let reverse =
            substitution_edited_model(make_descending_reverse_strand_transcript(), "NM_000002.1");
        // The model without the edit takes the same notation.
        let unedited = make_test_transcript();
        assert_eq!(
            [
                hgvsc(&edited, 25_006_000, 25_006_001, b"AA", b"CG"),
                hgvsc(&edited, 24_999_999, 25_000_000, b"AA", b"GC"),
                hgvsc(&edited, 25_005_999, 25_006_001, b"AAA", b"CGC"),
                hgvsc(&edited, 25_005_999, 25_006_001, b"AAA", b"ACG"),
                hgvsc(&edited, 25_006_000, 25_006_001, b"AA", b"G"),
                hgvsc(&edited, 25_006_000, 25_006_001, b"AA", b"CA"),
                hgvsc(&reverse, 24_999_999, 25_000_000, b"AA", b"CG"),
                hgvsc(&unedited, 25_006_000, 25_006_001, b"AA", b"CG"),
            ],
            [
                "NM_000001.1:c.*1701A>C",
                "NM_000001.1:c.-50A>C",
                "NM_000001.1:c.*1700_*1701delinsCG",
                "NM_000001.1:c.*1701A>C",
                "NM_000001.1:c.*1701A>G",
                "NM_000001.1:c.*1701A>C",
                "NM_000002.1:c.*1701T>C",
                "ENST00000000001.1:c.*1701A>C",
            ]
        );
    }

    /// The misalignment offset of an `NM_` model stands where the whole
    /// alternate allele differs from the bases of the clamped span
    /// (`get_misalignment_offset`, TranscriptVariationAllele.pm 1517,
    /// 1594-1624, reached once 1474-1483 type a notation), whether or not an
    /// HGVSc follows. The fixture's alignment added `AAA` after cDNA 9
    /// (`_rna_edit` `10 9 AAA`), so its mapper numbers the last base 2604,
    /// c.*1701 past the coding end 903, and an allele from cDNA 10 on carries
    /// 3. `AA` to `AG` over the last base and the base after it changes no
    /// base of the model: no notation, offset 3. `AAA` to `A` over the last
    /// base and the two after it keeps the base inside, and its alternate
    /// equals that base: no notation and no offset, as Ensembl VEP writes
    /// neither. `AA` to `CA` writes `c.*1701A>C` beside the offset, the base
    /// after the last base alone lies outside the model and carries neither,
    /// and the same model under an `ENST` name has no offset. Ensembl VEP
    /// 116.3 writes the offset of the first record beside an `ins` of the
    /// alternate's bases past the last base, the whole alternate read against
    /// the base inside (1474-1483).
    #[test]
    fn edited_model_refseq_offset_stands_on_a_row_without_a_notation() {
        use vep_core::transcript::{Attribute, ExonCoordMapper, MapperPair};
        let insertion_edited = |stable_id: &str| {
            let mut model = make_test_transcript();
            model.stable_id = stable_id.into();
            model.attributes.push(Attribute {
                code: "_rna_edit".into(),
                value: "10 9 AAA".to_string(),
            });
            let mapper = model.vefc.as_mut().unwrap().mapper.as_mut().unwrap();
            mapper.cdna_coding_start = 54;
            mapper.cdna_coding_end = 903;
            let pair = |from_start, from_end, to_start, to_end| MapperPair {
                from_start,
                from_end,
                to_start,
                to_end,
                ori: 1,
            };
            mapper.exon_coord_mapper = ExonCoordMapper::new(vec![
                pair(1, 9, 25_000_000, 25_000_008),
                pair(13, 303, 25_000_009, 25_000_299),
                pair(304, 603, 25_002_000, 25_002_299),
                pair(604, 2604, 25_004_000, 25_006_000),
            ]);
            model.derived = Default::default();
            model
        };
        let hgvs = |model: &Transcript, start: u64, end: u64, reference: &[u8], alt: &[u8]| {
            let variant =
                InputVariant::new("21".into(), start, end, reference.to_vec(), alt.to_vec());
            let hgvs = generate_hgvs(&variant, model, None, SHIFTED);
            (hgvs.hgvsc.unwrap_or_default(), hgvs.refseq_offset)
        };
        let refseq = insertion_edited("NM_000003.1");
        let ensembl = insertion_edited("ENST00000000003");
        assert_eq!(
            [
                hgvs(&refseq, 25_006_000, 25_006_001, b"AA", b"AG"),
                hgvs(&refseq, 25_006_000, 25_006_002, b"AAA", b"A"),
                hgvs(&refseq, 25_006_000, 25_006_001, b"AA", b"CA"),
                hgvs(&refseq, 25_006_001, 25_006_001, b"A", b"G"),
                hgvs(&ensembl, 25_006_000, 25_006_001, b"AA", b"AG"),
            ],
            [
                (String::new(), Some(3)),
                (String::new(), None),
                ("NM_000003.1:c.*1701A>C".to_string(), Some(3)),
                (String::new(), None),
                (String::new(), None),
            ]
        );
    }

    /// The same insertion on a reverse-strand edited model: the bases `CTGGCA`
    /// read on the transcript's strand at c.202-c.207 are `TGCCAG` on the
    /// forward strand at genomic 5844-5849 (CDS position 6051 minus the genomic
    /// one), and the inserted `GA` of the transcript is the forward-strand `TC`
    /// between 5847 and 5848. The shift moves it one base towards the
    /// transcript's 3' end, genomic 5846-5847, offset -1, and the allele turns
    /// the other way on the forward strand (`CT`), `AG` on the transcript's: both
    /// runs write `c.204_205insAG`, where the allele as written spells
    /// `c.204_205insGA`.
    #[test]
    fn reverse_strand_edited_model_hgvsc_spells_a_shifted_insertion_from_the_rotated_allele() {
        use crate::concordance_tests::{
            bases_from, compact_on_chr1, hgvsc_fasta_fixture,
            make_descending_reverse_strand_transcript,
        };
        let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(5844, "TGCCAG"));
        let insertion = InputVariant::new("1".into(), 5848, 5847, b"-".to_vec(), b"TC".to_vec());
        let edited = substitution_edited_model(
            compact_on_chr1(make_descending_reverse_strand_transcript()),
            "NM_000002.1",
        );
        for options in [SHIFTED, SHIFT_3PRIME] {
            let notation = generate_hgvs(&insertion, &edited, Some(&fasta), options);
            assert_eq!(
                notation.hgvsc.as_deref(),
                Some("NM_000002.1:c.204_205insAG"),
                "{options:?}"
            );
            assert_eq!(notation.offset, Some(-1), "{options:?}");
        }
        assert_eq!(
            generate_hgvs(&insertion, &edited, Some(&fasta), UNSHIFTED)
                .hgvsc
                .as_deref(),
            Some("NM_000002.1:c.203_204insGA")
        );
    }

    /// A row in an untranslated region of an edited model carries the rotated
    /// bases as a coding row does, under `--shift_3prime 1` too. Ensembl VEP
    /// 116.3 keeps the bases as written there in both runs: `codon` returns
    /// before `shift_feature_seqs` on a row without a translation position
    /// (TranscriptVariationAllele.pm 815-820), and `hgvs_transcript` writes
    /// the unrotated `feature_seq` over the shifted notation (1513). The
    /// forward-strand fixture lays NM_006933.7's bases around c.*1349,
    /// `CTGTGT` (GRCh37 21:35471001-35471006), on c.*1347-c.*1352 (genomic 4399
    /// plus the distance past the stop codon): `GA` inserted after c.*1348
    /// shifts one base and reads `c.*1349_*1350insAG`, where Ensembl VEP writes
    /// `c.*1349_*1350insGA` for `21 35471002 T TGA` on that model. The
    /// reverse-strand fixture lays the forward-strand bases around
    /// NM_001395858.1's c.-57, `GGCAAT` (GRCh37 21:16031360-16031365), at
    /// genomic 6084-6089 of the 5' UTR (genomic 6050 plus the distance before
    /// the start codon): the forward-strand `TC` between 6086 and 6087 is the
    /// transcript's `GA` after c.-37, both runs write `c.-36_-35insAG`, and
    /// Ensembl VEP writes `c.-57_-56insGA` for `21 16031362 C CTC` on that
    /// model where the rotated bases read `c.-57_-56insAG`.
    #[test]
    fn edited_model_hgvsc_rotates_a_shifted_insertion_on_an_untranslated_row() {
        use crate::concordance_tests::{
            bases_from, compact_on_chr1, hgvsc_fasta_fixture, make_compact_test_transcript,
            make_descending_reverse_strand_transcript,
        };
        let forward = substitution_edited_model(make_compact_test_transcript(), "NM_000001.1");
        let reverse = substitution_edited_model(
            compact_on_chr1(make_descending_reverse_strand_transcript()),
            "NM_000002.1",
        );
        let three_prime = InputVariant::new("1".into(), 5748, 5747, b"-".to_vec(), b"GA".to_vec());
        let five_prime = InputVariant::new("1".into(), 6087, 6086, b"-".to_vec(), b"TC".to_vec());
        let rows = [
            (
                &three_prime,
                &forward,
                bases_from(5746, "CTGTGT"),
                "NM_000001.1:c.*1349_*1350insAG",
                1,
                "NM_000001.1:c.*1348_*1349insGA",
            ),
            (
                &five_prime,
                &reverse,
                bases_from(6084, "GGCAAT"),
                "NM_000002.1:c.-36_-35insAG",
                -1,
                "NM_000002.1:c.-37_-36insGA",
            ),
        ];
        for (insertion, edited, bases, shifted, offset, as_written) in rows {
            let (_dir, fasta) = hgvsc_fasta_fixture(&bases);
            for options in [SHIFTED, SHIFT_3PRIME] {
                let notation = generate_hgvs(insertion, edited, Some(&fasta), options);
                assert_eq!(notation.hgvsc.as_deref(), Some(shifted), "{options:?}");
                assert_eq!(notation.offset, Some(offset), "{options:?}");
            }
            assert_eq!(
                generate_hgvs(insertion, edited, Some(&fasta), UNSHIFTED)
                    .hgvsc
                    .as_deref(),
                Some(as_written)
            );
        }
    }

    /// The module documentation names what the HGVS port leaves out of Ensembl
    /// VEP's `hgvs_transcript` and `hgvs_protein`: the second shift along an
    /// edited RefSeq transcript's own sequence and its GenBank protein name;
    /// where Ensembl VEP itself returns before that shift, and where the port
    /// shifts along the genome and where along the model's own sequence;
    /// that an LRG reference name takes no version (`versioned_reference_name`);
    /// and what it covers of an edited model's HGVSc, the misalignment offset
    /// and the `shift_feature_seqs` rotation of a shifted insertion's bases.
    #[test]
    fn documented_subset_names_the_refseq_edited_transcript_paths() {
        let module_docs: String = include_str!("hgvs.rs")
            .lines()
            .take_while(|l| l.starts_with("//!") || l.starts_with("//") || l.trim().is_empty())
            .filter(|l| l.starts_with("//!"))
            .map(|l| l.trim_start_matches("//!").trim())
            .collect::<Vec<_>>()
            .join(" ");
        for phrase in [
            "Two paths of `hgvs_transcript` and `hgvs_protein` stay out: the second shift along a RefSeq transcript's own edited sequence (the tail of `_return_3prime`, 211-256) and the GenBank protein name of such a transcript (1750-1756).",
            "Away from a transcript's last bases and its edited bases Ensembl VEP returns before that shift, copying a shift it already has (182-187).",
            "The port shifts an insertion or a deletion on an edited transcript 3' along the genome (`hgvs_shift`), and along the model's own sequence only on or beside bases the model holds and the genome lacks (`place_3prime_on_model`).",
            "edited sequence",
            "`_return_3prime`",
            "GenBank",
            "An LRG reference name carries no version, as there",
            "cached identifiers",
            "misalignment offset",
            "`shift_feature_seqs`",
        ] {
            assert!(module_docs.contains(phrase), "{phrase}: {module_docs}");
        }
        assert!(
            !module_docs.contains("and the LRG reference names"),
            "the LRG rule is ported, not left out: {module_docs}"
        );
        assert!(
            !module_docs.contains("147-258"),
            "the second shift is lines 211-256 of the sub's tail: {module_docs}"
        );
    }

    /// A slice with no FASTA: `substr` reads nothing, so the reference comes from
    /// the caller and no duplication lookup can succeed.
    fn unreadable_slice() -> TranscriptSlice<'static> {
        TranscriptSlice {
            chr: "21",
            tr_start: 1,
            tr_end: 10_000,
            reverse: false,
            fasta: None,
        }
    }

    #[test]
    fn test_versioned_reference_name() {
        assert_eq!(
            versioned_reference_name("ENST00000000001", Some(1)),
            "ENST00000000001.1"
        );
        assert_eq!(
            versioned_reference_name("ENST00000000001", None),
            "ENST00000000001"
        );
        assert_eq!(
            versioned_reference_name("NM_000001.2", Some(3)),
            "NM_000001.2"
        );
        assert_eq!(versioned_reference_name("LRG_1t1", Some(1)), "LRG_1t1");
        assert_eq!(
            versioned_reference_name("ENSP00000000001", Some(7)),
            "ENSP00000000001.7"
        );
    }

    #[test]
    fn test_hgvs_variant_notation_types() {
        let slice = unreadable_slice();
        let kind = |r: &[u8], a: &[u8]| {
            let end = 100 + r.len() as i64 - 1;
            hgvs_variant_notation(&slice, r, a.to_vec(), 100, end, DupLookup::AfterFirst)
                .map(|n| (n.kind, n.start, n.end))
        };
        assert_eq!(kind(b"A", b"G"), Some((HgvscKind::Sub, 100, 100)));
        assert_eq!(kind(b"ACG", b""), Some((HgvscKind::Del, 100, 102)));
        assert_eq!(kind(b"ACG", b"CGT"), Some((HgvscKind::Inv, 100, 102)));
        assert_eq!(kind(b"ACG", b"TTT"), Some((HgvscKind::Delins, 100, 102)));
        assert_eq!(kind(b"AC", b"ACAC"), Some((HgvscKind::Dup, 100, 101)));
        assert_eq!(
            kind(b"AC", b"ACACAC"),
            Some((HgvscKind::Multiple(3), 100, 101))
        );
        assert_eq!(kind(b"AC", b"ACACAG"), Some((HgvscKind::Delins, 100, 101)));
        assert_eq!(kind(b"ACG", b"ACG"), None);
        // An insertion (`ref_end < ref_start`) lists the smaller coordinate first.
        assert_eq!(
            hgvs_variant_notation(&slice, b"", b"TT".to_vec(), 100, 99, DupLookup::AfterFirst)
                .map(|n| (n.kind, n.start, n.end)),
            Some((HgvscKind::Ins, 99, 100))
        );
    }

    #[test]
    fn test_hgvsc_clip_alleles_trims_shared_flanks_and_retypes() {
        let mut n = HgvscNotation {
            start: 10,
            end: 14,
            ref_seq: b"GATTC".to_vec(),
            alt_seq: b"GACTC".to_vec(),
            kind: HgvscKind::Delins,
        };
        hgvsc_clip_alleles(&mut n);
        assert_eq!((n.start, n.end), (12, 12));
        assert_eq!(
            (n.ref_seq.as_slice(), n.alt_seq.as_slice()),
            (&b"T"[..], &b"C"[..])
        );
        assert_eq!(n.kind, HgvscKind::Sub);

        let mut n = HgvscNotation {
            start: 10,
            end: 12,
            ref_seq: b"GAT".to_vec(),
            alt_seq: b"GATTC".to_vec(),
            kind: HgvscKind::Delins,
        };
        hgvsc_clip_alleles(&mut n);
        assert_eq!((n.start, n.end), (13, 12));
        assert_eq!(n.alt_seq, b"TC");
        assert_eq!(n.kind, HgvscKind::Ins);

        let mut n = HgvscNotation {
            start: 10,
            end: 14,
            ref_seq: b"GATTC".to_vec(),
            alt_seq: b"GTC".to_vec(),
            kind: HgvscKind::Delins,
        };
        hgvsc_clip_alleles(&mut n);
        assert_eq!((n.start, n.end), (11, 12));
        assert_eq!(n.ref_seq, b"AT");
        assert_eq!(n.kind, HgvscKind::Del);
    }

    /// An equal-length record whose alleles differ at two places with a base
    /// they share between them is one `delins` over its whole span, as Ensembl
    /// VEP writes it: `hgvs_variant_notation` classes alleles of one length,
    /// two bases or more and no inversion, as a `delins`
    /// (Utils/Sequence.pm:548-566), and `_clip_alleles` strips the bases they
    /// share from the two ends alone (TranscriptVariationAllele.pm:2203-2242),
    /// so the shared base stays in the notation. The protein notation follows
    /// the residues: a pair inside one codon is a substitution, and a pair
    /// that changes two neighbouring residues is one `delins` over both
    /// (TranscriptVariationAllele.pm:2075-2077).
    #[test]
    fn two_changes_a_base_apart_are_written_as_one_deletion_insertion() {
        let tx = make_test_transcript();
        for (place, start, reference, alt, hgvsc, hgvsp) in [
            (
                "bases 4 and 6 of intron 1",
                25_000_303,
                b"ACA",
                b"GCT",
                "ENST00000000001.1:c.250+4_250+6delinsGCT",
                None,
            ),
            (
                "bases 1 and 3 of codon 4",
                25_000_059,
                b"AAA",
                b"GAC",
                "ENST00000000001.1:c.10_12delinsGAC",
                Some("ENSP00000000001.1:p.Lys4Asp"),
            ),
            (
                "base 3 of codon 4 and base 2 of codon 5",
                25_000_061,
                b"ATT",
                b"CTC",
                "ENST00000000001.1:c.12_14delinsCTC",
                Some("ENSP00000000001.1:p.Lys4_Phe5delinsAsnSer"),
            ),
        ] {
            let variant = InputVariant::new(
                "21".into(),
                start,
                start + 2,
                reference.to_vec(),
                alt.to_vec(),
            );
            let notation = generate_hgvs(&variant, &tx, None, SHIFTED);
            assert_eq!(notation.hgvsc.as_deref(), Some(hgvsc), "{place}");
            assert_eq!(notation.hgvsp.as_deref(), hgvsp, "{place}");
        }
    }

    #[test]
    fn test_format_hgvs_string_bodies() {
        let n = |kind: HgvscKind, r: &[u8], a: &[u8]| HgvscNotation {
            start: 0,
            end: 0,
            ref_seq: r.to_vec(),
            alt_seq: a.to_vec(),
            kind,
        };
        let f =
            |kind, r: &[u8], a: &[u8], s, e| format_hgvs_string("T.1", 'c', s, e, &n(kind, r, a));
        assert_eq!(f(HgvscKind::Sub, b"A", b"G", "5", "5"), "T.1:c.5A>G");
        assert_eq!(f(HgvscKind::Del, b"AC", b"", "5", "6"), "T.1:c.5_6del");
        assert_eq!(f(HgvscKind::Del, b"A", b"", "5", "5"), "T.1:c.5del");
        assert_eq!(f(HgvscKind::Inv, b"AC", b"GT", "5", "6"), "T.1:c.5_6inv");
        assert_eq!(f(HgvscKind::Dup, b"", b"AC", "5", "6"), "T.1:c.5_6dup");
        assert_eq!(f(HgvscKind::Ins, b"", b"AC", "5", "6"), "T.1:c.5_6insAC");
        assert_eq!(
            f(HgvscKind::Delins, b"AC", b"TTT", "250+3", "251-2"),
            "T.1:c.250+3_251-2delinsTTT"
        );
        assert_eq!(
            f(HgvscKind::Multiple(3), b"AC", b"ACACAC", "5", "6"),
            "T.1:c.5_6[3]"
        );
    }

    #[test]
    fn test_hgvsp_missense() {
        let tx = make_test_transcript();
        // CDS pos 5, GCT -> GAT = Ala -> Asp
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let hgvsp = generate_hgvsp(&variant, &tx, None);
        assert!(hgvsp.is_some());
        let notation = hgvsp.unwrap();
        assert!(
            notation.contains("Ala2Asp"),
            "Missense should show Ala2Asp, got: {notation}"
        );
    }

    #[test]
    fn test_hgvsp_synonymous() {
        let tx = make_test_transcript();
        // CDS pos 6, GCT -> GCC = Ala -> Ala (synonymous)
        let variant = InputVariant::new(
            "21".into(),
            25_000_055,
            25_000_055,
            b"T".to_vec(),
            b"C".to_vec(),
        );
        let hgvsp = generate_hgvsp(&variant, &tx, None);
        assert!(hgvsp.is_some());
        let notation = hgvsp.unwrap();
        assert!(
            notation.contains("Ala2="),
            "Synonymous should show Ala2=, got: {notation}"
        );
    }

    #[test]
    fn test_hgvsp_stop_gained() {
        let tx = make_test_transcript();
        // CDS pos 7, GGA -> TGA = Gly -> Stop
        let variant = InputVariant::new(
            "21".into(),
            25_000_056,
            25_000_056,
            b"G".to_vec(),
            b"T".to_vec(),
        );
        let hgvsp = generate_hgvsp(&variant, &tx, None);
        assert!(hgvsp.is_some());
        let notation = hgvsp.unwrap();
        assert!(
            notation.contains("Gly3Ter"),
            "Stop gained should show Gly3Ter, got: {notation}"
        );
    }

    /// `_get_hgvs_protein_format` opens a parenthesis right after `p.` and closes
    /// it after the whole description when asked for the predicted format
    /// (TranscriptVariationAllele.pm 1908, 2040), which `--hgvsp_use_prediction`
    /// passes down (OutputFactory.pm 1769). A synonymous description keeps its `=`
    /// inside the parentheses here; the output layer escapes it afterwards
    /// (OutputFactory.pm 1773), so the text formats print `p.(Ala2%3D)`.
    #[test]
    fn hgvsp_prediction_format_wraps_the_description() {
        let tx = make_test_transcript();
        let predicted = HgvsOptions {
            hgvsp_use_prediction: true,
            ..SHIFTED
        };
        let hgvsp = |pos: u64, reference: &[u8], alt: &[u8], options: HgvsOptions| {
            let variant =
                InputVariant::new("21".into(), pos, pos, reference.to_vec(), alt.to_vec());
            generate_hgvs(&variant, &tx, None, options).hgvsp
        };
        assert_eq!(
            hgvsp(25_000_054, b"C", b"A", predicted).as_deref(),
            Some("ENSP00000000001.1:p.(Ala2Asp)")
        );
        assert_eq!(
            hgvsp(25_000_056, b"G", b"T", predicted).as_deref(),
            Some("ENSP00000000001.1:p.(Gly3Ter)")
        );
        assert_eq!(
            hgvsp(25_000_055, b"T", b"C", predicted).as_deref(),
            Some("ENSP00000000001.1:p.(Ala2=)")
        );
        assert_eq!(
            hgvsp(25_000_055, b"T", b"C", SHIFTED).as_deref(),
            Some("ENSP00000000001.1:p.Ala2=")
        );
    }

    /// `hgvs_protein` appends the translation's version to the reference name
    /// unless the version is undefined, the name already carries one, the name is
    /// an LRG, or `remove_hgvsp_version` is set (TranscriptVariationAllele.pm
    /// 1759-1760); the flag skips the append and strips nothing. `hgvs_transcript`
    /// (1499) has no such flag, so HGVSc keeps its version.
    #[test]
    fn remove_hgvsp_version_drops_the_translation_version_only() {
        let tx = make_test_transcript();
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let unversioned = HgvsOptions {
            remove_hgvsp_version: true,
            ..SHIFTED
        };
        let notation = generate_hgvs(&variant, &tx, None, unversioned);
        assert_eq!(notation.hgvsp.as_deref(), Some("ENSP00000000001:p.Ala2Asp"));
        assert_eq!(notation.hgvsc.as_deref(), Some("ENST00000000001.1:c.5C>A"));
        assert_eq!(
            generate_hgvs(&variant, &tx, None, SHIFTED).hgvsp.as_deref(),
            Some("ENSP00000000001.1:p.Ala2Asp")
        );
        assert_eq!(
            generate_hgvs(
                &variant,
                &tx,
                None,
                HgvsOptions {
                    hgvsp_use_prediction: true,
                    ..unversioned
                }
            )
            .hgvsp
            .as_deref(),
            Some("ENSP00000000001:p.(Ala2Asp)")
        );

        let mut versioned_id = make_test_transcript();
        versioned_id.protein_id = Some("NP_001234.2".into());
        assert_eq!(
            generate_hgvs(&variant, &versioned_id, None, unversioned)
                .hgvsp
                .as_deref(),
            Some("NP_001234.2:p.Ala2Asp")
        );
    }

    #[test]
    fn test_hgvsp_start_lost() {
        let tx = make_test_transcript();
        // ATG -> ACG = Met -> Thr (start_lost)
        let variant = InputVariant::new(
            "21".into(),
            25_000_051,
            25_000_051,
            b"T".to_vec(),
            b"C".to_vec(),
        );
        let hgvsp = generate_hgvsp(&variant, &tx, None);
        assert!(hgvsp.is_some());
        let notation = hgvsp.unwrap();
        assert!(
            notation.contains("Met1?"),
            "Start lost should show Met1?, got: {notation}"
        );
    }

    #[test]
    fn test_hgvsp_frameshift() {
        let tx = make_test_transcript();
        // 2 bp insertion after the start codon (CDS 3|4): the alternate frame reads
        // ATG AAG CTG ..., so residue 2 becomes Lys and no stop precedes the
        // N-padded 3' UTR. An insertion inside ATG is a start loss (`Met1?`).
        let variant = InputVariant::new(
            "21".into(),
            25_000_053,
            25_000_052,
            b"-".to_vec(),
            b"AA".to_vec(),
        );
        let hgvsp = generate_hgvsp(&variant, &tx, None);
        assert_eq!(hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ala2LysfsTer?"));
    }

    /// The fixture transcript with an 849-base CDS ending in a stop codon: 282
    /// residues (codons 7 to 282 are `GCT`, Ala) and `TAA` at CDS 847-849, genomic
    /// 25_004_296..=25_004_298, followed by the given 3' UTR, which the alternate
    /// CDS of a frameshift reads through.
    fn terminal_stop_transcript(three_prime_utr: &str) -> Transcript {
        let mut tx = make_test_transcript();
        let mut cds = String::from("ATGGCTGGAAAATTCGAT");
        while cds.len() < 846 {
            cds.push_str("GCT");
        }
        cds.push_str("TAA");
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.translateable_seq = Some(cds);
        vefc.three_prime_utr = Some(three_prime_utr.to_string());
        vefc.mapper.as_mut().unwrap().cdna_coding_end = 899;
        tx.cdna_coding_end = Some(899);
        tx.coding_region_end = Some(25_004_298);
        tx.translation_end = Some(25_004_298);
        tx.translation.as_mut().unwrap().end = 299;
        tx
    }

    /// The `fs` branch of `_get_hgvs_protein_format` (TranscriptVariationAllele.pm
    /// 2007-2010): deleting the last base of codon 282 and the first of the stop
    /// leaves Ala282 as `GCA`, so the first changed residue is Ter283, read as
    /// `AGC` (Ser); the new frame stops at residue 284 (`TAA` from the UTR),
    /// `_stop_loss_extra_AA` counts 284 - 282 = 2, and the extension drops the
    /// replaced stop from the count: `extTer1`. Rows of this shape in
    /// `tests/golden/116/GRCh37-hgvs`: `21:35743146-35743147 -` against
    /// ENST00000290310, `ENSP00000290310.2:p.Ter124IleextTer14`, and
    /// `21:45210898 -` against ENST00000483896, `ENSP00000426898.1:p.Ter46AsnextTer78`.
    #[test]
    fn test_hgvsp_frameshift_starting_at_stop_codon_is_extension() {
        let tx = terminal_stop_transcript("GCTAAGCATAA");
        let variant = InputVariant::new(
            "21".into(),
            25_004_295,
            25_004_296,
            b"TT".to_vec(),
            b"-".to_vec(),
        );
        assert_eq!(
            generate_hgvsp(&variant, &tx, None).as_deref(),
            Some("ENSP00000000001.1:p.Ter283SerextTer1")
        );
    }

    /// The same deletion with no stop in the new frame: `_stop_loss_extra_AA`
    /// returns `undef`, the count is `?`, and `?` is not decremented.
    #[test]
    fn test_hgvsp_frameshift_starting_at_stop_codon_without_new_stop_keeps_query() {
        let tx = terminal_stop_transcript("GCAGCAGCA");
        let variant = InputVariant::new(
            "21".into(),
            25_004_295,
            25_004_296,
            b"TT".to_vec(),
            b"-".to_vec(),
        );
        assert_eq!(
            generate_hgvsp(&variant, &tx, None).as_deref(),
            Some("ENSP00000000001.1:p.Ter283SerextTer?")
        );
    }

    /// A frameshift whose first changed residue is not the stop keeps `fsTer`
    /// with the full count: deleting CDS 845-846 turns Ala282 into `GTA` (Val),
    /// and the new frame stops at residue 284, so 284 - 281 = 3.
    #[test]
    fn test_hgvsp_frameshift_before_stop_codon_keeps_fs_ter_count() {
        let tx = terminal_stop_transcript("GCTAAGCATAA");
        let variant = InputVariant::new(
            "21".into(),
            25_004_294,
            25_004_295,
            b"CT".to_vec(),
            b"-".to_vec(),
        );
        assert_eq!(
            generate_hgvsp(&variant, &tx, None).as_deref(),
            Some("ENSP00000000001.1:p.Ala282ValfsTer3")
        );
    }

    /// `frameshift` (VariationEffect.pm 1554) returns 0 when the reference peptide
    /// starts with the stop, so a one-base deletion inside the stop codon is typed
    /// by its peptides: the two bases left of the stop codon translate to `X`,
    /// three-lettered `Xaa` and rewritten `Ter`, equal to the reference `Ter`, which
    /// `hgvs_protein` prints as the stop unchanged. The stop is lost: the whole
    /// protein is kept and residues follow it to the transcript's end without a stop, so
    /// the notation is an extension named at the stop. Rows of this shape in
    /// `tests/golden/116/GRCh37-hgvs`: `21:34955948 -` against ENST00000457359 and
    /// `21:35791518 -` against ENST00000450895.
    #[test]
    fn test_hgvsp_deletion_inside_stop_codon_is_not_a_frameshift() {
        let tx = terminal_stop_transcript("CCAGCAGCA");
        let variant = InputVariant::new(
            "21".into(),
            25_004_297,
            25_004_297,
            b"A".to_vec(),
            b"-".to_vec(),
        );
        assert_eq!(
            generate_hgvsp(&variant, &tx, None).as_deref(),
            Some("ENSP00000000001.1:p.Ter283TyrextTer?")
        );
    }

    /// A shifted allele keeps the two alleles as the record writes them: the
    /// row's own where the record's were not cut, the kept string otherwise.
    #[test]
    fn a_shifted_allele_carries_the_alleles_as_the_record_writes_them() {
        let tx = make_test_transcript();
        let insertion = InputVariant::new(
            "21".into(),
            25_000_059,
            25_000_058,
            b"-".to_vec(),
            b"TAANNN".to_vec(),
        );
        let moved = shifted_variant(
            &insertion,
            Strand::Forward,
            &tx,
            (25_000_061, 25_000_060),
            None,
        )
        .expect("an insertion has a shifted allele");
        assert_eq!(moved.original_allele_string.as_deref(), Some("-/TAANNN"));
        let mut cut = insertion.clone();
        cut.original_allele_string = Some("GGA/NNNTTAGGA".to_string());
        let moved = shifted_variant(&cut, Strand::Forward, &tx, (25_000_061, 25_000_060), None)
            .expect("an insertion has a shifted allele");
        assert_eq!(
            moved.original_allele_string.as_deref(),
            Some("GGA/NNNTTAGGA")
        );
    }
}
