// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! HGVS notation generation: HGVSc (coding DNA) and HGVSp (protein) for a
//! variant-transcript pair, and HGVSg (genomic) for a variant allele.
//!
//! The strings follow the Human Genome Variation Society nomenclature
//! conventions, matching Perl VEP's output.
//!
//! Perl citations name modules of ensembl-variation release/116
//! (`Bio/EnsEMBL/Variation/...`).
//!
//! Reference: <https://varnomen.hgvs.org/>

use std::borrow::Cow;

use vep_core::coordinate::Strand;
use vep_core::transcript::{MapperPair, Transcript};
use vep_core::variant::{InputVariant, VariantClass};

/// Convert a single-letter amino acid code to its three-letter abbreviation.
pub fn amino_acid_three_letter(one_letter: u8) -> &'static str {
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

/// How far Perl VEP's `perform_shift` (TranscriptVariationAllele.pm 291) can
/// move an insertion or deletion of `len` bases: it compares against 1,000 bases
/// of flank, so a pattern that fits the flank moves at most `1001 - len`
/// positions and a longer one at most 1,000, except that on the reverse strand
/// `len == 1001` gets a loop bound of zero and never moves.
fn hgvs_shift_limit(len: usize, reverse: bool) -> u64 {
    match len {
        0..=1000 => 1001 - len as u64,
        1001 if reverse => 0,
        _ => 1000,
    }
}

/// The 3'-shifted span of an insertion or deletion, Perl's `_return_3prime(1)`
/// (TranscriptVariationAllele.pm 109): along the genome in the transcript's 3'
/// direction, whatever `--shift_3prime` says, when a reference FASTA is present
/// and the shift moves the variant.
fn hgvs_shift(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<(u64, u64)> {
    let len = match variant.variant_class {
        VariantClass::Insertion => variant.alt_allele().len(),
        VariantClass::Deletion => variant.ref_allele.len(),
        _ => return None,
    };
    let fasta = reference_fasta?;
    let reverse = transcript.strand == Strand::Reverse;
    crate::consequences::shift_indel_3prime_coords(
        variant,
        transcript,
        fasta,
        hgvs_shift_limit(len, reverse),
    )
    .filter(|&(start, end)| (start, end) != (variant.start, variant.end))
}

/// The transcript's genomic span read in transcript orientation (Perl's
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

    /// Perl `substr($slice->seq, $start - 1, $len)`, or `None` when the FASTA is
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
}

/// Perl's HGVS variant type (`hgvs_variant_notation`, Utils/Sequence.pm 493),
/// shared by the transcript (`c.`/`n.`) and genomic (`g.`) descriptions.
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

/// Perl's `$hgvs_notation` hash for one description.
struct HgvscNotation {
    start: i64,
    end: i64,
    ref_seq: Vec<u8>,
    alt_seq: Vec<u8>,
    kind: HgvscKind,
}

/// `hgvs_variant_notation` (Utils/Sequence.pm 493) with the default lookup
/// order: the type and displayed span of `alt_seq` replacing reference positions
/// `ref_start..=ref_end`. `read(start, len)` returns the `len` reference bases
/// from 1-based `start`, `None` when they cannot be read; the reference allele
/// falls back to `fallback_ref` when given. `None` where Perl returns `undef`
/// (the alleles are equal) or the reference bases cannot be read.
fn hgvs_variant_notation(
    read: impl Fn(i64, usize) -> Option<Vec<u8>>,
    fallback_ref: Option<&[u8]>,
    alt_seq: Vec<u8>,
    ref_start: i64,
    ref_end: i64,
) -> Option<HgvscNotation> {
    let ref_length = (ref_end - ref_start + 1).max(0) as usize;
    let ref_seq = match read(ref_start, ref_length) {
        Some(s) => s,
        None => fallback_ref?.to_vec(),
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
    if alt_length == 0 {
        return Some(n);
    }
    if ref_length == alt_length {
        n.kind = if ref_length == 1 {
            HgvscKind::Sub
        } else if n.alt_seq == crate::coding::reverse_complement(&n.ref_seq) {
            HgvscKind::Inv
        } else {
            HgvscKind::Delins
        };
        return Some(n);
    }
    if ref_length == 0 {
        // Perl's `lookup_order` 1: the bases after the site are tried first.
        let after = read(ref_end + 1, alt_length);
        if after.as_deref() == Some(n.alt_seq.as_slice()) {
            n.end = n.start + alt_length as i64 - 1;
            n.kind = HgvscKind::Dup;
            return Some(n);
        }
        let before = read(ref_end - alt_length as i64 + 1, alt_length);
        if before.as_deref() == Some(n.alt_seq.as_slice()) {
            n.start = n.end - alt_length as i64 + 1;
            n.kind = HgvscKind::Dup;
            return Some(n);
        }
        (n.start, n.end) = (ref_end, ref_start);
        n.kind = HgvscKind::Ins;
        return Some(n);
    }
    if alt_length.is_multiple_of(ref_length) {
        let multiple = alt_length / ref_length;
        if n.alt_seq == n.ref_seq.repeat(multiple) {
            n.kind = if multiple == 2 {
                HgvscKind::Dup
            } else {
                HgvscKind::Multiple(multiple)
            };
            return Some(n);
        }
    }
    n.kind = HgvscKind::Delins;
    Some(n)
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

/// One exon of the transcript as `_get_cDNA_position` walks them: genomic span
/// and cDNA span, from the transcript mapper's pairs.
struct ExonSpan {
    start: i64,
    end: i64,
    cdna_start: i64,
    cdna_end: i64,
}

impl ExonSpan {
    /// The mapper's exon pairs sorted by genomic start, Perl's `_sorted_exons`.
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

/// `_get_cDNA_position` (TranscriptVariationAllele.pm 2755): the HGVS cDNA
/// coordinate of 1-based slice position `pos`, exonic as a plain coordinate,
/// intronic as the nearest exon boundary with a `+` or `-` distance (the
/// upstream exon on a tie), then made relative to the start codon (`-` before it)
/// and the stop codon (`*` after it) on a coding transcript.
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
            // Perl drops the coordinate and the `+` here: the stop codon's last base
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

/// The `(exon coordinate, intron offset)` pair Perl reads back from an HGVS
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

/// The HGVS output of one variant-transcript pair under `--hgvs`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct HgvsNotation {
    /// HGVSc, `generate_hgvsc`.
    pub hgvsc: Option<String>,
    /// HGVSp, `generate_hgvsp`.
    pub hgvsp: Option<String>,
    /// VEP's `HGVS_OFFSET`: the bases the insertion or deletion was shifted 3'
    /// along the transcript, negative on the reverse strand, present only when
    /// the shift moved it and at least one notation was produced.
    pub offset: Option<i64>,
}

/// HGVSc, HGVSp and the shift offset for a variant-transcript pair, the indel
/// shift computed once for both notations.
pub fn generate_hgvs(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> HgvsNotation {
    let shifted = hgvs_shift(variant, transcript, reference_fasta);
    let hgvsc = hgvsc_at(variant, transcript, reference_fasta, shifted);
    let hgvsp = hgvsp_at(variant, transcript, reference_fasta, shifted);
    let offset = shifted
        .filter(|_| hgvsc.is_some() || hgvsp.is_some())
        .map(|(start, _)| {
            let length = start.abs_diff(variant.start) as i64;
            if transcript.strand == Strand::Reverse {
                -length
            } else {
                length
            }
        });
    HgvsNotation {
        hgvsc,
        hgvsp,
        offset,
    }
}

/// Generate HGVSc notation for a variant-transcript pair: Perl VEP's
/// `hgvs_transcript` (TranscriptVariationAllele.pm 1375) for an Ensembl
/// transcript.
///
/// Returns `None` when the alternate allele carries a character outside
/// `ACGT-` or equals the reference, when the variant does not overlap the
/// transcript's genomic span (so an upstream or downstream variant has no HGVSc),
/// when the 3'-shifted span leaves it, or when the alleles agree at the shifted
/// position. A variant that overlaps the span but runs over either edge is
/// described over the bases inside it, each end clamped onto the transcript
/// (`_var2transcript_slice_coords`, 2732-2733). An insertion or deletion is
/// described at its most 3' position on the transcript strand when a FASTA is
/// present; the type and the displayed span follow `hgvs_variant_notation`,
/// shared bases are clipped, an exonic single-base substitution inside the CDS
/// takes its CDS coordinate and every other position goes through
/// `_get_cDNA_position`. Without a FASTA the reference bases come from the
/// variant and no duplication is detected.
pub fn generate_hgvsc(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<String> {
    let shifted = hgvs_shift(variant, transcript, reference_fasta);
    hgvsc_at(variant, transcript, reference_fasta, shifted)
}

/// `generate_hgvsc` with the 3' shift already computed.
fn hgvsc_at(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
    shifted: Option<(u64, u64)>,
) -> Option<String> {
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
    let offset = shifted.map_or(0, |(start, _)| start.abs_diff(variant.start) as i64);

    // `hgvs_transcript` reverse-complements an allele onto the transcript strand
    // and `hgvs_variant_notation` strips its gap characters before reading it.
    let in_transcript_orientation = |mut allele: Vec<u8>| {
        allele.make_ascii_uppercase();
        allele.retain(|&b| b != b'-');
        if reverse {
            allele = crate::coding::reverse_complement(&allele);
        }
        allele
    };
    // A shifted insertion is described by its rotated allele.
    let alt_seq = in_transcript_orientation(match variant.variant_class {
        VariantClass::Insertion => normalized_hgvs_alt_allele(variant, transcript, shifted)
            .unwrap_or_else(|| alt_allele.to_vec()),
        _ => alt_allele.to_vec(),
    });
    let fallback_ref = in_transcript_orientation(variant.ref_allele.clone());

    // OutputFactory.pm asks for HGVS only `within_feature`, the raw overlap of the
    // variant with the transcript, so an insertion butting against a transcript
    // edge has no notation even though the clamp below would place it. This gate
    // also covers the `return undef` of `_var2transcript_slice_coords`
    // (TranscriptVariationAllele.pm 2729) for a span wholly before or wholly past
    // the transcript.
    if variant.end < transcript.start || variant.start > transcript.end {
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
    // `_var2transcript_slice_coords`: on the reverse strand the ends swap.
    let (raw_start, raw_end) = if reverse {
        (slice.position(vf_end), slice.position(vf_start))
    } else {
        (slice.position(vf_start), slice.position(vf_end))
    };
    // `_var2transcript_slice_coords` (TranscriptVariationAllele.pm 2732-2733): an
    // end outside the slice is clamped onto its edge, so a deletion running over
    // an edge is described over the bases inside the transcript. `clamp` panics on
    // an empty slice, which Perl never checks.
    let tr_len = slice.len();
    if tr_len < 1 {
        return None;
    }
    let slice_start = raw_start.clamp(1, tr_len);
    let slice_end = raw_end.clamp(1, tr_len);
    if tr_len < slice_end + offset {
        return None;
    }
    // Without a FASTA the reference bases come from the variant: the bases the
    // clamp cut off are dropped from that fallback too.
    let front_cut = (slice_start - raw_start).max(0) as usize;
    let back_cut = (raw_end - slice_end).max(0) as usize;
    let fallback_ref = fallback_ref
        .get(front_cut..fallback_ref.len().saturating_sub(back_cut))
        .unwrap_or(&fallback_ref);

    let mut n = hgvs_variant_notation(
        |start, len| slice.substr(start, len),
        slice.fasta.is_none().then_some(fallback_ref),
        alt_seq,
        slice_start + offset,
        slice_end + offset,
    )?;
    if n.kind != HgvscKind::Dup {
        hgvsc_clip_alleles(&mut n);
    }

    let ref_name = versioned_reference_name(&transcript.stable_id, transcript.version);

    let coding = (transcript.translation.is_some() && mapper.cdna_coding_start > 0).then_some((
        mapper.cdna_coding_start as i64,
        mapper.cdna_coding_end as i64,
    ));
    let exons = ExonSpan::in_genomic_order(pairs);

    let same_pos = n.start == n.end;
    let is_snp = variant.ref_allele.len() == 1
        && alt_allele.len() == 1
        && variant.ref_allele != b"-"
        && alt_allele != b"-";
    let exonic = exons
        .iter()
        .any(|e| vf_start.min(vf_end) <= e.end && vf_start.max(vf_end) >= e.start);
    let cds_snp_position = if is_snp && exonic && coding.is_some() {
        crate::coding::perl_span(transcript, variant.start, variant.end)
            .and_then(|s| s.cds_start.zip(s.cds_end))
            .map(|(cds_start, _)| cds_start)
    } else {
        None
    };
    let (start, end) = match cds_snp_position {
        Some(cds_start) => (cds_start.to_string(), cds_start.to_string()),
        None => {
            let start = perl_cdna_position(&slice, &exons, coding, n.start)?;
            let end = if same_pos {
                start.clone()
            } else {
                perl_cdna_position(&slice, &exons, coding, n.end)?
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
    let (start, end) = if (exon_start > exon_end
        || (exon_start == exon_end && intron_start > intron_end))
        && !end.contains('*')
    {
        (end, start)
    } else {
        (start, end)
    };

    let numbering = if coding.is_some() { 'c' } else { 'n' };
    Some(format_hgvs_string(&ref_name, numbering, &start, &end, &n))
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
    if variant.variant_class != VariantClass::Insertion {
        return None;
    }

    let shifted_start = shifted?.0;
    let shift_steps = shifted_start.abs_diff(variant.start) as usize;
    if shift_steps == 0 {
        return None;
    }

    let alt_allele = variant.alt_allele();
    if alt_allele == b"-" || alt_allele.is_empty() {
        return None;
    }

    Some(rotate_sequence_for_transcript_3prime(
        alt_allele,
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
/// alternate allele carries a character outside `ACGT-` or equals the
/// reference, or Perl would return `undef` for the span. An insertion or
/// deletion is first shifted to its most 3' position on the transcript strand
/// (Perl's `_return_3prime(1)`, applied for HGVS whatever `--shift_3prime`
/// says) when a reference FASTA is available; the shifted alleles are the
/// deleted reference bases at the new position and the rotated insertion.
pub fn generate_hgvsp(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<String> {
    let shifted = hgvs_shift(variant, transcript, reference_fasta);
    hgvsp_at(variant, transcript, reference_fasta, shifted)
}

/// `generate_hgvsp` with the 3' shift already computed.
fn hgvsp_at(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
    shifted: Option<(u64, u64)>,
) -> Option<String> {
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let has_coding_model = transcript.translation.is_some()
        && mapper.cdna_coding_start > 0
        && mapper.cdna_coding_end >= mapper.cdna_coding_start;
    if !has_coding_model {
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

    let protein_id = transcript.protein_id.as_ref()?;
    let protein_prefix = format!(
        "{}:p.",
        versioned_reference_name(
            protein_id,
            transcript.translation.as_ref().and_then(|t| t.version)
        )
    );

    let shifted_variant = shifted.map(|(start, end)| {
        let ref_allele = normalized_hgvs_ref_allele(variant, shifted, reference_fasta)
            .unwrap_or_else(|| b"-".to_vec());
        let alt_allele = normalized_hgvs_alt_allele(variant, transcript, shifted)
            .unwrap_or_else(|| b"-".to_vec());
        InputVariant::new(variant.chr.clone(), start, end, ref_allele, alt_allele)
    });
    let notation = crate::coding::perl_hgvs_protein(
        variant,
        shifted_variant.as_ref().map(|v| (v, v.start, v.end)),
        transcript,
        reference_fasta,
    )?;
    Some(format!("{protein_prefix}{notation}"))
}

/// HGVSg, the genomic notation `--hgvsg` adds: Perl VEP's `hgvs_genomic`
/// (VariationFeature.pm 1877) on the variant's chromosome, as
/// `VariationFeatureOverlapAllele_to_output_hash` (OutputFactory.pm 1366-1395)
/// calls it for every sequence-variant row.
///
/// The reference name is the chromosome as the input names it (`21:g.`,
/// `chr21:g.`) and the numbering `m.` on `MT`. The alleles are those of the
/// whole record ([`InputVariant::record_allele_string`]); a record with several
/// ALTs describes the ALT at `ALLELE_NUM`, trimmed against the reference
/// (`trim_sequences`) to find its span and class, a bi-allelic record its one ALT
/// as written. An insertion or deletion is moved to its most 3' position on the
/// forward strand, whatever the transcript strand, within the flank Perl reads
/// (100 bases, more for an allele over 50; `_get_flank_seq`, VariationFeature.pm
/// 1805); the reference bases come from the FASTA, never from the input, so the
/// type follows `hgvs_variant_notation` on the sequence there and a `delins` is
/// clipped by `_clip_alleles` (VariationFeature.pm 2046).
///
/// `None` without a FASTA or when the FASTA lacks the chromosome, for a
/// structural variant (whose rows carry no HGVSg), for an allele with a
/// character outside `ACGT-`, for a span outside the chromosome, and where the
/// allele equals the reference sequence at its (shifted) position.
pub fn generate_hgvsg(
    variant: &InputVariant,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<String> {
    if variant.is_structural || variant.variant_class.is_structural() || variant.no_allele {
        return None;
    }
    let fasta = reference_fasta?;
    let (contig, chr_len) = fasta.contig(&variant.chr)?;
    let chr_len = chr_len as i64;

    // The VariationFeature's allele string, upper-cased by `validate_vf`.
    let record = variant.record_allele_string().to_ascii_uppercase();
    let all_alleles: Vec<&str> = record.split('/').collect();
    let is_multi_allelic = all_alleles.len() > 2;
    // OutputFactory.pm 1387-1390: the ALT at `allele_number` of a multi-allelic
    // record, else the second allele.
    let allele_number = if is_multi_allelic {
        variant
            .allele_num
            .map_or(variant.allele_index + 1, |n| n as usize)
    } else {
        1
    };
    let original_allele = *all_alleles.get(allele_number)?;

    let (vf_start, vf_end) = (variant.start as i64, variant.end as i64);
    // VariationFeature.pm 1904-1907: the span must lie on the chromosome.
    if vf_start < 1 || vf_end < 1 || vf_start > chr_len || vf_end > chr_len {
        return None;
    }

    // `_get_flank_seq` (VariationFeature.pm 1805-1849): 100 bases each side, or
    // 50 more than the length of the last allele longer than 50.
    let mut add_length: i64 = 100;
    for allele in &all_alleles {
        if allele.len() > 50 {
            add_length = allele.len() as i64 + 50;
        }
    }
    let mut seq_start = vf_start - add_length;
    let seq_end = vf_end + add_length;
    let mut ref_start = add_length;
    let mut ref_end = add_length + vf_end - vf_start;
    if seq_start < 0 {
        ref_start += seq_start;
        ref_end += seq_start;
        seq_start = 0;
    }
    let flank = padded_sequence(fasta, &variant.chr, chr_len, seq_start + 1, seq_end)?;

    // VariationFeature.pm 1942-1944: written tandems spelt out, then any allele
    // with a character outside `ACGT-` skipped.
    let allele = expand_tandems(original_allele).into_bytes();
    if allele
        .iter()
        .any(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'-'))
    {
        return None;
    }
    // VariationFeature.pm 1949-1957: an Ensembl-format line on the reverse strand
    // is described on the forward strand of the chromosome slice.
    let flip = variant.strand == Strand::Reverse;

    // VariationFeature.pm 1967-1976: an ALT of a multi-allelic record is trimmed
    // against the reference, narrowing the span. Perl took its `$check_allele`
    // before this, so an allele that is not an insertion or deletion is still
    // described by its untrimmed bases over the narrowed span (`GT/G/CT` gives
    // `delinsCT` at the `G`).
    let mut chr_start = vf_start;
    let (trimmed_ref, trimmed_alt): (Cow<[u8]>, Cow<[u8]>) = if is_multi_allelic {
        let (trimmed_ref, trimmed_alt, start, end) =
            trim_sequences(all_alleles[0].as_bytes(), &allele, vf_start, vf_end);
        ref_start += start - vf_start;
        ref_end += end - vf_end;
        chr_start = start;
        (
            dash_if_empty(trimmed_ref).into(),
            dash_if_empty(trimmed_alt).into(),
        )
    } else {
        (all_alleles[0].as_bytes().into(), allele.as_slice().into())
    };

    // `var_allele_class` (SO_variation_class, Utils/Sequence.pm 335): only an
    // `insertion` (`-/ACG`) or a `deletion` (`ACG/-`) is shifted.
    let base_run = |s: &[u8]| !s.is_empty() && s.iter().all(u8::is_ascii_uppercase);
    let insertion = *trimmed_ref == *b"-" && base_run(&trimmed_alt);
    let deletion = base_run(&trimmed_ref) && trimmed_alt.contains(&b'-');

    // VariationFeature.pm 1981-2012 (`DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME` is 1):
    // the deleted or inserted bases rotate along the flank after the variant.
    let (check_allele, offset) = if insertion || deletion {
        let mut seq_to_check = if deletion { trimmed_ref } else { trimmed_alt }.into_owned();
        if flip {
            seq_to_check = crate::coding::reverse_complement(&seq_to_check);
        }
        let (three_prime_allele, offset) = match flank.get(ref_end as usize..) {
            Some(downstream) => three_prime_seq_offset(seq_to_check, downstream),
            None => (seq_to_check, 0),
        };
        let check_allele = if insertion {
            three_prime_allele
        } else {
            b"-".to_vec()
        };
        (check_allele, offset)
    } else if flip {
        // Perl's `lookup_order` -1 matters only where the span is empty, an
        // insertion, which always takes the branch above.
        (crate::coding::reverse_complement(&allele), 0)
    } else {
        (allele, 0)
    };

    let alt_seq: Vec<u8> = check_allele.into_iter().filter(|&b| b != b'-').collect();
    let mut n = hgvs_variant_notation(
        |start, len| perl_substr(&flank, start - 1, len).map(<[u8]>::to_vec),
        None,
        alt_seq,
        ref_start + offset,
        ref_end + offset,
    )?;
    // The flank positions the notation carries become chromosome positions
    // (`$chr_start + $offset` and `$chr_end + $offset` are the display span).
    let to_chromosome = chr_start - ref_start;
    n.start += to_chromosome;
    n.end += to_chromosome;
    if n.kind == HgvscKind::Delins {
        genomic_clip_alleles(&mut n);
    }

    // VariationFeature.pm 2034: `m.` only when the slice is named `MT`; offline the
    // slice is the FASTA contig the input chromosome resolves to (BaseVEP.pm
    // `get_slice`), so a UCSC `chrM` FASTA gives `g.`.
    let numbering = if contig == "MT" { 'm' } else { 'g' };
    Some(format_hgvs_string(
        &variant.original_chr,
        numbering,
        &n.start.to_string(),
        &n.end.to_string(),
        &n,
    ))
}

/// Bases `start..=end` (1-based) of `chr`, read as Perl's FASTA-backed slice
/// reads them: positions past the end of the contig are `N` (`_raw_seq`,
/// FastaSequence.pm 627-636).
fn padded_sequence(
    fasta: &vep_fasta::IndexedFasta,
    chr: &str,
    chr_len: i64,
    start: i64,
    end: i64,
) -> Option<Vec<u8>> {
    let mut seq = fasta.sequence(chr, start as u64, end.min(chr_len) as u64)?;
    seq.resize((end - start + 1) as usize, b'N');
    Some(seq)
}

/// Perl's `substr($s, $offset, $len)`: a negative offset counts from the end,
/// a span partly outside the string yields the part inside it, and an offset
/// past the end yields `undef`.
fn perl_substr(s: &[u8], offset: i64, len: usize) -> Option<&[u8]> {
    let n = s.len() as i64;
    let start = if offset < 0 { n + offset } else { offset };
    if start > n {
        return None;
    }
    let end = (start + len as i64).clamp(0, n);
    let start = start.clamp(0, n);
    Some(&s[start as usize..end.max(start) as usize])
}

/// `expand` (ensembl Utils/Sequence.pm): `s/(\w*)\((\w+)\)(\d+)/$1.$2 x $3/eg`,
/// so `(CA)3` reads `CACACA`.
fn expand_tandems(allele: &str) -> String {
    if !allele.contains('(') {
        return allele.to_string();
    }
    let is_word = |b: u8| b.is_ascii_alphanumeric() || b == b'_';
    let bytes = allele.as_bytes();
    let mut out = String::with_capacity(allele.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'(' {
            let unit_end = (i + 1..bytes.len())
                .find(|&j| !is_word(bytes[j]))
                .unwrap_or(bytes.len());
            let count_end = (unit_end + 1..bytes.len())
                .find(|&j| !bytes[j].is_ascii_digit())
                .unwrap_or(bytes.len());
            if unit_end > i + 1 && bytes.get(unit_end) == Some(&b')') && count_end > unit_end + 1 {
                let count: usize = allele[unit_end + 1..count_end].parse().unwrap_or(0);
                out.push_str(&allele[i + 1..unit_end].repeat(count));
                i = count_end;
                continue;
            }
        }
        out.push(bytes[i] as char);
        i += 1;
    }
    out
}

/// `trim_sequences` (Utils/Sequence.pm 965) on the forward strand: the bases
/// `ref_allele` and `alt` share are trimmed from the left, moving `start`, then
/// from the right, moving `end`.
fn trim_sequences(
    ref_allele: &[u8],
    alt: &[u8],
    mut start: i64,
    mut end: i64,
) -> (Vec<u8>, Vec<u8>, i64, i64) {
    let (mut r, mut a) = (ref_allele, alt);
    while let (Some(x), Some(y)) = (r.first(), a.first()) {
        if x != y {
            break;
        }
        (r, a) = (&r[1..], &a[1..]);
        start += 1;
    }
    while let (Some(x), Some(y)) = (r.last(), a.last()) {
        if x != y {
            break;
        }
        (r, a) = (&r[..r.len() - 1], &a[..a.len() - 1]);
        end -= 1;
    }
    (r.to_vec(), a.to_vec(), start, end)
}

fn dash_if_empty(allele: Vec<u8>) -> Vec<u8> {
    if allele.is_empty() {
        b"-".to_vec()
    } else {
        allele
    }
}

/// `get_3prime_seq_offset` (Utils/Sequence.pm 792): `seq` rotated one base at a
/// time while its first base equals the next base of `post_seq`, and the number
/// of positions it moved.
fn three_prime_seq_offset(mut seq: Vec<u8>, post_seq: &[u8]) -> (Vec<u8>, i64) {
    let mut offset = 0;
    for &next in post_seq {
        if seq.first() != Some(&next) {
            break;
        }
        seq.rotate_left(1);
        offset += 1;
    }
    (seq, offset)
}

/// `_clip_alleles` of VariationFeature.pm (2046), which `hgvs_genomic` applies
/// to a `delins`: the shared leading and trailing bases are trimmed, moving
/// `start` and `end`, and the type is re-read from what remains. An emptied
/// reference becomes a `dup` of the trimmed leading bases when they end in the
/// inserted ones, else an `ins` between the base before and the one after; an
/// emptied alternate a `del`. The `=` type is unreachable: a `delins` has
/// differing alleles.
fn genomic_clip_alleles(n: &mut HgvscNotation) {
    let mut preseq: Vec<u8> = Vec::new();
    let (mut r, mut a): (&[u8], &[u8]) = (&n.ref_seq, &n.alt_seq);
    while let (Some(&x), Some(&y)) = (r.first(), a.first()) {
        if x != y {
            break;
        }
        preseq.push(x);
        (r, a) = (&r[1..], &a[1..]);
        n.start += 1;
    }
    while let (Some(x), Some(y)) = (r.last(), a.last()) {
        if x != y {
            break;
        }
        (r, a) = (&r[..r.len() - 1], &a[..a.len() - 1]);
        n.end -= 1;
    }
    let (r, a) = (dash_if_empty(r.to_vec()), dash_if_empty(a.to_vec()));
    if r != b"-" && r.len() == 1 && a.len() == 1 {
        n.kind = HgvscKind::Sub;
    } else if r == b"-" {
        if preseq.ends_with(&a) {
            n.kind = HgvscKind::Dup;
            n.start -= a.len() as i64;
        } else {
            n.kind = HgvscKind::Ins;
            n.start -= 1;
            n.end = n.start + 1;
        }
    } else if a == b"-" {
        n.kind = HgvscKind::Del;
    }
    n.ref_seq = r;
    n.alt_seq = a;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::make_test_transcript;

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
        // No reference sequence: the reference allele is the fallback.
        let unreadable = |_: i64, _: usize| None;
        let kind = |r: &[u8], a: &[u8]| {
            let end = 100 + r.len() as i64 - 1;
            hgvs_variant_notation(unreadable, Some(r), a.to_vec(), 100, end)
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
            hgvs_variant_notation(unreadable, Some(b""), b"TT".to_vec(), 100, 99)
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
    /// by its peptides: the two bases left of `TAA` translate to `X`,
    /// three-lettered `Xaa` and rewritten `Ter`, equal to the reference `Ter`, so
    /// the notation is `Ter283=`. Rows of this shape in
    /// `tests/golden/116/GRCh37-hgvs`: `21:34955948 -` against ENST00000457359,
    /// `ENSP00000415191.1:p.Ter211=`, and `21:35791518 -` against ENST00000450895,
    /// `ENSP00000397039.1:p.Ter140=`.
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
            Some("ENSP00000000001.1:p.Ter283=")
        );
    }

    /// Chromosome `21` of the HGVSg fixture: an `AAAA` run at 11-14, a `CACACA`
    /// run at 21-26, `GGATC` at 31-35 and `CC` at 35-36. Every expected HGVSg
    /// below is what Ensembl VEP 116's `hgvs_genomic` prints for the same variant
    /// on this sequence (a FASTA-backed whole-chromosome slice, as `--offline
    /// --fasta` builds it).
    const HGVSG_CHR21: &str = "CGTACGTACGAAAATCCGATCACACAGTTCGGATCCTAGGTCAGAGACTTGCATGCATGC";

    /// The fixture FASTA: `21` as above and a 14-base `MT`.
    fn hgvsg_fasta() -> (tempfile::TempDir, vep_fasta::IndexedFasta) {
        let dir = tempfile::tempdir().expect("tempdir");
        let fa = dir.path().join("ref.fa");
        std::fs::write(&fa, format!(">21\n{HGVSG_CHR21}\n>MT\nGATTACAGATTACA\n"))
            .expect("write FASTA");
        std::fs::write(
            dir.path().join("ref.fa.fai"),
            "21\t60\t4\t60\t61\nMT\t14\t69\t14\t15\n",
        )
        .expect("write FASTA index");
        let fasta = vep_fasta::IndexedFasta::from_path(&fa).expect("load indexed FASTA");
        (dir, fasta)
    }

    /// A bi-allelic variant as the parsers hold it (an insertion has `end` one
    /// below `start`).
    fn biallelic(chr: &str, start: u64, end: u64, ref_allele: &str, alt: &str) -> InputVariant {
        InputVariant::new(
            chr.into(),
            start,
            end,
            ref_allele.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        )
    }

    fn hgvsg(variant: &InputVariant, fasta: &vep_fasta::IndexedFasta) -> Option<String> {
        generate_hgvsg(variant, Some(fasta))
    }

    #[test]
    fn test_hgvsg_snv_is_named_by_the_input_chromosome() {
        let (_dir, fasta) = hgvsg_fasta();
        assert_eq!(
            hgvsg(&biallelic("21", 33, 33, "A", "G"), &fasta).as_deref(),
            Some("21:g.33A>G")
        );
        assert_eq!(
            hgvsg(&biallelic("chr21", 33, 33, "A", "G"), &fasta).as_deref(),
            Some("chr21:g.33A>G")
        );
    }

    /// The reference bases come from the FASTA, not the input: a mismatching REF
    /// is described with the FASTA's base, and an ALT equal to it has no HGVSg.
    #[test]
    fn test_hgvsg_reads_the_reference_from_the_fasta() {
        let (_dir, fasta) = hgvsg_fasta();
        assert_eq!(
            hgvsg(&biallelic("21", 33, 33, "C", "G"), &fasta).as_deref(),
            Some("21:g.33A>G")
        );
        assert_eq!(hgvsg(&biallelic("21", 33, 33, "G", "A"), &fasta), None);
    }

    #[test]
    fn test_hgvsg_deletion_outside_a_repeat_keeps_its_position() {
        let (_dir, fasta) = hgvsg_fasta();
        assert_eq!(
            hgvsg(&biallelic("21", 34, 34, "T", "-"), &fasta).as_deref(),
            Some("21:g.34del")
        );
    }

    /// A deletion in a repeat moves to the repeat's 3' end on the forward strand.
    #[test]
    fn test_hgvsg_deletion_in_a_repeat_is_shifted_3prime() {
        let (_dir, fasta) = hgvsg_fasta();
        assert_eq!(
            hgvsg(&biallelic("21", 11, 11, "A", "-"), &fasta).as_deref(),
            Some("21:g.14del")
        );
        assert_eq!(
            hgvsg(&biallelic("21", 21, 22, "CA", "-"), &fasta).as_deref(),
            Some("21:g.25_26del")
        );
    }

    #[test]
    fn test_hgvsg_insertion_repeating_no_neighbour_is_ins() {
        let (_dir, fasta) = hgvsg_fasta();
        assert_eq!(
            hgvsg(&biallelic("21", 35, 34, "-", "G"), &fasta).as_deref(),
            Some("21:g.34_35insG")
        );
    }

    /// An insertion of bases the reference already holds is a `dup` of the
    /// most 3' copy: in place after the `AAAA` run, and after shifting from the
    /// run's start or the `CACACA` run's start.
    #[test]
    fn test_hgvsg_insertion_repeating_its_neighbour_is_dup() {
        let (_dir, fasta) = hgvsg_fasta();
        assert_eq!(
            hgvsg(&biallelic("21", 15, 14, "-", "A"), &fasta).as_deref(),
            Some("21:g.14dup")
        );
        assert_eq!(
            hgvsg(&biallelic("21", 11, 10, "-", "A"), &fasta).as_deref(),
            Some("21:g.14dup")
        );
        assert_eq!(
            hgvsg(&biallelic("21", 21, 20, "-", "CA"), &fasta).as_deref(),
            Some("21:g.25_26dup")
        );
    }

    /// Equal-length alleles: a `delins`, an `inv` when the ALT is the reverse
    /// complement, and a `delins` clipped of the bases it shares with the
    /// reference (`_clip_alleles`).
    #[test]
    fn test_hgvsg_mnv_is_delins_or_inv() {
        let (_dir, fasta) = hgvsg_fasta();
        assert_eq!(
            hgvsg(&biallelic("21", 34, 35, "TC", "AG"), &fasta).as_deref(),
            Some("21:g.34_35delinsAG")
        );
        assert_eq!(
            hgvsg(&biallelic("21", 31, 33, "GGA", "TCC"), &fasta).as_deref(),
            Some("21:g.31_33inv")
        );
        assert_eq!(
            hgvsg(&biallelic("21", 31, 33, "GGA", "GCT"), &fasta).as_deref(),
            Some("21:g.32_33delinsCT")
        );
    }

    /// The record `TC/T/GC/TCC` at 34-35 (VCF `21 34 . TC T,GC,TCC`): each ALT,
    /// picked by its ALLELE_NUM, is trimmed against the reference first. `T` is a
    /// deletion of the `C` at 35, shifted to 36; `TCC` an insertion of `C`, a
    /// `dup` of 36; and `GC`, neither, keeps its untrimmed bases over the trimmed
    /// one-base span, as Perl does.
    #[test]
    fn test_hgvsg_multi_allelic_record_describes_each_alt() {
        let (_dir, fasta) = hgvsg_fasta();
        let record: std::sync::Arc<str> = "TC/T/GC/TCC".into();
        let expected = ["21:g.36del", "21:g.34delinsGC", "21:g.36dup"];
        for (index, (alt, expected)) in ["T", "GC", "TCC"].into_iter().zip(expected).enumerate() {
            let mut variant = biallelic("21", 34, 35, "TC", alt);
            variant.record_allele_string_multi = Some(record.clone());
            variant.allele_index = index;
            assert_eq!(hgvsg(&variant, &fasta).as_deref(), Some(expected), "{alt}");
        }
    }

    #[test]
    fn test_hgvsg_mitochondrion_uses_m_numbering() {
        let (_dir, fasta) = hgvsg_fasta();
        assert_eq!(
            hgvsg(&biallelic("MT", 2, 2, "A", "C"), &fasta).as_deref(),
            Some("MT:m.2A>C")
        );
        assert_eq!(
            hgvsg(&biallelic("chrM", 2, 2, "A", "C"), &fasta).as_deref(),
            Some("chrM:m.2A>C")
        );
    }

    /// Ensembl VEP picks `m.` by the FASTA contig's name, not the input's: a UCSC
    /// FASTA's `chrM` is numbered `g.` (Ensembl VEP 116 prints `chrM:g.2A>C`).
    #[test]
    fn test_hgvsg_ucsc_mitochondrion_uses_g_numbering() {
        let dir = tempfile::tempdir().expect("tempdir");
        let fa = dir.path().join("ucsc.fa");
        std::fs::write(&fa, ">chrM\nGATTACAGATTACA\n").expect("write FASTA");
        std::fs::write(dir.path().join("ucsc.fa.fai"), "chrM\t14\t6\t14\t15\n")
            .expect("write FASTA index");
        let fasta = vep_fasta::IndexedFasta::from_path(&fa).expect("load indexed FASTA");
        assert_eq!(
            hgvsg(&biallelic("chrM", 2, 2, "A", "C"), &fasta).as_deref(),
            Some("chrM:g.2A>C")
        );
    }

    /// The 3' shift stops at the flank Ensembl VEP reads: 100 bases, or 50 more
    /// than an allele longer than 50 bases. Contig `1` is `CG` x 100, then 150 `A`s
    /// (201-350), then 50 `C`s. Ensembl VEP 116.2 prints `1:g.301del`,
    /// `1:g.291_350del` and `1:g.251_305dup` for these records.
    #[test]
    fn test_hgvsg_shift_stops_at_the_flank() {
        let dir = tempfile::tempdir().expect("tempdir");
        let fa = dir.path().join("flank.fa");
        let seq = format!("{}{}{}", "CG".repeat(100), "A".repeat(150), "C".repeat(50));
        std::fs::write(&fa, format!(">1\n{seq}\n")).expect("write FASTA");
        std::fs::write(dir.path().join("flank.fa.fai"), "1\t400\t3\t400\t401\n")
            .expect("write FASTA index");
        let fasta = vep_fasta::IndexedFasta::from_path(&fa).expect("load indexed FASTA");
        // VCF `1 200 GA G`
        assert_eq!(
            hgvsg(&biallelic("1", 201, 201, "A", "-"), &fasta).as_deref(),
            Some("1:g.301del")
        );
        // VCF `1 200 G(A x 60) G`
        assert_eq!(
            hgvsg(&biallelic("1", 201, 260, &"A".repeat(60), "-"), &fasta).as_deref(),
            Some("1:g.291_350del")
        );
        // VCF `1 200 G G(A x 55)`
        assert_eq!(
            hgvsg(&biallelic("1", 201, 200, "-", &"A".repeat(55)), &fasta).as_deref(),
            Some("1:g.251_305dup")
        );
    }

    /// A deletion shifted up to the contig's last base (the flank past the end
    /// reads as `N`), and Ensembl VEP's `0_1ins` for an insertion before base 1
    /// in a multi-allelic record (`22:g.17del`, `21:g.0_1insT`, `21:g.1C>G` in
    /// Ensembl VEP 116.2).
    #[test]
    fn test_hgvsg_at_the_contig_ends() {
        let dir = tempfile::tempdir().expect("tempdir");
        let fa = dir.path().join("ends.fa");
        std::fs::write(&fa, ">22\nACGTACGTACGTAAAAA\n").expect("write FASTA");
        std::fs::write(dir.path().join("ends.fa.fai"), "22\t17\t4\t17\t18\n")
            .expect("write FASTA index");
        let ends = vep_fasta::IndexedFasta::from_path(&fa).expect("load indexed FASTA");
        // VCF `22 12 TA T`
        assert_eq!(
            hgvsg(&biallelic("22", 13, 13, "A", "-"), &ends).as_deref(),
            Some("22:g.17del")
        );

        // VCF `21 1 C TC,G`
        let (_dir, fasta) = hgvsg_fasta();
        let record: std::sync::Arc<str> = "C/TC/G".into();
        for (index, (alt, expected)) in [("TC", "21:g.0_1insT"), ("G", "21:g.1C>G")]
            .into_iter()
            .enumerate()
        {
            let mut variant = biallelic("21", 1, 1, "C", alt);
            variant.record_allele_string_multi = Some(record.clone());
            variant.allele_index = index;
            assert_eq!(hgvsg(&variant, &fasta).as_deref(), Some(expected), "{alt}");
        }
    }

    /// An Ensembl-format line on the reverse strand is described on the forward
    /// strand: its alleles are complemented, and a deletion is shifted only when
    /// its complemented bases repeat downstream.
    #[test]
    fn test_hgvsg_reverse_strand_input_is_described_on_the_forward_strand() {
        let (_dir, fasta) = hgvsg_fasta();
        let mut snv = biallelic("21", 34, 34, "T", "G");
        snv.strand = Strand::Reverse;
        assert_eq!(hgvsg(&snv, &fasta).as_deref(), Some("21:g.34T>C"));
        let mut deletion = biallelic("21", 21, 22, "CA", "-");
        deletion.strand = Strand::Reverse;
        assert_eq!(hgvsg(&deletion, &fasta).as_deref(), Some("21:g.21_22del"));
    }

    /// Ensembl VEP leaves HGVSg empty without a reference sequence, for an ALT
    /// outside `ACGT-`, for a span off the chromosome and for a structural
    /// variant; vep-rs leaves it empty too when the FASTA lacks the chromosome.
    #[test]
    fn test_hgvsg_is_empty_where_vep_leaves_it_empty() {
        let (_dir, fasta) = hgvsg_fasta();
        assert_eq!(
            generate_hgvsg(&biallelic("21", 33, 33, "A", "G"), None),
            None
        );
        assert_eq!(hgvsg(&biallelic("22", 33, 33, "A", "G"), &fasta), None);
        assert_eq!(hgvsg(&biallelic("21", 33, 33, "A", "N"), &fasta), None);
        assert_eq!(hgvsg(&biallelic("21", 61, 60, "-", "A"), &fasta), None);
        let mut structural = biallelic("21", 30, 40, "N", "<DEL>");
        structural.variant_class = VariantClass::StructuralDeletion;
        structural.is_structural = true;
        assert_eq!(hgvsg(&structural, &fasta), None);
    }
}
