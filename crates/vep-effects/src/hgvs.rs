// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! HGVS notation generation for variant-transcript consequences.
//!
//! Produces HGVSc (coding DNA) and HGVSp (protein) notation strings following
//! the Human Genome Variation Society nomenclature conventions, matching Perl
//! VEP's output.
//!
//! Reference: <https://varnomen.hgvs.org/>

use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;
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

/// How far Perl VEP's `perform_shift` (TranscriptVariationAllele.pm 290) can
/// move an insertion or deletion of `len` bases: it compares against 1,000 bases
/// of flank, so a pattern that fits the flank moves at most `1001 - len`
/// positions and a longer one at most 1,000, except that on the reverse strand
/// a 1,001-base pattern gets a loop bound of zero and never moves.
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

/// Perl's HGVS variant type (`hgvs_variant_notation`, Utils/Sequence.pm 493).
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

/// Perl's `$hgvs_notation` hash for a transcript-level description.
struct HgvscNotation {
    start: i64,
    end: i64,
    ref_seq: Vec<u8>,
    alt_seq: Vec<u8>,
    kind: HgvscKind,
}

/// `hgvs_variant_notation` (Utils/Sequence.pm 493) with the default lookup
/// order: the type and displayed span of `alt_seq` replacing slice positions
/// `ref_start..=ref_end`. `None` where Perl returns `undef` (the alleles are
/// equal) or the reference bases cannot be read.
fn hgvs_variant_notation(
    slice: &TranscriptSlice<'_>,
    fallback_ref: &[u8],
    alt_seq: Vec<u8>,
    ref_start: i64,
    ref_end: i64,
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
        // The bases after the site first, then the bases before it: a match is a
        // duplication of those bases.
        let after = slice.substr(ref_end + 1, alt_length);
        if after.as_deref() == Some(n.alt_seq.as_slice()) {
            n.end = n.start + alt_length as i64 - 1;
            n.kind = HgvscKind::Dup;
            return Some(n);
        }
        let before = slice.substr(ref_end - alt_length as i64 + 1, alt_length);
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

/// `_clip_alleles` (TranscriptVariationAllele.pm 2118) as `hgvs_transcript`
/// calls it, before `numbering` is set: the shared residues are trimmed from the
/// front and the back, moving `start` and `end`, and the type is re-read from
/// what remains; the `=` and `dup` re-sets, which need `numbering`, never fire.
fn hgvsc_clip_alleles(n: &mut HgvscNotation) {
    let mut ref_seq: &[u8] = &n.ref_seq;
    let mut alt_seq: &[u8] = &n.alt_seq;
    for _ in 0..n.ref_seq.len() {
        match (ref_seq.first(), alt_seq.first()) {
            (Some(r), Some(a)) if r == a => {
                n.start += 1;
                ref_seq = &ref_seq[1..];
                alt_seq = &alt_seq[1..];
            }
            _ => break,
        }
    }
    for _ in 0..ref_seq.len() {
        match (ref_seq.last(), alt_seq.last()) {
            (Some(r), Some(a)) if r == a => {
                ref_seq = &ref_seq[..ref_seq.len() - 1];
                alt_seq = &alt_seq[..alt_seq.len() - 1];
                n.end -= 1;
            }
            _ => break,
        }
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
    let (ref_seq, alt_seq) = (ref_seq.to_vec(), alt_seq.to_vec());
    n.ref_seq = ref_seq;
    n.alt_seq = alt_seq;
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

/// `_get_cDNA_position` (TranscriptVariationAllele.pm 2683): the HGVS cDNA
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
            // Perl clears the coordinate, prefixes `*` and strips the `+` from the
            // offset, so the stop codon's last base plus an intronic distance
            // prints as `*` followed by the bare distance.
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
    let body = match &n.kind {
        HgvscKind::Sub => format!("{start}{}>{alt}", String::from_utf8_lossy(&n.ref_seq)),
        HgvscKind::Inv if n.ref_seq.len() == 1 => {
            format!("{start}{}>{alt}", String::from_utf8_lossy(&n.ref_seq))
        }
        HgvscKind::Del => format!("{coordinates}del"),
        HgvscKind::Inv => format!("{coordinates}inv"),
        HgvscKind::Dup => format!("{coordinates}dup"),
        HgvscKind::Delins => format!("{coordinates}delins{alt}"),
        HgvscKind::Ins => format!("{coordinates}ins{alt}"),
        HgvscKind::Multiple(m) => format!("{coordinates}[{m}]"),
    };
    format!("{ref_name}:{numbering}.{body}")
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
/// `hgvs_transcript` (TranscriptVariationAllele.pm 1311) for an Ensembl
/// transcript.
///
/// Returns `None` when the alternate allele carries a character outside
/// `ACGT-` or equals the reference, when the variant lies outside the
/// transcript's genomic span (so an upstream or downstream variant has no HGVSc),
/// when the 3'-shifted span leaves it, or when the alleles agree at the shifted
/// position. An insertion or deletion is described at its most 3' position on
/// the transcript strand when a FASTA is present; the type and the displayed
/// span follow `hgvs_variant_notation`, shared bases are clipped, an exonic
/// single-base substitution inside the CDS takes its CDS coordinate and every
/// other position goes through `_get_cDNA_position`. Without a FASTA the
/// reference bases come from the variant and no duplication is detected.
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

    // The alternate sequence in transcript orientation: the rotated allele of a
    // shifted insertion, dashes removed, complemented on the reverse strand.
    let mut alt_seq: Vec<u8> = match variant.variant_class {
        VariantClass::Insertion => normalized_hgvs_alt_allele(variant, transcript, shifted)
            .unwrap_or_else(|| alt_allele.to_ascii_uppercase()),
        _ => alt_allele.to_ascii_uppercase(),
    };
    alt_seq.retain(|&b| b != b'-');
    if reverse {
        alt_seq = crate::coding::reverse_complement(&alt_seq);
    }
    let mut fallback_ref: Vec<u8> = variant.ref_allele.to_ascii_uppercase();
    fallback_ref.retain(|&b| b != b'-');
    if reverse {
        fallback_ref = crate::coding::reverse_complement(&fallback_ref);
    }

    let slice = TranscriptSlice {
        chr: &variant.chr,
        tr_start: transcript.start as i64,
        tr_end: transcript.end as i64,
        reverse,
        fasta: reference_fasta,
    };
    let (vf_start, vf_end) = (variant.start as i64, variant.end as i64);
    let (slice_start, slice_end) = if reverse {
        (slice.tr_end - vf_end + 1, slice.tr_end - vf_start + 1)
    } else {
        (vf_start - slice.tr_start + 1, vf_end - slice.tr_start + 1)
    };
    let tr_len = slice.len();
    if slice_start < 1 || slice_end < 1 || slice_start > tr_len || slice_end > tr_len {
        return None;
    }
    if tr_len < slice_end + offset {
        return None;
    }

    let mut n = hgvs_variant_notation(
        &slice,
        &fallback_ref,
        alt_seq,
        slice_start + offset,
        slice_end + offset,
    )?;
    if n.kind != HgvscKind::Dup {
        hgvsc_clip_alleles(&mut n);
    }

    let ref_name = match transcript.version {
        Some(version) if !transcript.stable_id.contains("LRG") => {
            format!("{}.{version}", transcript.stable_id)
        }
        _ => transcript.stable_id.to_string(),
    };

    let coding = (transcript.translation.is_some() && mapper.cdna_coding_start > 0).then_some((
        mapper.cdna_coding_start as i64,
        mapper.cdna_coding_end as i64,
    ));
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
    let protein_prefix = match transcript.translation.as_ref().and_then(|t| t.version) {
        Some(version) => format!("{protein_id}.{version}:p."),
        None => format!("{protein_id}:p."),
    };

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

    /// The fixture transcript as `_get_cDNA_position` sees it: the mapper's
    /// three exon pairs and the CDS bounds cDNA 51..900.
    fn fixture_slice_and_exons(tx: &Transcript) -> (TranscriptSlice<'static>, Vec<ExonSpan>) {
        let mapper = tx.vefc.as_ref().unwrap().mapper.as_ref().unwrap();
        let exons = mapper
            .exon_coord_mapper
            .pairs
            .iter()
            .map(|p| ExonSpan {
                start: p.to_start as i64,
                end: p.to_end as i64,
                cdna_start: p.from_start as i64,
                cdna_end: p.from_end as i64,
            })
            .collect();
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
    fn test_hgvs_variant_notation_types() {
        let slice = unreadable_slice();
        let kind = |r: &[u8], a: &[u8]| {
            let end = 100 + r.len() as i64 - 1;
            hgvs_variant_notation(&slice, r, a.to_vec(), 100, end).map(|n| (n.kind, n.start, n.end))
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
            hgvs_variant_notation(&slice, b"", b"TT".to_vec(), 100, 99)
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
}
