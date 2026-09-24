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

/// Generate HGVSc notation for a variant-transcript pair.
///
/// Returns notation like:
/// - `ENST00000366667.4:c.803T>C` (coding SNV)
/// - `ENST00000366667.4:c.803delT` (coding deletion)
/// - `ENST00000366667.4:c.803_804insAA` (coding insertion)
/// - `ENST00000366667.4:c.-14T>C` (5' UTR)
/// - `ENST00000366667.4:c.*42T>C` (3' UTR)
/// - `ENST00000366667.4:c.803+2T>C` (intronic near donor)
/// - `ENST00000366667.4:n.803T>C` (non-coding transcript)
pub fn generate_hgvsc(
    variant: &InputVariant,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<String> {
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let pairs = &mapper.exon_coord_mapper.pairs;
    if pairs.is_empty() {
        return None;
    }

    // Perl VEP 3'-shifts every indel for HGVS: to its most 3' position on the
    // transcript strand.
    let shifted = reference_fasta.and_then(|fasta| {
        crate::consequences::shift_indel_3prime_coords(variant, transcript, fasta, 2000)
    });
    let (var_start, var_end) = if let Some((s, e)) = shifted {
        (s, e)
    } else {
        (variant.start, variant.end)
    };

    let normalized_ref_allele = normalized_hgvs_ref_allele(variant, shifted, reference_fasta);
    let normalized_alt_allele = normalized_hgvs_alt_allele(variant, transcript, shifted);
    let ref_allele = normalized_ref_allele
        .as_deref()
        .unwrap_or(variant.ref_allele.as_slice());
    let alt_allele_raw = normalized_alt_allele
        .as_deref()
        .unwrap_or_else(|| variant.alt_allele());

    let has_coding_model = transcript.translation.is_some()
        && mapper.cdna_coding_start > 0
        && mapper.cdna_coding_end >= mapper.cdna_coding_start;

    let prefix_type = if has_coding_model { "c" } else { "n" };

    let transcript_prefix = if let Some(version) = transcript.version {
        format!("{}.{}:{}", transcript.stable_id, version, prefix_type)
    } else {
        format!("{}:{}", transcript.stable_id, prefix_type)
    };

    let cdna_start = genomic_to_cdna_pos(var_start, pairs);
    let cdna_end = if var_end != var_start {
        genomic_to_cdna_pos(var_end, pairs)
    } else {
        cdna_start
    };

    let start_intronic =
        cdna_start.is_none() && var_start >= transcript.start && var_start <= transcript.end;
    let end_intronic =
        cdna_end.is_none() && var_end >= transcript.start && var_end <= transcript.end;

    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele_raw == b"-" || alt_allele_raw.is_empty();

    // `Cow` avoids a heap allocation for the common "-" case.
    let (display_ref, display_alt): (std::borrow::Cow<'_, str>, std::borrow::Cow<'_, str>) =
        if transcript.strand == Strand::Reverse {
            (
                if ref_is_dash {
                    std::borrow::Cow::Borrowed("-")
                } else {
                    std::borrow::Cow::Owned(reverse_complement_string(ref_allele))
                },
                if alt_is_dash {
                    std::borrow::Cow::Borrowed("-")
                } else {
                    std::borrow::Cow::Owned(reverse_complement_string(alt_allele_raw))
                },
            )
        } else {
            (
                if ref_is_dash {
                    std::borrow::Cow::Borrowed("-")
                } else {
                    // Alleles are always ASCII nucleotides; use from_utf8 directly
                    std::borrow::Cow::Owned(
                        std::str::from_utf8(ref_allele)
                            .unwrap_or("-")
                            .to_uppercase(),
                    )
                },
                if alt_is_dash {
                    std::borrow::Cow::Borrowed("-")
                } else {
                    std::borrow::Cow::Owned(
                        std::str::from_utf8(alt_allele_raw)
                            .unwrap_or("-")
                            .to_uppercase(),
                    )
                },
            )
        };

    let is_insertion = variant.variant_class == VariantClass::Insertion;
    let is_deletion = variant.variant_class == VariantClass::Deletion;

    if start_intronic || end_intronic {
        if is_insertion {
            if let Some(dup_desc) =
                check_duplication_cdna(variant, var_start, var_end, transcript, reference_fasta)
            {
                return Some(format!("{transcript_prefix}{dup_desc}"));
            }
        }

        let ctx = IntronicHgvsContext {
            transcript,
            pairs,
            mapper,
            has_coding_model,
        };
        return generate_hgvsc_intronic_shifted(
            variant,
            var_start,
            var_end,
            &ctx,
            &transcript_prefix,
            &display_ref,
            &display_alt,
        );
    }

    let cdna_s = cdna_start?;
    let cdna_e = cdna_end.unwrap_or(cdna_s);

    let cdna_lo = cdna_s.min(cdna_e);
    let cdna_hi = cdna_s.max(cdna_e);

    let pos_start = format_cdna_position(cdna_lo, mapper, has_coding_model);
    let pos_end = if cdna_hi != cdna_lo {
        Some(format_cdna_position(cdna_hi, mapper, has_coding_model))
    } else {
        None
    };

    let variant_desc = if is_insertion {
        // Perl VEP writes an insertion between its two flanking positions. A
        // VEP-convention insertion has start > end, so cdna_s maps from the
        // higher genomic position; HGVS needs ascending order.
        let (ins_start, ins_end) = if cdna_s == cdna_e {
            // Both endpoints map to the same cDNA position, so use its flanks.
            let lo = cdna_s.saturating_sub(1).max(1);
            (
                format_cdna_position(lo, mapper, has_coding_model),
                format_cdna_position(cdna_s, mapper, has_coding_model),
            )
        } else {
            let lo = cdna_s.min(cdna_e);
            let hi = cdna_s.max(cdna_e);
            (
                format_cdna_position(lo, mapper, has_coding_model),
                format_cdna_position(hi, mapper, has_coding_model),
            )
        };

        if let Some(dup_desc) =
            check_duplication_cdna(variant, var_start, var_end, transcript, reference_fasta)
        {
            dup_desc
        } else {
            format!(".{ins_start}_{ins_end}ins{display_alt}")
        }
    } else if is_deletion {
        if let Some(ref pos_e) = pos_end {
            format!(".{pos_start}_{pos_e}del")
        } else {
            format!(".{pos_start}del")
        }
    } else if is_hgvs_inversion(ref_allele, alt_allele_raw) {
        if let Some(ref pos_e) = pos_end {
            format!(".{pos_start}_{pos_e}inv")
        } else {
            format!(".{pos_start}inv")
        }
    } else if ref_allele.len() > 1 && alt_allele_raw.len() > 1 && !ref_is_dash && !alt_is_dash {
        if ref_allele.len() == 1 {
            format!(".{pos_start}delins{display_alt}")
        } else if let Some(ref pos_e) = pos_end {
            format!(".{pos_start}_{pos_e}delins{display_alt}")
        } else {
            format!(".{pos_start}delins{display_alt}")
        }
    } else if ref_allele.len() == 1 && alt_allele_raw.len() == 1 && !ref_is_dash && !alt_is_dash {
        format!(".{pos_start}{display_ref}>{display_alt}")
    } else if ref_is_dash && !alt_is_dash {
        let prev_pos =
            format_cdna_position(cdna_s.saturating_sub(1).max(1), mapper, has_coding_model);
        format!(".{prev_pos}_{pos_start}ins{display_alt}")
    } else if !ref_is_dash && alt_is_dash {
        if let Some(ref pos_e) = pos_end {
            format!(".{pos_start}_{pos_e}del")
        } else {
            format!(".{pos_start}del")
        }
    } else {
        if let Some(ref pos_e) = pos_end {
            format!(".{pos_start}_{pos_e}delins{display_alt}")
        } else {
            format!(".{pos_start}delins{display_alt}")
        }
    };

    Some(format!("{transcript_prefix}{variant_desc}"))
}

/// Format a cDNA position as an HGVS coding position.
///
/// Returns positions like:
/// - `"42"`: coding region (c.42)
/// - `"-14"`: 5' UTR (c.-14)
/// - `"*42"`: 3' UTR (c.*42)
fn format_cdna_position(
    cdna_pos: u64,
    mapper: &vep_core::transcript::TranscriptMapper,
    has_coding_model: bool,
) -> String {
    if !has_coding_model {
        return cdna_pos.to_string();
    }

    let coding_start = mapper.cdna_coding_start;
    let coding_end = mapper.cdna_coding_end;
    let start_phase_offset = u64::try_from(mapper.start_phase).unwrap_or(0);

    if cdna_pos < coding_start {
        let offset = coding_start - cdna_pos;
        format!("-{offset}")
    } else if cdna_pos > coding_end {
        let offset = cdna_pos - coding_end;
        format!("*{offset}")
    } else {
        let cds_pos = cdna_pos - coding_start + 1 + start_phase_offset;
        cds_pos.to_string()
    }
}

/// Shared context for intronic HGVS coordinate calculations.
struct IntronicHgvsContext<'a> {
    transcript: &'a Transcript,
    pairs: &'a [vep_core::transcript::MapperPair],
    mapper: &'a vep_core::transcript::TranscriptMapper,
    has_coding_model: bool,
}

/// Generate HGVSc for intronic variants using (possibly shifted) coordinates.
///
/// Uses offset notation relative to nearest exon boundary:
/// - `c.803+2T>C` (donor/5' splice site)
/// - `c.804-3T>C` (acceptor/3' splice site)
fn generate_hgvsc_intronic_shifted(
    variant: &InputVariant,
    var_start: u64,
    var_end: u64,
    ctx: &IntronicHgvsContext<'_>,
    transcript_prefix: &str,
    display_ref: &str,
    display_alt: &str,
) -> Option<String> {
    let vefc = ctx.transcript.vefc.as_ref()?;
    let pairs = ctx.pairs;

    let introns = if !vefc.introns.is_empty() {
        &vefc.introns
    } else {
        &ctx.transcript.introns
    };

    for (intron_idx, intron) in introns.iter().enumerate() {
        let in_intron = var_start >= intron.start && var_start <= intron.end;
        if !in_intron {
            continue;
        }

        let (dist_to_prev_exon, dist_to_next_exon) = match ctx.transcript.strand {
            Strand::Forward => {
                let d_donor = var_start - intron.start + 1;
                let d_acceptor = intron.end - var_start + 1;
                (d_donor, d_acceptor)
            }
            Strand::Reverse => {
                let d_donor = intron.end - var_start + 1;
                let d_acceptor = var_start - intron.start + 1;
                (d_donor, d_acceptor)
            }
        };

        let use_donor = dist_to_prev_exon <= dist_to_next_exon;

        let (base_cdna, offset, sign) = if use_donor {
            let preceding_exon_idx = intron_idx;
            if preceding_exon_idx < pairs.len() {
                let pair = &pairs[preceding_exon_idx];
                let exon_cdna_end = pair.from_end;
                (exon_cdna_end, dist_to_prev_exon as i64, "+")
            } else {
                return None;
            }
        } else {
            let following_exon_idx = intron_idx + 1;
            if following_exon_idx < pairs.len() {
                let pair = &pairs[following_exon_idx];
                let exon_cdna_start = pair.from_start;
                (exon_cdna_start, -(dist_to_next_exon as i64), "-")
            } else {
                return None;
            }
        };

        // Perl VEP does not apply start_phase_offset to the intronic reference
        // position. Counteract the offset that format_cdna_position will add by
        // pre-subtracting it from the cDNA coordinate.
        let spo = u64::try_from(ctx.mapper.start_phase).unwrap_or(0);
        let adjusted_cdna = if ctx.has_coding_model && spo > 0 && base_cdna >= spo {
            base_cdna - spo
        } else {
            base_cdna
        };
        let base_pos = format_cdna_position(adjusted_cdna, ctx.mapper, ctx.has_coding_model);
        let offset_abs = offset.unsigned_abs();

        let is_snv = variant.variant_class == VariantClass::Snv;

        let variant_desc = if is_snv {
            format!(".{base_pos}{sign}{offset_abs}{display_ref}>{display_alt}")
        } else if variant.variant_class == VariantClass::Deletion {
            if var_start == var_end {
                format!(".{base_pos}{sign}{offset_abs}del")
            } else {
                let end_offset = compute_intronic_end_offset(var_end, intron, intron_idx, ctx);
                if let Some((end_base_pos, end_sign, end_offset_abs)) = end_offset {
                    // Ensure HGVS ascending genomic order:
                    // For donor (+): smaller offset first (closer to exon)
                    // For acceptor (-): larger offset first (further from exon)
                    let needs_swap = base_pos == end_base_pos
                        && sign == end_sign
                        && ((sign == "+" && offset_abs > end_offset_abs)
                            || (sign == "-" && offset_abs < end_offset_abs));
                    let (p1, s1, o1, p2, s2, o2) = if needs_swap {
                        (
                            &end_base_pos,
                            end_sign,
                            end_offset_abs,
                            &base_pos,
                            sign,
                            offset_abs,
                        )
                    } else {
                        (
                            &base_pos,
                            sign,
                            offset_abs,
                            &end_base_pos,
                            end_sign,
                            end_offset_abs,
                        )
                    };
                    format!(".{p1}{s1}{o1}_{p2}{s2}{o2}del")
                } else {
                    format!(".{base_pos}{sign}{offset_abs}del")
                }
            }
        } else if variant.variant_class == VariantClass::Insertion {
            let prev_offset = if offset > 0 { offset - 1 } else { offset };
            let next_offset = offset;
            if prev_offset == 0 {
                format!(".{base_pos}_{base_pos}{sign}{next_offset}ins{display_alt}")
            } else {
                let prev_abs = prev_offset.unsigned_abs();
                format!(".{base_pos}{sign}{prev_abs}_{base_pos}{sign}{offset_abs}ins{display_alt}")
            }
        } else {
            if var_start != var_end {
                let end_offset = compute_intronic_end_offset(var_end, intron, intron_idx, ctx);
                if let Some((end_base_pos, end_sign, end_offset_abs)) = end_offset {
                    format!(".{base_pos}{sign}{offset_abs}_{end_base_pos}{end_sign}{end_offset_abs}delins{display_alt}")
                } else {
                    format!(".{base_pos}{sign}{offset_abs}delins{display_alt}")
                }
            } else {
                format!(".{base_pos}{sign}{offset_abs}delins{display_alt}")
            }
        };

        return Some(format!("{transcript_prefix}{variant_desc}"));
    }

    None
}

/// Compute intronic end offset for multi-base intronic variants.
fn compute_intronic_end_offset(
    end_pos: u64,
    intron: &vep_core::transcript::Intron,
    intron_idx: usize,
    ctx: &IntronicHgvsContext<'_>,
) -> Option<(String, &'static str, u64)> {
    if end_pos >= intron.start && end_pos <= intron.end {
        let (dist_to_prev, dist_to_next) = match ctx.transcript.strand {
            Strand::Forward => {
                let d_donor = end_pos - intron.start + 1;
                let d_acceptor = intron.end - end_pos + 1;
                (d_donor, d_acceptor)
            }
            Strand::Reverse => {
                let d_donor = intron.end - end_pos + 1;
                let d_acceptor = end_pos - intron.start + 1;
                (d_donor, d_acceptor)
            }
        };

        let use_donor = dist_to_prev <= dist_to_next;
        let spo = u64::try_from(ctx.mapper.start_phase).unwrap_or(0);
        if use_donor {
            let preceding_exon_idx = intron_idx;
            if preceding_exon_idx < ctx.pairs.len() {
                let pair = &ctx.pairs[preceding_exon_idx];
                let cdna = if ctx.has_coding_model && spo > 0 && pair.from_end >= spo {
                    pair.from_end - spo
                } else {
                    pair.from_end
                };
                let base_pos = format_cdna_position(cdna, ctx.mapper, ctx.has_coding_model);
                return Some((base_pos, "+", dist_to_prev));
            }
        } else {
            let following_exon_idx = intron_idx + 1;
            if following_exon_idx < ctx.pairs.len() {
                let pair = &ctx.pairs[following_exon_idx];
                let cdna = if ctx.has_coding_model && spo > 0 && pair.from_start >= spo {
                    pair.from_start - spo
                } else {
                    pair.from_start
                };
                let base_pos = format_cdna_position(cdna, ctx.mapper, ctx.has_coding_model);
                return Some((base_pos, "-", dist_to_next));
            }
        }
    }

    None
}

/// Convert genomic position to cDNA position using mapper pairs.
fn genomic_to_cdna_pos(
    genomic_pos: u64,
    pairs: &[vep_core::transcript::MapperPair],
) -> Option<u64> {
    for pair in pairs {
        if genomic_pos >= pair.to_start && genomic_pos <= pair.to_end {
            let cdna_pos = if pair.ori == 1 {
                pair.from_start + (genomic_pos - pair.to_start)
            } else {
                pair.from_start + (pair.to_end - genomic_pos)
            };
            return Some(cdna_pos);
        }
    }
    None
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

fn is_hgvs_inversion(ref_allele: &[u8], alt_allele: &[u8]) -> bool {
    ref_allele.len() > 1
        && ref_allele.len() == alt_allele.len()
        && !ref_allele.eq_ignore_ascii_case(alt_allele)
        && reverse_complement_string(ref_allele)
            .as_bytes()
            .eq_ignore_ascii_case(alt_allele)
}

/// Check if an insertion is actually a duplication of the preceding sequence.
///
/// When `reference_fasta` is available, uses the genomic reference sequence
/// for duplication detection (works in UTR, intronic, intergenic contexts).
/// Falls back to CDS-only check when FASTA is unavailable.
fn check_duplication_cdna(
    variant: &InputVariant,
    var_start: u64,
    var_end: u64,
    transcript: &Transcript,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<String> {
    let alt_allele = variant.alt_allele();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();
    if alt_is_dash {
        return None;
    }
    let alt_len = alt_allele.len();
    if alt_len == 0 {
        return None;
    }

    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let pairs = &mapper.exon_coord_mapper.pairs;

    let has_coding_model = transcript.translation.is_some()
        && mapper.cdna_coding_start > 0
        && mapper.cdna_coding_end >= mapper.cdna_coding_start;

    // After 3' shifting the alt allele is rotated, so the reference bases
    // adjacent to the insertion are compared against it by cyclic rotation.
    // Forward strand shifts right, so the duplicated region precedes the
    // insertion point; reverse shifts left, so it follows.
    if let Some(fasta) = reference_fasta {
        let al = alt_len as u64;
        // For VEP-convention insertions, var_end < var_start.
        if transcript.strand == Strand::Forward {
            let ins_pos = var_end;
            if ins_pos >= al {
                let dup_start = ins_pos - al + 1;
                let dup_end = ins_pos;
                if let Some(ref_seg) = fasta.sequence(&variant.chr, dup_start, dup_end) {
                    if ref_seg.len() == alt_len && is_rotation(&ref_seg, alt_allele) {
                        return build_dup_notation(
                            dup_start,
                            dup_end,
                            transcript,
                            pairs,
                            mapper,
                            has_coding_model,
                        );
                    }
                }
            }
        } else {
            // Reverse strand: duplicated region is after the insertion point (higher coords)
            let dup_start = var_start;
            let dup_end = var_start + al - 1;
            if let Some(ref_seg) = fasta.sequence(&variant.chr, dup_start, dup_end) {
                if ref_seg.len() == alt_len && is_rotation(&ref_seg, alt_allele) {
                    return build_dup_notation(
                        dup_start,
                        dup_end,
                        transcript,
                        pairs,
                        mapper,
                        has_coding_model,
                    );
                }
            }
        }
        return None;
    }

    // Without a FASTA, only the CDS (translateable_seq) can be checked.
    let translateable_seq = vefc.translateable_seq.as_ref()?;
    let cds = translateable_seq.as_bytes();

    if !has_coding_model {
        return None;
    }

    let cdna_pos = genomic_to_cdna_pos(var_start, pairs)?;
    let coding_start = mapper.cdna_coding_start;
    let coding_end = mapper.cdna_coding_end;

    if cdna_pos < coding_start || cdna_pos > coding_end {
        return None;
    }

    let start_phase_offset = u64::try_from(mapper.start_phase).unwrap_or(0);
    let cds_pos = cdna_pos - coding_start + 1 + start_phase_offset;
    let idx = (cds_pos - 1) as usize;

    if idx < alt_len {
        return None;
    }

    let preceding_start = idx - alt_len;
    if preceding_start + alt_len > cds.len() {
        return None;
    }
    let preceding = &cds[preceding_start..preceding_start + alt_len];

    let alt_mrna = if transcript.strand == Strand::Reverse {
        crate::coding::reverse_complement_pub(alt_allele)
    } else {
        alt_allele.to_vec()
    };

    if preceding.eq_ignore_ascii_case(&alt_mrna) {
        let dup_start = format_cdna_position(
            coding_start + (preceding_start as u64) - start_phase_offset,
            mapper,
            has_coding_model,
        );
        if alt_len == 1 {
            Some(format!(".{dup_start}dup"))
        } else {
            let dup_end = format_cdna_position(
                coding_start + (preceding_start + alt_len - 1) as u64 - start_phase_offset,
                mapper,
                has_coding_model,
            );
            Some(format!(".{dup_start}_{dup_end}dup"))
        }
    } else {
        None
    }
}

const HGVS_POSITION_SCALE: i64 = 1_000_000;

struct HgvsPosition {
    display: String,
    sort_key: i64,
}

/// Build duplication notation from genomic coordinates of the duplicated region.
fn build_dup_notation(
    genomic_dup_start: u64,
    genomic_dup_end: u64,
    transcript: &Transcript,
    pairs: &[vep_core::transcript::MapperPair],
    mapper: &vep_core::transcript::TranscriptMapper,
    has_coding_model: bool,
) -> Option<String> {
    let start_pos = genomic_to_hgvs_position(
        genomic_dup_start,
        transcript,
        pairs,
        mapper,
        has_coding_model,
    )?;
    let end_pos =
        genomic_to_hgvs_position(genomic_dup_end, transcript, pairs, mapper, has_coding_model)?;

    let (first, second) = if start_pos.sort_key <= end_pos.sort_key {
        (start_pos, end_pos)
    } else {
        (end_pos, start_pos)
    };

    if first.sort_key == second.sort_key {
        Some(format!(".{}dup", first.display))
    } else {
        Some(format!(".{}_{}dup", first.display, second.display))
    }
}

fn genomic_to_hgvs_position(
    genomic_pos: u64,
    transcript: &Transcript,
    pairs: &[vep_core::transcript::MapperPair],
    mapper: &vep_core::transcript::TranscriptMapper,
    has_coding_model: bool,
) -> Option<HgvsPosition> {
    if let Some(cdna_pos) = genomic_to_cdna_pos(genomic_pos, pairs) {
        return Some(HgvsPosition {
            display: format_cdna_position(cdna_pos, mapper, has_coding_model),
            sort_key: i64::try_from(cdna_pos).ok()? * HGVS_POSITION_SCALE,
        });
    }

    let vefc = transcript.vefc.as_ref()?;
    let introns = if !vefc.introns.is_empty() {
        &vefc.introns
    } else {
        &transcript.introns
    };

    let spo = u64::try_from(mapper.start_phase).unwrap_or(0);

    for (intron_idx, intron) in introns.iter().enumerate() {
        if genomic_pos < intron.start || genomic_pos > intron.end {
            continue;
        }

        let (dist_to_prev, dist_to_next) = match transcript.strand {
            Strand::Forward => (genomic_pos - intron.start + 1, intron.end - genomic_pos + 1),
            Strand::Reverse => (intron.end - genomic_pos + 1, genomic_pos - intron.start + 1),
        };

        let use_donor = dist_to_prev <= dist_to_next;
        if use_donor {
            let pair = pairs.get(intron_idx)?;
            let adjusted_cdna = if has_coding_model && spo > 0 && pair.from_end >= spo {
                pair.from_end - spo
            } else {
                pair.from_end
            };
            let offset = i64::try_from(dist_to_prev).ok()?;
            return Some(HgvsPosition {
                display: format!(
                    "{}+{}",
                    format_cdna_position(adjusted_cdna, mapper, has_coding_model),
                    dist_to_prev
                ),
                sort_key: i64::try_from(pair.from_end).ok()? * HGVS_POSITION_SCALE + offset,
            });
        }

        let pair = pairs.get(intron_idx + 1)?;
        let adjusted_cdna = if has_coding_model && spo > 0 && pair.from_start >= spo {
            pair.from_start - spo
        } else {
            pair.from_start
        };
        let offset = i64::try_from(dist_to_next).ok()?;
        return Some(HgvsPosition {
            display: format!(
                "{}-{}",
                format_cdna_position(adjusted_cdna, mapper, has_coding_model),
                dist_to_next
            ),
            sort_key: i64::try_from(pair.from_start).ok()? * HGVS_POSITION_SCALE - offset,
        });
    }

    None
}

/// Check if `a` is a cyclic rotation of `b` (case-insensitive).
fn is_rotation(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() || a.is_empty() {
        return false;
    }
    // (b ++ b) contains a iff a is a rotation of b
    let mut doubled: Vec<u8> = b.iter().map(|c| c.to_ascii_uppercase()).collect();
    doubled.extend(b.iter().map(|c| c.to_ascii_uppercase()));
    let needle: Vec<u8> = a.iter().map(|c| c.to_ascii_uppercase()).collect();
    doubled
        .windows(needle.len())
        .any(|w| w == needle.as_slice())
}

fn reverse_complement_string(seq: &[u8]) -> String {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'N',
        })
        .collect()
}

/// Generate HGVSp notation for a variant-transcript pair: Perl VEP's
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

    let shifted = reference_fasta.and_then(|fasta| {
        crate::consequences::shift_indel_3prime_coords(variant, transcript, fasta, 2000)
    });
    let notation = match shifted {
        Some((start, end)) if (start, end) != (variant.start, variant.end) => {
            let ref_allele = normalized_hgvs_ref_allele(variant, shifted, reference_fasta)
                .unwrap_or_else(|| b"-".to_vec());
            let alt_allele = normalized_hgvs_alt_allele(variant, transcript, shifted)
                .unwrap_or_else(|| b"-".to_vec());
            let shifted_variant =
                InputVariant::new(variant.chr.clone(), start, end, ref_allele, alt_allele);
            crate::coding::perl_hgvs_protein(
                &shifted_variant,
                transcript,
                start,
                end,
                reference_fasta,
            )?
        }
        _ => crate::coding::perl_hgvs_protein(
            variant,
            transcript,
            variant.start,
            variant.end,
            reference_fasta,
        )?,
    };
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

    #[test]
    fn test_format_cdna_position_coding() {
        let tx = make_test_transcript();
        let vefc = tx.vefc.as_ref().unwrap();
        let mapper = vefc.mapper.as_ref().unwrap();

        // CDS pos 1 = cDNA 51
        assert_eq!(format_cdna_position(51, mapper, true), "1");
        // CDS pos 10 = cDNA 60
        assert_eq!(format_cdna_position(60, mapper, true), "10");
    }

    #[test]
    fn test_format_cdna_position_five_prime_utr() {
        let tx = make_test_transcript();
        let vefc = tx.vefc.as_ref().unwrap();
        let mapper = vefc.mapper.as_ref().unwrap();

        // cDNA 50 is 1bp before coding start (51) -> c.-1
        assert_eq!(format_cdna_position(50, mapper, true), "-1");
        // cDNA 1 is 50bp before coding start -> c.-50
        assert_eq!(format_cdna_position(1, mapper, true), "-50");
    }

    #[test]
    fn test_format_cdna_position_three_prime_utr() {
        let tx = make_test_transcript();
        let vefc = tx.vefc.as_ref().unwrap();
        let mapper = vefc.mapper.as_ref().unwrap();

        // cDNA 901 is 1bp after coding end (900) -> c.*1
        assert_eq!(format_cdna_position(901, mapper, true), "*1");
        // cDNA 910 is 10bp after coding end -> c.*10
        assert_eq!(format_cdna_position(910, mapper, true), "*10");
    }

    #[test]
    fn test_format_cdna_position_non_coding() {
        let tx = make_test_transcript();
        let vefc = tx.vefc.as_ref().unwrap();
        let mapper = vefc.mapper.as_ref().unwrap();

        // Non-coding: just the raw cDNA position
        assert_eq!(format_cdna_position(42, mapper, false), "42");
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
        let hgvsc = generate_hgvsc(&variant, &tx, None);
        assert!(hgvsc.is_some(), "HGVSc should be generated for coding SNV");
        let notation = hgvsc.unwrap();
        assert!(
            notation.starts_with("ENST00000000001.1:c."),
            "Should have transcript prefix, got: {notation}"
        );
        assert!(
            notation.contains(">"),
            "SNV should use substitution notation, got: {notation}"
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
        let hgvsc = generate_hgvsc(&variant, &tx, None);
        assert!(hgvsc.is_some());
        let notation = hgvsc.unwrap();
        assert!(
            notation.contains("c.-"),
            "5' UTR should use negative notation, got: {notation}"
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
        let hgvsc = generate_hgvsc(&variant, &tx, None);
        assert!(hgvsc.is_some());
        let notation = hgvsc.unwrap();
        assert!(
            notation.contains("c.*"),
            "3' UTR should use * notation, got: {notation}"
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
        let hgvsc = generate_hgvsc(&variant, &tx, None);
        assert!(hgvsc.is_some());
        let notation = hgvsc.unwrap();
        assert!(
            notation.contains("+"),
            "Intronic donor side should use + notation, got: {notation}"
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
        let hgvsc = generate_hgvsc(&variant, &tx, None);
        assert!(hgvsc.is_some());
        let notation = hgvsc.unwrap();
        assert!(
            notation.contains("del"),
            "Deletion should use del notation, got: {notation}"
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
        let hgvsc = generate_hgvsc(&variant, &tx, None);
        assert!(hgvsc.is_some());
        let notation = hgvsc.unwrap();
        assert!(
            notation.contains(":n."),
            "Non-coding transcript should use n. prefix, got: {notation}"
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

    #[test]
    fn test_reverse_complement_string() {
        assert_eq!(reverse_complement_string(b"ATCG"), "CGAT");
        assert_eq!(reverse_complement_string(b"A"), "T");
        assert_eq!(reverse_complement_string(b"AAAA"), "TTTT");
    }
}
