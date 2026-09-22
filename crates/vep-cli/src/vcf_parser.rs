// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VCF line parsing and allele processing.
//!
//! Handles parsing of VCF data lines into `InputVariant`s, including
//! multi-allelic splitting, structural variant classification, BND
//! notation parsing, and allele trimming/normalization.
//!
//! Perl citations name modules of ensembl-vep release/115 (`Bio/EnsEMBL/VEP/...`);
//! `Sequence.pm` and `StructuralVariationOverlap.pm` are in ensembl-variation
//! release/115 and `BaseVCF4.pm` is in ensembl-io release/115.

use anyhow::{bail, Context};

use vep_core::allele::trim_alleles;
use vep_core::variant::InputVariant;

/// Parse a VCF data line into one or more `InputVariant`s.
///
/// Multi-allelic sites are split into separate records, one per ALT allele.
/// This is a minimal parser; the noodles-based parser is in the vep-io crate.
///
/// `capture_raw_input` controls whether each produced variant retains a copy of
/// the source line in `InputVariant::raw_input`. Only the VCF output formatter
/// reads that field (it re-emits the original record with a CSQ INFO field
/// appended); every other output format ignores it, and the formatter has a
/// complete fallback that rebuilds a minimal record from the variant's own
/// coordinates and alleles. Capturing unconditionally would allocate one
/// line-length `String` per variant only to drop it. Callers pass `true` only
/// when the output format needs it; see `RAW_INPUT_OUTPUT_FORMATS` in `runner.rs`.
pub fn parse_vcf_line(line: &str, capture_raw_input: bool) -> anyhow::Result<Vec<InputVariant>> {
    let fields = split_fixed_fields(line)?;

    let chr = fields[0].to_string();
    let pos: u64 = fields[1]
        .parse()
        .with_context(|| format!("Invalid POS: '{}'", fields[1]))?;
    let id = fields[2];
    let ref_allele = fields[3];
    let alt_field = fields[4];

    if alt_field == "." || alt_field == "*" {
        bail!("No alternate allele");
    }

    let alts: Vec<&str> = alt_field.split(',').collect();
    let valid_alts: Vec<(usize, &str)> = alts
        .iter()
        .enumerate()
        .filter_map(|(i, alt)| {
            // Perl VEP skips `*` alleles (neither the indel nor the SV path fires)
            // and `.`. `<*>` is not filtered: Perl treats it as a non-SV allele at
            // POS (its SV regex excludes `<*>`). `<NON_REF>` and `<CPX>` are not
            // filtered either: absent from `%SO_TERMS`, they are still annotated in
            // tab output because `vep_skip` is checked only by File-based
            // annotation sources and VCF/JSON output.
            if *alt == "*" || *alt == "." {
                None
            } else {
                Some((i, *alt))
            }
        })
        .collect();
    if valid_alts.is_empty() {
        bail!("No alternate allele");
    }

    let info_field = fields[7];

    let mut symbolic_alts: Vec<(usize, &str)> = Vec::with_capacity(valid_alts.len());
    let mut explicit_alts: Vec<(usize, &str)> = Vec::with_capacity(valid_alts.len());
    for (allele_idx, alt) in &valid_alts {
        if is_symbolic_alt(alt) {
            symbolic_alts.push((*allele_idx, alt));
        } else {
            explicit_alts.push((*allele_idx, alt));
        }
    }

    let mut variants = Vec::with_capacity(valid_alts.len());
    let is_mixed = !explicit_alts.is_empty() && !symbolic_alts.is_empty();

    // In a mixed record (`C  G,<INV>`) Perl VEP processes only the symbolic
    // allele: VCF.pm builds separate VariationFeature and
    // StructuralVariationFeature objects from one line, and the explicit path
    // fires only when no symbolic allele is present.
    if !explicit_alts.is_empty() && symbolic_alts.is_empty() {
        let trim_alts: Vec<&str> = explicit_alts.iter().map(|(_, alt)| *alt).collect();
        let (trimmed_ref, trimmed_alts, prefix_trimmed) =
            trim_multi_alleles(ref_allele, &trim_alts);
        let (vep_start, vep_end) =
            vep_coords_from_trimmed(pos, prefix_trimmed, trimmed_ref.as_bytes());

        // The allele string VEP names an ID-less record by: the raw REF and every
        // ALT joined by `/`. VEP keeps the raw string as `nontrimmed_allele_string`
        // and promotes it to `original_allele_string` for any record whose alleles
        // differ in length; a record without an indel allele is never trimmed, so
        // its raw and working strings coincide (Parser/VCF.pm, Parser.pm
        // `post_process_vfs` and `minimise_alleles`, OutputFactory.pm).
        let uploaded_alleles = {
            let mut s = ref_allele.to_string();
            for a in &trim_alts {
                s.push('/');
                s.push_str(a);
            }
            s
        };

        for (trim_idx, (allele_idx, _alt)) in explicit_alts.iter().enumerate() {
            let vep_alt = trimmed_alts
                .get(trim_idx)
                .with_context(|| format!("missing normalized ALT for trim index {trim_idx}"))?;

            let mut variant = InputVariant::new(
                chr.clone(),
                vep_start,
                vep_end,
                trimmed_ref.as_bytes().to_vec(),
                vep_alt.as_bytes().to_vec(),
            );

            // Secondary minimisation of a bi-allelic indel. Perl VEP
            // (Parser.pm `post_process_vfs`) minimises a bi-allelic record whose
            // alleles differ in length whatever the flags: `trim_sequences`
            // (Sequence.pm) trims from the left while both alleles are non-empty,
            // then from the right, and maps an emptied allele to `-`, so `T/TT`
            // at pos P becomes `-/T` at P+1. `minimise_alleles` returns a
            // multi-allelic string (`.+/.+/.+`) unchanged: without `--minimal` a
            // multi-allelic record keeps its anchor-chopped alleles, and every
            // position, codon and distance is computed on that representation.
            if explicit_alts.len() == 1
                && variant.ref_allele != b"-"
                && variant.alt_alleles[0] != b"-"
                && variant.ref_allele.len() != variant.alt_alleles[0].len()
            {
                let (min_ref, min_alt, extra_prefix) =
                    trim_alleles(&variant.ref_allele, &variant.alt_alleles[0]);
                if extra_prefix > 0
                    || min_ref != variant.ref_allele
                    || min_alt != variant.alt_alleles[0]
                {
                    let total_prefix = prefix_trimmed + extra_prefix;
                    let (new_start, new_end) = vep_coords_from_trimmed(pos, total_prefix, &min_ref);
                    variant.original_allele_string = Some(variant.allele_string.clone());
                    variant.allele_string = format!(
                        "{}/{}",
                        String::from_utf8_lossy(&min_ref),
                        String::from_utf8_lossy(&min_alt),
                    );
                    variant.ref_allele = min_ref;
                    variant.alt_alleles = vec![min_alt];
                    variant.start = new_start;
                    variant.end = new_end;
                    variant.variant_class = vep_core::variant::classify_variant(
                        &variant.ref_allele,
                        &variant.alt_alleles[0],
                    );
                    variant.minimised = true;
                }
            }

            variant.allele_index = *allele_idx;
            variant.uploaded_allele_string = Some(uploaded_alleles.clone());
            if explicit_alts.len() > 1 {
                let mut chopped = trimmed_ref.clone();
                for a in &trimmed_alts {
                    chopped.push('/');
                    chopped.push_str(a);
                }
                variant.record_allele_string_multi = Some(chopped);
            }
            if capture_raw_input {
                variant.raw_input = Some(line.to_string());
            }
            if id != "." && !id.is_empty() {
                variant.id = Some(id.to_string());
            }
            variants.push(variant);
        }
    }

    // Perl VEP (Parser.pm:get_SO_term) joins multi-allelic ALTs with "/": for
    // "<CN0>,<CN2>" the joined "<CN0>/<CN2>" matches no single-allele regex and
    // falls through to generic "copy_number_variation" as one SVF, so multiple
    // CN-type alleles coalesce into a single CopyNumberVariation variant.
    let coalesce_cn = symbolic_alts.len() > 1
        && explicit_alts.is_empty()
        && symbolic_alts.iter().all(|(_, alt)| {
            let upper = alt.to_uppercase();
            upper.starts_with("<CN") && !upper.starts_with("<CNV") && upper.ends_with('>')
        });

    let symbolic_alts_to_process: Vec<(usize, &str)> = if coalesce_cn {
        // ALT becomes <CNV> so `parse_copy_number()` in cnv.rs returns Generic
        // rather than the first allele's Deletion/Duplication.
        vec![(symbolic_alts[0].0, "<CNV>")]
    } else {
        symbolic_alts.clone()
    };

    for (allele_idx, alt) in &symbolic_alts_to_process {
        let svtype = get_info_value(info_field, "SVTYPE");
        let sv_class = if coalesce_cn {
            vep_core::variant::VariantClass::CopyNumberVariation
        } else if is_mixed {
            if let Some(svtype) = svtype {
                let fallback_alt = format!("<{svtype}>");
                classify_symbolic_alt(&fallback_alt, Some(svtype))
            } else {
                classify_symbolic_alt(alt, svtype)
            }
        } else {
            classify_symbolic_alt(alt, svtype)
        };
        // Perl VEP (ensembl-io BaseVCF4.pm get_end()) universally prioritises
        // SVLEN over INFO/END for all SV types:
        //   if(defined($info->{SVLEN})) { $end = $self->get_start + abs($svlen) - 1; }
        //   elsif(defined($info->{END})) { $end = $info->{END}; }
        // Since get_start() = POS+1 for SVs, this is: end = POS + abs(SVLEN).
        let mut sv_end_val = get_info_value(info_field, "SVLEN")
            .and_then(|v| v.parse::<i64>().ok())
            .map(|len| pos + len.unsigned_abs())
            .or_else(|| get_info_value(info_field, "END").and_then(|v| v.parse::<u64>().ok()))
            .unwrap_or(pos);

        // Perl VEP uses POS+1 for all SV types (the VCF POS is the anchor base;
        // the actual structural event starts at POS+1).
        let sv_start = if sv_class.is_span_type() {
            pos + 1
        } else {
            pos
        };

        // Perl VEP uses SVLEN (not END) for tandem repeat end coordinates.
        // VCF.pm:502-512: $end = $start + max(@svlen) - 1 when $so_term =~ /tandem/
        if sv_class == vep_core::variant::VariantClass::TandemRepeat {
            if let Some(svlen_str) = get_info_value(info_field, "SVLEN") {
                if let Ok(svlen) = svlen_str.parse::<i64>() {
                    sv_end_val = sv_start + svlen.unsigned_abs() - 1;
                }
            }
        }

        // For a symbolic <BND> Perl VEP spans POS+1 to POS+1+SVLEN from SVLEN,
        // not INFO/END (typically POS+1 on BND records), then synthesizes paired
        // and single-breakend forms.
        let is_symbolic_bnd =
            sv_class == vep_core::variant::VariantClass::Translocation && alt.starts_with('<');
        if is_symbolic_bnd {
            if let Some(svlen_str) = get_info_value(info_field, "SVLEN") {
                if let Ok(svlen) = svlen_str.parse::<i64>() {
                    sv_end_val = sv_start + svlen.unsigned_abs() - 1;
                }
            }
        }

        // A point-type SV (BND, ME, INS) without END/SVLEN defaults to POS; the
        // event starts at POS+1, so the clamp keeps it a single position.
        let sv_end_val = sv_end_val.max(sv_start);
        let mut variant = InputVariant::new(
            chr.clone(),
            sv_start,
            sv_end_val,
            ref_allele.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        );
        variant.variant_class = sv_class;
        variant.is_structural = true;
        // VEP's allele string for a structural variant record is every ALT of
        // the record joined by `/`, explicit alleles of a mixed record included,
        // never the REF base.
        let joined_alts = valid_alts
            .iter()
            .map(|(_, a)| *a)
            .collect::<Vec<_>>()
            .join("/");
        variant.vep_skip =
            !vep_supports_sv_type(&joined_alts, get_info_value(info_field, "SVTYPE"));
        variant.uploaded_allele_string = Some(joined_alts);
        variant.sv_end = Some(sv_end_val);
        variant.sv_type = get_info_value(info_field, "SVTYPE").map(|s| s.to_string());
        variant.sv_len = get_info_value(info_field, "SVLEN").and_then(|v| v.parse().ok());
        if sv_class == vep_core::variant::VariantClass::TandemRepeat {
            variant.tr_alt_bases = tandem_repeat_alt_bases(info_field);
        }
        variant.ci_pos = parse_ci_field(get_info_value(info_field, "CIPOS"));
        variant.ci_end = parse_ci_field(get_info_value(info_field, "CIEND"));
        variant.mate_id = get_info_value(info_field, "MATEID")
            .or_else(|| get_info_value(info_field, "PARID"))
            .map(|s| s.to_string());
        if sv_class == vep_core::variant::VariantClass::Translocation {
            if let Some((mate_chr, mate_pos)) = parse_bnd_mate(alt) {
                variant.mate_chr = Some(mate_chr);
                variant.mate_pos = Some(mate_pos);
            }
        }
        variant.allele_index = *allele_idx;
        if capture_raw_input {
            variant.raw_input = Some(line.to_string());
        }
        if id != "." && !id.is_empty() {
            variant.id = Some(id.to_string());
        }
        let is_paired_bnd = is_paired_breakend_alt(alt);
        let is_native_single_breakend = is_native_single_breakend_alt(alt);
        variant.is_single_breakend = is_native_single_breakend;
        // Perl VEP creates two allele entries per BND record, the paired notation
        // ("T[21:33033339[") and the single-breakend form ("T." or ".T"), each
        // annotated independently against all overlapping transcripts.
        let is_bnd = sv_class == vep_core::variant::VariantClass::Translocation;
        // A symbolic <BND> with CHR2/END2 enters dual emission even when SVLEN
        // gives end == start: Perl VEP creates both forms regardless of span.
        let has_chr2_end2 = get_info_value(info_field, "CHR2").is_some()
            && (get_info_value(info_field, "END2").is_some()
                || get_info_value(info_field, "POS2").is_some());
        if is_bnd && is_symbolic_bnd && (variant.end > variant.start || has_chr2_end2) {
            let ref_str = std::str::from_utf8(ref_allele.as_bytes()).unwrap_or("N");

            // Perl VCF.pm:483-498 stashes the raw `INFO/END` as `$incorrect_end`,
            // deletes it so it cannot be read as a span end, and with CHR2 present
            // resolves `my $breakend_pos = $info->{END2} || $incorrect_end;`. The
            // END2-absent fallback is therefore the raw `INFO/END`, which Manta and
            // gnomAD use for the mate position although VCF 4.4 does not sanction
            // it. POS2 is never consulted, so it is not in the chain below.
            // The allele string keeps CHR2 as the input spells it (`chr21`): VEP
            // builds `N[chr21:pos[` from the INFO value verbatim, while the mate
            // lookup below uses the cache's chromosome name.
            let chr2_raw = get_info_value(info_field, "CHR2").map(str::to_string);
            let chr2 = chr2_raw
                .as_deref()
                .map(vep_core::coordinate::normalize_chromosome);
            let end2 = get_info_value(info_field, "END2")
                .and_then(|v| v.parse::<u64>().ok())
                .or_else(|| get_info_value(info_field, "END").and_then(|v| v.parse::<u64>().ok()));

            let (mate_chr_val, mate_chr_label, mate_pos) =
                if let (Some(ref c2), Some(ref raw), Some(e2)) = (&chr2, &chr2_raw, end2) {
                    (c2.clone(), raw.clone(), e2)
                } else {
                    // No CHR2: Perl leaves the allele symbolic, so the SVLEN-derived
                    // span end is the only available mate anchor.
                    (chr.clone(), chr.clone(), variant.end + 1)
                };

            let paired_alt = format!("{}[{}:{}[", ref_str, mate_chr_label, mate_pos);
            variant.alt_alleles = vec![paired_alt.as_bytes().to_vec()];
            variant.allele_string = format!("{}/{}", ref_str, paired_alt);
            variant.mate_chr = Some(mate_chr_val);
            variant.mate_pos = Some(mate_pos);
            variant.is_single_breakend = false;
            // A set `mate_id` makes `display_allele()` render the literal allele
            // string rather than the SO term; gnomAD <BND> records carry no MATEID.
            if variant.mate_id.is_none() {
                let synthetic_id = variant
                    .id
                    .clone()
                    .unwrap_or_else(|| format!("{}:{}", chr, variant.start));
                variant.mate_id = Some(synthetic_id);
            }

            let single_alt = format!("{}.", ref_str);
            let mut sb_variant = variant.clone();
            sb_variant.alt_alleles = vec![single_alt.as_bytes().to_vec()];
            sb_variant.allele_string = format!("{}/{}", ref_str, single_alt);
            sb_variant.mate_chr = None;
            sb_variant.mate_pos = None;
            sb_variant.is_single_breakend = true;

            // The record's own breakend precedes its mates, as in Perl's
            // StructuralVariationOverlap.pm (`for ($vf, @$breakends)`).
            variants.push(sb_variant);
            variants.push(variant);
        } else if is_bnd && is_paired_bnd {
            let ref_str = std::str::from_utf8(ref_allele.as_bytes()).unwrap_or("N");
            // Orientation from the paired ALT: `]p]t` / `[p[t` put the ref base
            // after the break (".t"); `t[p[` / `t]p]` put it before ("t.").
            let single_breakend_alt = if alt.starts_with('[') || alt.starts_with(']') {
                format!(".{}", ref_str)
            } else {
                format!("{}.", ref_str)
            };
            let mut sb_variant = variant.clone();
            sb_variant.alt_alleles = vec![single_breakend_alt.as_bytes().to_vec()];
            sb_variant.allele_string = format!("{}/{}", ref_str, single_breakend_alt);
            sb_variant.mate_chr = None;
            sb_variant.mate_pos = None;
            sb_variant.is_single_breakend = true;
            variants.push(sb_variant);
            variants.push(variant);
        } else {
            variants.push(variant);
        }
    }

    Ok(variants)
}

pub(crate) fn is_paired_breakend_alt(alt: &str) -> bool {
    alt.contains('[') || alt.contains(']')
}

pub(crate) fn is_native_single_breakend_alt(alt: &str) -> bool {
    alt.len() > 1 && (alt.ends_with('.') || alt.starts_with('.')) && !is_paired_breakend_alt(alt)
}

/// Returns true if an ALT allele is a symbolic (structural variant) allele.
pub(crate) fn is_symbolic_alt(alt: &str) -> bool {
    if alt.is_empty() {
        return false;
    }
    // Perl VEP's SV regex `[<\[\]][^\*]+[>\]\[]` excludes `<*>`, treating it as
    // a non-SV allele at POS.
    if alt.starts_with('<') && alt.ends_with('>') && alt != "<*>" {
        return true;
    }
    if is_paired_breakend_alt(alt) {
        return true;
    }
    if is_native_single_breakend_alt(alt) {
        return true;
    }
    false
}

/// Classify a symbolic ALT allele into a VariantClass, with optional SVTYPE fallback.
pub(crate) fn classify_symbolic_alt(
    alt: &str,
    svtype: Option<&str>,
) -> vep_core::variant::VariantClass {
    use vep_core::variant::VariantClass;
    let upper = alt.to_uppercase();

    // Mobile elements must be checked before generic DEL/INS: <DEL:ME> starts with <DEL.
    if upper.starts_with("<INS:ME") {
        return VariantClass::MobileElementInsertion;
    }
    if upper.starts_with("<DEL:ME") {
        return VariantClass::MobileElementDeletion;
    }
    if upper.starts_with("<DEL") {
        return VariantClass::StructuralDeletion;
    }
    if upper.starts_with("<INS") {
        return VariantClass::StructuralInsertion;
    }
    if upper.starts_with("<DUP:TANDEM") {
        return VariantClass::TandemDuplication;
    }
    if upper.starts_with("<DUP") {
        return VariantClass::Duplication;
    }
    if upper.starts_with("<INV") {
        return VariantClass::Inversion;
    }
    if upper.starts_with("<CNV:TR") {
        return VariantClass::TandemRepeat;
    }
    // Perl VEP maps CN0 -> DEL ("deletion"), CN2 -> DUP ("duplication")
    if upper.starts_with("<CN") && !upper.starts_with("<CNV") {
        let inner = &upper[3..upper.len().saturating_sub(1)]; // strip "<CN" and ">"
        let num_str = inner.strip_prefix('=').unwrap_or(inner);
        if let Ok(cn) = num_str.parse::<u32>() {
            return match cn {
                0 => VariantClass::StructuralDeletion,
                2 => VariantClass::Duplication,
                _ => VariantClass::CopyNumberVariation,
            };
        }
    }
    if upper.starts_with("<CNV") || upper.starts_with("<CN") {
        return VariantClass::CopyNumberVariation;
    }
    if upper.starts_with("<CPX") {
        return VariantClass::ComplexStructural;
    }
    // <NON_REF> is a gVCF reference-confidence allele. Perl VEP treats it as an
    // SV with class_SO_term "NON_REF" and overlap-only consequences (no
    // feature_truncation/ablation), which CopyNumberVariation reproduces.
    if upper == "<NON_REF>" {
        return VariantClass::CopyNumberVariation;
    }
    if alt.contains('[')
        || alt.contains(']')
        || (alt.len() > 1 && (alt.ends_with('.') || alt.starts_with('.')))
    {
        return VariantClass::Translocation;
    }
    if let Some(svt) = svtype {
        match svt.to_uppercase().as_str() {
            "DEL" => return VariantClass::StructuralDeletion,
            "INS" => return VariantClass::StructuralInsertion,
            "DUP" => return VariantClass::Duplication,
            "INV" => return VariantClass::Inversion,
            "BND" => return VariantClass::Translocation,
            "CNV" => return VariantClass::CopyNumberVariation,
            _ => {}
        }
    }
    VariantClass::ComplexStructural
}

/// Whether VEP knows the structural variant type of a record: the port of
/// `Parser::get_SO_term` (ensembl-vep) applied to every ALT joined by `/`, with
/// the abbreviation looked up in `%SO_TERMS` (ensembl-variation `Utils/Config.pm`).
/// A record whose type is unknown is parsed and, in the default format,
/// annotated, but VEP's VCF and JSON writers leave it without consequences.
pub(crate) fn vep_supports_sv_type(joined_alts: &str, svtype: Option<&str>) -> bool {
    const SO_TERMS: [&str; 18] = [
        "INS",
        "INS_ME",
        "INS_ALU",
        "INS_HERV",
        "INS_LINE1",
        "INS_SVA",
        "DEL",
        "DEL_ME",
        "DEL_ALU",
        "DEL_HERV",
        "DEL_LINE1",
        "DEL_SVA",
        "TREP",
        "TDUP",
        "DUP",
        "CNV",
        "INV",
        "BND",
    ];
    let upper = joined_alts.to_ascii_uppercase();
    let follows_vcf44 = {
        let t = upper.strip_prefix('<').unwrap_or(&upper);
        ["DEL", "INS", "DUP", "INV", "CNV"]
            .iter()
            .any(|p| t.starts_with(p))
            || (t.starts_with("CN")
                && t[2..]
                    .strip_prefix('=')
                    .unwrap_or(&t[2..])
                    .starts_with(|c: char| c.is_ascii_digit()))
    };
    let kind = match svtype {
        Some(st) if !follows_vcf44 => st.to_ascii_uppercase(),
        _ => upper,
    };
    let abbrev = if let Some(pos) = kind.find("INS:ME").or_else(|| kind.find("DEL:ME")) {
        let base = &kind[pos..pos + 3];
        let rest = kind[pos + 6..].trim_start_matches(':');
        let element: String = rest
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric())
            .collect();
        let element = if element == "L1" {
            "LINE1".to_string()
        } else {
            element
        };
        let subtype = if ["ALU", "HERV", "LINE1", "SVA"].contains(&element.as_str()) {
            element
        } else {
            "ME".to_string()
        };
        format!("{base}_{subtype}")
    } else if kind.contains("DUP:TANDEM") {
        "TDUP".to_string()
    } else if kind.contains("CNV:TR") {
        "TREP".to_string()
    } else if kind.contains("CNV")
        || kind.contains("CN=")
        || kind
            .find("CN")
            .is_some_and(|i| kind[i + 2..].starts_with(|c: char| c.is_ascii_digit()))
    {
        "CNV".to_string()
    } else if kind.contains(['[', ']']) || kind.starts_with('.') || kind.ends_with('.') {
        "BND".to_string()
    } else if kind.starts_with('<') || kind.ends_with('>') {
        let stripped: String = kind.chars().filter(|c| *c != '<' && *c != '>').collect();
        stripped.split(':').next().unwrap_or("").to_string()
    } else {
        kind.clone()
    };
    SO_TERMS.contains(&abbrev.as_str())
}

/// Extract a value from the VCF INFO field by key.
/// The alternate allele's length in bases for a `<CNV:TR>` record: INFO/RB, else RUC times
/// the repeat unit's length (RUS's length, else RUL). The first value of each field is read,
/// the record having one tandem-repeat allele; a fractional RUC is truncated the way Perl's
/// `x` operator truncates its count (`Parser/VCF.pm`, `_expand_tandem_repeat_allele_string`).
fn tandem_repeat_alt_bases(info: &str) -> Option<u64> {
    let first = |key: &str| get_info_value(info, key).and_then(|v| v.split(',').next());
    if let Some(rb) = first("RB").and_then(|v| v.parse::<f64>().ok()) {
        return (rb >= 0.0).then_some(rb as u64);
    }
    let ruc = first("RUC")?.parse::<f64>().ok()?;
    let unit_len = match first("RUS") {
        Some(rus) if rus != "." => rus.len() as f64,
        _ => first("RUL")?.parse::<f64>().ok()?,
    };
    (ruc >= 0.0).then_some((ruc.trunc() * unit_len) as u64)
}

pub(crate) fn get_info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    if info == "." || info.is_empty() {
        return None;
    }
    for field in info.split(';') {
        if let Some((k, v)) = field.split_once('=') {
            if k == key {
                return Some(v);
            }
        } else if field == key {
            return Some("");
        }
    }
    None
}

/// Parse a VCF confidence interval INFO field (CIPOS or CIEND).
///
/// Format: two comma-separated integers, e.g., "-100,50".
/// Returns `Some((lo, hi))` on success, `None` if missing or malformed.
pub(crate) fn parse_ci_field(value: Option<&str>) -> Option<(i64, i64)> {
    let s = value?;
    let (a, b) = s.split_once(',')?;
    Some((a.parse().ok()?, b.parse().ok()?))
}

/// Parse the mate chromosome and position from a BND ALT notation.
///
/// VCF BND notation encodes the mate position in the ALT field:
/// - `N[chr:pos[`: forward to forward join
/// - `]chr:pos]N`: reverse to forward join
/// - `N]chr:pos]`: forward to reverse join
/// - `[chr:pos[N`: reverse to reverse join
///
/// Returns `(chromosome, position)` with the chromosome normalized (no `chr` prefix).
pub(crate) fn parse_bnd_mate(alt: &str) -> Option<(String, u64)> {
    // Find the chr:pos portion inside brackets. A BND ALT delimits it with a
    // matched pair of either `[` or `]`; `[` wins when both appear, and a single
    // unmatched delimiter is not a BND.
    let delim = if alt.contains('[') { '[' } else { ']' };
    let start = alt.find(delim)?;
    let end = alt[start + 1..].find(delim).map(|e| start + 1 + e)?;
    let bracket_content = &alt[start + 1..end];

    let (chr_raw, pos_str) = bracket_content.split_once(':')?;
    let pos: u64 = pos_str.parse().ok()?;

    let chr = vep_core::coordinate::normalize_chromosome(chr_raw);

    Some((chr, pos))
}

/// The eight fixed columns of a VCF data line (CHROM through INFO). The split
/// stops after INFO: FORMAT and the sample columns are never read, and a
/// multisample line can carry thousands of them.
fn split_fixed_fields(line: &str) -> anyhow::Result<[&str; 8]> {
    let mut parts = line.splitn(9, '\t');
    let mut fields = [""; 8];
    for field in fields.iter_mut() {
        match parts.next() {
            Some(part) => *field = part,
            None => bail!("VCF line has fewer than 8 fields"),
        }
    }
    Ok(fields)
}

pub(crate) fn vcf_line_is_non_variant(line: &str) -> bool {
    // ALT='.' is the canonical VCF representation of a non-variant record.
    let mut fields = line.split('\t');
    let _chrom = fields.next();
    let _pos = fields.next();
    let _id = fields.next();
    let _ref = fields.next();
    matches!(fields.next(), Some("."))
}

/// Compute VEP coordinates from already-trimmed alleles.
pub(crate) fn vep_coords_from_trimmed(
    pos: u64,
    prefix_trimmed: usize,
    trimmed_ref: &[u8],
) -> (u64, u64) {
    let start = pos + prefix_trimmed as u64;
    let end = if trimmed_ref == b"-" {
        // VEP insertion convention: end < start.
        start - 1
    } else {
        start + trimmed_ref.len() as u64 - 1
    };
    (start, end)
}

/// Trim shared prefix/suffix across REF and all ALT alleles together.
///
/// This matches Perl VEP behavior for multi-allelic records where all ALT
/// alleles share a common normalized coordinate window.
pub(crate) fn trim_multi_alleles(
    ref_allele: &str,
    alt_alleles: &[&str],
) -> (String, Vec<String>, usize) {
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes: Vec<&[u8]> = alt_alleles.iter().map(|a| a.as_bytes()).collect();

    // Perl VEP behavior (Parser/VCF.pm lines 294-336):
    // 1. Check if any ALT has a different length from REF ($is_indel flag)
    // 2. For same-length substitutions (MNPs): no trimming at all
    // 3. For indels: trim only the shared first base (no suffix trimming)
    let is_indel = alt_bytes.iter().any(|a| a.len() != ref_bytes.len());

    if !is_indel {
        let trimmed_alts = alt_alleles.iter().map(|a| a.to_string()).collect();
        return (ref_allele.to_string(), trimmed_alts, 0);
    }

    // Indel: only the shared first base is trimmed across all alleles; per-allele
    // suffix trimming is the secondary step in parse_vcf_line() (Parser.pm:834-855).
    let mut prefix_len = 0usize;
    if !ref_bytes.is_empty() {
        let first = ref_bytes[0].to_ascii_uppercase();
        if alt_bytes
            .iter()
            .all(|a| !a.is_empty() && a[0].to_ascii_uppercase() == first)
        {
            prefix_len = 1;
        }
    }

    let ref_after = &ref_bytes[prefix_len..];
    let trimmed_ref = if ref_after.is_empty() {
        "-".to_string()
    } else {
        String::from_utf8_lossy(ref_after).to_string()
    };

    let mut trimmed_alts = Vec::with_capacity(alt_bytes.len());
    for alt in &alt_bytes {
        let alt_after = &alt[prefix_len..];
        trimmed_alts.push(if alt_after.is_empty() {
            "-".to_string()
        } else {
            String::from_utf8_lossy(alt_after).to_string()
        });
    }

    (trimmed_ref, trimmed_alts, prefix_len)
}

/// Convert VCF padded alleles to VEP minimal representation.
///
/// VCF: POS=100, REF=ACGT, ALT=A -> VEP: start=101, end=103, ref=CGT, alt=-
/// VCF: POS=100, REF=A, ALT=ACGT -> VEP: start=101, end=100, ref=-, alt=CGT
/// VCF: POS=100, REF=A, ALT=G -> VEP: start=100, end=100, ref=A, alt=G
#[cfg(test)]
fn vcf_to_vep_alleles(
    ref_allele: &str,
    alt_allele: &str,
    pos: u64,
) -> (Vec<u8>, Vec<u8>, u64, u64) {
    let (trimmed_ref, trimmed_alt, prefix_trimmed) =
        trim_alleles(ref_allele.as_bytes(), alt_allele.as_bytes());
    let (start, end) = vep_coords_from_trimmed(pos, prefix_trimmed, &trimmed_ref);
    (trimmed_ref, trimmed_alt, start, end)
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::variant::VariantClass;

    #[test]
    fn test_vcf_to_vep_snv() {
        let (r, a, s, e) = vcf_to_vep_alleles("A", "G", 100);
        assert_eq!(r, b"A");
        assert_eq!(a, b"G");
        assert_eq!(s, 100);
        assert_eq!(e, 100);
    }

    #[test]
    fn test_vcf_to_vep_deletion() {
        let (r, a, s, e) = vcf_to_vep_alleles("ACGT", "A", 100);
        assert_eq!(r, b"CGT");
        assert_eq!(a, b"-");
        assert_eq!(s, 101);
        assert_eq!(e, 103);
    }

    #[test]
    fn test_vcf_to_vep_insertion() {
        let (r, a, s, e) = vcf_to_vep_alleles("A", "ACGT", 100);
        assert_eq!(r, b"-");
        assert_eq!(a, b"CGT");
        assert_eq!(s, 101);
        assert_eq!(e, 100);
    }

    #[test]
    fn test_vcf_to_vep_complex_indel() {
        let (r, a, s, e) = vcf_to_vep_alleles("ACG", "AT", 100);
        assert_eq!(r, b"CG");
        assert_eq!(a, b"T");
        assert_eq!(s, 101);
        assert_eq!(e, 102);
    }

    #[test]
    fn test_vcf_to_vep_prefix_and_suffix_trim() {
        let (r, a, s, e) = vcf_to_vep_alleles("ACGT", "ATGT", 100);
        assert_eq!(r, b"C");
        assert_eq!(a, b"T");
        assert_eq!(s, 101);
        assert_eq!(e, 101);
    }

    /// The gate actually elides the capture. Without this, a change that
    /// hardcoded `capture_raw_input = true` would pass every other test while
    /// silently restoring one line-length `String` allocation per variant.
    #[test]
    fn raw_input_not_captured_when_gated_off() {
        let line = "21\t25000100\trs123\tA\tG\t50\tPASS\t.";
        let variants = parse_vcf_line(line, false).unwrap();
        assert_eq!(variants.len(), 1);
        assert!(
            variants[0].raw_input.is_none(),
            "raw_input must be absent when capture is gated off"
        );
    }

    /// The VCF writer needs the verbatim source line, so capture must retain it
    /// byte-for-byte (it re-emits the record with a CSQ INFO field appended).
    #[test]
    fn raw_input_captured_verbatim_when_gated_on() {
        let line = "21\t25000100\trs123\tA\tG\t50\tPASS\tDP=42";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants[0].raw_input.as_deref(), Some(line));
    }

    /// Multi-allelic records split into one variant per ALT, and the gate must
    /// apply to every split product.
    #[test]
    fn raw_input_gate_applies_to_every_split_allele() {
        let line = "21\t25000100\t.\tA\tG,T\t50\tPASS\t.";
        let off = parse_vcf_line(line, false).unwrap();
        assert_eq!(off.len(), 2);
        assert!(off.iter().all(|v| v.raw_input.is_none()));

        let on = parse_vcf_line(line, true).unwrap();
        assert_eq!(on.len(), 2);
        assert!(on.iter().all(|v| v.raw_input.as_deref() == Some(line)));
    }

    /// Symbolic-allele (SV) records take a separate assignment site in the
    /// parser from explicit alleles, so the gate must hold on both.
    #[test]
    fn raw_input_gate_applies_to_symbolic_alleles() {
        let line = "21\t25000100\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=25001000";
        let off = parse_vcf_line(line, false).unwrap();
        assert!(!off.is_empty());
        assert!(off.iter().all(|v| v.raw_input.is_none()));

        let on = parse_vcf_line(line, true).unwrap();
        assert!(on.iter().all(|v| v.raw_input.as_deref() == Some(line)));
    }

    /// FORMAT and the sample columns are never read, so a multisample line
    /// must parse exactly as its eight-column prefix does, `raw_input` aside.
    #[test]
    fn multisample_line_parses_like_its_eight_field_prefix() {
        let prefix = "21\t25000100\trs123\tA\tG,T\t50\tPASS\tDP=42;AF=0.5";
        let mut line = format!("{prefix}\tGT:DP");
        for i in 0..3_000 {
            line.push_str(if i % 2 == 0 { "\t0/1:12" } else { "\t1/1:7" });
        }
        let from_prefix = parse_vcf_line(prefix, false).unwrap();
        let from_line = parse_vcf_line(&line, false).unwrap();
        assert_eq!(from_line.len(), 2);
        assert_eq!(
            serde_json::to_string(&from_line).unwrap(),
            serde_json::to_string(&from_prefix).unwrap()
        );
        let captured = parse_vcf_line(&line, true).unwrap();
        assert!(captured
            .iter()
            .all(|v| v.raw_input.as_deref() == Some(line.as_str())));
    }

    /// A line with an INFO column but no FORMAT is the eight-field minimum.
    #[test]
    fn eight_field_line_reads_info() {
        let line = "21\t25000100\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=25001000";
        let variants = parse_vcf_line(line, false).unwrap();
        assert_eq!(variants[0].sv_type.as_deref(), Some("DEL"));
        assert_eq!(variants[0].sv_end, Some(25001000));
    }

    #[test]
    fn seven_field_line_keeps_its_error_text() {
        let err = parse_vcf_line("21\t25000100\trs123\tA\tG\t50\tPASS", false).unwrap_err();
        assert_eq!(err.to_string(), "VCF line has fewer than 8 fields");
        let err = parse_vcf_line("", false).unwrap_err();
        assert_eq!(err.to_string(), "VCF line has fewer than 8 fields");
    }

    #[test]
    fn test_parse_vcf_line_snv() {
        let line = "21\t25000100\trs123\tA\tG\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].chr, "21");
        assert_eq!(variants[0].start, 25000100);
        assert_eq!(variants[0].end, 25000100);
        assert_eq!(variants[0].ref_allele, b"A");
        assert_eq!(variants[0].alt_allele(), b"G");
        assert_eq!(variants[0].id.as_deref(), Some("rs123"));
    }

    #[test]
    fn test_parse_vcf_line_multiallelic() {
        let line = "21\t25000100\t.\tA\tG,T\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].display_allele(), "G");
        assert_eq!(variants[1].display_allele(), "T");
        assert_eq!(variants[0].allele_index, 0);
        assert_eq!(variants[1].allele_index, 1);
    }

    #[test]
    fn test_parse_vcf_line_multiallelic_indel_shared_window() {
        let line = "1\t100\t.\tCTTT\tC,CTT\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        // ALT C: prefix C trimmed -> TTT/- (deletion, alt is "-" -> no secondary trim)
        assert_eq!(variants[0].start, 101);
        assert_eq!(variants[0].end, 103);
        assert_eq!(variants[0].ref_allele, b"TTT");
        assert_eq!(variants[0].alt_allele(), b"-");

        // ALT CTT: prefix C chopped -> TTT/TT and kept there: a multi-allelic
        // record is never minimised further, so VEP computes this allele's
        // positions and distances on the three-base window.
        assert_eq!(variants[1].start, 101);
        assert_eq!(variants[1].end, 103);
        assert_eq!(variants[1].ref_allele, b"TTT");
        assert_eq!(variants[1].alt_allele(), b"TT");
        assert!(!variants[1].minimised);
        assert_eq!(variants[1].uploaded_variation(), "1_101_CTTT/C/CTT");
        assert_eq!(variants[1].record_allele_string(), "TTT/-/TT");
    }

    #[test]
    fn test_parse_vcf_line_ignores_star_for_joint_trimming() {
        // * (spanning deletion) is filtered: Perl VEP does not annotate it.
        let line = "1\t100\t.\tATAAT\tA,*\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].start, 101);
        assert_eq!(variants[0].end, 104);
        assert_eq!(variants[0].ref_allele, b"TAAT");
        assert_eq!(variants[0].alt_allele(), b"-");
    }

    #[test]
    fn test_parse_vcf_line_no_alt() {
        let line = "21\t25000100\t.\tA\t.\t50\tPASS\t.";
        assert!(parse_vcf_line(line, true).is_err());
    }

    #[test]
    fn test_parse_vcf_line_deletion() {
        let line = "21\t25000100\t.\tACGT\tA\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].start, 25000101);
        assert_eq!(variants[0].end, 25000103);
        assert_eq!(variants[0].ref_allele, b"CGT");
        assert_eq!(variants[0].alt_allele(), b"-");
    }

    #[test]
    fn test_parse_vcf_line_mixed_symbolic_uses_svtype_override() {
        let line = "21\t1000\t.\tC\tG,<INV>\t.\t.\tSVTYPE=DEL;END=2000;SVLEN=-1000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            1,
            "mixed records should emit only the symbolic ALT"
        );

        let variant = &variants[0];
        assert_eq!(variant.allele_index, 1);
        assert_eq!(variant.alt_allele(), b"<INV>");
        assert_eq!(variant.variant_class, VariantClass::StructuralDeletion);
        assert!(variant.is_structural);
        assert_eq!(variant.start, 1001);
        assert_eq!(variant.sv_end, Some(2000));
    }

    #[test]
    fn test_parse_vcf_line_paired_bnd_dual_emission() {
        let line = "13\t500\tbnd_1\tA\tA[21:100[\t.\t.\tSVTYPE=BND;MATEID=bnd_2";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(paired.variant_class, VariantClass::Translocation);
        assert_eq!(paired.alt_allele(), b"A[21:100[");
        assert_eq!(paired.start, 501);
        assert_eq!(paired.end, 501);
        assert_eq!(paired.mate_chr.as_deref(), Some("21"));
        assert_eq!(paired.mate_pos, Some(100));
        assert!(!paired.is_single_breakend);

        let single = &variants[0];
        assert_eq!(single.variant_class, VariantClass::Translocation);
        assert_eq!(single.alt_allele(), b"A.");
        assert!(single.is_single_breakend);
        assert!(single.mate_chr.is_none());
        assert!(single.mate_pos.is_none());
    }

    #[test]
    fn test_parse_vcf_line_native_single_breakend_not_duplicated() {
        let line = "13\t500\tbnd_1\tA\tA.\t.\t.\tSVTYPE=BND";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);

        let variant = &variants[0];
        assert_eq!(variant.variant_class, VariantClass::Translocation);
        assert_eq!(variant.alt_allele(), b"A.");
        assert!(variant.is_single_breakend);
        assert!(variant.mate_chr.is_none());
        assert!(variant.mate_pos.is_none());
    }

    #[test]
    fn test_parse_vcf_line_symbolic_bnd_with_svlen_dual_emission() {
        // <BND> symbolic alleles with SVLEN (gnomAD format): Perl VEP computes
        // span from POS+1 to POS+1+SVLEN-1 and creates paired + single-breakend.
        let line = "21\t10870051\tgnomAD_BND\tN\t<BND>\t.\t.\tSVTYPE=BND;END=10870052;SVLEN=52221";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            2,
            "symbolic <BND> with SVLEN should emit 2 variants"
        );

        let paired = &variants[1];
        assert_eq!(paired.variant_class, VariantClass::Translocation);
        assert_eq!(paired.start, 10870052); // POS+1
        assert_eq!(paired.end, 10922272); // POS+1 + 52221 - 1
        assert!(!paired.is_single_breakend);
        assert_eq!(paired.mate_chr.as_deref(), Some("21"));
        assert_eq!(paired.mate_pos, Some(10922273)); // end + 1
        assert_eq!(paired.sv_type.as_deref(), Some("BND"));
        assert_eq!(paired.sv_end, Some(10922272));

        let single = &variants[0];
        assert_eq!(single.variant_class, VariantClass::Translocation);
        assert_eq!(single.alt_allele(), b"N.");
        assert!(single.is_single_breakend);
        assert!(single.mate_chr.is_none());
        assert!(single.mate_pos.is_none());
        assert_eq!(single.start, 10870052);
        assert_eq!(single.end, 10922272);
        assert_eq!(single.sv_end, Some(10922272));
    }

    #[test]
    fn test_parse_vcf_line_symbolic_bnd_without_svlen_single_variant() {
        // <BND> without SVLEN: END defaults to POS+1, no span to create dual emission.
        let line = "21\t100\ttest_bnd\tN\t<BND>\t.\t.\tSVTYPE=BND;END=101";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            1,
            "symbolic <BND> without SVLEN should emit 1 variant"
        );
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::Translocation);
        assert_eq!(v.start, 101); // POS+1
        assert_eq!(v.end, 101); // END clamped to max(sv_start)
    }

    #[test]
    fn test_parse_vcf_line_bracket_bnd_takes_dual_emission_path() {
        // Bracket BND notation takes the dual-emission path, not the symbolic one.
        let line = "13\t500\tbnd_1\tA\tA[21:100[\t.\t.\tSVTYPE=BND;MATEID=bnd_2;SVLEN=5000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(paired.alt_allele(), b"A[21:100[");
        assert_eq!(paired.start, 501); // POS+1
                                       // Bracket BND: END falls back to POS+abs(SVLEN)=5500, clamped to max(sv_start)
        assert_eq!(paired.end, 5500);
        assert!(!paired.is_single_breakend);

        let single = &variants[0];
        assert_eq!(single.alt_allele(), b"A.");
        assert!(single.is_single_breakend);
    }

    #[test]
    fn test_parse_vcf_line_symbolic_bnd_with_chr2_end2() {
        // gnomAD <BND> with CHR2/END2: Perl VCF.pm uses END2 as the mate position
        // instead of the SVLEN-derived span end. This converts <BND> to a standard
        // paired BND (e.g., N[21:42458165[).
        let line = "21\t15727643\tBND_21_56725\tN\t<BND>\t.\t.\tSVTYPE=BND;END=15727644;SVLEN=26730521;CHR2=21;END2=42458165";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            2,
            "symbolic <BND> with CHR2/END2 should emit 2 variants"
        );

        let paired = &variants[1];
        assert_eq!(paired.variant_class, VariantClass::Translocation);
        assert_eq!(paired.start, 15727644); // POS+1
        assert!(!paired.is_single_breakend);
        assert_eq!(paired.mate_chr.as_deref(), Some("21"));
        assert_eq!(paired.mate_pos, Some(42458165)); // END2, not end+1
        let alt_str = std::str::from_utf8(paired.alt_allele()).unwrap();
        assert!(
            alt_str.contains("42458165"),
            "paired ALT should use END2 position: {alt_str}"
        );

        let single = &variants[0];
        assert!(single.is_single_breakend);
        assert_eq!(single.alt_allele(), b"N.");
        assert!(single.mate_chr.is_none());
    }

    /// With CHR2 but no END2, the mate position is the raw `INFO/END`.
    ///
    /// Perl's breakpoint branch (VCF.pm:483-498) stashes `INFO/END` as
    /// `$incorrect_end`, deletes it so it cannot be read as a span end, and then
    /// resolves `my $breakend_pos = $info->{END2} || $incorrect_end;`. POS2 is
    /// never consulted, so a record carrying POS2 but not END2 still resolves to
    /// `INFO/END`.
    #[test]
    fn test_parse_vcf_line_symbolic_bnd_pos2_no_end2_uses_info_end() {
        let line = "21\t5000000\tBND_test\tN\t<BND>\t.\t.\tSVTYPE=BND;END=5000001;SVLEN=100000;CHR2=21;POS2=5100001";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(
            paired.mate_pos,
            Some(5000001),
            "Perl resolves END2 || INFO/END and never reads POS2 on this branch"
        );
        assert_eq!(paired.mate_chr.as_deref(), Some("21"));
    }

    #[test]
    fn test_parse_vcf_line_symbolic_bnd_inter_chromosomal_chr2() {
        let line =
            "21\t1000000\tBND_inter\tN\t<BND>\t.\t.\tSVTYPE=BND;END=1000001;SVLEN=500;CHR2=chrX;END2=50000000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(paired.mate_chr.as_deref(), Some("X")); // Normalized from chrX
        assert_eq!(paired.mate_pos, Some(50000000));
        // The allele string spells the mate chromosome as the input did.
        assert_eq!(paired.alt_allele(), b"N[chrX:50000000[");
        assert_eq!(paired.allele_string, "N/N[chrX:50000000[");
    }

    /// CHR2 with neither END2 nor POS2 resolves to the raw `INFO/END`.
    ///
    /// CHR2 and END present, END2 absent: the mate position is the raw INFO/END
    /// (5000001 for the record below), not POS + |SVLEN| (5010001).
    #[test]
    fn test_parse_vcf_line_symbolic_bnd_chr2_without_end2_uses_info_end() {
        let line =
            "21\t5000000\tBND_noend2\tN\t<BND>\t.\t.\tSVTYPE=BND;END=5000001;SVLEN=10000;CHR2=21";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(
            paired.mate_pos,
            Some(5000001),
            "Perl uses the raw INFO/END as the mate position, not POS + |SVLEN|"
        );
        assert_eq!(paired.mate_chr.as_deref(), Some("21"));
        let alt_str = std::str::from_utf8(paired.alt_allele()).unwrap();
        assert!(
            alt_str.contains("21:5000001"),
            "paired ALT must reference CHR2:END: {alt_str}"
        );
    }

    /// Guard: END2 still wins over `INFO/END` when both are present.
    ///
    /// `$info->{END2} || $incorrect_end` is an ordered fallback, not a swap: a
    /// record that carries END2 must keep using it, so the rule cannot degenerate
    /// into "always use INFO/END".
    #[test]
    fn test_parse_vcf_line_symbolic_bnd_end2_preferred_over_info_end() {
        let line = "21\t5000000\tBND_both\tN\t<BND>\t.\t.\tSVTYPE=BND;END=5000001;SVLEN=500;CHR2=21;END2=5900000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);

        let paired = &variants[1];
        assert_eq!(
            paired.mate_pos,
            Some(5900000),
            "END2 must take precedence over INFO/END"
        );
    }

    #[test]
    fn test_parse_vcf_line_non_ref_is_annotated() {
        // <NON_REF> is a gVCF reference-confidence allele that Perl VEP annotates
        // as an SV with overlap-only consequences, hence CopyNumberVariation.
        let line = "21\t100\t.\tA\t<NON_REF>\t.\t.\tEND=200";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::CopyNumberVariation);
        assert!(v.is_structural);
        assert_eq!(v.sv_end, Some(200));
        assert_eq!(v.display_allele(), "<NON_REF>");
    }

    #[test]
    fn vep_skip_marks_unsupported_types_only() {
        for (alts, svtype, supported) in [
            ("<DEL>", Some("DEL"), true),
            ("<INS:ME:ALU>", Some("INS"), true),
            ("<DEL:ME:L1>", None, true),
            ("<DUP:TANDEM>", None, true),
            ("<CNV:TR>", None, true),
            ("<CN0>", None, true),
            ("<CN0>/<CN2>", None, true),
            ("<CNV>", None, true),
            ("<INV>", None, true),
            ("<BND>", Some("BND"), true),
            ("N[21:100[", Some("BND"), true),
            ("<DUP:INT>", None, true),
            ("<CPX>", Some("CPX"), false),
            ("<NON_REF>", None, false),
            ("<TRA>", Some("TRA"), false),
        ] {
            assert_eq!(
                vep_supports_sv_type(alts, svtype),
                supported,
                "{alts} {svtype:?}"
            );
        }
        let cpx = parse_vcf_line("21\t100\t.\tA\t<CPX>\t.\t.\tSVTYPE=CPX;END=500", true).unwrap();
        assert!(cpx[0].vep_skip);
        let del = parse_vcf_line("21\t100\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=500", true).unwrap();
        assert!(!del[0].vep_skip);
    }

    #[test]
    fn test_parse_vcf_line_cpx_is_annotated() {
        // <CPX> variants are annotated by Perl VEP (despite vep_skip flag) with
        // overlap-only consequences. Allele displayed as "CPX".
        let line = "21\t100\t.\tA\t<CPX>\t.\t.\tSVTYPE=CPX;END=500";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::ComplexStructural);
        assert!(v.is_structural);
        assert_eq!(v.sv_end, Some(500));
        assert_eq!(v.display_allele(), "CPX");
    }

    #[test]
    fn test_parse_vcf_line_copy_number_classification() {
        // CN0 → DEL, CN2 → DUP, CN3 → generic CNV (matching Perl single-allele behavior)
        let mut cn0_variants =
            parse_vcf_line("21\t300\t.\tA\t<CN0>\t.\t.\tSVTYPE=CNV;END=800", true).unwrap();
        let cn0 = cn0_variants.remove(0);
        assert_eq!(cn0.variant_class, VariantClass::StructuralDeletion);

        let mut cn2_variants =
            parse_vcf_line("21\t300\t.\tA\t<CN2>\t.\t.\tSVTYPE=CNV;END=800", true).unwrap();
        let cn2 = cn2_variants.remove(0);
        assert_eq!(cn2.variant_class, VariantClass::Duplication);

        let mut cn3_variants =
            parse_vcf_line("21\t300\t.\tA\t<CN3>\t.\t.\tSVTYPE=CNV;END=800", true).unwrap();
        let cn3 = cn3_variants.remove(0);
        assert_eq!(cn3.variant_class, VariantClass::CopyNumberVariation);
    }

    #[test]
    fn test_parse_vcf_line_multi_allelic_cn_coalesced() {
        // Perl VEP joins "<CN0>,<CN2>" into "<CN0>/<CN2>", which matches no
        // single-allele regex and becomes one generic copy_number_variation SVF.
        let line = "21\t14504804\tDUP_gs_CNV\tC\t<CN0>,<CN2>\t.\t.\tSVTYPE=CNV;END=14530597";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            1,
            "multi-allelic CN should produce 1 variant, not 2"
        );
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::CopyNumberVariation);
        assert!(v.is_structural);
        assert_eq!(v.display_allele(), "copy_number_variation");
        // Alt allele must be <CNV> (not <CN0>) so cnv.rs parse_copy_number()
        // returns Generic, avoiding transcript_ablation/feature_elongation overcalls.
        assert_eq!(
            std::str::from_utf8(v.alt_allele()).unwrap(),
            "<CNV>",
            "coalesced CN alt should be <CNV>, not the first allele"
        );
    }

    #[test]
    fn test_parse_vcf_line_multi_allelic_cn_four_alleles() {
        let line = "21\t15244513\t.\tT\t<CN0>,<CN2>,<CN3>,<CN4>\t.\t.\tSVTYPE=CNV;END=15253042";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].variant_class, VariantClass::CopyNumberVariation);
        assert_eq!(
            std::str::from_utf8(variants[0].alt_allele()).unwrap(),
            "<CNV>",
            "coalesced 4-allele CN should have <CNV> alt"
        );
    }

    #[test]
    fn test_parse_vcf_line_single_cn0_still_deletion() {
        let line = "21\t300\t.\tA\t<CN0>\t.\t.\tSVTYPE=CNV;END=800";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].variant_class, VariantClass::StructuralDeletion);
    }

    #[test]
    fn test_parse_vcf_line_single_cn2_still_duplication() {
        let line = "21\t300\t.\tA\t<CN2>\t.\t.\tSVTYPE=CNV;END=800";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].variant_class, VariantClass::Duplication);
    }

    #[test]
    fn test_parse_vcf_line_too_few_fields() {
        let line = "21\t25000100\t.\tA";
        assert!(parse_vcf_line(line, true).is_err());
    }

    #[test]
    fn test_parse_vcf_line_tandem_repeat_uses_svlen_for_end() {
        // VCF.pm:502-512: $end = $start + max(@svlen) - 1 when $so_term =~ /tandem/
        let line = "21\t27253937\tsynth_cnv_tr_0001\tC\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;END=27256025;SVLEN=307";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::TandemRepeat);
        // sv_start = POS+1 = 27253938, sv_end = sv_start + SVLEN - 1 = 27254244
        assert_eq!(v.start, 27253938);
        assert_eq!(v.sv_end, Some(27254244));
        // Must not use the END field value
        assert_ne!(v.sv_end, Some(27256025));
    }

    #[test]
    fn test_parse_vcf_line_cnv_svlen_over_end() {
        // SVLEN=307 and END=27256025: SVLEN wins (BaseVCF4.pm get_end);
        // end = POS + abs(SVLEN) = 27253937 + 307 = 27254244
        let line = "21\t27253937\t.\tC\t<CNV>\t.\tPASS\tSVTYPE=CNV;END=27256025;SVLEN=307";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::CopyNumberVariation);
        assert_eq!(v.sv_end, Some(27254244));
    }

    #[test]
    fn test_tandem_repeat_alt_bases_reads_rb_then_ruc_times_unit() {
        // RB wins; without it RUC times the unit length, the unit read from RUS and from RUL
        // when RUS is missing; a fractional RUC truncates as Perl's `x` operator does.
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RB=42;RUC=7;RUS=CAG"),
            Some(42)
        );
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RUC=7;RUS=CAG"),
            Some(21)
        );
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RUC=7.9;RUS=CAG"),
            Some(21)
        );
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RUC=7;RUS=.;RUL=4"),
            Some(28)
        );
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RB=42,60;RUS=CAG,AT"),
            Some(42)
        );
        assert_eq!(tandem_repeat_alt_bases("SVTYPE=CNV;RUC=7"), None);
        assert_eq!(tandem_repeat_alt_bases("SVTYPE=CNV;SVLEN=300"), None);
    }

    #[test]
    fn test_parse_vcf_line_tandem_repeat_carries_alt_bases() {
        let line = "21\t27253937\t.\tC\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;END=27254237;SVLEN=300;RN=1;RUS=CAG;RUL=3;RUC=120;RB=360";
        let v = &parse_vcf_line(line, true).unwrap()[0];
        assert_eq!(v.variant_class, VariantClass::TandemRepeat);
        assert_eq!(v.tr_alt_bases, Some(360));
        // The field is read for tandem repeats only.
        let cnv = "21\t27253937\t.\tC\t<CNV>\t.\tPASS\tSVTYPE=CNV;END=27254237;SVLEN=300;RB=360";
        assert_eq!(parse_vcf_line(cnv, true).unwrap()[0].tr_alt_bases, None);
    }

    #[test]
    fn test_parse_vcf_line_tandem_repeat_svlen_only() {
        let line = "21\t27253937\t.\tC\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;SVLEN=500";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::TandemRepeat);
        // sv_start = POS+1 = 27253938, sv_end = 27253938 + 500 - 1 = 27254437
        assert_eq!(v.sv_end, Some(27254437));
    }

    #[test]
    fn test_parse_vcf_line_complex_substitution_suffix_trimmed() {
        // VCF: TAA -> GCTCA (3bp ref, 5bp alt: complex substitution)
        // Perl VEP trims trailing A: ref=TA, alt=GCTC, end=27256844
        let line = "21\t27256843\t.\tTAA\tGCTCA\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ref_allele, b"TA");
        assert_eq!(v.alt_alleles[0], b"GCTC");
        assert_eq!(v.start, 27256843);
        assert_eq!(v.end, 27256844); // 27256843 + 2 - 1
        assert_eq!(v.allele_string, "TA/GCTC");
        assert!(v.minimised);
    }

    #[test]
    fn test_parse_vcf_line_simple_deletion_no_extra_suffix_trim() {
        // VCF: ACGT -> A (simple deletion: prefix trim removes A, leaving CGT/-)
        // Already minimal after prefix trim, no suffix trimming needed
        let line = "1\t100\t.\tACGT\tA\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ref_allele, b"CGT");
        assert_eq!(v.alt_alleles[0], b"-");
        assert_eq!(v.start, 101);
        assert_eq!(v.end, 103);
    }

    #[test]
    fn test_parse_vcf_line_mnp_no_trimming() {
        // Same-length MNP: no trimming at all (ref.len() == alt.len())
        let line = "1\t100\t.\tAC\tGT\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ref_allele, b"AC");
        assert_eq!(v.alt_alleles[0], b"GT");
        assert_eq!(v.start, 100);
        assert_eq!(v.end, 101);
        assert!(!v.minimised);
    }

    #[test]
    fn test_parse_vcf_line_insertion_with_shared_suffix() {
        // VCF: AG -> TCCAG (2bp ref, 5bp alt)
        // Multi-allele trim: A vs T no match -> prefix_len=0, ref=AG, alt=TCCAG
        // Secondary minimisation: no prefix match, suffix G matches G -> trim,
        // then suffix A matches A -> trim -> ref=-, alt=TCC
        let line = "1\t100\t.\tAG\tTCCAG\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ref_allele, b"-");
        assert_eq!(v.alt_alleles[0], b"TCC");
        assert!(v.minimised);
    }

    /// A bi-allelic record whose REF is a strict prefix of its ALT must be
    /// re-minimised to a pure insertion.
    ///
    /// Perl runs two minimisation stages. `VEP/Parser/VCF.pm:328-335` strips one
    /// shared leading base, then `post_process_vfs` (`VEP/Parser.pm:834-856`) calls
    /// `minimise_alleles` -> `trim_sequences` (`Sequence.pm:1004-1035`), which is
    /// not gated on `--minimal` and trims from the left while both alleles are
    /// non-empty. With `empty_to_dash` set, `T/TT` therefore becomes `-/T`; Perl
    /// does not keep a base on the shorter allele.
    #[test]
    fn test_parse_vcf_line_biallelic_prefix_ref_minimises_to_insertion() {
        // AT -> ATT: stage 1 strips the shared leading A leaving T/TT, then
        // trim_sequences strips the shared leading T, emptying ref.
        let line = "1\t100\t.\tAT\tATT\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(
            v.ref_allele,
            b"-",
            "Perl's trim_sequences empties the ref and empty_to_dash maps it to \
             '-'; got ref={:?} alt={:?}",
            String::from_utf8_lossy(&v.ref_allele),
            String::from_utf8_lossy(&v.alt_alleles[0])
        );
        assert_eq!(v.alt_alleles[0], b"T");
        // VEP insertion convention: end == start - 1.
        assert_eq!(v.start, 102);
        assert_eq!(v.end, 101);
        assert!(v.minimised);
    }

    /// Guard: the same shape with a second ALT must not be re-minimised.
    ///
    /// `minimise_alleles` returns any allele string matching `.+/.+/.+` unchanged
    /// (`VEP/Parser.pm:951-957`), so a multi-allelic record never reaches
    /// `trim_sequences` and keeps its anchor base.
    #[test]
    fn test_parse_vcf_line_multiallelic_prefix_ref_keeps_anchor_base() {
        let line = "1\t100\t.\tAT\tATT,A\t50\tPASS\t.";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        let ins = &variants[0];
        assert_eq!(
            ins.ref_allele,
            b"T",
            "a 2-ALT record returns early at Parser.pm:951 and keeps its anchor \
             base; got ref={:?}",
            String::from_utf8_lossy(&ins.ref_allele)
        );
        assert_eq!(ins.alt_alleles[0], b"TT");
        assert_eq!(ins.start, 101);
        assert_eq!(ins.end, 101);
    }

    #[test]
    fn test_parse_ci_field() {
        assert_eq!(parse_ci_field(Some("-100,50")), Some((-100, 50)));
        assert_eq!(parse_ci_field(Some("0,0")), Some((0, 0)));
        assert_eq!(parse_ci_field(Some("-500,500")), Some((-500, 500)));
        assert_eq!(parse_ci_field(None), None);
        assert_eq!(parse_ci_field(Some("")), None);
        assert_eq!(parse_ci_field(Some("abc")), None);
        assert_eq!(parse_ci_field(Some("-100")), None);
    }

    #[test]
    fn test_parse_vcf_line_imprecise_sv_cipos_ciend() {
        let line =
            "21\t27253938\tid1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=27260000;CIPOS=-100,50;CIEND=-50,100";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ci_pos, Some((-100, 50)));
        assert_eq!(v.ci_end, Some((-50, 100)));
        assert!(v.is_structural);
        assert_eq!(v.sv_end, Some(27260000));
    }

    #[test]
    fn test_parse_vcf_line_sv_without_cipos_ciend() {
        let line = "21\t27253938\tid1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=27260000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.ci_pos, None);
        assert_eq!(v.ci_end, None);
    }

    #[test]
    fn test_parse_vcf_line_star_ref_block_not_filtered() {
        // Perl VEP does not filter <*>: the SV regex explicitly excludes it.
        // <*> is treated as a non-SV allele at position POS.
        let line = "21\t27253474\trefblock1\tT\t<*>\t.\tPASS\tEND=27255568;LEN=2094";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            1,
            "<*> should produce a variant, not be filtered"
        );
        let v = &variants[0];
        assert_eq!(v.start, 27253474);
        assert_eq!(v.end, 27253474); // single position, not SV span
        assert!(!v.is_structural, "<*> should not be structural");
        assert_eq!(v.alt_allele(), b"<*>");
        assert_eq!(v.allele_string, "T/<*>");
    }

    #[test]
    fn test_parse_vcf_line_star_allele_still_filtered() {
        // Plain * (spanning deletion) should still be filtered
        let line = "21\t100\t.\tA\t*\t.\tPASS\t.";
        assert!(parse_vcf_line(line, true).is_err());
    }

    #[test]
    fn test_parse_vcf_line_parid_as_mate_id() {
        // VCF 4.5 uses PARID instead of MATEID for BND partner identification;
        // PARID is read when MATEID is absent.
        let line = "21\t100\tbnd_1\tA\tA[21:200[\t.\t.\tSVTYPE=BND;PARID=bnd_2";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        let paired = &variants[1];
        assert_eq!(paired.mate_id.as_deref(), Some("bnd_2"));
        let single = &variants[0];
        assert_eq!(single.mate_id.as_deref(), Some("bnd_2"));
        assert!(single.is_single_breakend);
        assert_eq!(single.display_allele(), "A.");
    }

    #[test]
    fn test_parse_vcf_line_mateid_preferred_over_parid() {
        // When both MATEID and PARID are present, MATEID takes precedence
        let line = "21\t100\tbnd_1\tA\tA[21:200[\t.\t.\tSVTYPE=BND;MATEID=mate1;PARID=par1";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants[0].mate_id.as_deref(), Some("mate1"));
    }

    #[test]
    fn test_gnomad_grch38_del_positive_svlen() {
        // gnomAD v4.1 GRCh38 DEL: chr21, positive SVLEN, with END and CHR2
        // Real example: chr21:5033802 SVLEN=157 END=5033959
        let line = "chr21\t5033802\tgnomAD-SV_v3_DEL_chr21_test\tN\t<DEL>\t.\t.\tALGORITHMS=manta;CHR2=chr21;END=5033959;SVLEN=157;SVTYPE=DEL";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.chr, "21");
        assert_eq!(v.original_chr, "chr21");
        // sv_start = POS + 1 (DEL is span_type)
        assert_eq!(v.start, 5033803);
        // sv_end = POS + unsigned_abs(SVLEN) = 5033802 + 157 = 5033959
        assert_eq!(v.end, 5033959);
        assert_eq!(v.variant_class, VariantClass::StructuralDeletion);
        assert!(v.is_structural);
        assert_eq!(v.location(), "chr21:5033803-5033959");
        assert_eq!(v.display_allele(), "deletion");
    }

    #[test]
    fn test_gnomad_grch38_dup() {
        let line = "chr21\t10000000\tgnomAD_DUP_test\tN\t<DUP>\t.\t.\tCHR2=chr21;END=10000500;SVLEN=500;SVTYPE=DUP";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.chr, "21");
        assert_eq!(v.original_chr, "chr21");
        assert_eq!(v.start, 10000001);
        assert_eq!(v.end, 10000500);
        assert_eq!(v.variant_class, VariantClass::Duplication);
        assert_eq!(v.location(), "chr21:10000001-10000500");
        assert_eq!(v.display_allele(), "duplication");
    }

    #[test]
    fn test_gnomad_grch38_ins_me_alu() {
        let line = "chr21\t20000000\tgnomAD_ME_ALU_test\tN\t<INS:ME:ALU>\t.\t.\tCHR2=chr21;END=20000001;SVLEN=300;SVTYPE=INS:ME:ALU";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.chr, "21");
        // ME insertions are point-type symbolic alleles; `is_span_type()` is true for every SV class.
        assert_eq!(v.variant_class, VariantClass::MobileElementInsertion);
        assert_eq!(v.display_allele(), "Alu_insertion");
    }

    #[test]
    fn test_gnomad_grch38_bnd_svlen_minus1() {
        // BND with SVLEN=-1 and CHR2+END2: dual emission triggers on CHR2 plus
        // END2/POS2 even when SVLEN yields end == start, so 2 variants result.
        let line = "chr21\t10838490\tgnomAD_BND_test\tN\t<BND>\t.\t.\tCHR2=chr2;END=10838491;END2=89831875;SVLEN=-1;SVTYPE=BND";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            2,
            "BND SVLEN=-1 with CHR2+END2: should produce paired + single-breakend"
        );
        let paired = &variants[1];
        assert_eq!(paired.chr, "21");
        assert_eq!(paired.variant_class, VariantClass::Translocation);
        assert!(!paired.is_single_breakend);
        assert_eq!(paired.mate_chr.as_deref(), Some("2"));
        assert_eq!(paired.mate_pos, Some(89831875));
        let sb = &variants[0];
        assert_eq!(sb.chr, "21");
        assert_eq!(sb.variant_class, VariantClass::Translocation);
        assert!(sb.is_single_breakend);
        assert!(sb.mate_chr.is_none());
    }

    #[test]
    fn test_gnomad_grch38_bnd_large_svlen_with_chr2_end2() {
        let line = "chr21\t10838490\tgnomAD_BND_large_test\tN\t<BND>\t.\t.\tCHR2=chr2;END=10873831;END2=89831875;SVLEN=35340;SVTYPE=BND";
        let variants = parse_vcf_line(line, true).unwrap();
        // Large SVLEN: sv_start = 10838491, sv_end = 10838491 + 35340 - 1 = 10873830
        // end > start, so dual-emission gate passes -> 2 variants
        assert_eq!(variants.len(), 2);
        let paired = &variants[1];
        // CHR2/END2 present: mate_chr = normalize("chr2") = "2", mate_pos = 89831875
        assert_eq!(paired.mate_chr.as_deref(), Some("2"));
        assert_eq!(paired.mate_pos, Some(89831875));
        let single = &variants[0];
        assert!(single.is_single_breakend);
    }

    #[test]
    fn test_gnomad_grch38_ctx() {
        // <CTX> has no explicit handling in classify_symbolic_alt and SVTYPE=CTX
        // is in no match arm, so it lands in ComplexStructural.
        let line = "chr21\t10000000\tgnomAD_CTX_test\tN\t<CTX>\t.\t.\tCHR2=chr3;END=10000001;END2=50000000;SVLEN=-1;SVTYPE=CTX";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::ComplexStructural);
        // Perl VEP treats CTX as a translocation/BND; vep-rs does not.
    }
}
