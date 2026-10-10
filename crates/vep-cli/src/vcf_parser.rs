// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VCF line parsing and allele processing.
//!
//! Handles parsing of VCF data lines into `InputVariant`s, including
//! multi-allelic splitting, structural variant classification, BND
//! notation parsing, and allele trimming/normalization.
//!
//! Perl citations name modules of ensembl-vep release/116 (`Bio/EnsEMBL/VEP/...`);
//! `Sequence.pm`, `StructuralVariationFeature.pm` and `StructuralVariationOverlap.pm`
//! are in ensembl-variation release/116 and `BaseVCF4.pm` is in ensembl-io release/116.

use anyhow::{bail, Context};

use vep_core::allele::trim_alleles;
use vep_core::coordinate::{parse_coordinate, MAX_COORD};
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
///
/// Tandem-repeat lengths are read under Ensembl VEP's default `--max_sv_size`;
/// [`parse_vcf_line_with_max_sv_size`] takes the run's own limit.
pub fn parse_vcf_line(line: &str, capture_raw_input: bool) -> anyhow::Result<Vec<InputVariant>> {
    parse_vcf_line_with_max_sv_size(
        line,
        capture_raw_input,
        Some(crate::args::DEFAULT_MAX_SV_SIZE as u64),
    )
}

/// [`parse_vcf_line`] under the run's `--max_sv_size` (`None` for `-1`), which
/// decides whether a `<CNV:TR>` record's alternate allele is expanded to a length
/// or kept symbolic, as `_expand_tandem_repeat_allele_string` (`Parser/VCF.pm`)
/// decides under `max_sv_size`.
pub fn parse_vcf_line_with_max_sv_size(
    line: &str,
    capture_raw_input: bool,
    max_sv_size: Option<u64>,
) -> anyhow::Result<Vec<InputVariant>> {
    let fields = split_fixed_fields(line)?;

    let chr = fields[0].to_string();
    let pos =
        parse_coordinate(fields[1]).with_context(|| format!("Invalid POS: '{}'", fields[1]))?;
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
            // POS (its SV regex excludes `<*>`). `<CPX>` is not filtered either:
            // absent from `%SO_TERMS`, it is still annotated in tab output because
            // `vep_skip` is checked only by File-based annotation sources and
            // VCF/JSON output. `<NON_REF>` is not filtered: the record is written
            // as one intergenic row (`VariantClass::ReferenceBlock`).
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

    // A record with any symbolic allele is one structural variant (Parser/VCF.pm
    // `create_VariationFeatures` routes the whole line to
    // `create_StructuralVariationFeatures`), so the sequence path below fires
    // only when no symbolic allele is present; a mixed record's sequence ALTs
    // are alleles of the structural variant.
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

    // VEP's allele string for a structural variant record is every ALT of the
    // record joined by `/`, explicit alleles of a mixed record included
    // (Parser/VCF.pm `create_StructuralVariationFeatures`); the REF base it
    // prepends for a breakend record alone is added by `record_allele_string()`.
    let joined_alts = valid_alts
        .iter()
        .map(|(_, a)| *a)
        .collect::<Vec<_>>()
        .join("/");
    let svtype = get_info_value(info_field, "SVTYPE");
    // One structural variant record has one class and one SO term, read from the
    // joined ALTs with INFO/SVTYPE standing in for an ALT outside the VCF 4.4 forms
    // (Parser.pm `get_SO_term`, Parser/VCF.pm:472-475), so `N.` with `SVTYPE=DEL`
    // is a deletion; the term labels the row of every ALT, each ALT keeping its own
    // bytes. A record whose joined type is in no SO table is parsed with each ALT's
    // own class and the `vep_skip` mark; a breakend record keeps one variant per
    // bracket ALT.
    let record_term = sv_so_term_and_class(&sv_type_abbreviation(&joined_alts, svtype));
    let per_alt_classes: Vec<(usize, &str, vep_core::variant::VariantClass)> = match record_term {
        Some((_, class))
            if !symbolic_alts.is_empty()
                && class != vep_core::variant::VariantClass::Translocation =>
        {
            valid_alts
                .iter()
                .map(|&(idx, alt)| (idx, alt, class))
                .collect()
        }
        _ => symbolic_alts
            .iter()
            .map(|&(idx, alt)| {
                let class = match svtype {
                    // In a mixed record VEP reads the type from INFO/SVTYPE, since
                    // the joined ALTs open with a sequence allele.
                    Some(svtype) if is_mixed => {
                        classify_symbolic_alt(&format!("<{svtype}>"), Some(svtype))
                    }
                    _ => classify_symbolic_alt(alt, svtype),
                };
                (idx, alt, class)
            })
            .collect(),
    };
    let tandem_repeat_alt_count = |alts: &[(usize, &str, vep_core::variant::VariantClass)]| {
        alts.iter()
            .filter(|(_, a, _)| a.to_ascii_uppercase().starts_with("<CNV:TR"))
            .count()
    };
    let record_tandem_repeat_alts = tandem_repeat_alt_count(&per_alt_classes);

    for (k, &(allele_idx, alt, sv_class)) in per_alt_classes.iter().enumerate() {
        if sv_class == vep_core::variant::VariantClass::ReferenceBlock {
            tracing::warn!(
                "{}:{} {}: NON_REF is not a supported structural variant type; \
                 a reference-confidence block is written as one intergenic row",
                chr,
                pos,
                alt
            );
        }
        // Perl VEP uses POS+1 for all SV types (the VCF POS is the anchor base;
        // the actual structural event starts at POS+1).
        let sv_start = if sv_class.is_span_type() {
            pos + 1
        } else {
            pos
        };

        // Perl VEP (ensembl-io BaseVCF4.pm get_end()) universally prioritises
        // SVLEN over INFO/END for all SV types, and a symbolic record with
        // neither is the single base at get_start():
        //   if(defined($info->{SVLEN})) { $end = $self->get_start + abs($svlen) - 1; }
        //   elsif(defined($info->{END})) { $end = $info->{END}; }
        //   elsif(SVTYPE or a symbolic ALT) { $end = $self->get_start; }
        // Since get_start() = POS+1 for SVs, SVLEN gives end = POS + abs(SVLEN),
        // and SVLEN=0 or END=POS gives end = POS, one base before the start. A
        // comma-separated SVLEN is read as its first value (BaseVCF4.pm:304,
        // `(split(',',$info->{SVLEN}))[0]`); a tandem repeat or a tandem duplication
        // (`$so_term =~ /tandem/`, Parser/VCF.pm:502) takes the largest
        // (`:505-511`, `$end = $start + max(@svlen) - 1`), signed.
        let svlen_values: Vec<i64> = get_info_value(info_field, "SVLEN")
            .map(|v| v.split(',').filter_map(|s| s.parse().ok()).collect())
            .unwrap_or_default();
        let is_bnd = sv_class == vep_core::variant::VariantClass::Translocation;
        let is_tandem_repeat = sv_class == vep_core::variant::VariantClass::TandemRepeat;
        let is_tandem =
            is_tandem_repeat || sv_class == vep_core::variant::VariantClass::TandemDuplication;
        let svlen_span = if is_tandem {
            svlen_values
                .iter()
                .max()
                .map(|&len| pos.saturating_add(len.unsigned_abs()).min(MAX_COORD))
        } else {
            svlen_values
                .first()
                .map(|&len| pos.saturating_add(len.unsigned_abs()).min(MAX_COORD))
        };
        let info_end = get_info_value(info_field, "END");
        // A breakpoint record's INFO/END never ends its span: Parser/VCF.pm:487-489
        // deletes it and reads the end again, so a `<BND>` spans POS+1 to the SVLEN
        // end or is the single base at POS+1; the deleted END serves only as the
        // mate position below. Every other class reads END only without SVLEN
        // (BaseVCF4.pm:302-308), and an END it reads that is not a coordinate fails
        // the record, as Ensembl VEP skips it (Parser.pm `validate_vf`: `start >
        // end+1` for a negative END, the size test for one past the limit, "Invalid
        // start or end coordinate" for text). A tandem repeat or a tandem duplication
        // whose largest SVLEN is negative ends before its start and is skipped the
        // same way.
        let largest_svlen = svlen_values.iter().max().copied();
        let sv_end_val = if is_bnd {
            svlen_span.unwrap_or(sv_start)
        } else if is_tandem && largest_svlen.is_some_and(|len| len < 0) {
            bail!(
                "Invalid SVLEN for a tandem repeat or a tandem duplication: '{}'",
                largest_svlen.unwrap_or_default()
            );
        } else if let Some(span) = svlen_span {
            span
        } else {
            match info_end {
                Some(text) => {
                    parse_coordinate(text).with_context(|| format!("Invalid END: '{text}'"))?
                }
                None => sv_start,
            }
        };

        let is_symbolic_bnd = is_bnd && alt.starts_with('<');

        // An insertion whose END is its position or whose SVLEN is 0 stays the
        // pair `start = POS + 1, end = POS` Perl carries: the insertion between
        // the two bases, with the predicates of a sequence insertion and the
        // Location of its two flanks. Every other class keeps a span of at least
        // one base.
        let is_insertion_pair = matches!(
            sv_class,
            vep_core::variant::VariantClass::StructuralInsertion
                | vep_core::variant::VariantClass::MobileElementInsertion
        ) && sv_end_val.checked_add(1) == Some(sv_start);
        let sv_end_val = if is_insertion_pair {
            sv_end_val
        } else {
            sv_end_val.max(sv_start)
        };
        let mut variant = InputVariant::new(
            chr.clone(),
            sv_start,
            sv_end_val,
            ref_allele.as_bytes().to_vec(),
            alt.as_bytes().to_vec(),
        );
        variant.variant_class = sv_class;
        variant.is_structural = true;
        variant.vep_skip = record_term.is_none();
        variant.record_so_term = record_term.map(|(term, _)| term.to_string());
        variant.uploaded_allele_string = Some(joined_alts.clone());
        variant.sv_end = Some(sv_end_val);
        variant.sv_type = svtype.map(|s| s.to_string());
        variant.sv_len = svlen_values.first().copied();
        if is_tandem_repeat {
            variant.tr_alt_bases = tandem_repeat_alt_bases(
                info_field,
                tandem_repeat_alt_count(&per_alt_classes[..k]),
                record_tandem_repeat_alts,
                max_sv_size,
            );
        }
        variant.ci_pos = parse_ci_field(get_info_value(info_field, "CIPOS"));
        variant.ci_end = parse_ci_field(get_info_value(info_field, "CIEND"));
        if is_bnd {
            if let Some((mate_chr, mate_pos)) = parse_bnd_mate(alt) {
                variant.mate_chr = Some(mate_chr);
                variant.mate_pos = Some(mate_pos);
            }
        }
        variant.allele_index = allele_idx;
        if capture_raw_input {
            variant.raw_input = Some(line.to_string());
        }
        if id != "." && !id.is_empty() {
            variant.id = Some(id.to_string());
        }
        let is_paired_bnd = is_paired_breakend_alt(alt);
        // A translocation allele that names no mate is a single breakend: the native
        // `N.` and `.N` forms, and a symbolic `<BND>` without a mate, which Perl
        // holds as `N/<BND>` with no breakend string (Parser/VCF.pm:499,
        // `StructuralVariationFeature.pm` `_parse_breakends`).
        variant.is_single_breakend =
            is_bnd && (is_native_single_breakend_alt(alt) || is_symbolic_bnd);
        // Perl VEP creates two allele entries per BND record, the paired notation
        // ("T[21:33033339[") and the single-breakend form ("T." or ".T"), each
        // annotated independently against all overlapping transcripts.
        //
        // A symbolic `<BND>` names a mate only through INFO/CHR2 with END2, else the
        // deleted INFO/END (Parser/VCF.pm:491-497, `$info->{END2} || $incorrect_end`,
        // a Manta and gnomAD habit VCF 4.4 does not sanction); POS2 is never read,
        // and a span alone names none. The allele string keeps CHR2 as the input
        // spells it (`chr21`): VEP builds `N[chr21:pos[` from the INFO value verbatim,
        // while the mate lookup below uses the cache's chromosome name.
        let symbolic_mate = if is_symbolic_bnd {
            // Perl's `||` reads END only for an END2 that is absent, empty or `0`;
            // any other END2 is the mate position, and one that is not a coordinate
            // names no mate (`StructuralVariationFeature.pm:1276` matches `[0-9]+`,
            // so a signed `+5`, which `u64::from_str` accepts, names none), END or
            // not. A `CHR2=` or a bare `CHR2` names none either: `get_info`
            // (BaseVCF4.pm:582) leaves its value undef, so `defined $info->{CHR2}`
            // (Parser/VCF.pm:491) is false.
            let bracket_pos = |text: &str| {
                text.bytes()
                    .all(|b| b.is_ascii_digit())
                    .then(|| parse_coordinate(text))
                    .flatten()
            };
            let mate_pos = match get_info_value(info_field, "END2") {
                Some(end2) if !end2.is_empty() && end2 != "0" => bracket_pos(end2),
                _ => info_end.and_then(bracket_pos),
            };
            get_info_value(info_field, "CHR2")
                .filter(|chr2| !chr2.is_empty())
                .zip(mate_pos)
        } else {
            None
        };
        if let Some((chr2_raw, mate_pos)) = symbolic_mate {
            let ref_str = std::str::from_utf8(ref_allele.as_bytes()).unwrap_or("N");
            let paired_alt = format!("{}[{}:{}[", ref_str, chr2_raw, mate_pos);
            variant.alt_alleles = vec![paired_alt.as_bytes().to_vec()];
            variant.allele_string = format!("{}/{}", ref_str, paired_alt);
            variant.mate_chr = Some(vep_core::coordinate::normalize_chromosome(chr2_raw));
            variant.mate_pos = Some(mate_pos);
            variant.is_single_breakend = false;
            mark_paired_breakend_record(&mut variant, &chr, info_field);

            let single_alt = format!("{}.", ref_str);
            // The record's own breakend precedes its mates, as in Perl's
            // StructuralVariationOverlap.pm (`for ($vf, @$breakends)`).
            variants.push(own_breakend(&variant, ref_str, &single_alt));
            variants.push(variant);
        } else if is_bnd && is_paired_bnd {
            mark_paired_breakend_record(&mut variant, &chr, info_field);
            // One own breakend per record, ahead of its first mate: `_parse_breakends`
            // (`StructuralVariationFeature.pm:1274-1321`) names it once from every
            // bracket part of the allele string, `.t` when any mate is left-placed
            // (`]p]t`, `[p[t`) and `t.` otherwise, and `StructuralVariationOverlap::new`
            // annotates `($vf, @$breakends)`: that breakend, then one entry per bracket.
            let first_bracket_alt = !per_alt_classes[..k]
                .iter()
                .any(|(_, earlier, _)| is_paired_breakend_alt(earlier));
            if first_bracket_alt {
                let ref_str = std::str::from_utf8(ref_allele.as_bytes()).unwrap_or("N");
                let any_left_breakend = valid_alts
                    .iter()
                    .any(|(_, a)| a.starts_with('[') || a.starts_with(']'));
                let single_breakend_alt = if any_left_breakend {
                    format!(".{}", ref_str)
                } else {
                    format!("{}.", ref_str)
                };
                variants.push(own_breakend(&variant, ref_str, &single_breakend_alt));
            }
            variants.push(variant);
        } else {
            variants.push(variant);
        }
    }

    let line_name = input_line_name(&chr, pos, ref_allele, alt_field);
    for variant in &mut variants {
        variant.input_line_name = Some(line_name.clone());
    }

    Ok(variants)
}

/// The record's own breakend: a copy of its mate `variant` whose only allele is
/// `single_alt` (`t.` or `.t`) over `ref_str`, flagged single, with no mate coordinates.
fn own_breakend(variant: &InputVariant, ref_str: &str, single_alt: &str) -> InputVariant {
    let mut own = variant.clone();
    own.alt_alleles = vec![single_alt.as_bytes().to_vec()];
    own.allele_string = format!("{ref_str}/{single_alt}");
    own.mate_chr = None;
    own.mate_pos = None;
    own.is_single_breakend = true;
    own
}

/// Names the mate of a breakend record that has one: INFO/MATEID, else PARID, else
/// the record's ID or `chr:start` (gnomAD `<BND>` records carry no MATEID). A set
/// `mate_id` is how `display_allele()` and the breakend predicates tell the record's
/// own breakend, written `N.` or `.N` (`StructuralVariationFeature.pm`
/// `_parse_breakends` names it from the bracket or synthesised mate allele alone),
/// from a record whose only ALT is a single breakend, whose Allele is the bare class
/// term `chromosome_breakpoint` (`OutputFactory.pm:2029`,
/// `$vfoa->{breakend}->{string} || $svf->class_SO_term`). Ensembl VEP never reads
/// MATEID, so one on a record without a bracket mate names nothing and is not kept.
fn mark_paired_breakend_record(variant: &mut InputVariant, chr: &str, info: &str) {
    let mate_id = get_info_value(info, "MATEID")
        .or_else(|| get_info_value(info, "PARID"))
        .map(str::to_string)
        .or_else(|| variant.id.clone())
        .unwrap_or_else(|| format!("{}:{}", chr, variant.start));
    variant.mate_id = Some(mate_id);
}

/// The name VEP gives a VCF record whose ID column cannot name it
/// (OutputFactory.pm `VariationFeature_to_output_hash`, the VCF `_line` branch):
/// CHROM, POS, REF and ALT as written, the ALT commas as `/`, REF alone when the
/// ALT column is empty.
fn input_line_name(chr: &str, pos: u64, ref_allele: &str, alt_field: &str) -> String {
    if alt_field.is_empty() {
        format!("{chr}_{pos}_{ref_allele}")
    } else {
        format!("{chr}_{pos}_{ref_allele}/{}", alt_field.replace(',', "/"))
    }
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
    if upper.starts_with("<BND") {
        return VariantClass::Translocation;
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
    // A `<NON_REF>` with a supported SVTYPE took that class above, as in VEP, where
    // INFO/SVTYPE replaces an ALT outside the VCF 4.4 forms before the SO lookup
    // (Parser/VCF.pm `create_StructuralVariationFeatures`).
    if upper == "<NON_REF>" {
        return VariantClass::ReferenceBlock;
    }
    VariantClass::ComplexStructural
}

/// The structural variant type abbreviation VEP derives for a record: the port of
/// `Parser::get_SO_term` (ensembl-vep) applied to every ALT joined by `/`, with
/// INFO/SVTYPE replacing a joined string that opens outside the VCF 4.4 forms
/// (`Parser/VCF.pm`). Identical parts collapse to one (`<INS>/<INS>` reads `INS`),
/// and a string naming both DEL and DUP, `<DEL>/<DUP:TANDEM>` included, is `CNV`.
pub(crate) fn sv_type_abbreviation(joined_alts: &str, svtype: Option<&str>) -> String {
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
    let mut kind = match svtype {
        Some(st) if !follows_vcf44 => st.to_ascii_uppercase(),
        _ => upper,
    };
    if kind.contains(['/', ',']) {
        let parts: Vec<&str> = kind.split(['/', ',']).collect();
        fn normalized(part: &str) -> &str {
            part.trim().trim_start_matches('<').trim_end_matches('>')
        }
        let first = normalized(parts[0]);
        if parts.iter().all(|p| normalized(p) == first) {
            kind = parts[0].to_string();
        }
    }
    // The leftmost match of Perl's `/(INS|DEL):(ME):?([A-Z0-9]+)?/i` (Parser.pm:685).
    let mobile_element = ["INS:ME", "DEL:ME"]
        .into_iter()
        .filter_map(|needle| kind.find(needle))
        .min();
    if let Some(pos) = mobile_element {
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
    } else if kind.contains("DEL") && kind.contains("DUP") {
        "CNV".to_string()
    } else if kind.contains("DUP:TANDEM") {
        "TDUP".to_string()
    } else if kind.contains("CNV:TR") {
        "TREP".to_string()
    } else if kind.contains("CN=")
        || kind
            .find("CN")
            .is_some_and(|i| kind[i + 2..].starts_with(|c: char| c.is_ascii_digit()))
    {
        // Perl's `/^<?CN=?0>?$/` and `/^<?CN=?2>?$/` match the whole string, so a
        // list of copy numbers is `CNV`.
        let bare = kind.trim_start_matches('<').trim_end_matches('>');
        let copies = bare
            .strip_prefix("CN")
            .map(|c| c.strip_prefix('=').unwrap_or(c));
        match copies {
            Some("0") => "DEL".to_string(),
            Some("2") => "DUP".to_string(),
            _ => "CNV".to_string(),
        }
    } else if kind.contains("CNV") {
        "CNV".to_string()
    } else if kind.contains(['[', ']']) || kind.starts_with('.') || kind.ends_with('.') {
        "BND".to_string()
    } else if kind.starts_with('<') || kind.ends_with('>') {
        let stripped: String = kind.chars().filter(|c| *c != '<' && *c != '>').collect();
        stripped.split(':').next().unwrap_or("").to_string()
    } else {
        kind
    }
}

/// VEP's `%SO_TERMS` (ensembl-variation `Utils/Config.pm`): the SO term of a
/// record whose type abbreviation is one of its keys, with the class that
/// annotates the term, `None` for an unsupported type. Ensembl VEP 116 drops a
/// record of unknown type before annotation (`Parser.pm` `validate_svf`) and
/// writes nothing for it in any format; vep-rs parses and annotates it, carrying
/// its VCF line without consequences and writing an `input`-only JSON object
/// (`docs/intended-divergences.md`).
pub(crate) fn sv_so_term_and_class(
    abbrev: &str,
) -> Option<(&'static str, vep_core::variant::VariantClass)> {
    use vep_core::variant::VariantClass::*;
    Some(match abbrev {
        "INS" => ("insertion", StructuralInsertion),
        "INS_ME" => ("mobile_element_insertion", MobileElementInsertion),
        "INS_ALU" => ("Alu_insertion", MobileElementInsertion),
        "INS_HERV" => ("HERV_insertion", MobileElementInsertion),
        "INS_LINE1" => ("LINE1_insertion", MobileElementInsertion),
        "INS_SVA" => ("SVA_insertion", MobileElementInsertion),
        "DEL" => ("deletion", StructuralDeletion),
        "DEL_ME" => ("mobile_element_deletion", MobileElementDeletion),
        "DEL_ALU" => ("Alu_deletion", MobileElementDeletion),
        "DEL_HERV" => ("HERV_deletion", MobileElementDeletion),
        "DEL_LINE1" => ("LINE1_deletion", MobileElementDeletion),
        "DEL_SVA" => ("SVA_deletion", MobileElementDeletion),
        "TREP" => ("tandem_repeat", TandemRepeat),
        "TDUP" => ("tandem_duplication", TandemDuplication),
        "DUP" => ("duplication", Duplication),
        "CNV" => ("copy_number_variation", CopyNumberVariation),
        "INV" => ("inversion", Inversion),
        "BND" => ("chromosome_breakpoint", Translocation),
        _ => return None,
    })
}

/// The alternate allele's length in bases for the `ordinal`th `<CNV:TR>` allele of a
/// record carrying `alleles` of them: the length of the allele
/// `_expand_tandem_repeat_allele_string` (`Parser/VCF.pm:385-428`) builds. INFO/RN gives
/// each tandem-repeat allele its number of repeat sequences (one each when absent,
/// `:394-397`; an allele past the list's end has none, `:408`), consumed in order from
/// the RUS, RUC and RB lists (`:412-414`); a `.` unit is the one base `N` (`:413`); each
/// sequence's count is its RUC entry whenever the record carries RUC, else its RB entry
/// divided by the unit's length (`:414`), and the count is truncated the way Perl's `x`
/// operator truncates (`:419`). RUL is never read. `None` without RUS (`:385`), when an
/// entry the allele needs is missing, or is `.` in RN, RUC or RB, when the RN counts
/// through the allele sum past `usize::MAX` or its length past `u64::MAX`, and for every
/// allele of the record when any repeat sequence of the record is longer than
/// `max_sv_size` (`:403`, `:417-418`, `:428`), which keeps the symbolic reading: under
/// `None`, the flag's `-1`, no sequence is ever short enough.
fn tandem_repeat_alt_bases(
    info: &str,
    ordinal: usize,
    alleles: usize,
    max_sv_size: Option<u64>,
) -> Option<u64> {
    // Perl's `$info->{RUC} || ''` (`:394-400`): an empty or `0` value is no list.
    let list = |key: &str| -> Vec<&str> {
        get_info_value(info, key)
            .filter(|v| !v.is_empty() && *v != "0")
            .map(|v| v.split(',').collect())
            .unwrap_or_default()
    };
    let rus = list("RUS");
    if rus.is_empty() {
        return None;
    }
    let rn = list("RN");
    // One materialised allele per RN entry (`:408`): an allele past the list's end
    // has no sequences, and every allele has one when RN is absent.
    let sequences = |allele: usize| -> Option<usize> {
        if rn.is_empty() {
            Some(1)
        } else {
            rn.get(allele)?.parse().ok()
        }
    };
    let ruc = list("RUC");
    let rb = list("RB");
    // The `j`th repeat sequence's unit length and count, `(length($seq), $num)`.
    let sequence = |j: usize| -> Option<(f64, f64)> {
        let unit_len = match *rus.get(j)? {
            "." => 1.0,
            unit => unit.len() as f64,
        };
        let count: f64 = if ruc.is_empty() {
            rb.get(j)?.parse::<f64>().ok()? / unit_len
        } else {
            ruc.get(j)?.parse().ok()?
        };
        Some((unit_len, count))
    };
    // `$self->param('max_sv_size') || 5000` (`:403`): `-1` passes `||` unchanged, so
    // every sequence is oversized under it, and `0` falls to 5000.
    let cap = match max_sv_size {
        None => return None,
        Some(0) => 5000.0,
        Some(n) => n as f64,
    };
    let record_alleles = if rn.is_empty() { alleles } else { rn.len() };
    let mut total = 0usize;
    for allele in 0..record_alleles {
        total = total.checked_add(sequences(allele).unwrap_or(0))?;
    }
    for j in 0..total.min(rus.len()) {
        if sequence(j).is_some_and(|(unit_len, count)| unit_len * count > cap) {
            return None;
        }
    }
    let mut first = 0usize;
    for earlier in 0..ordinal {
        first = first.checked_add(sequences(earlier)?)?;
    }
    let last = first.checked_add(sequences(ordinal)?)?;
    let mut bases = 0u64;
    for j in first..last {
        let (unit_len, count) = sequence(j)?;
        let length = (count >= 0.0).then_some((count.trunc() * unit_len) as u64)?;
        bases = bases.checked_add(length)?;
    }
    Some(bases)
}

/// The value of `key` in a VCF INFO column, `""` for a flag, `None` when absent.
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
    // suffix trimming is the secondary step in parse_vcf_line() (Parser.pm:860-881).
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

    /// `tandem_repeat_alt_bases` for a record of `ordinal + 1` tandem-repeat alleles
    /// under Ensembl VEP's default `--max_sv_size`.
    fn alt_bases(info: &str, ordinal: usize) -> Option<u64> {
        let limit = Some(crate::args::DEFAULT_MAX_SV_SIZE as u64);
        tandem_repeat_alt_bases(info, ordinal, ordinal + 1, limit)
    }

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
        assert_eq!(variants[1].uploaded_variation(), "1_100_CTTT/C/CTT");
        assert_eq!(variants[1].record_allele_string(), "TTT/-/TT");
    }

    /// `Uploaded_variation` of every allele parsed from `line`.
    fn uploaded_variations(line: &str) -> Vec<String> {
        parse_vcf_line(line, false)
            .unwrap()
            .iter()
            .map(|v| v.uploaded_variation())
            .collect()
    }

    /// OutputFactory.pm `VariationFeature_to_output_hash`: the ID column names
    /// the record (the first of `;`-joined IDs, as Parser/VCF.pm takes it).
    #[test]
    fn uploaded_variation_is_the_id_when_the_record_carries_one() {
        assert_eq!(
            uploaded_variations("21\t9937918\trs867932594\tGGAGGGTCTGAATGTTTGTCCC\tG\t.\t.\t."),
            ["rs867932594"]
        );
        assert_eq!(
            uploaded_variations("21\t34971477\trs1;rs2\tC\tA\t.\t.\t."),
            ["rs1"]
        );
    }

    /// Ensembl VEP 116.2 on the GRCh37-hgvs corpus names `21 34971477 . C A`
    /// `21_34971477_C/A`.
    #[test]
    fn dot_id_names_an_snv_by_its_line() {
        assert_eq!(
            uploaded_variations("21\t34971477\t.\tC\tA\t.\t.\t."),
            ["21_34971477_C/A"]
        );
    }

    /// Ensembl VEP 116.2 names `21 10463009 . T TG` (OutputFactory.pm, the line-name branch)
    /// `21_10463009_T/TG` and `21 30677516 . TCACA T` (GRCh37-hgvs corpus)
    /// `21_30677516_TCACA/T`: the POS as written, not the anchor-trimmed start.
    #[test]
    fn dot_id_names_an_indel_by_its_line() {
        let variants = parse_vcf_line("21\t10463009\t.\tT\tTG\t.\t.\t.", false).unwrap();
        assert_eq!(variants[0].start, 10463010);
        assert_eq!(variants[0].uploaded_variation(), "21_10463009_T/TG");
        let variants = parse_vcf_line("21\t30677516\t.\tTCACA\tT\t.\t.\t.", false).unwrap();
        assert_eq!(variants[0].start, 30677517);
        assert_eq!(variants[0].allele_string, "CACA/-");
        assert_eq!(variants[0].uploaded_variation(), "21_30677516_TCACA/T");
    }

    /// Parser.pm `minimise_alleles` copies the record's line onto the minimised
    /// copy, so `21 100 . ATT AT`, minimised to `T/-` at 102, is named
    /// `21_100_ATT/AT` by OutputFactory.pm.
    #[test]
    fn minimised_copy_is_named_by_its_line() {
        let variants = parse_vcf_line("21\t100\t.\tATT\tAT\t.\t.\t.", false).unwrap();
        assert!(variants[0].minimised);
        assert_eq!(variants[0].start, 102);
        assert_eq!(variants[0].allele_string, "T/-");
        assert_eq!(variants[0].uploaded_variation(), "21_100_ATT/AT");
    }

    /// Ensembl VEP 116.2 on the GRCh37 corpus names both alleles of
    /// `21 9489854 . T TAA,TAC` `21_9489854_T/TAA/TAC`; the ALT column is taken
    /// as written, a `*` allele included (OutputFactory.pm turns its commas to `/`).
    #[test]
    fn dot_id_names_a_multi_allelic_record_by_its_line() {
        assert_eq!(
            uploaded_variations("21\t9489854\t.\tT\tTAA,TAC\t.\t.\t."),
            ["21_9489854_T/TAA/TAC", "21_9489854_T/TAA/TAC"]
        );
        assert_eq!(
            uploaded_variations("21\t100\t.\tA\tG,*\t.\t.\t."),
            ["21_100_A/G/*"]
        );
    }

    /// Ensembl VEP 116.2 on the GRCh37 corpus names `21 33867341 . C <CN2>`
    /// `21_33867341_C/<CN2>`: the REF base and the POS, where the structural
    /// variant itself starts at 33867342.
    #[test]
    fn dot_id_names_a_symbolic_record_by_its_line() {
        let variants =
            parse_vcf_line("21\t33867341\t.\tC\t<CN2>\t.\t.\tEND=34000018", false).unwrap();
        assert_eq!(variants[0].start, 33867342);
        assert_eq!(variants[0].uploaded_variation(), "21_33867341_C/<CN2>");
    }

    /// OutputFactory.pm: an ID equal to the `chr_start_alleles` name Parser.pm
    /// `validate_vf` generates for an empty ID column is not an identifier, so the
    /// record is named by its line like an ID-less one.
    #[test]
    fn an_id_equal_to_the_generated_name_yields_the_line_name() {
        assert_eq!(
            uploaded_variations("21\t30677516\t21_30677517_CACA/-\tTCACA\tT\t.\t.\t."),
            ["21_30677516_TCACA/T"]
        );
        assert_eq!(
            uploaded_variations("21\t33867341\t21_33867342_<CN2>\tC\t<CN2>\t.\t.\tEND=34000018"),
            ["21_33867341_C/<CN2>"]
        );
        // An ID that merely resembles the generated name is kept.
        assert_eq!(
            uploaded_variations("21\t30677516\t21_30677516_CACA/-\tTCACA\tT\t.\t.\t."),
            ["21_30677516_CACA/-"]
        );
    }

    /// A bracket breakend record's generated name is built from the one allele
    /// string Parser/VCF.pm `create_StructuralVariationFeatures` gives the record,
    /// REF joined to the bracket ALT (`N/N[21:9412900[`) at the start POS + 1, so
    /// an ID equal to it yields the line name on the record's own breakend row and
    /// on the mate row alike (OutputFactory.pm `VariationFeature_to_output_hash`
    /// names the record once); a native single breakend keeps its bare `.N`.
    #[test]
    fn breakend_id_equal_to_the_generated_name_yields_the_line_name_on_both_rows() {
        let variants = parse_vcf_line(
            "21\t9412500\t21_9412501_N/N[21:9412900[\tN\tN[21:9412900[\t.\t.\tSVTYPE=BND",
            false,
        )
        .unwrap();
        assert_eq!(variants.len(), 2);
        assert!(variants[0].is_single_breakend);
        for v in &variants {
            assert_eq!(v.record_allele_string(), "N/N[21:9412900[");
            assert_eq!(v.uploaded_variation(), "21_9412500_N/N[21:9412900[");
        }
        assert_eq!(
            uploaded_variations("21\t9412700\t21_9412701_.N\tN\t.N\t.\t.\tSVTYPE=BND"),
            ["21_9412700_N/.N"]
        );
        // An ID that merely resembles the generated name (no REF) is kept on both rows.
        assert_eq!(
            uploaded_variations(
                "21\t9412500\t21_9412501_N[21:9412900[\tN\tN[21:9412900[\t.\t.\tSVTYPE=BND"
            ),
            ["21_9412501_N[21:9412900[", "21_9412501_N[21:9412900["]
        );
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

    /// A sequence ALT beside a symbolic one: VEP 116.2 reads the record's type from
    /// INFO/SVTYPE (the joined ALTs open with a sequence allele) and writes one row
    /// per ALT, each with the class term as its Allele (the `synth_multi_mix` records
    /// of `tests/golden/116/GRCh37`, `deletion` twice per transcript).
    #[test]
    fn test_parse_vcf_line_mixed_symbolic_uses_svtype_override() {
        let line = "21\t1000\t.\tC\tG,<INV>\t.\t.\tSVTYPE=DEL;END=2000;SVLEN=-1000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(
            variants.len(),
            2,
            "one variant per ALT, the sequence ALT included"
        );
        for (variant, idx, alt) in [
            (&variants[0], 0, b"G" as &[u8]),
            (&variants[1], 1, b"<INV>"),
        ] {
            assert_eq!(variant.allele_index, idx);
            assert_eq!(variant.alt_allele(), alt);
            assert_eq!(variant.variant_class, VariantClass::StructuralDeletion);
            assert_eq!(variant.display_allele(), "deletion");
            assert!(variant.is_structural);
            assert!(!variant.vep_skip);
            assert_eq!(variant.start, 1001);
            assert_eq!(variant.sv_end, Some(2000));
            assert_eq!(variant.uploaded_allele_string.as_deref(), Some("G/<INV>"));
        }
    }

    /// A record whose ALTs are all sequence alleles is never a structural variant,
    /// whatever INFO/SVTYPE says (Parser/VCF.pm `create_VariationFeatures` routes a
    /// line to the structural path only on a symbolic or breakend ALT).
    #[test]
    fn all_sequence_multi_allelic_record_ignores_svtype() {
        let line = "21\t1000\t.\tC\tG,T\t.\t.\tSVTYPE=DEL;END=2000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        for v in &variants {
            assert_eq!(v.variant_class, VariantClass::Snv);
            assert!(!v.is_structural);
        }
    }

    /// A mixed record without INFO/SVTYPE has no type VEP knows (`G/<INV>` is in no
    /// SO table), so the symbolic ALT alone is parsed, with the `vep_skip` mark.
    #[test]
    fn test_parse_vcf_line_mixed_symbolic_without_svtype_is_unsupported() {
        let line = "21\t1000\t.\tC\tG,<INV>\t.\t.\tEND=2000";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].alt_allele(), b"<INV>");
        assert_eq!(variants[0].variant_class, VariantClass::Inversion);
        assert!(variants[0].vep_skip);
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
        // One allele, no bracket mate: VEP prints the class term as its Allele.
        assert_eq!(variant.display_allele(), "chromosome_breakpoint");
    }

    /// INFO/MATEID on a record whose only ALT is a single breakend names no mate:
    /// Ensembl VEP never reads MATEID, `_parse_breakends`
    /// (`StructuralVariationFeature.pm:1279-1283`) finds no bracket in `.N`, and the
    /// row's Allele falls to `class_SO_term` (`OutputFactory.pm:2029`),
    /// `chromosome_breakpoint`, as for the record without MATEID.
    #[test]
    fn single_breakend_with_mateid_is_parsed_as_without() {
        let with =
            parse_vcf_line("21\t100\t.\tN\t.N\t.\t.\tSVTYPE=BND;MATEID=bnd_x", false).unwrap();
        let without = parse_vcf_line("21\t100\t.\tN\t.N\t.\t.\tSVTYPE=BND", false).unwrap();
        assert_eq!(format!("{with:?}"), format!("{without:?}"));
        assert_eq!(with.len(), 1);
        assert!(with[0].is_single_breakend);
        assert_eq!(with[0].mate_id, None);
        assert_eq!(with[0].display_allele(), "chromosome_breakpoint");
    }

    /// A bracket breakend record names its own breakend `N.` or `.N` whether or not
    /// INFO carries a MATEID (`StructuralVariationFeature.pm` `_parse_breakends`
    /// builds the string from the bracket alleles alone): Ensembl VEP 116.2 writes
    /// `N.` on both rows of `bnd_genic_with_intergenic_mate` (no MATEID) in the
    /// GRCh38-release116 golden corpus.
    #[test]
    fn bracket_breakend_without_mateid_names_its_own_allele_n_dot() {
        let line =
            "21\t30425650\tbnd_genic_with_intergenic_mate\tN\tN[21:11600000[\t.\t.\tSVTYPE=BND";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        assert!(variants[0].is_single_breakend);
        assert_eq!(variants[0].display_allele(), "N.");
        assert_eq!(variants[1].display_allele(), "N[21:11600000[");

        let left =
            parse_vcf_line("21\t11500000\t.\tN\t]21:11000000]N\t.\t.\tSVTYPE=BND", true).unwrap();
        assert_eq!(left[0].display_allele(), ".N");
    }

    /// A record listing two bracket ALTs has one own breakend and two mates, three
    /// variants: `_parse_breakends` (`StructuralVariationFeature.pm:1274-1321`) pushes
    /// one breakend per bracket part of `N/N[21:11500000[/N]21:11600000]` and names
    /// the record's own string once, `N.` when every mate is right-placed and `.N`
    /// when any is left-placed (`:1298`, `:1321`), and `StructuralVariationOverlap::new`
    /// annotates `($vf, @$breakends)`, the own breakend first (`:90`).
    #[test]
    fn two_bracket_alts_yield_one_own_breakend_and_two_mates() {
        let line = "21\t11000000\tmb2\tN\tN[21:11500000[,N]21:11600000]\t.\t.\tSVTYPE=BND";
        let variants = parse_vcf_line(line, true).unwrap();
        let alleles: Vec<String> = variants.iter().map(|v| v.display_allele()).collect();
        assert_eq!(alleles, ["N.", "N[21:11500000[", "N]21:11600000]"]);
        assert!(variants[0].is_single_breakend);
        assert_eq!(variants[0].allele_index, 0);
        assert_eq!(variants[1].mate_pos, Some(11500000));
        assert_eq!(variants[2].mate_pos, Some(11600000));
        assert_eq!(variants[2].allele_index, 1);
        // The record string of every variant is the whole record, REF and both ALTs,
        // as `Parser/VCF.pm` joins the ALTs (`:469`) and prepends REF (`:499`); JSON
        // `allele_string` prints it (`OutputFactory/JSON.pm:241`).
        for v in &variants {
            assert_eq!(v.record_allele_string(), "N/N[21:11500000[/N]21:11600000]");
        }

        let mixed = "21\t11000000\tmb_mixed\tN\tN[21:11500000[,]21:11600000]N\t.\t.\tSVTYPE=BND";
        let variants = parse_vcf_line(mixed, true).unwrap();
        let alleles: Vec<String> = variants.iter().map(|v| v.display_allele()).collect();
        assert_eq!(alleles, [".N", "N[21:11500000[", "]21:11600000]N"]);
    }

    /// A `<BND>` with SVLEN and no CHR2 names no mate: Parser/VCF.pm:491-497 builds
    /// a bracket mate only under `defined $info->{CHR2}`, so the record stays
    /// `N/<BND>` (`:499`) spanning POS+1 to POS+|SVLEN| (BaseVCF4.pm:302-305) with
    /// no breakend string, one allele whose Allele is `chromosome_breakpoint`
    /// (OutputFactory.pm:2029); Ensembl VEP 116.2 writes `21:10870052-10922272
    /// chromosome_breakpoint` for it. INFO/END is deleted on a breakpoint record
    /// (`:487-489`) and never ends the span.
    #[test]
    fn symbolic_bnd_with_svlen_and_no_chr2_is_one_single_breakend_over_the_span() {
        let line = "21\t10870051\tgnomAD_BND\tN\t<BND>\t.\t.\tSVTYPE=BND;END=10870052;SVLEN=52221";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::Translocation);
        assert_eq!((v.start, v.end), (10870052, 10922272));
        assert_eq!(v.sv_end, Some(10922272));
        assert!(v.is_single_breakend);
        assert_eq!(v.mate_id, None);
        assert!(v.mate_chr.is_none() && v.mate_pos.is_none());
        assert_eq!(v.alt_allele(), b"<BND>");
        assert_eq!(v.display_allele(), "chromosome_breakpoint");
        assert_eq!(v.record_allele_string(), "N/<BND>");
    }

    /// A breakpoint record's INFO/END never ends its span (Parser/VCF.pm:487-489
    /// deletes it before the end is read again): `<BND>` with `END=30425750` and a
    /// bracket breakend with the same END are the single base at POS+1, as Ensembl
    /// VEP 116.2 writes `21:30425651` for both; with CHR2 the deleted END is the
    /// mate position (`:493`), the span still the point.
    #[test]
    fn a_breakpoint_record_s_end_never_ends_its_span() {
        let symbolic = parse_vcf_line(
            "21\t30425650\tend_span\tN\t<BND>\t.\t.\tSVTYPE=BND;END=30425750",
            true,
        )
        .unwrap();
        assert_eq!(symbolic.len(), 1);
        assert_eq!((symbolic[0].start, symbolic[0].end), (30425651, 30425651));
        assert_eq!(symbolic[0].display_allele(), "chromosome_breakpoint");
        let bracket = parse_vcf_line(
            "21\t30425650\tbracket_end\tN\tN[21:30425700[\t.\t.\tSVTYPE=BND;END=30425750",
            true,
        )
        .unwrap();
        assert_eq!(bracket.len(), 2);
        for v in &bracket {
            assert_eq!((v.start, v.end), (30425651, 30425651));
        }
        let chr2_end = parse_vcf_line(
            "21\t30425650\tchr2_far_end\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=2;END=40000000",
            true,
        )
        .unwrap();
        assert_eq!(chr2_end.len(), 2);
        assert_eq!((chr2_end[1].start, chr2_end[1].end), (30425651, 30425651));
        assert_eq!(chr2_end[1].mate_chr.as_deref(), Some("2"));
        assert_eq!(chr2_end[1].mate_pos, Some(40000000));
        assert_eq!(chr2_end[1].display_allele(), "N[2:40000000[");
    }

    /// A `<BND>` names a mate through CHR2 with END2, else the deleted INFO/END
    /// (Parser/VCF.pm:493, `$info->{END2} || $incorrect_end`); POS2 is never read
    /// and an END2 that is not a coordinate names nothing, so Ensembl VEP 116.2
    /// writes `chromosome_breakpoint` with no mate row for `CHR2=2;POS2=40000000`
    /// and for `CHR2=21;END2=-7`. Perl's `||` falls to END only for an END2 that is
    /// empty or `0`: `END2=-7` beside `END=30425700` still names no mate, while
    /// `END2=0` reads END. A signed `+30425700` as END2 or END names none, as
    /// `_parse_breakends` reads the position as `[0-9]+`
    /// (`StructuralVariationFeature.pm:1276`). A `CHR2=` or a bare `CHR2` names
    /// none, as `split('=')` leaves the value undef and `defined $info->{CHR2}`
    /// (`:491`) is false: Ensembl VEP 116.2 writes `chromosome_breakpoint` for
    /// `CHR2=;END2=30425800`.
    #[test]
    fn pos2_or_an_unreadable_end2_names_no_mate() {
        for info in [
            "SVTYPE=BND;CHR2=2;POS2=40000000",
            "SVTYPE=BND;CHR2=21;END2=-7",
            "SVTYPE=BND;CHR2=21;END2=-7;END=30425700",
            "SVTYPE=BND;CHR2=21;END2=abc;END=30425700",
            "SVTYPE=BND;CHR2=21;END2=+30425700",
            "SVTYPE=BND;CHR2=21;END=+30425700",
            "SVTYPE=BND;CHR2=;END2=30425800",
            "SVTYPE=BND;CHR2;END2=30425800",
        ] {
            let variants =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t<BND>\t.\t.\t{info}"), true).unwrap();
            assert_eq!(variants.len(), 1, "{info}");
            assert!(variants[0].is_single_breakend, "{info}");
            assert_eq!(variants[0].mate_pos, None, "{info}");
            assert_eq!(
                variants[0].display_allele(),
                "chromosome_breakpoint",
                "{info}"
            );
        }
        let with_end2 = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700",
            true,
        )
        .unwrap();
        assert_eq!(with_end2.len(), 2);
        assert_eq!(with_end2[1].mate_pos, Some(30425700));
        let zero_end2 = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=0;END=30425700",
            true,
        )
        .unwrap();
        assert_eq!(zero_end2.len(), 2);
        assert_eq!(zero_end2[1].mate_pos, Some(30425700));
    }

    /// A comma-separated SVLEN is read as its first value for the span
    /// (BaseVCF4.pm:304, `(split(',',$info->{SVLEN}))[0]`, so `-300,400` ends at
    /// POS+300) and as its largest for a tandem repeat (Parser/VCF.pm:505-511,
    /// `$end = $start + max(@svlen) - 1`, so `300,400` ends at POS+400), as Ensembl
    /// VEP 116.2 writes `21:30425601-30425900` and `21:30425601-30426000`.
    #[test]
    fn a_comma_separated_svlen_reads_its_first_value_and_a_tandem_repeat_its_largest() {
        let cnv = parse_vcf_line(
            "21\t30425600\tx\tN\t<DEL>,<DUP>\t.\t.\tSVTYPE=CNV;SVLEN=-300,400",
            true,
        )
        .unwrap();
        assert_eq!(cnv.len(), 2);
        for v in &cnv {
            assert_eq!((v.start, v.end), (30425601, 30425900));
            assert_eq!(v.sv_len, Some(-300));
        }
        let repeat = parse_vcf_line(
            "21\t30425600\tx\tN\t<CNV:TR>,<CNV:TR>\t.\t.\tSVTYPE=CNV;SVLEN=300,400;RUS=AT,AT;RUC=100,50",
            true,
        )
        .unwrap();
        assert_eq!(repeat.len(), 2);
        for v in &repeat {
            assert_eq!((v.start, v.end), (30425601, 30426000));
        }
    }

    /// A tandem duplication reads SVLEN as a tandem repeat does: Parser/VCF.pm:502
    /// routes every record whose SO term matches `/tandem/` (`tandem_duplication`
    /// included) to `$end = $start + max(@svlen) - 1` (`:505-511`), so Ensembl VEP
    /// 116.2 writes `21:30425601-30426000 tandem_duplication` for `<DUP:TANDEM>` with
    /// `SVLEN=300,400` and for `<DUP:TANDEM>,<DUP:TANDEM>` with the same SVLEN, and
    /// skips `SVLEN=-300` (`start > end+1`) with or without END, while a `<DUP>` with
    /// `SVLEN=300,400` or `SVLEN=-300` is `21:30425601-30425900 duplication`.
    #[test]
    fn a_tandem_duplication_reads_svlen_as_a_tandem_repeat_does() {
        let tandem = parse_vcf_line(
            "21\t30425600\tx\tN\t<DUP:TANDEM>\t.\t.\tSVTYPE=DUP;SVLEN=300,400",
            true,
        )
        .unwrap();
        assert_eq!(tandem.len(), 1);
        assert_eq!(tandem[0].variant_class, VariantClass::TandemDuplication);
        assert_eq!((tandem[0].start, tandem[0].end), (30425601, 30426000));
        assert_eq!(tandem[0].sv_len, Some(300));
        let two = parse_vcf_line(
            "21\t30425600\tx\tN\t<DUP:TANDEM>,<DUP:TANDEM>\t.\t.\tSVTYPE=DUP;SVLEN=300,400",
            true,
        )
        .unwrap();
        assert_eq!(two.len(), 2);
        for v in &two {
            assert_eq!(v.variant_class, VariantClass::TandemDuplication);
            assert_eq!((v.start, v.end), (30425601, 30426000));
        }
        for info in [
            "SVTYPE=DUP;SVLEN=-300",
            "SVTYPE=DUP;SVLEN=-300;END=30425900",
        ] {
            let err = parse_vcf_line(
                &format!("21\t30425600\tx\tN\t<DUP:TANDEM>\t.\t.\t{info}"),
                true,
            )
            .unwrap_err()
            .to_string();
            assert!(err.contains("Invalid SVLEN"), "{info}: {err}");
        }
        for info in ["SVTYPE=DUP;SVLEN=300,400", "SVTYPE=DUP;SVLEN=-300"] {
            let dup =
                parse_vcf_line(&format!("21\t30425600\tx\tN\t<DUP>\t.\t.\t{info}"), true).unwrap();
            assert_eq!(dup[0].variant_class, VariantClass::Duplication, "{info}");
            assert_eq!((dup[0].start, dup[0].end), (30425601, 30425900), "{info}");
        }
    }

    /// `N.` with `SVTYPE=DEL` is a deletion: the ALT is outside the VCF 4.4 forms, so
    /// INFO/SVTYPE is the type (Parser/VCF.pm:472-475) and `get_SO_term` reads
    /// `deletion`; Ensembl VEP 116.2 writes `21:30425601-30425899 deletion` with
    /// `frameshift_variant,feature_truncation`. A `<BND>` without SVTYPE is a
    /// breakpoint from its own ALT.
    #[test]
    fn svtype_names_the_class_of_a_dot_alt_and_a_bare_bnd_is_a_breakpoint() {
        let deletion = parse_vcf_line(
            "21\t30425600\tgnomad_del\tN\tN.\t.\t.\tSVTYPE=DEL;END=30425899",
            true,
        )
        .unwrap();
        assert_eq!(deletion.len(), 1);
        assert_eq!(deletion[0].variant_class, VariantClass::StructuralDeletion);
        assert_eq!((deletion[0].start, deletion[0].end), (30425601, 30425899));
        assert!(!deletion[0].is_single_breakend);
        assert_eq!(deletion[0].display_allele(), "deletion");
        let breakpoint =
            parse_vcf_line("21\t30425650\tnosvtype\tN\t<BND>\t.\t.\tEND=30425651", true).unwrap();
        assert_eq!(breakpoint[0].variant_class, VariantClass::Translocation);
        assert_eq!(breakpoint[0].display_allele(), "chromosome_breakpoint");
    }

    #[test]
    fn test_parse_vcf_line_symbolic_bnd_without_svlen_single_variant() {
        // <BND> without SVLEN or CHR2: the single base at POS+1, naming no mate.
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
        assert_eq!(v.end, 101); // INFO/END never ends a breakpoint's span
    }

    /// A coordinate is read up to `MAX_COORD`, 2^62 - 1: a POS past it fails the line,
    /// as does a negative POS, `-0` included; an END past it (`i64::MAX` included), a
    /// negative END or text fails the line of every class but a breakpoint, as Ensembl
    /// VEP skips the record (Parser.pm `validate_vf`: `start > end+1`, the size test,
    /// "Invalid start or end coordinate"), while a breakpoint record ignores END for
    /// its span and is the single base at POS + 1 and a record with SVLEN never reads
    /// END (BaseVCF4.pm:302-308); a tandem repeat whose SVLEN is negative ends before
    /// its start and fails the line as Ensembl VEP skips it; and a span from SVLEN
    /// saturates at the ceiling, so a POS just below it with SVLEN `i64::MIN` ends
    /// there whatever the class, and a `<BND>` with CHR2 and such an END2 anchors its
    /// mate there.
    #[test]
    fn coordinates_are_read_up_to_the_ceiling() {
        let past = [MAX_COORD + 1, i64::MAX as u64, u64::MAX];
        for pos in past.iter().map(|n| n.to_string()).chain(["-0".to_string()]) {
            let err = parse_vcf_line(&format!("21\t{pos}\t.\tA\tG\t.\t.\t."), true)
                .unwrap_err()
                .to_string();
            assert!(err.contains("Invalid POS"), "{pos}: {err}");
        }
        for end in past
            .iter()
            .map(|n| n.to_string())
            .chain(["-5".to_string(), "abc".to_string()])
        {
            let err = parse_vcf_line(
                &format!("21\t100\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END={end}"),
                true,
            )
            .unwrap_err()
            .to_string();
            assert!(err.contains("Invalid END"), "{end}: {err}");
            let breakpoint = parse_vcf_line(
                &format!("21\t100\t.\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END={end}"),
                true,
            )
            .unwrap();
            assert_eq!(breakpoint.len(), 1, "{end}");
            assert_eq!(
                (breakpoint[0].start, breakpoint[0].end),
                (101, 101),
                "{end}"
            );
            assert_eq!(
                breakpoint[0].display_allele(),
                "chromosome_breakpoint",
                "{end}"
            );
            let with_svlen = parse_vcf_line(
                &format!("21\t100\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=100;END={end}"),
                true,
            )
            .unwrap();
            assert_eq!(
                (with_svlen[0].start, with_svlen[0].end),
                (101, 200),
                "{end}"
            );
        }
        let err = parse_vcf_line(
            "21\t100\t.\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;SVLEN=-300;RUS=AT;RUC=1",
            true,
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("Invalid SVLEN"), "{err}");
        let kept = parse_vcf_line(
            &format!("21\t100\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END={MAX_COORD}"),
            true,
        )
        .unwrap();
        assert_eq!(kept[0].end, MAX_COORD);
        let svlen = format!("SVLEN={}", i64::MIN);
        let below = MAX_COORD - 1;
        let widest = parse_vcf_line(
            &format!("21\t{below}\t.\tN\t<DUP>\t.\t.\tSVTYPE=DUP;{svlen}"),
            true,
        )
        .unwrap();
        assert_eq!(widest[0].sv_end, Some(MAX_COORD));
        let insertion = parse_vcf_line(
            &format!("21\t{below}\t.\tN\t<INS>\t.\t.\tSVTYPE=INS;{svlen}"),
            true,
        )
        .unwrap();
        assert_eq!(
            (insertion[0].start, insertion[0].end),
            (MAX_COORD, MAX_COORD)
        );
        let breakend = parse_vcf_line(
            &format!(
                "21\t{}\t.\tN\t<BND>\t.\t.\tSVTYPE=BND;{svlen};CHR2=21;END2={MAX_COORD}",
                MAX_COORD - 2
            ),
            true,
        )
        .unwrap();
        assert_eq!(breakend.len(), 2);
        assert_eq!(breakend[1].end, MAX_COORD);
        assert_eq!(breakend[1].mate_pos, Some(MAX_COORD));
    }

    /// A bare `<BND>`, with neither CHR2 nor a span, names no mate: `Parser/VCF.pm`
    /// leaves the ALT as written without CHR2 (the `:491-498` block that rewrites it
    /// is skipped) and prepends REF (`:499`), `_parse_breakends` finds no bracket in
    /// `N/<BND>` (`StructuralVariationFeature.pm:1279-1283`), and the row's Allele
    /// falls to `class_SO_term` (`OutputFactory.pm:2029`), `chromosome_breakpoint`,
    /// with or without INFO/MATEID.
    #[test]
    fn bare_symbolic_bnd_prints_the_class_term() {
        for info in [
            "SVTYPE=BND",
            "SVTYPE=BND;MATEID=bnd_x",
            "SVTYPE=BND;END=101",
        ] {
            let line = format!("21\t100\tbare\tN\t<BND>\t.\t.\t{info}");
            let variants = parse_vcf_line(&line, true).unwrap();
            assert_eq!(variants.len(), 1, "{info}");
            assert_eq!(variants[0].mate_id, None, "{info}");
            assert_eq!(
                variants[0].display_allele(),
                "chromosome_breakpoint",
                "{info}"
            );
        }
        let with_mate = parse_vcf_line(
            "21\t100\tpaired\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=500",
            true,
        )
        .unwrap();
        assert_eq!(with_mate[0].display_allele(), "N.");
        assert_eq!(with_mate[1].display_allele(), "N[21:500[");
    }

    /// ensembl-io `BaseVCF4::get_end`: SVLEN first (`start + |SVLEN| - 1`), then END,
    /// then `get_start` for a symbolic record with neither. A symbolic insertion whose
    /// END is its position (1000 Genomes `<INS:MT>`) or whose SVLEN is 0 (`<INS:ME:ALU>`)
    /// is therefore the pair `start = POS + 1, end = POS`, printed as the two flanks;
    /// one with neither field is the single base at POS + 1.
    #[test]
    fn test_parse_vcf_line_symbolic_insertion_end_at_pos_keeps_the_pair() {
        let line = "11\t7309883\t.\tT\t<INS:MT>\t.\t.\tCIEND=-1,0;CIPOS=0,1;END=7309883;SVTYPE=INS;IMPRECISE";
        let v = &parse_vcf_line(line, true).unwrap()[0];
        assert_eq!(v.variant_class, VariantClass::StructuralInsertion);
        assert_eq!(
            (v.start, v.end, v.sv_end),
            (7309884, 7309883, Some(7309883))
        );
        assert_eq!(v.location(), "11:7309883-7309884");

        let line = "17\t78056052\t.\tG\t<INS:ME:ALU>\t.\t.\tMEINFO=AluY,0,0,-;SVLEN=0;SVTYPE=ALU;TSD=CAGGCGTG";
        let v = &parse_vcf_line(line, true).unwrap()[0];
        assert_eq!(v.variant_class, VariantClass::MobileElementInsertion);
        assert_eq!((v.start, v.end), (78056053, 78056052));
        assert_eq!(v.location(), "17:78056052-78056053");

        let line = "21\t25880572\tsynth_ins_sva_0001\tT\t<INS:ME:SVA>\t.\tPASS\tSVTYPE=INS;MEINFO=SVA_E,64,2808,-";
        let v = &parse_vcf_line(line, true).unwrap()[0];
        assert_eq!((v.start, v.end), (25880573, 25880573));
        assert_eq!(v.location(), "21:25880573");

        // A ranged insertion: SVLEN gives end = POS + SVLEN.
        let line = "21\t1000\t.\tN\t<INS>\t.\t.\tSVTYPE=INS;SVLEN=300";
        let v = &parse_vcf_line(line, true).unwrap()[0];
        assert_eq!((v.start, v.end), (1001, 1300));
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

    /// A `<NON_REF>` record is a reference-confidence block: one structural record
    /// spanning POS+1..END that no consequence path annotates (VEP's own row for it
    /// is `21:101-200 <NON_REF> - - - intergenic_variant`), and whose VCF and JSON
    /// output VEP leaves without consequences (`vep_skip`).
    #[test]
    fn test_parse_vcf_line_non_ref_is_a_reference_block() {
        let line = "21\t100\t.\tA\t<NON_REF>\t.\t.\tEND=200";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        assert_eq!(v.variant_class, VariantClass::ReferenceBlock);
        assert!(v.is_structural);
        assert!(v.vep_skip);
        assert_eq!((v.start, v.end), (101, 200));
        assert_eq!(v.sv_end, Some(200));
        assert_eq!(v.location(), "21:101-200");
        assert_eq!(v.display_allele(), "<NON_REF>");
    }

    /// Only `<NON_REF>` (in any case) is the reference block; the other symbolic
    /// alleles keep their classes.
    #[test]
    fn classify_symbolic_alt_reference_block() {
        assert_eq!(
            classify_symbolic_alt("<NON_REF>", None),
            VariantClass::ReferenceBlock
        );
        assert_eq!(
            classify_symbolic_alt("<non_ref>", None),
            VariantClass::ReferenceBlock
        );
        assert_eq!(
            classify_symbolic_alt("<NON_REF>", Some("NON_REF")),
            VariantClass::ReferenceBlock
        );
        // A supported SVTYPE names the record's type (Parser/VCF.pm replaces an ALT
        // outside the VCF 4.4 forms with INFO/SVTYPE before the SO lookup).
        assert_eq!(
            classify_symbolic_alt("<NON_REF>", Some("CNV")),
            VariantClass::CopyNumberVariation
        );
        assert_eq!(
            classify_symbolic_alt("<CNV>", None),
            VariantClass::CopyNumberVariation
        );
        assert_eq!(
            classify_symbolic_alt("<CPX>", Some("CPX")),
            VariantClass::ComplexStructural
        );
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
                sv_so_term_and_class(&sv_type_abbreviation(alts, svtype)).is_some(),
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

    /// `<CN0>,<CN2>` is one `copy_number_variation` record (Parser.pm `get_SO_term`:
    /// the joined `<CN0>/<CN2>` is neither `<CN0>` nor `<CN2>` alone) written as one
    /// row per ALT, both labelled `copy_number_variation` (record `multi_cn0_cn2` of the
    /// GRCh38-release116 golden corpus: two rows identical except ALLELE_NUM). Each
    /// variant keeps its own ALT bytes and carries the record's term, which is what
    /// gives `<CN0>` the record's copy-number reading rather than a deletion's.
    #[test]
    fn test_parse_vcf_line_multi_allelic_cn_one_variant_per_alt() {
        let line = "21\t14504804\tDUP_gs_CNV\tC\t<CN0>,<CN2>\t.\t.\tSVTYPE=CNV;END=14530597";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2, "one variant per ALT");
        for (idx, v) in variants.iter().enumerate() {
            assert_eq!(v.allele_index, idx);
            assert_eq!(v.variant_class, VariantClass::CopyNumberVariation);
            assert!(v.is_structural);
            assert!(!v.vep_skip);
            assert_eq!(v.display_allele(), "copy_number_variation");
            assert_eq!(v.record_so_term.as_deref(), Some("copy_number_variation"));
            assert_eq!(v.uploaded_allele_string.as_deref(), Some("<CN0>/<CN2>"));
            assert_eq!(v.start, 14504805);
            assert_eq!(v.sv_end, Some(14530597));
        }
        assert_eq!(variants[0].alt_allele(), b"<CN0>");
        assert_eq!(variants[1].alt_allele(), b"<CN2>");
    }

    #[test]
    fn test_parse_vcf_line_multi_allelic_cn_four_alleles() {
        let line = "21\t15244513\t.\tT\t<CN0>,<CN2>,<CN3>,<CN4>\t.\t.\tSVTYPE=CNV;END=15253042";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 4);
        for (v, alt) in variants.iter().zip(["<CN0>", "<CN2>", "<CN3>", "<CN4>"]) {
            assert_eq!(v.variant_class, VariantClass::CopyNumberVariation);
            assert_eq!(v.alt_allele(), alt.as_bytes());
            assert_eq!(v.record_so_term.as_deref(), Some("copy_number_variation"));
        }
    }

    /// A `<DEL>,<DUP>` record is one `copy_number_variation` (Parser.pm `get_SO_term`:
    /// DEL beside DUP) whose two variants keep the ALT bytes the line wrote, `<DEL>`
    /// and `<DUP>`, under the record's term; its JSON `allele_string` is the ALTs
    /// joined by `/`. Ensembl VEP 116.2 writes two `copy_number_variation` rows with
    /// `3_prime_UTR_variant` and `IMPACT=MODIFIER` for `21 30165700 multi_del_dup N
    /// <DEL>,<DUP>` against ENST00000286808 in the GRCh38-release116 golden corpus,
    /// the reading of neither a deletion nor a duplication.
    #[test]
    fn del_beside_dup_keeps_each_alt_under_the_record_term() {
        for (alts, expected_alts) in [
            ("<DEL>,<DUP>", ["<DEL>", "<DUP>"]),
            ("<DUP>,<DEL>", ["<DUP>", "<DEL>"]),
            ("<DEL>,<DUP:TANDEM>", ["<DEL>", "<DUP:TANDEM>"]),
        ] {
            let line = format!("21\t30165700\tmulti_del_dup\tN\t{alts}\t.\t.\tEND=30165900");
            let variants = parse_vcf_line(&line, true).unwrap();
            assert_eq!(variants.len(), 2, "{alts}");
            for (v, alt) in variants.iter().zip(expected_alts) {
                assert_eq!(v.alt_allele(), alt.as_bytes(), "{alts}");
                assert_eq!(v.variant_class, VariantClass::CopyNumberVariation, "{alts}");
                assert_eq!(
                    v.record_so_term.as_deref(),
                    Some("copy_number_variation"),
                    "{alts}"
                );
                assert_eq!(v.display_allele(), "copy_number_variation", "{alts}");
                assert_eq!(
                    v.record_allele_string(),
                    alts.replace(',', "/"),
                    "{alts}: JSON allele_string"
                );
            }
        }
    }

    /// Every multi-allelic symbolic shape of ensembl-vep PR 1998 with its type and
    /// row count under 115.2 (the joined string looked up as written: one allele,
    /// and unsupported unless it was a copy-number list) and under 116.2 (identical
    /// parts collapse, DEL beside DUP is a copy-number variation, one allele per
    /// ALT, every row's Allele the record's `class_SO_term`). The 116.2 column is
    /// the GRCh38-release116 golden corpus's `multi_*` records; `<INS>,<DEL>` and
    /// `<INV>,<DUP>` are skipped by 116.2 and stay annotated here with the
    /// `vep_skip` mark, as 115.2 left them. The mobile-element shapes follow
    /// Parser.pm `get_SO_term` (the `(INS|DEL):ME` match runs once over the joined
    /// string and before the DEL-and-DUP test) and `%SO_TERMS` (`Utils/Config.pm`):
    /// `<INS:ME:ALU>,<INS:ME:LINE1>` is `INS_ALU`, `Alu_insertion` on both rows
    /// (OutputFactory.pm `BaseStructuralVariationOverlapAllele_to_output_hash`).
    #[test]
    fn multi_allelic_symbolic_shapes_take_the_record_class() {
        use VariantClass::*;
        // (ALTs, 115.2 type and alleles, 116.2 classes per variant, supported,
        // 116.2 Allele column per variant)
        let shapes = [
            (
                "<INS>,<INS>",
                "unsupported INS/INS, 1 allele",
                vec![StructuralInsertion; 2],
                true,
                vec!["insertion"; 2],
            ),
            (
                "<DEL>,<DEL>",
                "unsupported DEL/DEL, 1 allele",
                vec![StructuralDeletion; 2],
                true,
                vec!["deletion"; 2],
            ),
            (
                "<DEL>,<DUP>",
                "unsupported DEL/DUP, 1 allele",
                vec![CopyNumberVariation; 2],
                true,
                vec!["copy_number_variation"; 2],
            ),
            (
                "<DUP>,<DEL>",
                "unsupported DUP/DEL, 1 allele",
                vec![CopyNumberVariation; 2],
                true,
                vec!["copy_number_variation"; 2],
            ),
            (
                "<DEL>,<DUP:TANDEM>",
                "tandem_duplication, 1 allele",
                vec![CopyNumberVariation; 2],
                true,
                vec!["copy_number_variation"; 2],
            ),
            (
                "<CN0>,<CN2>",
                "copy_number_variation, 1 allele",
                vec![CopyNumberVariation; 2],
                true,
                vec!["copy_number_variation"; 2],
            ),
            (
                "<INS>,<DEL>",
                "unsupported INS/DEL, 1 allele",
                vec![StructuralInsertion, StructuralDeletion],
                false,
                vec!["insertion", "deletion"],
            ),
            (
                "<INV>,<DUP>",
                "unsupported INV/DUP, 1 allele",
                vec![Inversion, Duplication],
                false,
                vec!["inversion", "duplication"],
            ),
            (
                "G,<INV>",
                "deletion, 1 allele",
                vec![StructuralDeletion; 2],
                true,
                vec!["deletion"; 2],
            ),
            (
                "<INS:ME:ALU>,<INS:ME:ALU>",
                "Alu_insertion, 1 allele",
                vec![MobileElementInsertion; 2],
                true,
                vec!["Alu_insertion"; 2],
            ),
            (
                "<INS:ME:ALU>,<INS:ME:LINE1>",
                "Alu_insertion, 1 allele",
                vec![MobileElementInsertion; 2],
                true,
                vec!["Alu_insertion"; 2],
            ),
            (
                "<INS:ME:ALU>,<INS>",
                "Alu_insertion, 1 allele",
                vec![MobileElementInsertion; 2],
                true,
                vec!["Alu_insertion"; 2],
            ),
            (
                "<INS>,<INS:ME:ALU>",
                "Alu_insertion, 1 allele",
                vec![MobileElementInsertion; 2],
                true,
                vec!["Alu_insertion"; 2],
            ),
            (
                "<DEL:ME:ALU>,<DUP>",
                "Alu_deletion, 1 allele",
                vec![MobileElementDeletion; 2],
                true,
                vec!["Alu_deletion"; 2],
            ),
        ];
        for (alts, release_115, classes, supported, labels) in shapes {
            let line = format!("21\t30165700\t.\tN\t{alts}\t.\t.\tSVTYPE=DEL;END=30165900");
            let variants = parse_vcf_line(&line, true).unwrap();
            let got: Vec<VariantClass> = variants.iter().map(|v| v.variant_class).collect();
            assert_eq!(got, classes, "{alts} (115.2: {release_115})");
            assert!(
                variants.iter().all(|v| v.vep_skip != supported),
                "{alts}: vep_skip must be {}",
                !supported
            );
            let alleles: Vec<usize> = variants.iter().map(|v| v.allele_index).collect();
            assert_eq!(alleles, (0..variants.len()).collect::<Vec<_>>(), "{alts}");
            let shown: Vec<String> = variants.iter().map(|v| v.display_allele()).collect();
            assert_eq!(shown, labels, "{alts}");
            let bytes: Vec<&[u8]> = variants.iter().map(|v| v.alt_allele()).collect();
            let written: Vec<&[u8]> = alts.split(',').map(str::as_bytes).collect();
            assert_eq!(bytes, written, "{alts}: each variant keeps its own ALT");
        }
    }

    /// A record carrying both a mobile-element deletion and a mobile-element
    /// insertion takes the leftmost one: Parser.pm `get_SO_term` matches
    /// `/(INS|DEL):(ME):?([A-Z0-9]+)?/i` once over the joined string (`:685`), so
    /// `<DEL:ME:ALU>,<INS:ME:LINE1>` is `DEL_ALU`, `Alu_deletion` on both rows, and
    /// the reversed order `INS_LINE1`, `LINE1_insertion` on both.
    #[test]
    fn mixed_mobile_element_record_takes_the_leftmost_alt() {
        use VariantClass::*;
        for (alts, abbrev, class, label) in [
            (
                "<DEL:ME:ALU>,<INS:ME:LINE1>",
                "DEL_ALU",
                MobileElementDeletion,
                "Alu_deletion",
            ),
            (
                "<INS:ME:LINE1>,<DEL:ME:ALU>",
                "INS_LINE1",
                MobileElementInsertion,
                "LINE1_insertion",
            ),
        ] {
            assert_eq!(
                sv_type_abbreviation(&alts.replace(',', "/"), Some("DEL")),
                abbrev,
                "{alts}"
            );
            let line = format!("21\t30165700\tme_mix\tN\t{alts}\t.\t.\tEND=30165900;SVTYPE=DEL");
            let variants = parse_vcf_line(&line, true).unwrap();
            assert_eq!(variants.len(), 2, "{alts}");
            for v in &variants {
                assert_eq!(v.variant_class, class, "{alts}");
                assert_eq!(v.record_so_term.as_deref(), Some(label), "{alts}");
                assert_eq!(v.display_allele(), label, "{alts}");
            }
        }
    }

    /// A single mobile-element ALT keeps its subtype label from the same table:
    /// `<INS:ME:ALU>` is `Alu_insertion`, `<INS:ME>` `mobile_element_insertion`,
    /// `<INS:ME:L1>` `LINE1_insertion` (Parser.pm `get_SO_term` maps `L1` to
    /// `LINE1`), `<DEL:ME:SVA>` `SVA_deletion`, and an element outside Parser.pm's
    /// list (`<INS:ME:LINE>`) the bare `mobile_element_insertion`. Ensembl VEP
    /// 116.2 prints `Alu_insertion` for `21 10524521 gnomAD-SV_v3_INS_chr21_1d284ae3
    /// N <INS:ME:ALU>` and `mobile_element_insertion` for `21 31659650
    /// synth_ins_me_0117 A <INS:ME>` in the GRCh38 golden corpus, and
    /// `mobile_element_insertion` for `21 27284085 synth_ins_line_0011 G
    /// <INS:ME:LINE>` and `mobile_element_deletion` for `21 33041214
    /// synth_del_me_0080 G <DEL:ME>` in the GRCh37 golden corpus.
    #[test]
    fn single_mobile_element_alts_keep_their_subtype_labels() {
        for (alt, class, label) in [
            (
                "<INS:ME:ALU>",
                VariantClass::MobileElementInsertion,
                "Alu_insertion",
            ),
            (
                "<INS:ME>",
                VariantClass::MobileElementInsertion,
                "mobile_element_insertion",
            ),
            (
                "<INS:ME:L1>",
                VariantClass::MobileElementInsertion,
                "LINE1_insertion",
            ),
            (
                "<INS:ME:LINE>",
                VariantClass::MobileElementInsertion,
                "mobile_element_insertion",
            ),
            (
                "<DEL:ME:SVA>",
                VariantClass::MobileElementDeletion,
                "SVA_deletion",
            ),
            (
                "<DEL:ME>",
                VariantClass::MobileElementDeletion,
                "mobile_element_deletion",
            ),
        ] {
            let line = format!("21\t10524521\t.\tN\t{alt}\t.\t.\tSVTYPE=INS;END=10524572");
            let variants = parse_vcf_line(&line, true).unwrap();
            assert_eq!(variants.len(), 1, "{alt}");
            assert_eq!(variants[0].variant_class, class, "{alt}");
            assert_eq!(variants[0].display_allele(), label, "{alt}");
            assert_eq!(variants[0].record_so_term.as_deref(), Some(label), "{alt}");
        }
    }

    /// `<CN0>,<CN0>` collapses to `<CN0>`, a deletion, under 116.2's normalisation
    /// (Parser.pm `get_SO_term`: identical parts reduce to one).
    #[test]
    fn repeated_cn0_collapses_to_a_deletion() {
        let variants =
            parse_vcf_line("21\t100\t.\tA\t<CN0>,<CN0>\t.\t.\tSVTYPE=CNV;END=800", true).unwrap();
        assert_eq!(variants.len(), 2);
        for v in &variants {
            assert_eq!(v.variant_class, VariantClass::StructuralDeletion);
            assert_eq!(v.alt_allele(), b"<CN0>");
            assert!(!v.vep_skip);
        }
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

    /// `_expand_tandem_repeat_allele_string` (`Parser/VCF.pm`) builds `RUS x RUC` when
    /// the record carries RUC and `RUS x (RB / length(RUS))` otherwise (`:414`), so
    /// `RB=42;RUC=7;RUS=CAG` is `CAG x 7`, 21 bases, and `RB=307;RUS=AT` is
    /// `AT x 153.5`, which Perl's `x` truncates to 306 bases (`:419`); a `.` unit is
    /// the one base `N` and RUL is never read (`:413`), so `RUC=7;RUS=.;RUL=4` is
    /// `NNNNNNN`, 7 bases; an empty or `0` RUC is no RUC (`:399`, Perl's `||`); without
    /// RUS nothing expands (`:385`).
    #[test]
    fn tandem_repeat_alt_bases_reads_ruc_before_rb_and_never_rul() {
        assert_eq!(alt_bases("SVTYPE=CNV;RB=42;RUC=7;RUS=CAG", 0), Some(21));
        assert_eq!(alt_bases("SVTYPE=CNV;RUC=7;RUS=CAG", 0), Some(21));
        assert_eq!(alt_bases("SVTYPE=CNV;RUC=7.9;RUS=CAG", 0), Some(21));
        assert_eq!(alt_bases("SVTYPE=CNV;RB=42;RUS=CAG", 0), Some(42));
        assert_eq!(alt_bases("SVTYPE=CNV;RB=307;RUS=AT", 0), Some(306));
        assert_eq!(
            alt_bases("SVTYPE=CNV;RUC=153.9;RB=307;RUS=AT", 0),
            Some(306)
        );
        assert_eq!(alt_bases("SVTYPE=CNV;RUC=7;RUS=.;RUL=4", 0), Some(7));
        assert_eq!(alt_bases("SVTYPE=CNV;RUC=;RB=42;RUS=CAG", 0), Some(42));
        assert_eq!(alt_bases("SVTYPE=CNV;RUC=0;RB=42;RUS=CAG", 0), Some(42));
        assert_eq!(alt_bases("SVTYPE=CNV;RB=42,60;RUS=CAG,AT", 0), Some(42));
        assert_eq!(alt_bases("SVTYPE=CNV;RUC=7;RUL=4", 0), None);
        assert_eq!(alt_bases("SVTYPE=CNV;RUC=7", 0), None);
        assert_eq!(alt_bases("SVTYPE=CNV;SVLEN=300", 0), None);
    }

    /// `_expand_tandem_repeat_allele_string` (Parser/VCF.pm) builds one alternate
    /// allele per RN entry (one per `<CNV:TR>` allele when RN is absent) by
    /// concatenating that many RUS x RUC repeat sequences in order, so the allele's
    /// length is their sum: `RN=2;RUS=CT,GA;RUC=2,1` expands to `CTCTGA`, 6 bases,
    /// as does `RB=4,2`; with two `<CNV:TR>` alleles and `RUS=CT,CT;RUC=2,8` the
    /// first allele is `CTCT` (4 bases) and the second `CT` x 8 (16 bases); a
    /// second allele whose entries the lists lack has no length.
    #[test]
    fn tandem_repeat_alt_bases_sums_an_allele_s_repeat_sequences() {
        assert_eq!(
            alt_bases("SVTYPE=CNV;SVLEN=10;RN=2;RUS=CT,GA;RUC=2,1", 0),
            Some(6)
        );
        assert_eq!(
            alt_bases("SVTYPE=CNV;SVLEN=10;RN=2;RUS=CT,GA;RB=4,2", 0),
            Some(6)
        );
        let two_alleles = "SVTYPE=CNV;SVLEN=10;RUS=CT,CT;RUC=2,8";
        assert_eq!(alt_bases(two_alleles, 0), Some(4));
        assert_eq!(alt_bases(two_alleles, 1), Some(16));
        assert_eq!(
            alt_bases("SVTYPE=CNV;SVLEN=10;RN=1,2;RUS=CT,GA,TTA;RUC=2,1,2", 1),
            Some(8)
        );
        assert_eq!(alt_bases("SVTYPE=CNV;RUS=CT;RUC=2", 1), None);
        assert_eq!(alt_bases("SVTYPE=CNV;RN=.;RUS=CT;RUC=2", 0), None);

        let line = "21\t25002099\tmulti_tr\tA\t<CNV:TR>,<CNV:TR>\t.\t.\tSVTYPE=CNV;END=25002109;SVLEN=10;RUS=CT,CT;RUC=2,8";
        let variants = parse_vcf_line(line, true).unwrap();
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].tr_alt_bases, Some(4));
        assert_eq!(variants[1].tr_alt_bases, Some(16));
    }

    /// An RN count that pushes the running sum past the machine word states no
    /// length for any allele: the first allele's `usize::MAX` sequences exhaust the
    /// two RUS entries, the second allele's one sequence takes the sum past
    /// `usize::MAX`, and a third allele's offset lies past it. An RN of `usize::MAX`
    /// on its own sums without overflow: the oversize check reads only the one
    /// sequence RUS lists, and the allele's second sequence finds no unit, so it
    /// states no length. Under a limit of `u64::MAX`, which admits a `1e19`-base
    /// sequence, one such sequence is its length and two sum past `u64::MAX`.
    #[test]
    fn tandem_repeat_alt_bases_refuses_an_rn_sum_past_the_machine_word() {
        let info = "SVTYPE=CNV;SVLEN=10;RN=18446744073709551615,1;RUS=AT,AT;RUC=100,100";
        assert_eq!(alt_bases(info, 0), None);
        assert_eq!(alt_bases(info, 1), None);
        assert_eq!(alt_bases(info, 2), None);
        let line = format!("21\t25002099\tmulti_tr\tA\t<CNV:TR>,<CNV:TR>\t.\t.\t{info}");
        let variants = parse_vcf_line(&line, true).unwrap();
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].tr_alt_bases, None);
        assert_eq!(variants[1].tr_alt_bases, None);
        assert_eq!(
            alt_bases("SVTYPE=CNV;RN=18446744073709551615;RUS=AT;RUC=100", 0),
            None
        );
        let unlimited = Some(u64::MAX);
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RUS=A;RUC=1e19", 0, 1, unlimited),
            Some(10_000_000_000_000_000_000)
        );
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RN=2;RUS=A,A;RUC=1e19,1e19", 0, 1, unlimited),
            None
        );
    }

    /// `_expand_tandem_repeat_allele_string` keeps the symbolic reading when any repeat
    /// sequence is longer than `max_sv_size` (`Parser/VCF.pm:417-418`, `:428`), read
    /// as `$self->param('max_sv_size') || 5000` (`:403`), so `-1` passes `||` unchanged
    /// and nothing expands under it: `RUS=AT;RUC=100` is 200 bases under the
    /// 10,000,000 default and under a cap of 200, symbolic under 199 and under `-1`;
    /// `RUS=CAG;RUC=4000000`, 12,000,000 bases, is symbolic under the default; one
    /// oversized allele keeps every allele of its record symbolic; a cap of `0` is
    /// Perl's 5000.
    #[test]
    fn tandem_repeat_expansion_stops_at_max_sv_size() {
        let info = "SVTYPE=CNV;RUS=AT;RUC=100;RB=300";
        assert_eq!(
            tandem_repeat_alt_bases(info, 0, 1, Some(10_000_000)),
            Some(200)
        );
        assert_eq!(tandem_repeat_alt_bases(info, 0, 1, Some(200)), Some(200));
        assert_eq!(tandem_repeat_alt_bases(info, 0, 1, Some(199)), None);
        assert_eq!(tandem_repeat_alt_bases(info, 0, 1, None), None);
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RUS=AT;RUC=4000000", 0, 1, Some(10_000_000)),
            Some(8_000_000)
        );
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RUS=CAG;RUC=4000000", 0, 1, Some(10_000_000)),
            None
        );
        let two = "SVTYPE=CNV;RUS=AT,CAG;RUC=10,4000000";
        assert_eq!(tandem_repeat_alt_bases(two, 0, 2, Some(10_000_000)), None);
        assert_eq!(
            tandem_repeat_alt_bases(two, 0, 2, Some(20_000_000)),
            Some(20)
        );
        assert_eq!(
            tandem_repeat_alt_bases(two, 1, 2, Some(20_000_000)),
            Some(12_000_000)
        );
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RUS=AT;RUC=2600", 0, 1, Some(0)),
            None
        );
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RUS=AT;RUC=2500", 0, 1, Some(0)),
            Some(5000)
        );

        let line = "21\t27253937\ttr\tC\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=27254237;SVLEN=300;RUS=AT;RUC=100;RB=300";
        assert_eq!(
            parse_vcf_line(line, true).unwrap()[0].tr_alt_bases,
            Some(200)
        );
        for limit in [None, Some(100)] {
            let variants = parse_vcf_line_with_max_sv_size(line, true, limit).unwrap();
            assert_eq!(variants[0].tr_alt_bases, None, "{limit:?}");
        }
    }

    /// `_expand_tandem_repeat_allele_string` builds one allele per RN entry (`:408`,
    /// `for my $records (@RN)`), so a `<CNV:TR>` ALT past a short RN list has no
    /// materialised allele and states no length, and the oversize scan covers only the
    /// sequences the list materialises: `RN=1;RUS=AT,CAGA;RUC=10,4000000` on two ALTs
    /// is 20 bases for the first and nothing for the second, whose 16,000,000-base
    /// sequence is never built; `RN=1;RUS=AT,AT;RUC=10,10` likewise leaves the second
    /// ALT without a length, where `RN=1,1` gives both 20.
    #[test]
    fn an_alt_past_the_rn_list_states_no_length() {
        let limit = Some(10_000_000);
        let record = |info: &str| {
            let line = format!("21\t30425600\tx\tN\t<CNV:TR>,<CNV:TR>\t.\t.\tEND=30425899;{info}");
            let variants = parse_vcf_line_with_max_sv_size(&line, true, limit).unwrap();
            (variants[0].tr_alt_bases, variants[1].tr_alt_bases)
        };
        assert_eq!(
            record("SVTYPE=CNV;RN=1;RUS=AT,CAGA;RUC=10,4000000"),
            (Some(20), None)
        );
        assert_eq!(
            record("SVTYPE=CNV;RN=1;RUS=AT,AT;RUC=10,10"),
            (Some(20), None)
        );
        assert_eq!(
            record("SVTYPE=CNV;RN=1,1;RUS=AT,AT;RUC=10,10"),
            (Some(20), Some(20))
        );
        assert_eq!(
            record("SVTYPE=CNV;RUS=AT,AT;RUC=10,10"),
            (Some(20), Some(20))
        );
        assert_eq!(
            tandem_repeat_alt_bases("SVTYPE=CNV;RN=1;RUS=AT,CAGA;RUC=10,4000000", 1, 2, limit),
            None
        );
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
    /// shared leading base, then `post_process_vfs` (`VEP/Parser.pm:860-882`) calls
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
    /// (`VEP/Parser.pm:978-984`), so a multi-allelic record never reaches
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
            "a 2-ALT record returns early at Parser.pm:978 and keeps its anchor \
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
        // BND with SVLEN=-1 and CHR2+END2: CHR2 with END2 names the mate whatever
        // the span, so the record yields its own breakend and the mate even though
        // SVLEN leaves end == start.
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
        // CHR2 with END2 names the mate -> 2 variants
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
