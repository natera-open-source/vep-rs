// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! VCF line parsing and allele processing.
//!
//! Parses a VCF data line into `InputVariant`s: the line is routed through
//! `vep_core::sv_routing`, a sequence record is split per ALT with its alleles
//! trimmed and minimised here, and a structural record's alleles are built by
//! `vep_core::sv_record::structural_variants`.
//!
//! Perl citations name modules of ensembl-vep release/116 (`Bio/EnsEMBL/VEP/...`);
//! `Sequence.pm`, `StructuralVariationFeature.pm` and `StructuralVariationOverlap.pm`
//! are in ensembl-variation release/116 and `BaseVCF4.pm` is in ensembl-io release/116.

use std::sync::Arc;

use anyhow::{bail, Context};

use vep_core::allele::trim_alleles;
use vep_core::coordinate::parse_pos;
use vep_core::sv_record::{structural_variants, DEFAULT_MAX_SV_SIZE};
use vep_core::sv_routing::{perl_true, routes_to_structural, variant_alts, written_alts};
use vep_core::variant::{InputVariant, VariantClass};
use vep_core::vcf_info::info_value;

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
    parse_vcf_line_with_max_sv_size(line, capture_raw_input, Some(DEFAULT_MAX_SV_SIZE))
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
    let pos = parse_pos(fields[1]).with_context(|| format!("Invalid POS: '{}'", fields[1]))?;
    let id = fields[2];
    let ref_allele = fields[3];
    let alt_field = fields[4];

    if alt_field == "." {
        bail!("No alternate allele");
    }

    let info_field = fields[7];
    let svtype = info_value(info_field, "SVTYPE").filter(|v| perl_true(v));

    // The ALT column as Perl's list-context `split` reads it, a gVCF variant site read
    // as its real alleles (`variant_alts`), each ALT keeping its written position.
    let written = written_alts(alt_field);
    let mut alts = written.clone();
    variant_alts(&mut alts, svtype);
    // `*` and `.` are written as no allele here; both stay in the record's allele
    // string, which Ensembl VEP builds from every ALT (Parser/VCF.pm:469). Ensembl VEP
    // writes no row for `*` (neither the indel nor the SV path fires) and reads a `.`
    // beside another ALT as an allele of its own. `<*>` is not filtered: Perl treats
    // it as a non-SV allele at POS unless INFO/SVTYPE routes the record. `<CPX>` is
    // not filtered: absent from `%SO_TERMS`, it is parsed with the `vep_skip` mark
    // (`sv_so_term_and_class`).
    let valid_alts: Vec<(usize, &str)> = alts
        .iter()
        .copied()
        .filter(|(_, alt)| *alt != "*" && *alt != ".")
        .collect();

    let raw_input: Option<Arc<str>> = capture_raw_input.then(|| Arc::from(line));
    let explicit_id = (id != "." && !id.is_empty()).then(|| id.to_string());
    let line_name: Arc<str> = input_line_name(&chr, pos, ref_allele, alt_field).into();

    if valid_alts.is_empty() {
        if written.is_empty() || written.iter().any(|(_, alt)| *alt != "*") {
            bail!("No alternate allele");
        }
        let mut variant =
            InputVariant::without_allele(chr, pos, ref_allele, written.iter().map(|(_, alt)| *alt));
        // The record's raw string, as on the sequence path: a VCF record's JSON `id`
        // keeps the ID column's `.` (`variation_name`), told by this field being set.
        variant.uploaded_allele_string = Some(variant.allele_string.as_str().into());
        variant.raw_input = raw_input;
        variant.id = explicit_id;
        variant.input_line_name = Some(line_name);
        return Ok(vec![variant]);
    }

    let mut variants = Vec::with_capacity(valid_alts.len());

    if !routes_to_structural(ref_allele, &alts, svtype) {
        let alt_strings: Vec<&str> = alts.iter().map(|(_, alt)| *alt).collect();
        let (mut trimmed_ref, mut trimmed_alts, prefix_trimmed) =
            trim_multi_alleles(ref_allele, &alt_strings);
        // `validate_vf` upper-cases the allele string after the anchor chop and before
        // `post_process_vfs` minimises it (Parser.pm:554, `tr/[a-z]/[A-Z]/`).
        trimmed_ref.make_ascii_uppercase();
        trimmed_alts
            .iter_mut()
            .for_each(|alt| alt.make_ascii_uppercase());
        let (vep_start, vep_end) =
            vep_coords_from_trimmed(pos, prefix_trimmed, trimmed_ref.as_bytes());

        // The allele string VEP names an ID-less record by: the raw REF and every
        // ALT joined by `/`, `*` and `.` at their places. VEP keeps the raw string as
        // `nontrimmed_allele_string` and promotes it to `original_allele_string` for
        // any record whose alleles differ in length; a record without an indel allele
        // is never trimmed, so its raw and working strings coincide (Parser/VCF.pm,
        // Parser.pm `post_process_vfs` and `minimise_alleles`, OutputFactory.pm).
        let uploaded_alleles = allele_string(ref_allele, alts.iter().map(|(_, alt)| *alt));
        // Perl's working string of a record with several ALTs, `*` and `.` counted
        // (`scalar @$alts > 1`, Parser/VCF.pm:301-322): every allele behind the chopped
        // anchor, a string `minimise_alleles` leaves as it is.
        let chopped = (alts.len() > 1)
            .then(|| allele_string(&trimmed_ref, trimmed_alts.iter().map(String::as_str)));

        for ((allele_idx, alt), vep_alt) in alts.iter().zip(&trimmed_alts) {
            if *alt == "*" || *alt == "." {
                continue;
            }
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
            // multi-allelic string (`.+/.+/.+`) unchanged, a `*` part counted
            // (Parser.pm:978-984): without `--minimal` a multi-allelic record keeps
            // its anchor-chopped alleles, and every position, codon and distance is
            // computed on that representation.
            if alts.len() == 1
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
            variant.record_allele_string_multi = chopped.clone();
            variant.raw_input = raw_input.clone();
            variant.id = explicit_id.clone();
            variants.push(variant);
        }
    } else {
        variants = structural_variants(
            &chr,
            pos,
            ref_allele,
            &written,
            info_field,
            max_sv_size,
            explicit_id.as_deref(),
        )?;
        if variants
            .iter()
            .any(|v| v.variant_class == VariantClass::ReferenceBlock)
        {
            tracing::warn!(
                "{}:{} {}: NON_REF is not a supported structural variant type; \
                 a reference-confidence block is written as one intergenic row",
                chr,
                pos,
                alt_field
            );
        }
        for variant in &mut variants {
            variant.raw_input = raw_input.clone();
            variant.id = explicit_id.clone();
        }
    }

    for variant in &mut variants {
        variant.input_line_name = Some(line_name.clone());
    }

    Ok(variants)
}

/// `REF/ALT1/ALT2`: the reference and every ALT joined by `/`, one buffer for every
/// allele of the record.
fn allele_string<'a>(reference: &str, alts: impl IntoIterator<Item = &'a str>) -> Arc<str> {
    let mut s = reference.to_string();
    for alt in alts {
        s.push('/');
        s.push_str(alt);
    }
    s.into()
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

/// The anchor chop of Parser/VCF.pm:294-335 over every written ALT, `*` and `.` included:
/// the record is an indel when any ALT's length differs from REF's (a `*` of one base
/// against a two-base REF counts, so `AG AC,*` is chopped where `AG AC` is not, :296-298),
/// and the first base is chopped when REF and every ALT other than `*` share it (:306-308
/// with several ALTs, :328 with one), `*` kept as written and an emptied allele reading
/// `-` (:315-316). Same-length alleles are never trimmed; the secondary minimisation of a
/// bi-allelic indel is `parse_vcf_line`'s (Parser.pm:860-881). Returns the trimmed REF, the
/// trimmed ALTs in their written order and the number of bases chopped.
pub(crate) fn trim_multi_alleles(
    ref_allele: &str,
    alt_alleles: &[&str],
) -> (String, Vec<String>, usize) {
    let ref_bytes = ref_allele.as_bytes();
    let is_indel = alt_alleles.iter().any(|a| a.len() != ref_bytes.len());
    let shares_first_base = |allele: &&str| {
        allele
            .as_bytes()
            .first()
            .is_some_and(|b| b.eq_ignore_ascii_case(&ref_bytes[0]))
    };
    // With one ALT its own first base decides (`<*>` against `N` is not chopped); with
    // several, the alleles carrying a `*` are left out of the comparison.
    let anchored = is_indel
        && !ref_bytes.is_empty()
        && match alt_alleles {
            [single] => shares_first_base(single),
            _ => alt_alleles
                .iter()
                .filter(|a| !a.contains('*'))
                .all(shares_first_base),
        };
    let prefix_len = usize::from(anchored);
    let chop = |allele: &str| match &allele.as_bytes()[prefix_len..] {
        [] => "-".to_string(),
        rest => String::from_utf8_lossy(rest).into_owned(),
    };
    let trimmed_alts = alt_alleles
        .iter()
        .map(|alt| {
            if alt.contains('*') {
                alt.to_string()
            } else {
                chop(alt)
            }
        })
        .collect();
    (chop(ref_allele), trimmed_alts, prefix_len)
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
    use vep_core::coordinate::MAX_COORD;
    use vep_core::sv_record::{parse_ci_field, tandem_repeat_alt_bases};
    use vep_core::sv_routing::{classify_symbolic_alt, sv_so_term_and_class, sv_type_abbreviation};
    use vep_core::variant::VariantClass;

    /// `tandem_repeat_alt_bases` for a record of `ordinal + 1` tandem-repeat alleles
    /// under Ensembl VEP's default `--max_sv_size`.
    fn alt_bases(info: &str, ordinal: usize) -> Option<u64> {
        tandem_repeat_alt_bases(info, ordinal, ordinal + 1, Some(DEFAULT_MAX_SV_SIZE))
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

    /// A `*`-only record is named like every other VCF record: Parser/VCF.pm:347 gives
    /// its feature the ID column as written (`.` for a dot), the `id` of its JSON object
    /// (OutputFactory/JSON.pm:236), and OutputFactory.pm's line branch names it
    /// `CHROM_POS_REF/ALT`; the alleles that build that name are the record's own.
    #[test]
    fn dot_id_names_a_star_only_record_by_its_line() {
        let variants = parse_vcf_line("21\t30425650\t.\tN\t*\t.\t.\t.", false).unwrap();
        assert!(variants[0].no_allele);
        assert_eq!(variants[0].variation_name(), ".");
        assert_eq!(variants[0].uploaded_variation(), "21_30425650_N/*");
        assert_eq!(variants[0].uploaded_allele_string.as_deref(), Some("N/*"));
        let named = parse_vcf_line("21\t30425650\tstar\tN\t*,*\t.\t.\t.", false).unwrap();
        assert_eq!(named[0].variation_name(), "star");
        assert_eq!(named[0].uploaded_allele_string.as_deref(), Some("N/*/*"));
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

    /// A record whose REF and ALTs are a run of `ACGT` is never a structural variant,
    /// whatever INFO/SVTYPE says (Parser/VCF.pm:238, `join("", $ref, @$alts) !~
    /// /^[ACGT]+$/`).
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
    /// (`StructuralVariationFeature.pm:1276`). A bare `CHR2` names none (undef to
    /// `defined $info->{CHR2}`, `:491`), and `CHR2=` is defined and empty, so it builds
    /// `N[:30425800[`, which the bracket regex refuses: Ensembl VEP 116.2 writes
    /// `chromosome_breakpoint` for `CHR2=;END2=30425800`.
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

    /// The mate INFO/CHR2 names is written `N[CHR2:pos[` whatever the reference base
    /// (Parser/VCF.pm:495-496, `"N"` for a lone `<BND>`), and the own breakend takes the
    /// reference base (`_parse_breakends`, the `A` part of `A/N[21:30425700[`): Ensembl
    /// VEP 116.2 writes `A.` and `N[21:30425700[` with ALLELE_NUM 0 and 1 for REF `A`,
    /// allele string `A/N[21:30425700[`, over the SVLEN span when there is one. A
    /// `<BND:X>` keeps its part (`N/<BND:X>/N[21:30425700[`) and `<BND>,<BND>` both of
    /// its (`N/<BND>/<BND>/N[21:30425700[`), each with one own breakend and one mate.
    #[test]
    fn a_synthesised_mate_is_a_literal_n_and_the_own_breakend_the_reference_base() {
        let variants = parse_vcf_line(
            "21\t30425650\tx\tA\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700",
            true,
        )
        .unwrap();
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].display_allele(), "A.");
        assert_eq!(variants[1].display_allele(), "N[21:30425700[");
        assert_eq!(variants[1].mate_pos, Some(30425700));
        assert_eq!(variants[1].mate_chr.as_deref(), Some("21"));
        for v in &variants {
            assert_eq!(v.record_allele_string(), "A/N[21:30425700[");
            assert_eq!((v.start, v.end), (30425651, 30425651));
        }
        let spanned = parse_vcf_line(
            "21\t30425650\tx\tA\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700;SVLEN=100",
            true,
        )
        .unwrap();
        assert_eq!(spanned.len(), 2);
        assert!(spanned
            .iter()
            .all(|v| (v.start, v.end) == (30425651, 30425750)));
        for (alt, record) in [
            ("<BND:X>", "N/<BND:X>/N[21:30425700["),
            ("<BND>,<BND>", "N/<BND>/<BND>/N[21:30425700["),
        ] {
            let variants = parse_vcf_line(
                &format!("21\t30425650\tx\tN\t{alt}\t.\t.\tSVTYPE=BND;CHR2=21;END2=30425700"),
                true,
            )
            .unwrap();
            assert_eq!(variants.len(), 2, "{alt}");
            assert_eq!(variants[0].display_allele(), "N.", "{alt}");
            assert_eq!(variants[1].display_allele(), "N[21:30425700[", "{alt}");
            assert_eq!(variants[0].record_allele_string(), record, "{alt}");
        }
    }

    /// The bracket text keeps INFO's values as written and the record's string carries
    /// the text whether or not a mate is parsed from it: Ensembl VEP 116.2 writes the
    /// Allele `N[21:00030425700[` for `END2=00030425700` (the mate at 30425700), the
    /// allele string `N/N[21:-7[`, `N/N[21:+30425700[` and `N/N[21:abc[` with
    /// `chromosome_breakpoint` rows for an END2 the regex refuses, `N/N[:30425800[` for
    /// `CHR2=` (defined and empty, Perl's split keeps the field), and `N/<BND>` for a bare
    /// `CHR2` flag or an END2 written empty with no END.
    #[test]
    fn a_mate_text_follows_the_info_values_as_written() {
        let padded = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=00030425700",
            true,
        )
        .unwrap();
        assert_eq!(padded.len(), 2);
        assert_eq!(padded[1].display_allele(), "N[21:00030425700[");
        assert_eq!(padded[1].mate_pos, Some(30425700));
        assert_eq!(padded[0].record_allele_string(), "N/N[21:00030425700[");
        for (info, record) in [
            ("SVTYPE=BND;CHR2=21;END2=-7", "N/N[21:-7["),
            ("SVTYPE=BND;CHR2=21;END2=+30425700", "N/N[21:+30425700["),
            ("SVTYPE=BND;CHR2=21;END2=abc", "N/N[21:abc["),
            ("SVTYPE=BND;CHR2=;END2=30425800", "N/N[:30425800["),
            ("SVTYPE=BND;CHR2;END2=30425800", "N/<BND>"),
            ("SVTYPE=BND;CHR2=21;END2=", "N/<BND>"),
        ] {
            let variants =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t<BND>\t.\t.\t{info}"), true).unwrap();
            assert_eq!(variants.len(), 1, "{info}");
            assert_eq!(
                variants[0].display_allele(),
                "chromosome_breakpoint",
                "{info}"
            );
            assert_eq!(variants[0].record_allele_string(), record, "{info}");
            assert!(
                variants[0].is_single_breakend && variants[0].mate_id.is_none(),
                "{info}"
            );
        }
    }

    /// A dot-bounded ALT beside a bracket ALT lends the own breakend its base and gets no
    /// allele of its own (`_parse_breakends`, StructuralVariationFeature.pm:1280-1283
    /// against :1322): Ensembl VEP 116.2 writes `N.` and `N[21:30425700[` for
    /// `.N,N[21:30425700[` with REF `N` or `A`, and `.N` and `]21:30425700]N` for
    /// `N.,]21:30425700]N` (left-placed mate), whose string is `N/N./]21:30425700]N`.
    #[test]
    fn a_dot_alt_beside_a_bracket_alt_lends_the_own_breakend_its_base() {
        for reference in ["N", "A"] {
            let variants = parse_vcf_line(
                &format!("21\t30425650\tx\t{reference}\t.N,N[21:30425700[\t.\t.\tSVTYPE=BND"),
                true,
            )
            .unwrap();
            assert_eq!(variants.len(), 2, "{reference}");
            assert_eq!(variants[0].display_allele(), "N.", "{reference}");
            assert_eq!(
                variants[1].display_allele(),
                "N[21:30425700[",
                "{reference}"
            );
            assert_eq!(
                variants[0].record_allele_string(),
                ".N/N[21:30425700[",
                "{reference}"
            );
        }
        let left = parse_vcf_line(
            "21\t30425650\tx\tN\tN.,]21:30425700]N\t.\t.\tSVTYPE=BND",
            true,
        )
        .unwrap();
        assert_eq!(left.len(), 2);
        assert_eq!(left[0].display_allele(), ".N");
        assert_eq!(left[1].display_allele(), "]21:30425700]N");
        assert_eq!(left[0].record_allele_string(), "N/N./]21:30425700]N");
        assert_eq!(left[1].mate_pos, Some(30425700));
    }

    /// SVLEN is read as Perl's `abs()` numifies it (BaseVCF4.pm:304-305): Ensembl VEP
    /// 116.2 writes the pair `21:30425650-30425651` for `<DUP>`, `<INV>`, `<CNV>`, `<INS>`,
    /// `<BND>` and `<INS:ME:ALU>` with `SVLEN=0`, `SVLEN=abc` or `SVLEN=` (an END beside
    /// them unread), and `21:30425651-30425950` for `SVLEN=300abc`, `SVLEN= 300` and
    /// `SVLEN=3e2`. A `<DEL>` with the pair is a record Ensembl VEP skips ("deletion looks
    /// incomplete"), parsed here as the pair.
    #[test]
    fn svlen_is_read_as_perl_numifies_it() {
        for alt in [
            "<DUP>",
            "<INV>",
            "<CNV>",
            "<INS>",
            "<BND>",
            "<INS:ME:ALU>",
            "<DEL>",
        ] {
            for svlen in ["0", "abc", ""] {
                let variants = parse_vcf_line(
                    &format!("21\t30425650\tx\tN\t{alt}\t.\t.\tSVLEN={svlen};END=30425899"),
                    true,
                )
                .unwrap();
                assert_eq!(variants.len(), 1, "{alt} {svlen:?}");
                assert_eq!(
                    (variants[0].start, variants[0].end),
                    (30425651, 30425650),
                    "{alt} {svlen:?}"
                );
                assert_eq!(variants[0].sv_len, Some(0), "{alt} {svlen:?}");
            }
        }
        for svlen in ["300abc", " 300", "3e2", "+300", "-300"] {
            let variants = parse_vcf_line(
                &format!("21\t30425650\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;SVLEN={svlen}"),
                true,
            )
            .unwrap();
            assert_eq!(
                (variants[0].start, variants[0].end),
                (30425651, 30425950),
                "{svlen:?}"
            );
        }
    }

    /// `validate_vf` (Parser.pm:483-488) skips a record whose end lies more than one base
    /// before its start, and Ensembl VEP 116.2 does so for `END=0` and `END=30425000` at
    /// POS 30425650 for every class, a `<BND>` included, whose deleted END is re-read only
    /// when Perl reads it as true (Parser/VCF.pm:489), so `END=0` keeps the 0 and `END=`
    /// the empty string ("Invalid start or end coordinate"). `END=30425650` is the pair
    /// and `END=30425651` the single base at POS + 1.
    #[test]
    fn an_end_before_pos_is_refused_and_an_end_at_pos_is_the_pair() {
        for (alt, info) in [
            ("<DUP>", "SVTYPE=DUP;END=0"),
            ("<DUP>", "SVTYPE=DUP;END=30425000"),
            ("<DEL>", "SVTYPE=DEL;END=0"),
            ("<BND>", "SVTYPE=BND;CHR2=21;END=0"),
            ("<BND>", "SVTYPE=BND;END=0"),
        ] {
            let err = parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true)
                .unwrap_err()
                .to_string();
            assert!(err.contains("start > end+1"), "{alt} {info}: {err}");
        }
        for (alt, info) in [("<BND>", "SVTYPE=BND;END="), ("<DUP>", "SVTYPE=DUP;END=")] {
            let err = parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true)
                .unwrap_err()
                .to_string();
            assert!(err.contains("Invalid END: ''"), "{alt} {info}: {err}");
        }
        let pair = parse_vcf_line(
            "21\t30425650\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=30425650",
            true,
        )
        .unwrap();
        assert_eq!((pair[0].start, pair[0].end), (30425651, 30425650));
        let point = parse_vcf_line(
            "21\t30425650\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=30425651",
            true,
        )
        .unwrap();
        assert_eq!((point[0].start, point[0].end), (30425651, 30425651));
        let kept_bnd = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;SVLEN=100;END=0",
            true,
        )
        .unwrap();
        assert_eq!((kept_bnd[0].start, kept_bnd[0].end), (30425651, 30425750));
    }

    /// A repeated INFO key reads its last value, as the hash Ensembl VEP's parser builds
    /// does (BaseVCF4.pm:581-583): Ensembl VEP 116.2 writes `21:30425651-30425950` for
    /// `SVLEN=-100;SVLEN=-300`, `21:30425651-30425900` for `END=30425700;END=30425900`,
    /// the mate `N[21:30425800[` for `END2=30425700;END2=30425800` and for
    /// `CHR2=2;CHR2=21`, and `duplication` for `N.` with `SVTYPE=DEL;SVTYPE=DUP`. A bare
    /// key written last is undef, so `CHR2=21;CHR2` names no mate.
    #[test]
    fn a_repeated_info_key_reads_its_last_value() {
        let svlen = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-100;SVLEN=-300",
            true,
        )
        .unwrap();
        assert_eq!((svlen[0].start, svlen[0].end), (30425651, 30425950));
        let end = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425700;END=30425900",
            true,
        )
        .unwrap();
        assert_eq!(end[0].end, 30425900);
        for info in [
            "SVTYPE=BND;CHR2=21;END2=30425700;END2=30425800",
            "SVTYPE=BND;CHR2=2;CHR2=21;END2=30425800",
        ] {
            let variants =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t<BND>\t.\t.\t{info}"), true).unwrap();
            assert_eq!(variants.len(), 2, "{info}");
            assert_eq!(variants[1].display_allele(), "N[21:30425800[", "{info}");
            assert_eq!(variants[1].mate_chr.as_deref(), Some("21"), "{info}");
        }
        let svtype = parse_vcf_line(
            "21\t30425650\tx\tN\tN.\t.\t.\tSVTYPE=DEL;SVTYPE=DUP;END=30425899",
            true,
        )
        .unwrap();
        assert_eq!(svtype[0].variant_class, VariantClass::Duplication);
        let bare_last = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;CHR2;END2=30425800",
            true,
        )
        .unwrap();
        assert_eq!(bare_last.len(), 1);
        assert_eq!(bare_last[0].record_allele_string(), "N/<BND>");
    }

    /// A signed coordinate is a number to Perl (`looks_like_number`, Parser.pm:475;
    /// `$start++` on the raw POS): Ensembl VEP 116.2 annotates `END=+30425899` as the end
    /// 30425899 and `SVLEN=+300` as 300, as a POS `+30425650` is 30425650; the sign it
    /// prints inside the Location (`21:30425651-+30425899`) is not reproduced.
    #[test]
    fn a_signed_coordinate_reads_as_its_number() {
        let end = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=+30425899",
            true,
        )
        .unwrap();
        assert_eq!((end[0].start, end[0].end), (30425651, 30425899));
        let pos = parse_vcf_line(
            "21\t+30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425899",
            true,
        )
        .unwrap();
        assert_eq!((pos[0].start, pos[0].end), (30425651, 30425899));
        let svlen = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=+300",
            true,
        )
        .unwrap();
        assert_eq!(svlen[0].end, 30425950);
    }

    /// A record with a parsed breakend is one breakpoint whatever its other ALTs: Ensembl
    /// VEP 116.2 writes `N.` and `N[21:30425700[` at the point `21:30425651` for
    /// `<DEL>,N[21:30425700[` with `SVTYPE=BND;END=30425899`, the `<DEL>` part lending
    /// nothing (`_parse_breakends` reads no base from it) and the END unread.
    #[test]
    fn a_record_with_a_parsed_breakend_is_one_breakpoint_whatever_its_other_alts() {
        let mixed = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>,N[21:30425700[\t.\t.\tSVTYPE=BND;END=30425899",
            true,
        )
        .unwrap();
        assert_eq!(mixed.len(), 2);
        for v in &mixed {
            assert_eq!(v.variant_class, VariantClass::Translocation);
            assert_eq!((v.start, v.end), (30425651, 30425651));
            assert_eq!(v.record_allele_string(), "N/<DEL>/N[21:30425700[");
        }
        assert_eq!(mixed[0].display_allele(), "N.");
        assert_eq!(mixed[1].display_allele(), "N[21:30425700[");
    }

    /// A bare `SVLEN` or `END` flag is undef to Perl (`get_end` tests both with `defined`,
    /// BaseVCF4.pm:302-307), so `<DUP>` with a bare `SVLEN` beside `END=30425899` ends at
    /// 30425899 and one with a bare `END` alone is the single base at POS + 1.
    #[test]
    fn a_bare_svlen_or_end_flag_is_undef() {
        let svlen_flag = parse_vcf_line(
            "21\t30425650\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;SVLEN;END=30425899",
            true,
        )
        .unwrap();
        assert_eq!(
            (svlen_flag[0].start, svlen_flag[0].end),
            (30425651, 30425899)
        );
        assert_eq!(svlen_flag[0].sv_len, None);
        let end_flag =
            parse_vcf_line("21\t30425650\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END", true).unwrap();
        assert_eq!((end_flag[0].start, end_flag[0].end), (30425651, 30425651));
    }

    /// Parser/VCF.pm:237-245 routes a line to the structural path when REF and the ALTs
    /// together are not a run of `ACGT` and INFO/SVTYPE is true or an ALT is written in
    /// a structural form: Ensembl VEP 116.2 writes `21:30425651-30425949 deletion` for
    /// `N G` and for `<*>` with `SVTYPE=DEL;END=30425949`, and for `A g` (a lowercase
    /// base is outside `ACGT`), while `A G` with the same INFO is the SNV
    /// `21:30425650 G missense_variant` and `N G` alone the SNV
    /// `coding_sequence_variant`. `SVTYPE=` and `SVTYPE=0` are false to Perl and route
    /// nothing; `<*>` alone is a sequence allele at POS.
    #[test]
    fn svtype_routes_a_record_whose_alleles_are_not_plain_acgt() {
        for (reference, alt) in [("N", "G"), ("N", "<*>"), ("A", "g")] {
            let v = parse_vcf_line(
                &format!("21\t30425650\tx\t{reference}\t{alt}\t.\t.\tSVTYPE=DEL;END=30425949"),
                true,
            )
            .unwrap();
            assert_eq!(v.len(), 1, "{reference} {alt}");
            assert_eq!(
                v[0].variant_class,
                VariantClass::StructuralDeletion,
                "{reference} {alt}"
            );
            assert_eq!(
                (v[0].start, v[0].end),
                (30425651, 30425949),
                "{reference} {alt}"
            );
            assert_eq!(v[0].display_allele(), "deletion", "{reference} {alt}");
        }
        for (reference, alt, info) in [
            ("A", "G", "SVTYPE=DEL;END=30425949"),
            ("N", "G", "."),
            ("N", "G", "SVTYPE=;END=30425949"),
            ("N", "G", "SVTYPE=0;END=30425949"),
            ("N", "<*>", "."),
        ] {
            let v = parse_vcf_line(
                &format!("21\t30425650\tx\t{reference}\t{alt}\t.\t.\t{info}"),
                true,
            )
            .unwrap();
            assert_eq!(v.len(), 1, "{reference} {alt} {info}");
            assert!(!v[0].is_structural, "{reference} {alt} {info}");
            assert_eq!(
                (v[0].start, v[0].end),
                (30425650, 30425650),
                "{reference} {alt} {info}"
            );
            assert_eq!(v[0].display_allele(), alt, "{reference} {alt} {info}");
        }
    }

    /// Parser/VCF.pm:473 tests the VCF 4.4 forms without `/i`, so a lowercase ALT is
    /// outside them and SVTYPE names its type; `get_SO_term`'s fall-through branches
    /// (Parser.pm:712-717) keep the type as written and `%SO_TERMS` has upper-case keys
    /// only. Ensembl VEP 116.2 writes `insertion` for `<ins:me:alu>` with `SVTYPE=INS`,
    /// `duplication` for `<del>` with `SVTYPE=DUP`, `copy_number_variation` for `<cnv:tr>`
    /// with `SVTYPE=CNV` and for `<DEL>,<dup>` with `SVTYPE=DEL`, and skips `<del>` alone
    /// and `<BND>` with `SVTYPE=bnd` as unsupported types.
    #[test]
    fn a_lowercase_symbolic_alt_is_outside_the_vcf44_forms() {
        for (alt, info, term, class) in [
            (
                "<ins:me:alu>",
                "SVTYPE=INS",
                "insertion",
                VariantClass::StructuralInsertion,
            ),
            (
                "<del>",
                "SVTYPE=DUP;END=30425949",
                "duplication",
                VariantClass::Duplication,
            ),
            (
                "<cnv:tr>",
                "SVTYPE=CNV;END=30425949",
                "copy_number_variation",
                VariantClass::CopyNumberVariation,
            ),
            (
                "<DEL>,<dup>",
                "SVTYPE=DEL;END=30425949",
                "copy_number_variation",
                VariantClass::CopyNumberVariation,
            ),
            (
                "<DEL>,<del>",
                "SVTYPE=DEL;END=30425949",
                "deletion",
                VariantClass::StructuralDeletion,
            ),
        ] {
            let v =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap();
            for variant in &v {
                assert_eq!(
                    variant.record_so_term.as_deref(),
                    Some(term),
                    "{alt} {info}"
                );
                assert_eq!(variant.variant_class, class, "{alt} {info}");
                assert!(!variant.vep_skip, "{alt} {info}");
            }
        }
        for (alt, info) in [
            ("<del>", "."),
            ("<BND>", "SVTYPE=bnd"),
            ("<del>", "SVTYPE=del;END=30425949"),
        ] {
            let v =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap();
            assert_eq!(v.len(), 1, "{alt} {info}");
            assert!(v[0].vep_skip, "{alt} {info}");
            assert_eq!(v[0].record_so_term, None, "{alt} {info}");
        }
        assert_eq!(sv_type_abbreviation("<del>", None), "del");
        assert_eq!(sv_type_abbreviation("<BND>", Some("bnd")), "bnd");
        assert_eq!(
            sv_type_abbreviation("<cn0>", None),
            "CNV",
            "`/^<?CN=?0>?$/` is upper case only"
        );
        assert_eq!(sv_type_abbreviation("<CN0>", None), "DEL");
    }

    /// A gVCF variant site lists `<NON_REF>` beside its real alleles; Ensembl VEP 116.2
    /// drops the line (`G/NON_REF is not a supported structural variant type`), and
    /// vep-rs reads the real alleles, numbered as they stand: `N G,<NON_REF>` is the SNV
    /// `G` at POS. A record whose only ALT is `<NON_REF>` is a reference block.
    #[test]
    fn a_gvcf_variant_site_is_its_real_alleles() {
        let site = parse_vcf_line("21\t30425650\tx\tA\tG,<NON_REF>\t.\t.\t.", true).unwrap();
        assert_eq!(site.len(), 1);
        assert_eq!(
            (
                site[0].display_allele(),
                site[0].allele_index,
                site[0].is_structural
            ),
            ("G".to_string(), 0, false)
        );
        let indel = parse_vcf_line("21\t30425650\tx\tN\tNA,<NON_REF>\t.\t.\t.", true).unwrap();
        assert_eq!(
            (indel.len(), indel[0].display_allele()),
            (1, "A".to_string())
        );
        let block =
            parse_vcf_line("21\t30425650\tx\tN\t<NON_REF>\t.\t.\tEND=30425700", true).unwrap();
        assert_eq!(
            (block.len(), block[0].variant_class),
            (1, VariantClass::ReferenceBlock)
        );
    }

    /// One StructuralVariationFeature carries the record (Parser/VCF.pm:560-574), so a
    /// record whose class is the breakpoint's gives every ALT the breakpoint's class and
    /// point: Ensembl VEP 116.2 writes `21:30425651 chromosome_breakpoint` for each part of
    /// `<BND>,<DEL>` with `SVTYPE=BND;END=30425949`, the END unread, with the allele string
    /// `N/<BND>/<DEL>`. A `*` allele stays in the string (`*/<DEL>`) and names no variant;
    /// a `.` ALT names none whatever SVTYPE says.
    #[test]
    fn a_breakpoint_record_puts_every_alt_at_the_point_and_keeps_star_in_its_string() {
        let v = parse_vcf_line(
            "21\t30425650\tx\tN\t<BND>,<DEL>\t.\t.\tSVTYPE=BND;END=30425949",
            true,
        )
        .unwrap();
        assert_eq!(v.len(), 2);
        for variant in &v {
            assert_eq!(variant.variant_class, VariantClass::Translocation);
            assert_eq!((variant.start, variant.end), (30425651, 30425651));
            assert_eq!(variant.display_allele(), "chromosome_breakpoint");
            assert_eq!(variant.record_allele_string(), "N/<BND>/<DEL>");
        }
        let star = parse_vcf_line(
            "21\t30425650\tx\tN\t*,<DEL>\t.\t.\tSVTYPE=DEL;END=30425949",
            true,
        )
        .unwrap();
        assert_eq!(star.len(), 1);
        assert_eq!(star[0].record_allele_string(), "*/<DEL>");
        assert_eq!(star[0].allele_index, 1);
        for info in ["SVTYPE=BND", "SVTYPE=DEL;END=30425949", "."] {
            let err = parse_vcf_line(&format!("21\t30425650\tx\tN\t.\t.\t.\t{info}"), true)
                .unwrap_err()
                .to_string();
            assert!(err.contains("No alternate allele"), "{info}: {err}");
        }
        let err = parse_vcf_line("21\t30425650\tx\tN\tG\t.\t.\tSVTYPE=FOO;END=30425949", true)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("FOO is not a supported structural variant type"),
            "{err}"
        );
    }

    /// ALLELE_NUM numbers the parts of Ensembl VEP's allele string from 1
    /// (StructuralVariationOverlap.pm:64-77): a bare `<BND>` is `N/<BND>`, so its row
    /// prints 2 (116.2 writes `ALLELE_NUM=2` beside the reference allele's 1); `.N` is
    /// not prefixed and prints 1; a record with parsed breakends prints 0 for its own
    /// breakend and 1 upward for its mates; every other class prints the ALT's position.
    #[test]
    fn allele_num_follows_ensembl_veps_allele_string() {
        let bare = parse_vcf_line("21\t30425650\tx\tN\t<BND>\t.\t.\tSVTYPE=BND", true).unwrap();
        assert_eq!(bare[0].allele_num, Some(2));
        let dot = parse_vcf_line("21\t30425650\tx\tN\t.N\t.\t.\tSVTYPE=BND", true).unwrap();
        assert_eq!(dot[0].allele_num, Some(1));
        let paired = parse_vcf_line(
            "21\t30425650\tx\tN\tN[21:30425700[,N]22:100]\t.\t.\tSVTYPE=BND",
            true,
        )
        .unwrap();
        assert_eq!(
            paired.iter().map(|v| v.allele_num).collect::<Vec<_>>(),
            [Some(0), Some(1), Some(2)]
        );
        let del = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>,<DUP>\t.\t.\tSVTYPE=CNV;END=30425949",
            true,
        )
        .unwrap();
        assert_eq!(
            del.iter().map(|v| v.allele_num).collect::<Vec<_>>(),
            [None, None]
        );
    }

    /// `get_SO_term`'s regexes decide the type of a shape at their edges: `/CN=?[0-9]/i`
    /// needs a digit, the list normalisation strips one bracket a side, the mobile-element
    /// pattern takes one colon before the element, and `s/\:.+//g` needs a character after
    /// the colon. Ensembl VEP 116.2 skips `<CN=X>`, `<<DEL>>,<DEL>` and `<DEL:>` as
    /// unsupported types (`CN=X`, `DEL/DEL`, `DEL:`), writes `copy_number_variation` for
    /// `<<CN0>>` and `mobile_element_insertion` for `<INS:ME::ALU>` with `SVTYPE=INS`, and
    /// skips `<NON_REF>,<NON_REF>` as `NON_REF`, which vep-rs reads as a reference block
    /// of two alleles, one intergenic row each; under `SVTYPE=DEL` it writes two `deletion`
    /// rows per transcript for `G,<NON_REF>`, the block an allele of the record.
    #[test]
    fn the_type_follows_get_so_terms_regexes_at_their_edges() {
        let parse = |alt: &str, info: &str| {
            parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap()
        };
        for (alt, info, abbrev) in [
            ("<CN=X>", "END=30425949", "CN=X"),
            ("<<DEL>>,<DEL>", "END=30425949", "DEL/DEL"),
            ("<DEL:>", "END=30425949", "DEL:"),
        ] {
            let v = parse(alt, info);
            assert!(v.iter().all(|v| v.vep_skip), "{alt}");
            assert!(v.iter().all(|v| v.record_so_term.is_none()), "{alt}");
            assert_eq!(sv_type_abbreviation(&alt.replace(',', "/"), None), abbrev);
        }
        for (alt, info, term) in [
            ("<CN=X>", "SVTYPE=DEL;END=30425949", "deletion"),
            ("<<DEL>>", "END=30425949", "deletion"),
            ("<<CN0>>", "END=30425949", "copy_number_variation"),
            ("<CN0>", "END=30425949", "deletion"),
            ("<INS:ME::ALU>", "SVTYPE=INS", "mobile_element_insertion"),
            ("<INS:ME:ALU>", "SVTYPE=INS", "Alu_insertion"),
            ("<DEL:X>", "END=30425949", "deletion"),
            ("<CNX>,<CN3>", "END=30425949", "copy_number_variation"),
        ] {
            let v = parse(alt, info);
            assert!(
                !v.is_empty() && v.iter().all(|v| !v.vep_skip),
                "{alt} {info}"
            );
            assert!(
                v.iter().all(|v| v.record_so_term.as_deref() == Some(term)),
                "{alt} {info}: {:?}",
                v[0].record_so_term
            );
        }
        let site = parse("G,<NON_REF>", "SVTYPE=DEL;END=30425949");
        assert_eq!(
            site.len(),
            2,
            "a true SVTYPE keeps the block as a deletion allele"
        );
        for variant in &site {
            assert_eq!(variant.record_so_term.as_deref(), Some("deletion"));
            assert_eq!(variant.record_allele_string(), "G/<NON_REF>");
        }
        let blocks = parse("<NON_REF>,<NON_REF>", "END=30425700");
        assert_eq!(blocks.len(), 2);
        for (variant, index) in blocks.iter().zip(0..) {
            assert_eq!(variant.variant_class, VariantClass::ReferenceBlock);
            assert!(variant.vep_skip);
            assert_eq!(variant.allele_index, index);
            assert_eq!(variant.record_allele_string(), "<NON_REF>/<NON_REF>");
        }
    }

    /// Perl's list-context `split` drops a trailing empty ALT field (BaseVCF4.pm:489) and a
    /// trailing empty type part (Parser.pm:669): Ensembl VEP 116.2 writes the one allele `G`
    /// for `A G,`, a `deletion` over END for `<DEL>,` with `SVTYPE=DEL`, and a `deletion`
    /// for `N.` with `SVTYPE=DEL/`.
    #[test]
    fn a_trailing_empty_alt_field_names_no_allele() {
        let seq = parse_vcf_line("21\t30425650\tx\tA\tG,\t.\t.\t.", true).unwrap();
        assert_eq!(seq.len(), 1);
        assert_eq!(
            (seq[0].display_allele(), seq[0].allele_index),
            ("G".to_string(), 0)
        );
        assert_eq!(seq[0].uploaded_allele_string.as_deref(), Some("A/G"));
        for (alt, info) in [
            ("<DEL>,", "SVTYPE=DEL;END=30425949"),
            ("N.", "SVTYPE=DEL/;END=30425949"),
        ] {
            let v =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap();
            assert_eq!(v.len(), 1, "{alt} {info}");
            assert_eq!(
                v[0].record_so_term.as_deref(),
                Some("deletion"),
                "{alt} {info}"
            );
            assert!(!v[0].vep_skip, "{alt} {info}");
            assert_eq!((v[0].start, v[0].end), (30425651, 30425949), "{alt} {info}");
        }
    }

    /// A gVCF site's real allele keeps its written position whichever side of `<NON_REF>`
    /// it stands on, so `<NON_REF>,G` prints ALLELE_NUM 2 for `G`, the number a genotype
    /// `0/2` names.
    #[test]
    fn a_gvcf_sites_allele_keeps_its_written_position() {
        let after = parse_vcf_line("21\t30425650\tx\tA\t<NON_REF>,G\t.\t.\t.", true).unwrap();
        assert_eq!(after.len(), 1);
        assert_eq!(
            (after[0].display_allele(), after[0].allele_index),
            ("G".to_string(), 1)
        );
        let before = parse_vcf_line("21\t30425650\tx\tA\tG,<NON_REF>\t.\t.\t.", true).unwrap();
        assert_eq!(before[0].allele_index, 0);
    }

    /// BaseVCF4.pm:581-583 builds the INFO hash key by key as written (`$info_data{$key} =
    /// $value`), the last occurrence winning: a bare `END` written last leaves END undef, and
    /// `end=` or ` END=` is not `END`, so a deletion carrying only those has no end and is
    /// Ensembl VEP's `deletion looks incomplete` (Parser/VCF.pm:555-558); vep-rs ends such a
    /// record at its start.
    #[test]
    fn an_info_key_matches_exactly_and_a_bare_key_written_last_is_undef() {
        let dup =
            parse_vcf_line("21\t1000\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=1100;END", true).unwrap();
        assert_eq!((dup[0].start, dup[0].end), (1001, 1001));
        for info in [
            "svtype=DEL;end=30425949",
            "SVTYPE=DEL; END=30425949",
            "SVTYPE=DEL;END =30425949",
        ] {
            let v =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t<DEL>\t.\t.\t{info}"), true).unwrap();
            assert_eq!((v[0].start, v[0].end), (30425651, 30425651), "{info}");
            assert_eq!(v[0].record_so_term.as_deref(), Some("deletion"), "{info}");
        }
        let exact = parse_vcf_line(
            "21\t30425650\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425949",
            true,
        )
        .unwrap();
        assert_eq!(exact[0].end, 30425949);
    }

    /// `_parse_breakends`' match tolerates whitespace around the brackets, the name and the
    /// colon, and BaseVCF4.pm:582 cuts an INFO value at a second `=`: Ensembl VEP 116.2
    /// writes the own breakend `N.` (ALLELE_NUM 0) and a mate row for `N[ 21 : 30425700 [`,
    /// for `<BND>` with `CHR2= 21;END2=30425700` (Allele `N[ 21:30425700[`) and for `<BND>`
    /// with `CHR2=21;END2=30425700=5` (Allele `N[21:30425700[`).
    #[test]
    fn a_breakend_mate_is_read_through_the_whitespace_perl_allows() {
        for (alt, info, mate_text) in [
            ("N[ 21 : 30425700 [", "SVTYPE=BND", "N[ 21 : 30425700 ["),
            (
                "<BND>",
                "SVTYPE=BND;CHR2= 21;END2=30425700",
                "N[ 21:30425700[",
            ),
            (
                "<BND>",
                "SVTYPE=BND;CHR2=21;END2=30425700=5",
                "N[21:30425700[",
            ),
        ] {
            let v =
                parse_vcf_line(&format!("21\t30425650\tx\tN\t{alt}\t.\t.\t{info}"), true).unwrap();
            assert_eq!(v.len(), 2, "{alt} {info}");
            assert_eq!(
                (v[0].display_allele(), v[0].allele_num),
                ("N.".to_string(), Some(0)),
                "{alt}"
            );
            assert_eq!(
                (v[1].display_allele(), v[1].allele_num),
                (mate_text.to_string(), Some(1)),
                "{alt}"
            );
            assert_eq!(
                (v[1].mate_chr.as_deref(), v[1].mate_pos),
                (Some("21"), Some(30425700)),
                "{alt}"
            );
        }
    }

    /// A unit of no bases under RUC is the empty allele (`'' x 5`, Parser/VCF.pm:419), no
    /// bases and so the loss of the run; under RB it divides by zero (`:414`) and keeps the
    /// symbolic reading; a count is read as Perl numifies it (`RUC=2abc` repeats twice, a
    /// loss of three bases on a nine-base run).
    #[test]
    fn a_zero_length_repeat_unit_is_read_as_perl_expands_it() {
        let run = "SVTYPE=CNV;END=30425949;SVLEN=300;";
        for (info, bases) in [("RUS=,AC;RUC=5", Some(0)), ("RUS=,AC;RB=12", None)] {
            let v = parse_vcf_line(
                &format!("21\t30425650\tx\tN\t<CNV:TR>\t.\t.\t{run}{info}"),
                true,
            )
            .unwrap();
            assert_eq!(v[0].tr_alt_bases, bases, "{info}");
        }
        let numified = parse_vcf_line(
            "21\t30425587\tx\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=2abc",
            true,
        )
        .unwrap();
        assert_eq!(numified[0].tr_alt_bases, Some(6));
    }

    /// The record-wide strings (the uploaded allele string, the record's working string,
    /// the line name and the raw line) are one buffer per record shared by its alleles.
    #[test]
    fn the_record_strings_are_shared_by_its_alleles() {
        let shared = |variants: &[InputVariant]| {
            let first = &variants[0];
            variants.iter().skip(1).all(|v| {
                Arc::ptr_eq(
                    v.uploaded_allele_string.as_ref().unwrap(),
                    first.uploaded_allele_string.as_ref().unwrap(),
                ) && Arc::ptr_eq(
                    v.input_line_name.as_ref().unwrap(),
                    first.input_line_name.as_ref().unwrap(),
                ) && Arc::ptr_eq(
                    v.raw_input.as_ref().unwrap(),
                    first.raw_input.as_ref().unwrap(),
                ) && match (
                    &v.record_allele_string_multi,
                    &first.record_allele_string_multi,
                ) {
                    (Some(a), Some(b)) => Arc::ptr_eq(a, b),
                    (None, None) => true,
                    _ => false,
                }
            })
        };
        let seq = parse_vcf_line("21\t1000\tx\tA\tG,T,C\t.\t.\t.", true).unwrap();
        assert!(seq.len() == 3 && shared(&seq));
        let sv = parse_vcf_line(
            "21\t1000\tx\tN\t<DEL>,<DUP>\t.\t.\tSVTYPE=CNV;END=1100",
            true,
        )
        .unwrap();
        assert!(sv.len() == 2 && shared(&sv));
        let bnd = parse_vcf_line(
            "21\t1000\tx\tN\tN[21:1700[,N]22:100]\t.\t.\tSVTYPE=BND",
            true,
        )
        .unwrap();
        assert!(bnd.len() == 3 && shared(&bnd));
    }

    /// Both VCF readers read a record through the shared `vep_core` modules (INFO, the ALT
    /// list, routing, the structural record), so over the shapes that decide the path and
    /// the alleles (an `N` reference base, `<*>`, a gVCF site in either order, a block list,
    /// a lowercase ALT, a mixed list, a breakpoint list, bracket and CHR2 mates, dot
    /// breakends, a bare `<BND>` with MATEID, repeated and bare INFO keys, a trailing
    /// empty ALT, POS 0, a POS spelled `1e3` or `1000.0`, a `*`-only record, a `*` beside
    /// sequence alleles in either order, `<CNV:TR>`
    /// counts) every variant's path, class, skip mark, term, span, allele index and
    /// number, Allele, mate, record string, allele-less mark and repeat length agree; a
    /// POS past the ceiling, and a structural record's POS at it, are refused by both.
    #[test]
    fn the_vep_io_reader_routes_as_this_parser_does() {
        use vep_io::input::vcf::VcfParser as Reader;
        use vep_io::InputParser;
        let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let summary = |v: &InputVariant| {
            (
                (
                    v.is_structural,
                    v.variant_class,
                    v.vep_skip,
                    v.record_so_term.clone(),
                ),
                (v.start, v.end, v.allele_index, v.allele_num),
                (
                    v.display_allele(),
                    v.mate_id.is_some(),
                    v.is_single_breakend,
                ),
                (v.mate_chr.clone(), v.mate_pos, v.tr_alt_bases, v.no_allele),
                // The record strings of the sequence path are the command-line parser's
                // alone; the structural ones come from the shared record builder.
                (
                    v.is_structural
                        .then(|| v.uploaded_allele_string.as_deref().map(str::to_string)),
                    v.is_structural.then(|| v.record_allele_string()),
                    v.input_line_name.as_deref().map(str::to_string),
                ),
            )
        };
        let from_reader = |record: &str| {
            let data = format!("{header}{record}\n");
            let mut reader =
                Reader::new(Box::new(std::io::Cursor::new(data.into_bytes()))).unwrap();
            let mut variants = Vec::new();
            while let Some(v) = reader.next_variant() {
                variants.push(v.map(|v| summary(&v)).map_err(|e| e.to_string()));
            }
            variants
        };
        for record in [
            "21\t1000\tx\tN\tG\t.\t.\tSVTYPE=DEL;END=1100",
            "21\t1000\tx\tA\tG\t.\t.\tSVTYPE=DEL;END=1100",
            "21\t1000\tx\tN\tG\t.\t.\t.",
            "21\t1000\tx\tN\t<*>\t.\t.\t.",
            "21\t1000\tx\tN\t<*>\t.\t.\tSVTYPE=DEL;END=1100",
            "21\t1000\tx\tA\tG,<NON_REF>\t.\t.\t.",
            "21\t1000\tx\tA\t<NON_REF>,G\t.\t.\t.",
            "21\t1000\tx\tN\tG,<NON_REF>\t.\t.\tSVTYPE=DEL;END=1100",
            "21\t1000\tx\tN\t<NON_REF>\t.\t.\tEND=1100",
            "21\t1000\tx\tN\t<NON_REF>,<NON_REF>\t.\t.\tEND=1100",
            "21\t1000\tx\tN\t<del>\t.\t.\tEND=1100",
            "21\t1000\tx\tN\t<del>\t.\t.\tSVTYPE=DUP;END=1100",
            "21\t1000\tx\tA\tG,<DEL>\t.\t.\tSVTYPE=DEL;END=1100",
            "21\t1000\tx\tC\tG,<INV>\t.\t.\tSVTYPE=DEL;END=2000;SVLEN=-1000",
            "21\t1000\tx\tN\t<BND>,<DEL>\t.\t.\tSVTYPE=BND;END=1100",
            "21\t1000\tx\tN\t<INS>,<DEL>\t.\t.\tEND=1100",
            "21\t1000\tx\tN\t<INS:ME::ALU>\t.\t.\tSVTYPE=INS",
            "21\t1000\tx\tN\t<CN=X>\t.\t.\tEND=1100",
            "21\t1000\tx\tN\t<<DEL>>,<DEL>\t.\t.\tEND=1100",
            "21\t1000\tx\tN\t<DEL:>\t.\t.\tEND=1100",
            "21\t1000\tx\tN\t<CNX>,<CN3>\t.\t.\tEND=1100",
            "21\t1000\tx\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=1009;SVLEN=9;RUS=AGC;RUC=4",
            "21\t1000\tx\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=1009;SVLEN=9;RUS=AGC;RUC=2abc",
            "21\t1000\tx\tN\t<CNV:TR>\t.\t.\tSVTYPE=CNV;END=1009;SVLEN=9;RUS=,AC;RUC=5",
            "21\t1000\tx\tN\tN[21:1700[,N]22:100]\t.\t.\tSVTYPE=BND",
            "21\t1000\tx\tN\tN[ 21 : 1700 [\t.\t.\tSVTYPE=BND",
            "21\t1000\tx\tA\t<BND>\t.\t.\tSVTYPE=BND;CHR2= 21;END2=1700=5",
            "21\t1000\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;MATEID=m1",
            "21\t1000\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;SVLEN=100",
            "21\t1000\tx\tN\t.N\t.\t.\tSVTYPE=BND;MATEID=m1",
            "21\t1000\tx\tN\tN.,N[21:1700[\t.\t.\tSVTYPE=BND",
            "21\t1000\tx\tN\t*,<DEL>\t.\t.\tSVTYPE=DEL;END=1100",
            "21\t1000\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=1100;END",
            "21\t1000\tx\tN\t<DUP>\t.\t.\tSVTYPE=DUP;END=1100;END=1200",
            "21\t1000\tx\tN\t<DEL>\t.\t.\tsvtype=DEL;end=1100",
            "21\t1000\tx\tA\tG,\t.\t.\t.",
            "21\t1000\tx\tN\t<DEL>,\t.\t.\tSVTYPE=DEL;END=1100",
            "21\t1000\tx\tN\tN.\t.\t.\tSVTYPE=DEL/;END=1100",
            "21\t0\tx\tA\tAT\t.\t.\t.",
            "21\t1e3\tx\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=1100",
            "21\t1000.0\tx\tA\tG\t.\t.\t.",
            "21\t1000\tx\tA\tg\t.\t.\t.",
            "21\t1000\tx\tagc\ta\t.\t.\t.",
            "21\t1000\tx\tN\t*\t.\t.\t.",
            "21\t1000\tx\tNAC\t*,*\t.\t.\tSVTYPE=DEL;END=1100",
            "21\t1000\tx\tN\tN[21:1700[,<NON_REF>\t.\t.\t.",
            "21\t1000\tx\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;END2=1700 ",
            "21\t30425650\tx\tA\tG,*\t.\t.\t.",
            "21\t30425650\tx\tAG\tAC,*\t.\t.\t.",
            "21\t30425650\tx\tAG\tTC,*\t.\t.\t.",
            "21\t30425650\tx\tAGC\tATT,*\t.\t.\t.",
            "21\t30425650\tx\tAG\tA,*\t.\t.\t.",
            "21\t30425650\tx\tA\tAG,*\t.\t.\t.",
            "21\t30425650\tx\tAG\t*,AC\t.\t.\t.",
            "21\t30425650\tx\tAG\tAC,AT,*\t.\t.\t.",
            "21\t30425650\tx\tA\t*,G\t.\t.\t.",
            "21\t30425650\tx\tAG\t*,A\t.\t.\t.",
            "21\t30425650\tx\tA\tG,T,*\t.\t.\t.",
        ] {
            let cli: Vec<Result<_, String>> = parse_vcf_line(record, false)
                .unwrap()
                .iter()
                .map(|v| Ok(summary(v)))
                .collect();
            assert_eq!(from_reader(record), cli, "{record}");
        }
        // A record whose type has no term and whose ALTs are all sequence alleles is an
        // error here and a skipped record there, both writing nothing.
        let record = "21\t1000\tx\tN\tG\t.\t.\tSVTYPE=FOO;END=1100";
        assert!(parse_vcf_line(record, false).is_err());
        assert!(from_reader(record).is_empty());
        // A position past the ceiling is refused by both, as is a structural record's
        // position at it.
        for (pos, alt, info) in [
            (MAX_COORD + 1, "G", "."),
            (u64::MAX, "G", "."),
            (MAX_COORD, "<DEL>", "SVTYPE=DEL;SVLEN=-300"),
        ] {
            let record = format!("21\t{pos}\tx\tA\t{alt}\t.\t.\t{info}");
            assert!(parse_vcf_line(&record, false).is_err(), "{pos}");
            let read = from_reader(&record);
            assert!(read.len() == 1 && read[0].is_err(), "{pos}: {read:?}");
        }
        assert!(parse_vcf_line("21\t18446744073709551616\tx\tA\tG\t.\t.\t.", false).is_err());
        assert_eq!(
            from_reader("21\t18446744073709551616\tx\tA\tG\t.\t.\t.").len(),
            1
        );
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
    /// there whatever the class. A breakpoint's deleted END names its mate with CHR2
    /// (Parser/VCF.pm:493): a run of digits past the ceiling is a mate whose text keeps
    /// the digits (Ensembl VEP 116.2 writes `N[21:18446744073709551615[`) and whose
    /// position saturates at the ceiling; `-5` or `abc` names none and the record's
    /// string is `N/N[21:-5[`.
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
            let digits = end.bytes().all(|b| b.is_ascii_digit());
            assert_eq!(breakpoint.len(), if digits { 2 } else { 1 }, "{end}");
            assert_eq!(
                (breakpoint[0].start, breakpoint[0].end),
                (101, 101),
                "{end}"
            );
            assert_eq!(
                breakpoint[0].record_allele_string(),
                format!("N/N[21:{end}["),
                "{end}"
            );
            if digits {
                assert_eq!(breakpoint[0].display_allele(), "N.", "{end}");
                assert_eq!(breakpoint[1].display_allele(), format!("N[21:{end}["));
                assert_eq!(breakpoint[1].mate_pos, Some(MAX_COORD), "{end}");
            } else {
                assert_eq!(
                    breakpoint[0].display_allele(),
                    "chromosome_breakpoint",
                    "{end}"
                );
            }
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
        // Perl VEP does not filter <*>: the SV regex explicitly excludes it, so without
        // INFO/SVTYPE it is a non-SV allele at position POS.
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

    /// The anchor chop counts a `*` as an allele of one base (Parser/VCF.pm:296-298) and
    /// leaves it out of the first-base comparison (:306), so Ensembl VEP 116.2 chops the
    /// shared first base of a record whose other alleles have the REF's length: at 30425650
    /// it writes `AG AC,*` as `21:30425651 C` (allele string `G/C/*`), `AGC ATT,*` as
    /// `21:30425651-30425652 TT` (`GC/TT/*`), `AG *,AC` as `21:30425651 C` (ALLELE_NUM 2,
    /// `G/*/C`) and `AG AC,AT,*` as `21:30425651 C` and `T` (`G/C/T/*`); the chop needs the
    /// shared base (`AG TC,*` stays `21:30425650-30425651 TC`), a one-base REF has nothing
    /// to chop (`A G,*`, `A *,G`, `A G,T,*` stay at 30425650) and without the `*` the
    /// same-length record is not an indel (`AG AC,AT` stays `21:30425650-30425651`); an
    /// indel beside a `*` is chopped as before (`AG A,*` `21:30425651 -`, `A AG,*`
    /// `21:30425650-30425651 G`, `AG *,A` `21:30425651 -` with ALLELE_NUM 2). A `.` beside the
    /// alleles is compared and never shares the base, so `A AT,.` keeps its anchor
    /// (`21:30425650 AT`, as Ensembl VEP 116.2 writes it beside its own rows under `.`).
    #[test]
    fn a_star_counts_as_an_allele_of_one_base_in_the_anchor_chop() {
        for (reference, alts, allele_string, rows) in [
            ("A", "AT,.", "A/AT/.", vec![("21:30425650", "AT", 0)]),
            (
                "ATT",
                "AT,.",
                "ATT/AT/.",
                vec![("21:30425650-30425652", "AT", 0)],
            ),
            ("A", "G,*", "A/G/*", vec![("21:30425650", "G", 0)]),
            ("AG", "AC,*", "G/C/*", vec![("21:30425651", "C", 0)]),
            (
                "AG",
                "TC,*",
                "AG/TC/*",
                vec![("21:30425650-30425651", "TC", 0)],
            ),
            (
                "AGC",
                "ATT,*",
                "GC/TT/*",
                vec![("21:30425651-30425652", "TT", 0)],
            ),
            ("AG", "A,*", "G/-/*", vec![("21:30425651", "-", 0)]),
            ("A", "AG,*", "-/G/*", vec![("21:30425650-30425651", "G", 0)]),
            ("AG", "*,AC", "G/*/C", vec![("21:30425651", "C", 1)]),
            (
                "AG",
                "AC,AT,*",
                "G/C/T/*",
                vec![("21:30425651", "C", 0), ("21:30425651", "T", 1)],
            ),
            ("A", "*,G", "A/*/G", vec![("21:30425650", "G", 1)]),
            ("AG", "*,A", "G/*/-", vec![("21:30425651", "-", 1)]),
            (
                "A",
                "G,T,*",
                "A/G/T/*",
                vec![("21:30425650", "G", 0), ("21:30425650", "T", 1)],
            ),
            (
                "AG",
                "AC,AT",
                "AG/AC/AT",
                vec![
                    ("21:30425650-30425651", "AC", 0),
                    ("21:30425650-30425651", "AT", 1),
                ],
            ),
        ] {
            let record = format!("21\t30425650\t.\t{reference}\t{alts}\t.\t.\t.");
            let variants = parse_vcf_line(&record, false).unwrap();
            let read: Vec<(String, String, usize)> = variants
                .iter()
                .map(|v| (v.location(), v.display_allele(), v.allele_index))
                .collect();
            let expected: Vec<(String, String, usize)> = rows
                .iter()
                .map(|(location, allele, index)| (location.to_string(), allele.to_string(), *index))
                .collect();
            assert_eq!(read, expected, "{reference} {alts}");
            for v in &variants {
                assert_eq!(
                    v.record_allele_string(),
                    allele_string,
                    "{reference} {alts}"
                );
                assert!(!v.minimised, "{reference} {alts}");
            }
        }
    }

    /// A `*` beside an indel keeps the record multi-allelic: `minimise_alleles`
    /// (Parser.pm:978-984) leaves a three-part string as it is, so Ensembl VEP 116.2 writes
    /// `21 30425650 x ATT AT,*` at `21:30425651-30425652` with the Allele `T` (the
    /// anchor-chopped `TT/T/*`), `AT ATT,*` at `21:30425651` with `TT`, and `ATTG ATC,*` at
    /// `21:30425651-30425653` with `TC`; without the `*` the bi-allelic indel is minimised
    /// from both ends (`21:30425652 -`, `21:30425651-30425652 T`, `21:30425652-30425653 C`).
    #[test]
    fn a_star_beside_an_indel_keeps_the_anchor_chopped_alleles() {
        let read = |record: &str| {
            let v = parse_vcf_line(record, false).unwrap();
            assert_eq!(v.len(), 1, "{record}");
            (
                v[0].location(),
                v[0].display_allele(),
                v[0].record_allele_string(),
                v[0].minimised,
            )
        };
        for (record, location, allele, string) in [
            (
                "21\t30425650\tx\tATT\tAT,*\t.\t.\t.",
                "21:30425651-30425652",
                "T",
                "TT/T/*",
            ),
            (
                "21\t30425650\tx\tAT\tATT,*\t.\t.\t.",
                "21:30425651",
                "TT",
                "T/TT/*",
            ),
            (
                "21\t30425650\tx\tATTG\tATC,*\t.\t.\t.",
                "21:30425651-30425653",
                "TC",
                "TTG/TC/*",
            ),
        ] {
            let (loc, al, st, minimised) = read(record);
            assert_eq!(
                (loc.as_str(), al.as_str(), st.as_str()),
                (location, allele, string),
                "{record}"
            );
            assert!(!minimised, "{record}");
        }
        for (record, location, allele) in [
            ("21\t30425650\tx\tATT\tAT\t.\t.\t.", "21:30425652", "-"),
            (
                "21\t30425650\tx\tAT\tATT\t.\t.\t.",
                "21:30425651-30425652",
                "T",
            ),
            (
                "21\t30425650\tx\tATTG\tATC\t.\t.\t.",
                "21:30425652-30425653",
                "C",
            ),
        ] {
            let (loc, al, _, minimised) = read(record);
            assert_eq!((loc.as_str(), al.as_str()), (location, allele), "{record}");
            assert!(minimised, "{record}");
        }
    }

    /// A line whose every ALT is `*` names no allele: Ensembl VEP 116.2 writes
    /// `21 30425650 x N * . . .` and `N *,*` with no row, the VCF line without CSQ and a
    /// JSON object (`"allele_string":"N/*"`, start and end 30425650,
    /// `"most_severe_consequence":"?"`), with or without a true SVTYPE. The one variant
    /// carries the record and is never annotated; a `*` beside an allele stays filtered.
    #[test]
    fn a_star_only_record_is_carried_as_the_allele_less_feature() {
        for (line, allele_string) in [
            ("21\t100\t.\tA\t*\t.\tPASS\t.", "A/*"),
            ("21\t100\t.\tNAC\t*,*\t.\tPASS\t.", "NAC/*/*"),
            ("21\t100\t.\tN\t*\t.\tPASS\tSVTYPE=DEL;END=400", "N/*"),
        ] {
            let variants = parse_vcf_line(line, true).unwrap();
            assert_eq!(variants.len(), 1, "{line}");
            let v = &variants[0];
            assert!(v.no_allele, "{line}");
            assert_eq!(v.allele_string, allele_string, "{line}");
            assert_eq!(
                (v.start, v.end),
                (100, 99 + v.ref_allele.len() as u64),
                "{line}"
            );
            assert!(!v.is_structural && v.transcript_consequences.is_empty());
            assert_eq!(v.record_allele_string(), allele_string, "{line}");
        }
        assert_eq!(
            parse_vcf_line("21\t100\t.\tA\t*\t.\tPASS\t.", true).unwrap()[0].uploaded_variation(),
            "21_100_A/*"
        );
        let beside = parse_vcf_line("21\t100\t.\tA\tG,*\t.\tPASS\t.", true).unwrap();
        assert_eq!(beside.len(), 1);
        assert!(!beside[0].no_allele);
        assert_eq!(beside[0].alt_allele(), b"G");
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
