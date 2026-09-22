// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Co-located variant matching and frequency extraction.
//!
//! Matches input variants against cached known variants (e.g., dbSNP, gnomAD)
//! by genomic position, extracts allele-specific population frequencies, and
//! builds the Extra output fields (AF, CLIN_SIG, PUBMED, etc.) for VEP output.
//! Implements the `--check_existing` pipeline stage.

use std::collections::HashMap;

use rustc_hash::FxHashMap;
use vep_core::variant::{ColocatedVariant, InputVariant};
use vep_core::variation::CachedVariation;

/// Find co-located known variants for an input variant.
///
/// Uses simple position-based matching: any cached variant at the same
/// start position is reported as co-located. This matches the default
/// Perl VEP behaviour with `--check_existing`.
pub fn find_colocated_variants(
    variant: &InputVariant,
    variations: &[CachedVariation],
    position_index: &FxHashMap<u64, Vec<usize>>,
) -> Vec<ColocatedVariant> {
    let Some(indices) = position_index.get(&variant.start) else {
        return Vec::new();
    };

    let mut colocated = Vec::new();
    for &idx in indices {
        let cv = &variations[idx];
        if cv.failed != 0 {
            continue;
        }
        colocated.push(cached_to_colocated(cv, variant));
    }

    colocated
}

/// Convert a CachedVariation to a ColocatedVariant, extracting allele-specific
/// frequencies for the input variant's alt allele.
///
/// When doing position-based matching (no allele check), the input allele may
/// not appear in the cached variant's frequency strings. In that case the first
/// non-reference allele frequency available is used.
fn cached_to_colocated(cv: &CachedVariation, variant: &InputVariant) -> ColocatedVariant {
    let alt_allele = variant.display_allele();

    let mut frequencies = HashMap::new();
    for (pop, freq_str) in &cv.frequencies {
        let freq = extract_allele_frequency(freq_str, &alt_allele)
            .or_else(|| extract_first_frequency(freq_str));
        if let Some(f) = freq {
            frequencies.insert(pop.clone(), f);
        }
    }

    let clin_sig = cv
        .clin_sig
        .as_deref()
        .filter(|s| !s.is_empty())
        .map(|s| s.split(',').map(|t| t.trim().to_string()).collect())
        .unwrap_or_default();

    let pubmed = cv
        .pubmed
        .as_deref()
        .filter(|s| !s.is_empty())
        .map(|s| s.split(',').map(|t| t.trim().to_string()).collect())
        .unwrap_or_default();

    ColocatedVariant {
        id: cv.variation_name.clone(),
        start: cv.start,
        end: cv.end,
        allele_string: Some(cv.allele_string.clone()).filter(|s| !s.is_empty()),
        strand: cv.strand,
        somatic: cv.somatic != 0,
        clin_sig,
        phenotype: cv.phenotype_or_disease != 0,
        frequencies,
        pubmed,
    }
}

/// Parse a frequency string like "T:0.003,C:0.997" or "T:0.003" and extract
/// the frequency for the given allele.
///
/// Returns `None` if the allele is not found in the frequency string.
pub fn extract_allele_frequency(freq_str: &str, allele: &str) -> Option<f64> {
    for entry in freq_str.split(',') {
        let entry = entry.trim();
        if let Some((allele_part, freq_part)) = entry.split_once(':') {
            if allele_part == allele {
                return freq_part.parse::<f64>().ok();
            }
        }
    }
    None
}

/// Extract the first available frequency from a frequency string,
/// regardless of allele. Used as a fallback when the specific input
/// allele is not found (e.g., due to strand differences).
fn extract_first_frequency(freq_str: &str) -> Option<f64> {
    for entry in freq_str.split(',') {
        let entry = entry.trim();
        if let Some((_allele_part, freq_part)) = entry.split_once(':') {
            if let Ok(f) = freq_part.parse::<f64>() {
                return Some(f);
            }
        }
    }
    None
}

/// Render `value` as decimal ASCII into `buf`, returning the populated slice.
///
/// `buf` must hold at least 20 bytes (`u64::MAX` is 20 digits). Byte-identical
/// to `u64`'s `Display`, but avoids the `core::fmt` machinery on a path that
/// runs once per (variant x transcript consequence).
fn format_u64(mut value: u64, buf: &mut [u8; 20]) -> &[u8] {
    if value == 0 {
        buf[0] = b'0';
        return &buf[..1];
    }
    let mut i = buf.len();
    while value > 0 {
        i -= 1;
        buf[i] = b'0' + (value % 10) as u8;
        value /= 10;
    }
    &buf[i..]
}

/// Write the `Extra` column of one default-format row: `KEY=VALUE` pairs joined
/// by `;`, keys in VEP's order (IMPACT, DISTANCE, STRAND, FLAGS, then the active
/// flag fields, then plugin keys alphabetically), values escaping only `;` as
/// `%3B`, lists joined by `,`, and a lone `-` when no field applies. `tc` is
/// `None` for an intergenic row, which carries IMPACT and no STRAND.
///
/// This is the hot path for the default `-o vep` format: it runs once per
/// (variant x transcript consequence). With no optional field group active and no
/// plugin data it writes the four fixed keys directly without allocating.
pub fn write_extra_fields(
    variant: &InputVariant,
    tc: Option<&vep_core::consequence::TranscriptConsequence>,
    plan: &vep_io::output::fields::ExtraFieldsPlan,
    writer: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    let has_plugin_data =
        !variant.plugin_data.is_empty() || tc.is_some_and(|t| !t.plugin_data.is_empty());
    if !plan.any_optional_group && !has_plugin_data && !variant.is_structural {
        return write_user_group(tc, writer);
    }

    let options = &plan.options;
    let mut first = true;
    for key in vep_io::output::fields::extra_keys_with(&plan.flag_fields, variant, tc, options) {
        let value = vep_io::output::fields::field_value(key, variant, tc, options, ",");
        if value.is_empty() {
            continue;
        }
        if !first {
            writer.write_all(b";")?;
        }
        first = false;
        writer.write_all(key.as_bytes())?;
        writer.write_all(b"=")?;
        if value.contains(';') {
            writer.write_all(value.replace(';', "%3B").as_bytes())?;
        } else {
            writer.write_all(value.as_bytes())?;
        }
    }
    if first {
        writer.write_all(b"-")?;
    }
    Ok(())
}

/// The four unconditional keys, written as bytes. IMPACT is always present;
/// DISTANCE only on an up/downstream row; STRAND and FLAGS only against a
/// transcript, FLAGS only when the transcript carries a `cds_` attribute.
fn write_user_group(
    tc: Option<&vep_core::consequence::TranscriptConsequence>,
    writer: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    let Some(tc) = tc else {
        return writer.write_all(b"IMPACT=MODIFIER");
    };
    writer.write_all(b"IMPACT=")?;
    writer.write_all(tc.impact.as_str().as_bytes())?;
    if let Some(dist) = tc.distance {
        let mut buf = [0u8; 20];
        writer.write_all(b";DISTANCE=")?;
        writer.write_all(format_u64(dist, &mut buf))?;
    }
    writer.write_all(b";STRAND=")?;
    writer.write_all(if tc.strand < 0 { b"-1" } else { b"1" })?;
    let mut wrote_flags = false;
    for flag in tc.flags.iter().filter(|f| f.starts_with("cds_")) {
        writer.write_all(if wrote_flags { b"," } else { b";FLAGS=" })?;
        writer.write_all(flag.as_bytes())?;
        wrote_flags = true;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_allele_frequency_single() {
        assert_eq!(extract_allele_frequency("T:0.003", "T"), Some(0.003));
        assert_eq!(extract_allele_frequency("T:0.003", "A"), None);
    }

    #[test]
    fn test_extract_allele_frequency_multiple() {
        assert_eq!(
            extract_allele_frequency("T:0.003,C:0.997", "T"),
            Some(0.003)
        );
        assert_eq!(
            extract_allele_frequency("T:0.003,C:0.997", "C"),
            Some(0.997)
        );
        assert_eq!(extract_allele_frequency("T:0.003,C:0.997", "G"), None);
    }

    #[test]
    fn test_extract_allele_frequency_zero() {
        assert_eq!(extract_allele_frequency("T:0", "T"), Some(0.0));
    }

    #[test]
    fn test_extract_allele_frequency_scientific() {
        assert_eq!(
            extract_allele_frequency("T:1.886e-05", "T"),
            Some(1.886e-05)
        );
    }

    #[test]
    fn test_find_colocated_basic() {
        let var = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());

        let mut freqs = HashMap::new();
        freqs.insert("AFR".to_string(), "G:0.003".to_string());

        let cached = vec![CachedVariation {
            variation_name: "rs123".to_string(),
            start: 100,
            end: 100,
            allele_string: "A/G".to_string(),
            strand: 1,
            frequencies: freqs,
            ..Default::default()
        }];

        let mut pos_idx = FxHashMap::default();
        pos_idx.insert(100u64, vec![0usize]);

        let colocated = find_colocated_variants(&var, &cached, &pos_idx);
        assert_eq!(colocated.len(), 1);
        assert_eq!(colocated[0].id, "rs123");
        assert_eq!(colocated[0].frequencies.get("AFR"), Some(&0.003));
    }

    #[test]
    fn test_find_colocated_skips_failed() {
        let var = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());

        let cached = vec![CachedVariation {
            variation_name: "rs_failed".to_string(),
            start: 100,
            end: 100,
            failed: 1,
            ..Default::default()
        }];

        let mut pos_idx = FxHashMap::default();
        pos_idx.insert(100u64, vec![0usize]);

        let colocated = find_colocated_variants(&var, &cached, &pos_idx);
        assert!(colocated.is_empty());
    }

    #[test]
    fn test_find_colocated_no_match() {
        let var = InputVariant::new("21".into(), 200, 200, b"A".to_vec(), b"G".to_vec());

        let cached = vec![CachedVariation {
            variation_name: "rs123".to_string(),
            start: 100,
            end: 100,
            ..Default::default()
        }];

        let mut pos_idx = FxHashMap::default();
        pos_idx.insert(100u64, vec![0usize]);

        let colocated = find_colocated_variants(&var, &cached, &pos_idx);
        assert!(colocated.is_empty());
    }

    // The `-o vep` path streams the Extra column straight into the writer, and
    // the end-to-end concordance comparison exercises only the no-flag workload,
    // so these tests pin the bytes of every branch: VEP's key order, the `;`
    // separator, the `-` placeholder, and the `;`-only escaping.

    use vep_core::consequence::{Impact, TranscriptConsequence};
    use vep_io::output::fields::{ExtraFieldsPlan, FieldOptions};

    /// Render `write_extra_fields` into a `String` for byte-level assertions.
    fn extra(
        variant: &InputVariant,
        tc: Option<&TranscriptConsequence>,
        options: &FieldOptions,
    ) -> String {
        let mut buf: Vec<u8> = Vec::new();
        let plan = ExtraFieldsPlan::new(options.clone());
        write_extra_fields(variant, tc, &plan, &mut buf).expect("write_extra_fields");
        String::from_utf8(buf).expect("Extra field is valid UTF-8")
    }

    #[test]
    fn test_format_u64_matches_display() {
        // format_u64 replaces `u64`'s Display on the per-consequence writer
        // path, so it must agree byte-for-byte over the interesting range:
        // zero, every single digit, every power-of-ten boundary (+/- 1), and
        // u64::MAX (the 20-digit case that sizes the buffer).
        let mut cases: Vec<u64> = (0..=1000).collect();
        let mut p: u64 = 1;
        for _ in 0..19 {
            p = p.saturating_mul(10);
            cases.extend_from_slice(&[p - 1, p, p + 1]);
        }
        cases.extend_from_slice(&[u64::MAX - 1, u64::MAX, 2500, 12_345_678_901_234_567_890]);

        for v in cases {
            let mut buf = [0u8; 20];
            let got = std::str::from_utf8(format_u64(v, &mut buf)).unwrap();
            assert_eq!(got, v.to_string(), "format_u64 mismatch for {v}");
        }
    }

    #[test]
    fn default_row_carries_impact_and_strand() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            impact: Impact::MODERATE,
            strand: -1,
            ..Default::default()
        };
        assert_eq!(
            extra(&variant, Some(&tc), &FieldOptions::default()),
            "IMPACT=MODERATE;STRAND=-1"
        );
    }

    #[test]
    fn distance_precedes_strand_and_flags_follow_it() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let tc = TranscriptConsequence {
            impact: Impact::MODIFIER,
            strand: 1,
            distance: Some(2500),
            flags: std::sync::Arc::from([
                "cds_start_NF".to_string(),
                "cds_end_NF".to_string(),
                "gencode_basic".to_string(),
            ]),
            ..Default::default()
        };
        assert_eq!(
            extra(&variant, Some(&tc), &FieldOptions::default()),
            "IMPACT=MODIFIER;DISTANCE=2500;STRAND=1;FLAGS=cds_start_NF,cds_end_NF"
        );
    }

    #[test]
    fn intergenic_row_carries_impact_alone() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        assert_eq!(
            extra(&variant, None, &FieldOptions::default()),
            "IMPACT=MODIFIER"
        );
    }

    #[test]
    fn optional_groups_follow_vep_order() {
        let mut variant = InputVariant::new("21".into(), 100, 101, b"AC".to_vec(), b"-".to_vec());
        let mut frequencies = std::collections::HashMap::new();
        frequencies.insert("AF".to_string(), 0.003);
        variant.colocated_variants.push(ColocatedVariant {
            id: "rs123".to_string(),
            clin_sig: vec!["pathogenic".to_string(), "likely_pathogenic".to_string()],
            somatic: true,
            phenotype: false,
            frequencies,
            ..Default::default()
        });
        let tc = TranscriptConsequence {
            impact: Impact::HIGH,
            strand: 1,
            hgvsc: Some("ENST00000366667.4:c.428C>T".to_string()),
            hgvsp: Some("ENSP00000355627.3:p.Ala143=".to_string()),
            sift: Some("deleterious(0.01)".to_string()),
            ..Default::default()
        };
        let options = FieldOptions {
            variant_class: true,
            hgvs: true,
            sift: true,
            af: true,
            check_existing: true,
            ..Default::default()
        };
        assert_eq!(
            extra(&variant, Some(&tc), &options),
            "IMPACT=HIGH;STRAND=1;VARIANT_CLASS=deletion;SIFT=deleterious(0.01);HGVSc=ENST00000366667.4:c.428C>T;HGVSp=ENSP00000355627.3:p.Ala143%3D;AF=0.003;CLIN_SIG=pathogenic,likely_pathogenic;SOMATIC=1;PHENO=0"
        );
    }

    #[test]
    fn plugin_keys_follow_flag_fields_alphabetically_and_escape_semicolons_only() {
        let variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let mut tc = TranscriptConsequence {
            impact: Impact::LOW,
            strand: 1,
            ..Default::default()
        };
        tc.plugin_data
            .insert("ZPlugin".to_string(), "a=b;c".to_string());
        tc.plugin_data
            .insert("APlugin".to_string(), "x".to_string());
        assert_eq!(
            extra(&variant, Some(&tc), &FieldOptions::default()),
            "IMPACT=LOW;STRAND=1;APlugin=x;ZPlugin=a=b%3Bc"
        );
    }
}
