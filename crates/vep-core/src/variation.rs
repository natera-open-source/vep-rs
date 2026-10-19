// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Variation records of the Ensembl VEP variation cache.

use std::collections::HashMap;

/// One known variant read from the variation cache, the hash Ensembl VEP's
/// `BaseCacheVariation.pm::parse_variation` builds at release/116.2.
///
/// The named fields are the cache's fixed columns; every other column of the
/// cache's `variation_cols` is a population frequency and lands in
/// [`frequencies`](Self::frequencies) under its own name, so the record carries
/// whatever populations the cache declares.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CachedVariation {
    /// The chromosome column of the start-indexed layout, the file's own
    /// chromosome; `None` on the per-region layout, which has no such column.
    pub chr: Option<String>,
    /// Variation name / identifier (e.g., "rs699").
    pub variation_name: String,
    /// Failed status: 0 = not failed.
    pub failed: u8,
    /// Somatic flag: 0 = germline, 1 = somatic.
    pub somatic: u8,
    /// 1-based start position.
    pub start: u64,
    /// 1-based end position (inclusive); an insertion has `end = start - 1`.
    pub end: u64,
    /// Allele string (e.g., "A/G", "C/A/T"); a null allele has no `/`.
    pub allele_string: String,
    /// Strand: 1 or -1.
    pub strand: i8,
    /// Clinical significance terms (e.g., "benign,likely_benign").
    pub clin_sig: Option<String>,
    /// Phenotype or disease association flag: 0 or 1.
    pub phenotype_or_disease: u8,
    /// The ClinVar allele and somatic fields, present on the few records that
    /// carry any of them; boxed so the records without them pay one pointer.
    pub clinvar: Option<Box<ClinVarDetail>>,
    /// PubMed IDs (comma-separated).
    pub pubmed: Option<String>,
    /// Variant synonyms (encoded string from cache).
    pub var_synonyms: Option<String>,
    /// Population allele frequencies.
    /// Key = population column name, value = "allele:freq,allele:freq" encoded string.
    pub frequencies: HashMap<String, String>,
}

impl CachedVariation {
    /// [`ClinVarDetail::clin_sig_allele_entries`] of the record's ClinVar
    /// fields; empty when it carries none.
    pub fn clin_sig_allele_entries(&self) -> Vec<(&str, &str)> {
        self.clinvar
            .as_deref()
            .map(ClinVarDetail::clin_sig_allele_entries)
            .unwrap_or_default()
    }
}

/// The ClinVar fields of a variation record beyond `clin_sig`.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct ClinVarDetail {
    /// Significance per risk allele, `allele:term` pairs joined by `;`.
    pub clin_sig_allele: Option<String>,
    /// The reference allele behind `clin_sig_allele`.
    pub clin_sig_ref_allele: Option<String>,
    /// Somatic classification, `phenotype(classification)` entries joined by `;`.
    pub clinical_impact: Option<String>,
}

impl ClinVarDetail {
    /// The `(allele, term)` pairs of `clin_sig_allele`, in the cache's order.
    ///
    /// An entry is `allele:term`; the allele is the text before its first `:`
    /// and the term the text between the first and the second (`OutputFactory.pm:1106-1108`
    /// at release/116.2), so an entry with no `:` has an empty term.
    pub fn clin_sig_allele_entries(&self) -> Vec<(&str, &str)> {
        let Some(entries) = self.clin_sig_allele.as_deref() else {
            return Vec::new();
        };
        entries
            .split(';')
            .map(|entry| {
                let (allele, rest) = entry.split_once(':').unwrap_or((entry, ""));
                let term = rest.split(':').next().unwrap_or("");
                (allele, term)
            })
            .collect()
    }
}

/// `reverse_comp` of ensembl's `Utils/Sequence.pm`: the sequence reversed and
/// complemented over `ACGT`, the IUPAC codes `RYMKSWHBVDNX` and their lower-case
/// forms, every other character kept. `S`, `W`, `N` and `X` are their own
/// complement, so they need no arm of their own.
pub fn reverse_comp(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            't' => 'a',
            'r' => 'y',
            'y' => 'r',
            'm' => 'k',
            'k' => 'm',
            'h' => 'd',
            'b' => 'v',
            'v' => 'b',
            'd' => 'h',
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'R' => 'Y',
            'Y' => 'R',
            'M' => 'K',
            'K' => 'M',
            'H' => 'D',
            'B' => 'V',
            'V' => 'B',
            'D' => 'H',
            other => other,
        })
        .collect()
}

/// A string's value as a Perl numeric context reads it: its longest numeric
/// prefix, 0 when it has none.
pub fn perl_numeric(text: &str) -> f64 {
    let text = text.trim_start();
    (1..=text.len())
        .rev()
        .filter(|&len| text.is_char_boundary(len))
        .find_map(|len| text[..len].parse::<f64>().ok())
        .unwrap_or(0.0)
}

/// A number as Perl stringifies it (`sv_2pv_flags` in `sv.c`): `0` for a zero of
/// either sign, `Inf`, `-Inf` or `NaN` for a non-finite value (a cache value such
/// as `1e999` reads as infinite), and `%.15g` otherwise.
pub fn perl_number_string(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value == f64::INFINITY {
        return "Inf".to_string();
    }
    if value == f64::NEG_INFINITY {
        return "-Inf".to_string();
    }
    let scientific = format!("{value:.14e}");
    let (mantissa, exponent) = scientific.split_once('e').expect("exponent form");
    let exponent: i32 = exponent.parse().expect("exponent digits");
    if !(-4..15).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!(
            "{}e{sign}{:02}",
            trim_fraction_zeros(mantissa),
            exponent.abs()
        )
    } else {
        let decimals = (14 - exponent).max(0) as usize;
        trim_fraction_zeros(&format!("{value:.decimals$}"))
    }
}

fn trim_fraction_zeros(text: &str) -> String {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        text.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perl_number_string_prints_fifteen_significant_digits() {
        assert_eq!(perl_number_string(1.0 - 0.7), "0.3");
        assert_eq!(perl_number_string(1.0 - 0.9999), "9.9999999999989e-05");
        assert_eq!(perl_number_string(0.0), "0");
        assert_eq!(perl_number_string(1.0), "1");
        assert_eq!(perl_number_string(0.25), "0.25");
        assert_eq!(perl_number_string(1e-05), "1e-05");
        assert_eq!(perl_number_string(1e15), "1e+15");
        assert_eq!(perl_number_string(123456789012345.0), "123456789012345");
        assert_eq!(perl_number_string(1.0 - perl_numeric("1e999")), "-Inf");
        assert_eq!(perl_number_string(f64::INFINITY), "Inf");
        assert_eq!(perl_number_string(f64::NAN), "NaN");
        assert_eq!(perl_numeric("0.7abc"), 0.7);
        assert_eq!(perl_numeric("abc"), 0.0);
        assert_eq!(perl_numeric("1e-05"), 1e-05);
    }

    #[test]
    fn test_cached_variation_default() {
        let v = CachedVariation::default();
        assert!(v.chr.is_none());
        assert_eq!(v.variation_name, "");
        assert_eq!(v.failed, 0);
        assert_eq!(v.somatic, 0);
        assert_eq!(v.start, 0);
        assert_eq!(v.end, 0);
        assert_eq!(v.allele_string, "");
        assert_eq!(v.strand, 0);
        assert!(v.clin_sig.is_none());
        assert_eq!(v.phenotype_or_disease, 0);
        assert!(v.clinvar.is_none());
        assert!(v.pubmed.is_none());
        assert!(v.frequencies.is_empty());
        assert!(v.var_synonyms.is_none());
        assert!(v.clin_sig_allele_entries().is_empty());
    }

    /// The `rs142566111` fixture is that record's `clin_sig_allele` in the
    /// release-116 GRCh38 cache.
    #[test]
    fn clin_sig_allele_entries_split_on_semicolon_then_colon() {
        let v = CachedVariation {
            variation_name: "rs142566111".to_string(),
            allele_string: "T/A/G".to_string(),
            clinvar: Some(Box::new(ClinVarDetail {
                clin_sig_allele: Some("G:likely_benign;A:likely_benign".to_string()),
                ..Default::default()
            })),
            ..Default::default()
        };
        assert_eq!(
            v.clin_sig_allele_entries(),
            vec![("G", "likely_benign"), ("A", "likely_benign")]
        );
        let odd = CachedVariation {
            clinvar: Some(Box::new(ClinVarDetail {
                clin_sig_allele: Some(":pathogenic;C:a:b;G".to_string()),
                ..Default::default()
            })),
            ..Default::default()
        };
        assert_eq!(
            odd.clin_sig_allele_entries(),
            vec![("", "pathogenic"), ("C", "a"), ("G", "")]
        );
    }

    /// `reverse_comp` of ensembl's `Utils/Sequence.pm` over the bases, the IUPAC
    /// codes in both cases and characters outside the table, which are kept and
    /// reversed: `ACGTN-` is `-NACGT`, `RYMKSWHBVDX` is `XHBVDWSMKRY`, lower case
    /// likewise, and `COSMIC_MUTATION` is `NOIATAUK_GIKSOG`.
    #[test]
    fn reverse_comp_matches_ensembl() {
        assert_eq!(reverse_comp("ACGTN-"), "-NACGT");
        assert_eq!(reverse_comp("acgtn"), "nacgt");
        assert_eq!(reverse_comp("RYMKSWHBVDX"), "XHBVDWSMKRY");
        assert_eq!(reverse_comp("rymkswhbvdx"), "xhbvdwsmkry");
        assert_eq!(reverse_comp("COSMIC_MUTATION"), "NOIATAUK_GIKSOG");
        assert_eq!(reverse_comp("GAT-"), "-ATC");
    }
}
