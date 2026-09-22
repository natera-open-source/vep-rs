// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Allele manipulation utilities: trimming, normalization, and Ts/Tv classification.
//!
//! [`trim_alleles`] converts VCF-style allele representations (with shared prefix/suffix)
//! to VEP-style minimized alleles (e.g., `"ACGT"/"ACG"` becomes `"T"/"-"` with a
//! coordinate offset). Also provides IUPAC ambiguity code lookup and transition/
//! transversion classification.

/// Transition/Transversion classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TsTv {
    Transition,
    Transversion,
}

impl TsTv {
    pub fn short_name(self) -> &'static str {
        match self {
            TsTv::Transition => "Ts",
            TsTv::Transversion => "Tv",
        }
    }
}

/// Classify an SNV as transition or transversion.
/// Returns None if not a valid single-base substitution.
pub fn classify_ts_tv(ref_base: u8, alt_base: u8) -> Option<TsTv> {
    let r = ref_base.to_ascii_uppercase();
    let a = alt_base.to_ascii_uppercase();
    match (r, a) {
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') => Some(TsTv::Transition),
        (b'A', b'C')
        | (b'C', b'A')
        | (b'G', b'T')
        | (b'T', b'G')
        | (b'A', b'T')
        | (b'T', b'A')
        | (b'C', b'G')
        | (b'G', b'C') => Some(TsTv::Transversion),
        _ => None,
    }
}

/// IUPAC ambiguity codes for pairs of bases.
pub fn ambiguity_code(bases: &[u8]) -> Option<u8> {
    let mut sorted: Vec<u8> = bases.iter().map(|b| b.to_ascii_uppercase()).collect();
    sorted.sort_unstable();
    sorted.dedup();
    match sorted.as_slice() {
        [b'A'] => Some(b'A'),
        [b'C'] => Some(b'C'),
        [b'G'] => Some(b'G'),
        [b'T'] => Some(b'T'),
        [b'A', b'G'] => Some(b'R'),
        [b'C', b'T'] => Some(b'Y'),
        [b'A', b'C'] => Some(b'M'),
        [b'G', b'T'] => Some(b'K'),
        [b'C', b'G'] => Some(b'S'),
        [b'A', b'T'] => Some(b'W'),
        [b'A', b'C', b'G'] => Some(b'V'),
        [b'A', b'C', b'T'] => Some(b'H'),
        [b'A', b'G', b'T'] => Some(b'D'),
        [b'C', b'G', b'T'] => Some(b'B'),
        [b'A', b'C', b'G', b'T'] => Some(b'N'),
        _ => None,
    }
}

/// Trim shared prefix and suffix from ref and alt alleles.
///
/// Returns (trimmed_ref, trimmed_alt, bases_trimmed_from_start).
/// Converts VCF-style representation (e.g., "ACGT" -> "ACG" becomes "T" -> "-")
/// to VEP-style ("T/-" with adjusted coordinates).
pub fn trim_alleles(ref_allele: &[u8], alt_allele: &[u8]) -> (Vec<u8>, Vec<u8>, usize) {
    // An index rather than `Vec::remove(0)`, which makes the trim O(n^2).
    let prefix_len = ref_allele
        .iter()
        .zip(alt_allele.iter())
        .take_while(|(r, a)| r.eq_ignore_ascii_case(a))
        .count();
    let mut ref_a = &ref_allele[prefix_len..];
    let mut alt_a = &alt_allele[prefix_len..];

    while !ref_a.is_empty()
        && !alt_a.is_empty()
        && ref_a
            .last()
            .unwrap()
            .eq_ignore_ascii_case(alt_a.last().unwrap())
    {
        ref_a = &ref_a[..ref_a.len() - 1];
        alt_a = &alt_a[..alt_a.len() - 1];
    }

    // VEP uses "-" for empty alleles
    let trimmed_ref = if ref_a.is_empty() {
        b"-".to_vec()
    } else {
        ref_a.to_vec()
    };
    let trimmed_alt = if alt_a.is_empty() {
        b"-".to_vec()
    } else {
        alt_a.to_vec()
    };

    (trimmed_ref, trimmed_alt, prefix_len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ts_tv() {
        assert_eq!(classify_ts_tv(b'A', b'G'), Some(TsTv::Transition));
        assert_eq!(classify_ts_tv(b'C', b'T'), Some(TsTv::Transition));
        assert_eq!(classify_ts_tv(b'A', b'C'), Some(TsTv::Transversion));
        assert_eq!(classify_ts_tv(b'G', b'T'), Some(TsTv::Transversion));
        assert_eq!(classify_ts_tv(b'A', b'N'), None);
    }

    #[test]
    fn test_trim_alleles_snv() {
        let (r, a, trim) = trim_alleles(b"A", b"G");
        assert_eq!(r, b"A");
        assert_eq!(a, b"G");
        assert_eq!(trim, 0);
    }

    #[test]
    fn test_trim_alleles_insertion() {
        // VCF: pos=100, ref=A, alt=ATT -> VEP: pos=100-101, ref="-", alt="TT"
        let (r, a, trim) = trim_alleles(b"A", b"ATT");
        assert_eq!(r, b"-");
        assert_eq!(a, b"TT");
        assert_eq!(trim, 1);
    }

    #[test]
    fn test_trim_alleles_deletion() {
        // VCF: pos=100, ref=ATT, alt=A -> VEP: pos=101-102, ref="TT", alt="-"
        let (r, a, trim) = trim_alleles(b"ATT", b"A");
        assert_eq!(r, b"TT");
        assert_eq!(a, b"-");
        assert_eq!(trim, 1);
    }

    #[test]
    fn test_trim_alleles_complex() {
        let (r, a, trim) = trim_alleles(b"ACGT", b"ATGT");
        assert_eq!(r, b"C");
        assert_eq!(a, b"T");
        assert_eq!(trim, 1);
    }

    #[test]
    fn test_ambiguity_code() {
        assert_eq!(ambiguity_code(b"AG"), Some(b'R'));
        assert_eq!(ambiguity_code(b"CT"), Some(b'Y'));
        assert_eq!(ambiguity_code(b"ACGT"), Some(b'N'));
        assert_eq!(ambiguity_code(b"A"), Some(b'A'));
    }
}
