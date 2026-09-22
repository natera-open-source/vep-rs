// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Genomic coordinate types: strand, position, range, and chromosome normalization.
//!
//! All coordinates are 1-based inclusive, matching VEP's convention.
//! Provides [`normalize_chromosome`] for converting UCSC-style `chr` prefixes
//! to Ensembl-style names at input time.

use std::fmt;

/// Strand orientation of a genomic feature (forward/+1 or reverse/-1).
///
/// Affects how upstream/downstream distances are computed and how alternate
/// alleles are reverse-complemented for coding effect analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    /// Numeric representation matching VEP convention (1 or -1).
    pub fn as_i8(self) -> i8 {
        match self {
            Strand::Forward => 1,
            Strand::Reverse => -1,
        }
    }

    pub fn from_i8(val: i8) -> Option<Self> {
        match val {
            1 => Some(Strand::Forward),
            -1 => Some(Strand::Reverse),
            _ => None,
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_i8())
    }
}

/// A position on a chromosome, 1-based.
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub struct GenomicPosition {
    pub chr: String,
    pub pos: u64,
    pub strand: Strand,
}

impl fmt::Display for GenomicPosition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.chr, self.pos)
    }
}

/// A range on a chromosome, 1-based inclusive on both ends.
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub struct GenomicRange {
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

impl GenomicRange {
    pub fn new(chr: String, start: u64, end: u64, strand: Strand) -> Self {
        Self {
            chr,
            start,
            end,
            strand,
        }
    }

    /// Length of the range in base pairs.
    pub fn len(&self) -> u64 {
        if self.end >= self.start {
            self.end - self.start + 1
        } else {
            0
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Check if two ranges on the same chromosome overlap.
    pub fn overlaps(&self, other: &GenomicRange) -> bool {
        self.chr == other.chr && self.start <= other.end && other.start <= self.end
    }

    /// Compute the cache region index for a given position.
    /// Region index `s` covers coordinates `[s * region_size + 1, (s+1) * region_size]`.
    pub fn cache_region_index(pos: u64, region_size: u64) -> u64 {
        if pos == 0 {
            0
        } else {
            (pos - 1) / region_size
        }
    }
}

impl fmt::Display for GenomicRange {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.start == self.end {
            write!(f, "{}:{}", self.chr, self.start)
        } else {
            write!(f, "{}:{}-{}", self.chr, self.start, self.end)
        }
    }
}

/// Normalize a chromosome name to match Ensembl cache convention.
///
/// Strips the `chr` prefix (case-insensitive) used by UCSC/GRCh38 VCFs so that
/// chromosome names match the Ensembl cache (`1`, `2`, ..., `X`, `Y`, `MT`).
/// Also converts `chrM` to `MT` for mitochondrial consistency.
///
/// This is applied at input time so all downstream matching (transcript lookup,
/// variation lookup, overlap detection) works with a single representation.
pub fn normalize_chromosome(chr: &str) -> String {
    let stripped = chr
        .strip_prefix("chr")
        .or_else(|| chr.strip_prefix("Chr"))
        .or_else(|| chr.strip_prefix("CHR"))
        .unwrap_or(chr);

    if stripped.eq_ignore_ascii_case("M") {
        "MT".to_string()
    } else {
        stripped.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_range_display() {
        let snp = GenomicRange::new("21".into(), 25000100, 25000100, Strand::Forward);
        assert_eq!(snp.to_string(), "21:25000100");

        let del = GenomicRange::new("21".into(), 25000100, 25000105, Strand::Forward);
        assert_eq!(del.to_string(), "21:25000100-25000105");
    }

    #[test]
    fn test_range_overlap() {
        let a = GenomicRange::new("21".into(), 100, 200, Strand::Forward);
        let b = GenomicRange::new("21".into(), 150, 250, Strand::Forward);
        let c = GenomicRange::new("21".into(), 201, 300, Strand::Forward);
        let d = GenomicRange::new("22".into(), 100, 200, Strand::Forward);

        assert!(a.overlaps(&b));
        assert!(!a.overlaps(&c));
        assert!(!a.overlaps(&d));
    }

    #[test]
    fn test_cache_region_index() {
        assert_eq!(GenomicRange::cache_region_index(1, 1_000_000), 0);
        assert_eq!(GenomicRange::cache_region_index(1_000_000, 1_000_000), 0);
        assert_eq!(GenomicRange::cache_region_index(1_000_001, 1_000_000), 1);
        assert_eq!(GenomicRange::cache_region_index(25_000_001, 1_000_000), 25);
    }

    #[test]
    fn test_normalize_chromosome_no_prefix() {
        assert_eq!(normalize_chromosome("21"), "21");
        assert_eq!(normalize_chromosome("X"), "X");
        assert_eq!(normalize_chromosome("MT"), "MT");
    }

    #[test]
    fn test_normalize_chromosome_strip_chr() {
        assert_eq!(normalize_chromosome("chr21"), "21");
        assert_eq!(normalize_chromosome("chrX"), "X");
        assert_eq!(normalize_chromosome("Chr1"), "1");
        assert_eq!(normalize_chromosome("CHR22"), "22");
    }

    #[test]
    fn test_normalize_chromosome_mito() {
        assert_eq!(normalize_chromosome("chrM"), "MT");
        assert_eq!(normalize_chromosome("M"), "MT");
        assert_eq!(normalize_chromosome("MT"), "MT");
    }
}
