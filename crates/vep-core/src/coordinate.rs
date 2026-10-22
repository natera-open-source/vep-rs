// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Genomic coordinate types: strand, position, range, and chromosome normalization.
//!
//! All coordinates are 1-based inclusive, matching VEP's convention.
//! Provides [`normalize_chromosome`] for converting UCSC-style `chr` prefixes
//! to Ensembl-style names at input time.

use std::fmt;

/// The largest coordinate an input reader accepts, 2^62 - 1: above every chromosome,
/// and below `i64::MAX` by enough that a coordinate plus a span, a flanking base or a
/// distance stays representable in `u64` and in the `i64` the effects layer converts
/// to.
pub const MAX_COORD: u64 = (1 << 62) - 1;

/// Reads a coordinate column: a non-negative integer no larger than [`MAX_COORD`],
/// `None` for anything else.
pub fn parse_coordinate(text: &str) -> Option<u64> {
    text.parse::<u64>().ok().filter(|&n| n <= MAX_COORD)
}

/// Reads the POS column as Ensembl VEP numifies it (ensembl-io BaseVCF4.pm:263,
/// `$start++`): a coordinate in any spelling that names an exact integer, so `3e7`,
/// `30425949.0`, `30425949.` and ` 30425949` are coordinates where `30425949.5`,
/// `30_425_949`, `0x1D04A1D`, `.` and `-5` are not.
pub fn parse_pos(text: &str) -> Option<u64> {
    parse_coordinate(text).or_else(|| {
        let trimmed = text.trim_ascii();
        if trimmed.starts_with('-') {
            return None;
        }
        let value: f64 = trimmed.parse().ok()?;
        // `MAX_COORD as f64` rounds up to 2^62, so the bound is strict.
        (value.is_finite() && value.fract() == 0.0 && value < MAX_COORD as f64)
            .then_some(value as u64)
    })
}

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

/// Whether `name` can name a chromosome's directory or file of a cache: one
/// ordinary path component, so a path joined from it stays inside its parent. An
/// empty name, `.`, `..`, an absolute path and a name holding a separator of
/// either convention are refused.
pub fn is_cache_path_component(name: &str) -> bool {
    let mut components = std::path::Path::new(name).components();
    matches!(
        (components.next(), components.next()),
        (Some(std::path::Component::Normal(_)), None)
    ) && !name.contains(['/', '\\'])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A chromosome name is one ordinary path component: the parent and current
    /// directories, an empty name, a rooted path and a name with a separator of
    /// either convention would place a cache file outside its chromosome's directory.
    #[test]
    fn a_cache_path_component_is_one_ordinary_name() {
        for name in ["21", "X", "MT", "LRG_485", "HSCHR6_MHC_COX", "chr21", "a.b"] {
            assert!(is_cache_path_component(name), "{name:?}");
        }
        for name in [
            "", ".", "..", "../21", "21/22", "/tmp", "21\\22", "..\\21", "./21",
        ] {
            assert!(!is_cache_path_component(name), "{name:?}");
        }
    }

    /// `MAX_COORD` is read and one past it is not; a negative number, `-0` included,
    /// and a non-integer are not coordinates.
    #[test]
    fn parse_coordinate_reads_up_to_the_ceiling() {
        assert_eq!(parse_coordinate("0"), Some(0));
        assert_eq!(parse_coordinate(&MAX_COORD.to_string()), Some(MAX_COORD));
        assert_eq!(parse_coordinate(&(MAX_COORD + 1).to_string()), None);
        assert_eq!(parse_coordinate(&u64::MAX.to_string()), None);
        for text in ["-1", "-0", "1e3", "", "."] {
            assert_eq!(parse_coordinate(text), None, "{text:?}");
        }
    }

    /// Ensembl VEP 116.2 reads `21 3e7 x N <DEL> . . SVTYPE=DEL;END=30425949` as the
    /// deletion `21:30000001-30425949` and `21 30425949.0 x A G` at 30425949, numifying
    /// the POS text; `30_425_949` numifies to 30 and `0x1D04A1D` and `.` to 0, so those
    /// spellings name no coordinate here, nor does a fraction or a negative number.
    #[test]
    fn parse_pos_reads_an_exact_integer_in_any_spelling() {
        for (text, pos) in [
            ("3e7", 30_000_000),
            ("30425949.0", 30_425_949),
            ("30425949.", 30_425_949),
            ("3.0425949E7", 30_425_949),
            ("1e3", 1000),
            (" 30425949", 30_425_949),
            ("30425949 ", 30_425_949),
            ("+30425949", 30_425_949),
            ("030425949", 30_425_949),
            ("0", 0),
        ] {
            assert_eq!(parse_pos(text), Some(pos), "{text:?}");
        }
        for text in [
            "30425949.5",
            "30_425_949",
            "0x1D04A1D",
            "\u{ff13}\u{ff10}",
            ".",
            "",
            "-5",
            "-0.0",
            "NaN",
            "Inf",
            "1e400",
            "4611686018427387904",
            "4611686018427387904.0",
            "4.611686018427388e18",
        ] {
            assert_eq!(parse_pos(text), None, "{text:?}");
        }
    }

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
