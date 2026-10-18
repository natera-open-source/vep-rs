// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Build ExonCoordMapper pairs from GFF3 exon/CDS data.
//!
//! The mapper pairs define the cDNA <-> genomic coordinate mapping for each
//! transcript. The pairs
//! must exactly match Perl's TranscriptMapper output.

use crate::gff3::{CdsRecord, ExonRecord};

/// A computed mapper pair for the JSON cache.
#[derive(Debug, Clone)]
pub struct MapperPair {
    pub from_start: u64,
    pub from_end: u64,
    pub to_start: u64,
    pub to_end: u64,
    pub ori: i8,
}

/// Result of building the mapper for a transcript.
#[derive(Debug)]
pub struct MapperResult {
    pub pairs: Vec<MapperPair>,
    pub cdna_coding_start: Option<u64>,
    pub cdna_coding_end: Option<u64>,
    /// Exon 1's Ensembl phase, the value `TranscriptMapper.pm:110` records.
    pub start_phase: i8,
}

/// Build ExonCoordMapper pairs for a transcript.
///
/// For each transcript:
/// 1. Sort exons in 5'->3' order (ascending for forward, descending for reverse)
/// 2. Accumulate cDNA positions while traversing exons
/// 3. Each exon becomes a mapper pair: from (cDNA) -> to (genomic)
/// 4. Compute cdna_coding_start/end by mapping the first/last CDS positions
///
/// `first_exon_phase` is the Ensembl phase of exon 1, recorded as the mapper's
/// `start_phase` (`TranscriptMapper.pm:110`); a transcript without exons records -1.
pub fn build_mapper(
    exons: &[ExonRecord],
    cds_regions: &[CdsRecord],
    strand: i8,
    first_exon_phase: i8,
) -> MapperResult {
    if exons.is_empty() {
        return MapperResult {
            pairs: Vec::new(),
            cdna_coding_start: None,
            cdna_coding_end: None,
            start_phase: -1,
        };
    }

    let mut sorted_exons: Vec<&ExonRecord> = exons.iter().collect();
    if strand >= 0 {
        sorted_exons.sort_by_key(|e| e.start);
    } else {
        sorted_exons.sort_by_key(|e| std::cmp::Reverse(e.start));
    }

    let ori = if strand >= 0 { 1i8 } else { -1i8 };

    let mut pairs = Vec::with_capacity(sorted_exons.len());
    let mut cdna_pos: u64 = 1;

    for exon in &sorted_exons {
        let exon_len = exon.end - exon.start + 1;
        let from_start = cdna_pos;
        let from_end = cdna_pos + exon_len - 1;

        pairs.push(MapperPair {
            from_start,
            from_end,
            to_start: exon.start,
            to_end: exon.end,
            ori,
        });

        cdna_pos = from_end + 1;
    }

    let (cdna_coding_start, cdna_coding_end) = compute_cdna_coding(&pairs, cds_regions, strand);

    MapperResult {
        pairs,
        cdna_coding_start,
        cdna_coding_end,
        start_phase: first_exon_phase,
    }
}

/// Compute cdna_coding_start and cdna_coding_end from CDS regions.
///
/// The CDS defines the coding portion of the transcript; the first and last CDS
/// bases are located in cDNA coordinates. Both are `None` without a CDS.
fn compute_cdna_coding(
    pairs: &[MapperPair],
    cds_regions: &[CdsRecord],
    strand: i8,
) -> (Option<u64>, Option<u64>) {
    let (Some(cds_genomic_start), Some(cds_genomic_end)) = (
        cds_regions.iter().map(|c| c.start).min(),
        cds_regions.iter().map(|c| c.end).max(),
    ) else {
        return (None, None);
    };

    let (first_cds_pos, last_cds_pos) = if strand >= 0 {
        (cds_genomic_start, cds_genomic_end)
    } else {
        (cds_genomic_end, cds_genomic_start)
    };

    (
        genomic_to_cdna(pairs, first_cds_pos, strand),
        genomic_to_cdna(pairs, last_cds_pos, strand),
    )
}

/// Map a genomic position to its cDNA coordinate using the mapper pairs.
fn genomic_to_cdna(pairs: &[MapperPair], genomic_pos: u64, strand: i8) -> Option<u64> {
    for pair in pairs {
        if genomic_pos >= pair.to_start && genomic_pos <= pair.to_end {
            let offset = if strand >= 0 {
                genomic_pos - pair.to_start
            } else {
                pair.to_end - genomic_pos
            };
            return Some(pair.from_start + offset);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gff3::{CdsRecord, ExonRecord};

    fn make_exon(start: u64, end: u64, strand: i8) -> ExonRecord {
        ExonRecord {
            exon_id: None,
            start,
            end,
            strand,
            rank: None,
            ensembl_phase: None,
            ensembl_end_phase: None,
        }
    }

    fn make_cds(start: u64, end: u64, strand: i8, phase: i8) -> CdsRecord {
        CdsRecord {
            start,
            end,
            strand,
            phase,
            protein_id: None,
            protein_version: None,
        }
    }

    #[test]
    fn test_forward_strand_two_exons() {
        let exons = vec![make_exon(100, 200, 1), make_exon(300, 500, 1)];
        let cds = vec![];
        let result = build_mapper(&exons, &cds, 1, -1);

        assert_eq!(result.pairs.len(), 2);

        // Exon 1: genomic 100-200 -> cDNA 1-101
        assert_eq!(result.pairs[0].from_start, 1);
        assert_eq!(result.pairs[0].from_end, 101);
        assert_eq!(result.pairs[0].to_start, 100);
        assert_eq!(result.pairs[0].to_end, 200);
        assert_eq!(result.pairs[0].ori, 1);

        // Exon 2: genomic 300-500 -> cDNA 102-302
        assert_eq!(result.pairs[1].from_start, 102);
        assert_eq!(result.pairs[1].from_end, 302);
        assert_eq!(result.pairs[1].to_start, 300);
        assert_eq!(result.pairs[1].to_end, 500);
        assert_eq!(result.pairs[1].ori, 1);
    }

    #[test]
    fn test_reverse_strand_two_exons() {
        let exons = vec![make_exon(100, 200, -1), make_exon(300, 500, -1)];
        let cds = vec![];
        let result = build_mapper(&exons, &cds, -1, -1);

        assert_eq!(result.pairs.len(), 2);

        // Reverse strand: 5'->3' means descending position.
        // Exon at 300-500 comes first in cDNA.
        assert_eq!(result.pairs[0].from_start, 1);
        assert_eq!(result.pairs[0].from_end, 201);
        assert_eq!(result.pairs[0].to_start, 300);
        assert_eq!(result.pairs[0].to_end, 500);
        assert_eq!(result.pairs[0].ori, -1);

        // Exon at 100-200 is second in cDNA.
        assert_eq!(result.pairs[1].from_start, 202);
        assert_eq!(result.pairs[1].from_end, 302);
        assert_eq!(result.pairs[1].to_start, 100);
        assert_eq!(result.pairs[1].to_end, 200);
        assert_eq!(result.pairs[1].ori, -1);
    }

    #[test]
    fn test_forward_strand_with_cds() {
        let exons = vec![make_exon(100, 200, 1), make_exon(300, 500, 1)];
        let cds = vec![make_cds(150, 200, 1, 0), make_cds(300, 480, 1, 2)];
        let result = build_mapper(&exons, &cds, 1, -1);

        // CDS starts at genomic 150, which is in exon 1 (100-200).
        // Offset = 150 - 100 = 50, so cDNA coding start = 1 + 50 = 51.
        assert_eq!(result.cdna_coding_start, Some(51));

        // CDS ends at genomic 480, which is in exon 2 (300-500).
        // Offset = 480 - 300 = 180, so cDNA coding end = 102 + 180 = 282.
        assert_eq!(result.cdna_coding_end, Some(282));
    }

    #[test]
    fn test_reverse_strand_with_cds() {
        // Reverse strand transcript: exons at 100-200 and 300-500.
        // CDS from 120-200 and 300-450.
        let exons = vec![make_exon(100, 200, -1), make_exon(300, 500, -1)];
        let cds = vec![make_cds(120, 200, -1, 1), make_cds(300, 450, -1, 0)];
        let result = build_mapper(&exons, &cds, -1, -1);

        // Reverse strand: first CDS in 5'->3' is the one with highest genomic end.
        // CDS at 300-450: genomic end 450 is the 5' start.
        // In mapper pairs: pair[0] = 300-500 -> cDNA 1-201
        // Offset of genomic 450 from end: 500 - 450 = 50, cDNA = 1 + 50 = 51
        assert_eq!(result.cdna_coding_start, Some(51));

        // Last CDS in 5'->3' is CDS at 120-200: genomic start 120 is the 3' end.
        // In mapper pairs: pair[1] = 100-200 -> cDNA 202-302
        // Offset of genomic 120 from end: 200 - 120 = 80, cDNA = 202 + 80 = 282
        assert_eq!(result.cdna_coding_end, Some(282));
    }

    #[test]
    fn test_empty_exons() {
        let result = build_mapper(&[], &[], 1, 2);
        assert!(result.pairs.is_empty());
        assert!(result.cdna_coding_start.is_none());
        assert!(result.cdna_coding_end.is_none());
        assert_eq!(result.start_phase, -1);
    }

    /// The mapper's `start_phase` is exon 1's Ensembl phase (`TranscriptMapper.pm:110`),
    /// here -1 for a transcript whose first exon is untranslated, not the GFF3 phase
    /// column of the 5'-most CDS segment (2 here) and not the translation start exon's
    /// phase; a `cds_start_NF` transcript whose exon 1 begins mid-codon records that
    /// exon's phase, 1.
    #[test]
    fn start_phase_is_ensembl_phase() {
        let exons = vec![make_exon(100, 200, 1), make_exon(300, 500, 1)];
        let cds = vec![make_cds(300, 480, 1, 2)];
        let result = build_mapper(&exons, &cds, 1, -1);
        assert_eq!(result.start_phase, -1);
        assert_eq!(result.cdna_coding_start, Some(102));

        let partial = build_mapper(&exons, &[make_cds(100, 200, 1, 2)], 1, 1);
        assert_eq!(partial.start_phase, 1);
    }
}
