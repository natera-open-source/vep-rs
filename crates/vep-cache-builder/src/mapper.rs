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
    pub start_phase: i8,
}

/// Build ExonCoordMapper pairs for a transcript.
///
/// For each transcript:
/// 1. Sort exons in 5'->3' order (ascending for forward, descending for reverse)
/// 2. Accumulate cDNA positions while traversing exons
/// 3. Each exon becomes a mapper pair: from (cDNA) -> to (genomic)
/// 4. Compute cdna_coding_start/end by mapping the first/last CDS positions
pub fn build_mapper(exons: &[ExonRecord], cds_regions: &[CdsRecord], strand: i8) -> MapperResult {
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

    let (cdna_coding_start, cdna_coding_end, start_phase) = if cds_regions.is_empty() {
        (None, None, -1i8)
    } else {
        compute_cdna_coding(&pairs, cds_regions, strand)
    };

    MapperResult {
        pairs,
        cdna_coding_start,
        cdna_coding_end,
        start_phase,
    }
}

/// Compute cdna_coding_start and cdna_coding_end from CDS regions.
///
/// The CDS defines the coding portion of the transcript; the first and last CDS
/// bases are located in cDNA coordinates.
fn compute_cdna_coding(
    pairs: &[MapperPair],
    cds_regions: &[CdsRecord],
    strand: i8,
) -> (Option<u64>, Option<u64>, i8) {
    if cds_regions.is_empty() {
        return (None, None, -1);
    }

    let cds_genomic_start = cds_regions.iter().map(|c| c.start).min().unwrap();
    let cds_genomic_end = cds_regions.iter().map(|c| c.end).max().unwrap();

    let (first_cds_pos, last_cds_pos) = if strand >= 0 {
        (cds_genomic_start, cds_genomic_end)
    } else {
        (cds_genomic_end, cds_genomic_start)
    };

    let cdna_coding_start = genomic_to_cdna(pairs, first_cds_pos, strand);
    let cdna_coding_end = genomic_to_cdna(pairs, last_cds_pos, strand);

    // Start phase comes from the first CDS feature in 5'->3' order.
    let start_phase = if strand >= 0 {
        cds_regions
            .iter()
            .min_by_key(|c| c.start)
            .map(|c| c.phase)
            .unwrap_or(-1)
    } else {
        cds_regions
            .iter()
            .max_by_key(|c| c.end)
            .map(|c| c.phase)
            .unwrap_or(-1)
    };

    (cdna_coding_start, cdna_coding_end, start_phase)
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
        }
    }

    fn make_cds(start: u64, end: u64, strand: i8, phase: i8) -> CdsRecord {
        CdsRecord {
            start,
            end,
            strand,
            phase,
        }
    }

    #[test]
    fn test_forward_strand_two_exons() {
        let exons = vec![make_exon(100, 200, 1), make_exon(300, 500, 1)];
        let cds = vec![];
        let result = build_mapper(&exons, &cds, 1);

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
        let result = build_mapper(&exons, &cds, -1);

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
        let result = build_mapper(&exons, &cds, 1);

        // CDS starts at genomic 150, which is in exon 1 (100-200).
        // Offset = 150 - 100 = 50, so cDNA coding start = 1 + 50 = 51.
        assert_eq!(result.cdna_coding_start, Some(51));

        // CDS ends at genomic 480, which is in exon 2 (300-500).
        // Offset = 480 - 300 = 180, so cDNA coding end = 102 + 180 = 282.
        assert_eq!(result.cdna_coding_end, Some(282));

        assert_eq!(result.start_phase, 0);
    }

    #[test]
    fn test_reverse_strand_with_cds() {
        // Reverse strand transcript: exons at 100-200 and 300-500.
        // CDS from 120-200 and 300-450.
        let exons = vec![make_exon(100, 200, -1), make_exon(300, 500, -1)];
        let cds = vec![make_cds(120, 200, -1, 1), make_cds(300, 450, -1, 0)];
        let result = build_mapper(&exons, &cds, -1);

        // Reverse strand: first CDS in 5'->3' is the one with highest genomic end.
        // CDS at 300-450: genomic end 450 is the 5' start.
        // In mapper pairs: pair[0] = 300-500 -> cDNA 1-201
        // Offset of genomic 450 from end: 500 - 450 = 50, cDNA = 1 + 50 = 51
        assert_eq!(result.cdna_coding_start, Some(51));

        // Last CDS in 5'->3' is CDS at 120-200: genomic start 120 is the 3' end.
        // In mapper pairs: pair[1] = 100-200 -> cDNA 202-302
        // Offset of genomic 120 from end: 200 - 120 = 80, cDNA = 202 + 80 = 282
        assert_eq!(result.cdna_coding_end, Some(282));

        // Start phase from CDS with highest end (300-450, phase=0).
        assert_eq!(result.start_phase, 0);
    }

    #[test]
    fn test_empty_exons() {
        let result = build_mapper(&[], &[], 1);
        assert!(result.pairs.is_empty());
        assert!(result.cdna_coding_start.is_none());
        assert!(result.cdna_coding_end.is_none());
        assert_eq!(result.start_phase, -1);
    }
}
