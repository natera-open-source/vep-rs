// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Transcript-variant overlap detection.
//!
//! Determines which transcripts a variant could affect by comparing the
//! variant's genomic interval against each transcript's extended region
//! (span +/- upstream/downstream distance, adjusted for strand).
//! Used as the first step in the annotation pipeline before coordinate
//! mapping and consequence calculation.

use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

/// Find all transcripts that a variant could affect.
///
/// Checks if the variant's genomic interval overlaps the transcript's
/// extended region (transcript span +/- upstream/downstream distance,
/// adjusted for strand).
pub fn find_overlapping_transcripts<'a>(
    variant: &InputVariant,
    transcripts: &'a [Transcript],
    upstream_distance: u64,
    downstream_distance: u64,
) -> Vec<&'a Transcript> {
    let variant_start = variant.start.min(variant.end);
    let variant_end = variant.start.max(variant.end);
    transcripts
        .iter()
        .filter(|tx| {
            if tx.chr != variant.chr {
                return false;
            }
            let (region_start, region_end) =
                extended_region(tx, upstream_distance, downstream_distance);
            variant_start <= region_end && variant_end >= region_start
        })
        .collect()
}

/// Compute the extended genomic region for a transcript, accounting for strand.
///
/// Forward strand: upstream is before start, downstream is after end.
/// Reverse strand: upstream is after end, downstream is before start.
fn extended_region(tx: &Transcript, upstream_dist: u64, downstream_dist: u64) -> (u64, u64) {
    match tx.strand {
        Strand::Forward => {
            let start = tx.start.saturating_sub(upstream_dist);
            let end = tx.end.saturating_add(downstream_dist);
            (start, end)
        }
        Strand::Reverse => {
            let start = tx.start.saturating_sub(downstream_dist);
            let end = tx.end.saturating_add(upstream_dist);
            (start, end)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use vep_core::coordinate::Strand;

    fn make_transcript(chr: &str, start: u64, end: u64, strand: Strand) -> Transcript {
        Transcript {
            stable_id: "ENST00000000001".into(),
            version: None,
            db_id: None,
            gene_stable_id: "ENSG00000000001".into(),
            chr: chr.into(),
            start,
            end,
            strand,
            biotype: "protein_coding".into(),
            source: "Ensembl".into(),
            description: None,
            gene_symbol: None,
            gene_symbol_source: None,
            hgnc_id: None,
            gene_phenotype: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            refseq: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            exons: vec![],
            introns: vec![],
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            translation_start: None,
            translation_end: None,
            translation: None,
            cdna_sequence: None,
            protein_sequence: None,
            flags: Arc::from([]),
            gencode_primary: false,
            attributes: vec![],
            vefc: None,
            derived: Default::default(),
        }
    }

    #[test]
    fn test_overlap_within_transcript() {
        let txs = [make_transcript(
            "21",
            25_000_000,
            25_010_000,
            Strand::Forward,
        )];
        let variant = InputVariant::new(
            "21".into(),
            25_005_000,
            25_005_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = find_overlapping_transcripts(&variant, &txs, 5000, 5000);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_overlap_upstream_forward() {
        let txs = [make_transcript(
            "21",
            25_000_000,
            25_010_000,
            Strand::Forward,
        )];
        // 3000bp upstream of forward-strand transcript
        let variant = InputVariant::new(
            "21".into(),
            24_997_000,
            24_997_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = find_overlapping_transcripts(&variant, &txs, 5000, 5000);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_no_overlap_too_far_upstream() {
        let txs = [make_transcript(
            "21",
            25_000_000,
            25_010_000,
            Strand::Forward,
        )];
        let variant = InputVariant::new(
            "21".into(),
            24_990_000,
            24_990_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = find_overlapping_transcripts(&variant, &txs, 5000, 5000);
        assert!(result.is_empty());
    }

    #[test]
    fn test_no_overlap_different_chr() {
        let txs = [make_transcript(
            "22",
            25_000_000,
            25_010_000,
            Strand::Forward,
        )];
        let variant = InputVariant::new(
            "21".into(),
            25_005_000,
            25_005_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = find_overlapping_transcripts(&variant, &txs, 5000, 5000);
        assert!(result.is_empty());
    }

    #[test]
    fn test_overlap_downstream_reverse_strand() {
        // For reverse strand, "downstream" (3') is before the transcript start
        let txs = [make_transcript(
            "21",
            25_000_000,
            25_010_000,
            Strand::Reverse,
        )];
        let variant = InputVariant::new(
            "21".into(),
            24_996_000,
            24_996_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = find_overlapping_transcripts(&variant, &txs, 5000, 5000);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_multiple_transcripts() {
        let transcripts = vec![
            make_transcript("21", 25_000_000, 25_010_000, Strand::Forward),
            make_transcript("21", 25_100_000, 25_110_000, Strand::Forward),
            make_transcript("21", 25_005_000, 25_015_000, Strand::Reverse),
        ];
        let variant = InputVariant::new(
            "21".into(),
            25_005_000,
            25_005_000,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let result = find_overlapping_transcripts(&variant, &transcripts, 5000, 5000);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_insertion_interval_overlap_with_reversed_coordinates() {
        let txs = [make_transcript(
            "21",
            25_000_000,
            25_010_000,
            Strand::Forward,
        )];
        // VEP-style insertion coordinates use start > end.
        let variant = InputVariant::new(
            "21".into(),
            25_000_000,
            24_999_999,
            b"-".to_vec(),
            b"T".to_vec(),
        );
        let result = find_overlapping_transcripts(&variant, &txs, 0, 0);
        assert_eq!(result.len(), 1);
    }
}
