// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Consequence calculation engine for VEP.
//!
//! This crate implements the biological logic for determining the effect
//! of genomic variants on transcripts, including consequence assignment,
//! coordinate mapping, codon translation, and amino acid comparison.

pub mod coding;
pub mod consequences;
pub mod display;
pub mod hgvs;
pub mod mapper;
pub mod overlap;
pub mod sv;

pub use consequences::{calculate_consequences, EffectsConfig};
pub use overlap::find_overlapping_transcripts;

/// Test helpers for building mock transcripts. Available in tests only.
#[cfg(test)]
pub(crate) mod test_helpers {
    use std::sync::Arc;
    use vep_core::codon::reverse_complement;
    use vep_core::coordinate::Strand;
    use vep_core::transcript::*;
    use vep_core::variant::InputVariant;

    /// Build a test transcript on chr21, forward strand, with 3 exons and 2 introns.
    ///
    /// Layout (all 1-based, forward strand):
    ///   Exon 1:   25_000_000 - 25_000_299  (300bp, cDNA 1-300)
    ///   Intron 1: 25_000_300 - 25_001_999  (1700bp)
    ///   Exon 2:   25_002_000 - 25_002_299  (300bp, cDNA 301-600)
    ///   Intron 2: 25_002_300 - 25_003_999  (1700bp)
    ///   Exon 3:   25_004_000 - 25_006_000  (2001bp, cDNA 601-2601)
    ///
    /// CDS: cDNA positions 51-900 (850bp of coding sequence)
    ///   5' UTR: cDNA 1-50 (first 50bp of exon 1)
    ///   CDS:    cDNA 51-900 (last 250bp of exon 1 + all of exon 2 + first 300bp of exon 3)
    ///   3' UTR: cDNA 901-2601 (remaining 1701bp of exon 3)
    ///
    /// Translateable sequence: 850bp starting with ATG GCT GGA ... (Met Ala Gly ...)
    pub fn make_test_transcript() -> Transcript {
        // 850bp = 283 codons + 1 extra base; the first six codons are known.
        let mut cds = String::with_capacity(850);
        cds.push_str("ATG"); // codon 1: Met (M)
        cds.push_str("GCT"); // codon 2: Ala (A)
        cds.push_str("GGA"); // codon 3: Gly (G)
        cds.push_str("AAA"); // codon 4: Lys (K)
        cds.push_str("TTC"); // codon 5: Phe (F)
        cds.push_str("GAT"); // codon 6: Asp (D)
                             // GCT (Ala) fills the rest; up to 2bp overshoot is truncated.
        while cds.len() < 850 {
            cds.push_str("GCT");
        }
        cds.truncate(850);

        let exons = vec![
            Exon {
                stable_id: Some("ENSE00000000001".into()),
                start: 25_000_000,
                end: 25_000_299,
                rank: 1,
                phase: -1,
                end_phase: 0,
            },
            Exon {
                stable_id: Some("ENSE00000000002".into()),
                start: 25_002_000,
                end: 25_002_299,
                rank: 2,
                phase: 0,
                end_phase: 0,
            },
            Exon {
                stable_id: Some("ENSE00000000003".into()),
                start: 25_004_000,
                end: 25_006_000,
                rank: 3,
                phase: 0,
                end_phase: -1,
            },
        ];

        let introns = vec![
            Intron {
                start: 25_000_300,
                end: 25_001_999,
                rank: 1,
            },
            Intron {
                start: 25_002_300,
                end: 25_003_999,
                rank: 2,
            },
        ];

        let mapper_pairs = vec![
            MapperPair {
                from_start: 1,
                from_end: 300,
                to_start: 25_000_000,
                to_end: 25_000_299,
                ori: 1,
            },
            MapperPair {
                from_start: 301,
                from_end: 600,
                to_start: 25_002_000,
                to_end: 25_002_299,
                ori: 1,
            },
            MapperPair {
                from_start: 601,
                from_end: 2601,
                to_start: 25_004_000,
                to_end: 25_006_000,
                ori: 1,
            },
        ];

        let vefc = TranscriptVEFC {
            codon_table: 1,
            five_prime_utr: None,
            three_prime_utr: None,
            translateable_seq: Some(cds),
            peptide: None,
            introns: introns.clone(),
            sorted_exons: exons.clone(),
            mapper: Some(TranscriptMapper {
                start_phase: 0,
                cdna_coding_start: 51,
                cdna_coding_end: 900,
                exon_coord_mapper: ExonCoordMapper::new(mapper_pairs),
            }),
            protein_features: vec![],
            protein_function_predictions: None,
            seq_edits: vec![],
        };

        Transcript {
            stable_id: "ENST00000000001".into(),
            version: Some(1),
            db_id: None,
            gene_stable_id: "ENSG00000000001".into(),
            chr: "21".into(),
            start: 25_000_000,
            end: 25_006_000,
            strand: Strand::Forward,
            biotype: "protein_coding".into(),
            source: "Ensembl".into(),
            description: None,
            gene_symbol: Some("TEST1".into()),
            gene_symbol_source: Some("HGNC".into()),
            hgnc_id: Some("HGNC:0001".into()),
            gene_phenotype: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: Some("ENSP00000000001".into()),
            refseq: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            exons: exons.clone(),
            introns,
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(900),
            coding_region_start: Some(25_000_050),
            coding_region_end: Some(25_004_299),
            translation_start: Some(25_000_050),
            translation_end: Some(25_004_299),
            translation: Some(Translation {
                stable_id: "ENSP00000000001".into(),
                version: Some(1),
                db_id: None,
                start: 51,
                end: 300,
                start_exon_index: 0,
                end_exon_index: 2,
                seq: None,
            }),
            cdna_sequence: None,
            protein_sequence: None,
            flags: Arc::from([]),
            gencode_primary: false,
            attributes: vec![],
            vefc: Some(vefc),
            derived: Default::default(),
        }
    }

    /// The 5' UTR length of [`make_single_exon_transcript`].
    pub const SINGLE_EXON_UTR5: u64 = 30;

    /// A single-exon coding transcript on chromosome 21: a 30-base 5' UTR, the
    /// coding sequence `cds` (transcript orientation, ending in its stop codon) and
    /// the 3' UTR `utr3`, on `strand`. The genomic span starts at 40,000,001; on the
    /// reverse strand the transcript's first base is the span's last.
    pub fn make_single_exon_transcript(strand: Strand, cds: &str, utr3: &str) -> Transcript {
        const START: u64 = 40_000_001;
        let cds_len = cds.len() as u64;
        let cdna_len = SINGLE_EXON_UTR5 + cds_len + utr3.len() as u64;
        let end = START + cdna_len - 1;
        let (ori, crs, cre) = if strand == Strand::Reverse {
            (
                -1,
                end - (SINGLE_EXON_UTR5 + cds_len - 1),
                end - SINGLE_EXON_UTR5,
            )
        } else {
            (
                1,
                START + SINGLE_EXON_UTR5,
                START + SINGLE_EXON_UTR5 + cds_len - 1,
            )
        };
        let exons = vec![Exon {
            stable_id: Some("ENSE00000000001".into()),
            start: START,
            end,
            rank: 1,
            phase: -1,
            end_phase: -1,
        }];
        let pairs = vec![MapperPair {
            from_start: 1,
            from_end: cdna_len,
            to_start: START,
            to_end: end,
            ori,
        }];
        let vefc = TranscriptVEFC {
            codon_table: 1,
            five_prime_utr: None,
            three_prime_utr: Some(utr3.to_string()),
            translateable_seq: Some(cds.to_string()),
            peptide: None,
            introns: vec![],
            sorted_exons: exons.clone(),
            mapper: Some(TranscriptMapper {
                start_phase: 0,
                cdna_coding_start: SINGLE_EXON_UTR5 + 1,
                cdna_coding_end: SINGLE_EXON_UTR5 + cds_len,
                exon_coord_mapper: ExonCoordMapper::new(pairs),
            }),
            protein_features: vec![],
            protein_function_predictions: None,
            seq_edits: vec![],
        };
        Transcript {
            stable_id: "ENST00000000001".into(),
            version: Some(1),
            db_id: None,
            gene_stable_id: "ENSG00000000001".into(),
            chr: "21".into(),
            start: START,
            end,
            strand,
            biotype: "protein_coding".into(),
            source: "Ensembl".into(),
            description: None,
            gene_symbol: Some("TEST1".into()),
            gene_symbol_source: Some("HGNC".into()),
            hgnc_id: Some("HGNC:0001".into()),
            gene_phenotype: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: Some("ENSP00000000001".into()),
            refseq: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            exons,
            introns: vec![],
            cdna_coding_start: Some(SINGLE_EXON_UTR5 + 1),
            cdna_coding_end: Some(SINGLE_EXON_UTR5 + cds_len),
            coding_region_start: Some(crs),
            coding_region_end: Some(cre),
            translation_start: Some(crs),
            translation_end: Some(cre),
            translation: Some(Translation {
                stable_id: "ENSP00000000001".into(),
                version: Some(1),
                db_id: None,
                start: SINGLE_EXON_UTR5 + 1,
                end: SINGLE_EXON_UTR5 + cds_len,
                start_exon_index: 0,
                end_exon_index: 0,
                seq: None,
            }),
            cdna_sequence: None,
            protein_sequence: None,
            flags: Arc::from([]),
            gencode_primary: false,
            attributes: vec![],
            vefc: Some(vefc),
            derived: Default::default(),
        }
    }

    /// A 69-base coding sequence: `ATG`, twenty `GCT`, then `penultimate` and `stop`,
    /// so the stop codon is CDS 67..69 and the penultimate codon CDS 64..66.
    pub fn cds_ending_in(penultimate: &str, stop: &str) -> String {
        let mut cds = String::from("ATG");
        cds.push_str(&"GCT".repeat(20));
        cds.push_str(penultimate);
        cds.push_str(stop);
        cds
    }

    /// The genomic coordinate of CDS position `cds_pos` on a single-exon transcript.
    pub fn genomic_of_cds(tx: &Transcript, cds_pos: u64) -> u64 {
        let cdna = SINGLE_EXON_UTR5 + cds_pos;
        if tx.strand == Strand::Reverse {
            tx.end - (cdna - 1)
        } else {
            tx.start + cdna - 1
        }
    }

    /// An insertion of `seq` (transcript orientation) between CDS positions `after`
    /// and `after + 1`, as the input record carries it on the forward strand.
    pub fn insertion_after_cds(tx: &Transcript, after: u64, seq: &str) -> InputVariant {
        let a = genomic_of_cds(tx, after);
        let b = genomic_of_cds(tx, after + 1);
        let anchor = a.min(b);
        let alt = if tx.strand == Strand::Reverse {
            reverse_complement(seq.as_bytes())
        } else {
            seq.as_bytes().to_vec()
        };
        InputVariant::new("21".into(), anchor + 1, anchor, b"-".to_vec(), alt)
    }

    /// A deletion of CDS positions `from..=to` (transcript orientation; positions
    /// past the CDS end run into the 3' UTR), `seq` being the deleted transcript bases.
    pub fn deletion_of_cds(tx: &Transcript, from: u64, to: u64, seq: &str) -> InputVariant {
        let a = genomic_of_cds(tx, from);
        let b = genomic_of_cds(tx, to);
        let r = if tx.strand == Strand::Reverse {
            reverse_complement(seq.as_bytes())
        } else {
            seq.as_bytes().to_vec()
        };
        InputVariant::new("21".into(), a.min(b), a.max(b), r, b"-".to_vec())
    }

    /// Build the standard test transcript with transcript flags applied.
    pub fn make_test_transcript_with_flags(flags: &[&str]) -> Transcript {
        let mut tx = make_test_transcript();
        tx.flags = flags
            .iter()
            .map(|flag| (*flag).to_string())
            .collect::<Vec<_>>()
            .into();
        tx
    }

    /// Assert that a transcript consequence carries exactly the given Sequence
    /// Ontology term set, in any order.
    ///
    /// The concordance comparator scores the whole `consequence_set` of a tuple,
    /// so a test that only checks `contains` can pass while the tuple still
    /// mismatches VEP on a term it did not look at. `expected` takes the exact
    /// VEP strings (`"5_prime_UTR_variant"`); the failure message shows both
    /// sets sorted so the missing and extra terms read off directly.
    pub fn assert_consequence_set_eq(
        tc: &vep_core::consequence::TranscriptConsequence,
        expected: &[&str],
    ) {
        let mut got: Vec<String> = tc.consequences.iter().map(|c| c.to_string()).collect();
        got.sort();
        got.dedup();
        let mut want: Vec<String> = expected.iter().map(|s| (*s).to_string()).collect();
        want.sort();
        want.dedup();
        assert_eq!(
            got,
            want,
            "consequence set mismatch for {}: missing {:?}, extra {:?}",
            tc.transcript_id,
            want.iter().filter(|t| !got.contains(t)).collect::<Vec<_>>(),
            got.iter().filter(|t| !want.contains(t)).collect::<Vec<_>>(),
        );
    }
}

#[cfg(test)]
mod concordance_tests;
