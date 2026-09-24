// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Benchmark for consequence calculation.
//!
//! Builds a realistic test transcript and benchmarks the core
//! `calculate_consequences` function for different variant types.

use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use vep_core::coordinate::Strand;
use vep_core::transcript::*;
use vep_core::variant::InputVariant;
use vep_effects::{calculate_consequences, find_overlapping_transcripts, EffectsConfig};

/// Build a test transcript on chr21 with 3 exons.
///
/// Layout mirrors the test_helpers in vep-effects (forward strand):
///   Exon 1:   25_000_000 - 25_000_299
///   Intron 1: 25_000_300 - 25_001_999
///   Exon 2:   25_002_000 - 25_002_299
///   Intron 2: 25_002_300 - 25_003_999
///   Exon 3:   25_004_000 - 25_006_000
///   CDS: cDNA 51-900
fn build_bench_transcript() -> Transcript {
    let mut cds = String::with_capacity(850);
    cds.push_str("ATG");
    cds.push_str("GCT");
    cds.push_str("GGA");
    cds.push_str("AAA");
    cds.push_str("TTC");
    cds.push_str("GAT");
    while cds.len() < 849 {
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
        flags: vec![].into(),
        gencode_primary: false,
        attributes: vec![],
        vefc: Some(vefc),
        derived: Default::default(),
    }
}

fn consequence_benchmark(c: &mut Criterion) {
    let transcript = build_bench_transcript();
    let config = EffectsConfig::default();

    let snv = InputVariant::new(
        "21".into(),
        25_000_060,
        25_000_060,
        b"G".to_vec(),
        b"T".to_vec(),
    );

    let deletion = InputVariant::new(
        "21".into(),
        25_000_295,
        25_000_305,
        b"ACGTACGTACG".to_vec(),
        b"-".to_vec(),
    );

    let insertion = InputVariant::new(
        "21".into(),
        25_002_100,
        25_002_099,
        b"-".to_vec(),
        b"AAA".to_vec(),
    );

    let upstream = InputVariant::new(
        "21".into(),
        24_998_000,
        24_998_000,
        b"A".to_vec(),
        b"G".to_vec(),
    );

    c.bench_function("calculate_consequences_snv", |b| {
        b.iter(|| calculate_consequences(black_box(&snv), black_box(&transcript), &config))
    });

    c.bench_function("calculate_consequences_deletion", |b| {
        b.iter(|| calculate_consequences(black_box(&deletion), black_box(&transcript), &config))
    });

    c.bench_function("calculate_consequences_insertion", |b| {
        b.iter(|| calculate_consequences(black_box(&insertion), black_box(&transcript), &config))
    });

    c.bench_function("calculate_consequences_upstream", |b| {
        b.iter(|| calculate_consequences(black_box(&upstream), black_box(&transcript), &config))
    });

    let transcripts = vec![transcript];
    c.bench_function("find_overlapping_transcripts", |b| {
        b.iter(|| {
            find_overlapping_transcripts(
                black_box(&snv),
                black_box(&transcripts),
                config.upstream_distance,
                config.downstream_distance,
            )
        })
    });
}

criterion_group!(benches, consequence_benchmark);
criterion_main!(benches);
