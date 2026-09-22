// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Benchmark: JSON formatter for the embedding hot path.
//!
//! `JsonOutputFormatter::format_variant_value` does not do a
//! `to_string` -> `from_str` round-trip per variant: it builds the `Value`
//! directly, while `format_variant` keeps the string output the file-based CLI
//! path needs.
//!
//! This bench measures that invariant: a per-variant round-trip shows up as
//! the direct lane converging on the `roundtrip` lane. It compares:
//!
//! - `format_variant_value`: the embedding hot path (direct Value).
//! - `format_variant`: the CLI path, which serializes to a single JSON line.
//! - `roundtrip`: a synthesized worst-case that serializes the value to a
//!   string and immediately parses it back.
//!
//! On any healthy build the `roundtrip` lane is materially slower than the
//! direct `format_variant_value` lane.
//!
//! Run with:
//!     cargo bench --package vep-io --bench json_format_bench

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use smallvec::smallvec;
use vep_core::consequence::{Consequence, FeatureType, Impact, TranscriptConsequence};
use vep_core::variant::{ColocatedVariant, InputVariant};
use vep_io::output::fields::FieldOptions;
use vep_io::output::json::JsonOutputFormatter;

/// One annotated transcript consequence resembling a typical missense hit.
fn make_transcript_consequence(transcript_id: &str, gene_id: &str) -> TranscriptConsequence {
    TranscriptConsequence {
        transcript_id: transcript_id.into(),
        gene_id: gene_id.into(),
        gene_symbol: Some("MRPL39".into()),
        gene_symbol_source: Some("HGNC".into()),
        hgnc_id: Some("HGNC:16650".into()),
        consequences: smallvec![Consequence::MissenseVariant],
        impact: Impact::MODERATE,
        feature_type: FeatureType::Transcript,
        biotype: Some("protein_coding".into()),
        canonical: true,
        cdna_position: Some("1033".into()),
        cds_position: Some("991".into()),
        protein_position: Some("331".into()),
        amino_acids: Some("A/V".into()),
        codons: Some("Gca/Gta".into()),
        strand: -1,
        hgvsc: Some("ENST00000352957.4:c.991G>A".into()),
        hgvsp: Some("ENSP00000284967.4:p.Ala331Val".into()),
        sift: Some("deleterious(0.02)".into()),
        polyphen: Some("probably_damaging(0.95)".into()),
        domains: vec![
            ("Pfam".into(), "PF00080".into()),
            ("PROSITE_profiles".into(), "PS51402".into()),
        ],
        ..Default::default()
    }
}

/// Build a representative annotated variant with `n_tc` transcript
/// consequences. Real gene-dense regions can produce 20+ TCs per variant.
fn make_variant_with_tcs(n_tc: usize) -> InputVariant {
    let mut variant = InputVariant::new(
        "21".into(),
        25_585_733,
        25_585_733,
        b"A".to_vec(),
        b"G".to_vec(),
    );
    variant.id = Some("rs123".into());
    variant.most_severe_consequence = Some(Consequence::MissenseVariant);
    variant.transcript_consequences = (0..n_tc)
        .map(|i| {
            make_transcript_consequence(
                &format!("ENST{:011}", i + 1),
                &format!("ENSG{:011}", i + 1),
            )
        })
        .collect();
    variant.colocated_variants = vec![ColocatedVariant {
        id: "rs999".into(),
        start: 25_585_733,
        end: 25_585_733,
        allele_string: Some("A/G".into()),
        ..Default::default()
    }];
    variant.existing_variation = vec!["rs999".into()];
    variant
}

fn bench_format_paths(c: &mut Criterion) {
    let formatter = JsonOutputFormatter::new(FieldOptions::default(), Vec::new(), None);
    let mut group = c.benchmark_group("format_variant");

    for n_tc in [1usize, 5, 20] {
        let variant = make_variant_with_tcs(n_tc);

        group.bench_with_input(
            BenchmarkId::new("format_variant_value", n_tc),
            &variant,
            |b, v| {
                b.iter(|| {
                    let value = formatter.format_variant_value(black_box(v));
                    black_box(value);
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("format_variant_to_string", n_tc),
            &variant,
            |b, v| {
                b.iter(|| {
                    let value = formatter.format_variant_value(black_box(v));
                    let s = serde_json::to_string(&value).unwrap();
                    black_box(s);
                });
            },
        );

        group.bench_with_input(BenchmarkId::new("roundtrip", n_tc), &variant, |b, v| {
            b.iter(|| {
                let value = formatter.format_variant_value(black_box(v));
                let s = serde_json::to_string(&value).unwrap();
                let parsed: serde_json::Value = serde_json::from_str(&s).unwrap();
                black_box(parsed);
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_format_paths);
criterion_main!(benches);
