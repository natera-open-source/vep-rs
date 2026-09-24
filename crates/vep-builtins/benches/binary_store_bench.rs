// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Benchmark: binary annotation store vs tabix lookup performance.
//!
//! Generates synthetic annotation data, writes it to a binary store,
//! and benchmarks batch and single-position queries.
//!
//! Run with `cargo bench --package vep-builtins --bench binary_store_bench`.

use std::sync::Arc;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;

use vep_builtins::annotation_store::AnnotationStore;
use vep_builtins::binary_store::{BinaryAnnotator, BinaryStoreWriter};
use vep_builtins::tabix::TabixRecord;

/// Generate synthetic records: `record_count` records across `chr_count` chromosomes,
/// evenly spaced starting at position 1000.
fn generate_records(
    chr_count: usize,
    records_per_chr: usize,
    col_count: usize,
) -> (Vec<String>, Vec<(String, Vec<TabixRecord>)>) {
    let column_names: Vec<String> = (0..col_count).map(|i| format!("col_{i}")).collect();

    let header = Arc::new(column_names.clone());
    let mut grouped = Vec::with_capacity(chr_count);
    for c in 0..chr_count {
        let chr = format!("{}", c + 1);
        let mut records = Vec::with_capacity(records_per_chr);
        for i in 0..records_per_chr {
            let pos = 1000 + (i as u64) * 3; // positions: 1000, 1003, 1006, ...
            let columns: Vec<String> = (0..col_count).map(|j| format!("{}.{}", pos, j)).collect();
            records.push(TabixRecord {
                chr: chr.clone(),
                start: pos,
                end: pos,
                columns,
                header: Some(Arc::clone(&header)),
            });
        }
        grouped.push((chr, records));
    }

    (column_names, grouped)
}

fn bench_binary_store_queries(c: &mut Criterion) {
    let mut group = c.benchmark_group("binary_store");

    for &records_per_chr in &[1_000, 10_000, 100_000] {
        let chr_count = 3;
        let col_count = 6;

        let (column_names, grouped) = generate_records(chr_count, records_per_chr, col_count);

        let dir = tempfile::tempdir().unwrap();
        let vpd_path = dir.path().join("bench.vpd");
        let writer = BinaryStoreWriter::new(vpd_path.clone());
        writer.write(&column_names, &grouped).unwrap();

        let store = BinaryAnnotator::open(&vpd_path).unwrap();

        group.bench_with_input(
            BenchmarkId::new("single_query", records_per_chr),
            &records_per_chr,
            |b, _| {
                b.iter(|| {
                    let result = store.query(black_box("1"), black_box(1000), black_box(1000));
                    black_box(result.unwrap());
                });
            },
        );

        let batch_size = 5000.min(records_per_chr);
        let regions: Vec<(&str, u64, u64)> = (0..batch_size)
            .map(|i| {
                let pos = 1000 + (i as u64) * 3;
                ("1", pos, pos)
            })
            .collect();

        group.bench_with_input(
            BenchmarkId::new("batch_query_5000", records_per_chr),
            &records_per_chr,
            |b, _| {
                b.iter(|| {
                    let result = store.query_batch(black_box(&regions));
                    black_box(result.unwrap());
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_binary_store_queries);
criterion_main!(benches);
