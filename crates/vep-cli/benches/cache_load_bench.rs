// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Benchmark: JSON cache loading performance.
//!
//! Measures transcript and variation cache loading from converted JSON files.
//! Uses the full cache directory as the benchmark fixture.
//!
//! To run:
//!   cargo bench --package vep-cli --bench cache_load_bench
//!
//! The cache root is `$VEP_JSON_CACHE_DIR` when set, else
//! `tmp/concordance/cache_json/` under the workspace root. `$VEP_BENCH_CHR`
//! (default `21`) names the chromosome the per-chromosome benchmark parses,
//! which is the unit `LazyTranscriptIndexes::prewarm` loads at run time.
//! Benchmarks are silently skipped if the cache is not available.

use std::path::Path;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};

use vep_cli::json_cache;

/// Root of the JSON cache (relative to workspace root).
const CACHE_DIR: &str = "tmp/concordance/cache_json";

fn cache_dir() -> Option<String> {
    let candidates = [
        std::env::var("VEP_JSON_CACHE_DIR").unwrap_or_default(),
        CACHE_DIR.to_string(),
        format!(
            "{}/{}",
            env!("CARGO_MANIFEST_DIR").trim_end_matches("/crates/vep-cli"),
            CACHE_DIR
        ),
    ];
    for c in &candidates {
        if !c.is_empty() && Path::new(c).join("transcripts").is_dir() {
            return Some(c.clone());
        }
    }
    None
}

fn bench_load_transcripts_for_chr(c: &mut Criterion) {
    let Some(dir) = cache_dir() else {
        eprintln!("SKIP: JSON cache not found at {CACHE_DIR}; run vep-cache-converter first");
        return;
    };
    let chr = std::env::var("VEP_BENCH_CHR").unwrap_or_else(|_| "21".to_string());
    let shards = json_cache::enumerate_transcript_shards(&dir).expect("shard enumeration failed");
    let Some(paths) = shards.get(&chr) else {
        eprintln!("SKIP: chromosome {chr} not in cache {dir}");
        return;
    };

    let mut group = c.benchmark_group("cache_loading");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(20));
    group.bench_function(format!("transcripts_chr{chr}"), |b| {
        b.iter(|| {
            json_cache::load_transcripts_for_chr(&chr, paths).expect("transcript loading failed");
        });
    });
    group.finish();
}

fn bench_load_transcripts(c: &mut Criterion) {
    let Some(dir) = cache_dir() else {
        return;
    };

    let mut group = c.benchmark_group("cache_loading");
    group.sample_size(10);
    group.bench_function("transcripts", |b| {
        b.iter(|| {
            json_cache::load_all_transcripts(&dir).expect("transcript loading failed");
        });
    });
    group.finish();
}

fn bench_load_variations(c: &mut Criterion) {
    let Some(dir) = cache_dir() else {
        return;
    };

    let var_dir = Path::new(&dir).join("variations");
    if !var_dir.is_dir() {
        eprintln!("SKIP: No variations directory in cache");
        return;
    }

    c.bench_function("cache_loading/variations", |b| {
        b.iter(|| {
            json_cache::load_all_variations(&dir).expect("variation loading failed");
        });
    });
}

criterion_group!(
    benches,
    bench_load_transcripts_for_chr,
    bench_load_transcripts,
    bench_load_variations
);
criterion_main!(benches);
