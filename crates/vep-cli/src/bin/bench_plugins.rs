// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Plugin overhead benchmark binary.
//!
//! Runs the VEP annotation pipeline on real data and emits structured
//! JSON timing reports decomposed by pipeline phase. Designed to measure
//! plugin overhead across multiple dimensions:
//!
//! - Per-plugin marginal cost
//! - Parallel vs sequential prefetch
//! - Tabix vs binary store backends
//! - Buffer size scaling
//!
//! # Usage
//!
//! ```bash
//! cargo build --release --bin bench_plugins
//!
//! target/release/bench_plugins \
//!   --input <input.vcf.gz> \
//!   --json-cache <json_cache_dir> \
//!   --plugin CADD,snv=<cadd_snvs.tsv.gz>,indels=<cadd_indels.tsv.gz> \
//!   --plugin REVEL,file=<revel.tsv.gz> \
//!   --variants 50000 \
//!   --runs 3 \
//!   --output-json results.json
//! ```

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::Parser;
use serde::Serialize;

use vep_builtins::{BuiltinRegistry, PluginTimings};
use vep_cli::json_cache;
use vep_cli::transcript_index::TranscriptBinIndex;
use vep_cli::variation_matcher;
use vep_core::allele::trim_alleles;
use vep_core::consequence::Consequence;
use vep_core::variant::InputVariant;

#[derive(Parser, Debug)]
#[command(
    name = "bench_plugins",
    about = "Benchmark plugin overhead in the VEP annotation pipeline"
)]
struct Args {
    /// Input VCF file path.
    #[arg(short = 'i', long)]
    input: String,

    /// Path to converted JSON cache directory.
    #[arg(long = "json-cache")]
    json_cache: String,

    /// Plugin specifications (repeatable). Same format as vep --plugin.
    #[arg(long = "plugin", action = clap::ArgAction::Append)]
    plugin: Vec<String>,

    /// Maximum number of variants to load from the input file.
    #[arg(long, default_value_t = 50000)]
    variants: usize,

    /// Number of timed runs (results are aggregated).
    #[arg(long, default_value_t = 3)]
    runs: usize,

    /// Buffer size for batch processing.
    #[arg(long = "buffer-size", default_value_t = 5000)]
    buffer_size: usize,

    /// Number of threads for consequence annotation.
    #[arg(long, default_value_t = 1)]
    threads: usize,

    /// Force sequential plugin execution (no parallel prefetch).
    #[arg(long, default_value_t = false)]
    sequential: bool,

    /// Output JSON report path. Prints to stdout if not specified.
    #[arg(long = "output-json")]
    output_json: Option<String>,

    /// Label for this benchmark run (included in the report).
    #[arg(long)]
    label: Option<String>,
}

#[derive(Serialize)]
struct BenchReport {
    config: BenchConfig,
    phases: PhaseTimings,
    plugin_detail: HashMap<String, PluginDetail>,
    runs: Vec<RunTimings>,
    derived: DerivedMetrics,
}

#[derive(Serialize)]
struct BenchConfig {
    input: String,
    variants_loaded: usize,
    plugins: Vec<String>,
    buffer_size: usize,
    threads: usize,
    runs: usize,
    mode: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    label: Option<String>,
}

#[derive(Serialize)]
struct PhaseTimings {
    cache_load_ms: u64,
    parse_ms: u64,
    consequence_ms: AggregatedTiming,
    plugin_prefetch_ms: AggregatedTiming,
    plugin_annotate_ms: AggregatedTiming,
    plugin_total_ms: AggregatedTiming,
    total_ms: AggregatedTiming,
}

#[derive(Serialize, Default)]
struct AggregatedTiming {
    mean: f64,
    min: u64,
    max: u64,
}

#[derive(Serialize, Clone)]
struct PluginDetail {
    prefetch_ms: AggregatedTimingF64,
    annotate_ms: AggregatedTimingF64,
}

#[derive(Serialize, Clone, Default)]
struct AggregatedTimingF64 {
    mean: f64,
    min: f64,
    max: f64,
}

#[derive(Serialize, Clone)]
struct RunTimings {
    run: usize,
    consequence_ms: u64,
    plugin_prefetch_ms: u64,
    plugin_annotate_ms: u64,
    plugin_total_ms: u64,
    total_ms: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    plugin_timings: Option<PluginTimings>,
}

#[derive(Serialize)]
struct DerivedMetrics {
    variants_per_second: f64,
    plugin_overhead_pct: f64,
    plugin_overhead_ms: f64,
    /// Only present if baseline (no-plugin) run data is embedded.
    #[serde(skip_serializing_if = "Option::is_none")]
    consequence_only_variants_per_second: Option<f64>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let subscriber = tracing_subscriber::fmt()
        .with_max_level(tracing::Level::INFO)
        .with_writer(std::io::stderr)
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("failed to set tracing subscriber");

    // A scoped pool, not the global one, so this binary coexists with library
    // code that owns a pool per Annotator.
    let num_threads = args.threads.max(1);
    let pool = vep_cli::annotator::build_thread_pool(num_threads)?;

    eprintln!("=== Plugin Overhead Benchmark ===");
    eprintln!(
        "Input: {}, Variants: {}, Runs: {}, Threads: {}, Buffer: {}, Mode: {}",
        args.input,
        args.variants,
        args.runs,
        num_threads,
        args.buffer_size,
        if args.sequential {
            "sequential"
        } else {
            "parallel"
        }
    );

    let cache_start = Instant::now();
    let effects_config = Arc::new(vep_effects::EffectsConfig {
        upstream_distance: 5000,
        downstream_distance: 5000,
        reference_fasta: None,
        enable_indel_3prime_shift: false,
        max_indel_3prime_shift: 2000,
        populate_loftee_context: false,
        compute_hgvs: false,
        compute_exon_intron_numbers: true,
    });

    let raw_transcripts = json_cache::load_all_transcripts(&args.json_cache)?;
    let mut indexed: HashMap<String, TranscriptBinIndex> = HashMap::new();
    for (chr, txs) in raw_transcripts {
        indexed.insert(
            chr,
            TranscriptBinIndex::new(
                txs,
                effects_config.upstream_distance,
                effects_config.downstream_distance,
            ),
        );
    }
    let transcripts = Arc::new(indexed);

    let variation_data: Option<Arc<json_cache::VariationData>> = {
        let vd = json_cache::load_all_variations(&args.json_cache)?;
        if vd.variations.is_empty() {
            None
        } else {
            Some(Arc::new(vd))
        }
    };

    let cache_load_ms = cache_start.elapsed().as_millis() as u64;
    eprintln!("Cache loaded in {cache_load_ms} ms");

    let mut registry = BuiltinRegistry::new();
    if !args.plugin.is_empty() {
        let unresolved = registry.load_from_args(&args.plugin)?;
        if !unresolved.is_empty() {
            eprintln!("Warning: unresolved plugins (skipped): {:?}", unresolved);
        }
    }
    let plugin_names = registry.plugin_names();
    eprintln!(
        "Plugins: {} ({})",
        plugin_names.len(),
        if plugin_names.is_empty() {
            "none".to_string()
        } else {
            plugin_names.join(", ")
        }
    );

    let parse_start = Instant::now();
    let variants = load_variants(&args.input, args.variants)?;
    let parse_ms = parse_start.elapsed().as_millis() as u64;
    let variant_count = variants.len();
    eprintln!("Parsed {} variants in {} ms", variant_count, parse_ms);

    if variant_count == 0 {
        bail!("No variants loaded from input file");
    }

    let use_parallel = num_threads > 1;
    let mut run_results: Vec<RunTimings> = Vec::with_capacity(args.runs);

    for run_idx in 0..args.runs {
        eprintln!("--- Run {}/{} ---", run_idx + 1, args.runs);

        let mut batch = variants.clone();

        let run_start = Instant::now();

        let csq_start = Instant::now();
        pool.install(|| {
            annotate_consequences(
                &mut batch,
                &transcripts,
                &variation_data,
                &effects_config,
                use_parallel,
            );
        });
        let csq_ms = csq_start.elapsed().as_millis() as u64;

        let plugin_timings = if registry.plugin_count() > 0 {
            if args.sequential {
                Some(pool.install(|| registry.run_batch_sequential(&mut batch))?)
            } else {
                Some(pool.install(|| registry.run_batch_timed(&mut batch))?)
            }
        } else {
            None
        };

        let total_ms = run_start.elapsed().as_millis() as u64;
        let plugin_prefetch_ms = plugin_timings
            .as_ref()
            .map(|t| t.total_prefetch_ms)
            .unwrap_or(0);
        let plugin_annotate_ms = plugin_timings
            .as_ref()
            .map(|t| t.total_annotate_ms)
            .unwrap_or(0);
        let plugin_total_ms = plugin_timings.as_ref().map(|t| t.total_ms).unwrap_or(0);

        eprintln!(
            "  consequence: {} ms, plugin_prefetch: {} ms, plugin_annotate: {} ms, total: {} ms",
            csq_ms, plugin_prefetch_ms, plugin_annotate_ms, total_ms
        );

        run_results.push(RunTimings {
            run: run_idx + 1,
            consequence_ms: csq_ms,
            plugin_prefetch_ms,
            plugin_annotate_ms,
            plugin_total_ms,
            total_ms,
            plugin_timings,
        });
    }

    let report = build_report(
        &args,
        variant_count,
        cache_load_ms,
        parse_ms,
        &plugin_names,
        &run_results,
    );

    let json = serde_json::to_string_pretty(&report)?;

    if let Some(ref path) = args.output_json {
        std::fs::write(path, &json).with_context(|| format!("failed to write report to {path}"))?;
        eprintln!("Report written to {path}");
    } else {
        println!("{json}");
    }

    eprintln!(
        "=== Done: {:.0} variants/sec, plugin overhead {:.1}% ===",
        report.derived.variants_per_second, report.derived.plugin_overhead_pct
    );

    Ok(())
}

/// Parse VCF variants from an input file.
///
/// Supports both plain `.vcf` and compressed `.vcf.gz`/`.bgz` files.
/// Compressed files are decompressed via `zcat` / `gzcat` subprocess.
fn load_variants(path: &str, max_variants: usize) -> Result<Vec<InputVariant>> {
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") || path.ends_with(".bgz") {
        // Decompress via subprocess (avoids needing noodles::bgzf in the binary).
        let child = std::process::Command::new("zcat")
            .arg(path)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .or_else(|_| {
                // macOS: try gzcat as fallback
                std::process::Command::new("gzcat")
                    .arg(path)
                    .stdout(std::process::Stdio::piped())
                    .spawn()
            })
            .with_context(|| format!("failed to decompress {path} (tried zcat/gzcat)"))?;
        Box::new(BufReader::new(child.stdout.unwrap()))
    } else {
        let file = std::fs::File::open(path).with_context(|| format!("failed to open {path}"))?;
        Box::new(BufReader::new(file))
    };

    let mut variants = Vec::new();
    for line_result in reader.lines() {
        if variants.len() >= max_variants {
            break;
        }
        let line = line_result.context("failed to read line")?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        match parse_vcf_line(trimmed) {
            Ok(vs) => variants.extend(vs),
            Err(_) => continue,
        }
    }

    Ok(variants)
}

/// Minimal VCF line parser: the allele-trimming subset of `vcf_parser::parse_vcf_line`,
/// without symbolic-allele or raw-input handling.
fn parse_vcf_line(line: &str) -> Result<Vec<InputVariant>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        bail!("VCF line has fewer than 8 fields");
    }

    let chr = fields[0].to_string();
    let pos: u64 = fields[1].parse().context("invalid POS")?;
    let ref_allele = fields[3];
    let alt_field = fields[4];

    if alt_field == "." || alt_field == "*" {
        bail!("No alternate allele");
    }

    let alts: Vec<&str> = alt_field.split(',').collect();
    let mut variants = Vec::new();

    for (i, alt) in alts.iter().enumerate() {
        if *alt == "*" || *alt == "." {
            continue;
        }

        let (trimmed_ref, trimmed_alt, prefix_trimmed) =
            trim_alleles(ref_allele.as_bytes(), alt.as_bytes());

        let start = pos + prefix_trimmed as u64;
        let end = if trimmed_ref == b"-" {
            start.saturating_sub(1)
        } else {
            start + trimmed_ref.len() as u64 - 1
        };

        let mut variant = InputVariant::new(chr.clone(), start, end, trimmed_ref, trimmed_alt);

        variant.allele_index = i;
        let var_id = fields[2];
        if var_id != "." {
            variant.id = Some(var_id.to_string());
        }

        variants.push(variant);
    }

    Ok(variants)
}

/// Run consequence annotation on all variants (mirrors annotate_batch in runner.rs).
fn annotate_consequences(
    variants: &mut [InputVariant],
    transcripts: &Arc<HashMap<String, TranscriptBinIndex>>,
    variation_data: &Option<Arc<json_cache::VariationData>>,
    effects_config: &Arc<vep_effects::EffectsConfig>,
    use_parallel: bool,
) {
    use rayon::prelude::*;

    let annotate_one = |variant: &mut InputVariant| {
        if let Some(chr_transcripts) = transcripts.get(&variant.chr) {
            let variant_start = variant.start.min(variant.end);
            let variant_end = variant.start.max(variant.end);
            chr_transcripts.for_each_overlapping(variant_start, variant_end, |transcript| {
                if let Some(tc) =
                    vep_effects::calculate_consequences(variant, transcript, effects_config)
                {
                    variant.transcript_consequences.push(tc);
                }
            });
        }

        if let Some(ref vd) = variation_data {
            if let Some(chr_vars) = vd.variations.get(&variant.chr) {
                if let Some(chr_idx) = vd.position_index.get(&variant.chr) {
                    let colocated =
                        variation_matcher::find_colocated_variants(variant, chr_vars, chr_idx);
                    for cv in &colocated {
                        variant.existing_variation.push(cv.id.clone());
                    }
                    variant.colocated_variants = colocated;
                }
            }
        }

        if variant.transcript_consequences.is_empty() && variant.most_severe_consequence.is_none() {
            variant.most_severe_consequence = Some(Consequence::IntergenicVariant);
        }
    };

    if use_parallel {
        variants.par_iter_mut().for_each(annotate_one);
    } else {
        variants.iter_mut().for_each(annotate_one);
    }
}

fn build_report(
    args: &Args,
    variant_count: usize,
    cache_load_ms: u64,
    parse_ms: u64,
    plugin_names: &[String],
    runs: &[RunTimings],
) -> BenchReport {
    let csq_times: Vec<u64> = runs.iter().map(|r| r.consequence_ms).collect();
    let prefetch_times: Vec<u64> = runs.iter().map(|r| r.plugin_prefetch_ms).collect();
    let annotate_times: Vec<u64> = runs.iter().map(|r| r.plugin_annotate_ms).collect();
    let plugin_total_times: Vec<u64> = runs.iter().map(|r| r.plugin_total_ms).collect();
    let total_times: Vec<u64> = runs.iter().map(|r| r.total_ms).collect();

    let phases = PhaseTimings {
        cache_load_ms,
        parse_ms,
        consequence_ms: aggregate(&csq_times),
        plugin_prefetch_ms: aggregate(&prefetch_times),
        plugin_annotate_ms: aggregate(&annotate_times),
        plugin_total_ms: aggregate(&plugin_total_times),
        total_ms: aggregate(&total_times),
    };

    let mut plugin_detail: HashMap<String, PluginDetail> = HashMap::new();
    if !plugin_names.is_empty() {
        for name in plugin_names {
            let prefetch_vals: Vec<f64> = runs
                .iter()
                .filter_map(|r| {
                    r.plugin_timings.as_ref().and_then(|t| {
                        t.prefetch_ms
                            .iter()
                            .find(|(n, _)| n == name)
                            .map(|(_, ms)| *ms as f64)
                    })
                })
                .collect();
            let annotate_vals: Vec<f64> = runs
                .iter()
                .filter_map(|r| {
                    r.plugin_timings.as_ref().and_then(|t| {
                        t.annotate_ms
                            .iter()
                            .find(|(n, _)| n == name)
                            .map(|(_, ms)| *ms as f64)
                    })
                })
                .collect();

            plugin_detail.insert(
                name.clone(),
                PluginDetail {
                    prefetch_ms: aggregate_f64(&prefetch_vals),
                    annotate_ms: aggregate_f64(&annotate_vals),
                },
            );
        }
    }

    let mean_total = phases.total_ms.mean;
    let mean_csq = phases.consequence_ms.mean;
    let mean_plugin_total = phases.plugin_total_ms.mean;

    let variants_per_second = if mean_total > 0.0 {
        (variant_count as f64) / (mean_total / 1000.0)
    } else {
        0.0
    };

    let plugin_overhead_pct = if mean_csq > 0.0 {
        (mean_plugin_total / mean_csq) * 100.0
    } else {
        0.0
    };

    BenchReport {
        config: BenchConfig {
            input: args.input.clone(),
            variants_loaded: variant_count,
            plugins: plugin_names.to_vec(),
            buffer_size: args.buffer_size,
            threads: args.threads,
            runs: args.runs,
            mode: if args.sequential {
                "sequential".to_string()
            } else {
                "parallel".to_string()
            },
            label: args.label.clone(),
        },
        phases,
        plugin_detail,
        runs: runs.to_vec(),
        derived: DerivedMetrics {
            variants_per_second,
            plugin_overhead_pct,
            plugin_overhead_ms: mean_plugin_total,
            consequence_only_variants_per_second: if mean_csq > 0.0 {
                Some((variant_count as f64) / (mean_csq / 1000.0))
            } else {
                None
            },
        },
    }
}

fn aggregate(values: &[u64]) -> AggregatedTiming {
    if values.is_empty() {
        return AggregatedTiming::default();
    }
    let sum: u64 = values.iter().sum();
    AggregatedTiming {
        mean: sum as f64 / values.len() as f64,
        min: *values.iter().min().unwrap(),
        max: *values.iter().max().unwrap(),
    }
}

fn aggregate_f64(values: &[f64]) -> AggregatedTimingF64 {
    if values.is_empty() {
        return AggregatedTimingF64::default();
    }
    let sum: f64 = values.iter().sum();
    AggregatedTimingF64 {
        mean: sum / values.len() as f64,
        min: values.iter().cloned().fold(f64::INFINITY, f64::min),
        max: values.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    }
}
