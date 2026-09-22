// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Pipeline runner that orchestrates input parsing, annotation, and output.
//!
//! Reads VCF variants, loads transcripts from a converted JSON cache,
//! computes consequence annotations using the vep-effects engine,
//! and writes annotated output. The stages run as the three-thread pipeline
//! in [`crate::pipeline`]; this module owns the setup and `annotate_batch`.
//!
//! Annotation of each batch is parallelised across the resolved `--fork` count
//! (default: every logical CPU; `--fork 1` runs serially) using rayon's `par_iter_mut`.
//!
//! Perl citations name modules of ensembl-vep release/115 (`Bio/EnsEMBL/VEP/...`).

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{bail, Context};
use noodles::bgzf;
use rayon::prelude::*;
use rayon::ThreadPool;
use tracing::{debug, info};

use crate::annotator::{build_thread_pool, LoadedPlugins};
use crate::config::Config;
use crate::json_cache;
use crate::pick::FilterConfig;
use crate::pipeline;
use crate::transcript_index::{LazyTranscriptIndexes, TranscriptIndex};
use crate::variation_matcher;
use vep_builtins::BuiltinRegistry;
use vep_core::consequence::{Consequence, Impact};
use vep_core::variant::{InputVariant, VariantClass};
use vep_io::output::fields::{ExtraFieldsPlan, FieldOptions};
use vep_io::output::json::JsonOutputFormatter;
use vep_io::output::parquet::ParquetOutputFormatter;
use vep_io::output::tab::TabOutputFormatter;
use vep_io::output::vcf_output::VcfOutputFormatter;
use vep_io::OutputFormatter;
use vep_plugin::PluginRegistry;

/// Read buffer over the input file, under the bgzf block reader when there is
/// one.
const INPUT_BUFFER_BYTES: usize = 1 << 20;

/// Write buffer over the output file: annotated rows come out in batches of
/// several megabytes, so the default 8 KiB would flush hundreds of times per batch.
const OUTPUT_BUFFER_BYTES: usize = 4 << 20;

/// Atomic counters for pipeline statistics.
pub(crate) struct PipelineStats {
    variants_processed: AtomicU64,
    consequences_calculated: AtomicU64,
}

/// Pipeline runner for VEP.
pub struct Runner {
    config: Config,
    /// The cache's `info.json`, when the JSON cache carries one; it supplies the
    /// header's `## <source> version` lines and the JSON `assembly_name`.
    cache_info: Option<vep_core::cache_info::CacheInfo>,
    /// The invocation as VEP prints it in `## VEP command-line:`.
    command_line: String,
}

/// Renders the invocation as VEP prints it in `## VEP command-line:`
/// (Config.pm `full_command`): `vep` then every option sorted by name, a value
/// after its flag, and each path collapsed to `[PATH]/<basename>`. Short flags
/// print under their long names.
fn vep_command_line(args: impl Iterator<Item = String>) -> String {
    let args: Vec<String> = args.collect();
    let mut options: Vec<(String, Option<String>)> = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        let (name, inline_value) = match arg.split_once('=') {
            Some((n, v)) if n.starts_with('-') => (n.to_string(), Some(v.to_string())),
            _ => (arg.clone(), None),
        };
        if !name.starts_with('-') {
            i += 1;
            continue;
        }
        let long = match name.as_str() {
            "-i" => "input_file",
            "-o" => "output_file",
            "-a" => "assembly",
            "--force" => "force_overwrite",
            other => other.trim_start_matches('-'),
        }
        .to_string();
        let value = match inline_value {
            Some(v) => Some(v),
            None if i + 1 < args.len() && !args[i + 1].starts_with('-') => {
                i += 1;
                Some(args[i].clone())
            }
            None => None,
        };
        options.push((long, value.map(|v| collapse_paths(&v))));
        i += 1;
    }
    options.sort_by(|a, b| a.0.cmp(&b.0));
    let mut out = String::from("vep");
    for (flag, value) in options {
        out.push_str(" --");
        out.push_str(&flag);
        if let Some(v) = value {
            out.push(' ');
            out.push_str(&v);
        }
    }
    out
}

/// `[PATH]/<basename>` for each path in a value (`,`-separated lists included).
fn collapse_paths(value: &str) -> String {
    value
        .split(',')
        .map(|part| match part.rfind('/') {
            Some(idx) if idx + 1 < part.len() => format!("[PATH]/{}", &part[idx + 1..]),
            _ => part.to_string(),
        })
        .collect::<Vec<_>>()
        .join(",")
}

/// VEP's `--max_sv_size` check (Parser.pm `validate_vf`): a structural variant
/// whose span exceeds the limit keeps its VCF line without consequences and
/// is absent from the JSON output.
pub(crate) fn mark_oversize_sv(variant: &mut InputVariant, max_sv_size: u64) {
    if variant.is_structural && variant.end.saturating_sub(variant.start) > max_sv_size {
        variant.oversize_sv = true;
    }
}

impl Runner {
    pub fn new(config: Config) -> Self {
        let cache_info = config
            .json_cache
            .as_deref()
            .and_then(|dir| std::fs::read_to_string(Path::new(dir).join("info.json")).ok())
            .and_then(|text| serde_json::from_str(&text).ok());
        let command_line = vep_command_line(std::env::args().skip(1));
        Self {
            config,
            cache_info,
            command_line,
        }
    }

    /// Which optional output fields the configuration switches on, shared by
    /// every output format.
    fn field_options(&self) -> FieldOptions {
        let c = &self.config;
        FieldOptions {
            allele_number: c.allele_number,
            show_ref_allele: c.show_ref_allele,
            uploaded_allele: false,
            include_pick: c.flag_pick
                || c.flag_pick_allele
                || c.flag_pick_allele_gene
                || c.pick
                || c.pick_allele
                || c.pick_allele_gene,
            variant_class: c.variant_class,
            minimal: c.minimal,
            symbol: c.symbol,
            biotype: c.biotype,
            canonical: c.canonical,
            mane_select: c.mane_select,
            mane: c.mane,
            tsl: c.tsl,
            gencode_primary: false,
            appris: c.appris,
            ccds: c.ccds,
            protein: c.protein,
            uniprot: c.uniprot,
            xref_refseq: c.xref_refseq,
            gene_phenotype: c.gene_phenotype,
            sift: c.sift.is_some(),
            polyphen: c.polyphen.is_some(),
            numbers: c.numbers,
            domains: c.domains,
            mirna: c.mirna,
            hgvs: c.hgvs,
            hgvsg: c.hgvsg,
            af: c.af,
            af_1kg: c.af_1kg,
            af_gnomade: c.af_gnomade,
            af_gnomadg: c.af_gnomadg,
            max_af: c.max_af,
            check_existing: c.check_existing,
            pubmed: c.pubmed,
            dont_skip: c.dont_skip,
            overlaps: false,
            regulatory: c.regulatory,
            no_escape: c.output_format == "json",
        }
    }

    /// The assembly VEP names in JSON output: the flag, else the cache's.
    /// The run's identity for the Parquet footer: version, assembly, cache and
    /// the command line as the headers print it.
    fn parquet_kv_metadata(&self) -> Vec<(String, String)> {
        let mut kv = vec![
            (
                "vep_rs_version".to_string(),
                env!("CARGO_PKG_VERSION").to_string(),
            ),
            (
                "vep_api_version".to_string(),
                vep_core::VEP_VERSION.to_string(),
            ),
            ("cache".to_string(), self.config.cache_dir_path()),
            ("command_line".to_string(), self.command_line.clone()),
        ];
        if let Some(asm) = self.assembly_name() {
            kv.push(("assembly".to_string(), asm));
        }
        if let Some(info) = &self.cache_info {
            for (k, v) in &info.source_versions {
                kv.push((format!("source_{k}"), v.clone()));
            }
        }
        kv
    }

    fn assembly_name(&self) -> Option<String> {
        self.config
            .assembly
            .clone()
            .or_else(|| self.cache_info.as_ref().map(|i| i.assembly.clone()))
            .filter(|a| !a.is_empty())
    }

    /// Execute the full VEP pipeline.
    pub fn run(&self) -> anyhow::Result<()> {
        let start = Instant::now();

        if self.config.output_file != "STDOUT" {
            let path = Path::new(&self.config.output_file);
            if path.exists() && !self.config.force_overwrite {
                bail!(
                    "Output file '{}' already exists. Use --force_overwrite to overwrite.",
                    self.config.output_file
                );
            }
        }

        info!("input: {}", self.config.input_file);
        info!("output: {}", self.config.output_file);
        info!("output format: {}", self.config.output_format);

        // A scoped pool, not the global one, so embedding callers can own their
        // own pool without contention.
        let num_threads = self.config.fork.max(1);
        let pool = build_thread_pool(num_threads)?;
        if num_threads > 1 {
            info!("Using {} threads for annotation", num_threads);
        }

        if self.config.shift_3prime && self.config.fasta.is_none() {
            bail!("--shift_3prime requires --fasta <path> (indexed FASTA: .fa + .fai)");
        }

        let reference_fasta = if let Some(ref fasta_path) = self.config.fasta {
            info!("Using reference FASTA: {}", fasta_path);
            Some(Arc::new(
                vep_fasta::IndexedFasta::from_path(fasta_path)
                    .with_context(|| format!("Failed to open indexed FASTA '{fasta_path}'"))?,
            ))
        } else {
            None
        };

        let needs_loftee_context = self.config.plugin.iter().any(|arg| {
            arg.split(',')
                .next()
                .is_some_and(|name| name.eq_ignore_ascii_case("LoFTEE"))
        });

        // Sequence-aware plugins (LoFTEE) share the FASTA Arc, so the FASTA is
        // not opened twice.
        let plugin_fasta = reference_fasta.clone();

        let effects_config = Arc::new(vep_effects::EffectsConfig {
            upstream_distance: self.config.distance.0,
            downstream_distance: self.config.distance.1,
            reference_fasta,
            enable_indel_3prime_shift: self.config.shift_3prime,
            // This primarily needs to cover simple repeat-induced shifting around splice windows.
            max_indel_3prime_shift: 2000,
            populate_loftee_context: needs_loftee_context,
            compute_hgvs: needs_hgvs_computation(&self.config),
            compute_exon_intron_numbers: needs_exon_intron_numbers(&self.config),
        });

        let LoadedPlugins { builtin, dylib } = self.load_plugins(plugin_fasta)?;
        let builtin_plugins = Arc::new(builtin);
        let dylib_plugins = Arc::new(dylib);
        if builtin_plugins.plugin_count() > 0 {
            info!(
                "Loaded {} built-in plugin(s)",
                builtin_plugins.plugin_count()
            );
        }
        if dylib_plugins.plugin_count() > 0 {
            info!("Loaded {} dylib plugin(s)", dylib_plugins.plugin_count());
        }

        // Indexed lazily per chromosome: the key set comes from a `read_dir`
        // walk, so the Perl-parity presence predicates are complete immediately,
        // while a chromosome's JSON parses only when a variant queries it.
        let transcripts: Arc<LazyTranscriptIndexes> =
            if let Some(ref json_cache_path) = self.config.json_cache {
                let index_impl = self.config.transcript_index_impl;
                info!(
                    "Transcript cache: {} (lazy per-chromosome, index impl={})",
                    json_cache_path,
                    index_impl.name()
                );
                Arc::new(LazyTranscriptIndexes::new(
                    json_cache_path,
                    index_impl,
                    effects_config.upstream_distance,
                    effects_config.downstream_distance,
                )?)
            } else {
                Arc::new(LazyTranscriptIndexes::empty())
            };

        let variation_data: Option<Arc<json_cache::VariationData>> = if self.config.check_existing
            || self.config.af
            || self.config.af_1kg
            || self.config.af_gnomade
            || self.config.af_gnomadg
            || self.config.max_af
            || self.config.pubmed
        {
            if let Some(ref json_cache_path) = self.config.json_cache {
                info!("Loading variations from JSON cache: {}", json_cache_path);
                let data = json_cache::load_all_variations(json_cache_path)?;
                Some(Arc::new(data))
            } else {
                None
            }
        } else {
            None
        };

        // `MultithreadedReader` when decompression_threads > 1: block
        // decompression is the largest single wall-time lever on .vcf.gz input.
        // Both bgzf readers are already `BufRead` over their decompressed block,
        // so they take no `BufReader` of their own; the buffer sits under them,
        // where the block reader would otherwise issue two small reads per
        // compressed block.
        let reader: Box<dyn BufRead + Send> = if self.config.input_file == "STDIN" {
            Box::new(BufReader::with_capacity(INPUT_BUFFER_BYTES, io::stdin()))
        } else {
            let f = File::open(&self.config.input_file).with_context(|| {
                format!("Failed to open input file '{}'", self.config.input_file)
            })?;
            let f = BufReader::with_capacity(INPUT_BUFFER_BYTES, f);
            if self.config.input_file.ends_with(".gz") {
                match std::num::NonZeroUsize::new(self.config.decompression_threads) {
                    Some(n) if n.get() > 1 => {
                        Box::new(bgzf::MultithreadedReader::with_worker_count(n, f))
                    }
                    _ => Box::new(bgzf::Reader::new(f)),
                }
            } else {
                Box::new(f)
            }
        };

        // Parquet: a TSV intermediate at `<output_file>.tsv.tmp`, converted by a
        // DuckDB subprocess after the write loop.
        let parquet_intermediate_path: Option<String> =
            if self.config.output_format == "parquet" && self.config.output_file != "STDOUT" {
                // A run that cannot be finalized must fail before annotating.
                check_duckdb()?;
                Some(format!("{}.tsv.tmp", self.config.output_file))
            } else {
                None
            };

        let mut writer: Box<dyn Write + Send> = if self.config.output_file == "STDOUT" {
            Box::new(BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, io::stdout()))
        } else if let Some(ref intermediate) = parquet_intermediate_path {
            let f = File::create(intermediate).with_context(|| {
                format!(
                    "Failed to create Parquet intermediate TSV '{}'",
                    intermediate
                )
            })?;
            Box::new(BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, f))
        } else {
            let f = File::create(&self.config.output_file).with_context(|| {
                format!("Failed to create output file '{}'", self.config.output_file)
            })?;
            Box::new(BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, f))
        };

        let format = self.detect_format(&self.config.format);
        debug!("detected input format: {}", format);

        let filter_config = FilterConfig {
            pick: self.config.pick,
            pick_allele: self.config.pick_allele,
            per_gene: self.config.per_gene,
            pick_allele_gene: self.config.pick_allele_gene,
            most_severe: self.config.most_severe,
            summary: self.config.summary,
            flag_pick: self.config.flag_pick,
            flag_pick_allele: self.config.flag_pick_allele,
            flag_pick_allele_gene: self.config.flag_pick_allele_gene,
        };
        let filters_active = filter_config.any_active();

        let output_format = self.config.output_format.as_str();
        let mut vcf_formatter: Option<VcfOutputFormatter> = None;
        let mut json_formatter: Option<JsonOutputFormatter> = None;
        let mut parquet_formatter: Option<ParquetOutputFormatter> = None;
        let mut tab_formatter: Option<TabOutputFormatter> = None;
        let plugin_fields = builtin_plugins
            .all_header_info()
            .into_iter()
            .map(|(name, _)| name)
            .collect::<Vec<_>>();
        let options = self.field_options();
        let extra_plan = ExtraFieldsPlan::new(options.clone());

        match output_format {
            "vcf" => {
                vcf_formatter = Some(
                    VcfOutputFormatter::new(
                        self.config.vcf_info_field.clone(),
                        options,
                        plugin_fields,
                    )
                    .with_run_info(
                        &self.config.cache_dir_path(),
                        self.cache_info.as_ref(),
                        &self.command_line,
                    ),
                );
            }
            "json" => {
                json_formatter = Some(JsonOutputFormatter::new(
                    options,
                    plugin_fields,
                    self.assembly_name(),
                ));
            }
            "tab" => {
                let header = vep_io::output::fields::tab_header(
                    &options,
                    &plugin_fields,
                    &self.config.cache_dir_path(),
                    self.cache_info.as_ref(),
                    &self.command_line,
                );
                tab_formatter = Some(TabOutputFormatter::new(
                    options,
                    plugin_fields,
                    header,
                    self.config.no_headers,
                ));
            }
            "parquet" => {
                parquet_formatter = Some(ParquetOutputFormatter::new(options, plugin_fields));
            }
            _ => {
                // The default VEP text format is written by `write_record`.
            }
        }

        // The Parquet TSV intermediate is read by DuckDB with `header=true`, so
        // its header row is mandatory regardless of `--no_headers`, a
        // VEP-default-tab flag.
        if output_format == "parquet" {
            if let Some(ref mut fmt) = parquet_formatter {
                fmt.write_header(&mut writer)?;
            }
        } else if !self.config.no_headers {
            match output_format {
                "vcf" => {
                    // The header block is written once the input's own header
                    // lines are known: at its #CHROM line or first data line.
                }
                "json" => {
                    // JSON has no header
                }
                "tab" => {
                    if let Some(ref mut fmt) = tab_formatter {
                        fmt.write_header(&mut writer)?;
                    }
                }
                _ => {
                    self.write_header(&mut writer)?;
                }
            }
        }

        let stats = PipelineStats {
            variants_processed: AtomicU64::new(0),
            consequences_calculated: AtomicU64::new(0),
        };

        let batch_size = self.config.buffer_size;
        let use_parallel = num_threads > 1;
        let stats_enabled = !self.config.quiet;
        let capture_raw_input = needs_raw_input_capture(&self.config);

        let annotation_resources = AnnotationResources {
            transcripts,
            variation_data,
            effects_config,
            builtin_plugins,
            dylib_plugins,
            prediction_config: PredictionConfig {
                sift: self.config.sift.clone(),
                polyphen: self.config.polyphen.clone(),
            },
        };

        // VCF is read line by line because header passthrough needs the lines;
        // every other format goes through its parser.
        let input = match format {
            "vcf" => pipeline::InputSource::Vcf(reader),
            _ => {
                let parser: Box<dyn vep_io::InputParser> = match format {
                    "ensembl" => Box::new(vep_io::input::ensembl::EnsemblParser::new(reader)),
                    "region" => Box::new(vep_io::input::region::RegionParser::new(reader)),
                    "hgvs" => Box::new(vep_io::input::hgvs::HgvsParser::new(reader)),
                    other => {
                        bail!(
                            "Input format '{}' is not supported. Supported formats: vcf, ensembl, region, hgvs.",
                            other
                        );
                    }
                };
                if output_format == "vcf" && !self.config.no_headers {
                    // A non-VCF input has no header lines of its own.
                    if let Some(ref fmt) = vcf_formatter {
                        for h in fmt.header_lines(&[]) {
                            writeln!(writer, "{}", h)?;
                        }
                    }
                }
                pipeline::InputSource::Parser(parser)
            }
        };

        let prewarm_order = if self.config.input_file.ends_with(".gz") {
            pipeline::tabix_contig_order(&self.config.input_file)
        } else {
            Vec::new()
        };
        let reader_cfg = pipeline::ReaderConfig {
            batch_size,
            capture_raw_input,
            allow_non_variant: self.config.allow_non_variant,
            dont_skip: self.config.dont_skip,
            max_sv_size: self.config.max_sv_size,
            output_format,
            no_headers: self.config.no_headers,
            vcf_formatter: vcf_formatter.as_ref(),
            transcripts: &annotation_resources.transcripts,
            prewarm_order: &prewarm_order,
        };
        let write_ctx = pipeline::WriteContext {
            filter_config: &filter_config,
            filters_active,
            output_format,
            vcf_formatter: vcf_formatter.as_ref(),
            json_formatter: json_formatter.as_ref(),
            parquet_formatter: parquet_formatter.as_ref(),
            tab_formatter: tab_formatter.as_ref(),
            extra_plan: &extra_plan,
        };
        let coordinator = pipeline::Coordinator {
            pool: &pool,
            use_parallel,
            resources: &annotation_resources,
            stats: &stats,
            stats_enabled,
            write_ctx: &write_ctx,
            capture_raw_input,
            dont_skip: self.config.dont_skip,
            max_sv_size: self.config.max_sv_size,
            batch_size,
            progress: !self.config.quiet,
        };

        // The writer thread flushes and closes the output before this returns.
        let variant_count = pipeline::run(input, &reader_cfg, &coordinator, writer)?;

        if !self.config.quiet && variant_count >= 10000 {
            eprint!("\r");
        }

        if let Some(intermediate) = parquet_intermediate_path {
            finalize_parquet_output(
                &intermediate,
                &self.config.output_file,
                &self.config.parquet_shape,
                parquet_formatter
                    .as_ref()
                    .map(|f| f.header_columns())
                    .unwrap_or_default(),
                &self.parquet_kv_metadata(),
                self.config.parquet_row_group_size,
            )?;
            if !self.config.keep_intermediate {
                std::fs::remove_file(&intermediate).with_context(|| {
                    format!(
                        "Failed to remove Parquet intermediate '{}' after DuckDB success",
                        intermediate
                    )
                })?;
            }
        }

        let elapsed = start.elapsed();
        if !self.config.quiet {
            let conseqs = stats.consequences_calculated.load(Ordering::Relaxed);
            eprintln!(
                "Processed {} variant{} ({} consequence{}) in {:.1}s{}",
                variant_count,
                if variant_count == 1 { "" } else { "s" },
                conseqs,
                if conseqs == 1 { "" } else { "s" },
                elapsed.as_secs_f64(),
                if num_threads > 1 {
                    format!(" using {} threads", num_threads)
                } else {
                    String::new()
                },
            );
        }

        Ok(())
    }

    /// Load the `--plugin` CLI arguments: each as a built-in when one claims the
    /// name, otherwise as a C-ABI dylib resolved by path or under `--dir_plugins`.
    /// A name that resolves to neither is warned about and skipped, which is
    /// Ensembl VEP's default without `--safe`.
    ///
    /// The resolved `--assembly` is passed through so assembly-sensitive plugins
    /// (REVEL's dual coordinate columns, LoFTEE's per-assembly GERP format) can
    /// select the right column and reader, and so plugins carrying an assembly
    /// marker in their data file can reject a mismatched one.
    fn load_plugins(
        &self,
        reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    ) -> anyhow::Result<LoadedPlugins> {
        crate::annotator::load_plugins(
            &self.config.plugin,
            self.config.dir_plugins.as_deref(),
            reference_fasta,
            crate::annotator::resolve_assembly(self.config.assembly.as_deref()),
        )
    }

    /// Detect the input format from the `--format` flag or the input file's extension.
    fn detect_format<'a>(&self, format: &'a str) -> &'a str {
        if format != "guess" {
            return format;
        }
        let input = &self.config.input_file;
        if input.ends_with(".vcf") || input.ends_with(".vcf.gz") || input.ends_with(".bcf") {
            return "vcf";
        }
        if input.ends_with(".hgvs") {
            return "hgvs";
        }
        "vcf"
    }

    /// Write the default-format header lines VEP writes, then the column line.
    fn write_header(&self, writer: &mut dyn Write) -> anyhow::Result<()> {
        let options = self.field_options();
        for line in vep_io::output::fields::default_format_header(
            &options,
            &self.config.cache_dir_path(),
            self.cache_info.as_ref(),
            &self.command_line,
        ) {
            writeln!(writer, "{}", line)?;
        }
        Ok(())
    }
}

/// The oldest DuckDB release whose `COPY ... (FORMAT PARQUET)` accepts every
/// option the Parquet writer uses (`PARQUET_VERSION`, `WRITE_BLOOM_FILTER`,
/// `BLOOM_FILTER_FALSE_POSITIVE_RATIO`, `DICTIONARY_SIZE_LIMIT`,
/// `DICTIONARY_COMPRESSION_RATIO_THRESHOLD`, `WRITE_PARTITION_COLUMNS`,
/// `KV_METADATA`).
const MIN_DUCKDB_VERSION: (u32, u32) = (1, 2);

/// Parses `duckdb --version` output such as `v1.5.5 (Variegata) d8cdaa33fd`.
fn parse_duckdb_version(text: &str) -> Option<(u32, u32, u32)> {
    let token = text.split_whitespace().next()?.trim_start_matches('v');
    let mut parts = token.split('.').map(|p| p.parse::<u32>().ok());
    Some((
        parts.next()??,
        parts.next()??,
        parts.next().flatten().unwrap_or(0),
    ))
}

/// Checks that a `duckdb` CLI new enough for the Parquet writer is on PATH.
fn check_duckdb() -> anyhow::Result<()> {
    let output = std::process::Command::new("duckdb")
        .arg("--version")
        .output()
        .ok()
        .filter(|o| o.status.success())
        .ok_or_else(|| {
            anyhow::anyhow!(
                "--output_format=parquet requires the `duckdb` CLI on PATH. \
                 Install via `brew install duckdb` (macOS) or the release binary \
                 at https://duckdb.org/docs/installation/"
            )
        })?;
    let text = String::from_utf8_lossy(&output.stdout);
    match parse_duckdb_version(&text) {
        Some((major, minor, _)) if (major, minor) >= MIN_DUCKDB_VERSION => Ok(()),
        Some((major, minor, patch)) => anyhow::bail!(
            "--output_format=parquet needs duckdb {}.{} or newer; found {}.{}.{}",
            MIN_DUCKDB_VERSION.0,
            MIN_DUCKDB_VERSION.1,
            major,
            minor,
            patch
        ),
        None => anyhow::bail!("could not parse `duckdb --version` output: {}", text.trim()),
    }
}

/// DuckDB writes row groups in whole vectors of this many rows, so a requested
/// row group size is rounded up to a multiple of it. The dictionary size limit
/// must equal the rounded size: a column with a distinct value per row in a
/// full row group otherwise exceeds the limit and loses its dictionary and,
/// with it, its Bloom filter.
const DUCKDB_VECTOR_SIZE: u64 = 2_048;

/// Rows per Parquet row group as DuckDB will write them.
fn parquet_row_group_rows(requested: u64) -> u64 {
    requested.max(1).div_ceil(DUCKDB_VECTOR_SIZE) * DUCKDB_VECTOR_SIZE
}

/// Key columns that lead the flat schema and key the nested shape's groups.
const PARQUET_KEY_COLUMNS: [&str; 5] = ["chrom", "pos", "end", "ref", "alt"];

/// Columns stored as integers rather than text; everything else stays VARCHAR.
fn parquet_column_sql(name: &str) -> String {
    match name {
        "pos" | "end" => format!("TRY_CAST({0} AS BIGINT) AS {0}", sql_ident(name)),
        "DISTANCE" => format!("TRY_CAST({0} AS INTEGER) AS {0}", sql_ident(name)),
        "STRAND" => format!("TRY_CAST({0} AS TINYINT) AS {0}", sql_ident(name)),
        _ => sql_ident(name),
    }
}

/// Columns of the nested shape that describe the input record rather than one
/// consequence; they stay scalar and join the grouping key.
const PARQUET_PER_VARIANT_COLUMNS: [&str; 4] = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Existing_variation",
];

/// Subprocesses the `duckdb` CLI to turn the TSV intermediate into a
/// Hive-partitioned Parquet directory (`<out>/chrom=<c>/*.parquet`).
///
/// Schema: the five key columns then exactly the `--tab` columns
/// (`ParquetOutputFormatter`); `pos`/`end` are BIGINT, `DISTANCE` INTEGER,
/// `STRAND` TINYINT, everything else VARCHAR with `-` stored as NULL. The nested
/// shape groups rows by (chrom, pos, ref, alt) and stores the per-consequence
/// columns as LIST, in the row order of the flat shape.
///
/// Writer options favour the reader: Parquet v2 pages, ZSTD level 9, row groups
/// of `row_group_size` rows (rounded up to whole DuckDB vectors) sorted by
/// (chrom, pos, ref, alt) so the `pos` statistics are tight, a dictionary on
/// every column (the size limit equals the row group and the
/// compression-ratio threshold is off) because
/// DuckDB writes a Bloom filter only for dictionary-encoded columns, Bloom
/// filters at a 0.1% false-positive rate, the partition column kept inside the
/// files, and the run's identity in the footer's key-value metadata. LIST
/// columns cannot carry Bloom filters; the scalar key columns still do.
fn finalize_parquet_output(
    tsv_path: &str,
    out_path: &str,
    shape: &str,
    header_columns: Vec<&str>,
    kv_metadata: &[(String, String)],
    row_group_size: u64,
) -> anyhow::Result<()> {
    use std::process::Command;

    let tab_columns: Vec<&str> = header_columns
        .iter()
        .skip(PARQUET_KEY_COLUMNS.len())
        .copied()
        .collect();

    // `all_varchar` keeps every value textual until the explicit casts below:
    // DuckDB's type sniffing otherwise reads `cDNA_position` as BIGINT from early
    // rows and fails on the first `204-205`. The intermediate carries no quoting.
    let source = format!(
        "read_csv('{}', delim='\\t', header=true, all_varchar=true, quote='', escape='')",
        escape_sql_string(tsv_path)
    );
    let select_sql = match shape {
        "flat" => {
            let cols = header_columns
                .iter()
                .map(|c| parquet_column_sql(c))
                .collect::<Vec<_>>()
                .join(", ");
            format!("SELECT {cols} FROM {source} ORDER BY chrom, pos, ref, alt")
        }
        _ => {
            // Groups are keyed by the variant key and the per-variant columns,
            // so two input records at one site stay two rows.
            let mut cols: Vec<String> = vec![
                "chrom".to_string(),
                "TRY_CAST(pos AS BIGINT) AS pos".to_string(),
                "TRY_CAST(\"end\" AS BIGINT) AS \"end\"".to_string(),
                "ref".to_string(),
                "alt".to_string(),
            ];
            let mut group_by: Vec<String> = vec![
                "chrom".to_string(),
                "TRY_CAST(pos AS BIGINT)".to_string(),
                "TRY_CAST(\"end\" AS BIGINT)".to_string(),
                "ref".to_string(),
                "alt".to_string(),
            ];
            for c in &tab_columns {
                if PARQUET_PER_VARIANT_COLUMNS.contains(c) {
                    cols.push(sql_ident(c));
                    group_by.push(sql_ident(c));
                } else {
                    // `list` keeps the flat row order because the source is
                    // ordered before grouping.
                    let inner = match *c {
                        "DISTANCE" => format!("TRY_CAST({} AS INTEGER)", sql_ident(c)),
                        "STRAND" => format!("TRY_CAST({} AS TINYINT)", sql_ident(c)),
                        _ => sql_ident(c),
                    };
                    cols.push(format!(
                        "list({inner} ORDER BY row_order) AS {}",
                        sql_ident(c)
                    ));
                }
            }
            format!(
                "SELECT {} FROM (SELECT *, row_number() OVER () AS row_order FROM {source}) \
                 GROUP BY {} ORDER BY chrom, pos, ref, alt",
                cols.join(", "),
                group_by.join(", ")
            )
        }
    };

    let kv = kv_metadata
        .iter()
        .map(|(k, v)| format!("'{}': '{}'", escape_sql_string(k), escape_sql_string(v)))
        .collect::<Vec<_>>()
        .join(", ");
    let sql = format!(
        "COPY ({select_sql}) TO '{}' (FORMAT PARQUET, PARQUET_VERSION V2, COMPRESSION 'zstd', \
         COMPRESSION_LEVEL 9, ROW_GROUP_SIZE {rg}, DICTIONARY_SIZE_LIMIT {rg}, \
         DICTIONARY_COMPRESSION_RATIO_THRESHOLD -1, WRITE_BLOOM_FILTER true, \
         BLOOM_FILTER_FALSE_POSITIVE_RATIO 0.001, PARTITION_BY (chrom), \
         WRITE_PARTITION_COLUMNS true, KV_METADATA {{{kv}}}, OVERWRITE_OR_IGNORE true);",
        escape_sql_string(out_path),
        rg = parquet_row_group_rows(row_group_size),
    );

    info!(
        "Parquet finalize: duckdb from '{}' -> '{}' (shape={})",
        tsv_path, out_path, shape,
    );

    let output = Command::new("duckdb")
        .arg("-c")
        .arg(&sql)
        .output()
        .context("Failed to spawn `duckdb` for Parquet finalize")?;
    if !output.status.success() {
        anyhow::bail!(
            "duckdb COPY to Parquet failed (exit {}): stderr=\n{}\nsql=\n{}",
            output.status,
            String::from_utf8_lossy(&output.stderr),
            sql,
        );
    }
    Ok(())
}

/// Escape a string for embedding in a single-quoted SQL literal.
/// Replaces `'` with `''` (standard SQL escape).
fn escape_sql_string(s: &str) -> String {
    s.replace('\'', "''")
}

/// Quote a column identifier for DuckDB (double-quotes, standard SQL).
/// VEP CSQ field names like "cDNA_position" need quoting because SQL
/// identifiers are typically lowercased/case-insensitive otherwise.
fn sql_ident(name: &str) -> String {
    format!("\"{}\"", name.replace('"', "\"\""))
}

/// Shared annotation resources passed to each batch.
pub(crate) struct AnnotationResources {
    transcripts: Arc<LazyTranscriptIndexes>,
    variation_data: Option<Arc<json_cache::VariationData>>,
    effects_config: Arc<vep_effects::EffectsConfig>,
    builtin_plugins: Arc<BuiltinRegistry>,
    dylib_plugins: Arc<PluginRegistry>,
    prediction_config: PredictionConfig,
}

/// Whether HGVSc/HGVSp must be computed for this run: only under `--hgvs`.
/// VEP fills the HGVS fields of every format, the VCF CSQ columns included, only
/// when the flag is set, so without it the columns stay empty and the computation
/// (which needs the reference sequence) is skipped.
pub(crate) fn needs_hgvs_computation(config: &Config) -> bool {
    config.hgvs
}

/// Whether the `EXON`/`INTRON` ordinals must be computed for this run: under
/// `--numbers` (which `--vcf` and `--everything` imply), the only flag any
/// output format prints them for, and whenever a plugin is loaded, since LoFTEE
/// reads them from the consequence and a dylib plugin receives the whole
/// consequence as JSON.
pub(crate) fn needs_exon_intron_numbers(config: &Config) -> bool {
    config.numbers || !config.plugin.is_empty()
}

/// Output formats that read `InputVariant::raw_input`.
///
/// The VCF writer re-emits each original record with a CSQ INFO field appended,
/// and the JSON writer carries the line as VEP's `input` key. Every other format
/// builds its rows from the variant's own fields. The VCF writer also has a
/// fallback that reconstructs a minimal record from `original_chr`/`start`/`id`/
/// `ref_allele`/`alt_alleles` when `raw_input` is `None`, so a missing capture
/// degrades silently to a rebuilt line rather than erroring, which is why the
/// capture must be gated on the format list here rather than on a heuristic.
const RAW_INPUT_OUTPUT_FORMATS: [&str; 2] = ["vcf", "json"];

/// Whether the parser must retain the source VCF line on each variant.
///
/// See [`RAW_INPUT_OUTPUT_FORMATS`] and `vcf_parser::parse_vcf_line`.
pub(crate) fn needs_raw_input_capture(config: &Config) -> bool {
    RAW_INPUT_OUTPUT_FORMATS
        .iter()
        .any(|fmt| config.output_format == *fmt)
}

fn structural_variant_bounds(variant: &InputVariant) -> (u64, u64) {
    let variant_start = variant.start.min(variant.end);
    let variant_end = variant
        .sv_end
        .unwrap_or(variant.end)
        .max(variant.start.max(variant.end));
    (variant_start, variant_end)
}

/// Annotate a batch of variants, optionally in parallel using rayon.
///
/// When `use_parallel` is true, uses `par_iter_mut` to distribute annotation
/// across the rayon thread pool. The shared data (transcripts, variation_data,
/// effects_config) is behind `Arc` and only read, so no locking is needed.
///
/// After consequence assignment, built-in plugins are invoked via `run_batch`
/// to annotate the entire buffer in one pass (enabling batched tabix queries).
pub(crate) fn annotate_batch(
    batch: &mut Vec<InputVariant>,
    resources: &AnnotationResources,
    pool: &ThreadPool,
    use_parallel: bool,
    stats: &PipelineStats,
    stats_enabled: bool,
) -> anyhow::Result<()> {
    let batch_len = batch.len();
    let annotation_phase_start = Instant::now();

    // Materialize every chromosome this batch touches before fanning out.
    //
    // Three reasons this must not happen inside the parallel closure:
    //  1. `load_transcripts_for_chr` parses its shards with rayon. Calling it
    //     from a worker already inside `pool.install(...)` blocks that worker on
    //     a job that needs the pool it is occupying: a deadlock.
    //  2. Even without the deadlock, N workers racing the same cold chromosome
    //     each parse it and all but one discard the result.
    //  3. A parse failure has no error channel from inside the per-variant
    //     closure. Here it aborts the run. Otherwise a truncated shard would annotate every variant
    //     on that chromosome as `intergenic_variant` and still exit 0, because
    //     `contains_key` answers from the directory listing and so still reports
    //     the chromosome as present.
    //
    // Pre-warming serially here means each chromosome is parsed exactly once,
    // outside the pool, and every `get` in the closure below is a pure lookup.
    for variant in batch.iter() {
        resources.transcripts.prewarm(&variant.chr)?;
        if let Some(mate_chr) = variant.mate_chr.as_deref() {
            resources.transcripts.prewarm(mate_chr)?;
        }
    }

    let annotate_one = |variant: &mut InputVariant| -> u64 {
        // Annotate with transcript consequences. Every variant, a breakend spanning
        // more than VEP's 10 Mb `--max_sv_size` included, is annotated against the
        // chromosome's complete transcript set, so the output of one record never
        // depends on which other records share its batch. VEP marks such breakends
        // skipped at cache-region load and then annotates them against whatever
        // transcripts their batch-mates happened to load; that batch dependence is a
        // VEP defect and is not reproduced here.
        let is_paired_breakend = variant.variant_class == VariantClass::Translocation
            && !variant.is_single_breakend
            && variant.mate_chr.is_some()
            && variant.mate_pos.is_some();
        let is_interchrom_paired_breakend = is_paired_breakend
            && variant
                .mate_chr
                .as_deref()
                .map(|mate_chr| mate_chr != variant.chr)
                .unwrap_or(false);
        if !is_interchrom_paired_breakend {
            if let Some(chr_transcripts) = resources.transcripts.get(&variant.chr) {
                annotate_variant(
                    variant,
                    chr_transcripts,
                    &resources.effects_config,
                    &resources.prediction_config,
                );
            }
        }

        append_bnd_mate_consequences(
            variant,
            resources.transcripts.as_ref(),
            &resources.effects_config,
        );
        let csq_count = variant.transcript_consequences.len() as u64;

        if let Some(ref vd) = resources.variation_data {
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

        // Intergenic only when the chromosome has transcripts in the cache: Perl
        // VEP silently drops variants on chromosomes it never loaded. A BND never
        // gets the intergenic fallback: Perl emits consequences only via local
        // and mate annotation, never bare intergenic with Feature="-".
        let is_bnd = variant.variant_class == VariantClass::Translocation;
        if variant.transcript_consequences.is_empty()
            && variant.most_severe_consequence.is_none()
            && resources.transcripts.contains_key(&variant.chr)
            && !is_bnd
        {
            variant.most_severe_consequence = Some(Consequence::IntergenicVariant);
        }

        csq_count
    };

    if use_parallel {
        if stats_enabled {
            let (vars, csqs) = pool.install(|| {
                batch
                    .par_iter_mut()
                    .map(|v| (1u64, annotate_one(v)))
                    .reduce(|| (0u64, 0u64), |a, b| (a.0 + b.0, a.1 + b.1))
            });
            stats.variants_processed.fetch_add(vars, Ordering::Relaxed);
            stats
                .consequences_calculated
                .fetch_add(csqs, Ordering::Relaxed);
        } else {
            pool.install(|| {
                batch.par_iter_mut().for_each(|v| {
                    annotate_one(v);
                });
            });
        }
    } else if stats_enabled {
        let mut vars = 0u64;
        let mut csqs = 0u64;
        for v in batch.iter_mut() {
            csqs += annotate_one(v);
            vars += 1;
        }
        stats.variants_processed.fetch_add(vars, Ordering::Relaxed);
        stats
            .consequences_calculated
            .fetch_add(csqs, Ordering::Relaxed);
    } else {
        batch.iter_mut().for_each(|v| {
            annotate_one(v);
        });
    }
    let annotation_elapsed = annotation_phase_start.elapsed();
    debug!(
        batch_size = batch_len,
        annotation_ms = annotation_elapsed.as_millis() as u64,
        "batch annotation phase complete"
    );

    // Plugins run after consequence assignment so they see transcript
    // consequences and gene symbols; one batch call enables batched tabix queries.
    if !resources.builtin_plugins.is_empty() {
        let plugin_phase_start = Instant::now();
        if let Err(e) = resources.builtin_plugins.run_batch(batch) {
            tracing::warn!("plugin execution error: {}", e);
        }
        let plugin_elapsed = plugin_phase_start.elapsed();
        debug!(
            batch_size = batch_len,
            plugin_ms = plugin_elapsed.as_millis() as u64,
            "batch plugin phase complete"
        );
    }
    crate::dylib_plugins::run(&resources.dylib_plugins, batch)?;

    Ok(())
}

/// Configuration for SIFT/PolyPhen prediction lookup.
#[derive(Debug, Clone, Default)]
pub struct PredictionConfig {
    /// SIFT output mode: "p" (prediction), "s" (score), "b" (both), or None.
    pub sift: Option<String>,
    /// PolyPhen output mode: "p" (prediction), "s" (score), "b" (both), or None.
    pub polyphen: Option<String>,
}

/// Annotate a variant with transcript consequences and optional SIFT/PolyPhen predictions.
fn annotate_variant(
    variant: &mut InputVariant,
    transcripts: &TranscriptIndex,
    config: &vep_effects::EffectsConfig,
    prediction_config: &PredictionConfig,
) {
    annotate_variant_with_predictions(variant, transcripts, config, prediction_config);
}

fn annotate_variant_with_predictions(
    variant: &mut InputVariant,
    transcripts: &TranscriptIndex,
    config: &vep_effects::EffectsConfig,
    prediction_config: &PredictionConfig,
) {
    let (variant_start, variant_end) = structural_variant_bounds(variant);
    let has_predictions = prediction_config.sift.is_some() || prediction_config.polyphen.is_some();

    // For structural variants, use sv_end to query the full SV span so all
    // overlapping transcripts are found (not just those at the 1bp anchor).
    transcripts.for_each_overlapping(variant_start, variant_end, |transcript| {
        if let Some(mut tc) = vep_effects::calculate_consequences(variant, transcript, config) {
            if has_predictions && tc.consequences.contains(&Consequence::MissenseVariant) {
                lookup_predictions(&mut tc, transcript, prediction_config);
            }
            // LoFTEE's structural context is captured here because `transcript` is
            // out of scope in the plugin phase; gated so other runs pay nothing.
            if config.populate_loftee_context {
                tc.loftee_ctx = Some(Box::new(
                    vep_core::consequence::LofteeContext::from_transcript(transcript),
                ));
            }
            variant.transcript_consequences.push(tc);
        }
    });

    // Post-process <*> (unspecified ALT) variants: Perl VEP treats <*> as a
    // non-SV allele but assigns coding_sequence_variant instead of frameshift
    // because the allele length is unknown. Replace frameshift→coding_sequence_variant.
    if variant.alt_allele() == b"<*>" {
        for tc in &mut variant.transcript_consequences {
            if tc.consequences.contains(&Consequence::FrameshiftVariant) {
                tc.consequences
                    .retain(|c| *c != Consequence::FrameshiftVariant);
                tc.consequences.push(Consequence::CodingSequenceVariant);
                tc.consequences.sort_by_key(|c| c.rank());
                tc.impact = tc
                    .consequences
                    .iter()
                    .map(|c| c.impact())
                    .min_by_key(|i| match i {
                        Impact::HIGH => 0,
                        Impact::MODERATE => 1,
                        Impact::LOW => 2,
                        Impact::MODIFIER => 3,
                    })
                    .unwrap_or(Impact::MODIFIER);
            }
        }
    }

    if let Some(most_severe) = variant
        .transcript_consequences
        .iter()
        .flat_map(|tc| tc.consequences.iter())
        .copied()
        .min_by_key(|c| c.rank())
    {
        variant.most_severe_consequence = Some(most_severe);
    }
}

/// Look up SIFT/PolyPhen predictions for a missense consequence.
///
/// Parses the protein position and amino acids from the transcript consequence,
/// then queries the prediction matrix stored in the transcript's VEFC.
fn lookup_predictions(
    tc: &mut vep_core::consequence::TranscriptConsequence,
    transcript: &vep_core::transcript::Transcript,
    prediction_config: &PredictionConfig,
) {
    use vep_core::prediction::{format_prediction, AnalysisType};

    let position: usize = tc
        .protein_position
        .as_deref()
        .and_then(|s| s.split('-').next())
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    if position == 0 {
        return;
    }

    let alt_aa: u8 = tc
        .amino_acids
        .as_deref()
        .and_then(|s| s.split('/').nth(1))
        .and_then(|s| s.as_bytes().first().copied())
        .unwrap_or(0);
    if alt_aa == 0 || alt_aa == b'-' {
        return;
    }

    let predictions = match transcript
        .vefc
        .as_ref()
        .and_then(|v| v.protein_function_predictions.as_ref())
    {
        Some(p) => p,
        None => return,
    };

    if let Some(ref mode) = prediction_config.sift {
        if let Some(ref matrix) = predictions.sift {
            if let Some(pred) = matrix.lookup(position, alt_aa, AnalysisType::Sift) {
                tc.sift = Some(format_prediction(&pred, mode));
            }
        }
    }

    // PolyPhen (default: humvar, matching Perl VEP)
    if let Some(ref mode) = prediction_config.polyphen {
        if let Some(ref matrix) = predictions.polyphen_humvar {
            if let Some(pred) = matrix.lookup(position, alt_aa, AnalysisType::PolyPhenHumVar) {
                tc.polyphen = Some(format_prediction(&pred, mode));
            }
        }
    }
}

/// Append the mate-side consequences of a paired BND, in place.
///
/// The single implementation: `annotator.rs` calls this rather than carrying its
/// own, because two implementations of one Perl behaviour diverge (dedup against
/// locally annotated transcripts, `calculate_paired_mate` versus the generic
/// `annotate_variant`, `local_span_is_ranged`).
pub(crate) fn append_bnd_mate_consequences(
    variant: &mut InputVariant,
    transcripts_by_chr: &LazyTranscriptIndexes,
    config: &vep_effects::EffectsConfig,
) {
    // Perl VEP annotates BND alleles against transcripts at the mate position.
    // This produces:
    // - feature_truncation: when mate falls within a transcript body
    // - intergenic_variant: when mate is near but outside a transcript
    // Only applies to paired BND alleles (not single-breakend forms).
    if variant.variant_class != vep_core::variant::VariantClass::Translocation
        || variant.is_single_breakend
    {
        return;
    }

    let mate_chr = match variant.mate_chr.as_deref() {
        Some(chr) => chr,
        None => return,
    };
    let mate_pos = match variant.mate_pos {
        Some(pos) => pos,
        None => return,
    };

    // The mate's own chromosome having loaded transcripts is the only requirement,
    // and `get` below enforces it. No gate on the local chromosome: Perl loads a
    // region because a variant overlaps it, and the mate coordinate is such a
    // variant; nothing in `AnnotationSource::Cache` conditions the mate's region
    // on the local chromosome being cached.
    let chr_transcripts = match transcripts_by_chr.get(mate_chr) {
        Some(idx) => idx,
        None => return,
    };

    // Create a temporary variant at the mate position for annotation.
    // calculate_paired_mate uses variant.start as the breakpoint.
    let orig_start = variant.start;
    let orig_end = variant.end;
    variant.start = mate_pos;
    variant.end = mate_pos;

    let upstream = config.upstream_distance;
    let downstream = config.downstream_distance;

    // Whether the original local span covered more than one base, captured
    // before the overwrite above. Perl's mate-side context predicates read the
    // local variation feature while `feature_truncation` reads the mate, so a
    // point local span gets bare `feature_truncation`; see
    // `calculate_paired_mate`.
    let local_span_is_ranged = orig_end > orig_start;

    let mut mate_consequences = Vec::new();
    chr_transcripts.for_each_overlapping(mate_pos, mate_pos, |transcript| {
        if let Some(tc) = vep_effects::sv::breakend::calculate_paired_mate(
            variant,
            transcript,
            upstream,
            downstream,
            local_span_is_ranged,
        ) {
            mate_consequences.push(tc);
        }
    });

    variant.start = orig_start;
    variant.end = orig_end;

    let mut appended = 0usize;
    for tc in mate_consequences {
        let already_has = variant
            .transcript_consequences
            .iter()
            .any(|existing| existing.transcript_id == tc.transcript_id);
        if !already_has {
            variant.transcript_consequences.push(tc);
            appended += 1;
        }
    }

    // `annotate_variant` computes `most_severe_consequence` from the local
    // consequences before this runs, so a more severe mate-side term (or a
    // mate-only BND, where the field is None) needs a recompute; the downstream
    // `intergenic_variant` fallback fires only on an empty list, which this append
    // has just made non-empty. Conditional, so a call that appended nothing leaves
    // the field as `annotate_variant` set it.
    if appended > 0 {
        if let Some(most_severe) = variant
            .transcript_consequences
            .iter()
            .flat_map(|tc| tc.consequences.iter())
            .copied()
            .min_by_key(|c| c.rank())
        {
            variant.most_severe_consequence = Some(most_severe);
        }
    }
}

/// Annotate one VCF line through the CLI path and return per-variant
/// `transcript_id:so_terms` strings, sorted.
///
/// Exists so `annotator.rs`'s tests can compare the CLI path against the annotator path
/// without duplicating this file's `AnnotationResources` and thread-pool scaffolding.
/// Test-only, and deliberately outside `mod tests` because a test module is not visible
/// from another file.
#[cfg(test)]
pub(crate) fn annotate_batch_for_test(
    vcf_line: &str,
    transcripts: Arc<LazyTranscriptIndexes>,
) -> Vec<Vec<String>> {
    let resources = AnnotationResources {
        transcripts,
        variation_data: None,
        effects_config: Arc::new(vep_effects::EffectsConfig {
            upstream_distance: 5_000,
            downstream_distance: 5_000,
            ..Default::default()
        }),
        builtin_plugins: Arc::new(BuiltinRegistry::new()),
        dylib_plugins: Arc::new(PluginRegistry::new()),
        prediction_config: PredictionConfig::default(),
    };
    let stats = PipelineStats {
        variants_processed: AtomicU64::new(0),
        consequences_calculated: AtomicU64::new(0),
    };
    let pool = build_thread_pool(1).unwrap();
    let mut batch = crate::vcf_parser::parse_vcf_line(vcf_line, true).expect("parse test VCF line");
    annotate_batch(&mut batch, &resources, &pool, false, &stats, false).expect("annotate");
    batch
        .iter()
        .map(|v| {
            let mut ids: Vec<String> = v
                .transcript_consequences
                .iter()
                .map(|tc| {
                    format!(
                        "{}:{}",
                        tc.transcript_id,
                        tc.consequences
                            .iter()
                            .map(|c| c.so_term())
                            .collect::<Vec<_>>()
                            .join(",")
                    )
                })
                .collect();
            ids.sort();
            ids
        })
        .collect()
}

/// Two-exon test transcript on `chr` at `start` (1 kb span), shared by the runner
/// and annotator tests so both batch paths are pinned against the same fixture.
#[cfg(test)]
pub(crate) fn make_runner_test_transcript(
    stable_id: &str,
    gene_id: &str,
    chr: &str,
    start: u64,
) -> vep_core::transcript::Transcript {
    use vep_core::coordinate::Strand;
    use vep_core::transcript::{Exon, Intron, Transcript};
    let exons = vec![
        Exon {
            stable_id: None,
            start,
            end: start + 199,
            rank: 1,
            phase: -1,
            end_phase: -1,
        },
        Exon {
            stable_id: None,
            start: start + 800,
            end: start + 999,
            rank: 2,
            phase: -1,
            end_phase: -1,
        },
    ];
    let introns = vec![Intron {
        start: start + 200,
        end: start + 799,
        rank: 1,
    }];

    Transcript {
        stable_id: stable_id.into(),
        version: None,
        db_id: None,
        gene_stable_id: gene_id.into(),
        chr: chr.into(),
        start,
        end: start + 999,
        strand: Strand::Forward,
        biotype: "lncRNA".into(),
        source: "Ensembl".into(),
        description: None,
        gene_symbol: Some(stable_id.into()),
        gene_symbol_source: Some("TEST".into()),
        hgnc_id: None,
        gene_phenotype: None,
        canonical: true,
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
        exons,
        introns,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript_index::TranscriptIndexImpl;
    use crate::vcf_parser::parse_vcf_line;
    use std::collections::HashMap;
    use std::sync::atomic::AtomicU64;

    #[test]
    fn parquet_row_groups_are_whole_duckdb_vectors() {
        assert_eq!(parquet_row_group_rows(122_880), 122_880);
        assert_eq!(parquet_row_group_rows(5_000), 6_144);
        assert_eq!(parquet_row_group_rows(1), 2_048);
        assert_eq!(parquet_row_group_rows(0), 2_048);
    }

    #[test]
    fn command_line_prints_like_vep() {
        let args = [
            "-i",
            "/input/sample.vcf",
            "--offline",
            "--json_cache=/caches/grch37",
            "-o",
            "out.txt",
            "--force",
            "--species",
            "homo_sapiens",
        ];
        assert_eq!(
            vep_command_line(args.iter().map(|s| s.to_string())),
            "vep --force_overwrite --input_file [PATH]/sample.vcf --json_cache [PATH]/grch37 \
             --offline --output_file out.txt --species homo_sapiens"
        );
    }
    use std::sync::Arc;
    use vep_core::variant::VariantClass;

    /// Build a `Config` from real CLI flags so the test tracks flag resolution
    /// and never breaks when unrelated `Config` fields are added.
    fn config_from(flags: &[&str]) -> Config {
        use clap::Parser;
        let mut argv = vec!["vep"];
        argv.extend_from_slice(flags);
        Config::from_args(crate::args::Args::parse_from(argv))
    }

    /// The capture gate must be derived from the output format, not guessed.
    /// `vcf` is the only format whose writer reads `raw_input`; because that
    /// writer silently falls back to a reconstructed record when the field is
    /// absent, a wrong gate here degrades output quietly instead of erroring.
    #[test]
    fn raw_input_is_captured_for_the_formats_that_re_emit_the_line() {
        for fmt in ["vcf", "json"] {
            assert!(
                needs_raw_input_capture(&config_from(&["--output_format", fmt])),
                "{fmt} re-emits the input line"
            );
        }
        assert!(needs_raw_input_capture(&config_from(&["--vcf"])));
        assert!(needs_raw_input_capture(&config_from(&["--json"])));
        assert!(!needs_raw_input_capture(&config_from(&[])));
        for fmt in ["vep", "tab", "parquet"] {
            assert!(
                !needs_raw_input_capture(&config_from(&["--output_format", fmt])),
                "{fmt} does not read raw_input"
            );
        }
    }

    #[test]
    fn hgvs_is_computed_only_under_the_flag() {
        for fmt in ["vep", "tab", "vcf", "json", "parquet"] {
            assert!(!needs_hgvs_computation(&config_from(&[
                "--output_format",
                fmt
            ])));
            assert!(needs_hgvs_computation(&config_from(&[
                "--hgvs",
                "--output_format",
                fmt
            ])));
        }
        assert!(!needs_hgvs_computation(&config_from(&["--vcf"])));
        assert!(needs_hgvs_computation(&config_from(&["--hgvs"])));
    }

    /// The exon/intron ordinals reach the output only under `--numbers`; `--vcf`
    /// and `--everything` imply it, and a loaded plugin may read them off the
    /// consequence, so those runs compute them too.
    #[test]
    fn exon_intron_numbers_are_computed_only_when_something_reads_them() {
        for fmt in ["vep", "tab", "json", "parquet"] {
            assert!(!needs_exon_intron_numbers(&config_from(&[
                "--output_format",
                fmt
            ])));
            assert!(needs_exon_intron_numbers(&config_from(&[
                "--numbers",
                "--output_format",
                fmt
            ])));
        }
        assert!(!needs_exon_intron_numbers(&config_from(&["--hgvs"])));
        assert!(needs_exon_intron_numbers(&config_from(&["--vcf"])));
        assert!(needs_exon_intron_numbers(&config_from(&[
            "--output_format",
            "vcf"
        ])));
        assert!(needs_exon_intron_numbers(&config_from(&["--everything"])));
        assert!(needs_exon_intron_numbers(&config_from(&[
            "--plugin",
            "LoFTEE,loftee_path:/x"
        ])));
    }

    /// Write a JSON cache holding one 1,000 bp transcript per `(chr, start)`, in
    /// the array format `vep-cache-builder` emits. Self-cleaning on drop.
    fn write_json_cache(tag: &str, entries: &[(&str, u64)]) -> tempfile::TempDir {
        let dir = tempfile::Builder::new()
            .prefix(&format!("vep_runner_{tag}_"))
            .tempdir()
            .unwrap();
        for (chr, start) in entries {
            let chr_dir = dir.path().join("transcripts").join(chr);
            std::fs::create_dir_all(&chr_dir).unwrap();
            std::fs::write(
                chr_dir.join("0-1000.json"),
                format!(
                    r#"[{{"stable_id":"ENST{chr}","gene_stable_id":"ENSG{chr}","chr":"{chr}","start":{start},"end":{end},"strand":1,"biotype":"lncRNA","source":"Ensembl"}}]"#,
                    end = start + 999
                ),
            )
            .unwrap();
        }
        dir
    }

    fn resources_over(transcripts: Arc<LazyTranscriptIndexes>) -> AnnotationResources {
        AnnotationResources {
            transcripts,
            variation_data: None,
            effects_config: Arc::new(vep_effects::EffectsConfig {
                upstream_distance: 5_000,
                downstream_distance: 5_000,
                ..Default::default()
            }),
            builtin_plugins: Arc::new(BuiltinRegistry::new()),
            dylib_plugins: Arc::new(PluginRegistry::new()),
            prediction_config: PredictionConfig::default(),
        }
    }

    fn zero_stats() -> PipelineStats {
        PipelineStats {
            variants_processed: AtomicU64::new(0),
            consequences_calculated: AtomicU64::new(0),
        }
    }

    /// A malformed shard must abort `annotate_batch` with an error, even under
    /// parallel annotation.
    ///
    /// This is the pre-warm loop's own test. Every other
    /// `annotate_batch` test builds its index with `from_indexes`, which pre-fills
    /// every slot, so `get` never misses and deleting the loop entirely leaves them
    /// all green. Here the index is built lazily from a cache with a truncated
    /// shard, so the loop's `prewarm(...)?` is the only thing that can surface the
    /// failure; without it, every variant on chr21 would be annotated
    /// `intergenic_variant` and the run would exit 0, because `contains_key`
    /// answers from the directory listing.
    ///
    /// `use_parallel = true` with a 4-thread pool also exercises the deadlock
    /// ordering: pre-warming happens before `pool.install`, so nothing parses from
    /// inside a rayon worker.
    #[test]
    fn test_annotate_batch_propagates_malformed_shard_error() {
        let dir = write_json_cache("malformed", &[("21", 1_000)]);
        std::fs::write(
            dir.path()
                .join("transcripts")
                .join("21")
                .join("0-1000.json"),
            "not json",
        )
        .unwrap();

        let transcripts = Arc::new(
            LazyTranscriptIndexes::new(
                dir.path().to_str().unwrap(),
                TranscriptIndexImpl::Bin,
                5_000,
                5_000,
            )
            .unwrap(),
        );
        assert!(
            transcripts.contains_key("21"),
            "the chromosome reads as present; that is the trap the loop closes"
        );

        let resources = resources_over(transcripts);
        let stats = zero_stats();
        let pool = build_thread_pool(4).unwrap();
        let mut batch =
            parse_vcf_line("21\t1300\trs1\tA\tG\t.\t.\t.", true).expect("parse test variant");

        let err = annotate_batch(&mut batch, &resources, &pool, true, &stats, false)
            .expect_err("a truncated shard must abort the batch, not annotate as intergenic");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("chr21"),
            "error must name the chromosome, got: {msg}"
        );
    }

    /// Inter-chromosomal BND annotation must be identical lazy vs eager.
    ///
    /// The mate chromosome is reached through `prewarm(mate_chr)` in the pre-warm
    /// loop and `get(mate_chr)` in `append_bnd_mate_consequences`, two different
    /// code paths that must agree on the same chromosome name. If they ever
    /// disagreed (e.g. one normalized `chr21` and the other did not), `get` would
    /// return `None` and the mate breakpoint would be silently unannotated while
    /// the run still exited 0.
    ///
    /// Both cache shapes are covered, and both must annotate the mate. Whether the
    /// local chromosome is cached has no bearing on whether the mate coordinate falls
    /// inside a transcript, so a chr21-only cache and a chr1+chr21 cache must agree.
    /// Naming the transcript in both shapes is what makes this meaningful: the
    /// lazy-vs-eager equality would hold for two identically-empty results.
    #[test]
    fn test_annotate_batch_interchrom_bnd_lazy_matches_eager() {
        let bnd_line = "1\t1300\tbnd_1\tA\tA[21:25001300[\t.\t.\tSVTYPE=BND;MATEID=bnd_2";

        let annotate_with = |transcripts: Arc<LazyTranscriptIndexes>| -> Vec<Vec<String>> {
            let resources = resources_over(transcripts);
            let stats = zero_stats();
            let pool = build_thread_pool(1).unwrap();
            let mut batch = parse_vcf_line(bnd_line, true).expect("parse BND");
            annotate_batch(&mut batch, &resources, &pool, false, &stats, false).expect("annotate");
            batch
                .iter()
                .map(|v| {
                    let mut ids: Vec<String> = v
                        .transcript_consequences
                        .iter()
                        .map(|tc| tc.transcript_id.to_string())
                        .collect();
                    ids.sort();
                    ids
                })
                .collect()
        };

        let compare = |tag: &str, chrs: &[(&str, u64)]| -> Vec<Vec<String>> {
            let dir = write_json_cache(tag, chrs);
            let path = dir.path().to_str().unwrap();
            let lazy = annotate_with(Arc::new(
                LazyTranscriptIndexes::new(path, TranscriptIndexImpl::Bin, 5_000, 5_000).unwrap(),
            ));
            let eager_raw = crate::json_cache::load_all_transcripts(path).unwrap();
            let eager = annotate_with(Arc::new(LazyTranscriptIndexes::from_loaded(
                eager_raw,
                TranscriptIndexImpl::Bin,
                5_000,
                5_000,
            )));
            assert_eq!(
                lazy, eager,
                "{tag}: inter-chromosomal BND consequences must be identical lazy vs eager"
            );
            lazy
        };

        // Full-genome shape: chr1 present, so the chr21 mate is annotated. This is
        // the case that exercises `prewarm(mate_chr)` + `get(mate_chr)` agreeing on
        // the chromosome name.
        let both = compare("bnd_both", &[("1", 1_000), ("21", 25_000_000)]);
        assert!(
            both.iter().any(|ids| ids.iter().any(|id| id == "ENST21")),
            "the chr21 mate transcript must be annotated when chr1 is cached, got: {both:?}"
        );

        // Single-chromosome shape: chr1 absent. The mate is still annotated, because
        // whether the local chromosome is cached has no bearing on whether the mate
        // coordinate falls inside a transcript.
        //
        // Asserting both cache shapes is the point. The lazy-vs-eager equality inside
        // `compare` would hold for two identically-empty results, so without naming the
        // transcript the check would pass on a mate that is never annotated at all.
        let mate_only = compare("bnd_mate_only", &[("21", 25_000_000)]);
        assert!(
            mate_only
                .iter()
                .any(|ids| ids.iter().any(|id| id == "ENST21")),
            "the chr21 mate transcript must be annotated even when the local chr1 is \
             uncached, got: {mate_only:?}"
        );
    }

    /// Annotation must be identical whether `--fork` fans out or runs serially.
    ///
    /// `scripts/concordance/run_clone_measurement.sh` passes the same `--fork` to
    /// every engine on every cell, and an SV cell's wall time and its concordance
    /// F1 come from that one invocation, so a `--fork` that could move the output
    /// would make each of those a different measurement than the accompanying
    /// manuscript reports.
    ///
    /// Order is asserted, not just set equality. `par_iter_mut` annotates in
    /// whatever order workers pick up variants, and output is emitted by walking
    /// the batch afterwards, so a change that collected results instead of
    /// mutating in place would reorder the file while leaving every tuple present
    /// which is invisible to a set comparison and to any F1, which deduplicates.
    ///
    /// The batch spans two chromosomes and four input lines deliberately: a
    /// single-variant batch cannot detect a reordering, and a single chromosome
    /// would not exercise the pre-warm loop's multi-chromosome path.
    ///
    /// Breakends rather than SNVs, for two reasons. They are what this fixture
    /// transcript annotates (`make_runner_test_transcript` is a CDS-less lncRNA
    /// and the SNV path yields nothing on it, which the emptiness guard below
    /// would catch). And they do strictly more per-variant work inside the
    /// fan-out: `annotate_one` takes both sides of its
    /// `is_interchrom_paired_breakend` split across these four lines, and
    /// `append_bnd_mate_consequences` reaches into the other chromosome's index
    /// while a worker is annotating this one: the cross-variant index access
    /// where a parallel-vs-serial divergence would most plausibly hide.
    #[test]
    fn test_annotate_batch_fork_1_matches_fork_4() {
        // Transcript 1000-1999 on each chromosome: exon 1 1000-1199, intron
        // 1200-1799, exon 2 1800-1999. Lines 1 and 2 are same-chromosome paired
        // breakends (local annotation kept); lines 3 and 4 are inter-chromosomal
        // in both directions (local suppressed, mate annotated). Each line parses
        // to a paired plus a single-breakend variant, so the batch is 8 variants.
        let lines = [
            "1\t1100\tbnd_a\tA\tA[1:1800[\t.\t.\tSVTYPE=BND;MATEID=bnd_a2",
            "21\t1300\tbnd_b\tA\tA[21:1850[\t.\t.\tSVTYPE=BND;MATEID=bnd_b2",
            "1\t1150\tbnd_c\tA\tA[21:1400[\t.\t.\tSVTYPE=BND;MATEID=bnd_c2",
            "21\t1900\tbnd_d\tA\tA[1:1200[\t.\t.\tSVTYPE=BND;MATEID=bnd_d2",
        ];

        // (chr, start, allele, transcript, consequence terms) per consequence, in
        // exactly the order the writer would walk.
        let annotate_at = |threads: usize| -> Vec<(String, u64, String, String, String)> {
            let mut by_chr = HashMap::new();
            for chr in ["1", "21"] {
                let tx = make_runner_test_transcript(
                    &format!("ENST{chr}"),
                    &format!("ENSG{chr}"),
                    chr,
                    1_000,
                );
                by_chr.insert(
                    chr.to_string(),
                    TranscriptIndex::new(TranscriptIndexImpl::Bin, vec![tx], 5_000, 5_000),
                );
            }
            let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(by_chr));
            let resources = resources_over(transcripts);
            let stats = zero_stats();
            let pool = build_thread_pool(threads).unwrap();
            let mut batch = Vec::new();
            for line in lines {
                batch.extend(parse_vcf_line(line, true).expect("parse test variant"));
            }
            annotate_batch(&mut batch, &resources, &pool, threads > 1, &stats, false)
                .expect("annotate");
            batch
                .iter()
                .flat_map(|v| {
                    v.transcript_consequences.iter().map(move |tc| {
                        (
                            v.chr.to_string(),
                            v.start,
                            v.allele_string.clone(),
                            tc.transcript_id.to_string(),
                            tc.consequences
                                .iter()
                                .map(|c| c.so_term())
                                .collect::<Vec<_>>()
                                .join(","),
                        )
                    })
                })
                .collect()
        };

        let serial = annotate_at(1);
        let parallel = annotate_at(4);

        // Guard against the two-empty-results false pass: a fork-invariance check
        // over a batch nothing annotated holds trivially, and this fixture's SNV
        // path does produce exactly that.
        assert!(
            serial.len() >= lines.len(),
            "expected at least one consequence per input line, got {}: {serial:?}",
            serial.len()
        );
        assert!(
            serial.iter().any(|c| c.3 == "ENST1") && serial.iter().any(|c| c.3 == "ENST21"),
            "both chromosomes' transcripts must be annotated, got: {serial:?}"
        );

        assert_eq!(
            serial, parallel,
            "annotation differs between --fork 1 and --fork 4; --fork is not a pure \
             performance knob and no cross-fork measurement comparison is sound"
        );
    }

    #[test]
    fn test_append_bnd_mate_consequences_interchrom_annotates_mate() {
        // Perl VEP annotates inter-chromosomal BND mate-side transcripts when the
        // mate chromosome's cache is loaded; here every loaded transcript is
        // available, so the annotation is unconditional.
        let local_tx = make_runner_test_transcript("ENSTLOCAL", "ENSGLOCAL", "1", 1_000);
        let mate_tx = make_runner_test_transcript("ENSTMATE", "ENSGMATE", "21", 25_000_000);

        let mut by_chr = HashMap::new();
        by_chr.insert(
            "1".to_string(),
            TranscriptIndex::new(TranscriptIndexImpl::Bin, vec![local_tx], 5_000, 5_000),
        );
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(TranscriptIndexImpl::Bin, vec![mate_tx], 5_000, 5_000),
        );

        let config = vep_effects::EffectsConfig::default();
        let mut variant = InputVariant::new(
            "1".into(),
            1_300,
            1_300,
            b"A".to_vec(),
            b"A[21:25001300[".to_vec(),
        );
        variant.variant_class = VariantClass::Translocation;
        variant.is_structural = true;
        variant.sv_end = Some(1_300);
        variant.mate_chr = Some("21".into());
        variant.mate_pos = Some(25_001_300);

        annotate_variant(
            &mut variant,
            by_chr.get("1").unwrap(),
            &config,
            &PredictionConfig::default(),
        );
        assert_eq!(variant.transcript_consequences.len(), 1);

        append_bnd_mate_consequences(
            &mut variant,
            &LazyTranscriptIndexes::from_indexes(by_chr),
            &config,
        );

        assert_eq!(variant.transcript_consequences.len(), 2);
        let ids: Vec<&str> = variant
            .transcript_consequences
            .iter()
            .map(|tc| &*tc.transcript_id)
            .collect();
        assert!(ids.contains(&"ENSTLOCAL"));
        assert!(ids.contains(&"ENSTMATE"));
    }

    #[test]
    fn test_append_bnd_mate_consequences_skips_single_breakend() {
        let local_tx = make_runner_test_transcript("ENSTLOCAL", "ENSGLOCAL", "1", 1_000);
        let mate_tx = make_runner_test_transcript("ENSTMATE", "ENSGMATE", "21", 25_000_000);

        let mut by_chr = HashMap::new();
        by_chr.insert(
            "1".to_string(),
            TranscriptIndex::new(TranscriptIndexImpl::Bin, vec![local_tx], 5_000, 5_000),
        );
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(TranscriptIndexImpl::Bin, vec![mate_tx], 5_000, 5_000),
        );

        let config = vep_effects::EffectsConfig::default();
        let mut variant =
            InputVariant::new("1".into(), 1_300, 1_300, b"A".to_vec(), b"A.".to_vec());
        variant.variant_class = VariantClass::Translocation;
        variant.is_structural = true;
        variant.is_single_breakend = true;
        variant.sv_end = Some(1_300);
        variant.mate_chr = Some("21".into());
        variant.mate_pos = Some(25_001_300);

        annotate_variant(
            &mut variant,
            by_chr.get("1").unwrap(),
            &config,
            &PredictionConfig::default(),
        );
        append_bnd_mate_consequences(
            &mut variant,
            &LazyTranscriptIndexes::from_indexes(by_chr),
            &config,
        );

        assert_eq!(variant.transcript_consequences.len(), 1);
        assert_eq!(
            &*variant.transcript_consequences[0].transcript_id,
            "ENSTLOCAL"
        );
    }

    #[test]
    fn test_append_bnd_mate_consequences_interchrom_no_transcripts_noop() {
        let local_tx = make_runner_test_transcript("ENSTLOCAL", "ENSGLOCAL", "1", 1_000);

        let mut by_chr = HashMap::new();
        by_chr.insert(
            "1".to_string(),
            TranscriptIndex::new(TranscriptIndexImpl::Bin, vec![local_tx], 5_000, 5_000),
        );

        let config = vep_effects::EffectsConfig::default();
        let mut variant = InputVariant::new(
            "1".into(),
            1_300,
            1_300,
            b"A".to_vec(),
            b"A[21:25001300[".to_vec(),
        );
        variant.variant_class = VariantClass::Translocation;
        variant.is_structural = true;
        variant.sv_end = Some(1_300);
        variant.mate_chr = Some("21".into());
        variant.mate_pos = Some(25_001_300);

        annotate_variant(
            &mut variant,
            by_chr.get("1").unwrap(),
            &config,
            &PredictionConfig::default(),
        );
        assert_eq!(variant.transcript_consequences.len(), 1);

        append_bnd_mate_consequences(
            &mut variant,
            &LazyTranscriptIndexes::from_indexes(by_chr),
            &config,
        );

        assert_eq!(variant.transcript_consequences.len(), 1);
        assert_eq!(
            &*variant.transcript_consequences[0].transcript_id,
            "ENSTLOCAL"
        );
    }

    #[test]
    fn test_annotate_batch_interchrom_bnd_mate_annotation() {
        // Inter-chromosomal paired BND: local annotation suppressed (Perl parity),
        // but mate-side annotation enabled via append_bnd_mate_consequences.
        // Single-breakend form gets local annotation but not mate annotation.
        let local_tx = make_runner_test_transcript("ENSTLOCAL", "ENSGLOCAL", "1", 1_000);
        let mate_tx = make_runner_test_transcript("ENSTMATE", "ENSGMATE", "21", 25_000_000);

        let mut by_chr = HashMap::new();
        by_chr.insert(
            "1".to_string(),
            TranscriptIndex::new(TranscriptIndexImpl::Bin, vec![local_tx], 5_000, 5_000),
        );
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(TranscriptIndexImpl::Bin, vec![mate_tx], 5_000, 5_000),
        );

        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(by_chr));
        let variation_data: Option<Arc<json_cache::VariationData>> = None;
        let effects_config = Arc::new(vep_effects::EffectsConfig::default());
        let builtin_plugins = Arc::new(BuiltinRegistry::new());
        let resources = AnnotationResources {
            transcripts,
            variation_data,
            effects_config,
            builtin_plugins,
            dylib_plugins: Arc::new(PluginRegistry::new()),
            prediction_config: PredictionConfig::default(),
        };
        let stats = PipelineStats {
            variants_processed: AtomicU64::new(0),
            consequences_calculated: AtomicU64::new(0),
        };

        let mut batch = parse_vcf_line(
            "1\t1300\tbnd_1\tA\tA[21:25001300[\t.\t.\tSVTYPE=BND;MATEID=bnd_2",
            true,
        )
        .unwrap();
        let pool = build_thread_pool(1).unwrap();
        annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();

        let paired = batch.iter().find(|v| !v.is_single_breakend).unwrap();
        let single = batch.iter().find(|v| v.is_single_breakend).unwrap();

        let paired_ids: Vec<&str> = paired
            .transcript_consequences
            .iter()
            .map(|tc| &*tc.transcript_id)
            .collect();
        let single_ids: Vec<&str> = single
            .transcript_consequences
            .iter()
            .map(|tc| &*tc.transcript_id)
            .collect();

        assert!(!paired_ids.contains(&"ENSTLOCAL"));
        assert!(paired_ids.contains(&"ENSTMATE"));
        assert!(single_ids.contains(&"ENSTLOCAL"));
        assert!(!single_ids.contains(&"ENSTMATE"));
    }

    #[test]
    fn test_annotate_batch_same_chrom_paired_breakend_keeps_local_annotation() {
        let local_tx = make_runner_test_transcript("ENSTLOCAL", "ENSGLOCAL", "21", 1_000);

        let mut by_chr = HashMap::new();
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(TranscriptIndexImpl::Bin, vec![local_tx], 5_000, 5_000),
        );

        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(by_chr));
        let variation_data: Option<Arc<json_cache::VariationData>> = None;
        let effects_config = Arc::new(vep_effects::EffectsConfig::default());
        let builtin_plugins = Arc::new(BuiltinRegistry::new());
        let resources = AnnotationResources {
            transcripts,
            variation_data,
            effects_config,
            builtin_plugins,
            dylib_plugins: Arc::new(PluginRegistry::new()),
            prediction_config: PredictionConfig::default(),
        };
        let stats = PipelineStats {
            variants_processed: AtomicU64::new(0),
            consequences_calculated: AtomicU64::new(0),
        };

        let mut batch = parse_vcf_line(
            "21\t1300\tbnd_1\tA\tA[21:1800[\t.\t.\tSVTYPE=BND;MATEID=bnd_2",
            true,
        )
        .unwrap();
        let pool = build_thread_pool(1).unwrap();
        annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();

        let paired = batch.iter().find(|v| !v.is_single_breakend).unwrap();
        let single = batch.iter().find(|v| v.is_single_breakend).unwrap();

        assert_eq!(paired.transcript_consequences.len(), 1);
        assert_eq!(single.transcript_consequences.len(), 1);
        assert_eq!(
            &*paired.transcript_consequences[0].transcript_id,
            "ENSTLOCAL"
        );
        assert_eq!(
            &*single.transcript_consequences[0].transcript_id,
            "ENSTLOCAL"
        );
    }

    /// A giant breakend (span above VEP's 10 Mb `--max_sv_size` default) is annotated
    /// against the chromosome's complete transcript set, and its output does not depend
    /// on which other records share its batch. VEP marks such breakends skipped at
    /// cache-region load and annotates them only against transcripts their batch-mates
    /// loaded. This pins batch independence: the same tuples alone as with a neighbour,
    /// and the breakend's own distant transcripts are among them.
    #[test]
    fn test_annotate_batch_large_bnd_is_batch_independent() {
        fn run(with_neighbor: bool) -> Vec<(bool, Vec<String>)> {
            let buffered_tx =
                make_runner_test_transcript("ENSTBUFFERED", "ENSGBUFFERED", "21", 1_000);
            let distant_tx =
                make_runner_test_transcript("ENSTDISTANT", "ENSGDISTANT", "21", 20_000_000);
            let mut by_chr = HashMap::new();
            by_chr.insert(
                "21".to_string(),
                TranscriptIndex::new(
                    TranscriptIndexImpl::Bin,
                    vec![buffered_tx, distant_tx],
                    5_000,
                    5_000,
                ),
            );
            let resources = AnnotationResources {
                transcripts: Arc::new(LazyTranscriptIndexes::from_indexes(by_chr)),
                variation_data: None,
                effects_config: Arc::new(vep_effects::EffectsConfig::default()),
                builtin_plugins: Arc::new(BuiltinRegistry::new()),
                dylib_plugins: Arc::new(PluginRegistry::new()),
                prediction_config: PredictionConfig::default(),
            };
            let stats = PipelineStats {
                variants_processed: AtomicU64::new(0),
                consequences_calculated: AtomicU64::new(0),
            };
            let mut batch = Vec::new();
            if with_neighbor {
                batch.extend(parse_vcf_line("21\t1500\tnearby\tA\tG\t.\t.\t.", true).unwrap());
            }
            batch.extend(
                parse_vcf_line(
                    "21\t999\tbig_bnd\tN\t<BND>\t.\t.\tSVTYPE=BND;END=1000;SVLEN=20000000;CHR2=21;END2=20001000",
                    true,
                )
                .unwrap(),
            );
            let pool = build_thread_pool(1).unwrap();
            annotate_batch(&mut batch, &resources, &pool, false, &stats, false).unwrap();
            batch
                .iter()
                .filter(|v| v.id.as_deref() == Some("big_bnd"))
                .map(|v| {
                    let mut ids: Vec<String> = v
                        .transcript_consequences
                        .iter()
                        .map(|tc| tc.transcript_id.to_string())
                        .collect();
                    ids.sort();
                    ids.dedup();
                    (v.is_single_breakend, ids)
                })
                .collect()
        }

        let alone = run(false);
        let with_neighbor = run(true);
        assert_eq!(
            alone, with_neighbor,
            "a giant breakend's annotation must not depend on its batch-mates"
        );
        let single = alone
            .iter()
            .find(|(is_single, _)| *is_single)
            .expect("the BND record yields a single-breakend entry");
        assert!(
            single.1.iter().any(|id| id == "ENSTDISTANT"),
            "the breakend is annotated against its own span's transcripts, not only its neighbours': {:?}",
            single.1
        );
    }

    #[test]
    fn test_annotate_imprecise_sv_uses_nominal_coordinates() {
        // Perl VEP offline mode does not use outer_start/outer_end (CIPOS/CIEND
        // expansion) for transcript overlap. Only nominal coordinates are used.
        let tx_inside = make_runner_test_transcript("ENST_INSIDE", "GENE_INSIDE", "21", 27255000);
        // tx_edge spans 27240000-27240999, well outside nominal DEL start (27254100)
        // plus upstream distance (5000bp). Even with CIPOS=-14500 it should not
        // be annotated because Perl VEP offline mode uses nominal coordinates.
        let tx_edge = make_runner_test_transcript("ENST_EDGE", "GENE_EDGE", "21", 27240000);
        let transcripts = TranscriptIndex::new(
            TranscriptIndexImpl::Bin,
            vec![tx_inside, tx_edge],
            5_000,
            5_000,
        );
        let config = vep_effects::EffectsConfig::default();

        let mut variant = InputVariant::new(
            "21".to_string(),
            27254100,
            27260000,
            b"N".to_vec(),
            b"<DEL>".to_vec(),
        );
        variant.variant_class = VariantClass::StructuralDeletion;
        variant.is_structural = true;
        variant.sv_end = Some(27260000);
        variant.ci_pos = Some((-14500, 0));

        annotate_variant(
            &mut variant,
            &transcripts,
            &config,
            &PredictionConfig::default(),
        );
        let features: Vec<&str> = variant
            .transcript_consequences
            .iter()
            .map(|tc| tc.transcript_id.as_ref())
            .collect();
        assert!(
            features.contains(&"ENST_INSIDE"),
            "tx within nominal range must be annotated"
        );
        assert!(
            !features.contains(&"ENST_EDGE"),
            "tx outside nominal range should not be annotated even with CIPOS"
        );
    }
}
