// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reusable annotation engine for file-based and service-based entrypoints.

use std::sync::Arc;
use std::time::Instant;

use anyhow::Context;
use rayon::prelude::*;
use rayon::ThreadPool;
use thiserror::Error;
use tracing::{debug, info};

use crate::config::Config;
use crate::json_cache;
use crate::transcript_index::{
    LazyTranscriptIndexes, MaterializedTranscriptIndexes, TranscriptIndex, TranscriptIndexImpl,
};
use crate::variation_matcher;
use crate::vcf_parser::parse_vcf_line;
use vep_builtins::BuiltinRegistry;
use vep_core::consequence::{Consequence, Impact};
use vep_core::variant::{InputVariant, VariantClass};
use vep_plugin::PluginRegistry;

/// Shared annotation resources that can be reused across batches.
pub struct AnnotationResources {
    pub transcripts: Arc<LazyTranscriptIndexes>,
    pub variation_data: Option<Arc<json_cache::VariationData>>,
    pub effects_config: Arc<vep_effects::EffectsConfig>,
    pub builtin_plugins: Arc<BuiltinRegistry>,
    /// C-ABI plugins loaded from the `--plugin` arguments no built-in claimed.
    pub dylib_plugins: Arc<PluginRegistry>,
    pub prediction_config: PredictionConfig,
}

/// Build a scoped rayon thread pool sized for one annotator.
///
/// `rayon::ThreadPoolBuilder::build_global()` is avoided: concurrent embedding
/// consumers (long-lived services, request handlers) would have every in-flight
/// request fighting for one shared global pool of N threads, oversubscribing the
/// CPU. With a scoped pool each annotator's parallel work runs inside
/// `pool.install(...)` and stays bounded to its configured fork count.
pub fn build_thread_pool(num_threads: usize) -> anyhow::Result<Arc<ThreadPool>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .thread_name(|i| format!("vep-annotator-{i}"))
        .build()
        .context("failed to build rayon thread pool for annotator")?;
    Ok(Arc::new(pool))
}

/// Configuration for SIFT/PolyPhen prediction lookup.
#[derive(Debug, Clone, Default)]
pub struct PredictionConfig {
    /// SIFT output mode: "p" (prediction), "s" (score), "b" (both), or None.
    pub sift: Option<String>,
    /// PolyPhen output mode: "p" (prediction), "s" (score), "b" (both), or None.
    pub polyphen: Option<String>,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct AnnotationStats {
    pub variants_processed: u64,
    pub consequences_calculated: u64,
}

/// Per-row error returned by [`Annotator::try_annotate_rows`].
///
/// Embedding consumers want to map a single bad row to a row-level error
/// rather than failing the whole batch. The
/// taxonomy distinguishes input/configuration errors (caller's fault, return
/// to caller as a row error) from engine errors (annotator's fault, but still
/// row-scoped, never panicking the whole batch).
#[derive(Debug, Clone, Error)]
pub enum AnnotationError {
    /// Row failed pre-annotation validation (empty allele, ALT='.', etc.).
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// Assembly string did not match a loaded annotator.
    #[error("unsupported assembly: {0}")]
    UnsupportedAssembly(String),
    /// Engine raised an error mid-annotation. Unused: the annotation hot path
    /// is infallible, so the variant exists for API stability.
    #[error("annotation engine error: {0}")]
    Annotation(String),
}

/// A simple-allele row: chr, 1-based pos, ref, alt. The shape a batch service
/// consumer produces after splitting a `[row_id, asm, chr, pos, ref, alt]`
/// payload.
///
/// This is the per-row library input. For richer VCF features (multi-allelic,
/// symbolic SVs, breakends) callers should construct `InputVariant` themselves
/// and use [`Annotator::annotate_batch`] / [`Annotator::try_annotate_batch`]
/// directly.
#[derive(Debug, Clone)]
pub struct RawRow {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
}

/// Validate a [`RawRow`] against the simple-allele constraints embedding
/// consumers rely on.
pub fn validate_raw_row(row: &RawRow) -> Result<(), AnnotationError> {
    if row.pos == 0 {
        return Err(AnnotationError::InvalidInput(
            "position must be 1-based and > 0".into(),
        ));
    }
    if row.chrom.trim().is_empty() {
        return Err(AnnotationError::InvalidInput(
            "chromosome must not be empty".into(),
        ));
    }
    if row.ref_allele.trim().is_empty() || row.alt_allele.trim().is_empty() {
        return Err(AnnotationError::InvalidInput(
            "ref and alt alleles must not be empty".into(),
        ));
    }
    if row.alt_allele == "." {
        return Err(AnnotationError::InvalidInput(
            "ALT='.' is not supported by the simple-allele API".into(),
        ));
    }
    if row.ref_allele.contains(',') || row.alt_allele.contains(',') {
        return Err(AnnotationError::InvalidInput(
            "multi-allelic rows are not supported; split alleles before calling".into(),
        ));
    }
    if row.alt_allele.contains('<') || row.alt_allele.contains('[') || row.alt_allele.contains(']')
    {
        return Err(AnnotationError::InvalidInput(
            "symbolic and breakend alleles are not supported by the simple-allele API".into(),
        ));
    }
    Ok(())
}

/// Convert a [`RawRow`] into a single [`InputVariant`] using the same VCF
/// parser the file pipeline uses, after running [`validate_raw_row`].
pub fn row_to_input_variant(row: &RawRow) -> Result<InputVariant, AnnotationError> {
    validate_raw_row(row)?;
    let line = format!(
        "{}\t{}\t.\t{}\t{}\t.\tPASS\t.",
        row.chrom, row.pos, row.ref_allele, row.alt_allele
    );
    // `raw_input` capture is unconditional: embedders may read it, the
    // synthesized line is the only record of the row, and this path annotates
    // one row at a time, so the per-variant cost is irrelevant.
    let mut variants = parse_vcf_line(&line, true).map_err(|err| {
        AnnotationError::InvalidInput(format!(
            "VCF parse error for {}:{} {}>{}: {err}",
            row.chrom, row.pos, row.ref_allele, row.alt_allele
        ))
    })?;
    if variants.len() != 1 {
        return Err(AnnotationError::InvalidInput(format!(
            "{}:{} {}>{} produced {} alleles; multi-allelic inputs are not supported",
            row.chrom,
            row.pos,
            row.ref_allele,
            row.alt_allele,
            variants.len()
        )));
    }
    Ok(variants.remove(0))
}

/// Reusable annotator bound to one loaded assembly/cache.
///
/// Holds a scoped rayon thread pool so concurrent annotators (e.g. one per
/// assembly inside an HTTP service) don't oversubscribe a single global pool.
/// Parallel work inside `annotate_batch` runs through `pool.install(...)`.
///
/// `Clone` is cheap: all state is wrapped in `Arc`. Callers can hand
/// out `Annotator` clones across HTTP request handlers, worker threads, etc.
/// without rewrapping in `Arc<Annotator>` themselves.
#[derive(Clone)]
pub struct Annotator {
    inner: Arc<AnnotatorInner>,
}

struct AnnotatorInner {
    resources: AnnotationResources,
    pool: Arc<ThreadPool>,
}

impl Annotator {
    /// Build an annotator over `resources`.
    ///
    /// Fails when `resources.transcripts` is not fully materialized. This path
    /// annotates inside `pool.install` with no pre-warm pass and `get` never
    /// parses, so a lazily-built index would silently annotate every variant on an
    /// unmaterialized chromosome as `intergenic_variant` and exit 0. Checking once
    /// here (rather than per batch, and rather than in a `debug_assert!` that
    /// release builds drop) is what makes the invariant hold in production.
    /// `load_transcript_indexes` produces a conforming index via `from_loaded`.
    pub fn new(resources: AnnotationResources, pool: Arc<ThreadPool>) -> anyhow::Result<Self> {
        // Verified at construction; the wrapper is dropped because the batch path
        // needs the `Arc` itself, and the check cannot be invalidated afterwards:
        // materialization only ever moves from empty to filled.
        MaterializedTranscriptIndexes::try_from_lazy(Arc::clone(&resources.transcripts))?;
        Ok(Self {
            inner: Arc::new(AnnotatorInner { resources, pool }),
        })
    }

    pub fn from_config(config: &Config) -> anyhow::Result<Self> {
        let pool = build_thread_pool(config.fork)?;
        Self::from_config_with_pool(config, pool)
    }

    pub fn from_config_with_pool(config: &Config, pool: Arc<ThreadPool>) -> anyhow::Result<Self> {
        if config.shift_3prime && config.fasta.is_none() {
            anyhow::bail!("--shift_3prime requires --fasta <path> (indexed FASTA: .fa + .fai)");
        }

        let reference_fasta = config
            .fasta
            .as_deref()
            .map(load_reference_fasta)
            .transpose()?;
        let effects_config = build_effects_config(
            config.distance,
            reference_fasta.clone(),
            config.shift_3prime,
            crate::runner::needs_hgvs_computation(config),
            crate::runner::needs_exon_intron_numbers(config),
        );
        let LoadedPlugins { builtin, dylib } = load_plugins(
            &config.plugin,
            config.dir_plugins.as_deref(),
            reference_fasta,
            resolve_assembly(config.assembly.as_deref()),
        )?;
        let builtin_plugins = Arc::new(builtin);
        let dylib_plugins = Arc::new(dylib);
        let transcripts = if let Some(ref json_cache_path) = config.json_cache {
            load_transcript_indexes(
                json_cache_path,
                &effects_config,
                config.transcript_index_impl,
            )?
        } else {
            Arc::new(LazyTranscriptIndexes::empty())
        };
        let variation_data = if config.check_existing
            || config.af
            || config.af_1kg
            || config.af_gnomade
            || config.af_gnomadg
            || config.max_af
            || config.pubmed
        {
            config
                .json_cache
                .as_deref()
                .map(load_variation_data)
                .transpose()?
        } else {
            None
        };

        Self::new(
            AnnotationResources {
                transcripts,
                variation_data,
                effects_config,
                builtin_plugins,
                dylib_plugins,
                prediction_config: PredictionConfig {
                    sift: config.sift.clone(),
                    polyphen: config.polyphen.clone(),
                },
            },
            pool,
        )
    }

    pub fn resources(&self) -> &AnnotationResources {
        &self.inner.resources
    }

    pub fn pool(&self) -> &Arc<ThreadPool> {
        &self.inner.pool
    }

    pub fn plugin_count(&self) -> usize {
        self.inner.resources.builtin_plugins.plugin_count()
            + self.inner.resources.dylib_plugins.plugin_count()
    }

    pub fn plugin_header_fields(&self) -> Vec<String> {
        self.inner
            .resources
            .builtin_plugins
            .all_header_info()
            .into_iter()
            .map(|(name, _)| name)
            .collect()
    }

    pub fn annotate_batch(&self, batch: &mut Vec<InputVariant>) -> AnnotationStats {
        self.annotate_batch_with_parallelism(batch, true)
    }

    /// Annotate a batch while explicitly controlling intra-batch rayon use.
    ///
    /// Embedders can use this to run many concurrent batches
    /// serially without rebuilding the annotator or changing the pool size.
    /// The work still runs inside the annotator's scoped pool so any nested
    /// rayon code stays on the intended worker set.
    pub fn annotate_batch_with_parallelism(
        &self,
        batch: &mut Vec<InputVariant>,
        use_parallel: bool,
    ) -> AnnotationStats {
        let resources = &self.inner.resources;
        let pool = &self.inner.pool;
        let use_parallel = use_parallel && pool.current_num_threads() > 1;
        // No pre-warm pass: `annotate_batch` calls `transcripts.get` inside the
        // rayon closure and `get` never parses, so the index must already be
        // materialized; `Annotator::new` enforces that, hence no per-batch check.
        pool.install(|| annotate_batch(batch, resources, use_parallel))
    }

    pub fn annotate_variants(&self, batch: &mut Vec<InputVariant>) -> AnnotationStats {
        self.annotate_batch(batch)
    }

    pub fn annotate_variants_with_parallelism(
        &self,
        batch: &mut Vec<InputVariant>,
        use_parallel: bool,
    ) -> AnnotationStats {
        self.annotate_batch_with_parallelism(batch, use_parallel)
    }

    /// Annotate a batch with per-row error semantics.
    ///
    /// Mutates each variant in place (same as [`annotate_batch`]) and returns
    /// a parallel `Vec<Result<(), AnnotationError>>` where the i-th entry maps
    /// to `batch[i]`. The engine is infallible, so every entry is `Ok`; embedding
    /// consumers code against this signature (see [`AnnotationError::Annotation`]
    /// for the reserved error path).
    pub fn try_annotate_batch(
        &self,
        batch: &mut Vec<InputVariant>,
    ) -> Vec<Result<(), AnnotationError>> {
        self.try_annotate_batch_with_parallelism(batch, true)
    }

    pub fn try_annotate_batch_with_parallelism(
        &self,
        batch: &mut Vec<InputVariant>,
        use_parallel: bool,
    ) -> Vec<Result<(), AnnotationError>> {
        let len = batch.len();
        self.annotate_batch_with_parallelism(batch, use_parallel);
        (0..len).map(|_| Ok(())).collect()
    }

    /// Annotate a batch of [`RawRow`]s, producing one annotated `InputVariant`
    /// per row or an [`AnnotationError`] for rows that fail validation/parse.
    ///
    /// This is the high-level per-row entry point intended for embedding consumers
    /// that already speak the simple-allele input shape (chr, pos, ref, alt).
    /// Output is index-aligned with `rows`: `result[i]` corresponds to
    /// `rows[i]`.
    pub fn try_annotate_rows(&self, rows: &[RawRow]) -> Vec<Result<InputVariant, AnnotationError>> {
        self.try_annotate_rows_with_parallelism(rows, true)
    }

    pub fn try_annotate_rows_with_parallelism(
        &self,
        rows: &[RawRow],
        use_parallel: bool,
    ) -> Vec<Result<InputVariant, AnnotationError>> {
        let mut errors: Vec<Option<AnnotationError>> = Vec::with_capacity(rows.len());
        let mut batch: Vec<InputVariant> = Vec::with_capacity(rows.len());
        for row in rows {
            match row_to_input_variant(row) {
                Ok(variant) => {
                    errors.push(None);
                    batch.push(variant);
                }
                Err(err) => errors.push(Some(err)),
            }
        }

        self.annotate_batch_with_parallelism(&mut batch, use_parallel);

        let mut batch_iter = batch.into_iter();
        errors
            .into_iter()
            .map(|maybe_err| match maybe_err {
                Some(err) => Err(err),
                None => Ok(batch_iter
                    .next()
                    .expect("annotated batch length must equal success-row count")),
            })
            .collect()
    }
}

/// Resolve `--assembly` into an [`Assembly`](vep_core::assembly::Assembly), warning
/// when a value was supplied that is not a recognized human assembly.
///
/// Shared by both plugin-loading entry points so they cannot disagree about what a
/// given `--assembly` string means, or about whether an unrecognized one is silent.
pub fn resolve_assembly(raw: Option<&str>) -> Option<vep_core::assembly::Assembly> {
    let assembly = raw.and_then(vep_core::assembly::Assembly::parse);
    if let Some(raw) = raw {
        if assembly.is_none() {
            tracing::warn!(
                "--assembly '{}' is not a recognized human assembly; \
                 assembly-sensitive plugins will run without an assembly hint",
                raw
            );
        }
    }
    assembly
}

/// The plugins one run loads: the built-ins claimed by name, and the C-ABI dylibs
/// the remaining `--plugin` arguments resolved to.
pub(crate) struct LoadedPlugins {
    pub(crate) builtin: BuiltinRegistry,
    pub(crate) dylib: PluginRegistry,
}

/// Load the `--plugin` arguments: each as a built-in when one claims the name,
/// otherwise as a dylib resolved by path or under `dir_plugins`
/// ([`crate::dylib_plugins::load_unresolved`], which skips a name that resolves
/// to no file).
///
/// The assembly is load-bearing rather than advisory: REVEL picks which of its two
/// coordinate columns tabix indexed from it, LoFTEE picks between a bigWig and a
/// tabix-TSV GERP reader, and CADD / dbscSNV use it to refuse a data file whose own
/// declared assembly contradicts the run. Passing `None` leaves every one of those
/// guards inert, so a wrong-assembly file annotates nothing and still exits 0.
pub(crate) fn load_plugins(
    plugin_args: &[String],
    dir_plugins: Option<&str>,
    reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    assembly: Option<vep_core::assembly::Assembly>,
) -> anyhow::Result<LoadedPlugins> {
    let mut builtin = BuiltinRegistry::new();
    if plugin_args.is_empty() {
        return Ok(LoadedPlugins {
            builtin,
            dylib: PluginRegistry::new(),
        });
    }

    let unresolved = builtin
        .load_from_args_with_context(plugin_args, reference_fasta, assembly)
        .map_err(|e| anyhow::anyhow!("plugin loading failed: {}", e))?;
    let dylib = crate::dylib_plugins::load_unresolved(&unresolved, dir_plugins)?;

    Ok(LoadedPlugins { builtin, dylib })
}

pub fn load_reference_fasta(path: &str) -> anyhow::Result<Arc<vep_fasta::IndexedFasta>> {
    info!("Using reference FASTA: {}", path);
    Ok(Arc::new(
        vep_fasta::IndexedFasta::from_path(path)
            .with_context(|| format!("Failed to open indexed FASTA '{path}'"))?,
    ))
}

pub fn build_effects_config(
    distance: (u64, u64),
    reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    shift_3prime: bool,
    compute_hgvs: bool,
    compute_exon_intron_numbers: bool,
) -> Arc<vep_effects::EffectsConfig> {
    Arc::new(vep_effects::EffectsConfig {
        upstream_distance: distance.0,
        downstream_distance: distance.1,
        reference_fasta,
        enable_indel_3prime_shift: shift_3prime,
        // This primarily needs to cover simple repeat-induced shifting around splice windows.
        max_indel_3prime_shift: 2000,
        // LoFTEE context is attached by the main runner path when the plugin is
        // active; this shared helper is not on that path.
        populate_loftee_context: false,
        compute_hgvs,
        compute_exon_intron_numbers,
    })
}

/// Load and index every chromosome up front.
///
/// The lazy path (`LazyTranscriptIndexes::new`) is right for one-shot file runs,
/// which typically touch a fraction of the cache. This eager path is right for
/// a long-lived consumer: it loads once and then serves many requests over
/// arbitrary chromosomes, so deferring per-chromosome parsing only moves the
/// same total cost into the first request that happens to touch each chromosome.
pub fn load_transcript_indexes(
    json_cache_path: &str,
    effects_config: &vep_effects::EffectsConfig,
    index_impl: TranscriptIndexImpl,
) -> anyhow::Result<Arc<LazyTranscriptIndexes>> {
    info!("Loading transcripts from JSON cache: {}", json_cache_path);
    let raw = json_cache::load_all_transcripts(json_cache_path)?;
    info!("Building transcript index with impl={}", index_impl.name());
    Ok(Arc::new(LazyTranscriptIndexes::from_loaded(
        raw,
        index_impl,
        effects_config.upstream_distance,
        effects_config.downstream_distance,
    )))
}

pub fn load_variation_data(
    json_cache_path: &str,
) -> anyhow::Result<Arc<json_cache::VariationData>> {
    info!("Loading variations from JSON cache: {}", json_cache_path);
    let data = json_cache::load_all_variations(json_cache_path)?;
    Ok(Arc::new(data))
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
pub fn annotate_batch(
    batch: &mut Vec<InputVariant>,
    resources: &AnnotationResources,
    use_parallel: bool,
) -> AnnotationStats {
    let batch_len = batch.len();
    let annotation_phase_start = Instant::now();
    let annotate_one = |variant: &mut InputVariant| {
        // Every variant, a breakend spanning more than VEP's 10 Mb `--max_sv_size`
        // included, is annotated against the chromosome's complete transcript set, so
        // one record's output never depends on which records share its batch. This is
        // the same rule as the runner's `annotate_batch`; the two paths must not
        // diverge on a giant breakend.
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

        // The runner's implementation is shared, not duplicated: a second copy
        // here could let the worker and the CLI disagree on the same BND.
        crate::runner::append_bnd_mate_consequences(
            variant,
            resources.transcripts.as_ref(),
            &resources.effects_config,
        );

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

        let is_bnd = variant.variant_class == VariantClass::Translocation;
        if variant.transcript_consequences.is_empty()
            && variant.most_severe_consequence.is_none()
            && resources.transcripts.contains_key(&variant.chr)
            && !is_bnd
        {
            variant.most_severe_consequence = Some(Consequence::IntergenicVariant);
        }
    };

    let stats = if use_parallel {
        let (variants_processed, consequences_calculated) = batch
            .par_iter_mut()
            .map(|variant| {
                annotate_one(variant);
                (1u64, variant.transcript_consequences.len() as u64)
            })
            .reduce(|| (0u64, 0u64), |a, b| (a.0 + b.0, a.1 + b.1));
        AnnotationStats {
            variants_processed,
            consequences_calculated,
        }
    } else {
        let mut stats = AnnotationStats::default();
        for variant in batch.iter_mut() {
            annotate_one(variant);
            stats.variants_processed += 1;
            stats.consequences_calculated += variant.transcript_consequences.len() as u64;
        }
        stats
    };
    let annotation_elapsed = annotation_phase_start.elapsed();
    debug!(
        batch_size = batch_len,
        annotation_ms = annotation_elapsed.as_millis() as u64,
        "batch annotation phase complete"
    );

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
    // No Result to return here: a dylib plugin failure is logged and the batch goes on.
    if let Err(e) = crate::dylib_plugins::run(&resources.dylib_plugins, batch) {
        tracing::error!("dylib plugin execution error: {:#}", e);
    }

    stats
}

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

    transcripts.for_each_overlapping(variant_start, variant_end, |transcript| {
        if let Some(mut tc) = vep_effects::calculate_consequences(variant, transcript, config) {
            if has_predictions && tc.consequences.contains(&Consequence::MissenseVariant) {
                lookup_predictions(&mut tc, transcript, prediction_config);
            }
            variant.transcript_consequences.push(tc);
        }
    });

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

    if let Some(ref mode) = prediction_config.polyphen {
        if let Some(ref matrix) = predictions.polyphen_humvar {
            if let Some(pred) = matrix.lookup(position, alt_aa, AnalysisType::PolyPhenHumVar) {
                tc.polyphen = Some(format_prediction(&pred, mode));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn annotator_over(transcripts: Arc<LazyTranscriptIndexes>) -> anyhow::Result<Annotator> {
        let pool = build_thread_pool(1).unwrap();
        let resources = AnnotationResources {
            transcripts,
            variation_data: None,
            effects_config: Arc::new(vep_effects::EffectsConfig::default()),
            builtin_plugins: Arc::new(BuiltinRegistry::new()),
            dylib_plugins: Arc::new(PluginRegistry::new()),
            prediction_config: PredictionConfig::default(),
        };
        Annotator::new(resources, pool)
    }

    fn empty_annotator() -> Annotator {
        annotator_over(Arc::new(LazyTranscriptIndexes::empty())).unwrap()
    }

    /// The library path annotates a giant breakend against the chromosome's complete
    /// transcript set, whatever else shares its batch, and agrees with the CLI runner
    /// on the same record. VEP restricts such a breakend to the transcripts its
    /// batch-mates loaded; that batch dependence is a VEP defect neither path reproduces.
    #[test]
    fn test_annotator_batch_large_bnd_is_batch_independent() {
        use crate::runner::make_runner_test_transcript;
        use crate::transcript_index::{TranscriptIndex, TranscriptIndexImpl};
        use std::collections::HashMap;

        fn two_transcript_index() -> Arc<LazyTranscriptIndexes> {
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
            Arc::new(LazyTranscriptIndexes::from_indexes(by_chr))
        }
        const BND: &str =
            "21\t999\tbig_bnd\tN\t<BND>\t.\t.\tSVTYPE=BND;END=1000;SVLEN=20000000;CHR2=21;END2=20001000";

        fn run(with_neighbor: bool) -> Vec<(bool, Vec<String>)> {
            let annotator = annotator_over(two_transcript_index()).unwrap();
            let mut batch = Vec::new();
            if with_neighbor {
                batch.extend(parse_vcf_line("21\t1500\tnearby\tA\tG\t.\t.\t.", true).unwrap());
            }
            batch.extend(parse_vcf_line(BND, true).unwrap());
            let _ = annotator.annotate_batch(&mut batch);
            batch
                .iter()
                .filter(|v| v.id.as_deref() == Some("big_bnd"))
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
            single.1.iter().any(|id| id.starts_with("ENSTDISTANT:")),
            "the breakend is annotated against its own span's transcripts: {:?}",
            single.1
        );

        // The two batch paths must agree transcript for transcript on the same record.
        let mut cli: Vec<Vec<String>> =
            crate::runner::annotate_batch_for_test(BND, two_transcript_index());
        for ids in &mut cli {
            ids.sort();
        }
        let mut lib: Vec<Vec<String>> = alone.iter().map(|(_, ids)| ids.clone()).collect();
        lib.sort();
        cli.sort();
        assert_eq!(
            lib, cli,
            "annotator and runner disagree on a giant breakend"
        );
    }

    /// `Annotator::new` must reject a lazily-built index.
    ///
    /// This path fans out to rayon with no pre-warm pass and `get` never parses, so
    /// an unmaterialized chromosome yields `None` while `contains_key` still answers
    /// `true` from the directory listing: every variant on it would be annotated
    /// `intergenic_variant` and the run would exit 0. The check must be a real
    /// error, not a `debug_assert!`: release builds compile assertions out, which
    /// is exactly where the silent mis-annotation would ship.
    #[test]
    fn annotator_rejects_unmaterialized_transcript_index() {
        let dir = tmp_json_cache("annotator_guard");
        let lazy = LazyTranscriptIndexes::new(
            dir.path().to_str().unwrap(),
            TranscriptIndexImpl::Bin,
            5_000,
            5_000,
        )
        .unwrap();
        assert!(
            lazy.contains_key("21"),
            "chromosome reads as present pre-parse; that is the trap"
        );
        assert!(
            lazy.get("21").is_none(),
            "and `get` misses without a prewarm"
        );

        let err = match annotator_over(Arc::new(lazy)) {
            Err(err) => err,
            Ok(_) => panic!("a lazily-built index must not reach the annotator"),
        };
        let msg = format!("{err:#}");
        assert!(
            msg.contains("not fully materialized"),
            "error must name the invariant, got: {msg}"
        );

        let raw = crate::json_cache::load_all_transcripts(dir.path().to_str().unwrap()).unwrap();
        let eager = LazyTranscriptIndexes::from_loaded(raw, TranscriptIndexImpl::Bin, 5_000, 5_000);
        assert!(eager.get("21").is_some());
        annotator_over(Arc::new(eager)).expect("an eagerly-built index is accepted");
    }

    /// The CLI path and the annotator path must produce the same mate annotation.
    ///
    /// Both paths call the runner's `append_bnd_mate_consequences`, so the assertion
    /// is not that two copies agree but that no second copy exists: a divergent copy
    /// (no dedup against locally annotated transcripts, mate annotation through the
    /// generic `annotate_variant` rather than `calculate_paired_mate`, no
    /// `local_span_is_ranged`) is what breaks this equality.
    ///
    /// The inter-chromosomal BND is the case that matters. The variant sits on chr1 and
    /// its mate on chr21, so every consequence in the result is mate-side, which means an
    /// implementation that annotates no mate at all yields an empty set on both sides and
    /// would satisfy a bare equality. Naming the transcript is what prevents that.
    #[test]
    fn annotator_and_cli_agree_on_interchromosomal_bnd_mate_annotation() {
        let dir = tempfile::Builder::new()
            .prefix("vep_annotator_bnd_")
            .tempdir()
            .unwrap();
        for (chr, start) in [("1", 1_000u64), ("21", 25_000_000u64)] {
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
        let path = dir.path().to_str().unwrap();
        let bnd_line = "1\t1300\tbnd_1\tA\tA[21:25001300[\t.\t.\tSVTYPE=BND;MATEID=bnd_2";

        let ids_of = |v: &vep_core::variant::InputVariant| -> Vec<String> {
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
        };

        let eager_raw = crate::json_cache::load_all_transcripts(path).unwrap();
        let eager = Arc::new(LazyTranscriptIndexes::from_loaded(
            eager_raw,
            TranscriptIndexImpl::Bin,
            5_000,
            5_000,
        ));
        let annotator = annotator_over(eager).expect("eager index is accepted");
        let mut worker_batch =
            crate::vcf_parser::parse_vcf_line(bnd_line, true).expect("parse BND");
        let _ = annotator.annotate_batch(&mut worker_batch);

        let lazy = Arc::new(
            LazyTranscriptIndexes::new(path, TranscriptIndexImpl::Bin, 5_000, 5_000).unwrap(),
        );
        let cli_ids = crate::runner::annotate_batch_for_test(bnd_line, lazy);

        let worker_ids: Vec<Vec<String>> = worker_batch.iter().map(ids_of).collect();
        assert_eq!(
            worker_ids, cli_ids,
            "the annotator and CLI paths must agree on an inter-chromosomal BND mate"
        );
        assert!(
            worker_ids
                .iter()
                .any(|ids| ids.iter().any(|s| s.starts_with("ENST21:"))),
            "the chr21 mate transcript must be annotated, or this equality compares two \
             empty sets: {worker_ids:?}"
        );
    }

    /// Minimal one-chromosome JSON cache in a self-cleaning temp dir.
    fn tmp_json_cache(tag: &str) -> tempfile::TempDir {
        let dir = tempfile::Builder::new()
            .prefix(&format!("vep_annotator_{tag}_"))
            .tempdir()
            .unwrap();
        let chr_dir = dir.path().join("transcripts").join("21");
        std::fs::create_dir_all(&chr_dir).unwrap();
        std::fs::write(
            chr_dir.join("0-1000.json"),
            r#"[{"stable_id":"ENST1","gene_stable_id":"ENSG1","chr":"21","start":100,"end":200,"strand":1,"biotype":"protein_coding","source":"ensembl"}]"#,
        )
        .unwrap();
        dir
    }

    #[test]
    fn validate_raw_row_rejects_zero_pos() {
        let row = RawRow {
            chrom: "21".into(),
            pos: 0,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
        };
        let err = validate_raw_row(&row).unwrap_err();
        assert!(matches!(err, AnnotationError::InvalidInput(_)));
        assert!(err.to_string().contains("1-based"));
    }

    #[test]
    fn validate_raw_row_rejects_symbolic_alt() {
        let row = RawRow {
            chrom: "21".into(),
            pos: 1,
            ref_allele: "A".into(),
            alt_allele: "<DEL>".into(),
        };
        let err = validate_raw_row(&row).unwrap_err();
        assert!(err.to_string().contains("symbolic"));
    }

    #[test]
    fn try_annotate_rows_preserves_input_order_and_errors() {
        let annotator = empty_annotator();
        let rows = vec![
            RawRow {
                chrom: "21".into(),
                pos: 100,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
            },
            RawRow {
                chrom: "".into(),
                pos: 100,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
            },
            RawRow {
                chrom: "21".into(),
                pos: 0,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
            },
            RawRow {
                chrom: "X".into(),
                pos: 200,
                ref_allele: "C".into(),
                alt_allele: "T".into(),
            },
        ];

        let results = annotator.try_annotate_rows(&rows);
        assert_eq!(results.len(), 4);
        assert!(results[0].is_ok());
        assert!(matches!(
            results[1].as_ref().unwrap_err(),
            AnnotationError::InvalidInput(_)
        ));
        assert!(matches!(
            results[2].as_ref().unwrap_err(),
            AnnotationError::InvalidInput(_)
        ));
        assert!(results[3].is_ok());
    }

    #[test]
    fn try_annotate_batch_returns_one_result_per_input() {
        let annotator = empty_annotator();
        let mut batch: Vec<InputVariant> = vec![
            InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec()),
            InputVariant::new("X".into(), 200, 200, b"C".to_vec(), b"T".to_vec()),
        ];
        let results = annotator.try_annotate_batch(&mut batch);
        assert_eq!(results.len(), 2);
        assert!(results.iter().all(|r| r.is_ok()));
    }

    #[test]
    fn explicit_serial_annotation_preserves_batch_contract() {
        let annotator = empty_annotator();
        let mut batch: Vec<InputVariant> = vec![
            InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec()),
            InputVariant::new("X".into(), 200, 200, b"C".to_vec(), b"T".to_vec()),
        ];

        let stats = annotator.annotate_batch_with_parallelism(&mut batch, false);

        assert_eq!(stats.variants_processed, 2);
        assert_eq!(batch.len(), 2);
        assert_eq!(batch[0].chr, "21");
        assert_eq!(batch[1].chr, "X");
    }
}
