// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reusable annotation engine for file-based and service-based entrypoints.

use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;

use anyhow::Context;
use rayon::prelude::*;
use rayon::ThreadPool;
use thiserror::Error;
use tracing::{debug, info};

use crate::config::{CliError, Config};
use crate::filter_set::FilterSet;
use crate::json_cache;
use crate::nearest::NearestType;
use crate::sv_source::{self, SvSource};
use crate::transcript_index::{
    LazyTranscriptIndexes, MaterializedTranscriptIndexes, TranscriptSet,
};
use crate::variation_matcher::{self, ColocationRules};
use crate::variation_source::{ReadFailure, VariationSource};
use crate::vcf_parser::parse_vcf_line;
use vep_builtins::BuiltinRegistry;
use vep_core::consequence::Consequence;
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::{InputVariant, VariantClass};
use vep_plugin::PluginRegistry;

/// The refusal of `--hgvs` or `--hgvsg` with no reference FASTA, Ensembl VEP's
/// text (Runner.pm 759), under the `fasta` component.
pub(crate) const HGVS_NEEDS_FASTA: &str = "Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode without a FASTA file (see --fasta)";

/// Shared annotation resources that can be reused across batches.
pub struct AnnotationResources {
    pub transcripts: Arc<LazyTranscriptIndexes>,
    pub variation_source: Option<Arc<VariationSource>>,
    /// The known structural variants of `--check_svs <file>`, `None` without
    /// the flag.
    pub sv_source: Option<Arc<SvSource>>,
    /// Which known variants of a window co-locate with an input variant.
    pub colocation: ColocationRules,
    pub effects_config: Arc<vep_effects::EffectsConfig>,
    pub builtin_plugins: Arc<BuiltinRegistry>,
    /// C-ABI plugins loaded from the `--plugin` arguments no built-in claimed.
    pub dylib_plugins: Arc<PluginRegistry>,
    pub prediction_config: PredictionConfig,
    /// The identifier `--nearest` reports for each record's nearest
    /// protein-coding transcription start sites; `None` leaves
    /// [`InputVariant::nearest`] unset.
    pub nearest: Option<NearestType>,
    /// The compiled `--transcript_filter`, which every transcript passes before
    /// its consequences are computed, as in the `vep` binary; `None` filters
    /// nothing.
    pub transcript_filter: Option<Arc<FilterSet>>,
    /// The chromosome synonyms every input chromosome is resolved through
    /// before the cache lookups, as the `vep` binary resolves them: Ensembl VEP
    /// fetches a record's regions, and a breakend mate's, by the name
    /// `get_source_chr_name` gives (`AnnotationSource.pm:243-245` and
    /// `:251-254`, through `:179`; `BaseVEP.pm:689-738`). An empty table still
    /// resolves a `chr` prefix.
    pub synonyms: Arc<ChromosomeSynonyms>,
}

/// Build a scoped rayon thread pool sized for one annotator.
///
/// `rayon::ThreadPoolBuilder::build_global()` is avoided: concurrent embedding
/// consumers (long-lived services, request handlers) would have every in-flight
/// request fighting for one shared global pool of N threads, oversubscribing the
/// CPU. With a scoped pool each annotator's parallel work runs inside
/// `pool.install(...)` and stays bounded to its configured fork count.
///
/// The workers are `vep-annot-<i>`: Linux keeps 15 bytes of a thread name, so
/// the prefix leaves room for the index of every worker a machine can run.
pub fn build_thread_pool(num_threads: usize) -> anyhow::Result<Arc<ThreadPool>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .thread_name(|i| format!("vep-annot-{i}"))
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
    /// Read PolyPhen from the HumDiv matrix rather than HumVar (`--humdiv`).
    pub humdiv: bool,
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
    /// A variation file of the cache or the `--check_svs` file could not be
    /// read, so the batch's co-located or overlapping known variants are
    /// unknown; every row of that batch carries it.
    #[error("annotation engine error: {0}")]
    Annotation(String),
}

/// A simple-allele row: chr, 1-based pos, ref, alt, one allele per row.
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
/// Holds a scoped rayon thread pool so concurrent annotators (one per assembly,
/// say) don't oversubscribe a single global pool. Parallel work inside
/// `annotate_batch` runs through `pool.install(...)`.
///
/// `Clone` is cheap: all state is wrapped in `Arc`. Callers can hand out
/// `Annotator` clones across threads without rewrapping in `Arc<Annotator>`
/// themselves.
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

    /// Resolves the cache and the FASTA the configuration names and settles the
    /// HGVS flags, as the `vep` binary does, so a configuration with no cache is
    /// an error here too rather than an empty index, `--fasta_dir` is read, and
    /// `--hgvs` without a FASTA is refused.
    pub fn from_config_with_pool(config: &Config, pool: Arc<ThreadPool>) -> anyhow::Result<Self> {
        let cache = crate::cache_resolve::resolve(config)?;
        // After the cache and before the notices and the FASTA, as the binary
        // orders them (`Runner::new` resolves the cache; `Runner::run` compiles
        // the filter, prints the configuration notices, then opens the FASTA and
        // runs its checks): with a missing cache, a bad filter and a missing
        // FASTA together, both paths report them in that order.
        let transcript_filter =
            crate::runner::compile_transcript_filter(&config.transcript_filter)?;
        log_notices(&config.notices);
        let predictions = crate::runner::prediction_checks(config, &cache.info)?;
        log_notices(&predictions.notices);
        let reference_fasta = crate::runner::open_fasta(config, &cache)?;
        let checks = crate::runner::fasta_setup_checks(config, reference_fasta.is_some())?;
        log_notices(&checks.notices);

        let mut effects = new_effects_config(
            config.distance,
            reference_fasta.clone(),
            config.shift_3prime,
            hgvs_options(config),
            checks.hgvs,
            crate::runner::needs_exon_intron_numbers(config),
            config.protein_version,
        );
        effects.hgvsg = hgvsg_options(config, &cache.synonyms);
        effects.spdi = spdi_options(config, &cache.synonyms);
        effects.domains = config.domains;
        effects.use_transcript_ref = config.use_transcript_ref;
        effects.process_ref_homs = config.process_ref_homs;
        let effects_config = Arc::new(effects);
        let LoadedPlugins { builtin, dylib } = load_plugins(
            &config.plugin,
            &config.dir_plugins,
            config.safe,
            reference_fasta,
            resolve_assembly(cache.assembly.as_deref()),
        )?;
        let builtin_plugins = Arc::new(builtin);
        let dylib_plugins = Arc::new(dylib);
        let transcripts = load_transcript_indexes(
            &cache.dir,
            &effects_config,
            TranscriptSet::from_config(config),
        )?;
        let variation_source = crate::runner::open_variation_source(config, &cache)?;
        let sv_source = crate::runner::open_sv_source(config)?;

        Self::new(
            AnnotationResources {
                transcripts,
                variation_source,
                sv_source,
                colocation: ColocationRules::from_config(config),
                effects_config,
                builtin_plugins,
                dylib_plugins,
                prediction_config: PredictionConfig {
                    sift: predictions.sift,
                    polyphen: predictions.polyphen,
                    humdiv: config.humdiv,
                },
                nearest: config.nearest,
                transcript_filter,
                synonyms: Arc::clone(&cache.synonyms),
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
        let (stats, failure) = self.run_batch(batch, use_parallel);
        // This entry point returns statistics, not a result, so a read failure
        // of a variation or structural-variant file is logged under the error
        // contract's component; `try_annotate_batch` returns it per row.
        if let Some((component, message)) = failure {
            tracing::error!("{}", CliError::new(component, message));
        }
        stats
    }

    /// Annotates the batch on the annotator's pool; the second value is the
    /// batch's read failure of a variation or structural-variant file, if one
    /// occurred, with the error component it is reported under.
    fn run_batch(
        &self,
        batch: &mut Vec<InputVariant>,
        use_parallel: bool,
    ) -> (AnnotationStats, Option<(&'static str, String)>) {
        let resources = &self.inner.resources;
        let pool = &self.inner.pool;
        let use_parallel = use_parallel && pool.current_num_threads() > 1;
        let failure = ReadFailure::default();
        // No pre-warm pass: `annotate_batch` calls `transcripts.get` inside the
        // rayon closure and `get` never parses, so the index must already be
        // materialized; `Annotator::new` enforces that, hence no per-batch check.
        let stats = pool.install(|| annotate_batch(batch, resources, use_parallel, &failure));
        (stats, failure.take_with_component())
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
    /// to `batch[i]`. Every entry is `Ok` unless a variation file of the cache
    /// or the `--check_svs` file could not be read, when every entry carries
    /// that failure as [`AnnotationError::Annotation`]: the rows were annotated
    /// without their co-located or overlapping known variants.
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
        let (_, failure) = self.run_batch(batch, use_parallel);
        let row_result = match failure {
            Some((_component, message)) => Err(AnnotationError::Annotation(message)),
            None => Ok(()),
        };
        vec![row_result; len]
    }

    /// Annotate a batch of [`RawRow`]s, producing one annotated `InputVariant`
    /// per row or an [`AnnotationError`] for rows that fail validation/parse,
    /// and for every parsed row when a variation file of the cache or the
    /// `--check_svs` file could not be read.
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

        let (_, read_failure) = self.run_batch(&mut batch, use_parallel);

        let mut batch_iter = batch.into_iter();
        errors
            .into_iter()
            .map(|maybe_err| match maybe_err {
                Some(err) => Err(err),
                None => {
                    let variant = batch_iter
                        .next()
                        .expect("annotated batch length must equal success-row count");
                    match &read_failure {
                        Some((_component, message)) => {
                            Err(AnnotationError::Annotation(message.clone()))
                        }
                        None => Ok(variant),
                    }
                }
            })
            .collect()
    }
}

/// Logs the `INFO:` lines the binary prints to stderr; the tracing level
/// stands in for the prefix.
fn log_notices(notices: &[&str]) {
    for notice in notices {
        info!("{}", notice.trim_start_matches("INFO: "));
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

/// Load the `--plugin` arguments as `get_all_Plugins` does (Runner.pm:901-960):
/// a plugins directory that does not exist is warned about when a plugin is named,
/// an error under `safe` (`:911-918`); each argument is `Name,param,param`
/// (`:928`), loaded as a built-in when one claims the name, otherwise as a dylib
/// resolved by path or under the directory ([`crate::dylib_plugins::load_unresolved`]);
/// a name that resolves to neither, or a dylib that does not open, is warned
/// about and dropped, an error under `safe` (`:929-939`), as is a built-in or a
/// dylib whose `init` fails (`Failed to instantiate plugin <Name>: <reason>`,
/// `:946-954`). The arguments load one at a time in their order, so the
/// warnings follow it as Ensembl VEP's single loop writes them (`:925`).
///
/// The assembly is load-bearing rather than advisory: REVEL picks which of its two
/// coordinate columns tabix indexed from it, LoFTEE picks between a bigWig and a
/// tabix-TSV GERP reader, and CADD / dbscSNV use it to refuse a data file whose own
/// declared assembly contradicts the run. Passing `None` leaves every one of those
/// guards inert, so a wrong-assembly file annotates nothing and still exits 0.
pub(crate) fn load_plugins(
    plugin_args: &[String],
    dir_plugins: &str,
    safe: bool,
    reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    assembly: Option<vep_core::assembly::Assembly>,
) -> anyhow::Result<LoadedPlugins> {
    let mut builtin = BuiltinRegistry::new();
    let mut dylib = PluginRegistry::new();
    if plugin_args.is_empty() {
        return Ok(LoadedPlugins { builtin, dylib });
    }

    if !std::path::Path::new(dir_plugins).is_dir() {
        plugin_warning(
            format!("Plugins directory '{dir_plugins}' not found."),
            safe,
        )?;
    }

    for arg in plugin_args {
        let outcome = builtin.load_from_args_with_context(
            std::slice::from_ref(arg),
            reference_fasta.clone(),
            assembly,
        );
        for (name, error) in &outcome.failed {
            // `Init`'s `Display` adds a prefix Ensembl VEP's line does not carry, so
            // its text alone is the reason.
            let reason = match error {
                vep_builtins::PluginError::Init(text) => text.clone(),
                other => other.to_string(),
            };
            plugin_warning(
                format!("Failed to instantiate plugin {name}: {reason}"),
                safe,
            )?;
        }
        for unresolved in &outcome.unresolved {
            crate::dylib_plugins::load_unresolved(&mut dylib, unresolved, dir_plugins, safe)?;
        }
    }

    Ok(LoadedPlugins { builtin, dylib })
}

/// A plugin-loading warning (Runner.pm:911-918, `:929-939`, `:946-954`): the
/// plugin component's error carrying `message` under `safe`, else the warning,
/// and the run goes on.
pub(crate) fn plugin_warning(message: String, safe: bool) -> anyhow::Result<()> {
    if safe {
        return Err(crate::config::CliError::new("plugin", message).into());
    }
    crate::warnings::global().warn(&message);
    Ok(())
}

/// The effects configuration of an embedding caller, without HGVSg; `hgvs`
/// carries the HGVS shift, the two HGVSp forms and the `--shift_3prime` the
/// HGVSp predicates read, which agrees with `shift_3prime` here.
pub fn build_effects_config(
    distance: (u64, u64),
    reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    shift_3prime: bool,
    hgvs: vep_effects::hgvs::HgvsOptions,
    compute_hgvs: bool,
    compute_exon_intron_numbers: bool,
    protein_version: bool,
) -> Arc<vep_effects::EffectsConfig> {
    Arc::new(new_effects_config(
        distance,
        reference_fasta,
        shift_3prime,
        hgvs,
        compute_hgvs,
        compute_exon_intron_numbers,
        protein_version,
    ))
}

/// The effects configuration without HGVSg or SPDI; `Annotator::from_config_with_pool`
/// sets `hgvsg` and `spdi` on the result.
fn new_effects_config(
    distance: (u64, u64),
    reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    shift_3prime: bool,
    hgvs: vep_effects::hgvs::HgvsOptions,
    compute_hgvs: bool,
    compute_exon_intron_numbers: bool,
    protein_version: bool,
) -> vep_effects::EffectsConfig {
    vep_effects::EffectsConfig {
        upstream_distance: distance.0,
        downstream_distance: distance.1,
        reference_fasta,
        enable_indel_3prime_shift: shift_3prime,
        hgvs,
        // LoFTEE context is attached by the main runner path when the plugin is
        // active; this shared helper is not on that path.
        populate_loftee_context: false,
        compute_hgvs,
        compute_exon_intron_numbers,
        protein_version,
        hgvsg: None,
        spdi: None,
        domains: false,
        use_transcript_ref: false,
        process_ref_homs: false,
    }
}

/// The HGVSc and HGVSp options of a configuration: the shift under `--shift_hgvs`,
/// the two HGVSp forms and whether `--shift_3prime 1` moved the consequences.
pub(crate) fn hgvs_options(config: &Config) -> vep_effects::hgvs::HgvsOptions {
    vep_effects::hgvs::HgvsOptions {
        shift_hgvs: config.shift_hgvs,
        hgvsp_use_prediction: config.hgvsp_use_prediction,
        remove_hgvsp_version: config.remove_hgvsp_version,
        shift_3prime: config.shift_3prime,
    }
}

/// The HGVSg options of a configuration: `None` without `--hgvsg`. Insertions
/// and deletions are described at their most 3' position unless `--shift_hgvs 0`,
/// the one setting that gates every HGVS shift (`shift_hgvs_variants_3prime`,
/// VariationFeature.pm 1981-1988); the chromosome is named as the input names
/// it, or under `--hgvsg_use_accession` by the first `NC_` name among the input
/// name's own synonyms (OutputFactory.pm 1367-1377).
pub(crate) fn hgvsg_options(
    config: &Config,
    synonyms: &ChromosomeSynonyms,
) -> Option<vep_effects::hgvsg::HgvsgOptions> {
    config.hgvsg.then(|| vep_effects::hgvsg::HgvsgOptions {
        shift: config.shift_hgvs,
        accessions: config.hgvsg_use_accession.then(|| nc_accessions(synonyms)),
    })
}

/// The SPDI options of a configuration: `None` without `--spdi` or `--ga4gh_vrs`;
/// the GA4GH VRS allele is built under `--ga4gh_vrs` from the `NC_` accession
/// each chromosome's synonyms carry.
pub(crate) fn spdi_options(
    config: &Config,
    synonyms: &ChromosomeSynonyms,
) -> Option<vep_effects::spdi::SpdiOptions> {
    (config.spdi || config.ga4gh_vrs).then(|| vep_effects::spdi::SpdiOptions {
        vrs_accessions: config.ga4gh_vrs.then(|| nc_accessions(synonyms)),
    })
}

/// Every name of the synonyms table with an `NC_` synonym, mapped to the first
/// such synonym in the table's order ([`vep_effects::hgvsg::nc_accession`]).
fn nc_accessions(synonyms: &ChromosomeSynonyms) -> HashMap<String, String> {
    synonyms
        .names()
        .filter_map(|name| {
            let accession = vep_effects::hgvsg::nc_accession(
                synonyms.synonyms_of(name).iter().map(String::as_str),
            )?;
            Some((name.to_string(), accession.to_string()))
        })
        .collect()
}

/// Load and index every chromosome up front, keeping the transcripts of
/// `transcript_set` alone.
///
/// The lazy path (`LazyTranscriptIndexes::new`) is right for one-shot file runs,
/// which typically touch a fraction of the cache. This eager path is right for
/// a long-lived consumer: it loads once and then serves many requests over
/// arbitrary chromosomes, so deferring per-chromosome parsing only moves the
/// same total cost into the first request that happens to touch each chromosome.
pub fn load_transcript_indexes(
    json_cache_path: &str,
    effects_config: &vep_effects::EffectsConfig,
    transcript_set: TranscriptSet,
) -> anyhow::Result<Arc<LazyTranscriptIndexes>> {
    info!("Loading transcripts from JSON cache: {}", json_cache_path);
    let mut raw = json_cache::load_all_transcripts(json_cache_path)?;
    for transcripts in raw.values_mut() {
        transcripts.retain(|transcript| transcript_set.accepts(transcript));
    }
    info!("Building transcript index");
    Ok(Arc::new(LazyTranscriptIndexes::from_loaded(
        raw,
        effects_config.upstream_distance,
        effects_config.downstream_distance,
    )))
}

/// Annotate a batch of variants, optionally in parallel using rayon.
///
/// A failure to read a variation file of the cache or the `--check_svs` file
/// lands in `failure`, the caller's per-batch slot: this function returns
/// statistics, and the lookups inside the closure have no error channel of
/// their own. A caller that drops the slot unread annotates the batch without
/// its co-located or overlapping known variants and never learns it.
pub fn annotate_batch(
    batch: &mut Vec<InputVariant>,
    resources: &AnnotationResources,
    use_parallel: bool,
    failure: &ReadFailure,
) -> AnnotationStats {
    let batch_len = batch.len();
    let annotation_phase_start = Instant::now();
    for variant in batch.iter_mut() {
        crate::pipeline::name_source_chromosomes(
            variant,
            &resources.synonyms,
            &resources.transcripts,
        );
    }
    if let Some(source) = &resources.variation_source {
        if let Err(error) = source.prefetch(batch.iter().map(|v| (v.chr.as_str(), v.start, v.end)))
        {
            failure.record(format!("{error:#}"));
        }
    }
    if let Some(source) = &resources.sv_source {
        if let Err(error) = source.prefetch(sv_source::batch_spans(batch)) {
            failure.record_under("check_svs", format!("{error:#}"));
        }
    }
    let annotate_one = |variant: &mut InputVariant| {
        crate::nearest::annotate(variant, resources.transcripts.as_ref(), resources.nearest);
        if variant.no_allele {
            return;
        }
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
        let transcript_filter = resources.transcript_filter.as_deref();
        if !is_interchrom_paired_breakend {
            if let Some(chr_transcripts) = resources.transcripts.get(&variant.chr) {
                crate::runner::annotate_variant(
                    variant,
                    chr_transcripts,
                    &resources.effects_config,
                    &resources.prediction_config,
                    transcript_filter,
                );
            }
        }

        // The runner's implementation is shared, not duplicated: a second copy
        // here could let the two entry points disagree on the same BND.
        crate::runner::append_bnd_mate_consequences(
            variant,
            resources.transcripts.as_ref(),
            &resources.effects_config,
            transcript_filter,
        );

        if let Some(source) = &resources.variation_source {
            variation_matcher::annotate_colocated(variant, source, &resources.colocation, failure);
        }
        if let Some(source) = &resources.sv_source {
            sv_source::annotate_overlapping_svs(variant, source, failure);
        }

        // A breakend's intergenic verdict is its record's: `apply_record_rules`.
        let is_bnd = variant.variant_class == VariantClass::Translocation;
        if variant.transcript_consequences.is_empty()
            && variant.most_severe_consequence.is_none()
            && resources.transcripts.contains_key(&variant.chr)
            && !is_bnd
        {
            variant.most_severe_consequence = Some(Consequence::IntergenicVariant);
            variant.intergenic_shift =
                crate::runner::intergenic_shift(variant, &resources.effects_config);
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
    crate::runner::apply_record_rules(batch, Some(resources.transcripts.as_ref()));
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

#[cfg(test)]
mod tests {
    use super::*;

    fn annotator_over(transcripts: Arc<LazyTranscriptIndexes>) -> anyhow::Result<Annotator> {
        annotator_with_nearest(transcripts, None)
    }

    fn annotator_with_nearest(
        transcripts: Arc<LazyTranscriptIndexes>,
        nearest: Option<NearestType>,
    ) -> anyhow::Result<Annotator> {
        let pool = build_thread_pool(1).unwrap();
        let resources = AnnotationResources {
            transcripts,
            variation_source: None,
            sv_source: None,
            colocation: ColocationRules::default(),
            effects_config: Arc::new(vep_effects::EffectsConfig::default()),
            builtin_plugins: Arc::new(BuiltinRegistry::new()),
            dylib_plugins: Arc::new(PluginRegistry::new()),
            prediction_config: PredictionConfig::default(),
            nearest,
            transcript_filter: None,
            synonyms: Arc::new(ChromosomeSynonyms::new()),
        };
        Annotator::new(resources, pool)
    }

    fn empty_annotator() -> Annotator {
        annotator_over(Arc::new(LazyTranscriptIndexes::empty())).unwrap()
    }

    /// Every pool thread carries a distinct `vep-annot-<i>` name that fits the
    /// 15 bytes Linux keeps of a thread name, so a per-thread profile of a
    /// 16-worker run attributes each worker separately.
    #[test]
    fn pool_threads_carry_distinct_names_within_the_comm_limit() {
        let pool = build_thread_pool(16).unwrap();
        let mut names: Vec<String> =
            pool.broadcast(|_| std::thread::current().name().unwrap().to_string());
        names.sort();
        names.dedup();
        assert_eq!(names.len(), 16);
        for name in &names {
            assert!(name.starts_with("vep-annot-"), "{name}");
            assert!(name.len() <= 15, "{name}");
        }
    }

    /// The library path annotates a giant breakend against the chromosome's complete
    /// transcript set, whatever else shares its batch, and agrees with the CLI runner
    /// on the same record. VEP restricts such a breakend to the transcripts its
    /// batch-mates loaded; that batch dependence is a VEP defect neither path reproduces.
    #[test]
    fn test_annotator_batch_large_bnd_is_batch_independent() {
        use crate::runner::make_runner_test_transcript;
        use crate::transcript_index::TranscriptIndex;
        use std::collections::HashMap;

        fn two_transcript_index() -> Arc<LazyTranscriptIndexes> {
            let buffered_tx =
                make_runner_test_transcript("ENSTBUFFERED", "ENSGBUFFERED", "21", 1_000);
            let distant_tx =
                make_runner_test_transcript("ENSTDISTANT", "ENSGDISTANT", "21", 20_000_000);
            let mut by_chr = HashMap::new();
            by_chr.insert(
                "21".to_string(),
                TranscriptIndex::new(vec![buffered_tx, distant_tx], 5_000, 5_000),
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
        let lazy = LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap();
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
        let eager = LazyTranscriptIndexes::from_loaded(raw, 5_000, 5_000);
        assert!(eager.get("21").is_some());
        annotator_over(Arc::new(eager)).expect("an eagerly-built index is accepted");
    }

    /// The CLI path and the annotator path must produce the same mate annotation.
    ///
    /// Both paths call the runner's `append_bnd_mate_consequences`, so the assertion
    /// is not that two copies agree but that no second copy exists: a divergent copy
    /// (no dedup against locally annotated transcripts, mate annotation through the
    /// generic `annotate_variant` rather than `calculate_paired_mate`) is what breaks
    /// this equality.
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
        let eager = Arc::new(LazyTranscriptIndexes::from_loaded(eager_raw, 5_000, 5_000));
        let annotator = annotator_over(eager).expect("eager index is accepted");
        let mut worker_batch =
            crate::vcf_parser::parse_vcf_line(bnd_line, true).expect("parse BND");
        let _ = annotator.annotate_batch(&mut worker_batch);

        let lazy = Arc::new(LazyTranscriptIndexes::new(path, 5_000, 5_000).unwrap());
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

    /// An annotator with no transcripts over the region files under
    /// `<dir>/variations/<chr>/`, read by the five-column layout, on a pool of
    /// `threads`.
    fn annotator_over_variation_files(dir: &std::path::Path, threads: usize) -> Annotator {
        let info = vep_core::cache_info::CacheInfo {
            variation_cols: ["variation_name", "failed", "start", "end", "allele_string"]
                .iter()
                .map(|name| name.to_string())
                .collect(),
            ..Default::default()
        };
        let source = crate::variation_source::VariationSource::open(
            dir,
            &info,
            false,
            Arc::new(vep_core::cache_manifest::CacheManifest::absent(dir)),
        )
        .unwrap()
        .expect("a region file is a variation source");
        let resources = AnnotationResources {
            transcripts: Arc::new(LazyTranscriptIndexes::empty()),
            variation_source: Some(Arc::new(source)),
            sv_source: None,
            colocation: ColocationRules::default(),
            effects_config: Arc::new(vep_effects::EffectsConfig::default()),
            builtin_plugins: Arc::new(BuiltinRegistry::new()),
            dylib_plugins: Arc::new(PluginRegistry::new()),
            prediction_config: PredictionConfig::default(),
            nearest: None,
            transcript_filter: None,
            synonyms: Arc::new(ChromosomeSynonyms::new()),
        };
        Annotator::new(resources, build_thread_pool(threads).unwrap()).unwrap()
    }

    /// A variation file that cannot be read reaches the per-row API as an
    /// error on every row of the batch; the statistics API logs it instead.
    #[test]
    fn try_annotate_batch_reports_a_variation_read_failure_on_every_row() {
        let dir = tempfile::tempdir().unwrap();
        let region = dir
            .path()
            .join("variations")
            .join("21")
            .join("1-1000000_var.gz");
        std::fs::create_dir_all(region.parent().unwrap()).unwrap();
        std::fs::write(&region, "not a gzip body").unwrap();
        let annotator = annotator_over_variation_files(dir.path(), 1);
        let mut batch: Vec<InputVariant> = vec![
            InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec()),
            InputVariant::new("21".into(), 200, 200, b"C".to_vec(), b"T".to_vec()),
        ];
        let results = annotator.try_annotate_batch(&mut batch);
        assert_eq!(results.len(), 2);
        for result in &results {
            match result {
                Err(AnnotationError::Annotation(message)) => assert!(
                    message.starts_with(&format!("reading {}: ", region.display())),
                    "{message}"
                ),
                other => panic!("expected the read failure, got {other:?}"),
            }
        }
        let rows = [RawRow {
            chrom: "21".into(),
            pos: 100,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
        }];
        assert!(matches!(
            annotator.try_annotate_rows(&rows).as_slice(),
            [Err(AnnotationError::Annotation(_))]
        ));
        // The statistics API returns no result, so it logs the failure instead.
        annotator.annotate_batch(&mut batch);
    }

    /// The slot a batch records its read failure in is that batch's own: two
    /// clones annotating at once, one batch on a corrupt region file and the
    /// other on a sound one, each report their own result, whatever the
    /// interleaving of their prefetches and lookups on the shared source.
    #[test]
    fn concurrent_batches_report_only_their_own_variation_read_failure() {
        let dir = tempfile::tempdir().unwrap();
        let chromosome = dir.path().join("variations").join("21");
        std::fs::create_dir_all(&chromosome).unwrap();
        let corrupt = chromosome.join("1-1000000_var.gz");
        std::fs::write(&corrupt, "not a gzip body").unwrap();
        let mut sound = flate2::write::GzEncoder::new(
            std::fs::File::create(chromosome.join("1000001-2000000_var.gz")).unwrap(),
            flate2::Compression::default(),
        );
        std::io::Write::write_all(&mut sound, b"rs2 0 1500000 1500000 A/G\n").unwrap();
        sound.finish().unwrap();
        let annotator = annotator_over_variation_files(dir.path(), 4);
        let snv = |pos: u64| InputVariant::new("21".into(), pos, pos, b"A".to_vec(), b"G".to_vec());
        for _ in 0..50 {
            let barrier = Arc::new(std::sync::Barrier::new(2));
            let run = |positions: Vec<u64>| {
                let annotator = annotator.clone();
                let barrier = Arc::clone(&barrier);
                std::thread::spawn(move || {
                    let mut batch: Vec<InputVariant> = positions.into_iter().map(snv).collect();
                    barrier.wait();
                    let results = annotator.try_annotate_batch(&mut batch);
                    (results, batch)
                })
            };
            let (on_corrupt, on_sound) = (run(vec![500, 700]), run(vec![1_500_000]));
            let (results, _) = on_corrupt.join().unwrap();
            assert!(
                results.iter().all(|result| matches!(
                    result,
                    Err(AnnotationError::Annotation(message))
                        if message.starts_with(&format!("reading {}: ", corrupt.display()))
                )),
                "{results:?}"
            );
            let (results, batch) = on_sound.join().unwrap();
            assert!(results.iter().all(Result::is_ok), "{results:?}");
            assert_eq!(batch[0].existing_variation, ["rs2"]);
        }
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

    /// The library and the command-line paths annotate through one function, so a
    /// record whose rows depend on a path-specific step cannot differ between them:
    /// a one-base `<*>` deletion with `SVTYPE=DEL` inside ENST00000390690's coding
    /// sequence is `frameshift_variant, feature_truncation`, where the sequence-allele
    /// reading of a bare `<*>` would write `coding_sequence_variant`; the `<DEL>`
    /// record with the same INFO is the control. Both paths run the records against
    /// the release 116 corpus cache.
    #[test]
    fn the_library_path_annotates_through_the_command_line_path() {
        const EXPECTED_ROW: &str = "ENST00000390690:frameshift_variant,feature_truncation";
        let cache = format!(
            "{}/../../tests/golden/116/GRCh38-release116/json_cache",
            env!("CARGO_MANIFEST_DIR")
        );
        let raw = crate::json_cache::load_all_transcripts(&cache).unwrap();
        let transcripts = Arc::new(LazyTranscriptIndexes::from_loaded(raw, 5_000, 5_000));
        let annotator = annotator_over(Arc::clone(&transcripts)).unwrap();
        let rows = |variant: &InputVariant| -> Vec<String> {
            let mut ids: Vec<String> = variant
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
        for record in [
            "21\t30425650\tstar_del1\tN\t<*>\t.\t.\tSVTYPE=DEL;END=30425651",
            "21\t30425650\tdel1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=30425651",
        ] {
            let mut library_batch = parse_vcf_line(record, true).unwrap();
            let _ = annotator.annotate_batch(&mut library_batch);
            let library: Vec<Vec<String>> = library_batch.iter().map(rows).collect();
            let command_line =
                crate::runner::annotate_batch_for_test(record, Arc::clone(&transcripts));
            assert_eq!(library, command_line, "{record}");
            assert!(
                library[0].iter().any(|row| row == EXPECTED_ROW),
                "{record}: {:?}",
                library[0]
            );
        }
    }

    /// The library path sets NEAREST on a record whose only ALT is `*`, the
    /// fixture and value of `runner::tests::nearest_is_set_on_a_record_with_no_allele`,
    /// so the two paths agree on the record Ensembl VEP annotates against nothing.
    #[test]
    fn the_library_path_sets_nearest_on_a_record_with_no_allele() {
        use crate::runner::make_runner_test_transcript;
        use crate::transcript_index::TranscriptIndex;
        use std::collections::HashMap;

        let mut coding = make_runner_test_transcript("ENSTCODING", "ENSGCODING", "21", 10_000);
        coding.biotype = "protein_coding".into();
        let mut by_chr = HashMap::new();
        by_chr.insert(
            "21".to_string(),
            TranscriptIndex::new(vec![coding], 5_000, 5_000),
        );
        let transcripts = Arc::new(LazyTranscriptIndexes::from_indexes(by_chr));
        let annotator = annotator_with_nearest(transcripts, Some(NearestType::Gene)).unwrap();

        let mut batch = parse_vcf_line("21\t20000\tr_star\tA\t*\t.\t.\t.", true).unwrap();
        batch.extend(parse_vcf_line("21\t20000\trs1\tA\tG\t.\t.\t.", true).unwrap());
        let _ = annotator.annotate_batch(&mut batch);

        assert!(batch[0].no_allele);
        assert!(batch[0].transcript_consequences.is_empty());
        assert!(batch[0].most_severe_consequence.is_none());
        for variant in &batch {
            assert_eq!(
                variant.nearest,
                Some(vec!["ENSGCODING".to_string()]),
                "{:?}",
                variant.id
            );
        }
    }

    /// The eager loader keeps the transcript set's members alone, the same
    /// transcripts the lazy path keeps chromosome by chromosome
    /// (AnnotationType/Transcript.pm:223-230): on the GENCODE corpus the primary
    /// set holds fewer transcripts than the basic set, which holds fewer than
    /// the cache.
    #[test]
    fn eager_indexes_keep_the_transcript_set_alone() {
        const CHROMOSOMES: [&str; 2] = ["21", "22"];
        let cache = format!(
            "{}/../../tests/golden/116/GRCh38-gencode-basic/json_cache",
            env!("CARGO_MANIFEST_DIR")
        );
        let effects = vep_effects::EffectsConfig {
            upstream_distance: 5_000,
            downstream_distance: 5_000,
            ..Default::default()
        };
        let indexed = |indexes: &LazyTranscriptIndexes| {
            CHROMOSOMES
                .iter()
                .map(|chr| indexes.get(chr).map_or(0, |index| index.len()))
                .sum::<usize>()
        };
        let count = |set: TranscriptSet| {
            let eager = load_transcript_indexes(&cache, &effects, set).unwrap();
            let lazy = LazyTranscriptIndexes::new(&cache, 5_000, 5_000)
                .unwrap()
                .with_transcript_set(set);
            for chr in CHROMOSOMES {
                lazy.prewarm(chr).unwrap();
            }
            let eager_count = indexed(&eager);
            assert_eq!(eager_count, indexed(&lazy), "{set:?}");
            eager_count
        };
        let every = count(TranscriptSet::default());
        let basic = count(TranscriptSet {
            gencode_basic: true,
            gencode_primary: false,
        });
        let primary = count(TranscriptSet {
            gencode_basic: false,
            gencode_primary: true,
        });
        assert!(
            primary > 0 && primary < basic && basic < every,
            "{primary} {basic} {every}"
        );
    }

    /// `from_config` resolves the cache as the `vep` binary does: a root in Ensembl
    /// VEP's layout and a `--json_cache` leaf both load the leaf's transcripts, and
    /// a root holding no cache is an error, not an empty index.
    #[test]
    fn from_config_resolves_the_cache_as_the_binary_does() {
        let golden = format!("{}/../../tests/golden/116", env!("CARGO_MANIFEST_DIR"));
        let root = format!("{golden}/GRCh38-run-control/cache_root");
        let leaf = format!("{golden}/GRCh38-release116/json_cache");
        let config = |flags: &[&str]| {
            let mut argv = vec!["vep", "--fork", "1"];
            argv.extend_from_slice(flags);
            Config::from_argv(argv).unwrap()
        };
        for flags in [&["--dir_cache", &root][..], &["--json_cache", &leaf][..]] {
            let annotator = Annotator::from_config(&config(flags)).unwrap();
            assert!(
                annotator.resources().transcripts.contains_key("21"),
                "{flags:?}"
            );
        }
        let empty = tempfile::tempdir().unwrap();
        let err = Annotator::from_config(&config(&["--dir_cache", empty.path().to_str().unwrap()]))
            .err()
            .expect("an empty cache root is refused");
        assert!(err.to_string().starts_with("ERROR: [cache] "), "{err}");
    }

    /// The library path resolves an input chromosome through the cache's
    /// synonyms as the binary does: Ensembl VEP fetches a record's regions
    /// through `get_source_chr_name` (`AnnotationSource.pm:243-245`, then
    /// `:179`; `BaseVEP.pm:689-738`), so on the GRCh38-synonyms corpus, whose
    /// `chr_synonyms.txt` lists `NC_000021.9` and `CM000683.2` for `21`, a row
    /// on either name gets the rows the `21` row gets and keeps its own name
    /// for the output.
    #[test]
    fn from_config_resolves_input_chromosomes_through_the_synonyms() {
        let cache = format!(
            "{}/../../tests/golden/116/GRCh38-synonyms/json_cache",
            env!("CARGO_MANIFEST_DIR")
        );
        let annotator = Annotator::from_config(
            &Config::from_argv(["vep", "--fork", "1", "--json_cache", &cache]).unwrap(),
        )
        .unwrap();
        let row = |chrom: &str| RawRow {
            chrom: chrom.to_string(),
            pos: 9_548_815,
            ref_allele: "T".to_string(),
            alt_allele: "C".to_string(),
        };
        let names = ["21", "NC_000021.9", "CM000683.2", "chr21"];
        let annotated: Vec<InputVariant> = annotator
            .try_annotate_rows(&names.map(row))
            .into_iter()
            .map(|result| result.unwrap())
            .collect();
        let rows = |variant: &InputVariant| -> Vec<(String, Vec<Consequence>)> {
            variant
                .transcript_consequences
                .iter()
                .map(|tc| (tc.transcript_id.to_string(), tc.consequences.to_vec()))
                .collect()
        };
        let expected = rows(&annotated[0]);
        assert!(expected.len() > 1, "{expected:?}");
        for (name, variant) in names.iter().zip(&annotated) {
            assert_eq!(rows(variant), expected, "{name}");
            assert_eq!(variant.chr, "21", "{name}");
            assert_eq!(variant.original_chr, *name);
        }
    }

    /// The library path resolves a paired breakend's mate chromosome through the
    /// cache's synonyms as it does the record's: Ensembl VEP fetches each
    /// breakend's region through `get_source_chr_name` (`AnnotationSource.pm:251-254`,
    /// then `:179`; `BaseVEP.pm:689-738`), so a record written `chr21` with a
    /// `]chr21:...]` mate gets the rows the bare `21` spelling gets, under the
    /// cache's name for both. On the GRCh38-synonyms corpus, whose `chr_synonyms.txt`
    /// lists `chr21` for `21`, Ensembl VEP 116.2 prints 60 rows for this record
    /// under either spelling, 30 transcripts for each of its two breakend
    /// alleles, as the `vep` binary does over the corpus's `json_cache`.
    #[test]
    fn from_config_resolves_a_breakend_mate_through_the_synonyms() {
        let cache = format!(
            "{}/../../tests/golden/116/GRCh38-synonyms/json_cache",
            env!("CARGO_MANIFEST_DIR")
        );
        let annotator = Annotator::from_config(
            &Config::from_argv(["vep", "--fork", "1", "--json_cache", &cache]).unwrap(),
        )
        .unwrap();
        let annotated = |chrom: &str| -> Vec<InputVariant> {
            let line = format!(
                "{chrom}\t31659665\tbnd\tA\t]{chrom}:31659730]A\t.\tPASS\tSVTYPE=BND;MATEID=mate"
            );
            let mut batch = parse_vcf_line(&line, false).unwrap();
            let _ = annotator.annotate_batch(&mut batch);
            batch
        };
        let rows = |variant: &InputVariant| -> Vec<(String, Vec<Consequence>)> {
            variant
                .transcript_consequences
                .iter()
                .map(|tc| (tc.transcript_id.to_string(), tc.consequences.to_vec()))
                .collect()
        };
        let bare = annotated("21");
        let prefixed = annotated("chr21");
        assert_eq!(bare.len(), 2, "the record's own breakend and its mate");
        assert_eq!(prefixed.len(), 2);
        let total = |batch: &[InputVariant]| -> usize {
            batch.iter().map(|variant| rows(variant).len()).sum()
        };
        assert_eq!(total(&bare), 60);
        assert_eq!(total(&prefixed), total(&bare));
        for (bare_variant, prefixed_variant) in bare.iter().zip(&prefixed) {
            let bare_rows = rows(bare_variant);
            assert_eq!(
                rows(prefixed_variant),
                bare_rows,
                "{}",
                prefixed_variant.display_allele()
            );
            assert!(!bare_rows.is_empty(), "{}", bare_variant.display_allele());
            assert_eq!(prefixed_variant.chr, "21");
            assert_eq!(prefixed_variant.original_chr, "chr21");
            assert_eq!(prefixed_variant.mate_chr, bare_variant.mate_chr);
        }
        let paired = prefixed.iter().find(|v| !v.is_single_breakend).unwrap();
        assert_eq!(paired.mate_chr.as_deref(), Some("21"));
    }

    /// `from_config` carries the co-location flags into the library path: on the
    /// GRCh38-colocated corpus, `21 5230609 A C` gains the failed null allele
    /// `rs1979841086` beside `rs373402393` under `--failed 1`, as Ensembl VEP 116.2
    /// prints it, and `rs373402393` alone without the flag.
    #[test]
    fn from_config_carries_the_colocation_rules() {
        let cache = format!(
            "{}/../../tests/golden/116/GRCh38-colocated/json_cache",
            env!("CARGO_MANIFEST_DIR")
        );
        let existing = |flags: &[&str]| -> Vec<String> {
            let mut argv = vec![
                "vep",
                "--fork",
                "1",
                "--json_cache",
                &cache,
                "--check_existing",
            ];
            argv.extend_from_slice(flags);
            let annotator = Annotator::from_config(&Config::from_argv(argv).unwrap()).unwrap();
            let mut batch = parse_vcf_line("21\t5230609\t.\tA\tC\t.\t.\t.", false).unwrap();
            let _ = annotator.annotate_batch(&mut batch);
            batch[0].existing_variation.clone()
        };
        assert_eq!(existing(&[]), ["rs373402393"]);
        assert_eq!(
            existing(&["--failed", "1"]),
            ["rs1979841086", "rs373402393"]
        );
    }

    /// The HGVSg options follow the flags: the shift is off under `--shift_hgvs 0`,
    /// and `--hgvsg_use_accession` maps every name of the synonyms table that has
    /// an `NC_` synonym to the first one (`21` and `MT` here, `chr21` to nothing,
    /// as its one synonym is `21`). The SPDI options carry the same map for the
    /// VRS object under `--ga4gh_vrs` and none under `--spdi` alone.
    #[test]
    fn hgvsg_and_spdi_options_take_the_accessions_from_the_synonyms() {
        let table = ChromosomeSynonyms::parse(
            "21\tNC_000021.9\n21\tchr21\n21\tCM000683.2\nMT\tchrM\nMT\tNC_012920.1\n",
        );
        let accessions = HashMap::from([
            ("21".to_string(), "NC_000021.9".to_string()),
            ("MT".to_string(), "NC_012920.1".to_string()),
        ]);
        let config = Config::from_argv(["vep", "--hgvsg", "--hgvsg_use_accession"]).unwrap();
        assert_eq!(
            hgvsg_options(&config, &table),
            Some(vep_effects::hgvsg::HgvsgOptions {
                shift: true,
                accessions: Some(accessions.clone()),
            })
        );
        let config = Config::from_argv(["vep", "--hgvsg", "--shift_hgvs", "0"]).unwrap();
        assert_eq!(
            hgvsg_options(&config, &table),
            Some(vep_effects::hgvsg::HgvsgOptions {
                shift: false,
                accessions: None,
            })
        );
        let config = Config::from_argv(["vep", "--hgvsg_use_accession"]).unwrap();
        assert_eq!(hgvsg_options(&config, &table), None);

        let config = Config::from_argv(["vep", "--ga4gh_vrs", "--tab"]).unwrap();
        assert_eq!(
            spdi_options(&config, &table),
            Some(vep_effects::spdi::SpdiOptions {
                vrs_accessions: Some(accessions),
            })
        );
        let config = Config::from_argv(["vep", "--spdi"]).unwrap();
        assert_eq!(
            spdi_options(&config, &table),
            Some(vep_effects::spdi::SpdiOptions {
                vrs_accessions: None,
            })
        );
        assert_eq!(
            spdi_options(&Config::from_argv(["vep"]).unwrap(), &table),
            None
        );
        let empty = ChromosomeSynonyms::new();
        let config = Config::from_argv(["vep", "--ga4gh_vrs", "--tab"]).unwrap();
        assert_eq!(
            spdi_options(&config, &empty).unwrap().vrs_accessions,
            Some(HashMap::new())
        );
    }

    /// `from_config` refuses `--hgvsg` with no FASTA as the binary does, with
    /// Ensembl VEP's message (Runner.pm 759) under the `fasta` component.
    #[test]
    fn from_config_refuses_hgvsg_without_a_fasta() {
        let leaf = format!(
            "{}/../../tests/golden/116/GRCh38-release116/json_cache",
            env!("CARGO_MANIFEST_DIR")
        );
        let config =
            Config::from_argv(["vep", "--fork", "1", "--json_cache", &leaf, "--hgvsg"]).unwrap();
        let err = Annotator::from_config(&config)
            .err()
            .expect("--hgvsg without a FASTA is refused");
        assert_eq!(
            CliError::from_run(&err),
            format!("ERROR: [fasta] {HGVS_NEEDS_FASTA}")
        );
    }

    /// `from_config` compiles the `--transcript_filter` before it opens the
    /// FASTA and runs its checks, as the binary does (`Runner::run`): `--hgvs`
    /// with no FASTA and a filter Ensembl VEP would refuse reports the `filter`
    /// component, not `fasta`.
    #[test]
    fn from_config_refuses_a_bad_transcript_filter_before_a_missing_fasta() {
        let leaf = format!(
            "{}/../../tests/golden/116/GRCh38-release116/json_cache",
            env!("CARGO_MANIFEST_DIR")
        );
        let config = Config::from_argv([
            "vep",
            "--fork",
            "1",
            "--json_cache",
            &leaf,
            "--hgvs",
            "--transcript_filter",
            "biotype equals x",
        ])
        .unwrap();
        let err = Annotator::from_config(&config)
            .err()
            .expect("a filter Ensembl VEP would refuse is refused");
        assert_eq!(
            CliError::from_run(&err),
            "ERROR: [filter] No such operator \"equals\""
        );
    }

    /// `from_config` compiles the `--transcript_filter` and the batch path applies
    /// it as the binary does: over the transcript-filter corpus under
    /// `biotype is protein_coding`, the embedding caller's rows, by record, allele
    /// and feature, are exactly the rows of Ensembl VEP 116.2's reference output
    /// for that flag set, intergenic rows for the records every transcript of
    /// which fails included, and fewer than the unfiltered rows.
    #[test]
    fn from_config_applies_the_transcript_filter_as_the_binary_does() {
        use std::io::Read as _;
        let corpus = format!(
            "{}/../../tests/golden/116/GRCh38-transcript-filter",
            env!("CARGO_MANIFEST_DIR")
        );
        let cache = format!("{corpus}/json_cache");
        let rows = |filter: &[&str]| -> Vec<(String, String, String)> {
            let mut argv = vec!["vep", "--fork", "1", "--json_cache", &cache];
            argv.extend_from_slice(filter);
            let annotator = Annotator::from_config(&Config::from_argv(argv).unwrap()).unwrap();
            let vcf = std::fs::read_to_string(format!("{corpus}/variants.vcf")).unwrap();
            let mut batch: Vec<InputVariant> = vcf
                .lines()
                .filter(|l| !l.starts_with('#'))
                .flat_map(|l| parse_vcf_line(l, false).unwrap())
                .collect();
            let _ = annotator.annotate_batch(&mut batch);
            let mut out = Vec::new();
            for v in &batch {
                for tc in &v.transcript_consequences {
                    out.push((
                        v.uploaded_variation(),
                        v.display_allele(),
                        tc.transcript_id.to_string(),
                    ));
                }
                if v.transcript_consequences.is_empty() && v.most_severe_consequence.is_some() {
                    out.push((v.uploaded_variation(), v.display_allele(), "-".to_string()));
                }
            }
            out.sort();
            out
        };
        let filtered = rows(&["--transcript_filter", "biotype is protein_coding"]);
        let unfiltered = rows(&[]);
        assert!(filtered.len() < unfiltered.len());

        let mut reference = String::new();
        flate2::read::GzDecoder::new(
            std::fs::File::open(format!("{corpus}/expected_filter-biotype/default.txt.gz"))
                .unwrap(),
        )
        .read_to_string(&mut reference)
        .unwrap();
        let mut expected: Vec<(String, String, String)> = reference
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f[0].to_string(), f[2].to_string(), f[4].to_string())
            })
            .collect();
        expected.sort();
        assert_eq!(filtered, expected);
    }

    /// `from_config` lands every HGVS-related field of the effects configuration
    /// as the binary's run computes it from the same configuration, on the
    /// GRCh38-identifiers corpus with its reference: `--hgvsg` gives the shifted
    /// HGVSg options, `--shift_hgvs 0` the unshifted ones with the HGVSc and
    /// HGVSp options beside them, and `--hgvs` alone computes HGVSc and HGVSp
    /// with no HGVSg.
    #[test]
    fn from_config_carries_the_hgvs_and_hgvsg_options() {
        let corpus = format!(
            "{}/../../tests/golden/116/GRCh38-identifiers",
            env!("CARGO_MANIFEST_DIR")
        );
        let cache = format!("{corpus}/json_cache");
        let fasta = format!("{corpus}/reference.fa.gz");
        let effects = |flags: &[&str]| {
            let mut argv = vec![
                "vep",
                "--fork",
                "1",
                "--json_cache",
                &cache,
                "--fasta",
                &fasta,
            ];
            argv.extend_from_slice(flags);
            let config = Config::from_argv(argv).unwrap();
            let annotator = Annotator::from_config(&config).unwrap();
            let effects = Arc::clone(&annotator.resources().effects_config);
            assert_eq!(effects.hgvs, hgvs_options(&config), "{flags:?}");
            // No flag here asks for an accession, so the table does not enter
            // the options.
            assert_eq!(
                effects.hgvsg,
                hgvsg_options(&config, &ChromosomeSynonyms::new()),
                "{flags:?}"
            );
            assert_eq!(
                effects.compute_hgvs,
                crate::runner::fasta_setup_checks(&config, true)
                    .unwrap()
                    .hgvs,
                "{flags:?}"
            );
            assert_eq!(
                effects.compute_exon_intron_numbers,
                crate::runner::needs_exon_intron_numbers(&config),
                "{flags:?}"
            );
            assert_eq!(effects.enable_indel_3prime_shift, config.shift_3prime);
            assert_eq!(effects.protein_version, config.protein_version);
            assert_eq!(
                (effects.upstream_distance, effects.downstream_distance),
                config.distance
            );
            assert!(effects.reference_fasta.is_some(), "{flags:?}");
            effects
        };

        let shifted = effects(&["--hgvsg"]);
        assert_eq!(shifted.hgvsg.as_ref().map(|o| o.shift), Some(true));
        assert!(shifted.hgvs.shift_hgvs && !shifted.compute_hgvs);

        let unshifted = effects(&["--hgvsg", "--shift_hgvs", "0"]);
        assert_eq!(unshifted.hgvsg.as_ref().map(|o| o.shift), Some(false));
        assert!(!unshifted.hgvs.shift_hgvs && !unshifted.enable_indel_3prime_shift);

        let hgvs = effects(&[
            "--hgvs",
            "--hgvsp_use_prediction",
            "--remove_hgvsp_version",
            "--protein_version",
            "--numbers",
            "--distance",
            "2000,1000",
        ]);
        assert!(hgvs.compute_hgvs && hgvs.hgvsg.is_none());
        assert!(hgvs.hgvs.hgvsp_use_prediction && hgvs.hgvs.remove_hgvsp_version);
        assert!(hgvs.protein_version && hgvs.compute_exon_intron_numbers);
        assert_eq!(
            (hgvs.upstream_distance, hgvs.downstream_distance),
            (2000, 1000)
        );
    }

    /// `--shift_hgvs 0` is the one setting that turns every HGVS shift off
    /// (`shift_hgvs_variants_3prime`, read by `hgvs_genomic` at VariationFeature.pm
    /// 1981-1988 as by `hgvs_transcript` and `hgvs_protein`), so the HGVSg options
    /// and the HGVSc and HGVSp options read the same flag; the two HGVSp forms
    /// travel beside it.
    #[test]
    fn hgvs_and_hgvsg_options_follow_the_flags() {
        let config = |flags: &[&str]| {
            let mut argv = vec!["vep", "--fork", "1", "--hgvsg"];
            argv.extend_from_slice(flags);
            Config::from_argv(argv).unwrap()
        };
        let no_table = ChromosomeSynonyms::new();
        let default = config(&[]);
        assert_eq!(
            hgvsg_options(&default, &no_table).map(|o| o.shift),
            Some(true)
        );
        assert_eq!(
            hgvs_options(&default),
            vep_effects::hgvs::HgvsOptions::default()
        );

        let unshifted = config(&["--shift_hgvs", "0"]);
        assert_eq!(
            hgvsg_options(&unshifted, &no_table).map(|o| o.shift),
            Some(false)
        );
        assert!(!hgvs_options(&unshifted).shift_hgvs);
        assert!(hgvs_options(&config(&["--shift_3prime", "1"])).shift_3prime);

        let forms = config(&["--hgvsp_use_prediction", "--remove_hgvsp_version"]);
        assert_eq!(
            hgvs_options(&forms),
            vep_effects::hgvs::HgvsOptions {
                shift_hgvs: true,
                hgvsp_use_prediction: true,
                remove_hgvsp_version: true,
                shift_3prime: false,
            }
        );
    }

    /// An embedding caller passes the whole `HgvsOptions`, so the two HGVSp
    /// forms reach it as `--shift_hgvs` does; HGVSg stays off on this path.
    #[test]
    fn build_effects_config_carries_the_hgvs_options() {
        let options = vep_effects::hgvs::HgvsOptions {
            shift_hgvs: false,
            hgvsp_use_prediction: true,
            remove_hgvsp_version: true,
            shift_3prime: false,
        };
        let effects =
            build_effects_config((5_000, 5_000), None, false, options, true, false, false);
        assert_eq!(effects.hgvs, options);
        assert!(effects.compute_hgvs && effects.hgvsg.is_none());
    }

    /// The annotator path reaches the same record-level breakend verdict as the CLI:
    /// a breakend record with no transcript within range of either breakend is one
    /// `intergenic_variant` per allele, and a record whose mate reaches a transcript
    /// gets no intergenic row.
    #[test]
    fn annotator_applies_the_breakend_record_verdict() {
        let dir = tmp_json_cache("bnd_record");
        let lazy = Arc::new(
            LazyTranscriptIndexes::new(dir.path().to_str().unwrap(), 5_000, 5_000).unwrap(),
        );
        lazy.prewarm("21").unwrap();
        let annotator = annotator_over(lazy).unwrap();

        // The one transcript of `tmp_json_cache` spans 21:100-200.
        let mut intergenic = parse_vcf_line(
            "21\t1000000\tbnd_a\tN\tN[21:1000500[\t.\t.\tSVTYPE=BND;MATEID=bnd_b",
            true,
        )
        .unwrap();
        let _ = annotator.annotate_batch(&mut intergenic);
        assert_eq!(intergenic.len(), 2);
        for allele in &intergenic {
            assert!(allele.transcript_consequences.is_empty());
            assert_eq!(
                allele.most_severe_consequence,
                Some(Consequence::IntergenicVariant)
            );
        }

        let mut genic_mate = parse_vcf_line(
            "21\t1000000\tbnd_a\tN\tN[21:150[\t.\t.\tSVTYPE=BND;MATEID=bnd_b",
            true,
        )
        .unwrap();
        let _ = annotator.annotate_batch(&mut genic_mate);
        assert!(genic_mate[0].transcript_consequences.is_empty());
        assert_eq!(genic_mate[0].most_severe_consequence, None);
        assert_eq!(genic_mate[1].transcript_consequences.len(), 1);
    }

    /// A plugins directory that does not exist is fatal under `--safe` when a
    /// plugin is named (Runner.pm:911-916), the plugin component naming the
    /// directory with Ensembl VEP's text; with no plugin named it is not looked at.
    #[test]
    fn missing_plugins_dir_dies_under_safe() {
        let err = load_plugins(
            &["NoSuchPlugin".to_string()],
            "/no/such/dir",
            true,
            None,
            None,
        )
        .err()
        .expect("a missing plugins directory is fatal under --safe");
        assert_eq!(
            err.downcast_ref::<crate::config::CliError>()
                .unwrap()
                .to_string(),
            "ERROR: [plugin] Plugins directory '/no/such/dir' not found."
        );
        let loaded = load_plugins(&[], "/no/such/dir", true, None, None).unwrap();
        assert_eq!(loaded.builtin.plugin_count(), 0);
    }

    /// `Name,param,param` splits on every comma (Runner.pm:928: `split /,/`), the
    /// first part the plugin, the rest its parameters in order: LoFtool reads its
    /// scores file from the first parameter, so the file it fails to open names
    /// the parameter the split gave it. The failure is fatal under `--safe` alone
    /// (`:946-954`: `throw($msg) if $self->param('safe')`, else a warning and
    /// `next`), the plugin component carrying Ensembl VEP's text; without it the
    /// run goes on with no plugin loaded.
    #[test]
    fn plugin_name_and_params_split_as_perl_does() {
        let dir = tempfile::tempdir().unwrap();
        let dir_path = dir.path().to_str().unwrap();
        for (arg, first_param) in [
            ("LoFtool,/no/such/scores.txt", "/no/such/scores.txt"),
            ("LoFtool,/no/such/first,second,third", "/no/such/first"),
        ] {
            let err = load_plugins(&[arg.to_string()], dir_path, true, None, None)
                .err()
                .expect("a scores file that does not exist fails init under --safe");
            let message = CliError::from_run(&err);
            assert!(
                message.starts_with(&format!(
                    "ERROR: [plugin] Failed to instantiate plugin LoFtool: LoFtool: failed to open {first_param}:"
                )),
                "{arg}: {message}"
            );
            let loaded = load_plugins(&[arg.to_string()], dir_path, false, None, None).unwrap();
            assert_eq!(loaded.builtin.plugin_count(), 0);
        }
        let err = load_plugins(&["LoFtool".to_string()], dir_path, true, None, None)
            .err()
            .expect("LoFtool without a parameter fails init under --safe");
        assert_eq!(
            CliError::from_run(&err),
            "ERROR: [plugin] Failed to instantiate plugin LoFtool: LoFtool requires a path to the scores file"
        );
    }
}
