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

use crate::config::{CliError, Config, PluginDirs};
use crate::dylib_plugins::DylibOutcome;
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
use vep_builtins::custom::CustomSources;
use vep_builtins::traits::{OutputFormat, RunContext};
use vep_builtins::BuiltinRegistry;
use vep_core::consequence::Consequence;
use vep_core::synonyms::ChromosomeSynonyms;
use vep_core::variant::{InputVariant, VariantClass};
use vep_effects::LofteeContextRows;
use vep_plugin::PluginRegistry;

/// The refusal of `--hgvs` or `--hgvsg` with no reference FASTA, Ensembl VEP's
/// text (Runner.pm 759), under the `fasta` component.
pub(crate) const HGVS_NEEDS_FASTA: &str = "Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode without a FASTA file (see --fasta)";

/// Shared annotation resources that can be reused across batches.
pub struct AnnotationResources {
    pub transcripts: Arc<LazyTranscriptIndexes>,
    pub variation_source: Option<Arc<VariationSource>>,
    pub regulatory_source: Option<Arc<crate::regulatory::RegulatorySource>>,
    /// The epigenomes `--cell_type` named, resolved against the regulatory
    /// source's list; `None` without the flag or without a source.
    pub cell_types: Option<Arc<crate::regulatory::CellTypes>>,
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
    /// The `--custom` sources, which annotate each batch after the plugins;
    /// `None` without a specification.
    pub custom_sources: Option<Arc<CustomSources>>,
    /// What the output headers print for each custom source, in order
    /// ([`Annotator::custom_source_infos`]); empty without a source.
    pub custom_source_infos: Vec<vep_io::output::fields::CustomSourceInfo>,
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
    /// [`Annotator::bam_edited`].
    bam_edited: bool,
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
    ///
    /// The regulatory source, when the configuration has one, is loaded whole
    /// here for the same reason: a corrupt regulatory shard fails construction,
    /// as a corrupt transcript shard does, instead of leaving every variant of
    /// its chromosome without regulatory rows.
    ///
    /// With a reference sequence in `resources.effects_config`, a stop codon
    /// an `_rna_edit` completes is paired with the bases the genome holds of
    /// it, as `load_transcript_indexes` pairs it: an index the caller built
    /// without the reference sequence is replaced, in the annotator's
    /// resources, by a copy holding the paired transcripts, and one a load
    /// already paired is kept as given.
    pub fn new(resources: AnnotationResources, pool: Arc<ThreadPool>) -> anyhow::Result<Self> {
        Self::with_bam_edited(resources, pool, false)
    }

    /// [`Self::new`] with the value [`Self::bam_edited`] returns.
    fn with_bam_edited(
        mut resources: AnnotationResources,
        pool: Arc<ThreadPool>,
        bam_edited: bool,
    ) -> anyhow::Result<Self> {
        // Verified at construction; the wrapper is dropped because the batch path
        // needs the `Arc` itself, and the check cannot be invalidated afterwards:
        // materialization only ever moves from empty to filled.
        MaterializedTranscriptIndexes::try_from_lazy(Arc::clone(&resources.transcripts))?;
        if let Some(reference) = &resources.effects_config.reference_fasta {
            if let Some(paired) = resources.transcripts.paired_with(reference) {
                resources.transcripts = Arc::new(paired);
            }
        }
        if let Some(source) = &resources.regulatory_source {
            source.prewarm_all()?;
        }
        Ok(Self {
            inner: Arc::new(AnnotatorInner {
                resources,
                pool,
                bam_edited,
            }),
        })
    }

    /// The pool is `--fork`'s count by the host's rule alone: an embedding
    /// consumer's variants arrive without an input file whose shape could lower
    /// it.
    pub fn from_config(config: &Config) -> anyhow::Result<Self> {
        let pool = build_thread_pool(config.fork.host_count())?;
        Self::from_config_with_pool(config, pool)
    }

    /// Resolves the cache and the FASTA the configuration names and settles the
    /// HGVS flags, as the `vep` binary does, so a configuration with no cache is
    /// an error here too rather than an empty index, `--fasta_dir` is read, and
    /// `--hgvs` without a FASTA is refused.
    pub fn from_config_with_pool(config: &Config, pool: Arc<ThreadPool>) -> anyhow::Result<Self> {
        // The configuration's notice comes first, as the binary prints it
        // (`Runner::new`, Runner.pm:117-120), before a missing cache is refused
        // (CacheDir.pm:109).
        log_notices(&config.notices);
        // The resolution settles the cache's decisions on a copy of the
        // configuration, the one read from here on.
        let mut settled = config.clone();
        let cache = crate::cache_resolve::resolve(&mut settled)?;
        let config = &settled;
        // After the cache and before the cache's notices and the FASTA, as the
        // binary orders them (`Runner::new` resolves the cache; `Runner::run`
        // compiles the filter, prints the cache's notices, then opens the FASTA
        // and runs its checks): with a missing cache, a bad filter and a missing
        // FASTA together, both paths report them in that order.
        let transcript_filter =
            crate::runner::compile_transcript_filter(&config.transcript_filter)?;
        log_notices(&cache.notices);
        let predictions = crate::runner::prediction_checks(config, &cache.info)?;
        log_notices(&predictions.notices);
        let custom_specs = crate::runner::parse_custom_specs(config)?;
        let reference_fasta = crate::runner::open_fasta(config, &cache)?;
        let custom_sources = crate::runner::open_custom_sources(
            custom_specs,
            reference_fasta.as_deref(),
            &cache.synonyms,
            crate::runner::custom_run_options(config),
        )?;
        let checks = match crate::runner::fasta_setup_checks(config, reference_fasta.is_some()) {
            Ok(checks) => checks,
            Err(refusal) => {
                log_notices(&refusal.notices);
                return Err(refusal.error.into());
            }
        };
        log_notices(&checks.notices);

        let LoadedPlugins { builtin, dylib } = load_plugins(
            &config.plugin,
            &config.dir_plugins,
            config.safe,
            reference_fasta.clone(),
            cache.assembly.as_deref(),
            run_context(config, checks.hgvs),
        )?;
        let mut effects = new_effects_config(
            config.distance,
            reference_fasta,
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
        effects.loftee_context = loftee_context_rows(&builtin);
        effects.context_sequences = context_sequence_rows(&builtin);
        let effects_config = Arc::new(effects);
        let builtin_plugins = Arc::new(builtin);
        let dylib_plugins = Arc::new(dylib);
        let regulatory = regulatory_enabled(config, &dylib_plugins);
        let transcripts = load_transcript_indexes_with_exon_ids(
            &cache.dir,
            &effects_config,
            TranscriptSet::from_config(config),
            exon_ids_wanted(&config.plugin),
        )?;
        let variation_source = crate::runner::open_variation_source(config, &cache)?;
        let sv_source = crate::runner::open_sv_source(config)?;
        let regulatory_source =
            crate::regulatory::RegulatorySource::open_for(config, &cache, regulatory)?;
        let cell_types = crate::regulatory::CellTypes::resolve(
            config.cell_type.as_deref(),
            regulatory_source.as_deref(),
        )?;

        Self::with_bam_edited(
            AnnotationResources {
                transcripts,
                variation_source,
                regulatory_source,
                cell_types,
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
                custom_source_infos: crate::runner::custom_source_infos(
                    custom_sources.as_deref(),
                    crate::runner::custom_run_options(config),
                ),
                custom_sources,
            },
            pool,
            config.bam_edited,
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

    /// The names of every loaded plugin's fields: the built-in plugins' first,
    /// then the dylibs', each group in load order with every plugin's keys sorted.
    pub fn plugin_header_fields(&self) -> Vec<String> {
        let resources = &self.inner.resources;
        resources
            .builtin_plugins
            .all_header_info()
            .into_iter()
            .chain(resources.dylib_plugins.all_header_info())
            .map(|(name, _)| name)
            .collect()
    }

    /// What the output headers print for each `--custom` source, in order: an
    /// embedder puts them in its `FieldOptions::custom_sources` so the formats
    /// list the sources' keys after the plugin fields. Empty without a source.
    pub fn custom_source_infos(&self) -> Vec<vep_io::output::fields::CustomSourceInfo> {
        self.inner.resources.custom_source_infos.clone()
    }

    /// Whether each row is analysed against its transcript's own reference
    /// allele: `--use_transcript_ref`, which a BAM-edited cache also switches
    /// on, unless `--use_given_ref` was given, as [`Self::from_config`]
    /// resolves it, on the annotator's copy of the configuration and not on
    /// the caller's. An embedder puts it in its
    /// `FieldOptions::use_transcript_ref` so the formats write the `GIVEN_REF`
    /// and `USED_REF` fields as the `vep` binary does.
    pub fn use_transcript_ref(&self) -> bool {
        self.inner.resources.effects_config.use_transcript_ref
    }

    /// Whether [`Self::from_config`] resolved a BAM-edited cache without
    /// `--use_given_ref`, the case in which the `vep` binary lists the
    /// `BAM_EDIT` field: an embedder puts it in its `FieldOptions::bam_edited`.
    /// `false` for an annotator built by [`Self::new`], which resolves no cache.
    pub fn bam_edited(&self) -> bool {
        self.inner.bam_edited
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
    /// occurred, with the error component it is reported under. A built-in
    /// plugin's failures on the batch are written as the batch ends
    /// ([`crate::warnings::Sink::write_plugin_failures`]), one line per plugin:
    /// a program has no run end for the sink to hold them to.
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
        // materialized; `Annotator::new` enforces that, and loads the regulatory
        // source whole, hence no per-batch check of either.
        let stats = pool.install(|| annotate_batch(batch, resources, use_parallel, &failure));
        crate::warnings::global().write_plugin_failures();
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

/// The assembly hint the built-in plugins read, from the run's assembly
/// (`--assembly`, else the cache's own): one of the two assemblies
/// [`Assembly`](vep_core::assembly::Assembly) names, or `None` for no assembly
/// and for a name [`Assembly::parse`](vep_core::assembly::Assembly::parse) reads
/// as neither (`GRCm39`). `load_plugins` warns about an assembly that gives no
/// hint, once per built-in it loads, so a run without plugins says nothing
/// about it.
pub fn resolve_assembly(raw: Option<&str>) -> Option<vep_core::assembly::Assembly> {
    raw.and_then(vep_core::assembly::Assembly::parse)
}

/// The plugins one run loads: the built-ins claimed by name, and the C-ABI dylibs
/// the remaining `--plugin` arguments resolved to.
pub(crate) struct LoadedPlugins {
    pub(crate) builtin: BuiltinRegistry,
    pub(crate) dylib: PluginRegistry,
}

/// Whether the run annotates regulatory features: `--regulatory` (given or
/// implied by `--everything`), or a loaded plugin that declares regulatory or
/// motif features, which Ensembl VEP switches the option on for
/// (`Runner.pm:998-1002`).
pub(crate) fn regulatory_enabled(config: &Config, dylib: &PluginRegistry) -> bool {
    if config.regulatory {
        return true;
    }
    let wanting = dylib.plugins_wanting_regulatory_features();
    for plugin in &wanting {
        info!("Fetching regulatory features for plugin: {plugin}");
    }
    !wanting.is_empty()
}

/// Load the `--plugin` arguments as `get_all_Plugins` does (Runner.pm:901-960):
/// each match of the plugins directory value that is not a directory is warned
/// about when a plugin is named, an error under `safe` (`:911-918`); each argument
/// is `Name,param,param` (`:928`), loaded as a built-in when one claims the name,
/// otherwise as a dylib resolved by path or under the directories, the last match
/// first ([`crate::dylib_plugins::load_unresolved`]). A name that resolves to
/// neither is settled in this order: an Ensembl VEP plugin vep-rs does not
/// provide ([`vep_builtins::unported_reason`]) is the plugin component's error
/// naming the reason and the support table, with or without `safe`; a module
/// `<dir>/<Name>.pm` under the directories is refused as a Perl plugin, with or
/// without `safe`, as is a module given by path; anything else is warned about
/// and dropped, an error under `safe` (`:929-939`), as is a dylib that does not
/// open and a built-in or a dylib whose `init` fails (`Failed to instantiate
/// plugin <Name>: <reason>`, `:946-954`). The arguments load one at a time in
/// their order, so the warnings follow it as Ensembl VEP's single loop writes
/// them (`:925`).
///
/// The assembly hint ([`resolve_assembly`] of `assembly`, the run's assembly name)
/// is load-bearing rather than advisory: REVEL picks which of its two coordinate
/// columns tabix indexed from it, LoFTEE picks between a bigWig and a tabix-TSV
/// GERP reader, and CADD / dbscSNV use it to refuse a data file whose own declared
/// assembly contradicts the run. A hint of `None` leaves every one of those guards
/// inert, so a wrong-assembly file annotates nothing and still exits 0; an assembly
/// name that gives no hint is therefore warned about as each built-in loads under
/// it, naming the plugin, and not at all when no plugin loads. `run` carries the
/// output format and the transcript set Ensembl VEP's plugins read from
/// `$self->{config}` as they construct; the directories of `dir_plugins` reach
/// the built-ins through it, so LoFtool and pLI find their bundled score files.
pub(crate) fn load_plugins(
    plugin_args: &[String],
    dir_plugins: &PluginDirs,
    safe: bool,
    reference_fasta: Option<Arc<vep_fasta::IndexedFasta>>,
    assembly: Option<&str>,
    run: RunContext<'_>,
) -> anyhow::Result<LoadedPlugins> {
    let mut builtin = BuiltinRegistry::new();
    let mut dylib = PluginRegistry::new();
    if plugin_args.is_empty() {
        return Ok(LoadedPlugins { builtin, dylib });
    }
    let run = RunContext {
        plugin_dirs: &dir_plugins.dirs,
        ..run
    };
    let hint = resolve_assembly(assembly);
    let unrecognised_assembly = assembly.filter(|_| hint.is_none());

    for missing in &dir_plugins.not_found {
        plugin_warning(format!("Plugins directory '{missing}' not found."), safe)?;
    }

    for arg in plugin_args {
        let outcome = builtin.load_from_args_for_run(
            std::slice::from_ref(arg),
            reference_fasta.clone(),
            hint,
            run,
        );
        let builtin_loaded = outcome.failed.is_empty() && outcome.unresolved.is_empty();
        if let Some(raw) = unrecognised_assembly {
            if builtin_loaded {
                tracing::warn!(
                    "--assembly '{raw}' is not an assembly {} recognises; it runs without an assembly hint",
                    vep_core::sv_routing::split_plugin_arg(arg).0
                );
            }
        }
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
            let dylib_outcome =
                crate::dylib_plugins::load_unresolved(&mut dylib, unresolved, dir_plugins, safe)?;
            if dylib_outcome == DylibOutcome::NotFound {
                let name = vep_core::sv_routing::split_plugin_arg(unresolved).0;
                settle_unresolved_name(name, dir_plugins, safe)?;
            }
        }
    }

    Ok(LoadedPlugins { builtin, dylib })
}

/// The sentence every unported-plugin refusal ends with.
const SUPPORT_TABLE: &str =
    "The plugin support table in docs/plugins.md lists every Ensembl VEP plugin and its status.";

/// Settles a `--plugin` name that is neither a built-in nor a file: the plugin
/// component's error for an Ensembl VEP plugin vep-rs does not provide, or for
/// a Perl module of that name under the plugins directories, whatever `safe`
/// says; otherwise the `Failed to load plugin` warning, an error under `safe`
/// (Runner.pm:929-939).
fn settle_unresolved_name(name: &str, dir_plugins: &PluginDirs, safe: bool) -> anyhow::Result<()> {
    if let Some(reason) = vep_builtins::unported_reason(name) {
        return Err(CliError::new(
            "plugin",
            format!(
                "{name} is an Ensembl VEP plugin vep-rs does not provide: {reason} {SUPPORT_TABLE}"
            ),
        )
        .into());
    }
    if let Some(module) = crate::dylib_plugins::perl_module_under(dir_plugins, name) {
        return Err(crate::dylib_plugins::perl_plugin_refusal(&module).into());
    }
    let searched = dir_plugins.searched();
    plugin_warning(
        format!(
            "Failed to load plugin {name}: {name} is not a built-in plugin and no lib{name}.so or lib{name}.dylib is under {searched}"
        ),
        safe,
    )
}

/// The run's settings the built-in plugins read before `init`: the output format
/// as the writers name it, the `--refseq` / `--merged` transcript set and
/// whether the run computes the HGVS notations (`hgvs`: the `--hgvs` the FASTA
/// checks kept); [`load_plugins`] adds the plugins directories it is handed.
pub(crate) fn run_context(config: &Config, hgvs: bool) -> RunContext<'static> {
    RunContext {
        output_format: OutputFormat::from_name(&config.output_format),
        refseq: config.refseq,
        merged: config.merged,
        hgvs,
        plugin_dirs: &[],
    }
}

/// The built-in plugins that read the transcript model on every row.
const EVERY_ROW_CONTEXT_PLUGINS: [&str; 4] = ["LoFTEE", "MaxEntScan", "NMD", "NearestExonJB"];

/// Whether the cache loads keep each exon's stable id ([`vep_core::transcript::Transcript::exon_ids`]):
/// for a run whose `--plugin` arguments name a plugin that prints them
/// (NearestExonJB). The arguments are read, not the loaded plugins: the
/// binary builds the transcript indexes before the plugins load.
pub(crate) fn exon_ids_wanted(plugin_args: &[String]) -> bool {
    plugin_args
        .iter()
        .any(|arg| vep_core::sv_routing::split_plugin_arg(arg).0 == "NearestExonJB")
}

/// The consequence rows that carry the transcript model a loaded built-in reads:
/// every row for a plugin `EVERY_ROW_CONTEXT_PLUGINS` names, the
/// `5_prime_UTR_variant` rows for UTRAnnotator, none otherwise. The loaded
/// plugins' names are read, not the arguments': `--plugin LoF` loads LoFTEE too.
pub(crate) fn loftee_context_rows(builtin: &BuiltinRegistry) -> LofteeContextRows {
    let names = builtin.plugin_names();
    let is_loaded = |plugin: &str| names.iter().any(|name| name == plugin);
    if EVERY_ROW_CONTEXT_PLUGINS
        .iter()
        .any(|plugin| is_loaded(plugin))
    {
        LofteeContextRows::Every
    } else if is_loaded("UTRAnnotator") {
        LofteeContextRows::WithTerm(Consequence::FivePrimeUtrVariant)
    } else {
        LofteeContextRows::None
    }
}

/// The consequence rows that carry the transcript's coding, 3' UTR and protein
/// sequences: the `frameshift_variant` rows for Downstream, which translates
/// the changed coding sequence, none otherwise.
pub(crate) fn context_sequence_rows(builtin: &BuiltinRegistry) -> LofteeContextRows {
    if builtin
        .plugin_names()
        .iter()
        .any(|name| name == "Downstream")
    {
        LofteeContextRows::WithTerm(Consequence::FrameshiftVariant)
    } else {
        LofteeContextRows::None
    }
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

/// Runs the built-in plugins over `batch` and warns about each one that failed
/// ([`report_plugin_failures`]).
pub(crate) fn run_builtin_plugins(registry: &BuiltinRegistry, batch: &mut [InputVariant]) {
    report_plugin_failures(registry.run_batch(batch));
}

/// Reports each plugin that failed on a batch to the warning sink
/// ([`crate::warnings::Sink::plugin_failure`]), which holds them for one line
/// per plugin, written at the run's end or, on the library path, as the batch
/// ends ([`Annotator::run_batch`]). The reason is a `Run` error's own text and
/// any other error's description: `run_plugins` warns `Plugin '<Class>' went
/// wrong: $@` when a plugin dies, and the next plugin runs
/// (OutputFactory.pm:2243-2245).
pub(crate) fn report_plugin_failures(outcome: vep_builtins::BatchOutcome) {
    for failure in outcome.failures {
        let reason = match failure.error {
            vep_builtins::PluginError::Run(reason) => reason,
            other => other.to_string(),
        };
        crate::warnings::global().plugin_failure(
            &failure.plugin,
            failure.instance,
            &reason,
            failure.failed_variants,
        );
    }
}

/// The effects configuration of an embedding caller, without HGVSg; `hgvs`
/// carries the HGVS shift, the two HGVSp forms and the `--shift_3prime` the
/// HGVSc gate and the HGVSp predicates read, which agrees with `shift_3prime`
/// here.
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
        loftee_context: LofteeContextRows::None,
        context_sequences: LofteeContextRows::None,
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
/// name's own synonyms (OutputFactory.pm 1367-1377) and, where it has none,
/// among the synonyms of the chromosome the name resolves to; a sequence with
/// no `NC_` name takes its RefSeq name the same way
/// ([`vep_effects::hgvsg::refseq_name`]).
pub(crate) fn hgvsg_options(
    config: &Config,
    synonyms: &ChromosomeSynonyms,
) -> Option<vep_effects::hgvsg::HgvsgOptions> {
    config.hgvsg.then(|| vep_effects::hgvsg::HgvsgOptions {
        shift: config.shift_hgvs,
        accessions: config.hgvsg_use_accession.then(|| refseq_names(synonyms)),
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
        cache_name_on_every_row: config.spdi_on_the_cache_name,
    })
}

/// Every name of the synonyms table with an `NC_` synonym, mapped to the first
/// such synonym in the table's order ([`vep_effects::hgvsg::nc_accession`]), in
/// the latest form the table lists that accession in
/// ([`vep_effects::hgvsg::latest_form`]), the form HGVSg takes.
fn nc_accessions(synonyms: &ChromosomeSynonyms) -> HashMap<String, String> {
    synonyms
        .names()
        .filter_map(|name| {
            let listed = synonyms.synonyms_of(name);
            let accession = vep_effects::hgvsg::nc_accession(listed.iter().map(String::as_str))?;
            let accession = vep_effects::hgvsg::latest_form(
                accession,
                listed
                    .iter()
                    .map(String::as_str)
                    .chain(std::iter::once(name)),
            );
            Some((name.to_string(), accession.to_string()))
        })
        .collect()
}

/// Every name of the synonyms table with a RefSeq synonym, mapped to its `NC_`
/// accession or, where it has none, to its first synonym in the form of a RefSeq
/// accession ([`vep_effects::hgvsg::refseq_name`]); a name that is a later
/// form of that accession is mapped to itself.
fn refseq_names(synonyms: &ChromosomeSynonyms) -> HashMap<String, String> {
    synonyms
        .names()
        .filter_map(|name| {
            let refseq = vep_effects::hgvsg::refseq_name(
                synonyms.synonyms_of(name).iter().map(String::as_str),
            )?;
            let refseq = vep_effects::hgvsg::later_form(refseq, name).unwrap_or(refseq);
            Some((name.to_string(), refseq.to_string()))
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
    load_transcript_indexes_with_exon_ids(json_cache_path, effects_config, transcript_set, false)
}

/// [`load_transcript_indexes`] keeping each exon's stable id on its transcript
/// under `keep_exon_ids` ([`vep_core::transcript::Transcript::exon_ids`]), for a
/// run whose plugins print them (`exon_ids_wanted`).
pub fn load_transcript_indexes_with_exon_ids(
    json_cache_path: &str,
    effects_config: &vep_effects::EffectsConfig,
    transcript_set: TranscriptSet,
    keep_exon_ids: bool,
) -> anyhow::Result<Arc<LazyTranscriptIndexes>> {
    info!("Loading transcripts from JSON cache: {}", json_cache_path);
    let mut raw = json_cache::load_all_transcripts_in_set_with(
        json_cache_path,
        Some(&transcript_set),
        keep_exon_ids,
    )?;
    if let Some(reference) = effects_config.reference_fasta.as_deref() {
        for transcripts in raw.values_mut() {
            json_cache::pair_completed_stop_codons(transcripts, reference);
        }
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
    let source_chromosomes = crate::pipeline::source_chromosomes(
        resources.custom_sources.as_deref(),
        resources.regulatory_source.as_deref(),
    );
    for variant in batch.iter_mut() {
        crate::pipeline::name_source_chromosomes(
            variant,
            &resources.synonyms,
            &resources.transcripts,
            &source_chromosomes,
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
        if let Some(source) = &resources.regulatory_source {
            source.annotate_variant(
                variant,
                &resources.effects_config,
                resources.cell_types.as_deref(),
            );
        }

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
                (
                    1u64,
                    variant.transcript_consequences.len() as u64
                        + u64::from(variant.intergenic_overlap),
                )
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
            stats.consequences_calculated += variant.transcript_consequences.len() as u64
                + u64::from(variant.intergenic_overlap);
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
        run_builtin_plugins(&resources.builtin_plugins, batch);
        let plugin_elapsed = plugin_phase_start.elapsed();
        debug!(
            batch_size = batch_len,
            plugin_ms = plugin_elapsed.as_millis() as u64,
            "batch plugin phase complete"
        );
    }
    // A dylib plugin's failure is a warning inside `run`; the error is the host's
    // own serialisation of a row, logged, and the batch goes on.
    if let Err(e) = crate::dylib_plugins::run(&resources.dylib_plugins, batch) {
        tracing::error!("dylib plugin execution error: {:#}", e);
    }
    // A custom source's read failure lands in the batch's slot, as a variation
    // file's does, under the `custom` component.
    if let Some(sources) = &resources.custom_sources {
        if let Err(error) = sources.run_batch(batch) {
            failure.record_under("custom", error.to_string());
        }
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
        annotator_with(transcripts, nearest, BuiltinRegistry::new())
    }

    fn annotator_with(
        transcripts: Arc<LazyTranscriptIndexes>,
        nearest: Option<NearestType>,
        builtin_plugins: BuiltinRegistry,
    ) -> anyhow::Result<Annotator> {
        let pool = build_thread_pool(1).unwrap();
        let resources = AnnotationResources {
            transcripts,
            variation_source: None,
            regulatory_source: None,
            cell_types: None,
            sv_source: None,
            colocation: ColocationRules::default(),
            effects_config: Arc::new(vep_effects::EffectsConfig::default()),
            builtin_plugins: Arc::new(builtin_plugins),
            dylib_plugins: Arc::new(PluginRegistry::new()),
            prediction_config: PredictionConfig::default(),
            nearest,
            transcript_filter: None,
            synonyms: Arc::new(ChromosomeSynonyms::new()),
            custom_sources: None,
            custom_source_infos: Vec::new(),
        };
        Annotator::new(resources, pool)
    }

    fn empty_annotator() -> Annotator {
        annotator_over(Arc::new(LazyTranscriptIndexes::empty())).unwrap()
    }

    /// The messages of the `WARN` events `work` emits on this thread.
    fn warnings_of(work: impl FnOnce()) -> Vec<String> {
        messages_of(tracing::Level::WARN, work)
    }

    /// The messages of the events of `level` that `work` emits on this thread.
    fn messages_of(level: tracing::Level, work: impl FnOnce()) -> Vec<String> {
        use std::sync::Mutex;
        use tracing::field::{Field, Visit};
        use tracing_subscriber::layer::{Context, Layer, SubscriberExt};

        struct Capture {
            level: tracing::Level,
            lines: Arc<Mutex<Vec<String>>>,
        }
        struct Message(String);
        impl Visit for Message {
            fn record_debug(&mut self, field: &Field, value: &dyn std::fmt::Debug) {
                if field.name() == "message" {
                    self.0 = format!("{value:?}");
                }
            }
        }
        impl<S: tracing::Subscriber> Layer<S> for Capture {
            fn on_event(&self, event: &tracing::Event<'_>, _ctx: Context<'_, S>) {
                if *event.metadata().level() != self.level {
                    return;
                }
                let mut message = Message(String::new());
                event.record(&mut message);
                self.lines.lock().unwrap().push(message.0);
            }
        }
        let lines = Arc::new(Mutex::new(Vec::new()));
        let subscriber = tracing_subscriber::registry().with(Capture {
            level,
            lines: Arc::clone(&lines),
        });
        tracing::subscriber::with_default(subscriber, work);
        let captured = lines.lock().unwrap().clone();
        captured
    }

    /// A REVEL data file of one row, bgzipped and tabix-indexed on its `hg19_pos`
    /// column, the layout the GRCh37 preparation gives it.
    fn revel_fixture(dir: &std::path::Path) -> std::path::PathBuf {
        tabix_fixture(
            dir,
            "revel.tsv.gz",
            "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid",
            "21\t100200\t500600\tA\tG\tM\tV\t0.742\tENST00000400000",
        )
    }

    /// An AlphaMissense data file of one row, bgzipped and tabix-indexed on its
    /// `POS` column, as the published preparation indexes it.
    fn alphamissense_fixture(dir: &std::path::Path) -> std::path::PathBuf {
        tabix_fixture(
            dir,
            "alphamissense.tsv.gz",
            "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class",
            "chr21\t100200\tA\tG\thg19\tQ9NQ38\tENST00000400000.1\tM1V\t0.9123\tlikely_pathogenic",
        )
    }

    /// A one-row tabix file of `header` and `row`, indexed on the second column.
    fn tabix_fixture(
        dir: &std::path::Path,
        name: &str,
        header: &str,
        row: &str,
    ) -> std::path::PathBuf {
        tabix_fixture_in_blocks(dir, name, header, &[row], false).0
    }

    /// [`tabix_fixture`] over `rows`, in position order; under `a_block_per_row`
    /// the header and each row end their own bgzf block. The second value is
    /// the byte offset of the block holding each row.
    fn tabix_fixture_in_blocks(
        dir: &std::path::Path,
        name: &str,
        header: &str,
        rows: &[&str],
        a_block_per_row: bool,
    ) -> (std::path::PathBuf, Vec<usize>) {
        use noodles::bgzf;
        use noodles::core::Position;
        use noodles::csi::binning_index::index::header::Builder as HeaderBuilder;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        use std::io::Write as _;

        let path = dir.join(name);
        let mut writer = bgzf::io::Writer::new(std::fs::File::create(&path).unwrap());
        writeln!(writer, "{header}").unwrap();
        let mut indexer = noodles::tabix::index::Indexer::default();
        indexer.set_header(
            HeaderBuilder::default()
                .set_reference_sequence_name_index(0)
                .set_start_position_index(1)
                .build(),
        );
        let mut block_offsets = Vec::with_capacity(rows.len());
        for row in rows {
            if a_block_per_row {
                writer.flush().unwrap();
            }
            let start = writer.virtual_position();
            writeln!(writer, "{row}").unwrap();
            let end = writer.virtual_position();
            let mut fields = row.split('\t');
            let chromosome = fields.next().unwrap();
            let position: usize = fields.next().unwrap().parse().unwrap();
            let position = Position::try_from(position).unwrap();
            indexer
                .add_record(chromosome, position, position, Chunk::new(start, end))
                .unwrap();
            block_offsets.push(start.compressed() as usize);
        }
        writer.finish().unwrap();
        noodles::tabix::fs::write(format!("{}.tbi", path.display()), &indexer.build()).unwrap();
        (path, block_offsets)
    }

    /// An assembly name that gives the built-in plugins no hint (`GRCm39`) is
    /// warned about once per built-in that loads under it, naming the plugin, and
    /// not at all when no plugin loads, when the plugin fails to load, or when the
    /// assembly is one the plugins recognise; so a run without `--plugin` under
    /// such an assembly prints no line about it. REVEL refuses to load without
    /// a hint, as `REVEL.pm:133` dies without an assembly it knows, so its run
    /// under `GRCm39` is a failed load with no hint line.
    #[test]
    fn the_assembly_hint_warning_fires_only_when_a_plugin_loads() {
        let _guard = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let dir_plugins = &PluginDirs::expand(dir.path().to_str().unwrap());
        let alphamissense = format!(
            "AlphaMissense,file={}",
            alphamissense_fixture(dir.path()).display()
        );
        let revel = format!(
            "REVEL,file={},no_match=1",
            revel_fixture(dir.path()).display()
        );
        let hint_lines = |lines: Vec<String>| -> Vec<String> {
            lines
                .into_iter()
                .filter(|line| line.starts_with("--assembly"))
                .collect()
        };

        let none = warnings_of(|| {
            load_plugins(
                &[],
                dir_plugins,
                false,
                None,
                Some("GRCm39"),
                RunContext::default(),
            )
            .unwrap();
        });
        assert!(none.is_empty(), "{none:?}");

        let loaded = warnings_of(|| {
            let plugins = load_plugins(
                std::slice::from_ref(&alphamissense),
                dir_plugins,
                false,
                None,
                Some("GRCm39"),
                RunContext::default(),
            )
            .unwrap();
            assert_eq!(plugins.builtin.plugin_count(), 1);
        });
        assert_eq!(
            hint_lines(loaded),
            ["--assembly 'GRCm39' is not an assembly AlphaMissense recognises; it runs without an assembly hint"]
        );

        let recognised = warnings_of(|| {
            let plugins = load_plugins(
                &[alphamissense.clone(), revel.clone()],
                dir_plugins,
                false,
                None,
                Some("GRCh37"),
                RunContext::default(),
            )
            .unwrap();
            assert_eq!(plugins.builtin.plugin_count(), 2);
        });
        assert!(hint_lines(recognised).is_empty());

        for failing in [revel, "REVEL,file=/no/such/revel.tsv.gz".to_string()] {
            let failed = warnings_of(|| {
                let plugins = load_plugins(
                    std::slice::from_ref(&failing),
                    dir_plugins,
                    false,
                    None,
                    Some("GRCm39"),
                    RunContext::default(),
                )
                .unwrap();
                assert_eq!(plugins.builtin.plugin_count(), 0);
            });
            assert!(hint_lines(failed).is_empty());
        }
    }

    /// A registry that loaded the built-ins `names` name, every name claimed.
    fn registry_loading(names: &[&str]) -> BuiltinRegistry {
        let mut registry = BuiltinRegistry::new();
        let args: Vec<String> = names.iter().map(|name| name.to_string()).collect();
        assert!(registry.load_from_args(&args).unwrap().is_empty());
        registry
    }

    /// The context rows follow the names of the built-ins the registry loaded,
    /// not the arguments: `LoF` loads the built-in named `LoFTEE`, so it asks
    /// for every row as `LoFTEE` and `NMD` do; UTRAnnotator alone asks for the
    /// `5_prime_UTR_variant` rows, and beside LoFTEE for every row; an empty
    /// registry asks for none. LoFTEE loads under the percentile rule, the one
    /// that opens no GERP file.
    #[test]
    fn loftee_context_rows_follow_the_loaded_plugin_names() {
        let rows_under = |names: &[&str]| {
            let registry = registry_loading(names);
            assert_eq!(registry.plugin_count(), names.len());
            loftee_context_rows(&registry)
        };
        assert_eq!(rows_under(&[]), LofteeContextRows::None);
        assert_eq!(
            rows_under(&["LoF,use_gerp_end_trunc=0"]),
            LofteeContextRows::Every
        );
        assert_eq!(
            rows_under(&["LoFTEE,use_gerp_end_trunc=0"]),
            LofteeContextRows::Every
        );
        assert_eq!(rows_under(&["NMD"]), LofteeContextRows::Every);
        assert_eq!(rows_under(&["NearestExonJB"]), LofteeContextRows::Every);
        assert_eq!(
            rows_under(&["UTRAnnotator"]),
            LofteeContextRows::WithTerm(Consequence::FivePrimeUtrVariant)
        );
        assert_eq!(
            rows_under(&["UTRAnnotator", "LoF,use_gerp_end_trunc=0"]),
            LofteeContextRows::Every
        );
        assert_eq!(
            rows_under(&["Downstream", "TSSDistance"]),
            LofteeContextRows::None
        );
    }

    /// The transcript's sequences travel on the `frameshift_variant` rows when
    /// Downstream is loaded and on no row otherwise.
    #[test]
    fn context_sequence_rows_follow_the_downstream_plugin() {
        let rows_under = |names: &[&str]| context_sequence_rows(&registry_loading(names));
        assert_eq!(rows_under(&[]), LofteeContextRows::None);
        assert_eq!(
            rows_under(&["NMD", "LoF,use_gerp_end_trunc=0"]),
            LofteeContextRows::None
        );
        assert_eq!(
            rows_under(&["Downstream"]),
            LofteeContextRows::WithTerm(Consequence::FrameshiftVariant)
        );
        assert_eq!(
            rows_under(&["Blosum62", "Downstream"]),
            LofteeContextRows::WithTerm(Consequence::FrameshiftVariant)
        );
    }

    /// The exon ids load for a run whose `--plugin` arguments name
    /// NearestExonJB, bare or with parameters, read as the registry reads a
    /// name (the text before the first comma), and for no other run.
    #[test]
    fn exon_ids_load_only_for_a_run_naming_nearest_exon_jb() {
        let wanted = |args: &[&str]| {
            let args: Vec<String> = args.iter().map(|arg| arg.to_string()).collect();
            exon_ids_wanted(&args)
        };
        assert!(wanted(&["NearestExonJB"]));
        assert!(wanted(&[
            "Downstream",
            "NearestExonJB,intronic=1,max_range=50"
        ]));
        assert!(!wanted(&[]));
        assert!(!wanted(&["Downstream", "TSSDistance"]));
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
                                tc.transcript_id(),
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
                        tc.transcript_id(),
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
            regulatory_source: None,
            cell_types: None,
            sv_source: None,
            colocation: ColocationRules::default(),
            effects_config: Arc::new(vep_effects::EffectsConfig::default()),
            builtin_plugins: Arc::new(BuiltinRegistry::new()),
            dylib_plugins: Arc::new(PluginRegistry::new()),
            prediction_config: PredictionConfig::default(),
            nearest: None,
            transcript_filter: None,
            synonyms: Arc::new(ChromosomeSynonyms::new()),
            custom_sources: None,
            custom_source_infos: Vec::new(),
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
                        tc.transcript_id(),
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

    /// A built-in plugin's failure on a batch a program hands over is written as
    /// the batch ends, one `went wrong` line carrying the variants it failed on:
    /// the sink holds a command-line run's plugin failures to the run's end,
    /// which a program has none of. REVEL over a file whose second data block
    /// cannot be inflated loads on the first (`last_preamble_line` reads through
    /// the first data line) and fails the variant whose span is read from the
    /// second.
    #[test]
    fn the_library_path_writes_a_plugins_failure_as_the_batch_ends() {
        use crate::warnings::Target;

        let _guard = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let (path, row_blocks) = tabix_fixture_in_blocks(
            dir.path(),
            "revel.tsv.gz",
            "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid",
            &[
                "21\t100200\t500600\tA\tG\tM\tV\t0.742\tENST00000400000",
                "21\t100300\t500700\tC\tT\tR\tW\t0.311\tENST00000400000",
            ],
            true,
        );
        let mut revel = std::fs::read(&path).unwrap();
        // The second row's block runs to the 28-byte end-of-file block, its
        // deflate payload between the 18-byte block header and the 8-byte
        // trailer. Zeroed, the payload is a stored block whose length check fails.
        let payload = row_blocks[1] + 18..revel.len() - 28 - 8;
        revel[payload].fill(0);
        std::fs::write(&path, &revel).unwrap();

        let plugins = load_plugins(
            &[format!("REVEL,file={}", path.display())],
            &PluginDirs::expand(dir.path().to_str().unwrap()),
            false,
            None,
            Some("GRCh37"),
            RunContext::default(),
        )
        .unwrap();
        assert_eq!(plugins.builtin.plugin_count(), 1);
        let annotator = annotator_with(
            Arc::new(LazyTranscriptIndexes::empty()),
            None,
            plugins.builtin,
        )
        .unwrap();

        let sink = crate::warnings::global();
        sink.configure(Target::Memory, true, false);
        let mut batch = parse_vcf_line("21\t100300\tin_block\tC\tT\t.\t.\t.", true).unwrap();
        let _ = annotator.annotate_batch(&mut batch);
        // Another test's warning can land in the same process-wide sink while
        // the window is open, so only the plugin lines are counted.
        let texts: Vec<String> = sink
            .take_collected()
            .into_iter()
            .map(|(_, text)| text)
            .filter(|text| text.starts_with("WARNING: Plugin '"))
            .collect();
        sink.configure(Target::Stderr, false, false);
        assert_eq!(texts.len(), 1, "{texts:#?}");
        assert!(
            texts[0].starts_with("WARNING: Plugin 'REVEL' went wrong: "),
            "{}",
            texts[0]
        );
        assert!(
            texts[0].ends_with(
                "; the plugin annotated nothing on the variants it failed on, the first at \
                 21:100300; 1 variant affected\n"
            ),
            "{}",
            texts[0]
        );
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
            ..TranscriptSet::default()
        });
        let primary = count(TranscriptSet {
            gencode_primary: true,
            ..TranscriptSet::default()
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
                .map(|tc| (tc.transcript_id().to_string(), tc.consequences.to_vec()))
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
                .map(|tc| (tc.transcript_id().to_string(), tc.consequences.to_vec()))
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
                cache_name_on_every_row: false,
            })
        );
        let config = Config::from_argv(["vep", "--spdi"]).unwrap();
        assert_eq!(
            spdi_options(&config, &table),
            Some(vep_effects::spdi::SpdiOptions {
                vrs_accessions: None,
                cache_name_on_every_row: false,
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

    /// Ensembl VEP 116.2 prints the `--everything` notice (Runner.pm:752,
    /// BaseVEP.pm:785) before it throws for `--hgvsg` without a FASTA (:759), so
    /// `from_config` logs it at `INFO` before it returns that refusal, as the
    /// binary prints it before the error.
    #[test]
    fn from_config_logs_the_everything_notice_before_the_hgvsg_refusal() {
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
            "--everything",
            "--hgvsg",
        ])
        .unwrap();
        let mut refusal = None;
        let logged = messages_of(tracing::Level::INFO, || {
            refusal = Annotator::from_config(&config).err();
        });
        let err = refusal.expect("--everything --hgvsg without a FASTA is refused");
        assert_eq!(
            CliError::from_run(&err),
            format!("ERROR: [fasta] {HGVS_NEEDS_FASTA}")
        );
        assert_eq!(
            logged
                .iter()
                .filter(|m| *m == "Disabling --hgvs; using --offline and no FASTA file found")
                .count(),
            1,
            "{logged:?}"
        );
    }

    /// Ensembl VEP 116.2 prints the `--var_synonyms` notice as `init` begins
    /// (Runner.pm:117-120), before the annotation sources are built (:128),
    /// where a missing cache directory is refused (CacheDir.pm:109) and the
    /// transcript filter is compiled (AnnotationSource.pm:87-90), so
    /// `from_config` logs it at `INFO` once before it returns either refusal.
    #[test]
    fn from_config_logs_the_var_synonyms_notice_before_a_cache_or_filter_refusal() {
        const NOTICE: &str = "Colocated variants will be provided when using --var_synonyms";
        let leaf = format!(
            "{}/../../tests/golden/116/GRCh38-release116/json_cache",
            env!("CARGO_MANIFEST_DIR")
        );
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("missing_cache").display().to_string();
        for (argv, error) in [
            (
                vec![
                    "vep",
                    "--fork",
                    "1",
                    "--json_cache",
                    &missing,
                    "--var_synonyms",
                ],
                format!("ERROR: [cache] Directory {missing} does not exist"),
            ),
            (
                vec![
                    "vep",
                    "--fork",
                    "1",
                    "--json_cache",
                    &leaf,
                    "--var_synonyms",
                    "--transcript_filter",
                    "biotype equals x",
                ],
                "ERROR: [filter] No such operator \"equals\"".to_string(),
            ),
        ] {
            let config = Config::from_argv(argv.clone()).unwrap();
            let mut refusal = None;
            let logged = messages_of(tracing::Level::INFO, || {
                refusal = Annotator::from_config(&config).err();
            });
            let err = refusal.unwrap_or_else(|| panic!("{argv:?} was not refused"));
            assert_eq!(CliError::from_run(&err), error, "{argv:?}");
            assert_eq!(
                logged.iter().filter(|m| *m == NOTICE).count(),
                1,
                "{argv:?}: {logged:?}"
            );
        }
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
                        tc.transcript_id().to_string(),
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

    /// `from_config` gates the per-row plugin context as the binary does: with
    /// Downstream, HGVSIntronOffset and NearestExonJB loaded under `--hgvs` on
    /// the GRCh38-plugin-computed corpus, every row carries the transcript
    /// model, the `frameshift_variant` rows carry the transcript's sequences
    /// and the HGVS notations are computed, so the record `21 36896882 TG/T`
    /// gets on ENST00000336648 the fields the corpus documents for the row:
    /// `DownstreamProtein=PQRCVSWMG`, which starts at the codon the deletion
    /// changes where the reference row's `QRCVSWMG` starts one codon later
    /// (`Downstream.pm:113`), and the reference row's
    /// `ProteinLengthChange=-242`,
    /// `NearestExonJB=ENSE00001282176|23|start|272` and both intron offsets `0`.
    #[test]
    fn from_config_gates_the_plugin_context_as_the_binary_does() {
        let corpus = format!(
            "{}/../../tests/golden/116/GRCh38-plugin-computed",
            env!("CARGO_MANIFEST_DIR")
        );
        let cache = format!("{corpus}/json_cache");
        let fasta = format!("{corpus}/reference.fa.gz");
        let config = Config::from_argv([
            "vep",
            "--fork",
            "1",
            "--json_cache",
            &cache,
            "--fasta",
            &fasta,
            "--hgvs",
            "--plugin",
            "Downstream",
            "--plugin",
            "HGVSIntronOffset",
            "--plugin",
            "NearestExonJB",
        ])
        .unwrap();
        let annotator = Annotator::from_config(&config).unwrap();
        assert_eq!(annotator.plugin_count(), 3);
        let effects = &annotator.resources().effects_config;
        assert_eq!(effects.loftee_context, LofteeContextRows::Every);
        assert_eq!(
            effects.context_sequences,
            LofteeContextRows::WithTerm(Consequence::FrameshiftVariant)
        );
        assert!(effects.compute_hgvs);

        let mut batch = parse_vcf_line("21\t36896882\t1069237\tTG\tT\t.\t.\t.", false).unwrap();
        let _ = annotator.annotate_batch(&mut batch);
        let row = batch[0]
            .transcript_consequences
            .iter()
            .find(|tc| tc.transcript_id() == "ENST00000336648")
            .expect("the frameshift row");
        assert!(row.consequences.contains(&Consequence::FrameshiftVariant));
        assert_eq!(row.plugin_data["DownstreamProtein"], "PQRCVSWMG");
        assert_eq!(row.plugin_data["ProteinLengthChange"], "-242");
        assert_eq!(
            row.plugin_data["NearestExonJB"],
            "ENSE00001282176|23|start|272"
        );
        assert_eq!(row.plugin_data["HGVS_IntronStartOffset"], "0");
        assert_eq!(row.plugin_data["HGVS_IntronEndOffset"], "0");
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
        let dirs = PluginDirs::expand("/no/such/dir");
        let err = load_plugins(
            &["NoSuchPlugin".to_string()],
            &dirs,
            true,
            None,
            None,
            RunContext::default(),
        )
        .err()
        .expect("a missing plugins directory is fatal under --safe");
        assert_eq!(
            err.downcast_ref::<crate::config::CliError>()
                .unwrap()
                .to_string(),
            "ERROR: [plugin] Plugins directory '/no/such/dir' not found."
        );
        let loaded = load_plugins(&[], &dirs, true, None, None, RunContext::default()).unwrap();
        assert_eq!(loaded.builtin.plugin_count(), 0);
    }

    /// `Name,param,param` splits on every comma (Runner.pm:928: `split /,/`), the
    /// first part the plugin, the rest its parameters in order: LoFtool reads its
    /// scores file from the first parameter, so the file it does not find names
    /// the parameter the split gave it. The failure is fatal under `--safe` alone
    /// (`:946-954`: `throw($msg) if $self->param('safe')`, else a warning and
    /// `next`), the plugin component carrying Ensembl VEP's text; without it the
    /// run goes on with no plugin loaded. Without a parameter the plugin reads
    /// its bundled file from the plugins directory, which reaches it through the
    /// run context.
    #[test]
    fn plugin_name_and_params_split_as_perl_does() {
        // The LoFtool load below warns through the process-wide sink.
        let _guard = crate::warnings::sink_test_guard();
        let dir = tempfile::tempdir().unwrap();
        let dirs = PluginDirs::expand(dir.path().to_str().unwrap());
        for (arg, first_param) in [
            ("LoFtool,/no/such/scores.txt", "/no/such/scores.txt"),
            ("LoFtool,/no/such/first,second,third", "/no/such/first"),
        ] {
            let err = load_plugins(
                &[arg.to_string()],
                &dirs,
                true,
                None,
                None,
                RunContext::default(),
            )
            .err()
            .expect("a scores file that does not exist fails init under --safe");
            assert_eq!(
                CliError::from_run(&err),
                format!(
                    "ERROR: [plugin] Failed to instantiate plugin LoFtool: LoFtool: scores file {first_param} not found"
                ),
                "{arg}"
            );
            let loaded = load_plugins(
                &[arg.to_string()],
                &dirs,
                false,
                None,
                None,
                RunContext::default(),
            )
            .unwrap();
            assert_eq!(loaded.builtin.plugin_count(), 0);
        }
        let err = load_plugins(
            &["LoFtool".to_string()],
            &dirs,
            true,
            None,
            None,
            RunContext::default(),
        )
        .err()
        .expect("LoFtool without a parameter fails init under --safe");
        assert_eq!(
            CliError::from_run(&err),
            format!(
                "ERROR: [plugin] Failed to instantiate plugin LoFtool: LoFtool: scores file {} not found",
                dir.path().join("LoFtool_scores.txt").display()
            )
        );
    }
}
