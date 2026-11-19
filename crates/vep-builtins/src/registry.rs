// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Plugin registry: discovers and manages built-in plugins by name.
//!
//! The registry resolves `--plugin Name,param1,param2` arguments to built-in
//! plugin instances. If a name is not found in the registry, the caller
//! falls back to the dylib plugin loader.
//!
//! Plugin execution uses a two-phase protocol when possible:
//! 1. **Prefetch** (parallel): all plugins that support it run I/O in
//!    parallel via rayon, reading from independent data files. A plugin whose
//!    prefetch reads the records' coordinates alone may be prefetched before
//!    the batch is annotated ([`BuiltinRegistry::prefetch_ahead`]); one that
//!    reads the consequence rows is prefetched once they exist.
//! 2. **Annotate** (sequential): prefetched results are applied to variants.
//!
//! Plugins that don't implement prefetch fall back to sequential `run_batch`.

use std::time::Instant;

use rayon::prelude::*;
use tracing::{debug, info};

use crate::plugins;
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData, RunContext};
use vep_core::consequence::RowKind;
use vep_core::sv_routing::split_plugin_arg;
use vep_core::variant::InputVariant;

/// What one plugin's prefetch returned.
type Prefetched = Result<Option<PrefetchData>, PluginError>;

/// One batch's prefetch results, one slot per loaded plugin in load order: a
/// filled slot holds what [`BuiltinPlugin::prefetch`] returned for the batch,
/// an empty one marks a plugin whose prefetch still has to run, which
/// [`BuiltinRegistry::annotate_batch`] does once the batch is annotated.
pub struct PluginPrefetch(Vec<Option<Prefetched>>);

impl PluginPrefetch {
    /// Every plugin's prefetch still to run.
    fn unread(plugins: usize) -> Self {
        Self((0..plugins).map(|_| None).collect())
    }
}

/// Per-plugin timing measurements from a single `run_batch_timed` call.
#[derive(Debug, Clone, serde::Serialize)]
pub struct PluginTimings {
    /// Per-plugin prefetch durations: `(plugin_name, milliseconds)`.
    /// For plugins that don't support prefetch, the value is 0.
    pub prefetch_ms: Vec<(String, u64)>,
    /// Per-plugin annotate/run_batch durations: `(plugin_name, milliseconds)`.
    pub annotate_ms: Vec<(String, u64)>,
    /// Total wall-clock time for the parallel prefetch phase.
    pub total_prefetch_ms: u64,
    /// Total wall-clock time for the sequential annotate phase.
    pub total_annotate_ms: u64,
    /// Total wall-clock time for the entire plugin batch (prefetch + annotate).
    pub total_ms: u64,
}

/// What loading the `--plugin` arguments did not load: the arguments no built-in
/// claims, as written, and the built-ins whose `init` failed, each with the name
/// the argument gave and the error.
#[derive(Debug, Default)]
pub struct LoadOutcome {
    pub unresolved: Vec<String>,
    pub failed: Vec<(String, PluginError)>,
}

/// What one batch's plugin phase could not do: every plugin whose `prefetch`,
/// `annotate` or `run_batch` failed on a variant after the variant-by-variant
/// retry, with the name it reports and its error, in load order, one entry per
/// plugin per batch. A plugin that fails costs its own fields alone, on the
/// variants it failed on (the whole variant, or the rows it failed on, as the
/// error's note states): the plugins after it run as if it were absent, as
/// Ensembl VEP's `run_plugins` catches each plugin's failure per line and goes
/// on to the next plugin (OutputFactory.pm:2212-2256). Empty when every plugin
/// annotated every variant.
#[derive(Debug, Default)]
#[must_use = "a failed plugin is named here and goes unreported when dropped"]
pub struct BatchOutcome {
    pub failures: Vec<PluginFailure>,
}

/// One plugin's failure on a batch: the plugin's name and its place in load
/// order (a built-in loaded twice is two plugins of one name, each with its own
/// failures), the error the first variant it failed on met, with the note of
/// the phase's [`Footprint`] naming that variant, and the number of the batch's
/// variants it failed on, which the run's warning sums over the batches into
/// the one line it writes for the plugin.
#[derive(Debug)]
pub struct PluginFailure {
    pub plugin: String,
    pub instance: usize,
    pub error: PluginError,
    pub failed_variants: usize,
}

/// What a variant a plugin failed on keeps of the plugin's fields, by the
/// phase that failed; the note the plugin's failure carries states it.
#[derive(Debug, Clone, Copy)]
enum Footprint {
    /// Nothing: a failed `prefetch` runs no `annotate`, and an `annotate`
    /// that fails writes nothing on the variant
    /// ([`BuiltinPlugin::annotate`]).
    Variant,
    /// Every row but the ones the plugin failed on: `run_batch` runs every
    /// row and reports the first failure after the loop
    /// ([`BuiltinPlugin::run_batch`]), as does an `annotate` whose plugin
    /// declares its failures per row ([`BuiltinPlugin::annotate_fails_per_row`]).
    Rows,
}

impl Footprint {
    /// The note a failure carries after its error text, naming `first`, the
    /// location of the first variant of the batch the plugin failed on.
    fn note(self, first: &str) -> String {
        match self {
            Footprint::Variant => format!(
                "the plugin annotated nothing on the variants it failed on, the first at {first}"
            ),
            Footprint::Rows => format!(
                "the plugin annotated every row but the ones it failed on, the first of them \
                 on the variant at {first}"
            ),
        }
    }
}

/// Phase two for one plugin over `variants`: the data its prefetch returned
/// applied, or `run_batch` when it prefetched nothing; a failed prefetch is the
/// failure. An error comes with the footprint of the phase that failed, which
/// for `annotate` is the rows' when the plugin declares its failures per row
/// ([`BuiltinPlugin::annotate_fails_per_row`]).
fn apply(
    plugin: &dyn BuiltinPlugin,
    variants: &mut [InputVariant],
    prefetched: Result<Option<PrefetchData>, PluginError>,
) -> Result<(), (Footprint, PluginError)> {
    match prefetched {
        Ok(Some(data)) => plugin.annotate(variants, data).map_err(|error| {
            let footprint = if plugin.annotate_fails_per_row() {
                Footprint::Rows
            } else {
                Footprint::Variant
            };
            (footprint, error)
        }),
        Ok(None) => plugin
            .run_batch(variants)
            .map_err(|error| (Footprint::Rows, error)),
        Err(error) => Err((Footprint::Variant, error)),
    }
}

/// Registry of built-in plugins.
///
/// Manages the lifecycle of built-in VEP plugins: initialization from CLI
/// arguments, two-phase batch execution (parallel prefetch + sequential
/// annotate), and timing measurement for benchmarking.
pub struct BuiltinRegistry {
    plugins: Vec<Box<dyn BuiltinPlugin>>,
}

impl BuiltinRegistry {
    /// Create an empty registry.
    pub fn new() -> Self {
        Self {
            plugins: Vec::new(),
        }
    }

    /// Phase two for one plugin: applies the data its prefetch returned, or runs
    /// it over the batch when it prefetched nothing. A failure of the batch is
    /// retried variant by variant, each variant prefetched and annotated on its
    /// own, as Ensembl VEP's plugins read each line's own span and lose the
    /// lines whose span fails alone (OutputFactory.pm:2225-2226, REVEL.pm:168);
    /// a batch of one variant is its own retry. After a `run_batch` failure
    /// the batch pass wrote every row but the failed ones ([`Footprint::Rows`]),
    /// so the variants after the first that fails on its own have nothing to
    /// gain from a retry and are run again on a copy, for the count of the
    /// variants the plugin failed on alone. The plugin is recorded in `outcome`
    /// once under its name, with the error the first failing variant met, the
    /// [`Footprint`] of the phase it failed in, the location of that variant
    /// and the number of the batch's variants it failed on; a retry that
    /// annotates every variant records nothing.
    fn apply_prefetched(
        outcome: &mut BatchOutcome,
        instance: usize,
        plugin: &dyn BuiltinPlugin,
        variants: &mut [InputVariant],
        prefetched: Result<Option<PrefetchData>, PluginError>,
    ) {
        let Err(batch_failure) = apply(plugin, variants, prefetched) else {
            return;
        };
        let (first_failure, failed_variants) = if variants.len() < 2 {
            let first = variants
                .first()
                .map(|variant| (variant.location(), batch_failure));
            (first, 1)
        } else {
            let count_after_first = matches!(batch_failure.0, Footprint::Rows);
            let mut failed_variants = 0usize;
            let mut first_failure: Option<(String, (Footprint, PluginError))> = None;
            for variant in variants.iter_mut() {
                // Past the first failure under `Footprint::Rows` the batch pass's
                // rows stand, so the variant is run on a copy, for the count alone.
                let mut copy;
                let retried_variant = if count_after_first && first_failure.is_some() {
                    copy = variant.clone();
                    &mut copy
                } else {
                    variant
                };
                let one = std::slice::from_mut(retried_variant);
                let prefetched = plugin.prefetch(one);
                if let Err(failure) = apply(plugin, one, prefetched) {
                    failed_variants += 1;
                    if first_failure.is_none() {
                        first_failure = Some((one[0].location(), failure));
                    }
                }
            }
            debug!(
                plugin = plugin.name(),
                batch_size = variants.len(),
                failed_variants,
                "failed plugin batch retried variant by variant"
            );
            (first_failure, failed_variants)
        };
        let Some((location, (footprint, error))) = first_failure else {
            return;
        };
        outcome.failures.push(PluginFailure {
            plugin: plugin.name().to_string(),
            instance,
            error: error.with_note(&footprint.note(&location)),
            failed_variants,
        });
    }

    /// Load plugins from `--plugin` argument strings.
    ///
    /// Each argument has the format `"PluginName,param1,param2"`.
    /// Returns names of plugins that were not found as builtins (for dylib
    /// fallback); the first built-in whose `init` fails is the error.
    pub fn load_from_args(&mut self, plugin_args: &[String]) -> Result<Vec<String>, PluginError> {
        let outcome = self.load_from_args_with_context(plugin_args, None, None);
        match outcome.failed.into_iter().next() {
            Some((_, error)) => Err(error),
            None => Ok(outcome.unresolved),
        }
    }

    /// Load plugins with the shared reference FASTA (from `--fasta`) and the
    /// target assembly (from `--assembly`), under the default run context: the
    /// `vep` output format and the Ensembl transcript set.
    pub fn load_from_args_with_context(
        &mut self,
        plugin_args: &[String],
        reference_fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>,
        assembly: Option<vep_core::assembly::Assembly>,
    ) -> LoadOutcome {
        self.load_from_args_for_run(
            plugin_args,
            reference_fasta,
            assembly,
            RunContext::default(),
        )
    }

    /// Load plugins with the full pipeline context: the shared reference FASTA
    /// (from `--fasta`), the target assembly (from `--assembly`) and the run's
    /// output format and transcript set (`--refseq`, `--merged`).
    ///
    /// The assembly, the format and the transcript set are injected **before**
    /// `init()` so a plugin can use them while opening its data files, as Ensembl
    /// VEP's plugins read `$self->{config}` in their constructors. The assembly's
    /// ordering is load-bearing: REVEL must know the assembly to pick which
    /// coordinate column tabix indexed, and LoFTEE derives its default END_TRUNC
    /// cutoff from it. The FASTA is injected after `init()`. A built-in whose
    /// `init` fails is left out and reported in the outcome with its name; the
    /// arguments after it still load.
    pub fn load_from_args_for_run(
        &mut self,
        plugin_args: &[String],
        reference_fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>,
        assembly: Option<vep_core::assembly::Assembly>,
        run: RunContext<'_>,
    ) -> LoadOutcome {
        let mut outcome = LoadOutcome::default();

        for arg in plugin_args {
            let (name, params) = split_plugin_arg(arg);
            let Some(mut plugin) = create_builtin(name) else {
                outcome.unresolved.push(arg.clone());
                continue;
            };
            if let Err(error) = load_builtin(
                plugin.as_mut(),
                &params,
                reference_fasta.clone(),
                assembly,
                run,
            ) {
                outcome.failed.push((name.to_string(), error));
                continue;
            }
            info!(name = plugin.name(), "loaded builtin plugin");
            self.plugins.push(plugin);
        }

        outcome
    }

    /// Run all plugins on a buffer of variants using two-phase execution.
    ///
    /// **Phase 1 (parallel prefetch):** Calls [`prefetch`] on every plugin
    /// in parallel via rayon. Each plugin reads from its own data file(s)
    /// using only a shared `&[InputVariant]` reference. Plugins that return
    /// `None` from prefetch are handled in phase 2 via `run_batch`.
    ///
    /// **Phase 2 (sequential annotate):** For each plugin, either applies
    /// the prefetched data via [`annotate`], or falls back to [`run_batch`].
    /// This phase needs `&mut [InputVariant]` so it runs sequentially.
    ///
    /// A plugin whose prefetch, annotate or run_batch fails on a variant after the
    /// variant-by-variant retry is named in the returned [`BatchOutcome`] and the
    /// plugins after it still run.
    ///
    /// [`prefetch`]: BuiltinPlugin::prefetch
    /// [`annotate`]: BuiltinPlugin::annotate
    /// [`run_batch`]: BuiltinPlugin::run_batch
    pub fn run_batch(&self, variants: &mut [InputVariant]) -> BatchOutcome {
        self.annotate_batch(variants, PluginPrefetch::unread(self.plugins.len()))
    }

    /// Phase 1 for the plugins whose prefetch reads no consequence row
    /// ([`BuiltinPlugin::prefetch_reads_consequences`] false), in parallel on
    /// the current rayon pool, over a batch that need not be annotated yet;
    /// the other plugins' slots are left for [`Self::annotate_batch`] to fill.
    pub fn prefetch_ahead(&self, variants: &[InputVariant]) -> PluginPrefetch {
        PluginPrefetch(
            self.plugins
                .par_iter()
                .map(|plugin| {
                    if plugin.prefetch_reads_consequences() {
                        None
                    } else {
                        Some(plugin.prefetch(variants))
                    }
                })
                .collect(),
        )
    }

    /// Phase 1 for the plugins whose slot `prefetch` left empty, in parallel, then
    /// phase 2 plugin by plugin in load order; every plugin whose prefetch,
    /// annotate or run_batch failed on a variant after the variant-by-variant
    /// retry is named in the returned [`BatchOutcome`].
    pub fn annotate_batch(
        &self,
        variants: &mut [InputVariant],
        prefetch: PluginPrefetch,
    ) -> BatchOutcome {
        let mut outcome = BatchOutcome::default();
        if self.plugins.is_empty() {
            return outcome;
        }
        debug_assert_eq!(prefetch.0.len(), self.plugins.len());

        let batch_size = variants.len();
        let overall_start = Instant::now();

        // The mutable slice is reborrowed as shared for prefetch; no plugin
        // writes to variants in that phase.
        let variants_shared: &[InputVariant] = variants;
        let prefetch_start = Instant::now();

        let prefetched: Vec<Prefetched> = self
            .plugins
            .par_iter()
            .zip(prefetch.0)
            .map(|(plugin, ahead)| ahead.unwrap_or_else(|| plugin.prefetch(variants_shared)))
            .collect();

        let prefetch_elapsed = prefetch_start.elapsed();
        debug!(
            batch_size,
            prefetch_ms = prefetch_elapsed.as_millis() as u64,
            plugins = self.plugins.len(),
            "parallel prefetch phase complete"
        );

        for (instance, (plugin, prefetch_result)) in self.plugins.iter().zip(prefetched).enumerate()
        {
            let started = Instant::now();
            Self::apply_prefetched(
                &mut outcome,
                instance,
                plugin.as_ref(),
                variants,
                prefetch_result,
            );
            let elapsed = started.elapsed();
            debug!(
                plugin = plugin.name(),
                batch_size,
                annotate_ms = elapsed.as_millis() as u64,
                "plugin annotate phase completed"
            );
        }

        let overall_elapsed = overall_start.elapsed();
        debug!(
            batch_size,
            total_plugin_ms = overall_elapsed.as_millis() as u64,
            "all plugins completed"
        );

        outcome
    }

    /// Like [`Self::run_batch`] but returns structured per-plugin timing data
    /// beside the outcome.
    ///
    /// Used by benchmarking tools to measure plugin overhead. The pipeline
    /// behavior is identical to `run_batch`.
    pub fn run_batch_timed(&self, variants: &mut [InputVariant]) -> (PluginTimings, BatchOutcome) {
        let mut timings = PluginTimings {
            prefetch_ms: Vec::with_capacity(self.plugins.len()),
            annotate_ms: Vec::with_capacity(self.plugins.len()),
            total_prefetch_ms: 0,
            total_annotate_ms: 0,
            total_ms: 0,
        };
        let mut outcome = BatchOutcome::default();

        if self.plugins.is_empty() {
            return (timings, outcome);
        }

        let overall_start = Instant::now();

        let variants_shared: &[InputVariant] = variants;
        let prefetch_wall_start = Instant::now();

        let prefetched: Vec<(Result<Option<PrefetchData>, PluginError>, u64)> = self
            .plugins
            .par_iter()
            .map(|plugin| {
                let t = Instant::now();
                let result = plugin.prefetch(variants_shared);
                let elapsed_ms = t.elapsed().as_millis() as u64;
                (result, elapsed_ms)
            })
            .collect();

        timings.total_prefetch_ms = prefetch_wall_start.elapsed().as_millis() as u64;

        for (plugin, (_result, ms)) in self.plugins.iter().zip(prefetched.iter()) {
            timings.prefetch_ms.push((plugin.name().to_string(), *ms));
        }

        let annotate_wall_start = Instant::now();

        for (instance, (plugin, (prefetch_result, _ms))) in
            self.plugins.iter().zip(prefetched).enumerate()
        {
            let t = Instant::now();
            Self::apply_prefetched(
                &mut outcome,
                instance,
                plugin.as_ref(),
                variants,
                prefetch_result,
            );
            let elapsed_ms = t.elapsed().as_millis() as u64;
            timings
                .annotate_ms
                .push((plugin.name().to_string(), elapsed_ms));
        }

        timings.total_annotate_ms = annotate_wall_start.elapsed().as_millis() as u64;
        timings.total_ms = overall_start.elapsed().as_millis() as u64;

        (timings, outcome)
    }

    /// Run all plugins **sequentially** (no parallel prefetch).
    ///
    /// Used by benchmarks to measure the speedup from parallel prefetch.
    /// Each plugin's prefetch + annotate runs one at a time; a failing plugin
    /// is named in the outcome as [`Self::run_batch`] names it.
    pub fn run_batch_sequential(
        &self,
        variants: &mut [InputVariant],
    ) -> (PluginTimings, BatchOutcome) {
        let mut timings = PluginTimings {
            prefetch_ms: Vec::with_capacity(self.plugins.len()),
            annotate_ms: Vec::with_capacity(self.plugins.len()),
            total_prefetch_ms: 0,
            total_annotate_ms: 0,
            total_ms: 0,
        };
        let mut outcome = BatchOutcome::default();

        if self.plugins.is_empty() {
            return (timings, outcome);
        }

        let overall_start = Instant::now();
        let mut total_prefetch = 0u64;
        let mut total_annotate = 0u64;

        for (instance, plugin) in self.plugins.iter().enumerate() {
            let t = Instant::now();
            let prefetch_result = plugin.prefetch(variants);
            let prefetch_ms = t.elapsed().as_millis() as u64;
            timings
                .prefetch_ms
                .push((plugin.name().to_string(), prefetch_ms));
            total_prefetch += prefetch_ms;

            let t = Instant::now();
            Self::apply_prefetched(
                &mut outcome,
                instance,
                plugin.as_ref(),
                variants,
                prefetch_result,
            );
            let annotate_ms = t.elapsed().as_millis() as u64;
            timings
                .annotate_ms
                .push((plugin.name().to_string(), annotate_ms));
            total_annotate += annotate_ms;
        }

        timings.total_prefetch_ms = total_prefetch;
        timings.total_annotate_ms = total_annotate;
        timings.total_ms = overall_start.elapsed().as_millis() as u64;

        (timings, outcome)
    }

    /// Returns the names of all loaded plugins.
    pub fn plugin_names(&self) -> Vec<String> {
        self.plugins.iter().map(|p| p.name().to_string()).collect()
    }

    /// Returns the number of loaded plugins.
    pub fn plugin_count(&self) -> usize {
        self.plugins.len()
    }

    /// Returns `true` if no plugins are loaded.
    pub fn is_empty(&self) -> bool {
        self.plugins.is_empty()
    }

    /// Returns header info from all loaded plugins: plugin by plugin in load
    /// order, each plugin's fields sorted by name, the order Ensembl VEP's
    /// `OutputFactory.pm` `get_plugin_headers` gives its plugin columns.
    pub fn all_header_info(&self) -> Vec<(String, String)> {
        let mut headers = Vec::new();
        for plugin in &self.plugins {
            let mut fields = plugin.header_info();
            fields.sort_by(|a, b| a.0.cmp(&b.0));
            headers.extend(fields);
        }
        headers
    }

    /// Each loaded plugin's fields with the kinds of row the plugin runs on
    /// ([`BuiltinPlugin::feature_types`]), for the writers that print a field on
    /// the rows of those kinds alone.
    pub fn field_row_kinds(&self) -> Vec<(String, Vec<RowKind>)> {
        self.plugins
            .iter()
            .flat_map(|plugin| {
                plugin
                    .header_info()
                    .into_iter()
                    .map(|(name, _)| (name, plugin.feature_types().to_vec()))
            })
            .collect()
    }
}

impl Default for BuiltinRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// Hands a created built-in the run's context in the order Ensembl VEP's plugins
/// read it: the assembly, the output format, the transcript set and the plugins
/// directories before `init` (their constructors read `$self->{config}`), the
/// shared FASTA after it. A parameter naming a remote data file under the
/// plugin's own grammar ([`BuiltinPlugin::remote_data_file`]) is refused
/// before `init` as the instantiation failure
/// `remote data files are not supported; download <url> and index it locally`:
/// Ensembl VEP's tabix plugins stream a `://` path through htslib
/// (`BaseVepTabixPlugin.pm` `check_file`), which vep-rs does not do.
fn load_builtin(
    plugin: &mut dyn BuiltinPlugin,
    params: &[String],
    reference_fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>,
    assembly: Option<vep_core::assembly::Assembly>,
    run: RunContext<'_>,
) -> Result<(), PluginError> {
    if let Some(url) = params
        .iter()
        .find_map(|param| plugin.remote_data_file(param))
    {
        return Err(PluginError::Init(format!(
            "remote data files are not supported; download {url} and index it locally"
        )));
    }
    plugin.set_assembly(assembly);
    plugin.set_output_format(run.output_format);
    plugin.set_transcript_set(run.refseq, run.merged);
    plugin.set_plugin_dirs(run.plugin_dirs);
    plugin.init(params)?;
    if plugin.requires_reference_fasta() && reference_fasta.is_none() {
        return Err(PluginError::Init(format!(
            "{}: a reference FASTA is needed; pass --fasta",
            plugin.name()
        )));
    }
    if plugin.requires_hgvs() && !run.hgvs {
        return Err(PluginError::Init(format!(
            "{}: the HGVS notations are needed; pass --hgvs",
            plugin.name()
        )));
    }
    plugin.set_reference_fasta(reference_fasta);
    Ok(())
}

/// Attempt to create a built-in plugin by name.
///
/// Returns `None` if the name doesn't match any built-in plugin. `LoF`, the
/// name of LoFTEE's Ensembl VEP module (`LoF.pm`), creates the LoFTEE built-in,
/// so an Ensembl VEP argument string loads it unchanged.
fn create_builtin(name: &str) -> Option<Box<dyn BuiltinPlugin>> {
    match name {
        "CADD" => Some(Box::new(plugins::cadd::CaddPlugin::new())),
        "gnomADc" => Some(Box::new(plugins::gnomadc::GnomadcPlugin::new())),
        "REVEL" => Some(Box::new(plugins::revel::RevelPlugin::new())),
        "SpliceAI" => Some(Box::new(plugins::spliceai::SpliceAiPlugin::new())),
        "AlphaMissense" => Some(Box::new(plugins::alpha_missense::AlphaMissensePlugin::new())),
        "AncestralAllele" => Some(Box::new(
            plugins::ancestral_allele::AncestralAllelePlugin::new(),
        )),
        "dbNSFP" => Some(Box::new(plugins::dbnsfp::DbNsfpPlugin::new())),
        "dbscSNV" => Some(Box::new(plugins::dbscsnv::DbscSnvPlugin::new())),
        "DosageSensitivity" => Some(Box::new(
            plugins::dosage_sensitivity::DosageSensitivityPlugin::new(),
        )),
        "GWAS" => Some(Box::new(plugins::gwas::GwasPlugin::new())),
        "LOEUF" => Some(Box::new(plugins::loeuf::LoeufPlugin::new())),
        "LoFtool" => Some(Box::new(plugins::loftool::LoFtoolPlugin::new())),
        "LoFTEE" | "LoF" => Some(Box::new(plugins::loftee::LofteePlugin::new())),
        "MaxEntScan" => Some(Box::new(plugins::maxentscan::MaxEntScanPlugin::new())),
        "NMD" => Some(Box::new(plugins::nmd::NmdPlugin::new())),
        "OpenTargets" => Some(Box::new(plugins::open_targets::OpenTargetsPlugin::new())),
        "pLI" => Some(Box::new(plugins::pli::PliPlugin::new())),
        "PrimateAI" => Some(Box::new(plugins::primate_ai::PrimateAiPlugin::new())),
        "ClinPred" => Some(Box::new(plugins::clinpred::ClinPredPlugin::new())),
        "Conservation" => Some(Box::new(plugins::conservation::ConservationPlugin::new())),
        "BayesDel" => Some(Box::new(plugins::bayesdel::BayesDelPlugin::new())),
        "UTRAnnotator" => Some(Box::new(plugins::utr_annotator::UtrAnnotatorPlugin::new())),
        "AVI" => Some(Box::new(plugins::avi::AviPlugin::new())),
        "TranscriptAnnotator" => Some(Box::new(
            plugins::transcript_annotator::TranscriptAnnotatorPlugin::new(),
        )),
        "gnomADMt" => Some(Box::new(plugins::gnomadmt::GnomadMtPlugin::new())),
        "SpliceVault" => Some(Box::new(plugins::splicevault::SpliceVaultPlugin::new())),
        "Blosum62" => Some(Box::new(plugins::blosum62::Blosum62Plugin::new())),
        "Downstream" => Some(Box::new(plugins::downstream::DownstreamPlugin::new())),
        "HGVSIntronOffset" => Some(Box::new(
            plugins::hgvs_intron_offset::HgvsIntronOffsetPlugin::new(),
        )),
        "NearestExonJB" => Some(Box::new(
            plugins::nearest_exon_jb::NearestExonJbPlugin::new(),
        )),
        "TSSDistance" => Some(Box::new(plugins::tss_distance::TssDistancePlugin::new())),
        _ => None,
    }
}

/// The Ensembl VEP plugins vep-rs does not provide, each with the reason a
/// `--plugin` naming it is refused with: the modules of the release 116 plugin
/// collection that no built-in claims, the external `Gwava` of its web manifest,
/// and `PostGAP` and `StructuralVariantOverlap`, which release 116 removed from
/// the collection. A dylib of the same name, by path or under the plugins
/// directories, loads before the table is consulted. Every reason is one or two
/// sentences ending in a full stop; the refusal appends the pointer to the
/// support table.
pub const UNPORTED_ENSEMBL_PLUGINS: &[(&str, &str)] = &[
    (
        "AVADA",
        "Literature-evidence PubMed identifiers from one GRCh37-only file; absent from the release 116 web manifest.",
    ),
    (
        "CAPICE",
        "A GRCh37-only pathogenicity score from a single data release.",
    ),
    (
        "CCDSFilter",
        "An example filter plugin absent from the release 116 web manifest; offline, the cached transcripts carry no database cross-references for it to filter on, and the vep filter tool covers the use.",
    ),
    (
        "CSN",
        "Clinical Sequencing Nomenclature strings computed from the transcript model; rarely used, and absent from the release 116 web manifest.",
    ),
    (
        "Carol",
        "A combined SIFT and PolyPhen score; REVEL and AlphaMissense cover the use.",
    ),
    (
        "Condel",
        "A consensus of SIFT and PolyPhen scores that REVEL and AlphaMissense cover; the --condel flag is accepted without effect.",
    ),
    (
        "DAS",
        "An obsolete remote protocol queried per variant; absent from the release 116 web manifest.",
    ),
    (
        "DeNovo",
        "A filter plugin that writes per-run report files from a pedigree and the record's genotypes; a vep-rs plugin drops no output line and writes no file.",
    ),
    (
        "Draw",
        "Writes a PNG image per output row through the GD library.",
    ),
    (
        "EVE",
        "Missense-only evolutionary scores, popEVE included since release 116, distributed without stated data terms; AlphaMissense covers missense pathogenicity.",
    ),
    (
        "Enformer",
        "Over five thousand regulatory expression-impact features per variant from a lifted-over score file; rarely used.",
    ),
    (
        "FATHMM",
        "Runs an external Python script against a MySQL database.",
    ),
    ("FATHMM_MKL", "GRCh37-only coding and non-coding scores."),
    (
        "FlagLRG",
        "Flags Locus Reference Genomic transcripts; the LRG resource is retired and the --lrg flag is refused.",
    ),
    (
        "FunMotifs",
        "Tissue BED files served from a preprint's data server, GRCh37 only; rarely used.",
    ),
    (
        "G2P",
        "Diagnostic allelic-requirement filtering that aggregates per individual, reads VCF frequency collections and writes HTML and text reports.",
    ),
    (
        "GO",
        "The gene ontology file is built from the Ensembl core database on first use; vep-rs connects to no database.",
    ),
    (
        "GXA",
        "Queries the Gene Expression Atlas web service per gene; absent from the release 116 web manifest.",
    ),
    (
        "GeneBe",
        "Queries the GeneBe web service per variant under academic-only, rate-limited terms.",
    ),
    (
        "GeneSplicer",
        "Runs the external GeneSplicer binary per variant; MaxEntScan covers splice-site scoring.",
    ),
    (
        "Geno2MP",
        "Rare-variant genotype counts linked to HPO terms in GRCh37 coordinates; rarely used.",
    ),
    (
        "Gwava",
        "An external plugin outside the Ensembl VEP collection whose GRCh37 data source is retired.",
    ),
    (
        "HGVSReferenceBase",
        "The reference base of the long HGVS deletion or insertion form, computed from the transcript model; absent from the release 116 web manifest.",
    ),
    (
        "IntAct",
        "Molecular-interaction effects from the IntAct mutation tables across seven species; rarely used.",
    ),
    (
        "LD",
        "Linkage disequilibrium computed through the Ensembl variation database and genotype VCFs; vep-rs connects to no database.",
    ),
    (
        "LOVD",
        "Queries the LOVD web service per variant, which the plugin itself asks users not to do intensively.",
    ),
    (
        "LocalID",
        "Serves identifier input from a SQLite database built from the variation database; vep-rs reads identifier input from the variation cache with --format id.",
    ),
    (
        "MPC",
        "GRCh37-only regional missense constraint scores.",
    ),
    ("MTR", "A GRCh37-only missense tolerance ratio."),
    (
        "Mastermind",
        "Literature citation counts and identifiers from a vendor file behind registration, under commercial terms.",
    ),
    (
        "MaveDB",
        "Deep mutational scanning scores whose licence differs per dataset; rarely used.",
    ),
    (
        "MechPredict",
        "Gene-level probabilities of dominant-negative, gain-of-function and loss-of-function mechanisms, assembled from three OSF datasets; rarely used.",
    ),
    (
        "NearestGene",
        "Reads the Ensembl core database by default; the --nearest flag covers the nearest-gene lookup from the cache.",
    ),
    (
        "NonSynonymousFilter",
        "An example filter plugin absent from the release 116 web manifest; a vep-rs plugin drops no output line, and the vep filter tool covers the use.",
    ),
    (
        "PON_P2",
        "Queries its web service per variant, or reads a pre-computed GRCh37 file derived from COSMIC under its restricted terms.",
    ),
    (
        "Paralogues",
        "Variants at paralogous aligned positions, read from the Compara database, a web service or a GRCh38 matches file; rarely used.",
    ),
    (
        "PhenotypeOrthologous",
        "Mouse and rat orthologue phenotypes from a GRCh38 GFF3; rarely used.",
    ),
    (
        "Phenotypes",
        "The phenotype file is generated from the Ensembl variation database on first use; vep-rs connects to no database.",
    ),
    (
        "PolyPhen_SIFT",
        "SIFT and PolyPhen predictions from a SQLite database built from the Ensembl database; the vep-rs cache carries both natively.",
    ),
    (
        "PostGAP",
        "Removed from the Ensembl VEP plugin collection in release 116.",
    ),
    (
        "PromoterAI",
        "Promoter expression-impact scores under a commercial licence agreement.",
    ),
    (
        "ProtVar",
        "Structure stability, pocket and interface context from a SQLite database; vep-rs embeds no SQLite.",
    ),
    (
        "ProteinSeqs",
        "Writes two FASTA files of reference and mutated protein sequences as a side effect.",
    ),
    (
        "RankFilter",
        "An example filter plugin absent from the release 116 web manifest; --most_severe, --pick and the vep filter tool cover the use.",
    ),
    (
        "RefSeqHGVS",
        "Reads the Ensembl otherfeatures database; --refseq with --hgvs yields RefSeq HGVS notation directly, and the plugin is absent from the release 116 web manifest.",
    ),
    (
        "ReferenceQuality",
        "Reference-assembly issue annotations from GRC GFF3 files merged by the user; rarely used.",
    ),
    (
        "RiboseqORFs",
        "Recomputes consequences over Ribo-seq open reading frame models; absent from the release 116 web manifest.",
    ),
    (
        "SameCodon",
        "Reads the Ensembl variation database for other variants in the same codon.",
    ),
    (
        "SingleLetterAA",
        "Rewrites HGVSp with one-letter amino-acid codes, a form of the protein notation vep-rs does not compute.",
    ),
    (
        "SpliceRegion",
        "Its splice_donor_5th_base_variant, splice_donor_region_variant and splice_polypyrimidine_tract_variant terms are core consequence terms of both engines; only its Extended mode adds a key of its own.",
    ),
    (
        "StructuralVariantOverlap",
        "Removed from the Ensembl VEP plugin collection in release 116; a --custom structural-variant VCF source covers the use.",
    ),
    (
        "SubsetVCF",
        "All alternate alleles and INFO fields of the overlapping records of a VCF; a --custom VCF source covers the common case, and the plugin is absent from the release 116 web manifest.",
    ),
    ("VARITY", "A GRCh37-only rare-missense predictor."),
    (
        "mutfunc",
        "Protein destabilisation and motif effects from a SQLite database; vep-rs embeds no SQLite.",
    ),
    (
        "neXtProt",
        "Queries the neXtProt SPARQL endpoint per variant and needs a database connection; absent from the release 116 web manifest.",
    ),
    (
        "satMutMPRA",
        "Saturation mutagenesis MPRA effects for twenty-one regulatory elements; rarely used.",
    ),
];

/// The reason `name` is refused with when it is an Ensembl VEP plugin vep-rs
/// does not provide ([`UNPORTED_ENSEMBL_PLUGINS`]); `None` for a built-in or
/// a name the collection does not carry.
pub fn unported_reason(name: &str) -> Option<&'static str> {
    UNPORTED_ENSEMBL_PLUGINS
        .iter()
        .find(|(plugin, _)| *plugin == name)
        .map(|(_, reason)| *reason)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::OutputFormat;
    use indexmap::IndexMap;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, Condvar, Mutex};
    use std::time::Duration;
    use vep_core::consequence::{FeatureType, TranscriptConsequence};

    #[test]
    fn new_registry_is_empty() {
        let registry = BuiltinRegistry::new();
        assert!(registry.is_empty());
        assert_eq!(registry.plugin_count(), 0);
    }

    #[test]
    fn unknown_plugin_is_unresolved() {
        let mut registry = BuiltinRegistry::new();
        let unresolved = registry
            .load_from_args(&["UnknownPlugin,foo".to_string()])
            .unwrap();
        assert_eq!(unresolved, vec!["UnknownPlugin,foo"]);
        assert!(registry.is_empty());
    }

    /// A built-in whose `init` fails is left out of the registry and named in the
    /// outcome, and the arguments after it are still read; `load_from_args`
    /// returns that failure as its error.
    #[test]
    fn failed_init_is_reported_and_the_rest_still_load() {
        let mut registry = BuiltinRegistry::new();
        let outcome = registry.load_from_args_with_context(
            &[
                "LoFtool,/no/such/scores.txt".to_string(),
                "UnknownPlugin".to_string(),
            ],
            None,
            None,
        );
        assert!(registry.is_empty());
        assert_eq!(outcome.unresolved, ["UnknownPlugin"]);
        let [(name, error)] = &outcome.failed[..] else {
            panic!("one failure, got {:?}", outcome.failed);
        };
        assert_eq!(name, "LoFtool");
        assert_eq!(
            error.to_string(),
            "plugin initialization failed: LoFtool: scores file /no/such/scores.txt not found"
        );
        assert!(registry
            .load_from_args(&["LoFtool,/no/such/scores.txt".to_string()])
            .is_err());
    }

    /// `LoF` is the name of LoFTEE's Ensembl VEP module, so a `--plugin LoF,...`
    /// argument written for Ensembl VEP creates the built-in named `LoFTEE`.
    #[test]
    fn create_builtin_accepts_lof_as_loftee() {
        let plugin = create_builtin("LoF").expect("LoF names the LoFTEE built-in");
        assert_eq!(plugin.name(), "LoFTEE");
        assert_eq!(create_builtin("LoFTEE").unwrap().name(), "LoFTEE");
        assert!(create_builtin("lof").is_none());
    }

    /// Every plugin name of Ensembl VEP's release 116: the 83 modules of the
    /// plugin collection, the external `Gwava` of its web manifest, `PostGAP` and
    /// `StructuralVariantOverlap`, which release 116 removed from the collection,
    /// and LoFTEE's `LoF`.
    const ENSEMBL_PLUGIN_NAMES: &[&str] = &[
        "AlphaMissense",
        "AncestralAllele",
        "AVADA",
        "AVI",
        "BayesDel",
        "Blosum62",
        "CADD",
        "CAPICE",
        "Carol",
        "CCDSFilter",
        "ClinPred",
        "Condel",
        "Conservation",
        "CSN",
        "DAS",
        "dbNSFP",
        "dbscSNV",
        "DeNovo",
        "DosageSensitivity",
        "Downstream",
        "Draw",
        "Enformer",
        "EVE",
        "FATHMM",
        "FATHMM_MKL",
        "FlagLRG",
        "FunMotifs",
        "G2P",
        "GeneBe",
        "GeneSplicer",
        "Geno2MP",
        "gnomADc",
        "gnomADMt",
        "GO",
        "GWAS",
        "Gwava",
        "GXA",
        "HGVSIntronOffset",
        "HGVSReferenceBase",
        "IntAct",
        "LD",
        "LocalID",
        "LOEUF",
        "LoF",
        "LoFtool",
        "LOVD",
        "Mastermind",
        "MaveDB",
        "MaxEntScan",
        "MechPredict",
        "MPC",
        "MTR",
        "mutfunc",
        "NearestExonJB",
        "NearestGene",
        "neXtProt",
        "NMD",
        "NonSynonymousFilter",
        "OpenTargets",
        "Paralogues",
        "PhenotypeOrthologous",
        "Phenotypes",
        "pLI",
        "PolyPhen_SIFT",
        "PON_P2",
        "PostGAP",
        "PrimateAI",
        "PromoterAI",
        "ProteinSeqs",
        "ProtVar",
        "RankFilter",
        "ReferenceQuality",
        "RefSeqHGVS",
        "REVEL",
        "RiboseqORFs",
        "SameCodon",
        "satMutMPRA",
        "SingleLetterAA",
        "SpliceAI",
        "SpliceRegion",
        "SpliceVault",
        "StructuralVariantOverlap",
        "SubsetVCF",
        "TranscriptAnnotator",
        "TSSDistance",
        "UTRAnnotator",
        "VARITY",
    ];

    /// Some built-ins that compute their fields from the transcript model and
    /// take no data file; the unported table refuses none of them.
    const COMPUTED_FROM_THE_TRANSCRIPT_MODEL: &[&str] = &[
        "Blosum62",
        "Downstream",
        "HGVSIntronOffset",
        "NearestExonJB",
        "TSSDistance",
    ];

    /// The built-in names and the unported table partition Ensembl VEP's plugin
    /// names: no name is both, every unported entry is one of the 87, and every
    /// one of the 87 is a built-in or an unported entry, the listed computed
    /// plugins among the built-ins.
    #[test]
    fn unported_and_builtin_names_partition_the_collection() {
        assert_eq!(ENSEMBL_PLUGIN_NAMES.len(), 87);
        for (name, reason) in UNPORTED_ENSEMBL_PLUGINS {
            assert!(
                ENSEMBL_PLUGIN_NAMES.contains(name),
                "{name} is not an Ensembl VEP plugin name"
            );
            assert!(
                create_builtin(name).is_none(),
                "{name} is a built-in and an unported entry"
            );
            assert!(
                !COMPUTED_FROM_THE_TRANSCRIPT_MODEL.contains(name),
                "{name} is computed from the transcript model and refused"
            );
            assert!(
                reason.ends_with('.') && reason.matches(". ").count() <= 1,
                "{name}: the reason is one or two sentences ending in a full stop: {reason}"
            );
            assert!(
                reason
                    .split(|c: char| !c.is_ascii_digit())
                    .filter(|digits| digits.len() >= 3)
                    .all(|digits| digits == "116"),
                "{name}: release 116 is the only version number a reason names: {reason}"
            );
            assert_eq!(unported_reason(name), Some(*reason));
        }
        let mut names: Vec<&str> = UNPORTED_ENSEMBL_PLUGINS
            .iter()
            .map(|(name, _)| *name)
            .collect();
        names.sort_unstable();
        names.dedup();
        assert_eq!(
            names.len(),
            UNPORTED_ENSEMBL_PLUGINS.len(),
            "duplicate entry"
        );

        for name in ENSEMBL_PLUGIN_NAMES {
            let claimed = create_builtin(name).is_some()
                || unported_reason(name).is_some()
                || COMPUTED_FROM_THE_TRANSCRIPT_MODEL.contains(name);
            assert!(
                claimed,
                "{name} is neither a built-in nor an unported entry"
            );
        }
        assert_eq!(unported_reason("REVEL"), None);
        assert_eq!(unported_reason("NoSuchPlugin"), None);
    }

    /// A built-in whose data file is a URL is refused before `init`, as the
    /// instantiation failure naming the file, with or without a `file=` key.
    #[test]
    fn a_remote_data_file_is_an_instantiation_failure() {
        let mut registry = BuiltinRegistry::new();
        let outcome = registry.load_from_args_with_context(
            &[
                "REVEL,file=https://example.invalid/x.tsv.gz".to_string(),
                "CADD,ftp://example.invalid/synthetic_scores.tsv.gz".to_string(),
            ],
            None,
            None,
        );
        assert!(registry.is_empty());
        assert!(outcome.unresolved.is_empty());
        let failures: Vec<(String, String)> = outcome
            .failed
            .iter()
            .map(|(name, error)| (name.clone(), error.to_string()))
            .collect();
        assert_eq!(
            failures,
            [
                (
                    "REVEL".to_string(),
                    "plugin initialization failed: remote data files are not supported; download https://example.invalid/x.tsv.gz and index it locally".to_string()
                ),
                (
                    "CADD".to_string(),
                    "plugin initialization failed: remote data files are not supported; download ftp://example.invalid/synthetic_scores.tsv.gz and index it locally".to_string()
                ),
            ]
        );
        let revel = create_builtin("REVEL").unwrap();
        assert_eq!(revel.remote_data_file("a=b"), None);
        assert_eq!(revel.remote_data_file("x"), None);
    }

    /// The refusal reads a word under the plugin's own grammar. LoFTEE's
    /// `LoF.pm` keeps a word only when `split /:/` yields exactly two fields
    /// (lines 68-69), and a URL's `://` makes three or more, so no LoFTEE word
    /// names a remote data file: `gerp_bigwig:https://...` is dropped, logged,
    /// and the plugin loads without it, as Ensembl VEP's does; `gerp_bigwig=` with
    /// a URL is the two-field word with the unknown key `gerp_bigwig=https`,
    /// stored and never read there, ignored here. Under the `key=value` grammar
    /// a word without `=` is a positional path, so REVEL's `file:https://...`
    /// stays the local path it names; a positional URL is named whole, query
    /// and all, and a `key=` URL is cut at the first `=` of its query, as
    /// `params_to_hash` cuts the value it hands the plugin.
    #[test]
    fn a_remote_data_file_is_read_under_the_plugins_grammar() {
        const URL: &str = "https://example.invalid/gerp.bw";
        let corpus = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh38-plugin-loftee/plugin_data");
        let loftee_path = format!("loftee_path:{}/", corpus.display());
        let bigwig = corpus
            .join("gerp_conservation_scores.homo_sapiens.GRCh38.bw")
            .display()
            .to_string();
        let load = |arg: String| {
            let mut registry = BuiltinRegistry::new();
            let outcome = registry.load_from_args_with_context(
                &[arg],
                None,
                Some(vep_core::assembly::Assembly::Grch38),
            );
            let failures: Vec<(String, String)> = outcome
                .failed
                .iter()
                .map(|(name, error)| (name.clone(), error.to_string()))
                .collect();
            (registry.plugin_names(), failures)
        };

        let query = format!("{URL}?to=https://other.invalid/g.bw");
        for word in [
            format!("gerp_bigwig:{URL}"),
            format!("gerp_bigwig={URL}"),
            format!("gerp_bigwig:{bigwig}"),
            URL.to_string(),
        ] {
            assert_eq!(
                load(format!("LoF,{loftee_path},{word}")),
                (vec!["LoFTEE".to_string()], vec![]),
                "{word}"
            );
        }

        let loftee = create_builtin("LoF").unwrap();
        for word in [
            format!("gerp_bigwig:{URL}"),
            format!("gerp_bigwig:{URL}::"),
            format!("gerp_bigwig:{query}"),
            format!("gerp_bigwig={URL}"),
            URL.to_string(),
            query.clone(),
            format!("gerp_bigwig:{bigwig}"),
            format!("gerp_bigwig:/data/x={URL}"),
            "min_intron_size:20:".to_string(),
        ] {
            assert_eq!(loftee.remote_data_file(&word), None, "{word}");
        }
        let revel = create_builtin("REVEL").unwrap();
        assert_eq!(revel.remote_data_file(&format!("file:{URL}")), None);
        assert_eq!(revel.remote_data_file(&format!("file={URL}")), Some(URL));
        assert_eq!(revel.remote_data_file(&query), Some(query.as_str()));
        let cut = format!("{URL}?to");
        assert_eq!(
            revel.remote_data_file(&format!("file={query}")),
            Some(cut.as_str())
        );
        assert_eq!(
            load(format!("REVEL,file={query}")),
            (
                vec![],
                vec![(
                    "REVEL".to_string(),
                    format!("plugin initialization failed: remote data files are not supported; download {cut} and index it locally")
                )]
            )
        );
    }

    /// A plugin that records the order of the context calls around `init`.
    struct ContextOrderPlugin {
        calls: Vec<String>,
    }

    impl BuiltinPlugin for ContextOrderPlugin {
        fn name(&self) -> &str {
            "ContextOrder"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            self.calls.push("init".to_string());
            Ok(())
        }
        fn set_assembly(&mut self, assembly: Option<vep_core::assembly::Assembly>) {
            self.calls.push(format!("assembly {assembly:?}"));
        }
        fn set_output_format(&mut self, format: OutputFormat) {
            self.calls.push(format!("format {format:?}"));
        }
        fn set_transcript_set(&mut self, refseq: bool, merged: bool) {
            self.calls.push(format!("set {refseq} {merged}"));
        }
        fn set_plugin_dirs(&mut self, dirs: &[String]) {
            self.calls.push(format!("dirs {dirs:?}"));
        }
        fn set_reference_fasta(&mut self, _fasta: Option<Arc<vep_fasta::IndexedFasta>>) {
            self.calls.push("fasta".to_string());
        }
    }

    /// Ensembl VEP's plugins read the output format, the transcript set and the
    /// plugins directory from `$self->{config}` in their constructors, so the
    /// registry hands all three to a plugin before `init`, after the assembly;
    /// the FASTA still follows `init`.
    #[test]
    fn output_format_and_transcript_set_reach_a_plugin_before_init() {
        let mut plugin = ContextOrderPlugin { calls: Vec::new() };
        let dirs = ["/plugins".to_string()];
        let run = RunContext {
            output_format: OutputFormat::Json,
            refseq: true,
            merged: false,
            hgvs: false,
            plugin_dirs: &dirs,
        };
        load_builtin(
            &mut plugin,
            &[],
            None,
            Some(vep_core::assembly::Assembly::Grch38),
            run,
        )
        .unwrap();
        assert_eq!(
            plugin.calls,
            [
                "assembly Some(Grch38)",
                "format Json",
                "set true false",
                "dirs [\"/plugins\"]",
                "init",
                "fasta"
            ]
        );
    }

    /// A built-in that reads the HGVSc coordinates loads only when the run
    /// computes the HGVS notations; otherwise it is left out and named in the
    /// outcome with the flag to pass.
    #[test]
    fn a_plugin_needing_the_hgvs_notations_is_refused_without_them() {
        let args = ["HGVSIntronOffset".to_string()];
        let mut registry = BuiltinRegistry::new();
        let outcome = registry.load_from_args_for_run(&args, None, None, RunContext::default());
        assert!(registry.is_empty());
        let [(name, error)] = &outcome.failed[..] else {
            panic!("one failure, got {:?}", outcome.failed);
        };
        assert_eq!(name, "HGVSIntronOffset");
        assert_eq!(
            error.to_string(),
            "plugin initialization failed: HGVSIntronOffset: the HGVS notations are needed; \
             pass --hgvs"
        );
        let run = RunContext {
            hgvs: true,
            ..RunContext::default()
        };
        let outcome = registry.load_from_args_for_run(&args, None, None, run);
        assert!(outcome.failed.is_empty() && outcome.unresolved.is_empty());
        assert_eq!(registry.plugin_names(), ["HGVSIntronOffset"]);
    }

    /// A mock plugin that supports the two-phase prefetch/annotate protocol.
    /// Tracks call counts to verify the parallel execution path is used.
    struct MockPrefetchPlugin {
        name: String,
        prefetch_count: Arc<AtomicUsize>,
        annotate_count: Arc<AtomicUsize>,
    }

    impl BuiltinPlugin for MockPrefetchPlugin {
        fn name(&self) -> &str {
            &self.name
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn prefetch(
            &self,
            _variants: &[InputVariant],
        ) -> Result<Option<PrefetchData>, PluginError> {
            self.prefetch_count.fetch_add(1, Ordering::Relaxed);
            Ok(Some(PrefetchData::new(42u64)))
        }
        fn annotate(
            &self,
            variants: &mut [InputVariant],
            data: PrefetchData,
        ) -> Result<(), PluginError> {
            let marker: u64 = data.downcast()?;
            assert_eq!(marker, 42);
            self.annotate_count.fetch_add(1, Ordering::Relaxed);
            for v in variants.iter_mut() {
                v.plugin_data
                    .insert(format!("{}_ran", self.name), "true".to_string());
            }
            Ok(())
        }
        fn run(
            &self,
            _consequence: &TranscriptConsequence,
            _variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            Ok(IndexMap::new())
        }
    }

    /// A mock plugin that does not support prefetch (returns None),
    /// forcing the registry to fall back to run_batch.
    struct MockFallbackPlugin {
        run_batch_count: Arc<AtomicUsize>,
    }

    impl BuiltinPlugin for MockFallbackPlugin {
        fn name(&self) -> &str {
            "Fallback"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
            self.run_batch_count.fetch_add(1, Ordering::Relaxed);
            for v in variants.iter_mut() {
                v.plugin_data
                    .insert("fallback_ran".to_string(), "true".to_string());
            }
            Ok(())
        }
        fn run(
            &self,
            _consequence: &TranscriptConsequence,
            _variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            Ok(IndexMap::new())
        }
    }

    #[test]
    fn two_phase_execution_calls_prefetch_and_annotate() {
        let prefetch_count = Arc::new(AtomicUsize::new(0));
        let annotate_count = Arc::new(AtomicUsize::new(0));

        let plugin = MockPrefetchPlugin {
            name: "TestPlugin".to_string(),
            prefetch_count: Arc::clone(&prefetch_count),
            annotate_count: Arc::clone(&annotate_count),
        };

        let mut registry = BuiltinRegistry::new();
        registry.plugins.push(Box::new(plugin));

        let mut variants = vec![InputVariant::new(
            "1".into(),
            100,
            100,
            b"A".to_vec(),
            b"G".to_vec(),
        )];

        assert!(registry.run_batch(&mut variants).failures.is_empty());

        assert_eq!(prefetch_count.load(Ordering::Relaxed), 1);
        assert_eq!(annotate_count.load(Ordering::Relaxed), 1);
        assert_eq!(
            variants[0].plugin_data.get("TestPlugin_ran").unwrap(),
            "true"
        );
    }

    #[test]
    fn fallback_to_run_batch_when_prefetch_returns_none() {
        let run_batch_count = Arc::new(AtomicUsize::new(0));

        let plugin = MockFallbackPlugin {
            run_batch_count: Arc::clone(&run_batch_count),
        };

        let mut registry = BuiltinRegistry::new();
        registry.plugins.push(Box::new(plugin));

        let mut variants = vec![InputVariant::new(
            "1".into(),
            100,
            100,
            b"A".to_vec(),
            b"G".to_vec(),
        )];

        assert!(registry.run_batch(&mut variants).failures.is_empty());

        assert_eq!(run_batch_count.load(Ordering::Relaxed), 1);
        assert_eq!(variants[0].plugin_data.get("fallback_ran").unwrap(), "true");
    }

    #[test]
    fn mixed_prefetch_and_fallback_plugins() {
        let prefetch_count_a = Arc::new(AtomicUsize::new(0));
        let annotate_count_a = Arc::new(AtomicUsize::new(0));
        let prefetch_count_b = Arc::new(AtomicUsize::new(0));
        let annotate_count_b = Arc::new(AtomicUsize::new(0));
        let run_batch_count = Arc::new(AtomicUsize::new(0));

        let plugin_a = MockPrefetchPlugin {
            name: "PluginA".to_string(),
            prefetch_count: Arc::clone(&prefetch_count_a),
            annotate_count: Arc::clone(&annotate_count_a),
        };
        let plugin_b = MockPrefetchPlugin {
            name: "PluginB".to_string(),
            prefetch_count: Arc::clone(&prefetch_count_b),
            annotate_count: Arc::clone(&annotate_count_b),
        };
        let fallback = MockFallbackPlugin {
            run_batch_count: Arc::clone(&run_batch_count),
        };

        let mut registry = BuiltinRegistry::new();
        registry.plugins.push(Box::new(plugin_a));
        registry.plugins.push(Box::new(fallback));
        registry.plugins.push(Box::new(plugin_b));

        let mut variants = vec![InputVariant::new(
            "1".into(),
            100,
            100,
            b"A".to_vec(),
            b"G".to_vec(),
        )];

        assert!(registry.run_batch(&mut variants).failures.is_empty());

        assert_eq!(prefetch_count_a.load(Ordering::Relaxed), 1);
        assert_eq!(annotate_count_a.load(Ordering::Relaxed), 1);
        assert_eq!(prefetch_count_b.load(Ordering::Relaxed), 1);
        assert_eq!(annotate_count_b.load(Ordering::Relaxed), 1);
        assert_eq!(run_batch_count.load(Ordering::Relaxed), 1);

        assert_eq!(variants[0].plugin_data.get("PluginA_ran").unwrap(), "true");
        assert_eq!(variants[0].plugin_data.get("PluginB_ran").unwrap(), "true");
        assert_eq!(variants[0].plugin_data.get("fallback_ran").unwrap(), "true");
    }

    /// `count` substitutions at distinct positions of one chromosome.
    fn batch_of(count: u64) -> Vec<InputVariant> {
        (0..count)
            .map(|i| {
                let pos = 1_000 + 100 * i;
                InputVariant::new("21".into(), pos, pos, b"A".to_vec(), b"G".to_vec())
            })
            .collect()
    }

    /// Ensembl VEP's tabix plugins read their data once per variant: every
    /// `run` call asks `get_data` for the variant's own coordinates, and the
    /// base class answers from a results cache keyed on those coordinates or a
    /// region cache of one region per chromosome, which a miss refills with a
    /// tabix query over the variant's span widened by `expand_right`
    /// (`BaseVepTabixPlugin.pm:112-114, 350-449`). The registry reads once per
    /// batch instead: a plugin's `prefetch` and `annotate` each run exactly
    /// once for a batch of 200 variants and once more per further batch, on
    /// every batch entry point, and the batch is annotated whole.
    #[test]
    fn prefetch_runs_once_per_batch_not_per_variant() {
        let prefetches = Arc::new(AtomicUsize::new(0));
        let annotates = Arc::new(AtomicUsize::new(0));
        let mut registry = BuiltinRegistry::new();
        registry.plugins.push(Box::new(MockPrefetchPlugin {
            name: "Once".to_string(),
            prefetch_count: Arc::clone(&prefetches),
            annotate_count: Arc::clone(&annotates),
        }));

        let mut variants = batch_of(200);
        assert!(registry.run_batch(&mut variants).failures.is_empty());
        assert_eq!(prefetches.load(Ordering::Relaxed), 1);
        assert_eq!(annotates.load(Ordering::Relaxed), 1);
        assert!(variants
            .iter()
            .all(|v| v.plugin_data.get("Once_ran").map(String::as_str) == Some("true")));

        let mut timed = batch_of(200);
        assert!(registry.run_batch_timed(&mut timed).1.failures.is_empty());
        let mut sequential = batch_of(200);
        assert!(registry
            .run_batch_sequential(&mut sequential)
            .1
            .failures
            .is_empty());
        let mut ahead = batch_of(200);
        let prefetch = registry.prefetch_ahead(&ahead);
        assert_eq!(prefetches.load(Ordering::Relaxed), 4);
        assert!(registry
            .annotate_batch(&mut ahead, prefetch)
            .failures
            .is_empty());
        assert_eq!(prefetches.load(Ordering::Relaxed), 4);
        assert_eq!(annotates.load(Ordering::Relaxed), 4);
        for batch in [&timed, &sequential, &ahead] {
            assert_eq!(batch.len(), 200);
            assert!(batch.iter().all(|v| v.plugin_data.contains_key("Once_ran")));
        }
    }

    /// What the plugins of one registry see of each other's prefetch: how
    /// many prefetches are in flight, the order the prefetches finished in
    /// and the order the annotate phase ran the plugins in.
    #[derive(Default)]
    struct PrefetchLog {
        in_flight: usize,
        most_in_flight: usize,
        prefetched: Vec<&'static str>,
        annotated: Vec<&'static str>,
    }

    /// A rendezvous of every plugin's prefetch: `wait_until` blocks on the
    /// shared log until a condition holds or `WAIT` passes, so a registry
    /// that ran the prefetches one after another fails the test instead of
    /// hanging it.
    #[derive(Default)]
    struct Rendezvous {
        log: Mutex<PrefetchLog>,
        changed: Condvar,
    }

    const WAIT: Duration = Duration::from_secs(30);

    impl Rendezvous {
        fn wait_until(&self, holds: impl Fn(&PrefetchLog) -> bool) -> bool {
            let guard = self.log.lock().unwrap();
            let (_guard, wait) = self
                .changed
                .wait_timeout_while(guard, WAIT, |log| !holds(log))
                .unwrap();
            !wait.timed_out()
        }

        fn update(&self, change: impl FnOnce(&mut PrefetchLog)) {
            change(&mut self.log.lock().unwrap());
            self.changed.notify_all();
        }
    }

    /// Plugin number `index` of `total` sharing one [`Rendezvous`]: its
    /// prefetch waits until every plugin's prefetch is in flight, then until
    /// the plugins after it in load order have finished theirs, so the
    /// prefetches finish in reverse load order; its annotate writes `name`.
    struct RendezvousPlugin {
        name: &'static str,
        index: usize,
        total: usize,
        rendezvous: Arc<Rendezvous>,
    }

    impl BuiltinPlugin for RendezvousPlugin {
        fn name(&self) -> &str {
            self.name
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn prefetch(
            &self,
            _variants: &[InputVariant],
        ) -> Result<Option<PrefetchData>, PluginError> {
            self.rendezvous.update(|log| {
                log.in_flight += 1;
                log.most_in_flight = log.most_in_flight.max(log.in_flight);
            });
            // The high-water mark, not the live count: a plugin that finished
            // its prefetch has decremented `in_flight`, and a waiter checking
            // the live count after that blocks until the timeout.
            if !self
                .rendezvous
                .wait_until(|log| log.most_in_flight == self.total)
            {
                return Err(PluginError::Run(format!(
                    "{}: the other prefetches never ran beside this one",
                    self.name
                )));
            }
            let later = self.total - 1 - self.index;
            if !self
                .rendezvous
                .wait_until(|log| log.prefetched.len() == later)
            {
                return Err(PluginError::Run(format!(
                    "{}: the plugins after it never finished their prefetch",
                    self.name
                )));
            }
            self.rendezvous.update(|log| {
                log.prefetched.push(self.name);
                log.in_flight -= 1;
            });
            Ok(Some(PrefetchData::new(self.index)))
        }
        fn annotate(
            &self,
            variants: &mut [InputVariant],
            data: PrefetchData,
        ) -> Result<(), PluginError> {
            assert_eq!(data.downcast::<usize>()?, self.index);
            self.rendezvous.update(|log| log.annotated.push(self.name));
            for v in variants.iter_mut() {
                v.plugin_data
                    .insert(self.name.to_string(), self.index.to_string());
            }
            Ok(())
        }
    }

    /// The prefetch phase runs every plugin at once on the rayon pool, so the
    /// three prefetches are in flight together and finish in the order their
    /// data arrives (here the reverse of load order); the annotate phase still
    /// runs the plugins in load order, the order Ensembl VEP's `run_plugins`
    /// walks `@{$self->{plugins}}` (`OutputFactory.pm:2212-2256`), so the
    /// fields land in load order whatever the prefetches' timing. Each
    /// prefetch holds its thread until the others are in flight, so the test
    /// runs on its own pool of one thread per plugin; a pool with fewer
    /// threads than plugins never has every prefetch in flight at once.
    #[test]
    fn plugins_prefetch_in_parallel_and_annotate_in_registry_order() {
        const NAMES: [&str; 3] = ["First", "Second", "Third"];
        let rendezvous = Arc::new(Rendezvous::default());
        let mut registry = BuiltinRegistry::new();
        for (index, name) in NAMES.iter().enumerate() {
            registry.plugins.push(Box::new(RendezvousPlugin {
                name,
                index,
                total: NAMES.len(),
                rendezvous: Arc::clone(&rendezvous),
            }));
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(NAMES.len())
            .build()
            .unwrap();

        let mut variants = batch_of(5);
        let outcome = pool.install(|| registry.run_batch(&mut variants));
        assert!(outcome.failures.is_empty(), "{:?}", outcome.failures);

        let log = rendezvous.log.lock().unwrap();
        assert_eq!(log.most_in_flight, NAMES.len());
        assert_eq!(log.prefetched, ["Third", "Second", "First"]);
        assert_eq!(log.annotated, NAMES);
        for variant in &variants {
            assert_eq!(variant.plugin_data.keys().collect::<Vec<_>>(), NAMES);
        }
    }

    /// A two-phase plugin that prefetches data for a batch holding a record on
    /// chromosome 21 and, for every other batch, opts out of the two-phase
    /// protocol as the tabix built-ins do for a batch with no region to query;
    /// its `run_batch` writes a different field from its `annotate`.
    struct OptOutPlugin {
        prefetches: Arc<AtomicUsize>,
        annotates: Arc<AtomicUsize>,
        run_batches: Arc<AtomicUsize>,
    }

    impl BuiltinPlugin for OptOutPlugin {
        fn name(&self) -> &str {
            "OptOut"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
            self.prefetches.fetch_add(1, Ordering::Relaxed);
            if variants.iter().any(|v| v.chr == "21") {
                Ok(Some(PrefetchData::new(())))
            } else {
                Ok(None)
            }
        }
        fn annotate(
            &self,
            variants: &mut [InputVariant],
            _data: PrefetchData,
        ) -> Result<(), PluginError> {
            self.annotates.fetch_add(1, Ordering::Relaxed);
            for v in variants.iter_mut() {
                v.plugin_data
                    .insert("annotated".to_string(), "1".to_string());
            }
            Ok(())
        }
        fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
            self.run_batches.fetch_add(1, Ordering::Relaxed);
            for v in variants.iter_mut() {
                v.plugin_data
                    .insert("run_batch".to_string(), "1".to_string());
            }
            Ok(())
        }
    }

    /// `Ok(None)` from a plugin's prefetch is the plugin's choice, per batch,
    /// not to use the two-phase protocol: the registry calls its `run_batch`
    /// for that batch and never its `annotate`, and uses `annotate` again on
    /// the next batch whose prefetch returns data, through each batch entry
    /// point, the prefetch run ahead of annotation included (the `None`
    /// travels in the plugin's slot).
    #[test]
    fn a_plugin_returning_none_from_prefetch_falls_back_to_run_batch() {
        for way in 0..4 {
            let prefetches = Arc::new(AtomicUsize::new(0));
            let annotates = Arc::new(AtomicUsize::new(0));
            let run_batches = Arc::new(AtomicUsize::new(0));
            let mut registry = BuiltinRegistry::new();
            registry.plugins.push(Box::new(OptOutPlugin {
                prefetches: Arc::clone(&prefetches),
                annotates: Arc::clone(&annotates),
                run_batches: Arc::clone(&run_batches),
            }));
            let run = |variants: &mut [InputVariant]| match way {
                0 => registry.run_batch(variants),
                1 => registry.run_batch_timed(variants).1,
                2 => registry.run_batch_sequential(variants).1,
                _ => {
                    let ahead = registry.prefetch_ahead(variants);
                    registry.annotate_batch(variants, ahead)
                }
            };

            let mut elsewhere = vec![
                InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec()),
                InputVariant::new("X".into(), 200, 200, b"C".to_vec(), b"T".to_vec()),
            ];
            assert!(run(&mut elsewhere).failures.is_empty(), "way {way}");
            assert_eq!(prefetches.load(Ordering::Relaxed), 1, "way {way}");
            assert_eq!(run_batches.load(Ordering::Relaxed), 1, "way {way}");
            assert_eq!(annotates.load(Ordering::Relaxed), 0, "way {way}");
            for variant in &elsewhere {
                assert_eq!(
                    variant.plugin_data.keys().collect::<Vec<_>>(),
                    ["run_batch"],
                    "way {way}"
                );
            }

            let mut on_21 = batch_of(3);
            assert!(run(&mut on_21).failures.is_empty(), "way {way}");
            assert_eq!(prefetches.load(Ordering::Relaxed), 2, "way {way}");
            assert_eq!(run_batches.load(Ordering::Relaxed), 1, "way {way}");
            assert_eq!(annotates.load(Ordering::Relaxed), 1, "way {way}");
            for variant in &on_21 {
                assert_eq!(
                    variant.plugin_data.keys().collect::<Vec<_>>(),
                    ["annotated"],
                    "way {way}"
                );
            }
        }
    }

    /// A two-phase plugin whose prefetch reads the consequence rows, so the
    /// registry prefetches it once the batch is annotated.
    struct MockRowReadingPlugin {
        inner: MockPrefetchPlugin,
    }

    impl BuiltinPlugin for MockRowReadingPlugin {
        fn name(&self) -> &str {
            self.inner.name()
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn prefetch_reads_consequences(&self) -> bool {
            true
        }
        fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
            self.inner.prefetch(variants)
        }
        fn annotate(
            &self,
            variants: &mut [InputVariant],
            data: PrefetchData,
        ) -> Result<(), PluginError> {
            self.inner.annotate(variants, data)
        }
    }

    /// The prefetch run ahead of annotation covers the plugins whose prefetch
    /// reads no consequence row and leaves the others, whose slots the annotate
    /// phase fills before applying every plugin in load order; the batch ends
    /// as `run_batch` leaves it.
    #[test]
    fn prefetch_ahead_leaves_row_reading_plugins_to_the_annotate_phase() {
        let ahead_prefetches = Arc::new(AtomicUsize::new(0));
        let ahead_annotates = Arc::new(AtomicUsize::new(0));
        let later_prefetches = Arc::new(AtomicUsize::new(0));
        let later_annotates = Arc::new(AtomicUsize::new(0));
        let fallbacks = Arc::new(AtomicUsize::new(0));
        let mut registry = BuiltinRegistry::new();
        registry.plugins.push(Box::new(MockPrefetchPlugin {
            name: "Ahead".to_string(),
            prefetch_count: Arc::clone(&ahead_prefetches),
            annotate_count: Arc::clone(&ahead_annotates),
        }));
        registry.plugins.push(Box::new(MockRowReadingPlugin {
            inner: MockPrefetchPlugin {
                name: "Later".to_string(),
                prefetch_count: Arc::clone(&later_prefetches),
                annotate_count: Arc::clone(&later_annotates),
            },
        }));
        registry.plugins.push(Box::new(MockFallbackPlugin {
            run_batch_count: Arc::clone(&fallbacks),
        }));
        let variant = || InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec());

        let mut variants = vec![variant()];
        let prefetch = registry.prefetch_ahead(&variants);
        assert_eq!(
            prefetch
                .0
                .iter()
                .map(|slot| slot.is_some())
                .collect::<Vec<_>>(),
            [true, false, true]
        );
        assert_eq!(ahead_prefetches.load(Ordering::Relaxed), 1);
        assert_eq!(later_prefetches.load(Ordering::Relaxed), 0);
        assert!(variants[0].plugin_data.is_empty());

        assert!(registry
            .annotate_batch(&mut variants, prefetch)
            .failures
            .is_empty());
        assert_eq!(ahead_prefetches.load(Ordering::Relaxed), 1);
        assert_eq!(later_prefetches.load(Ordering::Relaxed), 1);
        assert_eq!(ahead_annotates.load(Ordering::Relaxed), 1);
        assert_eq!(later_annotates.load(Ordering::Relaxed), 1);
        assert_eq!(fallbacks.load(Ordering::Relaxed), 1);

        let mut whole = vec![variant()];
        assert!(registry.run_batch(&mut whole).failures.is_empty());
        assert_eq!(variants[0].plugin_data, whole[0].plugin_data);
        assert_eq!(
            variants[0].plugin_data.keys().collect::<Vec<_>>(),
            ["Ahead_ran", "Later_ran", "fallback_ran"]
        );
    }

    /// Which step of a test plugin fails.
    #[derive(Clone, Copy, PartialEq)]
    enum Fails {
        Nothing,
        Prefetch,
        Annotate,
    }

    /// A two-phase plugin named `name` that fails at the step `fails` names with
    /// the error text `boom` and otherwise writes `<NAME>=1` on every variant.
    struct StepPlugin {
        name: &'static str,
        fails: Fails,
        annotate_calls: Arc<AtomicUsize>,
    }

    impl BuiltinPlugin for StepPlugin {
        fn name(&self) -> &str {
            self.name
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn prefetch(
            &self,
            _variants: &[InputVariant],
        ) -> Result<Option<PrefetchData>, PluginError> {
            if self.fails == Fails::Prefetch {
                return Err(PluginError::Io(std::io::Error::other("boom")));
            }
            Ok(Some(PrefetchData::new(())))
        }
        fn annotate(
            &self,
            variants: &mut [InputVariant],
            _data: PrefetchData,
        ) -> Result<(), PluginError> {
            self.annotate_calls.fetch_add(1, Ordering::Relaxed);
            if self.fails == Fails::Annotate {
                return Err(PluginError::Run("boom".to_string()));
            }
            for v in variants.iter_mut() {
                v.plugin_data
                    .insert(self.name.to_uppercase(), "1".to_string());
            }
            Ok(())
        }
        fn run(
            &self,
            _consequence: &TranscriptConsequence,
            _variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            Ok(IndexMap::new())
        }
    }

    /// A registry of `plugins` run over two variants through each of the four
    /// batch entry points, the prefetch run ahead of annotation then
    /// `annotate_batch` among them; returns each way's variants and outcome.
    fn run_every_way(plugins: &[(&'static str, Fails)]) -> Vec<(Vec<InputVariant>, BatchOutcome)> {
        (0..4)
            .map(|way| {
                let mut registry = BuiltinRegistry::new();
                for (name, fails) in plugins {
                    registry.plugins.push(Box::new(StepPlugin {
                        name,
                        fails: *fails,
                        annotate_calls: Arc::new(AtomicUsize::new(0)),
                    }));
                }
                let mut variants = vec![
                    InputVariant::new("1".into(), 100, 100, b"A".to_vec(), b"G".to_vec()),
                    InputVariant::new("1".into(), 200, 200, b"C".to_vec(), b"T".to_vec()),
                ];
                let outcome = match way {
                    0 => registry.run_batch(&mut variants),
                    1 => registry.run_batch_timed(&mut variants).1,
                    2 => registry.run_batch_sequential(&mut variants).1,
                    _ => {
                        let ahead = registry.prefetch_ahead(&variants);
                        registry.annotate_batch(&mut variants, ahead)
                    }
                };
                (variants, outcome)
            })
            .collect()
    }

    /// The note a failure in `prefetch` or `annotate` carries after its error
    /// text, naming the first variant of the batch the plugin failed on.
    fn failed_on_note(first: &str) -> String {
        format!("the plugin annotated nothing on the variants it failed on, the first at {first}")
    }

    /// The note a failure in `run_batch` carries after its error text, naming
    /// the variant of the first row the plugin failed on.
    fn failed_rows_note(first: &str) -> String {
        format!(
            "the plugin annotated every row but the ones it failed on, the first of them on the \
             variant at {first}"
        )
    }

    /// A per-row plugin whose `run` fails with the error text `boom` on a row
    /// whose `exon` reads `fails` and otherwise writes `ROW=1` on the row;
    /// counts its `run` calls.
    struct RowPlugin {
        runs: Arc<AtomicUsize>,
    }

    impl BuiltinPlugin for RowPlugin {
        fn name(&self) -> &str {
            "Row"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn run(
            &self,
            consequence: &TranscriptConsequence,
            _variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            self.runs.fetch_add(1, Ordering::Relaxed);
            if consequence.exon.as_deref() == Some("fails") {
                return Err(PluginError::Run("boom".to_string()));
            }
            Ok(IndexMap::from([("ROW".to_string(), "1".to_string())]))
        }
    }

    /// A transcript row [`RowPlugin`] `fails` on, or writes on.
    fn transcript_row(fails: bool) -> TranscriptConsequence {
        TranscriptConsequence {
            exon: fails.then(|| "fails".to_string()),
            ..row(FeatureType::Transcript)
        }
    }

    /// Ensembl VEP's `run_plugins` evals each row's plugin call on its own
    /// (OutputFactory.pm:2225-2244), so a row the plugin fails on loses the
    /// plugin's fields alone: the rows before and after it and the other
    /// variant's rows keep the field, in a batch of two (retried variant by
    /// variant) and in a batch of one (its own retry), and the outcome's note
    /// states that footprint, naming the variant the first failing row is on,
    /// with the one variant failed on counted. The batch pass wrote every other
    /// row, so the variants after the first failing one are run on a copy for
    /// the count alone: three rows per variant, the batch pass runs six, the
    /// retry three per variant.
    #[test]
    fn a_failing_row_loses_the_field_alone_and_the_note_states_the_rows_footprint() {
        for (starts, failing, runs_expected) in [
            (&[100, 200][..], 100, 6 + 3 + 3),
            (&[100, 200][..], 200, 6 + 3 + 3),
            (&[100][..], 100, 3),
        ] {
            let runs = Arc::new(AtomicUsize::new(0));
            let mut registry = BuiltinRegistry::new();
            registry.plugins.push(Box::new(RowPlugin {
                runs: Arc::clone(&runs),
            }));
            let mut variants: Vec<InputVariant> = starts.iter().map(|&s| snv_at(s)).collect();
            for variant in &mut variants {
                variant.transcript_consequences = vec![
                    transcript_row(false),
                    transcript_row(variant.start == failing),
                    transcript_row(false),
                ];
            }
            let outcome = registry.run_batch(&mut variants);
            let [PluginFailure {
                plugin,
                error: PluginError::Run(reason),
                failed_variants,
                ..
            }] = &outcome.failures[..]
            else {
                panic!("one Run failure, got {:?}", outcome.failures);
            };
            assert_eq!(plugin, "Row");
            assert_eq!(*failed_variants, 1);
            assert_eq!(
                *reason,
                format!("boom; {}", failed_rows_note(&format!("1:{failing}")))
            );
            let rows: Vec<Vec<Option<&str>>> = variants
                .iter()
                .map(|variant| {
                    variant
                        .transcript_consequences
                        .iter()
                        .map(|tc| tc.plugin_data.get("ROW").map(String::as_str))
                        .collect()
                })
                .collect();
            let expected: Vec<Vec<Option<&str>>> = variants
                .iter()
                .map(|variant| {
                    if variant.start == failing {
                        vec![Some("1"), None, Some("1")]
                    } else {
                        vec![Some("1"); 3]
                    }
                })
                .collect();
            assert_eq!(rows, expected, "{starts:?} failing {failing}");
            assert_eq!(
                runs.load(Ordering::Relaxed),
                runs_expected,
                "{starts:?} failing {failing}"
            );
        }
    }

    /// `run_plugins` catches a plugin's failure per line and runs the next plugin
    /// (OutputFactory.pm:2243-2245, 2250): a plugin whose annotate fails on every
    /// variant costs its own fields alone, the plugin after it in load order still
    /// annotates every variant, and the outcome names the failure once, with the
    /// plugin's name, the error, the first variant it failed on and the count of
    /// the batch's variants it failed on, on each of the four entry points.
    #[test]
    fn a_failing_plugin_does_not_stop_the_plugins_after_it() {
        let plugins = [("Boom", Fails::Annotate), ("Ok", Fails::Nothing)];
        for (variants, outcome) in run_every_way(&plugins) {
            for variant in &variants {
                assert_eq!(variant.plugin_data.get("OK").unwrap(), "1");
                assert!(!variant.plugin_data.contains_key("BOOM"));
            }
            let [PluginFailure {
                plugin,
                error: PluginError::Run(reason),
                failed_variants,
                ..
            }] = &outcome.failures[..]
            else {
                panic!("one Run failure, got {:?}", outcome.failures);
            };
            assert_eq!(plugin, "Boom");
            assert_eq!(*reason, format!("boom; {}", failed_on_note("1:100")));
            assert_eq!(*failed_variants, variants.len());
        }
    }

    /// A plugin whose prefetch fails is skipped in the annotate phase, named in
    /// the outcome with the prefetch error, and the plugins before and after it
    /// annotate as if it were absent; a prefetch that failed ahead of annotation
    /// travels in its slot and is named the same way.
    #[test]
    fn a_failing_prefetch_skips_only_that_plugin() {
        let plugins = [
            ("First", Fails::Nothing),
            ("Broken", Fails::Prefetch),
            ("Last", Fails::Nothing),
        ];
        for (variants, outcome) in run_every_way(&plugins) {
            for variant in &variants {
                assert_eq!(variant.plugin_data.get("FIRST").unwrap(), "1");
                assert_eq!(variant.plugin_data.get("LAST").unwrap(), "1");
                assert!(!variant.plugin_data.contains_key("BROKEN"));
            }
            let [PluginFailure {
                plugin,
                error: PluginError::Io(error),
                failed_variants,
                ..
            }] = &outcome.failures[..]
            else {
                panic!("one Io failure, got {:?}", outcome.failures);
            };
            assert_eq!(plugin, "Broken");
            assert_eq!(*failed_variants, variants.len());
            assert_eq!(
                error.to_string(),
                format!("boom; {}", failed_on_note("1:100"))
            );
        }
        let mut registry = BuiltinRegistry::new();
        let annotate_calls = Arc::new(AtomicUsize::new(0));
        registry.plugins.push(Box::new(StepPlugin {
            name: "Broken",
            fails: Fails::Prefetch,
            annotate_calls: Arc::clone(&annotate_calls),
        }));
        let mut variants = vec![InputVariant::new(
            "1".into(),
            100,
            100,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        assert_eq!(registry.run_batch(&mut variants).failures.len(), 1);
        assert_eq!(annotate_calls.load(Ordering::Relaxed), 0);
    }

    /// A two-phase plugin whose prefetch fails with the error text `bad block`
    /// on every batch `fails` accepts and otherwise writes `SPAN=1` on every
    /// variant of the batch; counts its prefetch calls.
    struct SpanPlugin {
        fails: fn(&[InputVariant]) -> bool,
        prefetches: Arc<AtomicUsize>,
    }

    impl BuiltinPlugin for SpanPlugin {
        fn name(&self) -> &str {
            "Span"
        }
        fn header_info(&self) -> Vec<(String, String)> {
            vec![]
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        fn prefetch(&self, variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
            self.prefetches.fetch_add(1, Ordering::Relaxed);
            if (self.fails)(variants) {
                return Err(PluginError::Run("bad block".to_string()));
            }
            Ok(Some(PrefetchData::new(())))
        }
        fn annotate(
            &self,
            variants: &mut [InputVariant],
            _data: PrefetchData,
        ) -> Result<(), PluginError> {
            for v in variants.iter_mut() {
                v.plugin_data.insert("SPAN".to_string(), "1".to_string());
            }
            Ok(())
        }
    }

    /// A single-base variant of chromosome 1 at `start`.
    fn snv_at(start: u64) -> InputVariant {
        InputVariant::new("1".into(), start, start, b"A".to_vec(), b"G".to_vec())
    }

    /// A one-plugin registry of a [`SpanPlugin`] failing on `fails` run over
    /// `starts`' variants through each of the four batch entry points; returns
    /// each way's variants, outcome and prefetch count.
    fn span_plugin_every_way(
        fails: fn(&[InputVariant]) -> bool,
        starts: &[u64],
    ) -> Vec<(Vec<InputVariant>, BatchOutcome, usize)> {
        (0..4)
            .map(|way| {
                let prefetches = Arc::new(AtomicUsize::new(0));
                let mut registry = BuiltinRegistry::new();
                registry.plugins.push(Box::new(SpanPlugin {
                    fails,
                    prefetches: Arc::clone(&prefetches),
                }));
                let mut variants: Vec<InputVariant> = starts.iter().map(|&s| snv_at(s)).collect();
                let outcome = match way {
                    0 => registry.run_batch(&mut variants),
                    1 => registry.run_batch_timed(&mut variants).1,
                    2 => registry.run_batch_sequential(&mut variants).1,
                    _ => {
                        let ahead = registry.prefetch_ahead(&variants);
                        registry.annotate_batch(&mut variants, ahead)
                    }
                };
                (variants, outcome, prefetches.load(Ordering::Relaxed))
            })
            .collect()
    }

    /// The `Run` reason of an outcome's one failure, which names `Span`, and
    /// the count of variants failed on.
    fn span_failure(outcome: &BatchOutcome) -> (String, usize) {
        let [PluginFailure {
            plugin,
            error: PluginError::Run(reason),
            failed_variants,
            ..
        }] = &outcome.failures[..]
        else {
            panic!("one Run failure, got {:?}", outcome.failures);
        };
        assert_eq!(plugin, "Span");
        (reason.clone(), *failed_variants)
    }

    /// Ensembl VEP 116.2's plugins read each line's own span inside `run_plugins`'s
    /// per-line `eval` (REVEL.pm:168 `get_data` per variation feature,
    /// OutputFactory.pm:2225-2226), so a span that cannot be read costs its own
    /// line the field and no other: a batch whose prefetch fails is retried
    /// variant by variant, the variants whose own span reads keep the field, the
    /// failing one has none, and the outcome names the plugin once with the error
    /// and the variant it annotated nothing on, on each of the four entry points.
    #[test]
    fn a_failing_batch_is_retried_variant_by_variant_and_the_failing_span_alone_loses_the_field() {
        let fails = |variants: &[InputVariant]| variants.iter().any(|v| v.start == 200);
        for (variants, outcome, prefetches) in span_plugin_every_way(fails, &[100, 200, 300]) {
            let spans: Vec<Option<&str>> = variants
                .iter()
                .map(|v| v.plugin_data.get("SPAN").map(String::as_str))
                .collect();
            assert_eq!(spans, [Some("1"), None, Some("1")]);
            assert_eq!(
                span_failure(&outcome),
                (format!("bad block; {}", failed_on_note("1:200")), 1)
            );
            assert_eq!(prefetches, 1 + 3, "the batch, then each variant");
        }
    }

    /// Each batch's record names the first variant of that batch the plugin
    /// failed on and counts the batch's failed variants; the warning sink keeps
    /// the first batch's text and sums the counts into the one line it writes
    /// for the run (Ensembl VEP's `warning_msg` writes each text once,
    /// BaseVEP.pm:810-813); the variants whose spans read keep the field in
    /// every batch.
    #[test]
    fn every_batch_of_a_run_names_the_first_variant_the_plugin_failed_on() {
        let mut registry = BuiltinRegistry::new();
        registry.plugins.push(Box::new(SpanPlugin {
            fails: |variants| variants.iter().any(|v| v.start == 200 || v.start == 400),
            prefetches: Arc::new(AtomicUsize::new(0)),
        }));
        let mut first_batch = vec![snv_at(100), snv_at(200), snv_at(300)];
        let mut second_batch = vec![snv_at(400), snv_at(500)];
        let first = registry.run_batch(&mut first_batch);
        let second = registry.run_batch(&mut second_batch);
        let with_field = |variants: &[InputVariant]| -> Vec<u64> {
            variants
                .iter()
                .filter(|v| v.plugin_data.contains_key("SPAN"))
                .map(|v| v.start)
                .collect()
        };
        assert_eq!(with_field(&first_batch), [100, 300]);
        assert_eq!(with_field(&second_batch), [500]);
        assert_eq!(
            span_failure(&first),
            (format!("bad block; {}", failed_on_note("1:200")), 1)
        );
        assert_eq!(
            span_failure(&second),
            (format!("bad block; {}", failed_on_note("1:400")), 1)
        );
    }

    /// A built-in loaded twice is two plugins of one name, each with its own
    /// data file, as Ensembl VEP loads one instance per `--plugin` argument
    /// (Runner.pm:925-954): each instance's record names the first variant that
    /// instance failed on, in load order.
    #[test]
    fn two_instances_of_one_plugin_name_each_name_their_own_first_variant() {
        let mut registry = BuiltinRegistry::new();
        for fails in [
            (|variants: &[InputVariant]| variants.iter().any(|v| v.start == 300))
                as fn(&[InputVariant]) -> bool,
            |variants: &[InputVariant]| variants.iter().any(|v| v.start == 100),
        ] {
            registry.plugins.push(Box::new(SpanPlugin {
                fails,
                prefetches: Arc::new(AtomicUsize::new(0)),
            }));
        }
        let mut variants = vec![snv_at(100), snv_at(200), snv_at(300)];
        let outcome = registry.run_batch(&mut variants);
        let instances: Vec<usize> = outcome.failures.iter().map(|f| f.instance).collect();
        assert_eq!(instances, [0, 1]);
        let reasons: Vec<String> = outcome
            .failures
            .iter()
            .map(|failure| {
                assert_eq!(failure.plugin, "Span");
                assert_eq!(failure.failed_variants, 1);
                failure.error.to_string()
            })
            .collect();
        assert_eq!(
            reasons,
            [
                format!(
                    "plugin execution failed: bad block; {}",
                    failed_on_note("1:300")
                ),
                format!(
                    "plugin execution failed: bad block; {}",
                    failed_on_note("1:100")
                ),
            ]
        );
    }

    /// A batch failure that every variant's own retry recovers from (a read
    /// across several variants' spans that fails where each span alone reads)
    /// leaves every variant its field and records nothing, as Ensembl VEP, which
    /// reads one span at a time, meets no failure there.
    #[test]
    fn a_batch_failure_the_retry_recovers_everywhere_is_not_recorded() {
        let fails = |variants: &[InputVariant]| variants.len() > 1;
        for (variants, outcome, prefetches) in span_plugin_every_way(fails, &[100, 200, 300]) {
            for variant in &variants {
                assert_eq!(variant.plugin_data.get("SPAN").unwrap(), "1");
            }
            assert!(outcome.failures.is_empty(), "{:?}", outcome.failures);
            assert_eq!(prefetches, 1 + 3);
        }
    }

    /// A one-variant batch that fails is its own span: the failure is recorded
    /// with the batch's error and the variant's location, and the plugin is not
    /// called again.
    #[test]
    fn a_one_variant_batch_fails_once_without_a_retry() {
        let fails = |_: &[InputVariant]| true;
        for (variants, outcome, prefetches) in span_plugin_every_way(fails, &[100]) {
            assert!(!variants[0].plugin_data.contains_key("SPAN"));
            assert_eq!(
                span_failure(&outcome),
                (format!("bad block; {}", failed_on_note("1:100")), 1)
            );
            assert_eq!(prefetches, 1);
        }
    }

    /// A plugin whose header fields are those it declares, in the order it
    /// declares them.
    struct MockHeaderPlugin {
        name: &'static str,
        fields: Vec<(&'static str, &'static str)>,
    }

    impl BuiltinPlugin for MockHeaderPlugin {
        fn name(&self) -> &str {
            self.name
        }
        fn header_info(&self) -> Vec<(String, String)> {
            self.fields
                .iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        }
        fn init(&mut self, _params: &[String]) -> Result<(), PluginError> {
            Ok(())
        }
        /// Marks the row with the first declared field.
        fn run(
            &self,
            _consequence: &TranscriptConsequence,
            _variant: &InputVariant,
        ) -> Result<IndexMap<String, String>, PluginError> {
            let mut result = IndexMap::new();
            if let Some((field, _)) = self.fields.first() {
                result.insert(field.to_string(), "1".to_string());
            }
            Ok(result)
        }
    }

    fn row(feature_type: FeatureType) -> TranscriptConsequence {
        TranscriptConsequence {
            feature_type,
            ..TranscriptConsequence::default()
        }
    }

    /// Every built-in whose Ensembl VEP module declares `feature_types` of
    /// `['Transcript']`, or does not override the default (`BaseVepPlugin.pm:99`),
    /// runs on transcript rows alone: the default `run_batch` writes on a
    /// record's transcript row and on neither of its regulatory rows, and the
    /// registry reports the plugin's fields as transcript-row fields.
    #[test]
    fn plugins_run_on_transcript_rows_only_by_default() {
        for name in [
            "REVEL",
            "SpliceAI",
            "AlphaMissense",
            "dbNSFP",
            "dbscSNV",
            "LoFtool",
            "LoFTEE",
            "pLI",
        ] {
            let plugin = create_builtin(name).unwrap();
            assert_eq!(plugin.feature_types(), [RowKind::Transcript], "{name}");
        }
        let mut registry = BuiltinRegistry::new();
        registry.plugins.push(Box::new(MockHeaderPlugin {
            name: "Marker",
            fields: vec![("Marker", "marks the rows it runs on")],
        }));
        assert_eq!(
            registry.field_row_kinds(),
            [("Marker".to_string(), vec![RowKind::Transcript])]
        );
        let mut variants = vec![InputVariant::new(
            "21".into(),
            100,
            100,
            b"A".to_vec(),
            b"G".to_vec(),
        )];
        variants[0].transcript_consequences.extend([
            row(FeatureType::RegulatoryFeature),
            row(FeatureType::Transcript),
            row(FeatureType::MotifFeature),
        ]);
        assert!(registry.run_batch(&mut variants).failures.is_empty());
        let written: Vec<bool> = variants[0]
            .transcript_consequences
            .iter()
            .map(|tc| tc.plugin_data.contains_key("Marker"))
            .collect();
        assert_eq!(written, [false, true, false]);
    }

    /// Every built-in whose Ensembl VEP module declares `feature_types` of
    /// `['Feature', 'Intergenic']` runs on every kind of row, so its fields
    /// print on regulatory, motif and intergenic rows too. None of these modules
    /// names a motif or regulatory kind by itself, so none switches regulatory
    /// annotation on (`Runner.pm:999` matches the declared words, `Feature`
    /// among them not).
    #[test]
    fn plugins_declaring_feature_and_intergenic_run_on_every_kind_of_row() {
        for name in [
            "CADD",
            "gnomADc",
            "GWAS",
            "BayesDel",
            "PrimateAI",
            "gnomADMt",
            "AncestralAllele",
            "AVI",
            "OpenTargets",
            "Conservation",
        ] {
            let plugin = create_builtin(name).unwrap();
            assert_eq!(plugin.feature_types(), RowKind::ALL, "{name}");
            for feature_type in [
                FeatureType::Transcript,
                FeatureType::RegulatoryFeature,
                FeatureType::MotifFeature,
            ] {
                assert!(plugin.annotates_row(&row(feature_type)), "{name}");
            }
        }
        let mut registry = BuiltinRegistry::new();
        registry.plugins.push(create_builtin("CADD").unwrap());
        for (_, kinds) in registry.field_row_kinds() {
            assert_eq!(kinds, RowKind::ALL);
        }
    }

    /// `OutputFactory.pm` `get_plugin_headers` walks the plugins in load order and
    /// each plugin's header hash by `sort keys`, Perl's byte-order `cmp`, so an
    /// upper-case name precedes a lower-case one and the second plugin's fields
    /// follow the first's whatever the plugin or field names.
    #[test]
    fn header_info_is_sorted_within_each_plugin_in_load_order() {
        let mut registry = BuiltinRegistry::new();
        registry.plugins.push(Box::new(MockHeaderPlugin {
            name: "Zulu",
            fields: vec![
                ("zeta", "last letter"),
                ("alpha", "first letter"),
                ("Mu", "upper"),
            ],
        }));
        registry.plugins.push(Box::new(MockHeaderPlugin {
            name: "Alfa",
            fields: vec![("beta", "second letter"), ("Alpha", "upper first")],
        }));

        let headers = registry.all_header_info();

        let pairs: Vec<(&str, &str)> = headers
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
        assert_eq!(
            pairs,
            [
                ("Mu", "upper"),
                ("alpha", "first letter"),
                ("zeta", "last letter"),
                ("Alpha", "upper first"),
                ("beta", "second letter"),
            ]
        );
    }
}
