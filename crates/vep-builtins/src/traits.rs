// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Plugin trait definition for built-in VEP plugins.
//!
//! Built-in plugins operate on native Rust types without JSON serialization
//! or C FFI overhead. The plugin execution model has three layers:
//!
//! - [`prefetch`](BuiltinPlugin::prefetch) + [`annotate`](BuiltinPlugin::annotate):
//!   Two-phase execution that separates I/O (parallelizable) from annotation
//!   (sequential). Tabix plugins implement these to enable parallel prefetch
//!   across independent data files.
//! - [`run_batch`](BuiltinPlugin::run_batch): Single-phase fallback for plugins
//!   that don't benefit from parallel I/O (e.g., in-memory lookups).
//! - [`run`](BuiltinPlugin::run): Per-transcript-consequence annotation for
//!   simple per-transcript plugins.

use std::sync::Arc;

use indexmap::IndexMap;
use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;

/// Errors from builtin plugin execution.
#[derive(Debug, thiserror::Error)]
pub enum PluginError {
    #[error("plugin initialization failed: {0}")]
    Init(String),
    #[error("plugin execution failed: {0}")]
    Run(String),
    #[error("tabix I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Type-erased container for data returned by [`BuiltinPlugin::prefetch`].
///
/// Each plugin stores its own concrete type (e.g., `BatchQueryResult`) inside
/// the `Box<dyn Any>` and downcasts it in [`BuiltinPlugin::annotate`]. This
/// keeps the trait free of backend-specific types.
pub struct PrefetchData(pub Box<dyn std::any::Any + Send>);

impl PrefetchData {
    /// Create a new `PrefetchData` wrapping a concrete value.
    pub fn new<T: std::any::Any + Send + 'static>(data: T) -> Self {
        Self(Box::new(data))
    }

    /// Downcast the inner data to a concrete type.
    ///
    /// Returns `Err(PluginError::Run)` if the downcast fails, which indicates
    /// a bug in the plugin (prefetch/annotate type mismatch).
    pub fn downcast<T: 'static>(self) -> Result<T, PluginError> {
        self.0
            .downcast::<T>()
            .map(|b| *b)
            .map_err(|_| PluginError::Run("prefetch data type mismatch".into()))
    }
}

/// Trait that all built-in VEP plugins implement.
///
/// Plugins are compiled directly into the binary and operate on native Rust
/// types. There are three execution strategies (the registry picks the best
/// one automatically):
///
/// 1. **Two-phase (prefetch + annotate):** The registry calls [`prefetch`]
///    on all plugins in parallel (via rayon), then calls [`annotate`]
///    sequentially to write results into variants. Tabix plugins use this
///    to overlap I/O from independent data files.
///
/// 2. **Single-phase (run_batch):** Fallback for plugins that return `None`
///    from [`prefetch`]. Called sequentially with `&mut` access to variants.
///
/// 3. **Per-transcript (run):** Default implementation of [`run_batch`]
///    dispatches to [`run`] for each transcript consequence.
///
/// [`prefetch`]: BuiltinPlugin::prefetch
/// [`annotate`]: BuiltinPlugin::annotate
pub trait BuiltinPlugin: Send + Sync {
    /// Plugin name matching the `--plugin` CLI argument (e.g., "CADD").
    fn name(&self) -> &str;

    /// Header info for output columns: `(field_name, description)` pairs.
    fn header_info(&self) -> Vec<(String, String)>;

    /// Initialize the plugin with CLI parameters.
    ///
    /// Parameters come from `--plugin Name,param1,param2,...` after the name.
    fn init(&mut self, params: &[String]) -> Result<(), PluginError>;

    /// Inject the pipeline's shared reference FASTA, if one was provided via
    /// `--fasta`. Called by the registry after [`init`](BuiltinPlugin::init).
    ///
    /// Default is a no-op; only sequence-aware consequence-filtering plugins
    /// (LoFTEE, which reads intron boundary motifs and ancestral context) use
    /// it. Sharing the already-loaded `Arc` avoids opening the FASTA twice.
    fn set_reference_fasta(&mut self, _fasta: Option<Arc<vep_fasta::IndexedFasta>>) {}

    /// Inject the run's target assembly, resolved from `--assembly`. Called by
    /// the registry **before** [`init`](BuiltinPlugin::init), so `init` can use
    /// it while opening data files.
    ///
    /// Default is a no-op. Plugins that need it fall into two groups:
    ///
    /// - **Dual-coordinate files**: REVEL ships one file holding both
    ///   `hg19_pos` and `grch38_pos`, indexed on a different column per
    ///   assembly. Without the assembly the plugin cannot know which column
    ///   holds the coordinates tabix indexed, and reads the wrong one.
    /// - **Per-assembly bundles**: LoFTEE's upstream GERP data is a tabix TSV
    ///   on GRCh37 and a bigWig on GRCh38, with different END_TRUNC cutoffs.
    ///
    /// Plugins keyed on gene symbol (LoFtool, pLI) genuinely ignore this: one file
    /// serves both assemblies.
    ///
    /// Several plugins whose databases are pre-split per assembly override it, to
    /// *validate* that the supplied file matches the requested build rather than
    /// annotating from the wrong genome. A pre-split file is not self-evidently
    /// correct: the two copies share a column layout, so pointing a run at the wrong
    /// one resolves contigs, matches nothing, and exits 0. Each uses whichever
    /// in-file signal that database carries: CADD a `##`-preamble marker
    /// ([`tabix::detect_assembly_marker`](crate::tabix::detect_assembly_marker)),
    /// dbscSNV its coordinate column names, AlphaMissense its per-row `genome`
    /// column ([`tabix::first_data_row_field`](crate::tabix::first_data_row_field)).
    ///
    /// gnomADc deliberately does not: its two files are different
    /// releases (v2.1 genomes on GRCh37, v4.0 exomes on GRCh38) whose headers the
    /// plugin already reads to discover its output fields, and neither carries a
    /// build marker. There is no in-file signal to validate against, so a check here
    /// would have to infer the build from a contig-name convention, exactly the
    /// heuristic that would break on the next release.
    fn set_assembly(&mut self, _assembly: Option<vep_core::assembly::Assembly>) {}

    /// I/O phase: fetch annotation data for a batch of variants.
    ///
    /// This method takes a **shared** reference to the variant buffer, so the
    /// registry can call it on multiple plugins in parallel. The returned
    /// [`PrefetchData`] is an opaque container that will be passed to
    /// [`annotate`](BuiltinPlugin::annotate) in the sequential phase.
    ///
    /// Return `Ok(None)` to opt out of two-phase execution; the registry
    /// will fall back to calling [`run_batch`](BuiltinPlugin::run_batch).
    fn prefetch(&self, _variants: &[InputVariant]) -> Result<Option<PrefetchData>, PluginError> {
        Ok(None)
    }

    /// Annotation phase: apply prefetched data to variants.
    ///
    /// Called sequentially after all [`prefetch`](BuiltinPlugin::prefetch)
    /// calls complete. The `data` argument is the value returned from
    /// `prefetch` for this plugin.
    fn annotate(
        &self,
        _variants: &mut [InputVariant],
        _data: PrefetchData,
    ) -> Result<(), PluginError> {
        Ok(())
    }

    /// Batch-annotate a buffer of variants (single-phase fallback).
    ///
    /// Called once per buffer (typically 5000 variants) when the plugin does
    /// not implement the two-phase prefetch/annotate protocol. Tabix plugins
    /// should prefer implementing [`prefetch`](BuiltinPlugin::prefetch) +
    /// [`annotate`](BuiltinPlugin::annotate) instead.
    ///
    /// Default implementation calls [`run`](BuiltinPlugin::run) for each
    /// variant/consequence pair.
    fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        for variant in variants.iter_mut() {
            for tc_idx in 0..variant.transcript_consequences.len() {
                let result = self.run(&variant.transcript_consequences[tc_idx], variant)?;
                for (k, v) in result {
                    variant.transcript_consequences[tc_idx]
                        .plugin_data
                        .insert(k, v);
                }
            }
        }
        Ok(())
    }

    /// Annotate a single transcript consequence (per-transcript plugins).
    ///
    /// Returns key-value pairs to add to the Extra column. Return an empty
    /// map if the plugin has no annotation for this consequence.
    fn run(
        &self,
        _consequence: &TranscriptConsequence,
        _variant: &InputVariant,
    ) -> Result<IndexMap<String, String>, PluginError> {
        Ok(IndexMap::new())
    }
}
