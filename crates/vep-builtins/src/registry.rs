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
//!    parallel via rayon, reading from independent data files.
//! 2. **Annotate** (sequential): prefetched results are applied to variants.
//!
//! Plugins that don't implement prefetch fall back to sequential `run_batch`.

use std::time::Instant;

use rayon::prelude::*;
use tracing::{debug, info};

use crate::plugins;
use crate::traits::{BuiltinPlugin, PluginError, PrefetchData};
use vep_core::variant::InputVariant;

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

    /// Load plugins from `--plugin` argument strings.
    ///
    /// Each argument has the format `"PluginName,param1,param2"`.
    /// Returns names of plugins that were not found as builtins (for dylib fallback).
    pub fn load_from_args(&mut self, plugin_args: &[String]) -> Result<Vec<String>, PluginError> {
        self.load_from_args_with_context(plugin_args, None, None)
    }

    /// Load plugins with the full pipeline context: the shared reference FASTA
    /// (from `--fasta`) and the target assembly (from `--assembly`).
    ///
    /// The assembly is injected **before** `init()` so a plugin can use it while
    /// opening its data files. That ordering is load-bearing: REVEL must know the
    /// assembly to pick which coordinate column tabix indexed, and LoFTEE must
    /// know it to choose between a bigWig and a tabix-TSV GERP reader. The FASTA
    /// is injected after `init()`.
    pub fn load_from_args_with_context(
        &mut self,
        plugin_args: &[String],
        reference_fasta: Option<std::sync::Arc<vep_fasta::IndexedFasta>>,
        assembly: Option<vep_core::assembly::Assembly>,
    ) -> Result<Vec<String>, PluginError> {
        let mut unresolved = Vec::new();

        for arg in plugin_args {
            let parts: Vec<&str> = arg.splitn(2, ',').collect();
            let name = parts[0];
            let params: Vec<String> = if parts.len() > 1 {
                parts[1].split(',').map(String::from).collect()
            } else {
                Vec::new()
            };

            if let Some(mut plugin) = create_builtin(name) {
                plugin.set_assembly(assembly);
                plugin.init(&params)?;
                plugin.set_reference_fasta(reference_fasta.clone());
                info!(name = plugin.name(), "loaded builtin plugin");
                self.plugins.push(plugin);
            } else {
                unresolved.push(arg.clone());
            }
        }

        Ok(unresolved)
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
    /// [`prefetch`]: BuiltinPlugin::prefetch
    /// [`annotate`]: BuiltinPlugin::annotate
    /// [`run_batch`]: BuiltinPlugin::run_batch
    pub fn run_batch(&self, variants: &mut [InputVariant]) -> Result<(), PluginError> {
        if self.plugins.is_empty() {
            return Ok(());
        }

        let batch_size = variants.len();
        let overall_start = Instant::now();

        // The mutable slice is reborrowed as shared for prefetch; no plugin
        // writes to variants in that phase.
        let variants_shared: &[InputVariant] = variants;
        let prefetch_start = Instant::now();

        let prefetched: Vec<Result<Option<PrefetchData>, PluginError>> = self
            .plugins
            .par_iter()
            .map(|plugin| plugin.prefetch(variants_shared))
            .collect();

        let prefetch_elapsed = prefetch_start.elapsed();
        debug!(
            batch_size,
            prefetch_ms = prefetch_elapsed.as_millis() as u64,
            plugins = self.plugins.len(),
            "parallel prefetch phase complete"
        );

        for (plugin, prefetch_result) in self.plugins.iter().zip(prefetched) {
            let started = Instant::now();
            match prefetch_result? {
                Some(data) => {
                    plugin.annotate(variants, data)?;
                }
                None => {
                    plugin.run_batch(variants)?;
                }
            }
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

        Ok(())
    }

    /// Like [`run_batch`] but returns structured per-plugin timing data.
    ///
    /// Used by benchmarking tools to measure plugin overhead. The pipeline
    /// behavior is identical to `run_batch`.
    pub fn run_batch_timed(
        &self,
        variants: &mut [InputVariant],
    ) -> Result<PluginTimings, PluginError> {
        let mut timings = PluginTimings {
            prefetch_ms: Vec::with_capacity(self.plugins.len()),
            annotate_ms: Vec::with_capacity(self.plugins.len()),
            total_prefetch_ms: 0,
            total_annotate_ms: 0,
            total_ms: 0,
        };

        if self.plugins.is_empty() {
            return Ok(timings);
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

        for (plugin, (prefetch_result, _ms)) in self.plugins.iter().zip(prefetched) {
            let t = Instant::now();
            match prefetch_result? {
                Some(data) => {
                    plugin.annotate(variants, data)?;
                }
                None => {
                    plugin.run_batch(variants)?;
                }
            }
            let elapsed_ms = t.elapsed().as_millis() as u64;
            timings
                .annotate_ms
                .push((plugin.name().to_string(), elapsed_ms));
        }

        timings.total_annotate_ms = annotate_wall_start.elapsed().as_millis() as u64;
        timings.total_ms = overall_start.elapsed().as_millis() as u64;

        Ok(timings)
    }

    /// Run all plugins **sequentially** (no parallel prefetch).
    ///
    /// Used by benchmarks to measure the speedup from parallel prefetch.
    /// Each plugin's prefetch + annotate runs one at a time.
    pub fn run_batch_sequential(
        &self,
        variants: &mut [InputVariant],
    ) -> Result<PluginTimings, PluginError> {
        let mut timings = PluginTimings {
            prefetch_ms: Vec::with_capacity(self.plugins.len()),
            annotate_ms: Vec::with_capacity(self.plugins.len()),
            total_prefetch_ms: 0,
            total_annotate_ms: 0,
            total_ms: 0,
        };

        if self.plugins.is_empty() {
            return Ok(timings);
        }

        let overall_start = Instant::now();
        let mut total_prefetch = 0u64;
        let mut total_annotate = 0u64;

        for plugin in &self.plugins {
            let t = Instant::now();
            let prefetch_result = plugin.prefetch(variants)?;
            let prefetch_ms = t.elapsed().as_millis() as u64;
            timings
                .prefetch_ms
                .push((plugin.name().to_string(), prefetch_ms));
            total_prefetch += prefetch_ms;

            let t = Instant::now();
            match prefetch_result {
                Some(data) => {
                    plugin.annotate(variants, data)?;
                }
                None => {
                    plugin.run_batch(variants)?;
                }
            }
            let annotate_ms = t.elapsed().as_millis() as u64;
            timings
                .annotate_ms
                .push((plugin.name().to_string(), annotate_ms));
            total_annotate += annotate_ms;
        }

        timings.total_prefetch_ms = total_prefetch;
        timings.total_annotate_ms = total_annotate;
        timings.total_ms = overall_start.elapsed().as_millis() as u64;

        Ok(timings)
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

    /// Returns header info from all loaded plugins.
    pub fn all_header_info(&self) -> Vec<(String, String)> {
        let mut headers = Vec::new();
        for plugin in &self.plugins {
            headers.extend(plugin.header_info());
        }
        headers
    }
}

impl Default for BuiltinRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// Attempt to create a built-in plugin by name.
///
/// Returns `None` if the name doesn't match any built-in plugin.
fn create_builtin(name: &str) -> Option<Box<dyn BuiltinPlugin>> {
    match name {
        "CADD" => Some(Box::new(plugins::cadd::CaddPlugin::new())),
        "gnomADc" => Some(Box::new(plugins::gnomadc::GnomadcPlugin::new())),
        "REVEL" => Some(Box::new(plugins::revel::RevelPlugin::new())),
        "SpliceAI" => Some(Box::new(plugins::spliceai::SpliceAiPlugin::new())),
        "AlphaMissense" => Some(Box::new(plugins::alpha_missense::AlphaMissensePlugin::new())),
        "dbNSFP" => Some(Box::new(plugins::dbnsfp::DbNsfpPlugin::new())),
        "dbscSNV" => Some(Box::new(plugins::dbscsnv::DbscSnvPlugin::new())),
        "GWAS" => Some(Box::new(plugins::gwas::GwasPlugin::new())),
        "LoFtool" => Some(Box::new(plugins::loftool::LoFtoolPlugin::new())),
        "LoFTEE" => Some(Box::new(plugins::loftee::LofteePlugin::new())),
        "pLI" => Some(Box::new(plugins::pli::PliPlugin::new())),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use indexmap::IndexMap;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;
    use vep_core::consequence::TranscriptConsequence;

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

        registry.run_batch(&mut variants).unwrap();

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

        registry.run_batch(&mut variants).unwrap();

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

        registry.run_batch(&mut variants).unwrap();

        assert_eq!(prefetch_count_a.load(Ordering::Relaxed), 1);
        assert_eq!(annotate_count_a.load(Ordering::Relaxed), 1);
        assert_eq!(prefetch_count_b.load(Ordering::Relaxed), 1);
        assert_eq!(annotate_count_b.load(Ordering::Relaxed), 1);
        assert_eq!(run_batch_count.load(Ordering::Relaxed), 1);

        assert_eq!(variants[0].plugin_data.get("PluginA_ran").unwrap(), "true");
        assert_eq!(variants[0].plugin_data.get("PluginB_ran").unwrap(), "true");
        assert_eq!(variants[0].plugin_data.get("fallback_ran").unwrap(), "true");
    }
}
