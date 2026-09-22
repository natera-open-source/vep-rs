// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Plugin registry: discovers, loads, and manages plugin dylibs.

use std::path::{Path, PathBuf};

use crate::loader::LoadedPlugin;
use crate::PluginError;

/// Registry of loaded plugins.
#[derive(Debug)]
pub struct PluginRegistry {
    plugins: Vec<LoadedPlugin>,
}

impl PluginRegistry {
    /// Create an empty plugin registry.
    pub fn new() -> Self {
        Self {
            plugins: Vec::new(),
        }
    }

    /// Load plugins from `--plugin` argument strings.
    ///
    /// Each argument has the format `"PluginName,param1,param2"` or
    /// `"/path/to/plugin.so,param1,param2"`. The first comma-separated token
    /// is the plugin name or path; the rest are init parameters.
    ///
    /// `plugin_dir` is an optional directory where named plugins are searched
    /// for (e.g. `~/.vep/plugins/`).
    pub fn load_from_args(
        &mut self,
        plugin_args: &[String],
        plugin_dir: Option<&str>,
    ) -> Result<(), PluginError> {
        for arg in plugin_args {
            let parts: Vec<&str> = arg.splitn(2, ',').collect();
            let name_or_path = parts[0];
            let params: Vec<String> = if parts.len() > 1 {
                parts[1].split(',').map(String::from).collect()
            } else {
                Vec::new()
            };

            let path = resolve_plugin_path(name_or_path, plugin_dir)?;
            let mut plugin = LoadedPlugin::load(&path)?;

            if !params.is_empty() {
                plugin.init(&params)?;
            }

            tracing::info!(
                name = plugin.name(),
                version = plugin.version(),
                "loaded plugin"
            );
            self.plugins.push(plugin);
        }
        Ok(())
    }

    /// Returns the number of loaded plugins.
    pub fn plugin_count(&self) -> usize {
        self.plugins.len()
    }

    /// Returns `true` if no plugins are loaded.
    pub fn is_empty(&self) -> bool {
        self.plugins.is_empty()
    }

    /// Returns an iterator over the loaded plugins.
    pub fn plugins(&self) -> &[LoadedPlugin] {
        &self.plugins
    }

    /// Run all loaded plugins on the given consequence/variant pair (both as JSON).
    ///
    /// Returns a vec of `(plugin_name, result_json)` for each plugin.
    pub fn run_all(
        &self,
        consequence_json: &str,
        variant_json: &str,
    ) -> Result<Vec<(String, String)>, PluginError> {
        let mut results = Vec::with_capacity(self.plugins.len());
        for plugin in &self.plugins {
            let result = plugin.run(consequence_json, variant_json)?;
            results.push((plugin.name().to_owned(), result));
        }
        Ok(results)
    }
}

impl Default for PluginRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// Resolve a plugin name or path to a concrete file path.
///
/// If `name` is already an existing file path, it is returned directly.
/// Otherwise `plugin_dir` is searched for `lib{name}.so` and `lib{name}.dylib`.
fn resolve_plugin_path(name: &str, plugin_dir: Option<&str>) -> Result<PathBuf, PluginError> {
    let path = Path::new(name);
    if path.exists() {
        return Ok(path.to_path_buf());
    }

    if let Some(dir) = plugin_dir {
        let dir = Path::new(dir);
        for pattern in [format!("lib{name}.so"), format!("lib{name}.dylib")] {
            let candidate = dir.join(&pattern);
            if candidate.exists() {
                return Ok(candidate);
            }
        }
    }

    Err(PluginError::Load(format!("plugin '{name}' not found")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_registry_is_empty() {
        let registry = PluginRegistry::new();
        assert!(registry.is_empty());
        assert_eq!(registry.plugin_count(), 0);
    }

    #[test]
    fn default_registry_is_empty() {
        let registry = PluginRegistry::default();
        assert!(registry.is_empty());
    }

    #[test]
    fn resolve_nonexistent_plugin_without_dir() {
        let result = resolve_plugin_path("nonexistent_plugin", None);
        assert!(result.is_err());
        match result.unwrap_err() {
            PluginError::Load(msg) => {
                assert!(
                    msg.contains("nonexistent_plugin"),
                    "error should mention plugin name: {msg}"
                );
            }
            other => panic!("expected PluginError::Load, got: {other:?}"),
        }
    }

    #[test]
    fn resolve_nonexistent_plugin_with_dir() {
        let result = resolve_plugin_path("nonexistent_plugin", Some("/tmp/no_such_dir"));
        assert!(result.is_err());
    }

    #[test]
    fn resolve_existing_file_path() {
        // Use Cargo.toml as a stand-in for "an existing file"
        let manifest = env!("CARGO_MANIFEST_DIR");
        let cargo_toml = format!("{manifest}/Cargo.toml");
        let result = resolve_plugin_path(&cargo_toml, None);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().to_str().unwrap(), cargo_toml);
    }

    #[test]
    fn load_from_args_with_no_args() {
        let mut registry = PluginRegistry::new();
        let result = registry.load_from_args(&[], None);
        assert!(result.is_ok());
        assert!(registry.is_empty());
    }

    #[test]
    fn load_from_args_with_nonexistent_plugin() {
        let mut registry = PluginRegistry::new();
        let args = vec!["NoSuchPlugin".to_string()];
        let result = registry.load_from_args(&args, None);
        assert!(result.is_err());
    }

    #[test]
    fn plugins_slice_is_empty_initially() {
        let registry = PluginRegistry::new();
        assert!(registry.plugins().is_empty());
    }
}
