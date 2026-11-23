// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Plugin registry: discovers, loads, and manages plugin dylibs.

use std::ffi::CString;
use std::path::{Path, PathBuf};

use crate::loader::LoadedPlugin;
use crate::{FeatureType, PluginError};
use vep_core::sv_routing::split_plugin_arg;

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
    /// is the plugin name or path; the rest are init parameters, read as Ensembl
    /// VEP reads the list ([`split_plugin_arg`]) and passed to `vep_plugin_init`
    /// on every load, an empty list included.
    ///
    /// `plugin_dirs` are the directories a name is searched for in, in order,
    /// the first holding the library winning.
    pub fn load_from_args(
        &mut self,
        plugin_args: &[String],
        plugin_dirs: &[String],
    ) -> Result<(), PluginError> {
        for arg in plugin_args {
            let (name_or_path, params) = split_plugin_arg(arg);

            let path = resolve_plugin_path(name_or_path, plugin_dirs)?;
            let plugin = LoadedPlugin::load(&path, &params)?;

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

    /// Returns the loaded plugins as a slice, in load order.
    pub fn plugins(&self) -> &[LoadedPlugin] {
        &self.plugins
    }

    /// Returns header info from all loaded plugins: plugin by plugin in load
    /// order, each plugin's fields sorted by name, the order Ensembl VEP's
    /// `OutputFactory.pm` `get_plugin_headers` gives its plugin columns.
    pub fn all_header_info(&self) -> Vec<(String, String)> {
        self.plugins
            .iter()
            .flat_map(|plugin| plugin.header_info().iter().cloned())
            .collect()
    }

    /// Whether any loaded plugin runs on a row of kind `row`
    /// ([`LoadedPlugin::runs_on`]).
    pub fn runs_on(&self, row: FeatureType) -> bool {
        self.plugins.iter().any(|plugin| plugin.runs_on(row))
    }

    /// The names of the loaded plugins that list `RegulatoryFeature` or
    /// `MotifFeature`, each of which switches regulatory annotation on for the
    /// run, as Ensembl VEP's `Runner.pm:998-1002` does for a plugin whose
    /// declared feature types name motif or regulatory features.
    pub fn plugins_wanting_regulatory_features(&self) -> Vec<&str> {
        self.plugins
            .iter()
            .filter(|plugin| {
                let types = plugin.feature_types();
                types.contains(&FeatureType::RegulatoryFeature)
                    || types.contains(&FeatureType::MotifFeature)
            })
            .map(|plugin| plugin.name())
            .collect()
    }

    /// Each loaded plugin's fields with the kinds of row the plugin runs on
    /// ([`LoadedPlugin::feature_types`]), for the writers that print a field on
    /// the rows of those kinds alone.
    pub fn field_row_kinds(&self) -> Vec<(String, Vec<FeatureType>)> {
        self.plugins
            .iter()
            .flat_map(|plugin| {
                plugin
                    .header_info()
                    .iter()
                    .map(|(name, _)| (name.clone(), plugin.feature_types().to_vec()))
            })
            .collect()
    }

    /// Run all loaded plugins, whatever their feature types, on the given
    /// consequence/variant pair (both as JSON).
    ///
    /// Returns a vec of `(plugin_name, result_json)` for each plugin.
    pub fn run_all(
        &self,
        consequence_json: &str,
        variant_json: &str,
    ) -> Result<Vec<(String, String)>, PluginError> {
        let consequence_json = CString::new(consequence_json)
            .map_err(|e| PluginError::Run(format!("invalid consequence JSON: {e}")))?;
        let variant_json = CString::new(variant_json)
            .map_err(|e| PluginError::Run(format!("invalid variant JSON: {e}")))?;
        let mut results = Vec::with_capacity(self.plugins.len());
        for plugin in &self.plugins {
            let result = plugin.run(&consequence_json, &variant_json)?;
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

/// Whether a `--plugin` token is written as a file, a path holding a separator
/// or a file name carrying a suffix (`libecho.so`, `Name.pm`), rather than as a
/// plugin's bare name. Ensembl VEP's `use $module` (Runner.pm:931-933) compiles
/// a bare word alone, so only a bare word can name an Ensembl VEP plugin.
pub fn names_a_file(token: &str) -> bool {
    has_separator(token) || Path::new(token).extension().is_some()
}

fn has_separator(token: &str) -> bool {
    token.chars().any(std::path::is_separator)
}

/// Resolve a plugin name or path to a library file.
///
/// A token written as a file ([`names_a_file`]) that exists is the plugin file:
/// a path as given, a file name in the working directory as `./<name>`. A bare
/// name is the library's stem: each of `plugin_dirs` in turn is searched for
/// `lib{name}.so` and `lib{name}.dylib`, and the first regular file found is the
/// plugin. The working directory is never searched for a bare name: Ensembl VEP
/// finds the module through `@INC`, onto which only the plugins directories are
/// pushed (Runner.pm:921, `:931-933`), so an entry named like the plugin in the
/// working directory is not a plugin.
fn resolve_plugin_path(name: &str, plugin_dirs: &[String]) -> Result<PathBuf, PluginError> {
    if names_a_file(name) {
        // The dynamic loader reads a name without a separator against its own
        // search path, which on Linux never holds the working directory, so a
        // file name is handed over with its directory.
        let path = if has_separator(name) {
            PathBuf::from(name)
        } else {
            Path::new(".").join(name)
        };
        if path.exists() {
            return Ok(path);
        }
    }

    for dir in plugin_dirs {
        let dir = Path::new(dir);
        for pattern in [format!("lib{name}.so"), format!("lib{name}.dylib")] {
            let candidate = dir.join(&pattern);
            if candidate.is_file() {
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
        let result = resolve_plugin_path("nonexistent_plugin", &[]);
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
        let result = resolve_plugin_path("nonexistent_plugin", &["/tmp/no_such_dir".to_string()]);
        assert!(result.is_err());
    }

    #[test]
    fn resolve_existing_file_path() {
        // Use Cargo.toml as a stand-in for "an existing file"
        let manifest = env!("CARGO_MANIFEST_DIR");
        let cargo_toml = format!("{manifest}/Cargo.toml");
        let result = resolve_plugin_path(&cargo_toml, &[]);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().to_str().unwrap(), cargo_toml);
    }

    /// A token holding a path separator or a suffix is written as a file; a bare
    /// word is a plugin name, the only token Ensembl VEP 116.3's `use $module`
    /// compiles (Runner.pm:931-933).
    #[test]
    fn names_a_file_tells_a_path_or_a_suffixed_name_from_a_bare_name() {
        for token in [
            "./GeneSplicer",
            "/opt/plugins/libecho.so",
            "plugins/x",
            "libecho.so",
            "Foo.pm",
        ] {
            assert!(names_a_file(token), "{token}");
        }
        for token in ["GeneSplicer", "echo_plugin", "dbNSFP", "FATHMM_MKL"] {
            assert!(!names_a_file(token), "{token}");
        }
    }

    /// A bare name is never read against the working directory, where Ensembl VEP
    /// 116.3 never looks for a module (`@INC` holds the plugins directories
    /// alone, Runner.pm:921): the crate's own `src` directory, an entry of the
    /// test process's working directory, resolves to nothing, while a file name
    /// written with its suffix beside it resolves to that file under `./`, the
    /// form the dynamic loader opens as a file.
    #[test]
    fn a_bare_name_matching_a_working_directory_entry_is_not_found() {
        assert!(
            Path::new("src").is_dir() && Path::new("Cargo.toml").is_file(),
            "the test runs in the crate directory"
        );
        let err = resolve_plugin_path("src", &[]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "plugin load failed: plugin 'src' not found"
        );
        assert_eq!(
            resolve_plugin_path("Cargo.toml", &[]).unwrap(),
            PathBuf::from("./Cargo.toml")
        );
        assert_eq!(
            resolve_plugin_path("./Cargo.toml", &[]).unwrap(),
            PathBuf::from("./Cargo.toml")
        );
    }

    #[test]
    fn load_from_args_with_no_args() {
        let mut registry = PluginRegistry::new();
        let result = registry.load_from_args(&[], &[]);
        assert!(result.is_ok());
        assert!(registry.is_empty());
    }

    #[test]
    fn load_from_args_with_nonexistent_plugin() {
        let mut registry = PluginRegistry::new();
        let args = vec!["NoSuchPlugin".to_string()];
        let result = registry.load_from_args(&args, &[]);
        assert!(result.is_err());
    }

    #[test]
    fn plugins_slice_is_empty_initially() {
        let registry = PluginRegistry::new();
        assert!(registry.plugins().is_empty());
    }
}
