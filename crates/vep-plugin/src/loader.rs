// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Dynamic library loader for VEP plugins.
//!
//! [`LoadedPlugin`] wraps a `libloading::Library` and provides safe access to
//! the plugin's FFI functions.

use std::ffi::{c_char, c_void, CStr, CString};
use std::path::Path;

use libloading::{Library, Symbol};

use crate::abi::{self, ffi};
use crate::{FeatureType, PluginError};

/// A loaded plugin dylib with its instance pointer and function pointers.
///
/// The library handle is kept alive so that all symbol references remain valid.
/// When dropped, the plugin instance is destroyed via the plugin's destroy function
/// before the library is unloaded.
pub struct LoadedPlugin {
    _library: Library,
    instance: *mut c_void,
    name: String,
    version: String,
    destroy_fn: ffi::PluginDestroyFn,
    init_fn: ffi::PluginInitFn,
    run_fn: ffi::PluginRunFn,
    free_fn: ffi::PluginFreeFn,
    feature_types_fn: Option<ffi::PluginFeatureTypesFn>,
    feature_types: Vec<FeatureType>,
    declares_regulatory: bool,
}

impl std::fmt::Debug for LoadedPlugin {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LoadedPlugin")
            .field("name", &self.name)
            .field("version", &self.version)
            .field("feature_types", &self.feature_types)
            .finish_non_exhaustive()
    }
}

// The plugin FFI contract requires that the instance is Send + Sync.
// Plugin authors guarantee thread safety by implementing the C API accordingly.
unsafe impl Send for LoadedPlugin {}
unsafe impl Sync for LoadedPlugin {}

impl LoadedPlugin {
    /// Load a plugin from a dylib file at `path`.
    ///
    /// This checks the ABI version, creates an instance, and reads the plugin's
    /// name and version strings.
    pub fn load(path: &Path) -> Result<Self, PluginError> {
        unsafe {
            let lib = Library::new(path)
                .map_err(|e| PluginError::Load(format!("{}: {}", path.display(), e)))?;

            let abi_version_fn: Symbol<ffi::PluginAbiVersionFn> = lib
                .get(abi::symbol_names::ABI_VERSION)
                .map_err(|e| PluginError::Load(format!("missing vep_plugin_abi_version: {e}")))?;
            let abi_version = abi_version_fn();
            if abi_version != abi::VEP_PLUGIN_ABI_VERSION {
                return Err(PluginError::AbiMismatch {
                    expected: abi::VEP_PLUGIN_ABI_VERSION,
                    got: abi_version,
                });
            }

            let create_fn: Symbol<ffi::PluginCreateFn> = lib
                .get(abi::symbol_names::CREATE)
                .map_err(|e| PluginError::Load(format!("missing vep_plugin_create: {e}")))?;
            let destroy_fn: ffi::PluginDestroyFn = *lib
                .get::<ffi::PluginDestroyFn>(abi::symbol_names::DESTROY)
                .map_err(|e| PluginError::Load(format!("missing vep_plugin_destroy: {e}")))?;
            let name_fn: Symbol<ffi::PluginNameFn> = lib
                .get(abi::symbol_names::NAME)
                .map_err(|e| PluginError::Load(format!("missing vep_plugin_name: {e}")))?;
            let version_fn: Symbol<ffi::PluginVersionFn> = lib
                .get(abi::symbol_names::VERSION)
                .map_err(|e| PluginError::Load(format!("missing vep_plugin_version: {e}")))?;
            let init_fn: ffi::PluginInitFn = *lib
                .get::<ffi::PluginInitFn>(abi::symbol_names::INIT)
                .map_err(|e| PluginError::Load(format!("missing vep_plugin_init: {e}")))?;
            let run_fn: ffi::PluginRunFn = *lib
                .get::<ffi::PluginRunFn>(abi::symbol_names::RUN)
                .map_err(|e| PluginError::Load(format!("missing vep_plugin_run: {e}")))?;
            let free_fn: ffi::PluginFreeFn = *lib
                .get::<ffi::PluginFreeFn>(abi::symbol_names::FREE)
                .map_err(|e| PluginError::Load(format!("missing vep_plugin_free: {e}")))?;
            let feature_types_fn: Option<ffi::PluginFeatureTypesFn> = lib
                .get::<ffi::PluginFeatureTypesFn>(abi::symbol_names::FEATURE_TYPES)
                .ok()
                .map(|symbol| *symbol);

            let instance = create_fn();
            if instance.is_null() {
                return Err(PluginError::Init("plugin create returned null".into()));
            }

            let name = read_c_string(name_fn(instance), "name")?;
            let version = read_c_string(version_fn(instance), "version")?;

            Ok(Self {
                _library: lib,
                instance,
                name,
                version,
                destroy_fn,
                init_fn,
                run_fn,
                free_fn,
                feature_types_fn,
                feature_types: vec![FeatureType::Transcript],
                declares_regulatory: false,
            })
        }
    }

    /// Reads the feature types the plugin declares, after its parameters have
    /// been given; a plugin without the symbol keeps `Transcript` alone.
    pub fn read_feature_types(&mut self) -> Result<(), PluginError> {
        let Some(feature_types_fn) = self.feature_types_fn else {
            return Ok(());
        };
        let text = unsafe { read_c_string(feature_types_fn(self.instance), "feature types")? };
        self.feature_types = FeatureType::parse_list(&text).map_err(|e| match e {
            PluginError::Load(reason) => {
                PluginError::Load(format!("plugin '{}': {reason}", self.name))
            }
            other => other,
        })?;
        // Ensembl VEP reads the declared words, not the classes they cover
        // (Runner.pm:998 greps them for motif or regulatory), so `Feature`
        // alone switches nothing on.
        self.declares_regulatory = text
            .split(',')
            .map(str::trim)
            .any(|word| matches!(word, "RegulatoryFeature" | "MotifFeature"));
        Ok(())
    }

    /// Whether the plugin's declared list names `RegulatoryFeature` or
    /// `MotifFeature` by word, which switches regulatory annotation on for the
    /// run (`Runner.pm:998-1002`); `Feature`, which covers both, does not.
    pub fn declares_regulatory_features(&self) -> bool {
        self.declares_regulatory
    }

    /// The feature types the plugin annotates ([`Self::read_feature_types`]).
    pub fn feature_types(&self) -> &[FeatureType] {
        &self.feature_types
    }

    /// Whether the plugin runs on a row of `kind`, as Ensembl VEP runs a plugin
    /// only on the feature classes it declares (`OutputFactory.pm` `run_plugins`
    /// through `check_feature_type`).
    pub fn annotates(&self, kind: vep_core::consequence::FeatureType) -> bool {
        self.feature_types.iter().any(|t| t.row_kind() == kind)
    }

    /// Initialize the plugin with the given parameters.
    pub fn init(&mut self, params: &[String]) -> Result<(), PluginError> {
        let c_params: Vec<CString> = params
            .iter()
            .map(|p| {
                CString::new(p.as_str())
                    .map_err(|e| PluginError::Init(format!("invalid parameter string: {e}")))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let c_param_ptrs: Vec<*const c_char> = c_params.iter().map(|s| s.as_ptr()).collect();

        let rc =
            unsafe { (self.init_fn)(self.instance, c_param_ptrs.as_ptr(), c_param_ptrs.len()) };
        if rc != 0 {
            return Err(PluginError::Init(format!(
                "plugin '{}' init returned error code {rc}",
                self.name
            )));
        }
        Ok(())
    }

    /// Run the plugin on a consequence/variant pair, both serialized as JSON.
    ///
    /// Returns the plugin's result as a JSON string.
    pub fn run(&self, consequence_json: &str, variant_json: &str) -> Result<String, PluginError> {
        let c_consequence = CString::new(consequence_json)
            .map_err(|e| PluginError::Run(format!("invalid consequence JSON: {e}")))?;
        let c_variant = CString::new(variant_json)
            .map_err(|e| PluginError::Run(format!("invalid variant JSON: {e}")))?;

        unsafe {
            let result_ptr =
                (self.run_fn)(self.instance, c_consequence.as_ptr(), c_variant.as_ptr());
            if result_ptr.is_null() {
                return Err(PluginError::Run(format!(
                    "plugin '{}' run returned null",
                    self.name
                )));
            }

            let result = CStr::from_ptr(result_ptr)
                .to_str()
                .map_err(|e| PluginError::Run(format!("invalid UTF-8 in plugin result: {e}")))?
                .to_owned();

            (self.free_fn)(result_ptr);

            Ok(result)
        }
    }

    /// Returns the plugin name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the plugin version.
    pub fn version(&self) -> &str {
        &self.version
    }
}

impl Drop for LoadedPlugin {
    fn drop(&mut self) {
        if !self.instance.is_null() {
            unsafe {
                (self.destroy_fn)(self.instance);
            }
            self.instance = std::ptr::null_mut();
        }
    }
}

/// Read a null-terminated C string from a plugin, returning a Rust `String`.
unsafe fn read_c_string(ptr: *const c_char, field: &str) -> Result<String, PluginError> {
    if ptr.is_null() {
        return Err(PluginError::Load(format!(
            "plugin returned null for {field}"
        )));
    }
    CStr::from_ptr(ptr)
        .to_str()
        .map(|s| s.to_owned())
        .map_err(|e| PluginError::Load(format!("invalid UTF-8 in plugin {field}: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn load_nonexistent_dylib_returns_error() {
        let path = PathBuf::from("/nonexistent/plugin.so");
        let result = LoadedPlugin::load(&path);
        assert!(result.is_err());
        match result.unwrap_err() {
            PluginError::Load(msg) => {
                assert!(
                    msg.contains("/nonexistent/plugin.so"),
                    "error should contain path: {msg}"
                );
            }
            other => panic!("expected PluginError::Load, got: {other:?}"),
        }
    }
}
