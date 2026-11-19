// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Dynamic library loader for VEP plugins.
//!
//! [`LoadedPlugin`] wraps a `libloading::Library` and provides safe access to
//! the plugin's FFI functions.

use std::collections::BTreeMap;
use std::ffi::{c_char, c_void, CStr, CString};
use std::path::Path;

use libloading::Library;

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
    header_info: Vec<(String, String)>,
    feature_types: Option<Vec<FeatureType>>,
    destroy_fn: ffi::PluginDestroyFn,
    run_fn: ffi::PluginRunFn,
    free_fn: ffi::PluginFreeFn,
}

impl std::fmt::Debug for LoadedPlugin {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LoadedPlugin")
            .field("name", &self.name)
            .field("version", &self.version)
            .field("header_info", &self.header_info)
            .field("feature_types", &self.feature_types)
            .finish_non_exhaustive()
    }
}

// The plugin FFI contract requires that the instance is Send + Sync.
// Plugin authors guarantee thread safety by implementing the C API accordingly.
unsafe impl Send for LoadedPlugin {}
unsafe impl Sync for LoadedPlugin {}

impl LoadedPlugin {
    /// Load a plugin from a dylib file at `path` and initialise it with `params`,
    /// the `--plugin` parameters after the name.
    ///
    /// This checks the ABI version, resolves the symbols, creates an instance,
    /// reads its name and version, calls `vep_plugin_init` once (with zero or
    /// more parameters), then reads the header info and, when the library
    /// exports it, the feature types.
    pub fn load(path: &Path, params: &[String]) -> Result<Self, PluginError> {
        unsafe {
            let lib = Library::new(path).map_err(|e| {
                PluginError::Load(format!("{}: {}", path.display(), loader_message(&e)))
            })?;

            let abi_version_fn =
                required::<ffi::PluginAbiVersionFn>(&lib, abi::symbol_names::ABI_VERSION)?;
            let abi_version = abi_version_fn();
            if abi_version != abi::VEP_PLUGIN_ABI_VERSION {
                return Err(PluginError::AbiMismatch {
                    expected: abi::VEP_PLUGIN_ABI_VERSION,
                    got: abi_version,
                });
            }

            let create_fn = required::<ffi::PluginCreateFn>(&lib, abi::symbol_names::CREATE)?;
            let destroy_fn = required::<ffi::PluginDestroyFn>(&lib, abi::symbol_names::DESTROY)?;
            let name_fn = required::<ffi::PluginNameFn>(&lib, abi::symbol_names::NAME)?;
            let version_fn = required::<ffi::PluginVersionFn>(&lib, abi::symbol_names::VERSION)?;
            let init_fn = required::<ffi::PluginInitFn>(&lib, abi::symbol_names::INIT)?;
            let header_info_fn =
                required::<ffi::PluginHeaderInfoFn>(&lib, abi::symbol_names::HEADER_INFO)?;
            let feature_types_fn = lib
                .get::<ffi::PluginFeatureTypesFn>(abi::symbol_names::FEATURE_TYPES)
                .ok()
                .map(|symbol| *symbol);
            let run_fn = required::<ffi::PluginRunFn>(&lib, abi::symbol_names::RUN)?;
            let free_fn = required::<ffi::PluginFreeFn>(&lib, abi::symbol_names::FREE)?;

            let instance = create_fn();
            if instance.is_null() {
                return Err(PluginError::Init("plugin create returned null".into()));
            }
            // From here the instance is owned, so an error below destroys it.
            let mut plugin = Self {
                _library: lib,
                instance,
                name: String::new(),
                version: String::new(),
                header_info: Vec::new(),
                feature_types: None,
                destroy_fn,
                run_fn,
                free_fn,
            };
            plugin.name = read_c_string(name_fn(instance), "name")?;
            plugin.version = read_c_string(version_fn(instance), "version")?;
            plugin.init(init_fn, params)?;
            let header_info_text = plugin.take_text(header_info_fn(instance), "header info")?;
            plugin.header_info = parse_header_info(&header_info_text)?;
            if let Some(feature_types_fn) = feature_types_fn {
                let feature_types_text =
                    plugin.take_text(feature_types_fn(instance), "feature types")?;
                plugin.feature_types = Some(plugin.parse_feature_types(&feature_types_text)?);
            }
            Ok(plugin)
        }
    }

    /// Calls `vep_plugin_init` with `params`.
    fn init(&self, init_fn: ffi::PluginInitFn, params: &[String]) -> Result<(), PluginError> {
        let c_params: Vec<CString> = params
            .iter()
            .map(|p| {
                CString::new(p.as_str())
                    .map_err(|e| PluginError::Init(format!("invalid parameter string: {e}")))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let c_param_ptrs: Vec<*const c_char> = c_params.iter().map(|s| s.as_ptr()).collect();

        let rc = unsafe { init_fn(self.instance, c_param_ptrs.as_ptr(), c_param_ptrs.len()) };
        if rc != 0 {
            return Err(PluginError::Init(format!(
                "plugin '{}' init returned error code {rc}",
                self.name
            )));
        }
        Ok(())
    }

    /// Copies a string the plugin allocated for `what` and releases it through
    /// `vep_plugin_free`.
    ///
    /// # Safety
    /// `ptr` is null or a string the plugin returned and has not freed.
    unsafe fn take_text(&self, ptr: *mut c_char, what: &str) -> Result<String, PluginError> {
        if ptr.is_null() {
            return Err(PluginError::Load(format!(
                "plugin '{}' returned null for {what}",
                self.name
            )));
        }
        let text = CStr::from_ptr(ptr)
            .to_str()
            .map(str::to_owned)
            .map_err(|e| {
                PluginError::Load(format!(
                    "invalid UTF-8 in plugin '{}' {what}: {e}",
                    self.name
                ))
            });
        (self.free_fn)(ptr);
        text
    }

    /// The feature types a plugin lists, each a name [`FeatureType::parse`] knows.
    fn parse_feature_types(&self, text: &str) -> Result<Vec<FeatureType>, PluginError> {
        let names: Vec<String> = serde_json::from_str(text).map_err(|e| {
            PluginError::Load(format!(
                "plugin '{}' feature types are not a JSON array of strings: {e}",
                self.name
            ))
        })?;
        names
            .iter()
            .map(|name| {
                FeatureType::parse(name).ok_or_else(|| {
                    PluginError::Load(format!(
                        "plugin '{}' lists an unknown feature type '{name}'",
                        self.name
                    ))
                })
            })
            .collect()
    }

    /// Run the plugin on a consequence/variant pair, both serialized as JSON.
    ///
    /// `consequence_json` is `null` on an intergenic row. Returns the plugin's
    /// result as a JSON string.
    pub fn run(&self, consequence_json: &CStr, variant_json: &CStr) -> Result<String, PluginError> {
        unsafe {
            let result_ptr = (self.run_fn)(
                self.instance,
                consequence_json.as_ptr(),
                variant_json.as_ptr(),
            );
            if result_ptr.is_null() {
                return Err(PluginError::Run(format!(
                    "plugin '{}' run returned null",
                    self.name
                )));
            }

            let result = CStr::from_ptr(result_ptr)
                .to_str()
                .map_err(|e| PluginError::Run(format!("invalid UTF-8 in plugin result: {e}")))
                .map(str::to_owned);

            (self.free_fn)(result_ptr);

            result
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

    /// The fields the plugin writes, as `(key, description)` pairs sorted by key,
    /// the order Ensembl VEP's `get_plugin_headers` gives one plugin's keys.
    pub fn header_info(&self) -> &[(String, String)] {
        &self.header_info
    }

    /// The feature types the plugin runs on: the kinds it listed, or
    /// `Transcript` alone when the library exports no `vep_plugin_feature_types`,
    /// Ensembl VEP's default for a plugin that does not override `feature_types`
    /// (`BaseVepPlugin.pm:99`).
    pub fn feature_types(&self) -> &[FeatureType] {
        self.feature_types
            .as_deref()
            .unwrap_or(&FeatureType::DEFAULT)
    }

    /// Whether the plugin runs on a row of kind `row` (Ensembl VEP's
    /// `check_feature_type`, `BaseVepPlugin.pm:310-313`).
    pub fn runs_on(&self, row: FeatureType) -> bool {
        self.feature_types().contains(&row)
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

/// The function pointer behind a required symbol; a missing one is a load error
/// naming the symbol.
///
/// # Safety
/// `T` is the function-pointer type the symbol was exported with.
unsafe fn required<T: Copy>(lib: &Library, symbol: &[u8]) -> Result<T, PluginError> {
    lib.get::<T>(symbol).map(|function| *function).map_err(|e| {
        let name = String::from_utf8_lossy(symbol);
        PluginError::Load(format!(
            "missing {}: {}",
            name.trim_end_matches('\0'),
            loader_message(&e)
        ))
    })
}

/// The text of a library error: the dynamic loader's own message (dlopen's or
/// dlsym's) where `libloading` carries one as the error's source, else the
/// library's summary such as `dlopen failed`.
fn loader_message(error: &libloading::Error) -> String {
    match std::error::Error::source(error) {
        Some(loader_error) => loader_error.to_string(),
        None => error.to_string(),
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

/// The header info a plugin returned: a JSON array of `[key, description]` pairs,
/// or `null`, read into pairs sorted by key with the last description of a
/// repeated key kept.
fn parse_header_info(text: &str) -> Result<Vec<(String, String)>, PluginError> {
    let pairs: Option<Vec<(String, String)>> = serde_json::from_str(text).map_err(|e| {
        PluginError::Load(format!(
            "plugin header info is not a JSON array of [key, description] pairs: {e}"
        ))
    })?;
    let sorted: BTreeMap<String, String> = pairs.unwrap_or_default().into_iter().collect();
    Ok(sorted.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn load_nonexistent_dylib_returns_error() {
        let path = PathBuf::from("/nonexistent/plugin.so");
        let result = LoadedPlugin::load(&path, &[]);
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

    /// The dynamic loader's own words reach the message: "cannot open shared
    /// object file: No such file or directory" from dlopen on Linux, "no such
    /// file" on macOS; the library's summary alone is "dlopen failed".
    #[test]
    fn a_missing_library_is_reported_with_the_loaders_reason() {
        let path = PathBuf::from("/nonexistent/plugin.so");
        let PluginError::Load(msg) = LoadedPlugin::load(&path, &[]).unwrap_err() else {
            panic!("expected PluginError::Load");
        };
        let lower = msg.to_ascii_lowercase();
        assert!(
            lower.contains("no such file") || lower.contains("cannot open"),
            "{msg}"
        );
    }

    #[test]
    fn header_info_is_sorted_by_key_and_null_is_none() {
        let pairs = parse_header_info(r#"[["B_KEY","second"],["A_KEY","first"]]"#).unwrap();
        assert_eq!(
            pairs,
            vec![
                ("A_KEY".to_string(), "first".to_string()),
                ("B_KEY".to_string(), "second".to_string()),
            ]
        );
        assert!(parse_header_info("null").unwrap().is_empty());
        assert!(parse_header_info("[]").unwrap().is_empty());
        let err = parse_header_info(r#"{"A_KEY":"first"}"#).unwrap_err();
        assert!(matches!(err, PluginError::Load(msg) if msg.contains("[key, description]")));
    }
}
