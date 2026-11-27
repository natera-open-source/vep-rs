// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! C ABI interface that plugin dylibs must export.
//!
//! Plugin authors compile their plugin as a `cdylib` crate and export functions
//! matching the signatures defined in [`ffi`]. The loader looks up these symbols
//! by name via `libloading`.

/// ABI version for plugin compatibility checking.
///
/// Bumped whenever the FFI contract changes in a backwards-incompatible way.
/// Version 2 adds the required `vep_plugin_header_info` symbol and the optional
/// `vep_plugin_feature_types` symbol, and calls `vep_plugin_init` on every load.
pub const VEP_PLUGIN_ABI_VERSION: u32 = 2;

/// Symbol names used when looking up exported functions from a plugin dylib.
pub mod symbol_names {
    pub const ABI_VERSION: &[u8] = b"vep_plugin_abi_version\0";
    pub const CREATE: &[u8] = b"vep_plugin_create\0";
    pub const DESTROY: &[u8] = b"vep_plugin_destroy\0";
    pub const NAME: &[u8] = b"vep_plugin_name\0";
    pub const VERSION: &[u8] = b"vep_plugin_version\0";
    pub const INIT: &[u8] = b"vep_plugin_init\0";
    pub const HEADER_INFO: &[u8] = b"vep_plugin_header_info\0";
    pub const FEATURE_TYPES: &[u8] = b"vep_plugin_feature_types\0";
    pub const RUN: &[u8] = b"vep_plugin_run\0";
    pub const FREE: &[u8] = b"vep_plugin_free\0";
}

/// C function signatures that a plugin dylib must export.
///
/// These are looked up by name via `libloading`.
pub mod ffi {
    use std::ffi::{c_char, c_void};

    /// Returns the ABI version the plugin was compiled against.
    pub type PluginAbiVersionFn = unsafe extern "C" fn() -> u32;

    /// Creates a new plugin instance. Returns an opaque pointer.
    pub type PluginCreateFn = unsafe extern "C" fn() -> *mut c_void;

    /// Destroys a plugin instance created by [`PluginCreateFn`].
    pub type PluginDestroyFn = unsafe extern "C" fn(*mut c_void);

    /// Returns the plugin name as a null-terminated C string.
    /// The returned pointer must remain valid for the lifetime of the instance.
    pub type PluginNameFn = unsafe extern "C" fn(*const c_void) -> *const c_char;

    /// Returns the plugin version as a null-terminated C string.
    /// The returned pointer must remain valid for the lifetime of the instance.
    pub type PluginVersionFn = unsafe extern "C" fn(*const c_void) -> *const c_char;

    /// Initialize the plugin with parameters. Returns 0 on success, non-zero on error.
    /// Called exactly once per instance, with the `--plugin` parameters after the
    /// name, zero or more, before any other instance call but `vep_plugin_name`
    /// and `vep_plugin_version`. At a count of 0 the parameter pointer is not
    /// null but must not be dereferenced.
    pub type PluginInitFn = unsafe extern "C" fn(*mut c_void, *const *const c_char, usize) -> i32;

    /// Returns the fields the plugin writes, as a JSON array of `[key, description]`
    /// pairs (`[]` or `null` for none), which the headers print and the column
    /// formats take their plugin columns from. The caller frees the string via
    /// [`PluginFreeFn`].
    pub type PluginHeaderInfoFn = unsafe extern "C" fn(*const c_void) -> *mut c_char;

    /// Returns the feature types the plugin runs on, as a JSON array of
    /// `"Transcript"`, `"RegulatoryFeature"`, `"MotifFeature"` and `"Intergenic"`.
    /// Optional: a library that does not export it runs on `Transcript` rows alone, Ensembl VEP's
    /// default (`BaseVepPlugin.pm:99`). The caller frees the string via [`PluginFreeFn`].
    pub type PluginFeatureTypesFn = unsafe extern "C" fn(*const c_void) -> *mut c_char;

    /// Run the plugin on a consequence.
    ///
    /// - First `*const c_void` is the plugin instance.
    /// - Second `*const c_char` is the consequence as a JSON string, the JSON text
    ///   `null` (never a null pointer) on an intergenic row.
    /// - Third `*const c_char` is the variant as a JSON string.
    /// - Returns a JSON string that the caller must free via [`PluginFreeFn`].
    ///
    /// The host may call it on one instance from several threads at once, so
    /// the plugin must be thread-safe.
    pub type PluginRunFn =
        unsafe extern "C" fn(*const c_void, *const c_char, *const c_char) -> *mut c_char;

    /// Free a string returned by [`PluginRunFn`], [`PluginHeaderInfoFn`] or
    /// [`PluginFeatureTypesFn`].
    pub type PluginFreeFn = unsafe extern "C" fn(*mut c_char);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn abi_version_is_defined() {
        assert_eq!(VEP_PLUGIN_ABI_VERSION, 2);
    }

    #[test]
    fn symbol_names_are_null_terminated() {
        // All symbol name byte slices must end with \0 for libloading
        for name in [
            symbol_names::ABI_VERSION,
            symbol_names::CREATE,
            symbol_names::DESTROY,
            symbol_names::NAME,
            symbol_names::VERSION,
            symbol_names::INIT,
            symbol_names::HEADER_INFO,
            symbol_names::FEATURE_TYPES,
            symbol_names::RUN,
            symbol_names::FREE,
        ] {
            assert_eq!(
                name.last().copied(),
                Some(0),
                "symbol name must be null-terminated"
            );
        }
    }
}
