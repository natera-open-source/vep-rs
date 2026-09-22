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
pub const VEP_PLUGIN_ABI_VERSION: u32 = 1;

/// Symbol names used when looking up exported functions from a plugin dylib.
pub mod symbol_names {
    pub const ABI_VERSION: &[u8] = b"vep_plugin_abi_version\0";
    pub const CREATE: &[u8] = b"vep_plugin_create\0";
    pub const DESTROY: &[u8] = b"vep_plugin_destroy\0";
    pub const NAME: &[u8] = b"vep_plugin_name\0";
    pub const VERSION: &[u8] = b"vep_plugin_version\0";
    pub const INIT: &[u8] = b"vep_plugin_init\0";
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
    pub type PluginInitFn = unsafe extern "C" fn(*mut c_void, *const *const c_char, usize) -> i32;

    /// Run the plugin on a consequence.
    ///
    /// - First `*const c_void` is the plugin instance.
    /// - Second `*const c_char` is the consequence as a JSON string.
    /// - Third `*const c_char` is the variant as a JSON string.
    /// - Returns a JSON string that the caller must free via [`PluginFreeFn`].
    pub type PluginRunFn =
        unsafe extern "C" fn(*const c_void, *const c_char, *const c_char) -> *mut c_char;

    /// Free a string returned by [`PluginRunFn`].
    pub type PluginFreeFn = unsafe extern "C" fn(*mut c_char);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn abi_version_is_defined() {
        assert_eq!(VEP_PLUGIN_ABI_VERSION, 1);
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
