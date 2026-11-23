// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The plugin interface, held as literals of this file: the interface number,
//! the ten symbol names the loader looks up and the C signature of each symbol.
//!
//! A plugin library is built against these values and never against this crate,
//! so a change to any of them breaks every library already built. The tests
//! compare them with the public items of `vep_plugin::abi`.

use std::ffi::{c_char, c_void};

use vep_plugin::abi::{ffi, symbol_names, VEP_PLUGIN_ABI_VERSION};

/// The number a plugin library returns from `vep_plugin_abi_version`.
const INTERFACE_NUMBER: u32 = 2;

/// The ten names a plugin library exports, each as the loader hands it to the
/// dynamic loader: the bytes of the C name, then the terminating NUL.
const EXPORTED_NAMES: [&[u8]; 10] = [
    b"vep_plugin_abi_version\0",
    b"vep_plugin_create\0",
    b"vep_plugin_destroy\0",
    b"vep_plugin_name\0",
    b"vep_plugin_version\0",
    b"vep_plugin_init\0",
    b"vep_plugin_header_info\0",
    b"vep_plugin_feature_types\0",
    b"vep_plugin_run\0",
    b"vep_plugin_free\0",
];

/// Implemented for a type by that type alone, so a bound `A: SameAs<B>` holds
/// only where `A` and `B` are one type; no coercion between function pointers
/// satisfies it.
trait SameAs<T> {}

impl<T> SameAs<T> for T {}

fn same_type<A: SameAs<B>, B>() {}

/// The loader refuses a library whose `vep_plugin_abi_version` returns any other
/// number.
#[test]
fn the_interface_number_is_two() {
    assert_eq!(VEP_PLUGIN_ABI_VERSION, INTERFACE_NUMBER);
}

#[test]
fn the_loader_looks_up_the_ten_names_a_plugin_library_exports() {
    let looked_up: [&[u8]; 10] = [
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
    ];
    for (looked_up, exported) in looked_up.into_iter().zip(EXPORTED_NAMES) {
        assert_eq!(
            looked_up,
            exported,
            "the loader looks up {:?} where a plugin library exports {:?}",
            String::from_utf8_lossy(looked_up),
            String::from_utf8_lossy(exported)
        );
    }
}

/// Each function-pointer type of `abi::ffi` beside the C signature a plugin
/// library exports under that symbol. The compiler makes the comparison: a
/// signature that differs stops the build of this file.
#[test]
fn each_symbol_keeps_its_c_signature() {
    same_type::<ffi::PluginAbiVersionFn, unsafe extern "C" fn() -> u32>();
    same_type::<ffi::PluginCreateFn, unsafe extern "C" fn() -> *mut c_void>();
    same_type::<ffi::PluginDestroyFn, unsafe extern "C" fn(*mut c_void)>();
    same_type::<ffi::PluginNameFn, unsafe extern "C" fn(*const c_void) -> *const c_char>();
    same_type::<ffi::PluginVersionFn, unsafe extern "C" fn(*const c_void) -> *const c_char>();
    same_type::<
        ffi::PluginInitFn,
        unsafe extern "C" fn(*mut c_void, *const *const c_char, usize) -> i32,
    >();
    same_type::<ffi::PluginHeaderInfoFn, unsafe extern "C" fn(*const c_void) -> *mut c_char>();
    same_type::<ffi::PluginFeatureTypesFn, unsafe extern "C" fn(*const c_void) -> *mut c_char>();
    same_type::<
        ffi::PluginRunFn,
        unsafe extern "C" fn(*const c_void, *const c_char, *const c_char) -> *mut c_char,
    >();
    same_type::<ffi::PluginFreeFn, unsafe extern "C" fn(*mut c_char)>();
}
