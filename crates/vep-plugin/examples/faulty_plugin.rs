// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Example C-ABI plugin that fails on purpose, the way a defective plugin would:
//! its first `run` returns `[1]`, a JSON value that is not an object, its second
//! returns NULL, and every call after those writes `FAULTY_OK=1`. The host warns
//! once about each failure, keeps the row and goes on, so a run with this plugin
//! loaded exercises that contract; a plugin that fails the same way on every row
//! draws one warning, since the host writes a warning text once.
//!
//! Built as a `cdylib`. The plugin lists `Transcript` as its one feature type,
//! so it runs on transcript rows alone. The host may run batches on several
//! threads, so the call counter is atomic and the two failures are the first
//! two calls in time, not the first two rows of the output. Every string handed
//! to the host is a `CString::into_raw`, which `vep_plugin_free` reclaims with
//! `CString::from_raw`.

use std::ffi::{c_char, c_void, CString};
use std::ptr;
use std::sync::atomic::{AtomicU64, Ordering};

struct FaultyPlugin {
    name: CString,
    version: CString,
    runs: AtomicU64,
}

#[no_mangle]
pub extern "C" fn vep_plugin_abi_version() -> u32 {
    2
}

#[no_mangle]
pub extern "C" fn vep_plugin_create() -> *mut c_void {
    let plugin = FaultyPlugin {
        name: CString::new("faulty").expect("no interior NUL"),
        version: CString::new("0.1.0").expect("no interior NUL"),
        runs: AtomicU64::new(0),
    };
    Box::into_raw(Box::new(plugin)).cast()
}

/// # Safety
/// `instance` must be a pointer returned by [`vep_plugin_create`] that has not
/// been destroyed.
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_destroy(instance: *mut c_void) {
    if !instance.is_null() {
        drop(Box::from_raw(instance.cast::<FaultyPlugin>()));
    }
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_name(instance: *const c_void) -> *const c_char {
    (*instance.cast::<FaultyPlugin>()).name.as_ptr()
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_version(instance: *const c_void) -> *const c_char {
    (*instance.cast::<FaultyPlugin>()).version.as_ptr()
}

/// The plugin takes no parameters, so every call succeeds.
///
/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_init(
    _instance: *mut c_void,
    _params: *const *const c_char,
    _count: usize,
) -> i32 {
    0
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_header_info(_instance: *const c_void) -> *mut c_char {
    CString::new(r#"[["FAULTY_OK","1 on every row after the plugin's two failing calls"]]"#)
        .expect("no interior NUL")
        .into_raw()
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_feature_types(_instance: *const c_void) -> *mut c_char {
    CString::new(r#"["Transcript"]"#)
        .expect("no interior NUL")
        .into_raw()
}

/// The first call returns a JSON array, the second NULL, every later call the
/// `FAULTY_OK` field.
///
/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_run(
    instance: *const c_void,
    _consequence_json: *const c_char,
    _variant_json: *const c_char,
) -> *mut c_char {
    let plugin = &*instance.cast::<FaultyPlugin>();
    let text = match plugin.runs.fetch_add(1, Ordering::SeqCst) {
        0 => "[1]",
        1 => return ptr::null_mut(),
        _ => r#"{"FAULTY_OK":"1"}"#,
    };
    CString::new(text).expect("no interior NUL").into_raw()
}

/// # Safety
/// `text` must be a pointer returned by [`vep_plugin_run`],
/// [`vep_plugin_header_info`] or [`vep_plugin_feature_types`] that has not been
/// freed.
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_free(text: *mut c_char) {
    if !text.is_null() {
        drop(CString::from_raw(text));
    }
}
