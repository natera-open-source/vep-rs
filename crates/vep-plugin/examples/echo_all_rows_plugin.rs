// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Example C-ABI plugin that runs on every row: it writes the kind of row it was
//! given as `ECHO_ROW`, the consequence's `feature_type` (`Transcript`,
//! `RegulatoryFeature`, `MotifFeature`) or `intergenic` where the consequence
//! document is `null`.
//!
//! Built as a `cdylib`. The library's `vep_plugin_feature_types` lists all four
//! kinds, which is what makes the host hand it every row, the intergenic ones
//! included; compiled with the `no-feature-types` feature it exports no such
//! symbol, the fixture for the host's `Transcript` default. Every string handed
//! to the host is a `CString::into_raw`, which `vep_plugin_free` reclaims with
//! `CString::from_raw`.

use std::ffi::{c_char, c_void, CStr, CString};
use std::ptr;

struct EchoAllRowsPlugin {
    name: CString,
    version: CString,
}

#[no_mangle]
pub extern "C" fn vep_plugin_abi_version() -> u32 {
    2
}

#[no_mangle]
pub extern "C" fn vep_plugin_create() -> *mut c_void {
    let plugin = EchoAllRowsPlugin {
        name: CString::new("echo_all_rows").expect("no interior NUL"),
        version: CString::new("0.1.0").expect("no interior NUL"),
    };
    Box::into_raw(Box::new(plugin)).cast()
}

/// # Safety
/// `instance` must be a pointer returned by [`vep_plugin_create`] that has not
/// been destroyed.
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_destroy(instance: *mut c_void) {
    if !instance.is_null() {
        drop(Box::from_raw(instance.cast::<EchoAllRowsPlugin>()));
    }
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_name(instance: *const c_void) -> *const c_char {
    (*instance.cast::<EchoAllRowsPlugin>()).name.as_ptr()
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_version(instance: *const c_void) -> *const c_char {
    (*instance.cast::<EchoAllRowsPlugin>()).version.as_ptr()
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
    CString::new(r#"[["ECHO_ROW","The kind of row: the feature type, or intergenic"]]"#)
        .expect("no interior NUL")
        .into_raw()
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[cfg(not(feature = "no-feature-types"))]
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_feature_types(_instance: *const c_void) -> *mut c_char {
    CString::new(r#"["Transcript","RegulatoryFeature","MotifFeature","Intergenic"]"#)
        .expect("no interior NUL")
        .into_raw()
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`]; the
/// consequence JSON must be a valid NUL-terminated string.
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_run(
    _instance: *const c_void,
    consequence_json: *const c_char,
    _variant_json: *const c_char,
) -> *mut c_char {
    let Ok(consequence_text) = CStr::from_ptr(consequence_json).to_str() else {
        return ptr::null_mut();
    };
    let kind = if consequence_text == "null" {
        "intergenic"
    } else {
        feature_type(consequence_text).unwrap_or("unknown")
    };
    match CString::new(format!("{{\"ECHO_ROW\":\"{kind}\"}}")) {
        Ok(text) => text.into_raw(),
        Err(_) => ptr::null_mut(),
    }
}

/// The consequence's `feature_type`, one of three fixed names, so a scan to the
/// closing quote reads it.
fn feature_type(json: &str) -> Option<&str> {
    let (_, after_key) = json.split_once("\"feature_type\":\"")?;
    let (value, _) = after_key.split_once('"')?;
    Some(value)
}

/// # Safety
/// `text` must be a pointer returned by [`vep_plugin_run`],
/// [`vep_plugin_header_info`] or `vep_plugin_feature_types` that has not been
/// freed.
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_free(text: *mut c_char) {
    if !text.is_null() {
        drop(CString::from_raw(text));
    }
}
