// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Example C-ABI plugin: echoes the consequence's transcript id, and the first
//! `--plugin` parameter when one was given, back as Extra-column fields.
//!
//! Built as a `cdylib`, so the ten `vep_plugin_*` symbols below are exactly what
//! `LoadedPlugin::load` looks up, with the signatures in `vep_plugin::abi::ffi`.
//! The plugin lists `Transcript` as its one feature type, so the host hands it
//! transcript rows and no other, unless a `feature_types=A;B` parameter names
//! other kinds (semicolons between them, since the `--plugin` list splits on
//! commas); an `undeclared` parameter adds a field the header info does not
//! declare. Every string handed to the host is a `CString::into_raw`, which
//! `vep_plugin_free` reclaims with `CString::from_raw`; the name and version
//! pointers stay valid because their `CString`s live inside the boxed instance.
//! The consequence JSON is read with a string scan for the one field the plugin
//! echoes, so the example adds no dependency to the crate.
//!
//! The host calls `vep_plugin_init` exactly once per instance, with zero or more
//! parameters, before any `vep_plugin_run`; the plugin relies on that order and
//! treats a departure from it as a defect of the host.

use std::ffi::{c_char, c_void, CStr, CString};
use std::ptr;

struct EchoPlugin {
    name: CString,
    version: CString,
    init_calls: u32,
    param: Option<String>,
    feature_types: String,
    undeclared: bool,
}

#[no_mangle]
pub extern "C" fn vep_plugin_abi_version() -> u32 {
    if cfg!(fixture_abi_one) {
        1
    } else {
        2
    }
}

#[no_mangle]
pub extern "C" fn vep_plugin_create() -> *mut c_void {
    let plugin = EchoPlugin {
        name: CString::new("echo").expect("no interior NUL"),
        version: CString::new("0.2.0").expect("no interior NUL"),
        init_calls: 0,
        param: None,
        feature_types: r#"["Transcript"]"#.to_string(),
        undeclared: false,
    };
    Box::into_raw(Box::new(plugin)).cast()
}

/// # Safety
/// `instance` must be a pointer returned by [`vep_plugin_create`] that has not
/// been destroyed.
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_destroy(instance: *mut c_void) {
    if !instance.is_null() {
        drop(Box::from_raw(instance.cast::<EchoPlugin>()));
    }
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_name(instance: *const c_void) -> *const c_char {
    (*instance.cast::<EchoPlugin>()).name.as_ptr()
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_version(instance: *const c_void) -> *const c_char {
    (*instance.cast::<EchoPlugin>()).version.as_ptr()
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`], and
/// `params` must point to `count` valid NUL-terminated strings.
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_init(
    instance: *mut c_void,
    params: *const *const c_char,
    count: usize,
) -> i32 {
    let plugin = &mut *instance.cast::<EchoPlugin>();
    plugin.init_calls += 1;
    if plugin.init_calls > 1 {
        return 2;
    }
    for i in 0..count {
        let Ok(param) = CStr::from_ptr(*params.add(i)).to_str() else {
            return 1;
        };
        if let Some(kinds) = param.strip_prefix("feature_types=") {
            let listed: Vec<String> = kinds
                .split(';')
                .filter(|kind| !kind.is_empty())
                .map(|kind| format!("\"{kind}\""))
                .collect();
            plugin.feature_types = format!("[{}]", listed.join(","));
        } else if param == "undeclared" {
            plugin.undeclared = true;
        } else if plugin.param.is_none() {
            plugin.param = Some(param.to_owned());
        }
    }
    0
}

/// The two fields the plugin declares, listed in any order: the host sorts a
/// plugin's keys for the header and the columns. The field the `undeclared`
/// parameter adds is left out on purpose.
///
/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_header_info(_instance: *const c_void) -> *mut c_char {
    let pairs = r#"[["ECHO_PARAM","The first plugin parameter"],["ECHO_FEATURE","The transcript id the consequence names"]]"#;
    CString::new(pairs).expect("no interior NUL").into_raw()
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`].
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_feature_types(instance: *const c_void) -> *mut c_char {
    let plugin = &*instance.cast::<EchoPlugin>();
    match CString::new(plugin.feature_types.as_str()) {
        Ok(text) => text.into_raw(),
        Err(_) => ptr::null_mut(),
    }
}

/// # Safety
/// `instance` must be a live pointer returned by [`vep_plugin_create`], and both
/// JSON arguments must be valid NUL-terminated strings.
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_run(
    instance: *const c_void,
    consequence_json: *const c_char,
    _variant_json: *const c_char,
) -> *mut c_char {
    let plugin = &*instance.cast::<EchoPlugin>();
    if plugin.init_calls != 1 {
        return ptr::null_mut();
    }
    let Ok(consequence_text) = CStr::from_ptr(consequence_json).to_str() else {
        return ptr::null_mut();
    };
    let transcript_id = json_string_field(consequence_text, "transcript_id").unwrap_or_default();

    let mut result = format!("{{\"ECHO_FEATURE\":{}", json_quote(&transcript_id));
    if let Some(param) = &plugin.param {
        result.push_str(&format!(",\"ECHO_PARAM\":{}", json_quote(param)));
    }
    if plugin.undeclared {
        result.push_str(",\"ECHO_UNDECLARED\":\"u\"");
    }
    result.push('}');
    match CString::new(result) {
        Ok(text) => text.into_raw(),
        Err(_) => ptr::null_mut(),
    }
}

/// The string value of a top-level `"key":"value"` pair, with JSON escapes decoded.
fn json_string_field(json: &str, key: &str) -> Option<String> {
    let needle = format!("\"{key}\":");
    let start = json.find(&needle)? + needle.len();
    let rest = json[start..].trim_start();
    let body = rest.strip_prefix('"')?;
    let mut out = String::new();
    let mut chars = body.chars();
    while let Some(c) = chars.next() {
        match c {
            '"' => return Some(out),
            '\\' => match chars.next()? {
                'n' => out.push('\n'),
                't' => out.push('\t'),
                'u' => {
                    let code: String = chars.by_ref().take(4).collect();
                    out.push(char::from_u32(u32::from_str_radix(&code, 16).ok()?)?);
                }
                other => out.push(other),
            },
            other => out.push(other),
        }
    }
    None
}

/// `s` as a JSON string literal.
fn json_quote(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// # Safety
/// `text` must be a pointer returned by [`vep_plugin_run`],
/// [`vep_plugin_header_info`] or [`vep_plugin_feature_types`] that has not been
/// freed.
#[cfg(not(fixture_missing_free))]
#[no_mangle]
pub unsafe extern "C" fn vep_plugin_free(text: *mut c_char) {
    if !text.is_null() {
        drop(CString::from_raw(text));
    }
}
