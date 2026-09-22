// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Round trip through a real dylib: the `echo_plugin` example, built beside this
//! test by `cargo test`, is loaded, initialised and run on a serialised
//! consequence and variant, and its fields come back through the registry.

use std::path::PathBuf;

use vep_plugin::{LoadedPlugin, PluginRegistry};

const DYLIB_EXT: &str = if cfg!(target_os = "macos") {
    "dylib"
} else {
    "so"
};

/// The profile's `examples/` directory: this test binary runs from
/// `<profile>/deps/`, and cargo writes example artifacts to `<profile>/examples/`.
fn examples_dir() -> PathBuf {
    let exe = std::env::current_exe().expect("current_exe");
    let deps = exe.parent().expect("test binary has a parent directory");
    let profile = deps
        .parent()
        .expect("deps directory has a parent directory");
    profile.join("examples")
}

fn echo_plugin_path() -> PathBuf {
    let path = examples_dir().join(format!("libecho_plugin.{DYLIB_EXT}"));
    assert!(
        path.is_file(),
        "echo_plugin example dylib not found at {}; `cargo test -p vep-plugin` builds it, \
         a `--test echo_roundtrip` filter alone does not",
        path.display()
    );
    path
}

/// The host serialises a `TranscriptConsequence` with serde; the plugin reads only the
/// `transcript_id` field, so the fixture is that serialisation's shape, written out.
fn consequence_json(transcript_id: &str) -> String {
    format!(
        "{{\"transcript_id\":\"{transcript_id}\",\"gene_id\":\"ENSG00000000001\",\"consequences\":[\"missense_variant\"],\"plugin_data\":{{}}}}"
    )
}

fn variant_json() -> String {
    "{\"chr\":\"21\",\"start\":25000000,\"end\":25000000,\"ref_allele\":\"A\",\"alt_allele\":\"G\"}"
        .to_string()
}

/// The plugin returns a flat JSON object of string values; this reads it without a JSON
/// library so the crate carries no test-only dependency.
fn parse_object(text: &str) -> std::collections::BTreeMap<String, String> {
    let body = text
        .trim()
        .strip_prefix('{')
        .and_then(|t| t.strip_suffix('}'))
        .unwrap_or_else(|| panic!("plugin returned a non-object: {text}"));
    let mut out = std::collections::BTreeMap::new();
    for pair in body.split(',').filter(|p| !p.trim().is_empty()) {
        let (k, v) = pair
            .split_once(':')
            .unwrap_or_else(|| panic!("bad pair {pair} in {text}"));
        out.insert(
            k.trim().trim_matches('"').to_string(),
            v.trim().trim_matches('"').to_string(),
        );
    }
    out
}

#[test]
fn loaded_plugin_reports_name_and_version_and_echoes_after_init() {
    let mut plugin = LoadedPlugin::load(&echo_plugin_path()).expect("load echo plugin");
    assert_eq!(plugin.name(), "echo");
    assert_eq!(plugin.version(), "0.1.0");

    plugin.init(&["hello".to_string()]).expect("init");
    let result = plugin
        .run(&consequence_json("ENST00000123456"), &variant_json())
        .expect("run");
    let fields = parse_object(&result);
    assert_eq!(fields["ECHO_FEATURE"], "ENST00000123456");
    assert_eq!(fields["ECHO_PARAM"], "hello");
    assert_eq!(fields.len(), 2);
}

#[test]
fn plugin_without_init_params_omits_the_param_field() {
    let plugin = LoadedPlugin::load(&echo_plugin_path()).expect("load echo plugin");
    let result = plugin
        .run(&consequence_json("ENST00000000001"), &variant_json())
        .expect("run");
    let fields = parse_object(&result);
    assert_eq!(fields["ECHO_FEATURE"], "ENST00000000001");
    assert!(!fields.contains_key("ECHO_PARAM"));
}

#[test]
fn registry_loads_by_file_path_with_params() {
    let arg = format!("{},first,second", echo_plugin_path().display());
    let mut registry = PluginRegistry::new();
    registry.load_from_args(&[arg], None).expect("load by path");
    assert_eq!(registry.plugin_count(), 1);
    assert_eq!(registry.plugins()[0].name(), "echo");

    let results = registry
        .run_all(&consequence_json("ENST00000999999"), &variant_json())
        .expect("run_all");
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].0, "echo");
    let fields = parse_object(&results[0].1);
    assert_eq!(fields["ECHO_FEATURE"], "ENST00000999999");
    assert_eq!(fields["ECHO_PARAM"], "first");
}

#[test]
fn registry_resolves_a_bare_name_under_the_plugin_dir() {
    let dir = echo_plugin_path().parent().unwrap().to_path_buf();
    let mut registry = PluginRegistry::new();
    registry
        .load_from_args(&["echo_plugin".to_string()], dir.to_str())
        .expect("load by name under --dir_plugins");
    assert_eq!(registry.plugin_count(), 1);
    let results = registry
        .run_all(&consequence_json("ENST00000000002"), &variant_json())
        .expect("run_all");
    let fields = parse_object(&results[0].1);
    assert_eq!(fields["ECHO_FEATURE"], "ENST00000000002");
    assert!(!fields.contains_key("ECHO_PARAM"));
}

#[test]
fn registry_reports_a_bare_name_missing_from_the_plugin_dir_as_not_found() {
    let dir = examples_dir();
    let mut registry = PluginRegistry::new();
    let err = registry
        .load_from_args(&["no_such_plugin".to_string()], dir.to_str())
        .unwrap_err();
    assert_eq!(
        err.to_string(),
        "plugin load failed: plugin 'no_such_plugin' not found"
    );
    assert!(registry.is_empty());
}
