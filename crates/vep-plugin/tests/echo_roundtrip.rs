// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Round trip through real dylibs: the `echo_plugin` and `echo_all_rows_plugin`
//! examples, built beside this test by `cargo test`, are loaded, initialised and
//! run on a serialised consequence and variant, and their fields, header info and
//! feature types come back through the loader and the registry.

use std::ffi::CString;
use std::path::PathBuf;
use std::process::Command;

use vep_plugin::{FeatureType, LoadedPlugin, PluginError, PluginRegistry};

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

fn example_path(name: &str) -> PathBuf {
    let path = examples_dir().join(format!("lib{name}.{DYLIB_EXT}"));
    assert!(
        path.is_file(),
        "{name} example dylib not found at {}; `cargo test -p vep-plugin` builds it, \
         a `--test echo_roundtrip` filter alone does not",
        path.display()
    );
    path
}

fn echo_plugin_path() -> PathBuf {
    example_path("echo_plugin")
}

fn echo_all_rows_plugin_path() -> PathBuf {
    example_path("echo_all_rows_plugin")
}

/// The host serialises a `TranscriptConsequence` with serde; the plugins read only
/// the `transcript_id` and `feature_type` fields, so the fixture is that
/// serialisation's shape, written out.
fn consequence_json(transcript_id: &str) -> CString {
    CString::new(format!(
        "{{\"transcript_id\":\"{transcript_id}\",\"gene_id\":\"ENSG00000000001\",\"consequences\":[\"missense_variant\"],\"feature_type\":\"Transcript\",\"plugin_data\":{{}}}}"
    ))
    .unwrap()
}

fn variant_json() -> CString {
    CString::new(
        "{\"chr\":\"21\",\"start\":25000000,\"end\":25000000,\"ref_allele\":\"A\",\"alt_allele\":\"G\"}",
    )
    .unwrap()
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

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|s| s.to_string()).collect()
}

#[test]
fn loaded_plugin_reports_name_and_version_and_echoes_after_init() {
    let plugin =
        LoadedPlugin::load(&echo_plugin_path(), &strings(&["hello"])).expect("load echo plugin");
    assert_eq!(plugin.name(), "echo");
    assert_eq!(plugin.version(), "0.2.0");

    let result = plugin
        .run(&consequence_json("ENST00000123456"), &variant_json())
        .expect("run");
    let fields = parse_object(&result);
    assert_eq!(fields["ECHO_FEATURE"], "ENST00000123456");
    assert_eq!(fields["ECHO_PARAM"], "hello");
    assert_eq!(fields.len(), 2);
}

/// The plugin lists its keys with `ECHO_PARAM` first; the loader sorts them.
#[test]
fn loaded_plugin_reports_its_header_info() {
    let plugin = LoadedPlugin::load(&echo_plugin_path(), &[]).expect("load echo plugin");
    assert_eq!(
        plugin.header_info(),
        &[
            (
                "ECHO_FEATURE".to_string(),
                "The transcript id the consequence names".to_string()
            ),
            (
                "ECHO_PARAM".to_string(),
                "The first plugin parameter".to_string()
            ),
        ]
    );
}

/// A load without parameters still calls `init`, once, with zero parameters: the
/// example returns nothing from `run` before its one `init` and refuses a second,
/// so a successful run without `ECHO_PARAM` is `init_calls == 1` and `count == 0`.
#[test]
fn init_is_called_with_zero_params() {
    let plugin = LoadedPlugin::load(&echo_plugin_path(), &[]).expect("load without parameters");
    let result = plugin
        .run(&consequence_json("ENST00000000001"), &variant_json())
        .expect("run after a parameterless init");
    let fields = parse_object(&result);
    assert_eq!(fields["ECHO_FEATURE"], "ENST00000000001");
    assert!(!fields.contains_key("ECHO_PARAM"));
    assert_eq!(fields.len(), 1);
}

/// A library's `vep_plugin_feature_types` names the kinds of row it runs on; a
/// library that exports no such symbol runs on transcript rows alone, Ensembl
/// VEP's default `feature_types` (`BaseVepPlugin.pm:99`). The
/// `fixture_no_feature_types` fixture is `echo_all_rows_plugin.rs` compiled
/// without the export.
#[test]
fn feature_types_export_is_optional_and_read() {
    let echo = LoadedPlugin::load(&echo_plugin_path(), &[]).expect("load echo plugin");
    assert_eq!(echo.feature_types(), [FeatureType::Transcript]);
    assert!(echo.runs_on(FeatureType::Transcript));
    assert!(!echo.runs_on(FeatureType::RegulatoryFeature));
    assert!(!echo.runs_on(FeatureType::Intergenic));

    let all_rows =
        LoadedPlugin::load(&echo_all_rows_plugin_path(), &[]).expect("load echo_all_rows plugin");
    assert_eq!(all_rows.name(), "echo_all_rows");
    assert_eq!(all_rows.feature_types(), FeatureType::ALL);
    for &row in FeatureType::ALL {
        assert!(all_rows.runs_on(row), "{row:?}");
    }
    let on_row = all_rows
        .run(&consequence_json("ENST00000000001"), &variant_json())
        .expect("run on a transcript row");
    assert_eq!(parse_object(&on_row)["ECHO_ROW"], "Transcript");
    let intergenic = all_rows
        .run(c"null", &variant_json())
        .expect("run on an intergenic row");
    assert_eq!(parse_object(&intergenic)["ECHO_ROW"], "intergenic");

    let fixture = compile_fixture(
        "echo_all_rows_plugin.rs",
        "echo_all_rows_no_feature_types",
        "fixture_no_feature_types",
    );
    let unlisted = LoadedPlugin::load(&fixture, &[]).expect("load the fixture without the export");
    assert_eq!(unlisted.feature_types(), [FeatureType::Transcript]);
    assert!(unlisted.runs_on(FeatureType::Transcript));
    for row in [
        FeatureType::RegulatoryFeature,
        FeatureType::MotifFeature,
        FeatureType::Intergenic,
    ] {
        assert!(!unlisted.runs_on(row), "{row:?}");
    }
    std::fs::remove_file(&fixture).unwrap();
}

/// An example compiled by `rustc` alone with the `cfg` `fixture_cfg` set, into
/// `CARGO_TARGET_TMPDIR`: cargo builds each example once, without it, and the
/// other tests load that build. The build's `ECHO_DESTROY_MARKER` names the
/// fixture's path with `.destroyed` appended.
fn compile_fixture(example: &str, crate_name: &str, fixture_cfg: &str) -> PathBuf {
    let source = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join(example);
    let out_dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    std::fs::create_dir_all(&out_dir).unwrap();
    let fixture = out_dir.join(format!(
        "lib{crate_name}_{}.{DYLIB_EXT}",
        std::process::id()
    ));
    let rustc = std::env::var_os("RUSTC").unwrap_or_else(|| "rustc".into());
    let output = Command::new(rustc)
        .args(["--edition", "2021", "--crate-type", "cdylib"])
        .args(["--crate-name", crate_name])
        .args(["--cfg", fixture_cfg])
        .env("ECHO_DESTROY_MARKER", destroy_marker(&fixture))
        .arg("-o")
        .arg(&fixture)
        .arg(&source)
        .output()
        .expect("rustc runs");
    assert!(
        output.status.success(),
        "compiling the {fixture_cfg} fixture failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    fixture
}

/// The file a fixture's `vep_plugin_destroy` adds a line to.
fn destroy_marker(fixture: &std::path::Path) -> PathBuf {
    let mut marker = fixture.as_os_str().to_owned();
    marker.push(".destroyed");
    PathBuf::from(marker)
}

/// Each fixture is `echo_plugin.rs` compiled with `--cfg fixture_fail="<stage>"`,
/// failing that stage. Once `create` has returned an instance, every failure
/// destroys it, once, while the library is still loaded (`loader.rs`, the
/// instance owned from `Self` on); the marker file outlives the library.
#[test]
fn every_load_failure_after_create_destroys_the_instance_once() {
    let stages: [(&str, fn(&PluginError) -> bool, &str, usize); 6] = [
        (
            "create",
            |e| matches!(e, PluginError::Init(_)),
            "plugin create returned null",
            0,
        ),
        (
            "name",
            |e| matches!(e, PluginError::Load(_)),
            "invalid UTF-8 in plugin name",
            1,
        ),
        (
            "init",
            |e| matches!(e, PluginError::Init(_)),
            "plugin 'echo' init returned error code 1",
            1,
        ),
        (
            "header_info",
            |e| matches!(e, PluginError::Load(_)),
            "plugin header info is not a JSON array of [key, description] pairs",
            1,
        ),
        (
            "feature_types",
            |e| matches!(e, PluginError::Load(_)),
            "plugin 'echo' lists an unknown feature type 'Gene'",
            1,
        ),
        (
            "run",
            |e| matches!(e, PluginError::Run(_)),
            "invalid UTF-8 in plugin result",
            1,
        ),
    ];
    for (stage, variant, message, destroys) in stages {
        let fixture = compile_fixture(
            "echo_plugin.rs",
            &format!("echo_plugin_fail_{stage}"),
            &format!("fixture_fail=\"{stage}\""),
        );
        let marker = destroy_marker(&fixture);
        let err = match LoadedPlugin::load(&fixture, &[]) {
            Ok(plugin) => plugin
                .run(&consequence_json("ENST00000000001"), &variant_json())
                .unwrap_err(),
            Err(err) => err,
        };
        assert!(variant(&err), "{stage}: {err:?}");
        assert!(err.to_string().contains(message), "{stage}: {err}");
        let destroyed = std::fs::read_to_string(&marker).unwrap_or_default();
        assert_eq!(destroyed.lines().count(), destroys, "{stage}");
        let _ = std::fs::remove_file(&marker);
        std::fs::remove_file(&fixture).unwrap();
    }
}

/// A plugin that lists `RegulatoryFeature` or `MotifFeature` switches regulatory
/// annotation on for the run, as Ensembl VEP's `Runner.pm:998-1002` does for a
/// plugin whose declared types name motif or regulatory features; the default
/// `Transcript` list switches nothing on. The example takes its list from a
/// `feature_types=` parameter, read after `init`; the all-rows example lists
/// every kind, so it switches the annotation on and its field prints on every
/// kind of row.
#[test]
fn a_plugin_listing_regulatory_kinds_switches_regulatory_annotation_on() {
    let path = echo_plugin_path().display().to_string();
    let registry_for = |arg: String| {
        let mut registry = PluginRegistry::new();
        registry.load_from_args(&[arg], &[]).expect("load");
        registry
    };
    assert!(registry_for(path.clone())
        .plugins_wanting_regulatory_features()
        .is_empty());
    let motif = registry_for(format!(
        "{path},feature_types=RegulatoryFeature;MotifFeature"
    ));
    assert_eq!(
        motif.plugins()[0].feature_types(),
        [FeatureType::RegulatoryFeature, FeatureType::MotifFeature]
    );
    assert_eq!(motif.plugins_wanting_regulatory_features(), ["echo"]);
    assert!(
        registry_for(format!("{path},feature_types=Transcript;Intergenic"))
            .plugins_wanting_regulatory_features()
            .is_empty()
    );
    let mut all_rows = PluginRegistry::new();
    all_rows
        .load_from_args(&[echo_all_rows_plugin_path().display().to_string()], &[])
        .expect("load echo_all_rows plugin");
    assert_eq!(
        all_rows.plugins_wanting_regulatory_features(),
        ["echo_all_rows"]
    );
    assert_eq!(
        all_rows.field_row_kinds(),
        [("ECHO_ROW".to_string(), FeatureType::ALL.to_vec())]
    );
}

/// The v1 fixture is `echo_plugin.rs` compiled with `--cfg fixture_abi_one`,
/// under which `vep_plugin_abi_version` returns 1.
#[test]
fn a_v1_library_is_refused() {
    let fixture = compile_fixture("echo_plugin.rs", "echo_plugin_abi_one", "fixture_abi_one");

    let err = LoadedPlugin::load(&fixture, &[]).unwrap_err();
    assert!(
        matches!(
            err,
            PluginError::AbiMismatch {
                expected: 2,
                got: 1
            }
        ),
        "{err:?}"
    );
    assert_eq!(
        err.to_string(),
        "incompatible plugin ABI version: expected 2, got 1"
    );
    std::fs::remove_file(&fixture).unwrap();
}

/// The fixture is `echo_plugin.rs` compiled with `--cfg fixture_missing_free`,
/// which leaves `vep_plugin_free` out. The load error names the symbol, then
/// the dynamic loader's own words ("undefined symbol" from dlsym on Linux,
/// "symbol not found" on macOS), never the library's bare "dlsym failed".
#[test]
fn a_library_missing_a_required_symbol_names_the_symbol_and_the_loaders_reason() {
    let fixture = compile_fixture(
        "echo_plugin.rs",
        "echo_plugin_missing_free",
        "fixture_missing_free",
    );

    let err = LoadedPlugin::load(&fixture, &[]).unwrap_err();
    let PluginError::Load(msg) = &err else {
        panic!("expected PluginError::Load, got: {err:?}");
    };
    let reason = msg
        .strip_prefix("missing vep_plugin_free: ")
        .unwrap_or_else(|| panic!("{msg}"));
    assert!(reason.to_ascii_lowercase().contains("symbol"), "{msg}");
    std::fs::remove_file(&fixture).unwrap();
}

#[test]
fn registry_loads_by_file_path_with_params() {
    let arg = format!("{},first,second", echo_plugin_path().display());
    let mut registry = PluginRegistry::new();
    registry.load_from_args(&[arg], &[]).expect("load by path");
    assert_eq!(registry.plugin_count(), 1);
    assert_eq!(registry.plugins()[0].name(), "echo");

    let results = registry
        .run_all(
            consequence_json("ENST00000999999").to_str().unwrap(),
            variant_json().to_str().unwrap(),
        )
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
        .load_from_args(
            &strings(&["echo_plugin"]),
            &strings(&[dir.to_str().unwrap()]),
        )
        .expect("load by name under --dir_plugins");
    assert_eq!(registry.plugin_count(), 1);
    let results = registry
        .run_all(
            consequence_json("ENST00000000002").to_str().unwrap(),
            variant_json().to_str().unwrap(),
        )
        .expect("run_all");
    let fields = parse_object(&results[0].1);
    assert_eq!(fields["ECHO_FEATURE"], "ENST00000000002");
    assert!(!fields.contains_key("ECHO_PARAM"));
}

/// Plugins in load order, each plugin's keys sorted (`OutputFactory.pm`
/// `get_plugin_headers`).
#[test]
fn registry_lists_header_info_plugin_by_plugin_in_load_order() {
    let dir = examples_dir();
    let mut registry = PluginRegistry::new();
    registry
        .load_from_args(
            &strings(&["echo_all_rows_plugin", "echo_plugin,hello"]),
            &strings(&[dir.to_str().unwrap()]),
        )
        .expect("load both examples");
    let names: Vec<String> = registry
        .all_header_info()
        .into_iter()
        .map(|(name, _)| name)
        .collect();
    assert_eq!(names, strings(&["ECHO_ROW", "ECHO_FEATURE", "ECHO_PARAM"]));
    assert!(registry.runs_on(FeatureType::Intergenic));
    assert!(registry.runs_on(FeatureType::Transcript));
}

/// Under the plugins directories a bare name takes a regular `lib<name>` file
/// alone: a directory of that name, and a plain entry of the bare name itself,
/// are passed over for the library a later directory holds, or for the
/// not-found error, as Ensembl VEP 116.3's `use $module` takes only a module
/// file from `@INC` (Runner.pm:931-933).
#[test]
fn registry_resolves_a_bare_name_only_to_a_library_file_under_the_plugin_dirs() {
    let scratch = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("bare-name-resolution-{}", std::process::id()));
    let decoys = scratch.join("decoys");
    std::fs::create_dir_all(decoys.join(format!("libGeneSplicer.{DYLIB_EXT}"))).unwrap();
    std::fs::write(decoys.join("GeneSplicer"), "").unwrap();
    let decoys_dir = decoys.to_str().unwrap().to_string();

    let mut registry = PluginRegistry::new();
    let err = registry
        .load_from_args(
            &strings(&["GeneSplicer"]),
            std::slice::from_ref(&decoys_dir),
        )
        .unwrap_err();
    assert_eq!(
        err.to_string(),
        "plugin load failed: plugin 'GeneSplicer' not found"
    );
    assert!(registry.is_empty());

    let libraries = scratch.join("libraries");
    std::fs::create_dir_all(&libraries).unwrap();
    std::fs::copy(
        echo_plugin_path(),
        libraries.join(format!("libGeneSplicer.{DYLIB_EXT}")),
    )
    .unwrap();
    registry
        .load_from_args(
            &strings(&["GeneSplicer"]),
            &[decoys_dir, libraries.to_str().unwrap().to_string()],
        )
        .expect("the library under the later directory loads");
    assert_eq!(registry.plugin_count(), 1);
    assert_eq!(registry.plugins()[0].name(), "echo");
    std::fs::remove_dir_all(&scratch).unwrap();
}

#[test]
fn registry_reports_a_bare_name_missing_from_the_plugin_dir_as_not_found() {
    let dir = examples_dir();
    let mut registry = PluginRegistry::new();
    let err = registry
        .load_from_args(
            &strings(&["no_such_plugin"]),
            &strings(&[dir.to_str().unwrap()]),
        )
        .unwrap_err();
    assert_eq!(
        err.to_string(),
        "plugin load failed: plugin 'no_such_plugin' not found"
    );
    assert!(registry.is_empty());
}
