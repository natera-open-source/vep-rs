// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! C-ABI dylib plugins: loading the `--plugin` arguments no built-in claimed, and
//! running the loaded plugins over a batch so their JSON results land in
//! `TranscriptConsequence::plugin_data`.
//!
//! Perl citations name modules of ensembl-vep release/115 (`Bio/EnsEMBL/VEP/...`).

use anyhow::{bail, Context};
use serde_json::Value;

use vep_core::consequence::TranscriptConsequence;
use vep_core::variant::InputVariant;
use vep_plugin::{PluginError, PluginRegistry};

/// Loads every `--plugin` argument that matched no built-in as a dylib.
///
/// A name that resolves to no file, neither as a path nor as `lib<name>.so` /
/// `lib<name>.dylib` under `dir_plugins`, is warned about and dropped: Ensembl VEP
/// continues without a plugin module that fails to load unless `--safe` is given
/// (`Runner.pm` `get_all_Plugins`). A file that exists but does not open, reports
/// another ABI version, or fails `init` is an error, because the user named a
/// real plugin and it will not run.
pub fn load_unresolved(
    unresolved: &[String],
    dir_plugins: Option<&str>,
) -> anyhow::Result<PluginRegistry> {
    let mut registry = PluginRegistry::new();
    for arg in unresolved {
        let name = arg.split(',').next().unwrap_or(arg);
        match registry.load_from_args(std::slice::from_ref(arg), dir_plugins) {
            Ok(()) => {}
            // The registry reports an unresolvable name with exactly this message;
            // any other `Load` text comes from a file that exists.
            Err(PluginError::Load(msg)) if msg == format!("plugin '{name}' not found") => {
                tracing::warn!(
                    "plugin '{name}' not found as a built-in or as a dylib under --dir_plugins; skipped"
                );
            }
            Err(e) => {
                return Err(
                    anyhow::Error::new(e).context(format!("plugin '{name}' failed to load"))
                );
            }
        }
    }
    Ok(registry)
}

/// Runs every loaded dylib plugin on every transcript consequence in `batch`.
///
/// The variant is serialised once and that text is handed to every plugin for
/// every consequence of the variant; only the consequence is serialised per call.
pub fn run(registry: &PluginRegistry, batch: &mut [InputVariant]) -> anyhow::Result<()> {
    if registry.is_empty() {
        return Ok(());
    }
    for variant in batch.iter_mut() {
        let variant_json = serde_json::to_string(&*variant).with_context(|| {
            format!("serialising variant {} for plugins", variant.allele_string)
        })?;
        for tc in variant.transcript_consequences.iter_mut() {
            let consequence_json = serde_json::to_string(&*tc).with_context(|| {
                format!(
                    "serialising consequence on {} for plugins",
                    tc.transcript_id
                )
            })?;
            for plugin in registry.plugins() {
                let result = plugin
                    .run(&consequence_json, &variant_json)
                    .with_context(|| format!("plugin '{}' run failed", plugin.name()))?;
                merge_plugin_result(tc, plugin.name(), &result)?;
            }
        }
    }
    Ok(())
}

/// Merges one plugin's JSON result into `tc.plugin_data`.
///
/// The result must be a JSON object; each entry becomes one Extra-column field. A
/// string value is stored verbatim, every other value as its JSON text, so a
/// number stays `0.42` and a nested structure stays parseable.
pub fn merge_plugin_result(
    tc: &mut TranscriptConsequence,
    plugin_name: &str,
    result_json: &str,
) -> anyhow::Result<()> {
    let value: Value = serde_json::from_str(result_json)
        .with_context(|| format!("plugin '{plugin_name}' returned invalid JSON"))?;
    let Value::Object(entries) = value else {
        bail!(
            "plugin '{plugin_name}' returned a JSON {} where an object was expected",
            json_type_name(&value)
        );
    };
    for (key, value) in entries {
        let text = match value {
            Value::String(s) => s,
            other => other.to_string(),
        };
        tc.plugin_data.insert(key, text);
    }
    Ok(())
}

fn json_type_name(value: &Value) -> &'static str {
    match value {
        Value::Null => "null",
        Value::Bool(_) => "boolean",
        Value::Number(_) => "number",
        Value::String(_) => "string",
        Value::Array(_) => "array",
        Value::Object(_) => "object",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_plugin_result_inserts_strings_verbatim_and_stringifies_the_rest() {
        let mut tc = TranscriptConsequence::default();
        merge_plugin_result(&mut tc, "echo", r#"{"LABEL": "lof", "SCORE": 0.42}"#).unwrap();
        assert_eq!(tc.plugin_data.len(), 2);
        assert_eq!(tc.plugin_data["LABEL"], "lof");
        assert_eq!(tc.plugin_data["SCORE"], "0.42");
    }

    #[test]
    fn merge_plugin_result_keeps_existing_fields_and_overwrites_a_repeated_key() {
        let mut tc = TranscriptConsequence::default();
        tc.plugin_data
            .insert("REVEL_score".to_string(), "0.9".to_string());
        merge_plugin_result(&mut tc, "echo", r#"{"A": "1"}"#).unwrap();
        merge_plugin_result(&mut tc, "echo", r#"{"A": "2", "B": null}"#).unwrap();
        assert_eq!(tc.plugin_data["REVEL_score"], "0.9");
        assert_eq!(tc.plugin_data["A"], "2");
        assert_eq!(tc.plugin_data["B"], "null");
    }

    #[test]
    fn merge_plugin_result_rejects_a_non_object() {
        let mut tc = TranscriptConsequence::default();
        let err = merge_plugin_result(&mut tc, "echo", r#"["not", "an", "object"]"#).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("echo"), "{msg}");
        assert!(msg.contains("array"), "{msg}");
        assert!(tc.plugin_data.is_empty());

        let err = merge_plugin_result(&mut tc, "echo", "not json").unwrap_err();
        assert!(err.to_string().contains("echo"), "{err}");
    }

    #[test]
    fn load_unresolved_skips_an_unknown_name_and_returns_an_empty_registry() {
        let registry = load_unresolved(&["NoSuchPlugin,param".to_string()], None).unwrap();
        assert!(registry.is_empty());
        let registry =
            load_unresolved(&["NoSuchPlugin".to_string()], Some("/no/such/dir")).unwrap();
        assert!(registry.is_empty());
    }

    #[test]
    fn load_unresolved_fails_on_a_file_that_is_not_a_plugin() {
        let not_a_dylib = concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml");
        let err = load_unresolved(&[not_a_dylib.to_string()], None).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains(not_a_dylib), "{msg}");
        assert!(msg.contains("failed to load"), "{msg}");
    }

    #[test]
    fn run_over_an_empty_registry_touches_nothing() {
        let mut variant =
            InputVariant::new("1".to_string(), 100, 100, b"A".to_vec(), b"G".to_vec());
        variant
            .transcript_consequences
            .push(TranscriptConsequence::default());
        run(&PluginRegistry::new(), std::slice::from_mut(&mut variant)).unwrap();
        assert!(variant.transcript_consequences[0].plugin_data.is_empty());
    }
}
