// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! C-ABI dylib plugins: loading the `--plugin` arguments no built-in claimed, and
//! running the loaded plugins over a batch so their JSON results land in
//! `TranscriptConsequence::plugin_data`.
//!
//! Perl citations name modules of ensembl-vep release/116.2 (`Bio/EnsEMBL/VEP/...`).

use anyhow::{bail, Context};
use serde_json::Value;

use vep_core::consequence::TranscriptConsequence;
use vep_core::sv_routing::split_plugin_arg;
use vep_core::variant::InputVariant;
use vep_plugin::{PluginError, PluginRegistry};

/// Loads one `--plugin` argument that matched no built-in as a dylib into
/// `registry`.
///
/// Every failure drops the plugin with a warning and the run goes on without
/// it, or is the plugin component's error under `safe`, as Ensembl VEP treats
/// a plugin module it cannot use (`throw($msg) if $self->param('safe')`, else
/// `warning_msg` and `next`). A name that resolves to no file, neither as a
/// path nor as `lib<name>.so` / `lib<name>.dylib` under `dir_plugins`, is
/// `Failed to load plugin <name>`; a file that does not open as a shared
/// library, lacks a symbol or reports another ABI version is `Failed to compile
/// plugin <name>: <reason>`, the module that fails `use` (Runner.pm:929-939);
/// a file whose `init` returns non-zero, or whose constructor returns no
/// instance, is `Failed to instantiate plugin <name>: <reason>`, the
/// constructor that dies (`:946-954`).
pub fn load_unresolved(
    registry: &mut PluginRegistry,
    arg: &str,
    dir_plugins: &str,
    safe: bool,
) -> anyhow::Result<()> {
    let name = split_plugin_arg(arg).0;
    let message = match registry.load_from_args(&[arg.to_string()], Some(dir_plugins)) {
        Ok(()) => return Ok(()),
        // The registry reports an unresolvable name with exactly this message;
        // any other `Load` text comes from a file that exists.
        Err(PluginError::Load(msg)) if msg == format!("plugin '{name}' not found") => format!(
            "Failed to load plugin {name}: {name} is not a built-in plugin and no lib{name}.so or lib{name}.dylib is under {dir_plugins}"
        ),
        Err(PluginError::Init(reason)) => format!("Failed to instantiate plugin {name}: {reason}"),
        Err(PluginError::Load(reason)) => format!("Failed to compile plugin {name}: {reason}"),
        Err(other) => format!("Failed to compile plugin {name}: {other}"),
    };
    crate::annotator::plugin_warning(message, safe)
}

/// Runs every loaded dylib plugin on the consequences in `batch` whose feature
/// type it declares (`OutputFactory.pm:2223` `check_feature_type`; a plugin
/// declares `Transcript` unless its `vep_plugin_feature_types` says otherwise).
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
            let plugins: Vec<_> = registry
                .plugins()
                .iter()
                .filter(|plugin| plugin.annotates(tc.feature_type))
                .collect();
            if plugins.is_empty() {
                continue;
            }
            let consequence_json = serde_json::to_string(&*tc).with_context(|| {
                format!(
                    "serialising consequence on {} for plugins",
                    tc.transcript_id
                )
            })?;
            for plugin in plugins {
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

    /// An unknown name is dropped and the registry stays empty; under `safe` it is
    /// the plugin component's error naming the plugin and the directory searched
    /// (Runner.pm:929-939: `Failed to compile plugin`, thrown under `--safe`).
    #[test]
    fn load_unresolved_skips_an_unknown_name_and_dies_under_safe() {
        let mut registry = PluginRegistry::new();
        load_unresolved(&mut registry, "NoSuchPlugin,param", "/no/such/dir", false).unwrap();
        assert!(registry.is_empty());
        let err = load_unresolved(&mut registry, "NoSuchPlugin", "/no/such/dir", true).unwrap_err();
        assert_eq!(
            err.downcast_ref::<crate::config::CliError>().unwrap().to_string(),
            "ERROR: [plugin] Failed to load plugin NoSuchPlugin: NoSuchPlugin is not a built-in plugin and no libNoSuchPlugin.so or libNoSuchPlugin.dylib is under /no/such/dir"
        );
    }

    /// A file that exists but is no shared library is dropped as a module that
    /// fails to compile is (Runner.pm:929-939: `Failed to compile plugin
    /// <Name>: <reason>`, a warning and `next`, thrown under `--safe`), the
    /// reason the loader's own text naming the file.
    #[test]
    fn load_unresolved_drops_a_file_that_is_not_a_plugin_and_dies_under_safe() {
        let not_a_dylib = concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml");
        let mut registry = PluginRegistry::new();
        load_unresolved(&mut registry, not_a_dylib, "/no/such/dir", false).unwrap();
        assert!(registry.is_empty());
        let err = load_unresolved(&mut registry, not_a_dylib, "/no/such/dir", true).unwrap_err();
        let message = err
            .downcast_ref::<crate::config::CliError>()
            .unwrap()
            .to_string();
        assert!(
            message.starts_with(&format!(
                "ERROR: [plugin] Failed to compile plugin {not_a_dylib}: {not_a_dylib}: "
            )),
            "{message}"
        );
    }

    /// The consequence document a plugin receives carries these keys and no
    /// other: every field of `TranscriptConsequence` except `id_rank` and
    /// `record_span_in_cds`, with `cdna_total`, `cds_total`, `protein_total`,
    /// `plugin_data`, `loftee_ctx`, `cds_length`, `source_cache`,
    /// `own_span_positions`, `used_ref`, `shifted_allele` and `motif` present
    /// only when populated.
    #[test]
    fn consequence_json_carries_every_documented_key_and_nothing_else() {
        const ALWAYS: [&str; 46] = [
            "transcript_id",
            "gene_id",
            "gene_symbol",
            "gene_symbol_source",
            "hgnc_id",
            "consequences",
            "impact",
            "biotype",
            "canonical",
            "cdna_position",
            "cds_position",
            "protein_position",
            "amino_acids",
            "codons",
            "protein_id",
            "distance",
            "strand",
            "exon",
            "intron",
            "hgvsc",
            "hgvsp",
            "hgvs_offset",
            "shift_length",
            "sift",
            "polyphen",
            "domains",
            "feature_type",
            "flags",
            "tsl",
            "mane_select",
            "mane_plus_clinical",
            "appris",
            "ccds",
            "gencode_primary",
            "swissprot",
            "trembl",
            "uniparc",
            "uniprot_isoform",
            "refseq",
            "transcript_version",
            "gene_version",
            "gene_phenotype",
            "mirna",
            "cdna_length",
            "feature_start",
            "feature_end",
        ];
        let sorted_keys = |tc: &TranscriptConsequence| -> Vec<String> {
            let value = serde_json::to_value(tc).unwrap();
            let mut names: Vec<String> = value.as_object().unwrap().keys().cloned().collect();
            names.sort();
            names
        };
        let mut expected: Vec<String> = ALWAYS.iter().map(|k| k.to_string()).collect();
        expected.sort();

        let mut tc = TranscriptConsequence {
            id_rank: 7,
            record_span_in_cds: true,
            ..TranscriptConsequence::default()
        };
        assert_eq!(sorted_keys(&tc), expected);

        tc.plugin_data
            .insert("REVEL_score".to_string(), "0.9".to_string());
        tc.loftee_ctx = Some(Box::new(
            vep_core::consequence::LofteeContext::from_transcript(
                &crate::runner::make_runner_test_transcript("ENST1", "ENSG1", "1", 100),
            ),
        ));
        tc.cds_length = Some(300);
        tc.source_cache = Some(vep_core::consequence::SourceCache::Ensembl);
        tc.cdna_total = Some(2500);
        tc.cds_total = Some(300);
        tc.protein_total = Some(99);
        tc.own_span_positions = Some(Box::new(vep_core::consequence::OwnSpanPositions {
            cdna: Some("439".to_string()),
            cds: Some("87".to_string()),
            protein: Some("29".to_string()),
        }));
        tc.used_ref = Some(std::sync::Arc::from("A"));
        tc.shifted_allele = Some(std::sync::Arc::from("GAA"));
        tc.motif = Some(Box::new(vep_core::consequence::MotifFields {
            name: std::sync::Arc::from("ENSPFM0290"),
            pos: Some(7),
            factors: std::sync::Arc::from([std::sync::Arc::from("HNF4A")]),
        }));
        expected.extend(
            [
                "plugin_data",
                "loftee_ctx",
                "cds_length",
                "source_cache",
                "cdna_total",
                "cds_total",
                "protein_total",
                "own_span_positions",
                "used_ref",
                "shifted_allele",
                "motif",
            ]
            .map(str::to_string),
        );
        expected.sort();
        assert_eq!(sorted_keys(&tc), expected);
    }

    /// The variant document a plugin receives names the input line under
    /// `raw_input` alone: the bytes of a line that is not UTF-8, kept for the
    /// VCF output, are no key of it.
    #[test]
    fn variant_json_carries_the_input_line_as_text_only() {
        let mut variant = InputVariant::new("21".into(), 100, 100, b"A".to_vec(), b"G".to_vec());
        let line: &[u8] = b"21\t100\tid_\xE9\tA\tG\t.\t.\t.";
        variant.raw_input = Some(String::from_utf8_lossy(line).into());
        variant.raw_input_bytes = Some(line.into());
        let value = serde_json::to_value(&variant).unwrap();
        let object = value.as_object().unwrap();
        assert_eq!(
            object["raw_input"].as_str(),
            Some("21\t100\tid_\u{FFFD}\tA\tG\t.\t.\t.")
        );
        assert!(!object.contains_key("raw_input_bytes"), "{value}");
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
