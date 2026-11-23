// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! C-ABI dylib plugins: loading the `--plugin` arguments no built-in claimed, and
//! running the loaded plugins over a batch so their JSON results land in
//! `TranscriptConsequence::plugin_data`, or in `InputVariant::plugin_data` on an
//! intergenic row.
//!
//! Ensembl VEP source citations name modules of ensembl-vep release/116.3
//! (`Bio/EnsEMBL/VEP/...`).

use std::ffi::{CStr, CString};
use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context};
use indexmap::IndexMap;
use serde_json::Value;

use crate::config::{CliError, PluginDirs};
use vep_core::sv_routing::split_plugin_arg;
use vep_core::variant::InputVariant;
use vep_plugin::registry::names_a_file;
use vep_plugin::{FeatureType, LoadedPlugin, PluginError, PluginRegistry};

/// What [`load_unresolved`] did with a `--plugin` argument no built-in claimed.
#[derive(Debug, PartialEq, Eq)]
pub enum DylibOutcome {
    /// The argument named a file: the plugin loaded, or its failure was warned
    /// about and the run goes on without it (the plugin component's error under
    /// `safe`).
    Settled,
    /// No file answers to the name, as a path or as `lib<name>.so` /
    /// `lib<name>.dylib` under the plugins directories.
    NotFound,
}

/// Loads one `--plugin` argument that matched no built-in as a dylib into
/// `registry`.
///
/// A path to an Ensembl VEP plugin module ([`is_perl_module`]) is refused
/// before any open ([`perl_plugin_refusal`]), with or without `safe`. Every
/// other failure of a file drops the plugin with a warning and the run goes on
/// without it, or is the plugin component's error under `safe`, as Ensembl VEP
/// treats a plugin module it cannot use (`throw($msg) if $self->param('safe')`,
/// else `warning_msg` and `next`): a file that does not open as a shared
/// library, lacks a symbol or reports another ABI version is `Failed to compile
/// plugin <name>: <reason>`, the module that fails `use` (Runner.pm:929-939); a
/// file whose `init` returns non-zero, or whose constructor returns no
/// instance, is `Failed to instantiate plugin <name>: <reason>`, the
/// constructor that dies (`:946-954`). A name that resolves to no file, neither
/// as a path nor under one of the plugins directories, searched in their order,
/// is returned as [`DylibOutcome::NotFound`] for the caller to settle. A bare
/// name ([`names_a_file`] false) is read against the directories alone, never
/// the working directory, as `use $module` reads `@INC` alone (`:921`, `:931`).
pub fn load_unresolved(
    registry: &mut PluginRegistry,
    arg: &str,
    dir_plugins: &PluginDirs,
    safe: bool,
) -> anyhow::Result<DylibOutcome> {
    let name = split_plugin_arg(arg).0;
    let path = Path::new(name);
    if names_a_file(name) && is_perl_module(path) {
        bail!(perl_plugin_refusal(path));
    }
    let message = match registry.load_from_args(&[arg.to_string()], &dir_plugins.dirs) {
        Ok(()) => return Ok(DylibOutcome::Settled),
        // The registry reports an unresolvable name with exactly this message;
        // any other `Load` text comes from a file that exists.
        Err(PluginError::Load(msg)) if msg == format!("plugin '{name}' not found") => {
            return Ok(DylibOutcome::NotFound);
        }
        Err(PluginError::Init(reason)) => format!("Failed to instantiate plugin {name}: {reason}"),
        Err(PluginError::Load(reason)) => format!("Failed to compile plugin {name}: {reason}"),
        Err(other) => format!("Failed to compile plugin {name}: {other}"),
    };
    crate::annotator::plugin_warning(message, safe)?;
    Ok(DylibOutcome::Settled)
}

/// Whether `path` is a file holding an Ensembl VEP plugin module: its suffix is
/// `.pm`, or its first non-blank line begins `package `. Only the file's first
/// 4 KiB are read, so a shared library costs one short read.
pub fn is_perl_module(path: &Path) -> bool {
    if !path.is_file() {
        return false;
    }
    if path.extension().is_some_and(|suffix| suffix == "pm") {
        return true;
    }
    let Ok(file) = std::fs::File::open(path) else {
        return false;
    };
    let mut head = Vec::with_capacity(4096);
    if file.take(4096).read_to_end(&mut head).is_err() {
        return false;
    }
    String::from_utf8_lossy(&head)
        .lines()
        .map(str::trim_start)
        .find(|line| !line.is_empty())
        .is_some_and(|line| line.starts_with("package "))
}

/// `<dir>/<name>.pm` under the first plugins directory, in search order, that
/// holds it: the module an Ensembl VEP `--plugin <name>` would load.
pub fn perl_module_under(dir_plugins: &PluginDirs, name: &str) -> Option<PathBuf> {
    dir_plugins
        .dirs
        .iter()
        .map(|dir| Path::new(dir).join(format!("{name}.pm")))
        .find(|module| module.is_file())
}

/// The plugin component's error for an Ensembl VEP plugin module: vep-rs
/// embeds no Perl, so the module is refused by name, with or without `safe`.
pub fn perl_plugin_refusal(module: &Path) -> CliError {
    CliError::new(
        "plugin",
        format!(
            "{} is an Ensembl VEP Perl plugin; vep-rs runs no Perl plugin, see docs/plugins.md",
            module.display()
        ),
    )
}

/// Runs every loaded dylib plugin over `batch`: on each consequence whose
/// feature type the plugin listed (`Transcript` alone when it listed none), and
/// on the intergenic row of a variant that overlaps no transcript, its regulatory
/// and motif rows notwithstanding, for a plugin that listed `Intergenic`
/// (OutputFactory.pm:2229-2232, `check_feature_type(ref($bvfoa->feature) ||
/// 'Intergenic')`). A plugin's result lands in the consequence's `plugin_data`,
/// or the variant's on an intergenic row, where the Extra column and the plugin
/// columns read it on the rows of the kinds the plugin listed.
///
/// Each document is serialised at most once, and only when a loaded plugin runs
/// on one of its rows: the variant when some plugin takes one of its rows, a
/// consequence when some plugin takes it; the C string handed to every plugin is
/// made once per document. The intergenic row's consequence document is `null`.
///
/// A plugin's failure on a row is a warning and costs that row the plugin's
/// fields alone (`run_plugin_on_row`); the error is the host's own
/// serialisation of a row.
pub fn run(registry: &PluginRegistry, batch: &mut [InputVariant]) -> anyhow::Result<()> {
    if registry.is_empty() {
        return Ok(());
    }
    let intergenic_wanted = registry.runs_on(FeatureType::Intergenic);
    for variant in batch.iter_mut() {
        let intergenic_row = intergenic_wanted
            && (variant.intergenic_overlap || variant.transcript_consequences.is_empty());
        let some_row_wanted = intergenic_row
            || variant
                .transcript_consequences
                .iter()
                .any(|tc| registry.runs_on(tc.feature_type.into()));
        if !some_row_wanted {
            continue;
        }
        let variant_json = variant_c_json(variant)?;
        if intergenic_row {
            run_plugins_on_intergenic_row(registry, &variant_json, variant);
        }
        for tc in variant.transcript_consequences.iter_mut() {
            let row = FeatureType::from(tc.feature_type);
            if !registry.runs_on(row) {
                continue;
            }
            let consequence_json = c_json(&*tc).with_context(|| {
                format!(
                    "serialising consequence on {} for plugins",
                    tc.transcript_id()
                )
            })?;
            run_plugins_on_row(
                registry,
                row,
                &consequence_json,
                &variant_json,
                &mut tc.plugin_data,
            );
        }
    }
    Ok(())
}

/// Runs every plugin that takes a row of kind `row` on one consequence and
/// variant document pair, merging each result into `fields`, the row's
/// `plugin_data` ([`run_plugin_on_row`]).
fn run_plugins_on_row(
    registry: &PluginRegistry,
    row: FeatureType,
    consequence_json: &CStr,
    variant_json: &CStr,
    fields: &mut IndexMap<String, String>,
) {
    for plugin in registry
        .plugins()
        .iter()
        .filter(|plugin| plugin.runs_on(row))
    {
        run_plugin_on_row(plugin, consequence_json, variant_json, fields);
    }
}

/// The intergenic row's plugins, with a `null` consequence document. The
/// result lands in the variant's `plugin_data`, the intergenic row having no
/// consequence of its own, and its keys are noted on the variant
/// ([`InputVariant::intergenic_row_plugin_keys`]), so they print on that row
/// alone whatever the plugin declared: Ensembl VEP keeps a plugin's result in
/// the hash of the row it ran on (`OutputFactory.pm:2229-2244`), and the rows
/// of the other kinds carry their own results.
fn run_plugins_on_intergenic_row(
    registry: &PluginRegistry,
    variant_json: &CStr,
    variant: &mut InputVariant,
) {
    for plugin in registry
        .plugins()
        .iter()
        .filter(|plugin| plugin.runs_on(FeatureType::Intergenic))
    {
        let Some(keys) = run_plugin_on_row(plugin, c"null", variant_json, &mut variant.plugin_data)
        else {
            continue;
        };
        for key in keys {
            if !variant.intergenic_row_plugin_keys.contains(&key) {
                variant.intergenic_row_plugin_keys.push(key);
            }
        }
    }
}

/// Runs one plugin on one consequence and variant document pair and merges its
/// result into `fields`, as `run_plugins` does per plugin per line
/// (OutputFactory.pm:2235-2254): a result that is not a JSON object is
/// `Plugin '<name>' did not return a hashref, output ignored!`, a `run` that
/// returns NULL or text that is not JSON is `Plugin '<name>' went wrong:
/// <reason>`, each a warning after which the next plugin runs on the row. The
/// warning sink writes a text once, so a plugin failing the same way on every
/// row warns once (BaseVEP.pm:810-813). Returns the keys the plugin wrote, or
/// `None` when its result was discarded.
fn run_plugin_on_row(
    plugin: &LoadedPlugin,
    consequence_json: &CStr,
    variant_json: &CStr,
    fields: &mut IndexMap<String, String>,
) -> Option<Vec<String>> {
    let name = plugin.name();
    let merged = match plugin.run(consequence_json, variant_json) {
        Ok(result) => merge_plugin_result(fields, name, &result),
        Err(PluginError::Run(reason)) => Err(reason),
        Err(other) => Err(other.to_string()),
    };
    let warning = match merged {
        Ok(Merge::Fields(keys)) => return Some(keys),
        Ok(Merge::NotAnObject) => {
            format!("Plugin '{name}' did not return a hashref, output ignored!")
        }
        Err(reason) => format!("Plugin '{name}' went wrong: {reason}"),
    };
    crate::warnings::global().warn(&warning);
    None
}

fn variant_c_json(variant: &InputVariant) -> anyhow::Result<CString> {
    c_json(variant)
        .with_context(|| format!("serialising variant {} for plugins", variant.allele_string))
}

/// `value` as the NUL-terminated JSON text a plugin reads. serde_json escapes
/// every control character, so the text carries no interior NUL and
/// `CString::new` takes the serialised buffer as it is.
fn c_json<T: serde::Serialize>(value: &T) -> anyhow::Result<CString> {
    Ok(CString::new(serde_json::to_vec(value)?)?)
}

/// What a plugin's result text became: the keys of the fields merged into the
/// row, or nothing because the text was a JSON value other than an object.
#[derive(Debug, PartialEq, Eq)]
pub enum Merge {
    Fields(Vec<String>),
    NotAnObject,
}

/// Merges one plugin's JSON result into `fields`, a consequence's or a variant's
/// `plugin_data`, and returns the keys it wrote.
///
/// A JSON object merges: each entry becomes one Extra-column field, a string
/// value stored verbatim and every other value as its JSON text, so a number
/// stays `0.42` and a nested structure stays parseable. Any other JSON value is
/// [`Merge::NotAnObject`] and merges nothing; text that is not JSON is the
/// error, `plugin '<name>' returned invalid JSON: <reason>`.
pub fn merge_plugin_result(
    fields: &mut IndexMap<String, String>,
    plugin_name: &str,
    result_json: &str,
) -> Result<Merge, String> {
    let value: Value = serde_json::from_str(result_json)
        .map_err(|error| format!("plugin '{plugin_name}' returned invalid JSON: {error}"))?;
    let Value::Object(entries) = value else {
        return Ok(Merge::NotAnObject);
    };
    let mut keys = Vec::with_capacity(entries.len());
    for (key, value) in entries {
        let text = match value {
            Value::String(s) => s,
            other => other.to_string(),
        };
        fields.insert(key.clone(), text);
        keys.push(key);
    }
    Ok(Merge::Fields(keys))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runner::{edit_transcript, set_lengths};
    use vep_core::consequence::TranscriptConsequence;

    #[test]
    fn merge_plugin_result_inserts_strings_verbatim_and_stringifies_the_rest() {
        let mut tc = TranscriptConsequence::default();
        assert_eq!(
            merge_plugin_result(
                &mut tc.plugin_data,
                "echo",
                r#"{"LABEL": "lof", "SCORE": 0.42}"#,
            ),
            Ok(Merge::Fields(vec![
                "LABEL".to_string(),
                "SCORE".to_string()
            ]))
        );
        assert_eq!(tc.plugin_data.len(), 2);
        assert_eq!(tc.plugin_data["LABEL"], "lof");
        assert_eq!(tc.plugin_data["SCORE"], "0.42");
    }

    #[test]
    fn merge_plugin_result_keeps_existing_fields_and_overwrites_a_repeated_key() {
        let mut tc = TranscriptConsequence::default();
        tc.plugin_data
            .insert("REVEL_score".to_string(), "0.9".to_string());
        merge_plugin_result(&mut tc.plugin_data, "echo", r#"{"A": "1"}"#).unwrap();
        merge_plugin_result(&mut tc.plugin_data, "echo", r#"{"A": "2", "B": null}"#).unwrap();
        assert_eq!(tc.plugin_data["REVEL_score"], "0.9");
        assert_eq!(tc.plugin_data["A"], "2");
        assert_eq!(tc.plugin_data["B"], "null");
    }

    /// A JSON value that is not an object merges nothing and is told apart from
    /// text that is not JSON, which names the plugin and the parse error: the
    /// two warnings `run_plugins` writes (OutputFactory.pm:2244, 2253).
    #[test]
    fn merge_plugin_result_tells_a_non_object_from_invalid_json() {
        let mut tc = TranscriptConsequence::default();
        for text in [r#"["not", "an", "object"]"#, "1", "null", r#""text""#] {
            assert_eq!(
                merge_plugin_result(&mut tc.plugin_data, "echo", text),
                Ok(Merge::NotAnObject),
                "{text}"
            );
        }
        assert!(tc.plugin_data.is_empty());

        let reason = merge_plugin_result(&mut tc.plugin_data, "echo", "not json").unwrap_err();
        assert!(
            reason.starts_with("plugin 'echo' returned invalid JSON: expected"),
            "{reason}"
        );
    }

    /// A name that is no path and no library under the directories is returned
    /// as not found, with or without `safe`, and the registry stays empty; the
    /// caller settles it: a name the unported table does not carry is dropped
    /// with the `Failed to load plugin` warning naming the plugin and the
    /// directory searched, the plugin component's error under `safe`
    /// (Runner.pm:929-939, thrown under `--safe`).
    #[test]
    fn an_unknown_name_is_returned_not_found_and_the_caller_skips_it_or_dies_under_safe() {
        // The `Failed to load plugin` warning goes to the process-wide sink.
        let _guard = crate::warnings::sink_test_guard();
        let dirs = PluginDirs::expand("/no/such/dir");
        let mut registry = PluginRegistry::new();
        for safe in [false, true] {
            let outcome =
                load_unresolved(&mut registry, "NoSuchPlugin,param", &dirs, safe).unwrap();
            assert_eq!(outcome, DylibOutcome::NotFound);
        }
        assert!(registry.is_empty());

        let dir = tempfile::tempdir().unwrap();
        let dirs = PluginDirs::expand(dir.path().to_str().unwrap());
        let args = ["NoSuchPlugin,param".to_string()];
        let context = vep_builtins::traits::RunContext::default();
        let loaded =
            crate::annotator::load_plugins(&args, &dirs, false, None, None, context).unwrap();
        assert!(loaded.dylib.is_empty() && loaded.builtin.is_empty());
        let err = crate::annotator::load_plugins(&args, &dirs, true, None, None, context)
            .err()
            .expect("an unknown name is fatal under --safe");
        assert_eq!(
            err.downcast_ref::<CliError>().unwrap().to_string(),
            format!(
                "ERROR: [plugin] Failed to load plugin NoSuchPlugin: NoSuchPlugin is not a built-in plugin and no libNoSuchPlugin.so or libNoSuchPlugin.dylib is under {}",
                dir.path().display()
            )
        );
    }

    /// A `.pm` file, or a file whose first non-blank line begins `package `, is
    /// an Ensembl VEP plugin module; a comment-first file, a non-Perl file, a
    /// directory and a path that does not exist are not. A module given by path
    /// is refused before any open, whatever `safe` says.
    #[test]
    fn a_perl_module_is_recognised_by_suffix_or_package_line_and_refused() {
        let dir = tempfile::tempdir().unwrap();
        let by_suffix = dir.path().join("Foo.pm");
        std::fs::write(&by_suffix, "=head1 LICENSE\n").unwrap();
        let by_package = dir.path().join("Foo");
        std::fs::write(&by_package, "\n\n  package Foo;\n").unwrap();
        let comment_first = dir.path().join("script");
        std::fs::write(&comment_first, "# package Foo;\npackage Foo;\n").unwrap();
        assert!(is_perl_module(&by_suffix));
        assert!(is_perl_module(&by_package));
        assert!(!is_perl_module(&comment_first));
        assert!(!is_perl_module(Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/Cargo.toml"
        ))));
        assert!(!is_perl_module(dir.path()));
        assert!(!is_perl_module(&dir.path().join("absent.pm")));

        let dirs = PluginDirs::expand(dir.path().to_str().unwrap());
        assert_eq!(perl_module_under(&dirs, "Foo"), Some(by_suffix.clone()));
        assert_eq!(perl_module_under(&dirs, "Bar"), None);

        let mut registry = PluginRegistry::new();
        for safe in [false, true] {
            let arg = format!("{},param", by_suffix.display());
            let err = load_unresolved(&mut registry, &arg, &dirs, safe).unwrap_err();
            assert_eq!(
                err.downcast_ref::<CliError>().unwrap().to_string(),
                format!(
                    "ERROR: [plugin] {} is an Ensembl VEP Perl plugin; vep-rs runs no Perl plugin, see docs/plugins.md",
                    by_suffix.display()
                )
            );
        }
        assert!(registry.is_empty());
    }

    /// A file that exists but is no shared library is dropped as a module that
    /// fails to compile is (Runner.pm:929-939: `Failed to compile plugin
    /// <Name>: <reason>`, a warning and `next`, thrown under `--safe`), the
    /// reason the loader's own text naming the file.
    #[test]
    fn load_unresolved_drops_a_file_that_is_not_a_plugin_and_dies_under_safe() {
        // The `Failed to compile plugin` warning goes to the process-wide sink.
        let _guard = crate::warnings::sink_test_guard();
        let not_a_dylib = concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml");
        let dirs = PluginDirs::expand("/no/such/dir");
        let mut registry = PluginRegistry::new();
        load_unresolved(&mut registry, not_a_dylib, &dirs, false).unwrap();
        assert!(registry.is_empty());
        let err = load_unresolved(&mut registry, not_a_dylib, &dirs, true).unwrap_err();
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
    /// other: the transcript's values under their accessor names beside the
    /// row's fields, with `cdna_total`, `cds_total`, `protein_total`,
    /// `plugin_data`, `loftee_ctx`, `cds_length`, `source_cache`,
    /// `bam_edit_status`, `own_span_positions`, `used_ref`, `shifted_allele`,
    /// `refseq_match`, `refseq_offset`, `motif` and `cell_type` present only
    /// when populated.
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

        let mut tc = TranscriptConsequence::default();
        tc.id_rank = 7;
        tc.record_span_in_cds = true;
        assert_eq!(sorted_keys(&tc), expected);

        tc.plugin_data
            .insert("REVEL_score".to_string(), "0.9".to_string());
        tc.loftee_ctx = Some(Box::new(
            vep_core::consequence::LofteeContext::from_transcript(
                &crate::runner::make_runner_test_transcript("ENST1", "ENSG1", "1", 100),
            ),
        ));
        edit_transcript(&mut tc, |transcript| {
            set_lengths(transcript, Some(300), 2500);
            transcript.source_cache = Some(vep_core::consequence::SourceCache::Ensembl);
        });
        tc.bam_edit_status = Some(vep_core::transcript::BamEditStatus::Ok);
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
        tc.motif = Some(Box::new(vep_core::consequence::MotifFields::new(
            std::sync::Arc::from("ENSPFM0290"),
            Some(7),
            Some(true),
            Some(-0.089),
            std::sync::Arc::from([std::sync::Arc::from("HNF4A")]),
        )));
        tc.cell_type = Some(Box::new([vep_core::consequence::CellTypeActivity {
            name: std::sync::Arc::from("HUVEC"),
            activity: vep_core::consequence::Activity::Repressed,
        }]));
        expected.extend(
            [
                "plugin_data",
                "loftee_ctx",
                "cds_length",
                "source_cache",
                "bam_edit_status",
                "cdna_total",
                "cds_total",
                "protein_total",
                "own_span_positions",
                "used_ref",
                "shifted_allele",
                "motif",
                "cell_type",
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

    /// One run whose plugin documents a committed file holds.
    struct DocumentRun {
        /// The stem of the run's two files under
        /// `tests/fixtures/plugin-documents`: `<name>.vcf`, the input records,
        /// and `<name>.jsonl`, the documents.
        name: &'static str,
        /// The directory under `tests/golden/116` that `{corpus}` names in
        /// `flags`.
        corpus: &'static str,
        /// The command line after `vep --offline --fork 1 --species
        /// homo_sapiens --assembly GRCh38`.
        flags: &'static [&'static str],
    }

    /// The runs, each over records that fill keys of the two documents no
    /// other run fills: an Ensembl cache under the flags of `--everything`
    /// that read no known variant, with HGVS, SPDI, `--nearest` and a 3'
    /// shift; a RefSeq cache and a merged cache; the LoFTEE built-in;
    /// regulatory and motif rows with cell types, on sequence and structural
    /// records; known variants; and structural records under `--check_svs`.
    const DOCUMENT_RUNS: &[DocumentRun] = &[
        DocumentRun {
            name: "ensembl",
            corpus: "GRCh38-everything-fasta",
            flags: &[
                "--json_cache",
                "{corpus}/json_cache",
                "--fasta",
                "{corpus}/reference.fa.gz",
                "--sift",
                "b",
                "--polyphen",
                "b",
                "--ccds",
                "--hgvs",
                "--hgvsg",
                "--spdi",
                "--ga4gh_vrs",
                "--symbol",
                "--numbers",
                "--domains",
                "--regulatory",
                "--canonical",
                "--protein",
                "--biotype",
                "--uniprot",
                "--mane",
                "--tsl",
                "--appris",
                "--variant_class",
                "--gene_phenotype",
                "--mirna",
                "--total_length",
                "--nearest",
                "symbol",
                "--shift_3prime",
                "1",
            ],
        },
        DocumentRun {
            name: "refseq",
            corpus: "GRCh38-refseq",
            flags: &[
                "--cache",
                "--dir_cache",
                "{corpus}/cache_root",
                "--cache_version",
                "116",
                "--refseq",
                "--fasta",
                "{corpus}/reference.fa.gz",
                "--hgvs",
                "--numbers",
                "--total_length",
            ],
        },
        DocumentRun {
            name: "merged",
            corpus: "GRCh38-merged",
            flags: &[
                "--json_cache",
                "{corpus}/json_cache",
                "--merged",
                "--fasta",
                "{corpus}/reference.fa.gz",
                "--hgvs",
                "--numbers",
                "--total_length",
            ],
        },
        DocumentRun {
            name: "loftee",
            corpus: "GRCh38-plugin-loftee",
            flags: &[
                "--json_cache",
                "{corpus}/json_cache",
                "--fasta",
                "{corpus}/reference.fa.gz",
                "--dir_plugins",
                "{corpus}/plugin_data",
                "--plugin",
                "LoF,loftee_path:{corpus}/plugin_data,human_ancestor_fa:{corpus}/plugin_data/human_ancestor.fa.gz,gerp_bigwig:{corpus}/plugin_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw",
            ],
        },
        DocumentRun {
            name: "regulatory",
            corpus: "GRCh38-regulatory",
            flags: &[
                "--json_cache",
                "{corpus}/json_cache",
                "--fasta",
                "{corpus}/reference.fa.gz",
                "--regulatory",
                "--cell_type",
                "GM12878_(f),H1_(m),astrocyte",
            ],
        },
        DocumentRun {
            name: "colocated",
            corpus: "GRCh38-colocated",
            flags: &["--json_cache", "{corpus}/json_cache", "--check_existing"],
        },
        DocumentRun {
            name: "structural",
            corpus: "GRCh38-checksvs",
            flags: &[
                "--json_cache",
                "{corpus}/json_cache",
                "--check_svs",
                "{corpus}/structural_variants.vcf.gz",
            ],
        },
    ];

    /// The command that rewrites the committed documents.
    const REGENERATE: &str =
        "VEP_UPDATE_PLUGIN_DOCUMENTS=1 cargo test -p vep-cli --lib plugin_documents";

    /// Where `written` first departs from `committed`: the number of the line,
    /// counted from 1, and the length in bytes of that line in each; `None`
    /// for the same bytes. A line itself is never returned, since a document
    /// can hold sequence.
    fn first_difference(committed: &[u8], written: &[u8]) -> Option<(usize, usize, usize)> {
        if committed == written {
            return None;
        }
        let mut committed_lines = committed.split(|byte| *byte == b'\n');
        let mut written_lines = written.split(|byte| *byte == b'\n');
        let mut number = 1;
        loop {
            match (committed_lines.next(), written_lines.next()) {
                (Some(committed_line), Some(written_line)) if committed_line == written_line => {
                    number += 1;
                }
                (committed_line, written_line) => {
                    return Some((
                        number,
                        committed_line.map_or(0, <[u8]>::len),
                        written_line.map_or(0, <[u8]>::len),
                    ));
                }
            }
        }
    }

    /// The two documents a plugin receives are, byte for byte, the committed
    /// ones on the records of each run of `DOCUMENT_RUNS`: the library's
    /// annotator, built from the run's command line, annotates each record of
    /// `<name>.vcf` as one batch, and every variant of the batch is written
    /// through `variant_c_json` and each of its consequences through `c_json`,
    /// the two calls `run` makes, one document a line in batch order, the
    /// variant's line before its consequences'. No run loads a shared library
    /// or a `--custom` source, so the batch the annotator returns is the batch
    /// `run` is handed.
    ///
    /// The variant document writes the frequencies of a known variant
    /// (`frequency_strings`) in the iteration order of a `HashMap`, which
    /// differs from one process to the next once the map holds two
    /// populations. The `colocated` run, the one run that reads known
    /// variants, therefore takes records whose known variants carry no
    /// frequency, and the bytes of a populated map are held by no run.
    ///
    /// `VEP_UPDATE_PLUGIN_DOCUMENTS=1` rewrites the committed files instead of
    /// comparing them.
    #[test]
    fn plugin_documents_of_the_corpus_runs_equal_the_committed_bytes() {
        let _guard = crate::warnings::sink_test_guard();
        let update = std::env::var_os("VEP_UPDATE_PLUGIN_DOCUMENTS").is_some();
        let fixtures =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/plugin-documents");
        let mut differing = Vec::new();
        for run in DOCUMENT_RUNS {
            let corpus = format!(
                "{}/../../tests/golden/116/{}",
                env!("CARGO_MANIFEST_DIR"),
                run.corpus
            );
            let mut argv: Vec<String> = [
                "vep",
                "--offline",
                "--fork",
                "1",
                "--species",
                "homo_sapiens",
                "--assembly",
                "GRCh38",
            ]
            .map(str::to_string)
            .to_vec();
            argv.extend(
                run.flags
                    .iter()
                    .map(|flag| flag.replace("{corpus}", &corpus)),
            );
            let config = crate::config::Config::from_argv(argv).unwrap();
            let annotator = crate::annotator::Annotator::from_config(&config).unwrap();

            let input = std::fs::read_to_string(fixtures.join(format!("{}.vcf", run.name)))
                .unwrap_or_else(|error| panic!("{}.vcf: {error}", run.name));
            let records = input
                .lines()
                .filter(|line| !line.is_empty() && !line.starts_with('#'));
            let mut documents = Vec::new();
            for (index, line) in records.enumerate() {
                let record = index + 1;
                let mut batch = crate::vcf_parser::parse_vcf_line(line, true)
                    .unwrap_or_else(|_| panic!("{}: record {record} does not parse", run.name));
                annotator.annotate_batch(&mut batch);
                for variant in &batch {
                    let document = variant_c_json(variant).unwrap_or_else(|_| {
                        panic!(
                            "{}: a variant of record {record} does not serialise",
                            run.name
                        )
                    });
                    documents.extend_from_slice(document.as_bytes());
                    documents.push(b'\n');
                    for tc in &variant.transcript_consequences {
                        let document = c_json(tc).unwrap_or_else(|_| {
                            panic!(
                                "{}: a consequence of record {record} does not serialise",
                                run.name
                            )
                        });
                        documents.extend_from_slice(document.as_bytes());
                        documents.push(b'\n');
                    }
                }
            }
            assert!(!documents.is_empty(), "{}.vcf holds no record", run.name);

            let path = fixtures.join(format!("{}.jsonl", run.name));
            if update {
                std::fs::write(&path, &documents).unwrap();
                continue;
            }
            let committed = std::fs::read(&path).unwrap_or_default();
            if let Some((line, committed_len, written_len)) =
                first_difference(&committed, &documents)
            {
                differing.push(format!(
                    "{}: line {line} is {committed_len} bytes committed and {written_len} bytes \
                     written, the first line that differs",
                    path.display()
                ));
            }
        }
        assert!(
            differing.is_empty(),
            "{}\nregenerate with `{REGENERATE}`",
            differing.join("\n")
        );
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
