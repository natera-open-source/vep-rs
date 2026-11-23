// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The line of a built-in plugin that fails to load, through the `vep`
//! binary. Ensembl VEP writes `Failed to instantiate plugin <module>: <the
//! constructor's text>` (Runner.pm:950), and no text its OpenTargets module
//! dies with names the module (OpenTargets.pm:110, :118, :227, :244, :258), so
//! its line names the plugin once. A built-in plugin opens the texts it words
//! with its own name, and the line carries the text after that name: a
//! warning in the warnings file, and the plugin component's error under
//! `--safe`.

use std::path::{Path, PathBuf};
use std::process::Output;

mod common;

fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-plugin-opentargets")
}

/// Runs the binary on the corpus's records with one `--plugin` value under
/// `--quiet`, so standard error holds an error alone and the warnings go to
/// the file beside the output.
fn run(out: &Path, plugin: &str, more: &[&str]) -> Output {
    let corpus = corpus();
    common::isolated_vep_anywhere()
        .arg("-i")
        .arg(corpus.join("variants.vcf"))
        .arg("-o")
        .arg(out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .arg("--dir_plugins")
        .arg(&corpus)
        .args(["--plugin", plugin])
        .args(more)
        .output()
        .expect("vep binary runs")
}

fn warnings(out: &Path) -> String {
    let path = out.with_file_name(format!(
        "{}_warnings.txt",
        out.file_name().unwrap().to_string_lossy()
    ));
    std::fs::read_to_string(&path).unwrap_or_default()
}

/// One refused `--plugin` value: the run ends 0, writes its output, writes
/// nothing to standard error and holds `line` as its whole warnings file.
fn dropped_with(plugin: &str, line: &str) {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("out.txt");
    let output = run(&out, plugin, &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(0), "{plugin}: {stderr}");
    assert_eq!(stderr, "", "{plugin}: standard error under --quiet");
    assert!(out.exists(), "{plugin}: the run writes its output");
    assert_eq!(warnings(&out), format!("{line}\n"), "{plugin}");
}

/// Each refusal the OpenTargets plugin words for a file it cannot use names the
/// plugin once.
#[test]
fn a_failed_open_targets_load_names_the_plugin_once() {
    dropped_with(
        "OpenTargets",
        "WARNING: Failed to instantiate plugin OpenTargets: no file specified; give \
         file=<path to the Open Targets export>",
    );
    dropped_with(
        "OpenTargets,file=/no/such/export.tsv.gz",
        "WARNING: Failed to instantiate plugin OpenTargets: data file /no/such/export.tsv.gz \
         not found",
    );
    let export = corpus().join("open_targets_vep.slice.tsv.gz");
    dropped_with(
        &format!("OpenTargets,file={},cutoff=1", export.display()),
        "WARNING: Failed to instantiate plugin OpenTargets: unknown parameter 'cutoff'; \
         OpenTargets accepts file, cols, study_types, lead_variants_only",
    );
}

/// Under `--safe` the same line is the plugin component's error: exit 1, the
/// line alone on standard error, no output file and no warnings file.
#[test]
fn a_failed_open_targets_load_under_safe_names_the_plugin_once() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("out.txt");
    let output = run(&out, "OpenTargets", &["--safe"]);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "ERROR: [plugin] Failed to instantiate plugin OpenTargets: no file specified; give \
         file=<path to the Open Targets export>\n"
    );
    assert!(!out.exists(), "a refused run writes no output");
    assert_eq!(warnings(&out), "", "a refused run writes no warnings file");
}

/// The line is the loader's, so every built-in plugin whose text opens with
/// its own name is named once: two more plugins, each on a refusal of its own.
#[test]
fn other_built_in_plugins_are_named_once_as_well() {
    dropped_with(
        "LoFtool,/no/such/scores.txt",
        "WARNING: Failed to instantiate plugin LoFtool: scores file /no/such/scores.txt not found",
    );
    dropped_with(
        "CADD",
        "WARNING: Failed to instantiate plugin CADD: no CADD files specified; give snv=, \
         indels= or sv= a path, or list the files by position",
    );
}

/// A text that opens with another name than the `--plugin` value's, or with no
/// name, is written whole: the LoF module's texts name LoFTEE, and the refusal
/// of a remote data file names no plugin.
#[test]
fn a_text_that_does_not_open_with_the_plugin_name_is_written_whole() {
    dropped_with(
        "LoF,run_splice_predictions=1",
        "WARNING: Failed to instantiate plugin LoF: LoFTEE: run_splice_predictions=1 asks for \
         the splice-prediction layer, which vep-rs does not implement",
    );
    dropped_with(
        "REVEL,https://example.invalid/x.tsv.gz",
        "WARNING: Failed to instantiate plugin REVEL: remote data files are not supported; \
         download https://example.invalid/x.tsv.gz and index it locally",
    );
}
