// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The dylib plugin path through the `vep` binary on the pruned golden corpus: an
//! unresolvable `--plugin` name and a file that is not a plugin are each warned
//! about and skipped, the file an error under `--safe`, and the `echo_plugin`
//! example resolved under `--dir_plugins` lands its fields in the Extra column of
//! every transcript row.

use std::path::{Path, PathBuf};
use std::process::Output;

mod common;

/// The warning an unresolvable name writes, up to the plugins directory it names.
const SKIPPED_WARNING: &str =
    "WARNING: Failed to load plugin NoSuchPlugin: NoSuchPlugin is not a built-in plugin and no libNoSuchPlugin.so or libNoSuchPlugin.dylib is under ";

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37")
}

/// Runs `vep` in the default output format on the corpus without `--quiet`, so
/// warnings reach stderr.
fn vep_on_corpus(out: &Path, extra: &[&str]) -> Output {
    let corpus = corpus_dir();
    common::isolated_vep_anywhere()
        .arg("-i")
        .arg(corpus.join("variants.vcf"))
        .arg("-o")
        .arg(out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh37"])
        .args(["--force_overwrite", "--no_stats"])
        .args(extra)
        .output()
        .expect("vep binary runs")
}

fn data_lines(out: &Path) -> Vec<String> {
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

#[test]
fn unresolvable_plugin_name_is_warned_and_skipped() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let output = vep_on_corpus(&out, &["--plugin", "NoSuchPlugin,param1"]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(stderr.contains(SKIPPED_WARNING), "{stderr}");
    assert!(!data_lines(&out).is_empty(), "the run still annotates");
}

#[test]
fn unresolvable_plugin_name_is_skipped_with_a_plugin_dir_too() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let output = vep_on_corpus(
        &out,
        &["--plugin", "NoSuchPlugin", "--dir_plugins", "/no/such/dir"],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(stderr.contains(SKIPPED_WARNING), "{stderr}");
}

/// A file that is not a shared library is dropped with Ensembl VEP's warning for
/// a plugin module that fails to compile (Runner.pm:929-939) and the run goes on;
/// under `--safe` it is the plugin component's error and no output is written.
#[test]
fn a_file_that_is_not_a_plugin_is_dropped_and_fatal_under_safe() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let not_a_plugin = Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let not_a_plugin = not_a_plugin.to_str().unwrap();
    let plugins_dir = corpus_dir();
    let plugins_dir = plugins_dir.to_str().unwrap();
    let output = vep_on_corpus(
        &out,
        &["--plugin", not_a_plugin, "--dir_plugins", plugins_dir],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(
        stderr.contains(&format!(
            "WARNING: Failed to compile plugin {not_a_plugin}: {not_a_plugin}: "
        )),
        "{stderr}"
    );
    assert!(!data_lines(&out).is_empty(), "the run still annotates");
    std::fs::remove_file(&out).unwrap();

    let output = vep_on_corpus(
        &out,
        &[
            "--plugin",
            not_a_plugin,
            "--dir_plugins",
            plugins_dir,
            "--safe",
        ],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.contains(&format!(
            "ERROR: [plugin] Failed to compile plugin {not_a_plugin}: {not_a_plugin}: "
        )),
        "{stderr}"
    );
    assert!(
        !out.exists(),
        "no output is written when --safe refuses a plugin"
    );
}

/// The profile's `examples/` directory when the `echo_plugin` example dylib has
/// been built there. `cargo test` on a selection that includes `vep-plugin` (the
/// workspace, or `-p vep-plugin -p vep-cli`) builds it; `-p vep-cli` alone does not.
fn echo_plugin_dir() -> Option<PathBuf> {
    let exe = std::env::current_exe().unwrap();
    let examples = exe.parent().unwrap().parent().unwrap().join("examples");
    let ext = if cfg!(target_os = "macos") {
        "dylib"
    } else {
        "so"
    };
    let dylib = examples.join(format!("libecho_plugin.{ext}"));
    if dylib.is_file() {
        return Some(examples);
    }
    // Under CI the whole workspace is under test, so an absent example is a build
    // defect rather than a narrower selection.
    assert!(
        std::env::var_os("CI").is_none(),
        "echo_plugin example dylib not found at {}",
        dylib.display()
    );
    eprintln!(
        "Skipping: echo_plugin example dylib not found at {}",
        dylib.display()
    );
    None
}

#[test]
fn echo_plugin_under_dir_plugins_annotates_every_transcript_row() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let output = vep_on_corpus(
        &out,
        &[
            "--plugin",
            "echo_plugin,hello",
            "--dir_plugins",
            dir.to_str().unwrap(),
        ],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "exit {:?}\n{stderr}",
        output.status.code()
    );
    assert!(!stderr.contains("not found"), "{stderr}");

    let lines = data_lines(&out);
    let transcript_rows: Vec<Vec<&str>> = lines
        .iter()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|cols| cols[5] == "Transcript")
        .collect();
    assert!(!transcript_rows.is_empty(), "corpus has transcript rows");
    for cols in &transcript_rows {
        let (feature, extra) = (cols[4], cols[13]);
        let fields: Vec<&str> = extra.split(';').collect();
        assert!(
            fields.contains(&format!("ECHO_FEATURE={feature}").as_str()),
            "row {} lacks ECHO_FEATURE={feature}: {extra}",
            cols[0]
        );
        assert!(
            fields.contains(&"ECHO_PARAM=hello"),
            "row {} lacks ECHO_PARAM=hello: {extra}",
            cols[0]
        );
    }
}

/// A plugin runs on the rows whose feature type it declares, `Transcript` alone
/// without a `vep_plugin_feature_types` symbol or list, as Ensembl VEP runs a
/// plugin only on the feature classes its `feature_types` names
/// (`OutputFactory.pm:2223`, `BaseVepPlugin.pm` `_check_types`); and a plugin
/// that declares regulatory or motif features switches regulatory annotation
/// on for the run without `--regulatory` (`Runner.pm:998-1002`).
#[test]
fn plugins_run_on_the_rows_their_feature_types_declare() {
    let Some(dir) = echo_plugin_dir() else {
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let cache = common::regulatory_cache(tmp.path());
    let input = tmp.path().join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tinside\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    let run = |plugin: &str, regulatory: bool| -> Vec<(String, String, String)> {
        let out = tmp
            .path()
            .join(format!("out-{}-{regulatory}.txt", plugin.len()));
        let mut command = common::isolated_vep(tmp.path());
        command
            .arg("-i")
            .arg(&input)
            .arg("-o")
            .arg(&out)
            .args(["--json_cache", cache.to_str().unwrap()])
            .args(["--assembly", "GRCh38", "--no_fasta", "--no_stats"])
            .args(["--force_overwrite", "--plugin", plugin])
            .args(["--dir_plugins", dir.to_str().unwrap()]);
        if regulatory {
            command.arg("--regulatory");
        }
        let output = command.output().unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{stderr}");
        assert!(!stderr.contains("Failed to"), "{stderr}");
        data_lines(&out)
            .iter()
            .map(|line| {
                let cols: Vec<&str> = line.split('\t').collect();
                let echo = cols[13]
                    .split(';')
                    .find_map(|f| f.strip_prefix("ECHO_FEATURE="))
                    .unwrap_or("-")
                    .to_string();
                (cols[4].to_string(), cols[5].to_string(), echo)
            })
            .collect()
    };
    let rows = |echo_regulatory: &str, echo_motif: &str| {
        vec![
            (
                "ENSR_FIXTURE".to_string(),
                "RegulatoryFeature".to_string(),
                echo_regulatory.to_string(),
            ),
            (
                "ENSM_FORWARD".to_string(),
                "MotifFeature".to_string(),
                echo_motif.to_string(),
            ),
            ("-".to_string(), "-".to_string(), "-".to_string()),
        ]
    };
    // No declaration: transcript rows alone, and the record has none.
    assert_eq!(run("echo_plugin", true), rows("-", "-"));
    // `Feature` is every feature class.
    assert_eq!(
        run("echo_plugin,feature_types=Feature", true),
        rows("ENSR_FIXTURE", "ENSM_FORWARD")
    );
    // One class, declared without --regulatory, switches the annotation on.
    assert_eq!(
        run("echo_plugin,feature_types=MotifFeature", false),
        rows("-", "ENSM_FORWARD")
    );
    // Without a plugin that names them and without the flag, no regulatory row:
    // `Feature` covers the classes but names neither word (Runner.pm:998).
    for plugin in ["echo_plugin", "echo_plugin,feature_types=Feature"] {
        assert_eq!(
            run(plugin, false),
            vec![("-".to_string(), "-".to_string(), "-".to_string())],
            "{plugin}"
        );
    }
}
