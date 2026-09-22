// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The dylib plugin path through the `vep` binary on the pruned golden corpus: an
//! unresolvable `--plugin` name is warned about and skipped, a file that is not a
//! plugin is a hard error, and the `echo_plugin` example resolved under
//! `--dir_plugins` lands its fields in the Extra column of every transcript row.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const SKIPPED_WARNING: &str =
    "plugin 'NoSuchPlugin' not found as a built-in or as a dylib under --dir_plugins; skipped";

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/115/GRCh37")
}

/// Runs `vep` in the default output format on the corpus without `--quiet`, so
/// warnings reach stderr.
fn vep_on_corpus(out: &Path, extra: &[&str]) -> Output {
    let corpus = corpus_dir();
    Command::new(env!("CARGO_BIN_EXE_vep"))
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

#[test]
fn a_file_that_is_not_a_plugin_is_a_hard_error() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("default.txt");
    let not_a_plugin = Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
    let output = vep_on_corpus(&out, &["--plugin", not_a_plugin.to_str().unwrap()]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success(), "{stderr}");
    assert!(stderr.contains("failed to load"), "{stderr}");
    assert!(stderr.contains("Cargo.toml"), "{stderr}");
    assert!(
        !out.exists(),
        "no output is written when a named plugin fails"
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
