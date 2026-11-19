// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A file that is no bigWig handed to a bigWig reader, through the `vep`
//! binary: the Conservation plugin, LoFTEE's GERP source and a `--custom`
//! bigWig source. Standard error carries vep-rs's own refusal and the run's
//! status lines alone, and under `--quiet` nothing but a `--safe` error: the
//! bigWig library's path-taking opener prints its own `Error when opening`
//! line beside the error it returns, and vep-rs owns its standard error, so
//! every reader opens the file itself and hands the library the open handle.

mod common;

use std::path::{Path, PathBuf};
use std::process::Output;

/// A short text file is no bigWig: the header read ends before the magic
/// number, so the open fails with the reader's end-of-file error.
const REASON: &str = "failed to fill whole buffer";

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{name}"))
}

/// Writes the file that is no bigWig into `dir`, named as a bigWig.
fn text_bigwig(dir: &Path) -> String {
    let path = dir.join("scores.bw");
    std::fs::write(&path, "not a bigWig\n").unwrap();
    path.display().to_string()
}

/// Runs the binary on `corpus`'s records with `extra` after the fixed
/// arguments, writing to `out`; the corpus's reference FASTA is passed when the
/// corpus carries one.
fn run(corpus: &Path, out: &Path, extra: &[&str]) -> Output {
    let mut command = common::isolated_vep_anywhere();
    command
        .arg("-i")
        .arg(corpus.join("variants.vcf"))
        .arg("-o")
        .arg(out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats"]);
    let fasta = corpus.join("reference.fa.gz");
    if fasta.exists() {
        command.arg("--fasta").arg(fasta);
    }
    command.args(extra).output().expect("vep binary runs")
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn warnings(out: &Path) -> String {
    let path = out.with_file_name(format!(
        "{}_warnings.txt",
        out.file_name().unwrap().to_string_lossy()
    ));
    std::fs::read_to_string(&path).unwrap_or_default()
}

/// The run's own status lines: the banner, the `tracing` status messages and
/// the summary line. Anything else on standard error is a warning or an error.
fn is_status_line(line: &str) -> bool {
    line.starts_with("--------------------")
        || line.starts_with("Ensembl VEP (Rust)")
        || line.contains("INFO")
        || line.starts_with("Processed ")
}

/// A plugin given the text file: under `--quiet` standard error is empty and
/// the refusal goes to the warnings file; without `--quiet` the refusal is the
/// one line beside the status lines; under `--quiet --safe` the error line is
/// all standard error carries.
fn plugin_refusal_alone_on_standard_error(corpus_name: &str, plugin: &str, refusal: &str) {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let corpus = corpus(corpus_name);
    let data = corpus.join("plugin_data");
    let text_bw = text_bigwig(dir);
    let plugin = plugin
        .replace("{plugins}", &data.display().to_string())
        .replace("{text}", &text_bw);
    let refusal = refusal
        .replace("{text}", &text_bw)
        .replace("{reason}", REASON);
    let plugin_args = ["--dir_plugins", data.to_str().unwrap(), "--plugin", &plugin];

    let out = dir.join("quiet.txt");
    let quiet = run(&corpus, &out, &[&plugin_args[..], &["--quiet"]].concat());
    assert_eq!(quiet.status.code(), Some(0), "{}", text(&quiet.stderr));
    assert_eq!(
        text(&quiet.stderr),
        "",
        "--quiet leaves standard error empty"
    );
    assert_eq!(
        warnings(&out).lines().collect::<Vec<_>>(),
        [format!("WARNING: {refusal}")],
        "the refusal is the warnings file's one line"
    );
    assert!(
        std::fs::read_to_string(&out)
            .unwrap()
            .lines()
            .any(|line| !line.starts_with('#')),
        "the run annotates without the plugin"
    );

    let out = dir.join("loud.txt");
    let loud = run(&corpus, &out, &plugin_args);
    assert_eq!(loud.status.code(), Some(0), "{}", text(&loud.stderr));
    let stderr = text(&loud.stderr);
    let beside_status: Vec<&str> = stderr.lines().filter(|l| !is_status_line(l)).collect();
    assert_eq!(
        beside_status,
        [format!("WARNING: {refusal}")],
        "the refusal is the one line beside the status lines:\n{stderr}"
    );

    let out = dir.join("safe.txt");
    let safe = run(
        &corpus,
        &out,
        &[&plugin_args[..], &["--quiet", "--safe"]].concat(),
    );
    assert_eq!(safe.status.code(), Some(1), "{}", text(&safe.stderr));
    assert_eq!(
        text(&safe.stderr),
        format!("ERROR: [plugin] {refusal}\n"),
        "the error line is all standard error carries"
    );
}

#[test]
fn conservation_text_file_draws_the_refusal_alone() {
    plugin_refusal_alone_on_standard_error(
        "GRCh38-plugin-conservation",
        "Conservation,{text}",
        "Failed to instantiate plugin Conservation: Conservation: {text} is not a readable \
         bigWig: {reason}",
    );
}

#[test]
fn loftee_gerp_text_file_draws_the_refusal_alone() {
    plugin_refusal_alone_on_standard_error(
        "GRCh38-plugin-loftee",
        "LoF,loftee_path:{plugins},human_ancestor_fa:{plugins}/human_ancestor.fa.gz,\
         gerp_bigwig:{text}",
        "Failed to instantiate plugin LoF: LoFTEE: failed to open GERP bigWig '{text}': \
         {reason}",
    );
}

/// A `--custom` bigWig source given the text file is refused before output,
/// under `--quiet` with the error line alone on standard error.
#[test]
fn custom_bigwig_text_file_draws_the_refusal_alone() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let corpus = corpus("GRCh38-custom");
    let text_bw = text_bigwig(dir);
    let source = format!("file={text_bw},short_name=cons,format=bigwig,type=overlap");
    let refusal = format!("ERROR: [custom] cannot open the --custom source {text_bw}: {REASON}\n");

    let out = dir.join("quiet.txt");
    let quiet = run(&corpus, &out, &["--custom", &source, "--quiet"]);
    assert_eq!(quiet.status.code(), Some(1), "{}", text(&quiet.stderr));
    assert_eq!(
        text(&quiet.stderr),
        refusal,
        "the error line is all standard error carries"
    );
    assert!(!out.exists(), "nothing is written");

    let out = dir.join("loud.txt");
    let loud = run(&corpus, &out, &["--custom", &source]);
    assert_eq!(loud.status.code(), Some(1), "{}", text(&loud.stderr));
    let stderr = text(&loud.stderr);
    let beside_status: Vec<&str> = stderr.lines().filter(|l| !is_status_line(l)).collect();
    assert_eq!(
        beside_status,
        [refusal.trim_end()],
        "the error is the one line beside the status lines:\n{stderr}"
    );
}
