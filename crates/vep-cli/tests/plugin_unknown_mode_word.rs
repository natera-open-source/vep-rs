// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A mode word outside a plugin's set, through the `vep` binary: LOEUF's
//! `match_by`, pLI's second parameter and a word after Conservation's `MAX`.
//! Ensembl VEP stores the word, loads the plugin and writes its header line
//! with no value on any row (LOEUF.pm:146-149, :207-209; pLI.pm:109-111,
//! :218-227) or the other method's values (Conservation.pm:82-89); vep-rs
//! refuses the plugin when it loads, the documented class
//! `plugin_unknown_mode_word_refused`. The ordinary word of each plugin runs
//! beside it, so the fixture is shown to reach the plugin.

use std::path::{Path, PathBuf};
use std::process::Output;

mod common;

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{name}"))
}

/// Runs the binary on a corpus's records with one `--plugin` value, `--quiet`,
/// so the warnings go to the file beside the output alone.
fn run(corpus: &Path, out: &Path, plugin: &str) -> Output {
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
        .arg(corpus)
        .args(["--plugin", plugin])
        .output()
        .expect("vep binary runs")
}

/// The output's `## <field> :` header lines and the data rows carrying
/// `<field>=`, out of all its data rows.
fn header_and_cells(out: &Path, field: &str) -> (usize, usize, usize) {
    let text = std::fs::read_to_string(out).unwrap();
    let header = text
        .lines()
        .filter(|line| line.starts_with(&format!("## {field} : ")))
        .count();
    let rows: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
    let cells = rows
        .iter()
        .filter(|row| {
            row.split('\t').nth(13).is_some_and(|extra| {
                extra
                    .split(';')
                    .any(|item| item.starts_with(&format!("{field}=")))
            })
        })
        .count();
    (header, cells, rows.len())
}

fn warnings(out: &Path) -> String {
    let path = out.with_file_name(format!(
        "{}_warnings.txt",
        out.file_name().unwrap().to_string_lossy()
    ));
    std::fs::read_to_string(&path).unwrap_or_default()
}

/// One plugin: its ordinary word loads it and writes the field; the word
/// outside its set drops it with the named refusal, the rows written alike.
fn refused_where_ensembl_vep_loads(
    corpus_name: &str,
    field: &str,
    ordinary: &str,
    outside: &str,
    refusal: &str,
) {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let corpus = corpus(corpus_name);
    let data = corpus.to_string_lossy();

    let loaded = dir.join("ordinary.txt");
    let output = run(&corpus, &loaded, &ordinary.replace("{corpus}", &data));
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (header, cells, rows) = header_and_cells(&loaded, field);
    assert_eq!(header, 1, "{ordinary}: one `## {field}` header line");
    assert!(
        cells > 0,
        "{ordinary}: the field reaches the rows ({cells} of {rows})"
    );
    assert!(
        !warnings(&loaded).contains("Failed to instantiate"),
        "{ordinary} loads: {}",
        warnings(&loaded)
    );

    let dropped = dir.join("outside.txt");
    let output = run(&corpus, &dropped, &outside.replace("{corpus}", &data));
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let warned = warnings(&dropped);
    assert!(
        warned.contains(&format!("WARNING: Failed to instantiate plugin {refusal}")),
        "{outside}: the refusal names the word: {warned}"
    );
    assert_eq!(
        header_and_cells(&dropped, field),
        (0, 0, rows),
        "{outside}: no header line, no value, the same {rows} rows"
    );
}

#[test]
fn loeuf_match_by_outside_gene_and_transcript_is_refused() {
    refused_where_ensembl_vep_loads(
        "GRCh38-plugin-loeuf",
        "LOEUF",
        "LOEUF,file={corpus}/data/loeuf_dataset_grch38.tsv.gz,match_by=gene",
        "LOEUF,file={corpus}/data/loeuf_dataset_grch38.tsv.gz,match_by=exon",
        "LOEUF: LOEUF: match_by must be gene or transcript, not 'exon'",
    );
}

#[test]
fn pli_mode_outside_gene_and_transcript_is_refused() {
    refused_where_ensembl_vep_loads(
        "GRCh38-plugin-pli",
        "pLI_gene_value",
        "pLI,{corpus}/plugin_data/pLI_values.txt,gene",
        "pLI,{corpus}/plugin_data/pLI_values.txt,exon",
        "pLI: pLI: mode 'exon' is not gene or transcript",
    );
}

#[test]
fn conservation_word_after_max_is_refused() {
    refused_where_ensembl_vep_loads(
        "GRCh38-plugin-conservation",
        "Conservation",
        "Conservation,{corpus}/plugin_data/gerp.bw,MAX",
        "Conservation,{corpus}/plugin_data/gerp.bw,MAX,foo",
        "Conservation: Conservation: MAX must be the last parameter",
    );
}
