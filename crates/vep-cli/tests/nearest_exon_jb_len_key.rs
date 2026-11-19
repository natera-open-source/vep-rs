// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--plugin NearestExonJB,max_range=len` through the `vep` binary: the plugin
//! loads, writes its header line, warns once that the value is not a number,
//! and reports no exon on any row. Ensembl VEP keys every exon's length under
//! `len` in the hash it tests the range value against (NearestExonJB.pm:229-237,
//! :268), so it reports every overlapped exon with the word as its distance,
//! the documented class `nearest_exon_jb_len_key_collision`.

use std::path::{Path, PathBuf};

mod common;

/// The corpus whose cache carries ENST00000336648; its record `1069237`
/// deletes a base 23 bases into ENSE00001282176 (272 bases).
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-plugin-computed")
}

const RECORD: &str = "21\t36896882\t1069237\tTG\tT";

/// Runs the binary on the record with one `--plugin NearestExonJB` value,
/// `--quiet`, so the warnings go to the file beside the output alone, and
/// returns the output's lines and the warnings file's text, which a run with
/// no warning never writes; each run writes its own output name.
fn annotate(dir: &Path, name: &str, plugin: &str) -> (Vec<String>, String) {
    let input = dir.join("in.vcf");
    std::fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{RECORD}\t.\t.\t.\n"),
    )
    .unwrap();
    let output = dir.join(format!("{name}.txt"));
    let warnings_file = dir.join(format!("{name}.txt_warnings.txt"));
    let corpus = corpus();
    let run = common::isolated_vep_anywhere()
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .arg("--dir_plugins")
        .arg(&corpus)
        .args(["--plugin", plugin])
        .output()
        .expect("vep binary runs");
    assert_eq!(
        run.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let lines = std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect();
    let warnings = std::fs::read_to_string(&warnings_file).unwrap_or_default();
    (lines, warnings)
}

/// The `NearestExonJB` item of the row on `feature`, if the row carries one.
fn item(lines: &[String], feature: &str) -> Option<String> {
    let row = lines
        .iter()
        .find(|line| line.split('\t').nth(4) == Some(feature))
        .expect("the transcript's row");
    row.split('\t')
        .nth(13)?
        .split(';')
        .find_map(|pair| pair.strip_prefix("NearestExonJB=").map(str::to_owned))
}

#[test]
fn max_range_len_loads_the_plugin_and_reports_no_exon() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();

    let (default, warnings) = annotate(dir, "default", "NearestExonJB");
    assert_eq!(
        item(&default, "ENST00000336648").as_deref(),
        Some("ENSE00001282176|23|start|272"),
        "the default range reports the overlapped exon"
    );
    assert!(!warnings.contains("NearestExonJB"), "{warnings}");

    let (len, warnings) = annotate(dir, "len", "NearestExonJB,max_range=len");
    assert_eq!(
        len.iter()
            .filter(|line| line.starts_with("## NearestExonJB : "))
            .count(),
        1,
        "the plugin loads and writes its header line"
    );
    let rows: Vec<&String> = len.iter().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(rows.len(), 17, "the record's rows are written");
    assert!(
        rows.iter().all(|row| !row.contains("NearestExonJB=")),
        "no row carries an item"
    );
    assert_eq!(
        warnings
            .lines()
            .filter(|line| line.contains("max_range='len' is not a number; it is read as 0"))
            .count(),
        1,
        "{warnings}"
    );
    assert!(!warnings.contains("Failed to instantiate"), "{warnings}");
}
