// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A conversion of transcript files whose translations hold no `seq_edits`
//! list ends on one warning that says so and names `seq_edits` under
//! `left_out` in the `info.json` it writes; a conversion of files that hold
//! the list ends on its summary line and copies `info.json` byte for byte.

use std::fs;
use std::path::Path;
use std::process::Command;

/// The metadata of the made-up flat caches.
const INFO: &str = r#"{"species":"homo_sapiens","assembly":"GRCh38","cache_version":116}"#;

/// A transcript record with a translated peptide and no `seq_edits` list.
const WITHOUT_THE_LIST: &str =
    r#"{"stable_id":"T1","variation_effect_feature_cache":{"peptide":"MK"}}"#;

/// The same record with an empty list, as the Storable converter writes a
/// translation without edits.
const WITH_AN_EMPTY_LIST: &str =
    r#"{"stable_id":"T2","variation_effect_feature_cache":{"peptide":"MK","seq_edits":[]}}"#;

/// Writes at `cache` a flat cache of one region file of chromosome 21 that
/// holds `records`.
fn flat_cache(cache: &Path, records: &[&str]) {
    fs::create_dir_all(cache.join("transcripts")).unwrap();
    fs::write(cache.join("info.json"), INFO).unwrap();
    fs::write(
        cache.join("transcripts/21_1-1000000.json"),
        format!("[{}]", records.join(",")),
    )
    .unwrap();
}

/// Converts `cache` into `out` and returns the log, each line from its level
/// on.
fn conversion_log(cache: &Path, out: &Path) -> Vec<String> {
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .arg("--json-dir")
        .arg(cache)
        .arg("--output-dir")
        .arg(out)
        .output()
        .expect("vep-cache-converter runs");
    let log = String::from_utf8(output.stderr).unwrap();
    assert!(output.status.success(), "{log}");
    log.lines()
        .filter_map(|line| line.split_once(' '))
        .map(|(_, logged)| logged.trim_start().to_string())
        .collect()
}

#[test]
fn a_conversion_without_the_list_ends_on_the_line_and_records_the_name() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("flat");
    let out = dir.path().join("out");
    flat_cache(&cache, &[WITHOUT_THE_LIST, WITH_AN_EMPTY_LIST]);
    let log = conversion_log(&cache, &out);
    assert_eq!(
        log.last().map(String::as_str),
        Some(
            "WARN vep_cache_converter: No sequence edits (seq_edits) for 1 translation of the \
             cache: the transcript files hold no list of them, so a variant at codon 1 of a \
             transcript with an alternative initiation codon is called from the cached peptide; \
             a conversion of Ensembl VEP's cache with vep-storable-to-json carries them \
             (docs/cache-conversion.md)"
        ),
        "{log:?}"
    );
    let info: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out.join("info.json")).unwrap()).unwrap();
    assert_eq!(info["left_out"], serde_json::json!(["seq_edits"]));
    assert_eq!(info["species"], "homo_sapiens");
    assert_eq!(info["cache_version"], 116);
}

#[test]
fn a_conversion_with_the_list_ends_on_its_summary_and_copies_the_metadata() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("flat");
    let out = dir.path().join("out");
    flat_cache(&cache, &[WITH_AN_EMPTY_LIST]);
    let log = conversion_log(&cache, &out);
    let last = log.last().expect("the conversion logs its summary");
    assert!(
        last.starts_with("INFO vep_cache_converter: Conversion complete: "),
        "{log:?}"
    );
    assert!(
        !log.iter().any(|line| line.contains("seq_edits")),
        "{log:?}"
    );
    assert_eq!(fs::read_to_string(out.join("info.json")).unwrap(), INFO);
}
