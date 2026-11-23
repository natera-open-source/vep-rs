// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Golden corpora round-tripped through the `vep` binary when one is built: a corpus is
//! annotated once into the `--tab` format and once into Parquet, and the canonical forms of
//! the two outputs are equal, a line a row, the tab output read under the types of the
//! Parquet output. One run is a flag set of its corpus: under `ALL` the synthetic dbNSFP
//! fixture declares columns named `ref` and `alt`, as two key columns of a part are named.
//!
//! The binary is `VEP_BINARY`, else the release or debug build under the cargo target
//! directory of the repository (`CARGO_TARGET_DIR`, else `target/`). Without one the tests
//! pass without running, and say so; `VEP_ROUND_TRIPS_MUST_RUN=1` fails them instead.

use std::env;
use std::ffi::OsString;
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

fn repository() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn vep_binary() -> Option<PathBuf> {
    if let Some(named) = env::var_os("VEP_BINARY").filter(|named| !named.is_empty()) {
        return Some(PathBuf::from(named));
    }
    let target =
        env::var_os("CARGO_TARGET_DIR").map_or_else(|| repository().join("target"), PathBuf::from);
    ["release", "debug"]
        .iter()
        .map(|profile| target.join(profile).join("vep"))
        .find(|candidate| {
            fs::metadata(candidate)
                .is_ok_and(|file| file.is_file() && file.permissions().mode() & 0o111 != 0)
        })
}

/// Ends a test that cannot run: quietly, or as a failure where the run is required.
fn skip(reason: &str) {
    assert!(
        env::var_os("VEP_ROUND_TRIPS_MUST_RUN").is_none_or(|value| value.is_empty()),
        "the round trip must run: {reason}"
    );
    eprintln!("skipped: {reason}");
}

/// The canonical form of the output at `path`; a `--tab` output is read under the types of
/// the Parquet output `types_of`.
fn canonical(path: &Path, types_of: Option<&Path>) -> String {
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-parquet-to-tab"));
    command.arg("--canonical");
    if let Some(parquet) = types_of {
        command.arg("--types-of").arg(parquet);
    }
    let output = command.arg(path).output().expect("the binary runs");
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}: {}",
        path.display(),
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("UTF-8 output")
}

/// The flags of a run of the corpus at `directory`: the manifest's own or, for a flag set,
/// that set's. Inside a flag `{corpus}` stands for the corpus directory and `{plugins}` for
/// its `plugin_data` directory.
fn flags_of(
    manifest: &serde_json::Value,
    flag_set: Option<&str>,
    directory: &Path,
) -> Vec<OsString> {
    let run = flag_set.map_or(manifest, |name| {
        manifest["flag_sets"]
            .as_array()
            .and_then(|sets| sets.iter().find(|set| set["name"] == name))
            .unwrap_or_else(|| panic!("the manifest lists no flag set {name}"))
    });
    let corpus = directory.display().to_string();
    let plugins = directory.join("plugin_data").display().to_string();
    run["flags"]
        .as_array()
        .map_or(&[][..], Vec::as_slice)
        .iter()
        .map(|flag| {
            let flag = flag.as_str().expect("a flag is a string");
            OsString::from(
                flag.replace("{corpus}", &corpus)
                    .replace("{plugins}", &plugins),
            )
        })
        .collect()
}

/// Annotates a corpus under its own flags, or under those of its flag set `flag_set`, into
/// the `--tab` format and into Parquet, and holds the two canonical forms equal, the column
/// line included. Returns the column line, `None` for a skipped run.
fn a_golden_corpus_round_trips_through_the_writer(
    corpus: &str,
    flag_set: Option<&str>,
) -> Option<String> {
    let Some(vep) = vep_binary() else {
        skip("no vep binary: set VEP_BINARY or build one under the cargo target directory");
        return None;
    };
    let directory = repository().join("tests/golden/116").join(corpus);
    let Ok(manifest) = fs::read_to_string(directory.join("manifest.json")) else {
        skip(&format!("no golden corpus at {}", directory.display()));
        return None;
    };
    let manifest: serde_json::Value = serde_json::from_str(&manifest).expect("a JSON manifest");
    let assembly = manifest["assembly"].as_str().expect("an assembly");
    let flags = flags_of(&manifest, flag_set, &directory);
    let scratch = env::temp_dir().join(format!(
        "vep-parquet-to-tab-round-trip-{corpus}-{}-{}",
        flag_set.unwrap_or("default"),
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&scratch);
    fs::create_dir_all(&scratch).expect("a scratch directory");
    let mut common: Vec<OsString> = vec![
        "--offline".into(),
        "--json_cache".into(),
        directory.join("json_cache").into(),
        "--species".into(),
        "homo_sapiens".into(),
        "--assembly".into(),
        assembly.into(),
        "--buffer_size".into(),
        "5000".into(),
        "--force_overwrite".into(),
        "--no_stats".into(),
        "--quiet".into(),
        "-i".into(),
        directory.join("variants.vcf").into(),
    ];
    common.extend(flags);
    let tab = scratch.join("out.tab");
    let parquet = scratch.join("out.parquet");
    let annotate = |format: &[&str], out: &Path| {
        let status = Command::new(&vep)
            .args(&common)
            .args(format)
            .arg("-o")
            .arg(out)
            .status()
            .expect("the vep binary runs");
        assert!(status.success(), "{format:?}: {status}");
    };
    annotate(&["--tab"], &tab);
    annotate(&["--output_format", "parquet"], &parquet);
    let from_parquet = canonical(&parquet, None);
    let from_tab = canonical(&tab, Some(&parquet));
    assert_eq!(
        from_parquet.lines().next(),
        from_tab.lines().next(),
        "{corpus}: the column lines differ"
    );
    assert!(
        from_parquet == from_tab,
        "{corpus}: the canonical forms differ"
    );
    let rows = fs::read_to_string(&tab)
        .expect("the tab output")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .count();
    assert!(rows > 0);
    assert_eq!(from_parquet.lines().count(), 1 + rows);
    fs::remove_dir_all(&scratch).expect("the scratch directory goes");
    from_parquet.lines().next().map(str::to_string)
}

#[test]
fn a_golden_corpus_round_trips_through_the_writer_grch38() {
    a_golden_corpus_round_trips_through_the_writer("GRCh38", None);
}

#[test]
fn a_golden_corpus_round_trips_through_the_writer_grch38_everything() {
    a_golden_corpus_round_trips_through_the_writer("GRCh38-everything", None);
}

#[test]
fn a_golden_corpus_round_trips_through_the_writer_grch38_minor_allele_columns() {
    a_golden_corpus_round_trips_through_the_writer("GRCh38-minor-allele-columns", None);
}

/// The synthetic dbNSFP fixture declares the columns `ref` and `alt` under `ALL` (every
/// column of the file's header, `dbNSFP.pm:233-236`), which the writer stores as `tab_ref`
/// and `tab_alt` beside the key columns of those names. Each prints the cells of its own
/// field under the name the `--tab` output's column line gives it, and the rows are the
/// `--tab` rows.
#[test]
fn a_run_whose_plugin_names_columns_as_the_key_columns_round_trips_through_the_writer() {
    let Some(column_line) = a_golden_corpus_round_trips_through_the_writer(
        "GRCh38-plugin-dbnsfp-fixture",
        Some("dbnsfp-all"),
    ) else {
        return;
    };
    let columns: Vec<&str> = column_line.trim_start_matches('#').split('\t').collect();
    for key in ["ref", "alt"] {
        assert_eq!(
            columns.iter().filter(|column| **column == key).count(),
            1,
            "{key}"
        );
        assert!(!columns.contains(&format!("tab_{key}").as_str()), "{key}");
    }
}
