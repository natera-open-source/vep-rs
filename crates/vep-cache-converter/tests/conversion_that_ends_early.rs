// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a cache conversion leaves in `--output-dir` when it ends early. The
//! conversion writes its checksum manifest twice: before its first file as
//! the one line that lists no file, and last as the list. A conversion that
//! fails between the two leaves the first, whichever files it had written, a
//! directory `vep` refuses by name; converting into the same directory
//! completes the cache.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use vep_core::cache_manifest::{CacheManifest, ManifestError, MANIFEST_FILE};

/// Writes a flat cache of one empty region file on chromosome 21 into `src`.
fn write_source(src: &Path) {
    std::fs::create_dir_all(src.join("transcripts")).unwrap();
    std::fs::write(
        src.join("info.json"),
        r#"{"species": "homo_sapiens", "assembly": "GRCh38", "cache_version": 116}"#,
    )
    .unwrap();
    std::fs::write(src.join("transcripts/21_1-1000000.json"), "[]").unwrap();
}

/// Adds a regulatory tree to the flat cache in `src` whose metadata names no
/// `regulatory_schema`: the conversion refuses the tree, after it has written
/// the region files.
fn add_refused_regulatory_tree(src: &Path) {
    std::fs::create_dir_all(src.join("regulatory/21")).unwrap();
    std::fs::write(src.join("regulatory/21/1-1000000.json"), "[]").unwrap();
}

fn convert(src: &Path, out: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .arg("--json-dir")
        .arg(src)
        .arg("--output-dir")
        .arg(out)
        .output()
        .expect("vep-cache-converter runs")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// Every entry under `dir` that is not a directory, by its path relative to
/// it, `/` between components.
fn files(dir: &Path) -> BTreeMap<String, PathBuf> {
    fn collect(dir: &Path, root: &Path, files: &mut BTreeMap<String, PathBuf>) {
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                collect(&path, root, files);
            } else {
                let relative: Vec<String> = path
                    .strip_prefix(root)
                    .unwrap()
                    .components()
                    .map(|part| part.as_os_str().to_string_lossy().into_owned())
                    .collect();
                files.insert(relative.join("/"), path);
            }
        }
    }
    let mut files = BTreeMap::new();
    collect(dir, dir, &mut files);
    files
}

/// The line the manifest's loader and `vep` refuse `out` with.
fn refusal(out: &Path) -> String {
    format!(
        "{} holds a cache whose build did not finish: build it again",
        out.display()
    )
}

/// Asserts that `out` is a directory a conversion left early: loading its
/// manifest is refused by name, and `vep` ends on that refusal under the
/// `cache` component.
fn assert_refused_by_name(dir: &Path, out: &Path) {
    match CacheManifest::load(out) {
        Err(error @ ManifestError::Unfinished { .. }) => {
            assert_eq!(error.to_string(), refusal(out));
        }
        other => panic!("the manifest of {} loads as {other:?}", out.display()),
    }
    assert_eq!(
        vep_refusal(dir, out),
        format!("ERROR: [cache] {}", refusal(out))
    );
}

/// The line `vep` ends on when it is run against `cache`.
fn vep_refusal(dir: &Path, cache: &Path) -> String {
    let input = dir.join("input.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t1500\trecord\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.join("out.txt");
    let argv = [
        "vep".to_string(),
        "--offline".to_string(),
        "--json_cache".to_string(),
        cache.display().to_string(),
        "--no_stats".to_string(),
        "--quiet".to_string(),
        "--force_overwrite".to_string(),
        "-i".to_string(),
        input.display().to_string(),
        "-o".to_string(),
        out.display().to_string(),
    ];
    let config = vep_cli::Config::from_argv(argv).expect("the command line is accepted");
    let error = vep_cli::runner::Runner::new(config)
        .and_then(|runner| runner.run())
        .expect_err("vep refuses the directory");
    vep_cli::config::CliError::from_run(&error)
}

/// A regulatory tree the conversion refuses fails it after the region files
/// are written, exit 1. The output directory then holds the region file and a
/// manifest that lists no file, which `vep` refuses by name; a conversion of
/// the source without the refused tree into the same directory completes and
/// lists the directory.
#[test]
fn a_conversion_that_fails_after_its_region_files_leaves_no_cache_behind() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("src");
    let out = dir.path().join("cache");
    write_source(&src);
    add_refused_regulatory_tree(&src);

    let failed = convert(&src, &out);
    let err = stderr(&failed);
    assert_eq!(failed.status.code(), Some(1), "{err}");
    assert!(
        err.contains("declare no schema; this build reads schema 1"),
        "{err}"
    );
    assert_eq!(
        files(&out).keys().collect::<Vec<_>>(),
        [MANIFEST_FILE, "transcripts/21/1-1000000.json"],
        "{err}"
    );
    assert_refused_by_name(dir.path(), &out);

    std::fs::remove_dir_all(src.join("regulatory")).unwrap();
    let converted = convert(&src, &out);
    assert!(converted.status.success(), "{}", stderr(&converted));
    let manifest = CacheManifest::load(&out).unwrap();
    assert!(manifest.is_present());
    manifest.require_complete().unwrap();
    assert_eq!(
        files(&out).keys().collect::<Vec<_>>(),
        [MANIFEST_FILE, "info.json", "transcripts/21/1-1000000.json"]
    );
    assert_eq!(manifest.len(), 2);
}

/// A conversion whose last write, the list, fails has written every file of
/// the cache, `info.json` among them, and exits 1: the manifest it leaves is
/// still the line that lists no file, so `vep` refuses the directory by name
/// and never reads it as a cache without a manifest. The list cannot be
/// written here because the output directory holds a link to a path that does
/// not exist.
#[cfg(unix)]
#[test]
fn a_conversion_whose_list_cannot_be_written_leaves_no_cache_behind() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("src");
    let out = dir.path().join("cache");
    write_source(&src);
    std::fs::create_dir_all(&out).unwrap();
    std::os::unix::fs::symlink(dir.path().join("absent"), out.join("link")).unwrap();

    let failed = convert(&src, &out);
    let err = stderr(&failed);
    assert_eq!(failed.status.code(), Some(1), "{err}");
    assert!(err.contains("Failed to write CHECKSUMS.sha256"), "{err}");
    assert_eq!(
        files(&out).keys().collect::<Vec<_>>(),
        [
            MANIFEST_FILE,
            "info.json",
            "link",
            "transcripts/21/1-1000000.json"
        ],
        "{err}"
    );
    assert_refused_by_name(dir.path(), &out);

    std::fs::remove_file(out.join("link")).unwrap();
    let converted = convert(&src, &out);
    assert!(converted.status.success(), "{}", stderr(&converted));
    let manifest = CacheManifest::load(&out).unwrap();
    manifest.require_complete().unwrap();
    assert_eq!(manifest.len(), 2);
}

/// A conversion into a directory that holds a finished conversion replaces
/// that cache's list by the line that lists no file before it writes: when it
/// then fails, `vep` refuses the directory by name, not under the list of the
/// cache it held.
#[test]
fn a_failed_conversion_over_a_finished_one_leaves_no_list_behind() {
    let dir = tempfile::tempdir().unwrap();
    let src = dir.path().join("src");
    let out = dir.path().join("cache");
    write_source(&src);
    let converted = convert(&src, &out);
    assert!(converted.status.success(), "{}", stderr(&converted));
    assert_eq!(CacheManifest::load(&out).unwrap().len(), 2);

    add_refused_regulatory_tree(&src);
    let failed = convert(&src, &out);
    let err = stderr(&failed);
    assert_eq!(failed.status.code(), Some(1), "{err}");
    assert_eq!(
        files(&out).keys().collect::<Vec<_>>(),
        [MANIFEST_FILE, "transcripts/21/1-1000000.json"],
        "{err}"
    );
    assert_refused_by_name(dir.path(), &out);
}
