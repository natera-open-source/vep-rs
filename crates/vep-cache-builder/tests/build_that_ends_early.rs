// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a build leaves in `--output-dir`. The build writes its checksum
//! manifest twice: before its first file as the one line that lists no file,
//! and last as the list. A build that fails after its region files leaves the
//! first, a directory `vep` refuses by name and a build of the same directory
//! completes; a finished build leaves the list, every file of the directory
//! in it; and a directory that holds files keeps them, listed beside the
//! build's own.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use vep_core::cache_manifest::{CacheManifest, ManifestError, MANIFEST_FILE};

/// One gene of one single-exon transcript on chromosome 21.
const MODEL: &str = "##gff-version 3\n\
     21\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1;biotype=lncRNA\n\
     21\tx\tlnc_RNA\t1000\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=lncRNA\n\
     21\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript:T1\n";

/// A regulatory GFF the bake refuses at its first line: the `ID` decodes to
/// bytes that are not UTF-8.
const REFUSED_REGULATORY_GFF: &str = "21\tEnsembl\tenhancer\t1000\t1100\t.\t.\t.\tID=caf%FF\n";

/// Writes the model and a protein FASTA into `dir` and builds them into
/// `out`, with `regulatory` as `--regulatory-gff` when given.
fn build(dir: &Path, out: &Path, regulatory: Option<&Path>) -> Output {
    let model = dir.join("model.gff3");
    std::fs::write(&model, MODEL).unwrap();
    let peptides = dir.join("peptides.fa");
    std::fs::write(&peptides, ">ENSP1 pep transcript:T1 gene:G1\nMK\n").unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"));
    command
        .args(["--release", "116", "--assembly", "GRCh38", "--gff3"])
        .arg(&model)
        .arg("--fasta")
        .arg(&peptides)
        .arg("--output-dir")
        .arg(out);
    if let Some(regulatory) = regulatory {
        command.arg("--regulatory-gff").arg(regulatory);
    }
    command.output().expect("vep-cache-builder runs")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// Every file under `dir` by its path relative to it, `/` between components.
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

/// A regulatory GFF the bake refuses fails the build after the transcript
/// shards are written, exit 1. The output directory then holds the shards and
/// a manifest that lists no file: loading it is refused by name, `vep` ends on
/// that refusal under the `cache` component, and a build of the same
/// directory without the refused file completes and lists the directory.
#[test]
fn a_build_that_fails_after_the_region_files_leaves_no_cache_behind() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("cache");
    let regulatory = dir.path().join("refused.gff3");
    std::fs::write(&regulatory, REFUSED_REGULATORY_GFF).unwrap();

    let failed = build(dir.path(), &out, Some(&regulatory));
    let err = stderr(&failed);
    assert_eq!(failed.status.code(), Some(1), "{err}");
    let errors: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(errors.len(), 1, "{err}");
    assert!(
        errors[0].starts_with("ERROR: [cache-builder] Failed to build the regulatory tree: "),
        "{err}"
    );

    let left = files(&out);
    assert_eq!(
        left.keys().collect::<Vec<_>>(),
        [MANIFEST_FILE, "transcripts/21/1-1000000.json"],
        "{err}"
    );
    let refusal = format!(
        "{} holds a cache whose build did not finish: build it again",
        out.display()
    );
    let error = CacheManifest::load(&out).unwrap_err();
    assert!(matches!(error, ManifestError::Unfinished { .. }), "{error}");
    assert_eq!(error.to_string(), refusal);
    assert_eq!(
        vep_refusal(dir.path(), &out),
        format!("ERROR: [cache] {refusal}")
    );

    let rebuilt = build(dir.path(), &out, None);
    assert!(rebuilt.status.success(), "{}", stderr(&rebuilt));
    let manifest = CacheManifest::load(&out).unwrap();
    assert!(manifest.is_present());
    assert_eq!(manifest.len(), files(&out).len() - 1);
}

/// A finished build leaves the list in place of the mark: the manifest
/// loads, every file of the directory but the manifest is listed with its
/// digest, and nothing is written beside the directory.
#[test]
fn a_finished_build_leaves_the_list_of_its_files_and_no_mark() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("cache");
    let output = build(dir.path(), &out, None);
    assert!(output.status.success(), "{}", stderr(&output));

    let manifest = CacheManifest::load(&out).unwrap();
    assert!(manifest.is_present());
    manifest.require_complete().unwrap();
    let written = files(&out);
    assert_eq!(
        written.keys().collect::<Vec<_>>(),
        [MANIFEST_FILE, "info.json", "transcripts/21/1-1000000.json"]
    );
    assert_eq!(manifest.len(), written.len() - 1);
    for (relative, path) in &written {
        if relative != MANIFEST_FILE {
            manifest
                .verify_bytes(path, &std::fs::read(path).unwrap())
                .unwrap();
        }
    }
    let mut beside: Vec<String> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    beside.sort();
    assert_eq!(beside, ["cache", "model.gff3", "peptides.fa"]);
}

/// An output directory that holds files is built into as it stands: a file
/// of another name is kept byte for byte and listed beside the build's own,
/// and a manifest the directory held is replaced by the list of all of them.
#[test]
fn a_build_into_a_directory_that_holds_files_keeps_and_lists_them() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("cache");
    std::fs::create_dir_all(out.join("notes")).unwrap();
    std::fs::write(out.join("notes/kept.txt"), "kept\n").unwrap();
    assert_eq!(CacheManifest::write(&out).unwrap(), 1);

    let output = build(dir.path(), &out, None);
    assert!(output.status.success(), "{}", stderr(&output));

    assert_eq!(
        std::fs::read_to_string(out.join("notes/kept.txt")).unwrap(),
        "kept\n"
    );
    let manifest = CacheManifest::load(&out).unwrap();
    let written = files(&out);
    assert_eq!(
        written.keys().collect::<Vec<_>>(),
        [
            MANIFEST_FILE,
            "info.json",
            "notes/kept.txt",
            "transcripts/21/1-1000000.json"
        ]
    );
    assert_eq!(manifest.len(), 3);
    for relative in [
        "info.json",
        "notes/kept.txt",
        "transcripts/21/1-1000000.json",
    ] {
        manifest.require_listed(&out.join(relative)).unwrap();
    }
}
