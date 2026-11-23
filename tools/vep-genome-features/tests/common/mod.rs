// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests share: the committed species cache and copies of it a test rewrites,
//! transcript records built by hand, scratch directories and a runner of the built
//! binary.

#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use flate2::write::GzEncoder;
use flate2::Compression;
use serde_json::{json, Value};
use vep_genome_features::extract::{extract, Request};
use vep_genome_features::features::Features;

/// The committed cache of a made-up species: sequence `1` with one protein-coding
/// transcript of two exons on the standard code and one lncRNA of one exon, both on
/// the forward strand; one variation file with its columns named in `info.json`; no
/// synonyms file, no side file, no SIFT marker and no regulatory flag.
pub fn species_cache() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../crates/vep-cli/tests/fixtures/species/json_cache")
}

/// A directory of this test's own, removed when the value is dropped.
pub struct Scratch(PathBuf);

impl Scratch {
    pub fn new(label: &str) -> Scratch {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let dir = std::env::temp_dir().join(format!(
            "vep-genome-features-{label}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&dir).expect("a scratch directory");
        Scratch(dir)
    }

    pub fn path(&self) -> &Path {
        &self.0
    }

    /// A copy of the committed cache under this directory.
    pub fn species_copy(&self) -> PathBuf {
        let copy = self.0.join("json_cache");
        copy_tree(&species_cache(), &copy);
        copy
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).expect("a directory");
    for entry in fs::read_dir(from).expect("a readable directory") {
        let entry = entry.expect("a directory entry").path();
        let target = to.join(entry.file_name().expect("a file name"));
        if entry.is_dir() {
            copy_tree(&entry, &target);
        } else {
            fs::copy(&entry, &target).expect("a copied file");
        }
    }
}

/// A transcript record without a coding sequence, its exons in the order given.
pub fn transcript(id: &str, strand: i8, biotype: &str, exons: &[(u64, u64)]) -> Value {
    let start = exons.iter().map(|exon| exon.0).min().expect("an exon");
    let end = exons.iter().map(|exon| exon.1).max().expect("an exon");
    json!({
        "stable_id": id,
        "start": start,
        "end": end,
        "strand": strand,
        "biotype": biotype,
        "exons": exons.iter().map(|(start, end)| json!({"start": start, "end": end})).collect::<Vec<_>>(),
        "attributes": [],
    })
}

/// `record` with a coding sequence over its whole spliced length: the genetic code it
/// names (or none), the cached coding sequence and the cached peptide (or none).
pub fn coding(mut record: Value, code: Option<u8>, sequence: &str, peptide: Option<&str>) -> Value {
    record["cdna_coding_start"] = json!(1);
    record["cdna_coding_end"] = json!(sequence.len());
    let mut cache = json!({
        "translateable_seq": sequence,
        "mapper": {"start_phase": -1},
    });
    if let Some(code) = code {
        cache["codon_table"] = json!(code);
    }
    if let Some(peptide) = peptide {
        cache["peptide"] = json!(peptide);
    }
    record["variation_effect_feature_cache"] = cache;
    record
}

/// `record` with the named attributes.
pub fn with_attributes(mut record: Value, codes: &[&str]) -> Value {
    record["attributes"] = codes.iter().map(|code| json!({"code": code})).collect();
    record
}

/// Writes `records` as the plain shard `file` of `sequence`.
pub fn write_shard(cache: &Path, sequence: &str, file: &str, records: &[Value]) {
    let dir = cache.join("transcripts").join(sequence);
    fs::create_dir_all(&dir).expect("a sequence directory");
    fs::write(
        dir.join(file),
        serde_json::to_vec(records).expect("a shard's bytes"),
    )
    .expect("a written shard");
}

/// Writes `shard`, any JSON value, as the gzip shard `file` of `sequence`.
pub fn write_gzip_shard(cache: &Path, sequence: &str, file: &str, shard: &Value) {
    let dir = cache.join("transcripts").join(sequence);
    fs::create_dir_all(&dir).expect("a sequence directory");
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder
        .write_all(&serde_json::to_vec(shard).expect("a shard's bytes"))
        .expect("a compressed shard");
    fs::write(dir.join(file), encoder.finish().expect("a gzip stream")).expect("a written shard");
}

/// Rewrites `info.json` through `edit`.
pub fn edit_info(cache: &Path, edit: impl FnOnce(&mut Value)) {
    let path = cache.join("info.json");
    let mut info: Value =
        serde_json::from_slice(&fs::read(&path).expect("info.json")).expect("a JSON object");
    edit(&mut info);
    fs::write(
        &path,
        serde_json::to_vec_pretty(&info).expect("info.json's bytes"),
    )
    .expect("a written info.json");
}

/// The features of the cache at `dir`, read with the side file beside `info.json`
/// where there is one.
pub fn features(dir: &Path) -> Features {
    extract(dir, &Request::default()).expect("the cache's features")
}

/// How a run of the built binary ended.
pub struct Ended {
    pub code: Option<i32>,
    pub stdout: String,
    pub stderr: String,
}

/// Runs the built binary.
pub fn run(args: &[&str]) -> Ended {
    let output = Command::new(env!("CARGO_BIN_EXE_vep-genome-features"))
        .args(args)
        .output()
        .expect("the built binary runs");
    Ended {
        code: output.status.code(),
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}
