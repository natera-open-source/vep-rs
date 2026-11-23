// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The transcript-set flags through the `vep` binary: `--refseq` and `--merged`
//! select the suffixed species directory of the cache root and the transcript
//! set the flag names, under `--json_cache` too, and `--all_refseq` keeps every
//! model. Each test quotes the documentation sentence it proves and cites the
//! Ensembl VEP source (release/116.3) that sentence describes.

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// The release 116 GRCh38 corpus cache: an Ensembl cache, every stable id an
/// `ENST` accession, no `source_cache` stamp and no `source_type` key.
fn corpus_cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

/// Records of the corpus on both of its chromosomes, each inside a transcript.
const RECORDS: &str = "\
21\t30282169\tstop_ins4_twostops_rev\tC\tCTATT\t.\t.\t.
21\t30348409\tstop_del2_lost_rev\tCTA\tC\t.\t.\t.
21\t30430756\tstop_ins5_newstop_fwd\tT\tTCCTAA\t.\t.\t.
21\t34716478\tpenult_ins3_stop_fwd\tA\tAATA\t.\t.\t.
22\t20061573\tpartial_codon_del3_fwd\tGAAG\tG\t.\t.\t.
";

fn write_input(dir: &Path) -> PathBuf {
    let path = dir.join("in.vcf");
    fs::write(
        &path,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{RECORDS}"),
    )
    .unwrap();
    path
}

/// Runs the binary on `in.vcf` under `dir`, writing `<dir>/<label>.default`,
/// with `cache` naming the cache and `extra` the flags under test. Nothing is
/// asserted about the exit status.
fn run(dir: &Path, label: &str, cache: &[&str], extra: &[&str]) -> Output {
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(dir.join("in.vcf"))
        .arg("-o")
        .arg(dir.join(format!("{label}.default")))
        .args(cache)
        .args(["--assembly", "GRCh38", "--force_overwrite"])
        .args(extra)
        .stdin(Stdio::null());
    cmd.output().unwrap()
}

/// [`run`] asserting success; the output text and the standard error.
fn annotate(dir: &Path, label: &str, cache: &[&str], extra: &[&str]) -> (String, String) {
    let output = run(dir, label, cache, extra);
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{extra:?}: {stderr}");
    let text = fs::read_to_string(dir.join(format!("{label}.default"))).unwrap();
    assert!(
        data_lines(&text).len() >= RECORDS.lines().count(),
        "{extra:?}: one row per record at least:\n{text}"
    );
    (text, stderr)
}

/// [`annotate`] on the corpus cache under `--json_cache`.
fn annotate_corpus(dir: &Path, label: &str, extra: &[&str]) -> (String, String) {
    let cache = corpus_cache();
    annotate(
        dir,
        label,
        &["--json_cache", cache.to_str().unwrap()],
        extra,
    )
}

/// The output without the header lines that differ between any two runs: the
/// command line, which names the flags under test, the time of the run and the
/// cache path.
fn without_volatile_headers(text: &str) -> Vec<&str> {
    const VOLATILE: [&str; 3] = [
        "## VEP command-line:",
        "## Output produced at ",
        "## Using cache in ",
    ];
    text.lines()
        .filter(|line| !VOLATILE.iter().any(|prefix| line.starts_with(prefix)))
        .collect()
}

fn data_lines(text: &str) -> Vec<&str> {
    text.lines().filter(|line| !line.starts_with('#')).collect()
}

/// The header lines of the columns `--refseq` and `--merged` add
/// (`REFSEQ_MATCH`, `SOURCE`, `REFSEQ_OFFSET`; Constants.pm:97-98), which no
/// row of an Ensembl cache fills.
const SET_COLUMN_HEADERS: [&str; 3] = ["## REFSEQ_MATCH : ", "## SOURCE : ", "## REFSEQ_OFFSET : "];

/// [`without_volatile_headers`] less the set columns' header lines, which the
/// run under the flag carries and the plain run does not, so the two compare
/// on the rows and every other header; asserts the flagged run lists
/// `REFSEQ_MATCH`.
fn without_set_columns(flagged: &str) -> Vec<&str> {
    assert!(
        flagged
            .lines()
            .any(|line| line.starts_with(SET_COLUMN_HEADERS[0])),
        "the flagged run lists REFSEQ_MATCH:\n{flagged}"
    );
    without_volatile_headers(flagged)
        .into_iter()
        .filter(|line| {
            !SET_COLUMN_HEADERS
                .iter()
                .any(|prefix| line.starts_with(prefix))
        })
        .collect()
}

/// Asserts every data row of `text` is the row of a record no transcript
/// overlaps: `Feature` `-` and `Consequence` `intergenic_variant`, one row per
/// record.
fn assert_every_row_intergenic(text: &str, what: &str) {
    let rows = data_lines(text);
    assert_eq!(
        rows.len(),
        RECORDS.lines().count(),
        "{what}: one row per record:\n{text}"
    );
    for row in rows {
        let columns: Vec<&str> = row.split('\t').collect();
        assert_eq!(columns[4], "-", "{what}: Feature of {row}");
        assert_eq!(
            columns[6], "intergenic_variant",
            "{what}: Consequence of {row}"
        );
    }
}

/// A cache root holding the corpus cache as `<root>/<species_dir>/116_GRCh38`
/// alone.
fn cache_root(dir: &Path, species_dir: &str) -> PathBuf {
    let root = dir.join("cache_root");
    copy_dir(&corpus_cache(), &root.join(species_dir).join("116_GRCh38"));
    root
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// The cache flags naming `root` in Ensembl VEP's layout.
fn layout_flags(root: &Path) -> [&str; 7] {
    [
        "--cache",
        "--dir_cache",
        root.to_str().unwrap(),
        "--cache_version",
        "116",
        "--species",
        "homo_sapiens",
    ]
}

/// Asserts `flag` resolves the cache under `<root>/<species_dir>/116_GRCh38`
/// and that the plain species directory, which the root lacks, is the cache
/// error without the flag; returns the suffixed run's output and the plain
/// `--json_cache` run's for the caller's comparison.
fn assert_species_directory(dir: &Path, flag: &str, species_dir: &str) -> (String, String) {
    write_input(dir);
    let root = cache_root(dir, species_dir);
    let layout = layout_flags(&root);
    let (text, _) = annotate(dir, "suffixed", &layout, &[flag]);
    let leaf = root.join(species_dir).join("116_GRCh38");
    assert!(
        text.lines()
            .any(|line| line == format!("## Using cache in {}", leaf.display())),
        "{flag}: the header names the suffixed directory:\n{text}"
    );

    let output = run(dir, "unsuffixed", &layout, &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "without {flag}: {stderr}");
    assert!(
        stderr.lines().any(|line| {
            line == format!(
                "ERROR: [cache] Cache directory {} not found",
                root.join("homo_sapiens").display()
            )
        }),
        "without {flag}: {stderr}"
    );
    let (reference, _) = annotate_corpus(dir, "plain", &[]);
    (text, reference)
}

/// `--refseq` selects the `_refseq` species directory of the cache root and
/// annotates against the RefSeq set: of the cache it finds there, the RefSeq
/// accessions alone, so the corpus cache's `ENST` models off `MT` are all
/// dropped and every record, none on `MT`, is intergenic, and `--all_refseq`
/// keeps them.
///
/// Documented: "The cache is read under the species directory suffixed
/// _refseq, or under `--json_cache` as named, and of its models only the RefSeq
/// accessions are kept" (Filtering, picking and transcript sets, Transcript
/// sets). Ensembl VEP appends `_refseq` to the species directory name under the
/// flag (`CacheDir.pm:312`), dies when the directory is absent (`:318`), reads
/// the set from the flags (`:517`) and drops every model whose stable id is not
/// a RefSeq accession unless `all_refseq` is set
/// (`AnnotationType/Transcript.pm:236-248`).
#[test]
fn refseq_selects_the_refseq_species_directory() {
    let dir = tempfile::tempdir().unwrap();
    let (text, reference) = assert_species_directory(dir.path(), "--refseq", "homo_sapiens_refseq");
    assert_every_row_intergenic(&text, "--refseq on an Ensembl cache");
    assert_ne!(
        without_volatile_headers(&text),
        without_volatile_headers(&reference)
    );

    let root = dir.path().join("cache_root");
    let (every_model, _) = annotate(
        dir.path(),
        "all",
        &layout_flags(&root),
        &["--refseq", "--all_refseq"],
    );
    assert_eq!(
        without_set_columns(&every_model),
        without_volatile_headers(&reference),
        "--refseq --all_refseq annotates every model of the cache"
    );
}

/// `--merged` selects the `_merged` species directory of the cache root
/// (`CacheDir.pm:312`, `:318`) and annotates against the merged set, whose
/// RefSeq clause runs on the models stamped `RefSeq` alone
/// (`AnnotationType/Transcript.pm:240`): the corpus cache carries no stamp, so
/// every model is kept and the run annotates as the corpus cache does.
///
/// Documented: "The cache is read under the species directory suffixed
/// _merged, or under `--json_cache` as named; every Ensembl model is kept"
/// (Filtering, picking and transcript sets, Transcript sets).
#[test]
fn merged_selects_the_merged_species_directory() {
    let dir = tempfile::tempdir().unwrap();
    let (text, reference) = assert_species_directory(dir.path(), "--merged", "homo_sapiens_merged");
    assert_eq!(
        without_set_columns(&text),
        without_volatile_headers(&reference),
        "--merged keeps every unstamped model"
    );
}

/// The golden corpus `name` of release 116.
fn golden_corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The cache flags naming `root` in Ensembl VEP's layout for the mouse.
fn mouse_layout_flags(root: &Path) -> [&str; 9] {
    [
        "--cache",
        "--dir_cache",
        root.to_str().unwrap(),
        "--cache_version",
        "116",
        "--species",
        "mus_musculus",
        "--assembly",
        "GRCm39",
    ]
}

/// Runs the binary on a mouse corpus's own records and reference, writing
/// `<dir>/<label>.default`, with `cache` naming the cache and `flags` the
/// transcript set. Nothing is asserted about the exit status.
fn run_mouse(dir: &Path, label: &str, corpus: &Path, cache: &[&str], flags: &[&str]) -> Output {
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(corpus.join("variants.vcf"))
        .arg("-o")
        .arg(dir.join(format!("{label}.default")))
        .arg("--offline")
        .args(cache)
        .args(["--force_overwrite", "--no_stats"])
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
        .args(flags)
        .stdin(Stdio::null());
    cmd.output().unwrap()
}

/// Asserts the run of `flag` on `corpus` under the layout root `root` reads
/// the cache at `<root>/<species_dir>/116_GRCm39` and writes the corpus's
/// expected rows: the same Uploaded_variation, Location, Allele, Gene,
/// Feature and Feature_type on every row, and the whole row wherever the
/// manifest lists no documented divergence for its Location, Allele and
/// Feature.
fn assert_mouse_rows(dir: &Path, corpus: &Path, root: &Path, flag: &str, species_dir: &str) {
    let output = run_mouse(dir, "suffixed", corpus, &mouse_layout_flags(root), &[flag]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{flag}: {stderr}");
    let text = fs::read_to_string(dir.join("suffixed.default")).unwrap();
    let leaf = root.join(species_dir).join("116_GRCm39");
    assert!(
        text.lines()
            .any(|line| line == format!("## Using cache in {}", leaf.display())),
        "{flag}: the header names the suffixed directory:\n{text}"
    );

    let manifest: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(corpus.join("manifest.json")).unwrap()).unwrap();
    let documented: HashSet<(String, String, String)> = manifest["divergences"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|item| {
            let cell = |key: &str| item[key].as_str().unwrap().to_string();
            (cell("location"), cell("allele"), cell("feature"))
        })
        .collect();
    let expected = common::read_expected(corpus, "default.txt");
    let split = |text: &str| -> (Vec<String>, Vec<String>) {
        let mut keys = Vec::new();
        let mut whole = Vec::new();
        for row in data_lines(text) {
            let columns: Vec<&str> = row.split('\t').collect();
            keys.push(columns[..6].join("\t"));
            let key = (
                columns[1].to_string(),
                columns[2].to_string(),
                columns[4].to_string(),
            );
            if !documented.contains(&key) {
                whole.push(row.to_string());
            }
        }
        keys.sort();
        whole.sort();
        (keys, whole)
    };
    let (keys, whole) = split(&text);
    let (expected_keys, expected_whole) = split(&expected);
    assert!(!keys.is_empty(), "{flag}: the run writes rows");
    assert_eq!(
        keys, expected_keys,
        "{flag}: a row for every row of the reference"
    );
    assert_eq!(
        whole, expected_whole,
        "{flag}: the rows outside the documented divergences"
    );
}

/// `--refseq` selects the `_refseq` species directory of a species other than
/// human: under `--species mus_musculus` the run reads
/// `cache_root/mus_musculus_refseq/116_GRCm39` of the mouse RefSeq corpus and
/// writes that corpus's expected rows, and without the flag the plain species
/// directory, which the root lacks, is the cache error.
///
/// Documented: "Under it the species directory is `--species` (default
/// `homo_sapiens`, the Latin name only), suffixed `_refseq` or `_merged` under
/// those flags" (Cache setup, Naming the cache). Ensembl VEP appends the suffix
/// to whatever species the run names (`CacheDir.pm:312`) and dies when the
/// directory is absent (`:318`).
#[test]
fn refseq_selects_the_refseq_directory_of_another_species() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = golden_corpus("GRCm39-mouse-refseq");
    let root = corpus.join("cache_root");
    assert_mouse_rows(
        dir.path(),
        &corpus,
        &root,
        "--refseq",
        "mus_musculus_refseq",
    );

    let output = run_mouse(
        dir.path(),
        "unsuffixed",
        &corpus,
        &mouse_layout_flags(&root),
        &[],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "without --refseq: {stderr}");
    assert!(
        stderr.lines().any(|line| {
            line == format!(
                "ERROR: [cache] Cache directory {} not found",
                root.join("mus_musculus").display()
            )
        }),
        "without --refseq: {stderr}"
    );
}

/// `--merged` selects the `_merged` species directory of a species other than
/// human (`CacheDir.pm:312`): the mouse merged corpus's cache, placed as
/// `<root>/mus_musculus_merged/116_GRCm39`, is the cache the run reads under
/// `--species mus_musculus --merged`, and the run writes that corpus's
/// expected rows.
///
/// Documented: "Under it the species directory is `--species` (default
/// `homo_sapiens`, the Latin name only), suffixed `_refseq` or `_merged` under
/// those flags" (Cache setup, Naming the cache).
#[test]
fn merged_selects_the_merged_directory_of_another_species() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = golden_corpus("GRCm39-mouse-merged");
    let root = dir.path().join("cache_root");
    copy_dir(
        &corpus.join("json_cache"),
        &root.join("mus_musculus_merged").join("116_GRCm39"),
    );
    assert_mouse_rows(
        dir.path(),
        &corpus,
        &root,
        "--merged",
        "mus_musculus_merged",
    );
}

/// Under `--json_cache` the flags name the set of the cache given, as
/// `CacheDir.pm:517` reads the set from the flags whatever directory
/// `--full_cache_dir` names (`AnnotationSourceAdaptor.pm:128`): `--refseq` on
/// the Ensembl corpus cache, which declares no `source_type`, drops every
/// model off `MT` and annotates the records, none on `MT`, with intergenic
/// rows alone; `--merged` and `--refseq --all_refseq` keep them all; and the
/// same cache declaring `source_type` `ensembl` is refused under `--refseq`
/// before any input is read.
///
/// Documented: "Under `--json_cache` the flags name the set of the cache given,
/// so `--refseq` on a converted Ensembl cache without the `source_type` key
/// keeps only the transcripts that pass its test of a RefSeq id, as Ensembl
/// VEP does" and "A run whose `--refseq` or `--merged` disagrees with the
/// declared set is refused before any input is read, where Ensembl VEP filters
/// whatever models the cache holds by the flags alone" (Cache conversion,
/// RefSeq and merged caches).
#[test]
fn refseq_under_json_cache_applies_the_set() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let (reference, _) = annotate_corpus(dir.path(), "plain", &[]);
    let (refseq, stderr) = annotate_corpus(dir.path(), "refseq", &["--refseq"]);
    assert_every_row_intergenic(&refseq, "--refseq --json_cache");
    assert!(!stderr.contains("ERROR"), "{stderr}");
    for flags in [&["--merged"][..], &["--refseq", "--all_refseq"][..]] {
        let (text, _) = annotate_corpus(dir.path(), "kept", flags);
        assert_eq!(
            without_set_columns(&text),
            without_volatile_headers(&reference),
            "{flags:?} keeps every model of the corpus cache"
        );
    }

    let declared = dir.path().join("declared");
    copy_dir(&corpus_cache(), &declared);
    let info = declared.join("info.json");
    let mut metadata: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(&info).unwrap()).unwrap();
    metadata["source_type"] = serde_json::Value::String("ensembl".to_string());
    fs::write(&info, serde_json::to_string_pretty(&metadata).unwrap()).unwrap();
    vep_core::cache_manifest::CacheManifest::write(&declared).unwrap();
    let output = run(
        dir.path(),
        "declared",
        &["--json_cache", declared.to_str().unwrap()],
        &["--refseq"],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.lines().any(|line| {
            line == format!(
                "ERROR: [cache] --refseq names the RefSeq transcript set but {} declares source_type ensembl",
                info.display()
            )
        }),
        "{stderr}"
    );
    assert!(!dir.path().join("declared.default").exists());
    let output = run(
        dir.path(),
        "declared-plain",
        &["--json_cache", declared.to_str().unwrap()],
        &[],
    );
    assert!(
        output.status.success(),
        "the declared Ensembl set loads without the flag"
    );
}
