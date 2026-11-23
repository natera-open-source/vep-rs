// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the `vep` binary prints where a run names a GFF or GTF transcript
//! source and the resolution of its cache is refused. The file is read beside
//! a cache, so where no cache stands at the place the run names, the refusal
//! is the source's: under `custom`, with the cache's own message and both
//! routes, naming the cache or building one from the file. A cache that stands
//! there and is refused, a `--species` or `--assembly` value the resolution
//! refuses and a `--synonyms` file that cannot be read end the run on the line
//! they draw without the source, under `cache` or `config`.

mod common;

use std::fs;
use std::path::Path;
use std::process::Stdio;

use common::{corpora, Corpus};
use vep_cli::config::GFF_GTF_RECIPE;

/// The corpus whose GTF is read beside its cache.
fn corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-gene-model-beside-cache")
        .expect("the gene-model corpus is committed")
}

/// The two spellings of the corpus's GTF as a transcript source.
fn sources(corpus: &Corpus) -> [Vec<String>; 2] {
    let file = corpus.dir.join("chr21_genes.gtf.gz").display().to_string();
    [
        vec!["--gtf".to_string(), file.clone()],
        vec!["--custom".to_string(), format!("file={file},format=gtf")],
    ]
}

/// A `vep` command over the corpus's records and FASTA, writing to `out.txt`
/// of `dir`, with `cache` naming its cache and `source` its transcript source.
fn command(
    dir: &Path,
    corpus: &Corpus,
    cache: &[&str],
    source: &[String],
) -> std::process::Command {
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(corpus.dir.join("variants.vcf"))
        .arg("-o")
        .arg(dir.join("out.txt"))
        .arg("--fasta")
        .arg(corpus.dir.join("reference.fa.gz"))
        .args(["--offline", "--no_stats", "--force_overwrite", "--quiet"])
        .args(cache)
        .args(source)
        .stdin(Stdio::null());
    cmd
}

/// The one `ERROR:` line of a run that `cache` and `source` end with exit code
/// 1 and no output file.
fn error_line(dir: &Path, corpus: &Corpus, cache: &[&str], source: &[String]) -> String {
    let output = command(dir, corpus, cache, source).output().unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(
        output.status.code(),
        Some(1),
        "{cache:?} {source:?}: {stderr}"
    );
    let errors: Vec<&str> = stderr.lines().filter(|l| l.starts_with("ERROR:")).collect();
    assert_eq!(errors.len(), 1, "{cache:?} {source:?}: {stderr}");
    assert!(
        !dir.join("out.txt").exists(),
        "{cache:?} {source:?}: a refused run wrote its output file"
    );
    errors[0].to_string()
}

/// Holds a run under `cache` to `line` without a transcript source and to the
/// same line with either spelling of one.
fn assert_the_line_stands_beside_a_source(dir: &Path, corpus: &Corpus, cache: &[&str], line: &str) {
    assert_eq!(error_line(dir, corpus, cache, &[]), line, "{cache:?}");
    for source in sources(corpus) {
        assert_eq!(
            error_line(dir, corpus, cache, &source),
            line,
            "{cache:?} {source:?}"
        );
    }
}

/// Holds a run under `cache` to `ERROR: [cache] <message>` without a
/// transcript source, and with either spelling of one to the source's refusal
/// under `custom`, which carries `message` and both routes.
fn assert_the_source_is_refused_for_want_of_a_cache(
    dir: &Path,
    corpus: &Corpus,
    cache: &[&str],
    message: &str,
) {
    assert_eq!(
        error_line(dir, corpus, cache, &[]),
        format!("ERROR: [cache] {message}"),
        "{cache:?}"
    );
    for source in sources(corpus) {
        assert_eq!(
            error_line(dir, corpus, cache, &source),
            format!(
                "ERROR: [custom] a GFF or GTF transcript source is read beside a cache, and \
                 this run has none ({message}); name the cache of the run with --json_cache \
                 <dir>, or with --cache and --dir_cache, or {GFF_GTF_RECIPE}"
            ),
            "{cache:?} {source:?}"
        );
    }
}

/// A directory under `dir` named `name`, with the directories `children` in it.
fn directory(dir: &Path, name: &str, children: &[&str]) -> String {
    let made = dir.join(name);
    fs::create_dir_all(&made).unwrap();
    for child in children {
        fs::create_dir_all(made.join(child)).unwrap();
    }
    made.display().to_string()
}

/// The corpus's cache is there and its assembly is not the one `--assembly`
/// names: the run ends on the cache's own line, Ensembl VEP's text with the
/// file renamed (`CacheDir.pm:407-408`), with and without the source.
#[test]
fn a_cache_whose_assembly_disagrees_is_refused_as_it_is_without_the_source() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();
    let cache = corpus.cache_dir().display().to_string();
    assert_the_line_stands_beside_a_source(
        dir.path(),
        &corpus,
        &["--json_cache", &cache, "--assembly", "GRCh37"],
        "ERROR: [cache] Mismatch in assembly versions from config (GRCh37) and cache info.json \
         file (GRCh38)",
    );
}

/// The corpus's cache is there and the `--synonyms` file is not: the run ends
/// on the file's own line under `config`, with and without the source.
#[test]
fn a_synonyms_file_that_cannot_be_read_is_refused_as_it_is_without_the_source() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();
    let cache = corpus.cache_dir().display().to_string();
    let synonyms = dir.path().join("absent-synonyms.txt").display().to_string();
    assert_the_line_stands_beside_a_source(
        dir.path(),
        &corpus,
        &[
            "--json_cache",
            &cache,
            "--assembly",
            "GRCh38",
            "--synonyms",
            &synonyms,
        ],
        &format!(
            "ERROR: [config] Could not read synonyms file {synonyms}: No such file or directory"
        ),
    );
}

/// A `--species` value Ensembl VEP refuses as a species directory
/// (`CacheDir.pm:310-311`) is refused under `cache` before any directory is
/// looked for, with and without the source.
#[test]
fn a_suffixed_species_is_refused_as_it_is_without_the_source() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();
    let root = directory(dir.path(), "root", &[]);
    assert_the_line_stands_beside_a_source(
        dir.path(),
        &corpus,
        &["--species", "homo_sapiens_refseq", "--dir_cache", &root],
        "ERROR: [cache] Should not use homo_sapiens_refseq as --species. Try using flags \
         --refseq or --merged with --species homo_sapiens",
    );
}

/// Under the cache root's layout, a cache directory whose name carries another
/// assembly than `--assembly` names, and several assemblies of one version
/// with no `--assembly`, are refused with Ensembl VEP's texts, the first
/// without its line on `--offline` (`CacheDir.pm:347-353`, `:361-364`), with
/// and without the source.
#[test]
fn a_layout_whose_assembly_the_flags_do_not_settle_is_refused_as_it_is_without_the_source() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();
    let one = directory(dir.path(), "one", &["homo_sapiens/116_GRCh38"]);
    assert_the_line_stands_beside_a_source(
        dir.path(),
        &corpus,
        &[
            "--dir_cache",
            &one,
            "--cache_version",
            "116",
            "--assembly",
            "GRCh37",
        ],
        "ERROR: [cache] Cache assembly version (GRCh38) and database or selected assembly \
         version (GRCh37) do not match",
    );
    let two = directory(
        dir.path(),
        "two",
        &["homo_sapiens/116_GRCh37", "homo_sapiens/116_GRCh38"],
    );
    assert_the_line_stands_beside_a_source(
        dir.path(),
        &corpus,
        &["--dir_cache", &two, "--cache_version", "116"],
        "ERROR: [cache] Multiple assemblies found for cache version 116 (GRCh37, GRCh38) - \
         specify one using --assembly [assembly]",
    );
}

/// A cache directory that is there and whose own files refuse it ends the run
/// on the line it draws without the source: a transcript set the flags do not
/// select, metadata that is no JSON object and the manifest of a build that
/// did not finish.
#[test]
fn a_cache_its_own_files_refuse_is_refused_as_it_is_without_the_source() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();

    let other_set = directory(dir.path(), "other-set", &["transcripts"]);
    fs::write(
        Path::new(&other_set).join("info.json"),
        r#"{"assembly":"GRCh38","source_type":"refseq"}"#,
    )
    .unwrap();
    assert_the_line_stands_beside_a_source(
        dir.path(),
        &corpus,
        &["--json_cache", &other_set],
        &format!(
            "ERROR: [cache] neither --refseq nor --merged was given, selecting the Ensembl \
             transcript set, but {other_set}/info.json declares source_type refseq"
        ),
    );

    let no_object = directory(dir.path(), "no-object", &["transcripts"]);
    fs::write(Path::new(&no_object).join("info.json"), "[]").unwrap();
    assert_the_line_stands_beside_a_source(
        dir.path(),
        &corpus,
        &["--json_cache", &no_object],
        &format!("ERROR: [cache] {no_object}/info.json: not a JSON object"),
    );

    let unfinished = directory(dir.path(), "unfinished", &[]);
    vep_core::cache_manifest::CacheManifest::mark_unfinished(&unfinished).unwrap();
    let without = error_line(dir.path(), &corpus, &["--json_cache", &unfinished], &[]);
    assert!(
        without.starts_with("ERROR: [cache] ") && without.contains(&unfinished),
        "{without}"
    );
    assert_the_line_stands_beside_a_source(
        dir.path(),
        &corpus,
        &["--json_cache", &unfinished],
        &without,
    );
}

/// A species path under the cache root that is a file cannot be listed, and
/// the run ends on the line that says so, with and without the source.
#[test]
fn a_species_path_that_cannot_be_listed_is_refused_as_it_is_without_the_source() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();
    let no_listing = directory(dir.path(), "no-listing", &[]);
    fs::write(Path::new(&no_listing).join("homo_sapiens"), "").unwrap();
    let without = error_line(dir.path(), &corpus, &["--dir_cache", &no_listing], &[]);
    assert!(
        without.starts_with(&format!("ERROR: [cache] {no_listing}/homo_sapiens: ")),
        "{without}"
    );
    assert_the_line_stands_beside_a_source(
        dir.path(),
        &corpus,
        &["--dir_cache", &no_listing],
        &without,
    );
}

/// No cache stands where the run names one: the source's refusal under
/// `custom` carries the cache's own message, Ensembl VEP's text where it has
/// one (`CacheDir.pm:109`, `:318`, `:336`, `:369`), and both routes.
#[test]
fn a_run_with_no_cache_where_it_names_one_is_told_both_routes() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();

    let absent = dir.path().join("absent").display().to_string();
    assert_the_source_is_refused_for_want_of_a_cache(
        dir.path(),
        &corpus,
        &["--json_cache", &absent],
        &format!("Directory {absent} does not exist"),
    );

    let empty_root = directory(dir.path(), "empty-root", &[]);
    assert_the_source_is_refused_for_want_of_a_cache(
        dir.path(),
        &corpus,
        &["--dir_cache", &empty_root],
        &format!("Cache directory {empty_root}/homo_sapiens not found"),
    );

    let no_version = directory(dir.path(), "no-version", &["homo_sapiens"]);
    assert_the_source_is_refused_for_want_of_a_cache(
        dir.path(),
        &corpus,
        &["--dir_cache", &no_version, "--cache_version", "116"],
        "No cache found for homo_sapiens, version 116",
    );

    let two = directory(
        dir.path(),
        "two",
        &["homo_sapiens/116_GRCh37", "homo_sapiens/116_GRCh38"],
    );
    assert_the_source_is_refused_for_want_of_a_cache(
        dir.path(),
        &corpus,
        &[
            "--dir_cache",
            &two,
            "--cache_version",
            "116",
            "--assembly",
            "CHM13",
        ],
        "No cache found for homo_sapiens, version 116, assembly CHM13",
    );

    let nothing = directory(dir.path(), "nothing", &[]);
    assert_the_source_is_refused_for_want_of_a_cache(
        dir.path(),
        &corpus,
        &["--json_cache", &nothing],
        &format!("no cache found at {nothing}: no transcripts directory"),
    );

    let no_metadata = directory(dir.path(), "no-metadata", &["transcripts"]);
    assert_the_source_is_refused_for_want_of_a_cache(
        dir.path(),
        &corpus,
        &["--json_cache", &no_metadata],
        &format!("no cache found at {no_metadata}: no info.json"),
    );
}

/// The cache flags alone reach these refusals: the same source beside the
/// corpus's cache, under the cache's own assembly, annotates, and rows of the
/// output carry the source's name.
#[test]
fn the_same_source_beside_the_corpus_cache_annotates() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();
    let cache = corpus.cache_dir().display().to_string();
    for source in sources(&corpus) {
        let output = command(
            dir.path(),
            &corpus,
            &["--json_cache", &cache, "--assembly", "GRCh38"],
            &source,
        )
        .output()
        .unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{source:?}: {stderr}");
        let written = fs::read_to_string(dir.path().join("out.txt")).unwrap();
        assert!(
            written
                .lines()
                .any(|l| !l.starts_with('#') && l.contains("SOURCE=chr21_genes.gtf.gz")),
            "{source:?}: no row carries the source's name"
        );
    }
}
