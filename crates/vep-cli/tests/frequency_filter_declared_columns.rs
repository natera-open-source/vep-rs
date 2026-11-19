// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `--check_frequency` filter on a cache whose `info.json` declares
//! `variation_cols` while its variation files are absent, held to what Ensembl
//! VEP 116.2 writes on the same cache.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const INCLUDE: &[&str] = &[
    "--check_frequency",
    "--freq_pop",
    "1KG_ALL",
    "--freq_freq",
    "0.01",
    "--freq_gt_lt",
    "gt",
    "--freq_filter",
    "include",
];
const EXCLUDE: &[&str] = &[
    "--check_frequency",
    "--freq_pop",
    "1KG_ALL",
    "--freq_freq",
    "0.01",
    "--freq_gt_lt",
    "gt",
    "--freq_filter",
    "exclude",
];
const FILTER_COMMON: &[&str] = &["--filter_common"];

/// A committed corpus whose cache declares `variation_cols`, with the records
/// and rows a `--check_existing` run writes for its `variants.vcf`.
struct Corpus {
    name: &'static str,
    species: &'static str,
    assembly: &'static str,
    fasta: Option<&'static str>,
    records: usize,
    rows: usize,
}

const HUMAN: Corpus = Corpus {
    name: "GRCh38-freqfilter",
    species: "homo_sapiens",
    assembly: "GRCh38",
    fasta: None,
    records: 15,
    rows: 291,
};
const YEAST: Corpus = Corpus {
    name: "R64-1-1-yeast",
    species: "saccharomyces_cerevisiae",
    assembly: "R64-1-1",
    fasta: Some("reference.fa.gz"),
    records: 263,
    rows: 1733,
};

fn corpus_dir(corpus: &Corpus) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(corpus.name)
}

/// Copies the directory tree at `from` to `to`, leaving out the top-level
/// entries named in `skip`.
fn copy_tree(from: &Path, to: &Path, skip: &[&str]) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        if skip.iter().any(|name| entry.file_name() == *name) {
            continue;
        }
        let target = to.join(entry.file_name());
        if entry.path().is_dir() {
            copy_tree(&entry.path(), &target, &[]);
        } else {
            std::fs::copy(entry.path(), &target).unwrap();
        }
    }
}

/// The corpus's cache copied to `dir` without its variation files and without
/// the checksum manifest that lists them; `info.json` keeps `variation_cols`.
fn cache_without_variation_files(dir: &Path, corpus: &Corpus) -> PathBuf {
    let cache = dir.join(format!("{}-declared", corpus.name));
    copy_tree(
        &corpus_dir(corpus).join("json_cache"),
        &cache,
        &["variations", "CHECKSUMS.sha256"],
    );
    assert!(!cache.join("variations").exists());
    let info = read_info(&cache);
    assert!(
        info["variation_cols"]
            .as_array()
            .is_some_and(|c| !c.is_empty()),
        "{} declares variation_cols",
        corpus.name
    );
    cache
}

/// [`cache_without_variation_files`] with `variation_cols` and `var_type`
/// taken out of `info.json`.
fn cache_without_declaration(dir: &Path, corpus: &Corpus) -> PathBuf {
    let declared = cache_without_variation_files(dir, corpus);
    let cache = dir.join(format!("{}-undeclared", corpus.name));
    copy_tree(&declared, &cache, &[]);
    let mut info = read_info(&cache);
    let fields = info.as_object_mut().unwrap();
    fields.remove("variation_cols");
    fields.remove("var_type");
    std::fs::write(cache.join("info.json"), info.to_string()).unwrap();
    cache
}

fn read_info(cache: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(cache.join("info.json")).unwrap()).unwrap()
}

/// Runs `vep` on `input`, a file of the corpus or an absolute path, against
/// `cache` with `flags` and `stats` (`--no_stats` or `--stats_text`) in the
/// default format; the run exits 0. Returns the output file.
fn run(
    dir: &Path,
    corpus: &Corpus,
    cache: &Path,
    input: &str,
    stats: &str,
    flags: &[&str],
) -> PathBuf {
    let out = dir.join("out.txt");
    let _ = std::fs::remove_file(&out);
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
    cmd.stdin(Stdio::null())
        .arg("-i")
        .arg(corpus_dir(corpus).join(input))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache)
        .args(["--offline", "--cache_version", "116"])
        .args(["--species", corpus.species, "--assembly", corpus.assembly])
        .args(["--format", "vcf", "--buffer_size", "5000", "--fork", "1"])
        .args([stats, "--quiet", "--force_overwrite"])
        .args(flags);
    if let Some(fasta) = corpus.fasta {
        cmd.arg("--fasta").arg(corpus_dir(corpus).join(fasta));
    }
    let output = cmd.output().unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{} {flags:?}: {}",
        corpus.name,
        String::from_utf8_lossy(&output.stderr)
    );
    out
}

/// The data rows of a run on the corpus's `input` against `cache` with `flags`.
fn rows_of(dir: &Path, corpus: &Corpus, cache: &Path, input: &str, flags: &[&str]) -> Vec<String> {
    let out = run(dir, corpus, cache, input, "--no_stats", flags);
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// [`rows_of`] the corpus's `variants.vcf`.
fn rows(dir: &Path, corpus: &Corpus, cache: &Path, flags: &[&str]) -> Vec<String> {
    rows_of(dir, corpus, cache, "variants.vcf", flags)
}

/// The `Variants filtered out` count in the text summary of a run on the
/// corpus's `input` against `cache` with `flags`.
fn filtered_out(dir: &Path, corpus: &Corpus, cache: &Path, input: &str, flags: &[&str]) -> String {
    let out = run(dir, corpus, cache, input, "--stats_text", flags);
    let summary = std::fs::read_to_string(format!("{}_summary.txt", out.display())).unwrap();
    summary
        .lines()
        .find_map(|line| line.strip_prefix("Variants filtered out\t"))
        .expect("the summary counts the filtered variants")
        .to_string()
}

/// The number of input records `rows` stand for, by their first column.
fn records(rows: &[String]) -> usize {
    let mut names: Vec<&str> = rows
        .iter()
        .map(|row| row.split('\t').next().unwrap())
        .collect();
    names.dedup();
    names.len()
}

/// The `--check_existing` rows of the corpus on `cache`, held to the corpus's
/// measured record and row counts.
fn control_rows(dir: &Path, corpus: &Corpus, cache: &Path) -> Vec<String> {
    let control = rows(dir, corpus, cache, &["--check_existing"]);
    assert_eq!(
        (records(&control), control.len()),
        (corpus.records, corpus.rows),
        "{}: the control run's records and rows",
        corpus.name
    );
    control
}

/// On a cache that declares `variation_cols` and holds no variation file the
/// filter runs, as Ensembl VEP builds its variation source from the declared
/// columns (CacheDir.pm:179) and reads no known variant from an absent file
/// (VariationTabix.pm:156 and 270): no record has frequency data, so `include` drops
/// every record (BaseCacheVariation.pm:319-322), each counted as filtered out
/// (:362-373), and `exclude`, which `--filter_common` sets, keeps every record
/// as the control writes it (:340-345).
fn assert_filter_runs_without_variation_files(corpus: &Corpus) {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_without_variation_files(dir.path(), corpus);
    let control = control_rows(dir.path(), corpus, &cache);
    let included = rows(dir.path(), corpus, &cache, INCLUDE);
    assert_eq!(
        (records(&included), included.len()),
        (0, 0),
        "{}: --freq_filter include keeps no record of {} without frequency data",
        corpus.name,
        corpus.records
    );
    assert_eq!(
        filtered_out(dir.path(), corpus, &cache, "variants.vcf", INCLUDE),
        corpus.records.to_string(),
        "{}: the summary counts every dropped record as filtered out",
        corpus.name
    );
    for flags in [EXCLUDE, FILTER_COMMON] {
        assert_eq!(
            rows(dir.path(), corpus, &cache, flags),
            control,
            "{} {flags:?}: every record is kept as the control writes it",
            corpus.name
        );
    }
}

/// Ensembl VEP 116.2 writes 0 of the 15 records of the GRCh38-freqfilter corpus
/// under `include`, its summary counting 15 filtered out, and all 15 (291 rows)
/// under `exclude` and `--filter_common`.
#[test]
fn frequency_filter_runs_on_declared_columns_without_variation_files_human() {
    assert_filter_runs_without_variation_files(&HUMAN);
}

/// Ensembl VEP 116.2 writes 0 of the 263 records of the yeast corpus under
/// `include`, its summary counting 263 filtered out, and all 263 (1,733 rows)
/// under `exclude` and `--filter_common`.
#[test]
fn frequency_filter_runs_on_declared_columns_without_variation_files_yeast() {
    assert_filter_runs_without_variation_files(&YEAST);
}

/// `flags` after `--minimal`.
fn minimal<'a>(flags: &[&'a str]) -> Vec<&'a str> {
    let mut all = vec!["--minimal"];
    all.extend_from_slice(flags);
    all
}

/// Under `--minimal` each alternate allele of a multi-allelic record is a
/// feature of its own, judged alone (BaseCacheVariation.pm:315), so the filter
/// judges a split record without a variation source too. On the
/// GRCh38-freqfilter corpus's `minimal.vcf` (15 records, two of them with two
/// alternate alleles) Ensembl VEP 116.2 writes 341 rows under
/// `--check_existing`, none under `include`, its summary counting 17 features
/// filtered out, and the 341 under `exclude` and `--filter_common`.
#[test]
fn frequency_filter_judges_split_records_without_variation_files() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_without_variation_files(dir.path(), &HUMAN);
    let rows_under =
        |flags: &[&str]| rows_of(dir.path(), &HUMAN, &cache, "minimal.vcf", &minimal(flags));
    let control = rows_under(&["--check_existing"]);
    assert_eq!((records(&control), control.len()), (15, 341));
    assert_eq!(
        rows_under(INCLUDE),
        Vec::<String>::new(),
        "--minimal --freq_filter include keeps no feature without frequency data"
    );
    assert_eq!(
        filtered_out(dir.path(), &HUMAN, &cache, "minimal.vcf", &minimal(INCLUDE)),
        "17",
        "the summary counts each of the 17 split features as filtered out"
    );
    for flags in [EXCLUDE, FILTER_COMMON] {
        assert_eq!(
            rows_under(flags),
            control,
            "--minimal {flags:?}: every record is kept as the control writes it"
        );
    }
}

/// A split record whose first alternate allele writes no row (`*` before an
/// allele that writes rows) has a first feature no allele stands for, judged
/// with the others; without a variation source none has frequency data. On
/// one such record Ensembl VEP 116.2 writes 16 rows under `--check_existing`,
/// none under `include`, its summary counting 2 features filtered out, and the
/// 16 under `exclude` and `--filter_common`.
#[test]
fn frequency_filter_judges_a_rowless_first_feature_without_variation_files() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache_without_variation_files(dir.path(), &HUMAN);
    let input = dir.path().join("star_first.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t31148973\t.\tGG\t*,TG\t.\t.\t.\n",
    )
    .unwrap();
    let input = input.to_str().unwrap();
    let rows_under = |flags: &[&str]| rows_of(dir.path(), &HUMAN, &cache, input, &minimal(flags));
    let control = rows_under(&["--check_existing"]);
    assert_eq!((records(&control), control.len()), (1, 16));
    assert_eq!(
        rows_under(INCLUDE),
        Vec::<String>::new(),
        "--minimal --freq_filter include keeps no feature without frequency data"
    );
    assert_eq!(
        filtered_out(dir.path(), &HUMAN, &cache, input, &minimal(INCLUDE)),
        "2",
        "the summary counts the rowless feature and the allele's as filtered out"
    );
    for flags in [EXCLUDE, FILTER_COMMON] {
        assert_eq!(
            rows_under(flags),
            control,
            "--minimal {flags:?}: the record is kept as the control writes it"
        );
    }
}

/// On a cache that declares no `variation_cols` Ensembl VEP builds no variation
/// source (CacheDir.pm:179), so no filter runs: `include`, `exclude` and
/// `--filter_common` each keep all 15 records (291 rows) of the
/// GRCh38-freqfilter corpus and all 263 (1,733 rows) of the yeast corpus, as
/// Ensembl VEP 116.2 writes them.
#[test]
fn frequency_filter_is_inactive_without_declared_columns() {
    let dir = tempfile::tempdir().unwrap();
    for corpus in [&HUMAN, &YEAST] {
        let cache = cache_without_declaration(dir.path(), corpus);
        let control = control_rows(dir.path(), corpus, &cache);
        for flags in [INCLUDE, EXCLUDE, FILTER_COMMON] {
            assert_eq!(
                rows(dir.path(), corpus, &cache, flags),
                control,
                "{} {flags:?}: the filter is inactive and every record is kept",
                corpus.name
            );
        }
    }
}
