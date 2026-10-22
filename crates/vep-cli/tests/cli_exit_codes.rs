// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The exit-status contract of the `vep` binary: 0 on success, 1 with one
//! `ERROR: [<component>] <message>` line on a refusal or a failed run, 2 with
//! clap's `error:` on an argument syntax error.

use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

mod common;

fn vep() -> Command {
    let mut cmd = common::isolated_vep_anywhere();
    cmd.stdin(Stdio::null());
    cmd
}

fn corpus_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-release116")
}

fn corpus_args(out: &Path) -> Vec<String> {
    let dir = corpus_dir();
    vec![
        "-i".to_string(),
        dir.join("variants.vcf").display().to_string(),
        "-o".to_string(),
        out.display().to_string(),
        "--offline".to_string(),
        "--json_cache".to_string(),
        dir.join("json_cache").display().to_string(),
        "--assembly".to_string(),
        "GRCh38".to_string(),
        "--force_overwrite".to_string(),
        "--no_stats".to_string(),
        "--max_sv_size".to_string(),
        "-1".to_string(),
        "--quiet".to_string(),
    ]
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

#[test]
fn success_is_0() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = vep().args(corpus_args(&out)).output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
    assert!(out.is_file());
    assert!(!stderr(&output).contains("ERROR:"));
}

/// A run that fails after parsing exits 1 with one contract line: here the
/// input file does not exist (Parser.pm:239).
#[test]
fn run_error_is_1() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let mut args = corpus_args(&out);
    args.extend(["-i".to_string(), "/nonexistent/input.vcf".to_string()]);
    let output = vep().args(args).output().unwrap();
    assert_eq!(output.status.code(), Some(1));
    let err = stderr(&output);
    let lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(lines.len(), 1, "{err}");
    assert_eq!(
        lines[0], "ERROR: [input] File \"/nonexistent/input.vcf\" does not exist",
        "{err}"
    );
    assert!(!out.exists(), "no output is written for a failed run");
}

/// An output file that cannot be created ends the run under the `output`
/// component with Ensembl VEP's text (BaseRunner.pm:226, `Could not write to
/// output file <file>`), before any variant is annotated; an empty output
/// name is such a file, named by the empty string.
#[test]
fn unwritable_output_is_1() {
    let dir = tempfile::tempdir().unwrap();
    for out in [
        dir.path().join("no_such_dir").join("out.txt"),
        PathBuf::new(),
    ] {
        let output = vep()
            .args(corpus_args(&out))
            .current_dir(dir.path())
            .output()
            .unwrap();
        let err = stderr(&output);
        assert_eq!(output.status.code(), Some(1), "{err}");
        let lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
        assert_eq!(
            lines,
            [format!(
                "ERROR: [output] Could not write to output file {}",
                out.display()
            )],
            "{err}"
        );
        assert!(output.stdout.is_empty());
    }
}

fn duckdb_available() -> bool {
    Command::new("duckdb")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// A Parquet run writes its rows to `<output>.tsv.tmp` beside the output file
/// first, so an output location that cannot be written fails the same way as
/// a plain output, naming that file, and a run whose output name is empty is
/// refused before any intermediate is written in the working directory; the
/// check needs the DuckDB CLI the format converts with.
#[test]
fn unwritable_parquet_output_is_1() {
    if !duckdb_available() {
        eprintln!("skipping: duckdb CLI not installed");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("no_such_dir").join("out.parquet");
    let mut args = corpus_args(&out);
    args.extend(["--output_format".to_string(), "parquet".to_string()]);
    let output = vep().args(args).output().unwrap();
    let err = stderr(&output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(
        lines,
        [format!(
            "ERROR: [output] Could not write to output file {}.tsv.tmp",
            out.display()
        )],
        "{err}"
    );
    assert!(output.stdout.is_empty());

    let mut args = corpus_args(Path::new(""));
    args.extend(["--output_format".to_string(), "parquet".to_string()]);
    let output = vep().args(args).current_dir(dir.path()).output().unwrap();
    let err = stderr(&output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(
        lines,
        ["ERROR: [output] Could not write to output file "],
        "{err}"
    );
    assert!(
        !dir.path().join(".tsv.tmp").exists(),
        "no intermediate is written for a refused output"
    );
}

/// A `--check_existing` run and the paths its assertions name: the output
/// file, the cache's `info.json` and its one region file.
struct VariationCacheRun {
    output: Output,
    out_file: PathBuf,
    info_json: PathBuf,
    region_file: PathBuf,
}

/// Runs `vep --check_existing` over one variant on chromosome 21 against a
/// cache under `dir` whose only variation file is
/// `variations/21/1-1000000_var.gz`, holding `region_body`, and whose
/// `info.json` declares `variation_cols` (a JSON array).
fn run_against_variation_cache(
    dir: &Path,
    region_body: &[u8],
    variation_cols: &str,
) -> VariationCacheRun {
    let cache = dir.join("json_cache");
    // The cache resolves only with a transcript shard beside the variation files.
    let shard = cache.join("transcripts").join("21").join("1-1000000.json");
    std::fs::create_dir_all(shard.parent().unwrap()).unwrap();
    std::fs::write(&shard, b"[]").unwrap();
    let region_file = cache.join("variations").join("21").join("1-1000000_var.gz");
    std::fs::create_dir_all(region_file.parent().unwrap()).unwrap();
    std::fs::write(&region_file, region_body).unwrap();
    let info_json = cache.join("info.json");
    std::fs::write(
        &info_json,
        format!(r#"{{"assembly":"GRCh38","variation_cols":{variation_cols}}}"#),
    )
    .unwrap();
    let input = dir.join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t500\t.\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    let out_file = dir.join("out.txt");
    let mut args = corpus_args(&out_file);
    args.extend([
        "-i".to_string(),
        input.display().to_string(),
        "--json_cache".to_string(),
        cache.display().to_string(),
        "--check_existing".to_string(),
    ]);
    let output = vep().args(args).output().unwrap();
    VariationCacheRun {
        output,
        out_file,
        info_json,
        region_file,
    }
}

/// A variation file of the cache that cannot be read ends the run under the
/// `cache` component before the batch that met it is written: no variant
/// goes out without its co-located variants.
#[test]
fn unreadable_variation_file_is_1() {
    let dir = tempfile::tempdir().unwrap();
    let run = run_against_variation_cache(
        dir.path(),
        b"not a gzip body",
        r#"["variation_name","failed","start","end","allele_string"]"#,
    );
    let err = stderr(&run.output);
    assert_eq!(run.output.status.code(), Some(1), "{err}");
    let lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(lines.len(), 1, "{err}");
    assert!(
        lines[0].starts_with(&format!(
            "ERROR: [cache] reading {}: ",
            run.region_file.display()
        )),
        "{err}"
    );
    let written = std::fs::read_to_string(&run.out_file).unwrap_or_default();
    assert!(
        written.lines().all(|line| line.starts_with('#')),
        "no variant row is written for the batch that met the failure:\n{written}"
    );
}

/// A variation layout the cache's `info.json` contradicts is refused when the
/// source opens, under the `cache` component, before any output is written.
#[test]
fn contradicted_variation_layout_is_1() {
    let dir = tempfile::tempdir().unwrap();
    let run = run_against_variation_cache(
        dir.path(),
        b"",
        r#"["variation_name","failed","end","allele_string"]"#,
    );
    let err = stderr(&run.output);
    assert_eq!(run.output.status.code(), Some(1), "{err}");
    let lines: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(
        lines,
        [format!(
            "ERROR: [cache] variation_cols of {} names no start column",
            run.info_json.display()
        )],
        "{err}"
    );
    assert!(
        !run.out_file.exists(),
        "no output is written for a run refused at open"
    );
}

/// An argument syntax error is clap's: `error:` prefix, exit 2.
#[test]
fn usage_error_is_2() {
    for argv in [
        vec!["--no_such_flag"],
        vec!["--fork"],
        vec!["--output_format", "xml"],
        vec!["--buffer_size", "many"],
    ] {
        let output = vep().args(&argv).output().unwrap();
        assert_eq!(output.status.code(), Some(2), "{argv:?}");
        let err = stderr(&output);
        assert!(err.starts_with("error: "), "{argv:?}: {err}");
        assert!(!err.contains("ERROR: ["), "{argv:?}: {err}");
    }
}

/// A flag the error table refuses exits 1 with its one line, before the input is
/// read; the clap usage path (exit 2) is not taken.
#[test]
fn refused_flag_is_1() {
    let output = vep()
        .args(["--database", "-i", "/nonexistent/input.vcf", "--quiet"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        stderr(&output),
        "ERROR: [cli] --database is not supported: vep-rs reads a local cache only and never connects to an Ensembl database; name the cache with --json_cache, or --cache with --dir_cache\n"
    );
    assert!(output.stdout.is_empty());
}
