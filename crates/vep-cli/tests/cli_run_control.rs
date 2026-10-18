// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Run control through the `vep` binary: help and version, the output file and
//! `STDOUT`, the cache resolved before anything runs, the warnings file and the
//! skipped-variants file, `--quiet`, `--verbose`, `--safe`, the plugin loading
//! rules, and the `runs` every golden-corpus manifest declares.

mod common;

use std::collections::BTreeSet;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use common::{corpora, execute_manifest_run, manifest_runs, warning_entries, Corpus, WarningEntry};

fn vep() -> Command {
    Command::new(env!("CARGO_BIN_EXE_vep"))
}

fn corpus_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-release116")
}

/// The run-control corpus: the release116 records plus one on a chromosome the
/// cache lacks and one with the alleles `X/Y`, which Ensembl VEP skips.
fn run_control_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-run-control")
}

/// A corpus's own arguments, output path left to the caller; the cache is named as
/// the corpus lays it out (`Corpus::cache_args`).
fn args_for(dir: &Path) -> Vec<String> {
    let corpus = corpora()
        .into_iter()
        .find(|corpus| corpus.dir.canonicalize().ok() == dir.canonicalize().ok())
        .unwrap_or_else(|| panic!("{} is not a corpus under tests/golden", dir.display()));
    let mut args = vec![
        "-i".to_string(),
        dir.join("variants.vcf").display().to_string(),
        "--offline".to_string(),
    ];
    args.extend(corpus.cache_args());
    args.extend(
        [
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--max_sv_size",
            "-1",
        ]
        .map(str::to_string),
    );
    args
}

/// The golden corpus's own arguments, output path left to the caller.
fn corpus_args() -> Vec<String> {
    args_for(&corpus_dir())
}

/// The golden corpus's arguments with `input` read in place of its own records.
fn corpus_args_with_input(input: &Path) -> Vec<String> {
    let mut args = corpus_args();
    let input_flag = args.iter().position(|a| a == "-i").unwrap();
    args[input_flag + 1] = input.display().to_string();
    args
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// Data lines of an output, the header dropped.
fn data_lines(output: &str) -> Vec<&str> {
    output.lines().filter(|l| !l.starts_with('#')).collect()
}

/// Runs the binary on the run-control corpus with `extra` after the corpus
/// arguments, in `cwd`.
fn run_control(extra: &[&str], cwd: &Path) -> Output {
    vep()
        .args(args_for(&run_control_dir()))
        .args(extra)
        .current_dir(cwd)
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

/// The warnings Ensembl VEP writes for the run-control corpus, as committed
/// beside its expected outputs.
fn expected_warnings() -> String {
    std::fs::read_to_string(run_control_dir().join("expected/default.txt_warnings.txt")).unwrap()
}

fn warning_lines(stderr: &str) -> Vec<&str> {
    stderr
        .lines()
        .filter(|l| l.starts_with("WARNING: "))
        .collect()
}

/// A command line with no argument prints the help and exits 0, as the `vep`
/// script does (vep:29 tests the argument count alone), whether standard input
/// is a terminal, a closed pipe or an open one.
#[test]
fn empty_argv_prints_help_and_exits_0() {
    let closed = vep().stdin(Stdio::null()).output().unwrap();
    assert_eq!(closed.status.code(), Some(0));
    let stdout = text(&closed.stdout);
    assert!(stdout.starts_with("vep "), "{stdout}");
    assert!(
        stdout.contains("Predict functional effects of genomic variants"),
        "{stdout}"
    );

    let mut child = vep()
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    // The child never reads its input and may have exited before the write
    // lands, which leaves a broken pipe and nothing else.
    if let Err(e) = child
        .stdin
        .take()
        .unwrap()
        .write_all(b"21\t30165801\trs_ok\tA\tG\t.\t.\t.\n")
    {
        assert_eq!(e.kind(), std::io::ErrorKind::BrokenPipe, "{e}");
    }
    let piped = child.wait_with_output().unwrap();
    assert_eq!(piped.status.code(), Some(0));
    assert_eq!(
        text(&piped.stdout),
        stdout,
        "the help, not an annotation of the pipe"
    );
    let help = vep().arg("--help").output().unwrap();
    assert_eq!(help.status.code(), Some(0));
    assert_eq!(text(&help.stdout), stdout, "--help prints the same text");
}

/// The usage and the version win over a flag refused before parsing (vep:29
/// before vep:30): `--af_exac` beside `--help` or `--version`, in either order,
/// exits 0 with the help or the version.
#[test]
fn help_and_version_win_over_a_flag_refused_before_parsing() {
    for (argv, alone) in [
        (&["--af_exac", "--help"][..], "--help"),
        (&["--help", "--af_exac"][..], "--help"),
        (&["--af_exac=1", "-h"][..], "-h"),
        (&["--af_exac", "--version"][..], "--version"),
        (&["--version", "--af_exac"][..], "--version"),
        (&["-V", "--af_exac"][..], "-V"),
    ] {
        let expected = text(&vep().arg(alone).output().unwrap().stdout);
        let output = vep().args(argv).output().unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{argv:?}: {}",
            text(&output.stderr)
        );
        assert_eq!(text(&output.stdout), expected, "{argv:?}");
    }
}

/// `--version` prints the crate version and exits 0.
#[test]
fn version_prints_crate_version_and_exits_0() {
    let output = vep().arg("--version").output().unwrap();
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(
        text(&output.stdout),
        format!("vep {}\n", env!("CARGO_PKG_VERSION"))
    );
    assert!(output.stderr.is_empty());
}

fn run_to(output_args: &[&str], cwd: &Path) -> Output {
    vep()
        .args(corpus_args())
        .args(output_args)
        .current_dir(cwd)
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

/// `STDOUT` is matched in any case (BaseRunner.pm:213 `uc eq 'STDOUT'`): every
/// spelling writes the annotation to standard output and no file.
#[test]
fn stdout_literal_is_case_insensitive() {
    let dir = tempfile::tempdir().unwrap();
    let upper = run_to(&["-o", "STDOUT", "--quiet"], dir.path());
    assert_eq!(upper.status.code(), Some(0), "{}", text(&upper.stderr));
    let upper_rows = data_lines(&text(&upper.stdout)).len();
    assert!(upper_rows > 30, "{upper_rows} rows");
    for spelling in ["stdout", "StdOut", "Stdout"] {
        let output = run_to(&["-o", spelling, "--quiet"], dir.path());
        assert_eq!(
            output.status.code(),
            Some(0),
            "{spelling}: {}",
            text(&output.stderr)
        );
        assert_eq!(
            data_lines(&text(&output.stdout)),
            data_lines(&text(&upper.stdout)),
            "{spelling}"
        );
        assert!(
            !dir.path().join(spelling).exists(),
            "{spelling}: a file was written instead"
        );
    }
}

/// `-o` given twice keeps the last value and warns, as the `vep` script does
/// (vep:32-33: `WARNING: More than one output file command-line flag provided.
/// Using the last value provided`), under `--quiet` too, since the script warns
/// before any option is read; the first path is never written.
#[test]
fn repeated_output_file_keeps_last_with_warning() {
    let dir = tempfile::tempdir().unwrap();
    let first = dir.path().join("first.txt");
    let second = dir.path().join("second.txt");
    let warning = format!(
        "WARNING: More than one output file command-line flag provided. Using the last value provided: '{}'",
        second.display()
    );
    for quiet in [&[][..], &["--quiet"][..]] {
        let mut args = vec![
            "-o",
            first.to_str().unwrap(),
            "-o",
            second.to_str().unwrap(),
        ];
        args.extend_from_slice(quiet);
        let output = run_to(&args, dir.path());
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        assert!(second.is_file() && !first.exists());
        let stderr = text(&output.stderr);
        assert!(stderr.lines().any(|l| l == warning), "{quiet:?}: {stderr}");
        std::fs::remove_file(&second).unwrap();
    }
    let once = run_to(&["-o", first.to_str().unwrap()], dir.path());
    assert!(
        !text(&once.stderr).contains("More than one output file"),
        "a single -o is not warned about"
    );
}

/// `-o stdout` in any case streams the annotation to standard output, the status
/// lines stay on standard error, and the only file written is the warnings file
/// the corpus's skipped records call for, `<spelling>_warnings.txt` in the
/// working directory (BaseVEP.pm:890 names it after the output file as given).
#[test]
fn stdout_any_case_writes_stdout() {
    let dir = tempfile::tempdir().unwrap();
    let output = run_to(&["-o", "StdOut"], dir.path());
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let stdout = text(&output.stdout);
    assert!(
        stdout.starts_with("## ENSEMBL VARIANT EFFECT PREDICTOR"),
        "{stdout}"
    );
    assert!(data_lines(&stdout).len() > 30);
    let stderr = text(&output.stderr);
    assert!(
        stderr.contains("Ensembl VEP (Rust)"),
        "status lines on stderr: {stderr}"
    );
    assert!(
        !stdout.contains("Ensembl VEP (Rust)"),
        "no status line on stdout"
    );
    let written: Vec<String> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert_eq!(
        written,
        ["StdOut_warnings.txt"],
        "the warnings file alone is written to the working directory"
    );
}

/// An output file that exists is refused before the run with Ensembl VEP's text
/// (BaseRunner.pm:203-204) under the output component; `--force_overwrite` lifts it.
#[test]
fn existing_output_is_refused_with_ensembl_veps_text() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("taken.txt");
    std::fs::write(&out, "kept\n").unwrap();
    let mut args = corpus_args();
    args.retain(|a| a != "--force_overwrite");
    let refused = vep()
        .args(&args)
        .args(["-o", out.to_str().unwrap(), "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        format!(
            "ERROR: [output] Output file {} already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite",
            out.display()
        )
    );
    assert_eq!(std::fs::read_to_string(&out).unwrap(), "kept\n");
    let forced = run_to(&["-o", out.to_str().unwrap(), "--quiet"], dir.path());
    assert_eq!(forced.status.code(), Some(0), "{}", text(&forced.stderr));
    assert!(data_lines(&std::fs::read_to_string(&out).unwrap()).len() > 30);
}

/// A run whose cache resolves to nothing exits 1 under the cache component before
/// the input is read or the output created: a cache root with no species
/// directory (`CacheDir.pm:318`), a species directory with no cache of the
/// version (`:336`), and a `--json_cache` directory that is no cache.
#[test]
fn no_cache_is_fatal_before_first_variant() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let input = dir.path().join("variants.vcf");
    std::fs::write(&input, "21\t30165801\trs_ok\tA\tG\t.\t.\t.\n").unwrap();
    let missing_species = dir.path().join("empty-root");
    std::fs::create_dir_all(&missing_species).unwrap();
    let missing_version = dir.path().join("root");
    std::fs::create_dir_all(missing_version.join("homo_sapiens")).unwrap();
    let no_cache = dir.path().join("not-a-cache");
    std::fs::create_dir_all(&no_cache).unwrap();
    let cases = [
        (
            ["--dir_cache", missing_species.to_str().unwrap()],
            format!(
                "ERROR: [cache] Cache directory {}/homo_sapiens not found",
                missing_species.display()
            ),
        ),
        (
            ["--dir_cache", missing_version.to_str().unwrap()],
            "ERROR: [cache] No cache found for homo_sapiens, version 116".to_string(),
        ),
        (
            ["--json_cache", no_cache.to_str().unwrap()],
            format!(
                "ERROR: [cache] no cache found at {}: no transcripts directory",
                no_cache.display()
            ),
        ),
    ];
    for (cache_flags, expected) in cases {
        let output = vep()
            .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
            .args(["--offline", "--assembly", "GRCh38", "--quiet"])
            .args(cache_flags)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1), "{cache_flags:?}");
        assert_eq!(text(&output.stderr).trim_end(), expected, "{cache_flags:?}");
        assert!(output.stdout.is_empty(), "{cache_flags:?}");
        assert!(!out.exists(), "{cache_flags:?}: an output file was created");
    }
    // A run with no input file at all is refused for its cache first.
    let output = vep()
        .args(["-i", "/no/such/input.vcf", "--json_cache"])
        .arg(&no_cache)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(
        text(&output.stderr).contains("ERROR: [cache] no cache found at"),
        "{}",
        text(&output.stderr)
    );
}

/// `--show_cache_info` prints the cache's source versions as `<key>\t<value>`
/// lines sorted by key and exits 0 (vep:39-43 prints them in hash order), with
/// the cache resolved by the same flags as a run; a named input file that does
/// not exist is refused first (BaseRunner.pm:406, Parser.pm:239), and an
/// existing one is not read.
#[test]
fn show_cache_info_prints_info_json_and_exits_0() {
    let dir = corpus_dir();
    let info: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(dir.join("json_cache").join("info.json")).unwrap(),
    )
    .unwrap();
    let mut expected: Vec<String> = info["source_versions"]
        .as_object()
        .unwrap()
        .iter()
        .map(|(k, v)| format!("{k}\t{}", v.as_str().unwrap()))
        .collect();
    expected.sort();
    assert!(expected.len() > 5, "{expected:?}");
    let expected = expected.join("\n") + "\n";
    let output = vep()
        .args(corpus_args())
        .arg("--show_cache_info")
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(text(&output.stdout), expected);
    let mut args = corpus_args();
    let input_flag = args.iter().position(|a| a == "-i").unwrap();
    args[input_flag + 1] = "/no/such/input.vcf".to_string();
    let refused = vep()
        .args(&args)
        .args(["--show_cache_info", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [run] File \"/no/such/input.vcf\" does not exist"
    );
    assert!(refused.stdout.is_empty());
    // Under the Ensembl layout the same flags that select a cache print its info.
    let root = tempfile::tempdir().unwrap();
    let leaf = root.path().join("homo_sapiens").join("116_GRCh38");
    std::fs::create_dir_all(leaf.join("transcripts")).unwrap();
    std::fs::write(
        leaf.join("info.json"),
        r#"{"assembly":"GRCh38","source_versions":{"dbSNP":"156","assembly":"GRCh38.p14"}}"#,
    )
    .unwrap();
    let output = vep()
        .args(["--show_cache_info", "--cache", "--dir_cache"])
        .arg(root.path())
        .args(["--species", "homo_sapiens", "--cache_version", "116"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(text(&output.stdout), "assembly\tGRCh38.p14\ndbSNP\t156\n");
    let refused = vep()
        .args(["--show_cache_info", "--quiet", "--dir_cache"])
        .arg(root.path())
        .args(["--assembly", "GRCh37"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [cache] Cache assembly version (GRCh38) and database or selected assembly version (GRCh37) do not match"
    );
    assert!(refused.stdout.is_empty());
    // A refused --transcript_filter ends the run before the info prints: Ensembl
    // VEP compiles it while building the sources the info is read from
    // (BaseRunner.pm:415, AnnotationSource.pm:87-89).
    let refused = vep()
        .args(["--show_cache_info", "--quiet", "--cache", "--dir_cache"])
        .arg(root.path())
        .args(["--transcript_filter", "biotype equals protein_coding"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [filter] No such operator \"equals\""
    );
    assert!(refused.stdout.is_empty());
}

/// Parquet is a directory and cannot stream to STDOUT: the pair is refused
/// before anything runs.
#[test]
fn parquet_to_stdout_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let output = run_to(
        &["-o", "stdout", "--output_format", "parquet", "--quiet"],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        text(&output.stderr).trim_end(),
        "ERROR: [output] --output_format parquet needs a file path: Parquet cannot be written to STDOUT"
    );
    assert!(output.stdout.is_empty());
}

/// Every warning of a run is written to `<output_file>_warnings.txt` beside the
/// output (BaseVEP.pm:890), as `WARNING: line N skipped (<line>): <description>`
/// for a skipped record (`:840-863`), and echoed on standard error; the file is
/// Ensembl VEP's for the run-control corpus. A run without a warning leaves no
/// file (`:881`: the handle is opened at the first warning).
#[test]
fn warnings_file_is_written_by_default() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = run_control(&["-o", out.to_str().unwrap()], dir.path());
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    assert_eq!(warnings, expected_warnings());
    assert_eq!(
        warning_lines(&text(&output.stderr)),
        expected_warnings().lines().collect::<Vec<_>>(),
        "every warning is echoed on standard error"
    );

    let clean = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38");
    let mut args = args_for(&clean);
    args.extend([
        "-o".to_string(),
        dir.path().join("clean.txt").display().to_string(),
    ]);
    let output = vep().args(&args).stdin(Stdio::null()).output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(
        !dir.path().join("clean.txt_warnings.txt").exists(),
        "a run without a warning writes no warnings file"
    );
}

/// `--gencode_primary` or `--flag_gencode_primary` with an `--assembly` other
/// than GRCh38 is warned about in Ensembl VEP's words as the configuration
/// resolves (Config.pm:722-729), so the warning reaches the warnings file and
/// standard error before the cache is opened; here the GRCh38 cache then
/// refuses the GRCh37 assembly (CacheDir.pm:331-369) and the run exits 1.
#[test]
fn gencode_primary_assembly_warning_precedes_the_run() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let corpus = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-gencode-basic");
    let mut args = args_for(&corpus);
    let assembly = args.iter().position(|a| a == "--assembly").unwrap();
    args[assembly + 1] = "GRCh37".to_string();
    args.extend(["--gencode_primary", "-o", out.to_str().unwrap()].map(str::to_string));
    let output = vep().args(&args).stdin(Stdio::null()).output().unwrap();
    let warning =
        "WARNING: --gencode_primary option is currently only available for human on the GRCh38 assembly\n";
    assert_eq!(output.status.code(), Some(1), "{}", text(&output.stderr));
    assert_eq!(
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap(),
        warning
    );
    let stderr = text(&output.stderr);
    assert_eq!(warning_lines(&stderr), [warning.trim_end()]);
    assert!(
        stderr
            .lines()
            .last()
            .unwrap()
            .starts_with("ERROR: [cache] Mismatch in assembly versions"),
        "{stderr}"
    );
}

/// `--warning_file` names the file; the literal `stderr` in any case
/// (BaseVEP.pm:884 `/^stderr$/i`) keeps the warnings on standard error, once
/// each, and writes no file, `--quiet` or not.
#[test]
fn warning_file_stderr_literal() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let named = dir.path().join("named_warnings.txt");
    let output = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--warning_file",
            named.to_str().unwrap(),
        ],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        std::fs::read_to_string(&named).unwrap(),
        expected_warnings()
    );
    assert!(!dir.path().join("out.txt_warnings.txt").exists());

    for (spelling, quiet) in [("STDERR", false), ("stderr", true), ("StdErr", false)] {
        let mut extra = vec!["-o", out.to_str().unwrap(), "--warning_file", spelling];
        if quiet {
            extra.push("--quiet");
        }
        let output = run_control(&extra, dir.path());
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        assert_eq!(
            warning_lines(&text(&output.stderr)),
            expected_warnings().lines().collect::<Vec<_>>(),
            "{spelling} quiet={quiet}: the warnings go to standard error once each"
        );
        assert!(
            !dir.path().join(spelling).exists(),
            "{spelling}: a file was written"
        );
        assert!(!dir.path().join("out.txt_warnings.txt").exists());
    }
}

/// A warning seen before is written once (BaseVEP.pm:810-813): the same unknown
/// plugin named twice is one line of the warnings file. The dedup key sets a
/// trailing `line N` aside, which `warnings::tests` pins on the sink itself.
#[test]
fn warnings_deduplicate_ignoring_line_number() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let corpus = run_control_dir();
    let output = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir_plugins",
            corpus.to_str().unwrap(),
            "--plugin",
            "NoSuchPlugin",
            "--plugin",
            "NoSuchPlugin",
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    let plugin_lines: Vec<&str> = warnings
        .lines()
        .filter(|l| l.contains("Failed to load plugin NoSuchPlugin"))
        .collect();
    assert_eq!(plugin_lines.len(), 1, "{warnings}");
}

/// Under `--quiet` the unknown-plugin warning is in the warnings file and not on
/// standard error (BaseVEP.pm:817-821: the file is written, the `warn` copy is
/// not).
#[test]
fn quiet_keeps_the_warning_in_the_warnings_file() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let corpus = run_control_dir();
    let output = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir_plugins",
            corpus.to_str().unwrap(),
            "--plugin",
            "NoSuchPlugin",
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(output.stderr.is_empty(), "{}", text(&output.stderr));
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    assert!(
        warnings
            .lines()
            .any(|l| l.starts_with("WARNING: Failed to load plugin NoSuchPlugin: ")),
        "{warnings}"
    );
}

/// Under `--dont_skip` a record on a chromosome the cache does not carry is
/// kept (Parser/VCF.pm:199-201) with its warning still written, and, having no
/// feature, it is intergenic as Ensembl VEP writes a kept record with none.
#[test]
fn dont_skip_keeps_an_absent_chromosome_record_as_intergenic() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = run_control(
        &["-o", out.to_str().unwrap(), "--dont_skip", "--quiet"],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let written = std::fs::read_to_string(&out).unwrap();
    let row = data_lines(&written)
        .into_iter()
        .find(|l| l.starts_with("unknown_chromosome\t"))
        .expect("the kept record has a row");
    let columns: Vec<&str> = row.split('\t').collect();
    assert_eq!(columns[1], "Z:1000");
    assert_eq!(columns[6], "intergenic_variant");
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    let expected = "WARNING: line 40 skipped (Z 1000 unknown_chromosome A G . . .): Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features";
    assert!(warnings.lines().any(|l| l == expected), "{warnings}");
}

/// A skip in the ensembl, region and HGVS formats names the record's number as
/// Ensembl VEP counts it, one per record read (Parser/VEP_input.pm:119,
/// Parser/Region.pm:167, Parser/HGVS.pm:151), so two records skipped for one
/// reason are two warnings. The quoted line Ensembl VEP adds from its parser's
/// `current_block` (BaseVEP.pm:846-851) is absent: these parsers expose no line
/// text.
#[test]
fn non_vcf_skips_carry_the_record_number() {
    let dir = tempfile::tempdir().unwrap();
    for (format, records) in [
        (
            "ensembl",
            "21 30165801 30165801 A/G 1\n\
             Z 1000 1000 A/G 1\n\
             Z 2000 2000 A/G 1\n\
             21 30165811 30165811 A/G 1\n",
        ),
        (
            "region",
            "21:30165801-30165801:1/G\n\
             Z:1000-1000:1/G\n\
             Z:2000-2000:1/G\n\
             21:30165811-30165811:1/G\n",
        ),
        (
            "hgvs",
            "21:g.30165801A>G\n\
             Z:g.1000A>G\n\
             Z:g.2000A>G\n\
             21:g.30165811A>G\n",
        ),
    ] {
        let input = dir.path().join(format!("in.{format}"));
        std::fs::write(&input, records).unwrap();
        let out = dir.path().join(format!("out.{format}"));
        let output = vep()
            .args(corpus_args_with_input(&input))
            .args(["--format", format, "-o", out.to_str().unwrap(), "--quiet"])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{format}: {}",
            text(&output.stderr)
        );
        assert_eq!(
            std::fs::read_to_string(format!("{}_warnings.txt", out.display()))
                .unwrap()
                .lines()
                .collect::<Vec<_>>(),
            [
                "WARNING: line 2 skipped Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features",
                "WARNING: line 3 skipped Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features",
            ],
            "{format}"
        );
        let written = std::fs::read_to_string(&out).unwrap();
        let rows = data_lines(&written);
        assert!(
            !rows.is_empty(),
            "{format}: the records on chromosome 21 are annotated"
        );
        assert!(
            rows.iter().all(|l| !l.contains("\tZ:")),
            "{format}: a record on chromosome Z is written: {rows:?}"
        );
    }
}

/// A record on a chromosome the cache does not carry is skipped with
/// Parser.pm:542-546's text and gets no row, line or object in any format: no
/// intergenic verdict stands in for the annotation.
#[test]
fn chromosome_absent_is_warned_not_intergenic() {
    let dir = tempfile::tempdir().unwrap();
    let expected = "WARNING: line 40 skipped (Z 1000 unknown_chromosome A G . . .): Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features";
    for (format, flag) in [
        ("default", None),
        ("tab", Some("--tab")),
        ("vcf", Some("--vcf")),
        ("json", Some("--json")),
    ] {
        let out = dir.path().join(format!("out.{format}"));
        let mut extra = vec!["-o", out.to_str().unwrap(), "--quiet"];
        extra.extend(flag);
        let output = run_control(&extra, dir.path());
        assert_eq!(
            output.status.code(),
            Some(0),
            "{format}: {}",
            text(&output.stderr)
        );
        let written = std::fs::read_to_string(&out).unwrap();
        assert!(
            !written.contains("unknown_chromosome") && !written.contains("Z:1000"),
            "{format}: the record on chromosome Z is written:\n{written}"
        );
        let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
        assert!(
            warnings.lines().any(|l| l == expected),
            "{format}: {warnings}"
        );
    }
}

/// A line the parser refuses and a record `validate_vf` refuses are each one
/// `WARNING: line N skipped (<line>): <description>` on standard error and in the
/// warnings file, N the data record's ordinal (Parser/VCF.pm:228), the line
/// collapsed to single spaces; the records before and after them are annotated.
#[test]
fn skipped_line_writes_the_vep_warning_line_to_stderr_and_the_warnings_file() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30165801\tfirst\tA\tG\t.\t.\t.\n\
         21\t-1\tbad_pos\tA\tG\t.\t.\t.\n\
         21\t30165901\tinvalid_alleles\tX\tY\t.\t.\t.\n\
         21\t30165811\tlast\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.path().join("out.txt");
    let output = vep()
        .args(corpus_args_with_input(&input))
        .args(["-o", out.to_str().unwrap()])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let expected = [
        "WARNING: line 2 skipped (21 -1 bad_pos A G . . .): Invalid POS: '-1'",
        "WARNING: line 3 skipped (21 30165901 invalid_alleles X Y . . .): Invalid allele string X/Y or possible parsing error",
    ];
    assert_eq!(warning_lines(&text(&output.stderr)), expected);
    assert_eq!(
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt"))
            .unwrap()
            .lines()
            .collect::<Vec<_>>(),
        expected
    );
    let ids: BTreeSet<String> = data_lines(&std::fs::read_to_string(&out).unwrap())
        .iter()
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(
        ids,
        BTreeSet::from(["first".to_string(), "last".to_string()])
    );
}

/// The allele-string check reads the working string after the anchor chop and the
/// minimisation (Parser.pm:554-556 runs on the post-processed feature): an indel
/// whose shared leading base is not a nucleotide minimises to `X/-`, which carries
/// a `-` and is annotated, while `X/Y` is skipped.
#[test]
fn allele_string_check_reads_the_minimised_string() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30165801\tkept_indel\tXX\tX\t.\t.\t.\n\
         21\t30165901\tinvalid_alleles\tX\tY\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.path().join("out.txt");
    let output = vep()
        .args(corpus_args_with_input(&input))
        .args(["-o", out.to_str().unwrap(), "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt"))
            .unwrap()
            .lines()
            .collect::<Vec<_>>(),
        ["WARNING: line 2 skipped (21 30165901 invalid_alleles X Y . . .): Invalid allele string X/Y or possible parsing error"]
    );
    let written = std::fs::read_to_string(&out).unwrap();
    let rows = data_lines(&written);
    assert!(!rows.is_empty(), "the indel is annotated");
    assert!(
        rows.iter().all(|l| l.starts_with("kept_indel\t")),
        "a row other than the indel's is written: {rows:?}"
    );
}

/// A structural record's skips come in Ensembl VEP's order: the type `get_SO_term`
/// refuses (Parser.pm:725), then `deletion looks incomplete` for a type string
/// reading `/del/i`, the raw type when `%SO_TERMS` lacks it, over a span with no
/// extent (Parser/VCF.pm:469-481, :555-558), then the span `validate_vf` refuses
/// (Parser.pm:484-486).
#[test]
fn deletion_looks_incomplete_follows_the_type_and_precedes_the_span() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30165900\tdel\tN\t<DEL>\t.\t.\tEND=30165850\n\
         21\t30165910\todd\tN\t<FOO:DEL>\t.\t.\t.\n\
         21\t30165920\tseq\tN\tN\t.\t.\tSVTYPE=XDEL\n\
         21\t30165930\todd_end\tN\t<FOO:DEL>\t.\t.\tEND=30165800\n\
         21\t30165940\tinv\tN\t<INV>\t.\t.\tEND=30165800\n",
    )
    .unwrap();
    let out = dir.path().join("out.txt");
    let output = vep()
        .args(corpus_args_with_input(&input))
        .args(["-o", out.to_str().unwrap()])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let expected = [
        "WARNING: line 1 skipped (21 30165900 del N <DEL> . . END=30165850): deletion looks incomplete",
        "WARNING: line 1 skipped (21 30165900 del N <DEL> . . END=30165850): start > end+1 : (START=30165901, END=30165850)",
        "WARNING: line 2 skipped (21 30165910 odd N <FOO:DEL> . . .): FOO is not a supported structural variant type",
        "WARNING: line 2 skipped (21 30165910 odd N <FOO:DEL> . . .): deletion looks incomplete",
        "WARNING: line 3 skipped (21 30165920 seq N N . . SVTYPE=XDEL): XDEL is not a supported structural variant type",
        "WARNING: line 3 skipped (21 30165920 seq N N . . SVTYPE=XDEL): deletion looks incomplete",
        "WARNING: line 4 skipped (21 30165930 odd_end N <FOO:DEL> . . END=30165800): FOO is not a supported structural variant type",
        "WARNING: line 4 skipped (21 30165930 odd_end N <FOO:DEL> . . END=30165800): deletion looks incomplete",
        "WARNING: line 4 skipped (21 30165930 odd_end N <FOO:DEL> . . END=30165800): start > end+1 : (START=30165931, END=30165800)",
        "WARNING: line 5 skipped (21 30165940 inv N <INV> . . END=30165800): start > end+1 : (START=30165941, END=30165800)",
    ];
    assert_eq!(warning_lines(&text(&output.stderr)), expected);
}

/// `--skipped_variants_file` lists every skipped record in Ensembl VEP's layout
/// (BaseRunner.pm:364-380): a header naming the input, then
/// `Line %-5s\t%-40.40s\t%s` per skip; the rows are Ensembl VEP's for the
/// run-control corpus.
#[test]
fn skipped_variants_file_lists_every_skip_in_perl_layout() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let skipped = dir.path().join("skipped.txt");
    let output = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--skipped_variants_file",
            skipped.to_str().unwrap(),
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let written = std::fs::read_to_string(&skipped).unwrap();
    let expected = std::fs::read_to_string(run_control_dir().join("expected/skipped.txt")).unwrap();
    let mut lines = written.lines();
    assert_eq!(
        lines.next().unwrap(),
        format!(
            "[VEP skipped the following variants from {}]",
            run_control_dir().join("variants.vcf").display()
        )
    );
    assert_eq!(
        lines.collect::<Vec<_>>(),
        expected.lines().skip(1).collect::<Vec<_>>()
    );
}

/// `--quiet` drops the standard-error copy of each warning and nothing else: the
/// warnings file is still written in full (BaseVEP.pm:817-821).
#[test]
fn quiet_keeps_warnings_file() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = run_control(&["-o", out.to_str().unwrap(), "--quiet"], dir.path());
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(output.stderr.is_empty(), "{}", text(&output.stderr));
    assert_eq!(
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap(),
        expected_warnings()
    );
}

/// `--quiet` with `--verbose` is Ensembl VEP's incompatible pair (Config.pm:620,
/// :961-967), refused under the cli component; `--safe` skips the check (`:961`).
#[test]
fn quiet_and_verbose_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = run_control(
        &["-o", out.to_str().unwrap(), "--quiet", "--verbose"],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        text(&output.stderr).trim_end(),
        "ERROR: [cli] Can't use --quiet and --verbose together"
    );
    assert!(!out.exists());
}

/// `--safe` skips the incompatible-pair check, so `--quiet --verbose --safe` runs.
#[test]
fn safe_skips_incompatible_check() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--quiet",
            "--verbose",
            "--safe",
        ],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(data_lines(&std::fs::read_to_string(&out).unwrap()).len() > 30);
}

/// `--verbose` prints the status lines of the run, the detected input format
/// among them; without it the format is not reported.
#[test]
fn verbose_prints_the_detected_format() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let verbose = run_control(&["-o", out.to_str().unwrap(), "--verbose"], dir.path());
    assert_eq!(verbose.status.code(), Some(0), "{}", text(&verbose.stderr));
    assert!(
        text(&verbose.stderr).contains("detected input format: vcf"),
        "{}",
        text(&verbose.stderr)
    );
    let plain = run_control(&["-o", out.to_str().unwrap()], dir.path());
    assert!(!text(&plain.stderr).contains("detected input format"));
}

/// An unknown plugin name is a warnings-file line and the run goes on without it
/// (Runner.pm:929-939 warns on a compile failure); under `--safe` it is the
/// plugin component's error, before any output is written.
#[test]
fn safe_makes_unknown_plugin_fatal() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let corpus = run_control_dir();
    let warned = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir_plugins",
            corpus.to_str().unwrap(),
            "--plugin",
            "NoSuchPlugin",
        ],
        dir.path(),
    );
    assert_eq!(warned.status.code(), Some(0), "{}", text(&warned.stderr));
    let line = format!(
        "WARNING: Failed to load plugin NoSuchPlugin: NoSuchPlugin is not a built-in plugin and no libNoSuchPlugin.so or libNoSuchPlugin.dylib is under {}",
        corpus.display()
    );
    assert!(
        text(&warned.stderr).lines().any(|l| l == line),
        "{}",
        text(&warned.stderr)
    );
    assert!(
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt"))
            .unwrap()
            .lines()
            .any(|l| l == line)
    );
    assert!(data_lines(&std::fs::read_to_string(&out).unwrap()).len() > 30);
    std::fs::remove_file(&out).unwrap();

    let fatal = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir_plugins",
            corpus.to_str().unwrap(),
            "--plugin",
            "NoSuchPlugin",
            "--safe",
        ],
        dir.path(),
    );
    assert_eq!(fatal.status.code(), Some(1));
    let error = line.replacen("WARNING: ", "ERROR: [plugin] ", 1);
    assert!(
        text(&fatal.stderr).lines().any(|l| l == error),
        "{}",
        text(&fatal.stderr)
    );
    assert!(
        !out.exists(),
        "no output is written when --safe refuses a plugin"
    );
}

/// The plugins directory is `--dir_plugins`, else `<--dir>/Plugins`
/// (Runner.pm:907); one that does not exist is `Plugins directory '<dir>' not
/// found.` in the warnings file when a plugin is named (`:911-918`), nothing when
/// none is, and the plugin component's error under `--safe`.
#[test]
fn missing_plugins_dir_warns_only_when_a_plugin_is_named() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let vep_dir = dir.path().join("vepdata");
    std::fs::create_dir(&vep_dir).unwrap();
    let line = format!(
        "WARNING: Plugins directory '{}/Plugins' not found.",
        vep_dir.display()
    );
    let with_plugin = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir",
            vep_dir.to_str().unwrap(),
            "--plugin",
            "NoSuchPlugin",
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(
        with_plugin.status.code(),
        Some(0),
        "{}",
        text(&with_plugin.stderr)
    );
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    assert!(warnings.lines().any(|l| l == line), "{warnings}");
    assert!(
        warnings.contains(&format!(
            "no libNoSuchPlugin.so or libNoSuchPlugin.dylib is under {}/Plugins",
            vep_dir.display()
        )),
        "{warnings}"
    );

    let without_plugin = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir",
            vep_dir.to_str().unwrap(),
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(without_plugin.status.code(), Some(0));
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    assert!(!warnings.contains("Plugins directory"), "{warnings}");

    let safe = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir_plugins",
            "/nonexistent",
            "--plugin",
            "NoSuchPlugin",
            "--safe",
        ],
        dir.path(),
    );
    assert_eq!(safe.status.code(), Some(1));
    assert!(
        text(&safe.stderr)
            .lines()
            .any(|l| l == "ERROR: [plugin] Plugins directory '/nonexistent' not found."),
        "{}",
        text(&safe.stderr)
    );
}

/// An existing `--skipped_variants_file` is refused before the run with Ensembl
/// VEP's text (BaseRunner.pm:337-340, which Ensembl VEP checks at the end of the
/// run) under the output component; `--force_overwrite` lifts it.
#[test]
fn existing_skipped_file_is_refused_without_force() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let skipped = dir.path().join("skipped.txt");
    std::fs::write(&skipped, "kept\n").unwrap();
    let mut args = args_for(&run_control_dir());
    args.retain(|a| a != "--force_overwrite");
    let refused = vep()
        .args(&args)
        .args([
            "-o",
            out.to_str().unwrap(),
            "--skipped_variants_file",
            skipped.to_str().unwrap(),
            "--quiet",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        format!(
            "ERROR: [output] Skipped variants file {} already exists. Specify a different output file using --skipped_variants_file or overwrite existing file with --force_overwrite",
            skipped.display()
        )
    );
    assert_eq!(std::fs::read_to_string(&skipped).unwrap(), "kept\n");
    assert!(!out.exists(), "the refusal comes before the run");
    let forced = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--skipped_variants_file",
            skipped.to_str().unwrap(),
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(forced.status.code(), Some(0), "{}", text(&forced.stderr));
    assert!(std::fs::read_to_string(&skipped)
        .unwrap()
        .starts_with("[VEP skipped the following variants from "));
}

/// A warnings file that cannot be created fails the run with Ensembl VEP's text
/// (BaseVEP.pm:894) under the output component, exit 1; under `--quiet` the
/// error line is all that reaches standard error, no warning is echoed.
#[test]
fn unwritable_warnings_file_fails_the_run() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--warning_file",
            "/nonexistent-dir/w.txt",
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        text(&output.stderr).trim_end(),
        "ERROR: [output] Could not write to warnings file /nonexistent-dir/w.txt"
    );
}

/// `warning_entries` reads every line shape of a warnings file and a
/// skipped-variants file to one (line number, description) entry: a skip with
/// its number and one without (`variant skipped`, the `--input_data` shape;
/// BaseVEP.pm:854-855), a skip quoting no input line (`:847`), a
/// skipped-variants row (BaseRunner.pm:375), a warning whose trailing `line N`
/// (Parser.pm:284) or ` on line N` (Haplo/Parser/VCF.pm:220) is its number, a
/// plain warning, and a line of no known shape. The quoted input line is dropped,
/// a path token is masked, and the header line of a skipped-variants file, a
/// blank line and a repeated entry are left out.
#[test]
fn warning_entries_reads_every_line_shape() {
    let lines = "\
[VEP skipped the following variants from /data/variants.vcf]
WARNING: line 41 skipped (21 30165901 invalid_alleles X Y . . .): Invalid allele string X/Y or possible parsing error
WARNING: variant skipped (21 1 . X Y): Invalid allele string X/Y or possible parsing error
WARNING: line 8 skipped deletion looks incomplete
Line 13   \t21 30165800 multi_ins_del N <INS>,<DEL> \tINS/DEL is not a supported structural variant type
Line 40   \tZ 1000 unknown_chromosome A G . . .     \tChromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features
Line 41   \t21 30165901 invalid_alleles X Y . . .   \tInvalid allele string X/Y or possible parsing error
WARNING: Skipped 2 empty lines from line 7
WARNING: Chromosome X not found in annotation sources on line 9
WARNING: Failed to load plugin NoSuchPlugin: NoSuchPlugin is not a built-in plugin and no libNoSuchPlugin.so or libNoSuchPlugin.dylib is under /data/vep/Plugins
WARNING: text
Line 5 of the input is odd

some other line
";
    let expected: BTreeSet<WarningEntry> = [
        (Some(41), "Invalid allele string X/Y or possible parsing error"),
        (None, "Invalid allele string X/Y or possible parsing error"),
        (Some(8), "deletion looks incomplete"),
        (Some(13), "INS/DEL is not a supported structural variant type"),
        (
            Some(40),
            "Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features",
        ),
        (Some(7), "Skipped 2 empty lines from"),
        (Some(9), "Chromosome X not found in annotation sources"),
        (
            None,
            "Failed to load plugin NoSuchPlugin: NoSuchPlugin is not a built-in plugin and no libNoSuchPlugin.so or libNoSuchPlugin.dylib is under <path>",
        ),
        (None, "text"),
        (None, "Line 5 of the input is odd"),
        (None, "some other line"),
    ]
    .into_iter()
    .map(|(number, description)| (number, description.to_string()))
    .collect();
    assert_eq!(warning_entries(lines), expected);
}

/// Every `runs` entry of every corpus manifest behaves as its `expect` block says.
#[test]
fn manifest_runs_behave_as_expected() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    let mut executed = 0;
    for corpus in corpora() {
        for run in manifest_runs(&corpus) {
            executed += 1;
            let report = execute_manifest_run(&corpus, &run, tmp.path());
            if !report.is_empty() {
                failures.push(report);
            }
        }
    }
    assert!(executed > 0, "no manifest declares runs");
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

/// `manifest_runs` refuses the `runs` entries that would certify a manifest's
/// `covers` with a run that expects nothing, or that are not the shape the
/// helper documents, naming the entry and the reason: a name that is `default`,
/// repeated or not a plain word; an unknown key at either level; a value of the
/// wrong type; a `stdin` or `output` outside the corpus; an unknown format; and
/// an `expect` block with nothing beyond `exit`, an empty text counting as
/// nothing. The documented shape runs.
#[test]
fn manifest_runs_refuses_the_entries_that_expect_nothing_or_are_malformed() {
    let source = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-run-control")
        .expect("the GRCh38-run-control corpus is committed");
    let label = format!("{}/{}", source.release, source.name);
    let with_runs = |runs: serde_json::Value| {
        let mut manifest = source.manifest.clone();
        manifest["runs"] = runs;
        Corpus {
            release: source.release.clone(),
            name: source.name.clone(),
            assembly: source.assembly.clone(),
            dir: source.dir.clone(),
            manifest,
        }
    };
    let help = serde_json::json!({"exit": 0, "stdout_starts_with": "vep "});
    let named = |name: &str| serde_json::json!({"name": name, "argv": ["--help"], "expect": help});
    let expecting = |expect: serde_json::Value| serde_json::json!([{"name": "x", "argv": ["--help"], "expect": expect}]);
    let not_a_name =
        |name: &str| format!("{label}: run name {name:?} is not a plain word other than default");
    let refused = [
        (serde_json::json!([named("default")]), not_a_name("default")),
        (serde_json::json!([named("")]), not_a_name("")),
        (serde_json::json!([named("a/b")]), not_a_name("a/b")),
        (
            serde_json::json!([named("x"), named("x")]),
            format!("{label}: run x is named twice"),
        ),
        (
            serde_json::json!([{"name": "x", "argv": [], "expect": help, "cwd": "."}]),
            format!(
                "{label}: run key \"cwd\" is not one of {:?}",
                ["name", "argv", "env", "stdin", "expect"]
            ),
        ),
        (
            serde_json::json!([{"name": "x", "argv": "--help", "expect": help}]),
            format!("{label} run x: argv is a list of strings"),
        ),
        (
            serde_json::json!([{"name": "x", "argv": [], "env": "A=1", "expect": help}]),
            format!("{label} run x: env is an object"),
        ),
        (
            serde_json::json!([{"name": "x", "argv": [], "env": {"A": 1}, "expect": help}]),
            format!("{label} run x: env A is a string"),
        ),
        (
            serde_json::json!([{"name": "x", "argv": [], "stdin": "../variants.vcf", "expect": help}]),
            format!("{label} run x: stdin \"../variants.vcf\" is not a path inside the corpus"),
        ),
        (
            serde_json::json!([{"name": "x", "argv": []}]),
            format!("{label} run x: expect is an object"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "stdout": "vep "})),
            format!(
                "{label} run x: expect key \"stdout\" is not one of {:?}",
                [
                    "exit",
                    "output",
                    "format",
                    "stdout_starts_with",
                    "stdout_contains",
                    "stderr_contains",
                    "warnings_file",
                    "skipped_file"
                ]
            ),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "skipped_file": "skipped.txt"})),
            format!("{label} run x: skipped_file needs --skipped_variants_file in argv"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "warnings_file": "../w.txt"})),
            format!("{label} run x: warnings_file \"../w.txt\" is not a path inside the corpus"),
        ),
        (
            expecting(serde_json::json!({"stdout_starts_with": "vep "})),
            format!("{label} run x: expect.exit is an integer"),
        ),
        (
            expecting(serde_json::json!({"exit": "0", "stdout_starts_with": "vep "})),
            format!("{label} run x: expect.exit is an integer"),
        ),
        (
            expecting(serde_json::json!({"exit": 0})),
            format!("{label} run x: expects nothing beyond the exit status"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "output": null, "stderr_contains": null})),
            format!("{label} run x: expects nothing beyond the exit status"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "stderr_contains": ""})),
            format!("{label} run x: expect.stderr_contains is empty"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "stdout_starts_with": 1})),
            format!("{label} run x: expect.stdout_starts_with is a string"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "output": "../default.txt"})),
            format!("{label} run x: output \"../default.txt\" is not a path inside the corpus"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "output": "default.txt", "format": "parquet"})),
            format!("{label} run x: unknown format parquet"),
        ),
    ];
    let mut wrong = Vec::new();
    for (runs, reason) in refused {
        let corpus = with_runs(runs.clone());
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            manifest_runs(&corpus).len()
        }));
        let message = match &outcome {
            Ok(count) => format!("accepted, {count} runs"),
            Err(panic) => panic
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| panic.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_default(),
        };
        if message != reason {
            wrong.push(format!("{runs}: expected {reason:?}, got {message:?}"));
        }
    }
    assert!(wrong.is_empty(), "\n{}", wrong.join("\n"));
    let corpus = with_runs(serde_json::json!([
        {"name": "x1_-", "argv": ["{golden}", "-o", "{out}"], "expect": {"exit": 0, "output": "default.txt"}},
        {"name": "y", "argv": [], "env": {"A": "1"}, "stdin": "variants.vcf", "expect": {"exit": 0, "format": "tab", "stderr_contains": "x"}}
    ]));
    let accepted = manifest_runs(&corpus);
    assert_eq!(accepted.len(), 2);
    assert_eq!(accepted[0].expect.format, "default");
    assert_eq!(accepted[1].expect.format, "tab");
    assert_eq!(accepted[1].stdin.as_deref(), Some("variants.vcf"));
}

/// A manifest's `cache_layout` selects the cache directory and the arguments
/// naming it: absent, `json_cache/` through `--json_cache`; `ensembl`, the
/// committed leaf `cache_root/homo_sapiens/<release>_<assembly>` through
/// `--cache --dir_cache <root> --cache_version <release>`. Any other value,
/// `null` included, is refused with the corpus named rather than read as either.
#[test]
fn cache_layout_selects_the_cache_directory_and_refuses_another_value() {
    let source = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-run-control")
        .expect("the GRCh38-run-control corpus is committed");
    let with_layout = |layout: Option<serde_json::Value>| {
        let mut manifest = source.manifest.clone();
        match layout {
            Some(value) => manifest["cache_layout"] = value,
            None => {
                manifest.as_object_mut().unwrap().remove("cache_layout");
            }
        }
        Corpus {
            release: source.release.clone(),
            name: source.name.clone(),
            assembly: source.assembly.clone(),
            dir: source.dir.clone(),
            manifest,
        }
    };
    let root = source.dir.join("cache_root");
    let ensembl = with_layout(Some(serde_json::json!("ensembl")));
    assert!(ensembl.ensembl_cache_layout());
    assert_eq!(
        ensembl.cache_dir(),
        root.join("homo_sapiens").join("116_GRCh38")
    );
    assert!(ensembl.cache_dir().join("info.json").is_file());
    assert_eq!(
        ensembl.cache_args(),
        [
            "--cache",
            "--dir_cache",
            root.to_str().unwrap(),
            "--cache_version",
            "116"
        ]
    );
    let json_cache = source.dir.join("json_cache");
    let json = with_layout(None);
    assert!(!json.ensembl_cache_layout());
    assert_eq!(json.cache_dir(), json_cache);
    assert_eq!(
        json.cache_args(),
        ["--json_cache", json_cache.to_str().unwrap()]
    );
    for (value, shown) in [
        (serde_json::json!("json"), "\"json\""),
        (serde_json::json!("Ensembl"), "\"Ensembl\""),
        (serde_json::Value::Null, "null"),
    ] {
        let corpus = with_layout(Some(value));
        let outcome =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| corpus.cache_args()));
        let message = outcome
            .err()
            .and_then(|panic| panic.downcast_ref::<String>().cloned())
            .unwrap_or_default();
        assert_eq!(
            message,
            format!(
                "{}/{}: cache_layout is {shown}, not \"ensembl\"",
                source.release, source.name
            )
        );
    }
}
