// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Run control through the `vep` binary: help and version, the output file and
//! `STDOUT`, the cache resolved before anything runs, and the `runs` every
//! golden-corpus manifest declares.

mod common;

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use common::{corpora, execute_manifest_run, manifest_runs, Corpus};

fn vep() -> Command {
    Command::new(env!("CARGO_BIN_EXE_vep"))
}

fn corpus_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-release116")
}

/// The golden corpus's own arguments, output path left to the caller.
fn corpus_args() -> Vec<String> {
    let dir = corpus_dir();
    vec![
        "-i".to_string(),
        dir.join("variants.vcf").display().to_string(),
        "--offline".to_string(),
        "--json_cache".to_string(),
        dir.join("json_cache").display().to_string(),
        "--assembly".to_string(),
        "GRCh38".to_string(),
        "--force_overwrite".to_string(),
        "--no_stats".to_string(),
        "--max_sv_size".to_string(),
        "-1".to_string(),
    ]
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// Data lines of an output, the header dropped.
fn data_lines(output: &str) -> Vec<&str> {
    output.lines().filter(|l| !l.starts_with('#')).collect()
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
/// lines stay on standard error, and nothing is written to disk.
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
    assert_eq!(
        std::fs::read_dir(dir.path()).unwrap().count(),
        0,
        "nothing written to the working directory"
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
                    "stderr_contains"
                ]
            ),
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
