// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Run control through the `vep` binary: help and version, the output file and
//! `STDOUT`, the cache resolved before anything runs, the FASTA and the HGVS
//! flags settled before the output opens, the warnings file and the
//! skipped-variants file, `--quiet`, `--verbose`, `--no_progress`, the banner,
//! the progress and summary lines, `--safe`, the plugin loading rules, and the
//! `runs` every golden-corpus manifest declares. What a run outside the `vep`
//! binary writes on standard error comes from the library run in this test
//! binary, started again as a child process.

mod common;

use std::collections::BTreeSet;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, ExitStatus, Output, Stdio};
use std::time::{Duration, Instant};

use noodles_bgzf::gzi;

use common::{
    check_html_report, compare_summary_sections, corpora, execute_manifest_run, manifest_runs,
    stats_html_report_path, stats_report_path, stats_text_file, warning_entries, Corpus,
    WarningEntry,
};

fn vep() -> Command {
    common::isolated_vep_anywhere()
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

/// A corpus's arguments with `input` read in place of its own records.
fn with_input(mut args: Vec<String>, input: &Path) -> Vec<String> {
    let input_flag = args.iter().position(|a| a == "-i").unwrap();
    args[input_flag + 1] = input.display().to_string();
    args
}

/// The golden corpus's arguments with `input` read in place of its own records.
fn corpus_args_with_input(input: &Path) -> Vec<String> {
    with_input(corpus_args(), input)
}

/// A corpus's arguments without `--force_overwrite`, for a run that must find
/// an existing file in its way.
fn without_force_overwrite(mut args: Vec<String>) -> Vec<String> {
    args.retain(|a| a != "--force_overwrite");
    args
}

/// The GRCh37-hgvs corpus's reference FASTA decompressed into `dir` with its
/// index beside it: an indexed FASTA a run can name, for a format that needs one.
fn decompressed_fasta(dir: &Path) -> PathBuf {
    let corpus =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37-hgvs");
    let compressed = std::fs::File::open(corpus.join("reference.fa.gz")).unwrap();
    let mut data = Vec::new();
    // The reference is bgzipped, one gzip member per block.
    flate2::read::MultiGzDecoder::new(compressed)
        .read_to_end(&mut data)
        .unwrap();
    let target = dir.join("reference.fa");
    std::fs::write(&target, data).unwrap();
    std::fs::copy(
        corpus.join("reference.fa.fai"),
        dir.join("reference.fa.fai"),
    )
    .unwrap();
    target
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

/// Runs `vep` with `args` writing to `out`, standard input closed.
fn run_with(args: &[String], extra: &[&str], out: &Path) -> Output {
    vep()
        .args(args)
        .args(extra)
        .args(["-o", out.to_str().unwrap()])
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

/// The usage wins over every other token of the command line (vep:29 before
/// vep:30, under Getopt::Long's `pass_through`, Config.pm:62): `--help` beside a
/// flag refused before parsing (`--af_exac`), an unknown flag, a bare word, a
/// refused value (`--fork x`), an integer option left without its value (`--fork --help`,
/// whose `--help` Getopt::Long pushes back) or a value flag of this version's
/// own left without its value (`--pipeline_depth --help`, which Ensembl VEP
/// leaves in `@ARGV` whole and prints its usage after), in either order, exits 0
/// with the text of `--help` alone, as Ensembl VEP 116.3 prints its usage and
/// exits 0 on each of these lines. `--help=1` is refused on both engines
/// (vep:30; exit 2 here, in clap's words). A string option takes `--help` as
/// its value on both engines (`input_file=s`, `cell_type=s`, Config.pm:209),
/// so `--input_file --help` and `--cell_type --help` print no usage and the
/// run goes on to its next refusal, the cache it cannot find here (exit 1).
#[test]
fn help_and_version_win_over_a_flag_refused_before_parsing() {
    let expected = text(&vep().arg("--help").output().unwrap().stdout);
    for argv in [
        &["--af_exac", "--help"][..],
        &["--help", "--af_exac"][..],
        &["--bogus", "--help"][..],
        &["extra", "--help"][..],
        &["--fork", "x", "--help"][..],
        &["--fork", "--help"][..],
        &["--buffer_size", "--help"][..],
        &["--pipeline_depth", "--help"][..],
        &["--writer_strategy", "--help"][..],
        &["--help", "--input_file"][..],
    ] {
        let output = vep().args(argv).output().unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{argv:?}: {}",
            text(&output.stderr)
        );
        assert_eq!(text(&output.stdout), expected, "{argv:?}");
    }
    let output = vep().arg("--help=1").output().unwrap();
    assert_eq!(output.status.code(), Some(2), "{}", text(&output.stderr));
    assert!(text(&output.stderr).starts_with("error:"));
    assert!(output.stdout.is_empty(), "no usage");
    for argv in [
        &["--input_file", "--help"][..],
        &["--cell_type", "--help"][..],
    ] {
        let output = vep().args(argv).output().unwrap();
        assert_eq!(
            output.status.code(),
            Some(1),
            "{argv:?}: {}",
            text(&output.stderr)
        );
        assert!(
            text(&output.stderr).contains("ERROR: [cache] Cache directory"),
            "{argv:?}: {}",
            text(&output.stderr)
        );
        assert!(
            !text(&output.stdout).contains("Usage"),
            "{argv:?}: no usage"
        );
    }
}

/// The usage error `Config::from_argv` returns for a command line is what the
/// binary ends on for it: the error's text is the text the binary writes, on
/// standard output for the help and the version and on standard error for a
/// syntax error, with nothing on the other stream, and the error's exit code is
/// the binary's exit status.
#[test]
fn library_usage_error_carries_the_text_and_the_status_the_binary_ends_with() {
    use vep_cli::config::ParseFailure;

    for (flags, asks_for_help_or_version, status) in [
        (&["--help"][..], true, 0),
        (&["--version"][..], true, 0),
        (&["--help=1"][..], false, 2),
        (&["--fork", "x"][..], false, 2),
        (&["--bogus"][..], false, 2),
    ] {
        let output = vep().args(flags).stdin(Stdio::null()).output().unwrap();
        let usage = match vep_cli::Config::from_argv(["vep"].iter().chain(flags)) {
            Err(ParseFailure::Usage(usage)) => usage,
            Err(other) => panic!("{flags:?} is refused under the error contract: {other}"),
            Ok(_) => panic!("{flags:?} parses"),
        };
        assert_eq!(
            usage.is_help_or_version(),
            asks_for_help_or_version,
            "{flags:?}"
        );
        assert_eq!(usage.exit_code(), status, "{flags:?}");
        assert_eq!(output.status.code(), Some(status), "{flags:?}");
        let (written, silent) = if asks_for_help_or_version {
            (&output.stdout, &output.stderr)
        } else {
            (&output.stderr, &output.stdout)
        };
        assert_eq!(text(written), usage.to_string(), "{flags:?}");
        assert!(silent.is_empty(), "{flags:?}: {}", text(silent));
    }
}

/// This version's own usage spellings, which Ensembl VEP 116.3 refuses or lacks
/// (`-h` is a prefix of eight of its options, left in `@ARGV` and refused at
/// vep:30; it has no `--version`): `-h` prints the help, `--version` and `-V` the
/// version, each exit 0. `-h` and `--version` win over the rest of the line as
/// `--help` does, `--version` ahead of `--help` included, where the help is
/// printed; `-V` is the version as the only argument.
#[test]
fn short_help_and_version_are_this_versions_own_and_win_the_same_way() {
    let help = text(&vep().arg("--help").output().unwrap().stdout);
    let version = text(&vep().arg("--version").output().unwrap().stdout);
    for (argv, expected) in [
        (&["--af_exac=1", "-h"][..], &help),
        (&["-h", "--bogus"][..], &help),
        (&["--af_exac", "--version"][..], &version),
        (&["--version", "--af_exac"][..], &version),
        (&["-V"][..], &version),
        (&["--bogus", "--version"][..], &version),
        (&["--fork", "--version"][..], &version),
        (&["--version", "--help"][..], &help),
    ] {
        let output = vep().args(argv).output().unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{argv:?}: {}",
            text(&output.stderr)
        );
        assert_eq!(&text(&output.stdout), expected, "{argv:?}");
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
    let upper_stdout = text(&upper.stdout);
    let upper_rows = data_lines(&upper_stdout);
    assert!(upper_rows.len() > 30, "{} rows", upper_rows.len());
    for spelling in ["stdout", "StdOut", "Stdout"] {
        let output = run_to(&["-o", spelling, "--quiet"], dir.path());
        assert_eq!(
            output.status.code(),
            Some(0),
            "{spelling}: {}",
            text(&output.stderr)
        );
        let stdout = text(&output.stdout);
        assert_eq!(data_lines(&stdout), upper_rows, "{spelling}");
        assert!(
            !dir.path().join(spelling).exists(),
            "{spelling}: a file was written instead"
        );
    }
}

/// `-o -` writes the annotation to standard output and no file, as Ensembl VEP's
/// two-argument `open(">-")` does (BaseRunner.pm:226), while the warnings file is
/// named after the value as spelled (`-_warnings.txt`, BaseVEP.pm:890); under
/// `--compress_output` the shell redirection `>-` (BaseRunner.pm:217) writes a
/// file named `-` on both engines.
#[test]
fn dash_output_file_is_standard_output_on_the_plain_path() {
    let dir = tempfile::tempdir().unwrap();
    let reference = run_to(&["-o", "STDOUT", "--quiet"], dir.path());
    assert_eq!(
        reference.status.code(),
        Some(0),
        "{}",
        text(&reference.stderr)
    );
    let reference_stdout = text(&reference.stdout);
    let reference_rows = data_lines(&reference_stdout);
    let output = run_to(&["-o", "-", "--quiet"], dir.path());
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let stdout = text(&output.stdout);
    let rows = data_lines(&stdout);
    assert_eq!(rows, reference_rows);
    assert!(rows.len() > 30);
    assert!(!dir.path().join("-").exists(), "a file was written instead");
    assert!(dir.path().join("-_warnings.txt").is_file());

    let compressed = run_to(
        &["-o", "-", "--quiet", "--compress_output", "gzip"],
        dir.path(),
    );
    assert_eq!(
        compressed.status.code(),
        Some(0),
        "{}",
        text(&compressed.stderr)
    );
    assert!(
        compressed.stdout.is_empty(),
        "the compressed output is a file"
    );
    let mut file_text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(dir.path().join("-")).unwrap())
        .read_to_string(&mut file_text)
        .unwrap();
    assert_eq!(
        data_lines(&file_text),
        reference_rows,
        "the file named `-` holds the rows"
    );
}

/// The existing-file check runs on the `-o` value as spelled, before it is read
/// as the `STDOUT` literal (BaseRunner.pm:203 before :213 and :222): with a file
/// named `STDOUT` in the working directory, `-o STDOUT` is refused with Ensembl
/// VEP's text, exit 1, and with `--force_overwrite` the rows go to standard
/// output while the file keeps its bytes.
#[test]
fn output_file_named_stdout_is_checked_as_spelled() {
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("STDOUT");
    std::fs::write(&file, "kept\n").unwrap();
    let refused = vep()
        .args(without_force_overwrite(corpus_args()))
        .args(["-o", "STDOUT", "--quiet"])
        .current_dir(dir.path())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(refused.status.code(), Some(1), "{}", text(&refused.stderr));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [output] Output file STDOUT already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite"
    );
    assert!(refused.stdout.is_empty());
    assert_eq!(std::fs::read_to_string(&file).unwrap(), "kept\n");

    let forced = run_to(&["-o", "STDOUT", "--quiet"], dir.path());
    assert_eq!(forced.status.code(), Some(0), "{}", text(&forced.stderr));
    assert!(data_lines(&text(&forced.stdout)).len() > 30);
    assert_eq!(std::fs::read_to_string(&file).unwrap(), "kept\n");
}

/// `-o` given twice keeps the last value and warns, as the `vep` script does
/// (vep:32-33: `warn("WARNING: More than one output file command-line flag
/// provided. Using the last value provided: '...' \n")`, a space before the
/// newline), under `--quiet` too, since the script warns before any option is
/// read; the first path is never written.
#[test]
fn repeated_output_file_keeps_last_with_warning() {
    let dir = tempfile::tempdir().unwrap();
    let first = dir.path().join("first.txt");
    let second = dir.path().join("second.txt");
    let warning = format!(
        "WARNING: More than one output file command-line flag provided. Using the last value provided: '{}' \n",
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
        assert!(
            stderr.starts_with(&warning),
            "{quiet:?}: expected stderr to open with {warning:?}, got {stderr:?}"
        );
        std::fs::remove_file(&second).unwrap();
    }
    let once = run_to(&["-o", first.to_str().unwrap()], dir.path());
    assert!(
        !text(&once.stderr).contains("More than one output file"),
        "a single -o is not warned about"
    );
}

/// `-o stdout` in any case streams the annotation to standard output, the status
/// lines stay on standard error, and the files written are the warnings file the
/// corpus's skipped records call for and the statistics report written by
/// default, `<spelling>_warnings.txt` and `<spelling>_summary.html` in the
/// working directory (BaseVEP.pm:890 and BaseRunner.pm:254 name them after the
/// output file as given).
#[test]
fn stdout_any_case_writes_stdout() {
    let dir = tempfile::tempdir().unwrap();
    let output = vep()
        .args(without_no_stats(corpus_args()))
        .args(["-o", "StdOut"])
        .current_dir(dir.path())
        .stdin(Stdio::null())
        .output()
        .unwrap();
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
        files_in(dir.path()),
        ["StdOut_summary.html", "StdOut_warnings.txt"],
        "the report and the warnings file are written to the working directory"
    );
    assert!(
        std::fs::read_to_string(dir.path().join("StdOut_summary.html"))
            .unwrap()
            .starts_with("<html>\n")
    );
}

/// An output file that exists is refused before the run with Ensembl VEP's text
/// (BaseRunner.pm:203-204) under the output component; `--force_overwrite` lifts it.
#[test]
fn existing_output_is_refused_with_ensembl_veps_text() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("taken.txt");
    std::fs::write(&out, "kept\n").unwrap();
    let refused = vep()
        .args(without_force_overwrite(corpus_args()))
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

/// A one-record VCF under the `##fileformat` lines `headers`, before its `#CHROM`
/// line, written to `dir` as `name`.
fn vcf_with_fileformat(dir: &Path, name: &str, headers: &[&str]) -> PathBuf {
    let path = dir.join(name);
    let mut text = headers.iter().map(|h| format!("{h}\n")).collect::<String>();
    text.push_str(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30282169\tstop_ins4_twostops_rev\tC\tCTATT\t.\t.\t.\n",
    );
    std::fs::write(&path, text).unwrap();
    path
}

/// Ensembl VEP's ensembl-io parser reads the `##fileformat` value of every such
/// header line (BaseVCF4.pm:93-94, :116-130) and ends the run on a version whose
/// major number is not its own 4: `confess "The VCF file format version
/// $f_version is not compatible with the parser version (VCF v4.2)"` (`:125`),
/// `$f_version` the first `<digits>.<digits>` of the value as written, so
/// `VCFv003.3` reports `003.3`; and on a value with no such pair: `die "The
/// script can't read the VCF file format version of 'fileformat'"` (`:129`). The
/// key follows `##` and any blanks (`:93`), and a second `##fileformat` line is
/// read like the first. Ensembl VEP 116.3 exits 25 on each of these inputs with
/// that text and, dying in `Runner::init` (BaseRunner.pm:406) before its output
/// opens, writes no output file and leaves an existing one as it was; vep-rs
/// refuses them under the error contract with the output file untouched.
#[test]
fn vcf_fileformat_major_version_other_than_4_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let incompatible = |version: &str| {
        format!(
            "ERROR: [input] The VCF file format version {version} is not compatible with the parser version (VCF v4.2)"
        )
    };
    let unreadable =
        "ERROR: [input] The script can't read the VCF file format version of 'fileformat'"
            .to_string();
    let cases: [(&[&str], String); 7] = [
        (&["##fileformat=VCFv3.3"], incompatible("3.3")),
        (&["##fileformat=VCFv5.0"], incompatible("5.0")),
        (&["##fileformat=VCFv003.3"], incompatible("003.3")),
        (&["## fileformat=VCFv3.3"], incompatible("3.3")),
        (
            &["##fileformat=VCFv4.2", "##fileformat=VCFv3.3"],
            incompatible("3.3"),
        ),
        (&["##fileformat=VCFv4"], unreadable.clone()),
        (&["##fileformat=VCF"], unreadable),
    ];
    for (i, (headers, expected)) in cases.iter().enumerate() {
        let input = vcf_with_fileformat(dir.path(), &format!("in{i}.vcf"), headers);
        let out = dir.path().join(format!("out{i}.txt"));
        let output = run_with(&corpus_args_with_input(&input), &["--quiet"], &out);
        assert_eq!(output.status.code(), Some(1), "{headers:?}");
        assert_eq!(text(&output.stderr).trim_end(), *expected, "{headers:?}");
        assert!(output.stdout.is_empty(), "{headers:?}");
        assert!(!out.exists(), "{headers:?}: an output file was written");
    }
    // An existing output file survives the refusal under the corpus arguments'
    // `--force_overwrite`.
    let input = vcf_with_fileformat(dir.path(), "kept.vcf", &["##fileformat=VCFv3.3"]);
    let out = dir.path().join("kept.txt");
    std::fs::write(&out, "an earlier result\n").unwrap();
    let output = run_with(&corpus_args_with_input(&input), &["--quiet"], &out);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        "an earlier result\n"
    );
}

/// The minor version is not read (BaseVCF4.pm:126 is commented out), the major
/// number is compared as a number, only the first `<digits>.<digits>` of the
/// value counts, `read_metadata` matches a `##key=value` line with a non-empty
/// value and the key `fileformat` as spelled (`:93`, `:116`), and a file without
/// the line is not checked: `VCFv4.5`, `VCFv4.5.1`, `VCFv04.2`, a bare `4.2`, an
/// empty value, `##FILEFORMAT=VCFv3.3` and no `##fileformat` line each annotate
/// to the rows a `VCFv4.2` header gives, with no warning, as Ensembl VEP 116.3
/// annotates them (exit 0, the same row on each).
#[test]
fn vcf_fileformat_4x_and_a_missing_line_annotate_as_vcfv42() {
    let dir = tempfile::tempdir().unwrap();
    let run = |name: &str, headers: &[&str]| -> (Vec<String>, String) {
        let input = vcf_with_fileformat(dir.path(), &format!("{name}.vcf"), headers);
        let out = dir.path().join(format!("{name}.txt"));
        let output = run_with(&corpus_args_with_input(&input), &[], &out);
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        let rows = data_lines(&std::fs::read_to_string(&out).unwrap())
            .into_iter()
            .map(str::to_string)
            .collect();
        (rows, text(&output.stderr))
    };
    let (expected, _) = run("v42", &["##fileformat=VCFv4.2"]);
    assert!(!expected.is_empty());
    let accepted: [&[&str]; 7] = [
        &["##fileformat=VCFv4.5"],
        &["##fileformat=VCFv4.5.1"],
        &["##fileformat=VCFv04.2"],
        &["##fileformat=4.2"],
        &["##fileformat="],
        &["##FILEFORMAT=VCFv3.3"],
        &[],
    ];
    for (i, headers) in accepted.iter().enumerate() {
        let (rows, stderr) = run(&format!("accepted{i}"), headers);
        assert_eq!(rows, expected, "{headers:?}");
        assert!(!stderr.contains("WARNING"), "{headers:?}: {stderr}");
    }
}

/// `--show_cache_info` prints the cache's source versions as `<key>\t<value>`
/// lines sorted by key and exits 0 (vep:39-43 prints them in hash order), with
/// the cache resolved by the same flags as a run; a named input file that does
/// not exist is refused first (BaseRunner.pm:406, Parser.pm:239), and an
/// existing one is not read. The info prints before `post_setup_checks` would
/// run (vep:39-46 exits before `run()`), so `--hgvs` with no FASTA in reach does
/// not stop it.
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
        .args(["--show_cache_info", "--hgvs"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(text(&output.stdout), expected);
    let refused = vep()
        .args(corpus_args_with_input(Path::new("/no/such/input.vcf")))
        .args(["--show_cache_info", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [input] File \"/no/such/input.vcf\" does not exist"
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
    // The cache carries no checksum manifest, so each run warns and the warnings
    // file lands in the working directory, kept out of the source tree here.
    let output = vep()
        .current_dir(root.path())
        .args(["--show_cache_info", "--cache", "--dir_cache"])
        .arg(root.path())
        .args(["--species", "homo_sapiens", "--cache_version", "116"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(text(&output.stdout), "assembly\tGRCh38.p14\ndbSNP\t156\n");
    let refused = vep()
        .current_dir(root.path())
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
        .current_dir(root.path())
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

/// `--show_cache_info` prints and exits 0 before the refusals of `Runner::init`
/// (Runner.pm:109-140), which vep:39-43 exits ahead of: an unreadable
/// `--synonyms` file (read at :131), a plugin that cannot be loaded under
/// `--safe` (:125), a `--fasta` that is not there (:134) and an output path that
/// cannot be written (opened in `run()`, BaseRunner.pm:226) stop a run and not
/// the print.
#[test]
fn show_cache_info_prints_before_the_init_stage_refusals() {
    let dir = tempfile::tempdir().unwrap();
    let expected = text(&run_control(&["--show_cache_info"], dir.path()).stdout);
    assert!(expected.lines().count() > 5, "{expected}");
    for extra in [
        &["--synonyms", "/no/such/synonyms.txt"][..],
        &["--plugin", "NoSuchPlugin", "--safe"][..],
        &["--fasta", "/no/such/reference.fa"][..],
        &["-o", "/no/such/dir/out.txt"][..],
    ] {
        let mut with_info = extra.to_vec();
        with_info.push("--show_cache_info");
        let output = run_control(&with_info, dir.path());
        assert_eq!(
            output.status.code(),
            Some(0),
            "{extra:?}: {}",
            text(&output.stderr)
        );
        assert_eq!(text(&output.stdout), expected, "{extra:?}");
        let refused = run_control(extra, dir.path());
        assert_eq!(
            refused.status.code(),
            Some(1),
            "{extra:?} stops a run: {}",
            text(&refused.stderr)
        );
    }
}

/// Under `--regulatory_gff` the print carries a `regulatory_gff\t<file>` line
/// beside the cache's source versions, the first file of the flag: Ensembl VEP
/// 116.3 prints the version data of every annotation source, first source first
/// (BaseRunner.pm:415-417), and a regulatory GFF source's is its file
/// (RegFeat.pm:674-677). The sources are built before the print, so a cache
/// holding regulatory data of its own draws the replacement warning into the
/// warnings file there too (CacheDir.pm:159-165).
#[test]
fn show_cache_info_lists_the_regulatory_gff_file() {
    let dir = tempfile::tempdir().unwrap();
    let plain = run_control(&["--show_cache_info"], dir.path());
    assert_eq!(plain.status.code(), Some(0), "{}", text(&plain.stderr));
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/regulatory-gff");
    let regulatory_features = fixtures.join("regulatory_features.gff3.gz");
    let gff = format!(
        "file={},emars={},motifs={}",
        regulatory_features.display(),
        fixtures.join("EMARs.gff.gz").display(),
        fixtures.join("motif_features.gff3.gz").display()
    );
    let output = run_control(
        &["--show_cache_info", "--quiet", "--regulatory_gff", &gff],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let mut expected: Vec<String> = text(&plain.stdout).lines().map(str::to_owned).collect();
    expected.push(format!("regulatory_gff\t{}", regulatory_features.display()));
    expected.sort();
    assert_eq!(text(&output.stdout).lines().collect::<Vec<_>>(), expected);
    let warnings =
        std::fs::read_to_string(dir.path().join("variant_effect_output.txt_warnings.txt"))
            .unwrap_or_default();
    assert!(
        warnings.contains("WARNING: Using regulatory features from --regulatory_gff; ignoring the regulatory data in this cache"),
        "{warnings}"
    );
}

/// A `--regulatory_gff` file that is not there ends a `--show_cache_info` run
/// before anything prints: Ensembl VEP 116.3 builds the regulatory GFF sources
/// with the cache's, inside the chain the print reads from
/// (AnnotationSourceAdaptor.pm:102, :264-267), and refuses the file by key and
/// name.
#[test]
fn show_cache_info_refuses_a_missing_regulatory_gff_before_printing() {
    let dir = tempfile::tempdir().unwrap();
    let refused = run_control(
        &[
            "--show_cache_info",
            "--quiet",
            "--regulatory_gff",
            "file=/no/such/reg.gff3.gz",
        ],
        dir.path(),
    );
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [cli] --regulatory_gff file=/no/such/reg.gff3.gz not found"
    );
    assert!(refused.stdout.is_empty());
}

/// The custom corpus's indexed BED source, with its `.tbi` beside it.
fn indexed_bed() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-custom/custom/regions.bed.gz")
}

/// A `--custom` file that cannot be opened ends a `--show_cache_info` run
/// before anything prints, as a run without the flag refuses it: Ensembl VEP
/// 116.3 builds the custom sources with the cache's inside the chain the print
/// reads from (vep:40, BaseRunner.pm:406, Runner.pm:827,
/// AnnotationSourceAdaptor.pm:103) and opens each file's tabix index for its
/// chromosome names there (BaseRunner.pm:453, File.pm:371), so a file that is
/// missing, bgzipped without an index, or not bgzipped dies with nothing on
/// standard output (ensembl-io TabixParser.pm:52, :63). The refusal is
/// vep-rs's own line, exit 1, and no warnings file is written.
#[test]
fn show_cache_info_refuses_an_unopenable_custom_file_before_printing() {
    let dir = tempfile::tempdir().unwrap();
    let unindexed = dir.path().join("unindexed.bed.gz");
    std::fs::copy(indexed_bed(), &unindexed).unwrap();
    let plain = dir.path().join("plain.bed");
    std::fs::write(&plain, "21\t11000000\t11000100\tregion\n").unwrap();
    for file in [
        "/no/such/regions.bed.gz".to_string(),
        unindexed.display().to_string(),
        plain.display().to_string(),
    ] {
        let custom = format!("file={file},short_name=regions,format=bed");
        let refused = run_control(
            &["--show_cache_info", "--quiet", "--custom", &custom],
            dir.path(),
        );
        assert_eq!(refused.status.code(), Some(1), "{custom}");
        assert_eq!(
            text(&refused.stderr).trim_end(),
            format!(
                "ERROR: [custom] cannot open the --custom source {file}: \
                 tabix index not found for {file}. Expected {file}.csi or {file}.tbi"
            ),
            "{custom}"
        );
        assert!(
            refused.stdout.is_empty(),
            "{custom}: the cache info was printed"
        );
        assert!(
            !dir.path()
                .join("variant_effect_output.txt_warnings.txt")
                .exists(),
            "{custom}: the warnings file was created"
        );
    }
}

/// `--show_cache_info` opens the `--custom` sources as a run does and prints
/// the cache's lines alone: a custom source's info is its `custom_info` entry,
/// which the version data skips (File.pm:170, BaseRunner.pm:417). The FASTA is
/// reached only through a `###CHR###` file, whose chromosome list it supplies
/// (AnnotationSourceAdaptor.pm:440, BaseVEP.pm:568, FastaSequence.pm:187): a
/// `--fasta` that is not there stops the print under such a file and not
/// otherwise, and a pattern matching no file is refused with the text of
/// AnnotationSourceAdaptor.pm:452.
#[test]
fn show_cache_info_prints_the_cache_lines_beside_an_opened_custom_source() {
    let dir = tempfile::tempdir().unwrap();
    let plain = run_control(&["--show_cache_info", "--quiet"], dir.path());
    assert_eq!(plain.status.code(), Some(0), "{}", text(&plain.stderr));
    let expected = text(&plain.stdout);
    assert!(expected.lines().count() > 5, "{expected}");
    let indexed = format!(
        "file={},short_name=regions,format=bed",
        indexed_bed().display()
    );
    std::fs::copy(indexed_bed(), dir.path().join("21.bed.gz")).unwrap();
    std::fs::copy(
        format!("{}.tbi", indexed_bed().display()),
        dir.path().join("21.bed.gz.tbi"),
    )
    .unwrap();
    let per_chromosome = format!(
        "file={},short_name=regions,format=bed",
        dir.path().join("###CHR###.bed.gz").display()
    );
    let no_such_fasta = ["--fasta", "/no/such/reference.fa"];
    for extra in [
        &["--custom", &indexed][..],
        &["--custom", &indexed, no_such_fasta[0], no_such_fasta[1]][..],
        &["--custom", &per_chromosome][..],
    ] {
        let mut with_info = vec!["--show_cache_info", "--quiet"];
        with_info.extend_from_slice(extra);
        let output = run_control(&with_info, dir.path());
        assert_eq!(
            output.status.code(),
            Some(0),
            "{extra:?}: {}",
            text(&output.stderr)
        );
        assert_eq!(text(&output.stdout), expected, "{extra:?}");
    }
    let refused = run_control(
        &[
            "--show_cache_info",
            "--quiet",
            "--custom",
            &per_chromosome,
            no_such_fasta[0],
            no_such_fasta[1],
        ],
        dir.path(),
    );
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [fasta] Specified FASTA file/directory /no/such/reference.fa not found"
    );
    assert!(refused.stdout.is_empty());
    let unmatched = dir.path().join("###CHR###.vcf.gz");
    let refused = run_control(
        &[
            "--show_cache_info",
            "--quiet",
            "--custom",
            &format!("file={},short_name=sites,format=vcf", unmatched.display()),
        ],
        dir.path(),
    );
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        format!(
            "ERROR: [custom] No files with pattern {} were found",
            unmatched.display()
        )
    );
    assert!(refused.stdout.is_empty());
}

/// A warning drawn before the print goes to the warnings file, and one that
/// cannot be written ends a `--show_cache_info` run with nothing printed, as
/// Ensembl VEP dies on the warning itself (BaseVEP.pm:817, :894); without a
/// warning the same unwritable path still prints, as
/// `show_cache_info_prints_before_the_init_stage_refusals` holds. The
/// `--fasta_dir` with no FASTA for the species warns through the `###CHR###`
/// file's chromosome list (BaseVEP.pm:516, :568), the fractional
/// `num_records` through the source's open, `--regulatory_gff` over a cache
/// holding regulatory data while the sources are built (CacheDir.pm:159-165),
/// and a cache without a checksum manifest at its resolution, a warning of
/// vep-rs's own. The regulatory warning comes before any custom file is
/// opened, so its unwritten file wins over a missing custom file, as the
/// order of Ensembl VEP's chain has it.
#[test]
fn show_cache_info_refuses_an_unwritable_warnings_file_once_a_warning_is_drawn() {
    let dir = tempfile::tempdir().unwrap();
    let empty_fasta_dir = dir.path().join("no-fasta");
    std::fs::create_dir(&empty_fasta_dir).unwrap();
    std::fs::copy(indexed_bed(), dir.path().join("21.bed.gz")).unwrap();
    std::fs::copy(
        format!("{}.tbi", indexed_bed().display()),
        dir.path().join("21.bed.gz.tbi"),
    )
    .unwrap();
    let per_chromosome = format!(
        "file={},short_name=regions,format=bed",
        dir.path().join("###CHR###.bed.gz").display()
    );
    let fractional = format!(
        "file={},short_name=regions,format=bed,num_records=2.5",
        indexed_bed().display()
    );
    let regulatory_gff = format!(
        "file={}",
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/regulatory-gff/regulatory_features.gff3.gz")
            .display()
    );
    let unwritten =
        "ERROR: [output] Could not write to warnings file /no/such/dir/out.txt_warnings.txt";
    for extra in [
        &[
            "--custom",
            &per_chromosome,
            "--fasta_dir",
            empty_fasta_dir.to_str().unwrap(),
        ][..],
        &["--custom", &fractional][..],
        &[
            "--regulatory_gff",
            &regulatory_gff,
            "--custom",
            "file=/no/such/regions.bed.gz,short_name=regions,format=bed",
        ][..],
    ] {
        let mut args = vec!["--show_cache_info", "--quiet", "-o", "/no/such/dir/out.txt"];
        args.extend_from_slice(extra);
        let refused = run_control(&args, dir.path());
        assert_eq!(refused.status.code(), Some(1), "{extra:?}");
        assert_eq!(text(&refused.stderr).trim_end(), unwritten, "{extra:?}");
        assert!(
            refused.stdout.is_empty(),
            "{extra:?}: the cache info was printed"
        );
    }
    let root = tempfile::tempdir().unwrap();
    let leaf = root.path().join("homo_sapiens").join("116_GRCh38");
    std::fs::create_dir_all(leaf.join("transcripts")).unwrap();
    std::fs::write(
        leaf.join("info.json"),
        r#"{"assembly":"GRCh38","source_versions":{"dbSNP":"156"}}"#,
    )
    .unwrap();
    let refused = vep()
        .current_dir(root.path())
        .args(["--show_cache_info", "--quiet", "--cache", "--dir_cache"])
        .arg(root.path())
        .args(["--species", "homo_sapiens", "--cache_version", "116"])
        .args(["-o", "/no/such/dir/out.txt"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(text(&refused.stderr).trim_end(), unwritten);
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

/// The footer of a Parquet output's first part: its key-value metadata and
/// whether its row groups declare a sort.
#[cfg(feature = "parquet")]
fn parquet_footer(out: &Path) -> (std::collections::BTreeMap<String, String>, Vec<bool>) {
    use parquet::file::reader::FileReader;
    let parts = vep_io::output::parquet::part_files(out).unwrap();
    let reader = parquet::file::serialized_reader::SerializedFileReader::new(
        std::fs::File::open(&parts[0]).unwrap(),
    )
    .unwrap();
    let metadata = reader.metadata();
    let kv = metadata
        .file_metadata()
        .key_value_metadata()
        .unwrap()
        .iter()
        .map(|kv| (kv.key.clone(), kv.value.clone().unwrap_or_default()))
        .collect();
    let sorted = metadata
        .row_groups()
        .iter()
        .map(|row_group| row_group.sorting_columns().is_some())
        .collect();
    (kv, sorted)
}

/// A VCF of intergenic records on the release corpus's chromosomes, one per
/// `(chromosome, position)`.
fn intergenic_vcf(dir: &Path, name: &str, records: &[(&str, u64)]) -> PathBuf {
    let path = dir.join(name);
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (i, (chrom, pos)) in records.iter().enumerate() {
        text.push_str(&format!("{chrom}\t{pos}\tr{}\tA\tG\t.\t.\t.\n", i + 1));
    }
    std::fs::write(&path, text).unwrap();
    path
}

/// A Parquet run of the release corpus's arguments over `input`, `extra` flags
/// added, written to `out`.
#[cfg(feature = "parquet")]
fn parquet_run(input: &Path, extra: &[&str], out: &Path) -> Output {
    let mut args = corpus_args_with_input(input);
    args.extend(["--output_format", "parquet", "--quiet"].map(str::to_string));
    run_with(&args, extra, out)
}

/// A Parquet run refuses the first record out of position order on its
/// chromosome, naming the record's line and the two positions, exits 1 and
/// leaves no file; a text format annotates the same input in input order.
#[cfg(feature = "parquet")]
#[test]
fn parquet_refuses_an_unsorted_input_naming_the_line() {
    let dir = tempfile::tempdir().unwrap();
    let input = intergenic_vcf(dir.path(), "in.vcf", &[("21", 100), ("21", 50)]);
    let out = dir.path().join("out.parquet");
    let output = parquet_run(&input, &[], &out);
    assert_eq!(output.status.code(), Some(1), "{}", text(&output.stderr));
    assert_eq!(
        text(&output.stderr).trim_end(),
        "ERROR: [parquet] input is not sorted by chromosome and position at line 2 (21:50 after 21:100); sort the input or pass --no_check_variants_order"
    );
    assert!(!out.exists(), "the refused run leaves no file");

    let tab = dir.path().join("out.tab");
    let mut args = corpus_args_with_input(&input);
    args.extend(["--tab", "--quiet"].map(str::to_string));
    let output = run_with(&args, &[], &tab);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let written = std::fs::read_to_string(&tab).unwrap();
    let locations: Vec<&str> = data_lines(&written)
        .iter()
        .map(|l| l.split('\t').nth(1).unwrap())
        .collect();
    assert_eq!(locations, ["21:100", "21:50"]);
}

/// The input lines the warnings file at `path` names as skipped, in file
/// order, and nothing when the file is not there.
#[cfg(any(feature = "parquet", target_os = "linux"))]
fn skipped_lines(path: &Path) -> Vec<u64> {
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .map(|line| {
            let number = line.strip_prefix("WARNING: line ").expect(line);
            number.split(' ').next().unwrap().parse().unwrap()
        })
        .collect()
}

/// A streaming Parquet run refused for a record out of order leaves what the
/// batched pipeline mode leaves, at any worker count: no output, and the
/// warnings of the whole batch that holds the refused record, in input order.
/// Of 96 lines in one batch, line 20 is out of order and lines 3, 22, 28 and
/// 40 are skipped. One worker cuts the batch into tasks of 16 lines and more
/// workers into tasks of 8, so the tasks that hold lines 28 and 40 come due
/// after the refused one, in whatever order they finished: the writer reports
/// the rest of the batch before it ends the run, and every run names all four.
#[cfg(feature = "parquet")]
#[test]
fn a_refused_parquet_run_leaves_the_whole_batchs_warnings_at_any_worker_count() {
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<(&str, u64)> = (1..=96u64)
        .map(|line| match line {
            3 | 22 | 28 | 40 => ("Z", 1_000 + line),
            20 => ("21", 5),
            _ => ("21", 100 + line),
        })
        .collect();
    let input = intergenic_vcf(dir.path(), "in.vcf", &records);
    let out = dir.path().join("out.parquet");
    let warnings = dir.path().join("out.parquet_warnings.txt");
    let left = |extra: &[&str]| {
        let _ = std::fs::remove_file(&warnings);
        let output = parquet_run(&input, extra, &out);
        assert_eq!(output.status.code(), Some(1), "{}", text(&output.stderr));
        assert!(!out.exists(), "the refused run leaves no output");
        (text(&output.stderr), skipped_lines(&warnings))
    };
    let batched = left(&["--pipeline_mode", "batched"]);
    assert_eq!(
        batched.0.trim_end(),
        "ERROR: [parquet] input is not sorted by chromosome and position at line 20 (21:5 after 21:119); sort the input or pass --no_check_variants_order"
    );
    assert_eq!(batched.1, [3, 22, 28, 40]);
    for fork in ["1", "2", "8", "32"] {
        for run in 0..6 {
            assert_eq!(
                left(&["--pipeline_mode", "streaming", "--fork", fork]),
                batched,
                "--fork {fork}, run {run}"
            );
        }
    }
    for run in 0..6 {
        assert_eq!(
            left(&["--pipeline_mode", "streaming"]),
            batched,
            "run {run}"
        );
    }
}

/// The flags of every run a run that ends early is repeated under: both
/// pipeline modes at 1, 2, 8 and 32 workers, `repeats` times each.
#[cfg(any(feature = "parquet", target_os = "linux"))]
fn modes_and_worker_counts(repeats: usize) -> Vec<[&'static str; 4]> {
    let mut runs = Vec::new();
    for mode in ["streaming", "batched"] {
        for fork in ["1", "2", "8", "32"] {
            for _ in 0..repeats {
                runs.push(["--pipeline_mode", mode, "--fork", fork]);
            }
        }
    }
    runs
}

/// A Parquet run of several batches refused for a record out of order leaves,
/// in each pipeline mode and at any worker count, the warnings of the batches
/// through the one that holds the refused record and of none after it, in
/// input order. Of 400 records in batches of 100, each batch holds one skipped
/// line (30, 130, 230 and 330) and line 150, in the second, is out of order:
/// every run names lines 30 and 130, whatever the batches it had read beyond.
#[cfg(feature = "parquet")]
#[test]
fn a_refused_parquet_run_of_several_batches_leaves_the_stated_warnings_in_each_mode() {
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<(&str, u64)> = (1..=400u64)
        .map(|line| match line {
            30 | 130 | 230 | 330 => ("Z", 1_000 + line),
            150 => ("21", 5),
            _ => ("21", 100 + line),
        })
        .collect();
    let input = intergenic_vcf(dir.path(), "in.vcf", &records);
    let out = dir.path().join("out.parquet");
    let warnings = dir.path().join("out.parquet_warnings.txt");
    for flags in modes_and_worker_counts(6) {
        let _ = std::fs::remove_file(&warnings);
        let mut extra = vec!["--buffer_size", "100"];
        extra.extend(flags);
        let output = parquet_run(&input, &extra, &out);
        assert_eq!(output.status.code(), Some(1), "{flags:?}");
        assert_eq!(
            text(&output.stderr).trim_end(),
            "ERROR: [parquet] input is not sorted by chromosome and position at line 150 (21:5 after 21:249); sort the input or pass --no_check_variants_order",
            "{flags:?}"
        );
        assert!(!out.exists(), "{flags:?}: the refused run leaves no output");
        assert_eq!(skipped_lines(&warnings), [30, 130], "{flags:?}");
    }
}

/// A Parquet run refused for a record out of order leaves the warnings its
/// reader found, a run of empty lines and, in a format the reader parses, a
/// skipped line, through the batch that holds the refused record and no
/// further, each in its place in input order, in both pipeline modes at any
/// worker count: the reader hands them on with their batch, however far ahead
/// of the writer it has read. Of 400 records in batches of 50, record 70, in
/// the second batch, is out of order, and a run of two empty lines follows the
/// record that fills that batch, so it is read once the batch is full and
/// travels with the third: record 100 of the VCF, and record 102 of the
/// Ensembl-format lines, whose batches count the records kept. The VCF skips
/// records 7 and 27, around the empty line after record 20, and record 57; the
/// Ensembl-format lines skip records 7 and 57, ahead of the empty line after
/// record 80. Both skip record 107, which the run does not reach.
#[cfg(feature = "parquet")]
#[test]
fn a_refused_parquet_run_leaves_the_readers_warnings_through_its_batch_in_input_order() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.parquet");
    let warnings = dir.path().join("out.parquet_warnings.txt");
    // 400 records with one empty line after record `empty_after` and a run of
    // two after record `run_after`.
    let lines = |skipped: &[u64],
                 empty_after: u64,
                 run_after: u64,
                 record: &dyn Fn(&str, u64) -> String| {
        let mut text = String::new();
        for line in 1..=400u64 {
            text.push_str(&match line {
                70 => record("21", 5),
                _ if skipped.contains(&line) => record("Z", 1_000 + line),
                _ => record("21", 100 + line),
            });
            if line == empty_after {
                text.push('\n');
            }
            if line == run_after {
                text.push_str("\n\n");
            }
        }
        text
    };
    let vcf = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}",
        lines(&[7, 27, 57, 107], 20, 100, &|chrom, pos| {
            format!("{chrom}\t{pos}\t.\tA\tG\t.\t.\t.\n")
        })
    );
    let ensembl = lines(&[7, 57, 107], 80, 102, &|chrom, pos| {
        format!("{chrom}\t{pos}\t{pos}\tA/G\t1\n")
    });
    // A skipped line as `line N`, any other warning as its text.
    let warned = || -> Vec<String> {
        std::fs::read_to_string(&warnings)
            .unwrap_or_default()
            .lines()
            .map(|line| {
                let text = line.strip_prefix("WARNING: ").expect(line);
                match text.split_once(" skipped (") {
                    Some((line_number, _)) => line_number.to_string(),
                    None => text.to_string(),
                }
            })
            .collect()
    };
    for (format, text_in, expected) in [
        (
            "vcf",
            &vcf,
            vec![
                "line 7",
                "Skipped 1 empty line from line 21",
                "line 27",
                "line 57",
            ],
        ),
        (
            "ensembl",
            &ensembl,
            vec!["line 7", "line 57", "Skipped 1 empty line from line 81"],
        ),
    ] {
        let input = dir.path().join(format!("in.{format}"));
        std::fs::write(&input, text_in).unwrap();
        for flags in modes_and_worker_counts(6) {
            let _ = std::fs::remove_file(&warnings);
            let mut extra = vec!["--format", format, "--buffer_size", "50"];
            extra.extend(flags);
            let output = parquet_run(&input, &extra, &out);
            assert_eq!(output.status.code(), Some(1), "{format} {flags:?}");
            assert_eq!(
                text(&output.stderr).trim_end(),
                "ERROR: [parquet] input is not sorted by chromosome and position at line 70 (21:5 after 21:169); sort the input or pass --no_check_variants_order",
                "{format} {flags:?}"
            );
            assert!(!out.exists(), "{format} {flags:?}");
            assert_eq!(warned(), expected, "{format} {flags:?}");
        }
    }
}

/// A text output that cannot be written ends the run with one error and one
/// warnings file in both pipeline modes at any worker count: the warnings of
/// every batch through the one whose write failed, in input order, and of no
/// batch after it. Of 400 records in batches of 50 every fifth is skipped, and
/// without its header the output is 40 rows of 65 to 67 bytes a batch: 7,954
/// bytes for three batches and 10,634 for four. The output is Linux's
/// `/dev/full`, which fails every write, so the buffered writer's first write,
/// of a buffer of 9,000 bytes, fails inside the fourth batch and the run names
/// the skipped lines through line 200.
#[cfg(target_os = "linux")]
#[test]
fn a_text_output_that_cannot_be_written_leaves_one_warnings_file_at_any_worker_count() {
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<(&str, u64)> = (1..=400u64)
        .map(|line| match line % 5 {
            3 => ("Z", 1_000 + line),
            _ => ("21", 100 + line),
        })
        .collect();
    let input = intergenic_vcf(dir.path(), "in.vcf", &records);
    let warnings = dir.path().join("warnings.txt");
    let mut args = corpus_args_with_input(&input);
    args.extend(
        [
            "--quiet",
            "--no_headers",
            "--buffer_size",
            "50",
            "--output_buffer_size",
            "9000",
            "--warning_file",
            warnings.to_str().unwrap(),
        ]
        .map(str::to_string),
    );
    let through_the_fourth_batch: Vec<u64> = (1..=200).filter(|line| line % 5 == 3).collect();
    for flags in modes_and_worker_counts(6) {
        let _ = std::fs::remove_file(&warnings);
        let output = run_with(&args, &flags, Path::new("/dev/full"));
        assert_eq!(output.status.code(), Some(1), "{flags:?}");
        assert_eq!(
            text(&output.stderr).trim_end(),
            "ERROR: [run] Failed to write output: No space left on device (os error 28)",
            "{flags:?}"
        );
        assert_eq!(
            skipped_lines(&warnings),
            through_the_fourth_batch,
            "{flags:?}"
        );
    }
}

/// Chromosomes may come in any order and each resets the position, as long as
/// a chromosome's records stay together, so a chromosome that returns after
/// another is the refusal; each chromosome's rows fill a part of their own, the
/// parts numbered in input order and each in `contig_index` order. The check
/// reads the position as written, so a deletion followed by a substitution at
/// its POS, whose annotated start precedes the deletion's, passes and the part
/// still sorts by `pos`.
#[cfg(feature = "parquet")]
#[test]
fn parquet_chromosomes_stay_together_in_any_order() {
    let dir = tempfile::tempdir().unwrap();
    let input = intergenic_vcf(
        dir.path(),
        "any_order.vcf",
        &[("22", 100), ("21", 500), ("21", 500)],
    );
    let out = dir.path().join("any_order.parquet");
    let output = parquet_run(&input, &[], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let parts = vep_io::output::parquet::part_files(&out).unwrap();
    assert_eq!(
        parts,
        [
            out.join("part-00000.parquet"),
            out.join("part-00001.parquet")
        ]
    );
    assert_eq!(
        vep_io::output::parquet::sort_keys(&parts[0]).unwrap(),
        [(1, 100)]
    );
    assert_eq!(
        vep_io::output::parquet::sort_keys(&parts[1]).unwrap(),
        [(0, 500), (0, 500)]
    );
    let (kv, sorted) = parquet_footer(&out);
    assert_eq!(kv["sort"], "contig_index,pos");
    assert!(sorted.iter().all(|declared| *declared));

    let input = intergenic_vcf(
        dir.path(),
        "returns.vcf",
        &[("21", 100), ("22", 5), ("21", 600)],
    );
    let out = dir.path().join("returns.parquet");
    let output = parquet_run(&input, &[], &out);
    assert_eq!(output.status.code(), Some(1), "{}", text(&output.stderr));
    assert_eq!(
        text(&output.stderr).trim_end(),
        "ERROR: [parquet] input is not sorted by chromosome and position at line 3 (21:600 after 22:5); sort the input or pass --no_check_variants_order"
    );
    assert!(!out.exists());

    let input = dir.path().join("anchor.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t100\tdel\tAT\tA\t.\t.\t.\n\
         21\t100\tsnv\tA\tG\t.\t.\t.\n\
         21\t101\tsnv2\tT\tC\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.path().join("anchor.parquet");
    let output = parquet_run(&input, &[], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let keys = vep_io::output::parquet::sort_keys(&out).unwrap();
    assert_eq!(
        keys.iter().map(|k| k.1).collect::<Vec<_>>(),
        [100, 101, 101]
    );
    let (kv, _) = parquet_footer(&out);
    assert_eq!(kv["sort"], "contig_index,pos");
}

/// `contig_index` ranks the cache's chromosomes first, in its order, then the
/// reference FASTA's other sequences in its index's order, so two records
/// `--dont_skip` keeps on contigs the cache lacks take distinct ranks; the
/// footer lists the order.
#[cfg(feature = "parquet")]
#[test]
fn parquet_contig_order_follows_the_cache_then_the_fasta_index() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("ref.fa");
    let mut sequences = String::new();
    let mut fai = String::new();
    for name in ["Z", "21", "Q"] {
        sequences.push_str(&format!(">{name}\n"));
        fai.push_str(&format!("{name}\t60\t{}\t60\t61\n", sequences.len()));
        sequences.push_str(&"A".repeat(60));
        sequences.push('\n');
    }
    std::fs::write(&fasta, sequences).unwrap();
    std::fs::write(dir.path().join("ref.fa.fai"), fai).unwrap();
    let input = intergenic_vcf(dir.path(), "in.vcf", &[("21", 10), ("Z", 5), ("Q", 7)]);
    let out = dir.path().join("out.parquet");
    let output = parquet_run(
        &input,
        &["--dont_skip", "--fasta", fasta.to_str().unwrap()],
        &out,
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let (kv, _) = parquet_footer(&out);
    assert_eq!(kv["contig_order"], "21,22,Z,Q");
    assert_eq!(kv["sort"], "contig_index,pos");
    let keys = vep_io::output::parquet::sort_keys(&out).unwrap();
    assert_eq!(keys, [(0, 10), (2, 5), (3, 7)]);

    let output = parquet_run(
        &input,
        &["--dont_skip"],
        &dir.path().join("no_fasta.parquet"),
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let (kv, _) = parquet_footer(&dir.path().join("no_fasta.parquet"));
    assert_eq!(kv["contig_order"], "21,22");
}

/// A reference sequence the cache holds under another name adds no entry to
/// the footer's `contig_order`. The GRCh38-custom cache's `chr_synonyms.txt`
/// links `NC_000021.9` and `CM000683.2` to its chromosome `21`, so a reference
/// index that names both accessions and a sequence the cache lacks lists `21`
/// once, the other sequence after it and neither accession, and a record
/// written with an accession ranks as `21`.
#[cfg(feature = "parquet")]
#[test]
fn parquet_contig_order_folds_a_reference_synonym_of_a_cache_chromosome() {
    let corpus =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-custom");
    let synonyms = std::fs::read_to_string(corpus.join("json_cache/chr_synonyms.txt")).unwrap();
    for accession in ["NC_000021.9", "CM000683.2"] {
        assert!(
            synonyms
                .lines()
                .any(|line| line == format!("21\t{accession}")),
            "the cache links {accession} to 21"
        );
    }
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("ref.fa");
    let mut sequences = String::new();
    let mut fai = String::new();
    for name in ["NC_000021.9", "Q", "CM000683.2"] {
        sequences.push_str(&format!(">{name}\n"));
        fai.push_str(&format!("{name}\t60\t{}\t60\t61\n", sequences.len()));
        sequences.push_str(&"A".repeat(60));
        sequences.push('\n');
    }
    std::fs::write(&fasta, sequences).unwrap();
    std::fs::write(dir.path().join("ref.fa.fai"), fai).unwrap();
    let input = intergenic_vcf(dir.path(), "in.vcf", &[("NC_000021.9", 10), ("Q", 7)]);
    let out = dir.path().join("out.parquet");
    let mut args = with_input(args_for(&corpus), &input);
    args.extend(["--output_format", "parquet", "--quiet"].map(str::to_string));
    let output = run_with(
        &args,
        &["--dont_skip", "--fasta", fasta.to_str().unwrap()],
        &out,
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let (kv, _) = parquet_footer(&out);
    assert_eq!(kv["contig_order"], "21,Q");
    assert!(
        !kv["contig_order"]
            .split(',')
            .any(|entry| entry == "NC_000021.9" || entry == "CM000683.2"),
        "{}",
        kv["contig_order"]
    );
    let keys = vep_io::output::parquet::sort_keys(&out).unwrap();
    assert_eq!(keys, [(0, 10), (1, 7)]);
}

/// The entries of a Parquet output directory as paths under it, sorted.
#[cfg(feature = "parquet")]
fn directory_entries(out: &Path) -> Vec<String> {
    fn walk(dir: &Path, root: &Path, entries: &mut Vec<String>) {
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                walk(&path, root, entries);
            } else {
                entries.push(
                    path.strip_prefix(root)
                        .unwrap()
                        .to_string_lossy()
                        .replace('\\', "/"),
                );
            }
        }
    }
    let mut entries = Vec::new();
    walk(out, out, &mut entries);
    entries.sort();
    entries
}

/// A Parquet run writes a directory of parts, one per chromosome, and a
/// manifest; a second run into it is refused as any existing output is, and
/// under `--force_overwrite` the previous run's parts and manifest go and the
/// directory holds the new run alone.
#[cfg(feature = "parquet")]
#[test]
fn parquet_force_overwrite_removes_only_vep_rs_parts() {
    use vep_io::output::parquet::Manifest;
    let dir = tempfile::tempdir().unwrap();
    let three = intergenic_vcf(dir.path(), "three.vcf", &[("21", 10), ("Z", 5), ("Q", 7)]);
    let out = dir.path().join("out.parquet");
    let output = parquet_run(&three, &["--dont_skip"], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        directory_entries(&out),
        [
            "_vep_rs_manifest.json",
            "part-00000.parquet",
            "part-00001.parquet",
            "part-00002.parquet"
        ]
    );
    let previous = Manifest::read(&out).unwrap();
    assert_eq!(previous.rows, 3);
    assert_eq!(
        previous
            .parts
            .iter()
            .map(|p| (p.path.as_str(), p.rows))
            .collect::<Vec<_>>(),
        [
            ("part-00000.parquet", 1),
            ("part-00001.parquet", 1),
            ("part-00002.parquet", 1)
        ]
    );
    assert_eq!(previous.vep_rs_version, env!("CARGO_PKG_VERSION"));

    let one = intergenic_vcf(dir.path(), "one.vcf", &[("21", 10), ("21", 20)]);
    let mut args = without_force_overwrite(corpus_args_with_input(&one));
    args.extend(["--output_format", "parquet", "--quiet"].map(str::to_string));
    let refused = run_with(&args, &[], &out);
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        format!(
            "ERROR: [output] Output file {} already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite",
            out.display()
        )
    );
    assert_eq!(
        directory_entries(&out).len(),
        4,
        "the refusal removes nothing"
    );

    let output = parquet_run(&one, &[], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        directory_entries(&out),
        ["_vep_rs_manifest.json", "part-00000.parquet"]
    );
    let manifest = Manifest::read(&out).unwrap();
    assert_eq!(manifest.rows, 2);
    assert_ne!(manifest.run_id, previous.run_id);
    assert_eq!(
        vep_io::output::parquet::sort_keys(&out).unwrap(),
        [(0, 10), (0, 20)]
    );
}

/// A directory holding anything vep-rs did not write, a text file or a Parquet
/// file without `vep_rs_version` in its footer, refuses `--force_overwrite`
/// with every file left in place.
#[cfg(feature = "parquet")]
#[test]
fn parquet_force_overwrite_refuses_a_foreign_file() {
    use vep_io::output::fields::FieldOptions;
    use vep_io::output::parquet::{InputOrder, ParquetDirectoryWriter, ParquetOutputFormatter};
    let dir = tempfile::tempdir().unwrap();
    let input = intergenic_vcf(dir.path(), "in.vcf", &[("21", 10), ("22", 5)]);
    let out = dir.path().join("out.parquet");
    let output = parquet_run(&input, &[], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let written = directory_entries(&out);
    assert_eq!(written.len(), 3);
    let digests = |entries: &[String]| -> Vec<(String, u64)> {
        entries
            .iter()
            .map(|e| (e.clone(), std::fs::metadata(out.join(e)).unwrap().len()))
            .collect()
    };
    let before = digests(&written);

    std::fs::write(out.join("notes.txt"), "mine\n").unwrap();
    let refused = parquet_run(&input, &[], &out);
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [parquet] output directory contains files vep-rs did not write"
    );
    let mut expected = written.clone();
    expected.push("notes.txt".to_string());
    expected.sort();
    assert_eq!(directory_entries(&out), expected, "nothing is deleted");
    assert_eq!(digests(&written), before);
    assert_eq!(
        std::fs::read_to_string(out.join("notes.txt")).unwrap(),
        "mine\n"
    );
    std::fs::remove_file(out.join("notes.txt")).unwrap();

    // A Parquet file of the output's own schema whose footer lacks the key.
    let formatter = ParquetOutputFormatter::new(FieldOptions::default(), vec![], vec![]);
    let unkeyed = dir.path().join("unkeyed");
    ParquetDirectoryWriter::create(&unkeyed, &formatter, vec![], InputOrder::Unchecked, None)
        .unwrap()
        .close()
        .unwrap();
    let foreign = out.join("part-00009.parquet");
    std::fs::copy(unkeyed.join("part-00000.parquet"), &foreign).unwrap();
    assert!(!vep_io::output::parquet::is_vep_rs_part(&foreign));
    assert!(vep_io::output::parquet::is_vep_rs_part(
        &out.join("part-00000.parquet")
    ));
    let refused = parquet_run(&input, &[], &out);
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [parquet] output directory contains files vep-rs did not write"
    );
    let mut expected = written.clone();
    expected.push("part-00009.parquet".to_string());
    expected.sort();
    assert_eq!(directory_entries(&out), expected, "nothing is deleted");
    assert_eq!(digests(&written), before);
}

/// The files of a Parquet output directory with their bytes, as paths under
/// it, sorted.
#[cfg(feature = "parquet")]
fn directory_contents(out: &Path) -> Vec<(String, Vec<u8>)> {
    directory_entries(out)
        .into_iter()
        .map(|entry| {
            let bytes = std::fs::read(out.join(&entry)).unwrap();
            (entry, bytes)
        })
        .collect()
}

/// The temporary directories of Parquet runs, and the outputs set aside while
/// one is moved into place, among the entries of `dir`.
#[cfg(feature = "parquet")]
fn temporary_entries(dir: &Path) -> Vec<String> {
    files_in(dir)
        .into_iter()
        .filter(|name| name.contains(".vep-rs-tmp-") || name.contains(".vep-rs-previous-"))
        .collect()
}

/// A two-chromosome Parquet output at `<dir>/out.parquet`, `extra` flags
/// added, with its files and their bytes.
#[cfg(feature = "parquet")]
fn earlier_parquet_output(dir: &Path, extra: &[&str]) -> (PathBuf, Vec<(String, Vec<u8>)>) {
    let sorted = intergenic_vcf(dir, "sorted.vcf", &[("21", 10), ("22", 5)]);
    let out = dir.join("out.parquet");
    let output = parquet_run(&sorted, extra, &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let contents = directory_contents(&out);
    assert_eq!(contents.len(), 3, "two parts and the manifest");
    (out, contents)
}

/// A Parquet run into an earlier output that is refused for its input's order
/// leaves the earlier output as it was, byte for byte, and nothing beside it.
#[cfg(feature = "parquet")]
#[test]
fn a_refused_parquet_replacement_leaves_the_earlier_output_whole() {
    let dir = tempfile::tempdir().unwrap();
    let (out, earlier) = earlier_parquet_output(dir.path(), &[]);
    let unsorted = intergenic_vcf(dir.path(), "unsorted.vcf", &[("21", 100), ("21", 50)]);
    let refused = parquet_run(&unsorted, &[], &out);
    assert_eq!(refused.status.code(), Some(1), "{}", text(&refused.stderr));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [parquet] input is not sorted by chromosome and position at line 2 (21:50 after 21:100); sort the input or pass --no_check_variants_order"
    );
    assert!(
        out.is_dir(),
        "the refused replacement leaves the earlier output in place"
    );
    assert_eq!(directory_contents(&out), earlier);
    assert_eq!(temporary_entries(dir.path()), Vec::<String>::new());
}

/// Runs a Parquet run into `out` that fails once its rows are rendered: one
/// record draws a warning, and the warnings file is Linux's `/dev/full`, a
/// device, which the run does not test before its first record and which
/// fails the write of the warning.
#[cfg(all(feature = "parquet", target_os = "linux"))]
fn parquet_run_failing_on_its_warnings_file(dir: &Path, mode: &str, out: &Path) {
    let warns = intergenic_vcf(dir, "warns.vcf", &[("21", 10), ("Z", 5), ("22", 7)]);
    let warnings = Path::new("/dev/full");
    let failed = parquet_run(
        &warns,
        &[
            "--pipeline_mode",
            mode,
            "--warning_file",
            warnings.to_str().unwrap(),
        ],
        out,
    );
    assert_eq!(failed.status.code(), Some(1), "{}", text(&failed.stderr));
    assert_eq!(
        text(&failed.stderr).trim_end(),
        format!(
            "ERROR: [output] Could not write to warnings file {}",
            warnings.display()
        ),
        "{mode}"
    );
}

/// A Parquet run into an earlier output that fails once rows are rendered
/// leaves the earlier output as it was, in both pipeline modes.
#[cfg(all(feature = "parquet", target_os = "linux"))]
#[test]
fn a_parquet_replacement_that_fails_after_rendering_rows_leaves_the_earlier_output_whole() {
    let dir = tempfile::tempdir().unwrap();
    let (out, earlier) = earlier_parquet_output(dir.path(), &[]);
    for mode in ["batched", "streaming"] {
        parquet_run_failing_on_its_warnings_file(dir.path(), mode, &out);
        assert_eq!(
            directory_contents(&out),
            earlier,
            "{mode}: the failed replacement leaves the earlier output whole"
        );
        assert_eq!(
            temporary_entries(dir.path()),
            Vec::<String>::new(),
            "{mode}"
        );
    }
}

/// A Parquet run into a new path that fails once rows are rendered leaves no
/// directory there and nothing beside it, in both pipeline modes.
#[cfg(all(feature = "parquet", target_os = "linux"))]
#[test]
fn a_failed_parquet_run_into_a_new_path_leaves_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("new.parquet");
    for mode in ["batched", "streaming"] {
        parquet_run_failing_on_its_warnings_file(dir.path(), mode, &out);
        assert!(
            !out.exists(),
            "{mode}: a failed run into a new path leaves no directory"
        );
        assert_eq!(files_in(dir.path()), ["warns.vcf"], "{mode}");
    }
}

/// A Parquet run that fails before its first record, on a statistics report it
/// cannot create, leaves the earlier output as it was.
#[cfg(feature = "parquet")]
#[test]
fn a_parquet_run_that_fails_before_its_first_record_leaves_the_earlier_output_whole() {
    let dir = tempfile::tempdir().unwrap();
    let (out, earlier) = earlier_parquet_output(dir.path(), &[]);
    let stats = dir.path().join("missing").join("stats.html");
    let mut args = without_no_stats(corpus_args_with_input(&dir.path().join("sorted.vcf")));
    args.extend(["--output_format", "parquet", "--quiet"].map(str::to_string));
    let failed = run_with(&args, &["--stats_file", stats.to_str().unwrap()], &out);
    assert_eq!(failed.status.code(), Some(1), "{}", text(&failed.stderr));
    assert_eq!(
        text(&failed.stderr).trim_end(),
        format!(
            "ERROR: [output] Could not write to stats file {}",
            stats.display()
        )
    );
    assert_eq!(
        directory_contents(&out),
        earlier,
        "the failed replacement leaves the earlier output whole"
    );
    assert_eq!(temporary_entries(dir.path()), Vec::<String>::new());
}

/// The release corpus's arguments for a Parquet run that reads VCF records on
/// standard input, `extra` flags added.
#[cfg(feature = "parquet")]
fn parquet_stdin_args(extra: &[&str]) -> Vec<String> {
    let mut args = corpus_args_with_input(Path::new("STDIN"));
    args.extend(["--format", "vcf", "--output_format", "parquet", "--quiet"].map(str::to_string));
    args.extend(extra.iter().map(|flag| flag.to_string()));
    args
}

/// Runs `vep` with `args` writing a Parquet output to `out`, its records on
/// standard input, which stays open until the run has made its temporary
/// directory: all that while what is at `out` is as it was before the run.
/// `while_running` then runs, the input ends, and the finished run is returned.
#[cfg(feature = "parquet")]
fn parquet_run_held_open(
    args: &[String],
    records: &[(&str, u64)],
    out: &Path,
    while_running: impl FnOnce(),
) -> Output {
    let parent = out.parent().unwrap();
    let before = out.is_dir().then(|| directory_contents(out));
    let input = intergenic_vcf(parent, "held.vcf", records);
    let mut child = vep()
        .args(args)
        .args(["-o", out.to_str().unwrap()])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    stdin.write_all(&std::fs::read(&input).unwrap()).unwrap();
    stdin.flush().unwrap();
    let deadline = Instant::now() + Duration::from_secs(120);
    loop {
        assert_eq!(
            out.is_dir().then(|| directory_entries(out)),
            before
                .as_ref()
                .map(|files| files.iter().map(|(name, _)| name.clone()).collect()),
            "the output path is as it was while the replacement runs"
        );
        if !temporary_entries(parent).is_empty() {
            break;
        }
        assert!(
            child.try_wait().unwrap().is_none(),
            "the run ended before its input did"
        );
        assert!(
            Instant::now() < deadline,
            "the run made no temporary directory"
        );
        std::thread::sleep(Duration::from_millis(5));
    }
    assert_eq!(
        out.is_dir().then(|| directory_contents(out)),
        before,
        "the earlier output is whole while the replacement runs"
    );
    while_running();
    drop(stdin);
    child.wait_with_output().unwrap()
}

/// A `--parquet_partition chrom` output replaced by a run over fewer
/// chromosomes holds the new run alone, the directory of the chromosome the
/// new run lacks gone; the earlier output stays whole until the new one is
/// complete.
#[cfg(feature = "parquet")]
#[test]
fn a_partitioned_parquet_output_is_replaced_by_one_with_fewer_chromosomes() {
    use vep_io::output::parquet::Manifest;
    let dir = tempfile::tempdir().unwrap();
    let (out, earlier) = earlier_parquet_output(dir.path(), &["--parquet_partition", "chrom"]);
    assert_eq!(
        earlier.iter().map(|(name, _)| name).collect::<Vec<_>>(),
        [
            "_vep_rs_manifest.json",
            "chrom=21/part-00000.parquet",
            "chrom=22/part-00001.parquet"
        ]
    );
    let earlier_run = Manifest::read(&out).unwrap().run_id;
    let output = parquet_run_held_open(
        &parquet_stdin_args(&["--parquet_partition", "chrom"]),
        &[("21", 10), ("21", 20)],
        &out,
        || {},
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        directory_entries(&out),
        ["_vep_rs_manifest.json", "chrom=21/part-00000.parquet"]
    );
    assert!(!out.join("chrom=22").exists());
    let manifest = Manifest::read(&out).unwrap();
    assert_eq!(manifest.rows, 2);
    assert_ne!(manifest.run_id, earlier_run);
    assert_eq!(temporary_entries(dir.path()), Vec::<String>::new());
}

/// A partitioned output replaced by a plain one, and a plain one by a
/// partitioned one, holds the new run's layout alone.
#[cfg(feature = "parquet")]
#[test]
fn a_parquet_output_is_replaced_across_the_partitioned_and_the_plain_layout() {
    let dir = tempfile::tempdir().unwrap();
    let (out, _) = earlier_parquet_output(dir.path(), &["--parquet_partition", "chrom"]);
    let records = [("21", 10), ("22", 5)];
    let output = parquet_run_held_open(&parquet_stdin_args(&[]), &records, &out, || {});
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        directory_entries(&out),
        [
            "_vep_rs_manifest.json",
            "part-00000.parquet",
            "part-00001.parquet"
        ]
    );
    let partitioned = parquet_stdin_args(&["--parquet_partition", "chrom"]);
    let output = parquet_run_held_open(&partitioned, &records, &out, || {});
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        directory_entries(&out),
        [
            "_vep_rs_manifest.json",
            "chrom=21/part-00000.parquet",
            "chrom=22/part-00001.parquet"
        ]
    );
    assert_eq!(
        vep_io::output::parquet::sort_keys(&out).unwrap(),
        [(0, 10), (1, 5)]
    );
    assert_eq!(temporary_entries(dir.path()), Vec::<String>::new());
}

/// A populated output replaced by a run with no rows holds one empty part and
/// a manifest of 0 rows, every directory of the earlier output gone.
#[cfg(feature = "parquet")]
#[test]
fn a_populated_parquet_output_is_replaced_by_an_empty_run() {
    use vep_io::output::parquet::Manifest;
    let dir = tempfile::tempdir().unwrap();
    let (out, _) = earlier_parquet_output(dir.path(), &["--parquet_partition", "chrom"]);
    let partitioned = parquet_stdin_args(&["--parquet_partition", "chrom"]);
    let output = parquet_run_held_open(&partitioned, &[("Z", 5)], &out, || {});
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        directory_entries(&out),
        ["_vep_rs_manifest.json", "part-00000.parquet"]
    );
    let manifest = Manifest::read(&out).unwrap();
    assert_eq!(manifest.rows, 0);
    assert_eq!(
        manifest
            .parts
            .iter()
            .map(|p| (p.path.as_str(), p.rows))
            .collect::<Vec<_>>(),
        [("part-00000.parquet", 0)]
    );
    assert_eq!(
        vep_io::output::parquet::sort_keys(&out).unwrap(),
        Vec::<(i32, i64)>::new()
    );
    assert_eq!(temporary_entries(dir.path()), Vec::<String>::new());
}

/// The output path is checked again once the run has succeeded. A file vep-rs
/// did not write that arrives in the earlier output while the run reads its
/// input refuses the replacement, and a directory that arrives at a new path
/// refuses a run without `--force_overwrite`; either way everything stays in
/// place and the run's own directory is removed.
#[cfg(feature = "parquet")]
#[test]
fn a_parquet_output_path_is_checked_again_before_the_run_replaces_it() {
    let dir = tempfile::tempdir().unwrap();
    let (out, earlier) = earlier_parquet_output(dir.path(), &[]);
    let notes = out.join("notes.txt");
    let records = [("21", 10)];
    let refused = parquet_run_held_open(&parquet_stdin_args(&[]), &records, &out, || {
        std::fs::write(&notes, "mine\n").unwrap();
    });
    assert_eq!(refused.status.code(), Some(1), "{}", text(&refused.stderr));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        "ERROR: [parquet] output directory contains files vep-rs did not write"
    );
    assert_eq!(std::fs::read_to_string(&notes).unwrap(), "mine\n");
    std::fs::remove_file(&notes).unwrap();
    assert_eq!(directory_contents(&out), earlier);
    assert_eq!(temporary_entries(dir.path()), Vec::<String>::new());

    let new_path = dir.path().join("new.parquet");
    let unforced = without_force_overwrite(parquet_stdin_args(&[]));
    let refused = parquet_run_held_open(&unforced, &records, &new_path, || {
        std::fs::create_dir(&new_path).unwrap();
    });
    assert_eq!(refused.status.code(), Some(1), "{}", text(&refused.stderr));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        format!(
            "ERROR: [output] Output file {} already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite",
            new_path.display()
        )
    );
    assert_eq!(std::fs::read_dir(&new_path).unwrap().count(), 0);
    assert_eq!(temporary_entries(dir.path()), Vec::<String>::new());
}

/// A Parquet output directory keeps its permission bits when a run replaces
/// it.
#[cfg(all(feature = "parquet", unix))]
#[test]
fn a_parquet_replacement_keeps_the_output_directorys_permission_bits() {
    use std::os::unix::fs::PermissionsExt;
    let dir = tempfile::tempdir().unwrap();
    let (out, _) = earlier_parquet_output(dir.path(), &[]);
    std::fs::set_permissions(&out, std::fs::Permissions::from_mode(0o700)).unwrap();
    let output = parquet_run(&dir.path().join("sorted.vcf"), &[], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        std::fs::metadata(&out).unwrap().permissions().mode() & 0o7777,
        0o700,
        "the replaced directory keeps its permission bits"
    );
    assert_eq!(directory_entries(&out).len(), 3);
    assert_eq!(temporary_entries(dir.path()), Vec::<String>::new());
}

/// A failed Parquet run that cannot remove its temporary directory, here
/// because the output's directory lost its write permission while the run
/// read its input, names that directory in its error line.
#[cfg(all(feature = "parquet", unix))]
#[test]
fn a_failed_parquet_run_names_a_temporary_directory_it_cannot_remove() {
    use std::os::unix::fs::PermissionsExt;
    let dir = tempfile::tempdir().unwrap();
    let parent = dir.path().join("results");
    std::fs::create_dir(&parent).unwrap();
    let out = parent.join("out.parquet");
    let read_only = std::cell::Cell::new(false);
    let failed = parquet_run_held_open(&parquet_stdin_args(&[]), &[("21", 10)], &out, || {
        std::fs::set_permissions(&parent, std::fs::Permissions::from_mode(0o555)).unwrap();
        read_only.set(std::fs::write(parent.join("probe"), "").is_err());
    });
    std::fs::set_permissions(&parent, std::fs::Permissions::from_mode(0o755)).unwrap();
    if !read_only.get() {
        // This user's writes ignore the mode, so nothing fails.
        return;
    }
    let left = temporary_entries(&parent);
    assert_eq!(left.len(), 1, "{left:?}");
    assert_eq!(failed.status.code(), Some(1), "{}", text(&failed.stderr));
    assert_eq!(
        text(&failed.stderr).trim_end(),
        format!(
            "ERROR: [output] Could not write to output file {}; the temporary directory {} could not be removed",
            out.display(),
            parent.join(&left[0]).display()
        )
    );
    assert!(!out.exists());
}

/// A symbolic link at the output path, or inside an earlier output, refuses
/// the run by name before it starts, and what the link points at keeps its
/// files; an output path under a linked directory is written as any other.
#[cfg(all(feature = "parquet", unix))]
#[test]
fn a_symbolic_link_at_or_inside_a_parquet_output_refuses_the_run() {
    use std::os::unix::fs::symlink;
    let dir = tempfile::tempdir().unwrap();
    let (target, earlier) = earlier_parquet_output(dir.path(), &[]);
    let input = dir.path().join("sorted.vcf");

    let link = dir.path().join("link.parquet");
    symlink(&target, &link).unwrap();
    let refused = parquet_run(&input, &[], &link);
    assert_eq!(refused.status.code(), Some(1), "{}", text(&refused.stderr));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        format!(
            "ERROR: [parquet] output path {} is a symbolic link",
            link.display()
        )
    );
    assert!(link.symlink_metadata().unwrap().file_type().is_symlink());
    assert_eq!(
        directory_contents(&target),
        earlier,
        "the link's target keeps its files"
    );

    let dangling = dir.path().join("dangling.parquet");
    symlink(dir.path().join("nowhere"), &dangling).unwrap();
    let mut unforced = without_force_overwrite(corpus_args_with_input(&input));
    unforced.extend(["--output_format", "parquet", "--quiet"].map(str::to_string));
    for refused in [
        parquet_run(&input, &[], &dangling),
        run_with(&unforced, &[], &dangling),
    ] {
        assert_eq!(refused.status.code(), Some(1), "{}", text(&refused.stderr));
        assert_eq!(
            text(&refused.stderr).trim_end(),
            format!(
                "ERROR: [parquet] output path {} is a symbolic link",
                dangling.display()
            )
        );
        assert!(!dir.path().join("nowhere").exists());
    }

    let elsewhere = dir.path().join("elsewhere");
    std::fs::create_dir(&elsewhere).unwrap();
    std::fs::copy(
        target.join("part-00000.parquet"),
        elsewhere.join("part-00002.parquet"),
    )
    .unwrap();
    let inside = target.join("chrom=22");
    symlink(&elsewhere, &inside).unwrap();
    let refused = parquet_run(&input, &[], &target);
    assert_eq!(refused.status.code(), Some(1), "{}", text(&refused.stderr));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        format!(
            "ERROR: [parquet] output directory entry {} is a symbolic link",
            inside.display()
        )
    );
    assert!(
        elsewhere.join("part-00002.parquet").exists(),
        "the link's target keeps its files"
    );
    std::fs::remove_file(&inside).unwrap();
    assert_eq!(directory_contents(&target), earlier);

    let linked_parent = dir.path().join("linked");
    symlink(&elsewhere, &linked_parent).unwrap();
    let output = parquet_run(&input, &[], &linked_parent.join("out.parquet"));
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(directory_entries(&elsewhere.join("out.parquet")).len(), 3);
    assert_eq!(temporary_entries(dir.path()), Vec::<String>::new());
    assert_eq!(temporary_entries(&elsewhere), Vec::<String>::new());
}

/// `--no_check_variants_order` writes the unsorted input in input order: no
/// row group declares a sort, the footer's `sort` is `none`, `contig_index` is
/// still set on every row, and the run warns once, into the warnings file and
/// on standard error. An identifier, HGVS, SPDI or region input is written the
/// same way without the flag, since its records carry no position to sort by.
#[cfg(feature = "parquet")]
#[test]
fn no_check_variants_order_relaxes_the_parquet_layout() {
    let dir = tempfile::tempdir().unwrap();
    let input = intergenic_vcf(dir.path(), "in.vcf", &[("21", 100), ("21", 50), ("22", 7)]);
    let out = dir.path().join("out.parquet");
    let mut args = corpus_args_with_input(&input);
    args.extend(["--output_format", "parquet"].map(str::to_string));
    let output = run_with(&args, &["--no_check_variants_order"], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let warning = "WARNING: input order not checked: Parquet row groups are not sorted and parts are not rolled at contig boundaries";
    assert_eq!(warning_lines(&text(&output.stderr)), [warning]);
    let warnings = std::fs::read_to_string(dir.path().join("out.parquet_warnings.txt")).unwrap();
    assert_eq!(warnings.lines().collect::<Vec<_>>(), [warning]);
    let keys = vep_io::output::parquet::sort_keys(&out).unwrap();
    assert_eq!(
        keys,
        [(0, 100), (0, 50), (1, 7)],
        "rows in input order, ranked"
    );
    let (kv, sorted) = parquet_footer(&out);
    assert_eq!(kv["sort"], "none");
    assert_eq!(sorted, [false]);

    let region = dir.path().join("in.region");
    std::fs::write(
        &region,
        "21:30165811-30165811:1/G\n21:30165801-30165801:1/G\n",
    )
    .unwrap();
    let fasta = decompressed_fasta(dir.path());
    let out = dir.path().join("region.parquet");
    let mut args = corpus_args_with_input(&region);
    args.extend(
        [
            "--output_format",
            "parquet",
            "--format",
            "region",
            "--quiet",
        ]
        .map(str::to_string),
    );
    let output = run_with(&args, &["--fasta", fasta.to_str().unwrap()], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(
        text(&output.stderr).is_empty(),
        "quiet: {}",
        text(&output.stderr)
    );
    let warnings = std::fs::read_to_string(dir.path().join("region.parquet_warnings.txt")).unwrap();
    assert_eq!(warnings.lines().collect::<Vec<_>>(), [warning]);
    let keys = vep_io::output::parquet::sort_keys(&out).unwrap();
    assert_eq!(
        keys.iter().map(|k| k.1).collect::<Vec<_>>(),
        [30165811, 30165801],
        "input order"
    );
    assert_eq!(parquet_footer(&out).0["sort"], "none");
}

/// Under `--parquet_partition chrom` a part closes at every change of spelling
/// whether or not the order is checked, so a run that skips the check warns
/// with the partitioned layout's own text and never that its parts are not
/// rolled at contig boundaries: the unsorted input is written in input order
/// as one part per run of a spelling, each in its spelling's directory.
#[cfg(feature = "parquet")]
#[test]
fn the_unchecked_order_warning_is_true_of_a_partitioned_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = intergenic_vcf(dir.path(), "in.vcf", &[("21", 100), ("21", 50), ("22", 7)]);
    let out = dir.path().join("out.parquet");
    let mut args = corpus_args_with_input(&input);
    args.extend(["--output_format", "parquet"].map(str::to_string));
    let output = run_with(
        &args,
        &["--no_check_variants_order", "--parquet_partition", "chrom"],
        &out,
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let warning = "WARNING: input order not checked: Parquet row groups are not sorted and a part closes at every change of chromosome spelling";
    let stderr = text(&output.stderr);
    assert_eq!(
        stderr.lines().filter(|line| *line == warning).count(),
        1,
        "{stderr}"
    );
    assert!(
        !stderr.contains("parts are not rolled at contig boundaries"),
        "{stderr}"
    );
    let warnings = std::fs::read_to_string(dir.path().join("out.parquet_warnings.txt")).unwrap();
    assert_eq!(
        warnings.lines().filter(|line| *line == warning).count(),
        1,
        "{warnings}"
    );
    assert!(
        !warnings.contains("parts are not rolled at contig boundaries"),
        "{warnings}"
    );
    assert_eq!(
        directory_entries(&out),
        [
            "_vep_rs_manifest.json",
            "chrom=21/part-00000.parquet",
            "chrom=22/part-00001.parquet",
        ]
    );
    let keys = vep_io::output::parquet::sort_keys(&out).unwrap();
    assert_eq!(keys, [(0, 100), (0, 50), (1, 7)], "rows in input order");
    assert_eq!(parquet_footer(&out).0["sort"], "none");

    let returning = intergenic_vcf(
        dir.path(),
        "returning.vcf",
        &[("21", 100), ("22", 7), ("21", 50)],
    );
    let out = dir.path().join("returning.parquet");
    let mut args = corpus_args_with_input(&returning);
    args.extend(["--output_format", "parquet", "--quiet"].map(str::to_string));
    let output = run_with(
        &args,
        &["--no_check_variants_order", "--parquet_partition", "chrom"],
        &out,
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        directory_entries(&out),
        [
            "_vep_rs_manifest.json",
            "chrom=21/part-00000.parquet",
            "chrom=21/part-00002.parquet",
            "chrom=22/part-00001.parquet",
        ]
    );
}

/// Every warning of a run is written to `<output_file>_warnings.txt` beside the
/// output (BaseVEP.pm:890), as `WARNING: line N skipped (<line>): <description>`
/// for a skipped record (`:840-863`), and echoed on standard error; the file is
/// Ensembl VEP's for the run-control corpus. A run without a warning leaves no
/// file (`:881`: the handle is opened at the first warning); the attributes
/// corpus is such a run, where the release corpus's tandem duplication, read
/// without a FASTA, draws Ensembl VEP's warning (Parser/VCF.pm:521-524).
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

    let clean =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-attributes");
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

/// A cache without variation data, asked for co-located variants, warns once
/// that the co-location columns stay empty, naming the flag that asks for
/// them and the caches that hold the data, where Ensembl VEP annotates
/// silently (CacheDir.pm:179 attaches no variation source). The line is a
/// warning of the run like any other (BaseVEP.pm:802-822): written to the
/// warnings file and echoed on standard error, and under `--quiet` written to
/// the file alone (`:821` drops the standard-error copy). The release corpus's
/// cache holds transcripts and no variation files.
#[test]
fn cache_without_variation_data_warns_into_the_warnings_file() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let args = corpus_args();
    let cache_flag = args.iter().position(|a| a == "--json_cache").unwrap();
    let cache = &args[cache_flag + 1];
    let warning = format!(
        "WARNING: [cache] no variation data under {cache}: Existing_variation and the frequency \
         columns stay empty under --everything; a cache from vep-cache-builder --auto or \
         --include-variations, or a conversion of Ensembl VEP's cache, holds it \
         (docs/cache-variation-data.md)"
    );
    let output = run_with(&args, &["--everything"], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    assert!(
        warnings.lines().any(|l| l == warning),
        "the warnings file lacks the cache warning:\n{warnings}"
    );
    let stderr = text(&output.stderr);
    assert_eq!(
        warning_lines(&stderr)
            .iter()
            .filter(|l| **l == warning)
            .count(),
        1,
        "{stderr}"
    );
    let mut echoed = warning_lines(&stderr);
    echoed.sort_unstable();
    let mut written: Vec<&str> = warnings.lines().collect();
    written.sort_unstable();
    assert_eq!(
        echoed, written,
        "every warning of the run is the warnings file's"
    );

    // A second output: a warnings file is created at a run's first warning, so
    // a reused path could still hold the first run's file.
    let quiet_out = dir.path().join("quiet.txt");
    let quiet = run_with(&args, &["--everything", "--quiet"], &quiet_out);
    assert_eq!(quiet.status.code(), Some(0), "{}", text(&quiet.stderr));
    assert!(
        warning_lines(&text(&quiet.stderr)).is_empty(),
        "{}",
        text(&quiet.stderr)
    );
    let warnings = std::fs::read_to_string(dir.path().join("quiet.txt_warnings.txt")).unwrap();
    assert!(
        warnings.lines().any(|l| l == warning),
        "--quiet drops the standard-error copy alone:\n{warnings}"
    );
}

/// `--warning_file` names the file; the literal `stderr` in any case
/// (BaseVEP.pm:885 `/^stderr$/i`) keeps the warnings on standard error and
/// writes no file: each warning twice, as `print $fh` (:819) then `warn` (:821)
/// write it, and once under `--quiet`, which drops the `warn` copy alone.
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
        let copies = if quiet { 1 } else { 2 };
        let expected_text = expected_warnings();
        let expected: Vec<&str> = expected_text
            .lines()
            .flat_map(|line| std::iter::repeat_n(line, copies))
            .collect();
        assert_eq!(
            warning_lines(&text(&output.stderr)),
            expected,
            "{spelling} quiet={quiet}: {copies} copies of each warning on standard error"
        );
        assert!(
            !dir.path().join(spelling).exists(),
            "{spelling}: a file was written"
        );
        assert!(!dir.path().join("out.txt_warnings.txt").exists());
    }
}

/// A `--warning_file` or `--skipped_variants_file` value is a file name as
/// given (`target_for` in `warnings.rs` builds a `PathBuf` from every value but
/// `stderr`, `reserve_skipped_variants_file` from every value): `-`, ` sp.txt`
/// and `>pw.txt` are created under those names in the working directory, the
/// file whose name follows the `>` is left alone and nothing reaches standard
/// output. Ensembl VEP opens both with Perl's two-argument `open` (BaseVEP.pm:894
/// `$fh->open(">".$file)`, BaseRunner.pm:344 `$fh->open(">$file_name")`), which
/// sends `-` to standard output, writes `' sp.txt'` to `sp.txt` and appends
/// `'>pw.txt'` to `pw.txt`, the readings `docs/cli/compatibility.md` states.
#[test]
fn warning_file_and_skipped_file_values_are_file_names_as_given() {
    for flag in ["--warning_file", "--skipped_variants_file"] {
        for value in ["-", " sp.txt", ">pw.txt"] {
            let dir = tempfile::tempdir().unwrap();
            let out = dir.path().join("out.txt");
            std::fs::write(dir.path().join("pw.txt"), "kept\n").unwrap();
            let output = run_control(
                &["-o", out.to_str().unwrap(), flag, value, "--quiet"],
                dir.path(),
            );
            let case = format!("{flag} {value:?}");
            assert_eq!(
                output.status.code(),
                Some(0),
                "{case}: {}",
                text(&output.stderr)
            );
            assert_eq!(text(&output.stdout), "", "{case}: standard output");
            let written = std::fs::read_to_string(dir.path().join(value))
                .unwrap_or_else(|e| panic!("{case}: no file named as given: {e}"));
            if flag == "--warning_file" {
                assert_eq!(written, expected_warnings(), "{case}");
            } else {
                assert!(
                    written.starts_with("[VEP skipped the following variants from "),
                    "{case}: {written}"
                );
            }
            assert_eq!(
                std::fs::read_to_string(dir.path().join("pw.txt")).unwrap(),
                "kept\n",
                "{case}: the file named after the `>` is left alone"
            );
            assert!(
                !dir.path().join("sp.txt").exists(),
                "{case}: no file under the trimmed name"
            );
        }
    }
}

/// A `--warning_file` whose value is false in Perl, the empty string or `0`,
/// is no file: `$file ||= ...'_warnings.txt'` (BaseVEP.pm:890) falls back to
/// the default beside the output, and nothing named `0` is written.
#[test]
fn warning_file_false_value_is_the_default_file() {
    for value in ["", "0"] {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.txt");
        let output = run_control(
            &["-o", out.to_str().unwrap(), "--warning_file", value],
            dir.path(),
        );
        assert_eq!(
            output.status.code(),
            Some(0),
            "{value:?}: {}",
            text(&output.stderr)
        );
        assert_eq!(
            std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap(),
            expected_warnings(),
            "{value:?}"
        );
        assert!(
            !dir.path().join("0").exists(),
            "{value:?}: a file `0` was written"
        );
    }
}

/// A `--skipped_variants_file` value that is false in Perl is a file name here
/// (`reserve_skipped_variants_file` in `warnings.rs` reserves every value): `0`
/// is written as a file named `0` and the empty value is refused before the
/// first record is read, with no output written, with
/// `Could not write to skipped variants file` (exit 1), where Ensembl VEP's
/// `$self->dump_skipped_variants if $self->param('skipped_variants_file')`
/// (Runner.pm:176) writes no file for either and exits 0, the reading
/// `docs/cli/compatibility.md` states.
#[test]
fn skipped_variants_file_false_value_is_a_file_name() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let zero = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--skipped_variants_file",
            "0",
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(zero.status.code(), Some(0), "{}", text(&zero.stderr));
    assert!(std::fs::read_to_string(dir.path().join("0"))
        .unwrap()
        .starts_with("[VEP skipped the following variants from "));
    let written = std::fs::read(&out).unwrap();

    let empty = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--skipped_variants_file",
            "",
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(empty.status.code(), Some(1));
    assert_eq!(
        text(&empty.stderr),
        "ERROR: [output] Could not write to skipped variants file \n",
        "the empty name follows the message's trailing space"
    );
    assert!(
        std::fs::read(&out).unwrap() == written,
        "the refused run wrote over the earlier output"
    );
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
/// reason are two warnings, each quoting its line as Ensembl VEP quotes its
/// parser's `current_block` (BaseVEP.pm:846-851). The region and HGVS runs name
/// a FASTA, which both formats need (Parser/Region.pm:115; the HGVS notation is
/// resolved against it), so a chromosome the FASTA lacks fails an HGVS notation
/// before the cache is consulted, with the parser's own warning, as the slice
/// fetch fails it in Ensembl VEP.
#[test]
fn non_vcf_skips_carry_the_record_number() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = decompressed_fasta(dir.path());
    let not_in_cache = |line: &str| {
        format!("({line}): Chromosome Z not found in annotation sources or synonyms; chromosome Z does not overlap any features")
    };
    let not_in_fasta = |notation: &str| {
        format!("({notation}): Unable to parse HGVS notation '{notation}': The reference FASTA has no sequence named Z")
    };
    for (format, records, warnings) in [
        (
            "ensembl",
            "21\t30165801\t30165801\tA/G\t1\n\
             Z\t1000\t1000\tA/G\t1\n\
             Z\t2000\t2000\tA/G\t1\n\
             21\t30165811\t30165811\tA/G\t1\n",
            [
                not_in_cache("Z 1000 1000 A/G 1"),
                not_in_cache("Z 2000 2000 A/G 1"),
            ],
        ),
        (
            "region",
            "21:30165801-30165801:1/G\n\
             Z:1000-1000:1/G\n\
             Z:2000-2000:1/G\n\
             21:30165811-30165811:1/G\n",
            [
                not_in_cache("Z:1000-1000:1/G"),
                not_in_cache("Z:2000-2000:1/G"),
            ],
        ),
        (
            "hgvs",
            "21:g.30165801A>G\n\
             Z:g.1000A>G\n\
             Z:g.2000A>G\n\
             21:g.30165811T>G\n",
            [not_in_fasta("Z:g.1000A>G"), not_in_fasta("Z:g.2000A>G")],
        ),
    ] {
        let input = dir.path().join(format!("in.{format}"));
        std::fs::write(&input, records).unwrap();
        let out = dir.path().join(format!("out.{format}"));
        // The corpus lifts the structural-variant limit; a genomic HGVS notation
        // under the lifted limit is skipped, as Ensembl VEP skips it.
        let mut args = corpus_args_with_input(&input);
        if format == "hgvs" {
            let limit = args.iter().position(|a| a == "--max_sv_size").unwrap();
            args.drain(limit..limit + 2);
        }
        let mut command = vep();
        command
            .args(args)
            .args(["--format", format, "-o", out.to_str().unwrap(), "--quiet"]);
        if format != "ensembl" {
            command.arg("--fasta").arg(&fasta);
        }
        let output = command.stdin(Stdio::null()).output().unwrap();
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
                format!("WARNING: line 2 skipped {}", warnings[0]),
                format!("WARNING: line 3 skipped {}", warnings[1]),
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
/// a `-` and is annotated, while `X/Y` is skipped. The format is named: a first
/// data line whose REF is `XX` is of no shape `check_format` knows.
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
        .args(["-o", out.to_str().unwrap(), "--quiet", "--format", "vcf"])
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

/// A structural record whose INFO/END is no number to Perl is skipped with two
/// lines: `deletion looks incomplete` (Parser/VCF.pm:555-558: `abc` numifies to 0,
/// below the start) and then `validate_vf`'s coordinate test over the start Ensembl VEP
/// reads, `POS + 1`, and the END text as written (Parser.pm:475-477, the text kept
/// by BaseVCF4.pm:307-308). Both lines are Ensembl VEP 116.3's for the record,
/// which writes no row for it and reaches no cache.
#[test]
fn end_that_is_no_number_is_refused_with_the_coordinate_text() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t100\t.\tN\t<DEL>\t.\t.\tEND=abc\n",
    )
    .unwrap();
    let out = dir.path().join("out.txt");
    let cache = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-vcf-edge-cases/json_cache");
    let output = vep()
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache)
        .args(["--assembly", "GRCh38", "--force_overwrite", "--no_stats"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let expected = [
        "WARNING: line 1 skipped (21 100 . N <DEL> . . END=abc): deletion looks incomplete",
        "WARNING: line 1 skipped (21 100 . N <DEL> . . END=abc): Invalid start '101' or end 'abc' coordinate",
    ];
    assert_eq!(warning_lines(&text(&output.stderr)), expected);
    assert_eq!(
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt"))
            .unwrap()
            .lines()
            .collect::<Vec<_>>(),
        expected
    );
    assert!(data_lines(&std::fs::read_to_string(&out).unwrap()).is_empty());
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

/// Runs the run-control cache over `records` (VCF, with `columns` after INFO in
/// the header line and the first record a sequence variant the format is
/// detected from, so `records` starts at line 2) in `dir`, `extra` flags last,
/// and returns the process output with the warnings file's text.
fn run_control_over(records: &str, columns: &str, extra: &[&str], dir: &Path) -> (Output, String) {
    let input = dir.join("records.vcf");
    let samples = columns.matches('\t').count().saturating_sub(1);
    let first = format!(
        "21\t30282169\tstop_ins4_twostops_rev\tC\tCTATT\t.\t.\t.{}{}\n",
        if samples > 0 { "\tGT" } else { "" },
        "\t0/1".repeat(samples)
    );
    std::fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO{columns}\n{first}{records}"),
    )
    .unwrap();
    let out = dir.join("out.txt");
    let output = vep()
        .args(with_input(args_for(&run_control_dir()), &input))
        .args(["-o", out.to_str().unwrap()])
        .args(extra)
        .current_dir(dir)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let warnings = std::fs::read_to_string(dir.join("out.txt_warnings.txt")).unwrap_or_default();
    (output, warnings)
}

/// Under `--individual all`, each sample's copy of a refused record is validated
/// apart (VCF.pm:199 over the copies of :595-666), so one record is two skips:
/// the identical warning text is written once, as `_seen_warnings` keys the whole
/// text (BaseVEP.pm:810-813), on standard error and in the warnings file alike,
/// while `skipped_variants` records every copy (:858), so the skipped-variants
/// file lists the line twice.
#[test]
fn individual_copies_of_a_refused_record_warn_once_and_are_listed_twice() {
    let dir = tempfile::tempdir().unwrap();
    let skipped = dir.path().join("skipped.txt");
    let (output, warnings) = run_control_over(
        "21\t30165901\tinvalid_alleles\tX\tY\t.\t.\t.\tGT\t0/1\t0/1\n",
        "\tFORMAT\tS1\tS2",
        &[
            "--individual",
            "all",
            "--skipped_variants_file",
            skipped.to_str().unwrap(),
        ],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let line = "WARNING: line 2 skipped (21 30165901 invalid_alleles X Y . . . GT 0/1 0/1): Invalid allele string X/Y or possible parsing error";
    assert_eq!(warning_lines(&text(&output.stderr)), [line]);
    assert_eq!(warnings, format!("{line}\n"));
    let listed = std::fs::read_to_string(&skipped).unwrap();
    let rows: Vec<&str> = listed.lines().skip(1).collect();
    assert_eq!(rows.len(), 2, "{listed}");
    for row in rows {
        assert_eq!(
            row,
            "Line 2    \t21 30165901 invalid_alleles X Y . . . GT\tInvalid allele string X/Y or possible parsing error"
        );
    }
}

/// The line a skip warning quotes is cut by bytes, as `length` and `substr`
/// count them on the undecoded line Ensembl VEP reads (BaseVEP.pm:848-850,
/// `$maxChar = 50`): a 62-byte record whose identifier is twenty two-byte
/// characters is quoted through byte 46, the first seventeen of them, then
/// `...`, the bytes Ensembl VEP 116.3 writes for this line.
#[test]
fn skip_warning_quotes_the_line_cut_by_bytes() {
    let dir = tempfile::tempdir().unwrap();
    let id = "é".repeat(20);
    let record = format!("21\t30165901\t{id}\tX\tY\t.\t.\t.\n");
    assert_eq!(record.len(), 63, "62 bytes and the newline");
    let (output, warnings) = run_control_over(&record, "", &[], dir.path());
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let line = format!(
        "WARNING: line 2 skipped (21 30165901 {}...): Invalid allele string X/Y or possible parsing error",
        "é".repeat(17)
    );
    assert_eq!(warning_lines(&text(&output.stderr)), [line.as_str()]);
    assert_eq!(warnings, format!("{line}\n"));
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

/// `--safe` skips the incompatible-pair check (Config.pm:961), so `--quiet
/// --verbose --safe` runs, and `--quiet` wins: `status_msg` returns on it before
/// printing (BaseVEP.pm:782), so standard error stays empty, the status lines
/// `--verbose` would add included, while the warnings file is still written.
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
    assert_eq!(text(&output.stderr), "", "--quiet silences --verbose");
    assert_eq!(
        std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap(),
        expected_warnings()
    );
}

/// The pick-summary corpus: multi-allelic records over several transcripts, the
/// fixture of the pick family's flag sets.
fn pick_summary_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-pick-summary")
}

/// Runs the binary on the pick-summary corpus with `extra` after the corpus
/// arguments, writing `out`.
fn pick_summary(extra: &[&str], out: &Path) -> Output {
    run_with(&args_for(&pick_summary_dir()), extra, out)
}

/// The rows of a default-format output for one record, each as its first seven
/// columns (the identifier through the consequence).
fn record_rows(output: &str, id: &str) -> Vec<String> {
    let prefix = format!("{id}\t");
    data_lines(output)
        .into_iter()
        .filter(|l| l.starts_with(&prefix))
        .map(|l| l.split('\t').take(7).collect::<Vec<_>>().join("\t"))
        .collect()
}

/// A repeated `--pick_order` keeps its last value, as Getopt::Long keeps the last
/// value of the scalar `pick_order=s` (Config.pm:117, not in `@ALLOW_MULTIPLE`,
/// :340-346): under `--pick --pick_order length --pick_order rank` Ensembl VEP
/// 116.3 picks `C ENST00000339775`, the row `rank` alone picks, for the record
/// `21 17604928 G A,C`; `length` alone would pick another transcript.
#[test]
fn repeated_pick_order_keeps_the_last_list() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let repeated = pick_summary(
        &["--pick", "--pick_order", "length", "--pick_order", "rank"],
        &out,
    );
    assert_eq!(
        repeated.status.code(),
        Some(0),
        "{}",
        text(&repeated.stderr)
    );
    let rows = record_rows(&std::fs::read_to_string(&out).unwrap(), "ps_17604928");
    assert_eq!(
        rows,
        ["ps_17604928\t21:17604928\tC\tENSG00000154640\tENST00000339775\tTranscript\tmissense_variant"]
    );
    let length_alone = pick_summary(&["--pick", "--pick_order", "length"], &out);
    assert_eq!(length_alone.status.code(), Some(0));
    assert_ne!(
        record_rows(&std::fs::read_to_string(&out).unwrap(), "ps_17604928"),
        rows,
        "the last value decides the pick"
    );
}

/// `--pick_order ""` is the empty order: Config.pm:761 splits the value to
/// nothing, `_is_flag_active` reads the empty list as inactive (:902-906) so the
/// value check never sees it, and `pick_worst_VariationFeatureOverlapAllele`
/// compares no category and returns the first candidate (OutputFactory.pm:808).
/// Ensembl VEP 116.3 picks `A ENST00000339775 synonymous_variant`, its first row,
/// for the record `21 17604928 G A,C`; the default order picks another row.
#[test]
fn an_empty_pick_order_picks_the_first_row() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let empty = pick_summary(&["--pick", "--pick_order", ""], &out);
    assert_eq!(empty.status.code(), Some(0), "{}", text(&empty.stderr));
    assert!(
        !text(&empty.stderr).contains("WARNING"),
        "{}",
        text(&empty.stderr)
    );
    let rows = record_rows(&std::fs::read_to_string(&out).unwrap(), "ps_17604928");
    assert_eq!(
        rows,
        ["ps_17604928\t21:17604928\tA\tENSG00000154640\tENST00000339775\tTranscript\tsynonymous_variant"]
    );
    let default_order = pick_summary(&["--pick"], &out);
    assert_eq!(default_order.status.code(), Some(0));
    assert_ne!(
        record_rows(&std::fs::read_to_string(&out).unwrap(), "ps_17604928"),
        rows,
        "the default order compares its criteria"
    );
}

/// `--safe` skips the check that refuses `--summary` with `--most_severe`
/// (Config.pm:961), and the pair then lists every term: `summary_only` routes on
/// either flag (OutputFactory.pm:387) and tests `summary` first (`:436`), the
/// most-severe term alone being its `else` (`:446-447`). The data rows are the
/// `--summary` run's.
#[test]
fn safe_lets_summary_win_over_most_severe() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let written_rows = || -> Vec<String> {
        data_lines(&std::fs::read_to_string(&out).unwrap())
            .into_iter()
            .map(str::to_string)
            .collect()
    };
    let pair = pick_summary(&["--safe", "--summary", "--most_severe"], &out);
    assert_eq!(pair.status.code(), Some(0), "{}", text(&pair.stderr));
    let pair_rows = written_rows();
    let summary = pick_summary(&["--summary"], &out);
    assert_eq!(summary.status.code(), Some(0));
    let summary_rows = written_rows();
    assert_eq!(pair_rows, summary_rows);
    assert!(
        pair_rows
            .iter()
            .any(|row| row.split('\t').nth(6).is_some_and(|c| c.contains(','))),
        "a summary lists several terms on some record"
    );
}

/// `--verbose` prints the status lines of the run, the detected input format
/// among them in Ensembl VEP's words (Parser.pm:147); without it the format is
/// not reported.
#[test]
fn verbose_prints_the_detected_format() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let verbose = run_control(&["-o", out.to_str().unwrap(), "--verbose"], dir.path());
    assert_eq!(verbose.status.code(), Some(0), "{}", text(&verbose.stderr));
    assert!(
        text(&verbose.stderr).contains("No input file format specified - detected vcf format"),
        "{}",
        text(&verbose.stderr)
    );
    let plain = run_control(&["-o", out.to_str().unwrap()], dir.path());
    assert!(!text(&plain.stderr).contains("detected vcf format"));
}

/// Under `--verbose`, a run with `--decompression_threads auto` reports the
/// worker count it resolved and the shape it read it from, and a run with
/// `--fork 0` the annotator count it resolved beside the inflaters: the
/// run-control input is a plain VCF, so the shape is unknown, the worker count
/// the `0` rule's and the annotator count the host's rule with no inflater; a
/// run with `--decompression_threads 0` and a pinned `--fork` reports neither.
#[test]
fn verbose_prints_the_auto_decompression_resolution() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let auto = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--verbose",
            "--decompression_threads",
            "auto",
            "--fork",
            "0",
        ],
        dir.path(),
    );
    assert_eq!(auto.status.code(), Some(0), "{}", text(&auto.stderr));
    let workers =
        vep_cli::config::host_parallelism().min(vep_cli::config::DECOMPRESSION_THREADS_AUTO_CAP);
    assert!(
        text(&auto.stderr).contains(&format!(
            "--decompression_threads auto resolved {workers} workers from an input of unknown shape"
        )),
        "{}",
        text(&auto.stderr)
    );
    let annotators = vep_cli::config::resolve_fork(0);
    assert!(
        text(&auto.stderr).contains(&format!(
            "--fork 0 resolved {annotators} annotators beside 0 inflaters from an input of unknown shape"
        )),
        "{}",
        text(&auto.stderr)
    );
    let host_capped = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--verbose",
            "--decompression_threads",
            "0",
            "--fork",
            "2",
        ],
        dir.path(),
    );
    assert_eq!(
        host_capped.status.code(),
        Some(0),
        "{}",
        text(&host_capped.stderr)
    );
    assert!(!text(&host_capped.stderr).contains("--decompression_threads auto resolved"));
    assert!(!text(&host_capped.stderr).contains("--fork 0 resolved"));
}

/// The status lines of `stderr` at `level` (`INFO` or `DEBUG`), each without
/// its timestamp: `<LEVEL> <target>: <message>`. The `INFO:` notices and the
/// summary line carry no timestamp and are not among them.
fn status_lines(stderr: &str, level: &str) -> Vec<String> {
    stderr
        .lines()
        .filter_map(|line| {
            let (stamp, rest) = line.split_once(' ')?;
            let timestamp = stamp.contains('T') && stamp.ends_with('Z');
            timestamp.then(|| rest.trim_start().to_string())
        })
        .filter(|line| line.starts_with(level) && line[level.len()..].starts_with(' '))
        .collect()
}

/// A `--verbose` run of the run-control corpus on one annotation thread,
/// writing to `out` with `NO_COLOR` unset: its INFO status lines are the same
/// on every host but for their timestamps.
fn verbose_run_on_one_thread(out: &Path) -> Command {
    let mut command = vep();
    command
        .args(args_for(&run_control_dir()))
        .args(["--verbose", "--fork", "1", "-o"])
        .arg(out)
        .env_remove("NO_COLOR")
        .stdin(Stdio::null());
    command
}

/// Asserts that `stderr`, the standard error of [`verbose_run_on_one_thread`]
/// writing to `out`, is plain text: no escape character anywhere, every INFO
/// status line whole after its timestamp, and the fields of a DEBUG status
/// line as `name=value`.
fn assert_plain_status_lines(stderr: &str, out: &Path) {
    assert!(!stderr.contains('\x1b'), "{stderr:?}");
    let args = args_for(&run_control_dir());
    let value_of = |flag: &str| &args[args.iter().position(|a| a == flag).unwrap() + 1];
    assert_eq!(
        status_lines(stderr, "INFO"),
        [
            format!("INFO vep_cli::runner: input: {}", value_of("-i")),
            format!("INFO vep_cli::runner: output: {}", out.display()),
            "INFO vep_cli::runner: output format: vep".to_string(),
            format!(
                "INFO vep_cli::runner: Transcript cache: {}/homo_sapiens/116_GRCh38 \
                 (lazy per-chromosome)",
                value_of("--dir_cache")
            ),
            "INFO vep_cli::transcript_index:   chromosome 21: 2276 transcripts".to_string(),
            "INFO vep_cli::transcript_index:   chromosome 22: 90 transcripts".to_string(),
        ],
        "{stderr}"
    );
    let tasks: Vec<String> = status_lines(stderr, "DEBUG")
        .iter()
        .filter_map(|line| {
            line.strip_prefix("DEBUG vep_cli::pipeline: batch handed to the pool seq=0 tasks=")
        })
        .map(str::to_string)
        .collect();
    assert_eq!(tasks.len(), 1, "{stderr}");
    assert!(tasks[0].parse::<u32>().is_ok(), "{stderr}");
}

/// The status lines a run sends to a file are plain text, `<timestamp> <LEVEL>
/// <target>: <message>`: the status layer colours its lines on a terminal
/// alone, and Ensembl VEP's status lines carry no colour code wherever they go
/// (BaseVEP.pm:785).
#[test]
fn status_lines_sent_to_a_file_are_plain_text() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let log = dir.path().join("run.log");
    let status = verbose_run_on_one_thread(&out)
        .stdout(Stdio::null())
        .stderr(std::fs::File::create(&log).unwrap())
        .status()
        .unwrap();
    let stderr = std::fs::read_to_string(&log).unwrap();
    assert_eq!(status.code(), Some(0), "{stderr}");
    assert_plain_status_lines(&stderr, &out);
}

/// The status lines a run sends to a pipe are the plain text a file takes.
#[test]
fn status_lines_sent_to_a_pipe_are_plain_text() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = verbose_run_on_one_thread(&out).output().unwrap();
    let stderr = text(&output.stderr);
    assert_eq!(output.status.code(), Some(0), "{stderr}");
    assert_plain_status_lines(&stderr, &out);
}

/// A status line names a chromosome as the cache spells it, with nothing put
/// before the name, and counts one transcript as `1 transcript`. The fly
/// corpus's cache names its sequences `4`, `2L`, `rDNA`, `mitochondrion_genome`
/// and by scaffold names; its own run on one annotation thread loads them in
/// the order of the input's first record on each, `MT` read through the
/// synonyms file. Each count is the number of distinct transcript ids in the
/// shards under that name.
#[test]
fn status_lines_name_a_chromosome_as_the_cache_spells_it() {
    let corpus =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/BDGP6.54-fly");
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = vep()
        .current_dir(&corpus)
        .args([
            "-i",
            "variants.vcf",
            "--offline",
            "--json_cache",
            "json_cache",
            "--species",
            "drosophila_melanogaster",
            "--assembly",
            "BDGP6.54",
            "--fasta",
            "reference.fa.gz",
            "--synonyms",
            "synonyms.txt",
            "--no_stats",
            "--fork",
            "1",
            "-o",
        ])
        .arg(&out)
        .env_remove("NO_COLOR")
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let stderr = text(&output.stderr);
    assert_eq!(output.status.code(), Some(0), "{stderr}");
    assert_eq!(
        status_lines(&stderr, "INFO"),
        [
            "INFO vep_cli::runner: Using reference FASTA: reference.fa.gz".to_string(),
            "INFO vep_cli::runner: input: variants.vcf".to_string(),
            format!("INFO vep_cli::runner: output: {}", out.display()),
            "INFO vep_cli::runner: output format: vep".to_string(),
            "INFO vep_cli::runner: Transcript cache: json_cache (lazy per-chromosome)".to_string(),
            "INFO vep_cli::transcript_index:   chromosome 4: 109 transcripts".to_string(),
            "INFO vep_cli::transcript_index:   chromosome 211000022278760: 1 transcript"
                .to_string(),
            "INFO vep_cli::transcript_index:   chromosome 211000022278978: 2 transcripts"
                .to_string(),
            "INFO vep_cli::transcript_index:   chromosome 211000022279188: 1 transcript"
                .to_string(),
            "INFO vep_cli::transcript_index:   chromosome 211000022280328: 1 transcript"
                .to_string(),
            "INFO vep_cli::transcript_index:   chromosome Y: 1 transcript".to_string(),
            "INFO vep_cli::transcript_index:   chromosome mitochondrion_genome: 38 transcripts"
                .to_string(),
            "INFO vep_cli::transcript_index:   chromosome Unmapped_Scaffold_8_D1580_D1567: \
             4 transcripts"
                .to_string(),
            "INFO vep_cli::transcript_index:   chromosome rDNA: 2 transcripts".to_string(),
            "INFO vep_cli::transcript_index:   chromosome Unmapped_Scaffold_60_D1601: \
             1 transcript"
                .to_string(),
            "INFO vep_cli::transcript_index:   chromosome Y_mapped_Scaffold_9_D1573: \
             1 transcript"
                .to_string(),
            "INFO vep_cli::transcript_index:   chromosome 2L: 567 transcripts".to_string(),
        ],
        "{stderr}"
    );
}

/// The counts of the `\rProcessed <count> variants...` progress lines of
/// `stderr`, in order; the summary line, `Processed <count> variants (`, is not
/// one.
fn progress_counts(stderr: &str) -> Vec<u64> {
    stderr
        .split('\r')
        .skip(1)
        .filter_map(|segment| {
            let (count, tail) = segment.strip_prefix("Processed ")?.split_once(' ')?;
            tail.starts_with("variants...")
                .then(|| count.parse().ok())
                .flatten()
        })
        .collect()
}

/// The last line of `stderr`, whole: the summary line of a run that printed
/// one, which stands on a line of its own after the progress lines.
fn summary_line(stderr: &str) -> &str {
    stderr.lines().last().unwrap_or("")
}

/// `output` without the two header lines that differ between two runs of one
/// input: the time and the command line.
fn without_run_lines(output: &str) -> Vec<&str> {
    output
        .lines()
        .filter(|l| {
            !l.starts_with("## Output produced at ") && !l.starts_with("## VEP command-line: ")
        })
        .collect()
}

/// `--no_progress` (Config.pm:161, declared and never read: Ensembl VEP's `vep`
/// prints no progress bar) switches off the progress line and nothing else. The
/// line, `\rProcessed <count> variants...`, is written on standard error each
/// time the count of annotated variants passes another 10,000, so 20,000
/// records print a count of at least 10,000. With the flag no `\rProcessed` is
/// written, while the `INFO:` notice,
/// the status lines and the summary line stay and the output is the same but
/// for the header lines naming the time and the command line.
#[test]
fn no_progress_silences_the_progress_line_and_keeps_info() {
    const NOTICE: &str = "INFO: --shift_hgvs has been set to zero, setting shift_3prime to 0";
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<(&str, u64)> = (0..20_000).map(|i| ("21", 100 + i)).collect();
    let input = intergenic_vcf(dir.path(), "many.vcf", &records);
    let args = corpus_args_with_input(&input);
    let out = dir.path().join("out.txt");

    let with_progress = run_with(&args, &["--shift_hgvs", "0"], &out);
    assert_eq!(
        with_progress.status.code(),
        Some(0),
        "{}",
        text(&with_progress.stderr)
    );
    let stderr = text(&with_progress.stderr);
    assert!(
        progress_counts(&stderr)
            .iter()
            .any(|&count| count >= 10_000),
        "a progress line with a count of at least 10,000:\n{stderr}"
    );
    let output = std::fs::read_to_string(&out).unwrap();

    let silenced = run_with(&args, &["--shift_hgvs", "0", "--no_progress"], &out);
    assert_eq!(
        silenced.status.code(),
        Some(0),
        "{}",
        text(&silenced.stderr)
    );
    let silenced_stderr = text(&silenced.stderr);
    assert!(
        !silenced_stderr.contains("\rProcessed"),
        "{silenced_stderr}"
    );
    for stderr in [&stderr, &silenced_stderr] {
        assert!(stderr.lines().any(|l| l == NOTICE), "{stderr}");
        assert!(!status_lines(stderr, "INFO").is_empty(), "{stderr}");
        assert!(
            summary_line(stderr).starts_with("Processed 20000 variants ("),
            "{stderr}"
        );
    }
    assert_eq!(
        without_run_lines(&std::fs::read_to_string(&out).unwrap()),
        without_run_lines(&output),
        "the output is unchanged"
    );
}

/// A progress line is due each time the count of annotated variants passes
/// another 10,000: 25,000 records print exactly `\rProcessed 10000 variants...`
/// and `\rProcessed 20000 variants...` on both coordinators, the batched one as
/// its count grows by whole batches of 5,000 (`progress_due`), the streaming
/// one as its tasks finish (`multiple_passed`), the second line after its last
/// batch is handed out. 12,000 records, all handed out before 10,000 of them
/// are done, print the one line in both modes. `--quiet` prints nothing at all
/// and `--no_progress` neither a progress line nor the carriage return that
/// opens one, in either mode.
#[test]
fn progress_line_appears_every_ten_thousand_variants_unless_quiet_or_no_progress() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let run = |args: &[String], extra: &[&str]| {
        let output = run_with(args, extra, &out);
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        text(&output.stderr)
    };

    let records: Vec<(&str, u64)> = (0..25_000).map(|i| ("21", 100 + i)).collect();
    let args = corpus_args_with_input(&intergenic_vcf(dir.path(), "many.vcf", &records));
    for mode in ["batched", "streaming"] {
        let stderr = run(&args, &["--pipeline_mode", mode]);
        assert_eq!(
            progress_counts(&stderr),
            [10_000, 20_000],
            "{mode}:\n{stderr}"
        );
        assert!(
            summary_line(&stderr).starts_with("Processed 25000 variants ("),
            "{mode}:\n{stderr}"
        );
        assert_eq!(
            run(&args, &["--pipeline_mode", mode, "--quiet"]),
            "",
            "{mode}"
        );
        let silenced = run(&args, &["--pipeline_mode", mode, "--no_progress"]);
        assert!(!silenced.contains('\r'), "{mode}:\n{silenced}");
        assert!(
            summary_line(&silenced).starts_with("Processed 25000 variants ("),
            "{silenced}"
        );
    }

    let records: Vec<(&str, u64)> = (0..12_000).map(|i| ("21", 100 + i)).collect();
    let args = corpus_args_with_input(&intergenic_vcf(dir.path(), "fewer.vcf", &records));
    for mode in ["batched", "streaming"] {
        let stderr = run(&args, &["--pipeline_mode", mode]);
        assert_eq!(progress_counts(&stderr), [10_000], "{mode}:\n{stderr}");
        assert!(
            summary_line(&stderr).starts_with("Processed 12000 variants ("),
            "{mode}:\n{stderr}"
        );
    }
}

/// A line a run writes after a progress line opens a line of its own. The
/// progress line is written after a carriage return and has no line end, so a
/// line end closes it before any other line: the `ERROR:` line of a run that
/// fails, the echo of a warning, a status line and the summary line. Three
/// inputs show it in both pipeline modes, the first two of 25,000 records. A
/// compressed input cut inside its deflate stream ends the run with exit 1. A
/// record 15,000 on a chromosome the cache lacks is echoed, and a warnings
/// file that fails at its first write then ends the run with exit 1: Linux's
/// `/dev/full`, a device, which the run does not test before its first
/// record, so that run is made on Linux alone. Three records
/// on a second chromosome follow 35,000 on the first, and at a pipeline depth
/// of 1 at most four batches stand between the reader and the batches done,
/// so that chromosome's status line follows a progress line. The batched
/// coordinator reaches each failure between whole batches of 5,000, so its
/// progress lines are the ones before the failure.
#[test]
fn a_line_written_after_a_progress_line_opens_a_line_of_its_own() {
    const SKIPPED: &str = "WARNING: line 15000 skipped (Z 15099 r15000 A G . . .): \
        Chromosome Z not found in annotation sources or synonyms; \
        chromosome Z does not overlap any features";
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let run = |args: &[String], extra: &[&str]| {
        let output = run_with(args, extra, &out);
        (output.status.code(), text(&output.stderr))
    };
    let last_line = |stderr: &str| stderr.lines().last().unwrap_or("").to_string();

    let many: Vec<(&str, u64)> = (0..25_000).map(|i| ("21", 100 + i)).collect();
    let whole = std::fs::read(intergenic_vcf(dir.path(), "many.vcf", &many)).unwrap();
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(&whole).unwrap();
    let mut compressed = encoder.finish().unwrap();
    compressed.truncate(compressed.len() * 9 / 10);
    let cut = dir.path().join("cut.vcf.gz");
    std::fs::write(&cut, compressed).unwrap();
    let cut_args = corpus_args_with_input(&cut);

    let warned: Vec<(&str, u64)> = (0..25_000)
        .map(|i| (if i == 14_999 { "Z" } else { "21" }, 100 + i))
        .collect();
    let warned_args = corpus_args_with_input(&intergenic_vcf(dir.path(), "warned.vcf", &warned));
    let full_device = Path::new("/dev/full");
    let unwritable = format!(
        "ERROR: [output] Could not write to warnings file {}",
        full_device.display()
    );

    let second: Vec<(&str, u64)> = (0..35_003)
        .map(|i| (if i < 35_000 { "21" } else { "22" }, 100 + i))
        .collect();
    let second_args = corpus_args_with_input(&intergenic_vcf(dir.path(), "second.vcf", &second));

    for mode in ["batched", "streaming"] {
        let (code, stderr) = run(&cut_args, &["--pipeline_mode", mode]);
        assert_eq!(code, Some(1), "{mode}:\n{stderr}");
        assert!(
            last_line(&stderr).starts_with("ERROR: [run] Failed to read input line: "),
            "{mode}:\n{stderr}"
        );
        if mode == "batched" {
            assert_eq!(progress_counts(&stderr), [10_000, 20_000], "{stderr}");
            assert!(
                stderr.contains("\rProcessed 20000 variants...\nERROR: [run] "),
                "{stderr}"
            );
        }

        if cfg!(target_os = "linux") {
            let (code, stderr) = run(
                &warned_args,
                &[
                    "--pipeline_mode",
                    mode,
                    "--warning_file",
                    full_device.to_str().unwrap(),
                ],
            );
            assert_eq!(code, Some(1), "{mode}:\n{stderr}");
            assert!(stderr.lines().any(|l| l == SKIPPED), "{mode}:\n{stderr}");
            assert_eq!(last_line(&stderr), unwritable, "{mode}:\n{stderr}");
            if mode == "batched" {
                assert_eq!(progress_counts(&stderr), [10_000], "{stderr}");
                assert!(
                    stderr.contains("\rProcessed 10000 variants...\nWARNING: line 15000 skipped "),
                    "{stderr}"
                );
            }
        }

        let (code, stderr) = run(&warned_args, &["--pipeline_mode", mode]);
        assert_eq!(code, Some(0), "{mode}:\n{stderr}");
        assert!(stderr.lines().any(|l| l == SKIPPED), "{mode}:\n{stderr}");
        assert!(!progress_counts(&stderr).is_empty(), "{mode}:\n{stderr}");
        assert!(
            last_line(&stderr).starts_with("Processed 24999 variants ("),
            "{mode}:\n{stderr}"
        );

        let (code, stderr) = run(
            &second_args,
            &["--pipeline_mode", mode, "--pipeline_depth", "1"],
        );
        assert_eq!(code, Some(0), "{mode}:\n{stderr}");
        let loaded: Vec<&str> = stderr
            .lines()
            .filter(|l| l.contains("INFO vep_cli::transcript_index:   chromosome 22: "))
            .collect();
        assert_eq!(loaded.len(), 1, "{mode}:\n{stderr}");
        assert!(!loaded[0].contains('\r'), "{mode}:\n{stderr}");
        let (before, _) = stderr.split_once(loaded[0]).unwrap();
        assert!(before.ends_with(" variants...\n"), "{mode}:\n{stderr}");
        assert!(
            last_line(&stderr).starts_with("Processed 35003 variants ("),
            "{mode}:\n{stderr}"
        );
    }
}

/// A progress line that cannot be written is dropped and the run goes on. The
/// reader of standard error goes once the first progress line has arrived, so
/// every later line fails to write: either coordinator still writes all
/// 60,000 records to the output, and the run ends in success past the
/// line end and the summary line it cannot write. A streaming task
/// that panicked on the failed write would never send its chunk, which the
/// writer waits for.
#[test]
fn a_progress_line_that_cannot_be_written_leaves_the_output_complete() {
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<(&str, u64)> = (0..60_000).map(|i| ("21", 100 + i)).collect();
    let args = corpus_args_with_input(&intergenic_vcf(dir.path(), "many.vcf", &records));
    for mode in ["streaming", "batched"] {
        let out = dir.path().join(format!("{mode}.txt"));
        let status = status_once_the_reader_goes(&args, mode, &out);
        assert!(status.success(), "{mode}: {status}");
        let output = std::fs::read_to_string(&out).unwrap();
        let written: BTreeSet<&str> = data_lines(&output)
            .iter()
            .map(|row| row.split('\t').next().unwrap())
            .collect();
        assert_eq!(written.len(), 60_000, "{mode}");
    }
}

/// The exit status of a run of `args` in pipeline mode `mode` writing to
/// `out`, whose reader of standard error goes once the first progress line
/// has arrived.
fn status_once_the_reader_goes(args: &[String], mode: &str, out: &Path) -> ExitStatus {
    const FIRST: &[u8] = b"\rProcessed 10000 variants...";
    let mut child = vep()
        .args(args)
        .args(["--pipeline_mode", mode, "-o", out.to_str().unwrap()])
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut stderr = child.stderr.take().unwrap();
    let mut seen = Vec::new();
    let mut byte = [0u8; 1];
    while !seen.ends_with(FIRST) {
        assert_eq!(
            stderr.read(&mut byte).unwrap(),
            1,
            "standard error ended before the first progress line:\n{}",
            text(&seen)
        );
        seen.push(byte[0]);
    }
    drop(stderr);
    status_before_the_deadline(child)
}

/// The exit status of `child`, a run whose standard error cannot be written:
/// one still running after two minutes is killed and fails the test.
fn status_before_the_deadline(mut child: Child) -> ExitStatus {
    let deadline = Instant::now() + Duration::from_secs(120);
    loop {
        if let Some(status) = child.try_wait().unwrap() {
            break status;
        }
        if Instant::now() > deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("the run did not end once standard error could not be written");
        }
        std::thread::sleep(Duration::from_millis(20));
    }
}

/// The streams that take no bytes on this host: a device that refuses every
/// write, where the host has one, and a pipe whose reader has gone.
fn unwritable_streams() -> &'static [&'static str] {
    if Path::new("/dev/full").exists() {
        &["/dev/full", "a pipe whose reader has gone"]
    } else {
        &["a pipe whose reader has gone"]
    }
}

/// A standard error that takes no bytes, as [`unwritable_streams`] names it.
fn unwritable(stream: &str) -> Stdio {
    if stream == "/dev/full" {
        return std::fs::File::options()
            .write(true)
            .open(stream)
            .unwrap()
            .into();
    }
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    writer.into()
}

/// Standard error on a stream that takes no bytes never ends a run, a device
/// that refuses every write (`/dev/full`) or a pipe whose reader has gone.
/// Ensembl VEP's one write there on a run that completes, the `warn` of
/// `warning_msg` (BaseVEP.pm:821), fails on such a device without ending the
/// run. In both pipeline modes a run of 12,000 records with a skipped line
/// (the banner, the status lines, a warning's echo, a progress line and the
/// summary line all refused) exits 0 with the output and the warnings file
/// of the same run on a readable stream, `--warning_file STDERR` included. A
/// refused run keeps its exit status of 1 though its `ERROR:` line is lost.
#[test]
fn an_unwritable_error_stream_never_fails_a_run() {
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<(&str, u64)> = (0..12_000)
        .map(|i| {
            if i == 7 {
                ("Z", 1_000)
            } else {
                ("21", 100 + i)
            }
        })
        .collect();
    let args = corpus_args_with_input(&intergenic_vcf(dir.path(), "many.vcf", &records));
    let run = |args: &[String], extra: &[&str], out: &Path, stderr: Stdio| {
        vep()
            .args(args)
            .args(extra)
            .args(["-o", out.to_str().unwrap()])
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(stderr)
            .status()
            .unwrap()
    };
    let left = |out: &Path| {
        let warnings = PathBuf::from(format!("{}_warnings.txt", out.display()));
        (
            std::fs::read_to_string(out)
                .map(|text| without_run_lines(&text).join("\n"))
                .ok(),
            std::fs::read_to_string(warnings).ok(),
        )
    };
    let streams = unwritable_streams();
    let mut runs = 0;
    let mut fresh = || {
        runs += 1;
        dir.path().join(format!("out{runs}.txt"))
    };
    for mode in ["streaming", "batched"] {
        for target in [&[][..], &["--warning_file", "STDERR"][..]] {
            let extra = [&["--pipeline_mode", mode][..], target].concat();
            let readable = fresh();
            let status = run(&args, &extra, &readable, Stdio::null());
            assert!(status.success(), "{mode} {target:?}: {status}");
            let expected = left(&readable);
            assert!(expected.0.is_some(), "{mode} {target:?}");
            assert_eq!(expected.1.is_some(), target.is_empty(), "{mode}");
            for stream in streams {
                let out = fresh();
                let status = run(&args, &extra, &out, unwritable(stream));
                assert!(status.success(), "{mode} {target:?}, {stream}: {status}");
                assert_eq!(left(&out), expected, "{mode} {target:?}, {stream}");
            }
        }
    }

    let refused = corpus_args_with_input(&dir.path().join("no_such.vcf"));
    for stream in streams {
        let out = dir.path().join("refused.txt");
        let status = run(&refused, &[], &out, unwritable(stream));
        assert_eq!(status.code(), Some(1), "{stream}: {status}");
    }
}

/// A line a plugin prints on standard error is dropped by such a stream too,
/// as a bare `warn` in an Ensembl VEP plugin fails there without ending the
/// run. dbNSFP's bare-version warning at its load (dbNSFP.pm:208) and its
/// per-row list errors under `transcript_match=1` (dbNSFP.pm:379, :384) are
/// refused, and each run exits 0 with the output and the warnings file of the
/// same run on a readable stream.
#[test]
fn an_unwritable_error_stream_never_fails_a_run_whose_plugin_warns() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-plugin-dbnsfp-fixture");
    let bare = dir.path().join("dbNSFP4.9_grch38_synthetic.gz");
    for suffix in ["", ".tbi"] {
        std::fs::copy(
            corpus.join(format!(
                "plugin_data/dbNSFP4.9a_grch38_synthetic.gz{suffix}"
            )),
            dir.path()
                .join(format!("dbNSFP4.9_grch38_synthetic.gz{suffix}")),
        )
        .unwrap();
    }
    let lists = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/dbnsfp-lists");
    let cases = [
        (
            "the load",
            args_for(&corpus),
            format!("dbNSFP,{},SIFT_score", bare.display()),
            "WARNING: Specific suffix 'a' or 'c' not found for version 4.9. Proceeding with \
             general version.",
        ),
        (
            "a row",
            with_input(args_for(&corpus), &lists.join("variants.vcf")),
            format!(
                "dbNSFP,{},transcript_match=1,SIFT_score,SIFT_pred",
                lists.join("dbNSFP4.9a_grch38_lists_synthetic.gz").display()
            ),
            "ERROR: Transcript index out of range for field SIFT_score",
        ),
    ];
    let run = |args: &[String], plugin: &str, out: &Path, stderr: Stdio| {
        vep()
            .args(args)
            .args(["--quiet", "--plugin", plugin])
            .args(["-o", out.to_str().unwrap()])
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(stderr)
            .spawn()
            .unwrap()
    };
    let left = |out: &Path| {
        let warnings = PathBuf::from(format!("{}_warnings.txt", out.display()));
        (
            std::fs::read_to_string(out)
                .map(|text| without_run_lines(&text).join("\n"))
                .ok(),
            std::fs::read_to_string(warnings).ok(),
        )
    };
    let mut runs = 0;
    let mut fresh = || {
        runs += 1;
        dir.path().join(format!("out{runs}.txt"))
    };
    for (case, args, plugin, line) in &cases {
        let readable = fresh();
        let printed = run(args, plugin, &readable, Stdio::piped())
            .wait_with_output()
            .unwrap();
        assert!(printed.status.success(), "{case}: {}", printed.status);
        let stderr = String::from_utf8_lossy(&printed.stderr);
        assert!(
            stderr.lines().any(|l| l == *line),
            "{case}: the plugin's line on a readable stream:\n{stderr}"
        );
        let expected = left(&readable);
        assert!(expected.0.is_some(), "{case}");
        for stream in unwritable_streams() {
            let out = fresh();
            let status = status_before_the_deadline(run(args, plugin, &out, unwritable(stream)));
            assert!(status.success(), "{case}, {stream}: {status}");
            assert_eq!(left(&out), expected, "{case}, {stream}");
        }
    }
}

/// The variable that carries a library run's arguments, one per line, to the
/// process [`library_run`] starts.
const LIBRARY_RUN: &str = "LIBRARY_RUN_ARGUMENTS";

/// The body of the process [`library_run`] starts: the library runs the
/// arguments [`LIBRARY_RUN`] carries, in a process that never calls
/// `vep_cli::stderr::drop_refused_writes`. Without the variable, as in a run
/// of these tests, it does nothing.
#[test]
fn library_run_from_the_environment() {
    let Ok(args) = std::env::var(LIBRARY_RUN) else {
        return;
    };
    let argv = std::iter::once("vep").chain(args.lines());
    let config = vep_cli::Config::from_argv(argv).unwrap();
    vep_cli::runner::Runner::new(config).unwrap().run().unwrap();
}

/// This test binary started again to run `args` through the library, in the
/// environment [`vep`] gives the binary: its harness runs
/// [`library_run_from_the_environment`] alone, under the flags `harness`.
fn library_run(args: &[String], harness: &[&str]) -> Command {
    let mut cmd = Command::new(std::env::current_exe().unwrap());
    for (key, value) in vep().get_envs() {
        match value {
            Some(value) => cmd.env(key, value),
            None => cmd.env_remove(key),
        };
    }
    cmd.env_remove("RUST_TEST_NOCAPTURE")
        .env(LIBRARY_RUN, args.join("\n"))
        .args(["--exact", "library_run_from_the_environment"])
        .args(harness)
        .stdin(Stdio::null());
    cmd
}

/// A process that runs the library and never calls
/// `vep_cli::stderr::drop_refused_writes` writes the run's lines as `eprint!`
/// does, so a test harness captures them. A batched run of 12,000 records
/// with a skipped line puts a warning, a progress line and the summary
/// line on standard error under a harness that does not capture
/// (`--nocapture`), and none of them under one that does.
#[test]
fn a_library_runs_lines_are_captured_by_a_test_harness() {
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<(&str, u64)> = (0..12_000)
        .map(|i| {
            if i == 7 {
                ("Z", 1_000)
            } else {
                ("21", 100 + i)
            }
        })
        .collect();
    let out = dir.path().join("out.txt");
    let mut args = corpus_args_with_input(&intergenic_vcf(dir.path(), "many.vcf", &records));
    args.extend(["--pipeline_mode", "batched", "-o", out.to_str().unwrap()].map(str::to_string));
    let stderr_under = |harness: &[&str]| {
        let output = library_run(&args, harness).output().unwrap();
        assert!(
            output.status.success(),
            "{harness:?}:\n{}",
            text(&output.stdout)
        );
        text(&output.stderr)
    };

    let shown = stderr_under(&["--nocapture"]);
    assert!(shown.lines().any(|l| l.starts_with("WARNING: ")), "{shown}");
    assert_eq!(progress_counts(&shown).len(), 1, "{shown}");
    assert!(
        summary_line(&shown).starts_with("Processed 11999 variants ("),
        "{shown}"
    );

    let captured = stderr_under(&[]);
    assert!(
        !captured.contains("WARNING: ") && !captured.contains("Processed "),
        "{captured}"
    );
}

/// A process that runs the library and never configures the warning sink
/// reads each warning once on standard error: the sink writes its one copy
/// there and echoes none beside it. A library run of three records, the
/// second on a chromosome the cache lacks, puts that record's skip line on
/// standard error exactly once under a harness that does not capture. The two
/// copies are the binary's under `--warning_file STDERR`
/// ([`warning_file_stderr_literal`]).
#[test]
fn a_library_run_writes_each_warning_once_on_standard_error() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let records = [("21", 100), ("Z", 1_000), ("21", 300)];
    let mut args = corpus_args_with_input(&intergenic_vcf(dir.path(), "few.vcf", &records));
    args.extend(["-o", out.to_str().unwrap()].map(str::to_string));
    let output = library_run(&args, &["--nocapture"]).output().unwrap();
    assert!(output.status.success(), "{}", text(&output.stdout));
    let stderr = text(&output.stderr);
    let warnings = warning_lines(&stderr);
    assert_eq!(warnings.len(), 1, "{stderr}");
    assert!(
        warnings[0].starts_with("WARNING: line 2 skipped (Z 1000 r2 "),
        "{stderr}"
    );
}

/// The streaming coordinator's progress line is dropped when standard error
/// refuses it in every process, the `vep` binary or another: a task that
/// panicked on that write would never send its chunk, and the output would
/// end before its records. A library run of 12,000 records in a process that
/// never calls `vep_cli::stderr::drop_refused_writes`, its harness not
/// capturing, prints its one progress line on a readable stream; with
/// standard error a pipe whose reader has gone it writes every record to the
/// output all the same. Its next line there, the summary line behind the line
/// end that closes the progress line, is an `eprint!` the stream refuses,
/// whose panic fails the child's one test.
#[test]
fn a_refused_progress_line_leaves_a_library_runs_output_complete() {
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<(&str, u64)> = (0..12_000).map(|i| ("21", 100 + i)).collect();
    let input = intergenic_vcf(dir.path(), "many.vcf", &records);
    let run_to = |out: &Path| {
        let mut args = corpus_args_with_input(&input);
        args.extend(
            ["--pipeline_mode", "streaming", "-o", out.to_str().unwrap()].map(str::to_string),
        );
        library_run(&args, &["--nocapture"])
    };

    let readable = run_to(&dir.path().join("readable.txt")).output().unwrap();
    assert!(readable.status.success(), "{}", text(&readable.stdout));
    assert_eq!(progress_counts(&text(&readable.stderr)), [10_000]);

    let out = dir.path().join("out.txt");
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    let child = run_to(&out)
        .stdout(Stdio::null())
        .stderr(writer)
        .spawn()
        .unwrap();
    let status = status_before_the_deadline(child);
    assert!(!status.success(), "{status}");
    let output = std::fs::read_to_string(&out).unwrap();
    let written: BTreeSet<&str> = data_lines(&output)
        .iter()
        .map(|row| row.split('\t').next().unwrap())
        .collect();
    assert_eq!(written.len(), 12_000);
}

/// The frequency-filter corpus: records with cache frequencies, the fixture of
/// `--filter_common` and `--check_frequency`.
fn freqfilter_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-freqfilter")
}

/// The summary line, this engine's own, closes every run on standard error:
/// `Processed <N> variants (<M> consequences, <F> filtered out) in <X.X>s using
/// <T> threads`, where N counts the records annotated, M the transcript
/// consequences computed (the `Transcript` rows of a run without a pick) and F
/// the records a frequency filter dropped; the filter clause goes when F is 0,
/// the thread clause under one thread, and the nouns are singular at 1. On the
/// frequency-filter corpus `--filter_common --freq_pop AF` keeps 11 of the 15
/// records with 251 transcript rows, so N, M and F are read back from the input
/// and the output. `--quiet` withholds the line (BaseVEP.pm:782).
#[test]
fn final_summary_line_states_variants_consequences_and_filtered_count() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let args = args_for(&freqfilter_dir());
    let filtered_run = run_with(
        &args,
        &["--filter_common", "--freq_pop", "AF", "--fork", "4"],
        &out,
    );
    assert_eq!(
        filtered_run.status.code(),
        Some(0),
        "{}",
        text(&filtered_run.stderr)
    );
    let records = std::fs::read_to_string(freqfilter_dir().join("variants.vcf")).unwrap();
    let records = data_lines(&records).len();
    let output = std::fs::read_to_string(&out).unwrap();
    let rows = data_lines(&output);
    let kept: BTreeSet<&str> = rows.iter().map(|r| r.split('\t').next().unwrap()).collect();
    let transcript_rows = rows
        .iter()
        .filter(|r| r.split('\t').nth(5) == Some("Transcript"))
        .count();
    let filtered = records - kept.len();
    assert_eq!((kept.len(), transcript_rows, filtered), (11, 251, 4));
    let stderr = text(&filtered_run.stderr);
    let (head, tail) = summary_line(&stderr).split_once(" in ").unwrap();
    assert_eq!(
        head,
        format!(
            "Processed {} variants ({transcript_rows} consequences, {filtered} filtered out)",
            kept.len()
        ),
        "{stderr}"
    );
    let (seconds, threads) = tail.split_once('s').unwrap();
    assert!(
        seconds.parse::<f64>().is_ok()
            && seconds.split_once('.').is_some_and(|(_, f)| f.len() == 1),
        "{tail}"
    );
    assert_eq!(threads, " using 4 threads");

    let one = dir.path().join("one.vcf");
    std::fs::write(
        &one,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30289183\tstop_ins1_lost_rev\tC\tCG\t.\t.\t.\n",
    )
    .unwrap();
    let single = run_with(
        &with_input(args_for(&run_control_dir()), &one),
        &["--fork", "1"],
        &out,
    );
    assert_eq!(single.status.code(), Some(0), "{}", text(&single.stderr));
    let stderr = text(&single.stderr);
    assert_eq!(data_lines(&std::fs::read_to_string(&out).unwrap()).len(), 1);
    let (head, tail) = summary_line(&stderr).split_once(" in ").unwrap();
    assert_eq!(head, "Processed 1 variant (1 consequence)", "{stderr}");
    assert!(tail.ends_with('s') && !tail.contains("thread"), "{tail}");

    let quiet = run_with(
        &with_input(args_for(&run_control_dir()), &one),
        &["--quiet"],
        &out,
    );
    assert_eq!(quiet.status.code(), Some(0));
    assert!(quiet.stderr.is_empty(), "{}", text(&quiet.stderr));
}

/// The banner opens standard error on every run, three lines before the first
/// status line: a rule of three dashed columns, `Ensembl VEP (Rust)   v116.3`
/// (the Ensembl VEP release the cache targets, not the crate version) and the
/// rule again. Ensembl VEP prints its own banner only in its usage text
/// (vep:55-58). `--quiet` withholds the banner with every other status line
/// (BaseVEP.pm:782), and nothing of it reaches standard output.
#[test]
fn banner_prints_three_lines_to_stderr_unless_quiet() {
    const RULE: &str = "-------------------- -------------------- --------------------";
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = run_control(&["-o", out.to_str().unwrap()], dir.path());
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let stderr = text(&output.stderr);
    let banner: Vec<&str> = stderr.lines().take(3).collect();
    assert_eq!(
        banner,
        [RULE, "Ensembl VEP (Rust)   v116.3", RULE],
        "{stderr}"
    );
    assert!(output.stdout.is_empty());

    let quiet = run_control(&["-o", out.to_str().unwrap(), "--quiet"], dir.path());
    assert_eq!(quiet.status.code(), Some(0));
    assert!(quiet.stderr.is_empty(), "{}", text(&quiet.stderr));
}

/// Every status line goes to standard error, so `-o STDOUT` carries the output
/// alone: the banner, the `INFO:` notices of the FASTA checks in Ensembl VEP's
/// order (`Disabling --hgvs` at Runner.pm:752 before `--shift_hgvs` at :793),
/// the timestamped status lines and the summary line are on standard error,
/// while standard output is the file a run with `-o <file>` writes, but for the
/// two header lines naming the time and the command line.
#[test]
fn status_lines_go_to_stderr_with_perls_info_notices() {
    const HGVS: &str = "INFO: Disabling --hgvs; using --offline and no FASTA file found";
    const SHIFT: &str = "INFO: --shift_hgvs has been set to zero, setting shift_3prime to 0";
    let dir = tempfile::tempdir().unwrap();
    let flags = ["--everything", "--shift_hgvs", "0"];
    let to_stdout = run_control(&[&flags[..], &["-o", "STDOUT"]].concat(), dir.path());
    assert_eq!(
        to_stdout.status.code(),
        Some(0),
        "{}",
        text(&to_stdout.stderr)
    );
    let stderr = text(&to_stdout.stderr);
    let notices: Vec<&str> = stderr.lines().filter(|l| l.starts_with("INFO: ")).collect();
    assert_eq!(notices, [HGVS, SHIFT], "{stderr}");
    assert!(stderr.starts_with("--------------------"), "{stderr}");
    assert!(!status_lines(&stderr, "INFO").is_empty(), "{stderr}");
    assert!(summary_line(&stderr).starts_with("Processed "), "{stderr}");
    let stdout = text(&to_stdout.stdout);
    assert!(
        stdout
            .lines()
            .all(|l| l.starts_with('#') || l.split('\t').count() == 14),
        "a line that is neither a header nor a row:\n{stdout}"
    );

    let out = dir.path().join("out.txt");
    let to_file = run_control(
        &[&flags[..], &["-o", out.to_str().unwrap()]].concat(),
        dir.path(),
    );
    assert_eq!(to_file.status.code(), Some(0), "{}", text(&to_file.stderr));
    assert!(to_file.stdout.is_empty());
    assert_eq!(
        without_run_lines(&stdout),
        without_run_lines(&std::fs::read_to_string(&out).unwrap())
    );
}

/// `-v` is `--verbose` (Config.pm:159, `verbose|v`): the DEBUG status lines the
/// long spelling prints, the detected input format among them, print under the
/// short one too, and a plain run prints none.
#[test]
fn short_v_is_verbose() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let debug_lines = |extra: &[&str]| {
        let output = run_control(
            &[&["-o", out.to_str().unwrap()], extra].concat(),
            dir.path(),
        );
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        status_lines(&text(&output.stderr), "DEBUG")
    };
    let long = debug_lines(&["--verbose"]);
    assert!(
        long.iter()
            .any(|l| l.ends_with("No input file format specified - detected vcf format")),
        "{long:?}"
    );
    assert_eq!(debug_lines(&["-v"]), long);
    assert!(debug_lines(&[]).is_empty());
}

/// `--buffer_size` (Config.pm:118, default 5000 at :300) sets how many records
/// are annotated per batch and changes nothing in the output: 1, 7 and the
/// default write byte-identical rows for the run-control corpus. Ensembl VEP's
/// buffer never spans two chromosomes (InputBuffer.pm:181-267); a vep-rs batch
/// may, which is likewise invisible in the output.
#[test]
fn buffer_size_changes_no_output() {
    let dir = tempfile::tempdir().unwrap();
    let mut outputs = Vec::new();
    for size in ["1", "7", "5000"] {
        let out = dir.path().join(format!("out-{size}.txt"));
        let output = run_control(
            &["-o", out.to_str().unwrap(), "--buffer_size", size],
            dir.path(),
        );
        assert_eq!(
            output.status.code(),
            Some(0),
            "{size}: {}",
            text(&output.stderr)
        );
        let rows = data_lines(&std::fs::read_to_string(&out).unwrap())
            .into_iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        assert!(rows.len() > 30, "{size}: the run annotates");
        outputs.push(rows);
    }
    assert_eq!(
        outputs[0], outputs[2],
        "--buffer_size 1 against the default"
    );
    assert_eq!(
        outputs[1], outputs[2],
        "--buffer_size 7 against the default"
    );
}

/// `--decompression_threads` sets the BGZF block reader's worker count and
/// changes nothing in the output: a BGZF copy of the run-control input annotates
/// to the same rows with one worker, four, and the `0` rule's count, and to the
/// uncompressed file's rows.
#[test]
fn decompression_threads_changes_no_output() {
    let dir = tempfile::tempdir().unwrap();
    let plain = dir.path().join("plain.txt");
    let output = run_control(&["-o", plain.to_str().unwrap()], dir.path());
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let reference = std::fs::read_to_string(&plain).unwrap();
    let reference = data_lines(&reference);

    let mut writer = noodles::bgzf::io::Writer::new(Vec::new());
    writer
        .write_all(&std::fs::read(run_control_dir().join("variants.vcf")).unwrap())
        .unwrap();
    let bgzf = dir.path().join("variants.vcf.gz");
    std::fs::write(&bgzf, writer.finish().unwrap()).unwrap();
    for threads in ["1", "4", "0"] {
        let out = dir.path().join(format!("out-{threads}.txt"));
        let output = vep()
            .args(with_input(args_for(&run_control_dir()), &bgzf))
            .args([
                "-o",
                out.to_str().unwrap(),
                "--decompression_threads",
                threads,
            ])
            .current_dir(dir.path())
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{threads}: {}",
            text(&output.stderr)
        );
        let rows = std::fs::read_to_string(&out).unwrap();
        assert_eq!(
            data_lines(&rows),
            reference,
            "--decompression_threads {threads}"
        );
    }
}

/// The line `--var_synonyms` without `--check_existing` prints (Runner.pm:117-120).
const VAR_SYNONYMS_NOTICE: &str =
    "INFO: Colocated variants will be provided when using --var_synonyms";

/// `--var_synonyms` without `--check_existing` prints `INFO: Colocated variants
/// will be provided when using --var_synonyms` on standard error as a line of its
/// own and switches co-location on (Runner.pm:117-120), so the output carries the
/// `VAR_SYNONYMS` key and Existing_variation; `--quiet` withholds the notice as
/// `status_msg` does (BaseVEP.pm:782), and `--check_existing` given beside the
/// flag leaves nothing to switch on and no notice. A refused command line prints
/// its error alone: `Config.pm:769` refuses it before the runner's implication.
#[test]
fn var_synonyms_prints_the_notice_and_switches_check_existing_on() {
    let colocated =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-colocated");
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let run = |extra: &[&str]| {
        vep()
            .args(args_for(&colocated))
            .args(["-o", out.to_str().unwrap()])
            .args(extra)
            .current_dir(dir.path())
            .stdin(Stdio::null())
            .output()
            .unwrap()
    };
    let implied = run(&["--var_synonyms"]);
    assert_eq!(implied.status.code(), Some(0), "{}", text(&implied.stderr));
    let stderr = text(&implied.stderr);
    assert!(
        stderr.lines().any(|l| l == VAR_SYNONYMS_NOTICE),
        "the notice is a line of its own:\n{stderr}"
    );
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(output.contains("## VAR_SYNONYMS : List of known variation synonyms and their sources"));
    assert!(
        output.contains("VAR_SYNONYMS=ClinVar::RCV000886041,VCV000713837"),
        "{output}"
    );
    assert!(
        data_lines(&output)
            .iter()
            .any(|l| l.split('\t').nth(12) == Some("rs374325784")),
        "co-location is on"
    );

    // The configuration's notice precedes the lines of `post_setup_checks`, as
    // `init` prints them (Runner.pm:117-120, then :138).
    let ordered = run(&["--var_synonyms", "--shift_hgvs", "0"]);
    assert_eq!(ordered.status.code(), Some(0), "{}", text(&ordered.stderr));
    let stderr = text(&ordered.stderr);
    let lines: Vec<&str> = stderr.lines().collect();
    let notice_at = lines
        .iter()
        .position(|l| *l == VAR_SYNONYMS_NOTICE)
        .unwrap_or_else(|| panic!("no var_synonyms notice:\n{stderr}"));
    let shift_at = lines
        .iter()
        .position(|l| *l == "INFO: --shift_hgvs has been set to zero, setting shift_3prime to 0")
        .unwrap_or_else(|| panic!("no shift_hgvs line:\n{stderr}"));
    assert!(notice_at < shift_at, "{stderr}");

    let quiet = run(&["--var_synonyms", "--quiet"]);
    assert_eq!(quiet.status.code(), Some(0), "{}", text(&quiet.stderr));
    assert!(!text(&quiet.stderr).contains("INFO: Colocated"));
    assert!(std::fs::read_to_string(&out)
        .unwrap()
        .contains("VAR_SYNONYMS=ClinVar::RCV000886041,VCV000713837"));

    let explicit = run(&["--var_synonyms", "--check_existing"]);
    assert_eq!(
        explicit.status.code(),
        Some(0),
        "{}",
        text(&explicit.stderr)
    );
    assert!(!text(&explicit.stderr).contains("INFO: Colocated"));

    let refused = vep()
        .args(["--var_synonyms", "--json", "--vcf"])
        .current_dir(dir.path())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr),
        "ERROR: [cli] Can't use --json and --vcf together\n"
    );
}

/// `var_synonyms 0` in a `--config` file or `vep.ini`, and `VEP_VAR_SYNONYMS=0`,
/// set the flag defined but false: `Config.pm:1101` and `:1158` store the `0`,
/// Runner.pm:117 tests the flag with `defined`, so the notice prints and
/// `check_existing` is switched on, while OutputFactory.pm:1090 reads the `0` as
/// false and writes no `VAR_SYNONYMS`. On `21 14967531 T C` Ensembl VEP 116.3
/// prints the notice, the row with `rs139263261`, `CLIN_SIG=benign;PHENO=1`, no
/// `## VAR_SYNONYMS` header line and `--var_synonyms 0` in the command-line
/// header, where `full_command` prints a stored value other than `1` after its
/// key (Config.pm:745); a quoted empty value, `var_synonyms ""`, is stored as the
/// empty string (`:1092` strips the quotes) and prints `--var_synonyms` with an
/// empty value, the trailing blank trimmed at the end of the line (`:749`).
/// `var_synonyms 1` from the file is the flag itself, printed bare, and the bare
/// record, with no flag from any source, prints neither the notice nor a
/// co-located variant.
#[test]
fn var_synonyms_zero_from_a_source_is_defined_but_false() {
    #[derive(Debug, PartialEq)]
    struct Observed {
        notice: bool,
        header: bool,
        existing: BTreeSet<String>,
        clinical: bool,
        flag_printed: String,
    }
    let colocated =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-colocated");
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("one.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t14967531\t.\tT\tC\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.path().join("out.txt");
    let run = |extra: &[&str], env: &[(&str, &str)]| {
        let mut command = vep();
        command
            .args(args_for(&colocated))
            .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
            .args(extra)
            .current_dir(dir.path())
            .stdin(Stdio::null());
        for (key, value) in env {
            command.env(key, value);
        }
        let output = command.output().unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{extra:?} {env:?}: {}",
            text(&output.stderr)
        );
        let notice = text(&output.stderr)
            .lines()
            .any(|l| l == VAR_SYNONYMS_NOTICE);
        let written = std::fs::read_to_string(&out).unwrap();
        let header = written.lines().any(|l| l.starts_with("## VAR_SYNONYMS"));
        let command_line = written
            .lines()
            .find(|l| l.starts_with("## VEP command-line:"))
            .unwrap();
        // `var_synonyms` sorts last among these flags, so the line ends with
        // the flag and, when it has one, its value.
        let flag_printed = command_line
            .find("--var_synonyms")
            .map_or(String::new(), |at| command_line[at..].to_string());
        assert!(
            !command_line.ends_with(' '),
            "no trailing blank (Config.pm:749): {command_line:?}"
        );
        let rows = data_lines(&written);
        let existing: BTreeSet<String> = rows
            .iter()
            .map(|l| l.split('\t').nth(12).unwrap().to_string())
            .collect();
        let clinical = rows.iter().all(|l| l.contains("CLIN_SIG=benign;PHENO=1"));
        Observed {
            notice,
            header,
            existing,
            clinical,
            flag_printed,
        }
    };
    let known = |name: &str| BTreeSet::from([name.to_string()]);
    let defined_off = || Observed {
        notice: true,
        header: false,
        existing: known("rs139263261"),
        clinical: true,
        flag_printed: "--var_synonyms 0".to_string(),
    };
    let config = dir.path().join("flags.ini");
    std::fs::write(&config, "var_synonyms 0\n").unwrap();
    assert_eq!(
        run(&["--config", config.to_str().unwrap()], &[]),
        defined_off(),
        "a 0 from the config file"
    );
    assert_eq!(
        run(&[], &[("VEP_VAR_SYNONYMS", "0")]),
        defined_off(),
        "a 0 from the environment"
    );
    let ini_dir = dir.path().join("ini");
    std::fs::create_dir(&ini_dir).unwrap();
    std::fs::write(ini_dir.join("vep.ini"), "var_synonyms 0\n").unwrap();
    assert_eq!(
        run(&["--dir", ini_dir.to_str().unwrap()], &[]),
        defined_off(),
        "a 0 from <dir>/vep.ini"
    );
    std::fs::write(&config, "var_synonyms \"\"\n").unwrap();
    assert_eq!(
        run(&["--config", config.to_str().unwrap()], &[]),
        Observed {
            flag_printed: "--var_synonyms".to_string(),
            ..defined_off()
        },
        "a quoted empty value from the config file: defined, printed with no value"
    );
    std::fs::write(&config, "var_synonyms 1\n").unwrap();
    assert_eq!(
        run(&["--config", config.to_str().unwrap()], &[]),
        Observed {
            header: true,
            flag_printed: "--var_synonyms".to_string(),
            ..defined_off()
        },
        "a 1 from the config file is the flag"
    );
    assert_eq!(
        run(&[], &[]),
        Observed {
            notice: false,
            header: false,
            existing: known("-"),
            clinical: false,
            flag_printed: String::new(),
        },
        "no flag from any source"
    );
}

/// A run that cannot open its FASTA or its input has printed the `--var_synonyms`
/// notice by then: `init` prints it (Runner.pm:117-120) before `fasta_db` (:134)
/// and the parser's open (:136), so the notice is the line before the error, the
/// error is the last line, and no output file is created.
#[test]
fn var_synonyms_notice_precedes_a_failed_fasta_or_input_open() {
    let colocated =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-colocated");
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let missing_fasta = dir.path().join("absent.fa");
    let missing_input = dir.path().join("absent.vcf");
    for (args, extra, error) in [
        (
            args_for(&colocated),
            vec!["--var_synonyms", "--fasta", missing_fasta.to_str().unwrap()],
            format!(
                "ERROR: [fasta] Specified FASTA file/directory {} not found",
                missing_fasta.display()
            ),
        ),
        (
            with_input(args_for(&colocated), &missing_input),
            vec!["--var_synonyms"],
            format!(
                "ERROR: [input] File \"{}\" does not exist",
                missing_input.display()
            ),
        ),
    ] {
        let refused = run_with(&args, &extra, &out);
        let stderr = text(&refused.stderr);
        assert_eq!(refused.status.code(), Some(1), "{extra:?}:\n{stderr}");
        let lines: Vec<&str> = stderr.lines().collect();
        let notice_at = lines
            .iter()
            .position(|l| *l == VAR_SYNONYMS_NOTICE)
            .unwrap_or_else(|| panic!("{extra:?}: no var_synonyms notice:\n{stderr}"));
        assert!(
            lines
                .get(notice_at + 1)
                .is_some_and(|l| l.starts_with(&error)),
            "{extra:?}: the notice is not the line before the error:\n{stderr}"
        );
        assert_eq!(notice_at + 2, lines.len(), "{extra:?}:\n{stderr}");
        assert!(refused.stdout.is_empty(), "{extra:?}");
        assert!(!out.exists(), "{extra:?}: an output file was created");
    }
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

/// A plugin whose constructor fails is dropped with the warning `Failed to
/// instantiate plugin <Name>: <reason>` and the run completes, or is the plugin
/// component's error under `--safe` (Runner.pm:946-954: `throw($msg) if
/// $self->param('safe'); $self->warning_msg($msg); next;`), with no output
/// written.
#[test]
fn plugin_init_failure_warns_and_is_fatal_under_safe() {
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
            "LoFtool,/no/such/scores.txt",
        ],
        dir.path(),
    );
    assert_eq!(warned.status.code(), Some(0), "{}", text(&warned.stderr));
    let prefix =
        "WARNING: Failed to instantiate plugin LoFtool: scores file /no/such/scores.txt not found";
    assert!(
        text(&warned.stderr).lines().any(|l| l == prefix),
        "{}",
        text(&warned.stderr)
    );
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    assert!(warnings.lines().any(|l| l == prefix), "{warnings}");
    let output = std::fs::read_to_string(&out).unwrap();
    let rows = data_lines(&output);
    assert!(rows.len() > 30);
    assert!(
        rows.iter().all(|row| !row.contains("LoFtool")),
        "the plugin is dropped"
    );
    std::fs::remove_file(&out).unwrap();

    let fatal = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir_plugins",
            corpus.to_str().unwrap(),
            "--plugin",
            "LoFtool,/no/such/scores.txt",
            "--safe",
        ],
        dir.path(),
    );
    assert_eq!(fatal.status.code(), Some(1));
    let stderr = text(&fatal.stderr);
    let errors: Vec<&str> = stderr
        .lines()
        .filter(|l| l.starts_with("ERROR: "))
        .collect();
    assert_eq!(errors.len(), 1, "{stderr}");
    assert_eq!(
        errors[0],
        prefix.replacen("WARNING: ", "ERROR: [plugin] ", 1)
    );
    assert!(
        !out.exists(),
        "no output is written when --safe refuses a plugin"
    );
}

/// UTRAnnotator reads the 5' UTR and coding sequences from the reference, so a
/// run under `--no_fasta` refuses it at load as a plugin whose constructor
/// fails: the warning `Failed to instantiate plugin UTRAnnotator: <reason>` on
/// standard error and in the warnings file, a complete run whose header and
/// rows name none of the plugin's five fields, and the plugin component's
/// error under `--safe` with no output written. The same run with the corpus's
/// reference declares the five fields. Ensembl VEP's module prints a notice
/// and runs on (UTRAnnotator.pm:88), the difference
/// `docs/plugins/transcript-model.md` states.
#[test]
fn utrannotator_without_a_reference_warns_and_is_fatal_under_safe() {
    const FIELDS: [&str; 5] = [
        "5UTR_consequence",
        "5UTR_annotation",
        "Existing_uORFs",
        "Existing_OutOfFrame_oORFs",
        "Existing_InFrame_oORFs",
    ];
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let corpus = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-plugin-utr-annotator");
    let args = args_for(&corpus);

    let warned = run_with(&args, &["--no_fasta", "--plugin", "UTRAnnotator"], &out);
    assert_eq!(warned.status.code(), Some(0), "{}", text(&warned.stderr));
    let prefix = "WARNING: Failed to instantiate plugin UTRAnnotator: a reference FASTA is needed; pass --fasta";
    assert!(
        text(&warned.stderr).lines().any(|l| l == prefix),
        "{}",
        text(&warned.stderr)
    );
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    assert!(warnings.lines().any(|l| l == prefix), "{warnings}");
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(!data_lines(&output).is_empty());
    for field in FIELDS {
        assert!(
            output.lines().all(|line| !line.contains(field)),
            "{field} is named though the plugin is dropped"
        );
    }
    std::fs::remove_file(&out).unwrap();

    let fatal = run_with(
        &args,
        &["--no_fasta", "--plugin", "UTRAnnotator", "--safe"],
        &out,
    );
    assert_eq!(fatal.status.code(), Some(1));
    let stderr = text(&fatal.stderr);
    let errors: Vec<&str> = stderr
        .lines()
        .filter(|l| l.starts_with("ERROR: "))
        .collect();
    assert_eq!(errors.len(), 1, "{stderr}");
    assert_eq!(
        errors[0],
        prefix.replacen("WARNING: ", "ERROR: [plugin] ", 1)
    );
    assert!(
        !out.exists(),
        "no output is written when --safe refuses a plugin"
    );

    let fasta = corpus.join("reference.fa.gz");
    let loaded = run_with(
        &args,
        &[
            "--fasta",
            fasta.to_str().unwrap(),
            "--plugin",
            "UTRAnnotator",
        ],
        &out,
    );
    assert_eq!(loaded.status.code(), Some(0), "{}", text(&loaded.stderr));
    let output = std::fs::read_to_string(&out).unwrap();
    for field in FIELDS {
        assert!(
            output
                .lines()
                .any(|line| line.starts_with(&format!("## {field} : "))),
            "{field} is declared when the plugin loads"
        );
    }
}

/// A built-in the run refuses for what the run lacks never reaches its
/// `init`, so its refusal is the one warning the run writes, on standard error
/// and in the warnings file: no notice of a `file=` it was not given and none
/// of a parameter it does not read. Under `--safe` the refusal is the run's
/// one error and no warning precedes it. UTRAnnotator and MaxEntScan need the
/// reference, HGVSIntronOffset the HGVS notations.
#[test]
fn a_plugin_refused_for_what_the_run_lacks_writes_its_refusal_alone() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let corpus = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-plugin-utr-annotator");
    let args = args_for(&corpus);
    let no_reference = "a reference FASTA is needed; pass --fasta";
    let no_hgvs = "the HGVS notations are needed; pass --hgvs";
    let cases: [(&[&str], &str, &str); 5] = [
        (&["--no_fasta"], "UTRAnnotator", no_reference),
        (&["--no_fasta"], "UTRAnnotator,stray", no_reference),
        (&["--no_fasta"], "MaxEntScan", no_reference),
        (&["--no_fasta"], "MaxEntScan,/no/such/dir", no_reference),
        (&[], "HGVSIntronOffset,stray", no_hgvs),
    ];
    for (flags, plugin, need) in cases {
        let name = plugin.split(',').next().unwrap();
        let refusal = format!("Failed to instantiate plugin {name}: {need}");
        let warning = format!("WARNING: {refusal}");
        let mut extra = flags.to_vec();
        extra.extend(["--plugin", plugin]);

        let warned = run_with(&args, &extra, &out);
        let stderr = text(&warned.stderr);
        assert_eq!(warned.status.code(), Some(0), "{plugin}: {stderr}");
        assert_eq!(warning_lines(&stderr), [warning.as_str()], "{plugin}");
        let warnings_file = dir.path().join("out.txt_warnings.txt");
        let warnings = std::fs::read_to_string(&warnings_file).unwrap();
        assert_eq!(warnings, format!("{warning}\n"), "{plugin}");
        // A run that writes no warning leaves an earlier run's warnings file in place.
        std::fs::remove_file(&warnings_file).unwrap();
        std::fs::remove_file(&out).unwrap();

        extra.push("--safe");
        let fatal = run_with(&args, &extra, &out);
        let stderr = text(&fatal.stderr);
        assert_eq!(fatal.status.code(), Some(1), "{plugin}: {stderr}");
        assert!(warning_lines(&stderr).is_empty(), "{plugin}: {stderr}");
        let errors: Vec<&str> = stderr
            .lines()
            .filter(|l| l.starts_with("ERROR: "))
            .collect();
        assert_eq!(errors, [format!("ERROR: [plugin] {refusal}")], "{plugin}");
        assert!(!out.exists(), "{plugin}: --safe writes no output");
    }
}

/// `LoF`, the name of LoFTEE's Ensembl VEP module, loads the LoFTEE built-in:
/// the run under `--plugin LoF` writes the `## LoF` header lines and the `LoF`
/// keys of the run under `--plugin LoFTEE`, byte for byte, so an Ensembl VEP
/// argument string loads it unchanged, and a value asking for the
/// splice-prediction layer is refused at initialization. With the GERP rule on
/// and no GERP word, the plugin reads the bigWig Ensembl VEP's plugin reads,
/// `loftee_path` followed by `gerp_conservation_scores.homo_sapiens.GRCh38.bw`
/// (`LoF.pm` line 130), so `GERP_DIST` reaches `LoF_info`; the bare word is
/// refused at initialization naming the file under `/vep/loftee/` (`LoF.pm`
/// line 83), as Ensembl VEP refuses it for its own files there: a warning and
/// no `LoF` field, the plugin component's error under `--safe`.
#[test]
fn lof_names_the_loftee_plugin() {
    let dir = tempfile::tempdir().unwrap();
    let mut outputs = Vec::new();
    for name in ["LoFTEE", "LoF"] {
        let out = dir.path().join(format!("{name}.txt"));
        let run = run_control(
            &[
                "-o",
                out.to_str().unwrap(),
                "--plugin",
                &format!("{name},min_intron_size:20,use_gerp_end_trunc:0"),
            ],
            dir.path(),
        );
        assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
        let output = std::fs::read_to_string(&out).unwrap();
        assert!(output.contains("\n## LoF : "), "{name}: no LoF header line");
        assert!(
            data_lines(&output).iter().any(|row| row.contains(";LoF=")),
            "{name}: no row carries a LoF call"
        );
        // The run's own lines (the time, the command line naming the plugin
        // and the output file) differ between the two runs by construction.
        let shared: Vec<&str> = output
            .lines()
            .filter(|l| {
                !l.starts_with("## Output produced at ") && !l.starts_with("## VEP command-line:")
            })
            .collect();
        outputs.push(shared.join("\n"));
    }
    assert_eq!(outputs[0], outputs[1]);

    let out = dir.path().join("refused.txt");
    let refused = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir_plugins",
            dir.path().to_str().unwrap(),
            "--plugin",
            "LoF,run_splice_predictions:1",
            "--safe",
        ],
        dir.path(),
    );
    assert_eq!(refused.status.code(), Some(1));
    assert!(
        text(&refused.stderr).lines().any(|l| l == "ERROR: [plugin] Failed to instantiate plugin LoF: LoFTEE: run_splice_predictions=1 asks for the splice-prediction layer, which vep-rs does not implement"),
        "{}",
        text(&refused.stderr)
    );
    assert!(!out.exists());

    let loftee_path = dir.path().join("loftee");
    std::fs::create_dir(&loftee_path).unwrap();
    let bigwig = "gerp_conservation_scores.homo_sapiens.GRCh38.bw";
    std::fs::copy(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116/GRCh38-plugin-loftee/plugin_data")
            .join(bigwig),
        loftee_path.join(bigwig),
    )
    .unwrap();
    let out = dir.path().join("derived.txt");
    let derived = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--plugin",
            &format!("LoF,loftee_path:{}/", loftee_path.display()),
        ],
        dir.path(),
    );
    assert_eq!(derived.status.code(), Some(0), "{}", text(&derived.stderr));
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(
        data_lines(&output)
            .iter()
            .any(|row| row.contains(";LoF=") && row.contains("GERP_DIST:")),
        "no LoF row read the derived GERP bigWig:\n{output}"
    );

    // The bare word reads the plugin's own default directory, so its refusal is
    // asserted only on a host without one.
    if std::path::Path::new("/vep/loftee/").exists() {
        return;
    }
    let refusal = "Failed to instantiate plugin LoF: LoFTEE: GERP data file not found: \
                   '/vep/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw'; \
                   Ensembl VEP's LoF plugin derives this GERP source from loftee_path when no \
                   gerp_bigwig is given; pass the GERP file, or use_gerp_end_trunc=0 for the \
                   percentile rule";
    let out = dir.path().join("bare.txt");
    let bare = run_control(
        &["-o", out.to_str().unwrap(), "--plugin", "LoF"],
        dir.path(),
    );
    assert_eq!(bare.status.code(), Some(0), "{}", text(&bare.stderr));
    assert!(
        text(&bare.stderr)
            .lines()
            .any(|l| l == format!("WARNING: {refusal}")),
        "{}",
        text(&bare.stderr)
    );
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(
        !output.contains("\n## LoF : ")
            && data_lines(&output).iter().all(|row| !row.contains("LoF=")),
        "the plugin is dropped, so no LoF header or field is written:\n{output}"
    );

    let out = dir.path().join("bare-safe.txt");
    let bare_safe = run_control(
        &["-o", out.to_str().unwrap(), "--plugin", "LoF", "--safe"],
        dir.path(),
    );
    assert_eq!(bare_safe.status.code(), Some(1));
    assert!(
        text(&bare_safe.stderr)
            .lines()
            .any(|l| l == format!("ERROR: [plugin] {refusal}")),
        "{}",
        text(&bare_safe.stderr)
    );
    assert!(!out.exists());
}

/// LoFTEE's `conservation_file` names the database Ensembl VEP's plugin reads
/// for the PhyloCSF items of `LoF_flags` and `LoF_info` (`LoF.pm` lines
/// 274-285), so a value naming no file is refused at initialization, and so
/// is the word `mysql`, under which the plugin opens a database connection
/// (line 122): a warning and no `LoF` field, the plugin component's error
/// under `--safe`. Without the parameter the corpus's three stop-gained rows
/// carry their `LoF` fields, and under the off word `false` (line 274) the
/// output is the one written without the parameter, byte for byte.
#[test]
fn loftee_refuses_a_conservation_file_it_cannot_read_and_runs_unchanged_without_one() {
    let dir = tempfile::tempdir().unwrap();
    let run = |name: &str, plugin: &str, extra: &[&str]| {
        let out = dir.path().join(format!("{name}.txt"));
        let mut args = vec!["-o", out.to_str().unwrap(), "--plugin", plugin];
        args.extend_from_slice(extra);
        (run_control(&args, dir.path()), out)
    };

    let (without, out) = run("without", "LoF,use_gerp_end_trunc:0", &[]);
    assert_eq!(without.status.code(), Some(0), "{}", text(&without.stderr));
    let without = std::fs::read_to_string(&out).unwrap();
    let lof_fields: Vec<String> = data_lines(&without)
        .iter()
        .filter_map(|row| {
            let fields: Vec<&str> = row
                .rsplit('\t')
                .next()?
                .split(';')
                .filter(|field| field.starts_with("LoF"))
                .collect();
            (!fields.is_empty()).then(|| fields.join(";"))
        })
        .collect();
    assert_eq!(
        lof_fields,
        [
            "LoF=LC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;\
             LoF_info=PERCENTILE:0.993589743589744",
            "LoF=LC;LoF_filter=END_TRUNC;LoF_info=PERCENTILE:0.993019197207679",
            "LoF=LC;LoF_filter=END_TRUNC;LoF_info=PERCENTILE:0.993019197207679",
        ]
    );

    let (off, out) = run(
        "off",
        "LoF,use_gerp_end_trunc:0,conservation_file:false",
        &[],
    );
    assert_eq!(off.status.code(), Some(0), "{}", text(&off.stderr));
    assert_eq!(
        without_run_lines(&std::fs::read_to_string(&out).unwrap()),
        without_run_lines(&without)
    );

    let missing = dir.path().join("phylocsf.synthetic.sqlite");
    let missing = missing.to_str().unwrap();
    for (name, plugin, refusal) in [
        (
            "missing",
            format!("LoF,use_gerp_end_trunc:0,conservation_file:{missing}"),
            format!(
                "Failed to instantiate plugin LoF: LoFTEE: conservation_file database not \
                 found: '{missing}'"
            ),
        ),
        (
            "connection",
            "LoF,use_gerp_end_trunc:0,conservation_file:mysql".to_string(),
            "Failed to instantiate plugin LoF: LoFTEE: conservation_file=mysql asks for a \
             database connection, which vep-rs does not open: pass the database file's path \
             under conservation_file"
                .to_string(),
        ),
    ] {
        let (dropped, out) = run(&format!("{name}-dropped"), &plugin, &[]);
        assert_eq!(dropped.status.code(), Some(0), "{}", text(&dropped.stderr));
        assert!(
            text(&dropped.stderr)
                .lines()
                .any(|l| l == format!("WARNING: {refusal}")),
            "{}",
            text(&dropped.stderr)
        );
        let output = std::fs::read_to_string(&out).unwrap();
        assert!(
            !output.contains("\n## LoF : ")
                && data_lines(&output).iter().all(|row| !row.contains("LoF=")),
            "the plugin is dropped, so no LoF header or field is written:\n{output}"
        );

        let (refused, out) = run(&format!("{name}-refused"), &plugin, &["--safe"]);
        assert_eq!(refused.status.code(), Some(1));
        assert!(
            text(&refused.stderr)
                .lines()
                .any(|l| l == format!("ERROR: [plugin] {refusal}")),
            "{}",
            text(&refused.stderr)
        );
        assert!(!out.exists());
    }
    assert!(
        !std::path::Path::new(missing).exists(),
        "no database file is created"
    );
}

/// Ensembl VEP's LoF plugin stores every `key:value` word and reads the keys
/// its own branch names, printing nothing for a word it never reads (`LoF.pm`
/// lines 67-72 on `master` and on `grch38`). `use_gerp_end_trunc`,
/// `gerp_end_trunc_cutoff` and `filter_position` are read on `grch38` alone,
/// `skip_lof_info` and `tabix_path` on `master` alone. A run that passes the
/// other assembly's keys therefore writes the output and the warnings of the
/// run without them: the GRCh37 frameshift at `21:10970009-10970012` keeps
/// its GERP items and the GRCh38 stop-gained rows keep `LoF_info`. Each such
/// key is named in one status line on standard error.
#[test]
fn loftee_keys_of_the_other_assembly_change_no_row_and_raise_no_warning() {
    let dir = tempfile::tempdir().unwrap();
    let warnings_of =
        |out: &Path| std::fs::read_to_string(format!("{}_warnings.txt", out.display())).ok();
    let lof_fields = |output: &str, location: Option<&str>| -> Vec<String> {
        data_lines(output)
            .iter()
            .map(|row| row.split('\t').collect::<Vec<_>>())
            .filter(|cells| location.is_none_or(|location| cells[1] == location))
            .filter_map(|cells| {
                let fields: Vec<&str> = cells
                    .last()?
                    .split(';')
                    .filter(|field| field.starts_with("LoF"))
                    .collect();
                (!fields.is_empty()).then(|| format!("{} {}", cells[4], fields.join(";")))
            })
            .collect()
    };
    // The two runs agree in every output line but the two that name the run,
    // in the warnings file and in the warnings on standard error, and the run
    // with the keys names each in one status line.
    let same_run = |without: (&Output, &Path), with: (&Output, &Path), notes: &[String]| {
        for (run, _) in [without, with] {
            assert_eq!(run.status.code(), Some(0), "{}", text(&run.stderr));
        }
        let (without_stderr, with_stderr) = (text(&without.0.stderr), text(&with.0.stderr));
        assert_eq!(
            without_run_lines(&std::fs::read_to_string(with.1).unwrap()),
            without_run_lines(&std::fs::read_to_string(without.1).unwrap())
        );
        assert_eq!(warnings_of(with.1), warnings_of(without.1));
        assert_eq!(
            warning_lines(&with_stderr),
            warning_lines(&without_stderr),
            "{with_stderr}"
        );
        for note in notes {
            let named = |stderr: &str| stderr.lines().filter(|l| l.contains(note)).count();
            assert_eq!(named(&with_stderr), 1, "{note}\n{with_stderr}");
            assert_eq!(named(&without_stderr), 0, "{note}\n{without_stderr}");
        }
    };

    let gerp = gerp_tabix_tsv(dir.path(), &[("21", 10970009, "2.5")]);
    let grch37 = format!("LoF,gerp_file:{},run_splice_predictions:0", gerp.display());
    let run = |name: &str, plugin: &str| {
        let out = dir.path().join(format!("{name}.txt"));
        (
            run_with(&hgvs_corpus_args(), &["--plugin", plugin], &out),
            out,
        )
    };
    let (without, without_out) = run("grch37_without", &grch37);
    let (with, with_out) = run(
        "grch37_with",
        &format!("{grch37},gerp_end_trunc_cutoff:-1000,use_gerp_end_trunc:0,filter_position:0.9"),
    );
    let notes = [
        "gerp_end_trunc_cutoff",
        "use_gerp_end_trunc",
        "filter_position",
    ]
    .map(|key| {
        format!(
            "LoFTEE: ignoring '{key}' on GRCh37: Ensembl VEP's LoF plugin for that assembly \
             never reads it and judges END_TRUNC by the GERP rule at its own cutoff of 180"
        )
    });
    same_run((&without, &without_out), (&with, &with_out), &notes);
    assert_eq!(warnings_of(&with_out), None);
    assert_eq!(
        lof_fields(
            &std::fs::read_to_string(&with_out).unwrap(),
            Some("21:10970009-10970012")
        ),
        [
            "ENST00000298232 LoF=HC;LoF_info=GERP_DIST:2.5,BP_DIST:1466,\
             PERCENTILE:0.0742821473158552,DIST_FROM_LAST_EXON:1331,50_BP_RULE:PASS",
            "ENST00000342420 LoF=HC;LoF_info=GERP_DIST:2.5,BP_DIST:1407,\
             PERCENTILE:0.0771725032425422,DIST_FROM_LAST_EXON:1272,50_BP_RULE:PASS",
            "ENST00000361285 LoF=HC;LoF_info=GERP_DIST:2.5,BP_DIST:1519,\
             PERCENTILE:0.0718599033816425,DIST_FROM_LAST_EXON:1384,50_BP_RULE:PASS",
        ]
    );

    let run = |name: &str, plugin: &str| {
        let out = dir.path().join(format!("{name}.txt"));
        let args = ["-o", out.to_str().unwrap(), "--plugin", plugin];
        (run_control(&args, dir.path()), out)
    };
    let grch38 = "LoF,use_gerp_end_trunc:0";
    let (without, without_out) = run("grch38_without", grch38);
    let (with, with_out) = run(
        "grch38_with",
        &format!("{grch38},skip_lof_info:1,tabix_path:/usr/local/bin/tabix"),
    );
    let notes = ["skip_lof_info", "tabix_path"].map(|key| {
        format!("LoFTEE: ignoring '{key}': Ensembl VEP's LoF plugin reads it on GRCh37 alone")
    });
    same_run((&without, &without_out), (&with, &with_out), &notes);
    let rows = lof_fields(&std::fs::read_to_string(&with_out).unwrap(), None);
    let rows: Vec<&str> = rows
        .iter()
        .map(|row| row.split_once(' ').unwrap().1)
        .collect();
    assert_eq!(
        rows,
        [
            "LoF=LC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;\
             LoF_info=PERCENTILE:0.993589743589744",
            "LoF=LC;LoF_filter=END_TRUNC;LoF_info=PERCENTILE:0.993019197207679",
            "LoF=LC;LoF_filter=END_TRUNC;LoF_info=PERCENTILE:0.993019197207679",
        ]
    );
}

/// A made-up tabix-indexed GERP TSV holding `rows`, under `dir`: chromosome,
/// position and score, the score in the third field, where the `master`
/// branch of LoFTEE reads it (`gerp_dist.pl` lines 96 and 121).
fn gerp_tabix_tsv(dir: &Path, rows: &[(&str, usize, &str)]) -> PathBuf {
    use noodles::csi::binning_index::index::header::Builder as HeaderBuilder;

    let lines: Vec<String> = rows
        .iter()
        .map(|(chromosome, position, score)| format!("{chromosome}\t{position}\t{score}"))
        .collect();
    let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
    common::indexed_lines(
        dir,
        "GERP_scores.synthetic.txt.gz",
        HeaderBuilder::default()
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(None)
            .build(),
        &["#chrom\tpos\tGERP"],
        &lines,
        |fields| {
            let base = fields[1].parse::<usize>().unwrap();
            (base, base)
        },
    )
}

/// master's `LoF.pm`, the LoFTEE plugin of a GRCh37 run, runs its
/// splice-prediction layer unless `run_splice_predictions` is `0` (lines 98
/// and 186), a layer vep-rs does not implement, so a GRCh37 run that leaves
/// the key out is refused at initialization: a warning and no `LoF` field, the
/// plugin component's error under `--safe`. Under `run_splice_predictions:0`
/// the GRCh37-hgvs corpus's rows carry their `LoF` fields, the frameshift at
/// `21:10970009-10970012` with the GERP total of the named file's one row.
#[test]
fn loftee_refuses_a_grch37_run_that_leaves_splice_predictions_on() {
    let dir = tempfile::tempdir().unwrap();
    let gerp = gerp_tabix_tsv(dir.path(), &[("21", 10970009, "2.5")]);
    let run = |name: &str, plugin: &str, extra: &[&str]| {
        let out = dir.path().join(format!("{name}.txt"));
        let mut args = vec!["--plugin", plugin];
        args.extend_from_slice(extra);
        (run_with(&hgvs_corpus_args(), &args, &out), out)
    };

    let left_on = format!("LoF,gerp_file:{}", gerp.display());
    let refusal = "Failed to instantiate plugin LoF: LoFTEE: this GRCh37 run asks for the \
                   splice-prediction layer, which vep-rs does not implement: Ensembl VEP's LoF \
                   plugin runs it on that assembly unless run_splice_predictions is 0; pass \
                   run_splice_predictions:0";
    let (dropped, out) = run("dropped", &left_on, &[]);
    assert_eq!(dropped.status.code(), Some(0), "{}", text(&dropped.stderr));
    assert!(
        text(&dropped.stderr)
            .lines()
            .any(|l| l == format!("WARNING: {refusal}")),
        "{}",
        text(&dropped.stderr)
    );
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(
        !output.contains("\n## LoF : ")
            && data_lines(&output).iter().all(|row| !row.contains("LoF=")),
        "the plugin is dropped, so no LoF header or field is written:\n{output}"
    );

    let (refused, out) = run("refused", &left_on, &["--safe"]);
    assert_eq!(refused.status.code(), Some(1));
    assert!(
        text(&refused.stderr)
            .lines()
            .any(|l| l == format!("ERROR: [plugin] {refusal}")),
        "{}",
        text(&refused.stderr)
    );
    assert!(!out.exists());

    let (off, out) = run("off", &format!("{left_on},run_splice_predictions:0"), &[]);
    assert_eq!(off.status.code(), Some(0), "{}", text(&off.stderr));
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(output.contains("\n## LoF : "), "no LoF header line");
    let lof_rows: Vec<Vec<&str>> = data_lines(&output)
        .iter()
        .map(|row| row.split('\t').collect::<Vec<_>>())
        .filter(|cells| cells.last().is_some_and(|extra| extra.contains("LoF=")))
        .collect();
    assert_eq!(lof_rows.len(), 529);
    let frameshift: Vec<(&str, String)> = lof_rows
        .iter()
        .filter(|cells| cells[1] == "21:10970009-10970012")
        .map(|cells| {
            let fields: Vec<&str> = cells
                .last()
                .unwrap()
                .split(';')
                .filter(|field| field.starts_with("LoF"))
                .collect();
            (cells[4], fields.join(";"))
        })
        .collect();
    let row = |transcript: &'static str, bp_dist: u32, percentile: &str, from_last_exon: u32| {
        (
            transcript,
            format!(
                "LoF=HC;LoF_info=GERP_DIST:2.5,BP_DIST:{bp_dist},PERCENTILE:{percentile},\
                 DIST_FROM_LAST_EXON:{from_last_exon},50_BP_RULE:PASS"
            ),
        )
    };
    assert_eq!(
        frameshift,
        [
            row("ENST00000298232", 1466, "0.0742821473158552", 1331),
            row("ENST00000342420", 1407, "0.0771725032425422", 1272),
            row("ENST00000361285", 1519, "0.0718599033816425", 1384),
        ]
    );
}

/// master's `LoF.pm` opens a database connection for its GERP scores when
/// `gerp_file` is the word `mysql` (lines 134-136), which vep-rs does not
/// open, so a GRCh37 run naming the word is refused at initialization under a
/// line that says so, whatever file another GERP key names: a warning and no
/// `LoF` field, the plugin component's error under `--safe`.
#[test]
fn loftee_refuses_a_gerp_database_connection_on_grch37() {
    let dir = tempfile::tempdir().unwrap();
    let gerp = gerp_tabix_tsv(dir.path(), &[("21", 10970009, "2.5")]);
    let run = |name: &str, plugin: &str, extra: &[&str]| {
        let out = dir.path().join(format!("{name}.txt"));
        let mut args = vec!["--plugin", plugin];
        args.extend_from_slice(extra);
        (run_with(&hgvs_corpus_args(), &args, &out), out)
    };
    let refusal = "Failed to instantiate plugin LoF: LoFTEE: gerp_file=mysql asks for a \
                   database connection, which vep-rs does not open: pass the GERP table's path \
                   under gerp_file or gerp_tabix";
    for (name, plugin) in [
        (
            "word",
            "LoF,gerp_file:mysql,run_splice_predictions:0".to_string(),
        ),
        (
            "beside-a-file",
            format!(
                "LoF,gerp_tabix={},gerp_file:mysql,run_splice_predictions:0",
                gerp.display()
            ),
        ),
    ] {
        let (dropped, out) = run(&format!("{name}-dropped"), &plugin, &[]);
        assert_eq!(dropped.status.code(), Some(0), "{}", text(&dropped.stderr));
        assert!(
            text(&dropped.stderr)
                .lines()
                .any(|l| l == format!("WARNING: {refusal}")),
            "{}",
            text(&dropped.stderr)
        );
        let output = std::fs::read_to_string(&out).unwrap();
        assert!(
            !output.contains("\n## LoF : ")
                && data_lines(&output).iter().all(|row| !row.contains("LoF=")),
            "the plugin is dropped, so no LoF header or field is written:\n{output}"
        );

        let (refused, out) = run(&format!("{name}-refused"), &plugin, &["--safe"]);
        assert_eq!(refused.status.code(), Some(1));
        assert!(
            text(&refused.stderr)
                .lines()
                .any(|l| l == format!("ERROR: [plugin] {refusal}")),
            "{}",
            text(&refused.stderr)
        );
        assert!(!out.exists());
    }
}

/// A `--plugin` file that exists but does not open as a shared library is
/// dropped with `Failed to compile plugin <Name>: <reason>`, as Ensembl VEP drops
/// a module `use` cannot load (Runner.pm:929-939), and the run completes; under
/// `--safe` it is the plugin component's error. The warnings follow the
/// `--plugin` arguments' order, as Ensembl VEP's one loop over them writes
/// them (`:925`): a failing library, then a failing built-in, then a name
/// nothing resolves.
#[test]
fn plugin_load_failures_warn_in_argument_order_and_a_bad_library_is_fatal_under_safe() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let corpus = run_control_dir();
    let not_a_library = corpus.join("manifest.json");
    let not_a_library = not_a_library.to_str().unwrap();
    let warned = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir_plugins",
            corpus.to_str().unwrap(),
            "--plugin",
            not_a_library,
            "--plugin",
            "LoFtool,/no/such/scores.txt",
            "--plugin",
            "NoSuchPlugin",
        ],
        dir.path(),
    );
    assert_eq!(warned.status.code(), Some(0), "{}", text(&warned.stderr));
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    let plugin_lines: Vec<&str> = warnings
        .lines()
        .filter(|l| l.starts_with("WARNING: Failed to "))
        .collect();
    assert_eq!(plugin_lines.len(), 3, "{warnings}");
    assert!(
        plugin_lines[0].starts_with(&format!(
            "WARNING: Failed to compile plugin {not_a_library}: {not_a_library}: "
        )),
        "{warnings}"
    );
    assert_eq!(
        plugin_lines[1],
        "WARNING: Failed to instantiate plugin LoFtool: scores file /no/such/scores.txt not found",
        "{warnings}"
    );
    assert_eq!(
        plugin_lines[2],
        format!(
            "WARNING: Failed to load plugin NoSuchPlugin: NoSuchPlugin is not a built-in plugin and no libNoSuchPlugin.so or libNoSuchPlugin.dylib is under {}",
            corpus.display()
        ),
        "{warnings}"
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
            not_a_library,
            "--safe",
        ],
        dir.path(),
    );
    assert_eq!(fatal.status.code(), Some(1));
    let stderr = text(&fatal.stderr);
    let errors: Vec<&str> = stderr
        .lines()
        .filter(|l| l.starts_with("ERROR: "))
        .collect();
    assert_eq!(errors.len(), 1, "{stderr}");
    assert!(
        errors[0].starts_with(&format!(
            "ERROR: [plugin] Failed to compile plugin {not_a_library}: {not_a_library}: "
        )),
        "{}",
        errors[0]
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

/// `--dir_plugins` is expanded as `glob` expands it (Runner.pm:908): over `p1/`,
/// `p2/` and `f.txt`, `<tmp>/*` makes the two directories the search path, last
/// match first, and the file is `Plugins directory '<tmp>/f.txt' not found.` in
/// the warnings file when `--plugin` names a plugin (`:911-916`), nothing when
/// none is named, and the plugin component's error under `--safe`.
#[test]
fn dir_plugins_non_directory_match_warns_only_with_a_plugin() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let plugins = dir.path().join("plugins");
    for name in ["p1", "p2"] {
        std::fs::create_dir_all(plugins.join(name)).unwrap();
    }
    std::fs::write(plugins.join("f.txt"), "").unwrap();
    let pattern = format!("{}/*", plugins.display());
    let not_found = format!(
        "WARNING: Plugins directory '{}/f.txt' not found.",
        plugins.display()
    );

    let with_plugin = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir_plugins",
            &pattern,
            "--plugin",
            "echo_plugin",
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
    let directory_lines: Vec<&str> = warnings
        .lines()
        .filter(|l| l.contains("Plugins directory"))
        .collect();
    assert_eq!(directory_lines, [not_found.as_str()], "{warnings}");
    assert!(
        warnings.contains(&format!(
            "no libecho_plugin.so or libecho_plugin.dylib is under {0}/p2, {0}/p1",
            plugins.display()
        )),
        "{warnings}"
    );

    let without_plugin = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--dir_plugins",
            &pattern,
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
            &pattern,
            "--plugin",
            "echo_plugin",
            "--safe",
        ],
        dir.path(),
    );
    assert_eq!(safe.status.code(), Some(1));
    assert!(
        text(&safe.stderr)
            .lines()
            .any(|l| l == format!("ERROR: [plugin] {}", &not_found["WARNING: ".len()..])),
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
    let refused = vep()
        .args(without_force_overwrite(args_for(&run_control_dir())))
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

/// `--vcf_info_field` takes a VCF INFO key alone (`[A-Za-z_][0-9A-Za-z_.]*`, or
/// the reserved `1000G`): the empty string, `0` and a name with a `;` are argument
/// syntax errors (exit 2) naming the flag. Ensembl VEP 116.3 accepts them, prints
/// the given name in the `##INFO` line (OutputFactory/VCF.pm:268) and keys the
/// entry `CSQ` for the two Perl reads as false (`:219`, `:312`). A dotted key,
/// which the grammar admits, and `1000G` are accepted.
#[test]
fn vcf_info_field_refuses_a_value_that_is_not_an_info_key() {
    let run = |value: &str| {
        vep()
            .args([
                "--vcf_info_field",
                value,
                "-i",
                "/nonexistent/input.vcf",
                "--quiet",
            ])
            .stdin(Stdio::null())
            .output()
            .unwrap()
    };
    for value in ["", "0", "A;B", "1KG", "a b"] {
        let output = run(value);
        let stderr = text(&output.stderr);
        assert_eq!(output.status.code(), Some(2), "{value:?}: {stderr}");
        assert!(
            stderr.contains(&format!(
                "error: invalid value '{value}' for '--vcf_info_field <VCF_INFO_FIELD>': a VCF \
                 INFO key is expected: a letter or underscore, then letters, digits, underscores \
                 or dots, or 1000G"
            )),
            "{value:?}: {stderr}"
        );
    }
    for value in ["A.N", "1000G"] {
        let output = run(value);
        let stderr = text(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{value}: {stderr}");
        assert!(stderr.starts_with("ERROR: ["), "{value}: {stderr}");
    }
}

/// A warnings file that `--warning_file` names and that cannot be created
/// fails the run with Ensembl VEP's text (BaseVEP.pm:894) under the output
/// component, exit 1, before the first record is read: no output is written
/// and no warning is echoed, with and without `--quiet`. Ensembl VEP 116.3
/// opens the file at its first warning and dies there, with the header alone
/// written and no warning echoed, exit 2; `docs/cli/compatibility.md` states
/// the difference.
#[test]
fn unwritable_warnings_file_fails_the_run() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let unwritten = "ERROR: [output] Could not write to warnings file /nonexistent-dir/w.txt";
    for quiet in [&["--quiet"][..], &[]] {
        let mut args = vec![
            "-o",
            out.to_str().unwrap(),
            "--warning_file",
            "/nonexistent-dir/w.txt",
        ];
        args.extend(quiet);
        let output = run_control(&args, dir.path());
        assert_eq!(output.status.code(), Some(1), "{quiet:?}");
        let stderr = text(&output.stderr);
        assert_eq!(
            stderr
                .lines()
                .filter(|line| line.starts_with("ERROR:"))
                .collect::<Vec<_>>(),
            [unwritten],
            "{quiet:?}: {stderr}"
        );
        assert_eq!(warning_lines(&stderr), Vec::<&str>::new(), "{quiet:?}");
        assert!(!out.exists(), "{quiet:?}: the refused run wrote an output");
    }
    let quiet = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--warning_file",
            "/nonexistent-dir/w.txt",
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(text(&quiet.stderr).trim_end(), unwritten);
}

/// A warnings file that fails at its first write ends the run after the batch
/// that warned, its rows written (here the corpus is one batch, so every row
/// of a run without the flag); under `--quiet` the error line is all that
/// reaches standard error, and without it each warning is echoed on standard
/// error before the error line. The file is Linux's `/dev/full`, a device,
/// which the run does not test before its first record.
#[cfg(target_os = "linux")]
#[test]
fn a_warnings_file_that_fails_at_its_first_write_ends_the_run_after_the_batch_that_warned() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let unwritten = "ERROR: [output] Could not write to warnings file /dev/full";
    let output = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--warning_file",
            "/dev/full",
            "--quiet",
        ],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(text(&output.stderr).trim_end(), unwritten);
    let row_count = data_lines(&std::fs::read_to_string(&out).unwrap()).len();
    let plain = dir.path().join("plain.txt");
    assert!(
        run_control(&["-o", plain.to_str().unwrap(), "--quiet"], dir.path())
            .status
            .success()
    );
    assert_eq!(
        row_count,
        data_lines(&std::fs::read_to_string(&plain).unwrap()).len()
    );

    let output = run_control(
        &["-o", out.to_str().unwrap(), "--warning_file", "/dev/full"],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(1));
    let stderr = text(&output.stderr);
    assert_eq!(
        warning_lines(&stderr),
        expected_warnings().lines().collect::<Vec<_>>(),
        "{stderr}"
    );
    assert!(stderr.trim_end().ends_with(unwritten), "{stderr}");
}

/// A Parquet run that holds a record out of order and cannot write its
/// warnings file ends on one error line in both pipeline modes at any worker
/// count. The warnings file's failure is the error when the batch that holds
/// the refused record, or a batch before it, drew a warning, wherever in the
/// batch the skipped line stands, as a `--warning_file` that fails at the
/// warning is raised after the batch that warned; the order refusal is the
/// error when the first skipped line stands in a later batch, whose warnings
/// a refused run never reports. Of 400 records in batches of 50, one is
/// skipped and one is out of order. The warnings file is Linux's `/dev/full`,
/// a device, which the run does not test before its first record.
#[cfg(all(feature = "parquet", target_os = "linux"))]
#[test]
fn a_refused_batch_that_could_not_write_a_warning_ends_on_the_warnings_file_in_each_mode() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.parquet");
    let warnings = Path::new("/dev/full");
    let unwritten = format!(
        "ERROR: [output] Could not write to warnings file {}",
        warnings.display()
    );
    let out_of_order_at_line_20 = "ERROR: [parquet] input is not sorted by chromosome and position at line 20 (21:5 after 21:119); sort the input or pass --no_check_variants_order";
    for (skipped, out_of_order, error) in [
        (3, 20, unwritten.as_str()),
        (40, 20, unwritten.as_str()),
        (3, 70, unwritten.as_str()),
        (70, 20, out_of_order_at_line_20),
    ] {
        let records: Vec<(&str, u64)> = (1..=400u64)
            .map(|line| match line {
                _ if line == skipped => ("Z", 1_000 + line),
                _ if line == out_of_order => ("21", 5),
                _ => ("21", 100 + line),
            })
            .collect();
        let input = intergenic_vcf(dir.path(), "in.vcf", &records);
        for flags in modes_and_worker_counts(3) {
            let mut extra = vec![
                "--buffer_size",
                "50",
                "--warning_file",
                warnings.to_str().unwrap(),
            ];
            extra.extend(flags);
            let output = parquet_run(&input, &extra, &out);
            let run =
                format!("line {skipped} skipped, line {out_of_order} out of order, {flags:?}");
            assert_eq!(output.status.code(), Some(1), "{run}");
            assert_eq!(text(&output.stderr).trim_end(), error, "{run}");
            assert!(!out.exists(), "{run}: the run leaves no output");
        }
    }
}

/// A text output that takes no bytes and a warnings file that cannot be
/// written end the run on the warnings file's error in both pipeline modes at
/// any worker count, whatever the output's buffer holds. One batch of 96
/// records, four of them skipped, writes 92 rows without its header, 6,064
/// bytes: with a buffer of 1,000 bytes a write of the batch fails, and with
/// one of 9,000 the buffer holds the whole output and the flush that completes
/// it is the first write to fail. The output and the warnings file are
/// Linux's `/dev/full`, a device, which the run does not test before its
/// first record.
#[cfg(target_os = "linux")]
#[test]
fn a_batch_that_could_not_write_a_warning_ends_on_the_warnings_file_at_any_output_buffer_size() {
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<(&str, u64)> = (1..=96u64)
        .map(|line| match line {
            3 | 22 | 28 | 40 => ("Z", 1_000 + line),
            _ => ("21", 100 + line),
        })
        .collect();
    let input = intergenic_vcf(dir.path(), "in.vcf", &records);
    let warnings = Path::new("/dev/full");
    let unwritten = format!(
        "ERROR: [output] Could not write to warnings file {}",
        warnings.display()
    );
    let mut args = corpus_args_with_input(&input);
    args.extend(
        [
            "--quiet",
            "--no_headers",
            "--warning_file",
            warnings.to_str().unwrap(),
        ]
        .map(str::to_string),
    );
    for buffer in ["1000", "9000"] {
        for flags in modes_and_worker_counts(3) {
            let mut extra = vec!["--output_buffer_size", buffer];
            extra.extend(flags);
            let output = run_with(&args, &extra, Path::new("/dev/full"));
            let run = format!("a buffer of {buffer} bytes, {flags:?}");
            assert_eq!(output.status.code(), Some(1), "{run}");
            assert_eq!(text(&output.stderr).trim_end(), unwritten, "{run}");
        }
    }
}

/// A comma-joined or trailing-comma value of a clap-enumerated `%VALID` flag is an
/// argument syntax error (exit 2) in clap's words. Ensembl VEP 116.3 splits the
/// value on commas before its `%VALID` test (Config.pm:946) and reads `b,` by its
/// first letter (OutputFactory.pm:1817), annotating `SIFT=tolerated(0.19)` on
/// `21 14216193 C G`; `docs/cli/compatibility.md` states the difference.
#[test]
fn a_comma_joined_enumerated_value_is_a_syntax_error() {
    for (flag, value) in [
        ("--sift", "b,"),
        ("--polyphen", "b,p"),
        ("--format", "vcf,"),
    ] {
        let output = vep()
            .args([flag, value, "-i", "/nonexistent/input.vcf", "--quiet"])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        let stderr = text(&output.stderr);
        assert_eq!(output.status.code(), Some(2), "{flag} {value}: {stderr}");
        assert!(
            stderr.contains(&format!("error: invalid value '{value}' for '{flag}")),
            "{flag} {value}: {stderr}"
        );
    }
}

/// A value Perl reads as false, `0` or the empty string, given to a
/// clap-enumerated `%VALID` flag (`--sift` and `--polyphen` take `s`, `p` or
/// `b` and `--nearest` takes `transcript`, `gene` or `symbol`, the
/// `value_parser` lists of `args.rs`) is an argument syntax error (exit 2) in
/// clap's words, `0` an invalid value and the empty string a missing one,
/// before any run and with no output file written. Ensembl VEP 116.3 reads such
/// a value as the flag unset:
/// `_is_flag_active` (Config.pm:902-905, `return 0 unless
/// defined($config->{$flag}); return ref($config->{$flag}) eq 'ARRAY' ? scalar
/// @{$config->{$flag}} : $config->{$flag};`) keeps it out of the `%VALID` test
/// (`:945`, `foreach my $flag(grep {$self->_is_flag_active($config, $_)} keys
/// %VALID) {`), and the column gates are truthiness tests too,
/// OutputFactory.pm:1814 (`if (my $opt = $self->{$lc_tool}) {`) for SIFT and
/// PolyPhen and AnnotationType/Transcript.pm:170 (`if(my $nearest_type =
/// $self->{nearest}) {`) for NEAREST, so `--sift 0`, `--sift ''`, `--polyphen 0`
/// and `--nearest 0` each exit 0 there and annotate without the column;
/// `--sift x` exits 2 in both. `docs/cli/compatibility.md` states the
/// difference.
#[test]
fn a_false_value_of_an_enumerated_flag_is_a_syntax_error() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    for (flag, value, first_line) in [
        (
            "--sift",
            "0",
            "error: invalid value '0' for '--sift <SIFT>'",
        ),
        (
            "--sift",
            "",
            "error: a value is required for '--sift <SIFT>' but none was supplied",
        ),
        (
            "--polyphen",
            "0",
            "error: invalid value '0' for '--polyphen <POLYPHEN>'",
        ),
        (
            "--nearest",
            "0",
            "error: invalid value '0' for '--nearest <TYPE>'",
        ),
    ] {
        let output = vep()
            .args([
                flag,
                value,
                "-i",
                "/nonexistent/input.vcf",
                "-o",
                out.to_str().unwrap(),
                "--quiet",
            ])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        let stderr = text(&output.stderr);
        assert_eq!(output.status.code(), Some(2), "{flag} {value:?}: {stderr}");
        assert_eq!(
            stderr.lines().next(),
            Some(first_line),
            "{flag} {value:?}: {stderr}"
        );
        assert!(
            !out.exists(),
            "{flag} {value:?}: an output file was created"
        );
    }
}

/// The `-o` value names the output file as given: `-o ' lead.txt'` writes a file
/// named ` lead.txt`, where Ensembl VEP's two-argument `open` (BaseRunner.pm:226)
/// strips the whitespace and writes `lead.txt`; `docs/cli/compatibility.md` states
/// the difference.
#[test]
fn output_file_name_is_taken_as_given() {
    let dir = tempfile::tempdir().unwrap();
    let output = run_control(&["-o", " lead.txt", "--quiet"], dir.path());
    assert!(output.status.success(), "{}", text(&output.stderr));
    assert!(dir.path().join(" lead.txt").is_file());
    assert!(!dir.path().join("lead.txt").exists());
    assert!(
        !data_lines(&std::fs::read_to_string(dir.path().join(" lead.txt")).unwrap()).is_empty()
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

/// The statistics report Ensembl VEP wrote for the statistics corpus's default run,
/// the fixture a `stats_file` expectation compares with.
fn statistics_reference() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-statistics/expected/default.summary.log");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// `text` with the row `key` of `[section]` carrying `value` instead.
fn with_summary_row(text: &str, section: &str, key: &str, value: &str) -> String {
    let mut current = String::new();
    let mut replaced = false;
    let lines: Vec<String> = text
        .lines()
        .map(|line| {
            if let Some(title) = line.strip_prefix('[').and_then(|l| l.strip_suffix(']')) {
                current = title.to_string();
            } else if current == section && line.split('\t').next() == Some(key) {
                replaced = true;
                return format!("{key}\t{value}");
            }
            line.to_string()
        })
        .collect();
    assert!(replaced, "no row {key:?} in [{section}]");
    lines.join("\n") + "\n"
}

/// A `stats_file` expectation compares the run's `--stats_text` report with the
/// reference's as `tools/vep-stats-compare` does, its rules for the two
/// documented differences excepted: the seven
/// volatile rows of `[VEP run statistics]` (the command line, the two times, the
/// run time, the input and output paths) compare by key alone, the paths of
/// `Annotation sources` are masked, a chart section compares as a set of
/// `(key, count)` pairs whatever its order, and everything else compares
/// verbatim: a `[General statistics]` count, a `[Data version]` row, a section
/// present on one side, the order of the tabular rows and a row printed twice are
/// each reported. The report's path follows Ensembl VEP's `--stats_file` rule
/// (BaseRunner.pm:254-271).
#[test]
fn stats_file_expectation_applies_the_comparators_masks() {
    let reference = statistics_reference();
    let run_stats = "VEP run statistics";
    let mut candidate = reference.clone();
    for (key, value) in [
        (
            "Annotation sources",
            "Cache: /home/user/vep-rs/tests/golden/116/GRCh38-statistics/json_cache",
        ),
        (
            "Command line options",
            "vep --assembly GRCh38 --json_cache [PATH]/json_cache --offline --stats_text",
        ),
        ("Start time", "2030-01-01 00:00:00"),
        ("End time", "2030-01-01 00:00:01"),
        ("Run time", "1 seconds"),
        (
            "Input file",
            "/home/user/vep-rs/tests/golden/116/GRCh38-statistics/variants.vcf",
        ),
        (
            "Output file",
            "/home/user/out/GRCh38-statistics-run-stats.out",
        ),
    ] {
        candidate = with_summary_row(&candidate, run_stats, key, value);
    }
    // The two seven-count classes in key order, as a deterministic engine sorts a tie.
    candidate = candidate.replace(
        "duplication\t1\nsequence_alteration\t7\nindel\t7\n",
        "duplication\t1\nindel\t7\nsequence_alteration\t7\n",
    );
    assert_ne!(candidate, reference);
    assert_eq!(compare_summary_sections(&reference, &candidate), "");
    let masked_cache = with_summary_row(
        &reference,
        run_stats,
        "Annotation sources",
        "Cache: [PATH]/116_GRCh38",
    );
    assert_eq!(compare_summary_sections(&reference, &masked_cache), "");
    let custom = with_summary_row(
        &reference,
        run_stats,
        "Annotation sources",
        "Cache: /work/cache/homo_sapiens/116_GRCh38; Custom: /work/corpus/custom/sites.vcf.gz (vcf)",
    );
    assert_eq!(
        compare_summary_sections(&reference, &custom),
        "  [VEP run statistics]: Annotation sources: expected \"Cache: [PATH]\", vep-rs \"Cache: [PATH]; Custom: [PATH] (vcf)\"\n"
    );
    let species = with_summary_row(&reference, run_stats, "Species", "mus_musculus");
    assert_eq!(
        compare_summary_sections(&reference, &species),
        "  [VEP run statistics]: Species: expected \"homo_sapiens\", vep-rs \"mus_musculus\"\n"
    );
    let genes = with_summary_row(&reference, "General statistics", "Overlapped genes", "97");
    assert_eq!(
        compare_summary_sections(&reference, &genes),
        "  [General statistics]: Overlapped genes: expected \"98\", vep-rs \"97\"\n"
    );
    let version = with_summary_row(&reference, "Data version", "dbSNP", "157");
    assert_eq!(
        compare_summary_sections(&reference, &version),
        "  [Data version]: dbSNP: expected \"156\", vep-rs \"157\"\n"
    );
    let count = with_summary_row(&reference, "Consequences (all)", "intron_variant", "964");
    assert_eq!(
        compare_summary_sections(&reference, &count),
        "  [Consequences (all)]: MISSING in vep-rs: intron_variant\t963\n  [Consequences (all)]: UNEXPECTED in vep-rs: intron_variant\t964\n"
    );
    let without_protein = reference
        .split("\n[Position in protein]\n")
        .next()
        .unwrap()
        .to_string()
        + "\n";
    assert_eq!(
        compare_summary_sections(&reference, &without_protein),
        "  MISSING in vep-rs: section [Position in protein]\n"
    );
    assert_eq!(
        compare_summary_sections(&without_protein, &reference),
        "  UNEXPECTED in vep-rs: section [Position in protein]\n"
    );
    let reordered = reference.replacen(
        "VEP version (API)\t 116.3 (116)\nAnnotation sources\tCache: /work/cache/homo_sapiens/116_GRCh38\n",
        "Annotation sources\tCache: /work/cache/homo_sapiens/116_GRCh38\nVEP version (API)\t 116.3 (116)\n",
        1,
    );
    let report = compare_summary_sections(&reference, &reordered);
    assert!(
        report.starts_with("  [VEP run statistics]: rows in another order: expected [\"VEP version (API)\", \"Annotation sources\", \"Species\""),
        "{report}"
    );
    let doubled = reference.replacen("SNV\t243\n", "SNV\t243\nSNV\t1\n", 1);
    assert_eq!(
        compare_summary_sections(&reference, &doubled),
        "  [Variant classes]: SNV printed twice in vep-rs\n  [Variant classes]: UNEXPECTED in vep-rs: SNV\t1\n"
    );
    let repeated = reference.replacen(
        "Overlapped genes\t98\n",
        "Overlapped genes\t98\nOverlapped genes\t1\n",
        1,
    );
    assert_eq!(
        compare_summary_sections(&reference, &repeated),
        "  [General statistics]: Overlapped genes printed twice in vep-rs\n"
    );
    let twice = format!("{reference}\n[Variant classes]\nSNV\t243\n");
    assert_eq!(
        compare_summary_sections(&reference, &twice),
        "  section [Variant classes] appears twice in vep-rs\n"
    );
    assert_eq!(
        compare_summary_sections(&reference, ""),
        reference
            .lines()
            .filter(|l| l.starts_with('['))
            .map(|l| format!("  MISSING in vep-rs: section {l}\n"))
            .collect::<String>()
    );
    for (root, file) in [
        ("out.txt_summary", "out.txt_summary.txt"),
        ("/run/out.txt.stats.html", "/run/out.txt.stats.txt"),
        ("stats.TXT", "stats.TXT"),
        ("report.foo", "report.foo.txt"),
        ("report", "report.txt"),
        ("/data/a.b/report", "/data/a.b/report.txt"),
    ] {
        assert_eq!(stats_text_file(root), PathBuf::from(file), "{root}");
    }
    // The run's report: `--stats_file` or `--sf` through the rule above, else beside the
    // output, which for a standard-output run is the `-o` spelling in the working directory.
    let out = Path::new("/data/run/x.out");
    let dir = Path::new("/data/run");
    let argv = |words: &[&str]| words.iter().map(|w| w.to_string()).collect::<Vec<_>>();
    assert_eq!(
        stats_report_path(&argv(&["-i", "in.vcf", "-o", "/data/run/x.out"]), out, dir),
        PathBuf::from("/data/run/x.out_summary.txt")
    );
    assert_eq!(
        stats_report_path(&argv(&["--stats_file", "/data/run/x.stats.html"]), out, dir),
        PathBuf::from("/data/run/x.stats.txt")
    );
    assert_eq!(
        stats_report_path(&argv(&["--sf", "/data/run/named"]), out, dir),
        PathBuf::from("/data/run/named.txt")
    );
    assert_eq!(
        stats_report_path(&argv(&["-o", "STDOUT", "--stats_text"]), out, dir),
        PathBuf::from("/data/run/STDOUT_summary.txt")
    );
    assert_eq!(
        stats_report_path(&argv(&["--output_file", "stdout"]), out, dir),
        PathBuf::from("/data/run/stdout_summary.txt")
    );
    // The HTML report takes the same root with the other extension.
    for (words, file) in [
        (
            &["-i", "in.vcf", "-o", "/data/run/x.out"][..],
            "/data/run/x.out_summary.html",
        ),
        (
            &["--stats_file", "/data/run/x.stats.txt"][..],
            "/data/run/x.stats.html",
        ),
        (
            &["--sf", "/data/run/named.HTML"][..],
            "/data/run/named.HTML",
        ),
        (
            &["--sf", "/data/run/report.foo"][..],
            "/data/run/report.foo.html",
        ),
        (&["-o", "StdOut"][..], "/data/run/StdOut_summary.html"),
    ] {
        assert_eq!(
            stats_html_report_path(&argv(words), out, dir),
            PathBuf::from(file),
            "{words:?}"
        );
    }
    let lacking = check_html_report("<html><title>VEP summary</title><script></script>");
    assert!(
        lacking.contains("  no element with id masthead\n"),
        "{lacking}"
    );
    assert!(
        lacking.contains("  the page holds a <script element\n"),
        "{lacking}"
    );
    assert!(!lacking.contains("remote"), "{lacking}");
}

/// The GRCh37-hgvs corpus, the one committed with a reference FASTA.
fn hgvs_corpus_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37-hgvs")
}

/// The GRCh37-hgvs corpus's bgzipped FASTA, indexed with its `.fai` and `.gzi`.
fn hgvs_corpus_fasta() -> PathBuf {
    hgvs_corpus_dir().join("reference.fa.gz")
}

/// A directory holding the corpus FASTA and its indexes under the name Ensembl
/// gives a chromosome file, which carries the species and the assembly.
fn ensembl_named_fasta_dir(parent: &Path) -> PathBuf {
    let dir = parent.join("fasta_dir");
    std::fs::create_dir_all(&dir).unwrap();
    let source = hgvs_corpus_fasta();
    for suffix in ["", ".fai", ".gzi"] {
        std::fs::copy(
            format!("{}{suffix}", source.display()),
            dir.join(format!(
                "Homo_sapiens.GRCh37.dna.chromosome.21.fa.gz{suffix}"
            )),
        )
        .unwrap();
    }
    dir
}

/// The GRCh37-hgvs corpus's own arguments without its FASTA and its `--hgvs`,
/// output path left to the caller.
fn hgvs_corpus_args() -> Vec<String> {
    let dir = hgvs_corpus_dir();
    vec![
        "-i".to_string(),
        dir.join("variants.vcf").display().to_string(),
        "--offline".to_string(),
        "--json_cache".to_string(),
        dir.join("json_cache").display().to_string(),
        "--assembly".to_string(),
        "GRCh37".to_string(),
        "--force_overwrite".to_string(),
        "--no_stats".to_string(),
    ]
}

/// The line a run asking for HGVS notations with no FASTA ends with.
const HGVS_REFUSAL: &str = "ERROR: [fasta] Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode without a FASTA file (see --fasta)";

/// `--hgvs` or `--hgvsg` with no FASTA in reach (none named, none in the cache
/// directory) is refused under the fasta component with Ensembl VEP's text
/// (Runner.pm:757-759) before the output is created, and before an existing
/// output file is refused (Runner.pm:200-202 opens the FASTA and runs its checks
/// in init(), ahead of the output handle), as a `--fasta` path that is not there
/// is (Utils/FastaSequence.pm:187); an input file that is not there is refused
/// ahead of the HGVS check, as init() opens the input between the two
/// (Runner.pm:134-138). With a FASTA the same flags annotate.
#[test]
fn hgvs_without_fasta_is_fatal() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    for flags in [
        &["--hgvs"][..],
        &["--hgvsg"][..],
        &["--hgvs", "--hgvsg"][..],
    ] {
        let mut extra = flags.to_vec();
        extra.push("--quiet");
        let refused = run_with(&corpus_args(), &extra, &out);
        assert_eq!(refused.status.code(), Some(1), "{flags:?}");
        assert_eq!(text(&refused.stderr).trim_end(), HGVS_REFUSAL, "{flags:?}");
        assert!(refused.stdout.is_empty(), "{flags:?}");
        assert!(!out.exists(), "{flags:?}: an output file was created");
    }

    std::fs::write(&out, "kept\n").unwrap();
    let args = without_force_overwrite(corpus_args());
    let refused = run_with(&args, &["--hgvs", "--quiet"], &out);
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(text(&refused.stderr).trim_end(), HGVS_REFUSAL);
    assert_eq!(std::fs::read_to_string(&out).unwrap(), "kept\n");
    let missing = dir.path().join("absent.fa");
    let refused = run_with(
        &args,
        &["--quiet", "--fasta", missing.to_str().unwrap()],
        &out,
    );
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        format!(
            "ERROR: [fasta] Specified FASTA file/directory {} not found",
            missing.display()
        )
    );
    assert_eq!(std::fs::read_to_string(&out).unwrap(), "kept\n");
    let absent_input = dir.path().join("absent.vcf");
    let refused = run_with(
        &corpus_args_with_input(&absent_input),
        &["--hgvs", "--quiet"],
        &out,
    );
    assert_eq!(refused.status.code(), Some(1));
    assert!(
        text(&refused.stderr).starts_with(&format!(
            "ERROR: [input] File \"{}\" does not exist",
            absent_input.display()
        )),
        "{}",
        text(&refused.stderr)
    );

    let fasta = hgvs_corpus_fasta().display().to_string();
    let annotated = run_with(
        &hgvs_corpus_args(),
        &["--hgvs", "--hgvsg", "--quiet", "--fasta", &fasta],
        &out,
    );
    assert_eq!(
        annotated.status.code(),
        Some(0),
        "{}",
        text(&annotated.stderr)
    );
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(output.contains("## HGVSc : HGVS coding sequence name"));
    assert!(
        data_lines(&output).iter().any(|l| l.contains("HGVSc=")),
        "no HGVSc value in the output"
    );
}

/// `--everything` without a FASTA prints `INFO: Disabling --hgvs; using --offline
/// and no FASTA file found` on standard error and drops the HGVS notations it
/// switched on, so the output carries no HGVSc or HGVSp (Runner.pm:745-754);
/// `--quiet` withholds the notice as status_msg does (BaseVEP.pm:782); `--hgvsg`
/// is not in the set the notice clears, so `--everything --hgvsg` is still
/// refused; with a FASTA nothing is disabled.
#[test]
fn everything_without_fasta_prints_info_and_disables_hgvs() {
    const NOTICE: &str = "INFO: Disabling --hgvs; using --offline and no FASTA file found";
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let disabled = run_with(&corpus_args(), &["--everything"], &out);
    assert_eq!(
        disabled.status.code(),
        Some(0),
        "{}",
        text(&disabled.stderr)
    );
    let stderr = text(&disabled.stderr);
    assert!(
        stderr.lines().any(|l| l == NOTICE),
        "the notice is a line of its own:\n{stderr}"
    );
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(!output.contains("HGVSc"), "HGVSc survived the notice");
    assert!(!output.contains("HGVSp"), "HGVSp survived the notice");
    assert!(output.contains("## SYMBOL"), "the rest of the set is on");

    let quiet = run_with(&corpus_args(), &["--everything", "--quiet"], &out);
    assert_eq!(quiet.status.code(), Some(0), "{}", text(&quiet.stderr));
    assert!(!text(&quiet.stderr).contains("INFO: Disabling"));

    let refused = run_with(
        &corpus_args(),
        &["--everything", "--hgvsg", "--quiet"],
        &out,
    );
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(text(&refused.stderr).trim_end(), HGVS_REFUSAL);

    let fasta = hgvs_corpus_fasta().display().to_string();
    let with_fasta = run_with(
        &hgvs_corpus_args(),
        &["--everything", "--fasta", &fasta],
        &out,
    );
    assert_eq!(
        with_fasta.status.code(),
        Some(0),
        "{}",
        text(&with_fasta.stderr)
    );
    assert!(!text(&with_fasta.stderr).contains("INFO: Disabling"));
    let output = std::fs::read_to_string(&out).unwrap();
    assert!(output.contains("## HGVSc : HGVS coding sequence name"));
}

/// Ensembl VEP 116.3 prints `INFO: Disabling --hgvs; using --offline and no
/// FASTA file found` as it reaches it (Runner.pm:752, BaseVEP.pm:785) and then
/// throws for a flag `--everything` does not clear that needs a FASTA
/// (:759-762), so `--everything` with `--hgvsg`, `--check_ref`, `--lookup_ref`
/// or `--use_transcript_ref` and no FASTA prints the notice before the `ERROR:`
/// line, the error last; `--quiet` withholds the notice and leaves the error
/// alone (BaseVEP.pm:782).
#[test]
fn the_everything_notice_prints_before_each_fasta_refusal() {
    const NOTICE: &str = "INFO: Disabling --hgvs; using --offline and no FASTA file found";
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    for (flag, error) in [
        ("--hgvsg", HGVS_REFUSAL),
        (
            "--check_ref",
            "ERROR: [fasta] Cannot check reference sequences (--check_ref) without a FASTA file (see --fasta)",
        ),
        (
            "--lookup_ref",
            "ERROR: [fasta] Cannot look up reference sequences (--lookup_ref) without a FASTA file (see --fasta)",
        ),
        (
            "--use_transcript_ref",
            "ERROR: [fasta] Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref",
        ),
    ] {
        let refused = run_with(&corpus_args(), &["--everything", flag], &out);
        assert_eq!(refused.status.code(), Some(1), "{flag}");
        let stderr = text(&refused.stderr);
        let lines: Vec<&str> = stderr.lines().collect();
        let info_lines: Vec<&str> = lines
            .iter()
            .copied()
            .filter(|l| l.starts_with("INFO: "))
            .collect();
        assert_eq!(info_lines, [NOTICE], "{flag}: the INFO lines of:\n{stderr}");
        let notice_at = lines.iter().position(|l| *l == NOTICE).unwrap();
        let error_at = lines
            .iter()
            .position(|l| *l == error)
            .unwrap_or_else(|| panic!("{flag}: no error line in:\n{stderr}"));
        assert!(notice_at < error_at, "{flag}: the error precedes the notice:\n{stderr}");
        assert_eq!(error_at, lines.len() - 1, "{flag}: the error is not last:\n{stderr}");

        let quiet = run_with(&corpus_args(), &["--everything", flag, "--quiet"], &out);
        assert_eq!(quiet.status.code(), Some(1), "{flag}");
        assert_eq!(text(&quiet.stderr).trim_end(), error, "{flag}");
    }
}

/// Ensembl VEP 116.3 prints `INFO: Colocated variants will be provided when
/// using --var_synonyms` as `init` begins (Runner.pm:117-120), before the
/// annotation sources are built (:128), where a missing cache directory is
/// refused (CacheDir.pm:109) and the transcript filter is compiled
/// (AnnotationSource.pm:87-90), so `--var_synonyms` with a missing cache
/// directory or a filter Ensembl VEP refuses prints the notice before the
/// `ERROR:` line, the error last; `--quiet` withholds the notice
/// (BaseVEP.pm:782), and `--show_cache_info` runs no `init` (vep:37-43), so
/// it prints none.
#[test]
fn the_var_synonyms_notice_prints_before_a_cache_or_filter_refusal() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let missing = dir.path().join("missing_cache");
    let mut missing_cache = corpus_args();
    let cache_flag = missing_cache
        .iter()
        .position(|a| a == "--json_cache")
        .unwrap();
    missing_cache[cache_flag + 1] = missing.display().to_string();
    let cache_error = format!(
        "ERROR: [cache] Directory {} does not exist",
        missing.display()
    );
    let with_cache = corpus_args();
    let filter_error = "ERROR: [filter] No such operator \"equals\"";
    for (args, extra, error) in [
        (
            &missing_cache,
            &["--var_synonyms"][..],
            cache_error.as_str(),
        ),
        (
            &with_cache,
            &["--var_synonyms", "--transcript_filter", "biotype equals x"][..],
            filter_error,
        ),
    ] {
        let refused = run_with(args, extra, &out);
        assert_eq!(refused.status.code(), Some(1), "{extra:?}");
        let stderr = text(&refused.stderr);
        let lines: Vec<&str> = stderr.lines().collect();
        let info_lines: Vec<&str> = lines
            .iter()
            .copied()
            .filter(|l| l.starts_with("INFO: "))
            .collect();
        assert_eq!(
            info_lines,
            [VAR_SYNONYMS_NOTICE],
            "{extra:?}: the INFO lines of:\n{stderr}"
        );
        let notice_at = lines
            .iter()
            .position(|l| *l == VAR_SYNONYMS_NOTICE)
            .unwrap();
        let error_at = lines
            .iter()
            .position(|l| *l == error)
            .unwrap_or_else(|| panic!("{extra:?}: no error line in:\n{stderr}"));
        assert!(
            notice_at < error_at,
            "{extra:?}: the error precedes the notice:\n{stderr}"
        );
        assert_eq!(
            error_at,
            lines.len() - 1,
            "{extra:?}: the error is not last:\n{stderr}"
        );

        let mut quiet = extra.to_vec();
        quiet.push("--quiet");
        let quiet = run_with(args, &quiet, &out);
        assert_eq!(quiet.status.code(), Some(1), "{extra:?}");
        assert_eq!(text(&quiet.stderr).trim_end(), error, "{extra:?}");
    }

    let info = vep()
        .args(&missing_cache)
        .args(["--var_synonyms", "--show_cache_info"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(info.status.code(), Some(1));
    let stderr = text(&info.stderr);
    assert!(!stderr.contains(VAR_SYNONYMS_NOTICE), "{stderr}");
    assert_eq!(
        stderr.lines().last(),
        Some(cache_error.as_str()),
        "{stderr}"
    );
}

/// A `--fasta` given without its indexes is indexed beside itself, as Ensembl
/// VEP's reader indexes it (`Variation/Utils/FastaSequence.pm:323-358`, through
/// `Bio::DB::HTS::Faidx`): the bgzipped corpus FASTA gets the `.fai` and `.gzi`
/// `samtools faidx` wrote for it, its plain form a `.fai`, and either run
/// writes the HGVS output of the run on the indexed corpus FASTA.
#[test]
fn a_fasta_without_its_indexes_is_indexed_beside_itself() {
    let dir = tempfile::tempdir().unwrap();
    let source = hgvs_corpus_fasta();
    let indexed_out = dir.path().join("indexed.txt");
    let indexed = run_with(
        &hgvs_corpus_args(),
        &["--hgvs", "--quiet", "--fasta", source.to_str().unwrap()],
        &indexed_out,
    );
    assert!(indexed.status.success(), "{}", text(&indexed.stderr));
    // The header names the run's time and its FASTA; the rows follow it.
    let rows = |path: &Path| -> String {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with("##"))
            .map(|line| format!("{line}\n"))
            .collect()
    };
    let expected = rows(&indexed_out);
    assert!(expected.lines().count() > 1, "{expected}");

    let bgzf = dir.path().join("bgzf.fa.gz");
    std::fs::copy(&source, &bgzf).unwrap();
    let mut plain_bases = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&source).unwrap())
        .read_to_end(&mut plain_bases)
        .unwrap();
    let plain = dir.path().join("plain.fa");
    std::fs::write(&plain, plain_bases).unwrap();

    for (fasta, indexes) in [(&bgzf, &["fai", "gzi"][..]), (&plain, &["fai"][..])] {
        let out = dir.path().join("out.txt");
        let run = run_with(
            &hgvs_corpus_args(),
            &["--hgvs", "--quiet", "--fasta", fasta.to_str().unwrap()],
            &out,
        );
        assert!(
            run.status.success(),
            "{}: {}",
            fasta.display(),
            text(&run.stderr)
        );
        assert!(run.stderr.is_empty(), "{}", text(&run.stderr));
        assert_eq!(rows(&out), expected, "{}", fasta.display());
        for extension in indexes {
            assert!(
                Path::new(&format!("{}.{extension}", fasta.display())).is_file(),
                "{}.{extension}",
                fasta.display()
            );
        }
        assert!(!Path::new(&format!("{}.vep.lock", fasta.display())).exists());
    }
    for extension in ["fai", "gzi"] {
        assert_eq!(
            std::fs::read(format!("{}.{extension}", bgzf.display())).unwrap(),
            std::fs::read(format!("{}.{extension}", source.display())).unwrap(),
            "{extension}"
        );
    }
}

/// The file `--fasta_dir` picks is indexed beside itself as the file `--fasta`
/// names is: the bgzipped corpus FASTA under a name that carries the species
/// and the assembly, given without its indexes, gets the `.fai` and `.gzi`
/// `samtools faidx` wrote for it, no lock file is left, and the run writes the
/// HGVS output of the run on the indexed corpus FASTA.
#[test]
fn the_file_fasta_dir_picks_is_indexed_beside_itself() {
    let dir = tempfile::tempdir().unwrap();
    let source = hgvs_corpus_fasta();
    let indexed_out = dir.path().join("indexed.txt");
    let indexed = run_with(
        &hgvs_corpus_args(),
        &["--hgvs", "--quiet", "--fasta", source.to_str().unwrap()],
        &indexed_out,
    );
    assert!(indexed.status.success(), "{}", text(&indexed.stderr));
    let expected = std::fs::read_to_string(&indexed_out).unwrap();
    assert!(data_lines(&expected).len() > 1, "{expected}");

    let fasta_dir = dir.path().join("fasta_dir");
    std::fs::create_dir(&fasta_dir).unwrap();
    let picked = fasta_dir.join("Homo_sapiens.GRCh37.dna.chromosome.21.fa.gz");
    std::fs::copy(&source, &picked).unwrap();
    let out = dir.path().join("out.txt");
    let run = run_with(
        &hgvs_corpus_args(),
        &[
            "--hgvs",
            "--quiet",
            "--fasta_dir",
            fasta_dir.to_str().unwrap(),
        ],
        &out,
    );
    assert!(run.status.success(), "{}", text(&run.stderr));
    assert!(run.stderr.is_empty(), "{}", text(&run.stderr));
    let actual = std::fs::read_to_string(&out).unwrap();
    assert_eq!(data_lines(&actual), data_lines(&expected));
    for extension in ["fai", "gzi"] {
        assert_eq!(
            std::fs::read(format!("{}.{extension}", picked.display())).unwrap(),
            std::fs::read(format!("{}.{extension}", source.display())).unwrap(),
            "{extension}"
        );
    }
    assert!(!Path::new(&format!("{}.vep.lock", picked.display())).exists());
}

/// A `--fasta` without its index in a directory the run cannot write ends the
/// run with one line under the fasta component: Ensembl VEP's sentence for
/// the lock file it cannot write (`Variation/Utils/FastaSequence.pm:334`),
/// then the cause and the two ways on.
#[cfg(unix)]
#[test]
fn a_fasta_whose_lock_file_cannot_be_written_ends_the_run_with_one_line() {
    use std::os::unix::fs::PermissionsExt;
    let dir = tempfile::tempdir().unwrap();
    let closed = dir.path().join("closed");
    std::fs::create_dir(&closed).unwrap();
    let fasta = closed.join("made_up.fa");
    std::fs::write(&fasta, ">21\nACGTACGT\n").unwrap();
    std::fs::set_permissions(&closed, std::fs::Permissions::from_mode(0o555)).unwrap();
    // Root writes into the directory anyway, and the state cannot be made.
    let unwritable = std::fs::write(closed.join("probe"), "").is_err();
    let run = run_with(
        &hgvs_corpus_args(),
        &["--hgvs", "--quiet", "--fasta", fasta.to_str().unwrap()],
        &dir.path().join("out.txt"),
    );
    std::fs::set_permissions(&closed, std::fs::Permissions::from_mode(0o755)).unwrap();
    if !unwritable {
        return;
    }
    assert_eq!(run.status.code(), Some(1), "{}", text(&run.stderr));
    assert_eq!(
        text(&run.stderr),
        format!(
            "ERROR: [fasta] Could not write to FASTA lock file {}.vep.lock: Permission denied (os error 13); index the FASTA with samtools faidx, or link it into a directory the run can write\n",
            fasta.display()
        )
    );
}

/// A `.gzi` that misstates one block's length is refused at open; a damaged
/// block, whose ISIZE field and the index agree while its payload inflates to
/// another length, opens but makes every read of that block "no reference", so
/// the run reports it once and exits 1 under the fasta component once its reads
/// are done, and the output written meanwhile is not taken for the file's
/// annotation. Both fixtures are the committed FASTA with the block whose
/// sequence most of the corpus's variants read (entry 544 ends it) altered.
#[test]
fn misdescribed_block_fails_the_run() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let source = hgvs_corpus_fasta();
    let staged = dir.path().join("reference.fa.gz");
    std::fs::copy(&source, &staged).unwrap();
    std::fs::copy(
        format!("{}.fai", source.display()),
        format!("{}.fai", staged.display()),
    )
    .unwrap();
    let source_gzi = format!("{}.gzi", source.display());
    let staged_gzi = format!("{}.gzi", staged.display());
    let index: Vec<(u64, u64)> = gzi::fs::read(&source_gzi).unwrap().as_ref().to_vec();
    assert!(
        index.len() > 600,
        "the corpus index has {} entries",
        index.len()
    );
    let flags = ["--hgvs", "--quiet", "--fasta", staged.to_str().unwrap()];

    let mut misstated = index.clone();
    misstated[544].1 += 1;
    gzi::fs::write(&staged_gzi, &gzi::Index::from(misstated)).unwrap();
    let refused = run_with(&hgvs_corpus_args(), &flags, &out);
    assert_eq!(refused.status.code(), Some(1));
    let stderr = text(&refused.stderr);
    assert!(
        stderr.starts_with(&format!(
            "ERROR: [fasta] BGZF index {staged_gzi} does not describe {} (the block at byte {} ends at byte",
            staged.display(),
            index[543].0
        )),
        "{stderr}"
    );
    assert!(!out.exists(), "refused before the output opened");

    let mut damaged = std::fs::read(&source).unwrap();
    let isize_at = index[544].0 as usize - 4;
    let isize = u32::from_le_bytes(damaged[isize_at..isize_at + 4].try_into().unwrap());
    damaged[isize_at..isize_at + 4].copy_from_slice(&(isize + 1).to_le_bytes());
    std::fs::write(&staged, &damaged).unwrap();
    let mut agreeing = index.clone();
    for entry in &mut agreeing[544..] {
        entry.1 += 1;
    }
    gzi::fs::write(&staged_gzi, &gzi::Index::from(agreeing)).unwrap();
    let failed = run_with(&hgvs_corpus_args(), &flags, &out);
    assert_eq!(failed.status.code(), Some(1));
    let stderr = text(&failed.stderr);
    assert!(
        stderr.trim_end().ends_with(&format!(
            "ERROR: [fasta] BGZF index {staged_gzi} does not describe {} (a block did not read as the index says): index the file again with samtools faidx",
            staged.display()
        )),
        "{stderr}"
    );
    assert!(
        out.exists(),
        "the output was written before the reads were judged"
    );

    std::fs::copy(&source, &staged).unwrap();
    std::fs::copy(&source_gzi, &staged_gzi).unwrap();
    let annotated = run_with(&hgvs_corpus_args(), &flags, &out);
    assert_eq!(
        annotated.status.code(),
        Some(0),
        "{}",
        text(&annotated.stderr)
    );
}

/// `--fasta_dir` supplies the FASTA whose name carries the species and the
/// assembly (BaseVEP.pm:495-513), so a run with it annotates as the same run
/// with `--fasta` does; a directory holding no such file (the corpus directory,
/// whose `reference.fa.gz` names neither; one holding another species' file)
/// warns and leaves the run without a FASTA (:516), and one that cannot be
/// opened is passed over without a word (:496).
#[test]
fn fasta_dir_selects_the_species_and_assembly_file() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let corpus = hgvs_corpus_dir().display().to_string();
    let fasta = hgvs_corpus_fasta().display().to_string();
    let by_file = run_with(
        &hgvs_corpus_args(),
        &["--hgvs", "--quiet", "--fasta", &fasta],
        &out,
    );
    assert_eq!(by_file.status.code(), Some(0), "{}", text(&by_file.stderr));
    let expected = std::fs::read_to_string(&out).unwrap();

    let named = ensembl_named_fasta_dir(dir.path()).display().to_string();
    let by_dir = run_with(
        &hgvs_corpus_args(),
        &["--hgvs", "--quiet", "--fasta_dir", &named],
        &out,
    );
    assert_eq!(by_dir.status.code(), Some(0), "{}", text(&by_dir.stderr));
    let actual = std::fs::read_to_string(&out).unwrap();
    assert_eq!(data_lines(&actual), data_lines(&expected));
    assert!(data_lines(&actual).iter().any(|l| l.contains("HGVSc=")));

    let warned = run_with(
        &hgvs_corpus_args(),
        &["--hgvs", "--fasta_dir", &corpus],
        &out,
    );
    assert_eq!(warned.status.code(), Some(1));
    let stderr = text(&warned.stderr);
    assert!(
        stderr.contains(&format!("No compatible FASTA files found in {corpus}")),
        "{stderr}"
    );
    assert!(
        stderr.contains("ERROR: [fasta] Cannot generate HGVS coordinates"),
        "{stderr}"
    );

    let other_species = run_with(
        &hgvs_corpus_args(),
        &["--hgvs", "--fasta_dir", &named, "--species", "mus_musculus"],
        &out,
    );
    assert_eq!(other_species.status.code(), Some(1));
    assert!(
        text(&other_species.stderr)
            .contains(&format!("No compatible FASTA files found in {named}")),
        "{}",
        text(&other_species.stderr)
    );

    let unreadable = dir.path().join("absent").display().to_string();
    let passed_over = run_with(
        &hgvs_corpus_args(),
        &["--hgvs", "--fasta_dir", &unreadable, "--fasta", &fasta],
        &out,
    );
    assert_eq!(
        passed_over.status.code(),
        Some(0),
        "{}",
        text(&passed_over.stderr)
    );
    assert!(!text(&passed_over.stderr).contains("No compatible FASTA files"));
    assert_eq!(
        data_lines(&std::fs::read_to_string(&out).unwrap()),
        data_lines(&expected)
    );
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

/// The run-control corpus as `corpora` reads it, manifest included.
fn run_control_corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-run-control")
        .expect("the GRCh38-run-control corpus is committed")
}

/// `source` with `manifest` in place of its own.
fn with_manifest(source: &Corpus, manifest: serde_json::Value) -> Corpus {
    Corpus {
        release: source.release.clone(),
        name: source.name.clone(),
        assembly: source.assembly.clone(),
        species: source.species.clone(),
        dir: source.dir.clone(),
        manifest,
    }
}

/// The message of a caught panic: a formatted `panic!` carries a `String`, a
/// literal one a `&str`; any other payload reads as empty.
fn panic_message(panic: &(dyn std::any::Any + Send)) -> String {
    panic
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| panic.downcast_ref::<&str>().map(|s| s.to_string()))
        .unwrap_or_default()
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
    let source = run_control_corpus();
    let label = format!("{}/{}", source.release, source.name);
    let with_runs = |runs: serde_json::Value| {
        let mut manifest = source.manifest.clone();
        manifest["runs"] = runs;
        with_manifest(&source, manifest)
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
                    "skipped_file",
                    "stats_file",
                    "stats_html"
                ]
            ),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "stats_html": "yes"})),
            format!("{label} run x: expect.stats_html is a boolean"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "stats_html": false})),
            format!("{label} run x: expects nothing beyond the exit status"),
        ),
        (
            serde_json::json!([{"name": "x", "argv": ["{golden}", "-o", "{out}"], "expect": {"exit": 0, "stats_html": true}}]),
            format!("{label} run x: stats_html needs no --no_stats or {{golden}} in argv, and --stats_html beside --stats_text"),
        ),
        (
            serde_json::json!([{"name": "x", "argv": ["-i", "{input}", "--stats_text", "-o", "{out}"], "expect": {"exit": 0, "stats_html": true}}]),
            format!("{label} run x: stats_html needs no --no_stats or {{golden}} in argv, and --stats_html beside --stats_text"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "skipped_file": "skipped.txt"})),
            format!("{label} run x: skipped_file needs --skipped_variants_file in argv"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "stats_file": "default.summary.log"})),
            format!("{label} run x: stats_file needs --stats_text and no --no_stats or {{golden}} in argv"),
        ),
        (
            serde_json::json!([{"name": "x", "argv": ["{golden}", "--stats_text", "-o", "{out}"], "expect": {"exit": 0, "stats_file": "default.summary.log"}}]),
            format!("{label} run x: stats_file needs --stats_text and no --no_stats or {{golden}} in argv"),
        ),
        (
            expecting(serde_json::json!({"exit": 0, "stats_file": "../s.log"})),
            format!("{label} run x: stats_file \"../s.log\" is not a path inside the corpus"),
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
            Err(panic) => panic_message(panic.as_ref()),
        };
        if message != reason {
            wrong.push(format!("{runs}: expected {reason:?}, got {message:?}"));
        }
    }
    assert!(wrong.is_empty(), "\n{}", wrong.join("\n"));
    let corpus = with_runs(serde_json::json!([
        {"name": "x1_-", "argv": ["{golden}", "-o", "{out}"], "expect": {"exit": 0, "output": "default.txt"}},
        {"name": "y", "argv": [], "env": {"A": "1"}, "stdin": "variants.vcf", "expect": {"exit": 0, "format": "tab", "stderr_contains": "x"}},
        {"name": "s", "argv": ["-i", "{input}", "--stats_text", "-o", "{out}"], "expect": {"exit": 0, "stats_file": "default.summary.log"}},
        {"name": "h", "argv": ["-i", "{input}", "-o", "{out}"], "expect": {"exit": 0, "stats_html": true}},
        {"name": "hs", "argv": ["-i", "{input}", "--stats_text", "--stats_html", "-o", "{out}"], "expect": {"exit": 0, "stats_html": true}}
    ]));
    let accepted = manifest_runs(&corpus);
    assert_eq!(accepted.len(), 5);
    assert_eq!(accepted[0].expect.format, "default");
    assert_eq!(accepted[1].expect.format, "tab");
    assert_eq!(accepted[1].stdin.as_deref(), Some("variants.vcf"));
    assert_eq!(
        accepted[2].expect.stats_file.as_deref(),
        Some("default.summary.log")
    );
    assert!(!accepted[2].expect.stats_html);
    assert!(accepted[3].expect.stats_html && accepted[4].expect.stats_html);
}

/// A manifest's `cache_layout` selects the cache directory and the arguments
/// naming it: absent, `json_cache/` through `--json_cache`; `ensembl`, the
/// committed leaf `cache_root/homo_sapiens/<release>_<assembly>` through
/// `--cache --dir_cache <root> --cache_version <release>`, the species
/// directory suffixed `_refseq` or `_merged` when the manifest's flags carry
/// `--refseq` or `--merged`. Any other value, `null` included, is refused with
/// the corpus named rather than read as either.
#[test]
fn cache_layout_selects_the_cache_directory_and_refuses_another_value() {
    let source = run_control_corpus();
    let with_layout = |layout: Option<serde_json::Value>| {
        let mut manifest = source.manifest.clone();
        match layout {
            Some(value) => manifest["cache_layout"] = value,
            None => {
                manifest.as_object_mut().unwrap().remove("cache_layout");
            }
        }
        with_manifest(&source, manifest)
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
    // The species directory follows the manifest's flags as the binary follows
    // the run's (`CacheDir.pm:312`): `--refseq` reads `<species>_refseq`,
    // `--merged` `<species>_merged`, and the root passed stays the same.
    for (flag, suffix) in [("--refseq", "_refseq"), ("--merged", "_merged")] {
        let mut manifest = ensembl.manifest.clone();
        manifest["flags"] = serde_json::json!([flag, "--hgvs"]);
        let set = with_manifest(&source, manifest);
        assert_eq!(set.species_dir(), format!("homo_sapiens{suffix}"));
        assert_eq!(
            set.cache_dir(),
            root.join(format!("homo_sapiens{suffix}"))
                .join("116_GRCh38")
        );
        assert_eq!(set.cache_args(), ensembl.cache_args());
    }
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
            .map(|panic| panic_message(panic.as_ref()))
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

/// A `--fasta` value that is false in Perl, `0` or the empty value, is a path
/// here (`fasta: args.fasta` in `config.rs` takes it as given, and `open_fasta`
/// in `runner.rs` opens it), refused before any output is written with
/// vep-fasta's `Specified FASTA file/directory <value> not found` (exit 1),
/// where Ensembl VEP 116.3 reads the flag as unset: CacheDir.pm:235
/// `unless($self->param('fasta') || $self->param('no_fasta')) {` auto-detects the
/// cache directory's FASTA and BaseVEP.pm:521 `if(my $fasta_file =
/// $self->param('fasta')) {` opens none otherwise, the reading
/// `docs/cli/compatibility.md` and `docs/cache-setup.md` state.
#[test]
fn fasta_false_value_is_a_path() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    for value in ["0", ""] {
        let output = run_control(
            &["-o", out.to_str().unwrap(), "--fasta", value, "--quiet"],
            dir.path(),
        );
        assert_eq!(output.status.code(), Some(1), "{value:?}");
        assert_eq!(
            text(&output.stderr),
            format!("ERROR: [fasta] Specified FASTA file/directory {value} not found\n"),
            "{value:?}"
        );
        assert!(!out.exists(), "{value:?}: an output was written");
    }
}

/// `--fasta_dir` filters its candidates on the assembly the cache resolves to
/// (`fasta_source` in `runner.rs` passes the resolved cache's assembly, here the
/// corpus `info.json`'s `GRCh37`), so under `--json_cache` or its alias
/// `--full_cache_dir` with no `--assembly`, a directory holding a GRCh38-named
/// `toplevel` file beside the GRCh37 chromosome file opens the GRCh37 file and
/// annotates as `--fasta` does. Ensembl VEP 116.3 filters on `param('assembly')`
/// (BaseVEP.pm:497), which `--offline --full_cache_dir` without `--assembly` or
/// `--grch37` leaves unset (CacheDir.pm:356 fills it inside the layout scan of
/// :303 alone, and BaseRunner.pm:139 from a database `--offline` never opens),
/// so the `/$assembly/i` test of BaseVEP.pm:508 is an empty pattern, which Perl
/// reads as the last successful match, the species test of :507, and every
/// species file passes; :502 then ranks the `toplevel` name first and :513 makes
/// it the run's FASTA, which :523 opens. On such a directory Ensembl VEP 116.3
/// warns `Use of uninitialized value $assembly in regexp compilation` (reported
/// at :506, the grep's first expression) and annotates against the GRCh38-named
/// file, where with `--assembly GRCh37` it opens the GRCh37 file, as this run
/// does. `docs/cache-setup.md` states the difference.
#[test]
fn fasta_dir_filters_on_the_cache_assembly_without_assembly_flag() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let corpus = hgvs_corpus_dir();
    let input = corpus.join("variants.vcf").display().to_string();
    let cache = corpus.join("json_cache").display().to_string();
    let args = |cache_flag: &str| -> Vec<String> {
        [
            "-i",
            &input,
            "--offline",
            cache_flag,
            &cache,
            "--force_overwrite",
            "--no_stats",
            "--hgvs",
        ]
        .into_iter()
        .map(str::to_string)
        .collect()
    };
    let fasta = hgvs_corpus_fasta().display().to_string();
    let by_file = run_with(&args("--json_cache"), &["--quiet", "--fasta", &fasta], &out);
    assert_eq!(by_file.status.code(), Some(0), "{}", text(&by_file.stderr));
    let expected = std::fs::read_to_string(&out).unwrap();

    let mixed = ensembl_named_fasta_dir(dir.path());
    std::fs::write(
        mixed.join("homo_sapiens.GRCh38.dna.toplevel.fa"),
        ">21\nACGT\n",
    )
    .unwrap();
    let grch37 = mixed.join("Homo_sapiens.GRCh37.dna.chromosome.21.fa.gz");
    for cache_flag in ["--json_cache", "--full_cache_dir"] {
        let by_dir = run_with(
            &args(cache_flag),
            &["--fasta_dir", mixed.to_str().unwrap()],
            &out,
        );
        let stderr = text(&by_dir.stderr);
        assert_eq!(by_dir.status.code(), Some(0), "{cache_flag}: {stderr}");
        assert!(
            stderr.contains(&format!("Using reference FASTA: {}", grch37.display())),
            "{cache_flag}: {stderr}"
        );
        let actual = std::fs::read_to_string(&out).unwrap();
        assert_eq!(data_lines(&actual), data_lines(&expected), "{cache_flag}");
        assert!(
            data_lines(&actual).iter().any(|l| l.contains("HGVSc=")),
            "{cache_flag}"
        );
    }
}

/// The statistics corpus: a few hundred chromosome 21 records with Ensembl VEP's
/// `_summary.txt` beside each run's expected output.
fn statistics_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-statistics")
}

/// A corpus's arguments without the `--no_stats` every corpus run passes, for a
/// run that writes the statistics report.
fn without_no_stats(mut args: Vec<String>) -> Vec<String> {
    args.retain(|a| a != "--no_stats");
    args
}

/// The files `dir` holds, sorted.
fn files_in(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// The rows of one `[section]` of a statistics report.
fn summary_section(text: &str, title: &str) -> Vec<String> {
    text.split("\n\n")
        .find(|block| block.starts_with(&format!("[{title}]\n")))
        .unwrap_or_else(|| panic!("no section {title} in {text}"))
        .lines()
        .skip(1)
        .map(str::to_string)
        .collect()
}

/// `--no_stats` wins over `--stats_text` and `--stats_html` (Runner.pm:144
/// tests it first; Stats.pm:230, 277 count nothing): no report is written, the
/// annotation is the same, and an existing report is left alone.
#[test]
fn no_stats_wins_over_stats_text_and_stats_html() {
    let dir = tempfile::tempdir().unwrap();
    let args = args_for(&statistics_dir());
    let out = dir.path().join("out.txt");
    std::fs::write(dir.path().join("out.txt_summary.txt"), "kept\n").unwrap();
    let output = run_with(&args, &["--quiet", "--stats_text", "--stats_html"], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        files_in(dir.path()),
        ["out.txt", "out.txt_summary.txt", "out.txt_warnings.txt"]
    );
    assert_eq!(
        std::fs::read_to_string(dir.path().join("out.txt_summary.txt")).unwrap(),
        "kept\n"
    );
    let plain = dir.path().join("plain.txt");
    let output = run_with(&args, &["--quiet"], &plain);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        data_lines(&std::fs::read_to_string(&out).unwrap()),
        data_lines(&std::fs::read_to_string(&plain).unwrap())
    );
}

/// The rows of the table with element id `<id>_table` of an HTML report, each
/// as `key<TAB>count`, the header row dropped.
fn html_table_rows(html: &str, id: &str) -> Vec<String> {
    let open = format!("<table id=\"{id}_table\"");
    let start = html
        .find(&open)
        .unwrap_or_else(|| panic!("no table {id}_table in {html}"));
    let table = &html[start..];
    let table = &table[..table.find("</table>").unwrap()];
    table
        .split("<tr>")
        .skip(1)
        .filter(|row| row.starts_with("<td>"))
        .map(|row| {
            row.trim_end_matches("</tr>")
                .trim_start_matches("<td>")
                .trim_end_matches("</td>")
                .replace("</td><td>", "\t")
        })
        .collect()
}

/// `--stats_text` writes `<output>_summary.txt` beside the output
/// (BaseRunner.pm:254, 305-309) with the corpus's counts, and `--stats_html`
/// beside it writes `<output>_summary.html` (:311-315) with the same rows; the
/// annotation output is the same byte for byte outside its header.
#[test]
fn stats_text_and_stats_html_write_two_files() {
    let dir = tempfile::tempdir().unwrap();
    let args = without_no_stats(args_for(&statistics_dir()));
    let out = dir.path().join("out.txt");
    let output = run_with(&args, &["--quiet", "--stats_text", "--stats_html"], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        files_in(dir.path()),
        [
            "out.txt",
            "out.txt_summary.html",
            "out.txt_summary.txt",
            "out.txt_warnings.txt"
        ]
    );
    let summary = std::fs::read_to_string(dir.path().join("out.txt_summary.txt")).unwrap();
    let html = std::fs::read_to_string(dir.path().join("out.txt_summary.html")).unwrap();
    for (id, title) in [
        ("var_class", "Variant classes"),
        ("var_cons", "Consequences (most severe)"),
        ("consequences", "Consequences (all)"),
        ("coding", "Coding consequences"),
        ("chr", "Variants by chromosome"),
    ] {
        assert_eq!(
            html_table_rows(&html, id),
            summary_section(&summary, title),
            "{title}"
        );
    }
    assert!(html.contains("<td>Variants processed</td><td>339</td>"));
    assert!(html.contains("<td>Overlapped transcripts</td><td>843</td>"));
    assert!(summary.starts_with("[VEP run statistics]\nVEP version (API)\t 116.3 (116)\n"));
    assert_eq!(
        summary_section(&summary, "General statistics"),
        [
            "Lines of input read\t340",
            "Variants processed\t339",
            "Variants filtered out\t0",
            "Novel / existing variants\t-",
            "Overlapped genes\t98",
            "Overlapped transcripts\t843",
            "Overlapped regulatory features\t-",
        ]
    );
    let expected: Vec<String> = summary_section(
        &std::fs::read_to_string(statistics_dir().join("expected/default.summary.log")).unwrap(),
        "Consequences (all)",
    );
    assert_eq!(summary_section(&summary, "Consequences (all)"), expected);
    let silent = dir.path().join("silent.txt");
    let output = run_with(&args, &["--quiet", "--no_stats"], &silent);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        data_lines(&std::fs::read_to_string(&out).unwrap()),
        data_lines(&std::fs::read_to_string(&silent).unwrap())
    );
}

/// Asserts that `html` is a self-contained statistics report as the harness
/// checks one ([`check_html_report`]).
fn assert_html_report(html: &str) {
    let lacks = check_html_report(html);
    assert!(lacks.is_empty(), "{lacks}\n{html}");
}

/// With neither `--stats_text` nor `--stats_html`, the HTML report is written
/// (Runner.pm:787-790 turns `stats_html` on) as `<output>_summary.html` and no
/// text report is; it is checked for an existing file before the output
/// (Runner.pm:144-148 before :202) and refused without `--force_overwrite`, and
/// a refused output leaves no report behind.
#[test]
fn stats_html_is_the_default_report() {
    let dir = tempfile::tempdir().unwrap();
    let args = without_no_stats(args_for(&statistics_dir()));
    let out = dir.path().join("out.txt");
    let output = run_with(&args, &["--quiet"], &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(
        files_in(dir.path()),
        ["out.txt", "out.txt_summary.html", "out.txt_warnings.txt"]
    );
    let html = std::fs::read_to_string(dir.path().join("out.txt_summary.html")).unwrap();
    assert_html_report(&html);
    assert!(html.contains("<td>Variants processed</td><td>339</td>"));
    assert!(
        html.contains(" id=\"chr\""),
        "the corpus has records, so the bar chart prints"
    );
    assert!(html.contains(" id=\"protein\""));
    assert!(!html.contains(" id=\"sift\"") && !html.contains(" id=\"chr_21\""));
    let expected: Vec<String> = summary_section(
        &std::fs::read_to_string(statistics_dir().join("expected/default.summary.log")).unwrap(),
        "Consequences (all)",
    );
    assert_eq!(html_table_rows(&html, "consequences"), expected);

    let unforced = without_force_overwrite(args.clone());
    let summary = dir.path().join("out.txt_summary.html");
    let refused = run_with(&unforced, &["--quiet"], &out);
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        format!(
            "ERROR: [output] Stats file {} already exists. Specify a different output file with --stats_file, overwrite existing file with --force_overwrite or disable stats with --no_stats",
            summary.display()
        )
    );
    assert_eq!(std::fs::read_to_string(&summary).unwrap(), html);
    std::fs::remove_file(&summary).unwrap();
    let refused = run_with(&unforced, &["--quiet"], &out);
    assert_eq!(refused.status.code(), Some(1));
    assert!(
        text(&refused.stderr).starts_with(&format!(
            "ERROR: [output] Output file {} already exists.",
            out.display()
        )),
        "{}",
        text(&refused.stderr)
    );
    assert!(
        !summary.exists(),
        "a refused output leaves no report behind"
    );
    let silent = dir.path().join("silent.txt");
    let output = run_with(&args, &["--quiet", "--no_stats"], &silent);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(!dir.path().join("silent.txt_summary.html").exists());
    assert_eq!(
        data_lines(&std::fs::read_to_string(&out).unwrap()),
        data_lines(&std::fs::read_to_string(&silent).unwrap())
    );
}

/// A run that writes a report without a FASTA says once, in an `INFO:` line on
/// standard error, that the `Distribution of variants on chromosome` charts are
/// left out, whichever report it writes. The line reaches no file of the run,
/// whose report is the one a `--quiet` run writes, and `--quiet` withholds it
/// as it withholds every notice. A run under `--no_stats` and a run with a
/// FASTA, whose report has the chart, print none.
#[test]
fn a_report_without_a_fasta_says_the_distribution_charts_are_left_out() {
    const NOTICE: &str = "INFO: No FASTA file found; the statistics report leaves out the \"Distribution of variants on chromosome\" charts, which need the chromosome lengths (see --fasta)";
    fn notices(output: &Output) -> Vec<String> {
        text(&output.stderr)
            .lines()
            .filter(|line| line.starts_with("INFO: "))
            .map(str::to_string)
            .collect()
    }
    // The report from its general statistics on: every count and chart, without
    // the rows that name the run's own time, command line and output path.
    fn counts_and_charts(dir: &Path, report: &str) -> String {
        let html = std::fs::read_to_string(dir.join(report)).unwrap();
        html[html.find("<h3 id=\"gen_stats\">").unwrap()..].to_string()
    }
    let dir = tempfile::tempdir().unwrap();
    let args = without_no_stats(args_for(&statistics_dir()));
    for (name, flags) in [
        ("html", &[][..]),
        ("text", &["--stats_text"][..]),
        ("both", &["--stats_text", "--stats_html"][..]),
    ] {
        let output = run_with(&args, flags, &dir.path().join(format!("{name}.txt")));
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        assert_eq!(notices(&output), [NOTICE], "{name}");
        assert!(output.stdout.is_empty(), "{name}");
    }
    let written = files_in(dir.path());
    assert_eq!(
        written,
        [
            "both.txt",
            "both.txt_summary.html",
            "both.txt_summary.txt",
            "both.txt_warnings.txt",
            "html.txt",
            "html.txt_summary.html",
            "html.txt_warnings.txt",
            "text.txt",
            "text.txt_summary.txt",
            "text.txt_warnings.txt",
        ]
    );
    for file in &written {
        let content = std::fs::read_to_string(dir.path().join(file)).unwrap();
        assert!(!content.contains("No FASTA file found"), "{file}");
        assert!(
            !content.contains("Distribution of variants on chromosome"),
            "{file}"
        );
    }

    let quiet = run_with(&args, &["--quiet"], &dir.path().join("quiet.txt"));
    assert_eq!(quiet.status.code(), Some(0), "{}", text(&quiet.stderr));
    assert!(text(&quiet.stderr).is_empty(), "{}", text(&quiet.stderr));
    assert_eq!(
        counts_and_charts(dir.path(), "html.txt_summary.html"),
        counts_and_charts(dir.path(), "quiet.txt_summary.html")
    );

    let no_report = run_with(&args, &["--no_stats"], &dir.path().join("none.txt"));
    assert_eq!(no_report.status.code(), Some(0));
    assert!(notices(&no_report).is_empty(), "{:?}", notices(&no_report));

    let fasta = dir.path().join("ref.fa");
    std::fs::write(&fasta, format!(">21\n{}\n", "N".repeat(60))).unwrap();
    std::fs::write(dir.path().join("ref.fa.fai"), "21\t60\t4\t60\t61\n").unwrap();
    let with_fasta = run_with(
        &args,
        &["--fasta", fasta.to_str().unwrap()],
        &dir.path().join("fasta.txt"),
    );
    assert_eq!(
        with_fasta.status.code(),
        Some(0),
        "{}",
        text(&with_fasta.stderr)
    );
    assert!(
        notices(&with_fasta).is_empty(),
        "{:?}",
        notices(&with_fasta)
    );
    let html = std::fs::read_to_string(dir.path().join("fasta.txt_summary.html")).unwrap();
    assert!(html.contains("<h3 id=\"chr_21\">Distribution of variants on chromosome 21</h3>"));
}

/// `-o stdout` writes the report as `stdout_summary.html` in the working
/// directory, the output name as spelled (BaseRunner.pm:254 reads
/// `output_file` as given, ahead of the `STDOUT` comparison), beside the
/// warnings file named the same way.
#[test]
fn stdout_output_writes_stdout_summary_html_as_spelled() {
    let dir = tempfile::tempdir().unwrap();
    let args = without_no_stats(args_for(&statistics_dir()));
    let output = vep()
        .args(&args)
        .args(["--quiet", "-o", "stdout"])
        .current_dir(dir.path())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(data_lines(&text(&output.stdout)).len() > 300);
    assert_eq!(
        files_in(dir.path()),
        ["stdout_summary.html", "stdout_warnings.txt"]
    );
    let html = std::fs::read_to_string(dir.path().join("stdout_summary.html")).unwrap();
    assert_html_report(&html);
    assert!(html.contains("<td>Output file</td><td><pre>stdout</pre></td>"));
}

/// An existing report is refused before the run with Ensembl VEP's text
/// (BaseRunner.pm:275-276) under the output component, exit 1, the output not
/// created; `--force_overwrite` lifts it. A second `-o STDOUT --stats_text` run
/// in one directory is refused the same way, as Ensembl VEP refuses it.
#[test]
fn existing_summary_file_is_refused_without_force() {
    let dir = tempfile::tempdir().unwrap();
    let args = without_force_overwrite(without_no_stats(args_for(&statistics_dir())));
    let out = dir.path().join("out.txt");
    let summary = dir.path().join("out.txt_summary.txt");
    std::fs::write(&summary, "kept\n").unwrap();
    let refused = run_with(&args, &["--quiet", "--stats_text"], &out);
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(
        text(&refused.stderr).trim_end(),
        format!(
            "ERROR: [output] Stats file {} already exists. Specify a different output file with --stats_file, overwrite existing file with --force_overwrite or disable stats with --no_stats",
            summary.display()
        )
    );
    assert_eq!(std::fs::read_to_string(&summary).unwrap(), "kept\n");
    assert!(!out.exists(), "the refusal comes before the output opens");
    let forced = run_with(
        &args,
        &["--quiet", "--stats_text", "--force_overwrite"],
        &out,
    );
    assert_eq!(forced.status.code(), Some(0), "{}", text(&forced.stderr));
    assert!(std::fs::read_to_string(&summary)
        .unwrap()
        .starts_with("[VEP run statistics]\n"));

    let first = vep()
        .args(&args)
        .args(["--quiet", "--stats_text", "-o", "STDOUT"])
        .current_dir(dir.path())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(first.status.code(), Some(0), "{}", text(&first.stderr));
    assert!(dir.path().join("STDOUT_summary.txt").is_file());
    let second = vep()
        .args(&args)
        .args(["--quiet", "--stats_text", "-o", "STDOUT"])
        .current_dir(dir.path())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(second.status.code(), Some(1));
    assert!(
        text(&second.stderr)
            .starts_with("ERROR: [output] Stats file STDOUT_summary.txt already exists."),
        "{}",
        text(&second.stderr)
    );
    assert!(second.stdout.is_empty());
}

/// With both the report and the output in the way, the refusal names the report
/// (Runner.pm:144-148 opens it in `init`, before `run` opens the output at :202),
/// and neither file is touched.
#[test]
fn existing_summary_and_output_name_the_stats_file_first() {
    let dir = tempfile::tempdir().unwrap();
    let args = without_force_overwrite(without_no_stats(args_for(&statistics_dir())));
    let out = dir.path().join("out.txt");
    let summary = dir.path().join("out.txt_summary.txt");
    std::fs::write(&out, "output kept\n").unwrap();
    std::fs::write(&summary, "summary kept\n").unwrap();
    let refused = run_with(&args, &["--quiet", "--stats_text"], &out);
    assert_eq!(refused.status.code(), Some(1));
    assert!(
        text(&refused.stderr).starts_with(&format!(
            "ERROR: [output] Stats file {} already exists.",
            summary.display()
        )),
        "{}",
        text(&refused.stderr)
    );
    assert_eq!(std::fs::read_to_string(&out).unwrap(), "output kept\n");
    assert_eq!(std::fs::read_to_string(&summary).unwrap(), "summary kept\n");
    // The output alone in the way: its refusal, and no empty report left behind.
    std::fs::remove_file(&summary).unwrap();
    let refused = run_with(&args, &["--quiet", "--stats_text"], &out);
    assert_eq!(refused.status.code(), Some(1));
    assert!(
        text(&refused.stderr).starts_with(&format!(
            "ERROR: [output] Output file {} already exists.",
            out.display()
        )),
        "{}",
        text(&refused.stderr)
    );
    assert!(
        !summary.exists(),
        "a refused output leaves no report behind"
    );
}

/// `--quiet` silences the lines on standard error and nothing else: the report
/// is written with every count (Stats.pm keys the counters on `--no_stats`
/// alone), the same as without the flag.
#[test]
fn quiet_keeps_the_summary_file() {
    let dir = tempfile::tempdir().unwrap();
    let args = without_no_stats(args_for(&statistics_dir()));
    let quiet = dir.path().join("quiet.txt");
    let output = run_with(&args, &["--quiet", "--stats_text"], &quiet);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(output.stderr.is_empty(), "{}", text(&output.stderr));
    let loud = dir.path().join("loud.txt");
    let output = run_with(&args, &["--stats_text"], &loud);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(text(&output.stderr).contains("Processed "));
    let quiet_summary = std::fs::read_to_string(dir.path().join("quiet.txt_summary.txt")).unwrap();
    let loud_summary = std::fs::read_to_string(dir.path().join("loud.txt_summary.txt")).unwrap();
    for title in [
        "General statistics",
        "Variant classes",
        "Consequences (all)",
    ] {
        assert_eq!(
            summary_section(&quiet_summary, title),
            summary_section(&loud_summary, title),
            "{title}"
        );
    }
    assert_eq!(
        summary_section(&quiet_summary, "General statistics")[1],
        "Variants processed\t339"
    );
}

/// `--force_overwrite` covers the three guarded files: the output
/// (BaseRunner.pm:203), the report (:275) and the skipped-variants file (:337);
/// without it each existing file is refused in that order of checks, report
/// before output before skipped file.
#[test]
fn force_overwrite_covers_output_stats_and_skipped_files() {
    let dir = tempfile::tempdir().unwrap();
    let args = without_force_overwrite(without_no_stats(args_for(&statistics_dir())));
    let out = dir.path().join("out.txt");
    let summary = dir.path().join("out.txt_summary.txt");
    let skipped = dir.path().join("skipped.txt");
    for file in [&out, &summary, &skipped] {
        std::fs::write(file, "kept\n").unwrap();
    }
    let extra = [
        "--quiet",
        "--stats_text",
        "--skipped_variants_file",
        skipped.to_str().unwrap(),
    ];
    let refused = run_with(&args, &extra, &out);
    assert!(text(&refused.stderr).starts_with("ERROR: [output] Stats file "));
    std::fs::remove_file(&summary).unwrap();
    let refused = run_with(&args, &extra, &out);
    assert!(text(&refused.stderr).starts_with("ERROR: [output] Output file "));
    std::fs::remove_file(&out).unwrap();
    let refused = run_with(&args, &extra, &out);
    assert!(text(&refused.stderr).starts_with("ERROR: [output] Skipped variants file "));
    assert!(!summary.exists() && !out.exists());
    for file in [&out, &summary, &skipped] {
        std::fs::write(file, "kept\n").unwrap();
    }
    let mut forced = extra.to_vec();
    forced.push("--force_overwrite");
    let output = run_with(&args, &forced, &out);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(data_lines(&std::fs::read_to_string(&out).unwrap()).len() > 30);
    assert!(std::fs::read_to_string(&summary)
        .unwrap()
        .starts_with("[VEP run statistics]\n"));
    assert!(std::fs::read_to_string(&skipped)
        .unwrap()
        .starts_with("[VEP skipped the following variants from "));
}
