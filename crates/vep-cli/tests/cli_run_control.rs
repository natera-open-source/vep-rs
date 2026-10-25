// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Run control through the `vep` binary: help and version, the output file and
//! `STDOUT`, the cache resolved before anything runs, the FASTA and the HGVS
//! flags settled before the output opens, the warnings file and the
//! skipped-variants file, `--quiet`, `--verbose`, `--safe`, the plugin loading
//! rules, and the `runs` every golden-corpus manifest declares.

mod common;

use std::collections::BTreeSet;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use noodles_bgzf::gzi;

use common::{corpora, execute_manifest_run, manifest_runs, warning_entries, Corpus, WarningEntry};

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
/// flag refused before parsing (`--af_exac`), an unknown flag, an option Ensembl
/// VEP has and this version lacks (`--stats_html`), a bare word, a refused value
/// (`--fork x`), an integer option left without its value (`--fork --help`,
/// whose `--help` Getopt::Long pushes back) or a value flag of this version's
/// own left without its value (`--pipeline_depth --help`, which Ensembl VEP
/// leaves in `@ARGV` whole and prints its usage after), in either order, exits 0
/// with the text of `--help` alone, as Ensembl VEP 116.2 prints its usage and
/// exits 0 on each of these lines. `--help=1` is refused on both engines
/// (vep:30; exit 2 here, in clap's words), and so is `--gff --help`: `gff=s`
/// (Config.pm:260) takes `--help` as its value on both engines, so no usage
/// prints, and this version refuses the option by name (exit 2).
#[test]
fn help_and_version_win_over_a_flag_refused_before_parsing() {
    let expected = text(&vep().arg("--help").output().unwrap().stdout);
    for argv in [
        &["--af_exac", "--help"][..],
        &["--help", "--af_exac"][..],
        &["--bogus", "--help"][..],
        &["--stats_html", "--help"][..],
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
    for argv in [&["--help=1"][..], &["--gff", "--help"][..]] {
        let output = vep().args(argv).output().unwrap();
        assert_eq!(
            output.status.code(),
            Some(2),
            "{argv:?}: {}",
            text(&output.stderr)
        );
        assert!(text(&output.stderr).starts_with("error:"), "{argv:?}");
        assert!(output.stdout.is_empty(), "{argv:?}: no usage");
    }
}

/// This version's own usage spellings, which Ensembl VEP 116.2 refuses or lacks
/// (`-h` is a prefix of eight of its options, left in `@ARGV` and refused at
/// vep:30; it has no `--version`): `-h` prints the help, `--version` and `-V` the
/// version, each exit 0 and each winning over the rest of the line as `--help`
/// does, `--version` ahead of `--help` included, where the help is printed.
#[test]
fn short_help_and_version_are_this_versions_own_and_win_the_same_way() {
    let help = text(&vep().arg("--help").output().unwrap().stdout);
    let version = text(&vep().arg("--version").output().unwrap().stdout);
    for (argv, expected) in [
        (&["--af_exac=1", "-h"][..], &help),
        (&["-h", "--bogus"][..], &help),
        (&["--af_exac", "--version"][..], &version),
        (&["--version", "--af_exac"][..], &version),
        (&["-V", "--af_exac"][..], &version),
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
/// read like the first. Ensembl VEP 116.2 exits 25 on each of these inputs with
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
/// to the rows a `VCFv4.2` header gives, with no warning, as Ensembl VEP 116.2
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
/// that the co-location columns stay empty, where Ensembl VEP annotates
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
         columns stay empty"
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
/// is written as a file named `0` and the empty value fails the run's end with
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
/// below the start) and then `validate_vf`'s coordinate test over the start Perl
/// reads, `POS + 1`, and the END text as written (Parser.pm:475-477, the text kept
/// by BaseVCF4.pm:307-308). Both lines are Ensembl VEP 116.2's for the record,
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
/// `...`, the bytes Ensembl VEP 116.2 writes for this line.
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
/// 116.2 picks `C ENST00000339775`, the row `rank` alone picks, for the record
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
/// Ensembl VEP 116.2 picks `A ENST00000339775 synonymous_variant`, its first row,
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
/// worker count it resolved and the shape it read it from: the run-control
/// input is a plain VCF, so the shape is unknown and the count the `0` rule's;
/// a run with `0` reports no resolution.
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
    let host_capped = run_control(
        &[
            "-o",
            out.to_str().unwrap(),
            "--verbose",
            "--decompression_threads",
            "0",
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
/// false and writes no `VAR_SYNONYMS`. On `21 14967531 T C` Ensembl VEP 116.2
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
    let prefix = "WARNING: Failed to instantiate plugin LoFtool: LoFtool: failed to open /no/such/scores.txt:";
    assert!(
        text(&warned.stderr).lines().any(|l| l.starts_with(prefix)),
        "{}",
        text(&warned.stderr)
    );
    let warnings = std::fs::read_to_string(dir.path().join("out.txt_warnings.txt")).unwrap();
    assert!(
        warnings.lines().any(|l| l.starts_with(prefix)),
        "{warnings}"
    );
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
    assert!(
        errors[0].starts_with(&prefix.replacen("WARNING: ", "ERROR: [plugin] ", 1)),
        "{}",
        errors[0]
    );
    assert!(
        !out.exists(),
        "no output is written when --safe refuses a plugin"
    );
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
    assert!(
        plugin_lines[1].starts_with(
            "WARNING: Failed to instantiate plugin LoFtool: LoFtool: failed to open /no/such/scores.txt:"
        ),
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
/// syntax errors (exit 2) naming the flag. Ensembl VEP 116.2 accepts them, prints
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

/// A warnings file that cannot be created fails the run with Ensembl VEP's text
/// (BaseVEP.pm:894) under the output component, exit 1; under `--quiet` the
/// error line is all that reaches standard error, no warning is echoed. The run
/// ends after the batch that warned, its rows written (here the corpus is one
/// batch, so every row of a run without the flag), and without `--quiet` each
/// warning is echoed on standard error before the error line. Ensembl VEP 116.2
/// dies inside the first warning, with the header alone written and no warning
/// echoed, exit 2; `docs/cli/compatibility.md` states the difference.
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
        &[
            "-o",
            out.to_str().unwrap(),
            "--warning_file",
            "/nonexistent-dir/w.txt",
        ],
        dir.path(),
    );
    assert_eq!(output.status.code(), Some(1));
    let stderr = text(&output.stderr);
    assert_eq!(
        warning_lines(&stderr),
        expected_warnings().lines().collect::<Vec<_>>(),
        "{stderr}"
    );
    assert!(
        stderr
            .trim_end()
            .ends_with("ERROR: [output] Could not write to warnings file /nonexistent-dir/w.txt"),
        "{stderr}"
    );
}

/// A comma-joined or trailing-comma value of a clap-enumerated `%VALID` flag is an
/// argument syntax error (exit 2) in clap's words. Ensembl VEP 116.2 splits the
/// value on commas before its `%VALID` test (Config.pm:946) and reads `b,` by its
/// first letter (OutputFactory.pm:1808), annotating `SIFT=tolerated(0.19)` on
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
/// before any run and with no output file written. Ensembl VEP 116.2 reads such
/// a value as the flag unset:
/// `_is_flag_active` (Config.pm:902-905, `return 0 unless
/// defined($config->{$flag}); return ref($config->{$flag}) eq 'ARRAY' ? scalar
/// @{$config->{$flag}} : $config->{$flag};`) keeps it out of the `%VALID` test
/// (`:945`, `foreach my $flag(grep {$self->_is_flag_active($config, $_)} keys
/// %VALID) {`), and the column gates are truthiness tests too,
/// OutputFactory.pm:1805 (`if (my $opt = $self->{$lc_tool}) {`) for SIFT and
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
            Err(panic) => panic_message(panic.as_ref()),
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

/// A `duckdb` CLI below 1.3.2 on PATH refuses `--output_format parquet` before
/// anything is written: the refusal names the floor and the version found, exit
/// 1, and no output directory, intermediate or Parquet file exists afterwards. A
/// 1.2 release's Parquet COPY rejects `WRITE_BLOOM_FILTER` yet its CLI exits 0,
/// so without the floor such a run would report success and write nothing. The
/// CLI is a stub answering `--version` alone, so the test needs no DuckDB.
#[cfg(unix)]
#[test]
fn parquet_output_refuses_a_duckdb_below_the_floor_before_writing() {
    use std::os::unix::fs::PermissionsExt;
    let dir = tempfile::tempdir().unwrap();
    let bin = dir.path().join("bin");
    std::fs::create_dir(&bin).unwrap();
    let stub = bin.join("duckdb");
    std::fs::write(
        &stub,
        "#!/bin/sh\ncase \"$1\" in --version) echo 'v1.2.2 (Histrionicus) 7c039464e4';; *) exit 0;; esac\n",
    )
    .unwrap();
    std::fs::set_permissions(&stub, std::fs::Permissions::from_mode(0o755)).unwrap();
    let out = dir.path().join("out.parquet");
    let output = vep()
        .args(corpus_args())
        .args([
            "-o",
            out.to_str().unwrap(),
            "--output_format",
            "parquet",
            "--quiet",
        ])
        .env("PATH", &bin)
        .current_dir(dir.path())
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1), "{}", text(&output.stderr));
    assert_eq!(
        text(&output.stderr).trim_end(),
        "ERROR: [run] --output_format=parquet needs duckdb 1.3.2 or newer; found 1.2.2"
    );
    assert!(!out.exists(), "an output directory was created");
    assert!(
        !dir.path().join("out.parquet.tsv.tmp").exists(),
        "an intermediate was written"
    );
    let written: Vec<String> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name != "bin")
        .collect();
    assert!(written.is_empty(), "{written:?}");
}

/// A `--fasta` value that is false in Perl, `0` or the empty value, is a path
/// here (`fasta: args.fasta` in `config.rs` takes it as given, and `open_fasta`
/// in `runner.rs` opens it), refused before any output is written with
/// vep-fasta's `Specified FASTA file/directory <value> not found` (exit 1),
/// where Ensembl VEP 116.2 reads the flag as unset: CacheDir.pm:235
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
/// annotates as `--fasta` does. Ensembl VEP 116.2 filters on `param('assembly')`
/// (BaseVEP.pm:497), which `--offline --full_cache_dir` without `--assembly` or
/// `--grch37` leaves unset (CacheDir.pm:356 fills it inside the layout scan of
/// :303 alone, and BaseRunner.pm:139 from a database `--offline` never opens),
/// so the `/$assembly/i` test of BaseVEP.pm:508 is an empty pattern, which Perl
/// reads as the last successful match, the species test of :507, and every
/// species file passes; :502 then ranks the `toplevel` name first and :513 makes
/// it the run's FASTA, which :523 opens. On such a directory Ensembl VEP 116.2
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
