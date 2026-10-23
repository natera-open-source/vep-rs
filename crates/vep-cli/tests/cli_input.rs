// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary opens its input: a gzip or BGZF file is read whatever it
//! is called and however it is given, `--input_data` stands in for the file, the
//! format comes from the first data line unless `--format` names it, `STDIN` is
//! matched in any case and needs `--format`, and a file that does not exist is
//! refused with Ensembl VEP's text; and how the output's command-line header
//! names the run's flags, whichever source gave them.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Sequence records of the release-116 corpus on both of its chromosomes.
const RECORDS: &str = "\
21\t30282169\tstop_ins4_twostops_rev\tC\tCTATT\t.\t.\t.
21\t30348409\tstop_del2_lost_rev\tCTA\tC\t.\t.\t.
21\t30430756\tstop_ins5_newstop_fwd\tT\tTCCTAA\t.\t.\t.
21\t34716478\tpenult_ins3_stop_fwd\tA\tAATA\t.\t.\t.
22\t20061573\tpartial_codon_del3_fwd\tGAAG\tG\t.\t.\t.
";

/// The same five records in Ensembl's format (Parser/VEP_input.pm:109-168): the
/// alleles without the VCF anchor base, an insertion's start one past its end,
/// strand `1`, the VCF ID as the identifier, space-delimited.
const ENSEMBL_RECORDS: &str = "\
21 30282170 30282169 -/TATT 1 stop_ins4_twostops_rev
21 30348410 30348411 TA/- 1 stop_del2_lost_rev
21 30430757 30430756 -/CCTAA 1 stop_ins5_newstop_fwd
21 34716479 34716478 -/ATA 1 penult_ins3_stop_fwd
22 20061574 20061576 AAG/- 1 partial_codon_del3_fwd
";

fn vcf_text() -> String {
    format!("{VCF_HEADER}{RECORDS}")
}

fn write(dir: &Path, name: &str, bytes: &[u8]) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, bytes).unwrap();
    path
}

fn plain_gzip(text: &str) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(text.as_bytes()).unwrap();
    encoder.finish().unwrap()
}

fn bgzf(text: &str) -> Vec<u8> {
    let mut writer = noodles::bgzf::io::Writer::new(Vec::new());
    writer.write_all(text.as_bytes()).unwrap();
    writer.finish().unwrap()
}

/// The binary's arguments after the input: the release-116 cache and the
/// output path, with `extra` last.
fn args(out: &Path, extra: &[&str]) -> Vec<String> {
    let mut args: Vec<String> = vec![
        "-o".into(),
        out.display().to_string(),
        "--json_cache".into(),
        cache().display().to_string(),
        "--assembly".into(),
        "GRCh38".into(),
        "--force_overwrite".into(),
        "--no_stats".into(),
        "--quiet".into(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    args
}

fn vep(input: &[&str], out: &Path, extra: &[&str]) -> Output {
    common::isolated_vep_anywhere()
        .args(input)
        .args(args(out, extra))
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

fn vep_with_stdin(input: &[&str], stdin: &[u8], out: &Path, extra: &[&str]) -> Output {
    let mut child = common::isolated_vep_anywhere()
        .args(input)
        .args(args(out, extra))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut pipe = child.stdin.take().unwrap();
    let bytes = stdin.to_vec();
    let writer = std::thread::spawn(move || pipe.write_all(&bytes));
    let output = child.wait_with_output().unwrap();
    // A child that refused its input before reading it leaves the writer a
    // broken pipe, which the child's own standard error explains.
    let _ = writer.join().unwrap();
    output
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The output's data lines: every row of the default format.
fn rows(out: &Path) -> Vec<String> {
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The rows of a successful run on `input` with `extra` flags.
fn annotate(dir: &Path, input: &[&str], extra: &[&str]) -> Vec<String> {
    let out = output_path(dir);
    let output = vep(input, &out, extra);
    assert!(
        output.status.success(),
        "{input:?} {extra:?}: {}",
        stderr(&output)
    );
    rows(&out)
}

/// An output path no earlier run of the process wrote.
fn output_path(dir: &Path) -> PathBuf {
    static COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    dir.join(format!("out-{n}.txt"))
}

/// The rows the plain VCF file yields: the reference every other way of giving
/// the same records is compared with.
fn reference_rows(dir: &Path) -> Vec<String> {
    let plain = write(dir, "plain.vcf", vcf_text().as_bytes());
    let rows = annotate(dir, &["-i", plain.to_str().unwrap()], &[]);
    assert!(rows.len() >= 5, "the plain VCF annotates: {rows:?}");
    rows
}

/// Exit 1 with `line` as the one contract line on standard error.
fn assert_refused(output: &Output, line: &str) {
    let err = stderr(output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let errors: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(errors, [line], "{err}");
}

/// A plain gzip file (one member, no BGZF extra field) is decompressed as
/// Ensembl VEP's `get_compressed_filehandle` decompresses any binary input
/// (Parser.pm:242, Utils.pm:429-453), under any name: here `.vcf.gz`, `.gz`
/// and a name with no suffix at all.
#[test]
fn plain_gzip_input_is_read() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path());
    let bytes = plain_gzip(&vcf_text());
    for name in ["plain.vcf.gz", "records.gz", "records"] {
        let path = write(dir.path(), name, &bytes);
        assert_eq!(
            annotate(dir.path(), &["-i", path.to_str().unwrap()], &[]),
            reference,
            "{name}"
        );
    }
    // Bytes after the last member are ignored, as `gzip -dc` ignores them.
    for (name, trailer) in [
        ("trailer.vcf.gz", &b"GARBAGE TRAILER\n"[..]),
        ("zeros.vcf.gz", &[0u8; 512][..]),
    ] {
        let path = write(
            dir.path(),
            name,
            &[bytes.clone(), trailer.to_vec()].concat(),
        );
        assert_eq!(
            annotate(dir.path(), &["-i", path.to_str().unwrap()], &[]),
            reference,
            "{name}"
        );
    }
}

/// A BGZF file is recognised by its `BC` extra subfield, not its name: `.bgz`,
/// `.vcf.gz` and a suffix-less name read the same, with one decompression
/// thread or several.
#[test]
fn bgzf_input_without_gz_suffix_is_read() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path());
    let bytes = bgzf(&vcf_text());
    for (name, threads) in [
        ("records.bgz", "1"),
        ("records.vcf.gz", "4"),
        ("records", "0"),
    ] {
        let path = write(dir.path(), name, &bytes);
        assert_eq!(
            annotate(
                dir.path(),
                &["-i", path.to_str().unwrap()],
                &["--decompression_threads", threads]
            ),
            reference,
            "{name} with {threads} threads"
        );
    }
}

/// A file concatenating a BGZF stream, its end-of-file block included, and a
/// plain gzip member (`cat records.vcf.bgz more.vcf.gz`) is read whole in either
/// order, and so is one whose blocks resume after the plain member, under one
/// decompression thread and several, as `gzip -dc` reads every member of a
/// stream whatever its kind (Parser.pm:242, Utils.pm:441).
#[test]
fn bgzf_and_plain_gzip_members_of_one_stream_are_all_read() {
    let dir = tempfile::tempdir().unwrap();
    let mut expected = reference_rows(dir.path());
    expected.sort();
    let records: Vec<&str> = RECORDS.lines().collect();
    let first_two = format!("{}\n{}\n", records[0], records[1]);
    let last_three = format!("{}\n{}\n{}\n", records[2], records[3], records[4]);
    let streams = [
        (
            "blocks-then-member.vcf.gz",
            [
                bgzf(&format!("{VCF_HEADER}{first_two}")),
                plain_gzip(&last_three),
            ]
            .concat(),
        ),
        (
            "member-then-blocks.vcf.gz",
            [
                plain_gzip(&format!("{VCF_HEADER}{last_three}")),
                bgzf(&first_two),
            ]
            .concat(),
        ),
        (
            "blocks-member-blocks.vcf.gz",
            [bgzf(VCF_HEADER), plain_gzip(&last_three), bgzf(&first_two)].concat(),
        ),
    ];
    for (name, bytes) in &streams {
        let path = write(dir.path(), name, bytes);
        for threads in ["1", "4"] {
            let mut rows = annotate(
                dir.path(),
                &["-i", path.to_str().unwrap()],
                &["--decompression_threads", threads],
            );
            rows.sort();
            assert_eq!(rows, expected, "{name} with {threads} threads");
        }
    }
}

/// A BGZF stream cut inside a block is a read failure, exit 1 with one line,
/// under one decompression worker and under several: the block reader with
/// workers would otherwise end the stream at the cut and the run would write
/// the records before it as a complete output.
#[test]
fn truncated_bgzf_input_is_a_read_error() {
    let dir = tempfile::tempdir().unwrap();
    let records: Vec<&str> = RECORDS.lines().collect();
    let first_stream = bgzf(&format!("{VCF_HEADER}{}\n{}\n", records[0], records[1]));
    let second_stream = bgzf(&format!("{}\n{}\n{}\n", records[2], records[3], records[4]));
    // The second stream less its 28-byte end-of-file block and the last ten
    // bytes of its data block.
    let cut = second_stream.len() - 28 - 10;
    let mut bytes = first_stream;
    bytes.extend_from_slice(&second_stream[..cut]);
    let path = write(dir.path(), "truncated.vcf.gz", &bytes);
    for threads in ["1", "4", "0"] {
        let out = output_path(dir.path());
        let output = vep(
            &["-i", path.to_str().unwrap()],
            &out,
            &["--decompression_threads", threads],
        );
        let err = stderr(&output);
        assert_eq!(output.status.code(), Some(1), "{threads}: {err}");
        let errors: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
        assert_eq!(errors.len(), 1, "{threads}: {err}");
        assert!(
            errors[0].starts_with("ERROR: [run] Failed to read input"),
            "{threads}: {err}"
        );
    }
}

/// Standard input is sniffed like a file: gzip and BGZF bytes on the pipe are
/// decompressed, with `--format` naming the format as a pipe always must.
#[test]
fn gzip_on_stdin_is_read() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path());
    for (label, bytes) in [
        ("plain text", vcf_text().into_bytes()),
        ("gzip", plain_gzip(&vcf_text())),
        ("bgzf", bgzf(&vcf_text())),
    ] {
        let out = output_path(dir.path());
        let output = vep_with_stdin(&["-i", "STDIN"], &bytes, &out, &["--format", "vcf"]);
        assert!(output.status.success(), "{label}: {}", stderr(&output));
        assert_eq!(rows(&out), reference, "{label}");
    }
}

/// `--input_data` (Config.pm:79, `input_data|id`) is read in place of the input
/// file (Runner.pm:817-821), VCF records or Ensembl records alike, each line a
/// record, with the format detected from its first line; `--id` is the same
/// flag; and the command-line header leaves the data out (Config.pm:732).
#[test]
fn input_data_string_is_annotated() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path());
    let out = output_path(dir.path());
    let output = vep(&["--input_data", RECORDS], &out, &[]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(rows(&out), reference);
    let header = std::fs::read_to_string(&out).unwrap();
    let command_line = header
        .lines()
        .find(|l| l.starts_with("## VEP command-line:"))
        .unwrap();
    assert!(
        !command_line.contains("input_data") && !command_line.contains("stop_ins4"),
        "{command_line}"
    );

    let out = output_path(dir.path());
    let output = vep(&["--id", ENSEMBL_RECORDS], &out, &[]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(rows(&out), reference, "Ensembl records through --id");

    // Given both, the string is read and the file is not opened.
    let out = output_path(dir.path());
    let output = vep(
        &["-i", "/nonexistent/input.vcf", "--input_data", RECORDS],
        &out,
        &[],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(rows(&out), reference, "--input_data beside -i");

    // An empty string and `0` are no input data (Runner.pm:818 reads the string
    // only when it is true in Perl), so the input file stands.
    for no_data in ["", "0"] {
        let output = vep(
            &["-i", "/nonexistent/input.vcf", "--input_data", no_data],
            &out,
            &[],
        );
        assert_refused(
            &output,
            "ERROR: [input] File \"/nonexistent/input.vcf\" does not exist",
        );
    }
}

/// The command-line header is `full_command` over the resolved flags
/// (Config.pm:730-749): a flag of the `--config` file and one of a `VEP_*`
/// variable print beside the command line's own, each under its long name,
/// sorted, and a value of exactly `1` prints as the flag alone (:745). Ensembl
/// VEP 116.2 prints `--chr 21,22 --config [PATH]/flags.ini ... --fork ...
/// --symbol` for this config file with `VEP_FORK=1` in the environment.
#[test]
fn config_file_and_environment_flags_print_in_the_command_line_header() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path());
    let plain = dir.path().join("plain.vcf");
    let config = write(dir.path(), "flags.ini", b"symbol 1\nchr 21,22\n");
    let out = output_path(dir.path());
    let output = common::isolated_vep_anywhere()
        .args(["-i", plain.to_str().unwrap()])
        .args(args(&out, &["--config", config.to_str().unwrap()]))
        .env("VEP_FORK", "1")
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    let header = std::fs::read_to_string(&out).unwrap();
    let command_line = header
        .lines()
        .find(|l| l.starts_with("## VEP command-line:"))
        .unwrap();
    assert_eq!(
        command_line,
        format!(
            "## VEP command-line: vep --assembly GRCh38 --chr 21,22 --config [PATH]/flags.ini \
             --force_overwrite --fork --input_file [PATH]/plain.vcf --json_cache [PATH]/json_cache \
             --no_stats --output_file [PATH]/{} --quiet --symbol",
            out.file_name().unwrap().to_str().unwrap()
        )
    );
    // The file's flags acted: the rows carry the symbols it asked for, which
    // the plain run's rows lack.
    let rows = rows(&out);
    assert_eq!(rows.len(), reference.len());
    assert!(reference.iter().all(|r| !r.contains("SYMBOL=")));
    let symbol_rows = rows.iter().filter(|r| r.contains("SYMBOL=")).count();
    assert!(symbol_rows > 1, "{symbol_rows} rows carry a symbol");
}

/// A command line holding an argument that is not UTF-8 completes. The header is
/// rendered from the parsed tokens, so a program name of such bytes, which no
/// option reads, leaves the run and its header untouched; a non-UTF-8 value of
/// an option is clap's usage error, exit 2, before any output is written.
#[cfg(unix)]
#[test]
fn non_utf8_argument_does_not_abort_the_run() {
    use std::ffi::OsString;
    use std::os::unix::ffi::OsStringExt;
    use std::os::unix::process::CommandExt;

    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path());
    let plain = dir.path().join("plain.vcf");
    let out = output_path(dir.path());
    let output = common::isolated_vep_anywhere()
        .arg0(OsString::from_vec(b"v\xE9p".to_vec()))
        .args(["-i", plain.to_str().unwrap()])
        .args(args(&out, &[]))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(rows(&out), reference);
    let header = std::fs::read_to_string(&out).unwrap();
    assert!(
        header
            .lines()
            .any(|l| l.starts_with("## VEP command-line: vep --assembly GRCh38 ")),
        "{header}"
    );

    let mut non_utf8_out = dir.path().as_os_str().to_owned().into_vec();
    non_utf8_out.extend_from_slice(b"/out-\xE9.txt");
    let non_utf8_out = OsString::from_vec(non_utf8_out);
    // `args` renders its output path as text, which would replace the byte, so
    // the path goes in as its bytes in place of the `-o <path>` pair `args` opens with.
    let shared_args = args(&out, &[]);
    let output = common::isolated_vep_anywhere()
        .args(["-i", plain.to_str().unwrap(), "-o"])
        .arg(&non_utf8_out)
        .args(&shared_args[2..])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let err = stderr(&output);
    assert_eq!(output.status.code(), Some(2), "{err}");
    assert!(err.contains("invalid UTF-8"), "{err}");
    assert!(!Path::new(&non_utf8_out).exists());
}

/// `--format guess`, the default, reads the first data line as Ensembl VEP's
/// `detect_format` does (Parser.pm:368-435): an Ensembl-format file with no
/// telling suffix is parsed as Ensembl, a VCF with a `.txt` suffix as VCF, and
/// `--verbose` reports the detection in Ensembl VEP's words (Parser.pm:147).
#[test]
fn format_is_detected_from_the_first_data_line() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path());
    let ensembl = write(dir.path(), "records.txt", ENSEMBL_RECORDS.as_bytes());
    assert_eq!(
        annotate(dir.path(), &["-i", ensembl.to_str().unwrap()], &[]),
        reference
    );
    let vcf = write(dir.path(), "records.tsv", vcf_text().as_bytes());
    let out = output_path(dir.path());
    let mut verbose = args(&out, &["--verbose"]);
    verbose.retain(|a| a != "--quiet");
    let output = common::isolated_vep_anywhere()
        .args(["-i", vcf.to_str().unwrap()])
        .args(&verbose)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(
        stderr(&output).contains("No input file format specified - detected vcf format"),
        "{}",
        stderr(&output)
    );
    assert_eq!(rows(&out), reference);
}

/// An input whose first data line is of no known shape, or that has no data
/// line at all, is refused with Ensembl VEP's text (Parser.pm:150) before any
/// output is written.
#[test]
fn undetectable_input_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    for (name, text) in [
        ("three-columns.txt", "21\t25587759\t25587769\n"),
        ("header-only.vcf", VCF_HEADER),
        ("empty.vcf", ""),
    ] {
        let path = write(dir.path(), name, text.as_bytes());
        let out = output_path(dir.path());
        let output = vep(&["-i", path.to_str().unwrap()], &out, &[]);
        assert_refused(&output, "ERROR: [input] Can't detect input format");
        assert!(!out.exists(), "{name}: no output for a refused input");
        // Named, the format reads the same file.
        let output = vep(&["-i", path.to_str().unwrap()], &out, &["--format", "vcf"]);
        assert!(output.status.success(), "{name}: {}", stderr(&output));
    }
}

/// A line longer than one line of input may hold is refused under the `input`
/// component, naming its line, by every path that reads it: the VCF reader
/// (the third line of a VCF, the line after its two header lines), format
/// detection (the same file without `--format`, whose first data line is the
/// long one), and the identifier list (`--format id`, line 2). The refusal
/// prints once, and a line of exactly the cap is read.
#[test]
fn line_past_the_cap_is_refused_by_every_reader() {
    let cap = vep_io::input::MAX_LINE_BYTES;
    let refusal = |line: u64| {
        format!(
            "ERROR: [input] line {line} is longer than the {cap}-byte limit for one line of input"
        )
    };
    let dir = tempfile::tempdir().unwrap();
    let out = output_path(dir.path());
    let record_prefix = "21\t30282169\trs1\tC\tA\t.\t.\t";
    let vcf_with_padded_record = |info_bytes: usize| {
        let mut bytes = format!("{VCF_HEADER}{record_prefix}").into_bytes();
        bytes.extend(std::iter::repeat_n(b'x', info_bytes));
        bytes.push(b'\n');
        bytes
    };

    let path = write(
        dir.path(),
        "long-third-line.vcf",
        &vcf_with_padded_record(cap),
    );
    let output = vep(&["-i", path.to_str().unwrap()], &out, &["--format", "vcf"]);
    assert_refused(&output, &refusal(3));
    let output = vep(&["-i", path.to_str().unwrap()], &out, &[]);
    assert_refused(&output, &refusal(3));

    let mut long_id = b"rs1\n".to_vec();
    long_id.extend(std::iter::repeat_n(b'x', cap + 1));
    long_id.push(b'\n');
    let path = write(dir.path(), "long-second-line.txt", &long_id);
    let output = vep(&["-i", path.to_str().unwrap()], &out, &["--format", "id"]);
    assert_refused(&output, &refusal(2));

    // Exactly the cap: an INFO value of the cap less the record's other bytes.
    let path = write(
        dir.path(),
        "line-at-the-cap.vcf",
        &vcf_with_padded_record(cap - record_prefix.len()),
    );
    let output = vep(&["-i", path.to_str().unwrap()], &out, &["--format", "vcf"]);
    assert!(output.status.success(), "{}", stderr(&output));
}

/// `STDIN` names standard input in any case (Parser.pm:235, `uc($file) eq
/// 'STDIN'`).
#[test]
fn stdin_literal_is_case_insensitive() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path());
    for literal in ["stdin", "Stdin", "STDIN"] {
        let out = output_path(dir.path());
        let output = vep_with_stdin(
            &["-i", literal],
            vcf_text().as_bytes(),
            &out,
            &["--format", "vcf"],
        );
        assert!(output.status.success(), "{literal}: {}", stderr(&output));
        assert_eq!(rows(&out), reference, "{literal}");
    }
}

/// A named input file that does not exist is refused with Ensembl VEP's text
/// (Parser.pm:239), under the input component, and nothing is written.
#[test]
fn missing_input_file_is_fatal() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = vep(&["-i", "/nonexistent/input.vcf"], &out, &[]);
    assert_refused(
        &output,
        "ERROR: [input] File \"/nonexistent/input.vcf\" does not exist",
    );
    assert!(!out.exists());
}

/// An input name Perl reads as false, the empty string or `0`, is refused with
/// the text of the parser's first test (Parser.pm:132 `die("ERROR: No input
/// file given\n") unless $hashref->{file}`), which comes before the existence
/// test (:239), and nothing is written.
#[test]
fn false_input_file_name_is_no_input_file() {
    let dir = tempfile::tempdir().unwrap();
    for name in ["", "0"] {
        let out = output_path(dir.path());
        let output = vep(&["-i", name], &out, &[]);
        assert_refused(&output, "ERROR: [input] No input file given");
        assert!(!out.exists(), "{name:?}");
    }
}

/// A region input with no FASTA in reach (none named, none in the cache
/// directory) is refused under the input component with Ensembl VEP's text
/// (Parser/Region.pm:115), whether `--format region` names the format or the
/// default `--format guess` detects it from the line, and nothing is written.
#[test]
fn region_input_without_fasta_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let region = write(dir.path(), "region.txt", b"21:30165901-30165901:1/A\n");
    for format in [&["--format", "region"][..], &[][..]] {
        let out = output_path(dir.path());
        let output = vep(&["-i", region.to_str().unwrap()], &out, format);
        assert_refused(
            &output,
            "ERROR: [input] Cannot use Region format in offline mode without a FASTA file",
        );
        assert!(!out.exists(), "{format:?}: no output for a refused input");
    }
}

/// Standard input under `--format guess` is refused before anything is read,
/// with Ensembl VEP's text (Parser.pm:377), whether the input is named `STDIN`
/// or left to the default.
#[test]
fn stdin_without_format_is_fatal() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let output = vep_with_stdin(&["-i", "STDIN"], vcf_text().as_bytes(), &out, &[]);
    assert_refused(
        &output,
        "ERROR: [input] Cannot detect format from STDIN - specify format with --format [format]",
    );
    assert!(!out.exists());
    let output = common::isolated_vep_anywhere()
        .args(args(&out, &[]))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_refused(
        &output,
        "ERROR: [input] Cannot detect format from STDIN - specify format with --format [format]",
    );
}

/// A BGZF file followed by bytes that open no gzip member, zero padding or a
/// text trailer of at least a block header's length, is read whole, as `gzip -dc`
/// reads it for Ensembl VEP (Utils.pm:441; GNU gzip ends at trailing zeros and
/// at trailing garbage after the last member), under one decompression worker
/// and several; the block reader alone would refuse the trailer as a frame.
#[test]
fn bgzf_input_with_a_trailer_is_read_whole() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_rows(dir.path());
    let bytes = bgzf(&vcf_text());
    for (name, trailer) in [
        ("zeros.vcf.gz", &[0u8; 512][..]),
        ("trailer.vcf.gz", &b"GARBAGE TRAILER BYTES\n"[..]),
    ] {
        assert!(
            trailer.len() >= 18,
            "{name}: at least a block header's length"
        );
        let path = write(
            dir.path(),
            name,
            &[bytes.clone(), trailer.to_vec()].concat(),
        );
        for threads in ["1", "4"] {
            assert_eq!(
                annotate(
                    dir.path(),
                    &["-i", path.to_str().unwrap()],
                    &["--decompression_threads", threads]
                ),
                reference,
                "{name} with {threads} threads"
            );
        }
    }
}

/// `--input_data` given with a Perl-false value (the empty string or `0`) beside
/// `-i` reads the file (Runner.pm:818 tests the string's truth) but is a defined
/// parameter to `skipped_variant_msg` (BaseVEP.pm:854-855 tests `defined`), so
/// a skipped record is reported as `variant skipped` with no line number, while
/// the skipped-variants file still carries the record's line number (:858-862).
#[test]
fn input_data_false_value_drops_the_line_number_from_skip_warnings() {
    let dir = tempfile::tempdir().unwrap();
    let records = format!("{VCF_HEADER}{RECORDS}21\t30165901\tinvalid_alleles\tX\tY\t.\t.\t.\n");
    let input = write(dir.path(), "with-skip.vcf", records.as_bytes());
    let reference = annotate(dir.path(), &["-i", input.to_str().unwrap()], &[]);
    let described = "(21 30165901 invalid_alleles X Y . . .): Invalid allele string X/Y or possible parsing error";
    for value in ["", "0"] {
        let out = output_path(dir.path());
        let skipped = dir.path().join(format!("skipped-{}.txt", value.len()));
        let output = vep(
            &["-i", input.to_str().unwrap(), "--input_data", value],
            &out,
            &["--skipped_variants_file", skipped.to_str().unwrap()],
        );
        assert!(output.status.success(), "{value:?}: {}", stderr(&output));
        assert_eq!(rows(&out), reference, "{value:?}: the file is read");
        let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
        assert_eq!(
            warnings,
            format!("WARNING: variant skipped {described}\n"),
            "{value:?}"
        );
        let listed = std::fs::read_to_string(&skipped).unwrap();
        let row = listed.lines().nth(1).unwrap_or_default();
        assert!(row.starts_with("Line 6    \t"), "{value:?}: {listed}");
    }
    let out = output_path(dir.path());
    let output = vep(&["-i", input.to_str().unwrap()], &out, &[]);
    assert!(output.status.success(), "{}", stderr(&output));
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert_eq!(
        warnings,
        format!("WARNING: line 6 skipped {described}\n"),
        "without the flag the line is numbered"
    );
}

/// `--full_cache_dir`, Ensembl VEP's own option (Config.pm:243) that this version
/// carries as an alias of `--json_cache`, prints under its own key in the
/// command-line header, `## VEP command-line:` and the VCF `##VEP-command-line=`
/// alike, whether it came from the command line, a `--config` file or
/// `VEP_FULL_CACHE_DIR`: `full_command` walks the stored keys (Config.pm:734-745)
/// and Ensembl VEP 116.2 prints `--full_cache_dir [PATH]/json_cache` from all
/// three sources.
#[test]
fn full_cache_dir_keeps_its_own_key_in_the_command_line_header() {
    let dir = tempfile::tempdir().unwrap();
    let plain = write(dir.path(), "plain.vcf", vcf_text().as_bytes());
    let cache = cache().display().to_string();
    let config = write(
        dir.path(),
        "flags.ini",
        format!("full_cache_dir {cache}\n").as_bytes(),
    );
    let header_lines = |out: &Path| -> Vec<String> {
        std::fs::read_to_string(out)
            .unwrap()
            .lines()
            .filter(|l| {
                l.starts_with("## VEP command-line:") || l.starts_with("##VEP-command-line=")
            })
            .map(str::to_string)
            .collect()
    };
    for (source, flags, env) in [
        (
            "command line",
            vec!["--full_cache_dir", cache.as_str()],
            None,
        ),
        (
            "config file",
            vec!["--config", config.to_str().unwrap()],
            None,
        ),
        (
            "environment",
            vec![],
            Some(("VEP_FULL_CACHE_DIR", cache.as_str())),
        ),
    ] {
        for format in [&[][..], &["--vcf"][..]] {
            let out = output_path(dir.path());
            let mut command = common::isolated_vep_anywhere();
            command
                .args(["-i", plain.to_str().unwrap(), "-o", out.to_str().unwrap()])
                .args([
                    "--assembly",
                    "GRCh38",
                    "--force_overwrite",
                    "--no_stats",
                    "--quiet",
                ])
                .args(&flags)
                .args(format)
                .stdin(Stdio::null());
            if let Some((key, value)) = env {
                command.env(key, value);
            }
            let output = command.output().unwrap();
            assert!(
                output.status.success(),
                "{source} {format:?}: {}",
                stderr(&output)
            );
            let lines = header_lines(&out);
            assert_eq!(lines.len(), 1, "{source} {format:?}: {lines:?}");
            assert!(
                lines[0].contains("--full_cache_dir [PATH]/json_cache"),
                "{source} {format:?}: {}",
                lines[0]
            );
            assert!(
                !lines[0].contains("--json_cache"),
                "{source} {format:?}: {}",
                lines[0]
            );
            assert!(
                rows(&out).len() >= 5,
                "{source} {format:?}: the cache resolved"
            );
        }
    }
}

/// Under `--full_cache_dir` with no `--assembly`, the JSON `assembly_name` is the
/// cache's own, read from `info.json`. Ensembl VEP 116.2 prints `null` there:
/// JSON.pm:242 reads the `assembly` parameter, set by `--assembly` or by the
/// cache directory's name (CacheDir.pm:356), which the `--full_cache_dir` branch
/// (:303) never reaches, or by `cache_assembly`, declared at JSON.pm:168 and set
/// nowhere. The cache's assembly is the one the run annotates against, so it
/// is printed rather than the unset parameter.
#[test]
fn full_cache_dir_json_names_the_caches_assembly_without_the_flag() {
    let dir = tempfile::tempdir().unwrap();
    let plain = write(dir.path(), "plain.vcf", vcf_text().as_bytes());
    let info: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(cache().join("info.json")).unwrap()).unwrap();
    let assembly = info["assembly"].as_str().unwrap();
    assert!(!assembly.is_empty(), "{info}");
    let out = output_path(dir.path());
    let output = common::isolated_vep_anywhere()
        .args(["-i", plain.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--full_cache_dir", cache().to_str().unwrap()])
        .args(["--json", "--force_overwrite", "--no_stats", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    let records: Vec<serde_json::Value> = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(records.len(), 5, "one object per record");
    for record in &records {
        assert_eq!(record["assembly_name"], assembly, "{record}");
    }
}
