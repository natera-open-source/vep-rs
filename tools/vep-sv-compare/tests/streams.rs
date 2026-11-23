// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary on a stdout or stderr that cannot be written, and on one file for both: the
//! exit code, what the other stream carries, whether the reports are written, and where the
//! lines of stderr stand among the text of stdout.
//!
//! The comparisons read the committed fixture `tests/fixtures/reference_release`, or 320
//! copies of its one input under long names and without outputs where stdout has to outgrow
//! its first block. A file-size limit is set by a shell that then runs the binary.

use std::ffi::OsString;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use vep_sv_compare::cli::streams::{block_bytes, MIN_BLOCK_BYTES};
use vep_sv_compare::cli::{usage, HELP};
use vep_sv_compare::threads::THREADS_ENV;
use vep_sv_compare::TempDir;

/// The four report files of a comparison.
const REPORT_FILES: [&str; 4] = [
    "concordance_report.json",
    "concordance_report.md",
    "discordant.tsv",
    "discordant_adjusted.tsv",
];

/// What the interpreter writes to stderr as it ends with a stdout whose reader has gone.
const STDOUT_NOT_FLUSHED: &str =
    "Exception ignored while flushing sys.stdout:\nBrokenPipeError: [Errno 32] Broken pipe\n";

/// What the interpreter writes to stderr as it ends with a stdout file at the file-size
/// limit.
const STDOUT_TOO_LARGE: &str =
    "Exception ignored while flushing sys.stdout:\nOSError: [Errno 27] File too large\n";

/// The line of a run that fails on a write past the file-size limit.
const FILE_TOO_LARGE: &str = "ERROR: [vep-sv-compare] [Errno 27] File too large\n";

/// The line of a run that fails on a write to a stdout whose reader has gone.
const BROKEN_PIPE: &str = "ERROR: [vep-sv-compare] [Errno 32] Broken pipe\n";

/// What the comparison prints first under the fixture's cache.
const CACHE_LINE: &str = "cross-chromosome mask active over 2 chromosomes (13 transcripts)\n";

/// The bytes of the first block of stdout on the 320 copies: the chunks of text, each handed
/// on once it holds 8,192 bytes, before the first that does not fit in 131,072.
const FIRST_BLOCK: usize = 123_688;

/// The bytes of stdout on the 320 copies and the fixture's input at the refusal.
const HELD_AT_REFUSAL: usize = 150_832;

/// Where a stream of the binary goes.
#[derive(Clone, Copy)]
enum To {
    /// To the test, which reads it.
    Read,
    /// To a pipe whose reader has gone.
    Gone,
}

impl To {
    fn stdio(self) -> Stdio {
        match self {
            To::Read => Stdio::piped(),
            To::Gone => {
                let (reader, writer) = std::io::pipe().expect("a pipe");
                drop(reader);
                Stdio::from(writer)
            }
        }
    }
}

/// How a run of the binary ended; a stream that was not read is empty.
struct Finished {
    code: Option<i32>,
    stdout: String,
    stderr: String,
}

/// The binary on `arguments`, in the directory `at` where one is given.
fn tool(at: Option<&Path>, arguments: &[OsString]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-sv-compare"));
    command
        .args(arguments)
        .env(THREADS_ENV, "2")
        .stdin(Stdio::null());
    if let Some(directory) = at {
        command.current_dir(directory);
    }
    command
}

fn run(arguments: &[OsString], stdout: To, stderr: To) -> Finished {
    run_at(None, arguments, stdout, stderr)
}

fn run_at(at: Option<&Path>, arguments: &[OsString], stdout: To, stderr: To) -> Finished {
    let output = tool(at, arguments)
        .stdout(stdout.stdio())
        .stderr(stderr.stdio())
        .output()
        .expect("the binary runs");
    Finished {
        code: output.status.code(),
        stdout: String::from_utf8(output.stdout).expect("stdout is UTF-8"),
        stderr: String::from_utf8(output.stderr).expect("stderr is UTF-8"),
    }
}

/// Runs the binary with one file for stdout and stderr: the exit code, the bytes of the
/// file and the block a stream on that file holds.
fn run_shared(at: Option<&Path>, arguments: &[OsString]) -> (Option<i32>, Vec<u8>, usize) {
    use std::os::unix::fs::MetadataExt;
    let scratch = TempDir::new("vep_sv_compare_streams_").expect("a scratch directory");
    let path = scratch.path.join("both.txt");
    let file = File::create(&path).expect("the shared file");
    let block = block_bytes(file.metadata().expect("the shared file's size").blksize());
    let twin = file
        .try_clone()
        .expect("a second handle of the shared file");
    let status = tool(at, arguments)
        .stdout(Stdio::from(file))
        .stderr(Stdio::from(twin))
        .status()
        .expect("the binary runs");
    let bytes = std::fs::read(&path).expect("the shared file is read");
    (status.code(), bytes, block)
}

/// The binary on `arguments` under a file-size limit of `blocks` blocks, set by `ulimit -f`
/// of the shell that runs it: a block is 512 bytes, or 1,024 in some shells.
fn limited(at: Option<&Path>, blocks: usize, arguments: &[OsString]) -> Command {
    let mut command = Command::new("sh");
    command
        .arg("-c")
        .arg("ulimit -f \"$0\" && exec \"$@\"")
        .arg(blocks.to_string())
        .arg(env!("CARGO_BIN_EXE_vep-sv-compare"))
        .args(arguments)
        .env(THREADS_ENV, "2")
        .stdin(Stdio::null());
    if let Some(directory) = at {
        command.current_dir(directory);
    }
    command
}

/// Whether `bytes` are the first `blocks` blocks of `whole`, of either size of a block.
fn cut_at(bytes: &[u8], whole: &[u8], blocks: usize) -> bool {
    [512, 1024]
        .iter()
        .any(|block| bytes.len() == blocks * block)
        && whole.starts_with(bytes)
}

/// The committed fixture.
fn committed() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/reference_release")
}

fn golden(expected: &str, name: &str) -> String {
    std::fs::read_to_string(committed().join(expected).join(name)).expect("a golden file")
}

fn strings(texts: &[&str]) -> Vec<OsString> {
    texts.iter().map(OsString::from).collect()
}

/// The arguments of a comparison against the 116.2 reference of the directories under
/// `inputs`, the reference output under `perl`, the fixture's cache and the reports under
/// `out`.
fn arguments(inputs: &Path, perl: &str, out: &Path) -> Vec<OsString> {
    let mut arguments = strings(&["--reference-release", "116.2", "--assembly", "grch37"]);
    for (option, path) in [
        ("--input-dir", inputs.join("in")),
        ("--perl-dir", inputs.join(perl)),
        ("--rust-dir", inputs.join("rust")),
        ("--vep-rs-cache", committed().join("cache")),
        ("--output-dir", out.to_path_buf()),
    ] {
        arguments.push(option.into());
        arguments.push(path.into());
    }
    arguments
}

/// The reports under `out`, a file that is not there as `None`.
fn reports(out: &Path) -> [Option<String>; 4] {
    REPORT_FILES.map(|name| std::fs::read_to_string(out.join(name)).ok())
}

/// The bytes of the reports under `out`, a file that is not there as `None`.
fn report_bytes(out: &Path) -> [Option<Vec<u8>>; 4] {
    REPORT_FILES.map(|name| std::fs::read(out.join(name)).ok())
}

fn scratch() -> TempDir {
    TempDir::new("vep_sv_compare_streams_").expect("a scratch directory")
}

/// A directory holding `count` copies of the fixture's input and of its three outputs,
/// under the names `f0000` and on.
fn copies(count: usize) -> TempDir {
    let directory = scratch();
    let committed = committed();
    for (sub, name) in [
        ("in", "sv.vcf"),
        ("perl", "sv.txt"),
        ("perl_116", "sv.txt"),
        ("rust", "sv.txt"),
    ] {
        let bytes = std::fs::read(committed.join(sub).join(name)).expect("a fixture file");
        let ending = name.strip_prefix("sv").expect("a file named sv");
        std::fs::create_dir(directory.path.join(sub)).expect("a directory of copies");
        for index in 0..count {
            let copy = directory
                .path
                .join(sub)
                .join(format!("f{index:04}{ending}"));
            std::fs::write(copy, &bytes).expect("a copy");
        }
    }
    directory
}

#[test]
fn a_run_that_succeeds_with_no_reader_of_stdout_writes_its_reports_and_exits_120() {
    let out = scratch();
    let finished = run(
        &arguments(&committed(), "perl_116", &out.path),
        To::Gone,
        To::Read,
    );
    assert_eq!(finished.code, Some(120));
    assert_eq!(finished.stderr, STDOUT_NOT_FLUSHED);
    let expected = REPORT_FILES.map(|name| Some(golden("expected_116_2", name)));
    assert_eq!(reports(&out.path), expected);
}

#[test]
fn a_run_that_succeeds_with_no_reader_of_stderr_exits_0() {
    let out = scratch();
    let finished = run(
        &arguments(&committed(), "perl_116", &out.path),
        To::Read,
        To::Gone,
    );
    assert_eq!(finished.code, Some(0));
    let expected = golden("expected_116_2", "stdout.txt").replace(
        "{{OUTPUT_DIR}}",
        out.path.to_str().expect("a UTF-8 scratch path"),
    );
    assert_eq!(finished.stdout, expected);
    assert!(reports(&out.path).iter().all(Option::is_some));
}

#[test]
fn an_early_refusal_exits_120_where_a_stream_holds_what_it_cannot_write() {
    let directory = scratch();
    let path = |name: &str| directory.path.join(name);
    std::fs::create_dir(path("empty_in")).expect("an input directory");
    std::fs::create_dir_all(path("empty_cache/transcripts")).expect("a cache directory");
    let option = |name: &str, value: PathBuf| [OsString::from(name), value.into()];
    let out = option("--output-dir", path("out"));
    let cache = option("--vep-rs-cache", committed().join("cache"));
    // The arguments, the text stdout holds at the refusal, the refusal, and the exit code
    // with no reader of stdout.
    let cases: [(Vec<OsString>, &str, String, i32); 4] = [
        (
            [option("--input-dir", path("absent")), out.clone()].concat(),
            "",
            format!(
                "ERROR: [compare_sv_concordance] input-dir not found: {}\n",
                path("absent").display()
            ),
            1,
        ),
        (
            [
                option("--input-dir", committed().join("in")),
                out.clone(),
                option("--vep-rs-cache", path("empty_cache")),
            ]
            .concat(),
            "",
            format!(
                "ERROR: [compare_sv_concordance] --vep-rs-cache {} yielded no transcripts; the cross-chromosome mask cannot adjudicate and a zero count would read as a clean run\n",
                path("empty_cache").display()
            ),
            1,
        ),
        (
            [option("--input-dir", path("empty_in")), out.clone()].concat(),
            "",
            format!(
                "ERROR: [compare_sv_concordance] no VCF files found in {}\n",
                path("empty_in").display()
            ),
            1,
        ),
        (
            [option("--input-dir", path("empty_in")), out.clone(), cache].concat(),
            CACHE_LINE,
            format!(
                "ERROR: [compare_sv_concordance] no VCF files found in {}\n",
                path("empty_in").display()
            ),
            120,
        ),
    ];
    for (arguments, held, refusal, code_without_stdout) in cases {
        let read = run(&arguments, To::Read, To::Read);
        assert_eq!(read.code, Some(1), "{refusal}");
        assert_eq!((read.stdout.as_str(), &read.stderr), (held, &refusal));

        let no_stdout = run(&arguments, To::Gone, To::Read);
        assert_eq!(no_stdout.code, Some(code_without_stdout), "{refusal}");
        let closing = if held.is_empty() {
            ""
        } else {
            STDOUT_NOT_FLUSHED
        };
        assert_eq!(no_stdout.stderr, format!("{refusal}{closing}"));

        let no_stderr = run(&arguments, To::Read, To::Gone);
        assert_eq!(no_stderr.code, Some(120), "{refusal}");
        assert_eq!(no_stderr.stdout, held);

        let neither = run(&arguments, To::Gone, To::Gone);
        assert_eq!(neither.code, Some(120), "{refusal}");
    }
}

#[test]
fn the_skipped_record_refusal_exits_120_with_either_reader_gone() {
    let out = scratch();
    let arguments = arguments(&committed(), "perl", &out.path);
    let lines = golden("expected_116_2_refused", "stderr.txt");
    let held = golden("expected_116_2_refused", "stdout.txt");

    let no_stdout = run(&arguments, To::Gone, To::Read);
    assert_eq!(no_stdout.code, Some(120));
    assert_eq!(no_stdout.stderr, format!("{lines}{STDOUT_NOT_FLUSHED}"));

    let no_stderr = run(&arguments, To::Read, To::Gone);
    assert_eq!(no_stderr.code, Some(120));
    assert_eq!(no_stderr.stdout, held);

    assert_eq!(run(&arguments, To::Gone, To::Gone).code, Some(120));
    assert_eq!(reports(&out.path), [None, None, None, None]);
}

#[test]
fn a_usage_error_and_the_help_exit_120_where_their_stream_has_no_reader() {
    let usage_error = strings(&["--assembly", "bad"]);
    let read = run(&usage_error, To::Gone, To::Read);
    assert_eq!(read.code, Some(2));
    assert_eq!(
        read.stderr,
        format!(
            "{}vep-sv-compare: error: argument --assembly: invalid choice: 'bad' (choose from 'grch37', 'grch38')\n",
            usage()
        )
    );
    assert_eq!(run(&usage_error, To::Read, To::Gone).code, Some(120));
    assert_eq!(run(&usage_error, To::Gone, To::Gone).code, Some(120));

    let help = strings(&["--help"]);
    let read = run(&help, To::Read, To::Gone);
    assert_eq!((read.code, read.stdout.as_str()), (Some(0), HELP));
    let no_stdout = run(&help, To::Gone, To::Read);
    assert_eq!(no_stdout.code, Some(120));
    assert_eq!(no_stdout.stderr, STDOUT_NOT_FLUSHED);
    assert_eq!(run(&help, To::Gone, To::Gone).code, Some(120));
}

/// The fixture with a row of two columns after the reference's output, and the line of the
/// failure it ends in.
fn cut_output() -> (TempDir, String) {
    let directory = copies(1);
    let output = directory.path.join("perl_116/f0000.txt");
    let mut text = std::fs::read_to_string(&output).expect("the output copied");
    text.push_str("cut\t21:5\n");
    std::fs::write(&output, text).expect("the output with a cut row");
    let failure = format!(
        "ERROR: [vep-sv-compare] {}: line 8 has 2 columns, expected at least 7. A truncated or malformed engine output must fail the comparison rather than silently reduce the tuple set.\n",
        output.display()
    );
    (directory, failure)
}

#[test]
fn a_failure_exits_120_with_either_reader_gone() {
    let (inputs, failure) = cut_output();
    let out = scratch();
    let arguments = arguments(&inputs.path, "perl_116", &out.path);
    let held = format!("{CACHE_LINE}Found 1 input VCF files\n  Processing f0000...\n");

    let read = run(&arguments, To::Read, To::Read);
    assert_eq!(read.code, Some(1));
    assert_eq!((&read.stdout, &read.stderr), (&held, &failure));

    let no_stdout = run(&arguments, To::Gone, To::Read);
    assert_eq!(no_stdout.code, Some(120));
    assert_eq!(no_stdout.stderr, format!("{failure}{STDOUT_NOT_FLUSHED}"));

    let no_stderr = run(&arguments, To::Read, To::Gone);
    assert_eq!(no_stderr.code, Some(120));
    assert_eq!(no_stderr.stdout, held);
}

/// A directory of 320 copies of the fixture's input under names of 100 characters, with
/// empty directories for the two engines' outputs: every line the comparison prints for
/// them names no path outside the directory when it is run there on [`relative`]
/// arguments. With `refused`, the fixture's input and its outputs of the 115.2 shape stand
/// after the copies.
fn long_names(refused: bool) -> TempDir {
    let directory = scratch();
    let committed = committed();
    for sub in ["in", "perl", "rust"] {
        std::fs::create_dir(directory.path.join(sub)).expect("a directory of the comparison");
    }
    let input = std::fs::read(committed.join("in/sv.vcf")).expect("the fixture's input");
    for index in 0..320 {
        let name = format!("f{index:04}{}.vcf", "x".repeat(95));
        std::fs::write(directory.path.join("in").join(name), &input).expect("a copy");
    }
    if refused {
        for name in ["in/sv.vcf", "perl/sv.txt", "rust/sv.txt"] {
            std::fs::copy(committed.join(name), directory.path.join(name)).expect("a fixture file");
        }
    }
    directory
}

/// The arguments of a comparison against the 116.2 reference run in a directory of
/// [`long_names`].
fn relative() -> Vec<OsString> {
    strings(&[
        "--reference-release",
        "116.2",
        "--input-dir",
        "in",
        "--perl-dir",
        "perl",
        "--rust-dir",
        "rust",
        "--output-dir",
        "out",
    ])
}

#[test]
fn stdout_written_in_mid_run_to_no_reader_fails_the_run_before_its_reports() {
    let inputs = long_names(false);
    let at = Some(inputs.path.as_path());
    let no_stdout = run_at(at, &relative(), To::Gone, To::Read);
    assert_eq!(no_stdout.code, Some(120));
    assert_eq!(
        no_stdout.stderr,
        format!("{BROKEN_PIPE}{STDOUT_NOT_FLUSHED}")
    );
    assert_eq!(reports(&inputs.path.join("out")), [None, None, None, None]);
    // With a reader the same run writes all of its text and its reports.
    let read = run_at(at, &relative(), To::Read, To::Gone);
    assert_eq!(read.code, Some(0));
    assert_eq!(read.stdout.len(), 151_867);
    assert!(reports(&inputs.path.join("out"))
        .iter()
        .all(Option::is_some));
}

#[test]
fn one_file_for_both_streams_carries_a_refusal_before_the_text_stdout_held() {
    let out = scratch();
    let arguments = arguments(&committed(), "perl", &out.path);
    let (code, shared, _) = run_shared(None, &arguments);
    assert_eq!(code, Some(1));
    let expected = [
        golden("expected_116_2_refused", "stderr.txt"),
        golden("expected_116_2_refused", "stdout.txt"),
    ]
    .concat();
    assert_eq!(String::from_utf8(shared).expect("UTF-8 text"), expected);
}

#[test]
fn one_file_for_both_streams_carries_a_failure_after_the_text_of_stdout() {
    let (inputs, failure) = cut_output();
    let out = scratch();
    let (code, shared, _) = run_shared(None, &arguments(&inputs.path, "perl_116", &out.path));
    assert_eq!(code, Some(1));
    assert_eq!(
        String::from_utf8(shared).expect("UTF-8 text"),
        format!("{CACHE_LINE}Found 1 input VCF files\n  Processing f0000...\n{failure}")
    );
}

#[test]
fn one_file_for_both_streams_carries_a_refusal_after_the_first_block_of_stdout() {
    let inputs = long_names(true);
    let at = Some(inputs.path.as_path());
    let read = run_at(at, &relative(), To::Read, To::Read);
    assert_eq!(read.code, Some(1));
    assert_eq!(read.stdout.len(), HELD_AT_REFUSAL);
    assert_eq!(read.stderr, golden("expected_116_2_refused", "stderr.txt"));
    let (code, shared, block) = run_shared(at, &relative());
    assert_eq!(code, Some(1));
    // A file that states a larger block of its own holds all of this stdout to the end.
    let first = if block == MIN_BLOCK_BYTES {
        FIRST_BLOCK
    } else {
        0
    };
    let expected = [
        &read.stdout[..first],
        read.stderr.as_str(),
        &read.stdout[first..],
    ]
    .concat();
    assert!(
        String::from_utf8(shared).expect("UTF-8 text") == expected,
        "the lines of the refusal stand after the first {first} bytes of stdout"
    );
}

#[cfg(target_os = "linux")]
#[test]
fn a_device_with_no_space_left_fails_the_flush_in_python_s_words() {
    let full = || {
        let device = File::options()
            .write(true)
            .open("/dev/full")
            .expect("the full device");
        Stdio::from(device)
    };
    let out = scratch();
    let arguments = arguments(&committed(), "perl_116", &out.path);
    let output = tool(None, &arguments)
        .stdout(full())
        .stderr(Stdio::piped())
        .output()
        .expect("the binary runs");
    assert_eq!(output.status.code(), Some(120));
    assert_eq!(
        String::from_utf8(output.stderr).expect("stderr is UTF-8"),
        "Exception ignored while flushing sys.stdout:\nOSError: [Errno 28] No space left on device\n"
    );
    assert!(reports(&out.path).iter().all(Option::is_some));
    // A refusal that stderr cannot take: stdout is written whole, and the exit code is 120.
    let refused = tool(None, &self::arguments(&committed(), "perl", &out.path))
        .stdout(Stdio::piped())
        .stderr(full())
        .output()
        .expect("the binary runs");
    assert_eq!(refused.status.code(), Some(120));
    assert_eq!(
        String::from_utf8(refused.stdout).expect("stdout is UTF-8"),
        golden("expected_116_2_refused", "stdout.txt")
    );
}

#[test]
fn a_stdout_file_at_the_file_size_limit_ends_the_run_with_120() {
    let directory = scratch();
    let file = directory.path.join("help.txt");
    for blocks in [0, 1] {
        let output = limited(None, blocks, &strings(&["--help"]))
            .stdout(File::create(&file).expect("the file of stdout"))
            .stderr(Stdio::piped())
            .output()
            .expect("the shell runs");
        assert_eq!(output.status.code(), Some(120), "{blocks} blocks");
        assert_eq!(
            String::from_utf8(output.stderr).expect("stderr is UTF-8"),
            STDOUT_TOO_LARGE
        );
        let written = std::fs::read(&file).expect("the file of stdout is read");
        assert!(
            cut_at(&written, HELP.as_bytes(), blocks),
            "{} bytes under {blocks} blocks",
            written.len()
        );
    }
}

#[test]
fn a_stderr_file_at_the_file_size_limit_ends_a_refusal_with_120() {
    let out = scratch();
    let file = out.path.join("stderr.txt");
    let output = limited(
        None,
        0,
        &arguments(&committed(), "perl", &out.path.join("out")),
    )
    .stdout(Stdio::piped())
    .stderr(File::create(&file).expect("the file of stderr"))
    .output()
    .expect("the shell runs");
    assert_eq!(output.status.code(), Some(120));
    assert_eq!(
        String::from_utf8(output.stdout).expect("stdout is UTF-8"),
        golden("expected_116_2_refused", "stdout.txt")
    );
    assert_eq!(
        std::fs::read(&file).expect("the file of stderr is read"),
        b""
    );
}

#[test]
fn a_report_past_the_file_size_limit_fails_the_run_and_is_cut_at_the_limit() {
    let out = scratch();
    let arguments = arguments(&committed(), "perl_116", &out.path);
    let whole = golden("expected_116_2", "concordance_report.json");
    // What stdout holds when the first report is written: the lines of the cache and of
    // the one input.
    let held: String = golden("expected_116_2", "stdout.txt")
        .split_inclusive('\n')
        .take(4)
        .collect();
    for blocks in [0, 2] {
        let output = limited(None, blocks, &arguments)
            .output()
            .expect("the shell runs");
        assert_eq!(output.status.code(), Some(1), "{blocks} blocks");
        assert_eq!(
            String::from_utf8(output.stderr).expect("stderr is UTF-8"),
            FILE_TOO_LARGE
        );
        assert_eq!(
            String::from_utf8(output.stdout).expect("stdout is UTF-8"),
            held
        );
        let [json, markdown, raw, adjusted] = report_bytes(&out.path);
        let json = json.expect("the JSON report is made");
        assert!(
            cut_at(&json, whole.as_bytes(), blocks),
            "{} bytes under {blocks} blocks",
            json.len()
        );
        assert_eq!((markdown, raw, adjusted), (None, None, None));
    }
}

#[test]
fn a_table_past_the_file_size_limit_leaves_the_reports_before_it_whole() {
    // One record whose reference output holds 600 rows the scored engine has none of: the
    // table of one-sided tuples outgrows the JSON report and the markdown.
    let directory = scratch();
    let write = |name: &str, text: &str| {
        let path = directory.path.join(name);
        std::fs::create_dir_all(path.parent().expect("a file under the directory"))
            .expect("a directory of the comparison");
        std::fs::write(path, text).expect("a file of the comparison");
    };
    let head = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n";
    let rows: String = (0..600)
        .map(|index| format!("a1\t21:100\tT\tG\tENST{index:011}\tTranscript\tintron_variant\n"))
        .collect();
    write(
        "in/one.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t100\ta1\tA\tT\t.\t.\t.\n",
    );
    write("perl/one.txt", &format!("{head}{rows}"));
    write("rust/one.txt", head);
    let at = Some(directory.path.as_path());
    let out = directory.path.join("out");

    let unlimited = run_at(at, &relative(), To::Read, To::Read);
    assert_eq!(unlimited.code, Some(0), "stderr: {}", unlimited.stderr);
    let [json, markdown, raw, adjusted] = report_bytes(&out);
    let (json, markdown) = (
        json.expect("the JSON report"),
        markdown.expect("the markdown"),
    );
    let (raw, adjusted) = (
        raw.expect("the raw table"),
        adjusted.expect("the adjusted table"),
    );
    let blocks = 40;
    assert!(json.len() < blocks * 512 && markdown.len() < blocks * 512);
    assert!(raw.len() > blocks * 1024 && adjusted == raw);
    std::fs::remove_dir_all(&out).expect("the reports are removed");

    let output = limited(at, blocks, &relative())
        .output()
        .expect("the shell runs");
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8(output.stderr).expect("stderr is UTF-8"),
        FILE_TOO_LARGE
    );
    let [cut_json, cut_markdown, cut_raw, cut_adjusted] = report_bytes(&out);
    assert_eq!((cut_json, cut_markdown), (Some(json), Some(markdown)));
    let cut_raw = cut_raw.expect("the raw table is made");
    assert!(
        cut_at(&cut_raw, &raw, blocks),
        "{} bytes under {blocks} blocks",
        cut_raw.len()
    );
    assert_eq!(cut_adjusted, None);
}

#[test]
fn stdout_past_the_file_size_limit_in_mid_run_fails_the_run_before_its_reports() {
    let inputs = long_names(false);
    let at = Some(inputs.path.as_path());
    let out = inputs.path.join("out");
    let whole = run_at(at, &relative(), To::Read, To::Read);
    assert_eq!(whole.code, Some(0), "stderr: {}", whole.stderr);
    std::fs::remove_dir_all(&out).expect("the reports are removed");

    // The limit lies inside the first block of stdout: the file takes a part of the block.
    let blocks = 64;
    let file = inputs.path.join("stdout.txt");
    let output = limited(at, blocks, &relative())
        .stdout(File::create(&file).expect("the file of stdout"))
        .stderr(Stdio::piped())
        .output()
        .expect("the shell runs");
    assert_eq!(output.status.code(), Some(120));
    assert_eq!(
        String::from_utf8(output.stderr).expect("stderr is UTF-8"),
        format!("{FILE_TOO_LARGE}{STDOUT_TOO_LARGE}")
    );
    let written = std::fs::read(&file).expect("the file of stdout is read");
    assert!(
        cut_at(&written, whole.stdout.as_bytes(), blocks),
        "{} bytes under {blocks} blocks",
        written.len()
    );
    assert_eq!(reports(&out), [None, None, None, None]);
}
