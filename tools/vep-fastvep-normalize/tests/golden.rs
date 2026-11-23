// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over the fixtures and compares the output file, stderr and the exit code
//! with the goldens, byte for byte and at several thread counts.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_vep-fastvep-normalize");
const ROOT: &str = env!("CARGO_MANIFEST_DIR");
const THREADS_ENV: &str = "VEP_FASTVEP_NORMALIZE_THREADS";
const USAGE: &str = "usage: vep-fastvep-normalize <in> <out>\n";
const ERROR: &str = "ERROR: [vep-fastvep-normalize] ";

/// The thread bounds every fixture is run under; `None` leaves the host's cores.
const THREADS: [Option<&str>; 4] = [Some("1"), Some("2"), Some("7"), None];

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-fastvep-normalize-golden-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn fixture(name: &str) -> PathBuf {
    Path::new(ROOT).join("tests/fixtures").join(name)
}

fn bytes(name: &str) -> Vec<u8> {
    fs::read(fixture(name)).unwrap()
}

fn command(threads: Option<&str>) -> Command {
    let mut command = Command::new(BIN);
    match threads {
        Some(bound) => command.env(THREADS_ENV, bound),
        None => command.env_remove(THREADS_ENV),
    };
    command
}

fn normalise(input: &Path, output: &Path, threads: Option<&str>) -> Output {
    command(threads).arg(input).arg(output).output().unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The fixture `<name>.tab` gives `expected_<name>.tab` and `expected_<name>.stderr` under
/// every thread bound, with exit code 0 and nothing on stdout.
fn held(name: &str) {
    let scratch = Scratch::new(name);
    let written = scratch.path("output.tab");
    for threads in THREADS {
        let output = normalise(&fixture(&format!("{name}.tab")), &written, threads);
        assert_eq!(output.status.code(), Some(0), "{name} {threads:?}");
        assert!(output.stdout.is_empty(), "{name} {threads:?}");
        assert_eq!(
            output.stderr,
            bytes(&format!("expected_{name}.stderr")),
            "{name} {threads:?}: {}",
            stderr(&output)
        );
        assert!(
            fs::read(&written).unwrap() == bytes(&format!("expected_{name}.tab")),
            "{name} {threads:?}: the output differs from the golden"
        );
    }
}

#[test]
fn intergenic_rows_are_rewritten_and_every_other_row_is_kept() {
    held("raw");
}

#[test]
fn every_line_ending_is_written_as_a_line_feed() {
    held("line_endings");
    held("carriage_returns");
}

#[test]
fn a_last_line_that_starts_with_a_hash_keeps_its_missing_ending() {
    held("hash_line_last");
}

#[test]
fn an_input_without_rows_counts_none() {
    held("empty");
    held("header_only");
}

#[test]
fn a_byte_order_mark_makes_the_first_line_a_row() {
    held("byte_order_mark");
    held("byte_order_mark_row");
}

#[test]
fn cells_outside_ascii_pass_through_and_end_no_line() {
    held("beyond_ascii");
}

#[test]
fn an_input_that_is_not_utf_8_exits_1_after_the_lines_before_its_line() {
    let scratch = Scratch::new("not-utf8");
    let written = scratch.path("output.tab");
    for (name, first) in [("not_utf8.tab", 0xff), ("not_utf8_at_the_end.tab", 0xe2)] {
        let input = fixture(name);
        let offset = bytes(name).iter().position(|byte| *byte == first).unwrap();
        for threads in THREADS {
            let output = normalise(&input, &written, threads);
            assert_eq!(output.status.code(), Some(1), "{name} {threads:?}");
            assert!(output.stdout.is_empty());
            assert_eq!(
                stderr(&output),
                format!(
                    "{ERROR}{}: the byte at offset {offset} is not UTF-8\n",
                    input.display()
                )
            );
            assert!(
                fs::read(&written).unwrap() == bytes("expected_raw.tab"),
                "{name} {threads:?}: the output is not the rewrite of the lines before"
            );
        }
    }
}

#[test]
fn any_argument_count_but_two_prints_the_usage_and_exits_1() {
    let scratch = Scratch::new("usage");
    let input = fixture("raw.tab");
    let written = scratch.path("output.tab");
    let (input, output) = (input.as_path(), written.as_path());
    let argument_lists: [&[&Path]; 4] = [
        &[],
        &[input],
        &[input, output, output],
        &[Path::new("--help")],
    ];
    for arguments in argument_lists {
        let output = command(None).args(arguments).output().unwrap();
        assert_eq!(output.status.code(), Some(1), "{arguments:?}");
        assert_eq!(stderr(&output), USAGE, "{arguments:?}");
        assert!(output.stdout.is_empty());
        assert!(!written.exists(), "{arguments:?}");
    }
}

#[test]
fn an_argument_that_looks_like_an_option_is_a_file_name() {
    let scratch = Scratch::new("option-names");
    let output = command(None)
        .current_dir(&scratch.0)
        .args(["-h", "--out"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(stderr(&output).starts_with(&format!("{ERROR}-h: ")));
    fs::copy(fixture("raw.tab"), scratch.path("-h")).unwrap();
    let output = command(None)
        .current_dir(&scratch.0)
        .args(["-h", "--out"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(output.stderr, bytes("expected_raw.stderr"));
    assert!(fs::read(scratch.path("--out")).unwrap() == bytes("expected_raw.tab"));
}

#[test]
fn an_input_that_cannot_be_read_exits_1_and_makes_no_output() {
    let scratch = Scratch::new("no-input");
    let written = scratch.path("output.tab");
    for input in [scratch.path("missing.tab"), scratch.0.clone()] {
        let output = normalise(&input, &written, None);
        assert_eq!(output.status.code(), Some(1), "{}", input.display());
        assert!(
            stderr(&output).starts_with(&format!("{ERROR}{}: ", input.display())),
            "{}",
            stderr(&output)
        );
        assert_eq!(stderr(&output).lines().count(), 1);
        assert!(output.stdout.is_empty());
        assert!(!written.exists(), "{}", input.display());
    }
}

#[test]
fn an_output_that_cannot_be_made_or_written_exits_1() {
    let scratch = Scratch::new("no-output");
    let mut outputs = vec![scratch.path("missing/output.tab"), scratch.0.clone()];
    // A device that refuses every write, where the host has one.
    if Path::new("/dev/full").exists() {
        outputs.push(PathBuf::from("/dev/full"));
    }
    for written in outputs {
        let output = normalise(&fixture("raw.tab"), &written, None);
        assert_eq!(output.status.code(), Some(1), "{}", written.display());
        assert!(
            stderr(&output).starts_with(&format!("{ERROR}{}: ", written.display())),
            "{}",
            stderr(&output)
        );
        assert_eq!(stderr(&output).lines().count(), 1);
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn an_output_that_exists_is_replaced() {
    let scratch = Scratch::new("replaced");
    let written = scratch.path("output.tab");
    fs::write(
        &written,
        [bytes("expected_raw.tab"), bytes("raw.tab")].concat(),
    )
    .unwrap();
    let output = normalise(&fixture("raw.tab"), &written, None);
    assert_eq!(output.status.code(), Some(0));
    assert!(fs::read(&written).unwrap() == bytes("expected_raw.tab"));
}

#[test]
fn one_name_for_the_input_and_the_output_leaves_the_file_empty() {
    let scratch = Scratch::new("one-name");
    let both = scratch.path("raw.tab");
    fs::copy(fixture("raw.tab"), &both).unwrap();
    let output = normalise(&both, &both, None);
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(output.stderr, bytes("expected_empty.stderr"));
    assert!(fs::read(&both).unwrap().is_empty());
}

#[test]
fn an_input_that_is_a_pipe_is_read_as_a_file_is() {
    if !Path::new("/dev/stdin").exists() {
        return;
    }
    let scratch = Scratch::new("pipe");
    let written = scratch.path("output.tab");
    let mut child = command(Some("2"))
        .arg("/dev/stdin")
        .arg(&written)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&bytes("line_endings.tab"))
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", stderr(&output));
    assert_eq!(output.stderr, bytes("expected_line_endings.stderr"));
    assert!(fs::read(&written).unwrap() == bytes("expected_line_endings.tab"));
}

#[test]
fn a_thread_bound_that_is_no_count_exits_1_and_makes_no_output() {
    let scratch = Scratch::new("threads");
    let written = scratch.path("output.tab");
    for bound in ["0", "many", "-2"] {
        let output = normalise(&fixture("raw.tab"), &written, Some(bound));
        assert_eq!(output.status.code(), Some(1), "{bound}");
        assert_eq!(
            stderr(&output),
            format!(
                "{ERROR}{THREADS_ENV}={bound}: expected a whole number of threads, 1 or more\n"
            )
        );
        assert!(!written.exists(), "{bound}");
    }
    let output = normalise(&fixture("raw.tab"), &written, Some(""));
    assert_eq!(output.status.code(), Some(0));
    assert!(fs::read(&written).unwrap() == bytes("expected_raw.tab"));
}

/// The rewrite written out the long way: universal newlines first, then every line split
/// into its cells.
fn reference(input: &[u8]) -> (Vec<u8>, u64, u64) {
    let mut text = Vec::with_capacity(input.len());
    let mut at = 0;
    while at < input.len() {
        if input[at] == b'\r' {
            text.push(b'\n');
            if input.get(at + 1) == Some(&b'\n') {
                at += 1;
            }
        } else {
            text.push(input[at]);
        }
        at += 1;
    }
    let mut out = Vec::with_capacity(text.len() + 1);
    let (mut rows, mut rewritten) = (0, 0);
    for line in text.split_inclusive(|byte| *byte == b'\n') {
        if line.starts_with(b"#") {
            out.extend_from_slice(line);
            continue;
        }
        rows += 1;
        let cells = line.strip_suffix(b"\n").unwrap_or(line);
        let mut cells: Vec<&[u8]> = cells.split(|byte| *byte == b'\t').collect();
        if cells.len() >= 6 && cells[3] == b"-" && cells[4] == b"-" && cells[5] == b"Transcript" {
            cells[5] = b"-";
            rewritten += 1;
        }
        out.extend_from_slice(&cells.join(&b'\t'));
        out.push(b'\n');
    }
    (out, rows, rewritten)
}

/// About seven runs of rows under every line ending, with `#` lines between them, one row
/// longer than two runs and a last row without an ending.
fn generated() -> Vec<u8> {
    let mut input = b"## generated\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tExtra\n".to_vec();
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut row = 0u64;
    while input.len() < 7 << 20 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        row += 1;
        let cells = match state % 5 {
            0 | 1 => "-\t-\tTranscript",
            2 => "-\t-\t-",
            3 => "ENSG00000000001\tENST00000000011\tTranscript",
            _ => "-\tENSR00000000051\tRegulatoryFeature",
        };
        let extra = "x".repeat((state >> 8) as usize % 300);
        let ending = match (state >> 20) % 16 {
            0 => "\r\n",
            1 => "\r",
            2 => "\n\n",
            3 => "\n#a line between rows\n",
            _ => "\n",
        };
        input.extend_from_slice(
            format!(
                "var{row}\t21:{row}\tG\t{cells}\tterm{}\tnote={extra}{ending}",
                state % 7
            )
            .as_bytes(),
        );
        if row == 9_000 {
            let long = "y".repeat((2 << 20) + 4_321);
            input.extend_from_slice(
                format!("long\t21:0\tG\t-\t-\tTranscript\t{long}\r\n").as_bytes(),
            );
        }
    }
    input.extend_from_slice(b"last\t21:1\tG\t-\t-\tTranscript");
    input
}

#[test]
fn a_generated_input_of_several_runs_is_rewritten_alike_under_every_thread_bound() {
    let scratch = Scratch::new("generated");
    let input = scratch.path("generated.tab");
    let written = scratch.path("output.tab");
    let generated = generated();
    let (expected, rows, rewritten) = reference(&generated);
    assert!(rewritten > 10_000 && rows > rewritten);
    fs::write(&input, &generated).unwrap();
    for threads in THREADS {
        let output = normalise(&input, &written, threads);
        assert_eq!(output.status.code(), Some(0), "{threads:?}");
        assert_eq!(
            stderr(&output),
            format!("normalize: rewrote {rewritten}/{rows} intergenic Feature_type=Transcript \u{2192} -\n"),
            "{threads:?}"
        );
        assert!(
            fs::read(&written).unwrap() == expected,
            "{threads:?}: the output differs from the reference rewrite"
        );
    }
}

#[test]
fn the_reference_rewrite_gives_the_goldens() {
    for name in [
        "raw",
        "line_endings",
        "carriage_returns",
        "hash_line_last",
        "empty",
        "header_only",
        "byte_order_mark",
        "byte_order_mark_row",
        "beyond_ascii",
    ] {
        let (out, rows, rewritten) = reference(&bytes(&format!("{name}.tab")));
        assert!(out == bytes(&format!("expected_{name}.tab")), "{name}");
        assert_eq!(
            format!("normalize: rewrote {rewritten}/{rows} intergenic Feature_type=Transcript \u{2192} -\n").into_bytes(),
            bytes(&format!("expected_{name}.stderr")),
            "{name}"
        );
    }
}
