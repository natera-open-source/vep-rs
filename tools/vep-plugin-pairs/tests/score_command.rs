// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `score` command over small outputs, on the tool's own registry and class table: the
//! refusals of the command line and of the checks before any input is scored, each one
//! error line and exit status 2 with nothing left behind; the report's writing; and the
//! class table a command line that names none reads.

use flate2::write::GzEncoder;
use flate2::Compression;
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::fs;
use std::io::{Read, Write};
use std::os::unix::fs::{FileTypeExt, OpenOptionsExt, PermissionsExt};
use std::path::{Path, PathBuf};
use vep_plugin_pairs::registry::BuiltInRegistry;
use vep_plugin_pairs::score::api::BuiltIn;
use vep_plugin_pairs::score::run;

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
const ERROR: &str = "ERROR: [score_plugin_tuples] ";
/// The tree's class table, which the program holds for a command line that names none.
const COMMITTED_CLASS_TABLE: &str =
    include_str!("../../../scripts/concordance/plugin_class_table.tsv");

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-score-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_string()
    }

    fn file(&self, name: &str, content: impl AsRef<[u8]>) -> String {
        fs::write(self.0.join(name), content).unwrap();
        self.path(name)
    }

    /// A default-format output: the column header and `rows`.
    fn output(&self, name: &str, rows: &[String]) -> String {
        self.file(name, format!("{HEADER}{}", rows.concat()))
    }

    /// A class table that declares no class.
    fn no_class(&self) -> String {
        self.file(
            "classes.tsv",
            "#class\tplugin\tfield names\tmet or not met\n",
        )
    }

    fn exists(&self, name: &str) -> bool {
        self.0.join(name).exists()
    }

    fn is_empty(&self, name: &str) -> bool {
        fs::read_dir(self.0.join(name)).map_or(true, |mut entries| entries.next().is_none())
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn row(location: &str, allele: &str, feature: &str, extra: &str) -> String {
    format!("v\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\tmissense_variant\t-\t-\t-\t-\t-\t-\t{extra}\n")
}

fn reference_rows() -> Vec<String> {
    vec![
        row("1:100", "T", "ENST1", "LOEUF=0.31;IMPACT=LOW"),
        row("1:200", "G", "ENST2", "LOEUF=0.5;Custom=id1"),
    ]
}

fn gzip(content: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(content).unwrap();
    encoder.finish().unwrap()
}

/// The exit status, the lines of standard output and standard error of one command line.
fn run_main(arguments: &[&str]) -> (i32, Vec<String>, String) {
    let arguments: Vec<String> = arguments.iter().map(|word| word.to_string()).collect();
    let (mut out, mut error) = (Vec::new(), Vec::new());
    let code = run(&arguments, &BuiltInRegistry, &BuiltIn, &mut out, &mut error);
    let out = String::from_utf8(out).unwrap();
    (
        code,
        out.lines().map(str::to_string).collect(),
        String::from_utf8(error).unwrap(),
    )
}

/// Asserts one refusal: exit status 2, no summary, one error line that holds `message`.
fn assert_refused(done: &(i32, Vec<String>, String), message: &str) {
    let (code, lines, error) = done;
    assert_eq!((*code, lines.len()), (2, 0), "{error}");
    assert!(
        error.starts_with(ERROR) && error.matches('\n').count() == 1,
        "{error}"
    );
    assert!(error.contains(message), "{error}");
}

#[test]
fn class_options_are_refused_with_one_error_line() {
    let scratch = Scratch::new("class-options");
    let reference = scratch.output("reference.txt", &reference_rows());
    let table = scratch.no_class();
    let missing = scratch.path("missing.tsv");
    let engine = format!("engine={reference}");
    let (report, work) = (scratch.path("report.json"), scratch.path("work"));
    let cases: Vec<(Vec<&str>, String)> = vec![
        (
            vec!["--classes", "planted_engine_zero"],
            "--classes names planted_engine_zero, a class the class table does not".to_string(),
        ),
        (
            vec!["--class-table", &missing],
            format!("not a file: {missing}"),
        ),
        (
            vec!["--class-population", "planted_engine_zero=many"],
            "takes CLASS=N".to_string(),
        ),
        (
            vec!["--class-population", "no_such_class=1"],
            "--class-population names no_such_class".to_string(),
        ),
        (
            vec!["--plugins", "LOEUF,Nope"],
            "--plugins names Nope, no plugin of the registry".to_string(),
        ),
        (vec!["--plugins", ""], "no source to score".to_string()),
        (
            vec!["--engine-binary", "engine=abc"],
            "--engine-binary takes ID=SHA256, not engine=abc".to_string(),
        ),
        (
            vec!["--engine-commit", "engine=not-hex"],
            "--engine-commit takes ID=COMMIT, not engine=not-hex".to_string(),
        ),
        (
            vec!["--engine-commit", "other=0123abc"],
            "other is no engine of this run".to_string(),
        ),
    ];
    for (arguments, message) in cases {
        let mut line = vec![
            "--reference",
            &reference,
            "--engine",
            &engine,
            "--report",
            &report,
            "--work-dir",
            &work,
            "--plugins",
            "LOEUF",
            "--class-table",
            &table,
        ];
        line.extend(arguments);
        assert_refused(&run_main(&line), &message);
        assert!(!scratch.exists("report.json") && !scratch.exists("work"));
    }
}

#[test]
fn run_command_that_no_shell_reads_is_refused() {
    let scratch = Scratch::new("run-command");
    let reference = scratch.output("reference.txt", &reference_rows());
    let engine = format!("engine={reference}");
    let (code, lines, error) = run_main(&[
        "--reference",
        &reference,
        "--engine",
        &engine,
        "--plugins",
        "LOEUF",
        "--report",
        &scratch.path("report.json"),
        "--class-table",
        &scratch.no_class(),
        "--run-command",
        "vep 'a file",
    ]);
    assert_eq!((code, lines.len()), (2, 0));
    assert_eq!(
        error,
        "ERROR: [score_plugin_tuples] --run-command takes the run's command line, not \"vep 'a file\"\n"
    );
}

#[test]
fn every_input_is_checked_before_any_is_scored() {
    let scratch = Scratch::new("first-lines");
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let narrow = scratch.file("narrow.tsv", "1:100\tT\tENST1\n");
    let cut = scratch.file(
        "cut.txt.gz",
        &gzip(fs::read(&reference).unwrap().as_slice())[..20],
    );
    let headed = scratch.file(
        "headed.txt",
        format!("##a comment\n#Uploaded_variation\tLocation\n{}", rows[0]),
    );
    let table = scratch.no_class();
    let (report, work, residuals) = (
        scratch.path("report.json"),
        scratch.path("work"),
        scratch.path("residuals"),
    );
    let scored = |first: &str, later: &str| {
        let engine = format!("engine={reference}");
        let later = format!("later={later}");
        run_main(&[
            "--reference",
            first,
            "--engine",
            &engine,
            "--engine",
            &later,
            "--plugins",
            "LOEUF",
            "--class-table",
            &table,
            "--report",
            &report,
            "--work-dir",
            &work,
            "--residuals",
            &residuals,
        ])
    };
    assert_refused(
        &scored(&reference, &headed),
        "headed.txt: line 2, the column header, names 2 columns, not 14\n",
    );
    assert_refused(
        &scored(&reference, &narrow),
        "narrow.tsv: line 1 has 3 columns; a default-format row has 14, a tuple file line 5\n",
    );
    assert_refused(
        &scored(&reference, &cut),
        "cut.txt.gz: cannot read the gzip stream: Compressed file ended before the end-of-stream marker was reached\n",
    );
    assert_refused(
        &scored(&narrow, &reference),
        "narrow.tsv: line 1 has 3 columns; a default-format row has 14\n",
    );
    assert!(!scratch.exists("work") && !scratch.exists("residuals"));
}

#[test]
fn earlier_report_that_cannot_be_replaced_is_refused_before_scoring() {
    // SAFETY: `geteuid` reads the process's user and cannot fail.
    if unsafe { libc::geteuid() } == 0 {
        // The superuser writes in any directory.
        return;
    }
    for protected in ["the report", "its directory"] {
        let scratch = Scratch::new("locked-report");
        let reference = scratch.output("reference.txt", &reference_rows());
        fs::create_dir(scratch.0.join("reports")).unwrap();
        let report = scratch.file("reports/report.json", b"an earlier report\n");
        let locked = match protected {
            "the report" => PathBuf::from(&report),
            _ => scratch.0.join("reports"),
        };
        let mode = fs::metadata(&locked).unwrap().permissions().mode() & 0o777;
        fs::set_permissions(&locked, fs::Permissions::from_mode(mode & 0o555)).unwrap();
        let engine = format!("engine={reference}");
        let done = run_main(&[
            "--reference",
            &reference,
            "--engine",
            &engine,
            "--plugins",
            "LOEUF",
            "--class-table",
            &scratch.no_class(),
            "--work-dir",
            &scratch.path("work"),
            "--report",
            &report,
        ]);
        fs::set_permissions(&locked, fs::Permissions::from_mode(mode)).unwrap();
        assert_eq!((done.0, done.1.len()), (2, 0), "{protected}");
        assert_eq!(
            done.2,
            format!("{ERROR}cannot write the report {report}\n"),
            "{protected}"
        );
        assert_eq!(fs::read(&report).unwrap(), b"an earlier report\n");
        assert!(!scratch.exists("work"));
    }
}

#[test]
fn unwritable_output_path_is_refused_before_scoring() {
    for (unwritable, shape) in [
        ("--report", "below a file"),
        ("--report", "a directory"),
        ("--residuals", "below a file"),
    ] {
        let scratch = Scratch::new("unwritable");
        let reference =
            scratch.output("reference.txt", &[row("1:100", "T", "ENST1", "Custom=id1")]);
        let below = format!("{reference}/below-a-file");
        let directory = scratch.0.to_str().unwrap().to_string();
        let taken = match shape {
            "below a file" => &below,
            _ => &directory,
        };
        let (mut report, mut residuals) = (scratch.path("report.json"), scratch.path("residuals"));
        match unwritable {
            "--report" => report = taken.clone(),
            _ => residuals = taken.clone(),
        }
        let engine = format!("own={reference}");
        let (code, lines, error) = run_main(&[
            "--reference",
            &reference,
            "--engine",
            &engine,
            "--custom-source",
            "Custom=Custom",
            "--class-table",
            &scratch.no_class(),
            "--work-dir",
            &scratch.path("work"),
            "--report",
            &report,
            "--residuals",
            &residuals,
        ]);
        assert_eq!((code, lines.len()), (2, 0), "{unwritable} {shape}: {error}");
        assert!(
            error.starts_with(ERROR) && error.matches('\n').count() == 1,
            "{error}"
        );
        assert!(!scratch.exists("report.json"));
        assert!(!scratch.exists("work") || unwritable == "--residuals");
        assert!(scratch.is_empty("work"));
    }
}

#[test]
fn command_line_refusals_are_one_error_line() {
    let scratch = Scratch::new("refusals");
    let reference = scratch.output("reference.txt", &[row("1:100", "T", "ENST1", "Custom=id1")]);
    let missing = scratch.path("missing.txt");
    let report = scratch.path("report.json");
    let own = format!("own={reference}");
    let named = format!("a={reference}");
    let cases: Vec<(&str, Vec<&str>, String)> = vec![
        (
            &missing,
            vec!["--engine", &own],
            format!("not a file: {missing}"),
        ),
        (
            &reference,
            vec!["--engine", &reference],
            "--engine takes ID=PATH".to_string(),
        ),
        (
            &reference,
            vec!["--engine", &named, "--engine", &named],
            "--engine names a twice".to_string(),
        ),
        (
            &reference,
            vec!["--engine", &named, "--custom-source", "Custom"],
            "--custom-source takes NAME=PREFIX".to_string(),
        ),
    ];
    for (target, arguments, message) in cases {
        let mut line = vec!["--reference", target, "--report", &report];
        line.extend(arguments);
        assert_refused(&run_main(&line), &message);
        assert!(!scratch.exists("report.json"));
    }
}

#[test]
fn input_that_is_no_vcf_or_cannot_be_read_is_refused_with_nothing_left_behind() {
    let scratch = Scratch::new("inputs");
    let reference = scratch.output("reference.txt", &reference_rows());
    let table = scratch.no_class();
    let empty = scratch.file("empty.vcf", format!("{VCF_HEADER}1\t5\trs1\tC\n"));
    let whole = format!("{VCF_HEADER}1\t100\trs1\tC\tT\t.\t.\t.\n1\t200\trs2\tA\tG\t.\t.\t.\n");
    let packed = gzip(whole.as_bytes());
    let cut = scratch.file("cut.vcf.gz", &packed[..packed.len() - 12]);
    let nowhere = scratch.path("nowhere.vcf");
    let engine = format!("engine={reference}");
    let report = scratch.path("report.json");
    for (name, path, message) in [
        ("missing", &nowhere, format!("not a file: {nowhere}")),
        (
            "empty",
            &empty,
            format!("{empty}: no VCF record, a data line of five or more columns"),
        ),
        ("cut", &cut, format!("{cut}: cannot read the gzip stream: ")),
    ] {
        let (work, residuals) = (format!("work-{name}"), name);
        let done = run_main(&[
            "--reference",
            &reference,
            "--engine",
            &engine,
            "--plugins",
            "LOEUF",
            "--class-table",
            &table,
            "--work-dir",
            &scratch.path(&work),
            "--report",
            &report,
            "--residuals",
            &scratch.path(residuals),
            "--input",
            path,
        ]);
        assert_refused(&done, &format!("{ERROR}{message}"));
        assert!(!scratch.exists("report.json"));
        assert!(
            scratch.is_empty(residuals) && scratch.is_empty(&work),
            "{name}"
        );
    }
}

#[test]
fn threads_that_are_no_whole_number_from_one_are_refused_where_the_passes_start() {
    let scratch = Scratch::new("threads");
    let reference = scratch.output("reference.txt", &reference_rows());
    let engine = format!("engine={reference}");
    for count in ["0", "-5", "99999999999999999999999999"] {
        let done = run_main(&[
            "--reference",
            &reference,
            "--engine",
            &engine,
            "--plugins",
            "LOEUF",
            "--class-table",
            &scratch.no_class(),
            "--work-dir",
            &scratch.path("work"),
            "--report",
            &scratch.path("report.json"),
            "--pairing-threads",
            count,
        ]);
        assert_eq!(
            done.2,
            format!("{ERROR}argument --threads: invalid value: '{count}'\n")
        );
        assert_eq!(done.0, 2);
        assert!(scratch.exists("work") && scratch.is_empty("work"));
    }
}

#[test]
fn command_line_that_cannot_be_read_ends_in_the_message_line_after_the_usage() {
    let (code, lines, error) = run_main(&["--reference", "r", "--bogus"]);
    assert_eq!((code, lines.len()), (2, 0));
    let told: Vec<&str> = error.lines().collect();
    assert_eq!(told.len(), 2, "{error}");
    assert!(told[0].starts_with("usage: vep-plugin-pairs score "));
    assert_eq!(
        told[1],
        "vep-plugin-pairs score: error: the following arguments are required: --engine, --report"
    );
}

#[test]
fn help_prints_the_usage_and_every_option() {
    for word in ["-h", "--help"] {
        let (code, lines, error) = run_main(&[word]);
        assert_eq!((code, error.as_str()), (0, ""));
        assert!(lines[0].starts_with("usage: vep-plugin-pairs score --reference "));
        for option in [
            "--reference",
            "--engine",
            "--engine-binary",
            "--engine-commit",
            "--plugins",
            "--run-command",
            "--custom-source",
            "--work-dir",
            "--report",
            "--residuals",
            "--class-table",
            "--gwas-catalog",
            "--revel-table",
            "--revel-no-match",
            "--utrannotator-sweep",
            "--json-cache",
            "--spliceai-table",
            "--input",
            "--classes",
            "--class-population",
            "--pairing-tool",
            "--pairing-threads",
            "--sort-memory",
            "--sort-threads",
        ] {
            assert!(
                lines
                    .iter()
                    .any(|line| line.trim_start().starts_with(&format!("{option} "))),
                "{option}"
            );
        }
    }
}

#[test]
fn binary_runs_the_score_command_by_its_first_word() {
    let scratch = Scratch::new("binary");
    let done = std::process::Command::new(env!("CARGO_BIN_EXE_vep-plugin-pairs"))
        .args(["score", "--reference"])
        .arg(Path::new(&scratch.path("missing.txt")))
        .output()
        .unwrap();
    assert_eq!(done.status.code(), Some(2));
    assert!(done.stdout.is_empty());
    let error = String::from_utf8(done.stderr).unwrap();
    assert!(
        error.ends_with(
            "vep-plugin-pairs score: error: the following arguments are required: --engine, --report\n"
        ),
        "{error}"
    );
}

/// The report's entry of one source, engine and projection.
fn entry<'a>(report: &'a Value, source: &str, engine: &str, projection: &str) -> &'a Value {
    report["scores"]
        .as_array()
        .unwrap()
        .iter()
        .find(|entry| {
            entry["source"] == source
                && entry["engine"] == engine
                && entry["projection"] == projection
        })
        .unwrap_or_else(|| panic!("no entry of {source}, {engine}, {projection}"))
}

/// An entry's raw and adjusted counts: the reference's, the engine's, the intersection.
fn counts(entry: &Value) -> [i64; 6] {
    [
        "perl",
        "rust",
        "intersection",
        "adj_perl",
        "adj_rust",
        "adj_intersection",
    ]
    .map(|key| entry[key].as_i64().unwrap_or(-1))
}

fn read_report(path: &str) -> Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}

fn residual_lines(directory: &str, source: &str, engine: &str, projection: &str) -> Vec<String> {
    let path = format!("{directory}/{source}.{engine}.{projection}.tsv");
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

fn names(directory: impl AsRef<Path>) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(directory)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    names
}

/// The committed class table's rows as their cells, its header and blank lines left out.
fn committed_rows() -> Vec<Vec<&'static str>> {
    let rows: Vec<Vec<&str>> = COMMITTED_CLASS_TABLE
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .map(|line| line.split('\t').collect())
        .collect();
    assert!(
        !rows.is_empty(),
        "the committed class table declares no class"
    );
    rows
}

fn sha256(content: &[u8]) -> String {
    Sha256::digest(content)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

#[test]
fn command_line_writes_the_report() {
    let scratch = Scratch::new("writes");
    let reference = scratch.output(
        "reference.txt",
        &[row("1:100", "T", "ENST1", "Custom=id1;Custom_SIG=p")],
    );
    let engine = scratch.file("engine.tsv", "1:100\tT\tENST1\tCustom\tid1\n");
    let report_path = scratch.path("reports/report.json");
    let residuals = scratch.path("kept/residuals");
    let (code, lines, error) = run_main(&[
        "--reference",
        &reference,
        "--engine",
        &format!("own={reference}"),
        "--engine",
        &format!("other={engine}"),
        "--custom-source",
        "Custom=Custom",
        "--work-dir",
        &scratch.path("work"),
        "--report",
        &report_path,
        "--sort-memory",
        "16M",
        "--sort-threads",
        "1",
        "--residuals",
        &residuals,
    ]);
    assert_eq!((code, error.as_str()), (0, ""));
    let report = read_report(&report_path);
    let ids: Vec<&str> = report["engines"]
        .as_array()
        .unwrap()
        .iter()
        .map(|engine| engine["id"].as_str().unwrap())
        .collect();
    assert_eq!(ids, ["own", "other"]);
    assert_eq!(report["class_table"]["path"], "plugin_class_table.tsv");
    assert_eq!(
        counts(entry(&report, "Custom", "own", "feature")),
        [2, 2, 2, 2, 2, 2]
    );
    assert_eq!(
        counts(entry(&report, "Custom", "other", "allele")),
        [2, 1, 1, 2, 1, 1]
    );
    assert!(lines.contains(
        &"Custom\tother\tallele\tperl=2 rust=1 intersection=1 adj_perl=2 adj_rust=1 adj_intersection=1 text_only_differences=0 residual_rows=1 shared_keys=0"
            .to_string()
    ));
    assert!(scratch.is_empty("work"));
    assert_eq!(
        names(&residuals),
        [
            "Custom.other.allele.tsv",
            "Custom.other.feature.tsv",
            "Custom.own.allele.tsv",
            "Custom.own.feature.tsv"
        ]
    );
    assert_eq!(
        residual_lines(&residuals, "Custom", "other", "feature"),
        ["reference\t1:100\tT\tENST1\tCustom_SIG\tp"]
    );
    assert_eq!(
        residual_lines(&residuals, "Custom", "other", "allele"),
        ["reference\t1:100\tT\tCustom_SIG\tp"]
    );
    assert!(residual_lines(&residuals, "Custom", "own", "feature").is_empty());
}

/// A run of one LOEUF output against itself, with `more` after the options every run
/// takes; the report's path is `report`.
fn own_run(scratch: &Scratch, report: &str, more: &[&str]) -> (i32, Vec<String>, String) {
    let reference = scratch.output("reference.txt", &reference_rows());
    let engine = format!("engine={reference}");
    let work = scratch.path("work");
    let mut line = vec![
        "--reference",
        &reference,
        "--engine",
        &engine,
        "--plugins",
        "LOEUF",
        "--work-dir",
        &work,
        "--sort-memory",
        "16M",
        "--report",
        report,
    ];
    line.extend_from_slice(more);
    run_main(&line)
}

const OWN_FIRST_LINE: &str = "LOEUF\tengine\tfeature\tperl=2 rust=2 intersection=2 ";

#[test]
fn report_replaces_an_earlier_one_and_takes_its_mode() {
    let scratch = Scratch::new("replaces");
    let report = scratch.file("report.json", b"an earlier report\n");
    fs::set_permissions(&report, fs::Permissions::from_mode(0o600)).unwrap();
    let (code, lines, error) = own_run(&scratch, &report, &[]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert!(lines[0].starts_with(OWN_FIRST_LINE), "{}", lines[0]);
    assert_eq!(
        read_report(&report)["plugins"],
        serde_json::json!(["LOEUF"])
    );
    assert_eq!(
        fs::metadata(&report).unwrap().permissions().mode() & 0o777,
        0o600
    );
    assert_eq!(names(&scratch.0), ["reference.txt", "report.json", "work"]);
}

#[test]
fn report_path_that_is_no_regular_file_is_written_as_it_is() {
    for shape in ["a link", "a pipe"] {
        let scratch = Scratch::new("shape");
        let report = scratch.path("report.json");
        if shape == "a link" {
            scratch.file("linked.json", b"an earlier report\n");
            std::os::unix::fs::symlink("linked.json", &report).unwrap();
        } else {
            let name = std::ffi::CString::new(report.as_bytes()).unwrap();
            // SAFETY: `name` is a NUL-terminated path that outlives the call.
            assert_eq!(unsafe { libc::mkfifo(name.as_ptr(), 0o644) }, 0);
        }
        // A pipe takes a writer only while a reader holds it open, and this reader never blocks.
        let mut reader = fs::OpenOptions::new()
            .read(true)
            .custom_flags(libc::O_NONBLOCK)
            .open(&report)
            .unwrap();
        let (code, lines, error) = own_run(&scratch, &report, &[]);
        let mut written = Vec::new();
        let _ = reader.read_to_end(&mut written);
        assert_eq!((code, error.as_str()), (0, ""), "{shape}");
        assert!(lines[0].starts_with(OWN_FIRST_LINE), "{shape}");
        let stated: Value = serde_json::from_slice(&written).unwrap();
        assert_eq!(stated["plugins"], serde_json::json!(["LOEUF"]), "{shape}");
        let kind = fs::symlink_metadata(&report).unwrap().file_type();
        assert!(
            if shape == "a link" {
                kind.is_symlink()
            } else {
                kind.is_fifo()
            },
            "{shape}"
        );
    }
}

#[test]
fn report_that_cannot_be_written_is_named_and_leaves_no_residual_file() {
    // The staged file's name is the report's and a suffix: a name near the longest a
    // directory takes leaves no room for it.
    let long = format!("{}.json", "r".repeat(245));
    for earlier in [None, Some(b"an earlier report\n")] {
        let scratch = Scratch::new("unwritten");
        let report = scratch.path(&long);
        if let Some(content) = earlier {
            fs::write(&report, content).unwrap();
        }
        let residuals = scratch.path("residuals");
        let (code, lines, error) = own_run(&scratch, &report, &["--residuals", &residuals]);
        assert_eq!((code, lines.len()), (2, 0));
        assert_eq!(
            error,
            format!("{ERROR}cannot write the report {report}: File name too long\n")
        );
        assert!(scratch.is_empty("residuals"));
        let mut kept = vec!["reference.txt".to_string(), "residuals".to_string()];
        if earlier.is_some() {
            kept.push(long.clone());
        }
        kept.push("work".to_string());
        kept.sort();
        assert_eq!(names(&scratch.0), kept);
        if let Some(content) = earlier {
            assert_eq!(fs::read(&report).unwrap(), content);
        }
    }
    // A report written as it is, to a device that takes no byte.
    let scratch = Scratch::new("full");
    let report = scratch.path("report.json");
    std::os::unix::fs::symlink("/dev/full", &report).unwrap();
    let residuals = scratch.path("residuals");
    let (code, lines, error) = own_run(&scratch, &report, &["--residuals", &residuals]);
    assert_eq!((code, lines.len()), (2, 0));
    assert_eq!(
        error,
        format!("{ERROR}cannot write the report {report}: No space left on device\n")
    );
    assert!(scratch.is_empty("residuals"));
}

#[test]
fn run_goes_through_off_the_main_thread() {
    let scratch = Scratch::new("thread");
    let report = scratch.path("report.json");
    let done = std::thread::scope(|scope| {
        scope
            .spawn(|| own_run(&scratch, &report, &[]))
            .join()
            .unwrap()
    });
    assert_eq!((done.0, done.2.as_str()), (0, ""));
    assert!(scratch.exists("report.json"));
}

#[test]
fn run_command_and_engine_facts_are_carried_into_the_report() {
    let scratch = Scratch::new("carried");
    let report = scratch.path("report.json");
    let (code, _, error) = own_run(&scratch, &report, &[]);
    assert_eq!((code, error.as_str()), (0, ""));
    let stated = read_report(&report);
    assert!(stated["run_command"].is_null());
    assert!(stated["engines"][0].get("binary_sha256").is_none());
    let digest = "ab".repeat(32);
    let (code, _, error) = own_run(
        &scratch,
        &report,
        &[
            "--run-command",
            "vep --plugin 'LOEUF,a file' --offline",
            "--engine-commit",
            "engine=0123abc",
            "--engine-binary",
            &format!("engine={digest}"),
        ],
    );
    assert_eq!((code, error.as_str()), (0, ""));
    let stated = read_report(&report);
    assert_eq!(
        stated["run_command"],
        serde_json::json!(["vep", "--plugin", "LOEUF,a file", "--offline"])
    );
    let engine = &stated["engines"][0];
    let keys: Vec<&String> = engine.as_object().unwrap().keys().collect();
    assert_eq!(
        keys,
        [
            "id",
            "binary_sha256",
            "commit",
            "path",
            "sha256",
            "format",
            "rows",
            "unclaimed_keys"
        ]
    );
    assert_eq!(engine["binary_sha256"], digest.as_str());
    assert_eq!(engine["commit"], "0123abc");
    assert!(stated["reference"].get("binary_sha256").is_none());
}

#[test]
fn committed_class_table_names_classes_of_the_tool_alone() {
    use vep_plugin_pairs::score::api::Rules;
    let unknown: Vec<&str> = committed_rows()
        .iter()
        .map(|cells| cells[0])
        .filter(|name| *name != "order_only" && !BuiltIn.knows(name))
        .collect();
    assert_eq!(unknown, Vec::<&str>::new());
}

#[test]
fn committed_class_table_declares_no_class_twice() {
    let mut seen = std::collections::BTreeSet::new();
    let twice: Vec<&str> = committed_rows()
        .iter()
        .map(|cells| cells[0])
        .filter(|name| !seen.insert(*name))
        .collect();
    assert_eq!(twice, Vec::<&str>::new());
}

#[test]
fn committed_class_table_marks_every_row_met_or_not_met() {
    let header: Vec<&str> = COMMITTED_CLASS_TABLE
        .lines()
        .next()
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(
        header[..4],
        ["#class", "plugin", "field names", "met or not met"]
    );
    for cells in committed_rows() {
        assert_eq!(cells.len(), header.len(), "{}", cells[0]);
        assert!(matches!(cells[3], "met" | "not met"), "{}", cells[0]);
    }
}

/// The report of a run that scores one LOEUF row against itself under the tool's own
/// classes and registry, with `more` on its command line.
fn default_table_run(scratch: &Scratch, more: &[&str]) -> Value {
    let report = scratch.path("report.json");
    let (code, _, error) = own_run(scratch, &report, more);
    assert_eq!((code, error.as_str()), (0, ""));
    read_report(&report)
}

fn applied_classes(report: &Value) -> std::collections::BTreeSet<String> {
    report["classes"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|item| item["applied"] == true)
        .map(|item| item["class"].as_str().unwrap().to_string())
        .collect()
}

const OTHER_CLASS_TABLE: &str =
    "loeuf_zero_minimum\tLOEUF\tLOEUF\norder_only\tNearestExonJB\tNearestExonJB\tnot met\t-\n";

fn named(names: &[&str]) -> std::collections::BTreeSet<String> {
    names.iter().map(|name| name.to_string()).collect()
}

#[test]
fn command_line_reads_the_committed_class_table_by_default() {
    let scratch = Scratch::new("default-table");
    let report = default_table_run(&scratch, &[]);
    assert_eq!(report["class_table"]["path"], "plugin_class_table.tsv");
    assert_eq!(
        report["class_table"]["sha256"],
        sha256(COMMITTED_CLASS_TABLE.as_bytes()).as_str()
    );
    let declared: Vec<&str> = report["classes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|item| item["class"].as_str().unwrap())
        .collect();
    let committed: Vec<&str> = committed_rows().iter().map(|cells| cells[0]).collect();
    assert_eq!(declared, committed);
}

#[test]
fn command_line_class_table_and_classes_decide_over_the_default() {
    let scratch = Scratch::new("other-table");
    let other = scratch.file("other.tsv", OTHER_CLASS_TABLE);
    let report = default_table_run(
        &scratch,
        &["--class-table", &other, "--classes", "order_only"],
    );
    assert_eq!(report["class_table"]["path"], other.as_str());
    let declared: Vec<&str> = report["classes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|item| item["class"].as_str().unwrap())
        .collect();
    assert_eq!(declared, ["loeuf_zero_minimum", "order_only"]);
    assert_eq!(applied_classes(&report), named(&["order_only"]));
    let none = default_table_run(&scratch, &["--class-table", &other, "--classes", ""]);
    assert!(applied_classes(&none).is_empty());
    let rows = committed_rows();
    let not_met = rows
        .iter()
        .find(|cells| cells[3] == "not met")
        .unwrap_or(rows.last().unwrap())[0];
    let chosen = named(&[rows[0][0], not_met]);
    let listed: Vec<&str> = chosen.iter().map(String::as_str).collect();
    let report = default_table_run(&scratch, &["--classes", &listed.join(",")]);
    assert_eq!(report["class_table"]["path"], "plugin_class_table.tsv");
    assert_eq!(applied_classes(&report), chosen);
    assert!(applied_classes(&default_table_run(&scratch, &["--classes", ""])).is_empty());
}

#[test]
fn command_line_applies_by_default_every_row_but_those_that_read_not_met() {
    let scratch = Scratch::new("met-alone");
    let met: std::collections::BTreeSet<String> = committed_rows()
        .iter()
        .filter(|cells| cells[3] != "not met")
        .map(|cells| cells[0].to_string())
        .collect();
    assert_eq!(applied_classes(&default_table_run(&scratch, &[])), met);
    let other = scratch.file("other.tsv", OTHER_CLASS_TABLE);
    assert_eq!(
        applied_classes(&default_table_run(&scratch, &["--class-table", &other])),
        named(&["loeuf_zero_minimum"])
    );
}

/// A file of the tool's fixture directory.
fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

#[test]
fn a_table_with_a_csi_index_alone_is_read() {
    let scratch = Scratch::new("csi-alone");
    let table = scratch.path("revel.synthetic.tsv.gz");
    fs::copy(fixture("revel.synthetic.tsv.gz"), &table).unwrap();
    fs::copy(
        fixture("revel.synthetic.tsv.gz.csi"),
        format!("{table}.csi"),
    )
    .unwrap();
    let report = scratch.path("report.json");
    let (code, lines, error) = own_run(&scratch, &report, &["--revel-table", &table]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert!(lines[0].starts_with(OWN_FIRST_LINE), "{}", lines[0]);
    assert_eq!(
        read_report(&report)["revel_table"],
        serde_json::json!({
            "path": table,
            "sha256": sha256(&fs::read(&table).unwrap()),
            "no_match": false,
        })
    );
    assert_eq!(
        names(&scratch.0),
        [
            "reference.txt",
            "report.json",
            "revel.synthetic.tsv.gz",
            "revel.synthetic.tsv.gz.csi",
            "work"
        ]
    );
}

#[test]
fn a_table_without_an_index_is_refused_with_the_names_looked_for() {
    let scratch = Scratch::new("no-index");
    let table = scratch.path("revel.synthetic.tsv.gz");
    fs::copy(fixture("revel.synthetic.tsv.gz"), &table).unwrap();
    let stem = scratch.path("revel.synthetic.tsv");
    let report = scratch.path("report.json");
    let done = own_run(&scratch, &report, &["--revel-table", &table]);
    assert_eq!(
        done,
        (
            2,
            vec![],
            format!(
                "{ERROR}--revel-table: cannot read the tabix index of {table}: none of \
                 {table}.csi, {stem}.csi, {table}.tbi, {stem}.tbi exists\n"
            )
        )
    );
    assert!(!scratch.exists("report.json"));
}
