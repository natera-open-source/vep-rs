// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The files a Parquet run writes beside its output, named at the output path
//! or inside the output directory, through the `vep` binary and, for a
//! program that calls the runner itself, through the library. The run writes
//! its parts into a temporary directory and renames it to the output path as
//! its last step, after a check that takes only what the Parquet writer
//! wrote, so a statistics report, a warnings file or a skipped-variants file
//! placed there would be refused once every record is annotated. Such a run
//! is refused before its input is opened, and for the warnings file before a
//! warning can create it: one line names the file and the flag that places
//! it, and nothing at the output path or inside it is created or changed. A
//! run whose files stand beside the directory completes, and the default
//! names stand there however the output path ends: a path written with a
//! trailing separator, or with a `.` step after it, names its report and its
//! warnings file as the same path written without them does. A text output
//! keeps Ensembl VEP's names on the path as spelled, and the Parquet page's
//! sentences on these files are held to the lines the binary prints.

#![cfg(feature = "parquet")]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// The transcript cache of the release corpus of `assembly`.
fn cache(assembly: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(assembly)
        .join("json_cache")
}

/// Copies the directory `from` to `to`, every file but the one named
/// `except`.
fn copy_without(from: &Path, to: &Path, except: &str) {
    fs::create_dir(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_without(&entry.path(), &target, except);
        } else if entry.file_name() != except {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// A VCF of one made-up substitution per `(chromosome, position)`, written
/// under `dir`.
fn vcf(dir: &Path, name: &str, records: &[(&str, u64)]) -> PathBuf {
    let path = dir.join(name);
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (i, (chrom, pos)) in records.iter().enumerate() {
        text.push_str(&format!("{chrom}\t{pos}\tr{}\tA\tG\t.\t.\t.\n", i + 1));
    }
    fs::write(&path, text).unwrap();
    path
}

/// A scratch directory and, inside it, the directory the runs of one test
/// read and write; the binary's home stands beside that directory.
struct Scratch {
    root: tempfile::TempDir,
}

impl Scratch {
    fn new() -> Self {
        let root = tempfile::tempdir().unwrap();
        fs::create_dir(root.path().join("work")).unwrap();
        Self { root }
    }

    fn work(&self) -> PathBuf {
        self.root.path().join("work")
    }

    fn at(&self, name: &str) -> String {
        self.work().join(name).display().to_string()
    }

    /// A Parquet run over `input` to `output` as spelled, started in the work
    /// directory, with `extra` after its own flags. The statistics report is
    /// on unless `extra` switches it off.
    fn run(&self, input: &Path, output: &str, extra: &[&str]) -> Output {
        self.run_over(&cache("GRCh38"), "GRCh38", input, output, extra)
    }

    /// [`Self::run`] over the transcript cache `cache` of `assembly`.
    fn run_over(
        &self,
        cache: &Path,
        assembly: &str,
        input: &Path,
        output: &str,
        extra: &[&str],
    ) -> Output {
        common::isolated_vep(self.root.path())
            .current_dir(self.work())
            .args(["--json_cache", cache.to_str().unwrap()])
            .args(["--assembly", assembly, "--no_fasta", "--quiet"])
            .args(["--output_format", "parquet"])
            .args(["--input_file", input.to_str().unwrap()])
            .args(["--output_file", output])
            .args(extra)
            .stdin(Stdio::null())
            .output()
            .unwrap()
    }

    /// A completed output of two parts at `name` in the work directory.
    fn earlier_output(&self, input: &Path, name: &str) {
        let output = self.run(input, &self.at(name), &["--no_stats"]);
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    }

    /// Every file under the work directory with its bytes, and every
    /// directory, as paths under it, sorted. A symbolic link is named with
    /// its target and read through it.
    fn contents(&self) -> Vec<(String, Option<Vec<u8>>)> {
        fn walk(dir: &Path, root: &Path, entries: &mut Vec<(String, Option<Vec<u8>>)>) {
            for entry in fs::read_dir(dir).unwrap() {
                let path = entry.unwrap().path();
                let mut name = path.strip_prefix(root).unwrap().display().to_string();
                if let Ok(target) = fs::read_link(&path) {
                    name = format!("{name} -> {}", target.display());
                }
                if path.is_dir() && !path.is_symlink() {
                    entries.push((name, None));
                    walk(&path, root, entries);
                } else {
                    entries.push((name, fs::read(&path).ok()));
                }
            }
        }
        let mut entries = Vec::new();
        walk(&self.work(), &self.work(), &mut entries);
        entries.sort();
        entries
    }

    /// The names under the work directory, sorted.
    fn names(&self) -> Vec<String> {
        self.contents().into_iter().map(|(name, _)| name).collect()
    }
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// The `ERROR:` lines of a run's standard error.
fn error_lines(output: &Output) -> Vec<String> {
    text(&output.stderr)
        .lines()
        .filter(|line| line.starts_with("ERROR:"))
        .map(str::to_string)
        .collect()
}

/// The refusal of `file`, a `what` placed by `flag`, inside `output`.
fn inside(what: &str, file: &str, output: &str, flag: &str) -> String {
    format!(
        "ERROR: [parquet] {what} {file} lies inside the output directory {output}; name it outside that directory with {flag}"
    )
}

/// The refusal of `file`, a `what` placed by `flag`, named as the output path.
fn at_the_output_path(what: &str, file: &str, flag: &str) -> String {
    format!(
        "ERROR: [parquet] {what} {file} is the output path itself; name another file with {flag}"
    )
}

/// `output` is refused with exit 1 and `line` as its one `ERROR:` line, and
/// the work directory holds what it held in `before`, byte for byte.
#[track_caller]
fn assert_refused_leaving(
    scratch: &Scratch,
    output: &Output,
    line: &str,
    before: &[(String, Option<Vec<u8>>)],
) {
    assert_eq!(output.status.code(), Some(1), "{}", text(&output.stderr));
    assert_eq!(error_lines(output), [line]);
    assert_eq!(
        scratch.contents(),
        before,
        "the refused run created, changed or removed a file"
    );
}

/// A statistics report named inside the output directory is refused by name:
/// at a new path, over an earlier output without `--force_overwrite` and with
/// it, before a record out of order is met and before an input that does not
/// exist is looked for. With both reports asked for the line names the text
/// report, the first the run would create.
#[test]
fn a_report_named_inside_the_output_directory_is_refused_before_the_input_is_opened() {
    let scratch = Scratch::new();
    let input = vcf(&scratch.work(), "in.vcf", &[("21", 10), ("22", 5)]);
    let unsorted = vcf(&scratch.work(), "unsorted.vcf", &[("21", 100), ("21", 50)]);
    let out = scratch.at("out.parquet");
    let report = scratch.at("out.parquet/summary.html");
    let line = inside("statistics report", &report, &out, "--stats_file");

    let before = scratch.contents();
    let new_path = scratch.run(&input, &out, &["--stats_file", &report]);
    assert_refused_leaving(&scratch, &new_path, &line, &before);

    scratch.earlier_output(&input, "out.parquet");
    let before = scratch.contents();
    assert_eq!(
        scratch.names(),
        [
            "in.vcf",
            "out.parquet",
            "out.parquet/_vep_rs_manifest.json",
            "out.parquet/part-00000.parquet",
            "out.parquet/part-00001.parquet",
            "unsorted.vcf"
        ]
    );
    let unforced = scratch.run(&input, &out, &["--stats_file", &report]);
    assert_refused_leaving(&scratch, &unforced, &line, &before);
    let forced = scratch.run(
        &input,
        &out,
        &["--force_overwrite", "--stats_file", &report],
    );
    assert_refused_leaving(&scratch, &forced, &line, &before);
    let out_of_order = scratch.run(
        &unsorted,
        &out,
        &["--force_overwrite", "--stats_file", &report],
    );
    assert_refused_leaving(&scratch, &out_of_order, &line, &before);
    let absent = scratch.work().join("absent.vcf");
    let unopened = scratch.run(
        &absent,
        &out,
        &["--force_overwrite", "--stats_file", &report],
    );
    assert_refused_leaving(&scratch, &unopened, &line, &before);

    let root = scratch.at("out.parquet/summary");
    let both = scratch.run(
        &input,
        &out,
        &[
            "--force_overwrite",
            "--stats_text",
            "--stats_html",
            "--stats_file",
            &root,
        ],
    );
    let text_report = format!("{root}.txt");
    assert_refused_leaving(
        &scratch,
        &both,
        &inside("statistics report", &text_report, &out, "--stats_file"),
        &before,
    );
}

/// The warnings file of a run over the input [`skips`] writes, whole.
const SKIP_WARNING: &str = "WARNING: line 2 skipped (Z 5 r2 A G . . .): Chromosome Z not found in \
     annotation sources or synonyms; chromosome Z does not overlap any features\n";

/// The names of a completed output of one part at `out.parquet`.
const ONE_PART: [&str; 3] = [
    "out.parquet",
    "out.parquet/_vep_rs_manifest.json",
    "out.parquet/part-00000.parquet",
];

/// A VCF of one record the cache annotates and one on a chromosome it lacks,
/// which draws one warning.
fn skips(dir: &Path) -> PathBuf {
    vcf(dir, "skips.vcf", &[("21", 10), ("Z", 5)])
}

/// The names under the work directory but the inputs', `in.vcf` and
/// `skips.vcf`, sorted.
fn written(scratch: &Scratch) -> Vec<String> {
    let mut names = scratch.names();
    names.retain(|name| !matches!(name.as_str(), "in.vcf" | "skips.vcf"));
    names
}

/// `output` is a completed run's: exit 0 and no `ERROR:` line.
#[track_caller]
fn assert_completed(output: &Output) {
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(error_lines(output), [""; 0]);
}

/// An output path written with a trailing separator, with two, or with a `.`
/// step after one names its default report and warnings file as the same
/// path written without them does: the run completes, the directory holds
/// the Parquet output alone and both files stand beside it, the warnings
/// file byte for byte the one the plain path writes. A relative path does
/// the same from the working directory.
#[test]
fn an_output_path_with_a_trailing_separator_names_its_default_files_beside_the_directory() {
    let beside = [
        "out.parquet",
        "out.parquet/_vep_rs_manifest.json",
        "out.parquet/part-00000.parquet",
        "out.parquet_summary.html",
        "out.parquet_warnings.txt",
    ];
    let completed = |output_as_spelled: &dyn Fn(&Scratch) -> String| {
        let scratch = Scratch::new();
        let input = skips(&scratch.work());
        let output = scratch.run(&input, &output_as_spelled(&scratch), &[]);
        assert_completed(&output);
        assert_eq!(written(&scratch), beside, "{}", output_as_spelled(&scratch));
        let read = |name: &str| fs::read_to_string(scratch.work().join(name)).unwrap();
        assert!(
            read("out.parquet_summary.html").starts_with("<html>\n<head>\n  <title>VEP summary"),
            "the report is the run's"
        );
        read("out.parquet_warnings.txt")
    };

    let plain = completed(&|scratch| scratch.at("out.parquet"));
    assert_eq!(plain, SKIP_WARNING);
    for ending in ["/", "//", "/.", "/./"] {
        let ended = completed(&|scratch| format!("{}{ending}", scratch.at("out.parquet")));
        assert_eq!(ended, plain, "out.parquet{ending}");
    }
    let relative = completed(&|_| "out.parquet/".to_string());
    assert_eq!(relative, plain);
}

/// A Parquet directory named `0`, a name Perl reads as false, names its
/// warnings file `vep_warnings.txt` (BaseVEP.pm:890) and its report
/// `0_summary.html`, with a separator after the name as without one.
#[test]
fn a_directory_named_zero_names_its_default_files_with_a_separator_as_without_one() {
    for spelled in ["0", "0/", "0/."] {
        let scratch = Scratch::new();
        let input = skips(&scratch.work());
        let output = scratch.run(&input, spelled, &[]);
        assert_completed(&output);
        assert_eq!(
            written(&scratch),
            [
                "0",
                "0/_vep_rs_manifest.json",
                "0/part-00000.parquet",
                "0_summary.html",
                "vep_warnings.txt"
            ],
            "{spelled}"
        );
        assert_eq!(
            fs::read_to_string(scratch.work().join("vep_warnings.txt")).unwrap(),
            SKIP_WARNING,
            "{spelled}"
        );
    }
}

/// A text output keeps Ensembl VEP's default names on the output path as
/// spelled (BaseRunner.pm:254): written with a trailing separator it names
/// its report inside a directory that does not exist, so the run is refused
/// on the report and creates nothing.
#[test]
fn a_text_output_path_with_a_trailing_separator_names_its_report_as_spelled() {
    let scratch = Scratch::new();
    let input = skips(&scratch.work());
    let out = format!("{}/", scratch.at("out.txt"));
    let before = scratch.contents();
    // The last `--output_format` value replaces the `parquet` the helper gives.
    let as_text = scratch.run(&input, &out, &["--output_format", "vep"]);
    assert_refused_leaving(
        &scratch,
        &as_text,
        &format!("ERROR: [output] Could not write to stats file {out}_summary.html"),
        &before,
    );
}

/// Each default file of an output path written with a trailing separator
/// follows its own flag: no warnings file without a warning, no report under
/// `--no_stats`, the text report under `--stats_text`, both under both
/// flags, and a file a flag names elsewhere stands where the flag names it
/// while the other keeps its default name beside the directory. A
/// skipped-variants file named beside the directory leaves both default
/// names as they are.
#[test]
fn each_default_file_of_a_path_with_a_trailing_separator_follows_its_own_flag() {
    // The files a run over `records` with `extra` writes beside the output of
    // one part, and its standard error.
    let run = |records: &[(&str, u64)], extra: &[&str]| {
        let scratch = Scratch::new();
        let input = vcf(&scratch.work(), "in.vcf", records);
        let out = format!("{}/", scratch.at("out.parquet"));
        let extra: Vec<String> = extra
            .iter()
            .map(|word| word.replace("<work>", &scratch.work().display().to_string()))
            .collect();
        let extra: Vec<&str> = extra.iter().map(String::as_str).collect();
        let output = scratch.run(&input, &out, &extra);
        assert_completed(&output);
        let mut names = written(&scratch);
        assert_eq!(names[..3], ONE_PART, "{extra:?}");
        names.drain(..3);
        for warnings in ["out.parquet_warnings.txt", "w.txt"] {
            if let Ok(warnings) = fs::read_to_string(scratch.work().join(warnings)) {
                assert_eq!(warnings, SKIP_WARNING, "{extra:?}");
            }
        }
        (names, text(&output.stderr))
    };
    let annotated: &[(&str, u64)] = &[("21", 10)];
    let skipping: &[(&str, u64)] = &[("21", 10), ("Z", 5)];
    let none: [&str; 0] = [];

    let (names, _) = run(annotated, &[]);
    assert_eq!(names, ["out.parquet_summary.html"]);
    // This command line names no file of its own and draws no warning.
    let (names, stderr) = run(annotated, &["--no_stats"]);
    assert_eq!(names, none);
    assert_eq!(stderr, "");
    let (names, _) = run(skipping, &["--no_stats"]);
    assert_eq!(names, ["out.parquet_warnings.txt"]);
    let (names, _) = run(skipping, &["--stats_text"]);
    assert_eq!(
        names,
        ["out.parquet_summary.txt", "out.parquet_warnings.txt"]
    );
    let (names, _) = run(skipping, &["--stats_text", "--stats_html"]);
    assert_eq!(
        names,
        [
            "out.parquet_summary.html",
            "out.parquet_summary.txt",
            "out.parquet_warnings.txt"
        ]
    );
    let (names, stderr) = run(skipping, &["--warning_file", "STDERR"]);
    assert_eq!(names, ["out.parquet_summary.html"]);
    assert_eq!(stderr, SKIP_WARNING);
    let (names, stderr) = run(skipping, &["--no_stats", "--warning_file", "STDERR"]);
    assert_eq!(names, none);
    assert_eq!(stderr, SKIP_WARNING);
    let (names, _) = run(skipping, &["--stats_file", "<work>/report.html"]);
    assert_eq!(names, ["out.parquet_warnings.txt", "report.html"]);
    let (names, _) = run(skipping, &["--warning_file", "<work>/w.txt"]);
    assert_eq!(names, ["out.parquet_summary.html", "w.txt"]);
    let skipped = [
        "--skipped_variants_file",
        "<work>/skipped.txt",
        "--stats_text",
    ];
    let (names, _) = run(skipping, &skipped);
    assert_eq!(
        names,
        [
            "out.parquet_summary.txt",
            "out.parquet_warnings.txt",
            "skipped.txt"
        ]
    );
}

/// Over an earlier output a path written with a trailing separator is
/// refused as the plain path is, on the path as spelled and with nothing
/// changed or created beside the directory; under `--force_overwrite` the
/// run replaces the earlier output of two parts with its one part and writes
/// both default files beside it, and a further replacement under
/// `--no_stats` leaves the report as it is. A run that only shows the cache
/// writes nothing.
#[test]
fn a_path_with_a_trailing_separator_replaces_an_earlier_output_as_the_plain_path_does() {
    let scratch = Scratch::new();
    let input = vcf(&scratch.work(), "in.vcf", &[("21", 10), ("22", 5)]);
    let skipping = skips(&scratch.work());
    scratch.earlier_output(&input, "out.parquet");
    let out = format!("{}/", scratch.at("out.parquet"));
    let before = scratch.contents();

    let unforced = scratch.run(&skipping, &out, &[]);
    assert_refused_leaving(
        &scratch,
        &unforced,
        &format!(
            "ERROR: [output] Output file {out} already exists. Specify a different output file \
             with --output_file or overwrite existing file with --force_overwrite"
        ),
        &before,
    );

    let shown = scratch.run(&skipping, &out, &["--show_cache_info"]);
    assert_completed(&shown);
    assert_eq!(
        scratch.contents(),
        before,
        "a run that shows the cache wrote a file"
    );

    let forced = scratch.run(&skipping, &out, &["--force_overwrite"]);
    assert_completed(&forced);
    assert_eq!(
        written(&scratch),
        [
            "out.parquet",
            "out.parquet/_vep_rs_manifest.json",
            "out.parquet/part-00000.parquet",
            "out.parquet_summary.html",
            "out.parquet_warnings.txt"
        ]
    );
    assert_eq!(
        fs::read_to_string(scratch.work().join("out.parquet_warnings.txt")).unwrap(),
        SKIP_WARNING
    );

    // The manifest carries the identifier of the run that wrote it.
    let read = |name: &str| fs::read(scratch.work().join(name)).unwrap();
    let names = written(&scratch);
    let report = read("out.parquet_summary.html");
    let manifest = read("out.parquet/_vep_rs_manifest.json");
    let unreported = scratch.run(&skipping, &out, &["--force_overwrite", "--no_stats"]);
    assert_completed(&unreported);
    assert_eq!(written(&scratch), names);
    assert_eq!(text(&read("out.parquet_warnings.txt")), SKIP_WARNING);
    assert_eq!(read("out.parquet_summary.html"), report);
    assert_ne!(read("out.parquet/_vep_rs_manifest.json"), manifest);
}

/// A file that a flag names inside a directory written with a trailing
/// separator stays refused by name, the other file's default name beside
/// the directory notwithstanding: the report, then the warnings file, each
/// under the name a default made from the path as spelled would have.
#[test]
fn a_file_a_flag_names_inside_a_path_with_a_trailing_separator_stays_refused_by_name() {
    let scratch = Scratch::new();
    let input = skips(&scratch.work());
    let out = format!("{}/", scratch.at("out.parquet"));
    let before = scratch.contents();

    let report = format!("{out}_summary.html");
    let reported = scratch.run(&input, &out, &["--stats_file", &report]);
    assert_refused_leaving(
        &scratch,
        &reported,
        &inside("statistics report", &report, &out, "--stats_file"),
        &before,
    );

    let warnings = format!("{out}_warnings.txt");
    let warned = scratch.run(&input, &out, &["--warning_file", &warnings]);
    assert_refused_leaving(
        &scratch,
        &warned,
        &inside("warnings file", &warnings, &out, "--warning_file"),
        &before,
    );
}

/// A default name that a symbolic link leads inside the output directory is
/// refused as a named file is, in the one line that names the file, the
/// directory and the flag: for a path written with a trailing separator and
/// for the plain path alike.
#[cfg(unix)]
#[test]
fn a_default_name_a_link_leads_inside_the_output_directory_is_refused_by_name() {
    use std::os::unix::fs::symlink;
    let scratch = Scratch::new();
    let input = skips(&scratch.work());
    symlink(
        "out.parquet/warnings.txt",
        scratch.work().join("out.parquet_warnings.txt"),
    )
    .unwrap();
    let warnings = scratch.at("out.parquet_warnings.txt");
    let before = scratch.contents();
    for ending in ["/", ""] {
        let out = format!("{}{ending}", scratch.at("out.parquet"));
        let linked = scratch.run(&input, &out, &["--no_stats"]);
        assert_refused_leaving(
            &scratch,
            &linked,
            &inside("warnings file", &warnings, &out, "--warning_file"),
            &before,
        );
    }
}

/// A regular file and a symbolic link written with a trailing separator are
/// no output directory, and the default names beside them change neither
/// refusal; the link is refused before any file of the run is created. The
/// file is met when the run tests its output, before it creates its report,
/// so nothing is created, with the report on as under `--no_stats` and with
/// `--force_overwrite` as without it, and the file keeps its bytes.
#[cfg(unix)]
#[test]
fn a_file_or_a_link_written_with_a_trailing_separator_is_no_output_directory() {
    use std::os::unix::fs::symlink;
    let scratch = Scratch::new();
    let input = vcf(&scratch.work(), "in.vcf", &[("21", 10), ("22", 5)]);
    let skipping = skips(&scratch.work());
    scratch.earlier_output(&input, "real.parquet");
    symlink("real.parquet", scratch.work().join("link")).unwrap();
    fs::write(scratch.work().join("kept"), "kept text\n").unwrap();
    let before = scratch.contents();

    let link = scratch.at("link");
    let linked = scratch.run(&skipping, &format!("{link}/"), &["--force_overwrite"]);
    assert_refused_leaving(
        &scratch,
        &linked,
        &format!("ERROR: [parquet] output path {link} is a symbolic link"),
        &before,
    );

    let kept = format!("{}/", scratch.at("kept"));
    let unwritable = format!("ERROR: [output] Could not write to output file {kept}");
    let unreported = scratch.run(&skipping, &kept, &["--no_stats"]);
    assert_refused_leaving(&scratch, &unreported, &unwritable, &before);
    let forced = scratch.run(&skipping, &kept, &["--force_overwrite", "--no_stats"]);
    assert_refused_leaving(&scratch, &forced, &unwritable, &before);

    let reported = scratch.run(&skipping, &kept, &[]);
    assert_refused_leaving(&scratch, &reported, &unwritable, &before);
}

/// A warnings file and a skipped-variants file named inside the output
/// directory are refused as the report is, over an input that would draw a
/// warning and a skip; of several files named there the line names the
/// warnings file, then the report.
#[test]
fn a_warnings_file_or_a_skipped_variants_file_named_inside_the_output_is_refused_by_name() {
    let scratch = Scratch::new();
    let input = vcf(&scratch.work(), "in.vcf", &[("21", 10), ("22", 5)]);
    let skips = vcf(&scratch.work(), "skips.vcf", &[("21", 10), ("Z", 5)]);
    scratch.earlier_output(&input, "out.parquet");
    let out = scratch.at("out.parquet");
    let report = scratch.at("out.parquet/summary.html");
    let warnings = scratch.at("out.parquet/warnings.txt");
    let skipped = scratch.at("out.parquet/skipped.txt");
    let before = scratch.contents();

    let warned = scratch.run(
        &skips,
        &out,
        &[
            "--force_overwrite",
            "--no_stats",
            "--warning_file",
            &warnings,
        ],
    );
    assert_refused_leaving(
        &scratch,
        &warned,
        &inside("warnings file", &warnings, &out, "--warning_file"),
        &before,
    );

    let skipping = scratch.run(
        &skips,
        &out,
        &[
            "--force_overwrite",
            "--no_stats",
            "--warning_file",
            "STDERR",
            "--skipped_variants_file",
            &skipped,
        ],
    );
    assert_refused_leaving(
        &scratch,
        &skipping,
        &inside(
            "skipped-variants file",
            &skipped,
            &out,
            "--skipped_variants_file",
        ),
        &before,
    );

    let every_file = [
        "--force_overwrite",
        "--stats_file",
        &report,
        "--skipped_variants_file",
        &skipped,
        "--warning_file",
    ];
    let all = scratch.run(
        &skips,
        &out,
        &[&every_file[..], &[warnings.as_str()][..]].concat(),
    );
    assert_refused_leaving(
        &scratch,
        &all,
        &inside("warnings file", &warnings, &out, "--warning_file"),
        &before,
    );
    let without_the_warnings_file =
        scratch.run(&skips, &out, &[&every_file[..], &["STDERR"][..]].concat());
    assert_refused_leaving(
        &scratch,
        &without_the_warnings_file,
        &inside("statistics report", &report, &out, "--stats_file"),
        &before,
    );
}

/// A warning raised before the run starts creates no warnings file inside the
/// output: the warnings file is refused before the command line is resolved,
/// where a flag draws a warning, and before the cache is resolved, where a
/// cache with no checksum manifest draws one.
#[test]
fn a_warning_before_the_run_starts_creates_no_warnings_file_inside_the_output() {
    let scratch = Scratch::new();
    let input = vcf(&scratch.work(), "in.vcf", &[("21", 10), ("22", 5)]);
    scratch.earlier_output(&input, "out.parquet");
    let unchecked = scratch.root.path().join("cache");
    copy_without(&cache("GRCh38"), &unchecked, "CHECKSUMS.sha256");
    let out = scratch.at("out.parquet");
    let warnings = scratch.at("out.parquet/warnings.txt");
    let named = [
        "--force_overwrite",
        "--no_stats",
        "--warning_file",
        &warnings,
    ];
    let line = inside("warnings file", &warnings, &out, "--warning_file");
    let before = scratch.contents();

    let flagged = scratch.run_over(
        &cache("GRCh37"),
        "GRCh37",
        &input,
        &out,
        &[&named[..], &["--gencode_primary"][..]].concat(),
    );
    assert_refused_leaving(&scratch, &flagged, &line, &before);

    let unmanifested = scratch.run_over(&unchecked, "GRCh38", &input, &out, &named);
    assert_refused_leaving(&scratch, &unmanifested, &line, &before);

    // Beside the directory the run over that cache completes and writes its
    // warning.
    let beside = scratch.at("out.parquet_warnings.txt");
    let named = ["--force_overwrite", "--no_stats", "--warning_file", &beside];
    let completed = scratch.run_over(&unchecked, "GRCh38", &input, &out, &named);
    assert_eq!(
        completed.status.code(),
        Some(0),
        "{}",
        text(&completed.stderr)
    );
    assert!(
        fs::read_to_string(&beside)
            .unwrap()
            .contains("carries no checksum manifest"),
        "the warning of the cache is in the warnings file"
    );
}

/// A program that calls the runner itself and gives the warning sink a file
/// inside the output directory is refused by the runner, with the line the
/// binary prints: for a file under a name of its own, and for one inside a
/// directory the output path writes with a trailing separator.
#[test]
fn a_program_that_gives_the_sink_a_file_inside_the_output_is_refused_by_the_runner() {
    use vep_cli::warnings::{self, Target};
    let scratch = Scratch::new();
    let input = vcf(&scratch.work(), "in.vcf", &[("21", 10), ("22", 5)]);
    let refusal_of = |output: &str, warnings_file: &str| {
        warnings::global().configure(Target::File(PathBuf::from(warnings_file)), true, false);
        let config = vep_cli::Config::from_argv([
            "vep",
            "--json_cache",
            cache("GRCh38").to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--no_fasta",
            "--quiet",
            "--no_stats",
            "--output_format",
            "parquet",
            "--input_file",
            input.to_str().unwrap(),
            "--output_file",
            output,
        ])
        .unwrap_or_else(|_| panic!("the command line resolves"));
        let refusal = vep_cli::runner::Runner::new(config)
            .unwrap()
            .run()
            .unwrap_err();
        refusal.to_string()
    };
    let before = scratch.contents();

    let out = scratch.at("out.parquet");
    let named = scratch.at("out.parquet/warnings.txt");
    assert_eq!(
        refusal_of(&out, &named),
        inside("warnings file", &named, &out, "--warning_file")
    );

    let out = format!("{out}/");
    let under_the_separator = format!("{out}_warnings.txt");
    assert_eq!(
        refusal_of(&out, &under_the_separator),
        inside(
            "warnings file",
            &under_the_separator,
            &out,
            "--warning_file"
        )
    );
    assert_eq!(scratch.contents(), before, "a refused run wrote a file");
}

/// A file is placed by where its name leads and not by how it is spelled: a
/// relative name under an output named another way, a name through a
/// symbolic link to the output directory and a link whose target lies inside
/// it are refused, each named as the command line spells it.
#[cfg(unix)]
#[test]
fn a_file_another_spelling_places_inside_the_output_directory_is_refused() {
    use std::os::unix::fs::symlink;
    let scratch = Scratch::new();
    let input = vcf(&scratch.work(), "in.vcf", &[("21", 10), ("22", 5)]);
    scratch.earlier_output(&input, "out.parquet");
    fs::create_dir(scratch.work().join("sub")).unwrap();
    symlink("out.parquet", scratch.work().join("results")).unwrap();
    symlink(
        "out.parquet/summary.html",
        scratch.work().join("report.html"),
    )
    .unwrap();
    let before = scratch.contents();

    let relative = scratch.run(
        &input,
        "out.parquet",
        &[
            "--force_overwrite",
            "--stats_file",
            "./sub/../out.parquet/summary.html",
        ],
    );
    assert_refused_leaving(
        &scratch,
        &relative,
        &inside(
            "statistics report",
            "./sub/../out.parquet/summary.html",
            "out.parquet",
            "--stats_file",
        ),
        &before,
    );

    let out = scratch.at("out.parquet");
    for spelling in ["results/summary.html", "report.html"] {
        let report = scratch.at(spelling);
        let linked = scratch.run(
            &input,
            &out,
            &["--force_overwrite", "--stats_file", &report],
        );
        assert_refused_leaving(
            &scratch,
            &linked,
            &inside("statistics report", &report, &out, "--stats_file"),
            &before,
        );
    }
}

/// A report named as the output path itself is refused by name, at a new path
/// and under `--force_overwrite`, and so is one named as a regular file that
/// the output path spells with a trailing separator.
#[test]
fn a_report_named_as_the_output_path_is_refused_by_name() {
    let scratch = Scratch::new();
    let input = vcf(&scratch.work(), "in.vcf", &[("21", 10), ("22", 5)]);
    let out = scratch.at("out.html");
    let line = at_the_output_path("statistics report", &out, "--stats_file");
    let before = scratch.contents();
    for force in [&[][..], &["--force_overwrite"][..]] {
        let refused = scratch.run(&input, &out, &[force, &["--stats_file", &out][..]].concat());
        assert_refused_leaving(&scratch, &refused, &line, &before);
    }

    // A regular file named with a trailing separator is no output directory,
    // and a report named as that file would be written over it.
    let kept = scratch.at("kept.html");
    fs::write(&kept, "kept text\n").unwrap();
    let before = scratch.contents();
    let over_a_file = scratch.run(
        &input,
        &format!("{kept}/"),
        &["--force_overwrite", "--stats_file", &kept],
    );
    assert_refused_leaving(
        &scratch,
        &over_a_file,
        &at_the_output_path("statistics report", &kept, "--stats_file"),
        &before,
    );
}

/// A run whose files stand beside the output directory replaces an earlier
/// output and writes each of them: under a name that opens with the
/// directory's, in a directory whose name opens with it, and under a name
/// that steps into the directory and out again.
#[test]
fn a_run_whose_files_stand_beside_the_output_directory_completes() {
    let scratch = Scratch::new();
    let input = vcf(&scratch.work(), "in.vcf", &[("21", 10), ("22", 5)]);
    let skips = vcf(&scratch.work(), "skips.vcf", &[("21", 10), ("Z", 5)]);
    scratch.earlier_output(&input, "out.parquet");
    fs::create_dir(scratch.work().join("out.parquet.files")).unwrap();
    let out = scratch.at("out.parquet");

    let completed = scratch.run(
        &skips,
        &out,
        &[
            "--force_overwrite",
            "--stats_text",
            "--stats_file",
            &scratch.at("out.parquet/../out.parquet_summary.txt"),
            "--warning_file",
            &scratch.at("out.parquet_warnings.txt"),
            "--skipped_variants_file",
            &scratch.at("out.parquet.files/skipped.txt"),
        ],
    );
    assert_eq!(
        completed.status.code(),
        Some(0),
        "{}",
        text(&completed.stderr)
    );
    assert_eq!(
        scratch.names(),
        [
            "in.vcf",
            "out.parquet",
            "out.parquet.files",
            "out.parquet.files/skipped.txt",
            "out.parquet/_vep_rs_manifest.json",
            "out.parquet/part-00000.parquet",
            "out.parquet_summary.txt",
            "out.parquet_warnings.txt",
            "skips.vcf"
        ]
    );
    let read = |name: &str| fs::read_to_string(scratch.work().join(name)).unwrap();
    assert!(
        read("out.parquet_summary.txt").starts_with("[VEP run statistics]\n"),
        "the report is the run's"
    );
    assert!(
        read("out.parquet_warnings.txt").starts_with("WARNING: line 2 skipped "),
        "the warnings file is the run's"
    );
    assert!(
        read("out.parquet.files/skipped.txt").starts_with("[VEP skipped the following variants "),
        "the skipped-variants file is the run's"
    );
}

/// A report the run does not write is not refused: under `--no_stats` a
/// `--stats_file` inside the output directory names no file, and the run
/// replaces the earlier output of two parts with its one part.
#[test]
fn a_report_named_inside_the_output_that_the_run_does_not_write_is_not_refused() {
    let scratch = Scratch::new();
    let input = vcf(&scratch.work(), "in.vcf", &[("21", 10), ("22", 5)]);
    let replacement = vcf(&scratch.work(), "replacement.vcf", &[("21", 10)]);
    scratch.earlier_output(&input, "out.parquet");
    let out = scratch.at("out.parquet");
    let report = scratch.at("out.parquet/summary.html");
    let unreported = scratch.run(
        &replacement,
        &out,
        &["--force_overwrite", "--no_stats", "--stats_file", &report],
    );
    assert_eq!(
        unreported.status.code(),
        Some(0),
        "{}",
        text(&unreported.stderr)
    );
    assert_eq!(
        scratch.names(),
        [
            "in.vcf",
            "out.parquet",
            "out.parquet/_vep_rs_manifest.json",
            "out.parquet/part-00000.parquet",
            "replacement.vcf"
        ]
    );
}

/// The Parquet page quotes both refusals as the binary prints them and says
/// when the run is refused, what another spelling changes and what a trailing
/// separator leaves as it is. A page that puts the refusal at the run's last
/// step, or that has a trailing separator place the default names inside the
/// directory or refuse the run, is refused.
#[test]
fn the_parquet_page_states_the_refusal_of_a_file_named_inside_the_output() {
    let page = fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/output/parquet.md"),
    )
    .unwrap();
    let quoted_inside = inside("statistics report", "<file>", "<path>", "--stats_file");
    let quoted_at = at_the_output_path("statistics report", "<file>", "--stats_file");
    let (_, quoted_at) = quoted_at.split_once("<file> ").unwrap();
    let stated = [
        "Name the warnings file, the skipped-variants file and the statistics report outside \
         the output directory. A run whose warnings file, skipped-variants file or report would \
         be written inside it is refused before it opens its input, with exit 1"
            .to_string(),
        format!("The line names the file and the flag that places it: `{quoted_inside}`."),
        format!("A file named as the output path is refused as `{quoted_at}`."),
        "The run reads each name as the file system resolves it, through a symbolic link and a \
         `..` step, so a file that only another spelling places in the directory is refused as \
         well."
            .to_string(),
        "An output path written with a trailing separator (`out.parquet/`), or with a `.` step \
         after it, names its default report and warnings file as the same path without them \
         does, beside the directory: `out.parquet_summary.html` and `out.parquet_warnings.txt`."
            .to_string(),
    ];
    for sentence in &stated {
        assert!(
            page.contains(sentence.as_str()),
            "docs/output/parquet.md lacks {sentence:?}"
        );
    }
    for replaced in [
        "refuses that replacement at its last step",
        "places their default names inside it",
        "places the default names",
        "is refused the same way",
        "write the path without one",
    ] {
        assert!(
            !page.contains(replaced),
            "docs/output/parquet.md states {replaced:?}"
        );
    }
}
