// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The files a run names for itself, through the `vep` binary in every output
//! format: the statistics report, the output, the warnings file and the
//! skipped-variants file. Before one is created or emptied and before the
//! first record is read, the run tests that it can write the report, the
//! output, a warnings file that `--warning_file` names and the
//! skipped-variants file, and refuses the first that fails by name, in the
//! order it would meet them, with exit 1 and everything in its directory as
//! it was. A warnings file is still created by the first warning alone, one
//! drawn before the files are tested included, and a run whose files can all
//! be written completes. The warnings file under its default name is not
//! tested: the first warning meets it, as Ensembl VEP's first warning does
//! (BaseVEP.pm:817 and 890-894), so a run that draws no warning completes
//! where that file could not be created, from a directory closed to writing
//! too. `vep recode` tests its `--warning_file` before it annotates and
//! names no warnings file by default. A warnings file outside a Parquet
//! output takes a warning raised before a refusal of a file named inside it.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// A golden corpus of release 116.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The output formats and the flags that ask for each.
fn formats() -> Vec<(&'static str, &'static [&'static str])> {
    let mut formats: Vec<(&str, &[&str])> = vec![
        ("default", &[]),
        ("tab", &["--tab"]),
        ("vcf", &["--vcf"]),
        ("json", &["--json"]),
    ];
    if cfg!(feature = "parquet") {
        formats.push(("parquet", &["--output_format", "parquet"]));
    }
    formats
}

/// A scratch directory and, inside it, the directory the runs of one test
/// read and write, each run started there so that its files are named by
/// relative paths; the binary's home stands beside that directory.
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

    /// A VCF of one made-up substitution per `(chromosome, position)`,
    /// written beside the work directory. The release cache holds no
    /// chromosome `Z`, so a record there is skipped with a warning.
    fn vcf(&self, name: &str, records: &[(&str, u64)]) -> PathBuf {
        let path = self.root.path().join(name);
        let mut text =
            String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (i, (chrom, pos)) in records.iter().enumerate() {
            text.push_str(&format!("{chrom}\t{pos}\tr{}\tA\tG\t.\t.\t.\n", i + 1));
        }
        fs::write(&path, text).unwrap();
        path
    }

    /// An input of three records, the second of them skipped with a warning.
    fn input_that_warns(&self) -> PathBuf {
        self.vcf("warns.vcf", &[("21", 100), ("Z", 5), ("21", 300)])
    }

    /// An input of two records that draws no warning.
    fn input_without_a_warning(&self) -> PathBuf {
        self.vcf("clean.vcf", &[("21", 100), ("21", 300)])
    }

    /// A run of `format` over `input`, started in the work directory, with
    /// `extra` after its own flags.
    fn run(&self, format: &[&str], input: &Path, extra: &[&str]) -> Output {
        common::isolated_vep(self.root.path())
            .current_dir(self.work())
            .arg("--json_cache")
            .arg(corpus("GRCh38").join("json_cache"))
            .args(["--assembly", "GRCh38", "--no_fasta"])
            .args(format)
            .arg("--input_file")
            .arg(input)
            .args(extra)
            .stdin(Stdio::null())
            .output()
            .unwrap()
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
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// The lines of a run's standard error that open with `prefix`.
fn lines_opening(output: &Output, prefix: &str) -> Vec<String> {
    text(&output.stderr)
        .lines()
        .filter(|line| line.starts_with(prefix))
        .map(str::to_string)
        .collect()
}

/// `output` is refused with exit 1 and `line` as its one `ERROR:` line, it
/// echoes no warning, and the work directory holds what it held in `before`,
/// byte for byte.
#[track_caller]
fn assert_refused_leaving(
    scratch: &Scratch,
    output: &Output,
    line: &str,
    before: &[(String, Option<Vec<u8>>)],
    run: &str,
) {
    let stderr = text(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{run}: {stderr}");
    assert_eq!(lines_opening(output, "ERROR:"), [line], "{run}: {stderr}");
    assert_eq!(
        lines_opening(output, "WARNING:"),
        Vec::<String>::new(),
        "{run}: a record was read before the refusal"
    );
    assert_eq!(
        scratch.contents(),
        before,
        "{run}: the refused run created, changed or removed a file"
    );
}

/// A warnings file that `--warning_file` names and that cannot be written is
/// refused by name before the first record is read, whether or not the input
/// would draw a warning: one named in a directory that does not exist, and
/// one named where a directory stands. No output is written.
#[test]
fn a_warnings_file_that_cannot_be_written_is_refused_before_the_first_record() {
    for (format, flags) in formats() {
        let scratch = Scratch::new();
        fs::create_dir(scratch.work().join("w.dir")).unwrap();
        let inputs = [
            ("a warning", scratch.input_that_warns()),
            ("no warning", scratch.input_without_a_warning()),
        ];
        for (drawn, input) in &inputs {
            for name in ["absent/w.txt", "w.dir"] {
                let before = scratch.contents();
                let output = scratch.run(
                    flags,
                    input,
                    &["-o", "out", "--no_stats", "--warning_file", name],
                );
                assert_refused_leaving(
                    &scratch,
                    &output,
                    &format!("ERROR: [output] Could not write to warnings file {name}"),
                    &before,
                    &format!("{format}, {drawn}, {name}"),
                );
            }
        }
    }
}

/// The data lines of an output: every line that is not a header line.
fn data_lines(output: &[u8]) -> Vec<String> {
    text(output)
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The warnings file under its default name, beside the output, is not
/// tested before the first record: the first warning creates it, as Ensembl
/// VEP's first warning does (BaseVEP.pm:817 calls `warning_fh`, which names
/// the file at 890 and opens it at 894). Where a directory stands at that
/// name a run that draws no warning completes with its output written, and a
/// run that draws a warning ends once the batch that warned is written, on
/// the line of a warnings file that cannot be written, the warning echoed
/// before it. A Parquet run that ends so leaves nothing at its output path.
#[test]
fn the_default_warnings_file_is_first_met_by_the_first_warning() {
    let skipped = "WARNING: line 2 skipped (Z 5 r2 A G . . .): Chromosome Z not found in \
                   annotation sources or synonyms; chromosome Z does not overlap any features";
    for (format, flags) in formats() {
        let scratch = Scratch::new();
        let work = scratch.work();
        for output in ["clean", "out"] {
            fs::create_dir(work.join(format!("{output}_warnings.txt"))).unwrap();
        }
        let clean = scratch.run(
            flags,
            &scratch.input_without_a_warning(),
            &["-o", "clean", "--no_stats"],
        );
        let stderr = text(&clean.stderr);
        assert_eq!(clean.status.code(), Some(0), "{format}: {stderr}");
        assert_eq!(
            lines_opening(&clean, "ERROR:"),
            Vec::<String>::new(),
            "{format}: {stderr}"
        );
        assert!(work.join("clean").exists(), "{format}: no output");

        // The same input with every warning on standard error gives the rows
        // the batch that warned holds.
        let plain = scratch.run(
            flags,
            &scratch.input_that_warns(),
            &["-o", "plain", "--no_stats", "--warning_file", "STDERR"],
        );
        assert_eq!(plain.status.code(), Some(0), "{format}");
        let warned = scratch.run(
            flags,
            &scratch.input_that_warns(),
            &["-o", "out", "--no_stats"],
        );
        let stderr = text(&warned.stderr);
        assert_eq!(warned.status.code(), Some(1), "{format}: {stderr}");
        assert_eq!(
            lines_opening(&warned, "ERROR:"),
            ["ERROR: [output] Could not write to warnings file out_warnings.txt"],
            "{format}: {stderr}"
        );
        assert_eq!(
            lines_opening(&warned, "WARNING:"),
            [skipped],
            "{format}: {stderr}"
        );
        assert!(
            stderr.find(skipped) < stderr.find("ERROR:"),
            "{format}: the warning is not echoed before the error line: {stderr}"
        );
        if format == "parquet" {
            assert!(
                !work.join("out").exists(),
                "{format}: the failed run left something at its output path"
            );
        } else {
            let rows = data_lines(&fs::read(work.join("out")).unwrap());
            assert_eq!(rows.len(), 2, "{format}: {rows:?}");
            assert_eq!(
                rows,
                data_lines(&fs::read(work.join("plain")).unwrap()),
                "{format}: the batch that warned is not written whole"
            );
        }
        for output in ["clean", "out"] {
            let stands = work.join(format!("{output}_warnings.txt"));
            assert_eq!(
                fs::read_dir(&stands).unwrap().count(),
                0,
                "{format}: something was written where the directory {} stands",
                stands.display()
            );
        }
    }
}

/// From a directory closed to writing, a run that writes to standard output,
/// asks for no report and draws no warning completes: nothing tests the
/// warnings file its default name would put there (`STDOUT_warnings.txt`),
/// which the first warning alone would create. The rows are those of the
/// same run from a directory it can write, and the directory stays empty.
/// A `--warning_file` value that is false in Perl, `0` or empty, names no
/// file (BaseVEP.pm:890), so that run completes as well. `vep recode`, which
/// names no warnings file by default, completes there with and without a
/// warning.
#[cfg(unix)]
#[test]
fn a_run_that_draws_no_warning_completes_from_a_directory_closed_to_writing() {
    use std::os::unix::fs::PermissionsExt;
    let scratch = Scratch::new();
    let work = scratch.work();
    let input = scratch.input_without_a_warning();
    let files = ["-o", "STDOUT", "--no_stats"];
    let written: Vec<_> = formats()
        .into_iter()
        .filter(|(format, _)| *format != "parquet")
        .map(|(format, flags)| (format, flags, scratch.run(flags, &input, &files)))
        .collect();
    assert!(scratch.contents().is_empty(), "a run left a file behind");

    fs::set_permissions(&work, fs::Permissions::from_mode(0o555)).unwrap();
    for (format, flags, open) in &written {
        assert_eq!(open.status.code(), Some(0), "{format}");
        let closed = scratch.run(flags, &input, &files);
        let stderr = text(&closed.stderr);
        assert_eq!(closed.status.code(), Some(0), "{format}: {stderr}");
        assert_eq!(
            lines_opening(&closed, "ERROR:"),
            Vec::<String>::new(),
            "{format}: {stderr}"
        );
        let rows = data_lines(&closed.stdout);
        assert_eq!(rows.len(), 2, "{format}: {rows:?}");
        assert_eq!(rows, data_lines(&open.stdout), "{format}");
        if *format == "default" {
            assert_eq!(
                rows,
                [
                    "r1\t21:100\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
                    "r2\t21:300\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
                ]
            );
        }
    }

    let (_, _, default_format) = &written[0];
    for false_value in ["0", ""] {
        let closed = scratch.run(
            &[],
            &input,
            &["-o", "STDOUT", "--no_stats", "--warning_file", false_value],
        );
        let stderr = text(&closed.stderr);
        assert_eq!(closed.status.code(), Some(0), "{false_value:?}: {stderr}");
        assert_eq!(
            data_lines(&closed.stdout),
            data_lines(&default_format.stdout),
            "{false_value:?}"
        );
    }

    let recoder = corpus("GRCh38-recoder");
    for (id, array) in [
        (
            "rs781261918",
            "[{\"T\":{\"hgvsg\":[\"21:g.25702257C>T\"],\"input\":\"rs781261918\"}}]",
        ),
        (
            "rs_nope",
            "[{\"warnings\":[\"No variant found with ID 'rs_nope'\"]}]",
        ),
    ] {
        let recoded = common::isolated_vep(scratch.root.path())
            .current_dir(&work)
            .arg("recode")
            .arg("--json_cache")
            .arg(recoder.join("json_cache"))
            .arg("--fasta")
            .arg(recoder.join("reference.fa.gz"))
            .args(["--fields", "hgvsg", "--id", id])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        let stderr = text(&recoded.stderr);
        assert_eq!(recoded.status.code(), Some(0), "{id}: {stderr}");
        assert_eq!(text(&recoded.stdout), array, "{id}");
    }
    assert!(scratch.contents().is_empty(), "a run left a file behind");
    fs::set_permissions(&work, fs::Permissions::from_mode(0o755)).unwrap();
}

/// A warning drawn before the files are tested is written to the warnings
/// file at once: a run that names a `--fasta_dir` holding no FASTA and an
/// output it cannot write is refused for the output, with that warning in the
/// warnings file in place of the file's earlier bytes and no report created.
#[test]
fn a_warning_drawn_before_the_files_are_tested_is_written_to_the_warnings_file() {
    let scratch = Scratch::new();
    let work = scratch.work();
    fs::create_dir(work.join("fasta.dir")).unwrap();
    fs::write(work.join("w.txt"), "earlier warnings\n").unwrap();
    let output = scratch.run(
        &[],
        &scratch.input_without_a_warning(),
        &[
            "--fasta_dir",
            "fasta.dir",
            "-o",
            "absent/out",
            "--stats_file",
            "s.html",
            "--warning_file",
            "w.txt",
        ],
    );
    let stderr = text(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert_eq!(
        lines_opening(&output, "ERROR:"),
        ["ERROR: [output] Could not write to output file absent/out"],
        "{stderr}"
    );
    assert_eq!(
        text(&fs::read(work.join("w.txt")).unwrap()),
        "WARNING: No compatible FASTA files found in fasta.dir\n"
    );
    assert!(
        !work.join("s.html").exists(),
        "the refused run created its report"
    );
}

/// A skipped-variants file that cannot be written is refused by name before
/// the first record is read, though the run writes the file only once every
/// record is annotated: one named in a directory that does not exist, one
/// named where a directory stands, under `--force_overwrite`, and the empty
/// name. Neither the output nor the statistics report is written.
#[test]
fn a_skipped_variants_file_that_cannot_be_written_is_refused_before_the_first_record() {
    for (format, flags) in formats() {
        let scratch = Scratch::new();
        fs::create_dir(scratch.work().join("skipped.dir")).unwrap();
        let input = scratch.input_that_warns();
        for name in ["absent/skipped.txt", "skipped.dir", ""] {
            let before = scratch.contents();
            let output = scratch.run(
                flags,
                &input,
                &[
                    "-o",
                    "out",
                    "--stats_file",
                    "s.html",
                    "--force_overwrite",
                    "--skipped_variants_file",
                    name,
                ],
            );
            assert_refused_leaving(
                &scratch,
                &output,
                &format!("ERROR: [output] Could not write to skipped variants file {name}"),
                &before,
                &format!("{format}, {name:?}"),
            );
        }
    }
}

/// A statistics report is neither emptied nor created for a run whose output
/// cannot be written: an earlier report keeps its bytes under
/// `--force_overwrite`, and a run with no earlier report leaves none. With
/// both reports asked for, the text report keeps its bytes when the HTML
/// report cannot be written.
#[test]
fn a_report_is_neither_emptied_nor_created_for_a_run_that_is_refused() {
    for (format, flags) in formats() {
        let scratch = Scratch::new();
        let input = scratch.input_without_a_warning();
        let unwritable_output = "ERROR: [output] Could not write to output file absent/out";

        let before = scratch.contents();
        let output = scratch.run(
            flags,
            &input,
            &["-o", "absent/out", "--stats_file", "s.html"],
        );
        assert_refused_leaving(
            &scratch,
            &output,
            unwritable_output,
            &before,
            &format!("{format}, no earlier report"),
        );

        fs::write(scratch.work().join("s.html"), "an earlier report").unwrap();
        let before = scratch.contents();
        let output = scratch.run(
            flags,
            &input,
            &[
                "-o",
                "absent/out",
                "--stats_file",
                "s.html",
                "--force_overwrite",
            ],
        );
        assert_refused_leaving(
            &scratch,
            &output,
            unwritable_output,
            &before,
            &format!("{format}, an earlier report"),
        );

        fs::write(scratch.work().join("both.txt"), "an earlier text report").unwrap();
        fs::create_dir(scratch.work().join("both.html")).unwrap();
        let before = scratch.contents();
        let output = scratch.run(
            flags,
            &input,
            &[
                "-o",
                "out",
                "--stats_text",
                "--stats_html",
                "--stats_file",
                "both",
                "--force_overwrite",
            ],
        );
        assert_refused_leaving(
            &scratch,
            &output,
            "ERROR: [output] Could not write to stats file both.html",
            &before,
            &format!("{format}, both reports"),
        );
    }
}

/// An output named under a regular file, by a trailing separator, and an
/// output path with no final name under `--force_overwrite` cannot be
/// written, and the run is refused for the output with no report created and
/// an earlier report at its bytes: the output is tested before a report is
/// created or emptied.
#[test]
fn an_output_under_a_file_or_with_no_final_name_is_refused_with_no_report_created() {
    for (format, flags) in formats() {
        for earlier_report in [false, true] {
            let scratch = Scratch::new();
            let input = scratch.input_without_a_warning();
            fs::write(scratch.work().join("kept"), "kept").unwrap();
            if earlier_report {
                fs::write(scratch.work().join("s.html"), "an earlier report").unwrap();
            }
            let before = scratch.contents();
            for output in ["kept/", "."] {
                let refused = scratch.run(
                    flags,
                    &input,
                    &[
                        "-o",
                        output,
                        "--force_overwrite",
                        "--stats_file",
                        "s.html",
                        "--warning_file",
                        "w.txt",
                    ],
                );
                assert_refused_leaving(
                    &scratch,
                    &refused,
                    &format!("ERROR: [output] Could not write to output file {output}"),
                    &before,
                    &format!("{format}, {output}, an earlier report: {earlier_report}"),
                );
            }
        }
    }
}

/// Under `--force_overwrite` an earlier output, report and warnings file keep
/// their bytes when the run's skipped-variants file cannot be written. The
/// earlier output of a Parquet run is a completed output; of every other
/// format, a file of bytes the run would replace.
#[test]
fn earlier_files_keep_their_bytes_when_a_later_file_cannot_be_written() {
    for (format, flags) in formats() {
        let scratch = Scratch::new();
        let input = scratch.input_that_warns();
        if format == "parquet" {
            let earlier = scratch.run(flags, &input, &["-o", "out", "--no_stats", "--quiet"]);
            assert_eq!(earlier.status.code(), Some(0), "{}", text(&earlier.stderr));
            fs::remove_file(scratch.work().join("out_warnings.txt")).unwrap();
        } else {
            fs::write(scratch.work().join("out"), "an earlier output").unwrap();
        }
        fs::write(scratch.work().join("s.html"), "an earlier report").unwrap();
        fs::write(scratch.work().join("w.txt"), "earlier warnings\n").unwrap();
        let before = scratch.contents();
        let output = scratch.run(
            flags,
            &input,
            &[
                "-o",
                "out",
                "--force_overwrite",
                "--stats_file",
                "s.html",
                "--warning_file",
                "w.txt",
                "--skipped_variants_file",
                "absent/skipped.txt",
            ],
        );
        assert_refused_leaving(
            &scratch,
            &output,
            "ERROR: [output] Could not write to skipped variants file absent/skipped.txt",
            &before,
            format,
        );
    }
}

/// Of several files that cannot be written the refusal names the first in the
/// order the run would meet them: the statistics report, the output, the
/// warnings file, the skipped-variants file. No file that can be written is
/// created for a run that names a later one that cannot.
#[test]
fn the_refusal_names_the_first_file_in_the_order_the_run_meets_them() {
    for (format, flags) in formats() {
        let scratch = Scratch::new();
        let input = scratch.input_that_warns();
        let before = scratch.contents();
        for (report, out, warnings, line) in [
            (
                "absent/s.html",
                "absent/out",
                "absent/w.txt",
                "ERROR: [output] Could not write to stats file absent/s.html",
            ),
            (
                "s.html",
                "absent/out",
                "absent/w.txt",
                "ERROR: [output] Could not write to output file absent/out",
            ),
            (
                "s.html",
                "out",
                "absent/w.txt",
                "ERROR: [output] Could not write to warnings file absent/w.txt",
            ),
            (
                "s.html",
                "out",
                "w.txt",
                "ERROR: [output] Could not write to skipped variants file absent/skipped.txt",
            ),
        ] {
            let output = scratch.run(
                flags,
                &input,
                &[
                    "-o",
                    out,
                    "--stats_file",
                    report,
                    "--warning_file",
                    warnings,
                    "--skipped_variants_file",
                    "absent/skipped.txt",
                ],
            );
            assert_refused_leaving(
                &scratch,
                &output,
                line,
                &before,
                &format!("{format}, {report}, {out}, {warnings}"),
            );
        }
    }
}

/// A run whose files can all be written completes, and the test of each
/// leaves nothing behind: a run that draws no warning writes its output, its
/// report and its skipped-variants file and no warnings file, and leaves an
/// earlier warnings file its bytes. A run that draws a warning writes the
/// warnings file.
#[test]
fn a_run_whose_files_can_be_written_completes_and_a_warning_alone_creates_the_warnings_file() {
    for (format, flags) in formats() {
        let scratch = Scratch::new();
        let work = scratch.work();
        let files = [
            "-o",
            "out",
            "--stats_file",
            "s.html",
            "--warning_file",
            "w.txt",
            "--skipped_variants_file",
            "skipped.txt",
            "--force_overwrite",
        ];
        let clean = scratch.run(flags, &scratch.input_without_a_warning(), &files);
        assert_eq!(
            clean.status.code(),
            Some(0),
            "{format}: {}",
            text(&clean.stderr)
        );
        assert!(work.join("out").exists(), "{format}: no output");
        assert!(
            text(&fs::read(work.join("s.html")).unwrap()).contains("</html>"),
            "{format}: the report is not written"
        );
        assert!(
            text(&fs::read(work.join("skipped.txt")).unwrap())
                .starts_with("[VEP skipped the following variants from "),
            "{format}"
        );
        assert!(
            !work.join("w.txt").exists(),
            "{format}: a run that drew no warning left a warnings file"
        );

        fs::write(work.join("w.txt"), "earlier warnings\n").unwrap();
        let clean = scratch.run(flags, &scratch.input_without_a_warning(), &files);
        assert_eq!(
            clean.status.code(),
            Some(0),
            "{format}: {}",
            text(&clean.stderr)
        );
        assert_eq!(
            text(&fs::read(work.join("w.txt")).unwrap()),
            "earlier warnings\n",
            "{format}: the test of the warnings file changed it"
        );

        let warned = scratch.run(flags, &scratch.input_that_warns(), &files);
        assert_eq!(
            warned.status.code(),
            Some(0),
            "{format}: {}",
            text(&warned.stderr)
        );
        assert_eq!(
            text(&fs::read(work.join("w.txt")).unwrap()),
            "WARNING: line 2 skipped (Z 5 r2 A G . . .): Chromosome Z not found in annotation \
             sources or synonyms; chromosome Z does not overlap any features\n",
            "{format}"
        );
    }
}

/// A warnings file named through a symbolic link that leads nowhere is tested
/// at the name the link gives: a run that draws no warning leaves the link as
/// it was and creates nothing, and a warning is written through it.
#[cfg(unix)]
#[test]
fn a_warnings_file_named_through_a_link_that_leads_nowhere_is_tested_where_it_leads() {
    let scratch = Scratch::new();
    let work = scratch.work();
    std::os::unix::fs::symlink("w.target", work.join("w.link")).unwrap();
    let files = ["-o", "out", "--no_stats", "--warning_file", "w.link"];
    let clean = scratch.run(&[], &scratch.input_without_a_warning(), &files);
    assert_eq!(clean.status.code(), Some(0), "{}", text(&clean.stderr));
    assert!(
        !work.join("w.target").exists(),
        "the test left a file behind"
    );
    assert_eq!(
        fs::read_link(work.join("w.link")).unwrap(),
        Path::new("w.target")
    );

    let files = [&files[..], &["--force_overwrite"]].concat();
    let warned = scratch.run(&[], &scratch.input_that_warns(), &files);
    assert_eq!(warned.status.code(), Some(0), "{}", text(&warned.stderr));
    assert!(text(&fs::read(work.join("w.target")).unwrap()).starts_with("WARNING: line 2 skipped"));

    std::os::unix::fs::symlink("absent/w.target", work.join("lost.link")).unwrap();
    let before = scratch.contents();
    let refused = scratch.run(
        &[],
        &scratch.input_that_warns(),
        &[
            "-o",
            "out",
            "--no_stats",
            "--force_overwrite",
            "--warning_file",
            "lost.link",
        ],
    );
    assert_refused_leaving(
        &scratch,
        &refused,
        "ERROR: [output] Could not write to warnings file lost.link",
        &before,
        "a link into a directory that does not exist",
    );
}

/// `vep recode` refuses a `--warning_file` it cannot write before it
/// annotates, whether or not the input draws a warning, though it writes the
/// file only once the annotation has ended and a warning stands to be
/// written. The file is tested where the annotation run tests its own files,
/// so a FASTA that is not there is refused first. A file it can write is
/// created by a warning alone, and an earlier file keeps its bytes through a
/// run that draws none. A device is first met by that write: Linux's
/// `/dev/full` ends a run that draws a warning on the line of a warnings file
/// that cannot be written, with no output, and a run that draws none
/// completes.
#[test]
fn vep_recode_refuses_a_warnings_file_it_cannot_write_before_it_annotates() {
    let scratch = Scratch::new();
    let work = scratch.work();
    let recoder = corpus("GRCh38-recoder");
    let recode_over = |fasta: &Path, id: &str, warnings: &str| {
        common::isolated_vep(scratch.root.path())
            .current_dir(&work)
            .arg("recode")
            .arg("--json_cache")
            .arg(recoder.join("json_cache"))
            .arg("--fasta")
            .arg(fasta)
            .args(["--fields", "hgvsg", "--id", id])
            .args(["--warning_file", warnings])
            .stdin(Stdio::null())
            .output()
            .unwrap()
    };
    let recode =
        |id: &str, warnings: &str| recode_over(&recoder.join("reference.fa.gz"), id, warnings);
    let before = scratch.contents();
    for id in ["rs_nope", "rs781261918"] {
        let refused = recode(id, "absent/w.txt");
        assert_eq!(refused.status.code(), Some(1), "{id}");
        assert_eq!(
            lines_opening(&refused, "ERROR:"),
            ["ERROR: [recode] Could not write to warnings file absent/w.txt"],
            "{id}: {}",
            text(&refused.stderr)
        );
        assert!(refused.stdout.is_empty(), "{id}: an output was written");
        assert_eq!(scratch.contents(), before, "{id}");
    }
    let no_fasta = recode_over(&work.join("absent.fa"), "rs_nope", "absent/w.txt");
    assert_eq!(no_fasta.status.code(), Some(1));
    let errors = lines_opening(&no_fasta, "ERROR:");
    assert!(
        matches!(errors.as_slice(), [line] if line.starts_with("ERROR: [fasta] ")),
        "{errors:?}"
    );

    let completed = recode("rs781261918", "w.txt");
    assert_eq!(
        completed.status.code(),
        Some(0),
        "{}",
        text(&completed.stderr)
    );
    assert!(
        !work.join("w.txt").exists(),
        "a run that drew no warning left a warnings file"
    );
    fs::write(work.join("w.txt"), "earlier warnings\n").unwrap();
    let completed = recode("rs781261918", "w.txt");
    assert_eq!(
        completed.status.code(),
        Some(0),
        "{}",
        text(&completed.stderr)
    );
    assert_eq!(
        text(&fs::read(work.join("w.txt")).unwrap()),
        "earlier warnings\n"
    );
    let warned = recode("rs_nope", "w.txt");
    assert_eq!(warned.status.code(), Some(0), "{}", text(&warned.stderr));
    assert_eq!(
        text(&fs::read(work.join("w.txt")).unwrap()),
        "WARNING: No variant found with ID 'rs_nope'\n"
    );

    if cfg!(target_os = "linux") {
        let full = recode("rs_nope", "/dev/full");
        assert_eq!(full.status.code(), Some(1), "{}", text(&full.stderr));
        assert_eq!(
            lines_opening(&full, "ERROR:"),
            ["ERROR: [recode] Could not write to warnings file /dev/full"]
        );
        assert!(full.stdout.is_empty(), "an output was written");
        let unmet = recode("rs781261918", "/dev/full");
        assert_eq!(unmet.status.code(), Some(0), "{}", text(&unmet.stderr));
    }
}

/// A warnings file outside a Parquet output takes a warning raised before the
/// run is refused for a file named inside that output: a cache with no
/// checksum manifest draws its warning while the cache is resolved, the
/// report named inside the output directory refuses the run before its input
/// is opened, and the warnings file holds the warning alone.
#[cfg(feature = "parquet")]
#[test]
fn a_warnings_file_outside_a_parquet_output_takes_a_warning_raised_before_the_refusal() {
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
    let scratch = Scratch::new();
    let work = scratch.work();
    let unchecked = scratch.root.path().join("cache");
    copy_without(
        &corpus("GRCh38").join("json_cache"),
        &unchecked,
        "CHECKSUMS.sha256",
    );
    let refused = common::isolated_vep(scratch.root.path())
        .current_dir(&work)
        .arg("--json_cache")
        .arg(&unchecked)
        .args(["--assembly", "GRCh38", "--no_fasta"])
        .args(["--output_format", "parquet"])
        .arg("--input_file")
        .arg(scratch.input_without_a_warning())
        .args(["-o", "out.parquet"])
        .args(["--stats_file", "out.parquet/s.html"])
        .args(["--warning_file", "w.txt"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let stderr = text(&refused.stderr);
    assert_eq!(refused.status.code(), Some(1), "{stderr}");
    assert_eq!(
        lines_opening(&refused, "ERROR:"),
        [
            "ERROR: [parquet] statistics report out.parquet/s.html lies inside the output \
          directory out.parquet; name it outside that directory with --stats_file"
        ],
        "{stderr}"
    );
    assert_eq!(
        text(&fs::read(work.join("w.txt")).unwrap()),
        format!(
            "WARNING: cache {} carries no checksum manifest\n",
            unchecked.display()
        )
    );
    assert!(
        !work.join("out.parquet").exists(),
        "the refused run left something at its output path"
    );
}

/// The Parquet page says when the report beside an output is emptied: once
/// the tests the run makes before its first record have passed, of the
/// report, the output, a file `--warning_file` names and the skipped-variants
/// file. It names no test of the default warnings file, which the run does
/// not make.
#[test]
fn the_parquet_page_says_when_the_report_beside_an_output_is_emptied() {
    let page = fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/output/parquet.md"),
    )
    .unwrap();
    let page = page.split_whitespace().collect::<Vec<_>>().join(" ");
    assert!(
        page.contains(
            "The statistics report beside it is a file of its own, emptied before the run reads \
             its first record, once the tests it makes of the report, the output, a \
             `--warning_file` and the skipped-variants file have passed."
        ),
        "docs/output/parquet.md does not say when the report is emptied"
    );
    for older in [
        "emptied as a replacement starts",
        "tested that it can write every file it names",
    ] {
        assert!(
            !page.contains(older),
            "docs/output/parquet.md states {older:?}"
        );
    }
}
