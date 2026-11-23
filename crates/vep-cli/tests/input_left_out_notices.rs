// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the `vep` binary says when it leaves part of what it was given to read
//! out of a run: a line of a list-based input that begins with `#`, read as a
//! comment; a FASTA file passed over for want of its index; under `--hgvs` or
//! `--hgvsg`, a chromosome the FASTA holds no sequence for; under `--hgvs`, an
//! insertion or a deletion on a transcript edited away from the genome; and a
//! VCF record whose ALT is a lone `.` under a true INFO/SVTYPE. A warning is
//! written once, to the warnings file and, unless `--quiet`, to standard error;
//! an `INFO:` notice once, to standard error alone, and not under `--quiet`.
//! Neither changes a row of the output. `vep recode`, whose annotation run is
//! quiet, writes the notices of that run itself, before its own line, in a run
//! that completes.

use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

/// The identifiers corpus: a chromosome 21 cache slice and a reference FASTA
/// over its records.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// One run of the binary over `input` with the corpus cache and FASTA.
struct Run {
    output: Output,
    /// The data rows of the output file.
    rows: Vec<String>,
    /// The warnings file, empty when the run wrote none.
    warnings: String,
}

impl Run {
    /// The lines of standard error that open with `WARNING` and hold `about`,
    /// the words of the warning a test reads; any other line is not its
    /// concern.
    fn warned_on_standard_error(&self, about: &str) -> Vec<String> {
        text(&self.output.stderr)
            .lines()
            .filter(|line| line.starts_with("WARNING") && line.contains(about))
            .map(str::to_string)
            .collect()
    }

    fn warning_lines(&self) -> Vec<&str> {
        self.warnings.lines().collect()
    }
}

/// Runs the binary on `input`, writing `<name>.txt` and its warnings file under
/// `dir`, and requires exit status 0.
fn annotate(dir: &Path, name: &str, input: &Path, extra: &[&str]) -> Run {
    let out = dir.join(format!("{name}.txt"));
    let mut command = common::isolated_vep(dir);
    command
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus().join("json_cache"))
        .arg("--fasta")
        .arg(corpus().join("reference.fa.gz"))
        .args([
            "--assembly",
            "GRCh38",
            "--offline",
            "--force_overwrite",
            "--no_stats",
        ])
        .args(extra);
    let output = command.output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let rows = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect();
    let warnings =
        std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap_or_default();
    Run {
        output,
        rows,
        warnings,
    }
}

/// A line that begins with `#` in Ensembl-format, region or HGVS input is a
/// comment here, where ensembl-io's list parsers name no metadata line
/// (ListBasedParser.pm:60, Parser/VEP_input.pm:61) and Ensembl VEP reads a
/// record. The run says so once, at the end of the input, with the number of
/// such lines and the first one's line number, with and without `--dont_skip`;
/// the line stands in the warnings file under `--quiet` and on standard error
/// as well without it. The rows are those of the same input without its `#`
/// lines, and an input without one writes no warnings file.
#[test]
fn hash_lines_read_as_comments_are_counted_once_a_run() {
    let dir = tempfile::tempdir().unwrap();
    let formats: [(&str, [&str; 2]); 3] = [
        (
            "ensembl",
            [
                "21\t10538885\t10538885\tA/C\t1\tfirst",
                "21\t14221107\t14221107\tC/A\t1\tsecond",
            ],
        ),
        (
            "region",
            ["21:10538885-10538885:1/C", "21:14221107-14221107:1/A"],
        ),
        ("hgvs", ["21:g.10538885A>C", "21:g.14221107C>A"]),
    ];
    for (format, [first, second]) in formats {
        let plain = write(
            dir.path(),
            &format!("plain.{format}"),
            &format!("{first}\n{second}\n"),
        );
        let expected = annotate(
            dir.path(),
            &format!("plain_{format}"),
            &plain,
            &["--format", format, "--quiet"],
        );
        assert!(!expected.rows.is_empty(), "{format}");
        assert_eq!(expected.warnings, "", "{format}");

        let commented = write(
            dir.path(),
            &format!("commented.{format}"),
            &format!("{first}\n# a header\n#another\n{second}\n# the last line\n"),
        );
        let notice = "WARNING: Read 3 lines that begin with '#' as comments, the first at line 2";
        for flags in [&[][..], &["--dont_skip"][..]] {
            let mut extra = vec!["--format", format, "--quiet"];
            extra.extend(flags);
            let quiet = annotate(
                dir.path(),
                &format!("commented_{format}_{}", flags.len()),
                &commented,
                &extra,
            );
            assert_eq!(quiet.rows, expected.rows, "{format} {flags:?}");
            assert_eq!(quiet.warning_lines(), [notice], "{format} {flags:?}");
            assert!(
                quiet.warned_on_standard_error("with '#' as").is_empty(),
                "{format} {flags:?}: {}",
                text(&quiet.output.stderr)
            );
        }
        let loud = annotate(
            dir.path(),
            &format!("loud_{format}"),
            &commented,
            &["--format", format],
        );
        assert_eq!(loud.rows, expected.rows, "{format}");
        assert_eq!(loud.warning_lines(), [notice], "{format}");
        assert_eq!(
            loud.warned_on_standard_error("with '#' as"),
            [notice],
            "{format}"
        );

        let single = write(
            dir.path(),
            &format!("single.{format}"),
            &format!("# a header\n{first}\n{second}\n"),
        );
        let one = annotate(
            dir.path(),
            &format!("single_{format}"),
            &single,
            &["--format", format, "--quiet"],
        );
        assert_eq!(one.rows, expected.rows, "{format}");
        assert_eq!(
            one.warning_lines(),
            ["WARNING: Read 1 line that begins with '#' as a comment, at line 1"],
            "{format}"
        );
    }
}

/// The `INFO:` lines of standard error that say what the run left out of
/// what it was given to read: a FASTA file passed over, a chromosome the
/// FASTA lacks, a notation on an edited transcript. Any other line of the
/// run is no concern of these tests.
fn notices(output: &Output) -> Vec<String> {
    text(&output.stderr)
        .lines()
        .filter(|line| {
            line.starts_with("INFO: FASTA file ")
                || line.starts_with("INFO: The FASTA holds no sequence ")
                || (line.starts_with("INFO: ") && line.contains("edited away from the genome"))
        })
        .map(str::to_string)
        .collect()
}

/// Copies the files and directories of `from` into the new directory `to`.
fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// A sequence of 120 bases named `22` in a plain FASTA at `path`, with its
/// `.fai` beside it when `indexed`.
fn write_small_fasta(path: &Path, indexed: bool) {
    std::fs::write(
        path,
        format!(">22\n{}\n{}\n", "ACGT".repeat(15), "ACGT".repeat(15)),
    )
    .unwrap();
    if indexed {
        std::fs::write(format!("{}.fai", path.display()), "22\t120\t4\t60\t61\n").unwrap();
    }
}

/// Runs the binary on the identifiers corpus's records with `flags`, writing
/// `<name>.txt` under `dir`; exit status 0 is required.
fn run_records(dir: &Path, name: &str, flags: &[&str]) -> Run {
    run_input(dir, name, &corpus().join("variants.vcf"), flags)
}

/// Runs the binary on `input` with `flags`, writing `<name>.txt` under `dir`;
/// exit status 0 is required.
fn run_input(dir: &Path, name: &str, input: &Path, flags: &[&str]) -> Run {
    let out = dir.join(format!("{name}.txt"));
    let mut command = common::isolated_vep(dir);
    command
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .args([
            "--assembly",
            "GRCh38",
            "--offline",
            "--force_overwrite",
            "--no_stats",
        ])
        .args(flags);
    let output = command.output().unwrap();
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let rows = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect();
    let warnings =
        std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap_or_default();
    Run {
        output,
        rows,
        warnings,
    }
}

/// A `.fa` or `.fa.gz` in the cache directory without the index it needs is
/// passed over by the detection that finds the run's FASTA there, where
/// Ensembl VEP takes the file and writes its index
/// (`CacheDir.pm:235-251`, `Utils/FastaSequence.pm:308-358`). The run names
/// each such file once on standard error, with the index it lacks, what
/// writes it, and the run's reference FASTA without the file: none, or the
/// first indexed file. The line is a notice, withheld under `--quiet` and absent
/// from the warnings file, and the rows are those of the same run under
/// `--no_fasta`. A directory whose FASTA files all have their indexes prints
/// no such line.
#[test]
fn a_cache_directory_fasta_without_its_index_is_named_once() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("cache");
    copy_dir(&corpus().join("json_cache"), &cache);
    let json_cache = cache.to_str().unwrap();
    let passed_over = |file: &str, lacks: &str, reads: &str| {
        format!(
            "INFO: FASTA file {} in the cache directory was passed over ({lacks}); {reads}",
            cache.join(file).display()
        )
    };
    let no_fai = |file: &str| {
        format!(
            "FASTA index {}.fai not found: index the FASTA with samtools faidx",
            cache.join(file).display()
        )
    };

    let expected = run_records(
        dir.path(),
        "no_fasta",
        &["--json_cache", json_cache, "--no_fasta"],
    );
    assert!(!expected.rows.is_empty());
    assert!(notices(&expected.output).is_empty());

    write_small_fasta(&cache.join("a.fa"), false);
    let alone = run_records(dir.path(), "alone", &["--json_cache", json_cache]);
    assert_eq!(
        notices(&alone.output),
        [passed_over(
            "a.fa",
            &no_fai("a.fa"),
            "the run has no reference FASTA"
        )]
    );
    assert_eq!(alone.rows, expected.rows);
    assert_eq!(alone.warnings, "");
    let quiet = run_records(
        dir.path(),
        "quiet",
        &["--json_cache", json_cache, "--quiet"],
    );
    assert!(notices(&quiet.output).is_empty());
    assert_eq!(quiet.rows, expected.rows);
    assert_eq!(quiet.warnings, "");
    let named = run_records(
        dir.path(),
        "named",
        &["--json_cache", json_cache, "--no_fasta"],
    );
    assert!(
        notices(&named.output).is_empty(),
        "--no_fasta looks for no FASTA in the cache directory"
    );

    // A cache without its checksum manifest reads an indexed FASTA that the
    // manifest does not list.
    std::fs::remove_file(cache.join("CHECKSUMS.sha256")).unwrap();
    write_small_fasta(&cache.join("b.fa"), true);
    let beside = run_records(dir.path(), "beside", &["--json_cache", json_cache]);
    assert_eq!(
        notices(&beside.output),
        [passed_over(
            "a.fa",
            &no_fai("a.fa"),
            &format!(
                "the run's reference FASTA is {}",
                cache.join("b.fa").display()
            )
        )]
    );
    assert_eq!(beside.rows, expected.rows);

    std::fs::remove_file(cache.join("a.fa")).unwrap();
    let indexed = run_records(dir.path(), "indexed", &["--json_cache", json_cache]);
    assert!(
        notices(&indexed.output).is_empty(),
        "{}",
        text(&indexed.output.stderr)
    );
    assert_eq!(indexed.rows, expected.rows);
}

/// The ancestral-allele corpus: a chromosome 21 cache slice and a slice of
/// the ancestral FASTA the AncestralAllele plugin reads.
fn ancestral_corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-plugin-ancestral")
}

/// The line for a file passed over in the cache directory speaks of the
/// run's reference FASTA, which a plugin's own FASTA is not: a run whose
/// AncestralAllele plugin reads its ancestral FASTA has no reference FASTA
/// when the one file of the cache directory lacks its index, the line says
/// so, and the plugin's values are those of the same run on the cache
/// without that file.
#[test]
fn a_plugin_that_reads_its_own_fasta_gives_the_run_no_reference_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("cache");
    copy_dir(&ancestral_corpus().join("json_cache"), &cache);
    let variants = ancestral_corpus().join("variants.vcf");
    let ancestral = ancestral_corpus().join("ancestral/homo_sapiens_ancestor_21.slice.fa.gz");
    let plugin = format!("AncestralAllele,{}", ancestral.display());
    let flags = [
        "--json_cache",
        cache.to_str().unwrap(),
        "--plugin",
        plugin.as_str(),
    ];

    let expected = run_input(dir.path(), "indexed", &variants, &flags);
    assert!(notices(&expected.output).is_empty());
    let read_from_the_plugin_fasta = |row: &String| {
        row.split("AA=")
            .nth(1)
            .and_then(|value| value.chars().next())
            .is_some_and(|base| base.is_ascii_alphabetic())
    };
    assert!(
        expected.rows.iter().any(read_from_the_plugin_fasta),
        "{:?}",
        expected.rows.first()
    );

    write_small_fasta(&cache.join("a.fa"), false);
    let beside = run_input(dir.path(), "beside", &variants, &flags);
    assert_eq!(
        notices(&beside.output),
        [format!(
            "INFO: FASTA file {a} in the cache directory was passed over (FASTA index {a}.fai not found: index the FASTA with samtools faidx); the run has no reference FASTA",
            a = cache.join("a.fa").display()
        )]
    );
    assert_eq!(beside.rows, expected.rows);
}

/// The line names the run's reference FASTA before the file is opened, so it
/// holds when the cache's checksum manifest then refuses that file: an
/// indexed FASTA the manifest does not list, beside one passed over, is
/// named as the run's reference FASTA, and the run ends on the manifest's
/// refusal of it with exit status 1 and no output.
#[test]
fn a_cache_fasta_the_checksum_manifest_refuses_is_named_before_the_refusal() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("cache");
    copy_dir(&corpus().join("json_cache"), &cache);
    assert!(cache.join("CHECKSUMS.sha256").is_file());
    write_small_fasta(&cache.join("a.fa"), false);
    write_small_fasta(&cache.join("b.fa"), true);
    let out = dir.path().join("refused.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(corpus().join("variants.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(&cache)
        .args([
            "--assembly",
            "GRCh38",
            "--offline",
            "--force_overwrite",
            "--no_stats",
        ])
        .output()
        .unwrap();
    let stderr = text(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    let named = format!(
        "INFO: FASTA file {a} in the cache directory was passed over (FASTA index {a}.fai not found: index the FASTA with samtools faidx); the run's reference FASTA is {b}",
        a = cache.join("a.fa").display(),
        b = cache.join("b.fa").display()
    );
    assert_eq!(notices(&output), [named.as_str()]);
    let refusal = format!(
        "ERROR: [cache] checksum mismatch for b.fa in {}",
        cache.display()
    );
    let lines: Vec<&str> = stderr
        .lines()
        .filter(|line| *line == named || line.starts_with("ERROR:"))
        .collect();
    assert_eq!(lines, [named.as_str(), refusal.as_str()], "{stderr}");
    assert!(!out.exists());
}

/// A `--fasta` directory is read as its indexed `.fa` and `.fa.gz` files; a
/// file in it without the index it needs is passed over, and the run names
/// each such file once on standard error with the index it lacks: the `.fai`
/// of a plain file, the `.gzi` of a BGZF file. The rows are those of the run
/// given the indexed file alone, and a directory of indexed files prints no
/// such line.
#[test]
fn a_fasta_directory_file_without_its_index_is_named_once() {
    let dir = tempfile::tempdir().unwrap();
    let fastas = dir.path().join("fastas");
    std::fs::create_dir(&fastas).unwrap();
    for suffix in ["", ".fai", ".gzi"] {
        std::fs::copy(
            corpus().join(format!("reference.fa.gz{suffix}")),
            fastas.join(format!("reference.fa.gz{suffix}")),
        )
        .unwrap();
    }
    let json_cache = corpus().join("json_cache");
    let json_cache = json_cache.to_str().unwrap();
    let reference = corpus().join("reference.fa.gz");
    let expected = run_records(
        dir.path(),
        "file",
        &[
            "--json_cache",
            json_cache,
            "--fasta",
            reference.to_str().unwrap(),
            "--hgvs",
        ],
    );
    assert!(!expected.rows.is_empty());
    assert!(notices(&expected.output).is_empty());

    let flags = [
        "--json_cache",
        json_cache,
        "--fasta",
        fastas.to_str().unwrap(),
        "--hgvs",
    ];
    let indexed = run_records(dir.path(), "indexed", &flags);
    assert!(
        notices(&indexed.output).is_empty(),
        "{}",
        text(&indexed.output.stderr)
    );
    assert_eq!(indexed.rows, expected.rows);

    write_small_fasta(&fastas.join("plain.fa"), false);
    std::fs::copy(&reference, fastas.join("second.fa.gz")).unwrap();
    std::fs::copy(
        corpus().join("reference.fa.gz.fai"),
        fastas.join("second.fa.gz.fai"),
    )
    .unwrap();
    let shown = |file: &str| fastas.join(file).display().to_string();
    let lines = [
        format!(
            "INFO: FASTA file {} in the --fasta directory was passed over (FASTA index {} not found: index the FASTA with samtools faidx); its sequences are not read",
            shown("plain.fa"),
            shown("plain.fa.fai")
        ),
        format!(
            "INFO: FASTA file {} in the --fasta directory was passed over (BGZF index {} not found: index the compressed FASTA with samtools faidx, which writes the .fai and the .gzi); its sequences are not read",
            shown("second.fa.gz"),
            shown("second.fa.gz.gzi")
        ),
    ];
    let passed_over = run_records(dir.path(), "passed_over", &flags);
    assert_eq!(notices(&passed_over.output), lines);
    assert_eq!(passed_over.rows, expected.rows);
    assert_eq!(passed_over.warnings, expected.warnings);
    let mut quietly = flags.to_vec();
    quietly.push("--quiet");
    let quiet = run_records(dir.path(), "quiet", &quietly);
    assert!(notices(&quiet.output).is_empty());
    assert_eq!(quiet.rows, expected.rows);
}

/// Under `--hgvs` a substitution or a deletion on a chromosome the FASTA holds
/// no sequence for gets no HGVSc, where Ensembl VEP reads `N` per base there
/// and prints one (`Utils/FastaSequence.pm:597-601`; its warning at `:599` is
/// raised above the verbosity a run sets, `BaseRunner.pm:491`, so neither
/// engine warns). The run names the chromosome in one notice at its end,
/// withheld under `--quiet` and absent from the warnings file; the rows are
/// the same with and without the notice, and a FASTA that holds the
/// chromosome, or a run under neither `--hgvs` nor `--hgvsg`, prints none.
#[test]
fn a_chromosome_the_fasta_lacks_is_named_once_under_hgvs() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let json_cache = corpus().join("json_cache");
    let json_cache = json_cache.to_str().unwrap();
    let lacking = [
        "--json_cache",
        json_cache,
        "--fasta",
        other.to_str().unwrap(),
    ];
    let notice = NO_SEQUENCE_FOR_21;

    let mut hgvs = lacking.to_vec();
    hgvs.push("--hgvs");
    let loud = run_records(dir.path(), "loud", &hgvs);
    assert_eq!(notices(&loud.output), [notice]);
    assert_eq!(loud.warnings, "");
    assert!(!loud.rows.is_empty());
    let mut quietly = hgvs.clone();
    quietly.push("--quiet");
    let quiet = run_records(dir.path(), "quiet", &quietly);
    assert!(notices(&quiet.output).is_empty());
    assert_eq!(quiet.rows, loud.rows);
    assert_eq!(quiet.warnings, "");

    let without_hgvs = run_records(dir.path(), "without_hgvs", &lacking);
    assert!(notices(&without_hgvs.output).is_empty());

    let reference = corpus().join("reference.fa.gz");
    let held = run_records(
        dir.path(),
        "held",
        &[
            "--json_cache",
            json_cache,
            "--fasta",
            reference.to_str().unwrap(),
            "--hgvs",
        ],
    );
    assert!(
        notices(&held.output).is_empty(),
        "{}",
        text(&held.output.stderr)
    );
}

/// The identifiers corpus's one-base substitution in an intron of nine
/// transcripts (line 6 of its `variants.vcf`), alone in a VCF under `dir`.
fn intron_substitution(dir: &Path) -> PathBuf {
    write(
        dir,
        "intron_snv.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t10538885\tintron_snv_1\tA\tC\t.\t.\t.\n",
    )
}

/// Under `--hgvsg` a record on a chromosome the FASTA holds no sequence for
/// gets no HGVSg: both ends of a record must be bases the FASTA holds. Ensembl
/// VEP reads `N` per base there (`Utils/FastaSequence.pm:597-601`) and, for a
/// record that takes a transcript's slice (`AnnotationType/Transcript.pm:122`),
/// cuts the reference allele of its notation from those bases
/// (`VariationFeature.pm:1927`, `Utils/Sequence.pm:526`): over a FASTA holding
/// chromosome 22 alone Ensembl VEP 116.3 prints `21:g.10538885N>C` on each of
/// this record's nine rows. Over the corpus reference both engines print
/// `21:g.10538885A>C` there (`expected_noshift/default.txt.gz` holds Ensembl
/// VEP's rows). No warning is written. `docs/cache-setup.md` states the
/// difference.
#[test]
fn a_contig_the_fasta_lacks_gives_no_hgvsg() {
    let dir = tempfile::tempdir().unwrap();
    let input = intron_substitution(dir.path());
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    let rows_over = |name: &str, fasta: &Path| {
        run_input(
            dir.path(),
            name,
            &input,
            &[
                "--json_cache",
                json_cache.to_str().unwrap(),
                "--fasta",
                fasta.to_str().unwrap(),
                "--hgvsg",
                "--tab",
                "--fields",
                "Uploaded_variation,Location,Allele,Feature,HGVSg",
            ],
        )
    };
    let expected = |hgvsg: &str| -> Vec<String> {
        [
            "ENST00000427445",
            "ENST00000612746",
            "ENST00000612957",
            "ENST00000618007",
            "ENST00000622113",
            "ENST00001091160",
            "ENST00001091161",
            "ENST00001091162",
            "ENST00001131865",
        ]
        .iter()
        .map(|transcript| format!("intron_snv_1\t21:10538885\tC\t{transcript}\t{hgvsg}"))
        .collect()
    };

    let lacking = rows_over("lacking", &other);
    assert_eq!(lacking.rows, expected("-"));
    assert_eq!(lacking.warnings, "");
    let held = rows_over("held", &reference);
    assert_eq!(held.rows, expected("21:g.10538885A>C"));
    assert_eq!(held.warnings, "");
}

/// The notice of a run under `--hgvsg` alone whose FASTA holds no sequence for
/// chromosome 21.
const NO_SEQUENCE_FOR_21_HGVSG: &str = "INFO: The FASTA holds no sequence for chromosome 21, \
    which the input has records on: a record there gets no HGVSg";

/// The notice of a run under `--hgvs` and `--hgvsg` whose FASTA holds no
/// sequence for chromosome 21.
const NO_SEQUENCE_FOR_21_BOTH: &str = "INFO: The FASTA holds no sequence for chromosome 21, \
    which the input has records on: a record there gets no HGVSg, and a substitution or a \
    deletion no HGVSc, nor does an insertion written outside a transcript";

/// A run under `--hgvsg` names a chromosome the FASTA holds no sequence for
/// in the notice a run under `--hgvs` prints at its end: under `--hgvsg`
/// alone the line names HGVSg, and under both flags the one line names both
/// notations. It is withheld under `--quiet` and absent from the warnings
/// file, the rows are the same with and without it, and a FASTA that holds
/// the chromosome prints none.
#[test]
fn a_chromosome_the_fasta_lacks_is_named_once_under_hgvsg() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let json_cache = corpus().join("json_cache");
    let json_cache = json_cache.to_str().unwrap();
    let lacking = [
        "--json_cache",
        json_cache,
        "--fasta",
        other.to_str().unwrap(),
    ];
    let with = |more: &[&'static str]| {
        let mut flags = lacking.to_vec();
        flags.extend_from_slice(more);
        flags
    };
    let genomic = |run: &Run| run.rows.iter().filter(|row| row.contains("HGVSg=")).count();

    let loud = run_records(dir.path(), "loud", &with(&["--hgvsg"]));
    assert_eq!(notices(&loud.output), [NO_SEQUENCE_FOR_21_HGVSG]);
    assert_eq!(loud.warnings, "");
    assert!(!loud.rows.is_empty());
    assert_eq!(genomic(&loud), 0);
    let quiet = run_records(dir.path(), "quiet", &with(&["--hgvsg", "--quiet"]));
    assert!(notices(&quiet.output).is_empty());
    assert_eq!(quiet.rows, loud.rows);
    assert_eq!(quiet.warnings, "");

    let both = run_records(dir.path(), "both", &with(&["--hgvs", "--hgvsg"]));
    assert_eq!(notices(&both.output), [NO_SEQUENCE_FOR_21_BOTH]);
    assert_eq!(genomic(&both), 0);
    let coding = run_records(dir.path(), "coding", &with(&["--hgvs"]));
    assert_eq!(notices(&coding.output), [NO_SEQUENCE_FOR_21]);

    let reference = corpus().join("reference.fa.gz");
    let held = run_records(
        dir.path(),
        "held",
        &[
            "--json_cache",
            json_cache,
            "--fasta",
            reference.to_str().unwrap(),
            "--hgvsg",
        ],
    );
    assert!(
        notices(&held.output).is_empty(),
        "{}",
        text(&held.output.stderr)
    );
    assert_eq!(genomic(&held), held.rows.len());
    assert_eq!(held.rows.len(), loud.rows.len());
}

/// One run of the binary under the HGVS flags: its name, the notation flags,
/// the flags that settle what its output writes, and the line it prints for
/// a chromosome the FASTA lacks, `None` where it prints none.
type NotationCase = (
    &'static str,
    &'static [&'static str],
    &'static [&'static str],
    Option<&'static str>,
);

/// Whether the rows of `run` hold an HGVSg (`:g.`) and whether they hold an
/// HGVSc (`:c.` or `:n.`).
fn notations_written(run: &Run) -> (bool, bool) {
    let holds = |mark: &str| run.rows.iter().any(|row| row.contains(mark));
    (holds(":g."), holds(":c.") || holds(":n."))
}

/// The line for a chromosome the FASTA holds no sequence for names a notation
/// only where the run's output writes it. A tab or a VCF output writes the
/// fields its `--fields` list names, and `--most_severe` and `--summary`
/// write one row a record before any row that carries a notation is built
/// (`OutputFactory.pm:387`; HGVSg is set at `:1393` and HGVSc at `:1763`), so
/// the line leaves out the notation such a run does not write, and a run that
/// writes neither prints no line. The default and JSON outputs read no
/// `--fields` list and write both. Each case is held to what its flags write
/// over the corpus FASTA, which holds the chromosome.
#[test]
fn a_chromosome_the_fasta_lacks_is_named_for_the_notations_the_output_writes() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    let over = |name: &str, fasta: &Path, more: &[&str]| {
        let mut flags = vec![
            "--json_cache",
            json_cache.to_str().unwrap(),
            "--fasta",
            fasta.to_str().unwrap(),
        ];
        flags.extend_from_slice(more);
        run_records(dir.path(), name, &flags)
    };
    const BOTH_FLAGS: [&str; 2] = ["--hgvs", "--hgvsg"];
    let cases: [NotationCase; 14] = [
        (
            "tab_without_hgvsc",
            &["--hgvs"],
            &["--tab", "--fields", "Uploaded_variation,Feature,HGVSp"],
            None,
        ),
        (
            "vcf_without_hgvsc",
            &["--hgvs"],
            &["--vcf", "--fields", "Allele,Feature,HGVSp"],
            None,
        ),
        (
            "tab_without_hgvsg",
            &["--hgvsg"],
            &["--tab", "--fields", "Uploaded_variation,Feature"],
            None,
        ),
        (
            "tab_hgvsg_alone",
            &BOTH_FLAGS,
            &["--tab", "--fields", "Uploaded_variation,Feature,HGVSg"],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        (
            "tab_hgvsc_alone",
            &BOTH_FLAGS,
            &["--tab", "--fields", "Uploaded_variation,Feature,HGVSc"],
            Some(NO_SEQUENCE_FOR_21),
        ),
        (
            "vcf_hgvsc_alone",
            &BOTH_FLAGS,
            &["--vcf", "--fields", "Allele,Feature,HGVSc"],
            Some(NO_SEQUENCE_FOR_21),
        ),
        (
            "tab_both_fields",
            &BOTH_FLAGS,
            &["--tab", "--fields", "Uploaded_variation,HGVSc,HGVSg"],
            Some(NO_SEQUENCE_FOR_21_BOTH),
        ),
        (
            "tab_every_field",
            &BOTH_FLAGS,
            &["--tab"],
            Some(NO_SEQUENCE_FOR_21_BOTH),
        ),
        ("most_severe", &BOTH_FLAGS, &["--most_severe"], None),
        ("summary", &BOTH_FLAGS, &["--summary"], None),
        (
            "most_severe_tab",
            &BOTH_FLAGS,
            &["--most_severe", "--tab"],
            None,
        ),
        ("summary_json", &BOTH_FLAGS, &["--summary", "--json"], None),
        (
            "default_reads_no_field_list",
            &BOTH_FLAGS,
            &["--fields", "Uploaded_variation,Feature"],
            Some(NO_SEQUENCE_FOR_21_BOTH),
        ),
        (
            "json",
            &BOTH_FLAGS,
            &["--json"],
            Some(NO_SEQUENCE_FOR_21_BOTH),
        ),
    ];
    for (name, notation_flags, output_flags, line) in cases {
        let flags = [notation_flags, output_flags].concat();
        let held = over(&format!("{name}_held"), &reference, &flags);
        assert!(notices(&held.output).is_empty(), "{name}");
        let written = match line {
            None => (false, false),
            Some(NO_SEQUENCE_FOR_21_HGVSG) => (true, false),
            Some(NO_SEQUENCE_FOR_21) => (false, true),
            Some(_) => (true, true),
        };
        assert_eq!(notations_written(&held), written, "{name}");
        let lacking = over(&format!("{name}_lacking"), &other, &flags);
        let expected: Vec<&str> = line.into_iter().collect();
        assert_eq!(notices(&lacking.output), expected, "{name}");
        assert_eq!(lacking.warnings, "", "{name}");
        assert_eq!(lacking.rows.len(), held.rows.len(), "{name}");
    }
}

/// One run of the binary whose flags settle the rows its output keeps: its
/// name, its flags, and the line it prints for a chromosome the FASTA lacks,
/// `None` where it prints none.
type KeptRowsCase = (&'static str, &'static [&'static str], Option<&'static str>);

/// How many HGVSg (`:g.`) and how many HGVSc (`:c.` or `:n.`) the rows of
/// `run` hold.
fn notations_held(run: &Run) -> (usize, usize) {
    let held = |mark: &str| -> usize { run.rows.iter().map(|row| row.matches(mark).count()).sum() };
    (held(":g."), held(":c.") + held(":n."))
}

/// The filter no transcript of the identifiers corpus passes.
const NO_TRANSCRIPT: [&str; 2] = ["--transcript_filter", "stable_id is ENST_none"];

/// The line for a chromosome the FASTA holds no sequence for names the
/// chromosome, and a notation, only where a row the output keeps on that
/// chromosome is without the notation for want of the sequence. `--chr`
/// leaves the records of another chromosome out before any row is built; a
/// `--transcript_filter` no transcript passes leaves each record its
/// intergenic row, which holds an HGVSg and never an HGVSc; `--no_intergenic`
/// then leaves no row; and `--coding_only` keeps rows that hold an HGVSc. Each
/// case is held to what the same flags write over the corpus FASTA, which
/// holds the chromosome: the line names a notation exactly where those rows
/// hold more of it.
#[test]
fn a_chromosome_the_fasta_lacks_is_named_where_a_kept_row_is_without_the_notation() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    let over = |name: &str, fasta: &Path, more: &[&str]| {
        let mut flags = vec![
            "--json_cache",
            json_cache.to_str().unwrap(),
            "--fasta",
            fasta.to_str().unwrap(),
        ];
        flags.extend_from_slice(more);
        run_records(dir.path(), name, &flags)
    };
    let cases: [KeptRowsCase; 11] = [
        (
            "chr_leaves_the_records_out",
            &["--hgvs", "--hgvsg", "--chr", "22"],
            None,
        ),
        (
            "no_transcript_row_under_both",
            &["--hgvs", "--hgvsg", NO_TRANSCRIPT[0], NO_TRANSCRIPT[1]],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        (
            "no_transcript_row_under_hgvs",
            &["--hgvs", NO_TRANSCRIPT[0], NO_TRANSCRIPT[1]],
            None,
        ),
        (
            "no_transcript_row_under_hgvsg",
            &["--hgvsg", NO_TRANSCRIPT[0], NO_TRANSCRIPT[1]],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        (
            "no_transcript_row_in_json",
            &[
                "--hgvs",
                "--hgvsg",
                NO_TRANSCRIPT[0],
                NO_TRANSCRIPT[1],
                "--json",
            ],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        (
            "no_transcript_row_in_vcf",
            &[
                "--hgvs",
                "--hgvsg",
                NO_TRANSCRIPT[0],
                NO_TRANSCRIPT[1],
                "--vcf",
            ],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        (
            "no_transcript_row_in_tab",
            &[
                "--hgvs",
                "--hgvsg",
                NO_TRANSCRIPT[0],
                NO_TRANSCRIPT[1],
                "--tab",
            ],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        (
            "no_row_at_all",
            &[
                "--hgvs",
                "--hgvsg",
                NO_TRANSCRIPT[0],
                NO_TRANSCRIPT[1],
                "--no_intergenic",
            ],
            None,
        ),
        (
            "coding_rows_under_hgvs",
            &["--hgvs", "--coding_only"],
            Some(NO_SEQUENCE_FOR_21),
        ),
        (
            "coding_rows_under_both",
            &["--hgvs", "--hgvsg", "--coding_only"],
            Some(NO_SEQUENCE_FOR_21_BOTH),
        ),
        (
            "picked_rows_under_both",
            &["--hgvs", "--hgvsg", "--pick"],
            Some(NO_SEQUENCE_FOR_21_BOTH),
        ),
    ];
    for (name, flags, line) in cases {
        let held = over(&format!("{name}_held"), &reference, flags);
        assert!(notices(&held.output).is_empty(), "{name}");
        let lacking = over(&format!("{name}_lacking"), &other, flags);
        assert_eq!(lacking.rows.len(), held.rows.len(), "{name}");
        assert_eq!(lacking.warnings, "", "{name}");
        let (held_genomic, held_coding) = notations_held(&held);
        let (genomic, coding) = notations_held(&lacking);
        let without = match line {
            None => (false, false),
            Some(NO_SEQUENCE_FOR_21_HGVSG) => (true, false),
            Some(NO_SEQUENCE_FOR_21) => (false, true),
            Some(_) => (true, true),
        };
        assert_eq!(
            (genomic < held_genomic, coding < held_coding),
            without,
            "{name}"
        );
        let expected: Vec<&str> = line.into_iter().collect();
        assert_eq!(notices(&lacking.output), expected, "{name}");
    }
}

/// The identifiers corpus's records that add bases to the reference and take
/// none away (a one-base REF that opens the ALT), alone in a VCF under `dir`,
/// with how many there are.
fn insertions_alone(dir: &Path) -> (PathBuf, usize) {
    let records = std::fs::read_to_string(corpus().join("variants.vcf")).unwrap();
    let kept: Vec<&str> = records
        .lines()
        .filter(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            line.starts_with('#')
                || (cells[3].len() == 1
                    && cells[4].len() > 1
                    && cells[4].starts_with(cells[3])
                    && !cells[4].contains(','))
        })
        .collect();
    let insertions = kept.iter().filter(|line| !line.starts_with('#')).count();
    (write(dir, "insertions.vcf", &lines_of(&kept)), insertions)
}

/// An insertion names no reference base, so on a chromosome the FASTA holds
/// no sequence for it has its HGVSc on a transcript it lies within, and no
/// HGVSg. A run under `--hgvs` that keeps such a row alone prints no line,
/// and under `--hgvsg` as well its line names HGVSg alone. A row of bases
/// written before a transcript's first base holds an HGVSc only where the 3'
/// rule places the bases on the transcript, along bases such a run cannot
/// read: a run that keeps such a row names HGVSc too, as the corpus's
/// insertions do, 55 of whose rows lie upstream of their transcript, though
/// each of their rows holds the HGVSc it holds over the corpus FASTA.
#[test]
fn insertions_draw_the_line_for_hgvsc_where_a_kept_row_lies_before_its_transcript() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let (input, insertions) = insertions_alone(dir.path());
    assert!(insertions > 0);
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    let over = |name: &str, fasta: &Path, more: &[&str]| {
        let mut flags = vec![
            "--json_cache",
            json_cache.to_str().unwrap(),
            "--fasta",
            fasta.to_str().unwrap(),
        ];
        flags.extend_from_slice(more);
        run_input(dir.path(), name, &input, &flags)
    };

    // The rows of the three insertions that lie within ENST00001131865.
    let within = [
        "--transcript_filter",
        "stable_id is ENST00001131865",
        "--no_intergenic",
    ];
    let flags = |more: &[&'static str]| [&within[..], more].concat();
    let held = over("within_held", &reference, &flags(&["--hgvs", "--hgvsg"]));
    assert_eq!(held.rows.len(), 3);
    assert!(held
        .rows
        .iter()
        .all(|row| !row.contains("stream_gene_variant")));
    let (held_genomic, held_coding) = notations_held(&held);
    assert_eq!(
        (held_genomic, held_coding),
        (held.rows.len(), held.rows.len())
    );
    let coding = over("within_coding", &other, &flags(&["--hgvs"]));
    assert_eq!(notations_held(&coding), (0, held_coding));
    assert!(
        notices(&coding.output).is_empty(),
        "{}",
        text(&coding.output.stderr)
    );
    let both = over("within_both", &other, &flags(&["--hgvs", "--hgvsg"]));
    assert_eq!(notations_held(&both), (0, held_coding));
    assert_eq!(notices(&both.output), [NO_SEQUENCE_FOR_21_HGVSG]);

    // Every row of the fifteen records, the upstream rows among them.
    let held = over("held", &reference, &["--hgvs", "--hgvsg"]);
    let (held_genomic, held_coding) = notations_held(&held);
    assert!(held_genomic > 0 && held_coding > 0);
    let upstream = |run: &Run| {
        run.rows
            .iter()
            .filter(|row| row.contains("upstream_gene_variant"))
            .count()
    };
    assert_eq!(upstream(&held), 55);
    let coding = over("coding", &other, &["--hgvs"]);
    assert_eq!(notations_held(&coding), (0, held_coding));
    assert_eq!(upstream(&coding), 55);
    assert_eq!(notices(&coding.output), [NO_SEQUENCE_FOR_21]);
    let both = over("both", &other, &["--hgvs", "--hgvsg"]);
    assert_eq!(notations_held(&both), (0, held_coding));
    assert_eq!(notices(&both.output), [NO_SEQUENCE_FOR_21_BOTH]);
}

/// Two rules decide an HGVSc from bases beside the record, at a transcript's
/// end: bases added before a transcript's first base, or directly after its
/// last, take the notation of their last placement where the run of their
/// placements reaches the transcript, and an allele that keeps bases of its
/// reference and deletes others over a transcript's first or last base has a
/// notation where the deleted bases can be written wholly inside the
/// transcript. On a chromosome the FASTA holds no sequence for neither rule
/// reads a base, so the row has no HGVSc, and whether it holds one over the
/// chromosome's own bases cannot be told: the line names HGVSc for every
/// such row some bases would give a notation, those that hold one over the
/// corpus FASTA and those that hold none there alike. A row that no bases
/// give a notation draws no line: bases added past the base after the last
/// base, a deletion the record's own bases hold over the last base, and
/// bases added before the first base under `--shift_hgvs 0` or
/// `--shift_3prime 1`, or in a record whose class does not shift. Each run
/// keeps one transcript's rows alone, and is held to what the same flags
/// write over the corpus FASTA.
#[test]
fn a_row_whose_hgvsc_turns_on_bases_beside_the_record_draws_the_line() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    for (name, record, transcript, more, held_coding, named) in [
        // `GT` before `n.1` of a forward snRNA whose first bases repeat it:
        // `n.2_3dup` over the corpus FASTA.
        (
            "added_before_the_first_base",
            "21\t13621576\tbefore_dup_fwd\tT\tTGT\t.\t.\t.\n",
            "ENST00000384745",
            &[][..],
            1,
            true,
        ),
        (
            "added_before_the_first_base_unshifted",
            "21\t13621576\tbefore_dup_fwd\tT\tTGT\t.\t.\t.\n",
            "ENST00000384745",
            &["--shift_hgvs", "0"][..],
            0,
            false,
        ),
        // The same bases in a record with a `*` allele, whose class does not
        // shift (`vep_effects::shift_motif`): they stay before the first
        // base over any sequence.
        (
            "added_before_the_first_base_in_a_record_that_does_not_shift",
            "21\t13621576\tbefore_dup_fwd\tT\tTGT,*\t.\t.\t.\n",
            "ENST00000384745",
            &[][..],
            0,
            false,
        ),
        // `A` directly after the last base of a forward transcript, an `A`
        // that an `A` follows: no notation over the corpus FASTA.
        (
            "added_after_the_last_base",
            "21\t17614113\ton_fwd\tA\tAA\t.\t.\t.\n",
            "ENST00000630155",
            &[][..],
            0,
            true,
        ),
        // Under `--shift_3prime 1` the record's own shift reads the same
        // bases: over the corpus FASTA the row is annotated on the exon the
        // shift carries the bases to, with its notation, and here where the
        // record writes them, outside the transcript, asked for none.
        (
            "added_before_the_first_base_under_the_record_shift",
            "21\t13621576\tbefore_dup_fwd\tT\tTGT\t.\t.\t.\n",
            "ENST00000384745",
            &["--shift_3prime", "1"][..],
            1,
            false,
        ),
        // `A` between the last base of a forward transcript and the base
        // after it, in a record that keeps the last two bases and that base:
        // `n.296dup` over the corpus FASTA.
        (
            "added_after_the_last_base_to_a_kept_reference",
            "21\t10576585\tover_one_fwd\tTAAT\tTAAAT,TA\t.\t.\t.\n",
            "ENST00000616920",
            &[][..],
            1,
            true,
        ),
        (
            "added_past_the_base_after_the_last",
            "21\t17614114\tbeyond_fwd\tA\tAA\t.\t.\t.\n",
            "ENST00000630155",
            &[][..],
            0,
            false,
        ),
        // A run of three `A` over the base before a reverse transcript and
        // its first two bases: `n.1del` over the corpus FASTA.
        (
            "deleted_over_the_first_base_in_a_run",
            "21\t13980158\trev_run\tAAAG\tAG,AAAGA\t.\t.\t.\n",
            "ENST00000766862",
            &[][..],
            1,
            true,
        ),
        // The last two bases and the base after them, beside a kept base
        // after them: the base before the three is not read.
        (
            "deleted_over_the_last_base_beside_a_base_not_read",
            "21\t13645822\tfwd_last_unread\tTCAA\tA,TCAAC\t.\t.\t.\n",
            "ENST00000299443",
            &[][..],
            0,
            true,
        ),
        // The same three bases beside a kept base before them, which is no
        // copy of the last of the three.
        (
            "deleted_over_the_last_base_beside_a_kept_base",
            "21\t13645820\tfwd_last_before\tCTTCA\tCT,CTTCAC\t.\t.\t.\n",
            "ENST00000299443",
            &[][..],
            0,
            false,
        ),
    ] {
        let input = write(
            dir.path(),
            &format!("{name}.vcf"),
            &format!("{VCF_HEADER}{record}"),
        );
        let filter = format!("stable_id is {transcript}");
        let over = |kind: &str, fasta: &Path| {
            let mut flags = vec![
                "--json_cache",
                json_cache.to_str().unwrap(),
                "--fasta",
                fasta.to_str().unwrap(),
                "--hgvs",
                "--transcript_filter",
                &filter,
                "--no_intergenic",
            ];
            flags.extend_from_slice(more);
            run_input(dir.path(), &format!("{name}_{kind}"), &input, &flags)
        };
        let held = over("held", &reference);
        assert!(!held.rows.is_empty(), "{name}");
        assert_eq!(notations_held(&held), (0, held_coding), "{name}");
        assert!(notices(&held.output).is_empty(), "{name}");
        let lacking = over("lacking", &other);
        assert_eq!(lacking.rows.len(), held.rows.len(), "{name}");
        assert_eq!(notations_held(&lacking), (0, 0), "{name}");
        let expected: &[&str] = if named { &[NO_SEQUENCE_FOR_21] } else { &[] };
        assert_eq!(notices(&lacking.output), expected, "{name}");
    }
}

/// Under `--process_ref_homs` the reference row of a record with one alternate
/// allele holds that allele's HGVSg (`OutputFactory.pm:1386-1394`), so it is
/// without one on a chromosome the FASTA holds no sequence for. `--pick` and
/// `--per_gene` keep that row alone of the record's rows, and the line names
/// HGVSg for it. Each run is held to what the same flags write over the corpus
/// FASTA, which holds the chromosome.
#[test]
fn a_reference_row_kept_alone_draws_the_line_for_the_hgvsg_it_is_without() {
    let dir = tempfile::tempdir().unwrap();
    let input = intron_substitution(dir.path());
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    let over = |name: &str, fasta: &Path, pick: &str| {
        run_input(
            dir.path(),
            name,
            &input,
            &[
                "--json_cache",
                json_cache.to_str().unwrap(),
                "--fasta",
                fasta.to_str().unwrap(),
                "--hgvsg",
                "--process_ref_homs",
                pick,
                "--tab",
                "--fields",
                "Allele,HGVSg",
            ],
        )
    };
    let record = std::fs::read_to_string(&input).unwrap();
    let reference_allele = record.lines().last().unwrap().split('\t').nth(3).unwrap();
    for (name, pick) in [("pick", "--pick"), ("per_gene", "--per_gene")] {
        let held = over(&format!("{name}_held"), &reference, pick);
        assert!(notices(&held.output).is_empty(), "{name}");
        let [row] = held.rows.as_slice() else {
            panic!("{name}: {:?}", held.rows)
        };
        let (allele, notation) = row.split_once('\t').unwrap();
        assert_eq!(allele, reference_allele, "{name}");
        assert!(notation.contains(":g."), "{name}: {row}");
        let lacking = over(&format!("{name}_lacking"), &other, pick);
        assert_eq!(lacking.rows, [format!("{reference_allele}\t-")], "{name}");
        assert_eq!(lacking.warnings, "", "{name}");
        assert_eq!(
            notices(&lacking.output),
            [NO_SEQUENCE_FOR_21_HGVSG],
            "{name}"
        );
    }
}

/// A run with a plugin annotates each batch whole and renders it in a pass
/// of its own, and reads its line for a chromosome the FASTA holds no
/// sequence for from the rows the record step leaves there too. Under
/// `--plugin NMD` the corpus's records draw the line for both notations,
/// their intergenic rows alone the line for HGVSg, and `--chr 22`, which
/// leaves no record, none. Each case is held to what the same flags write
/// over the corpus FASTA, which holds the chromosome.
#[test]
fn a_run_a_plugin_keeps_batched_names_the_chromosome_from_the_rows_the_output_keeps() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    let over = |name: &str, fasta: &Path, more: &[&str]| {
        let mut flags = vec![
            "--json_cache",
            json_cache.to_str().unwrap(),
            "--fasta",
            fasta.to_str().unwrap(),
            "--plugin",
            "NMD",
            "--hgvs",
            "--hgvsg",
        ];
        flags.extend_from_slice(more);
        run_records(dir.path(), name, &flags)
    };
    let cases: [KeptRowsCase; 3] = [
        ("every_row", &[], Some(NO_SEQUENCE_FOR_21_BOTH)),
        (
            "no_transcript_row",
            &[NO_TRANSCRIPT[0], NO_TRANSCRIPT[1]],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        ("chr_leaves_the_records_out", &["--chr", "22"], None),
    ];
    for (name, flags, line) in cases {
        let held = over(&format!("{name}_held"), &reference, flags);
        assert!(notices(&held.output).is_empty(), "{name}");
        let lacking = over(&format!("{name}_lacking"), &other, flags);
        assert_eq!(lacking.rows.len(), held.rows.len(), "{name}");
        let (held_genomic, held_coding) = notations_held(&held);
        let (genomic, coding) = notations_held(&lacking);
        let without = match line {
            None => (false, false),
            Some(NO_SEQUENCE_FOR_21_HGVSG) => (true, false),
            Some(_) => (true, true),
        };
        assert_eq!(
            (genomic < held_genomic, coding < held_coding),
            without,
            "{name}"
        );
        let expected: Vec<&str> = line.into_iter().collect();
        assert_eq!(notices(&lacking.output), expected, "{name}");
    }
}

/// The RefSeq corpus: a chromosome 21 slice of the RefSeq cache, in Ensembl
/// VEP's directory layout, and a reference FASTA over its records.
fn refseq_corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-refseq")
}

/// Under `--hgvs` an insertion or a deletion is shifted 3' along the genome
/// before its notation is written, also on a RefSeq model whose cached
/// sequence carries an `_rna_edit` and so differs from the genome; a shift
/// along the model's own sequence is made only on bases the genome lacks
/// (Ensembl VEP's own, `TranscriptVariationAllele.pm:211-256`, stands after
/// a return its `_return_3prime` takes at `:182-187` away from a
/// transcript's last and edited bases). A run on a RefSeq cache
/// says so once at its end, with the number of insertion and deletion alleles
/// whose written rows hold a notation on such a model: the corpus's
/// `--refseq --hgvs` run, which writes every row and
/// whose seven such alleles are those with an `HGVSc` on one of the 18 cached
/// models carrying an `_rna_edit` in the stored Ensembl VEP rows
/// (`expected_hgvs/default.txt.gz`), and a run on one of its insertions alone,
/// which lies on `NM_030582.4`. Under `--individual` a record is annotated
/// once for each individual that carries an allele of it
/// (`Parser/VCF.pm:627-662`), and an allele two individuals carry is counted
/// once: the count over the corpus's records, each carried by two
/// individuals, is the same seven, and a record given a second insertion
/// allele that one of the two carries alone counts two, or one for the
/// individual that carries the first alone. The
/// notice is withheld under `--quiet`, changes no row, and is not printed
/// under `--shift_hgvs 0`, without `--hgvs`, or on a cache of the Ensembl set.
#[test]
fn indels_on_edited_refseq_models_are_counted_once_under_hgvs() {
    let dir = tempfile::tempdir().unwrap();
    let cache_root = refseq_corpus().join("cache_root");
    let reference = refseq_corpus().join("reference.fa.gz");
    let refseq = [
        "--cache",
        "--dir_cache",
        cache_root.to_str().unwrap(),
        "--cache_version",
        "116",
        "--species",
        "homo_sapiens",
        "--refseq",
        "--fasta",
        reference.to_str().unwrap(),
    ];
    let with = |more: &[&'static str]| -> Vec<&str> {
        let mut flags = refseq.to_vec();
        flags.extend(more);
        flags
    };
    let variants = refseq_corpus().join("variants.vcf");
    let about_edited_models = |output: &Output| -> Vec<String> {
        notices(output)
            .into_iter()
            .filter(|line| line.contains("edited away from the genome"))
            .collect()
    };
    let tail = "an HGVS notation on a transcript whose sequence is edited away from the genome: the notation is shifted 3' along the genome, and along the transcript's own sequence only on bases the genome lacks";

    let loud = run_input(dir.path(), "loud", &variants, &with(&["--hgvs"]));
    assert_eq!(
        about_edited_models(&loud.output),
        [format!("INFO: 7 insertion or deletion alleles have {tail}")]
    );
    let quiet = run_input(
        dir.path(),
        "quiet",
        &variants,
        &with(&["--hgvs", "--quiet"]),
    );
    assert!(about_edited_models(&quiet.output).is_empty());
    assert_eq!(quiet.rows, loud.rows);
    assert_eq!(quiet.warnings, loud.warnings);

    let records = std::fs::read_to_string(&variants).unwrap();
    let insertion: String = records
        .lines()
        .filter(|line| line.starts_with('#') || line.split('\t').nth(2) == Some("rs749541171"))
        .map(|line| format!("{line}\n"))
        .collect();
    let one = write(dir.path(), "insertion.vcf", &insertion);
    let single = run_input(dir.path(), "single", &one, &with(&["--hgvs"]));
    assert!(
        single
            .rows
            .iter()
            .any(|row| row.contains("NM_030582.4") && row.contains("HGVSc=")),
        "{:?}",
        single.rows
    );
    assert_eq!(
        about_edited_models(&single.output),
        [format!("INFO: 1 insertion or deletion allele has {tail}")]
    );

    let genotyped: String = records
        .lines()
        .map(|line| match line {
            header if header.starts_with("##") => format!("{header}\n"),
            columns if columns.starts_with('#') => format!("{columns}\tFORMAT\tS1\tS2\n"),
            record => format!("{record}\tGT\t0/1\t1/1\n"),
        })
        .collect();
    let carried_twice = write(dir.path(), "genotyped.vcf", &genotyped);
    let copies = run_input(
        dir.path(),
        "copies",
        &carried_twice,
        &with(&["--hgvs", "--individual", "all"]),
    );
    for individual in ["IND=S1", "IND=S2"] {
        assert!(
            copies.rows.iter().any(|row| row.contains("NM_030582.4")
                && row.contains("HGVSc=")
                && row.contains(individual)),
            "{individual}"
        );
    }
    assert_eq!(
        about_edited_models(&copies.output),
        [format!("INFO: 7 insertion or deletion alleles have {tail}")]
    );
    let one_individual = run_input(
        dir.path(),
        "one_individual",
        &carried_twice,
        &with(&["--hgvs", "--individual", "S2"]),
    );
    assert_eq!(
        about_edited_models(&one_individual.output),
        [format!("INFO: 7 insertion or deletion alleles have {tail}")]
    );

    // The insertion's record with a second allele, its first with the last
    // base written once more: `S1` carries the first, `S2` both.
    let two_alleles: String = records
        .lines()
        .filter_map(|line| {
            if line.starts_with("##") {
                return Some(format!("{line}\n"));
            }
            if line.starts_with('#') {
                return Some(format!("{line}\tFORMAT\tS1\tS2\n"));
            }
            let mut columns: Vec<String> = line.split('\t').map(str::to_string).collect();
            if columns[2] != "rs749541171" {
                return None;
            }
            let first = columns[4].clone();
            columns[4] = format!("{first},{first}{}", &first[first.len() - 1..]);
            Some(format!("{}\tGT\t0/1\t1/2\n", columns.join("\t")))
        })
        .collect();
    let second_allele = write(dir.path(), "two_alleles.vcf", &two_alleles);
    for (individuals, counted) in [
        ("all", "2 insertion or deletion alleles have"),
        ("S1", "1 insertion or deletion allele has"),
        ("S2", "2 insertion or deletion alleles have"),
    ] {
        let run = run_input(
            dir.path(),
            &format!("two_alleles_{individuals}"),
            &second_allele,
            &with(&["--hgvs", "--individual", individuals]),
        );
        assert_eq!(
            about_edited_models(&run.output),
            [format!("INFO: {counted} {tail}")],
            "--individual {individuals}"
        );
    }

    for (name, flags) in [
        ("unshifted", with(&["--hgvs", "--shift_hgvs", "0"])),
        ("no_hgvs", with(&[])),
    ] {
        let run = run_input(dir.path(), name, &variants, &flags);
        assert!(
            about_edited_models(&run.output).is_empty(),
            "{name}: {}",
            text(&run.output.stderr)
        );
    }
    let ensembl_cache = corpus().join("json_cache");
    let ensembl_reference = corpus().join("reference.fa.gz");
    let ensembl_set = run_records(
        dir.path(),
        "ensembl_set",
        &[
            "--json_cache",
            ensembl_cache.to_str().unwrap(),
            "--fasta",
            ensembl_reference.to_str().unwrap(),
            "--hgvs",
        ],
    );
    assert!(
        notices(&ensembl_set.output).is_empty(),
        "{}",
        text(&ensembl_set.output.stderr)
    );
}

/// The line of standard error that counts `alleles` insertion or deletion
/// alleles with a notation on an edited model, whole.
fn edited_models_line(alleles: usize) -> String {
    let counted = if alleles == 1 {
        "1 insertion or deletion allele has".to_string()
    } else {
        format!("{alleles} insertion or deletion alleles have")
    };
    format!("INFO: {counted} an HGVS notation on a transcript whose sequence is edited away from the genome: the notation is shifted 3' along the genome, and along the transcript's own sequence only on bases the genome lacks")
}

/// The lines of a run's standard error about edited models.
fn about_edited_models(output: &Output) -> Vec<String> {
    notices(output)
        .into_iter()
        .filter(|line| line.contains("edited away from the genome"))
        .collect()
}

/// Runs the binary under `--refseq --hgvs` and `flags` over the RefSeq
/// corpus's cache and FASTA, writing `<name>.txt` under `dir`: on the
/// corpus's insertion `rs749541171` alone when `insertion_alone`, on every
/// record of the corpus otherwise.
fn refseq_hgvs_run(dir: &Path, name: &str, insertion_alone: bool, flags: &[&str]) -> Run {
    let cache_root = refseq_corpus().join("cache_root");
    let reference = refseq_corpus().join("reference.fa.gz");
    let variants = refseq_corpus().join("variants.vcf");
    let input = if insertion_alone {
        let insertion: String = std::fs::read_to_string(&variants)
            .unwrap()
            .lines()
            .filter(|line| line.starts_with('#') || line.split('\t').nth(2) == Some("rs749541171"))
            .map(|line| format!("{line}\n"))
            .collect();
        write(dir, "insertion.vcf", &insertion)
    } else {
        variants
    };
    let mut all = vec![
        "--cache",
        "--dir_cache",
        cache_root.to_str().unwrap(),
        "--cache_version",
        "116",
        "--species",
        "homo_sapiens",
        "--refseq",
        "--fasta",
        reference.to_str().unwrap(),
        "--hgvs",
    ];
    all.extend(flags);
    run_input(dir, name, &input, &all)
}

/// Whether a row of `run` holds an HGVSc on one of the two edited models
/// the insertion `rs749541171` overlaps, `NM_030582.4` and `NM_130444.3`:
/// every format writes the notation as the model's own id before `:c.`.
fn writes_an_edited_model_notation(run: &Run) -> bool {
    run.rows.iter().any(|row| {
        ["NM_030582.4:c.", "NM_130444.3:c."]
            .iter()
            .any(|notation| row.contains(notation))
    })
}

/// Holds a run under `flags` that writes no HGVSc on an edited model for the
/// insertion `rs749541171` to printing no line about edited models: the run
/// writes `rows` rows, `with_hgvsc` of them with a coding notation and none
/// on an edited model. Over the corpus's records the same flags print the
/// line for `corpus_alleles` alleles, or no line where that is `None`.
fn assert_counts_no_dropped_notation(
    flags: &[&str],
    rows: usize,
    with_hgvsc: usize,
    corpus_alleles: Option<usize>,
) {
    let dir = tempfile::tempdir().unwrap();
    let one = refseq_hgvs_run(dir.path(), "insertion", true, flags);
    assert_eq!(one.rows.len(), rows, "{:?}", one.rows);
    assert_eq!(
        one.rows.iter().filter(|row| row.contains(":c.")).count(),
        with_hgvsc,
        "{:?}",
        one.rows
    );
    assert!(!writes_an_edited_model_notation(&one), "{:?}", one.rows);
    assert!(
        about_edited_models(&one.output).is_empty(),
        "{}",
        text(&one.output.stderr)
    );
    let corpus = refseq_hgvs_run(dir.path(), "corpus", false, flags);
    assert_eq!(
        about_edited_models(&corpus.output),
        corpus_alleles
            .map(edited_models_line)
            .into_iter()
            .collect::<Vec<_>>()
    );
}

/// The count of the edited-model line is taken from the rows the output
/// keeps: an allele counts where a row the run writes holds an `HGVSc` on a
/// model edited away from the genome. The corpus's insertion `rs749541171`
/// gets an `HGVSc` on six RefSeq transcripts, two of them edited models.
/// `--coding_only` and `--flag_pick` keep a row on an edited model, as does
/// a `--pick` run whose `--transcript_filter` leaves `NM_030582.4` alone, a
/// tab or a VCF output whose `--fields` list names `HGVSc`, and the JSON
/// output: each prints the line for one allele, and writes the same rows
/// under `--quiet` without it. Over the corpus's records, of the seven
/// alleles the unfiltered run counts, four keep an `HGVSc` on one of the 18
/// cached models carrying an `_rna_edit` under `--coding_only`.
#[test]
fn indels_on_edited_models_are_counted_from_the_rows_the_output_keeps() {
    let dir = tempfile::tempdir().unwrap();
    for (name, flags, rows) in [
        ("all_rows", &[][..], 6),
        ("coding_only", &["--coding_only"][..], 3),
        ("flag_pick", &["--flag_pick"][..], 6),
        (
            "pick_of_an_edited_model",
            &["--pick", "--transcript_filter", "stable_id is NM_030582.4"][..],
            1,
        ),
        (
            "tab_fields_with_hgvsc",
            &["--tab", "--fields", "Uploaded_variation,Feature,HGVSc"][..],
            6,
        ),
        (
            "vcf_fields_with_hgvsc",
            &["--vcf", "--fields", "Allele,Feature,HGVSc"][..],
            1,
        ),
        ("json", &["--json"][..], 1),
    ] {
        let kept = refseq_hgvs_run(dir.path(), name, true, flags);
        assert_eq!(kept.rows.len(), rows, "{name}: {:?}", kept.rows);
        assert!(
            writes_an_edited_model_notation(&kept),
            "{name}: {:?}",
            kept.rows
        );
        assert_eq!(
            about_edited_models(&kept.output),
            [edited_models_line(1)],
            "{name}"
        );
        let quiet = refseq_hgvs_run(
            dir.path(),
            &format!("{name}_quiet"),
            true,
            &[flags, &["--quiet"][..]].concat(),
        );
        assert_eq!(quiet.rows, kept.rows, "{name}");
        assert!(about_edited_models(&quiet.output).is_empty(), "{name}");
    }
    let coding = refseq_hgvs_run(dir.path(), "corpus_coding_only", false, &["--coding_only"]);
    assert_eq!(about_edited_models(&coding.output), [edited_models_line(4)]);
}

/// `--pick` keeps one row of the insertion, on an unedited model, and the
/// run prints no line; over the corpus two alleles keep their picked row on
/// an edited model.
#[test]
fn pick_counts_the_indels_whose_kept_row_is_on_an_edited_model() {
    assert_counts_no_dropped_notation(&["--pick"], 1, 1, Some(2));
}

/// `--pick_allele` keeps one row of the insertion's one allele, on an
/// unedited model: no line, and two alleles over the corpus.
#[test]
fn pick_allele_counts_the_indels_whose_kept_row_is_on_an_edited_model() {
    assert_counts_no_dropped_notation(&["--pick_allele"], 1, 1, Some(2));
}

/// `--per_gene` keeps one row a gene of the insertion, both on unedited
/// models: no line, and two alleles over the corpus.
#[test]
fn per_gene_counts_the_indels_whose_kept_row_is_on_an_edited_model() {
    assert_counts_no_dropped_notation(&["--per_gene"], 2, 2, Some(2));
}

/// `--pick_allele_gene` keeps one row a gene of the insertion's one allele,
/// both on unedited models: no line, and two alleles over the corpus.
#[test]
fn pick_allele_gene_counts_the_indels_whose_kept_row_is_on_an_edited_model() {
    assert_counts_no_dropped_notation(&["--pick_allele_gene"], 2, 2, Some(2));
}

/// `--most_severe` writes one row a record with no transcript and no
/// notation (`OutputFactory.pm:426-459`): no line, for the insertion and
/// over the corpus.
#[test]
fn most_severe_writes_no_notation_and_counts_no_indel_on_an_edited_model() {
    assert_counts_no_dropped_notation(&["--most_severe"], 1, 0, None);
}

/// `--summary` writes one row a record with no transcript and no notation
/// (`OutputFactory.pm:426-459`): no line, for the insertion and over the
/// corpus.
#[test]
fn summary_writes_no_notation_and_counts_no_indel_on_an_edited_model() {
    assert_counts_no_dropped_notation(&["--summary"], 1, 0, None);
}

/// A tab or a VCF output whose `--fields` list names no `HGVSc` writes no
/// notation, and the run prints no line, for the insertion and over the
/// corpus.
#[test]
fn a_field_list_without_hgvsc_counts_no_indel_on_an_edited_model() {
    assert_counts_no_dropped_notation(
        &[
            "--tab",
            "--fields",
            "Uploaded_variation,Feature,Consequence",
        ],
        6,
        0,
        None,
    );
    assert_counts_no_dropped_notation(
        &["--vcf", "--fields", "Allele,Feature,Consequence"],
        1,
        0,
        None,
    );
}

/// A run with a plugin annotates each batch whole and renders it in a pass
/// of its own, and counts from the rows the record step leaves there too:
/// under `--plugin NMD` the insertion keeps its rows on the two edited models
/// and the run prints the line for one allele, seven over the corpus, and
/// under `--pick` beside it the insertion's kept row is on an unedited model
/// and the run prints no line, two alleles over the corpus.
#[test]
fn a_run_a_plugin_keeps_batched_counts_from_the_rows_the_output_keeps() {
    let dir = tempfile::tempdir().unwrap();
    let batched = ["--plugin", "NMD"];
    let kept = refseq_hgvs_run(dir.path(), "insertion", true, &batched);
    assert!(writes_an_edited_model_notation(&kept), "{:?}", kept.rows);
    assert_eq!(about_edited_models(&kept.output), [edited_models_line(1)]);
    let corpus = refseq_hgvs_run(dir.path(), "corpus", false, &batched);
    assert_eq!(about_edited_models(&corpus.output), [edited_models_line(7)]);
    assert_counts_no_dropped_notation(&["--plugin", "NMD", "--pick"], 1, 1, Some(2));
}

/// A record over `--max_sv_size` keeps its rows in the default and tab
/// outputs, while the VCF output carries its line without them and the JSON
/// output leaves the record out. The corpus's deletion `rs2037126010` lies on
/// the two edited models of the insertion `rs749541171` and spans 30 bases as
/// the size test measures it: every format prints the line for one allele,
/// and under `--max_sv_size 5` the default and tab outputs still write the
/// notation and print the line, where the VCF and JSON outputs write no
/// notation and the run prints no line.
#[test]
fn an_indel_over_max_sv_size_is_counted_where_the_output_writes_its_rows() {
    let dir = tempfile::tempdir().unwrap();
    let cache_root = refseq_corpus().join("cache_root");
    let reference = refseq_corpus().join("reference.fa.gz");
    let deletion: String = std::fs::read_to_string(refseq_corpus().join("variants.vcf"))
        .unwrap()
        .lines()
        .filter(|line| line.starts_with('#') || line.split('\t').nth(2) == Some("rs2037126010"))
        .map(|line| format!("{line}\n"))
        .collect();
    let input = write(dir.path(), "deletion.vcf", &deletion);
    let run = |name: &str, flags: &[&str]| {
        let mut all = vec![
            "--cache",
            "--dir_cache",
            cache_root.to_str().unwrap(),
            "--cache_version",
            "116",
            "--species",
            "homo_sapiens",
            "--refseq",
            "--fasta",
            reference.to_str().unwrap(),
            "--hgvs",
        ];
        all.extend(flags);
        run_input(dir.path(), name, &input, &all)
    };
    for (name, format, rows_written) in [
        ("default", &[][..], true),
        ("tab", &["--tab"][..], true),
        ("vcf", &["--vcf"][..], false),
        ("json", &["--json"][..], false),
    ] {
        let within = run(&format!("{name}_within"), format);
        assert!(
            writes_an_edited_model_notation(&within),
            "{name}: {:?}",
            within.rows
        );
        assert_eq!(
            about_edited_models(&within.output),
            [edited_models_line(1)],
            "{name}"
        );
        let over = run(
            &format!("{name}_over"),
            &[format, &["--max_sv_size", "5"][..]].concat(),
        );
        assert_eq!(
            writes_an_edited_model_notation(&over),
            rows_written,
            "{name}: {:?}",
            over.rows
        );
        let counted: Vec<String> = rows_written
            .then(|| edited_models_line(1))
            .into_iter()
            .collect();
        assert_eq!(about_edited_models(&over.output), counted, "{name}");
    }
}

/// The line a completed `vep recode` run writes last over a cache read as the
/// Ensembl transcript set.
const RECODE_ENSEMBL_SET: &str = "INFO: this cache is read as the Ensembl transcript set, and \
    RefSeq notation is written from a merged cache; no variant is mapped to an LRG, which needs \
    the Ensembl database";

/// The line a completed `vep recode` run writes last over a merged cache.
const RECODE_MERGED_SET: &str =
    "INFO: no variant is mapped to an LRG, which needs the Ensembl database";

/// The line a `vep recode` run whose fields name `id`, one of its default
/// fields, writes before its own over a cache that holds no variation data,
/// as the caches of these tests hold none; `cache` is the directory as the
/// command names it.
fn no_id_under(cache: &Path) -> String {
    format!(
        "INFO: no variation data under {}: id, read from a cache's known variants, stays \
         unwritten; a cache from vep-cache-builder --auto or --include-variations, or a \
         conversion of Ensembl VEP's cache, holds it (docs/cache-variation-data.md)",
        cache.display()
    )
}

/// The notice of a run under `--hgvs` alone whose FASTA holds no sequence for
/// chromosome 21.
const NO_SEQUENCE_FOR_21: &str = "INFO: The FASTA holds no sequence for chromosome 21, which the \
    input has records on: a substitution or a deletion there gets no HGVSc, nor does an \
    insertion written outside a transcript";

/// One `vep recode` run with `args`, under a home of the run's own: its exit
/// status, its whole standard error and its output.
fn recode(args: &[&str]) -> (Option<i32>, String, Vec<u8>) {
    let home = tempfile::tempdir().unwrap();
    let output = common::isolated_vep(home.path())
        .arg("recode")
        .args(args)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    (output.status.code(), text(&output.stderr), output.stdout)
}

/// [`recode`] for a run that must complete: its whole standard error and its
/// output.
fn recoded(args: &[&str]) -> (String, Vec<u8>) {
    let (code, stderr, stdout) = recode(args);
    assert_eq!(code, Some(0), "{args:?}: {stderr}");
    (stderr, stdout)
}

/// `lines`, each closed by a line feed: the whole standard error of a run.
fn lines_of(lines: &[&str]) -> String {
    lines.iter().map(|line| format!("{line}\n")).collect()
}

/// `vep recode` runs its annotation quiet, so the line of a FASTA file passed
/// over is the tool's to write: once a file, before the tool's own line, in a
/// run that completes, and the file changes no byte of the output. A cache
/// directory whose one FASTA lacks its index leaves the default fields, which
/// need a FASTA, refused with the error alone; `--fields id` needs none and
/// completes. A directory of indexed files draws no line of a file. Each run
/// that completes names `id` among its fields, so the tool's line of a cache
/// without variation data stands before its own.
#[test]
fn vep_recode_names_a_fasta_file_passed_over_before_its_own_line() {
    let dir = tempfile::tempdir().unwrap();
    let variants = corpus().join("variants.vcf");
    let variants = variants.to_str().unwrap();

    let cache = dir.path().join("cache");
    copy_dir(&corpus().join("json_cache"), &cache);
    let in_cache = ["-i", variants, "--json_cache", cache.to_str().unwrap()];
    let ids = [&in_cache[..], &["--fields", "id"][..]].concat();
    let no_id = no_id_under(&cache);
    let (stderr, expected_ids) = recoded(&ids);
    assert_eq!(stderr, lines_of(&[&no_id, RECODE_ENSEMBL_SET]));

    write_small_fasta(&cache.join("a.fa"), false);
    let a = cache.join("a.fa").display().to_string();
    let passed_over = |reads: &str| {
        format!(
            "INFO: FASTA file {a} in the cache directory was passed over (FASTA index {a}.fai \
             not found: index the FASTA with samtools faidx); {reads}"
        )
    };
    let (stderr, stdout) = recoded(&ids);
    assert_eq!(
        stderr,
        lines_of(&[
            &passed_over("the run has no reference FASTA"),
            &no_id,
            RECODE_ENSEMBL_SET
        ])
    );
    assert_eq!(stdout, expected_ids);
    let (code, stderr, stdout) = recode(&in_cache);
    assert_eq!(code, Some(1));
    assert_eq!(
        stderr,
        "ERROR: [fasta] Cannot generate HGVS coordinates (--hgvs and --hgvsg) in offline mode \
         without a FASTA file (see --fasta)\n"
    );
    assert!(stdout.is_empty());

    // A cache without its checksum manifest reads an indexed FASTA that the
    // manifest does not list; this one holds no sequence for the records'
    // chromosome, which the run says at its end.
    std::fs::remove_file(cache.join("CHECKSUMS.sha256")).unwrap();
    write_small_fasta(&cache.join("b.fa"), true);
    let reads_b = passed_over(&format!(
        "the run's reference FASTA is {}",
        cache.join("b.fa").display()
    ));
    let (stderr, beside) = recoded(&in_cache);
    assert_eq!(
        stderr,
        lines_of(&[
            &reads_b,
            NO_SEQUENCE_FOR_21_BOTH,
            &no_id,
            RECODE_ENSEMBL_SET
        ])
    );
    std::fs::remove_file(cache.join("a.fa")).unwrap();
    let (stderr, indexed) = recoded(&in_cache);
    assert_eq!(
        stderr,
        lines_of(&[NO_SEQUENCE_FOR_21_BOTH, &no_id, RECODE_ENSEMBL_SET])
    );
    assert_eq!(beside, indexed);

    let fastas = dir.path().join("fastas");
    std::fs::create_dir(&fastas).unwrap();
    for suffix in ["", ".fai", ".gzi"] {
        std::fs::copy(
            corpus().join(format!("reference.fa.gz{suffix}")),
            fastas.join(format!("reference.fa.gz{suffix}")),
        )
        .unwrap();
    }
    let json_cache = corpus().join("json_cache");
    let from_directory = [
        "-i",
        variants,
        "--json_cache",
        json_cache.to_str().unwrap(),
        "--fasta",
        fastas.to_str().unwrap(),
    ];
    let no_id = no_id_under(&json_cache);
    let (stderr, expected) = recoded(&from_directory);
    assert_eq!(stderr, lines_of(&[&no_id, RECODE_ENSEMBL_SET]));
    write_small_fasta(&fastas.join("plain.fa"), false);
    let plain = fastas.join("plain.fa").display().to_string();
    let (stderr, stdout) = recoded(&from_directory);
    assert_eq!(
        stderr,
        lines_of(&[
            &format!(
                "INFO: FASTA file {plain} in the --fasta directory was passed over (FASTA index \
                 {plain}.fai not found: index the FASTA with samtools faidx); its sequences are \
                 not read"
            ),
            &no_id,
            RECODE_ENSEMBL_SET
        ])
    );
    assert_eq!(stdout, expected);
}

/// The default fields of `vep recode` name `hgvsg` and `hgvsc`, so a FASTA
/// that holds no sequence for a chromosome of the input leaves every allele
/// there without the first, and a substitution or a deletion without the
/// second. The tool writes the run's notice of it once, before its own line,
/// and the line names the notations the fields ask for: both under the default
/// fields, `hgvsg` alone or `hgvsc` alone under either field. A FASTA that
/// holds the chromosome, or fields that ask for no HGVS notation, draw no
/// line of a chromosome. Where the fields name `id`, as the default fields
/// do, the tool's line of a cache without variation data stands before its
/// own.
#[test]
fn vep_recode_names_a_chromosome_the_fasta_lacks_before_its_own_line() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let variants = corpus().join("variants.vcf");
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    let over = |fasta: &Path, more: &[&'static str]| {
        let mut args = vec![
            "-i".to_string(),
            variants.display().to_string(),
            "--json_cache".to_string(),
            json_cache.display().to_string(),
            "--fasta".to_string(),
            fasta.display().to_string(),
        ];
        args.extend(more.iter().map(|flag| flag.to_string()));
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        recoded(&args)
    };

    let no_id = no_id_under(&json_cache);
    let (stderr, lacking) = over(&other, &[]);
    assert_eq!(
        stderr,
        lines_of(&[NO_SEQUENCE_FOR_21_BOTH, &no_id, RECODE_ENSEMBL_SET])
    );
    let (stderr, held) = over(&reference, &[]);
    assert_eq!(stderr, lines_of(&[&no_id, RECODE_ENSEMBL_SET]));
    let notations = |output: &[u8], key: &str| text(output).matches(key).count();
    assert!(
        notations(&lacking, "\"hgvsc\"") < notations(&held, "\"hgvsc\""),
        "the FASTA that lacks the chromosome leaves alleles without an hgvsc"
    );
    assert_eq!(notations(&lacking, "\"hgvsg\""), 0);
    assert!(notations(&held, "\"hgvsg\"") > 0);
    let (stderr, genomic) = over(&other, &["--fields", "hgvsg"]);
    assert_eq!(
        stderr,
        lines_of(&[NO_SEQUENCE_FOR_21_HGVSG, RECODE_ENSEMBL_SET])
    );
    assert_eq!(notations(&genomic, "\"hgvsg\""), 0);
    let (stderr, _) = over(&other, &["--fields", "hgvsc"]);
    assert_eq!(stderr, lines_of(&[NO_SEQUENCE_FOR_21, RECODE_ENSEMBL_SET]));
    let (stderr, _) = over(&other, &["--fields", "id"]);
    assert_eq!(stderr, lines_of(&[&no_id, RECODE_ENSEMBL_SET]));
}

/// `vep recode` writes an `hgvsg` or an `hgvsc` where its fields name the key
/// or `--mane_select` writes its triples, and its line for a chromosome the
/// FASTA holds no sequence for names those notations alone. Fields that need
/// the HGVS flag and write no `hgvsc` (`hgvsp`, and `hgvs`, which reads a key
/// no entry carries) draw the tool's line alone; beside `hgvsg` they draw the
/// line for HGVSg, beside `hgvsc` the line for HGVSc, and under
/// `--mane_select`, whose triples hold both, the line for both. Each case is
/// held to the keys its fields write over the corpus FASTA, which holds the
/// chromosome.
#[test]
fn vep_recode_names_a_chromosome_the_fasta_lacks_for_the_notations_it_writes() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let variants = corpus().join("variants.vcf");
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    let over = |fasta: &Path, more: &[&str]| {
        let mut args = vec![
            "-i".to_string(),
            variants.display().to_string(),
            "--json_cache".to_string(),
            json_cache.display().to_string(),
            "--fasta".to_string(),
            fasta.display().to_string(),
        ];
        args.extend(more.iter().map(|flag| flag.to_string()));
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        recoded(&args)
    };
    let writes = |output: &[u8], key: &str| text(output).contains(key);
    let cases: [(&[&str], Option<&str>); 7] = [
        (&["--fields", "hgvsp"], None),
        (&["--fields", "hgvs"], None),
        (&["--fields", "hgvsp,spdi"], None),
        (&["--fields", "hgvsp,hgvsg"], Some(NO_SEQUENCE_FOR_21_HGVSG)),
        (&["--fields", "hgvsp,hgvsc"], Some(NO_SEQUENCE_FOR_21)),
        (&["--fields", "hgvs,hgvsc"], Some(NO_SEQUENCE_FOR_21)),
        (
            &["--fields", "hgvsp", "--mane_select"],
            Some(NO_SEQUENCE_FOR_21_BOTH),
        ),
    ];
    for (fields, line) in cases {
        let (stderr, held) = over(&reference, fields);
        assert_eq!(stderr, lines_of(&[RECODE_ENSEMBL_SET]), "{fields:?}");
        let written = match line {
            None => (false, false),
            Some(NO_SEQUENCE_FOR_21_HGVSG) => (true, false),
            Some(NO_SEQUENCE_FOR_21) => (false, true),
            Some(_) => (true, true),
        };
        assert_eq!(
            (writes(&held, "\"hgvsg\""), writes(&held, "\"hgvsc\"")),
            written,
            "{fields:?}"
        );
        let (stderr, _) = over(&other, fields);
        let expected: Vec<&str> = line.into_iter().chain([RECODE_ENSEMBL_SET]).collect();
        assert_eq!(stderr, lines_of(&expected), "{fields:?}");
    }
}

/// `vep recode` names a chromosome the FASTA holds no sequence for, and a
/// notation, only where an entry it writes there is without the notation for
/// want of the sequence. Under a `--transcript_filter` no transcript passes
/// every entry is an intergenic one, which holds an `hgvsg` and never an
/// `hgvsc`: fields that write `hgvsg` draw the line for HGVSg alone, and
/// fields that write `hgvsc` alone, or `--mane_select` beside `hgvsp`, whose
/// triples are those of MANE Select entries, draw the tool's line alone. Over
/// records that are insertions alone, each `hgvsc` stands over either FASTA
/// and the line names HGVSg alone. Each case is held to what the same
/// arguments write over the corpus FASTA, which holds the chromosome: the
/// line names a notation exactly where that output holds more of it. A case
/// that names no `--fields` reads the default fields, `id` among them, and
/// ends on the tool's line of a cache without variation data and its own.
#[test]
fn vep_recode_names_a_chromosome_the_fasta_lacks_where_an_entry_is_without_the_notation() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let variants = corpus().join("variants.vcf");
    let (insertions, _) = insertions_alone(dir.path());
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    let over = |input: &Path, fasta: &Path, more: &[&str]| {
        let mut args = vec![
            "-i".to_string(),
            input.display().to_string(),
            "--json_cache".to_string(),
            json_cache.display().to_string(),
            "--fasta".to_string(),
            fasta.display().to_string(),
        ];
        args.extend(more.iter().map(|flag| flag.to_string()));
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        recoded(&args)
    };
    let held_in = |output: &[u8]| {
        let output = text(output);
        let held = |mark: &str| output.matches(mark).count();
        (held(":g."), held(":c.") + held(":n."))
    };
    let cases: [(&Path, &[&str], Option<&str>); 8] = [
        (
            &variants,
            &[NO_TRANSCRIPT[0], NO_TRANSCRIPT[1]],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        (
            &variants,
            &["--fields", "hgvsg", NO_TRANSCRIPT[0], NO_TRANSCRIPT[1]],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        (
            &variants,
            &["--fields", "hgvsc", NO_TRANSCRIPT[0], NO_TRANSCRIPT[1]],
            None,
        ),
        (
            &variants,
            &[
                "--fields",
                "hgvsp",
                "--mane_select",
                NO_TRANSCRIPT[0],
                NO_TRANSCRIPT[1],
            ],
            None,
        ),
        (
            &variants,
            &[
                "--fields",
                "hgvsg",
                "--mane_select",
                NO_TRANSCRIPT[0],
                NO_TRANSCRIPT[1],
            ],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        (&insertions, &["--fields", "hgvsc"], None),
        (
            &insertions,
            &["--fields", "hgvsc,hgvsg"],
            Some(NO_SEQUENCE_FOR_21_HGVSG),
        ),
        (&insertions, &[], Some(NO_SEQUENCE_FOR_21_HGVSG)),
    ];
    let no_id = no_id_under(&json_cache);
    for (input, more, line) in cases {
        let own: Vec<&str> = (!more.contains(&"--fields"))
            .then_some(no_id.as_str())
            .into_iter()
            .chain([RECODE_ENSEMBL_SET])
            .collect();
        let (stderr, held) = over(input, &reference, more);
        assert_eq!(stderr, lines_of(&own), "{more:?}");
        let (stderr, lacking) = over(input, &other, more);
        let (held_genomic, held_coding) = held_in(&held);
        let (genomic, coding) = held_in(&lacking);
        let without = match line {
            None => (false, false),
            Some(NO_SEQUENCE_FOR_21_HGVSG) => (true, false),
            Some(NO_SEQUENCE_FOR_21) => (false, true),
            Some(_) => (true, true),
        };
        assert_eq!(
            (genomic < held_genomic, coding < held_coding),
            without,
            "{more:?}"
        );
        let expected: Vec<&str> = line.into_iter().chain(own).collect();
        assert_eq!(stderr, lines_of(&expected), "{more:?}");
    }
}

/// `vep recode` writes the `hgvsc` of an allele that keeps bases of its
/// reference and deletes others over a transcript's first base where the
/// deleted bases can be written wholly inside the transcript, which turns on
/// a base beside the record. On a chromosome the FASTA holds no sequence for
/// the entry is without it, and the line names HGVSc before the tool's own
/// line; the same deletion beside a kept base that settles it has no `hgvsc`
/// over either FASTA and draws no line. Each run keeps one transcript's
/// entries.
#[test]
fn vep_recode_names_a_chromosome_where_an_entry_turns_on_bases_beside_its_record() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    write_small_fasta(&other, true);
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    for (name, record, transcript, held_coding, line) in [
        (
            "deleted_over_the_first_base_in_a_run",
            "21\t13980158\trev_run\tAAAG\tAG,AAAGA\t.\t.\t.\n",
            "ENST00000766862",
            1,
            Some(NO_SEQUENCE_FOR_21),
        ),
        (
            "deleted_over_the_last_base_beside_a_kept_base",
            "21\t13645820\tfwd_last_before\tCTTCA\tCT,CTTCAC\t.\t.\t.\n",
            "ENST00000299443",
            0,
            None,
        ),
    ] {
        let input = write(
            dir.path(),
            &format!("{name}.vcf"),
            &format!("{VCF_HEADER}{record}"),
        );
        let filter = format!("stable_id is {transcript}");
        let over = |fasta: &Path| {
            recoded(&[
                "-i",
                input.to_str().unwrap(),
                "--json_cache",
                json_cache.to_str().unwrap(),
                "--fasta",
                fasta.to_str().unwrap(),
                "--fields",
                "hgvsc",
                "--transcript_filter",
                &filter,
            ])
        };
        let coding = |output: &[u8]| {
            let output = text(output);
            output.matches(":c.").count() + output.matches(":n.").count()
        };
        let (stderr, held) = over(&reference);
        assert_eq!(stderr, lines_of(&[RECODE_ENSEMBL_SET]), "{name}");
        assert_eq!(coding(&held), held_coding, "{name}");
        let (stderr, lacking) = over(&other);
        assert_eq!(coding(&lacking), 0, "{name}");
        let expected: Vec<&str> = line.into_iter().chain([RECODE_ENSEMBL_SET]).collect();
        assert_eq!(stderr, lines_of(&expected), "{name}");
    }
}

/// The merged corpus: a chromosome 21 slice of the merged cache, converted,
/// and a reference FASTA over its records.
fn merged_corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-merged")
}

/// `vep recode` reads a merged cache, whose RefSeq models can carry an
/// `_rna_edit`, and its default fields name `hgvsc`: the tool writes the
/// run's count of insertion and deletion alleles with a notation on such a
/// model once, before its own line. Over the merged corpus the count is 12,
/// the alleles whose rows carry an `HGVSc` on one of the slice's 28 models
/// with an `_rna_edit`, and the line is the one `vep --merged --hgvs` prints
/// over the same cache, FASTA and records. `--shift_hgvs 0`, which shifts no
/// notation, and fields that ask for no HGVS notation draw no such count.
/// The fields of each run name `id` and the cache holds no variation data,
/// so the tool's line of the unwritten `id` stands before its own.
#[test]
fn vep_recode_counts_the_indels_on_edited_models_before_its_own_line() {
    let dir = tempfile::tempdir().unwrap();
    let variants = merged_corpus().join("variants.vcf");
    let json_cache = merged_corpus().join("json_cache");
    let reference = merged_corpus().join("reference.fa.gz");
    let resources = [
        "--json_cache",
        json_cache.to_str().unwrap(),
        "--fasta",
        reference.to_str().unwrap(),
    ];
    let counted = "INFO: 12 insertion or deletion alleles have an HGVS notation on a transcript \
        whose sequence is edited away from the genome: the notation is shifted 3' along the \
        genome, and along the transcript's own sequence only on bases the genome lacks";

    let annotated = run_input(
        dir.path(),
        "annotated",
        &variants,
        &[&resources[..], &["--merged", "--hgvs"][..]].concat(),
    );
    assert_eq!(notices(&annotated.output), [counted]);

    let input = ["-i", variants.to_str().unwrap()];
    let recode_with = |more: &[&str]| recoded(&[&input[..], &resources[..], more].concat());
    let no_id = no_id_under(&json_cache);
    let (stderr, _) = recode_with(&[]);
    assert_eq!(stderr, lines_of(&[counted, &no_id, RECODE_MERGED_SET]));
    let (stderr, _) = recode_with(&["--shift_hgvs", "0"]);
    assert_eq!(stderr, lines_of(&[&no_id, RECODE_MERGED_SET]));
    let (stderr, _) = recode_with(&["--fields", "id"]);
    assert_eq!(stderr, lines_of(&[&no_id, RECODE_MERGED_SET]));
}

/// The transcripts whose `hgvsc` a `vep recode` output writes, by the id that
/// opens each notation: under the `hgvsc` key of an allele's object, and in
/// its `mane_select` triples.
fn written_hgvsc_ids(output: &[u8]) -> (Vec<String>, Vec<String>) {
    let id_of = |notation: &serde_json::Value| -> Option<String> {
        notation
            .as_str()
            .and_then(|notation| notation.split_once(':'))
            .map(|(id, _)| id.to_string())
    };
    let array: serde_json::Value = serde_json::from_slice(output).unwrap();
    let (mut under_the_key, mut in_triples) = (Vec::new(), Vec::new());
    let alleles = array
        .as_array()
        .unwrap()
        .iter()
        .filter_map(serde_json::Value::as_object)
        .flat_map(|object| object.values());
    for allele in alleles {
        let list = |key: &str| allele[key].as_array().cloned().unwrap_or_default();
        under_the_key.extend(list("hgvsc").iter().filter_map(id_of));
        in_triples.extend(
            list("mane_select")
                .iter()
                .filter_map(|triple| id_of(&triple["hgvsc"])),
        );
    }
    under_the_key.sort();
    in_triples.sort();
    (under_the_key, in_triples)
}

/// `vep recode` counts the insertion and deletion alleles whose notation on
/// an edited model it writes: the `hgvsc` of every consequence entry where
/// the fields name `hgvsc`, and of the MANE Select entries alone under
/// `--mane_select` without that field. The deletion `rs754050130` of the
/// merged corpus overlaps three models with an `_rna_edit`, none of them
/// MANE Select. `--fields hgvsc` writes their notations and draws the line
/// for one allele, with `--mane_select` as without. `--fields hgvsp` and
/// `--fields hgvs`, a name no consequence entry has a key of, write no
/// `hgvsc`, and `--mane_select --fields hgvsp` writes the notations of the
/// two MANE Select entries, on unedited models: each draws the tool's line
/// alone. Over the corpus's records, fields that write no `hgvsc` draw no
/// count, and neither does `--mane_select --fields hgvsp`, since no row of
/// the annotation run is on a MANE Select model with a `BAM_EDIT`. Fields
/// that write no `hgvsc` and name `id`, over either input, draw the tool's
/// line of a cache without variation data and then its own; the others
/// draw the tool's line alone.
#[test]
fn vep_recode_counts_the_indels_whose_notation_on_an_edited_model_it_writes() {
    const EDITED: [&str; 3] = ["NM_001001713.1", "NM_001317741.1", "NM_001317742.1"];
    const MANE_SELECT: [&str; 2] = ["ENST00000333634.10", "NM_007341.3"];
    let dir = tempfile::tempdir().unwrap();
    let variants = merged_corpus().join("variants.vcf");
    let deletion: String = std::fs::read_to_string(&variants)
        .unwrap()
        .lines()
        .filter(|line| line.starts_with('#') || line.split('\t').nth(2) == Some("rs754050130"))
        .map(|line| format!("{line}\n"))
        .collect();
    let deletion = write(dir.path(), "deletion.vcf", &deletion);
    let json_cache = merged_corpus().join("json_cache");
    let reference = merged_corpus().join("reference.fa.gz");
    let resources = [
        "--json_cache",
        json_cache.to_str().unwrap(),
        "--fasta",
        reference.to_str().unwrap(),
    ];
    let over = |input: &Path, more: &[&str]| {
        let (stderr, output) =
            recoded(&[&["-i", input.to_str().unwrap()][..], &resources[..], more].concat());
        (stderr, written_hgvsc_ids(&output))
    };
    let on_edited_models = |ids: &[String]| -> Vec<String> {
        ids.iter()
            .filter(|id| EDITED.contains(&id.as_str()))
            .cloned()
            .collect()
    };
    let counted_once = lines_of(&[&edited_models_line(1), RECODE_MERGED_SET]);
    let own_line_alone = lines_of(&[RECODE_MERGED_SET]);
    let own_line_after_no_id = lines_of(&[&no_id_under(&json_cache), RECODE_MERGED_SET]);
    let uncounted = |fields: &str| {
        if fields.split(',').any(|field| field == "id") {
            &own_line_after_no_id
        } else {
            &own_line_alone
        }
    };

    let (stderr, (under_the_key, in_triples)) = over(&deletion, &["--fields", "hgvsc"]);
    assert_eq!(on_edited_models(&under_the_key), EDITED);
    assert!(in_triples.is_empty());
    assert_eq!(stderr, counted_once);
    let (stderr, (under_the_key, in_triples)) =
        over(&deletion, &["--mane_select", "--fields", "hgvsc"]);
    assert_eq!(on_edited_models(&under_the_key), EDITED);
    assert_eq!(in_triples, MANE_SELECT);
    assert_eq!(stderr, counted_once);

    for fields in ["hgvsp", "hgvs", "id,hgvsp,spdi"] {
        let (stderr, (under_the_key, in_triples)) = over(&deletion, &["--fields", fields]);
        assert!(
            under_the_key.is_empty() && in_triples.is_empty(),
            "{fields}"
        );
        assert_eq!(&stderr, uncounted(fields), "{fields}");
    }
    let (stderr, (under_the_key, in_triples)) =
        over(&deletion, &["--mane_select", "--fields", "hgvsp"]);
    assert!(under_the_key.is_empty());
    assert_eq!(in_triples, MANE_SELECT);
    assert_eq!(stderr, own_line_alone);

    for fields in ["hgvsp", "hgvs", "id,hgvsp,spdi"] {
        let (stderr, (under_the_key, in_triples)) = over(&variants, &["--fields", fields]);
        assert!(
            under_the_key.is_empty() && in_triples.is_empty(),
            "{fields}"
        );
        assert_eq!(&stderr, uncounted(fields), "{fields}");
    }
    // Each row of the annotation run by its model, its MANE Select partner
    // and its edit status: no MANE Select row is on an edited model.
    let rows = run_input(
        dir.path(),
        "rows",
        &variants,
        &[
            &resources[..],
            &[
                "--merged",
                "--hgvs",
                "--mane_select",
                "--tab",
                "--fields",
                "Feature,MANE_SELECT,BAM_EDIT",
            ][..],
        ]
        .concat(),
    )
    .rows;
    let cells = |row: &str| -> Vec<String> { row.split('\t').map(str::to_string).collect() };
    assert!(rows.iter().any(|row| cells(row)[1] != "-"));
    assert!(rows.iter().any(|row| cells(row)[2] == "OK"));
    assert!(!rows
        .iter()
        .any(|row| cells(row)[1] != "-" && cells(row)[2] != "-"));
    let (stderr, (under_the_key, in_triples)) =
        over(&variants, &["--mane_select", "--fields", "hgvsp"]);
    assert!(under_the_key.is_empty() && !in_triples.is_empty());
    assert_eq!(stderr, own_line_alone);
}

/// The lines a `vep recode` run keeps stand in the order `vep` prints them:
/// the FASTA file passed over, which a run names at its start; the line of a
/// `--transcript_filter` on a key a cache built from a GTF leaves undefined;
/// the chromosome the FASTA lacks, which a run names at its end; then the
/// tool's line of a cache without variation data, the default fields naming
/// `id`, and its own line. No transcript passes the filter, so every entry
/// the run writes is an intergenic one, which holds an `hgvsg` and no `hgvsc`
/// over any FASTA: the line for the chromosome names HGVSg alone.
#[test]
fn vep_recode_writes_the_lines_its_run_kept_in_the_order_vep_prints_them() {
    let dir = tempfile::tempdir().unwrap();
    let built =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-gtf-builder");
    let cache = dir.path().join("cache");
    copy_dir(&built.join("json_cache"), &cache);
    std::fs::remove_file(cache.join("CHECKSUMS.sha256")).unwrap();
    write_small_fasta(&cache.join("a.fa"), false);
    write_small_fasta(&cache.join("b.fa"), true);
    let variants = built.join("variants.vcf");
    let (stderr, _) = recoded(&[
        "-i",
        variants.to_str().unwrap(),
        "--json_cache",
        cache.to_str().unwrap(),
        "--transcript_filter",
        "dbID exists",
    ]);
    let a = cache.join("a.fa").display().to_string();
    assert_eq!(
        stderr,
        lines_of(&[
            &format!(
                "INFO: FASTA file {a} in the cache directory was passed over (FASTA index \
                 {a}.fai not found: index the FASTA with samtools faidx); the run's reference FASTA \
                 is {}",
                cache.join("b.fa").display()
            ),
            "INFO: --transcript_filter names dbID, which no transcript of this cache carries: \
             a cache built from a GTF or GFF3 leaves it undefined, so exists and a comparison on \
             it fail for each of its transcripts",
            NO_SEQUENCE_FOR_21_HGVSG,
            &no_id_under(&cache),
            RECODE_ENSEMBL_SET,
        ])
    );
}

/// The two warnings of an input that leaves part of itself out stand where
/// `vep recode` carries every warning of its annotation run, in the `warnings`
/// member of the output: a `#` line of a list input read as a comment, and a
/// VCF record whose ALT is a lone `.` under a true INFO/SVTYPE. Standard
/// error holds the tool's line alone.
#[test]
fn vep_recode_carries_the_two_input_warnings_in_its_output() {
    let dir = tempfile::tempdir().unwrap();
    let json_cache = corpus().join("json_cache");
    let reference = corpus().join("reference.fa.gz");
    let warnings_of = |input: &Path| -> serde_json::Value {
        let (stderr, stdout) = recoded(&[
            "-i",
            input.to_str().unwrap(),
            "--json_cache",
            json_cache.to_str().unwrap(),
            "--fasta",
            reference.to_str().unwrap(),
            "--fields",
            "hgvsg",
        ]);
        assert_eq!(stderr, lines_of(&[RECODE_ENSEMBL_SET]));
        let objects: serde_json::Value = serde_json::from_slice(&stdout).unwrap();
        assert_eq!(objects.as_array().map(Vec::len), Some(1), "{objects}");
        objects[0]["warnings"].clone()
    };

    let list = write(
        dir.path(),
        "list.txt",
        "#a comment line\n21 10538885 10538885 A/C +\n",
    );
    assert_eq!(
        warnings_of(&list),
        serde_json::json!(["Read 1 line that begins with '#' as a comment, at line 1"])
    );

    let records = write(
        dir.path(),
        "records.vcf",
        &format!(
            "{VCF_HEADER}21\t10538885\tabsent\tA\t.\t.\t.\tSVTYPE=DEL;END=10538900\n\
             21\t10538885\tgiven\tA\tC\t.\t.\t.\n"
        ),
    );
    assert_eq!(
        warnings_of(&records),
        serde_json::json!([
            "line 1 (21 10538885 absent A . . . SVTYPE=DEL;END=10538900): ALT is '.', which \
             names no allele, so the record is read as a non-variant record and its INFO/SVTYPE \
             is not annotated"
        ])
    );
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A VCF record whose ALT column is a lone `.` names no allele and is a
/// non-variant record: dropped, or carried without consequences under
/// `--allow_non_variant`. Under a true INFO/SVTYPE Ensembl VEP tests the
/// structural shape first (Parser/VCF.pm:237-245, before the missing ALT of
/// `:259`) and annotates the record from INFO alone; here it stays a
/// non-variant record, and the run warns about it once by its line, with and
/// without the flag. A lone `.` without `SVTYPE`, or under one Perl reads as
/// false (`SVTYPE=0`), draws no warning, as in Ensembl VEP, and the rows are
/// those of the input without the non-variant records.
#[test]
fn a_missing_alt_under_svtype_is_warned_about_by_its_line() {
    let dir = tempfile::tempdir().unwrap();
    let first = "21\t10538885\tfirst\tA\tC\t.\t.\t.\n";
    let second = "21\t14221107\tsecond\tC\tA\t.\t.\t.\n";
    let plain = write(
        dir.path(),
        "plain.vcf",
        &format!("{VCF_HEADER}{first}{second}"),
    );
    let dotted = write(
        dir.path(),
        "dotted.vcf",
        &format!(
            "{VCF_HEADER}{first}\
             21\t10538900\tdel\tN\t.\t.\t.\tSVTYPE=DEL;END=10539200\n\
             21\t10538950\tplain\tA\t.\t.\t.\t.\n\
             21\t10538960\tfalse\tA\t.\t.\t.\tSVTYPE=0\n\
             21\t10538970\tbnd\tN\t.\t.\t.\tSVTYPE=BND\n\
             {second}"
        ),
    );
    let warned = [
        "WARNING: line 2 (21 10538900 del N . . . SVTYPE=DEL;END=10539200): ALT is '.', which names no allele, so the record is read as a non-variant record and its INFO/SVTYPE is not annotated",
        "WARNING: line 5 (21 10538970 bnd N . . . SVTYPE=BND): ALT is '.', which names no allele, so the record is read as a non-variant record and its INFO/SVTYPE is not annotated",
    ];
    for flags in [&["--quiet"][..], &["--quiet", "--allow_non_variant"][..]] {
        let name = format!("flags{}", flags.len());
        let expected = annotate(dir.path(), &format!("plain_{name}"), &plain, flags);
        assert!(!expected.rows.is_empty());
        assert_eq!(expected.warnings, "", "{flags:?}");
        let run = annotate(dir.path(), &format!("dotted_{name}"), &dotted, flags);
        assert_eq!(run.rows, expected.rows, "{flags:?}");
        assert_eq!(run.warning_lines(), warned, "{flags:?}");
    }
    let loud = annotate(dir.path(), "loud", &dotted, &[]);
    assert_eq!(loud.warning_lines(), warned);
    assert_eq!(loud.warned_on_standard_error("ALT is '.'"), warned);
}

/// The warning for a lone `.` under a true INFO/SVTYPE is drawn once a record
/// whatever the kept record writes: the reference allele's rows under
/// `--process_ref_homs`, a copy per individual under `--individual`, which
/// switches `--allow_non_variant` on, the one record under
/// `--individual_zyg`, its reference row where the one individual named is
/// reference-homozygous, and nothing where no listed individual has a
/// genotype on the record. A lone `.` without `SVTYPE` draws none under any
/// of them.
#[test]
fn a_missing_alt_under_svtype_is_warned_about_whatever_the_record_writes() {
    let dir = tempfile::tempdir().unwrap();
    let genotyped = write(
        dir.path(),
        "genotyped.vcf",
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
         21\t10538885\tfirst\tA\tC\t.\t.\t.\tGT\t0/1\t0/1\n\
         21\t10538900\tnone\tN\t.\t.\t.\tSVTYPE=DEL\tGT\t./.\t./.\n\
         21\t10538910\tone\tN\t.\t.\t.\tSVTYPE=DEL\tGT\t0/0\t./.\n\
         21\t10538920\tboth\tN\t.\t.\t.\tSVTYPE=DEL\tGT\t0/0\t0/0\n\
         21\t10538950\tplain\tA\t.\t.\t.\t.\tGT\t0/0\t0/0\n\
         21\t14221107\tsecond\tC\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
    );
    let warned = [
        "WARNING: line 2 (21 10538900 none N . . . SVTYPE=DEL GT ./. ./.): ALT is '.', which names no allele, so the record is read as a non-variant record and its INFO/SVTYPE is not annotated",
        "WARNING: line 3 (21 10538910 one N . . . SVTYPE=DEL GT 0/0 ./.): ALT is '.', which names no allele, so the record is read as a non-variant record and its INFO/SVTYPE is not annotated",
        "WARNING: line 4 (21 10538920 both N . . . SVTYPE=DEL GT 0/0 0/0): ALT is '.', which names no allele, so the record is read as a non-variant record and its INFO/SVTYPE is not annotated",
    ];
    for (name, flags) in [
        ("dropped", &["--quiet"][..]),
        ("kept", &["--quiet", "--allow_non_variant"][..]),
        (
            "reference_rows",
            &["--quiet", "--allow_non_variant", "--process_ref_homs"][..],
        ),
        ("copies", &["--quiet", "--individual", "all"][..]),
        (
            "copies_reference_rows",
            &["--quiet", "--individual", "all", "--process_ref_homs"][..],
        ),
        ("zygosities", &["--quiet", "--individual_zyg", "all"][..]),
        ("one_zygosity", &["--quiet", "--individual_zyg", "S1"][..]),
    ] {
        let run = annotate(dir.path(), name, &genotyped, flags);
        assert!(!run.rows.is_empty(), "{flags:?}");
        assert_eq!(run.warning_lines(), warned, "{flags:?}");
    }
}

/// The page under `docs/` at `path`, its white space collapsed.
fn page(path: &str) -> String {
    std::fs::read_to_string(format!("{}/../../docs/{path}", env!("CARGO_MANIFEST_DIR")))
        .unwrap()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

/// Each page that states what a run leaves out of what it was given to read
/// states the line the run prints for it, in the words the tests above read
/// from the binary.
#[test]
fn the_pages_state_the_line_each_run_prints() {
    for (path, phrases) in [
        (
            "cache-setup.md",
            &[
                "`INFO: FASTA file <file> in the cache directory was passed over (<what the file lacks and what mends it>); the run has no reference FASTA`, ending `; the run's reference FASTA is <file>` when another file there is indexed",
                "`INFO: FASTA file <file> in the --fasta directory was passed over (<the same>); its sequences are not read`",
                "`INFO: The FASTA holds no sequence for chromosome <name>, which the input has records on: a substitution or a deletion there gets no HGVSc, nor does an insertion written outside a transcript`",
            ][..],
        ),
        (
            "output/refseq-and-merged.md",
            &[
                "An insertion or a deletion on such a model is shifted 3' along the genome before its notation is written, and along the model's own sequence only on bases the genome lacks",
                "A run that writes such a notation counts those alleles in one `INFO:` line at its end, withheld under `--quiet`",
                "Ensembl VEP's source holds a second shift along the edited sequence (`TranscriptVariationAllele.pm` lines 211-256)",
                "Away from a transcript's last bases and its edited bases it returns before that shift (lines 182-187), copying a shift it already has",
                "Its GenBank protein name (lines 1750-1756) needs the transcript's analysis record, which a cache does not carry, so both engines name the cached protein identifier",
                "An LRG name takes no version in either engine",
            ][..],
        ),
        (
            "divergences/defect-j.md",
            &[
                "`WARNING: line <n> (<the record's line>): ALT is '.', which names no allele, so the record is read as a non-variant record and its INFO/SVTYPE is not annotated`",
            ][..],
        ),
    ] {
        let text = page(path);
        for phrase in phrases {
            assert!(text.contains(phrase), "docs/{path} lacks: {phrase}");
        }
    }
}
