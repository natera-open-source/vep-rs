// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Flags vep-rs declares and refuses, through the `vep` binary: each exits 1
//! with exactly one line `ERROR: [cli] <message>` on standard error, before any
//! input is read. The messages are the error table's (`vep_cli::config`). The
//! inputs refused by content (BCF, a CAid line) exit 1 the same way under the
//! `input` component, before any record is annotated.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Stdio;

use vep_cli::config::{error_table_entry, Treatment, ERROR_TABLE, GFF_GTF_RECIPE};

mod common;

/// Runs the binary with `flags` ahead of an input that does not exist, so a flag
/// that were accepted would fail later, on the input, under another component.
fn refusal(flags: &[&str]) -> (Option<i32>, String) {
    let output = common::isolated_vep_anywhere()
        .args(flags)
        .args(["-i", "/nonexistent/input.vcf", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

fn assert_refused(flags: &[&str], message: &str) {
    let (code, stderr) = refusal(flags);
    assert_eq!(code, Some(1), "{flags:?}: {stderr}");
    assert_eq!(stderr, format!("ERROR: [cli] {message}\n"), "{flags:?}");
}

/// Ensembl VEP resolves an ambiguous HGVSp to every variant through its
/// database (Parser/HGVS.pm:178-186, 224); vep-rs resolves a genomic,
/// transcript or protein notation to one variant and skips an ambiguous
/// protein notation with Ensembl VEP's message (VariationFeatureAdaptor.pm:1956).
#[test]
fn ambiguous_hgvs_is_refused() {
    assert_refused(
        &["--ambiguous_hgvs"],
        "--ambiguous_hgvs is not supported: it resolves an HGVSp input to every matching variant through the Ensembl database; HGVS input is genomic (g., m.), transcript (c., n.) and protein (p.) notation, each resolved to one variant, a protein change with several equally short nucleotide paths skipped with Ensembl VEP's message",
    );
}

/// A bare `--` is an extra command-line parameter: Getopt::Long runs with
/// `pass_through` (Config.pm:62), which leaves it and every later token in
/// `@ARGV`, and vep dies on what is left (vep:30, `ERROR: Unexpected extra
/// command-line parameter(s): --`, the tokens joined by `, `), unless `--help`
/// was read ahead of the `--` (vep:29); a `--help` behind it is one more extra
/// parameter. An empty value, `--pick_order ''`, is not an extra parameter
/// either: it is the empty order (Config.pm:761 splits it to nothing) and the
/// run goes on to its input, while a middle empty field, `canonical,,rank`, is
/// a value `%VALID` lacks (:946), clap's syntax error here, exit 2.
#[test]
fn a_bare_double_dash_is_an_extra_parameter() {
    let run = |argv: &[&str]| {
        let output = common::isolated_vep_anywhere()
            .args(argv)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        (
            output.status.code(),
            String::from_utf8_lossy(&output.stderr).into_owned(),
        )
    };
    let refused = |argv: &[&str], extra: &str| {
        let (code, stderr) = run(argv);
        assert_eq!(code, Some(1), "{argv:?}: {stderr}");
        assert_eq!(
            stderr,
            format!("ERROR: [cli] Unexpected extra command-line parameter(s): {extra}\n"),
            "{argv:?}"
        );
    };
    refused(&["-i", "in.vcf", "--"], "--");
    refused(&["-i", "in.vcf", "--", "--hgvs", "x"], "--, --hgvs, x");
    refused(&["--", "-i", "in.vcf"], "--, -i, in.vcf");
    refused(&["--", "--help"], "--, --help");
    let (code, stderr) = run(&["--help", "--"]);
    assert_eq!(code, Some(0), "{stderr}");
    let (code, stderr) = refusal(&["--pick_order", ""]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(
        stderr.starts_with("ERROR: [") && !stderr.contains("pick_order"),
        "the run passes the parser and fails on its inputs: {stderr}"
    );
    let (code, stderr) = refusal(&["--pick_order", "canonical,,rank"]);
    assert_eq!(code, Some(2), "{stderr}");
    assert!(
        stderr.starts_with(
            "error: a value is required for '--pick_order <PICK_ORDER>' but none was supplied"
        ),
        "{stderr}"
    );
}

/// `--convert` is refused with Ensembl VEP's deprecation text (Config.pm:977-985:
/// `ERROR: --convert has been deprecated`), the `ERROR: ` prefix the contract's.
#[test]
fn convert_is_deprecated_error() {
    assert_refused(&["--convert", "vcf"], "--convert has been deprecated");
    assert_refused(&["--convert=ensembl"], "--convert has been deprecated");
}

/// The database connection flags (Config.pm:86-92) are refused one by one, the
/// short forms `-u` (user|u) and `--pass` (password|pass) with their long
/// flag's text.
#[test]
fn database_flags_are_refused() {
    let refused_message = |flag: &str| -> &'static str {
        match error_table_entry(flag).map(|entry| entry.treatment) {
            Some(Treatment::Refused { message, .. }) => message,
            other => panic!("{flag}: {other:?}"),
        }
    };
    assert_refused(&["--database"], refused_message("database"));
    assert_refused(
        &["--host", "ensembldb.ensembl.org"],
        refused_message("host"),
    );
    assert_refused(&["--port", "3337"], refused_message("port"));
    assert_refused(&["--user", "anonymous"], refused_message("user"));
    assert_refused(&["-u", "anonymous"], refused_message("user"));
    assert_refused(&["--password", "x"], refused_message("password"));
    assert_refused(&["--pass", "x"], refused_message("password"));
    assert_refused(&["--registry", "registry.pm"], refused_message("registry"));
    // A refused flag is refused however many flags accompany it, and before
    // the option rules: --database with --offline is not the incompatibility.
    assert_refused(
        &["--offline", "--cache", "--database"],
        refused_message("database"),
    );
}

/// Every refused entry of the error table names its flag first, says what the
/// flag needs or does after a colon, and gives the route; every entry cites its
/// Ensembl VEP source; and the binary prints exactly the table's line.
#[test]
fn every_error_table_entry_names_flag_need_and_route() {
    assert!(ERROR_TABLE.len() >= 24, "{} entries", ERROR_TABLE.len());
    for entry in ERROR_TABLE {
        let flag = format!("--{}", entry.flag);
        assert!(
            entry.ensembl_vep.contains(".pm:") || entry.ensembl_vep.contains("vep:"),
            "{flag}: the Ensembl VEP column cites a source line, got {:?}",
            entry.ensembl_vep
        );
        assert!(!entry.route.is_empty(), "{flag}: a route");
        match entry.treatment {
            Treatment::Refused { component, message } => {
                assert_eq!(component, "cli", "{flag}");
                assert!(message.starts_with(&flag), "{flag}: {message}");
                assert!(
                    message.contains(": ") || message.ends_with("deprecated"),
                    "{flag}: the message says what the flag needs: {message}"
                );
                let argv: Vec<&str> = if entry.declared && takes_value(entry.flag) {
                    vec![flag.as_str(), "x"]
                } else {
                    vec![flag.as_str()]
                };
                assert_refused(&argv, message);
            }
            Treatment::AcceptedWithoutEffect { reason } => {
                assert!(!reason.is_empty(), "{flag}");
                assert!(entry.declared, "{flag}: an accepted flag is declared");
            }
        }
    }
}

/// Whether a declared flag takes a value, asked of the parser.
fn takes_value(flag: &str) -> bool {
    use clap::CommandFactory;
    let mut command = vep_cli::args::Args::command();
    command.build();
    let arg = command
        .get_arguments()
        .find(|a| a.get_long() == Some(flag))
        .unwrap_or_else(|| panic!("--{flag} is declared"));
    !matches!(arg.get_action(), clap::ArgAction::SetTrue)
}

/// `--db_version` selects an Ensembl database version (Config.pm:92).
#[test]
fn db_version_is_refused() {
    assert_refused(
        &["--db_version", "116"],
        "--db_version is not supported: it selects an Ensembl database version; --cache_version selects a cache",
    );
}

/// `--af_exac` left Ensembl VEP's option list with release 107, the release
/// that added `--af_gnomade`; Ensembl VEP 116.3 names it only in
/// `%INCOMPATIBLE{database}` (Config.pm:618). It is refused before parsing, in
/// both spellings, naming its replacement.
#[test]
fn af_exac_is_refused_as_removed_upstream() {
    let message = "--af_exac was removed in Ensembl VEP release 107: use --af_gnomade";
    assert_refused(&["--af_exac"], message);
    assert_refused(&["--af_exac=1"], message);
    assert_refused(&["--everything", "--af_exac"], message);
}

/// `--genomes` selects the Ensembl Genomes host (Config.pm:390-396).
#[test]
fn genomes_is_refused() {
    assert_refused(
        &["--genomes"],
        "--genomes is not supported: it selects the Ensembl Genomes database host; vep-rs reads a local cache only",
    );
}

/// `--gp` reads coordinates from INFO/GP (Parser/VCF.pm:271-285).
#[test]
fn gp_is_refused() {
    assert_refused(
        &["--gp"],
        "--gp is not supported: it reads a GRCh37 chr:pos from the INFO/GP field that some 1000 Genomes VCFs carry",
    );
}

/// `--is_multispecies` describes a multi-species database (Config.pm:101).
#[test]
fn is_multispecies_is_refused() {
    assert_refused(
        &["--is_multispecies", "1"],
        "--is_multispecies is not supported: it describes a multi-species Ensembl database; vep-rs reads a local cache only",
    );
}

/// `--lrg` needs the core database; Ensembl VEP refuses it offline too
/// (Runner.pm:766 `Cannot map to LRGs in offline mode`).
#[test]
fn lrg_is_refused() {
    assert_refused(
        &["--lrg"],
        "--lrg is not supported: LRG mapping needs the Ensembl core database (Ensembl VEP: Cannot map to LRGs in offline mode)",
    );
}

/// Asserts that a transcript-source flag, given a file in both value spellings
/// and with an output path, is refused with the error table's line, which
/// opens with `opening` and ends with the cache recipe, that no output file is
/// written, and that `--help` lists the flag with its value; the message.
fn assert_transcript_source_flag_refused(flag: &str, opening: &str) -> &'static str {
    let message = match error_table_entry(flag).map(|entry| entry.treatment) {
        Some(Treatment::Refused { message, .. }) => message,
        other => panic!("--{flag}: {other:?}"),
    };
    assert!(message.starts_with(opening), "--{flag}: {message}");
    assert!(message.ends_with(GFF_GTF_RECIPE), "--{flag}: {message}");
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let out = out.to_str().unwrap();
    let separate = format!("--{flag}");
    let attached = format!("--{flag}=models.gz");
    for argv in [
        vec![separate.as_str(), "models.gz", "-o", out],
        vec![attached.as_str(), "-o", out],
    ] {
        assert_refused(&argv, message);
        assert!(!Path::new(out).exists(), "{argv:?}: no output is written");
    }
    let help = common::isolated_vep_anywhere()
        .arg("--help")
        .output()
        .unwrap();
    let help = String::from_utf8(help.stdout).unwrap();
    assert!(
        help.contains(&format!("--{flag} <FILE>")),
        "--help lists --{flag}"
    );
    message
}

/// A file of the corpus whose runs read a gene-model file beside the cache:
/// its cache slice, its GTF and GFF3, and the FASTA holding their sequence.
fn gene_model_corpus(file: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-gene-model-beside-cache")
        .join(file)
}

/// Asserts that a transcript-source flag naming `file`, in both value
/// spellings, runs beside the corpus's cache and writes rows carrying
/// `SOURCE=<the file's name>`, the source's short name (`File.pm:154-157`),
/// under the header lines Ensembl VEP writes for the column and the source,
/// and that `--help` lists the flag with its value.
fn assert_transcript_source_flag_read_beside_a_cache(flag: &str, file: &str) {
    let dir = tempfile::tempdir().unwrap();
    let path = gene_model_corpus(file);
    let path = path.to_str().unwrap();
    let separate = format!("--{flag}");
    let attached = format!("--{flag}={path}");
    for argv in [vec![separate.as_str(), path], vec![attached.as_str()]] {
        let out = dir.path().join("out.txt");
        let output = common::isolated_vep(dir.path())
            .arg("-i")
            .arg(gene_model_corpus("variants.vcf"))
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(gene_model_corpus("json_cache"))
            .arg("--fasta")
            .arg(gene_model_corpus("reference.fa.gz"))
            .args(["--assembly", "GRCh38", "--force_overwrite", "--quiet"])
            .args(&argv)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{argv:?}: {stderr}");
        let text = std::fs::read_to_string(&out).unwrap();
        for header in [
            "## SOURCE : Source of transcript".to_string(),
            format!("## {file} : [PATH]/{file}"),
        ] {
            assert!(
                text.lines().any(|line| line == header),
                "{argv:?}: {header}"
            );
        }
        let item = format!("SOURCE={file}");
        let rows: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
        let from_file = rows
            .iter()
            .filter(|row| {
                row.split('\t')
                    .next_back()
                    .unwrap()
                    .split(';')
                    .any(|i| i == item)
            })
            .count();
        assert!(
            from_file > 0 && from_file < rows.len(),
            "{argv:?}: {from_file} of {} rows carry {item}",
            rows.len()
        );
    }
    let help = common::isolated_vep_anywhere()
        .arg("--help")
        .output()
        .unwrap();
    let help = String::from_utf8(help.stdout).unwrap();
    assert!(
        help.contains(&format!("--{flag} <FILE>")),
        "--help lists --{flag}"
    );
}

/// Asserts that a transcript-source flag in a run with no cache where its
/// flags name one, here the default cache root of a home that holds no
/// cache, is refused under `custom` with the cache's own message and both
/// routes, naming the cache or building one from the file, and writes no
/// output.
fn assert_transcript_source_flag_refused_without_a_cache(flag: &str, file: &str) {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.txt");
    let out = out.to_str().unwrap();
    let separate = format!("--{flag}");
    let attached = format!("--{flag}={file}");
    for argv in [
        vec![separate.as_str(), file, "-o", out],
        vec![attached.as_str(), "-o", out],
    ] {
        let (code, stderr) = refusal(&argv);
        assert_eq!(code, Some(1), "{argv:?}: {stderr}");
        assert!(
            stderr.starts_with(
                "ERROR: [custom] a GFF or GTF transcript source is read beside a cache, and this \
                 run has none ("
            ),
            "{argv:?}: {stderr}"
        );
        assert!(
            stderr.ends_with(&format!(
                "; name the cache of the run with --json_cache <dir>, or with --cache and \
                 --dir_cache, or {GFF_GTF_RECIPE}\n"
            )),
            "{argv:?}: {stderr}"
        );
        assert_eq!(stderr.lines().count(), 1, "{argv:?}: {stderr}");
        assert!(!Path::new(out).exists(), "{argv:?}: no output is written");
    }
}

/// `--gff` is Ensembl VEP's shortcut to a `--custom file=<gff>,format=gff`
/// transcript source (Config.pm:260; the option set at :532-536), built after
/// the cache's sources and annotating beside them
/// (AnnotationSourceAdaptor.pm:94-105).
#[test]
fn gff_is_read_beside_a_cache() {
    assert_transcript_source_flag_read_beside_a_cache("gff", "novel_models.gff3.gz");
}

/// `--gtf` is Ensembl VEP's shortcut to a `--custom file=<gtf>,format=gtf`
/// transcript source (Config.pm:261; the option set at :538-541), built after
/// the cache's sources and annotating beside them
/// (AnnotationSourceAdaptor.pm:94-105).
#[test]
fn gtf_is_read_beside_a_cache() {
    assert_transcript_source_flag_read_beside_a_cache("gtf", "chr21_genes.gtf.gz");
}

/// A run naming `--gff` with no cache where its flags name one: the file is
/// read beside a cache, and a file alone is annotated against a cache built
/// from it, so the refusal names both routes.
#[test]
fn gff_without_a_cache_is_refused_with_both_routes() {
    assert_transcript_source_flag_refused_without_a_cache("gff", "genes.gff3.gz");
}

/// A run naming `--gtf` with no cache where its flags name one is refused
/// as `--gff` is, with both routes.
#[test]
fn gtf_without_a_cache_is_refused_with_both_routes() {
    assert_transcript_source_flag_refused_without_a_cache("gtf", "genes.gtf.gz");
}

/// `--bam` (Config.pm:253) hands an alignment file to the cache, GFF, GTF and
/// database transcript sources (CacheDir.pm:153; AnnotationSourceAdaptor.pm:162,
/// 402), whose RefSeq models Ensembl VEP corrects against it as they are read
/// (AnnotationType/Transcript.pm:421-631 applies the edits as SeqEdits), and
/// its option set switches on `--use_transcript_ref` (Config.pm:573-579).
/// The RefSeq cache Ensembl distributes was dumped with those edits applied, the
/// cache vep-rs reads under `--refseq`, so the refusal names that route.
#[test]
fn bam_is_refused_naming_the_refseq_cache_route() {
    assert_transcript_source_flag_refused(
        "bam",
        "--bam is not supported: Ensembl VEP applies a BAM's alignments to RefSeq transcript models at run time; the RefSeq cache Ensembl distributes already carries those edits, which vep-rs reads from a converted <species>_refseq cache under --refseq",
    );
}

/// `--bam` opens a BAM file of transcript alignments through Bio::DB::HTS
/// (AnnotationType/Transcript.pm:391-399), edits each transcript the input
/// overlaps to the sequence of its alignment (:119, :421-631) and switches on
/// `--use_transcript_ref` (Config.pm:573-579); the refusal is the same line in
/// the two-dash, one-dash and attached spellings.
#[test]
fn bam_is_refused() {
    let message = format!(
        "--bam is not supported: Ensembl VEP applies a BAM's alignments to RefSeq transcript models at run time; the RefSeq cache Ensembl distributes already carries those edits, which vep-rs reads from a converted <species>_refseq cache under --refseq; a GFF or GTF transcript set of your own is read without the BAM's edits; {GFF_GTF_RECIPE}"
    );
    assert_refused(&["--bam", "edits.bam"], &message);
    assert_refused(&["-bam", "edits.bam"], &message);
    assert_refused(&["--bam=edits.bam"], &message);
}

/// Runs the binary on `input` against the release-116 cache, so a refusal of
/// the input's content comes after the cache resolves and before any row is
/// written; `extra` flags follow the fixed ones.
fn input_refusal(input: &Path, out: &Path, extra: &[&str]) -> (Option<i32>, String) {
    let cache = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache");
    let output = common::isolated_vep_anywhere()
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(out)
        .arg("--json_cache")
        .arg(cache)
        .args(["--assembly", "GRCh38", "--quiet"])
        .args(extra)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

fn assert_input_refused(input: &Path, out: &Path, message: &str) {
    let (code, stderr) = input_refusal(input, out, &[]);
    assert_eq!(code, Some(1), "{}: {stderr}", input.display());
    assert_eq!(
        stderr,
        format!("ERROR: [input] {message}\n"),
        "{}",
        input.display()
    );
    assert!(!out.exists(), "{}: no output is written", input.display());
}

/// `--format` takes Ensembl VEP's `%VALID{format}` values (Config.pm:597) as a
/// clap enumeration: any other value is a usage error, exit 2, listing them, as
/// Ensembl VEP lists them after `is not a valid value for --format`
/// (Config.pm:949-955).
#[test]
fn format_invalid_value_lists_valid_values() {
    let (code, stderr) = refusal(&["--format", "foo"]);
    assert_eq!(code, Some(2), "{stderr}");
    assert!(
        stderr.contains("error: invalid value 'foo' for '--format <FORMAT>'"),
        "{stderr}"
    );
    assert!(
        stderr.contains("[possible values: ensembl, vcf, hgvs, id, spdi, region, guess]"),
        "{stderr}"
    );
}

/// A first data line of the CAid shape (`^CA\d+$`, ensembl-variation
/// Utils/VEP.pm:412-415) detects the `caid` format, whose records Ensembl VEP
/// looks up in its variation database (Parser/CAID.pm:53); its parser refuses
/// the format offline with `ERROR: Cannot use CAid format in offline mode`
/// (Parser/CAID.pm:103), the text vep-rs refuses it with under the `input`
/// component. `--format caid` is not a value the flag takes, as
/// `%VALID{format}` (Config.pm:597) lacks it.
#[test]
fn format_caid_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("caid.txt");
    std::fs::write(&input, "CA9985736\n").unwrap();
    assert_input_refused(
        &input,
        &dir.path().join("out.txt"),
        "Cannot use CAid format in offline mode",
    );
    let (code, stderr) = refusal(&["--format", "caid"]);
    assert_eq!(code, Some(2), "{stderr}");
    assert!(
        stderr.contains("error: invalid value 'caid' for '--format <FORMAT>'"),
        "{stderr}"
    );
}

/// A BCF stream (the magic `BCF` and version 2, raw or inside its BGZF blocks)
/// is refused on its first bytes, before any line is read: Ensembl VEP's
/// `%FORMAT_MAP` (Parser.pm:88-96) names text formats only.
#[test]
fn bcf_input_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let message =
        "BCF input is not supported: write a VCF with bcftools view -Oz and pass it with -i";
    let bcf = b"BCF\x02\x02\x10\x00\x00\x00##fileformat=VCFv4.2\n";
    let raw = dir.path().join("raw.bcf");
    std::fs::write(&raw, bcf).unwrap();
    assert_input_refused(&raw, &dir.path().join("raw.txt"), message);

    let mut writer = noodles::bgzf::io::Writer::new(Vec::new());
    writer.write_all(bcf).unwrap();
    let compressed = dir.path().join("records.vcf");
    std::fs::write(&compressed, writer.finish().unwrap()).unwrap();
    assert_input_refused(&compressed, &dir.path().join("compressed.txt"), message);
}

/// An option-rule violation is refused under the same contract with Ensembl
/// VEP's own text (Config.pm:964, :972).
#[test]
fn option_rule_violations_are_refused() {
    assert_refused(
        &["--quiet", "--verbose"],
        "Can't use --quiet and --verbose together",
    );
    assert_refused(&["--vcf", "--tab"], "Can't use --vcf and --tab together");
    assert_refused(
        &["--most_severe", "--everything"],
        "Can't use --most_severe and --biotype together",
    );
}

/// `--transcript_index_impl` is not an option of vep-rs: the overlap index has
/// one implementation (UCSC binning) and Ensembl VEP has no such flag, so clap
/// refuses the name as an unknown argument (exit 2) in both value spellings,
/// and a run that passes it stops before any input is read.
#[test]
fn transcript_index_impl_is_rejected() {
    for flags in [
        &["--transcript_index_impl", "bin"][..],
        &["--transcript_index_impl=coitrees"][..],
    ] {
        let (code, stderr) = refusal(flags);
        assert_eq!(code, Some(2), "{flags:?}: {stderr}");
        assert!(
            stderr.starts_with("error: unexpected argument '--transcript_index_impl'"),
            "{flags:?}: {stderr}"
        );
    }
}

/// The Parquet writer has one row shape, so `--parquet_shape` is not an option
/// of vep-rs: clap refuses the name as an unknown argument (exit 2) in both
/// value spellings, and a run that passes it stops before any input is read.
#[test]
fn parquet_shape_is_rejected() {
    for flags in [
        &["--parquet_shape", "flat"][..],
        &["--parquet_shape=nested"][..],
    ] {
        let (code, stderr) = refusal(flags);
        assert_eq!(code, Some(2), "{flags:?}: {stderr}");
        assert!(
            stderr.starts_with("error: unexpected argument '--parquet_shape'"),
            "{flags:?}: {stderr}"
        );
    }
}

/// The Parquet writer writes its file directly, with no intermediate to keep,
/// so `--keep_intermediate` is not an option of vep-rs: clap refuses the name
/// as an unknown argument (exit 2), as a switch and with a value.
#[test]
fn keep_intermediate_is_rejected() {
    for flags in [&["--keep_intermediate"][..], &["--keep_intermediate=1"][..]] {
        let (code, stderr) = refusal(flags);
        assert_eq!(code, Some(2), "{flags:?}: {stderr}");
        assert!(
            stderr.starts_with("error: unexpected argument '--keep_intermediate'"),
            "{flags:?}: {stderr}"
        );
    }
}

/// The Parquet writer bounds its row groups by encoded size, so
/// `--parquet_row_group_size` is not an option of vep-rs: clap refuses the name
/// as an unknown argument (exit 2) in both value spellings.
#[test]
fn parquet_row_group_size_is_rejected() {
    for flags in [
        &["--parquet_row_group_size", "2048"][..],
        &["--parquet_row_group_size=2048"][..],
    ] {
        let (code, stderr) = refusal(flags);
        assert_eq!(code, Some(2), "{flags:?}: {stderr}");
        assert!(
            stderr.starts_with("error: unexpected argument '--parquet_row_group_size'"),
            "{flags:?}: {stderr}"
        );
    }
}

/// `af_exac 1` from a `--config` file or a `--dir`'s `vep.ini` reaches the
/// refusal `--af_exac` reaches typed, with the identical line and exit 1, as
/// every flag read from a source does: `key_storage` in `argv.rs` stores a key
/// naming a flag the error table refuses before parsing under its own name, and
/// `tokens_of` lays it out as the switch `--af_exac` for the scan in
/// `config::parse_command_line` to refuse. A `0` value adds nothing, the rule
/// for a switch set to a false value from a source. Ensembl VEP 116.3 drops the
/// key from a file or the environment, since `GetOptions` run against the
/// `VEP_PARAMS` list (Config.pm:1046) matches it to no option and
/// `is_valid_param` (:1047) reads the empty result as invalid, so
/// `read_config_from_file` skips the line (:1096, :1100); typed, it dies on it
/// as an unexpected parameter (vep:30).
#[test]
fn af_exac_is_refused_from_a_configuration_file() {
    let message = "--af_exac was removed in Ensembl VEP release 107: use --af_gnomade";
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("flags.ini");
    std::fs::write(&file, "af_exac 1\n").unwrap();
    assert_refused(&["--config", file.to_str().unwrap()], message);

    let ini_dir = dir.path().join("dir");
    std::fs::create_dir_all(&ini_dir).unwrap();
    std::fs::write(ini_dir.join("vep.ini"), "af_exac 1\n").unwrap();
    assert_refused(&["--dir", ini_dir.to_str().unwrap()], message);

    let off = dir.path().join("off.ini");
    std::fs::write(&off, "af_exac 0\n").unwrap();
    let (code, stderr) = refusal(&["--config", off.to_str().unwrap()]);
    assert_eq!(code, Some(1), "the run fails later, on its cache: {stderr}");
    assert!(
        !stderr.contains("af_exac"),
        "a false value adds no token: {stderr}"
    );
}
