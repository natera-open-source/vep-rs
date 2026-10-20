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

use vep_cli::config::{error_table_entry, Treatment, ERROR_TABLE};

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
    assert!(ERROR_TABLE.len() >= 27, "{} entries", ERROR_TABLE.len());
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

/// `--af_exac` left `@VEP_PARAMS` with release 115 and survives only in
/// `%INCOMPATIBLE{database}` (Config.pm:618); it is refused before parsing, in
/// both spellings, naming its replacement.
#[test]
fn af_exac_is_refused_as_removed_upstream() {
    let message = "--af_exac was removed in Ensembl VEP release 115: use --af_gnomade";
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
        "--gp is not supported: it reads 1000 Genomes pilot coordinates from INFO/GP, which no current VCF carries",
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
/// Utils/VEP.pm:412-415) detects the `caid` format, which resolves through the
/// ClinGen allele registry; Ensembl VEP refuses it offline (Parser/CAID.pm:103)
/// and vep-rs refuses it by name. `--format caid` is not a value the flag
/// takes, as `%VALID{format}` (Config.pm:597) lacks it.
#[test]
fn format_caid_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("caid.txt");
    std::fs::write(&input, "CA9985736\n").unwrap();
    assert_input_refused(
        &input,
        &dir.path().join("out.txt"),
        "CAid input (CA...) resolves through the ClinGen allele registry, a remote service; vep-rs does not resolve it",
    );
    let (code, stderr) = refusal(&["--format", "caid"]);
    assert_eq!(code, Some(2), "{stderr}");
    assert!(
        stderr.contains("error: invalid value 'caid' for '--format <FORMAT>'"),
        "{stderr}"
    );
}

/// An HGVS input with an r. notation is refused by name when the line is read,
/// however many lines follow it, and no record is annotated: Ensembl VEP's
/// resolver throws on the type (ensembl-variation
/// DBSQL/VariationFeatureAdaptor.pm:1989, `HGVS type r not recognised`), which
/// its parser warns about and skips (Parser/HGVS.pm:204-216), and HGVS input
/// itself is refused offline (Parser/HGVS.pm:104); here the run stops. The run
/// carries the reference FASTA the format needs, so the refusal is the type's,
/// not the format's.
#[test]
fn hgvs_rna_input_is_refused() {
    let expected_stderr =
        "ERROR: [input] HGVS r. notation is not supported (Ensembl VEP refuses it too: HGVS type r not recognised)\n";
    let fasta = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-input-hgvs/reference.fa.gz");
    let with_fasta = ["--fasta", fasta.to_str().unwrap()];
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rna.hgvs");
    std::fs::write(&input, "ENST00000366667:r.803C>T\n21:g.30282169C>A\n").unwrap();
    let out = dir.path().join("out.txt");
    let (code, stderr) = input_refusal(&input, &out, &with_fasta);
    assert_eq!(code, Some(1), "{stderr}");
    assert_eq!(stderr, expected_stderr);
    assert_no_rows(&out);

    // The refusal ends a run whose earlier records were read.
    let later = dir.path().join("rna_later.hgvs");
    std::fs::write(&later, "21:g.30282169C>A\nENST00000366667:r.803C>T\n").unwrap();
    let later_out = dir.path().join("later.txt");
    let (code, stderr) = input_refusal(&later, &later_out, &with_fasta);
    assert_eq!(code, Some(1), "{stderr}");
    assert_eq!(stderr, expected_stderr);
    assert_no_rows(&later_out);
}

/// Asserts that `out` holds no data row: the run wrote no file, or its header
/// alone.
fn assert_no_rows(out: &Path) {
    let written = std::fs::read_to_string(out).unwrap_or_default();
    let rows = written.lines().filter(|l| !l.starts_with('#')).count();
    assert_eq!(rows, 0, "no record is annotated: {written:?}");
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
