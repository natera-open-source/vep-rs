// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Flags vep-rs declares and refuses, through the `vep` binary: each exits 1
//! with exactly one line `ERROR: [cli] <message>` on standard error, before any
//! input is read. The messages are the error table's (`vep_cli::config`).

use std::process::{Command, Stdio};

use vep_cli::config::{error_table_entry, Treatment, ERROR_TABLE};

/// Runs the binary with `flags` ahead of an input that does not exist, so a flag
/// that were accepted would fail later, on the input, under another component.
fn refusal(flags: &[&str]) -> (Option<i32>, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_vep"))
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

/// Ensembl VEP resolves an ambiguous HGVSp through its database
/// (Parser/HGVS.pm:178-186, 224); vep-rs reads genomic HGVS only.
#[test]
fn ambiguous_hgvs_is_refused() {
    assert_refused(
        &["--ambiguous_hgvs"],
        "--ambiguous_hgvs is not supported: it resolves an HGVSp input to every matching variant through the Ensembl database; HGVS input is genomic (g.) only",
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
