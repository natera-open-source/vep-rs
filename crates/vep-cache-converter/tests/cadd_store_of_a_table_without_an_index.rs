// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The conversion that CADD's refusal names, run as the refusal names it.
//!
//! Beside a table with no tabix index CADD refuses a store built on other
//! columns than the table's header names, and the line that ends its load
//! names the flags that rebuild the store and the flag the conversion needs
//! on such a table. Each test stages a committed table of made-up scores
//! with no index, converts it on other flags, reads the refusal a `vep` run
//! writes, runs `convert-plugin-data` with what that line names, and holds
//! the conversion's exit status and the store it writes: the next run loads
//! the plugin from that store and writes the rows the table gives through
//! its index.

use std::path::{Path, PathBuf};
use std::process::Command;

use vep_cli::warnings::{self, Target};

fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{name}"))
}

/// Copies a corpus's table into `dir`, with its `.tbi` where `indexed`;
/// returns the copy's path.
fn staged_table(dir: &Path, table: &Path, indexed: bool) -> PathBuf {
    let name = table.file_name().unwrap().to_string_lossy().into_owned();
    let staged = dir.join(&name);
    std::fs::copy(table, &staged).unwrap();
    if indexed {
        std::fs::copy(
            table.with_file_name(format!("{name}.tbi")),
            dir.join(format!("{name}.tbi")),
        )
        .unwrap();
    }
    staged
}

/// The store a plugin looks for beside `table`.
fn store_beside(table: &Path) -> PathBuf {
    PathBuf::from(table.to_str().unwrap().replace(".tsv.gz", ".tsv.vpd"))
}

/// `vep-cache-converter convert-plugin-data` on `table` with `flags`, writing
/// the store beside it; returns the exit code and what the run printed on
/// standard error.
fn convert(table: &Path, flags: &[&str]) -> (Option<i32>, String) {
    let run = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .arg("convert-plugin-data")
        .arg("--input")
        .arg(table)
        .arg("--output")
        .arg(store_beside(table))
        .args(flags)
        .output()
        .unwrap();
    (
        run.status.code(),
        String::from_utf8_lossy(&run.stderr).into_owned(),
    )
}

/// What one run left: its warnings file whole (empty when it wrote none) and
/// its data rows.
struct Run {
    warnings: String,
    rows: Vec<String>,
}

/// Runs `vep` through the library on a corpus's records and cache with one
/// `--plugin` value. The run's warnings go to the file beside its output
/// alone, where the `vep` binary puts them under `--quiet`.
fn annotate(dir: &Path, label: &str, corpus: &Path, plugin: &str) -> Run {
    let out = dir.join(format!("{label}.txt"));
    let warnings_file = dir.join(format!("{label}.txt_warnings.txt"));
    let input = corpus.join("variants.vcf");
    let json_cache = corpus.join("json_cache");
    let argv = [
        "vep",
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--offline",
        "--json_cache",
        json_cache.to_str().unwrap(),
        "--species",
        "homo_sapiens",
        "--assembly",
        "GRCh38",
        "--force_overwrite",
        "--no_stats",
        "--quiet",
        "--dir_plugins",
        dir.to_str().unwrap(),
        "--plugin",
        plugin,
    ];
    warnings::global().configure(Target::File(warnings_file.clone()), true, false);
    let config = vep_cli::Config::from_argv(argv).unwrap();
    vep_cli::runner::Runner::new(config).unwrap().run().unwrap();
    let warnings = std::fs::read_to_string(&warnings_file).unwrap_or_default();
    let rows = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_owned)
        .collect();
    Run { warnings, rows }
}

/// One table: a store of it converted on `other_flags`, which the table's
/// header does not name, is refused with no index beside the table; the
/// refusal names `named_flags`; the conversion with them ends with exit code
/// 1 for want of the index and with exit code 0 under the flag the refusal
/// names; and the store it then writes gives the table's rows.
fn the_conversion_the_refusal_names_writes_a_store_the_plugin_reads(
    corpus_name: &str,
    table: &str,
    key: &str,
    other_flags: &[&str],
    named_flags: &str,
) {
    let corpus = corpus(corpus_name);
    let tmp = tempfile::tempdir().unwrap();
    let with_index = tmp.path().join("with-index");
    let without_index = tmp.path().join("without-index");
    std::fs::create_dir_all(&with_index).unwrap();
    std::fs::create_dir_all(&without_index).unwrap();

    let indexed = staged_table(&with_index, &corpus.join(table), true);
    let through_index = annotate(
        &with_index,
        "table",
        &corpus,
        &format!("CADD,{key}={}", indexed.display()),
    );
    assert!(
        through_index.rows.iter().any(|row| row.contains("CADD_")),
        "{corpus_name}: the table scores its records"
    );

    let table = staged_table(&without_index, &corpus.join(table), false);
    let store = store_beside(&table);
    let value = format!("CADD,{key}={}", table.display());
    let (code, stderr) = convert(&table, &[other_flags, &["--no-verify"]].concat());
    assert_eq!(code, Some(0), "{corpus_name}: {stderr}");

    let refused = annotate(&without_index, "refused", &corpus, &value);
    let refusal = refused.warnings.lines().next().unwrap_or_default();
    let built_on = if other_flags.is_empty() {
        "--chr-col 0 --start-col 1".to_string()
    } else {
        format!("--chr-col 0 --start-col 1 {}", other_flags.join(" "))
    };
    assert_eq!(
        refusal,
        format!(
            "WARNING: Failed to instantiate plugin CADD: tabix index {}.tbi not found; a plugin \
             data file needs its .tbi index beside it, so index it with tabix without -C; the \
             store {} beside the table is refused: store columns `{built_on}` differ from the \
             table's `{named_flags}`; convert the table again with vep-cache-converter \
             convert-plugin-data and the table's flags; a data file with no index is converted \
             with --no-verify",
            table.display(),
            store.display()
        ),
        "{corpus_name}"
    );
    assert!(
        refused.rows.iter().all(|row| !row.contains("CADD_")),
        "{corpus_name}: the plugin is not loaded"
    );

    // The command the line names: the table's flags, read from the line.
    let named = refusal
        .split("differ from the table's `")
        .nth(1)
        .and_then(|rest| rest.split('`').next())
        .expect("the refusal names the table's flags");
    let flags: Vec<&str> = named.split(' ').collect();
    let (code, stderr) = convert(&table, &flags);
    assert_eq!(code, Some(1), "{corpus_name}: {stderr}");
    assert!(
        stderr.contains("failed to open the table for verification, which --no-verify skips"),
        "{corpus_name}: {stderr}"
    );
    let (code, stderr) = convert(&table, &[flags.as_slice(), &["--no-verify"]].concat());
    assert_eq!(code, Some(0), "{corpus_name}: {stderr}");
    assert!(store.exists() && store.with_extension("vpdi").exists());

    let read = annotate(&without_index, "read", &corpus, &value);
    assert_eq!(
        read.warnings, through_index.warnings,
        "{corpus_name}: the plugin loads from the store"
    );
    assert_eq!(read.rows, through_index.rows, "{corpus_name}");
}

/// The committed table of single positions, whose header names a `Pos`
/// column, and the committed table of spans, whose header names a `Start`
/// and an `End` column. The runs go one after another, since each points the
/// one warnings sink the process holds at its own file.
#[test]
fn the_conversion_cadds_refusal_names_on_a_table_with_no_index_writes_a_store_it_reads() {
    the_conversion_the_refusal_names_writes_a_store_the_plugin_reads(
        "GRCh38-plugin-cadd-fixture",
        "plugin_data/whole_genome_SNVs_synthetic.tsv.gz",
        "snv",
        &["--zero-based"],
        "--chr-col 0 --start-col 1",
    );
    the_conversion_the_refusal_names_writes_a_store_the_plugin_reads(
        "GRCh38-plugin-cadd-sv-fixture",
        "plugin_data/1000G_phase3_SVs_synthetic.tsv.gz",
        "sv",
        &[],
        "--chr-col 0 --start-col 1 --end-col 2 --zero-based",
    );
}
