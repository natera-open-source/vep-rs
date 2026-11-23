// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The conversion that a plugin's refusal of a store names, run as the
//! refusal names it on a table with no tabix index.
//!
//! Beside a table with no index a plugin whose store is refused is not
//! loaded, and the line that ends its load names the conversion that rebuilds
//! the store and the flag that conversion needs on such a table. Each case
//! stages a committed table of made-up values with no index beside a store
//! the plugin refuses, reads the refusal a `vep` run writes, runs
//! `convert-plugin-data` as that line names it, and holds the conversion's
//! exit status and the store it writes: without the flag the conversion ends
//! with exit code 1 and leaves the refused store as it stood, with the flag
//! it ends with exit code 0, and the next run loads the plugin from the store
//! and writes the rows the table gives through its index. A table whose
//! index stands under the name of the file's stem is refused for its `.tbi`
//! all the same; its line names no flag, and the conversion it names ends
//! with exit code 0 as named.
//!
//! REVEL checks its table before it opens a store and TranscriptAnnotator
//! opens its store first; CADD's own cases stand in
//! `cadd_store_of_a_table_without_an_index.rs`.

use std::path::{Path, PathBuf};
use std::process::Command;

use vep_cli::warnings::{self, Target};

/// What a refusal that names the conversion ends with beside a table with no
/// index.
const NO_INDEX_CONVERSION: &str = "a data file with no index is converted with --no-verify";

/// What the converter prints when it cannot open the table for its check.
const NO_TABLE_FOR_THE_CHECK: &str =
    "failed to open the table for verification, which --no-verify skips";

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

/// The bytes of a store's two files, the data file first.
fn store_bytes(store: &Path) -> (Vec<u8>, Vec<u8>) {
    (
        std::fs::read(store).unwrap(),
        std::fs::read(store.with_extension("vpdi")).unwrap(),
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

/// One plugin's table and the layout flags its conversion takes.
struct Case {
    corpus: &'static str,
    table: &'static str,
    plugin: &'static str,
    /// The text that opens each value the plugin writes into a row.
    field: &'static str,
    flags: &'static [&'static str],
}

/// The refusal that ends a plugin's load beside `table`, which has no index,
/// and `store`, which is refused with `store_refusal`.
fn refusal_beside(case: &Case, table: &Path, store: &Path, store_refusal: &str) -> String {
    format!(
        "WARNING: Failed to instantiate plugin {}: tabix index {}.tbi not found; a plugin data \
         file needs its .tbi index beside it, so index it with tabix without -C; the store {} \
         beside the table is refused: {store_refusal}; {NO_INDEX_CONVERSION}",
        case.plugin,
        table.display(),
        store.display()
    )
}

/// `case`'s table with no index, beside the store `refused` leaves: the run
/// writes `store_refusal` in the line that ends the plugin's load; the
/// conversion with `flags` ends with exit code 1 for want of the index and
/// leaves the refused store as it stood; under the flag the line names it
/// ends with exit code 0, and the next run writes the rows of `through_index`.
fn the_named_conversion_writes_a_store_the_plugin_reads(
    case: &Case,
    dir: &Path,
    through_index: &Run,
    refused: impl FnOnce(&Path),
    store_refusal: &str,
    flags: &[&str],
) {
    std::fs::create_dir_all(dir).unwrap();
    let corpus = corpus(case.corpus);
    let table = staged_table(dir, &corpus.join(case.table), false);
    let store = store_beside(&table);
    let value = format!("{},file={}", case.plugin, table.display());
    refused(&table);

    let refused_run = annotate(dir, "refused", &corpus, &value);
    assert_eq!(
        refused_run.warnings.lines().next().unwrap_or_default(),
        refusal_beside(case, &table, &store, store_refusal),
        "{}",
        case.plugin
    );
    assert!(
        refused_run.rows.iter().all(|row| !row.contains(case.field)),
        "{}: the plugin is not loaded",
        case.plugin
    );

    let before = store_bytes(&store);
    let (code, stderr) = convert(&table, flags);
    assert_eq!(code, Some(1), "{}: {stderr}", case.plugin);
    assert!(
        stderr.contains(NO_TABLE_FOR_THE_CHECK),
        "{}: {stderr}",
        case.plugin
    );
    assert!(
        store_bytes(&store) == before,
        "{}: a conversion that ends before it reads the table leaves the store as it stood",
        case.plugin
    );

    let (code, stderr) = convert(&table, &[flags, &["--no-verify"]].concat());
    assert_eq!(code, Some(0), "{}: {stderr}", case.plugin);

    let read = annotate(dir, "read", &corpus, &value);
    assert_eq!(
        read.warnings, through_index.warnings,
        "{}: the plugin loads from the store",
        case.plugin
    );
    assert_eq!(read.rows, through_index.rows, "{}", case.plugin);
}

/// Converts `table` with `flags` and no check, then sets the version its
/// store's index carries to 4, one below the version a plugin reads.
fn store_of_an_earlier_index_version(table: &Path, flags: &[&str]) {
    let (code, stderr) = convert(table, &[flags, &["--no-verify"]].concat());
    assert_eq!(code, Some(0), "{stderr}");
    let index_path = store_beside(table).with_extension("vpdi");
    let mut index = std::fs::read(&index_path).unwrap();
    index[4..8].copy_from_slice(&4u32.to_le_bytes());
    std::fs::write(&index_path, index).unwrap();
}

/// A store refused for its index version, for REVEL and for
/// TranscriptAnnotator, and a store REVEL refuses for its columns, built on
/// the converter's default flags where the plugin reads the third column as
/// the position. The line of the first names the conversion, which is run
/// with the layout flags of the table; the line of the second names the
/// flags itself, which the conversion is run with. Last, REVEL's table with
/// its index under the name of the file's stem, beside a store of an earlier
/// index version: the line names the conversion and no flag, and the
/// conversion runs as named. The runs go one after another, since each
/// points the one warnings sink the process holds at its own file.
#[test]
fn the_conversion_a_refusal_names_on_a_table_with_no_index_writes_a_store_the_plugin_reads() {
    const RECONVERT: &str = "convert the table again with vep-cache-converter convert-plugin-data";
    let cases = [
        Case {
            corpus: "GRCh38-plugin-revel",
            table: "plugin_data/revel.tsv.gz",
            plugin: "REVEL",
            field: "REVEL=",
            flags: &["--start-col", "2"],
        },
        Case {
            corpus: "GRCh38-plugin-transcript-annotator",
            table: "data/transcript_annotator.tsv.gz",
            plugin: "TranscriptAnnotator",
            field: "transcript_annotator_",
            flags: &[],
        },
    ];
    let tmp = tempfile::tempdir().unwrap();
    for case in &cases {
        let root = tmp.path().join(case.plugin);
        let with_index = root.join("with-index");
        std::fs::create_dir_all(&with_index).unwrap();
        let corpus = corpus(case.corpus);
        let indexed = staged_table(&with_index, &corpus.join(case.table), true);
        let through_index = annotate(
            &with_index,
            "table",
            &corpus,
            &format!("{},file={}", case.plugin, indexed.display()),
        );
        assert!(
            through_index
                .rows
                .iter()
                .any(|row| row.contains(case.field)),
            "{}: the table answers",
            case.plugin
        );

        the_named_conversion_writes_a_store_the_plugin_reads(
            case,
            &root.join("index-version"),
            &through_index,
            |table| store_of_an_earlier_index_version(table, case.flags),
            &format!("unsupported index version: 4 (expected 5); {RECONVERT}"),
            case.flags,
        );
    }

    // REVEL's table converted on the default flags, which place a row by its
    // second column: the refusal names the flags that rebuild the store.
    let revel = &cases[0];
    let root = tmp.path().join("REVEL-columns");
    let with_index = root.join("with-index");
    std::fs::create_dir_all(&with_index).unwrap();
    let corpus = corpus(revel.corpus);
    let indexed = staged_table(&with_index, &corpus.join(revel.table), true);
    let through_index = annotate(
        &with_index,
        "table",
        &corpus,
        &format!("REVEL,file={}", indexed.display()),
    );
    let named = "--chr-col 0 --start-col 2";
    the_named_conversion_writes_a_store_the_plugin_reads(
        revel,
        &root.join("columns"),
        &through_index,
        |table| {
            let (code, stderr) = convert(table, &["--no-verify"]);
            assert_eq!(code, Some(0), "{stderr}");
        },
        &format!(
            "store columns `--chr-col 0 --start-col 1` differ from the table's `{named}`; \
             {RECONVERT} and the table's flags"
        ),
        &named.split(' ').collect::<Vec<_>>(),
    );

    // REVEL's table with its index under the name of the file's stem, beside
    // a store of an earlier index version: the plugin refuses the table for
    // its `.tbi`, the conversion's check reads the index where it stands, so
    // the line ends on the conversion and the conversion runs as named.
    let root = tmp.path().join("REVEL-index-on-the-stem");
    std::fs::create_dir_all(&root).unwrap();
    let table = staged_table(&root, &corpus.join(revel.table), true);
    let tbi = PathBuf::from(format!("{}.tbi", table.display()));
    std::fs::rename(&tbi, table.with_extension("tbi")).unwrap();
    let store = store_beside(&table);
    let value = format!("REVEL,file={}", table.display());
    store_of_an_earlier_index_version(&table, revel.flags);

    let refused_run = annotate(&root, "refused", &corpus, &value);
    assert_eq!(
        refused_run.warnings.lines().next().unwrap_or_default(),
        format!(
            "WARNING: Failed to instantiate plugin REVEL: tabix index {} not found; a plugin \
             data file needs its .tbi index beside it, so index it with tabix without -C; the \
             store {} beside the table is refused: unsupported index version: 4 (expected 5); \
             {RECONVERT}",
            tbi.display(),
            store.display()
        ),
        "an index on the stem: the line names no flag"
    );
    assert!(
        refused_run
            .rows
            .iter()
            .all(|row| !row.contains(revel.field)),
        "an index on the stem: the plugin is not loaded"
    );

    let (code, stderr) = convert(&table, revel.flags);
    assert_eq!(code, Some(0), "an index on the stem: {stderr}");
    let read = annotate(&root, "read", &corpus, &value);
    assert_eq!(
        read.warnings, through_index.warnings,
        "an index on the stem: the plugin loads from the store"
    );
    assert_eq!(read.rows, through_index.rows, "an index on the stem");
}
