// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A plugin's table beside a `.vpd` store that is refused for its index
//! version, through the `vep` binary. A store whose index carries another
//! version than the one the reader takes is refused whatever it holds, an
//! earlier version as a later one, and the table is then all that is left to
//! read. With `<table>.tbi` beside it the run warns of the store and reads the
//! table. Without it the plugin is not loaded, under one line that names the
//! missing index, the store, the version its index carries, the conversion
//! that rebuilds it and the flag that conversion needs on a table with no
//! index. REVEL checks its table before it opens a store
//! (`check_table_or_store`) and TranscriptAnnotator opens its store first
//! (`open_backend`); both write the same lines, and so does CADD, which reads
//! its coordinate columns before it opens a store.
//!
//! Beside a store that is accepted, CADD's table needs no index either: with
//! no index the plugin takes its coordinate columns from the names in the
//! table's header. It accepts such a store only where it was built on those
//! columns, and refuses any other by the flags that rebuild it.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use vep_builtins::binary_store::{BinaryStoreWriter, ContentIdentity};
use vep_builtins::tabix::{TabixConfig, TabixRecord};

mod common;

/// What the refusal of a store's index version ends with.
const RECONVERT: &str = "convert the table again with vep-cache-converter convert-plugin-data";

/// What a plugin's refusal of a store beside a table with no index ends with,
/// where it names the conversion: the conversion checks a store through the
/// table's index and ends without one unless it is given this flag.
const NO_INDEX_CONVERSION: &str = "a data file with no index is converted with --no-verify";

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{name}"))
}

/// Copies a corpus's table and its `.tbi` into `dir`; returns the copy's path.
fn staged_table(dir: &Path, table: &Path) -> PathBuf {
    let name = table.file_name().unwrap().to_string_lossy().into_owned();
    let staged = dir.join(&name);
    fs::copy(table, &staged).unwrap();
    fs::copy(
        table.with_file_name(format!("{name}.tbi")),
        dir.join(format!("{name}.tbi")),
    )
    .unwrap();
    staged
}

/// Writes beside `table`, under the name a plugin looks for, a store of one
/// made-up row bound to the table's content, then sets the version its index
/// carries to `version`; returns the `.vpd` path. A plugin refuses an index
/// for its version before it reads the store's columns or content, so the row
/// is never reached.
fn store_of_index_version(table: &Path, version: u32) -> PathBuf {
    let store = PathBuf::from(table.to_str().unwrap().replace(".tsv.gz", ".tsv.vpd"));
    let header: Vec<String> = ["chrom", "pos", "name"].map(String::from).to_vec();
    let record = TabixRecord::new(
        "21".to_string(),
        100,
        100,
        ["21", "100", "B"].map(String::from).to_vec(),
        Some(Arc::new(header.clone())),
    );
    let mut config = TabixConfig::default();
    config.file_path = table.to_path_buf();
    config.chr_col = 0;
    config.start_col = 1;
    config.end_col = None;
    config.ref_col = None;
    config.alt_col = None;
    config.zero_based = false;
    BinaryStoreWriter::new(store.clone())
        .write(
            &header,
            &[("21".to_string(), vec![record])],
            &config,
            &ContentIdentity::of_file(table).unwrap(),
        )
        .unwrap();
    let index_path = store.with_extension("vpdi");
    let mut index = fs::read(&index_path).unwrap();
    index[4..8].copy_from_slice(&version.to_le_bytes());
    fs::write(&index_path, index).unwrap();
    store
}

/// Writes beside `table`, under the name a plugin looks for, a store of one
/// made-up row with an index of the version a plugin reads, placed by the
/// first column and the column `start_col` and bound to the content of
/// `built_from`; returns the `.vpd` path. A plugin refuses such a store when
/// it opens it, for its columns where the plugin reads other columns and for
/// its content where `built_from` is not the table, before it reads the row.
fn current_store_built_on(table: &Path, start_col: usize, built_from: &Path) -> PathBuf {
    let store = PathBuf::from(table.to_str().unwrap().replace(".tsv.gz", ".tsv.vpd"));
    let header: Vec<String> = ["chrom", "first", "second"].map(String::from).to_vec();
    let record = TabixRecord::new(
        "21".to_string(),
        100,
        100,
        ["21", "100", "100"].map(String::from).to_vec(),
        Some(Arc::new(header.clone())),
    );
    let mut config = TabixConfig::default();
    config.file_path = table.to_path_buf();
    config.start_col = start_col;
    BinaryStoreWriter::new(store.clone())
        .write(
            &header,
            &[("21".to_string(), vec![record])],
            &config,
            &ContentIdentity::of_file(built_from).unwrap(),
        )
        .unwrap();
    store
}

/// Writes beside `table`, under the name a plugin looks for, a store of the
/// table's own rows, bound to the table's content and placed by its first two
/// columns, as the converter's default flags place them; returns the `.vpd`
/// path. The column names are the last `#` line before the first row.
fn store_of_the_tables_rows(table: &Path) -> PathBuf {
    store_of_the_tables_rows_built_on(table, None, false)
}

/// [`store_of_the_tables_rows`] with the rows placed as the converter places
/// them under `--end-col <end_col>` and, where `zero_based`, `--zero-based`:
/// a zero-based start is the base before the row's first, and without an end
/// column a row is the one base at its start.
fn store_of_the_tables_rows_built_on(
    table: &Path,
    end_col: Option<usize>,
    zero_based: bool,
) -> PathBuf {
    use std::io::Read;

    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(table).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    let names: Vec<String> = text
        .lines()
        .take_while(|line| line.starts_with('#'))
        .last()
        .expect("the table has a header line")
        .trim_start_matches('#')
        .split('\t')
        .map(String::from)
        .collect();
    let header = Arc::new(names.clone());
    let mut chromosomes: Vec<(String, Vec<TabixRecord>)> = Vec::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let columns: Vec<String> = line.split('\t').map(String::from).collect();
        let start = columns[1].parse::<u64>().unwrap() + u64::from(zero_based);
        let end = end_col.map_or(start, |col| columns[col].parse().unwrap());
        let record = TabixRecord::new(
            columns[0].clone(),
            start,
            end,
            columns,
            Some(Arc::clone(&header)),
        );
        match chromosomes.last_mut() {
            Some((chr, records)) if *chr == record.chr => records.push(record),
            _ => chromosomes.push((record.chr.clone(), vec![record])),
        }
    }
    let store = PathBuf::from(table.to_str().unwrap().replace(".tsv.gz", ".tsv.vpd"));
    let mut config = TabixConfig::default();
    config.file_path = table.to_path_buf();
    config.end_col = end_col;
    config.zero_based = zero_based;
    BinaryStoreWriter::new(store.clone())
        .write(
            &names,
            &chromosomes,
            &config,
            &ContentIdentity::of_file(table).unwrap(),
        )
        .unwrap();
    store
}

/// What one run left: its exit code, its warnings file whole (empty when it
/// wrote none), its data rows and the number of them that carry `<field>=`.
struct Run {
    code: Option<i32>,
    warnings: String,
    rows: Vec<String>,
    cells: usize,
}

/// Runs the binary on a corpus's records and cache with one `--plugin` value
/// under `--quiet`, so the warnings go to the file beside the output alone.
fn run(dir: &Path, label: &str, corpus: &Path, plugin: &str, field: &str) -> Run {
    run_on(
        dir,
        label,
        &corpus.join("variants.vcf"),
        corpus,
        plugin,
        field,
    )
}

/// [`run`] on the records of `input` in place of the corpus's own.
fn run_on(dir: &Path, label: &str, input: &Path, corpus: &Path, plugin: &str, field: &str) -> Run {
    let out = dir.join(format!("{label}.txt"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .args(["--plugin", plugin])
        .output()
        .expect("vep binary runs");
    let warnings =
        fs::read_to_string(dir.join(format!("{label}.txt_warnings.txt"))).unwrap_or_default();
    let rows: Vec<String> = fs::read_to_string(&out)
        .unwrap_or_default()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_owned)
        .collect();
    let cells = rows
        .iter()
        .filter(|row| {
            row.split('\t').nth(13).is_some_and(|extra| {
                extra
                    .split(';')
                    .any(|item| item.starts_with(&format!("{field}=")))
            })
        })
        .count();
    Run {
        code: output.status.code(),
        warnings,
        rows,
        cells,
    }
}

/// One plugin on its corpus's table, beside a store of an earlier and of a
/// later index version: with the table's `.tbi` the store's refusal is one
/// warning and the table is read; without it the plugin's load ends on one
/// line that names the index, the store, what rebuilds the store and the
/// flag that conversion needs on a table with no index.
fn table_is_read_with_its_index_and_named_with_its_store_without(
    corpus_name: &str,
    table: &str,
    plugin: &str,
    field: &str,
) {
    let corpus = corpus(corpus_name);
    for version in [4u32, 6] {
        let tmp = tempfile::tempdir().unwrap();
        let dir = tmp.path();
        let table = staged_table(dir, &corpus.join(table));
        let store = store_of_index_version(&table, version);
        let value = format!("{plugin},file={}", table.display());
        let store_refusal =
            format!("unsupported index version: {version} (expected 5); {RECONVERT}");

        let read = run(dir, "with-index", &corpus, &value, field);
        assert_eq!(read.code, Some(0), "{plugin}, version {version}");
        assert_eq!(
            read.warnings,
            format!(
                "WARNING: binary store refused, the table is read instead store={} table={} \
                 reason={store_refusal}\n",
                store.display(),
                table.display()
            ),
            "{plugin}, version {version}: the store's refusal alone"
        );
        assert!(
            read.cells > 0,
            "{plugin}, version {version}: the table answers"
        );

        let tbi = PathBuf::from(format!("{}.tbi", table.display()));
        fs::remove_file(&tbi).unwrap();
        let refused = run(dir, "without-index", &corpus, &value, field);
        assert_eq!(refused.code, Some(0), "{plugin}, version {version}");
        assert_eq!(
            refused.warnings,
            format!(
                "WARNING: Failed to instantiate plugin {plugin}: tabix index {} not found; a \
                 plugin data file needs its .tbi index beside it, so index it with tabix without \
                 -C; the store {} beside the table is refused: {store_refusal}; \
                 {NO_INDEX_CONVERSION}\n",
                tbi.display(),
                store.display()
            ),
            "{plugin}, version {version}: one line names the index, the store and the flag"
        );
        assert_eq!(
            refused.cells, 0,
            "{plugin}, version {version}: the plugin is not loaded"
        );
    }
}

#[test]
fn a_plugin_that_checks_its_table_first_names_the_refused_store() {
    table_is_read_with_its_index_and_named_with_its_store_without(
        "GRCh38-plugin-revel",
        "plugin_data/revel.tsv.gz",
        "REVEL",
        "REVEL",
    );
}

#[test]
fn a_plugin_that_opens_its_store_first_names_the_refused_store() {
    table_is_read_with_its_index_and_named_with_its_store_without(
        "GRCh38-plugin-transcript-annotator",
        "data/transcript_annotator.tsv.gz",
        "TranscriptAnnotator",
        "transcript_annotator_SIFT_score",
    );
}

#[test]
fn a_plugin_that_reads_its_columns_before_its_store_names_the_refused_store() {
    table_is_read_with_its_index_and_named_with_its_store_without(
        "GRCh38-plugin-cadd-fixture",
        "plugin_data/whole_genome_SNVs_synthetic.tsv.gz",
        "CADD",
        "CADD_PHRED",
    );
}

/// REVEL beside a store whose index is of the version it reads, which the
/// plugin's check of its table passes and which is refused when it is
/// opened: for its columns, where the store places a row by the second
/// column and the plugin reads the third as the position, and for its
/// content, where the store is bound to another file's content. With the
/// table's `.tbi` either refusal is one warning, which names the conversion,
/// and the table is read. Without it the plugin is not loaded, and the line
/// that ends its load names the conversion and the flag it needs on a table
/// with no index.
#[test]
fn a_store_refused_when_it_is_opened_names_the_flag_beside_a_table_without_an_index() {
    let corpus = corpus("GRCh38-plugin-revel");
    for refused_for in ["columns", "content"] {
        let tmp = tempfile::tempdir().unwrap();
        let dir = tmp.path();
        let table = staged_table(dir, &corpus.join("plugin_data/revel.tsv.gz"));
        let tbi = PathBuf::from(format!("{}.tbi", table.display()));
        let (store, store_refusal) = if refused_for == "columns" {
            (
                current_store_built_on(&table, 1, &table),
                format!(
                    "store columns `--chr-col 0 --start-col 1` differ from the table's \
                     `--chr-col 0 --start-col 2`; {RECONVERT} and the table's flags"
                ),
            )
        } else {
            let (holds, records) = (
                fs::metadata(&table).unwrap().len(),
                fs::metadata(&tbi).unwrap().len(),
            );
            assert_ne!(holds, records, "the two files differ in length");
            (
                current_store_built_on(&table, 2, &tbi),
                format!(
                    "store source differs: the table {} holds {holds} bytes and its store's \
                     index records {records}; {RECONVERT}",
                    table.display()
                ),
            )
        };
        let value = format!("REVEL,file={}", table.display());

        let read = run(dir, "with-index", &corpus, &value, "REVEL");
        assert_eq!(read.code, Some(0), "{refused_for}");
        assert_eq!(
            read.warnings,
            format!(
                "WARNING: binary store refused, the table is read instead store={} table={} \
                 reason={store_refusal}\n",
                store.display(),
                table.display()
            ),
            "{refused_for}: the store's refusal alone, with no flag"
        );
        assert!(read.cells > 0, "{refused_for}: the table answers");

        fs::remove_file(&tbi).unwrap();
        let refused = run(dir, "without-index", &corpus, &value, "REVEL");
        assert_eq!(refused.code, Some(0), "{refused_for}");
        assert_eq!(
            refused.warnings,
            format!(
                "WARNING: Failed to instantiate plugin REVEL: tabix index {} not found; a \
                 plugin data file needs its .tbi index beside it, so index it with tabix without \
                 -C; the store {} beside the table is refused: {store_refusal}; \
                 {NO_INDEX_CONVERSION}\n",
                tbi.display(),
                store.display()
            ),
            "{refused_for}: one line names the index, the store and the flag"
        );
        assert_eq!(refused.cells, 0, "{refused_for}: the plugin is not loaded");
    }
}

/// CADD takes its coordinate columns from its table's tabix index and, where
/// the table has none, from the names in its header, and reads the store
/// beside it where the store was built on them. With no index beside a store
/// of the table's own rows the plugin loads, the warnings file stays empty and
/// every row is the row the table gives through tabix. The table then replaced under its name, the
/// store is refused for its content and the plugin is not loaded, under one
/// line that names the missing index, the store and what differs.
#[test]
fn cadd_reads_an_accepted_store_beside_a_table_without_an_index() {
    let corpus = corpus("GRCh38-plugin-cadd-fixture");
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let table = staged_table(
        dir,
        &corpus.join("plugin_data/whole_genome_SNVs_synthetic.tsv.gz"),
    );
    let value = format!("CADD,snv={}", table.display());

    let through_tabix = run(dir, "table", &corpus, &value, "CADD_PHRED");
    assert_eq!(through_tabix.code, Some(0));
    assert_eq!(through_tabix.warnings, "");
    assert!(through_tabix.cells > 0, "the table scores its records");

    let store = store_of_the_tables_rows(&table);
    let tbi = PathBuf::from(format!("{}.tbi", table.display()));
    fs::remove_file(&tbi).unwrap();
    let through_store = run(dir, "store-without-index", &corpus, &value, "CADD_PHRED");
    assert_eq!(through_store.code, Some(0));
    assert_eq!(
        through_store.warnings, "",
        "the plugin loads from the store"
    );
    assert_eq!(through_store.rows, through_tabix.rows);

    let built_from = fs::metadata(&table).unwrap().len();
    fs::copy(
        corpus.join("plugin_data/gnomad.genomes.r4.0.indel.synthetic.tsv.gz"),
        &table,
    )
    .unwrap();
    let replaced = fs::metadata(&table).unwrap().len();
    assert_ne!(replaced, built_from, "the two tables differ in length");
    let refused = run(dir, "replaced-table", &corpus, &value, "CADD_PHRED");
    assert_eq!(refused.code, Some(0));
    assert_eq!(
        refused.warnings,
        format!(
            "WARNING: Failed to instantiate plugin CADD: tabix index {} not found; a plugin \
             data file needs its .tbi index beside it, so index it with tabix without -C; the \
             store {} beside the table is refused: store source differs: the table {} holds \
             {replaced} bytes and its store's index records {built_from}; {RECONVERT}; \
             {NO_INDEX_CONVERSION}\n",
            tbi.display(),
            store.display(),
            table.display()
        )
    );
    assert_eq!(refused.cells, 0, "the plugin is not loaded");
}

/// The line that ends CADD's load beside a table with no index and a store
/// built on other columns than the table's header names: the missing index,
/// the store, the flags it was built on, the flags that rebuild it and the
/// flag the conversion needs on a table with no index.
fn refused_for_its_columns(table: &Path, store: &Path, built_on: &str, named: &str) -> String {
    format!(
        "WARNING: Failed to instantiate plugin CADD: tabix index {}.tbi not found; a plugin data \
         file needs its .tbi index beside it, so index it with tabix without -C; the store {} \
         beside the table is refused: store columns `{built_on}` differ from the table's \
         `{named}`; {RECONVERT} and the table's flags; {NO_INDEX_CONVERSION}\n",
        table.display(),
        store.display()
    )
}

/// A store of the committed table of single positions, built under
/// `--zero-based`, holds every row one base along. Beside the table's index
/// CADD refuses it for the index's columns and reads the table. With no index
/// the columns are the ones the table's header names, a `Pos` column of
/// one-based positions (CADD.pm:209, 266), so the store is refused by them and
/// the plugin is not loaded, where the store's own flags would score other
/// records than the table does.
#[test]
fn cadd_refuses_a_store_of_shifted_positions_beside_a_table_without_an_index() {
    let corpus = corpus("GRCh38-plugin-cadd-fixture");
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let table = staged_table(
        dir,
        &corpus.join("plugin_data/whole_genome_SNVs_synthetic.tsv.gz"),
    );
    let value = format!("CADD,snv={}", table.display());
    let through_tabix = run(dir, "table", &corpus, &value, "CADD_PHRED");
    assert!(through_tabix.cells > 0, "the table scores its records");

    let store = store_of_the_tables_rows_built_on(&table, None, true);
    let (built_on, named) = (
        "--chr-col 0 --start-col 1 --zero-based",
        "--chr-col 0 --start-col 1",
    );
    let beside_index = run(dir, "with-index", &corpus, &value, "CADD_PHRED");
    assert_eq!(beside_index.code, Some(0));
    assert_eq!(
        beside_index.warnings,
        format!(
            "WARNING: binary store refused, the table is read instead store={} table={} \
             reason=store columns `{built_on}` differ from the table's `{named}`; {RECONVERT} \
             and the table's flags\n",
            store.display(),
            table.display()
        )
    );
    assert_eq!(beside_index.rows, through_tabix.rows);

    fs::remove_file(format!("{}.tbi", table.display())).unwrap();
    let refused = run(dir, "without-index", &corpus, &value, "CADD_PHRED");
    assert_eq!(refused.code, Some(0));
    assert_eq!(
        refused.warnings,
        refused_for_its_columns(&table, &store, built_on, named)
    );
    assert_eq!(refused.cells, 0, "the plugin is not loaded");
}

/// The committed table of spans names a `Start` and an `End` column. The
/// plugin matches a structural variant at its start less one against the
/// `Start` cell (CADD.pm:220, 243, 249, 266), so the cell is a zero-based
/// start, and with no index CADD reads a store of zero-based spans alone: its
/// rows are the rows the table gives through its index, on the corpus's
/// records and on an insertion whose lookup ends at its own POS, which the
/// table does not score. A store of one-based spans, or of points at the
/// `Start` cell, would score that insertion; each is refused by the flags that
/// rebuild it.
#[test]
fn cadd_reads_a_table_of_spans_without_an_index_from_a_store_of_zero_based_spans_alone() {
    let corpus = corpus("GRCh38-plugin-cadd-sv-fixture");
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let table = staged_table(
        dir,
        &corpus.join("plugin_data/1000G_phase3_SVs_synthetic.tsv.gz"),
    );
    let value = format!("CADD,sv={}", table.display());

    // The corpus's insertion alone, with `END=` its own POS added to its INFO.
    let records = fs::read_to_string(corpus.join("variants.vcf")).unwrap();
    let insertion: Vec<String> = records
        .lines()
        .filter_map(|line| {
            if line.starts_with('#') {
                return Some(line.to_owned());
            }
            let mut cells: Vec<&str> = line.split('\t').collect();
            let ended = format!("{};END={}", cells[7], cells[1]);
            (cells[2] == "insertion_sv").then(|| {
                cells[7] = &ended;
                cells.join("\t")
            })
        })
        .collect();
    assert_eq!(insertion.iter().filter(|l| !l.starts_with('#')).count(), 1);
    let one_record = dir.join("one_record.vcf");
    fs::write(&one_record, insertion.join("\n") + "\n").unwrap();
    let inputs = [corpus.join("variants.vcf"), one_record];

    let through_tabix: Vec<Run> = inputs
        .iter()
        .enumerate()
        .map(|(i, input)| {
            run_on(
                dir,
                &format!("table-{i}"),
                input,
                &corpus,
                &value,
                "CADD_PHRED",
            )
        })
        .collect();
    assert!(through_tabix[0].cells > 0, "the table scores its records");
    assert_eq!(
        through_tabix[1].cells, 0,
        "the table gives the insertion that ends at its POS no score"
    );
    fs::remove_file(format!("{}.tbi", table.display())).unwrap();

    let named = "--chr-col 0 --start-col 1 --end-col 2 --zero-based";
    for (end_col, built_on) in [
        (None, "--chr-col 0 --start-col 1"),
        (Some(2), "--chr-col 0 --start-col 1 --end-col 2"),
    ] {
        let store = store_of_the_tables_rows_built_on(&table, end_col, false);
        for (i, (input, table_run)) in inputs.iter().zip(&through_tabix).enumerate() {
            let label = format!("refused-{}-{i}", end_col.is_some());
            let refused = run_on(dir, &label, input, &corpus, &value, "CADD_PHRED");
            assert_eq!(refused.code, Some(0), "{built_on}");
            // The refusal, then the warnings the records draw in any run.
            assert_eq!(
                refused.warnings,
                refused_for_its_columns(&table, &store, built_on, named) + &table_run.warnings
            );
            assert_eq!(refused.cells, 0, "{built_on}: the plugin is not loaded");
        }
    }

    store_of_the_tables_rows_built_on(&table, Some(2), true);
    for (i, (input, table_run)) in inputs.iter().zip(&through_tabix).enumerate() {
        let label = format!("store-{i}");
        let read = run_on(dir, &label, input, &corpus, &value, "CADD_PHRED");
        assert_eq!(read.code, Some(0));
        assert_eq!(
            read.warnings, table_run.warnings,
            "the plugin loads from the store"
        );
        assert_eq!(read.rows, table_run.rows);
    }
}
