// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The Parquet type of a column that a plugin or a `--custom` source declares
//! under the name of a typed column, read from the parts the binary writes
//! over two committed corpora: a string holding the `--tab` cell, whatever
//! the name, while the column a flag alone gives keeps the type of its name.
#![cfg(feature = "parquet")]

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

use vep_io::output::parquet::{
    canonical_tab_text, column_type, part_files, tab_column_name, tab_rows, KEY_COLUMNS,
};

mod common;

/// The tab columns of an output's first part, in order, each under the name the
/// `--tab` format gives it (`tab_column_name`) with the type its footer stores
/// it under: a primitive by its physical type, a group as a list or a struct.
fn stored_types(out: &Path) -> Vec<(String, &'static str)> {
    use parquet::basic::{ConvertedType, Type as Physical};
    use parquet::file::reader::FileReader;
    let parts = part_files(out).unwrap();
    let reader = parquet::file::serialized_reader::SerializedFileReader::new(
        fs::File::open(&parts[0]).unwrap(),
    )
    .unwrap();
    reader
        .metadata()
        .file_metadata()
        .schema()
        .get_fields()
        .iter()
        .skip(KEY_COLUMNS.len())
        .map(|field| {
            let stored = if field.is_group() {
                match field.get_basic_info().converted_type() {
                    ConvertedType::LIST => "LIST",
                    _ => "STRUCT",
                }
            } else {
                match field.get_physical_type() {
                    Physical::BYTE_ARRAY => "STRING",
                    Physical::INT32 => "INT32",
                    Physical::DOUBLE => "DOUBLE",
                    Physical::BOOLEAN => "BOOLEAN",
                    other => panic!("{}: stored as {other}", field.name()),
                }
            };
            (tab_column_name(field.name()).to_string(), stored)
        })
        .collect()
}

/// The stored type of one tab column of an output.
fn stored_type(out: &Path, column: &str) -> &'static str {
    stored_types(out)
        .into_iter()
        .find(|(name, _)| name == column)
        .unwrap_or_else(|| panic!("{} has no column {column}", out.display()))
        .1
}

/// The columns of a `--tab` output, from its column line.
fn tab_columns(tab: &str) -> Vec<String> {
    tab.lines()
        .find(|l| l.starts_with('#') && !l.starts_with("##"))
        .expect("the tab column line")
        .trim_start_matches('#')
        .split('\t')
        .map(str::to_string)
        .collect()
}

/// The data rows of a `--tab` output as its Parquet output holds them: a cell
/// of a column in `declared` as written, every other cell as the type of its
/// column's name prints it.
fn rows_as_stored(tab: &str, declared: &BTreeSet<String>) -> Vec<String> {
    let columns = tab_columns(tab);
    let mut rows: Vec<String> = tab
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|line| {
            line.split('\t')
                .zip(&columns)
                .map(|(cell, column)| {
                    if declared.contains(column) {
                        cell.to_string()
                    } else {
                        canonical_tab_text(column_type(column), cell)
                    }
                })
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect();
    rows.sort_unstable();
    rows
}

/// The rows of a Parquet output read back through the cell contract, sorted.
fn stored_rows(out: &Path) -> Vec<String> {
    let mut rows = tab_rows(out).unwrap();
    rows.sort_unstable();
    rows
}

/// The cell of `column` on each row of `rows` at `location`, by the row's
/// `Feature`.
fn cells_by_feature(
    rows: &[String],
    columns: &[String],
    location: &str,
    column: &str,
) -> BTreeMap<String, String> {
    let index = |name: &str| columns.iter().position(|c| c == name).unwrap();
    let (at, feature, cell) = (index("Location"), index("Feature"), index(column));
    rows.iter()
        .map(|row| row.split('\t').collect::<Vec<_>>())
        .filter(|cells| cells[at] == location)
        .map(|cells| (cells[feature].to_string(), cells[cell].to_string()))
        .collect()
}

/// The run's Parquet output under `extra`, written to `name` in `out_dir`.
fn write_parquet(run: &common::Run, out_dir: &Path, name: &str, extra: &[&str]) -> PathBuf {
    let path = out_dir.join(name);
    let mut cmd = common::vep_command(run, out_dir);
    cmd.arg("-o")
        .arg(&path)
        .args(["--output_format", "parquet"])
        .args(extra);
    let out = common::run_output(run, &mut cmd);
    assert!(
        out.status.success(),
        "vep --output_format parquet failed on {}:\n{}",
        run.label(),
        String::from_utf8_lossy(&out.stderr)
    );
    path
}

/// The corpus `name` among `corpora`.
fn corpus<'a>(corpora: &'a [common::Corpus], name: &str) -> &'a common::Corpus {
    corpora
        .iter()
        .find(|corpus| corpus.name == name)
        .unwrap_or_else(|| panic!("the corpus {name}"))
}

/// The flag set `dbnsfp-all` of the dbNSFP fixture corpus declares the made-up
/// table's `TSL` column, the name of the INT32 column `--tsl` gives, and its
/// cell on the row of `variants.vcf` line 18 and ENST00000290200 is a
/// `,`-joined list. Every field the plugin declares is a string column holding
/// the `--tab` cell, with `--tsl` and without it; under the flag the rows the
/// plugin left without a value carry the cache's level in the same column.
/// Without the plugin the flag's column is INT32.
#[test]
fn a_plugin_field_named_as_a_typed_column_is_stored_as_the_string_the_plugin_wrote() {
    const RECORD: &str = "21:33296341";
    const PLUGIN_ROW: &str = "ENST00000290200";
    const OTHER_ROW: &str = "ENST00000422891";
    const FILE_TSL: &str = "5,1,2,3,2,NA";

    let corpora = common::corpora();
    let corpus = corpus(&corpora, "GRCh38-plugin-dbnsfp-fixture");
    let runs = common::runs(corpus);
    let plugin_run = runs
        .iter()
        .find(|run| run.name == "dbnsfp-all")
        .expect("the flag set dbnsfp-all");
    let plain_run = common::default_run(corpus);
    let tmp = tempfile::tempdir().unwrap();

    let plain_columns: BTreeSet<String> =
        tab_columns(&common::run_vep(&plain_run, "tab", tmp.path(), &[]))
            .into_iter()
            .collect();
    let declared: BTreeSet<String> =
        tab_columns(&common::run_vep(plugin_run, "tab", tmp.path(), &[]))
            .into_iter()
            .filter(|column| !plain_columns.contains(column))
            .collect();
    assert!(declared.contains("TSL") && declared.contains("APPRIS"));

    for (extra, flag_level) in [(&[][..], "-"), (&["--tsl"][..], "3")] {
        let label = format!("{} {extra:?}", plugin_run.label());
        let tab = common::run_vep(plugin_run, "tab", tmp.path(), extra);
        let columns = tab_columns(&tab);
        let out = write_parquet(plugin_run, tmp.path(), "plugin.parquet", extra);

        let expected_types: Vec<(String, &str)> = columns
            .iter()
            .map(|column| {
                let stored = if declared.contains(column) {
                    "STRING"
                } else {
                    column_type(column).stored_type_name()
                };
                (column.clone(), stored)
            })
            .collect();
        assert_eq!(stored_types(&out), expected_types, "{label}");
        assert_eq!(stored_type(&out, "TSL"), "STRING", "{label}");
        assert_eq!(stored_type(&out, "STRAND"), "INT32", "{label}");

        let rows = stored_rows(&out);
        assert_eq!(rows, rows_as_stored(&tab, &declared), "{label}");
        let levels = cells_by_feature(&rows, &columns, RECORD, "TSL");
        assert_eq!(levels[PLUGIN_ROW], FILE_TSL, "{label}");
        assert_eq!(levels[OTHER_ROW], flag_level, "{label}");
    }

    let out = write_parquet(&plain_run, tmp.path(), "plain.parquet", &["--tsl"]);
    assert_eq!(stored_type(&out, "TSL"), "INT32");
}

/// A `--custom` VCF source with the short name `gnomADg` and the field `AF`
/// declares the key `gnomADg_AF`, the name of the frequency column
/// `--af_gnomadg` gives. The source's keys are string columns holding the
/// `--tab` cells, with the flag and without it; without the source the flag's
/// column is a list.
#[test]
fn a_custom_key_that_spells_a_frequency_column_is_a_string_column() {
    let corpora = common::corpora();
    let corpus = corpus(&corpora, "GRCh38-custom");
    let run = common::default_run(corpus);
    let tmp = tempfile::tempdir().unwrap();
    let source = format!(
        "file={},short_name=gnomADg,format=vcf,type=exact,fields=AF%CLASS",
        corpus.dir.join("custom/sites.vcf.gz").display()
    );
    let declared: BTreeSet<String> = ["gnomADg", "gnomADg_AF", "gnomADg_CLASS"]
        .map(str::to_string)
        .into();

    for flags in [&[][..], &["--af_gnomadg"][..]] {
        let extra: Vec<&str> = ["--custom", source.as_str()]
            .into_iter()
            .chain(flags.iter().copied())
            .collect();
        let tab = common::run_vep(&run, "tab", tmp.path(), &extra);
        let out = write_parquet(&run, tmp.path(), "custom.parquet", &extra);
        for key in &declared {
            assert_eq!(stored_type(&out, key), "STRING", "{key} {flags:?}");
        }
        assert_eq!(
            stored_rows(&out),
            rows_as_stored(&tab, &declared),
            "{flags:?}"
        );
    }

    let out = write_parquet(&run, tmp.path(), "flag.parquet", &["--af_gnomadg"]);
    assert_eq!(stored_type(&out, "gnomADg_AF"), "LIST");
}

/// The Parquet page states the rule the two tests above hold the parts to.
#[test]
fn the_parquet_page_says_a_declared_field_is_a_string_whatever_its_name() {
    let page = fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/output/parquet.md"),
    )
    .unwrap();
    for sentence in [
        "A field that a plugin or a `--custom` source declares is a STRING whatever its name.",
        "A plugin's `TSL` field is a string column with `--tsl` and without it, \
         and so is a `--custom` key that spells a frequency column, such as `gnomADg_AF`.",
    ] {
        assert!(
            page.contains(sentence),
            "docs/output/parquet.md: {sentence}"
        );
    }
}
