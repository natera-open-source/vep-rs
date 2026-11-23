// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The names of the columns of a Parquet part in any letter case, read from
//! the parts the binary writes over three committed corpora: no two columns of
//! a part differ in letter case alone; a tab column whose name an earlier
//! column carries in other letters is stored behind `tab_`; and the footer's
//! `tab_columns` key states the `--tab` column names, which the round trip to
//! the tab format prints.
#![cfg(feature = "parquet")]

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

use parquet::file::reader::FileReader;
use parquet::file::serialized_reader::SerializedFileReader;
use vep_io::output::parquet::{
    canonical_tab_text, part_files, tab_column_types, tab_rows, KEY_COLUMNS,
};

mod common;

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

/// The names of a part's columns, in order, and the value of its footer's
/// `tab_columns` key.
fn names_and_tab_columns(part: &Path) -> (Vec<String>, Option<String>) {
    let reader = SerializedFileReader::new(fs::File::open(part).unwrap()).unwrap();
    let footer = reader.metadata().file_metadata();
    let names = footer
        .schema()
        .get_fields()
        .iter()
        .map(|field| field.name().to_string())
        .collect();
    let stated = footer
        .key_value_metadata()
        .into_iter()
        .flatten()
        .find(|entry| entry.key == "tab_columns")
        .and_then(|entry| entry.value.clone());
    (names, stated)
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

/// Holds every part of `out` to the six key columns, then `stored`, with no
/// two names equal in any letter case and the tab columns of `tab` under the
/// footer's `tab_columns` key; then holds the output read back through the
/// cell contract to the tab output's column names and rows.
fn assert_part_names(out: &Path, tab: &str, stored: &[String], label: &str) {
    let columns = tab_columns(tab);
    let parts = part_files(out).unwrap();
    assert!(!parts.is_empty(), "{label}: no part");
    for part in &parts {
        let (names, stated) = names_and_tab_columns(part);
        let part = part.display();
        assert_eq!(names[..KEY_COLUMNS.len()], KEY_COLUMNS, "{label} {part}");
        assert_eq!(names[KEY_COLUMNS.len()..], *stored, "{label} {part}");
        let distinct: BTreeSet<String> = names.iter().map(|name| name.to_lowercase()).collect();
        assert_eq!(distinct.len(), names.len(), "{label} {part}");
        assert_eq!(
            stated.as_deref(),
            Some(columns.join("\t").as_str()),
            "{label} {part}"
        );
    }

    let types = tab_column_types(out).unwrap();
    let read_back: Vec<&str> = types.iter().map(|(name, _)| name.as_str()).collect();
    assert_eq!(read_back, columns, "{label}");
    let mut expected: Vec<String> = tab
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            line.split('\t')
                .zip(&types)
                .map(|(cell, (_, column_type))| canonical_tab_text(*column_type, cell))
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect();
    expected.sort_unstable();
    let mut rows = tab_rows(out).unwrap();
    rows.sort_unstable();
    assert!(
        rows == expected,
        "{label}: the rows read back differ from the tab output's"
    );
}

/// A `--custom` source with the short name `REF` declares the columns `REF`,
/// `REF_AF` and `REF_CLASS`. `REF` is the key column `ref` to a reader that
/// compares names without regard to case, so the part stores it as `tab_REF`;
/// the two other columns and every column of a run that names no such source
/// keep their names.
#[test]
fn a_custom_source_named_as_a_key_column_in_capitals_is_stored_behind_the_prefix() {
    let corpora = common::corpora();
    let corpus = corpus(&corpora, "GRCh38-custom");
    let run = common::default_run(corpus);
    let tmp = tempfile::tempdir().unwrap();
    let source = format!(
        "file={},short_name=REF,format=vcf,type=exact,fields=AF%CLASS",
        corpus.dir.join("custom/sites.vcf.gz").display()
    );
    let extra = ["--custom", source.as_str()];

    let tab = common::run_vep(&run, "tab", tmp.path(), &extra);
    let columns = tab_columns(&tab);
    let place = columns.iter().position(|column| column == "REF").unwrap();
    assert_eq!(columns[place..], ["REF", "REF_AF", "REF_CLASS"]);
    let mut stored = columns.clone();
    stored[place] = "tab_REF".to_string();
    let out = write_parquet(&run, tmp.path(), "custom.parquet", &extra);
    assert_part_names(&out, &tab, &stored, "a source named REF");

    let tab = common::run_vep(&run, "tab", tmp.path(), &[]);
    let out = write_parquet(&run, tmp.path(), "plain.parquet", &[]);
    assert_part_names(&out, &tab, &tab_columns(&tab), "no source");
}

/// The CADD plugin declares `CADD_PHRED` and `CADD_RAW`, and the dbNSFP plugin
/// given the columns `CADD_phred` and `CADD_raw` declares those two: loaded
/// together, over the made-up tables of the two plugins' fixture corpora, they
/// give four tab columns that differ in letter case alone, pair by pair. The
/// part stores the second of each pair behind `tab_`.
#[test]
fn two_plugins_columns_that_differ_in_letter_case_alone_are_stored_apart() {
    let corpora = common::corpora();
    let dbnsfp = corpus(&corpora, "GRCh38-plugin-dbnsfp-fixture");
    let cadd = corpus(&corpora, "GRCh38-plugin-cadd-fixture");
    let run = common::default_run(dbnsfp);
    let tmp = tempfile::tempdir().unwrap();
    let cadd_plugin = format!(
        "CADD,snv={}",
        cadd.dir
            .join("plugin_data/whole_genome_SNVs_synthetic.tsv.gz")
            .display()
    );
    let dbnsfp_plugin = format!(
        "dbNSFP,{},CADD_phred,CADD_raw",
        dbnsfp
            .dir
            .join("plugin_data/dbNSFP4.9a_grch38_synthetic.gz")
            .display()
    );
    let extra = [
        "--plugin",
        cadd_plugin.as_str(),
        "--plugin",
        dbnsfp_plugin.as_str(),
    ];

    let tab = common::run_vep(&run, "tab", tmp.path(), &extra);
    let columns = tab_columns(&tab);
    let place = columns.len() - 4;
    assert_eq!(
        columns[place..],
        ["CADD_PHRED", "CADD_RAW", "CADD_phred", "CADD_raw"]
    );
    let mut stored = columns.clone();
    stored[place + 2] = "tab_CADD_phred".to_string();
    stored[place + 3] = "tab_CADD_raw".to_string();
    let out = write_parquet(&run, tmp.path(), "plugins.parquet", &extra);
    assert_part_names(&out, &tab, &stored, "CADD beside dbNSFP");
}

/// The format comparison's README says under which name its Parquet reader
/// prints a column the writer stored behind `tab_`.
#[test]
fn the_format_comparison_readme_says_where_a_stored_column_takes_its_name() {
    let readme = fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/vep-formats-compare/README.md"),
    )
    .unwrap();
    let text = readme.split_whitespace().collect::<Vec<_>>().join(" ");
    let sentence = "the columns after them are the tab columns, one the writer stored behind \
                    `tab_` (a name an earlier column carries in any letter case) under the name \
                    the file's footer states under `tab_columns`,";
    assert!(
        text.contains(sentence),
        "tools/vep-formats-compare/README.md: {sentence}"
    );
}

/// The Parquet page states the rule the two tests above hold the parts to,
/// names the footer key and says which name the round trip to the tab format
/// prints.
#[test]
fn the_parquet_page_says_how_names_that_differ_in_letter_case_are_stored() {
    let page = fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/output/parquet.md"),
    )
    .unwrap();
    for sentence in [
        "No two columns share a name, in any letter case.",
        "A tab column whose name an earlier column carries in other letters is stored behind \
         `tab_` too, as many times as give a name no earlier column carries.",
        "A `--custom` source named `REF` is stored as `tab_REF`, and the `CADD_phred` column of \
         the committed synthetic dbNSFP fixture, loaded after the CADD plugin's `CADD_PHRED`, \
         as `tab_CADD_phred`.",
        "DuckDB and Spark compare column names without regard to case.",
        "The footer's [`tab_columns`](#footer-metadata) key states the `--tab` name of every \
         column, and the [round trip](#round-trip-to---tab) prints each column under it.",
        "| `tab_columns`                           | the `--tab` column names in column order, \
         joined by tabs: the name each column after the six key columns prints under",
    ] {
        assert!(
            page.contains(sentence),
            "docs/output/parquet.md: {sentence}"
        );
    }
    assert!(
        !page.contains("No two columns share a name. "),
        "docs/output/parquet.md states the rule for names equal letter for letter alone"
    );
}
