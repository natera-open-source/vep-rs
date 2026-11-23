// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the test files share: a scratch directory, the lines of a fabricated tab output
//! and a run of the binary over one.

#![allow(dead_code)]

use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

pub const BIN: &str = env!("CARGO_BIN_EXE_vep-fastvep-plugin-tuples");

/// The field map a run without `--field-map` reads.
pub const FIELD_MAP: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/data/fastvep_field_map.tsv");

/// The header line of a field map.
pub const MAP_HEADER: &str =
    "fastvep_column\tpart\tvep_field\tsource\tkind\tnote\tlist_separator\treference_set";

/// What opens every refusal.
pub const ERROR: &str = "ERROR: [fastvep_plugin_tuples] ";

/// The end of a summary line that counts nothing left out.
pub const QUIET: &str = "values outside --fields: 0; parts not written word(score): 0; \
                         cells with several entries: 0; unmapped columns: none";

/// The built-in columns of the fabricated tab output.
pub const BUILT_IN: [&str; 14] = [
    "#Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "Extra",
];

/// A directory of the system's temporary directory, removed when the test ends.
pub struct Scratch(pub PathBuf);

impl Scratch {
    pub fn new() -> Scratch {
        // Tests run side by side, so a directory's name is its own.
        static DIRECTORIES: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "vep-fastvep-plugin-tuples-{}-{}",
            std::process::id(),
            DIRECTORIES.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    pub fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    /// Writes a file of the directory and returns its path as text.
    pub fn write(&self, name: &str, bytes: impl AsRef<[u8]>) -> String {
        let path = self.path(name);
        fs::write(&path, bytes).unwrap();
        text(&path)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// A path of the tests as text.
pub fn text(path: &Path) -> String {
    path.to_str().expect("a test path is UTF-8").to_string()
}

/// One run of the binary over an input of a scratch directory.
pub struct Run {
    pub code: i32,
    pub stdout: Vec<u8>,
    pub stderr: String,
    pub tuple_file: PathBuf,
}

impl Run {
    /// The tuple file's lines, each split at its tabs.
    pub fn tuples(&self) -> Vec<Vec<String>> {
        let text = String::from_utf8(self.tuple_bytes()).unwrap();
        let mut lines: Vec<&str> = text.split('\n').collect();
        assert_eq!(lines.pop(), Some(""), "the tuple file ends at a line feed");
        lines
            .iter()
            .map(|line| line.split('\t').map(str::to_string).collect())
            .collect()
    }

    pub fn tuple_bytes(&self) -> Vec<u8> {
        fs::read(&self.tuple_file).unwrap()
    }

    pub fn wrote(&self) -> bool {
        self.tuple_file.exists()
    }
}

/// One row of fourteen built-in cells under gene `ENSG1`, then the supplementary cells.
pub fn row(location: &str, allele: &str, feature: &str, cells: &[&str]) -> String {
    row_of_gene("ENSG1", location, allele, feature, cells)
}

pub fn row_of_gene(
    gene: &str,
    location: &str,
    allele: &str,
    feature: &str,
    cells: &[&str],
) -> String {
    let mut all = vec![
        "v",
        location,
        allele,
        gene,
        feature,
        "Transcript",
        "missense_variant",
        "10",
        "10",
        "4",
        "A/T",
        "gcc/Acc",
        "-",
        "IMPACT=MODERATE",
    ];
    all.extend_from_slice(cells);
    all.join("\t")
}

/// The header line over the built-in columns and `columns`, then the rows.
pub fn table(columns: &[&str], rows: &[String]) -> String {
    let mut text = [&BUILT_IN[..], columns].concat().join("\t");
    text.push('\n');
    for row in rows {
        text.push_str(row);
        text.push('\n');
    }
    text
}

/// The binary run with `args`.
pub fn binary<S: AsRef<OsStr>>(args: &[S]) -> Output {
    Command::new(BIN).args(args).output().unwrap()
}

pub fn exit_code(output: &Output) -> i32 {
    output.status.code().expect("the binary exits with a code")
}

/// The binary over `bytes` written as `raw.tab` of the directory, the tuple file
/// `tuples.tsv` beside it, then the options.
pub fn run_bytes(scratch: &Scratch, bytes: &[u8], options: &[&str]) -> Run {
    let raw = scratch.write("raw.tab", bytes);
    let tuple_file = scratch.path("tuples.tsv");
    let mut args = vec![raw, text(&tuple_file)];
    args.extend(options.iter().map(|option| option.to_string()));
    let output = binary(&args);
    Run {
        code: exit_code(&output),
        stdout: output.stdout,
        stderr: String::from_utf8(output.stderr).unwrap(),
        tuple_file,
    }
}

pub fn run_text(scratch: &Scratch, text: &str, options: &[&str]) -> Run {
    run_bytes(scratch, text.as_bytes(), options)
}

pub fn run(scratch: &Scratch, columns: &[&str], rows: &[String], options: &[&str]) -> Run {
    run_text(scratch, &table(columns, rows), options)
}

/// The binary over one row at 21:1000, allele G, feature ENST1, with one supplementary
/// column holding `cell`.
pub fn one_cell(scratch: &Scratch, column: &str, cell: &str, options: &[&str]) -> Run {
    run(
        scratch,
        &[column],
        &[row("21:1000", "G", "ENST1", &[cell])],
        options,
    )
}

/// The tuples of `pairs` at one location, allele and feature.
pub fn at_place(
    pairs: &[(&str, &str)],
    location: &str,
    allele: &str,
    feature: &str,
) -> Vec<Vec<String>> {
    pairs
        .iter()
        .map(|(field, value)| {
            [location, allele, feature, field, value]
                .iter()
                .map(|cell| cell.to_string())
                .collect()
        })
        .collect()
}

/// The tuples of `pairs` at 21:1000, allele G, feature ENST1.
pub fn at(pairs: &[(&str, &str)]) -> Vec<Vec<String>> {
    at_place(pairs, "21:1000", "G", "ENST1")
}

/// The tuples of `pairs` at 21:1000, allele G and `feature`.
pub fn at_feature(pairs: &[(&str, &str)], feature: &str) -> Vec<Vec<String>> {
    at_place(pairs, "21:1000", "G", feature)
}

/// One tuple as its five cells.
pub fn tuple(cells: [&str; 5]) -> Vec<String> {
    cells.iter().map(|cell| cell.to_string()).collect()
}

/// A map of the scratch directory: the header line, then `rows`.
pub fn custom_map(scratch: &Scratch, rows: &str) -> String {
    scratch.write("map.tsv", format!("{MAP_HEADER}\n{rows}"))
}
