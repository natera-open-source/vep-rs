// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Table directories the tests build for themselves, and a runner of the built binary.
//! The header and the rows are spelled out here, apart from the library's own constants.

#![allow(dead_code)]

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

pub const HEADER: &str =
    "genome\tassembly\ttranscripts\tset\tcache_version\tseed\trecords\toutcome\tsha256";
pub const OUTCOMES_HEADER: &str = "outcome\tmeaning";
pub const SETS_HEADER: &str = "set\trelease\tcache_version\tcaches";

pub const UNCLASSIFIED: &str = "differs:unclassified\tThe outputs differ, with no name.";
pub const NO_TRANSCRIPTS: &str = "excluded:no-transcripts\tThe cache holds no transcript.";

pub const BACTERIUM_KEY: &str =
    "bacteria_0_collection___example_bacterium_gca_000000001__ASM1v1__ensembl";
pub const ARCHAEON_KEY: &str =
    "bacteria_0_collection__example_archaeon_gca_000000002__ASM2v1__ensembl";
pub const CEREVISIAE_KEY: &str = "fungi__saccharomyces_cerevisiae__R64-1-1__ensembl";
pub const POMBE_KEY: &str = "fungi__schizosaccharomyces_pombe__ASM294v2__ensembl";
pub const HUMAN_KEY: &str = "vertebrates__homo_sapiens__GRCh38__ensembl";
pub const HUMAN_MERGED_KEY: &str = "vertebrates__homo_sapiens__GRCh38__merged";

/// A made-up digest: one hex digit, 64 times.
pub fn digest_of(digit: char) -> String {
    digit.to_string().repeat(64)
}

/// The two rows of the bacterial file, in key order: a genome whose name begins with an
/// underscore, then an excluded one.
pub fn bacteria_rows() -> Vec<String> {
    vec![
        format!(
            "_example_bacterium_gca_000000001\tASM1v1\tensembl\tbacteria_0_collection\t63\t11\t28\texact\t{}",
            digest_of('1')
        ),
        "example_archaeon_gca_000000002\tASM2v1\tensembl\tbacteria_0_collection\t63\t11\t0\texcluded:no-transcripts\t-"
            .to_owned(),
    ]
}

/// The four rows of the plain file, in key order: two of `fungi`, two of `vertebrates`.
pub fn genome_rows() -> Vec<String> {
    vec![
        format!(
            "saccharomyces_cerevisiae\tR64-1-1\tensembl\tfungi\t63\t11\t30\texact\t{}",
            digest_of('2')
        ),
        format!(
            "schizosaccharomyces_pombe\tASM294v2\tensembl\tfungi\t63\t11\t29\tdiffers:unclassified\t{}",
            digest_of('3')
        ),
        format!(
            "homo_sapiens\tGRCh38\tensembl\tvertebrates\t116\t11\t31\texact\t{}",
            digest_of('4')
        ),
        format!(
            "homo_sapiens\tGRCh38\tmerged\tvertebrates\t116\t11\t31\texact\t{}",
            digest_of('5')
        ),
    ]
}

/// The sets file that matches the six rows.
pub fn set_rows() -> Vec<String> {
    vec![
        "bacteria_0_collection\t63\t63\t2".to_owned(),
        "fungi\t63\t63\t2".to_owned(),
        "vertebrates\t116\t116\t2".to_owned(),
    ]
}

pub fn outcome_rows() -> Vec<String> {
    vec![UNCLASSIFIED.to_owned(), NO_TRANSCRIPTS.to_owned()]
}

/// A header line and `rows`, each ended by a newline.
pub fn text(header: &str, rows: &[String]) -> String {
    let mut text = format!("{header}\n");
    for row in rows {
        text.push_str(row);
        text.push('\n');
    }
    text
}

pub fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

pub fn gunzip(bytes: &[u8]) -> String {
    let mut text = String::new();
    MultiGzDecoder::new(bytes)
        .read_to_string(&mut text)
        .unwrap();
    text
}

/// A fresh directory for the test `name`, under the build's own scratch directory.
pub fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// A table directory for the test `name` with the four files as given.
pub fn table(
    name: &str,
    genomes: &[String],
    bacteria: &[String],
    outcomes: &[String],
    sets: &[String],
) -> PathBuf {
    let dir = scratch(name);
    fs::write(dir.join("genomes.tsv"), text(HEADER, genomes)).unwrap();
    fs::write(
        dir.join("bacteria.tsv.gz"),
        gzip(text(HEADER, bacteria).as_bytes()),
    )
    .unwrap();
    fs::write(dir.join("outcomes.tsv"), text(OUTCOMES_HEADER, outcomes)).unwrap();
    fs::write(dir.join("sets.tsv"), text(SETS_HEADER, sets)).unwrap();
    dir
}

/// The valid table of six rows across three sets.
pub fn valid_table(name: &str) -> PathBuf {
    table(
        name,
        &genome_rows(),
        &bacteria_rows(),
        &outcome_rows(),
        &set_rows(),
    )
}

pub struct Run {
    pub code: i32,
    pub stdout: String,
    pub stderr: String,
}

/// Runs the built binary with `args`.
pub fn tool(args: &[&str]) -> Run {
    let output = Command::new(env!("CARGO_BIN_EXE_vep-genome-table"))
        .args(args)
        .output()
        .unwrap();
    Run {
        code: output.status.code().unwrap(),
        stdout: String::from_utf8(output.stdout).unwrap(),
        stderr: String::from_utf8(output.stderr).unwrap(),
    }
}

/// The one line a failed `check` of `dir` prints for `file` at `line`.
pub fn failure(dir: &Path, file: &str, line: usize, what: &str) -> String {
    format!(
        "ERROR: [vep-genome-table] {}:{line}: {what}\n",
        dir.join(file).display()
    )
}
