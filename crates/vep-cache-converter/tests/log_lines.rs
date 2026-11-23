// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The conversion's log where standard error is not a terminal: every line is
//! plain text, `<timestamp> <level> <target>: <message>`, with no colour code,
//! and a count of one stands beside its noun in the singular.

use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::Command;

/// The whole log of [`conversion`], each line without its timestamp and the
/// space after it.
const EXPECTED: [&str; 9] = [
    " WARN vep_cache_converter: Skipping file with unexpected name: .._1-1000000",
    " INFO vep_cache_converter: transcripts: chr 21 region 1-1000000: 0 entries",
    " INFO vep_cache_converter: Total transcripts: 1 region, 0 entries",
    " WARN vep_cache_converter: Skipping file with unexpected name: .._all_vars.gz",
    " WARN vep_cache_converter: Skipping file with unexpected name: .._all_vars.gz.tbi",
    " INFO vep_cache_converter: Total variations: 0 files on 0 chromosomes",
    " INFO vep_cache_converter: Cache info: species=homo_sapiens, assembly=GRCh37",
    " INFO vep_cache_converter: Copied info.json",
    " INFO vep_cache_converter: Conversion complete: 1 transcript region (0 transcripts), \
     0 variation files on 0 chromosomes, 0 regulatory regions, 2 files in the checksum manifest",
];

/// The whole log of the conversion of [`cache_of_one_chromosome`], each line
/// without its timestamp and the space after it.
const EXPECTED_OF_ONE_CHROMOSOME: [&str; 7] = [
    " INFO vep_cache_converter: transcripts: chr 21 region 1-1000000: 1 entry",
    " INFO vep_cache_converter: Total transcripts: 1 region, 1 entry",
    " INFO vep_cache_converter: variations: chr 21 1-1000000_var.gz",
    " INFO vep_cache_converter: Total variations: 1 file on 1 chromosome",
    " INFO vep_cache_converter: Cache info: species=homo_sapiens, assembly=GRCh37",
    " INFO vep_cache_converter: Copied info.json",
    " INFO vep_cache_converter: Conversion complete: 1 transcript region (1 transcript), \
     1 variation file on 1 chromosome, 0 regulatory regions, 3 files in the checksum manifest",
];

/// Whether `text` is a UTC timestamp to the microsecond,
/// `YYYY-MM-DDTHH:MM:SS.ffffffZ`.
fn is_timestamp(text: &str) -> bool {
    text.len() == 27
        && text.bytes().enumerate().all(|(at, byte)| match at {
            4 | 7 => byte == b'-',
            10 => byte == b'T',
            13 | 16 => byte == b':',
            19 => byte == b'.',
            26 => byte == b'Z',
            _ => byte.is_ascii_digit(),
        })
}

/// The committed flat cache.
fn committed_cache() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/parent_named_shards")
}

/// The conversion of the flat cache `cache` into `out`, with `NO_COLOR`
/// unset.
fn conversion_of(cache: &Path, out: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"));
    command
        .arg("--json-dir")
        .arg(cache)
        .arg("--output-dir")
        .arg(out)
        .env_remove("NO_COLOR");
    command
}

/// The conversion of the committed flat cache into `out`, with `NO_COLOR`
/// unset.
fn conversion(out: &Path) -> Command {
    conversion_of(&committed_cache(), out)
}

/// Writes at `cache` a made-up flat cache of chromosome 21 alone: the
/// committed cache's `info.json`, one region file holding one entry and one
/// variation region file. The conversion counts a region file's entries,
/// reads each for a translation without its `seq_edits` list alone, and
/// copies a variation file byte for byte, so an empty object stands for the
/// transcript and an empty file for the variation file.
fn cache_of_one_chromosome(cache: &Path) {
    fs::create_dir_all(cache.join("transcripts")).unwrap();
    fs::create_dir_all(cache.join("variations/21")).unwrap();
    fs::copy(committed_cache().join("info.json"), cache.join("info.json")).unwrap();
    fs::write(cache.join("transcripts/21_1-1000000.json"), "[{}]\n").unwrap();
    fs::write(cache.join("variations/21/1-1000000_var.gz"), "").unwrap();
}

/// Asserts that `log` is the whole log `expected`: each line a timestamp, one
/// space and the literal rest of the line.
fn assert_plain_log(log: &str, expected: &[&str]) {
    let rests: Vec<&str> = log
        .lines()
        .map(|line| {
            let (stamp, rest) = line.split_once(' ').unwrap_or((line, ""));
            assert!(is_timestamp(stamp), "{line:?} opens with no timestamp");
            rest
        })
        .collect();
    assert_eq!(rests, expected, "{log:?}");
}

#[test]
fn a_conversion_log_sent_to_a_file_is_plain_text() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("conversion.log");
    let status = conversion(&dir.path().join("out"))
        .stderr(File::create(&path).unwrap())
        .status()
        .expect("vep-cache-converter runs");
    let log = std::fs::read_to_string(&path).unwrap();
    assert!(status.success(), "{log}");
    assert_plain_log(&log, &EXPECTED);
}

#[test]
fn a_conversion_log_sent_to_a_pipe_is_plain_text() {
    let dir = tempfile::tempdir().unwrap();
    let output = conversion(&dir.path().join("out"))
        .output()
        .expect("vep-cache-converter runs");
    let log = String::from_utf8(output.stderr).unwrap();
    assert!(output.status.success(), "{log}");
    assert_plain_log(&log, &EXPECTED);
}

#[test]
fn a_count_of_one_stands_beside_its_noun_in_the_singular() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("cache");
    cache_of_one_chromosome(&cache);
    let output = conversion_of(&cache, &dir.path().join("out"))
        .output()
        .expect("vep-cache-converter runs");
    let log = String::from_utf8(output.stderr).unwrap();
    assert!(output.status.success(), "{log}");
    assert_plain_log(&log, &EXPECTED_OF_ONE_CHROMOSOME);
}
