// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` file and a cache's variation file written as two bgzip streams
//! joined, read through the `vep` binary. The first stream's end marker, a
//! block that holds no bytes, stands inside such a file, and `tabix` gives the
//! row that opens the second stream that block's offset as its chunk start
//! (htslib `bgzf_read_block` skips empty blocks), so a query at that row
//! enters the file at the empty block.

mod common;

use std::path::{Path, PathBuf};
use std::process::Output;

/// The flags every run passes.
const FLAGS: &[&str] = &[
    "--offline",
    "--species",
    "homo_sapiens",
    "--assembly",
    "GRCh38",
    "--buffer_size",
    "5000",
    "--force_overwrite",
    "--no_stats",
    "--quiet",
    "--no_headers",
];

/// The flags of the co-located corpus's `colocated` flag set.
const COLOCATED_FLAGS: &[&str] = &[
    "--check_existing",
    "--af",
    "--af_1kg",
    "--af_gnomade",
    "--af_gnomadg",
    "--max_af",
    "--pubmed",
    "--var_synonyms",
    "--clinvar_somatic_classification",
    "--failed",
    "1",
];

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{name}"))
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// Copies the directory `from` and everything under it to `to`.
fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// Runs `vep` on `input` against `cache` with [`FLAGS`] and `more`; the
/// default-format output follows the process.
fn run(out_dir: &Path, input: &Path, cache: &Path, more: &[&str]) -> (Output, String) {
    let out = out_dir.join("out.txt");
    let output = common::isolated_vep(out_dir)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache)
        .args(FLAGS)
        .args(more)
        .output()
        .unwrap();
    let written = std::fs::read_to_string(&out).unwrap_or_default();
    (output, written)
}

/// The score file of `tests/fixtures/spliceai-two-streams` as a `--custom`
/// VCF source: the record at the row that opens the second stream carries the
/// source's field with that row's value, as the records in either stream
/// carry theirs.
#[test]
fn a_custom_vcf_of_two_joined_streams_is_read_at_the_row_that_opens_the_second() {
    let fixture = fixture("spliceai-two-streams");
    let corpus = corpus("GRCh38-plugin-spliceai");
    let tmp = tempfile::tempdir().unwrap();
    let custom = format!(
        "file={},short_name=TwoStreams,format=vcf,type=exact,coords=0,fields=SpliceAI",
        fixture
            .join("scores_two_streams.synthetic.vcf.gz")
            .display()
    );
    let fasta = corpus.join("reference.fa.gz");
    let (output, written) = run(
        tmp.path(),
        &fixture.join("variants.vcf"),
        &corpus.join("json_cache"),
        &["--fasta", fasta.to_str().unwrap(), "--custom", &custom],
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let rows = |record: &str| -> Vec<&str> {
        written
            .lines()
            .filter(|row| row.split('\t').next() == Some(record))
            .collect()
    };
    let opening = rows("opens_second_stream");
    assert_eq!(opening.len(), 1);
    assert!(
        opening[0].contains("TwoStreams_SpliceAI=G|TMPRSS15|0.47|0.45|1.00|0.35|-36|26|28|20"),
        "{}",
        opening[0]
    );
    for (record, count) in [
        ("first_stream_one_gene", 2),
        ("first_stream_two_genes", 4),
        ("second_stream_two_alternates", 12),
    ] {
        let rows = rows(record);
        assert_eq!(rows.len(), count, "{record}");
        assert!(
            rows.iter().all(|row| row.contains("TwoStreams_SpliceAI=")),
            "{record}"
        );
    }
}

/// `tests/fixtures/variation-two-streams` in place of the co-located corpus's
/// variation file: the corpus's input under its `colocated` flags reads the
/// known variants at the row that opens the second stream and after it, so
/// the output is the one the corpus's own cache gives.
#[test]
fn a_variation_file_of_two_joined_streams_is_read_at_the_row_that_opens_the_second() {
    let corpus = corpus("GRCh38-colocated");
    let input = corpus.join("variants.vcf");
    let tmp = tempfile::tempdir().unwrap();
    let (output, expected) = run(
        &tmp.path().join("own"),
        &input,
        &corpus.join("json_cache"),
        COLOCATED_FLAGS,
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(expected.contains("rs1555830109"));

    let cache = tmp.path().join("cache");
    copy_dir(
        &corpus.join("json_cache/transcripts"),
        &cache.join("transcripts"),
    );
    std::fs::copy(corpus.join("json_cache/info.json"), cache.join("info.json")).unwrap();
    let variations = cache.join("variations/21");
    std::fs::create_dir_all(&variations).unwrap();
    for name in ["all_vars.gz", "all_vars.gz.tbi"] {
        std::fs::copy(
            fixture("variation-two-streams").join(name),
            variations.join(name),
        )
        .unwrap();
    }
    let (output, actual) = run(&tmp.path().join("joined"), &input, &cache, COLOCATED_FLAGS);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert!(
        actual == expected,
        "{}",
        common::compare_headerless(&expected, &actual)
    );
}
