// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--transcript_filter` on the two keys a RefSeq or merged cache's transcripts
//! carry and an Ensembl cache's lack: `_gene_symbol_id`, the EntrezGene
//! identifier the symbol of a RefSeq gene is taken from
//! (`AnnotationSource/Database/Transcript.pm:563` at release/116.2), and
//! `_bam_edit_status`, the outcome of editing a model against the cache's
//! alignment file (`AnnotationType/Transcript.pm:601, 621`). Ensembl VEP's
//! filter reads both as exact keys of the transcript hash (`FilterSet.pm:538`)
//! and finds a key the transcript lacks undefined, so `exists` and a unit with
//! a value drop the transcript and `nex` keeps it (`FilterSet.pm:481-483,
//! 704-711`). Run through the `vep` binary on the GRCh38 RefSeq golden corpus
//! and on the GRCh38 transcript-filter corpus, in the four output formats.

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// A release 116 golden corpus.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{name}"))
}

/// The record within reach of five RefSeq models of two genes: `NM_003225.3` of
/// `TFF1` (EntrezGene 7031) and four predicted `XR_` models of `LOC105372815`
/// (EntrezGene 105372815).
const TFF1: &str = "21\t42362510\t21:42362510:G:T\tG\tT\t.\t.\t.";

/// The record within reach of `NR_109925.1` alone, a model the cache's
/// alignment file edited (`bam_edit_status` `ok`).
const EDITED: &str = "21\t18737924\trs975049326\tC\tT\t.\t.\t.";

/// The chromosome 21 record of the transcript-filter corpus within reach of 41
/// Ensembl transcripts.
const ENSEMBL: &str = "21\t17608184\ttf_17608184\tA\tG\t.\t.\t.";

const FIVE: [&str; 5] = [
    "NM_003225.3",
    "XR_007067875.1",
    "XR_007067876.1",
    "XR_007067877.1",
    "XR_937755.3",
];

/// The output formats and the flag that selects each.
const FORMATS: [(&str, &[&str]); 4] = [
    ("default", &[]),
    ("tab", &["--tab"]),
    ("vcf", &["--vcf"]),
    ("json", &["--json"]),
];

/// Runs the binary on `record` with the corpus cache flags `cache` and
/// `filter`, in format `format`, and returns the output file's text.
fn annotate(dir: &Path, record: &str, cache: &[&str], format: &str, filter: &[&str]) -> String {
    let input = dir.join("in.vcf");
    fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{record}\n"),
    )
    .unwrap();
    let out = dir.join(format!("out.{format}"));
    let (_, format_flags) = FORMATS.iter().find(|(name, _)| *name == format).unwrap();
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args([
            "--offline",
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh38",
        ])
        .args(cache)
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .args(*format_flags)
        .args(filter)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{filter:?} {format}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read_to_string(&out).unwrap()
}

/// The models of the rows of `text` in `format`: the `Feature` column of a
/// default or tab row (`-` on an intergenic row), the `Feature` field of each
/// CSQ entry, or the `transcript_id` of each JSON transcript consequence, with
/// `-` for a JSON line holding intergenic consequences alone.
fn models(text: &str, format: &str) -> Vec<String> {
    match format {
        "default" | "tab" => text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').nth(4).unwrap().to_string())
            .collect(),
        "vcf" => {
            let header = text
                .lines()
                .find(|line| line.starts_with("##INFO=<ID=CSQ"))
                .unwrap();
            let fields: Vec<&str> = header
                .split("Format: ")
                .nth(1)
                .unwrap()
                .trim_end_matches("\">")
                .split('|')
                .collect();
            let feature = fields.iter().position(|f| *f == "Feature").unwrap();
            text.lines()
                .filter(|line| !line.starts_with('#'))
                .flat_map(|line| {
                    let info = line.split('\t').nth(7).unwrap();
                    let csq = info
                        .split(';')
                        .find_map(|item| item.strip_prefix("CSQ="))
                        .unwrap();
                    csq.split(',')
                        .map(|entry| {
                            let value = entry.split('|').nth(feature).unwrap();
                            if value.is_empty() { "-" } else { value }.to_string()
                        })
                        .collect::<Vec<_>>()
                })
                .collect()
        }
        "json" => text
            .lines()
            .flat_map(|line| {
                let record: serde_json::Value = serde_json::from_str(line).unwrap();
                match record.get("transcript_consequences") {
                    Some(serde_json::Value::Array(entries)) => entries
                        .iter()
                        .map(|entry| entry["transcript_id"].as_str().unwrap().to_string())
                        .collect::<Vec<_>>(),
                    _ => {
                        assert!(record.get("intergenic_consequences").is_some(), "{line}");
                        vec!["-".to_string()]
                    }
                }
            })
            .collect(),
        other => panic!("{other}"),
    }
}

/// The GRCh38 RefSeq corpus cache under `--refseq`, with its reference.
fn refseq_cache() -> Vec<String> {
    let corpus = corpus("GRCh38-refseq");
    vec![
        "--cache".to_string(),
        "--dir_cache".to_string(),
        corpus.join("cache_root").display().to_string(),
        "--cache_version".to_string(),
        "116".to_string(),
        "--refseq".to_string(),
        "--fasta".to_string(),
        corpus.join("reference.fa.gz").display().to_string(),
    ]
}

/// `_gene_symbol_id` reads the EntrezGene identifier on every model of the
/// RefSeq cache: `exists` keeps the five models, `is 7031` the TFF1 model
/// alone, `gene_symbol_id` reaches the key through the `_` rule
/// (`FilterSet.pm:562-564`), and `_gene_symbol_i`, a prefix of the one key,
/// keeps the five models too, where Ensembl VEP 116.2 keeps four: its prefix
/// branch reads nothing on the first transcript, the documented divergence
/// `transcript_filter_abbreviated_fields`
/// (`docs/divergences/transcript-filter-fields.md`, "An abbreviated field on
/// the first transcript"); a numeric comparison reads the identifier as a
/// number. The kept rows are the unfiltered run's own rows.
#[test]
fn gene_symbol_id_is_read_on_a_refseq_cache() {
    let dir = tempfile::tempdir().unwrap();
    let cache = refseq_cache();
    let cache: Vec<&str> = cache.iter().map(String::as_str).collect();
    for (format, _) in FORMATS {
        let plain = annotate(dir.path(), TFF1, &cache, format, &[]);
        assert_eq!(models(&plain, format), FIVE, "{format}");
        let kept = annotate(
            dir.path(),
            TFF1,
            &cache,
            format,
            &["--transcript_filter", "_gene_symbol_id exists"],
        );
        assert_eq!(models(&kept, format), FIVE, "{format}");
        fn rows(text: &str) -> Vec<&str> {
            text.lines().filter(|line| !line.starts_with('#')).collect()
        }
        assert_eq!(rows(&kept), rows(&plain), "{format}");
        for filter in ["_gene_symbol_id is 7031", "gene_symbol_id < 8000"] {
            let one = annotate(
                dir.path(),
                TFF1,
                &cache,
                format,
                &["--transcript_filter", filter],
            );
            assert_eq!(models(&one, format), ["NM_003225.3"], "{format} {filter}");
        }
        for filter in ["_gene_symbol_i exists", "gene_symbol_id > 7000"] {
            let all = annotate(
                dir.path(),
                TFF1,
                &cache,
                format,
                &["--transcript_filter", filter],
            );
            assert_eq!(models(&all, format), FIVE, "{format} {filter}");
        }
        let none = annotate(
            dir.path(),
            TFF1,
            &cache,
            format,
            &["--transcript_filter", "_gene_symbol_id nex"],
        );
        assert_eq!(models(&none, format), ["-"], "{format}");
    }
}

/// `_bam_edit_status` reads the status the cache stores on an edited model:
/// `exists` and `is ok` keep `NR_109925.1`, whose `BAM_EDIT` column prints
/// `OK`; `is OK` drops it, since `eq` compares the stored lower-case text
/// (`FilterSet.pm:704`), while `match OK` keeps it (714, case-insensitive);
/// `nex` drops it.
#[test]
fn bam_edit_status_is_read_on_an_edited_model() {
    let dir = tempfile::tempdir().unwrap();
    let cache = refseq_cache();
    let cache: Vec<&str> = cache.iter().map(String::as_str).collect();
    for (format, _) in FORMATS {
        let plain = annotate(dir.path(), EDITED, &cache, format, &[]);
        assert_eq!(models(&plain, format), ["NR_109925.1"], "{format}");
        if format == "default" {
            assert!(plain.contains("BAM_EDIT=OK"), "{plain}");
        }
        for filter in [
            "_bam_edit_status exists",
            "_bam_edit_status is ok",
            "bam_edit_status match OK",
            "_b exists",
        ] {
            let kept = annotate(
                dir.path(),
                EDITED,
                &cache,
                format,
                &["--transcript_filter", filter],
            );
            assert_eq!(models(&kept, format), ["NR_109925.1"], "{format} {filter}");
        }
        for filter in ["_bam_edit_status is OK", "_bam_edit_status nex"] {
            let dropped = annotate(
                dir.path(),
                EDITED,
                &cache,
                format,
                &["--transcript_filter", filter],
            );
            assert_eq!(models(&dropped, format), ["-"], "{format} {filter}");
        }
    }
}

/// An Ensembl cache's transcripts carry neither key, so a filter naming one
/// reads it undefined, as Ensembl VEP reads a key the transcript hash lacks:
/// `exists` and a unit with a value drop every transcript, leaving one
/// `intergenic_variant` row, and `nex` or a negated unit keeps all 41.
#[test]
fn an_ensembl_cache_leaves_both_keys_undefined() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus("GRCh38-transcript-filter");
    let json_cache = corpus.join("json_cache").display().to_string();
    let cache = ["--json_cache", json_cache.as_str()];
    for (format, _) in FORMATS {
        let plain = annotate(dir.path(), ENSEMBL, &cache, format, &[]);
        assert_eq!(models(&plain, format).len(), 41, "{format}");
        for filter in [
            "_gene_symbol_id exists",
            "_gene_symbol_id is 7031",
            "_bam_edit_status exists",
            "_bam_edit_status is ok",
        ] {
            let dropped = annotate(
                dir.path(),
                ENSEMBL,
                &cache,
                format,
                &["--transcript_filter", filter],
            );
            assert_eq!(models(&dropped, format), ["-"], "{format} {filter}");
        }
        for filter in [
            "_gene_symbol_id nex",
            "not _gene_symbol_id is 7031",
            "_bam_edit_status nex",
        ] {
            let kept = annotate(
                dir.path(),
                ENSEMBL,
                &cache,
                format,
                &["--transcript_filter", filter],
            );
            assert_eq!(
                models(&kept, format),
                models(&plain, format),
                "{format} {filter}"
            );
        }
    }
}
