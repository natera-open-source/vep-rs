// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A cache transcript whose gene has no stable id, through the `vep` binary.
//!
//! A RefSeq cache of a species other than human can hold cDNA alignments whose
//! gene carries no stable id: the transcript has its own `stable_id` and an
//! undefined `_gene_stable_id`, and the converted cache writes its record
//! without the `gene_stable_id` key. Ensembl VEP keeps such a transcript (its
//! filter and its merge step test the transcript's own id alone,
//! `AnnotationType/Transcript.pm:220` and `:280`), annotates it like any other
//! and copies the undefined value into the row's gene
//! (`OutputFactory.pm:1477`): the default and tab formats print `-`, the VCF
//! format leaves the field empty and the JSON format writes no `gene_id` key.
//!
//! The cache here is one such record of the release 116 rabbit RefSeq cache
//! (`oryctolagus_cuniculus_refseq`, OryCun2.0) as the converter writes it, and
//! each expected row is the row Ensembl VEP 116.3 writes for that transcript
//! on the one record, on that cache under the same flags.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

use serde_json::json;

mod common;

/// The transcript's record: a four-exon cDNA alignment on the reverse strand
/// with no translation and no `gene_stable_id` key.
fn transcript_record() -> serde_json::Value {
    let exons = [
        (1661834, 1662176),
        (1647614, 1647766),
        (1645713, 1646183),
        (1642373, 1643588),
    ];
    let exon = |(start, end): (u64, u64)| json!({"start": start, "end": end, "phase": -1, "end_phase": -1});
    json!({
        "stable_id": "NM_001205004.1",
        "dbID": 93202,
        "start": 1642373,
        "end": 1662176,
        "strand": -1,
        "biotype": "cDNA",
        "source": "ensembl",
        "is_canonical": 1,
        "ccds": "-",
        "gene_phenotype": 0,
        "attributes": [],
        "exons": exons.iter().copied().map(exon).collect::<Vec<_>>(),
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "translateable_seq": "",
            "seq_edits": [],
            "introns": [
                {"start": 1647767, "end": 1661833},
                {"start": 1646184, "end": 1647613},
                {"start": 1643589, "end": 1645712},
            ],
            "sorted_exons": exons.iter().rev().copied().map(exon).collect::<Vec<_>>(),
            "mapper": {
                "pair_count": 4,
                "start_phase": -1,
                "pairs": [
                    {"from_start": 1, "from_end": 343, "to_start": 1661834, "to_end": 1662176, "ori": -1},
                    {"from_start": 344, "from_end": 496, "to_start": 1647614, "to_end": 1647766, "ori": -1},
                    {"from_start": 497, "from_end": 967, "to_start": 1645713, "to_end": 1646183, "ori": -1},
                    {"from_start": 968, "from_end": 2183, "to_start": 1642373, "to_end": 1643588, "ori": -1},
                ],
            },
        },
    })
}

/// Writes the one-transcript RefSeq cache and the one record under `dir`;
/// the cache and the input.
fn stage(dir: &Path) -> (PathBuf, PathBuf) {
    let cache = dir.join("cache");
    let shard = cache.join("transcripts/1/1000001-2000000.json");
    fs::create_dir_all(shard.parent().unwrap()).unwrap();
    fs::write(&shard, json!([transcript_record()]).to_string()).unwrap();
    fs::write(
        cache.join("info.json"),
        json!({
            "species": "oryctolagus_cuniculus",
            "assembly": "OryCun2.0",
            "cache_version": 116,
            "source_type": "refseq",
        })
        .to_string(),
    )
    .unwrap();
    let input = dir.join("one.vcf");
    fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t1645800\tno_gene_id_1\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    (cache, input)
}

/// Runs `vep --refseq` on the record with `flags` and returns the output's
/// data lines (those not opening with `#`).
fn rows(flags: &[&str]) -> Vec<String> {
    let dir = tempfile::tempdir().unwrap();
    let (cache, input) = stage(dir.path());
    let out = dir.path().join("out.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["--offline", "--json_cache"])
        .arg(&cache)
        .args([
            "--species",
            "oryctolagus_cuniculus",
            "--assembly",
            "OryCun2.0",
            "--refseq",
            "--no_stats",
            "--force_overwrite",
        ])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{flags:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The default format prints `-` for the gene and the transcript's own row.
#[test]
fn default_format_prints_a_dash_for_the_gene() {
    assert_eq!(
        rows(&[]),
        [
            "no_gene_id_1\t1:1645800\tG\t-\tNM_001205004.1\tTranscript\t\
             non_coding_transcript_exon_variant\t880\t-\t-\t-\t-\t-\t\
             IMPACT=MODIFIER;STRAND=-1"
        ]
    );
}

/// The tab format prints `-` in its `Gene` column.
#[test]
fn tab_format_prints_a_dash_for_the_gene() {
    assert_eq!(
        rows(&["--tab"]),
        [
            "no_gene_id_1\t1:1645800\tG\t-\tNM_001205004.1\tTranscript\t\
             non_coding_transcript_exon_variant\t880\t-\t-\t-\t-\t-\t\
             MODIFIER\t-\t-1\t-\t-\t-"
        ]
    );
}

/// The VCF format leaves the `Gene` field of the `CSQ` entry empty.
#[test]
fn vcf_format_leaves_the_gene_field_empty() {
    let rows = rows(&["--vcf"]);
    assert_eq!(rows.len(), 1, "{rows:?}");
    assert_eq!(
        rows[0].split('\t').nth(7).unwrap(),
        "CSQ=G|non_coding_transcript_exon_variant|MODIFIER|||Transcript|NM_001205004.1|cDNA|\
         3/4||||880|||||||-1|||||"
    );
}

/// The JSON format writes the transcript's consequence without a `gene_id`
/// key.
#[test]
fn json_format_writes_no_gene_id_key() {
    let rows = rows(&["--json"]);
    assert_eq!(rows.len(), 1, "{rows:?}");
    let record: serde_json::Value = serde_json::from_str(&rows[0]).unwrap();
    assert_eq!(
        record["transcript_consequences"],
        json!([{
            "cdna_end": 880,
            "cdna_start": 880,
            "consequence_terms": ["non_coding_transcript_exon_variant"],
            "impact": "MODIFIER",
            "strand": -1,
            "transcript_id": "NM_001205004.1",
            "variant_allele": "G",
        }])
    );
}
