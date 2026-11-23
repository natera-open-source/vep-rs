// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A translation the cache holds without a stable id, through the `vep`
//! binary on a synthetic cache.
//!
//! A ready-made Ensembl VEP cache can hold a translation that carries a
//! database id and no stable id: in the RefSeq transcript set of several
//! species the models with the source `ensembl` do, T-cell receptor segments
//! and models of the biotype `cDNA` among them. The cache converter writes
//! such a translation without the `stable_id` key. Ensembl VEP reads a cached
//! translation's stable id in one place, the reference name of an HGVSp string
//! (`TranscriptVariationAllele.pm:1748` at ensembl-variation release/116), so
//! on the default columns such a model is annotated like any other coding
//! model.
//!
//! The segment model's transcript carries the source `ensembl` and a stable id
//! that is no RefSeq accession, so a run on the RefSeq set drops it unless
//! `--all_refseq` is given (`AnnotationType/Transcript.pm:236-246` at
//! release/116.3).

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

/// The first genomic base of the segment model's one exon.
const SEGMENT_START: u64 = 2001;

/// The whole codons of the segment's coding sequence; one more base follows
/// them, which the exon's end phase of 1 leaves for the next segment.
const SEGMENT_CODONS: u64 = 15;

/// The codon the record changes, counted from 1.
const CHANGED_CODON: u64 = 5;

/// The segment model: one forward exon that is coding sequence throughout,
/// fifteen alanine codons and one base, with a translation that carries a
/// database id and no stable id, and no protein id on the transcript.
fn segment_model() -> serde_json::Value {
    let coding = format!("{}G", "GCT".repeat(SEGMENT_CODONS as usize));
    let length = coding.len() as u64;
    let end = SEGMENT_START + length - 1;
    json!({
        "stable_id": "SEGJ1.1",
        "version": 1,
        "dbID": 9001,
        "gene_stable_id": "SEGJ1",
        "gene_version": 1,
        "chr": "1",
        "start": SEGMENT_START,
        "end": end,
        "strand": 1,
        "biotype": "TR_J_gene",
        "source": "ensembl",
        "is_canonical": 1,
        "exons": [{"stable_id": "SEGJ1.1.1", "start": SEGMENT_START, "end": end,
            "phase": 0, "end_phase": 1}],
        "cdna_coding_start": 1,
        "cdna_coding_end": length,
        "translation": {"dbID": 7001, "start": 1, "end": length},
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "translateable_seq": coding,
            "peptide": "A".repeat(SEGMENT_CODONS as usize),
            "introns": [],
            "mapper": {"start_phase": 0, "cdna_coding_start": 1,
                "cdna_coding_end": length,
                "pairs": [{"from_start": 1, "from_end": length,
                    "to_start": SEGMENT_START, "to_end": end, "ori": 1}]},
        },
    })
}

/// A model named by a RefSeq accession whose one intron holds the segment.
fn accession_model() -> serde_json::Value {
    json!({
        "stable_id": "XR_000000001.1",
        "version": 1,
        "dbID": 9002,
        "gene_stable_id": "900000001",
        "gene_version": 1,
        "chr": "1",
        "start": 1001,
        "end": 3000,
        "strand": 1,
        "biotype": "lncRNA",
        "source": "Gnomon",
        "is_canonical": 1,
        "exons": [
            {"stable_id": "XR_000000001.1.1", "start": 1001, "end": 1100,
                "phase": -1, "end_phase": -1},
            {"stable_id": "XR_000000001.1.2", "start": 2901, "end": 3000,
                "phase": -1, "end_phase": -1},
        ],
        "variation_effect_feature_cache": {
            "introns": [{"start": 1101, "end": 2900}],
            "mapper": {"pairs": [
                {"from_start": 1, "from_end": 100,
                    "to_start": 1001, "to_end": 1100, "ori": 1},
                {"from_start": 101, "from_end": 200,
                    "to_start": 2901, "to_end": 3000, "ori": 1},
            ]},
        },
    })
}

/// The genomic position of the second base of [`CHANGED_CODON`].
fn record_position() -> u64 {
    SEGMENT_START + 3 * (CHANGED_CODON - 1) + 1
}

/// Writes under `dir` a cache of the RefSeq set that holds the two models and
/// one record that changes the second base of the segment's fifth codon, runs
/// the binary with `set_flags`, and returns the default-format rows sorted.
fn annotate(dir: &Path, set_flags: &[&str]) -> Vec<String> {
    let cache = dir.join("cache");
    let shard = cache.join("transcripts").join("1").join("1-1000000.json");
    fs::create_dir_all(shard.parent().unwrap()).unwrap();
    fs::write(
        &shard,
        json!([segment_model(), accession_model()]).to_string(),
    )
    .unwrap();
    fs::write(
        cache.join("info.json"),
        json!({"assembly": "GRCh38", "source_type": "refseq"}).to_string(),
    )
    .unwrap();

    let input = dir.join("record.vcf");
    fs::write(
        &input,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             1\t{}\tsegment_codon\tC\tA\t.\t.\t.\n",
            record_position()
        ),
    )
    .unwrap();
    let out = dir.join("record.out");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(&cache)
        .args(["--assembly", "GRCh38", "--offline"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .args(set_flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{set_flags:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut rows: Vec<String> = fs::read_to_string(&out)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    rows.sort();
    rows
}

/// The row of the segment model: the record changes the second base of codon
/// 5, base 14 of the coding sequence and of the transcript, from cytosine to
/// adenine, so the alanine codon becomes an aspartate codon.
const SEGMENT_ROW: &str = "segment_codon\t1:2014\tA\tSEGJ1\tSEGJ1.1\tTranscript\t\
    missense_variant\t14\t14\t5\tA/D\tgCt/gAt\t-\tIMPACT=MODERATE;STRAND=1";

/// The row of the accession model, whose intron holds the record.
const ACCESSION_ROW: &str = "segment_codon\t1:2014\tA\t900000001\tXR_000000001.1\tTranscript\t\
    intron_variant,non_coding_transcript_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1";

/// A cache whose shard holds a translation without a stable id loads, and
/// with every model of the RefSeq set kept the segment model's row carries the
/// coding consequence, positions, amino acids and codons of its own coding
/// sequence, beside the row of the model named by an accession.
#[test]
fn a_model_whose_translation_has_no_stable_id_is_annotated_with_all_refseq() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        annotate(dir.path(), &["--refseq", "--all_refseq"]),
        [ACCESSION_ROW, SEGMENT_ROW]
    );
}

/// On the RefSeq set without `--all_refseq` the same cache loads and the
/// record's one row is that of the model named by an accession.
#[test]
fn the_refseq_set_alone_reads_the_shard_and_drops_the_segment_model() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(annotate(dir.path(), &["--refseq"]), [ACCESSION_ROW]);
}
