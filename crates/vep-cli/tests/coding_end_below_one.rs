// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A cached model whose coding end lies below 1, through the `vep` binary on
//! a synthetic cache.
//!
//! A ready-made Ensembl VEP cache can hold a coding model whose
//! `cdna_coding_end` is a negative number, on the transcript and in its
//! mapper: ensembl core takes the coding end of a model with sequence edits
//! from its `_transl_end` attribute as it stands (`Transcript.pm:1042-1044`).
//! The cache converter writes the value as the cache holds it, a string
//! holding the number. Ensembl VEP reads the cache and annotates every record
//! of the sequence; no coding position maps to such a model
//! (`TranscriptMapper.pm:441`, where every cDNA position lies past the coding
//! end).
//!
//! The made-up cache holds two models of the RefSeq set on one sequence in one
//! shard: the model with the negative coding end, and an ordinary coding
//! model far from it.

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

/// The first base of the model with the negative coding end.
const ODD_START: u64 = 2001;

/// The first base of the ordinary coding model, more than the upstream and
/// downstream distance away from the other model.
const ORDINARY_START: u64 = 20_001;

/// The model with the negative coding end: two forward exons of 30 and 60
/// bases around an intron of 70, a coding start at the fourth base of the
/// transcript and a coding end of -20, held as a string on the transcript and
/// in its mapper. Its sequence fields are the ones ensembl core cuts for such a
/// pair: a 3' UTR of the last 20 bases and a coding sequence that runs from the
/// coding start to 23 bases before the transcript's end.
fn model_with_a_negative_coding_end() -> serde_json::Value {
    let spliced = format!("GGG{}", "GCT".repeat(29));
    assert_eq!(spliced.len(), 90);
    let coding = &spliced[3..spliced.len() - 23];
    json!({
        "stable_id": "NM_000000001.1",
        "version": 1,
        "dbID": 9001,
        "gene_stable_id": "900000001",
        "gene_version": 1,
        "chr": "1",
        "start": ODD_START,
        "end": ODD_START + 159,
        "strand": 1,
        "biotype": "protein_coding",
        "source": "BestRefSeq",
        "is_canonical": 1,
        "exons": [
            {"stable_id": "NM_000000001.1.1", "start": ODD_START, "end": ODD_START + 29,
                "phase": -1, "end_phase": 0},
            {"stable_id": "NM_000000001.1.2", "start": ODD_START + 100,
                "end": ODD_START + 159, "phase": 0, "end_phase": -1},
        ],
        "cdna_coding_start": 4,
        "cdna_coding_end": "-20",
        "translation": {"stable_id": "NP_000000001.1", "version": 1, "dbID": 7001,
            "start": 4, "end": 1},
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "five_prime_utr": &spliced[..3],
            "translateable_seq": coding,
            "three_prime_utr": &spliced[spliced.len() - 20..],
            "peptide": "A".repeat(coding.len() / 3),
            "introns": [{"start": ODD_START + 30, "end": ODD_START + 99}],
            "mapper": {"start_phase": -1, "cdna_coding_start": 4,
                "cdna_coding_end": "-20",
                "pairs": [
                    {"from_start": 1, "from_end": 30,
                        "to_start": ODD_START, "to_end": ODD_START + 29, "ori": 1},
                    {"from_start": 31, "from_end": 90,
                        "to_start": ODD_START + 100, "to_end": ODD_START + 159, "ori": 1},
                ]},
        },
    })
}

/// The ordinary model: one forward exon that is coding sequence throughout, a
/// start codon, eight alanine codons and a stop codon.
fn ordinary_model() -> serde_json::Value {
    let coding = format!("ATG{}TGA", "GCT".repeat(8));
    let length = coding.len() as u64;
    let end = ORDINARY_START + length - 1;
    json!({
        "stable_id": "NM_000000002.1",
        "version": 1,
        "dbID": 9002,
        "gene_stable_id": "900000002",
        "gene_version": 1,
        "chr": "1",
        "start": ORDINARY_START,
        "end": end,
        "strand": 1,
        "biotype": "protein_coding",
        "source": "BestRefSeq",
        "is_canonical": 1,
        "exons": [{"stable_id": "NM_000000002.1.1", "start": ORDINARY_START, "end": end,
            "phase": 0, "end_phase": 0}],
        "cdna_coding_start": 1,
        "cdna_coding_end": length,
        "translation": {"stable_id": "NP_000000002.1", "version": 1, "dbID": 7002,
            "start": 1, "end": length},
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "translateable_seq": coding,
            "peptide": format!("M{}", "A".repeat(8)),
            "introns": [],
            "mapper": {"start_phase": 0, "cdna_coding_start": 1,
                "cdna_coding_end": length,
                "pairs": [{"from_start": 1, "from_end": length,
                    "to_start": ORDINARY_START, "to_end": end, "ori": 1}]},
        },
    })
}

/// Writes under `dir` a cache of the RefSeq set whose one shard holds the two
/// models and the three records, runs the binary under `--refseq`, and returns
/// the default-format rows sorted.
fn annotate(dir: &Path) -> Vec<String> {
    let cache = dir.join("cache");
    let shard = cache.join("transcripts").join("1").join("1-1000000.json");
    fs::create_dir_all(shard.parent().unwrap()).unwrap();
    fs::write(
        &shard,
        json!([model_with_a_negative_coding_end(), ordinary_model()]).to_string(),
    )
    .unwrap();
    fs::write(
        cache.join("info.json"),
        json!({"assembly": "GRCh38", "source_type": "refseq"}).to_string(),
    )
    .unwrap();

    let input = dir.join("records.vcf");
    fs::write(
        &input,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             1\t{}\tbefore_the_model\tA\tG\t.\t.\t.\n\
             1\t{}\tin_its_intron\tA\tG\t.\t.\t.\n\
             1\t{}\tin_the_ordinary_model\tC\tA\t.\t.\t.\n",
            ODD_START - 51,
            ODD_START + 59,
            ORDINARY_START + 4,
        ),
    )
    .unwrap();
    let out = dir.join("records.out");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(&cache)
        .args(["--assembly", "GRCh38", "--offline", "--refseq"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
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

/// The row of the ordinary model: the record changes the second base of codon
/// 2, base 5 of the coding sequence and of the transcript, from cytosine to
/// adenine, so the alanine codon becomes an aspartate codon.
const ORDINARY_ROW: &str = "in_the_ordinary_model\t1:20005\tA\t900000002\tNM_000000002.1\t\
    Transcript\tmissense_variant\t5\t5\t2\tA/D\tgCt/gAt\t-\tIMPACT=MODERATE;STRAND=1";

/// The row of a record 51 bases before the model with the negative coding end.
const BEFORE_ROW: &str = "before_the_model\t1:1950\tG\t900000001\tNM_000000001.1\tTranscript\t\
    upstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=51;STRAND=1";

/// The row of a record inside that model's intron.
const INTRON_ROW: &str = "in_its_intron\t1:2060\tG\t900000001\tNM_000000001.1\tTranscript\t\
    intron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1";

/// A shard that holds a coding model with a negative coding end loads, as
/// Ensembl VEP reads the region file: a record in another model of the
/// sequence gets that model's row, and a record before the model or in its
/// intron gets the row Ensembl VEP writes for the model there, which needs no
/// coding position.
#[test]
fn a_shard_holding_a_model_with_a_negative_coding_end_is_read() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(annotate(dir.path()), [BEFORE_ROW, INTRON_ROW, ORDINARY_ROW]);
}
