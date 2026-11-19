// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary writes HGVSc for a deletion that covers a whole
//! transcript: the `transcript_ablation` row carries the notation clamped onto
//! the transcript's first and last bases.
//! Every expected value is Ensembl VEP 116.2's output for the same record on
//! the release-116 GRCh38 cache and the chromosome 21 reference the corpus
//! carries.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The identifiers corpus: a chromosome 21 cache slice and the bgzipped
/// reference over its records.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The data lines of a `--hgvs` run on `records`.
fn annotate(dir: &Path, name: &str, records: &str) -> Vec<String> {
    let input = dir.join(format!("{name}.vcf"));
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join(format!("{name}.out"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus().join("json_cache"))
        .arg("--fasta")
        .arg(corpus().join("reference.fa.gz"))
        .args(["--assembly", "GRCh38", "--offline", "--hgvs"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// `(Consequence, cDNA_position, Extra)` of `record`'s row on `transcript`.
fn row_cells(rows: &[String], record: &str, transcript: &str) -> [String; 3] {
    let row = rows
        .iter()
        .find(|row| row.starts_with(record) && row.contains(transcript))
        .unwrap_or_else(|| panic!("no row of {record} on {transcript}: {rows:?}"));
    let c: Vec<&str> = row.split('\t').collect();
    [c[6], c[7], c[13]].map(str::to_string)
}

/// A deletion covering a whole transcript is `within_feature`
/// (BaseVariationFeatureOverlapAllele.pm:480), so Ensembl VEP 116.2 asks
/// `hgvs_transcript` for its `transcript_ablation` row (OutputFactory.pm:1721,
/// 1750-1754): `_var2transcript_slice_coords` refuses only a span with both
/// ends outside on one side and clamps each end onto the transcript
/// (TranscriptVariationAllele.pm:2729-2733), so at a zero 3' shift the row
/// reads `n.1_<length>del`, or `c.-N_*Mdel` on a coding transcript, and no
/// HGVSp (`hgvs_protein` finds no peptide position, :1730-1738). The three
/// records delete the snRNA ENST00000384745 (forward, 107 bases), the miRNA
/// ENST00000581787 (reverse, 76 bases) and the single-exon coding
/// ENST00000286808 (reverse, 1,241 bases) whole, each with a first deleted
/// base that differs from the base after the deletion, so nothing shifts.
#[test]
fn a_deletion_covering_a_whole_transcript_writes_hgvsc_on_its_ablation_row() {
    let dir = tempfile::tempdir().unwrap();
    let records = concat!(
        "21\t13621568\tsnrna_fwd\t",
        "TAGGGCTTTGTGCCTGCTTTGGCGGCACATATCCTAAAATTGGAACAATACAGAGAAAGTTAGCATGG",
        "CTTCTGCATAAGGAGGCAGCACAACTCTTTGAGGCATTCCATATTTTGTG",
        "\tT\t.\t.\t.\n",
        "21\t13644766\tmirna_rev\t",
        "GCTGTAAATCACACTAGAAGAATTTTCATAATGCAGTCACATGGTGATTGTGTGTGATTGCATTATGA",
        "ATATTATTGTATGTGTGGGC",
        "\tG\t.\t.\t.\n",
        "21\t30165560\tcoding_rev\t",
        "TTACATTTTAGAATAGGTATATGTTTTATTGTTATAGAATATCAATTTCAACTTTAGATTCTTATCAG",
        "ATAATACATGATGGTCAAAAAGTAATAAGTAGTCCATCATTGTACACTTGCTAATGATTATGGGTATA",
        "CTTGAATTTAGAATATGAAAATTGGAGCAGCGATTGGGTCCAACAATAGTCAATGCATTTAATACACC",
        "TGAAGTCAGAATGTCCAAGTAAGTCATTTTCATTCCTACCACCAGTAACAAATTGAAACTTTCGTATC",
        "AATGTAAATCTAGACATAAAGCAAGTTCTCCTGCCTCACATACTTACAGTTTCTAGCAGGACTTTATA",
        "AAAAACATTGACCATAGTCCATACTTGGAGCCAAAAGGAGGCATTAGACATAACTGGTGGAGGTCTTA",
        "CTAAGCATTGTCGTATTTCTTCGCTTATCTGTGTGTGGCACACGGTAGCCAGGCACTGGATATCTGTA",
        "CCCTTGCTTCTTTCTGTTGCAGCAGCAAAATCCACAAAGCAGACCCCCTCCAATGAAGAGGACAGCAG",
        "CGCTTGCCCAGCCAAGGAAAAGTGCTGCTCCCAGCTCTCGTTTCTGACCTATGTGGATGGCTGGGTTG",
        "TAGAAATCTCTGATGATTATATTGGCTGTCCAGCTCACCGGAATCAGAACGAAGATGCCCGTCAGGAT",
        "GAAGAGGACTCCTGAAGTTCCCAGAAGGTATGCTTTGGCCCTCTCGTTAGAGCCTGTGCACTGGACCT",
        "GCTTCATGCCACAGATGCCAATAAGCAGGGCGATCAAGGAGAGAGCAACAGCCACACACATGAGGGCC",
        "CGGGCTGTTTCCAGGGCAGGCGGGAGAGCCAACAAGGAGCTATAGAACTTGCATTGCAACCGGACCCT",
        "GGCTTGTCGGATGCAATTCATCCAGAGCCCTTCCCAGAGCCTCTCAAAGACAATAATGTTGCTGCCAA",
        "CAAAAGCTGATACTCTCCACTGAGGCAGAAGGGTTGTGGCAAGAGTCCCCACCATGCCAAGGAACCCA",
        "AGAACCAGCCCAGCAATTTGCAAGGGATAAAATGCCATTGCCTTTACTTTGAAGACTGGTTCAATTCG",
        "GAGTGGTAGAAGATGCTCAAGATGTAGAACGTGGAAGCCTTTTGAGTCTAGAGAAGTCTTTGATGATG",
        "TTCTCCTTGGCTGCTGCCCATAACGTTTTAGTCCAAATCAGTAGCTGTGACTGAACTTGGCCTAACTA",
        "GAAGTACCTGTTGTAAATGCATGTTG",
        "\tT\t.\t.\t.\n",
    );
    let rows = annotate(dir.path(), "whole", records);
    assert_eq!(
        row_cells(&rows, "snrna_fwd", "ENST00000384745"),
        [
            "transcript_ablation",
            "-",
            "IMPACT=HIGH;STRAND=1;HGVSc=ENST00000384745.1:n.1_107del",
        ]
    );
    assert_eq!(
        row_cells(&rows, "mirna_rev", "ENST00000581787"),
        [
            "transcript_ablation",
            "-",
            "IMPACT=HIGH;STRAND=-1;HGVSc=ENST00000581787.1:n.1_76del",
        ]
    );
    assert_eq!(
        row_cells(&rows, "coding_rev", "ENST00000286808"),
        [
            "transcript_ablation",
            "-",
            "IMPACT=HIGH;STRAND=-1;HGVSc=ENST00000286808.5:c.-188_*378del",
        ]
    );
    // The same deletion inside a neighbouring transcript is clamped by neither
    // end and reads as any intronic deletion.
    assert_eq!(
        row_cells(&rows, "snrna_fwd", "ENST00000299443"),
        [
            "intron_variant",
            "-",
            "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000299443.6:c.1056-1236_1056-1120del",
        ]
    );
}
