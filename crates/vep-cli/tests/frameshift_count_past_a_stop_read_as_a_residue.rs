// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The count to the new stop in the protein notation of a frameshift whose
//! first changed residue lies after a stop codon that the transcript's
//! peptide reads as a residue, through the `vep` binary on a made-up cache
//! and reference: the first shape of the documented difference
//! `hgvsp_frameshift_count_past_a_stop_read_as_a_residue`.
//!
//! The model is one exon: six bases of 5' UTR, a coding sequence of 24 codons
//! and nine bases of 3' UTR. Codon 10 is a stop codon that the peptide reads
//! as a selenocysteine, through the `_selenocysteine` sequence edit of its
//! cache record. Three records each delete one or two coding bases. Ensembl
//! VEP 116.3 counts to the first stop of the alternate sequence's
//! translation, which is codon 10 for a change after it
//! (`TranscriptVariationAllele.pm:2492`, `:2503`, `:2508` and `:2519`), and
//! writes the unknown count on both records after that codon.

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

const CHROMOSOME: &str = "1";

/// The genomic position of the model's first base.
const MODEL_START: usize = 101;

const UTR5: &str = "GTCACA";

/// The 3' UTR: its third to fifth bases are the same stop codon as codon 10.
const UTR3: [&str; 3] = ["CA", "TGA", "TTTT"];

/// Codons 21 to 24 of the coding sequence: a threonine, a leucine and a
/// serine codon and the stop codon that closes it. Read one base on, the
/// bases of the first three and the next spell a proline codon and a stop
/// codon that is another codon than codon 10.
const LAST_CODONS: [&str; 4] = ["ACC", "CTA", "AGC", "TAG"];

/// The residue number of the stop codon that the peptide reads as a residue.
const SELENOCYSTEINE: usize = 10;

/// The model's coding sequence: a start codon, eight alanine codons, the
/// stop codon read as a selenocysteine, ten alanine codons and the four
/// codons of `LAST_CODONS`: 72 bases.
fn coding() -> String {
    format!(
        "ATG{}TGA{}{}",
        "GCT".repeat(8),
        "GCT".repeat(10),
        LAST_CODONS.concat()
    )
}

/// The model's peptide: 23 residues, the tenth a selenocysteine.
fn peptide() -> String {
    format!("M{}U{}TLS", "A".repeat(8), "A".repeat(10))
}

fn cdna() -> String {
    format!("{UTR5}{}{}", coding(), UTR3.concat())
}

/// The genomic position of coding base `coding`.
fn genomic(coding: usize) -> usize {
    MODEL_START + UTR5.len() + coding - 1
}

/// The model's cache record, with the sequence edit that makes residue 10 a
/// selenocysteine.
fn model() -> serde_json::Value {
    let cdna = cdna();
    let end = MODEL_START + cdna.len() - 1;
    let coding_start = UTR5.len() + 1;
    let coding_end = UTR5.len() + coding().len();
    json!({
        "stable_id": "ENST00000000001",
        "version": 1,
        "gene_stable_id": "ENSG00000000001",
        "chr": CHROMOSOME,
        "start": MODEL_START,
        "end": end,
        "strand": 1,
        "biotype": "protein_coding",
        "source": "ensembl",
        "protein": "ENSP00000000001.1",
        "exons": [{"start": MODEL_START, "end": end, "rank": 1, "phase": -1, "end_phase": -1}],
        "cdna_coding_start": coding_start,
        "cdna_coding_end": coding_end,
        "attributes": [],
        "translation": {"stable_id": "ENSP00000000001", "version": 1, "start": coding_start,
            "end": coding_end},
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "five_prime_utr": UTR5,
            "three_prime_utr": UTR3.concat(),
            "spliced_seq": cdna,
            "translateable_seq": coding(),
            "peptide": peptide(),
            "seq_edits": [{"code": "_selenocysteine", "name": "Selenocysteine",
                "start": SELENOCYSTEINE, "end": SELENOCYSTEINE, "alt_seq": "U",
                "length_diff": 0}],
            "mapper": {"start_phase": -1, "cdna_coding_start": coding_start,
                "cdna_coding_end": coding_end,
                "pairs": [{"from_start": 1, "from_end": cdna.len(), "ori": 1,
                    "to_start": MODEL_START, "to_end": end}]},
        },
    })
}

/// The made-up chromosome: the model's bases at their place, `C` elsewhere.
fn chromosome() -> Vec<u8> {
    let cdna = cdna();
    let mut bases = vec![b'C'; 300];
    let from = MODEL_START - 1;
    bases[from..from + cdna.len()].copy_from_slice(cdna.as_bytes());
    bases
}

/// A record that deletes the coding bases `first` to `last`, written from
/// the base before them.
fn deletion(name: &str, first: usize, last: usize) -> String {
    let coding = coding();
    let reference = &coding[first - 2..last];
    format!(
        "{CHROMOSOME}\t{}\t{name}\t{reference}\t{}\t.\t.\t.\n",
        genomic(first - 1),
        &reference[..1]
    )
}

/// The lines of the default output of `vep --hgvs` on the three records,
/// from its column line on.
fn rows(dir: &Path) -> Vec<String> {
    let cache = dir.join("cache");
    let shard = cache
        .join("transcripts")
        .join(CHROMOSOME)
        .join("1-1000000.json");
    fs::create_dir_all(shard.parent().unwrap()).unwrap();
    fs::write(&shard, json!([model()]).to_string()).unwrap();
    fs::write(
        cache.join("info.json"),
        json!({"assembly": "GRCh38"}).to_string(),
    )
    .unwrap();

    let bases = chromosome();
    let fasta = dir.join("reference.fa");
    let header = format!(">{CHROMOSOME}\n");
    let mut text = header.clone().into_bytes();
    text.extend_from_slice(&bases);
    text.push(b'\n');
    fs::write(&fasta, text).unwrap();
    fs::write(
        dir.join("reference.fa.fai"),
        format!(
            "{CHROMOSOME}\t{len}\t{offset}\t{len}\t{width}\n",
            len = bases.len(),
            offset = header.len(),
            width = bases.len() + 1
        ),
    )
    .unwrap();

    let input = dir.join("records.vcf");
    fs::write(
        &input,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}{}{}",
            deletion("before_the_selenocysteine_codon", 13, 13),
            deletion("to_a_stop_of_another_codon", 43, 43),
            deletion("to_the_selenocysteine_codon", 43, 44),
        ),
    )
    .unwrap();
    let out = dir.join("records.txt");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(&cache)
        .arg("--fasta")
        .arg(&fasta)
        .args(["--assembly", "GRCh38", "--offline"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .arg("--hgvs")
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with("##"))
        .map(str::to_string)
        .collect()
}

/// The first record deletes coding base 13, the first base of codon 5,
/// before the selenocysteine codon: the new frame reads that codon's bases
/// across two codons and meets its first stop at residue 22, 18 residues
/// counted from the first changed one, as both engines write it.
///
/// The second deletes coding base 43, the first base of codon 15. The new
/// frame reads leucine codons from the alanine codons 15 to 19, a leucine
/// codon over the end of codon 20, a proline codon and, at residue 22, a stop
/// codon that is another codon than codon 10: eight residues counted from
/// the first changed one. Ensembl VEP's alternate translation has its first
/// stop at residue 10, 10 - 14 is not positive, and it writes
/// `p.Ala15LeufsTer?`.
///
/// The third deletes coding bases 43 and 44. The new frame's first stop, at
/// residue 25, is the third to fifth bases of the 3' UTR, the same stop codon
/// as the one the peptide reads as a selenocysteine. Whether that codon ends
/// the protein there is not known from the model, and both engines write the
/// unknown count.
#[test]
fn a_frameshift_past_a_stop_read_as_a_residue_counts_to_a_stop_of_another_codon() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows(dir.path()),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "before_the_selenocysteine_codon\t1:119\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t19\t13\t5\tA/X\tGct/ct\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.13del;HGVSp=ENSP00000000001.1:p.Ala5LeufsTer18",
            "to_a_stop_of_another_codon\t1:149\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t49\t43\t15\tA/X\tGct/ct\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.43del;HGVSp=ENSP00000000001.1:p.Ala15LeufsTer8",
            "to_the_selenocysteine_codon\t1:149-150\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t49-50\t43-44\t15\tA/X\tGCt/t\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.43_44del;HGVSp=ENSP00000000001.1:p.Ala15CysfsTer?",
        ]
    );
}
