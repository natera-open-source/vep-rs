// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation of a frameshift that the end of the transcript cuts
//! short, through the `vep` binary on a made-up cache and reference: the
//! `HGVSp` pair of the documented difference
//! `stop_codon_overlap_without_stop_codon`.
//!
//! Two models of one exon each hold the same eight whole codons, the last of
//! which names a residue, and end on their last coding base, with no 3' UTR.
//! The first carries no flag; the second is flagged `cds_end_NF`. A deletion
//! that shifts the frame in the last codon leaves fewer than three bases from
//! the first changed codon on, so the alternate translation ends before the
//! first changed residue. Ensembl VEP 116.3 then words the row's reference
//! residues as deleted (`TranscriptVariationAllele.pm:2342-2346`) and, on the
//! model without the flag, where its lost-stop test holds
//! (`VariationEffect.pm:1373-1381`), adds an extension of unknown length
//! (`TranscriptVariationAllele.pm:1916-1936`).

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

const CHROMOSOME: &str = "1";

/// The genomic position of the first base of the model without a flag.
const UNFLAGGED_START: usize = 101;

/// The genomic position of the first base of the model flagged `cds_end_NF`,
/// out of the other model's reach.
const FLAGGED_START: usize = 10_201;

const UTR5: &str = "GTCACA";

/// The coding sequence of both models: a start codon, five alanine codons, a
/// proline codon at coding bases 19 to 21 and an aspartate codon at 22 to 24.
fn eight_codons() -> String {
    format!("ATG{}CCGGAC", "GCT".repeat(5))
}

/// The protein of both models: eight residues.
const PEPTIDE: &str = "MAAAAAPD";

/// The genomic position of coding base `coding` of the model that starts at
/// `model_start`.
fn genomic(model_start: usize, coding: usize) -> usize {
    model_start + UTR5.len() + coding - 1
}

/// A cache record of one exon that ends on the last base of `coding`: the 5'
/// UTR, the coding sequence and no 3' UTR, carrying `flags`.
fn model(number: usize, start: usize, coding: &str, flags: &[&str]) -> serde_json::Value {
    let cdna = format!("{UTR5}{coding}");
    let end = start + cdna.len() - 1;
    let coding_start = UTR5.len() + 1;
    let coding_end = cdna.len();
    let attributes: Vec<serde_json::Value> = flags
        .iter()
        .map(|flag| json!({"code": flag, "value": "1"}))
        .collect();
    json!({
        "stable_id": format!("ENST0000000000{number}"),
        "version": 1,
        "gene_stable_id": format!("ENSG0000000000{number}"),
        "chr": CHROMOSOME,
        "start": start,
        "end": end,
        "strand": 1,
        "biotype": "protein_coding",
        "source": "ensembl",
        "protein": format!("ENSP0000000000{number}.1"),
        "exons": [{"start": start, "end": end, "rank": 1, "phase": -1, "end_phase": -1}],
        "cdna_coding_start": coding_start,
        "cdna_coding_end": coding_end,
        "attributes": attributes,
        "translation": {"stable_id": format!("ENSP0000000000{number}"), "version": 1,
            "start": coding_start, "end": coding_end},
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "five_prime_utr": UTR5,
            "spliced_seq": cdna,
            "translateable_seq": coding,
            "peptide": PEPTIDE,
            "seq_edits": [],
            "mapper": {"start_phase": -1, "cdna_coding_start": coding_start,
                "cdna_coding_end": coding_end,
                "pairs": [{"from_start": 1, "from_end": cdna.len(), "ori": 1,
                    "to_start": start, "to_end": end}]},
        },
    })
}

/// The made-up chromosome: each model's bases at their place and `C`
/// elsewhere.
fn chromosome() -> Vec<u8> {
    let mut bases = vec![b'C'; 10_400];
    for (start, coding) in [
        (UNFLAGGED_START, eight_codons()),
        (FLAGGED_START, eight_codons()),
    ] {
        let cdna = format!("{UTR5}{coding}");
        bases[start - 1..start - 1 + cdna.len()].copy_from_slice(cdna.as_bytes());
    }
    bases
}

/// A record that deletes the coding bases `first` to `last` of the model that
/// starts at `model_start`, written from the base before them.
fn deletion(name: &str, model_start: usize, first: usize, last: usize) -> String {
    let bases = chromosome();
    let from = genomic(model_start, first - 1) - 1;
    let reference = std::str::from_utf8(&bases[from..genomic(model_start, last)]).unwrap();
    format!(
        "{CHROMOSOME}\t{}\t{name}\t{reference}\t{}\t.\t.\t.\n",
        from + 1,
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
    fs::write(
        &shard,
        json!([
            model(1, UNFLAGGED_START, &eight_codons(), &[]),
            model(2, FLAGGED_START, &eight_codons(), &["cds_end_NF"]),
        ])
        .to_string(),
    )
    .unwrap();
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
            deletion("first_base_of_the_last_codon", UNFLAGGED_START, 22, 22),
            deletion("an_earlier_codon_into_the_last", UNFLAGGED_START, 19, 22),
            deletion(
                "first_base_of_the_last_codon_flagged",
                FLAGGED_START,
                22,
                22
            ),
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

/// The first record deletes coding base 22 of the model without a flag, the
/// first base of codon 8. The changed sequence holds seven whole codons, the
/// reference's first seven residues, and two bases of an eighth; nothing
/// follows them. Ensembl VEP writes `p.Asp8delextTer?`: the residue deleted,
/// and a stop lost where none stands. The frame shifts at residue 8, and the
/// row names it as the first residue affected, `p.Asp8fs`.
///
/// The second deletes coding bases 19 to 22, the whole of codon 7 and the
/// first base of codon 8: six whole codons and two bases are left. Ensembl
/// VEP writes `p.Pro7_Asp8delextTer?`; the row names the first of the two
/// residues, `p.Pro7fs`.
///
/// The third deletes coding base 22 of the model flagged `cds_end_NF`.
/// Ensembl VEP's lost-stop test does not hold under the flag
/// (`VariationEffect.pm:1373`), so it writes the deletion alone, and the row
/// keeps that form, `p.Asp8del`.
#[test]
fn a_frameshift_cut_short_at_an_unclosed_coding_end_is_named_at_its_first_residue() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows(dir.path()),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "first_base_of_the_last_codon\t1:128\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t28\t22\t8\tD/X\tGac/ac\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.22del;HGVSp=ENSP00000000001.1:p.Asp8fs",
            "an_earlier_codon_into_the_last\t1:125-128\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t25-28\t19-22\t7-8\tPD/X\tCCGGac/ac\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.19_22del;HGVSp=ENSP00000000001.1:p.Pro7fs",
            "first_base_of_the_last_codon_flagged\t1:10228\t-\tENSG00000000002\tENST00000000002\tTranscript\tframeshift_variant\t28\t22\t8\tD/X\tGac/ac\t-\tIMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000000002.1:c.22del;HGVSp=ENSP00000000002.1:p.Asp8del",
        ]
    );
}
