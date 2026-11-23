// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation of a frameshift whose alternate translation ends
//! before a residue differs from the reference protein, through the `vep`
//! binary on a made-up cache and reference: the documented difference
//! `hgvsp_frameshift_with_no_changed_residue`.
//!
//! Two models of one exon each end on their last coding base, with no 3'
//! UTR. The first holds eight whole codons and one base of a ninth and is
//! flagged `cds_end_NF`; the second holds the same eight codons and a stop
//! codon. Ensembl VEP 116.3 walks the reference and the alternate translation
//! from the row's residue for the first one that differs
//! (`TranscriptVariationAllele.pm:2348`); where the alternate translation
//! ends first it returns the pair it read last, two equal residues, at the
//! position after them (`:2351-2352`, `:2361` and `:2364`), and writes a
//! frameshift that names a residue by the one before it, changed to itself.

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

const CHROMOSOME: &str = "1";

/// The genomic position of the first base of the model no stop codon closes.
const OPEN_START: usize = 101;

/// The genomic position of the first base of the model a stop codon closes,
/// out of the other model's reach.
const CLOSED_START: usize = 10_201;

const UTR5: &str = "GTCACA";

/// Codons 1 to 8 of both models: a start codon, five alanine codons, a
/// proline codon and an aspartate codon, 24 bases.
fn eight_codons() -> String {
    format!("ATG{}CCGGAC", "GCT".repeat(5))
}

/// The protein of both models: eight residues.
const PEPTIDE: &str = "MAAAAAPD";

/// The coding sequence that no stop codon closes: the eight codons and one
/// base of a ninth, 25 bases.
fn open_coding() -> String {
    format!("{}T", eight_codons())
}

/// The coding sequence that its stop codon closes: the eight codons and a
/// stop codon, 27 bases.
fn closed_coding() -> String {
    format!("{}TAA", eight_codons())
}

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
    for (start, coding) in [(OPEN_START, open_coding()), (CLOSED_START, closed_coding())] {
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

/// A record that writes `inserted` between the coding bases `after` and
/// `after + 1` of the model that starts at `model_start`.
fn insertion(name: &str, model_start: usize, after: usize, inserted: &str) -> String {
    let bases = chromosome();
    let at = genomic(model_start, after);
    let anchor = bases[at - 1] as char;
    format!("{CHROMOSOME}\t{at}\t{name}\t{anchor}\t{anchor}{inserted}\t.\t.\t.\n")
}

/// The lines of the default output of `vep --hgvs` on the six records, from
/// its column line on.
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
            model(1, OPEN_START, &open_coding(), &["cds_end_NF"]),
            model(2, CLOSED_START, &closed_coding(), &[]),
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
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}{}{}{}{}{}",
            deletion("changed_residue_found", OPEN_START, 21, 21),
            deletion("next_residue_left_two_bases", OPEN_START, 21, 22),
            deletion("next_residue_left_no_base", OPEN_START, 21, 24),
            insertion("last_codon_kept_by_an_insertion", OPEN_START, 23, "T"),
            deletion("last_codon_kept_by_a_deletion", OPEN_START, 24, 24),
            deletion("stop_codon_broken", CLOSED_START, 24, 24),
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

/// The first record deletes coding base 21 of the open model, the third base
/// of codon 7, which equals the base after it. Codon 7 still codes proline
/// and codon 8 reads the last three coding bases, a threonine codon where the
/// reference holds an aspartate: the walk finds that residue, and both
/// engines write the frameshift at residue 8.
///
/// The second deletes coding bases 21 and 22. Codon 7 still codes proline and
/// two bases follow it, which spell a leucine codon whatever base comes
/// third; a protein is read over whole codons, so no residue differs.
/// Ensembl VEP writes `p.Pro8ProfsTer?`.
///
/// The third deletes coding bases 21 to 24. Codon 7 closes on the last coding
/// base and still codes proline, and no base of residue 8 is left. Ensembl
/// VEP writes `p.Pro8ProfsTer?`, naming the reference's aspartate by the
/// proline before it.
///
/// The fourth writes one base between coding bases 23 and 24, and the fifth
/// deletes coding base 24: codon 8 still codes aspartate on both, the
/// protein's eight residues are unchanged, and the reference holds no residue
/// 9. Ensembl VEP writes `p.Asp9AspfsTer?` on both.
///
/// The sixth deletes coding base 24 of the closed model. Codon 8 closes on
/// the stop codon's first base and still codes aspartate; the stop codon's
/// last two bases follow it and are no codon. Ensembl VEP writes
/// `p.Asp9AspfsTer?`, naming the stop by the residue before it.
///
/// No residue the transcript's bases give changes on the last five records,
/// and each row carries no HGVSp.
#[test]
fn a_frameshift_that_changes_no_residue_the_transcript_holds_has_no_protein_notation() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows(dir.path()),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "changed_residue_found\t1:127\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t27\t21\t7\tP/X\tccG/cc\t-\tIMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000000001.1:c.22del;HGVSp=ENSP00000000001.1:p.Asp8ThrfsTer?;HGVS_OFFSET=1",
            "next_residue_left_two_bases\t1:127-128\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t27-28\t21-22\t7-8\tPD/PX\tccGGac/ccac\t-\tIMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000000001.1:c.21_22del",
            "next_residue_left_no_base\t1:127-130\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t27-30\t21-24\t7-8\tPD/X\tccGGAC/cc\t-\tIMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000000001.1:c.21_24del",
            "last_codon_kept_by_an_insertion\t1:129-130\tT\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t29-30\t23-24\t8\tD/DX\tgac/gaTc\t-\tIMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000000001.1:c.23_24insT",
            "last_codon_kept_by_a_deletion\t1:130\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t30\t24\t8\tD/X\tgaC/ga\t-\tIMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF;HGVSc=ENST00000000001.1:c.24del",
            "stop_codon_broken\t1:10230\t-\tENSG00000000002\tENST00000000002\tTranscript\tframeshift_variant\t30\t24\t8\tD/X\tgaC/ga\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000002.1:c.24del",
        ]
    );
}
