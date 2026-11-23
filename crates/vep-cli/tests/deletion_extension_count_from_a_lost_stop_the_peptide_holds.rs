// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The count of an extension in the protein notation of a deletion of whole
//! residues that takes a stop, on a transcript whose peptide holds a stop and
//! names the residues after it, through the `vep` binary on a made-up cache
//! and reference: the eleventh shape of the documented difference
//! `hgvsp_stop_codon_protein_comparison`.
//!
//! Two models of one exon each, as a polymorphic pseudogene's cache record
//! gives them, with six bases of 5' UTR and no sequence edit. The first has a
//! coding sequence of 24 codons, the last a stop codon, a peptide that holds
//! stops at residues 10 and 15, an aspartate at residue 9 and a proline at
//! residue 11, and eleven bases of 3' UTR. The second has 23 codons that no
//! stop codon closes, a peptide that holds a stop at residue 10, and no 3'
//! UTR. Ensembl VEP 116.3 subtracts one and the peptide's length from the
//! place of the first stop of the alternate sequence's translation
//! (`TranscriptVariationAllele.pm:2503` and `:2513`) and writes the unknown
//! count on all four records.

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

const CHROMOSOME: &str = "1";

const UTR5: &str = "GTCACA";

/// The 3' UTR of the first model: a proline codon, an alanine codon, a stop
/// codon and two more bases.
const UTR3: &str = "CCTGCTTGATT";

/// The genomic position of the first model's first base, and of the
/// second's, which lies too far from the first for a record of one to reach
/// the other.
const FIRST_MODEL: usize = 101;
const SECOND_MODEL: usize = 10_301;

/// The three bases of the made-up chromosome right after the second model's
/// last base: a stop codon in the model's frame, which no transcript holds.
const AFTER_THE_SECOND_MODEL: &str = "TAA";

/// A coding sequence of 23 codons and the peptide that is its translation: a
/// start codon, alanine codons, a stop codon at each residue number of
/// `stops`, which the peptide holds as a stop, and the codon and residue of
/// `others` at their residue numbers.
fn coding(stops: &[usize], others: &[(usize, &str, char)]) -> (String, String) {
    let mut bases = String::new();
    let mut peptide = String::new();
    for residue in 1..=23 {
        let other = others.iter().find(|(at, _, _)| *at == residue);
        let (codon, letter) = if stops.contains(&residue) {
            ("TGA", '*')
        } else if residue == 1 {
            ("ATG", 'M')
        } else if let Some((_, codon, letter)) = other {
            (*codon, *letter)
        } else {
            ("GCT", 'A')
        };
        bases.push_str(codon);
        peptide.push(letter);
    }
    (bases, peptide)
}

/// The first model's coding sequence, closed by a stop codon, and its
/// peptide.
fn closed() -> (String, String) {
    let (bases, peptide) = coding(&[10, 15], &[(9, "GAT", 'D'), (11, "CCT", 'P')]);
    (format!("{bases}TAG"), peptide)
}

/// A model's cache record: `number` names it, it starts at `start`, and its
/// peptide holds its stops with no sequence edit to read one as a residue.
fn model(
    number: usize,
    start: usize,
    coding: &str,
    peptide: &str,
    utr3: &str,
) -> serde_json::Value {
    let cdna = format!("{UTR5}{coding}{utr3}");
    let end = start + cdna.len() - 1;
    let coding_start = UTR5.len() + 1;
    let coding_end = UTR5.len() + coding.len();
    json!({
        "stable_id": format!("ENST0000000000{number}"),
        "version": 1,
        "gene_stable_id": format!("ENSG0000000000{number}"),
        "chr": CHROMOSOME,
        "start": start,
        "end": end,
        "strand": 1,
        "biotype": "polymorphic_pseudogene",
        "source": "ensembl",
        "protein": format!("ENSP0000000000{number}.1"),
        "exons": [{"start": start, "end": end, "rank": 1, "phase": -1, "end_phase": -1}],
        "cdna_coding_start": coding_start,
        "cdna_coding_end": coding_end,
        "attributes": [],
        "translation": {"stable_id": format!("ENSP0000000000{number}"), "version": 1,
            "start": coding_start, "end": coding_end},
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "five_prime_utr": UTR5,
            "three_prime_utr": utr3,
            "spliced_seq": cdna,
            "translateable_seq": coding,
            "peptide": peptide,
            "seq_edits": [],
            "mapper": {"start_phase": -1, "cdna_coding_start": coding_start,
                "cdna_coding_end": coding_end,
                "pairs": [{"from_start": 1, "from_end": cdna.len(), "ori": 1,
                    "to_start": start, "to_end": end}]},
        },
    })
}

/// The made-up chromosome: each model's bases at their place, the bases of
/// `AFTER_THE_SECOND_MODEL` right after the second, and `C` elsewhere.
fn chromosome() -> Vec<u8> {
    let (first, _) = closed();
    let (second, _) = coding(&[10], &[]);
    let mut bases = vec![b'C'; 10_500];
    for (start, laid) in [
        (FIRST_MODEL, format!("{UTR5}{first}{UTR3}")),
        (
            SECOND_MODEL,
            format!("{UTR5}{second}{AFTER_THE_SECOND_MODEL}"),
        ),
    ] {
        bases[start - 1..start - 1 + laid.len()].copy_from_slice(laid.as_bytes());
    }
    bases
}

/// A record that deletes the coding bases `from` to `to` of the model that
/// starts at `start`, written with the base before them.
fn deletion(name: &str, start: usize, from: usize, to: usize) -> String {
    let bases = chromosome();
    let before = start + UTR5.len() + from - 2;
    let last = start + UTR5.len() + to - 1;
    let reference = String::from_utf8_lossy(&bases[before - 1..last]).into_owned();
    let kept = &reference[..1];
    format!("{CHROMOSOME}\t{before}\t{name}\t{reference}\t{kept}\t.\t.\t.\n")
}

/// The lines of the default output of `vep --hgvs` on the four records, from
/// its column line on.
fn rows(dir: &Path) -> Vec<String> {
    let cache = dir.join("cache");
    let shard = cache
        .join("transcripts")
        .join(CHROMOSOME)
        .join("1-1000000.json");
    fs::create_dir_all(shard.parent().unwrap()).unwrap();
    let (first, first_peptide) = closed();
    let (second, second_peptide) = coding(&[10], &[]);
    fs::write(
        &shard,
        json!([
            model(1, FIRST_MODEL, &first, &first_peptide, UTR3),
            model(2, SECOND_MODEL, &second, &second_peptide, ""),
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
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}{}{}{}",
            deletion("through_a_held_stop", FIRST_MODEL, 25, 30),
            deletion("opens_at_a_held_stop", FIRST_MODEL, 43, 48),
            deletion("closing_stop_behind_held_stops", FIRST_MODEL, 67, 72),
            deletion("no_stop_in_the_transcript", SECOND_MODEL, 28, 33),
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

/// The first record deletes the codons of residues 9 and 10 of the first
/// model, an aspartate and a stop the peptide holds. The proline of residue
/// 11 stands at residue 9, the first changed residue, and the stop the
/// peptide holds at residue 15 stands at residue 13: three places after the
/// lost stop's. The first stop of Ensembl VEP's alternate translation is that
/// stop, 13 - 1 - 23 is not above nought, and it writes
/// `p.Asp9_Ter10delextTer?`.
///
/// The second deletes the codons of residues 15 and 16, so the description
/// opens at the stop: the stop codon that closes the coding sequence stands
/// at residue 22, seven places on. Ensembl VEP's first stop is residue 10,
/// before the change.
///
/// The third deletes the codons of residue 23 and of the closing stop,
/// residue 24: the proline codon of the 3' UTR stands at residue 23 and the
/// stop codon of the 3' UTR at residue 25, one place on. Ensembl VEP's first
/// stop is again residue 10.
///
/// The fourth deletes the codons of residues 10 and 11 of the second model,
/// whose coding sequence no stop codon closes and which has no 3' UTR. The
/// stop codon that stands in frame on the chromosome right after the model's
/// last base is no base of the transcript: the stop is sought through the 3'
/// UTR and no farther, and both engines write the unknown count.
#[test]
fn a_deletion_through_a_stop_the_peptide_holds_counts_to_the_first_stop_after_the_change() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows(dir.path()),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "through_a_held_stop\t1:131-136\t-\tENSG00000000001\tENST00000000001\tTranscript\tstop_lost,inframe_deletion\t31-36\t25-30\t9-10\tD*/-\tGATTGA/-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.25_30del;HGVSp=ENSP00000000001.1:p.Asp9_Ter10delextTer3",
            "opens_at_a_held_stop\t1:149-154\t-\tENSG00000000001\tENST00000000001\tTranscript\tstop_lost,inframe_deletion\t49-54\t43-48\t15-16\t*A/-\tTGAGCT/-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.43_48del;HGVSp=ENSP00000000001.1:p.Ter15_Ala16delextTer7",
            "closing_stop_behind_held_stops\t1:173-178\t-\tENSG00000000001\tENST00000000001\tTranscript\tstop_lost,inframe_deletion\t73-78\t67-72\t23-24\tA*/-\tGCTTAG/-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.67_72del;HGVSp=ENSP00000000001.1:p.Ala23_Ter24delextTer1",
            "no_stop_in_the_transcript\t1:10334-10339\t-\tENSG00000000002\tENST00000000002\tTranscript\tstop_lost,inframe_deletion\t34-39\t28-33\t10-11\t*A/-\tTGAGCT/-\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000002.1:c.28_33del;HGVSp=ENSP00000000002.1:p.Ter10_Ala11delextTer?",
        ]
    );
}
