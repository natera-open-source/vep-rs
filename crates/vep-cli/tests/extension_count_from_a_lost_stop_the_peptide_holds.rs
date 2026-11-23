// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The count of an extension in the protein notation of a substitution that
//! changes a stop, on a transcript whose peptide holds a stop and names the
//! residues after it, through the `vep` binary on a made-up cache and
//! reference: the eleventh shape of the documented difference
//! `hgvsp_stop_codon_protein_comparison`.
//!
//! Two models of one exon each, as a polymorphic pseudogene's cache record
//! gives them, with six bases of 5' UTR and no sequence edit. The first has a
//! coding sequence of 24 codons, the last a stop codon, a peptide that holds
//! stops at residues 10 and 15, and eleven bases of 3' UTR. The second has 23
//! codons that no stop codon closes, a peptide that holds a stop at residue
//! 10, and no 3' UTR. Ensembl VEP 116.3 subtracts the peptide's length from
//! the place of the first stop of the alternate sequence's translation
//! (`TranscriptVariationAllele.pm:2503` and `:2513`) and writes the unknown
//! count on all four records.

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

const CHROMOSOME: &str = "1";

const UTR5: &str = "GTCACA";

/// The 3' UTR of the first model: two alanine codons, a stop codon and two
/// more bases.
const UTR3: &str = "GCTGCTTGATT";

/// The genomic position of the first model's first base, and of the
/// second's, which lies too far from the first for a record of one to reach
/// the other.
const FIRST_MODEL: usize = 101;
const SECOND_MODEL: usize = 10_301;

/// The three bases of the made-up chromosome right after the second model's
/// last base: a stop codon in the model's frame, which no transcript holds.
const AFTER_THE_SECOND_MODEL: &str = "TAA";

/// A coding sequence of 23 codons and the peptide that is its translation: a
/// start codon, alanine codons, and a stop codon at each residue number of
/// `stops`, which the peptide holds as a stop.
fn coding(stops: &[usize]) -> (String, String) {
    let mut bases = String::new();
    let mut peptide = String::new();
    for residue in 1..=23 {
        let (codon, letter) = if stops.contains(&residue) {
            ("TGA", '*')
        } else if residue == 1 {
            ("ATG", 'M')
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
    let (bases, peptide) = coding(&[10, 15]);
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
    let (second, _) = coding(&[10]);
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

/// A record that puts `base` in the place of coding base `coding` of the
/// model that starts at `start`.
fn substitution(name: &str, start: usize, coding: usize, base: char) -> String {
    let position = start + UTR5.len() + coding - 1;
    let reference = chromosome()[position - 1] as char;
    format!("{CHROMOSOME}\t{position}\t{name}\t{reference}\t{base}\t.\t.\t.\n")
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
    let (second, second_peptide) = coding(&[10]);
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
            substitution("to_the_next_held_stop", FIRST_MODEL, 30, 'G'),
            substitution("to_the_closing_stop", FIRST_MODEL, 45, 'G'),
            substitution("closing_stop_to_the_utr", FIRST_MODEL, 70, 'C'),
            substitution("no_stop_in_the_transcript", SECOND_MODEL, 30, 'G'),
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

/// The first record changes the third base of codon 10 of the first model, a
/// stop the peptide holds, to a tryptophan codon: the frame's next stop is
/// the one the peptide holds at residue 15, five residues on. The first stop
/// of Ensembl VEP's alternate translation is that residue 15, 15 - 1 - 23 is
/// not positive, and it writes `p.Ter10TrpextTer?`.
///
/// The second changes the third base of codon 15: the frame reads on to the
/// stop codon that closes the coding sequence, residue 24, nine residues on.
/// Ensembl VEP's first stop is residue 10, before the change.
///
/// The third changes the first base of the closing stop codon to a glutamine
/// codon: the frame reads two alanine codons of the 3' UTR and its stop
/// codon, residue 27, three residues on. Ensembl VEP's first stop is again
/// residue 10.
///
/// The fourth changes the stop of residue 10 of the second model, whose
/// coding sequence no stop codon closes and which has no 3' UTR. The stop
/// codon that stands in frame on the chromosome right after the model's last
/// base is no base of the transcript: the stop is sought through the 3' UTR
/// and no farther, and both engines write the unknown count.
#[test]
fn a_lost_stop_on_a_peptide_that_holds_a_stop_counts_to_the_next_stop_of_the_frame() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows(dir.path()),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "to_the_next_held_stop\t1:136\tG\tENSG00000000001\tENST00000000001\tTranscript\tstop_lost\t36\t30\t10\t*/W\ttgA/tgG\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.30A>G;HGVSp=ENSP00000000001.1:p.Ter10TrpextTer5",
            "to_the_closing_stop\t1:151\tG\tENSG00000000001\tENST00000000001\tTranscript\tstop_lost\t51\t45\t15\t*/W\ttgA/tgG\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.45A>G;HGVSp=ENSP00000000001.1:p.Ter15TrpextTer9",
            "closing_stop_to_the_utr\t1:176\tC\tENSG00000000001\tENST00000000001\tTranscript\tstop_lost\t76\t70\t24\t*/Q\tTag/Cag\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.70T>C;HGVSp=ENSP00000000001.1:p.Ter24GlnextTer3",
            "no_stop_in_the_transcript\t1:10336\tG\tENSG00000000002\tENST00000000002\tTranscript\tstop_lost\t36\t30\t10\t*/W\ttgA/tgG\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000002.1:c.30A>G;HGVSp=ENSP00000000002.1:p.Ter10TrpextTer?",
        ]
    );
}
