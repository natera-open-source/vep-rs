// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The count of the extension that names a frameshift whose first changed
//! residue is a stop, on a transcript whose peptide holds stops and reads on
//! past them, through the `vep` binary on a made-up cache and reference: the
//! third shape of the documented difference
//! `hgvsp_frameshift_count_past_a_stop_read_as_a_residue`.
//!
//! The model is one exon, as a polymorphic pseudogene's cache record gives
//! it: six bases of 5' UTR, a coding sequence of 24 codons and nine bases of
//! 3' UTR, with no sequence edit. Codons 5 and 10 are stop codons, and the
//! peptide has a stop at residues 5 and 10 and names the residues after each.
//! Five records each insert or delete one or two coding bases right before a
//! stop codon, so that the stop is the first residue that changes, and both
//! engines write the change as an extension named at that stop
//! (`TranscriptVariationAllele.pm:2007` to `:2010`). Ensembl VEP 116.3 counts
//! to the first stop of the alternate sequence's translation, which is codon
//! 5 for a change after it (`:2492`, `:2503`, `:2508` and `:2519`), and
//! writes the unknown count on all four records after that codon.

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

const CHROMOSOME: &str = "1";

/// The genomic position of the model's first base.
const MODEL_START: usize = 101;

const UTR5: &str = "GTCACA";

/// The 3' UTR: its third to fifth bases are a stop codon.
const UTR3: [&str; 3] = ["CA", "TGA", "TTTT"];

/// Codons 21 to 24 of the coding sequence: a threonine, a leucine and a
/// serine codon and the stop codon that closes it.
const LAST_CODONS: [&str; 4] = ["ACC", "CTA", "AGC", "TAG"];

/// The four bases of the made-up chromosome right after the model's last
/// base. Read in the frame that two bases inserted before the closing stop
/// open, the last two bases of the 3' UTR and the first of these spell a
/// leucine codon and the other three a stop codon, which no transcript holds.
const AFTER_THE_MODEL: &str = "GTAA";

/// The residue numbers of the two stop codons that the peptide holds as
/// stops.
const HELD_STOPS: [usize; 2] = [5, 10];

/// The model's coding sequence: a start codon, three alanine codons, a stop
/// codon, four alanine codons, a second stop codon, ten alanine codons and
/// the four codons of `LAST_CODONS`: 72 bases.
fn coding() -> String {
    format!(
        "ATG{}TGA{}TAG{}{}",
        "GCT".repeat(3),
        "GCT".repeat(4),
        "GCT".repeat(10),
        LAST_CODONS.concat()
    )
}

/// The model's peptide: 23 letters, the fifth and the tenth a stop.
fn peptide() -> String {
    let peptide = format!("M{}*{}*{}TLS", "A".repeat(3), "A".repeat(4), "A".repeat(10));
    for held in HELD_STOPS {
        assert_eq!(peptide.as_bytes()[held - 1], b'*');
    }
    peptide
}

fn cdna() -> String {
    format!("{UTR5}{}{}", coding(), UTR3.concat())
}

/// The genomic position of coding base `coding`.
fn genomic(coding: usize) -> usize {
    MODEL_START + UTR5.len() + coding - 1
}

/// The model's cache record: the peptide holds the stops of codons 5 and 10
/// and no sequence edit reads either as a residue.
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
        "biotype": "polymorphic_pseudogene",
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
            "seq_edits": [],
            "mapper": {"start_phase": -1, "cdna_coding_start": coding_start,
                "cdna_coding_end": coding_end,
                "pairs": [{"from_start": 1, "from_end": cdna.len(), "ori": 1,
                    "to_start": MODEL_START, "to_end": end}]},
        },
    })
}

/// The made-up chromosome: the model's bases at their place, the bases of
/// `AFTER_THE_MODEL` right after them and `C` elsewhere.
fn chromosome() -> Vec<u8> {
    let on_the_chromosome = format!("{}{AFTER_THE_MODEL}", cdna());
    let mut bases = vec![b'C'; 300];
    let from = MODEL_START - 1;
    bases[from..from + on_the_chromosome.len()].copy_from_slice(on_the_chromosome.as_bytes());
    bases
}

/// A record that inserts `inserted` between the coding bases `after` and
/// `after + 1`, written from the base before the insertion.
fn insertion(name: &str, after: usize, inserted: &str) -> String {
    let coding = coding();
    let before = &coding[after - 1..after];
    format!(
        "{CHROMOSOME}\t{}\t{name}\t{before}\t{before}{inserted}\t.\t.\t.\n",
        genomic(after)
    )
}

/// A record that deletes the coding base `base`, written from the base
/// before it.
fn deletion(name: &str, base: usize) -> String {
    let coding = coding();
    let reference = &coding[base - 2..base];
    format!(
        "{CHROMOSOME}\t{}\t{name}\t{reference}\t{}\t.\t.\t.\n",
        genomic(base - 1),
        &reference[..1]
    )
}

/// The lines of the default output of `vep --hgvs` on the five records, from
/// its column line on.
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
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}{}{}{}{}",
            insertion("at_the_first_held_stop", 12, "C"),
            insertion("at_a_later_held_stop", 27, "C"),
            deletion("to_a_stop_inside_the_coding_sequence", 27),
            insertion("at_the_closing_stop", 69, "C"),
            insertion("no_stop_in_the_transcript", 69, "CC"),
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

/// The first record inserts one base between codon 4 and the stop codon of
/// residue 5, the first stop the peptide holds: that stop is the first
/// changed residue, a leucine, and the new frame meets its first stop at
/// residue 26, the third to fifth bases of the 3' UTR, 21 residues from the
/// lost stop, as both engines write it: no stop stands before the change.
///
/// The second inserts one base between codon 9 and the stop codon of residue
/// 10. That stop becomes a leucine and the new frame meets the same stop of
/// the 3' UTR, sixteen residues on. Ensembl VEP's alternate translation has
/// its first stop at residue 5, 5 - 9 is not positive, and it writes
/// `p.Ter10LeuextTer?`.
///
/// The third deletes coding base 27, the last base of codon 9. That codon
/// takes the first base of the stop codon and stays an alanine codon, the
/// stop of residue 10 becomes an arginine, and the new frame meets its first
/// stop at residue 22, inside the coding sequence: twelve residues on, where
/// Ensembl VEP writes `p.Ter10ArgextTer?`. The base deleted is the same as
/// the first base of the stop codon, so the notation of the coding sequence
/// names base 28.
///
/// The fourth inserts one base between codon 23 and the stop codon that
/// closes the coding sequence, behind both held stops. The closing stop
/// becomes a leucine and the new frame meets the stop of the 3' UTR two
/// residues on, where Ensembl VEP writes `p.Ter24LeuextTer?`. The inserted
/// base repeats the base before it.
///
/// The fifth inserts two bases at the same place. The new frame reads four
/// codons to the end of the 3' UTR and meets no stop in the transcript; the
/// stop codon that stands in that frame on the chromosome, right after the
/// transcript's last base, is no base of the transcript. The stop is sought
/// through the 3' UTR and no farther, and both engines write the unknown
/// count.
#[test]
fn a_frameshift_that_changes_a_stop_the_peptide_holds_counts_its_extension_from_that_stop() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows(dir.path()),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "at_the_first_held_stop\t1:118-119\tC\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t18-19\t12-13\t4-5\t-/X\t-/C\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.12_13insC;HGVSp=ENSP00000000001.1:p.Ter5LeuextTer21",
            "at_a_later_held_stop\t1:133-134\tC\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t33-34\t27-28\t9-10\t-/X\t-/C\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.27_28insC;HGVSp=ENSP00000000001.1:p.Ter10LeuextTer16",
            "to_a_stop_inside_the_coding_sequence\t1:133\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t33\t27\t9\tA/X\tgcT/gc\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.28del;HGVSp=ENSP00000000001.1:p.Ter10ArgextTer12;HGVS_OFFSET=1",
            "at_the_closing_stop\t1:175-176\tC\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t75-76\t69-70\t23-24\t-/X\t-/C\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.69dup;HGVSp=ENSP00000000001.1:p.Ter24LeuextTer2",
            "no_stop_in_the_transcript\t1:175-176\tCC\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t75-76\t69-70\t23-24\t-/X\t-/CC\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.69_70insCC;HGVSp=ENSP00000000001.1:p.Ter24ProextTer?",
        ]
    );
}
