// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows of an insertion or a deletion whose length change is no multiple
//! of three and whose first base lies inside the codon of a stop that the
//! transcript's peptide holds and reads on past, through the `vep` binary on
//! a made-up cache and reference: the documented difference
//! `frameshift_inside_a_held_stop_codon`.
//!
//! The model is one exon, as a polymorphic pseudogene's cache record gives
//! it: six bases of 5' UTR, a coding sequence of 24 codons and nine bases of
//! 3' UTR, with no sequence edit. Codon 10 is a stop codon, and the peptide
//! has a stop at residue 10 and names the 13 residues after it. Ensembl VEP's
//! `frameshift` returns 0 where the reference peptide of the allele's codons
//! begins with a stop (`Utils/VariationEffect.pm:1554`), and its walk for the
//! first changed residue of a frameshift ends on a stop that both
//! translations hold (`TranscriptVariationAllele.pm:2354` to `:2357`). Both
//! hold for the stop codon no residue follows.

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
/// base.
const AFTER_THE_MODEL: &str = "GTAA";

/// The residue number of the stop codon that the peptide holds as a stop.
const HELD_STOP: usize = 10;

/// The model's coding sequence: a start codon, eight alanine codons, a stop
/// codon, ten alanine codons and the four codons of `LAST_CODONS`: 72 bases.
fn coding() -> String {
    format!(
        "ATG{}TGA{}{}",
        "GCT".repeat(8),
        "GCT".repeat(10),
        LAST_CODONS.concat()
    )
}

/// The model's peptide: 23 letters, the tenth a stop.
fn peptide() -> String {
    let peptide = format!("M{}*{}TLS", "A".repeat(8), "A".repeat(10));
    assert_eq!(peptide.as_bytes()[HELD_STOP - 1], b'*');
    peptide
}

fn cdna() -> String {
    format!("{UTR5}{}{}", coding(), UTR3.concat())
}

/// The genomic position of coding base `coding`.
fn genomic(coding: usize) -> usize {
    MODEL_START + UTR5.len() + coding - 1
}

/// The model's cache record: the peptide holds the stop of codon 10 and no
/// sequence edit reads it as a residue.
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

/// The lines of the default output of `vep --hgvs` on the six records, from
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
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}{}{}{}{}{}",
            insertion("written_before_the_stop_codon", 27, "T"),
            insertion("written_inside_the_stop_codon", 28, "T"),
            deletion("a_base_lost_inside_the_stop_codon", 30, 30),
            deletion("five_bases_lost_from_inside_the_stop_codon", 29, 33),
            insertion("four_bases_added_and_the_stop_kept", 29, "AATG"),
            deletion("a_base_lost_and_the_stop_kept", 29, 29),
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

/// The last base of codon 9 and the first base of the stop codon of residue
/// 10 are the same base. The first record adds a copy of it between the two
/// and the second between the stop codon's first and second bases: one and
/// the same sequence, with one HGVSc. Residue 10 becomes a leucine and the
/// new frame meets its first stop at residue 26, the third to fifth bases
/// of the 3' UTR, sixteen residues from the lost stop. The first record's
/// first base lies in codon 9, and both engines write the frameshift and
/// its extension. The second lies inside the stop codon, where Ensembl
/// VEP's rule gives `coding_sequence_variant` and a deletion-insertion.
///
/// The third deletes the last base of the stop codon: the stop becomes a
/// tryptophan and the new frame meets a stop at residue 22, twelve residues
/// on. Ensembl VEP's rule gives `coding_sequence_variant` and the stop
/// unchanged.
///
/// The fourth deletes five bases from the stop codon's second base: the
/// stop becomes a cysteine and the new frame meets the same stop of the 3'
/// UTR, at residue 24, fourteen residues on. The alternate codon string is
/// the first base of the reference's, so Ensembl VEP's rule gives
/// `inframe_deletion` (`Utils/VariationEffect.pm:1178`).
///
/// The fifth adds four bases after the stop codon's second base, which
/// leave a stop codon at residue 10. Residue 11, an alanine, is the first
/// that changes, to a methionine, and the new frame's first stop is residue
/// 27, the seventeenth residue counted from the changed one. Ensembl VEP's
/// rule gives `inframe_insertion,coding_sequence_variant`.
///
/// The sixth deletes the middle base of the stop codon, which leaves a stop
/// codon at residue 10 as well: residue 11 becomes a leucine and the new
/// frame's first stop is residue 22, the twelfth residue counted from the
/// changed one. Ensembl VEP's rule gives `coding_sequence_variant` and the
/// stop unchanged.
#[test]
fn a_frameshift_length_allele_inside_a_stop_the_peptide_holds_is_a_frameshift() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows(dir.path()),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "written_before_the_stop_codon\t1:133-134\tT\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t33-34\t27-28\t9-10\t-/X\t-/T\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.28dup;HGVSp=ENSP00000000001.1:p.Ter10LeuextTer16;HGVS_OFFSET=1",
            "written_inside_the_stop_codon\t1:134-135\tT\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t34-35\t28-29\t10\t*/LX\ttga/tTga\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.28dup;HGVSp=ENSP00000000001.1:p.Ter10LeuextTer16",
            "a_base_lost_inside_the_stop_codon\t1:136\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t36\t30\t10\t*/X\ttgA/tg\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.30del;HGVSp=ENSP00000000001.1:p.Ter10TrpextTer12",
            "five_bases_lost_from_inside_the_stop_codon\t1:135-139\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t35-39\t29-33\t10-11\t*A/X\ttGAGCT/t\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.30_34del;HGVSp=ENSP00000000001.1:p.Ter10CysextTer14;HGVS_OFFSET=1",
            "four_bases_added_and_the_stop_kept\t1:135-136\tAATG\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t35-36\t29-30\t10\t*/*MX\ttga/tgAATGa\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.30_31insATGA;HGVSp=ENSP00000000001.1:p.Ala11MetfsTer17;HGVS_OFFSET=1",
            "a_base_lost_and_the_stop_kept\t1:135\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t35\t29\t10\t*/X\ttGa/ta\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000001.1:c.29del;HGVSp=ENSP00000000001.1:p.Ala11LeufsTer12",
        ]
    );
}
