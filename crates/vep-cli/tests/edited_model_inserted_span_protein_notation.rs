// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation of an insertion that the 3' shift of `--hgvs` moves
//! onto bases an edited RefSeq model holds and the genome lacks, inside its
//! coding sequence, through the `vep` binary on a made-up cache and
//! reference: the shape of the documented difference
//! `refseq_edited_model_inserted_span_hgvsp`.
//!
//! The model is the genome plus nine bases after its coding position 14. An
//! `_rna_edit` puts them at coding positions 15 to 23, and the model's two
//! mapper pairs meet there. The nine bases end in the three genomic bases
//! before them, coding positions 12 to 14, and the base after them is not
//! their first. The record inserts 18 bases three bases before the nine,
//! between coding positions 11 and 12: turned by three, they are the nine
//! added bases twice, so the shift carries the insertion three bases on and
//! stops, and the genome that carries the record is the model with coding
//! positions 15 to 23 repeated.
//!
//! Beside that row, one test holds the `code` row of the entry's record to
//! the functions that implement the entry.

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

const CHROMOSOME: &str = "1";

/// The genomic position of the model's first base.
const MODEL_START: usize = 101;

const UTR5: &str = "GTCACA";

const UTR3: &str = "GGTCAT";

/// The model's coding sequence as codons with their residues: twelve sense
/// codons and a stop codon, 39 bases. Coding positions 15 to 23 are the bases
/// the genome lacks: the last base of codon 5, codons 6 and 7 and the first
/// two bases of codon 8.
const CODING: [(&str, char); 13] = [
    ("ATG", 'M'),
    ("GAA", 'E'),
    ("CCT", 'P'),
    ("CTT", 'L'),
    ("GCA", 'A'),
    ("CAG", 'Q'),
    ("GAT", 'D'),
    ("GCC", 'A'),
    ("TGG", 'W'),
    ("CAT", 'H'),
    ("AAA", 'K'),
    ("GTT", 'V'),
    ("TAA", '*'),
];

/// The first and the last coding position the genome lacks.
const ADDED: (usize, usize) = (15, 23);

fn coding() -> String {
    CODING.iter().map(|(codon, _)| *codon).collect()
}

/// The model's cDNA as the cache holds it, with the added bases.
fn cdna() -> String {
    format!("{UTR5}{}{UTR3}", coding())
}

/// The cDNA position of coding position `coding`.
fn cdna_position(coding: usize) -> usize {
    UTR5.len() + coding
}

/// The bases the genome lacks, in the model's order.
fn added_bases() -> String {
    coding()[ADDED.0 - 1..ADDED.1].to_string()
}

/// The model's bases that the genome holds, in order: the cDNA before the
/// added bases and the cDNA after them.
fn held() -> (String, String) {
    let cdna = cdna();
    (
        cdna[..cdna_position(ADDED.0) - 1].to_string(),
        cdna[cdna_position(ADDED.1)..].to_string(),
    )
}

/// The model's cache record: an `_rna_edit` that inserts the added bases
/// after the last cDNA base before them, its span counted in the unedited
/// cDNA, and one mapper pair on each side of them.
fn model() -> serde_json::Value {
    let (before, after) = held();
    let last_before = before.len();
    let first_after = cdna_position(ADDED.1) + 1;
    let end = MODEL_START + before.len() + after.len() - 1;
    let coding_start = cdna_position(1);
    let coding_end = cdna_position(coding().len());
    let peptide: String = CODING[..CODING.len() - 1]
        .iter()
        .map(|(_, residue)| *residue)
        .collect();
    let edit = format!("{} {last_before} {}", last_before + 1, added_bases());
    json!({
        "stable_id": "NM_000004.1",
        "gene_stable_id": "4",
        "chr": CHROMOSOME,
        "start": MODEL_START,
        "end": end,
        "strand": 1,
        "biotype": "protein_coding",
        "source": "BestRefSeq",
        "bam_edit_status": "ok",
        "protein": "NP_000004.1",
        "exons": [{"start": MODEL_START, "end": end, "rank": 1, "phase": -1, "end_phase": -1}],
        "cdna_coding_start": coding_start,
        "cdna_coding_end": coding_end,
        "attributes": [{"code": "_rna_edit", "value": edit}],
        "translation": {"stable_id": "NP_000004.1", "version": 1, "start": coding_start,
            "end": coding_end - added_bases().len()},
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "five_prime_utr": UTR5,
            "three_prime_utr": UTR3,
            "spliced_seq": cdna(),
            "translateable_seq": coding(),
            "peptide": peptide,
            "mapper": {"start_phase": -1, "cdna_coding_start": coding_start,
                "cdna_coding_end": coding_end,
                "pairs": [
                    {"from_start": 1, "from_end": last_before, "ori": 1,
                        "to_start": MODEL_START, "to_end": MODEL_START + last_before - 1},
                    {"from_start": first_after, "from_end": cdna().len(), "ori": 1,
                        "to_start": MODEL_START + last_before, "to_end": end},
                ]},
        },
    })
}

/// The made-up chromosome: the bases the genome holds of the model at their
/// place, `C` elsewhere.
fn chromosome() -> Vec<u8> {
    let (before, after) = held();
    let mut bases = vec![b'C'; 300];
    let from = MODEL_START - 1;
    let on_genome = format!("{before}{after}");
    bases[from..from + on_genome.len()].copy_from_slice(on_genome.as_bytes());
    bases
}

/// The record: after the base of coding position 11, the three genomic bases
/// before the added bases and the first six added bases, twice.
fn record() -> String {
    let coding = coding();
    let added = added_bases();
    let anchor = &coding[10..11];
    let inserted = format!("{}{}", &coding[11..14], &added[..6]).repeat(2);
    let position = MODEL_START + cdna_position(11) - 1;
    format!("{CHROMOSOME}\t{position}\tinserted\t{anchor}\t{anchor}{inserted}\t.\t.\t.\n")
}

/// The lines of the default output of `vep --refseq --hgvs` on the record,
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

    let input = dir.join("record.vcf");
    fs::write(
        &input,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}",
            record()
        ),
    )
    .unwrap();
    let out = dir.join("record.txt");
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
        .args(["--refseq", "--hgvs"])
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

/// The 3' shift carries the insertion three bases on, onto the nine bases
/// the model holds and the genome lacks, and the row prints `HGVS_OFFSET=3`.
/// vep-rs reads the reference codons of the protein notation from the
/// model's coding sequence as the cache holds it, at the span the shift moves
/// the insertion to: the reference peptide is the model's residues 5 to 8,
/// the alternate peptide repeats residues 6 to 8 after them, and the `HGVSp`
/// is the `dup` of residues 6 to 8, the model's own residues. Ensembl VEP
/// 116.3 writes the record's reference allele, which is empty, over that span
/// before it reads the reference codons of a model that carries an
/// `_rna_edit` (`TranscriptVariationAllele.pm:872-874` and `:877`), so its
/// reference peptide is read with the nine bases cut out. The row's position
/// columns, its consequence term, `Codons` and `Amino_acids` are read at the
/// point the record writes, inside residue 4.
#[test]
fn an_insertion_shifted_onto_bases_an_edited_model_adds_names_the_models_own_residues() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows(dir.path()),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "inserted\t1:117-118\tTGCACAGGATGCACAGGA\t4\tNM_000004.1\tTranscript\tinframe_insertion\t17-18\t11-12\t4\tL/LAQDAQD\tctt/ctTGCACAGGATGCACAGGAt\t-\tIMPACT=MODERATE;STRAND=1;HGVSc=NM_000004.1:c.15_23dup;HGVSp=NP_000004.1:p.Gln6_Ala8dup;HGVS_OFFSET=3;BAM_EDIT=OK",
        ]
    );
}

/// The items the `code` row of the entry's record names, in the record's
/// order: the strings of the `code` list of the `[implementation]` table of
/// `scripts/divergences/reproductions/refseq_edited_model_inserted_span_hgvsp.toml`.
fn code_row() -> Vec<String> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../scripts/divergences/reproductions")
        .join("refseq_edited_model_inserted_span_hgvsp.toml");
    let record = fs::read_to_string(&path).unwrap();
    let list = record
        .split_once("\n[implementation]\n")
        .and_then(|(_, table)| table.split_once("code = ["))
        .and_then(|(_, rest)| rest.split_once(']'))
        .unwrap_or_else(|| panic!("{}: no code row", path.display()))
        .0;
    list.split('"')
        .skip(1)
        .step_by(2)
        .map(str::to_string)
        .collect()
}

/// The record names the two functions that keep the reference peptide of the
/// protein notation on the model's own residues. `PerlCodingEval::codon`
/// reads the reference codons from the model's coding sequence as the cache
/// holds it, without the write of `TranscriptVariationAllele.pm:872-874`.
/// Where the model's cache record gives no `seq_edits`, as on the model of
/// this file, `overlay_seq_edits` then writes the cached peptide's residue
/// over each residue of the translated read that differs from it, except at
/// position 1 and on a partial terminal codon. Each of the two keeps the row
/// of the test above without the other: a read with the added bases cut out
/// takes the model's residues back from the cached peptide, and a read of the
/// model's own bases needs no residue from it. The row is Ensembl VEP's
/// `delins` only with the cut in `codon` and no overlay.
#[test]
fn the_record_names_both_functions_that_keep_the_reference_peptide_on_the_models_residues() {
    assert_eq!(
        code_row(),
        [
            "crates/vep-effects/src/coding.rs::PerlCodingEval::codon",
            "crates/vep-effects/src/coding.rs::overlay_seq_edits",
        ]
    );
}
