// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The position inside the HGVSc of a row past the own edit of a RefSeq
//! model whose alignment correction failed (`BAM_EDIT=FAILED`), through the
//! `vep` binary on a made-up cache and reference: the shape of the documented
//! difference `refseq_failed_edit_model_hgvsc_offset`.
//!
//! The model's coding sequence reads one genomic base twice. The cache
//! writes the double reading as one extra position, coding position 8, by an
//! `_rna_edit` that writes the one cDNA position 14 as two bases, and its two
//! mapper pairs meet there: past the edit a cached cDNA number is the count
//! along the exon plus one. The extra position is no base of the transcript,
//! so vep-rs numbers every notation along the exon. Ensembl VEP 116.3
//! numbers a one-base coding substitution at the row's CDS start, which
//! counts the extra position (`TranscriptVariationAllele.pm:1519-1522`), and
//! any other notation along the exon plus the two bases the edit writes
//! (`:1524` and `:1617`).

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

const CHROMOSOME: &str = "1";

/// The genomic position of the model's first base.
const MODEL_START: usize = 101;

const UTR5: &str = "GTCACA";

const UTR3: &str = "GGTCATCGATGCTAGC";

/// The model's coding sequence as the cache holds it, as codons with their
/// residues: eight sense codons and a stop codon, 27 bases. Coding position
/// 8, the middle base of codon 3, is the position the genome lacks.
const CODING: [(&str, char); 9] = [
    ("ATG", 'M'),
    ("GCA", 'A'),
    ("AAA", 'K'),
    ("CCT", 'P'),
    ("GAT", 'D'),
    ("TGG", 'W'),
    ("CAT", 'H'),
    ("GTT", 'V'),
    ("TAA", '*'),
];

/// The coding position the cache holds and the genome lacks.
const EXTRA: usize = 8;

fn coding() -> String {
    CODING.iter().map(|(codon, _)| *codon).collect()
}

/// The cached cDNA position of coding position `coding`.
fn cdna_position(coding: usize) -> usize {
    UTR5.len() + coding
}

/// The model's bases as the genome holds them: the cached cDNA without the
/// extra position.
fn on_genome() -> String {
    let cached = format!("{UTR5}{}{UTR3}", coding());
    let extra = cdna_position(EXTRA);
    format!("{}{}", &cached[..extra - 1], &cached[extra..])
}

/// The genomic position of a cached cDNA position past the extra one.
fn genomic_past_the_edit(cdna: usize) -> usize {
    MODEL_START + cdna - 2
}

/// The model's cache record: status `failed`, no cached cDNA, an `_rna_edit`
/// that writes the cDNA position before which the extra one stands as two
/// bases, its span counted along the exon, and one mapper pair on each side
/// of the extra position.
fn model() -> serde_json::Value {
    let extra = cdna_position(EXTRA);
    let cached_length = UTR5.len() + coding().len() + UTR3.len();
    let end = MODEL_START + on_genome().len() - 1;
    let coding_start = cdna_position(1);
    let coding_end = cdna_position(coding().len());
    let peptide: String = CODING[..CODING.len() - 1]
        .iter()
        .map(|(_, residue)| *residue)
        .collect();
    let read_twice = &coding()[EXTRA..=EXTRA];
    let edit = format!("{extra} {extra} {read_twice}{read_twice}");
    json!({
        "stable_id": "NM_000007.1",
        "gene_stable_id": "7",
        "chr": CHROMOSOME,
        "start": MODEL_START,
        "end": end,
        "strand": 1,
        "biotype": "protein_coding",
        "source": "BestRefSeq",
        "bam_edit_status": "failed",
        "protein": "NP_000007.1",
        "exons": [{"start": MODEL_START, "end": end, "rank": 1, "phase": -1, "end_phase": -1}],
        "cdna_coding_start": coding_start,
        "cdna_coding_end": coding_end,
        "attributes": [{"code": "_rna_edit", "value": edit}],
        "translation": {"stable_id": "NP_000007.1", "version": 1, "start": coding_start,
            "end": coding_end - 1},
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "five_prime_utr": UTR5,
            "three_prime_utr": UTR3,
            "translateable_seq": coding(),
            "peptide": peptide,
            "mapper": {"start_phase": -1, "cdna_coding_start": coding_start,
                "cdna_coding_end": coding_end,
                "pairs": [
                    {"from_start": 1, "from_end": extra - 1, "ori": 1,
                        "to_start": MODEL_START, "to_end": MODEL_START + extra - 2},
                    {"from_start": extra + 1, "from_end": cached_length, "ori": 1,
                        "to_start": MODEL_START + extra - 1, "to_end": end},
                ]},
        },
    })
}

/// The made-up chromosome: the model's bases at their place, `C` elsewhere.
fn chromosome() -> Vec<u8> {
    let mut bases = vec![b'C'; 300];
    let from = MODEL_START - 1;
    let on_genome = on_genome();
    bases[from..from + on_genome.len()].copy_from_slice(on_genome.as_bytes());
    bases
}

/// Six records, each written where the bases beside it differ, so that no 3'
/// shift moves it: a deletion of cached coding position 5, before the edit; a
/// deletion of cached coding position 14, a substitution of position 17 and
/// an insertion between positions 19 and 20, past it; and in the 3' UTR a
/// substitution of its fifth base and an insertion between its seventh and
/// eighth.
fn records() -> String {
    let coding = coding();
    let base = |coding_position: usize| &coding[coding_position - 1..coding_position];
    let utr = |position: usize| &UTR3[position - 1..position];
    let coding_end = cdna_position(coding.len());
    let lines = [
        format!(
            "{}\tcoding_deletion_before_the_edit\t{}{}\t{}",
            MODEL_START + cdna_position(4) - 1,
            base(4),
            base(5),
            base(4)
        ),
        format!(
            "{}\tcoding_deletion_past_the_edit\t{}{}\t{}",
            genomic_past_the_edit(cdna_position(13)),
            base(13),
            base(14),
            base(13)
        ),
        format!(
            "{}\tcoding_substitution_past_the_edit\t{}\tC",
            genomic_past_the_edit(cdna_position(17)),
            base(17)
        ),
        format!(
            "{}\tcoding_insertion_past_the_edit\t{}\t{}G",
            genomic_past_the_edit(cdna_position(19)),
            base(19),
            base(19)
        ),
        format!(
            "{}\tutr_substitution\t{}\tG",
            genomic_past_the_edit(coding_end + 5),
            utr(5)
        ),
        format!(
            "{}\tutr_insertion\t{}\t{}A",
            genomic_past_the_edit(coding_end + 7),
            utr(7),
            utr(7)
        ),
    ];
    lines
        .iter()
        .map(|line| format!("{CHROMOSOME}\t{line}\t.\t.\t.\n"))
        .collect()
}

/// The rows of the default output of `vep --refseq --hgvs` on the records:
/// each row's record, `cDNA_position`, `CDS_position` and the `Extra` items
/// named in `items`, joined by ` | `.
fn rows(dir: &Path, items: &[&str]) -> Vec<String> {
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
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}",
            records()
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
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            let extra = |name: &str| {
                cells[13]
                    .split(';')
                    .find_map(|item| item.strip_prefix(&format!("{name}=")))
                    .unwrap_or("-")
            };
            let mut row = vec![cells[0], cells[6], cells[7], cells[8]];
            row.extend(items.iter().map(|name| extra(name)));
            row.join(" | ")
        })
        .collect()
}

/// Each notation stands at the count along the exon, the transcript's own
/// bases. A coding one past the edit is one below the row's `CDS_position`,
/// which counts the cache's extra position, and a 3' UTR one is the count of
/// bases from the stop codon's last base, which the edit before the stop
/// codon does not change. `cDNA_position` and `CDS_position` stay the cached
/// numbers. Before the edit the notation is the row's `CDS_position`. Every
/// row prints the model's status, and a row past the edit the offset Ensembl
/// VEP sums for it, 2.
#[test]
fn a_notation_past_the_edit_of_a_failed_edit_model_is_counted_along_the_exon() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        rows(dir.path(), &["BAM_EDIT", "REFSEQ_OFFSET", "HGVSc"]),
        [
            "coding_deletion_before_the_edit | frameshift_variant | 11 | 5 | FAILED | - | NM_000007.1:c.5del",
            "coding_deletion_past_the_edit | frameshift_variant | 20 | 14 | FAILED | 2 | NM_000007.1:c.13del",
            "coding_substitution_past_the_edit | missense_variant | 23 | 17 | FAILED | 2 | NM_000007.1:c.16G>C",
            "coding_insertion_past_the_edit | frameshift_variant | 25-26 | 19-20 | FAILED | 2 | NM_000007.1:c.18_19insG",
            "utr_substitution | 3_prime_UTR_variant | 38 | - | FAILED | 2 | NM_000007.1:c.*5A>G",
            "utr_insertion | 3_prime_UTR_variant | 40-41 | - | FAILED | 2 | NM_000007.1:c.*7_*8insA",
        ]
    );
}
