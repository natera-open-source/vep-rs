// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reference codons of an insertion written exactly where an edited
//! RefSeq model holds bases the genome lacks, inside its coding sequence,
//! through the `vep` binary on a made-up cache and reference: the two
//! frameshift rows of the documented difference
//! `refseq_edited_model_inserted_span_codons`.
//!
//! Each of two made-up models is the genome of its own chromosome plus the
//! bases an `_rna_edit` puts after its coding position 14, where its two
//! mapper pairs meet. The first model holds six such bases, coding positions
//! 15 to 20, and its record inserts one base at that point. The second holds
//! three, coding positions 15 to 17, and its record inserts two. Each row
//! prints the span of the added bases as its `CDS_position` and is a
//! `frameshift_variant`: against the model, the record's genome has the
//! added bases replaced by the inserted ones.

use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

/// The genomic position of each model's first base on its chromosome.
const MODEL_START: usize = 101;

/// The first coding position the genome lacks, in both models.
const FIRST_ADDED: usize = 15;

/// One made-up edited model and the record written on its added bases.
struct Model {
    chromosome: &'static str,
    gene: &'static str,
    stable_id: &'static str,
    protein: &'static str,
    utr5: &'static str,
    utr3: &'static str,
    /// The coding sequence as codons with their residues, the stop codon last.
    codons: &'static [(&'static str, char)],
    /// How many coding positions from [`FIRST_ADDED`] on the genome lacks.
    added: usize,
    /// The record's name and the bases it inserts at the point.
    record: (&'static str, &'static str),
}

/// Six added bases: the last base of codon 5, codon 6 and the first two
/// bases of codon 7. The record inserts one base.
const SIX_ADDED: Model = Model {
    chromosome: "1",
    gene: "4",
    stable_id: "NM_000004.1",
    protein: "NP_000004.1",
    utr5: "GTCACA",
    utr3: "GGTCAT",
    codons: &[
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
    ],
    added: 6,
    record: ("one_base_on_six", "T"),
};

/// Three added bases: the last base of codon 5 and the first two bases of
/// codon 6. The record inserts two bases.
const THREE_ADDED: Model = Model {
    chromosome: "2",
    gene: "5",
    stable_id: "NM_000005.1",
    protein: "NP_000005.1",
    utr5: "TGACCT",
    utr3: "GACTTC",
    codons: &[
        ("ATG", 'M'),
        ("TCA", 'S'),
        ("GGT", 'G'),
        ("AAC", 'N'),
        ("CTG", 'L'),
        ("AGC", 'S'),
        ("TTT", 'F'),
        ("GAC", 'D'),
        ("CCA", 'P'),
        ("ATT", 'I'),
        ("GAG", 'E'),
        ("TGA", '*'),
    ],
    added: 3,
    record: ("two_bases_on_three", "CA"),
};

impl Model {
    fn coding(&self) -> String {
        self.codons.iter().map(|(codon, _)| *codon).collect()
    }

    /// The model's cDNA as the cache holds it, with the added bases.
    fn cdna(&self) -> String {
        format!("{}{}{}", self.utr5, self.coding(), self.utr3)
    }

    /// The cDNA position of coding position `coding`.
    fn cdna_position(&self, coding: usize) -> usize {
        self.utr5.len() + coding
    }

    /// The last coding position the genome lacks.
    fn last_added(&self) -> usize {
        FIRST_ADDED + self.added - 1
    }

    /// The bases the genome lacks, in the model's order.
    fn added_bases(&self) -> String {
        self.coding()[FIRST_ADDED - 1..self.last_added()].to_string()
    }

    /// The model's bases that the genome holds, in order: the cDNA before the
    /// added bases and the cDNA after them.
    fn held(&self) -> (String, String) {
        let cdna = self.cdna();
        (
            cdna[..self.cdna_position(FIRST_ADDED) - 1].to_string(),
            cdna[self.cdna_position(self.last_added())..].to_string(),
        )
    }

    /// The model's cache record: an `_rna_edit` that inserts the added bases
    /// after the last cDNA base before them, its span counted in the unedited
    /// cDNA, and one mapper pair on each side of them.
    fn cache_record(&self) -> serde_json::Value {
        let (before, after) = self.held();
        let last_before = before.len();
        let first_after = self.cdna_position(self.last_added()) + 1;
        let end = MODEL_START + before.len() + after.len() - 1;
        let coding_start = self.cdna_position(1);
        let coding_end = self.cdna_position(self.coding().len());
        let peptide: String = self.codons[..self.codons.len() - 1]
            .iter()
            .map(|(_, residue)| *residue)
            .collect();
        let edit = format!("{} {last_before} {}", last_before + 1, self.added_bases());
        json!({
            "stable_id": self.stable_id,
            "gene_stable_id": self.gene,
            "chr": self.chromosome,
            "start": MODEL_START,
            "end": end,
            "strand": 1,
            "biotype": "protein_coding",
            "source": "BestRefSeq",
            "bam_edit_status": "ok",
            "protein": self.protein,
            "exons": [{"start": MODEL_START, "end": end, "rank": 1, "phase": -1, "end_phase": -1}],
            "cdna_coding_start": coding_start,
            "cdna_coding_end": coding_end,
            "attributes": [{"code": "_rna_edit", "value": edit}],
            "translation": {"stable_id": self.protein, "version": 1, "start": coding_start,
                "end": coding_end - self.added},
            "variation_effect_feature_cache": {
                "codon_table": 1,
                "five_prime_utr": self.utr5,
                "three_prime_utr": self.utr3,
                "spliced_seq": self.cdna(),
                "translateable_seq": self.coding(),
                "peptide": peptide,
                "mapper": {"start_phase": -1, "cdna_coding_start": coding_start,
                    "cdna_coding_end": coding_end,
                    "pairs": [
                        {"from_start": 1, "from_end": last_before, "ori": 1,
                            "to_start": MODEL_START, "to_end": MODEL_START + last_before - 1},
                        {"from_start": first_after, "from_end": self.cdna().len(), "ori": 1,
                            "to_start": MODEL_START + last_before, "to_end": end},
                    ]},
            },
        })
    }

    /// The made-up chromosome: the bases the genome holds of the model at
    /// their place, `C` elsewhere.
    fn chromosome_bases(&self) -> Vec<u8> {
        let (before, after) = self.held();
        let mut bases = vec![b'C'; 300];
        let from = MODEL_START - 1;
        let on_genome = format!("{before}{after}");
        bases[from..from + on_genome.len()].copy_from_slice(on_genome.as_bytes());
        bases
    }

    /// The record: after the base of the last coding position before the
    /// added bases, the inserted bases.
    fn record_line(&self) -> String {
        let coding = self.coding();
        let anchor = &coding[FIRST_ADDED - 2..FIRST_ADDED - 1];
        let position = MODEL_START + self.cdna_position(FIRST_ADDED - 1) - 1;
        let (name, inserted) = self.record;
        format!(
            "{}\t{position}\t{name}\t{anchor}\t{anchor}{inserted}\t.\t.\t.\n",
            self.chromosome
        )
    }
}

/// Writes the two models' cache, the reference of the two chromosomes and the
/// two records into `dir`, and returns the cache, the reference and the input.
fn stage(dir: &Path) -> (std::path::PathBuf, std::path::PathBuf, std::path::PathBuf) {
    let models = [&SIX_ADDED, &THREE_ADDED];
    let cache = dir.join("cache");
    let mut reference = Vec::new();
    let mut index = String::new();
    let mut records = String::new();
    for model in models {
        let shard = cache
            .join("transcripts")
            .join(model.chromosome)
            .join("1-1000000.json");
        fs::create_dir_all(shard.parent().unwrap()).unwrap();
        fs::write(&shard, json!([model.cache_record()]).to_string()).unwrap();

        let bases = model.chromosome_bases();
        reference.extend_from_slice(format!(">{}\n", model.chromosome).as_bytes());
        index.push_str(&format!(
            "{}\t{len}\t{offset}\t{len}\t{width}\n",
            model.chromosome,
            len = bases.len(),
            offset = reference.len(),
            width = bases.len() + 1
        ));
        reference.extend_from_slice(&bases);
        reference.push(b'\n');
        records.push_str(&model.record_line());
    }
    fs::write(
        cache.join("info.json"),
        json!({"assembly": "GRCh38"}).to_string(),
    )
    .unwrap();
    let fasta = dir.join("reference.fa");
    fs::write(&fasta, reference).unwrap();
    fs::write(dir.join("reference.fa.fai"), index).unwrap();
    let input = dir.join("records.vcf");
    fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"),
    )
    .unwrap();
    (cache, fasta, input)
}

/// The lines of the default output of `vep --refseq` with `flags` on the two
/// records, from its column line on.
fn rows(dir: &Path, name: &str, flags: &[&str]) -> Vec<String> {
    let (cache, fasta, input) = stage(dir);
    let out = dir.join(name);
    let mut command = common::isolated_vep(dir);
    command
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(&cache)
        .args(["--assembly", "GRCh38", "--offline"])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .arg("--refseq");
    for flag in flags {
        if *flag == "--fasta" {
            command.arg(flag).arg(&fasta);
        } else {
            command.arg(flag);
        }
    }
    let output = command.stdin(Stdio::null()).output().unwrap();
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

/// The two rows as vep-rs writes them.
///
/// On the first model the inserted base stands on coding positions 15 to 20,
/// residues 5 to 7, so the reference codons are the model's coding positions
/// 13 to 21 and the reference amino acids its residues 5 to 7. The alternate
/// codons are four bases: coding positions 13 and 14, the inserted base and
/// coding position 21, one whole codon and one base of the next.
///
/// On the second model the two inserted bases stand on coding positions 15 to
/// 17, residues 5 and 6, so the reference codons are coding positions 13 to
/// 18 and the reference amino acids residues 5 and 6. The alternate codons
/// are five bases: coding positions 13 and 14, the two inserted bases and
/// coding position 18.
///
/// Ensembl VEP 116.3 writes the record's reference allele, which is empty,
/// over the span of the added bases before it reads the reference codons of
/// a model that carries an `_rna_edit` (`TranscriptVariationAllele.pm:872-874`
/// and `:877`). Its reference codons are then coding positions 13, 14 and 21
/// to 27 on the first model and 13, 14 and 18 to 21 on the second, and its
/// reference amino acids their translation, which names a residue of neither
/// model at the row's last position.
const ROWS: [&str; 3] = [
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
    "one_base_on_six\t1:120-121\tT\t4\tNM_000004.1\tTranscript\tframeshift_variant\t21-26\t15-20\t5-7\tAQD/AX\tgcacaggat/gcTt\t-\tIMPACT=HIGH;STRAND=1;BAM_EDIT=OK",
    "two_bases_on_three\t2:120-121\tCA\t5\tNM_000005.1\tTranscript\tframeshift_variant\t21-23\t15-17\t5-6\tLS/LX\tctgagc/ctCAc\t-\tIMPACT=HIGH;STRAND=1;BAM_EDIT=OK",
];

/// The run of the entry's record: `--refseq --use_given_ref`, which a cache of
/// edited models takes without a reference.
#[test]
fn an_insertion_on_bases_an_edited_model_adds_reads_the_models_own_reference_codons() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(rows(dir.path(), "given.txt", &["--use_given_ref"]), ROWS);
}

/// The same two rows on the reference, the run a cache of edited models
/// takes without `--use_given_ref`: an insertion's reference allele is empty
/// either way.
#[test]
fn the_same_rows_are_written_on_the_reference() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(rows(dir.path(), "reference.txt", &["--fasta"]), ROWS);
}
