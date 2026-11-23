// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The count to the new stop codon in the protein notation of a frameshift
//! whose first changed residue is residue 1, through the `vep` binary: the
//! documented difference `hgvsp_frameshift_count_at_the_first_residue`.
//!
//! Ensembl VEP 116.3 makes no count for such a frameshift on any transcript:
//! the frameshift form hands its count routine the first changed residue's
//! number less one (`TranscriptVariationAllele.pm:2003`), the routine returns
//! nothing for the number 0 (`:2483`), and the form writes the unknown count
//! (`:2006`). vep-rs counts where no start codon stands at the coding start,
//! on a transcript flagged `cds_start_NF` or one whose record shows no start
//! codon, and keeps the unknown count where a start codon stands.
//!
//! The first test runs three made-up models, one on each of three made-up
//! chromosomes, that share a coding body of fifteen codons and a stop codon.
//! Read from its third base the body meets a stop codon at the fourth codon,
//! and read from the second base behind a start codon at the seventh. The
//! second test runs the first record of the entry's reproduction on the
//! committed GRCh37 cache.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

use serde_json::json;

mod common;

/// The genomic position of each model's first base on its chromosome.
const MODEL_START: usize = 101;

/// The first seven codons of the coding body the models share: a glutamate,
/// a proline, a valine, an alanine, a serine, a leucine and a threonine
/// codon. Read from the body's third base they spell a threonine and two
/// cysteine codons and a stop codon.
const BODY: [&str; 7] = ["GAA", "CCT", "GTT", "GCT", "AGC", "CTA", "ACT"];

/// The 5' UTR of the model whose start codon stands: its last base is the
/// first base of the start codon.
const UTR5: [&str; 2] = ["CCT", "TCA"];

const UTR3: [&str; 2] = ["CAT", "CC"];

/// The coding body: the codons of `BODY`, eight alanine codons and a stop
/// codon, 48 bases.
fn body() -> String {
    format!("{}{}TAA", BODY.concat(), "GCT".repeat(8))
}

/// The body's residues, without the stop.
fn body_peptide() -> String {
    format!("EPVASLT{}", "A".repeat(8))
}

/// One made-up model on a chromosome of its own.
struct Model {
    chromosome: &'static str,
    number: usize,
    utr5: String,
    coding: String,
    peptide: String,
    flagged: bool,
}

impl Model {
    /// The body alone, flagged `cds_start_NF`: no 5' UTR and no start codon.
    fn flagged() -> Self {
        Self {
            chromosome: "1",
            number: 1,
            utr5: String::new(),
            coding: body(),
            peptide: body_peptide(),
            flagged: true,
        }
    }

    /// The body alone without the flag: its record shows no start codon.
    fn without_a_start_codon() -> Self {
        Self {
            chromosome: "2",
            number: 2,
            flagged: false,
            ..Self::flagged()
        }
    }

    /// The body behind a 5' UTR and a start codon.
    fn with_a_start_codon() -> Self {
        Self {
            chromosome: "3",
            number: 3,
            utr5: UTR5.concat(),
            coding: format!("ATG{}", body()),
            peptide: format!("M{}", body_peptide()),
            flagged: false,
        }
    }

    fn cdna(&self) -> String {
        format!("{}{}{}", self.utr5, self.coding, UTR3.concat())
    }

    /// The genomic position of coding base `coding`.
    fn genomic(&self, coding: usize) -> usize {
        MODEL_START + self.utr5.len() + coding - 1
    }

    /// The chromosome: the model's bases at their place, `C` elsewhere.
    fn bases(&self) -> Vec<u8> {
        let cdna = self.cdna();
        let mut bases = vec![b'C'; 300];
        let from = MODEL_START - 1;
        bases[from..from + cdna.len()].copy_from_slice(cdna.as_bytes());
        bases
    }

    /// The model's cache record. A model without a 5' UTR opens in phase 0.
    fn record(&self) -> serde_json::Value {
        let cdna = self.cdna();
        let end = MODEL_START + cdna.len() - 1;
        let coding_start = self.utr5.len() + 1;
        let coding_end = self.utr5.len() + self.coding.len();
        let phase = if self.utr5.is_empty() { 0 } else { -1 };
        let attributes = if self.flagged {
            json!([{"code": "cds_start_NF", "name": "CDS start not found", "value": "1"}])
        } else {
            json!([])
        };
        json!({
            "stable_id": format!("ENST0000000000{}", self.number),
            "version": 1,
            "gene_stable_id": format!("ENSG0000000000{}", self.number),
            "chr": self.chromosome,
            "start": MODEL_START,
            "end": end,
            "strand": 1,
            "biotype": "protein_coding",
            "source": "ensembl",
            "protein": format!("ENSP0000000000{}.1", self.number),
            "exons": [{"start": MODEL_START, "end": end, "rank": 1, "phase": phase,
                "end_phase": -1}],
            "cdna_coding_start": coding_start,
            "cdna_coding_end": coding_end,
            "attributes": attributes,
            "translation": {"stable_id": format!("ENSP0000000000{}", self.number),
                "version": 1, "start": coding_start, "end": coding_end},
            "variation_effect_feature_cache": {
                "codon_table": 1,
                "five_prime_utr": self.utr5,
                "three_prime_utr": UTR3.concat(),
                "spliced_seq": cdna,
                "translateable_seq": self.coding,
                "peptide": self.peptide,
                "seq_edits": [],
                "mapper": {"start_phase": phase, "cdna_coding_start": coding_start,
                    "cdna_coding_end": coding_end,
                    "pairs": [{"from_start": 1, "from_end": cdna.len(), "ori": 1,
                        "to_start": MODEL_START, "to_end": end}]},
            },
        })
    }

    /// A record that deletes the coding bases `first` to `last`, written from
    /// the base before them.
    fn deletion(&self, name: &str, first: usize, last: usize) -> String {
        let bases = self.bases();
        let before = self.genomic(first) - 1;
        let reference = String::from_utf8(bases[before - 1..self.genomic(last)].to_vec()).unwrap();
        format!(
            "{}\t{before}\t{name}\t{reference}\t{}\t.\t.\t.\n",
            self.chromosome,
            &reference[..1]
        )
    }

    /// A record that inserts a copy of coding base 1 before it, written on
    /// the last base of the 5' UTR.
    fn insertion_before_the_coding_start(&self, name: &str) -> String {
        let bases = self.bases();
        let before = self.genomic(1) - 1;
        let reference = bases[before - 1] as char;
        format!(
            "{}\t{before}\t{name}\t{reference}\t{reference}{}\t.\t.\t.\n",
            self.chromosome,
            &self.coding[..1]
        )
    }
}

/// The lines of the default output of `vep --hgvs` on `records` against
/// `cache` and `fasta`, from its column line on.
fn rows(dir: &Path, cache: &Path, fasta: &Path, assembly: &str, records: &str) -> Vec<String> {
    let input = dir.join("records.vcf");
    fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"),
    )
    .unwrap();
    let out = dir.join("records.txt");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache)
        .arg("--fasta")
        .arg(fasta)
        .args(["--assembly", assembly, "--offline"])
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

/// The made-up cache and reference of `models` under `dir`.
fn made_up(dir: &Path, models: &[Model]) -> (PathBuf, PathBuf) {
    let cache = dir.join("cache");
    let mut reference = Vec::new();
    let mut index = String::new();
    for model in models {
        let shard = cache
            .join("transcripts")
            .join(model.chromosome)
            .join("1-1000000.json");
        fs::create_dir_all(shard.parent().unwrap()).unwrap();
        fs::write(&shard, json!([model.record()]).to_string()).unwrap();
        let bases = model.bases();
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
    }
    fs::write(
        cache.join("info.json"),
        json!({"assembly": "GRCh38"}).to_string(),
    )
    .unwrap();
    let fasta = dir.join("reference.fa");
    fs::write(&fasta, reference).unwrap();
    fs::write(dir.join("reference.fa.fai"), index).unwrap();
    (cache, fasta)
}

/// On the flagged model, a deletion of coding bases 1 and 2 changes residue
/// 1, a glutamate, to a threonine, and the new frame reads two cysteine
/// codons and a stop codon, its fourth codon: a count of 4. A deletion of
/// coding bases 4 and 5 keeps residue 1, changes residue 2, a proline, to a
/// cysteine and meets the same stop codon: a count of 3, as Ensembl VEP
/// writes it. On the model without the flag whose record shows no start
/// codon, the first deletion is counted as on the flagged model and takes no
/// start term.
///
/// On the model whose start codon stands, no count is written. A copy of the
/// start codon's first base inserted before it is described after that base,
/// and a deletion of the start codon's first base, the base the 5' UTR ends
/// on, leaves the start codon in place: after each the edited transcript
/// holds the whole coding sequence behind a start codon, and the notation is
/// the first residue unchanged. A deletion of the start codon's first two
/// bases loses the start, and the notation is the lost start's.
#[test]
fn a_frameshift_named_at_residue_1_is_counted_only_where_no_start_codon_stands() {
    let dir = tempfile::tempdir().unwrap();
    let flagged = Model::flagged();
    let unflagged = Model::without_a_start_codon();
    let standing = Model::with_a_start_codon();
    let records = [
        flagged.deletion("flagged_first_codon", 1, 2),
        flagged.deletion("flagged_second_codon", 4, 5),
        unflagged.deletion("no_start_codon_shown_first_codon", 1, 2),
        standing.insertion_before_the_coding_start("start_codon_stands_insertion_before_it"),
        standing.deletion("start_codon_stands_first_base", 1, 1),
        standing.deletion("start_codon_stands_first_two_bases", 1, 2),
    ]
    .concat();
    let (cache, fasta) = made_up(dir.path(), &[flagged, unflagged, standing]);
    assert_eq!(
        rows(dir.path(), &cache, &fasta, "GRCh38", &records),
        [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
            "flagged_first_codon\t1:101-102\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t1-2\t1-2\t1\tE/X\tGAa/a\t-\tIMPACT=HIGH;STRAND=1;FLAGS=cds_start_NF;HGVSc=ENST00000000001.1:c.1_2del;HGVSp=ENSP00000000001.1:p.Glu1ThrfsTer4",
            "flagged_second_codon\t1:104-105\t-\tENSG00000000001\tENST00000000001\tTranscript\tframeshift_variant\t4-5\t4-5\t2\tP/X\tCCt/t\t-\tIMPACT=HIGH;STRAND=1;FLAGS=cds_start_NF;HGVSc=ENST00000000001.1:c.4_5del;HGVSp=ENSP00000000001.1:p.Pro2CysfsTer3",
            "no_start_codon_shown_first_codon\t2:101-102\t-\tENSG00000000002\tENST00000000002\tTranscript\tframeshift_variant\t1-2\t1-2\t1\tE/X\tGAa/a\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000002.1:c.1_2del;HGVSp=ENSP00000000002.1:p.Glu1ThrfsTer4",
            "start_codon_stands_insertion_before_it\t3:106-107\tA\tENSG00000000003\tENST00000000003\tTranscript\t5_prime_UTR_variant\t6-7\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000000003.1:c.1dup;HGVSp=ENSP00000000003.1:p.Met1%3D;HGVS_OFFSET=1",
            "start_codon_stands_first_base\t3:107\t-\tENSG00000000003\tENST00000000003\tTranscript\tframeshift_variant,start_retained_variant\t7\t1\t1\tM/X\tAtg/tg\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000003.1:c.1del;HGVSp=ENSP00000000003.1:p.Met1%3D",
            "start_codon_stands_first_two_bases\t3:107-108\t-\tENSG00000000003\tENST00000000003\tTranscript\tframeshift_variant,start_lost\t7-8\t1-2\t1\tM/X\tATg/g\t-\tIMPACT=HIGH;STRAND=1;HGVSc=ENST00000000003.1:c.1_2del;HGVSp=ENSP00000000003.1:p.Met1?",
        ]
    );
}

/// The first record of the entry's reproduction deletes `21:35247664-35247665`,
/// the two bases that open the coding sequence of `ENST00000415023`, which the
/// committed GRCh37 cache flags `cds_start_NF`. The first codon of the new
/// frame names a lysine where the record's first residue is a glutamate, and
/// its fourth codon, coding bases 12 to 14, is a stop codon: a count of 4.
#[test]
fn two_bases_deleted_from_the_first_codon_of_a_flagged_transcript_count_to_the_new_stop_codon() {
    let dir = tempfile::tempdir().unwrap();
    let corpus =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37-hgvs");
    let written = rows(
        dir.path(),
        &corpus.join("json_cache"),
        &corpus.join("reference.fa.gz"),
        "GRCh37",
        "21\t35247663\tfirst_codon_two_bases_deleted\tGGA\tG\t.\t.\t.\n",
    );
    let row: Vec<&String> = written
        .iter()
        .filter(|line| line.split('\t').nth(4) == Some("ENST00000415023"))
        .collect();
    assert_eq!(
        row,
        ["first_codon_two_bases_deleted\t21:35247664-35247665\t-\tENSG00000205726\tENST00000415023\tTranscript\tframeshift_variant\t1-2\t1-2\t1\tE/X\tGAa/a\t-\tIMPACT=HIGH;STRAND=1;FLAGS=cds_start_NF,cds_end_NF;HGVSc=ENST00000415023.1:c.1_2del;HGVSp=ENSP00000409800.1:p.Glu1LysfsTer4"]
    );
}
