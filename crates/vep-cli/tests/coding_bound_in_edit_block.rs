// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Records on a transcript whose first or last coding base lies in bases an
//! `_rna_edit` adds to the model and the genome lacks, through the `vep`
//! binary on a synthetic cache and reference.
//!
//! A RefSeq model can carry a block of bases before its first genomic base or
//! after its last, and its coding sequence can begin or end inside the block.
//! No mapper pair holds that coding bound. Ensembl VEP numbers a coding
//! position from the cDNA coding start alone (ensembl core
//! `TranscriptMapper.pm:421` and `:447-448`) and takes the genomic coding
//! region from the translation's exons (ensembl core
//! `Transcript.pm:1087-1154`), so a record on a coding base the genome holds
//! carries its `CDS_position`, `Protein_position`, `Amino_acids` and `Codons`.
//! Each record here reads on the edited model as it reads on a twin whose
//! whole cDNA lies on the genome.

use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;

mod common;

const CHROMOSOME: &str = "1";
const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

/// The coding sequence every model shares, as codons with their residues: a
/// start codon, nine sense codons and a stop codon, 33 bases.
const CODING: [(&str, char); 11] = [
    ("ATG", 'M'),
    ("GCT", 'A'),
    ("CTG", 'L'),
    ("CCA", 'P'),
    ("GAG", 'E'),
    ("GCT", 'A'),
    ("CTG", 'L'),
    ("CCA", 'P'),
    ("GAG", 'E'),
    ("AAA", 'K'),
    ("TAA", '*'),
];

fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        _ => 'A',
    }
}

/// A single-exon coding transcript of the synthetic chromosome.
struct Model {
    id: String,
    /// 1 or -1.
    strand: i8,
    /// The lowest genomic base the model occupies.
    start: u64,
    /// The 5' UTR, empty or three bases.
    utr5: &'static str,
    /// The 3' UTR, empty or six bases.
    utr3: &'static str,
    /// How many bases at the cDNA's 5' end the genome lacks.
    prepended: usize,
    /// How many bases at the cDNA's 3' end the genome lacks.
    appended: usize,
}

impl Model {
    fn coding(&self) -> String {
        CODING.iter().map(|(codon, _)| *codon).collect()
    }

    fn cdna(&self) -> String {
        format!("{}{}{}", self.utr5, self.coding(), self.utr3)
    }

    /// The cDNA bases the genome holds, in transcript order.
    fn held(&self) -> String {
        let cdna = self.cdna();
        cdna[self.prepended..cdna.len() - self.appended].to_string()
    }

    fn end(&self) -> u64 {
        self.start + self.held().len() as u64 - 1
    }

    /// The bases the model occupies on the chromosome, forward strand.
    fn genomic_bases(&self) -> String {
        if self.strand == 1 {
            self.held()
        } else {
            self.held().chars().rev().map(complement).collect()
        }
    }

    /// The genomic position of the 1-based cDNA position `cdna`, a base the
    /// genome holds.
    fn genomic(&self, cdna: usize) -> u64 {
        let offset = (cdna - self.prepended - 1) as u64;
        if self.strand == 1 {
            self.start + offset
        } else {
            self.end() - offset
        }
    }

    /// The first and last coding cDNA positions the genome holds.
    fn held_coding(&self) -> (usize, usize) {
        let coding_start = self.utr5.len() + 1;
        let coding_end = self.utr5.len() + self.coding().len();
        (
            coding_start.max(self.prepended + 1),
            coding_end.min(self.cdna().len() - self.appended),
        )
    }

    /// The model's cache record: one mapper pair over the bases the genome
    /// holds, and an `_rna_edit` inserting each block the genome lacks, its
    /// span counted in the unedited cDNA.
    fn record(&self) -> serde_json::Value {
        let cdna = self.cdna();
        let held = self.held().len();
        let coding_start = self.utr5.len() + 1;
        let coding_end = self.utr5.len() + self.coding().len();
        let mut attributes = Vec::new();
        if self.prepended > 0 {
            let block = &cdna[..self.prepended];
            attributes.push(json!({"code": "_rna_edit", "value": format!("1 0 {block}")}));
        }
        if self.appended > 0 {
            let block = &cdna[cdna.len() - self.appended..];
            let edit = format!("{} {held} {block}", held + 1);
            attributes.push(json!({"code": "_rna_edit", "value": edit}));
        }
        let peptide: String = CODING[..CODING.len() - 1]
            .iter()
            .map(|(_, residue)| *residue)
            .collect();
        let phase = if self.utr5.is_empty() { 0 } else { -1 };
        let mut cache = json!({
            "codon_table": 1,
            "translateable_seq": self.coding(),
            "peptide": peptide,
            "mapper": {"start_phase": phase, "cdna_coding_start": coding_start,
                "cdna_coding_end": coding_end,
                "pairs": [{"from_start": self.prepended + 1,
                    "from_end": cdna.len() - self.appended,
                    "to_start": self.start, "to_end": self.end(), "ori": self.strand}]},
        });
        if !self.utr5.is_empty() {
            cache["five_prime_utr"] = json!(self.utr5);
        }
        if !self.utr3.is_empty() {
            cache["three_prime_utr"] = json!(self.utr3);
        }
        json!({
            "stable_id": self.id,
            "gene_stable_id": format!("gene_{}", self.id),
            "chr": CHROMOSOME,
            "start": self.start,
            "end": self.end(),
            "strand": self.strand,
            "biotype": "protein_coding",
            "source": "ensembl",
            "exons": [{"start": self.start, "end": self.end(), "rank": 1,
                "phase": phase, "end_phase": -1}],
            "cdna_coding_start": coding_start,
            "cdna_coding_end": coding_end,
            "attributes": attributes,
            "translation": {"stable_id": format!("protein_{}", self.id),
                "start": coding_start.saturating_sub(self.prepended).max(1),
                "end": coding_end.min(cdna.len() - self.appended) - self.prepended},
            "variation_effect_feature_cache": cache,
        })
    }
}

/// Five edited models, each followed by its twin with the same cDNA wholly on
/// the genome: the coding start in a prepended block on either strand, the
/// coding end in an appended block on either strand, and both at once.
fn models() -> Vec<(Model, Model)> {
    let shapes: [(&str, i8, &'static str, &'static str, usize, usize); 5] = [
        ("start_forward", 1, "", "GGTCAT", 12, 0),
        ("start_reverse", -1, "", "GGTCAT", 7, 0),
        ("end_forward", 1, "TCG", "", 0, 5),
        ("end_reverse", -1, "TCG", "", 0, 5),
        ("both_forward", 1, "", "", 12, 5),
    ];
    shapes
        .into_iter()
        .enumerate()
        .map(|(i, (name, strand, utr5, utr3, prepended, appended))| {
            let model = |id: String, start: u64, prepended, appended| Model {
                id,
                strand,
                start,
                utr5,
                utr3,
                prepended,
                appended,
            };
            let at = 1001 + i as u64 * 400;
            (
                model(name.to_string(), at, prepended, appended),
                model(format!("whole_{name}"), at + 200, 0, 0),
            )
        })
        .collect()
}

/// The synthetic chromosome: every model's bases at its place, `C` elsewhere.
fn chromosome(models: &[(Model, Model)]) -> Vec<u8> {
    let mut bases = vec![b'C'; 3200];
    for model in models.iter().flat_map(|(edited, whole)| [edited, whole]) {
        let from = model.start as usize - 1;
        let held = model.genomic_bases();
        bases[from..from + held.len()].copy_from_slice(held.as_bytes());
    }
    bases
}

/// The records of `model` at the coding bases its edited twin's genome holds,
/// `first..=last` in cDNA positions, as `(shape, position, reference,
/// alternate)`: a substitution of the first, the middle and the last such
/// base, and a one-base deletion and insertion at the middle one. Each change
/// is the same change of the transcript on both twins.
fn records(
    model: &Model,
    (first, last): (usize, usize),
    chromosome: &[u8],
) -> Vec<(String, u64, String, String)> {
    let base = |at: u64| chromosome[at as usize - 1] as char;
    let cdna = model.cdna();
    let on_genome = |transcript_base: char| {
        if model.strand == 1 {
            transcript_base
        } else {
            complement(transcript_base)
        }
    };
    let middle = (first + last) / 2;
    let mut records = Vec::new();
    for (label, at) in [("first", first), ("middle", middle), ("last", last)] {
        let position = model.genomic(at);
        let reference = cdna.as_bytes()[at - 1] as char;
        let alternate = *BASES.iter().find(|b| **b != reference).unwrap();
        records.push((
            format!("{label}_substituted"),
            position,
            base(position).to_string(),
            on_genome(alternate).to_string(),
        ));
    }
    let deleted = model.genomic(middle);
    records.push((
        "middle_deleted".to_string(),
        deleted - 1,
        format!("{}{}", base(deleted - 1), base(deleted)),
        base(deleted - 1).to_string(),
    ));
    // Between the middle base and the next one of the transcript.
    let anchor = model.genomic(middle).min(model.genomic(middle + 1));
    records.push((
        "middle_inserted".to_string(),
        anchor,
        base(anchor).to_string(),
        format!("{}{}", base(anchor), on_genome('T')),
    ));
    records
}

/// Consequence, cDNA_position, CDS_position, Protein_position, Amino_acids and
/// Codons.
type Cells = [String; 6];

/// Each record's cells on its own model, keyed by model id and shape.
type Rows = HashMap<(String, String), Cells>;

/// Writes the cache, the reference and every record under `dir`, runs the
/// binary, and returns the models and each record's cells on its own model,
/// keyed by model id and shape.
fn annotate(dir: &Path) -> (Vec<(Model, Model)>, Rows) {
    let models = models();
    let bases = chromosome(&models);

    let cache = dir.join("cache");
    let shard = cache
        .join("transcripts")
        .join(CHROMOSOME)
        .join("1-1000000.json");
    fs::create_dir_all(shard.parent().unwrap()).unwrap();
    let transcripts: Vec<serde_json::Value> = models
        .iter()
        .flat_map(|(edited, whole)| [edited.record(), whole.record()])
        .collect();
    fs::write(&shard, serde_json::Value::Array(transcripts).to_string()).unwrap();
    fs::write(
        cache.join("info.json"),
        json!({"assembly": "GRCh38"}).to_string(),
    )
    .unwrap();

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

    let mut vcf =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (edited, whole) in &models {
        for model in [edited, whole] {
            for (shape, position, reference, alternate) in
                records(model, edited.held_coding(), &bases)
            {
                vcf.push_str(&format!(
                    "{CHROMOSOME}\t{position}\t{}:{shape}\t{reference}\t{alternate}\t.\t.\t.\n",
                    model.id
                ));
            }
        }
    }
    let input = dir.join("records.vcf");
    fs::write(&input, vcf).unwrap();
    let out = dir.join("records.out");
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
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let mut rows = HashMap::new();
    for line in fs::read_to_string(&out).unwrap().lines() {
        let cells: Vec<&str> = line.split('\t').collect();
        let Some((model, shape)) = cells[0].split_once(':') else {
            continue;
        };
        if cells[4] == model {
            let kept: Cells = std::array::from_fn(|i| cells[6 + i].to_string());
            rows.insert((model.to_string(), shape.to_string()), kept);
        }
    }
    (models, rows)
}

/// Every record on a coding base the genome holds reads on the edited model
/// as on its twin, in the term, the three positions, the residues and the
/// codons; the coding position counts the bases of the block before it
/// (`cds = cdna - cdna_coding_start + 1`, TranscriptMapper.pm:447-448).
#[test]
fn a_coding_bound_in_an_edit_block_keeps_the_coding_columns() {
    let dir = tempfile::tempdir().unwrap();
    let (models, rows) = annotate(dir.path());
    assert_eq!(rows.len(), 2 * 5 * 5);

    let mut differing = Vec::new();
    for (edited, whole) in &models {
        for ((model, shape), cells) in &rows {
            if *model != edited.id {
                continue;
            }
            let twin = &rows[&(whole.id.clone(), shape.clone())];
            if cells != twin {
                differing.push(format!("{model}:{shape} {cells:?} twin {twin:?}"));
            }
        }
    }
    differing.sort();
    assert!(
        differing.is_empty(),
        "{} of 25 records differ from their twin's:\n{}",
        differing.len(),
        differing.join("\n")
    );

    // CDS_position and Protein_position of the first and last coding bases
    // the genome holds, by the model: the block's bases count.
    let expected = [
        ("start_forward", "first", "13", "5"),
        ("start_forward", "last", "33", "11"),
        ("start_reverse", "first", "8", "3"),
        ("end_forward", "first", "1", "1"),
        ("end_forward", "last", "28", "10"),
        ("end_reverse", "last", "28", "10"),
        ("both_forward", "first", "13", "5"),
        ("both_forward", "last", "28", "10"),
    ];
    for (model, base, cds, protein) in expected {
        let cells = &rows[&(model.to_string(), format!("{base}_substituted"))];
        assert_eq!(
            [cells[2].as_str(), cells[3].as_str()],
            [cds, protein],
            "{model} {base}"
        );
        assert!(
            cells[4] != "-" && cells[5] != "-",
            "{model} {base}: residues and codons are written"
        );
    }
}
