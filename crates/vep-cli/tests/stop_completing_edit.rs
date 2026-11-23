// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Records at the stop bases of a transcript whose stop codon an `_rna_edit`
//! completes, through the `vep` binary and the library's annotator on a
//! synthetic cache and reference.
//!
//! A mitochondrial RefSeq model can end one or two bases into its stop codon:
//! the genome holds the codon's first base or two and the transcript's 3' end
//! supplies the rest. The cache records that as an `_rna_edit` replacing the
//! exon's last one or two bases by the whole codon, and its mapper pairs the
//! replaced genomic bases with the codon's last bases: Ensembl's mapper reads an
//! edit that lengthens the cDNA as an insertion at the edit's first base
//! (TranscriptMapper.pm:209-211 after :179) and compares no base (:169). A run
//! with a reference sequence reads which bases of the codon the genome holds
//! and pairs them there, so under the default flags every record on such a
//! transcript reads, in `Consequence`, `cDNA_position`, `CDS_position`,
//! `Protein_position`, `Amino_acids` and `Codons`, as the same change reads on
//! a twin transcript with the same coding sequence whose stop codon is whole
//! on the genome. A run without a reference sequence keeps the cache's pairs.
//!
//! Those six cells under the default flags are what the twin measures. After
//! the held stop bases the twin's genome holds the rest of the codon and the
//! edited transcript's holds the chromosome's next bases, and `--hgvs` and
//! `--shift_3prime 1` carry an insertion 3' along the genome alone, as Ensembl
//! VEP does on these transcripts (TranscriptVariationAllele.pm:142). Under
//! those flags each transcript's row follows its own genome, and their tests
//! hold every insertion, on an edited transcript and on its twin, to that
//! reading.

use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;
use std::sync::Arc;

use serde_json::json;
use vep_cli::variation_source::ReadFailure;
use vep_cli::vcf_parser::parse_vcf_line;
use vep_cli::{AnnotationResources, Annotator, Config, InputVariant, LazyTranscriptIndexes};

mod common;

const CHROMOSOME: &str = "MT";
const STOP: &str = "TAA";
const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

/// The codons a model's coding sequence cycles through between its start
/// codon and its last sense codon, with their residues under the vertebrate
/// mitochondrial code.
const SENSE: [(&str, char); 3] = [("GCT", 'A'), ("CTG", 'L'), ("CCA", 'P')];

/// A single-exon forward-strand coding transcript of the synthetic chromosome.
struct Model {
    id: String,
    /// The transcript's first genomic base.
    start: u64,
    /// The number of sense codons, the start codon among them.
    sense_codons: u64,
    /// How many bases of the stop codon the genome holds: 3 on a twin, 1 or 2
    /// on a transcript whose stop codon an edit completes.
    held: u64,
}

impl Model {
    /// The sense codons with their residues: the start codon, the cycle of
    /// [`SENSE`], then a last codon ending in a base the stop codon lacks.
    fn sense(&self) -> Vec<(&'static str, char)> {
        let mut codons = vec![("ATG", 'M')];
        codons.extend((0..self.sense_codons - 2).map(|i| SENSE[i as usize % SENSE.len()]));
        codons.push(("GAG", 'E'));
        codons
    }

    fn sense_len(&self) -> u64 {
        self.sense_codons * 3
    }

    /// The genomic position of the stop codon's first base.
    fn stop_start(&self) -> u64 {
        self.start + self.sense_len()
    }

    /// The transcript's last genomic base.
    fn end(&self) -> u64 {
        self.stop_start() + self.held - 1
    }

    /// The bases the genome holds for the transcript.
    fn genomic_bases(&self) -> String {
        let sense: String = self.sense().iter().map(|(codon, _)| *codon).collect();
        sense + &STOP[..self.held as usize]
    }

    /// The model's cache record. A transcript whose stop codon an edit
    /// completes carries the edit and the two pairs Ensembl's mapper builds
    /// for it; its translateable sequence holds the whole codon either way.
    fn record(&self) -> serde_json::Value {
        let sense = self.sense();
        let coding: String = sense.iter().map(|(codon, _)| *codon).collect::<String>() + STOP;
        let peptide: String = sense.iter().map(|(_, residue)| *residue).collect();
        let coding_len = self.sense_len() + 3;
        let (pairs, attributes) = if self.held == 3 {
            let whole = json!({"from_start": 1, "from_end": coding_len,
                "to_start": self.start, "to_end": self.end(), "ori": 1});
            (json!([whole]), json!([]))
        } else {
            let before = json!({"from_start": 1, "from_end": self.sense_len(),
                "to_start": self.start, "to_end": self.stop_start() - 1, "ori": 1});
            let replaced = json!({"from_start": coding_len - self.held + 1,
                "from_end": coding_len,
                "to_start": self.stop_start(), "to_end": self.end(), "ori": 1});
            let edit = format!(
                "{} {} {STOP}",
                self.sense_len() + 1,
                self.sense_len() + self.held
            );
            (
                json!([before, replaced]),
                json!([{"code": "_rna_edit", "value": edit}]),
            )
        };
        json!({
            "stable_id": self.id,
            "gene_stable_id": format!("gene_{}", self.id),
            "chr": CHROMOSOME,
            "start": self.start,
            "end": self.end(),
            "strand": 1,
            "biotype": "protein_coding",
            "source": "ensembl",
            "is_canonical": 1,
            "protein": format!("protein_{}", self.id),
            "exons": [{"start": self.start, "end": self.end(), "rank": 1,
                "phase": 0, "end_phase": -1}],
            "cdna_coding_start": 1,
            "cdna_coding_end": coding_len,
            "attributes": attributes,
            "translation": {"stable_id": format!("protein_{}", self.id),
                "start": 1, "end": self.sense_len() + self.held},
            "variation_effect_feature_cache": {
                "codon_table": 2,
                "translateable_seq": coding,
                "peptide": peptide,
                "mapper": {"start_phase": 0, "cdna_coding_start": 1,
                    "cdna_coding_end": coding_len, "pairs": pairs},
            },
        })
    }
}

/// Six transcripts whose stop codon an edit completes, the first from two
/// genomic bases and the others from one, each followed by its twin: the same
/// coding sequence with the stop codon whole on the genome.
fn models() -> Vec<(Model, Model)> {
    (0..6u64)
        .map(|i| {
            let model = |name: &str, start: u64, held: u64| Model {
                id: format!("{name}{i}"),
                start,
                sense_codons: 3 + i,
                held,
            };
            (
                model("edited", 101 + i * 200, if i == 0 { 2 } else { 1 }),
                model("whole", 201 + i * 200, 3),
            )
        })
        .collect()
}

/// The synthetic chromosome: every model's bases at its place, `C` elsewhere.
fn chromosome(models: &[(Model, Model)]) -> Vec<u8> {
    let mut bases = vec![b'C'; 1400];
    for model in models.iter().flat_map(|(edited, whole)| [edited, whole]) {
        let from = model.start as usize - 1;
        let held = model.genomic_bases();
        bases[from..from + held.len()].copy_from_slice(held.as_bytes());
    }
    bases
}

/// Every run of one to three bases, 84 in all.
fn insertions() -> Vec<String> {
    let mut runs: Vec<String> = BASES.iter().map(char::to_string).collect();
    let mut all = runs.clone();
    for _ in 0..2 {
        runs = runs
            .iter()
            .flat_map(|run| BASES.iter().map(move |base| format!("{run}{base}")))
            .collect();
        all.extend(runs.iter().cloned());
    }
    all
}

/// The records of one model as `(shape, position, reference, alternate)`,
/// built from the chromosome's own bases. `held` is how many stop bases the
/// edited twin's genome holds, so both twins get the same shapes: each held
/// stop base substituted and deleted, the held stop bases deleted together,
/// deletions from the last sense base into the stop codon, every insertion
/// between the last sense base and the stop codon and between two held stop
/// bases, the last sense base substituted and deleted, and an insertion after
/// the transcript's last base.
fn records(model: &Model, held: u64, chromosome: &[u8]) -> Vec<(String, u64, String, String)> {
    let span = |from: u64, to: u64| {
        String::from_utf8(chromosome[from as usize - 1..to as usize].to_vec()).unwrap()
    };
    let base = |at: u64| span(at, at);
    let stop = model.stop_start();
    let mut records = Vec::new();
    let mut substitutions = |label: &str, at: u64| {
        for alternate in BASES.iter().filter(|b| b.to_string() != base(at)) {
            records.push((
                format!("{label}_to_{alternate}"),
                at,
                base(at),
                alternate.to_string(),
            ));
        }
    };
    for k in 0..held {
        substitutions(&format!("stop{k}"), stop + k);
    }
    substitutions("sense", stop - 1);
    let mut deletion = |label: &str, first: u64, last: u64| {
        records.push((
            label.to_string(),
            first - 1,
            span(first - 1, last),
            base(first - 1),
        ));
    };
    for k in 0..held {
        deletion(&format!("stop{k}_deleted"), stop + k, stop + k);
    }
    deletion("sense_deleted", stop - 1, stop - 1);
    deletion("seam_deleted", stop - 1, stop);
    if held == 2 {
        deletion("stop_deleted", stop, stop + 1);
        deletion("seam_and_stop_deleted", stop - 1, stop + 1);
    }
    let mut insertion = |label: &str, after: u64, inserted: &str| {
        records.push((
            format!("{label}_{inserted}"),
            after,
            base(after),
            base(after) + inserted,
        ));
    };
    for inserted in insertions() {
        insertion("seam_inserted", stop - 1, &inserted);
        if held == 2 {
            insertion("stop_inserted", stop, &inserted);
        }
    }
    for inserted in BASES {
        insertion("after_inserted", model.end(), &inserted.to_string());
    }
    records
}

/// The cells compared between twins: Consequence, cDNA_position, CDS_position,
/// Protein_position, Amino_acids and Codons.
type Cells = [String; 6];

/// Each record's cells on its own model, keyed by model id and shape.
type Rows = HashMap<(String, String), Cells>;

/// A record's row on its own model: the compared cells and the keys of its
/// `Extra` column.
struct Row {
    cells: Cells,
    extra: HashMap<String, String>,
}

/// Each record's row on its own model, keyed by model id and shape.
type FullRows = HashMap<(String, String), Row>;

/// [`annotate_under`] with no further flag, each row reduced to its cells.
fn annotate(dir: &Path, with_reference: bool) -> (Vec<(Model, Model)>, Rows) {
    let (models, rows) = annotate_under(dir, with_reference, &[]);
    let cells = rows.into_iter().map(|(key, row)| (key, row.cells));
    (models, cells.collect())
}

/// Writes the cache, the reference and every record of every model under
/// `dir`, runs the binary with `flags`, and with `--fasta` when
/// `with_reference`, and returns the models and each record's row on its own
/// model, keyed by model id and shape.
fn annotate_under(
    dir: &Path,
    with_reference: bool,
    flags: &[&str],
) -> (Vec<(Model, Model)>, FullRows) {
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

    let fasta: PathBuf = dir.join("reference.fa");
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
            for (shape, position, reference, alternate) in records(model, edited.held, &bases) {
                vcf.push_str(&format!(
                    "{CHROMOSOME}\t{position}\t{}:{shape}\t{reference}\t{alternate}\t.\t.\t.\n",
                    model.id
                ));
            }
        }
    }
    let name = if with_reference {
        "with_reference"
    } else {
        "without_reference"
    };
    let input = dir.join(format!("{name}.vcf"));
    fs::write(&input, vcf).unwrap();
    let out = dir.join(format!("{name}.out"));
    let mut command = common::isolated_vep(dir);
    command
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(&cache)
        .args(["--assembly", "GRCh38", "--offline"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .args(flags)
        .stdin(Stdio::null());
    if with_reference {
        command.arg("--fasta").arg(&fasta);
    }
    let output = command.output().unwrap();
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
            let extra = cells[13]
                .split(';')
                .filter_map(|pair| pair.split_once('='))
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect();
            let row = Row { cells: kept, extra };
            rows.insert((model.to_string(), shape.to_string()), row);
        }
    }
    (models, rows)
}

/// The insertions of [`records`] written inside `model`, between its last
/// sense base and its stop codon and between two held stop bases, as `(shape,
/// base the insertion follows, inserted bases)`.
fn inside_insertions(model: &Model, held: u64, chromosome: &[u8]) -> Vec<(String, u64, String)> {
    records(model, held, chromosome)
        .into_iter()
        .filter(|(shape, ..)| {
            shape.starts_with("seam_inserted_") || shape.starts_with("stop_inserted_")
        })
        .map(|(shape, after, _, alternate)| (shape, after, alternate[1..].to_string()))
        .collect()
}

/// Where the 3' rule leaves an insertion of `run` after genomic base `after`:
/// the base it then follows and the run as it reads there. The run moves one
/// base while its first base is the chromosome's next.
fn carried(chromosome: &[u8], after: u64, run: &str) -> (u64, String) {
    let mut run = run.as_bytes().to_vec();
    let mut after = after;
    while chromosome.get(after as usize) == Some(&run[0]) {
        run.rotate_left(1);
        after += 1;
    }
    (after, String::from_utf8(run).unwrap())
}

/// The `HGVSc` of an insertion of `run` after genomic base `after` of `model`
/// as the chromosome gives it, at the place the 3' rule leaves it when `shift`
/// and where the record wrote it otherwise; `None` where it has no coding
/// notation. A genomic base of the model has the cDNA position of its place in
/// the exon. An insertion that copies the bases before it, or the bases after
/// it inside the transcript, is a duplication of them; one left after the
/// transcript's last genomic base has no coordinate of its own.
fn coding_notation(
    model: &Model,
    chromosome: &[u8],
    after: u64,
    run: &str,
    shift: bool,
) -> Option<String> {
    let (at, run) = if shift {
        carried(chromosome, after, run)
    } else {
        (after, run.to_string())
    };
    let len = run.len() as u64;
    let copies = |first: u64, last: u64| {
        first >= model.start
            && last <= model.end()
            && &chromosome[first as usize - 1..last as usize] == run.as_bytes()
    };
    let position = |genomic: u64| genomic - model.start + 1;
    let duplication = |first: u64, last: u64| {
        if first == last {
            format!("{}:c.{}dup", model.id, position(first))
        } else {
            format!("{}:c.{}_{}dup", model.id, position(first), position(last))
        }
    };
    if copies(at + 1 - len, at) {
        Some(duplication(at + 1 - len, at))
    } else if copies(at + 1, at + len) {
        Some(duplication(at + 1, at + len))
    } else if at < model.end() {
        let (first, last) = (position(at), position(at) + 1);
        Some(format!("{}:c.{first}_{last}ins{run}", model.id))
    } else {
        None
    }
}

/// The Consequence cells of `rows` on the edited models whose shape starts
/// with `prefix`, counted.
fn consequences(rows: &Rows, prefix: &str) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for ((model, shape), cells) in rows {
        if model.starts_with("edited") && shape.starts_with(prefix) {
            *counts.entry(cells[0].clone()).or_default() += 1;
        }
    }
    counts
}

/// With a reference sequence every record on a transcript whose stop codon an
/// edit completes reads as the same record reads on its twin, at the position
/// the genome's base has in the coding sequence. A substitution of a stop base
/// the genome holds leaves no stop codon under the vertebrate mitochondrial
/// code, so all 21 are `stop_lost`; an insertion of one or two bases between
/// the last sense codon and the stop codon shifts the frame, one of three
/// bases is in frame, and the four that are themselves stop codons keep the
/// stop.
#[test]
fn a_record_at_a_completed_stop_codon_reads_as_on_a_whole_one() {
    let dir = tempfile::tempdir().unwrap();
    let (models, rows) = annotate(dir.path(), true);
    assert_eq!(rows.len(), 2 * 672);

    let mut differing = Vec::new();
    for ((model, shape), cells) in &rows {
        if let Some(pair) = model.strip_prefix("edited") {
            let twin = &rows[&(format!("whole{pair}"), shape.clone())];
            if cells != twin {
                differing.push(format!("{model}:{shape} {cells:?} twin {twin:?}"));
            }
        }
    }
    differing.sort();
    assert!(
        differing.is_empty(),
        "{} of 672 records differ from their twin's:\n{}",
        differing.len(),
        differing.join("\n")
    );

    for (edited, _) in &models {
        for k in 0..edited.held {
            for alternate in BASES
                .iter()
                .filter(|b| !STOP[k as usize..].starts_with(**b))
            {
                let cells = &rows[&(edited.id.clone(), format!("stop{k}_to_{alternate}"))];
                let position = (edited.sense_len() + 1 + k).to_string();
                let codon = (edited.sense_codons + 1).to_string();
                assert_eq!(
                    cells[..4],
                    ["stop_lost".to_string(), position.clone(), position, codon],
                    "{} stop{k}_to_{alternate}",
                    edited.id
                );
            }
        }
    }
    assert_eq!(
        consequences(&rows, "stop0_to_"),
        BTreeMap::from([("stop_lost".to_string(), 18)])
    );
    assert_eq!(
        consequences(&rows, "stop1_to_"),
        BTreeMap::from([("stop_lost".to_string(), 3)])
    );
    assert_eq!(
        consequences(&rows, "seam_inserted_"),
        BTreeMap::from([
            ("frameshift_variant".to_string(), 120),
            ("inframe_insertion".to_string(), 360),
            ("inframe_insertion,stop_retained_variant".to_string(), 24),
        ])
    );
}

/// Without a reference sequence nothing shows which bases of the codon the
/// genome holds, and the cache's pairs stand: a stop base the genome holds is
/// read at the codon's last base on a transcript holding one, and one position
/// on from its own on the transcript holding two.
#[test]
fn without_a_reference_sequence_the_cache_pairs_stand() {
    let dir = tempfile::tempdir().unwrap();
    let (models, rows) = annotate(dir.path(), false);
    for (edited, _) in &models {
        for k in 0..edited.held {
            let alternate = "G";
            let cells = &rows[&(edited.id.clone(), format!("stop{k}_to_{alternate}"))];
            let position = (edited.sense_len() + 4 - edited.held + k).to_string();
            assert_eq!(cells[2], position, "{} stop{k}", edited.id);
        }
    }
}

/// Each record's `HGVSp` on its own model, keyed by model id and shape.
type Notations = HashMap<(String, String), Option<String>>;

/// The models and each record's `HGVSp` under `--hgvs`: the binary reads
/// again, with the flag, the cache, the reference and the records
/// [`annotate`] writes under `dir`.
fn protein_notations(dir: &Path) -> (Vec<(Model, Model)>, Notations) {
    let (models, _) = annotate(dir, true);
    let out = dir.join("protein_notations.out");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(dir.join("with_reference.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(dir.join("cache"))
        .arg("--fasta")
        .arg(dir.join("reference.fa"))
        .args(["--assembly", "GRCh38", "--offline", "--hgvs"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut notations = HashMap::new();
    for line in fs::read_to_string(&out).unwrap().lines() {
        let cells: Vec<&str> = line.split('\t').collect();
        let Some((model, shape)) = cells[0].split_once(':') else {
            continue;
        };
        if cells[4] == model {
            let notation = cells[13]
                .split(';')
                .find_map(|pair| pair.strip_prefix("HGVSp="))
                .map(str::to_string);
            notations.insert((model.to_string(), shape.to_string()), notation);
        }
    }
    (models, notations)
}

/// The genomic base an insertion of `run` after base `after` follows where
/// the 3' rule leaves it: the run moves one base while its first base is the
/// chromosome's next.
fn base_before_shifted(chromosome: &[u8], after: u64, run: &str) -> u64 {
    let mut run = run.as_bytes().to_vec();
    let mut after = after;
    while chromosome.get(after as usize) == Some(&run[0]) {
        run.rotate_left(1);
        after += 1;
    }
    after
}

/// Under `--hgvs` an insertion's protein notation is read where the 3' rule
/// leaves it on the genome, and a place after the transcript's last genomic
/// base has no translation coordinates (TranscriptVariationAllele.pm:1730-1737,
/// TranscriptMapper.pm:441-443). A transcript whose stop codon an edit
/// completes ends on the genome a base or two before its coding sequence
/// does, and the rule leaves 100 of its insertions directly after that last
/// genomic base: none carries an `HGVSp`, as none of its twins' 11 does. The 6
/// that insert one base are the allele a record writes after the last base,
/// and both spellings read alike.
#[test]
fn an_insertion_left_after_the_last_genomic_base_has_no_protein_notation() {
    let dir = tempfile::tempdir().unwrap();
    let (models, notations) = protein_notations(dir.path());
    let bases = chromosome(&models);
    let mut left_after_the_last_base = [0, 0];
    let mut written_both_ways = 0;
    let mut noted_inside = 0;
    let mut noted_outside = Vec::new();
    for (edited, whole) in &models {
        for (side, model) in [edited, whole].into_iter().enumerate() {
            for (shape, after, _, alternate) in records(model, edited.held, &bases) {
                if !shape.contains("_inserted_") {
                    continue;
                }
                let run = &alternate[1..];
                let notation = &notations[&(model.id.clone(), shape.clone())];
                let at = base_before_shifted(&bases, after, run);
                if at < model.end() {
                    noted_inside += usize::from(notation.is_some());
                    continue;
                }
                if notation.is_some() {
                    noted_outside.push(format!("{}:{shape}", model.id));
                }
                if after < model.end() && at == model.end() {
                    left_after_the_last_base[side] += 1;
                    if side == 0 && run.len() == 1 {
                        let written_after = (model.id.clone(), format!("after_inserted_{run}"));
                        written_both_ways += usize::from(notation == &notations[&written_after]);
                    }
                }
            }
        }
    }
    noted_outside.sort();
    assert!(
        noted_outside.is_empty(),
        "{} insertions left after the last genomic base carry an HGVSp:\n{}",
        noted_outside.len(),
        noted_outside.join("\n")
    );
    assert_eq!(left_after_the_last_base, [100, 11]);
    assert_eq!(written_both_ways, 6);
    assert!(noted_inside > 0);
}

/// The library loads every chromosome when its annotator is built, the binary
/// each chromosome as a record first names it; both pair a completed stop
/// codon from the reference sequence, so the library reads a substitution of
/// each stop base the genome holds at the coding position the binary prints.
#[test]
fn the_library_reads_a_held_stop_base_where_the_binary_does() {
    let dir = tempfile::tempdir().unwrap();
    let (models, rows) = annotate(dir.path(), true);
    let bases = chromosome(&models);
    let cache = dir.path().join("cache");
    let fasta = dir.path().join("reference.fa");
    let config = Config::from_argv([
        "vep",
        "--json_cache",
        cache.to_str().unwrap(),
        "--fasta",
        fasta.to_str().unwrap(),
        "--assembly",
        "GRCh38",
        "--offline",
        "--fork",
        "1",
    ])
    .unwrap();
    let annotator = Annotator::from_config(&config).unwrap();
    for (edited, _) in &models {
        for k in 0..edited.held {
            let position = edited.stop_start() + k;
            let reference = bases[position as usize - 1] as char;
            let line = format!("{CHROMOSOME}\t{position}\t.\t{reference}\tG\t.\t.\t.");
            let mut record = parse_vcf_line(&line, false).unwrap();
            annotator.annotate_batch(&mut record);
            let row = record[0]
                .transcript_consequences
                .iter()
                .find(|row| row.transcript_id() == edited.id)
                .unwrap();
            let binary = &rows[&(edited.id.clone(), format!("stop{k}_to_G"))];
            assert_eq!(
                row.cds_position.as_deref(),
                Some(binary[2].as_str()),
                "{} stop{k}",
                edited.id
            );
            assert_eq!(binary[2], (edited.sense_len() + 1 + k).to_string());
        }
    }
}

/// Under `--hgvs` an insertion's notation is read where the 3' rule leaves it
/// on the genome, on an edited transcript and on its twin alike. One left
/// inside the transcript is a duplication of the bases before it or an
/// insertion between two coding positions. One left directly after the
/// transcript's last genomic base is a duplication of the transcript's last
/// bases or has no notation, since a place after the last exon has no
/// coordinate (TranscriptVariationAllele.pm:1524-1527 with :2834), and one
/// carried further has none (:1472). An edited transcript's genome holds one
/// or two stop bases, so 111 of the edited transcripts' insertions have no
/// notation and 3 of their twins'.
///
/// Under `--shift_hgvs 0` an insertion stays where the record wrote it and is
/// a duplication when it copies the bases before it or the bases after it
/// inside the transcript. The 12 insertions that copy the stop bases an edit
/// supplies are duplications on a twin and insertions on an edited
/// transcript, whose genome lacks those bases.
#[test]
fn under_hgvs_an_insertion_is_noted_where_the_genome_carries_it() {
    let dir = tempfile::tempdir().unwrap();
    let (models, shifted) = annotate_under(dir.path(), true, &["--hgvs"]);
    let (_, unshifted) = annotate_under(dir.path(), true, &["--hgvs", "--shift_hgvs", "0"]);
    let bases = chromosome(&models);
    let mut without_notation = [0, 0];
    let mut insertions_the_twin_duplicates = 0;
    for (edited, whole) in &models {
        for (side, model) in [edited, whole].into_iter().enumerate() {
            for (shape, after, run) in inside_insertions(model, edited.held, &bases) {
                let key = (model.id.clone(), shape.clone());
                let (at, _) = carried(&bases, after, &run);
                let notation = coding_notation(model, &bases, after, &run, true);
                let offset = notation
                    .as_ref()
                    .filter(|_| at > after)
                    .map(|_| (at - after).to_string());
                let row = &shifted[&key];
                assert_eq!(row.extra.get("HGVSc"), notation.as_ref(), "{key:?}");
                assert_eq!(row.extra.get("HGVS_OFFSET"), offset.as_ref(), "{key:?}");
                without_notation[side] += usize::from(notation.is_none());

                let notation = coding_notation(model, &bases, after, &run, false);
                let row = &unshifted[&key];
                assert_eq!(row.extra.get("HGVSc"), notation.as_ref(), "{key:?}");
                assert_eq!(row.extra.get("HGVS_OFFSET"), None, "{key:?}");
            }
        }
        for (shape, ..) in inside_insertions(edited, edited.held, &bases) {
            let notation =
                |model: &Model| &unshifted[&(model.id.clone(), shape.clone())].extra["HGVSc"];
            if notation(edited).contains("ins") && notation(whole).ends_with("dup") {
                insertions_the_twin_duplicates += 1;
            }
        }
    }
    assert_eq!(without_notation, [111, 3]);
    assert_eq!(insertions_the_twin_duplicates, 12);
}

/// Under `--shift_3prime 1` a row describes the place the 3' rule carries its
/// record to on the genome: the span moves by the shift before the overlap
/// with the transcript is tested (BaseVariationFeatureOverlapAllele.pm:460-462
/// and :480). An insertion left at or after the transcript's last genomic
/// base reads `downstream_gene_variant`, 131 on the edited transcripts and 12
/// on their twins; on 119 of the 131 the twin's genome holds the bases the
/// shift reads and its row keeps a coding term. An insertion left inside an
/// edited transcript reads as on its twin.
#[test]
fn under_shift_3prime_an_insertion_carried_past_the_last_base_is_downstream() {
    let dir = tempfile::tempdir().unwrap();
    let (models, rows) = annotate_under(dir.path(), true, &["--shift_3prime", "1"]);
    let bases = chromosome(&models);
    let mut downstream = [0, 0];
    let mut coding_on_the_twin = 0;
    for (edited, whole) in &models {
        for (side, model) in [edited, whole].into_iter().enumerate() {
            for (shape, after, run) in inside_insertions(model, edited.held, &bases) {
                let key = (model.id.clone(), shape);
                let (at, _) = carried(&bases, after, &run);
                let is_downstream = rows[&key].cells[0] == "downstream_gene_variant";
                assert_eq!(is_downstream, at >= model.end(), "{key:?}");
                downstream[side] += usize::from(is_downstream);
            }
        }
        for (shape, after, run) in inside_insertions(edited, edited.held, &bases) {
            let row = &rows[&(edited.id.clone(), shape.clone())];
            let twin = &rows[&(whole.id.clone(), shape.clone())];
            if carried(&bases, after, &run).0 < edited.end() {
                assert_eq!(row.cells, twin.cells, "{} {shape}", edited.id);
            } else if twin.cells[0] != "downstream_gene_variant" {
                coding_on_the_twin += 1;
            }
        }
    }
    assert_eq!(downstream, [131, 12]);
    assert_eq!(coding_on_the_twin, 119);
}

/// The reference sequence at `fasta`, as a library caller opens it.
fn reference(fasta: &Path) -> Arc<vep_fasta::IndexedFasta> {
    Arc::new(vep_fasta::IndexedFasta::from_path(fasta).unwrap())
}

/// The resources a library caller assembles by hand: its own index and, in
/// the effects configuration, the reference sequence at `fasta`.
fn caller_resources(transcripts: LazyTranscriptIndexes, fasta: &Path) -> AnnotationResources {
    AnnotationResources::new(
        Arc::new(transcripts),
        vep_cli::annotator::build_effects_config(
            (5_000, 5_000),
            Some(reference(fasta)),
            false,
            Default::default(),
            false,
            false,
            false,
        ),
    )
}

/// An annotator over `transcripts` as an embedder assembles one, from
/// [`caller_resources`].
fn embedder_annotator(transcripts: LazyTranscriptIndexes, fasta: &Path) -> Annotator {
    let resources = caller_resources(transcripts, fasta);
    Annotator::new(resources, vep_cli::build_thread_pool(1).unwrap()).unwrap()
}

/// `cells` with the consequence terms in alphabetical order, the order
/// [`library_cells`] gives them.
fn with_sorted_terms(cells: &Cells) -> Cells {
    let mut terms: Vec<&str> = cells[0].split(',').collect();
    terms.sort_unstable();
    let mut sorted = cells.clone();
    sorted[0] = terms.join(",");
    sorted
}

/// The compared cells of the VCF record `line` on the transcript `id` as
/// `annotate` reads them, an absent value as the `-` the binary prints and
/// the consequence terms in alphabetical order.
fn library_cells(annotate: impl Fn(&mut Vec<InputVariant>), id: &str, line: &str) -> Cells {
    let mut record = parse_vcf_line(line, false).unwrap();
    annotate(&mut record);
    let row = record[0]
        .transcript_consequences
        .iter()
        .find(|row| row.transcript_id() == id)
        .unwrap();
    let mut terms: Vec<&str> = row.consequences.iter().map(|c| c.so_term()).collect();
    terms.sort_unstable();
    let cell = |value: &Option<String>| value.clone().unwrap_or_else(|| "-".to_string());
    [
        terms.join(","),
        cell(&row.cdna_position),
        cell(&row.cds_position),
        cell(&row.protein_position),
        cell(&row.amino_acids),
        cell(&row.codons),
    ]
}

/// An embedder that builds the transcript index itself, lazily with the
/// chromosome warmed or from transcripts it loaded, and hands its annotator
/// the reference sequence reads every record on an edited transcript as the
/// binary reads it with `--fasta`: the annotator pairs a completed stop codon
/// from the reference sequence it is given, whichever way its index was built.
#[test]
fn an_embedder_with_the_reference_reads_every_record_as_the_binary_does() {
    let dir = tempfile::tempdir().unwrap();
    let (models, rows) = annotate(dir.path(), true);
    let bases = chromosome(&models);
    let cache = dir.path().join("cache");
    let cache = cache.to_str().unwrap();
    let fasta = dir.path().join("reference.fa");

    let warmed = LazyTranscriptIndexes::new(cache, 5_000, 5_000).unwrap();
    warmed.prewarm(CHROMOSOME).unwrap();
    let loaded = vep_cli::json_cache::load_all_transcripts(cache).unwrap();
    let loaded = LazyTranscriptIndexes::from_loaded(loaded, 5_000, 5_000);
    let mut differing = Vec::new();
    for (construction, transcripts) in [("warmed", warmed), ("loaded", loaded)] {
        let annotator = embedder_annotator(transcripts, &fasta);
        let mut compared = 0;
        let mut rows_differing = Vec::new();
        for (edited, _) in &models {
            for (shape, position, reference, alternate) in records(edited, edited.held, &bases) {
                let line =
                    format!("{CHROMOSOME}\t{position}\t.\t{reference}\t{alternate}\t.\t.\t.");
                let annotate = |record: &mut Vec<InputVariant>| {
                    annotator.annotate_batch(record);
                };
                let cells = library_cells(annotate, &edited.id, &line);
                let binary = with_sorted_terms(&rows[&(edited.id.clone(), shape.clone())]);
                compared += 1;
                if cells != binary {
                    rows_differing
                        .push(format!("{}:{shape} {cells:?} binary {binary:?}", edited.id));
                }
            }
        }
        assert_eq!(compared, 672);
        if !rows_differing.is_empty() {
            differing.push(format!(
                "{construction}: {} of 672 records differ from the binary's:\n{}",
                rows_differing.len(),
                rows_differing.join("\n")
            ));
        }
    }
    assert!(differing.is_empty(), "{}", differing.join("\n"));
}

/// `--total_length` prints the cDNA's length beside `cDNA_position`: the
/// exons and the bases an `_rna_edit` adds together, the length of the
/// sequence the position counts along. An edited transcript's total is its
/// twin's, on these transcripts the coding sequence's length, and its row
/// reads as its twin's in every compared cell. Ensembl VEP's total is the
/// exons' length alone, since the transcript its cache stores carries no
/// `edits_enabled` and `Transcript::length` then leaves the edits out
/// (OutputFactory.pm:1737; Transcript.pm:1675-1678).
#[test]
fn total_length_counts_the_bases_an_edit_adds() {
    let dir = tempfile::tempdir().unwrap();
    let (models, rows) = annotate_under(dir.path(), true, &["--total_length"]);
    let bases = chromosome(&models);
    let mut with_total = [0, 0];
    let mut differing = Vec::new();
    for (edited, whole) in &models {
        let length = (edited.sense_len() + 3).to_string();
        for (side, model) in [edited, whole].into_iter().enumerate() {
            for (shape, ..) in records(model, edited.held, &bases) {
                let cells = &rows[&(model.id.clone(), shape.clone())].cells;
                if let Some((_, total)) = cells[1].split_once('/') {
                    with_total[side] += 1;
                    if total != length {
                        differing.push(format!("{}:{shape} {} of {length}", model.id, cells[1]));
                    }
                }
                let twin = &rows[&(whole.id.clone(), shape.clone())].cells;
                if cells != twin {
                    differing.push(format!("{}:{shape} {cells:?} twin {twin:?}", model.id));
                }
            }
        }
    }
    assert_eq!(with_total, [648, 648]);
    differing.sort();
    assert!(
        differing.is_empty(),
        "{} rows differ:\n{}",
        differing.len(),
        differing.join("\n")
    );
}

/// A library caller that builds its own transcript index and annotates
/// through it without an [`Annotator`] hands the index the reference sequence
/// with `with_reference`, and reads every record on an edited transcript as
/// the binary reads it with `--fasta`: on a set handed the reference before
/// its chromosome loads, on one handed it afterwards and on one built from
/// transcripts the caller loaded. The same index without the reference
/// sequence keeps the cache's pairs, under which 624 of the 672 records read
/// otherwise.
#[test]
fn a_caller_that_hands_its_index_the_reference_reads_as_the_binary_does() {
    let dir = tempfile::tempdir().unwrap();
    let (models, rows) = annotate(dir.path(), true);
    let bases = chromosome(&models);
    let cache = dir.path().join("cache");
    let cache = cache.to_str().unwrap();
    let fasta = dir.path().join("reference.fa");
    let lazy = || LazyTranscriptIndexes::new(cache, 5_000, 5_000).unwrap();
    let loaded = || {
        let transcripts = vep_cli::json_cache::load_all_transcripts(cache).unwrap();
        LazyTranscriptIndexes::from_loaded(transcripts, 5_000, 5_000)
    };
    let warmed = |transcripts: LazyTranscriptIndexes| {
        transcripts.prewarm(CHROMOSOME).unwrap();
        transcripts
    };
    let handed =
        |transcripts: LazyTranscriptIndexes| transcripts.with_reference(Some(reference(&fasta)));
    let differing = |transcripts: LazyTranscriptIndexes| {
        let resources = caller_resources(transcripts, &fasta);
        let annotate = |record: &mut Vec<InputVariant>| {
            let failure = ReadFailure::default();
            vep_cli::annotator::annotate_batch(record, &resources, false, &failure);
        };
        let mut differing = 0;
        for (edited, _) in &models {
            for (shape, position, reference, alternate) in records(edited, edited.held, &bases) {
                let line =
                    format!("{CHROMOSOME}\t{position}\t.\t{reference}\t{alternate}\t.\t.\t.");
                let cells = library_cells(annotate, &edited.id, &line);
                let binary = with_sorted_terms(&rows[&(edited.id.clone(), shape)]);
                differing += usize::from(cells != binary);
            }
        }
        differing
    };
    assert_eq!(
        differing(warmed(handed(lazy()))),
        0,
        "handed before its load"
    );
    assert_eq!(
        differing(handed(warmed(lazy()))),
        0,
        "handed after its load"
    );
    assert_eq!(
        differing(handed(loaded())),
        0,
        "built from loaded transcripts"
    );
    assert_eq!(differing(warmed(lazy())), 624, "never handed the reference");
    assert_eq!(
        differing(loaded()),
        624,
        "loaded and never handed the reference"
    );
}
