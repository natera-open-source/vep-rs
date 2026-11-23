// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A record inside the span of a transcript model and in none of its exons
//! or introns, through the `vep` binary on a synthetic cache.
//!
//! The ready-made caches of some plant genomes hold organelle models spliced
//! from parts on both strands, their exons stored in an order that runs
//! against the sequence. Ensembl builds the intron between two neighbouring
//! exons from the first exon's strand (`Intron.pm:70-75`), so the intron
//! between two such parts has its start above its end, and Ensembl VEP passes
//! over an intron of that shape (`BaseTranscriptVariation.pm:894` and
//! `:1014`). A record between the parts overlaps the model's span, so the
//! model gets a row (`AnnotationType/Transcript.pm:103-108` and `:155-166`),
//! no predicate holds for it, and its allele takes the default term
//! (`BaseVariationFeatureOverlapAllele.pm:285`), `intergenic_variant`
//! (`Utils/Constants.pm:636-649`).
//!
//! The expected rows have the cells Ensembl VEP 116.3 writes for such a
//! record on the mitochondrial sequence of a ready-made plant cache: the
//! model's gene, its id and its strand, the default term, no position, and
//! the impact of that term.
//!
//! Ensembl VEP 116.3 raises no warning for such a record: it passes over the
//! intron and assigns the default term, at the lines above, without one. The
//! run's warnings file therefore holds no line for either record, with
//! `--quiet` and without it.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

use serde_json::json;

mod common;

/// The sequence region of the cache.
const REGION: &str = "1";

/// One part of the model: its genomic span and the strand it lies on.
struct Part {
    start: u64,
    end: u64,
    strand: i8,
}

/// The model's parts in the order the cache stores them: two on the forward
/// strand, the second before the first on the sequence, then two on the
/// reverse strand.
const PARTS: [Part; 4] = [
    Part {
        start: 8001,
        end: 8090,
        strand: 1,
    },
    Part {
        start: 7001,
        end: 7060,
        strand: 1,
    },
    Part {
        start: 2001,
        end: 2090,
        strand: -1,
    },
    Part {
        start: 1001,
        end: 1030,
        strand: -1,
    },
];

/// The intron Ensembl builds between two neighbouring parts, as its start and
/// its end (`Intron.pm:70-75`): the bases after the first part up to the
/// second, read on the first part's strand.
fn intron_between(first: &Part, second: &Part) -> (u64, u64) {
    if first.strand == -1 {
        (second.end + 1, first.start - 1)
    } else {
        (first.end + 1, second.start - 1)
    }
}

/// The cache record of the model: a protein-coding transcript on the forward
/// strand whose coding sequence is its whole cDNA, a start codon, alanine
/// codons and a stop codon.
fn record() -> serde_json::Value {
    let length: u64 = PARTS.iter().map(|part| part.end - part.start + 1).sum();
    let codons = (length / 3) as usize;
    let coding = format!("ATG{}TAA", "GCT".repeat(codons - 2));
    let peptide = format!("M{}", "A".repeat(codons - 2));

    let mut exons = Vec::new();
    let mut pairs = Vec::new();
    let mut cdna = 0;
    for (index, part) in PARTS.iter().enumerate() {
        let part_length = part.end - part.start + 1;
        exons.push(json!({"start": part.start, "end": part.end,
            "rank": index + 1, "phase": 0, "end_phase": 0}));
        pairs.push(
            json!({"from_start": cdna + 1, "from_end": cdna + part_length,
            "to_start": part.start, "to_end": part.end, "ori": part.strand}),
        );
        cdna += part_length;
    }
    let introns: Vec<serde_json::Value> = PARTS
        .windows(2)
        .map(|pair| {
            let (start, end) = intron_between(&pair[0], &pair[1]);
            json!({"start": start, "end": end})
        })
        .collect();

    json!({
        "stable_id": "spliced_model",
        "gene_stable_id": "spliced_gene",
        "chr": REGION,
        "start": PARTS.iter().map(|part| part.start).min().unwrap(),
        "end": PARTS.iter().map(|part| part.end).max().unwrap(),
        "strand": 1,
        "biotype": "protein_coding",
        "source": "ensembl",
        "exons": exons,
        "cdna_coding_start": 1,
        "cdna_coding_end": length,
        "translation": {"stable_id": "spliced_protein", "start": 1, "end": length},
        "variation_effect_feature_cache": {
            "codon_table": 1,
            "translateable_seq": coding,
            "peptide": peptide,
            "introns": introns,
            "mapper": {"start_phase": 0, "cdna_coding_start": 1,
                "cdna_coding_end": length, "pairs": pairs},
        },
    })
}

/// The input: one substitution between the forward parts and the reverse
/// parts, and one between the two forward parts, each further from every part
/// than the splice windows reach.
const INPUT: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
    1\t5000\tbetween_the_strands\tA\tG\t.\t.\t.\n\
    1\t7500\tagainst_the_sequence\tA\tG\t.\t.\t.\n";

/// The default-format rows of both records, whole.
const EXPECTED: &str = "between_the_strands\t1:5000\tG\tspliced_gene\tspliced_model\tTranscript\t\
    intergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1\n\
    against_the_sequence\t1:7500\tG\tspliced_gene\tspliced_model\tTranscript\t\
    intergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1\n";

/// Writes the model's cache and the two records under `dir`, runs the binary
/// on them under `flags` beside the run's own, and returns its output with the
/// path of its output file.
fn annotate(dir: &Path, flags: &[&str]) -> (Output, PathBuf) {
    let cache = dir.join("cache");
    let shard = cache
        .join("transcripts")
        .join(REGION)
        .join("1-1000000.json");
    fs::create_dir_all(shard.parent().unwrap()).unwrap();
    fs::write(&shard, json!([record()]).to_string()).unwrap();
    fs::write(
        cache.join("info.json"),
        json!({"assembly": "GRCh38"}).to_string(),
    )
    .unwrap();
    let input = dir.join("records.vcf");
    fs::write(&input, INPUT).unwrap();
    let out = dir.join("records.out");

    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(&cache)
        .args(["--assembly", "GRCh38", "--offline"])
        .args(["--force_overwrite", "--no_stats", "--no_headers"])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    (output, out)
}

/// A record inside the model's span that no exon and no intron of the model
/// holds gets the model's row with the default term, and no other row.
#[test]
fn record_between_the_parts_of_a_spliced_model_takes_the_default_term_on_its_row() {
    let backwards = PARTS
        .windows(2)
        .map(|pair| intron_between(&pair[0], &pair[1]))
        .filter(|(start, end)| start > end)
        .count();
    assert_eq!(
        backwards, 2,
        "two introns of the model start above their end"
    );

    let dir = tempfile::tempdir().unwrap();
    let (output, out) = annotate(dir.path(), &["--quiet"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(fs::read_to_string(&out).unwrap(), EXPECTED);
}

/// The run raises no warning for either record. The cache carries no checksum
/// manifest, and the warnings file holds that warning's line alone, with
/// `--quiet` and without it; without it no line of standard error names the
/// model.
#[test]
fn record_between_the_parts_of_a_spliced_model_draws_no_warning() {
    for flags in [&[][..], &["--quiet"][..]] {
        let dir = tempfile::tempdir().unwrap();
        let (output, out) = annotate(dir.path(), flags);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{flags:?}: {stderr}");
        assert_eq!(fs::read_to_string(&out).unwrap(), EXPECTED, "{flags:?}");

        let warnings = fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
        assert_eq!(
            warnings,
            format!(
                "WARNING: cache {} carries no checksum manifest\n",
                dir.path().join("cache").display()
            ),
            "{flags:?}"
        );
        let naming_the_model: Vec<&str> = stderr
            .lines()
            .filter(|line| line.contains("spliced_model"))
            .collect();
        assert!(
            naming_the_model.is_empty(),
            "{flags:?}: standard error names the model: {naming_the_model:?}"
        );
    }
}
