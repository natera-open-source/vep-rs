// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Two sequence regions of one cache that name different genetic codes,
//! through the `vep` binary on a synthetic cache.
//!
//! Ensembl VEP keeps the codon table on each transcript's record, taken from
//! its sequence region's `codon_table` attribute
//! (`AnnotationSource/Database/Transcript.pm:359-363`), and translates an
//! allele's codon under the table of its own transcript
//! (`BaseTranscriptVariation.pm:1146-1164`, `TranscriptVariationAllele.pm:712`
//! and `:736`). A cache that holds one region on one code and another region
//! on another therefore writes, for the same substitution in the same coding
//! sequence, the amino acids of each region's own code.
//!
//! Every codon of the coding sequence is computed with `vep_core::codon` over
//! the 64 codons in the order of `Bio::Tools::CodonTable`'s strings, so each
//! expected amino acid is the one that module's table gives.
//!
//! Genetic code 26, the nuclear code of Pachysolen tannophilus, is no table of
//! that module: its test pins whole rows written out by hand from NCBI's
//! genetic code table.

use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::process::Stdio;

use serde_json::json;
use vep_core::codon::{is_start_codon, is_ter_codon, translate_codon_with_table};

mod common;

type Codon = [u8; 3];

/// The table id of the standard genetic code.
const STANDARD: u8 = 1;

/// The two sequence regions of the cache.
const REGIONS: [&str; 2] = ["1", "2"];

/// The first genomic base of each region's transcript.
const TRANSCRIPT_START: u64 = 1001;

/// The 64 codons in the order of `Bio::Tools::CodonTable`'s strings: the first
/// base the most significant, the bases counted thymine, cytosine, adenine,
/// guanine.
fn codons_in_table_order() -> Vec<Codon> {
    let bases = *b"TCAG";
    let mut codons = Vec::with_capacity(64);
    for first in bases {
        for second in bases {
            for third in bases {
                codons.push([first, second, third]);
            }
        }
    }
    codons
}

/// The first codon in the module's order that satisfies `wanted`.
fn first_codon(wanted: impl Fn(&Codon) -> bool) -> Option<Codon> {
    codons_in_table_order()
        .into_iter()
        .find(|codon| wanted(codon))
}

/// The amino acid of `codon` under `table`.
fn amino_acid(codon: &Codon, table: u8) -> char {
    translate_codon_with_table(codon, table) as char
}

/// The first codon in the module's order that `table` and `other` translate
/// to different amino acids; none when the two translate all 64 codons alike.
fn first_codon_read_differently(table: u8, other: u8) -> Option<Codon> {
    first_codon(|codon| amino_acid(codon, table) != amino_acid(codon, other))
}

/// A one-base substitution in a codon of the shared coding sequence.
struct Substitution {
    /// The codon the coding sequence holds.
    reference: Codon,
    /// The codon the substitution leaves, one base away.
    alternate: Codon,
}

impl Substitution {
    /// The place in the codon, from 0, of the base that changes.
    fn offset(&self) -> usize {
        (0..3)
            .find(|at| self.reference[*at] != self.alternate[*at])
            .unwrap()
    }

    /// The `Amino_acids` cell under `table`: the reference and alternate
    /// residues around a slash, or the one residue where the two are equal
    /// (`TranscriptVariationAllele.pm:618-630`).
    fn amino_acids(&self, table: u8) -> String {
        let reference = amino_acid(&self.reference, table);
        let alternate = amino_acid(&self.alternate, table);
        if reference == alternate {
            reference.to_string()
        } else {
            format!("{reference}/{alternate}")
        }
    }

    /// The `Codons` cell: both codons in lower case with the changed base in
    /// upper case (`TranscriptVariationAllele.pm:666-681` and `:901-932`).
    fn codons(&self) -> String {
        let display = |codon: &Codon| -> String {
            codon
                .iter()
                .enumerate()
                .map(|(at, base)| {
                    if at == self.offset() {
                        base.to_ascii_uppercase() as char
                    } else {
                        base.to_ascii_lowercase() as char
                    }
                })
                .collect()
        };
        format!("{}/{}", display(&self.reference), display(&self.alternate))
    }
}

/// The substitutions made on both regions of a cache on `first` and `second`.
/// The first leaves the first codon in the module's order that the two codes
/// translate to different amino acids. Where the first codon a code
/// translates unlike the standard code is another codon, a substitution
/// leaving that codon follows, so a region read under the standard code in
/// place of its own writes other amino acids on some row. Each reference
/// codon is the first codon in the module's order one base away from the
/// codon left that neither code nor the standard code reads as a stop codon.
fn substitutions(first: u8, second: u8) -> Vec<Substitution> {
    let mut written = vec![first_codon_read_differently(first, second)
        .expect("the two codes translate some codon differently")];
    for table in [first, second] {
        let unlike_standard =
            first_codon_read_differently(table, STANDARD).filter(|codon| !written.contains(codon));
        written.extend(unlike_standard);
    }
    written
        .into_iter()
        .map(|alternate| {
            let reference = first_codon(|codon| {
                let changed = (0..3).filter(|at| codon[*at] != alternate[*at]).count();
                changed == 1
                    && [first, second, STANDARD]
                        .iter()
                        .all(|table| !is_ter_codon(codon, *table))
            })
            .expect("a sense codon one base away");
            Substitution {
                reference,
                alternate,
            }
        })
        .collect()
}

/// The codons of the coding sequence both regions carry: a codon both codes
/// mark as a start codon, each substitution's reference codon in turn, and a
/// codon both codes read as a stop codon, each end the first such codon in
/// the module's order.
fn coding_codons(first: u8, second: u8, substitutions: &[Substitution]) -> Vec<Codon> {
    let shared = |marked: fn(&[u8], u8) -> bool| {
        first_codon(|codon| marked(codon, first) && marked(codon, second))
            .expect("a codon both codes mark alike")
    };
    let mut codons = vec![shared(is_start_codon)];
    codons.extend(
        substitutions
            .iter()
            .map(|substitution| substitution.reference),
    );
    codons.push(shared(is_ter_codon));
    codons
}

/// The peptide Ensembl's `Transcript::translate` gives `codons` under `table`
/// (ensembl core `Transcript.pm:2152-2165`): each codon's amino acid, the
/// last one dropped where its codon is a stop codon of the table, and `M`
/// first where the first codon is a start codon of the table.
fn peptide(codons: &[Codon], table: u8) -> String {
    let mut residues: Vec<char> = codons
        .iter()
        .map(|codon| amino_acid(codon, table))
        .collect();
    if codons
        .last()
        .is_some_and(|codon| is_ter_codon(codon, table))
    {
        residues.pop();
    }
    if codons
        .first()
        .is_some_and(|codon| is_start_codon(codon, table))
    {
        if let Some(first) = residues.first_mut() {
            *first = 'M';
        }
    }
    residues.into_iter().collect()
}

/// The stable id of the transcript on `region`.
fn transcript_id(region: &str) -> String {
    format!("transcript_of_{region}")
}

/// The name of the input line making substitution `index` on `region`.
fn line_name(region: &str, index: usize) -> String {
    format!("region_{region}_substitution_{index}")
}

/// The cache record of the transcript on `region`: one forward exon holding
/// `codons` and nothing else, naming `table`.
fn record(region: &str, codons: &[Codon], table: u8) -> serde_json::Value {
    let coding: String = codons.iter().flatten().map(|base| *base as char).collect();
    let length = coding.len() as u64;
    let end = TRANSCRIPT_START + length - 1;
    json!({
        "stable_id": transcript_id(region),
        "gene_stable_id": format!("gene_of_{region}"),
        "chr": region,
        "start": TRANSCRIPT_START,
        "end": end,
        "strand": 1,
        "biotype": "protein_coding",
        "source": "ensembl",
        "exons": [{"start": TRANSCRIPT_START, "end": end, "rank": 1,
            "phase": 0, "end_phase": 0}],
        "cdna_coding_start": 1,
        "cdna_coding_end": length,
        "translation": {"stable_id": format!("protein_of_{region}"),
            "start": 1, "end": length},
        "variation_effect_feature_cache": {
            "codon_table": table,
            "translateable_seq": coding,
            "peptide": peptide(codons, table),
            "mapper": {"start_phase": 0, "cdna_coding_start": 1,
                "cdna_coding_end": length,
                "pairs": [{"from_start": 1, "from_end": length,
                    "to_start": TRANSCRIPT_START, "to_end": end, "ori": 1}]},
        },
    })
}

/// The `Feature`, `Amino_acids` and `Codons` cells of one output row.
#[derive(Debug)]
struct Row {
    feature: String,
    amino_acids: String,
    codons: String,
}

/// Writes under `dir` a cache whose two regions carry `codons`, the first on
/// `tables[0]` and the second on `tables[1]`, and one input line per region
/// and substitution, runs the binary, and returns each default-format row
/// under its input line's name.
fn annotate(
    dir: &Path,
    tables: [u8; 2],
    codons: &[Codon],
    substitutions: &[Substitution],
) -> HashMap<String, Row> {
    let cache = dir.join("cache");
    for (region, table) in REGIONS.into_iter().zip(tables) {
        let shard = cache
            .join("transcripts")
            .join(region)
            .join("1-1000000.json");
        fs::create_dir_all(shard.parent().unwrap()).unwrap();
        fs::write(&shard, json!([record(region, codons, table)]).to_string()).unwrap();
    }
    fs::write(
        cache.join("info.json"),
        json!({"assembly": "GRCh38"}).to_string(),
    )
    .unwrap();

    let mut vcf =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for region in REGIONS {
        for (index, substitution) in substitutions.iter().enumerate() {
            let offset = substitution.offset();
            // The start codon precedes the reference codons.
            let position = TRANSCRIPT_START + 3 * (index as u64 + 1) + offset as u64;
            vcf.push_str(&format!(
                "{region}\t{position}\t{}\t{}\t{}\t.\t.\t.\n",
                line_name(region, index),
                substitution.reference[offset] as char,
                substitution.alternate[offset] as char
            ));
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

    let text = fs::read_to_string(&out).unwrap();
    let mut rows = HashMap::new();
    for line in text.lines() {
        let cells: Vec<&str> = line.split('\t').collect();
        let row = Row {
            feature: cells[4].to_string(),
            amino_acids: cells[10].to_string(),
            codons: cells[11].to_string(),
        };
        let earlier = rows.insert(cells[0].to_string(), row);
        assert!(earlier.is_none(), "two rows for {}:\n{text}", cells[0]);
    }
    rows
}

/// One run on a cache of two sequence regions that name two genetic codes,
/// both carrying the same coding sequence, with the same substitutions made
/// on each: every row's `Amino_acids` is the text its own region's table
/// gives, while the `Codons` cells of the two regions are equal. The pairs are
/// the standard code beside each other genetic code the ready-made Ensembl
/// VEP caches name that translates some codon unlike it, and two non-standard
/// codes side by side. The two regions' amino acids differ on the first
/// substitution, and each code has a substitution on which the standard code
/// gives other amino acids, except the standard code itself and table 11,
/// which translates all 64 codons as the standard code does and is run beside
/// table 4.
#[test]
fn rows_of_two_sequences_follow_their_own_genetic_codes() {
    let beside_standard = [2, 3, 4, 5, 6, 9, 12, 13].map(|other| [STANDARD, other]);
    let both_non_standard = [[11, 4], [12, 4]];
    for tables in beside_standard.into_iter().chain(both_non_standard) {
        let [first, second] = tables;
        let pair = format!("tables {first} and {second}");
        let substitutions = substitutions(first, second);
        assert_ne!(
            substitutions[0].amino_acids(first),
            substitutions[0].amino_acids(second),
            "{pair}"
        );
        for table in tables {
            let translates_as_standard = codons_in_table_order()
                .iter()
                .all(|codon| amino_acid(codon, table) == amino_acid(codon, STANDARD));
            assert_eq!(
                translates_as_standard,
                table == STANDARD || table == 11,
                "{pair}: table {table}"
            );
            let told_from_standard = substitutions.iter().any(|substitution| {
                substitution.amino_acids(table) != substitution.amino_acids(STANDARD)
            });
            assert_ne!(
                told_from_standard, translates_as_standard,
                "{pair}: table {table}"
            );
        }

        let dir = tempfile::tempdir().unwrap();
        let codons = coding_codons(first, second, &substitutions);
        let rows = annotate(dir.path(), tables, &codons, &substitutions);
        assert_eq!(rows.len(), 2 * substitutions.len(), "{pair}: {rows:?}");
        for (region, table) in REGIONS.into_iter().zip(tables) {
            for (index, substitution) in substitutions.iter().enumerate() {
                let name = line_name(region, index);
                let row = &rows[&name];
                assert_eq!(row.feature, transcript_id(region), "{pair}: {name}");
                assert_eq!(
                    row.amino_acids,
                    substitution.amino_acids(table),
                    "{pair}: {name}"
                );
                assert_eq!(row.codons, substitution.codons(), "{pair}: {name}");
            }
        }
    }
}

/// A cache whose transcript names genetic code 26 loads, and its rows read
/// `CTG` as alanine, as NCBI's genetic code table gives that code, where the
/// standard code reads leucine. `Bio::Tools::CodonTable` in Ensembl VEP's
/// BioPerl holds no table 26: its `id` warns and falls to id 0
/// (`CodonTable.pm:328-338`), under which `translate` reads the module's last
/// table, the standard translation (`CodonTable.pm:412-460`), so Ensembl VEP
/// writes leucine for the codon, and a dumper on that module stores the
/// peptide the same way, as the record here carries it (`MALLE`).
///
/// The coding sequence is `ATG GCT CTG CTA GAA TAA` from base 1001. The first
/// record leaves `CTA` in place of `CTG` (alanine to leucine), the second
/// leaves `CTG` in place of `CTA` (leucine to alanine), and the third changes
/// a codon both codes read alike. Each expected row is whole.
#[test]
fn rows_of_a_transcript_on_genetic_code_26_read_ctg_as_alanine() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("cache");
    let shard = cache.join("transcripts").join("1").join("1-1000000.json");
    fs::create_dir_all(shard.parent().unwrap()).unwrap();
    let transcript = json!({
        "stable_id": "transcript_on_code_26",
        "gene_stable_id": "gene_on_code_26",
        "chr": "1",
        "start": 1001,
        "end": 1018,
        "strand": 1,
        "biotype": "protein_coding",
        "source": "ensembl",
        "exons": [{"start": 1001, "end": 1018, "rank": 1, "phase": 0, "end_phase": 0}],
        "cdna_coding_start": 1,
        "cdna_coding_end": 18,
        "translation": {"stable_id": "protein_on_code_26", "start": 1, "end": 18},
        "variation_effect_feature_cache": {
            "codon_table": "26",
            "translateable_seq": "ATGGCTCTGCTAGAATAA",
            "peptide": "MALLE",
            "seq_edits": [],
            "mapper": {"start_phase": 0, "cdna_coding_start": 1, "cdna_coding_end": 18,
                "pairs": [{"from_start": 1, "from_end": 18,
                    "to_start": 1001, "to_end": 1018, "ori": 1}]},
        },
    });
    fs::write(&shard, json!([transcript]).to_string()).unwrap();
    fs::write(
        cache.join("info.json"),
        json!({"assembly": "GRCh38"}).to_string(),
    )
    .unwrap();

    let input = dir.path().join("records.vcf");
    fs::write(
        &input,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t1009\tctg_to_cta\tG\tA\t.\t.\t.\n\
         1\t1012\tcta_to_ctg\tA\tG\t.\t.\t.\n\
         1\t1013\tgaa_to_caa\tG\tC\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.path().join("records.out");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(&cache)
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

    let rows = fs::read_to_string(&out).unwrap();
    let rows: Vec<&str> = rows.lines().collect();
    assert_eq!(
        rows,
        [
            "ctg_to_cta\t1:1009\tA\tgene_on_code_26\ttranscript_on_code_26\tTranscript\t\
             missense_variant\t9\t9\t3\tA/L\tctG/ctA\t-\tIMPACT=MODERATE;STRAND=1",
            "cta_to_ctg\t1:1012\tG\tgene_on_code_26\ttranscript_on_code_26\tTranscript\t\
             missense_variant\t12\t12\t4\tL/A\tctA/ctG\t-\tIMPACT=MODERATE;STRAND=1",
            "gaa_to_caa\t1:1013\tC\tgene_on_code_26\ttranscript_on_code_26\tTranscript\t\
             missense_variant\t13\t13\t5\tE/Q\tGaa/Caa\t-\tIMPACT=MODERATE;STRAND=1",
        ]
    );
}
