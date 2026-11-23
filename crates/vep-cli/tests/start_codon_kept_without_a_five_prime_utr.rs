// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads an insertion or deletion that overlaps the
//! start codon of a transcript with no 5' UTR and leaves a start codon at the
//! start of the coding sequence.
//!
//! Two transcripts of one gene share a start codon at 21:32298972-32298974:
//! `WITH_UTR` holds 117 bases before it and `WITHOUT_UTR` opens on it. Each
//! of the four records of the reproduction changes the length of the coding
//! sequence at that codon and leaves `ATG` as its first three bases, so both
//! rows carry `start_retained_variant` (SO:0002019) and no `start_lost`, with
//! the IMPACT of their terms and, under `--hgvs`, the protein notation of the
//! change itself.
//!
//! Ensembl VEP 116.3 writes `start_retained_variant` for each of the four on
//! the transcript with the 5' UTR. On the other it makes no kept-start test:
//! the test stands inside `if ($utr)` (VariationEffect.pm:1061-1065), a
//! deletion reads as an altered start (:1069) and an insertion by the
//! comparison of :1071. Its row there carries no `start_retained_variant`,
//! and `start_lost` where the allele shifts the frame (:870) or the peptide
//! test holds (:878-882), with the protein notation of a lost start
//! (TranscriptVariationAllele.pm:2161-2164). That is the documented
//! difference `start_codon_kept_without_five_prime_utr`.
//!
//! The records are made up on the chromosome 21 reference the corpus carries
//! and read on its release-116 GRCh38 cache. Each test runs one input line of
//! the reproduction and holds the rows of the two transcripts in the columns
//! the reproduction compares.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The reproduction record of the documented difference.
const REPRODUCTION: &str =
    "scripts/divergences/reproductions/start_codon_kept_without_five_prime_utr.toml";

/// The corpus whose chromosome 21 cache slice and reference the reproduction
/// names.
const CORPUS: &str = "tests/golden/116/GRCh38-plugin-loftee-fixture";

/// A transcript whose coding sequence starts at its cDNA position 118.
const WITH_UTR: &str = "ENST00000303645";
/// A transcript whose coding sequence starts at its first base.
const WITHOUT_UTR: &str = "ENST00000339944";

/// Four bases written between the first and the second base of the start
/// codon. They open with the codon's last two bases, so the codon stands
/// whole before the other two and the frame shifts after it.
const FRAMESHIFT_INSERTION: (&str, [&str; 2]) = (
    "frameshift_insertion_inside_the_codon",
    [
        "frameshift_insertion_inside_the_codon | ENST00000303645 | frameshift_variant,start_retained_variant | 118-119 | 1-2 | HIGH | ENSP00000306697.5:p.Ala2ProfsTer15",
        "frameshift_insertion_inside_the_codon | ENST00000339944 | frameshift_variant,start_retained_variant | 1-2 | 1-2 | HIGH | ENSP00000343661.4:p.Ala2ProfsTer15",
    ],
);

/// Three bases written at the same place. They open with the codon's last two
/// bases, so the codon stands whole before a new leucine codon.
const WHOLE_CODON_INSERTION: (&str, [&str; 2]) = (
    "whole_codon_insertion_inside_the_codon",
    [
        "whole_codon_insertion_inside_the_codon | ENST00000303645 | inframe_insertion,start_retained_variant | 118-119 | 1-2 | MODERATE | ENSP00000306697.5:p.Met1_Ala2insLeu",
        "whole_codon_insertion_inside_the_codon | ENST00000339944 | inframe_insertion,start_retained_variant | 1-2 | 1-2 | MODERATE | ENSP00000343661.4:p.Met1_Ala2insLeu",
    ],
);

/// The start codon's third base removed. The base after it is the same base
/// and takes its place, and the frame shifts after the codon.
const THIRD_BASE_DELETION: (&str, [&str; 2]) = (
    "third_base_removed_before_its_repeat",
    [
        "third_base_removed_before_its_repeat | ENST00000303645 | frameshift_variant,start_retained_variant | 120 | 3 | HIGH | ENSP00000306697.5:p.Ala2ProfsTer30",
        "third_base_removed_before_its_repeat | ENST00000339944 | frameshift_variant,start_retained_variant | 3 | 3 | HIGH | ENSP00000343661.4:p.Ala2ProfsTer30",
    ],
);

/// Nine bases removed from the start codon's third base on. The base after
/// them is the same base as the first of the nine, so the codon stands whole
/// and the protein lacks its residues 2 to 4.
const THREE_CODON_DELETION: (&str, [&str; 2]) = (
    "three_codons_removed_from_the_third_base",
    [
        "three_codons_removed_from_the_third_base | ENST00000303645 | inframe_deletion,start_retained_variant | 120-128 | 3-11 | MODERATE | ENSP00000306697.5:p.Ala2_Gly4del",
        "three_codons_removed_from_the_third_base | ENST00000339944 | inframe_deletion,start_retained_variant | 3-11 | 3-11 | MODERATE | ENSP00000343661.4:p.Ala2_Gly4del",
    ],
);

/// A file or directory of the repository, by its path from the root.
fn committed(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(path)
}

/// The lines of the TOML multi-line string or array of the reproduction that
/// opens on the line `opening`, up to the line that is `closing`.
fn block(opening: &str, closing: &str) -> Vec<String> {
    let path = committed(REPRODUCTION);
    let record = std::fs::read_to_string(&path).unwrap();
    let lines: Vec<String> = record
        .lines()
        .skip_while(|line| *line != opening)
        .skip(1)
        .take_while(|line| *line != closing)
        .map(str::to_string)
        .collect();
    assert!(!lines.is_empty(), "{}: {opening}", path.display());
    lines
}

/// The line of the reproduction's input that carries `id` in its ID column: a
/// VCF line without its header. The record is read here so that each input
/// line is spelt in one file.
fn record_line(id: &str) -> String {
    let lines: Vec<String> = block("text = \"\"\"", "\"\"\"")
        .into_iter()
        .filter(|line| line.split('\t').nth(2) == Some(id))
        .collect();
    assert_eq!(lines.len(), 1, "{REPRODUCTION}: {id}");
    format!("{}\n", lines[0])
}

/// The value of `key` in the `Extra` cell of a row of the default output.
fn extra<'a>(cell: &'a str, key: &str) -> &'a str {
    cell.split(';')
        .find_map(|item| item.strip_prefix(key)?.strip_prefix('='))
        .unwrap_or_else(|| panic!("no {key} in {cell}"))
}

/// The rows `vep --hgvs` writes for the record `id` on the two transcripts,
/// each in the columns the reproduction compares, joined as it joins them:
/// the record's ID, the transcript, the terms, the cDNA and the CDS position,
/// the IMPACT and the HGVSp. The rows are sorted, as the reproduction's are.
fn compared_rows(dir: &Path, id: &str) -> Vec<String> {
    let input = dir.join(format!("{id}.vcf"));
    std::fs::write(&input, format!("{VCF_HEADER}{}", record_line(id))).unwrap();
    let out = dir.join(format!("{id}.txt"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(committed(CORPUS).join("json_cache"))
        .arg("--fasta")
        .arg(committed(CORPUS).join("reference.fa.gz"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--format", "vcf", "--hgvs"])
        .args(["--force_overwrite", "--quiet", "--no_stats"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{id}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut rows: Vec<String> = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|cells| cells[0] == id && [WITH_UTR, WITHOUT_UTR].contains(&cells[4]))
        .map(|cells| {
            [
                cells[0],
                cells[4],
                cells[6],
                cells[7],
                cells[8],
                extra(cells[13], "IMPACT"),
                extra(cells[13], "HGVSp"),
            ]
            .join(" | ")
        })
        .collect();
    rows.sort();
    rows
}

/// A frameshift insertion written inside the start codon that leaves the
/// codon whole: `start_retained_variant` beside `frameshift_variant` on both
/// transcripts, with the frameshift notation. Ensembl VEP 116.3 writes
/// `frameshift_variant,start_lost` and `p.Ala2_?1` on the transcript without
/// a 5' UTR (VariationEffect.pm:1071, :870).
#[test]
fn a_frameshift_insertion_inside_the_start_codon_keeps_the_start_without_a_five_prime_utr() {
    let dir = tempfile::tempdir().unwrap();
    let (id, rows) = FRAMESHIFT_INSERTION;
    assert_eq!(compared_rows(dir.path(), id), rows);
}

/// An insertion of one whole codon written inside the start codon that leaves
/// the codon whole: `start_retained_variant` beside `inframe_insertion` on
/// both transcripts. Ensembl VEP 116.3 writes `inframe_insertion` alone on
/// the transcript without a 5' UTR (VariationEffect.pm:1071).
#[test]
fn a_whole_codon_insertion_inside_the_start_codon_keeps_the_start_without_a_five_prime_utr() {
    let dir = tempfile::tempdir().unwrap();
    let (id, rows) = WHOLE_CODON_INSERTION;
    assert_eq!(compared_rows(dir.path(), id), rows);
}

/// The start codon's third base removed before a base that repeats it:
/// `start_retained_variant` beside `frameshift_variant` on both transcripts,
/// with the frameshift notation. Ensembl VEP 116.3 writes
/// `frameshift_variant,start_lost` and `p.Ala2?` on the transcript without a
/// 5' UTR (VariationEffect.pm:1069, :870).
#[test]
fn the_third_base_removed_before_its_repeat_keeps_the_start_without_a_five_prime_utr() {
    let dir = tempfile::tempdir().unwrap();
    let (id, rows) = THIRD_BASE_DELETION;
    assert_eq!(compared_rows(dir.path(), id), rows);
}

/// Three codons removed from the start codon's third base on, before a base
/// that repeats it: `start_retained_variant` beside `inframe_deletion` on
/// both transcripts, IMPACT `MODERATE`, with a `del` notation. Ensembl VEP
/// 116.3 writes `start_lost,inframe_deletion`, `HIGH` and `p.AlaAsnGly2_?4`
/// on the transcript without a 5' UTR (VariationEffect.pm:1069, :878-882),
/// and both start terms on the other, the documented difference
/// `start_cooccurrence_swap`.
#[test]
fn three_codons_removed_from_the_third_base_keep_the_start_without_a_five_prime_utr() {
    let dir = tempfile::tempdir().unwrap();
    let (id, rows) = THREE_CODON_DELETION;
    assert_eq!(compared_rows(dir.path(), id), rows);
}

/// The reproduction's `vep_rs` list is the eight rows the four tests above
/// hold, in the order of its input.
#[test]
fn the_reproduction_records_the_rows_vep_rs_writes() {
    let recorded: Vec<String> = block("vep_rs = [", "]")
        .into_iter()
        .map(|line| {
            line.trim()
                .trim_end_matches(',')
                .trim_matches('"')
                .to_string()
        })
        .collect();
    let held: Vec<&str> = [
        FRAMESHIFT_INSERTION,
        WHOLE_CODON_INSERTION,
        THIRD_BASE_DELETION,
        THREE_CODON_DELETION,
    ]
    .iter()
    .flat_map(|(_, rows)| rows.iter().copied())
    .collect();
    assert_eq!(recorded, held);
}
