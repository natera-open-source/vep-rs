// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows the `vep` binary writes under `--hgvs` on the record of the
//! documented difference `start_cooccurrence_swap`: a deletion of the first
//! coding base of a transcript whose 5' UTR ends in the same base, so that
//! the start codon forms again one base earlier and the coding sequence
//! stands whole. The test annotates the input line of the entry's record
//! under `scripts/divergences/reproductions` over the cache and the reference
//! that record names and holds every row of the output. The row of the
//! transcript whose start codon the deletion overlaps carries
//! `start_retained_variant` and the first residue, unchanged, as its HGVSp;
//! Ensembl VEP writes `p.Met1?` there, its `start_lost` holding beside
//! `start_retained_variant` (TranscriptVariationAllele.pm:2161-2164).

use std::path::PathBuf;
use std::process::Stdio;

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

const ENTRY: &str = "start_cooccurrence_swap";

/// A file or directory of the repository, by its path from the root.
fn committed(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(path)
}

/// The line of the entry's reproduction record that carries `id` in its ID
/// column: a line of the `text` of the record's first `[input]` table, a VCF
/// without its header. The record is read here so that the input line is
/// spelt in one file.
fn record_line(id: &str) -> String {
    let path = committed(&format!("scripts/divergences/reproductions/{ENTRY}.toml"));
    let record = std::fs::read_to_string(&path).unwrap();
    let lines: Vec<&str> = record
        .lines()
        .skip_while(|line| *line != "[input]")
        .skip_while(|line| *line != "text = \"\"\"")
        .skip(1)
        .take_while(|line| *line != "\"\"\"")
        .filter(|line| line.split('\t').nth(2) == Some(id))
        .collect();
    assert_eq!(lines.len(), 1, "{}: {id}", path.display());
    format!("{}\n", lines[0])
}

/// Every row of the record's output under `--hgvs`. On the transcript whose
/// start codon the deletion overlaps, the coding sequence stands whole one
/// base earlier, and the notation is the first residue, unchanged; the
/// default output writes the `=` of a notation as `%3D`. On the four
/// transcripts whose coding sequence opens 72 bases upstream the same base
/// is coding base 73, and the deletion is the frameshift those rows name.
#[test]
fn a_deletion_that_leaves_the_coding_sequence_whole_behind_the_start_changes_no_residue() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("record.vcf");
    std::fs::write(&input, format!("{VCF_HEADER}{}", record_line("4629713"))).unwrap();
    let out = dir.path().join("out.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(committed("tests/golden/116/GRCh37/json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh37"])
        .args(["--format", "vcf"])
        .arg("--fasta")
        .arg(committed("tests/golden/116/GRCh37-hgvs/reference.fa.gz"))
        .args(["--force_overwrite", "--quiet", "--no_stats", "--hgvs"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let rows: Vec<String> = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect();
    let row = |feature: &str, cells: &str, extra: &str| {
        format!(
            "4629713\t21:36265246\t-\tENSG00000159216\t{feature}\tTranscript\t{cells}\t-\t{extra}"
        )
    };
    let upstream = |feature: &str, distance: u32, flags: &str| {
        row(
            feature,
            "upstream_gene_variant\t-\t-\t-\t-\t-",
            &format!("IMPACT=MODIFIER;DISTANCE={distance};STRAND=-1{flags}"),
        )
    };
    let frameshift_at_73 = |feature: &str, cdna: u32, flags: &str, hgvs: &str| {
        row(
            feature,
            &format!("frameshift_variant\t{cdna}\t73\t25\tM/X\tAtg/tg"),
            &format!("IMPACT=HIGH;STRAND=-1{flags};{hgvs}"),
        )
    };
    let expected = [
        frameshift_at_73(
            "ENST00000300305",
            518,
            "",
            "HGVSc=ENST00000300305.3:c.73del;HGVSp=ENSP00000300305.3:p.Met25Ter",
        ),
        upstream("ENST00000325074", 3212, ""),
        upstream("ENST00000344691", 4259, ""),
        upstream("ENST00000358356", 4259, ""),
        upstream("ENST00000399237", 3159, ";FLAGS=cds_end_NF"),
        frameshift_at_73(
            "ENST00000416754",
            322,
            ";FLAGS=cds_end_NF",
            "HGVSc=ENST00000416754.1:c.73del;HGVSp=ENSP00000405158.1:p.Met25Ter",
        ),
        frameshift_at_73(
            "ENST00000437180",
            263,
            "",
            "HGVSc=ENST00000437180.1:c.73del;HGVSp=ENSP00000409227.1:p.Met25Ter",
        ),
        row(
            "ENST00000455571",
            "intron_variant\t-\t-\t-\t-\t-",
            "IMPACT=MODIFIER;STRAND=-1;FLAGS=cds_end_NF;HGVSc=ENST00000455571.1:c.59-5853del",
        ),
        frameshift_at_73(
            "ENST00000475045",
            1507,
            ";FLAGS=cds_end_NF",
            "HGVSc=ENST00000475045.2:c.73del;HGVSp=ENSP00000477072.1:p.Met25Ter",
        ),
        row(
            "ENST00000482318",
            "intron_variant,NMD_transcript_variant\t-\t-\t-\t-\t-",
            "IMPACT=MODIFIER;STRAND=-1;HGVSc=ENST00000482318.1:c.59-12236del",
        ),
        row(
            "ENST00000486278",
            "frameshift_variant,start_retained_variant\t217\t1\t1\tM/X\tAtg/tg",
            "IMPACT=HIGH;STRAND=-1;HGVSc=ENST00000486278.2:c.1del;\
             HGVSp=ENSP00000438019.1:p.Met1%3D",
        ),
    ];
    assert_eq!(rows, expected);
}
