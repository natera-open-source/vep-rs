// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads an alternate allele that holds unknown bases
//! (`N`) after bases that settle a stop codon.
//!
//! A record can spell an inserted sequence known at its ends and unknown in
//! between. Where the whole codons before the first codon that holds an
//! unknown base hold a stop codon, the protein ends at that stop whatever the
//! unknown bases are, so the row has one reading: the residues through the
//! stop, `stop_gained` beside the term the length gives, the IMPACT of those
//! terms and, under `--hgvs`, the protein notation at the allele's last
//! placement with its HGVS_OFFSET. That is the row every completion of the
//! unknown bases shares, and the tests below hold each row to its completions.
//! Ensembl VEP 116.3 translates no allele that holds a letter outside A, C, G
//! and T (TranscriptVariationAllele.pm:701 and 1683-1684), so it writes no
//! residues, no `stop_gained` and no HGVSp for such an allele.
//! The coding notation stays empty: the completions share its positions and
//! not its inserted bases. An allele whose known bases settle no stop codon
//! before its first unknown base keeps the row of an unknown change.
//! The records are made up on the chromosome 21 reference the corpus carries
//! and read on its release-116 GRCh38 cache.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The identifiers corpus: a chromosome 21 cache slice and its reference.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The base before the insertion point, 21:10567700, the last base of codon
/// 141 of [`TRANSCRIPT`]. The codons after it are those of Ile142 and Leu143.
const POSITION: u64 = 10_567_700;
const SHARED_BASE: &str = "C";

/// A forward-strand coding transcript that holds the insertion point inside
/// one exon, between `c.423` and `c.424`.
const TRANSCRIPT: &str = "ENST00000622113";
const GENE: &str = "ENSG00000274391";
const PROTEIN: &str = "ENSP00000482040.1";

/// Four unknown bases.
const UNKNOWN: &str = "NNNN";

/// A sense codon and a stop codon, then the unknown bases: ten inserted
/// bases, one over whole codons.
const FRAMESHIFT: (&str, &str, &str) = ("stop_then_unknown_frameshift", "GAATAG", "");
/// The same with two known bases after the unknown ones: twelve inserted
/// bases, whole codons.
const INFRAME: (&str, &str, &str) = ("stop_then_unknown_inframe", "GAATAG", "CA");
/// An insertion whose first four bases equal the four reference bases after
/// the insertion point, so its last placement is four bases 3', inside the
/// codon of Leu143, where the first inserted bases complete a stop codon.
const SHIFTED: (&str, &str, &str) = ("stop_then_unknown_shifted", "ATTTAG", "CA");
/// Two sense codons before the unknown bases: no stop codon is settled.
const NO_STOP: (&str, &str, &str) = ("unknown_before_any_stop", "GAAGCT", "CA");

/// The inserted bases of a record, its unknown bases spelled as `middle`.
fn inserted(record: (&str, &str, &str), middle: &str) -> String {
    format!("{}{middle}{}", record.1, record.2)
}

fn vcf_line(id: &str, inserted: &str) -> String {
    format!("21\t{POSITION}\t{id}\t{SHARED_BASE}\t{SHARED_BASE}{inserted}\t.\t.\t.\n")
}

/// The data rows `vep --tab --hgvs` writes for `records` on the corpus.
fn annotate(dir: &Path, name: &str, records: &str) -> Vec<Vec<String>> {
    let input = dir.join(format!("{name}.vcf"));
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join(format!("{name}.tab"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus().join("json_cache"))
        .arg("--fasta")
        .arg(corpus().join("reference.fa.gz"))
        .args([
            "--tab",
            "--hgvs",
            "--no_stats",
            "--quiet",
            "--force_overwrite",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

/// Columns of the tab output under `--tab --hgvs`.
const ID: usize = 0;
const FEATURE: usize = 4;
const CONSEQUENCE: usize = 6;
const CDNA_POSITION: usize = 7;
const CDS_POSITION: usize = 8;
const PROTEIN_POSITION: usize = 9;
const AMINO_ACIDS: usize = 10;
const CODONS: usize = 11;
const IMPACT: usize = 13;
const HGVSC: usize = 17;
const HGVSP: usize = 18;
const HGVS_OFFSET: usize = 19;

/// The whole row of [`TRANSCRIPT`] for a record.
fn row(
    record: (&str, &str, &str),
    consequence: &str,
    amino_acids: &str,
    impact: &str,
    hgvsp: &str,
    hgvs_offset: &str,
) -> String {
    let allele = inserted(record, UNKNOWN);
    [
        record.0,
        &format!("21:{POSITION}-{}", POSITION + 1),
        &allele,
        GENE,
        TRANSCRIPT,
        "Transcript",
        consequence,
        "775-776",
        "423-424",
        "141-142",
        amino_acids,
        &format!("-/{allele}"),
        "-",
        impact,
        "-",
        "1",
        "-",
        "-",
        hgvsp,
        hgvs_offset,
    ]
    .join("\t")
}

fn rows_of<'a>(rows: &'a [Vec<String>], id: &str) -> Vec<&'a Vec<String>> {
    rows.iter().filter(|row| row[ID] == id).collect()
}

#[test]
fn a_stop_codon_settled_before_the_unknown_bases_is_a_gained_stop() {
    let dir = tempfile::tempdir().unwrap();
    let records: String = [FRAMESHIFT, INFRAME, SHIFTED]
        .iter()
        .map(|record| vcf_line(record.0, &inserted(*record, UNKNOWN)))
        .collect();
    let rows = annotate(dir.path(), "settled", &records);
    let written = |record: (&str, &str, &str)| -> String {
        rows_of(&rows, record.0)
            .into_iter()
            .find(|row| row[FEATURE] == TRANSCRIPT)
            .unwrap_or_else(|| panic!("{}: no row of {TRANSCRIPT}", record.0))
            .join("\t")
    };
    // A codon that holds an unknown base reads `X`; the partial codon that
    // ends a frameshift reads `X` as for any allele.
    assert_eq!(
        written(FRAMESHIFT),
        row(
            FRAMESHIFT,
            "stop_gained,frameshift_variant",
            "-/E*XX",
            "HIGH",
            &format!("{PROTEIN}:p.Ile142GlufsTer2"),
            "-",
        )
    );
    assert_eq!(
        written(INFRAME),
        row(
            INFRAME,
            "stop_gained,inframe_insertion",
            "-/E*XX",
            "HIGH",
            &format!("{PROTEIN}:p.Asn141_Ile142insGluTer"),
            "-",
        )
    );
    assert_eq!(
        written(SHIFTED),
        row(
            SHIFTED,
            "stop_gained,inframe_insertion",
            "-/I*XX",
            "HIGH",
            &format!("{PROTEIN}:p.Ile142_Leu143insTer"),
            "4",
        )
    );
}

#[test]
fn unknown_bases_before_any_stop_codon_keep_the_row_of_an_unknown_change() {
    let dir = tempfile::tempdir().unwrap();
    let rows = annotate(
        dir.path(),
        "unsettled",
        &vcf_line(NO_STOP.0, &inserted(NO_STOP, UNKNOWN)),
    );
    let written = rows_of(&rows, NO_STOP.0)
        .into_iter()
        .find(|row| row[FEATURE] == TRANSCRIPT)
        .unwrap()
        .join("\t");
    assert_eq!(
        written,
        row(
            NO_STOP,
            "coding_sequence_variant",
            "-",
            "MODIFIER",
            "-",
            "-"
        )
    );
}

/// The alternate half of `Amino_acids` through its first stop.
fn through_the_stop(amino_acids: &str) -> Option<&str> {
    let alternate = amino_acids.split_once('/')?.1;
    let stop = alternate.find('*')?;
    Some(&alternate[..=stop])
}

#[test]
fn each_row_is_the_row_every_completion_of_the_unknown_bases_shares() {
    let dir = tempfile::tempdir().unwrap();
    // Each unknown base as each base in turn, then completions whose first
    // unknown codon is a stop codon and completions of mixed bases.
    let completions = [
        "AAAA", "CCCC", "GGGG", "TTTT", "TAGT", "TGAC", "GATC", "CTAA",
    ];
    let mut records = String::new();
    for record in [FRAMESHIFT, INFRAME, SHIFTED] {
        records.push_str(&vcf_line(record.0, &inserted(record, UNKNOWN)));
        for (n, completion) in completions.iter().enumerate() {
            records.push_str(&vcf_line(
                &format!("{}_completed_{n}", record.0),
                &inserted(record, completion),
            ));
        }
    }
    let rows = annotate(dir.path(), "completions", &records);
    let shared = [
        CONSEQUENCE,
        CDNA_POSITION,
        CDS_POSITION,
        PROTEIN_POSITION,
        IMPACT,
        HGVSP,
        HGVS_OFFSET,
    ];
    let mut coding_rows = 0;
    for record in [FRAMESHIFT, INFRAME, SHIFTED] {
        let of_the_record = rows_of(&rows, record.0);
        assert_eq!(of_the_record.len(), 9, "{}", record.0);
        for unknown in of_the_record {
            let feature = &unknown[FEATURE];
            let coding = unknown[CDS_POSITION] != "-";
            coding_rows += usize::from(coding);
            for n in 0..completions.len() {
                let id = format!("{}_completed_{n}", record.0);
                let completed = rows_of(&rows, &id)
                    .into_iter()
                    .find(|row| &row[FEATURE] == feature)
                    .unwrap_or_else(|| panic!("{id}: no row of {feature}"));
                for column in shared {
                    assert_eq!(
                        unknown[column], completed[column],
                        "{id} {feature} column {column}"
                    );
                }
                if coding {
                    assert_eq!(
                        through_the_stop(&unknown[AMINO_ACIDS]),
                        through_the_stop(&completed[AMINO_ACIDS]),
                        "{id} {feature}"
                    );
                    assert!(through_the_stop(&unknown[AMINO_ACIDS]).is_some());
                    // The completions share the positions of the coding
                    // notation and not its bases: the row writes none.
                    assert_eq!(unknown[HGVSC], "-", "{id} {feature}");
                    assert_ne!(completed[HGVSC], "-", "{id} {feature}");
                    assert_ne!(unknown[CODONS], completed[CODONS], "{id} {feature}");
                }
            }
        }
    }
    // Eight coding transcripts hold the insertion point.
    assert_eq!(coding_rows, 24);
}
