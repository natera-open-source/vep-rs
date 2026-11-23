// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HGVSp the `vep` binary writes for an insertion with unknown bases
//! (`N`) that its record writes in an intron and the HGVS 3' rule places in
//! the coding sequence.
//!
//! An allele written a few bases before an exon can open with the bases that
//! follow it, so the 3' rule describes it past them, inside the exon. Where
//! its known bases spell a stop codon in the reading frame there before its
//! first unknown base, the protein ends at that stop whatever the unknown
//! bases are, and every completion of them has one HGVSp and one
//! HGVS_OFFSET. The row carries both, beside the terms of the intron the
//! record writes the allele in, as the row of an allele of known bases does.
//! The tests below pin those rows and hold them to their completions.
//! Ensembl VEP 116.3 writes no notation for an allele that holds a letter
//! outside A, C, G and T (TranscriptVariationAllele.pm:1397 and 1684), and
//! none for an allele of known bases whose written span is not coding
//! (TranscriptVariationAllele.pm:1731). The coding notation stays empty: the
//! completions share its positions and not its inserted bases.
//! The records are made up on the chromosome 21 reference the corpus carries
//! and read on its release-116 GRCh38 cache.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The identifiers corpus: a chromosome 21 cache slice and its reference.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The four unknown bases of every record, as the record writes them.
const UNKNOWN: &str = "NNNN";

/// A record written in an intron of one transcript, three bases before the
/// next coding exon in the transcript's direction.
struct Record {
    id: &'static str,
    /// The base before the insertion point on the forward strand.
    position: u64,
    shared_base: &'static str,
    /// The inserted bases on the forward strand.
    inserted: &'static str,
    gene: &'static str,
    transcript: &'static str,
    strand: &'static str,
    /// The notation with its protein, and the bases the 3' rule moves the
    /// insertion by, signed by the strand; `-` twice where none is written.
    hgvsp: &'static str,
    hgvs_offset: &'static str,
}

/// 21:10542391 lies in an intron of the forward-strand ENST00000622113; the
/// next exon opens at 21:10542395 with the last base of codon 22, a serine.
/// The inserted bases open with the four reference bases after the point,
/// `TAGT`, so the 3' rule places the insertion after 21:10542395, between
/// codons 22 and 23. There the transcript reads `GAA`, `TAG` and the unknown
/// bases: a glutamate and a stop codon.
const FORWARD_INFRAME: Record = Record {
    id: "forward_inframe",
    position: 10_542_391,
    shared_base: "C",
    inserted: "TAGTGAATAGNNNNC",
    gene: "ENSG00000274391",
    transcript: "ENST00000622113",
    strand: "1",
    hgvsp: "ENSP00000482040.1:p.Ser22_Pro23insGluTer",
    hgvs_offset: "4",
};

/// The same without the last base: fourteen inserted bases, two over whole
/// codons.
const FORWARD_FRAMESHIFT: Record = Record {
    id: "forward_frameshift",
    inserted: "TAGTGAATAGNNNN",
    hgvsp: "ENSP00000482040.1:p.Pro23GlufsTer2",
    ..FORWARD_INFRAME
};

/// The unknown bases follow the glutamate codon, so no stop codon stands
/// before them.
const FORWARD_UNSETTLED: Record = Record {
    id: "forward_unsettled",
    inserted: "TAGTGAANNNNTAGC",
    hgvsp: "-",
    hgvs_offset: "-",
    ..FORWARD_INFRAME
};

/// 21:14166454 and 21:14166455 lie in an intron of the reverse-strand
/// ENST00000681601; the next exon in its direction opens at 21:14166451 with
/// the second base of codon 215, a glycine. On the transcript the inserted
/// bases open with the four bases after the point, so the 3' rule places the
/// insertion between 21:14166450 and 21:14166451, after that base. There the
/// transcript reads `G`, which completes the codon as a glycine, then `GAA`,
/// `TAG` and the unknown bases.
const REVERSE_INFRAME: Record = Record {
    id: "reverse_inframe",
    position: 14_166_454,
    shared_base: "G",
    inserted: "NNNNCTATTCCCCTG",
    gene: "ENSG00000188992",
    transcript: "ENST00000681601",
    strand: "-1",
    hgvsp: "ENSP00000505323.1:p.Gly215_Leu216insGluTer",
    hgvs_offset: "-4",
};

/// The same with one more base after the unknown ones on the transcript:
/// sixteen inserted bases, one over whole codons.
const REVERSE_FRAMESHIFT: Record = Record {
    id: "reverse_frameshift",
    inserted: "GNNNNCTATTCCCCTG",
    hgvsp: "ENSP00000505323.1:p.Leu216GlufsTer2",
    ..REVERSE_INFRAME
};

/// The unknown bases follow the glutamate codon on the transcript.
const REVERSE_UNSETTLED: Record = Record {
    id: "reverse_unsettled",
    inserted: "CTANNNNTTCCCCTG",
    hgvsp: "-",
    hgvs_offset: "-",
    ..REVERSE_INFRAME
};

const SETTLED: [&Record; 4] = [
    &FORWARD_INFRAME,
    &FORWARD_FRAMESHIFT,
    &REVERSE_INFRAME,
    &REVERSE_FRAMESHIFT,
];
const UNSETTLED: [&Record; 2] = [&FORWARD_UNSETTLED, &REVERSE_UNSETTLED];

/// Each unknown base as each base in turn, two mixed fills, and a stop codon
/// repeated from each of its phases on either strand.
const COMPLETIONS: [&str; 12] = [
    "AAAA", "CCCC", "GGGG", "TTTT", "ACGT", "TGCA", "TAAT", "AATA", "ATAA", "TTAT", "TATT", "ATTA",
];

fn vcf_line(record: &Record, id: &str, inserted: &str) -> String {
    format!(
        "21\t{}\t{id}\t{base}\t{base}{inserted}\t.\t.\t.\n",
        record.position,
        base = record.shared_base,
    )
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
const IMPACT: usize = 13;
const HGVSC: usize = 17;
const HGVSP: usize = 18;
const HGVS_OFFSET: usize = 19;

/// The whole row of the record's transcript: the terms and the impact of the
/// intron the record writes the allele in, no position, residue or codon, no
/// HGVSc, and the record's notation and offset.
fn row(record: &Record) -> String {
    [
        record.id,
        &format!("21:{}-{}", record.position, record.position + 1),
        record.inserted,
        record.gene,
        record.transcript,
        "Transcript",
        "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
        "-",
        "-",
        "-",
        "-",
        "-",
        "-",
        "LOW",
        "-",
        record.strand,
        "-",
        "-",
        record.hgvsp,
        record.hgvs_offset,
    ]
    .join("\t")
}

fn rows_of<'a>(rows: &'a [Vec<String>], id: &str) -> Vec<&'a Vec<String>> {
    rows.iter().filter(|row| row[ID] == id).collect()
}

fn written(rows: &[Vec<String>], record: &Record) -> String {
    rows_of(rows, record.id)
        .into_iter()
        .find(|row| row[FEATURE] == record.transcript)
        .unwrap_or_else(|| panic!("{}: no row of {}", record.id, record.transcript))
        .join("\t")
}

#[test]
fn a_stop_codon_settled_at_the_shifted_placement_gives_the_intron_row_its_notation() {
    let dir = tempfile::tempdir().unwrap();
    let records: String = SETTLED
        .iter()
        .map(|record| vcf_line(record, record.id, record.inserted))
        .collect();
    let rows = annotate(dir.path(), "settled", &records);
    for record in SETTLED {
        assert_eq!(written(&rows, record), row(record));
    }
}

#[test]
fn unknown_bases_before_any_stop_codon_leave_the_intron_row_without_a_notation() {
    let dir = tempfile::tempdir().unwrap();
    let mut records = String::new();
    for record in UNSETTLED {
        records.push_str(&vcf_line(record, record.id, record.inserted));
        for (n, completion) in COMPLETIONS.iter().enumerate() {
            records.push_str(&vcf_line(
                record,
                &format!("{}_completed_{n}", record.id),
                &record.inserted.replace(UNKNOWN, completion),
            ));
        }
    }
    let rows = annotate(dir.path(), "unsettled", &records);
    for record in UNSETTLED {
        assert_eq!(written(&rows, record), row(record));
        // The completions share no notation: one that spells a stop codon
        // among the unknown bases ends the protein there.
        let notations: BTreeSet<&str> = (0..COMPLETIONS.len())
            .map(|n| {
                let completed_id = format!("{}_completed_{n}", record.id);
                rows_of(&rows, &completed_id)
                    .into_iter()
                    .find(|row| row[FEATURE] == record.transcript)
                    .unwrap_or_else(|| panic!("{completed_id}: no row of {}", record.transcript))
                    [HGVSP]
                    .as_str()
            })
            .collect();
        assert!(notations.len() > 1, "{}: {notations:?}", record.id);
    }
}

#[test]
fn each_notation_is_the_one_every_completion_of_the_unknown_bases_shares() {
    let dir = tempfile::tempdir().unwrap();
    let mut records = String::new();
    for record in SETTLED {
        records.push_str(&vcf_line(record, record.id, record.inserted));
        for (n, completion) in COMPLETIONS.iter().enumerate() {
            records.push_str(&vcf_line(
                record,
                &format!("{}_completed_{n}", record.id),
                &record.inserted.replace(UNKNOWN, completion),
            ));
        }
    }
    let rows = annotate(dir.path(), "completions", &records);
    let shared = [CONSEQUENCE, IMPACT, HGVSP];
    let mut rows_with_a_notation = 0;
    for record in SETTLED {
        for unknown in rows_of(&rows, record.id) {
            let feature = &unknown[FEATURE];
            let noted = unknown[HGVSP] != "-";
            rows_with_a_notation += usize::from(noted);
            for n in 0..COMPLETIONS.len() {
                let completed_id = format!("{}_completed_{n}", record.id);
                let completed = rows_of(&rows, &completed_id)
                    .into_iter()
                    .find(|row| &row[FEATURE] == feature)
                    .unwrap_or_else(|| panic!("{completed_id}: no row of {feature}"));
                for column in shared {
                    assert_eq!(
                        unknown[column], completed[column],
                        "{completed_id} {feature} column {column}"
                    );
                }
                // The completions share the positions of the coding notation
                // and not its bases: the row writes none, and writes the
                // offset beside its protein notation alone.
                assert_eq!(unknown[HGVSC], "-", "{completed_id} {feature}");
                if noted {
                    assert_ne!(completed[HGVSC], "-", "{completed_id} {feature}");
                    assert_eq!(
                        unknown[HGVS_OFFSET], completed[HGVS_OFFSET],
                        "{completed_id} {feature}"
                    );
                } else {
                    assert_eq!(unknown[HGVS_OFFSET], "-", "{completed_id} {feature}");
                }
            }
        }
    }
    // Seven coding transcripts hold the forward-strand exon and six the
    // reverse-strand one; each is read for two records.
    assert_eq!(rows_with_a_notation, 26);
}
