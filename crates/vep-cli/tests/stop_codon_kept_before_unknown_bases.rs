// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads an insertion whose known bases spell a stop
//! codon at the transcript's old stop position and whose later bases are
//! unknown (`N`).
//!
//! An insertion before or inside the stop codon can begin with bases that
//! leave a stop codon where the reference's stood, behind residues that are
//! all the reference's. The protein is then the reference protein whatever
//! the unknown bases are, and every completion of them has one row:
//! `inframe_insertion` beside `stop_retained_variant`, the residues through
//! the stop and, under `--hgvs`, the protein notation of a kept stop codon
//! at the allele's last placement with its HGVS_OFFSET. The tests below pin
//! that row and hold it to its completions.
//! Ensembl VEP 116.3 translates no allele that holds a letter outside A, C, G
//! and T (TranscriptVariationAllele.pm:701 and 1683-1684), and
//! `inframe_insertion` needs both peptides (Utils/VariationEffect.pm:1113-1115),
//! so it writes no residues, no `inframe_insertion` and no HGVSp for such an
//! allele. The coding notation stays empty: the completions share its
//! positions and not its inserted bases.
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

/// Where a record inserts its bases: between the last sense codon and the
/// stop codon of one transcript.
struct Site {
    /// The base before the insertion point on the forward strand.
    position: u64,
    shared_base: &'static str,
    gene: &'static str,
    transcript: &'static str,
    protein: &'static str,
    strand: &'static str,
    cdna_position: &'static str,
    cds_position: &'static str,
    protein_position: &'static str,
    /// The stop codon's place in the protein.
    stop: u32,
    /// The bases the 3' rule moves the insertion by, signed by the strand.
    hgvs_offset: &'static str,
    /// Whether the record writes the reverse complement of the bases the
    /// transcript reads.
    reverse: bool,
}

/// 21:10605549 is the last base of codon 533, an aspartate, of the
/// forward-strand ENST00000622113; its stop codon `TAA` is 21:10605550-10605552.
const FORWARD: Site = Site {
    position: 10_605_549,
    shared_base: "T",
    gene: "ENSG00000274391",
    transcript: "ENST00000622113",
    protein: "ENSP00000482040.1",
    strand: "1",
    cdna_position: "1951-1952",
    cds_position: "1599-1600",
    protein_position: "533-534",
    stop: 534,
    hgvs_offset: "2",
    reverse: false,
};

/// 21:14108995 is the first base of the stop codon `TAA` of the
/// reverse-strand ENST00000681601, which runs down to 21:14108993; the base
/// after it on the forward strand, 21:14108996, is the last base of codon
/// 460, a threonine.
const REVERSE: Site = Site {
    position: 14_108_995,
    shared_base: "A",
    gene: "ENSG00000188992",
    transcript: "ENST00000681601",
    protein: "ENSP00000505323.1",
    strand: "-1",
    cdna_position: "1490-1491",
    cds_position: "1380-1381",
    protein_position: "460-461",
    stop: 461,
    hgvs_offset: "-2",
    reverse: true,
};

/// Four unknown bases.
const UNKNOWN: &str = "NNNN";

/// A stop codon, then the unknown bases: seven inserted bases, one over
/// whole codons. The fields are the bases before the unknown ones and after
/// them, as the transcript reads them.
const FRAMESHIFT: (&str, &str) = ("TAG", "");
/// The same with two known bases after the unknown ones: nine inserted
/// bases, whole codons.
const INFRAME: (&str, &str) = ("TAG", "CA");
/// The stop codon's own third base is unknown: no stop codon is spelt at
/// the old stop position.
const UNKNOWN_STOP: (&str, &str) = ("TA", "GCA");

/// The inserted bases as the transcript reads them, the unknown bases spelt
/// as `middle`.
fn inserted(record: (&str, &str), middle: &str) -> String {
    format!("{}{middle}{}", record.0, record.1)
}

fn reverse_complement(bases: &str) -> String {
    bases
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// The inserted bases as the record writes them.
fn written_bases(site: &Site, bases: &str) -> String {
    if site.reverse {
        reverse_complement(bases)
    } else {
        bases.to_owned()
    }
}

fn vcf_line(site: &Site, id: &str, bases: &str) -> String {
    format!(
        "21\t{}\t{id}\t{base}\t{base}{}\t.\t.\t.\n",
        site.position,
        written_bases(site, bases),
        base = site.shared_base,
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
const CDNA_POSITION: usize = 7;
const CDS_POSITION: usize = 8;
const PROTEIN_POSITION: usize = 9;
const AMINO_ACIDS: usize = 10;
const CODONS: usize = 11;
const IMPACT: usize = 13;
const HGVSC: usize = 17;
const HGVSP: usize = 18;
const HGVS_OFFSET: usize = 19;

/// The whole row of the site's transcript for a record of unknown bases.
/// `cells` are its Consequence, Amino_acids, IMPACT, HGVSp and HGVS_OFFSET.
fn row(site: &Site, id: &str, record: (&str, &str), cells: [&str; 5]) -> String {
    let [consequence, amino_acids, impact, hgvsp, hgvs_offset] = cells;
    let read = inserted(record, UNKNOWN);
    [
        id,
        &format!("21:{}-{}", site.position, site.position + 1),
        &written_bases(site, &read),
        site.gene,
        site.transcript,
        "Transcript",
        consequence,
        site.cdna_position,
        site.cds_position,
        site.protein_position,
        amino_acids,
        &format!("-/{read}"),
        "-",
        impact,
        "-",
        site.strand,
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

fn written(rows: &[Vec<String>], site: &Site, id: &str) -> String {
    rows_of(rows, id)
        .into_iter()
        .find(|row| row[FEATURE] == site.transcript)
        .unwrap_or_else(|| panic!("{id}: no row of {}", site.transcript))
        .join("\t")
}

#[test]
fn a_stop_codon_kept_before_the_unknown_bases_reads_its_insertion() {
    let dir = tempfile::tempdir().unwrap();
    let cases = [
        (&FORWARD, "forward_frameshift", FRAMESHIFT),
        (&FORWARD, "forward_inframe", INFRAME),
        (&REVERSE, "reverse_frameshift", FRAMESHIFT),
        (&REVERSE, "reverse_inframe", INFRAME),
    ];
    let records: String = cases
        .iter()
        .map(|(site, id, record)| vcf_line(site, id, &inserted(*record, UNKNOWN)))
        .collect();
    let rows = annotate(dir.path(), "kept", &records);
    for (site, id, record) in cases {
        // The stop codon reads `*`, a codon that holds an unknown base `X`,
        // and the base left after two whole codons `X` as on any such row.
        assert_eq!(
            written(&rows, site, id),
            row(
                site,
                id,
                record,
                [
                    "inframe_insertion,stop_retained_variant",
                    "-/*XX",
                    "MODERATE",
                    &format!("{}:p.Ter{}%3D", site.protein, site.stop),
                    site.hgvs_offset,
                ],
            )
        );
    }
}

#[test]
fn an_unknown_base_in_the_stop_codon_keeps_the_row_of_an_unknown_change() {
    let dir = tempfile::tempdir().unwrap();
    let cases = [
        (&FORWARD, "forward_unknown_stop"),
        (&REVERSE, "reverse_unknown_stop"),
    ];
    let records: String = cases
        .iter()
        .map(|(site, id)| vcf_line(site, id, &inserted(UNKNOWN_STOP, UNKNOWN)))
        .collect();
    let rows = annotate(dir.path(), "unsettled", &records);
    for (site, id) in cases {
        assert_eq!(
            written(&rows, site, id),
            row(
                site,
                id,
                UNKNOWN_STOP,
                ["coding_sequence_variant", "-", "MODIFIER", "-", "-"],
            )
        );
    }
}

/// The alternate half of `Amino_acids` through its first stop.
fn through_the_stop(amino_acids: &str) -> Option<&str> {
    let alternate = amino_acids.split_once('/')?.1;
    let stop = alternate.find('*')?;
    Some(&alternate[..=stop])
}

#[test]
fn each_kept_stop_row_is_the_row_every_completion_of_the_unknown_bases_shares() {
    let dir = tempfile::tempdir().unwrap();
    // Each unknown base as each base in turn, then completions whose first
    // unknown codon is a stop codon and completions of mixed bases.
    let completions = [
        "AAAA", "CCCC", "GGGG", "TTTT", "TAGT", "TGAC", "GATC", "CTAA",
    ];
    let cases = [
        (&FORWARD, "forward_frameshift", FRAMESHIFT),
        (&FORWARD, "forward_inframe", INFRAME),
        (&REVERSE, "reverse_frameshift", FRAMESHIFT),
        (&REVERSE, "reverse_inframe", INFRAME),
    ];
    let mut records = String::new();
    for (site, id, record) in cases {
        records.push_str(&vcf_line(site, id, &inserted(record, UNKNOWN)));
        for (n, completion) in completions.iter().enumerate() {
            records.push_str(&vcf_line(
                site,
                &format!("{id}_completed_{n}"),
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
    let mut kept_stop_rows = 0;
    for (_, id, _) in cases {
        for unknown in rows_of(&rows, id) {
            let feature = &unknown[FEATURE];
            let kept = unknown[CONSEQUENCE].contains("stop_retained_variant");
            kept_stop_rows += usize::from(kept);
            for n in 0..completions.len() {
                let completed_id = format!("{id}_completed_{n}");
                let completed = rows_of(&rows, &completed_id)
                    .into_iter()
                    .find(|row| &row[FEATURE] == feature)
                    .unwrap_or_else(|| panic!("{completed_id}: no row of {feature}"));
                // A transcript that holds the insertion point elsewhere
                // has no stop codon there to keep.
                if !kept {
                    continue;
                }
                for column in shared {
                    assert_eq!(
                        unknown[column], completed[column],
                        "{completed_id} {feature} column {column}"
                    );
                }
                assert_eq!(
                    through_the_stop(&unknown[AMINO_ACIDS]),
                    through_the_stop(&completed[AMINO_ACIDS]),
                    "{completed_id} {feature}"
                );
                assert!(through_the_stop(&unknown[AMINO_ACIDS]).is_some());
                // The completions share the positions of the coding
                // notation and not its bases: the row writes none.
                assert_eq!(unknown[HGVSC], "-", "{completed_id} {feature}");
                assert_ne!(completed[HGVSC], "-", "{completed_id} {feature}");
                assert_ne!(unknown[CODONS], completed[CODONS]);
            }
        }
    }
    // Six coding transcripts end on the forward-strand stop codon and five on
    // the reverse-strand one; each is read for two records.
    assert_eq!(kept_stop_rows, 22);
}
