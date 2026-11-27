// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation of an insertion Ensembl VEP left unrotated on a reverse-strand
//! transcript, on a record's only or first allele and on a later allele read at the first
//! allele's shifted pair, and a later allele's terms read from the peptides at that pair.

mod downstream_common;

use downstream_common::*;
use vep_formats_compare::classes::SHIFTED_ALLELE_ROTATION_BY_STRAND as ROTATION;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::Side;

const CSQ_HEADER: &str = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|STRAND|HGVSc|HGVSp|HGVS_OFFSET\">\n";

/// One allele's CSQ entry: Allele, Consequence, CDS_position, Protein_position,
/// Amino_acids, Codons, HGVSp and HGVS_OFFSET.
type Csq<'a> = (
    &'a str,
    &'a str,
    &'a str,
    &'a str,
    &'a str,
    &'a str,
    &'a str,
    i64,
);

/// The VCF line of `record` on the transcript `feature` of strand `strand`, one entry per
/// allele, each allele's HGVSc its own and the same on both sides.
fn vcf(record: &str, feature: &str, strand: &str, entries: &[Csq]) -> String {
    let csq: Vec<String> = entries
        .iter()
        .map(|(allele, terms, cds, protein, amino_acids, codons, hgvsp, offset)| {
            format!(
                "{allele}|{terms}|Transcript|{feature}|-|{cds}|{protein}|{amino_acids}|{codons}|{strand}|{feature}.1:c.ins{allele}|{hgvsp}|{offset}"
            )
        })
        .collect();
    let line = record.strip_suffix('.').expect("the INFO is a dot");
    format!("{CSQ_HEADER}{line}CSQ={}\n", csq.join(","))
}

/// The report of the two VCF texts of `record` under `classes`, beside `side`.
fn judged(reference: &str, candidate: &str, record: &str, side: &Side, classes: &[&str]) -> Value {
    let inputs = [record.to_string()];
    let mut documented = built_in(classes);
    compared_under(
        &parse_vcf(reference, &inputs).expect("the reference parses"),
        &parse_vcf(candidate, &inputs).expect("the candidate parses"),
        &mut documented,
        side,
        Some(&inputs),
    )
}

/// TranscriptVariationAllele.pm:1337-1339 turn a reverse-strand transcript's allele by its
/// length less the shift, so `GCT` shifted four bases is read unturned where vep-rs turns
/// it once: the HGVSp of `chr6 16327676 T TGCT` on ATXN1 (ENST00000244769), Ensembl VEP
/// 116.3's `p.His226delinsGlnAla` against `p.Gln225dup`. The class reads the rows alone.
/// A forward strand, a shift no longer than the allele or a multiple of its length, another
/// HGVSc, another CDS position and another protein stay differences.
#[test]
fn an_unrotated_reverse_strand_insertion_s_hgvsp_is_the_rotation_class_from_the_rows_alone() {
    let record = "21\t2036\tr\tT\tTGCT\t.\t.\t.";
    let row = |strand: &str, offset: i64, cds: &str, hgvsp: &str| {
        vcf(
            record,
            TOY_REVERSE,
            strand,
            &[(
                "GCT",
                "inframe_insertion",
                cds,
                "2-3",
                "-/Q",
                "-/CAG",
                hgvsp,
                offset,
            )],
        )
    };
    let (vep, vep_rs) = ("ENSP1.1:p.His3delinsGlnAla", "ENSP1.1:p.Gln2dup");
    let side = Side::default();
    let once = count_map(&[(ROTATION, 1)]);
    let accepted = judged(
        &row("-1", -4, "6-7", vep),
        &row("-1", -4, "6-7", vep_rs),
        record,
        &side,
        &[ROTATION],
    );
    assert_eq!(masked(&accepted, "HGVSp"), once);
    assert_eq!(verdict_of(&accepted), "green");
    for (strand, offset) in [("1", -4), ("-1", -3), ("-1", -2), ("-1", -6), ("-1", 0)] {
        let report = judged(
            &row(strand, offset, "6-7", vep),
            &row(strand, offset, "6-7", vep_rs),
            record,
            &side,
            &[ROTATION],
        );
        assert_eq!(differing(&report, "HGVSp"), 1, "{strand} {offset}");
    }
    let other_protein = judged(
        &row("-1", -4, "6-7", vep),
        &row("-1", -4, "6-7", "ENSP2.1:p.Gln2dup"),
        record,
        &side,
        &[ROTATION],
    );
    assert_eq!(differing(&other_protein, "HGVSp"), 1);
    let other_cds = judged(
        &row("-1", -4, "10-11", vep),
        &row("-1", -4, "6-7", vep_rs),
        record,
        &side,
        &[ROTATION],
    );
    assert_eq!(differing(&other_cds, "HGVSp"), 1);
    let other_hgvsc = judged(
        &row("-1", -4, "6-7", vep),
        &row("-1", -4, "6-7", vep_rs).replace("c.insGCT", "c.insTGC"),
        record,
        &side,
        &[ROTATION],
    );
    assert_eq!(differing(&other_hgvsc, "HGVSp"), 1);
}

/// A later allele's HGVSp is written at its own shift (TranscriptVariationAllele.pm:1724-1733),
/// so on a reverse-strand transcript it differs by the rotation alone while its CDS position
/// is the first allele's moved pair: `chr6 16327676 T TGCT,TGCTGCTGCT` on ATXN1 gives the
/// later allele `p.His226delinsGlnAlaAlaAla` at `677-678` against `p.Gln223_Gln225dup` at
/// `634-635`. The shared-shift class takes such an HGVSp; on a forward strand, where both
/// engines turn the allele alike, a differing HGVSp stays a difference.
#[test]
fn hgvs_shift_shared_across_alleles_takes_a_later_allele_s_unrotated_hgvsp() {
    let directory = scratch_of("shared-shift-rotation");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let cases = [
        (
            "21\t2036\tw\tT\tTGG,TGCT\t.\t.\t.",
            TOY_REVERSE,
            "-1",
            -4,
            true,
        ),
        (
            "21\t1006\tw\tA\tAGG,AGCT\t.\t.\t.",
            TOY_FORWARD,
            "1",
            4,
            false,
        ),
    ];
    for (record, feature, strand, offset, taken) in cases {
        let first = (
            "GG",
            "frameshift_variant",
            "6-7",
            "2-3",
            "-/X",
            "-/GG",
            "-",
            offset,
        );
        let text = |cds: &str, protein: &str, hgvsp: &str| {
            vcf(
                record,
                feature,
                strand,
                &[
                    first,
                    (
                        "GCT",
                        "inframe_insertion",
                        cds,
                        protein,
                        "-/X",
                        "-/X",
                        hgvsp,
                        offset,
                    ),
                ],
            )
        };
        let report = judged(
            &text("10-11", "4", "ENSP1.1:p.His3delinsGlnAla"),
            &text("6-7", "2-3", "ENSP1.1:p.Gln2dup"),
            record,
            &side,
            &[SHARED_SHIFT, ROTATION],
        );
        let once = count_map(&[(SHARED_SHIFT, 1)]);
        assert_eq!(masked(&report, "CDS_position"), once, "{record}");
        if taken {
            assert_eq!(masked(&report, "HGVSp"), once, "{record}");
            assert_eq!(verdict_of(&report), "green", "{record}");
        } else {
            assert_eq!(differing(&report, "HGVSp"), 1, "{record}");
        }
    }
}

/// Without statistics a later allele's terms are read with its row, from the peptides at the
/// first allele's moved pair (TranscriptVariationAllele.pm:798-830): an insertion whose
/// alternate peptide neither begins nor ends with the reference peptide there is
/// protein_altering_variant (Utils/VariationEffect.pm:375-393, :1125-1127) where vep-rs's own
/// pair gives inframe_insertion, as on the later allele of `chr6 16327676 T
/// TGCT,TGCTGCTGCT` on ATXN1 (`H/QAAA` against `Q/QQQQ`). The shared-shift class takes such a
/// Consequence beside the Amino_acids and Codons it takes; terms the two peptide pairs do not
/// give, or peptides that are not the moved window, stay differences.
#[test]
fn hgvs_shift_shared_across_alleles_takes_a_later_allele_s_terms_read_at_the_moved_pair() {
    let directory = scratch_of("shared-shift-terms");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let record = "21\t1006\tw\tA\tAAGC,AAGCAGA\t.\t.\t.";
    let first = (
        "AGC",
        "inframe_insertion",
        "6-7",
        "2-3",
        "-/S",
        "-/AGC",
        "-",
        4,
    );
    let text = |terms: &str, cds: &str, protein: &str, amino_acids: &str, codons: &str| {
        vcf(
            record,
            TOY_FORWARD,
            "1",
            &[
                first,
                ("AGCAGA", terms, cds, protein, amino_acids, codons, "-", 4),
            ],
        )
    };
    let candidate = text("inframe_insertion", "6-7", "2-3", "-/SR", "-/AGCAGA");
    let once = count_map(&[(SHARED_SHIFT, 1)]);
    let report = judged(
        &text(
            "protein_altering_variant",
            "10-11",
            "4",
            "L/QQI",
            "-/CAGCAGatt",
        ),
        &candidate,
        record,
        &side,
        &[SHARED_SHIFT],
    );
    for column in ["Consequence", "Amino_acids", "Codons", "CDS_position"] {
        assert_eq!(masked(&report, column), once, "{column}");
    }
    assert_eq!(verdict_of(&report), "green");
    let other_term = judged(
        &text("missense_variant", "10-11", "4", "L/QQI", "-/CAGCAGatt"),
        &candidate,
        record,
        &side,
        &[SHARED_SHIFT],
    );
    assert_eq!(differing(&other_term, "Consequence"), 1);
    let other_window = judged(
        &text(
            "protein_altering_variant",
            "10-11",
            "4",
            "L/QQV",
            "-/CAGCAGatt",
        ),
        &candidate,
        record,
        &side,
        &[SHARED_SHIFT],
    );
    assert_eq!(differing(&other_window, "Consequence"), 1);
}
