// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvs_shifted_insertion_clamped_end` on a model that carries an edit, where Ensembl VEP
//! writes its clamped notation's alleles over from the row: a substitution printed with `-`
//! as its reference base, and an `ins` of the row's whole allele, each against the
//! candidate's duplication. The class is asked as a comparison asks it, and two outputs are
//! compared in the default, tab, VCF and JSON formats.

mod compare_common;
mod shifted_edited_common;

use std::collections::BTreeMap;

use compare_common::*;
use shifted_edited_common::*;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::pyjson::Value;

const CLASS: &str = "hgvs_shifted_insertion_clamped_end";
const INPUT: &str = "21\t26393634\tfirst_base_copy\tC\tCA\t.\t.\t.";

/// The cells of an insertion's row that the class reads beside its two notations.
#[derive(Clone, Copy)]
struct Row<'a> {
    location: &'a str,
    allele: &'a str,
    feature: &'a str,
    cdna: &'a str,
    strand: &'a str,
    /// The HGVS_OFFSET, and the SHIFT_LENGTH where the row prints one.
    offset: &'a str,
    shift: &'a str,
    /// The BAM_EDIT, empty for a row that prints none.
    bam_edit: &'a str,
}

/// One base inserted before the first base of a forward-strand model, a copy of that base:
/// the shift carries it one base into the transcript.
const ONE_BASE: Row = Row {
    location: "21:26393634-26393635",
    allele: "A",
    feature: "NR_135516.1",
    cdna: "1-2",
    strand: "1",
    offset: "1",
    shift: "1",
    bam_edit: "OK",
};

/// The same on a reverse-strand model.
const ONE_BASE_REVERSE: Row = Row {
    location: "21:18759812-18759813",
    allele: "G",
    feature: "NR_109925.1",
    cdna: "1-2",
    strand: "-1",
    offset: "-1",
    shift: "1",
    bam_edit: "OK",
};

/// Two bases inserted before the first base of the reverse-strand model, carried three bases
/// into the transcript.
const TWO_BASES_REVERSE: Row = Row {
    location: "21:18759812-18759813",
    allele: "CG",
    feature: "NR_109925.1",
    cdna: "3-4",
    strand: "-1",
    offset: "-3",
    shift: "3",
    bam_edit: "OK",
};

impl Row<'_> {
    /// The row's line of a default-format output with the notation after `feature:`.
    fn line(&self, notation: &str) -> String {
        let mut extra = format!("IMPACT=MODIFIER;STRAND={}", self.strand);
        if !self.bam_edit.is_empty() {
            extra.push_str(&format!(";BAM_EDIT={}", self.bam_edit));
        }
        extra.push_str(&format!(
            ";HGVSc={}:{notation};HGVS_OFFSET={}",
            self.feature, self.offset
        ));
        if !self.shift.is_empty() {
            extra.push_str(&format!(";SHIFT_LENGTH={}", self.shift));
        }
        format!(
            "first_base_copy\t{}\t{}\tENSG1\t{}\tTranscript\tnon_coding_transcript_exon_variant\t{}\t-\t-\t-\t-\t-\t{extra}",
            self.location, self.allele, self.feature, self.cdna
        )
    }

    /// What the class says of the row's HGVSc, the reference's notation `reference` and
    /// the candidate's `candidate`.
    fn reads(&self, reference: &str, candidate: &str) -> bool {
        let (e, a) = (
            default_entry(&self.line(reference)),
            default_entry(&self.line(candidate)),
        );
        let value = |notation: &str| format!("{}:{notation}", self.feature);
        class_accepts(
            CLASS,
            &e,
            "HGVSc",
            &value(reference),
            &value(candidate),
            Some(&a),
        )
    }

    /// The verdict of the row's two default-format outputs under the class and what was
    /// applied.
    fn compared(&self, reference: &str, candidate: &str) -> (String, BTreeMap<String, i64>) {
        default_under(
            CLASS,
            &default_text(&[self.line(reference)]),
            &default_text(&[self.line(candidate)]),
        )
    }
}

/// The row of `ONE_BASE` in the tab, VCF and JSON formats with the notation `notation`.
fn in_format(fmt: &str, notation: &str) -> Vec<Entry> {
    let inputs = vec![INPUT.to_string()];
    let hgvsc = format!("{}:{notation}", ONE_BASE.feature);
    match fmt {
        "tab" => {
            let header = tab_header(&[
                "IMPACT",
                "STRAND",
                "BAM_EDIT",
                "HGVSc",
                "HGVS_OFFSET",
                "SHIFT_LENGTH",
            ]);
            let row = format!(
                "first_base_copy\t{}\t{}\tENSG1\t{}\tTranscript\tnon_coding_transcript_exon_variant\t{}\t-\t-\t-\t-\t-\tMODIFIER\t1\tOK\t{hgvsc}\t1\t1",
                ONE_BASE.location, ONE_BASE.allele, ONE_BASE.feature, ONE_BASE.cdna
            );
            tab_entries(&format!("{header}\n{row}\n"))
        }
        "vcf" => {
            let header = csq_header(
                "Allele|Consequence|Feature|Feature_type|cDNA_position|BAM_EDIT|HGVSc|HGVS_OFFSET|SHIFT_LENGTH",
            );
            let csq = [
                ONE_BASE.allele,
                "non_coding_transcript_exon_variant",
                ONE_BASE.feature,
                "Transcript",
                ONE_BASE.cdna,
                "OK",
                &hgvsc,
                "1",
                "1",
            ]
            .join("|");
            vcf_entries(
                &format!("{header}{}CSQ={csq}\n", record_columns(INPUT)),
                &inputs,
            )
        }
        _ => {
            let record = [
                ("input", text(INPUT)),
                ("id", text("first_base_copy")),
                ("allele_string", text("-/A")),
                (
                    "most_severe_consequence",
                    text("non_coding_transcript_exon_variant"),
                ),
            ];
            let consequence = object(&[
                ("variant_allele", text(ONE_BASE.allele)),
                ("transcript_id", text(ONE_BASE.feature)),
                (
                    "consequence_terms",
                    Value::List(vec![text("non_coding_transcript_exon_variant")]),
                ),
                ("cdna_start", Value::Int(1)),
                ("cdna_end", Value::Int(2)),
                ("bam_edit", text("OK")),
                ("hgvsc", text(&hgvsc)),
                ("hgvs_offset", Value::Int(1)),
                ("shift_length", Value::Int(1)),
            ]);
            json_entries(&json_line(&record, consequence), &inputs)
        }
    }
}

fn red() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), nothing())
}

/// The reference's substitution at P with `-` as its reference base against the candidate's
/// duplication of the base before P is the class's, on either strand and in the 5' UTR's
/// numbering, and on a row that prints no BAM_EDIT, as a mitochondrial model's does; so is,
/// on a row carrying `BAM_EDIT` `OK`, the reference's `ins` of the row's whole allele between
/// P-1 and P against the candidate's duplication of as many bases ending at P-1, the bases
/// read on the transcript's strand.
#[test]
fn the_respelled_clamped_notations_are_read_against_the_duplication() {
    assert!(ONE_BASE.reads("n.2->A", "n.1dup"));
    assert!(ONE_BASE_REVERSE.reads("n.2->C", "n.1dup"));
    let in_the_utr = Row {
        location: "21:45455521-45455522",
        feature: "NM_030582.4",
        ..ONE_BASE
    };
    assert!(in_the_utr.reads("c.-8->A", "c.-9dup"));
    let no_status = Row {
        bam_edit: "",
        ..ONE_BASE
    };
    assert!(no_status.reads("n.2->A", "n.1dup"));
    assert!(TWO_BASES_REVERSE.reads("n.3_4insCG", "n.2_3dup"));
    let two_forward = Row {
        allele: "AC",
        strand: "1",
        offset: "3",
        ..TWO_BASES_REVERSE
    };
    assert!(two_forward.reads("n.3_4insAC", "n.2_3dup"));
    let two_reverse = Row {
        allele: "AC",
        ..TWO_BASES_REVERSE
    };
    assert!(two_reverse.reads("n.3_4insGT", "n.2_3dup"));
    for (row, reference, candidate) in [
        (ONE_BASE, "n.2->A", "n.1dup"),
        (ONE_BASE_REVERSE, "n.2->C", "n.1dup"),
        (TWO_BASES_REVERSE, "n.3_4insCG", "n.2_3dup"),
    ] {
        assert_eq!(
            row.compared(reference, candidate),
            ("green".to_string(), once(CLASS)),
            "{reference}"
        );
    }
}

/// Through the comparison the one-base row is set aside once under the class in the tab, VCF
/// and JSON formats as in the default format, and without the class it is a difference.
#[test]
fn the_respelled_substitution_is_set_aside_once_a_row_in_every_format() {
    assert_eq!(
        ONE_BASE.compared("n.2->A", "n.1dup"),
        ("green".to_string(), once(CLASS))
    );
    for fmt in ["tab", "vcf", "json"] {
        assert_eq!(
            under_class(CLASS, &in_format(fmt, "n.2->A"), &in_format(fmt, "n.1dup")),
            ("green".to_string(), once(CLASS)),
            "{fmt}"
        );
    }
    assert_eq!(
        verdict_of(
            &default_text(&[ONE_BASE.line("n.2->A")]),
            &default_text(&[ONE_BASE.line("n.1dup")]),
            &[]
        ),
        "red"
    );
}

/// The sides exchanged, a reference duplication against a candidate's respelled notation,
/// stay a difference in every format.
#[test]
fn the_respelled_notations_on_the_candidates_side_stay_a_difference() {
    assert!(ONE_BASE.reads("n.2->A", "n.1dup"));
    assert!(TWO_BASES_REVERSE.reads("n.3_4insCG", "n.2_3dup"));
    assert!(!ONE_BASE.reads("n.1dup", "n.2->A"));
    assert!(!TWO_BASES_REVERSE.reads("n.2_3dup", "n.3_4insCG"));
    assert_eq!(ONE_BASE.compared("n.1dup", "n.2->A"), red());
    for fmt in ["tab", "vcf", "json"] {
        assert_eq!(
            under_class(CLASS, &in_format(fmt, "n.1dup"), &in_format(fmt, "n.2->A")),
            red(),
            "{fmt}"
        );
    }
}

/// A row off the shape stays a difference: a `-` reference base before two bases or at
/// another position than the one after the duplication, a row that prints no SHIFT_LENGTH,
/// an `ins` whose bases are not the row's allele on the transcript's strand, hold more bases
/// than it or stand elsewhere, a whole-allele `ins` against a candidate `ins`, and a
/// whole-allele `ins` on a row that does not carry `BAM_EDIT` `OK`.
#[test]
fn a_respelled_clamped_notation_off_the_shape_stays_a_difference() {
    assert!(ONE_BASE.reads("n.2->A", "n.1dup"));
    assert!(TWO_BASES_REVERSE.reads("n.3_4insCG", "n.2_3dup"));
    for (reference, candidate) in [
        ("n.2->AC", "n.1dup"),
        ("n.3->A", "n.1dup"),
        ("n.2->A", "n.2dup"),
        ("n.2-->A", "n.1dup"),
        ("n.2->-", "n.1dup"),
        ("c.2->A", "n.1dup"),
    ] {
        assert!(
            !ONE_BASE.reads(reference, candidate),
            "{reference} {candidate}"
        );
        assert_eq!(
            ONE_BASE.compared(reference, candidate),
            red(),
            "{reference} {candidate}"
        );
    }
    let unlabelled = Row {
        shift: "",
        ..ONE_BASE
    };
    assert!(!unlabelled.reads("n.2->A", "n.1dup"));
    for (reference, candidate) in [
        ("n.3_4insGC", "n.2_3dup"),
        ("n.3_4insCGC", "n.2_3dup"),
        ("n.3_4insCG", "n.3_4dup"),
        ("n.4_5insCG", "n.2_3dup"),
        ("n.3_4insCG", "n.3_4insGC"),
        ("n.3_4insCG", "n.3dup"),
    ] {
        assert!(
            !TWO_BASES_REVERSE.reads(reference, candidate),
            "{reference} {candidate}"
        );
    }
    for bam_edit in ["", "FAILED"] {
        let row = Row {
            bam_edit,
            ..TWO_BASES_REVERSE
        };
        assert!(!row.reads("n.3_4insCG", "n.2_3dup"), "{bam_edit}");
        assert_eq!(row.compared("n.3_4insCG", "n.2_3dup"), red(), "{bam_edit}");
    }
    let two_forward = Row {
        allele: "AC",
        strand: "1",
        offset: "3",
        ..TWO_BASES_REVERSE
    };
    assert!(!two_forward.reads("n.3_4insGT", "n.2_3dup"));
    assert!(!two_forward.reads("n.3_4insCA", "n.2_3dup"));
}
