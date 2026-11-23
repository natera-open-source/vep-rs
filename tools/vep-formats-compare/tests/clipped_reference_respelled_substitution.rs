// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvs_clipped_reference_against_whole_alternate` on an edited model's edge, where two
//! bases, one inside the model and one outside it, are replaced by one: the reference types
//! a substitution of the base inside and then writes the row's whole reference allele over
//! it, two bases at one position, and the candidate names the base inside alone. The class
//! is asked as a comparison asks it, and two outputs are compared in the default, tab, VCF
//! and JSON formats.

mod compare_common;
mod shifted_edited_common;

use std::collections::BTreeMap;

use compare_common::*;
use shifted_edited_common::*;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::pyjson::Value;

const CLIPPED: &str = "hgvs_clipped_reference_against_whole_alternate";
const UNMINIMISED: &str = "refseq_edited_model_unminimised_allele_hgvsc";
const INPUT: &str = "21\t42297244\tlast_base\tAC\tG\t.\t.\t.";

/// The cells of a row that the class reads beside its two notations.
#[derive(Clone, Copy)]
struct Row<'a> {
    location: &'a str,
    feature: &'a str,
    cdna: &'a str,
    strand: &'a str,
    /// The row's reference allele as the genome and the model spell it.
    used_ref: &'a str,
    refseq_offset: &'a str,
    bam_edit: &'a str,
}

/// The last base of a forward-strand model and the base past it replaced by one base.
const LAST_BASE: Row = Row {
    location: "21:42297244-42297245",
    feature: "NM_207628.1",
    cdna: "3037-?",
    strand: "1",
    used_ref: "AC",
    refseq_offset: "",
    bam_edit: "OK",
};

impl Row<'_> {
    /// The row's line of a default-format output with the notation after `feature:`.
    fn line(&self, notation: &str) -> String {
        let mut extra = format!("IMPACT=MODIFIER;STRAND={}", self.strand);
        let hgvsc = format!("{}:{notation}", self.feature);
        for (name, value) in [
            ("REFSEQ_OFFSET", self.refseq_offset),
            ("GIVEN_REF", self.used_ref),
            ("USED_REF", self.used_ref),
            ("BAM_EDIT", self.bam_edit),
            ("HGVSc", hgvsc.as_str()),
        ] {
            if !value.is_empty() {
                extra.push_str(&format!(";{name}={value}"));
            }
        }
        format!(
            "last_base\t{}\tG\tENSG1\t{}\tTranscript\t3_prime_UTR_variant\t{}\t-\t-\t-\t-\t-\t{extra}",
            self.location, self.feature, self.cdna
        )
    }

    /// What `klass` says of the row's HGVSc, the reference's notation `reference` and the
    /// candidate's `candidate`.
    fn read_by(&self, klass: &str, reference: &str, candidate: &str) -> bool {
        let (e, a) = (
            default_entry(&self.line(reference)),
            default_entry(&self.line(candidate)),
        );
        let value = |notation: &str| format!("{}:{notation}", self.feature);
        class_accepts(
            klass,
            &e,
            "HGVSc",
            &value(reference),
            &value(candidate),
            Some(&a),
        )
    }

    fn reads(&self, reference: &str, candidate: &str) -> bool {
        self.read_by(CLIPPED, reference, candidate)
    }

    /// The verdict of the row's two default-format outputs under the class and what was
    /// applied.
    fn compared(&self, reference: &str, candidate: &str) -> (String, BTreeMap<String, i64>) {
        default_under(
            CLIPPED,
            &default_text(&[self.line(reference)]),
            &default_text(&[self.line(candidate)]),
        )
    }
}

/// The row of `LAST_BASE` in the tab, VCF and JSON formats with the notation `notation`.
fn in_format(fmt: &str, notation: &str) -> Vec<Entry> {
    let inputs = vec![INPUT.to_string()];
    let hgvsc = format!("{}:{notation}", LAST_BASE.feature);
    match fmt {
        "tab" => {
            let header = tab_header(&["IMPACT", "STRAND", "BAM_EDIT", "HGVSc"]);
            let row = format!(
                "last_base\t{}\tG\tENSG1\t{}\tTranscript\t3_prime_UTR_variant\t{}\t-\t-\t-\t-\t-\tMODIFIER\t1\tOK\t{hgvsc}",
                LAST_BASE.location, LAST_BASE.feature, LAST_BASE.cdna
            );
            tab_entries(&format!("{header}\n{row}\n"))
        }
        "vcf" => {
            let header =
                csq_header("Allele|Consequence|Feature|Feature_type|cDNA_position|BAM_EDIT|HGVSc");
            let csq = [
                "G",
                "3_prime_UTR_variant",
                LAST_BASE.feature,
                "Transcript",
                LAST_BASE.cdna,
                "OK",
                &hgvsc,
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
                ("id", text("last_base")),
                ("allele_string", text("AC/G")),
                ("most_severe_consequence", text("3_prime_UTR_variant")),
            ];
            let consequence = object(&[
                ("variant_allele", text("G")),
                ("transcript_id", text(LAST_BASE.feature)),
                (
                    "consequence_terms",
                    Value::List(vec![text("3_prime_UTR_variant")]),
                ),
                ("cdna_start", Value::Int(3037)),
                ("bam_edit", text("OK")),
                ("hgvsc", text(&hgvsc)),
            ]);
            json_entries(&json_line(&record, consequence), &inputs)
        }
    }
}

fn red() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), nothing())
}

/// The reference's substitution naming the row's two reference bases at one position against
/// the candidate's substitution of the base inside is the class's: on either strand, in `c.`
/// and `n.` numbering, with the candidate's coordinate one above the reference's as the
/// model's numbering puts it, and with the reference's coordinate slid past the last base by
/// the row's REFSEQ_OFFSET. The class of the unclipped pair, which reads a candidate spelt
/// without bases the two alleles share, refuses it.
#[test]
fn the_respelled_reference_allele_is_read_against_the_base_inside() {
    assert!(LAST_BASE.reads("c.*852AC>G", "c.*852A>G"));
    let reverse = Row {
        location: "21:38367260-38367261",
        feature: "NM_001243432.2",
        cdna: "1600-?",
        strand: "-1",
        used_ref: "AC",
        ..LAST_BASE
    };
    assert!(reverse.reads("c.*374GT>C", "c.*374G>C"));
    let non_coding = Row {
        location: "21:42849225-42849226",
        feature: "NR_048535.1",
        cdna: "2042-?",
        strand: "-1",
        used_ref: "GG",
        ..LAST_BASE
    };
    assert!(non_coding.reads("n.2042CC>T", "n.2042C>T"));
    let renumbered = Row {
        location: "21:31653586-31653587",
        feature: "NR_187558.1",
        cdna: "2632-?",
        strand: "-1",
        used_ref: "TA",
        ..LAST_BASE
    };
    assert!(renumbered.reads("n.2631TA>G", "n.2632T>G"));
    let slid = Row {
        location: "21:37223419-37223420",
        feature: "NM_001331018.1",
        cdna: "3458-?",
        strand: "-1",
        used_ref: "AT",
        refseq_offset: "19",
        ..LAST_BASE
    };
    assert!(slid.reads("c.*2123+19AT>G", "c.*2124A>G"));
    assert_eq!(
        LAST_BASE.compared("c.*852AC>G", "c.*852A>G"),
        ("green".to_string(), once(CLIPPED))
    );
    assert_eq!(
        slid.compared("c.*2123+19AT>G", "c.*2124A>G"),
        ("green".to_string(), once(CLIPPED))
    );
    assert!(!LAST_BASE.read_by(UNMINIMISED, "c.*852AC>G", "c.*852A>G"));
}

/// Through the comparison the row is set aside once under the class in the tab, VCF and JSON
/// formats as in the default format, and without the class it is a difference.
#[test]
fn the_respelled_reference_allele_is_set_aside_once_a_row_in_every_format() {
    assert_eq!(
        LAST_BASE.compared("c.*852AC>G", "c.*852A>G"),
        ("green".to_string(), once(CLIPPED))
    );
    for fmt in ["tab", "vcf", "json"] {
        assert_eq!(
            under_class(
                CLIPPED,
                &in_format(fmt, "c.*852AC>G"),
                &in_format(fmt, "c.*852A>G")
            ),
            ("green".to_string(), once(CLIPPED)),
            "{fmt}"
        );
    }
    assert_eq!(
        verdict_of(
            &default_text(&[LAST_BASE.line("c.*852AC>G")]),
            &default_text(&[LAST_BASE.line("c.*852A>G")]),
            &[]
        ),
        "red"
    );
}

/// The sides exchanged, a reference that names the base inside against a candidate that
/// names two reference bases at one position, stay a difference in every format.
#[test]
fn the_respelled_reference_allele_on_the_candidates_side_stays_a_difference() {
    assert!(LAST_BASE.reads("c.*852AC>G", "c.*852A>G"));
    assert!(!LAST_BASE.reads("c.*852A>G", "c.*852AC>G"));
    assert_eq!(LAST_BASE.compared("c.*852A>G", "c.*852AC>G"), red());
    for fmt in ["tab", "vcf", "json"] {
        assert_eq!(
            under_class(
                CLIPPED,
                &in_format(fmt, "c.*852A>G"),
                &in_format(fmt, "c.*852AC>G")
            ),
            red(),
            "{fmt}"
        );
    }
}

/// A row off the shape stays a difference: a row that does not carry `BAM_EDIT` `OK`, a
/// closed cDNA position, a reference that writes two alternate bases or a letter that names
/// no base, a candidate at a coordinate further than the span's length plus one from the
/// reference's, and another transcript's or numbering's notation.
#[test]
fn a_respelled_reference_allele_off_the_shape_stays_a_difference() {
    assert!(LAST_BASE.reads("c.*852AC>G", "c.*852A>G"));
    for bam_edit in ["", "FAILED"] {
        let row = Row {
            bam_edit,
            ..LAST_BASE
        };
        assert!(!row.reads("c.*852AC>G", "c.*852A>G"), "{bam_edit}");
        assert_eq!(row.compared("c.*852AC>G", "c.*852A>G"), red(), "{bam_edit}");
    }
    let closed = Row {
        cdna: "3037",
        ..LAST_BASE
    };
    assert!(!closed.reads("c.*852AC>G", "c.*852A>G"));
    for (reference, candidate) in [
        ("c.*852AC>GT", "c.*852A>G"),
        ("c.*852AR>G", "c.*852A>G"),
        ("c.*852AC>-", "c.*852A>G"),
        ("c.*852AC>G", "c.*856A>G"),
        ("c.*852AC>G", "c.*848A>G"),
        ("c.*852AC>G", "n.3037A>G"),
        ("c.852AC>G", "c.*852A>G"),
    ] {
        assert!(
            !LAST_BASE.reads(reference, candidate),
            "{reference} {candidate}"
        );
        assert_eq!(
            LAST_BASE.compared(reference, candidate),
            red(),
            "{reference} {candidate}"
        );
    }
    let other = format!("{}:c.*852A>G", "NM_207627.2");
    let (e, a) = (
        default_entry(&LAST_BASE.line("c.*852AC>G")),
        default_entry(&LAST_BASE.line("c.*852A>G")),
    );
    assert!(!class_accepts(
        CLIPPED,
        &e,
        "HGVSc",
        "NM_207628.1:c.*852AC>G",
        &other,
        Some(&a)
    ));
}
