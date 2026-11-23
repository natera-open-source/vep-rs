// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvs_clipped_reference_against_whole_alternate` on an edited model's last base: the
//! reference writes the whole alternate against the one base its clamp left, at a coordinate
//! the model's REFSEQ_OFFSET slid past that base, and the candidate writes the change inside
//! at the transcript's own coordinate. The class is asked as a comparison asks it, and two
//! outputs are compared in the default, tab, VCF and JSON formats.

mod compare_common;
mod shifted_edited_common;

use std::collections::BTreeMap;

use compare_common::*;
use shifted_edited_common::*;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::pyjson::Value;

const CLIPPED: &str = "hgvs_clipped_reference_against_whole_alternate";
const OFFSET: &str = "refseq_edited_model_hgvsc_offset";
const TAIL_TERM: &str = "refseq_edited_model_tail_term";
const FEATURE: &str = "NM_001331018.1";
const INPUT: &str = "21\t37223419\tlast_base\tAT\tAA\t.\t.\t.";
/// The reference's notation and the candidate's.
const WHOLE_ALTERNATE: &str = "NM_001331018.1:c.*2123+19delinsTT";
const BASE_INSIDE: &str = "NM_001331018.1:c.*2124A>T";

/// The cells of a row that the class reads beside its two notations.
#[derive(Clone, Copy)]
struct Row<'a> {
    location: &'a str,
    allele: &'a str,
    cdna: &'a str,
    refseq_offset: &'a str,
    bam_edit: &'a str,
}

/// Two bases replaced by two over the last base of a reverse-strand model that holds one
/// inserted base and eighteen appended ones: the base inside changes and the base past it
/// does not.
const LAST_BASE: Row = Row {
    location: "21:37223419-37223420",
    allele: "AA",
    cdna: "3458-?",
    refseq_offset: "19",
    bam_edit: "OK",
};

impl Row<'_> {
    /// The row's line of a default-format output under `terms` with the notation `hgvsc`.
    fn line_of(&self, terms: &str, hgvsc: &str) -> String {
        let mut extra = String::from("IMPACT=MODIFIER;STRAND=-1");
        for (name, value) in [
            ("REFSEQ_OFFSET", self.refseq_offset),
            ("BAM_EDIT", self.bam_edit),
            ("HGVSc", hgvsc),
        ] {
            if !value.is_empty() {
                extra.push_str(&format!(";{name}={value}"));
            }
        }
        format!(
            "last_base\t{}\t{}\t10311\t{FEATURE}\tTranscript\t{terms}\t{}\t-\t-\t-\t-\t-\t{extra}",
            self.location, self.allele, self.cdna
        )
    }

    fn line(&self, hgvsc: &str) -> String {
        self.line_of("3_prime_UTR_variant", hgvsc)
    }

    /// What the class says of the row's HGVSc, the reference's `reference` and the
    /// candidate's `candidate`.
    fn reads(&self, reference: &str, candidate: &str) -> bool {
        let (e, a) = (
            default_entry(&self.line(reference)),
            default_entry(&self.line(candidate)),
        );
        class_accepts(CLIPPED, &e, "HGVSc", reference, candidate, Some(&a))
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

/// The row of `LAST_BASE` in the tab, VCF and JSON formats with the notation `hgvsc`.
fn in_format(fmt: &str, hgvsc: &str) -> Vec<Entry> {
    let inputs = vec![INPUT.to_string()];
    match fmt {
        "tab" => {
            let header = tab_header(&["IMPACT", "STRAND", "REFSEQ_OFFSET", "BAM_EDIT", "HGVSc"]);
            let row = format!(
                "last_base\t{}\t{}\t10311\t{FEATURE}\tTranscript\t3_prime_UTR_variant\t{}\t-\t-\t-\t-\t-\tMODIFIER\t-1\t{}\t{}\t{hgvsc}",
                LAST_BASE.location,
                LAST_BASE.allele,
                LAST_BASE.cdna,
                LAST_BASE.refseq_offset,
                LAST_BASE.bam_edit
            );
            tab_entries(&format!("{header}\n{row}\n"))
        }
        "vcf" => {
            let header = csq_header(
                "Allele|Consequence|Feature|Feature_type|cDNA_position|REFSEQ_OFFSET|BAM_EDIT|HGVSc",
            );
            let csq = [
                LAST_BASE.allele,
                "3_prime_UTR_variant",
                FEATURE,
                "Transcript",
                LAST_BASE.cdna,
                LAST_BASE.refseq_offset,
                LAST_BASE.bam_edit,
                hgvsc,
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
                ("allele_string", text("AT/AA")),
                ("most_severe_consequence", text("3_prime_UTR_variant")),
            ];
            let consequence = object(&[
                ("variant_allele", text(LAST_BASE.allele)),
                ("transcript_id", text(FEATURE)),
                (
                    "consequence_terms",
                    Value::List(vec![text("3_prime_UTR_variant")]),
                ),
                ("cdna_start", Value::Int(3458)),
                ("refseq_offset", Value::Int(19)),
                ("bam_edit", text(LAST_BASE.bam_edit)),
                ("hgvsc", text(hgvsc)),
            ]);
            json_entries(&json_line(&record, consequence), &inputs)
        }
    }
}

fn red() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), nothing())
}

/// The reference's `delins` of the whole alternate at the slid coordinate against the
/// candidate's substitution of the base inside is the class's: the record that changes the
/// base inside alone, the record that changes both bases, and the same pair on a model whose
/// cDNA numbers differ. The offset class, which reads two notations of one kind, refuses it.
#[test]
fn the_whole_alternate_at_a_slid_coordinate_is_read_against_the_base_inside() {
    assert!(LAST_BASE.reads(WHOLE_ALTERNATE, BASE_INSIDE));
    assert_eq!(
        LAST_BASE.compared(WHOLE_ALTERNATE, BASE_INSIDE),
        ("green".to_string(), once(CLIPPED))
    );
    let both_change = Row {
        allele: "CA",
        ..LAST_BASE
    };
    assert!(both_change.reads("NM_001331018.1:c.*2123+19delinsTG", BASE_INSIDE));
    let other_model = Row {
        cdna: "3689-?",
        ..LAST_BASE
    };
    assert!(other_model.reads(WHOLE_ALTERNATE, BASE_INSIDE));
    let (e, a) = (
        default_entry(&LAST_BASE.line(WHOLE_ALTERNATE)),
        default_entry(&LAST_BASE.line(BASE_INSIDE)),
    );
    assert!(!class_accepts(
        OFFSET,
        &e,
        "HGVSc",
        WHOLE_ALTERNATE,
        BASE_INSIDE,
        Some(&a)
    ));
}

/// Through the comparison the row is set aside once under the class in the default, tab, VCF
/// and JSON formats; without the class it is a difference, and with the terms each engine
/// writes it is green under the class and the class of the model's last base.
#[test]
fn the_slid_coordinate_row_is_set_aside_once_in_every_format() {
    assert_eq!(
        LAST_BASE.compared(WHOLE_ALTERNATE, BASE_INSIDE),
        ("green".to_string(), once(CLIPPED))
    );
    for fmt in ["tab", "vcf", "json"] {
        assert_eq!(
            under_class(
                CLIPPED,
                &in_format(fmt, WHOLE_ALTERNATE),
                &in_format(fmt, BASE_INSIDE)
            ),
            ("green".to_string(), once(CLIPPED)),
            "{fmt}"
        );
    }
    let reference = default_text(&[LAST_BASE.line(WHOLE_ALTERNATE)]);
    let candidate = default_text(&[LAST_BASE.line(BASE_INSIDE)]);
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    assert_eq!(verdict_of(&reference, &candidate, &[OFFSET]), "red");
    let (report, documented) = compare_default(
        &default_text(&[LAST_BASE.line_of("intergenic_variant", WHOLE_ALTERNATE)]),
        &default_text(&[LAST_BASE.line_of("3_prime_UTR_variant", BASE_INSIDE)]),
        &[CLIPPED, TAIL_TERM],
    );
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(CLIPPED, 1), (TAIL_TERM, 1)])
    );
}

/// The sides exchanged, a reference substitution at the transcript's own coordinate against
/// a candidate `delins` at the slid one, stay a difference in the default and tab formats,
/// whose rows state the span.
#[test]
fn the_slid_coordinate_on_the_candidates_side_stays_a_difference() {
    assert!(LAST_BASE.reads(WHOLE_ALTERNATE, BASE_INSIDE));
    assert!(!LAST_BASE.reads(BASE_INSIDE, WHOLE_ALTERNATE));
    assert_eq!(LAST_BASE.compared(BASE_INSIDE, WHOLE_ALTERNATE), red());
    assert_eq!(
        under_class(
            CLIPPED,
            &in_format("tab", BASE_INSIDE),
            &in_format("tab", WHOLE_ALTERNATE)
        ),
        red()
    );
}

/// A row off the shape stays a difference: a distance that is not the row's REFSEQ_OFFSET, a
/// row without one or with a negative one, a row that does not carry `BAM_EDIT` `OK`, a cDNA
/// position that is closed or open at its start, and a candidate coordinate further than the
/// span's length plus one from the stretch between the exon base the reference's coordinate
/// counts from and that coordinate.
#[test]
fn a_slid_coordinate_off_the_shape_stays_a_difference() {
    assert!(LAST_BASE.reads(WHOLE_ALTERNATE, BASE_INSIDE));
    for reference in [
        "NM_001331018.1:c.*2123+18delinsTT",
        "NM_001331018.1:c.*2123+20delinsTT",
        "NM_001331018.1:c.*2123-19delinsTT",
        "NM_001331018.1:c.*2142delinsTT",
    ] {
        assert!(!LAST_BASE.reads(reference, BASE_INSIDE), "{reference}");
        assert_eq!(
            LAST_BASE.compared(reference, BASE_INSIDE),
            red(),
            "{reference}"
        );
    }
    let off = [
        Row {
            refseq_offset: "",
            ..LAST_BASE
        },
        Row {
            refseq_offset: "18",
            ..LAST_BASE
        },
        Row {
            refseq_offset: "-19",
            ..LAST_BASE
        },
        Row {
            bam_edit: "",
            ..LAST_BASE
        },
        Row {
            bam_edit: "FAILED",
            ..LAST_BASE
        },
        Row {
            cdna: "3458",
            ..LAST_BASE
        },
        Row {
            cdna: "?-3458",
            ..LAST_BASE
        },
    ];
    for row in off {
        let label = format!("{} {} {}", row.refseq_offset, row.bam_edit, row.cdna);
        assert!(!row.reads(WHOLE_ALTERNATE, BASE_INSIDE), "{label}");
        assert_eq!(row.compared(WHOLE_ALTERNATE, BASE_INSIDE), red(), "{label}");
    }
    for (candidate, read) in [
        ("NM_001331018.1:c.*2119A>T", false),
        ("NM_001331018.1:c.*2120A>T", true),
        ("NM_001331018.1:c.*2123A>T", true),
        ("NM_001331018.1:c.*2142A>T", true),
        ("NM_001331018.1:c.*2145A>T", true),
        ("NM_001331018.1:c.*2146A>T", false),
        ("NM_001331018.1:c.2124A>T", false),
        ("NM_001331018.1:n.3458A>T", false),
        ("NM_001331021.1:c.*2124A>T", false),
    ] {
        assert_eq!(
            LAST_BASE.reads(WHOLE_ALTERNATE, candidate),
            read,
            "{candidate}"
        );
    }
}
