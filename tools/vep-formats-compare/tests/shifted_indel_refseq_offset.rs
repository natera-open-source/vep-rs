// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `REFSEQ_OFFSET` of `hgvs_shifted_indel_written_span_gate`: on a row carrying
//! `BAM_EDIT` `OK` whose HGVSc and HGVS_OFFSET the class reads, the reference writing
//! neither, the offset of the model's edits that the candidate alone writes is set aside with
//! them. The class is asked as a comparison asks it, and two outputs are compared in the
//! default, tab, VCF and JSON formats.

mod compare_common;
mod shifted_edited_common;

use std::collections::BTreeMap;

use compare_common::*;
use shifted_edited_common::*;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::downstream::field_of;

const GATE: &str = "hgvs_shifted_indel_written_span_gate";
const FEATURE: &str = "NM_001331018.1";
const INPUT: &str = "21\t37223419\tlast_two\tATA\tA\t.\t.\t.";

/// What one side of a row writes of the fields the class reads.
#[derive(Clone, Copy)]
struct Written<'a> {
    hgvsc: &'a str,
    hgvs_offset: &'a str,
    refseq_offset: &'a str,
    bam_edit: &'a str,
}

/// The reference's side of the deletion of an edited model's last two bases: the 3' shift
/// carries the deletion over the model's last base, and the row has no notation, no
/// HGVS_OFFSET and no REFSEQ_OFFSET.
const REFERENCE: Written = Written {
    hgvsc: "",
    hgvs_offset: "",
    refseq_offset: "",
    bam_edit: "OK",
};

/// The candidate's side of that row: the shifted deletion, its shift, and the one base the
/// model's edits put before the allele.
const CANDIDATE: Written = Written {
    hgvsc: "NM_001331018.1:c.*2124del",
    hgvs_offset: "-1",
    refseq_offset: "1",
    bam_edit: "OK",
};

/// The cells of a row that both sides write alike.
#[derive(Clone, Copy)]
struct Row<'a> {
    location: &'a str,
    allele: &'a str,
    terms: &'a str,
    cdna: &'a str,
    strand: &'a str,
}

const DELETION: Row = Row {
    location: "21:37223420-37223421",
    allele: "-",
    terms: "3_prime_UTR_variant",
    cdna: "3457-3458",
    strand: "-1",
};

/// A made-up row of the same model read on the forward strand: one base duplicated behind an
/// edit that deleted bases.
const DUPLICATION: Row = Row {
    location: "21:45513719-45513720",
    allele: "A",
    terms: "3_prime_UTR_variant",
    cdna: "6586-6587",
    strand: "1",
};

impl Row<'_> {
    /// The row's line of a default-format output with what `side` writes.
    fn line(&self, side: Written) -> String {
        let mut extra = format!("IMPACT=MODIFIER;STRAND={}", self.strand);
        for (name, value) in [
            ("REFSEQ_OFFSET", side.refseq_offset),
            ("BAM_EDIT", side.bam_edit),
            ("HGVSc", side.hgvsc),
            ("HGVS_OFFSET", side.hgvs_offset),
        ] {
            if !value.is_empty() {
                extra.push_str(&format!(";{name}={value}"));
            }
        }
        format!(
            "last_two\t{}\t{}\t10311\t{FEATURE}\tTranscript\t{}\t{}\t-\t-\t-\t-\t-\t{extra}",
            self.location, self.allele, self.terms, self.cdna
        )
    }

    fn entry(&self, side: Written) -> Entry {
        default_entry(&self.line(side))
    }

    /// What the class says of the row's REFSEQ_OFFSET, the reference's side `reference` and
    /// the candidate's `candidate`.
    fn reads(&self, reference: Written, candidate: Written) -> bool {
        let (e, a) = (self.entry(reference), self.entry(candidate));
        let name = "REFSEQ_OFFSET";
        class_accepts(
            GATE,
            &e,
            name,
            &field_of(Some(&e), &[name]),
            &field_of(Some(&a), &[name]),
            Some(&a),
        )
    }

    /// The verdict of the row's two default-format outputs under the class and what was
    /// applied.
    fn compared(&self, reference: Written, candidate: Written) -> (String, BTreeMap<String, i64>) {
        default_under(
            GATE,
            &default_text(&[self.line(reference)]),
            &default_text(&[self.line(candidate)]),
        )
    }
}

/// The deletion's row in the tab, VCF and JSON formats with what `side` writes.
fn in_format(fmt: &str, side: Written) -> Vec<Entry> {
    let inputs = vec![INPUT.to_string()];
    let cell = |value: &str| if value.is_empty() { "-" } else { value }.to_string();
    match fmt {
        "tab" => {
            let header = tab_header(&[
                "IMPACT",
                "STRAND",
                "REFSEQ_OFFSET",
                "BAM_EDIT",
                "HGVSc",
                "HGVS_OFFSET",
            ]);
            let cells = [
                "MODIFIER".to_string(),
                DELETION.strand.to_string(),
                cell(side.refseq_offset),
                cell(side.bam_edit),
                cell(side.hgvsc),
                cell(side.hgvs_offset),
            ];
            let row = format!(
                "last_two\t{}\t-\t10311\t{FEATURE}\tTranscript\t{}\t{}\t-\t-\t-\t-\t-\t{}",
                DELETION.location,
                DELETION.terms,
                DELETION.cdna,
                cells.join("\t")
            );
            tab_entries(&format!("{header}\n{row}\n"))
        }
        "vcf" => {
            let header = csq_header(
                "Allele|Consequence|Feature|Feature_type|cDNA_position|REFSEQ_OFFSET|BAM_EDIT|HGVSc|HGVS_OFFSET",
            );
            let csq = [
                "-",
                DELETION.terms,
                FEATURE,
                "Transcript",
                DELETION.cdna,
                side.refseq_offset,
                side.bam_edit,
                side.hgvsc,
                side.hgvs_offset,
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
                ("id", text("last_two")),
                ("allele_string", text("TA/-")),
                ("most_severe_consequence", text(DELETION.terms)),
            ];
            let mut consequence = vec![
                ("variant_allele", text("-")),
                ("transcript_id", text(FEATURE)),
                ("consequence_terms", Value::List(vec![text(DELETION.terms)])),
                ("cdna_start", Value::Int(3457)),
                ("cdna_end", Value::Int(3458)),
            ];
            let whole = |value: &str| Value::Int(value.parse().expect("a whole number"));
            if !side.refseq_offset.is_empty() {
                consequence.push(("refseq_offset", whole(side.refseq_offset)));
            }
            if !side.bam_edit.is_empty() {
                consequence.push(("bam_edit", text(side.bam_edit)));
            }
            if !side.hgvsc.is_empty() {
                consequence.push(("hgvsc", text(side.hgvsc)));
            }
            if !side.hgvs_offset.is_empty() {
                consequence.push(("hgvs_offset", whole(side.hgvs_offset)));
            }
            json_entries(&json_line(&record, object(&consequence)), &inputs)
        }
    }
}

fn red() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), nothing())
}

/// The REFSEQ_OFFSET the candidate alone writes is the class's on the row whose notation and
/// HGVS_OFFSET it reads: the deletion of an edited model's last two bases on a reverse-strand
/// model, and a duplication beside the negative offset of a model whose edit deleted bases.
#[test]
fn the_written_span_gate_reads_the_refseq_offset_the_candidate_alone_writes() {
    assert!(DELETION.reads(REFERENCE, CANDIDATE));
    assert_eq!(
        DELETION.compared(REFERENCE, CANDIDATE),
        ("green".to_string(), once(GATE))
    );
    let duplicated = Written {
        hgvsc: "NM_001331018.1:c.*1322dup",
        hgvs_offset: "1",
        refseq_offset: "-9",
        ..CANDIDATE
    };
    assert!(DUPLICATION.reads(REFERENCE, duplicated));
    assert_eq!(
        DUPLICATION.compared(REFERENCE, duplicated),
        ("green".to_string(), once(GATE))
    );
}

/// Through the comparison the row is set aside once under the class in the tab, VCF and JSON
/// formats as in the default format, and without the class it is a difference.
#[test]
fn the_refseq_offset_is_set_aside_with_the_notation_once_a_row_in_every_format() {
    assert_eq!(
        DELETION.compared(REFERENCE, CANDIDATE),
        ("green".to_string(), once(GATE))
    );
    for fmt in ["tab", "vcf", "json"] {
        assert_eq!(
            under_class(GATE, &in_format(fmt, REFERENCE), &in_format(fmt, CANDIDATE)),
            ("green".to_string(), once(GATE)),
            "{fmt}"
        );
    }
    assert_eq!(
        verdict_of(
            &default_text(&[DELETION.line(REFERENCE)]),
            &default_text(&[DELETION.line(CANDIDATE)]),
            &[]
        ),
        "red"
    );
}

/// The sides exchanged, a reference that alone writes the notation, the HGVS_OFFSET and the
/// REFSEQ_OFFSET, stay a difference in every format.
#[test]
fn the_refseq_offset_the_reference_alone_writes_stays_a_difference() {
    assert!(DELETION.reads(REFERENCE, CANDIDATE));
    assert!(!DELETION.reads(CANDIDATE, REFERENCE));
    assert_eq!(DELETION.compared(CANDIDATE, REFERENCE), red());
    for fmt in ["tab", "vcf", "json"] {
        assert_eq!(
            under_class(GATE, &in_format(fmt, CANDIDATE), &in_format(fmt, REFERENCE)),
            red(),
            "{fmt}"
        );
    }
    // The reference writes the offset alone, beside a candidate row of the class's shape.
    let offset_alone = Written {
        refseq_offset: "1",
        ..REFERENCE
    };
    let without = Written {
        refseq_offset: "",
        ..CANDIDATE
    };
    assert!(!DELETION.reads(offset_alone, without));
    assert_eq!(DELETION.compared(offset_alone, without).0, "red");
}

/// A REFSEQ_OFFSET off the shape stays a difference: two different values, a value on a row
/// the class does not read (a reference notation, a candidate without an HGVS_OFFSET or with
/// a substitution's notation, a substitution's row, a row annotated downstream), a row that
/// does not carry `BAM_EDIT` `OK` on both sides, and a value that is no non-zero integer.
#[test]
fn a_refseq_offset_off_the_shape_stays_a_difference() {
    assert!(DELETION.reads(REFERENCE, CANDIDATE));
    let other_value = Written {
        refseq_offset: "19",
        ..REFERENCE
    };
    assert!(!DELETION.reads(other_value, CANDIDATE));
    assert_eq!(DELETION.compared(other_value, CANDIDATE).0, "red");

    let noted = Written {
        hgvsc: "NM_001331018.1:c.*2123del",
        ..REFERENCE
    };
    let unread: [(&str, Written, Written); 6] = [
        ("a reference notation", noted, CANDIDATE),
        (
            "no candidate HGVS_OFFSET",
            REFERENCE,
            Written {
                hgvs_offset: "",
                ..CANDIDATE
            },
        ),
        (
            "a zero HGVS_OFFSET",
            REFERENCE,
            Written {
                hgvs_offset: "0",
                ..CANDIDATE
            },
        ),
        (
            "no candidate notation",
            REFERENCE,
            Written {
                hgvsc: "",
                ..CANDIDATE
            },
        ),
        (
            "a substitution's notation",
            REFERENCE,
            Written {
                hgvsc: "NM_001331018.1:c.*2124A>T",
                ..CANDIDATE
            },
        ),
        (
            "another transcript's notation",
            REFERENCE,
            Written {
                hgvsc: "NM_001331021.1:c.*2124del",
                ..CANDIDATE
            },
        ),
    ];
    for (label, reference, candidate) in unread {
        assert!(!DELETION.reads(reference, candidate), "{label}");
        assert_eq!(DELETION.compared(reference, candidate).0, "red", "{label}");
    }
    let substitution = Row {
        location: "21:37223420",
        allele: "T",
        cdna: "3458",
        ..DELETION
    };
    assert!(!substitution.reads(REFERENCE, CANDIDATE));
    let downstream = Row {
        terms: "downstream_gene_variant",
        cdna: "-",
        ..DELETION
    };
    assert!(!downstream.reads(REFERENCE, CANDIDATE));

    for (reference_edit, candidate_edit) in [("", ""), ("FAILED", "FAILED"), ("OK", ""), ("", "OK")]
    {
        let reference = Written {
            bam_edit: reference_edit,
            ..REFERENCE
        };
        let candidate = Written {
            bam_edit: candidate_edit,
            ..CANDIDATE
        };
        let label = format!("{reference_edit} {candidate_edit}");
        assert!(!DELETION.reads(reference, candidate), "{label}");
        assert_eq!(DELETION.compared(reference, candidate).0, "red", "{label}");
    }
    for value in ["0", "+1", "1.5", "one"] {
        let candidate = Written {
            refseq_offset: value,
            ..CANDIDATE
        };
        assert!(!DELETION.reads(REFERENCE, candidate), "{value}");
    }
}
