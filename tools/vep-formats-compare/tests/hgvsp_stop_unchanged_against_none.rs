// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tenth shape of `hgvsp_stop_codon_protein_comparison`: the reference writes no HGVSp
//! and the candidate writes the stop as unchanged at the residue after the one codon the
//! row names, on a row whose terms hold `stop_retained_variant`, whose coding positions lie
//! inside that codon and whose amino acids are the codon's residue against the same residue
//! and then a stop. The class is asked as a comparison asks it, and two outputs are compared
//! in the default, tab, VCF and JSON formats.

mod compare_common;

use compare_common::*;
use indexmap::IndexMap;
use std::collections::BTreeSet;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::{field_class_accepts, RuleContext, Side};

const STOP_CODON_CLASS: &str = "hgvsp_stop_codon_protein_comparison";
const PROTEIN: &str = "ENSP00000000001.1";
const INPUT: &str = "21\t25004296\tstop_row\tA\tACTGAGG\t.\t.\t.";
const KEPT_STOP_TERMS: &str = "inframe_insertion,stop_retained_variant";

/// The cells of a row that the shape reads beside its terms and its FLAGS.
#[derive(Clone, Copy)]
struct Row {
    /// The `CDS_position` cell: one coding base, or the two an insertion lies between.
    cds: &'static str,
    /// The `Protein_position` cell.
    protein: &'static str,
    /// The `Amino_acids` cell.
    amino_acids: &'static str,
}

/// Six bases inserted between the second and third base of codon 365, the last residue's,
/// of a protein whose stop is residue 366: the codon's residue is kept and a stop follows.
const SIX_BASES_FORWARD: Row = Row {
    cds: "1094-1095",
    protein: "365",
    amino_acids: "V/V*G",
};
const KEPT_STOP: &str = "Ter366=";

/// An HGVSp value: the protein, `:p.` and the description, empty for no description; the
/// default, tab and VCF formats escape the `=` of a description, the JSON format does not.
fn hgvsp(description: &str, escaped: bool) -> String {
    if description.is_empty() {
        return String::new();
    }
    let text = format!("{PROTEIN}:p.{description}");
    if escaped {
        text.replace('=', "%3D")
    } else {
        text
    }
}

/// What the class says of a row whose reference HGVSp is the description `reference` and
/// whose candidate HGVSp is `candidate` (empty for none), asked as a comparison without side
/// tables asks it.
fn reads_on(row: Row, terms: &str, flags: &str, reference: &str, candidate: &str) -> bool {
    let mut held = vec![
        ("CDS_position", row.cds),
        ("Protein_position", row.protein),
        ("Amino_acids", row.amino_acids),
    ];
    if !flags.is_empty() {
        held.push(("FLAGS", flags));
    }
    let e = entry(
        "21:25004296-25004297|CTGAGG|ENST00000000001|Transcript",
        None,
        Some("21:25004296-25004297"),
        "CTGAGG",
        "ENST00000000001",
        terms,
        &held,
    );
    let side = Side::default();
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    field_class_accepts(
        STOP_CODON_CLASS,
        &e,
        "HGVSp",
        &hgvsp(reference, true),
        &hgvsp(candidate, true),
        None,
        &context,
    )
    .expect("the class reads the row")
}

/// `reads_on` a row whose terms hold `stop_retained_variant`, that writes no FLAGS and whose
/// reference writes no HGVSp.
fn reads(row: Row, candidate: &str) -> bool {
    reads_on(row, KEPT_STOP_TERMS, "", "", candidate)
}

fn texts(items: &[&str]) -> Value {
    Value::List(
        items
            .iter()
            .map(|item| Value::Str(item.to_string()))
            .collect(),
    )
}

fn object(members: Vec<(&str, Value)>) -> Value {
    Value::Dict(
        members
            .into_iter()
            .map(|(name, member)| (name.to_string(), member))
            .collect::<IndexMap<String, Value>>(),
    )
}

/// The verdict of one comparison under the class; a green one has set the one row aside
/// under the class.
fn checked(report: &Value, applied: &vep_formats_compare::manifest::Documented) -> String {
    let verdict = verdict(report).to_string();
    assert!(verdict == "green" || verdict == "red");
    if verdict == "green" {
        assert_eq!(
            counts(&applied.applied),
            count_map(&[(STOP_CODON_CLASS, 1)])
        );
    }
    verdict
}

/// The verdicts of one made row under the class in the default, tab, VCF and JSON formats:
/// the row's HGVSp is the description `reference` on the one side and `candidate` on the
/// other (empty for no notation). Without the class the default-format row is a difference.
fn verdicts(row: Row, terms: &str, reference: &str, candidate: &str) -> BTreeSet<String> {
    let mut verdicts = BTreeSet::new();
    let sides = [hgvsp(reference, true), hgvsp(candidate, true)];

    let default = |hgvsp: &str| {
        let mut extra = String::from("IMPACT=MODERATE;STRAND=1");
        if !hgvsp.is_empty() {
            extra += &format!(";HGVSp={hgvsp}");
        }
        let cells = [
            "stop_row",
            "21:25004296-25004297",
            "CTGAGG",
            "ENSG1",
            "ENST00000000001",
            "Transcript",
            terms,
            "1181-1182",
            row.cds,
            row.protein,
            row.amino_acids,
            "-",
            "-",
            &extra,
        ];
        default_text(&[cells.join("\t")])
    };
    let (report, documented) = compare_default(
        &default(&sides[0]),
        &default(&sides[1]),
        &[STOP_CODON_CLASS],
    );
    verdicts.insert(checked(&report, &documented));
    assert_eq!(
        verdict_of(&default(&sides[0]), &default(&sides[1]), &[]),
        "red"
    );

    let tab_header = DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tHGVSp");
    let tab = |hgvsp: &str| {
        let cells = [
            "stop_row",
            "21:25004296-25004297",
            "CTGAGG",
            "ENSG1",
            "ENST00000000001",
            "Transcript",
            terms,
            "1181-1182",
            row.cds,
            row.protein,
            row.amino_acids,
            "-",
            "-",
            "MODERATE",
            "1",
            if hgvsp.is_empty() { "-" } else { hgvsp },
        ];
        format!("{tab_header}\n{}\n", cells.join("\t"))
    };
    let mut documented = documented_with(&[STOP_CODON_CLASS]);
    let report = compare(
        &parse_tab(&tab(&sides[0])).expect("the reference parses"),
        &parse_tab(&tab(&sides[1])).expect("the candidate parses"),
        &mut documented,
    );
    verdicts.insert(checked(&report, &documented));

    let inputs = vec![INPUT.to_string()];
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|CDS_position|Protein_position|Amino_acids|HGVSp\">\n";
    let vcf = |hgvsp: &str| {
        format!(
            "{csq_header}{}CSQ=CTGAGG|{}|ENST00000000001|Transcript|{}|{}|{}|{hgvsp}\n",
            INPUT.strip_suffix('.').expect("an empty INFO column"),
            terms.replace(',', "&"),
            row.cds,
            row.protein,
            row.amino_acids,
        )
    };
    let mut documented = documented_with(&[STOP_CODON_CLASS]);
    let report = compare(
        &parse_vcf(&vcf(&sides[0]), &inputs).expect("the reference parses"),
        &parse_vcf(&vcf(&sides[1]), &inputs).expect("the candidate parses"),
        &mut documented,
    );
    verdicts.insert(checked(&report, &documented));

    let held_terms: Vec<&str> = terms.split(',').collect();
    let bound = |cell: &str, last: bool| {
        let (first, second) = cell.split_once('-').unwrap_or((cell, cell));
        let text = if last { second } else { first };
        Value::Int(text.parse().expect("a position"))
    };
    let json = |description: &str| {
        let mut entry = vec![
            ("variant_allele", Value::Str("CTGAGG".to_string())),
            ("transcript_id", Value::Str("ENST00000000001".to_string())),
            ("consequence_terms", texts(&held_terms)),
            ("cds_start", bound(row.cds, false)),
            ("cds_end", bound(row.cds, true)),
            ("protein_start", bound(row.protein, false)),
            ("protein_end", bound(row.protein, true)),
            ("amino_acids", Value::Str(row.amino_acids.to_string())),
        ];
        if !description.is_empty() {
            entry.push(("hgvsp", Value::Str(hgvsp(description, false))));
        }
        let record = object(vec![
            ("input", Value::Str(INPUT.to_string())),
            ("id", Value::Str("stop_row".to_string())),
            ("allele_string", Value::Str("-/CTGAGG".to_string())),
            (
                "most_severe_consequence",
                Value::Str(held_terms[0].to_string()),
            ),
            ("transcript_consequences", Value::List(vec![object(entry)])),
        ]);
        pyjson::dumps(&record, &pyjson::DumpOptions::default()) + "\n"
    };
    let mut documented = documented_with(&[STOP_CODON_CLASS]);
    let report = compare(
        &parse_json(&json(reference), &inputs).expect("the reference parses"),
        &parse_json(&json(candidate), &inputs).expect("the candidate parses"),
        &mut documented,
    );
    verdicts.insert(checked(&report, &documented));
    verdicts
}

fn only(verdict: &str) -> BTreeSet<String> {
    BTreeSet::from([verdict.to_string()])
}

/// The rows of made-up insertions inside the last residue's codon that keep the protein, on
/// a forward-strand transcript whose stop is residue 366 and a reverse-strand one whose stop
/// is residue 109: a whole number of codons on either strand, five and eight bases on the
/// reverse strand, and six bases written with a base the two alleles share, whose coding
/// position is one base. Each is read under the class.
#[test]
fn the_stop_unchanged_against_no_notation_is_read_on_either_strand_at_every_length() {
    let forward = |amino_acids| Row {
        cds: "1094-1095",
        protein: "365",
        amino_acids,
    };
    let reverse = |amino_acids| Row {
        cds: "323-324",
        protein: "108",
        amino_acids,
    };
    for amino_acids in ["V/V*G", "V/V*GV"] {
        assert!(reads(forward(amino_acids), KEPT_STOP), "{amino_acids}");
    }
    for amino_acids in ["A/A*X", "A/A*G", "A/A*GX", "A/A*GV"] {
        assert!(reads(reverse(amino_acids), "Ter109="), "{amino_acids}");
    }
    let shared_base = Row {
        cds: "1095",
        protein: "365",
        amino_acids: "V/V*G",
    };
    assert!(reads(shared_base, KEPT_STOP));
    // An insertion between the codon's first and second base.
    let after_the_first_base = Row {
        cds: "1093-1094",
        protein: "365",
        amino_acids: "V/V*G",
    };
    assert!(reads(after_the_first_base, KEPT_STOP));
}

/// Through the comparison the row is set aside once in the default, tab, VCF and JSON
/// formats, and the candidate's value inside prediction parentheses is read the same way.
#[test]
fn the_stop_unchanged_against_no_notation_is_set_aside_once_a_row_in_every_format() {
    assert_eq!(
        verdicts(SIX_BASES_FORWARD, KEPT_STOP_TERMS, "", KEPT_STOP),
        only("green")
    );
    assert_eq!(
        verdicts(
            SIX_BASES_FORWARD,
            KEPT_STOP_TERMS,
            "",
            &format!("({KEPT_STOP})")
        ),
        only("green")
    );
}

/// The sides exchanged, the reference's stop unchanged against a candidate that writes no
/// HGVSp, stay a difference in every format.
#[test]
fn the_stop_unchanged_against_no_notation_stays_a_difference_with_the_sides_exchanged() {
    assert!(reads(SIX_BASES_FORWARD, KEPT_STOP));
    assert!(!reads_on(
        SIX_BASES_FORWARD,
        KEPT_STOP_TERMS,
        "",
        KEPT_STOP,
        ""
    ));
    assert_eq!(
        verdicts(SIX_BASES_FORWARD, KEPT_STOP_TERMS, KEPT_STOP, ""),
        only("red")
    );
}

/// A row off the shape stays a difference: terms without `stop_retained_variant`, FLAGS
/// that name `cds_end_NF`, a residue number that is not the one after the row's codon, a row
/// that names two codons or none, coding positions outside the codon, and amino acids that
/// do not read the codon's residue and then a stop.
#[test]
fn an_empty_reference_value_off_the_shape_stays_a_difference() {
    assert!(reads(SIX_BASES_FORWARD, KEPT_STOP));
    for terms in [
        "inframe_insertion",
        "inframe_insertion,coding_sequence_variant",
        "stop_gained,frameshift_variant",
    ] {
        assert!(
            !reads_on(SIX_BASES_FORWARD, terms, "", "", KEPT_STOP),
            "{terms}"
        );
        assert_eq!(
            verdicts(SIX_BASES_FORWARD, terms, "", KEPT_STOP),
            only("red"),
            "{terms}"
        );
    }
    assert!(!reads_on(
        SIX_BASES_FORWARD,
        KEPT_STOP_TERMS,
        "cds_end_NF",
        "",
        KEPT_STOP
    ));
    for candidate in ["Ter365=", "Ter367=", "Ter1=", "Ter0366="] {
        assert!(!reads(SIX_BASES_FORWARD, candidate), "{candidate}");
    }
    for protein in ["365-366", "364-365", "-", "?"] {
        let row = Row {
            protein,
            ..SIX_BASES_FORWARD
        };
        assert!(!reads(row, KEPT_STOP), "{protein}");
    }
    for cds in ["1092-1093", "1095-1096", "1096", "1092", "-", "?-1095"] {
        let row = Row {
            cds,
            ..SIX_BASES_FORWARD
        };
        assert!(!reads(row, KEPT_STOP), "{cds}");
    }
    for amino_acids in [
        "V/VG*", "V/A*G", "V/V", "V/VG", "-/*G", "*/*", "X/X*", "VA/V*G", "V", "-",
    ] {
        let row = Row {
            amino_acids,
            ..SIX_BASES_FORWARD
        };
        assert!(!reads(row, KEPT_STOP), "{amino_acids}");
        assert_eq!(
            verdicts(row, KEPT_STOP_TERMS, "", KEPT_STOP),
            only("red"),
            "{amino_acids}"
        );
    }
}

/// Against an empty reference value the candidate's value is the stop unchanged alone: an
/// extension named at the stop, a frameshift, an insertion and an unclosed parenthesis stay
/// differences.
#[test]
fn an_empty_reference_value_against_another_description_stays_a_difference() {
    assert!(reads(SIX_BASES_FORWARD, KEPT_STOP));
    for candidate in [
        "Ter366LysextTer3",
        "Val365AlafsTer2",
        "Val365_Ter366insLys",
        "Val365=",
        "(Ter366=",
    ] {
        assert!(!reads(SIX_BASES_FORWARD, candidate), "{candidate}");
        assert_eq!(
            verdicts(SIX_BASES_FORWARD, KEPT_STOP_TERMS, "", candidate),
            only("red"),
            "{candidate}"
        );
    }
}
