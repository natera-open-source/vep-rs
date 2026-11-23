// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The ninth shape of `hgvsp_stop_codon_protein_comparison`: the reference writes one residue
//! and a stop as inserted between two neighbouring residues, the residue the right flank's,
//! and the candidate writes the stop as unchanged at the residue after the right flank, on a
//! row whose terms hold `stop_retained_variant`. The class is asked as a comparison asks it,
//! its function is called where one side writes no notation, and two outputs are compared in
//! the default, tab, VCF and JSON formats.

mod compare_common;

use compare_common::*;
use indexmap::IndexMap;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use vep_formats_compare::error::Error;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::hgvsp_stop::stop_codon_protein_comparison;
use vep_formats_compare::rules::{field_class_accepts, RuleContext, Side};

const STOP_CODON_CLASS: &str = "hgvsp_stop_codon_protein_comparison";
const PROTEIN: &str = "ENSP00000000001.1";
const INPUT: &str = "21\t25004296\tstop_row\tAC\tA\t.\t.\t.";

/// The shape on a made-up protein whose stop is residue 284: the reference's description is
/// the residue at 283 and a stop inserted between residues 282 and 283, the candidate's the
/// stop unchanged at 284, and the row's terms hold `stop_retained_variant`.
const INSERTED_RESIDUE_AND_STOP: &str = "Asp282_Trp283insTrpTer";
const KEPT_STOP: &str = "Ter284=";
const KEPT_STOP_TERMS: &str = "inframe_insertion,stop_retained_variant,NMD_transcript_variant";

/// An HGVSp value: the protein, `:p.` and the description, empty for no description; the
/// default, tab and VCF formats escape the `=` of a description, the JSON format does not.
fn hgvsp(description: &str, protein: &str, escaped: bool) -> String {
    if description.is_empty() {
        return String::new();
    }
    let text = format!("{protein}:p.{description}");
    if escaped {
        text.replace('=', "%3D")
    } else {
        text
    }
}

/// What the class says of the differing HGVSp values `ev` and `av` of a row whose terms are
/// `terms` and whose FLAGS are `flags`, asked as a comparison without side tables asks it.
fn class_reads(terms: &str, flags: &str, ev: &str, av: &str) -> Result<bool, Error> {
    let held: &[(&str, &str)] = if flags.is_empty() {
        &[]
    } else {
        &[("FLAGS", flags)]
    };
    let e = entry(
        "21:25004297-25004298|-|ENST00000000001|Transcript",
        None,
        Some("21:25004297-25004298"),
        "-",
        "ENST00000000001",
        terms,
        held,
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
    field_class_accepts(STOP_CODON_CLASS, &e, "HGVSp", ev, av, None, &context)
}

/// The class on one protein's two descriptions after `p.`, the `=` escaped, on a row whose
/// terms are `terms` and whose FLAGS are `flags`.
fn reads_on(reference: &str, candidate: &str, terms: &str, flags: &str) -> bool {
    class_reads(
        terms,
        flags,
        &hgvsp(reference, PROTEIN, true),
        &hgvsp(candidate, PROTEIN, true),
    )
    .expect("the class reads the row")
}

/// `reads_on` a row whose terms hold `stop_retained_variant` and that writes no FLAGS.
fn reads(reference: &str, candidate: &str) -> bool {
    reads_on(reference, candidate, KEPT_STOP_TERMS, "")
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
fn checked(report: &Value, documented: &Documented) -> String {
    let verdict = verdict(report).to_string();
    assert!(verdict == "green" || verdict == "red");
    if verdict == "green" {
        assert_eq!(
            counts(&documented.applied),
            count_map(&[(STOP_CODON_CLASS, 1)])
        );
    }
    verdict
}

/// The verdicts of one made row under the class in the default, tab, VCF and JSON formats:
/// the row's HGVSp is the description `reference` after `p.` on the one side and `candidate`
/// on the other (empty for no notation), its terms `consequence` and its FLAGS `flags`.
/// Without the class the default-format row is a difference.
fn verdicts_under(
    reference: &str,
    candidate: &str,
    consequence: &str,
    flags: &str,
    candidate_protein: &str,
) -> BTreeSet<String> {
    let mut verdicts = BTreeSet::new();
    let sides = [
        hgvsp(reference, PROTEIN, true),
        hgvsp(candidate, candidate_protein, true),
    ];

    let default = |hgvsp: &str| {
        let mut extra = String::from("IMPACT=HIGH;STRAND=1");
        if !flags.is_empty() {
            extra += &format!(";FLAGS={flags}");
        }
        extra += ";HGVSc=ENST00000000001.1:c.850del";
        if !hgvsp.is_empty() {
            extra += &format!(";HGVSp={hgvsp}");
        }
        default_text(&[default_row(
            "stop_row",
            "21:25004297",
            "-",
            "ENST00000000001",
            consequence,
            &extra,
        )])
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

    let tab_header = DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tFLAGS\tHGVSc\tHGVSp");
    let tab = |hgvsp: &str| {
        let cells = [
            "stop_row",
            "21:25004297",
            "-",
            "ENSG1",
            "ENST00000000001",
            "Transcript",
            consequence,
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
            "HIGH",
            "1",
            if flags.is_empty() { "-" } else { flags },
            "ENST00000000001.1:c.850del",
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
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|FLAGS|HGVSc|HGVSp\">\n";
    let vcf = |hgvsp: &str| {
        format!(
            "{csq_header}{}CSQ=-|{}|ENST00000000001|Transcript|{}|ENST00000000001.1:c.850del|{hgvsp}\n",
            INPUT.strip_suffix('.').expect("an empty INFO column"),
            consequence.replace(',', "&"),
            flags.replace(',', "&"),
        )
    };
    let mut documented = documented_with(&[STOP_CODON_CLASS]);
    let report = compare(
        &parse_vcf(&vcf(&sides[0]), &inputs).expect("the reference parses"),
        &parse_vcf(&vcf(&sides[1]), &inputs).expect("the candidate parses"),
        &mut documented,
    );
    verdicts.insert(checked(&report, &documented));

    let terms: Vec<&str> = consequence.split(',').collect();
    let json = |description: &str, protein: &str| {
        let mut entry = vec![
            ("variant_allele", Value::Str("-".to_string())),
            ("transcript_id", Value::Str("ENST00000000001".to_string())),
            ("consequence_terms", texts(&terms)),
            (
                "hgvsc",
                Value::Str("ENST00000000001.1:c.850del".to_string()),
            ),
        ];
        if !flags.is_empty() {
            entry.push(("flags", texts(&flags.split(',').collect::<Vec<_>>())));
        }
        if !description.is_empty() {
            entry.push(("hgvsp", Value::Str(hgvsp(description, protein, false))));
        }
        let record = object(vec![
            ("input", Value::Str(INPUT.to_string())),
            ("id", Value::Str("stop_row".to_string())),
            ("allele_string", Value::Str("C/-".to_string())),
            ("most_severe_consequence", Value::Str(terms[0].to_string())),
            ("transcript_consequences", Value::List(vec![object(entry)])),
        ]);
        pyjson::dumps(&record, &pyjson::DumpOptions::default()) + "\n"
    };
    let mut documented = documented_with(&[STOP_CODON_CLASS]);
    let report = compare(
        &parse_json(&json(reference, PROTEIN), &inputs).expect("the reference parses"),
        &parse_json(&json(candidate, candidate_protein), &inputs).expect("the candidate parses"),
        &mut documented,
    );
    verdicts.insert(checked(&report, &documented));
    verdicts
}

/// The verdicts of a row under the class, its candidate naming the reference's protein.
fn stop_codon_verdicts(
    reference: &str,
    candidate: &str,
    consequence: &str,
    flags: &str,
) -> BTreeSet<String> {
    verdicts_under(reference, candidate, consequence, flags, PROTEIN)
}

fn only(verdict: &str) -> BTreeSet<String> {
    BTreeSet::from([verdict.to_string()])
}

/// The classes a column's differences were set aside under, each with its count.
fn masked(report: &Value, column: &str) -> BTreeMap<String, i64> {
    at(report, &["columns", column, "masked"])
        .as_dict()
        .expect("the masked counts are an object")
        .iter()
        .map(|(klass, count)| match count {
            Value::Int(count) => (klass.clone(), *count),
            other => panic!("{klass} is counted {other:?}"),
        })
        .collect()
}

/// In the default format, a row of `protein` whose reference writes the description
/// `reference` and whose candidate writes `candidate` is set aside once, the HGVSp and the
/// Extra text each counting once under the class; without the class the row is one differing
/// HGVSp.
fn set_aside_once_in_the_default_format(protein: &str, reference: &str, candidate: &str) {
    let text = |description: &str| {
        let extra = format!(
            "IMPACT=MODERATE;STRAND=-1;HGVSp={}",
            hgvsp(description, protein, true)
        );
        default_text(&[default_row(
            "stop_row",
            "21:25004297-25004298",
            "-",
            "ENST00000000001",
            KEPT_STOP_TERMS,
            &extra,
        )])
    };
    let (reference, candidate) = (text(reference), text(candidate));
    let (report, documented) = compare_default(&reference, &candidate, &[STOP_CODON_CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(STOP_CODON_CLASS, 1)])
    );
    for column in ["HGVSp", "Extra"] {
        assert_eq!(
            masked(&report, column),
            count_map(&[(STOP_CODON_CLASS, 1)]),
            "{column}"
        );
        assert_eq!(
            number(&report, &["columns", column, "differing"]),
            0,
            "{column}"
        );
    }
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    assert!(reasons(&without).contains(&"field HGVSp: 1 differing entries"));
}

/// The class reads one residue and a stop inserted between two neighbouring residues, the
/// residue the right flank's, against the stop unchanged at the residue after the right
/// flank, on a row whose terms hold `stop_retained_variant`: beside other terms or alone,
/// the `=` escaped or not, at the protein's first residues, and with both descriptions
/// inside one pair of prediction parentheses.
#[test]
fn a_residue_and_a_stop_inserted_before_the_last_residue_is_read_against_the_stop_unchanged() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    assert!(reads_on(
        INSERTED_RESIDUE_AND_STOP,
        KEPT_STOP,
        "stop_retained_variant",
        ""
    ));
    assert_eq!(
        class_reads(
            KEPT_STOP_TERMS,
            "",
            &hgvsp(INSERTED_RESIDUE_AND_STOP, PROTEIN, false),
            &hgvsp(KEPT_STOP, PROTEIN, false),
        ),
        Ok(true)
    );
    assert!(reads("Met1_Trp2insTrpTer", "Ter3="));
    assert!(reads(
        &format!("({INSERTED_RESIDUE_AND_STOP})"),
        &format!("({KEPT_STOP})")
    ));
}

/// Through the comparison, one residue and a stop inserted between two neighbouring residues
/// against the candidate's stop unchanged at the residue after the right flank is set aside
/// once per row in the default, tab, VCF and JSON formats; in the default format the HGVSp
/// and the Extra text each count once under the class. Without the class the row is one
/// differing HGVSp.
#[test]
fn a_residue_and_a_stop_inserted_before_the_last_residue_is_set_aside_once_a_row() {
    assert_eq!(
        stop_codon_verdicts(INSERTED_RESIDUE_AND_STOP, KEPT_STOP, KEPT_STOP_TERMS, ""),
        only("green")
    );
    set_aside_once_in_the_default_format(PROTEIN, INSERTED_RESIDUE_AND_STOP, KEPT_STOP);
}

/// The stop unchanged in the reference against the candidate's inserted residue and stop
/// stays a difference in every format.
#[test]
fn an_inserted_residue_and_stop_stays_a_difference_with_the_sides_exchanged() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    assert!(!reads(KEPT_STOP, INSERTED_RESIDUE_AND_STOP));
    assert_eq!(
        stop_codon_verdicts(KEPT_STOP, INSERTED_RESIDUE_AND_STOP, KEPT_STOP_TERMS, ""),
        only("red")
    );
}

/// An inserted residue that is not the right flank's stays a difference: the protein's
/// residue before the inserted stop is then another than the reference protein's.
#[test]
fn an_inserted_residue_other_than_the_right_flanks_stays_a_difference() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    for reference in ["Asp282_Trp283insAspTer", "Asp282_Trp283insGlyTer"] {
        assert!(!reads(reference, KEPT_STOP), "{reference}");
        assert_eq!(
            stop_codon_verdicts(reference, KEPT_STOP, KEPT_STOP_TERMS, ""),
            only("red"),
            "{reference}"
        );
    }
}

/// Flanks that are not two neighbouring residues stay a difference: residue numbers more
/// than one apart, equal or falling, and a flank that is a stop.
#[test]
fn an_inserted_residue_and_stop_stays_a_difference_between_flanks_that_are_not_neighbours() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    for reference in [
        "Asp281_Trp283insTrpTer",
        "Asp283_Trp283insTrpTer",
        "Asp284_Trp283insTrpTer",
        "Ter282_Trp283insTrpTer",
        "Asp282_Ter283insTerTer",
    ] {
        assert!(!reads(reference, KEPT_STOP), "{reference}");
        assert_eq!(
            stop_codon_verdicts(reference, KEPT_STOP, KEPT_STOP_TERMS, ""),
            only("red"),
            "{reference}"
        );
    }
}

/// The candidate's stop unchanged at a residue other than the one after the right flank
/// stays a difference: the right flank's own number, the one before it and two past it.
#[test]
fn the_stop_unchanged_stays_a_difference_at_another_residue_than_the_one_after_the_right_flank() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    for candidate in ["Ter283=", "Ter282=", "Ter285="] {
        assert!(!reads(INSERTED_RESIDUE_AND_STOP, candidate), "{candidate}");
        assert_eq!(
            stop_codon_verdicts(INSERTED_RESIDUE_AND_STOP, candidate, KEPT_STOP_TERMS, ""),
            only("red"),
            "{candidate}"
        );
    }
}

/// More than one residue before the inserted stop stays a difference, whatever residue the
/// candidate names the stop at: the row does not name the reference residues past the right
/// flank. So do an inserted stop alone, an inserted residue with no stop and two stops.
#[test]
fn more_than_one_residue_before_the_inserted_stop_stays_a_difference() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    for (reference, candidate) in [
        ("Asp282_Trp283insTrpGlyTer", "Ter285="),
        ("Asp282_Trp283insTrpGlyTer", "Ter284="),
        ("Asp282_Trp283insTrpTrpTer", "Ter285="),
        ("Asp282_Trp283insTer", "Ter283="),
        ("Asp282_Trp283insTer", "Ter284="),
        ("Asp282_Trp283insTrp", "Ter284="),
        ("Asp282_Trp283insTrpTerTer", "Ter284="),
    ] {
        assert!(
            !reads(reference, candidate),
            "{reference} against {candidate}"
        );
        assert_eq!(
            stop_codon_verdicts(reference, candidate, KEPT_STOP_TERMS, ""),
            only("red"),
            "{reference} against {candidate}"
        );
    }
}

/// An inserted residue and stop on a row whose terms hold no `stop_retained_variant` stays a
/// difference: the insertion's term alone, a gained stop and a lost stop.
#[test]
fn an_inserted_residue_and_stop_stays_a_difference_on_a_row_without_the_retained_stop_term() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    for terms in [
        "inframe_insertion",
        "inframe_insertion,stop_gained",
        "stop_lost,inframe_insertion",
        "inframe_insertion,NMD_transcript_variant",
    ] {
        assert!(
            !reads_on(INSERTED_RESIDUE_AND_STOP, KEPT_STOP, terms, ""),
            "{terms}"
        );
        assert_eq!(
            stop_codon_verdicts(INSERTED_RESIDUE_AND_STOP, KEPT_STOP, terms, ""),
            only("red"),
            "{terms}"
        );
    }
}

/// An inserted residue and stop on a row whose FLAGS name `cds_end_NF` stays a difference,
/// with or without the start flag beside it.
#[test]
fn an_inserted_residue_and_stop_stays_a_difference_on_a_coding_sequence_no_stop_codon_closes() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    for flags in ["cds_end_NF", "cds_start_NF,cds_end_NF"] {
        assert!(
            !reads_on(INSERTED_RESIDUE_AND_STOP, KEPT_STOP, KEPT_STOP_TERMS, flags),
            "{flags}"
        );
        assert_eq!(
            stop_codon_verdicts(INSERTED_RESIDUE_AND_STOP, KEPT_STOP, KEPT_STOP_TERMS, flags),
            only("red"),
            "{flags}"
        );
    }
}

/// The two descriptions on two proteins stay a difference.
#[test]
fn an_inserted_residue_and_stop_stays_a_difference_on_another_protein() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    assert_eq!(
        class_reads(
            KEPT_STOP_TERMS,
            "",
            &hgvsp(INSERTED_RESIDUE_AND_STOP, PROTEIN, false),
            "ENSP00000000002.1:p.Ter284%3D",
        ),
        Ok(false)
    );
    assert_eq!(
        verdicts_under(
            INSERTED_RESIDUE_AND_STOP,
            KEPT_STOP,
            KEPT_STOP_TERMS,
            "",
            "ENSP00000000002.1"
        ),
        only("red")
    );
}

/// Residue numbers longer than the comparison converts end the comparison of an inserted
/// residue and stop with the comparator's error, as every malformed number does.
#[test]
fn an_inserted_residue_and_stop_at_a_residue_number_past_conversion_ends_the_comparison() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    let number = "7".repeat(5000);
    assert_eq!(
        class_reads(
            KEPT_STOP_TERMS,
            "",
            &hgvsp(&format!("Asp{number}_Trp{number}insTrpTer"), PROTEIN, true),
            &hgvsp(&format!("Ter{number}="), PROTEIN, true),
        ),
        Err(Error::Comparator(format!(
            "HGVSp residue number '{}...' is not an integer",
            "7".repeat(24)
        )))
    );
}

/// Through the comparison, an inserted residue and stop against the stop unchanged, each
/// description inside a pair of prediction parentheses, is set aside once per row in the
/// default, tab, VCF and JSON formats. With the protein named without its version, the
/// function reads the pair, and in the default format the HGVSp and the Extra text each count
/// once under the class. Without the class the row is one differing HGVSp.
#[test]
fn an_inserted_residue_and_stop_inside_prediction_parentheses_is_set_aside_once_a_row() {
    let predicted = format!("({INSERTED_RESIDUE_AND_STOP})");
    let kept = format!("({KEPT_STOP})");
    assert_eq!(
        stop_codon_verdicts(&predicted, &kept, KEPT_STOP_TERMS, ""),
        only("green")
    );
    let protein = "ENSP00000000001";
    assert_eq!(
        stop_codon_protein_comparison(
            KEPT_STOP_TERMS,
            "",
            &hgvsp(&predicted, protein, false),
            &hgvsp(&kept, protein, true),
        ),
        Ok(true)
    );
    set_aside_once_in_the_default_format(protein, &predicted, &kept);
}

/// One description inside prediction parentheses against the other outside them stays a
/// difference in every format, whichever side has the parentheses; so do two descriptions
/// after a parenthesis that nothing closes.
#[test]
fn an_inserted_residue_and_stop_stays_a_difference_with_one_value_alone_in_parentheses() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    for (reference, candidate) in [
        (
            format!("({INSERTED_RESIDUE_AND_STOP})"),
            KEPT_STOP.to_string(),
        ),
        (
            INSERTED_RESIDUE_AND_STOP.to_string(),
            format!("({KEPT_STOP})"),
        ),
        (
            format!("({INSERTED_RESIDUE_AND_STOP}"),
            format!("({KEPT_STOP}"),
        ),
    ] {
        assert!(
            !reads(&reference, &candidate),
            "{reference} against {candidate}"
        );
        assert_eq!(
            stop_codon_verdicts(&reference, &candidate, KEPT_STOP_TERMS, ""),
            only("red"),
            "{reference} against {candidate}"
        );
    }
}

/// An inserted residue and stop in the reference against a candidate that writes no HGVSp
/// stays a difference in every format; so does the candidate's inserted residue and stop
/// against a reference that writes none.
#[test]
fn an_inserted_residue_and_stop_stays_a_difference_against_no_notation() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    let inserted = hgvsp(INSERTED_RESIDUE_AND_STOP, PROTEIN, false);
    for none in ["", "-"] {
        assert_eq!(
            stop_codon_protein_comparison(KEPT_STOP_TERMS, "", &inserted, none),
            Ok(false),
            "{none:?}"
        );
        assert_eq!(
            stop_codon_protein_comparison(KEPT_STOP_TERMS, "", none, &inserted),
            Ok(false),
            "{none:?}"
        );
    }
    assert_eq!(
        stop_codon_verdicts(INSERTED_RESIDUE_AND_STOP, "", KEPT_STOP_TERMS, ""),
        only("red")
    );
    assert_eq!(
        stop_codon_verdicts("", INSERTED_RESIDUE_AND_STOP, KEPT_STOP_TERMS, ""),
        only("red")
    );
}

/// An inserted residue and stop in the reference stays a difference against a candidate
/// that writes anything but the stop unchanged: the stop replaced by a stop, an extension
/// named at the stop, the stop deleted, a frameshift or a duplication at the right flank, the
/// right flank unchanged, and the same insertion without its stop.
#[test]
fn an_inserted_residue_and_stop_stays_a_difference_against_anything_but_the_stop_unchanged() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    for candidate in [
        "Ter284delinsTer",
        "Ter284TrpextTer3",
        "Ter284del",
        "Trp283GlyfsTer2",
        "Trp283dup",
        "Trp283=",
        "Asp282_Trp283insTrp",
    ] {
        assert!(!reads(INSERTED_RESIDUE_AND_STOP, candidate), "{candidate}");
        assert_eq!(
            stop_codon_verdicts(INSERTED_RESIDUE_AND_STOP, candidate, KEPT_STOP_TERMS, ""),
            only("red"),
            "{candidate}"
        );
    }
}

/// The stop unchanged in the reference against a candidate that writes no HGVSp stays a
/// difference in every format, outside and inside prediction parentheses, on a row whose
/// terms hold `stop_retained_variant` and on one whose terms hold a gained stop and a
/// frameshift: the class reads two descriptions of one protein, never one against none.
#[test]
fn the_stop_unchanged_stays_a_difference_against_no_notation() {
    assert!(reads(INSERTED_RESIDUE_AND_STOP, KEPT_STOP));
    for terms in [KEPT_STOP_TERMS, "stop_gained,frameshift_variant"] {
        for reference in [KEPT_STOP.to_string(), format!("({KEPT_STOP})")] {
            let unchanged = hgvsp(&reference, PROTEIN, true);
            for none in ["", "-"] {
                assert_eq!(
                    stop_codon_protein_comparison(terms, "", &unchanged, none),
                    Ok(false),
                    "{terms}, {reference}, {none:?}"
                );
            }
            assert_eq!(
                stop_codon_verdicts(&reference, "", terms, ""),
                only("red"),
                "{terms}, {reference}"
            );
        }
    }
}
