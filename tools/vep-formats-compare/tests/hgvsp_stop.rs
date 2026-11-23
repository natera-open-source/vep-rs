// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HGVSp stop-codon classes through the comparison of two outputs: what each class sets
//! aside, once a row, in the default, tab, VCF and JSON formats, and the nearest shapes that
//! stay differences.

mod compare_common;

use compare_common::*;
use indexmap::IndexMap;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::hgvsp_stop::stop_named_on_incomplete_coding_sequence;
use vep_formats_compare::shared::ABSENT;

const STOP_CODON_CLASS: &str = "hgvsp_stop_codon_protein_comparison";
const INCOMPLETE_CLASS: &str = "hgvsp_stop_named_on_incomplete_coding_sequence";
const PROTEIN: &str = "ENSP00000000001.1";
const INPUT: &str = "21\t25004296\tstop_row\tAC\tA\t.\t.\t.";

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

/// The verdict of one comparison under `klass`; a green one has set the one row aside
/// under the class.
fn checked(report: &Value, documented: &Documented, klass: &str) -> String {
    let verdict = verdict(report).to_string();
    assert!(verdict == "green" || verdict == "red");
    if verdict == "green" {
        assert_eq!(counts(&documented.applied), count_map(&[(klass, 1)]));
    }
    verdict
}

/// The verdicts of one synthetic row under the class `klass` in the default, tab, VCF and
/// JSON formats: the row's HGVSp is the description `reference` after `p.` on the one side
/// and `candidate` on the other (empty for no notation), its terms `consequence` and its
/// FLAGS `flags`. Without the class the default-format row is a difference.
fn verdicts_under(
    klass: &str,
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
    let (report, documented) = compare_default(&default(&sides[0]), &default(&sides[1]), &[klass]);
    verdicts.insert(checked(&report, &documented, klass));
    assert_eq!(
        verdict_of(&default(&sides[0]), &default(&sides[1]), &[]),
        "red"
    );

    let tab_header = DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tFLAGS\tHGVSc\tHGVSp");
    let or_dash = |text: &str| {
        if text.is_empty() {
            "-".to_string()
        } else {
            text.to_string()
        }
    };
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
            &or_dash(flags),
            "ENST00000000001.1:c.850del",
            &or_dash(hgvsp),
        ];
        format!("{tab_header}\n{}\n", cells.join("\t"))
    };
    let mut documented = documented_with(&[klass]);
    let report = compare(
        &parse_tab(&tab(&sides[0])).expect("the reference parses"),
        &parse_tab(&tab(&sides[1])).expect("the candidate parses"),
        &mut documented,
    );
    verdicts.insert(checked(&report, &documented, klass));

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
    let mut documented = documented_with(&[klass]);
    let report = compare(
        &parse_vcf(&vcf(&sides[0]), &inputs).expect("the reference parses"),
        &parse_vcf(&vcf(&sides[1]), &inputs).expect("the candidate parses"),
        &mut documented,
    );
    verdicts.insert(checked(&report, &documented, klass));

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
    let mut documented = documented_with(&[klass]);
    let report = compare(
        &parse_json(&json(reference, PROTEIN), &inputs).expect("the reference parses"),
        &parse_json(&json(candidate, candidate_protein), &inputs).expect("the candidate parses"),
        &mut documented,
    );
    verdicts.insert(checked(&report, &documented, klass));
    verdicts
}

/// The verdicts of a row under the stop-codon class, its candidate naming the reference's
/// protein.
fn stop_codon_verdicts(
    reference: &str,
    candidate: &str,
    consequence: &str,
    flags: &str,
) -> BTreeSet<String> {
    verdicts_under(
        STOP_CODON_CLASS,
        reference,
        candidate,
        consequence,
        flags,
        PROTEIN,
    )
}

fn only(verdict: &str) -> BTreeSet<String> {
    BTreeSet::from([verdict.to_string()])
}

/// Each row of `rows`, a reference description, a candidate description and the row's
/// terms, is set aside in every format.
fn set_aside(rows: &[(&str, &str, &str)]) {
    for (reference, candidate, terms) in rows {
        assert_eq!(
            stop_codon_verdicts(reference, candidate, terms, ""),
            only("green"),
            "{reference} against {candidate}"
        );
    }
}

/// The stop replaced by a stop in the reference, where the protein does not change and the
/// candidate writes the stop as unchanged at the same residue, is set aside once per row in
/// the default, tab, VCF and JSON formats, and is a difference without the class.
#[test]
fn the_stop_replaced_by_a_stop_is_set_aside_against_the_stop_unchanged() {
    set_aside(&[
        (
            "Ter284delinsTer",
            "Ter284=",
            "inframe_insertion,stop_retained_variant",
        ),
        ("Ter291delinsTer", "Ter291=", "stop_retained_variant"),
    ]);
}

/// The stop alone replaced by the first residues that follow the protein, closed by `Ter`
/// or not, against the candidate's extension named at the stop with the first of them: a
/// count above the residues the reference names, an unknown count, and the number of the
/// residues named where no `Ter` closes them.
#[test]
fn the_stop_replaced_by_the_first_residues_is_set_aside_against_an_extension() {
    set_aside(&[
        ("Ter284delinsIleTer", "Ter284IleextTer22", "stop_lost"),
        ("Ter284delinsIleTer", "Ter284IleextTer?", "stop_lost"),
        (
            "Ter284delinsIleLeuValProTer",
            "Ter284IleextTer?",
            "stop_lost",
        ),
        (
            "Ter284delinsIleLeuValPro",
            "Ter284IleextTer9",
            "stop_lost,inframe_insertion",
        ),
        ("Ter284delinsIleLeuValPro", "Ter284IleextTer4", "stop_lost"),
    ]);
}

/// The stop replaced by every added residue and then the stop the inserted bases bring,
/// against the candidate's extension whose count is the number of residues named.
#[test]
fn the_stop_replaced_through_an_inserted_stop_is_set_aside_against_an_extension() {
    set_aside(&[
        (
            "Ter284delinsAsxSerGlxTrpArgPheGlyAspValHisTer",
            "Ter284AsxextTer10",
            "stop_lost,inframe_insertion",
        ),
        ("Ter284delinsAsxTer", "Ter284AsxextTer1", "stop_lost"),
    ]);
}

/// Residues inserted between the last residue and the stop against the candidate's extension
/// named at the stop with the first of them, its count their number.
#[test]
fn residues_inserted_before_the_stop_are_set_aside_against_an_extension() {
    set_aside(&[
        (
            "Asp283_Ter284insHis",
            "Ter284HisextTer1",
            "inframe_insertion",
        ),
        (
            "Asp283_Ter284insHisHisAsxGly",
            "Ter284HisextTer4",
            "inframe_insertion",
        ),
    ]);
}

/// The stop written as unchanged on a row whose terms hold `stop_lost` against the
/// candidate's extension named at the stop, with a count or an unknown one.
#[test]
fn a_lost_stop_written_as_unchanged_is_set_aside_against_an_extension() {
    set_aside(&[
        ("Ter284=", "Ter284LysextTer47", "stop_lost"),
        (
            "Ter284=",
            "Ter284LysextTer?",
            "stop_lost,NMD_transcript_variant",
        ),
    ]);
}

/// A range from a residue to the stop replaced by residues, none of them a stop, on a row
/// whose terms hold `stop_lost`, against the candidate's frameshift named at that residue
/// with the first of them: an unknown count, and a count above the residues the reference
/// names that leaves the new stop away from the stop's own residue.
#[test]
fn a_range_to_the_stop_replaced_is_set_aside_against_a_frameshift() {
    set_aside(&[
        (
            "Trp283_Ter284delinsAsx",
            "Trp283AsxfsTer4",
            "stop_lost,protein_altering_variant,NMD_transcript_variant",
        ),
        (
            "Trp283_Ter284delinsAsx",
            "Trp283AsxfsTer?",
            "stop_lost,protein_altering_variant",
        ),
        (
            "Trp281_Ter284delinsAsxGlyHis",
            "Trp281AsxfsTer24",
            "stop_lost",
        ),
        (
            "Trp281_Ter284delinsAsxGlyHis",
            "Trp281AsxfsTer9",
            "stop_lost,inframe_deletion",
        ),
    ]);
}

/// An extension with an unknown count on a row whose terms hold `stop_lost` against the
/// candidate's same extension with a count.
#[test]
fn an_unknown_extension_count_is_set_aside_against_a_count() {
    set_aside(&[
        ("Ter557AlaextTer?", "Ter557AlaextTer35", "stop_lost"),
        (
            "Ter591AlaextTer?",
            "Ter591AlaextTer1",
            "stop_lost,NMD_transcript_variant",
        ),
    ]);
}

/// The stop written as deleted with an extension count and no residue in the stop's place,
/// on a row whose terms hold `stop_lost`, against the candidate's extension at the same
/// residue with the same count.
#[test]
fn the_stop_deleted_with_an_extension_count_is_set_aside_against_an_extension() {
    set_aside(&[
        (
            "Ter284delextTer12",
            "Ter284LysextTer12",
            "stop_lost,inframe_deletion",
        ),
        (
            "Ter230delextTer11",
            "Ter230ArgextTer11",
            "stop_lost,NMD_transcript_variant",
        ),
        (
            "Ter220delextTer?",
            "Ter220ArgextTer?",
            "stop_lost,inframe_deletion",
        ),
    ]);
}

/// The nearest shapes stay differences in every format: each shape with the sides exchanged,
/// another residue number or first residue, a count the reference's residues do not bear out,
/// an insertion not beside the stop, the stop unchanged on a row whose terms hold no
/// `stop_lost`, a row whose FLAGS name `cds_end_NF`, a stop named at residue 1, a reference
/// without a notation, the stop written as deleted with no extension count, with another count
/// or residue number than the extension's or on a row without `stop_lost`, a range that ends
/// at the stop against a frameshift named at another residue, with another first residue, with
/// a count the reference's residues do not bear out or one that leaves the new stop at the
/// stop's own residue, on a row whose terms hold no `stop_lost`, or closed by a stop or
/// opening at the stop's own residue, two extensions that differ, an unknown count against a
/// count with the sides exchanged, with another first residue or residue number or on a row
/// whose terms hold no `stop_lost`, and the same notation on another protein.
#[test]
fn the_nearest_shapes_of_the_stop_codon_class_stay_differences() {
    #[rustfmt::skip]
    let rows: &[(&str, &str, &str, &str, &str)] = &[
        ("Ter284=", "Ter284delinsTer", "inframe_insertion,stop_retained_variant", "", "the sides exchanged"),
        ("Ter284IleextTer22", "Ter284delinsIleTer", "stop_lost", "", "the sides exchanged"),
        ("Ter284AsxextTer10", "Ter284delinsAsxSerGlxTrpArgPheGlyAspValHisTer", "stop_lost,inframe_insertion", "", "the sides exchanged"),
        ("Ter284HisextTer1", "Asp283_Ter284insHis", "inframe_insertion", "", "the sides exchanged"),
        ("Ter284LysextTer47", "Ter284=", "stop_lost", "", "the sides exchanged"),
        ("Ter284delinsTer", "Ter285=", "stop_retained_variant", "", "another residue number"),
        ("Ter284delinsIleTer", "Ter284=", "stop_lost", "", "the stop unchanged against added residues"),
        ("Ter284delinsIleTer", "Ter284LeuextTer22", "stop_lost", "", "another first residue"),
        ("Ter284delinsIleTer", "Ter285IleextTer22", "stop_lost", "", "another residue number"),
        ("Ter284delinsIleLeuValProTer", "Ter284IleextTer3", "stop_lost", "", "a count below the residues the reference names"),
        ("Ter284delinsAsxSerGlxTrpArgPheGlyAspValHisTer", "Ter284AsxextTer9", "stop_lost", "", "a count below the residues the reference names"),
        ("Asp283_Ter284insHisHisAsxGly", "Ter284HisextTer5", "inframe_insertion", "", "a count that is not the number of inserted residues"),
        ("Asp283_Ter284insHisHisAsxGly", "Ter284GlnextTer4", "inframe_insertion", "", "another first residue"),
        ("Asp282_Ter284insHisHisAsxGly", "Ter284HisextTer4", "inframe_insertion", "", "an insertion not beside the stop"),
        ("Ter284=", "Ter284LysextTer47", "splice_acceptor_variant,stop_retained_variant", "", "the stop unchanged on a row whose terms hold no stop_lost"),
        ("Ter284=", "Ter284LysextTer?", "stop_lost", "cds_end_NF", "a coding sequence that no stop codon closes"),
        ("Ter284=", "Ter284LysextTer?", "stop_lost", "cds_start_NF,cds_end_NF", "a coding sequence that no stop codon closes"),
        ("Ter1=", "Ter1LysextTer5", "stop_lost", "cds_start_NF", "a stop named at residue 1 of an incomplete first codon"),
        ("", "Ter284LysextTer12", "stop_lost", "", "no protein notation in the reference"),
        ("Asp283_Ter284delinsPhePro", "Asp283ProfsTer9", "stop_lost", "", "a frameshift whose new residue is not the first the reference names"),
        ("Ter284LysextTer12", "Ter284LysextTer13", "stop_lost", "", "two extensions that differ"),
        ("Ter284delinsIleTerTer", "Ter284IleextTer2", "stop_lost", "", "a stop inside the named residues"),
        ("Trp283AsxfsTer4", "Trp283_Ter284delinsAsx", "stop_lost", "", "the sides exchanged"),
        ("Trp283_Ter284delinsAsx", "Trp282AsxfsTer4", "stop_lost", "", "a frameshift named at another residue number"),
        ("Trp283_Ter284delinsAsx", "Pro283AsxfsTer4", "stop_lost", "", "a frameshift named at another residue"),
        ("Trp281_Ter284delinsAsxGlyHis", "Trp281GlyfsTer24", "stop_lost", "", "a frameshift whose new residue is not the first the reference names"),
        ("Trp281_Ter284delinsAsxGlyHis", "Trp281AsxfsTer3", "stop_lost", "", "a count the reference's residues do not bear out"),
        ("Trp283_Ter284delinsAsx", "Trp283AsxfsTer4", "inframe_deletion,stop_retained_variant", "", "a range ending at the stop on a row whose terms hold no stop_lost"),
        ("Trp283_Ter284delinsAsxTer", "Trp283AsxfsTer2", "stop_lost", "", "a stop after the residues the reference names"),
        ("Trp283_Ter284delinsAsx", "Trp283AsxfsTer4", "stop_lost", "cds_end_NF", "a coding sequence that no stop codon closes"),
        ("Trp284_Ter284delinsAsx", "Trp284AsxfsTer4", "stop_lost", "", "a range that opens at the stop's own residue"),
        ("Trp281_Ter284delinsAsxGlyHis", "Trp281AsxfsTer4", "stop_lost,inframe_deletion", "", "a count that leaves the new stop at the stop's own residue"),
        ("Trp283_Ter284delinsAsx", "Trp283AsxfsTer2", "stop_lost", "", "a count that leaves the new stop at the stop's own residue"),
        ("Ter557AlaextTer35", "Ter557AlaextTer?", "stop_lost", "", "a count against an unknown count"),
        ("Ter557AlaextTer?", "Ter557GlxextTer35", "stop_lost", "", "an unknown count against a count with another first residue"),
        ("Ter557AlaextTer?", "Ter558AlaextTer35", "stop_lost", "", "an unknown count against a count at another residue number"),
        ("Ter557AlaextTer?", "Ter557AlaextTer35", "frameshift_variant", "", "an unknown count against a count on a row whose terms hold no stop_lost"),
        ("Ter557AlaextTer?", "Ter557AlaextTer35", "stop_lost", "cds_end_NF", "an unknown count against a count on a coding sequence that no stop codon closes"),
        ("Ter230ArgextTer11", "Ter230delextTer11", "stop_lost", "", "an extension against the stop written as deleted with its count"),
        ("Ter230delextTer11", "Ter230ArgextTer12", "stop_lost", "", "the stop written as deleted with another count than the extension's"),
        ("Ter230delextTer?", "Ter230ArgextTer11", "stop_lost", "", "the stop written as deleted with an unknown count against a count"),
        ("Ter230delextTer11", "Ter231ArgextTer11", "stop_lost", "", "the stop written as deleted at another residue number"),
        ("Ter230delextTer11", "Ter230ArgextTer11", "inframe_deletion", "", "the stop written as deleted with an extension count on a row whose terms hold no stop_lost"),
        ("Ter230del", "Ter230ArgextTer11", "stop_lost", "", "the stop written as deleted with no extension count"),
    ];
    for (reference, candidate, terms, flags, shape) in rows {
        assert_eq!(
            stop_codon_verdicts(reference, candidate, terms, flags),
            only("red"),
            "{shape}"
        );
    }
    assert_eq!(
        verdicts_under(
            STOP_CODON_CLASS,
            "Ter284=",
            "Ter284LysextTer47",
            "stop_lost",
            "",
            "ENSP00000000002.1"
        ),
        only("red")
    );
}

/// A reference description that names the stop alone, as unchanged, deleted, replaced by one
/// residue or extended, where the candidate writes no HGVSp, is set aside once per row in the
/// default, tab, VCF and JSON formats on a row whose FLAGS name `cds_end_NF`, whose terms
/// hold `incomplete_terminal_codon_variant`, or whose description is `Ter1=` under FLAGS
/// naming `cds_start_NF`; it is a difference without the class.
#[test]
fn a_stop_named_on_an_incomplete_coding_sequence_is_set_aside() {
    #[rustfmt::skip]
    let rows: &[(&str, &str, &str)] = &[
        // The stop unchanged, the flag and the term.
        ("Ter305=", "incomplete_terminal_codon_variant,coding_sequence_variant", "cds_end_NF"),
        // The stop unchanged, the term alone.
        ("Ter305=", "incomplete_terminal_codon_variant,coding_sequence_variant", ""),
        // The stop deleted.
        ("Ter185del", "protein_altering_variant,incomplete_terminal_codon_variant", "cds_end_NF"),
        // The stop deleted, the flag alone.
        ("Ter185del", "frameshift_variant", "cds_start_NF,cds_end_NF"),
        // The stop replaced by one residue.
        ("Ter224Val", "incomplete_terminal_codon_variant,coding_sequence_variant", ""),
        // An extension of unknown count.
        ("Ter176HisextTer?", "frameshift_variant", "cds_end_NF"),
        // An extension with a count.
        ("Ter176HisextTer12", "frameshift_variant", "cds_end_NF"),
        // The stop unchanged at residue 1.
        ("Ter1=", "coding_sequence_variant", "cds_start_NF"),
        // The stop unchanged at residue 1, both flags.
        ("Ter1=", "coding_sequence_variant", "cds_start_NF,cds_end_NF"),
    ];
    for (reference, terms, flags) in rows {
        assert_eq!(
            verdicts_under(INCOMPLETE_CLASS, reference, "", terms, flags, PROTEIN),
            only("green"),
            "{reference} under {terms} and {flags}"
        );
    }
}

/// The nearest shapes stay differences in every format: a row that shows no incomplete coding
/// sequence, the start flag at another residue or on another description, a range deleted up to
/// the stop, a deletion-insertion of the stop alone, an insertion that names a stop at residue
/// 1, a notation that names no stop alone, a candidate that writes a notation, and the sides
/// exchanged.
#[test]
fn the_nearest_shapes_of_the_incomplete_coding_sequence_class_stay_differences() {
    #[rustfmt::skip]
    let rows: &[(&str, &str, &str, &str, &str)] = &[
        ("Ter305=", "", "coding_sequence_variant", "", "no flag and no term"),
        ("Ter1=", "", "coding_sequence_variant", "", "the stop unchanged at residue 1 on an unflagged transcript"),
        ("Ter5=", "", "coding_sequence_variant", "cds_start_NF", "the start flag at another residue than 1"),
        ("Ter1Val", "", "coding_sequence_variant", "cds_start_NF", "the start flag on another description at residue 1"),
        ("Asp223_Ter224del", "", "frameshift_variant", "cds_end_NF", "a range deleted up to the stop"),
        ("Ter224delinsValHis", "", "protein_altering_variant,incomplete_terminal_codon_variant", "cds_end_NF", "a deletion-insertion of the stop alone"),
        ("Ter1_Val2insHis", "", "inframe_insertion", "cds_start_NF", "an insertion that names a stop at residue 1"),
        ("Asp223ValfsTer4", "", "frameshift_variant", "cds_end_NF", "a notation that names no stop alone"),
        ("Ter305=", "Ter305HisextTer?", "stop_lost", "cds_end_NF", "a candidate that writes a notation"),
        ("", "Ter305=", "incomplete_terminal_codon_variant,coding_sequence_variant", "cds_end_NF", "the sides exchanged"),
    ];
    for (reference, candidate, terms, flags, shape) in rows {
        assert_eq!(
            verdicts_under(
                INCOMPLETE_CLASS,
                reference,
                candidate,
                terms,
                flags,
                PROTEIN
            ),
            only("red"),
            "{shape}"
        );
    }
}

/// A row of the incomplete coding sequence class under `--hgvsp_use_prediction`: the
/// reference's description is the stop written as unchanged inside one pair of prediction
/// parentheses, and the row's terms and FLAGS show the incomplete coding sequence.
const PREDICTED_STOP: &str = "(Ter305=)";
const INCOMPLETE_TERMS: &str = "incomplete_terminal_codon_variant,coding_sequence_variant";
const INCOMPLETE_FLAGS: &str = "cds_end_NF";

/// `stop_named_on_incomplete_coding_sequence` on one protein's two descriptions after `p.`,
/// the `=` escaped; an empty description is no notation.
fn incomplete_sequence_reads(reference: &str, candidate: &str, terms: &str, flags: &str) -> bool {
    let value = |description: &str| {
        if description.is_empty() {
            "-".to_string()
        } else {
            hgvsp(description, PROTEIN, true)
        }
    };
    stop_named_on_incomplete_coding_sequence(terms, flags, &value(reference), &value(candidate))
}

/// `incomplete_sequence_reads` on a row whose terms and FLAGS show the incomplete coding
/// sequence.
fn marked_row_reads(reference: &str, candidate: &str) -> bool {
    incomplete_sequence_reads(reference, candidate, INCOMPLETE_TERMS, INCOMPLETE_FLAGS)
}

/// The verdicts of a row under the incomplete coding sequence class, its candidate naming
/// the reference's protein.
fn incomplete_verdicts(
    reference: &str,
    candidate: &str,
    consequence: &str,
    flags: &str,
) -> BTreeSet<String> {
    verdicts_under(
        INCOMPLETE_CLASS,
        reference,
        candidate,
        consequence,
        flags,
        PROTEIN,
    )
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

/// The function reads each description of the class inside one pair of prediction
/// parentheses as it reads it without them, the `=` escaped or not, where the candidate
/// writes no notation.
#[test]
fn the_incomplete_coding_sequence_rule_reads_a_description_inside_prediction_parentheses() {
    assert!(marked_row_reads(PREDICTED_STOP, ""));
    let rows: &[(&str, &str, &str)] = &[
        ("Ter305=", INCOMPLETE_TERMS, "cds_end_NF"),
        ("Ter305=", INCOMPLETE_TERMS, ""),
        ("Ter185del", "frameshift_variant", "cds_start_NF,cds_end_NF"),
        ("Ter224Val", INCOMPLETE_TERMS, ""),
        ("Ter176HisextTer?", "frameshift_variant", "cds_end_NF"),
        ("Ter176HisextTer12", "frameshift_variant", "cds_end_NF"),
        ("Ter1=", "coding_sequence_variant", "cds_start_NF"),
    ];
    for (reference, terms, flags) in rows {
        assert!(
            incomplete_sequence_reads(reference, "", terms, flags),
            "{reference}"
        );
        let predicted = format!("({reference})");
        assert!(
            incomplete_sequence_reads(&predicted, "", terms, flags),
            "{reference}"
        );
        assert!(
            stop_named_on_incomplete_coding_sequence(
                terms,
                flags,
                &format!("{PROTEIN}:p.{predicted}"),
                ABSENT
            ),
            "{reference}"
        );
    }
}

/// Through the comparison, a reference description that names the stop alone inside one
/// pair of prediction parentheses, where the candidate writes no HGVSp, is set aside once per
/// row in the default, tab, VCF and JSON formats; in the default format the HGVSp and the
/// Extra text each count once under the class. It is a difference without the class.
#[test]
fn a_stop_named_inside_prediction_parentheses_is_set_aside_on_an_incomplete_coding_sequence() {
    assert_eq!(
        incomplete_verdicts(PREDICTED_STOP, "", INCOMPLETE_TERMS, INCOMPLETE_FLAGS),
        only("green")
    );
    let text = |hgvsp: &str| {
        let mut extra = format!("IMPACT=LOW;STRAND=1;FLAGS={INCOMPLETE_FLAGS}");
        if !hgvsp.is_empty() {
            extra += &format!(";HGVSp={hgvsp}");
        }
        default_text(&[default_row(
            "stop_row",
            "21:25004297",
            "-",
            "ENST00000000001",
            INCOMPLETE_TERMS,
            &extra,
        )])
    };
    let (reference, candidate) = (text("ENSP00000000001.1:p.(Ter305%3D)"), text(""));
    let (report, documented) = compare_default(&reference, &candidate, &[INCOMPLETE_CLASS]);
    assert_eq!(verdict(&report), "green");
    let once = count_map(&[(INCOMPLETE_CLASS, 1)]);
    assert_eq!(counts(&documented.applied), once);
    for column in ["HGVSp", "Extra"] {
        assert_eq!(masked(&report, column), once, "{column}");
        assert_eq!(
            number(&report, &["columns", column, "differing"]),
            0,
            "{column}"
        );
    }
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
}

/// A candidate that writes a notation, inside prediction parentheses or not, stays a
/// difference against the stop named alone inside them.
#[test]
fn a_candidate_notation_stays_a_difference_beside_prediction_parentheses() {
    assert!(marked_row_reads(PREDICTED_STOP, ""));
    for candidate in ["(Ter305HisextTer?)", "Ter305HisextTer?", "Ter305="] {
        assert!(!marked_row_reads(PREDICTED_STOP, candidate), "{candidate}");
    }
    assert_eq!(
        incomplete_verdicts(
            PREDICTED_STOP,
            "(Ter305HisextTer?)",
            INCOMPLETE_TERMS,
            INCOMPLETE_FLAGS
        ),
        only("red")
    );
}

/// A description that names more than the stop alone inside the pair stays a difference: a
/// range deleted up to the stop, a deletion-insertion of the stop, an insertion that names a
/// stop at residue 1 and a notation that names no stop alone, on the row's own terms and
/// FLAGS.
#[test]
fn more_than_the_stop_alone_stays_a_difference_inside_prediction_parentheses() {
    assert!(marked_row_reads(PREDICTED_STOP, ""));
    for reference in [
        "(Asp223_Ter224del)",
        "(Ter224delinsValHis)",
        "(Ter1_Val2insHis)",
        "(Asp223ValfsTer4)",
    ] {
        assert!(!marked_row_reads(reference, ""), "{reference}");
    }
    assert_eq!(
        incomplete_verdicts("(Asp223_Ter224del)", "", INCOMPLETE_TERMS, INCOMPLETE_FLAGS),
        only("red")
    );
}

/// A row that carries neither the flag nor the term stays a difference with its description
/// inside prediction parentheses, and so does the start flag at another residue than 1.
#[test]
fn an_unmarked_row_stays_a_difference_inside_prediction_parentheses() {
    assert!(marked_row_reads(PREDICTED_STOP, ""));
    let unmarked = "coding_sequence_variant";
    assert!(!incomplete_sequence_reads(PREDICTED_STOP, "", unmarked, ""));
    assert!(!incomplete_sequence_reads("(Ter1=)", "", unmarked, ""));
    assert!(!incomplete_sequence_reads(
        "(Ter5=)",
        "",
        unmarked,
        "cds_start_NF"
    ));
    assert_eq!(
        incomplete_verdicts(PREDICTED_STOP, "", unmarked, ""),
        only("red")
    );
}

/// One pair of parentheses is set aside and nothing else: an opening parenthesis alone, a
/// closing one alone and a second pair around the first stay differences.
#[test]
fn one_parenthesis_without_its_pair_stays_a_difference() {
    assert!(marked_row_reads(PREDICTED_STOP, ""));
    for reference in ["(Ter305=", "Ter305=)", "((Ter305=))"] {
        assert!(!marked_row_reads(reference, ""), "{reference}");
        assert_eq!(
            incomplete_verdicts(reference, "", INCOMPLETE_TERMS, INCOMPLETE_FLAGS),
            only("red"),
            "{reference}"
        );
    }
}
