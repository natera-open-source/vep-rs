// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvsp_frameshift_with_no_changed_residue` on the HGVSp of a frameshift's row: the
//! reference writes a frameshift with the unknown count that names one residue on both sides
//! of its number, and the candidate writes no HGVSp.
//! The class takes the field from the row alone, once a row, in the default, tab, VCF and
//! JSON formats, and leaves every other pair a difference.

mod compare_common;
mod support;

use compare_common::*;
use indexmap::IndexMap;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::path::Path;
use support::{args_of, compact, member, report, run};
use vep_formats_compare::cache_classes;
use vep_formats_compare::classes;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::frameshift_no_changed_residue::frameshift_with_no_changed_residue;
use vep_formats_compare::shared::ABSENT;

const CLASS: &str = "hgvsp_frameshift_with_no_changed_residue";
const PROTEIN: &str = "ENSP00000000001.1";
const HGVSC: &str = "ENST00000000001.1:c.32del";
const INPUT: &str = "21\t25004296\tfs_row\tAC\tA\t.\t.\t.";
const FRAMESHIFT: &str = "frameshift_variant";

/// The made-up description of the class: a frameshift with the unknown count that names a
/// leucine on both sides of residue 11.
const UNCHANGED: &str = "Leu11LeufsTer?";

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

/// The verdict of one comparison under the class: a green one has set the one row aside
/// under the class, and a red one has set nothing aside.
fn checked(report: &Value, documented: &Documented) -> String {
    let verdict = verdict(report).to_string();
    assert!(verdict == "green" || verdict == "red");
    let applied = if verdict == "green" {
        count_map(&[(CLASS, 1)])
    } else {
        count_map(&[])
    };
    assert_eq!(counts(&documented.applied), applied);
    verdict
}

/// A default-format output of one Transcript row whose terms are `consequence` and whose
/// HGVSp is `hgvsp`, empty for none.
fn default_output(hgvsp: &str, consequence: &str) -> String {
    let mut extra = format!("IMPACT=HIGH;STRAND=1;HGVSc={HGVSC}");
    if !hgvsp.is_empty() {
        extra += &format!(";HGVSp={hgvsp}");
    }
    default_text(&[default_row(
        "fs_row",
        "21:25004297",
        "-",
        "ENST00000000001",
        consequence,
        &extra,
    )])
}

/// The verdicts of one synthetic row under the class in the default, tab, VCF and JSON
/// formats: the row's HGVSp is the description `reference` after `p.` on the one side and
/// `candidate` on the other (empty for no notation), and its terms are `consequence` on
/// both sides. Without the class the default-format row is a difference.
fn verdicts(reference: &str, candidate: &str, consequence: &str) -> BTreeSet<String> {
    let mut verdicts = BTreeSet::new();
    let sides = [hgvsp(reference, true), hgvsp(candidate, true)];

    let default = |hgvsp: &str| default_output(hgvsp, consequence);
    let (report, documented) = compare_default(&default(&sides[0]), &default(&sides[1]), &[CLASS]);
    verdicts.insert(checked(&report, &documented));
    assert_eq!(
        verdict_of(&default(&sides[0]), &default(&sides[1]), &[]),
        "red"
    );

    let tab_header = DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tHGVSc\tHGVSp");
    let tab = |hgvsp: &str| {
        let cells = [
            "fs_row",
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
            HGVSC,
            if hgvsp.is_empty() { "-" } else { hgvsp },
        ];
        format!("{tab_header}\n{}\n", cells.join("\t"))
    };
    let mut documented = documented_with(&[CLASS]);
    let report = compare(
        &parse_tab(&tab(&sides[0])).expect("the reference parses"),
        &parse_tab(&tab(&sides[1])).expect("the candidate parses"),
        &mut documented,
    );
    verdicts.insert(checked(&report, &documented));

    let inputs = vec![INPUT.to_string()];
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|HGVSc|HGVSp\">\n";
    let vcf = |hgvsp: &str| {
        format!(
            "{csq_header}{}CSQ=-|{}|ENST00000000001|Transcript|{HGVSC}|{hgvsp}\n",
            INPUT.strip_suffix('.').expect("an empty INFO column"),
            consequence.replace(',', "&"),
        )
    };
    let mut documented = documented_with(&[CLASS]);
    let report = compare(
        &parse_vcf(&vcf(&sides[0]), &inputs).expect("the reference parses"),
        &parse_vcf(&vcf(&sides[1]), &inputs).expect("the candidate parses"),
        &mut documented,
    );
    verdicts.insert(checked(&report, &documented));

    let terms: Vec<&str> = consequence.split(',').collect();
    let json = |description: &str| {
        let mut entry = vec![
            ("variant_allele", Value::Str("-".to_string())),
            ("transcript_id", Value::Str("ENST00000000001".to_string())),
            ("consequence_terms", texts(&terms)),
            ("hgvsc", Value::Str(HGVSC.to_string())),
        ];
        if !description.is_empty() {
            entry.push(("hgvsp", Value::Str(hgvsp(description, false))));
        }
        let record = object(vec![
            ("input", Value::Str(INPUT.to_string())),
            ("id", Value::Str("fs_row".to_string())),
            ("allele_string", Value::Str("C/-".to_string())),
            ("most_severe_consequence", Value::Str(terms[0].to_string())),
            ("transcript_consequences", Value::List(vec![object(entry)])),
        ]);
        pyjson::dumps(&record, &pyjson::DumpOptions::default()) + "\n"
    };
    let mut documented = documented_with(&[CLASS]);
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

/// The walk of `_get_fs_peptides` leaves by its bound where the alternate translation ends
/// before a residue differs, and returns the pair of residues read last, two equal ones
/// (TranscriptVariationAllele.pm:2348, :2351-2352, :2361 and :2364). Such a description,
/// where the candidate writes no HGVSp, is set aside once per row in the default, tab, VCF
/// and JSON formats on a row whose terms hold `frameshift_variant`, alone or beside other
/// terms, and inside a pair of prediction parentheses as well.
#[test]
fn a_frameshift_that_names_one_residue_on_both_sides_is_set_aside_against_no_value() {
    #[rustfmt::skip]
    let rows: &[(&str, &str)] = &[
        (UNCHANGED, FRAMESHIFT),
        ("(Leu11LeufsTer?)", FRAMESHIFT),
        ("Ser204SerfsTer?", "frameshift_variant,splice_region_variant"),
        ("Sec5SecfsTer?", "splice_region_variant,frameshift_variant"),
    ];
    for (reference, terms) in rows {
        assert_eq!(
            verdicts(reference, "", terms),
            only("green"),
            "{reference} under {terms}"
        );
    }
}

/// In the default format the HGVSp and the Extra text of the row each count once under the
/// class and no column is left differing; without the class the row is a difference.
#[test]
fn the_default_format_counts_the_field_and_the_extra_text_once_each() {
    let reference = default_output(&hgvsp(UNCHANGED, true), FRAMESHIFT);
    let candidate = default_output("", FRAMESHIFT);
    let (report, documented) = compare_default(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    let once = count_map(&[(CLASS, 1)]);
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

/// A candidate that writes any notation stays a difference in every format, and so do the
/// two sides exchanged: the reference's description against a frameshift with a changed
/// residue, the same residues with a count, the residue unchanged and a notation inside
/// prediction parentheses, and no reference value against the candidate's description.
#[test]
fn a_candidate_notation_and_the_sides_exchanged_stay_differences() {
    #[rustfmt::skip]
    let rows: &[(&str, &str, &str)] = &[
        ("", UNCHANGED, "the sides exchanged"),
        (UNCHANGED, "Leu11ProfsTer5", "a candidate frameshift with a changed residue"),
        (UNCHANGED, "Leu11LeufsTer4", "a candidate that writes the same residues with a count"),
        (UNCHANGED, "Leu11=", "a candidate that writes the residue unchanged"),
        ("(Leu11LeufsTer?)", "(Leu11ProfsTer5)", "a candidate notation inside prediction parentheses"),
    ];
    for (reference, candidate, shape) in rows {
        assert_eq!(
            verdicts(reference, candidate, FRAMESHIFT),
            only("red"),
            "{shape}"
        );
    }
}

/// Any other reference description stays a difference in every format where the candidate
/// writes no HGVSp: two different residue names, a count that is a number, the stop named on
/// both sides, a residue number that is no residue's, a text after the unknown count, and a
/// description that is no frameshift.
#[test]
fn another_reference_description_stays_a_difference() {
    #[rustfmt::skip]
    let rows: &[(&str, &str)] = &[
        ("Leu11ProfsTer?", "two different residue names"),
        ("Leu11LeufsTer4", "a count that is a number"),
        ("Ter11TerfsTer?", "the stop named on both sides"),
        ("Leu0LeufsTer?", "the residue number zero"),
        ("Leu011LeufsTer?", "a residue number that opens with a zero"),
        ("Leu11LeufsTer?Val", "a text after the unknown count"),
        ("Leu11Leufs", "a frameshift with no count"),
        ("Leu11=", "the residue unchanged"),
        ("Leu11del", "a deletion"),
        ("Leu11_Val12insLeu", "an insertion"),
        ("Ter11LeuextTer?", "an extension of unknown count"),
    ];
    for (reference, shape) in rows {
        assert_eq!(verdicts(reference, "", FRAMESHIFT), only("red"), "{shape}");
    }
}

/// A row whose terms hold no `frameshift_variant` stays a difference in every format.
#[test]
fn a_row_whose_terms_hold_no_frameshift_stays_a_difference() {
    for terms in [
        "inframe_deletion",
        "coding_sequence_variant",
        "stop_retained_variant,splice_region_variant",
    ] {
        assert_eq!(verdicts(UNCHANGED, "", terms), only("red"), "{terms}");
    }
}

/// One pair of parentheses is set aside and nothing else: an opening parenthesis alone, a
/// closing one alone and a second pair around the first stay differences.
#[test]
fn one_parenthesis_without_its_pair_stays_a_difference() {
    for reference in ["(Leu11LeufsTer?", "Leu11LeufsTer?)", "((Leu11LeufsTer?))"] {
        assert_eq!(
            verdicts(reference, "", FRAMESHIFT),
            only("red"),
            "{reference}"
        );
    }
}

/// The rule reads the reference's value, the candidate's and the terms of both rows: a
/// candidate that writes no HGVSp in each spelling of no value, the terms joined by `,` or
/// by `&`, and no candidate row. A row of either side whose terms hold no
/// `frameshift_variant`, a reference value with no protein before `p.`, a candidate value
/// and the sides exchanged are no row of the class.
#[test]
fn the_rule_reads_the_two_values_and_the_terms_of_both_rows() {
    let value = hgvsp(UNCHANGED, true);
    let both = Some(FRAMESHIFT);
    for none in ["", "-", ABSENT] {
        assert!(
            frameshift_with_no_changed_residue(FRAMESHIFT, both, &value, none),
            "{none}"
        );
    }
    assert!(frameshift_with_no_changed_residue(
        FRAMESHIFT, None, &value, "-"
    ));
    for terms in [
        "splice_region_variant,frameshift_variant",
        "frameshift_variant&splice_region_variant",
    ] {
        assert!(
            frameshift_with_no_changed_residue(terms, Some(terms), &value, "-"),
            "{terms}"
        );
    }
    let in_frame = "inframe_deletion";
    assert!(!frameshift_with_no_changed_residue(
        FRAMESHIFT,
        Some(in_frame),
        &value,
        "-"
    ));
    assert!(!frameshift_with_no_changed_residue(
        in_frame, both, &value, "-"
    ));
    for nameless in [UNCHANGED, "p.Leu11LeufsTer?"] {
        assert!(
            !frameshift_with_no_changed_residue(FRAMESHIFT, both, nameless, "-"),
            "{nameless}"
        );
    }
    assert!(!frameshift_with_no_changed_residue(
        FRAMESHIFT, both, &value, &value
    ));
    assert!(!frameshift_with_no_changed_residue(
        FRAMESHIFT, both, "-", &value
    ));
}

/// The class reads the HGVSp under either letter case and the row alone: it stands in no
/// table of the classes that read a cache, a record's rows or its alleles, no release
/// corrects it, and the --masks help names it.
#[test]
fn the_class_tables_and_the_help_carry_the_class() {
    let tables = classes::tables();
    assert_eq!(classes::HGVSP_FRAMESHIFT_WITH_NO_CHANGED_RESIDUE, CLASS);
    assert_eq!(
        tables.field_classes_of(CLASS),
        Some(&["HGVSp", "hgvsp"][..])
    );
    assert!(!tables.in_cache_reading_classes(CLASS));
    assert!(!tables.in_transcript_reading_classes(CLASS));
    assert!(!tables.in_record_reading_classes(CLASS));
    assert!(!tables.in_allele_reading_classes(CLASS));
    assert!(!tables.in_record_reference_classes(CLASS));
    assert!(cache_classes::classes()
        .iter()
        .all(|(class, _)| *class != CLASS));
    assert_eq!(tables.corrected_in(CLASS), None);
    let help: Vec<&str> = include_str!("../data/help.txt")
        .split_whitespace()
        .collect();
    assert!(help.join(" ").contains(&format!(
        "{CLASS} (under --hgvs, a reference HGVSp that is a frameshift with the unknown count and one residue name, no stop, on both sides of its residue number, where the candidate writes none, on a row whose terms hold frameshift_variant on both sides; read from the row alone),"
    )));
}

/// The command line names the class like any other and reads no file beside the two
/// outputs for it: the class is applied once, and without it the pair is a difference.
#[test]
fn the_command_line_names_the_class_and_reads_the_two_outputs_alone() {
    let directory = scratch("frameshift-with-no-changed-residue-command-line");
    let text = |path: &Path| path.to_str().expect("a path of text").to_string();
    let written = |name: &str, hgvsp: &str| {
        text(&planted(
            &directory,
            name,
            &default_output(hgvsp, FRAMESHIFT),
        ))
    };
    let reference = written("reference.txt", &hgvsp(UNCHANGED, true));
    let candidate = written("candidate.txt", "");
    let compared = |more: &[&str], name: &str| {
        let path = directory.join(name);
        let mut argv = args_of(&[
            "--format",
            "default",
            "--expected",
            &reference,
            "--actual",
            &candidate,
            "--report",
            &text(&path),
        ]);
        argv.extend(args_of(more));
        (run(&directory, &argv), path)
    };
    let (ran, path) = compared(&["--masks", CLASS], "named.json");
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(
        compact(member(&report(&path), &["masks", "applied"])),
        format!(r#"{{"{CLASS}": 1}}"#)
    );
    let (ran, path) = compared(&[], "unnamed.json");
    assert_eq!(ran.code, 1, "{}", ran.err());
    assert_eq!(compact(member(&report(&path), &["masks", "applied"])), "{}");
}
