// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The UTRAnnotator class of a run without a reference sequence: the reference loads the
//! plugin and writes its five fields from no sequence, the candidate refuses the plugin and
//! names none of them, read through the comparison of planted rows in the four formats.

mod compare_common;
mod support;
mod utrannotator_common;

use compare_common::*;
use utrannotator_common::*;

const CLASS: &str = UTR_WITHOUT_REFERENCE;

/// The plugin's five fields as the reference writes them on a 5' UTR row in a run without
/// a reference sequence, in the order its Extra carries them (UTRAnnotator.pm:291-298).
const PLACEHOLDERS: &str = "5UTR_annotation=-;5UTR_consequence=-;Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=0";

/// A made-up substitution's row on a 5' UTR, with `fields` after the row's own Extra pairs.
fn utr(fields: &str) -> String {
    let extra = if fields.is_empty() {
        "IMPACT=MODIFIER;STRAND=1".to_string()
    } else {
        format!("IMPACT=MODIFIER;STRAND=1;{fields}")
    };
    default_row(
        "utr_snv",
        "21:1000",
        "A",
        "ENST1",
        "5_prime_UTR_variant",
        &extra,
    )
}

/// The same record's row on an intron of another transcript, which the plugin leaves alone.
fn intron() -> String {
    default_row(
        "utr_snv",
        "21:1000",
        "A",
        "ENST2",
        "intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    )
}

/// The reference's and the candidate's default-format outputs of the two rows.
fn default_pair() -> (String, String) {
    (
        default_text(&[utr(PLACEHOLDERS), intron()]),
        default_text(&[utr(""), intron()]),
    )
}

/// In the default format the reference's Extra carries the five fields on the 5' UTR row
/// alone, so that row is the class's and the intron's row is equal.
#[test]
fn a_run_without_a_reference_is_set_aside_in_the_default_format() {
    let (reference, candidate) = default_pair();
    assert!(compared(&reference, &candidate, &[CLASS]).is("green", &[(CLASS, 1)]));
    assert!(compared(&reference, &candidate, &[]).is("red", &[]));
}

/// The tab format gives every row the five columns the reference's header declares and the
/// candidate's lacks, `-` on a row the plugin left alone, and the VCF format gives every
/// CSQ entry the five fields its `Format` names, empty there: each entry is the class's.
#[test]
fn a_run_without_a_reference_is_set_aside_on_every_tab_and_vcf_entry() {
    let columns =
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tSTRAND";
    let rows = [
        "utr_snv\t21:1000\tA\tENSG1\tENST1\tTranscript\t5_prime_UTR_variant\t1",
        "utr_snv\t21:1000\tA\tENSG1\tENST2\tTranscript\tintron_variant\t1",
    ];
    let reference = format!(
        "{columns}\t5UTR_annotation\t5UTR_consequence\tExisting_InFrame_oORFs\tExisting_OutOfFrame_oORFs\tExisting_uORFs\n{}\t-\t-\t0\t0\t0\n{}\t-\t-\t-\t-\t-\n",
        rows[0], rows[1]
    );
    let candidate = format!("{columns}\n{}\n{}\n", rows[0], rows[1]);
    let (reference, candidate) = (tab_entries(&reference), tab_entries(&candidate));
    let outcome = compared_entries(&reference, &candidate, &[CLASS], None, None);
    assert!(outcome.is("green", &[(CLASS, 2)]));
    assert!(compared_entries(&reference, &candidate, &[], None, None).is("red", &[]));

    let inputs = lines(&["21\t1000\tutr_snv\tG\tA\t.\t.\t."]);
    let entries = "Allele|Consequence|Feature|Feature_type|STRAND";
    let reference = format!(
        "{}{}CSQ=A|5_prime_UTR_variant|ENST1|Transcript|1|||0|0|0,A|intron_variant|ENST2|Transcript|1|||||\n",
        csq_header(&format!(
            "{entries}|5UTR_annotation|5UTR_consequence|Existing_InFrame_oORFs|Existing_OutOfFrame_oORFs|Existing_uORFs"
        )),
        record_columns(&inputs[0])
    );
    let candidate = format!(
        "{}{}CSQ=A|5_prime_UTR_variant|ENST1|Transcript|1,A|intron_variant|ENST2|Transcript|1\n",
        csq_header(entries),
        record_columns(&inputs[0])
    );
    let reference = vcf_entries(&reference, &inputs);
    let candidate = vcf_entries(&candidate, &inputs);
    let outcome = compared_entries(&reference, &candidate, &[CLASS], Some(&inputs), None);
    assert!(outcome.is("green", &[(CLASS, 2)]));
    let unnamed = compared_entries(&reference, &candidate, &[], Some(&inputs), None);
    assert!(unnamed.is("red", &[]));
}

/// The JSON format writes the three counts as numbers on the 5' UTR's consequence object
/// and leaves out a text field that holds no value.
#[test]
fn a_run_without_a_reference_is_set_aside_in_the_json_format() {
    let inputs = lines(&["21\t1000\tutr_snv\tG\tA\t.\t.\t."]);
    let members = r#""variant_allele": "A", "transcript_id": "ENST1", "consequence_terms": ["5_prime_UTR_variant"], "strand": 1"#;
    let counts =
        r#""existing_inframe_oorfs": 0, "existing_outofframe_oorfs": 0, "existing_uorfs": 0"#;
    let reference = json_entries(
        &json_line(
            &inputs[0],
            "utr_snv",
            "G/A",
            &format!("{members}, {counts}"),
        ),
        &inputs,
    );
    let candidate = json_entries(&json_line(&inputs[0], "utr_snv", "G/A", members), &inputs);
    let outcome = compared_entries(&reference, &candidate, &[CLASS], Some(&inputs), None);
    assert!(outcome.is("green", &[(CLASS, 1)]));
    let unnamed = compared_entries(&reference, &candidate, &[], Some(&inputs), None);
    assert!(unnamed.is("red", &[]));
}

/// The class takes a row only where the candidate names none of the five fields and the
/// reference wrote them from no sequence: a candidate that ran the plugin, a candidate that
/// carries one of the fields, a count the reference read as other than `0` and a term other
/// than the gained start each stay a difference. So do one empty count beside two of `0`, and a
/// count or a term the reference's row carries under a second letter case.
#[test]
fn a_row_the_candidate_fills_or_the_reference_annotates_stays_a_difference() {
    let counted = PLACEHOLDERS.replace("Existing_uORFs=0", "Existing_uORFs=1");
    let uncounted = PLACEHOLDERS.replace("Existing_uORFs=0", "Existing_uORFs=-");
    let recounted = format!("{PLACEHOLDERS};existing_uorfs=1");
    let renamed = format!("{PLACEHOLDERS};5utr_consequence={UTR_LOSS}");
    let reported = PLACEHOLDERS
        .replace(
            "5UTR_annotation=-",
            "5UTR_annotation=CapDistanceToStart=63:type=uORF",
        )
        .replace(
            "5UTR_consequence=-",
            &format!("5UTR_consequence={UTR_LOSS}"),
        );
    let cases: [(&str, &str); 7] = [
        (PLACEHOLDERS, &counted),
        (PLACEHOLDERS, "Existing_uORFs=0"),
        (&counted, ""),
        (&reported, ""),
        (&uncounted, ""),
        (&recounted, ""),
        (&renamed, ""),
    ];
    for (reference, candidate) in cases {
        let outcome = compared(
            &default_text(&[utr(reference), intron()]),
            &default_text(&[utr(candidate), intron()]),
            &[CLASS],
        );
        assert!(outcome.is("red", &[]), "{reference} against {candidate}");
    }
}

/// A substitution whose alternate allele spells a start codon gets the gained start from
/// the reference, written over no sequence beside three counts of `0`
/// (UTRAnnotator.pm:345-354): that row is the class's too, while the same term beside a
/// count the reference read, or beside no count, stays a difference.
#[test]
fn a_gained_start_from_the_alleles_own_bases_is_set_aside() {
    let annotation = "5UTR_annotation=CapDistanceToStart=23:type=inFrame_oORF";
    let gained = PLACEHOLDERS
        .replace("5UTR_annotation=-", annotation)
        .replace(
            "5UTR_consequence=-",
            &format!("5UTR_consequence={UTR_GAIN}"),
        );
    let candidate = default_text(&[utr(""), intron()]);
    let outcome = compared(
        &default_text(&[utr(&gained), intron()]),
        &candidate,
        &[CLASS],
    );
    assert!(outcome.is("green", &[(CLASS, 1)]));
    let counted = gained.replace("Existing_uORFs=0", "Existing_uORFs=1");
    let uncounted = format!("{annotation};5UTR_consequence={UTR_GAIN}");
    for fields in [&counted, &uncounted] {
        let outcome = compared(
            &default_text(&[utr(fields), intron()]),
            &candidate,
            &[CLASS],
        );
        assert!(outcome.is("red", &[]), "{fields}");
    }
}

/// The command line names the class like any other.
#[test]
fn the_command_line_names_the_class_of_a_run_without_a_reference() {
    let (reference, candidate) = default_pair();
    let directory = planted_files(
        "utrannotator-without-reference",
        &[("reference", &reference), ("candidate", &candidate)],
    );
    assert_eq!(exit_code(&directory, &[]), 1);
    assert_eq!(exit_code(&directory, &["--masks", CLASS]), 0);
}
