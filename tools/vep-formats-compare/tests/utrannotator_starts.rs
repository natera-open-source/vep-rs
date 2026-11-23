// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The UTRAnnotator classes of a start codon: the gained-start window, a lost start an
//! insertion only moves on and a gained start an insertion leaves in place, read through the
//! comparison of planted rows.

mod compare_common;
mod support;
mod utrannotator_common;

use compare_common::*;
use utrannotator_common::*;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::record_rows::RecordRows;
use vep_formats_compare::self_parity::reverse_complement;
use vep_formats_compare::shared::allele_is_indel;
use vep_formats_compare::shared::ABSENT;

/// The comparison of two one-row outputs under the gained-start window with `inputs` read.
fn window_with(reference: &str, candidate: &str, inputs: &[&str]) -> Outcome {
    compared_with(reference, candidate, &[UTR_WINDOW], Some(&lines(inputs)))
}

/// The insertion 21 25734954 A ACAT on ENST00000000033, cDNA 64-65: Ensembl VEP reports
/// the start lost at 63, the candidate a start gained at 64, the first inserted base's own
/// index, which the window of an empty reference allele does not hold. Under the record's
/// identifier the row is the class's where the input states the record's alleles, counted
/// once with its 5UTR_annotation. So are a deletion and a shorter replacement whose reference
/// reports a start past the allele's end, a longer replacement whose candidate reports one
/// past the reference allele's end, and a reference start at a negative CapDistanceToStart.
#[test]
fn the_gained_start_window_accepts_a_start_one_window_alone_holds() {
    let rows: [(Row, Values, Values); 6] = [
        (at_cdna(UTR_INSERTION, "64-65"), utr_lost(), gained(64)),
        (at_cdna(UTR_DELETION, "336-340"), gained(337), utr_shifted()),
        (
            ("shorter", "21:2000-2004", "GCA", "ENST3", "402-406"),
            gained(405),
            nothing(),
        ),
        (
            ("longer", "21:3000", "GCA", "ENST4", "112"),
            utr_shifted(),
            gained(113),
        ),
        (
            ("wrapped", "21:4000", "A", "ENST5", "2"),
            gained(-1),
            nothing(),
        ),
        (
            ("wrapped_del", "21:5000", "-", "ENST6", "1"),
            gained(-2),
            nothing(),
        ),
    ];
    let reference: Vec<String> = rows
        .iter()
        .map(|(row, own, _)| utr_row(*row, own))
        .collect();
    let candidate: Vec<String> = rows
        .iter()
        .map(|(row, _, own)| utr_row(*row, own))
        .collect();
    let (reference, candidate) = (default_text(&reference), default_text(&candidate));
    let outcome = window_with(&reference, &candidate, &UTR_INPUTS);
    assert!(outcome.is("green", &[(UTR_WINDOW, 6)]));
    assert!(outcome.masks("5UTR_consequence", &[(UTR_WINDOW, 6)]));
    assert!(outcome.masks("5UTR_annotation", &[(UTR_WINDOW, 6)]));
    let unnamed = compared_with(&reference, &candidate, &[], Some(&lines(&UTR_INPUTS)));
    assert_eq!(unnamed.verdict(), "red");
    // A two-base Location under an identifier is an insertion's or a two-base reference
    // allele's: without the input that row alone stays a difference.
    let outcome = compared(&reference, &candidate, &[UTR_WINDOW]);
    assert!(outcome.is("red", &[(UTR_WINDOW, 5)]));
    assert_eq!(outcome.differing("5UTR_consequence"), 1);
    // An identifier-less Uploaded_variation spells the record's alleles itself.
    let spelt = |text: &str| text.replace("utr_ins", "21_25734954_A/ACAT");
    let outcome = compared(&spelt(&reference), &spelt(&candidate), &[UTR_WINDOW]);
    assert!(outcome.is("green", &[(UTR_WINDOW, 6)]));
}

/// The reference's and the candidate's one-row outputs of the insertion 21 25734954 A ACAT:
/// a start lost at 63 against a start gained at 64.
fn insertion_pair() -> (String, String) {
    let row = at_cdna(UTR_INSERTION, "64-65");
    (
        one_row(utr_row(row, &utr_lost())),
        one_row(utr_row(row, &gained(64))),
    )
}

/// The rows of two records that share the identifier and the Location: the insertion's and
/// one of the allele `G`.
fn two_record_entries() -> Vec<Entry> {
    default_entries(&default_text(&[
        utr_row(at_cdna(UTR_INSERTION, "64-65"), &utr_lost()),
        utr_row(
            ("utr_ins", UTR_INSERTION.1, "G", UTR_INSERTION.3, "64-65"),
            &nothing(),
        ),
    ]))
}

/// The insertion 21 25734954 A ACAT under an identifier: where one input line carries the
/// identifier the row's record is that line; where several carry it the record is the one
/// line that writes the row at its Location with its allele, whether the others stand at
/// the same place with another allele, at another place or on another chromosome. Two lines
/// that both write the row, one of them by two of its alternates in two ways included, lines
/// of which none writes it and an identifier no line carries hand nothing, and the row stays
/// a difference; one line alone that spells another allele is handed as the record and the
/// class refuses it by its own test.
#[test]
fn a_row_under_an_identifier_several_input_lines_carry_finds_its_record_by_place_and_alleles() {
    let (reference, candidate) = insertion_pair();
    let own = UTR_INPUTS[0];
    let other_allele = "21\t25734954\tutr_ins\tA\tAG\t.\t.\t.";
    let elsewhere = "21\t25734990\tutr_ins\tA\tACAT\t.\t.\t.";
    let other_chromosome = "22\t25734954\tutr_ins\tA\tACAT\t.\t.\t.";
    // A two-base reference allele with the row's allele as one alternate and as the middle of another.
    let two_ways = "21\t25734954\tutr_ins\tAC\tCAT,ACATC\t.\t.\t.";
    let prefixed = format!("chr{own}");
    let renamed = own.replace("utr_ins", "other");
    let grid: [(&str, Vec<&str>, bool); 11] = [
        ("one line", vec![own], true),
        (
            "two lines at one place with different alleles",
            vec![other_allele, own],
            true,
        ),
        ("two lines at different places", vec![elsewhere, own], true),
        (
            "two lines on different chromosomes",
            vec![own, other_chromosome],
            true,
        ),
        (
            "a line with the prefixed chromosome among three",
            vec![other_allele, &prefixed, elsewhere],
            true,
        ),
        ("two lines that both write the row", vec![own, own], false),
        (
            "a line whose two alternates write the row in two ways, beside the insertion's line",
            vec![two_ways, own],
            false,
        ),
        (
            "the insertion's line beside a line whose two alternates write the row in two ways",
            vec![own, two_ways],
            false,
        ),
        (
            "two lines of which none writes the row",
            vec![other_allele, elsewhere],
            false,
        ),
        (
            "one line that spells another allele",
            vec![other_allele],
            false,
        ),
        ("an identifier no line carries", vec![&renamed], false),
    ];
    for (what, given, taken) in &grid {
        let outcome = window_with(&reference, &candidate, given);
        assert!(outcome.taken_is(*taken, &[(UTR_WINDOW, 1)]), "{what}");
    }
    // The rows of two records that share the identifier and the Location are each handed their
    // own record alone; one line of two alternates keeps both of its rows.
    let entries = two_record_entries();
    let two_lines = lines(&[own, other_allele]);
    let mut two_records = RecordRows::new(&entries, Some(&two_lines));
    for (entry, alternate) in entries.iter().zip(["ACAT", "AG"]) {
        assert_eq!(
            handed(&mut two_records, entry),
            Some(vec![
                (entry.allele.clone(), None),
                (String::new(), Some(alternate.to_string()))
            ]),
            "{alternate}"
        );
    }
    let one_line = lines(&["21\t25734954\tutr_ins\tA\tACAT,AG\t.\t.\t."]);
    let mut one_record = RecordRows::new(&entries, Some(&one_line));
    assert_eq!(
        handed(&mut one_record, &entries[0]),
        Some(vec![
            ("CAT".to_string(), None),
            ("G".to_string(), None),
            (String::new(), Some("ACAT,AG".to_string()))
        ])
    );
    let twice = lines(&[own, own]);
    let mut both = RecordRows::new(&entries, Some(&twice));
    assert_eq!(handed(&mut both, &entries[0]), None);
    assert!(!both.one_record(&entries[0]).expect("the lines are read"));
}

/// The insertion 21 25734954 A ACAT on a line whose ID column joins several identifiers:
/// both engines print the first of them, so the row's record is that line, whether it joins
/// two identifiers or five, and the entry that stands in for the line carries the first
/// identifier as its name. A line that carries the row's identifier as a later element, or
/// behind a first element that is a dot, empty or `0`, is no record of the row, which stays
/// a difference.
#[test]
fn a_row_under_the_first_identifier_of_a_several_identifier_line_finds_its_record() {
    let (reference, candidate) = insertion_pair();
    let own = UTR_INPUTS[0];
    let grid = [
        ("two identifiers", "utr_ins;rs2", true),
        ("five identifiers", "utr_ins;rs2;rs3;rs4;rs5", true),
        ("the identifier as a later element", "rs2;utr_ins", false),
        ("a first element that is a dot", ".;utr_ins", false),
        ("a first element that is empty", ";utr_ins", false),
        ("a first element that is zero", "0;utr_ins", false),
    ];
    for (what, name, taken) in grid {
        let line = own.replace("utr_ins", name);
        let outcome = window_with(&reference, &candidate, &[&line]);
        assert!(outcome.taken_is(taken, &[(UTR_WINDOW, 1)]), "{what}");
    }
    let entries = default_entries(&reference);
    let joined = vec![own.replace("utr_ins", "utr_ins;rs2;rs3;rs4;rs5")];
    let mut records = RecordRows::new(&entries, Some(&joined));
    assert!(records.one_record(&entries[0]).expect("the lines are read"));
    let line = records
        .line_of(&entries[0])
        .expect("the lines are read")
        .expect("the row has a line");
    assert_eq!(line[..5], ["21", "25734954", "utr_ins", "A", "ACAT"]);
    let found: Vec<(String, Option<String>, Option<String>)> = records
        .of(&entries[0])
        .expect("the record is read")
        .expect("the row has a record")
        .iter()
        .map(|row| {
            (
                row.allele.clone(),
                row.fields.get("Uploaded_variation").cloned(),
                row.fields.get("ALT").cloned(),
            )
        })
        .collect();
    assert_eq!(
        found,
        [
            ("CAT".to_string(), Some("utr_ins".to_string()), None),
            (
                String::new(),
                Some("utr_ins".to_string()),
                Some("ACAT".to_string())
            )
        ]
    );
}

/// A line of one identifier and a line that joins that identifier with another stand under
/// one name, as two joined lines with one first identifier do. A row under the name is handed
/// the one of them that writes it at its Location with its allele, whichever of the two it
/// is and whether the other stands elsewhere or spells another allele; where both write the
/// row nothing is handed. The rows of two such records at one Location are each handed their
/// own record alone.
#[test]
fn lines_that_share_a_first_identifier_are_told_apart_by_place_and_allele() {
    let (reference, candidate) = insertion_pair();
    let own = UTR_INPUTS[0];
    let other_allele = "21\t25734954\tutr_ins\tA\tAG\t.\t.\t.";
    let elsewhere = "21\t25734990\tutr_ins\tAC\tA\t.\t.\t.";
    let joined = |line: &str, second: &str| line.replace("utr_ins", &format!("utr_ins;{second}"));
    let grid: [(&str, Vec<String>, bool); 7] = [
        (
            "the joined line writes the row, the plain line stands elsewhere",
            vec![elsewhere.to_string(), joined(own, "rs2")],
            true,
        ),
        (
            "the joined line writes the row, the plain line spells another allele",
            vec![other_allele.to_string(), joined(own, "rs2")],
            true,
        ),
        (
            "the plain line writes the row, the joined line stands elsewhere",
            vec![own.to_string(), joined(elsewhere, "rs2")],
            true,
        ),
        (
            "the plain line writes the row, the joined line spells another allele",
            vec![own.to_string(), joined(other_allele, "rs2")],
            true,
        ),
        (
            "both lines write the row",
            vec![own.to_string(), joined(own, "rs2")],
            false,
        ),
        (
            "two joined lines of which one writes the row",
            vec![joined(elsewhere, "rs3"), joined(own, "rs2")],
            true,
        ),
        (
            "two joined lines that both write the row",
            vec![joined(own, "rs3"), joined(own, "rs2")],
            false,
        ),
    ];
    for (what, given, taken) in &grid {
        let outcome = compared_with(&reference, &candidate, &[UTR_WINDOW], Some(given));
        assert!(outcome.taken_is(*taken, &[(UTR_WINDOW, 1)]), "{what}");
    }
    let entries = two_record_entries();
    let two_lines = vec![joined(own, "rs2"), other_allele.to_string()];
    let mut two_records = RecordRows::new(&entries, Some(&two_lines));
    for (entry, alternate) in entries.iter().zip(["ACAT", "AG"]) {
        assert_eq!(
            handed(&mut two_records, entry),
            Some(vec![
                (entry.allele.clone(), None),
                (String::new(), Some(alternate.to_string()))
            ]),
            "{alternate}"
        );
    }
}

/// The same value pair stays a difference on a one-allele substitution's row, on a
/// two-base substitution's, on an insertion one base further on, whose reference window
/// holds the start, and at a start past the candidate's own window. So do a start the
/// engine with the shorter allele alone reports, a deletion's start before or past the
/// deleted bases, two other terms, two gained starts, one term with two annotations, a
/// cDNA position that is not the allele's own and a gained start that names no position.
#[test]
fn the_gained_start_window_refuses_a_start_the_two_windows_do_not_explain() {
    let insertion = |cdna: &'static str| at_cdna(UTR_INSERTION, cdna);
    let deletion = |cdna: &'static str| at_cdna(UTR_DELETION, cdna);
    let cases: [(Row, Values, Values); 14] = [
        (
            ("snv", "21:600", "T", "ENST1", "64"),
            utr_lost(),
            gained(64),
        ),
        (
            ("mnv", "21:500-501", "GT", "ENST1", "64-65"),
            utr_lost(),
            gained(64),
        ),
        (insertion("65-66"), utr_lost(), gained(64)),
        (insertion("64-65"), utr_lost(), gained(67)),
        (insertion("64-65"), gained(64), utr_lost()),
        (deletion("336-340"), utr_shifted(), gained(337)),
        (deletion("336-340"), gained(334), utr_shifted()),
        (deletion("336-340"), gained(340), utr_shifted()),
        (insertion("64-65"), utr_lost(), utr_shifted()),
        (insertion("64-65"), gained(63), gained(64)),
        (
            insertion("64-65"),
            utr_lost(),
            values(UTR_LOSS, "CapDistanceToStart=63:type=inframe_oORF"),
        ),
        (insertion("?-65"), utr_lost(), gained(64)),
        (deletion("336-338"), gained(337), nothing()),
        (
            insertion("64-65"),
            utr_lost(),
            values(UTR_GAIN, "type=uORF"),
        ),
    ];
    for (row, reference, candidate) in &cases {
        let outcome = window_with(
            &one_row(utr_row(*row, reference)),
            &one_row(utr_row(*row, candidate)),
            &UTR_INPUTS,
        );
        assert!(
            outcome.is("red", &[]),
            "{row:?} {reference:?} {candidate:?}"
        );
    }
}

/// A start codon the reference reports lost at index 63 or 62 stands whole three bases on
/// when the three inserted bases end with its first base or its first two: a candidate's
/// gained start at that index plus the allele's length names the reference's own codon, not
/// one the inserted bases make, and the row stays a difference; so does the one-base
/// insertion whose gained start stands one base past the lost codon. A gained start at any
/// other index of the allele's span is the class's beside the same lost start, and so is
/// the moved codon's index beside a reference that reports another term or nothing, which
/// states no lost codon to hold it against. The JSON and VCF entries read alike.
#[test]
fn the_gained_start_window_refuses_a_start_codon_the_insertion_only_moved_on() {
    let one_base: Row = ("one", "21:600-601", "T", "ENST00000000033", "64-65");
    let mut inputs = UTR_INPUTS.to_vec();
    inputs.push("21\t600\tone\tA\tAT\t.\t.\t.");
    let lost_two_below = values(
        UTR_LOSS,
        "CapDistanceToStart=62:DistanceToCDS=111:KozakStrength=Moderate:type=uORF",
    );
    let insertion = at_cdna(UTR_INSERTION, "64-65");
    let cases: [(Row, Values, Values, bool); 10] = [
        (insertion, utr_lost(), gained(66), false),
        (insertion, lost_two_below.clone(), gained(65), false),
        (one_base, utr_lost(), gained(64), false),
        (insertion, utr_lost(), gained(64), true),
        (insertion, utr_lost(), gained(65), true),
        (insertion, lost_two_below.clone(), gained(64), true),
        (insertion, lost_two_below, gained(66), true),
        (insertion, utr_shifted(), gained(66), true),
        (insertion, nothing(), gained(66), true),
        (
            insertion,
            values(
                UTR_LOSS,
                "CapDistanceToStart=63:type=uORF,CapDistanceToStart=60:type=uORF",
            ),
            gained(66),
            true,
        ),
    ];
    for (row, reference, candidate, taken) in &cases {
        let outcome = window_with(
            &one_row(utr_row(*row, reference)),
            &one_row(utr_row(*row, candidate)),
            &inputs,
        );
        assert!(
            outcome.taken_is(*taken, &[(UTR_WINDOW, 1)]),
            "{row:?} {reference:?} {candidate:?}"
        );
    }
    let header = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|5UTR_consequence|5UTR_annotation",
    );
    let lost = format!("{UTR_LOSS}|CapDistanceToStart=63:type=uORF");
    let vcf_reference = format!(
        "{header}{}CSQ=CAT|5_prime_UTR_variant|ENST00000000033|Transcript|64-65|{lost}\n",
        record_columns(UTR_INPUTS[0])
    );
    let record = lines(&UTR_INPUTS[..1]);
    let object = |utr: &str, position: i64| {
        json_line(
            UTR_INPUTS[0],
            "r0",
            "-/CAT",
            &format!(
                r#""variant_allele": "CAT", "transcript_id": "ENST00000000033", "consequence_terms": ["5_prime_UTR_variant"], "cdna_start": 65, "cdna_end": 64, "5utr_consequence": "{utr}", "5utr_annotation": {{"1": {{"CapDistanceToStart": {position}, "type": "uORF"}}}}"#
            ),
        )
    };
    for (at, taken) in [(66, false), (65, true)] {
        let vcf_candidate = vcf_reference.replace(
            &lost,
            &format!("{UTR_GAIN}|CapDistanceToStart={at}:type=uORF"),
        );
        let outcome = compared_entries(
            &vcf_entries(&vcf_reference, &record),
            &vcf_entries(&vcf_candidate, &record),
            &[UTR_WINDOW],
            None,
            None,
        );
        assert_eq!(outcome.green(), taken, "{at}");
        let outcome = compared_entries(
            &json_entries(&object(UTR_LOSS, 63), &record),
            &json_entries(&object(UTR_GAIN, at), &record),
            &[UTR_WINDOW],
            None,
            None,
        );
        assert_eq!(outcome.green(), taken, "{at}");
    }
}

/// Three insertions whose mutated sequence can also be written three bases nearer the cap,
/// each under a made-up identifier with its transcript's strand, and their input lines.
const EARLIER_ROWS: [(Row, &str); 3] = [
    (
        ("earlier_one", "21:2000-2001", "GTAATG", "ENST1", "498-499"),
        "1",
    ),
    (
        ("earlier_two", "21:3000-3001", "CTATG", "ENST1", "23-24"),
        "1",
    ),
    (
        (
            "earlier_three",
            "21:4000-4001",
            "CATAAC",
            "ENST1",
            "496-497",
        ),
        "-1",
    ),
];
const EARLIER_INPUTS: [&str; 3] = [
    "21\t2000\tearlier_one\tC\tCGTAATG\t.\t.\t.",
    "21\t3000\tearlier_two\tC\tCCTATG\t.\t.\t.",
    "21\t4000\tearlier_three\tC\tCCATAAC\t.\t.\t.",
];

/// What the reference and the candidate write on each row of `EARLIER_ROWS`: nothing, the
/// frameshift and the lost stop against a start gained three bases before the base after
/// the insertion.
fn earlier_pairs() -> [(Values, Values); 3] {
    [
        (
            nothing(),
            values(
                UTR_GAIN,
                "CapDistanceToStart=495:DistanceToCDS=500:DistanceToStop=105:KozakContext=GTTATGG:KozakStrength=Strong:type=uORF",
            ),
        ),
        (
            values(
                UTR_FRAMESHIFT,
                "Evidence=False:KozakContext=GTAATGG:KozakStrength=Strong:alt_type=uORF:alt_type_length=327:ref_StartDistanceToCDS=432:ref_type=uORF:ref_type_length=63",
            ),
            values(
                UTR_GAIN,
                "CapDistanceToStart=20:DistanceToCDS=437:DistanceToStop=324:KozakContext=GTAATGC:KozakStrength=Moderate:type=uORF",
            ),
        ),
        (
            values(
                UTR_STOP_LOSS,
                "AltStop=True:AltStopDistanceToCDS=299:Evidence=False:FrameWithCDS=outOfFrame:KozakContext=TTAATGG:KozakStrength=Moderate",
            ),
            values(
                UTR_GAIN,
                "CapDistanceToStart=493:DistanceToCDS=306:DistanceToStop=33:KozakContext=AAGATGG:KozakStrength=Strong:type=uORF",
            ),
        ),
    ]
}

/// The row of `EARLIER_ROWS` at `at` with the two fields `own`, on a transcript of
/// `strand` or of the row's own.
fn earlier_row(at: usize, own: &Values, strand: Option<&str>) -> String {
    let (row, own_strand) = EARLIER_ROWS[at];
    counted_row(row, own).replacen(
        "STRAND=1;",
        &format!("STRAND={};", strand.unwrap_or(own_strand)),
        1,
    )
}

/// The candidate names the first start codon of the changed bases of the first way of
/// writing the mutated sequence. For a record written three bases further on, that codon
/// begins three bases before the base after the insertion, where the allele, read on the
/// transcript's strand, spells it from the base three before its end: the row is the
/// class's whatever the reference writes but a gained start, counted once with both fields.
/// A three-base allele that spells the codon from its last base on to its first is the
/// class's at every index one allele's length apart before the insertion, on either strand,
/// beside the indices of its own window; without the input the two-base Location stays a
/// difference.
#[test]
fn the_gained_start_window_takes_a_start_codon_the_same_sequence_written_earlier_shows() {
    let pairs = earlier_pairs();
    let rows = |candidate: bool| -> String {
        let rows: Vec<String> = pairs
            .iter()
            .enumerate()
            .map(|(at, (reference, own))| {
                earlier_row(at, if candidate { own } else { reference }, None)
            })
            .collect();
        default_text(&rows)
    };
    let (reference, candidate) = (rows(false), rows(true));
    let outcome = window_with(&reference, &candidate, &EARLIER_INPUTS);
    assert!(outcome.is("green", &[(UTR_WINDOW, 3)]));
    assert!(outcome.masks("5UTR_consequence", &[(UTR_WINDOW, 3)]));
    assert!(outcome.masks("5UTR_annotation", &[(UTR_WINDOW, 3)]));
    let outcome = compared(&reference, &candidate, &[UTR_WINDOW]);
    assert!(outcome.is("red", &[]));
    assert_eq!(outcome.differing("5UTR_consequence"), 3);

    // The base after the insertion is index 146: `TGA` spells the codon from its last base.
    let inputs = |allele: &str| [format!("21\t1000\tround_ins\tG\tG{allele}\t.\t.\t.")];
    for (allele, strand) in [("TGA", "1"), ("TCA", "-1")] {
        let row: Row = ("round_ins", "21:1000-1001", allele, "ENST1", "146-147");
        let mut taken = Vec::new();
        for at in 131..=150 {
            let text = |own: &Values| one_row(stranded_row(row, own, Some(strand)));
            let outcome = compared_with(
                &text(&nothing()),
                &text(&gained(at)),
                &[UTR_WINDOW],
                Some(&inputs(allele)),
            );
            if outcome.green() {
                assert!(outcome.is("green", &[(UTR_WINDOW, 1)]), "{allele} {at}");
                taken.push(at);
            } else {
                assert!(outcome.is("red", &[]), "{allele} {at}");
            }
        }
        assert_eq!(taken, [133, 136, 139, 142, 146, 147, 148], "{allele}");
    }
}

/// A start the candidate alone reports before the insertion stays a difference where the
/// allele does not spell the codon at that offset, one or two bases before the insertion,
/// on a row read on the other strand or that states none, beside a reference that reports a
/// gained start of its own, and under a record whose reference allele is not empty, of one
/// base or of two. The same start reported by the reference alone stays a difference too.
#[test]
fn the_gained_start_window_refuses_a_start_before_the_insertion_the_allele_does_not_spell() {
    let pairs = earlier_pairs();
    let compared_at = |at: usize, reference: &Values, own: &Values, strand: Option<&str>| {
        let one = |own: &Values| one_row(earlier_row(at, own, strand));
        compared_with(
            &one(reference),
            &one(own),
            &[UTR_WINDOW],
            Some(&lines(&EARLIER_INPUTS)),
        )
    };
    // `GTAATG` before index 498: the codon at offset 3 stands 3 and 9 bases before its end.
    for (position, taken) in [
        (497, false),
        (496, false),
        (495, true),
        (494, false),
        (493, false),
        (492, false),
        (489, true),
        (488, false),
    ] {
        let outcome = compared_at(0, &nothing(), &gained(position), None);
        assert!(outcome.taken_is(taken, &[(UTR_WINDOW, 1)]), "{position}");
    }
    for (at, strand) in [(0, "-1"), (2, "1")] {
        let (reference, own) = &pairs[at];
        let outcome = compared_at(at, reference, own, Some(strand));
        assert!(outcome.is("red", &[]), "{at} {strand}");
    }
    let unstranded = |own: &Values| one_row(stranded_row(EARLIER_ROWS[0].0, own, None));
    let outcome = compared_with(
        &unstranded(&pairs[0].0),
        &unstranded(&pairs[0].1),
        &[UTR_WINDOW],
        Some(&lines(&EARLIER_INPUTS)),
    );
    assert!(outcome.is("red", &[]));
    let outcome = compared_at(0, &gained(498), &pairs[0].1, None);
    assert!(outcome.is("red", &[]));
    let replaced = lines(&["21\t2000\tearlier_one\tCA\tCGTAATG\t.\t.\t."]);
    let one = |own: &Values| one_row(earlier_row(0, own, None));
    let outcome = compared_with(
        &one(&pairs[0].0),
        &one(&pairs[0].1),
        &[UTR_WINDOW],
        Some(&replaced),
    );
    assert!(outcome.is("red", &[]));
    // Two reference bases, the first at index 497: the codon stands three bases before it.
    let longer = lines(&["21\t2000\tearlier_one\tCAG\tCGTAATG\t.\t.\t."]);
    let outcome = compared_with(
        &one(&pairs[0].0),
        &one(&gained(494)),
        &[UTR_WINDOW],
        Some(&longer),
    );
    assert!(outcome.is("red", &[]));
    let outcome = compared_at(0, &pairs[0].1, &nothing(), None);
    assert!(outcome.is("red", &[]));
}

/// The comparison of two one-row outputs under the lost-start class.
fn lost_after(reference: &str, candidate: &str) -> Outcome {
    compared(reference, candidate, &[UTR_LOST_AFTER])
}

/// An insertion between cDNA positions 146 and 147: Ensembl VEP reports the start codon at
/// 5' UTR index 146 lost, the codon that begins directly after the insertion and stands
/// whole. The row is the class's, counted once with its 5UTR_annotation, where the candidate
/// reports nothing, a frameshift, or a start gained inside its allele. It is read from the
/// row alone: under an identifier one input line carries, one that two lines carry and one
/// no line carries, and with the record's alleles spelt in the first column.
#[test]
fn the_lost_start_after_insertion_accepts_a_whole_codon_the_reference_reports_lost() {
    let rows: [(Row, i64); 3] = [
        (UTR_WHOLE, 146),
        (("whole_ins_b", "21:2000-2001", "C", "ENST2", "40-41"), 40),
        (
            ("whole_ins_c", "21:3000-3001", "CCATGC", "ENST3", "300-301"),
            300,
        ),
    ];
    let reference: Vec<String> = rows
        .iter()
        .map(|(row, first)| counted_row(*row, &lost_at(*first)))
        .collect();
    let reference = default_text(&reference);
    let candidate = default_text(&[
        counted_row(rows[0].0, &nothing()),
        counted_row(rows[1].0, &frame_shifted()),
        counted_row(rows[2].0, &gained(302)),
    ]);
    let outcome = lost_after(&reference, &candidate);
    assert!(outcome.is("green", &[(UTR_LOST_AFTER, 3)]));
    assert!(outcome.masks("5UTR_consequence", &[(UTR_LOST_AFTER, 3)]));
    assert!(outcome.masks("5UTR_annotation", &[(UTR_LOST_AFTER, 3)]));
    assert_eq!(compared(&reference, &candidate, &[]).verdict(), "red");
    let one_line = vec![UTR_WHOLE_INPUTS[0]];
    let two_lines = vec![UTR_WHOLE_INPUTS[0], "21\t1000\twhole_ins\tG\tGT\t.\t.\t."];
    for inputs in [&one_line, &two_lines, &UTR_WHOLE_INPUTS[1..].to_vec()] {
        let outcome = compared_with(
            &reference,
            &candidate,
            &[UTR_LOST_AFTER],
            Some(&lines(inputs)),
        );
        assert!(outcome.is("green", &[(UTR_LOST_AFTER, 3)]), "{inputs:?}");
    }
    let spelt = |text: &str| text.replace("whole_ins\t", "21_1000_G/GCCA\t");
    let outcome = lost_after(&spelt(&reference), &spelt(&candidate));
    assert!(outcome.is("green", &[(UTR_LOST_AFTER, 3)]));
}

/// The same insertion with three bases that end with neither the first base of a start
/// codon nor its first two.
const UTR_PARTED: Row = ("parted_ins", "21:1000-1001", "CCC", "ENST1", "146-147");

/// The same value pair stays a difference where the reported codon begins one or two
/// bases before the base after the insertion and the allele does not end with the codon's
/// first base or its first two, the insertion parting it, or past that base; where the
/// codon begins three bases before it; where the candidate reports the lost start, alone or
/// beside the reference with another annotation; where the reference's annotation names two
/// codons or none; where the cDNA positions are not an insertion's two flanking bases or are
/// not numbers; on a deletion's row, a one-base and a three-base substitution's and a
/// structural allele's; on a two-base substitution's and a longer replacement's whose
/// record's alleles the first column spells; and where the candidate's row carries none of
/// the plugin's counts.
#[test]
fn the_lost_start_after_insertion_refuses_a_codon_the_insertion_parts_and_every_other_allele() {
    let two = values(
        UTR_LOSS,
        "CapDistanceToStart=146:type=uORF,CapDistanceToStart=143:type=uORF",
    );
    let whole_at =
        |cdna: &'static str| -> Row { (UTR_WHOLE.0, UTR_WHOLE.1, UTR_WHOLE.2, UTR_WHOLE.3, cdna) };
    let cases: [(Row, Values, Values, &str); 23] = [
        (UTR_PARTED, lost_at(145), nothing(), UTR_COUNTS),
        (UTR_PARTED, lost_at(144), nothing(), UTR_COUNTS),
        (UTR_WHOLE, lost_at(144), nothing(), UTR_COUNTS),
        (UTR_WHOLE, lost_at(143), nothing(), UTR_COUNTS),
        (UTR_WHOLE, lost_at(147), nothing(), UTR_COUNTS),
        (UTR_PARTED, lost_at(145), gained(146), UTR_COUNTS),
        (
            UTR_WHOLE,
            lost_at(146),
            values(UTR_LOSS, "CapDistanceToStart=146:type=inframe_oORF"),
            UTR_COUNTS,
        ),
        (UTR_WHOLE, nothing(), lost_at(146), UTR_COUNTS),
        (UTR_WHOLE, frame_shifted(), lost_at(146), UTR_COUNTS),
        (UTR_WHOLE, two, nothing(), UTR_COUNTS),
        (
            UTR_WHOLE,
            values(UTR_LOSS, "type=uORF"),
            nothing(),
            UTR_COUNTS,
        ),
        (UTR_WHOLE, gained(146), nothing(), UTR_COUNTS),
        (whole_at("146"), lost_at(146), nothing(), UTR_COUNTS),
        (whole_at("146-148"), lost_at(146), nothing(), UTR_COUNTS),
        (whole_at("?-147"), lost_at(146), nothing(), UTR_COUNTS),
        (whole_at("146-?"), lost_at(146), nothing(), UTR_COUNTS),
        (
            ("del", "21:1000-1001", "-", "ENST1", "146-147"),
            lost_at(146),
            nothing(),
            UTR_COUNTS,
        ),
        (
            ("snv", "21:600", "T", "ENST1", "146"),
            lost_at(146),
            nothing(),
            UTR_COUNTS,
        ),
        (
            ("21_500_ACG/GTA", "21:500-502", "GTA", "ENST1", "146-148"),
            lost_at(146),
            nothing(),
            UTR_COUNTS,
        ),
        (
            ("sv", "21:1000-1001", "<INS>", "ENST1", "146-147"),
            lost_at(146),
            nothing(),
            UTR_COUNTS,
        ),
        (
            ("21_500_AC/GT", "21:500-501", "GT", "ENST1", "146-147"),
            lost_at(146),
            nothing(),
            UTR_COUNTS,
        ),
        (
            ("21_700_AC/GTTT", "21:700-701", "GTTT", "ENST1", "146-147"),
            lost_at(146),
            nothing(),
            UTR_COUNTS,
        ),
        (UTR_WHOLE, lost_at(146), nothing(), ""),
    ];
    for (row, reference, candidate, counts) in &cases {
        let outcome = lost_after(
            &one_row(counted_row(*row, reference)),
            &one_row(counted_row_with(*row, candidate, counts)),
        );
        assert!(
            outcome.is("red", &[]),
            "{row:?} {reference:?} {candidate:?} {counts:?}"
        );
    }
}

/// A made-up grid over an insertion between cDNA positions 146 and 147, the base after it
/// 5' UTR index 146. The reference reports a start codon lost at index 146, 145, 144, 143
/// or 147; the allele, read on the transcript's strand, is one to six bases long and ends
/// with a start codon's first base, with its first two, or with neither; the transcript is
/// on the forward strand, on the reverse strand, where the row's Allele is the reverse
/// complement of what the transcript reads, or the row states no strand; the candidate
/// reports nothing, a frameshift, a lost stop or a start gained inside its allele. The row
/// is the class's exactly where the codon stands whole after the insertion: reported at
/// 146, whatever the allele and the strand; at 145 where the row states its strand and the
/// allele so read ends with the codon's first base; at 144 where it ends with the codon's
/// first two. It is counted once with its 5UTR_annotation.
#[test]
fn the_lost_start_after_insertion_takes_a_codon_the_insertion_only_moves_on() {
    let stop_lost = values(
        UTR_STOP_LOSS,
        "AltStop=False:FrameWithCDS=inFrame:KozakStrength=Weak",
    );
    let read_alleles = [
        "A", "T", "AT", "CA", "CC", "CCA", "CAT", "CCC", "GCCA", "GCAT", "GCCC", "CCGCCA",
        "CCGCAT", "ATA", "TAT",
    ];
    let (mut taken, mut refused) = (0, 0);
    for lost_index in [146, 145, 144, 143, 147] {
        for read in read_alleles {
            for strand in [Some("1"), Some("-1"), None] {
                let written = if strand == Some("-1") {
                    reverse_complement(read)
                } else {
                    read.to_string()
                };
                let row: Row = ("grid_ins", "21:1000-1001", &written, "ENST1", "146-147");
                // Without a strand the allele is read at the codon's own index alone.
                let on_strand = strand.map(|_| read);
                let whole = lost_index == 146
                    || (lost_index == 145 && on_strand.is_some_and(|read| read.ends_with('A')))
                    || (lost_index == 144 && on_strand.is_some_and(|read| read.ends_with("AT")));
                let candidates = [
                    nothing(),
                    frame_shifted(),
                    stop_lost.clone(),
                    gained(146 + read.len() as i64 - 1),
                ];
                for candidate in &candidates {
                    let outcome = lost_after(
                        &one_row(stranded_row(row, &lost_at(lost_index), strand)),
                        &one_row(stranded_row(row, candidate, strand)),
                    );
                    let case = format!("{lost_index} {read} {strand:?} {}", candidate.0);
                    if whole {
                        assert!(outcome.is("green", &[(UTR_LOST_AFTER, 1)]), "{case}");
                        assert!(
                            outcome.masks("5UTR_consequence", &[(UTR_LOST_AFTER, 1)]),
                            "{case}"
                        );
                        assert!(
                            outcome.masks("5UTR_annotation", &[(UTR_LOST_AFTER, 1)]),
                            "{case}"
                        );
                        taken += 1;
                    } else {
                        assert!(outcome.is("red", &[]), "{case}");
                        refused += 1;
                    }
                }
            }
        }
    }
    // 15 alleles, 3 strands and 4 candidate values at each index: every case at 146; at 145 the
    // 6 alleles that end with the first base on the 2 stated strands; at 144 the 5 that end
    // with the first two.
    let whole = 15 * 3 * 4 + 6 * 2 * 4 + 5 * 2 * 4;
    assert_eq!((taken, refused), (whole, 5 * 15 * 3 * 4 - whole));
}

/// The fields of the class's row that it does not read stay differences beside it.
const OTHER_FIELDS: [&str; 6] = [
    "Existing_uORFs",
    "Existing_InFrame_oORFs",
    "Existing_OutOfFrame_oORFs",
    "LOEUF",
    "Consequence",
    "Extra",
];

/// At a codon reported one or two bases before the base after the insertion, with an
/// allele that ends with its first base or its first two, the row stays a difference where
/// the candidate reports the lost start, alone or beside the reference; on a deletion's
/// row, a one-base substitution's and a structural allele's; on a two-base substitution's
/// and a longer replacement's whose record's alleles the first column spells or the input
/// states; and where the candidate's row carries none of the plugin's counts. A two-base
/// Location under an identifier with no record handed is read by its positions and its
/// allele alone. The plugin's counts and another plugin's field on the class's row stay
/// differences, with and without the seeded runs.
#[test]
fn the_lost_start_after_insertion_refuses_the_side_the_kind_and_the_field_at_a_moved_codon() {
    let moved: Row = ("moved_ins", "21:1000-1001", "CCA", "ENST1", "146-147");
    let moved_two: Row = ("moved_ins", "21:1000-1001", "CAT", "ENST1", "146-147");
    for (row, lost_index) in [(moved, 145), (moved_two, 144)] {
        let cases: [(Values, Values, &str); 4] = [
            (nothing(), lost_at(lost_index), UTR_COUNTS),
            (frame_shifted(), lost_at(lost_index), UTR_COUNTS),
            (
                lost_at(lost_index),
                values(
                    UTR_LOSS,
                    &format!("CapDistanceToStart={lost_index}:type=inframe_oORF"),
                ),
                UTR_COUNTS,
            ),
            (lost_at(lost_index), nothing(), ""),
        ];
        for (reference, candidate, counts) in &cases {
            let outcome = lost_after(
                &one_row(counted_row(row, reference)),
                &one_row(counted_row_with(row, candidate, counts)),
            );
            assert!(
                outcome.is("red", &[]),
                "{row:?} {reference:?} {candidate:?} {counts:?}"
            );
        }
    }
    let others: [Row; 5] = [
        ("del", "21:1000-1001", "-", "ENST1", "146-147"),
        ("snv", "21:600", "A", "ENST1", "146"),
        ("sv", "21:1000-1001", "<INS>", "ENST1", "146-147"),
        ("21_500_AC/CA", "21:500-501", "CA", "ENST1", "146-147"),
        ("21_700_AC/GTCA", "21:700-701", "GTCA", "ENST1", "146-147"),
    ];
    for row in others {
        let outcome = lost_after(
            &one_row(counted_row(row, &lost_at(145))),
            &one_row(counted_row(row, &nothing())),
        );
        assert!(outcome.is("red", &[]), "{row:?}");
    }
    // Under an identifier the input's line says whether the reference allele is empty.
    let inputs = lines(&[
        "21\t1000\tmoved_ins\tG\tGCCA\t.\t.\t.",
        "21\t500\tmnv\tAC\tCA\t.\t.\t.",
        "21\t700\tlonger\tAC\tGTCA\t.\t.\t.",
    ]);
    let rows: [(Row, bool, bool); 3] = [
        (moved, true, true),
        (("mnv", "21:500-501", "CA", "ENST1", "146-147"), false, true),
        (
            ("longer", "21:700-701", "GTCA", "ENST1", "146-147"),
            false,
            true,
        ),
    ];
    let texts = |own: &Values| -> Vec<String> {
        rows.iter()
            .map(|(row, _, _)| counted_row(*row, own))
            .collect()
    };
    let reference = default_entries(&default_text(&texts(&lost_at(145))));
    let candidate = default_entries(&default_text(&texts(&nothing())));
    for (at, (row, handed, alone)) in rows.iter().enumerate() {
        let found = accepts_with_and_without_the_record(
            UTR_LOST_AFTER,
            &reference,
            &inputs,
            at,
            &candidate[at],
            UTR_LOSS,
            "-",
        );
        assert_eq!(found, [*handed, *handed, *alone, *alone], "{row:?}");
    }
    // An identifier that one line, two lines and no line of the input carry.
    let pair = (
        one_row(counted_row(moved, &lost_at(145))),
        one_row(counted_row(moved, &nothing())),
    );
    let second = inputs[0].replace("GCCA", "GT");
    let given: [Option<Vec<String>>; 4] = [
        Some(vec![inputs[0].clone()]),
        Some(vec![inputs[0].clone(), second]),
        Some(inputs[1..].to_vec()),
        None,
    ];
    for lines in &given {
        let outcome = compared_with(&pair.0, &pair.1, &[UTR_LOST_AFTER], lines.as_deref());
        assert!(outcome.is("green", &[(UTR_LOST_AFTER, 1)]), "{lines:?}");
    }
    let (e, a) = (&reference[0], &candidate[0]);
    for name in OTHER_FIELDS {
        assert!(
            !accepts(UTR_LOST_AFTER, e, name, "1", "2", Some(a), None, None),
            "{name}"
        );
    }
    let other = counted_row_with(
        moved,
        &nothing(),
        &UTR_COUNTS.replace("Existing_uORFs=1", "Existing_uORFs=2"),
    ) + ";LOEUF=0.4";
    let lost = one_row(counted_row(moved, &lost_at(145)) + ";LOEUF=0.3");
    for swept in [false, true] {
        let outcome = if swept {
            compared_swept(
                &lost,
                &one_row(other.clone()),
                Some(&[&lost, &lost]),
                &[UTR_LOST_AFTER, UTR_PRECEDENCE],
            )
        } else {
            lost_after(&lost, &one_row(other.clone()))
        };
        assert!(outcome.is("red", &[(UTR_LOST_AFTER, 1)]));
        assert_eq!(outcome.differing("Existing_uORFs"), 1);
        assert_eq!(outcome.differing("LOEUF"), 1);
        assert!(outcome.masks("5UTR_consequence", &[(UTR_LOST_AFTER, 1)]));
    }
}

/// A row two classes hold at a codon the insertion stands inside is counted under the
/// first in the order a row is put to the classes. Five bases after the first base of the
/// codon at index 145 that hold a start codon and end with the codon's first base: the
/// candidate's gained start at 147 is the window class's where the record's alleles are
/// read, and this class's where they are not; a gained start at 150, the lost codon's index
/// plus the allele's length, names the moved codon, which the window class refuses and this
/// class takes. A frameshift the seeded runs write on the row beside the lost start is the
/// precedence class's, and this class's without them. A lost start both engines report,
/// the reference's annotation in another order, is the pair-order class's alone.
#[test]
fn the_lost_start_after_insertion_shares_a_moved_codon_with_the_plugins_other_classes() {
    let making: Row = ("moved_ins", "21:1000-1001", "CATGA", "ENST1", "146-147");
    let inputs = lines(&["21\t1000\tmoved_ins\tG\tGCATGA\t.\t.\t."]);
    let two_lines = vec![inputs[0].clone(), inputs[0].replace("GCATGA", "GT")];
    let reference = one_row(counted_row(making, &lost_at(145)));
    let cases: [(i64, Inputs, &str); 5] = [
        (147, Some(&inputs), UTR_WINDOW),
        (147, Some(&two_lines), UTR_WINDOW),
        (147, None, UTR_LOST_AFTER),
        (150, Some(&inputs), UTR_LOST_AFTER),
        (150, None, UTR_LOST_AFTER),
    ];
    for (at, lines, credited) in cases {
        for classes in [[UTR_WINDOW, UTR_LOST_AFTER], [UTR_LOST_AFTER, UTR_WINDOW]] {
            let outcome = compared_with(
                &reference,
                &one_row(counted_row(making, &gained(at))),
                &classes,
                lines,
            );
            assert!(
                outcome.is("green", &[(credited, 1)]),
                "{at} {lines:?} {classes:?}"
            );
        }
    }
    let moved: Row = ("moved_ins", "21:1000-1001", "A", "ENST1", "146-147");
    let lost = one_row(counted_row(moved, &lost_at(145)));
    let shifted = one_row(counted_row(moved, &frame_shifted()));
    let nothing_written = one_row(counted_row(moved, &nothing()));
    let both_terms = [lost.as_str(), shifted.as_str()];
    let lost_alone = [lost.as_str(), lost.as_str()];
    let cases: [(&str, Seeded, Classes, &str); 6] = [
        (
            &shifted,
            Some(&both_terms),
            &[UTR_PRECEDENCE, UTR_LOST_AFTER],
            UTR_PRECEDENCE,
        ),
        (
            &shifted,
            Some(&both_terms),
            &[UTR_LOST_AFTER, UTR_PRECEDENCE],
            UTR_PRECEDENCE,
        ),
        (
            &shifted,
            Some(&lost_alone),
            &[UTR_PRECEDENCE, UTR_LOST_AFTER],
            UTR_LOST_AFTER,
        ),
        (
            &shifted,
            None,
            &[UTR_PRECEDENCE, UTR_LOST_AFTER],
            UTR_LOST_AFTER,
        ),
        (
            &shifted,
            Some(&both_terms),
            &[UTR_LOST_AFTER],
            UTR_LOST_AFTER,
        ),
        (
            &nothing_written,
            Some(&both_terms),
            &[UTR_PRECEDENCE, UTR_LOST_AFTER],
            UTR_LOST_AFTER,
        ),
    ];
    for (candidate, sweep, classes, credited) in cases {
        let outcome = compared_swept(&lost, candidate, sweep, classes);
        assert!(
            outcome.is("green", &[(credited, 1)]),
            "{classes:?} {}",
            sweep.is_none()
        );
    }
    let hash_order = one_row(counted_row(
        moved,
        &values(
            UTR_LOSS,
            "type=uORF:KozakStrength=Moderate:DistanceToCDS=110:CapDistanceToStart=145",
        ),
    ));
    let outcome = compared(&hash_order, &lost, &[PAIR_ORDER, UTR_LOST_AFTER]);
    assert!(outcome.is("green", &[(PAIR_ORDER, 1)]));
}

/// The tab format names the strand as a column, the VCF format as a CSQ field and the
/// JSON format as the consequence object's `strand`: on a reverse-strand transcript a
/// codon reported one below the lower cDNA position is the class's where the row's allele,
/// reverse-complemented, ends with the codon's first base, and stays a difference where it
/// does not or the format carries no strand.
#[test]
fn the_lost_start_after_insertion_reads_a_moved_codon_in_the_tab_vcf_and_json_formats() {
    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tSTRAND\t5UTR_consequence\t5UTR_annotation\tExisting_uORFs\n";
    let tab = |allele: &str, strand: &str, own: &Values| {
        format!(
            "{header}moved_ins\t21:1000-1001\t{allele}\tENSG1\tENST1\tTranscript\t5_prime_UTR_variant\t146-147\t{strand}\t{}\t{}\t1\n",
            own.0, own.1
        )
    };
    for (allele, strand, taken) in [
        ("TGG", "-1", true),
        ("CCA", "1", true),
        ("CCA", "-1", false),
        ("TGG", "1", false),
        ("TGG", "-", false),
    ] {
        let outcome = compared_entries(
            &tab_entries(&tab(allele, strand, &lost_at(145))),
            &tab_entries(&tab(allele, strand, &nothing())),
            &[UTR_LOST_AFTER],
            None,
            None,
        );
        assert!(
            outcome.taken_is(taken, &[(UTR_LOST_AFTER, 1)]),
            "{allele} {strand}"
        );
    }

    let inputs = lines(&["21\t1000\tmoved_ins\tC\tCTGG\t.\t.\t."]);
    let lost = format!("{UTR_LOSS}|CapDistanceToStart=145:type=uORF");
    for (fields, strand, taken) in [
        ("STRAND|", "-1|", true),
        ("STRAND|", "1|", false),
        ("", "", false),
    ] {
        let csq = csq_header(&format!(
            "Allele|Consequence|Feature|Feature_type|cDNA_position|{fields}5UTR_consequence|5UTR_annotation|Existing_uORFs"
        ));
        let vcf_reference = format!(
            "{csq}{}CSQ=TGG|5_prime_UTR_variant|ENST1|Transcript|146-147|{strand}{lost}|1\n",
            record_columns(&inputs[0])
        );
        let outcome = compared_entries(
            &vcf_entries(&vcf_reference, &inputs),
            &vcf_entries(&vcf_reference.replace(&lost, "|"), &inputs),
            &[UTR_LOST_AFTER],
            None,
            None,
        );
        assert!(
            outcome.taken_is(taken, &[(UTR_LOST_AFTER, 1)]),
            "{fields} {strand}"
        );
    }

    let object = |strand: Option<i64>, lost_start: bool| {
        let mut entry = r#""variant_allele": "TGG", "transcript_id": "ENST1", "consequence_terms": ["5_prime_UTR_variant"], "cdna_start": 147, "cdna_end": 146, "existing_uorfs": 1"#.to_string();
        if let Some(strand) = strand {
            entry.push_str(&format!(r#", "strand": {strand}"#));
        }
        if lost_start {
            entry.push_str(&format!(
                r#", "5utr_consequence": "{UTR_LOSS}", "5utr_annotation": {{"1": {{"CapDistanceToStart": 145, "type": "uORF"}}}}"#
            ));
        }
        json_line(&inputs[0], "r0", "-/TGG", &entry)
    };
    for (strand, taken) in [(Some(-1), true), (Some(1), false), (None, false)] {
        let outcome = compared_entries(
            &json_entries(&object(strand, true), &inputs),
            &json_entries(&object(strand, false), &inputs),
            &[UTR_LOST_AFTER],
            None,
            None,
        );
        assert!(
            outcome.taken_is(taken, &[(UTR_LOST_AFTER, 1)]),
            "{strand:?}"
        );
    }
}

/// A row under an identifier spells no alleles. Handed its record, the class refuses a
/// two-base substitution's row and a longer replacement's, whose reference allele is not
/// empty, and takes the insertion's; without one it reads the positions, whose equality with
/// the codon's index only an insertion's row carries: with a reference allele
/// UTRAnnotator.pm:769 keeps the codon's index at or below the allele's own, the lower cDNA
/// position less one. The plugin's three counts and another plugin's field on the class's
/// row stay differences.
#[test]
fn the_lost_start_after_insertion_reads_a_record_it_is_handed_and_the_positions_without_one() {
    let rows: [(Row, bool, bool); 3] = [
        (UTR_WHOLE, true, true),
        (("mnv", "21:500-501", "GT", "ENST1", "146-147"), false, true),
        (
            ("longer", "21:700-701", "GTTT", "ENST1", "146-147"),
            false,
            true,
        ),
    ];
    let texts = |own: &Values| -> Vec<String> {
        rows.iter()
            .map(|(row, _, _)| counted_row(*row, own))
            .collect()
    };
    let reference = default_entries(&default_text(&texts(&lost_at(146))));
    let candidate = default_entries(&default_text(&texts(&nothing())));
    let inputs = lines(&UTR_WHOLE_INPUTS);
    for (at, (row, handed, alone)) in rows.iter().enumerate() {
        let found = accepts_with_and_without_the_record(
            UTR_LOST_AFTER,
            &reference,
            &inputs,
            at,
            &candidate[at],
            UTR_LOSS,
            "-",
        );
        assert_eq!(found, [*handed, *handed, *alone, *alone], "{row:?}");
    }
    let (e, a) = (&reference[0], &candidate[0]);
    for name in OTHER_FIELDS {
        assert!(
            !accepts(UTR_LOST_AFTER, e, name, "1", "2", Some(a), None, None),
            "{name}"
        );
    }
    assert!(!accepts(
        UTR_LOST_AFTER,
        e,
        "5UTR_consequence",
        UTR_LOSS,
        "-",
        None,
        None,
        None
    ));
    let moved = counted_row_with(
        UTR_WHOLE,
        &nothing(),
        &UTR_COUNTS.replace("Existing_uORFs=1", "Existing_uORFs=2"),
    ) + ";LOEUF=0.4";
    let outcome = lost_after(
        &one_row(counted_row(UTR_WHOLE, &lost_at(146)) + ";LOEUF=0.3"),
        &one_row(moved),
    );
    assert!(outcome.is("red", &[(UTR_LOST_AFTER, 1)]));
    assert_eq!(outcome.differing("Existing_uORFs"), 1);
    assert_eq!(outcome.differing("LOEUF"), 1);
    assert!(outcome.masks("5UTR_consequence", &[(UTR_LOST_AFTER, 1)]));
}

/// A row two classes hold is counted under the first in the order a row is put to the
/// classes: a start the candidate gains inside its allele is the window class's where the
/// record's alleles are read, and a frameshift the seeded runs write on the row beside the
/// lost start is the precedence class's. Named alone, or where the other class does not hold
/// (no record of the row, no seeded run of it, a candidate without a term), the row is this
/// class's.
#[test]
fn the_lost_start_after_insertion_follows_the_plugins_other_two_classes() {
    assert_eq!(
        plugin_classes(),
        [
            PAIR_ORDER,
            UTR_WINDOW,
            UTR_PRECEDENCE,
            UTR_STRUCTURAL,
            UTR_LOST_AFTER,
            UTR_REMADE,
            UTR_BOUNDARY,
            UTR_BEFORE_STOP,
            UTR_WITHOUT_REFERENCE
        ]
    );
    let gaining: Row = ("whole_ins", "21:1000-1001", "CCATGC", "ENST1", "146-147");
    let reference = one_row(counted_row(gaining, &lost_at(146)));
    let candidate = one_row(counted_row(gaining, &gained(148)));
    let inputs = lines(&["21\t1000\twhole_ins\tG\tGCCATGC\t.\t.\t."]);
    let two_lines = vec![inputs[0].clone(), inputs[0].replace("GCCATGC", "GT")];
    let cases: [(Classes, Inputs, &str); 5] = [
        (&[UTR_WINDOW, UTR_LOST_AFTER], Some(&inputs), UTR_WINDOW),
        (&[UTR_LOST_AFTER, UTR_WINDOW], Some(&inputs), UTR_WINDOW),
        (&[UTR_WINDOW, UTR_LOST_AFTER], Some(&two_lines), UTR_WINDOW),
        (&[UTR_WINDOW, UTR_LOST_AFTER], None, UTR_LOST_AFTER),
        (&[UTR_LOST_AFTER], Some(&inputs), UTR_LOST_AFTER),
    ];
    for (classes, lines, credited) in cases {
        let outcome = compared_with(&reference, &candidate, classes, lines);
        assert!(
            outcome.is("green", &[(credited, 1)]),
            "{classes:?} {lines:?}"
        );
    }
    let shifted = one_row(counted_row(UTR_WHOLE, &frame_shifted()));
    let lost = one_row(counted_row(UTR_WHOLE, &lost_at(146)));
    let nothing_written = one_row(counted_row(UTR_WHOLE, &nothing()));
    let both_terms = [lost.as_str(), shifted.as_str()];
    let lost_alone = [lost.as_str(), lost.as_str()];
    let cases: [(&str, Seeded, Classes, &str); 7] = [
        (
            &shifted,
            Some(&both_terms),
            &[UTR_PRECEDENCE, UTR_LOST_AFTER],
            UTR_PRECEDENCE,
        ),
        (
            &shifted,
            Some(&both_terms),
            &[UTR_LOST_AFTER, UTR_PRECEDENCE],
            UTR_PRECEDENCE,
        ),
        (
            &shifted,
            Some(&lost_alone),
            &[UTR_PRECEDENCE, UTR_LOST_AFTER],
            UTR_LOST_AFTER,
        ),
        (
            &shifted,
            None,
            &[UTR_PRECEDENCE, UTR_LOST_AFTER],
            UTR_LOST_AFTER,
        ),
        (
            &shifted,
            Some(&both_terms),
            &[UTR_LOST_AFTER],
            UTR_LOST_AFTER,
        ),
        (
            &nothing_written,
            Some(&both_terms),
            &[UTR_PRECEDENCE, UTR_LOST_AFTER],
            UTR_LOST_AFTER,
        ),
        (&nothing_written, None, &[UTR_LOST_AFTER], UTR_LOST_AFTER),
    ];
    for (candidate, sweep, classes, credited) in cases {
        let outcome = compared_swept(&lost, candidate, sweep, classes);
        assert!(
            outcome.is("green", &[(credited, 1)]),
            "{classes:?} {}",
            sweep.is_none()
        );
    }
}

/// The command line names the lost-start class like any other.
#[test]
fn the_command_line_names_the_lost_start_after_insertion_class() {
    let lost = one_row(counted_row(UTR_WHOLE, &lost_at(146)));
    let nothing_written = one_row(counted_row(UTR_WHOLE, &nothing()));
    let directory = planted_files(
        "utrannotator-lost-start",
        &[("reference", &lost), ("candidate", &nothing_written)],
    );
    assert_eq!(exit_code(&directory, &[]), 1);
    assert_eq!(exit_code(&directory, &["--masks", UTR_LOST_AFTER]), 0);
}

/// The tab format names the two fields and the count as columns. The VCF and JSON
/// entries carry their record's alleles, so an insertion's row is the class's and the same
/// values on a two-base substitution's record stay a difference; a row without a term
/// writes two empty VCF fields and no JSON key. The JSON format writes an insertion's
/// cdna_start past its cdna_end and the annotation as an object.
#[test]
fn the_lost_start_after_insertion_reads_the_tab_vcf_and_json_formats() {
    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\t5UTR_consequence\t5UTR_annotation\tExisting_uORFs\n";
    let tab = |own: &Values| {
        format!(
            "{header}whole_ins\t21:1000-1001\tCCA\tENSG1\tENST1\tTranscript\t5_prime_UTR_variant\t146-147\t{}\t{}\t1\n",
            own.0, own.1
        )
    };
    let under_class = |reference: &[Entry], candidate: &[Entry]| {
        compared_entries(reference, candidate, &[UTR_LOST_AFTER], None, None)
    };
    let outcome = under_class(
        &tab_entries(&tab(&lost_at(146))),
        &tab_entries(&tab(&nothing())),
    );
    assert!(outcome.is("green", &[(UTR_LOST_AFTER, 1)]));
    let outcome = under_class(
        &tab_entries(&tab(&lost_at(144))),
        &tab_entries(&tab(&nothing())),
    );
    assert_eq!(outcome.verdict(), "red");

    let inputs = lines(&UTR_WHOLE_INPUTS[..2]);
    let csq = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|5UTR_consequence|5UTR_annotation|Existing_uORFs",
    );
    let lost = format!("{UTR_LOSS}|CapDistanceToStart=146:type=uORF");
    let vcf_reference = format!(
        "{csq}{}CSQ=CCA|5_prime_UTR_variant|ENST1|Transcript|146-147|{lost}|1\n{}CSQ=GT|5_prime_UTR_variant|ENST3|Transcript|146-147|{lost}|1\n",
        record_columns(&inputs[0]),
        record_columns(&inputs[1])
    );
    let outcome = under_class(
        &vcf_entries(&vcf_reference, &inputs),
        &vcf_entries(&vcf_reference.replace(&lost, "|"), &inputs),
    );
    assert!(outcome.is("red", &[(UTR_LOST_AFTER, 1)]));
    assert_eq!(outcome.differing("5UTR_consequence"), 1);

    let object = |at: usize,
                  allele_string: &str,
                  allele: &str,
                  feature: &str,
                  start: i64,
                  end: i64,
                  lost_start: bool| {
        let mut entry = format!(
            r#""variant_allele": "{allele}", "transcript_id": "{feature}", "consequence_terms": ["5_prime_UTR_variant"], "cdna_start": {start}, "cdna_end": {end}, "existing_uorfs": 1"#
        );
        if lost_start {
            entry.push_str(&format!(
                r#", "5utr_consequence": "{UTR_LOSS}", "5utr_annotation": {{"1": {{"CapDistanceToStart": 146, "type": "uORF"}}}}"#
            ));
        }
        json_line(&inputs[at], &format!("r{at}"), allele_string, &entry)
    };
    let both = |lost_start: bool| {
        object(0, "-/CCA", "CCA", "ENST1", 147, 146, lost_start)
            + &object(1, "AC/GT", "GT", "ENST3", 146, 147, lost_start)
    };
    let outcome = under_class(
        &json_entries(&both(true), &inputs),
        &json_entries(&both(false), &inputs),
    );
    assert!(outcome.is("red", &[(UTR_LOST_AFTER, 1)]));
    assert_eq!(outcome.differing("5utr_consequence"), 1);
}

/// The comparison of two one-row outputs under the remade-in-place class.
fn remade(reference: &str, candidate: &str) -> Outcome {
    compared(reference, candidate, &[UTR_REMADE])
}

/// A made-up grid over an insertion between cDNA positions 146 and 147, the base after it
/// 5' UTR index 146. The reference reports a start codon gained at index 146, 145, 144, 143
/// or 147; the allele, read on the transcript's strand, is one to six bases long and opens
/// with a start codon's last two bases, with its last base, or with neither; the transcript
/// is on the forward strand, on the reverse strand, where the row's Allele is the reverse
/// complement of what the transcript reads, or the row states no strand; the candidate
/// reports nothing, a frameshift, a lost stop, a start gained three bases further on, or the
/// same gained start. The row is the class's exactly where the inserted bases leave the
/// reference's codon at its own index and the candidate does not report it gained: at 145
/// where the row states its strand and the allele so read opens with the codon's last two
/// bases, at 144 where it opens with its last. It is counted once with its 5UTR_annotation.
#[test]
fn the_gained_start_remade_in_place_takes_a_codon_the_insertion_leaves_at_its_own_index() {
    let stop_lost = values(
        UTR_STOP_LOSS,
        "AltStop=False:FrameWithCDS=inFrame:KozakStrength=Weak",
    );
    let read_alleles = [
        "G", "T", "TG", "GC", "CC", "TGC", "GCC", "CCC", "TGCA", "GTGC", "CCCC", "TGCCGC",
        "GCCGCC", "GT", "TGG",
    ];
    let (mut taken, mut refused) = (0, 0);
    for gained_index in [146, 145, 144, 143, 147] {
        for read in read_alleles {
            for strand in [Some("1"), Some("-1"), None] {
                let written = if strand == Some("-1") {
                    reverse_complement(read)
                } else {
                    read.to_string()
                };
                let row: Row = ("grid_ins", "21:1000-1001", &written, "ENST1", "146-147");
                let on_strand = strand.map(|_| read);
                let in_place = on_strand.is_some_and(|read| {
                    (gained_index == 145 && read.starts_with("TG"))
                        || (gained_index == 144 && read.starts_with('G'))
                });
                let candidates = [
                    nothing(),
                    frame_shifted(),
                    stop_lost.clone(),
                    gained(gained_index + 3),
                    gained(gained_index),
                ];
                for candidate in &candidates {
                    let outcome = remade(
                        &one_row(stranded_row(row, &gained(gained_index), strand)),
                        &one_row(stranded_row(row, candidate, strand)),
                    );
                    let case = format!("{gained_index} {read} {strand:?} {candidate:?}");
                    if *candidate == gained(gained_index) {
                        assert!(outcome.is("green", &[]), "{case}");
                    } else if in_place {
                        assert!(outcome.is("green", &[(UTR_REMADE, 1)]), "{case}");
                        if candidate.0 != UTR_GAIN {
                            assert!(
                                outcome.masks("5UTR_consequence", &[(UTR_REMADE, 1)]),
                                "{case}"
                            );
                        }
                        assert!(
                            outcome.masks("5UTR_annotation", &[(UTR_REMADE, 1)]),
                            "{case}"
                        );
                        taken += 1;
                    } else {
                        assert!(outcome.is("red", &[]), "{case}");
                        refused += 1;
                    }
                }
            }
        }
    }
    // 15 alleles, 3 strands and 4 differing candidate values at each index: at 145 the 5
    // alleles that open with the last two bases on the 2 stated strands; at 144 the 6 that
    // open with the last base.
    let in_place = 5 * 2 * 4 + 6 * 2 * 4;
    assert_eq!((taken, refused), (in_place, 5 * 15 * 3 * 4 - in_place));
}

/// At a start codon reported one below the lower cDNA position, with an allele that opens
/// with the codon's last two bases, the row stays a difference where the candidate alone
/// reports the gained start; where the reference's annotation names two codons or none;
/// where the cDNA positions are not an insertion's two flanking bases; on a deletion's row,
/// a one-base substitution's and a structural allele's; on a two-base substitution's and a
/// longer replacement's whose record's alleles the first column spells or the input states;
/// and where the candidate's row carries none of the plugin's counts; and where the gained
/// start stands at a negative index. The plugin's counts and another plugin's field on the
/// class's row stay differences.
#[test]
fn the_gained_start_remade_in_place_refuses_the_side_the_kind_and_the_field() {
    let kept: Row = ("kept_ins", "21:1000-1001", "TGC", "ENST1", "146-147");
    let two = values(
        UTR_GAIN,
        "CapDistanceToStart=145:type=uORF,CapDistanceToStart=143:type=uORF",
    );
    // A gained start the reference's wrapped window reports before the UTR's first base, two
    // below the lower cDNA position of an insertion at the UTR's second base, is the window
    // class's and none of this class's.
    let wrapped: Row = ("wrapped_ins", "21:1000-1001", "GCC", "ENST1", "1-2");
    let pair = (
        one_row(counted_row(wrapped, &gained(-1))),
        one_row(counted_row(wrapped, &nothing())),
    );
    assert!(remade(&pair.0, &pair.1).is("red", &[]));
    let outcome = compared(&pair.0, &pair.1, &[UTR_WINDOW, UTR_REMADE]);
    assert!(outcome.is("green", &[(UTR_WINDOW, 1)]));
    let kept_at = |cdna: &'static str| -> Row { (kept.0, kept.1, kept.2, kept.3, cdna) };
    let cases: [(Row, Values, Values, &str); 15] = [
        (kept, nothing(), gained(145), UTR_COUNTS),
        (kept, frame_shifted(), gained(145), UTR_COUNTS),
        (kept, lost_at(145), nothing(), UTR_COUNTS),
        (kept, two, nothing(), UTR_COUNTS),
        (kept, values(UTR_GAIN, "type=uORF"), nothing(), UTR_COUNTS),
        (kept, gained(145), values(UTR_GAIN, "type=uORF"), UTR_COUNTS),
        (kept_at("146"), gained(145), nothing(), UTR_COUNTS),
        (kept_at("146-148"), gained(145), nothing(), UTR_COUNTS),
        (kept_at("?-147"), gained(145), nothing(), UTR_COUNTS),
        (
            ("del", "21:1000-1001", "-", "ENST1", "146-147"),
            gained(145),
            nothing(),
            UTR_COUNTS,
        ),
        (
            ("snv", "21:600", "T", "ENST1", "146"),
            gained(145),
            nothing(),
            UTR_COUNTS,
        ),
        (
            ("sv", "21:1000-1001", "<INS>", "ENST1", "146-147"),
            gained(145),
            nothing(),
            UTR_COUNTS,
        ),
        (
            ("21_500_AC/TG", "21:500-501", "TG", "ENST1", "146-147"),
            gained(145),
            nothing(),
            UTR_COUNTS,
        ),
        (
            ("21_700_AC/TGCA", "21:700-701", "TGCA", "ENST1", "146-147"),
            gained(145),
            nothing(),
            UTR_COUNTS,
        ),
        (kept, gained(145), nothing(), ""),
    ];
    for (row, reference, candidate, counts) in &cases {
        let outcome = remade(
            &one_row(counted_row(*row, reference)),
            &one_row(counted_row_with(*row, candidate, counts)),
        );
        assert!(
            outcome.is("red", &[]),
            "{row:?} {reference:?} {candidate:?} {counts:?}"
        );
    }
    // Under an identifier the input's line says whether the reference allele is empty; without
    // it a two-base Location is read by its positions and its allele alone.
    let inputs = lines(&[
        "21\t1000\tkept_ins\tG\tGTGC\t.\t.\t.",
        "21\t500\tmnv\tAC\tTG\t.\t.\t.",
        "21\t700\tlonger\tAC\tTGCA\t.\t.\t.",
    ]);
    let rows: [(Row, bool, bool); 3] = [
        (kept, true, true),
        (("mnv", "21:500-501", "TG", "ENST1", "146-147"), false, true),
        (
            ("longer", "21:700-701", "TGCA", "ENST1", "146-147"),
            false,
            true,
        ),
    ];
    let texts = |own: &Values| -> Vec<String> {
        rows.iter()
            .map(|(row, _, _)| counted_row(*row, own))
            .collect()
    };
    let reference = default_entries(&default_text(&texts(&gained(145))));
    let candidate = default_entries(&default_text(&texts(&nothing())));
    for (at, (row, handed, alone)) in rows.iter().enumerate() {
        let found = accepts_with_and_without_the_record(
            UTR_REMADE,
            &reference,
            &inputs,
            at,
            &candidate[at],
            UTR_GAIN,
            "-",
        );
        assert_eq!(found, [*handed, *handed, *alone, *alone], "{row:?}");
    }
    let pair = (
        one_row(counted_row(kept, &gained(145))),
        one_row(counted_row(kept, &nothing())),
    );
    let second = inputs[0].replace("GTGC", "GT");
    let given: [Option<Vec<String>>; 4] = [
        Some(vec![inputs[0].clone()]),
        Some(vec![inputs[0].clone(), second]),
        Some(inputs[1..].to_vec()),
        None,
    ];
    for lines in &given {
        let outcome = compared_with(&pair.0, &pair.1, &[UTR_REMADE], lines.as_deref());
        assert!(outcome.is("green", &[(UTR_REMADE, 1)]), "{lines:?}");
    }
    let (e, a) = (&reference[0], &candidate[0]);
    for name in OTHER_FIELDS {
        assert!(
            !accepts(UTR_REMADE, e, name, "1", "2", Some(a), None, None),
            "{name}"
        );
    }
    assert!(!accepts(
        UTR_REMADE,
        e,
        "5UTR_consequence",
        UTR_GAIN,
        "-",
        None,
        None,
        None
    ));
    let moved = counted_row_with(
        kept,
        &nothing(),
        &UTR_COUNTS.replace("Existing_uORFs=1", "Existing_uORFs=2"),
    ) + ";LOEUF=0.4";
    let outcome = remade(
        &one_row(counted_row(kept, &gained(145)) + ";LOEUF=0.3"),
        &one_row(moved),
    );
    assert!(outcome.is("red", &[(UTR_REMADE, 1)]));
    assert_eq!(outcome.differing("Existing_uORFs"), 1);
    assert_eq!(outcome.differing("LOEUF"), 1);
    assert!(outcome.masks("5UTR_consequence", &[(UTR_REMADE, 1)]));
}

/// Where Ensembl VEP's hash order keeps the frameshift on the row, its gained start at the
/// codon's own index is one its runs under other hash seeds write: the candidate's
/// frameshift, which lists the codon's own ORF beside the reference's, is the class's with
/// those runs and a difference without them, as it is where the runs write the gained start
/// at another index, at two, or not at all, and where the reference writes no term. A row
/// two classes hold is counted under the first in the order a row is put to the classes:
/// beside a gained start the reference writes, a frameshift the seeded runs write on the row
/// is the precedence class's, with its annotation where the runs write that too and this
/// class's where they do not; a start the candidate gains inside its allele beside the
/// reference's frameshift is the window class's where the record's alleles are read; the
/// reference's gained start in another order of its pairs is the pair-order class's.
#[test]
fn the_gained_start_remade_in_place_reads_the_seeded_runs_and_follows_the_plugins_other_classes() {
    let kept: Row = ("kept_ins", "21:1000-1001", "TG", "ENST1", "146-147");
    let own = "KozakStrength=Strong:alt_type=uORF:ref_type=uORF";
    let reference_shift = one_row(counted_row(kept, &frame_shifted()));
    let candidate_shift = one_row(counted_row(
        kept,
        &values(UTR_FRAMESHIFT, &format!("{},{own}", frame_shifted().1)),
    ));
    let reference_gain = one_row(counted_row(kept, &gained(145)));
    let at_two = one_row(counted_row(
        kept,
        &values(
            UTR_GAIN,
            "CapDistanceToStart=145:type=uORF,CapDistanceToStart=143:type=uORF",
        ),
    ));
    let gained_at = |at: i64| one_row(counted_row(kept, &gained(at)));
    let sweeps: [(Option<Vec<String>>, bool); 6] = [
        (
            Some(vec![reference_shift.clone(), reference_gain.clone()]),
            true,
        ),
        (None, false),
        (
            Some(vec![reference_shift.clone(), reference_shift.clone()]),
            false,
        ),
        (Some(vec![reference_shift.clone(), gained_at(146)]), false),
        (
            Some(vec![
                reference_shift.clone(),
                reference_gain.clone(),
                gained_at(144),
            ]),
            false,
        ),
        (Some(vec![reference_shift.clone(), at_two]), false),
    ];
    for (sweep, taken) in &sweeps {
        let outcome = match sweep {
            Some(sweep) => compared_swept(
                &reference_shift,
                &candidate_shift,
                Some(sweep),
                &[UTR_REMADE, UTR_PRECEDENCE],
            ),
            None => remade(&reference_shift, &candidate_shift),
        };
        assert!(
            outcome.taken_is(*taken, &[(UTR_REMADE, 1)]),
            "{}",
            sweep.is_none()
        );
    }
    let nothing_written = one_row(counted_row(kept, &nothing()));
    let outcome = compared_swept(
        &nothing_written,
        &candidate_shift,
        Some(&[&reference_shift, &reference_gain]),
        &[UTR_REMADE, UTR_PRECEDENCE],
    );
    assert!(outcome.is("red", &[]));
    // Beside the reference's own gained start: the precedence class first.
    let both_terms = [reference_gain.as_str(), reference_shift.as_str()];
    let cases: [(&str, Classes, Applied); 5] = [
        (
            &reference_shift,
            &[UTR_PRECEDENCE, UTR_REMADE],
            &[(UTR_PRECEDENCE, 1)],
        ),
        (
            &reference_shift,
            &[UTR_REMADE, UTR_PRECEDENCE],
            &[(UTR_PRECEDENCE, 1)],
        ),
        (
            &candidate_shift,
            &[UTR_PRECEDENCE, UTR_REMADE],
            &[(UTR_PRECEDENCE, 1), (UTR_REMADE, 1)],
        ),
        (&reference_shift, &[UTR_REMADE], &[(UTR_REMADE, 1)]),
        (
            &nothing_written,
            &[UTR_PRECEDENCE, UTR_REMADE],
            &[(UTR_REMADE, 1)],
        ),
    ];
    for (candidate, classes, applied) in cases {
        let outcome = compared_swept(&reference_gain, candidate, Some(&both_terms), classes);
        assert!(outcome.is("green", applied), "{classes:?} {applied:?}");
    }
    // A start codon the candidate gains inside its allele, beside the reference's frameshift.
    let making: Row = ("kept_ins", "21:1000-1001", "TGATGC", "ENST1", "146-147");
    let inputs = lines(&["21\t1000\tkept_ins\tG\tGTGATGC\t.\t.\t."]);
    let reference = one_row(counted_row(making, &frame_shifted()));
    let candidate = one_row(counted_row(making, &gained(148)));
    let seeded = [
        reference.clone(),
        one_row(counted_row(making, &gained(145))),
    ];
    let cases: [(Inputs, bool, Applied); 3] = [
        (Some(&inputs), true, &[(UTR_WINDOW, 1)]),
        (Some(&inputs), false, &[(UTR_WINDOW, 1)]),
        (None, true, &[(UTR_PRECEDENCE, 1), (UTR_REMADE, 1)]),
    ];
    for (lines, swept, applied) in cases {
        let outcome = compared_entries(
            &default_entries(&reference),
            &default_entries(&candidate),
            &[UTR_WINDOW, UTR_PRECEDENCE, UTR_REMADE],
            lines,
            swept.then(|| seed_sweep(&seeded)),
        );
        assert!(outcome.is("green", applied), "{lines:?} {swept}");
    }
    // Both engines name the made codon's neighbour at the old index and the further one: the
    // annotation alone differs and is this class's.
    let outcome = compared(
        &one_row(counted_row(making, &gained(145))),
        &candidate,
        &[UTR_REMADE, UTR_WINDOW],
    );
    assert!(outcome.is("green", &[(UTR_REMADE, 1)]));
    assert!(outcome.masks("5UTR_annotation", &[(UTR_REMADE, 1)]));
    let hash_order = one_row(counted_row(
        kept,
        &values(
            UTR_GAIN,
            "type=uORF:KozakStrength=Weak:DistanceToCDS=110:CapDistanceToStart=145",
        ),
    ));
    let outcome = compared(&hash_order, &reference_gain, &[PAIR_ORDER, UTR_REMADE]);
    assert!(outcome.is("green", &[(PAIR_ORDER, 1)]));
}

/// The command line names the remade-in-place class like any other.
#[test]
fn the_command_line_names_the_gained_start_remade_in_place_class() {
    let kept: Row = ("kept_ins", "21:1000-1001", "TG", "ENST1", "146-147");
    let reference_gain = one_row(counted_row(kept, &gained(145)));
    let nothing_written = one_row(counted_row(kept, &nothing()));
    let directory = planted_files(
        "utrannotator-remade",
        &[
            ("reference", &reference_gain),
            ("candidate", &nothing_written),
        ],
    );
    assert_eq!(exit_code(&directory, &[]), 1);
    assert_eq!(exit_code(&directory, &["--masks", UTR_REMADE]), 0);
}

/// The tab format names the strand as a column, the VCF format as a CSQ field and the JSON
/// format as the consequence object's `strand`: on a reverse-strand transcript a gained
/// start one below the lower cDNA position is the class's where the row's allele,
/// reverse-complemented, opens with the codon's last two bases, and stays a difference
/// where it does not or the format carries no strand.
#[test]
fn the_gained_start_remade_in_place_reads_the_tab_vcf_and_json_formats() {
    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tSTRAND\t5UTR_consequence\t5UTR_annotation\tExisting_uORFs\n";
    let tab = |allele: &str, strand: &str, own: &Values| {
        format!(
            "{header}kept_ins\t21:1000-1001\t{allele}\tENSG1\tENST1\tTranscript\t5_prime_UTR_variant\t146-147\t{strand}\t{}\t{}\t1\n",
            own.0, own.1
        )
    };
    for (allele, strand, taken) in [
        ("GCA", "-1", true),
        ("TGC", "1", true),
        ("TGC", "-1", false),
        ("GCA", "1", false),
        ("GCA", "-", false),
    ] {
        let outcome = compared_entries(
            &tab_entries(&tab(allele, strand, &gained(145))),
            &tab_entries(&tab(allele, strand, &nothing())),
            &[UTR_REMADE],
            None,
            None,
        );
        assert!(
            outcome.taken_is(taken, &[(UTR_REMADE, 1)]),
            "{allele} {strand}"
        );
    }

    let inputs = lines(&["21\t1000\tkept_ins\tC\tCGCA\t.\t.\t."]);
    let gain = format!("{UTR_GAIN}|CapDistanceToStart=145:type=uORF");
    for (fields, strand, taken) in [
        ("STRAND|", "-1|", true),
        ("STRAND|", "1|", false),
        ("", "", false),
    ] {
        let csq = csq_header(&format!(
            "Allele|Consequence|Feature|Feature_type|cDNA_position|{fields}5UTR_consequence|5UTR_annotation|Existing_uORFs"
        ));
        let vcf_reference = format!(
            "{csq}{}CSQ=GCA|5_prime_UTR_variant|ENST1|Transcript|146-147|{strand}{gain}|1\n",
            record_columns(&inputs[0])
        );
        let outcome = compared_entries(
            &vcf_entries(&vcf_reference, &inputs),
            &vcf_entries(&vcf_reference.replace(&gain, "|"), &inputs),
            &[UTR_REMADE],
            None,
            None,
        );
        assert!(
            outcome.taken_is(taken, &[(UTR_REMADE, 1)]),
            "{fields} {strand}"
        );
    }

    let object = |strand: Option<i64>, gains: bool| {
        let mut entry = r#""variant_allele": "GCA", "transcript_id": "ENST1", "consequence_terms": ["5_prime_UTR_variant"], "cdna_start": 147, "cdna_end": 146, "existing_uorfs": 1"#.to_string();
        if let Some(strand) = strand {
            entry.push_str(&format!(r#", "strand": {strand}"#));
        }
        if gains {
            entry.push_str(&format!(
                r#", "5utr_consequence": "{UTR_GAIN}", "5utr_annotation": {{"1": {{"CapDistanceToStart": 145, "type": "uORF"}}}}"#
            ));
        }
        json_line(&inputs[0], "r0", "-/GCA", &entry)
    };
    for (strand, taken) in [(Some(-1), true), (Some(1), false), (None, false)] {
        let outcome = compared_entries(
            &json_entries(&object(strand, true), &inputs),
            &json_entries(&object(strand, false), &inputs),
            &[UTR_REMADE],
            None,
            None,
        );
        assert!(outcome.taken_is(taken, &[(UTR_REMADE, 1)]), "{strand:?}");
    }
}

/// On an insertion's row, where the gained-start window reads every 5UTR_ field, a
/// 5UTR_annotation differing in order alone is counted under the order class.
#[test]
fn the_annotation_pair_order_is_counted_before_the_window_class() {
    // The seven pairs of a substitution in the 5' UTR of one transcript that destroys an
    // upstream start codon, in one order Ensembl VEP writes them and in the order of their keys.
    let hash_order = "DistanceToStop=24:CapDistanceToStart=407:KozakContext=ATAATGG:KozakStrength=Strong:DistanceToCDS=23:Evidence=False:type=OutOfFrame_oORF";
    let key_order = "CapDistanceToStart=407:DistanceToCDS=23:DistanceToStop=24:Evidence=False:KozakContext=ATAATGG:KozakStrength=Strong:type=OutOfFrame_oORF";
    let insertion = one_row(utr_row(
        ("ins", "21:25734954-25734955", "T", "ENST1", "12-13"),
        &values(UTR_LOSS, hash_order),
    ));
    let outcome = compared(
        &insertion,
        &insertion.replace(hash_order, key_order),
        &[UTR_WINDOW, PAIR_ORDER],
    );
    assert!(outcome.is("green", &[(PAIR_ORDER, 1)]));
    let e = &default_entries(&one_row(default_row(
        "rs1",
        "21:100",
        "A",
        "ENST1",
        "intron_variant",
        "IMPACT=MODIFIER",
    )))[0];
    for (value, accepted) in [("0", true), ("0.5", false)] {
        assert_eq!(
            accepts(
                "loeuf_zero_minimum",
                e,
                "LOEUF",
                ABSENT,
                value,
                None,
                None,
                None
            ),
            accepted
        );
    }
}

/// A record of substitutions whose alleles share their first base (AC to AT,AG) is written
/// as it stands, so its AT row is a substitution's; beside a deletion (AC to AT,A) the base
/// is dropped and the rows read T and -. A 5UTR_ field differing on the substitution's row
/// is therefore not the gained-start class's.
#[test]
fn allele_is_indel_keeps_the_shared_first_base_of_a_record_without_an_indel() {
    assert_eq!(
        allele_is_indel("AC", &["AT", "AG"], "AT", false),
        Some(false)
    );
    assert_eq!(
        allele_is_indel("AC", &["AT", "AG"], "AG", false),
        Some(false)
    );
    assert_eq!(allele_is_indel("AC", &["AT", "A"], "T", false), Some(false));
    assert_eq!(allele_is_indel("AC", &["AT", "A"], "-", false), Some(true));
    assert_eq!(
        allele_is_indel("GC", &["GT", "GA", "*"], "T", false),
        Some(false)
    );
    let reference = one_row(default_row(
        "21_100_AC/AT/AG",
        "21:100-101",
        "AT",
        "ENST1",
        "5_prime_UTR_variant",
        "IMPACT=MODIFIER;STRAND=1;5UTR_consequence=uAUG_lost;5UTR_annotation=type=uORF",
    ));
    let outcome = compared(
        &reference,
        &reference.replace("uAUG_lost", "uAUG_gained"),
        &[UTR_WINDOW],
    );
    assert!(outcome.is("red", &[]));
}
