// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvsp_start_codon_test_past_the_shifted_allele` on the HGVSp of an insertion or deletion
//! that the 3' shift moves onto the start codon and that leaves the whole coding sequence
//! standing: the reference writes a frameshift at residue 1 with the unknown count, or an
//! insertion between a residue 0 and residue 1, and the candidate the first residue
//! unchanged.
//! The rows are made-up rows of three alleles in the columns a run under
//! `--hgvs --shift_3prime 1 --shift_length` writes: an insertion of one base, an insertion of
//! three bases and a deletion of one base, each on a transcript of its own.

mod compare_common;
mod shifted_edited_common;

use compare_common::*;
use shifted_edited_common::*;
use std::collections::BTreeMap;
use vep_formats_compare::classes;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::downstream::field_of;

const CLASS: &str = "hgvsp_start_codon_test_past_the_shifted_allele";
const WRITTEN_SPAN: &str = "hgvsp_shifted_indel_written_span_start_test";

/// The two rows of one allele: the thirteen columns, the HGVSc and the SHIFT_LENGTH both
/// engines write alike, and each engine's HGVSp, the reference's first.
#[derive(Clone, Copy)]
struct Written<'a> {
    columns: &'a str,
    hgvsc: &'a str,
    hgvsp: [&'a str; 2],
    shift: &'a str,
}

/// The insertion of one base, moved two bases 3': the reference's notation is a frameshift.
const INSERTION: Written = Written {
    columns: "inserted_onto_the_start_codon\t21:1000-1001\tA\tENSG1\tENST1\tTranscript\tframeshift_variant,start_retained_variant\t102-103\t2-3\t1\t-\t-\t-",
    hgvsc: "ENST1.1:c.2dup",
    hgvsp: ["ENSP1.1:p.Met1TrpfsTer?", "ENSP1.1:p.Met1%3D"],
    shift: "2",
};

/// The insertion of three bases, moved three bases 3': the reference's notation is an
/// insertion between a residue 0 and residue 1.
const IN_FRAME_INSERTION: Written = Written {
    columns: "inserted_in_frame_onto_the_start_codon\t21:2000-2001\tTCA\tENSG2\tENST2\tTranscript\tinframe_insertion,start_retained_variant\t52-53\t1-2\t1\t-\t-\t-",
    hgvsc: "ENST2.3:c.-2_1dup",
    hgvsp: ["ENSP2.3:p.Ala0_Ala1insSer", "ENSP2.3:p.Met1%3D"],
    shift: "3",
};

/// The deletion of one base, moved one base 3': the reference's notation is a frameshift.
const DELETION: Written = Written {
    columns: "deleted_onto_the_start_codon\t21:3000\t-\tENSG3\tENST3\tTranscript\tframeshift_variant,start_retained_variant\t12\t2\t1\t-\t-\t-",
    hgvsc: "ENST3.2:c.2del",
    hgvsp: ["ENSP3.2:p.Met1ValfsTer?", "ENSP3.2:p.Met1%3D"],
    shift: "1",
};

const ALLELES: [Written; 3] = [INSERTION, IN_FRAME_INSERTION, DELETION];

/// One default-format row: `columns`, then an Extra of the three fields, each left out
/// where empty.
fn row(columns: &str, hgvsc: &str, hgvsp: &str, shift: &str) -> String {
    let mut extra = String::from("IMPACT=HIGH;STRAND=1");
    for (name, value) in [("SHIFT_LENGTH", shift), ("HGVSc", hgvsc), ("HGVSp", hgvsp)] {
        if !value.is_empty() {
            extra.push_str(&format!(";{name}={value}"));
        }
    }
    format!("{columns}\t{extra}")
}

impl Written<'_> {
    fn reference(&self) -> String {
        row(self.columns, self.hgvsc, self.hgvsp[0], self.shift)
    }

    fn candidate(&self) -> String {
        row(self.columns, self.hgvsc, self.hgvsp[1], self.shift)
    }

    /// The row of the engine at `side` with the description `description` after its
    /// protein's name.
    fn describing(&self, side: usize, description: &str) -> String {
        let name = self.hgvsp[side].split(":p.").next().unwrap_or_default();
        row(
            self.columns,
            self.hgvsc,
            &format!("{name}:p.{description}"),
            self.shift,
        )
    }

    /// The reference's row with the description `description`.
    fn reference_describing(&self, description: &str) -> String {
        self.describing(0, description)
    }

    /// The candidate's row with the description `description`.
    fn candidate_describing(&self, description: &str) -> String {
        self.describing(1, description)
    }

    /// The column at `place` of the thirteen.
    fn column(&self, place: usize) -> &str {
        self.columns.split('\t').nth(place).unwrap_or_default()
    }
}

/// Whether the class takes the HGVSp of a pair of default-format rows.
fn takes(reference: &str, candidate: &str) -> bool {
    let (e, a) = (default_entry(reference), default_entry(candidate));
    class_accepts(
        CLASS,
        &e,
        "HGVSp",
        &field_of(Some(&e), &["HGVSp"]),
        &field_of(Some(&a), &["HGVSp"]),
        Some(&a),
    )
}

/// The verdict of two entry lists under `classes` and what was applied.
fn under(
    classes: &[&str],
    expected: &[Entry],
    actual: &[Entry],
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let report = compare(expected, actual, &mut documented);
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// The verdict of two lists of default-format rows under `classes` and what was applied.
fn rows_under(
    classes: &[&str],
    reference: &[String],
    candidate: &[String],
) -> (String, BTreeMap<String, i64>) {
    let (report, documented) =
        compare_default(&default_text(reference), &default_text(candidate), classes);
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// The rows of the three alleles are the class, each row counted once, and no row of the
/// class of a start codon test read on the written span. Without the class the three HGVSp
/// values are differences.
#[test]
fn the_rows_of_the_three_alleles_are_the_class() {
    let reference: Vec<String> = ALLELES.iter().map(Written::reference).collect();
    let candidate: Vec<String> = ALLELES.iter().map(Written::candidate).collect();
    let thrice = ("green".to_string(), count_map(&[(CLASS, 3)]));
    assert_eq!(rows_under(&[CLASS], &reference, &candidate), thrice);
    assert_eq!(
        rows_under(&[WRITTEN_SPAN, CLASS], &reference, &candidate),
        thrice
    );
    assert_eq!(
        rows_under(&[WRITTEN_SPAN], &reference, &candidate),
        ("red".to_string(), nothing())
    );
    for allele in ALLELES {
        assert!(
            takes(&allele.reference(), &allele.candidate()),
            "{}",
            allele.columns
        );
        assert_eq!(
            rows_under(&[], &[allele.reference()], &[allele.candidate()]).0,
            "red"
        );
    }
    // A row of a transcript whose protein nonsense-mediated decay takes carries its term
    // beside the retained start's.
    let decayed = |text: String| {
        text.replace(
            "start_retained_variant",
            "start_retained_variant,NMD_transcript_variant",
        )
    };
    assert!(takes(
        &decayed(INSERTION.reference()),
        &decayed(INSERTION.candidate())
    ));
}

/// The tab, VCF and JSON formats carry the cells apart and count a row once; the JSON
/// format writes the unchanged residue's `=` as it is and the CDS position as a start and
/// an end.
#[test]
fn the_class_reads_the_tab_vcf_and_json_formats() {
    let once = ("green".to_string(), count_map(&[(CLASS, 1)]));
    let header = tab_header(&["IMPACT", "STRAND", "SHIFT_LENGTH", "HGVSc", "HGVSp"]);
    for allele in ALLELES {
        let tab = |hgvsp: &str| {
            tab_entries(&format!(
                "{header}\n{}\tHIGH\t1\t{}\t{}\t{hgvsp}\n",
                allele.columns, allele.shift, allele.hgvsc
            ))
        };
        assert_eq!(
            under(&[CLASS], &tab(allele.hgvsp[0]), &tab(allele.hgvsp[1])),
            once,
            "{}",
            allele.columns
        );
        assert_eq!(
            under(&[], &tab(allele.hgvsp[0]), &tab(allele.hgvsp[1])).0,
            "red"
        );
    }

    let csq = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|SHIFT_LENGTH|HGVSc|HGVSp",
    );
    let records = [
        "21\t1000\tinserted_onto_the_start_codon\tC\tCA\t.\t.\t.",
        "21\t2000\tinserted_in_frame_onto_the_start_codon\tG\tGTCA\t.\t.\t.",
        "21\t2999\tdeleted_onto_the_start_codon\tGA\tG\t.\t.\t.",
    ];
    for (allele, record) in ALLELES.iter().zip(records) {
        let input = vec![record.to_string()];
        let vcf = |hgvsp: &str| {
            vcf_entries(
                &format!(
                    "{csq}{}CSQ={}|{}|{}|Transcript|{}|{}|{}|{}|{}|{hgvsp}\n",
                    record_columns(&input[0]),
                    allele.column(2),
                    allele.column(6).replace(',', "&"),
                    allele.column(4),
                    allele.column(7),
                    allele.column(8),
                    allele.column(9),
                    allele.shift,
                    allele.hgvsc
                ),
                &input,
            )
        };
        assert_eq!(
            under(&[CLASS], &vcf(allele.hgvsp[0]), &vcf(allele.hgvsp[1])),
            once,
            "{record}"
        );
        assert_eq!(
            under(&[], &vcf(allele.hgvsp[0]), &vcf(allele.hgvsp[1])).0,
            "red"
        );
    }

    // The insertion of one base as a JSON object, its reference notation a frameshift.
    let inserted = vec![records[0].to_string()];
    let record = [
        ("input", text(&inserted[0])),
        ("id", text("inserted_onto_the_start_codon")),
        ("allele_string", text("-/A")),
        ("most_severe_consequence", text("frameshift_variant")),
    ];
    let consequence = object(&[
        ("variant_allele", text("A")),
        ("transcript_id", text("ENST1")),
        (
            "consequence_terms",
            Value::List(vec![
                text("frameshift_variant"),
                text("start_retained_variant"),
            ]),
        ),
        ("cdna_start", Value::Int(102)),
        ("cdna_end", Value::Int(103)),
        ("cds_start", Value::Int(2)),
        ("cds_end", Value::Int(3)),
        ("protein_start", Value::Int(1)),
        ("protein_end", Value::Int(1)),
        ("shift_length", Value::Int(2)),
        ("hgvsc", text(INSERTION.hgvsc)),
    ]);
    let entries_of = |consequence: &Value, hgvsp: &str| {
        json_entries(
            &json_line(
                &record,
                with_members(consequence, &[("hgvsp", text(hgvsp))]),
            ),
            &inserted,
        )
    };
    let unchanged = INSERTION.hgvsp[1].replace("%3D", "=");
    assert_eq!(
        under(
            &[CLASS],
            &entries_of(&consequence, INSERTION.hgvsp[0]),
            &entries_of(&consequence, &unchanged)
        ),
        once
    );
    assert_eq!(
        under(
            &[],
            &entries_of(&consequence, INSERTION.hgvsp[0]),
            &entries_of(&consequence, &unchanged)
        )
        .0,
        "red"
    );
    // An object that leaves out the CDS start, the CDS end or the shift length, or whose
    // shift length is zero, is no row of the class.
    for other in [
        without_member(&consequence, "cds_start"),
        without_member(&consequence, "cds_end"),
        without_member(&consequence, "shift_length"),
        with_members(&consequence, &[("shift_length", Value::Int(0))]),
    ] {
        assert_eq!(
            under(
                &[CLASS],
                &entries_of(&other, INSERTION.hgvsp[0]),
                &entries_of(&other, &unchanged)
            ),
            ("red".to_string(), nothing())
        );
    }

    // The insertion of three bases as a JSON object, its reference notation an insertion
    // between a residue 0 and residue 1.
    let inserted = vec![records[1].to_string()];
    let record = [
        ("input", text(&inserted[0])),
        ("id", text("inserted_in_frame_onto_the_start_codon")),
        ("allele_string", text("-/TCA")),
        ("most_severe_consequence", text("inframe_insertion")),
    ];
    let consequence = object(&[
        ("variant_allele", text("TCA")),
        ("transcript_id", text("ENST2")),
        (
            "consequence_terms",
            Value::List(vec![
                text("inframe_insertion"),
                text("start_retained_variant"),
            ]),
        ),
        ("cdna_start", Value::Int(52)),
        ("cdna_end", Value::Int(53)),
        ("cds_start", Value::Int(1)),
        ("cds_end", Value::Int(2)),
        ("protein_start", Value::Int(1)),
        ("protein_end", Value::Int(1)),
        ("shift_length", Value::Int(3)),
        ("hgvsc", text(IN_FRAME_INSERTION.hgvsc)),
    ]);
    let entries_of = |hgvsp: &str| {
        json_entries(
            &json_line(
                &record,
                with_members(&consequence, &[("hgvsp", text(hgvsp))]),
            ),
            &inserted,
        )
    };
    let unchanged = IN_FRAME_INSERTION.hgvsp[1].replace("%3D", "=");
    assert_eq!(
        under(
            &[CLASS],
            &entries_of(IN_FRAME_INSERTION.hgvsp[0]),
            &entries_of(&unchanged)
        ),
        once
    );
    assert_eq!(
        under(
            &[CLASS],
            &entries_of("ENSP2.3:p.(Ala0_Ala1insSerVal)"),
            &entries_of("ENSP2.3:p.(Met1=)")
        ),
        once
    );
    assert_eq!(
        under(
            &[CLASS],
            &entries_of("ENSP2.3:p.Ala0_Val1insSer"),
            &entries_of(&unchanged)
        ),
        ("red".to_string(), nothing())
    );
}

/// An empty value on either side, a count in the reference's frameshift, another residue
/// number, unequal flank names and a candidate of another form stay differences, and so do
/// the sides exchanged and two values that name different proteins.
#[test]
fn any_other_pair_of_notations_stays_a_difference() {
    for allele in ALLELES {
        let (reference, candidate) = (allele.reference(), allele.candidate());
        let without = |text: &str, hgvsp: &str| text.replace(&format!(";HGVSp={hgvsp}"), "");
        assert!(!takes(&without(&reference, allele.hgvsp[0]), &candidate));
        assert!(!takes(&reference, &without(&candidate, allele.hgvsp[1])));
        assert!(!takes(&candidate, &reference));
        let renamed = candidate.replace("HGVSp=ENSP", "HGVSp=ENSQ");
        assert!(!takes(&reference, &renamed), "{}", allele.columns);
        assert_eq!(
            rows_under(&[CLASS], &[allele.reference()], &[renamed]).1,
            nothing()
        );
    }
    for (allele, reference, candidate) in [
        // A count in the reference's frameshift.
        (INSERTION, "Met1TrpfsTer9", "Met1%3D"),
        (DELETION, "Met1ValfsTer12", "Met1%3D"),
        // Another residue number, on one side or on both.
        (INSERTION, "Met2TrpfsTer?", "Met1%3D"),
        (INSERTION, "Met1TrpfsTer?", "Met2%3D"),
        (INSERTION, "Met2TrpfsTer?", "Met2%3D"),
        (INSERTION, "Met11TrpfsTer?", "Met11%3D"),
        (IN_FRAME_INSERTION, "Ala1_Ala2insSer", "Met1%3D"),
        (IN_FRAME_INSERTION, "Ala0_Ala2insSer", "Met1%3D"),
        (IN_FRAME_INSERTION, "Ala00_Ala1insSer", "Met1%3D"),
        (IN_FRAME_INSERTION, "Ala0_Ala1insSer", "Met2%3D"),
        (IN_FRAME_INSERTION, "Ala0_Ala1insSer", "Met0%3D"),
        // Unequal flank names in the reference's insertion.
        (IN_FRAME_INSERTION, "Ala0_Val1insSer", "Met1%3D"),
        (IN_FRAME_INSERTION, "Ala0_Met1insSer", "Met1%3D"),
        // A frameshift whose first residue is not the candidate's.
        (INSERTION, "Val1TrpfsTer?", "Met1%3D"),
        (INSERTION, "Met1TrpfsTer?", "Val1%3D"),
        // A candidate of another form.
        (INSERTION, "Met1TrpfsTer?", "Met1?"),
        (INSERTION, "Met1TrpfsTer?", "Met1Trp"),
        (INSERTION, "Met1TrpfsTer?", "Met1TrpfsTer3"),
        (INSERTION, "Met1TrpfsTer?", "Met1_?3"),
        (INSERTION, "Met1TrpfsTer?", "Met1%3D%3D"),
        (INSERTION, "Met1TrpfsTer?", "Met1del"),
        (IN_FRAME_INSERTION, "Ala0_Ala1insSer", "Met1?"),
        (IN_FRAME_INSERTION, "Ala0_Ala1insSer", "Met1_Ala2insSer"),
        (IN_FRAME_INSERTION, "Ala0_Ala1insSer", "Met1dup"),
        // A reference of another form: an insertion of no residue or of a cut name, a
        // duplication, a replacement, a frameshift that names no second residue, the lost
        // start and a deletion.
        (IN_FRAME_INSERTION, "Ala0_Ala1ins", "Met1%3D"),
        (IN_FRAME_INSERTION, "Ala0_Ala1insSerV", "Met1%3D"),
        (IN_FRAME_INSERTION, "Ala0_Ala1dup", "Met1%3D"),
        (IN_FRAME_INSERTION, "Ala0_Ala1delinsSer", "Met1%3D"),
        (INSERTION, "Met1fsTer?", "Met1%3D"),
        (INSERTION, "Met1?", "Met1%3D"),
        (DELETION, "Met1del", "Met1%3D"),
        (DELETION, "Met1_Pro3del", "Met1%3D"),
    ] {
        assert!(
            !takes(
                &allele.reference_describing(reference),
                &allele.candidate_describing(candidate)
            ),
            "{reference} {candidate}"
        );
    }
    // Each pair above is one cell away from a pair the class takes.
    for (allele, reference, candidate) in [
        (INSERTION, "Met1TrpfsTer?", "Met1%3D"),
        (INSERTION, "Met1TrpfsTer?", "Met1="),
        (INSERTION, "Leu1TrpfsTer?", "Leu1%3D"),
        (DELETION, "Met1ValfsTer?", "Met1%3D"),
        (IN_FRAME_INSERTION, "Ala0_Ala1insSer", "Met1%3D"),
        (IN_FRAME_INSERTION, "Ala0_Ala1insSerVal", "Met1="),
        (IN_FRAME_INSERTION, "Met0_Met1insSer", "Met1%3D"),
        (IN_FRAME_INSERTION, "Ala0_Ala1insSer", "Leu1%3D"),
    ] {
        assert!(
            takes(
                &allele.reference_describing(reference),
                &allele.candidate_describing(candidate)
            ),
            "{reference} {candidate}"
        );
    }
}

/// A row whose terms hold `start_lost`, lack `start_retained_variant` or differ between the
/// engines, a row without one positive SHIFT_LENGTH on both sides, a row that is no
/// insertion's or deletion's and a row whose CDS position lacks its start or its end stay
/// differences; the class reads no other field and nothing without the candidate's row.
#[test]
fn any_other_row_stays_a_difference() {
    for allele in ALLELES {
        let (reference, candidate) = (allele.reference(), allele.candidate());
        // The lost start beside the retained one, and in its place, on both rows.
        for terms in ["start_lost,start_retained_variant", "start_lost"] {
            let tested = |text: &str| text.replace("start_retained_variant", terms);
            assert!(!takes(&tested(&reference), &tested(&candidate)), "{terms}");
            assert_eq!(
                rows_under(&[CLASS], &[tested(&reference)], &[tested(&candidate)]),
                ("red".to_string(), nothing())
            );
        }
        // No retained start, and terms that differ.
        let elsewhere = |text: &str| text.replace(",start_retained_variant", "");
        assert!(!takes(&elsewhere(&reference), &elsewhere(&candidate)));
        assert!(!takes(&reference, &elsewhere(&candidate)));
        assert!(!takes(&elsewhere(&reference), &candidate));
        // The shift length: none, a zero, two that differ, one on one side alone, and a
        // text that is no positive integer.
        let shift = format!(";SHIFT_LENGTH={}", allele.shift);
        let shifted = |text: &str, to: &str| text.replace(&shift, to);
        for (on_reference, on_candidate) in [
            ("", ""),
            (";SHIFT_LENGTH=0", ";SHIFT_LENGTH=0"),
            (shift.as_str(), ";SHIFT_LENGTH=7"),
            (shift.as_str(), ""),
            ("", shift.as_str()),
            (";SHIFT_LENGTH=-2", ";SHIFT_LENGTH=-2"),
            (";SHIFT_LENGTH=02", ";SHIFT_LENGTH=02"),
            (";SHIFT_LENGTH=2.0", ";SHIFT_LENGTH=2.0"),
        ] {
            assert!(
                !takes(
                    &shifted(&reference, on_reference),
                    &shifted(&candidate, on_candidate)
                ),
                "{on_reference} {on_candidate}"
            );
        }
        assert!(takes(
            &shifted(&reference, ";SHIFT_LENGTH=40"),
            &shifted(&candidate, ";SHIFT_LENGTH=40")
        ));
        // Any other field of the row, and a comparison that holds no candidate row.
        let (e, a) = (default_entry(&reference), default_entry(&candidate));
        for name in ["HGVSc", "SHIFT_LENGTH", "Consequence", "Extra"] {
            assert!(
                !class_accepts(CLASS, &e, name, allele.hgvsp[0], allele.hgvsp[1], Some(&a)),
                "{name}"
            );
        }
        assert!(class_accepts(
            CLASS,
            &e,
            "hgvsp",
            allele.hgvsp[0],
            allele.hgvsp[1],
            Some(&a)
        ));
        assert!(!class_accepts(
            CLASS,
            &e,
            "HGVSp",
            allele.hgvsp[0],
            allele.hgvsp[1],
            None
        ));
        // A CDS position with no start, with no end or with neither, on both rows or on
        // one alone.
        let positions = format!("\t{}\t{}\t", allele.column(7), allele.column(8));
        let placed = |text: &str, to: &str| {
            text.replacen(&positions, &format!("\t{}\t{to}\t", allele.column(7)), 1)
        };
        for to in ["-", "?-3", "2-?", "?", "?-3/1500"] {
            assert!(
                !takes(&placed(&reference, to), &placed(&candidate, to)),
                "{to}"
            );
        }
        assert!(!takes(&reference, &placed(&candidate, "?-3")));
        assert!(!takes(&placed(&reference, "?-3"), &candidate));
        // The two ends may be any bases, one number stands for both, and the column may
        // carry a total.
        for to in ["1-2", "7", "2-3/1500"] {
            assert!(
                takes(&placed(&reference, to), &placed(&candidate, to)),
                "{to}"
            );
        }
    }
    // One base replaced by another is no insertion and no deletion.
    let replaced = |text: String| text.replace("21:1000-1001", "21:1001");
    assert!(!takes(
        &replaced(INSERTION.reference()),
        &replaced(INSERTION.candidate())
    ));
}

/// A pair of prediction parentheses that both values carry is set aside; one value in
/// parentheses beside one without stays a difference.
#[test]
fn a_pair_of_prediction_parentheses_on_both_values_is_set_aside() {
    for (allele, reference, candidate) in [
        (INSERTION, "(Met1TrpfsTer?)", "(Met1%3D)"),
        (INSERTION, "(Met1TrpfsTer?)", "(Met1=)"),
        (IN_FRAME_INSERTION, "(Ala0_Ala1insSer)", "(Met1%3D)"),
        (DELETION, "(Met1ValfsTer?)", "(Met1%3D)"),
    ] {
        let predicted = allele.reference_describing(reference);
        assert!(takes(&predicted, &allele.candidate_describing(candidate)));
        assert!(!takes(&predicted, &allele.candidate()));
        assert!(!takes(
            &allele.reference(),
            &allele.candidate_describing(candidate)
        ));
    }
    assert!(!takes(
        &INSERTION.reference_describing("(Met1TrpfsTer?"),
        &INSERTION.candidate_describing("(Met1%3D")
    ));
}

/// Two inserted bases at a two-base Location read as a substitution from the row's text
/// under an identifier; the candidate's duplication or insertion notation marks the row an
/// insertion's, and a substitution's notation does not.
#[test]
fn a_two_base_insertion_under_an_identifier_is_read_from_the_candidates_notation() {
    let columns = INSERTION.columns.replace("\tA\t", "\tAT\t");
    let pair = |hgvsc: &str| {
        let written = |hgvsp: &str| row(&columns, hgvsc, hgvsp, "2");
        takes(&written(INSERTION.hgvsp[0]), &written(INSERTION.hgvsp[1]))
    };
    assert!(pair("ENST1.1:c.1_2dup"));
    assert!(pair("ENST1.1:c.2_3insAT"));
    assert!(!pair("ENST1.1:c.1_2delinsTA"));
    assert!(!pair("ENST9.1:c.1_2dup"));
    assert!(!pair(""));
    // A row whose first column spells the record's alleles says so itself.
    let spelt = columns.replace("inserted_onto_the_start_codon", "21_1001_-/AT");
    let written = |hgvsp: &str| row(&spelt, "", hgvsp, "2");
    assert!(takes(
        &written(INSERTION.hgvsp[0]),
        &written(INSERTION.hgvsp[1])
    ));
}

/// The class has its constant, reads the HGVSp under either letter case and stands in the
/// class table after the class of a start codon kept without a 5' UTR, before the class of
/// the LoFTEE flags and the two classes that close the table; it reads the two rows alone,
/// and the `--masks` help names it.
#[test]
fn the_class_tables_and_the_help_carry_the_class() {
    assert_eq!(
        classes::HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE,
        CLASS
    );
    assert!(classes::NAMES.iter().any(|(name, value)| {
        *name == "HGVSP_START_CODON_TEST_PAST_THE_SHIFTED_ALLELE" && *value == CLASS
    }));
    let tables = classes::tables();
    assert_eq!(
        tables.field_classes_of(CLASS),
        Some(&["HGVSp", "hgvsp"][..])
    );
    let names: Vec<&str> = tables
        .field_classes()
        .iter()
        .map(|(class, _)| *class)
        .collect();
    assert_eq!(
        names[names.len() - 5..],
        [
            classes::START_CODON_KEPT_WITHOUT_FIVE_PRIME_UTR,
            CLASS,
            classes::LOFTEE_NAGNAG_SITE_FIRST_TRANSCRIPT_STRAND,
            classes::PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN,
            classes::PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR,
        ]
    );
    for table in [
        tables.in_term_pair_following_classes(CLASS),
        tables.in_record_reading_classes(CLASS),
        tables.in_record_reference_classes(CLASS),
        tables.in_allele_reading_classes(CLASS),
        tables.in_transcript_reading_classes(CLASS),
        tables.in_sibling_allele_classes(CLASS),
    ] {
        assert!(!table);
    }
    let help = include_str!("../data/help.txt");
    assert_eq!(help.matches(CLASS).count(), 1);
}
