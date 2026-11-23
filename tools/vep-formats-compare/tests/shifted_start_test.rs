// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `hgvsp_shifted_indel_written_span_start_test` on the HGVSp of an insertion or deletion
//! that the 3' shift moves into the first codon: the reference writes a frameshift at
//! residue 1 with the unknown count, or a deletion of residues from the first, and the
//! candidate the notation of the allele written at its 3' place.
//! The rows are those both engines write under `--hgvs` for made-up records on the
//! committed corpus `tests/golden/116/GRCh38-custom`: an insertion written at the last base
//! before the coding start of `ENST00000400183`, and two deletions that begin before the
//! first base of `ENST00000652508`, a transcript with no 5' UTR.

mod compare_common;
mod shifted_edited_common;

use compare_common::*;
use shifted_edited_common::*;
use std::collections::BTreeMap;
use vep_formats_compare::classes;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::downstream::field_of;

const CLASS: &str = "hgvsp_shifted_indel_written_span_start_test";
const CLAMPED_SPAN: &str = "hgvs_shifted_deletion_clamped_span";

/// The two rows of one allele as its record writes it: the thirteen columns both engines
/// write alike, then each engine's HGVSc and HGVSp, the reference's first, and the
/// HGVS_OFFSET both write.
#[derive(Clone, Copy)]
struct Written<'a> {
    columns: &'a str,
    hgvsc: [&'a str; 2],
    hgvsp: [&'a str; 2],
    offset: &'a str,
}

/// The insertion of one base after the last base of the 5' UTR, moved one base 3'.
const INSERTION: Written = Written {
    columns: "inserted_before_the_coding_start\t21:15730393-15730394\tA\tENSG00000155313\tENST00000400183\tTranscript\t5_prime_UTR_variant\t412-413\t-\t-\t-\t-\t-",
    hgvsc: ["ENST00000400183.7:c.1dup", "ENST00000400183.7:c.1dup"],
    hgvsp: [
        "ENSP00000383044.2:p.Met1AsnfsTer?",
        "ENSP00000383044.2:p.Met1%3D",
    ],
    offset: "1",
};

/// The deletion of seven bases, one before the transcript's first base, moved one base 3'.
const FRAMESHIFT_DELETION: Written = Written {
    columns: "deleted_across_the_first_base\t21:46326287-46326293\t-\tENSG00000286224\tENST00000652508\tTranscript\tcoding_sequence_variant,5_prime_UTR_variant\t?-6\t?-6\t?-2\t-\t-\t-",
    hgvsc: ["ENST00000652508.1:c.2_7del", "ENST00000652508.1:c.1_7del"],
    hgvsp: [
        "ENSP00000498568.1:p.Met1ProfsTer?",
        "ENSP00000498568.1:p.Met1_?3",
    ],
    offset: "1",
};

/// The deletion of nine bases, one before the transcript's first base, moved two bases 3'.
const IN_FRAME_DELETION: Written = Written {
    columns: "deleted_in_frame_across_the_first_base\t21:46326287-46326295\t-\tENSG00000286224\tENST00000652508\tTranscript\tcoding_sequence_variant,5_prime_UTR_variant\t?-8\t?-8\t?-3\t-\t-\t-",
    hgvsc: ["ENST00000652508.1:c.3_10del", "ENST00000652508.1:c.2_10del"],
    hgvsp: [
        "ENSP00000498568.1:p.Met1_Pro3del",
        "ENSP00000498568.1:p.MetValPro1_?3",
    ],
    offset: "2",
};

/// One default-format row: `columns`, then an Extra of the three fields, each left out
/// where empty.
fn row(columns: &str, hgvsc: &str, hgvsp: &str, offset: &str) -> String {
    let mut extra = String::from("IMPACT=MODIFIER;STRAND=1");
    for (name, value) in [("HGVSc", hgvsc), ("HGVSp", hgvsp), ("HGVS_OFFSET", offset)] {
        if !value.is_empty() {
            extra.push_str(&format!(";{name}={value}"));
        }
    }
    format!("{columns}\t{extra}")
}

impl Written<'_> {
    fn reference(&self) -> String {
        row(self.columns, self.hgvsc[0], self.hgvsp[0], self.offset)
    }

    fn candidate(&self) -> String {
        row(self.columns, self.hgvsc[1], self.hgvsp[1], self.offset)
    }

    /// The reference's row with the description `description` after its protein's name.
    fn reference_describing(&self, description: &str) -> String {
        let name = self.hgvsp[0].split(":p.").next().unwrap_or_default();
        row(
            self.columns,
            self.hgvsc[0],
            &format!("{name}:p.{description}"),
            self.offset,
        )
    }

    /// The candidate's row with the description `description` after its protein's name.
    fn candidate_describing(&self, description: &str) -> String {
        let name = self.hgvsp[1].split(":p.").next().unwrap_or_default();
        row(
            self.columns,
            self.hgvsc[1],
            &format!("{name}:p.{description}"),
            self.offset,
        )
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

/// The three alleles as their records write them are the class, each row counted once;
/// the HGVSc of the two deletions is the clamped span's. Without the class the three
/// HGVSp values are differences.
#[test]
fn the_rows_of_the_three_alleles_as_written_are_the_class() {
    let alleles = [INSERTION, FRAMESHIFT_DELETION, IN_FRAME_DELETION];
    let reference: Vec<String> = alleles.iter().map(Written::reference).collect();
    let candidate: Vec<String> = alleles.iter().map(Written::candidate).collect();
    assert_eq!(
        rows_under(&[CLASS, CLAMPED_SPAN], &reference, &candidate),
        (
            "green".to_string(),
            count_map(&[(CLASS, 3), (CLAMPED_SPAN, 2)])
        )
    );
    assert_eq!(
        rows_under(&[CLAMPED_SPAN], &reference, &candidate),
        ("red".to_string(), count_map(&[(CLAMPED_SPAN, 2)]))
    );
    for allele in alleles {
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
    // beside the 5' UTR's.
    let decayed = |text: String| {
        text.replace(
            "5_prime_UTR_variant",
            "5_prime_UTR_variant,NMD_transcript_variant",
        )
    };
    assert!(takes(
        &decayed(INSERTION.reference()),
        &decayed(INSERTION.candidate())
    ));
}

/// The same alleles written at their 3' place are tested where they are written, so their
/// rows are not the class: the insertion's rows differ in their terms and in a notation
/// with no shift, and the deletions' rows agree.
#[test]
fn the_rows_of_the_alleles_written_at_their_3_prime_place_are_not_the_class() {
    let columns = |terms: &str| {
        format!(
            "inserted_after_coding_base_1\t21:15730394-15730395\tA\tENSG00000155313\tENST00000400183\tTranscript\t{terms}\t413-414\t1-2\t1\tM/NX\tatg/aAtg\t-"
        )
    };
    let reference = row(
        &columns("frameshift_variant,start_lost,start_retained_variant"),
        INSERTION.hgvsc[0],
        "ENSP00000383044.2:p.Met1?",
        "",
    );
    let candidate = row(
        &columns("frameshift_variant,start_retained_variant"),
        INSERTION.hgvsc[1],
        INSERTION.hgvsp[1],
        "",
    );
    assert!(!takes(&reference, &candidate));
    assert_eq!(
        rows_under(&[CLASS], &[reference], &[candidate]),
        ("red".to_string(), nothing())
    );
    let moved = row(
        "deleted_at_its_3_prime_place\t21:46326288-46326294\t-\tENSG00000286224\tENST00000652508\tTranscript\tframeshift_variant,start_lost\t1-7\t1-7\t1-3\tMVP/X\tATGGTCCcc/cc\t-",
        FRAMESHIFT_DELETION.hgvsc[1],
        FRAMESHIFT_DELETION.hgvsp[1],
        "",
    );
    let moved = [moved];
    assert_eq!(
        rows_under(&[CLASS], &moved, &moved),
        ("green".to_string(), nothing())
    );
}

/// The tab, VCF and JSON formats carry the cells apart and count a row once; the JSON
/// format writes the unchanged residue's `=` as it is, and leaves out the CDS start a
/// deletion's written span lacks.
#[test]
fn the_class_reads_the_tab_vcf_and_json_formats() {
    let once = ("green".to_string(), count_map(&[(CLASS, 1)]));
    let header = tab_header(&["IMPACT", "STRAND", "HGVSc", "HGVSp", "HGVS_OFFSET"]);
    for allele in [INSERTION, FRAMESHIFT_DELETION, IN_FRAME_DELETION] {
        let tab = |hgvsp: &str| {
            tab_entries(&format!(
                "{header}\n{}\tMODIFIER\t1\t{}\t{hgvsp}\t{}\n",
                allele.columns, allele.hgvsc[1], allele.offset
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
        "Allele|Consequence|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|HGVSc|HGVSp|HGVS_OFFSET",
    );
    let inserted =
        vec!["21\t15730393\tinserted_before_the_coding_start\tC\tCA\t.\t.\t.".to_string()];
    let vcf = |hgvsp: &str| {
        vcf_entries(
            &format!(
                "{csq}{}CSQ=A|5_prime_UTR_variant|ENST00000400183|Transcript|412-413|||{}|{hgvsp}|1\n",
                record_columns(&inserted[0]),
                INSERTION.hgvsc[1]
            ),
            &inserted,
        )
    };
    assert_eq!(
        under(&[CLASS], &vcf(INSERTION.hgvsp[0]), &vcf(INSERTION.hgvsp[1])),
        once
    );
    let deleted =
        vec!["21\t46326286\tdeleted_across_the_first_base\tTCATGGTC\tT\t.\t.\t.".to_string()];
    let vcf = |hgvsp: &str| {
        vcf_entries(
            &format!(
                "{csq}{}CSQ=-|coding_sequence_variant&5_prime_UTR_variant|ENST00000652508|Transcript|?-6|?-6|?-2|{}|{hgvsp}|1\n",
                record_columns(&deleted[0]),
                FRAMESHIFT_DELETION.hgvsc[1]
            ),
            &deleted,
        )
    };
    assert_eq!(
        under(
            &[CLASS],
            &vcf(FRAMESHIFT_DELETION.hgvsp[0]),
            &vcf(FRAMESHIFT_DELETION.hgvsp[1])
        ),
        once
    );

    let record = [
        ("input", text(&inserted[0])),
        ("id", text("inserted_before_the_coding_start")),
        ("allele_string", text("-/A")),
        ("most_severe_consequence", text("5_prime_UTR_variant")),
    ];
    let consequence = object(&[
        ("variant_allele", text("A")),
        ("transcript_id", text("ENST00000400183")),
        (
            "consequence_terms",
            Value::List(vec![text("5_prime_UTR_variant")]),
        ),
        ("cdna_start", Value::Int(413)),
        ("cdna_end", Value::Int(412)),
        ("hgvsc", text(INSERTION.hgvsc[1])),
        ("hgvs_offset", Value::Int(1)),
    ]);
    let object_of = |hgvsp: &str| {
        json_entries(
            &json_line(
                &record,
                with_members(&consequence, &[("hgvsp", text(hgvsp))]),
            ),
            &inserted,
        )
    };
    let unchanged = INSERTION.hgvsp[1].replace("%3D", "=");
    assert_eq!(
        under(
            &[CLASS],
            &object_of(INSERTION.hgvsp[0]),
            &object_of(&unchanged)
        ),
        once
    );
    assert_eq!(
        under(&[], &object_of(INSERTION.hgvsp[0]), &object_of(&unchanged)).0,
        "red"
    );
    let record = [
        ("input", text(&deleted[0])),
        ("id", text("deleted_across_the_first_base")),
        ("allele_string", text("CATGGTC/-")),
        ("most_severe_consequence", text("coding_sequence_variant")),
    ];
    let consequence = object(&[
        ("variant_allele", text("-")),
        ("transcript_id", text("ENST00000652508")),
        (
            "consequence_terms",
            Value::List(vec![
                text("coding_sequence_variant"),
                text("5_prime_UTR_variant"),
            ]),
        ),
        ("cdna_end", Value::Int(6)),
        ("cds_end", Value::Int(6)),
        ("protein_end", Value::Int(2)),
        ("hgvsc", text(FRAMESHIFT_DELETION.hgvsc[1])),
        ("hgvs_offset", Value::Int(1)),
    ]);
    let object_of = |members: &[(&str, Value)]| {
        json_entries(
            &json_line(&record, with_members(&consequence, members)),
            &deleted,
        )
    };
    let reference = [("hgvsp", text(FRAMESHIFT_DELETION.hgvsp[0]))];
    let candidate = [("hgvsp", text(FRAMESHIFT_DELETION.hgvsp[1]))];
    assert_eq!(
        under(&[CLASS], &object_of(&reference), &object_of(&candidate)),
        once
    );
    // An object that names a CDS start is a row tested where it is written.
    let started = |hgvsp: &str| object_of(&[("cds_start", Value::Int(1)), ("hgvsp", text(hgvsp))]);
    assert_eq!(
        under(
            &[CLASS],
            &started(FRAMESHIFT_DELETION.hgvsp[0]),
            &started(FRAMESHIFT_DELETION.hgvsp[1])
        ),
        ("red".to_string(), nothing())
    );
}

/// An empty value on either side, another residue number, a count in the reference's
/// frameshift and a candidate of another form stay differences, and so do the sides
/// exchanged and two values that name different proteins.
#[test]
fn any_other_pair_of_notations_stays_a_difference() {
    for allele in [INSERTION, FRAMESHIFT_DELETION, IN_FRAME_DELETION] {
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
        // Another residue number, on one side or on both.
        (INSERTION, "Met2AsnfsTer?", "Met1%3D"),
        (INSERTION, "Met1AsnfsTer?", "Met2%3D"),
        (INSERTION, "Met2AsnfsTer?", "Met2%3D"),
        (INSERTION, "Met11AsnfsTer?", "Met11%3D"),
        (FRAMESHIFT_DELETION, "Met2ProfsTer?", "Met2_?3"),
        (FRAMESHIFT_DELETION, "Met1ProfsTer?", "Met2_?3"),
        (IN_FRAME_DELETION, "Met2_Pro3del", "MetValPro1_?3"),
        (IN_FRAME_DELETION, "Met1_Pro3del", "MetValPro2_?3"),
        (IN_FRAME_DELETION, "Met1_Pro4del", "MetValPro1_?3"),
        (IN_FRAME_DELETION, "Met1_Pro3del", "MetValPro1_?4"),
        // A count in the reference's frameshift.
        (INSERTION, "Met1AsnfsTer5", "Met1%3D"),
        (FRAMESHIFT_DELETION, "Met1ProfsTer12", "Met1_?3"),
        // A candidate of another form, or one that names another residue.
        (INSERTION, "Met1AsnfsTer?", "Met1?"),
        (INSERTION, "Met1AsnfsTer?", "Met1Asn"),
        (INSERTION, "Met1AsnfsTer?", "Met1AsnfsTer3"),
        (INSERTION, "Met1AsnfsTer?", "Val1%3D"),
        (INSERTION, "Met1AsnfsTer?", "Met1_?3"),
        (INSERTION, "Met1AsnfsTer?", "Met1%3D%3D"),
        (FRAMESHIFT_DELETION, "Met1ProfsTer?", "Met1%3D"),
        (FRAMESHIFT_DELETION, "Met1ProfsTer?", "Val1_?3"),
        (FRAMESHIFT_DELETION, "Met1ProfsTer?", "Met1_?"),
        (FRAMESHIFT_DELETION, "Met1ProfsTer?", "Met1_Pro3del"),
        (FRAMESHIFT_DELETION, "Met1ProfsTer?", "MetValPro1_?3"),
        (IN_FRAME_DELETION, "Met1_Pro3del", "Met1_?3"),
        (IN_FRAME_DELETION, "Met1_Pro3del", "MetVal1_?3"),
        (IN_FRAME_DELETION, "Met1_Pro3del", "MetValProThr1_?3"),
        (IN_FRAME_DELETION, "Met1_Pro3del", "MetValAla1_?3"),
        (IN_FRAME_DELETION, "Met1_Pro3del", "ValValPro1_?3"),
        (IN_FRAME_DELETION, "Met1_Pro3del", "MetValPro1_Pro3del"),
        (IN_FRAME_DELETION, "Met1_Pro3del", "Met1%3D"),
        // Another reference: a deletion of one residue, of a range that starts later, a
        // replacement and the lost start itself.
        (IN_FRAME_DELETION, "Met1del", "Met1?"),
        (IN_FRAME_DELETION, "Val2_Pro3del", "MetValPro1_?3"),
        (IN_FRAME_DELETION, "Met1_Pro3delinsThr", "MetValPro1_?3"),
        (INSERTION, "Met1?", "Met1%3D"),
        (FRAMESHIFT_DELETION, "Met1?", "Met1_?3"),
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
        (INSERTION, "Met1AsnfsTer?", "Met1%3D"),
        (INSERTION, "Met1AsnfsTer?", "Met1="),
        (INSERTION, "Leu1AsnfsTer?", "Leu1%3D"),
        (FRAMESHIFT_DELETION, "Met1ProfsTer?", "Met1_?3"),
        (IN_FRAME_DELETION, "Met1_Pro3del", "MetValPro1_?3"),
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

/// A row whose terms hold a term of a start codon test, lack the 5' UTR's term or differ
/// between the engines, a row without one non-zero HGVS_OFFSET on both sides, a row that
/// is no insertion's or deletion's and a row whose CDS position has another form stay
/// differences; the class reads no other field and nothing without the candidate's row.
#[test]
fn any_other_row_stays_a_difference() {
    let alleles = [INSERTION, FRAMESHIFT_DELETION, IN_FRAME_DELETION];
    for allele in alleles {
        let (reference, candidate) = (allele.reference(), allele.candidate());
        // A term of a start codon test, on both rows.
        for term in ["start_lost", "start_retained_variant"] {
            let tested = |text: &str| {
                text.replace(
                    "5_prime_UTR_variant",
                    &format!("5_prime_UTR_variant,{term}"),
                )
            };
            assert!(!takes(&tested(&reference), &tested(&candidate)), "{term}");
            assert_eq!(
                rows_under(&[CLASS], &[tested(&reference)], &[tested(&candidate)]),
                ("red".to_string(), nothing())
            );
        }
        // No term of the 5' UTR, and terms that differ.
        let elsewhere = |text: &str| text.replace("5_prime_UTR_variant", "intron_variant");
        assert!(!takes(&elsewhere(&reference), &elsewhere(&candidate)));
        assert!(!takes(&reference, &elsewhere(&candidate)));
        assert!(!takes(&elsewhere(&reference), &candidate));
        // The offset: none, a zero, two that differ and one on one side alone.
        let offset = format!(";HGVS_OFFSET={}", allele.offset);
        let shifted = |text: &str, to: &str| text.replace(&offset, to);
        for (on_reference, on_candidate) in [
            ("", ""),
            (";HGVS_OFFSET=0", ";HGVS_OFFSET=0"),
            (offset.as_str(), ";HGVS_OFFSET=7"),
            (offset.as_str(), ""),
            ("", offset.as_str()),
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
            &shifted(&reference, ";HGVS_OFFSET=-3"),
            &shifted(&candidate, ";HGVS_OFFSET=-3")
        ));
        // Any other field of the row, and a comparison that holds no candidate row.
        let (e, a) = (default_entry(&reference), default_entry(&candidate));
        for name in ["HGVSc", "HGVS_OFFSET", "Consequence", "Extra"] {
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
    }
    // One base replaced by another is no insertion and no deletion.
    let replaced = |text: String| text.replace("21:15730393-15730394", "21:15730394");
    assert!(!takes(
        &replaced(INSERTION.reference()),
        &replaced(INSERTION.candidate())
    ));
    // An insertion's row with a CDS position is written at the first codon.
    let coding = |text: String| text.replace("\t412-413\t-\t", "\t413-414\t1-2\t");
    assert!(!takes(
        &coding(INSERTION.reference()),
        &coding(INSERTION.candidate())
    ));
    let open = |text: String| text.replace("\t412-413\t-\t", "\t?-413\t?-1\t");
    assert!(!takes(
        &open(INSERTION.reference()),
        &open(INSERTION.candidate())
    ));
    for allele in [FRAMESHIFT_DELETION, IN_FRAME_DELETION] {
        let (reference, candidate) = (allele.reference(), allele.candidate());
        // A deletion's row without the coding term, without a CDS position, with a CDS
        // start, with no CDS end, or with one of them on one side alone.
        let uncoded = |text: &str| text.replace("coding_sequence_variant,", "");
        assert!(!takes(&uncoded(&reference), &uncoded(&candidate)));
        let position = allele.columns.split('\t').nth(8).unwrap_or_default();
        let placed = |text: &str, to: &str| {
            text.replacen(
                &format!("\t{position}\t{position}\t"),
                &format!("\t{position}\t{to}\t"),
                1,
            )
        };
        for to in ["-", "1-6", "6", "2-?", "?"] {
            assert!(
                !takes(&placed(&reference, to), &placed(&candidate, to)),
                "{to}"
            );
        }
        assert!(!takes(&reference, &placed(&candidate, "1-6")));
        assert!(!takes(&placed(&reference, "1-6"), &candidate));
        // The end of the range may be any base, and the column may carry a total.
        for to in ["?-9", "?-6/1500"] {
            assert!(
                takes(&placed(&reference, to), &placed(&candidate, to)),
                "{to}"
            );
        }
    }
}

/// The lost start ends on the residue of the last coding base of the candidate's `c.`
/// deletion, and is written without a last residue where that base lies in the first
/// codon; another residue stays a difference. A notation that is no such deletion says
/// nothing of the residue.
#[test]
fn the_lost_start_ends_on_the_residue_of_the_deletions_last_base() {
    let pair = |hgvsc: &str, description: &str| {
        let reference = FRAMESHIFT_DELETION.reference();
        let candidate = FRAMESHIFT_DELETION
            .candidate_describing(description)
            .replace(FRAMESHIFT_DELETION.hgvsc[1], hgvsc);
        takes(&reference, &candidate)
    };
    for (hgvsc, description) in [
        ("ENST00000652508.1:c.1_7del", "Met1_?3"),
        ("ENST00000652508.1:c.1_5del", "Met1_?2"),
        ("ENST00000652508.1:c.1_4del", "Met1_?2"),
        ("ENST00000652508.1:c.2_17del", "Met1_?6"),
        ("ENST00000652508.1:c.-2_6del", "Met1_?2"),
        ("ENST00000652508.1:c.1_3del", "Met1?"),
        ("ENST00000652508.1:c.1_2del", "Met1?"),
        ("ENST00000652508.1:c.2del", "Met1?"),
        ("ENST00000652508.1:c.4del", "Met1_?2"),
        (
            "ENST00000652508.1:c.1_30000000000000000000000del",
            "Met1_?10000000000000000000000",
        ),
        // A span that ends past the coding sequence or in an intron names no coding base.
        ("ENST00000652508.1:c.1_*4del", "Met1_?250"),
        ("ENST00000652508.1:c.1_7+2del", "Met1_?3"),
        ("ENST00000652508.1:c.1_7delinsAC", "Met1_?3"),
        ("ENST00000652508.1:c.1_7delinsAC", "Met1?"),
    ] {
        assert!(pair(hgvsc, description), "{hgvsc} {description}");
    }
    for (hgvsc, description) in [
        ("ENST00000652508.1:c.1_7del", "Met1_?2"),
        ("ENST00000652508.1:c.1_7del", "Met1_?4"),
        ("ENST00000652508.1:c.1_7del", "Met1?"),
        ("ENST00000652508.1:c.1_6del", "Met1_?3"),
        ("ENST00000652508.1:c.1_3del", "Met1_?1"),
        ("ENST00000652508.1:c.1_3del", "Met1_?2"),
        ("ENST00000652508.1:c.2del", "Met1_?2"),
        ("ENST00000652508.1:c.1_7del", "Met1_?03"),
        ("ENST00000652508.1:c.1_7+2del", "Met1_?1"),
        ("ENST00000652508.1:c.1_7delinsAC", "Met1_?1"),
    ] {
        assert!(!pair(hgvsc, description), "{hgvsc} {description}");
    }
}

/// A pair of prediction parentheses that both values carry is set aside; one value in
/// parentheses beside one without stays a difference.
#[test]
fn a_pair_of_prediction_parentheses_on_both_values_is_set_aside() {
    for (allele, reference, candidate) in [
        (INSERTION, "(Met1AsnfsTer?)", "(Met1%3D)"),
        (INSERTION, "(Met1AsnfsTer?)", "(Met1=)"),
        (FRAMESHIFT_DELETION, "(Met1ProfsTer?)", "(Met1_?3)"),
        (IN_FRAME_DELETION, "(Met1_Pro3del)", "(MetValPro1_?3)"),
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
        &INSERTION.reference_describing("(Met1AsnfsTer?"),
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
        takes(
            &written("ENSP00000383044.2:p.Met1IlefsTer?"),
            &written(INSERTION.hgvsp[1]),
        )
    };
    assert!(pair("ENST00000400183.7:c.1_2dup"));
    assert!(pair("ENST00000400183.7:c.2_3insAT"));
    assert!(!pair("ENST00000400183.7:c.1_2delinsTA"));
    assert!(!pair("ENST00000999999.1:c.1_2dup"));
    assert!(!pair(""));
    // A row whose first column spells the record's alleles says so itself.
    let spelt = columns.replace("inserted_before_the_coding_start", "21_15730394_-/AT");
    let written = |hgvsp: &str| row(&spelt, "", hgvsp, "2");
    assert!(takes(
        &written("ENSP00000383044.2:p.Met1IlefsTer?"),
        &written(INSERTION.hgvsp[1])
    ));
}

/// The class reads the HGVSp under either letter case and stands in the class table
/// before the class of a start codon kept without a 5' UTR, the class of a start codon
/// test past the shifted allele, the class of the LoFTEE flags and the two classes that
/// close the table; the `--masks` help names it.
#[test]
fn the_class_tables_and_the_help_carry_the_class() {
    assert_eq!(classes::HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_START_TEST, CLASS);
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
    assert_eq!(names[names.len() - 6], CLASS);
    for table in [
        tables.in_term_pair_following_classes(CLASS),
        tables.in_record_reading_classes(CLASS),
        tables.in_transcript_reading_classes(CLASS),
        tables.in_sibling_allele_classes(CLASS),
    ] {
        assert!(!table);
    }
    let help = include_str!("../data/help.txt");
    assert_eq!(help.matches(CLASS).count(), 1);
}
