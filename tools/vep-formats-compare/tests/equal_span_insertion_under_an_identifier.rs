// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `refseq_edited_model_equal_span_insertion_term` on a default or tab row named by an
//! identifier, through the binary. Such a row's first column spells no alleles. The class is
//! one of `RECORD_REFERENCE_CLASSES`: a comparison that is given the input reads from its
//! records the bases the row's allele stands in for, and one that is not reads the row from
//! its Location and its allele.

mod support;

use std::path::{Path, PathBuf};

use support::args_of;
use support::compact;
use support::member;
use support::report;
use support::run;
use support::scratch;
use support::Ran;
use vep_formats_compare::classes;
use vep_formats_compare::pyjson::Value;

const CLASS: &str = classes::REFSEQ_EDITED_MODEL_EQUAL_SPAN_INSERTION_TERM;

/// A made-up record under an identifier and its row on one edited model: the VCF line's
/// position and alleles, the row's Allele, and a CDS span that holds as many bases as the
/// Allele. Every row's Location is the record's position and the one after it, an
/// insertion's after its shared first base and a replacement's of two bases alike.
struct Record {
    name: &'static str,
    position: u32,
    reference: &'static str,
    alternate: &'static str,
    allele: &'static str,
    cdna: &'static str,
    cds: &'static str,
    protein: &'static str,
}

/// Two bases inserted after the record's first base.
const TWO_INSERTED: Record = Record {
    name: "two_inserted",
    position: 73385942,
    reference: "A",
    alternate: "AGG",
    allele: "GG",
    cdna: "186-187",
    cds: "75-76",
    protein: "25-26",
};

/// Two bases that replace two.
const TWO_REPLACED: Record = Record {
    name: "two_replaced",
    position: 73385950,
    reference: "CT",
    alternate: "GA",
    allele: "GA",
    cdna: "196-197",
    cds: "85-86",
    protein: "29",
};

/// Three bases inserted after the record's first base.
const THREE_INSERTED: Record = Record {
    name: "three_inserted",
    position: 73385960,
    reference: "A",
    alternate: "AGGC",
    allele: "GGC",
    cdna: "206-208",
    cds: "95-97",
    protein: "32-33",
};

/// Three bases that replace two.
const THREE_FOR_TWO: Record = Record {
    name: "three_for_two",
    position: 73385970,
    reference: "CT",
    alternate: "GAC",
    allele: "GAC",
    cdna: "216-218",
    cds: "105-107",
    protein: "35-36",
};

const COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation";

impl Record {
    fn location(&self) -> String {
        format!("2:{}-{}", self.position, self.position + 1)
    }

    /// The record's row with the Consequence `term` and the IMPACT `impact`: the default
    /// format carries the IMPACT and the edit flag in `Extra`, the tab format as columns.
    fn row(&self, fmt: &str, term: &str, impact: &str) -> String {
        let columns = format!(
            "{}\t{}\t{}\t7840\tNM_015120.4\tTranscript\t{term}\t{}\t{}\t{}\t-\t-\t-",
            self.name,
            self.location(),
            self.allele,
            self.cdna,
            self.cds,
            self.protein
        );
        match fmt {
            "default" => format!("{columns}\tIMPACT={impact};STRAND=1;BAM_EDIT=OK"),
            _ => format!("{columns}\t{impact}\t1\tOK"),
        }
    }

    fn vcf_line(&self) -> String {
        format!(
            "2\t{}\t{}\t{}\t{}\t.\t.\t.",
            self.position, self.name, self.reference, self.alternate
        )
    }

    /// The row's key in a report.
    fn key(&self) -> String {
        format!("{}|{}|NM_015120.4|Transcript", self.location(), self.allele)
    }
}

/// An output of `records` in the format `fmt`, each row with `term` and `impact`.
fn output(fmt: &str, records: &[&Record], term: &str, impact: &str) -> String {
    let header = match fmt {
        "default" => format!("{COLUMNS}\tExtra"),
        _ => format!("{COLUMNS}\tIMPACT\tSTRAND\tBAM_EDIT"),
    };
    let mut text = format!("{header}\n");
    for record in records {
        text.push_str(&record.row(fmt, term, impact));
        text.push('\n');
    }
    text
}

fn planted(directory: &Path, name: &str, text: &str) -> PathBuf {
    let path = directory.join(name);
    std::fs::write(&path, text).expect("the file is written");
    path
}

fn path_text(path: &Path) -> &str {
    path.to_str().expect("a UTF-8 path")
}

/// The run of the binary with its report. The reference's rows of `records` carry
/// `intergenic_variant` with its IMPACT and the candidate's `missense_variant` with its
/// own, compared in the format `fmt` under the class alone, beside the input VCF `input`
/// where one is given.
fn compared(test: &str, fmt: &str, records: &[&Record], input: Option<&str>) -> (Ran, Value) {
    let directory = scratch(&format!(
        "{test}-{fmt}-{}",
        if input.is_some() {
            "beside-the-input"
        } else {
            "without-the-input"
        }
    ));
    let expected = planted(
        &directory,
        "expected.txt",
        &output(fmt, records, "intergenic_variant", "MODIFIER"),
    );
    let actual = planted(
        &directory,
        "actual.txt",
        &output(fmt, records, "missense_variant", "MODERATE"),
    );
    let report_path = directory.join("report.json");
    let mut arguments = vec![
        "--format",
        fmt,
        "--expected",
        path_text(&expected),
        "--actual",
        path_text(&actual),
        "--masks",
        CLASS,
        "--report",
        path_text(&report_path),
    ];
    let input = input.map(|text| planted(&directory, "input.vcf", text));
    if let Some(input) = &input {
        arguments.extend(["--input", path_text(input)]);
    }
    let ran = run(&directory, &args_of(&arguments));
    let read = report(&report_path);
    (ran, read)
}

/// An input VCF of the lines `lines`.
fn input_of(lines: &[String]) -> String {
    let mut text =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_string();
    for line in lines {
        text.push_str(line);
        text.push('\n');
    }
    text
}

/// The input VCF that holds the lines of `records`.
fn records_of(records: &[&Record]) -> String {
    let lines: Vec<String> = records.iter().map(|record| record.vcf_line()).collect();
    input_of(&lines)
}

/// The fields in which a row of the format `fmt` differs: the Consequence and the IMPACT,
/// and in the default format the `Extra` text that carries the IMPACT.
fn differing_fields(fmt: &str) -> &'static [&'static str] {
    match fmt {
        "default" => &["Consequence", "Extra", "IMPACT"],
        _ => &["Consequence", "IMPACT"],
    }
}

/// What the binary prints for a red comparison in which `count` rows differ.
fn red(fmt: &str, count: usize) -> String {
    let mut printed = format!("{fmt}: red\n");
    for field in differing_fields(fmt) {
        printed.push_str(&format!("  field {field}: {count} differing entries\n"));
    }
    printed
}

/// The count the report holds for the class under `masks.applied`, 0 where it holds none.
fn applied(read: &Value) -> i64 {
    match member(read, &["masks", "applied"]).get(CLASS) {
        Some(Value::Int(count)) => *count,
        Some(other) => panic!("masks.applied holds {other:?} for the class"),
        None => 0,
    }
}

/// Whether the report's examples of a differing Consequence name the row of `record`.
fn left_a_difference(read: &Value, record: &Record) -> bool {
    compact(member(read, &["columns", "Consequence"])).contains(&format!(
        "{}: expected 'intergenic_variant', got 'missense_variant'",
        record.key()
    ))
}

/// Beside the input, an insertion under an identifier is the class's whatever its length:
/// its record's reference allele is empty once the shared first base is dropped. Without
/// the input the row of two inserted bases cannot be told from two replaced bases and stays
/// a difference.
#[test]
fn an_insertion_under_an_identifier_is_read_from_the_records_of_the_input() {
    let records = [&TWO_INSERTED, &THREE_INSERTED];
    for fmt in ["default", "tab"] {
        let (ran, read) = compared("inserted", fmt, &records, Some(&records_of(&records)));
        assert_eq!(ran.out(), format!("{fmt}: green\n"), "{fmt}");
        assert_eq!((ran.code, ran.err()), (0, String::new()), "{fmt}");
        assert_eq!(applied(&read), 2, "{fmt}");
        assert_eq!(compact(member(&read, &["reasons"])), "[]", "{fmt}");
        assert_eq!(
            compact(member(&read, &["columns", "Consequence", "masked"])),
            format!("{{\"{CLASS}\": 2}}"),
            "{fmt}"
        );
        assert_eq!(compact(member(&read, &["masks", "unused"])), "[]", "{fmt}");

        let (ran, read) = compared("inserted", fmt, &records, None);
        assert_eq!(ran.out(), red(fmt, 1), "{fmt}");
        assert_eq!((ran.code, ran.err()), (1, String::new()), "{fmt}");
        assert_eq!(applied(&read), 1, "{fmt}");
        assert!(left_a_difference(&read, &TWO_INSERTED), "{fmt}");
        assert!(!left_a_difference(&read, &THREE_INSERTED), "{fmt}");
    }
}

/// Beside the input, bases a record replaces are no insertion, whatever the allele's length,
/// and the class takes no such row. Without the input two bases replaced by three have an
/// insertion's Location and allele and are taken with the insertions, the limit of a row
/// that spells no alleles.
#[test]
fn bases_a_record_replaces_under_an_identifier_stay_differences_beside_the_input() {
    let records = [&TWO_REPLACED, &THREE_FOR_TWO];
    for fmt in ["default", "tab"] {
        let (ran, read) = compared("replaced", fmt, &records, Some(&records_of(&records)));
        assert_eq!(ran.out(), red(fmt, 2), "{fmt}");
        assert_eq!((ran.code, ran.err()), (1, String::new()), "{fmt}");
        assert_eq!(applied(&read), 0, "{fmt}");
        assert_eq!(
            compact(member(&read, &["masks", "unused"])),
            format!("[\"{CLASS}\"]"),
            "{fmt}"
        );
        assert!(left_a_difference(&read, &TWO_REPLACED), "{fmt}");
        assert!(left_a_difference(&read, &THREE_FOR_TWO), "{fmt}");

        let (ran, read) = compared("replaced", fmt, &records, None);
        assert_eq!(ran.out(), red(fmt, 1), "{fmt}");
        assert_eq!((ran.code, ran.err()), (1, String::new()), "{fmt}");
        assert_eq!(applied(&read), 1, "{fmt}");
        assert!(left_a_difference(&read, &TWO_REPLACED), "{fmt}");
        assert!(!left_a_difference(&read, &THREE_FOR_TWO), "{fmt}");
    }
}

/// A row no record of the input places answers as it does without the input: the input
/// here holds the records of another chromosome alone.
#[test]
fn a_row_no_record_of_the_input_places_answers_from_its_location_and_allele() {
    let records = [&TWO_INSERTED, &THREE_INSERTED];
    let elsewhere: Vec<String> = records
        .iter()
        .map(|record| record.vcf_line().replacen("2\t", "3\t", 1))
        .collect();
    for fmt in ["default", "tab"] {
        let (ran, read) = compared("unplaced", fmt, &records, Some(&input_of(&elsewhere)));
        assert_eq!(ran.out(), red(fmt, 1), "{fmt}");
        assert_eq!((ran.code, ran.err()), (1, String::new()), "{fmt}");
        assert_eq!(applied(&read), 1, "{fmt}");
        assert!(left_a_difference(&read, &TWO_INSERTED), "{fmt}");
        assert!(!left_a_difference(&read, &THREE_INSERTED), "{fmt}");
    }
}

/// The class stands in the table of the classes that read a row's reference bases from the
/// records of the input.
#[test]
fn the_class_is_one_of_the_classes_that_read_the_records_of_the_input() {
    let tables = classes::tables();
    assert!(tables.in_record_reference_classes(CLASS));
}
