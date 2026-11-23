// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A gzip `--input` under a structural-variant class.
//!
//! The comparison reads every member of the file, zero bytes between members being
//! padding. The class reads the file as `zcat` writes it: up to the first zero byte after
//! a member, and no member whose header `zcat` refuses. Where the two readings differ, the
//! run ends as these cases have it.

mod support;

use flate2::write::GzEncoder;
use std::io::Write as _;
use std::path::Path;
use support::Ran;

const HEAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
/// A deletion the reference annotates.
const KEPT: &str = "21\t100\tdel1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=300\n";
/// An insertion the reference drops before annotation for the length its INFO states.
const DROPPED: &str = "21\t1000\tbig1\tN\t<INS>\t.\t.\tSVTYPE=INS;SVLEN=20000000\n";

const COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t\
                       Consequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\t\
                       Codons\tExisting_variation\tExtra\n";

const COUNTS_DIFFER: &str = "records the structural-variant comparator reads; a record needs \
                             CHROM, POS, ID, REF and ALT\n";

fn row(name: &str, location: &str, allele: &str, feature: &str) -> String {
    format!(
        "{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t\
         -\tIMPACT=MODIFIER;STRAND=1\n"
    )
}

fn member(text: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(text).expect("the text deflates");
    encoder.finish().expect("the member ends")
}

/// A member of `text` whose header carries the flag `flag`.
fn flagged(text: &[u8], flag: u8) -> Vec<u8> {
    let mut member = member(text);
    member[3] |= flag;
    member
}

/// A member of `text` with a header checksum: the header's own, or another where `wrong`.
fn with_header_checksum(text: &[u8], wrong: bool) -> Vec<u8> {
    let mut member = flagged(text, 0x02);
    let mut crc = flate2::Crc::new();
    crc.update(&member[..10]);
    let checksum = (crc.sum() as u16) ^ u16::from(wrong);
    member.splice(10..10, checksum.to_le_bytes());
    member
}

/// A default-format comparison under `reference_skipped_record` beside the input `input`:
/// the reference holds the row of the kept record, the candidate that row and two rows
/// of the dropped record.
fn compared(name: &str, input: &[u8]) -> (Ran, bool) {
    let directory = support::scratch(&format!("structural-gzip-input-{name}"));
    let plant = |file: &str, bytes: &[u8]| {
        std::fs::write(directory.join(file), bytes).expect("a file is planted");
    };
    let kept = row("del1", "21:101-300", "deletion", "ENST1");
    let dropped: String = ["ENST1", "ENST1b"]
        .iter()
        .map(|feature| row("big1", "21:1001-20001000", "insertion", feature))
        .collect();
    plant("e.txt", format!("{COLUMNS}{kept}").as_bytes());
    plant("a.txt", format!("{COLUMNS}{kept}{dropped}").as_bytes());
    plant("in.vcf.gz", input);
    let ran = support::run(
        &directory,
        &[
            "--format",
            "default",
            "--expected",
            "e.txt",
            "--actual",
            "a.txt",
            "--masks",
            "reference_skipped_record",
            "--input",
            "in.vcf.gz",
            "--report",
            "report.json",
        ],
    );
    let written = Path::new(&directory).join("report.json").exists();
    (ran, written)
}

fn first() -> Vec<u8> {
    member(format!("{HEAD}{KEPT}").as_bytes())
}

fn second() -> Vec<u8> {
    member(DROPPED.as_bytes())
}

/// How a run that is refused for its two counts of records ends.
fn refused_for(ran: &(Ran, bool), lines: usize, records: usize) {
    let (ran, written) = ran;
    assert_eq!(
        ran.err(),
        format!(
            "ERROR: [compare_vep_formats] in.vcf.gz: {lines} data lines, {records} {COUNTS_DIFFER}"
        )
    );
    assert_eq!((ran.code, ran.out().as_str(), *written), (2, "", false));
}

#[test]
fn members_are_read_one_after_another_and_past_zero_bytes_after_the_last() {
    let whole = member(format!("{HEAD}{KEPT}{DROPPED}").as_bytes());
    let inputs = [
        ("one", whole),
        ("two", [first(), second()].concat()),
        ("padded", [first(), second(), vec![0; 8]].concat()),
        // Header lines alone stand after the zero bytes: both readings count two records.
        (
            "header-lines",
            [first(), second(), vec![0; 3], member(b"#late\n")].concat(),
        ),
    ];
    for (name, input) in inputs {
        let (ran, written) = compared(name, &input);
        assert_eq!(
            (ran.code, ran.out().as_str(), ran.err().as_str(), written),
            (0, "default: green\n", "", true),
            "{name}"
        );
    }
}

#[test]
fn a_record_after_zero_bytes_is_not_read_and_the_two_counts_of_records_are_refused() {
    for (name, padding) in [("one-zero", 1), ("eight-zeros", 8)] {
        let input = [first(), vec![0; padding], second()].concat();
        refused_for(&compared(name, &input), 2, 1);
    }
}

#[test]
fn a_line_that_zero_bytes_cut_is_a_record_of_the_columns_before_them() {
    let cut = DROPPED.find("\t.\t.\t").expect("the columns after ALT");
    let input = [
        member(format!("{HEAD}{KEPT}{}", &DROPPED[..cut]).as_bytes()),
        vec![0; 2],
        member(&DROPPED.as_bytes()[cut..]),
    ]
    .concat();
    // Without its INFO the record is one the reference annotates, so its rows count.
    let (ran, written) = compared("cut-line", &input);
    assert_eq!(
        (ran.code, ran.out().as_str(), ran.err().as_str(), written),
        (
            1,
            "default: red\n  unexpected in the candidate: 2 entries\n",
            "",
            true
        )
    );
}

#[test]
fn a_character_that_zero_bytes_cut_ends_the_run() {
    let noted = format!("{}NOTE=caf\u{e9}\n", DROPPED.replace('\n', ";"));
    // The two bytes of the character before the line's end stand on either side of the
    // zero byte.
    let (before, after) = noted.as_bytes().split_at(noted.len() - 2);
    let input = [
        member(&[HEAD.as_bytes(), KEPT.as_bytes(), before].concat()),
        vec![0],
        member(after),
    ]
    .concat();
    let (ran, written) = compared("cut-character", &input);
    assert_eq!(
        ran.err(),
        "ERROR: [compare_vep_formats] 'utf-8' codec can't decode byte 0xc3 in position 0: \
         unexpected end of data\n"
    );
    assert_eq!((ran.code, ran.out().as_str(), written), (1, "", false));
}

#[test]
fn a_member_whose_header_is_refused_ends_the_reading_before_it() {
    let whole = format!("{HEAD}{KEPT}{DROPPED}");
    let (ran, written) = compared("checksum", &with_header_checksum(whole.as_bytes(), false));
    assert_eq!(
        (ran.code, ran.out().as_str(), ran.err().as_str(), written),
        (0, "default: green\n", "", true)
    );
    let wrong = with_header_checksum(whole.as_bytes(), true);
    refused_for(&compared("wrong-checksum", &wrong), 2, 0);
    let later = [first(), with_header_checksum(DROPPED.as_bytes(), true)].concat();
    refused_for(&compared("later-wrong-checksum", &later), 2, 1);
    for (name, flag) in [
        ("encrypted", 0x20),
        ("reserved", 0x40),
        ("reserved-high", 0x80),
    ] {
        refused_for(&compared(name, &flagged(whole.as_bytes(), flag)), 2, 0);
        let later = [first(), flagged(DROPPED.as_bytes(), flag)].concat();
        refused_for(&compared(&format!("later-{name}"), &later), 2, 1);
    }
}
