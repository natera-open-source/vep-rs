// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The input index and the VCF parser, on the committed Ensembl VEP outputs of
//! `tests/golden/116/GRCh38-release116` and on planted lines with their expected values.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::{Error, Uncaught};
use vep_formats_compare::input::{input_records, read_text};
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::vcf::*;

const RECORDS: usize = 36;
const ROWS: usize = 977;
const SKIPPED_RECORDS: [i64; 3] = [12, 13, 15];

const CSQ_HEADER: &str = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence \
                          annotations from Ensembl VEP. Format: \
                          Allele|Consequence|Feature_type|Feature\">";

fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-release116")
}

fn strings(items: &[&str]) -> Vec<String> {
    items.iter().map(|item| item.to_string()).collect()
}

fn keys(entries: &[Entry]) -> Vec<&str> {
    entries.iter().map(|entry| entry.key.as_str()).collect()
}

fn fields(entry: &Entry) -> Vec<(&str, &str)> {
    entry
        .fields
        .iter()
        .map(|(name, value)| (name.as_str(), value.as_str()))
        .collect()
}

fn comparator_error(result: Result<Vec<Entry>, Error>) -> String {
    match result {
        Err(Error::Comparator(text)) => text,
        Ok(entries) => panic!("{} entries and no error", entries.len()),
    }
}

/// What the run ends with on an input no reading of the comparison expects.
fn uncaught<T>(run: impl FnOnce() -> T + std::panic::UnwindSafe) -> String {
    match std::panic::catch_unwind(run) {
        Ok(_) => panic!("the run did not end"),
        Err(payload) => match payload.downcast_ref::<Uncaught>() {
            Some(Uncaught(what)) => what.clone(),
            None => std::panic::resume_unwind(payload),
        },
    }
}

#[test]
fn the_committed_reference_vcf_gives_a_record_entry_a_line_and_an_entry_a_csq_chunk() {
    let text = read_text(&corpus().join("expected/vcf.vcf.gz")).unwrap();
    let inputs = input_records(&corpus().join("variants.vcf")).unwrap();
    let entries = parse_vcf(&text, &inputs).unwrap();
    let (records, csq): (Vec<&Entry>, Vec<&Entry>) = entries
        .iter()
        .partition(|entry| entry.key.ends_with("|record"));
    assert_eq!(records.len(), RECORDS);
    assert_eq!(csq.len(), ROWS);
    let named: BTreeSet<i64> = records.iter().filter_map(|entry| entry.record).collect();
    let annotated: BTreeSet<i64> = (0..39)
        .filter(|record| !SKIPPED_RECORDS.contains(record))
        .collect();
    assert_eq!(named, annotated);
    assert!(records.iter().all(|entry| entry.record.is_some()));
    let ins = records
        .iter()
        .find(|entry| entry.fields["ID"] == "multi_ins_ins")
        .unwrap();
    assert_eq!(ins.key, "rec11|record");
    assert_eq!(ins.fields["INFO_without_CSQ"], "END=30165780");
    assert_eq!(ins.fields["CSQ_present"], "yes");
    let chunks: Vec<&&Entry> = csq
        .iter()
        .filter(|entry| entry.record == Some(11))
        .collect();
    assert_eq!(
        chunks
            .iter()
            .map(|entry| entry.key.as_str())
            .collect::<Vec<_>>(),
        [
            "rec11|insertion|ENST00000286808|Transcript",
            "rec11|insertion|ENST00000286808|Transcript#2",
        ]
    );
    assert_eq!(
        chunks[0].consequence,
        "feature_elongation&3_prime_UTR_variant"
    );
    assert_eq!(chunks[0].fields["EXON"], "1/1");
    assert_eq!(chunks[0].fields["cDNA_position"], "1025-1026");
    assert_eq!(chunks[0].fields["BIOTYPE"], "protein_coding");
    assert_eq!(chunks[0].fields["HGVSc"], "");
}

#[test]
fn an_input_line_gives_its_record_in_the_vcf_and_in_the_ensembl_format() {
    let record = |chromosome: &str, position: i64, reference: &str, alternates: &[&str]| {
        Some((
            chromosome.to_string(),
            position,
            reference.to_string(),
            strings(alternates),
        ))
    };
    assert_eq!(
        input_line_record("1\t100\trec1\tG\tGTAT,C\t.\t.\t."),
        record("1", 100, "G", &["GTAT", "C"])
    );
    assert_eq!(
        input_line_record("1\t100\t7\tg\tgtat"),
        record("1", 100, "G", &["GTAT"])
    );
    assert_eq!(
        input_line_record("1\t101\t100\t-/TAT\t+\trec1"),
        record("1", 101, "", &["TAT"])
    );
    assert_eq!(
        input_line_record("1 101 100 -/ATA -"),
        record("1", 101, "", &["TAT"])
    );
    assert_eq!(
        input_line_record("1\t100\t100\tC/ATAC/G\t-1\trec1"),
        record("1", 100, "G", &["GTAT", "C"])
    );
    assert_eq!(
        input_line_record("1\t100\t100\tg/gtat\t1"),
        record("1", 100, "G", &["GTAT"])
    );
    for line in [
        "1\t101\t100\t-/TAT",
        "1\tx\t100\t-/TAT\t+",
        "1\t101\tx\t-/TAT\t+",
        "1\t101\t100\t-/TAT\t0",
        "1\t100\trec1\t.\tGTAT",
        "",
    ] {
        assert_eq!(input_line_record(line), None, "{line:?}");
    }
    let inputs = strings(&[
        "1\t101\t100\t-/TAT\t+\trec1",
        "1\t200\t200\tG/GTAT\t+\trec2",
        "1\t300\trec3\tG\tGTAT",
    ]);
    let ordinals: Vec<Option<i64>> = ["100", "101", "99", "102", "x"]
        .iter()
        .map(|position| ensembl_line_ordinal(&inputs, 0, "1", position))
        .collect();
    assert_eq!(ordinals, [Some(0), Some(0), None, None, None]);
    assert_eq!(ensembl_line_ordinal(&inputs, 1, "1", "199"), Some(1));
    assert_eq!(ensembl_line_ordinal(&inputs, 1, "2", "199"), None);
    assert_eq!(ensembl_line_ordinal(&inputs, 2, "1", "300"), None);
    assert_eq!(ensembl_line_ordinal(&inputs, 3, "1", "300"), None);
}

#[test]
fn a_malformed_line_is_an_error_in_every_format() {
    let inputs = input_records(&corpus().join("variants.vcf")).unwrap();
    assert!(parse_default("#A\tB\n1\t2\t3\n").is_err());
    assert!(parse_vcf(
        "##INFO=<ID=CSQ,Description=\"Format: A|B\">\n1\t2\t.\tA\tT\t.\t.\tCSQ=x\n",
        &inputs
    )
    .is_err());
    assert!(parse_json("[1, 2]\n", &inputs).is_err());
}

#[test]
fn identical_input_lines_are_matched_in_order_of_appearance() {
    let inputs = strings(&["a\t1", "b\t2", "a\t1", "c"]);
    let mut index = InputIndex::new(&inputs, &mut vcf_line_key);
    assert_eq!(index.record("a\t1"), Some(0));
    assert_eq!(index.record("c"), Some(3));
    assert_eq!(index.record("a\t1"), Some(2));
    assert_eq!(index.record("a\t1"), None);
    assert_eq!(index.record("d"), None);
    assert_eq!(index.record("b\t2"), Some(1));
    let mut by_first_letter = InputIndex::new(&inputs, &mut |line: &str| line[..1].to_string());
    assert_eq!(by_first_letter.record("a"), Some(0));
    assert_eq!(by_first_letter.record("a"), Some(2));
    assert_eq!(by_first_letter.record("a\t1"), None);
}

#[test]
fn a_line_is_named_by_its_first_five_columns() {
    assert_eq!(vcf_line_key("1\t2\t3\t4\t5\t6\t7"), "1\t2\t3\t4\t5");
    assert_eq!(vcf_line_key("1\t2\t3\t4\t5"), "1\t2\t3\t4\t5");
    assert_eq!(vcf_line_key("1\t2\t3\t4\t5\t"), "1\t2\t3\t4\t5");
    assert_eq!(vcf_line_key("1\t2"), "1\t2");
    assert_eq!(vcf_line_key(""), "");
}

#[test]
fn the_csq_header_names_the_fields_of_its_format_list() {
    assert_eq!(
        csq_format_fields(CSQ_HEADER),
        Ok(strings(&[
            "Allele",
            "Consequence",
            "Feature_type",
            "Feature"
        ]))
    );
    assert_eq!(
        csq_format_fields("##INFO=<ID=CSQ,Description=\"x Format: A|B Format: C\">\">"),
        Ok(strings(&["A", "B Format: C"]))
    );
    assert_eq!(
        csq_format_fields("##INFO=<ID=CSQ,Description=\"Format: A|"),
        Ok(strings(&["A", ""]))
    );
    assert_eq!(
        csq_format_fields("##INFO=<ID=CSQ,Description=\"none\">"),
        Err(Error::Comparator(
            "CSQ header without a Format list: ##INFO=<ID=CSQ,Description=\"none\">".to_string()
        ))
    );
}

#[test]
fn an_input_line_on_the_forward_strand_gives_its_record() {
    assert_eq!(
        input_line_record("1\t100\trec1\tG\tGTAT,C\t.\t.\t."),
        Some((
            "1".to_string(),
            100,
            "G".to_string(),
            strings(&["GTAT", "C"])
        ))
    );
    assert_eq!(
        input_line_record("1\t101\t100\t-/tat\t+1\trec1"),
        Some(("1".to_string(), 101, String::new(), strings(&["TAT"])))
    );
    // A position in the digits of another script is the number it writes.
    assert_eq!(
        input_line_record("1\t\u{ff16}\u{ff10}\t\u{ff16}\u{ff10}\tC/G\t+"),
        Some(("1".to_string(), 60, "C".to_string(), strings(&["G"])))
    );
    for line in [
        "1\t101\t100\t-/TAT",
        "1\t101\t100\t-/TAT\t0",
        "1\t100\tid\t.\tT",
        "",
    ] {
        assert_eq!(input_line_record(line), None, "{line:?}");
    }
}

#[test]
fn an_output_line_is_written_from_the_ensembl_format_line_at_its_ordinal() {
    let inputs = strings(&[
        "1\t101\t100\t-/TAT\t+\trec1",
        "1\t200\t200\tG/GTAT\t+\trec2",
        "1\t300\trec3\tG\tGTAT",
    ]);
    assert_eq!(ensembl_line_ordinal(&inputs, 0, "1", "100"), Some(0));
    assert_eq!(ensembl_line_ordinal(&inputs, 0, "1", "101"), Some(0));
    assert_eq!(ensembl_line_ordinal(&inputs, 0, "1", "99"), None);
    assert_eq!(ensembl_line_ordinal(&inputs, 0, "1", "x"), None);
    assert_eq!(ensembl_line_ordinal(&inputs, 1, "2", "199"), None);
    assert_eq!(ensembl_line_ordinal(&inputs, 2, "1", "300"), None);
    assert_eq!(ensembl_line_ordinal(&inputs, 3, "1", "300"), None);
    assert_eq!(ensembl_line_ordinal(&inputs, -1, "1", "300"), None);
    // A position longer than `i64` holds is neither the line's start nor the one below it.
    assert_eq!(
        ensembl_line_ordinal(&inputs, 0, "1", "100000000000000000000"),
        None
    );
}

#[test]
fn a_position_the_interpreter_does_not_convert_ends_the_run() {
    let inputs = strings(&["1\t2\t2\tG/T\t+", "1\t\u{b2}\t2\tG/T\t+"]);
    assert_eq!(
        uncaught(|| ensembl_line_ordinal(&inputs, 0, "1", "\u{b2}")),
        "invalid literal for int() with base 10: '\u{b2}'"
    );
    assert_eq!(
        uncaught(|| ensembl_line_ordinal(&inputs, 1, "1", "2")),
        "invalid literal for int() with base 10: '\u{b2}'"
    );
    assert_eq!(
        uncaught(|| input_line_record("1\t\u{b2}\tid\tG\tT")),
        "invalid literal for int() with base 10: '\u{b2}'"
    );
    let long = "9".repeat(4301);
    assert_eq!(
        uncaught(|| ensembl_line_ordinal(&inputs, 0, "1", &long)),
        "Exceeds the limit (4300 digits) for integer string conversion: value has 4301 digits; \
         use sys.set_int_max_str_digits() to increase the limit"
    );
    assert_eq!(
        uncaught(|| input_line_record("1\t100000000000000000000\tid\tG\tT")),
        "invalid literal for int() with base 10: '100000000000000000000'"
    );
}

#[test]
fn a_line_of_fewer_than_eight_columns_is_an_error() {
    let lines = strings(&[CSQ_HEADER, "1\t100\trs1\tG\tGTAT\t.\t."]);
    assert_eq!(
        comparator_error(parse_vcf_lines(&lines, &[0], &[])),
        "VCF line has fewer than 8 columns: 1\t100\trs1\tG\tGTAT\t.\t."
    );
    assert_eq!(
        comparator_error(parse_vcf_lines(&strings(&[""]), &[0], &[])),
        "VCF line has fewer than 8 columns: "
    );
}

#[test]
fn a_csq_chunk_of_another_width_than_the_format_list_is_an_error() {
    let line = "1\tx\t.\tA\tT\t.\t.\tCSQ=x";
    assert_eq!(
        comparator_error(parse_vcf_lines(&strings(&[line]), &[4], &[])),
        "CSQ chunk width on line4: x"
    );
    let lines = strings(&[CSQ_HEADER, "1\tx\t.\tA\tT\t.\t.\tCSQ=T|m|Transcript|T1,T|m"]);
    assert_eq!(
        comparator_error(parse_vcf_lines(&lines, &[0], &strings(&["1\tx\t.\tA\tT"]))),
        "CSQ chunk width on rec0: T|m"
    );
}

#[test]
fn a_csq_header_without_a_format_list_is_an_error() {
    let lines = strings(&["##INFO=<ID=CSQ,Description=\"none\">"]);
    assert_eq!(
        comparator_error(parse_vcf_lines(&lines, &[], &[])),
        "CSQ header without a Format list: ##INFO=<ID=CSQ,Description=\"none\">"
    );
}

#[test]
fn a_data_line_without_an_ordinal_is_an_error() {
    let lines = strings(&[CSQ_HEADER, "1\t300\t.\tA\tT\t.\t.\t."]);
    assert_eq!(
        comparator_error(parse_vcf_lines(&lines, &[], &[])),
        "a data line has no ordinal in its output"
    );
}

#[test]
fn a_line_that_matches_no_input_record_is_named_by_its_ordinal_in_the_whole_output() {
    let inputs = strings(&[
        "1\t100\trs1\tG\tGTAT,C\t.\t.\t.",
        "1\t200\t.\tAT\tA\t.\t.\t.",
    ]);
    let lines = strings(&[
        "##fileformat=VCFv4.2",
        CSQ_HEADER,
        "1\t200\t.\tAT\tA\t.\t.\tAC=1",
        "1\t300\t.\tA\tT\t.\t.\t.",
        "#a comment between the lines",
        "1\t300\t.\tA\tT\t.\t.\t.",
    ]);
    let entries = parse_vcf_lines(&lines, &[4, 7, 9], &inputs).unwrap();
    assert_eq!(
        keys(&entries),
        ["rec1|record", "line7|record", "line9|record"]
    );
    assert_eq!(
        entries.iter().map(|entry| entry.record).collect::<Vec<_>>(),
        [Some(1), None, None]
    );
    let text = format!("{}\n", lines.join("\n"));
    let entries = parse_vcf(&text, &inputs).unwrap();
    assert_eq!(
        keys(&entries),
        ["rec1|record", "line1|record", "line2|record"]
    );
}

#[test]
fn a_line_gives_a_record_entry_and_an_entry_a_chunk_of_its_last_csq() {
    let inputs = strings(&[
        "1\t100\trs1\tG\tGTAT,C\t.\t.\t.",
        "1\t200\t.\tAT\tA\t.\t.\t.",
    ]);
    let text = [
        "##fileformat=VCFv4.2",
        CSQ_HEADER,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "1\t100\trs1\tG\tGTAT,C\t.\tPASS\tAC=1;CSQ=ignored;DP=3;CSQ=TAT|x_variant|Transcript|T1,\
         C|y_variant|Transcript|T1,TAT|x_variant|Transcript|T1",
        "1\t200\t.\tAT\tA\t.\t.\tCSQ=-|d_variant|Transcript|T2",
        "1\t300\t.\tA\tT\t.\t.\t.",
        "",
    ]
    .join("\n");
    let entries = parse_vcf(&text, &inputs).unwrap();
    assert_eq!(
        keys(&entries),
        [
            "rec0|record",
            "rec0|TAT|T1|Transcript",
            "rec0|C|T1|Transcript",
            "rec0|TAT|T1|Transcript#2",
            "rec1|record",
            "rec1|-|T2|Transcript",
            "line2|record",
        ]
    );
    assert_eq!(
        fields(&entries[0]),
        [
            ("CHROM", "1"),
            ("POS", "100"),
            ("ID", "rs1"),
            ("REF", "G"),
            ("ALT", "GTAT,C"),
            ("QUAL", "."),
            ("FILTER", "PASS"),
            ("INFO_without_CSQ", "AC=1;DP=3"),
            ("CSQ_present", "yes"),
        ]
    );
    assert_eq!(
        entries[0],
        Entry {
            key: "rec0|record".to_string(),
            record: Some(0),
            fields: entries[0].fields.clone(),
            ..Entry::default()
        }
    );
    assert_eq!(
        fields(&entries[1]),
        [
            ("Allele", "TAT"),
            ("Consequence", "x_variant"),
            ("Feature_type", "Transcript"),
            ("Feature", "T1"),
        ]
    );
    // The first base of the two alternates is not dropped, since one is as long as the
    // reference: the insertion's allele is spelt by neither.
    let shapes: Vec<_> = entries
        .iter()
        .map(|entry| {
            (
                entry.allele.as_str(),
                entry.consequence.as_str(),
                entry.indel,
                entry.insertion,
                entry.reference.as_deref(),
                entry.line,
                entry.allele_location.as_deref(),
            )
        })
        .collect();
    assert_eq!(
        shapes,
        [
            ("", "", None, None, None, None, None),
            ("TAT", "x_variant", None, Some(false), None, Some(0), None),
            (
                "C",
                "y_variant",
                Some(false),
                Some(false),
                Some("G"),
                Some(0),
                Some("1:100")
            ),
            ("TAT", "x_variant", None, Some(false), None, Some(0), None),
            ("", "", None, None, None, None, None),
            (
                "-",
                "d_variant",
                Some(true),
                Some(false),
                Some("T"),
                Some(1),
                Some("1:201")
            ),
            ("", "", None, None, None, None, None),
        ]
    );
    assert_eq!(entries[4].fields["INFO_without_CSQ"], "");
    assert_eq!(entries[6].fields["INFO_without_CSQ"], ".");
    assert_eq!(entries[6].fields["CSQ_present"], "no");
}

#[test]
fn a_field_the_format_list_names_twice_keeps_its_first_place_and_its_last_value() {
    let text =
        "##INFO=<ID=CSQ,Description=\"Format: Allele|Feature|Allele\">\">\n1\tx\t.\tA\tT\t.\t.\tCSQ=T|f1|G\n";
    let entries = parse_vcf(text, &[]).unwrap();
    assert_eq!(keys(&entries), ["line0|record", "line0|G|f1|"]);
    assert_eq!(fields(&entries[1]), [("Allele", "G"), ("Feature", "f1")]);
    assert_eq!(entries[1].allele, "G");
    assert_eq!(entries[1].reference, None);
    assert_eq!(entries[1].allele_location, None);
}

#[test]
fn a_line_that_names_no_input_record_carries_the_ensembl_format_line_it_is_written_from() {
    let inputs = strings(&["1\t101\t100\t-/TAT\t+\tr1", "1\t200\t201\tAT/A\t+"]);
    let text = format!(
        "{CSQ_HEADER}\n1\t100\tr1\tG\tGTAT\t.\t.\tCSQ=TAT|m|Transcript|T9\n\
         1\t199\t.\tGAT\tG\t.\t.\tCSQ=AT|m|Transcript|T9\n\
         1\t300\t.\tA\tT\t.\t.\tCSQ=T|m|Transcript|T9\n"
    );
    let entries = parse_vcf(&text, &inputs).unwrap();
    let read: Vec<_> = entries
        .iter()
        .map(|entry| {
            (
                entry.key.as_str(),
                entry.line,
                entry.reference.as_deref(),
                entry.allele_location.as_deref(),
                entry.indel,
                entry.insertion,
            )
        })
        .collect();
    assert_eq!(
        read,
        [
            ("line0|record", None, None, None, None, None),
            (
                "line0|TAT|T9|Transcript",
                Some(0),
                Some(""),
                Some("1:100-101"),
                Some(true),
                Some(true)
            ),
            ("line1|record", None, None, None, None, None),
            (
                "line1|AT|T9|Transcript",
                Some(1),
                None,
                None,
                Some(true),
                Some(false)
            ),
            ("line2|record", None, None, None, None, None),
            (
                "line2|T|T9|Transcript",
                None,
                Some("A"),
                Some("1:300"),
                Some(false),
                Some(false)
            ),
        ]
    );
    assert!(entries.iter().all(|entry| entry.record.is_none()));
}

#[test]
fn the_header_line_names_the_columns_of_a_default_or_tab_output() {
    assert_eq!(
        header_columns("## x\n#A\tB\n#C\n1\t2\n"),
        strings(&["A", "B"])
    );
    assert_eq!(header_columns("## x\r\n#A\tB\r\n"), strings(&["A", "B"]));
    assert_eq!(header_columns("## x\n"), Vec::<String>::new());
    assert_eq!(header_columns("#"), strings(&[""]));
}
