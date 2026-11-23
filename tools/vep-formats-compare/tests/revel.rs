// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The sibling-allele classes and the REVEL table classes, read through the comparison of
//! two entry lists on planted rows and on made-up tables with their tabix index, and the
//! binary on a made-up table under a CSI index alone.

mod compare_common;
mod support;

use compare_common::*;
use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::Path;
use std::path::PathBuf;
use vep_formats_compare::classes::tables;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Error;
use vep_formats_compare::error::Uncaught;
use vep_formats_compare::input::read_text;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::field_class_accepts;
use vep_formats_compare::rules::revel::revel_record;
use vep_formats_compare::rules::revel::trimmed_allele_keys;
use vep_formats_compare::rules::revel::AlleleKeys;
use vep_formats_compare::rules::revel::RevelRow;
use vep_formats_compare::rules::revel::RevelTable;
use vep_formats_compare::rules::RuleContext;
use vep_formats_compare::rules::Side;
use vep_formats_compare::shared::ABSENT;

const REVEL_TABLE_HEADER: &str =
    "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid";
const SEVEN_COLUMN_HEADER: &str = "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL";
const REVEL_OTHER: &str = "revel_other_reference_amino_acid_rows";
const REVEL_SIBLING: &str = "revel_sibling_allele_rows";
const ALPHAMISSENSE_SIBLING: &str = "alphamissense_sibling_allele_rows";
const TRANSCRIPT: &str = "ENST00000299443";

fn coding_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    amino_acids: &str,
    extra: &str,
) -> String {
    format!(
        "{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\tmissense_variant\t-\t-\t-\t{amino_acids}\t-\t-\t{extra}"
    )
}

/// A row of `am_rows`: its allele, feature, Amino_acids, am_class and am_pathogenicity.
type AmRow<'a> = (&'a str, &'a str, &'a str, &'a str, &'a str);

/// Default-format rows of the record 21 13610388 A C,T, its alleles spelt by the
/// Uploaded_variation `name`; an empty class leaves the two fields out.
fn am_rows(rows: &[AmRow], name: &str) -> String {
    let rows: Vec<String> = rows
        .iter()
        .map(|(allele, feature, amino_acids, am_class, score)| {
            let extra = if am_class.is_empty() {
                "IMPACT=MODERATE;STRAND=1".to_string()
            } else {
                format!("IMPACT=MODERATE;STRAND=1;am_class={am_class};am_pathogenicity={score}")
            };
            coding_row(name, "21:13610388", allele, feature, amino_acids, &extra)
        })
        .collect();
    default_text(&rows)
}

const TWO_ALTERNATES: &str = "21_13610388_A/C/T";
const ONE_ALTERNATE: &str = "21_13610388_A/T";

fn revel_extra(score: &str) -> String {
    if score.is_empty() {
        "IMPACT=MODERATE;STRAND=1".to_string()
    } else {
        format!("IMPACT=MODERATE;STRAND=1;REVEL={score}")
    }
}

/// Default-format rows of the record 21 13610388 A C,T on the transcript, its alleles
/// spelt by the Uploaded_variation `name`: (allele, Amino_acids, REVEL) each, an empty
/// score leaving the field out.
fn revel_rows(rows: &[(&str, &str, &str)], name: &str) -> String {
    let rows: Vec<String> = rows
        .iter()
        .map(|(allele, amino_acids, score)| {
            coding_row(
                name,
                "21:13610388",
                allele,
                TRANSCRIPT,
                amino_acids,
                &revel_extra(score),
            )
        })
        .collect();
    default_text(&rows)
}

/// One BGZF block of `data`: a gzip member whose `BC` extra field holds its size.
fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut deflater = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::new(6));
    deflater.write_all(data).expect("the block deflates");
    let body = deflater.finish().expect("the block deflates");
    let mut crc = flate2::Crc::new();
    crc.update(data);
    let mut block = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00".to_vec();
    block.extend_from_slice(&(body.len() as u16 + 25).to_le_bytes());
    block.extend_from_slice(&body);
    block.extend_from_slice(&crc.sum().to_le_bytes());
    block.extend_from_slice(&(data.len() as u32).to_le_bytes());
    block
}

fn gzipped(bytes: &[u8]) -> Vec<u8> {
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gzip.write_all(bytes).expect("the index deflates");
    gzip.finish().expect("the index deflates")
}

fn index_path(table: &Path) -> PathBuf {
    PathBuf::from(format!("{}.tbi", table.display()))
}

/// A made-up REVEL table of `lines` under `header`, each line a BGZF block of its own, with
/// a tabix index on the 1-based `position_column` that files every row of a chromosome in
/// its first bin.
fn write_revel_table(
    directory: &Path,
    lines: &[String],
    position_column: i32,
    header: &str,
    name: &str,
) -> PathBuf {
    let mut blocks = vec![bgzf_block(format!("{header}\n").as_bytes())];
    let mut spans: Vec<(String, u64, u64)> = Vec::new();
    let mut at = blocks[0].len() as u64;
    for line in lines {
        let block = bgzf_block(format!("{line}\n").as_bytes());
        let end = at + block.len() as u64;
        let chromosome = line.split('\t').next().unwrap_or("");
        match spans.iter_mut().find(|(held, _, _)| held == chromosome) {
            Some(span) => span.2 = end,
            None => spans.push((chromosome.to_string(), at, end)),
        }
        at = end;
        blocks.push(block);
    }
    blocks.push(bgzf_block(b""));
    let path = directory.join(name);
    std::fs::write(&path, blocks.concat()).expect("the table is written");
    let chunks: Vec<(String, u64, u64)> = spans
        .into_iter()
        .map(|(chromosome, first, end)| (chromosome, first << 16, end << 16))
        .collect();
    write_tabix_index(&path, &chunks, position_column);
    path
}

/// Writes the tabix index of `table` over `chunks`, each a sequence name with the virtual
/// offsets its one chunk of the first bin runs between, a block's file offset above the low
/// 16 bits and a place in its text in them, on the 1-based `position_column`.
fn write_tabix_index(table: &Path, chunks: &[(String, u64, u64)], position_column: i32) {
    let names: Vec<u8> = chunks
        .iter()
        .flat_map(|(chromosome, _, _)| chromosome.bytes().chain([0]))
        .collect();
    let mut index = b"TBI\x01".to_vec();
    for number in [
        chunks.len() as i32,
        0,
        1,
        position_column,
        position_column,
        i32::from(b'#'),
        0,
        names.len() as i32,
    ] {
        index.extend_from_slice(&number.to_le_bytes());
    }
    index.extend_from_slice(&names);
    for (_, first, end) in chunks {
        index.extend_from_slice(&1i32.to_le_bytes());
        index.extend_from_slice(&0u32.to_le_bytes());
        index.extend_from_slice(&1i32.to_le_bytes());
        index.extend_from_slice(&first.to_le_bytes());
        index.extend_from_slice(&end.to_le_bytes());
        index.extend_from_slice(&0i32.to_le_bytes());
    }
    std::fs::write(index_path(table), gzipped(&index)).expect("the index is written");
}

fn nine_column_table(directory: &Path, lines: &[&str]) -> PathBuf {
    let lines: Vec<String> = lines.iter().map(|line| line.to_string()).collect();
    write_revel_table(directory, &lines, 3, REVEL_TABLE_HEADER, "revel.tsv.gz")
}

/// A row of a made-up table: its alternate allele, aaref, aaalt, score and transcript ids,
/// ids of `None` leaving the ninth cell out.
type TableRow<'a> = (&'a str, &'a str, &'a str, &'a str, Option<&'a str>);

/// Where the rows of a made-up nine-column table stand and how the run read their ids.
#[derive(Clone, Copy)]
struct Placed {
    no_match: bool,
    position: i64,
}

const AT_THE_RECORD: Placed = Placed {
    no_match: false,
    position: 13610388,
};

/// The side tables of a comparison that reads a made-up nine-column table whose rows stand
/// at 21 `placed.position` with the reference allele A, written under `directory`.
fn revel_side_placed(directory: &Path, rows: &[TableRow], placed: Placed) -> Side {
    let lines: Vec<String> = rows
        .iter()
        .map(|(alternate, aaref, aaalt, score, ids)| {
            let position = placed.position.to_string();
            let mut cells = vec!["21", "999", &position, "A", alternate, aaref, aaalt, score];
            cells.extend(ids);
            cells.join("\t")
        })
        .collect();
    std::fs::create_dir_all(directory).expect("a directory for the table");
    let path = write_revel_table(directory, &lines, 3, REVEL_TABLE_HEADER, "revel.tsv.gz");
    side_of(RevelTable::new(&path, placed.no_match).expect("the table reads"))
}

fn revel_side(directory: &Path, rows: &[TableRow]) -> Side {
    revel_side_placed(directory, rows, AT_THE_RECORD)
}

fn side_of(table: RevelTable) -> Side {
    Side {
        revel_table: Some(table),
        ..Side::default()
    }
}

fn applied(documented: &Documented) -> BTreeMap<String, i64> {
    counts(&documented.applied)
}

/// The classes a column's values were set aside under, each with its count.
fn masked(report: &Value, column: &str) -> BTreeMap<String, i64> {
    at(report, &["columns", column, "masked"])
        .as_dict()
        .expect("the masked counts")
        .iter()
        .map(|(klass, counted)| match counted {
            Value::Int(counted) => (klass.clone(), *counted),
            other => panic!("a count is {other:?}"),
        })
        .collect()
}

fn lines_of(lines: &[&str]) -> Vec<String> {
    lines.iter().map(|line| line.to_string()).collect()
}

fn refs(entries: &[Entry]) -> Vec<&Entry> {
    entries.iter().collect()
}

/// The report of two entry lists under the built-in `classes`, and what was applied.
fn compare_under(
    expected: &[Entry],
    actual: &[Entry],
    side: &Side,
    classes: &[&str],
    lines: Option<&[String]>,
) -> (Value, Documented) {
    let mut documented = documented_with(classes);
    let options = CompareOptions {
        input_lines: lines,
        ..CompareOptions::new(side)
    };
    let report =
        compare_entries(expected, actual, &mut documented, &options).expect("the entries compare");
    (report, documented)
}

/// The verdict of two default-format texts read beside the table of `side` and the input
/// lines `lines`, and what was applied.
fn compare_revel_lines(
    reference: &str,
    candidate: &str,
    side: &Side,
    classes: &[&str],
    lines: Option<&[&str]>,
) -> (String, BTreeMap<String, i64>) {
    let lines = lines.map(lines_of);
    let (report, documented) = compare_under(
        &parse_default(reference).expect("the reference parses"),
        &parse_default(candidate).expect("the candidate parses"),
        side,
        classes,
        lines.as_deref(),
    );
    (verdict(&report).to_string(), applied(&documented))
}

fn compare_revel(
    reference: &str,
    candidate: &str,
    side: &Side,
    classes: &[&str],
) -> (String, BTreeMap<String, i64>) {
    compare_revel_lines(reference, candidate, side, classes, None)
}

fn green(klass: &str) -> (String, BTreeMap<String, i64>) {
    ("green".to_string(), count_map(&[(klass, 1)]))
}

fn red() -> (String, BTreeMap<String, i64>) {
    ("red".to_string(), count_map(&[]))
}

/// The JSON text of one record of the input line `input` under the identifier `id`: the
/// members `stated` of the record object and its transcript entries.
fn json_record(input: &str, id: &str, stated: &str, entries: &[String]) -> String {
    let input = input.replace('\t', "\\t");
    format!(
        "{{\"input\": \"{input}\", \"id\": \"{id}\", {stated}, \"most_severe_consequence\": \"missense_variant\", \"transcript_consequences\": [{}]}}\n",
        entries.join(", ")
    )
}

#[test]
fn the_alphamissense_sibling_class_accepts_the_prediction_of_the_sibling_allele() {
    // The record 21 13610388 A C,T on ENST00000299443, both alleles giving M/L: the
    // reference writes the C row's prediction on the T row, the candidate the table's row
    // for T or, where the table has none, nothing. The T row is set aside once in every
    // format.
    let klass = ALPHAMISSENSE_SIBLING;
    let one = count_map(&[(klass, 1)]);
    let c_row = ("C", TRANSCRIPT, "M/L", "likely_benign", "0.1079");
    let reference = am_rows(
        &[c_row, ("T", TRANSCRIPT, "M/L", "likely_benign", "0.1079")],
        TWO_ALTERNATES,
    );
    let candidate = am_rows(
        &[
            c_row,
            ("T", TRANSCRIPT, "M/L", "likely_pathogenic", "0.9000"),
        ],
        TWO_ALTERNATES,
    );
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    assert_eq!(masked(&report, "am_class"), one);
    assert_eq!(masked(&report, "am_pathogenicity"), one);
    assert!(texts(at(&report, &["mask_examples", klass]))[0]
        .starts_with("21:13610388|T|ENST00000299443|Transcript: "));
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    let without_row = am_rows(&[c_row, ("T", TRANSCRIPT, "M/L", "", "")], TWO_ALTERNATES);
    let (report, documented) = compare_default(&reference, &without_row, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    // Under an identifier a default-format row does not say whether its Location's rows
    // are one record's, so the class leaves it; the JSON entries below name their record.
    let named: Vec<String> = [("likely_benign", "0.1079"), ("likely_pathogenic", "0.9000")]
        .iter()
        .map(|(am_class, score)| {
            am_rows(
                &[c_row, ("T", TRANSCRIPT, "M/L", am_class, score)],
                "sibling_ct",
            )
        })
        .collect();
    assert_eq!(verdict_of(&named[0], &named[1], &[klass]), "red");
    // The JSON format carries the two fields as one alphamissense object.
    let inputs = lines_of(&["21\t13610388\tsibling_ct\tA\tC,T\t.\t.\t."]);
    let record = |t_class: &str, t_score: &str| {
        let entries: Vec<String> = [("C", "likely_benign", "0.1079"), ("T", t_class, t_score)]
            .iter()
            .map(|(allele, am_class, score)| {
                format!(
                    "{{\"variant_allele\": \"{allele}\", \"transcript_id\": \"{TRANSCRIPT}\", \"consequence_terms\": [\"missense_variant\"], \"amino_acids\": \"M/L\", \"alphamissense\": {{\"am_class\": \"{am_class}\", \"am_pathogenicity\": {score}}}}}"
                )
            })
            .collect();
        json_record(
            &inputs[0],
            "sibling_ct",
            "\"allele_string\": \"A/C/T\"",
            &entries,
        )
    };
    let (report, documented) = compare_under(
        &parse_json(&record("likely_benign", "0.1079"), &inputs).expect("the reference parses"),
        &parse_json(&record("likely_pathogenic", "0.9"), &inputs).expect("the candidate parses"),
        &Side::default(),
        &[klass],
        None,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
}

#[test]
fn the_alphamissense_sibling_class_refuses_a_value_no_sibling_allele_explains() {
    // A differing am_ field stays a difference on a one-allele record, beside a sibling
    // giving another amino acid or sitting on another feature, where the reference's value
    // is not the sibling row's, and where the candidate alone fills the field.
    let klass = ALPHAMISSENSE_SIBLING;
    let c_row = ("C", TRANSCRIPT, "M/L", "likely_benign", "0.1079");
    let t_reference = ("T", TRANSCRIPT, "M/L", "likely_benign", "0.1079");
    let t_candidate = ("T", TRANSCRIPT, "M/L", "likely_pathogenic", "0.9000");
    assert_eq!(
        verdict_of(
            &am_rows(&[t_reference], ONE_ALTERNATE),
            &am_rows(&[t_candidate], ONE_ALTERNATE),
            &[klass]
        ),
        "red"
    );
    for other in [
        ("C", TRANSCRIPT, "M/V", "likely_benign", "0.1079"),
        ("C", "ENST00000620442", "M/L", "likely_benign", "0.1079"),
        ("C", TRANSCRIPT, "M/L", "ambiguous", "0.4500"),
    ] {
        assert_eq!(
            verdict_of(
                &am_rows(&[other, t_reference], TWO_ALTERNATES),
                &am_rows(&[other, t_candidate], TWO_ALTERNATES),
                &[klass]
            ),
            "red",
            "{other:?}"
        );
    }
    let t_empty = ("T", TRANSCRIPT, "M/L", "", "");
    let (report, documented) = compare_default(
        &am_rows(&[c_row, t_empty], TWO_ALTERNATES),
        &am_rows(&[c_row, t_candidate], TWO_ALTERNATES),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), count_map(&[]));
    // A row that shares its sibling's class and carries a score of its own is no copy of
    // the sibling's table row: the two fields are read together, as the JSON object is.
    let own_score = ("T", TRANSCRIPT, "M/L", "likely_benign", "0.5000");
    let (report, documented) = compare_default(
        &am_rows(&[c_row, own_score], TWO_ALTERNATES),
        &am_rows(
            &[
                c_row,
                ("T", TRANSCRIPT, "M/L", "likely_pathogenic", "0.5000"),
            ],
            TWO_ALTERNATES,
        ),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), count_map(&[]));
}

#[test]
fn the_revel_sibling_class_accepts_the_score_of_the_sibling_allele() {
    // The record 21 13610388 A C,T, both alleles giving M/L: the reference writes the C
    // row's score on the T row, the candidate the table's row for T or nothing; the tab and
    // JSON formats carry the field as a column and as the revel member.
    let klass = REVEL_SIBLING;
    let one = count_map(&[(klass, 1)]);
    let leucine = ("C", "M/L", "0.025");
    let reference = revel_rows(&[leucine, ("T", "M/L", "0.025")], TWO_ALTERNATES);
    let candidate = revel_rows(&[leucine, ("T", "M/L", "0.027")], TWO_ALTERNATES);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    assert_eq!(masked(&report, "REVEL"), one);
    assert_eq!(verdict_of(&reference, &candidate, &[]), "red");
    assert_eq!(
        verdict_of(
            &reference,
            &revel_rows(&[leucine, ("T", "M/L", "")], TWO_ALTERNATES),
            &[klass]
        ),
        "green"
    );
    // Under an identifier the default format's rows may be two records': the class leaves
    // them.
    let named: Vec<String> = ["0.025", "0.027"]
        .iter()
        .map(|score| revel_rows(&[leucine, ("T", "M/L", score)], "two_leucines"))
        .collect();
    assert_eq!(verdict_of(&named[0], &named[1], &[klass]), "red");
    let tab_header = "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tAmino_acids\tREVEL\n";
    let tab_rows: String = ["C", "T"]
        .iter()
        .map(|allele| {
            format!("{TWO_ALTERNATES}\t21:13610388\t{allele}\t{TRANSCRIPT}\tTranscript\tmissense_variant\tM/L\t0.025\n")
        })
        .collect();
    let tab = format!("{tab_header}{tab_rows}");
    let t_row = format!("\tT\t{TRANSCRIPT}\tTranscript\tmissense_variant\tM/L\t");
    let (report, documented) = compare_under(
        &parse_tab(&tab).expect("the reference parses"),
        &parse_tab(&tab.replace(&format!("{t_row}0.025"), &format!("{t_row}0.027")))
            .expect("the candidate parses"),
        &Side::default(),
        &[klass],
        None,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    let inputs = lines_of(&["21\t13610388\ttwo_leucines\tA\tC,T\t.\t.\t."]);
    let record = |t_score: &str| {
        let entries: Vec<String> = [("C", "0.025"), ("T", t_score)]
            .iter()
            .map(|(allele, score)| {
                format!(
                    "{{\"variant_allele\": \"{allele}\", \"transcript_id\": \"{TRANSCRIPT}\", \"consequence_terms\": [\"missense_variant\"], \"amino_acids\": \"M/L\", \"revel\": {score}}}"
                )
            })
            .collect();
        json_record(
            &inputs[0],
            "two_leucines",
            "\"allele_string\": \"A/C/T\"",
            &entries,
        )
    };
    let (report, documented) = compare_under(
        &parse_json(&record("0.025"), &inputs).expect("the reference parses"),
        &parse_json(&record("0.027"), &inputs).expect("the candidate parses"),
        &Side::default(),
        &[klass],
        None,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
}

#[test]
fn the_revel_sibling_class_refuses_a_score_no_sibling_allele_explains() {
    // A differing REVEL stays a difference on a one-allele record, beside a sibling giving
    // another amino acid, where the reference's score is not the sibling row's, and where
    // the candidate alone fills the field.
    let klass = REVEL_SIBLING;
    assert_eq!(
        verdict_of(
            &revel_rows(&[("T", "M/L", "0.025")], ONE_ALTERNATE),
            &revel_rows(&[("T", "M/L", "0.027")], ONE_ALTERNATE),
            &[klass]
        ),
        "red"
    );
    let (valine, other_pair) = (("G", "M/V", "0.025"), "21_13610388_A/G/T");
    assert_eq!(
        verdict_of(
            &revel_rows(&[valine, ("T", "M/L", "0.025")], other_pair),
            &revel_rows(&[valine, ("T", "M/L", "0.027")], other_pair),
            &[klass]
        ),
        "red"
    );
    let other_score = ("C", "M/L", "0.043");
    assert_eq!(
        verdict_of(
            &revel_rows(&[other_score, ("T", "M/L", "0.025")], TWO_ALTERNATES),
            &revel_rows(&[other_score, ("T", "M/L", "0.027")], TWO_ALTERNATES),
            &[klass]
        ),
        "red"
    );
    let leucine = ("C", "M/L", "0.025");
    let (report, documented) = compare_default(
        &revel_rows(&[leucine, ("T", "M/L", "")], TWO_ALTERNATES),
        &revel_rows(&[leucine, ("T", "M/L", "0.027")], TWO_ALTERNATES),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), count_map(&[]));
}

fn row(cells: [&str; 7]) -> RevelRow {
    let cell = |text: &str| Some(text.to_string());
    RevelRow {
        position: cells[0].to_string(),
        reference: cell(cells[1]),
        alternate: cell(cells[2]),
        reference_amino_acid: cells[3].to_string(),
        alternate_amino_acid: cells[4].to_string(),
        score: cell(cells[5]),
        transcript_ids: cell(cells[6]),
    }
}

fn rows_at(table: &RevelTable, chromosome: &str, start: i64, end: i64) -> Vec<RevelRow> {
    table
        .rows(chromosome, start, end)
        .expect("the rows are read")
}

#[test]
fn the_table_reads_the_rows_of_a_span_in_file_order_through_the_index() {
    // The committed table slice of each REVEL corpus, written by bgzip and tabix: at every
    // position, and over spans, the reader returns the rows a sequential read of the file
    // holds there, in file order, on the position column the index names; a name with
    // `chr`, a chromosome the table lacks and an empty span are read as the plugin's query
    // reads them.
    for (corpus, column) in [("GRCh37-plugin-revel", 1), ("GRCh38-plugin-revel", 2)] {
        let path = golden().join(corpus).join("plugin_data/revel.tsv.gz");
        let table = RevelTable::new(&path, false).expect("the table reads");
        let mut held: BTreeMap<i64, Vec<RevelRow>> = BTreeMap::new();
        for line in data_lines(&read_text(&path).expect("the table inflates")) {
            let cells: Vec<&str> = line.split('\t').collect();
            let position: i64 = cells[column].parse().expect("a position");
            held.entry(position).or_default().push(row([
                cells[column],
                cells[3],
                cells[4],
                cells[5],
                cells[6],
                cells[7],
                cells[8],
            ]));
        }
        assert!(held.len() > 100, "{corpus}");
        let expected = |start: i64, end: i64| -> Vec<RevelRow> {
            held.range(start..=end)
                .flat_map(|(_, rows)| rows.iter().cloned())
                .collect()
        };
        let positions: Vec<i64> = held.keys().copied().collect();
        for at in &positions {
            assert_eq!(
                rows_at(&table, "21", *at, *at),
                expected(*at, *at),
                "{corpus} {at}"
            );
        }
        let last = positions[positions.len() - 1];
        for (start, end) in [
            (positions[0], last),
            (positions[3] - 40, positions[40] + 7),
            (last + 1, last + 9),
        ] {
            assert_eq!(
                rows_at(&table, "21", start, end),
                expected(start, end),
                "{corpus} {start} {end}"
            );
        }
        let fifth = positions[5];
        assert_eq!(
            rows_at(&table, "chr21", fifth, fifth),
            rows_at(&table, "21", fifth, fifth)
        );
        assert!(rows_at(&table, "22", fifth, fifth).is_empty());
        assert!(rows_at(&table, "21", fifth, fifth - 1).is_empty());
    }
}

/// Writes over the index of `table`, which `write_revel_table` wrote of `lines` under
/// `header`, a tabix index on the 1-based `position_column` that files the rows as `tabix`
/// files rows of one position each, in file order: every row in a chunk of the 16 kb bin of
/// its position, and per sequence a linear index that holds, for each 16 kb window up to
/// its last row's, the virtual offset of the first row at or after the window's first
/// position.
fn write_windowed_tabix_index(table: &Path, header: &str, lines: &[String], position_column: i32) {
    // Per sequence name, the chunk of each bin and the linear index.
    type Sequence = (String, Vec<(u32, u64, u64)>, Vec<u64>);
    let mut sequences: Vec<Sequence> = Vec::new();
    let mut at = bgzf_block(format!("{header}\n").as_bytes()).len() as u64;
    for line in lines {
        let end = at + bgzf_block(format!("{line}\n").as_bytes()).len() as u64;
        let cells: Vec<&str> = line.split('\t').collect();
        let position: u64 = cells[position_column as usize - 1]
            .parse()
            .expect("a position");
        let window = ((position - 1) >> 14) as usize;
        if sequences.last().is_none_or(|held| held.0 != cells[0]) {
            sequences.push((cells[0].to_string(), Vec::new(), Vec::new()));
        }
        let (_, bins, linear) = sequences.last_mut().expect("a sequence");
        let bin = 4681 + window as u32;
        match bins.last_mut() {
            Some(held) if held.0 == bin => held.2 = end << 16,
            _ => bins.push((bin, at << 16, end << 16)),
        }
        // A window without a row holds the offset of the first row after it.
        if linear.len() <= window {
            linear.resize(window + 1, at << 16);
        }
        at = end;
    }
    let names: Vec<u8> = sequences
        .iter()
        .flat_map(|(name, _, _)| name.bytes().chain([0]))
        .collect();
    let mut index = b"TBI\x01".to_vec();
    for number in [
        sequences.len() as i32,
        0,
        1,
        position_column,
        position_column,
        i32::from(b'#'),
        0,
        names.len() as i32,
    ] {
        index.extend_from_slice(&number.to_le_bytes());
    }
    index.extend_from_slice(&names);
    for (_, bins, linear) in &sequences {
        index.extend_from_slice(&(bins.len() as i32).to_le_bytes());
        for (bin, first, end) in bins {
            index.extend_from_slice(&bin.to_le_bytes());
            index.extend_from_slice(&1i32.to_le_bytes());
            index.extend_from_slice(&first.to_le_bytes());
            index.extend_from_slice(&end.to_le_bytes());
        }
        index.extend_from_slice(&(linear.len() as i32).to_le_bytes());
        for offset in linear {
            index.extend_from_slice(&offset.to_le_bytes());
        }
    }
    std::fs::write(index_path(table), gzipped(&index)).expect("the index is written");
}

#[test]
fn a_span_that_starts_at_or_below_zero_is_read_from_the_first_position() {
    // A made-up table under an index that files its rows as `tabix` files them: the rows of
    // 21 stand in three 16 kb windows, so its linear index holds three offsets, and 22 holds
    // one row. A span that starts at zero or below it is read from the sequence's first
    // position, as `tabix` reads it: it gives the rows the span from position 1 gives,
    // however far below zero it starts, and no row where it ends below position 1.
    let directory = scratch("revel-span-from-zero-or-below");
    let lines = lines_of(&[
        "21\t900\t100\tA\tG\tQ\tH\t0.1\tENST1",
        "21\t901\t5000\tA\tG\tQ\tH\t0.2\tENST1",
        "21\t902\t20000\tA\tG\tQ\tH\t0.3\tENST1",
        "21\t903\t40000\tA\tG\tQ\tH\t0.4\tENST1",
        "22\t904\t50\tA\tG\tQ\tH\t0.5\tENST2",
    ]);
    let name = "revel.synthetic.tsv.gz";
    let path = write_revel_table(&directory, &lines, 3, REVEL_TABLE_HEADER, name);
    write_windowed_tabix_index(&path, REVEL_TABLE_HEADER, &lines, 3);
    let table = RevelTable::new(&path, false).expect("the table reads");
    assert_eq!(table._linear["21"].len(), 3);
    assert_eq!(table._linear["22"].len(), 1);
    let positions = |chromosome: &str, start: i64, end: i64| -> Vec<String> {
        rows_at(&table, chromosome, start, end)
            .into_iter()
            .map(|row| row.position)
            .collect()
    };
    // From position 1 or above the linear index gives the offset of the span's own window.
    let from_the_first_position_on: [(i64, i64, &[&str]); 5] = [
        (1, 100, &["100"]),
        (1, 5000, &["100", "5000"]),
        (101, 4999, &[]),
        (16_385, 40_000, &["20000", "40000"]),
        (40_000, 40_000, &["40000"]),
    ];
    for (start, end, expected) in from_the_first_position_on {
        assert_eq!(positions("21", start, end), expected, "21 {start} to {end}");
    }
    let from_zero_or_below: [(&str, i64, i64, &[&str]); 10] = [
        ("21", 0, 100, &["100"]),
        ("21", 0, 5000, &["100", "5000"]),
        ("21", -5, 5000, &["100", "5000"]),
        ("21", -20_000, 100, &["100"]),
        (
            "21",
            -1_000_000_000_000,
            40_000,
            &["100", "5000", "20000", "40000"],
        ),
        ("21", i64::MIN, 100, &["100"]),
        ("22", 0, 50, &["50"]),
        ("22", -1_000_000_000_000, 50, &["50"]),
        ("21", 0, 0, &[]),
        ("21", -9, -3, &[]),
    ];
    for (chromosome, start, end, expected) in from_zero_or_below {
        assert_eq!(
            positions(chromosome, start, end),
            expected,
            "{chromosome} {start} to {end}"
        );
    }
}

fn refusal(path: &Path) -> String {
    match RevelTable::new(path, false).expect_err("the file is refused") {
        Error::Comparator(text) => text,
    }
}

#[test]
fn the_table_reads_each_layout_and_refuses_what_is_no_table() {
    // The header's column count fixes the layout (REVEL.pm:127-130, :209-241): seven
    // columns have no transcript ids and one position column, eight and nine have the ids
    // in the ninth cell of a row that carries one; a cell at or past a row's last filled
    // one is undefined. A file without its index, and one whose header is no REVEL header,
    // is refused by name.
    let directory = scratch("revel-layouts");
    let nine = nine_column_table(
        &directory,
        &[
            "21\t100\t500\tA\tG\tQ\tH\t0.742\tENST1;ENST2",
            "21\t100\t500\tA\tT\tQ\tL\t0.311",
            "21\t100\t500\tA",
            "21\t101\t501\tA\tG\tQ\tH\t\tENST1",
        ],
    );
    let rows = rows_at(
        &RevelTable::new(&nine, false).expect("the table reads"),
        "21",
        500,
        501,
    );
    assert_eq!(
        rows[0],
        row(["500", "A", "G", "Q", "H", "0.742", "ENST1;ENST2"])
    );
    assert_eq!(
        rows[1],
        RevelRow {
            transcript_ids: None,
            ..row(["500", "A", "T", "Q", "L", "0.311", ""])
        }
    );
    assert_eq!(
        rows[2],
        RevelRow {
            alternate: None,
            score: None,
            transcript_ids: None,
            ..row(["500", "A", "", "", "", "", ""])
        }
    );
    assert_eq!(rows[3], row(["501", "A", "G", "Q", "H", "", "ENST1"]));
    assert_eq!(rows.len(), 4);
    let seven = write_revel_table(
        &directory,
        &lines_of(&["21\t100\tA\tG\tQ\tH\t0.742"]),
        2,
        SEVEN_COLUMN_HEADER,
        "seven.tsv.gz",
    );
    assert_eq!(
        rows_at(
            &RevelTable::new(&seven, false).expect("the table reads"),
            "21",
            100,
            100
        ),
        vec![RevelRow {
            transcript_ids: None,
            ..row(["100", "A", "G", "Q", "H", "0.742", ""])
        }]
    );
    let grch37 = write_revel_table(
        &directory,
        &lines_of(&["21\t100\t500\tA\tG\tQ\tH\t0.742\tENST1"]),
        2,
        REVEL_TABLE_HEADER,
        "grch37.tsv.gz",
    );
    let grch37 = RevelTable::new(&grch37, false).expect("the table reads");
    let positions: Vec<String> = rows_at(&grch37, "21", 100, 100)
        .into_iter()
        .map(|row| row.position)
        .collect();
    assert_eq!(positions, ["100"]);
    assert!(rows_at(&grch37, "21", 500, 500).is_empty());
    let six = write_revel_table(
        &directory,
        &lines_of(&["21\t100\t500\tA\tG\t0.742"]),
        3,
        "#chr\thg19_pos\tgrch38_pos\tref\talt\tREVEL",
        "six.tsv.gz",
    );
    assert_eq!(
        refusal(&six),
        format!(
            "--revel-table: {} is not a REVEL table (a header of 7, 8 or 9 columns)",
            six.display()
        )
    );
    std::fs::remove_file(index_path(&nine)).expect("the index is removed");
    assert_eq!(
        refusal(&nine),
        format!(
            "--revel-table: cannot read the tabix index {}: No such file or directory",
            index_path(&nine).display()
        )
    );
}

fn keys(items: &[(&str, &str, i64)]) -> AlleleKeys {
    items
        .iter()
        .map(|(reference, alternate, start)| (reference.to_string(), alternate.to_string(), *start))
        .collect()
}

#[test]
fn the_trimmed_allele_keys_follow_the_trimming_orders_of_the_allele_match() {
    // Utils/Sequence.pm:965-1038 and :1190-1209: a one-base pair is its own key; a pair
    // that shares its first base starts one base on; shared last bases leave the start; a
    // pair emptied on one side is written `-`; both trimming orders are kept where they
    // differ.
    assert_eq!(trimmed_allele_keys("A", "G", 100), keys(&[("A", "G", 100)]));
    assert_eq!(
        trimmed_allele_keys("CA", "CG", 100),
        keys(&[("A", "G", 101)])
    );
    assert_eq!(
        trimmed_allele_keys("AC", "GC", 100),
        keys(&[("A", "G", 100)])
    );
    assert_eq!(
        trimmed_allele_keys("AC", "A", 100),
        keys(&[("C", "-", 101)])
    );
    assert_eq!(trimmed_allele_keys("C", "-", 101), keys(&[("C", "-", 101)]));
    assert_eq!(
        trimmed_allele_keys("CAC", "C", 100),
        keys(&[("AC", "-", 101), ("CA", "-", 100)])
    );
}

#[test]
fn a_deletions_start_past_the_widest_whole_number_is_the_greatest_position_held() {
    // The record 21 POS CA C starts one base past its POS, on the base it deletes: a POS
    // of more digits than any whole number held leaves the start, and the end of the one
    // deleted base, at the greatest position.
    let record = entry(
        "3|record",
        Some(3),
        None,
        "",
        "",
        "",
        &[
            ("CHROM", "21"),
            ("POS", &"9".repeat(40)),
            ("REF", "CA"),
            ("ALT", "C"),
        ],
    );
    let row = entry(
        "3|-|ENST1|Transcript",
        Some(3),
        None,
        "-",
        "ENST1",
        "frameshift_variant",
        &[],
    );
    let (chromosome, start, end, of_each_alternate, place) =
        revel_record(&row, Some(&[&record, &row]))
            .expect("the record is read")
            .expect("the record states its alleles");
    assert_eq!(
        (chromosome.as_str(), start, end, place),
        ("21", i64::MAX, i64::MAX, 0)
    );
    assert_eq!(of_each_alternate, [keys(&[("A", "-", i64::MAX)])]);
}

/// Whether `klass` accepts the field `name` of the reference's entry `e` beside the
/// candidate's `a`, the table of `side` and nothing else named.
fn accepts(
    klass: &str,
    e: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
    a: Option<&Entry>,
    side: &Side,
) -> bool {
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side,
    };
    field_class_accepts(klass, e, name, expected_value, actual_value, a, &context)
        .expect("the class is asked")
}

#[test]
fn the_other_reference_amino_acid_class_accepts_the_score_of_a_row_of_another_reference_amino_acid()
{
    // The record 21 13610388 A T, M/L on the transcript, against made-up tables. The
    // reference writes the score of the table row its walk reaches, which lists the
    // transcript and names the row's alternate amino acid beside another reference amino
    // acid; the candidate writes nothing. The class sets the value aside with the table and
    // nothing without it; it leaves a value the row of the row's own pair gives, a value no
    // row at the position gives, and a candidate that writes anything but what its own walk
    // gives.
    let directory = scratch("revel-other-reference-amino-acid");
    let ids = Some(TRANSCRIPT);
    let name = ONE_ALTERNATE;
    let reference = revel_rows(&[("T", "M/L", "0.300")], name);
    let empty = revel_rows(&[("T", "M/L", "")], name);
    let no_table = Side::default();
    let other = revel_side(&directory.join("other"), &[("T", "V", "L", "0.300", ids)]);
    assert_eq!(
        compare_revel(&reference, &empty, &other, &[REVEL_OTHER]),
        green(REVEL_OTHER)
    );
    assert_eq!(
        compare_revel(&reference, &empty, &no_table, &[REVEL_OTHER]),
        red()
    );
    assert_eq!(compare_revel(&reference, &empty, &other, &[]), red());
    assert_eq!(
        compare_revel(
            &reference,
            &revel_rows(&[("T", "M/L", "0.999")], name),
            &other,
            &[REVEL_OTHER]
        ),
        red()
    );
    let own_pair = revel_side(&directory.join("own"), &[("T", "M", "L", "0.300", ids)]);
    assert_eq!(
        compare_revel(&reference, &empty, &own_pair, &[REVEL_OTHER]),
        red()
    );
    assert_eq!(
        compare_revel(
            &revel_rows(&[("T", "M/L", "0.555")], name),
            &empty,
            &other,
            &[REVEL_OTHER]
        ),
        red()
    );
    let elsewhere = revel_side_placed(
        &directory.join("elsewhere"),
        &[("T", "V", "L", "0.300", ids)],
        Placed {
            position: 13610399,
            ..AT_THE_RECORD
        },
    );
    assert_eq!(
        compare_revel(&reference, &empty, &elsewhere, &[REVEL_OTHER]),
        red()
    );
    let other_allele = revel_side(&directory.join("allele"), &[("G", "V", "L", "0.300", ids)]);
    assert_eq!(
        compare_revel(&reference, &empty, &other_allele, &[REVEL_OTHER]),
        red()
    );
    let other_alternate_residue =
        revel_side(&directory.join("residue"), &[("T", "V", "F", "0.300", ids)]);
    assert_eq!(
        compare_revel(&reference, &empty, &other_alternate_residue, &[REVEL_OTHER]),
        red()
    );
    // The candidate's own walk goes on to the first row of the row's own pair: its score is
    // the candidate's value there, and an empty field is not.
    let later = revel_side(
        &directory.join("later"),
        &[("T", "V", "L", "0.300", ids), ("T", "M", "L", "0.400", ids)],
    );
    let own_score = revel_rows(&[("T", "M/L", "0.400")], name);
    assert_eq!(
        compare_revel(&reference, &own_score, &later, &[REVEL_OTHER]),
        green(REVEL_OTHER)
    );
    assert_eq!(
        compare_revel(&reference, &empty, &later, &[REVEL_OTHER]),
        red()
    );
    assert_eq!(
        compare_revel(&reference, &own_score, &other, &[REVEL_OTHER]),
        red()
    );
    // The two outputs spell the row's amino acids alike, as one pair.
    assert_eq!(
        compare_revel(
            &reference,
            &revel_rows(&[("T", "V/L", "")], name),
            &other,
            &[REVEL_OTHER]
        ),
        red()
    );
    // The cells are compared as written: a sign the row does not spell is another cell.
    let sign = revel_side(&directory.join("sign"), &[("T", "*", "L", "0.300", ids)]);
    assert_eq!(
        compare_revel(
            &revel_rows(&[("T", "U/L", "0.300")], name),
            &revel_rows(&[("T", "U/L", "")], name),
            &sign,
            &[REVEL_OTHER]
        ),
        green(REVEL_OTHER)
    );
    // Another plugin's field on the row is not the class's, and the JSON member is.
    let e = &parse_default(&reference).expect("the reference parses")[0];
    let a = &parse_default(&empty).expect("the candidate parses")[0];
    assert!(accepts(
        REVEL_OTHER,
        e,
        "REVEL",
        "0.300",
        ABSENT,
        Some(a),
        &other
    ));
    assert!(accepts(
        REVEL_OTHER,
        e,
        "revel",
        "0.300",
        ABSENT,
        Some(a),
        &other
    ));
    assert!(!accepts(
        REVEL_OTHER,
        e,
        "am_pathogenicity",
        "0.300",
        ABSENT,
        Some(a),
        &other
    ));
    assert!(!accepts(
        REVEL_OTHER,
        e,
        "REVEL",
        "0.300",
        ABSENT,
        None,
        &other
    ));
    assert!(!accepts(
        REVEL_OTHER,
        e,
        "REVEL",
        ABSENT,
        "0.300",
        Some(a),
        &other
    ));
    assert!(!accepts(
        REVEL_SIBLING,
        e,
        "REVEL",
        "0.300",
        ABSENT,
        Some(a),
        &other
    ));
}

#[test]
fn the_other_reference_amino_acid_class_reads_the_transcript_ids_as_the_run_did() {
    // REVEL.pm:191-197: a row whose ids cell is filled is the walk's only where it lists
    // the row's transcript, unless the run set `no_match`; a row with no ids is read on its
    // amino acid alone. Under `no_match` the reference's row is the first that names the
    // row's alternate amino acid, and the candidate's the first of the row's own pair.
    let directory = scratch("revel-transcript-ids");
    let unmatched = Placed {
        no_match: true,
        ..AT_THE_RECORD
    };
    let name = ONE_ALTERNATE;
    let reference = revel_rows(&[("T", "M/L", "0.300")], name);
    let empty = revel_rows(&[("T", "M/L", "")], name);
    let rows = [(
        "T",
        "V",
        "L",
        "0.300",
        Some("ENST00000620442;ENST00000620443"),
    )];
    assert_eq!(
        compare_revel(
            &reference,
            &empty,
            &revel_side(&directory.join("matched"), &rows),
            &[REVEL_OTHER]
        ),
        red()
    );
    assert_eq!(
        compare_revel(
            &reference,
            &empty,
            &revel_side_placed(&directory.join("unmatched"), &rows, unmatched),
            &[REVEL_OTHER]
        ),
        green(REVEL_OTHER)
    );
    let both = format!("ENST00000620442;{TRANSCRIPT}");
    let listed = [("T", "V", "L", "0.300", Some(both.as_str()))];
    assert_eq!(
        compare_revel(
            &reference,
            &empty,
            &revel_side(&directory.join("listed"), &listed),
            &[REVEL_OTHER]
        ),
        green(REVEL_OTHER)
    );
    let no_ids = [("T", "V", "L", "0.300", None)];
    assert_eq!(
        compare_revel(
            &reference,
            &empty,
            &revel_side(&directory.join("no-ids"), &no_ids),
            &[REVEL_OTHER]
        ),
        green(REVEL_OTHER)
    );
    let two = [
        ("T", "V", "L", "0.300", Some("ENST00000620442")),
        ("T", "M", "L", "0.400", Some("ENST00000620443")),
    ];
    let own_score = revel_rows(&[("T", "M/L", "0.400")], name);
    let two_unmatched = revel_side_placed(&directory.join("two-unmatched"), &two, unmatched);
    assert_eq!(
        compare_revel(&reference, &own_score, &two_unmatched, &[REVEL_OTHER]),
        green(REVEL_OTHER)
    );
    assert_eq!(
        compare_revel(&reference, &empty, &two_unmatched, &[REVEL_OTHER]),
        red()
    );
    assert_eq!(
        compare_revel(
            &reference,
            &own_score,
            &revel_side(&directory.join("two-matched"), &two),
            &[REVEL_OTHER]
        ),
        red()
    );
    // The seven-column layout carries no ids: every row is read on its amino acids.
    let seven = write_revel_table(
        &directory,
        &lines_of(&["21\t13610388\tA\tT\tV\tL\t0.300"]),
        2,
        SEVEN_COLUMN_HEADER,
        "seven.tsv.gz",
    );
    let seven = side_of(RevelTable::new(&seven, false).expect("the table reads"));
    assert_eq!(
        compare_revel(&reference, &empty, &seven, &[REVEL_OTHER]),
        green(REVEL_OTHER)
    );
}

/// The reference's and the candidate's row of a two-base substitution over the positions
/// 21:13610387-13610388 under the Uploaded_variation `name`, the reference's scored.
fn two_base_rows(name: &str) -> Vec<String> {
    ["0.029", ""]
        .iter()
        .map(|score| {
            default_text(&[coding_row(
                name,
                "21:13610387-13610388",
                "GT",
                TRANSCRIPT,
                "M/L",
                &revel_extra(score),
            )])
        })
        .collect()
}

#[test]
fn the_other_reference_amino_acid_class_reads_the_row_of_another_alternate_of_the_record() {
    // The record 21 13610388 A C,T: REVEL.pm:165 hands the match every alternate, so the
    // reference's value on the T row is the score of the first row that passes, a row of C
    // naming another reference amino acid; the candidate reads T's own rows, of which none
    // names the row's pair. The class takes the row where the record is stated, by an
    // Uploaded_variation that spells its alleles or by the one input line that carries its
    // identifier, and the sibling class, whose C row carries no value, does not.
    let directory = scratch("revel-another-alternate");
    let ids = Some(TRANSCRIPT);
    let table = revel_side(
        &directory.join("t"),
        &[("C", "V", "L", "0.016", ids), ("T", "I", "L", "0.029", ids)],
    );
    let valine = ("C", "M/V", "");
    let reference = revel_rows(&[valine, ("T", "M/L", "0.016")], TWO_ALTERNATES);
    let candidate = revel_rows(&[valine, ("T", "M/L", "")], TWO_ALTERNATES);
    assert_eq!(
        compare_revel(&reference, &candidate, &table, &[REVEL_OTHER]),
        green(REVEL_OTHER)
    );
    assert_eq!(
        compare_revel(&reference, &candidate, &table, &[REVEL_SIBLING]),
        red()
    );
    assert_eq!(
        compare_revel(
            &reference,
            &candidate,
            &table,
            &[REVEL_SIBLING, REVEL_OTHER]
        ),
        green(REVEL_OTHER)
    );
    // Without the correction the candidate writes its own alternate's row: a value the
    // corrected walk does not give stays a difference.
    assert_eq!(
        compare_revel(
            &reference,
            &revel_rows(&[valine, ("T", "M/L", "0.029")], TWO_ALTERNATES),
            &table,
            &[REVEL_OTHER]
        ),
        red()
    );
    let named: Vec<String> = ["0.016", ""]
        .iter()
        .map(|score| revel_rows(&[valine, ("T", "M/L", score)], "two_alternates"))
        .collect();
    let one_line = ["21\t13610388\ttwo_alternates\tA\tC,T\t.\t.\t."];
    assert_eq!(
        compare_revel_lines(
            &named[0],
            &named[1],
            &table,
            &[REVEL_OTHER],
            Some(&one_line)
        ),
        green(REVEL_OTHER)
    );
    // A record of one alternate cannot have met the C row: its value is no row's of its
    // own.
    assert_eq!(
        compare_revel(
            &revel_rows(&[("T", "M/L", "0.016")], ONE_ALTERNATE),
            &revel_rows(&[("T", "M/L", "")], ONE_ALTERNATE),
            &table,
            &[REVEL_OTHER]
        ),
        red()
    );
    // A substitution of two bases is matched on its trimmed alleles
    // (Utils/Sequence.pm:1190-1209).
    let two_bases = two_base_rows("21_13610387_GA/GT");
    assert_eq!(
        compare_revel(&two_bases[0], &two_bases[1], &table, &[REVEL_OTHER]),
        green(REVEL_OTHER)
    );
}

#[test]
fn the_other_reference_amino_acid_class_reads_a_one_base_row_without_its_record() {
    // Under an identifier that no input line, or more than one, carries, the rows do not
    // state their record. A one-base row is then read without it: the walk can have reached
    // the first passing row of the row's own allele or, before it, the first of another
    // allele, whichever alternates the record holds. The class takes the row where every
    // such row that carries the reference's value names another reference amino acid, and
    // leaves it where one of them names the row's own.
    let directory = scratch("revel-without-the-record");
    let ids = Some(TRANSCRIPT);
    let reference = revel_rows(&[("T", "M/L", "0.029")], "shared_name");
    let candidate = revel_rows(&[("T", "M/L", "")], "shared_name");
    let split = [
        "21\t13610388\tshared_name\tA\tC,T\t.\t.\t.",
        "21\t13610388\tshared_name\tA\tT\t.\t.\t.",
    ];
    let every_other = revel_side(
        &directory.join("every"),
        &[("C", "V", "L", "0.016", ids), ("T", "I", "L", "0.029", ids)],
    );
    let another_name = ["21\t13610388\tanother_name\tA\tT\t.\t.\t."];
    for lines in [None, Some(&split[..]), Some(&another_name[..])] {
        assert_eq!(
            compare_revel_lines(&reference, &candidate, &every_other, &[REVEL_OTHER], lines),
            green(REVEL_OTHER),
            "{lines:?}"
        );
    }
    // The same score on a row of the row's own pair, of another allele and earlier in the
    // file: whether the walk met it depends on the record, so the row is left.
    let either = revel_side(
        &directory.join("either"),
        &[("C", "M", "L", "0.029", ids), ("T", "I", "L", "0.029", ids)],
    );
    assert_eq!(
        compare_revel_lines(
            &reference,
            &candidate,
            &either,
            &[REVEL_OTHER],
            Some(&split)
        ),
        red()
    );
    assert_eq!(
        compare_revel_lines(
            &reference,
            &candidate,
            &either,
            &[REVEL_OTHER],
            Some(&["21\t13610388\tshared_name\tA\tT\t.\t.\t."])
        ),
        green(REVEL_OTHER)
    );
    assert_eq!(
        compare_revel_lines(
            &reference,
            &candidate,
            &either,
            &[REVEL_OTHER],
            Some(&["21\t13610388\tshared_name\tA\tC,T\t.\t.\t."])
        ),
        red()
    );
    // A row of another allele after the own allele's first passing row is never reached.
    let after = revel_side(
        &directory.join("after"),
        &[("T", "I", "L", "0.029", ids), ("C", "M", "L", "0.029", ids)],
    );
    assert_eq!(
        compare_revel_lines(&reference, &candidate, &after, &[REVEL_OTHER], Some(&split)),
        green(REVEL_OTHER)
    );
    // A row that is no one-base row at one position is not read without its record.
    let spanning = two_base_rows("shared_name");
    assert_eq!(
        compare_revel_lines(
            &spanning[0],
            &spanning[1],
            &every_other,
            &[REVEL_OTHER],
            Some(&split)
        ),
        red()
    );
}

#[test]
fn a_row_both_revel_classes_accept_is_the_sibling_classs() {
    // A made-up table in which the C row's score stands twice, under another reference
    // amino acid and under the row's own: both engines write it on the C row, and the
    // reference writes it on the T row as well. The T row's value is a sibling allele's and
    // the score of a row of another reference amino acid at once; named together, the class
    // that stands first among the field classes takes it.
    let directory = scratch("revel-both-classes");
    let ids = Some(TRANSCRIPT);
    let table = revel_side(
        &directory.join("t"),
        &[("C", "V", "L", "0.025", ids), ("C", "M", "L", "0.025", ids)],
    );
    let leucine = ("C", "M/L", "0.025");
    let reference = revel_rows(&[leucine, ("T", "M/L", "0.025")], TWO_ALTERNATES);
    let candidate = revel_rows(&[leucine, ("T", "M/L", "")], TWO_ALTERNATES);
    assert_eq!(
        compare_revel(&reference, &candidate, &table, &[REVEL_OTHER]),
        green(REVEL_OTHER)
    );
    assert_eq!(
        compare_revel(&reference, &candidate, &table, &[REVEL_SIBLING]),
        green(REVEL_SIBLING)
    );
    assert_eq!(
        compare_revel(
            &reference,
            &candidate,
            &table,
            &[REVEL_OTHER, REVEL_SIBLING]
        ),
        green(REVEL_SIBLING)
    );
    let classes: Vec<&str> = tables()
        .field_classes()
        .iter()
        .map(|(klass, _)| *klass)
        .collect();
    let place = |klass: &str| classes.iter().position(|held| *held == klass);
    assert_eq!(place(REVEL_SIBLING).map(|at| at + 1), place(REVEL_OTHER));
    assert!(place(REVEL_SIBLING).is_some());
}

#[test]
fn the_other_reference_amino_acid_class_reads_the_tab_and_json_formats() {
    // The tab format carries the field as a column and the JSON format as the `revel`
    // member; the JSON record states its alleles in `allele_string`.
    let directory = scratch("revel-tab-and-json");
    let ids = Some(TRANSCRIPT);
    let one = count_map(&[(REVEL_OTHER, 1)]);
    let table = revel_side(
        &directory.join("t"),
        &[("C", "V", "L", "0.016", ids), ("T", "I", "L", "0.029", ids)],
    );
    let tab = |score: &str| {
        format!(
            "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tAmino_acids\tREVEL\n{TWO_ALTERNATES}\t21:13610388\tT\t{TRANSCRIPT}\tTranscript\tmissense_variant\tM/L\t{score}\n"
        )
    };
    let (report, documented) = compare_under(
        &parse_tab(&tab("0.016")).expect("the reference parses"),
        &parse_tab(&tab("-")).expect("the candidate parses"),
        &table,
        &[REVEL_OTHER],
        None,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    let inputs = lines_of(&["21\t13610388\ttwo_alternates\tA\tC,T\t.\t.\t."]);
    let record = |score: Option<&str>| {
        let revel = score.map_or_else(String::new, |score| format!(", \"revel\": {score}"));
        let entry = format!(
            "{{\"variant_allele\": \"T\", \"transcript_id\": \"{TRANSCRIPT}\", \"consequence_terms\": [\"missense_variant\"], \"amino_acids\": \"M/L\"{revel}}}"
        );
        json_record(
            &inputs[0],
            "two_alternates",
            "\"seq_region_name\": \"21\", \"start\": 13610388, \"end\": 13610388, \"allele_string\": \"A/C/T\"",
            &[entry],
        )
    };
    let scored = parse_json(&record(Some("0.016")), &inputs).expect("the reference parses");
    let unscored = parse_json(&record(None), &inputs).expect("the candidate parses");
    let (report, documented) = compare_under(&scored, &unscored, &table, &[REVEL_OTHER], None);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    let (report, documented) =
        compare_under(&scored, &unscored, &Side::default(), &[REVEL_OTHER], None);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), count_map(&[]));
}

/// A made-up record of two alternates, C and T, on one transcript, both giving M/L: the
/// reference writes one value on both alleles' rows. The value is the allele's whose row
/// the candidate fills too, so with the candidate's C row filled the empty T row is the
/// class's, and with its T row filled the empty C row is. A candidate that leaves both
/// rows empty differs on the named allele's own row as well: no row of the candidate
/// carries the value, neither row has a sibling the table names, and both stay
/// differences. A caller that hands no rows of the candidate's gets no acceptance.
fn a_sibling_allele_class_refuses_the_named_alleles_own_row(klass: &str, value: &str) {
    let rows = |filled: &[&str]| {
        let rows: Vec<String> = ["C", "T"]
            .iter()
            .map(|allele| {
                let extra = if filled.contains(allele) {
                    format!("IMPACT=MODERATE;STRAND=1;{value}")
                } else {
                    "IMPACT=MODERATE;STRAND=1".to_string()
                };
                coding_row(
                    "21_5000100_A/C/T",
                    "21:5000100",
                    allele,
                    "ENST00000000041",
                    "M/L",
                    &extra,
                )
            })
            .collect();
        default_text(&rows)
    };
    let reference = rows(&["C", "T"]);
    for (kept, set_aside) in [("C", "T"), ("T", "C")] {
        let (report, documented) = compare_default(&reference, &rows(&[kept]), &[klass]);
        assert_eq!(verdict(&report), "green", "{kept}");
        assert_eq!(applied(&documented), count_map(&[(klass, 1)]), "{kept}");
        assert!(
            texts(at(&report, &["mask_examples", klass]))[0].starts_with(&format!(
                "21:5000100|{set_aside}|ENST00000000041|Transcript: "
            )),
            "{kept}"
        );
    }
    let (report, documented) = compare_default(&reference, &rows(&[]), &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), count_map(&[]));
    let (name, filled) = value
        .split(';')
        .next()
        .and_then(|first| first.split_once('='))
        .expect("a field and its value");
    let reference_rows = parse_default(&reference).expect("the reference parses");
    let candidate_rows = parse_default(&rows(&["C"])).expect("the candidate parses");
    let (sibling, own) = (&reference_rows[1], &reference_rows[0]);
    let side = Side::default();
    let asked = |e: &Entry, a: &Entry, candidate_rows: Option<&[Entry]>| {
        let record_rows = refs(&reference_rows);
        let candidate_rows = candidate_rows.map(refs);
        let context = RuleContext {
            record_rows: Some(&record_rows),
            candidate_rows: candidate_rows.as_deref(),
            declared: &[],
            assembly: None,
            dbnsfp_columns: &[],
            side: &side,
        };
        field_class_accepts(klass, e, name, filled, ABSENT, Some(a), &context)
            .expect("the class is asked")
    };
    assert!(asked(sibling, &candidate_rows[1], Some(&candidate_rows)));
    assert!(!asked(sibling, &candidate_rows[1], None));
    let emptied = parse_default(&rows(&[])).expect("the candidate parses");
    assert!(!asked(own, &emptied[0], Some(&emptied)));
}

#[test]
fn the_alphamissense_sibling_class_refuses_the_named_alleles_own_row() {
    a_sibling_allele_class_refuses_the_named_alleles_own_row(
        ALPHAMISSENSE_SIBLING,
        "am_class=likely_benign;am_pathogenicity=0.1079",
    );
}

#[test]
fn the_revel_sibling_class_refuses_the_named_alleles_own_row() {
    a_sibling_allele_class_refuses_the_named_alleles_own_row(REVEL_SIBLING, "REVEL=0.025");
}

#[test]
fn the_splicevault_sibling_class_refuses_the_named_alleles_own_row() {
    a_sibling_allele_class_refuses_the_named_alleles_own_row(
        "splicevault_sibling_allele_rows",
        "SpliceVault_site_type=donor;SpliceVault_site_pos=5",
    );
}

/// What a run says that ends on an input no reading of the comparison expects.
fn uncaught<T>(run: impl FnOnce() -> T) -> String {
    let payload = std::panic::catch_unwind(std::panic::AssertUnwindSafe(run))
        .err()
        .expect("the run ends");
    payload
        .downcast_ref::<Uncaught>()
        .expect("an ending of the run's own")
        .0
        .clone()
}

/// A nine-column table of one row at 21:500, written under a directory of `name`.
fn one_row_table(name: &str) -> PathBuf {
    nine_column_table(&scratch(name), &["21\t100\t500\tA\tG\tQ\tH\t0.742\tENST1"])
}

/// The index of `table` rewritten with the whole number `value` at the byte `at` of its
/// inflated text.
fn patch_index(table: &Path, at: usize, value: i32) {
    let mut raw = Vec::new();
    let packed = std::fs::read(index_path(table)).expect("the index is read");
    std::io::Read::read_to_end(&mut flate2::read::GzDecoder::new(&packed[..]), &mut raw)
        .expect("the index inflates");
    raw[at..at + 4].copy_from_slice(&value.to_le_bytes());
    std::fs::write(index_path(table), gzipped(&raw)).expect("the index is written");
}

/// The bytes of the index header that hold the sequence column, the position column and
/// the character that opens a header line.
const INDEX_SEQUENCE_COLUMN: usize = 12;
const INDEX_POSITION_COLUMN: usize = 16;
const INDEX_META: usize = 24;

#[test]
fn a_table_file_the_interpreter_does_not_open_ends_the_run_as_its_exception_does() {
    // What the index raises is the comparator's error; of what the table file raises only
    // a refused deflate stream is: a table that is missing or a directory beside its index
    // ends the run with the text of the `OSError`.
    let missing = one_row_table("revel-missing-table");
    std::fs::remove_file(&missing).expect("the table is removed");
    assert_eq!(
        uncaught(|| RevelTable::new(&missing, false)),
        format!(
            "[Errno 2] No such file or directory: '{}'",
            missing.display()
        )
    );
    std::fs::create_dir(&missing).expect("a directory in the table's place");
    assert_eq!(
        uncaught(|| RevelTable::new(&missing, false)),
        format!("[Errno 21] Is a directory: '{}'", missing.display())
    );
    // A table read once is opened again for every span its index names a chromosome for,
    // and for no other span.
    let read = one_row_table("revel-table-removed-after-reading");
    let table = RevelTable::new(&read, false).expect("the table reads");
    std::fs::remove_file(&read).expect("the table is removed");
    assert!(rows_at(&table, "22", 500, 500).is_empty());
    assert!(rows_at(&table, "21", 500, 499).is_empty());
    assert_eq!(
        uncaught(|| table.rows("21", 500, 500)),
        format!("[Errno 2] No such file or directory: '{}'", read.display())
    );
}

#[test]
fn a_block_whose_extra_field_or_size_is_cut_short_ends_the_run_as_the_interpreters_error_does() {
    // The size subfield of a block is read where the extra field ends before it, and a
    // block size smaller than the block's own header is a length the interpreter does not
    // read by: neither exception is the comparator's.
    let cut = one_row_table("revel-cut-size-subfield");
    let mut bytes = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x04\x00BC\x02\x00".to_vec();
    bytes.extend_from_slice(&[0; 16]);
    std::fs::write(&cut, &bytes).expect("the table is written");
    assert_eq!(
        uncaught(|| RevelTable::new(&cut, false)),
        "unpack_from requires a buffer of at least 6 bytes for unpacking 2 bytes \
         at offset 4 (actual buffer size is 4)"
    );
    let small = one_row_table("revel-small-block-size");
    let mut bytes = std::fs::read(&small).expect("the table is read");
    bytes[16..18].copy_from_slice(&5u16.to_le_bytes());
    std::fs::write(&small, &bytes).expect("the table is written");
    assert_eq!(
        uncaught(|| RevelTable::new(&small, false)),
        "read length must be non-negative or -1"
    );
    // A size as large as the header reads no byte, which is no deflate stream.
    bytes[16..18].copy_from_slice(&25u16.to_le_bytes());
    std::fs::write(&small, &bytes).expect("the table is written");
    assert_eq!(
        refusal(&small),
        format!(
            "--revel-table: cannot read {}: Error -5 while decompressing data: incomplete or \
             truncated stream",
            small.display()
        )
    );
    // A size one byte short of the header's reads to the end of the file, the block's
    // stream first, and the block after it is sought one byte past that size.
    bytes[16..18].copy_from_slice(&24u16.to_le_bytes());
    std::fs::write(&small, &bytes).expect("the table is written");
    assert_eq!(
        refusal(&small),
        format!("--revel-table: {} is not a BGZF file", small.display())
    );
    // A member without the size subfield is no block.
    bytes[12..14].copy_from_slice(b"XY");
    std::fs::write(&small, &bytes).expect("the table is written");
    assert_eq!(
        refusal(&small),
        format!("--revel-table: {} is not a BGZF file", small.display())
    );
}

#[test]
fn the_index_names_the_columns_and_the_header_character_as_the_interpreter_reads_them() {
    // A column the index numbers 0 is the row's last cell, as an index of -1 counts from
    // the end of a list; one further below ends the run where a row is too short for it.
    let last_cell = nine_column_table(
        &scratch("revel-last-cell-column"),
        &["21\t100\t500\tA\tG\tQ\tH\t0.742\t500"],
    );
    patch_index(&last_cell, INDEX_POSITION_COLUMN, 0);
    let table = RevelTable::new(&last_cell, false).expect("the table reads");
    assert_eq!(table._position_column, -1);
    assert_eq!(
        rows_at(&table, "21", 500, 500),
        vec![row(["500", "A", "G", "Q", "H", "0.742", "500"])]
    );
    let short_row = nine_column_table(&scratch("revel-column-below-the-row"), &["21"]);
    patch_index(&short_row, INDEX_POSITION_COLUMN, -3);
    let table = RevelTable::new(&short_row, false).expect("the table reads");
    assert_eq!(
        uncaught(|| table.rows("21", 500, 500)),
        "list index out of range"
    );
    // The header character is any code point; a number that is none ends the run.
    let no_character = one_row_table("revel-meta-outside-the-code-points");
    patch_index(&no_character, INDEX_META, 0x110000);
    assert_eq!(
        uncaught(|| RevelTable::new(&no_character, false)),
        "chr() arg not in range(0x110000)"
    );
    // Under another header character the `#` line is a row, so no line is a header.
    let other_character = one_row_table("revel-another-meta");
    patch_index(&other_character, INDEX_META, i32::from(b'@'));
    assert_eq!(
        refusal(&other_character),
        format!(
            "--revel-table: {} is not a REVEL table (a header of 7, 8 or 9 columns)",
            other_character.display()
        )
    );
    // A sequence column past a row's cells leaves the row out.
    let wide = one_row_table("revel-sequence-column-past-the-row");
    patch_index(&wide, INDEX_SEQUENCE_COLUMN, 40);
    let table = RevelTable::new(&wide, false).expect("the table reads");
    assert!(rows_at(&table, "21", 500, 500).is_empty());
}

#[test]
fn the_header_is_the_last_of_the_leading_header_lines() {
    // Every leading line that opens with the header character is read as the header in
    // turn, so a comment line between the column names and the first row is the header.
    let directory = scratch("revel-leading-header-lines");
    let lines = lines_of(&["21\t100\t500\tA\tG\tQ\tH\t0.742\tENST1"]);
    let two_lines = format!("##a table\n{REVEL_TABLE_HEADER}");
    let named_last = write_revel_table(&directory, &lines, 3, &two_lines, "named-last.tsv.gz");
    let table = RevelTable::new(&named_last, false).expect("the table reads");
    assert_eq!(table._columns, (3, 4, 5, 6, 7, Some(8)));
    assert_eq!(rows_at(&table, "21", 500, 500).len(), 1);
    let comment_last = format!("{REVEL_TABLE_HEADER}\n#a comment");
    let refused = write_revel_table(&directory, &lines, 3, &comment_last, "comment-last.tsv.gz");
    assert_eq!(
        refusal(&refused),
        format!(
            "--revel-table: {} is not a REVEL table (a header of 7, 8 or 9 columns)",
            refused.display()
        )
    );
}

#[test]
fn a_table_of_concatenated_streams_gives_the_rows_behind_the_first_stream() {
    // Three BGZF streams one after the other, as compressed tables joined byte for byte
    // are: each is one block of text and closes with an empty block. The rows behind a
    // closing block are read as the rows before it are, under an index that numbers the
    // block behind an empty block by the empty block's offset, as htslib 1.9 writes one,
    // and under an index that numbers every block by its own offset.
    let directory = scratch("revel-concatenated-streams");
    let lines = [
        "21\t100\t500\tA\tG\tQ\tH\t0.1\tENST1",
        "21\t101\t501\tA\tG\tQ\tH\t0.2\tENST1",
        "21\t102\t502\tA\tG\tQ\tH\t0.3\tENST1",
        "21\t103\t503\tA\tG\tQ\tH\t0.4\tENST1",
        "22\t104\t504\tA\tG\tQ\tH\t0.5\tENST2",
        "22\t105\t505\tA\tG\tQ\tH\t0.6\tENST2",
    ];
    let text =
        |lines: &[&str]| -> String { lines.iter().map(|line| format!("{line}\n")).collect() };
    let closing_block = bgzf_block(b"");
    let stream = |text: String| [bgzf_block(text.as_bytes()), closing_block.clone()].concat();
    let first = stream(text(&[REVEL_TABLE_HEADER, lines[0], lines[1]]));
    let second = stream(text(&lines[2..5]));
    let third = stream(text(&lines[5..]));
    let path = directory.join("revel.tsv.gz");
    std::fs::write(
        &path,
        [first.as_slice(), second.as_slice(), third.as_slice()].concat(),
    )
    .expect("the table is written");
    // The file offsets of the first stream's closing block, of the second stream's text and
    // of the file's last block; then the place in the first stream's text where its rows
    // begin and the place in the second stream's where its rows of 21 end.
    let closing = (first.len() - closing_block.len()) as u64;
    let behind = first.len() as u64;
    let last = (first.len() + second.len() + third.len() - closing_block.len()) as u64;
    let after_header = REVEL_TABLE_HEADER.len() as u64 + 1;
    let after_21 = text(&lines[2..4]).len() as u64;
    let rows_of = |lines: &[&str]| -> Vec<RevelRow> {
        lines
            .iter()
            .map(|line| {
                let cells: Vec<&str> = line.split('\t').collect();
                row([
                    cells[2], cells[3], cells[4], cells[5], cells[6], cells[7], cells[8],
                ])
            })
            .collect()
    };
    for (numbered, block) in [("by the empty block", closing), ("by its own", behind)] {
        let split = (block << 16) | after_21;
        let chunks = [
            ("21".to_string(), after_header, split),
            ("22".to_string(), split, last << 16),
        ];
        write_tabix_index(&path, &chunks, 3);
        let table = RevelTable::new(&path, false).expect("the table reads");
        assert_eq!(table._columns, (3, 4, 5, 6, 7, Some(8)), "{numbered}");
        for (chromosome, start, end, expected) in [
            ("21", 500, 503, &lines[..4]),
            ("21", 502, 503, &lines[2..4]),
            ("22", 504, 505, &lines[4..]),
            ("22", 505, 505, &lines[5..]),
        ] {
            assert_eq!(
                rows_at(&table, chromosome, start, end),
                rows_of(expected),
                "{numbered}: {chromosome} {start} to {end}"
            );
        }
    }
}

#[test]
fn a_chunk_that_ends_in_the_text_behind_an_empty_block_gives_no_row_past_its_end() {
    // A chunk that runs over the first stream's closing block and ends in the text behind
    // it gives the rows before its end and no row after it, under an index that numbers
    // that text by the closing block's offset and under an index that numbers it by its
    // own.
    let directory = scratch("revel-chunk-end-behind-an-empty-block");
    let lines = [
        "21\t100\t500\tA\tG\tQ\tH\t0.1\tENST1",
        "21\t101\t501\tA\tG\tQ\tH\t0.2\tENST1",
        "21\t102\t502\tA\tG\tQ\tH\t0.3\tENST1",
    ];
    let text =
        |lines: &[&str]| -> String { lines.iter().map(|line| format!("{line}\n")).collect() };
    let closing_block = bgzf_block(b"");
    let stream = |text: String| [bgzf_block(text.as_bytes()), closing_block.clone()].concat();
    let first = stream(text(&[REVEL_TABLE_HEADER, lines[0]]));
    let second = stream(text(&lines[1..]));
    let path = directory.join("revel.tsv.gz");
    std::fs::write(&path, [first.as_slice(), second.as_slice()].concat())
        .expect("the table is written");
    // The file offsets of the first stream's closing block and of the second stream's
    // text; then the place in the first stream's text where its row begins and the place
    // in the second stream's where its first row ends.
    let closing = (first.len() - closing_block.len()) as u64;
    let behind = first.len() as u64;
    let after_header = REVEL_TABLE_HEADER.len() as u64 + 1;
    let after_first = text(&lines[1..2]).len() as u64;
    for (numbered, block) in [("by the empty block", closing), ("by its own", behind)] {
        let chunks = [("21".to_string(), after_header, (block << 16) | after_first)];
        write_tabix_index(&path, &chunks, 3);
        let table = RevelTable::new(&path, false).expect("the table reads");
        let positions: Vec<String> = rows_at(&table, "21", 500, 502)
            .into_iter()
            .map(|row| row.position)
            .collect();
        assert_eq!(positions, ["500", "501"], "{numbered}");
    }
}

#[test]
fn the_index_is_inflated_whole_and_refused_in_the_interpreters_words() {
    // The index is one gzip text inflated in one call: its checksum is refused without the
    // two sums a file read through `gzip.open` names, zero bytes after it are padding, and
    // any other byte after it opens what is no gzip member.
    let table = one_row_table("revel-index-inflated-whole");
    let index = index_path(&table);
    let refused = |reason: &str| {
        format!(
            "--revel-table: cannot read the tabix index {}: {reason}",
            index.display()
        )
    };
    let packed = std::fs::read(&index).expect("the index is read");
    let with = |more: &[u8]| [&packed[..], more].concat();
    std::fs::write(&index, with(&[0, 0, 0])).expect("the index is written");
    assert!(RevelTable::new(&table, false).is_ok());
    std::fs::write(&index, with(b"xx")).expect("the index is written");
    assert_eq!(refusal(&table), refused("Not a gzipped file (b'xx')"));
    let mut wrong_sum = packed.clone();
    let sum = wrong_sum.len() - 8;
    wrong_sum[sum] ^= 0xff;
    std::fs::write(&index, &wrong_sum).expect("the index is written");
    assert_eq!(refusal(&table), refused("CRC check failed"));
    std::fs::write(&index, &packed[..packed.len() - 9]).expect("the index is written");
    assert_eq!(
        refusal(&table),
        refused("Compressed file ended before the end-of-stream marker was reached")
    );
    std::fs::write(&index, gzipped(b"BAI\x01")).expect("the index is written");
    assert_eq!(refusal(&table), refused("not a tabix index"));
    std::fs::write(&index, gzipped(b"TBI\x01\x01\x00")).expect("the index is written");
    assert_eq!(
        refusal(&table),
        refused(
            "unpack_from requires a buffer of at least 36 bytes for unpacking 32 bytes at \
             offset 4 (actual buffer size is 6)"
        )
    );
}

/// A CSI index, the layout `tabix -C` writes, of the table `write_revel_table` wrote at
/// `table` under `header` with every row on `chromosome` and the position in the 1-based
/// `position_column`: the bin that spans the sequence holds one chunk over every row.
fn csi_index(table: &Path, header: &str, chromosome: &str, position_column: i32) -> Vec<u8> {
    let first = bgzf_block(format!("{header}\n").as_bytes()).len() as u64;
    let length = std::fs::metadata(table).expect("the table is there").len();
    let end = length - bgzf_block(b"").len() as u64;
    let mut names = chromosome.as_bytes().to_vec();
    names.push(0);
    let mut index = b"CSI\x01".to_vec();
    // The smallest bin's shift and the levels below the top bin, then, their length before
    // them, the seven numbers of a tabix index's header and the names.
    for number in [
        14,
        5,
        28 + names.len() as i32,
        0,
        1,
        position_column,
        position_column,
        i32::from(b'#'),
        0,
        names.len() as i32,
    ] {
        index.extend_from_slice(&number.to_le_bytes());
    }
    index.extend_from_slice(&names);
    // One sequence of one bin: the bin's number, the offset its rows begin at, one chunk.
    index.extend_from_slice(&1i32.to_le_bytes());
    index.extend_from_slice(&1i32.to_le_bytes());
    index.extend_from_slice(&0u32.to_le_bytes());
    index.extend_from_slice(&(first << 16).to_le_bytes());
    index.extend_from_slice(&1i32.to_le_bytes());
    index.extend_from_slice(&(first << 16).to_le_bytes());
    index.extend_from_slice(&(end << 16).to_le_bytes());
    [bgzf_block(&index), bgzf_block(b"")].concat()
}

#[test]
fn a_table_under_a_csi_index_alone_is_refused_for_want_of_its_tbi_index() {
    // The reader looks for `<table>.tbi`, the index Ensembl VEP's check demands beside a
    // plugin's table (BaseVepTabixPlugin.pm:273), and for no other. A comparison that names
    // a made-up table with that index reads it. With a CSI index in its place the table is
    // refused, naming the `.tbi` it lacks: the comparison ends with exit code 2 and no
    // verdict, and so does the `field-class` command, whose `--spliceai-table` goes to the
    // same reader, which asks for the index before it reads a line of the table.
    let directory = scratch("revel-csi-index-alone");
    let lines = lines_of(&["21\t100\t500\tA\tG\tQ\tH\t0.742\tENST1"]);
    let name = "revel.synthetic.tsv.gz";
    let table = write_revel_table(&directory, &lines, 3, REVEL_TABLE_HEADER, name);
    let rows = revel_rows(&[("T", "Q/H", "0.742")], ONE_ALTERNATE);
    let expected = planted(&directory, "expected.txt", &rows);
    let actual = planted(&directory, "actual.txt", &rows);
    let ended = |args: &[&std::ffi::OsStr]| {
        let ran = support::run(&directory, args);
        (ran.code, ran.out(), ran.err())
    };
    let compared = || {
        ended(&[
            "--format".as_ref(),
            "default".as_ref(),
            "--expected".as_ref(),
            expected.as_os_str(),
            "--actual".as_ref(),
            actual.as_os_str(),
            "--revel-table".as_ref(),
            table.as_os_str(),
        ])
    };
    assert_eq!(
        compared(),
        (0, "default: green\n".to_string(), String::new())
    );

    let csi = PathBuf::from(format!("{}.csi", table.display()));
    std::fs::write(&csi, csi_index(&table, REVEL_TABLE_HEADER, "21", 3))
        .expect("the index is written");
    std::fs::remove_file(index_path(&table)).expect("the index is removed");
    let unread = format!(
        "cannot read the tabix index {}: No such file or directory",
        index_path(&table).display()
    );
    assert_eq!(refusal(&table), format!("--revel-table: {unread}"));
    let error_line = |option: &str| format!("ERROR: [compare_vep_formats] {option}: {unread}\n");
    assert_eq!(compared(), (2, String::new(), error_line("--revel-table")));
    let questions = planted(&directory, "questions.tsv", "");
    let answers = directory.join("answers.txt");
    let asked = ended(&[
        "field-class".as_ref(),
        "--class".as_ref(),
        "spliceai_multi_entry_cell".as_ref(),
        "--spliceai-table".as_ref(),
        table.as_os_str(),
        "--questions".as_ref(),
        questions.as_os_str(),
        "--answers".as_ref(),
        answers.as_os_str(),
    ]);
    assert_eq!(asked, (2, String::new(), error_line("--spliceai-table")));
    assert!(!answers.exists());
}
