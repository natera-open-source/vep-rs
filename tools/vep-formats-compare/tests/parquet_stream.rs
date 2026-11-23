// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A Parquet comparison through the stream of record groups: the rows both outputs hold
//! are kept back whatever their order, a group that differs is handed over, and the report
//! is the one the comparison of the whole outputs writes. Every row is written here.

mod compare_common;
mod parquet_writer;

use std::collections::BTreeSet;
use std::io::Read;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

use compare_common::*;
use indexmap::IndexMap;
use parquet_writer::plant;
use parquet_writer::plant_tab_rows;
use parquet_writer::Stored;
use parquet_writer::AS_THE_WRITER;
use vep_formats_compare::compare_files::compare_files_through;
use vep_formats_compare::error::Error;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::parquet::parquet_as_tab;
use vep_formats_compare::parse::parquet::CanonicalText;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::remainder::FrontEnd;
use vep_formats_compare::remainder::STREAM;
use vep_formats_compare::rules::Side;
use vep_formats_compare::stream::Text;

const COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tIMPACT";

/// One transcript row of record `record`, at a Location of the record's own.
fn planted_row(record: usize, feature: &str, terms: &str) -> String {
    let allele = &"ACGT"[record % 4..record % 4 + 1];
    format!(
        "rec{record}\t21:{}\t{allele}\tENSG1\t{feature}\tTranscript\t{terms}\tMODIFIER",
        1000 + record
    )
}

/// Two transcript rows for each of `records` records.
fn planted_rows(records: usize) -> Vec<String> {
    (0..records)
        .flat_map(|record| {
            ["ENST1", "ENST2"].map(|feature| planted_row(record, feature, "intron_variant"))
        })
        .collect()
}

fn tab_text(rows: &[String]) -> String {
    let mut text = format!("## ENSEMBL VARIANT EFFECT PREDICTOR\n{COLUMNS}\n");
    for row in rows {
        text.push_str(row);
        text.push('\n');
    }
    text
}

/// A Parquet output of one part file for each of `parts`, in their order.
fn planted_parts(directory: &Path, name: &str, parts: &[&[String]]) -> PathBuf {
    let out = directory.join(name);
    std::fs::create_dir_all(&out).expect("the output's directory is made");
    for (at, rows) in parts.iter().enumerate() {
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        plant_tab_rows(&out.join(format!("part-{at:05}.parquet")), COLUMNS, &rows);
    }
    out
}

/// The Parquet comparison of `expected` and `actual` under `documented`: through
/// `front_end`, or of the whole outputs without one.
fn compared(
    front_end: Option<&mut FrontEnd>,
    expected: &Path,
    actual: &Path,
    documented: &mut Documented,
) -> Result<Value, Error> {
    let masks: BTreeSet<String> = documented.masks.clone();
    compare_files_through(
        front_end,
        "parquet",
        expected,
        actual,
        None,
        documented,
        &masks,
        &IndexMap::new(),
        &IndexMap::new(),
        4,
        false,
        None,
        None,
        None,
        &mut Side::default(),
    )
}

/// The report of the comparison through the stream of record groups and what the stream
/// says it did, the report held to the report of the whole outputs byte for byte.
fn through_the_stream(expected: &Path, actual: &Path) -> (Value, String) {
    let whole =
        compared(None, expected, actual, &mut no_documented()).expect("the whole outputs compare");
    let mut front_end = FrontEnd::new(STREAM);
    let handed = compared(Some(&mut front_end), expected, actual, &mut no_documented())
        .expect("the outputs compare through the stream");
    let options = pyjson::DumpOptions::default();
    assert_eq!(
        pyjson::dumps(&handed, &options),
        pyjson::dumps(&whole, &options)
    );
    (handed, front_end.note)
}

/// Forty records whose two rows both outputs hold, the candidate's second half in its
/// first part file: no row is handed over, and the report counts them all.
#[test]
fn a_parquet_comparison_keeps_back_the_rows_both_outputs_hold() {
    let directory = scratch("parquet_stream_kept_back");
    let rows = planted_rows(40);
    let expected = planted(&directory, "tab.txt", &tab_text(&rows));
    let out = planted_parts(&directory, "out.parquet", &[&rows[40..], &rows[..40]]);
    let (report, note) = through_the_stream(&expected, &out);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(number(&report, &["counts", "compared"]), 80);
    assert_eq!(
        note,
        "80 of 80 reference rows were kept back by the stream of record groups, \
         0 of 40 record groups handed over"
    );
}

/// One row of one record differs: its record's group alone is handed over, and the report
/// names the row as the report of the whole outputs does.
#[test]
fn a_row_that_differs_hands_its_group_over_alone() {
    let directory = scratch("parquet_stream_one_group");
    let rows = planted_rows(40);
    let mut candidate = rows.clone();
    candidate[13] = planted_row(6, "ENST2", "splice_region_variant,intron_variant");
    let expected = planted(&directory, "tab.txt", &tab_text(&rows));
    let out = planted_parts(
        &directory,
        "out.parquet",
        &[&candidate[40..], &candidate[..40]],
    );
    let (report, note) = through_the_stream(&expected, &out);
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["counts", "compared"]), 80);
    assert_eq!(number(&report, &["columns", "Consequence", "differing"]), 1);
    assert_eq!(
        note,
        "78 of 80 reference rows were kept back by the stream of record groups, \
         1 of 40 record groups handed over"
    );
}

/// What a comparison came to: its report as the tool writes it, or its error's line.
fn outcome(result: Result<Value, Error>) -> Result<String, String> {
    result
        .map(|report| pyjson::dumps(&report, &pyjson::DumpOptions::default()))
        .map_err(|error| error.to_string())
}

/// The outcome of the comparison through the stream of record groups, held to the outcome
/// of the comparison of the whole outputs, and what the stream says it did.
fn both_paths(expected: &Path, actual: &Path) -> (Result<String, String>, String) {
    let whole = outcome(compared(None, expected, actual, &mut no_documented()));
    let mut front_end = FrontEnd::new(STREAM);
    let handed = outcome(compared(
        Some(&mut front_end),
        expected,
        actual,
        &mut no_documented(),
    ));
    assert_eq!(handed, whole);
    (handed, front_end.note)
}

fn kept_back(rows: usize, groups: usize, handed_over: usize) -> String {
    format!(
        "{} of {rows} reference rows were kept back by the stream of record groups, \
         {handed_over} of {groups} record groups handed over",
        rows - 2 * handed_over
    )
}

fn compared_whole(outside: &str) -> String {
    format!("the whole outputs were compared: {outside} is outside the front end")
}

/// The candidate holds each record's two rows in the other order: a record group whose
/// lines stand in another order is handed over, and the comparison, which reads no order
/// of a Parquet output, is green.
#[test]
fn rows_in_another_order_inside_a_record_are_handed_over_and_compare_green() {
    let directory = scratch("parquet_stream_reversed");
    let rows = planted_rows(12);
    let mut candidate = rows.clone();
    candidate.reverse();
    let expected = planted(&directory, "tab.txt", &tab_text(&rows));
    let out = planted_parts(&directory, "out.parquet", &[&candidate]);
    let (report, note) = through_the_stream(&expected, &out);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(number(&report, &["counts", "compared"]), 24);
    assert_eq!(note, kept_back(24, 12, 12));
}

/// Thirteen records differ in one row each, more than the examples a report shows. The
/// remainder is parsed in the order of the whole canonical form, the byte order of its
/// lines, so the examples are the ones the whole outputs give: `rec10` stands before
/// `rec2` there, and both outputs hold `rec2` first.
#[test]
fn the_examples_of_a_remainder_are_the_whole_canonical_form_s() {
    let directory = scratch("parquet_stream_examples");
    let rows = planted_rows(40);
    let mut candidate = rows.clone();
    for record in 2..15 {
        candidate[2 * record] = planted_row(record, "ENST1", "splice_region_variant");
    }
    let expected = planted(&directory, "tab.txt", &tab_text(&rows));
    let out = planted_parts(&directory, "out.parquet", &[&candidate]);
    let (report, note) = through_the_stream(&expected, &out);
    assert_eq!(verdict(&report), "red");
    assert_eq!(
        number(&report, &["columns", "Consequence", "differing"]),
        13
    );
    let examples = pyjson::dumps_compact(at(&report, &["columns", "Consequence", "examples"]));
    let places: Vec<usize> = ["21:1010|", "21:1011|", "21:1012|", "21:1013|", "21:1002|"]
        .iter()
        .map(|key| examples.find(key).unwrap_or(usize::MAX))
        .collect();
    assert!(
        places[..4].windows(2).all(|pair| pair[0] < pair[1]) && places[3] < usize::MAX,
        "{examples}"
    );
    assert_eq!(places[4], usize::MAX, "{examples}");
    assert_eq!(note, kept_back(80, 40, 13));
}

/// A row the candidate lacks, a record the reference lacks and two rows of one key that
/// differ are reported as the whole outputs report them: the three records' groups are
/// handed over, five of the reference's 21 rows.
#[test]
fn one_sided_rows_and_rows_of_one_key_are_reported_as_the_whole_outputs_report_them() {
    let directory = scratch("parquet_stream_one_sided");
    let rows = planted_rows(10);
    let mut reference = rows.clone();
    // Record 4 holds its second row twice in the reference.
    reference.insert(9, rows[9].clone());
    let mut candidate = rows.clone();
    candidate.insert(9, planted_row(4, "ENST2", "splice_region_variant"));
    // Record 7 lacks its first row in the candidate, which holds a record of its own.
    candidate.remove(15);
    candidate.push(planted_row(30, "ENST1", "intron_variant"));
    let expected = planted(&directory, "tab.txt", &tab_text(&reference));
    let out = planted_parts(
        &directory,
        "out.parquet",
        &[&candidate[8..], &candidate[..8]],
    );
    let (report, note) = through_the_stream(&expected, &out);
    assert_eq!(verdict(&report), "red");
    assert_eq!(
        at(&report, &["one_sided"]),
        &json(r#"{"expected_only": 1, "actual_only": 1}"#)
    );
    assert_eq!(
        note,
        "16 of 21 reference rows were kept back by the stream of record groups, \
         3 of 11 record groups handed over"
    );
}

/// The reference spells a frequency with an exponent or with trailing zeros and the
/// Parquet file stores the number: the text the stream reads of each is the canonical
/// form, so the rows agree and are kept back, and another value still differs.
#[test]
fn frequencies_agree_as_numbers_in_the_text_the_stream_reads() {
    let directory = scratch("parquet_stream_frequencies");
    let columns = "#Uploaded_variation\tLocation\tAllele\tConsequence\tAF\tgnomADe_AF\tMAX_AF";
    let row = |record: usize, greatest: &str| {
        format!(
            "rs{record}\t21:{}\tA\tintergenic_variant\t1.47e-05\t0.0010,0.2770\t{greatest}",
            100 + record
        )
    };
    let rows: Vec<String> = (0..6).map(|record| row(record, "0.2770")).collect();
    let text = format!(
        "## ENSEMBL VARIANT EFFECT PREDICTOR\n{columns}\n{}\n",
        rows.join("\n")
    );
    let expected = planted(&directory, "tab.txt", &text);
    let plant_candidate = |rows: &[String]| {
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        let out = directory.join("out.parquet");
        plant_tab_rows(&out, columns, &rows);
        out
    };
    let (report, note) = through_the_stream(&expected, &plant_candidate(&rows));
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(
        note,
        "6 of 6 reference rows were kept back by the stream of record groups, \
         0 of 6 record groups handed over"
    );
    let mut other = rows.clone();
    other[3] = row(3, "0.2771");
    let (report, note) = through_the_stream(&expected, &plant_candidate(&other));
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "MAX_AF", "differing"]), 1);
    assert_eq!(
        note,
        "5 of 6 reference rows were kept back by the stream of record groups, \
         1 of 6 record groups handed over"
    );
}

/// A candidate that cannot be read ends the comparison with the adapter's line through
/// both paths, also beside a reference one of whose rows the comparison may refuse, a row
/// that holds a letter outside ASCII. Beside a reference no row of which may end the
/// comparison, the stream holds no line of the reference; beside the other it hands every
/// row over, so that the reference's own refusal would come first.
#[test]
fn a_candidate_that_cannot_be_read_ends_both_paths_on_the_adapter_s_line() {
    let directory = scratch("parquet_stream_absent");
    let rows = planted_rows(6);
    let absent = directory.join("absent.parquet");
    let refusal = format!(
        "the Parquet adapter failed on {}: ERROR: [parquet_to_vep_tab] {}: no such file or directory",
        absent.display(),
        absent.display()
    );
    let expected = planted(&directory, "tab.txt", &tab_text(&rows));
    let (ended, note) = both_paths(&expected, &absent);
    assert_eq!(ended, Err(refusal.clone()));
    assert_eq!(
        note,
        "the comparison ended on a candidate that is not read: no line of the reference was held"
    );
    let mut marked = rows.clone();
    marked[5] = planted_row(2, "ENST2", "intron_variant").replace("ENSG1", "ENSG\u{e9}");
    let expected = planted(&directory, "marked.txt", &tab_text(&marked));
    let (ended, note) = both_paths(&expected, &absent);
    assert_eq!(ended, Err(refusal));
    assert_eq!(
        note,
        "0 of 12 reference rows were kept back by the stream of record groups, \
         0 of 0 record groups handed over"
    );
    // A file that is no Parquet file and no tab output.
    let empty = directory.join("empty.parquet");
    std::fs::create_dir_all(&empty).expect("a directory is made");
    let (ended, _) = both_paths(&expected, &empty);
    assert_eq!(
        ended,
        Err(format!(
            "the Parquet adapter failed on {}: ERROR: [parquet_to_vep_tab] {}: no .parquet file under the directory",
            empty.display(),
            empty.display()
        ))
    );
}

/// The whole outputs are compared where the stream cannot read an output as its lines: a
/// reference or a candidate whose column line follows a data line, which the sorted form
/// reads; a row that holds a carriage return; a reference that is refused, in the
/// adapter's words.
#[test]
fn an_output_the_stream_cannot_read_line_by_line_is_compared_whole() {
    let directory = scratch("parquet_stream_whole");
    let rows = planted_rows(6);
    let out = planted_parts(&directory, "out.parquet", &[&rows]);
    let late = format!("{}\n{COLUMNS}\n{}\n", rows[0], rows[1..].join("\n"));
    let late = planted(&directory, "late.txt", &late);
    let (ended, note) = both_paths(&late, &out);
    assert!(ended.is_ok_and(|report| report.contains("\"verdict\": \"green\"")));
    assert_eq!(
        note,
        compared_whole("a reference that is not read as its canonical form")
    );
    let (ended, note) = both_paths(&out, &late);
    assert!(ended.is_ok_and(|report| report.contains("\"verdict\": \"green\"")));
    assert_eq!(
        note,
        compared_whole("a candidate that is not read as its canonical form")
    );
    let mut split = rows.clone();
    split[3] = split[3].replace("ENSG1", "EN\rSG1");
    let (ended, note) = both_paths(&planted(&directory, "split.txt", &tab_text(&split)), &out);
    assert!(ended.is_err_and(|error| error.starts_with("column count on line: ")));
    assert_eq!(note, compared_whole("a row that holds a line end"));
    let headless = format!("{}\n", rows.join("\n"));
    let reference = planted(&directory, "headless.txt", &headless);
    let (ended, note) = both_paths(&reference, &out);
    assert_eq!(
        ended,
        Err(format!(
            "the Parquet adapter failed on {}: ERROR: [parquet_to_vep_tab] {}: no # header line; a tab output needs its header",
            reference.display(),
            reference.display()
        ))
    );
    assert_eq!(
        note,
        compared_whole("a reference that is not read as its canonical form")
    );
}

/// A Parquet cell that holds a line end is one row of the sorted canonical form and two
/// lines of a text: the whole outputs are compared.
#[test]
fn a_parquet_cell_that_holds_a_line_end_is_compared_whole() {
    let directory = scratch("parquet_stream_cell");
    let names: Vec<&str> = COLUMNS[1..].split('\t').collect();
    let rows = planted_rows(3);
    let mut cells: Vec<Vec<String>> = rows
        .iter()
        .map(|row| row.split('\t').map(str::to_string).collect())
        .collect();
    cells[2][3] = "ENSG\n1".to_string();
    let columns: Vec<(&str, Stored)> = names
        .iter()
        .enumerate()
        .map(|(at, name)| {
            let values: Vec<Option<&str>> =
                cells.iter().map(|row| Some(row[at].as_str())).collect();
            (*name, Stored::texts(&values))
        })
        .collect();
    let out = directory.join("out.parquet");
    // Row groups of four rows: the cell stands in the first.
    plant(&out, &columns, 4, AS_THE_WRITER);
    let expected = planted(&directory, "tab.txt", &tab_text(&rows));
    let (ended, note) = both_paths(&expected, &out);
    assert!(ended.is_err_and(|error| error.starts_with("column count on line: ")));
    assert_eq!(note, compared_whole("a row that holds a line end"));
}

fn gzipped(text: &str) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder
        .write_all(text.as_bytes())
        .expect("the text is deflated");
    encoder.finish().expect("the member is ended")
}

/// The stream reads a gzip reference, a reference that is a Parquet output and a
/// candidate that is a tab output as it reads a tab reference and a Parquet candidate; a
/// reference of `\r\n` line ends is the same rows.
#[test]
fn either_output_is_a_tab_output_plain_or_gzip_or_a_parquet_output() {
    let directory = scratch("parquet_stream_kinds");
    let rows = planted_rows(8);
    let text = tab_text(&rows);
    let plain = planted(&directory, "tab.txt", &text);
    let gzip = directory.join("tab.txt.gz");
    std::fs::write(&gzip, gzipped(&text)).expect("the gzip reference is written");
    let returned = planted(&directory, "crlf.txt", &text.replace('\n', "\r\n"));
    let out = planted_parts(&directory, "out.parquet", &[&rows[6..], &rows[..6]]);
    for (expected, actual) in [
        (&gzip, &out),
        (&returned, &out),
        (&out, &out),
        (&out, &plain),
        (&plain, &gzip),
    ] {
        let (report, note) = through_the_stream(expected, actual);
        assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
        assert_eq!(number(&report, &["counts", "compared"]), 16);
        assert_eq!(note, kept_back(16, 8, 0), "{expected:?} {actual:?}");
    }
}

/// The text of an output: its column line, then its rows in the output's own order. Each
/// reading writes the same bytes, and its lines sorted are the sorted canonical form's.
#[test]
fn a_canonical_text_is_the_sorted_form_s_lines_in_the_output_s_order_at_each_reading() {
    let directory = scratch("parquet_stream_text");
    // More rows than one step of a reading prints.
    let rows = planted_rows(5000);
    let out = planted_parts(&directory, "out.parquet", &[&rows[9000..], &rows[..9000]]);
    let tab = planted(&directory, "tab.txt", &tab_text(&rows));
    let read = |text: &CanonicalText| {
        let mut bytes = String::new();
        Text::open(text)
            .expect("the text opens")
            .read_to_string(&mut bytes)
            .expect("the text is read");
        bytes
    };
    for (path, order) in [
        (&out, [&rows[9000..], &rows[..9000]].concat()),
        (&tab, rows.clone()),
    ] {
        let text = CanonicalText::of(path);
        let first = read(&text);
        assert_eq!(read(&text), first);
        assert_eq!(first, format!("{COLUMNS}\n{}\n", order.join("\n")));
        assert!(!text.unread() && !text.split());
        let mut sorted = order.clone();
        sorted.sort();
        assert_eq!(
            parquet_as_tab(path, true).expect("the sorted form is read"),
            format!("{COLUMNS}\n{}\n", sorted.join("\n"))
        );
    }
    let (report, note) = through_the_stream(&tab, &out);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(note, kept_back(10000, 5000, 0));
}

/// A candidate that stores a column named as an INT32 column as strings, as it stores a
/// plugin's field of that name: the stream reads the tab reference beside the candidate,
/// so the cells stay as written on both sides and the rows are kept back. A candidate that
/// stores the column as integers reads the same reference red through both paths.
#[test]
fn a_string_column_of_the_candidate_stays_as_written_in_the_text_of_the_reference() {
    let directory = scratch("parquet_stream_strings");
    let columns = "#Uploaded_variation\tLocation\tAllele\tFeature\tTSL";
    let cells = ["05", "1,5,2", "7", "-"];
    let rows: Vec<String> = cells
        .iter()
        .enumerate()
        .map(|(record, cell)| format!("rs{record}\t21:{}\tA\tENST1\t{cell}", 100 + record))
        .collect();
    let text = format!("{columns}\n{}\n", rows.join("\n"));
    let expected = planted(&directory, "tab.txt", &text);
    let stored = |at: usize| -> Vec<Option<&str>> {
        rows.iter()
            .map(|row| row.split('\t').nth(at).filter(|cell| *cell != "-"))
            .collect()
    };
    let names: Vec<&str> = columns[1..].split('\t').collect();
    let as_strings: Vec<(&str, Stored)> = names
        .iter()
        .enumerate()
        .map(|(at, name)| (*name, Stored::texts(&stored(at))))
        .collect();
    let out = directory.join("strings.parquet");
    plant(&out, &as_strings, 3, AS_THE_WRITER);
    let (report, note) = through_the_stream(&expected, &out);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(
        note,
        "4 of 4 reference rows were kept back by the stream of record groups, \
         0 of 4 record groups handed over"
    );
    let mut as_integers = as_strings.clone();
    as_integers[4] = ("TSL", Stored::Int32(vec![Some(5), None, Some(7), None]));
    let out = directory.join("integers.parquet");
    plant(&out, &as_integers, 3, AS_THE_WRITER);
    let (report, note) = through_the_stream(&expected, &out);
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "TSL", "differing"]), 1);
    assert_eq!(
        note,
        "3 of 4 reference rows were kept back by the stream of record groups, \
         1 of 4 record groups handed over"
    );
}

/// A row longer than the bytes one step of a reading takes and a last row without its line
/// end are each one row of the text, and two outputs of no row compare green.
#[test]
fn a_long_row_a_last_row_without_its_line_end_and_no_row_at_all_are_read_as_they_stand() {
    let directory = scratch("parquet_stream_edges");
    let mut rows = planted_rows(4);
    // A gene name of three million letters, past the megabyte a step reads.
    rows[5] = rows[5].replace("ENSG1", &format!("ENSG{}", "x".repeat(3_000_000)));
    let text = tab_text(&rows);
    let ended = planted(&directory, "ended.txt", &text);
    let unended = planted(&directory, "unended.txt", text.trim_end_matches('\n'));
    let out = planted_parts(&directory, "out.parquet", &[&rows[4..], &rows[..4]]);
    for (expected, actual) in [(&unended, &out), (&ended, &unended), (&unended, &ended)] {
        let (report, note) = through_the_stream(expected, actual);
        assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
        assert_eq!(number(&report, &["counts", "compared"]), 8);
        assert_eq!(note, kept_back(8, 4, 0), "{expected:?} {actual:?}");
    }
    let no_row = planted(&directory, "empty.txt", &tab_text(&[]));
    let empty = planted_parts(&directory, "empty.parquet", &[&[]]);
    let (report, note) = through_the_stream(&no_row, &empty);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(number(&report, &["counts", "compared"]), 0);
    assert_eq!(note, kept_back(0, 0, 0));
}

/// A row of a Parquet output that opens with `#` is a header line to a reader of lines,
/// and the sorted canonical form holds it among other rows than the output's own order
/// does: the whole outputs are compared, and both paths write one report.
#[test]
fn a_parquet_row_that_opens_with_a_number_sign_is_compared_whole() {
    let directory = scratch("parquet_stream_header_row");
    let names: Vec<&str> = COLUMNS[1..].split('\t').collect();
    let mut rows = planted_rows(3);
    // A row that stands before the `#` row in byte order and after it in the file.
    rows[0] = rows[0].replace("rec0", " rec0");
    rows.insert(0, "#h1\tc2\tc3\tc4\tc5\tc6\tc7\tc8".to_string());
    let cells: Vec<Vec<&str>> = rows.iter().map(|row| row.split('\t').collect()).collect();
    let columns: Vec<(&str, Stored)> = names
        .iter()
        .enumerate()
        .map(|(at, name)| {
            let values: Vec<Option<&str>> = cells.iter().map(|row| Some(row[at])).collect();
            (*name, Stored::texts(&values))
        })
        .collect();
    let out = directory.join("out.parquet");
    plant(&out, &columns, 4, AS_THE_WRITER);
    let (ended, note) = both_paths(&out, &out);
    assert!(ended.is_ok_and(|report| report.contains("\"verdict\": \"green\"")));
    assert_eq!(note, compared_whole("a row that opens with #"));
}
