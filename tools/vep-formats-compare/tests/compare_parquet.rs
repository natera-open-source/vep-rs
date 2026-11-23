// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A Parquet candidate compared with a tab reference: both sides are read through the
//! adapter in its canonical form. The Parquet files are planted with `pyarrow` where the
//! interpreter has it, under its own defaults (snappy pages), and with the `parquet`
//! crate's writer under the same codec where it has not.

mod compare_common;
mod parquet_writer;

use std::collections::BTreeSet;
use std::ffi::OsString;
use std::path::Path;
use std::process::Command;

use compare_common::*;
use indexmap::IndexMap;
use parquet_writer::{plant, plant_tab_rows_as, Codec, Encoding, Stored, AS_THE_WRITER};
use vep_formats_compare::compare_files::compare_files_through;
use vep_formats_compare::parse::parquet::{parquet_as_tab, tab_beside_parquet, text_columns};
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::Side;

/// Writes the rows of a tab output as a Parquet file of the writer's layout: the key
/// columns, then each tab column under its declared type (`DISTANCE` and `STRAND` integers,
/// `Consequence` and `FLAGS` lists, the three position columns structs, every other column
/// a string with `-` as NULL except `Allele`). Its arguments: the output, the tab output's
/// column line, a file of the rows.
const WRITE_PARQUET: &str = r#"
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

path = Path(sys.argv[1])
columns = sys.argv[2][1:].split("\t")
rows = [line for line in Path(sys.argv[3]).read_text(encoding="utf-8").split("\n") if line]
cells = [row.split("\t") for row in rows]


def span(text):
    if text == "-":
        return None
    positions, _, total = text.partition("/")
    start, dash, end = positions.partition("-")
    member = lambda m: None if m in ("?", "") else int(m)
    return {"start": member(start), "end": member(end) if dash else member(start), "total": int(total) if total else None}


def column(name, values):
    if name in ("DISTANCE", "STRAND"):
        return pa.array([None if v == "-" else int(v) for v in values], pa.int32())
    if name == "MAX_AF":
        return pa.array([None if v == "-" else float(v) for v in values], pa.float64())
    if name in ("AF", "gnomADe_AF"):
        return pa.array([None if v == "-" else [float(t) for t in v.split(",")] for v in values], pa.list_(pa.float64()))
    if name in ("Consequence", "FLAGS"):
        return pa.array([None if v == "-" else v.split(",") for v in values], pa.list_(pa.string()))
    if name in ("cDNA_position", "CDS_position", "Protein_position"):
        return pa.array([span(v) for v in values], pa.struct([("start", pa.int32()), ("end", pa.int32()), ("total", pa.int32())]))
    if name == "Allele":
        return pa.array(values, pa.string())
    return pa.array([None if v == "-" else v for v in values], pa.string())


n = len(cells)
arrays = [
    pa.array([0] * n, pa.int32()),
    pa.array(["21"] * n, pa.string()),
    pa.array(list(range(1, n + 1)), pa.int64()),
    pa.array(list(range(1, n + 1)), pa.int64()),
    pa.array(["A"] * n, pa.string()),
    pa.array(["G"] * n, pa.string()),
]
arrays.extend(column(name, [row[i] for row in cells]) for i, name in enumerate(columns))
table = pa.table(arrays, names=["contig_index", "chrom", "pos", "end", "ref", "alt", *columns])
path.parent.mkdir(parents=True, exist_ok=True)
pq.write_table(table.replace_schema_metadata({"vep_rs_version": "0"}), path)
"#;

fn interpreter() -> OsString {
    OsString::from("python3")
}

/// Whether the interpreter can write a Parquet file.
fn pyarrow_plants() -> bool {
    static PLANTS: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *PLANTS.get_or_init(|| {
        Command::new(interpreter())
            .args(["-c", "import pyarrow.parquet"])
            .output()
            .is_ok_and(|done| done.status.success())
    })
}

/// Plants `rows` of the tab output `tab_text` as the Parquet file `path`.
fn write_parquet(path: &Path, tab_text: &str, rows: &[String]) {
    let column_line = header_lines(tab_text)
        .pop()
        .expect("the tab output has its column line");
    if !pyarrow_plants() {
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        let as_pyarrow = Encoding {
            codec: Codec::Snappy,
            data_page_v2: false,
            dictionary: true,
        };
        plant_tab_rows_as(path, &column_line, &rows, as_pyarrow);
        return;
    }
    let held = path.with_extension("rows");
    if let Some(directory) = held.parent() {
        std::fs::create_dir_all(directory).expect("the directory is made");
    }
    std::fs::write(&held, rows.join("\n") + "\n").expect("the rows are written");
    let done = Command::new(interpreter())
        .args(["-c", WRITE_PARQUET])
        .arg(path)
        .arg(&column_line)
        .arg(&held)
        .output()
        .expect("the interpreter runs");
    assert!(
        done.status.success(),
        "{}",
        String::from_utf8_lossy(&done.stderr)
    );
}

/// The report of a tab reference against a Parquet candidate.
fn parquet_report(expected: &Path, actual: &Path) -> Value {
    compare_files_through(
        None,
        "parquet",
        expected,
        actual,
        None,
        &mut no_documented(),
        &BTreeSet::new(),
        &IndexMap::new(),
        &IndexMap::new(),
        4,
        false,
        None,
        None,
        None,
        &mut Side::default(),
    )
    .expect("the files compare")
}

/// Both sides go through the adapter, the tab reference for its canonical form: a Parquet
/// candidate holding the reference's rows in the reverse order compares green.
#[test]
fn a_parquet_candidate_goes_through_the_adapter_as_the_tab_reference_does() {
    let directory = scratch("parquet_reordered");
    let tab = text_of("tab");
    let mut rows = data_lines(&tab);
    rows.reverse();
    let out = directory.join("out.parquet");
    write_parquet(&out, &tab, &rows);
    let expected_file = planted(&directory, "tab.txt", &tab);
    let report = parquet_report(&expected_file, &out);
    assert_eq!(verdict(&report), "green");
    assert_eq!(at(&report, &["format"]).as_str(), Some("parquet"));
    assert_eq!(number(&report, &["counts", "compared"]), ROWS);
}

/// The adapter on a directory of two part files holding the reference's rows under the
/// writer's types: the comparison is green, and a row the Parquet side lacks is red.
#[test]
fn a_parquet_candidate_of_part_files_is_read_through_the_adapter() {
    let directory = scratch("parquet_parts");
    let tab = text_of("tab");
    let rows = data_lines(&tab);
    let half = rows.len() / 2;
    let out = directory.join("out.parquet");
    write_parquet(&out.join("part-00000.parquet"), &tab, &rows[..half]);
    write_parquet(&out.join("part-00001.parquet"), &tab, &rows[half..]);
    let expected_file = planted(&directory, "tab.txt", &tab);
    let as_tab = parquet_as_tab(&out, true).expect("the adapter reads the directory");
    let column_line = header_lines(&tab).pop().expect("a column line");
    assert_eq!(header_lines(&as_tab), [column_line]);
    let sorted = |mut lines: Vec<String>| {
        lines.sort();
        lines
    };
    assert_eq!(sorted(data_lines(&as_tab)), sorted(rows.clone()));
    let report = parquet_report(&expected_file, &out);
    assert_eq!(verdict(&report), "green");
    assert_eq!(number(&report, &["counts", "compared"]), ROWS);
    write_parquet(
        &out.join("part-00001.parquet"),
        &tab,
        &rows[half..rows.len() - 1],
    );
    let report = parquet_report(&expected_file, &out);
    assert_eq!(verdict(&report), "red");
    assert_eq!(
        at(&report, &["one_sided"]),
        &json(r#"{"expected_only": 1, "actual_only": 0}"#)
    );
    let refused =
        parquet_as_tab(&directory.join("absent.parquet"), true).expect_err("no such output");
    assert!(refused.to_string().contains("the Parquet adapter failed"));
}

/// Ensembl VEP's tab spells a frequency with an exponent or with trailing zeros; the
/// Parquet file stores the number, which the adapter prints in one decimal form. Both sides
/// are compared in the canonical form, so the spelling agrees and a different value still
/// differs.
#[test]
fn a_parquet_comparison_reads_frequencies_as_numbers() {
    let directory = scratch("parquet_frequencies");
    let header = "#Uploaded_variation\tLocation\tAllele\tConsequence\tAF\tgnomADe_AF\tMAX_AF";
    let tab_text = format!("## ENSEMBL VARIANT EFFECT PREDICTOR\n{header}\n");
    let row = "rs1\t21:100\tA\tintergenic_variant\t1.47e-05\t0.0010,0.2770\t0.2770";
    let expected_file = planted(&directory, "tab.txt", &format!("{tab_text}{row}\n"));
    let out = directory.join("out.parquet");
    write_parquet(&out, &tab_text, &[row.to_string()]);
    let canonical = |path: &Path| -> Vec<String> {
        data_lines(&parquet_as_tab(path, true).expect("the adapter reads the output"))
    };
    assert_eq!(
        canonical(&out),
        ["rs1\t21:100\tA\tintergenic_variant\t0.0000147\t0.001,0.277\t0.277"]
    );
    assert_eq!(canonical(&expected_file), canonical(&out));
    let report = parquet_report(&expected_file, &out);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(number(&report, &["counts", "compared"]), 1);
    write_parquet(&out, &tab_text, &[row.replace("\t0.2770", "\t0.2771")]);
    let report = parquet_report(&expected_file, &out);
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "MAX_AF", "differing"]), 1);
}

/// A plugin's field named as the flag's INT32 column holds a list in the tab reference, and
/// a `--custom` source's field named as a frequency column a number spelt with zeros at its
/// end. A candidate that stores both columns as strings holds each cell as written and
/// compares green. One that stores the first as INT32 holds NULL in the list's place, and
/// one whose strings spell another list, the integer without its zeros or the number
/// without its zeros differs: each is red on those cells.
#[test]
fn a_parquet_comparison_reads_a_string_column_as_stored_whatever_its_name() {
    let directory = scratch("parquet_string_columns");
    let header = "#Uploaded_variation\tLocation\tAllele\tConsequence\tTSL\tgnomADg_AF";
    let tab_text = format!("## ENSEMBL VARIANT EFFECT PREDICTOR\n{header}\n");
    let rows = [
        "rs1\t21:100\tA\tintergenic_variant\t1,5,2\t0.0010",
        "rs2\t21:101\tA\tintergenic_variant\t007\t-",
        "rs3\t21:102\tA\tintergenic_variant\t-\t1e-3",
    ];
    let expected_file = planted(
        &directory,
        "tab.txt",
        &format!("{tab_text}{}\n", rows.join("\n")),
    );
    let canonical = |path: &Path| -> Vec<String> {
        data_lines(&parquet_as_tab(path, true).expect("the adapter reads the output"))
    };
    let beside = |parquet: &Path| -> Vec<String> {
        data_lines(&tab_beside_parquet(&expected_file, parquet).expect("the tab output is read"))
    };
    let plant_with = |path: &Path, levels: Stored, frequencies: [Option<&str>; 3]| {
        let texts = |values: [&str; 3]| Stored::texts(&values.map(Some));
        let term: &[Option<&str>] = &[Some("intergenic_variant")];
        plant(
            path,
            &[
                ("Uploaded_variation", texts(["rs1", "rs2", "rs3"])),
                ("Location", texts(["21:100", "21:101", "21:102"])),
                ("Allele", texts(["A"; 3])),
                ("Consequence", Stored::text_lists(&[Some(term); 3])),
                ("TSL", levels),
                ("gnomADg_AF", Stored::texts(&frequencies)),
            ],
            usize::MAX,
            AS_THE_WRITER,
        );
    };
    let as_written = [Some("0.0010"), None, Some("1e-3")];

    // The part that stores the plugin's column as INT32: the list is NULL there.
    let lost = directory.join("lost.parquet");
    plant_with(&lost, Stored::Int32(vec![None, Some(7), None]), as_written);
    assert!(!text_columns(&lost).contains(&"TSL".to_string()));
    let report = parquet_report(&expected_file, &lost);
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "TSL", "differing"]), 1);
    assert_eq!(number(&report, &["columns", "gnomADg_AF", "differing"]), 0);
    assert_eq!(
        canonical(&lost),
        [
            "rs1\t21:100\tA\tintergenic_variant\t-\t0.0010",
            "rs2\t21:101\tA\tintergenic_variant\t7\t-",
            "rs3\t21:102\tA\tintergenic_variant\t-\t1e-3",
        ]
    );
    assert_eq!(
        beside(&lost),
        [
            "rs1\t21:100\tA\tintergenic_variant\t1,5,2\t0.0010",
            "rs2\t21:101\tA\tintergenic_variant\t7\t-",
            "rs3\t21:102\tA\tintergenic_variant\t-\t1e-3",
        ]
    );

    // The part that stores both columns as strings holds the rows as the tab output
    // writes them, and the tab output read beside it is those rows.
    let stored = directory.join("stored.parquet");
    plant_with(
        &stored,
        Stored::texts(&[Some("1,5,2"), Some("007"), None]),
        as_written,
    );
    assert_eq!(text_columns(&stored)[4..], ["TSL", "gnomADg_AF"]);
    assert_eq!(canonical(&stored), rows);
    assert_eq!(beside(&stored), rows);
    let report = parquet_report(&expected_file, &stored);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(number(&report, &["counts", "compared"]), 3);
    // A tab output read with no Parquet output beside it is reprinted by its names alone.
    assert_eq!(
        canonical(&expected_file),
        [
            "rs1\t21:100\tA\tintergenic_variant\t1,5,2\t0.001",
            "rs2\t21:101\tA\tintergenic_variant\t7\t-",
            "rs3\t21:102\tA\tintergenic_variant\t-\t0.001",
        ]
    );

    // Strings that spell another list, the integer and the number otherwise are other cells.
    let other = directory.join("other.parquet");
    plant_with(
        &other,
        Stored::texts(&[Some("1,5,3"), Some("7"), None]),
        [Some("0.001"), None, Some("1e-3")],
    );
    let report = parquet_report(&expected_file, &other);
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "TSL", "differing"]), 2);
    assert_eq!(number(&report, &["columns", "gnomADg_AF", "differing"]), 1);

    // A column is a string column of an output where every file of it stores strings.
    let mixed = directory.join("mixed.parquet");
    plant_with(
        &mixed.join("part-00000.parquet"),
        Stored::texts(&[Some("1,5,2"), Some("007"), None]),
        as_written,
    );
    plant_with(
        &mixed.join("part-00001.parquet"),
        Stored::Int32(vec![None, Some(7), None]),
        as_written,
    );
    assert_eq!(text_columns(&mixed)[4..], ["gnomADg_AF"]);
    assert!(text_columns(&directory.join("absent.parquet")).is_empty());
}

/// A candidate that is neither a regular file nor a directory names no string column and
/// is not opened: a named pipe that no writer holds open keeps a reader waiting.
#[test]
fn a_candidate_that_is_no_file_and_no_directory_is_not_opened_for_its_columns() {
    let directory = scratch("parquet_pipe");
    let pipe = directory.join("pipe.parquet");
    let made = Command::new("mkfifo")
        .arg(&pipe)
        .status()
        .expect("mkfifo runs");
    assert!(made.success());
    let (sender, receiver) = std::sync::mpsc::channel();
    std::thread::spawn(move || sender.send(text_columns(&pipe)));
    let named = receiver.recv_timeout(std::time::Duration::from_secs(30));
    assert_eq!(named, Ok(Vec::new()));
}
