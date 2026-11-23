// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary over files written under the writer's schema: the inverse of the cell
//! contract per type family, the header, every input shape, the canonical form, the
//! `--output` option and every refusal.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const NULL: &str = "-";
const KEY_COLUMNS: [&str; 6] = ["contig_index", "chrom", "pos", "end", "ref", "alt"];

/// One tab column of `cells.parquet` per type family and spelling, with the tab cells its
/// stored values print; a column shorter than the longest is NULL on the other rows.
const CELLS: &[(&str, &[&str])] = &[
    ("Feature", &["ENST00000456328", "", NULL]),
    ("Allele", &["-", "A"]),
    ("Location", &["21:25000100-25000101"]),
    ("STRAND", &["-1", "1", NULL]),
    ("DISTANCE", &["4821"]),
    ("MAX_AF", &["0.095", "1", "0.000098", "0", NULL]),
    ("MOTIF_SCORE_CHANGE", &["-0.012", "0.1"]),
    ("CANONICAL", &["YES", NULL, NULL]),
    ("PICK", &["1", NULL]),
    ("GENCODE_PRIMARY", &["1", NULL]),
    ("HIGH_INF_POS", &["Y", "N", NULL]),
    ("AF", &["0.0002", "0.2766,0.095", "-,0.1"]),
    ("gnomADe_AF", &["0.095,0.1"]),
    ("gnomADg_NFE_AF", &["0", NULL]),
    ("SOMATIC", &["0,1", NULL]),
    ("PHENO", &["1"]),
    ("MANE", &["MANE_Select,MANE_Plus_Clinical"]),
    (
        "Consequence",
        &[
            "missense_variant,splice_region_variant",
            "intergenic_variant",
        ],
    ),
    ("CLIN_SIG", &["benign;likely_benign,pathogenic"]),
    ("PUBMED", &["12345,67890", "-,1"]),
    ("EXON", &["3/12", "3-4/12"]),
    ("INTRON", &["7/11"]),
    (
        "cDNA_position",
        &["44", "44-45", "44/1200", "-/1200", "44-45/1200"],
    ),
    ("CDS_position", &["146-?"]),
    ("Protein_position", &["?-52", "49-?/300", NULL]),
    // A plugin column is a string; a BOOLEAN column of no listed spelling prints a flag.
    ("CADD_PHRED", &["23.4"]),
    ("FLAGGED", &["1", "0"]),
];
const ROWS: usize = 5;

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn fixture(name: &str) -> String {
    fixtures()
        .join(name)
        .to_str()
        .expect("a UTF-8 path")
        .to_string()
}

fn output_of(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vep-parquet-to-tab"))
        .args(args)
        .env("VEP_PARQUET_TO_TAB_THREADS", "2")
        .output()
        .expect("the binary runs")
}

/// The standard output of a run that succeeds.
fn run(args: &[&str]) -> String {
    let output = output_of(args);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty(), "{args:?}");
    String::from_utf8(output.stdout).expect("UTF-8 output")
}

/// The error line of a run that is refused: exit code 2 and the component's prefix.
fn fails(args: &[&str]) -> String {
    let output = output_of(args);
    assert_eq!(output.status.code(), Some(2), "{args:?}");
    let stderr = String::from_utf8(output.stderr).expect("a UTF-8 error line");
    assert!(
        stderr.starts_with("ERROR: [parquet_to_vep_tab] "),
        "{args:?}: {stderr}"
    );
    assert_eq!(stderr.matches('\n').count(), 1, "{args:?}: {stderr}");
    stderr
}

/// A directory of the test's own, removed when the test ends.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let path =
            std::env::temp_dir().join(format!("vep-parquet-to-tab-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).expect("a scratch directory");
        Self(path)
    }

    fn path(&self, name: &str) -> String {
        self.0
            .join(name)
            .to_str()
            .expect("a UTF-8 path")
            .to_string()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn every_typed_cell_prints_its_tab_text() {
    assert_eq!(CELLS.iter().map(|(_, texts)| texts.len()).max(), Some(ROWS));
    let stdout = run(&[&fixture("cells.parquet")]);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 1 + ROWS);
    for (row, line) in lines[1..].iter().enumerate() {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), CELLS.len());
        for ((column, texts), cell) in CELLS.iter().zip(cells) {
            let expected = texts.get(row).copied().unwrap_or(NULL);
            assert_eq!(cell, expected, "{column} row {row}");
        }
    }
}

#[test]
fn the_header_is_the_tab_columns_after_the_key_columns() {
    let stdout = run(&[&fixture("cells.parquet")]);
    let header = stdout.lines().next().expect("a header line");
    let names: Vec<&str> = CELLS.iter().map(|(name, _)| *name).collect();
    assert_eq!(header, format!("#{}", names.join("\t")));
    assert!(!header.split('\t').any(|cell| KEY_COLUMNS.contains(&cell)));
}

#[test]
fn the_parts_of_a_directory_come_out_in_path_order() {
    let stdout = run(&[&fixture("parts.parquet")]);
    assert_eq!(
        stdout.lines().collect::<Vec<_>>(),
        ["#Feature", "first", "second"]
    );
}

#[test]
fn a_hive_directory_is_read_file_by_file_with_chrom_from_the_files() {
    use arrow_array::cast::AsArray;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let stdout = run(&[&fixture("hive.parquet")]);
    assert_eq!(
        stdout.lines().collect::<Vec<_>>(),
        ["#Location", "21:100", "X:100", "X:101"]
    );
    // The stored `chrom` is what a reader of the files gets, whatever the directory says.
    let part = fs::File::open(fixtures().join("hive.parquet/chrom=X/part-00000.parquet"))
        .expect("the part");
    let reader = ParquetRecordBatchReaderBuilder::try_new(part)
        .expect("a footer")
        .build()
        .expect("a reader");
    let mut stored = Vec::new();
    for batch in reader {
        let batch = batch.expect("a batch");
        let chrom = batch.column_by_name("chrom").expect("the chrom column");
        stored.extend(
            chrom
                .as_string::<i32>()
                .iter()
                .map(|value| value.expect("a chromosome").to_string()),
        );
    }
    assert_eq!(stored, ["X", "X"]);
}

#[test]
fn canonical_sorts_the_rows_and_reprints_every_double_cell() {
    let scratch = Scratch::new("canonical");
    let header = "#Feature\tMAX_AF\tAF";
    let parquet = fixture("canonical/out.parquet");
    let canonical = run(&["--canonical", &parquet]);
    assert_eq!(
        canonical.lines().collect::<Vec<_>>(),
        [
            header,
            "a\t0.000098\t0.1",
            "b\t0.095\t0.2766,0.095",
            "c\t-\t-"
        ]
    );
    let plain = run(&[&parquet]);
    let mut rows: Vec<&str> = plain.lines().skip(1).collect();
    rows.sort_unstable();
    assert_eq!(rows, canonical.lines().skip(1).collect::<Vec<_>>());

    let tab = scratch.path("out.tab");
    let text = format!(
        "## ENSEMBL VARIANT EFFECT PREDICTOR\n## Output produced at ...\n\
         {header}\nc\t-\t-\nb\t0.0950\t0.2766,0.0950\na\t9.8e-05\t0.100\n"
    );
    fs::write(&tab, &text).expect("the tab output");
    assert_eq!(
        fs::read(fixtures().join("canonical/out.tab")).expect("the fixture"),
        text.as_bytes()
    );
    assert_eq!(run(&["--canonical", &tab]), canonical);
    let gz = scratch.path("out.tab.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        fs::File::create(&gz).expect("the gzip output"),
        flate2::Compression::default(),
    );
    std::io::Write::write_all(&mut encoder, text.as_bytes()).expect("the gzip data");
    encoder.finish().expect("the gzip trailer");
    assert_eq!(run(&["--canonical", &gz]), canonical);
    assert_eq!(
        run(&["--canonical", &fixture("canonical/out.tab.gz")]),
        canonical
    );

    let again = scratch.path("again.tab");
    fs::write(&again, &canonical).expect("the canonical form");
    assert_eq!(run(&["--canonical", &again]), canonical);

    // A line ends at a newline alone: a cell holding another separator character stays
    // one row, and a CRLF file reads like an LF one.
    let separators = scratch.path("separators.tab");
    fs::write(
        &separators,
        format!("{header}\r\nb\u{2028}c\t0.0950\t0.1\r\n"),
    )
    .expect("the tab output");
    assert_eq!(
        run(&["--canonical", &separators]),
        format!("{header}\nb\u{2028}c\t0.095\t0.1\n")
    );
}

#[test]
fn the_output_option_writes_a_file() {
    let scratch = Scratch::new("output");
    let parquet = fixture("cells.parquet");
    let out = scratch.path("out.tab");
    assert_eq!(run(&["--input", &parquet, "--output", &out]), "");
    assert_eq!(
        fs::read_to_string(&out).expect("the output"),
        run(&[&parquet])
    );
}

#[test]
fn errors_are_reported_on_stderr_with_exit_2() {
    let scratch = Scratch::new("errors");
    assert!(fails(&[&scratch.path("absent.parquet")]).contains("no such file"));
    assert!(fails(&["--canonical", &scratch.path("absent.tab")]).contains("no such file"));
    let empty = scratch.path("empty");
    fs::create_dir(&empty).expect("an empty directory");
    assert!(fails(&[&empty]).contains("no .parquet file"));
    assert!(fails(&[&fixture("errors/mixed.parquet")]).contains("tab columns differ"));
    assert!(fails(&["--canonical", &fixture("errors/headless.tab")]).contains("no # header line"));
    assert!(fails(&[&fixture("errors/keys.parquet")]).contains("no tab columns"));
    // A regular file without the Parquet magic is refused before its footer is read; a
    // file with the magic but no footer, a tab file that is not UTF-8 and a truncated gzip
    // get the ERROR line too.
    assert!(fails(&[&fixture("errors/plain.tab")]).contains("not a Parquet file"));
    let corrupt = fixture("errors/corrupt.parquet");
    fails(&[&corrupt]);
    fails(&["--canonical", &corrupt]);
    fails(&["--canonical", &fixture("errors/binary.tab")]);
    fails(&["--canonical", &fixture("errors/truncated.tab.gz")]);
    let fifo = scratch.path("fifo");
    let name = std::ffi::CString::new(fifo.as_str()).expect("a path without a NUL");
    // SAFETY: `name` is a NUL-terminated path that lives through the call.
    assert_eq!(unsafe { libc::mkfifo(name.as_ptr(), 0o600) }, 0);
    assert!(fails(&[&fifo]).contains("neither a regular file nor a directory"));
    assert!(fails(&["--canonical", &fifo]).contains("neither a regular file nor a directory"));
}

#[test]
fn a_refused_run_leaves_the_named_output_empty() {
    let scratch = Scratch::new("refused");
    let out = scratch.path("out.tab");
    fails(&[&scratch.path("absent.parquet"), "--output", &out]);
    assert_eq!(fs::read(&out).expect("the output"), b"");
    // A token no number spells ends the run after the header line.
    let header_only = output_of(&[
        "--canonical",
        &fixture("errors/unconvertible.tab"),
        "--output",
        &out,
    ]);
    assert_eq!(header_only.status.code(), Some(2));
    assert_eq!(fs::read(&out).expect("the output"), b"#Feature\tAF\n");
    // A `--types-of` that names no Parquet output ends the run before the header line, and
    // a refusal of the tab input itself comes first.
    let absent = scratch.path("absent.parquet");
    let tab = fixture("canonical/out.tab");
    fails(&["--canonical", "--types-of", &absent, &tab, "--output", &out]);
    assert_eq!(fs::read(&out).expect("the output"), b"");
    let headless = fixture("errors/headless.tab");
    let refused = fails(&["--canonical", "--types-of", &absent, &headless]);
    assert!(refused.contains("no # header line"), "{refused}");
}

#[test]
fn a_file_that_cannot_be_opened_is_named_with_its_errno() {
    use std::os::unix::fs::PermissionsExt;

    // SAFETY: `geteuid` reads the process's own identity.
    if unsafe { libc::geteuid() } == 0 {
        eprintln!("skipped: the superuser opens a file of any mode");
        return;
    }
    let scratch = Scratch::new("unreadable");
    let tab = scratch.path("locked.tab");
    fs::write(&tab, "#Feature\na\n").expect("the tab output");
    fs::set_permissions(&tab, fs::Permissions::from_mode(0o000)).expect("a mode");
    assert_eq!(
        fails(&["--canonical", &tab]),
        format!("ERROR: [parquet_to_vep_tab] [Errno 13] Permission denied: '{tab}'\n")
    );
    let parts = scratch.path("parts");
    fs::create_dir(&parts).expect("a directory");
    let part = format!("{parts}/part-00000.parquet");
    fs::copy(fixtures().join("parts.parquet/part-00000.parquet"), &part).expect("a part");
    fs::set_permissions(&part, fs::Permissions::from_mode(0o000)).expect("a mode");
    assert_eq!(
        fails(&[&parts]),
        format!(
            "ERROR: [parquet_to_vep_tab] [Errno 13] Failed to open local file '{part}'. \
             Detail: [errno 13] Permission denied\n"
        )
    );
}

#[test]
fn a_decode_error_past_the_first_block_names_its_offset_in_the_file() {
    let scratch = Scratch::new("offset");
    let tab = scratch.path("long.tab");
    let line = "row\t0.5\n";
    let rows = (20 << 20) / line.len();
    let mut text = String::from("#Feature\tMAX_AF\n").into_bytes();
    text.extend(line.repeat(rows).bytes());
    let offset = text.len();
    text.extend_from_slice(b"bad\xff\n");
    fs::write(&tab, &text).expect("the tab output");
    assert_eq!(
        fails(&["--canonical", &tab]),
        format!(
            "ERROR: [parquet_to_vep_tab] 'utf-8' codec can't decode byte 0xff in position {}: \
             invalid start byte\n",
            offset + 3
        )
    );
    // The same rows without the broken line come out sorted, one block or many.
    text.truncate(offset);
    fs::write(&tab, &text).expect("the tab output");
    let canonical = run(&["--canonical", &tab]);
    assert_eq!(canonical.len(), text.len());
    assert!(canonical.starts_with("#Feature\tMAX_AF\nrow\t0.5\n"));
}

#[test]
fn a_reader_that_closes_the_pipe_ends_the_run_quietly() {
    use std::process::Stdio;

    let mut child = Command::new(env!("CARGO_BIN_EXE_vep-parquet-to-tab"))
        .arg(fixture("cells.parquet"))
        .env("VEP_PARQUET_TO_TAB_THREADS", "2")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("the binary runs");
    drop(child.stdout.take());
    let output = child.wait_with_output().expect("the binary ends");
    assert_eq!(output.status.code(), Some(0));
    assert!(output.stderr.is_empty());
}

#[test]
fn a_row_group_that_cannot_be_read_ends_the_run_after_the_ones_before_it() {
    use parquet::file::reader::{FileReader, SerializedFileReader};

    let scratch = Scratch::new("pages");
    let source = fixtures().join("groups.parquet/part-00000.parquet");
    let footer =
        SerializedFileReader::new(fs::File::open(&source).expect("the part")).expect("a footer");
    let groups = footer.metadata().row_groups();
    assert_eq!(groups.len(), 7);
    // The first page of the fifth row group's tab column, its header overwritten.
    let column = groups[4].columns().last().expect("the tab column");
    let (start, length) = column.byte_range();
    assert!(length > 12);
    let mut bytes = fs::read(&source).expect("the part's bytes");
    for byte in &mut bytes[start as usize + 2..start as usize + 12] {
        *byte ^= 0xff;
    }
    let broken = scratch.path("broken.parquet");
    fs::write(&broken, &bytes).expect("the broken part");
    let output = output_of(&[&broken]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        "#Feature\nrow0\nrow1\nrow2\nrow3\n"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    let error = stderr.lines().last().expect("an error line");
    assert!(
        error.starts_with(&format!("ERROR: [parquet_to_vep_tab] {broken}: ")),
        "{stderr}"
    );
    // The canonical form is written whole or not at all past its header.
    let output = output_of(&["--canonical", &broken]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(String::from_utf8_lossy(&output.stdout), "#Feature\n");
}

/// Writes a part of the writer's layout: the six key columns, then `columns`, a row an item.
fn write_part(path: &str, columns: Vec<(&str, arrow_array::ArrayRef)>) {
    write_part_stating(path, columns, None);
}

/// Writes a part as [`write_part`] does, its footer stating `tab_columns` under the key of
/// that name where one is given.
fn write_part_stating(
    path: &str,
    columns: Vec<(&str, arrow_array::ArrayRef)>,
    tab_columns: Option<&str>,
) {
    use arrow_array::{ArrayRef, Int32Array, Int64Array, RecordBatch, StringArray};
    use parquet::arrow::ArrowWriter;
    use parquet::file::metadata::KeyValue;
    use parquet::file::properties::WriterProperties;
    use std::sync::Arc;

    let rows = columns.first().map_or(0, |(_, values)| values.len());
    let places: Vec<i64> = (1..=rows as i64).collect();
    let mut fields: Vec<(&str, ArrayRef, bool)> = vec![
        (
            "contig_index",
            Arc::new(Int32Array::from(vec![0; rows])),
            true,
        ),
        ("chrom", Arc::new(StringArray::from(vec!["21"; rows])), true),
        ("pos", Arc::new(Int64Array::from(places.clone())), true),
        ("end", Arc::new(Int64Array::from(places)), true),
        ("ref", Arc::new(StringArray::from(vec!["A"; rows])), true),
        ("alt", Arc::new(StringArray::from(vec!["G"; rows])), true),
    ];
    fields.extend(
        columns
            .into_iter()
            .map(|(name, values)| (name, values, true)),
    );
    let batch = RecordBatch::try_from_iter_with_nullable(fields).expect("the rows of the part");
    let file = fs::File::create(path).expect("the part's file");
    let stated = tab_columns.map(|names| {
        let entry = KeyValue::new("tab_columns".to_string(), names.to_string());
        WriterProperties::builder()
            .set_key_value_metadata(Some(vec![entry]))
            .build()
    });
    let mut writer = ArrowWriter::try_new(file, batch.schema(), stated).expect("a Parquet writer");
    writer.write(&batch).expect("the rows are written");
    writer.close().expect("the footer is written");
}

/// A plugin's field named as the flag's INT32 column holds a list, and a `--custom`
/// source's field named as a frequency column a number spelt with zeros at its end. A part
/// that stores both columns as strings holds each cell as the `--tab` output writes it: its
/// canonical form is those rows, and so is the tab output's under `--types-of`. A part that
/// stores the first column as INT32 holds NULL in the list's place, and its canonical form
/// differs from the tab output's on that cell.
#[test]
fn the_canonical_form_keeps_a_string_column_as_stored_whatever_its_name() {
    use arrow_array::{ArrayRef, Int32Array, StringArray};
    use std::sync::Arc;

    let scratch = Scratch::new("string-columns");
    let header = "#Feature\tTSL\tgnomADg_AF";
    let features = || -> ArrayRef { Arc::new(StringArray::from(vec!["a", "b", "c"])) };
    let frequencies = || -> ArrayRef {
        Arc::new(StringArray::from(vec![
            Some("0.0010,x"),
            None,
            Some("1e-3"),
        ]))
    };
    let stored = scratch.path("stored.parquet");
    write_part(
        &stored,
        vec![
            ("Feature", features()),
            (
                "TSL",
                Arc::new(StringArray::from(vec![Some("1,5,2"), Some("007"), None])),
            ),
            ("gnomADg_AF", frequencies()),
        ],
    );
    let lost = scratch.path("lost.parquet");
    write_part(
        &lost,
        vec![
            ("Feature", features()),
            ("TSL", Arc::new(Int32Array::from(vec![None, Some(7), None]))),
            ("gnomADg_AF", frequencies()),
        ],
    );
    let tab = scratch.path("out.tab");
    fs::write(
        &tab,
        format!("{header}\nc\t-\t1e-3\nb\t007\t-\na\t1,5,2\t0.0010,x\n"),
    )
    .expect("the tab output");

    // The part that stores strings: its rows, its canonical form and the tab output read
    // under its types are one text.
    let as_written = format!("{header}\na\t1,5,2\t0.0010,x\nb\t007\t-\nc\t-\t1e-3\n");
    assert_eq!(run(&[&stored]), as_written);
    assert_eq!(run(&["--canonical", &stored]), as_written);
    assert_eq!(
        run(&["--canonical", "--types-of", &stored, &tab]),
        as_written
    );
    // The part that stores INT32: the list is NULL there, and the tab output read under
    // its types keeps the list.
    assert_eq!(
        run(&["--canonical", &lost]),
        format!("{header}\na\t-\t0.0010,x\nb\t7\t-\nc\t-\t1e-3\n")
    );
    assert_eq!(
        run(&["--canonical", "--types-of", &lost, &tab]),
        format!("{header}\na\t1,5,2\t0.0010,x\nb\t7\t-\nc\t-\t1e-3\n")
    );
    // A tab output read alone is reprinted by the names of its columns.
    assert_eq!(
        run(&["--canonical", &tab]),
        format!("{header}\na\t1,5,2\t0.001,x\nb\t7\t-\nc\t-\t0.001\n")
    );

    // `--types-of` names a Parquet output, and only a tab input under `--canonical` reads it.
    let absent = scratch.path("absent.parquet");
    assert_eq!(
        fails(&["--canonical", "--types-of", &absent, &tab]),
        format!("ERROR: [parquet_to_vep_tab] {absent}: no such file or directory\n")
    );
    assert_eq!(
        fails(&["--canonical", "--types-of", &tab, &tab]),
        format!("ERROR: [parquet_to_vep_tab] {tab}: not a Parquet file\n")
    );
    assert_eq!(
        run(&["--canonical", "--types-of", &absent, &stored]),
        as_written
    );
    assert_eq!(run(&["--types-of", &absent, &stored]), as_written);
}

/// A cell is a string where the field that prints it stores strings: a list of strings
/// counts, and each file is read under its own types. A stored tab still splits the line,
/// and `--types-of` takes a name only where every file stores it as strings.
#[test]
fn a_cell_is_a_string_where_the_field_that_prints_it_stores_strings() {
    use arrow_array::builder::{ListBuilder, StringBuilder};
    use arrow_array::{ArrayRef, Int32Array, StringArray};
    use std::sync::Arc;

    let scratch = Scratch::new("printing-field");
    let text = |cell: &str| -> ArrayRef { Arc::new(StringArray::from(vec![cell])) };
    let levels = |level: Option<i32>| -> ArrayRef { Arc::new(Int32Array::from(vec![level])) };
    let tab_of = |name: &str, rows: &str| {
        let tab = scratch.path(name);
        fs::write(&tab, rows).expect("the tab output");
        tab
    };

    // A list of strings under a frequency column's name stays as stored.
    let mut tokens = ListBuilder::new(StringBuilder::new());
    tokens.append_value([Some("0.0010"), Some("1e-3")]);
    let listed = scratch.path("listed.parquet");
    write_part(
        &listed,
        vec![
            ("Feature", text("a")),
            ("gnomADg_AF", Arc::new(tokens.finish())),
        ],
    );
    let rows = "#Feature\tgnomADg_AF\na\t0.0010,1e-3\n";
    assert_eq!(run(&["--canonical", &listed]), rows);
    let tab = tab_of("listed.tab", rows);
    assert_eq!(run(&["--canonical", "--types-of", &listed, &tab]), rows);

    // A stored tab splits the line, and each cell is read by the column at its place: the
    // string that stands under the INT32 column is reprinted, on both sides alike.
    let split = scratch.path("split.parquet");
    write_part(
        &split,
        vec![
            ("Feature", text("a\tb")),
            ("SYMBOL", text("007")),
            ("TSL", levels(None)),
        ],
    );
    let header = "#Feature\tSYMBOL\tTSL";
    let by_place = format!("{header}\na\tb\t7\t-\n");
    assert_eq!(run(&["--canonical", &split]), by_place);
    let tab = tab_of("split.tab", &format!("{header}\na\tb\t007\t-\n"));
    assert_eq!(run(&["--canonical", "--types-of", &split, &tab]), by_place);

    // Two parts store one column under two types: the cells of each part follow its own
    // field, and `--types-of` reprints the column by its name.
    let parts = scratch.path("parts.parquet");
    fs::create_dir(&parts).expect("a directory");
    write_part(
        &format!("{parts}/part-00000.parquet"),
        vec![("Feature", text("a")), ("TSL", levels(Some(7)))],
    );
    let strings = format!("{parts}/part-00001.parquet");
    write_part(&strings, vec![("Feature", text("a")), ("TSL", text("007"))]);
    assert_eq!(
        run(&["--canonical", &parts]),
        "#Feature\tTSL\na\t007\na\t7\n"
    );
    let tab = tab_of("parts.tab", "#Feature\tTSL\na\t007\n");
    assert_eq!(
        run(&["--canonical", "--types-of", &parts, &tab]),
        "#Feature\tTSL\na\t7\n"
    );
    assert_eq!(
        run(&["--canonical", "--types-of", &strings, &tab]),
        "#Feature\tTSL\na\t007\n"
    );
}

/// A tab column is read from its own field, the one at its place after the leading key
/// columns: a plugin's `ref`, `alt` and `pos` print the plugin's cells, each later cell
/// stands under its own column and the row keeps its last cell. The canonical form of the
/// part is then the canonical form of the tab output that holds the same rows.
#[test]
fn a_tab_column_named_as_a_key_column_prints_its_own_field() {
    use arrow_array::{ArrayRef, BooleanArray, Float64Array, Int32Array, StringArray};
    use std::sync::Arc;

    let scratch = Scratch::new("own-field");
    let texts = |cells: Vec<Option<&str>>| -> ArrayRef { Arc::new(StringArray::from(cells)) };
    let part = scratch.path("named.parquet");
    write_part(
        &part,
        vec![
            ("Feature", texts(vec![Some("a"), Some("b")])),
            ("ref", texts(vec![Some("007"), None])),
            (
                "MAX_AF",
                Arc::new(Float64Array::from(vec![Some(0.001), None])),
            ),
            ("alt", texts(vec![None, Some("plugin-alt")])),
            (
                "CANONICAL",
                Arc::new(BooleanArray::from(vec![Some(true), None])),
            ),
            ("pos", texts(vec![Some("12"), Some("plugin-pos")])),
            ("TSL", Arc::new(Int32Array::from(vec![Some(5), None]))),
        ],
    );
    let header = "#Feature\tref\tMAX_AF\talt\tCANONICAL\tpos\tTSL";
    let rows =
        format!("{header}\na\t007\t0.001\t-\tYES\t12\t5\nb\t-\t-\tplugin-alt\t-\tplugin-pos\t-\n");
    assert_eq!(run(&[&part]), rows);
    assert_eq!(run(&["--canonical", &part]), rows);

    // The tab output of the same rows, in another order and with its numbers as a run
    // writes them.
    let tab = scratch.path("named.tab");
    fs::write(
        &tab,
        format!(
            "{header}\nb\t-\t-\tplugin-alt\t-\tplugin-pos\t-\na\t007\t0.0010\t-\tYES\t12\t05\n"
        ),
    )
    .expect("the tab output");
    assert_eq!(run(&["--canonical", "--types-of", &part, &tab]), rows);
}

/// Tab columns of one name are fields of their own: each prints the cells of the field at
/// its place, and the canonical form keeps the cell of the one that stores strings as stored.
#[test]
fn tab_columns_of_one_name_print_each_its_own_field() {
    use arrow_array::{ArrayRef, Int32Array, StringArray};
    use std::sync::Arc;

    let scratch = Scratch::new("one-name");
    let text = |cell: &str| -> ArrayRef { Arc::new(StringArray::from(vec![cell])) };
    let levels = |level: i32| -> ArrayRef { Arc::new(Int32Array::from(vec![level])) };
    let part = scratch.path("twice.parquet");
    write_part(
        &part,
        vec![
            ("Feature", text("a")),
            ("TSL", levels(5)),
            ("TSL", text("007")),
            ("TSL", levels(6)),
        ],
    );
    let rows = "#Feature\tTSL\tTSL\tTSL\na\t5\t007\t6\n";
    assert_eq!(run(&[&part]), rows);
    assert_eq!(run(&["--canonical", &part]), rows);
}

/// The writer stores a tab column whose name an earlier column carries in other letters
/// behind `tab_` and states the tab columns' names in the footer: the header prints the
/// stated names, and a BOOLEAN column prints the spelling of its stated name. A file that
/// states no such list, or a list of another length, names each column by its field, one
/// prefix before a key column's name taken off.
#[test]
fn a_column_stored_behind_the_prefix_for_its_letter_case_prints_under_the_stated_name() {
    use arrow_array::{ArrayRef, BooleanArray, StringArray};
    use std::sync::Arc;

    let scratch = Scratch::new("stated-names");
    let text = |cell: &str| -> ArrayRef { Arc::new(StringArray::from(vec![cell])) };
    let columns = || -> Vec<(&str, ArrayRef)> {
        vec![
            ("Feature", text("a")),
            ("tab_REF", text("source")),
            ("canonical", text("no")),
            ("tab_CANONICAL", Arc::new(BooleanArray::from(vec![true]))),
            ("tab_alt", text("plugin")),
        ]
    };
    let part = scratch.path("stated.parquet");
    let stated = "Feature\tREF\tcanonical\tCANONICAL\talt";
    write_part_stating(&part, columns(), Some(stated));
    let rows = "#Feature\tREF\tcanonical\tCANONICAL\talt\na\tsource\tno\tYES\tplugin\n";
    assert_eq!(run(&[&part]), rows);
    assert_eq!(run(&["--canonical", &part]), rows);
    let tab = scratch.path("stated.tab");
    fs::write(&tab, rows).expect("the tab output");
    assert_eq!(run(&["--canonical", "--types-of", &part, &tab]), rows);

    let by_field = "#Feature\ttab_REF\tcanonical\ttab_CANONICAL\talt\na\tsource\tno\t1\tplugin\n";
    for other in [
        None,
        Some("Feature\tREF\tcanonical\tCANONICAL"),
        Some("Feature\tREF\tcanonical\tCANONICAL\talt\tTSL"),
    ] {
        let part = scratch.path("unstated.parquet");
        write_part_stating(&part, columns(), other);
        assert_eq!(run(&[&part]), by_field, "{other:?}");
    }
}

/// The README says which names the writer stores behind the prefix and where the header
/// takes a tab column's name.
#[test]
fn the_readme_says_the_header_names_a_column_as_the_footer_states_it() {
    let readme = fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("README.md"))
        .expect("the README");
    let text = readme.split_whitespace().collect::<Vec<_>>().join(" ");
    for sentence in [
        "It stores a tab column whose name an earlier column carries in other letters behind \
         `tab_` as well (`tab_REF` for a custom source named `REF`), and states the `--tab` \
         name of every tab column in the footer, under the key `tab_columns`, joined by tabs.",
        "The header line names each column as the footer states it, and the spelling of a \
         BOOLEAN column and the canonical form follow that name.",
        "For a file that states no such list, or a list of another length, a stored name that \
         is a key column's behind one or more `tab_` prints without the first, and every \
         other name prints as stored.",
    ] {
        assert!(text.contains(sentence), "the README lacks: {sentence}");
    }
    assert!(
        !text.contains("The header line names each column as the `--tab` format does:"),
        "the README names the columns by their stored names alone"
    );
}

/// The README says that a tab column is read from its own field, by its place after the
/// leading key columns, and holds no sentence that reads a column by its name.
#[test]
fn the_readme_says_a_tab_column_is_read_from_its_own_field() {
    let readme = fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("README.md"))
        .expect("the README");
    let text = readme.split_whitespace().collect::<Vec<_>>().join(" ");
    for sentence in [
        "**Columns by place.** A tab column is read from its own field: the field at its place \
         after the leading key columns, whatever its name.",
        "A tab column named as a key column the file leads with (`ref` or `alt` among a \
         plugin's own columns) prints the cells of its own field, and the row holds a cell for \
         every column of the header.",
        "the tab columns by their place and the row groups.",
    ] {
        assert!(text.contains(sentence), "the README lacks: {sentence}");
    }
    for replaced in [
        "Columns by name",
        "selects the key column too",
        "stands one column on",
        "the columns read by name",
    ] {
        assert!(!text.contains(replaced), "the README holds: {replaced}");
    }
}

/// The README's items on the canonical form say what the INT32 rule does with a cell that
/// is no INT32, that a string column is no number column whatever its name, and how a
/// `--tab` output is read under the types of a Parquet output.
#[test]
fn the_readme_says_how_the_canonical_form_reads_a_string_column() {
    let readme = fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("README.md"))
        .expect("the README");
    let text = readme.split_whitespace().collect::<Vec<_>>().join(" ");
    for sentence in [
        "a cell that is no INT32, an empty one included, stays as written.",
        "A string is no number, whatever its column's name: a field that a plugin or a \
         `--custom` source declares under one of these names is a string column of a Parquet \
         output.",
        "A cell of a Parquet input that is printed from a string column stays as stored.",
        "A `--tab` input does not state the types of its columns, so each is reprinted by its \
         name.",
        "With `--types-of PARQUET` a column stays as written where every file of that Parquet \
         output stores every column of its name as strings.",
        "A `--tab` cell that an INT32 column of a Parquet output could not hold is NULL there \
         and prints `-`, so the two forms differ on that cell.",
        "<($vep --canonical --types-of out.parquet out.tab)",
    ] {
        assert!(text.contains(sentence), "the README lacks: {sentence}");
    }
    assert!(!text.contains("as the writer stores it"));
}
