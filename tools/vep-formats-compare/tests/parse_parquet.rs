// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A Parquet output read as tab text: every typed cell, the layouts of an output, the
//! canonical form, the refusals and the Parquet cases of the command line, on files planted
//! here under the writer's schema, each held to the text stated for it.
//!
//! Two ignored tests read the files of cases under the `parquet` folder of the directory
//! `VEP_FORMATS_COMPARE_PYDIFF_DIR` names: `canonical<TAB>path<TAB>answer` lines written
//! by CPython over files it planted, the path as the lowercase hex of its bytes, the
//! answer `K` and the hex of the text or `E` and the hex of the comparator's error. One
//! runs the candidate engine `VEP_BINARY` names over three golden corpora, and one times
//! the reading of the output `VEP_FORMATS_COMPARE_PARQUET` names.

mod parquet_writer;

use flate2::write::GzEncoder;
use parquet::basic::{LogicalType, Repetition, Type as PhysicalType};
use parquet::file::reader::FileReader;
use parquet::file::serialized_reader::SerializedFileReader;
use parquet_writer::{
    group_field, plant, plant_rows, plant_tab_rows, text_field, value_field, write_file, Codec,
    Encoding, LeafData, Stored, Values, AS_THE_WRITER,
};
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use vep_formats_compare::error::Error;
use vep_formats_compare::input::{read_text, text_lines};
use vep_formats_compare::parse::parquet::{parquet_as_tab, tab_beside_parquet, text_columns};

const NULL: &str = "-";
const FAILED: &str = "the Parquet adapter failed on ";
const ADAPTER: &str = "ERROR: [parquet_to_vep_tab] ";

/// An empty directory of this test's own.
fn scratch(name: &str) -> PathBuf {
    let dir =
        Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("vep-formats-compare-parquet-{name}"));
    if dir.exists() {
        std::fs::remove_dir_all(&dir).expect("the old scratch directory is removed");
    }
    std::fs::create_dir_all(&dir).expect("the scratch directory is made");
    dir
}

fn tab(path: &Path, canonical: bool) -> String {
    parquet_as_tab(path, canonical).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The comparator's error on `path`, without the words that open it and name the path.
fn refusal(path: &Path, canonical: bool) -> String {
    let Err(Error::Comparator(text)) = parquet_as_tab(path, canonical) else {
        panic!("{} is read", path.display());
    };
    let opening = format!("{FAILED}{}: ", path.display());
    text.strip_prefix(&opening)
        .unwrap_or_else(|| panic!("{text} does not open with {opening}"))
        .to_string()
}

fn lines(text: &str) -> Vec<&str> {
    text.strip_suffix('\n')
        .expect("the text ends its last line")
        .split('\n')
        .collect()
}

fn data_lines(text: &str) -> Vec<String> {
    text_lines(text)
        .into_iter()
        .filter(|line| !line.starts_with('#'))
        .collect()
}

fn sorted(mut rows: Vec<String>) -> Vec<String> {
    rows.sort();
    rows
}

fn texts(values: &[Option<&str>]) -> Stored {
    Stored::texts(values)
}

fn doubles(values: &[Option<f64>]) -> Stored {
    Stored::Double(values.to_vec())
}

fn double_lists(values: &[Option<&[Option<f64>]>]) -> Stored {
    Stored::DoubleList(values.iter().map(|list| list.map(<[_]>::to_vec)).collect())
}

fn bool_lists(values: &[Option<&[Option<bool>]>]) -> Stored {
    Stored::BoolList(values.iter().map(|list| list.map(<[_]>::to_vec)).collect())
}

fn text_lists(values: &[Option<&[Option<&str>]>]) -> Stored {
    Stored::text_lists(values)
}

/// The `start`, `end` and `total` of a positions struct.
type Span = (Option<i32>, Option<i32>, Option<i32>);

fn positions(values: &[Option<Span>]) -> Stored {
    Stored::Positions(
        values
            .iter()
            .map(|members| members.map(|(start, end, total)| [start, end, total]))
            .collect(),
    )
}

/// One tab column per type family and spelling: the column, its stored values and the tab
/// cells those values print. A column shorter than the longest is NULL on the other rows.
fn cells() -> Vec<(&'static str, Stored, Vec<&'static str>)> {
    let span = |start, end, total| Some((start, end, total));
    vec![
        (
            "Feature",
            texts(&[Some("ENST00000456328"), Some(""), None]),
            vec!["ENST00000456328", "", NULL],
        ),
        ("Allele", texts(&[Some("-"), Some("A")]), vec!["-", "A"]),
        (
            "Location",
            texts(&[Some("21:25000100-25000101")]),
            vec!["21:25000100-25000101"],
        ),
        (
            "STRAND",
            Stored::Int32(vec![Some(-1), Some(1), None]),
            vec!["-1", "1", NULL],
        ),
        ("DISTANCE", Stored::Int32(vec![Some(4821)]), vec!["4821"]),
        (
            "MAX_AF",
            doubles(&[Some(0.095), Some(1.0), Some(9.8e-05), Some(0.0), None]),
            vec!["0.095", "1", "0.000098", "0", NULL],
        ),
        (
            "MOTIF_SCORE_CHANGE",
            doubles(&[Some(-0.012), Some(0.1)]),
            vec!["-0.012", "0.1"],
        ),
        (
            "CANONICAL",
            Stored::Bool(vec![Some(true), Some(false), None]),
            vec!["YES", NULL, NULL],
        ),
        (
            "PICK",
            Stored::Bool(vec![Some(true), Some(false)]),
            vec!["1", NULL],
        ),
        (
            "GENCODE_PRIMARY",
            Stored::Bool(vec![Some(true), Some(false)]),
            vec!["1", NULL],
        ),
        (
            "HIGH_INF_POS",
            Stored::Bool(vec![Some(true), Some(false), None]),
            vec!["Y", "N", NULL],
        ),
        (
            "AF",
            double_lists(&[
                Some(&[Some(0.0002)]),
                Some(&[Some(0.2766), Some(0.095)]),
                Some(&[None, Some(0.1)]),
            ]),
            vec!["0.0002", "0.2766,0.095", "-,0.1"],
        ),
        (
            "gnomADe_AF",
            double_lists(&[Some(&[Some(0.095), Some(0.1)])]),
            vec!["0.095,0.1"],
        ),
        (
            "gnomADg_NFE_AF",
            double_lists(&[Some(&[Some(0.0)]), None]),
            vec!["0", NULL],
        ),
        (
            "SOMATIC",
            bool_lists(&[Some(&[Some(false), Some(true)]), Some(&[None])]),
            vec!["0,1", NULL],
        ),
        ("PHENO", bool_lists(&[Some(&[Some(true)])]), vec!["1"]),
        (
            "MANE",
            text_lists(&[Some(&[Some("MANE_Select"), Some("MANE_Plus_Clinical")])]),
            vec!["MANE_Select,MANE_Plus_Clinical"],
        ),
        (
            "Consequence",
            text_lists(&[
                Some(&[Some("missense_variant"), Some("splice_region_variant")]),
                Some(&[Some("intergenic_variant")]),
            ]),
            vec![
                "missense_variant,splice_region_variant",
                "intergenic_variant",
            ],
        ),
        (
            "CLIN_SIG",
            text_lists(&[Some(&[Some("benign;likely_benign"), Some("pathogenic")])]),
            vec!["benign;likely_benign,pathogenic"],
        ),
        (
            "PUBMED",
            text_lists(&[
                Some(&[Some("12345"), Some("67890")]),
                Some(&[None, Some("1")]),
            ]),
            vec!["12345,67890", "-,1"],
        ),
        (
            "EXON",
            positions(&[
                span(Some(3), Some(3), Some(12)),
                span(Some(3), Some(4), Some(12)),
            ]),
            vec!["3/12", "3-4/12"],
        ),
        (
            "INTRON",
            positions(&[span(Some(7), Some(7), Some(11))]),
            vec!["7/11"],
        ),
        (
            "cDNA_position",
            positions(&[
                span(Some(44), Some(44), None),
                span(Some(44), Some(45), None),
                span(Some(44), Some(44), Some(1200)),
                span(None, None, Some(1200)),
                span(Some(44), Some(45), Some(1200)),
            ]),
            vec!["44", "44-45", "44/1200", "-/1200", "44-45/1200"],
        ),
        (
            "CDS_position",
            positions(&[span(Some(146), None, None)]),
            vec!["146-?"],
        ),
        (
            "Protein_position",
            positions(&[
                span(None, Some(52), None),
                span(Some(49), None, Some(300)),
                span(None, None, None),
            ]),
            vec!["?-52", "49-?/300", NULL],
        ),
        // A plugin column is a string; a BOOLEAN column of no listed spelling prints a flag.
        ("CADD_PHRED", texts(&[Some("23.4")]), vec!["23.4"]),
        (
            "FLAGGED",
            Stored::Bool(vec![Some(true), Some(false)]),
            vec!["1", "0"],
        ),
    ]
}

const ROWS: usize = 5;

fn planted_cells(path: &Path, group: usize, how: Encoding) {
    let columns: Vec<(&str, Stored)> = cells()
        .into_iter()
        .map(|(name, stored, _)| (name, stored))
        .collect();
    plant(path, &columns, group, how);
}

#[test]
fn every_typed_cell_prints_its_tab_text() {
    let dir = scratch("cells");
    let encodings = [
        AS_THE_WRITER,
        Encoding {
            codec: Codec::None,
            data_page_v2: false,
            dictionary: false,
        },
        Encoding {
            codec: Codec::Zstd,
            data_page_v2: false,
            dictionary: true,
        },
    ];
    for (at, how) in encodings.into_iter().enumerate() {
        for group in [ROWS, 2, 1] {
            let path = dir.join(format!("cells-{at}-{group}.parquet"));
            planted_cells(&path, group, how);
            let text = tab(&path, false);
            let lines = lines(&text);
            assert_eq!(lines.len(), 1 + ROWS);
            for (row, line) in lines[1..].iter().enumerate() {
                let printed: Vec<&str> = line.split('\t').collect();
                assert_eq!(printed.len(), cells().len());
                for ((column, _, texts), cell) in cells().iter().zip(printed) {
                    let expected = texts.get(row).copied().unwrap_or(NULL);
                    assert_eq!(cell, expected, "{column} row {row}");
                }
            }
        }
    }
}

/// The codec the pages of the first column of a file's first row group are written in.
fn codec_of(path: &Path) -> parquet::basic::Compression {
    let file = std::fs::File::open(path).expect("the planted file opens");
    let reader = SerializedFileReader::new(file).expect("the planted file is a Parquet file");
    reader.metadata().row_group(0).column(0).compression()
}

#[test]
fn a_file_is_read_under_every_page_codec_the_adapter_reads() {
    let dir = scratch("codecs");
    let plain = dir.join("plain.parquet");
    planted_cells(
        &plain,
        ROWS,
        Encoding {
            codec: Codec::None,
            data_page_v2: true,
            dictionary: true,
        },
    );
    let wanted = [tab(&plain, false), tab(&plain, true)];
    assert_eq!(lines(&wanted[0]).len(), 1 + ROWS);
    for codec in Codec::ALL {
        for data_page_v2 in [false, true] {
            let path = dir.join(format!("{codec:?}-{data_page_v2}.parquet"));
            planted_cells(
                &path,
                2,
                Encoding {
                    codec,
                    data_page_v2,
                    dictionary: !data_page_v2,
                },
            );
            assert_eq!(codec_of(&path), codec.compression(), "{codec:?}");
            assert_eq!(tab(&path, false), wanted[0], "{codec:?}");
            assert_eq!(tab(&path, true), wanted[1], "{codec:?}");
        }
    }
    // The parts of one output may each have a codec of their own.
    let mixed = dir.join("mixed.parquet");
    for (part, codec) in Codec::ALL.into_iter().enumerate() {
        let how = Encoding {
            codec,
            data_page_v2: true,
            dictionary: true,
        };
        let name = format!("part-{part:05}.parquet");
        plant(
            &mixed.join(name),
            &[("Feature", texts(&[Some(&format!("{codec:?}"))]))],
            1,
            how,
        );
    }
    assert_eq!(
        lines(&tab(&mixed, false))[1..],
        ["None", "Zstd", "Snappy", "Gzip", "Brotli", "Lz4", "Lz4Raw"]
    );
}

#[test]
fn the_header_is_the_tab_columns_after_the_key_columns() {
    let dir = scratch("header");
    let path = dir.join("cells.parquet");
    planted_cells(&path, ROWS, AS_THE_WRITER);
    let text = tab(&path, false);
    let header = lines(&text)[0];
    let names: Vec<&str> = cells().iter().map(|(name, _, _)| *name).collect();
    assert_eq!(header, format!("#{}", names.join("\t")));
    for key in ["contig_index", "chrom", "pos", "end", "ref", "alt"] {
        assert!(!header.split('\t').any(|name| name == key), "{key}");
    }
}

#[test]
fn cell_text_prints_each_value_kind() {
    let dir = scratch("kinds");
    let path = dir.join("kinds.parquet");
    let columns = [
        ("Feature", texts(&[None, Some("")])),
        ("Allele", texts(&[Some("-")])),
        ("STRAND", Stored::Int32(vec![Some(-1)])),
        ("MAX_AF", doubles(&[Some(0.095)])),
        ("CANONICAL", Stored::Bool(vec![Some(true), Some(false)])),
        ("HIGH_INF_POS", Stored::Bool(vec![Some(false)])),
        (
            "SOMATIC",
            bool_lists(&[Some(&[Some(false), Some(true), None]), Some(&[])]),
        ),
        (
            "AF",
            double_lists(&[Some(&[Some(0.0002), None]), Some(&[])]),
        ),
        (
            "Consequence",
            text_lists(&[Some(&[Some("a"), None]), Some(&[Some("")])]),
        ),
        ("EXON", positions(&[Some((Some(3), Some(4), Some(12)))])),
        ("cDNA_position", positions(&[Some((None, None, None))])),
    ];
    plant(&path, &columns, 9, AS_THE_WRITER);
    assert_eq!(
        lines(&tab(&path, false))[1..],
        [
            "-\t-\t-1\t0.095\tYES\tN\t0,1,-\t0.0002,-\ta,-\t3-4/12\t-",
            // An empty string is a cell, an empty list an empty cell, and a NULL struct
            // prints as the struct of three NULL members does.
            "\t-\t-\t-\t-\t-\t\t\t\t-\t-",
        ]
    );
}

#[test]
fn shortest_decimal_prints_the_writers_text() {
    let dir = scratch("decimals");
    let path = dir.join("decimals.parquet");
    let cases = [
        (0.095, "0.095"),
        (1.0, "1"),
        (9.8e-05, "0.000098"),
        (0.0, "0"),
        (-0.0, "-0"),
        (1e16, "10000000000000000"),
        (1e21, "1000000000000000000000"),
        (0.1 + 0.2, "0.30000000000000004"),
        (f64::NAN, "NaN"),
        (f64::INFINITY, "inf"),
        (f64::NEG_INFINITY, "-inf"),
        // Two shortest decimals equally close: the one of larger magnitude.
        (-126430982921026.62, "-126430982921026.63"),
        (126430982921026.62, "126430982921026.63"),
    ];
    let values: Vec<Option<f64>> = cases.iter().map(|(value, _)| Some(*value)).collect();
    let listed: Vec<Option<f64>> = values.clone();
    plant(
        &path,
        &[
            ("MAX_AF", doubles(&values)),
            ("Extra", doubles(&values)),
            ("AF", double_lists(&[Some(&listed)])),
        ],
        4,
        AS_THE_WRITER,
    );
    let expected: Vec<&str> = cases.iter().map(|(_, text)| *text).collect();
    for canonical in [false, true] {
        let text = tab(&path, canonical);
        let mut printed: Vec<(&str, &str)> = Vec::new();
        for line in &lines(&text)[1..] {
            let cells: Vec<&str> = line.split('\t').collect();
            printed.push((cells[0], cells[1]));
            if cells[2] != NULL {
                assert_eq!(cells[2], expected.join(","));
            }
        }
        let mut wanted: Vec<(&str, &str)> = expected.iter().map(|text| (*text, *text)).collect();
        if canonical {
            wanted.sort();
            printed.sort();
        }
        assert_eq!(printed, wanted, "canonical {canonical}");
    }
}

#[test]
fn the_parts_of_a_directory_come_out_in_path_order() {
    let dir = scratch("parts");
    let out = dir.join("out.parquet");
    plant(
        &out.join("part-00001.parquet"),
        &[("Feature", texts(&[Some("second")]))],
        1,
        AS_THE_WRITER,
    );
    plant(
        &out.join("part-00000.parquet"),
        &[("Feature", texts(&[Some("first")]))],
        1,
        AS_THE_WRITER,
    );
    std::fs::write(out.join("_vep_rs_manifest.json"), "{}").unwrap();
    assert_eq!(lines(&tab(&out, false)), ["#Feature", "first", "second"]);
}

#[test]
fn a_hive_directory_is_read_file_by_file() {
    let dir = scratch("hive");
    let out = dir.join("out.parquet");
    plant(
        &out.join("chrom=21").join("part-00000.parquet"),
        &[("Location", texts(&[Some("21:100")]))],
        1,
        AS_THE_WRITER,
    );
    plant(
        &out.join("chrom=X").join("part-00000.parquet"),
        &[("Location", texts(&[Some("X:100"), Some("X:101")]))],
        2,
        AS_THE_WRITER,
    );
    assert_eq!(
        lines(&tab(&out, false)),
        ["#Location", "21:100", "X:100", "X:101"]
    );
}

#[test]
fn the_files_of_a_directory_are_its_parquet_files_in_the_order_of_their_parts() {
    let dir = scratch("order");
    let out = dir.join("out");
    let names = [
        "a-b.parquet",
        "a/b.parquet",
        "a.parquet",
        "B.parquet",
        "a/.hidden.parquet",
        "a b.parquet",
        "\u{e9}.parquet",
        "a/c/d.parquet",
        "z.parquet/inner.parquet",
        "x.PARQUET",
        ".parquet",
        "q.parquet.tmp",
    ];
    for name in names {
        plant(
            &out.join(name),
            &[("Feature", texts(&[Some(name)]))],
            1,
            AS_THE_WRITER,
        );
    }
    std::os::unix::fs::symlink(out.join("a.parquet"), out.join("link.parquet")).unwrap();
    std::os::unix::fs::symlink(out.join("a"), out.join("linkdir")).unwrap();
    std::os::unix::fs::symlink(out.join("nowhere"), out.join("dangling.parquet")).unwrap();
    // A path sorts part by part, a file reached through a link is read and a directory
    // reached through one is not searched.
    assert_eq!(
        lines(&tab(&out, false))[1..],
        [
            ".parquet",
            "B.parquet",
            "a/.hidden.parquet",
            "a/b.parquet",
            "a/c/d.parquet",
            "a b.parquet",
            "a-b.parquet",
            "a.parquet",
            "a.parquet",
            "z.parquet/inner.parquet",
            "\u{e9}.parquet",
        ]
    );
}

fn gzipped(text: &str) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(text.as_bytes()).unwrap();
    encoder.finish().unwrap()
}

#[test]
fn canonical_sorts_the_rows_and_reprints_every_double_cell() {
    let dir = scratch("canonical");
    let parquet = dir.join("out.parquet");
    let columns = [
        ("Feature", texts(&[Some("b"), Some("a"), Some("c")])),
        ("MAX_AF", doubles(&[Some(0.095), Some(9.8e-05), None])),
        (
            "AF",
            double_lists(&[Some(&[Some(0.2766), Some(0.095)]), Some(&[Some(0.1)]), None]),
        ),
    ];
    plant(&parquet, &columns, 2, AS_THE_WRITER);
    let header = "#Feature\tMAX_AF\tAF";
    let canonical = tab(&parquet, true);
    assert_eq!(
        lines(&canonical),
        [
            header,
            "a\t0.000098\t0.1",
            "b\t0.095\t0.2766,0.095",
            "c\t-\t-"
        ]
    );
    let plain = tab(&parquet, false);
    assert_eq!(sorted(data_lines(&plain)), data_lines(&canonical));
    assert_eq!(lines(&plain)[1], "b\t0.095\t0.2766,0.095");

    let tab_text = format!(
        "## ENSEMBL VARIANT EFFECT PREDICTOR\n## Output produced at ...\n\
         {header}\nc\t-\t-\nb\t0.0950\t0.2766,0.0950\na\t9.8e-05\t0.100\n"
    );
    let tab_file = dir.join("out.tab");
    std::fs::write(&tab_file, &tab_text).unwrap();
    assert_eq!(tab(&tab_file, true), canonical);
    let gz = dir.join("out.tab.gz");
    std::fs::write(&gz, gzipped(&tab_text)).unwrap();
    assert_eq!(tab(&gz, true), canonical);

    let again = dir.join("again.tab");
    std::fs::write(&again, &canonical).unwrap();
    assert_eq!(tab(&again, true), canonical);

    // A line ends at a newline alone: a cell holding another separator character stays
    // one row, and a CRLF file reads like an LF one.
    let separators = dir.join("separators.tab");
    std::fs::write(
        &separators,
        format!("{header}\r\nb\u{1c}c\t0.0950\t0.1\r\n"),
    )
    .unwrap();
    assert_eq!(
        tab(&separators, true),
        format!("{header}\nb\u{1c}c\t0.095\t0.1\n")
    );
}

#[test]
fn errors_are_reported_in_the_adapters_words() {
    let dir = scratch("errors");
    let absent = dir.join("absent.parquet");
    for canonical in [false, true] {
        assert_eq!(
            refusal(&absent, canonical),
            format!("{ADAPTER}{}: no such file or directory", absent.display())
        );
    }
    let empty = dir.join("empty");
    std::fs::create_dir(&empty).unwrap();
    assert_eq!(
        refusal(&empty, false),
        format!(
            "{ADAPTER}{}: no .parquet file under the directory",
            empty.display()
        )
    );
    let mixed = dir.join("mixed.parquet");
    plant(
        &mixed.join("part-00000.parquet"),
        &[("Feature", texts(&[Some("a")]))],
        1,
        AS_THE_WRITER,
    );
    plant(
        &mixed.join("part-00001.parquet"),
        &[("Gene", texts(&[Some("g")]))],
        1,
        AS_THE_WRITER,
    );
    assert_eq!(
        refusal(&mixed, false),
        format!(
            "{ADAPTER}{0}/part-00001.parquet: its tab columns differ from {0}/part-00000.parquet's",
            mixed.display()
        )
    );
    let headless = dir.join("headless.tab");
    std::fs::write(&headless, "a\t1\n").unwrap();
    assert_eq!(
        refusal(&headless, true),
        format!(
            "{ADAPTER}{}: no # header line; a tab output needs its header",
            headless.display()
        )
    );
    let twice = dir.join("twice.tab");
    std::fs::write(&twice, "#a\nx\n#b\n").unwrap();
    assert_eq!(
        refusal(&twice, true),
        format!("{ADAPTER}{}: two header lines", twice.display())
    );
    let keys_only = dir.join("keys.parquet");
    plant_rows(&keys_only, &[], 1, 1, AS_THE_WRITER);
    assert_eq!(
        refusal(&keys_only, false),
        format!("{ADAPTER}no tab columns in the Parquet schema")
    );
    // A regular file without the Parquet magic is refused before it is opened as one; a
    // file with the magic but no footer, a tab file that is not UTF-8 and a gzip file cut
    // short are refused in one line too.
    let plain = dir.join("plain.tab");
    std::fs::write(&plain, "#Feature\ta\n").unwrap();
    assert_eq!(
        refusal(&plain, false),
        format!("{ADAPTER}{}: not a Parquet file", plain.display())
    );
    let magic = "Parquet magic bytes not found in footer. Either the file is corrupted or this is not a parquet file.";
    let corrupt = dir.join("corrupt.parquet");
    std::fs::write(&corrupt, b"PAR1\0\0\0\0\0\0\0\0").unwrap();
    for canonical in [false, true] {
        assert_eq!(refusal(&corrupt, canonical), format!("{ADAPTER}{magic}"));
    }
    let binary = dir.join("binary.tab");
    std::fs::write(&binary, b"\xff\xfe\x00").unwrap();
    assert_eq!(
        refusal(&binary, true),
        format!("{ADAPTER}'utf-8' codec can't decode byte 0xff in position 0: invalid start byte")
    );
    let truncated = dir.join("truncated.tab.gz");
    let whole = gzipped(&"#Feature\ta\n".repeat(100));
    std::fs::write(&truncated, &whole[..whole.len() - 8]).unwrap();
    assert_eq!(
        refusal(&truncated, true),
        format!("{ADAPTER}Compressed file ended before the end-of-stream marker was reached")
    );
    let device = Path::new("/dev/null");
    for canonical in [false, true] {
        assert_eq!(
            refusal(device, canonical),
            format!("{ADAPTER}/dev/null: neither a regular file nor a directory")
        );
    }
}

#[test]
fn a_file_that_is_no_parquet_file_is_refused_by_its_size_and_its_last_bytes() {
    let dir = scratch("footers");
    let planted = dir.join("planted.parquet");
    plant(
        &planted,
        &[("Feature", texts(&[Some("a")]))],
        1,
        AS_THE_WRITER,
    );
    let whole = std::fs::read(&planted).unwrap();
    let size = whole.len();
    let magic = "Parquet magic bytes not found in footer. Either the file is corrupted or this is not a parquet file.";
    let cases: Vec<(&str, Vec<u8>, String)> = vec![
        ("empty", Vec::new(), "Parquet file size is 0 bytes".to_string()),
        (
            "magic-only",
            b"PAR1".to_vec(),
            "Parquet file size is 4 bytes, smaller than the minimum file footer (8 bytes)".to_string(),
        ),
        ("text", b"hello world, not parquet at all".to_vec(), magic.to_string()),
        ("cut", whole[..size - 1].to_vec(), magic.to_string()),
        (
            "magic-magic",
            b"PAR1PAR1".to_vec(),
            "Parquet file size is 8 bytes, smaller than the size reported by footer's (827474256bytes)".to_string(),
        ),
        (
            "long-footer",
            [&whole[..size - 8], &1_000_000u32.to_le_bytes(), b"PAR1"].concat(),
            format!("Parquet file size is {size} bytes, smaller than the size reported by footer's (1000000bytes)"),
        ),
        (
            "encrypted",
            b"PAR1abcdefgh\x08\0\0\0PARE".to_vec(),
            "Could not read encrypted metadata, no decryption found in reader's properties".to_string(),
        ),
    ];
    for (name, bytes, wanted) in cases {
        // A file named alone opens with the magic, or it is no Parquet file at all.
        let holder = dir.join(name);
        std::fs::create_dir(&holder).unwrap();
        std::fs::write(holder.join("part.parquet"), &bytes).unwrap();
        assert_eq!(
            refusal(&holder, false),
            format!("{ADAPTER}{wanted}"),
            "{name}"
        );
        if bytes.starts_with(b"PAR1") {
            let alone = dir.join(format!("{name}.parquet"));
            std::fs::write(&alone, &bytes).unwrap();
            assert_eq!(
                refusal(&alone, true),
                format!("{ADAPTER}{wanted}"),
                "{name}"
            );
        }
    }
    // The files of a directory are all opened before the first is read.
    let second = dir.join("second");
    plant(
        &second.join("part-0.parquet"),
        &[("Feature", texts(&[Some("a")]))],
        1,
        AS_THE_WRITER,
    );
    std::fs::write(second.join("part-1.parquet"), b"junk").unwrap();
    assert_eq!(
        refusal(&second, false),
        format!(
            "{ADAPTER}Parquet file size is 4 bytes, smaller than the minimum file footer (8 bytes)"
        )
    );
}

#[test]
fn a_failure_is_one_line_and_the_text_reads_every_line_end_alike() {
    let dir = scratch("line-ends");
    // The adapter's words are joined by single spaces, the path the comparator names is not.
    let spaced = dir.join("two  spaces\u{a0}dir");
    std::fs::create_dir(&spaced).unwrap();
    let Err(Error::Comparator(text)) = parquet_as_tab(&spaced, false) else {
        panic!("an empty directory is read");
    };
    assert_eq!(
        text,
        format!(
            "{FAILED}{0}/two  spaces\u{a0}dir: {ADAPTER}{0}/two spaces dir: no .parquet file under the directory",
            dir.display()
        )
    );
    // A carriage return inside a cell is a line end of the text.
    let path = dir.join("returns.parquet");
    let cells = [Some("a\rb"), Some("c\r\nd"), Some("e\nf"), Some("g\th")];
    plant(
        &path,
        &[
            ("Feature", texts(&cells)),
            ("MAX_AF", texts(&[Some("1.50\n"), Some("x")])),
        ],
        9,
        AS_THE_WRITER,
    );
    assert_eq!(
        tab(&path, false),
        "#Feature\tMAX_AF\na\nb\t1.50\n\nc\nd\tx\ne\nf\t-\ng\th\t-\n"
    );
    // The canonical form keeps the cell of a string column as stored, whatever the
    // column's name, a line end inside it included.
    assert_eq!(
        tab(&path, true),
        "#Feature\tMAX_AF\na\nb\t1.50\n\nc\nd\tx\ne\nf\t-\ng\th\t-\n"
    );
    // It reads a cell that holds a tab as two cells, each in the place it comes to stand
    // in. In the place of a string column the piece stays as stored.
    let tabbed = dir.join("tabbed.parquet");
    plant(
        &tabbed,
        &[
            ("Feature", texts(&[Some("g\t0.50")])),
            ("MAX_AF", texts(&[Some("0.250")])),
        ],
        9,
        AS_THE_WRITER,
    );
    assert_eq!(tab(&tabbed, true), "#Feature\tMAX_AF\ng\t0.50\t0.250\n");
    // In the place of a DOUBLE column the piece is a number, also before a line end, and
    // the column's own cell stands past the last column.
    let numbered = dir.join("numbered.parquet");
    plant(
        &numbered,
        &[
            ("Feature", texts(&[Some("g\t0.50"), Some("h\t1.50\n")])),
            ("MAX_AF", Stored::Double(vec![Some(0.25), Some(0.5)])),
        ],
        9,
        AS_THE_WRITER,
    );
    assert_eq!(
        tab(&numbered, true),
        "#Feature\tMAX_AF\ng\t0.5\t0.25\nh\t1.5\t0.5\n"
    );
}

#[test]
fn a_path_the_adapters_parser_reads_as_an_option_is_refused_with_its_usage() {
    let Err(Error::Comparator(text)) = parquet_as_tab(Path::new("-x.parquet"), true) else {
        panic!("the path is read");
    };
    assert_eq!(
        text,
        "the Parquet adapter failed on -x.parquet: usage: parquet_to_vep_tab.py [-h] \
         [--input INPUT_PATH] [--output OUTPUT] [--canonical] [path] parquet_to_vep_tab.py: \
         error: unrecognized arguments: -x.parquet"
    );
    assert_eq!(
        refusal(Path::new("-1.parquet"), false),
        format!("{ADAPTER}-1.parquet: no such file or directory")
    );
}

/// The committed `--tab` output of Ensembl VEP on the release corpus.
fn corpus_tab() -> String {
    let corpus =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-release116");
    read_text(&corpus.join("expected/tab.txt.gz")).unwrap()
}

fn column_line(text: &str) -> String {
    text_lines(text)
        .into_iter()
        .rfind(|line| line.starts_with('#'))
        .expect("a column line")
}

#[test]
fn a_parquet_actual_is_read_through_the_adapter() {
    let dir = scratch("corpus");
    let text = corpus_tab();
    let rows = data_lines(&text);
    let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
    let header = column_line(&text);
    let out = dir.join("out.parquet");
    plant_tab_rows(
        &out.join("part-00000.parquet"),
        &header,
        &rows[..rows.len() / 2],
    );
    plant_tab_rows(
        &out.join("part-00001.parquet"),
        &header,
        &rows[rows.len() / 2..],
    );
    let as_tab = tab(&out, true);
    let headers: Vec<String> = text_lines(&as_tab)
        .into_iter()
        .filter(|line| line.starts_with('#'))
        .collect();
    assert_eq!(headers, std::slice::from_ref(&header));
    let wanted: Vec<String> = rows.iter().map(|row| row.to_string()).collect();
    assert_eq!(rows.len(), 977);
    assert_eq!(sorted(data_lines(&as_tab)), sorted(wanted.clone()));
    assert_eq!(data_lines(&tab(&out, false)), wanted);
    // The reference's own text in the canonical form holds the same rows.
    let expected_file = dir.join("tab.txt");
    std::fs::write(&expected_file, &text).unwrap();
    assert_eq!(tab(&expected_file, true), as_tab);
    let Err(Error::Comparator(refused)) = parquet_as_tab(&dir.join("absent.parquet"), true) else {
        panic!("an absent output is read");
    };
    assert!(refused.starts_with("the Parquet adapter failed"));
}

#[test]
fn a_parquet_comparison_reads_frequencies_as_numbers() {
    let dir = scratch("frequencies");
    let header = "#Uploaded_variation\tLocation\tAllele\tConsequence\tAF\tgnomADe_AF\tMAX_AF";
    let row = "rs1\t21:100\tA\tintergenic_variant\t1.47e-05\t0.0010,0.2770\t0.2770";
    let expected_file = dir.join("tab.txt");
    std::fs::write(
        &expected_file,
        format!("## ENSEMBL VARIANT EFFECT PREDICTOR\n{header}\n{row}\n"),
    )
    .unwrap();
    let out = dir.join("out.parquet");
    plant_tab_rows(&out, header, &[row]);
    assert_eq!(
        data_lines(&tab(&out, true)),
        ["rs1\t21:100\tA\tintergenic_variant\t0.0000147\t0.001,0.277\t0.277"]
    );
    assert_eq!(
        data_lines(&tab(&expected_file, true)),
        data_lines(&tab(&out, true))
    );
    plant_tab_rows(&out, header, &[&row.replace("\t0.2770", "\t0.2771")]);
    assert_ne!(
        data_lines(&tab(&expected_file, true)),
        data_lines(&tab(&out, true))
    );
}

#[test]
fn self_parity_reads_the_parquet_rows_in_file_order() {
    let dir = scratch("file-order");
    let header =
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tIND";
    let rows = [
        "rs1\t21:100\tA\tENSG1\tENST1\tTranscript\tintergenic_variant\tsample_b",
        "rs1\t21:100\tA\tENSG1\tENST1\tTranscript\tintergenic_variant\tsample_a",
    ];
    let out = dir.join("out.parquet");
    plant_tab_rows(&out, header, &rows);
    assert_eq!(data_lines(&tab(&out, false)), rows);
    assert_eq!(data_lines(&tab(&out, true)), [rows[1], rows[0]]);
}

#[test]
fn rows_across_row_groups_and_files_keep_their_order_or_are_sorted_as_one() {
    let dir = scratch("many");
    let out = dir.join("out.parquet");
    let count = 150_000;
    let label = |row: usize| format!("row{:06}", (row * 7919) % count);
    let mut wanted = Vec::new();
    for (part, group) in [(0usize, 40_000usize), (1, 7), (2, 100_000)] {
        let rows = [70_000, 1_000, 79_000][part];
        let from: usize = [0, 70_000, 71_000][part];
        let names: Vec<Option<String>> = (from..from + rows).map(|row| Some(label(row))).collect();
        let numbers: Vec<Option<f64>> = (from..from + rows)
            .map(|row| (row % 13 != 0).then_some((row % 1000) as f64 / 1000.0))
            .collect();
        wanted.extend((from..from + rows).map(|row| {
            let number = match row % 13 {
                0 => "-".to_string(),
                _ => format!("{}", (row % 1000) as f64 / 1000.0),
            };
            format!("{}\t{number}", label(row))
        }));
        plant(
            &out.join(format!("part-{part:05}.parquet")),
            &[
                ("Feature", Stored::Text(names)),
                ("MAX_AF", Stored::Double(numbers)),
            ],
            group,
            AS_THE_WRITER,
        );
    }
    assert_eq!(data_lines(&tab(&out, false)), wanted);
    assert_eq!(data_lines(&tab(&out, true)), sorted(wanted));
}

/// A string column of the bytes `values`, UTF-8 or not, in row groups of `group` rows; a
/// column named in `numbers` is a DOUBLE column of the numbers its values spell.
fn planted_bytes(path: &Path, columns: &[(&str, Vec<Vec<u8>>)], numbers: &[&str], group: usize) {
    let rows = columns[0].1.len();
    let fields = columns
        .iter()
        .map(|(name, _)| {
            if numbers.contains(name) {
                value_field(name, PhysicalType::DOUBLE, Repetition::OPTIONAL, None)
            } else {
                text_field(name)
            }
        })
        .collect();
    let number = |text: &Vec<u8>| -> f64 {
        std::str::from_utf8(text)
            .expect("the text of a number")
            .parse()
            .expect("a number")
    };
    let mut row_groups = Vec::new();
    for from in (0..rows).step_by(group) {
        let to = rows.min(from + group);
        row_groups.push(
            columns
                .iter()
                .map(|(name, values)| LeafData {
                    values: if numbers.contains(name) {
                        Values::Double(values[from..to].iter().map(number).collect())
                    } else {
                        Values::Bytes(values[from..to].to_vec())
                    },
                    def: Some(vec![1; to - from]),
                    rep: None,
                })
                .collect(),
        );
    }
    write_file(path, fields, row_groups, AS_THE_WRITER);
}

#[test]
fn the_failure_met_first_is_the_one_the_adapters_batches_meet_first() {
    let dir = scratch("first-failure");
    let not_utf8 =
        format!("{ADAPTER}'utf-8' codec can't decode byte 0xff in position 0: invalid start byte");
    let cut_short = format!(
        "{ADAPTER}'utf-8' codec can't decode byte 0xc3 in position 1: unexpected end of data"
    );
    let no_number = format!("{ADAPTER}could not convert string to float: '\u{131}nf'");
    let dotted = format!("{ADAPTER}could not convert string to float: '\u{130}nf'");
    // `MAX_AF` is a DOUBLE column, and the canonical form reads a number where a piece of
    // a cell comes to stand in its place: a text named for `MAX_AF` is planted after a tab
    // in the `Feature` cell of its row.
    let plant = |name: &str, rows: usize, group: usize, bad: &[(&str, usize, &[u8])]| {
        let columns: Vec<(&str, Vec<Vec<u8>>)> = ["Feature", "MAX_AF", "Gene"]
            .iter()
            .map(|column| {
                let plain: &[u8] = if *column == "MAX_AF" { b"0.5" } else { b"x" };
                let mut values = vec![plain.to_vec(); rows];
                for (named, row, value) in bad {
                    if *named == "MAX_AF" {
                        if *column == "Feature" {
                            values[*row] = [b"x\t", *value].concat();
                        }
                    } else if named == column {
                        values[*row] = value.to_vec();
                    }
                }
                (*column, values)
            })
            .collect();
        let path = dir.join(name);
        planted_bytes(&path, &columns, &["MAX_AF"], group);
        path
    };
    let dotless = "\u{131}nf".as_bytes();
    // A number of the first batch of 65536 rows is refused before a string of the second
    // is decoded; without the canonical form no number is read.
    let path = plant(
        "batches.parquet",
        70_000,
        70_000,
        &[("Gene", 66_000, b"\xff\xfe"), ("MAX_AF", 100, dotless)],
    );
    assert_eq!(refusal(&path, true), no_number);
    assert_eq!(refusal(&path, false), not_utf8);
    // Every string of a batch is decoded before a number of its rows is read, the batch
    // going on into the next row group.
    let path = plant(
        "one-batch.parquet",
        70_000,
        70_000,
        &[("Gene", 60_000, b"\xff\xfe"), ("MAX_AF", 100, dotless)],
    );
    assert_eq!(refusal(&path, true), not_utf8);
    let path = plant(
        "next-group.parquet",
        80_000,
        40_000,
        &[("Gene", 50_000, b"a\xc3"), ("MAX_AF", 10, dotless)],
    );
    assert_eq!(refusal(&path, true), cut_short);
    // The columns of a batch are decoded in their order, each over all its rows.
    let path = plant(
        "two-columns.parquet",
        1_000,
        300,
        &[("Gene", 5, b"\xff\xfe"), ("Feature", 900, b"a\xc3")],
    );
    assert_eq!(refusal(&path, false), cut_short);
    let path = plant(
        "two-batches.parquet",
        140_000,
        30_000,
        &[("Gene", 5, b"\xff\xfe"), ("Feature", 70_000, b"a\xc3")],
    );
    assert_eq!(refusal(&path, false), not_utf8);
    // The numbers of a batch are read row by row.
    let dotted_bytes = "\u{130}nf".as_bytes();
    let path = plant(
        "two-numbers.parquet",
        140_000,
        50_000,
        &[
            ("MAX_AF", 131_073, dotless),
            ("MAX_AF", 131_072, dotted_bytes),
        ],
    );
    assert_eq!(refusal(&path, true), dotted);
    assert_eq!(lines(&tab(&path, false)).len(), 140_001);
    // The files are read one after another.
    let both = dir.join("two-files");
    planted_bytes(
        &both.join("part-1.parquet"),
        &[
            ("Feature", vec![b"a".to_vec(), b"\xff\xfe".to_vec()]),
            ("MAX_AF", vec![b"1".to_vec(), b"2".to_vec()]),
        ],
        &["MAX_AF"],
        9,
    );
    planted_bytes(
        &both.join("part-0.parquet"),
        &[
            ("Feature", vec![b"a".to_vec(), [b"b\t", dotless].concat()]),
            ("MAX_AF", vec![b"1".to_vec(), b"2".to_vec()]),
        ],
        &["MAX_AF"],
        9,
    );
    assert_eq!(refusal(&both, true), no_number);
    assert_eq!(refusal(&both, false), not_utf8);
}

#[test]
fn a_list_and_a_struct_of_another_layout_print_as_the_writers_do() {
    let dir = scratch("layouts");
    let path = dir.join("layouts.parquet");
    let int32 = |name: &str, repetition| value_field(name, PhysicalType::INT32, repetition, None);
    let fields = vec![
        // A repeated value without a list group: a list that is never NULL.
        int32("Repeated", Repetition::REPEATED),
        // A list group whose repeated field is the element.
        group_field(
            "Older",
            Some(Repetition::OPTIONAL),
            Some(LogicalType::List),
            vec![int32("array", Repetition::REPEATED)],
        ),
        // A required list of required elements.
        group_field(
            "Required",
            Some(Repetition::REQUIRED),
            Some(LogicalType::List),
            vec![group_field(
                "list",
                Some(Repetition::REPEATED),
                None,
                vec![int32("element", Repetition::REQUIRED)],
            )],
        ),
        // A required struct of a required and an optional member, and of another member.
        group_field(
            "EXON",
            Some(Repetition::REQUIRED),
            None,
            vec![
                int32("other", Repetition::OPTIONAL),
                int32("total", Repetition::OPTIONAL),
                int32("start", Repetition::REQUIRED),
            ],
        ),
        // A struct with none of the three members.
        group_field(
            "INTRON",
            Some(Repetition::OPTIONAL),
            None,
            vec![int32("other", Repetition::OPTIONAL)],
        ),
        value_field("Flag", PhysicalType::BOOLEAN, Repetition::REQUIRED, None),
        value_field(
            "Unsigned",
            PhysicalType::INT32,
            Repetition::OPTIONAL,
            Some(LogicalType::integer(32, false)),
        ),
        value_field(
            "Wide",
            PhysicalType::INT64,
            Repetition::OPTIONAL,
            Some(LogicalType::integer(64, false)),
        ),
    ];
    let leaf = |values: Values, def: Option<Vec<i16>>, rep: Option<Vec<i16>>| LeafData {
        values,
        def,
        rep,
    };
    // Three rows.
    let leaves = vec![
        // [1, 2], [], [3]
        leaf(
            Values::Int32(vec![1, 2, 3]),
            Some(vec![1, 1, 0, 1]),
            Some(vec![0, 1, 0, 0]),
        ),
        // NULL, [], [4, 5]
        leaf(
            Values::Int32(vec![4, 5]),
            Some(vec![0, 1, 2, 2]),
            Some(vec![0, 0, 0, 1]),
        ),
        // [6], [], [7, 8]
        leaf(
            Values::Int32(vec![6, 7, 8]),
            Some(vec![1, 0, 1, 1]),
            Some(vec![0, 0, 0, 1]),
        ),
        leaf(Values::Int32(vec![9]), Some(vec![1, 0, 0]), None),
        leaf(Values::Int32(vec![12]), Some(vec![0, 1, 0]), None),
        leaf(Values::Int32(vec![3, 4, 5]), None, None),
        leaf(Values::Int32(vec![1]), Some(vec![2, 1, 0]), None),
        leaf(Values::Bool(vec![true, false, true]), None, None),
        leaf(Values::Int32(vec![-1, 7]), Some(vec![1, 0, 1]), None),
        leaf(Values::Int64(vec![-1, 7]), Some(vec![1, 1, 0]), None),
    ];
    write_file(
        &path,
        fields,
        vec![leaves],
        Encoding {
            codec: Codec::None,
            data_page_v2: false,
            dictionary: false,
        },
    );
    assert_eq!(
        lines(&tab(&path, false)),
        [
            "#Repeated\tOlder\tRequired\tEXON\tINTRON\tFlag\tUnsigned\tWide",
            "1,2\t-\t6\t3-?\t-\t1\t4294967295\t18446744073709551615",
            "\t\t\t4-?/12\t-\t0\t-\t7",
            "3\t4,5\t7,8\t5-?\t-\t1\t7\t-",
        ]
    );
}

#[test]
fn the_key_columns_are_the_ones_that_lead_the_schema() {
    let dir = scratch("keys");
    let path = dir.join("keys.parquet");
    // Two key columns in another order than the writer's, then the tab columns, one of
    // which has a key column's name.
    let names = ["alt", "chrom", "Feature", "pos", "Gene"];
    let fields = names.iter().map(|name| text_field(name)).collect();
    let leaves = names
        .iter()
        .map(|name| LeafData {
            values: Values::Bytes(vec![name.to_lowercase().into_bytes()]),
            def: Some(vec![1]),
            rep: None,
        })
        .collect();
    write_file(&path, fields, vec![leaves], AS_THE_WRITER);
    assert_eq!(
        lines(&tab(&path, false)),
        ["#Feature\tpos\tGene", "feature\tpos\tgene"]
    );
}

/// The writer stores a tab column named as a key column as `tab_<name>`, so that no two
/// columns of a file share a name. The column line names each tab column as the `--tab`
/// format does: one prefix goes where a key column's name stands behind it, and every
/// other name prints as stored.
#[test]
fn a_tab_column_stored_behind_the_prefix_prints_under_its_own_name() {
    let dir = scratch("stored-names");
    let path = dir.join("names.parquet");
    let names = [
        "chrom",
        "ref",
        "Feature",
        "tab_ref",
        "tab_tab_alt",
        "tab_Gene",
    ];
    let fields = names.iter().map(|name| text_field(name)).collect();
    let leaves = names
        .iter()
        .map(|name| LeafData {
            values: Values::Bytes(vec![name.to_lowercase().into_bytes()]),
            def: Some(vec![1]),
            rep: None,
        })
        .collect();
    write_file(&path, fields, vec![leaves], AS_THE_WRITER);
    assert_eq!(
        lines(&tab(&path, false)),
        [
            "#Feature\tref\ttab_alt\ttab_Gene",
            "feature\ttab_ref\ttab_tab_alt\ttab_gene"
        ]
    );
    // The string columns a `--tab` output is read beside are named as that output names them.
    assert_eq!(
        text_columns(&path),
        ["Feature", "ref", "tab_alt", "tab_Gene"]
    );
}

/// The writer stores a tab column whose name an earlier column carries in other letters
/// behind `tab_` and states the tab columns' names in the footer: the header prints the
/// stated names, a BOOLEAN column prints the spelling of its stated name, and the tab
/// output of the same row read beside the part has the part's canonical form. A file that
/// states no such list, or a list of another length, names each column by its field, one
/// prefix before a key column's name taken off.
#[test]
fn a_column_stored_behind_the_prefix_for_its_letter_case_prints_under_the_stated_name() {
    let dir = scratch("stated");
    let plant = |path: &Path, tab_columns: Option<&str>| {
        let cells = ["A", "G", "a", "source", "no"];
        let mut fields: Vec<_> = ["ref", "alt", "Feature", "tab_REF", "canonical"]
            .iter()
            .map(|name| text_field(name))
            .collect();
        fields.push(value_field(
            "tab_CANONICAL",
            PhysicalType::BOOLEAN,
            Repetition::OPTIONAL,
            None,
        ));
        let mut leaves: Vec<LeafData> = cells
            .iter()
            .map(|cell| LeafData {
                values: Values::Bytes(vec![cell.as_bytes().to_vec()]),
                def: Some(vec![1]),
                rep: None,
            })
            .collect();
        leaves.push(LeafData {
            values: Values::Bool(vec![true]),
            def: Some(vec![1]),
            rep: None,
        });
        parquet_writer::write_file_stating(path, fields, vec![leaves], AS_THE_WRITER, tab_columns);
    };

    let stated = dir.join("stated.parquet");
    plant(&stated, Some("Feature\tREF\tcanonical\tCANONICAL"));
    let rows = ["#Feature\tREF\tcanonical\tCANONICAL", "a\tsource\tno\tYES"];
    assert_eq!(lines(&tab(&stated, false)), rows);
    assert_eq!(lines(&tab(&stated, true)), rows);
    let tab_file = dir.join("stated.tab");
    std::fs::write(&tab_file, format!("{}\n", rows.join("\n"))).expect("the tab output");
    assert_eq!(
        tab_beside_parquet(&tab_file, &stated).expect("the tab output is read"),
        tab(&stated, true)
    );
    // The string columns a `--tab` output is read beside are named as the footer states them.
    assert_eq!(text_columns(&stated), ["Feature", "REF", "canonical"]);

    for other in [
        None,
        Some("Feature\tREF\tcanonical"),
        Some("Feature\tREF\tcanonical\tCANONICAL\tTSL"),
    ] {
        let path = dir.join("unstated.parquet");
        plant(&path, other);
        assert_eq!(
            lines(&tab(&path, false)),
            [
                "#Feature\ttab_REF\tcanonical\ttab_CANONICAL",
                "a\tsource\tno\t1"
            ],
            "{other:?}"
        );
    }
}

#[test]
fn a_stored_type_the_cell_contract_lacks_is_refused_by_its_column() {
    let dir = scratch("unread");
    let path = dir.join("bytes.parquet");
    let fields = vec![
        text_field("Feature"),
        value_field("Raw", PhysicalType::BYTE_ARRAY, Repetition::OPTIONAL, None),
    ];
    let leaves = vec![
        LeafData {
            values: Values::Bytes(vec![b"a".to_vec()]),
            def: Some(vec![1]),
            rep: None,
        },
        LeafData {
            values: Values::Bytes(vec![b"b".to_vec()]),
            def: Some(vec![1]),
            rep: None,
        },
    ];
    write_file(&path, fields, vec![leaves], AS_THE_WRITER);
    assert_eq!(
        refusal(&path, false),
        format!(
            "{ADAPTER}{}: the column Raw holds BYTE_ARRAY values, which is not read",
            path.display()
        )
    );
}

/// A golden corpus through the writer: the canonical form of the engine's Parquet output
/// is that of its `--tab` output read beside it. `VEP_BINARY` names the built candidate
/// engine.
#[test]
#[ignore = "runs the engine VEP_BINARY names"]
fn a_golden_corpus_round_trips_through_the_writer() {
    let vep = PathBuf::from(std::env::var_os("VEP_BINARY").expect("VEP_BINARY names an engine"));
    for corpus in ["GRCh38", "GRCh38-everything", "GRCh38-minor-allele-columns"] {
        let directory = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/golden/116")
            .join(corpus);
        let manifest =
            std::fs::read_to_string(directory.join("manifest.json")).expect("a manifest");
        let manifest = vep_formats_compare::pyjson::loads(&manifest).expect("the manifest parses");
        let assembly = manifest
            .get("assembly")
            .and_then(|name| name.as_str())
            .expect("an assembly");
        let flags: Vec<&str> = manifest
            .get("flags")
            .and_then(|flags| flags.as_list())
            .map(|flags| flags.iter().filter_map(|flag| flag.as_str()).collect())
            .unwrap_or_default();
        let dir = scratch(&format!("round-trip-{corpus}"));
        let run = |options: &[&str], output: &Path| {
            let done = std::process::Command::new(&vep)
                .args(["--offline", "--json_cache"])
                .arg(directory.join("json_cache"))
                .args(["--species", "homo_sapiens", "--assembly", assembly])
                .args([
                    "--buffer_size",
                    "5000",
                    "--force_overwrite",
                    "--no_stats",
                    "--quiet",
                    "-i",
                ])
                .arg(directory.join("variants.vcf"))
                .args(&flags)
                .args(options)
                .arg("-o")
                .arg(output)
                .status()
                .expect("the engine runs");
            assert!(done.success(), "{corpus}");
        };
        let (tab_file, parquet) = (dir.join("out.tab"), dir.join("out.parquet"));
        run(&["--tab"], &tab_file);
        run(&["--output_format", "parquet"], &parquet);
        let from_parquet = tab(&parquet, true);
        let from_tab = tab_beside_parquet(&tab_file, &parquet).expect("the tab output is read");
        assert_eq!(from_parquet, from_tab, "{corpus}");
        let rows = data_lines(&std::fs::read_to_string(&tab_file).unwrap()).len();
        assert!(rows > 0);
        assert_eq!(lines(&from_parquet).len(), 1 + rows, "{corpus}");
    }
}

fn unhex(hex: &str) -> String {
    let bytes: Vec<u8> = (0..hex.len() / 2)
        .map(|at| u8::from_str_radix(&hex[2 * at..2 * at + 2], 16).expect("two hex digits"))
        .collect();
    String::from_utf8(bytes).expect("the hex of UTF-8 bytes")
}

fn hex(text: &str) -> String {
    let mut out = String::with_capacity(2 * text.len());
    for byte in text.bytes() {
        let _ = write!(out, "{byte:02x}");
    }
    out
}

/// What `parquet_as_tab` makes of a path, spelt as the case files spell it.
fn answer(path: &Path, canonical: bool) -> String {
    match parquet_as_tab(path, canonical) {
        Ok(text) => format!("K{}", hex(&text)),
        Err(Error::Comparator(text)) => format!("E{}", hex(&text)),
    }
}

/// The first lines two texts differ in, for a failure's message.
fn first_difference(found: &str, wanted: &str) -> String {
    let (found, wanted) = (unhex(&found[1..]), unhex(&wanted[1..]));
    let mut lines = found.split('\n').zip(wanted.split('\n')).enumerate();
    match lines.find(|(_, (found, wanted))| found != wanted) {
        Some((at, (found, wanted))) => {
            format!("line {}\n found {found:?}\n wants {wanted:?}", at + 1)
        }
        None => format!(
            "{} lines found, {} wanted",
            found.split('\n').count(),
            wanted.split('\n').count()
        ),
    }
}

/// The cases of `file` that `parquet_as_tab` answers otherwise than written, and the count
/// of all its cases.
fn differing(file: &Path) -> (usize, Vec<String>) {
    let cases = std::fs::read_to_string(file).expect("the cases are read");
    let mut differing = Vec::new();
    let mut count = 0;
    for (at, line) in cases.lines().enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        let [canonical, path, wanted] = fields[..] else {
            panic!("{}:{}: not three fields", file.display(), at + 1);
        };
        let path = unhex(path);
        let found = answer(Path::new(&path), canonical == "1");
        if found != wanted {
            differing.push(format!(
                "{}:{}: {path} canonical {canonical}: {} against {}: {}",
                file.display(),
                at + 1,
                &found[..1],
                &wanted[..1],
                first_difference(&found, wanted)
            ));
        }
        count += 1;
    }
    (count, differing)
}

/// The files of cases under the `parquet` folder of the directory
/// `VEP_FORMATS_COMPARE_PYDIFF_DIR` names, `unread.tsv` apart.
fn case_files() -> (Vec<PathBuf>, PathBuf) {
    let dir = std::env::var_os("VEP_FORMATS_COMPARE_PYDIFF_DIR")
        .expect("VEP_FORMATS_COMPARE_PYDIFF_DIR names the directory of the cases");
    let folder = PathBuf::from(dir).join("parquet");
    let mut files: Vec<PathBuf> = std::fs::read_dir(&folder)
        .expect("the folder of the cases is listed")
        .map(|entry| entry.expect("an entry of the folder").path())
        .filter(|path| path.extension().is_some_and(|kind| kind == "tsv"))
        .filter(|path| path.file_name().is_some_and(|name| name != "unread.tsv"))
        .collect();
    files.sort();
    (files, folder.join("unread.tsv"))
}

#[test]
#[ignore = "reads the cases under VEP_FORMATS_COMPARE_PYDIFF_DIR"]
fn every_planted_output_reads_as_its_case_is_written() {
    let (files, _) = case_files();
    assert!(!files.is_empty());
    let mut all = Vec::new();
    for file in files {
        let (count, differing) = differing(&file);
        println!(
            "{}: {count} cases, {} differ",
            file.display(),
            differing.len()
        );
        all.extend(differing);
    }
    assert!(all.is_empty(), "{}", all.join("\n"));
}

/// A stored type outside the writer's cell contract is refused by name, where the adapter
/// prints the interpreter's text of the value.
#[test]
#[ignore = "reads the cases under VEP_FORMATS_COMPARE_PYDIFF_DIR"]
fn a_stored_type_outside_the_cell_contract_is_refused() {
    let (_, unread) = case_files();
    let cases = std::fs::read_to_string(&unread).expect("the cases are read");
    for line in cases.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let path = unhex(fields[1]);
        let found = unhex(&answer(Path::new(&path), false)[1..]);
        println!(
            "{path}\n the tool: {found}\n the adapter: {:?}",
            unhex(&fields[2][1..])
        );
        assert!(found.ends_with(", which is not read"), "{found}");
    }
}

/// The time `parquet_as_tab` takes on the output `VEP_FORMATS_COMPARE_PARQUET` names, in
/// file order and in the canonical form; the two texts are kept as `file-order.tab` and
/// `canonical.tab` in the directory `VEP_FORMATS_COMPARE_PARQUET_KEEP` names.
#[test]
#[ignore = "reads the output VEP_FORMATS_COMPARE_PARQUET names"]
fn the_time_to_read_an_output() {
    let path = PathBuf::from(
        std::env::var_os("VEP_FORMATS_COMPARE_PARQUET")
            .expect("VEP_FORMATS_COMPARE_PARQUET names an output"),
    );
    if path.is_file() {
        println!("pages of the first column: {:?}", codec_of(&path));
    }
    for canonical in [false, true] {
        let started = std::time::Instant::now();
        let text = parquet_as_tab(&path, canonical).expect("the output is read");
        println!(
            "canonical {canonical}: {} bytes, {} lines in {:.2} s",
            text.len(),
            text.matches('\n').count(),
            started.elapsed().as_secs_f64()
        );
        if let Some(keep) = std::env::var_os("VEP_FORMATS_COMPARE_PARQUET_KEEP") {
            let name = if canonical {
                "canonical.tab"
            } else {
                "file-order.tab"
            };
            std::fs::write(PathBuf::from(keep).join(name), text).expect("the text is kept");
        }
    }
}
