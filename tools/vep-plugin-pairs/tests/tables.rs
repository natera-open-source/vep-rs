// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rules that read a record, a table or a cache: the tabix reader and the shard reader on
//! small files and on the committed table slices, and the classes of AlphaMissense, REVEL,
//! SpliceVault, LoFTEE, NearestExonJB, MaxEntScan and the mitochondrial contig on small
//! outputs. The expected values are those `tools/vep-formats-compare` gives on the same
//! inputs; the dbNSFP class of a mitochondrial record is handed the assembly and the
//! field's plugin, which that tool reads from a manifest.

use flate2::read::MultiGzDecoder;
use flate2::write::{DeflateEncoder, GzEncoder};
use flate2::{Compression, Crc};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;
use vep_plugin_pairs::pytext::Json;
use vep_plugin_pairs::rules::cache::{shards, TranscriptShards};
use vep_plugin_pairs::rules::entry::{parse_default, Entry, Row, ABSENT};
use vep_plugin_pairs::rules::maxentscan::SpliceWindows;
use vep_plugin_pairs::rules::nearest_exon::{nearest_exon_items, reordered_items, NearestExonItem};
use vep_plugin_pairs::rules::revel::{trimmed_allele_keys, RevelRow, RevelTable};
use vep_plugin_pairs::rules::shared::first_alternate;
use vep_plugin_pairs::rules::tabix::{IndexedTable, TableError, TableRow};
use vep_plugin_pairs::rules::{BuiltIn, Consult, RecordRows, Refusal, Rules};

const TABLE_HEADER: &str =
    "#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid";

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-tables-{}-{tag}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn file(&self, name: &str, content: &[u8]) -> PathBuf {
        let path = self.0.join(name);
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, content).unwrap();
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn gzip(data: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

/// One BGZF block of `data`: a gzip member whose `BC` extra field holds its size.
fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut encoder = DeflateEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).unwrap();
    let body = encoder.finish().unwrap();
    let mut crc = Crc::new();
    crc.update(data);
    let mut block = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00".to_vec();
    block.extend_from_slice(&(body.len() as u16 + 25).to_le_bytes());
    block.extend_from_slice(&body);
    block.extend_from_slice(&crc.sum().to_le_bytes());
    block.extend_from_slice(&(data.len() as u32).to_le_bytes());
    block
}

/// The bytes of a tabix index over `spans`, each a sequence name with the file offsets its
/// one chunk of bin 0 runs between, on the 1-based `position_column`.
fn tabix_index(spans: &[(String, u64, u64)], position_column: i32) -> Vec<u8> {
    let chunks: Vec<(String, u64, u64)> = spans
        .iter()
        .map(|(name, first, end)| (name.clone(), first << 16, end << 16))
        .collect();
    tabix_index_of_chunks(&chunks, position_column)
}

/// `tabix_index` over `chunks`, each a sequence name with the virtual offsets its one chunk
/// runs between: a block's file offset above the low 16 bits, a place in its text in them.
fn tabix_index_of_chunks(chunks: &[(String, u64, u64)], position_column: i32) -> Vec<u8> {
    let names: Vec<u8> = chunks
        .iter()
        .flat_map(|(name, _, _)| name.bytes().chain([0]))
        .collect();
    let mut index = b"TBI\x01".to_vec();
    for field in [
        chunks.len() as i32,
        0,
        1,
        position_column,
        position_column,
        i32::from(b'#'),
        0,
        names.len() as i32,
    ] {
        index.extend_from_slice(&field.to_le_bytes());
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
    index
}

/// A made-up table of `lines` under `header`, each line a BGZF block of its own, with a
/// tabix index on the 1-based `position_column` that files every row of a chromosome in
/// its first bin.
fn write_table(
    scratch: &Scratch,
    name: &str,
    header: &str,
    lines: &[&str],
    position_column: i32,
) -> PathBuf {
    let mut bytes = bgzf_block(format!("{header}\n").as_bytes());
    let mut spans: Vec<(String, u64, u64)> = Vec::new();
    for line in lines {
        let at = bytes.len() as u64;
        bytes.extend_from_slice(&bgzf_block(format!("{line}\n").as_bytes()));
        let chromosome = line.split('\t').next().unwrap();
        match spans.iter_mut().find(|(name, _, _)| name == chromosome) {
            Some(span) => span.2 = bytes.len() as u64,
            None => spans.push((chromosome.to_string(), at, bytes.len() as u64)),
        }
    }
    bytes.extend_from_slice(&bgzf_block(b""));
    scratch.file(
        &format!("{name}.tbi"),
        &gzip(&tabix_index(&spans, position_column)),
    );
    scratch.file(name, &bytes)
}

fn cells_of(rows: &[TableRow]) -> Vec<Vec<String>> {
    rows.iter().map(|row| row.cells.clone()).collect()
}

fn table_refusal(read: Result<IndexedTable, TableError>) -> String {
    match read {
        Ok(_) => panic!("the table was read"),
        Err(TableError(message)) => message,
    }
}

fn golden(corpus: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(corpus)
        .join("plugin_data/revel.tsv.gz")
}

#[test]
fn an_indexed_table_gives_the_rows_of_a_span_in_file_order() {
    for (corpus, column) in [("GRCh37-plugin-revel", 1), ("GRCh38-plugin-revel", 2)] {
        let path = golden(corpus);
        let table = IndexedTable::open(&path).unwrap();
        assert_eq!(table.header().len(), 9, "{corpus}");
        let mut text = String::new();
        MultiGzDecoder::new(std::fs::File::open(&path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        let mut held: BTreeMap<i128, Vec<Vec<String>>> = BTreeMap::new();
        for line in text.lines().filter(|line| !line.starts_with('#')) {
            let cells: Vec<String> = line.split('\t').map(str::to_string).collect();
            held.entry(cells[column].parse().unwrap())
                .or_default()
                .push(cells);
        }
        assert!(held.len() > 100, "{corpus}");
        let expected = |start: i128, end: i128| -> Vec<Vec<String>> {
            held.range(start..=end)
                .flat_map(|(_, rows)| rows.clone())
                .collect()
        };
        let positions: Vec<i128> = held.keys().copied().collect();
        for at in &positions {
            let rows = table.rows("21", *at, *at).unwrap();
            assert_eq!(cells_of(&rows), expected(*at, *at), "{corpus} {at}");
            assert!(rows.iter().all(|row| row.position == at.to_string()));
        }
        let last = positions[positions.len() - 1];
        for (start, end) in [
            (positions[0], last),
            (positions[3] - 40, positions[40] + 7),
            (last + 1, last + 9),
        ] {
            assert_eq!(
                cells_of(&table.rows("21", start, end).unwrap()),
                expected(start, end),
                "{corpus} {start} {end}"
            );
        }
        let at = positions[5];
        assert_eq!(
            table.rows("chr21", at, at).unwrap(),
            table.rows("21", at, at).unwrap()
        );
        assert!(table.rows("22", at, at).unwrap().is_empty());
        assert!(table.rows("21", at, at - 1).unwrap().is_empty());
    }
}

#[test]
fn a_table_row_ends_at_its_last_filled_cell_and_stands_at_the_indexed_column() {
    let scratch = Scratch::new("layouts");
    let nine = write_table(
        &scratch,
        "nine.tsv.gz",
        TABLE_HEADER,
        &[
            "21\t100\t500\tA\tG\tQ\tH\t0.742\tENST1;ENST2",
            "21\t100\t500\tA\tT\tQ\tL\t0.311",
            "21\t100\t500\tA",
            "21\t101\t501\tA\tG\tQ\tH\t\tENST1",
            "21\t102\t502\tA\tG\t\t\t\t",
            "22\t100\t500\tA\tG\tQ\tH\t0.5\tENST3",
            "chrM\t7\t9\tA\tG\tQ\tH\t0.1\tENST4",
        ],
        3,
    );
    let table = IndexedTable::open(&nine).unwrap();
    assert_eq!(table.path(), nine);
    assert_eq!(
        table.header(),
        [
            "chr",
            "hg19_pos",
            "grch38_pos",
            "ref",
            "alt",
            "aaref",
            "aaalt",
            "REVEL",
            "Ensembl_transcriptid"
        ]
    );
    let rows = table.rows("21", 500, 502).unwrap();
    assert_eq!(
        rows.iter().map(|row| row.cells.len()).collect::<Vec<_>>(),
        [9, 8, 4, 9, 5]
    );
    assert_eq!(
        rows.iter()
            .map(|row| row.position.as_str())
            .collect::<Vec<_>>(),
        ["500", "500", "500", "501", "502"]
    );
    assert_eq!(rows[0].cell(8), Some("ENST1;ENST2"));
    assert_eq!((rows[1].cell(7), rows[1].cell(8)), (Some("0.311"), None));
    assert_eq!((rows[2].cell(3), rows[2].cell(4)), (Some("A"), None));
    assert_eq!(
        (rows[3].cell(7), rows[3].cell(8)),
        (Some(""), Some("ENST1"))
    );
    assert_eq!(
        cells_of(&table.rows("21", 501, 501).unwrap()),
        [rows[3].cells.clone()]
    );
    assert_eq!(table.rows("chr22", 500, 500).unwrap().len(), 1);
    assert_eq!(table.rows("M", 9, 9).unwrap().len(), 1);
    assert_eq!(table.rows("chrM", 9, 9).unwrap().len(), 1);
    assert!(table.rows("21", 100, 102).unwrap().is_empty());
    assert!(table.rows("21", 0, 0).unwrap().is_empty());
    let first_column = write_table(
        &scratch,
        "second.tsv.gz",
        TABLE_HEADER,
        &["21\t100\t500\tA\tG\tQ\tH\t0.742\tENST1"],
        2,
    );
    let table = IndexedTable::open(&first_column).unwrap();
    assert_eq!(table.rows("21", 100, 100).unwrap()[0].position, "100");
    assert!(table.rows("21", 500, 500).unwrap().is_empty());
    let bare = write_table(&scratch, "bare.tsv.gz", "21\t5\t6", &["21\t7\t8"], 2);
    let table = IndexedTable::open(&bare).unwrap();
    assert!(table.header().is_empty());
    assert_eq!(table.rows("21", 7, 7).unwrap().len(), 1);
}

#[test]
fn a_table_without_a_readable_index_or_block_is_refused_with_its_reason() {
    let scratch = Scratch::new("refusals");
    let line = "21\t100\t500\tA\tG\tQ\tH\t0.742\tENST1";
    let table = write_table(&scratch, "table.tsv.gz", TABLE_HEADER, &[line], 3);
    let index_path = scratch.0.join("table.tsv.gz.tbi");
    let index = tabix_index(&[("21".to_string(), 0, 100)], 3);
    let reason = |index: &[u8]| {
        std::fs::write(&index_path, index).unwrap();
        table_refusal(IndexedTable::open(&table))
            .strip_prefix(&format!(
                "cannot read the tabix index {}: ",
                index_path.display()
            ))
            .unwrap()
            .to_string()
    };
    assert_eq!(reason(b"TBI\x01"), "Not a gzipped file (b'TB')");
    assert_eq!(reason(&gzip(b"")), "not a tabix index");
    assert_eq!(reason(b""), "not a tabix index");
    assert_eq!(reason(&gzip(b"BAI\x01 and more")), "not a tabix index");
    assert_eq!(
        reason(&gzip(&index[..20])),
        "unpack_from requires a buffer of at least 36 bytes for unpacking 32 bytes at offset 4 (actual buffer size is 20)"
    );
    assert_eq!(
        reason(&gzip(&index[..index.len() - 6])),
        format!(
            "unpack_from requires a buffer of at least {} bytes for unpacking 16 bytes at offset {} (actual buffer size is {})",
            index.len() - 4,
            index.len() - 20,
            index.len() - 6
        )
    );
    let packed = gzip(&index);
    assert_eq!(
        reason(&packed[..packed.len() - 9]),
        "Compressed file ended before the end-of-stream marker was reached"
    );
    let mut wrong_sum = packed.clone();
    let at = wrong_sum.len() - 8;
    wrong_sum[at] ^= 0xff;
    assert_eq!(reason(&wrong_sum), "CRC check failed");
    let mut wrong_length = packed.clone();
    let at = wrong_length.len() - 1;
    wrong_length[at] ^= 0xff;
    assert_eq!(reason(&wrong_length), "Incorrect length of data produced");
    let mut padded = packed.clone();
    padded.extend_from_slice(&[0, 0, 0]);
    padded.extend_from_slice(&gzip(b""));
    std::fs::write(&index_path, &padded).unwrap();
    assert_eq!(
        IndexedTable::open(&table)
            .unwrap()
            .rows("21", 500, 500)
            .unwrap()
            .len(),
        1
    );
    let mut negative_count = index.clone();
    let at = 36 + 3 + 4 + 4;
    negative_count[at..at + 4].copy_from_slice(&(-1i32).to_le_bytes());
    assert_eq!(reason(&gzip(&negative_count)), "bad char in struct format");
    let mut odd_name = tabix_index(&[("2\u{e9}".to_string(), 0, 100)], 3);
    odd_name[37] = 0xff;
    assert_eq!(
        reason(&gzip(&odd_name)),
        "'utf-8' codec can't decode byte 0xff in position 1: invalid start byte"
    );
    std::fs::remove_file(&index_path).unwrap();
    assert_eq!(
        table_refusal(IndexedTable::open(&table)),
        format!(
            "cannot read the tabix index of {table}: none of {table}.csi, {stem}.csi, {table}.tbi, {stem}.tbi exists",
            table = table.display(),
            stem = table.with_extension("").display()
        )
    );

    let plain = scratch.file(
        "plain.tsv.gz",
        &gzip(format!("{TABLE_HEADER}\n{line}\n").as_bytes()),
    );
    scratch.file("plain.tsv.gz.tbi", &gzip(&index));
    assert_eq!(
        table_refusal(IndexedTable::open(&plain)),
        format!("{} is not a BGZF file", plain.display())
    );
    let mut cut = bgzf_block(format!("{TABLE_HEADER}\n").as_bytes());
    cut.truncate(cut.len() - 12);
    let cut = scratch.file("cut.tsv.gz", &cut);
    scratch.file("cut.tsv.gz.tbi", &gzip(&index));
    assert_eq!(
        table_refusal(IndexedTable::open(&cut)),
        format!(
            "cannot read {}: Error -5 while decompressing data: incomplete or truncated stream",
            cut.display()
        )
    );
    let empty = scratch.file("empty.tsv.gz", b"");
    scratch.file("empty.tsv.gz.tbi", &gzip(&index));
    let table = IndexedTable::open(&empty).unwrap();
    assert!(table.header().is_empty() && table.rows("21", 500, 500).unwrap().is_empty());
}

/// Three BGZF streams one after the other, as compressed tables joined byte for byte are:
/// each is one block of text and closes with an empty block. The rows behind a closing
/// block are read as the rows before it are, under an index that numbers the block behind
/// an empty block by the empty block's offset, as htslib 1.9 writes one, and under an
/// index that numbers every block by its own offset.
#[test]
fn a_table_of_concatenated_streams_gives_the_rows_behind_the_first_stream() {
    let scratch = Scratch::new("streams");
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
    let first = stream(text(&[TABLE_HEADER, lines[0], lines[1]]));
    let second = stream(text(&lines[2..5]));
    let third = stream(text(&lines[5..]));
    let table = scratch.file(
        "streams.tsv.gz",
        &[first.as_slice(), second.as_slice(), third.as_slice()].concat(),
    );
    // The file offsets of the first stream's closing block, of the second stream's text and
    // of the file's last block; then the place in the first stream's text where its rows
    // begin and the place in the second stream's where its rows of 21 end.
    let closing = (first.len() - closing_block.len()) as u64;
    let behind = first.len() as u64;
    let last = (first.len() + second.len() + third.len() - closing_block.len()) as u64;
    let after_header = TABLE_HEADER.len() as u64 + 1;
    let after_21 = text(&lines[2..4]).len() as u64;
    let cells = |lines: &[&str]| -> Vec<Vec<String>> {
        lines
            .iter()
            .map(|line| line.split('\t').map(str::to_string).collect())
            .collect()
    };
    for (numbered, block) in [("by the empty block", closing), ("by its own", behind)] {
        let split = (block << 16) | after_21;
        let chunks = [
            ("21".to_string(), after_header, split),
            ("22".to_string(), split, last << 16),
        ];
        scratch.file(
            "streams.tsv.gz.tbi",
            &gzip(&tabix_index_of_chunks(&chunks, 3)),
        );
        let read = IndexedTable::open(&table).unwrap();
        assert_eq!(read.header().len(), 9, "{numbered}");
        for (chromosome, start, end, expected) in [
            ("21", 500, 503, &lines[..4]),
            ("21", 502, 503, &lines[2..4]),
            ("22", 504, 505, &lines[4..]),
            ("22", 505, 505, &lines[5..]),
        ] {
            assert_eq!(
                cells_of(&read.rows(chromosome, start, end).unwrap()),
                cells(expected),
                "{numbered}: {chromosome} {start} to {end}"
            );
        }
    }
}

/// A chunk that runs over the first stream's closing block and ends in the text behind it
/// gives the rows before its end and no row after it, under an index that numbers that
/// text by the closing block's offset and under an index that numbers it by its own.
#[test]
fn a_chunk_that_ends_in_the_text_behind_an_empty_block_gives_no_row_past_its_end() {
    let scratch = Scratch::new("chunk-end");
    let lines = [
        "21\t100\t500\tA\tG\tQ\tH\t0.1\tENST1",
        "21\t101\t501\tA\tG\tQ\tH\t0.2\tENST1",
        "21\t102\t502\tA\tG\tQ\tH\t0.3\tENST1",
    ];
    let text =
        |lines: &[&str]| -> String { lines.iter().map(|line| format!("{line}\n")).collect() };
    let closing_block = bgzf_block(b"");
    let stream = |text: String| [bgzf_block(text.as_bytes()), closing_block.clone()].concat();
    let first = stream(text(&[TABLE_HEADER, lines[0]]));
    let second = stream(text(&lines[1..]));
    let table = scratch.file(
        "end.tsv.gz",
        &[first.as_slice(), second.as_slice()].concat(),
    );
    // The file offsets of the first stream's closing block and of the second stream's
    // text; then the place in the first stream's text where its row begins and the place
    // in the second stream's where its first row ends.
    let closing = (first.len() - closing_block.len()) as u64;
    let behind = first.len() as u64;
    let after_header = TABLE_HEADER.len() as u64 + 1;
    let after_first = text(&lines[1..2]).len() as u64;
    let before_the_end: Vec<Vec<String>> = lines[..2]
        .iter()
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    for (numbered, block) in [("by the empty block", closing), ("by its own", behind)] {
        let chunks = [("21".to_string(), after_header, (block << 16) | after_first)];
        scratch.file("end.tsv.gz.tbi", &gzip(&tabix_index_of_chunks(&chunks, 3)));
        let read = IndexedTable::open(&table).unwrap();
        assert_eq!(
            cells_of(&read.rows("21", 500, 502).unwrap()),
            before_the_end,
            "{numbered}"
        );
    }
}

/// A line goes on in the next block that holds text: the line the third block leaves open
/// takes the text behind the empty block, and the end of the file ends the lines.
#[test]
fn a_line_split_over_blocks_is_read_whole_and_an_empty_block_before_the_last_is_passed_over() {
    let scratch = Scratch::new("blocks");
    let mut bytes = bgzf_block(b"#a\tb\tc\n21\t5");
    bytes.extend_from_slice(&bgzf_block(b"00\tx\n21\t6"));
    bytes.extend_from_slice(&bgzf_block(b"00\ty\n21\t700\tz"));
    let end = bytes.len() as u64;
    bytes.extend_from_slice(&bgzf_block(b""));
    bytes.extend_from_slice(&bgzf_block(b"21\t800\tpast\n"));
    let table = scratch.file("split.tsv.gz", &bytes);
    scratch.file(
        "split.tsv.gz.tbi",
        &gzip(&tabix_index(&[("21".to_string(), 0, end + 200)], 2)),
    );
    let table = IndexedTable::open(&table).unwrap();
    assert_eq!(table.header(), ["a", "b", "c"]);
    assert_eq!(
        cells_of(&table.rows("21", 1, 1000).unwrap()),
        [
            vec!["21", "500", "x"],
            vec!["21", "600", "y"],
            vec!["21", "700", "z21", "800", "past"]
        ]
    );
    assert_eq!(
        cells_of(&table.rows("21", 600, 600).unwrap()),
        [["21", "600", "y"]]
    );
}

/// The made-up table of the fixture directory, of three sequence names: a row's second
/// column is the position `tabix` indexed it at and its third the row's end, so that some
/// rows lie over the edge of a small bin and are filed in a wider one.
const SYNTHETIC_TABLE: &str = "revel.synthetic.tsv.gz";
/// The index `tabix -s 1 -b 2 -e 3` wrote of the made-up table.
const SYNTHETIC_TBI: &str = "revel.synthetic.tsv.gz.tbi";
/// The index `tabix -C` wrote of it on the same columns.
const SYNTHETIC_CSI: &str = "revel.synthetic.tsv.gz.csi";
/// The index `tabix -C -m 12` wrote of a copy of it: smallest bins of 4 kb, one level more.
const SYNTHETIC_CSI_SHIFT_12: &str = "revel.synthetic.min_shift_12.tsv.gz.csi";

/// The bytes of a file of the tool's fixture directory.
fn fixture(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name);
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// A copy of the made-up table in the directory `directory` of `scratch`, with each of
/// `beside`, a file name and its bytes, in the same directory.
fn synthetic_table(scratch: &Scratch, directory: &str, beside: &[(&str, &[u8])]) -> PathBuf {
    for (name, content) in beside {
        scratch.file(&format!("{directory}/{name}"), content);
    }
    scratch.file(
        &format!("{directory}/{SYNTHETIC_TABLE}"),
        &fixture(SYNTHETIC_TABLE),
    )
}

fn gunzip(packed: &[u8]) -> Vec<u8> {
    let mut data = Vec::new();
    MultiGzDecoder::new(packed).read_to_end(&mut data).unwrap();
    data
}

/// The bytes of a CSI index of `min_shift` and `depth` over `spans`, each a sequence name
/// with the file offsets its one chunk of bin 0 runs between, on the 1-based
/// `position_column`; after bin 0 each sequence has the bin of its counts.
fn csi_index(
    spans: &[(String, u64, u64)],
    position_column: i32,
    min_shift: i32,
    depth: u32,
) -> Vec<u8> {
    let names: Vec<u8> = spans
        .iter()
        .flat_map(|(name, _, _)| name.bytes().chain([0]))
        .collect();
    let mut index = b"CSI\x01".to_vec();
    for field in [
        min_shift,
        depth as i32,
        28 + names.len() as i32,
        0,
        1,
        position_column,
        position_column,
        i32::from(b'#'),
        0,
        names.len() as i32,
    ] {
        index.extend_from_slice(&field.to_le_bytes());
    }
    index.extend_from_slice(&names);
    index.extend_from_slice(&(spans.len() as i32).to_le_bytes());
    let counts = (((1u64 << (3 * (depth + 1))) - 1) / 7 + 1) as u32;
    for (_, first, end) in spans {
        let (first, end) = (first << 16, end << 16);
        index.extend_from_slice(&2i32.to_le_bytes());
        index.extend_from_slice(&0u32.to_le_bytes());
        index.extend_from_slice(&first.to_le_bytes());
        index.extend_from_slice(&1i32.to_le_bytes());
        index.extend_from_slice(&first.to_le_bytes());
        index.extend_from_slice(&end.to_le_bytes());
        index.extend_from_slice(&counts.to_le_bytes());
        index.extend_from_slice(&0u64.to_le_bytes());
        index.extend_from_slice(&2i32.to_le_bytes());
        for word in [first, end, 1, 0] {
            index.extend_from_slice(&word.to_le_bytes());
        }
    }
    index
}

/// The made-up table opened under each of its three indexes, each copy with its index in
/// a directory of `scratch` of the layout's name.
fn synthetic_tables(scratch: &Scratch) -> [(&'static str, IndexedTable); 3] {
    [
        ("tbi", SYNTHETIC_TBI, "tbi"),
        ("csi", SYNTHETIC_CSI, "csi"),
        ("csi-shift-12", SYNTHETIC_CSI_SHIFT_12, "csi"),
    ]
    .map(|(layout, index, extension)| {
        let beside = format!("{SYNTHETIC_TABLE}.{extension}");
        let index = fixture(index);
        let table = synthetic_table(scratch, layout, &[(beside.as_str(), index.as_slice())]);
        (layout, IndexedTable::open(&table).unwrap())
    })
}

#[test]
fn a_table_under_a_csi_index_gives_the_rows_its_tbi_index_gives() {
    let scratch = Scratch::new("layouts-alike");
    let tables = synthetic_tables(&scratch);
    let text = String::from_utf8(gunzip(&fixture(SYNTHETIC_TABLE))).unwrap();
    let lines: Vec<Vec<String>> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    let at = |cells: &[String]| -> i128 { cells[1].parse().unwrap() };
    let expected = |name: &str, start: i128, end: i128| -> Vec<Vec<String>> {
        lines
            .iter()
            .filter(|cells| cells[0] == name && start <= at(cells) && at(cells) <= end)
            .cloned()
            .collect()
    };
    let (mut hits, mut misses) = (0, 0);
    for name in ["21", "22", "X"] {
        let positions: Vec<i128> = lines
            .iter()
            .filter(|cells| cells[0] == name)
            .map(|cells| at(cells))
            .collect();
        let (first, last) = (positions[0], positions[positions.len() - 1]);
        let mut spans = vec![
            (first, first),
            (last, last),
            (first, last),
            (1, 1 << 40),
            (last + 1, last + 9),
        ];
        for pair in positions.windows(2) {
            spans.extend([
                (pair[0], pair[0]),
                (pair[0], pair[1]),
                (pair[0] + 1, pair[1] - 1),
            ]);
        }
        for (start, end) in spans {
            let rows = expected(name, start, end);
            for (layout, table) in &tables {
                assert_eq!(
                    cells_of(&table.rows(name, start, end).unwrap()),
                    rows,
                    "{layout} {name} {start} {end}"
                );
            }
            if rows.is_empty() {
                misses += 1;
            } else {
                hits += 1;
            }
        }
    }
    assert_eq!(lines.len(), 25);
    assert!(hits >= 50 && misses >= 20, "{hits} {misses}");
    for (layout, table) in &tables {
        assert_eq!(table.header().len(), 9, "{layout}");
        assert_eq!(table.rows("21", 100, 100).unwrap().len(), 2, "{layout}");
        assert_eq!(
            table.rows("22", 16_384, 16_385).unwrap().len(),
            2,
            "{layout}"
        );
        assert_eq!(
            table.rows("chrX", 7, 7).unwrap(),
            table.rows("X", 7, 7).unwrap(),
            "{layout}"
        );
        assert_eq!(table.rows("X", 7, 7).unwrap().len(), 2, "{layout}");
        // A row is found at the position it begins at, not at one inside its span.
        assert_eq!(
            table.rows("X", 40_000, 40_000).unwrap().len(),
            1,
            "{layout}"
        );
        assert!(
            table.rows("X", 50_000, 50_000).unwrap().is_empty(),
            "{layout}"
        );
        assert!(table.rows("7", 100, 100).unwrap().is_empty(), "{layout}");
    }
}

#[test]
fn a_span_from_zero_or_below_is_read_from_the_first_position() {
    let scratch = Scratch::new("from-zero");
    for (layout, table) in &synthetic_tables(&scratch) {
        assert_eq!(table.rows("21", 0, 100).unwrap().len(), 2, "{layout}");
        assert_eq!(table.rows("X", 0, 7).unwrap().len(), 2, "{layout}");
        assert_eq!(
            cells_of(&table.rows("21", -5, 5000).unwrap()),
            cells_of(&table.rows("21", 1, 5000).unwrap()),
            "{layout}"
        );
        assert_eq!(table.rows("21", 1, 5000).unwrap().len(), 4, "{layout}");
        assert_eq!(
            table.rows("22", -1_000_000_000_000, 50).unwrap().len(),
            1,
            "{layout}"
        );
        assert!(table.rows("21", 0, 0).unwrap().is_empty(), "{layout}");
        assert!(table.rows("21", -9, -3).unwrap().is_empty(), "{layout}");
    }
}

#[test]
fn a_csi_index_that_is_cut_short_or_is_not_one_is_refused_with_its_path() {
    let scratch = Scratch::new("csi-refusals");
    let table = synthetic_table(&scratch, "refused", &[]);
    let index_path = PathBuf::from(format!("{}.csi", table.display()));
    let reason = |index: &[u8]| {
        std::fs::write(&index_path, index).unwrap();
        table_refusal(IndexedTable::open(&table))
            .strip_prefix(&format!(
                "cannot read the tabix index {}: ",
                index_path.display()
            ))
            .unwrap()
            .to_string()
    };
    let packed = fixture(SYNTHETIC_CSI);
    let raw = gunzip(&packed);
    let cut_short = "Compressed file ended before the end-of-stream marker was reached";
    assert_eq!(reason(&packed[..packed.len() - 40]), cut_short);
    assert_eq!(reason(&packed[..3]), cut_short);
    assert_eq!(
        reason(&gzip(&raw[..10])),
        "unpack_from requires a buffer of at least 16 bytes for unpacking 12 bytes at offset 4 (actual buffer size is 10)"
    );
    assert_eq!(
        reason(&gzip(&raw[..raw.len() - 20])),
        format!(
            "unpack_from requires a buffer of at least {} bytes for unpacking 32 bytes at offset {} (actual buffer size is {})",
            raw.len() - 8,
            raw.len() - 40,
            raw.len() - 20
        )
    );
    let mut other_magic = raw.clone();
    other_magic[3] = 2;
    assert_eq!(reason(&gzip(&other_magic)), "not a tabix index");
    assert_eq!(reason(b""), "not a tabix index");
    assert_eq!(reason(b"CSI\x01"), "Not a gzipped file (b'CS')");
    // The header's shift, depth and length of auxiliary data stand at bytes 4, 8 and 12,
    // the length of the names at 40, the first sequence's count of bins at 56 and its
    // first bin's count of chunks at 72.
    let with = |at: usize, value: i32| {
        let mut index = raw.clone();
        index[at..at + 4].copy_from_slice(&value.to_le_bytes());
        gzip(&index)
    };
    assert_eq!(
        reason(&with(8, 11)),
        "a CSI index of shift 14 and depth 11 is not read"
    );
    assert_eq!(
        reason(&with(8, -1)),
        "a CSI index of shift 14 and depth -1 is not read"
    );
    assert_eq!(
        reason(&with(4, -1)),
        "a CSI index of shift -1 and depth 8 is not read"
    );
    assert_eq!(
        reason(&with(4, 39)),
        "a CSI index of shift 39 and depth 8 is not read"
    );
    assert_eq!(reason(&with(12, 0)), "a CSI index without a tabix header");
    assert_eq!(reason(&with(12, 27)), "a CSI index without a tabix header");
    assert_eq!(reason(&with(12, -1)), "bad char in struct format");
    assert_eq!(
        reason(&with(40, 9)),
        "unpack_from requires a buffer of at least 37 bytes for unpacking 9 bytes at offset 28 (actual buffer size is 36)"
    );
    assert_eq!(reason(&with(56, -1)), "bad char in struct format");
    assert_eq!(reason(&with(72, -1)), "bad char in struct format");
    std::fs::write(&index_path, &packed).unwrap();
    assert_eq!(
        IndexedTable::open(&table)
            .unwrap()
            .rows("21", 100, 100)
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn an_index_is_looked_for_under_four_names_in_the_order_tabix_looks() {
    let scratch = Scratch::new("index-names");
    let (tbi, csi) = (fixture(SYNTHETIC_TBI), fixture(SYNTHETIC_CSI));
    let (tbi, csi, other): (&[u8], &[u8], &[u8]) = (&tbi, &csi, b"no index");
    let (file_csi, stem_csi, file_tbi, stem_tbi) = (
        "revel.synthetic.tsv.gz.csi",
        "revel.synthetic.tsv.csi",
        "revel.synthetic.tsv.gz.tbi",
        "revel.synthetic.tsv.tbi",
    );
    let rows = |table: &Path| {
        IndexedTable::open(table)
            .unwrap()
            .rows("22", 1, 20_000)
            .unwrap()
            .len()
    };
    let refused_at = |table: &Path, name: &str| {
        assert_eq!(
            table_refusal(IndexedTable::open(table)),
            format!(
                "cannot read the tabix index {}: Not a gzipped file (b'no')",
                table.with_file_name(name).display()
            )
        );
    };

    let every = [
        (file_csi, other),
        (stem_csi, csi),
        (file_tbi, tbi),
        (stem_tbi, tbi),
    ];
    refused_at(&synthetic_table(&scratch, "first", &every), file_csi);
    let first_alone = [
        (file_csi, csi),
        (stem_csi, other),
        (file_tbi, other),
        (stem_tbi, other),
    ];
    assert_eq!(
        rows(&synthetic_table(&scratch, "first-read", &first_alone)),
        3
    );

    let from_second = [(stem_csi, other), (file_tbi, tbi), (stem_tbi, tbi)];
    refused_at(&synthetic_table(&scratch, "second", &from_second), stem_csi);
    let second_alone = [(stem_csi, csi), (file_tbi, other), (stem_tbi, other)];
    assert_eq!(
        rows(&synthetic_table(&scratch, "second-read", &second_alone)),
        3
    );

    let from_third = [(file_tbi, other), (stem_tbi, tbi)];
    refused_at(&synthetic_table(&scratch, "third", &from_third), file_tbi);
    let third_alone = [(file_tbi, tbi), (stem_tbi, other)];
    assert_eq!(
        rows(&synthetic_table(&scratch, "third-read", &third_alone)),
        3
    );

    assert_eq!(
        rows(&synthetic_table(&scratch, "fourth", &[(stem_tbi, tbi)])),
        3
    );
    refused_at(
        &synthetic_table(&scratch, "fourth-unread", &[(stem_tbi, other)]),
        stem_tbi,
    );

    let none = synthetic_table(&scratch, "none", &[]);
    let beside = |name: &str| none.with_file_name(name).display().to_string();
    assert_eq!(
        table_refusal(IndexedTable::open(&none)),
        format!(
            "cannot read the tabix index of {}: none of {}, {}, {}, {} exists",
            none.display(),
            beside(file_csi),
            beside(stem_csi),
            beside(file_tbi),
            beside(stem_tbi)
        )
    );

    // A file name without an extension has no stem of its own to look under.
    let bare = scratch.file("bare/synthetic", &fixture(SYNTHETIC_TABLE));
    scratch.file("bare/synthetic.tsv.csi", csi);
    assert_eq!(
        table_refusal(IndexedTable::open(&bare)),
        format!(
            "cannot read the tabix index of {bare}: none of {bare}.csi, {bare}.tbi exists",
            bare = bare.display()
        )
    );
    scratch.file("bare/synthetic.tbi", tbi);
    assert_eq!(rows(&bare), 3);

    // An index is read in the layout its bytes name, whatever its file is called.
    let renamed = [(file_csi, tbi)];
    assert_eq!(rows(&synthetic_table(&scratch, "tbi-as-csi", &renamed)), 3);
    let renamed = [(file_tbi, csi)];
    assert_eq!(rows(&synthetic_table(&scratch, "csi-as-tbi", &renamed)), 3);
}

#[test]
fn a_span_ends_at_the_last_position_a_csi_index_holds() {
    let scratch = Scratch::new("csi-last-position");
    let mut bytes = bgzf_block(b"#name\tposition\tvalue\n");
    let first = bytes.len() as u64;
    bytes.extend_from_slice(&bgzf_block(
        b"21\t100\tone\n21\t16384\ttwo\n21\t16385\tthree\n21\t20000\tfour\n",
    ));
    let spans = [("21".to_string(), first, bytes.len() as u64)];
    bytes.extend_from_slice(&bgzf_block(b""));
    let positions = |table: &IndexedTable, start: i128, end: i128| -> Vec<String> {
        table
            .rows("21", start, end)
            .unwrap()
            .iter()
            .map(|row| row.position.clone())
            .collect()
    };

    // Bin 0 alone, of 16 kb: the last position the index holds is 16384.
    let shallow = scratch.file("shallow/table.tsv.gz", &bytes);
    scratch.file(
        "shallow/table.tsv.gz.csi",
        &gzip(&csi_index(&spans, 2, 14, 0)),
    );
    let shallow = IndexedTable::open(&shallow).unwrap();
    assert_eq!(shallow.header(), ["name", "position", "value"]);
    assert_eq!(positions(&shallow, 1, 16_384), ["100", "16384"]);
    assert_eq!(positions(&shallow, 16_384, 30_000), ["16384"]);
    assert!(positions(&shallow, 16_385, 30_000).is_empty());
    assert!(positions(&shallow, 20_000, 20_000).is_empty());

    // One level under bin 0: the last position is 131072, and the same lines are read.
    let deeper = scratch.file("deeper/table.tsv.gz", &bytes);
    scratch.file(
        "deeper/table.tsv.gz.csi",
        &gzip(&csi_index(&spans, 2, 14, 1)),
    );
    let deeper = IndexedTable::open(&deeper).unwrap();
    assert_eq!(positions(&deeper, 16_385, 30_000), ["16385", "20000"]);
    assert_eq!(
        positions(&deeper, 1, 1 << 40),
        ["100", "16384", "16385", "20000"]
    );

    // Under a TBI index the span is read as given.
    let tbi = scratch.file("tbi/table.tsv.gz", &bytes);
    scratch.file("tbi/table.tsv.gz.tbi", &gzip(&tabix_index(&spans, 2)));
    let tbi = IndexedTable::open(&tbi).unwrap();
    assert_eq!(positions(&tbi, 16_385, 30_000), ["16385", "20000"]);
}

fn refusal<T>(read: Result<T, Refusal>) -> String {
    match read {
        Ok(_) => panic!("the read was taken"),
        Err(Refusal(message)) => message,
    }
}

fn stable_ids(shards: &TranscriptShards, chromosome: Option<&str>) -> Result<Vec<String>, Refusal> {
    let mut ids = Vec::new();
    let mut note = |transcript: Json| {
        if let Some(Json::Str(id)) = transcript.get("stable_id") {
            ids.push(id.clone());
        }
        Ok(())
    };
    match chromosome {
        Some(chromosome) => shards.of_chromosome(chromosome, &mut note)?,
        None => shards.of_every_chromosome(&mut note)?,
    }
    Ok(ids)
}

#[test]
fn a_cache_yields_its_transcripts_shard_by_shard_in_reading_order() {
    let scratch = Scratch::new("cache");
    assert_eq!(
        refusal(TranscriptShards::new(&scratch.0).map(|_| ())),
        format!(
            "--json-cache: {} has no transcripts/ directory",
            scratch.0.display()
        )
    );
    scratch.file(
        "transcripts/21/9-10.json",
        br#"[{"stable_id": "T9"}, null, {"stable_id": "T9b"}]"#,
    );
    scratch.file("transcripts/21/10-20.json", br#"[{"stable_id": "T10"}]"#);
    scratch.file(
        "transcripts/21/1-8.json.gz",
        &gzip(br#"{"1-4": [{"stable_id": "G1"}], "5-8": [null, {"stable_id": "G5"}]}"#),
    );
    scratch.file("transcripts/21/.hidden.json", br#"[{"stable_id": "H"}]"#);
    scratch.file("transcripts/21/notes.txt", b"no shard");
    scratch.file("transcripts/21/upper.JSON", b"no shard");
    scratch.file("transcripts/3/1-2.json", br#"[{"stable_id": "T3"}]"#);
    scratch.file("transcripts/MT/1-2.json", br#"[{"stable_id": "TM"}]"#);
    scratch.file("transcripts/info.json", br#"[{"stable_id": "none"}]"#);
    let cache = TranscriptShards::new(&scratch.0).unwrap();
    assert_eq!(cache.root(), scratch.0.join("transcripts"));
    assert_eq!(
        cache.chromosomes().iter().collect::<Vec<_>>(),
        ["21", "3", "MT"]
    );
    let names: Vec<String> = shards(&scratch.0.join("transcripts/21"))
        .iter()
        .map(|shard| shard.file_name().unwrap().to_str().unwrap().to_string())
        .collect();
    assert_eq!(
        names,
        [".hidden.json", "10-20.json", "9-10.json", "1-8.json.gz"]
    );
    assert_eq!(
        stable_ids(&cache, Some("21")).unwrap(),
        ["H", "T10", "T9", "T9b", "G1", "G5"]
    );
    assert_eq!(stable_ids(&cache, Some("MT")).unwrap(), ["TM"]);
    assert!(stable_ids(&cache, Some("22")).unwrap().is_empty());
    assert_eq!(
        stable_ids(&cache, None).unwrap(),
        ["H", "T10", "T9", "T9b", "G1", "G5", "T3", "TM"]
    );
    let mut seen = 0;
    let stopped = cache.of_every_chromosome(|_| {
        seen += 1;
        if seen == 2 {
            return Err(Refusal("enough".to_string()));
        }
        Ok(())
    });
    assert_eq!((refusal(stopped), seen), ("enough".to_string(), 2));

    let broken = scratch.file("transcripts/3/0-1.json", b"[{");
    assert_eq!(
        refusal(stable_ids(&cache, Some("3"))).split(": ").next(),
        Some(format!("{} is not JSON", broken.display()).as_str())
    );
    std::fs::write(&broken, br#"{"1-2": {"stable_id": "T"}}"#).unwrap();
    assert_eq!(
        refusal(stable_ids(&cache, None)),
        format!(
            "{}: a range-keyed map value is not a transcript array",
            broken.display()
        )
    );
    std::fs::write(&broken, b"[3]").unwrap();
    assert_eq!(
        refusal(stable_ids(&cache, Some("3"))),
        format!("{}: a transcript is not an object", broken.display())
    );
    std::fs::remove_file(&broken).unwrap();
    std::fs::create_dir(scratch.0.join("transcripts/3/dir.json")).unwrap();
    assert_eq!(
        refusal(stable_ids(&cache, Some("3"))),
        format!(
            "cannot read {}: Is a directory",
            scratch.0.join("transcripts/3/dir.json").display()
        )
    );
}

#[test]
fn a_cache_read_on_threads_yields_the_same_transcripts_and_the_same_refusal() {
    let scratch = Scratch::new("cache-threads");
    // Forty shards of three transcripts each on one chromosome, a few on another.
    for shard in 0..40 {
        let text = format!(
            r#"[{{"stable_id": "T{shard}a"}}, null, {{"stable_id": "T{shard}b"}}, {{"stable_id": "T{shard}c"}}]"#
        );
        scratch.file(
            &format!("transcripts/21/{:03}-{:03}.json", shard, shard + 1),
            text.as_bytes(),
        );
    }
    scratch.file("transcripts/3/1-2.json", br#"[{"stable_id": "T3"}]"#);
    let cache = TranscriptShards::new(&scratch.0).unwrap();
    let expected: Vec<String> = (0..40)
        .flat_map(|shard| ["a", "b", "c"].map(|part| format!("T{shard}{part}")))
        .collect();
    let every: Vec<String> = expected.iter().cloned().chain(["T3".to_string()]).collect();
    // What the reading threads make of each transcript reaches the visitor in the same
    // order as the transcripts themselves.
    let read_ids = |chromosome: Option<&str>| -> Result<Vec<String>, Refusal> {
        let read = |transcript: Json| match transcript.get("stable_id") {
            Some(Json::Str(id)) if !id.ends_with('b') => Some(id.clone()),
            _ => None,
        };
        let mut ids = Vec::new();
        let note = |id: String| {
            ids.push(id);
            Ok(())
        };
        match chromosome {
            Some(chromosome) => cache.read_of_chromosome(chromosome, read, note)?,
            None => cache.read_of_every_chromosome(read, note)?,
        }
        Ok(ids)
    };
    let without_b: Vec<String> = every
        .iter()
        .filter(|id| !id.ends_with('b'))
        .cloned()
        .collect();
    for threads in [1, 2, 5, 16, 64] {
        cache.read_on(threads);
        assert_eq!(
            stable_ids(&cache, Some("21")).unwrap(),
            expected,
            "{threads}"
        );
        assert_eq!(stable_ids(&cache, None).unwrap(), every, "{threads}");
        assert_eq!(read_ids(None).unwrap(), without_b, "{threads}");
        assert_eq!(read_ids(Some("3")).unwrap(), ["T3"], "{threads}");
        // A visitor's refusal ends the walk at the transcript it refuses.
        let mut seen = Vec::new();
        let stopped = cache.of_chromosome("21", |transcript| {
            seen.push(transcript.get("stable_id").cloned());
            match seen.len() {
                50 => Err(Refusal("enough".to_string())),
                _ => Ok(()),
            }
        });
        assert_eq!((refusal(stopped), seen.len()), ("enough".to_string(), 50));
    }
    // Of two shards that cannot be read the first in reading order is the one refused,
    // after every transcript of the shards before it.
    let first = scratch.file("transcripts/21/017-018.json", b"[{");
    scratch.file("transcripts/21/031-032.json", b"[3]");
    for threads in [1, 2, 5, 16, 64] {
        cache.read_on(threads);
        let mut seen = Vec::new();
        let walked = cache.of_chromosome("21", |transcript| {
            if let Some(Json::Str(id)) = transcript.get("stable_id") {
                seen.push(id.clone());
            }
            Ok(())
        });
        assert_eq!(
            refusal(walked).split(": ").next(),
            Some(format!("{} is not JSON", first.display()).as_str()),
            "{threads}"
        );
        assert_eq!(seen, expected[..17 * 3], "{threads}");
        assert_eq!(
            refusal(read_ids(None)).split(": ").next(),
            Some(format!("{} is not JSON", first.display()).as_str()),
            "{threads}"
        );
    }
}

const HEADER_14: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";

/// One default-format row of a transcript, `row` its name, Location, allele, Feature,
/// Consequence and Extra text, with its cDNA position `cdna` and its `Amino_acids`.
fn output_row(row: [&str; 6], cdna: &str, amino_acids: &str) -> String {
    let [name, location, allele, feature, consequence, extra] = row;
    format!("{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t{consequence}\t{cdna}\t-\t-\t{amino_acids}\t-\t-\t{extra}")
}

fn default_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    consequence: &str,
    extra: &str,
) -> String {
    output_row(
        [name, location, allele, feature, consequence, extra],
        "-",
        "-",
    )
}

fn default_text<S: AsRef<str>>(rows: &[S]) -> String {
    let mut text = HEADER_14.to_string();
    for row in rows {
        text.push('\n');
        text.push_str(row.as_ref());
    }
    text + "\n"
}

/// What a class is handed beside the two rows: the engine's row or none, the readers, the
/// assembly the run names and the plugin the run gives the field.
#[derive(Default, Clone, Copy)]
struct Handed<'a> {
    revel_table: Option<&'a RevelTable>,
    splice_windows: Option<&'a SpliceWindows>,
    assembly: Option<&'a str>,
    source: Option<&'a str>,
}

/// The class's answer on the field `name` of the row `e`: `values` are the reference's and
/// the engine's, `rows` the engine's row, the reference's rows of the record and the
/// engine's.
fn consult(
    klass: &str,
    e: &Row,
    name: &str,
    values: (&str, &str),
    rows: (Option<&Row>, Option<&[Row]>, Option<&[Row]>),
    handed: Handed,
) -> bool {
    let (expected, actual) = values;
    let (a, record_rows, candidate_rows) = rows;
    let declared = BTreeSet::new();
    BuiltIn
        .accepts(&Consult {
            klass,
            e,
            name,
            expected,
            actual,
            a,
            record_rows,
            candidate_rows,
            declared: &declared,
            assembly: handed.assembly,
            source: handed.source,
            gwas_pairs: None,
            revel_table: handed.revel_table,
            utrannotator_sweep: None,
            downstream_sequences: None,
            splice_windows: handed.splice_windows,
            spliceai_table: None,
            exon_spans: None,
        })
        .unwrap()
}

/// Whether `klass` takes the field `name` on the reference's first row of `allele`, handed
/// the engine's row under the same key and each output's rows of the row's record where
/// the output states them.
fn takes_with(
    klass: &str,
    reference: &str,
    candidate: &str,
    allele: &str,
    name: &str,
    handed: Handed,
) -> bool {
    let reference = parse_default(reference).unwrap();
    let candidate = parse_default(candidate).unwrap();
    let no_lines: [&str; 0] = [];
    let e = reference.iter().find(|row| row.allele == allele).unwrap();
    let a = candidate.iter().find(|row| row.key == e.key);
    let record_rows = RecordRows::new(&reference, &no_lines).of(e).unwrap();
    let candidate_rows = a.and_then(|a| RecordRows::new(&candidate, &no_lines).of(a).unwrap());
    let value = |row: Option<&Row>| -> String {
        row.and_then(|row| row.fields.get(name))
            .unwrap_or(ABSENT)
            .to_string()
    };
    consult(
        klass,
        e,
        name,
        (&value(Some(e)), &value(a)),
        (a, record_rows.as_deref(), candidate_rows.as_deref()),
        handed,
    )
}

fn takes(klass: &str, reference: &str, candidate: &str, allele: &str, name: &str) -> bool {
    takes_with(klass, reference, candidate, allele, name, Handed::default())
}

const AM_CLASS: &str = "alphamissense_sibling_allele_rows";

/// Rows of the record 21 13610388 A C,T under `name`: (allele, feature, Amino_acids,
/// am_class, am_pathogenicity) each, an empty class leaving the two fields out.
fn am_rows(name: &str, rows: &[[&str; 5]]) -> String {
    let rows: Vec<String> = rows
        .iter()
        .map(|[allele, feature, amino_acids, class, score]| {
            let extra = if class.is_empty() {
                "IMPACT=MODERATE;STRAND=1".to_string()
            } else {
                format!("IMPACT=MODERATE;STRAND=1;am_class={class};am_pathogenicity={score}")
            };
            output_row(
                [
                    name,
                    "21:13610388",
                    allele,
                    feature,
                    "missense_variant",
                    &extra,
                ],
                "-",
                amino_acids,
            )
        })
        .collect();
    default_text(&rows)
}

#[test]
fn alphamissense_sibling_allele_rows_takes_the_sibling_alleles_prediction_and_no_other() {
    let two = "21_13610388_A/C/T";
    let c_row = ["C", "ENST00000299443", "M/L", "likely_benign", "0.1079"];
    let t_reference = ["T", "ENST00000299443", "M/L", "likely_benign", "0.1079"];
    let t_candidate = ["T", "ENST00000299443", "M/L", "likely_pathogenic", "0.9000"];
    let reference = am_rows(two, &[c_row, t_reference]);
    let candidate = am_rows(two, &[c_row, t_candidate]);
    for field in ["am_class", "am_pathogenicity"] {
        assert!(
            takes(AM_CLASS, &reference, &candidate, "T", field),
            "{field}"
        );
    }
    let without_row = am_rows(two, &[c_row, ["T", "ENST00000299443", "M/L", "", ""]]);
    assert!(takes(AM_CLASS, &reference, &without_row, "T", "am_class"));
    // Under an identifier a row does not say whether its Location's rows are one record's.
    let named = |t| am_rows("sibling_ct", &[c_row, t]);
    assert!(!takes(
        AM_CLASS,
        &named(t_reference),
        &named(t_candidate),
        "T",
        "am_class"
    ));
    let one = "21_13610388_A/T";
    assert!(!takes(
        AM_CLASS,
        &am_rows(one, &[t_reference]),
        &am_rows(one, &[t_candidate]),
        "T",
        "am_class"
    ));
    for sibling in [
        ["C", "ENST00000299443", "M/V", "likely_benign", "0.1079"],
        ["C", "ENST00000620442", "M/L", "likely_benign", "0.1079"],
        ["C", "ENST00000299443", "M/L", "ambiguous", "0.4500"],
    ] {
        let reference = am_rows(two, &[sibling, t_reference]);
        let candidate = am_rows(two, &[sibling, t_candidate]);
        assert!(
            !takes(AM_CLASS, &reference, &candidate, "T", "am_class"),
            "{sibling:?}"
        );
    }
    let t_empty = ["T", "ENST00000299443", "M/L", "", ""];
    assert!(!takes(
        AM_CLASS,
        &am_rows(two, &[c_row, t_empty]),
        &candidate,
        "T",
        "am_class"
    ));
    // The row's fields are read together: a score of its own is no copy of the sibling's row.
    let own_score = ["T", "ENST00000299443", "M/L", "likely_benign", "0.5000"];
    let other_class = ["T", "ENST00000299443", "M/L", "likely_pathogenic", "0.5000"];
    assert!(!takes(
        AM_CLASS,
        &am_rows(two, &[c_row, own_score]),
        &am_rows(two, &[c_row, other_class]),
        "T",
        "am_class"
    ));
}

const REVEL_SIBLING: &str = "revel_sibling_allele_rows";
const REVEL_OTHER: &str = "revel_other_reference_amino_acid_rows";
const TRANSCRIPT: &str = "ENST00000299443";

/// Rows of the record 21 13610388 A C,T on one transcript under `name`: (allele,
/// Amino_acids, REVEL) each, an empty score leaving the field out.
fn revel_rows(name: &str, rows: &[[&str; 3]]) -> String {
    let rows: Vec<String> = rows
        .iter()
        .map(|[allele, amino_acids, score]| {
            let extra = if score.is_empty() {
                "IMPACT=MODERATE;STRAND=1".to_string()
            } else {
                format!("IMPACT=MODERATE;STRAND=1;REVEL={score}")
            };
            output_row(
                [
                    name,
                    "21:13610388",
                    allele,
                    TRANSCRIPT,
                    "missense_variant",
                    &extra,
                ],
                "-",
                amino_acids,
            )
        })
        .collect();
    default_text(&rows)
}

#[test]
fn revel_sibling_allele_rows_takes_the_sibling_alleles_score_and_no_other() {
    let two = "21_13610388_A/C/T";
    let leucine = ["C", "M/L", "0.025"];
    let reference = revel_rows(two, &[leucine, ["T", "M/L", "0.025"]]);
    let candidate = revel_rows(two, &[leucine, ["T", "M/L", "0.027"]]);
    assert!(takes(REVEL_SIBLING, &reference, &candidate, "T", "REVEL"));
    let empty = revel_rows(two, &[leucine, ["T", "M/L", ""]]);
    assert!(takes(REVEL_SIBLING, &reference, &empty, "T", "REVEL"));
    let named = |score| revel_rows("two_leucines", &[leucine, ["T", "M/L", score]]);
    assert!(!takes(
        REVEL_SIBLING,
        &named("0.025"),
        &named("0.027"),
        "T",
        "REVEL"
    ));
    let one = |score| revel_rows("21_13610388_A/T", &[["T", "M/L", score]]);
    assert!(!takes(
        REVEL_SIBLING,
        &one("0.025"),
        &one("0.027"),
        "T",
        "REVEL"
    ));
    let valine = |score| {
        revel_rows(
            "21_13610388_A/G/T",
            &[["G", "M/V", "0.025"], ["T", "M/L", score]],
        )
    };
    assert!(!takes(
        REVEL_SIBLING,
        &valine("0.025"),
        &valine("0.027"),
        "T",
        "REVEL"
    ));
    let other_score = |score| revel_rows(two, &[["C", "M/L", "0.043"], ["T", "M/L", score]]);
    assert!(!takes(
        REVEL_SIBLING,
        &other_score("0.025"),
        &other_score("0.027"),
        "T",
        "REVEL"
    ));
    assert!(!takes(REVEL_SIBLING, &empty, &candidate, "T", "REVEL"));
}

#[test]
fn the_sibling_allele_classes_refuse_the_named_alleles_own_row() {
    for (klass, value) in [
        (AM_CLASS, "am_class=likely_benign;am_pathogenicity=0.1079"),
        (REVEL_SIBLING, "REVEL=0.025"),
        (
            "splicevault_sibling_allele_rows",
            "SpliceVault_site_type=donor;SpliceVault_site_pos=5",
        ),
    ] {
        let rows = |filled: &str| {
            let rows: Vec<String> = ["C", "T"]
                .iter()
                .map(|allele| {
                    let extra = if filled.contains(allele) {
                        format!("IMPACT=MODERATE;STRAND=1;{value}")
                    } else {
                        "IMPACT=MODERATE;STRAND=1".to_string()
                    };
                    output_row(
                        [
                            "21_5000100_A/C/T",
                            "21:5000100",
                            allele,
                            "ENST00000000041",
                            "missense_variant",
                            &extra,
                        ],
                        "-",
                        "M/L",
                    )
                })
                .collect();
            default_text(&rows)
        };
        let (name, filled) = value.split(';').next().unwrap().split_once('=').unwrap();
        let reference = rows("CT");
        for (kept, set_aside) in [("C", "T"), ("T", "C")] {
            assert!(
                takes(klass, &reference, &rows(kept), set_aside, name),
                "{klass} {kept}"
            );
            assert!(
                !takes(klass, &reference, &rows(kept), kept, name),
                "{klass} {kept}"
            );
        }
        assert!(!takes(klass, &reference, &rows(""), "C", name), "{klass}");
        assert!(!takes(klass, &reference, &rows(""), "T", name), "{klass}");
        let reference_rows = parse_default(&reference).unwrap();
        let candidate_rows = parse_default(&rows("C")).unwrap();
        let (own, sibling) = (&reference_rows[0], &reference_rows[1]);
        let handed = Handed::default();
        assert!(consult(
            klass,
            sibling,
            name,
            (filled, ABSENT),
            (
                Some(&candidate_rows[1]),
                Some(&reference_rows),
                Some(&candidate_rows)
            ),
            handed
        ));
        assert!(!consult(
            klass,
            sibling,
            name,
            (filled, ABSENT),
            (Some(&candidate_rows[1]), Some(&reference_rows), None),
            handed
        ));
        let emptied = parse_default(&rows("")).unwrap();
        assert!(!consult(
            klass,
            own,
            name,
            (filled, ABSENT),
            (Some(&emptied[0]), Some(&reference_rows), Some(&emptied)),
            handed
        ));
    }
}

#[test]
fn splicevault_sibling_allele_rows_reads_the_rows_of_the_record_it_is_handed() {
    let klass = "splicevault_sibling_allele_rows";
    let first = "SpliceVault_site_type=donor;SpliceVault_site_pos=5";
    let rows = |g_fields: &str, name: &str, alleles: &str| {
        let rows: Vec<String> = [("C", first), ("G", g_fields)]
            .iter()
            .filter(|(allele, _)| alleles.contains(allele))
            .map(|(allele, fields)| {
                let extra = if fields.is_empty() {
                    "IMPACT=HIGH;STRAND=1".to_string()
                } else {
                    format!("IMPACT=HIGH;STRAND=1;{fields}")
                };
                default_row(
                    name,
                    "21:100",
                    allele,
                    "ENST1",
                    "splice_donor_variant",
                    &extra,
                )
            })
            .collect();
        default_text(&rows)
    };
    let two = "21_100_T/C/G";
    let own = "SpliceVault_site_type=acceptor;SpliceVault_site_pos=8";
    let reference = rows(first, two, "CG");
    for candidate in [rows("", two, "CG"), rows(own, two, "CG")] {
        assert!(takes(
            klass,
            &reference,
            &candidate,
            "G",
            "SpliceVault_site_type"
        ));
    }
    assert!(!takes(
        klass,
        &rows("", two, "CG"),
        &reference,
        "G",
        "SpliceVault_site_type"
    ));
    for (expected, candidate) in [
        (rows(first, "21_100_T/G", "G"), rows("", "21_100_T/G", "G")),
        (rows(first, "21_100_T/G", "G"), rows(own, "21_100_T/G", "G")),
        (rows(first, "sibling", "CG"), rows("", "sibling", "CG")),
        (
            rows(
                "SpliceVault_site_type=acceptor;SpliceVault_site_pos=6",
                two,
                "CG",
            ),
            rows(own, two, "CG"),
        ),
    ] {
        assert!(!takes(
            klass,
            &expected,
            &candidate,
            "G",
            "SpliceVault_site_type"
        ));
    }
    let reference_rows = parse_default(&reference).unwrap();
    let entry = &reference_rows[1];
    let kept = parse_default(&rows("", two, "CG")).unwrap();
    let handed = Handed::default();
    let name = "SpliceVault_site_type";
    assert!(!consult(
        klass,
        entry,
        name,
        ("donor", "-"),
        (None, None, None),
        handed
    ));
    assert!(consult(
        klass,
        entry,
        name,
        ("donor", "-"),
        (None, Some(&reference_rows), Some(&kept)),
        handed
    ));
    assert!(!consult(
        klass,
        entry,
        name,
        ("donor", "-"),
        (None, Some(&reference_rows), None),
        handed
    ));
    // The rows of another record at the Location are no sibling's.
    let other_record = parse_default(&rows(first, "21_100_A/C/G", "CG")).unwrap();
    assert!(!consult(
        klass,
        entry,
        name,
        ("donor", "-"),
        (None, Some(&other_record), Some(&kept)),
        handed
    ));
}

#[test]
fn trimmed_allele_keys_follow_the_trimming_orders_of_the_allele_match() {
    let keys = |reference, alternate, position| -> Vec<(String, String, i128)> {
        trimmed_allele_keys(reference, alternate, position)
            .into_iter()
            .collect()
    };
    let key = |reference: &str, alternate: &str, position: i128| {
        (reference.to_string(), alternate.to_string(), position)
    };
    assert_eq!(keys("A", "G", 100), [key("A", "G", 100)]);
    assert_eq!(keys("CA", "CG", 100), [key("A", "G", 101)]);
    assert_eq!(keys("AC", "GC", 100), [key("A", "G", 100)]);
    assert_eq!(keys("AC", "A", 100), [key("C", "-", 101)]);
    assert_eq!(keys("C", "-", 101), [key("C", "-", 101)]);
    assert_eq!(
        keys("CAC", "C", 100),
        [key("AC", "-", 101), key("CA", "-", 100)]
    );
}

fn revel_cells(table: &RevelTable, rows: &[TableRow]) -> Vec<RevelRow<'static>> {
    // The cells are leaked so that a test may hold them past the rows.
    rows.iter()
        .map(|row| {
            let row: &'static TableRow = Box::leak(Box::new(row.clone()));
            table.row(row)
        })
        .collect()
}

fn revel_row<'a>(
    position: &'a str,
    reference: Option<&'a str>,
    alternate: Option<&'a str>,
    amino_acids: (&'a str, &'a str),
    score: Option<&'a str>,
    transcript_ids: Option<&'a str>,
) -> RevelRow<'a> {
    RevelRow {
        position,
        reference,
        alternate,
        reference_amino_acid: amino_acids.0,
        alternate_amino_acid: amino_acids.1,
        score,
        transcript_ids,
    }
}

#[test]
fn the_revel_table_reads_the_rows_of_a_span_in_file_order_through_the_index() {
    for (corpus, column) in [("GRCh37-plugin-revel", 1), ("GRCh38-plugin-revel", 2)] {
        let path = golden(corpus);
        let table = RevelTable::new(&path, false).unwrap();
        let mut text = String::new();
        MultiGzDecoder::new(std::fs::File::open(&path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        let mut held: BTreeMap<i128, Vec<Vec<String>>> = BTreeMap::new();
        for line in text.lines().filter(|line| !line.starts_with('#')) {
            let cells: Vec<String> = line.split('\t').map(str::to_string).collect();
            held.entry(cells[column].parse().unwrap())
                .or_default()
                .push(cells);
        }
        assert!(held.len() > 100, "{corpus}");
        let agree = |start: i128, end: i128| {
            let rows = table.rows("21", start, end).unwrap();
            let expected: Vec<&Vec<String>> =
                held.range(start..=end).flat_map(|(_, rows)| rows).collect();
            assert_eq!(rows.len(), expected.len(), "{corpus} {start} {end}");
            for (row, cells) in rows.iter().zip(expected) {
                let cell = |at: usize| Some(cells[at].as_str());
                assert_eq!(
                    table.row(row),
                    revel_row(
                        &cells[column],
                        cell(3),
                        cell(4),
                        (&cells[5], &cells[6]),
                        cell(7),
                        cell(8)
                    ),
                    "{corpus} {start} {end}"
                );
            }
        };
        let positions: Vec<i128> = held.keys().copied().collect();
        for at in &positions {
            agree(*at, *at);
        }
        let last = positions[positions.len() - 1];
        agree(positions[0], last);
        agree(positions[3] - 40, positions[40] + 7);
        agree(last + 1, last + 9);
        let at = positions[5];
        assert_eq!(
            table.rows("chr21", at, at).unwrap(),
            table.rows("21", at, at).unwrap()
        );
        assert!(table.rows("22", at, at).unwrap().is_empty());
        assert!(table.rows("21", at, at - 1).unwrap().is_empty());
    }
}

#[test]
fn the_revel_table_reads_each_layout_and_refuses_what_is_no_table() {
    let scratch = Scratch::new("revel-layouts");
    let nine = write_table(
        &scratch,
        "revel.tsv.gz",
        TABLE_HEADER,
        &[
            "21\t100\t500\tA\tG\tQ\tH\t0.742\tENST1;ENST2",
            "21\t100\t500\tA\tT\tQ\tL\t0.311",
            "21\t100\t500\tA",
            "21\t101\t501\tA\tG\tQ\tH\t\tENST1",
        ],
        3,
    );
    let table = RevelTable::new(&nine, false).unwrap();
    let rows = table.rows("21", 500, 501).unwrap();
    assert_eq!(
        revel_cells(&table, &rows),
        [
            revel_row(
                "500",
                Some("A"),
                Some("G"),
                ("Q", "H"),
                Some("0.742"),
                Some("ENST1;ENST2")
            ),
            revel_row("500", Some("A"), Some("T"), ("Q", "L"), Some("0.311"), None),
            revel_row("500", Some("A"), None, ("", ""), None, None),
            revel_row(
                "501",
                Some("A"),
                Some("G"),
                ("Q", "H"),
                Some(""),
                Some("ENST1")
            ),
        ]
    );
    let seven = write_table(
        &scratch,
        "seven.tsv.gz",
        "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL",
        &["21\t100\tA\tG\tQ\tH\t0.742"],
        2,
    );
    let table = RevelTable::new(&seven, false).unwrap();
    assert_eq!(
        revel_cells(&table, &table.rows("21", 100, 100).unwrap()),
        [revel_row(
            "100",
            Some("A"),
            Some("G"),
            ("Q", "H"),
            Some("0.742"),
            None
        )]
    );
    let grch37 = write_table(
        &scratch,
        "grch37.tsv.gz",
        TABLE_HEADER,
        &["21\t100\t500\tA\tG\tQ\tH\t0.742\tENST1"],
        2,
    );
    let table = RevelTable::new(&grch37, false).unwrap();
    assert_eq!(table.rows("21", 100, 100).unwrap()[0].position, "100");
    assert!(table.rows("21", 500, 500).unwrap().is_empty());
    let six = write_table(
        &scratch,
        "six.tsv.gz",
        "#chr\thg19_pos\tgrch38_pos\tref\talt\tREVEL",
        &["21\t100\t500\tA\tG\t0.742"],
        3,
    );
    assert_eq!(
        refusal(RevelTable::new(&six, false).map(|_| ())),
        format!(
            "--revel-table: {} is not a REVEL table (a header of 7, 8 or 9 columns)",
            six.display()
        )
    );
    std::fs::remove_file(scratch.0.join("revel.tsv.gz.tbi")).unwrap();
    assert_eq!(
        refusal(RevelTable::new(&nine, true).map(|_| ())),
        format!(
            "--revel-table: cannot read the tabix index of {table}: none of {table}.csi, {stem}.csi, {table}.tbi, {stem}.tbi exists",
            table = nine.display(),
            stem = nine.with_extension("").display()
        )
    );
}

/// A made-up nine-column table in its own directory whose rows stand at 21 `position` with
/// the reference allele `A`: (alternate allele, aaref, aaalt, score, transcript ids) each,
/// ids of `None` leaving the ninth cell out.
fn revel_table(
    scratch: &Scratch,
    directory: &str,
    rows: &[(&str, &str, &str, &str, Option<&str>)],
    no_match: bool,
    position: i128,
) -> RevelTable {
    let lines: Vec<String> = rows
        .iter()
        .map(|(alternate, aaref, aaalt, score, ids)| {
            let ids = ids.map_or(String::new(), |ids| format!("\t{ids}"));
            format!("21\t999\t{position}\tA\t{alternate}\t{aaref}\t{aaalt}\t{score}{ids}")
        })
        .collect();
    let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
    let path = write_table(
        scratch,
        &format!("{directory}/revel.tsv.gz"),
        TABLE_HEADER,
        &lines,
        3,
    );
    RevelTable::new(&path, no_match).unwrap()
}

#[test]
fn revel_other_reference_amino_acid_rows_takes_the_score_of_a_row_of_another_reference_amino_acid()
{
    let scratch = Scratch::new("revel-other");
    let at = 13610388;
    let name = "21_13610388_A/T";
    let reference = revel_rows(name, &[["T", "M/L", "0.300"]]);
    let empty = revel_rows(name, &[["T", "M/L", ""]]);
    let under = |table: &RevelTable, reference: &str, candidate: &str| {
        let handed = Handed {
            revel_table: Some(table),
            ..Handed::default()
        };
        takes_with(REVEL_OTHER, reference, candidate, "T", "REVEL", handed)
    };
    let ids = Some(TRANSCRIPT);
    let other = revel_table(
        &scratch,
        "other",
        &[("T", "V", "L", "0.300", ids)],
        false,
        at,
    );
    assert!(under(&other, &reference, &empty));
    assert!(!takes(REVEL_OTHER, &reference, &empty, "T", "REVEL"));
    assert!(!under(
        &other,
        &reference,
        &revel_rows(name, &[["T", "M/L", "0.999"]])
    ));
    let own_pair = revel_table(&scratch, "own", &[("T", "M", "L", "0.300", ids)], false, at);
    assert!(!under(&own_pair, &reference, &empty));
    assert!(!under(
        &other,
        &revel_rows(name, &[["T", "M/L", "0.555"]]),
        &empty
    ));
    let elsewhere = revel_table(
        &scratch,
        "elsewhere",
        &[("T", "V", "L", "0.300", ids)],
        false,
        13610399,
    );
    assert!(!under(&elsewhere, &reference, &empty));
    let other_allele = revel_table(
        &scratch,
        "allele",
        &[("G", "V", "L", "0.300", ids)],
        false,
        at,
    );
    assert!(!under(&other_allele, &reference, &empty));
    let other_residue = revel_table(
        &scratch,
        "residue",
        &[("T", "V", "F", "0.300", ids)],
        false,
        at,
    );
    assert!(!under(&other_residue, &reference, &empty));
    // The engine's own walk goes on to the first row of the row's own pair.
    let later = revel_table(
        &scratch,
        "later",
        &[("T", "V", "L", "0.300", ids), ("T", "M", "L", "0.400", ids)],
        false,
        at,
    );
    let own_score = revel_rows(name, &[["T", "M/L", "0.400"]]);
    assert!(under(&later, &reference, &own_score));
    assert!(!under(&later, &reference, &empty));
    assert!(!under(&other, &reference, &own_score));
    // The two outputs spell the row's amino acids alike, as one pair.
    assert!(!under(
        &other,
        &reference,
        &revel_rows(name, &[["T", "V/L", ""]])
    ));
    // The cells are compared as written.
    let sign = revel_table(
        &scratch,
        "sign",
        &[("T", "*", "L", "0.300", ids)],
        false,
        at,
    );
    assert!(under(
        &sign,
        &revel_rows(name, &[["T", "U/L", "0.300"]]),
        &revel_rows(name, &[["T", "U/L", ""]])
    ));
    // Another plugin's field on the row is not the class's, and the JSON member is.
    let e = &parse_default(&reference).unwrap()[0];
    let a = &parse_default(&empty).unwrap()[0];
    let handed = Handed {
        revel_table: Some(&other),
        ..Handed::default()
    };
    assert!(consult(
        REVEL_OTHER,
        e,
        "REVEL",
        ("0.300", ABSENT),
        (Some(a), None, None),
        handed
    ));
    assert!(consult(
        REVEL_OTHER,
        e,
        "revel",
        ("0.300", ABSENT),
        (Some(a), None, None),
        handed
    ));
    assert!(!consult(
        REVEL_OTHER,
        e,
        "am_pathogenicity",
        ("0.300", ABSENT),
        (Some(a), None, None),
        handed
    ));
    assert!(!consult(
        REVEL_OTHER,
        e,
        "REVEL",
        ("0.300", ABSENT),
        (None, None, None),
        handed
    ));
    assert!(!consult(
        REVEL_OTHER,
        e,
        "REVEL",
        (ABSENT, "0.300"),
        (Some(a), None, None),
        handed
    ));
    assert!(!consult(
        REVEL_SIBLING,
        e,
        "REVEL",
        ("0.300", ABSENT),
        (Some(a), None, None),
        handed
    ));

    // REVEL.pm:191-197: a row whose ids cell is filled is the walk's only where it lists
    // the row's transcript, unless the run set `no_match`; a row with no ids is read on its
    // amino acid alone.
    let unlisted = [(
        "T",
        "V",
        "L",
        "0.300",
        Some("ENST00000620442;ENST00000620443"),
    )];
    assert!(!under(
        &revel_table(&scratch, "matched", &unlisted, false, at),
        &reference,
        &empty
    ));
    assert!(under(
        &revel_table(&scratch, "unmatched", &unlisted, true, at),
        &reference,
        &empty
    ));
    let listed = [(
        "T",
        "V",
        "L",
        "0.300",
        Some("ENST00000620442;ENST00000299443"),
    )];
    assert!(under(
        &revel_table(&scratch, "listed", &listed, false, at),
        &reference,
        &empty
    ));
    let no_ids = [("T", "V", "L", "0.300", None)];
    assert!(under(
        &revel_table(&scratch, "no-ids", &no_ids, false, at),
        &reference,
        &empty
    ));
    let two = [
        ("T", "V", "L", "0.300", Some("ENST00000620442")),
        ("T", "M", "L", "0.400", Some("ENST00000620443")),
    ];
    let unmatched = revel_table(&scratch, "two-unmatched", &two, true, at);
    assert!(under(&unmatched, &reference, &own_score));
    assert!(!under(&unmatched, &reference, &empty));
    assert!(!under(
        &revel_table(&scratch, "two-matched", &two, false, at),
        &reference,
        &own_score
    ));
    let seven = write_table(
        &scratch,
        "seven.tsv.gz",
        "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL",
        &["21\t13610388\tA\tT\tV\tL\t0.300"],
        2,
    );
    assert!(under(
        &RevelTable::new(&seven, false).unwrap(),
        &reference,
        &empty
    ));
}

const MITOCHONDRIAL: &str = "plugin_contig_mitochondrial_spelling";

/// The row of the record <contig> 3460 G A on one transcript under an identifier, its Extra
/// text ending with `fields` when there are any.
fn mitochondrial_row(fields: &str, location: &str) -> String {
    let extra = if fields.is_empty() {
        "IMPACT=MODERATE;STRAND=1".to_string()
    } else {
        format!("IMPACT=MODERATE;STRAND=1;{fields}")
    };
    default_text(&[default_row(
        "rs199476118",
        location,
        "A",
        "ENST00000361390",
        "missense_variant",
        &extra,
    )])
}

#[test]
fn plugin_contig_mitochondrial_spelling_takes_a_tabix_plugins_field_the_engine_alone_fills() {
    let on = |reference: &str, candidate: &str, location: &str, name: &str| {
        takes(
            MITOCHONDRIAL,
            &mitochondrial_row(reference, location),
            &mitochondrial_row(candidate, location),
            "A",
            name,
        )
    };
    assert!(on("", "ClinPred=0.5", "MT:3460", "ClinPred"));
    assert!(on(
        "",
        "SpliceAI_pred=GENE1|0.10|0.20|0.30|0.40|-4|-3|2|1",
        "MT:3460",
        "SpliceAI_pred"
    ));
    assert!(on("", "gnomAD_mean_cov=31.2", "MT:3460", "gnomAD_mean_cov"));
    let alphamissense = "am_class=likely_pathogenic;am_pathogenicity=0.624";
    assert!(on("", alphamissense, "MT:3460", "am_class"));
    assert!(on("", alphamissense, "MT:3460", "am_pathogenicity"));
    for contig in ["chrM", "chrMT"] {
        assert!(
            on(
                "",
                "gnomAD_mean_cov=31.2",
                &format!("{contig}:3460"),
                "gnomAD_mean_cov"
            ),
            "{contig}"
        );
    }
    assert!(!on("", "ClinPred=0.5", "21:3460", "ClinPred"));
    assert!(!on("", "ClinPred=0.5", "M:3460", "ClinPred"));
    assert!(!on(
        "",
        "gnomAD_hap_AF_hom=0.01",
        "MT:3460",
        "gnomAD_hap_AF_hom"
    ));
    assert!(!on("", "LoFtool=0.5", "MT:3460", "LoFtool"));
    assert!(!on("ClinPred=0.5", "", "MT:3460", "ClinPred"));
    assert!(!on("ClinPred=0.5", "ClinPred=0.6", "MT:3460", "ClinPred"));
    // A text row names its chromosome itself; a row without a Location takes it from the
    // record entry among the rows of its record.
    let text_row = &parse_default(&mitochondrial_row("", "MT:3460")).unwrap()[0];
    let handed = Handed::default();
    assert!(consult(
        MITOCHONDRIAL,
        text_row,
        "ClinPred",
        (ABSENT, "0.5"),
        (None, None, None),
        handed
    ));
    assert!(consult(
        MITOCHONDRIAL,
        text_row,
        "am_class",
        (ABSENT, "likely_pathogenic"),
        (None, None, None),
        handed
    ));
    let placed = |key: &str, record: Option<i64>| -> Row {
        Rc::new(Entry {
            key: key.to_string(),
            record,
            allele: "A".to_string(),
            feature: "ENST00000361390".to_string(),
            consequence: "missense_variant".to_string(),
            ..Entry::default()
        })
    };
    let mut record_line = Entry::clone(&placed("0|record", Some(0)));
    record_line.fields.insert("CHROM", "MT");
    let record_line: Row = Rc::new(record_line);
    let csq_entry = placed("0|A|ENST00000361390|Transcript", Some(0));
    assert!(!consult(
        MITOCHONDRIAL,
        &csq_entry,
        "ClinPred",
        ("", "0.5"),
        (None, None, None),
        handed
    ));
    let rows = [Rc::clone(&record_line), Rc::clone(&csq_entry)];
    assert!(consult(
        MITOCHONDRIAL,
        &csq_entry,
        "ClinPred",
        ("", "0.5"),
        (None, Some(&rows), None),
        handed
    ));
}

const DBNSFP_NAMED_M: &str = "dbnsfp_mitochondrial_record_named_m";

/// The record MT 3460 G A against a made-up dbNSFP table naming the contig M: dbNSFP.pm:313
/// names the record M to the table and :335-338 match the row at its position, so the
/// reference writes the row's value where the engine, which reads the table for the record
/// only where its rows state the assembly's bases, writes nothing. The class takes such a
/// field on a record spelt MT or chrMT where the run names GRCh37 and gives the field to
/// dbNSFP, or the field carries the plugin's own prefix.
#[test]
fn dbnsfp_mitochondrial_record_named_m_takes_a_dbnsfp_field_the_reference_alone_fills_on_grch37() {
    let asked = |fields: (&str, &str), location: &str, name: &str, handed: Handed| {
        takes_with(
            DBNSFP_NAMED_M,
            &mitochondrial_row(fields.0, location),
            &mitochondrial_row(fields.1, location),
            "A",
            name,
            handed,
        )
    };
    let dbnsfp = |assembly: Option<&'static str>| Handed {
        assembly,
        source: Some("dbNSFP"),
        ..Handed::default()
    };
    let scored = ("SIFT_score=0.01", "");
    for assembly in ["GRCh37", "grch37", "hg19", "37", " GRCh37\n"] {
        for location in ["MT:3460", "chrMT:3460"] {
            assert!(
                asked(scored, location, "SIFT_score", dbnsfp(Some(assembly))),
                "{assembly} {location}"
            );
        }
        // On a record spelt chrM or M the module hands the name on as written, and another
        // chromosome is no record of the class.
        for location in ["chrM:3460", "M:3460", "21:3460"] {
            assert!(
                !asked(scored, location, "SIFT_score", dbnsfp(Some(assembly))),
                "{assembly} {location}"
            );
        }
    }
    for assembly in [None, Some("GRCh38"), Some("TAIR10"), Some("grch37.p13")] {
        assert!(
            !asked(scored, "MT:3460", "SIFT_score", dbnsfp(assembly)),
            "{assembly:?}"
        );
    }
    let grch37 = dbnsfp(Some("GRCh37"));
    // A field both outputs fill and one the engine alone fills stay differences.
    assert!(!asked(
        ("SIFT_score=0.01", "SIFT_score=0.02"),
        "MT:3460",
        "SIFT_score",
        grch37
    ));
    assert!(!asked(
        ("", "SIFT_score=0.01"),
        "MT:3460",
        "SIFT_score",
        grch37
    ));
    // A column is the class's field where the run gives it to dbNSFP; the prefix is the
    // plugin's own (dbNSFP.pm:159, :257-261) and needs no plugin named.
    for source in [None, Some("ClinPred"), Some("TranscriptAnnotator")] {
        let handed = Handed {
            assembly: Some("GRCh37"),
            source,
            ..Handed::default()
        };
        assert!(
            !asked(scored, "MT:3460", "SIFT_score", handed),
            "{source:?}"
        );
        assert!(
            !asked(("ClinPred=0.5", ""), "MT:3460", "ClinPred", handed),
            "{source:?}"
        );
        assert!(
            asked(
                ("dbNSFP_SIFT_score=0.01", ""),
                "MT:3460",
                "dbNSFP_SIFT_score",
                handed
            ),
            "{source:?}"
        );
    }
    // A text row names its chromosome itself; a row without a Location takes it from the
    // record entry among the rows of its record.
    let placed = |key: &str, record: Option<i64>| -> Row {
        Rc::new(Entry {
            key: key.to_string(),
            record,
            allele: "A".to_string(),
            feature: "ENST00000361390".to_string(),
            consequence: "missense_variant".to_string(),
            ..Entry::default()
        })
    };
    let mut record_line = Entry::clone(&placed("0|record", Some(0)));
    record_line.fields.insert("CHROM", "chrMT");
    let record_line: Row = Rc::new(record_line);
    let csq_entry = placed("0|A|ENST00000361390|Transcript", Some(0));
    assert!(!consult(
        DBNSFP_NAMED_M,
        &csq_entry,
        "SIFT_score",
        ("0.01", ""),
        (None, None, None),
        grch37
    ));
    let rows = [Rc::clone(&record_line), Rc::clone(&csq_entry)];
    assert!(consult(
        DBNSFP_NAMED_M,
        &csq_entry,
        "SIFT_score",
        ("0.01", ""),
        (None, Some(&rows), None),
        grch37
    ));
}

/// Rows of one record on one transcript: (allele, Consequence, the LoF pairs of the Extra
/// text) each.
fn lof_rows(name: &str, location: &str, feature: &str, rows: &[(&str, &str, &str)]) -> String {
    let rows: Vec<String> = rows
        .iter()
        .map(|(allele, consequence, lof)| {
            let extra = if lof.is_empty() {
                "IMPACT=HIGH;STRAND=1".to_string()
            } else {
                format!("IMPACT=HIGH;STRAND=1;{lof}")
            };
            default_row(name, location, allele, feature, consequence, &extra)
        })
        .collect();
    default_text(&rows)
}

#[test]
fn loftee_anc_allele_multiallelic_record_takes_the_filter_the_engine_alone_applies() {
    let klass = "loftee_anc_allele_multiallelic_record";
    let sibling = ("A", "synonymous_variant", "");
    let pair = |reference_lof: &'static str,
                candidate_lof: &'static str,
                siblings: &[(&'static str, &'static str, &'static str)],
                name: &str,
                location: &str,
                consequence: &'static str| {
        let side = |lof: &'static str| {
            let mut rows = vec![("T", consequence, lof)];
            rows.extend_from_slice(siblings);
            lof_rows(name, location, "ENST00000339944", &rows)
        };
        (side(reference_lof), side(candidate_lof))
    };
    let verdict = |reference_lof,
                   candidate_lof,
                   siblings: &[_],
                   name: &str,
                   location: &str,
                   consequence,
                   field: &str| {
        let (reference, candidate) = pair(
            reference_lof,
            candidate_lof,
            siblings,
            name,
            location,
            consequence,
        );
        takes(klass, &reference, &candidate, "T", field)
    };
    let (two, at, stop) = ("21_32314600_C/T/A", "21:32314600", "stop_gained");
    let (filtered, unfiltered) = (
        "LoF=LC;LoF_filter=END_TRUNC,ANC_ALLELE",
        "LoF=LC;LoF_filter=END_TRUNC",
    );
    assert!(verdict(
        unfiltered,
        filtered,
        &[sibling],
        two,
        at,
        stop,
        "LoF_filter"
    ));
    // Where the item is the only filter the row's LoF reads HC against LC.
    for field in ["LoF", "LoF_filter"] {
        assert!(
            verdict(
                "LoF=HC",
                "LoF=LC;LoF_filter=ANC_ALLELE",
                &[sibling],
                two,
                at,
                stop,
                field
            ),
            "{field}"
        );
    }
    assert!(!verdict(
        unfiltered,
        filtered,
        &[],
        "21_32314600_C/T",
        at,
        stop,
        "LoF_filter"
    ));
    assert!(!verdict(
        unfiltered,
        filtered,
        &[sibling],
        "two",
        at,
        stop,
        "LoF_filter"
    ));
    // The record's other allele on another feature alone is no sibling row of this one.
    let elsewhere = |lof: &str| {
        default_text(&[
            default_row(
                two,
                at,
                "T",
                "ENST00000339944",
                stop,
                &format!("IMPACT=HIGH;STRAND=1;{lof}"),
            ),
            default_row(
                two,
                at,
                "A",
                "ENST2",
                "synonymous_variant",
                "IMPACT=LOW;STRAND=1",
            ),
        ])
    };
    assert!(!takes(
        klass,
        &elsewhere(unfiltered),
        &elsewhere(filtered),
        "T",
        "LoF_filter"
    ));
    assert!(!verdict(
        filtered,
        unfiltered,
        &[sibling],
        two,
        at,
        stop,
        "LoF_filter"
    ));
    assert!(!verdict(
        unfiltered,
        "LoF=LC;LoF_filter=ANC_ALLELE",
        &[sibling],
        two,
        at,
        stop,
        "LoF_filter"
    ));
    assert!(!verdict(
        unfiltered,
        "LoF=LC;LoF_filter=END_TRUNC,ANC_ALLELE,ANC_ALLELE",
        &[sibling],
        two,
        at,
        stop,
        "LoF_filter"
    ));
    assert!(!verdict(
        unfiltered,
        filtered,
        &[("TA", "frameshift_variant", "")],
        "21_32314601_-/T/TA",
        "21:32314600-32314601",
        stop,
        "LoF_filter"
    ));
    assert!(!verdict(
        unfiltered,
        filtered,
        &[sibling],
        two,
        at,
        "missense_variant",
        "LoF_filter"
    ));
    assert!(!verdict(
        "LoF=HC;LoF_filter=END_TRUNC",
        "LoF=LC;LoF_filter=END_TRUNC",
        &[sibling],
        two,
        at,
        stop,
        "LoF"
    ));
    for field in ["LoF", "LoF_filter"] {
        assert!(
            !verdict(
                "LoF=LC",
                "LoF=HC;LoF_filter=ANC_ALLELE",
                &[sibling],
                two,
                at,
                stop,
                field
            ),
            "{field}"
        );
        assert!(
            !verdict(
                "LoF=HC",
                "LoF=HC;LoF_filter=ANC_ALLELE",
                &[sibling],
                two,
                at,
                stop,
                field
            ),
            "{field}"
        );
    }
}

const DONOR_FILTERED: &str =
    "LoF=LC;LoF_filter=GC_TO_GT_DONOR;LoF_flags=NON_CAN_SPLICE;LoF_info=INTRON_SIZE:261";
const DONOR_SPARED: &str = "LoF=HC;LoF_flags=NON_CAN_SPLICE;LoF_info=INTRON_SIZE:261";

#[test]
fn loftee_first_alt_donor_filter_takes_the_filter_read_from_the_first_alternate() {
    let klass = "loftee_first_alt_donor_filter";
    let donor = "splice_donor_variant";
    let donor_rows = |name: &str, consequence: &str, rows: &[(&str, &str)]| {
        let rows: Vec<(&str, &str, &str)> = rows
            .iter()
            .map(|(allele, lof)| (*allele, consequence, *lof))
            .collect();
        lof_rows(name, "21:46126243", "ENST00001028357", &rows)
    };
    let first_t = "21_46126243_C/T/A";
    let t_row = ("T", DONOR_FILTERED);
    let reference = donor_rows(first_t, donor, &[t_row, ("A", DONOR_FILTERED)]);
    let candidate = donor_rows(first_t, donor, &[t_row, ("A", DONOR_SPARED)]);
    for field in ["LoF", "LoF_filter"] {
        assert!(takes(klass, &reference, &candidate, "A", field), "{field}");
    }
    let first_a = "21_46126243_C/A/T";
    let spared_first = donor_rows(first_a, donor, &[("A", DONOR_SPARED), ("T", DONOR_SPARED)]);
    let filtered_own = donor_rows(
        first_a,
        donor,
        &[("A", DONOR_SPARED), ("T", DONOR_FILTERED)],
    );
    assert!(takes(
        klass,
        &spared_first,
        &filtered_own,
        "T",
        "LoF_filter"
    ));

    // The item stays a difference on a one-allele record, on rows under an identifier, on
    // the first alternate's own row, on a row that is not a splice-donor row, beside
    // another differing item, and where the row's LoF does not read LC beside the item
    // and HC without it.
    let one = "21_46126243_C/A";
    let acceptor = "splice_acceptor_variant";
    let other_filter = DONOR_SPARED.replace("LoF=HC", "LoF=LC;LoF_filter=END_TRUNC");
    let other_call = DONOR_SPARED.replace("LoF=HC", "LoF=OS");
    let filtered_high = DONOR_FILTERED.replace("LoF=LC", "LoF=HC");
    let spared_low = DONOR_SPARED.replace("LoF=HC", "LoF=LC");
    let other_size = DONOR_FILTERED.replace("INTRON_SIZE:261", "INTRON_SIZE:262");
    for (reference, candidate, allele, fields) in [
        (
            donor_rows(one, donor, &[("A", DONOR_FILTERED)]),
            donor_rows(one, donor, &[("A", DONOR_SPARED)]),
            "A",
            &["LoF", "LoF_filter"][..],
        ),
        (
            donor_rows("497375", donor, &[t_row, ("A", DONOR_FILTERED)]),
            donor_rows("497375", donor, &[t_row, ("A", DONOR_SPARED)]),
            "A",
            &["LoF", "LoF_filter"],
        ),
        (
            reference.clone(),
            donor_rows(
                first_t,
                donor,
                &[("T", DONOR_SPARED), ("A", DONOR_FILTERED)],
            ),
            "T",
            &["LoF", "LoF_filter"],
        ),
        (
            donor_rows(first_t, acceptor, &[t_row, ("A", DONOR_FILTERED)]),
            donor_rows(first_t, acceptor, &[t_row, ("A", DONOR_SPARED)]),
            "A",
            &["LoF", "LoF_filter"],
        ),
        (
            reference.clone(),
            donor_rows(first_t, donor, &[t_row, ("A", &other_filter)]),
            "A",
            &["LoF_filter"],
        ),
        (
            reference.clone(),
            donor_rows(first_t, donor, &[t_row, ("A", &other_call)]),
            "A",
            &["LoF", "LoF_filter"],
        ),
        (
            donor_rows(first_t, donor, &[t_row, ("A", &filtered_high)]),
            donor_rows(first_t, donor, &[t_row, ("A", &spared_low)]),
            "A",
            &["LoF", "LoF_filter"],
        ),
        (
            reference.clone(),
            donor_rows(first_t, donor, &[t_row, ("A", &other_size)]),
            "A",
            &["LoF_info"],
        ),
    ] {
        for field in fields {
            assert!(
                !takes(klass, &reference, &candidate, allele, field),
                "{allele} {field}"
            );
        }
    }
    // Rows under an identifier state no allele order: where a caller vouches that they are
    // one record's, neither the first alternate's own row nor a later allele's is accepted.
    let named = parse_default(&donor_rows(
        "497375",
        donor,
        &[t_row, ("A", DONOR_FILTERED)],
    ))
    .unwrap();
    for (index, lofs) in [
        [("T", DONOR_SPARED), ("A", DONOR_FILTERED)],
        [t_row, ("A", DONOR_SPARED)],
    ]
    .iter()
    .enumerate()
    {
        let e = &named[index];
        let a = &parse_default(&donor_rows("497375", donor, lofs)).unwrap()[index];
        assert_eq!(first_alternate(e, &named), None);
        let rows = (Some(a), Some(&named[..]), None);
        let handed = Handed::default();
        assert!(
            !consult(
                klass,
                e,
                "LoF_filter",
                ("GC_TO_GT_DONOR", ABSENT),
                rows,
                handed
            ),
            "{index}"
        );
        assert!(
            !consult(klass, e, "LoF", ("LC", "HC"), rows, handed),
            "{index}"
        );
    }
}

/// The row of a frameshift on one transcript: its LoF pairs and a LoF_info carrying the
/// GERP total, the bases walked and the distance from the last exon.
fn gerp_row(lof: &str, gerp: &str, walked: &str, reach: &str, consequence: &str) -> String {
    let info = format!("PERCENTILE:0.95,GERP_DIST:{gerp},BP_DIST:{walked},DIST_FROM_LAST_EXON:{reach},50_BP_RULE:FAIL");
    default_text(&[default_row(
        "2780984",
        "21:36759739",
        "-",
        "ENST00001017196",
        consequence,
        &format!("IMPACT=HIGH;STRAND=1;{lof};LoF_info={info}"),
    )])
}

#[test]
fn loftee_gerp_zoom_level_summary_takes_a_total_read_through_a_zoom_level() {
    let klass = "loftee_gerp_zoom_level_summary";
    let (filtered, clear) = ("LoF=LC;LoF_filter=END_TRUNC", "LoF=HC");
    let (low, high) = ("-69.0908176183701", "-51.9470007419586");
    let frameshift = "frameshift_variant";
    let row = |lof, gerp| gerp_row(lof, gerp, "118", "12", frameshift);
    let fields = ["LoF", "LoF_filter", "LoF_info"];
    for field in fields {
        assert!(
            takes(klass, &row(filtered, low), &row(clear, high), "-", field),
            "{field}"
        );
        assert!(
            takes(klass, &row(clear, high), &row(filtered, low), "-", field),
            "{field}"
        );
    }
    assert!(takes(
        klass,
        &row(filtered, "-100.5"),
        &row(filtered, "-90.25"),
        "-",
        "LoF_info"
    ));
    let stop_gain = |lof, gerp| gerp_row(lof, gerp, "108", "-16", "stop_gained");
    assert!(takes(
        klass,
        &stop_gain(filtered, "-64.1545345306396"),
        &stop_gain(clear, "-34.7560015916824"),
        "-",
        "LoF_filter"
    ));
    // An interval of 80 bases, the shortest a zoom level serves, walks 79.
    let walked = |lof, gerp, walked| gerp_row(lof, gerp, walked, "12", frameshift);
    assert!(takes(
        klass,
        &walked(filtered, low, "79"),
        &walked(clear, high, "79"),
        "-",
        "LoF_info"
    ));
    for (reference, candidate) in [
        (walked(filtered, low, "78"), walked(clear, high, "78")),
        (row(filtered, low), walked(clear, high, "119")),
        (
            gerp_row(filtered, low, "118", "12", "missense_variant"),
            gerp_row(clear, high, "118", "12", "missense_variant"),
        ),
        (row(filtered, "-69.50"), row(clear, "-69.5")),
    ] {
        for field in fields {
            assert!(!takes(klass, &reference, &candidate, "-", field), "{field}");
        }
    }
    // The total itself is the class's on these rows; the filter beside it is not.
    for (reference, candidate) in [
        (row(filtered, "-69.5"), row(clear, "-60.5")),
        (row(clear, low), row(filtered, high)),
        (
            gerp_row(filtered, low, "118", "51", frameshift),
            gerp_row(clear, high, "118", "51", frameshift),
        ),
        (
            row(filtered, low),
            row("LoF=LC;LoF_filter=ANC_ALLELE", high),
        ),
    ] {
        assert!(takes(klass, &reference, &candidate, "-", "LoF_info"));
        assert!(!takes(klass, &reference, &candidate, "-", "LoF_filter"));
    }
}

/// The row of one variant on one transcript, `variant` its name, Location, allele and
/// Feature: its Consequence and the LoF pairs of its Extra text, none for an empty `lof`.
fn walk_row(variant: [&str; 4], consequence: &str, lof: &str) -> String {
    let [name, location, allele, feature] = variant;
    lof_rows(name, location, feature, &[(allele, consequence, lof)])
}

const WALK_STOP_GAIN: [&str; 4] = ["first_exon_stop", "21:37442543", "T", "ENST00000290349"];
const WALK_FRAMESHIFT: [&str; 4] = [
    "frameshift_before_utr_exon",
    "21:31588081",
    "-",
    "ENST00000286809",
];
const WALKED: &str = "LoF=HC;LoF_info=GERP_DIST:705,BP_DIST:702,PERCENTILE:0.155875299760192,DIST_FROM_LAST_EXON:266,50_BP_RULE:PASS";
const WALKED_FILTERED: &str = "LoF=LC;LoF_filter=END_TRUNC;LoF_info=GERP_DIST:150,BP_DIST:148,PERCENTILE:0.91,DIST_FROM_LAST_EXON:40,50_BP_RULE:FAIL";
const WALKED_FLAGGED: &str =
    "LoF=HC;LoF_flags=NO_EXON_NUMBER;LoF_info=GERP_DIST:612,BP_DIST:610,PERCENTILE:0.31";

#[test]
fn loftee_gerp_walk_without_conservation_database_takes_the_fields_the_engine_alone_writes() {
    let klass = "loftee_gerp_walk_without_conservation_database";
    let on = |variant: [&str; 4], consequence: &str, reference: &str, candidate: &str| {
        let reference = walk_row(variant, consequence, reference);
        let candidate = walk_row(variant, consequence, candidate);
        ["LoF", "LoF_info", "LoF_filter", "LoF_flags"]
            .into_iter()
            .filter(|field| takes(klass, &reference, &candidate, variant[2], field))
            .collect::<Vec<&str>>()
    };
    let (stop_gain, frameshift) = ("stop_gained", "frameshift_variant");
    let shifted_at_splice_region = "frameshift_variant,splice_region_variant";
    // The fields the engine writes on a row the reference leaves without the plugin's
    // fields are taken, and no field the engine leaves out. The last three rows are a run's
    // under `skip_lof_info`, which writes no LoF_info.
    for (variant, consequence, candidate, written) in [
        (WALK_STOP_GAIN, stop_gain, WALKED, &["LoF", "LoF_info"][..]),
        (WALK_FRAMESHIFT, frameshift, WALKED, &["LoF", "LoF_info"]),
        (
            WALK_STOP_GAIN,
            stop_gain,
            WALKED_FILTERED,
            &["LoF", "LoF_info", "LoF_filter"],
        ),
        (
            WALK_FRAMESHIFT,
            shifted_at_splice_region,
            WALKED_FLAGGED,
            &["LoF", "LoF_info", "LoF_flags"],
        ),
        (WALK_STOP_GAIN, stop_gain, "LoF=HC", &["LoF"]),
        (
            WALK_STOP_GAIN,
            stop_gain,
            "LoF=LC;LoF_filter=END_TRUNC",
            &["LoF", "LoF_filter"],
        ),
        (
            WALK_FRAMESHIFT,
            shifted_at_splice_region,
            "LoF=HC;LoF_flags=NO_EXON_NUMBER",
            &["LoF", "LoF_flags"],
        ),
    ] {
        assert_eq!(
            on(variant, consequence, "", candidate),
            written,
            "{consequence} {candidate}"
        );
    }
    // Nothing is taken where the engine writes no LoF, where its LoF_info holds items and
    // no GERP_DIST among them, or on a row of neither term.
    let none: [&str; 0] = [];
    for (consequence, candidate) in [
        (
            stop_gain,
            "LoF_info=GERP_DIST:705,BP_DIST:702,PERCENTILE:0.155875299760192",
        ),
        (stop_gain, "LoF_filter=END_TRUNC"),
        (
            stop_gain,
            "LoF=HC;LoF_info=BP_DIST:702,PERCENTILE:0.155875299760192",
        ),
        ("splice_donor_variant", WALKED),
        ("splice_acceptor_variant", "LoF=HC"),
        ("missense_variant", WALKED),
    ] {
        assert_eq!(
            on(WALK_STOP_GAIN, consequence, "", candidate),
            none,
            "{consequence} {candidate}"
        );
    }
    // Nor where the reference writes any of the plugin's four fields.
    for reference in [
        "LoF=HC",
        "LoF_info=PERCENTILE:0.91",
        "LoF_filter=END_TRUNC",
        "LoF_flags=SINGLE_EXON",
    ] {
        for candidate in [WALKED, WALKED_FILTERED, WALKED_FLAGGED, "LoF=LC"] {
            assert_eq!(
                on(WALK_FRAMESHIFT, frameshift, reference, candidate),
                none,
                "{reference} {candidate}"
            );
        }
    }
    // The engine's row holds one of the two terms as well; an entry that states no terms,
    // as a tuple file's does, is read by the reference's.
    let unwritten = walk_row(WALK_STOP_GAIN, stop_gain, "");
    for (terms, taken) in [
        ("stop_gained,splice_region_variant", true),
        (frameshift, true),
        ("splice_donor_variant", false),
        ("", true),
    ] {
        let candidate = walk_row(WALK_STOP_GAIN, terms, "LoF=HC");
        assert_eq!(
            takes(klass, &unwritten, &candidate, "T", "LoF"),
            taken,
            "{terms}"
        );
    }
    // The engine's row is read: a value handed without it is not taken.
    let reference = parse_default(&walk_row(WALK_STOP_GAIN, stop_gain, "")).unwrap();
    let candidate = parse_default(&walk_row(WALK_STOP_GAIN, stop_gain, WALKED)).unwrap();
    let handed = Handed::default();
    let values = (ABSENT, "HC");
    assert!(consult(
        klass,
        &reference[0],
        "LoF",
        values,
        (Some(&candidate[0]), None, None),
        handed
    ));
    assert!(!consult(
        klass,
        &reference[0],
        "LoF",
        values,
        (None, None, None),
        handed
    ));
}

/// One row whose Extra text ends with the NearestExonJB value, left out when `value` is
/// empty, and whose cDNA_position is `cdna`.
fn nearest_exon_row(row: [&str; 6], value: &str, cdna: &str) -> String {
    let [name, location, allele, feature, consequence, extra] = row;
    let extra = if value.is_empty() {
        extra.to_string()
    } else {
        format!("{extra};NearestExonJB={value}")
    };
    default_text(&[output_row(
        [name, location, allele, feature, consequence, &extra],
        cdna,
        "-",
    )])
}

const SHIFT_PROBE: [&str; 6] = [
    "shift_probe",
    "21:37480826",
    "-",
    "ENST00000338785",
    "splice_donor_variant",
    "IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1",
];

#[test]
fn nearest_exon_jb_stretch_priming_takes_a_value_within_twelve_bases_outside_the_transcript() {
    let klass = "nearest_exon_jb_stretch_priming";
    let row = |extra: &'static str, consequence: &'static str| {
        [
            "nine_before_first_exon",
            "21:6354255",
            "T",
            "ENST00000624847",
            consequence,
            extra,
        ]
    };
    let upstream = "upstream_gene_variant";
    let nine = row("IMPACT=MODIFIER;DISTANCE=9;STRAND=1", upstream);
    let value = "ENSE00003757672|9|start|259";
    let on = |row: [&str; 6], reference: &str, candidate: &str| {
        takes(
            klass,
            &nearest_exon_row(row, reference, "-"),
            &nearest_exon_row(row, candidate, "-"),
            "T",
            "NearestExonJB",
        )
    };
    assert!(on(nine, "", value));
    let twelve = row("IMPACT=MODIFIER;DISTANCE=12;STRAND=1", upstream);
    assert!(on(twelve, "", "ENSE00003757672|12|start|259"));
    let farther = row("IMPACT=MODIFIER;DISTANCE=13;STRAND=1", upstream);
    assert!(!on(farther, "", "ENSE00003757672|13|start|259"));
    assert!(!on(
        row("IMPACT=MODIFIER;STRAND=1", "intron_variant"),
        "",
        value
    ));
    assert!(!on(nine, value, ""));
    assert!(!on(nine, value, "ENSE00003757672|9|start|260"));
}

#[test]
fn nearest_exon_jb_shifted_span_takes_the_exon_set_read_at_the_shifted_span() {
    let klass = "nearest_exon_jb_shifted_span";
    let on = |row: [&str; 6], reference: &str, candidate: &str| {
        takes(
            klass,
            &nearest_exon_row(row, reference, "-"),
            &nearest_exon_row(row, candidate, "-"),
            "-",
            "NearestExonJB",
        )
    };
    let value = "ENSE00003584551|0|end|189";
    assert!(on(SHIFT_PROBE, "", value));
    let onto_exon = [
        "onto_exon",
        "21:37480637",
        "-",
        "ENST00000338785",
        "frameshift_variant,splice_region_variant",
        "IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=1",
    ];
    assert!(on(onto_exon, "ENSE00003584551|1|start|189", ""));
    assert!(on(
        SHIFT_PROBE,
        "ENSE00003584551|1|start|189",
        "ENSE00003584551|1|start|189,ENSE00003675684|0|end|100"
    ));
    let extra = |extra: &'static str| {
        let mut row = SHIFT_PROBE;
        row[5] = extra;
        row
    };
    assert!(!on(extra("IMPACT=HIGH;STRAND=1;SHIFT_LENGTH=0"), "", value));
    assert!(!on(extra("IMPACT=HIGH;STRAND=1"), "", value));
    let regulatory = |extra: &str| {
        default_text(&[format!(
            "shift_probe\t21:37480826\t-\tENSG1\tENSR00000141251\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\t{extra}"
        )])
    };
    assert!(!takes(
        klass,
        &regulatory("IMPACT=MODIFIER;SHIFT_LENGTH=1"),
        &regulatory(&format!(
            "IMPACT=MODIFIER;SHIFT_LENGTH=1;NearestExonJB={value}"
        )),
        "-",
        "NearestExonJB"
    ));
    assert!(!on(SHIFT_PROBE, value, "ENSE00003584551|1|end|189"));
    // Two spellings of an empty value are not a value against none.
    let entry = &parse_default(&nearest_exon_row(SHIFT_PROBE, "", "-")).unwrap()[0];
    let handed = Handed::default();
    assert!(!consult(
        klass,
        entry,
        "NearestExonJB",
        ("", "-"),
        (None, None, None),
        handed
    ));
    assert!(consult(
        klass,
        entry,
        "NearestExonJB",
        ("", value),
        (None, None, None),
        handed
    ));
}

#[test]
fn nearest_exon_items_are_read_in_every_formats_spelling() {
    let item = |exon_id, distance, side, length| NearestExonItem {
        exon_id,
        distance,
        side,
        length,
    };
    assert_eq!(
        nearest_exon_items("ENSE00003584551+1+start+189&ENSE00003675684+0+end+100"),
        Some(vec![
            item("ENSE00003584551", "1", "start", "189"),
            item("ENSE00003675684", "0", "end", "100")
        ])
    );
    assert_eq!(
        nearest_exon_items("|106|end|135,|2780|start|109"),
        Some(vec![
            item("", "106", "end", "135"),
            item("", "2780", "start", "109")
        ])
    );
    assert_eq!(
        nearest_exon_items("ENSE00003595154|87|start_end|175"),
        Some(vec![item("ENSE00003595154", "87", "start_end", "175")])
    );
    assert_eq!(nearest_exon_items("ENSE00003584551|1|start"), None);
    assert_eq!(nearest_exon_items("-"), None);
}

#[test]
fn nearest_exon_jb_id_less_exon_fold_takes_one_item_folded_from_two_id_less_exons() {
    let klass = "nearest_exon_jb_id_less_exon_fold";
    let two_exons = [
        "idless_two_exons",
        "21:10578401-10578530",
        "-",
        "compmerge.2.Heart.chr21",
        "splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,non_coding_transcript_exon_variant,intron_variant",
        "IMPACT=HIGH;STRAND=1;SOURCE=RefSeq",
    ];
    let intronic = [
        "idless_intronic",
        "21:10521800",
        "A",
        "compmerge.2.Heart.chr21",
        "intron_variant,non_coding_transcript_variant",
        "IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq",
    ];
    let on = |row: [&str; 6], reference: &str, candidate: &str, cdna: &str| {
        takes(
            klass,
            &nearest_exon_row(row, reference, cdna),
            &nearest_exon_row(row, candidate, cdna),
            row[2],
            "NearestExonJB",
        )
    };
    let closed = "1266-1297";
    assert!(on(two_exons, "|64|start|89", "|64|start|82", closed));
    for folded in ["|106|end|109", "|2780|start|109"] {
        assert!(
            on(intronic, folded, "|106|end|135,|2780|start|109", "-"),
            "{folded}"
        );
    }
    let mut one_exon = two_exons;
    one_exon[4] = "non_coding_transcript_exon_variant";
    assert!(!on(one_exon, "|64|start|89", "|64|start|82", closed));
    for cdna in ["-", "?-1297", "1266-?"] {
        assert!(
            !on(two_exons, "|64|start|89", "|64|start|82", cdna),
            "{cdna}"
        );
    }
    for (reference, candidate) in [
        ("ENSE00001659869|64|start|89", "ENSE00001659869|64|start|82"),
        ("|64|start|89", "ENSE00001659869|64|start|82"),
        ("|106|end|135,|2780|start|109", "|106|end|109"),
        ("|65|start|89", "|64|start|82"),
        ("|64|end|89", "|64|start|82"),
        ("|64|start|89", ""),
        ("", "|64|start|82"),
    ] {
        assert!(
            !on(two_exons, reference, candidate, closed),
            "{reference} {candidate}"
        );
    }
}

#[test]
fn nearest_exon_jb_item_order_takes_the_same_items_in_another_order() {
    let (first, second) = ("ENSE00000000040|7|start|120", "ENSE00000000041|2|end|80");
    assert!(reordered_items(
        &format!("{second},{first}"),
        &format!("{first},{second}")
    ));
    // The same value twice is no reordering.
    assert!(!reordered_items(
        &format!("{second},{first}"),
        &format!("{second},{first}")
    ));
    // Three items in another order, one of them written twice on both sides.
    let third = "ENSE00000000042|3|end|90";
    assert!(reordered_items(
        &format!("{third},{first},{second},{first}"),
        &format!("{first},{first},{second},{third}")
    ));
    assert!(!reordered_items(
        &format!("{third},{first},{second}"),
        &format!("{first},{first},{second},{third}")
    ));
    assert!(reordered_items("b,a", "a,b") && reordered_items("b&a", "a&b"));
    // A value carrying both list separators is split by neither.
    assert!(!reordered_items("a,b&c", "c&a,b") && !reordered_items("b&c,a", "a,b&c"));
    let row = |value: &str| {
        nearest_exon_row(
            [
                "a",
                "21:100",
                "T",
                "ENST1",
                "intron_variant",
                "IMPACT=MODIFIER;STRAND=1",
            ],
            value,
            "-",
        )
    };
    let klass = "nearest_exon_jb_item_order";
    assert!(takes(
        klass,
        &row(&format!("{second},{first}")),
        &row(&format!("{first},{second}")),
        "T",
        "NearestExonJB"
    ));
    assert!(!takes(
        klass,
        &row(&format!("{second},{first}")),
        &row(first),
        "T",
        "NearestExonJB"
    ));
    // The class reads the NearestExonJB field alone, under either spelling of its name.
    let items = ("ENSE00000000031|3|start|117", "ENSE00000000041|3|end|207");
    let terms = "splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,intron_variant";
    let extra = format!("IMPACT=HIGH;STRAND=1;NearestExonJB={},{}", items.1, items.0);
    let text = default_text(&[default_row(
        "exon_deletion",
        "21:1000011-1001793",
        "-",
        "ENST00000000031",
        terms,
        &extra,
    )]);
    let entry = &parse_default(&text).unwrap()[0];
    let reference_text = format!("{},{}", items.1, items.0);
    let candidate_text = format!("{},{}", items.0, items.1);
    let asked = |name: &str, values: (&str, &str)| {
        consult(
            klass,
            entry,
            name,
            values,
            (None, None, None),
            Handed::default(),
        )
    };
    assert!(asked("NearestExonJB", (&reference_text, &candidate_text)));
    assert!(asked("nearestexonjb", (&reference_text, &candidate_text)));
    for name in [
        "5UTR_annotation",
        "5UTR_consequence",
        "DOMAINS",
        "FREQS",
        "ZYG",
        "TSSDistance",
        "SpliceAI_pred",
        "NearestExonJB_x",
        "Extra",
    ] {
        assert!(!asked(name, (&reference_text, &candidate_text)), "{name}");
    }
    assert!(!asked("NearestExonJB", ("a,b&c", "c&a,b")));
}

const MAXENTSCAN: &str = "maxentscan_overlapping_windows";
const DONOR_TEXT: &str = "ACACACACG";
const ACCEPTOR_TEXT: &str = "ACACACACACACACACACACACG";

/// A cache of one shard on chromosome 1 with three transcripts of one exon layout: an exon
/// of three bases, 201 to 203, between an intron from 101 to 200 and one from 204 to 300.
/// On ENST1 (forward) the first intron's acceptor window is 181 to 203 and the second's
/// donor window 201 to 209; on ENST2 (reverse) the second intron comes first, its acceptor
/// window 201 to 223, and the first intron's donor window is 195 to 203. ENST3 has no such
/// exon: one intron from 101 to 300.
fn maxentscan_cache(scratch: &Scratch) -> PathBuf {
    let exons =
        r#"[{"start": 1, "end": 100}, {"start": 201, "end": 203}, {"start": 301, "end": 400}]"#;
    let shard = format!(
        r#"[{{"stable_id": "ENST1", "strand": 1, "exons": {exons}}}, {{"stable_id": "ENST2", "strand": -1, "exons": {exons}}}, {{"stable_id": "ENST3", "strand": 1, "exons": [{{"start": 1, "end": 100}}, {{"start": 301, "end": 400}}]}}]"#
    );
    scratch.file("cache/transcripts/1/1-1000000.json", shard.as_bytes());
    scratch.0.join("cache")
}

/// One substitution's row as the reference and as the engine write it: the three scores of
/// two different windows and, under `verbose`, the two window texts, the reference's a
/// donor's and the engine's an acceptor's.
fn maxentscan_texts(location: &str, feature: &str, verbose: bool) -> (String, String) {
    let side = |scores: [&str; 3], window: &str| {
        let mut extra = format!(
            "IMPACT=LOW;STRAND=1;MaxEntScan_ref={};MaxEntScan_alt={};MaxEntScan_diff={}",
            scores[0], scores[1], scores[2]
        );
        if verbose {
            let altered = format!("{}T", &window[..window.len() - 1]);
            extra += &format!(";MaxEntScan_ref_seq={window};MaxEntScan_alt_seq={altered}");
        }
        default_text(&[default_row(
            "a",
            location,
            "G",
            feature,
            "splice_region_variant",
            &extra,
        )])
    };
    (
        side(["8.100", "2.300", "5.800"], DONOR_TEXT),
        side(["6.400", "1.100", "5.300"], ACCEPTOR_TEXT),
    )
}

#[test]
fn maxentscan_overlapping_windows_takes_a_row_the_cache_places_in_two_windows_and_no_other() {
    let scratch = Scratch::new("maxentscan");
    let cache = maxentscan_cache(&scratch);
    let windows = SpliceWindows::new(&cache).unwrap();
    let handed = Handed {
        splice_windows: Some(&windows),
        ..Handed::default()
    };
    let scores = ["MaxEntScan_ref", "MaxEntScan_alt", "MaxEntScan_diff"];
    let texts = ["MaxEntScan_ref_seq", "MaxEntScan_alt_seq"];
    for feature in ["ENST1", "ENST2", "ENST2.4"] {
        for verbose in [false, true] {
            let (reference, candidate) = maxentscan_texts("1:202", feature, verbose);
            let fields = scores.iter().chain(texts.iter().filter(|_| verbose));
            for field in fields {
                assert!(
                    takes_with(MAXENTSCAN, &reference, &candidate, "G", field, handed),
                    "{feature} {field}"
                );
                assert!(
                    !takes(MAXENTSCAN, &reference, &candidate, "G", field),
                    "{feature} {field}"
                );
            }
        }
    }
    for (location, feature) in [
        ("1:199", "ENST1"),
        ("1:205", "ENST1"),
        ("1:150", "ENST1"),
        ("1:199", "ENST2"),
        ("1:205", "ENST2"),
        ("1:202", "ENST3"),
        ("1:202", "ENST9"),
        ("1:201-203", "ENST1"),
        ("2:202", "ENST1"),
    ] {
        let (reference, candidate) = maxentscan_texts(location, feature, true);
        for field in scores.iter().chain(&texts) {
            assert!(
                !takes_with(MAXENTSCAN, &reference, &candidate, "G", field, handed),
                "{location} {feature}"
            );
        }
    }
    let held: Vec<bool> = [
        ("1", 200),
        ("1", 201),
        ("1", 203),
        ("1", 204),
        ("chr1", 202),
    ]
    .iter()
    .map(|(name, position)| {
        let (reference, _) = maxentscan_texts(&format!("{name}:{position}"), "ENST1", false);
        let fresh = SpliceWindows::new(&cache).unwrap();
        fresh.holds(&parse_default(&reference).unwrap()[0]).unwrap()
    })
    .collect();
    assert_eq!(held, [false, true, true, false, true]);
    // A window text is taken where the reference writes a donor's nine bases and the engine
    // an acceptor's twenty-three, and in no other shape.
    let (reference, candidate) = maxentscan_texts("1:202", "ENST1", true);
    let reversed: String = DONOR_TEXT.chars().rev().collect();
    for other in [
        candidate.replace(ACCEPTOR_TEXT, &reversed),
        candidate.replace(ACCEPTOR_TEXT, &format!("{ACCEPTOR_TEXT}A")),
        candidate.replace(&format!(";MaxEntScan_ref_seq={ACCEPTOR_TEXT}"), ""),
        candidate.replace(ACCEPTOR_TEXT, &ACCEPTOR_TEXT.to_lowercase()),
    ] {
        assert_ne!(other, candidate);
        assert!(!takes_with(
            MAXENTSCAN,
            &reference,
            &other,
            "G",
            "MaxEntScan_ref_seq",
            handed
        ));
    }
    assert!(!takes_with(
        MAXENTSCAN,
        &candidate,
        &reference,
        "G",
        "MaxEntScan_ref_seq",
        handed
    ));
    // A score is taken where both engines write a number.
    let (reference, candidate) = maxentscan_texts("1:202", "ENST1", false);
    assert!(!takes_with(
        MAXENTSCAN,
        &reference,
        &candidate.replace("MaxEntScan_alt=1.100;", ""),
        "G",
        "MaxEntScan_alt",
        handed
    ));
    assert!(!takes_with(
        MAXENTSCAN,
        &reference.replace("8.100", "abc"),
        &candidate,
        "G",
        "MaxEntScan_ref",
        handed
    ));
    // An entry without a Location is placed by no input line of a run.
    let unplaced = Entry {
        key: "0|G|ENST1|Transcript".to_string(),
        record: Some(0),
        allele: "G".to_string(),
        feature: "ENST1".to_string(),
        consequence: "splice_region_variant".to_string(),
        ..Entry::default()
    };
    assert_eq!(windows.holds(&unplaced), Ok(false));
    assert_eq!(windows.place(&unplaced), None);
    assert_eq!(
        refusal(SpliceWindows::new(&scratch.0.join("none")).map(|_| ())),
        format!(
            "--json-cache: {} has no transcripts/ directory",
            scratch.0.join("none").display()
        )
    );
}

#[test]
fn splice_windows_read_integers_written_as_text_and_leave_out_what_two_shards_read_differently() {
    let scratch = Scratch::new("windows-shards");
    let exons = r#"[{"start": "1", "end": "100"}, {"start": "201", "end": "203"}, {"start": "301", "end": "400"}]"#;
    scratch.file(
        "transcripts/1/1-1000000.json",
        format!(
            r#"{{"1-500": [{{"stable_id": "TEXT", "strand": "1", "exons": {exons}}}, {{"stable_id": "TWICE", "strand": 1, "exons": {exons}}}, {{"stable_id": "SAME", "strand": "-1", "exons": {exons}}}, {{"stable_id": "FLOAT", "strand": 1.0, "exons": {exons}}}, {{"stable_id": "BOOL", "strand": true, "exons": {exons}}}, {{"stable_id": "OPEN", "strand": 1, "exons": [{{"start": 1}}, {{"start": 201, "end": 203}}, {{"start": 301, "end": 400}}]}}]}}"#
        )
        .as_bytes(),
    );
    scratch.file(
        "transcripts/1/1000001-2000000.json.gz",
        &gzip(
            format!(
                r#"[{{"stable_id": "TWICE", "strand": -1, "exons": {exons}}}, {{"stable_id": "SAME", "strand": -1, "exons": {exons}}}, null]"#
            )
            .as_bytes(),
        ),
    );
    let windows = SpliceWindows::new(&scratch.0).unwrap();
    let held = |feature: &str| {
        let (reference, _) = maxentscan_texts("chr1:202", feature, false);
        windows
            .holds(&parse_default(&reference).unwrap()[0])
            .unwrap()
    };
    assert_eq!(
        ["TEXT", "TWICE", "SAME", "FLOAT", "BOOL", "OPEN"].map(held),
        [true, false, true, false, false, false]
    );
}
