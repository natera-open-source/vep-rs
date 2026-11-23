// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The header lines and the blocks of an input, the tuples of an input cut into blocks
//! of any size on any number of threads, and the summary line.

use super::*;
use crate::columns::place_columns;
use crate::map::read_field_map;
use std::path::PathBuf;

/// A file of the system's temporary directory, removed when the test ends.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str, bytes: &[u8]) -> Scratch {
        // Tests run side by side, so a file's name is its own.
        static FILES: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "vep-fastvep-plugin-tuples-unit-{}-{}-{name}",
            std::process::id(),
            FILES.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, bytes).unwrap();
        Scratch(path)
    }

    fn input(&self, block: usize) -> Input {
        Input::with_blocks_of(File::open(&self.0).unwrap(), block)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

fn names(header: Option<Vec<Vec<u8>>>) -> Option<Vec<String>> {
    header.map(|cells| {
        cells
            .iter()
            .map(|cell| String::from_utf8_lossy(cell).into_owned())
            .collect()
    })
}

/// The header of an input and the rest of it, block by block, read in blocks of `block`.
fn header_and_blocks(bytes: &[u8], block: usize) -> (Option<Vec<String>>, Vec<Vec<u8>>) {
    let file = Scratch::new(&format!("blocks-{block}"), bytes);
    let mut input = file.input(block);
    let header = names(input.header().unwrap());
    let mut blocks = Vec::new();
    while let Some(block) = input.next_block().unwrap() {
        blocks.push(block);
    }
    (header, blocks)
}

/// An input, the cells of its column header and the bytes from its first row on.
type HeaderCase = (
    &'static [u8],
    Option<&'static [&'static str]>,
    &'static [u8],
);

#[test]
fn the_header_is_the_last_line_of_one_hash_before_the_first_row() {
    let cases: [HeaderCase; 12] = [
        (b"", None, b""),
        (b"row\n", None, b"row\n"),
        (b"\n#a\n", None, b"\n#a\n"),
        (b"## meta\n", None, b""),
        (b"#a\tb\n", Some(&["a", "b"]), b""),
        (b"#a\tb", Some(&["a", "b"]), b""),
        (b"#a\r\n", Some(&["a"]), b""),
        (b"#\n", Some(&[""]), b""),
        (
            b"#first\n## meta\n#second\tx\n## more\nrow\n#late\n",
            Some(&["second", "x"]),
            b"row\n#late\n",
        ),
        (b"## meta\n#a\nrow", Some(&["a"]), b"row"),
        (b"#a\n\nrow\n", Some(&["a"]), b"\nrow\n"),
        (b"\xef\xbb\xbf#a\nrow\n", None, b"\xef\xbb\xbf#a\nrow\n"),
    ];
    for (bytes, header, rest) in cases {
        let header: Option<Vec<String>> =
            header.map(|cells| cells.iter().map(|cell| cell.to_string()).collect());
        for block in [1, 2, 3, 5, 64, 1 << 20] {
            let (read, blocks) = header_and_blocks(bytes, block);
            assert_eq!(
                read,
                header,
                "{:?} in blocks of {block}",
                String::from_utf8_lossy(bytes)
            );
            assert_eq!(
                blocks.concat(),
                rest,
                "{:?} in blocks of {block}",
                String::from_utf8_lossy(bytes)
            );
        }
    }
}

#[test]
fn a_block_ends_at_a_line_end_and_the_blocks_are_the_input() {
    let mut bytes = b"#h\n".to_vec();
    for row in 0..200 {
        bytes.extend_from_slice(format!("row {row}\t{}\n", "x".repeat(row % 37)).as_bytes());
    }
    bytes.extend_from_slice(b"the last row has no line feed");
    for block in [1, 7, 16, 100, 1000, 1 << 20] {
        let (header, blocks) = header_and_blocks(&bytes, block);
        assert_eq!(header, Some(vec!["h".to_string()]));
        assert_eq!(blocks.concat(), &bytes[3..], "blocks of {block}");
        let (last, whole) = blocks.split_last().unwrap();
        assert!(
            whole.iter().all(|block| block.ends_with(b"\n")),
            "blocks of {block}"
        );
        assert!(last.ends_with(b"line feed"));
        if block <= 100 {
            assert!(blocks.len() > 20, "blocks of {block}");
        }
    }
    // A line longer than a block is whole in one block.
    let long = format!("#h\n{}\nshort\n", "y".repeat(5000));
    let (_, blocks) = header_and_blocks(long.as_bytes(), 64);
    assert!(blocks[0].len() >= 5001 && blocks[0].ends_with(b"\n"));
    assert_eq!(blocks.concat(), &long.as_bytes()[3..]);
}

const MAP: &str = "\
fastvep_column\tpart\tvep_field\tsource\tkind\tnote\tlist_separator\treference_set
FV_A\t1\tA_ONE\tS\ttext\tn\t-\t-
FV_A\t2\tA_LIST\tS\ttext\tn\t|\t-
FV_B\t1.score\tB_SCORE\tS\tnumber\tn\t-\t-
FV_B\t2.word\tB_WORD\tS\ttext\tn\t-\t-
";

/// An input of `rows` rows over the columns FV_A, FV_PHYLOP and FV_B, with comment lines,
/// short rows, repeated cells, cells with several entries and carriage returns among them.
fn table(rows: usize) -> Vec<u8> {
    let mut bytes =
        b"## fastVEP output\n#v\tLocation\tAllele\tGene\tFeature\tFV_A\tFV_PHYLOP\tFV_B\n".to_vec();
    for row in 0..rows {
        let line = match row % 11 {
            3 => "# a comment among the rows\n".to_string(),
            5 => format!("v\t21:{row}\tG\n"),
            7 => format!("v\t21:{row}\tT\tg\tENST{row}\tT|a%20b|x&y|drop\t1.5\t-\r\n"),
            9 => format!("v\t21:{row}\t-\tg\t-\tC|one,C|two|z\t-\tC|plain|w(0.5)\n"),
            _ => format!(
                "v\t21:{}\tG\tg\tENST{row}\tG|v{}|\t0.1\tG|d({})|\n",
                row / 4,
                row / 4,
                row / 4
            ),
        };
        bytes.extend_from_slice(line.as_bytes());
    }
    bytes
}

/// The tuple file and the counts of `bytes` read in blocks of `block` on `threads` threads.
fn tuples(bytes: &[u8], kept: Option<&[&str]>, block: usize, threads: usize) -> (Vec<u8>, String) {
    let kept: Option<Vec<Vec<u8>>> =
        kept.map(|names| names.iter().map(|name| name.as_bytes().to_vec()).collect());
    let map = read_field_map("map.tsv", MAP.as_bytes(), kept.as_deref()).unwrap();
    let file = Scratch::new(&format!("tuples-in-{block}-{threads}"), bytes);
    let out = Scratch::new(
        &format!("tuples-out-{block}-{threads}"),
        b"an earlier run\n",
    );
    let mut input = file.input(block);
    let header = input.header().unwrap();
    let placed = place_columns(header.as_deref(), &map, kept.is_some()).unwrap();
    let counts = write_tuples(input, &placed, File::create(&out.0).unwrap(), threads).unwrap();
    (std::fs::read(&out.0).unwrap(), summary(&placed, &counts))
}

#[test]
fn the_rows_become_tuples_and_the_summary_counts_them() {
    let bytes = table(11);
    let (written, line) = tuples(&bytes, None, 1 << 20, 1);
    assert_eq!(
        String::from_utf8(written).unwrap(),
        "21:0\tG\tENST0\tA_ONE\tv0\n\
         21:0\tG\tENST0\tB_SCORE\t0\n\
         21:0\tG\tENST1\tA_ONE\tv0\n\
         21:0\tG\tENST1\tB_SCORE\t0\n\
         21:0\tG\tENST2\tA_ONE\tv0\n\
         21:0\tG\tENST2\tB_SCORE\t0\n\
         21:1\tG\tENST4\tA_ONE\tv1\n\
         21:1\tG\tENST4\tB_SCORE\t1\n\
         21:1\tG\tENST6\tA_ONE\tv1\n\
         21:1\tG\tENST6\tB_SCORE\t1\n\
         21:7\tT\tENST7\tA_ONE\ta b\n\
         21:7\tT\tENST7\tA_LIST\tx|y\n\
         21:2\tG\tENST8\tA_ONE\tv2\n\
         21:2\tG\tENST8\tB_SCORE\t2\n\
         21:9\t-\t-\tA_ONE\tone\n\
         21:9\t-\t-\tA_ONE\ttwo\n\
         21:9\t-\t-\tA_LIST\tz\n\
         21:9\t-\t-\tB_WORD\tw\n\
         21:2\tG\tENST10\tA_ONE\tv2\n\
         21:2\tG\tENST10\tB_SCORE\t2\n"
    );
    assert_eq!(
        line,
        "fastvep_plugin_tuples: wrote 20 tuples from 10 rows (1 too short); by field: \
         A_ONE=10, A_LIST=2, B_SCORE=7, B_WORD=1; mapped columns: FV_A, FV_B; dropped parts \
         with no Ensembl VEP field: 9 (FV_A:3=1, FV_B:1.word=7, FV_B:2.score=1); values \
         outside --fields: 0; parts not written word(score): 1; cells with several entries: 1; \
         unmapped columns: FV_PHYLOP"
    );
    let (written, line) = tuples(&bytes, Some(&["B_WORD", "A_LIST"]), 1 << 20, 1);
    assert_eq!(
        String::from_utf8(written).unwrap(),
        "21:7\tT\tENST7\tA_LIST\tx|y\n21:9\t-\t-\tA_LIST\tz\n21:9\t-\t-\tB_WORD\tw\n"
    );
    assert_eq!(
        line,
        "fastvep_plugin_tuples: wrote 3 tuples from 10 rows (1 too short); by field: \
         A_LIST=2, B_WORD=1; mapped columns: FV_A, FV_B; dropped parts with no Ensembl VEP \
         field: 9 (FV_A:3=1, FV_B:1.word=7, FV_B:2.score=1); values outside --fields: 17; \
         parts not written word(score): 1; cells with several entries: 1; unmapped columns: \
         FV_PHYLOP"
    );
}

#[test]
fn blocks_of_any_size_on_any_number_of_threads_write_the_same_tuples_and_counts() {
    let bytes = table(3000);
    let expected = tuples(&bytes, None, 1 << 20, 1);
    assert_eq!(
        expected.0.iter().filter(|byte| **byte == b'\n').count(),
        5452
    );
    for (block, threads) in [
        (1, 1),
        (64, 1),
        (64, 2),
        (500, 3),
        (4096, 8),
        (37, 16),
        (1 << 20, 16),
    ] {
        assert!(
            tuples(&bytes, None, block, threads) == expected,
            "blocks of {block} on {threads} threads"
        );
    }
    // An input without rows writes an empty tuple file on any number of threads.
    let header = &bytes[..bytes.iter().rposition(|byte| *byte == b'B').unwrap() + 2];
    let (written, line) = tuples(header, None, 16, 4);
    assert!(written.is_empty());
    assert!(line.starts_with("fastvep_plugin_tuples: wrote 0 tuples from 0 rows (0 too short); "));
}

#[test]
fn a_summary_without_drops_or_unmapped_columns_says_so() {
    let (_, line) = tuples(
        b"#v\tLocation\tAllele\tGene\tFeature\tFV_B\nv\t21:1\tG\tg\tENST1\tG|(1)|()\n",
        None,
        1 << 20,
        1,
    );
    assert_eq!(
        line,
        "fastvep_plugin_tuples: wrote 1 tuples from 1 rows (0 too short); by field: \
         B_SCORE=1, B_WORD=0; mapped columns: FV_B; dropped parts with no Ensembl VEP field: 0; \
         values outside --fields: 0; parts not written word(score): 0; cells with several \
         entries: 0; unmapped columns: none"
    );
}

#[test]
fn the_thread_bound_is_a_whole_number_of_one_or_more() {
    let bound = |value: &str| thread_bound(Some(std::ffi::OsStr::new(value)));
    assert_eq!(bound("3"), Ok(3));
    assert_eq!(bound(" 12 "), Ok(12));
    assert!(thread_bound(None).unwrap() >= 1);
    assert_eq!(bound(""), thread_bound(None));
    for refused in ["0", "-1", "two", "1.5"] {
        assert_eq!(
            bound(refused),
            Err(format!(
                "VEP_FASTVEP_PLUGIN_TUPLES_THREADS={refused}: expected a whole number of \
                 threads, 1 or more"
            ))
        );
    }
}
