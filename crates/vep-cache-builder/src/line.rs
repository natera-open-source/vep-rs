// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One line of a GFF3, a GTF or a regulatory GFF, as every reader of this
//! crate takes it.
//!
//! A line is read without its line end, under a length limit
//! (`MAX_LINE_BYTES`), and the byte-order marks that open it are taken off
//! (`without_byte_order_marks`). Its columns have two readings. The index
//! build of Ensembl VEP's release ends a column at a tab or at a NUL byte
//! (htslib 1.9 `tbx.c:89`, `tbx_parse1`): that reading decides whether the
//! build places the line, which it does when the line has four columns
//! (`index_columns`) and a start and an end the build parses
//! (`check_indexed_coordinates`). An annotation reader ends a column at a
//! tab alone (ensembl-io `ColumnBasedParser.pm:250`, `read_record`): that
//! reading decides what a record holds (`tab_columns`). A line whose
//! columns NUL bytes end is one the index build places and one in which an
//! annotation reader finds fewer columns than a record has, so a reader
//! passes over it as over any placed line of fewer than nine columns.

use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
pub(crate) use vep_core::regulatory_gff::index_columns;

use crate::gff3::check_skipped_coordinates;

/// The longest line read, 16 MiB. A feature line's widest column is its
/// attributes, a few kilobytes even for a gene with every cross-reference
/// listed, so the cap refuses a file with no line end in it long before the
/// line buffer grows to its size.
pub(crate) const MAX_LINE_BYTES: usize = 16 << 20;

/// A line as the file writes it, without its line end and without the
/// byte-order marks that open it ([`without_byte_order_marks`]): its first
/// `MAX_LINE_BYTES` bytes, and whether those are the whole line.
pub(crate) struct Line {
    pub(crate) written: Vec<u8>,
    pub(crate) whole: bool,
}

/// The next line behind `reader`, none at the file's end. The rest of a
/// line longer than `MAX_LINE_BYTES` is left unread.
pub(crate) fn read_line(reader: &mut impl BufRead, path: &Path) -> Result<Option<Line>> {
    // The cap plus the line end: a line of exactly the cap is read whole.
    let limit = MAX_LINE_BYTES as u64 + 1;
    let mut written = Vec::new();
    let read = (&mut *reader)
        .take(limit)
        .read_until(b'\n', &mut written)
        .with_context(|| format!("Failed to read {}", path.display()))?;
    if read == 0 {
        return Ok(None);
    }
    let ended = written.last() == Some(&b'\n');
    if ended {
        written.pop();
    }
    let whole = ended || read as u64 != limit;
    written.truncate(MAX_LINE_BYTES);
    Ok(Some(Line {
        written: without_byte_order_marks(written),
        whole,
    }))
}

/// [`read_line`] where a line may be a sequence's, of any length: the rest
/// of a line longer than `MAX_LINE_BYTES` is passed over.
pub(crate) fn read_line_of_any_length(
    reader: &mut impl BufRead,
    path: &Path,
) -> Result<Option<Line>> {
    let line = read_line(reader, path)?;
    if line.as_ref().is_some_and(|line| !line.whole) {
        reader
            .skip_until(b'\n')
            .with_context(|| format!("Failed to read {}", path.display()))?;
    }
    Ok(line)
}

/// The lines of the file at `path`, read off `reader` one at a time
/// ([`read_line`]): a line longer than `MAX_LINE_BYTES` is an error naming
/// the file and the line number, as is a read failure, and either ends the
/// iteration.
pub(crate) struct Lines<R> {
    reader: R,
    path: PathBuf,
    read: usize,
    done: bool,
}

/// The lines of the file at `path`, behind `reader` ([`Lines`]).
pub(crate) fn lines<R: BufRead>(reader: R, path: &Path) -> Lines<R> {
    Lines {
        reader,
        path: path.to_path_buf(),
        read: 0,
        done: false,
    }
}

impl<R: BufRead> Iterator for Lines<R> {
    type Item = Result<Vec<u8>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        self.read += 1;
        let item = match read_line(&mut self.reader, &self.path) {
            Ok(None) => None,
            Ok(Some(line)) if line.whole => Some(Ok(line.written)),
            Ok(Some(_)) => Some(Err(anyhow!(
                "{} line {} is longer than the {MAX_LINE_BYTES}-byte limit for one line",
                self.path.display(),
                self.read
            ))),
            Err(error) => Some(Err(error)),
        };
        self.done = !matches!(item, Some(Ok(_)));
        item
    }
}

/// U+FEFF, the byte-order mark, as UTF-8 writes it.
const BYTE_ORDER_MARK: &[u8] = "\u{feff}".as_bytes();

/// `line` without the byte-order marks that open it, alone or behind white
/// space. A mark at the start of a line changes nothing a reader reads: a
/// file with marks parses as the same file without them.
pub(crate) fn without_byte_order_marks(mut line: Vec<u8>) -> Vec<u8> {
    let lead = line
        .iter()
        .take_while(|byte| byte.is_ascii_whitespace())
        .count();
    let marks = line[lead..]
        .chunks_exact(BYTE_ORDER_MARK.len())
        .take_while(|chunk| *chunk == BYTE_ORDER_MARK)
        .count();
    line.drain(lead..lead + marks * BYTE_ORDER_MARK.len());
    line
}

/// Whether `line` is blank or a comment as a transcript reader reads one:
/// nothing stands in it, or a `#` opens it, once the white space before is
/// dropped.
pub(crate) fn is_blank_or_a_comment(line: &[u8]) -> bool {
    let text = String::from_utf8_lossy(line);
    let text = text.trim_start();
    text.is_empty() || text.starts_with('#')
}

/// Whether `line` is a comment as a transcript reader reads one: a `#` opens
/// it once the white space before it is dropped.
pub(crate) fn is_a_comment(line: &[u8]) -> bool {
    String::from_utf8_lossy(line).trim_start().starts_with('#')
}

/// The columns of `text` as an annotation reader ends them: at tabs alone,
/// so a NUL byte is a byte of the column it stands in (ensembl-io
/// `ColumnBasedParser.pm:250`).
pub(crate) fn tab_columns(text: &str) -> std::str::Split<'_, char> {
    text.split('\t')
}

/// What a transcript reader reads of `written`, a line that is neither
/// blank nor a comment: the line's text and its columns, ended at tabs
/// ([`tab_columns`]). The text is the line without the white space around
/// it, unless a tab opens the line: it is then read at the columns the file
/// writes, its first column empty, as the index build reads it. `Err` is the
/// count of the line's columns as the index build ends them
/// ([`index_columns`]) when that is under four: the line ends before its
/// start, and the index build fails it.
pub(crate) fn record_columns(written: &str) -> std::result::Result<(&str, Vec<&str>), usize> {
    let placed = index_columns(written.as_bytes()).count();
    if placed < 4 {
        return Err(placed);
    }
    let text = if written.starts_with('\t') {
        written.trim_end()
    } else {
        written.trim()
    };
    Ok((text, tab_columns(text).collect()))
}

/// Whether white space and then a tab open `written`, a line as the file
/// writes it: a first column the index build reads and a reader's trim
/// drops, so each column of the record sits one place further on to that
/// build than to the reader.
pub(crate) fn opens_with_a_blank_column(written: &str) -> bool {
    written
        .split_once('\t')
        .is_some_and(|(first, _)| !first.is_empty() && first.trim().is_empty())
}

/// Whether the index build reads the start and the end of `written` in
/// other columns than a reader may: white space and a tab open the line
/// ([`opens_with_a_blank_column`]), or a NUL byte, which ends a column to
/// that build alone, stands in it.
pub(crate) fn index_columns_may_differ(written: &str) -> bool {
    opens_with_a_blank_column(written) || written.contains('\0')
}

/// Holds the `start` and `end` of `written`, a line as the file writes it,
/// to the index build's reading of them: `parse`, the reader's own reading
/// of a column, runs on the first of the two the build fails, so its error
/// refuses the line. The two are the fourth and fifth columns as the index
/// build ends them ([`index_columns`]), and the rule is that of
/// [`check_skipped_coordinates`], whose reading of a line
/// ([`vep_core::regulatory_gff::unindexed_coordinate`]) ends a column at a
/// NUL byte as at a tab: `strtol` reads on past a tab into the next column
/// and stops at a NUL byte (htslib 1.9 `tbx.c:95` and `:105` hand it the
/// column as a C string), so a start or an end whose blanks run into a NUL
/// byte has no number and fails.
pub(crate) fn check_indexed_coordinates(
    written: &str,
    parse: impl Fn(&str, &str) -> Result<u64>,
) -> Result<()> {
    check_skipped_coordinates(written, parse)
}

/// Whether the index build of Ensembl VEP's release could place `line`, a
/// whole line as the file writes it: the line has a start column as that
/// build ends its columns ([`index_columns`]), and the build parses its
/// start and its end ([`check_indexed_coordinates`]). It places every line
/// that reads as a record, so the line that makes a `>` line a header is
/// never a record the section leaves unread.
pub(crate) fn index_build_places(line: &[u8]) -> bool {
    let written = String::from_utf8_lossy(line);
    let fails = |_: &str, _: &str| -> Result<u64> { Err(anyhow!("the index build fails it")) };
    index_columns(line).nth(3).is_some() && check_indexed_coordinates(&written, fails).is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every line of `text`, read through [`lines`].
    fn read(text: &[u8]) -> Vec<Vec<u8>> {
        lines(text, Path::new("t.gff3"))
            .collect::<Result<Vec<_>>>()
            .unwrap()
    }

    /// A line is read without its line end and without the byte-order marks
    /// that open it; a carriage return before the line end is a byte of the
    /// line, a last line needs no line end, and a blank line is a line.
    #[test]
    fn a_line_is_read_without_its_line_end_and_its_byte_order_marks() {
        assert_eq!(
            read(b"\xef\xbb\xbf21\tx\n\n \xef\xbb\xbf\xef\xbb\xbf#c\r\nlast"),
            [&b"21\tx"[..], b"", b" #c\r", b"last"]
        );
        assert!(read(b"").is_empty());
    }

    /// A line longer than the limit ends the iteration with an error naming
    /// the file and the line; one of exactly the limit is read whole.
    #[test]
    fn a_line_longer_than_the_limit_is_an_error_naming_it() {
        let mut text = b"first\n".to_vec();
        text.resize(text.len() + MAX_LINE_BYTES, b'x');
        text.extend_from_slice(b"\nlast\n");
        assert_eq!(read(&text).len(), 3);

        text.insert(6, b'x');
        let mut read = lines(&text[..], Path::new("t.gff3"));
        assert_eq!(read.next().unwrap().unwrap(), b"first");
        assert_eq!(
            read.next().unwrap().unwrap_err().to_string(),
            format!("t.gff3 line 2 is longer than the {MAX_LINE_BYTES}-byte limit for one line")
        );
        assert!(read.next().is_none());
    }

    /// The index build ends a column at a tab and at a NUL byte, an
    /// annotation reader at a tab alone: a line whose columns NUL bytes end
    /// is nine columns to the first and one to the second, and a NUL byte
    /// inside a tab-ended column is one more column to the index build.
    #[test]
    fn a_column_ends_at_a_tab_or_a_nul_byte_to_the_index_build_and_at_a_tab_to_a_reader() {
        let columns = |line: &str| {
            (
                index_columns(line.as_bytes()).count(),
                tab_columns(line).count(),
            )
        };
        for (line, expected) in [
            ("21\tx\texon\t1000\t1100\t.\t+\t.\tParent=T", (9, 9)),
            ("21\0x\0exon\x001000\x001100\0.\0+\0.\0Parent=T", (9, 1)),
            ("21\tx\texon\t1000\t1100\t.\t+\t.\tParent=T\0", (10, 9)),
            ("21\tx\tex\0on\t1000\t1100\t.\t+\t.\tParent=T", (10, 9)),
            ("\t21\tx\texon\t1000", (5, 5)),
            ("21\tx\texon\t1000\t", (5, 5)),
            ("21\tx\texon\t1000\0", (5, 4)),
            ("21\tx\texon\t1000\t1100\r", (5, 5)),
            ("", (1, 1)),
        ] {
            assert_eq!(columns(line), expected, "{line:?}");
        }
        let parts: Vec<&[u8]> = index_columns(b"a\0\tb\0").collect();
        assert_eq!(parts, [&b"a"[..], b"", b"b", b""]);
    }

    /// A transcript reader reads a line's columns at tabs, without the white
    /// space around the line unless a tab opens it, when the first column is
    /// empty and a blank last column is kept off; a line of under four
    /// columns as the index build ends them gives that count, so one whose
    /// columns NUL bytes end is read, as one column.
    #[test]
    fn a_reader_takes_the_columns_of_a_line_the_index_build_gives_four() {
        let columns = |written| record_columns(written).map(|(_, columns)| columns);
        assert_eq!(
            record_columns("  21\tx\texon\t1000\t1100 \r"),
            Ok((
                "21\tx\texon\t1000\t1100",
                vec!["21", "x", "exon", "1000", "1100"]
            ))
        );
        assert_eq!(
            columns("\t21\tx\texon\t1000\t "),
            Ok(vec!["", "21", "x", "exon", "1000"])
        );
        assert_eq!(
            columns("21\tx\texon\t1000\t"),
            Ok(vec!["21", "x", "exon", "1000"])
        );
        assert_eq!(
            columns("21\0x\0exon\x001000\x001100"),
            Ok(vec!["21\0x\0exon\x001000\x001100"])
        );
        assert_eq!(columns("21\tx\texon"), Err(3));
        assert_eq!(columns("21\0x\texon"), Err(3));
        assert_eq!(columns("21 x exon 1000 1100"), Err(1));
    }

    /// The index build places a line whose columns NUL bytes end as it
    /// places the line with tabs in their place: four columns, a start a
    /// number leads and an end a number leads that does not read below zero.
    /// A carriage return that ends the line is no part of its last column.
    #[test]
    fn the_index_build_places_a_line_whose_columns_nul_bytes_end() {
        for placed in [
            &b"21\0x\0exon\x001000\x001100\0.\0+\0.\0Parent=T"[..],
            b">21\0x\0exon\x001000\x001100\0.\0+\0.\0Parent=T",
            b"21\0x\0exon\x001000",
            b"21\tx\0exon\t1000\x001100\r",
            b"21\tx\texon\t1000\t1100\t.\t+\t.\tParent=T\0",
            b"21\tx\texon\t10\x0000\t1100\t.\t+\t.\tParent=T",
            b"21\tx\texon\t1000\0\t5",
            b"\0\0\x007",
        ] {
            assert!(index_build_places(placed), "{placed:?}");
        }
        for failed in [
            &b"21\0x\0exon"[..],
            b"21\0x\0exon\0",
            b"21\0x\0exon\0abc\x001100",
            b"21\0x\0exon\x001000\0-9",
            b"21\tx\0y\texon\t1000\t1100\t.\t+\t.\tParent=T",
            b"21\tx\tex\0on\t1000\t1100\t.\t+\t.\tParent=T",
        ] {
            assert!(!index_build_places(failed), "{failed:?}");
        }
    }

    /// `strtol` reads on past the tab of an empty or blank column into the
    /// next and stops at a NUL byte, so an empty start before a numeric end
    /// is placed when a tab ends it and failed when a NUL byte does, with
    /// the reader's own error for the column: the blanks before the NUL byte
    /// are the column's value.
    #[test]
    fn a_blank_column_a_nul_byte_ends_has_no_number() {
        let parse = |column: &str, value: &str| -> Result<u64> {
            Err(anyhow!("{column} {value:?} is not a number"))
        };
        for placed in [
            "21\tx\texon\t\t1100\t.",
            "21\0x\0exon\0\t1100\0.",
            "21\tx\texon\t1000\t\t5",
        ] {
            assert!(index_build_places(placed.as_bytes()), "{placed:?}");
            check_indexed_coordinates(placed, parse).unwrap();
        }
        for (failed, error) in [
            ("21\tx\texon\t\x001100\t.", "start \"\" is not a number"),
            ("21\tx\texon\t \t\x001100", "start \" \" is not a number"),
            ("21\tx\texon\t1000\t \x005", "end \" \" is not a number"),
            ("21\tx\texon\t1000\0\0", "end \"\" is not a number"),
        ] {
            assert!(!index_build_places(failed.as_bytes()), "{failed:?}");
            assert_eq!(
                check_indexed_coordinates(failed, parse)
                    .unwrap_err()
                    .to_string(),
                error,
                "{failed:?}"
            );
        }
    }

    /// The index build may read a line's start and end in other columns than
    /// a reader when white space and a tab open the line or a NUL byte
    /// stands in it, and in the reader's columns otherwise, a line a tab
    /// opens among them.
    #[test]
    fn the_index_build_reads_other_columns_behind_a_blank_column_or_a_nul_byte() {
        for line in [" \t21\tx", "\r\t21", "21\tx\0", "2\x001\tx\texon"] {
            assert!(index_columns_may_differ(line), "{line:?}");
        }
        for line in ["21\tx\texon", "\t21\tx", " 21\tx", "", "x\t \ty"] {
            assert!(!index_columns_may_differ(line), "{line:?}");
        }
    }
}
