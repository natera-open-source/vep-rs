// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A table dump read row by row as bytes, keeping the values of the columns
//! asked for and passing over every other.
//!
//! The reader of the largest tables of a variation database: one row is in
//! memory at a time, and of that row only the kept values, each in a buffer
//! the next row reuses. The values are read with the escapes `LOAD DATA`
//! reads: a backslash before `0`, `b`, `n`, `r`, `t` or `Z` is that control
//! byte and before any other byte that byte itself, so a tab, a line end and
//! a backslash inside a value are written behind a backslash and a row spans
//! as many lines as its values hold line ends; a value that is `\N` alone is
//! NULL.

use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;

/// The longest row read, 16 MiB, the cap the line reader of the core dumps
/// puts on a line: a file with no line end in it is refused long before a
/// kept value grows to its size.
const MAX_ROW_BYTES: usize = 16 << 20;

/// The kept values of one row.
pub struct DumpRow<'a> {
    values: &'a [Vec<u8>],
    nulls: &'a [bool],
    /// The row's number in the file, from 1.
    pub number: u64,
    width: usize,
}

impl DumpRow<'_> {
    /// The value of the `nth` kept column, `None` for NULL.
    pub fn value(&self, nth: usize) -> Option<&[u8]> {
        (!self.nulls[nth]).then(|| self.values[nth].as_slice())
    }

    /// The number of columns the row holds, which may pass the schema's.
    pub fn width(&self) -> usize {
        self.width
    }
}

/// Reads the dump at `path`, gzip when its name ends `.gz`, and hands `row`
/// the values of the columns `keep` names (0-based, ascending) of each row,
/// in file order. A row of fewer than `columns` values, the schema's count,
/// is an error naming the file and the row, as is a row past the row cap and
/// a file that ends inside a row. Returns the number of rows read.
pub(crate) fn read_rows(
    path: &Path,
    columns: usize,
    keep: &[usize],
    row: impl FnMut(&DumpRow<'_>) -> Result<()>,
) -> Result<u64> {
    read_rows_up_to(path, columns, keep, MAX_ROW_BYTES, row)
}

/// `read_rows` with a row cap of `max_row_bytes` in place of 16 MiB.
pub fn read_rows_up_to(
    path: &Path,
    columns: usize,
    keep: &[usize],
    max_row_bytes: usize,
    mut row: impl FnMut(&DumpRow<'_>) -> Result<()>,
) -> Result<u64> {
    let file =
        std::fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };

    // The slot of each column up to the last one kept.
    let mut slots: Vec<Option<usize>> = vec![None; keep.last().map_or(0, |last| last + 1)];
    for (slot, column) in keep.iter().enumerate() {
        slots[*column] = Some(slot);
    }
    let mut values: Vec<Vec<u8>> = vec![Vec::new(); keep.len()];
    let mut nulls = vec![false; keep.len()];

    // The row under way: the column the next byte belongs to and its slot.
    let mut column = 0usize;
    let mut slot = slots.first().copied().flatten();
    // A value that so far is the NULL marker and nothing else.
    let mut null = false;
    let mut value_started = false;
    let mut escaped = false;
    let mut row_bytes = 0usize;
    let mut row_number = 1u64;

    loop {
        let buffer = reader
            .fill_buf()
            .with_context(|| format!("Failed to read {}", path.display()))?;
        if buffer.is_empty() {
            break;
        }
        let read = buffer.len();
        let mut at = 0;
        while at < read {
            // A run of plain bytes is taken whole.
            let plain = if escaped {
                0
            } else {
                buffer[at..]
                    .iter()
                    .position(|byte| matches!(byte, b'\\' | b'\t' | b'\n'))
                    .unwrap_or(read - at)
            };
            if plain > 0 {
                if let Some(slot) = slot {
                    if null {
                        // `\N` followed by more is the letter, not the marker.
                        values[slot].push(b'N');
                    }
                    values[slot].extend_from_slice(&buffer[at..at + plain]);
                }
                null = false;
                value_started = true;
                row_bytes += plain;
                at += plain;
                continue;
            }
            let byte = buffer[at];
            at += 1;
            row_bytes += 1;
            if escaped {
                escaped = false;
                if byte == b'N' && !value_started {
                    null = true;
                } else {
                    let decoded = match byte {
                        b'0' => 0,
                        b'b' => 8,
                        b'n' => b'\n',
                        b'r' => b'\r',
                        b't' => b'\t',
                        b'Z' => 26,
                        other => other,
                    };
                    if let Some(slot) = slot {
                        if null {
                            values[slot].push(b'N');
                        }
                        values[slot].push(decoded);
                    }
                    null = false;
                }
                value_started = true;
            } else if byte == b'\\' {
                escaped = true;
            } else {
                // The end of a value: a tab, or the row's line end.
                if let Some(slot) = slot {
                    nulls[slot] = null;
                }
                null = false;
                value_started = false;
                if byte == b'\t' {
                    column += 1;
                    slot = slots.get(column).copied().flatten();
                } else {
                    if column + 1 < columns {
                        bail!(
                            "{} row {row_number} has {} columns, the schema has {columns}",
                            path.display(),
                            column + 1
                        );
                    }
                    row(&DumpRow {
                        values: &values,
                        nulls: &nulls,
                        number: row_number,
                        width: column + 1,
                    })?;
                    for value in &mut values {
                        value.clear();
                    }
                    nulls.fill(false);
                    column = 0;
                    slot = slots.first().copied().flatten();
                    row_bytes = 0;
                    row_number += 1;
                }
            }
        }
        if row_bytes > max_row_bytes {
            bail!(
                "{} row {row_number} is longer than the {max_row_bytes}-byte limit for one row",
                path.display()
            );
        }
        reader.consume(read);
    }
    if column != 0 || value_started || escaped || null {
        bail!(
            "{} ends inside row {row_number}, which has no line end",
            path.display()
        );
    }
    Ok(row_number - 1)
}

/// `value` as an unsigned 32-bit number, the type of an identifier column:
/// `None` for anything but ASCII digits that fit it.
pub(crate) fn unsigned(value: &[u8]) -> Option<u32> {
    if value.is_empty() {
        return None;
    }
    value.iter().try_fold(0u32, |number, byte| {
        let digit = byte.checked_sub(b'0').filter(|digit| *digit < 10)?;
        number.checked_mul(10)?.checked_add(u32::from(digit))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The rows `read_rows` hands over from `bytes` written as a plain dump:
    /// each kept value as text, `None` for NULL.
    fn rows_of(bytes: &[u8], columns: usize, keep: &[usize]) -> Result<Vec<Vec<Option<String>>>> {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("table.txt");
        std::fs::write(&path, bytes).unwrap();
        let mut rows = Vec::new();
        let read = read_rows(&path, columns, keep, |row| {
            rows.push(
                (0..keep.len())
                    .map(|nth| {
                        row.value(nth)
                            .map(|value| String::from_utf8(value.to_vec()).unwrap())
                    })
                    .collect(),
            );
            Ok(())
        })?;
        assert_eq!(read as usize, rows.len());
        Ok(rows)
    }

    fn text(value: &str) -> Option<String> {
        Some(value.to_string())
    }

    /// A value holds what its escapes spell: a tab and a line end written
    /// behind a backslash stay inside the value, so the row that holds the
    /// line end spans two lines and is one row; `\N` alone is NULL, while
    /// `\N` before more bytes and a plain `N` are the letter. A column that
    /// is not kept is passed over, its escaped tab and line end counting as
    /// neither a column end nor a row end.
    #[test]
    fn a_row_is_read_by_its_escapes_and_a_null_is_told_from_the_letter() {
        let dump = b"7\tfirst\\\ttabbed\t\\N\tkept\n\
            8\ttwo\\\nlines\tN\tplain\n\
            9\t\\Nine\\\\\\0\\Z\tskipped\\\tfield\\\nhere\t\\N\n";
        assert_eq!(
            rows_of(dump, 4, &[0, 1, 2, 3]).unwrap(),
            [
                vec![text("7"), text("first\ttabbed"), None, text("kept")],
                vec![text("8"), text("two\nlines"), text("N"), text("plain")],
                vec![
                    text("9"),
                    text("Nine\\\0\u{1a}"),
                    text("skipped\tfield\nhere"),
                    None
                ],
            ]
        );
        assert_eq!(
            rows_of(dump, 4, &[1, 3]).unwrap(),
            [
                vec![text("first\ttabbed"), text("kept")],
                vec![text("two\nlines"), text("plain")],
                vec![text("Nine\\\0\u{1a}"), None],
            ]
        );
    }

    /// A gzip dump of several members reads as the plain one, an empty dump
    /// holds no row, and a row with more values than the schema names is
    /// read by its leading columns.
    #[test]
    fn a_compressed_dump_an_empty_one_and_a_wider_row_are_read() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("table.txt.gz");
        let mut bytes = Vec::new();
        for member in ["1\ta\n2\t", "b\n3\tc\tmore\n"] {
            let mut encoder =
                flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            std::io::Write::write_all(&mut encoder, member.as_bytes()).unwrap();
            bytes.extend(encoder.finish().unwrap());
        }
        std::fs::write(&path, bytes).unwrap();
        let mut rows = Vec::new();
        let read = read_rows(&path, 2, &[1], |row| {
            rows.push((row.number, row.value(0).unwrap().to_vec()));
            Ok(())
        })
        .unwrap();
        assert_eq!(read, 3);
        assert_eq!(
            rows,
            [(1, b"a".to_vec()), (2, b"b".to_vec()), (3, b"c".to_vec())]
        );
        assert!(rows_of(b"", 3, &[0]).unwrap().is_empty());
    }

    /// A row short of the schema's columns and a file that ends inside a row
    /// are refused, each naming the row.
    #[test]
    fn a_short_row_and_an_unended_row_are_refused() {
        let short = rows_of(b"1\ta\tb\n2\ta\n", 3, &[0]).unwrap_err();
        assert!(
            format!("{short:#}").ends_with("table.txt row 2 has 2 columns, the schema has 3"),
            "{short:#}"
        );
        for unended in [
            &b"1\ta\tb\n2\ta\tb"[..],
            b"1\ta\tb\n2\ta\tb\\",
            b"1\ta\t\\N",
        ] {
            let error = rows_of(unended, 3, &[0]).unwrap_err();
            assert!(
                format!("{error:#}").contains("ends inside row"),
                "{error:#}"
            );
        }
    }

    /// An identifier is its digits; a sign, a space, an empty value and a
    /// number past 32 bits are none.
    #[test]
    fn an_identifier_is_ascii_digits_that_fit_32_bits() {
        assert_eq!(unsigned(b"0"), Some(0));
        assert_eq!(unsigned(b"4294967295"), Some(u32::MAX));
        for not_one in [&b""[..], b"4294967296", b"-1", b"1 ", b"1.0", b"\\N"] {
            assert_eq!(unsigned(not_one), None, "{not_one:?}");
        }
    }
}
