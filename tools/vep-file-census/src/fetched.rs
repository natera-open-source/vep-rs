// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Fetched files held to a list of key, bytes and ETag: every file under a root against
//! the size and the ETag an object store states for its object.
//!
//! The list is tab-separated lines `key<TAB>bytes<TAB>etag` without a header; a blank
//! line and a line that starts with `#` are skipped. An ETag is 32 hexadecimal digits,
//! the md5 of the file, or `<32 hex>-<parts>` for an object uploaded in parts: the md5 of
//! the parts' binary md5s joined, then the count of parts. A third cell of 64 hexadecimal
//! digits is the file's sha256.
//!
//! The verdicts, one a key, for the file `<root>/<key>`: `md5`, `sha256`,
//! `multipart-md5 <MiB>` (the ETag is that of its upload in parts of that size),
//! `bytes-alone` (the size holds and no part size tried gave the ETag) and `DIFFERS`
//! (missing, unreadable, of another size or of another digest, or an ETag of neither
//! shape). The part sizes tried are the whole MiB for which the size makes the stated
//! count of parts: 8, 16, 32, 64, 128, 256, 512 and 1024 MiB first, then every other from
//! 5 to 5120 MiB upwards, twelve at most, one alone where the object is one part.

use crate::pipeline::ordered;
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

const MIB: u128 = 1 << 20;
const USUAL_PART_MIB: [u128; 8] = [8, 16, 32, 64, 128, 256, 512, 1024];
const MOST_PART_SIZES: usize = 12;

pub const DIFFERS: &str = "DIFFERS";
pub const BYTES_ALONE: &str = "bytes-alone";

/// One line of the list: the key, the size as the report writes it and as a number where
/// it fits one, and the ETag.
pub struct Row {
    pub key: String,
    pub size_text: String,
    size: Option<u128>,
    etag: String,
}

fn is_hex(text: &str, length: usize) -> bool {
    text.len() == length && text.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f'))
}

/// The lines of `text` as `str.splitlines()` cuts them.
fn split_lines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut chars = text.char_indices().peekable();
    while let Some((at, c)) = chars.next() {
        let ends = matches!(
            c,
            '\n' | '\r'
                | '\u{b}'
                | '\u{c}'
                | '\u{1c}'
                | '\u{1d}'
                | '\u{1e}'
                | '\u{85}'
                | '\u{2028}'
                | '\u{2029}'
        );
        if !ends {
            continue;
        }
        lines.push(&text[start..at]);
        let mut next = at + c.len_utf8();
        if c == '\r' && chars.peek().is_some_and(|(_, c)| *c == '\n') {
            chars.next();
            next += 1;
        }
        start = next;
    }
    if start < text.len() {
        lines.push(&text[start..]);
    }
    lines
}

/// Whether `str.isspace()` holds for the character.
fn is_py_space(c: char) -> bool {
    matches!(
        c,
        '\u{9}'..='\u{d}'
            | '\u{1c}'..='\u{20}'
            | '\u{85}'
            | '\u{a0}'
            | '\u{1680}'
            | '\u{2000}'..='\u{200a}'
            | '\u{2028}'
            | '\u{2029}'
            | '\u{202f}'
            | '\u{205f}'
            | '\u{3000}'
    )
}

/// The rows of the list at `path`, or the words it cannot be taken for.
pub fn listed(path: &Path) -> Result<Vec<Row>, String> {
    let bytes = std::fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let text =
        String::from_utf8(bytes).map_err(|_| format!("{} is not UTF-8 text", path.display()))?;
    let mut rows = Vec::new();
    for (number, line) in split_lines(&text).into_iter().enumerate() {
        if line.trim_matches(is_py_space).is_empty() || line.starts_with('#') {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        if cells.len() != 3 {
            return Err(format!(
                "line {} of {} holds {} columns, not key, bytes and etag",
                number + 1,
                path.display(),
                cells.len()
            ));
        }
        if cells[1].is_empty() || !cells[1].bytes().all(|b| b.is_ascii_digit()) {
            return Err(format!(
                "line {} of {} states bytes that are no number: {}",
                number + 1,
                path.display(),
                cells[1]
            ));
        }
        let digits = cells[1].trim_start_matches('0');
        rows.push(Row {
            key: cells[0].to_string(),
            size_text: if digits.is_empty() { "0" } else { digits }.to_string(),
            size: digits.parse().ok().or(digits.is_empty().then_some(0)),
            etag: cells[2].to_string(),
        });
    }
    Ok(rows)
}

/// The part sizes to try, in MiB and in the order tried, for `size` bytes uploaded in
/// `parts` parts.
fn part_sizes(size: u128, parts: u128) -> Vec<u128> {
    let fits = |mib: &u128| size.div_ceil(mib * MIB) == parts;
    let mut ordered: Vec<u128> = USUAL_PART_MIB.iter().copied().filter(fits).collect();
    ordered.extend(
        (5..=5120)
            .filter(|mib| !USUAL_PART_MIB.contains(mib))
            .filter(fits),
    );
    // One part is the whole file under every part size that fits.
    ordered.truncate(if parts == 1 { 1 } else { MOST_PART_SIZES });
    ordered
}

fn md5_of(path: &Path) -> io::Result<String> {
    let mut file = File::open(path)?;
    let mut state = md5::Context::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        match file.read(&mut buf) {
            Ok(0) => return Ok(format!("{:x}", state.finalize())),
            Ok(count) => state.consume(&buf[..count]),
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
}

fn sha256_of(path: &Path) -> io::Result<String> {
    let mut file = File::open(path)?;
    let mut state = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        match file.read(&mut buf) {
            Ok(0) => {
                return Ok(state
                    .finalize()
                    .iter()
                    .map(|byte| format!("{byte:02x}"))
                    .collect())
            }
            Ok(count) => state.update(&buf[..count]),
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
}

/// The ETag an object store states for the file uploaded in parts of `part_bytes` bytes.
fn etag_of_parts(path: &Path, part_bytes: u128) -> io::Result<String> {
    let mut file = File::open(path)?;
    let mut buf = vec![0u8; 1 << 20];
    let mut joined = md5::Context::new();
    let mut parts = 0u64;
    loop {
        let mut part = md5::Context::new();
        let mut left = part_bytes;
        while left > 0 {
            let want = (buf.len() as u128).min(left) as usize;
            let count = match file.read(&mut buf[..want]) {
                Ok(count) => count,
                Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(error),
            };
            if count == 0 {
                break;
            }
            part.consume(&buf[..count]);
            left -= count as u128;
        }
        if left == part_bytes {
            break;
        }
        joined.consume(part.finalize().0);
        parts += 1;
    }
    Ok(format!("{:x}-{parts}", joined.finalize()))
}

fn verdict(root: &Path, row: &Row) -> String {
    let path: PathBuf = root.join(&row.key);
    let checked = (|| -> io::Result<String> {
        let meta = std::fs::metadata(&path)?;
        if !meta.is_file() || Some(u128::from(meta.len())) != row.size {
            return Ok(DIFFERS.to_string());
        }
        let size = u128::from(meta.len());
        let etag = row.etag.as_str();
        if is_hex(etag, 32) {
            return Ok(if md5_of(&path)? == etag {
                "md5"
            } else {
                DIFFERS
            }
            .to_string());
        }
        if is_hex(etag, 64) {
            return Ok(if sha256_of(&path)? == etag {
                "sha256"
            } else {
                DIFFERS
            }
            .to_string());
        }
        let Some((digest, parts)) = etag.split_once('-') else {
            return Ok(DIFFERS.to_string());
        };
        let counted = !parts.is_empty()
            && !parts.starts_with('0')
            && parts.bytes().all(|b| b.is_ascii_digit());
        if !is_hex(digest, 32) || !counted {
            return Ok(DIFFERS.to_string());
        }
        // A count past what a number holds is the count of no part size.
        let Ok(parts) = parts.parse::<u128>() else {
            return Ok(BYTES_ALONE.to_string());
        };
        for mib in part_sizes(size, parts) {
            if etag_of_parts(&path, mib * MIB)? == etag {
                return Ok(format!("multipart-md5 {mib}"));
            }
        }
        Ok(BYTES_ALONE.to_string())
    })();
    checked.unwrap_or_else(|_| DIFFERS.to_string())
}

/// The verdict of every row, in the list's order, `jobs` files at a time.
pub fn verdicts(root: &Path, rows: &[Row], jobs: usize) -> io::Result<Vec<String>> {
    let mut all = Vec::with_capacity(rows.len());
    ordered(
        jobs.min(rows.len().max(1)),
        |emit| {
            for row in rows {
                if !emit(row) {
                    break;
                }
            }
        },
        |row: &Row| verdict(root, row),
        |word| {
            all.push(word);
            Ok(())
        },
    )?;
    Ok(all)
}

/// The closing line: how many keys each verdict holds.
pub fn held_line(verdicts: &[String]) -> String {
    let count = |word: &str| {
        verdicts
            .iter()
            .filter(|verdict| verdict.split(' ').next() == Some(word))
            .count()
    };
    let by_sha256 = match count("sha256") {
        0 => String::new(),
        held => format!("{held} by sha256, "),
    };
    format!(
        "held: {} by md5, {} by multipart md5, {} by bytes alone, {by_sha256}{} differ of {}",
        count("md5"),
        count("multipart-md5"),
        count(BYTES_ALONE),
        count(DIFFERS),
        verdicts.len()
    )
}

#[cfg(test)]
#[path = "fetched_tests.rs"]
mod tests;
