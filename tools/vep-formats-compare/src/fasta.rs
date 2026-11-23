// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reference sequence read through its index.

use crate::argparse::path_name;
use crate::argparse::path_text;
use crate::error::Error;
use crate::error::Result;
use crate::pyfmt;
use crate::report::strerror;
use indexmap::IndexMap;
use std::io;
use std::io::Read;
use std::io::Seek;
use std::path::Path;
use std::path::PathBuf;

/// The error number of a read from a directory.
const EISDIR: i32 = 21;

/// The reference FASTA `hgvs_shift_sequence_edge_window` reads: an uncompressed FASTA
/// beside its `.fai` index, a region's length and its bases by 1-based inclusive
/// coordinates, upper-cased and clipped to the region.
#[derive(Debug)]
pub struct ReferenceSequence {
    /// A region's length, the offset of its first base, the bases a line holds and the bytes
    /// a line holds.
    pub regions: IndexMap<String, (i64, i64, i64, i64)>,
    pub handle: std::fs::File,
}

/// What the system says of an error, as the interpreter prints an `OSError`.
fn os_error(error: &io::Error) -> String {
    match error.raw_os_error() {
        Some(code) => format!("[Errno {code}] {}", strerror(error)),
        None => strerror(error),
    }
}

/// `path` opened for reading; the interpreter refuses to open a directory.
fn opened(path: &Path) -> io::Result<std::fs::File> {
    let handle = std::fs::File::open(path)?;
    if handle.metadata()?.is_dir() {
        return Err(io::Error::from_raw_os_error(EISDIR));
    }
    Ok(handle)
}

impl ReferenceSequence {
    pub fn new(path: &Path) -> Result<Self> {
        // `path.suffix`: the name's text from its last dot, where a character other than a
        // dot stands before that dot.
        let compressed = path_name(path)
            .strip_suffix(".gz")
            .is_some_and(|stem| stem.chars().any(|c| c != '.'));
        if compressed {
            return Err(Error::comparator(format!(
                "--reference-fasta: {} is compressed; an uncompressed FASTA beside its .fai is read",
                path_text(path)
            )));
        }
        let mut index = path.as_os_str().to_os_string();
        index.push(".fai");
        let index = PathBuf::from(index);
        let cannot_read = |path: &Path, error: &io::Error| {
            Error::comparator(format!(
                "--reference-fasta: cannot read {}: {}",
                path_text(path),
                strerror(error)
            ))
        };
        let mut bytes = Vec::new();
        opened(&index)
            .and_then(|mut handle| handle.read_to_end(&mut bytes))
            .map_err(|error| cannot_read(&index, &error))?;
        // The system's errors alone are the comparator's error: an index that is not UTF-8
        // ends the run.
        let text = match pyfmt::decode_utf8(&bytes) {
            Ok(text) => text,
            Err(error) => crate::uncaught!("{error}"),
        };
        let handle = opened(path).map_err(|error| cannot_read(path, &error))?;
        let mut regions = IndexMap::new();
        for line in pyfmt::split_lines(text) {
            let cols: Vec<&str> = line.split('\t').collect();
            let refused = || {
                Error::comparator(format!(
                    "{}: not a FASTA index line: {line}",
                    path_text(&index)
                ))
            };
            if cols.len() < 5 || !cols[1..5].iter().all(|col| pyfmt::is_digit(col)) {
                return Err(refused());
            }
            // `str.isdigit` takes digits `int` does not read (a superscript), and `int`
            // refuses a number longer than the interpreter converts: either ends the run.
            // A number `i64` does not hold is no index line.
            let whole = |digits: &str| -> Result<i64> {
                match pyfmt::int_checked(digits) {
                    Ok(number) => number.parse().map_err(|_| refused()),
                    Err(error) => crate::uncaught!("{error}"),
                }
            };
            let per_line = whole(cols[3])?;
            if per_line == 0 {
                return Err(refused());
            }
            let entry = (whole(cols[1])?, whole(cols[2])?, per_line, whole(cols[4])?);
            regions.insert(cols[0].to_string(), entry);
        }
        Ok(ReferenceSequence { regions, handle })
    }

    pub fn length(&self, region: &str) -> Option<i64> {
        self.regions.get(region).map(|entry| entry.0)
    }

    pub fn bases(&self, region: &str, start: i64, end: i64) -> String {
        let Some(&(length, offset, per_line, width)) = self.regions.get(region) else {
            crate::uncaught!("{}", pyfmt::repr_str(region))
        };
        let (start, end) = (start.max(1), end.min(length));
        if start > end {
            return String::new();
        }
        let wide = i128::from;
        let (per_line, width) = (wide(per_line), wide(width));
        let span = wide(end) - wide(start) + 1;
        let before = wide(start) - 1;
        let position = wide(offset) + before.div_euclid(per_line) * width + before % per_line;
        // The line ends of the span and of one line more on either side.
        let size = span + (span.div_euclid(per_line) + 2) * (width - per_line);
        // A seek takes an offset `i64` holds.
        let position = i64::try_from(position)
            .ok()
            .and_then(|at| u64::try_from(at).ok());
        let Some(position) = position else {
            crate::uncaught!("cannot fit 'int' into an offset-sized integer")
        };
        // The interpreter's `read` takes -1 for the whole of what is left and no other
        // negative size.
        let most = match size {
            -1 => u64::MAX,
            size if size < -1 => crate::uncaught!("read length must be non-negative or -1"),
            size => u64::try_from(size).unwrap_or(u64::MAX),
        };
        let mut raw = Vec::new();
        let mut handle = &self.handle;
        let read = handle
            .seek(io::SeekFrom::Start(position))
            .and_then(|_| handle.take(most).read_to_end(&mut raw));
        if let Err(error) = read {
            crate::uncaught!("{}", os_error(&error));
        }
        let mut bases: Vec<u8> = raw
            .into_iter()
            .filter(|byte| !matches!(byte, b'\r' | b'\n'))
            .take(usize::try_from(span).unwrap_or(usize::MAX))
            .collect();
        if let Some(at) = bases.iter().position(|byte| !byte.is_ascii()) {
            crate::uncaught!(
                "'ascii' codec can't decode byte {:#04x} in position {at}: ordinal not in range(128)",
                bases[at]
            );
        }
        bases.make_ascii_uppercase();
        String::from_utf8(bases).expect("ASCII is UTF-8")
    }
}
