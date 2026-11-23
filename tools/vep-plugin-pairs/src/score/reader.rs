// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The lines of a plain or gzip file and the sha256 of the file as stored, with a gzip
//! stream read as CPython's `gzip` reads one: every member, zero bytes after a member
//! skipped, and a stream it cannot read refused in that module's words at the read that
//! reaches the fault.

use super::paths::OsFailure;
use super::text::repr_bytes;
use super::Stop;
use crate::input::{DEFAULT_COLUMNS, TUPLE_COLUMNS};
use flate2::{Crc, Decompress, FlushDecompress, Status};
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{self, Read};

/// The stored bytes `gzip` takes at a time, and the bytes its buffered reader asks of it.
const UNIT: usize = 128 << 10;
/// The bytes of a plain file read at a time.
const PLAIN_UNIT: usize = 1 << 20;
const TRUNCATED: &str = "Compressed file ended before the end-of-stream marker was reached";

/// Why a file's bytes could not be read.
#[derive(Debug)]
pub enum Failure {
    /// The system refused the call.
    Os(OsFailure),
    /// The gzip stream is not one `gzip` reads; the reason in its words.
    Gzip(String),
}

impl From<io::Error> for Failure {
    fn from(error: io::Error) -> Failure {
        Failure::Os(OsFailure::of(&error, None))
    }
}

/// The file as stored, digested as it is read.
struct Stored {
    file: File,
    digest: Sha256,
}

impl Stored {
    /// Appends up to `count` bytes to `into`, fewer only at the file's end.
    fn read(&mut self, count: usize, into: &mut Vec<u8>) -> io::Result<usize> {
        let start = into.len();
        into.resize(start + count, 0);
        let mut got = 0;
        while got < count {
            match self.file.read(&mut into[start + got..]) {
                Ok(0) => break,
                Ok(read) => got += read,
                Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                Err(error) => {
                    into.truncate(start);
                    return Err(error);
                }
            }
        }
        into.truncate(start + got);
        self.digest.update(&into[start..]);
        Ok(got)
    }
}

/// The members of a gzip stream, inflated a read at a time.
struct Gunzip {
    inflater: Decompress,
    /// Stored bytes read and untaken: the inflater's input and, past a member's end,
    /// what follows the member.
    pending: Vec<u8>,
    at: usize,
    needs_input: bool,
    member_ended: bool,
    new_member: bool,
    crc: Crc,
    size: u64,
}

impl Gunzip {
    fn new() -> Gunzip {
        Gunzip {
            inflater: Decompress::new(false),
            pending: Vec::new(),
            at: 0,
            needs_input: true,
            member_ended: false,
            new_member: true,
            crc: Crc::new(),
            size: 0,
        }
    }

    /// Up to `count` stored bytes, fewer only at the file's end.
    fn take(&mut self, stored: &mut Stored, count: usize) -> io::Result<Vec<u8>> {
        let held = (self.pending.len() - self.at).min(count);
        let mut taken = self.pending[self.at..self.at + held].to_vec();
        self.at += held;
        if held < count {
            stored.read(count - held, &mut taken)?;
        }
        Ok(taken)
    }

    fn exactly(&mut self, stored: &mut Stored, count: usize) -> Result<Vec<u8>, Failure> {
        let taken = self.take(stored, count)?;
        if taken.len() < count {
            return Err(Failure::Gzip(TRUNCATED.to_string()));
        }
        Ok(taken)
    }

    /// Reads a member's header; `false` at the end of the file.
    fn header(&mut self, stored: &mut Stored) -> Result<bool, Failure> {
        let magic = self.take(stored, 2)?;
        if magic.is_empty() {
            return Ok(false);
        }
        if magic != [0x1f, 0x8b] {
            return Err(Failure::Gzip(format!(
                "Not a gzipped file ({})",
                repr_bytes(&magic)
            )));
        }
        let fixed = self.exactly(stored, 8)?;
        let (method, flags) = (fixed[0], fixed[1]);
        if method != 8 {
            return Err(Failure::Gzip("Unknown compression method".to_string()));
        }
        if flags & 4 != 0 {
            let length = self.exactly(stored, 2)?;
            self.exactly(
                stored,
                usize::from(u16::from_le_bytes([length[0], length[1]])),
            )?;
        }
        // The file name and the comment each end at a zero byte, or at the file's end.
        for flag in [8, 16] {
            if flags & flag != 0 {
                loop {
                    let byte = self.take(stored, 1)?;
                    if byte.is_empty() || byte[0] == 0 {
                        break;
                    }
                }
            }
        }
        if flags & 2 != 0 {
            self.exactly(stored, 2)?;
        }
        Ok(true)
    }

    /// Checks a member's trailer and skips the zero bytes after it.
    fn trailer(&mut self, stored: &mut Stored) -> Result<(), Failure> {
        let stated = self.exactly(stored, 8)?;
        let crc = u32::from_le_bytes([stated[0], stated[1], stated[2], stated[3]]);
        let size = u32::from_le_bytes([stated[4], stated[5], stated[6], stated[7]]);
        if crc != self.crc.sum() {
            return Err(Failure::Gzip(format!(
                "CRC check failed {crc:#x} != {:#x}",
                self.crc.sum()
            )));
        }
        if u64::from(size) != self.size & 0xffff_ffff {
            return Err(Failure::Gzip(
                "Incorrect length of data produced".to_string(),
            ));
        }
        loop {
            let byte = self.take(stored, 1)?;
            match byte.first() {
                Some(0) => {}
                Some(&other) => {
                    self.pending.drain(..self.at);
                    self.pending.insert(0, other);
                    self.at = 0;
                    break;
                }
                None => break,
            }
        }
        Ok(())
    }

    /// Appends up to `most` inflated bytes to `into` and returns their number: at least one
    /// unless the stream has ended.
    fn read(
        &mut self,
        stored: &mut Stored,
        most: usize,
        into: &mut Vec<u8>,
    ) -> Result<usize, Failure> {
        loop {
            if self.member_ended {
                self.trailer(stored)?;
                self.new_member = true;
                self.member_ended = false;
                self.inflater = Decompress::new(false);
                self.needs_input = true;
            }
            if self.new_member {
                self.crc = Crc::new();
                self.size = 0;
                if !self.header(stored)? {
                    return Ok(0);
                }
                self.new_member = false;
            }
            if self.needs_input {
                self.pending.drain(..self.at);
                self.at = 0;
                let held = self.pending.len();
                if held < UNIT {
                    stored.read(UNIT - held, &mut self.pending)?;
                }
            }
            let start = into.len();
            into.resize(start + most, 0);
            let (read, written) = (self.inflater.total_in(), self.inflater.total_out());
            let status = self
                .inflater
                .decompress(
                    &self.pending[self.at..],
                    &mut into[start..],
                    FlushDecompress::None,
                )
                .map_err(|_| {
                    Failure::Gzip(
                        "Error -3 while decompressing data: invalid deflate data".to_string(),
                    )
                })?;
            let consumed = (self.inflater.total_in() - read) as usize;
            let produced = (self.inflater.total_out() - written) as usize;
            self.at += consumed;
            into.truncate(start + produced);
            self.member_ended = status == Status::StreamEnd;
            self.needs_input = !self.member_ended && self.at == self.pending.len();
            if produced > 0 {
                self.crc.update(&into[start..]);
                self.size += produced as u64;
                return Ok(produced);
            }
            if !self.member_ended && consumed == 0 {
                // The inflater wants stored bytes the file does not hold.
                return Err(Failure::Gzip(TRUNCATED.to_string()));
            }
        }
    }
}

/// The lines of a plain or gzip file as bytes, and the sha256 of the file as stored.
pub struct Reader {
    stored: Stored,
    gunzip: Option<Gunzip>,
    /// Bytes read that wait to be handed out.
    block: Vec<u8>,
    at: usize,
}

impl Reader {
    /// Opens the file and tells a gzip stream by its first two bytes.
    pub fn open(path: &str) -> Result<Reader, OsFailure> {
        let file = File::open(path).map_err(|error| OsFailure::of(&error, Some(path)))?;
        let mut stored = Stored {
            file,
            digest: Sha256::new(),
        };
        let mut first = Vec::new();
        stored
            .read(PLAIN_UNIT, &mut first)
            .map_err(|error| OsFailure::of(&error, None))?;
        let mut reader = Reader {
            stored,
            gunzip: None,
            block: Vec::new(),
            at: 0,
        };
        if first.starts_with(&[0x1f, 0x8b]) {
            let mut gunzip = Gunzip::new();
            gunzip.pending = first;
            reader.gunzip = Some(gunzip);
        } else {
            reader.block = first;
        }
        Ok(reader)
    }

    /// Reads the next run of bytes; `false` at the end.
    fn fill(&mut self) -> Result<bool, Failure> {
        self.block.clear();
        self.at = 0;
        let read = match &mut self.gunzip {
            Some(gunzip) => gunzip.read(&mut self.stored, UNIT, &mut self.block)?,
            None => self.stored.read(PLAIN_UNIT, &mut self.block)?,
        };
        Ok(read > 0)
    }

    /// Puts the next line into `line`, its line feed included where the file writes one;
    /// `false` after the last line.
    pub fn line(&mut self, line: &mut Vec<u8>) -> Result<bool, Failure> {
        line.clear();
        loop {
            if self.at == self.block.len() && !self.fill()? {
                return Ok(!line.is_empty());
            }
            let rest = &self.block[self.at..];
            match memchr::memchr(b'\n', rest) {
                Some(end) => {
                    line.extend_from_slice(&rest[..=end]);
                    self.at += end + 1;
                    return Ok(true);
                }
                None => {
                    line.extend_from_slice(rest);
                    self.at = self.block.len();
                }
            }
        }
    }

    /// The file's first `most` bytes, fewer where it holds fewer.
    pub fn head(&mut self, most: usize) -> Result<Vec<u8>, Failure> {
        let mut head = Vec::new();
        while head.len() < most {
            if self.at == self.block.len() && !self.fill()? {
                break;
            }
            let taken = (self.block.len() - self.at).min(most - head.len());
            head.extend_from_slice(&self.block[self.at..self.at + taken]);
            self.at += taken;
        }
        Ok(head)
    }

    /// The sha256 of the file as stored, the unread bytes included.
    pub fn sha256(mut self) -> io::Result<String> {
        let mut rest = Vec::new();
        loop {
            rest.clear();
            if self.stored.read(PLAIN_UNIT, &mut rest)? == 0 {
                break;
            }
        }
        Ok(crate::order::hex(&self.stored.digest.finalize()))
    }
}

/// A line without its line end: every trailing carriage return and line feed.
pub fn without_line_end(mut line: &[u8]) -> &[u8] {
    while let [kept @ .., b'\r' | b'\n'] = line {
        line = kept;
    }
    line
}

/// Refuses an input that cannot be opened, whose column header names other than 14 columns,
/// or whose first data line is neither a default-format row of 14 columns nor, where one is
/// allowed and no `#` line stands before it, a tuple file line of five.
pub fn check_first_line(path: &str, tuple_file_allowed: bool) -> Result<(), Stop> {
    let mut reader = Reader::open(path).map_err(Stop::Os)?;
    let mut headed = false;
    let mut line = Vec::new();
    let mut number = 0u64;
    loop {
        match reader.line(&mut line) {
            Ok(true) => {}
            Ok(false) => return Ok(()),
            Err(Failure::Os(failure)) => return Err(Stop::Os(failure)),
            Err(Failure::Gzip(reason)) => {
                return Err(Stop::Refused(format!(
                    "{path}: cannot read the gzip stream: {reason}"
                )))
            }
        }
        number += 1;
        let row = without_line_end(&line);
        let columns = memchr::memchr_iter(b'\t', row).count() + 1;
        if row.starts_with(b"#") {
            headed = true;
            if !row.starts_with(b"##") && columns != DEFAULT_COLUMNS {
                return Err(Stop::Refused(format!(
                    "{path}: line {number}, the column header, names {columns} columns, not 14"
                )));
            }
        } else if !row.is_empty() {
            let tuples = columns == TUPLE_COLUMNS && tuple_file_allowed && !headed;
            if !tuples && columns != DEFAULT_COLUMNS {
                let allowed = if tuple_file_allowed {
                    ", a tuple file line 5"
                } else {
                    ""
                };
                return Err(Stop::Refused(format!(
                    "{path}: line {number} has {columns} columns; a default-format row has 14{allowed}"
                )));
            }
            return Ok(());
        }
    }
}

/// The sha256 of a file as stored.
pub fn digest(path: &str) -> Result<String, OsFailure> {
    let named = |error: io::Error| OsFailure::of(&error, Some(path));
    let mut file = File::open(path).map_err(named)?;
    let mut digest = Sha256::new();
    let mut block = vec![0u8; PLAIN_UNIT];
    loop {
        match file.read(&mut block) {
            Ok(0) => break,
            Ok(read) => digest.update(&block[..read]),
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(OsFailure::of(&error, None)),
        }
    }
    Ok(crate::order::hex(&digest.finalize()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;
    use std::path::PathBuf;

    fn gzip(content: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(content).unwrap();
        encoder.finish().unwrap()
    }

    struct Scratch(PathBuf);

    impl Scratch {
        fn new(tag: &str) -> Scratch {
            let path = std::env::temp_dir().join(format!(
                "vep-plugin-pairs-reader-{}-{tag}",
                std::process::id()
            ));
            let _ = std::fs::remove_dir_all(&path);
            std::fs::create_dir_all(&path).unwrap();
            Scratch(path)
        }

        fn file(&self, name: &str, content: &[u8]) -> String {
            let path = self.0.join(name);
            std::fs::write(&path, content).unwrap();
            path.to_str().unwrap().to_string()
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// The lines of the file, or the gzip refusal the reading ends in.
    fn lines(path: &str) -> Result<Vec<Vec<u8>>, String> {
        let mut reader = Reader::open(path).unwrap();
        let mut found = Vec::new();
        let mut line = Vec::new();
        loop {
            match reader.line(&mut line) {
                Ok(true) => found.push(line.clone()),
                Ok(false) => return Ok(found),
                Err(Failure::Gzip(reason)) => return Err(reason),
                Err(Failure::Os(failure)) => panic!("{}", failure.text()),
            }
        }
    }

    fn sha256(content: &[u8]) -> String {
        crate::order::hex(&Sha256::digest(content))
    }

    #[test]
    fn plain_and_gzip_files_give_the_same_lines_and_their_own_digest() {
        let scratch = Scratch::new("lines");
        let content = b"one\ntwo\r\n\nlast without an end";
        let plain = scratch.file("plain.txt", content);
        let packed_bytes = gzip(content);
        let packed = scratch.file("packed.txt.gz", &packed_bytes);
        let expected: Vec<Vec<u8>> = vec![
            b"one\n".to_vec(),
            b"two\r\n".to_vec(),
            b"\n".to_vec(),
            b"last without an end".to_vec(),
        ];
        assert_eq!(lines(&plain).unwrap(), expected);
        assert_eq!(lines(&packed).unwrap(), expected);
        assert_eq!(
            Reader::open(&plain).unwrap().sha256().unwrap(),
            sha256(content)
        );
        assert_eq!(
            Reader::open(&packed).unwrap().sha256().unwrap(),
            sha256(&packed_bytes)
        );
        assert_eq!(digest(&packed).unwrap(), sha256(&packed_bytes));
        assert_eq!(lines(&scratch.file("empty.txt", b"")).unwrap().len(), 0);
        assert_eq!(without_line_end(b"a\r\r\n\n"), b"a");
    }

    #[test]
    fn every_member_is_read_and_zero_padding_is_skipped() {
        let scratch = Scratch::new("members");
        let mut stream = gzip(b"first\n");
        stream.extend_from_slice(&[0, 0, 0]);
        stream.extend_from_slice(&gzip(b"second\n"));
        stream.extend_from_slice(&[0, 0]);
        let path = scratch.file("two.gz", &stream);
        assert_eq!(
            lines(&path).unwrap(),
            vec![b"first\n".to_vec(), b"second\n".to_vec()]
        );
    }

    #[test]
    fn long_stream_is_read_across_many_reads() {
        let scratch = Scratch::new("long");
        let content: Vec<u8> = (0..400_000u32)
            .flat_map(|number| format!("line {number}\n").into_bytes())
            .collect();
        let path = scratch.file("long.gz", &gzip(&content));
        let found = lines(&path).unwrap();
        assert_eq!(found.len(), 400_000);
        assert_eq!(found.concat(), content);
        let mut reader = Reader::open(&path).unwrap();
        assert_eq!(reader.head(1 << 16).unwrap(), &content[..1 << 16]);
    }

    #[test]
    fn stream_that_cannot_be_read_is_refused_in_the_words_of_gzip() {
        let scratch = Scratch::new("broken");
        let content: Vec<u8> = b"hello world\n".repeat(1000);
        let whole = gzip(&content);
        let cut = |end: usize| lines(&scratch.file("cut.gz", &whole[..end])).unwrap_err();
        for end in [2, 5, 20, whole.len() - 12, whole.len() - 3] {
            assert_eq!(cut(end), TRUNCATED, "{end}");
        }
        let with = |tail: &[u8]| {
            let mut stream = whole.clone();
            stream.extend_from_slice(tail);
            lines(&scratch.file("tail.gz", &stream))
        };
        assert_eq!(with(b"xy").unwrap_err(), "Not a gzipped file (b'xy')");
        assert_eq!(with(b"\0\0x").unwrap_err(), "Not a gzipped file (b'x')");
        assert_eq!(with(b"\x1f").unwrap_err(), "Not a gzipped file (b'\\x1f')");
        assert_eq!(with(b"\0\0\0").unwrap().len(), 1000);
        let mut wrong = whole.clone();
        let end = wrong.len();
        let length = wrong[end - 4..].to_vec();
        wrong[end - 8..end - 4].copy_from_slice(&[1, 2, 3, 4]);
        let refused = lines(&scratch.file("crc.gz", &wrong)).unwrap_err();
        assert!(
            refused.starts_with("CRC check failed 0x4030201 != 0x"),
            "{refused}"
        );
        let mut wrong = whole.clone();
        wrong[end - 4..].copy_from_slice(&[1, 2, 3, 4]);
        assert_eq!(
            lines(&scratch.file("length.gz", &wrong)).unwrap_err(),
            "Incorrect length of data produced"
        );
        assert_ne!(length, [1, 2, 3, 4]);
        let mut wrong = whole.clone();
        wrong[2] = 7;
        assert_eq!(
            lines(&scratch.file("method.gz", &wrong)).unwrap_err(),
            "Unknown compression method"
        );
    }

    #[test]
    fn lines_before_a_fault_are_read_until_the_read_that_reaches_it() {
        let scratch = Scratch::new("early");
        // The first lines of a stream cut far behind them are read; the fault is met later.
        let content: Vec<u8> = (0..200_000u32)
            .flat_map(|number| format!("row {number}\t{}\n", number * 7919).into_bytes())
            .collect();
        let whole = gzip(&content);
        let path = scratch.file("late.gz", &whole[..whole.len() - 100]);
        let mut reader = Reader::open(&path).unwrap();
        let mut line = Vec::new();
        assert!(reader.line(&mut line).unwrap());
        assert_eq!(line, b"row 0\t0\n");
        assert_eq!(lines(&path).unwrap_err(), TRUNCATED);
    }
}
