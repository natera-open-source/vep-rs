// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A file read as CPython reads it with `open` or `gzip.open`: the bytes of a plain file or
//! of every member of a gzip file, and the lines of that text under `encoding="utf-8"`,
//! `errors="replace"` and universal newlines. A file that cannot be opened or read gives the
//! reason in the words of the CPython exception.

mod bgzf;

use crate::error::{failure, Failure};
use crate::pytext::os_error_text;
use flate2::{Crc, Decompress, DecompressError, FlushDecompress, Status};
use std::fmt::Write;
use std::fs::File;
use std::io::Read;

/// The bytes read from a file at a time.
const CHUNK: usize = 128 * 1024;

/// `str(exc)` of the `EOFError` a gzip file raises when it ends inside a member.
pub(crate) const CUT_SHORT: &str =
    "Compressed file ended before the end-of-stream marker was reached";

/// The tool's own words for deflate data that does not inflate.
const CORRUPT: &str = "corrupt deflate data in a gzip member";

/// Why deflate data does not inflate: `str(exc)` of CPython's `zlib.error` where the
/// deflate library linked names its reason, as zlib does, and the tool's own words where it
/// names none.
pub(crate) fn corrupt(error: &DecompressError) -> ReadError {
    match error.message() {
        Some(message) => {
            ReadError::Cpython(format!("Error -3 while decompressing data: {message}"))
        }
        None => ReadError::Tool(CORRUPT.to_string()),
    }
}

/// Why a file could not be opened or read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadError {
    /// `str(exc)` of the `OSError`, `EOFError` or `zlib.error` CPython raises for it.
    Cpython(String),
    /// The tool's own reason, where CPython raises a `zlib.error` in words of zlib's that
    /// the deflate library linked does not give.
    Tool(String),
}

impl ReadError {
    /// The reason, without the sentence around it.
    pub fn reason(&self) -> &str {
        match self {
            ReadError::Cpython(reason) | ReadError::Tool(reason) => reason,
        }
    }

    /// The failure that ends a run on this error, `message` making the line's message from
    /// the reason: the program's own line for a reason in CPython's words, the tool's own
    /// line for a reason of the tool's.
    pub fn failure(&self, message: impl FnOnce(&str) -> String) -> Failure {
        match self {
            ReadError::Cpython(reason) => failure(message(reason)),
            ReadError::Tool(reason) => Failure::own(message(reason)),
        }
    }
}

/// The refusal of a reader of text lines: `<path> is not a readable text file (<reason>)`.
pub fn not_readable_text(path: &str, error: &ReadError) -> Failure {
    error.failure(|reason| format!("{path} is not a readable text file ({reason})"))
}

/// `repr(bytes)`: `b'...'`, in double quotes when the bytes hold `'` and no `"`; `\t`, `\n`,
/// `\r`, the backslash and the quote escaped; every byte outside `0x20..=0x7e` as `\x` and
/// two lowercase hex digits.
pub fn py_repr_bytes(bytes: &[u8]) -> String {
    let quote = if bytes.contains(&b'\'') && !bytes.contains(&b'"') {
        b'"'
    } else {
        b'\''
    };
    let mut out = String::from("b");
    out.push(char::from(quote));
    for byte in bytes {
        match byte {
            b'\\' => out.push_str("\\\\"),
            b'\t' => out.push_str("\\t"),
            b'\n' => out.push_str("\\n"),
            b'\r' => out.push_str("\\r"),
            _ if *byte == quote => {
                out.push('\\');
                out.push(char::from(quote));
            }
            b' '..=b'~' => out.push(char::from(*byte)),
            _ => write!(out, "\\x{byte:02x}").unwrap(),
        }
    }
    out.push(char::from(quote));
    out
}

/// `str(exc)` of the `OSError` a read of an open file raises: `[Errno <n>] <strerror>`.
fn read_error(error: &std::io::Error) -> ReadError {
    let text = match error.raw_os_error() {
        Some(code) => {
            let described = std::io::Error::from_raw_os_error(code).to_string();
            let strerror = described
                .strip_suffix(&format!(" (os error {code})"))
                .unwrap_or(&described);
            format!("[Errno {code}] {strerror}")
        }
        None => error.to_string(),
    };
    ReadError::Cpython(text)
}

/// One read of a file, repeated when a signal interrupts it.
fn read_some(file: &mut File, into: &mut [u8]) -> Result<usize, ReadError> {
    loop {
        match file.read(into) {
            Ok(count) => return Ok(count),
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(read_error(&error)),
        }
    }
}

/// The compressed bytes of a gzip file, read ahead a chunk at a time.
struct Compressed {
    file: File,
    bytes: Vec<u8>,
    start: usize,
    end: usize,
}

impl Compressed {
    /// The bytes read ahead, after one more read when there are none; empty at the end of
    /// the file.
    fn available(&mut self) -> Result<&[u8], ReadError> {
        if self.start == self.end {
            self.start = 0;
            self.end = read_some(&mut self.file, &mut self.bytes)?;
        }
        Ok(&self.bytes[self.start..self.end])
    }

    fn next_byte(&mut self) -> Result<Option<u8>, ReadError> {
        let byte = self.available()?.first().copied();
        self.start += usize::from(byte.is_some());
        Ok(byte)
    }

    /// The next `N` bytes; a file that ends before them is cut short.
    fn exact<const N: usize>(&mut self) -> Result<[u8; N], ReadError> {
        let mut bytes = [0u8; N];
        for byte in &mut bytes {
            *byte = self
                .next_byte()?
                .ok_or_else(|| ReadError::Cpython(CUT_SHORT.to_string()))?;
        }
        Ok(bytes)
    }

    fn skip_exact(&mut self, count: usize) -> Result<(), ReadError> {
        for _ in 0..count {
            self.exact::<1>()?;
        }
        Ok(())
    }

    /// Skips a zero-terminated field of a member header; the end of the file ends it too.
    fn skip_text(&mut self) -> Result<(), ReadError> {
        while self.next_byte()?.is_some_and(|byte| byte != 0) {}
        Ok(())
    }
}

/// Where the reading of a gzip file stands.
enum Member {
    /// Before a member's header, or at the end of the file.
    Header,
    /// Inside a member's deflate data.
    Body,
    /// Before a member's checksum and length.
    Trailer,
    /// After the last member.
    End,
}

/// A gzip file as `gzip.open(path, "rb")` reads it: every member in turn, the zero bytes
/// between members skipped, each member held to its checksum and its length.
struct Gzip {
    compressed: Compressed,
    member: Member,
    inflate: Decompress,
    crc: Crc,
}

impl Gzip {
    /// Reads a member header; `false` at the end of the file.
    fn header(&mut self) -> Result<bool, ReadError> {
        let Some(first) = self.compressed.next_byte()? else {
            return Ok(false);
        };
        let mut magic = vec![first];
        magic.extend(self.compressed.next_byte()?);
        if magic != [0x1f, 0x8b] {
            return Err(ReadError::Cpython(format!(
                "Not a gzipped file ({})",
                py_repr_bytes(&magic)
            )));
        }
        let [method, flags, ..] = self.compressed.exact::<8>()?;
        if method != 8 {
            return Err(ReadError::Cpython("Unknown compression method".to_string()));
        }
        if flags & 4 != 0 {
            let length = u16::from_le_bytes(self.compressed.exact::<2>()?);
            self.compressed.skip_exact(usize::from(length))?;
        }
        if flags & 8 != 0 {
            self.compressed.skip_text()?;
        }
        if flags & 16 != 0 {
            self.compressed.skip_text()?;
        }
        if flags & 2 != 0 {
            self.compressed.exact::<2>()?;
        }
        Ok(true)
    }

    /// Holds the member that ended to its checksum and length, and skips the zero bytes
    /// after it.
    fn trailer(&mut self) -> Result<(), ReadError> {
        let [c0, c1, c2, c3, s0, s1, s2, s3] = self.compressed.exact::<8>()?;
        let stored = u32::from_le_bytes([c0, c1, c2, c3]);
        if stored != self.crc.sum() {
            return Err(ReadError::Cpython(format!(
                "CRC check failed {stored:#x} != {:#x}",
                self.crc.sum()
            )));
        }
        if u32::from_le_bytes([s0, s1, s2, s3]) != self.crc.amount() {
            return Err(ReadError::Cpython(
                "Incorrect length of data produced".to_string(),
            ));
        }
        self.skip_padding()
    }

    /// Skips the zero bytes between a member and the next.
    fn skip_padding(&mut self) -> Result<(), ReadError> {
        while self.compressed.available()?.first() == Some(&0) {
            self.compressed.start += 1;
        }
        Ok(())
    }

    /// Inflates into `into`; 0 after the last member.
    fn read(&mut self, into: &mut [u8]) -> Result<usize, ReadError> {
        loop {
            match self.member {
                Member::End => return Ok(0),
                Member::Header => {
                    if !self.header()? {
                        self.member = Member::End;
                        return Ok(0);
                    }
                    self.inflate.reset(false);
                    self.crc.reset();
                    self.member = Member::Body;
                }
                Member::Body => {
                    let input = self.compressed.available()?;
                    let file_ended = input.is_empty();
                    let (read_before, written_before) =
                        (self.inflate.total_in(), self.inflate.total_out());
                    let status = self
                        .inflate
                        .decompress(input, into, FlushDecompress::None)
                        .map_err(|error| corrupt(&error))?;
                    let consumed = (self.inflate.total_in() - read_before) as usize;
                    let produced = (self.inflate.total_out() - written_before) as usize;
                    self.compressed.start += consumed;
                    self.crc.update(&into[..produced]);
                    if status == Status::StreamEnd {
                        self.member = Member::Trailer;
                    } else if produced == 0 && file_ended {
                        return Err(ReadError::Cpython(CUT_SHORT.to_string()));
                    } else if produced == 0 && consumed == 0 {
                        return Err(ReadError::Tool(CORRUPT.to_string()));
                    }
                    if produced > 0 {
                        return Ok(produced);
                    }
                }
                Member::Trailer => {
                    self.trailer()?;
                    self.member = Member::Header;
                }
            }
        }
    }
}

/// The bytes of a file: a plain file's own, a gzip file's members inflated.
pub struct ByteSource {
    file: Opened,
}

enum Opened {
    Plain(File),
    Gzip(Box<Gzip>),
    /// A gzip file read as BGZF blocks until a member that is none shows.
    Blocks(Box<bgzf::Bgzf>),
}

impl ByteSource {
    /// `opener(path)(path, "rb")`: `gzip.open` for a name ending in `.gz` or `.bgz` (BGZF is
    /// gzip), else `open`.
    pub fn open(path: &str) -> Result<ByteSource, ReadError> {
        if path.ends_with(".gz") || path.ends_with(".bgz") {
            ByteSource::open_gzip(path)
        } else {
            ByteSource::open_plain(path)
        }
    }

    /// `open(path, "rb")`.
    pub fn open_plain(path: &str) -> Result<ByteSource, ReadError> {
        Ok(ByteSource {
            file: Opened::Plain(open_file(path)?),
        })
    }

    /// `gzip.open(path, "rb")`: the file is opened; whether it is gzip shows at the first
    /// read. The blocks of a BGZF file are inflated on the threads of the run's bound.
    pub fn open_gzip(path: &str) -> Result<ByteSource, ReadError> {
        ByteSource::open_gzip_on(path, super::blocks::thread_count().unwrap_or(1))
    }

    /// `open_gzip` with the blocks of a BGZF file inflated on at most `threads` threads;
    /// the bytes are the same on every count.
    pub fn open_gzip_on(path: &str, threads: usize) -> Result<ByteSource, ReadError> {
        if threads > 1 {
            if let Some(blocks) = bgzf::Bgzf::open(path, threads)? {
                return Ok(ByteSource {
                    file: Opened::Blocks(Box::new(blocks)),
                });
            }
        }
        let gzip = Gzip {
            compressed: Compressed {
                file: open_file(path)?,
                bytes: vec![0; CHUNK],
                start: 0,
                end: 0,
            },
            member: Member::Header,
            inflate: Decompress::new(false),
            crc: Crc::new(),
        };
        Ok(ByteSource {
            file: Opened::Gzip(Box::new(gzip)),
        })
    }

    /// Reads into `into`; 0 at the end, and for an `into` without room.
    pub fn read(&mut self, into: &mut [u8]) -> Result<usize, ReadError> {
        loop {
            let (offset, after_member) = match &mut self.file {
                _ if into.is_empty() => return Ok(0),
                Opened::Plain(file) => return read_some(file, into),
                Opened::Gzip(gzip) => return gzip.read(into),
                Opened::Blocks(blocks) => match blocks.read(into) {
                    bgzf::Step::Bytes(count) => return Ok(count),
                    bgzf::Step::Members {
                        offset,
                        after_member,
                    } => (offset, after_member),
                },
            };
            // From here on the file is read member by member.
            if let Opened::Blocks(blocks) = &mut self.file {
                let gzip = blocks.members(offset, after_member)?;
                self.file = Opened::Gzip(Box::new(gzip));
            }
        }
    }

    /// Every byte from here to the end.
    pub fn read_to_end(mut self) -> Result<Vec<u8>, ReadError> {
        let mut bytes = Vec::new();
        let mut chunk = vec![0u8; CHUNK];
        loop {
            let count = self.read(&mut chunk)?;
            if count == 0 {
                return Ok(bytes);
            }
            bytes.extend_from_slice(&chunk[..count]);
        }
    }
}

/// Opens a file for reading; a directory is refused as `open` refuses it.
fn open_file(path: &str) -> Result<File, ReadError> {
    let file = File::open(path).map_err(|error| ReadError::Cpython(os_error_text(&error, path)))?;
    let is_directory = file
        .metadata()
        .map_err(|error| ReadError::Cpython(os_error_text(&error, path)))?
        .is_dir();
    if is_directory {
        let error = std::io::Error::from_raw_os_error(21);
        return Err(ReadError::Cpython(os_error_text(&error, path)));
    }
    Ok(file)
}

/// The lines of a text file as iterating `open(path, "rt", encoding="utf-8",
/// errors="replace")` yields them, each without its line end: a line ends at `\n`, `\r` or
/// `\r\n`, a last line needs no line end, and each run of bytes that is not UTF-8 reads as
/// U+FFFD.
pub struct TextLines {
    source: ByteSource,
    bytes: Vec<u8>,
    /// The first byte no line has taken.
    start: usize,
    /// The first byte not searched for a line end.
    scan: usize,
    /// One past the last byte read.
    end: usize,
    source_ended: bool,
    line: String,
}

impl TextLines {
    /// `opener(path)(path, "rt", encoding="utf-8", errors="replace")`.
    pub fn open(path: &str) -> Result<TextLines, ReadError> {
        ByteSource::open(path).map(TextLines::over)
    }

    /// `gzip.open(path, "rt", encoding="utf-8", errors="replace")`, whatever the name.
    pub fn open_gzip(path: &str) -> Result<TextLines, ReadError> {
        ByteSource::open_gzip(path).map(TextLines::over)
    }

    /// The lines of `source`.
    pub fn over(source: ByteSource) -> TextLines {
        TextLines {
            source,
            bytes: vec![0; 2 * CHUNK],
            start: 0,
            scan: 0,
            end: 0,
            source_ended: false,
            line: String::new(),
        }
    }

    /// Moves to the next line; `false` after the last. A file that cannot be read on gives
    /// the error in place of the line it would have ended, as the lines before it were
    /// given.
    pub fn advance(&mut self) -> Result<bool, ReadError> {
        loop {
            let found = memchr::memchr2(b'\n', b'\r', &self.bytes[self.scan..self.end]);
            let Some(at) = found.map(|offset| self.scan + offset) else {
                self.scan = self.end;
                if !self.source_ended {
                    self.fill()?;
                    continue;
                }
                if self.start == self.end {
                    return Ok(false);
                }
                self.keep(self.end, self.end);
                return Ok(true);
            };
            let mut next = at + 1;
            if self.bytes[at] == b'\r' {
                if next == self.end && !self.source_ended {
                    // Whether `\n` follows shows only with the next bytes.
                    self.scan = at;
                    self.fill()?;
                    continue;
                }
                if next < self.end && self.bytes[next] == b'\n' {
                    next += 1;
                }
            }
            self.keep(at, next);
            return Ok(true);
        }
    }

    /// The line `advance` moved to.
    pub fn line(&self) -> &str {
        &self.line
    }

    /// The next line, `None` after the last.
    pub fn next_line(&mut self) -> Result<Option<&str>, ReadError> {
        Ok(if self.advance()? {
            Some(self.line.as_str())
        } else {
            None
        })
    }

    /// Takes the bytes up to `line_end` as the line and moves past `next`.
    fn keep(&mut self, line_end: usize, next: usize) {
        let bytes = &self.bytes[self.start..line_end];
        self.line.clear();
        match std::str::from_utf8(bytes) {
            Ok(text) => self.line.push_str(text),
            Err(_) => self.line.push_str(&String::from_utf8_lossy(bytes)),
        }
        self.start = next;
        self.scan = next;
    }

    /// Reads more bytes behind the ones no line has taken.
    fn fill(&mut self) -> Result<(), ReadError> {
        if self.start > 0 {
            self.bytes.copy_within(self.start..self.end, 0);
            self.scan -= self.start;
            self.end -= self.start;
            self.start = 0;
        }
        if self.end == self.bytes.len() {
            self.bytes.resize(self.bytes.len() * 2, 0);
        }
        let count = self.source.read(&mut self.bytes[self.end..])?;
        self.source_ended = count == 0;
        self.end += count;
        Ok(())
    }
}

#[cfg(test)]
#[path = "text_tests.rs"]
mod text_tests;
