// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Python's reading of a file and its spelling of names: universal newlines, strict UTF-8, gzip members, the words of its failures, `repr()` of a text, and the name, stem and suffix of a path.

use std::ffi::{OsStr, OsString};
use std::fmt;
use std::fs::File;
use std::io::{self, Read};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};
use std::sync::LazyLock;

use flate2::{Crc, Decompress, FlushDecompress, Status};
use memchr::{memchr2, memchr2_iter, memrchr};
use regex::Regex;

use crate::failure::Failure;

/// About how many bytes of a file a reader hands to a worker at a time.
pub const BLOCK_BYTES: usize = 1 << 20;

/// The longest decimal text Python's `int()` converts; a longer one raises `ValueError`.
pub const INT_MAX_STR_DIGITS: usize = 4300;

/// The failure `int()` ends in on a decimal text of `digits` digits, more than
/// [`INT_MAX_STR_DIGITS`].
pub fn int_conversion_failure(digits: usize) -> Failure {
    Failure::value(format!(
        "Exceeds the limit ({INT_MAX_STR_DIGITS} digits) for integer string conversion: value has {digits} digits; use sys.set_int_max_str_digits() to increase the limit"
    ))
}

/// Why a line of an input file ends the comparison.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LineDefect {
    /// The line is not UTF-8.
    NotUtf8,
    /// A decimal text of this many digits, more than [`INT_MAX_STR_DIGITS`].
    IntegerDigits(usize),
    /// A row of an engine's output with this many columns, fewer than seven.
    ShortRow(usize),
}

/// The first line of a run of lines that ends the comparison: the line's number within
/// the run, from 1, and why.
pub type RunDefect = (u64, LineDefect);

/// The failure a defective line of `path` ends the comparison with, `line` counted from 1
/// over the whole file.
pub fn line_failure(path: &Path, line: u64, defect: LineDefect) -> Failure {
    match defect {
        LineDefect::NotUtf8 => Failure::decode(format!(
            "{}: line {line} is not UTF-8",
            path.display()
        )),
        LineDefect::IntegerDigits(digits) => int_conversion_failure(digits),
        LineDefect::ShortRow(columns) => Failure::value(format!(
            "{}: line {line} has {columns} columns, expected at least 7. A truncated or malformed engine output must fail the comparison rather than silently reduce the tuple set.",
            path.display()
        )),
    }
}

/// The first defective line of a file read as runs of lines, given each run's line count
/// and first defect in the order of the file: its number over the whole file and why.
pub fn first_defect(runs: impl IntoIterator<Item = (u64, Option<RunDefect>)>) -> Option<RunDefect> {
    let mut before = 0;
    for (lines, defect) in runs {
        if let Some((line, defect)) = defect {
            return Some((before + line, defect));
        }
        before += lines;
    }
    None
}

/// The length of the longest prefix of `buf` that ends on a line boundary: through its last
/// `\n`, or without one through its last `\r` that another byte follows.
fn line_cut(buf: &[u8]) -> Option<usize> {
    if let Some(i) = memrchr(b'\n', buf) {
        return Some(i + 1);
    }
    memrchr(b'\r', &buf[..buf.len().saturating_sub(1)]).map(|i| i + 1)
}

/// Reads `reader` to its end and calls `emit` with runs of whole lines, in order, of about
/// `block_bytes` each; only the last run may end without a line terminator.
pub fn read_blocks(
    mut reader: impl Read,
    block_bytes: usize,
    mut emit: impl FnMut(Vec<u8>),
) -> io::Result<()> {
    let block_bytes = block_bytes.max(1);
    let mut carry: Vec<u8> = Vec::new();
    loop {
        let mut buf = Vec::with_capacity(carry.len() + block_bytes);
        buf.append(&mut carry);
        let held = buf.len();
        buf.resize(held + block_bytes, 0);
        let mut filled = held;
        let mut eof = false;
        while filled < buf.len() {
            match reader.read(&mut buf[filled..]) {
                Ok(0) => {
                    eof = true;
                    break;
                }
                Ok(n) => filled += n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        buf.truncate(filled);
        if eof {
            if !buf.is_empty() {
                emit(buf);
            }
            return Ok(());
        }
        match line_cut(&buf) {
            Some(cut) => {
                carry.extend_from_slice(&buf[cut..]);
                buf.truncate(cut);
                emit(buf);
            }
            None => carry = buf,
        }
    }
}

/// The lines of a run of whole lines, each without its terminator: `\n`, `\r\n` and a lone
/// `\r` each end a line, as Python's text mode reads them, and a last line needs none.
pub fn lines(block: &[u8]) -> Lines<'_> {
    Lines { rest: block }
}

/// The iterator of [`lines`].
pub struct Lines<'a> {
    rest: &'a [u8],
}

impl<'a> Iterator for Lines<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<&'a [u8]> {
        if self.rest.is_empty() {
            return None;
        }
        let Some(i) = memchr2(b'\n', b'\r', self.rest) else {
            return Some(std::mem::take(&mut self.rest));
        };
        let line = &self.rest[..i];
        let mut next = i + 1;
        if self.rest[i] == b'\r' && self.rest.get(next) == Some(&b'\n') {
            next += 1;
        }
        self.rest = &self.rest[next..];
        Some(line)
    }
}

/// How many lines [`lines`] yields for a run of whole lines.
pub fn line_count(block: &[u8]) -> u64 {
    let mut count = 0;
    for i in memchr2_iter(b'\n', b'\r', block) {
        if block[i] == b'\n' || block.get(i + 1) != Some(&b'\n') {
            count += 1;
        }
    }
    if block.last().is_some_and(|b| *b != b'\n' && *b != b'\r') {
        count += 1;
    }
    count
}

/// A reader that ends, silently, where the reader inside it fails.
pub struct UntilFailure<R> {
    inner: R,
    ended: bool,
}

impl<R> UntilFailure<R> {
    pub fn new(inner: R) -> UntilFailure<R> {
        UntilFailure {
            inner,
            ended: false,
        }
    }
}

impl<R: Read> Read for UntilFailure<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while !self.ended {
            match self.inner.read(buf) {
                Ok(n) => return Ok(n),
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(_) => self.ended = true,
            }
        }
        Ok(0)
    }
}

/// The bytes Python's text reader asks its file for at a time.
const PY_CHUNK_BYTES: usize = 8192;

/// The system error number of a directory opened as a file.
const EISDIR: i32 = 21;

/// Opens a file for reading as Python's `open` does: a directory is refused where it is
/// opened, with the error a read of it would give.
pub fn open_file(path: &Path) -> io::Result<File> {
    let file = File::open(path)?;
    if file.metadata()?.is_dir() {
        return Err(io::Error::from_raw_os_error(EISDIR));
    }
    Ok(file)
}

/// `repr()` of the text Python holds for the bytes of a path: in single quotes, or in
/// double quotes when it holds a single quote and no double one; a backslash, the quote
/// in use, a tab, a line feed and a carriage return escaped; every other character that
/// is not printable as `\xNN`, `\uNNNN` or `\UNNNNNNNN`; a byte that is not UTF-8 as the
/// lone surrogate Python carries it as (`\udcNN`).
pub fn py_repr(bytes: &[u8]) -> String {
    static NOT_PRINTABLE: LazyLock<Regex> = LazyLock::new(|| {
        Regex::new(r"^[\p{Cc}\p{Cf}\p{Co}\p{Cn}\p{Zl}\p{Zp}\p{Zs}]$").expect("regex")
    });
    let single = bytes.contains(&b'\'');
    let quote = if single && !bytes.contains(&b'"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(bytes.len() + 2);
    out.push(quote);
    for chunk in bytes.utf8_chunks() {
        for c in chunk.valid().chars() {
            let mut buffer = [0u8; 4];
            match c {
                '\\' => out.push_str("\\\\"),
                '\t' => out.push_str("\\t"),
                '\n' => out.push_str("\\n"),
                '\r' => out.push_str("\\r"),
                c if c == quote => {
                    out.push('\\');
                    out.push(c);
                }
                ' '..='~' => out.push(c),
                c if c.is_ascii() || NOT_PRINTABLE.is_match(c.encode_utf8(&mut buffer)) => {
                    match c as u32 {
                        n @ 0..=0xff => out.push_str(&format!("\\x{n:02x}")),
                        n @ 0x100..=0xffff => out.push_str(&format!("\\u{n:04x}")),
                        n => out.push_str(&format!("\\U{n:08x}")),
                    }
                }
                c => out.push(c),
            }
        }
        for byte in chunk.invalid() {
            out.push_str(&format!("\\udc{byte:02x}"));
        }
    }
    out.push(quote);
    out
}

/// The text of an exception Python's own readers raise, carried inside an [`io::Error`].
#[derive(Debug)]
pub struct PyError(pub String);

impl fmt::Display for PyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for PyError {}

fn py_error(text: impl Into<String>) -> io::Error {
    io::Error::other(PyError(text.into()))
}

/// `repr()` of a Python `bytes` value.
fn py_bytes_repr(bytes: &[u8]) -> String {
    let quote = if bytes.contains(&b'\'') && !bytes.contains(&b'"') {
        b'"'
    } else {
        b'\''
    };
    let mut out = String::from("b");
    out.push(char::from(quote));
    for &byte in bytes {
        match byte {
            b'\\' => out.push_str("\\\\"),
            b'\t' => out.push_str("\\t"),
            b'\n' => out.push_str("\\n"),
            b'\r' => out.push_str("\\r"),
            byte if byte == quote => {
                out.push('\\');
                out.push(char::from(byte));
            }
            b' '..=b'~' => out.push(char::from(byte)),
            byte => out.push_str(&format!("\\x{byte:02x}")),
        }
    }
    out.push(char::from(quote));
    out
}

/// What the `gzip` module says of a file that ends inside a member.
const GZIP_ENDED: &str = "Compressed file ended before the end-of-stream marker was reached";

/// A gzip file as Python's `gzip` module reads it: every member in turn, zero bytes after
/// a member skipped, a file of no bytes holding nothing. Where the module raises, the read
/// fails with a [`PyError`] of its text: a file that ends inside a member, bytes that begin
/// no member, an unknown compression method, a checksum or a length that the data does not
/// have. Deflate data that cannot be decoded fails in this reader's own words.
pub struct PyGzip<R> {
    input: R,
    /// Bytes of the file read and not used, from `at` on.
    held: Vec<u8>,
    at: usize,
    /// The file has no byte after those held.
    ended: bool,
    in_member: bool,
    done: bool,
    inflater: Decompress,
    crc: Crc,
}

impl<R: Read> PyGzip<R> {
    pub fn new(input: R) -> PyGzip<R> {
        PyGzip {
            input,
            held: Vec::new(),
            at: 0,
            ended: false,
            in_member: false,
            done: false,
            inflater: Decompress::new(false),
            crc: Crc::new(),
        }
    }

    /// Reads more of the file behind the bytes held; false at its end.
    fn fill(&mut self) -> io::Result<bool> {
        if self.ended {
            return Ok(false);
        }
        self.held.drain(..self.at);
        self.at = 0;
        let before = self.held.len();
        self.held.resize(before + (1 << 16), 0);
        let read = loop {
            match self.input.read(&mut self.held[before..]) {
                Ok(read) => break read,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => {
                    self.held.truncate(before);
                    return Err(e);
                }
            }
        };
        self.held.truncate(before + read);
        self.ended = read == 0;
        Ok(read > 0)
    }

    /// The next `count` bytes of the file, fewer only where it ends.
    fn take(&mut self, count: usize) -> io::Result<Vec<u8>> {
        while self.held.len() - self.at < count && self.fill()? {}
        let found = count.min(self.held.len() - self.at);
        let bytes = self.held[self.at..self.at + found].to_vec();
        self.at += found;
        Ok(bytes)
    }

    /// The next `count` bytes of the file, which a member needs whole.
    fn exactly(&mut self, count: usize) -> io::Result<Vec<u8>> {
        let bytes = self.take(count)?;
        if bytes.len() < count {
            return Err(py_error(GZIP_ENDED));
        }
        Ok(bytes)
    }

    /// Reads the header of a member; false where the file ends before one begins.
    fn header(&mut self) -> io::Result<bool> {
        const HEADER_CHECKSUM: u8 = 2;
        const EXTRA: u8 = 4;
        const NAME: u8 = 8;
        const COMMENT: u8 = 16;
        let magic = self.take(2)?;
        if magic.is_empty() {
            return Ok(false);
        }
        if magic != [0x1f, 0x8b] {
            return Err(py_error(format!(
                "Not a gzipped file ({})",
                py_bytes_repr(&magic)
            )));
        }
        let fixed = self.exactly(8)?;
        if fixed[0] != 8 {
            return Err(py_error("Unknown compression method"));
        }
        let flags = fixed[1];
        if flags & EXTRA != 0 {
            let length = self.exactly(2)?;
            self.exactly(usize::from(u16::from_le_bytes([length[0], length[1]])))?;
        }
        // A name or a comment runs to its zero byte, or to the end of the file.
        for text in [NAME, COMMENT] {
            while flags & text != 0 && self.take(1)?.first().is_some_and(|byte| *byte != 0) {}
        }
        if flags & HEADER_CHECKSUM != 0 {
            self.exactly(2)?;
        }
        Ok(true)
    }

    /// Reads the checksum and the length that close a member, holds the data read to them,
    /// and skips the zero bytes after.
    fn trailer(&mut self) -> io::Result<()> {
        let end = self.exactly(8)?;
        let stored = u32::from_le_bytes([end[0], end[1], end[2], end[3]]);
        let length = u32::from_le_bytes([end[4], end[5], end[6], end[7]]);
        if stored != self.crc.sum() {
            return Err(py_error(format!(
                "CRC check failed {stored:#x} != {:#x}",
                self.crc.sum()
            )));
        }
        if length != self.crc.amount() {
            return Err(py_error("Incorrect length of data produced"));
        }
        loop {
            match self.take(1)?.first() {
                Some(0) => {}
                Some(_) => {
                    self.at -= 1;
                    return Ok(());
                }
                None => return Ok(()),
            }
        }
    }
}

impl<R: Read> Read for PyGzip<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        loop {
            if self.done {
                return Ok(0);
            }
            if !self.in_member {
                if !self.header()? {
                    self.done = true;
                    return Ok(0);
                }
                self.inflater.reset(false);
                self.crc.reset();
                self.in_member = true;
            }
            if self.at == self.held.len() && !self.fill()? {
                return Err(py_error(GZIP_ENDED));
            }
            let (read, written) = (self.inflater.total_in(), self.inflater.total_out());
            let status = self
                .inflater
                .decompress(&self.held[self.at..], out, FlushDecompress::None)
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "corrupt deflate data"))?;
            let read = (self.inflater.total_in() - read) as usize;
            let written = (self.inflater.total_out() - written) as usize;
            self.at += read;
            self.crc.update(&out[..written]);
            if status == Status::StreamEnd {
                self.trailer()?;
                self.in_member = false;
            } else if read == 0 && written == 0 && !self.fill()? {
                return Err(py_error(GZIP_ENDED));
            }
            if written > 0 {
                return Ok(written);
            }
        }
    }
}

/// Opens an input VCF as Python's readers of its records do: through the `gzip` module
/// where its name has the suffix `.gz` or `.bgz`, as a plain file otherwise.
pub fn open_records(path: &Path) -> io::Result<Box<dyn Read>> {
    let file = open_file(path)?;
    let gzipped = matches!(
        py_suffix(py_name(path.as_os_str().as_bytes())),
        b".gz" | b".bgz"
    );
    Ok(if gzipped {
        Box::new(PyGzip::new(file))
    } else {
        Box::new(file)
    })
}

/// The failure an operating-system error ends the comparison with, in the words of
/// Python's `OSError`: `[Errno N] description`, and after it `: ` and the `repr()` of the
/// path when the failed call named one. A [`PyError`] is its text alone. Any other error
/// keeps its own words, after the path when there is one.
pub fn os_failure(error: &io::Error, path: Option<&Path>) -> Failure {
    let Some(code) = error.raw_os_error() else {
        if let Some(PyError(text)) = error.get_ref().and_then(|inner| inner.downcast_ref()) {
            return Failure::io(text.clone());
        }
        return Failure::io(match path {
            Some(path) => format!("{}: {error}", path.display()),
            None => error.to_string(),
        });
    };
    let text = error.to_string();
    let description = text
        .strip_suffix(&format!(" (os error {code})"))
        .unwrap_or(&text);
    Failure::io(match path {
        Some(path) => format!(
            "[Errno {code}] {description}: {}",
            py_repr(path.as_os_str().as_bytes())
        ),
        None => format!("[Errno {code}] {description}"),
    })
}

/// `Path.mkdir(parents=True, exist_ok=True)`: makes `path` and the directories above it
/// that are missing. The failure names the directory whose creation the system refused.
pub fn make_dirs(path: &Path) -> Result<(), Failure> {
    fn make(path: &Path) -> Result<(), (PathBuf, io::Error)> {
        let refused = |error: io::Error| {
            if path.is_dir() {
                Ok(())
            } else {
                Err((path.to_path_buf(), error))
            }
        };
        match std::fs::create_dir(path) {
            Ok(()) => Ok(()),
            Err(error) if error.kind() == io::ErrorKind::NotFound => {
                // The directory above a name of one component is the working directory.
                let parent = match path.parent() {
                    Some(parent) if parent.as_os_str().is_empty() => Path::new("."),
                    Some(parent) => parent,
                    None => return Err((path.to_path_buf(), error)),
                };
                if parent == path {
                    return Err((path.to_path_buf(), error));
                }
                make(parent)?;
                match std::fs::create_dir(path) {
                    Ok(()) => Ok(()),
                    Err(error) if error.kind() == io::ErrorKind::NotFound => {
                        Err((path.to_path_buf(), error))
                    }
                    Err(error) => refused(error),
                }
            }
            Err(error) => refused(error),
        }
    }
    make(path).map_err(|(at, error)| os_failure(&error, Some(&at)))
}

/// Where CPython's UTF-8 decoder stops in `data`: the count of bytes it decodes, short of
/// a sequence the end of `data` cuts when more input may follow (`last` is false), or the
/// first and the one-past-last index of the bytes it refuses, with its reason.
fn py_utf8_decode(data: &[u8], last: bool) -> Result<usize, (usize, usize, &'static str)> {
    const START: &str = "invalid start byte";
    const CONTINUATION: &str = "invalid continuation byte";
    let continues = |byte: u8| (0x80..=0xbf).contains(&byte);
    let mut at = 0;
    while at < data.len() {
        let lead = data[at];
        let left = data.len() - at;
        let second = data.get(at + 1).copied();
        let width = match lead {
            0..=0x7f => 1,
            0x80..=0xc1 | 0xf5..=0xff => return Err((at, at + 1, START)),
            0xc2..=0xdf => {
                let Some(second) = second else { break };
                if !continues(second) {
                    return Err((at, at + 1, CONTINUATION));
                }
                2
            }
            0xe0..=0xef => {
                let Some(second) = second else { break };
                let refused = !continues(second)
                    || (lead == 0xe0 && second < 0xa0)
                    || (lead == 0xed && second >= 0xa0);
                if refused {
                    // A surrogate cut after its second byte waits for more input.
                    if left == 2 && !last && lead == 0xed && continues(second) {
                        break;
                    }
                    return Err((at, at + 1, CONTINUATION));
                }
                if left < 3 {
                    break;
                }
                if !continues(data[at + 2]) {
                    return Err((at, at + 2, CONTINUATION));
                }
                3
            }
            0xf0..=0xf4 => {
                let Some(second) = second else { break };
                if !continues(second)
                    || (lead == 0xf0 && second < 0x90)
                    || (lead == 0xf4 && second >= 0x90)
                {
                    return Err((at, at + 1, CONTINUATION));
                }
                if left < 3 {
                    break;
                }
                if !continues(data[at + 2]) {
                    return Err((at, at + 2, CONTINUATION));
                }
                if left < 4 {
                    break;
                }
                if !continues(data[at + 3]) {
                    return Err((at, at + 3, CONTINUATION));
                }
                4
            }
        };
        at += width;
    }
    if last && at < data.len() {
        return Err((at, data.len(), "unexpected end of data"));
    }
    Ok(at)
}

/// How many lines the terminators of `text` end, a carriage return carried over from the
/// text before it counted with them when `carried` is set.
fn terminated_lines(carried: bool, text: &[u8]) -> u64 {
    let mut count = u64::from(carried && text.first() != Some(&b'\n'));
    for i in memchr2_iter(b'\n', b'\r', text) {
        if text[i] == b'\n' || text.get(i + 1) != Some(&b'\n') {
            count += 1;
        }
    }
    count
}

/// The decoding failure Python's text reader ends on over the bytes of `reader`, with the
/// number of lines it hands out before: the text of its `UnicodeDecodeError`, whose
/// position counts from the start of the bytes the decoder was given, a chunk of 8,192
/// bytes of the stream after the bytes of a sequence the chunk before cut. `None` for a
/// stream that is UTF-8 throughout.
///
/// A line is handed out once its terminator is decoded, and a carriage return that ends a
/// chunk's text once the text after it is: the failure comes before every line whose
/// terminator stands in the failing chunk or after it.
pub fn py_decode_failure(mut reader: impl Read) -> io::Result<Option<(u64, String)>> {
    let mut data: Vec<u8> = Vec::with_capacity(PY_CHUNK_BYTES + 4);
    let mut carried = false;
    let mut lines = 0;
    loop {
        let pending = data.len();
        data.resize(pending + PY_CHUNK_BYTES, 0);
        let mut filled = pending;
        while filled < data.len() {
            match reader.read(&mut data[filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        data.truncate(filled);
        let last = filled == pending;
        let decoded = match py_utf8_decode(&data, last) {
            Ok(decoded) => decoded,
            Err((start, end, reason)) => {
                let message = if end == start + 1 {
                    format!(
                        "'utf-8' codec can't decode byte 0x{:02x} in position {start}: {reason}",
                        data[start]
                    )
                } else {
                    format!(
                        "'utf-8' codec can't decode bytes in position {start}-{}: {reason}",
                        end - 1
                    )
                };
                return Ok(Some((lines, message)));
            }
        };
        if last {
            return Ok(None);
        }
        let mut text = &data[..decoded];
        let prefixed = carried && !text.is_empty();
        if prefixed {
            carried = false;
        }
        // A carriage return that ends the text is held until the text after it shows
        // whether a line feed follows.
        if let Some((b'\r', before)) = text.split_last() {
            text = before;
            carried = true;
        }
        lines += terminated_lines(prefixed, text);
        data.drain(..decoded);
    }
}

/// The failure Python's reading of a file ends in, given the first defective line the
/// comparison found in it (`line`, counted from 1) and the bytes of the file as they are
/// read: the decoding failure of [`py_decode_failure`] when the reader meets it before it
/// hands out that line, else the failure of the line.
pub fn reading_failure(
    path: &Path,
    stream: io::Result<impl Read>,
    line: u64,
    defect: LineDefect,
) -> Failure {
    if let Ok(Ok(Some((before, message)))) = stream.map(py_decode_failure) {
        if line > before {
            return Failure::decode(message);
        }
    }
    line_failure(path, line, defect)
}

/// `f"{n:,}"`: the decimal digits of `n` in groups of three from the right, a comma
/// between groups.
pub fn py_grouped(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            out.push(',');
        }
        out.push(digit);
    }
    out
}

/// `str.isalpha()`: at least one character, every one a letter (general category `L`).
pub fn py_isalpha(text: &str) -> bool {
    static LETTERS: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^\p{L}+$").expect("regex"));
    LETTERS.is_match(text)
}

/// `str(Path(arg))` as a path.
pub fn py_path(arg: &OsStr) -> PathBuf {
    PathBuf::from(OsString::from_vec(crate::pyfmt::py_path_str(
        arg.as_bytes(),
    )))
}

/// `Path(directory) / name` for a `directory` already in Python's spelling and a `name` of
/// one component: no `./` before the name under the directory `.`.
pub fn py_join(directory: &Path, name: &OsStr) -> PathBuf {
    let directory = directory.as_os_str().as_bytes();
    let mut joined = Vec::with_capacity(directory.len() + 1 + name.len());
    if directory != b"." {
        joined.extend_from_slice(directory);
        if !directory.ends_with(b"/") {
            joined.push(b'/');
        }
    }
    joined.extend_from_slice(name.as_bytes());
    PathBuf::from(OsString::from_vec(joined))
}

/// `PurePath.name` of a path: its final component, `..` a component as any other and `.`
/// and an empty one none; empty for the root and for the working directory.
pub fn py_name(path: &[u8]) -> &[u8] {
    path.rsplit(|byte| *byte == b'/')
        .find(|part| !part.is_empty() && *part != b".")
        .unwrap_or(b"")
}

/// `PurePath.stem` of a final path component: the name without its last suffix, a name
/// whose text before the last dot is empty or all dots having no suffix.
pub fn py_stem(name: &[u8]) -> &[u8] {
    match memrchr(b'.', name) {
        Some(i) if name[..i].iter().any(|b| *b != b'.') => &name[..i],
        _ => name,
    }
}

/// `PurePath.suffix` of a final path component: what [`py_stem`] drops of the name, from
/// its last dot on; empty for a name without a suffix.
pub fn py_suffix(name: &[u8]) -> &[u8] {
    &name[py_stem(name).len()..]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::failure::FailureKind;
    use crate::tempdir::TempDir;

    fn lines_of(bytes: &[u8]) -> Vec<Vec<u8>> {
        lines(bytes).map(<[u8]>::to_vec).collect()
    }

    /// A reader that hands out at most `step` bytes a call.
    struct Trickle<'a> {
        bytes: &'a [u8],
        step: usize,
    }

    impl Read for Trickle<'_> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = self.step.min(buf.len()).min(self.bytes.len());
            buf[..n].copy_from_slice(&self.bytes[..n]);
            self.bytes = &self.bytes[n..];
            Ok(n)
        }
    }

    #[test]
    fn every_newline_convention_ends_a_line() {
        assert_eq!(
            lines_of(b"a\nb\r\nc\rd\r\r\ne\n\nf"),
            [&b"a"[..], b"b", b"c", b"d", b"", b"e", b"", b"f"]
        );
        assert_eq!(lines_of(b"a\r"), [b"a"]);
        assert_eq!(lines_of(b"\n"), [b""]);
        assert!(lines_of(b"").is_empty());
        for text in [
            &b"a\nb\r\nc\rd\r\r\ne\n\nf"[..],
            b"a\r",
            b"\n",
            b"",
            b"\r\n\r\n",
            b"\n\r\n\r",
            b"x",
        ] {
            assert_eq!(line_count(text), lines_of(text).len() as u64, "{text:?}");
        }
    }

    #[test]
    fn blocks_of_any_size_carry_the_lines_of_the_whole_input() {
        let text: &[u8] = b"one\ntwo\r\nthree\rfour\r\n\r\nfive\n\rsix\r\n\nseven\r\r\reight";
        let whole = lines_of(text);
        for block_bytes in 1..=text.len() + 2 {
            for step in [1, 2, 3, 7, 64] {
                let mut read: Vec<Vec<u8>> = Vec::new();
                let mut counted = 0;
                let mut blocks = 0;
                read_blocks(Trickle { bytes: text, step }, block_bytes, |block| {
                    blocks += 1;
                    counted += line_count(&block);
                    read.extend(lines_of(&block));
                })
                .unwrap();
                assert_eq!(read, whole, "block {block_bytes} step {step}");
                assert_eq!(counted, whole.len() as u64);
                assert!(blocks >= 1);
            }
        }
        let mut none = 0;
        read_blocks(&b""[..], 8, |_| none += 1).unwrap();
        assert_eq!(none, 0);
    }

    #[test]
    fn a_failing_reader_ends_where_it_fails() {
        struct Failing(usize);
        impl Read for Failing {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.0 == 0 {
                    return Err(io::Error::other("corrupt"));
                }
                self.0 -= 1;
                buf[..2].copy_from_slice(b"a\n");
                Ok(2)
            }
        }
        let mut read = Vec::new();
        read_blocks(UntilFailure::new(Failing(3)), 64, |block| {
            read.extend(block)
        })
        .unwrap();
        assert_eq!(read, b"a\na\na\n");
        assert!(read_blocks(Failing(1), 64, |_| {}).is_err());
    }

    /// `at` bytes of lines of 99 letters, then `bad`, then `tail`.
    fn stream(at: usize, bad: &[u8], tail: &[u8]) -> Vec<u8> {
        let mut bytes: Vec<u8> = (0..at)
            .map(|i| if i % 100 == 99 { b'\n' } else { b'a' })
            .collect();
        bytes.extend_from_slice(bad);
        bytes.extend_from_slice(tail);
        bytes
    }

    #[test]
    fn a_decoding_failure_reads_as_the_interpreter_writes_it() {
        type Case = (
            usize,
            &'static [u8],
            &'static [u8],
            Option<(u64, &'static str)>,
        );
        let cases: [Case; 21] = [
            (0, b"\xe9", b" x\n", Some((0, "'utf-8' codec can't decode byte 0xe9 in position 0: invalid continuation byte"))),
            (5, b"\xff", b"", Some((0, "'utf-8' codec can't decode byte 0xff in position 5: invalid start byte"))),
            (250, b"\xe2\x82\x28", b"\n", Some((0, "'utf-8' codec can't decode bytes in position 250-251: invalid continuation byte"))),
            (250, b"\xf0\x9f\x98\x28", b"\n", Some((0, "'utf-8' codec can't decode bytes in position 250-252: invalid continuation byte"))),
            (250, b"\xc0\xaf", b"\n", Some((0, "'utf-8' codec can't decode byte 0xc0 in position 250: invalid start byte"))),
            (250, b"\xed\xa0\x80", b"\n", Some((0, "'utf-8' codec can't decode byte 0xed in position 250: invalid continuation byte"))),
            // A sequence the end of the file cuts: the lines before it are handed out first.
            (250, b"\xe2\x82", b"", Some((2, "'utf-8' codec can't decode bytes in position 0-1: unexpected end of data"))),
            (250, b"\xf0", b"", Some((2, "'utf-8' codec can't decode byte 0xf0 in position 0: unexpected end of data"))),
            // A sequence a chunk cuts is decoded with the next chunk, from position 0.
            (8191, b"\xc3", b"\x28", Some((81, "'utf-8' codec can't decode byte 0xc3 in position 0: invalid continuation byte"))),
            (8191, b"\xc3\xa9", b"\n\xff", Some((81, "'utf-8' codec can't decode byte 0xff in position 3: invalid start byte"))),
            (8190, b"\xe2\x82", b"", Some((81, "'utf-8' codec can't decode bytes in position 0-1: unexpected end of data"))),
            (8192, b"\x80", b"", Some((81, "'utf-8' codec can't decode byte 0x80 in position 0: invalid start byte"))),
            // A carriage return that ends a chunk ends no line before the next chunk.
            (8191, b"\r", b"\xff", Some((81, "'utf-8' codec can't decode byte 0xff in position 0: invalid start byte"))),
            (8190, b"\rb", b"\xff", Some((82, "'utf-8' codec can't decode byte 0xff in position 0: invalid start byte"))),
            (8191, b"\r", b"\n\xff", Some((81, "'utf-8' codec can't decode byte 0xff in position 1: invalid start byte"))),
            // The first two bytes of a surrogate wait for more input and are refused with it.
            (8191, b"\xed\xa0", b"\x80", Some((81, "'utf-8' codec can't decode byte 0xed in position 0: invalid continuation byte"))),
            (8190, b"\xed\xa0", b"\x80", Some((81, "'utf-8' codec can't decode byte 0xed in position 0: invalid continuation byte"))),
            (8190, b"\xed\xa0", b"", Some((81, "'utf-8' codec can't decode byte 0xed in position 0: invalid continuation byte"))),
            (8189, b"\xf0\x9f\x98", b"\x28", Some((81, "'utf-8' codec can't decode bytes in position 0-2: invalid continuation byte"))),
            (100_000, b"\xe9", b"A", Some((983, "'utf-8' codec can't decode byte 0xe9 in position 1696: invalid continuation byte"))),
            (16_384, b"", b"no failure here\n", None),
        ];
        for (at, bad, tail, expected) in cases {
            let bytes = stream(at, bad, tail);
            for step in [1, 7, 8192, 100_000] {
                let failure = py_decode_failure(Trickle {
                    bytes: &bytes,
                    step,
                })
                .unwrap();
                let failure = failure
                    .as_ref()
                    .map(|(lines, message)| (*lines, message.as_str()));
                assert_eq!(failure, expected, "{at} {bad:?} {tail:?} step {step}");
            }
        }
    }

    #[test]
    fn the_failure_of_a_file_is_the_one_its_reading_meets_first() {
        let short = LineDefect::ShortRow(2);
        let path = Path::new("out.txt");
        let row = line_failure(path, 2, short);
        // The whole first chunk is decoded before its second line is handed out.
        let mut bytes = stream(9000, b"", b"");
        bytes[300] = 0xff;
        let failure = reading_failure(path, Ok(&bytes[..]), 2, short);
        assert_eq!(failure.kind, FailureKind::Decode);
        assert_eq!(
            failure.message,
            "'utf-8' codec can't decode byte 0xff in position 300: invalid start byte"
        );
        // A line of the first chunk is handed out before the second chunk is decoded.
        let mut bytes = stream(9000, b"", b"");
        bytes[8500] = 0xff;
        assert_eq!(reading_failure(path, Ok(&bytes[..]), 2, short), row);
        assert_eq!(
            reading_failure(path, Ok(&bytes[..]), 81, short).kind,
            FailureKind::Value
        );
        assert_eq!(
            reading_failure(path, Ok(&bytes[..]), 82, short).message,
            "'utf-8' codec can't decode byte 0xff in position 308: invalid start byte"
        );
        // A file that is UTF-8 throughout, and one that cannot be read again.
        let bytes = stream(9000, b"", b"");
        assert_eq!(reading_failure(path, Ok(&bytes[..]), 2, short), row);
        let unread: io::Result<&[u8]> = Err(io::Error::other("gone"));
        assert_eq!(reading_failure(path, unread, 2, short), row);
        assert_eq!(
            reading_failure(path, Ok(&b"a\n\xff\n"[..]), 2, LineDefect::NotUtf8).message,
            "'utf-8' codec can't decode byte 0xff in position 2: invalid start byte"
        );
    }

    #[test]
    fn a_path_is_quoted_as_python_quotes_a_text() {
        let cases: [(&[u8], &str); 18] = [
            (b"/a/b.txt", "'/a/b.txt'"),
            (b"it's", "\"it's\""),
            (b"a'b\"c", "'a\\'b\"c'"),
            (b"tab\there", "'tab\\there'"),
            (b"nl\n\r", "'nl\\n\\r'"),
            (b"back\\slash", "'back\\\\slash'"),
            ("\u{e9}".as_bytes(), "'\u{e9}'"),
            ("\u{a0}".as_bytes(), "'\\xa0'"),
            ("\u{200b}".as_bytes(), "'\\u200b'"),
            ("\u{e000}".as_bytes(), "'\\ue000'"),
            ("\u{378}".as_bytes(), "'\\u0378'"),
            ("\u{1f600}".as_bytes(), "'\u{1f600}'"),
            ("\u{e0001}".as_bytes(), "'\\U000e0001'"),
            (b"\x7f\x01", "'\\x7f\\x01'"),
            (b"caf\xe9", "'caf\\udce9'"),
            ("\u{2028}".as_bytes(), "'\\u2028'"),
            ("\u{3000}".as_bytes(), "'\\u3000'"),
            (b" ~", "' ~'"),
        ];
        for (bytes, expected) in cases {
            assert_eq!(py_repr(bytes), expected, "{bytes:?}");
        }
    }

    #[test]
    fn a_system_error_reads_as_python_writes_it() {
        let dir = TempDir::new("text_os_").unwrap();
        let absent = dir.path.join("absent.txt");
        let error = open_file(&absent).unwrap_err();
        let failure = os_failure(&error, Some(&absent));
        assert_eq!(failure.kind, FailureKind::Io);
        assert_eq!(
            failure.message,
            format!(
                "[Errno 2] No such file or directory: '{}'",
                absent.display()
            )
        );
        assert_eq!(
            os_failure(&error, None).message,
            "[Errno 2] No such file or directory"
        );
        let folder = open_file(&dir.path).unwrap_err();
        assert_eq!(
            os_failure(&folder, Some(&dir.path)).message,
            format!("[Errno 21] Is a directory: '{}'", dir.path.display())
        );
        let other = io::Error::other("corrupt");
        assert_eq!(
            os_failure(&other, Some(Path::new("a/b"))).message,
            "a/b: corrupt"
        );
        assert_eq!(os_failure(&other, None).message, "corrupt");
    }

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        use std::io::Write;
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    /// What the reader of gzip files gives for `bytes`: the data, or the failure it ends in.
    fn gunzip(bytes: &[u8], step: usize) -> Result<Vec<u8>, String> {
        let mut data = Vec::new();
        PyGzip::new(Trickle { bytes, step })
            .read_to_end(&mut data)
            .map(|_| data)
            .map_err(|error| os_failure(&error, Some(Path::new("in.vcf.gz"))).message)
    }

    #[test]
    fn a_gzip_file_is_read_and_refused_as_the_gzip_module_reads_and_refuses_it() {
        let payload: &[u8] = b"21\t100\tid\tN\t<DEL>\n";
        let whole = gzip(payload);
        let twice = [payload, payload].concat();
        let ended = "Compressed file ended before the end-of-stream marker was reached";
        let with = |tail: &[u8]| [&whole[..], tail].concat();
        let set = |range: std::ops::Range<usize>, bytes: &[u8]| {
            let mut changed = whole.clone();
            changed[range].copy_from_slice(bytes);
            changed
        };
        let end = whole.len();
        type Read = Result<Vec<u8>, String>;
        let cases: Vec<(Vec<u8>, Read)> = vec![
            (whole.clone(), Ok(payload.to_vec())),
            (with(&whole), Ok(twice.clone())),
            (Vec::new(), Ok(Vec::new())),
            (with(&[0; 7]), Ok(payload.to_vec())),
            ([&whole[..], &[0; 7], &whole[..]].concat(), Ok(twice)),
            (whole[..end - 12].to_vec(), Err(ended.to_string())),
            (whole[..end - 3].to_vec(), Err(ended.to_string())),
            (whole[..5].to_vec(), Err(ended.to_string())),
            (
                b"##fileformat\n".to_vec(),
                Err("Not a gzipped file (b'##')".to_string()),
            ),
            (b"#".to_vec(), Err("Not a gzipped file (b'#')".to_string())),
            (
                with(b"\x1f"),
                Err("Not a gzipped file (b'\\x1f')".to_string()),
            ),
            (
                with(b"'\""),
                Err("Not a gzipped file (b'\\'\"')".to_string()),
            ),
            (with(b"'x"), Err("Not a gzipped file (b\"'x\")".to_string())),
            (
                with(b"\\\n"),
                Err("Not a gzipped file (b'\\\\\\n')".to_string()),
            ),
            (
                with(b"\xff\x7f"),
                Err("Not a gzipped file (b'\\xff\\x7f')".to_string()),
            ),
            (
                set(end - 8..end - 4, &[0, 1, 2, 3]),
                Err("CRC check failed 0x3020100 != 0x6a77ffb8".to_string()),
            ),
            (
                set(end - 4..end, &[1, 0, 0, 0]),
                Err("Incorrect length of data produced".to_string()),
            ),
            (
                set(2..3, &[7]),
                Err("Unknown compression method".to_string()),
            ),
        ];
        for (bytes, expected) in cases {
            for step in [1, 5, 1 << 16] {
                assert_eq!(gunzip(&bytes, step), expected, "{bytes:?} step {step}");
            }
        }
        // Deflate data that cannot be decoded is refused in the reader's own words, the
        // file named.
        let mut corrupt = whole.clone();
        corrupt[10] = 0xff;
        assert_eq!(
            gunzip(&corrupt, 5),
            Err("in.vcf.gz: corrupt deflate data".to_string())
        );
    }

    #[test]
    fn an_input_is_opened_by_its_name() {
        let dir = TempDir::new("text_open_").unwrap();
        let read = |name: &str, bytes: &[u8]| {
            let path = dir.path.join(name);
            std::fs::write(&path, bytes).unwrap();
            let mut data = Vec::new();
            open_records(&path)
                .and_then(|mut reader| reader.read_to_end(&mut data))
                .map(|_| data)
                .map_err(|error| os_failure(&error, Some(&path)).message)
        };
        let record: &[u8] = b"21\t5\ta\tA\tT\n";
        assert_eq!(read("in.vcf", record), Ok(record.to_vec()));
        assert_eq!(read("in.vcf.gz", &gzip(record)), Ok(record.to_vec()));
        assert_eq!(read("in.vcf.bgz", &gzip(record)), Ok(record.to_vec()));
        assert_eq!(read("gz.vcf", &gzip(record)), Ok(gzip(record)));
        // A name of dots before `gz` or `bgz` has no suffix, and one that begins with dots
        // and has other text before its last dot has one.
        for plain in ["...gz", "..gz", ".gz", "...bgz", ".bgz", "in.vcf.gz."] {
            assert_eq!(read(plain, &gzip(record)), Ok(gzip(record)), "{plain}");
        }
        for gzipped in ["..a.gz", ".vcf.gz", "..vcf.bgz", "a..gz"] {
            assert_eq!(
                read(gzipped, &gzip(record)),
                Ok(record.to_vec()),
                "{gzipped}"
            );
        }
        assert_eq!(
            read("plain.vcf.gz", record),
            Err("Not a gzipped file (b'21')".to_string())
        );
        let folder = dir.path.join("folder.vcf.gz");
        std::fs::create_dir(&folder).unwrap();
        assert_eq!(
            open_records(&folder)
                .map(|_| ())
                .map_err(|error| os_failure(&error, Some(&folder)).message),
            Err(format!("[Errno 21] Is a directory: '{}'", folder.display()))
        );
    }

    #[test]
    fn directories_are_made_as_pathlib_makes_them() {
        let dir = TempDir::new("text_mkdir_").unwrap();
        let deep = dir.path.join("a/b/c");
        make_dirs(&deep).unwrap();
        assert!(deep.is_dir());
        make_dirs(&deep).unwrap();
        let file = dir.path.join("file");
        std::fs::write(&file, b"x").unwrap();
        assert_eq!(
            make_dirs(&file).unwrap_err().message,
            format!("[Errno 17] File exists: '{}'", file.display())
        );
        // The directory named is the one the system refused, here the deepest.
        let under = file.join("x/y");
        assert_eq!(
            make_dirs(&under).unwrap_err().message,
            format!("[Errno 20] Not a directory: '{}'", under.display())
        );
    }

    #[test]
    fn a_count_is_grouped_in_threes() {
        let cases = [
            (0, "0"),
            (999, "999"),
            (1000, "1,000"),
            (16731, "16,731"),
            (1_234_567, "1,234,567"),
            (100_000, "100,000"),
            (u64::MAX, "18,446,744,073,709,551,615"),
        ];
        for (n, expected) in cases {
            assert_eq!(py_grouped(n), expected);
        }
    }

    #[test]
    fn letters_are_the_letter_categories_alone() {
        assert!(py_isalpha("A"));
        assert!(py_isalpha("N\u{e9}\u{ff41}"));
        assert!(!py_isalpha(""));
        assert!(!py_isalpha("A1"));
        // A letter number and a combining mark are alphabetic and are no letters.
        assert!(!py_isalpha("\u{2167}"));
        assert!(!py_isalpha("a\u{301}"));
    }

    #[test]
    fn paths_are_spelled_as_pathlib_spells_them() {
        let join = |directory: &str, name: &str| {
            py_join(&py_path(OsStr::new(directory)), OsStr::new(name))
                .to_str()
                .unwrap()
                .to_string()
        };
        assert_eq!(join("", "sv.txt"), "sv.txt");
        assert_eq!(join(".", "sv.txt"), "sv.txt");
        assert_eq!(join("./", ".txt"), ".txt");
        assert_eq!(join("a//b/", "sv.txt"), "a/b/sv.txt");
        assert_eq!(join("/", "sv.txt"), "/sv.txt");
        assert_eq!(join("//", "sv.txt"), "//sv.txt");
        assert_eq!(join("///", "sv.txt"), "/sv.txt");
        assert_eq!(join("a/../b", "sv.txt"), "a/../b/sv.txt");
        assert_eq!(join("./a", "sv.txt"), "a/sv.txt");
    }

    #[test]
    fn a_stem_drops_the_last_suffix() {
        let cases = [
            ("01_snv.vcf.gz", "01_snv.vcf"),
            ("x.vcf", "x"),
            (".vcf.gz", ".vcf"),
            (".vcf", ".vcf"),
            ("...vcf", "...vcf"),
            ("..vcf.gz", "..vcf"),
            ("a.", "a"),
            ("a..", "a."),
            (".a.", ".a"),
            ("..", ".."),
            ("...", "..."),
            ("a", "a"),
            ("..a.vcf", "..a"),
        ];
        for (name, stem) in cases {
            assert_eq!(py_stem(name.as_bytes()), stem.as_bytes(), "{name}");
        }
    }

    #[test]
    fn a_name_is_the_last_component_and_a_suffix_what_the_stem_drops() {
        let names = [
            ("a/..", ".."),
            ("..", ".."),
            ("a/../..", ".."),
            ("a/.", "a"),
            ("a/b/./", "b"),
            (".", ""),
            ("", ""),
            ("/", ""),
            ("//", ""),
            ("a/", "a"),
            ("a//b.vcf", "b.vcf"),
            ("/x/...gz", "...gz"),
        ];
        for (path, name) in names {
            assert_eq!(py_name(path.as_bytes()), name.as_bytes(), "{path}");
        }
        let suffixes = [
            ("...gz", ""),
            ("..gz", ""),
            (".gz", ""),
            ("a.gz", ".gz"),
            ("..a.gz", ".gz"),
            ("...bgz", ""),
            ("a.vcf.bgz", ".bgz"),
            ("foo.", "."),
            ("foo.gz.", "."),
            (".vcf.gz", ".gz"),
            ("a..gz", ".gz"),
            ("....", ""),
            ("..", ""),
            ("vcf", ""),
            ("", ""),
        ];
        for (name, suffix) in suffixes {
            assert_eq!(py_suffix(name.as_bytes()), suffix.as_bytes(), "{name}");
        }
    }

    #[test]
    fn the_integer_limit_names_its_digits() {
        assert_eq!(
            int_conversion_failure(4301).message,
            "Exceeds the limit (4300 digits) for integer string conversion: value has 4301 digits; use sys.set_int_max_str_digits() to increase the limit"
        );
    }
}
