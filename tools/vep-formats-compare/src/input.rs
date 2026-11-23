// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reading a file as text, its lines and the input records.

use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::parse::vcf::record_prefix;
use crate::pyfmt;
use crate::pyjson::{self, Value};
use crate::report::strerror;
use crate::stream::gzip::{Failure, Gunzip, Open};
use crate::stream::utf8;
use flate2::{Crc, Decompress, FlushDecompress, Status};
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fmt::Write as _;
use std::fs::File;
use std::io::Read;
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};

/// The call depth `canonical_json_text` runs at where a comparison parses a
/// default-format output. An object is refused where the calls that write its canonical
/// form (`pyjson::compact_frames`) would pass the interpreter's recursion limit
/// (`pyjson::RECURSION_LIMIT`) less this depth: objects around a whole number are read
/// 495 deep and refused at 496.
pub const CANONICAL_CALLERS: usize = 6;

/// The call depth `numberified_json_text` runs at, in the self-parity mode. An object
/// ends the run where the calls that walk it and write it would pass the interpreter's
/// recursion limit less this depth: objects around a whole number are read 495 deep and
/// end the run at 496.
const NUMBERIFIED_CALLERS: usize = 6;

/// The room the inflater is given at a call.
const INFLATE_ROOM: usize = 1 << 20;

/// Less room than the inflater's fast path asks for (260 bytes).
const SLOW_INFLATE_ROOM: usize = 64;

/// What the inflater says of a stream whose error its fast path found: the fast path's
/// own reason is replaced by this one before the call returns.
const REASON_LOST: &str = "repeated call with bad state";

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const GZIP_DEFLATE: u8 = 8;
const FHCRC: u8 = 2;
const FEXTRA: u8 = 4;
const FNAME: u8 = 8;
const FCOMMENT: u8 = 16;
const GZIP_CUT_SHORT: &str = "Compressed file ended before the end-of-stream marker was reached";

/// The text of `path`, gunzipped when its name ends in `.gz`.
///
/// The text is read as the interpreter's text mode reads it: every `\r\n` and every other
/// `\r` is a `\n`. A gzip file is read whole before its text is decoded, so a file that is
/// both cut short and not UTF-8 is reported as cut short. A path that holds a NUL ends
/// the run and is no comparator's error.
pub fn read_text(path: &Path) -> Result<String> {
    let path = pathlib_form(path);
    if path.as_os_str().as_bytes().contains(&0) {
        crate::uncaught!("embedded null byte");
    }
    let fail = |detail: String| {
        Error::comparator(format!(
            "cannot read {}: {detail}",
            path_text(path.as_os_str(), usize::MAX)
        ))
    };
    let mut bytes = std::fs::read(&path).map_err(|error| fail(strerror(&error)))?;
    if path.as_os_str().as_bytes().ends_with(b".gz") {
        bytes = gunzip(&bytes).map_err(fail)?;
    }
    match String::from_utf8(bytes) {
        Ok(text) => Ok(universal_newlines(text)),
        Err(error) => Err(fail(
            pyfmt::decode_utf8(error.as_bytes())
                .err()
                .unwrap_or_default(),
        )),
    }
}

/// The bytes of the file at `path` as a stream, inflated where `gzip` is set as `gzip.open`
/// inflates them: member after member, a failure in the interpreter's words
/// (`stream::gzip`). Such a failure is an error whose text is those words, which
/// `report::strerror` gives back as they stand. The text inflated before a failure is
/// read first, as the interpreter hands a reader the text before a member's bad trailer
/// or before the place a file is cut.
pub(crate) fn byte_stream(path: &Path, gzip: bool) -> std::io::Result<Box<dyn Read>> {
    if !gzip {
        return Ok(Box::new(File::open(path)?));
    }
    let path = path.to_path_buf();
    let open: Open = Box::new(move || Ok(Box::new(File::open(&path)?)));
    Ok(Box::new(Gunzipped {
        gunzip: Gunzip::new(open)?,
        failure: None,
    }))
}

/// The text of a gzip file, read as any other stream is.
struct Gunzipped {
    gunzip: Gunzip,
    /// The failure met behind text already handed out.
    failure: Option<std::io::Error>,
}

impl Read for Gunzipped {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        if let Some(failure) = self.failure.take() {
            return Err(failure);
        }
        let mut filled = 0;
        let failure = match self.gunzip.read(out, &mut filled) {
            Ok(_) => return Ok(filled),
            Err(Failure::System(error)) => error,
            Err(Failure::Text(words)) => std::io::Error::other(words),
        };
        if filled == 0 {
            return Err(failure);
        }
        self.failure = Some(failure);
        Ok(filled)
    }
}

/// Whether `path` is a directory. A path the system cannot look up (a name longer than a
/// file name can be) is a malformed input and ends the run with the comparator's error,
/// whatever the interpreter's own `is_dir` makes of it.
pub fn is_directory(path: &Path) -> Result<bool> {
    let path = pathlib_form(path);
    let fail = |detail: &str| {
        Error::comparator(format!(
            "cannot read {}: {detail}",
            path_text(path.as_os_str(), 200)
        ))
    };
    if path.as_os_str().as_bytes().contains(&0) {
        return Err(fail("stat: embedded null character in path"));
    }
    match std::fs::metadata(&path) {
        Ok(found) => Ok(found.is_dir()),
        Err(error)
            if matches!(
                error.kind(),
                std::io::ErrorKind::NotFound | std::io::ErrorKind::NotADirectory
            ) =>
        {
            Ok(false)
        }
        Err(error) => Err(fail(&strerror(&error))),
    }
}

/// The lines of `text` without their terminators, no trailing empty line, borrowed.
pub fn line_slices(text: &str) -> impl Iterator<Item = &str> {
    let body = text.strip_suffix('\n').unwrap_or(text);
    body.split('\n')
        .take(if text.is_empty() { 0 } else { usize::MAX })
        .map(|line| line.strip_suffix('\r').unwrap_or(line))
}

/// The lines of `text` without their terminators, no trailing empty line.
pub fn text_lines(text: &str) -> Vec<String> {
    line_slices(text).map(str::to_string).collect()
}

/// Why the reading of a text's lines ended before its last line.
pub(crate) enum Stopped<E> {
    /// The text could not be read.
    Read(std::io::Error),
    /// What the visitor of a line returned.
    Visit(E),
}

/// What is called with each line of a text: the line without its end, the line with it,
/// and the place of the line's first byte in the text.
pub(crate) type LineVisitor<'a, E> =
    dyn FnMut(&[u8], &[u8], u64) -> std::result::Result<(), E> + 'a;

/// Calls `visit` with each line of the text `reader` gives, as `read_text` and
/// `text_lines` cut a text together: a line feed, a carriage return and the two together
/// each end a line, an empty line is a line, and what follows the last line end is a line
/// where it holds a byte. A line is given without its end, then with it, then with the
/// place of its first byte in the text.
pub(crate) fn text_lines_of<E>(
    reader: &mut dyn Read,
    visit: &mut LineVisitor<'_, E>,
) -> std::result::Result<(), Stopped<E>> {
    let mut buffer = vec![0u8; INFLATE_ROOM];
    let (mut start, mut end, mut ended) = (0, 0, false);
    // The place in the text of the buffer's first byte.
    let mut before = 0u64;
    loop {
        while let Some(found) = memchr::memchr2(b'\n', b'\r', &buffer[start..end]) {
            let at = start + found;
            // A carriage return that ends what is read so far may open a `\r\n`.
            if buffer[at] == b'\r' && at + 1 == end && !ended {
                break;
            }
            let both = buffer[at] == b'\r' && buffer.get(at + 1) == Some(&b'\n') && at + 1 < end;
            let next = at + 1 + usize::from(both);
            visit(
                &buffer[start..at],
                &buffer[start..next],
                before + start as u64,
            )
            .map_err(Stopped::Visit)?;
            start = next;
        }
        if ended {
            if start < end {
                visit(
                    &buffer[start..end],
                    &buffer[start..end],
                    before + start as u64,
                )
                .map_err(Stopped::Visit)?;
            }
            return Ok(());
        }
        if start > 0 {
            buffer.copy_within(start..end, 0);
            before += start as u64;
            end -= start;
            start = 0;
        }
        if end == buffer.len() {
            buffer.resize(buffer.len() * 2, 0);
        }
        match reader.read(&mut buffer[end..]) {
            Ok(0) => ended = true,
            Ok(read) => end += read,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
            Err(error) => return Err(Stopped::Read(error)),
        }
    }
}

/// The data lines of an input VCF, in order, each as far as the comparison reads it
/// (`record_prefix`): the columns after the fifth, which hold most of a line, are not kept.
///
/// The file is read as a stream and fails as `read_text` fails on it: a file that cannot
/// be read, or inflated, to its end is reported as that, and only a file read to its end
/// is reported for its first bytes that are no UTF-8, with their place in the whole text.
pub fn input_records(path: &Path) -> Result<Vec<String>> {
    let path = pathlib_form(path);
    if path.as_os_str().as_bytes().contains(&0) {
        crate::uncaught!("embedded null byte");
    }
    let fail = |detail: String| {
        Error::comparator(format!(
            "cannot read {}: {detail}",
            path_text(path.as_os_str(), usize::MAX)
        ))
    };
    let gzip = path.as_os_str().as_bytes().ends_with(b".gz");
    let mut reader = byte_stream(&path, gzip).map_err(|error| fail(strerror(&error)))?;
    let mut records = Vec::new();
    let mut refused: Option<String> = None;
    let read = text_lines_of(reader.as_mut(), &mut |line, with_end, at| {
        if refused.is_none() {
            match std::str::from_utf8(line) {
                Ok(line) if line.starts_with('#') => {}
                Ok(line) => records.push(record_prefix(line).to_string()),
                Err(_) => refused = utf8::failure(with_end, at),
            }
        }
        Ok::<(), std::convert::Infallible>(())
    });
    match (read, refused) {
        (Err(Stopped::Read(error)), _) => Err(fail(strerror(&error))),
        (_, Some(detail)) => Err(fail(detail)),
        _ => Ok(records),
    }
}

/// Suffixes the second and later occurrences of a key with `#n`.
///
/// A key is counted as it stands before its suffix, so a key that spells an earlier key's
/// suffixed form is a first occurrence and keeps its text.
pub fn dedup_keys(entries: &mut [Entry]) {
    let mut seen: HashMap<String, u64> = HashMap::new();
    for entry in entries {
        match seen.get_mut(&entry.key) {
            Some(count) => {
                *count += 1;
                write!(entry.key, "#{count}").expect("a String takes any text");
            }
            None => {
                seen.insert(entry.key.clone(), 1);
            }
        }
    }
}

/// A JSON object's text with its members sorted (`json_compact`), or `text` as it is
/// when it is not a JSON object. Ensembl VEP encodes the `GA4GH_VRS` object with its keys
/// in hash order, which differs from row to row, so the text formats compare the parsed
/// object, as the JSON format's nested object does. A text nested deeper than the
/// interpreter reads or walks is a malformed input and ends the run with the comparator's
/// error.
pub fn canonical_json_text(text: &str) -> Result<String> {
    canonical_json_text_under(text, CANONICAL_CALLERS)
}

/// `canonical_json_text` called with `callers` frames of the interpreter's call stack
/// below it: the depth an object is walked to is what is left of the interpreter's frames.
pub fn canonical_json_text_under(text: &str, callers: usize) -> Result<String> {
    let too_deep = |what: &str| {
        Error::comparator(format!(
            "a JSON {what} is nested too deep to read: {}...",
            pyfmt::prefix(text, 24)
        ))
    };
    let value = match pyjson::loads(text) {
        Ok(value) => value,
        Err(error) if error.is_recursion_error() => return Err(too_deep("value")),
        Err(error) if error.is_value_error() => return Ok(text.to_string()),
        // The interpreter reads the text; only an object goes on to its canonical form.
        Err(_) if !opens_an_object(text) => return Ok(text.to_string()),
        Err(_) => crate::uncaught!(
            "a JSON object holds half a surrogate pair: {}...",
            pyfmt::prefix(text, 24)
        ),
    };
    if value.as_dict().is_none() {
        return Ok(text.to_string());
    }
    if callers + 1 + pyjson::compact_frames(&value) > pyjson::RECURSION_LIMIT {
        return Err(too_deep("object"));
    }
    Ok(pyjson::dumps_compact(&value))
}

/// `canonical_json_text` with every string member that is a whole number made a number:
/// the `GA4GH_VRS` object of the text formats carries its `start` as a string, and the
/// JSON format's nested object has both bounds as numbers (`OutputFactory/JSON.pm:288`,
/// `numberify`).
///
/// Neither the interpreter's recursion limit nor its limit on the digits of a whole number
/// is the comparator's error: a text past either ends the run. A text nested deeper than
/// the reader reads ends it in the reader's words. The object is then walked member by
/// member, and of a member nested past the recursion limit and a string of more digits
/// than a whole number may have, the one the walk meets first ends the run; an object
/// the walk passes may still be nested too deep to write.
pub fn numberified_json_text(text: &str) -> String {
    let value = match pyjson::loads(text) {
        Ok(value) => value,
        Err(error) if error.is_value_error() => return text.to_string(),
        Err(error) if error.is_recursion_error() => crate::uncaught!("{error}"),
        Err(_) if !opens_an_object(text) => return text.to_string(),
        Err(_) => crate::uncaught!(
            "a JSON object holds half a surrogate pair: {}...",
            pyfmt::prefix(text, 24)
        ),
    };
    if value.as_dict().is_none() {
        return text.to_string();
    }
    // The walk of the object itself is the second call above this function's callers.
    let mut value = value;
    numberify(&mut value, NUMBERIFIED_CALLERS + 2);
    if NUMBERIFIED_CALLERS + 1 + pyjson::compact_frames(&value) > pyjson::RECURSION_LIMIT {
        crate::uncaught!("{RECURSION}");
    }
    pyjson::dumps_compact(&value)
}

/// What `RecursionError` says of a walk or a writing that nests too deep.
const RECURSION: &str = "maximum recursion depth exceeded";

/// Makes every string of `value` that is a whole number a number. The interpreter's walk
/// takes one call a value, the members of an object and the items of a list in their
/// order: `frame` is the depth of the call for `value`, a member's or an item's is one
/// more, and a value whose call would pass the recursion limit ends the run.
fn numberify(value: &mut Value, frame: usize) {
    let mut pending = vec![(value, frame)];
    while let Some((value, frame)) = pending.pop() {
        if frame > pyjson::RECURSION_LIMIT {
            crate::uncaught!("{RECURSION}");
        }
        match value {
            Value::Dict(members) => {
                pending.extend(members.values_mut().rev().map(|member| (member, frame + 1)));
            }
            Value::List(items) => {
                pending.extend(items.iter_mut().rev().map(|item| (item, frame + 1)));
            }
            Value::Str(text) if pyfmt::is_decimal(text) => {
                let whole = match pyfmt::int_checked(text) {
                    Ok(decimal) => Value::of_int(&decimal),
                    Err(error) => crate::uncaught!("{error}"),
                };
                *value = whole;
            }
            _ => {}
        }
    }
}

/// Whether the first character of `text` after the white space of JSON is `{`.
fn opens_an_object(text: &str) -> bool {
    text.trim_start_matches([' ', '\t', '\n', '\r'])
        .starts_with('{')
}

/// `path` as `pathlib.Path` holds it: no empty and no `.` component, no trailing `/`,
/// `.` for the empty path, a root of exactly two `/` kept and of more read as one.
fn pathlib_form(path: &Path) -> PathBuf {
    let bytes = path.as_os_str().as_bytes();
    let slashes = bytes.iter().take_while(|byte| **byte == b'/').count();
    let mut form: Vec<u8> = match slashes {
        0 => Vec::new(),
        2 => b"//".to_vec(),
        _ => b"/".to_vec(),
    };
    let root = form.len();
    for part in bytes.split(|byte| *byte == b'/') {
        if part.is_empty() || part == b"." {
            continue;
        }
        if form.len() > root {
            form.push(b'/');
        }
        form.extend_from_slice(part);
    }
    if form.is_empty() {
        form.push(b'.');
    }
    PathBuf::from(OsStr::from_bytes(&form))
}

/// The first `limit` characters of a path as an error line prints it: a byte that is not
/// UTF-8 is one character, printed as the escape `\udcXX` of the surrogate the interpreter
/// holds it as.
fn path_text(path: &OsStr, limit: usize) -> String {
    let mut text = String::new();
    let mut left = limit;
    for chunk in path.as_bytes().utf8_chunks() {
        for c in chunk.valid().chars() {
            if left == 0 {
                return text;
            }
            text.push(c);
            left -= 1;
        }
        for byte in chunk.invalid() {
            if left == 0 {
                return text;
            }
            write!(text, "\\udc{byte:02x}").expect("a String takes any text");
            left -= 1;
        }
    }
    text
}

/// `text` with every `\r\n` and every other `\r` a `\n`.
pub(crate) fn universal_newlines(text: String) -> String {
    if !text.contains('\r') {
        return text;
    }
    let mut out = String::with_capacity(text.len());
    let mut rest = text.as_str();
    while let Some(at) = rest.find('\r') {
        out.push_str(&rest[..at]);
        out.push('\n');
        rest = &rest[at + 1..];
        rest = rest.strip_prefix('\n').unwrap_or(rest);
    }
    out.push_str(rest);
    out
}

/// `repr()` of a `bytes`.
fn bytes_repr(bytes: &[u8]) -> String {
    let quote = if bytes.contains(&b'\'') && !bytes.contains(&b'"') {
        '"'
    } else {
        '\''
    };
    let mut text = String::from("b");
    text.push(quote);
    for &byte in bytes {
        match byte {
            b'\\' => text.push_str("\\\\"),
            b'\t' => text.push_str("\\t"),
            b'\n' => text.push_str("\\n"),
            b'\r' => text.push_str("\\r"),
            b'\'' if quote == '\'' => text.push_str("\\'"),
            0x20..=0x7e => text.push(char::from(byte)),
            _ => write!(text, "\\x{byte:02x}").expect("a String takes any text"),
        }
    }
    text.push(quote);
    text
}

/// The members of a gzip file, one after another, as `gzip.open` reads them: a zero-byte
/// file is an empty text, and zero bytes after a member are padding, skipped up to the
/// next member or the end. The error is the text of the exception `gzip.open` raises.
pub(crate) fn gunzip(mut rest: &[u8]) -> std::result::Result<Vec<u8>, String> {
    let mut bytes = Vec::new();
    let mut room = vec![0; INFLATE_ROOM];
    while !rest.is_empty() {
        rest = gunzip_member(rest, &mut bytes, &mut room)?;
        while let [0, tail @ ..] = rest {
            rest = tail;
        }
    }
    Ok(bytes)
}

/// The bytes after one member's header, read as `gzip.open` reads it: the header checksum
/// is skipped without a check, a flag bit the format reserves is read past, and a file
/// name or comment may run to the end of the file.
fn after_gzip_header(member: &[u8]) -> std::result::Result<&[u8], String> {
    let cut_short = || GZIP_CUT_SHORT.to_string();
    let magic = &member[..member.len().min(2)];
    if magic != GZIP_MAGIC {
        return Err(format!("Not a gzipped file ({})", bytes_repr(magic)));
    }
    let fixed = member.get(2..10).ok_or_else(cut_short)?;
    if fixed[0] != GZIP_DEFLATE {
        return Err("Unknown compression method".to_string());
    }
    let flags = fixed[1];
    let mut rest = &member[10..];
    if flags & FEXTRA != 0 {
        let length = rest.get(..2).ok_or_else(cut_short)?;
        let length = usize::from(u16::from_le_bytes([length[0], length[1]]));
        rest = rest.get(2 + length..).ok_or_else(cut_short)?;
    }
    for text in [FNAME, FCOMMENT] {
        if flags & text != 0 {
            rest = memchr::memchr(0, rest).map_or(&[][..], |end| &rest[end + 1..]);
        }
    }
    if flags & FHCRC != 0 {
        rest = rest.get(2..).ok_or_else(cut_short)?;
    }
    Ok(rest)
}

/// Inflates one member onto `bytes` and returns what follows its trailer, the checksum and
/// the length of its text checked against the trailer's, the checksum first.
fn gunzip_member<'a>(
    member: &'a [u8],
    bytes: &mut Vec<u8>,
    room: &mut [u8],
) -> std::result::Result<&'a [u8], String> {
    let stream = after_gzip_header(member)?;
    let mut rest = stream;
    let mut inflater = Decompress::new(false);
    let mut crc = Crc::new();
    loop {
        let (read, written) = (inflater.total_in(), inflater.total_out());
        let status = inflater
            .decompress(rest, room, FlushDecompress::None)
            .map_err(|error| {
                let reason = match error.message() {
                    Some(REASON_LOST) | None => refusal_reason(stream, written),
                    Some(reason) => reason.to_string(),
                };
                format!("Error -3 while decompressing data: {reason}")
            })?;
        let text = &room[..(inflater.total_out() - written) as usize];
        crc.update(text);
        bytes.extend_from_slice(text);
        rest = &rest[(inflater.total_in() - read) as usize..];
        if status == Status::StreamEnd {
            break;
        }
        if text.is_empty() && inflater.total_in() == read {
            return Err(GZIP_CUT_SHORT.to_string());
        }
    }
    let trailer = rest.get(..8).ok_or_else(|| GZIP_CUT_SHORT.to_string())?;
    let stored = |at: usize| {
        u32::from_le_bytes([
            trailer[at],
            trailer[at + 1],
            trailer[at + 2],
            trailer[at + 3],
        ])
    };
    if stored(0) != crc.sum() {
        return Err(format!(
            "CRC check failed {:#x} != {:#x}",
            stored(0),
            crc.sum()
        ));
    }
    // The stored length is the text's length modulo 2^32.
    if stored(4) != inflater.total_out() as u32 {
        return Err("Incorrect length of data produced".to_string());
    }
    Ok(&rest[8..])
}

/// zlib's reason for refusing the deflate stream `stream`, whose first `valid` bytes of
/// text an inflater wrote before the call that refused it.
///
/// The stream is inflated again without the fast path ever reaching the error: up to the
/// valid text the room ends with that text, and the fast path stops 260 bytes before the
/// end of its room; past it the room is smaller than the fast path takes.
pub(crate) fn refusal_reason(stream: &[u8], valid: u64) -> String {
    const NO_REASON: &str = "invalid input data";
    let mut inflater = Decompress::new(false);
    let mut room = vec![0; INFLATE_ROOM];
    let mut rest = stream;
    loop {
        let (read, written) = (inflater.total_in(), inflater.total_out());
        let size = match valid.saturating_sub(written) {
            0 => SLOW_INFLATE_ROOM,
            left => left.min(INFLATE_ROOM as u64) as usize,
        };
        match inflater.decompress(rest, &mut room[..size], FlushDecompress::None) {
            Err(error) => return error.message().unwrap_or(NO_REASON).to_string(),
            Ok(Status::StreamEnd) => return NO_REASON.to_string(),
            Ok(_) if inflater.total_in() == read && inflater.total_out() == written => {
                return NO_REASON.to_string()
            }
            Ok(_) => rest = &rest[(inflater.total_in() - read) as usize..],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Uncaught;
    use std::io::Write as _;

    const TEXT: &[u8] = b"#h\nrow\n";

    fn deflated(text: &[u8]) -> Vec<u8> {
        let mut deflate =
            flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
        deflate.write_all(text).unwrap();
        deflate.finish().unwrap()
    }

    /// One gzip member whose deflate stream is `stream` and whose trailer is that of
    /// `text`, under a header of `flags`, `fields` standing after the fixed header.
    fn member_of(stream: &[u8], text: &[u8], flags: u8, fields: &[u8]) -> Vec<u8> {
        let mut crc = Crc::new();
        crc.update(text);
        let mut member = vec![0x1f, 0x8b, 8, flags, 0, 0, 0, 0, 0, 0xff];
        member.extend_from_slice(fields);
        member.extend_from_slice(stream);
        member.extend_from_slice(&crc.sum().to_le_bytes());
        member.extend_from_slice(&(text.len() as u32).to_le_bytes());
        member
    }

    fn member(text: &[u8]) -> Vec<u8> {
        member_of(&deflated(text), text, 0, b"")
    }

    fn refusal(bytes: &[u8]) -> String {
        gunzip(bytes).expect_err("the file is refused")
    }

    /// A fresh directory of this test's own.
    fn scratch(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "vep-formats-compare-input-{}-{name}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    fn error_text<T: std::fmt::Debug>(result: Result<T>) -> String {
        match result.expect_err("the comparator's error") {
            Error::Comparator(text) => text,
        }
    }

    /// What a run ends with on an input no reading of the comparison expects.
    fn uncaught(run: impl FnOnce() -> String + std::panic::UnwindSafe) -> String {
        let payload = std::panic::catch_unwind(run).expect_err("the run ends");
        let Uncaught(what) = payload
            .downcast_ref::<Uncaught>()
            .expect("an uncaught exception");
        what.clone()
    }

    fn nested(depth: usize, leaf: &str) -> String {
        format!("{}{leaf}{}", "{\"a\":".repeat(depth), "}".repeat(depth))
    }

    #[test]
    fn a_line_ends_at_a_line_feed_and_loses_one_carriage_return_before_it() {
        let cases: [(&str, &[&str]); 12] = [
            ("", &[]),
            ("\n", &[""]),
            ("a", &["a"]),
            ("a\n", &["a"]),
            ("a\n\n", &["a", ""]),
            ("\n\na", &["", "", "a"]),
            ("a\r\nb\r", &["a", "b"]),
            ("a\r\r\n", &["a\r"]),
            ("a\rb\n", &["a\rb"]),
            ("\r", &[""]),
            ("a\n\r", &["a", ""]),
            (
                "a\u{b}b\u{c}c\u{85}d\u{2028}e\n",
                &["a\u{b}b\u{c}c\u{85}d\u{2028}e"],
            ),
        ];
        for (text, lines) in cases {
            assert_eq!(text_lines(text), lines, "{text:?}");
            assert_eq!(line_slices(text).collect::<Vec<_>>(), lines, "{text:?}");
        }
    }

    #[test]
    fn text_mode_reads_every_carriage_return_as_a_line_end() {
        let cases = [
            ("#h\r\nrow\r\n", "#h\nrow\n"),
            ("#h\rrow\rlast", "#h\nrow\nlast"),
            ("a\r\r\nb\n\rc\r", "a\n\nb\n\nc\n"),
            ("\r", "\n"),
            ("no carriage return\n", "no carriage return\n"),
            (
                "a\u{b}b\u{c}c\u{1c}d\u{85}e\u{2028}f\u{2029}g\n",
                "a\u{b}b\u{c}c\u{1c}d\u{85}e\u{2028}f\u{2029}g\n",
            ),
        ];
        for (text, read) in cases {
            assert_eq!(universal_newlines(text.to_string()), read, "{text:?}");
        }
    }

    #[test]
    fn a_key_is_counted_as_it_stands_before_its_suffix() {
        let mut entries: Vec<Entry> = ["k", "k", "k#2", "k", "k#2"]
            .iter()
            .map(|key| Entry {
                key: key.to_string(),
                ..Entry::default()
            })
            .collect();
        dedup_keys(&mut entries);
        let keys: Vec<&str> = entries.iter().map(|e| e.key.as_str()).collect();
        assert_eq!(keys, ["k", "k#2", "k#2", "k#3", "k#2#2"]);
        dedup_keys(&mut []);
    }

    #[test]
    fn a_path_is_held_as_pathlib_holds_it() {
        let cases = [
            ("", "."),
            (".", "."),
            ("./", "."),
            ("a", "a"),
            ("./a", "a"),
            ("a/", "a"),
            ("a//b/./c/", "a/b/c"),
            ("a/..", "a/.."),
            ("../a", "../a"),
            ("/", "/"),
            ("/a", "/a"),
            ("//", "//"),
            ("//a", "//a"),
            ("///a//b", "/a/b"),
            ("/./a", "/a"),
            ("a.gz/", "a.gz"),
        ];
        for (given, held) in cases {
            assert_eq!(pathlib_form(Path::new(given)), Path::new(held), "{given:?}");
        }
    }

    #[test]
    fn a_path_byte_that_is_not_utf8_prints_as_its_surrogate_escape() {
        let name = OsStr::from_bytes(b"w\xff\xc3\xa4\xe4\xb8.txt");
        assert_eq!(
            path_text(name, usize::MAX),
            "w\\udcff\u{e4}\\udce4\\udcb8.txt"
        );
        // Each such byte is one character of the 200 an error line prints.
        assert_eq!(path_text(name, 3), "w\\udcff\u{e4}");
        assert_eq!(path_text(name, 0), "");
        assert_eq!(path_text(OsStr::new("plain"), 200), "plain");
    }

    #[test]
    fn bytes_print_as_the_interpreter_prints_them() {
        assert_eq!(bytes_repr(b"#h"), "b'#h'");
        assert_eq!(bytes_repr(b"\x1f"), "b'\\x1f'");
        assert_eq!(bytes_repr(b"\x00\xff"), "b'\\x00\\xff'");
        assert_eq!(bytes_repr(b"'"), "b\"'\"");
        assert_eq!(bytes_repr(b"'\""), "b'\\'\"'");
        assert_eq!(bytes_repr(b"\"x"), "b'\"x'");
        assert_eq!(bytes_repr(b"\\\t"), "b'\\\\\\t'");
        assert_eq!(bytes_repr(b"\n\r"), "b'\\n\\r'");
        assert_eq!(bytes_repr(b" ~"), "b' ~'");
        assert_eq!(bytes_repr(b"\x7f"), "b'\\x7f'");
        assert_eq!(bytes_repr(b""), "b''");
    }

    #[test]
    fn gzip_members_are_read_one_after_another_past_their_padding() {
        assert_eq!(gunzip(b"").unwrap(), b"");
        assert_eq!(gunzip(&member(TEXT)).unwrap(), TEXT);
        assert_eq!(gunzip(&member(b"")).unwrap(), b"");
        let mut two = member(b"#h\n");
        two.extend(member(b"row\n"));
        assert_eq!(gunzip(&two).unwrap(), TEXT);
        let mut padded = member(b"#h\n");
        padded.extend_from_slice(&[0; 300]);
        padded.extend(member(b"row\n"));
        padded.extend_from_slice(&[0; 7]);
        assert_eq!(gunzip(&padded).unwrap(), TEXT);
    }

    #[test]
    fn a_gzip_header_is_read_past_without_a_check_of_its_checksum_or_reserved_flags() {
        let stream = deflated(TEXT);
        // A header checksum that is not the header's.
        let checked = member_of(&stream, TEXT, FHCRC, &[0x12, 0x34]);
        assert_eq!(gunzip(&checked).unwrap(), TEXT);
        // The three flag bits the format reserves.
        for reserved in [0x20, 0x40, 0x80, 0xe0] {
            assert_eq!(
                gunzip(&member_of(&stream, TEXT, reserved, b"")).unwrap(),
                TEXT
            );
        }
        // An extra field, a file name, a comment and a header checksum, in that order.
        let all = FEXTRA | FNAME | FCOMMENT | FHCRC;
        let fields = b"\x03\x00abcname\0a note\0\xff\xff";
        assert_eq!(
            gunzip(&member_of(&stream, TEXT, all, fields)).unwrap(),
            TEXT
        );
    }

    #[test]
    fn a_file_that_is_no_gzip_file_is_refused_with_its_first_two_bytes() {
        let whole = member(TEXT);
        let after = |tail: &[u8]| [whole.as_slice(), tail].concat();
        let cases: [(Vec<u8>, &str); 9] = [
            (TEXT.to_vec(), "b'#h'"),
            (vec![0x1f], "b'\\x1f'"),
            (b"'".to_vec(), "b\"'\""),
            (b"'\"rest".to_vec(), "b'\\'\"'"),
            (b"\\\tx".to_vec(), "b'\\\\\\t'"),
            (vec![0; 16], "b'\\x00\\x00'"),
            (after(b"not a member"), "b'no'"),
            (after(b"x"), "b'x'"),
            (after(b"\0\0\0zz"), "b'zz'"),
        ];
        for (bytes, magic) in cases {
            assert_eq!(refusal(&bytes), format!("Not a gzipped file ({magic})"));
        }
    }

    #[test]
    fn a_gzip_file_cut_short_is_refused_wherever_it_is_cut() {
        let whole = member(TEXT);
        let last = whole.len() - 1;
        let mut cut: Vec<Vec<u8>> = [2, 9, 10, 12, last, last - 7]
            .iter()
            .map(|at| whole[..*at].to_vec())
            .collect();
        // A file name that never ends, an extra field longer than the file and a header
        // checksum of one byte.
        cut.push(b"\x1f\x8b\x08\x08\0\0\0\0\0\xffa name without its end".to_vec());
        cut.push(member_of(&deflated(b""), b"", FEXTRA, b"\x09\x00abc")[..15].to_vec());
        cut.push(b"\x1f\x8b\x08\x02\0\0\0\0\0\xff\x01".to_vec());
        for bytes in cut {
            assert_eq!(refusal(&bytes), GZIP_CUT_SHORT, "{bytes:?}");
        }
    }

    #[test]
    fn a_gzip_trailer_that_is_not_its_texts_is_refused_the_checksum_first() {
        let whole = member(TEXT);
        let last = whole.len() - 1;
        let mut other_method = whole.clone();
        other_method[2] = 7;
        assert_eq!(refusal(&other_method), "Unknown compression method");
        let mut checksum = whole.clone();
        checksum[last - 7] ^= 1;
        assert_eq!(
            refusal(&checksum),
            "CRC check failed 0xcf388bcd != 0xcf388bcc"
        );
        let mut length = whole.clone();
        length[last] ^= 1;
        assert_eq!(refusal(&length), "Incorrect length of data produced");
        let mut both = checksum.clone();
        both[last] ^= 1;
        assert_eq!(refusal(&both), "CRC check failed 0xcf388bcd != 0xcf388bcc");
        let empty = member(b"");
        let mut small = empty[..empty.len() - 8].to_vec();
        small.extend_from_slice(&[5, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(refusal(&small), "CRC check failed 0x5 != 0x0");
    }

    #[test]
    fn a_deflate_stream_zlib_refuses_is_refused_with_its_reason() {
        let cases: [(&[u8], &str); 9] = [
            (b"\x07\x00", "invalid block type"),
            (b"\x01\x05\x00\x00\x00abcde", "invalid stored block lengths"),
            (
                b"\xfd\xff\xff\xff\xff\xff\xff\xff",
                "too many length or distance symbols",
            ),
            (
                b"\x05\xc0\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff",
                "invalid code lengths set",
            ),
            (b"\x03\x1f\xff\xff\xff", "invalid distance too far back"),
            (b"\x00\x03\x00\xfc\xffabc\x07\x00", "invalid block type"),
            // Streams long enough for the inflater's fast path, which finds the error.
            (
                b"\x1b\x03\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
                "invalid literal/length code",
            ),
            (
                b"\x03\x1f\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff",
                "invalid distance too far back",
            ),
            (
                b"\x4b\x04\x3e\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
                "invalid distance code",
            ),
        ];
        for (stream, reason) in cases {
            assert_eq!(
                refusal(&member_of(stream, b"", 0, b"")),
                format!("Error -3 while decompressing data: {reason}"),
                "{stream:?}"
            );
        }
    }

    #[test]
    fn a_deflate_error_past_the_first_room_keeps_its_reason() {
        // A stream of more text than one room, flushed to a byte, then a block of the
        // fixed code whose first symbol is no literal and no length.
        let text: Vec<u8> = (0..3 * INFLATE_ROOM / 2)
            .map(|at| (at % 251) as u8)
            .collect();
        let mut deflate =
            flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::fast());
        deflate.write_all(&text).unwrap();
        deflate.flush().unwrap();
        let mut stream = deflate.get_ref().clone();
        stream.extend_from_slice(b"\x1b\x03\0\0\0\0\0\0\0\0\0\0\0\0\0\0");
        assert_eq!(
            refusal(&member_of(&stream, &text, 0, b"")),
            "Error -3 while decompressing data: invalid literal/length code"
        );
    }

    #[test]
    fn read_text_reads_plain_and_gzip_files_in_text_mode() {
        let directory = scratch("read");
        let write = |name: &str, bytes: &[u8]| {
            let path = directory.join(name);
            std::fs::write(&path, bytes).unwrap();
            path
        };
        assert_eq!(
            read_text(&write("plain.txt", b"#h\r\nrow\rlast")).unwrap(),
            "#h\nrow\nlast"
        );
        assert_eq!(read_text(&write("empty.txt", b"")).unwrap(), "");
        assert_eq!(
            read_text(&write("bom.txt", b"\xef\xbb\xbf#h\n")).unwrap(),
            "\u{feff}#h\n"
        );
        assert_eq!(
            read_text(&write("crlf.txt.gz", &member(b"#h\r\nrow\rlast"))).unwrap(),
            "#h\nrow\nlast"
        );
        assert_eq!(read_text(&write("empty.txt.gz", b"")).unwrap(), "");
        // A file is gunzipped for its name alone: a member under another name is read as
        // it stands, and is no text.
        let named = write("member.txt", &member(TEXT));
        assert_eq!(
            error_text(read_text(&named)),
            format!(
                "cannot read {}: 'utf-8' codec can't decode byte 0x8b in position 1: invalid start byte",
                named.display()
            )
        );
        let records = write(
            "input.vcf",
            b"##fileformat=VCFv4.2\n#CHROM\tPOS\n21\t5\r\n\n21\t7",
        );
        assert_eq!(input_records(&records).unwrap(), ["21\t5", "", "21\t7"]);
        std::fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn a_file_that_cannot_be_read_is_an_error_in_the_interpreters_words() {
        let directory = scratch("refused");
        let shown = |name: &str| directory.join(name).display().to_string();
        let write = |name: &str, bytes: &[u8]| {
            let path = directory.join(name);
            std::fs::write(&path, bytes).unwrap();
            path
        };
        let cases: [(&str, &[u8], &str); 8] = [
            (
                "start.txt",
                b"ab\xffcd\n",
                "'utf-8' codec can't decode byte 0xff in position 2: invalid start byte",
            ),
            (
                "continuation.txt",
                b"ab\xe4\xb8xcd\n",
                "'utf-8' codec can't decode bytes in position 2-3: invalid continuation byte",
            ),
            (
                "end.txt",
                b"ab\xe4\xb8",
                "'utf-8' codec can't decode bytes in position 2-3: unexpected end of data",
            ),
            ("not_gzip.txt.gz", TEXT, "Not a gzipped file (b'#h')"),
            ("cut.txt.gz", &member(TEXT)[..12], GZIP_CUT_SHORT),
            (
                "bad_utf8.txt.gz",
                &member(b"ab\xffcd\n"),
                "'utf-8' codec can't decode byte 0xff in position 2: invalid start byte",
            ),
            // The position counts the bytes of every member before it.
            (
                "second_member.txt.gz",
                &[member(b"abc\n"), member(b"\xff\n")].concat(),
                "'utf-8' codec can't decode byte 0xff in position 4: invalid start byte",
            ),
            // A file both cut short and not UTF-8 is cut short: it is inflated whole first.
            (
                "bad_utf8_and_cut.txt.gz",
                &member(b"ab\xffcd\n")[..14],
                GZIP_CUT_SHORT,
            ),
        ];
        for (name, bytes, detail) in cases {
            let path = write(name, bytes);
            assert_eq!(
                error_text(read_text(&path)),
                format!("cannot read {}: {detail}", shown(name))
            );
        }
        for name in ["absent.txt", "absent.txt.gz"] {
            assert_eq!(
                error_text(read_text(&directory.join(name))),
                format!("cannot read {}: No such file or directory", shown(name))
            );
        }
        for name in ["a_directory", "a_directory.gz"] {
            std::fs::create_dir(directory.join(name)).unwrap();
            assert_eq!(
                error_text(read_text(&directory.join(name))),
                format!("cannot read {}: Is a directory", shown(name))
            );
        }
        assert_eq!(
            error_text(read_text(&directory.join("start.txt/below"))),
            format!("cannot read {}: Not a directory", shown("start.txt/below"))
        );
        let long = "x".repeat(300);
        assert_eq!(
            error_text(read_text(&directory.join(&long))),
            format!("cannot read {}: File name too long", shown(&long))
        );
        // The path is printed as pathlib holds it, and that path is the one opened.
        let spelt = format!("{}//./absent.txt/", directory.display());
        assert_eq!(
            error_text(read_text(Path::new(&spelt))),
            format!(
                "cannot read {}: No such file or directory",
                shown("absent.txt")
            )
        );
        let spelt = format!("{}//./end.txt/", directory.display());
        assert_eq!(
            error_text(read_text(Path::new(&spelt))),
            format!(
                "cannot read {}: 'utf-8' codec can't decode bytes in position 2-3: unexpected end of data",
                shown("end.txt")
            )
        );
        std::fs::remove_dir_all(&directory).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn a_file_without_read_permission_is_refused_as_the_system_words_it() {
        use std::os::unix::fs::PermissionsExt;
        let directory = scratch("permission");
        let path = directory.join("closed.txt");
        std::fs::write(&path, b"x\n").unwrap();
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o000)).unwrap();
        // A process allowed to read every file meets no refusal.
        if std::fs::read(&path).is_err() {
            assert_eq!(
                error_text(read_text(&path)),
                format!("cannot read {}: Permission denied", path.display())
            );
        }
        std::fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn is_directory_tells_a_directory_and_refuses_a_path_the_system_cannot_look_up() {
        let directory = scratch("directory");
        std::fs::write(directory.join("file"), b"").unwrap();
        assert_eq!(is_directory(&directory), Ok(true));
        assert_eq!(is_directory(&directory.join("file")), Ok(false));
        assert_eq!(is_directory(&directory.join("absent")), Ok(false));
        assert_eq!(is_directory(&directory.join("file/below")), Ok(false));
        assert_eq!(is_directory(&directory.join(".")), Ok(true));
        let long = directory.join("x".repeat(300));
        let shown: String = long.display().to_string().chars().take(200).collect();
        assert_eq!(
            error_text(is_directory(&long)),
            format!("cannot read {shown}: File name too long")
        );
        let nul = directory.join("a\0b");
        assert_eq!(
            error_text(is_directory(&nul)),
            format!(
                "cannot read {}: stat: embedded null character in path",
                nul.display()
            )
        );
        std::fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn a_json_object_is_written_with_its_members_sorted_and_any_other_text_as_it_is() {
        let perl = concat!(
            r#"{"type":"Allele","state":{"sequence":"C","type":"SequenceState"},"#,
            r#""location":{"type":"SequenceLocation","sequence_id":"refseq:NC_000021.9","#,
            r#""interval":{"type":"SimpleInterval","end":10538885,"start":"10538884"}}}"#
        );
        let rust = concat!(
            r#"{"location":{"type":"SequenceLocation","interval":{"type":"SimpleInterval","#,
            r#""start":"10538884","end":10538885},"sequence_id":"refseq:NC_000021.9"},"#,
            r#""state":{"type":"SequenceState","sequence":"C"},"type":"Allele"}"#
        );
        let canonical = concat!(
            r#"{"location":{"interval":{"end":10538885,"start":"10538884","type":"SimpleInterval"},"#,
            r#""sequence_id":"refseq:NC_000021.9","type":"SequenceLocation"},"#,
            r#""state":{"sequence":"C","type":"SequenceState"},"type":"Allele"}"#
        );
        assert_eq!(canonical_json_text(perl).unwrap(), canonical);
        assert_eq!(canonical_json_text(rust).unwrap(), canonical);
        for text in [
            "21:10538884:A:C",
            "-",
            "1e5",
            "",
            "[{\"b\":1,\"a\":2}]",
            "{\"b\":1,\"a\":2} x",
            "\"x\"",
        ] {
            assert_eq!(canonical_json_text(text).unwrap(), text);
        }
        assert_eq!(
            canonical_json_text(" {\"b\": 1, \"a\": 2}\n").unwrap(),
            r#"{"a":2,"b":1}"#
        );
        // The number form of the same object: the start a number, the end as it was.
        let numbers = canonical.replace(r#""start":"10538884""#, r#""start":10538884"#);
        assert_eq!(numberified_json_text(perl), numbers);
        assert_eq!(numberified_json_text("21:10538884:A:C"), "21:10538884:A:C");
        // A superscript digit is a digit and no decimal: it stays a string.
        assert_eq!(
            numberified_json_text("{\"start\":\"\u{b2}\"}"),
            "{\"start\":\"\u{b2}\"}"
        );
        assert_eq!(
            numberified_json_text(
                r#"{"b":["007",{"c":"٣"},"1.0","x"],"a":"99999999999999999999"}"#
            ),
            r#"{"a":99999999999999999999,"b":[7,{"c":3},"1.0","x"]}"#
        );
    }

    #[test]
    fn an_object_nested_past_the_interpreters_frames_is_refused() {
        let refused =
            |text: &str, callers: usize| error_text(canonical_json_text_under(text, callers));
        let object =
            "a JSON object is nested too deep to read: {\"a\":{\"a\":{\"a\":{\"a\":{\"a\"...";
        // A level of nesting takes two frames and the innermost value its own count.
        assert_eq!(
            canonical_json_text(&nested(495, "1")).unwrap(),
            nested(495, "1")
        );
        assert_eq!(error_text(canonical_json_text(&nested(496, "1"))), object);
        assert_eq!(
            canonical_json_text(&nested(494, "null")).unwrap(),
            nested(494, "null")
        );
        assert_eq!(
            error_text(canonical_json_text(&nested(495, "null"))),
            object
        );
        // One frame more below the call, as under `parse_tab`, is half a level less.
        let text = nested(495, "\"x\"");
        assert_eq!(
            canonical_json_text_under(&text, CANONICAL_CALLERS).unwrap(),
            text
        );
        assert_eq!(refused(&text, CANONICAL_CALLERS + 1), object);
        // Lists nested deeper than the interpreter reads are no value at all.
        let lists = format!("{}{}", "[".repeat(100_000), "]".repeat(100_000));
        assert_eq!(
            error_text(canonical_json_text(&lists)),
            "a JSON value is nested too deep to read: [[[[[[[[[[[[[[[[[[[[[[[[..."
        );
    }

    #[test]
    fn the_number_form_ends_the_run_past_the_recursion_and_the_digit_limits() {
        const DIGITS: &str = "Exceeds the limit (4300 digits) for integer string conversion: \
                              value has 4301 digits; use sys.set_int_max_str_digits() to \
                              increase the limit";
        assert_eq!(numberified_json_text(&nested(495, "1")), nested(495, "1"));
        assert_eq!(
            uncaught(|| numberified_json_text(&nested(496, "1"))),
            RECURSION
        );
        assert_eq!(
            uncaught(|| numberified_json_text(&nested(2000, "\"7\""))),
            RECURSION
        );
        // A text nested deeper than the reader reads ends the run in the reader's words.
        let lists = |depth: usize, inner: &str| {
            format!("{}{inner}{}", "[".repeat(depth), "]".repeat(depth))
        };
        assert_eq!(
            uncaught(|| numberified_json_text(&lists(100_000, ""))),
            "Stack overflow (used 10200 kB) while decoding a JSON array from a unicode string"
        );
        // A whole number of more digits than the interpreter converts.
        let digits = |count: usize| format!("{{\"start\":\"{}\"}}", "7".repeat(count));
        assert_eq!(
            numberified_json_text(&digits(4300)),
            digits(4300).replace('"', "").replace("start", "\"start\"")
        );
        assert_eq!(uncaught(|| numberified_json_text(&digits(4301))), DIGITS);
        // Of a member nested past the recursion limit and a string of too many digits, the
        // one the walk meets first ends the run.
        let long = format!("\"{}\"", "7".repeat(4301));
        let members = |first: &str, second: &str| format!("{{\"a\":{first},\"b\":{second}}}");
        let ended = |text: String| uncaught(|| numberified_json_text(&text));
        assert_eq!(ended(members(&long, &lists(1100, "1"))), DIGITS);
        assert_eq!(ended(members(&lists(1100, "1"), &long)), RECURSION);
        // The walk reaches a text inside 991 lists and an empty list inside 992.
        for (depth, inner, words) in [
            (991, "\"x\"", DIGITS),
            (992, "\"x\"", RECURSION),
            (992, "", DIGITS),
            (993, "", RECURSION),
        ] {
            assert_eq!(
                ended(members(&lists(depth, inner), &long)),
                words,
                "{depth} lists around {inner:?}"
            );
        }
    }

    #[test]
    fn half_a_surrogate_pair_ends_the_run_only_inside_an_object() {
        for text in ["\"\\ud83d\"", "[\"\\ud83d\"]", " [{\"a\":\"\\ud83d\"}]"] {
            assert_eq!(canonical_json_text(text).unwrap(), text);
            assert_eq!(numberified_json_text(text), text);
        }
        let object = "{\"a\":\"\\ud83d\"}";
        let what = "a JSON object holds half a surrogate pair: {\"a\":\"\\ud83d\"}...";
        assert_eq!(uncaught(|| canonical_json_text(object).unwrap()), what);
        assert_eq!(uncaught(|| numberified_json_text(object)), what);
        // A whole pair is one character.
        assert_eq!(
            canonical_json_text("{\"a\":\"\\ud83d\\ude00\"}").unwrap(),
            "{\"a\":\"\u{1f600}\"}"
        );
    }

    /// What `read_text` and `text_lines` give of a file's data lines, each as far as the
    /// comparison reads it; else the error.
    fn records_of_the_text(path: &Path) -> Result<Vec<String>> {
        let text = read_text(path)?;
        Ok(line_slices(&text)
            .filter(|line| !line.starts_with('#'))
            .map(|line| record_prefix(line).to_string())
            .collect())
    }

    #[test]
    fn the_input_is_read_as_a_stream_into_what_the_comparison_reads_of_each_line() {
        let directory = scratch("records");
        let write = |name: &str, bytes: &[u8]| {
            let path = directory.join(name);
            std::fs::write(&path, bytes).unwrap();
            path
        };
        let text: &[u8] =
            b"##fileformat=VCFv4.2\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
            21\t5\trs1\tA\tT\t.\tPASS\tAC=1;AF=0.5\r\n\
            \n\
            21\t7\t\t8\tA/T\t+\t.\tAC=2\r\
            #not a record\n\
            21 9 9 A/T - id\n\
            \xef\xbb\xbf21\t11\t.\tA\tT\t.\t.\tx\n\
            last\twithout\tits\tend";
        let kept = [
            "21\t5\trs1\tA\tT",
            "",
            "21\t7\t\t8\tA/T\t+",
            "21 9 9 A/T - id",
            "\u{feff}21\t11\t.\tA\tT",
            "last\twithout\tits\tend",
        ];
        for (name, bytes) in [
            ("input.vcf", text.to_vec()),
            ("input.vcf.gz", member(text)),
            (
                "two.vcf.gz",
                [member(&text[..40]), member(&text[40..])].concat(),
            ),
        ] {
            let path = write(name, &bytes);
            assert_eq!(input_records(&path).unwrap(), kept, "{name}");
            assert_eq!(records_of_the_text(&path).unwrap(), kept, "{name}");
        }
        // A long line, and a carriage return and its line feed in two reads of the file.
        let long = [vec![b'x'; INFLATE_ROOM - 1], b"\r\n21\t5\n".to_vec()].concat();
        let path = write("long.vcf", &long);
        assert_eq!(
            input_records(&path).unwrap(),
            records_of_the_text(&path).unwrap()
        );
        assert_eq!(input_records(&path).unwrap().len(), 2);
        assert!(input_records(&write("empty.vcf", b"")).unwrap().is_empty());
        assert!(input_records(&write("empty.vcf.gz", b""))
            .unwrap()
            .is_empty());
        assert_eq!(input_records(&write("end.vcf", b"\n")).unwrap(), [""]);
        std::fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn an_input_that_cannot_be_read_fails_as_the_whole_text_fails() {
        let directory = scratch("records-refused");
        let write = |name: &str, bytes: &[u8]| {
            let path = directory.join(name);
            std::fs::write(&path, bytes).unwrap();
            path
        };
        let good = b"21\t5\t.\tA\tT\n";
        let whole = member(b"#h\n21\t5\t\xff\tA\tT\n21\t7\n");
        let mut wrong_sum = whole.clone();
        let at = wrong_sum.len() - 8;
        wrong_sum[at] ^= 0xff;
        let cases: Vec<(&str, Vec<u8>)> = vec![
            // Bytes that are no UTF-8: the first of them, with its place in the whole text.
            ("start.vcf", b"#h\n21\t5\t\xff\tA\tT\n21\t\xfe\n".to_vec()),
            ("continuation.vcf", b"21\t5\t\xe2\x82\n21\t7\n".to_vec()),
            (
                "before-a-carriage-return.vcf",
                b"21\t5\t\xe2\r\n21\t7\n".to_vec(),
            ),
            ("end.vcf", b"21\t5\n21\t\xe2\x82".to_vec()),
            ("range.vcf", b"21\t5\n\xf0\x9f\x98\n".to_vec()),
            ("in-a-header.vcf", b"#\xc0h\n21\t5\n".to_vec()),
            ("not-utf8.vcf.gz", whole.clone()),
            // A file that is not read, or not inflated, to its end is reported as that,
            // whatever its text held before.
            ("wrong-sum.vcf.gz", wrong_sum),
            ("cut.vcf.gz", whole[..whole.len() - 3].to_vec()),
            ("after.vcf.gz", [whole.as_slice(), b"no"].concat()),
            ("plain.vcf.gz", good.to_vec()),
        ];
        for (name, bytes) in &cases {
            let path = write(name, bytes);
            let refused = error_text(input_records(&path));
            assert_eq!(refused, error_text(read_text(&path)), "{name}");
            assert!(refused.starts_with("cannot read "), "{name}: {refused}");
        }
        assert_eq!(
            error_text(input_records(&write("place.vcf", b"#h\n21\t5\t\xff\tA\tT\n"))),
            format!(
                "cannot read {}: 'utf-8' codec can't decode byte 0xff in position 8: invalid start byte",
                directory.join("place.vcf").display()
            )
        );
        for name in ["absent.vcf", "absent.vcf.gz"] {
            let path = directory.join(name);
            assert_eq!(
                error_text(input_records(&path)),
                error_text(read_text(&path)),
                "{name}"
            );
        }
        assert_eq!(
            error_text(input_records(&directory)),
            error_text(read_text(&directory))
        );
        assert_eq!(
            uncaught(|| error_text(input_records(Path::new("a\0b.vcf")))),
            "embedded null byte"
        );
        std::fs::remove_dir_all(&directory).unwrap();
    }
}
