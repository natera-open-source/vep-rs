// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One output read whole and cut into the lines the comparator's `read_text` and
//! `text_lines` give: UTF-8, gunzipped when the name ends in `.gz`, a line ended by `\n`,
//! `\r\n` or a lone `\r`, a last line without its end kept.

use flate2::{Crc, Decompress, FlushDecompress, Status};
use memchr::{memchr, memchr2};
use std::fs;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

/// A line that starts with `#`.
pub const HEADER: u8 = 1;
/// A header line that does not start with `##`: the column names.
pub const COLUMNS: u8 = 2;

#[derive(Clone, Copy)]
pub struct Line {
    pub start: u64,
    pub len: u32,
    pub flags: u8,
}

pub struct Side {
    pub bytes: Vec<u8>,
    pub lines: Vec<Line>,
}

impl Side {
    pub fn read(path: &Path, threads: usize) -> Result<Side, String> {
        let bytes = read_bytes(path)?;
        let lines = split_lines(&bytes, threads).map_err(|e| format!("{}: {e}", path.display()))?;
        Ok(Side { bytes, lines })
    }

    pub fn text(&self, line: &Line) -> &[u8] {
        let start = line.start as usize;
        &self.bytes[start..start + line.len as usize]
    }
}

fn read_bytes(path: &Path) -> Result<Vec<u8>, String> {
    let fail = |e: std::io::Error| format!("cannot read {}: {e}", path.display());
    if !path.as_os_str().as_bytes().ends_with(b".gz") {
        return fs::read(path).map_err(fail);
    }
    gunzip(&fs::read(path).map_err(fail)?).map_err(fail)
}

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const GZIP_DEFLATE: u8 = 8;
const FHCRC: u8 = 2;
const FEXTRA: u8 = 4;
const FNAME: u8 = 8;
const FCOMMENT: u8 = 16;
/// The room kept free for the inflater's next write.
const INFLATE_ROOM: usize = 1 << 20;

fn gzip_error(detail: &str) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, detail.to_string())
}

fn gzip_cut_short() -> std::io::Error {
    gzip_error("the gzip file ends inside a member")
}

/// The members of a gzip file, one after another, as `gzip.open` reads them: a zero-byte
/// file is an empty text, and zero bytes after a member are padding, skipped up to the
/// next member or the end.
fn gunzip(mut rest: &[u8]) -> std::io::Result<Vec<u8>> {
    let mut bytes = Vec::new();
    while !rest.is_empty() {
        rest = gunzip_member(rest, &mut bytes)?;
        while let [0, tail @ ..] = rest {
            rest = tail;
        }
    }
    Ok(bytes)
}

/// The bytes after one member's header, read as `gzip.open` reads it: the header checksum
/// is skipped without a check, a flag bit the format reserves is read past, and a file
/// name or comment may run to the end of the file.
fn after_gzip_header(member: &[u8]) -> std::io::Result<&[u8]> {
    if !member.starts_with(&GZIP_MAGIC) {
        return Err(gzip_error("not a gzip member"));
    }
    let fixed = member.get(2..10).ok_or_else(gzip_cut_short)?;
    if fixed[0] != GZIP_DEFLATE {
        return Err(gzip_error("a gzip member of an unknown compression method"));
    }
    let flags = fixed[1];
    let mut rest = &member[10..];
    if flags & FEXTRA != 0 {
        let length = rest.get(..2).ok_or_else(gzip_cut_short)?;
        let length = usize::from(u16::from_le_bytes([length[0], length[1]]));
        rest = rest.get(2 + length..).ok_or_else(gzip_cut_short)?;
    }
    for text in [FNAME, FCOMMENT] {
        if flags & text != 0 {
            rest = memchr(0, rest).map_or(&[][..], |end| &rest[end + 1..]);
        }
    }
    if flags & FHCRC != 0 {
        rest = rest.get(2..).ok_or_else(gzip_cut_short)?;
    }
    Ok(rest)
}

/// Inflates one member onto `bytes` and returns what follows its trailer, the CRC-32 and
/// the length of its text checked against the trailer's.
fn gunzip_member<'a>(member: &'a [u8], bytes: &mut Vec<u8>) -> std::io::Result<&'a [u8]> {
    let mut rest = after_gzip_header(member)?;
    let from = bytes.len();
    let mut inflater = Decompress::new(false);
    loop {
        if bytes.capacity() - bytes.len() < INFLATE_ROOM {
            bytes.reserve(INFLATE_ROOM);
        }
        let (read, written) = (inflater.total_in(), inflater.total_out());
        let status = inflater
            .decompress_vec(rest, bytes, FlushDecompress::None)
            .map_err(|e| gzip_error(&e.to_string()))?;
        rest = &rest[(inflater.total_in() - read) as usize..];
        if status == Status::StreamEnd {
            break;
        }
        if rest.is_empty() && inflater.total_out() == written {
            return Err(gzip_cut_short());
        }
    }
    let trailer = rest.get(..8).ok_or_else(gzip_cut_short)?;
    let text = &bytes[from..];
    let mut crc = Crc::new();
    crc.update(text);
    let stored = |at: usize| {
        u32::from_le_bytes([
            trailer[at],
            trailer[at + 1],
            trailer[at + 2],
            trailer[at + 3],
        ])
    };
    // The stored length is the text's length modulo 2^32.
    if stored(0) != crc.sum() || stored(4) != text.len() as u32 {
        return Err(gzip_error(
            "a gzip member's checksum or length is not its text's",
        ));
    }
    Ok(&rest[8..])
}

/// Offsets that cut `bytes` into about `parts` pieces, each ending after a `\n`, so no
/// piece ends inside a line or between the two bytes of a `\r\n`.
fn piece_bounds(bytes: &[u8], parts: usize) -> Vec<usize> {
    let mut bounds = vec![0];
    for k in 1..parts {
        let from = bytes.len() / parts * k;
        let last = *bounds.last().expect("bounds start at 0");
        if from <= last {
            continue;
        }
        match memchr(b'\n', &bytes[from..]) {
            Some(i) if from + i + 1 < bytes.len() => bounds.push(from + i + 1),
            _ => break,
        }
    }
    bounds.push(bytes.len());
    bounds
}

fn piece_lines(bytes: &[u8], base: usize) -> Result<Vec<Line>, String> {
    if std::str::from_utf8(bytes).is_err() {
        return Err("the text is not UTF-8".to_string());
    }
    let mut lines = Vec::new();
    let mut at = 0;
    while at < bytes.len() {
        let (end, next) = match memchr2(b'\n', b'\r', &bytes[at..]) {
            Some(i) => {
                let end = at + i;
                let crlf = bytes[end] == b'\r' && bytes.get(end + 1) == Some(&b'\n');
                (end, end + 1 + usize::from(crlf))
            }
            None => (bytes.len(), bytes.len()),
        };
        let len = u32::try_from(end - at).map_err(|_| "a line is longer than 4 GiB".to_string())?;
        let text = &bytes[at..end];
        let flags = match text {
            [b'#', b'#', ..] => HEADER,
            [b'#', ..] => HEADER | COLUMNS,
            _ => 0,
        };
        lines.push(Line {
            start: (base + at) as u64,
            len,
            flags,
        });
        at = next;
    }
    Ok(lines)
}

fn split_lines(bytes: &[u8], threads: usize) -> Result<Vec<Line>, String> {
    let bounds = piece_bounds(bytes, threads.max(1));
    let pieces: Vec<Result<Vec<Line>, String>> = std::thread::scope(|scope| {
        let workers: Vec<_> = bounds
            .windows(2)
            .map(|w| {
                let (from, to) = (w[0], w[1]);
                scope.spawn(move || piece_lines(&bytes[from..to], from))
            })
            .collect();
        workers
            .into_iter()
            .map(|w| w.join().expect("line splitter thread"))
            .collect()
    });
    let mut lines = Vec::new();
    for piece in pieces {
        lines.extend(piece?);
    }
    if u32::try_from(lines.len()).is_err() {
        return Err("more lines than this tool indexes".to_string());
    }
    Ok(lines)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn texts(bytes: &[u8], threads: usize) -> Vec<String> {
        split_lines(bytes, threads)
            .unwrap()
            .iter()
            .map(|l| {
                String::from_utf8(
                    bytes[l.start as usize..l.start as usize + l.len as usize].to_vec(),
                )
                .unwrap()
            })
            .collect()
    }

    #[test]
    fn every_line_end_of_text_mode_ends_a_line() {
        let bytes = b"#h\tx\r\nr1\t1\rr2\t2\nr3\t3";
        for threads in [1, 2, 7] {
            assert_eq!(texts(bytes, threads), ["#h\tx", "r1\t1", "r2\t2", "r3\t3"]);
        }
    }

    #[test]
    fn an_empty_line_before_the_end_stays_and_the_last_end_adds_none() {
        assert_eq!(texts(b"a\n\n", 1), ["a", ""]);
        assert_eq!(texts(b"a\r", 1), ["a"]);
        assert_eq!(texts(b"", 3), Vec::<String>::new());
        assert_eq!(texts(b"\n", 1), [""]);
    }

    #[test]
    fn pieces_never_cut_a_line() {
        let mut bytes = Vec::new();
        for i in 0..1000 {
            bytes.extend_from_slice(format!("row{i}\tvalue\r\n").as_bytes());
        }
        let one = texts(&bytes, 1);
        assert_eq!(one.len(), 1000);
        for threads in [2, 3, 16, 64] {
            assert_eq!(texts(&bytes, threads), one);
        }
    }

    #[test]
    fn header_lines_are_flagged() {
        let bytes = b"## note\n#a\tb\nx\ty\n";
        let flags: Vec<u8> = split_lines(bytes, 1)
            .unwrap()
            .iter()
            .map(|l| l.flags)
            .collect();
        assert_eq!(flags, [HEADER, HEADER | COLUMNS, 0]);
    }

    #[test]
    fn zero_padding_after_a_gzip_member_is_skipped_and_other_bytes_are_refused() {
        use flate2::write::GzEncoder;
        use std::io::Write;
        let member = |text: &[u8]| {
            let mut gz = GzEncoder::new(Vec::new(), flate2::Compression::default());
            gz.write_all(text).unwrap();
            gz.finish().unwrap()
        };
        let mut padded = member(b"#h\n");
        padded.extend_from_slice(&[0; 300]);
        padded.extend(member(b"row\n"));
        padded.extend_from_slice(&[0; 7]);
        assert_eq!(gunzip(&padded).unwrap(), b"#h\nrow\n");
        assert_eq!(gunzip(b"").unwrap(), b"");
        let mut trailing = member(b"#h\n");
        trailing.extend_from_slice(b"not a member");
        assert!(gunzip(&trailing).is_err());
        assert!(gunzip(&[0; 16]).is_err());
    }

    /// One gzip member holding `text` under a header of `flags`, `fields` standing between
    /// the fixed header and the deflate stream.
    fn member_with(flags: u8, fields: &[u8], text: &[u8]) -> Vec<u8> {
        use flate2::write::DeflateEncoder;
        use std::io::Write;
        let mut deflate = DeflateEncoder::new(Vec::new(), flate2::Compression::default());
        deflate.write_all(text).unwrap();
        let mut crc = Crc::new();
        crc.update(text);
        let mut member = vec![0x1f, 0x8b, 8, flags, 0, 0, 0, 0, 0, 0xff];
        member.extend_from_slice(fields);
        member.extend(deflate.finish().unwrap());
        member.extend_from_slice(&crc.sum().to_le_bytes());
        member.extend_from_slice(&(text.len() as u32).to_le_bytes());
        member
    }

    #[test]
    fn a_gzip_header_is_read_past_without_a_check_of_its_checksum_or_reserved_flags() {
        let text = b"#h\nrow\n";
        assert_eq!(gunzip(&member_with(0, b"", text)).unwrap(), text);
        // A header checksum that is not the header's.
        assert_eq!(
            gunzip(&member_with(FHCRC, &[0x12, 0x34], text)).unwrap(),
            text
        );
        // The three flag bits the format reserves.
        for reserved in [0x20, 0x40, 0x80, 0xe0] {
            assert_eq!(gunzip(&member_with(reserved, b"", text)).unwrap(), text);
        }
        // An extra field, a file name, a comment and a header checksum, in that order.
        let all = FEXTRA | FNAME | FCOMMENT | FHCRC;
        let fields = b"\x03\x00abcname\0a note\0\xff\xff";
        assert_eq!(gunzip(&member_with(all, fields, text)).unwrap(), text);
        let mut two = member_with(FNAME, b"one\0", b"#h\n");
        two.extend(member_with(0x20, b"", b"row\n"));
        assert_eq!(gunzip(&two).unwrap(), text);
    }

    #[test]
    fn a_gzip_member_that_gzip_open_refuses_is_refused() {
        let text = b"#h\nrow\n";
        let whole = member_with(0, b"", text);
        let last = whole.len() - 1;
        // A stored length and a stored checksum that are not the text's.
        let mut wrong = whole.clone();
        wrong[last] ^= 1;
        assert!(gunzip(&wrong).is_err());
        let mut wrong = whole.clone();
        wrong[last - 4] ^= 1;
        assert!(gunzip(&wrong).is_err());
        // Another compression method, a member cut inside its header, its stream or its
        // trailer, and a file name that never ends.
        let mut other = whole.clone();
        other[2] = 7;
        assert!(gunzip(&other).is_err());
        for cut in [1, 2, 9, 10, 12, last] {
            assert!(gunzip(&whole[..cut]).is_err(), "cut at {cut}");
        }
        assert!(gunzip(&member_with(FEXTRA, b"\x09\x00abc", b"")[..15]).is_err());
        let mut unended = vec![0x1f, 0x8b, 8, FNAME, 0, 0, 0, 0, 0, 0xff];
        unended.extend_from_slice(b"a name without its end");
        assert!(gunzip(&unended).is_err());
    }

    #[test]
    fn text_that_is_not_utf8_is_refused() {
        assert!(split_lines(b"#h\nr\t\xff\n", 1).is_err());
    }
}
