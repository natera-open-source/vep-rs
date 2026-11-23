// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The members of a gzip file, one after another, as `gzip.open` of the interpreter reads
//! them, a failure worded as the interpreter words it.
//!
//! A zero-byte file is an empty text, and zero bytes after a member are padding, skipped
//! up to the next member or the end. A member's header checksum is skipped without a
//! check, a flag bit the format reserves is read past, a file name or comment may run to
//! the end of the file, and the checksum and length of a member's text are held to its
//! trailer.
//!
//! The inflater loses the words of a failure its fast decoder meets. Such a deflate
//! stream is inflated again from its start, its last stretch through an output too small
//! for the fast decoder, and the failure takes the words met there.

use flate2::{Crc, Decompress, FlushDecompress, Status};
use std::io::{self, Read};

const FHCRC: u8 = 2;
const FEXTRA: u8 = 4;
const FNAME: u8 = 8;
const FCOMMENT: u8 = 16;
const INPUT: usize = 1 << 18;

const CUT_SHORT: &str = "Compressed file ended before the end-of-stream marker was reached";

/// Why a gzip file could not be read: the system's error, or the text of the exception
/// `gzip.open` raises.
#[derive(Debug)]
pub(crate) enum Failure {
    System(io::Error),
    Text(String),
}

impl From<io::Error> for Failure {
    fn from(error: io::Error) -> Failure {
        Failure::System(error)
    }
}

fn cut_short() -> Failure {
    Failure::Text(CUT_SHORT.to_string())
}

/// `repr()` of a `bytes` value.
fn bytes_repr(bytes: &[u8]) -> String {
    let quote = if bytes.contains(&b'\'') && !bytes.contains(&b'"') {
        '"'
    } else {
        '\''
    };
    let mut repr = format!("b{quote}");
    for &byte in bytes {
        match byte {
            b'\\' => repr.push_str("\\\\"),
            b'\t' => repr.push_str("\\t"),
            b'\n' => repr.push_str("\\n"),
            b'\r' => repr.push_str("\\r"),
            _ if char::from(byte) == quote => {
                repr.push('\\');
                repr.push(quote);
            }
            0x20..=0x7e => repr.push(char::from(byte)),
            _ => repr.push_str(&format!("\\x{byte:02x}")),
        }
    }
    repr.push(quote);
    repr
}

enum State {
    /// Before a member's header, or at the end of the file.
    Between,
    /// Inside a member's deflate stream.
    Inside,
    Ended,
}

/// The output's file from its start.
pub(crate) type Opened = Box<dyn Read + Send>;
pub(crate) type Open = Box<dyn Fn() -> io::Result<Opened> + Send>;

/// The words the inflater leaves where its fast decoder met the failure.
const WORDS_LOST: &str = "repeated call with bad state";
/// An output under the length the inflater's fast decoder asks for.
const SLOW_OUTPUT: usize = 200;

pub(crate) struct Gunzip {
    open: Open,
    file: Opened,
    /// The bytes of the file read so far.
    read: u64,
    /// Where in the file the deflate stream of the member in work starts.
    stream_start: u64,
    input: Vec<u8>,
    at: usize,
    filled: usize,
    file_ended: bool,
    state: State,
    inflater: Decompress,
    crc: Crc,
    length: u64,
}

impl Gunzip {
    pub(crate) fn new(open: Open) -> io::Result<Gunzip> {
        Ok(Gunzip {
            file: open()?,
            open,
            read: 0,
            stream_start: 0,
            input: vec![0; INPUT],
            at: 0,
            filled: 0,
            file_ended: false,
            state: State::Between,
            inflater: Decompress::new(false),
            crc: Crc::new(),
            length: 0,
        })
    }

    /// The words of the failure of the deflate stream in work, which fails once
    /// `written` bytes of its text are out: the words the inflater's slow decoder has
    /// for it.
    fn words_again(&self, written: u64) -> Option<String> {
        let mut file = (self.open)().ok()?;
        io::copy(&mut file.by_ref().take(self.stream_start), &mut io::sink()).ok()?;
        let mut inflater = Decompress::new(false);
        let mut input = vec![0; INPUT];
        let mut out = vec![0; INPUT];
        let (mut at, mut filled) = (0, 0);
        loop {
            if at == filled {
                filled = file.read(&mut input).ok()?;
                at = 0;
                if filled == 0 {
                    return None;
                }
            }
            let fast = inflater.total_out() + (out.len() as u64) < written;
            let room = if fast { out.len() } else { SLOW_OUTPUT };
            let (read, written_before) = (inflater.total_in(), inflater.total_out());
            match inflater.decompress(&input[at..filled], &mut out[..room], FlushDecompress::None) {
                Ok(Status::StreamEnd) => return None,
                Ok(_) => {
                    let read = (inflater.total_in() - read) as usize;
                    if read == 0 && inflater.total_out() == written_before {
                        return None;
                    }
                    at += read;
                }
                Err(error) => {
                    return error
                        .message()
                        .filter(|words| *words != WORDS_LOST)
                        .map(str::to_string)
                }
            }
        }
    }

    /// Reads more of the file when every byte read is used up; `false` at the file's end.
    fn refill(&mut self) -> io::Result<bool> {
        if self.at < self.filled {
            return Ok(true);
        }
        if self.file_ended {
            return Ok(false);
        }
        self.at = 0;
        self.filled = 0;
        loop {
            match self.file.read(&mut self.input) {
                Ok(0) => {
                    self.file_ended = true;
                    return Ok(false);
                }
                Ok(read) => {
                    self.filled = read;
                    self.read += read as u64;
                    return Ok(true);
                }
                Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                Err(error) => return Err(error),
            }
        }
    }

    fn byte(&mut self) -> io::Result<Option<u8>> {
        if !self.refill()? {
            return Ok(None);
        }
        self.at += 1;
        Ok(Some(self.input[self.at - 1]))
    }

    fn exactly<const N: usize>(&mut self) -> Result<[u8; N], Failure> {
        let mut bytes = [0; N];
        for byte in &mut bytes {
            *byte = self.byte()?.ok_or_else(cut_short)?;
        }
        Ok(bytes)
    }

    fn skip(&mut self, mut count: usize) -> Result<(), Failure> {
        while count > 0 {
            if !self.refill()? {
                return Err(cut_short());
            }
            let step = count.min(self.filled - self.at);
            self.at += step;
            count -= step;
        }
        Ok(())
    }

    /// Reads past a text that a zero byte or the end of the file ends.
    fn skip_text(&mut self) -> io::Result<()> {
        while let Some(byte) = self.byte()? {
            if byte == 0 {
                break;
            }
        }
        Ok(())
    }

    /// Reads a member's header; `false` at the end of the file.
    fn header(&mut self) -> Result<bool, Failure> {
        let Some(first) = self.byte()? else {
            return Ok(false);
        };
        let second = self.byte()?;
        if (first, second) != (0x1f, Some(0x8b)) {
            let magic: Vec<u8> = std::iter::once(first).chain(second).collect();
            return Err(Failure::Text(format!(
                "Not a gzipped file ({})",
                bytes_repr(&magic)
            )));
        }
        let [method, flags, ..] = self.exactly::<8>()?;
        if method != 8 {
            return Err(Failure::Text("Unknown compression method".to_string()));
        }
        if flags & FEXTRA != 0 {
            let length = u16::from_le_bytes(self.exactly::<2>()?);
            self.skip(usize::from(length))?;
        }
        if flags & FNAME != 0 {
            self.skip_text()?;
        }
        if flags & FCOMMENT != 0 {
            self.skip_text()?;
        }
        if flags & FHCRC != 0 {
            self.exactly::<2>()?;
        }
        Ok(true)
    }

    /// Holds the text of the member that just ended to its trailer and reads past the
    /// zero bytes after it.
    fn trailer(&mut self) -> Result<(), Failure> {
        let trailer = self.exactly::<8>()?;
        let stored = |at: usize| {
            u32::from_le_bytes([
                trailer[at],
                trailer[at + 1],
                trailer[at + 2],
                trailer[at + 3],
            ])
        };
        if stored(0) != self.crc.sum() {
            return Err(Failure::Text(format!(
                "CRC check failed {:#x} != {:#x}",
                stored(0),
                self.crc.sum()
            )));
        }
        // The stored length is the text's length modulo 2^32.
        if stored(4) != self.length as u32 {
            return Err(Failure::Text(
                "Incorrect length of data produced".to_string(),
            ));
        }
        while self.refill()? && self.input[self.at] == 0 {
            self.at += 1;
        }
        Ok(())
    }

    /// Inflates onto `out[*filled..]` until `out` is full or the file ends; `false` once
    /// the file has ended.
    pub(crate) fn read(&mut self, out: &mut [u8], filled: &mut usize) -> Result<bool, Failure> {
        loop {
            match self.state {
                State::Ended => return Ok(false),
                State::Between => {
                    if !self.header()? {
                        self.state = State::Ended;
                        return Ok(false);
                    }
                    self.inflater.reset(false);
                    self.crc.reset();
                    self.length = 0;
                    self.stream_start = self.read - (self.filled - self.at) as u64;
                    self.state = State::Inside;
                }
                State::Inside => {
                    if *filled == out.len() {
                        return Ok(true);
                    }
                    let more = self.refill()?;
                    let (read, written) = (self.inflater.total_in(), self.inflater.total_out());
                    let status = match self.inflater.decompress(
                        &self.input[self.at..self.filled],
                        &mut out[*filled..],
                        FlushDecompress::None,
                    ) {
                        Ok(status) => status,
                        Err(error) => {
                            let words = match error.message() {
                                Some(WORDS_LOST) | None => {
                                    self.words_again(self.inflater.total_out())
                                }
                                Some(words) => Some(words.to_string()),
                            };
                            return Err(Failure::Text(format!(
                                "Error -3 while decompressing data: {}",
                                words.as_deref().unwrap_or("invalid input data")
                            )));
                        }
                    };
                    let read = (self.inflater.total_in() - read) as usize;
                    let written = (self.inflater.total_out() - written) as usize;
                    self.at += read;
                    self.crc.update(&out[*filled..*filled + written]);
                    self.length += written as u64;
                    *filled += written;
                    if status == Status::StreamEnd {
                        self.trailer()?;
                        self.state = State::Between;
                    } else if read == 0 && written == 0 && (!more || *filled < out.len()) {
                        // The stream wants bytes the file does not hold.
                        return Err(cut_short());
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::{DeflateEncoder, GzEncoder};
    use std::io::Write;

    /// The text of `file`, read into pieces of `piece` bytes.
    fn gunzip_by(file: &[u8], piece: usize) -> Result<Vec<u8>, String> {
        let bytes = file.to_vec();
        let open: Open = Box::new(move || Ok(Box::new(io::Cursor::new(bytes.clone()))));
        let mut gunzip = Gunzip::new(open).unwrap();
        let mut text = Vec::new();
        loop {
            let mut out = vec![0; piece];
            let mut filled = 0;
            let more = gunzip
                .read(&mut out, &mut filled)
                .map_err(|failure| match failure {
                    Failure::Text(text) => text,
                    Failure::System(error) => error.to_string(),
                })?;
            text.extend_from_slice(&out[..filled]);
            if !more {
                return Ok(text);
            }
        }
    }

    fn gunzip(file: &[u8]) -> Result<Vec<u8>, String> {
        let whole = gunzip_by(file, 1 << 16);
        for piece in [1, 3, 4096] {
            assert_eq!(gunzip_by(file, piece), whole, "pieces of {piece} bytes");
        }
        whole
    }

    fn member(text: &[u8]) -> Vec<u8> {
        let mut gz = GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(text).unwrap();
        gz.finish().unwrap()
    }

    /// One gzip member holding `text` under a header of `flags`, `fields` standing between
    /// the fixed header and the deflate stream.
    fn member_with(flags: u8, fields: &[u8], text: &[u8]) -> Vec<u8> {
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
    fn zero_padding_after_a_gzip_member_is_skipped_and_other_bytes_are_refused() {
        let mut padded = member(b"#h\n");
        padded.extend_from_slice(&[0; 300]);
        padded.extend(member(b"row\n"));
        padded.extend_from_slice(&[0; 7]);
        assert_eq!(gunzip(&padded).unwrap(), b"#h\nrow\n");
        assert_eq!(gunzip(b"").unwrap(), b"");
        let mut trailing = member(b"#h\n");
        trailing.extend_from_slice(b"not a member");
        assert_eq!(gunzip(&trailing).unwrap_err(), "Not a gzipped file (b'no')");
        assert_eq!(
            gunzip(&[0; 16]).unwrap_err(),
            "Not a gzipped file (b'\\x00\\x00')"
        );
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
        assert_eq!(
            gunzip(&wrong).unwrap_err(),
            "Incorrect length of data produced"
        );
        let mut wrong = whole.clone();
        wrong[last - 4] ^= 1;
        let mut crc = Crc::new();
        crc.update(text);
        assert_eq!(
            gunzip(&wrong).unwrap_err(),
            format!(
                "CRC check failed {:#x} != {:#x}",
                crc.sum() ^ 1 << 24,
                crc.sum()
            )
        );
        // Another compression method, a member cut inside its header, its stream or its
        // trailer, and a file name that never ends.
        let mut other = whole.clone();
        other[2] = 7;
        assert_eq!(gunzip(&other).unwrap_err(), "Unknown compression method");
        assert_eq!(
            gunzip(&whole[..1]).unwrap_err(),
            "Not a gzipped file (b'\\x1f')"
        );
        for cut in [2, 9, 10, 12, last] {
            assert_eq!(
                gunzip(&whole[..cut]).unwrap_err(),
                CUT_SHORT,
                "cut at {cut}"
            );
        }
        assert_eq!(
            gunzip(&member_with(FEXTRA, b"\x09\x00abc", b"")[..15]).unwrap_err(),
            CUT_SHORT
        );
        let mut unended = vec![0x1f, 0x8b, 8, FNAME, 0, 0, 0, 0, 0, 0xff];
        unended.extend_from_slice(b"a name without its end");
        assert_eq!(gunzip(&unended).unwrap_err(), CUT_SHORT);
    }

    #[test]
    fn a_deflate_stream_zlib_refuses_is_refused_in_its_words() {
        // A final block of the reserved type.
        let mut member = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff, 0x07];
        member.extend_from_slice(&[0; 8]);
        assert_eq!(
            gunzip(&member).unwrap_err(),
            "Error -3 while decompressing data: invalid block type"
        );
    }

    #[test]
    fn the_first_bytes_after_a_member_are_named_as_the_interpreter_writes_bytes() {
        for (bytes, repr) in [
            (&b"'\""[..], "b'\\'\"'"),
            (b"'a", "b\"'a\""),
            (b"\\\n", "b'\\\\\\n'"),
            (b"\xff\t", "b'\\xff\\t'"),
            (b"\r\x7f", "b'\\r\\x7f'"),
            (b"n", "b'n'"),
        ] {
            assert_eq!(bytes_repr(bytes), repr);
        }
    }
}
