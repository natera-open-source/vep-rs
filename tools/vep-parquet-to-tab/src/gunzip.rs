// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A gzip stream read member by member, as Python's `gzip` module reads one: any number of
//! members, zero bytes between them skipped, and each way a stream can be broken reported
//! in that module's words.

use std::io::{self, BufRead};

use flate2::{Crc, Decompress, FlushDecompress, Status};

use crate::pyfmt::py_bytes_repr;

pub const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const DEFLATE: u8 = 8;
const FHCRC: u8 = 2;
const FEXTRA: u8 = 4;
const FNAME: u8 = 8;
const FCOMMENT: u8 = 16;
const TRUNCATED: &str = "Compressed file ended before the end-of-stream marker was reached";

/// Why a read of the input ended.
#[derive(Debug)]
pub enum ReadFailure {
    Io(io::Error),
    /// A broken gzip stream, in the words of Python's `gzip` module.
    Gzip(String),
    /// Deflate data no decompressor accepts.
    CorruptDeflate,
}

impl From<io::Error> for ReadFailure {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

fn truncated() -> ReadFailure {
    ReadFailure::Gzip(TRUNCATED.to_string())
}

enum State {
    /// Before a member's header, or at the end of the file.
    Header,
    Body,
    Trailer,
    End,
}

/// The decompressed bytes of a gzip file.
pub struct GzipReader<R> {
    inner: R,
    state: State,
    inflater: Decompress,
    crc: Crc,
}

impl<R: BufRead> GzipReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            state: State::Header,
            inflater: Decompress::new(false),
            crc: Crc::new(),
        }
    }

    /// The next byte, `None` at the end of the file.
    fn byte(&mut self) -> io::Result<Option<u8>> {
        let byte = self.inner.fill_buf()?.first().copied();
        if byte.is_some() {
            self.inner.consume(1);
        }
        Ok(byte)
    }

    /// Exactly `N` bytes, or the truncation failure.
    fn exact<const N: usize>(&mut self) -> Result<[u8; N], ReadFailure> {
        let mut bytes = [0u8; N];
        for byte in &mut bytes {
            *byte = self.byte()?.ok_or_else(truncated)?;
        }
        Ok(bytes)
    }

    /// Reads a member's header; `false` at the end of the file. A name or a comment that
    /// the file ends inside is taken as ended.
    fn header(&mut self) -> Result<bool, ReadFailure> {
        let Some(first) = self.byte()? else {
            return Ok(false);
        };
        let mut magic = vec![first];
        magic.extend(self.byte()?);
        if magic != GZIP_MAGIC {
            return Err(ReadFailure::Gzip(format!(
                "Not a gzipped file ({})",
                py_bytes_repr(&magic)
            )));
        }
        let [method, flags, ..] = self.exact::<8>()?;
        if method != DEFLATE {
            return Err(ReadFailure::Gzip("Unknown compression method".to_string()));
        }
        if flags & FEXTRA != 0 {
            let length = u16::from_le_bytes(self.exact::<2>()?);
            for _ in 0..length {
                self.exact::<1>()?;
            }
        }
        for text in [FNAME, FCOMMENT] {
            if flags & text != 0 {
                while self.byte()?.is_some_and(|byte| byte != 0) {}
            }
        }
        if flags & FHCRC != 0 {
            self.exact::<2>()?;
        }
        Ok(true)
    }

    /// Holds a finished member to its stored checksum and length, then skips the zero
    /// bytes after it.
    fn trailer(&mut self) -> Result<(), ReadFailure> {
        let trailer = self.exact::<8>()?;
        let stored = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
        let length = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]);
        if stored != self.crc.sum() {
            return Err(ReadFailure::Gzip(format!(
                "CRC check failed {stored:#x} != {:#x}",
                self.crc.sum()
            )));
        }
        if length != self.crc.amount() {
            return Err(ReadFailure::Gzip(
                "Incorrect length of data produced".to_string(),
            ));
        }
        while self.inner.fill_buf()?.first() == Some(&0) {
            self.inner.consume(1);
        }
        Ok(())
    }

    /// Fills `out` with decompressed bytes; 0 at the end of the last member.
    pub fn read(&mut self, out: &mut [u8]) -> Result<usize, ReadFailure> {
        loop {
            match self.state {
                State::End => return Ok(0),
                State::Header => {
                    if !self.header()? {
                        self.state = State::End;
                        return Ok(0);
                    }
                    self.inflater.reset(false);
                    self.crc.reset();
                    self.state = State::Body;
                }
                State::Body => {
                    let input = self.inner.fill_buf()?;
                    let at_end = input.is_empty();
                    let (consumed, produced) =
                        (self.inflater.total_in(), self.inflater.total_out());
                    let status = self
                        .inflater
                        .decompress(input, out, FlushDecompress::None)
                        .map_err(|_| ReadFailure::CorruptDeflate)?;
                    let consumed = (self.inflater.total_in() - consumed) as usize;
                    let produced = (self.inflater.total_out() - produced) as usize;
                    self.inner.consume(consumed);
                    self.crc.update(&out[..produced]);
                    if status == Status::StreamEnd {
                        self.state = State::Trailer;
                    } else if produced == 0 && at_end {
                        return Err(truncated());
                    }
                    if produced > 0 {
                        return Ok(produced);
                    }
                }
                State::Trailer => {
                    self.trailer()?;
                    self.state = State::Header;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    fn gzip(data: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(data).unwrap();
        encoder.finish().unwrap()
    }

    fn read_all(bytes: &[u8]) -> Result<Vec<u8>, ReadFailure> {
        let mut reader = GzipReader::new(bytes);
        let mut data = Vec::new();
        let mut buffer = [0u8; 97];
        loop {
            match reader.read(&mut buffer)? {
                0 => return Ok(data),
                n => data.extend_from_slice(&buffer[..n]),
            }
        }
    }

    fn message(bytes: &[u8]) -> String {
        match read_all(bytes) {
            Err(ReadFailure::Gzip(message)) => message,
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn members_and_zero_padding_read_as_one_stream() {
        let text = b"#Feature\ta\n".repeat(100);
        assert_eq!(read_all(&gzip(&text)).unwrap(), text);
        let mut two = gzip(b"#Feature\n");
        two.extend_from_slice(&[0; 7]);
        two.extend(gzip(b"a\nb\n"));
        two.extend_from_slice(&[0; 3]);
        assert_eq!(read_all(&two).unwrap(), b"#Feature\na\nb\n");
        assert_eq!(read_all(&gzip(b"")).unwrap(), b"");
    }

    #[test]
    fn header_fields_are_skipped() {
        let plain = gzip(b"rows\n");
        let mut with_fields = plain[..3].to_vec();
        with_fields.push(FEXTRA | FNAME | FCOMMENT | FHCRC);
        with_fields.extend_from_slice(&plain[4..10]);
        with_fields.extend_from_slice(&[3, 0, b'x', b'y', b'z']);
        with_fields.extend_from_slice(b"name\0comment\0");
        with_fields.extend_from_slice(&[0xaa, 0xbb]);
        with_fields.extend_from_slice(&plain[10..]);
        assert_eq!(read_all(&with_fields).unwrap(), b"rows\n");
    }

    #[test]
    fn a_broken_stream_is_reported_in_the_gzip_modules_words() {
        let good = gzip(&b"#Feature\ta\n".repeat(100));
        assert_eq!(message(&good[..good.len() - 8]), TRUNCATED);
        assert_eq!(message(&good[..good.len() - 1]), TRUNCATED);
        assert_eq!(message(&good[..5]), TRUNCATED);
        assert_eq!(message(&good[..10]), TRUNCATED);
        assert_eq!(message(&[0x1f, 0x8b]), TRUNCATED);
        let mut checksum = good.clone();
        let at = checksum.len() - 8;
        checksum[at] ^= 0xff;
        let stored = u32::from_le_bytes(checksum[at..at + 4].try_into().unwrap());
        let computed = u32::from_le_bytes(good[at..at + 4].try_into().unwrap());
        assert_eq!(
            message(&checksum),
            format!("CRC check failed {stored:#x} != {computed:#x}")
        );
        let mut length = good.clone();
        *length.last_mut().unwrap() ^= 0xff;
        assert_eq!(message(&length), "Incorrect length of data produced");
        let mut method = good.clone();
        method[2] = 7;
        assert_eq!(message(&method), "Unknown compression method");
        let mut junk = good.clone();
        junk.extend_from_slice(b"junk");
        assert_eq!(message(&junk), "Not a gzipped file (b'ju')");
        let mut one = good.clone();
        one.push(0x1f);
        assert_eq!(message(&one), "Not a gzipped file (b'\\x1f')");
        let mut data = good.clone();
        for byte in &mut data[20..23] {
            *byte ^= 0xff;
        }
        assert!(matches!(
            read_all(&data),
            Err(ReadFailure::CorruptDeflate | ReadFailure::Gzip(_))
        ));
    }
}
