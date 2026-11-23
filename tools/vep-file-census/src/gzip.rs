// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A gzip file as `gzip -dc` reads it, and the blocks of a BGZF file.
//!
//! `Members` reads the members one after another. A member is the two magic bytes, the
//! deflate method, a flag byte, six bytes that are not looked at, the extra field, the name
//! and the comment the flags announce, the header checksum where the flags announce one,
//! the deflate stream, then the CRC32 checksum and the length modulo 2^32 of what it
//! inflates to; the three checks are made. What a member holds is read before a failure
//! inside it is reported.
//!
//! After a member the file ends where no byte is left, and where every byte left is zero.
//! Any other bytes there that do not start a member are trailing garbage: the reading ends
//! and `ended_in_garbage` says so, the state `gzip` warns of and leaves with status 2. A
//! file whose first bytes are not a member's is no gzip file, a zero-byte file included.
//! Under `pass_through`, the reading of `zcat -f`, bytes that do not start a member are
//! read as they are from there to the end of the file instead, at the start of the file as
//! after a member.
//!
//! A BGZF block is a member whose header is the eighteen bytes `bgzip` writes, which carry
//! the block's length. `inflate_blocks` inflates a run of such blocks without reading them
//! as a stream. A member under that header that announces more bytes than a block holds
//! here is no block: `Members` reads it.

use flate2::{Decompress, FlushDecompress, Status};
use std::io::{self, BufRead, Read};

fn corrupt(detail: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("corrupt gzip stream: {detail}"),
    )
}

fn truncated() -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, "unexpected end of file")
}

const FHCRC: u8 = 2;
const FEXTRA: u8 = 4;
const FNAME: u8 = 8;
const FCOMMENT: u8 = 16;
const ENCRYPTED: u8 = 0x20;
const RESERVED: u8 = 0xc0;

enum State {
    /// The start of the file.
    Start,
    /// After a member.
    Between,
    /// Inside a member's deflate stream.
    Member,
    /// After a deflate stream, before its eight trailing bytes are checked.
    Trailer,
    /// Bytes that are no member, read as they are to the end of the file.
    Plain,
    Ended,
}

pub struct Members<R> {
    raw: R,
    state: State,
    pass_through: bool,
    garbage: bool,
    /// Bytes read off the file that the plain reading starts with.
    pending: Vec<u8>,
    inflater: Decompress,
    crc: crc32fast::Hasher,
    length: u32,
    /// The header bytes read so far, kept where the flags announce a header checksum.
    header: Option<crc32fast::Hasher>,
}

impl<R: BufRead> Members<R> {
    /// Reads `raw` from the start of a gzip file.
    pub fn new(raw: R, pass_through: bool) -> Members<R> {
        Members::starting(raw, State::Start, pass_through)
    }

    /// Reads `raw` from the end of one member of a gzip file on.
    pub fn after_member(raw: R, pass_through: bool) -> Members<R> {
        Members::starting(raw, State::Between, pass_through)
    }

    fn starting(raw: R, state: State, pass_through: bool) -> Members<R> {
        Members {
            raw,
            state,
            pass_through,
            garbage: false,
            pending: Vec::new(),
            inflater: Decompress::new(false),
            crc: crc32fast::Hasher::new(),
            length: 0,
            header: None,
        }
    }

    /// Whether the reading ended at trailing garbage after a member.
    pub fn ended_in_garbage(&self) -> bool {
        self.garbage
    }

    fn peek(&mut self) -> io::Result<Option<u8>> {
        Ok(self.raw.fill_buf()?.first().copied())
    }

    fn byte(&mut self) -> io::Result<u8> {
        let byte = self.peek()?.ok_or_else(truncated)?;
        self.raw.consume(1);
        if let Some(header) = &mut self.header {
            header.update(&[byte]);
        }
        Ok(byte)
    }

    fn skip(&mut self, count: usize) -> io::Result<()> {
        for _ in 0..count {
            self.byte()?;
        }
        Ok(())
    }

    fn skip_zero_terminated(&mut self) -> io::Result<()> {
        while self.byte()? != 0 {}
        Ok(())
    }

    /// Reads what follows a member's two magic bytes, up to its deflate stream.
    fn header(&mut self) -> io::Result<()> {
        self.header = None;
        let method = self.byte()?;
        if method != 8 {
            return Err(corrupt("unknown compression method"));
        }
        let flags = self.byte()?;
        if flags & ENCRYPTED != 0 {
            return Err(corrupt("the member is encrypted"));
        }
        if flags & RESERVED != 0 {
            return Err(corrupt("the member has reserved flags"));
        }
        let mut fixed = [0u8; 6];
        for byte in &mut fixed {
            *byte = self.byte()?;
        }
        if flags & FHCRC != 0 {
            let mut seen = crc32fast::Hasher::new();
            seen.update(&[0x1f, 0x8b, method, flags]);
            seen.update(&fixed);
            self.header = Some(seen);
        }
        if flags & FEXTRA != 0 {
            let length = usize::from(self.byte()?) | usize::from(self.byte()?) << 8;
            self.skip(length)?;
        }
        if flags & FNAME != 0 {
            self.skip_zero_terminated()?;
        }
        if flags & FCOMMENT != 0 {
            self.skip_zero_terminated()?;
        }
        if let Some(seen) = self.header.take() {
            let computed = seen.finalize() & 0xffff;
            let stored = u32::from(self.byte()?) | u32::from(self.byte()?) << 8;
            if stored != computed {
                return Err(corrupt("header checksum mismatch"));
            }
        }
        Ok(())
    }

    /// Where a member would start: reads its magic bytes and header, or says how the
    /// reading goes on instead.
    fn boundary(&mut self, first: bool) -> io::Result<State> {
        let Some(byte) = self.peek()? else {
            return if first && !self.pass_through {
                Err(truncated())
            } else {
                Ok(State::Ended)
            };
        };
        if byte == 0x1f {
            self.raw.consume(1);
            let second = self.peek()?;
            if second == Some(0x8b) {
                self.raw.consume(1);
                self.header()?;
                return Ok(State::Member);
            }
            if self.pass_through {
                self.pending.push(0x1f);
                return Ok(State::Plain);
            }
            if first {
                return Err(corrupt("not in gzip format"));
            }
            if second.is_none() {
                return Err(truncated());
            }
            self.garbage = true;
            return Ok(State::Ended);
        }
        if self.pass_through {
            return Ok(State::Plain);
        }
        if first {
            return Err(corrupt("not in gzip format"));
        }
        if byte == 0 {
            // Zero bytes to the end of the file are padding.
            loop {
                let ahead = self.raw.fill_buf()?;
                if ahead.is_empty() {
                    return Ok(State::Ended);
                }
                let zeros = ahead.iter().take_while(|byte| **byte == 0).count();
                let all = zeros == ahead.len();
                self.raw.consume(zeros);
                if !all {
                    break;
                }
            }
        } else {
            // One byte that is not zero and nothing after it.
            self.raw.consume(1);
            if self.peek()?.is_none() {
                return Err(truncated());
            }
        }
        self.garbage = true;
        Ok(State::Ended)
    }

    fn check_trailer(&mut self) -> io::Result<()> {
        let mut trailer = [0u8; 8];
        for byte in &mut trailer {
            *byte = self.byte()?;
        }
        let (crc, length) = trailer.split_at(4);
        let computed = std::mem::replace(&mut self.crc, crc32fast::Hasher::new()).finalize();
        if crc != computed.to_le_bytes() {
            return Err(corrupt("crc error"));
        }
        if length != self.length.to_le_bytes() {
            return Err(corrupt("length error"));
        }
        Ok(())
    }
}

impl<R: BufRead> Read for Members<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            match self.state {
                State::Ended => return Ok(0),
                State::Start | State::Between => {
                    let first = matches!(self.state, State::Start);
                    self.state = self.boundary(first)?;
                    if matches!(self.state, State::Member) {
                        self.inflater.reset(false);
                        self.length = 0;
                    }
                }
                State::Plain => {
                    if !self.pending.is_empty() {
                        let count = self.pending.len().min(buf.len());
                        buf[..count].copy_from_slice(&self.pending[..count]);
                        self.pending.drain(..count);
                        return Ok(count);
                    }
                    let count = self.raw.read(buf)?;
                    if count == 0 {
                        self.state = State::Ended;
                    }
                    return Ok(count);
                }
                State::Member => {
                    let input = self.raw.fill_buf()?;
                    if input.is_empty() {
                        self.state = State::Ended;
                        return Err(truncated());
                    }
                    let (before_in, before_out) =
                        (self.inflater.total_in(), self.inflater.total_out());
                    let status = self.inflater.decompress(input, buf, FlushDecompress::None);
                    let used = (self.inflater.total_in() - before_in) as usize;
                    let made = (self.inflater.total_out() - before_out) as usize;
                    self.raw.consume(used);
                    self.crc.update(&buf[..made]);
                    self.length = self.length.wrapping_add(made as u32);
                    match status {
                        Ok(Status::StreamEnd) => self.state = State::Trailer,
                        Ok(_) if used == 0 && made == 0 => {
                            self.state = State::Ended;
                            return Err(corrupt("the deflate stream cannot be read further"));
                        }
                        Ok(_) => {}
                        Err(error) if made == 0 => {
                            self.state = State::Ended;
                            return Err(corrupt(&error.to_string()));
                        }
                        // What was inflated before the failure is handed on first; the
                        // next read meets the failure again.
                        Err(_) => {}
                    }
                    if made > 0 {
                        return Ok(made);
                    }
                }
                State::Trailer => {
                    self.state = State::Ended;
                    self.check_trailer()?;
                    self.state = State::Between;
                }
            }
        }
    }
}

/// The bytes of a BGZF block's header.
pub const BGZF_HEADER: usize = 18;

/// The length of the block whose header starts `bytes`, where that header is the one
/// `bgzip` writes: the magic bytes, deflate, the extra-field flag alone, and one extra
/// subfield `BC` of two bytes that holds the block's length less one.
pub fn bgzf_block_len(bytes: &[u8]) -> Option<usize> {
    let header = bytes.get(..BGZF_HEADER)?;
    if header[..4] != [0x1f, 0x8b, 8, FEXTRA] || header[10..16] != [6, 0, b'B', b'C', 2, 0] {
        return None;
    }
    let length = usize::from(u16::from_le_bytes([header[16], header[17]])) + 1;
    (length >= BGZF_HEADER + 8).then_some(length)
}

/// The most a block may inflate to here. `bgzip` puts at most 64 KiB in a block.
const BLOCK_LIMIT: usize = 1 << 20;

/// Whether `block`, a whole member under a BGZF block's header, announces in its last four
/// bytes that it inflates to no more than a block may here.
pub fn bgzf_block_in_bound(block: &[u8]) -> bool {
    block
        .last_chunk::<4>()
        .is_some_and(|announced| u32::from_le_bytes(*announced) as usize <= BLOCK_LIMIT)
}

/// Appends to `out` what the whole BGZF blocks of `batch` inflate to, each block held to
/// its CRC32 checksum and its length. What the blocks before a failing one inflate to
/// stays in `out`, and so does what the failing block inflated to before its failure.
pub fn inflate_blocks(
    batch: &[u8],
    inflater: &mut Decompress,
    out: &mut Vec<u8>,
) -> io::Result<()> {
    let mut rest = batch;
    while !rest.is_empty() {
        let length = bgzf_block_len(rest)
            .filter(|length| *length <= rest.len())
            .ok_or_else(truncated)?;
        let (block, after) = rest.split_at(length);
        rest = after;
        let (deflated, trailer) = block[BGZF_HEADER..].split_at(length - BGZF_HEADER - 8);
        let crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
        let announced =
            u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]) as usize;
        let start = out.len();
        // One byte more than announced, so a stream that inflates to more is seen.
        out.resize(start + announced.min(BLOCK_LIMIT) + 1, 0);
        inflater.reset(false);
        let status = inflater.decompress(deflated, &mut out[start..], FlushDecompress::Finish);
        let made = inflater.total_out() as usize;
        out.truncate(start + made);
        match status {
            Ok(Status::StreamEnd) if inflater.total_in() == deflated.len() as u64 => {}
            Ok(Status::StreamEnd) => {
                return Err(corrupt("a block's deflate stream ends before its trailer"))
            }
            Ok(_) if made > announced.min(BLOCK_LIMIT) => return Err(corrupt("length error")),
            Ok(_) => return Err(truncated()),
            Err(error) => return Err(corrupt(&error.to_string())),
        }
        if crc32fast::hash(&out[start..]) != crc {
            return Err(corrupt("crc error"));
        }
        if made != announced {
            return Err(corrupt("length error"));
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "gzip_tests.rs"]
mod tests;
