// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A gzip file as CPython's `gzip` module reads it, and the blocks of a BGZF file.
//!
//! `GzipMembers` reads the members one after another. A member is the two magic bytes,
//! the deflate method, a flag byte, six bytes that are not looked at, the extra field, the
//! name and the comment the flags announce (a name or comment the file ends inside ends
//! there), two bytes of header checksum that are not checked, the deflate stream, then the
//! CRC32 checksum and the length modulo 2^32 of what it inflates to, both checked. The zero bytes
//! that pad the file after a member are skipped, and the file ends wherever a member would
//! start and no byte is left, a zero-byte file included. Any other byte where a member
//! would start is a corrupt stream, zero bytes before the first member among them.
//!
//! A BGZF block is a member whose header is the eighteen bytes `bgzip` writes, which carry
//! the block's length. `inflate_blocks` inflates a run of such blocks without reading them
//! as a stream, and says so where a block is not exactly what `GzipMembers` would read,
//! so that the caller reads those bytes through `GzipMembers` instead.

use flate2::{Decompress, FlushDecompress, Status};
use std::io::{self, BufRead, Read};

fn corrupt(detail: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("corrupt gzip stream: {detail}"),
    )
}

fn truncated() -> io::Error {
    io::Error::new(
        io::ErrorKind::UnexpectedEof,
        "compressed file ended before the end-of-stream marker was reached",
    )
}

fn read_byte(raw: &mut impl BufRead) -> io::Result<Option<u8>> {
    let byte = raw.fill_buf()?.first().copied();
    if byte.is_some() {
        raw.consume(1);
    }
    Ok(byte)
}

fn read_exact(raw: &mut impl BufRead, buf: &mut [u8]) -> io::Result<()> {
    raw.read_exact(buf).map_err(|e| match e.kind() {
        io::ErrorKind::UnexpectedEof => truncated(),
        _ => e,
    })
}

fn skip_exact(raw: &mut impl BufRead, mut count: usize) -> io::Result<()> {
    while count > 0 {
        let available = raw.fill_buf()?.len();
        if available == 0 {
            return Err(truncated());
        }
        let step = available.min(count);
        raw.consume(step);
        count -= step;
    }
    Ok(())
}

/// Skips a string that ends at its first zero byte or at the end of the file.
fn skip_zero_terminated(raw: &mut impl BufRead) -> io::Result<()> {
    loop {
        let available = raw.fill_buf()?;
        if available.is_empty() {
            return Ok(());
        }
        match memchr::memchr(0, available) {
            Some(at) => {
                raw.consume(at + 1);
                return Ok(());
            }
            None => {
                let all = available.len();
                raw.consume(all);
            }
        }
    }
}

fn skip_zero_bytes(raw: &mut impl BufRead) -> io::Result<()> {
    loop {
        let available = raw.fill_buf()?;
        let zeros = available.iter().take_while(|byte| **byte == 0).count();
        let more = zeros > 0 && zeros == available.len();
        raw.consume(zeros);
        if !more {
            return Ok(());
        }
    }
}

const FHCRC: u8 = 2;
const FEXTRA: u8 = 4;
const FNAME: u8 = 8;
const FCOMMENT: u8 = 16;

/// Reads a member's header; false where the file ends instead.
fn read_header(raw: &mut impl BufRead) -> io::Result<bool> {
    let Some(first) = read_byte(raw)? else {
        return Ok(false);
    };
    if first != 0x1f || read_byte(raw)? != Some(0x8b) {
        return Err(corrupt("not a gzip member"));
    }
    let mut fixed = [0u8; 8];
    read_exact(raw, &mut fixed)?;
    let [method, flags, ..] = fixed;
    if method != 8 {
        return Err(corrupt("unknown compression method"));
    }
    if flags & FEXTRA != 0 {
        let mut length = [0u8; 2];
        read_exact(raw, &mut length)?;
        skip_exact(raw, usize::from(u16::from_le_bytes(length)))?;
    }
    if flags & FNAME != 0 {
        skip_zero_terminated(raw)?;
    }
    if flags & FCOMMENT != 0 {
        skip_zero_terminated(raw)?;
    }
    if flags & FHCRC != 0 {
        skip_exact(raw, 2)?;
    }
    Ok(true)
}

enum State {
    /// The start of the file: a member or the end.
    Start,
    /// After a member: zero padding, then a member or the end.
    Between,
    /// Inside a member's deflate stream.
    Member,
    /// After a deflate stream, before its eight trailing bytes are checked.
    Trailer,
    Ended,
}

pub struct GzipMembers<R> {
    raw: R,
    state: State,
    inflater: Decompress,
    crc: crc32fast::Hasher,
    length: u32,
}

impl<R: BufRead> GzipMembers<R> {
    /// Reads `raw` from the start of a gzip file.
    pub fn new(raw: R) -> GzipMembers<R> {
        GzipMembers::starting(raw, State::Start)
    }

    /// Reads `raw` from the end of one member of a gzip file on.
    pub fn after_member(raw: R) -> GzipMembers<R> {
        GzipMembers::starting(raw, State::Between)
    }

    fn starting(raw: R, state: State) -> GzipMembers<R> {
        GzipMembers {
            raw,
            state,
            inflater: Decompress::new(false),
            crc: crc32fast::Hasher::new(),
            length: 0,
        }
    }

    fn check_trailer(&mut self) -> io::Result<()> {
        let mut trailer = [0u8; 8];
        read_exact(&mut self.raw, &mut trailer)?;
        let (crc, length) = trailer.split_at(4);
        let computed = std::mem::replace(&mut self.crc, crc32fast::Hasher::new()).finalize();
        if crc != computed.to_le_bytes() {
            return Err(corrupt("CRC check failed"));
        }
        if length != self.length.to_le_bytes() {
            return Err(corrupt("incorrect length of data produced"));
        }
        Ok(())
    }
}

impl<R: BufRead> Read for GzipMembers<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            match self.state {
                State::Ended => return Ok(0),
                State::Start | State::Between => {
                    if matches!(self.state, State::Between) {
                        skip_zero_bytes(&mut self.raw)?;
                    }
                    if read_header(&mut self.raw)? {
                        self.inflater.reset(false);
                        self.length = 0;
                        self.state = State::Member;
                    } else {
                        self.state = State::Ended;
                    }
                }
                State::Member => {
                    let input = self.raw.fill_buf()?;
                    if input.is_empty() {
                        return Err(truncated());
                    }
                    let (before_in, before_out) =
                        (self.inflater.total_in(), self.inflater.total_out());
                    let status = self
                        .inflater
                        .decompress(input, buf, FlushDecompress::None)
                        .map_err(|e| corrupt(&e.to_string()))?;
                    let used = (self.inflater.total_in() - before_in) as usize;
                    let made = (self.inflater.total_out() - before_out) as usize;
                    self.raw.consume(used);
                    self.crc.update(&buf[..made]);
                    self.length = self.length.wrapping_add(made as u32);
                    if status == Status::StreamEnd {
                        self.state = State::Trailer;
                    } else if used == 0 && made == 0 {
                        return Err(corrupt("the deflate stream cannot be read further"));
                    }
                    if made > 0 {
                        return Ok(made);
                    }
                }
                State::Trailer => {
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

/// The most a block may inflate to here; one that announces more is left to `GzipMembers`.
const BLOCK_LIMIT: usize = 1 << 20;

/// Appends to `out` what the whole BGZF blocks of `batch` inflate to. False where a block
/// is not one `GzipMembers` would read to the same bytes and the same end: a header that
/// is not `bgzip`'s, a deflate stream that does not end where the block's trailer starts,
/// or a trailer that does not hold the CRC32 checksum and the length of what was inflated.
pub fn inflate_blocks(batch: &[u8], inflater: &mut Decompress, out: &mut Vec<u8>) -> bool {
    let mut rest = batch;
    while !rest.is_empty() {
        let Some(length) = bgzf_block_len(rest).filter(|length| *length <= rest.len()) else {
            return false;
        };
        let (block, after) = rest.split_at(length);
        rest = after;
        let (deflated, trailer) = block[BGZF_HEADER..].split_at(length - BGZF_HEADER - 8);
        let crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
        let announced =
            u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]) as usize;
        if announced > BLOCK_LIMIT {
            return false;
        }
        let start = out.len();
        // One byte more than announced, so a stream that inflates to more is seen.
        out.resize(start + announced + 1, 0);
        inflater.reset(false);
        let ended = matches!(
            inflater.decompress(deflated, &mut out[start..], FlushDecompress::Finish),
            Ok(Status::StreamEnd)
        );
        if !ended
            || inflater.total_in() != deflated.len() as u64
            || inflater.total_out() != announced as u64
        {
            return false;
        }
        out.truncate(start + announced);
        if crc32fast::hash(&out[start..]) != crc {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{Compress, Compression, FlushCompress};

    fn deflate(data: &[u8]) -> Vec<u8> {
        let mut out = vec![0u8; data.len() + data.len() / 8 + 64];
        let mut compress = Compress::new(Compression::fast(), false);
        let status = compress
            .compress(data, &mut out, FlushCompress::Finish)
            .unwrap();
        assert_eq!(status, Status::StreamEnd);
        out.truncate(compress.total_out() as usize);
        out
    }

    fn trailer(data: &[u8]) -> Vec<u8> {
        [
            crc32fast::hash(data).to_le_bytes(),
            (data.len() as u32).to_le_bytes(),
        ]
        .concat()
    }

    /// A member with the header fields `flags` announces: an extra field of three bytes, a
    /// name, a comment and a header checksum that is wrong.
    fn member_with(flags: u8, data: &[u8]) -> Vec<u8> {
        let mut out = vec![0x1f, 0x8b, 8, flags, 1, 2, 3, 4, 0, 255];
        if flags & FEXTRA != 0 {
            out.extend_from_slice(&[3, 0, b'x', b'y', b'z']);
        }
        if flags & FNAME != 0 {
            out.extend_from_slice(b"name.vcf\0");
        }
        if flags & FCOMMENT != 0 {
            out.extend_from_slice(b"a comment\0");
        }
        if flags & FHCRC != 0 {
            out.extend_from_slice(&[0xde, 0xad]);
        }
        out.extend_from_slice(&deflate(data));
        out.extend_from_slice(&trailer(data));
        out
    }

    fn member(data: &[u8]) -> Vec<u8> {
        member_with(0, data)
    }

    fn block(data: &[u8]) -> Vec<u8> {
        let deflated = deflate(data);
        let length = BGZF_HEADER + deflated.len() + 8;
        let mut out = vec![0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0];
        out.extend_from_slice(&((length - 1) as u16).to_le_bytes());
        out.extend_from_slice(&deflated);
        out.extend_from_slice(&trailer(data));
        out
    }

    /// Reads `file` whole through a buffer of `capacity` bytes, so a member's end and the
    /// padding after it fall at every place in the buffer.
    fn read_all(file: &[u8], capacity: usize) -> io::Result<Vec<u8>> {
        let mut out = Vec::new();
        GzipMembers::new(io::BufReader::with_capacity(capacity, file)).read_to_end(&mut out)?;
        Ok(out)
    }

    #[test]
    fn members_are_read_one_after_another_and_zero_padding_after_a_member_is_skipped() {
        let files: [(Vec<u8>, &[u8]); 9] = [
            (Vec::new(), b""),
            (member(b""), b""),
            (member(b"one\n"), b"one\n"),
            ([member(b"one\n"), member(b"two\n")].concat(), b"one\ntwo\n"),
            ([member(b"one\n"), vec![0; 512]].concat(), b"one\n"),
            (
                [member(b"on"), vec![0; 7], member(b"e\ntwo\n"), vec![0]].concat(),
                b"one\ntwo\n",
            ),
            ([member(b""), vec![0; 3], member(b"")].concat(), b""),
            (
                [block(b"one\n"), block(b""), member(b"two\n"), block(b"")].concat(),
                b"one\ntwo\n",
            ),
            (
                [
                    member_with(FEXTRA | FNAME | FCOMMENT | FHCRC, b"one\n"),
                    member_with(FNAME, b"two\n"),
                    member_with(FHCRC | 0xe0, b"three\n"),
                ]
                .concat(),
                b"one\ntwo\nthree\n",
            ),
        ];
        for (file, expected) in &files {
            for capacity in [1, 2, 3, 5, 8, 64, 4096] {
                assert_eq!(
                    read_all(file, capacity).unwrap(),
                    *expected,
                    "capacity {capacity}, {} bytes",
                    file.len()
                );
            }
        }
    }

    #[test]
    fn any_other_byte_where_a_member_would_start_is_a_corrupt_stream() {
        let good = member(b"one\n");
        let mut wrong_crc = good.clone();
        let at = wrong_crc.len() - 8;
        wrong_crc[at] ^= 1;
        let mut wrong_length = good.clone();
        let at = wrong_length.len() - 4;
        wrong_length[at] ^= 1;
        let mut wrong_method = good.clone();
        wrong_method[2] = 7;
        let files = [
            vec![0u8; 4],
            [vec![0], good.clone()].concat(),
            [good.clone(), b"garbage".to_vec()].concat(),
            [good.clone(), vec![0; 9], b"x".to_vec()].concat(),
            [good.clone(), vec![0x1f]].concat(),
            good[..good.len() - 9].to_vec(),
            good[..good.len() - 1].to_vec(),
            good[..5].to_vec(),
            b"#plain text\n".to_vec(),
            wrong_crc,
            wrong_length,
            wrong_method,
            // A name the file ends inside ends there, and no deflate stream follows.
            vec![0x1f, 0x8b, 8, FNAME, 0, 0, 0, 0, 0, 255, b'n', b'a'],
        ];
        for file in &files {
            for capacity in [1, 4, 4096] {
                assert!(
                    read_all(file, capacity).is_err(),
                    "capacity {capacity}, {} bytes",
                    file.len()
                );
            }
        }
    }

    #[test]
    fn what_a_member_holds_is_read_before_its_trailer_fails() {
        let mut file = member(b"one\ntwo\n");
        let at = file.len() - 8;
        file[at] ^= 1;
        let mut reader = GzipMembers::new(io::BufReader::new(&file[..]));
        let mut text = [0u8; 64];
        assert_eq!(reader.read(&mut text).unwrap(), 8);
        assert_eq!(&text[..8], b"one\ntwo\n");
        assert!(reader.read(&mut text).is_err());
    }

    #[test]
    fn reading_goes_on_after_a_member_as_it_does_inside_a_file() {
        let rest = [vec![0u8; 5], member(b"two\n")].concat();
        let mut out = Vec::new();
        GzipMembers::after_member(io::BufReader::new(&rest[..]))
            .read_to_end(&mut out)
            .unwrap();
        assert_eq!(out, b"two\n");
        assert!(GzipMembers::new(io::BufReader::new(&rest[..]))
            .read_to_end(&mut Vec::new())
            .is_err());
    }

    #[test]
    fn blocks_inflate_to_what_the_member_reader_reads() {
        let long: Vec<u8> = (0..9_000u32)
            .flat_map(|i| format!("21\t{i}\t.\tA\tC\n").into_bytes())
            .collect();
        let batch = [
            block(b"one\n"),
            block(b""),
            block(&long[..65_280]),
            block(b"two\n"),
        ]
        .concat();
        let mut out = b"kept".to_vec();
        assert!(inflate_blocks(
            &batch,
            &mut Decompress::new(false),
            &mut out
        ));
        assert_eq!(&out[..4], b"kept");
        assert_eq!(out[4..], read_all(&batch, 4096).unwrap());
        assert_eq!(bgzf_block_len(&batch), Some(block(b"one\n").len()));
    }

    #[test]
    fn a_block_the_member_reader_would_read_otherwise_is_not_inflated() {
        let good = block(b"one\ntwo\n");
        let mut wrong_crc = good.clone();
        let at = wrong_crc.len() - 8;
        wrong_crc[at] ^= 1;
        let mut wrong_length = good.clone();
        let at = wrong_length.len() - 4;
        wrong_length[at] ^= 1;
        let mut longer = good.clone();
        longer[16] += 1;
        longer.insert(BGZF_HEADER + 3, 0);
        let refused = [
            wrong_crc,
            wrong_length,
            longer,
            good[..good.len() - 1].to_vec(),
            [good.clone(), member(b"x")].concat(),
            [good.clone(), vec![0; 28]].concat(),
        ];
        for batch in &refused {
            let mut out = Vec::new();
            assert!(
                !inflate_blocks(batch, &mut Decompress::new(false), &mut out),
                "{} bytes",
                batch.len()
            );
        }
        for header in [&member(b"one\n")[..], &good[..17], &[0u8; 18][..]] {
            assert_eq!(bgzf_block_len(header), None);
        }
        // The shortest block a header can announce is its header and its trailer.
        let mut short = good.clone();
        short[16..18].copy_from_slice(&24u16.to_le_bytes());
        assert_eq!(bgzf_block_len(&short), None);
    }
}
