// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The text of a gzip file as GNU gzip's `zcat` writes it, read in process.
//!
//! `zcat` is not the reader `input::gunzip` is. It ends its reading at the first zero byte
//! after a member, where `gzip.open` takes zero bytes for padding and reads the member
//! after them. It refuses a member whose header checksum is not the header's, or whose
//! header carries the encryption flag or a reserved flag, where `gzip.open` reads past
//! all three. And what it inflated before the place it stops stays written: the members
//! before a refused one, and the text of a member whose trailer is cut short or is not
//! the text's.
//!
//! Three files are read here otherwise than `zcat` reads them, and the comparison's own
//! reading of the input (`input::input_records`) refuses each before this reading is
//! asked for: a deflate stream that is cut short, of which `zcat` may write one byte
//! fewer; a deflate stream zlib refuses, some of which `zcat` inflates on; and a file of
//! another container that `zcat` reads under a gzip name (zip, compress, pack, lzh).

use flate2::{Crc, Decompress, FlushDecompress, Status};

/// The two bytes that open a member, and the ones an early gzip wrote.
const MAGICS: [[u8; 2]; 2] = [[0x1f, 0x8b], [0x1f, 0x9e]];
const DEFLATE: u8 = 8;
const FHCRC: u8 = 0x02;
const FEXTRA: u8 = 0x04;
const FNAME: u8 = 0x08;
const FCOMMENT: u8 = 0x10;
/// The encryption flag and the two reserved flags.
const REFUSED_FLAGS: u8 = 0xe0;
/// The header's bytes before its optional fields.
const FIXED_HEADER: usize = 10;
/// The checksum and the length that close a member.
const TRAILER: usize = 8;
/// The room the inflater is given at a call.
const INFLATE_ROOM: usize = 1 << 20;

/// What `zcat` writes to its standard output of a file of the bytes `file`.
///
/// Members are read one after another while each is whole and the byte after it opens
/// the next. A file that opens with no member has no text. A member whose deflate stream
/// is cut short or refused gives the text inflated before that place.
pub fn text(file: &[u8]) -> Vec<u8> {
    let mut text = Vec::new();
    let mut inflater = Decompress::new(false);
    let mut room = vec![0; INFLATE_ROOM];
    let mut rest = file;
    while let Some(stream) = after_header(rest) {
        inflater.reset(false);
        let Some(after) = member(stream, &mut inflater, &mut room, &mut text) else {
            break;
        };
        // A zero byte after a member ends the reading, whatever follows the zero bytes.
        if after.first().is_none_or(|byte| *byte == 0) {
            break;
        }
        rest = after;
    }
    text
}

/// The deflate stream after a member's header, `None` for a header `zcat` does not read:
/// one that is cut short, opens with other bytes, names another method than deflate,
/// carries a refused flag, or whose checksum is not that of the header's bytes before it.
fn after_header(member: &[u8]) -> Option<&[u8]> {
    let fixed = member.get(..FIXED_HEADER)?;
    if !MAGICS.contains(&[fixed[0], fixed[1]]) || fixed[2] != DEFLATE {
        return None;
    }
    let flags = fixed[3];
    if flags & REFUSED_FLAGS != 0 {
        return None;
    }
    let mut at = FIXED_HEADER;
    if flags & FEXTRA != 0 {
        let length = member.get(at..at + 2)?;
        at += 2 + usize::from(u16::from_le_bytes([length[0], length[1]]));
    }
    for text in [FNAME, FCOMMENT] {
        if flags & text != 0 {
            at += memchr::memchr(0, member.get(at..)?)? + 1;
        }
    }
    if flags & FHCRC != 0 {
        let stored = member.get(at..at + 2)?;
        let mut crc = Crc::new();
        crc.update(&member[..at]);
        // The header's checksum is the low 16 bits of the 32-bit CRC of the bytes before it.
        if u16::from_le_bytes([stored[0], stored[1]]) != crc.sum() as u16 {
            return None;
        }
        at += 2;
    }
    member.get(at..)
}

/// Inflates the deflate stream that opens `stream` onto `text` and returns what follows
/// the member's trailer, `None` where the stream does not end, the trailer is cut short,
/// or its checksum or length is not the text's. The text inflated stays on `text` in
/// every case.
fn member<'a>(
    stream: &'a [u8],
    inflater: &mut Decompress,
    room: &mut [u8],
    text: &mut Vec<u8>,
) -> Option<&'a [u8]> {
    let mut crc = Crc::new();
    let mut rest = stream;
    loop {
        let (read, written) = (inflater.total_in(), inflater.total_out());
        let status = inflater.decompress(rest, room, FlushDecompress::None);
        let piece = &room[..(inflater.total_out() - written) as usize];
        crc.update(piece);
        text.extend_from_slice(piece);
        rest = &rest[(inflater.total_in() - read) as usize..];
        match status {
            Ok(Status::StreamEnd) => break,
            Ok(_) if !piece.is_empty() || inflater.total_in() != read => {}
            _ => return None,
        }
    }
    let trailer = rest.get(..TRAILER)?;
    let stored = |at: usize| {
        u32::from_le_bytes([
            trailer[at],
            trailer[at + 1],
            trailer[at + 2],
            trailer[at + 3],
        ])
    };
    // The stored length is the text's length modulo 2^32.
    (stored(0) == crc.sum() && stored(4) == inflater.total_out() as u32).then(|| &rest[TRAILER..])
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::DeflateEncoder;
    use std::io::Write as _;

    const FIRST: &[u8] = b"#h\n21\t100\t.\tA\tG\n";
    const SECOND: &[u8] = b"21\t200\t.\tC\tT\n";

    fn deflated(text: &[u8]) -> Vec<u8> {
        let mut deflate = DeflateEncoder::new(Vec::new(), flate2::Compression::default());
        deflate.write_all(text).unwrap();
        deflate.finish().unwrap()
    }

    /// One member holding `text` under a header that opens with `magic`, carries `flags`
    /// and holds `fields` after its fixed bytes.
    fn member_with(magic: [u8; 2], flags: u8, fields: &[u8], text: &[u8]) -> Vec<u8> {
        let mut crc = Crc::new();
        crc.update(text);
        let mut member = vec![magic[0], magic[1], DEFLATE, flags, 0, 0, 0, 0, 0, 0xff];
        member.extend_from_slice(fields);
        member.extend_from_slice(&deflated(text));
        member.extend_from_slice(&crc.sum().to_le_bytes());
        member.extend_from_slice(&(text.len() as u32).to_le_bytes());
        member
    }

    fn member_of(text: &[u8]) -> Vec<u8> {
        member_with(MAGICS[0], 0, b"", text)
    }

    /// The low 16 bits of the 32-bit CRC of a header of `flags` and `fields`.
    fn header_checksum(flags: u8, fields: &[u8]) -> [u8; 2] {
        let mut crc = Crc::new();
        crc.update(&[0x1f, 0x8b, DEFLATE, flags, 0, 0, 0, 0, 0, 0xff]);
        crc.update(fields);
        (crc.sum() as u16).to_le_bytes()
    }

    fn joined(parts: &[&[u8]]) -> Vec<u8> {
        parts.concat()
    }

    #[test]
    fn members_are_read_one_after_another_up_to_the_first_zero_byte_after_one() {
        let (first, second) = (member_of(FIRST), member_of(SECOND));
        let both = joined(&[FIRST, SECOND]);
        assert_eq!(text(b""), b"");
        assert_eq!(text(&first), FIRST);
        assert_eq!(text(&joined(&[&first, &second])), both);
        assert_eq!(text(&joined(&[&first, &member_of(b""), &second])), both);
        // Zero bytes after the last member.
        assert_eq!(text(&joined(&[&first, &second, &[0; 8]])), both);
        // Zero bytes between two members, one or several.
        for padding in [1, 2, 8] {
            let padded = joined(&[&first, &vec![0; padding], &second]);
            assert_eq!(text(&padded), FIRST, "{padding} zero bytes");
        }
        let later = joined(&[&first, &second, &[0; 3], &second]);
        assert_eq!(text(&later), both);
        // Zero bytes and then bytes that open no member.
        assert_eq!(text(&joined(&[&first, &[0; 4], b"junk"])), FIRST);
    }

    #[test]
    fn bytes_that_open_no_member_end_the_reading_and_keep_the_text_before_them() {
        let first = member_of(FIRST);
        for tail in [&b"junk\n"[..], b"x", b"\x1f", b"\x1f\x8b", b"\x1f\x8b\x08"] {
            assert_eq!(text(&joined(&[&first, tail])), FIRST, "{tail:?}");
        }
        // A file that opens with no member has no text, whatever follows.
        assert_eq!(text(FIRST), b"");
        assert_eq!(text(&joined(&[FIRST, &first])), b"");
        assert_eq!(text(&[0; 16]), b"");
        assert_eq!(text(&joined(&[&[0], &first[..]])), b"");
        // Another method than deflate.
        let mut stored = member_of(SECOND);
        stored[2] = 0;
        assert_eq!(text(&stored), b"");
        assert_eq!(text(&joined(&[&first, &stored])), FIRST);
    }

    #[test]
    fn a_member_opens_with_either_magic_and_its_header_fields_are_read_past() {
        let old = member_with(MAGICS[1], 0, b"", SECOND);
        assert_eq!(text(&old), SECOND);
        let both = joined(&[FIRST, SECOND]);
        assert_eq!(text(&joined(&[&member_of(FIRST), &old])), both);
        assert_eq!(text(&member_with(MAGICS[0], 1, b"", FIRST)), FIRST);
        let extra = b"\x05\x00extra";
        assert_eq!(text(&member_with(MAGICS[0], FEXTRA, extra, FIRST)), FIRST);
        assert_eq!(
            text(&member_with(MAGICS[0], FNAME, b"name\0", FIRST)),
            FIRST
        );
        assert_eq!(
            text(&member_with(MAGICS[0], FCOMMENT, b"c\0", FIRST)),
            FIRST
        );
        // A field that runs past the end of the file.
        assert_eq!(
            text(b"\x1f\x8b\x08\x08\0\0\0\0\0\xffname without its end"),
            b""
        );
        assert_eq!(text(b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x32\x00short"), b"");
    }

    #[test]
    fn a_header_checksum_is_checked_and_the_encryption_and_reserved_flags_are_refused() {
        let first = member_of(FIRST);
        let checked = |flags: u8, fields: &[u8], checksum: [u8; 2]| {
            member_with(MAGICS[0], flags, &joined(&[fields, &checksum]), SECOND)
        };
        assert_eq!(
            text(&checked(FHCRC, b"", header_checksum(FHCRC, b""))),
            SECOND
        );
        let all = FHCRC | FEXTRA | FNAME | FCOMMENT;
        let fields = b"\x03\x00abcn\0c\0";
        assert_eq!(
            text(&checked(all, fields, header_checksum(all, fields))),
            SECOND
        );
        for wrong in [
            checked(FHCRC, b"", [0x12, 0x34]),
            checked(all, fields, [0, 0]),
        ] {
            assert_eq!(text(&wrong), b"");
            assert_eq!(text(&joined(&[&first, &wrong])), FIRST);
        }
        for flag in [0x20, 0x40, 0x80] {
            let refused = member_with(MAGICS[0], flag, b"", SECOND);
            assert_eq!(text(&refused), b"", "flag {flag:#x}");
            assert_eq!(text(&joined(&[&first, &refused])), FIRST, "flag {flag:#x}");
        }
    }

    #[test]
    fn the_text_inflated_before_a_cut_or_a_trailer_of_another_text_stays() {
        let (first, second) = (member_of(FIRST), member_of(SECOND));
        let both = joined(&[FIRST, SECOND]);
        // A trailer that is cut short, the whole of it included.
        for cut in [1, 4, TRAILER] {
            assert_eq!(
                text(&first[..first.len() - cut]),
                FIRST,
                "{cut} bytes short"
            );
            let two = joined(&[&first, &second[..second.len() - cut]]);
            assert_eq!(text(&two), both, "{cut} bytes short");
        }
        // A checksum or a length that is not the text's ends the reading after the text.
        for field in [TRAILER, TRAILER / 2] {
            let mut wrong = first.clone();
            let at = wrong.len() - field;
            wrong[at] ^= 1;
            assert_eq!(text(&joined(&[&wrong, &second])), FIRST);
            assert_eq!(
                text(&joined(&[&second, &wrong, &second])),
                joined(&[SECOND, FIRST])
            );
        }
        // A header that is cut short has no text.
        for cut in [1, 2, 3, 9, FIXED_HEADER] {
            assert_eq!(text(&first[..cut]), b"", "{cut} bytes");
            assert_eq!(
                text(&joined(&[&first, &second[..cut]])),
                FIRST,
                "{cut} bytes"
            );
        }
    }

    #[test]
    fn a_stored_block_that_is_cut_short_gives_the_bytes_before_the_cut() {
        let mut deflate = DeflateEncoder::new(Vec::new(), flate2::Compression::none());
        deflate.write_all(FIRST).unwrap();
        let stream = deflate.finish().unwrap();
        // A stored block: one byte of block header, two lengths of two bytes, the text.
        assert_eq!(&stream[5..5 + FIRST.len()], FIRST);
        let mut member = vec![0x1f, 0x8b, DEFLATE, 0, 0, 0, 0, 0, 0, 0xff];
        member.extend_from_slice(&stream);
        for kept in [0, 1, 7, FIRST.len()] {
            let cut = &member[..FIXED_HEADER + 5 + kept];
            assert_eq!(text(cut), &FIRST[..kept], "{kept} bytes of text");
        }
    }
}
