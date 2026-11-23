// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A gzip file held whole in memory: read as CPython's `gzip.decompress` reads one, with the
//! reason one that does not read gives in the words of the CPython exception, and written as
//! `gzip.compress` writes one at level 9 without a time stamp.

use crate::rows::text::{corrupt, py_repr_bytes, ReadError, CUT_SHORT};
use flate2::{Compress, Compression, Crc, Decompress, FlushCompress, FlushDecompress, Status};

/// The two bytes a gzip member starts with.
pub const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

fn cut_short() -> ReadError {
    ReadError::Cpython(CUT_SHORT.to_string())
}

/// The deflate data of the member `data` starts with, `None` for no data: the fixed header,
/// then the extra field, the name, the comment and the header checksum its flags name.
fn after_header(data: &[u8]) -> Result<Option<&[u8]>, ReadError> {
    if data.is_empty() {
        return Ok(None);
    }
    let magic = &data[..data.len().min(2)];
    if magic != GZIP_MAGIC {
        return Err(ReadError::Cpython(format!(
            "Not a gzipped file ({})",
            py_repr_bytes(magic)
        )));
    }
    let fixed = data.get(2..10).ok_or_else(cut_short)?;
    let (method, flags) = (fixed[0], fixed[1]);
    if method != 8 {
        return Err(ReadError::Cpython("Unknown compression method".to_string()));
    }
    let mut rest = &data[10..];
    if flags & 4 != 0 {
        let length = rest.get(..2).ok_or_else(cut_short)?;
        let length = usize::from(u16::from_le_bytes([length[0], length[1]]));
        rest = rest.get(2 + length..).ok_or_else(cut_short)?;
    }
    for text in [8, 16] {
        if flags & text != 0 {
            // A zero byte ends the text, and so does the end of the data.
            rest = match memchr::memchr(0, rest) {
                Some(end) => &rest[end + 1..],
                None => &[],
            };
        }
    }
    if flags & 2 != 0 {
        rest = rest.get(2..).ok_or_else(cut_short)?;
    }
    Ok(Some(rest))
}

/// Inflates the deflate stream `body` starts with onto `into` and returns the bytes after
/// the stream's end, `None` for a body that ends inside the stream.
fn inflate<'d>(body: &'d [u8], into: &mut Vec<u8>) -> Result<Option<&'d [u8]>, ReadError> {
    let mut inflate = Decompress::new(false);
    loop {
        if into.capacity() == into.len() {
            into.reserve(into.len().max(64 * 1024));
        }
        let (read, written) = (inflate.total_in(), inflate.total_out());
        let status = inflate
            .decompress_vec(&body[read as usize..], into, FlushDecompress::None)
            .map_err(|error| corrupt(&error))?;
        if status == Status::StreamEnd {
            return Ok(Some(&body[inflate.total_in() as usize..]));
        }
        let stalled = inflate.total_in() == read && inflate.total_out() == written;
        if stalled && into.capacity() > into.len() {
            return Ok(None);
        }
    }
}

/// `gzip.decompress(data)`: every member inflated in turn and held to its checksum and its
/// length, the zero bytes between members skipped.
pub fn decompress(data: &[u8]) -> Result<Vec<u8>, ReadError> {
    let mut whole = Vec::new();
    let mut rest = data;
    while let Some(body) = after_header(rest)? {
        let start = whole.len();
        let unused = inflate(body, &mut whole)?.ok_or_else(cut_short)?;
        let trailer = unused.get(..8).ok_or_else(cut_short)?;
        let member = &whole[start..];
        let mut crc = Crc::new();
        crc.update(member);
        let stored = |bytes: &[u8]| u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
        if stored(&trailer[..4]) != crc.sum() {
            return Err(ReadError::Cpython("CRC check failed".to_string()));
        }
        if stored(&trailer[4..]) != member.len() as u32 {
            return Err(ReadError::Cpython(
                "Incorrect length of data produced".to_string(),
            ));
        }
        let padding = unused[8..].iter().take_while(|byte| **byte == 0).count();
        rest = &unused[8 + padding..];
    }
    Ok(whole)
}

/// `gzip.compress(data, compresslevel=9, mtime=0)`: one member at level 9 whose header zlib
/// writes, with no name and no time stamp.
pub fn compress(data: &[u8]) -> Vec<u8> {
    let mut deflate = Compress::new_gzip(Compression::best(), 15);
    let mut stored = Vec::with_capacity(data.len() / 2 + 64);
    loop {
        let read = deflate.total_in() as usize;
        let status = deflate
            .compress_vec(&data[read..], &mut stored, FlushCompress::Finish)
            .expect("a deflate stream of bytes held in memory");
        if status == Status::StreamEnd {
            return stored;
        }
        stored.reserve(stored.capacity().max(4096));
    }
}

#[cfg(test)]
#[path = "gzip_tests.rs"]
mod tests;
