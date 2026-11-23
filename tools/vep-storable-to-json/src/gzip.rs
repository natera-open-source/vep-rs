// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `Compress::Zlib::memGzip`: the gzip member the conversion makes of a prediction matrix
//! that is not gzipped yet.
//!
//! `memGzip` runs zlib's deflate at level 9 with memory level 9 and no zlib header. Its
//! bytes are zlib's own choice of matches and block boundaries, which another deflate
//! implementation does not repeat, so this module calls zlib itself: the stock library,
//! built from the source the `libz-sys` crate carries.

use std::mem::MaybeUninit;
use std::os::raw::c_int;

use libz_sys as z;

/// The header `memGzip` writes: deflate, no flags, no time, no extra flags, and the
/// operating system byte 255.
const MINIMUM_HEADER: [u8; 10] = [0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
/// `memGzip`'s deflate settings: `Z_BEST_COMPRESSION`, `-MAX_WBITS` (a raw stream) and
/// `MAX_MEM_LEVEL`.
const LEVEL: c_int = 9;
const WINDOW_BITS: c_int = -15;
const MEMORY_LEVEL: c_int = 9;
/// The most input handed to one `deflate` call; zlib counts it in 32 bits.
const STEP: usize = 1 << 30;

/// The gzip member of `bytes`: the minimum header, the raw deflate stream, the CRC32
/// checksum and the length modulo 2^32.
pub fn mem_gzip(bytes: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(bytes.len() / 2 + 64);
    out.extend_from_slice(&MINIMUM_HEADER);
    let mut crc = 0;
    let mut stream = MaybeUninit::<z::z_stream>::zeroed();
    // SAFETY: zlib's documented protocol. `deflateInit2_` takes a zeroed stream, whose
    // null `zalloc` and `zfree` it replaces with its own allocator, so the stream is
    // read as a `z_stream` only after it returned `Z_OK`. Every `deflate` call is given
    // input and output pointers into live buffers with their exact lengths: `next_in`
    // into `bytes`, which zlib only reads, and `next_out` into the spare capacity of
    // `out`, whose length grows by what zlib reports it wrote before the buffer is
    // touched again. `deflateEnd` frees the stream's state once, on every path past
    // the successful init.
    unsafe {
        let status = z::deflateInit2_(
            stream.as_mut_ptr(),
            LEVEL,
            z::Z_DEFLATED,
            WINDOW_BITS,
            MEMORY_LEVEL,
            z::Z_DEFAULT_STRATEGY,
            z::zlibVersion(),
            size_of::<z::z_stream>() as c_int,
        );
        assert_eq!(status, z::Z_OK, "zlib deflateInit2");
        let stream = stream.assume_init_mut();
        let mut rest = bytes;
        loop {
            let (step, later) = rest.split_at(rest.len().min(STEP));
            let flush = if later.is_empty() {
                z::Z_FINISH
            } else {
                z::Z_NO_FLUSH
            };
            crc = z::crc32(crc, step.as_ptr(), step.len() as u32);
            stream.next_in = step.as_ptr().cast_mut();
            stream.avail_in = step.len() as u32;
            let finished = loop {
                out.reserve(64 * 1024);
                let spare = out.spare_capacity_mut();
                let room = spare.len().min(STEP);
                stream.next_out = spare.as_mut_ptr().cast();
                stream.avail_out = room as u32;
                let status = z::deflate(stream, flush);
                let written = room - stream.avail_out as usize;
                out.set_len(out.len() + written);
                match status {
                    z::Z_STREAM_END => break true,
                    // More output room is wanted, or this step's input is used up.
                    z::Z_OK | z::Z_BUF_ERROR if stream.avail_in > 0 || stream.avail_out == 0 => {}
                    z::Z_OK | z::Z_BUF_ERROR => break false,
                    other => panic!("zlib deflate returned {other}"),
                }
            };
            rest = later;
            if finished {
                break;
            }
        }
        z::deflateEnd(stream);
    }
    out.extend_from_slice(&(crc as u32).to_le_bytes());
    out.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
    out
}

#[cfg(test)]
#[path = "gzip_tests.rs"]
mod tests;
