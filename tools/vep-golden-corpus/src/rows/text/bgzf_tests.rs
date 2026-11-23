// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A BGZF file gives the same bytes, and where it fails the same reason after the same
//! bytes, on one thread and on several.

use super::super::{ByteSource, ReadError};
use crate::TempDir;
use flate2::write::GzEncoder;
use flate2::{Compress, Compression, Crc, FlushCompress};
use std::io::Write as _;

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-bgzf-").unwrap()
}

/// One BGZF block that holds `data`.
fn block(data: &[u8]) -> Vec<u8> {
    let mut deflate = Compress::new(Compression::new(6), false);
    let mut body = Vec::with_capacity(data.len() + 1024);
    deflate
        .compress_vec(data, &mut body, FlushCompress::Finish)
        .unwrap();
    let mut crc = Crc::new();
    crc.update(data);
    let length = (18 + body.len() + 8 - 1) as u16;
    let mut out = vec![
        0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
    ];
    out.extend(length.to_le_bytes());
    out.extend(body);
    out.extend(crc.sum().to_le_bytes());
    out.extend((data.len() as u32).to_le_bytes());
    out
}

/// Bytes that compress a little: lines of a counter and letters.
fn text(lines: usize) -> Vec<u8> {
    let mut out = Vec::new();
    let mut state = 17u32;
    for line in 0..lines {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
        out.extend(format!("21\t{line}\trs{}\t{:x}\n", state >> 8, state).as_bytes());
    }
    out
}

/// A BGZF file of `data` in blocks of `size` bytes, with the empty block that ends one.
fn bgzf(data: &[u8], size: usize) -> Vec<u8> {
    let mut out: Vec<u8> = data.chunks(size).flat_map(block).collect();
    out.extend(block(b""));
    out
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

/// Every byte a file gives on `threads` threads, read `room` bytes at a time, and the
/// error that ended the reading.
fn read(path: &str, threads: usize, room: usize) -> (Vec<u8>, Option<ReadError>) {
    let mut source = ByteSource::open_gzip_on(path, threads).unwrap();
    let mut bytes = Vec::new();
    let mut into = vec![0u8; room];
    loop {
        match source.read(&mut into) {
            Ok(0) => return (bytes, None),
            Ok(count) => bytes.extend(&into[..count]),
            Err(error) => return (bytes, Some(error)),
        }
    }
}

/// Holds a file to the same reading on one thread and on several, and returns it.
fn same_on_every_count(
    directory: &TempDir,
    name: &str,
    file: &[u8],
) -> (Vec<u8>, Option<ReadError>) {
    let path = directory.path.join(name);
    std::fs::write(&path, file).unwrap();
    let path = path.to_str().unwrap();
    let one = read(path, 1, 100_000);
    for (threads, room) in [(2, 100_000), (3, 7), (16, 1 << 20), (64, 4096)] {
        let several = read(path, threads, room);
        assert!(
            several == one,
            "{name} on {threads} threads: {} bytes and {:?} against {} bytes and {:?}",
            several.0.len(),
            several.1,
            one.0.len(),
            one.1
        );
    }
    one
}

#[test]
fn blocks_are_read_in_the_files_order_on_every_count_of_threads() {
    let directory = scratch();
    let data = text(120_000);
    assert!(data.len() > 2_000_000);
    for (name, size) in [("blocks.vcf.gz", 60_000), ("small-blocks.vcf.gz", 1_000)] {
        let file = bgzf(&data, size);
        assert!(
            same_on_every_count(&directory, name, &file) == (data.clone(), None),
            "{name}"
        );
    }
    // A file without the empty block that ends one, an empty file and one empty block.
    let open: Vec<u8> = data.chunks(50_000).flat_map(block).collect();
    assert!(same_on_every_count(&directory, "open.gz", &open) == (data.clone(), None));
    assert_eq!(
        same_on_every_count(&directory, "empty.gz", b""),
        (Vec::new(), None)
    );
    assert_eq!(
        same_on_every_count(&directory, "end.gz", &block(b"")),
        (Vec::new(), None)
    );
}

#[test]
fn members_that_are_no_blocks_are_read_in_their_turn() {
    let directory = scratch();
    let data = text(30_000);
    let (first, rest) = data.split_at(200_000);
    let (second, third) = rest.split_at(150_000);
    // Blocks, then a gzip member that states no length, zero bytes, and blocks again.
    let mut file: Vec<u8> = first.chunks(40_000).flat_map(block).collect();
    file.extend(gzip(second));
    file.extend([0u8; 37]);
    file.extend(third.chunks(40_000).flat_map(block));
    file.extend([0u8; 5]);
    assert!(same_on_every_count(&directory, "mixed.gz", &file) == (data.clone(), None));
    // A gzip file that is not BGZF at all, and zero bytes before any member.
    assert!(same_on_every_count(&directory, "plain.gz", &gzip(&data)) == (data.clone(), None));
    let mut padded = vec![0u8; 3];
    padded.extend(bgzf(&data, 40_000));
    let (bytes, error) = same_on_every_count(&directory, "padded.gz", &padded);
    assert_eq!(
        (bytes.len(), error),
        (
            0,
            Some(ReadError::Cpython(
                "Not a gzipped file (b'\\x00\\x00')".to_string()
            ))
        )
    );
}

#[test]
fn a_block_that_fails_gives_the_bytes_before_it_and_the_reason() {
    let directory = scratch();
    let data = text(30_000);
    let blocks: Vec<Vec<u8>> = data.chunks(40_000).map(block).collect();
    let before: usize = 3 * 40_000;
    let with = |change: &dyn Fn(&mut Vec<u8>)| -> Vec<u8> {
        let mut changed = blocks.clone();
        change(&mut changed[3]);
        let mut file = changed.concat();
        file.extend(block(b""));
        file
    };
    let cpython = |reason: &str| Some(ReadError::Cpython(reason.to_string()));
    // A checksum that is not the block's: its bytes are given, then the reason.
    let file = with(&|block| {
        let at = block.len() - 8;
        block[at] ^= 0xff;
    });
    let (bytes, error) = same_on_every_count(&directory, "crc.gz", &file);
    assert_eq!(bytes.len(), before + 40_000);
    assert!(
        matches!(error, Some(ReadError::Cpython(reason)) if reason.starts_with("CRC check failed 0x"))
    );
    // A length that is not the block's.
    let file = with(&|block| {
        let at = block.len() - 4;
        block[at] ^= 0x01;
    });
    let (bytes, error) = same_on_every_count(&directory, "length.gz", &file);
    assert_eq!(
        (bytes.len(), error),
        (
            before + 40_000,
            cpython("Incorrect length of data produced")
        )
    );
    // Deflate data that does not inflate.
    let file = with(&|block| block[18] = 0xff);
    let (bytes, error) = same_on_every_count(&directory, "deflate.gz", &file);
    assert_eq!(
        (bytes.len(), error),
        (
            before,
            cpython("Error -3 while decompressing data: invalid block type")
        )
    );
    // A file cut inside a block, and inside a header.
    let whole = with(&|_| {});
    let fourth: usize = blocks[..3].iter().map(Vec::len).sum();
    for (name, cut) in [
        ("cut-block.gz", fourth + 500),
        ("cut-header.gz", fourth + 11),
    ] {
        let (bytes, error) = same_on_every_count(&directory, name, &whole[..cut]);
        assert!(
            bytes.len() >= before && bytes.len() < before + 40_000,
            "{name}"
        );
        assert_eq!(
            error,
            cpython("Compressed file ended before the end-of-stream marker was reached"),
            "{name}"
        );
    }
    // Bytes behind the blocks that are no member.
    let mut trailed = whole.clone();
    trailed.extend(b"junk");
    let (bytes, error) = same_on_every_count(&directory, "trailed.gz", &trailed);
    assert_eq!(
        (bytes.len(), error),
        (data.len(), cpython("Not a gzipped file (b'ju')"))
    );
}
