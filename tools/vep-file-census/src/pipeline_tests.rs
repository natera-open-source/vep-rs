// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use flate2::write::GzEncoder;
use flate2::{Compress, Compression, FlushCompress};
use std::io::Write;
use std::path::PathBuf;

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path =
            std::env::temp_dir().join(format!("vep-file-census-unit-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn file(&self, name: &str, bytes: &[u8]) -> File {
        let path = self.0.join(name);
        std::fs::write(&path, bytes).unwrap();
        File::open(path).unwrap()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Lines of every length up to a few kilobytes, about `bytes` bytes of them; the last
/// line has no newline.
fn text(bytes: usize) -> Vec<u8> {
    let mut out = Vec::new();
    let mut at = 0usize;
    while out.len() < bytes {
        at += 1;
        writeln!(out, "{at}\t{}", "x".repeat(at * 7919 % 4000)).unwrap();
    }
    out.extend_from_slice(b"the last line");
    out
}

fn block_gzip(text: &[u8], block: usize) -> Vec<u8> {
    let mut out = Vec::new();
    for chunk in text.chunks(block).chain(std::iter::once(&[][..])) {
        let mut deflated = vec![0u8; chunk.len() + 1024];
        let mut compress = Compress::new(Compression::fast(), false);
        compress
            .compress(chunk, &mut deflated, FlushCompress::Finish)
            .unwrap();
        deflated.truncate(compress.total_out() as usize);
        let length = BGZF_HEADER + deflated.len() + 8;
        out.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0]);
        out.extend_from_slice(&((length - 1) as u16).to_le_bytes());
        out.extend_from_slice(&deflated);
        out.extend_from_slice(&crc32fast::hash(chunk).to_le_bytes());
        out.extend_from_slice(&(chunk.len() as u32).to_le_bytes());
    }
    out
}

fn gzip_members(text: &[u8], members: usize) -> Vec<u8> {
    let mut out = Vec::new();
    for member in text.chunks(text.len() / members + 1) {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(member).unwrap();
        out.extend_from_slice(&encoder.finish().unwrap());
    }
    out
}

fn all_bytes(file: &File, reading: Reading, threads: usize) -> (Vec<u8>, io::Result<Finish>) {
    let mut out = Vec::new();
    let finish = pieces(
        file,
        reading,
        threads,
        |bytes| bytes,
        |bytes| {
            out.extend_from_slice(&bytes);
            Ok(())
        },
    );
    (out, finish)
}

#[test]
fn what_the_workers_make_reaches_the_caller_in_the_sources_order() {
    for threads in [1, 2, 7] {
        let mut seen = Vec::new();
        let emitted = ordered(
            threads,
            |emit| {
                let mut count = 0;
                for item in 0..500u64 {
                    if !emit(item) {
                        break;
                    }
                    count += 1;
                }
                count
            },
            |item| {
                // The later items finish first.
                std::thread::sleep(std::time::Duration::from_micros(500 - item));
                item * 2
            },
            |made| {
                seen.push(made);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(emitted, 500);
        assert_eq!(seen, (0..500u64).map(|item| item * 2).collect::<Vec<_>>());
    }
}

#[test]
fn a_caller_that_fails_stops_the_source() {
    let failure = ordered(
        3,
        |emit| (0..1_000_000u64).take_while(|item| emit(*item)).count(),
        |item| item,
        |made| {
            if made == 10 {
                Err(io::Error::other("enough"))
            } else {
                Ok(())
            }
        },
    )
    .unwrap_err();
    assert_eq!(failure.to_string(), "enough");
}

#[test]
fn a_file_is_read_whole_and_in_order_in_every_packing() {
    let scratch = Scratch::new("pieces");
    let text = text(5 << 20);
    let plain = scratch.file("plain", &text);
    let blocks = scratch.file("blocks.gz", &block_gzip(&text, 30_000));
    let members = scratch.file("members.gz", &gzip_members(&text, 3));
    // Blocks first, then a member that is no block, then blocks again.
    let third = text.len() / 3;
    let mixed = [
        block_gzip(&text[..third], 30_000),
        gzip_members(&text[third..2 * third], 1),
        block_gzip(&text[2 * third..], 30_000),
    ]
    .concat();
    let mixed = scratch.file("mixed.gz", &mixed);
    for threads in [1, 2, 6] {
        for (file, reading) in [
            (&plain, Reading::Plain { from: 0 }),
            (&plain, Reading::Auto),
            (&blocks, Reading::Gunzip),
            (&blocks, Reading::Auto),
            (&members, Reading::Gunzip),
            (&members, Reading::Auto),
            (&mixed, Reading::Gunzip),
        ] {
            let (read, finish) = all_bytes(file, reading, threads);
            assert!(read == text, "{reading:?} on {threads} thread(s)");
            assert!(!finish.unwrap().garbage);
        }
        let (read, _) = all_bytes(&plain, Reading::Plain { from: 12_345 }, threads);
        assert!(read == text[12_345..]);
    }
}

// cat blocks.gz large.gz blocks.gz | gzip -dc, where large.gz is one member under a block's
// header that inflates to more than a block holds.
#[test]
fn a_member_under_a_blocks_header_that_is_no_block_is_read_as_a_member() {
    let scratch = Scratch::new("large");
    let text = text(1 << 20);
    let large = b"21\t100\t.\tfiller\n".repeat(140_000);
    let half = text.len() / 2;
    let member = block_gzip(&large, large.len());
    assert_eq!(bgzf_block_len(&member), Some(member.len() - 28));
    let cases = [
        // Blocks before it and after it.
        (
            [
                block_gzip(&text[..half], 30_000),
                member.clone(),
                block_gzip(&text[half..], 30_000),
            ]
            .concat(),
            [&text[..half], &large[..], &text[half..]].concat(),
        ),
        // No block before it.
        (
            [member, block_gzip(&text, 30_000)].concat(),
            [&large[..], &text[..]].concat(),
        ),
    ];
    for (at, (bytes, expected)) in cases.iter().enumerate() {
        let file = scratch.file(&format!("large-{at}.gz"), bytes);
        for threads in [1, 2, 6] {
            for reading in [Reading::Gunzip, Reading::Auto] {
                let (read, finish) = all_bytes(&file, reading, threads);
                assert!(
                    read == *expected,
                    "file {at}, {reading:?} on {threads} thread(s)"
                );
                assert!(!finish.unwrap().garbage);
            }
        }
    }
}

#[test]
fn what_was_read_before_a_failure_reaches_the_caller_before_it() {
    let scratch = Scratch::new("failures");
    let text = text(3 << 20);
    let blocks = block_gzip(&text, 30_000);
    let cut = scratch.file("cut.gz", &blocks[..blocks.len() / 2]);
    let (read, finish) = all_bytes(&cut, Reading::Gunzip, 4);
    assert!(finish.is_err());
    assert!(read.len() > text.len() / 3 && read.len() < text.len() / 2 + 30_000);
    assert!(text.starts_with(&read));

    let garbage = scratch.file(
        "garbage.gz",
        &[gzip_members(&text, 2), b"garbage".to_vec()].concat(),
    );
    let (read, finish) = all_bytes(&garbage, Reading::Gunzip, 4);
    assert!(read == text);
    assert!(finish.unwrap().garbage);
    let (read, finish) = all_bytes(&garbage, Reading::Auto, 4);
    assert!(read == [text.clone(), b"garbage".to_vec()].concat());
    assert!(!finish.unwrap().garbage);

    let not_gzip = scratch.file("plain.gz", &text);
    let (read, finish) = all_bytes(&not_gzip, Reading::Gunzip, 4);
    assert!(read.is_empty() && finish.is_err());
}

#[test]
fn lines_reach_the_pass_whole_and_in_order() {
    let scratch = Scratch::new("lines");
    let text = text(5 << 20);
    let plain = scratch.file("plain", &text);
    let blocks = scratch.file("blocks.gz", &block_gzip(&text, 30_000));
    // One line longer than several pieces, between short ones.
    let long = [
        b"first\n".to_vec(),
        vec![b'y'; 7 << 20],
        b"\nlast\n".to_vec(),
    ]
    .concat();
    let long_file = scratch.file("long", &long);
    for threads in [1, 3] {
        for (file, reading, expected) in [
            (&plain, Reading::Plain { from: 0 }, &text),
            (&blocks, Reading::Gunzip, &text),
            (&long_file, Reading::Plain { from: 0 }, &long),
        ] {
            let mut joined = Vec::new();
            let mut runs = 0;
            let finish = lines(
                file,
                reading,
                threads,
                Ends::Newline,
                |data| data.to_vec(),
                |run: Vec<u8>| {
                    // Every run but the file's last ends where a line ends.
                    assert!(joined.last().is_none_or(|byte| *byte == b'\n'));
                    joined.extend_from_slice(&run);
                    runs += 1;
                    Ok(())
                },
            );
            assert!(finish.is_ok());
            assert!(joined == *expected, "{reading:?} on {threads} thread(s)");
            assert!(runs > 2);
        }
    }
}

#[test]
fn a_return_ends_a_line_where_the_pass_asks_for_it() {
    let scratch = Scratch::new("returns");
    let text: Vec<u8> = (0..300_000u32)
        .flat_map(|at| format!("row {at}{}", ["\n", "\r\n", "\r"][(at % 3) as usize]).into_bytes())
        .collect();
    let file = scratch.file("returns", &text);
    let mut rows = 0;
    lines(
        &file,
        Reading::Plain { from: 0 },
        3,
        Ends::NewlineOrReturn,
        |data| {
            assert!(matches!(data.last(), Some(b'\n' | b'\r')));
            data.split(|byte| matches!(byte, b'\n' | b'\r'))
                .filter(|line| !line.is_empty())
                .count()
        },
        |count| {
            rows += count;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(rows, 300_000);
}

#[test]
fn the_thread_bound_is_read_from_the_environment() {
    // One test sets the variable, so no other test of this binary reads it meanwhile.
    std::env::remove_var(THREADS_ENV);
    assert!(thread_bound().unwrap() >= 1);
    std::env::set_var(THREADS_ENV, "3");
    assert_eq!(thread_bound(), Ok(3));
    for refused in ["0", "-2", "many", "1.5"] {
        std::env::set_var(THREADS_ENV, refused);
        assert!(thread_bound().is_err(), "{refused}");
    }
    std::env::set_var(THREADS_ENV, "");
    assert!(thread_bound().unwrap() >= 1);
    std::env::remove_var(THREADS_ENV);
}
