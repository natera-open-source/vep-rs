// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

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
/// name, a comment and the header checksum.
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
        let checksum = crc32fast::hash(&out) & 0xffff;
        out.extend_from_slice(&(checksum as u16).to_le_bytes());
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

/// What the reader reads of `file` through a buffer of `capacity` bytes, whether the
/// reading failed, and whether it ended in trailing garbage.
fn read_all(file: &[u8], capacity: usize, pass_through: bool) -> (Vec<u8>, bool, bool) {
    let mut out = Vec::new();
    let mut members = Members::new(io::BufReader::with_capacity(capacity, file), pass_through);
    let failed = members.read_to_end(&mut out).is_err();
    (out, failed, members.ended_in_garbage())
}

const CAPACITIES: [usize; 6] = [1, 2, 3, 5, 64, 4096];

// cat one.gz two.gz | gzip -dc
#[test]
fn members_are_read_one_after_another() {
    let files: [(Vec<u8>, &[u8]); 5] = [
        (member(b""), b""),
        (member(b"one\n"), b"one\n"),
        ([member(b"one\n"), member(b"two\n")].concat(), b"one\ntwo\n"),
        (
            [block(b"one\n"), block(b""), member(b"two\n"), block(b"")].concat(),
            b"one\ntwo\n",
        ),
        (
            [
                member_with(FEXTRA | FNAME | FCOMMENT | FHCRC, b"one\n"),
                member_with(FNAME, b"two\n"),
                member_with(FHCRC, b"three\n"),
            ]
            .concat(),
            b"one\ntwo\nthree\n",
        ),
    ];
    for (file, expected) in &files {
        for capacity in CAPACITIES {
            for pass_through in [false, true] {
                assert_eq!(
                    read_all(file, capacity, pass_through),
                    (expected.to_vec(), false, false),
                    "capacity {capacity}, {} bytes",
                    file.len()
                );
            }
        }
    }
}

// { cat one.gz; head -c 512 /dev/zero; } | gzip -dc: status 0
#[test]
fn zero_bytes_to_the_end_of_the_file_are_padding() {
    for zeros in [1, 2, 9, 512] {
        let file = [member(b"one\n"), vec![0; zeros]].concat();
        for capacity in CAPACITIES {
            assert_eq!(
                read_all(&file, capacity, false),
                (b"one\n".to_vec(), false, false)
            );
        }
    }
}

// { cat one.gz; printf garbage; } | gzip -dc: "trailing garbage ignored", status 2
#[test]
fn other_bytes_after_a_member_are_trailing_garbage_and_end_the_reading() {
    let good = member(b"one\n");
    let files = [
        [good.clone(), b"garbage".to_vec()].concat(),
        [good.clone(), vec![0; 9], b"x".to_vec()].concat(),
        // Zero bytes, then a member: the padding is not skipped to reach it.
        [good.clone(), vec![0; 2], member(b"two\n")].concat(),
        [good.clone(), vec![0x1f, 0x00]].concat(),
    ];
    for file in &files {
        for capacity in CAPACITIES {
            assert_eq!(
                read_all(file, capacity, false),
                (b"one\n".to_vec(), false, true),
                "capacity {capacity}, {} bytes",
                file.len()
            );
        }
    }
}

// { cat one.gz; printf x; } | gzip -dc: "unexpected end of file", status 1
#[test]
fn one_byte_after_a_member_that_is_not_zero_is_a_file_that_ends_early() {
    for last in [b'x', 0x1f] {
        let file = [member(b"one\n"), vec![last]].concat();
        for capacity in CAPACITIES {
            assert_eq!(
                read_all(&file, capacity, false),
                (b"one\n".to_vec(), true, false)
            );
        }
    }
}

// gzip -dc on a file that is no gzip file, a zero-byte file among them: status 1
#[test]
fn a_file_that_starts_with_no_member_is_no_gzip_file() {
    for file in [
        &b""[..],
        b"#plain text\n",
        &[0u8; 4],
        &[0x1f],
        &[0x1f, 0x8c, 8, 0],
    ] {
        for capacity in CAPACITIES {
            assert_eq!(read_all(file, capacity, false), (Vec::new(), true, false));
        }
    }
}

// zcat -f: what is no member is written as it is, at the start of the file and after a member
#[test]
fn under_pass_through_bytes_that_start_no_member_are_read_as_they_are() {
    let good = member(b"one\n");
    let files: [(Vec<u8>, Vec<u8>); 6] = [
        (Vec::new(), Vec::new()),
        (b"#plain text\n".to_vec(), b"#plain text\n".to_vec()),
        (vec![0x1f], vec![0x1f]),
        (vec![0x1f, b'x', b'y'], vec![0x1f, b'x', b'y']),
        (
            [good.clone(), b"garbage".to_vec()].concat(),
            b"one\ngarbage".to_vec(),
        ),
        (
            [good.clone(), vec![0x1f, 0, 7]].concat(),
            [b"one\n".to_vec(), vec![0x1f, 0, 7]].concat(),
        ),
    ];
    for (file, expected) in &files {
        for capacity in CAPACITIES {
            assert_eq!(
                read_all(file, capacity, true),
                (expected.clone(), false, false),
                "capacity {capacity}, {} bytes",
                file.len()
            );
        }
    }
}

#[test]
fn what_a_member_holds_is_read_before_its_failure() {
    let good = member(b"one\ntwo\n");
    let mut wrong_crc = good.clone();
    let at = wrong_crc.len() - 8;
    wrong_crc[at] ^= 1;
    let mut wrong_length = good.clone();
    let at = wrong_length.len() - 4;
    wrong_length[at] ^= 1;
    for file in [wrong_crc, wrong_length, good[..good.len() - 3].to_vec()] {
        for capacity in CAPACITIES {
            assert_eq!(
                read_all(&file, capacity, false),
                (b"one\ntwo\n".to_vec(), true, false)
            );
        }
    }
}

#[test]
fn a_header_the_reader_cannot_take_fails_the_member() {
    let good = member_with(FHCRC, b"one\n");
    let mut wrong_method = member(b"one\n");
    wrong_method[2] = 7;
    let mut wrong_header_checksum = good.clone();
    wrong_header_checksum[10] ^= 1;
    let files = [
        wrong_method,
        wrong_header_checksum,
        member_with(ENCRYPTED, b"one\n"),
        member_with(0x40, b"one\n"),
        good[..7].to_vec(),
        // A name the file ends inside.
        vec![0x1f, 0x8b, 8, FNAME, 0, 0, 0, 0, 0, 255, b'n', b'a'],
    ];
    for file in &files {
        for capacity in CAPACITIES {
            assert_eq!(read_all(file, capacity, false), (Vec::new(), true, false));
        }
    }
}

#[test]
fn reading_goes_on_after_a_member_as_it_does_inside_a_file() {
    let rest = member(b"two\n");
    let mut out = Vec::new();
    Members::after_member(io::BufReader::new(&rest[..]), false)
        .read_to_end(&mut out)
        .unwrap();
    assert_eq!(out, b"two\n");
    // Nothing after a member is the end of the file, where nothing at its start is none.
    let mut out = Vec::new();
    Members::after_member(io::BufReader::new(&b""[..]), false)
        .read_to_end(&mut out)
        .unwrap();
    assert!(out.is_empty());
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
    inflate_blocks(&batch, &mut Decompress::new(false), &mut out).unwrap();
    assert_eq!(&out[..4], b"kept");
    assert_eq!(out[4..], read_all(&batch, 4096, false).0);
    assert_eq!(bgzf_block_len(&batch), Some(block(b"one\n").len()));
}

#[test]
fn what_the_blocks_before_a_failing_one_hold_stays() {
    let good = block(b"one\ntwo\n");
    let mut wrong_crc = good.clone();
    let at = wrong_crc.len() - 8;
    wrong_crc[at] ^= 1;
    let mut wrong_length = good.clone();
    let at = wrong_length.len() - 4;
    wrong_length[at] ^= 1;
    for (failing, held) in [
        (wrong_crc, &b"first\none\ntwo\n"[..]),
        (wrong_length, b"first\none\ntwo\n"),
        (good[..good.len() - 1].to_vec(), b"first\n"),
        (member(b"x"), b"first\n"),
    ] {
        let batch = [block(b"first\n"), failing].concat();
        let mut out = Vec::new();
        assert!(inflate_blocks(&batch, &mut Decompress::new(false), &mut out).is_err());
        assert_eq!(out, held);
    }
    for header in [&member(b"one\n")[..], &good[..17], &[0u8; 18][..]] {
        assert_eq!(bgzf_block_len(header), None);
    }
    // The shortest block a header can announce is its header and its trailer.
    let mut short = good.clone();
    short[16..18].copy_from_slice(&24u16.to_le_bytes());
    assert_eq!(bgzf_block_len(&short), None);
}

// A member that gzip -dc reads whole, under the header of a block.
#[test]
fn a_member_that_announces_more_than_a_block_holds_is_no_block() {
    let data = b"21\t100\t.\tfiller\n".repeat(70_000);
    assert!(data.len() > BLOCK_LIMIT);
    let large = block(&data);
    assert_eq!(bgzf_block_len(&large), Some(large.len()));
    assert!(!bgzf_block_in_bound(&large));
    assert!(bgzf_block_in_bound(&block(&data[..BLOCK_LIMIT])));
    assert!(bgzf_block_in_bound(&block(b"")));
    assert!(!bgzf_block_in_bound(&[0, 0, 0]));
    // The block reader refuses it, and the member reader reads every byte of it.
    let mut out = Vec::new();
    assert!(inflate_blocks(&large, &mut Decompress::new(false), &mut out).is_err());
    for capacity in CAPACITIES {
        let (read, failed, garbage) = read_all(&large, capacity, false);
        assert!(read == data && !failed && !garbage, "capacity {capacity}");
    }
}
