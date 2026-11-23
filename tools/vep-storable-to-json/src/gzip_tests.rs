// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Expected bytes are `Compress::Zlib::memGzip`'s on Compress::Zlib 2.102 with zlib
//! 1.2.11.

use super::mem_gzip;

// perl -MCompress::Zlib -e 'my $g = Compress::Zlib::memGzip("VEP" . ("\x01\x02\x03\x04" x 5000));
//     print length($g), " ", unpack("H*", substr($g, 0, 12)), " ", unpack("H*", substr($g, -8)), "\n"'
#[test]
fn a_repeating_block_is_the_member_memgzip_writes() {
    let mut bytes = b"VEP".to_vec();
    for _ in 0..5000 {
        bytes.extend_from_slice(&[1, 2, 3, 4]);
    }
    let member = mem_gzip(&bytes);
    assert_eq!(member.len(), 63);
    assert_eq!(
        member[..12],
        [0x1f, 0x8b, 0x08, 0, 0, 0, 0, 0, 0, 0xff, 0xed, 0xc3]
    );
    assert_eq!(
        member[member.len() - 8..],
        [0xd4, 0x57, 0xdc, 0xa7, 0x23, 0x4e, 0, 0]
    );
}

#[test]
fn an_empty_blob_is_a_member_of_an_empty_stream() {
    // The minimum header, the two bytes of an empty final block, a zero CRC and length.
    assert_eq!(
        mem_gzip(b""),
        [0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    );
}

#[test]
fn the_member_holds_the_bytes_it_was_given() {
    use std::io::Read;
    let bytes: Vec<u8> = (0..200_000u32).map(|n| (n * 7919 % 251) as u8).collect();
    let member = mem_gzip(&bytes);
    let mut back = Vec::new();
    flate2::read::GzDecoder::new(&member[..])
        .read_to_end(&mut back)
        .expect("a gzip member");
    assert_eq!(back, bytes);
}
