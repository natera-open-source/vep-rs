// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A gzip file held in memory: the members read and the reason one that does not read
//! gives, each held to what CPython's `gzip.decompress` gives for the same bytes, and the
//! member written, held to the bytes of `zlib.compress(data, level=9, wbits=31)`.

use super::*;

fn bytes(hex: &str) -> Vec<u8> {
    let digits = hex.as_bytes().chunks(2);
    digits
        .map(|pair| u8::from_str_radix(std::str::from_utf8(pair).unwrap(), 16).unwrap())
        .collect()
}

fn hex(data: &[u8]) -> String {
    data.iter().map(|byte| format!("{byte:02x}")).collect()
}

/// One member that holds `first\n`.
const ONE: &str = "1f8b08000000000002034bcb2c2a2ee102002ab34ac706000000";

#[test]
fn every_member_is_read_and_the_zero_bytes_between_members_are_skipped() {
    assert_eq!(decompress(&bytes(ONE)).unwrap(), b"first\n");
    // A second member with an extra field, a name, a comment and a header checksum, zero
    // bytes before and after it.
    let two = format!(
        "{ONE}0000001f8b081e00000000020306004243020010006e616d65006e6f74650012342b4e4dcecf4be102007ec00f060700000000"
    );
    assert_eq!(decompress(&bytes(&two)).unwrap(), b"first\nsecond\n");
    assert_eq!(decompress(b"").unwrap(), b"");
}

#[test]
fn a_file_that_does_not_read_gives_the_reason_of_the_cpython_exception() {
    let cut_short = "Compressed file ended before the end-of-stream marker was reached";
    let cases = [
        // Bytes after a member that start no member, zero bytes before them or not.
        (format!("{ONE}7879"), "Not a gzipped file (b'xy')"),
        (format!("{ONE}78"), "Not a gzipped file (b'x')"),
        (format!("{ONE}00007827"), "Not a gzipped file (b\"x'\")"),
        (
            "1f8b07000000000002034bcb2c2a2ee102002ab34ac706000000".to_string(),
            "Unknown compression method",
        ),
        // A header that ends inside its fixed part, its extra field or its checksum.
        ("1f8b08000000".to_string(), cut_short),
        ("1f8b".to_string(), cut_short),
        ("1f8b080400000000020309".to_string(), cut_short),
        ("1f8b080400000000020309000102".to_string(), cut_short),
        ("1f8b080200000000020301".to_string(), cut_short),
        // A name without its closing zero byte ends at the end of the data, where the
        // deflate data is missing.
        ("1f8b08080000000002036e616d65".to_string(), cut_short),
        // Deflate data that ends before its last block, and a trailer of fewer than eight
        // bytes.
        ("1f8b08000000000002034bcb2c2a".to_string(), cut_short),
        ("1f8b0800000000000203".to_string(), cut_short),
        (ONE[..ONE.len() - 2].to_string(), cut_short),
        (
            "1f8b08000000000002034bcb2c2a2ee102000000000006000000".to_string(),
            "CRC check failed",
        ),
        (
            "1f8b08000000000002034bcb2c2a2ee102002ab34ac701000000".to_string(),
            "Incorrect length of data produced",
        ),
        (
            "1f8b080000000000020307cb2c2a2ee102002ab34ac706000000".to_string(),
            "Error -3 while decompressing data: invalid block type",
        ),
    ];
    for (data, reason) in cases {
        assert_eq!(
            decompress(&bytes(&data)),
            Err(ReadError::Cpython(reason.to_string())),
            "{data}"
        );
    }
}

#[test]
fn a_member_is_written_as_zlib_writes_one_at_level_9() {
    assert_eq!(
        hex(&compress(&[b'x'; 10])),
        "1f8b0800000000000203aba88001007747df8e0a000000"
    );
    assert_eq!(
        hex(&compress(b"")),
        "1f8b080000000000020303000000000000000000"
    );
    // Data that the first room for the member does not hold comes back whole.
    let long: Vec<u8> = (0..1_000_000u32)
        .map(|place| (place.wrapping_mul(2_654_435_761) >> 24) as u8)
        .collect();
    let stored = compress(&long);
    assert_eq!(stored[..10], [0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 2, 3]);
    assert_eq!(decompress(&stored).unwrap(), long);
}
