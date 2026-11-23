// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;

fn hex_of(algorithm: Algorithm, bytes: &[u8]) -> String {
    let mut state = State::new(algorithm);
    state.update(bytes);
    state.hex()
}

// printf abc | sha256sum; printf abc | md5sum; printf abc | xxh128sum; the same of nothing
#[test]
fn the_three_digests_are_those_of_the_digest_programs() {
    assert_eq!(
        hex_of(Algorithm::Sha256, b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    assert_eq!(
        hex_of(Algorithm::Md5, b"abc"),
        "900150983cd24fb0d6963f7d28e17f72"
    );
    assert_eq!(
        hex_of(Algorithm::Xxh128, b"abc"),
        "06b05ab6733a618578af5f94892f3950"
    );
    assert_eq!(
        hex_of(Algorithm::Sha256, b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        hex_of(Algorithm::Md5, b""),
        "d41d8cd98f00b204e9800998ecf8427e"
    );
    assert_eq!(
        hex_of(Algorithm::Xxh128, b""),
        "99aa06d3014798d86001c324468d497f"
    );
}

#[test]
fn a_digest_does_not_depend_on_how_the_bytes_arrive() {
    let bytes: Vec<u8> = (0..200_000u32).map(|at| (at * 31 % 251) as u8).collect();
    for algorithm in [Algorithm::Sha256, Algorithm::Md5, Algorithm::Xxh128] {
        let mut state = State::new(algorithm);
        for piece in bytes.chunks(4099) {
            state.update(piece);
        }
        assert_eq!(state.hex(), hex_of(algorithm, &bytes));
    }
    let digests = digest_stream(
        &bytes[..],
        &[Algorithm::Sha256, Algorithm::Md5, Algorithm::Xxh128],
    )
    .unwrap();
    assert_eq!(
        digests,
        [
            hex_of(Algorithm::Sha256, &bytes),
            hex_of(Algorithm::Md5, &bytes),
            hex_of(Algorithm::Xxh128, &bytes)
        ]
    );
}

// sha256sum 'a b' 'back\slash' $'new\nline'
#[test]
fn a_name_with_a_backslash_or_a_newline_is_escaped_and_its_line_starts_with_a_backslash() {
    assert_eq!(line("00", b"a b", true), b"00  a b\n");
    assert_eq!(line("00", b"back\\slash", true), b"\\00  back\\\\slash\n");
    assert_eq!(line("00", b"new\nline", true), b"\\00  new\\nline\n");
    assert_eq!(line("00", b"re\rturn", true), b"00  re\rturn\n");
    // xxh128sum writes the name as it is.
    assert_eq!(line("00", b"back\\slash", false), b"00  back\\slash\n");
}

// printf '##a\nb\n## c\n\n##\nd' | grep -a -v '^##'
#[test]
fn lines_that_start_with_the_text_are_dropped_and_the_last_line_gains_its_newline() {
    let kept =
        |data: &[u8], text: &[u8]| String::from_utf8(without_lines_starting(data, text)).unwrap();
    assert_eq!(kept(b"##a\nb\n## c\n\n##\nd", b"##"), "b\n\nd\n");
    assert_eq!(kept(b"b\n##a", b"##"), "b\n");
    assert_eq!(kept(b"##a\n##b\n", b"##"), "");
    assert_eq!(kept(b"", b"##"), "");
    assert_eq!(kept(b"a ##\n#b\n", b"##"), "a ##\n#b\n");
    assert_eq!(
        kept(
            b"## Output produced at x\n## other\nrow\n",
            b"## Output produced at"
        ),
        "## other\nrow\n"
    );
    assert_eq!(kept(b"a\nb\n", b""), "");
}

#[test]
fn the_leading_lines_end_at_the_first_line_that_starts_otherwise() {
    let dir = std::env::temp_dir().join(format!(
        "vep-file-census-unit-{}-leading",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("file");
    for (bytes, text, expected) in [
        (&b"##a\n##b\n#c\n##d\n"[..], &b"##"[..], 8),
        (b"##a\n##b", b"##", 7),
        (b"#c\n##d\n", b"##", 0),
        (b"", b"##", 0),
        (b"#", b"##", 0),
        (b"##", b"##", 2),
        (b"a\nb\n", b"", 4),
    ] {
        std::fs::write(&path, bytes).unwrap();
        assert_eq!(
            leading_lines(&File::open(&path).unwrap(), text).unwrap(),
            expected
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
