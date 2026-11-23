// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reader and the writer against `JSON::PP` 4.06 itself.
//!
//! `tests/fixtures/jsonpp/cases.txt` holds what the module does, every text in
//! hexadecimal (`-` for the empty text):
//!
//! - `d <ok|die> <text>`: `decode_json` takes the text or dies on it;
//! - `s <text between a string's quotes> <the characters it decodes to, as UTF-8>`;
//! - `e <a string's characters as UTF-8> <what encode_json writes for it, quotes included>`.

use super::{decode, string_text, write_string};

fn bytes(hex: &str) -> Vec<u8> {
    if hex == "-" {
        return Vec::new();
    }
    (0..hex.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&hex[at..at + 2], 16).expect("hexadecimal"))
        .collect()
}

#[test]
fn texts_are_taken_and_refused_as_the_module_does() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/jsonpp/cases.txt"
    );
    let cases = std::fs::read_to_string(path).expect("the cases");
    let (mut taken, mut refused, mut strings, mut written) = (0, 0, 0, 0);
    for line in cases.lines() {
        let columns: Vec<&str> = line.split(' ').collect();
        match columns.as_slice() {
            ["d", verdict, text] => {
                let text = bytes(text);
                let shown = String::from_utf8_lossy(&text[..text.len().min(60)]).into_owned();
                assert_eq!(decode(&text).is_ok(), *verdict == "ok", "{shown:?}");
                if *verdict == "ok" {
                    taken += 1;
                } else {
                    refused += 1;
                }
            }
            ["s", raw, characters] => {
                assert_eq!(string_text(&bytes(raw)), bytes(characters), "{line}");
                strings += 1;
            }
            ["e", characters, json] => {
                let mut out = Vec::new();
                write_string(&bytes(characters), &mut out);
                assert_eq!(out, bytes(json), "{line}");
                written += 1;
            }
            _ => panic!("a line of another shape: {line:?}"),
        }
    }
    assert_eq!((taken, refused, strings, written), (41, 76, 9, 7));
}
