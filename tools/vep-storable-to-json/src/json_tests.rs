// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Expected texts are JSON::XS 4.04's on perl 5.34.0:
//! `JSON->new->utf8->canonical->encode([$string])` for the strings and
//! `JSON->new->canonical->pretty->encode(...)` for the document.

use std::collections::BTreeMap;

use super::{print_bytes, write_base64, write_pretty, write_string, PerlStr, Pretty};

fn encoded(bytes: &[u8], utf8: bool) -> Vec<u8> {
    let mut out = Vec::new();
    write_string(PerlStr { bytes, utf8 }, &mut out).expect("a string JSON::XS writes");
    out
}

#[test]
fn a_byte_string_is_escaped_as_json_xs_escapes_it() {
    assert_eq!(encoded(b"plain", false), b"\"plain\"");
    // A slash is written as itself.
    assert_eq!(encoded(b"q\"b\\s/", false), b"\"q\\\"b\\\\s/\"");
    assert_eq!(
        encoded(b"tab\tnl\ncr\rbs\x08ff\x0c", false),
        b"\"tab\\tnl\\ncr\\rbs\\bff\\f\""
    );
    // Other control characters take \u00xx in lower case; DEL is written as itself.
    assert_eq!(
        encoded(b"ctl\x00\x01\x1f\x7f", false),
        b"\"ctl\\u0000\\u0001\\u001f\x7f\""
    );
}

#[test]
fn a_byte_above_127_is_the_character_of_that_number_in_utf8() {
    assert_eq!(
        encoded(b"hi\xe9\xff", false),
        "\"hi\u{e9}\u{ff}\"".as_bytes()
    );
    // The same characters under the UTF-8 flag are already UTF-8.
    assert_eq!(
        encoded("caf\u{e9}".as_bytes(), true),
        "\"caf\u{e9}\"".as_bytes()
    );
    assert_eq!(
        encoded("\u{263a}wide".as_bytes(), true),
        "\"\u{263a}wide\"".as_bytes()
    );
    assert_eq!(
        encoded("\u{1F600}".as_bytes(), true),
        "\"\u{1F600}\"".as_bytes()
    );
}

#[test]
fn bytes_under_the_utf8_flag_that_are_not_utf8_are_refused_with_json_xs_text() {
    let mut out = Vec::new();
    let error = write_string(
        PerlStr {
            bytes: b"ab\xff\xfecd",
            utf8: true,
        },
        &mut out,
    )
    .expect_err("malformed");
    assert_eq!(
        error.0,
        b"malformed or illegal unicode character in string [\xff\xfecd], cannot convert to JSON"
    );
    // An overlong form and a sequence cut short are malformed too.
    for bytes in [&b"\xc0\x80"[..], b"\xe2\x98", b"\xe0\x80\x80"] {
        assert!(write_string(PerlStr { bytes, utf8: true }, &mut out).is_err());
    }
}

#[test]
fn characters_and_length_follow_the_utf8_flag() {
    let latin = PerlStr::bytes(b"caf\xe9");
    assert_eq!(latin.characters().as_ref(), "caf\u{e9}".as_bytes());
    assert_eq!(latin.char_len(), 4);
    let flagged = PerlStr {
        bytes: "caf\u{e9}\u{263a}".as_bytes(),
        utf8: true,
    };
    assert_eq!(
        flagged.characters().as_ref(),
        "caf\u{e9}\u{263a}".as_bytes()
    );
    assert_eq!(flagged.char_len(), 5);
}

// {a => [], b => {}, c => [1, [2]], d => {x => JSON::true, y => "s\xe9"}, e => "z"}
#[test]
fn the_pretty_document_is_indented_by_three_with_spaces_around_the_colon() {
    let object = |members: Vec<(&str, Pretty)>| {
        Pretty::Object(
            members
                .into_iter()
                .map(|(key, value)| (key.as_bytes().to_vec(), value))
                .collect::<BTreeMap<_, _>>(),
        )
    };
    let document = object(vec![
        ("a", Pretty::Array(Vec::new())),
        ("b", object(Vec::new())),
        (
            "c",
            Pretty::Array(vec![
                Pretty::Integer(1),
                Pretty::Array(vec![Pretty::Integer(2)]),
            ]),
        ),
        (
            "d",
            object(vec![
                ("x", Pretty::True),
                ("y", Pretty::String("s\u{e9}".as_bytes().to_vec())),
            ]),
        ),
        ("e", Pretty::String(b"z".to_vec())),
    ]);
    let text = write_pretty(&document).expect("a document");
    let expected = "{\n   \"a\" : [],\n   \"b\" : {},\n   \"c\" : [\n      1,\n      [\n         2\n      ]\n   ],\n   \"d\" : {\n      \"x\" : true,\n      \"y\" : \"s\u{e9}\"\n   },\n   \"e\" : \"z\"\n}\n";
    assert_eq!(text, expected.as_bytes());
    // Perl prints the document one byte a character while every character is below 256.
    let printed = print_bytes(&text);
    assert!(printed.windows(3).any(|window| window == b"s\xe9\""));
    assert_eq!(printed.len(), text.len() - 1);
}

#[test]
fn a_character_past_255_makes_perl_print_the_text_as_utf8() {
    let text = "caf\u{e9} \u{263a}".as_bytes();
    assert_eq!(print_bytes(text), text);
    assert_eq!(print_bytes("caf\u{e9}".as_bytes()), b"caf\xe9");
    assert_eq!(print_bytes(b"plain"), b"plain");
}

#[test]
fn base64_is_padded_and_unbroken() {
    let encode = |bytes: &[u8]| {
        let mut out = Vec::new();
        write_base64(bytes, &mut out);
        out
    };
    assert_eq!(encode(b""), b"");
    assert_eq!(encode(b"V"), b"Vg==");
    assert_eq!(encode(b"VE"), b"VkU=");
    assert_eq!(encode(b"VEP"), b"VkVQ");
    assert_eq!(encode(&[0xfb, 0xff, 0xfe, 0x00]), b"+//+AA==");
    assert_eq!(encode(&[0u8; 60]).len(), 80);
}
