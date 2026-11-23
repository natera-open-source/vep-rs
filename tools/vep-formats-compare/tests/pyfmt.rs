// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `pyfmt` against the answers of the interpreter.
//!
//! Each file under `tests/data/pyfmt` holds `case<TAB>expected` lines written by CPython.
//! A text is its UTF-8 bytes, each as two lowercase hex digits, joined by `_`; a float is
//! the eight bytes of its bits, written the same way; a whole number of any length and the
//! text of a float are texts; the parts of a case are joined by `~`. The ignored test reads
//! the files of the same names, with many more lines, under the `pyfmt` folder of the
//! directory `VEP_FORMATS_COMPARE_PYDIFF_DIR` names.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use vep_formats_compare::pyfmt;

fn unhex(hex: &str) -> Vec<u8> {
    if hex.is_empty() {
        return Vec::new();
    }
    hex.split('_')
        .map(|byte| {
            assert_eq!(byte.len(), 2, "two hex digits a byte in {hex}");
            // `from_str_radix` also takes a sign and a capital digit.
            assert!(
                byte.bytes()
                    .all(|digit| matches!(digit, b'0'..=b'9' | b'a'..=b'f')),
                "lowercase hex digits in {hex}"
            );
            u8::from_str_radix(byte, 16).expect("two hex digits")
        })
        .collect()
}

fn text(hex: &str) -> String {
    String::from_utf8(unhex(hex)).expect("the hex of UTF-8 bytes")
}

fn hex_of(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(3 * bytes.len());
    for (at, byte) in bytes.iter().enumerate() {
        if at > 0 {
            out.push('_');
        }
        let _ = write!(out, "{byte:02x}");
    }
    out
}

fn hex(text: &str) -> String {
    hex_of(text.as_bytes())
}

fn float(bits: &str) -> f64 {
    let bytes: [u8; 8] = unhex(bits).try_into().expect("the eight bytes of a float");
    f64::from_bits(u64::from_be_bytes(bytes))
}

fn bits(value: f64) -> String {
    hex_of(&value.to_bits().to_be_bytes())
}

fn hexes(parts: &[&str]) -> String {
    parts
        .iter()
        .map(|part| hex(part))
        .collect::<Vec<_>>()
        .join(",")
}

/// A list of texts: its length, a space, the hex of each joined by commas.
fn items(parts: &[&str]) -> String {
    format!("{} {}", parts.len(), hexes(parts))
}

fn flags(text: &str) -> String {
    [
        pyfmt::is_digit(text),
        pyfmt::is_decimal(text),
        pyfmt::is_numeric(text),
        pyfmt::is_alpha(text),
        pyfmt::is_alnum(text),
        pyfmt::is_space(text),
        pyfmt::is_printable(text),
        pyfmt::is_lower(text),
        pyfmt::is_upper(text),
        pyfmt::is_ascii(text),
    ]
    .iter()
    .map(|flag| if *flag { '1' } else { '0' })
    .collect()
}

fn bound(text: &str) -> Option<i64> {
    (text != "-").then(|| text.parse().expect("a slice bound"))
}

fn int_answer(case: &str) -> String {
    let text = text(case);
    let checked = pyfmt::int_checked(&text);
    let whole = checked.as_ref().ok().cloned();
    assert_eq!(pyfmt::int_of_big(&text), whole, "int_of_big of {text:?}");
    assert_eq!(
        pyfmt::int_of(&text),
        whole.and_then(|number| number.parse().ok()),
        "int_of of {text:?}"
    );
    match checked {
        Ok(number) => format!("ok {}", hex(&number)),
        Err(error) => format!("err {}", hex(&error.to_string())),
    }
}

fn float_answer(case: &str) -> String {
    let text = text(case);
    match pyfmt::float_of(&text) {
        Some(value) if value.is_nan() => "ok nan".to_string(),
        Some(value) => format!("ok {}", bits(value)),
        None => format!("err {}", hex(&pyfmt::float_error(&text))),
    }
}

fn int_of_f64_answer(case: &str) -> String {
    match pyfmt::int_of_f64(float(case)) {
        Ok(number) => format!("ok {}", hex(&number)),
        Err(error) => {
            let class = match error {
                pyfmt::FloatToIntError::Nan => "ValueError",
                pyfmt::FloatToIntError::Infinity => "OverflowError",
            };
            format!("err {class} {}", hex(&error.to_string()))
        }
    }
}

fn utf8_answer(case: &str) -> String {
    let bytes = unhex(case);
    let strict = match pyfmt::decode_utf8(&bytes) {
        Ok(_) => "ok".to_string(),
        Err(message) => format!("err {}", hex(&message)),
    };
    format!("{strict},{}", hex(&pyfmt::decode_utf8_replace(&bytes)))
}

fn codepoint_answer(case: &str) -> String {
    let c =
        char::from_u32(u32::from_str_radix(case, 16).expect("a code point")).expect("a character");
    let one = c.to_string();
    assert_eq!(pyfmt::char_is_space(c), pyfmt::is_space(&one));
    format!(
        "{},{},{},{},{},{}",
        flags(&one),
        hex(&pyfmt::lower(&one)),
        hex(&pyfmt::upper(&one)),
        hex(&pyfmt::repr_str(&one)),
        hex(&pyfmt::lower(&format!("a{c}Σ"))),
        hex(&pyfmt::lower(&format!("{c}Σ{c}"))),
    )
}

fn answer(file: &str, case: &str) -> String {
    let fields: Vec<&str> = case.split('~').collect();
    match file {
        "repr_str" => hex(&pyfmt::repr_str(&text(case))),
        "repr_f64" => hex(&pyfmt::repr_f64(float(case))),
        "format_fixed" => hex(&pyfmt::format_fixed(
            float(fields[0]),
            fields[1].parse().expect("a count of digits"),
        )),
        "int" => int_answer(case),
        "float" => float_answer(case),
        "int_of_f64" => int_of_f64_answer(case),
        "f64_of_int" => pyfmt::f64_of_int(&text(case)).map_or("overflow".to_string(), bits),
        "predicates" => flags(&text(case)),
        "strip" => {
            let text = text(case);
            hexes(&[
                pyfmt::strip(&text),
                pyfmt::lstrip(&text),
                pyfmt::rstrip(&text),
            ])
        }
        "strip_chars" => {
            let (text, chars) = (text(fields[0]), text(fields[1]));
            hexes(&[
                pyfmt::strip_chars(&text, &chars),
                pyfmt::lstrip_chars(&text, &chars),
                pyfmt::rstrip_chars(&text, &chars),
            ])
        }
        "split_whitespace" => items(&pyfmt::split_whitespace(&text(case))),
        "split_lines" => items(&pyfmt::split_lines(&text(case))),
        "split" => {
            let (text, sep) = (text(fields[0]), text(fields[1]));
            let most: usize = fields[2].parse().expect("a count of splits");
            [
                items(&pyfmt::split(&text, &sep)),
                items(&pyfmt::split_n(&text, &sep, most)),
                items(&pyfmt::rsplit_n(&text, &sep, most)),
            ]
            .join("|")
        }
        "partition" => {
            let (text, sep) = (text(fields[0]), text(fields[1]));
            let (head, found, tail) = pyfmt::partition(&text, &sep);
            let (rhead, rfound, rtail) = pyfmt::rpartition(&text, &sep);
            hexes(&[head, found, tail, rhead, rfound, rtail])
        }
        "case" => {
            let text = text(case);
            hexes(&[&pyfmt::lower(&text), &pyfmt::upper(&text)])
        }
        "affix" => {
            let (text, affix) = (text(fields[0]), text(fields[1]));
            hexes(&[
                pyfmt::removeprefix(&text, &affix),
                pyfmt::removesuffix(&text, &affix),
            ])
        }
        "slice" => {
            let text = text(fields[0]);
            let chars: usize = fields[3].parse().expect("a count of code points");
            format!(
                "{},{},{}",
                pyfmt::len(&text),
                hex(pyfmt::prefix(&text, chars)),
                hex(pyfmt::slice(&text, bound(fields[1]), bound(fields[2])))
            )
        }
        "utf8" => utf8_answer(case),
        "codepoints" => codepoint_answer(case),
        other => panic!("no answer for {other}"),
    }
}

const FILES: [&str; 19] = [
    "repr_str",
    "repr_f64",
    "format_fixed",
    "int",
    "float",
    "int_of_f64",
    "f64_of_int",
    "predicates",
    "strip",
    "strip_chars",
    "split_whitespace",
    "split_lines",
    "split",
    "partition",
    "case",
    "affix",
    "slice",
    "utf8",
    "codepoints",
];

/// The count of cases of every file under `dir`, each of at least `least` cases and
/// every case answered as the file expects.
fn agree(dir: &Path, least: usize) -> Vec<(&'static str, usize)> {
    let mut counts = Vec::new();
    let mut differences = Vec::new();
    for file in FILES {
        let path = dir.join(format!("{file}.tsv"));
        let lines = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        let mut count = 0;
        for line in lines.lines() {
            let (case, expected) = line.split_once('\t').expect("a case and its answer");
            let got = answer(file, case);
            if got != expected && differences.len() < 20 {
                let shown = |text: &str| text.chars().take(300).collect::<String>();
                differences.push(format!(
                    "{file} {}: expected {}, got {}",
                    shown(case),
                    shown(expected),
                    shown(&got)
                ));
            }
            count += 1;
        }
        assert!(count >= least, "{file}: {count} cases");
        counts.push((file, count));
    }
    assert!(differences.is_empty(), "{}", differences.join("\n"));
    counts
}

#[test]
fn a_cell_that_is_not_pairs_of_lowercase_hex_digits_is_refused() {
    assert_eq!(unhex(""), [0u8; 0]);
    assert_eq!(unhex("00_7f_ff"), [0x00, 0x7f, 0xff]);
    for damaged in [
        "6", "616", "61_", "_61", "61__62", "61 62", "6g", "6A", "+a",
    ] {
        assert!(
            std::panic::catch_unwind(|| unhex(damaged)).is_err(),
            "{damaged}"
        );
    }
}

#[test]
fn the_committed_cases_are_answered_as_the_interpreter_answers_them() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/pyfmt");
    agree(&dir, 100);
}

#[test]
#[ignore = "reads the case files under the directory VEP_FORMATS_COMPARE_PYDIFF_DIR names"]
fn every_case_of_the_whole_files_is_answered_as_the_interpreter_answers_it() {
    let dir =
        std::env::var_os("VEP_FORMATS_COMPARE_PYDIFF_DIR").expect("VEP_FORMATS_COMPARE_PYDIFF_DIR");
    for (file, count) in agree(&PathBuf::from(dir).join("pyfmt"), 20_000) {
        println!("pyfmt {file}: {count} cases agree");
    }
}
