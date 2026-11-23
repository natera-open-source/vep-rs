// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `pyjson` against stored answers of CPython.
//!
//! Each file under `tests/data/pyjson` holds `case<TAB>expected` lines written by CPython.
//! A text is its UTF-8 bytes, each as two lowercase hex digits, joined by `_`; a float is
//! the eight bytes of its bits, written the same way; the text of a float is a text; the
//! parts of a case are joined by `~`. The ignored tests read the files of the same names,
//! with many more lines, under the `pyjson` folder of the directory
//! `VEP_FORMATS_COMPARE_PYDIFF_DIR` names, and a file of JSON objects, one a line, that
//! `VEP_FORMATS_COMPARE_PYDIFF_JSON` names.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;
use vep_formats_compare::pyjson::{self, DumpOptions, JsonErrorKind, Value};

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

/// A value with its types spelt out: `n`, `t`, `f`, `i<decimal>`, `d<bits>`, `s<hex>`, a
/// list in brackets and an object in braces with `<hex of the name>:<value>` members.
fn spelt(value: &Value, out: &mut String) {
    match value {
        Value::Null => out.push('n'),
        Value::Bool(true) => out.push('t'),
        Value::Bool(false) => out.push('f'),
        Value::Int(whole) => {
            let _ = write!(out, "i{whole}");
        }
        Value::BigInt(whole) => {
            let _ = write!(out, "i{whole}");
        }
        Value::Float(float) if float.is_nan() => out.push_str("dnan"),
        Value::Float(float) => {
            let _ = write!(out, "d{}", bits(*float));
        }
        Value::Str(text) => {
            let _ = write!(out, "s{}", hex(text));
        }
        Value::List(items) => {
            out.push('[');
            for (at, item) in items.iter().enumerate() {
                if at > 0 {
                    out.push(',');
                }
                spelt(item, out);
            }
            out.push(']');
        }
        Value::Dict(members) => {
            out.push('{');
            for (at, (name, item)) in members.iter().enumerate() {
                if at > 0 {
                    out.push(',');
                }
                let _ = write!(out, "{}:", hex(name));
                spelt(item, out);
            }
            out.push('}');
        }
    }
}

/// What `loads` makes of a text: the value spelt out, `D<depth>` for a value nested more
/// than 100 deep, `E<hex>` and `V<hex>` for the two kinds of `ValueError` with the
/// message, `R` for a `RecursionError`, `L` for a string with half a surrogate pair.
fn loads_answer(case: &str) -> String {
    match pyjson::loads(&text(case)) {
        Ok(value) => {
            let depth = pyjson::depth(&value);
            if depth > 100 {
                return format!("D{depth}");
            }
            let mut out = String::new();
            spelt(&value, &mut out);
            out
        }
        Err(error) => {
            assert_eq!(
                error.is_value_error(),
                matches!(error.kind, JsonErrorKind::Decode | JsonErrorKind::IntDigits)
            );
            match error.kind {
                JsonErrorKind::Decode => format!("E{}", hex(&error.to_string())),
                JsonErrorKind::IntDigits => format!("V{}", hex(&error.to_string())),
                JsonErrorKind::TooDeep => "R".to_string(),
                JsonErrorKind::LoneSurrogate => "L".to_string(),
            }
        }
    }
}

fn read(hex: &str) -> Value {
    let text = text(hex);
    pyjson::loads(&text).unwrap_or_else(|error| panic!("{text}: {error}"))
}

fn answer(file: &str, case: &str) -> String {
    let fields: Vec<&str> = case.split('~').collect();
    match file {
        "loads" | "loads_deep" => loads_answer(case),
        "dumps" => {
            let options = DumpOptions {
                indent: (fields[1] != "-").then(|| fields[1].parse().expect("an indent")),
                sort_keys: fields[2] == "1",
                ensure_ascii: fields[3] == "1",
            };
            hex(&pyjson::dumps(&read(fields[0]), &options))
        }
        "compact" => {
            let value = read(case);
            format!(
                "{},{},{},{},{},{}",
                hex(&pyjson::dumps_compact(&value)),
                hex(&pyjson::json_scalar(&value)),
                hex(&pyjson::py_repr(&value)),
                hex(&pyjson::py_str(&value)),
                u8::from(pyjson::truthy(&value)),
                pyjson::depth(&value)
            )
        }
        "dumps_str" => {
            let text = text(case);
            format!(
                "{},{}",
                hex(&pyjson::dumps_str(&text, true)),
                hex(&pyjson::dumps_str(&text, false))
            )
        }
        "json_float" => hex(&pyjson::json_float(float(case))),
        "py_float" => match pyjson::py_float(&read(case)) {
            Ok(value) if value.is_nan() => "ok nan".to_string(),
            Ok(value) => format!("ok {}", bits(value)),
            Err(error) => {
                let class = match error {
                    pyjson::FloatError::Type(_) => "TypeError",
                    pyjson::FloatError::Value(_) => "ValueError",
                    pyjson::FloatError::Overflow(_) => "OverflowError",
                };
                format!("err {class} {}", hex(&error.to_string()))
            }
        },
        "eq" => {
            let (a, b) = (read(fields[0]), read(fields[1]));
            assert_eq!(pyjson::py_eq(&a, &b), pyjson::py_eq(&b, &a));
            u8::from(pyjson::py_eq(&a, &b)).to_string()
        }
        other => panic!("no answer for {other}"),
    }
}

/// The count of cases of each of `files` under `dir`, each of at least `least` cases and
/// every case answered as the file expects.
fn agree(dir: &Path, files: &[&'static str], least: usize) -> Vec<(&'static str, usize)> {
    let mut counts = Vec::new();
    let mut differences = Vec::new();
    for file in files {
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
        counts.push((*file, count));
    }
    assert!(differences.is_empty(), "{}", differences.join("\n"));
    counts
}

const FILES: [&str; 7] = [
    "loads",
    "dumps",
    "compact",
    "dumps_str",
    "json_float",
    "eq",
    "py_float",
];

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
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/pyjson");
    agree(&dir, &FILES, 100);
}

#[test]
#[ignore = "reads the case files under the directory VEP_FORMATS_COMPARE_PYDIFF_DIR names"]
fn every_case_of_the_whole_files_is_answered_as_the_interpreter_answers_it() {
    let dir =
        std::env::var_os("VEP_FORMATS_COMPARE_PYDIFF_DIR").expect("VEP_FORMATS_COMPARE_PYDIFF_DIR");
    let dir = PathBuf::from(dir).join("pyjson");
    let mut counts = agree(&dir, &FILES, 20_000);
    counts.extend(agree(&dir, &["loads_deep"], 1));
    for (file, count) in counts {
        println!("pyjson {file}: {count} cases agree");
    }
}

#[test]
#[ignore = "times the reader and the compact writer on the file VEP_FORMATS_COMPARE_PYDIFF_JSON names"]
fn every_object_of_a_file_is_read_and_written_compact_and_the_pace_is_printed() {
    let path = std::env::var_os("VEP_FORMATS_COMPARE_PYDIFF_JSON")
        .expect("VEP_FORMATS_COMPARE_PYDIFF_JSON");
    let content = std::fs::read_to_string(&path).expect("a file of JSON lines");
    let megabytes = content.len() as f64 / 1e6;
    // One line at a time, as the stream reads an output: a value lives for its line.
    let (mut objects, mut written) = (0usize, 0usize);
    let (mut read, mut wrote, mut dropped) = (0.0, 0.0, 0.0);
    for line in content.lines() {
        let started = Instant::now();
        let value = pyjson::loads(line).expect("a JSON object");
        read += started.elapsed().as_secs_f64();
        let started = Instant::now();
        written += pyjson::dumps_compact(&value).len();
        wrote += started.elapsed().as_secs_f64();
        let started = Instant::now();
        drop(value);
        dropped += started.elapsed().as_secs_f64();
        objects += 1;
    }
    println!(
        "{objects} objects, {megabytes:.0} MB: loads {:.0} MB/s, dumps_compact {:.0} MB/s of input ({:.0} MB/s of output), dropping the values {:.0} MB/s of input",
        megabytes / read,
        megabytes / wrote,
        written as f64 / 1e6 / wrote,
        megabytes / dropped,
    );
}
