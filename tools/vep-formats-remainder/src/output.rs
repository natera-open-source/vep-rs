// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the comparator reads back: each side's remainder (every header line, then every
//! handed-over data line, in file order, each ended by `\n`) and `counters.json`, the
//! counts of what was kept back.

use crate::counters::populated;
use crate::groups::Split;
use crate::input::{Side, HEADER};
use crate::layout::Layout;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

struct Written {
    data_lines: u64,
    handed_over: u64,
    bytes_handed_over: u64,
}

fn write_remainder(path: &Path, side: &Side, kept: &[bool]) -> Result<Written, String> {
    let fail = |e: std::io::Error| format!("cannot write {}: {e}", path.display());
    let mut out = BufWriter::with_capacity(1 << 20, fs::File::create(path).map_err(fail)?);
    let mut written = Written {
        data_lines: 0,
        handed_over: 0,
        bytes_handed_over: 0,
    };
    for (line, &kept) in side.lines.iter().zip(kept) {
        let data = line.flags & HEADER == 0;
        written.data_lines += u64::from(data);
        if kept {
            continue;
        }
        out.write_all(side.text(line)).map_err(fail)?;
        out.write_all(b"\n").map_err(fail)?;
        if data {
            written.handed_over += 1;
            written.bytes_handed_over += u64::from(line.len) + 1;
        }
    }
    out.flush().map_err(fail)?;
    Ok(written)
}

/// Runs of handed-over reference data lines, each with the number of kept-back lines that
/// precede it in the reference: what moves a place among the remainder's lines to its
/// place among all lines.
fn kept_before(side: &Side, kept: &[bool]) -> Vec<(u64, u64)> {
    let mut runs: Vec<(u64, u64)> = Vec::new();
    let mut before = 0u64;
    for (line, &kept) in side.lines.iter().zip(kept) {
        if line.flags & HEADER != 0 {
            continue;
        }
        if kept {
            before += 1;
            continue;
        }
        match runs.last_mut() {
            Some(run) if run.1 == before => run.0 += 1,
            _ => runs.push((1, before)),
        }
    }
    runs
}

/// A JSON string of UTF-8 text: the two characters JSON escapes and the control
/// characters below U+0020 escaped, everything else as it is.
fn json_string(text: &[u8]) -> String {
    let text = String::from_utf8_lossy(text);
    let mut out = String::with_capacity(text.len() + 2);
    out.push('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn counters_json(
    split: &Split,
    kept_entries: u64,
    counts: &BTreeMap<Vec<u8>, u64>,
    runs: &[(u64, u64)],
    expected: &Written,
    actual: &Written,
) -> String {
    let whole = split
        .whole
        .map_or("null".to_string(), |w| json_string(w.text().as_bytes()));
    let populated: Vec<String> = counts
        .iter()
        .map(|(name, n)| format!("  {}: {n}", json_string(name)))
        .collect();
    let runs: Vec<String> = runs.iter().map(|(n, b)| format!("[{n}, {b}]")).collect();
    let side = |w: &Written| {
        format!(
            "{{\"lines\": {}, \"handed_over\": {}, \"bytes_handed_over\": {}}}",
            w.data_lines, w.handed_over, w.bytes_handed_over
        )
    };
    let block = |items: &[String]| {
        if items.is_empty() {
            String::new()
        } else {
            format!("\n{}\n ", items.join(",\n"))
        }
    };
    format!(
        "{{\n \"whole\": {whole},\n \"kept_entries\": {kept_entries},\n \"kept_populated\": {{{}}},\n \"kept_before\": [{}],\n \"groups\": {},\n \"groups_handed_over\": {},\n \"expected\": {},\n \"actual\": {}\n}}\n",
        block(&populated),
        runs.join(", "),
        split.groups,
        split.groups_handed_over,
        side(expected),
        side(actual),
    )
}

pub fn write(
    out_dir: &Path,
    expected: &Side,
    actual: &Side,
    split: &Split,
    threads: usize,
) -> Result<(), String> {
    fs::create_dir_all(out_dir).map_err(|e| format!("cannot write {}: {e}", out_dir.display()))?;
    let expected_path = out_dir.join("expected.txt");
    let actual_path = out_dir.join("actual.txt");
    let (we, wa) = std::thread::scope(|scope| {
        let we = scope.spawn(|| write_remainder(&expected_path, expected, &split.kept_expected));
        let wa = scope.spawn(|| write_remainder(&actual_path, actual, &split.kept_actual));
        (
            we.join().expect("reference writer thread"),
            wa.join().expect("candidate writer thread"),
        )
    });
    let (we, wa) = (we?, wa?);
    let kept_entries = split.kept_expected.iter().filter(|k| **k).count() as u64;
    let counts = match Layout::of(expected) {
        Ok(layout) if kept_entries > 0 => {
            populated(expected, &split.kept_expected, &layout, threads)
        }
        _ => BTreeMap::new(),
    };
    let runs = kept_before(expected, &split.kept_expected);
    let json = counters_json(split, kept_entries, &counts, &runs, &we, &wa);
    let path = out_dir.join("counters.json");
    fs::write(&path, json).map_err(|e| format!("cannot write {}: {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_name_is_written_as_a_json_string() {
        assert_eq!(json_string(b"IMPACT"), "\"IMPACT\"");
        assert_eq!(json_string(b"a\"b\\c"), "\"a\\\"b\\\\c\"");
        assert_eq!(json_string(b"tab\there"), "\"tab\\u0009here\"");
        assert_eq!(json_string("5\u{2032}UTR".as_bytes()), "\"5\u{2032}UTR\"");
        assert_eq!(json_string(b""), "\"\"");
    }
}
