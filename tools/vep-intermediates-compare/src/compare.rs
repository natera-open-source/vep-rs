// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The merge of the two sorted sides: the last row under each key, the comparison of the
//! keys both sides hold, and the report's lines.

use crate::rowcodec::{key_len, take_field, unorderable, SEQUENCE_DIGITS};
use std::borrow::Cow;
use std::cmp::{Ordering, Reverse};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

const READ_BUF: usize = 4 << 20;
/// The patterns the summary lists, as `Counter.most_common(20)` cuts them.
const TOP_PATTERNS: usize = 20;

/// What the merge counted, in the report's terms.
#[derive(Default)]
pub struct Tally {
    pub trace_keys: u64,
    pub output_keys: u64,
    pub keys: u64,
    compared: u64,
    discordant: u64,
    missing_in_output: u64,
    missing_in_trace: u64,
    /// Diverging fields with their counts, in the order each was first seen.
    fields: Vec<(&'static str, u64)>,
    /// Consequence patterns with their count and the order each was first seen.
    patterns: HashMap<Vec<u8>, (u64, usize)>,
    /// Two keys with no order between them: the report ends after its three counts.
    pub unorderable: bool,
    /// The first compared row whose trace lacks a value the comparison reads: the report
    /// ends before that row.
    pub row_failure: Option<String>,
}

/// One side's sorted rows, a key at a time.
struct Groups {
    reader: BufReader<File>,
    current: Vec<u8>,
    ahead: Vec<u8>,
    has_ahead: bool,
    key_len: usize,
}

impl Groups {
    fn open(path: &Path) -> io::Result<Groups> {
        Ok(Groups {
            reader: BufReader::with_capacity(READ_BUF, File::open(path)?),
            current: Vec::new(),
            ahead: Vec::new(),
            has_ahead: false,
            key_len: 0,
        })
    }

    fn read_ahead(&mut self) -> io::Result<bool> {
        self.ahead.clear();
        if self.reader.read_until(b'\n', &mut self.ahead)? == 0 {
            return Ok(false);
        }
        if self.ahead.last() == Some(&b'\n') {
            self.ahead.pop();
        }
        Ok(true)
    }

    /// Moves to the next key and keeps its last row, the one a `dict` filled in file
    /// order would hold. Returns `false` after the last key.
    fn advance(&mut self) -> io::Result<bool> {
        if !self.has_ahead && !self.read_ahead()? {
            return Ok(false);
        }
        std::mem::swap(&mut self.current, &mut self.ahead);
        self.key_len = key_len(&self.current);
        loop {
            if !self.read_ahead()? {
                self.has_ahead = false;
                return Ok(true);
            }
            if self.ahead.starts_with(&self.current[..self.key_len]) {
                std::mem::swap(&mut self.current, &mut self.ahead);
            } else {
                self.has_ahead = true;
                return Ok(true);
            }
        }
    }

    fn key(&self) -> &[u8] {
        &self.current[..self.key_len]
    }

    fn values(&self) -> &[u8] {
        &self.current[self.key_len + SEQUENCE_DIGITS..]
    }
}

/// Merges the two sorted sides, writes the row lines to `rows` and returns the counts.
pub fn merge(
    trace_sorted: &Path,
    output_sorted: &Path,
    discordant_only: bool,
    rows: &mut impl Write,
) -> io::Result<Tally> {
    let mut trace = Groups::open(trace_sorted)?;
    let mut output = Groups::open(output_sorted)?;
    let mut tally = Tally::default();
    let mut previous_key: Vec<u8> = Vec::new();
    let mut line = Vec::with_capacity(1024);
    let mut in_trace = trace.advance()?;
    let mut in_output = output.advance()?;
    while in_trace || in_output {
        let order = match (in_trace, in_output) {
            (true, true) => trace.key().cmp(output.key()),
            (true, false) => Ordering::Less,
            _ => Ordering::Greater,
        };
        let key = if order == Ordering::Greater {
            output.key()
        } else {
            trace.key()
        };
        if tally.keys > 0 && unorderable(&previous_key, key) {
            tally.unorderable = true;
        }
        previous_key.clear();
        previous_key.extend_from_slice(key);
        tally.keys += 1;
        match order {
            Ordering::Less => {
                tally.trace_keys += 1;
                tally.missing_in_output += 1;
                in_trace = trace.advance()?;
            }
            Ordering::Greater => {
                tally.output_keys += 1;
                tally.missing_in_trace += 1;
                in_output = output.advance()?;
            }
            Ordering::Equal => {
                tally.trace_keys += 1;
                tally.output_keys += 1;
                if tally.row_failure.is_none() {
                    line.clear();
                    match compare_row(&trace, &output, &mut tally, &mut line) {
                        Ok(discordant) if discordant || !discordant_only => {
                            rows.write_all(&line)?;
                        }
                        Ok(_) => {}
                        Err(failure) => tally.row_failure = Some(failure),
                    }
                }
                in_trace = trace.advance()?;
                in_output = output.advance()?;
            }
        }
    }
    Ok(tally)
}

/// `set(cons.replace(",", "&").split("&")) - {""}`, in code-point order.
fn terms(consequence: &[u8]) -> Vec<&[u8]> {
    let mut terms: Vec<&[u8]> = consequence
        .split(|b| matches!(b, b',' | b'&'))
        .filter(|term| !term.is_empty())
        .collect();
    terms.sort_unstable();
    terms.dedup();
    terms
}

/// Whether `a.upper() == b.upper()` under `str.upper()`, the full Unicode mapping.
fn same_upper(a: &[u8], b: &[u8]) -> bool {
    if a.is_ascii() && b.is_ascii() {
        return a.eq_ignore_ascii_case(b);
    }
    let upper = |bytes: &[u8]| {
        std::str::from_utf8(bytes)
            .expect("both inputs were decoded as UTF-8")
            .to_uppercase()
    };
    upper(a) == upper(b)
}

/// A trace value as an f-string prints it: `None` for a column the row did not reach.
fn printed<'a>(value: &'a Option<Cow<'_, [u8]>>) -> &'a [u8] {
    value.as_deref().unwrap_or(&b"None"[..])
}

/// A trace value that is present and not empty.
fn truthy<'a>(value: &'a Option<Cow<'_, [u8]>>) -> Option<&'a [u8]> {
    value.as_deref().filter(|bytes| !bytes.is_empty())
}

const FIELD_ORDER: [&str; 6] = [
    "consequence",
    "cds_start",
    "ref_codon",
    "alt_codon",
    "ref_peptide",
    "alt_peptide",
];

/// Compares the rows the two sides hold under one key, counts the outcome, writes the
/// row's report line to `line` and returns whether a field diverges. `Err` names a trace
/// value the comparison calls a `str` method on and the row lacks.
fn compare_row(
    trace: &Groups,
    output: &Groups,
    tally: &mut Tally,
    line: &mut Vec<u8>,
) -> Result<bool, String> {
    let mut key = trace.key();
    let mut key_field = || take_field(&mut key).expect("a key both sides hold has every field");
    let (location, allele, transcript) = (key_field(), key_field(), key_field());

    let mut values = trace.values();
    let mut trace_value = || take_field(&mut values);
    let (cds_start, cds_end) = (trace_value(), trace_value());
    let (ref_codon, alt_codon) = (trace_value(), trace_value());
    let (ref_peptide, alt_peptide) = (trace_value(), trace_value());
    let trace_consequence = trace_value();

    let mut values = output.values();
    let mut output_value = || take_field(&mut values).expect("an output row has every value");
    let (consequence, cds_position) = (output_value(), output_value());
    let (amino_acids, codons) = (output_value(), output_value());

    let lacks = |column: &str| {
        format!(
            "the trace row of {} {} {} is shorter than its header and has no {column}",
            String::from_utf8_lossy(&location),
            String::from_utf8_lossy(&allele),
            String::from_utf8_lossy(&transcript),
        )
    };

    tally.compared += 1;
    let Some(trace_consequence) = trace_consequence.as_deref() else {
        return Err(lacks("consequence_set"));
    };
    let trace_terms = terms(trace_consequence);
    let output_terms = terms(&consequence);
    let mut divergent = [false; FIELD_ORDER.len()];
    divergent[0] = trace_terms != output_terms;

    let output_cds: &[u8] = if &*cds_position == b"-" {
        b""
    } else {
        cds_position
            .split(|b| *b == b'-')
            .next()
            .unwrap_or_default()
    };
    if let Some(cds_start) = truthy(&cds_start) {
        divergent[1] = !output_cds.is_empty() && cds_start != output_cds;
    }

    if let Some(ref_codon) = truthy(&ref_codon) {
        let parts: Vec<&[u8]> = codons.split(|b| *b == b'/').collect();
        if &*codons != b"-" && parts.len() == 2 {
            divergent[2] = !same_upper(ref_codon, parts[0]);
            let Some(alt_codon) = alt_codon.as_deref() else {
                return Err(lacks("alt_codon"));
            };
            divergent[3] = !same_upper(alt_codon, parts[1]);
        }
    }

    if let Some(ref_peptide) = truthy(&ref_peptide) {
        let parts: Vec<&[u8]> = amino_acids.split(|b| *b == b'/').collect();
        if &*amino_acids != b"-" && parts.len() == 2 {
            divergent[4] = ref_peptide != parts[0];
            divergent[5] = alt_peptide.as_deref() != Some(parts[1]);
        }
    }

    let discordant = divergent.contains(&true);
    if discordant {
        tally.discordant += 1;
        for name in FIELD_ORDER
            .iter()
            .zip(divergent)
            .filter_map(|(n, d)| d.then_some(*n))
        {
            match tally.fields.iter_mut().find(|(seen, _)| *seen == name) {
                Some((_, count)) => *count += 1,
                None => tally.fields.push((name, 1)),
            }
        }
        if divergent[0] {
            let mut pattern = b"perl=".to_vec();
            pattern.extend_from_slice(&trace_terms.join(&b','));
            pattern.extend_from_slice(b" rust=");
            pattern.extend_from_slice(&output_terms.join(&b','));
            let first_seen = tally.patterns.len();
            tally.patterns.entry(pattern).or_insert((0, first_seen)).0 += 1;
        }
    }

    line.extend_from_slice(if discordant { b"DISCORD" } else { b"match" });
    for part in [&*location, &*allele, &*transcript] {
        line.push(b'\t');
        line.extend_from_slice(part);
    }
    line.push(b'\t');
    let names: Vec<&[u8]> = FIELD_ORDER
        .iter()
        .zip(divergent)
        .filter_map(|(n, d)| d.then_some(n.as_bytes()))
        .collect();
    line.extend_from_slice(&names.join(&b','));
    let parts: [&[u8]; 20] = [
        b"\tperl_cons=",
        trace_consequence,
        b"\trust_cons=",
        &consequence,
        b"\tperl_cds=",
        printed(&cds_start),
        b"-",
        printed(&cds_end),
        b"\tperl_pep=",
        printed(&ref_peptide),
        b"/",
        printed(&alt_peptide),
        b"\trust_aa=",
        &amino_acids,
        b"\tperl_codon=",
        printed(&ref_codon),
        b"/",
        printed(&alt_codon),
        b"\trust_codon=",
        &codons,
    ];
    for part in parts {
        line.extend_from_slice(part);
    }
    line.push(b'\n');
    Ok(discordant)
}

/// The three counts that open the report and the blank line after them.
pub fn write_counts(out: &mut impl Write, tally: &Tally) -> io::Result<()> {
    writeln!(out, "Perl rows: {}", tally.trace_keys)?;
    writeln!(out, "Rust rows: {}", tally.output_keys)?;
    writeln!(out, "Unique keys: {}", tally.keys)?;
    writeln!(out)
}

/// The summary that closes the report. Both lists keep `Counter.most_common`'s order:
/// the larger count first and, between equal counts, the one counted first.
pub fn write_summary(out: &mut impl Write, tally: &Tally) -> io::Result<()> {
    let rule = "=".repeat(60);
    writeln!(out, "\n{rule}\nSUMMARY\n{rule}")?;
    writeln!(out, "Total compared:    {}", tally.compared)?;
    writeln!(out, "Total discordant:  {}", tally.discordant)?;
    writeln!(out, "Missing in Rust:   {}", tally.missing_in_output)?;
    writeln!(out, "Missing in Perl:   {}", tally.missing_in_trace)?;
    writeln!(out, "\nField divergence counts:")?;
    let mut fields = tally.fields.clone();
    fields.sort_by_key(|field| Reverse(field.1));
    for (field, count) in fields {
        let percent = 100.0 * count as f64 / tally.discordant as f64;
        writeln!(out, "  {field:<20}  {count:>6}  ({percent:.1}%)")?;
    }
    writeln!(out, "\nTop consequence divergence patterns:")?;
    let mut patterns: Vec<(&Vec<u8>, u64, usize)> = tally
        .patterns
        .iter()
        .map(|(pattern, &(count, first_seen))| (pattern, count, first_seen))
        .collect();
    patterns.sort_by_key(|&(_, count, first_seen)| (Reverse(count), first_seen));
    for (pattern, count, _) in patterns.into_iter().take(TOP_PATTERNS) {
        write!(out, "  {count:>6}  ")?;
        out.write_all(pattern)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rowcodec::{push_field, push_sequence};

    #[test]
    fn consequence_terms_are_a_set_over_both_separators() {
        let set = |s: &'static str| terms(s.as_bytes());
        assert_eq!(set("b,a&b,,a"), [&b"a"[..], b"b"]);
        assert_eq!(set(""), Vec::<&[u8]>::new());
        assert_eq!(set("&,"), Vec::<&[u8]>::new());
        assert_eq!(set(" a ,a"), [&b" a "[..], b"a"]);
    }

    #[test]
    fn codons_compare_under_the_full_uppercase_mapping() {
        let same = |a: &str, b: &str| same_upper(a.as_bytes(), b.as_bytes());
        assert!(same("aGt", "AGT"));
        assert!(!same("aGt", "AGA"));
        assert!(same("stra\u{df}e", "STRASSE"));
        assert!(same("\u{149}", "\u{2bc}N"));
        assert!(same("\u{3c3}", "\u{3c2}"));
        assert!(!same("\u{e9}", "e"));
        assert!(!same("i", "\u{130}"));
    }

    #[test]
    fn percentages_round_half_to_even_on_the_binary_value() {
        let cases = [
            (49, 400, "12.2"),
            (1, 16, "6.2"),
            (3, 16, "18.8"),
            (1, 3, "33.3"),
            (2, 3, "66.7"),
            (1, 1, "100.0"),
            (1, 2000, "0.1"),
            (1, 2001, "0.0"),
        ];
        for (count, total, expected) in cases {
            let percent = 100.0 * count as f64 / total as f64;
            assert_eq!(format!("{percent:.1}"), expected, "{count} of {total}");
        }
    }

    /// Rows of one side: the three key fields, then the value fields.
    type Rows<'a> = [(&'a [Option<&'a str>; 3], &'a [Option<&'a str>])];

    /// Writes one side's sorted file.
    fn side(path: &Path, rows: &Rows) {
        let mut lines: Vec<Vec<u8>> = rows
            .iter()
            .enumerate()
            .map(|(sequence, (key, values))| {
                let mut line = Vec::new();
                for field in key.iter() {
                    push_field(&mut line, field.map(str::as_bytes));
                }
                push_sequence(&mut line, sequence as u64);
                for value in values.iter() {
                    push_field(&mut line, value.map(str::as_bytes));
                }
                line
            })
            .collect();
        lines.sort();
        std::fs::write(path, [lines.join(&b'\n'), b"\n".to_vec()].concat()).unwrap();
    }

    fn merged(tag: &str, trace: &Rows, output: &Rows, discordant_only: bool) -> (Tally, String) {
        let dir = std::env::temp_dir().join(format!(
            "vep-intermediates-compare-merge-{}-{tag}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let (trace_path, output_path) = (dir.join("trace"), dir.join("output"));
        side(&trace_path, trace);
        side(&output_path, output);
        let mut rows = Vec::new();
        let tally = merge(&trace_path, &output_path, discordant_only, &mut rows).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        (tally, String::from_utf8(rows).unwrap())
    }

    const K1: [Option<&str>; 3] = [Some("1:100"), Some("A"), Some("ENST1")];
    const K2: [Option<&str>; 3] = [Some("1:200"), Some("T"), Some("ENST2")];
    const K3: [Option<&str>; 3] = [Some("1:300"), Some("-"), Some("ENST3")];

    #[test]
    fn the_last_row_under_a_key_is_compared_and_one_sided_keys_are_counted() {
        let trace_first = [
            Some("9"),
            Some("9"),
            Some("ccc"),
            Some("ctc"),
            Some("P"),
            Some("L"),
            Some("stop_gained"),
        ];
        let trace_last = [
            Some("5"),
            Some("5"),
            Some("gCa"),
            Some("gTa"),
            Some("A"),
            Some("V"),
            Some("missense_variant,splice_region_variant"),
        ];
        let intron = [
            Some(""),
            Some(""),
            Some(""),
            Some(""),
            Some(""),
            Some(""),
            Some("intron_variant"),
        ];
        let (tally, rows) = merged(
            "last-row",
            &[(&K1, &trace_first), (&K2, &intron), (&K1, &trace_last)],
            &[
                (
                    &K1,
                    &[
                        Some("splice_region_variant&missense_variant"),
                        Some("5-7"),
                        Some("A/V"),
                        Some("GCA/GTA"),
                    ],
                ),
                (
                    &K3,
                    &[Some("intron_variant"), Some("-"), Some("-"), Some("-")],
                ),
            ],
            false,
        );
        assert_eq!((tally.trace_keys, tally.output_keys, tally.keys), (2, 2, 3));
        assert_eq!(
            (
                tally.compared,
                tally.discordant,
                tally.missing_in_output,
                tally.missing_in_trace
            ),
            (1, 0, 1, 1)
        );
        assert_eq!(
            rows,
            "match\t1:100\tA\tENST1\t\tperl_cons=missense_variant,splice_region_variant\t\
             rust_cons=splice_region_variant&missense_variant\tperl_cds=5-5\tperl_pep=A/V\t\
             rust_aa=A/V\tperl_codon=gCa/gTa\trust_codon=GCA/GTA\n"
        );
        assert!(!tally.unorderable && tally.row_failure.is_none());
    }

    #[test]
    fn each_diverging_field_is_named_in_the_fixed_order() {
        let trace = [
            Some("5"),
            Some("6"),
            Some("gCa"),
            Some("gTa"),
            Some("A"),
            Some("V"),
            Some("missense_variant"),
        ];
        let (tally, rows) = merged(
            "fields",
            &[(&K1, &trace)],
            &[(
                &K1,
                &[
                    Some("stop_gained"),
                    Some("7-8"),
                    Some("a/*"),
                    Some("gCc/tGa"),
                ],
            )],
            true,
        );
        assert_eq!(
            rows,
            "DISCORD\t1:100\tA\tENST1\tconsequence,cds_start,ref_codon,alt_codon,ref_peptide,alt_peptide\t\
             perl_cons=missense_variant\trust_cons=stop_gained\tperl_cds=5-6\tperl_pep=A/V\t\
             rust_aa=a/*\tperl_codon=gCa/gTa\trust_codon=gCc/tGa\n"
        );
        let mut summary = Vec::new();
        write_summary(&mut summary, &tally).unwrap();
        let summary = String::from_utf8(summary).unwrap();
        assert!(summary.contains("  consequence                1  (100.0%)\n  cds_start "));
        assert!(summary.ends_with(
            "Top consequence divergence patterns:\n       1  perl=missense_variant rust=stop_gained\n"
        ));
    }

    #[test]
    fn values_the_output_leaves_out_are_not_compared() {
        let trace = [
            Some("5"),
            Some("6"),
            Some("gCa"),
            Some("gTa"),
            Some("A"),
            Some("V"),
            Some("missense_variant"),
        ];
        for output in [
            [Some("missense_variant"), Some("-"), Some("-"), Some("-")],
            [Some("missense_variant"), Some("-9"), Some("A"), Some("gCa")],
            [
                Some("missense_variant"),
                Some("5-9"),
                Some("A/V/L"),
                Some("x/y/z"),
            ],
        ] {
            let (tally, rows) = merged("left-out", &[(&K1, &trace)], &[(&K1, &output)], true);
            assert_eq!((tally.compared, tally.discordant), (1, 0), "{output:?}");
            assert!(rows.is_empty());
        }
    }

    #[test]
    fn a_trace_row_without_its_consequence_ends_the_rows_and_keeps_the_counts() {
        let good = [
            Some(""),
            Some(""),
            Some(""),
            Some(""),
            Some(""),
            Some(""),
            Some("intron_variant"),
        ];
        let short = [Some("5"), None, None, None, None, None, None];
        let output = [Some("intron_variant"), Some("-"), Some("-"), Some("-")];
        let (tally, rows) = merged(
            "short-row",
            &[(&K1, &good), (&K2, &short), (&K3, &good)],
            &[(&K1, &output), (&K2, &output), (&K3, &output)],
            false,
        );
        assert_eq!(rows.lines().count(), 1);
        assert!(rows.starts_with("match\t1:100\t"));
        assert_eq!((tally.trace_keys, tally.output_keys, tally.keys), (3, 3, 3));
        assert_eq!(
            tally.row_failure.as_deref(),
            Some("the trace row of 1:200 T ENST2 is shorter than its header and has no consequence_set")
        );
    }

    #[test]
    fn absent_trace_values_print_as_none_and_an_absent_alt_peptide_diverges() {
        let trace = [
            None,
            None,
            Some(""),
            None,
            Some("A"),
            None,
            Some("missense_variant"),
        ];
        let (tally, rows) = merged(
            "none",
            &[(&K1, &trace)],
            &[(
                &K1,
                &[
                    Some("missense_variant"),
                    Some("5"),
                    Some("A/V"),
                    Some("gCa/gTa"),
                ],
            )],
            false,
        );
        assert_eq!(
            rows,
            "DISCORD\t1:100\tA\tENST1\talt_peptide\tperl_cons=missense_variant\t\
             rust_cons=missense_variant\tperl_cds=None-None\tperl_pep=A/None\trust_aa=A/V\t\
             perl_codon=/None\trust_codon=gCa/gTa\n"
        );
        assert!(tally.row_failure.is_none());
        // A reference codon with no alternative one is the row that cannot be compared.
        let trace = [
            None,
            None,
            Some("gCa"),
            None,
            None,
            None,
            Some("missense_variant"),
        ];
        let (tally, rows) = merged(
            "no-alt-codon",
            &[(&K1, &trace)],
            &[(
                &K1,
                &[
                    Some("missense_variant"),
                    Some("5"),
                    Some("A/V"),
                    Some("gCa/gTa"),
                ],
            )],
            false,
        );
        assert!(rows.is_empty());
        assert!(tally.row_failure.unwrap().ends_with("has no alt_codon"));
    }

    #[test]
    fn keys_without_an_order_are_reported_whichever_side_holds_them() {
        let good = [
            Some(""),
            Some(""),
            Some(""),
            Some(""),
            Some(""),
            Some(""),
            Some("intron_variant"),
        ];
        let output = [Some("intron_variant"), Some("-"), Some("-"), Some("-")];
        let short_key = [Some("1:100"), None, None];
        let lone_key = [Some("1:150"), None, None];
        let (tally, _) = merged(
            "unorderable",
            &[(&short_key, &good)],
            &[(&K1, &output)],
            false,
        );
        assert!(tally.unorderable);
        assert_eq!((tally.trace_keys, tally.output_keys, tally.keys), (1, 1, 2));
        let (tally, _) = merged(
            "orderable",
            &[(&lone_key, &good), (&K2, &good)],
            &[(&K1, &output)],
            false,
        );
        assert!(!tally.unorderable);
        assert_eq!(tally.missing_in_output, 2);
    }

    #[test]
    fn the_summary_orders_equal_counts_by_first_appearance_and_lists_twenty_patterns() {
        let mut tally = Tally {
            compared: 40,
            discordant: 8,
            missing_in_output: 2,
            missing_in_trace: 3,
            fields: vec![
                ("ref_codon", 2),
                ("consequence", 5),
                ("cds_start", 2),
                ("alt_peptide", 5),
            ],
            ..Tally::default()
        };
        for i in 0..25usize {
            let count = if i == 7 { 3 } else { 1 };
            tally
                .patterns
                .insert(format!("perl=p{i:02} rust=r").into_bytes(), (count, i));
        }
        let mut out = Vec::new();
        write_summary(&mut out, &tally).unwrap();
        let text = String::from_utf8(out).unwrap();
        let expected_head = format!(
            "\n{rule}\nSUMMARY\n{rule}\nTotal compared:    40\nTotal discordant:  8\n\
             Missing in Rust:   2\nMissing in Perl:   3\n\nField divergence counts:\n  \
             consequence                5  (62.5%)\n  alt_peptide                5  (62.5%)\n  \
             ref_codon                  2  (25.0%)\n  cds_start                  2  (25.0%)\n\n\
             Top consequence divergence patterns:\n       3  perl=p07 rust=r\n       1  perl=p00 rust=r\n",
            rule = "=".repeat(60)
        );
        assert!(text.starts_with(&expected_head), "{text}");
        let listed: Vec<&str> = text
            .lines()
            .filter(|line| line.contains("perl=p"))
            .collect();
        assert_eq!(listed.len(), 20);
        assert_eq!(listed[19], "       1  perl=p19 rust=r");
    }
}
