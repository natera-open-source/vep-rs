// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary over cells whose values carry escapes, lists, bytes outside UTF-8 and
//! empty parts.

mod common;

use common::{at, binary, custom_map, exit_code, one_cell, row, run, table, text, Scratch, QUIET};
use std::fs;
use std::process::Command;

#[test]
fn every_escape_decodes_and_only_a_bare_ampersand_splits_a_list() {
    let scratch = Scratch::new();
    let cell = "G|Likely%20benign%3B%3D%25%22|criteria_provided%2C_single_submitter\
                |One%7Ctwo&three%26four|a&b|SO:0001483";
    let run = one_cell(&scratch, "FV_CLINVAR", cell, &[]);
    assert_eq!(
        run.tuples(),
        at(&[
            ("ClinVar_CLNSIG", "Likely benign;=%\""),
            ("ClinVar_CLNREVSTAT", "criteria_provided,_single_submitter"),
            ("ClinVar_CLNDN", "One|two|three&four"),
            ("ClinVar_CLNVC", "a&b"),
            ("ClinVar_CLNVCSO", "SO:0001483"),
        ])
    );
}

#[test]
fn a_tab_or_line_breaking_character_stays_encoded_so_a_tuple_stays_one_line() {
    let scratch = Scratch::new();
    // The escape of each character, as `urllib.parse.quote` writes it, then `x`.
    let breaks = "\t\n\r\u{0b}\u{0c}\u{1c}\u{1d}\u{1e}\u{85}\u{2028}\u{2029}";
    let escapes: String = breaks
        .chars()
        .map(|c| {
            let mut escape: String = c
                .to_string()
                .bytes()
                .map(|byte| format!("%{byte:02X}"))
                .collect();
            escape.push('x');
            escape
        })
        .collect();
    assert!(escapes.starts_with("%09x%0Ax%0Dx%0Bx"));
    assert_eq!(escapes.matches('%').count(), 16);
    let run = one_cell(
        &scratch,
        "FV_CLINVAR",
        &format!("G|a{escapes}|raw\u{0b}break|||"),
        &[],
    );
    let written = String::from_utf8(run.tuple_bytes()).unwrap();
    assert_eq!(
        written,
        format!(
            "21:1000\tG\tENST1\tClinVar_CLNSIG\ta{escapes}\n\
             21:1000\tG\tENST1\tClinVar_CLNREVSTAT\traw%0Bbreak\n"
        )
    );
    // No character of the file but its two line feeds ends a line.
    let breaks = written
        .chars()
        .filter(|c| {
            matches!(
                c,
                '\n' | '\r' | '\u{0b}' | '\u{0c}' | '\u{1c}'
                    ..='\u{1e}' | '\u{85}' | '\u{2028}' | '\u{2029}'
            )
        })
        .count();
    assert_eq!(breaks, 2);
}

#[test]
fn an_escaped_space_decodes_to_the_space_the_default_format_writes() {
    // Ensembl VEP's default format keeps the spaces of a `--custom` value and a literal
    // `%20` as they are; its VCF format turns a run of spaces into one `_`.
    let scratch = Scratch::new();
    let run = one_cell(
        &scratch,
        "FV_CLINVAR",
        "G||||two%20words|a%20%20b%2520c",
        &[],
    );
    assert_eq!(
        run.tuples(),
        at(&[
            ("ClinVar_CLNVC", "two words"),
            ("ClinVar_CLNVCSO", "a  b%20c")
        ])
    );
}

#[test]
fn an_escaped_percent_sign_decodes_once() {
    let scratch = Scratch::new();
    let run = one_cell(&scratch, "FV_CLINVAR", "G||||50%2520off|", &[]);
    assert_eq!(run.tuples(), at(&[("ClinVar_CLNVC", "50%20off")]));
}

#[test]
fn values_are_utf8_whatever_the_locale_and_a_byte_that_is_not_passes_through() {
    let scratch = Scratch::new();
    let line = row("21:1000", "G", "ENST1", &["G|caf%C3%A9|BYTE|esc%FFaped||"]) + "\n";
    let mut bytes = table(&["FV_CLINVAR"], &[]).into_bytes();
    let (before, after) = line.split_once("BYTE").unwrap();
    bytes.extend_from_slice(before.as_bytes());
    bytes.extend_from_slice(b"byte\xff");
    bytes.extend_from_slice(after.as_bytes());
    let raw = scratch.write("raw.tab", &bytes);
    let tuple_file = scratch.path("tuples.tsv");
    let output = Command::new(common::BIN)
        .args([&raw, &text(&tuple_file)])
        .env("LC_ALL", "C")
        .output()
        .unwrap();
    assert_eq!(
        exit_code(&output),
        0,
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read(&tuple_file).unwrap(),
        b"21:1000\tG\tENST1\tClinVar_CLNSIG\tcaf\xc3\xa9\n\
          21:1000\tG\tENST1\tClinVar_CLNREVSTAT\tbyte\xff\n\
          21:1000\tG\tENST1\tClinVar_CLNDN\tesc\xffaped\n"
    );
}

#[test]
fn a_raw_byte_and_an_escape_read_as_one_utf8_text() {
    let scratch = Scratch::new();
    let line = row(
        "21:1000",
        "G",
        "ENST1",
        &["G|oneRAW%85two||||x%C3%A9|,G|||||xRAWER%A9|"],
    ) + "\n";
    let mut bytes = table(&["FV_CLINVAR"], &[]).into_bytes();
    // `RAWER` stands for the byte 0xC3 and `RAW` for the byte 0xC2.
    let mut rest = line.as_str();
    while let Some(at) = rest.find("RAW") {
        bytes.extend_from_slice(&rest.as_bytes()[..at]);
        if rest[at..].starts_with("RAWER") {
            bytes.push(0xc3);
            rest = &rest[at + 5..];
        } else {
            bytes.push(0xc2);
            rest = &rest[at + 3..];
        }
    }
    bytes.extend_from_slice(rest.as_bytes());
    let raw = scratch.write("raw.tab", &bytes);
    let tuple_file = scratch.path("tuples.tsv");
    let output = binary(&[&raw, &text(&tuple_file)]);
    assert_eq!(
        exit_code(&output),
        0,
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read(&tuple_file).unwrap(),
        b"21:1000\tG\tENST1\tClinVar_CLNSIG\tone%C2%85two\n\
          21:1000\tG\tENST1\tClinVar_CLNVCSO\tx\xc3\xa9\n"
    );
}

#[test]
fn a_list_takes_the_separator_its_map_row_names() {
    let scratch = Scratch::new();
    let custom = custom_map(
        &scratch,
        "FV_X\t1\tX_LIST\tX\ttext\tfabricated\t,\t-\n\
         FV_X\t2\tX_PLAIN\tX\ttext\tfabricated\t-\t-\n",
    );
    let run = one_cell(&scratch, "FV_X", "G|a&b&c|d&e", &["--field-map", &custom]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(run.tuples(), at(&[("X_LIST", "a,b,c"), ("X_PLAIN", "d&e")]));
}

#[test]
fn a_joined_range_and_a_score_half_are_percent_decoded() {
    let scratch = Scratch::new();
    let rows = [row(
        "21:1000",
        "G",
        "ENST1",
        &[
            "G|GENE%2DA|0.01|0.00|0.85|0.00|-12|3|0|27",
            "G|probably%5Fdamaging(0%2E998)|",
        ],
    )];
    let run = run(&scratch, &["SpliceAI", "FV_DBNSFP"], &rows, &[]);
    assert_eq!(
        run.tuples(),
        at(&[
            ("SpliceAI_pred", "GENE-A|0.01|0.00|0.85|0.00|-12|3|0|27"),
            ("SIFT_score", "0.998"),
        ])
    );
}

#[test]
fn an_empty_column_writes_an_empty_tuple_file() {
    let scratch = Scratch::new();
    let rows = [
        row("21:1000", "G", "ENST1", &["-"]),
        row("21:1000", "G", "ENST2", &[""]),
    ];
    let run = run(&scratch, &["FV_REVEL"], &rows, &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(run.tuple_bytes(), b"");
    assert_eq!(
        run.stderr,
        format!(
            "fastvep_plugin_tuples: wrote 0 tuples from 2 rows (0 too short); \
             by field: REVEL=0; mapped columns: FV_REVEL; \
             dropped parts with no Ensembl VEP field: 0; {QUIET}\n"
        )
    );
}

#[test]
fn a_header_without_rows_writes_an_empty_tuple_file() {
    let scratch = Scratch::new();
    let run = run(&scratch, &["FV_REVEL"], &[], &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(run.tuple_bytes(), b"");
    assert!(run.stderr.starts_with(
        "fastvep_plugin_tuples: wrote 0 tuples from 0 rows (0 too short); by field: REVEL=0; "
    ));
}

fn an_empty_or_missing_part_writes_no_tuple(column: &str, cell: &str, pairs: &[(&str, &str)]) {
    let scratch = Scratch::new();
    let run = one_cell(&scratch, column, cell, &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(run.tuples(), at(pairs));
    assert!(run
        .stderr
        .contains("dropped parts with no Ensembl VEP field: 0;"));
}

#[test]
fn an_empty_or_missing_part_writes_no_tuple_loeuf_empty() {
    an_empty_or_missing_part_writes_no_tuple(
        "FV_GNOMAD_GENE",
        "GENE1|0.9981|||",
        &[("pLI_gene_value", "0.9981")],
    );
}

#[test]
fn an_empty_or_missing_part_writes_no_tuple_pli_empty() {
    an_empty_or_missing_part_writes_no_tuple(
        "FV_GNOMAD_GENE",
        "GENE1||0.2345||",
        &[("LOEUF", "0.2345")],
    );
}

#[test]
fn an_empty_or_missing_part_writes_no_tuple_first_empty() {
    an_empty_or_missing_part_writes_no_tuple(
        "FV_ALPHAMISSENSE",
        "G||likely_benign",
        &[("am_class", "likely_benign")],
    );
}

#[test]
fn an_empty_or_missing_part_writes_no_tuple_short_entry() {
    an_empty_or_missing_part_writes_no_tuple(
        "FV_ALPHAMISSENSE",
        "G|0.5",
        &[("am_pathogenicity", "0.5")],
    );
}

#[test]
fn an_empty_or_missing_part_writes_no_tuple_allele_alone() {
    an_empty_or_missing_part_writes_no_tuple("FV_ALPHAMISSENSE", "G", &[]);
}

#[test]
fn an_empty_or_missing_part_writes_no_tuple_all_empty() {
    an_empty_or_missing_part_writes_no_tuple("SpliceAI", "G|||||||||", &[]);
}
