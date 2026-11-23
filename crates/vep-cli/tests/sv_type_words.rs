// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The symbolic type words `<INS:ME:*>`, `<DEL:ME:*>` and `<DUP:TANDEM>`
//! through the `vep` binary. Each test annotates the records of the GRCh38 and
//! GRCh37 golden corpora whose ALT opens with its word and holds their rows to
//! the class word `get_SO_term` gives the ALT (`Parser.pm:683-698`) and to the
//! features and consequence terms of the corpus's Ensembl VEP output.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;

mod common;

/// The golden corpora whose records carry the three words.
const CORPORA: [&str; 2] = ["GRCh38", "GRCh37"];

/// The tab columns that identify a row and state its class and consequence.
const COLUMNS: [&str; 7] = [
    "#Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
];

const LOCATION: usize = 1;
const ALLELE: usize = 2;
const FEATURE: usize = 4;

fn corpus(name: &str) -> common::Corpus {
    common::corpora()
        .into_iter()
        .find(|corpus| corpus.release == "116" && corpus.name == name)
        .unwrap_or_else(|| panic!("no golden corpus 116/{name}"))
}

/// The [`COLUMNS`] cells of every data row of a tab output, by the row's
/// `Uploaded_variation`.
fn rows_by_record(text: &str) -> BTreeMap<String, BTreeSet<Vec<String>>> {
    let header: Vec<&str> = text
        .lines()
        .find(|line| line.starts_with(COLUMNS[0]))
        .expect("a tab header line")
        .split('\t')
        .collect();
    let indexes = COLUMNS.map(|name| {
        header
            .iter()
            .position(|cell| *cell == name)
            .unwrap_or_else(|| panic!("the tab header declares {name}"))
    });
    let mut rows: BTreeMap<String, BTreeSet<Vec<String>>> = BTreeMap::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let cells: Vec<&str> = line.split('\t').collect();
        let row: Vec<String> = indexes.iter().map(|&i| cells[i].to_string()).collect();
        rows.entry(row[0].clone()).or_default().insert(row);
    }
    rows
}

/// The Location, Allele and Feature of every row the corpus manifest documents
/// as differing from Ensembl VEP's, which the golden comparison holds to the
/// documented difference.
fn documented_rows(corpus: &common::Corpus) -> BTreeSet<(String, String, String)> {
    corpus.manifest["divergences"]
        .as_array()
        .map(|items| {
            items
                .iter()
                .map(|item| {
                    let cell = |key: &str| item[key].as_str().unwrap_or_default().to_string();
                    (cell("location"), cell("allele"), cell("feature"))
                })
                .collect()
        })
        .unwrap_or_default()
}

/// Annotates, corpus by corpus, the records whose ALT opens with `word` and
/// whose ID names one record, and asserts that every row of a record carries
/// the class word `classes` gives its ALT and that its rows outside the
/// documented ones are those of the corpus's Ensembl VEP output. Every entry of
/// `classes` must name an ALT of some record.
fn assert_word(word: &str, classes: &[(&str, &str)]) {
    let dir = tempfile::tempdir().unwrap();
    let mut alts_seen = BTreeSet::new();
    for name in CORPORA {
        let corpus = corpus(name);
        let inputs = common::input_records(&corpus);
        let mut id_counts: BTreeMap<&str, usize> = BTreeMap::new();
        for line in &inputs {
            *id_counts
                .entry(line.split('\t').nth(2).unwrap())
                .or_default() += 1;
        }
        let selected: Vec<(&str, &str, &str)> = inputs
            .iter()
            .filter_map(|line| {
                let cells: Vec<&str> = line.split('\t').collect();
                (cells[4].starts_with(word) && cells[2] != "." && id_counts[cells[2]] == 1)
                    .then_some((line.as_str(), cells[2], cells[4]))
            })
            .collect();
        assert!(!selected.is_empty(), "{name}: no record of {word}");

        let input = dir.path().join(format!("{name}.vcf"));
        let lines: Vec<&str> = selected.iter().map(|(line, _, _)| *line).collect();
        fs::write(
            &input,
            format!(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}\n",
                lines.join("\n")
            ),
        )
        .unwrap();
        let mut run = common::default_run(&corpus);
        run.input = input.to_str().unwrap().to_string();
        let actual = rows_by_record(&common::run_vep(&run, "tab", dir.path(), &[]));
        let expected = rows_by_record(&run.expected("tab.txt"));
        let documented = documented_rows(&corpus);
        let undocumented = |rows: &BTreeSet<Vec<String>>| -> BTreeSet<Vec<String>> {
            rows.iter()
                .filter(|row| {
                    !documented.contains(&(
                        row[LOCATION].clone(),
                        row[ALLELE].clone(),
                        row[FEATURE].clone(),
                    ))
                })
                .cloned()
                .collect()
        };

        for (_, id, alt) in selected {
            let class = classes
                .iter()
                .find(|(text, _)| *text == alt)
                .unwrap_or_else(|| panic!("{name} {id}: no class word listed for {alt}"))
                .1;
            alts_seen.insert(alt.to_string());
            let rows = actual
                .get(id)
                .unwrap_or_else(|| panic!("{name} {id} ({alt}): the run writes no row"));
            assert!(
                rows.iter().all(|row| row[ALLELE] == class),
                "{name} {id} ({alt}): every row's Allele is {class}: {rows:?}"
            );
            let reference = expected
                .get(id)
                .unwrap_or_else(|| panic!("{name} {id} ({alt}): the reference has no row"));
            assert_eq!(
                undocumented(rows),
                undocumented(reference),
                "{name} {id} ({alt}): the rows are Ensembl VEP's"
            );
        }
    }
    for (alt, _) in classes {
        assert!(alts_seen.contains(*alt), "no corpus record of {alt}");
    }
}

/// `<INS:ME:*>`: an element of `get_SO_term`'s list (`ALU`, `HERV`, `LINE1`,
/// `SVA`) names its own insertion term, and a bare `<INS:ME>` or an element
/// outside the list (`LINE`) is a `mobile_element_insertion`
/// (`Parser.pm:683-694`).
///
/// Documented: "A mobile element takes one colon before its name:
/// `<INS:ME::ALU>` is a `mobile_element_insertion`, `<INS:ME:ALU>` an
/// `Alu_insertion`." (How structural variants are written, The type of a
/// symbolic allele).
#[test]
fn mobile_element_insertions_carry_their_class_word_and_ensembl_veps_rows() {
    assert_word(
        "<INS:ME",
        &[
            ("<INS:ME:ALU>", "Alu_insertion"),
            ("<INS:ME:LINE1>", "LINE1_insertion"),
            ("<INS:ME:SVA>", "SVA_insertion"),
            ("<INS:ME:LINE>", "mobile_element_insertion"),
            ("<INS:ME>", "mobile_element_insertion"),
        ],
    );
}

/// `<DEL:ME:*>`: the deletion terms follow the same list, so `LINE1` names a
/// `LINE1_deletion`, and a bare `<DEL:ME>` or an element outside the list
/// (`HERVK`) is a `mobile_element_deletion` (`Parser.pm:683-694`).
///
/// Documented: "The type is read through the regexes of Ensembl VEP's
/// `get_SO_term`, so their edges hold" (How structural variants are written,
/// The type of a symbolic allele).
#[test]
fn mobile_element_deletions_carry_their_class_word_and_ensembl_veps_rows() {
    assert_word(
        "<DEL:ME",
        &[
            ("<DEL:ME:LINE1>", "LINE1_deletion"),
            ("<DEL:ME:HERVK>", "mobile_element_deletion"),
            ("<DEL:ME>", "mobile_element_deletion"),
        ],
    );
}

/// `<DUP:TANDEM>` is a `tandem_duplication` (`Parser.pm:697-698`), whose span
/// the page states: "A tandem repeat or tandem duplication ends at POS plus
/// the largest signed SVLEN" (How structural variants are written, Spans).
#[test]
fn tandem_duplications_carry_their_class_word_and_ensembl_veps_rows() {
    assert_word("<DUP:TANDEM", &[("<DUP:TANDEM>", "tandem_duplication")]);
}
