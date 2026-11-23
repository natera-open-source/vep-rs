// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Every genetic code of `Bio::Tools::CodonTable`, the module Ensembl VEP
//! translates with, held codon by codon to the module's own answers.
//!
//! `tests/data/genetic_codes.tsv` is the output of
//! `scripts/data/dump_genetic_codes.pl` inside the Ensembl VEP container
//! image: for each id from 0 to 40 that holds a table, one line a codon in the
//! module's codon order with the table id, the codon, its amino acid and 1 or
//! 0 for each of `is_start_codon` and `is_ter_codon`. Its `#` lines name the
//! BioPerl release, the ids that hold no table and each id the module answers
//! as another id.
//!
//! Genetic code 26, the nuclear code of Pachysolen tannophilus, is no table of
//! the module. It is held to the two strings NCBI's genetic code table gives
//! for it, written out here.

use std::collections::{BTreeMap, BTreeSet};

use vep_core::codon::{
    is_start_codon, is_ter_codon, known_table, module_table, translate_codon_with_table,
    translate_sequence,
};

const GENERATED: &str = include_str!("data/genetic_codes.tsv");

/// The genetic codes the ready-made Ensembl VEP caches name.
const CACHE_CODES: [u8; 11] = [1, 2, 3, 4, 5, 6, 9, 11, 12, 13, 26];

/// The one known genetic code the module holds no table for.
const PACHYSOLEN_NUCLEAR: u8 = 26;

/// Genetic code 26 as NCBI's genetic code table (`gc.prt`, version 4.6) gives
/// it, in the module's codon order: `ncbieaa`, the amino acid of each codon,
/// and `sncbieaa`, `M` at a start codon and `*` at a stop codon.
const NCBI_CODE_26: (&[u8; 64], &[u8; 64]) = (
    b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    b"----------**--*----M---------------M----------------------------",
);

/// One line of the generated table: a codon with the module's amino acid and
/// its `is_start_codon` and `is_ter_codon` answers.
#[derive(Debug, PartialEq)]
struct Row {
    codon: [u8; 3],
    amino_acid: u8,
    start: bool,
    stop: bool,
}

/// The generated table: the rows of each id that holds a table, the first and
/// last id asked, the ids that hold no table, and each id the module answers
/// as another id.
struct Generated {
    tables: BTreeMap<u8, Vec<Row>>,
    asked: (u8, u8),
    without_a_table: Vec<u8>,
    answered_as_another: Vec<(u8, u8)>,
}

/// The 64 codons in the module's order: the first base the most significant,
/// the bases counted `T`, `C`, `A`, `G`.
fn codons_in_module_order() -> Vec<[u8; 3]> {
    let bases = *b"TCAG";
    let mut codons = Vec::with_capacity(64);
    for first in bases {
        for second in bases {
            for third in bases {
                codons.push([first, second, third]);
            }
        }
    }
    codons
}

fn id(text: &str, line: usize) -> u8 {
    text.parse()
        .unwrap_or_else(|_| panic!("line {line}: {text:?} is no table id"))
}

fn mark(text: &str, line: usize) -> bool {
    match text {
        "1" => true,
        "0" => false,
        other => panic!("line {line}: mark {other:?} is neither 1 nor 0"),
    }
}

/// Reads the generated table, refusing a line of another shape and a table
/// whose codon column is not the 64 codons in the module's order.
fn generated() -> Generated {
    let mut tables: BTreeMap<u8, Vec<Row>> = BTreeMap::new();
    let mut asked = None;
    let mut without_a_table = None;
    let mut answered_as_another = None;
    for (index, text) in GENERATED.lines().enumerate() {
        let line = index + 1;
        if let Some(header) = text.strip_prefix("# ") {
            let (key, value) = header
                .split_once(':')
                .unwrap_or_else(|| panic!("line {line}: a header without a colon"));
            let value = value.trim();
            match key {
                "ids_asked" => {
                    let (first, last) = value
                        .split_once('-')
                        .unwrap_or_else(|| panic!("line {line}: {value:?} is no id range"));
                    asked = Some((id(first, line), id(last, line)));
                }
                "ids_without_a_table" => {
                    without_a_table = Some(value.split_whitespace().map(|v| id(v, line)).collect());
                }
                "ids_answered_as_another" => {
                    answered_as_another = Some(
                        value
                            .split_whitespace()
                            .map(|pair| {
                                let (asked_id, answered) =
                                    pair.split_once('=').unwrap_or_else(|| {
                                        panic!("line {line}: {pair:?} is no <id>=<id>")
                                    });
                                (id(asked_id, line), id(answered, line))
                            })
                            .collect(),
                    );
                }
                _ => {}
            }
            continue;
        }
        let fields: Vec<&str> = text.split('\t').collect();
        let [table, codon, amino_acid, start, stop] = fields.as_slice() else {
            panic!("line {line}: not five tab-separated columns");
        };
        let codon: [u8; 3] = codon
            .as_bytes()
            .try_into()
            .unwrap_or_else(|_| panic!("line {line}: the codon is not three bases"));
        let [amino_acid] = amino_acid.as_bytes() else {
            panic!("line {line}: the amino acid is not one letter");
        };
        tables.entry(id(table, line)).or_default().push(Row {
            codon,
            amino_acid: *amino_acid,
            start: mark(start, line),
            stop: mark(stop, line),
        });
    }
    let order = codons_in_module_order();
    for (table, rows) in &tables {
        let codons: Vec<[u8; 3]> = rows.iter().map(|row| row.codon).collect();
        assert_eq!(codons, order, "the codon column of table {table}");
    }
    Generated {
        tables,
        asked: asked.expect("the ids_asked header"),
        without_a_table: without_a_table.expect("the ids_without_a_table header"),
        answered_as_another: answered_as_another.expect("the ids_answered_as_another header"),
    }
}

/// How many of the 64 codons `table` and table 1 read differently: in amino
/// acid, in start mark and in stop mark.
fn differences_from_the_standard(table: u8) -> (usize, usize, usize) {
    let codons = codons_in_module_order();
    let count =
        |differs: &dyn Fn(&[u8; 3]) -> bool| codons.iter().filter(|codon| differs(codon)).count();
    (
        count(&|codon| {
            translate_codon_with_table(codon, table) != translate_codon_with_table(codon, 1)
        }),
        count(&|codon| is_start_codon(codon, table) != is_start_codon(codon, 1)),
        count(&|codon| is_ter_codon(codon, table) != is_ter_codon(codon, 1)),
    )
}

/// Under every table the module holds, each of the 64 codons translates to the
/// module's amino acid and is a start codon and a stop codon exactly where
/// `is_start_codon` and `is_ter_codon` say so; the 64 codons translated as one
/// sequence give the same amino acids.
#[test]
fn every_table_translates_every_codon_as_the_generated_table() {
    let generated = generated();
    assert!(!generated.tables.is_empty());
    for (&table, rows) in &generated.tables {
        for row in rows {
            let codon = String::from_utf8_lossy(&row.codon);
            assert_eq!(
                translate_codon_with_table(&row.codon, table) as char,
                row.amino_acid as char,
                "table {table}, codon {codon}"
            );
            assert_eq!(
                is_start_codon(&row.codon, table),
                row.start,
                "start mark of table {table}, codon {codon}"
            );
            assert_eq!(
                is_ter_codon(&row.codon, table),
                row.stop,
                "stop mark of table {table}, codon {codon}"
            );
        }
        let sequence: Vec<u8> = rows.iter().flat_map(|row| row.codon).collect();
        let amino_acids: Vec<u8> = rows.iter().map(|row| row.amino_acid).collect();
        assert_eq!(
            translate_sequence(&sequence, table),
            amino_acids,
            "table {table}"
        );
    }
}

/// The ids from 1 to 255 that [`known_table`] accepts are the ids the module
/// holds a table for and 26. Any other id is unknown: every codon reads `X`
/// under it and none is a start or a stop. The module answers id 0 as id 1,
/// and so does every function here.
#[test]
fn the_known_tables_are_the_modules_tables_and_code_26() {
    let generated = generated();
    let mut held: BTreeSet<u8> = generated
        .tables
        .keys()
        .copied()
        .filter(|&table| table != 0)
        .collect();
    assert!(!held.contains(&PACHYSOLEN_NUCLEAR));
    held.insert(PACHYSOLEN_NUCLEAR);
    let known: BTreeSet<u8> = (1..=u8::MAX).filter(|&table| known_table(table)).collect();
    assert_eq!(held, known);

    let (first, last) = generated.asked;
    for table in first..=last {
        assert_ne!(
            generated.tables.contains_key(&table),
            generated.without_a_table.contains(&table),
            "id {table} is listed once, with a table or without one"
        );
    }
    assert!(generated.without_a_table.contains(&PACHYSOLEN_NUCLEAR));
    let unknown = generated
        .without_a_table
        .iter()
        .filter(|&&table| table != PACHYSOLEN_NUCLEAR);
    for &table in unknown {
        assert!(!known_table(table), "id {table}");
        for codon in codons_in_module_order() {
            assert_eq!(
                translate_codon_with_table(&codon, table),
                b'X',
                "id {table}"
            );
            assert!(!is_start_codon(&codon, table), "id {table}");
            assert!(!is_ter_codon(&codon, table), "id {table}");
        }
    }

    assert_eq!(generated.answered_as_another, [(0, 1)]);
    assert_eq!(generated.tables[&0], generated.tables[&1]);
    assert!(known_table(0));
    for codon in codons_in_module_order() {
        assert_eq!(
            translate_codon_with_table(&codon, 0),
            translate_codon_with_table(&codon, 1)
        );
        assert_eq!(is_start_codon(&codon, 0), is_start_codon(&codon, 1));
        assert_eq!(is_ter_codon(&codon, 0), is_ter_codon(&codon, 1));
    }
}

/// Each genetic code a ready-made Ensembl VEP cache names is known, with at
/// least one start codon and one stop codon: a table the module holds, or
/// code 26.
#[test]
fn the_codes_of_the_shipped_caches_are_known() {
    let generated = generated();
    for table in CACHE_CODES {
        assert!(known_table(table), "table {table}");
        if table == PACHYSOLEN_NUCLEAR {
            assert!(!generated.tables.contains_key(&table));
            let codons = codons_in_module_order();
            assert!(codons.iter().any(|codon| is_start_codon(codon, table)));
            assert!(codons.iter().any(|codon| is_ter_codon(codon, table)));
            continue;
        }
        let rows = generated
            .tables
            .get(&table)
            .unwrap_or_else(|| panic!("the generated table lacks table {table}"));
        assert!(
            rows.iter()
                .any(|row| row.start && is_start_codon(&row.codon, table)),
            "table {table} has no start codon"
        );
        assert!(
            rows.iter()
                .any(|row| row.stop && is_ter_codon(&row.codon, table)),
            "table {table} has no stop codon"
        );
    }
}

/// No table but table 1 reads all 64 codons as table 1 does: each differs in
/// an amino acid, a start mark or a stop mark, in the module's answers and in
/// the functions here. Table 11 translates every codon as table 1 does and has
/// the same stop codons: it differs in its start codons alone.
#[test]
fn every_code_but_the_standard_differs_from_it_somewhere() {
    let generated = generated();
    let mut others = 0;
    for table in (2..=u8::MAX).filter(|&table| known_table(table)) {
        let (amino_acids, starts, stops) = differences_from_the_standard(table);
        assert!(
            amino_acids + starts + stops > 0,
            "table {table} reads every codon as table 1 does"
        );
        if table == PACHYSOLEN_NUCLEAR {
            continue;
        }
        assert_ne!(
            generated.tables[&table], generated.tables[&1],
            "table {table}"
        );
        others += 1;
    }
    assert_eq!(others + 2, generated.tables.len());
    assert_eq!(differences_from_the_standard(1), (0, 0, 0));

    let (amino_acids, starts, stops) = differences_from_the_standard(11);
    assert_eq!((amino_acids, stops), (0, 0));
    assert!(starts > 0);
    let (eleven, standard) = (&generated.tables[&11], &generated.tables[&1]);
    for (row, standard_row) in eleven.iter().zip(standard) {
        assert_eq!(row.amino_acid, standard_row.amino_acid);
        assert_eq!(row.stop, standard_row.stop);
    }
    assert!(eleven
        .iter()
        .zip(standard)
        .any(|(row, standard_row)| row.start != standard_row.start));
}

/// Genetic code 26 reads each of the 64 codons as NCBI's genetic code table
/// gives it: the amino acid, a start codon where `sncbieaa` holds `M`, a stop
/// codon where `ncbieaa` holds `*`. It differs from the standard code in one
/// amino acid, `CTG`, which codes alanine where the standard code reads
/// leucine, and in one start mark, `TTG`, which the standard code alone marks;
/// the stop codons are the same.
///
/// The module holds no table at 26: its `id` warns and falls to id 0, under
/// which `translate` reads the last table, 24 (`CodonTable.pm:328-338` and
/// `:412-460`), so [`module_table`] answers 24 for it and for every id past
/// it, and each id up to 24 as itself.
#[test]
fn genetic_code_26_is_ncbis_pachysolen_nuclear_code() {
    let (amino_acids, starts) = NCBI_CODE_26;
    let codons = codons_in_module_order();
    for (at, codon) in codons.iter().enumerate() {
        let text = String::from_utf8_lossy(codon);
        assert_eq!(
            translate_codon_with_table(codon, PACHYSOLEN_NUCLEAR) as char,
            amino_acids[at] as char,
            "codon {text}"
        );
        assert_eq!(
            is_start_codon(codon, PACHYSOLEN_NUCLEAR),
            starts[at] == b'M',
            "start mark of codon {text}"
        );
        assert_eq!(
            is_ter_codon(codon, PACHYSOLEN_NUCLEAR),
            amino_acids[at] == b'*',
            "stop mark of codon {text}"
        );
        assert_eq!(starts[at] == b'*', amino_acids[at] == b'*', "codon {text}");
    }
    let sequence: Vec<u8> = codons.iter().flatten().copied().collect();
    assert_eq!(
        translate_sequence(&sequence, PACHYSOLEN_NUCLEAR),
        amino_acids.to_vec()
    );

    assert_eq!(translate_codon_with_table(b"CTG", PACHYSOLEN_NUCLEAR), b'A');
    assert_eq!(translate_codon_with_table(b"CTG", 1), b'L');
    assert!(is_start_codon(b"TTG", 1) && !is_start_codon(b"TTG", PACHYSOLEN_NUCLEAR));
    assert_eq!(differences_from_the_standard(PACHYSOLEN_NUCLEAR), (1, 1, 0));

    assert_eq!(module_table(PACHYSOLEN_NUCLEAR), 24);
    for table in 0..=24 {
        assert_eq!(module_table(table), table, "id {table}");
    }
    for table in 25..=u8::MAX {
        assert_eq!(module_table(table), 24, "id {table}");
    }
    let generated = generated();
    let fallback = &generated.tables[&module_table(PACHYSOLEN_NUCLEAR)];
    let read_differently: Vec<String> = fallback
        .iter()
        .filter(|row| row.amino_acid != translate_codon_with_table(&row.codon, PACHYSOLEN_NUCLEAR))
        .map(|row| String::from_utf8_lossy(&row.codon).into_owned())
        .collect();
    assert_eq!(read_differently, ["CTG"]);
}
