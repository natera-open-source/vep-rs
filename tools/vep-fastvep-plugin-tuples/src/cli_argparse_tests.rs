// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line against `tests/fixtures/command_lines.tsv`, the table of what CPython's
//! `argparse` reads from each of its command lines, and against arguments outside UTF-8.

use super::*;

const TABLE: &str = include_str!("../tests/fixtures/command_lines.tsv");

/// The command lines the table holds; a table cut short fails here, not silently.
const COMMAND_LINES: usize = 1207;

/// An option of the table: `-` where it is absent, `=` and its value where it is given.
fn option(cell: &str) -> Option<Vec<u8>> {
    match cell {
        "-" => None,
        _ => Some(
            cell.strip_prefix('=')
                .expect("an option cell is `-` or starts with `=`")
                .as_bytes()
                .to_vec(),
        ),
    }
}

#[test]
fn every_command_line_of_the_table_reads_as_argparse_reads_it() {
    let mut read = 0;
    for line in TABLE.lines() {
        let cells: Vec<&str> = line.split('\t').collect();
        let argv: Vec<Vec<u8>> = cells[5..]
            .iter()
            .map(|arg| arg.as_bytes().to_vec())
            .collect();
        let expected = match cells[0] {
            "run" => Parsed::Run(Args {
                input: cells[1].as_bytes().to_vec(),
                tuples: cells[2].as_bytes().to_vec(),
                field_map: option(cells[3]),
                fields: option(cells[4]),
            }),
            "help" => Parsed::Help,
            "usage" => Parsed::Usage(cells[1].to_string()),
            other => panic!("an outcome is run, help or usage, not {other:?}"),
        };
        assert_eq!(parse(&argv), expected, "{:?}", &cells[5..]);
        read += 1;
    }
    assert_eq!(read, COMMAND_LINES);
}

#[test]
fn an_argument_outside_utf8_is_a_value_and_is_named_as_standard_error_writes_it() {
    let raw = [
        b"in\xff".to_vec(),
        b"--fields=\xffa,b".to_vec(),
        b"--field-map".to_vec(),
        b"./m\xe9/".to_vec(),
        b"\xfe".to_vec(),
    ];
    assert_eq!(
        parse(&raw),
        Parsed::Run(Args {
            input: b"in\xff".to_vec(),
            tuples: b"\xfe".to_vec(),
            field_map: Some(b"./m\xe9/".to_vec()),
            fields: Some(b"\xffa,b".to_vec()),
        })
    );
    let extra = [b"a".to_vec(), b"b".to_vec(), b"c\xff".to_vec()];
    assert_eq!(
        parse(&extra),
        Parsed::Usage("unrecognized arguments: c\\udcff".to_string())
    );
    let unknown = [b"a".to_vec(), b"b".to_vec(), b"--\xff".to_vec()];
    assert_eq!(
        parse(&unknown),
        Parsed::Usage("unrecognized arguments: --\\udcff".to_string())
    );
    let ignored = [b"--help=\xff'".to_vec()];
    assert_eq!(
        parse(&ignored),
        Parsed::Usage("argument -h/--help: ignored explicit argument \"\\udcff'\"".to_string())
    );
    let ambiguous = [b"--f=\xff".to_vec()];
    assert_eq!(
        parse(&ambiguous),
        Parsed::Usage(
            "ambiguous option: --f=\\udcff could match --field-map, --fields".to_string()
        )
    );
}
