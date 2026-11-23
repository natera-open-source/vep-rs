// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;

fn parse_bytes(args: &[&[u8]]) -> Result<Command, String> {
    use std::os::unix::ffi::OsStringExt;
    let argv: Vec<OsString> = args
        .iter()
        .map(|arg| OsString::from_vec(arg.to_vec()))
        .collect();
    parse(&argv)
}

fn parse_strs(args: &[&str]) -> Result<Command, String> {
    let args: Vec<&[u8]> = args.iter().map(|arg| arg.as_bytes()).collect();
    parse_bytes(&args)
}

/// The options of a command line that is accepted, after the two required ones.
fn accepted(options: &[&str]) -> Args {
    let args = [&["--input-dir", "in", "--output-dir", "out"], options].concat();
    match parse_strs(&args) {
        Ok(Command::Run(args)) => args,
        other => panic!("{options:?}: {other:?}"),
    }
}

/// The message a command line is refused with.
fn refused(args: &[&str]) -> String {
    match parse_strs(args) {
        Err(message) => message,
        other => panic!("{args:?}: {other:?}"),
    }
}

/// The message the options are refused with after the two required ones.
fn refused_after_required(options: &[&str]) -> String {
    refused(&[&["--input-dir", "in", "--output-dir", "out"], options].concat())
}

fn asks_for_help(args: &[&str]) -> bool {
    matches!(parse_strs(args), Ok(Command::Help))
}

#[test]
fn the_defaults_are_the_bare_filter_the_auto_mode_and_nothing_else() {
    let args = accepted(&[]);
    assert_eq!(args.input_dir, b"in");
    assert_eq!(args.output_dir, b"out");
    assert_eq!(args.chrom_mode, ChromMode::Auto);
    assert_eq!(args.cache_dir, None);
    assert_eq!(args.max_variants, None);
    assert_eq!(args.assembly, None);
    assert_eq!(args.canonical_contigs, ContigsFlag::Bare);
}

#[test]
fn every_option_takes_its_value_from_the_next_argument_or_after_an_equals_sign() {
    let args = accepted(&[
        "--chrom-mode",
        "strip",
        "--cache-dir=/cache",
        "--max-variants",
        "0012",
        "--assembly=",
        "--canonical-contigs",
        "1,2",
    ]);
    assert_eq!(args.chrom_mode, ChromMode::Strip);
    assert_eq!(args.cache_dir.as_deref(), Some(&b"/cache"[..]));
    let max = args.max_variants.unwrap();
    assert_eq!((max.decimal.as_str(), max.count), ("12", 12));
    assert_eq!(args.assembly.as_deref(), Some(&b""[..]));
    assert_eq!(args.canonical_contigs, ContigsFlag::Named(b"1,2".to_vec()));
    for (name, mode) in [
        ("auto", ChromMode::Auto),
        ("keep", ChromMode::Keep),
        ("strip", ChromMode::Strip),
        ("add", ChromMode::Add),
    ] {
        assert_eq!(accepted(&["--chrom-mode", name]).chrom_mode, mode);
        assert_eq!(mode.name(), name);
    }
    // A value keeps the bytes it was given in.
    let Ok(Command::Run(args)) = parse_bytes(&[
        b"--input-dir",
        b"in\xff",
        b"--output-dir=out\xfe",
        b"--assembly",
        b"\xe9",
    ]) else {
        panic!("a value outside UTF-8 is a value");
    };
    assert_eq!(
        (args.input_dir, args.output_dir, args.assembly),
        (
            b"in\xff".to_vec(),
            b"out\xfe".to_vec(),
            Some(b"\xe9".to_vec())
        )
    );
}

#[test]
fn a_later_value_replaces_an_earlier_one() {
    let args = accepted(&[
        "--input-dir",
        "other",
        "--max-variants=1",
        "--max-variants",
        "2",
    ]);
    assert_eq!(args.input_dir, b"other");
    assert_eq!(args.max_variants.unwrap().decimal, "2");
    let flag = |options: &[&str]| accepted(options).canonical_contigs;
    assert_eq!(flag(&["--no-canonical-contigs"]), ContigsFlag::Off);
    assert_eq!(
        flag(&["--no-canonical-contigs", "--canonical-contigs"]),
        ContigsFlag::Bare
    );
    assert_eq!(
        flag(&["--canonical-contigs", "1", "--no-canonical-contigs"]),
        ContigsFlag::Off
    );
    assert_eq!(
        flag(&["--no-canonical-contigs", "--canonical-contigs=X"]),
        ContigsFlag::Named(b"X".to_vec())
    );
}

#[test]
fn the_contig_option_takes_the_next_argument_only_where_it_is_not_an_option() {
    let flag = |args: &[&str]| match parse_strs(args) {
        Ok(Command::Run(args)) => args.canonical_contigs,
        other => panic!("{args:?}: {other:?}"),
    };
    let dirs = ["--input-dir", "in", "--output-dir", "out"];
    assert_eq!(
        flag(&[
            "--canonical-contigs",
            "--input-dir",
            "in",
            "--output-dir",
            "out"
        ]),
        ContigsFlag::Bare
    );
    assert_eq!(
        flag(&[&["--canonical-contigs", "-5"][..], &dirs].concat()),
        ContigsFlag::Named(b"-5".to_vec())
    );
    assert_eq!(
        flag(&[&["--canonical-contigs", "-x y"][..], &dirs].concat()),
        ContigsFlag::Named(b"-x y".to_vec())
    );
    assert_eq!(
        flag(&[&["--canonical-contigs", "-"][..], &dirs].concat()),
        ContigsFlag::Named(b"-".to_vec())
    );
    assert_eq!(
        flag(&[&["--canonical-contigs="][..], &dirs].concat()),
        ContigsFlag::Named(Vec::new())
    );
    assert_eq!(
        flag(&[&["--canonical-contigs=--"][..], &dirs].concat()),
        ContigsFlag::Named(b"--".to_vec())
    );
    assert_eq!(
        flag(
            &[
                &dirs[..],
                &["--canonical-contigs", "--no-canonical-contigs"]
            ]
            .concat()
        ),
        ContigsFlag::Off
    );
    // An unknown option after it is not its value, and neither is the marker.
    assert_eq!(
        refused_after_required(&["--canonical-contigs", "-x"]),
        "unrecognized arguments: -x"
    );
    assert_eq!(
        refused_after_required(&["--canonical-contigs", "--", "1,2"]),
        "unrecognized arguments: -- 1,2"
    );
}

#[test]
fn an_option_is_named_by_any_prefix_that_fits_it_alone() {
    let Ok(Command::Run(args)) = parse_strs(&[
        "--i", "in", "--o=out", "--ch", "add", "--cac", "c", "--m=3", "--a", "x", "--n",
    ]) else {
        panic!("every prefix names one option");
    };
    assert_eq!(
        (args.input_dir, args.output_dir),
        (b"in".to_vec(), b"out".to_vec())
    );
    assert_eq!(args.chrom_mode, ChromMode::Add);
    assert_eq!(args.cache_dir, Some(b"c".to_vec()));
    assert_eq!(args.max_variants.unwrap().decimal, "3");
    assert_eq!(args.assembly, Some(b"x".to_vec()));
    assert_eq!(args.canonical_contigs, ContigsFlag::Off);
    assert_eq!(
        accepted(&["--can", "1"]).canonical_contigs,
        ContigsFlag::Named(b"1".to_vec())
    );
    assert_eq!(
        refused_after_required(&["--c", "x"]),
        "ambiguous option: --c could match --chrom-mode, --cache-dir, --canonical-contigs"
    );
    assert_eq!(
        refused_after_required(&["--ca=x"]),
        "ambiguous option: --ca=x could match --cache-dir, --canonical-contigs"
    );
    assert_eq!(
        refused_after_required(&["--=x"]),
        "ambiguous option: --=x could match --help, --input-dir, --output-dir, --chrom-mode, \
         --cache-dir, --max-variants, --assembly, --canonical-contigs, --no-canonical-contigs"
    );
    assert_eq!(
        refused_after_required(&["-=x"]),
        "ambiguous option: -=x could match -h, --help, --input-dir, --output-dir, --chrom-mode, \
         --cache-dir, --max-variants, --assembly, --canonical-contigs, --no-canonical-contigs"
    );
}

#[test]
fn errors_come_in_the_order_the_arguments_are_read_then_required_then_unread() {
    assert_eq!(
        refused(&[]),
        "the following arguments are required: --input-dir, --output-dir"
    );
    assert_eq!(
        refused(&["--output-dir", "out"]),
        "the following arguments are required: --input-dir"
    );
    assert_eq!(
        refused(&["stray", "--input-dir", "in", "-x"]),
        "the following arguments are required: --output-dir"
    );
    assert_eq!(
        refused(&[
            "one",
            "--input-dir",
            "in",
            "--frob=1",
            "--output-dir",
            "out",
            "two",
            ""
        ]),
        "unrecognized arguments: one --frob=1 two "
    );
    assert_eq!(
        refused_after_required(&["--", "--max-variants", "x"]),
        "unrecognized arguments: -- --max-variants x"
    );
    assert_eq!(
        refused(&["--", "--input-dir", "in", "--output-dir", "out"]),
        "the following arguments are required: --input-dir, --output-dir"
    );
    assert_eq!(
        refused(&["--max-variants", "x", "--chrom-mode", "y"]),
        "argument --max-variants: invalid int value: 'x'"
    );
    assert_eq!(
        refused(&["--chrom-mode", "y", "--max-variants", "x"]),
        "argument --chrom-mode: invalid choice: 'y' (choose from 'auto', 'keep', 'strip', 'add')"
    );
}

#[test]
fn a_value_is_refused_in_reprs_words() {
    assert_eq!(
        refused(&["--chrom-mode", "it's"]),
        "argument --chrom-mode: invalid choice: \"it's\" (choose from 'auto', 'keep', 'strip', 'add')"
    );
    assert_eq!(
        refused(&["--chrom-mode="]),
        "argument --chrom-mode: invalid choice: '' (choose from 'auto', 'keep', 'strip', 'add')"
    );
    assert_eq!(
        refused(&["--ch", "Auto"]),
        "argument --chrom-mode: invalid choice: 'Auto' (choose from 'auto', 'keep', 'strip', 'add')"
    );
    for (value, shown) in [
        ("many", "'many'"),
        ("1.5", "'1.5'"),
        ("", "''"),
        ("1__0", "'1__0'"),
        ("\u{1c}2", "'\\x1c2'"),
        ("a\tb", "'a\\tb'"),
    ] {
        assert_eq!(
            refused(&["--max-variants", value]),
            format!("argument --max-variants: invalid int value: {shown}")
        );
    }
    assert_eq!(
        parse_bytes(&[b"--max-variants", b"\xff"]).unwrap_err(),
        "argument --max-variants: invalid int value: '\\udcff'"
    );
    assert_eq!(
        refused(&["--no-canonical-contigs=x"]),
        "argument --no-canonical-contigs: ignored explicit argument 'x'"
    );
    assert_eq!(
        refused(&["--no="]),
        "argument --no-canonical-contigs: ignored explicit argument ''"
    );
}

#[test]
fn a_value_that_starts_like_an_option_is_not_taken() {
    assert_eq!(
        refused(&["--input-dir", "--output-dir", "out"]),
        "argument --input-dir: expected one argument"
    );
    assert_eq!(
        refused(&["--input-dir", "in", "--output-dir"]),
        "argument --output-dir: expected one argument"
    );
    assert_eq!(
        refused(&["--input-dir", "-x", "--output-dir", "out"]),
        "argument --input-dir: expected one argument"
    );
    assert_eq!(
        refused(&["--input-dir", "--", "--output-dir", "out"]),
        "argument --input-dir: expected one argument"
    );
    // `-` alone, a negative number's start in the digits of any script and an argument
    // with a space are values.
    for value in [
        "-",
        "-5",
        "-5x",
        "-.5",
        "-x y",
        "",
        "-\u{663}",
        "-.\u{ff15}x",
    ] {
        let Ok(Command::Run(args)) = parse_strs(&["--input-dir", value, "--output-dir", "out"])
        else {
            panic!("{value:?} is a value");
        };
        assert_eq!(args.input_dir, value.as_bytes());
    }
    for value in ["-3", "-\u{663}"] {
        assert_eq!(
            accepted(&["--max-variants", value])
                .max_variants
                .unwrap()
                .decimal,
            "-3",
            "{value}"
        );
    }
    // A numeral that is no decimal digit, and a dot that no digit follows, start no number.
    for value in ["-\u{b2}", "-..5", "-."] {
        assert_eq!(
            refused(&["--input-dir", value, "--output-dir", "out"]),
            "argument --input-dir: expected one argument",
            "{value}"
        );
    }
}

#[test]
fn help_is_asked_for_where_it_is_read_before_an_error() {
    for args in [
        &["-h"][..],
        &["--help"],
        &["--he"],
        &["--h"],
        &["-hh"],
        &["-hx"],
        &["-hhx"],
        &["-hx=y"],
        &["-h", "--max-variants", "x"],
        &["-h", "--c"],
        &["-x", "-h"],
        &["stray", "--input-dir", "in", "-h"],
    ] {
        assert!(asks_for_help(args), "{args:?}");
    }
    assert_eq!(
        refused(&["--max-variants", "x", "-h"]),
        "argument --max-variants: invalid int value: 'x'"
    );
    assert_eq!(
        refused(&["--c", "-h"]),
        "ambiguous option: --c could match --chrom-mode, --cache-dir, --canonical-contigs"
    );
    assert_eq!(
        refused(&["--", "-h"]),
        "the following arguments are required: --input-dir, --output-dir"
    );
    for (arg, value) in [
        ("--help=x", "'x'"),
        ("--help=", "''"),
        ("-h=", "''"),
        ("-h=x", "'x'"),
        ("-h-x", "'-x'"),
        ("-hh=x", "'x'"),
        ("-hh=", "''"),
    ] {
        assert_eq!(
            refused(&[arg]),
            format!("argument -h/--help: ignored explicit argument {value}"),
            "{arg}"
        );
    }
}
