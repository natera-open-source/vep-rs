// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line as `argparse` reads it: the options and defaults of each subcommand, the
//! spellings an option takes, and the error line of each refused command line after the
//! usage of the parser that refuses it.

use super::*;

fn read(args: &[&str]) -> Parsed {
    let args: Vec<String> = args.iter().map(|arg| arg.to_string()).collect();
    parse(&args)
}

/// The error line of a refused command line, without the usage block before it.
fn refused(args: &[&str]) -> String {
    let Parsed::Refused(text) = read(args) else {
        panic!("{args:?} is not refused");
    };
    let (usage, line) = text
        .trim_end_matches('\n')
        .rsplit_once('\n')
        .expect("a usage block and an error line");
    let parser = line.split(": error: ").next().unwrap();
    let help = match parser {
        "vep-golden-corpus" => TOP.help,
        "vep-golden-corpus select" => SELECT.help,
        "vep-golden-corpus classify" => CLASSIFY.help,
        "vep-golden-corpus check" => CHECK.help,
        "vep-golden-corpus identifiers" => IDENTIFIERS.help,
        "vep-golden-corpus covers" => COVERS.help,
        "vep-golden-corpus flag-set" => FLAG_SET.help,
        "vep-golden-corpus generated" => GENERATED.help,
        other => panic!("no parser {other:?}"),
    };
    assert!(
        help.starts_with(&format!("{usage}\n\n")),
        "{args:?} prints the usage of {parser}"
    );
    line.to_string()
}

fn ran(command: Command) -> Parsed {
    Parsed::Run(Box::new(command))
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| value.to_string()).collect()
}

#[test]
fn select_reads_every_option_and_its_default() {
    let minimal = read(&[
        "select",
        "--release",
        "116",
        "--assembly",
        "GRCh38",
        "--dataset",
        "s01=ref.txt=in.vcf",
        "--out",
        "corpus",
    ]);
    let defaults = SelectArgs {
        release: "116".to_string(),
        assembly: "GRCh38".to_string(),
        species: None,
        dataset: strings(&["s01=ref.txt=in.vcf"]),
        out: "corpus".to_string(),
        k: 3,
        seed: 1,
        max_records: 1500,
        multi_allelic_per_dataset: 4,
        max_span: 1_000_000,
        focus_term: Vec::new(),
        focus_k: None,
        include: Vec::new(),
        keep_records: false,
        cache: None,
        flank: 5000,
        covers: Vec::new(),
        family_columns: None,
        keep_samples: 0,
        genotype_exemplars: 8,
        require_genotype_class: Vec::new(),
        record_classes: None,
        record_class_exemplars: 2,
        matrix: None,
        drop_reference_statements: false,
        flag: Vec::new(),
        fasta_name: None,
    };
    assert_eq!(minimal, ran(Command::Select(Box::new(defaults))));
    let full = read(&[
        "select",
        "--release=116",
        "--assembly",
        "GRCh38",
        "--species",
        "mus_musculus",
        "--dataset",
        "a=b=c",
        "--dataset",
        "d=e=f",
        "--out",
        "first",
        "--out",
        "corpus",
        "--k",
        " 2 ",
        "--seed",
        "-7",
        "--max-records",
        "1_000",
        "--multi-allelic-per-dataset",
        "0",
        "--max-span",
        "+5",
        "--focus-term",
        "stop_lost",
        "--focus-term",
        "stop_gained",
        "--focus-k",
        "9",
        "--include",
        "a=rs1",
        "--keep-records",
        "--cache",
        "cache",
        "--flank",
        "10",
        "--covers",
        "symbol",
        "--family-columns",
        "HGVSc,IMPACT=",
        "--keep-samples",
        "2",
        "--genotype-exemplars",
        "3",
        "--require-genotype-class",
        "het",
        "--record-classes",
        "regulatory_class",
        "--record-class-exemplars",
        "4",
        "--matrix",
        "matrix.tsv",
        "--drop-reference-statements",
        "--flag=--hgvs",
        "--flag=--symbol",
        "--fasta-name",
        "reference.fa.gz",
    ]);
    let given = SelectArgs {
        release: "116".to_string(),
        assembly: "GRCh38".to_string(),
        species: Some("mus_musculus".to_string()),
        dataset: strings(&["a=b=c", "d=e=f"]),
        out: "corpus".to_string(),
        k: 2,
        seed: -7,
        max_records: 1000,
        multi_allelic_per_dataset: 0,
        max_span: 5,
        focus_term: strings(&["stop_lost", "stop_gained"]),
        focus_k: Some(9),
        include: strings(&["a=rs1"]),
        keep_records: true,
        cache: Some("cache".to_string()),
        flank: 10,
        covers: strings(&["symbol"]),
        family_columns: Some("HGVSc,IMPACT=".to_string()),
        keep_samples: 2,
        genotype_exemplars: 3,
        require_genotype_class: strings(&["het"]),
        record_classes: Some("regulatory_class".to_string()),
        record_class_exemplars: 4,
        matrix: Some("matrix.tsv".to_string()),
        drop_reference_statements: true,
        flag: strings(&["--hgvs", "--symbol"]),
        fasta_name: Some("reference.fa.gz".to_string()),
    };
    assert_eq!(full, ran(Command::Select(Box::new(given))));
}

#[test]
fn classify_and_check_read_their_options() {
    assert_eq!(
        read(&[
            "classify",
            "--corpus",
            "corpus",
            "--vep-default",
            "vep.txt",
            "--vep-rs-output",
            "vep_rs.txt",
            "--vep-vcf",
            "vep.vcf",
            "--vep-rs-vcf=vep_rs.vcf",
            "--vep-rs-j",
            "vep_rs.json",
            "--field",
            "HGVSc",
            "--field",
            "HGVSp",
            "--flag-set",
            "hgvs",
            "--spliceai-table",
            "snv.synthetic.vcf.gz",
            "--reference-fasta",
            "reference.fa",
            "--formats-comparator",
            "vep-formats-compare",
            "--flag=--hgvs",
            "--fasta-name",
            "reference.fa.gz",
        ]),
        ran(Command::Classify(ClassifyArgs {
            corpus: "corpus".to_string(),
            vep_default: "vep.txt".to_string(),
            vep_rs_output: "vep_rs.txt".to_string(),
            vep_tab: None,
            vep_vcf: Some("vep.vcf".to_string()),
            vep_json: None,
            vep_rs_tab: None,
            vep_rs_vcf: Some("vep_rs.vcf".to_string()),
            vep_rs_json: Some("vep_rs.json".to_string()),
            field: strings(&["HGVSc", "HGVSp"]),
            flag_set: Some("hgvs".to_string()),
            spliceai_table: strings(&["snv.synthetic.vcf.gz"]),
            reference_fasta: Some("reference.fa".to_string()),
            formats_comparator: Some("vep-formats-compare".to_string()),
            flag: strings(&["--hgvs"]),
            fasta_name: Some("reference.fa.gz".to_string()),
        }))
    );
    assert_eq!(
        read(&["check", "--corpus", "corpus", "--dataset", "a=b=c"]),
        ran(Command::Check(CheckArgs {
            corpus: "corpus".to_string(),
            dataset: strings(&["a=b=c"]),
            matrix: None,
            family_columns: None,
            require_genotype_class: Vec::new(),
        }))
    );
    assert_eq!(
        read(&[
            "check",
            "--c=corpus",
            "--d",
            "a=b=c",
            "--m",
            "",
            "--f",
            "HGVSc",
            "--r",
            "het",
            "--r",
            "star",
        ]),
        ran(Command::Check(CheckArgs {
            corpus: "corpus".to_string(),
            dataset: strings(&["a=b=c"]),
            matrix: Some(String::new()),
            family_columns: Some("HGVSc".to_string()),
            require_genotype_class: strings(&["het", "star"]),
        }))
    );
}

#[test]
fn covers_and_flag_set_read_their_corpus_among_their_options() {
    let listed = CoversArgs {
        corpus: "corpus".to_string(),
        add: Vec::new(),
        list: true,
        flag_set: None,
        matrix: "tests/parity/matrix.tsv".to_string(),
    };
    for args in [
        &["covers", "corpus", "--list"][..],
        &["covers", "--list", "corpus"],
        &["covers", "--li", "corpus"],
        &["covers", "--l", "--", "corpus"],
        &["--", "covers", "corpus", "--list", "--list"],
    ] {
        assert_eq!(read(args), ran(Command::Covers(listed.clone())), "{args:?}");
    }
    // A value that starts like a negative number or holds a space is taken as a value.
    assert_eq!(
        read(&[
            "covers",
            "-1",
            "--add",
            "-.5",
            "--add=",
            "--a",
            "-x y",
            "--add",
            "-",
            "--flag-set",
            "",
            "--matrix",
            "m.tsv",
        ]),
        ran(Command::Covers(CoversArgs {
            corpus: "-1".to_string(),
            add: strings(&["-.5", "", "-x y", "-"]),
            list: false,
            flag_set: Some(String::new()),
            matrix: "m.tsv".to_string(),
        }))
    );
    assert_eq!(
        read(&["flag-set", "corpus"]),
        ran(Command::FlagSet(FlagSetArgs {
            corpus: "corpus".to_string(),
            name: None,
            flag: Vec::new(),
            input: "variants.vcf".to_string(),
            stdin: false,
            unordered: false,
            perl_mode: "offline".to_string(),
            format: None,
            expected_dir: None,
            covers: Vec::new(),
            matrix: "tests/parity/matrix.tsv".to_string(),
            list: false,
            drop_reference_statements: false,
        }))
    );
    assert_eq!(
        read(&[
            "flag-set",
            "--name",
            "pick",
            "--flag=--pick",
            "--fl=--format=vcf",
            "--input",
            "variants.ensembl",
            "--st",
            "--u",
            "--perl",
            "gtf",
            "--fo",
            "default,json",
            "--e",
            "expected_ens",
            "--c",
            "pick",
            "--covers",
            "per_gene",
            "--m",
            "m.tsv",
            "--l",
            "--d",
            "corpus",
        ]),
        ran(Command::FlagSet(FlagSetArgs {
            corpus: "corpus".to_string(),
            name: Some("pick".to_string()),
            flag: strings(&["--pick", "--format=vcf"]),
            input: "variants.ensembl".to_string(),
            stdin: true,
            unordered: true,
            perl_mode: "gtf".to_string(),
            format: Some("default,json".to_string()),
            expected_dir: Some("expected_ens".to_string()),
            covers: strings(&["pick", "per_gene"]),
            matrix: "m.tsv".to_string(),
            list: true,
            drop_reference_statements: true,
        }))
    );
}

#[test]
fn identifiers_reads_its_paths_up_to_the_next_option() {
    let named = |path: &[&str], check: bool| {
        ran(Command::Identifiers(IdentifiersArgs {
            path: strings(path),
            check,
        }))
    };
    let cases: [(&[&str], &[&str], bool); 13] = [
        (&["identifiers", "d1"], &["d1"], false),
        (&["identifiers", "--check", "d1", "d2"], &["d1", "d2"], true),
        (&["identifiers", "d1", "d2", "--check"], &["d1", "d2"], true),
        (&["identifiers", "--che", "d1"], &["d1"], true),
        (&["identifiers", "--c", "d1"], &["d1"], true),
        (&["identifiers", "d1", "--check", "--check"], &["d1"], true),
        // The first `--` among the paths is dropped, and what follows it is a path.
        (&["identifiers", "--", "d1", "d2"], &["d1", "d2"], false),
        (&["identifiers", "d1", "--", "d2"], &["d1", "d2"], false),
        (&["identifiers", "d1", "--"], &["d1"], false),
        (&["identifiers", "--check", "--", "d1"], &["d1"], true),
        (
            &["identifiers", "d1", "--", "--check"],
            &["d1", "--check"],
            false,
        ),
        (
            &["identifiers", "d1", "--", "--", "d2"],
            &["d1", "--", "d2"],
            false,
        ),
        (&["identifiers", "-5", "a b"], &["-5", "a b"], false),
    ];
    for (args, path, check) in cases {
        assert_eq!(read(args), named(path, check), "{args:?}");
    }
    let own = "vep-golden-corpus identifiers: error: ";
    let top = "vep-golden-corpus: error: ";
    let refusals: [(&[&str], String); 8] = [
        (
            &["identifiers"],
            format!("{own}the following arguments are required: path"),
        ),
        (
            &["identifiers", "--check"],
            format!("{own}the following arguments are required: path"),
        ),
        (
            &["identifiers", "--"],
            format!("{own}the following arguments are required: path"),
        ),
        (
            &["identifiers", "--check=1", "d1"],
            format!("{own}argument --check: ignored explicit argument '1'"),
        ),
        // The paths end at the first option: one after it is taken by nothing.
        (
            &["identifiers", "d1", "--check", "d2"],
            format!("{top}unrecognized arguments: d2"),
        ),
        (
            &["identifiers", "d1", "-x"],
            format!("{top}unrecognized arguments: -x"),
        ),
        (
            &["identifiers", "-c", "d1"],
            format!("{top}unrecognized arguments: -c"),
        ),
        (
            &["identifiers", "--checks", "d1"],
            format!("{top}unrecognized arguments: --checks"),
        ),
    ];
    for (args, line) in refusals {
        assert_eq!(refused(args), line, "{args:?}");
    }
}

#[test]
fn help_is_shown_by_the_parser_that_reads_it() {
    assert_eq!(read(&["-h"]), Parsed::Help(TOP.help));
    assert_eq!(read(&["--help"]), Parsed::Help(TOP.help));
    assert_eq!(read(&["--he", "covers"]), Parsed::Help(TOP.help));
    // A letter behind `-h` that names no option is left unread and the help is shown.
    assert_eq!(read(&["-hx"]), Parsed::Help(TOP.help));
    assert_eq!(read(&["-hh"]), Parsed::Help(TOP.help));
    assert_eq!(read(&["covers", "-h"]), Parsed::Help(COVERS.help));
    assert_eq!(
        read(&["covers", "corpus", "--add", "x", "--help"]),
        Parsed::Help(COVERS.help)
    );
    assert_eq!(read(&["select", "-h"]), Parsed::Help(SELECT.help));
    assert_eq!(read(&["classify", "--h"]), Parsed::Help(CLASSIFY.help));
    assert_eq!(read(&["check", "-h", "--nope"]), Parsed::Help(CHECK.help));
    assert_eq!(
        read(&["identifiers", "--check", "-h", "corpus"]),
        Parsed::Help(IDENTIFIERS.help)
    );
    assert_eq!(
        read(&["identifiers", "corpus", "--help", "extra"]),
        Parsed::Help(IDENTIFIERS.help)
    );
    assert_eq!(
        read(&["flag-set", "--nope", "-h"]),
        Parsed::Help(FLAG_SET.help)
    );
    let parsers = [
        &TOP,
        &SELECT,
        &CLASSIFY,
        &CHECK,
        &IDENTIFIERS,
        &COVERS,
        &FLAG_SET,
        &GENERATED,
    ];
    for parser in parsers {
        let opening = format!("usage: {} [-h]", parser.prog);
        assert!(parser.usage().starts_with(&opening));
        assert!(parser.usage().ends_with('\n') && !parser.usage().contains("\n\n"));
    }
}

#[test]
fn a_refused_command_line_ends_on_the_line_argparse_writes() {
    let top = "vep-golden-corpus: error: ";
    let choices =
        "(choose from 'select', 'classify', 'check', 'identifiers', 'covers', 'flag-set', 'generated')";
    let required = "the following arguments are required: cmd";
    let ignored = "argument -h/--help: ignored explicit argument";
    let cases: Vec<(&[&str], String)> = vec![
        (&[], format!("{top}{required}")),
        (&["-x"], format!("{top}{required}")),
        (&["--"], format!("{top}{required}")),
        (&["-h=x"], format!("{top}{ignored} 'x'")),
        (&["-hh=x"], format!("{top}{ignored} 'x'")),
        (&["--help="], format!("{top}{ignored} ''")),
        (&["-h-"], format!("{top}{ignored} '-'")),
        (
            &["nonsense"],
            format!("{top}argument cmd: invalid choice: 'nonsense' {choices}"),
        ),
        (
            &["-1", "covers", "c"],
            format!("{top}argument cmd: invalid choice: '-1' {choices}"),
        ),
        (
            &["-a b"],
            format!("{top}argument cmd: invalid choice: '-a b' {choices}"),
        ),
        (
            &["--", "--", "covers", "c"],
            format!("{top}argument cmd: invalid choice: '--' {choices}"),
        ),
        (
            &["--foo", "covers", "c"],
            format!("{top}unrecognized arguments: --foo"),
        ),
        (
            &["covers", "c", "--foo"],
            format!("{top}unrecognized arguments: --foo"),
        ),
        (
            &["covers", "c", "extra", "more", "--list"],
            format!("{top}unrecognized arguments: extra more"),
        ),
        (
            &["covers", "c", "--nope=1", "--list"],
            format!("{top}unrecognized arguments: --nope=1"),
        ),
        (
            &["covers", "c", "--", "--list"],
            format!("{top}unrecognized arguments: --list"),
        ),
        (
            &["covers", "--", "--", "c"],
            format!("{top}unrecognized arguments: c"),
        ),
        (
            &["covers", "--list", "--", "c", "--"],
            format!("{top}unrecognized arguments: --"),
        ),
        (
            &["covers", "c", "-x"],
            format!("{top}unrecognized arguments: -x"),
        ),
        (
            &["covers", "c", "-"],
            format!("{top}unrecognized arguments: -"),
        ),
        (
            &["covers", "c", "--flag-sets", "x"],
            format!("{top}unrecognized arguments: --flag-sets x"),
        ),
        (
            &[
                "select",
                "--release",
                "1",
                "--assembly",
                "a",
                "--dataset",
                "s",
                "--out",
                "o",
                "--nope",
                "-q",
            ],
            format!("{top}unrecognized arguments: --nope -q"),
        ),
    ];
    for (args, expected) in cases {
        assert_eq!(refused(args), expected, "{args:?}");
    }
}

#[test]
fn a_subcommand_refuses_with_its_own_usage() {
    let cases: [(&[&str], &str); 35] = [
        (
            &["covers"],
            "covers: error: the following arguments are required: corpus",
        ),
        (
            &["--foo", "covers"],
            "covers: error: the following arguments are required: corpus",
        ),
        (
            &["covers", "c", "--add"],
            "covers: error: argument --add: expected one argument",
        ),
        (
            &["covers", "c", "--add", "--list"],
            "covers: error: argument --add: expected one argument",
        ),
        (
            &["covers", "c", "--add", "-h"],
            "covers: error: argument --add: expected one argument",
        ),
        (
            &["covers", "c", "--add", "-x"],
            "covers: error: argument --add: expected one argument",
        ),
        (
            &["covers", "c", "--add", "--"],
            "covers: error: argument --add: expected one argument",
        ),
        (
            &["covers", "c", "--list=x"],
            "covers: error: argument --list: ignored explicit argument 'x'",
        ),
        (
            &["covers", "c", "--list="],
            "covers: error: argument --list: ignored explicit argument ''",
        ),
        (
            &["covers", "c", "-=x"],
            "covers: error: ambiguous option: -=x could match -h, --help, --add, --list, --flag-set, --matrix",
        ),
        (
            &["covers", "c", "--=x"],
            "covers: error: ambiguous option: --=x could match --help, --add, --list, --flag-set, --matrix",
        ),
        (
            &["select"],
            "select: error: the following arguments are required: --release, --assembly, --dataset, --out",
        ),
        (
            &["select", "positional"],
            "select: error: the following arguments are required: --release, --assembly, --dataset, --out",
        ),
        (
            &["select", "--out", "o", "--dataset", "s"],
            "select: error: the following arguments are required: --release, --assembly",
        ),
        (
            &["select", "--release", "--assembly", "x"],
            "select: error: argument --release: expected one argument",
        ),
        (
            &["select", "--k", "1.5"],
            "select: error: argument --k: invalid int value: '1.5'",
        ),
        (
            &["select", "--k=", "--release", "1"],
            "select: error: argument --k: invalid int value: ''",
        ),
        (
            &["select", "--seed", "1__0"],
            "select: error: argument --seed: invalid int value: '1__0'",
        ),
        (
            &["select", "--k", "-1", "--seed", "x"],
            "select: error: argument --seed: invalid int value: 'x'",
        ),
        (
            &["select", "--k", "-x"],
            "select: error: argument --k: expected one argument",
        ),
        (
            &["select", "--keep-records", "--k", "it's"],
            "select: error: argument --k: invalid int value: \"it's\"",
        ),
        (
            &["select", "--k", "tab\there"],
            "select: error: argument --k: invalid int value: 'tab\\there'",
        ),
        (
            &["select", "--k", "\u{a0}x"],
            "select: error: argument --k: invalid int value: '\\xa0x'",
        ),
        (
            &["select", "--keep-records=1"],
            "select: error: argument --keep-records: ignored explicit argument '1'",
        ),
        (
            &["select", "--keep", "1"],
            "select: error: ambiguous option: --keep could match --keep-records, --keep-samples",
        ),
        (
            &["select", "--m", "5"],
            "select: error: ambiguous option: --m could match --max-records, --multi-allelic-per-dataset, --max-span, --matrix",
        ),
        (
            &["select", "--f", "5"],
            "select: error: ambiguous option: --f could match --focus-term, --focus-k, --flank, --family-columns, --flag, --fasta-name",
        ),
        (
            &["classify", "--fl", "x"],
            "classify: error: ambiguous option: --fl could match --flag-set, --flag",
        ),
        (
            &["classify", "--vep-rs", "x"],
            "classify: error: ambiguous option: --vep-rs could match --vep-rs-output, --vep-rs-tab, --vep-rs-vcf, --vep-rs-json",
        ),
        (
            &["classify", "--vep-d", "x"],
            "classify: error: the following arguments are required: --corpus, --vep-rs-output",
        ),
        (
            &["classify", "--s", "x", "--c"],
            "classify: error: argument --corpus: expected one argument",
        ),
        (
            &["check", "--dataset", "s"],
            "check: error: the following arguments are required: --corpus",
        ),
        (
            &["flag-set", "c", "--perl", "it's"],
            "flag-set: error: argument --perl-mode: invalid choice: \"it's\" (choose from 'offline', 'database', 'gtf', 'gff')",
        ),
        (
            &["flag-set", "c", "--f", "x"],
            "flag-set: error: ambiguous option: --f could match --flag, --format",
        ),
        (
            &["flag-set", "c", "--name", "-x"],
            "flag-set: error: argument --name: expected one argument",
        ),
    ];
    for (args, expected) in cases {
        assert_eq!(
            refused(args),
            format!("vep-golden-corpus {expected}"),
            "{args:?}"
        );
    }
}

#[test]
fn the_first_refusal_or_help_from_the_left_ends_the_reading() {
    // An option is acted on as it is read: a bad count before `-h` is refused, one after it
    // is never read, and a missing required option is named only when every argument is read.
    assert_eq!(
        refused(&["select", "--k", "x", "-h"]),
        "vep-golden-corpus select: error: argument --k: invalid int value: 'x'"
    );
    assert_eq!(
        read(&["select", "-h", "--k", "x"]),
        Parsed::Help(SELECT.help)
    );
    assert_eq!(
        refused(&["select", "--keep-r", "--k", "x"]),
        "vep-golden-corpus select: error: argument --k: invalid int value: 'x'"
    );
    // The subcommand's refusal comes before the arguments nothing took.
    assert_eq!(
        refused(&["--foo", "select", "--nope"]),
        "vep-golden-corpus select: error: the following arguments are required: --release, --assembly, --dataset, --out"
    );
}

#[test]
fn generated_reads_no_abbreviation_and_sets_the_words_of_its_command_aside() {
    let generated = |check: bool, command: Option<&[&str]>| {
        ran(Command::Generated(GeneratedArgs {
            corpus: "c".to_string(),
            check,
            command: command.map(strings),
        }))
    };
    assert_eq!(
        read(&["generated", "--corpus", "c", "--check"]),
        generated(true, None)
    );
    assert_eq!(read(&["generated", "--corpus=c"]), generated(false, None));
    // Every word after `--command` is the command's: an option of the subcommand, `-h`, a
    // `--` and a second `--command` among them.
    let words = ["a", "--", "b", "--check", "-h", "--command"];
    let mut args = vec!["generated", "--corpus", "c", "--command"];
    args.extend(words);
    assert_eq!(read(&args), generated(false, Some(&words)));
    assert_eq!(
        read(&["generated", "--corpus", "c", "--check", "--command"]),
        generated(true, Some(&[]))
    );
    // After `=` the option takes that one word, as it was given last, and the words set
    // aside stand in its place.
    assert_eq!(
        read(&["generated", "--corpus", "c", "--command=a", "--command=b"]),
        generated(false, Some(&["b"]))
    );
    assert_eq!(
        read(&["generated", "--corpus", "c", "--command="]),
        generated(false, Some(&[""]))
    );
    assert_eq!(
        read(&[
            "generated",
            "--corpus",
            "c",
            "--command=a",
            "--command",
            "b"
        ]),
        generated(false, Some(&["b"]))
    );
    // Behind a `--` before the subcommand's name nothing is set aside: the option takes
    // the arguments up to a `--`, options among them.
    assert_eq!(
        read(&[
            "--",
            "generated",
            "--corpus",
            "c",
            "--command",
            "a",
            "--check"
        ]),
        generated(false, Some(&["a", "--check"]))
    );
    assert_eq!(
        read(&["--", "generated", "--corpus", "c", "--command"]),
        generated(false, Some(&[]))
    );
    assert_eq!(
        read(&["generated", "-h", "--command", "x"]),
        Parsed::Help(GENERATED.help)
    );
    assert_eq!(
        read(&["generated", "--corpus", "c", "-hx"]),
        Parsed::Help(GENERATED.help)
    );
    let top = "vep-golden-corpus: error: unrecognized arguments:";
    let own = "vep-golden-corpus generated: error:";
    let required = "the following arguments are required: --corpus";
    let cases: Vec<(&[&str], String)> = vec![
        (&["generated"], format!("{own} {required}")),
        (
            &["generated", "--corpus"],
            format!("{own} argument --corpus: expected one argument"),
        ),
        (
            &["generated", "--corpus", "--command", "x"],
            format!("{own} argument --corpus: expected one argument"),
        ),
        // A prefix of an option names none.
        (
            &["generated", "--corpus", "c", "--comm", "a"],
            format!("{top} --comm a"),
        ),
        (
            &["generated", "--corpus=c", "--che"],
            format!("{top} --che"),
        ),
        (&["generated", "--cor", "c"], format!("{own} {required}")),
        (&["generated", "--corpus", "c", "-x"], format!("{top} -x")),
        (
            &["generated", "--corpus", "c", "--check", "--command=x", "y"],
            format!("{top} y"),
        ),
        (
            &["generated", "--corpus", "c", "--check=1"],
            format!("{own} argument --check: ignored explicit argument '1'"),
        ),
        (
            &["generated", "--corpus", "c", "extra"],
            format!("{top} extra"),
        ),
        (
            &[
                "generated",
                "--corpus",
                "c",
                "--check",
                "extra",
                "--command",
                "x",
            ],
            format!("{top} extra"),
        ),
        // The words set aside leave the arguments before them to be read alone.
        (
            &["generated", "--command", "a", "--corpus", "c"],
            format!("{own} {required}"),
        ),
        (
            &["generated", "--command", "x", "-h"],
            format!("{own} {required}"),
        ),
        (
            &["generated", "--corpus", "c", "--", "--command", "a"],
            format!("{top} --"),
        ),
        (
            &["generated", "--corpus", "c", "--", "extra"],
            format!("{top} -- extra"),
        ),
        (
            &[
                "--",
                "generated",
                "--corpus",
                "c",
                "--command",
                "a",
                "--",
                "b",
            ],
            format!("{top} -- b"),
        ),
        (
            &["--", "generated", "--command"],
            format!("{own} {required}"),
        ),
        (
            &["--", "generated", "--command", "--corpus", "c"],
            format!("{own} {required}"),
        ),
    ];
    for (args, line) in cases {
        assert_eq!(refused(args), line, "{args:?}");
    }
}

#[test]
fn the_option_that_drops_reference_statements_is_a_flag_of_select_and_flag_set() {
    assert_eq!(
        refused(&["select", "--drop-reference-statements=1", "--release", "1"]),
        "vep-golden-corpus select: error: argument --drop-reference-statements: ignored explicit argument '1'"
    );
    assert_eq!(
        refused(&["flag-set", "c", "--drop-reference-statements=x"]),
        "vep-golden-corpus flag-set: error: argument --drop-reference-statements: ignored explicit argument 'x'"
    );
    assert_eq!(
        refused(&["select", "--dr"]),
        "vep-golden-corpus select: error: the following arguments are required: --release, --assembly, --dataset, --out"
    );
    let Parsed::Run(command) = read(&["flag-set", "c", "--drop-ref", "--list"]) else {
        panic!("a prefix of the option is read");
    };
    let Command::FlagSet(args) = *command else {
        panic!("a flag-set run");
    };
    assert!(args.drop_reference_statements && args.list);
}
