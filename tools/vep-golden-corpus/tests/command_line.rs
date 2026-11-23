// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary with requests for help and with refused command lines, and holds the
//! exit code and the bytes of each.

mod common;

use common::{clean, refused, run, Run};
use std::ffi::OsStr;
use std::os::unix::ffi::OsStrExt;

const TOP_USAGE: &str = "\
usage: vep-golden-corpus [-h]
                         {select,classify,check,identifiers,covers,flag-set,generated} ...
";
const GENERATED_USAGE: &str = "\
usage: vep-golden-corpus generated [-h] --corpus CORPUS [--check]
                                   [--command ...]
";
const IDENTIFIERS_USAGE: &str =
    "usage: vep-golden-corpus identifiers [-h] [--check] path [path ...]\n";
const COVERS_USAGE: &str = "\
usage: vep-golden-corpus covers [-h] [--add ADD] [--list]
                                [--flag-set FLAG_SET] [--matrix MATRIX]
                                corpus
";
const CHECK_USAGE: &str = "\
usage: vep-golden-corpus check [-h] --corpus CORPUS --dataset DATASET
                               [--matrix MATRIX]
                               [--family-columns FAMILY_COLUMNS]
                               [--require-genotype-class CLASS]
";

fn help(name: &str) -> String {
    let path = format!("{}/src/help/{name}.txt", common::ROOT);
    std::fs::read_to_string(path).unwrap()
}

fn usage_error(usage: &str, line: &str) -> Run {
    Run {
        code: 2,
        stdout: String::new(),
        stderr: format!("{usage}{line}\n"),
    }
}

#[test]
fn help_goes_to_standard_output_in_the_layout_of_80_columns() {
    assert_eq!(run(&["--help"]), clean(0, &help("top")));
    assert_eq!(run(&["-h"]), clean(0, &help("top")));
    let subcommands = [
        "select",
        "classify",
        "check",
        "identifiers",
        "covers",
        "flag-set",
        "generated",
    ];
    for name in subcommands {
        let text = help(name);
        assert_eq!(run(&[name, "--help"]), clean(0, &text), "{name}");
        assert_eq!(run(&[name, "-h"]), clean(0, &text), "{name}");
    }
    let top = help("top");
    assert!(top.starts_with(&format!("{TOP_USAGE}\nBuild a golden corpus:")));
    assert!(top.ends_with(
        "\npositional arguments:\n  {select,classify,check,identifiers,covers,flag-set,generated}\n\noptions:\n  -h, --help            show this help message and exit\n"
    ));
    assert!(help("generated").starts_with(&format!(
        "{GENERATED_USAGE}\nWrite or check the generated block of the corpus manifest: the command that\n"
    )));
    assert!(help("identifiers").starts_with(&format!(
        "{IDENTIFIERS_USAGE}\nReplace every identifier under restricted redistribution terms with the\n"
    )));
    assert!(help("covers").starts_with(&format!(
        "{COVERS_USAGE}\npositional arguments:\n  corpus\n"
    )));
    assert!(help("check").starts_with(&format!("{CHECK_USAGE}\noptions:\n")));
    // The width of a terminal changes no byte.
    let wide = std::process::Command::new(common::BIN)
        .args(["covers", "-h"])
        .env("COLUMNS", "200")
        .output()
        .unwrap();
    assert_eq!(String::from_utf8(wide.stdout).unwrap(), help("covers"));
}

#[test]
fn a_refused_command_line_ends_with_the_usage_and_exit_code_2() {
    assert_eq!(
        run::<&str>(&[]),
        usage_error(
            TOP_USAGE,
            "vep-golden-corpus: error: the following arguments are required: cmd"
        )
    );
    assert_eq!(
        run(&["build"]),
        usage_error(
            TOP_USAGE,
            "vep-golden-corpus: error: argument cmd: invalid choice: 'build' (choose from 'select', 'classify', 'check', 'identifiers', 'covers', 'flag-set', 'generated')"
        )
    );
    assert_eq!(
        run(&["generated", "--check"]),
        usage_error(
            GENERATED_USAGE,
            "vep-golden-corpus generated: error: the following arguments are required: --corpus"
        )
    );
    assert_eq!(
        run(&["generated", "--corpus", "corpus", "--comm", "program"]),
        usage_error(
            TOP_USAGE,
            "vep-golden-corpus: error: unrecognized arguments: --comm program"
        )
    );
    assert_eq!(
        run(&["identifiers", "--check"]),
        usage_error(
            IDENTIFIERS_USAGE,
            "vep-golden-corpus identifiers: error: the following arguments are required: path"
        )
    );
    assert_eq!(
        run(&["identifiers", "corpus", "--check", "other"]),
        usage_error(
            TOP_USAGE,
            "vep-golden-corpus: error: unrecognized arguments: other"
        )
    );
    assert_eq!(
        run(&["covers", "corpus", "--nope"]),
        usage_error(
            TOP_USAGE,
            "vep-golden-corpus: error: unrecognized arguments: --nope"
        )
    );
    assert_eq!(
        run(&["covers"]),
        usage_error(
            COVERS_USAGE,
            "vep-golden-corpus covers: error: the following arguments are required: corpus"
        )
    );
    assert_eq!(
        run(&["covers", "corpus", "--add"]),
        usage_error(
            COVERS_USAGE,
            "vep-golden-corpus covers: error: argument --add: expected one argument"
        )
    );
    assert_eq!(
        run(&["check", "--dataset", "s01=ref.txt=in.vcf"]),
        usage_error(
            CHECK_USAGE,
            "vep-golden-corpus check: error: the following arguments are required: --corpus"
        )
    );
    assert_eq!(
        run(&["check", "--corpus", "corpus", "--dataset"]),
        usage_error(
            CHECK_USAGE,
            "vep-golden-corpus check: error: argument --dataset: expected one argument"
        )
    );
}

#[test]
fn an_argument_outside_utf8_is_refused_by_its_place() {
    let arguments = [
        OsStr::new("covers"),
        OsStr::from_bytes(b"corp\xffus"),
        OsStr::new("--list"),
    ];
    assert_eq!(
        run(&arguments),
        refused(1, "ERROR: [vep-golden-corpus] argument 2 is not UTF-8")
    );
}

#[test]
fn a_thread_bound_that_is_no_count_is_refused_before_a_file_is_read() {
    let scratch = common::scratch("thread-bound");
    let corpus = common::two_record_corpus(&scratch.path);
    let reference = common::fixture("two_record/ref.txt");
    let reference = reference.to_str().unwrap();
    let input = common::fixture("two_record/in.vcf");
    let dataset = format!("pool={reference}={}", input.to_str().unwrap());
    let commands: [&[&str]; 4] = [
        &[
            "select",
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--dataset",
            &dataset,
            "--out",
            &corpus,
        ],
        &["check", "--corpus", &corpus, "--dataset", &dataset],
        &[
            "classify",
            "--corpus",
            &corpus,
            "--vep-default",
            reference,
            "--vep-rs-output",
            reference,
        ],
        &["identifiers", "--check", &corpus],
    ];
    let bounded = |command: &[&str], bound: &str| {
        let output = std::process::Command::new(common::BIN)
            .args(command)
            .env("VEP_GOLDEN_CORPUS_THREADS", bound)
            .output()
            .unwrap();
        Run {
            code: output.status.code().unwrap(),
            stdout: String::from_utf8(output.stdout).unwrap(),
            stderr: String::from_utf8(output.stderr).unwrap(),
        }
    };
    let written = common::manifest_text(&corpus);
    for bound in ["0", "-1", "two", "1.5"] {
        let line = format!(
            "ERROR: [vep-golden-corpus] VEP_GOLDEN_CORPUS_THREADS={bound}: expected a whole number of threads, 1 or more"
        );
        for command in commands {
            assert_eq!(
                bounded(command, bound),
                refused(1, &line),
                "{} under {bound}",
                command[0]
            );
        }
    }
    assert_eq!(
        common::manifest_text(&corpus),
        written,
        "a refused run writes nothing"
    );
    // A bound of one or more reads on, and so does an empty one.
    let read = format!(
        "{corpus}: 2 files read, 0 hold an identifier in another form than the fabricated one\n"
    );
    for bound in ["1", " 3 ", ""] {
        assert_eq!(
            bounded(commands[1], bound),
            clean(0, "observed 2, covered by corpus 2, missing 0\n"),
            "{bound:?}"
        );
        assert_eq!(bounded(commands[3], bound), clean(0, &read), "{bound:?}");
    }
    // `check` reads no file for identifiers under a bound that is no count: a corpus that
    // holds an identifier, made up here from its letters and its digits, is refused for the
    // bound, and under a bound of one for the identifier.
    let held = format!("x|{}{}|y\n", "COSV", "12");
    common::write(std::path::Path::new(&corpus), "notes.txt", &held);
    assert_eq!(
        bounded(commands[1], "1"),
        common::failed(
            "notes.txt holds 1 identifier under restricted redistribution terms in another \
             form than the fabricated one; the identifiers subcommand replaces them"
        )
    );
    assert_eq!(
        bounded(commands[1], "two"),
        refused(
            1,
            "ERROR: [vep-golden-corpus] VEP_GOLDEN_CORPUS_THREADS=two: expected a whole number of threads, 1 or more"
        )
    );
}
