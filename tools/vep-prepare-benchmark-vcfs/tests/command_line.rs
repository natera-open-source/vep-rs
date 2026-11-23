// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary's exit code and bytes for help, refused command lines, a value that names no
//! contig, a thread bound that is no count, and runs that fail: an input directory without
//! inputs, an output that cannot be made, a gzip input that is cut off.

mod common;

use common::{
    assert_clean_exit, fixture, names_in, prepare, run, run_with, stderr, stdout, Scratch, USAGE,
};
use std::fs;
use std::path::Path;

fn refusal(message: &str) -> String {
    format!("{USAGE}vep-prepare-benchmark-vcfs: error: {message}\n")
}

/// Runs `options` after the two directories and returns stderr, holding the exit code to
/// 2, stdout to nothing and the output directory to absent.
fn refused(scratch: &Scratch, options: &[&str]) -> String {
    let output = prepare(&fixture("contigs"), &scratch.out(), options, None);
    assert_eq!(output.status.code(), Some(2), "{options:?}");
    assert!(output.stdout.is_empty(), "{options:?}");
    assert!(!scratch.out().exists(), "{options:?}");
    stderr(&output)
}

#[test]
fn help_is_the_usage_and_the_options() {
    for args in [
        &["--help"][..],
        &["-h"],
        &["--he"],
        &["-h", "--max-variants", "x"],
    ] {
        let output = run(args);
        assert_clean_exit(&output, "help");
        assert_eq!(
            stdout(&output),
            fs::read_to_string(fixture("expected_help.txt")).unwrap()
        );
    }
}

#[test]
fn a_refused_command_line_exits_2_with_the_usage_and_one_error_line() {
    let scratch = Scratch::new("refused");
    let output = run::<&str>(&[]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        stderr(&output),
        refusal("the following arguments are required: --input-dir, --output-dir")
    );
    for (options, message) in [
        (&["stray"][..], "unrecognized arguments: stray"),
        (&["--frobnicate"], "unrecognized arguments: --frobnicate"),
        (&["--", "--no-canonical-contigs"], "unrecognized arguments: -- --no-canonical-contigs"),
        (
            &["--c", "x"],
            "ambiguous option: --c could match --chrom-mode, --cache-dir, --canonical-contigs",
        ),
        (
            &["--chrom-mode", "both"],
            "argument --chrom-mode: invalid choice: 'both' (choose from 'auto', 'keep', 'strip', 'add')",
        ),
        (
            &["--max-variants", "many"],
            "argument --max-variants: invalid int value: 'many'",
        ),
        (
            &["--max-variants", "1.5"],
            "argument --max-variants: invalid int value: '1.5'",
        ),
        (&["--cache-dir"], "argument --cache-dir: expected one argument"),
        (
            &["--no-canonical-contigs=x"],
            "argument --no-canonical-contigs: ignored explicit argument 'x'",
        ),
        (&["--help=x"], "argument -h/--help: ignored explicit argument 'x'"),
    ] {
        assert_eq!(refused(&scratch, options), refusal(message), "{options:?}");
    }
    // An integer past the digits `int()` converts is no integer.
    let long = "1".repeat(4301);
    assert_eq!(
        refused(&scratch, &["--max-variants", &long]),
        refusal(&format!(
            "argument --max-variants: invalid int value: '{long}'"
        ))
    );
}

#[test]
fn abbreviated_options_run_as_the_whole_ones_do() {
    let scratch = Scratch::new("abbreviated");
    let input = fixture("contigs");
    let whole = prepare(
        &input,
        &scratch.out(),
        &[
            "--chrom-mode",
            "add",
            "--max-variants",
            "3",
            "--no-canonical-contigs",
        ],
        None,
    );
    assert_clean_exit(&whole, "whole options");
    let whole_files = names_in(&scratch.out());
    let whole_text = fs::read(scratch.out().join("ucsc.vcf")).unwrap();
    fs::remove_dir_all(scratch.out()).unwrap();
    let abbreviated = run(&[
        "--i",
        input.to_str().unwrap(),
        &format!("--o={}", scratch.out().display()),
        "--ch=add",
        "--m",
        "3",
        "--no",
    ]);
    assert_clean_exit(&abbreviated, "abbreviated options");
    assert_eq!(stdout(&abbreviated), stdout(&whole));
    assert_eq!(names_in(&scratch.out()), whole_files);
    assert_eq!(
        fs::read(scratch.out().join("ucsc.vcf")).unwrap(),
        whole_text
    );
}

#[test]
fn a_negative_cap_in_the_digits_of_another_script_is_a_value() {
    let scratch = Scratch::new("other-digits");
    let output = prepare(
        &fixture("contigs"),
        &scratch.out(),
        &["--max-variants", "-\u{663}"],
        None,
    );
    assert_clean_exit(&output, "a cap in the digits of another script");
    let manifest = fs::read_to_string(scratch.out().join("manifest.json")).unwrap();
    assert!(
        manifest.contains("\n  \"max_variants\": -3,\n"),
        "{manifest}"
    );
}

#[test]
fn a_value_that_names_no_contig_exits_1_with_one_line_before_anything_is_made() {
    let scratch = Scratch::new("no-contig");
    let listing = fixture("no_contigs.txt");
    let listing = listing.to_str().unwrap();
    for (value, shown) in [
        (",", "','".to_string()),
        ("", "''".to_string()),
        (" ,\t, ", "' ,\\t, '".to_string()),
        (listing, format!("'{listing}'")),
    ] {
        let output = prepare(
            Path::new("/nonexistent/inputs"),
            &scratch.out(),
            &["--canonical-contigs", value],
            None,
        );
        assert_eq!(output.status.code(), Some(1), "{value:?}");
        assert_eq!(
            stderr(&output),
            format!("error: --canonical-contigs {shown} names no contig\n")
        );
        assert!(output.stdout.is_empty());
        assert!(!scratch.out().exists());
    }
}

#[test]
fn a_contig_list_that_is_not_utf8_fails_the_run() {
    let scratch = Scratch::new("contig-bytes");
    let list = scratch.at("contigs.txt");
    fs::write(&list, b"1\n\xff\n").unwrap();
    let output = prepare(
        &fixture("contigs"),
        &scratch.out(),
        &["--canonical-contigs", list.to_str().unwrap()],
        None,
    );
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        stderr(&output),
        format!(
            "ERROR: [vep-prepare-benchmark-vcfs] {}: the contig list is not UTF-8 text\n",
            list.display()
        )
    );
    assert!(!scratch.out().exists());
}

#[test]
fn a_directory_without_inputs_fails_after_the_output_directory_is_made() {
    let scratch = Scratch::new("no-inputs");
    let empty = scratch.at("empty");
    fs::create_dir(&empty).unwrap();
    fs::write(scratch.at("a_file"), b"").unwrap();
    for (input, options, expected_stdout) in [
        (empty.clone(), &[][..], ""),
        (scratch.at("absent"), &[], ""),
        (scratch.at("a_file"), &[], ""),
        (
            empty.clone(),
            &["--assembly", "GRCm39"],
            "WARNING: --canonical-contigs requested but --assembly is GRCm39, which names no \
             canonical set; the filter is OFF for this run. Name the set with --canonical-contigs \
             <file or list> (a cache's info.json contig_order) or pass --no-canonical-contigs to \
             silence this warning.\n",
        ),
    ] {
        let _ = fs::remove_dir_all(scratch.out());
        let output = prepare(&input, &scratch.out(), options, None);
        assert_eq!(output.status.code(), Some(1), "{input:?}");
        assert_eq!(
            stderr(&output),
            format!(
                "ERROR: [vep-prepare-benchmark-vcfs] No .vcf or .vcf.gz files found in {}\n",
                input.display()
            )
        );
        assert_eq!(stdout(&output), expected_stdout);
        assert!(names_in(&scratch.out()).is_empty(), "{input:?}");
    }
}

#[test]
fn an_output_directory_that_cannot_be_made_fails_the_run() {
    let scratch = Scratch::new("no-output");
    fs::write(scratch.at("a_file"), b"").unwrap();
    for output_dir in [scratch.at("a_file"), scratch.at("a_file/below")] {
        let output = prepare(&fixture("contigs"), &output_dir, &[], None);
        assert_eq!(output.status.code(), Some(1));
        assert!(
            stderr(&output).starts_with(&format!(
                "ERROR: [vep-prepare-benchmark-vcfs] {}: ",
                output_dir.display()
            )),
            "{}",
            stderr(&output)
        );
        assert_eq!(stderr(&output).lines().count(), 1);
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn a_failed_input_ends_the_run_after_the_inputs_before_it_and_without_a_manifest() {
    let scratch = Scratch::new("cut-off");
    let input = scratch.at("in");
    fs::create_dir(&input).unwrap();
    let good = fs::read(fixture("contigs/ensembl.vcf.gz")).unwrap();
    fs::write(input.join("a_whole.vcf.gz"), &good).unwrap();
    fs::write(input.join("b_cut_off.vcf.gz"), &good[..good.len() - 5]).unwrap();
    fs::write(input.join("c_whole.vcf.gz"), &good).unwrap();
    fs::write(input.join("0_plain.vcf"), b"#h\n1\t1\n").unwrap();
    for threads in ["1", "3"] {
        let _ = fs::remove_dir_all(scratch.out());
        let output = prepare(&input, &scratch.out(), &[], Some(threads));
        assert_eq!(output.status.code(), Some(1));
        assert_eq!(
            stdout(&output),
            "prepared 0_plain.vcf -> 0_plain.vcf (mode=keep, variants=1, headers=1, canonical_dropped=0)\n\
             prepared a_whole.vcf.gz -> a_whole.vcf (mode=keep, variants=8, headers=3, canonical_dropped=7)\n"
        );
        assert_eq!(
            stderr(&output),
            format!(
                "ERROR: [vep-prepare-benchmark-vcfs] {}: compressed file ended before the \
                 end-of-stream marker was reached\n",
                input.join("b_cut_off.vcf.gz").display()
            )
        );
        let written = names_in(&scratch.out());
        assert!(written.contains(&"a_whole.vcf".to_string()));
        assert!(!written.contains(&"manifest.json".to_string()));
        assert_eq!(
            fs::read(scratch.out().join("a_whole.vcf")).unwrap(),
            fs::read(fixture("expected/contigs_default/ensembl.vcf")).unwrap()
        );
    }
    // What is not gzip under a gzip name fails as it is opened, before its output exists.
    fs::remove_file(input.join("b_cut_off.vcf.gz")).unwrap();
    fs::write(input.join("b_plain_text.vcf.gz"), b"#h\n1\t1\n").unwrap();
    let _ = fs::remove_dir_all(scratch.out());
    let output = prepare(&input, &scratch.out(), &[], Some("1"));
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        stderr(&output),
        format!(
            "ERROR: [vep-prepare-benchmark-vcfs] {}: corrupt gzip stream: not a gzip member\n",
            input.join("b_plain_text.vcf.gz").display()
        )
    );
    assert_eq!(names_in(&scratch.out()), ["0_plain.vcf", "a_whole.vcf"]);
}

#[test]
fn a_thread_bound_that_is_no_count_fails_the_run() {
    let scratch = Scratch::new("thread-bound");
    for value in ["0", "four", "-1", "1.5"] {
        let output = run_with(
            &[
                "--input-dir",
                fixture("contigs").to_str().unwrap(),
                "--output-dir",
                scratch.out().to_str().unwrap(),
            ],
            Some(value),
        );
        assert_eq!(output.status.code(), Some(1), "{value}");
        assert_eq!(
            stderr(&output),
            format!(
                "ERROR: [vep-prepare-benchmark-vcfs] VEP_PREPARE_BENCHMARK_VCFS_THREADS={value}: \
                 expected a whole number of threads, 1 or more\n"
            )
        );
        assert!(!scratch.out().exists());
    }
    // An empty value is the host's cores.
    let output = prepare(&fixture("contigs"), &scratch.out(), &[], Some(""));
    assert_clean_exit(&output, "an empty thread bound");
}
