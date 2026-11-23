// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;

fn read(args: &[&str]) -> Result<Command, String> {
    let argv: Vec<String> = args.iter().map(|a| a.to_string()).collect();
    parse(&argv)
}

fn run(args: &[&str]) -> Args {
    match read(args) {
        Ok(Command::Run(args)) => *args,
        other => panic!("{args:?}: {other:?}"),
    }
}

#[test]
fn every_option_is_read_with_its_default() {
    assert_eq!(run(&[]), Args::default());
    let defaults = Args::default();
    assert_eq!((defaults.jobs, defaults.timeout), (1, 1800));
    assert_eq!(defaults.reference, "116.3");
    let all = run(&[
        "--image",
        "img@sha256:0",
        "--ensembl-vep",
        "/usr/bin/env vep",
        "--vep",
        "./target//release/vep/",
        "--perl-cache-dir",
        "/caches/a",
        "--perl-cache-dir=b/",
        "--plugins-dir",
        "plugins",
        "--tool-dir",
        "tools/one",
        "--tool-dir",
        "tools/two",
        "--extra-dir",
        "loftee=/x//y/",
        "--entry",
        "one",
        "--entry",
        "two",
        "--jobs",
        "6",
        "--work-dir",
        "work/",
        "--keep",
        "--timeout",
        "30",
        "--expect-all-present",
        "--json",
        "out.json",
        "--check",
        "--reference",
        "116.2",
        "--page",
        "docs/page.md",
        "--reproductions",
        "records",
        "--list",
    ]);
    let text = |s: &str| Some(s.to_string());
    let list = |items: &[&str]| -> Vec<String> { items.iter().map(|i| i.to_string()).collect() };
    assert_eq!(
        all,
        Args {
            image: text("img@sha256:0"),
            ensembl_vep: text("/usr/bin/env vep"),
            vep: text("target/release/vep"),
            vep_names_a_directory: true,
            perl_cache_dir: list(&["/caches/a", "b"]),
            plugins_dir: text("plugins"),
            tool_dirs: list(&["tools/one", "tools/two"]),
            extra_dir: list(&["loftee=/x//y/"]),
            entry: list(&["one", "two"]),
            jobs: 6,
            work_dir: text("work"),
            keep: true,
            timeout: 30,
            expect_all_present: true,
            json: text("out.json"),
            check: true,
            reference: "116.2".to_string(),
            page: text("docs/page.md"),
            reproductions: text("records"),
            list: true,
        }
    );
}

#[test]
fn a_binary_written_with_a_directory_is_told_from_a_bare_name() {
    for (written, kept, with_a_directory) in [
        ("vep", "vep", false),
        ("./vep", "vep", true),
        ("bin/vep", "bin/vep", true),
        ("/opt/vep", "/opt/vep", true),
        ("", ".", false),
    ] {
        let args = run(&["--vep", written]);
        assert_eq!(args.vep.as_deref(), Some(kept), "{written:?}");
        assert_eq!(args.vep_names_a_directory, with_a_directory, "{written:?}");
    }
    // A later value replaces an earlier one whole.
    let args = run(&["--vep=./vep", "--vep", "vep"]);
    assert_eq!(args.vep.as_deref(), Some("vep"));
    assert!(!args.vep_names_a_directory);
}

#[test]
fn an_option_is_named_by_a_prefix_that_fits_it_alone() {
    let args = run(&["--ch", "--ref=116.2", "--pa", "p.md", "--rep", "r", "--l"]);
    assert!(args.check && args.list);
    assert_eq!(args.reference, "116.2");
    assert_eq!(args.page.as_deref(), Some("p.md"));
    assert_eq!(args.reproductions.as_deref(), Some("r"));
    let args = run(&["--jobs=-3", "--timeout", " 1_0 ", "--vep", ""]);
    assert_eq!((args.jobs, args.timeout), (-3, 10));
    assert_eq!(args.vep.as_deref(), Some("."));
    assert_eq!(read(&["-h"]), Ok(Command::Help));
    assert_eq!(read(&["--help", "--bogus"]), Ok(Command::Help));
    assert_eq!(read(&["--he"]), Ok(Command::Help));
}

#[test]
fn a_command_line_is_refused_in_argparse_words() {
    let refused = |args: &[&str]| read(args).expect_err("a refusal");
    assert_eq!(refused(&["--vep"]), "argument --vep: expected one argument");
    assert_eq!(
        refused(&["--check", "--page"]),
        "argument --page: expected one argument"
    );
    assert_eq!(
        refused(&["--vep", "--list"]),
        "argument --vep: expected one argument"
    );
    assert_eq!(
        refused(&["--jobs", "x"]),
        "argument --jobs: invalid int value: 'x'"
    );
    assert_eq!(
        refused(&["--timeout=1.5"]),
        "argument --timeout: invalid int value: '1.5'"
    );
    assert_eq!(
        refused(&["--keep=1"]),
        "argument --keep: ignored explicit argument '1'"
    );
    assert_eq!(
        refused(&["--list="]),
        "argument --list: ignored explicit argument ''"
    );
    assert_eq!(
        refused(&["--help=x"]),
        "argument -h/--help: ignored explicit argument 'x'"
    );
    assert_eq!(
        refused(&["--e", "x"]),
        "ambiguous option: --e could match --ensembl-vep, --extra-dir, --entry, --expect-all-present"
    );
    assert_eq!(
        refused(&["--re=1"]),
        "ambiguous option: --re=1 could match --reference, --reproductions"
    );
    assert_eq!(
        refused(&["--j", "2", "--l"]),
        "ambiguous option: --j could match --jobs, --json"
    );
    assert_eq!(refused(&["--bogus"]), "unrecognized arguments: --bogus");
    assert_eq!(
        refused(&["stray", "--list", "-x"]),
        "unrecognized arguments: stray -x"
    );
    assert_eq!(
        refused(&["--", "--list"]),
        "unrecognized arguments: -- --list"
    );
}
