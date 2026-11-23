// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary's exit code and bytes for what the goldens' table cannot name: shortened
//! options, refused command lines, help, the report's path, the two environment variables,
//! names outside UTF-8, links, and the inputs the tool refuses where Python reads them.

mod common;

use common::{fixtures, in_json, run, run_args, text, Bounds, Scratch, BUFFER_ENV, THREADS_ENV};
use std::ffi::{OsStr, OsString};
use std::fs;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::os::unix::fs::symlink;
use std::path::{Path, PathBuf};
use std::process::Command;

const USAGE: &str = "usage: vep-cache-diff [-h] --baseline-dir BASELINE_DIR --builder-dir BUILDER_DIR --output-report OUTPUT_REPORT\n";
const CACHE: &str = "caches/identical/cache";

fn golden(name: &str, report: &Path) -> (String, String) {
    let expected = fixtures().join("expected");
    let fixtures = fixtures();
    let stdout = fs::read_to_string(expected.join(format!("{name}.stdout")))
        .unwrap()
        .replace("{{FIXTURES}}", fixtures.to_str().unwrap())
        .replace("{{REPORT}}", report.to_str().unwrap());
    let json = fs::read_to_string(expected.join(format!("{name}.json")))
        .unwrap()
        .replace("{{FIXTURES}}", &in_json(&fixtures));
    (stdout, json)
}

#[test]
fn options_may_be_shortened_and_take_their_value_behind_an_equals_sign() {
    let scratch = Scratch::new("shortened");
    let report = scratch.path("report.json");
    let (stdout, json) = golden("identical", &report);
    let out = report.to_str().unwrap();
    for args in [
        vec!["--ba", CACHE, "--bu", CACHE, "--o", out],
        vec![
            &format!("--output-report={out}"),
            &format!("--builder-dir={CACHE}"),
            &format!("--baseline={CACHE}"),
        ],
        vec![
            "--baseline-dir",
            "caches/nowhere",
            "--builder-d",
            CACHE,
            "--out",
            out,
            "--baseline-dir",
            CACHE,
        ],
    ] {
        let _ = fs::remove_file(&report);
        let output = run_args(&args, Bounds::DEFAULT);
        assert_eq!(output.status.code(), Some(0), "{args:?}");
        assert_eq!(text(&output.stdout), stdout, "{args:?}");
        assert!(output.stderr.is_empty(), "{args:?}");
        assert_eq!(fs::read_to_string(&report).unwrap(), json, "{args:?}");
    }
}

#[test]
fn a_refused_command_line_exits_2_with_the_usage_and_one_error_line() {
    let cases: [(&[&str], &str); 9] = [
        (
            &[],
            "the following arguments are required: --baseline-dir, --builder-dir, --output-report",
        ),
        (
            &["--baseline-dir", CACHE, "--builder-dir", CACHE],
            "the following arguments are required: --output-report",
        ),
        (
            &["--b", CACHE, "--builder-dir", CACHE, "--output-report", "o"],
            "ambiguous option: --b could match --baseline-dir, --builder-dir",
        ),
        (
            &[
                "--baseline-dir",
                CACHE,
                "--builder-dir",
                CACHE,
                "--output-report",
                "o",
                "extra",
                "more",
            ],
            "unrecognized arguments: extra more",
        ),
        (
            &[
                "--baseline-dir",
                CACHE,
                "--builder-dir",
                CACHE,
                "--output-report",
                "o",
                "--",
            ],
            "unrecognized arguments: --",
        ),
        (
            &["--baseline-dir", "--builder-dir", CACHE],
            "argument --baseline-dir: expected one argument",
        ),
        (
            &["--baseline-dir", CACHE, "--builder-dir", CACHE, "--output-report"],
            "argument --output-report: expected one argument",
        ),
        (
            &["--help=x"],
            "argument -h/--help: ignored explicit argument 'x'",
        ),
        (
            &["--=x"],
            "ambiguous option: --=x could match --help, --baseline-dir, --builder-dir, --output-report",
        ),
    ];
    for (args, message) in cases {
        let output = run_args(args, Bounds::DEFAULT);
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert!(output.stdout.is_empty(), "{args:?}");
        assert_eq!(
            text(&output.stderr),
            format!("{USAGE}vep-cache-diff: error: {message}\n"),
            "{args:?}"
        );
    }
    assert!(!fixtures().join("o").exists());
}

#[test]
fn help_is_written_on_stdout_with_exit_0() {
    for args in [
        &["-h"][..],
        &["--help"],
        &["--he"],
        &["-hx"],
        &["--baseline-dir", "nowhere", "--help", "--bogus"],
    ] {
        let output = run_args(args, Bounds::DEFAULT);
        assert_eq!(output.status.code(), Some(0), "{args:?}");
        let stdout = text(&output.stdout);
        assert!(stdout.starts_with(USAGE), "{args:?}");
        for part in [
            "--baseline-dir BASELINE_DIR",
            "--builder-dir BUILDER_DIR",
            "--output-report OUTPUT_REPORT",
            THREADS_ENV,
            BUFFER_ENV,
        ] {
            assert!(stdout.contains(part), "{part}");
        }
        assert!(output.stderr.is_empty(), "{args:?}");
    }
}

#[test]
fn the_report_goes_where_python_resolves_its_path() {
    let scratch = Scratch::new("report");
    let cache = fixtures().join(CACHE);
    // Directories above the report are made; a path is resolved before it is printed.
    let asked = scratch.path("made/./here/../here/report.json");
    let resolved = scratch.path("made/here/report.json");
    let output = run(&cache, &cache, &asked, Bounds::DEFAULT);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let (stdout, json) = golden("identical", &resolved);
    assert_eq!(text(&output.stdout), stdout);
    assert_eq!(fs::read_to_string(&resolved).unwrap(), json);
    // A report that exists is replaced; one behind a link is written where the link leads.
    fs::write(
        &resolved,
        "an earlier report, longer than the next one\n".repeat(200),
    )
    .unwrap();
    symlink(&resolved, scratch.path("link.json")).unwrap();
    let output = run(&cache, &cache, &scratch.path("link.json"), Bounds::DEFAULT);
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(text(&output.stdout), stdout);
    assert_eq!(fs::read_to_string(&resolved).unwrap(), json);
    // A relative path is taken from the working directory.
    let output = Command::new(common::BIN)
        .current_dir(&scratch.0)
        .args(["--baseline-dir".as_ref(), cache.as_os_str()])
        .args(["--builder-dir".as_ref(), cache.as_os_str()])
        .args(["--output-report", "relative/report.json"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0));
    let relative = scratch.path("relative/report.json");
    assert_eq!(text(&output.stdout), golden("identical", &relative).0);
    assert_eq!(fs::read_to_string(&relative).unwrap(), json);
}

#[test]
fn a_report_that_cannot_be_written_fails_after_the_last_comparison() {
    let scratch = Scratch::new("unwritable");
    let cache = fixtures().join(CACHE);
    scratch.write("file", "not a directory\n");
    for asked in [
        scratch.path("file/under/report.json"),
        scratch.path("file/report.json"),
        scratch.0.clone(),
    ] {
        let output = run(&cache, &cache, &asked, Bounds::DEFAULT);
        assert_eq!(output.status.code(), Some(1), "{}", asked.display());
        let (stdout, _) = golden("identical", &asked);
        let through_the_comparisons = stdout
            .split_once("\nReport written to: ")
            .map(|(before, _)| before.to_string())
            .unwrap();
        assert_eq!(text(&output.stdout), through_the_comparisons);
        let stderr = text(&output.stderr);
        assert!(
            stderr.starts_with(&format!("ERROR: [vep-cache-diff] {}: ", asked.display()))
                && stderr.lines().count() == 1,
            "{stderr}"
        );
    }
    assert_eq!(
        fs::read_to_string(scratch.path("file")).unwrap(),
        "not a directory\n"
    );
}

#[test]
fn the_environment_variables_take_a_whole_number_of_1_or_more() {
    let scratch = Scratch::new("environment");
    let report = scratch.path("report.json");
    let args: Vec<OsString> = [
        "--baseline-dir",
        CACHE,
        "--builder-dir",
        CACHE,
        "--output-report",
    ]
    .iter()
    .map(OsString::from)
    .chain([report.clone().into_os_string()])
    .collect();
    let run_with = |name: &str, value: &OsStr| {
        Command::new(common::BIN)
            .current_dir(fixtures())
            .args(&args)
            .env_remove(THREADS_ENV)
            .env_remove(BUFFER_ENV)
            .env(name, value)
            .output()
            .unwrap()
    };
    let (stdout, json) = golden("identical", &report);
    for name in [THREADS_ENV, BUFFER_ENV] {
        for value in ["", "1", " 2 ", "64", "1000000"] {
            let _ = fs::remove_file(&report);
            let output = run_with(name, value.as_ref());
            assert_eq!(output.status.code(), Some(0), "{name}={value:?}");
            assert_eq!(text(&output.stdout), stdout, "{name}={value:?}");
            assert_eq!(
                fs::read_to_string(&report).unwrap(),
                json,
                "{name}={value:?}"
            );
        }
        let _ = fs::remove_file(&report);
        for value in [
            OsStr::new("0"),
            OsStr::new("-1"),
            OsStr::new("two"),
            OsStr::new("1.5"),
            OsStr::from_bytes(b"\xff"),
        ] {
            let output = run_with(name, value);
            assert_eq!(output.status.code(), Some(1), "{name}={value:?}");
            assert!(output.stdout.is_empty(), "{name}={value:?}");
            let what = if name == THREADS_ENV {
                "threads"
            } else {
                "bytes"
            };
            assert_eq!(
                text(&output.stderr),
                format!(
                    "ERROR: [vep-cache-diff] {name}={}: expected a whole number of {what}, 1 or more\n",
                    value.to_string_lossy()
                )
            );
            assert!(!report.exists());
        }
    }
}

fn outside_utf8(prefix: &str, suffix: &str) -> OsString {
    let mut bytes = prefix.as_bytes().to_vec();
    bytes.push(0xff);
    bytes.extend_from_slice(suffix.as_bytes());
    OsString::from_vec(bytes)
}

#[test]
fn names_outside_utf8_are_sorted_and_written_as_python_decodes_them() {
    let scratch = Scratch::new("bytes");
    if let Err(e) = fs::create_dir(scratch.0.join(outside_utf8("dir", ""))) {
        // A file system that takes UTF-8 names alone cannot hold these inputs.
        eprintln!("not run: {e}");
        return;
    }
    let record = r#"[{"stable_id": "T1"}]"#;
    let shard = |side: &str, name: &OsStr| {
        let directory = scratch.path(side).join("transcripts/21");
        fs::create_dir_all(&directory).unwrap();
        fs::write(directory.join(name), record).unwrap();
    };
    shard("a", &outside_utf8("r", ".json"));
    shard("a", "\u{e9}.json".as_ref());
    shard("a", "\u{e000}.json".as_ref());
    shard("a", "z.json".as_ref());
    fs::create_dir_all(scratch.path("b")).unwrap();
    let report = scratch.path("report.json");
    let output = run(
        &scratch.path("a"),
        &scratch.path("b"),
        &report,
        Bounds::DEFAULT,
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let written = fs::read_to_string(&report).unwrap();
    let regions: Vec<&str> = written
        .lines()
        .filter(|line| line.contains("\"region\""))
        .collect();
    // The stray byte stands as U+DCFF: behind `z`, before U+E000, though 0xff is the
    // largest byte.
    assert_eq!(
        regions,
        [
            r#"        "region": "21/r\udcff.json","#,
            r#"        "region": "21/z.json","#,
            r#"        "region": "21/\u00e9.json","#,
            r#"        "region": "21/\ue000.json","#,
        ]
    );
    // A cache directory and a report outside UTF-8: the bytes on stdout, the lone
    // surrogate in the report.
    let cache = scratch.0.join(outside_utf8("dir", ""));
    fs::create_dir_all(cache.join("transcripts")).unwrap();
    let report = scratch.0.join(outside_utf8("out", ".json"));
    let output = run(&cache, &scratch.path("b"), &report, Bounds::DEFAULT);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let line = |label: &str, path: &Path| {
        let mut line = label.as_bytes().to_vec();
        line.extend_from_slice(path.as_os_str().as_bytes());
        line.push(b'\n');
        line
    };
    let stdout = output.stdout;
    let baseline_line = line("  Baseline: ", &cache);
    assert!(stdout
        .windows(baseline_line.len())
        .any(|window| window == baseline_line));
    assert!(stdout.ends_with(&line("\nReport written to: ", &report)));
    let written = fs::read_to_string(&report).unwrap();
    let in_report = format!(
        "  \"baseline_dir\": \"{}/dir\\udcff\",",
        in_json(&scratch.0)
    );
    assert_eq!(written.lines().nth(1), Some(in_report.as_str()));
    // A directory that is not there is named on stderr with the byte as an escape.
    let nowhere = scratch.0.join(outside_utf8("nodir", ""));
    let output = run(&nowhere, &scratch.path("b"), &report, Bounds::DEFAULT);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        text(&output.stderr),
        format!(
            "ERROR: [diff_caches] baseline directory not found: {}/nodir\\udcff\n",
            scratch.0.display()
        )
    );
}

#[test]
fn links_are_followed_as_python_follows_them() {
    let scratch = Scratch::new("links");
    let record = |id: &str| format!(r#"[{{"stable_id": "{id}"}}]"#);
    scratch.write("real/transcripts/21/1.json", record("T1"));
    scratch.write("real/transcripts/22/1.json", record("T2"));
    scratch.write(
        "real/info.json",
        r#"{"variation_cols": ["variation_name"]}"#,
    );
    scratch.write_gzipped("real/variations/21/all_vars.gz", "rs1\nrs2\n");
    // A cache whose directories and files are links: `transcripts` and a shard are read
    // through their links, a linked directory under `transcripts` is not entered, a
    // linked chromosome under `variations` is, and a link that leads nowhere is an empty
    // shard and an empty variation file.
    let linked = scratch.path("linked");
    fs::create_dir_all(linked.join("variations/9")).unwrap();
    fs::create_dir_all(linked.join("inner/21")).unwrap();
    symlink(scratch.path("linked/inner"), linked.join("transcripts")).unwrap();
    symlink(
        scratch.path("real/transcripts/21/1.json"),
        linked.join("inner/21/1.json"),
    )
    .unwrap();
    symlink(scratch.path("real/transcripts/22"), linked.join("inner/22")).unwrap();
    symlink("nowhere", linked.join("inner/21/2.json")).unwrap();
    symlink(
        scratch.path("real/variations/21"),
        linked.join("variations/21"),
    )
    .unwrap();
    symlink("nowhere", linked.join("variations/9/all_vars.gz")).unwrap();
    symlink(scratch.path("real/info.json"), linked.join("info.json")).unwrap();
    symlink(&linked, scratch.path("cache")).unwrap();
    let report = scratch.path("report.json");
    let output = run(
        &scratch.path("real"),
        &scratch.path("cache"),
        &report,
        Bounds::DEFAULT,
    );
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let stdout = text(&output.stdout);
    let expected = format!(
        "Comparing caches:\n  Baseline: {root}/real\n  Builder:  {root}/linked\n  Comparing transcripts...\n    1 shared, 1 matched, 0 differed, rate=1.0000\n  Comparing variations...\n    2 shared, 2 matched, 0 differed, rate=1.0000\n  Comparing info.json...\n    matches=True\n\nReport written to: {root}/report.json\n",
        root = scratch.0.display()
    );
    assert_eq!(stdout, expected);
    let written = fs::read_to_string(&report).unwrap();
    assert!(written.contains("\"region\": \"22/1.json\""), "{written}");
    assert!(!written.contains("\"region\": \"21/"), "{written}");
}

#[test]
fn inputs_the_tool_refuses_where_python_reads_them_end_the_run() {
    let scratch = Scratch::new("bounds");
    let report = scratch.path("report.json");
    let cases = [
        (
            "id_number",
            "baseline/transcripts/1/1.json: a record's id is neither a string nor a value Python counts false",
            "  Comparing transcripts...\n",
        ),
        (
            "column_number",
            "baseline/info.json: variation_cols holds an item that is no string",
            "  Comparing variations...\n",
        ),
        (
            "depth",
            "baseline/transcripts/1/1.json: a kept value nested past 1000 levels at byte 1031",
            "  Comparing transcripts...\n",
        ),
    ];
    for (name, detail, last_line) in cases {
        let directory: PathBuf = fixtures().join("caches/bounds").join(name);
        let output = run(
            &directory.join("baseline"),
            &directory.join("builder"),
            &report,
            Bounds::DEFAULT,
        );
        assert_eq!(output.status.code(), Some(1), "{name}");
        assert_eq!(
            text(&output.stderr),
            format!("ERROR: [vep-cache-diff] {}/{detail}\n", directory.display()),
            "{name}"
        );
        assert!(text(&output.stdout).ends_with(last_line), "{name}");
        assert!(!report.exists(), "{name}");
    }
}
