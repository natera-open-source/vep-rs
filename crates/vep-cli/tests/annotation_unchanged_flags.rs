// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Flags whose runs keep the annotation rows of a run without them. On a cache
//! without regulatory data, `--regulatory` declares the regulatory columns,
//! fills none of them and writes one line to the warnings file, where Ensembl
//! VEP writes none; a cache with regulatory data annotates its features
//! (`regulatory.rs`). The statistics flags write the report beside the output
//! and change no byte of it.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

fn corpus_cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

/// Records of the GRCh38 release 116 corpus on both of its chromosomes.
const RECORDS: &str = "\
21\t30282169\tstop_ins4_twostops_rev\tC\tCTATT\t.\t.\t.
21\t30348409\tstop_del2_lost_rev\tCTA\tC\t.\t.\t.
21\t30430756\tstop_ins5_newstop_fwd\tT\tTCCTAA\t.\t.\t.
21\t34716478\tpenult_ins3_stop_fwd\tA\tAATA\t.\t.\t.
22\t20061573\tpartial_codon_del3_fwd\tGAAG\tG\t.\t.\t.
";

/// The columns `--regulatory` declares, the `regulatory` entry of `@FLAG_FIELDS`
/// (`Constants.pm:136`): `BIOTYPE`, filled only under `--biotype`
/// (`OutputFactory.pm:1589`), then the motif columns.
const REGULATORY_COLUMNS: [&str; 6] = [
    "BIOTYPE",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
];

/// The columns that stay empty under `--regulatory` even when `--biotype`, which
/// `--everything` and `--vcf` switch on, fills `BIOTYPE`.
const MOTIF_COLUMNS: [&str; 5] = {
    let [_biotype, motif_columns @ ..] = REGULATORY_COLUMNS;
    motif_columns
};

const FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

fn write_input(dir: &Path) -> PathBuf {
    let path = dir.join("in.vcf");
    fs::write(
        &path,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{RECORDS}"),
    )
    .unwrap();
    path
}

/// Runs the binary on `in.vcf` under `dir`, writing `<dir>/<label>.<format>`,
/// with `cache` naming the cache and `extra` the flags under test. Nothing is
/// asserted about the exit status.
fn run(dir: &Path, label: &str, format: &str, cache: &[&str], extra: &[&str]) -> Output {
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(dir.join("in.vcf"))
        .arg("-o")
        .arg(dir.join(format!("{label}.{format}")))
        .args(cache)
        .args(["--assembly", "GRCh38", "--force_overwrite"])
        .args(extra)
        .stdin(Stdio::null());
    if format != "default" {
        cmd.arg(format!("--{format}"));
    }
    cmd.output().unwrap()
}

/// [`run`] on the corpus cache, asserting success and at least one data line
/// per record of [`RECORDS`]; the output text and the standard error.
fn annotate(dir: &Path, label: &str, format: &str, extra: &[&str]) -> (String, String) {
    let cache = corpus_cache();
    let output = run(
        dir,
        label,
        format,
        &["--json_cache", cache.to_str().unwrap()],
        extra,
    );
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{extra:?} ({format}): {stderr}");
    let text = fs::read_to_string(dir.join(format!("{label}.{format}"))).unwrap();
    assert!(
        data_lines(&text).len() >= RECORDS.lines().count(),
        "{extra:?} ({format}): the run annotates"
    );
    (text, stderr)
}

/// The output without the header lines that differ between any two runs: the
/// command line, which names the flags under test, the time of the run and the
/// cache path.
fn without_volatile_headers(text: &str) -> Vec<&str> {
    const VOLATILE: [&str; 5] = [
        "## VEP command-line:",
        "## Output produced at ",
        "## Using cache in ",
        "##VEP-command-line=",
        "##VEP=",
    ];
    text.lines()
        .filter(|line| !VOLATILE.iter().any(|prefix| line.starts_with(prefix)))
        .collect()
}

fn data_lines(text: &str) -> Vec<&str> {
    text.lines().filter(|line| !line.starts_with('#')).collect()
}

/// A cache without regulatory data warns and writes empty motif columns.
#[test]
fn missing_regulatory_data_warns_and_writes_empty_motif_columns() {
    const NOTICE: &str = "--regulatory: this cache holds no regulatory data; no RegulatoryFeature or MotifFeature row is written";
    const EVERYTHING_NOTICE: &str = NOTICE;
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());

    let (text, stderr) = annotate(dir.path(), "everything", "tab", &["--everything"]);
    assert_eq!(
        stderr.matches(EVERYTHING_NOTICE).count(),
        1,
        "everything: one warning on standard error:\n{stderr}"
    );
    let warnings = fs::read_to_string(dir.path().join("everything.tab_warnings.txt")).unwrap();
    assert!(
        warnings.contains(EVERYTHING_NOTICE),
        "everything: {warnings}"
    );
    let header: Vec<&str> = text
        .lines()
        .find(|line| line.starts_with('#') && !line.starts_with("##"))
        .expect("a tab header line")
        .split('\t')
        .collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|cell| *cell == name)
            .unwrap_or_else(|| panic!("everything: the tab header declares {name}"))
    };
    let feature_type = column("Feature_type");
    let motif_indexes = MOTIF_COLUMNS.map(column);
    for row in data_lines(&text) {
        let cells: Vec<&str> = row.split('\t').collect();
        assert!(
            !["RegulatoryFeature", "MotifFeature"].contains(&cells[feature_type]),
            "everything: no regulatory row: {row}"
        );
        for &index in &motif_indexes {
            assert_eq!(
                cells[index], "-",
                "everything: {} is empty on {row}",
                header[index]
            );
        }
    }

    let descriptions: Vec<String> = REGULATORY_COLUMNS
        .iter()
        .map(|column| format!("## {column} : "))
        .collect();
    for format in FORMATS {
        let (without, _) = annotate(dir.path(), "plain", format, &[]);
        let (with, stderr) = annotate(dir.path(), "regulatory", format, &["--regulatory"]);
        assert_eq!(
            stderr.matches(NOTICE).count(),
            1,
            "{format}: one warning on standard error:\n{stderr}"
        );
        let warnings = dir.path().join(format!("regulatory.{format}_warnings.txt"));
        let warnings = fs::read_to_string(&warnings).unwrap();
        assert_eq!(
            warnings.lines().count(),
            1,
            "{format}: one warnings-file line:\n{warnings}"
        );
        assert!(warnings.contains(NOTICE), "{format}: {warnings}");
        assert!(
            !data_lines(&with)
                .iter()
                .any(|line| line.contains("RegulatoryFeature") || line.contains("MotifFeature")),
            "{format}: no regulatory row"
        );
        let without = without_volatile_headers(&without);
        let with = without_volatile_headers(&with);
        match format {
            "json" => assert_eq!(with, without, "{format}: the records are unchanged"),
            "default" | "tab" => {
                for description in &descriptions {
                    assert_eq!(
                        with.iter()
                            .filter(|line| line.starts_with(description.as_str()))
                            .count(),
                        1,
                        "{format}: {description} declared once"
                    );
                }
                let other_lines: Vec<&str> = with
                    .iter()
                    .copied()
                    .filter(|line| !descriptions.iter().any(|d| line.starts_with(d.as_str())))
                    .collect();
                if format == "default" {
                    // The Extra column omits an empty field, so the descriptions
                    // are the whole difference.
                    assert_eq!(other_lines, without, "{format}: nothing else changes");
                } else {
                    let column_headers = REGULATORY_COLUMNS
                        .iter()
                        .map(|column| format!("\t{column}"))
                        .collect::<String>();
                    let empty_cells = "\t-".repeat(REGULATORY_COLUMNS.len());
                    let stripped: Vec<&str> = other_lines
                        .iter()
                        .map(|line| {
                            if line.starts_with("#Uploaded_variation") {
                                line.strip_suffix(column_headers.as_str())
                                    .unwrap_or_else(|| {
                                        panic!("{format}: six column names on {line}")
                                    })
                            } else if line.starts_with('#') {
                                line
                            } else {
                                line.strip_suffix(empty_cells.as_str()).unwrap_or_else(|| {
                                    panic!("{format}: six empty cells on {line}")
                                })
                            }
                        })
                        .collect();
                    assert_eq!(stripped, without, "{format}: nothing else changes");
                }
            }
            _ => {
                let fields = MOTIF_COLUMNS
                    .iter()
                    .map(|column| format!("|{column}"))
                    .collect::<String>();
                let stripped: Vec<String> = with
                    .iter()
                    .map(|line| {
                        if line.starts_with("##INFO=<ID=CSQ,") {
                            let head = line
                                .strip_suffix(&format!("{fields}\">"))
                                .unwrap_or_else(|| panic!("vcf: five field names on {line}"));
                            format!("{head}\">")
                        } else if line.starts_with('#') {
                            line.to_string()
                        } else {
                            let (record, csq) = line.split_once("CSQ=").unwrap();
                            let entries: Vec<&str> = csq
                                .split(',')
                                .map(|entry| {
                                    entry.strip_suffix("|||||").unwrap_or_else(|| {
                                        panic!("vcf: five empty fields on {entry}")
                                    })
                                })
                                .collect();
                            format!("{record}CSQ={}", entries.join(","))
                        }
                    })
                    .collect();
                assert_eq!(stripped, without, "{format}: nothing else changes");
            }
        }
    }
}

/// The statistics flags change no byte of the annotation output: the report is a
/// file of its own beside it (`BaseRunner.pm:254`), skipped under `--no_stats`
/// (`Runner.pm:144-148`), named by `--stats_file` with the report's extension
/// (`:258-271`), written as text under `--stats_text` (`:305-309`) and as HTML
/// under `--stats_html` or neither flag (`:311-315`, `Runner.pm:787-790`).
#[test]
fn statistics_flags_leave_the_annotation_output_unchanged() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let named = dir.path().join("stats.html");
    let named = named.to_str().unwrap();
    let bare = dir.path().join("report");
    let bare = bare.to_str().unwrap();
    let (reference, _) = annotate(dir.path(), "plain", "default", &["--no_stats"]);
    let runs: [(&str, &[&str], &[&str]); 6] = [
        ("no_stats", &["--no_stats", "--stats_text"], &[]),
        (
            "stats_text",
            &["--stats_text"],
            &["stats_text.default_summary.txt"],
        ),
        (
            "stats_file",
            &["--stats_text", "--stats_file", named],
            &["stats.txt"],
        ),
        ("sf", &["--stats_text", "--sf", bare], &["report.txt"]),
        ("neither", &[], &["neither.default_summary.html"]),
        (
            "both",
            &["--stats_text", "--stats_html"],
            &["both.default_summary.html", "both.default_summary.txt"],
        ),
    ];
    for (label, flags, _) in &runs {
        let (text, stderr) = annotate(dir.path(), label, "default", flags);
        assert_eq!(
            without_volatile_headers(&text),
            without_volatile_headers(&reference),
            "{flags:?} changes the output"
        );
        assert!(!stderr.contains("WARNING"), "{flags:?}: {stderr}");
    }
    let written: BTreeSet<String> = fs::read_dir(dir.path())
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    let expected: BTreeSet<String> = ["home", "in.vcf", "plain.default"]
        .into_iter()
        .map(str::to_string)
        .chain(runs.iter().map(|(label, _, _)| format!("{label}.default")))
        .chain(
            runs.iter()
                .flat_map(|(_, _, reports)| reports.iter().map(|r| r.to_string())),
        )
        .collect();
    assert_eq!(written, expected, "one report per format each run asks for");
    for report in [
        "stats_text.default_summary.txt",
        "stats.txt",
        "report.txt",
        "both.default_summary.txt",
    ] {
        let text = fs::read_to_string(dir.path().join(report)).unwrap();
        assert!(
            text.starts_with("[VEP run statistics]\n"),
            "{report}: {text}"
        );
    }
    for report in ["neither.default_summary.html", "both.default_summary.html"] {
        let html = fs::read_to_string(dir.path().join(report)).unwrap();
        assert!(html.starts_with("<html>\n"), "{report}: {html}");
        assert!(!html.contains("<script"), "{report}");
    }
}
