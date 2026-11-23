// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The consequence sections of the `--stats_text` report of a run whose
//! `--custom` source is keyed `Consequence`, on line 9 of the custom corpus's
//! input, for whose site the corpus's `custom/sites.vcf.gz` holds the record
//! `exact_2` on one of the two alleles. The report counts the consequence
//! terms of the run's rows and nothing else, in every output format. Ensembl
//! VEP 116.3 counts every item of the row hash's `Consequence` list
//! (`Stats.pm:282`), onto which a source of that key has pushed its records'
//! names (`OutputFactory.pm:1301`), so its report of a default, tab or VCF run
//! lists the name as a consequence type
//! (`docs/divergences/stats-consequences-custom-key-named-consequence.md`).

use std::collections::BTreeSet;
use std::path::Path;
use std::process::Stdio;

use vep_core::consequence::Consequence;
use vep_core::terms::Vocabulary;

mod common;

/// The 1-based line, in the custom corpus's `variants.vcf`, of the record
/// whose two alleles lie in no transcript's reach, and the columns that open
/// it.
const RECORD_LINE: usize = 9;
const RECORD: &str = "21\t5380280\trs1444031311\t";

/// The line of a record of the same input whose first alternate, which the
/// source's record `exact_3` names, has nine transcript rows.
const NINE_ROW_LINE: usize = 16;
const NINE_ROW_RECORD: &str = "21\t10542415\tmulti_eq_ref_second\t";

/// The committed GRCh38-custom corpus.
fn corpus() -> common::Corpus {
    common::corpora()
        .into_iter()
        .find(|corpus| corpus.name == "GRCh38-custom")
        .expect("the GRCh38-custom corpus is committed")
}

/// The corpus's header and its one record at `line`, which opens with
/// `record`, written under `dir`.
fn record_file(dir: &Path, line: usize, record: &str) -> std::path::PathBuf {
    let text = std::fs::read_to_string(corpus().dir.join("variants.vcf")).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let found = lines[line - 1];
    assert!(
        found.starts_with(record),
        "line {line} of the corpus's input is another record"
    );
    let mut kept: Vec<&str> = lines
        .iter()
        .copied()
        .filter(|line| line.starts_with('#'))
        .collect();
    kept.push(found);
    let input = dir.join(format!("record_{line}.vcf"));
    std::fs::write(&input, format!("{}\n", kept.join("\n"))).unwrap();
    input
}

/// Runs `vep` on the record at [`RECORD_LINE`] with the corpus's site file as
/// an `exact` source keyed `short_name` and returns the output and the text
/// report. `format` is the output format's flag, empty for the default format.
fn annotate(dir: &Path, short_name: &str, format: &str) -> (String, String) {
    annotate_line(dir, RECORD_LINE, RECORD, short_name, format)
}

/// [`annotate`] on the record at `line` of the corpus's input.
fn annotate_line(
    dir: &Path,
    line: usize,
    record: &str,
    short_name: &str,
    format: &str,
) -> (String, String) {
    let corpus = corpus();
    let out = dir.join(format!("{line}.{short_name}.{format}.out"));
    let source = format!(
        "file={},short_name={short_name},format=vcf,type=exact,fields=AF%CLASS",
        corpus.dir.join("custom/sites.vcf.gz").display()
    );
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(record_file(dir, line, record))
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .args(corpus.cache_args())
        .args(["--assembly", &corpus.assembly])
        .args(["--force_overwrite", "--stats_text", "--quiet"])
        .args(["--custom", &source]);
    if !format.is_empty() {
        cmd.arg(format!("--{format}"));
    }
    let output = cmd.stdin(Stdio::null()).output().unwrap();
    assert!(
        output.status.success(),
        "{short_name} {format:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report = std::fs::read_to_string(format!("{}_summary.txt", out.display())).unwrap();
    (std::fs::read_to_string(&out).unwrap(), report)
}

/// The rows of the report's section `title`, in the report's order.
fn section(report: &str, title: &str) -> Vec<String> {
    let header = format!("[{title}]");
    assert!(
        report.lines().any(|line| line == header),
        "the report has no section {header}"
    );
    report
        .lines()
        .skip_while(|line| *line != header)
        .skip(1)
        .take_while(|line| !line.is_empty())
        .map(str::to_string)
        .collect()
}

/// The record's two alleles lie in no transcript's reach, so the run writes
/// two rows whose one term is `intergenic_variant`, and the report counts
/// that term twice over the rows and once as the record's most severe. The
/// name of the source's record is no consequence type and adds no line, in
/// any of the four formats. Ensembl VEP 116.3 writes `exact_2\t1` before the
/// term's line in the default, tab and VCF formats and the term's line alone
/// in JSON, whose factory takes the custom annotations off the variant before
/// a row is built (`OutputFactory/JSON.pm:250-251`).
#[test]
fn a_custom_key_named_consequence_adds_no_consequence_type_to_the_report() {
    let dir = tempfile::tempdir().unwrap();
    for format in ["", "tab", "vcf", "json"] {
        let (_, report) = annotate(dir.path(), "Consequence", format);
        assert_eq!(
            section(&report, "Consequences (most severe)"),
            ["intergenic_variant\t1"],
            "{format:?}"
        );
        assert_eq!(
            section(&report, "Consequences (all)"),
            ["intergenic_variant\t2"],
            "{format:?}"
        );
        assert!(
            section(&report, "Coding consequences").is_empty(),
            "{format:?}"
        );
    }
}

/// The tab row of the allele the source's record names lists the record's
/// name after the row's term in the `Consequence` column, as Ensembl VEP's
/// row does: the column holds the key's values, and the report counts the
/// row's terms.
#[test]
fn the_consequence_column_lists_the_records_name_the_report_does_not_count() {
    let dir = tempfile::tempdir().unwrap();
    let (output, report) = annotate(dir.path(), "Consequence", "tab");
    let cells: Vec<&str> = output
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').nth(6).unwrap())
        .collect();
    assert_eq!(cells, ["intergenic_variant,exact_2", "intergenic_variant"]);
    assert_eq!(
        section(&report, "Consequences (all)"),
        ["intergenic_variant\t2"]
    );
}

/// The key's name moves no count of the report: every chart section and the
/// `[General statistics]` section of the run under the key `Consequence`
/// equal those of the same source under the key `sites`, in each format.
#[test]
fn the_report_counts_what_the_same_source_gives_under_another_key() {
    let dir = tempfile::tempdir().unwrap();
    for format in ["", "tab", "vcf", "json"] {
        let (_, keyed) = annotate(dir.path(), "Consequence", format);
        let (_, plain) = annotate(dir.path(), "sites", format);
        for title in [
            "General statistics",
            "Variant classes",
            "Consequences (most severe)",
            "Consequences (all)",
            "Coding consequences",
            "Variants by chromosome",
        ] {
            assert_eq!(
                section(&keyed, title),
                section(&plain, title),
                "[{title}] under {format:?}"
            );
        }
    }
}

/// On an allele of several rows the name is on each row's cell and in no
/// count: the first alternate of the record at line 16 has nine transcript
/// rows, each of which ends its `Consequence` cell with `exact_3`, and the
/// report counts their terms, seven `missense_variant`, two `intron_variant`
/// and two `NMD_transcript_variant`. Ensembl VEP 116.3 writes `exact_3\t9`
/// before them, one count for each row (`OutputFactory.pm:401`).
#[test]
fn a_name_on_nine_rows_of_one_allele_is_counted_on_none() {
    let dir = tempfile::tempdir().unwrap();
    let (output, report) = annotate_line(
        dir.path(),
        NINE_ROW_LINE,
        NINE_ROW_RECORD,
        "Consequence",
        "tab",
    );
    let named = output
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').nth(6).unwrap())
        .filter(|cell| cell.ends_with(",exact_3"))
        .count();
    assert_eq!(named, 9);
    assert_eq!(
        section(&report, "Consequences (all)"),
        [
            "missense_variant\t7",
            "intron_variant\t2",
            "NMD_transcript_variant\t2"
        ]
    );
    assert_eq!(
        section(&report, "Consequences (most severe)"),
        ["missense_variant\t1"]
    );
    assert_eq!(
        section(&report, "Coding consequences"),
        ["missense_variant\t7"]
    );
}

/// The comparison of two statistics reports sets a key of
/// `[Consequences (all)]` aside under this entry's class only where it names
/// no consequence type, and reads the names of the consequence types from its
/// own table, `tools/vep-stats-compare/data/consequence_terms.tsv`: one row
/// for each consequence, its SO term, then its display and NCBI labels, a
/// dash where the vocabulary has none. The table holds every name the report
/// can count a consequence under and no other.
#[test]
fn the_statistics_comparison_knows_every_name_of_a_consequence_type() {
    let path = format!(
        "{}/../../tools/vep-stats-compare/data/consequence_terms.tsv",
        env!("CARGO_MANIFEST_DIR")
    );
    let table = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut rows = table.lines();
    assert_eq!(rows.next(), Some("SO_term\tdisplay_term\tNCBI_term"));
    let rows: Vec<Vec<&str>> = rows.map(|row| row.split('\t').collect()).collect();
    let label = |consequence: Consequence, vocabulary: Vocabulary| {
        let term = vocabulary.term(consequence);
        if term == consequence.so_term() {
            "-"
        } else {
            term
        }
    };
    let expected: Vec<Vec<&str>> = Consequence::ALL
        .iter()
        .map(|&consequence| {
            vec![
                consequence.so_term(),
                label(consequence, Vocabulary::Display),
                label(consequence, Vocabulary::Ncbi),
            ]
        })
        .collect();
    let sorted = |rows: &[Vec<&str>]| -> BTreeSet<String> {
        rows.iter().map(|row| row.join("\t")).collect()
    };
    assert_eq!(rows.len(), Consequence::ALL.len(), "{path}");
    assert_eq!(sorted(&rows), sorted(&expected), "{path}");
}
