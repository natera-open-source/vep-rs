// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A key that one `--custom` source declares twice, in the headers the binary
//! writes for one record: line 9 of the custom corpus's input. A source named
//! `S` that asks for the record field `count` and for the statistic `count`
//! declares the key `S_count` as a field and as a statistic; a field or a
//! statistic named twice declares its key twice as well. Ensembl VEP 116.3
//! tests a source's keys against the keys of the sources before it alone
//! (`OutputFactory.pm:2451`, read at `:2457`, `:2482` and `:2499`), so it
//! pushes a pair for each declaration (`:2493`, `:2503-2506`): two tab columns
//! whose description lines both carry the later pair's text
//! (`OutputFactory/Tab.pm:149`, `:153`, `:195`), two `CSQ` fields and two
//! `##INFO` lines of one ID (`OutputFactory/VCF.pm:276`, `:449`), and two key
//! lines of the default format (`OutputFactory/VEP_output.pm:167`). The cell
//! under the key holds the statistic where one is assigned and the field's
//! values where none is, so the tab and VCF headers list the key once with the
//! statistic's text and then the field's. The rows held here are the ones
//! Ensembl VEP writes for the record, its second listing of the key left out,
//! and the default format's key lines are Ensembl VEP's own. The key's one
//! listing and its text are the documented difference
//! `custom_key_named_as_a_plain_value_field`.

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The 1-based line of the record in the custom corpus's `variants.vcf`.
const RECORD_LINE: usize = 9;

/// The source file as a header prints it.
const FILE: &str = "[PATH]/statfield.vcf.gz";

/// The line of the source's own key in the tab and default headers.
const NAME_LINE: &str = "## S : [PATH]/statfield.vcf.gz";

/// The text of the key declared as the record field `count`.
const FIELD_TEXT: &str =
    "A made-up field named as a statistic. count field from [PATH]/statfield.vcf.gz";

/// The text of the key declared as the statistic `count`.
const STATISTIC_TEXT: &str = "count data from [PATH]/statfield.vcf.gz";

/// The text of the key's one line where it is declared as the field and as the
/// statistic: the statistic's text and then the field's, the order of their
/// values in a cell that holds both.
const BOTH_TEXTS: &str = "count data from [PATH]/statfield.vcf.gz; A made-up field named as a statistic. count field from [PATH]/statfield.vcf.gz";

/// The text of the key declared as the record field `CLASS`.
const CLASS_TEXT: &str = "The class of the corpus variant the record was written for: exact, overlap, multi_allelic or sv_span. CLASS field from [PATH]/statfield.vcf.gz";

/// The tab format's column line up to the columns of the source.
const TAB_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE";

/// The thirteen columns of the record's two rows, which no case changes.
const ROWS: [&str; 2] = [
    "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-",
    "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-",
];

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The record alone, written under `dir`: the header lines of the custom
/// corpus's input and its line [`RECORD_LINE`].
fn record_file(dir: &Path) -> PathBuf {
    let input = fs::read_to_string(corpus("GRCh38-custom").join("variants.vcf")).unwrap();
    let record = input
        .lines()
        .nth(RECORD_LINE - 1)
        .expect("the record's line");
    assert!(
        record.starts_with("21\t5380280\t"),
        "line {RECORD_LINE} of the custom corpus's input is the record"
    );
    let mut text: String = input
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    text.push_str(record);
    text.push('\n');
    let path = dir.join("record.vcf");
    fs::write(&path, text).unwrap();
    path
}

/// A site file written under `dir` whose one record stands at the record's
/// site with two INFO fields, `count`, the name of a statistic, and `CLASS`.
fn site_file(dir: &Path) -> PathBuf {
    common::indexed_vcf(
        dir,
        "statfield.vcf.gz",
        &[
            "##fileformat=VCFv4.2",
            "##INFO=<ID=CLASS,Number=1,Type=String,Description=\"The class of the corpus variant the record was written for: exact, overlap, multi_allelic or sv_span\">",
            "##INFO=<ID=count,Number=1,Type=Integer,Description=\"A made-up field named as a statistic\">",
            "##contig=<ID=21>",
        ],
        "21\t5380280\texact_2\tA\tC\t.\tPASS\tcount=7;CLASS=exact",
    )
}

/// The source `S` over `file`, whose records the variant overlaps, with the
/// `declared` options after it.
fn source(file: &Path, declared: &str) -> String {
    format!(
        "file={},short_name=S,format=vcf,type=overlap,{declared}",
        file.display()
    )
}

/// What a run's output states of the source's keys.
#[derive(Debug, PartialEq)]
struct Written {
    /// The header lines that describe a key of the source, in the output's
    /// order.
    key_lines: Vec<String>,
    /// The columns from the source's first on: the end of the tab format's
    /// column line or of the `CSQ` format; empty in the default format, whose
    /// column line names no key.
    columns: String,
    /// The rows, each from its fourteenth column on; a VCF record's `CSQ`
    /// entries from their third field on.
    cells: Vec<String>,
}

/// The run over the record in `format` (`default`, `tab` or `vcf`) under
/// `sources`, in their order; `label` names the run's output file.
fn written(dir: &Path, label: &str, format: &str, sources: &[String]) -> Written {
    let out = dir.join(format!("{label}.{format}.out"));
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(record_file(dir))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus("GRCh38-everything").join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"]);
    for source in sources {
        cmd.args(["--custom", source]);
    }
    if format != "default" {
        cmd.arg(format!("--{format}"));
    }
    let run = cmd.output().unwrap();
    assert!(
        run.status.success(),
        "vep {format} {label}:\n{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let text = fs::read_to_string(&out).unwrap();
    let key_lines = text
        .lines()
        .filter(|line| {
            ["## S ", "## S_", "##INFO=<ID=S,", "##INFO=<ID=S_"]
                .iter()
                .any(|opening| line.starts_with(opening))
        })
        .map(String::from)
        .collect();
    let columns = match format {
        "tab" => {
            let line = text
                .lines()
                .find(|line| line.starts_with("#Uploaded_variation"))
                .expect("the column line");
            line.strip_prefix(TAB_COLUMNS)
                .unwrap_or_else(|| panic!("{label}: the columns before the source's: {line}"))
                .to_string()
        }
        "vcf" => {
            let line = text
                .lines()
                .find(|line| line.starts_with("##INFO=<ID=CSQ,"))
                .expect("the CSQ line");
            let (_, after) = line
                .split_once("|HGNC_ID|SOURCE")
                .unwrap_or_else(|| panic!("{label}: the fields before the source's: {line}"));
            after.to_string()
        }
        _ => String::new(),
    };
    let rows: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
    let cells = match format {
        "vcf" => {
            assert_eq!(rows.len(), 1, "{label}: one VCF record");
            let (line, entries) = rows[0].split_once("CSQ=").expect("the CSQ entries");
            assert_eq!(line, "21\t5380280\trs1444031311\tA\tC,G\t.\t.\t", "{label}");
            entries
                .split(',')
                .zip(["C", "G"])
                .map(|(entry, allele)| {
                    let opening = format!("{allele}|intergenic_variant|MODIFIER{}", "|".repeat(21));
                    entry
                        .strip_prefix(&opening)
                        .unwrap_or_else(|| {
                            panic!("{label}: the fields before the source's: {entry}")
                        })
                        .to_string()
                })
                .collect()
        }
        _ => {
            assert_eq!(rows.len(), 2, "{label}: one row an allele");
            rows.iter()
                .zip(ROWS)
                .map(|(row, columns)| {
                    row.strip_prefix(columns)
                        .unwrap_or_else(|| panic!("{label}: the thirteen columns: {row}"))
                        .to_string()
                })
                .collect()
        }
    };
    Written {
        key_lines,
        columns,
        cells,
    }
}

/// A `## <key> : <text>` line of the tab or default header.
fn line(key: &str, text: &str) -> String {
    format!("## {key} : {text}")
}

/// The `##INFO` line the VCF header prints for a key of a source.
fn info(key: &str, text: &str) -> String {
    format!("##INFO=<ID={key},Number=.,Type=String,Description=\"{text}\">")
}

/// The expected statement of a run: its key lines, its columns and the cells
/// of both rows, which are equal on this record.
fn expected(key_lines: &[String], columns: &str, cells: &str) -> Written {
    Written {
        key_lines: key_lines.to_vec(),
        columns: columns.to_string(),
        cells: vec![cells.to_string(), cells.to_string()],
    }
}

/// A record field and a statistic of one name. The row's one value under the
/// key is the statistic, which Ensembl VEP assigns over the field's list
/// (`OutputFactory.pm:1303`, then `:1308`): `1`, the count of the records the
/// variant overlaps, where the field alone prints `7`. Ensembl VEP's tab
/// header gives both of its `S_count` lines the statistic's text, and its VCF
/// header a `##INFO` line with the field's text before the one with the
/// statistic's. The key is one tab column, one `CSQ` field and one `##INFO`
/// line, which carries both texts; the default format lists both declarations,
/// as Ensembl VEP does.
#[test]
fn a_key_declared_as_a_field_and_a_statistic_is_listed_once_with_both_texts() {
    let tmp = tempfile::tempdir().unwrap();
    let sources = [source(
        &site_file(tmp.path()),
        "fields=count,summary_stats=count",
    )];
    assert_eq!(
        written(tmp.path(), "field", "tab", &sources),
        expected(
            &[NAME_LINE.to_string(), line("S_count", BOTH_TEXTS)],
            "\tS\tS_count",
            "\tMODIFIER\t-\t-\t-\t-\texact_2\t1",
        )
    );
    assert_eq!(
        written(tmp.path(), "field", "vcf", &sources),
        expected(
            &[info("S", FILE), info("S_count", BOTH_TEXTS)],
            "|S|S_count\">",
            "|exact_2|1",
        )
    );
    assert_eq!(
        written(tmp.path(), "field", "default", &sources),
        expected(
            &[
                NAME_LINE.to_string(),
                line("S_count", FIELD_TEXT),
                line("S_count", STATISTIC_TEXT),
            ],
            "",
            "\tIMPACT=MODIFIER;S=exact_2;S_FILTER=PASS;S_count=1",
        )
    );
}

/// The same key with another field's key between its two declarations.
/// Ensembl VEP lists `S_count`, `S_CLASS`, `S_count`: the key keeps the place
/// of its first declaration.
#[test]
fn a_key_declared_twice_keeps_its_first_place_among_the_sources_keys() {
    let tmp = tempfile::tempdir().unwrap();
    let sources = [source(
        &site_file(tmp.path()),
        "fields=count%CLASS,summary_stats=count",
    )];
    assert_eq!(
        written(tmp.path(), "between", "tab", &sources),
        expected(
            &[
                NAME_LINE.to_string(),
                line("S_count", BOTH_TEXTS),
                line("S_CLASS", CLASS_TEXT),
            ],
            "\tS\tS_count\tS_CLASS",
            "\tMODIFIER\t-\t-\t-\t-\texact_2\t1\texact",
        )
    );
    assert_eq!(
        written(tmp.path(), "between", "vcf", &sources),
        expected(
            &[
                info("S", FILE),
                info("S_count", BOTH_TEXTS),
                info("S_CLASS", CLASS_TEXT),
            ],
            "|S|S_count|S_CLASS\">",
            "|exact_2|1|exact",
        )
    );
    assert_eq!(
        written(tmp.path(), "between", "default", &sources),
        expected(
            &[
                NAME_LINE.to_string(),
                line("S_count", FIELD_TEXT),
                line("S_CLASS", CLASS_TEXT),
                line("S_count", STATISTIC_TEXT),
            ],
            "",
            "\tIMPACT=MODIFIER;S=exact_2;S_CLASS=exact;S_FILTER=PASS;S_count=1",
        )
    );
}

/// A field named twice and a statistic named twice: both pairs Ensembl VEP
/// pushes carry one text, the key's own, with the file named once.
#[test]
fn a_field_or_a_statistic_named_twice_is_one_key_with_its_own_text() {
    let tmp = tempfile::tempdir().unwrap();
    let file = site_file(tmp.path());
    for (label, declared, text, cell, extra) in [
        (
            "fieldtwice",
            "fields=count%count",
            FIELD_TEXT,
            "7",
            "S_FILTER=PASS;S_count=7",
        ),
        (
            "stattwice",
            "summary_stats=count%count",
            STATISTIC_TEXT,
            "1",
            "S_count=1",
        ),
    ] {
        let sources = [source(&file, declared)];
        assert_eq!(
            written(tmp.path(), label, "tab", &sources),
            expected(
                &[NAME_LINE.to_string(), line("S_count", text)],
                "\tS\tS_count",
                &format!("\tMODIFIER\t-\t-\t-\t-\texact_2\t{cell}"),
            ),
            "{label}"
        );
        assert_eq!(
            written(tmp.path(), label, "vcf", &sources),
            expected(
                &[info("S", FILE), info("S_count", text)],
                "|S|S_count\">",
                &format!("|exact_2|{cell}"),
            ),
            "{label}"
        );
        assert_eq!(
            written(tmp.path(), label, "default", &sources),
            expected(
                &[
                    NAME_LINE.to_string(),
                    line("S_count", text),
                    line("S_count", text),
                ],
                "",
                &format!("\tIMPACT=MODIFIER;S=exact_2;{extra}"),
            ),
            "{label}"
        );
    }
}

/// A key that an earlier source listed takes the later source's `,<file>`
/// on the declaration that source repeats: the second source of the short
/// name `S` declares `S_count` as a field and as a statistic, so the file
/// joins the field's text and the statistic has a text of its own, the first
/// on the key's one line. Ensembl VEP appends the file to the field's text
/// once for each declaration (`OutputFactory.pm:2484`, `:2501`) and names the
/// statistic in no header, over a cell that holds it.
#[test]
fn a_key_an_earlier_source_listed_takes_the_file_for_each_later_declaration() {
    let tmp = tempfile::tempdir().unwrap();
    let file = site_file(tmp.path());
    let sources = [
        source(&file, "fields=count"),
        source(&file, "fields=count,summary_stats=count"),
    ];
    let name = format!("{FILE},{FILE}");
    let field = format!("{FIELD_TEXT},{FILE}");
    let key = format!("{STATISTIC_TEXT}; {field}");
    assert_eq!(
        written(tmp.path(), "twosources", "tab", &sources),
        expected(
            &[line("S", &name), line("S_count", &key)],
            "\tS\tS_count",
            "\tMODIFIER\t-\t-\t-\t-\texact_2,exact_2\t1",
        )
    );
    assert_eq!(
        written(tmp.path(), "twosources", "vcf", &sources),
        expected(
            &[info("S", &name), info("S_count", &key)],
            "|S|S_count\">",
            "|exact_2&exact_2|1",
        )
    );
    assert_eq!(
        written(tmp.path(), "twosources", "default", &sources),
        expected(
            &[
                line("S", &name),
                line("S_count", &field),
                line("S_count", STATISTIC_TEXT),
            ],
            "",
            "\tIMPACT=MODIFIER;S=exact_2,exact_2;S_FILTER=PASS,PASS;S_count=1",
        )
    );
}

/// A site file written under `dir` as `name` whose one record stands at the
/// record's site with an INFO field named `min`, the name of a statistic, and
/// with `qual` in its `QUAL` column, a VCF record's score.
fn site_file_with_a_min_field(dir: &Path, name: &str, qual: &str) -> PathBuf {
    common::indexed_vcf(
        dir,
        name,
        &[
            "##fileformat=VCFv4.2",
            "##INFO=<ID=min,Number=1,Type=Integer,Description=\"A made-up field named as a statistic\">",
            "##contig=<ID=21>",
        ],
        &format!("21\t5380280\texact_2\tA\tC\t{qual}\tPASS\tmin=7"),
    )
}

/// `min`, `max`, `mean` and `sum` are assigned only where a matched record
/// carries a score (`AnnotationSource/File.pm:424-434`, `:347-349`), a VCF
/// record's `QUAL`. Over a record with no score the cell under a key declared
/// as the field `min` and as the statistic `min` holds the field's value, `7`,
/// and over a record with the score 30 it holds the statistic, in both
/// engines; Ensembl VEP's tab header reads `min data from` on both of its
/// lines in either run. The key's one line carries both texts in both runs.
#[test]
fn a_key_holds_the_field_where_no_statistic_is_assigned_and_its_line_names_both() {
    let tmp = tempfile::tempdir().unwrap();
    for (name, qual, cell) in [("unscored.vcf.gz", ".", "7"), ("scored.vcf.gz", "30", "30")] {
        let file = format!("[PATH]/{name}");
        let field = format!("A made-up field named as a statistic. min field from {file}");
        let statistic = format!("min data from {file}");
        let both = format!("{statistic}; {field}");
        let sources = [source(
            &site_file_with_a_min_field(tmp.path(), name, qual),
            "fields=min,summary_stats=min",
        )];
        assert_eq!(
            written(tmp.path(), name, "tab", &sources),
            expected(
                &[line("S", &file), line("S_min", &both)],
                "\tS\tS_min",
                &format!("\tMODIFIER\t-\t-\t-\t-\texact_2\t{cell}"),
            ),
            "{name}"
        );
        assert_eq!(
            written(tmp.path(), name, "vcf", &sources),
            expected(
                &[info("S", &file), info("S_min", &both)],
                "|S|S_min\">",
                &format!("|exact_2|{cell}"),
            ),
            "{name}"
        );
        assert_eq!(
            written(tmp.path(), name, "default", &sources),
            expected(
                &[
                    line("S", &file),
                    line("S_min", &field),
                    line("S_min", &statistic),
                ],
                "",
                &format!("\tIMPACT=MODIFIER;S=exact_2;S_FILTER=PASS;S_min={cell}"),
            ),
            "{name}"
        );
    }
}
