// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A key that two `--custom` sources of one short name declare, in the headers
//! the binary writes for one record: line 9 of the custom corpus's input. Two
//! sources named `S` over two site files declare the key `S_count`, or
//! `S_min`, one as a record field and one as a statistic, or one of them as
//! both. Ensembl VEP 116.3 tests a later source's key against the keys the
//! earlier sources listed (`OutputFactory.pm:2451`, read at `:2482` and
//! `:2499`) and, where it finds the key, appends the later file to the pair
//! the key has, whatever that pair declares (`:2483-2484`, `:2500-2501`): a
//! field's text then stands alone over a cell that holds the later source's
//! statistic, and a statistic's text over a cell that holds a field's values.
//! Here a later source's file joins the pair of the declaration it repeats, a
//! declaration the key lacks is a pair of its own, and the key's one line in
//! the tab and VCF headers names the statistic first, the order of the values
//! in a cell that holds both. The cells held here are the ones Ensembl VEP
//! writes wherever it completes its run. The lines are the documented
//! difference `custom_repeated_short_name_header_lines`.

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The 1-based line of the record in the custom corpus's `variants.vcf`.
const RECORD_LINE: usize = 9;

/// The first source's file as a header prints it.
const FIRST: &str = "[PATH]/a.vcf.gz";

/// The second source's file as a header prints it.
const SECOND: &str = "[PATH]/b.vcf.gz";

/// The line of the shared short name: both files, in the sources' order.
const NAMES: &str = "[PATH]/a.vcf.gz,[PATH]/b.vcf.gz";

/// The tab format's column line up to the columns of the sources.
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

/// The two site files written under `dir`, each with one record at the
/// record's site that names the allele `C` and carries the INFO field
/// `field`, the name of a statistic: `a.vcf.gz` holds the record `exact_2`
/// with the value 7 and `b.vcf.gz` the record `second_1` with the value 9.
/// `quals` are the two records' `QUAL` cells, a VCF record's score.
fn site_files(dir: &Path, field: &str, quals: [&str; 2]) -> [PathBuf; 2] {
    let file = |name: &str, described: &str, record: &str, qual: &str, value: u8| {
        common::indexed_vcf(
            dir,
            name,
            &[
                "##fileformat=VCFv4.2",
                &format!("##INFO=<ID={field},Number=1,Type=Integer,Description=\"{described}\">"),
                "##contig=<ID=21>",
            ],
            &format!("21\t5380280\t{record}\tA\tC\t{qual}\tPASS\t{field}={value}"),
        )
    };
    [
        file(
            "a.vcf.gz",
            "The first file's made-up field",
            "exact_2",
            quals[0],
            7,
        ),
        file(
            "b.vcf.gz",
            "The second file's made-up field",
            "second_1",
            quals[1],
            9,
        ),
    ]
}

/// The text of the key `field` declared as a record field of the first file.
fn first_field(field: &str) -> String {
    format!("The first file's made-up field. {field} field from {FIRST}")
}

/// The text of the key `field` declared as a record field of the second file.
fn second_field(field: &str) -> String {
    format!("The second file's made-up field. {field} field from {SECOND}")
}

/// The text of the key declared as the statistic `stat` of the source over
/// `file`.
fn statistic(stat: &str, file: &str) -> String {
    format!("{stat} data from {file}")
}

/// Two sources named `S` of the type `kind` over `files`, the first with the
/// options `first` and the second with `second`.
fn sources(files: &[PathBuf; 2], kind: &str, first: &str, second: &str) -> [String; 2] {
    let source = |file: &Path, declared: &str| {
        format!(
            "file={},short_name=S,format=vcf,type={kind},{declared}",
            file.display()
        )
    };
    [source(&files[0], first), source(&files[1], second)]
}

/// What a run's output states of the sources' keys.
#[derive(Debug, PartialEq)]
struct Written {
    /// The header lines that describe a key of the sources, in the output's
    /// order.
    key_lines: Vec<String>,
    /// The columns from the sources' first on: the end of the tab format's
    /// column line or of the `CSQ` format; empty in the default format, whose
    /// column line names no key.
    columns: String,
    /// The rows, each from its fourteenth column on; a VCF record's `CSQ`
    /// entries from the sources' first field on.
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
                .unwrap_or_else(|| panic!("{label}: the columns before the sources': {line}"))
                .to_string()
        }
        "vcf" => {
            let line = text
                .lines()
                .find(|line| line.starts_with("##INFO=<ID=CSQ,"))
                .expect("the CSQ line");
            let (_, after) = line
                .split_once("|HGNC_ID|SOURCE")
                .unwrap_or_else(|| panic!("{label}: the fields before the sources': {line}"));
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
                            panic!("{label}: the fields before the sources': {entry}")
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

/// What the three formats state of the key `S_<field>` under two overlap
/// sources that declare it as `first` and `second` over files whose records
/// carry `quals`: one tab line and one `##INFO` line with `joined`, one
/// default line for each text of `pairs`, and under the key the one `cell` on
/// both rows, beside both records' names and a `PASS` for each of the
/// `filters` sources that name a field.
fn holds(
    label: &str,
    field: &str,
    quals: [&str; 2],
    (first, second): (&str, &str),
    joined: &str,
    pairs: &[String],
    (filters, cell): (usize, &str),
) {
    let tmp = tempfile::tempdir().unwrap();
    let files = site_files(tmp.path(), field, quals);
    let sources = sources(&files, "overlap", first, second);
    let key = format!("S_{field}");
    let both_rows = |cells: String| vec![cells.clone(), cells];
    let passes = vec!["PASS"; filters].join(",");
    assert_eq!(
        written(tmp.path(), label, "tab", &sources),
        Written {
            key_lines: vec![line("S", NAMES), line(&key, joined)],
            columns: format!("\tS\t{key}"),
            cells: both_rows(format!("\tMODIFIER\t-\t-\t-\t-\texact_2,second_1\t{cell}")),
        },
        "{label}: tab"
    );
    assert_eq!(
        written(tmp.path(), label, "vcf", &sources),
        Written {
            key_lines: vec![info("S", NAMES), info(&key, joined)],
            columns: format!("|S|{key}\">"),
            cells: both_rows(format!("|exact_2&second_1|{}", cell.replace(',', "&"))),
        },
        "{label}: VCF"
    );
    let mut default_lines = vec![line("S", NAMES)];
    default_lines.extend(pairs.iter().map(|text| line(&key, text)));
    assert_eq!(
        written(tmp.path(), label, "default", &sources),
        Written {
            key_lines: default_lines,
            columns: String::new(),
            cells: both_rows(format!(
                "\tIMPACT=MODIFIER;S=exact_2,second_1;S_FILTER={passes};{key}={cell}"
            )),
        },
        "{label}: default"
    );
}

/// An earlier source names the field `count` and a later source of the same
/// short name asks for the statistic `count`, alone or beside the field. The
/// cell holds the later source's statistic, `1`, which Ensembl VEP assigns
/// over the list of both sources' field values (`OutputFactory.pm:1303`, then
/// `:1308`); its one line for the key is the first source's field text with
/// the later file appended once for each declaration (`:2484`, `:2501`). The
/// statistic has a text of its own, with the file of the source that asks for
/// it, and the later file joins the field's text only where the later source
/// names the field.
#[test]
fn a_statistic_a_later_source_asks_for_is_described_beside_the_earlier_sources_field() {
    let field = first_field("count");
    let asked = statistic("count", SECOND);
    holds(
        "field-then-both",
        "count",
        [".", "."],
        ("fields=count", "fields=count,summary_stats=count"),
        &format!("{asked}; {field},{SECOND}"),
        &[format!("{field},{SECOND}"), asked.clone()],
        (2, "1"),
    );
    holds(
        "field-then-statistic",
        "count",
        [".", "."],
        ("fields=count", "summary_stats=count"),
        &format!("{asked}; {field}"),
        &[field, asked],
        (1, "1"),
    );
}

/// The declarations the other way round: an earlier source asks for the
/// statistic `count` and a later source names the field `count`, alone or
/// beside the statistic. The cell holds a source's statistic, and Ensembl
/// VEP's line is the statistic's text with the later file appended once for
/// each declaration (`OutputFactory.pm:2484`, `:2501`), though one of them is
/// the field's. The field has its own text, from the later file's own
/// description, after the statistic's, which takes the later file only where
/// the later source asks for the statistic.
#[test]
fn a_field_a_later_source_names_is_described_beside_the_earlier_sources_statistic() {
    let asked = statistic("count", FIRST);
    let field = second_field("count");
    holds(
        "statistic-then-field",
        "count",
        [".", "."],
        ("summary_stats=count", "fields=count"),
        &format!("{asked}; {field}"),
        &[asked.clone(), field.clone()],
        (1, "1"),
    );
    holds(
        "statistic-then-both",
        "count",
        [".", "."],
        ("summary_stats=count", "fields=count,summary_stats=count"),
        &format!("{asked},{SECOND}; {field}"),
        &[format!("{asked},{SECOND}"), field],
        (1, "1"),
    );
}

/// An earlier source declares the key twice, as the field and as the
/// statistic, and a later source declares it again. Ensembl VEP appends the
/// later file to the later of the key's two pairs, the statistic's
/// (`OutputFactory.pm:2452`, `:2483-2484`, `:2500-2501`), once for each
/// declaration. The later file joins the field's text where the later source
/// names the field and the statistic's where it asks for the statistic.
#[test]
fn a_later_sources_file_joins_the_declaration_it_repeats() {
    let field = first_field("count");
    let asked = statistic("count", FIRST);
    holds(
        "both-then-field",
        "count",
        [".", "."],
        ("fields=count,summary_stats=count", "fields=count"),
        &format!("{asked}; {field},{SECOND}"),
        &[format!("{field},{SECOND}"), asked.clone()],
        (2, "1"),
    );
    holds(
        "both-then-statistic",
        "count",
        [".", "."],
        ("fields=count,summary_stats=count", "summary_stats=count"),
        &format!("{asked},{SECOND}; {field}"),
        &[field.clone(), format!("{asked},{SECOND}")],
        (1, "1"),
    );
    holds(
        "both-then-both",
        "count",
        [".", "."],
        (
            "fields=count,summary_stats=count",
            "fields=count,summary_stats=count",
        ),
        &format!("{asked},{SECOND}; {field},{SECOND}"),
        &[format!("{field},{SECOND}"), format!("{asked},{SECOND}")],
        (2, "1"),
    );
}

/// `min` is assigned only where a matched record of the asking source carries
/// a score (`AnnotationSource/File.pm:424-434`), a VCF record's `QUAL`. With
/// scores of 30 and 12 the cell holds the asking source's statistic; with no
/// score it holds the field's values, of both files or of the one whose
/// source names the field. The key's one line carries the statistic's text
/// and the field's in every run, each with the files of the sources that
/// declare it, so one of its two texts is true of the cell in each. Ensembl
/// VEP's line for the key is the first source's text alone, the later file
/// appended to the first header line (the key holds no `o`,
/// `OutputFactory.pm:2452`).
#[test]
fn the_line_names_what_the_cell_holds_with_a_score_and_without_one() {
    let first = first_field("min");
    let second = second_field("min");
    for (scored, quals) in [(false, [".", "."]), (true, ["30", "12"])] {
        let run = |label: &str| format!("{label}-{}", if scored { "scored" } else { "unscored" });
        let cell = |with_a_score: &'static str, without: &'static str| {
            if scored {
                with_a_score
            } else {
                without
            }
        };
        let asked = statistic("min", SECOND);
        holds(
            &run("field-then-both"),
            "min",
            quals,
            ("fields=min", "fields=min,summary_stats=min"),
            &format!("{asked}; {first},{SECOND}"),
            &[format!("{first},{SECOND}"), asked.clone()],
            (2, cell("12", "7,9")),
        );
        holds(
            &run("field-then-statistic"),
            "min",
            quals,
            ("fields=min", "summary_stats=min"),
            &format!("{asked}; {first}"),
            &[first.clone(), asked],
            (1, cell("12", "7")),
        );
        let asked = statistic("min", FIRST);
        holds(
            &run("statistic-then-field"),
            "min",
            quals,
            ("summary_stats=min", "fields=min"),
            &format!("{asked}; {second}"),
            &[asked.clone(), second.clone()],
            (1, cell("30", "9")),
        );
        holds(
            &run("both-then-field"),
            "min",
            quals,
            ("fields=min,summary_stats=min", "fields=min"),
            &format!("{asked}; {first},{SECOND}"),
            &[format!("{first},{SECOND}"), asked],
            (2, cell("30", "7,9")),
        );
    }
}

/// Under `type=exact` each file's record names the allele `C`. The row of `C`
/// holds the statistic and then each source's field value, in the sources'
/// order, and the row of `G` the statistic alone: the order the key's one
/// line gives its texts. Ensembl VEP writes its header and ends its run at
/// the push of a field's value onto the statistic (`OutputFactory.pm:1303`),
/// the documented difference `custom_key_named_as_a_plain_value_field`.
#[test]
fn a_row_of_the_allele_both_records_name_holds_the_statistic_then_each_field() {
    let tmp = tempfile::tempdir().unwrap();
    let files = site_files(tmp.path(), "count", [".", "."]);
    let field = first_field("count");
    let asked = statistic("count", SECOND);
    let sources = sources(
        &files,
        "exact",
        "fields=count",
        "fields=count,summary_stats=count",
    );
    assert_eq!(
        written(tmp.path(), "exact", "tab", &sources),
        Written {
            key_lines: vec![
                line("S", NAMES),
                line("S_count", &format!("{asked}; {field},{SECOND}")),
            ],
            columns: "\tS\tS_count".to_string(),
            cells: vec![
                "\tMODIFIER\t-\t-\t-\t-\texact_2,second_1\t1,7,9".to_string(),
                "\tMODIFIER\t-\t-\t-\t-\t-\t1".to_string(),
            ],
        }
    );
}

/// A field that an earlier source names twice and a later source once more.
/// Ensembl VEP pushes the first source's two pairs and appends the later file
/// to the later of them (`OutputFactory.pm:2451`, `:2452`, `:2483-2484`), and
/// every text names the field, whose values the cell holds, so the lines stand
/// as Ensembl VEP's pairs give them: its default header's two lines, and both
/// texts on the key's one line of the tab and VCF headers, the later first.
#[test]
fn a_field_named_twice_and_once_more_by_a_later_source_keeps_its_lines() {
    let field = first_field("count");
    holds(
        "field-twice-then-once",
        "count",
        [".", "."],
        ("fields=count%count", "fields=count"),
        &format!("{field},{SECOND}; {field}"),
        &[field.clone(), format!("{field},{SECOND}")],
        (2, "7,9"),
    );
}
