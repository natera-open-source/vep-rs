// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` key that spells a statistic of a source, on the rows the
//! binary writes for one record: line 9 of the custom corpus's input, for
//! whose site the corpus's `custom/sites.vcf.gz` holds a record on one of the
//! two alleles. A source named `S` that asks for `summary_stats=count`
//! declares the key `S_count`, which a second source named `S_count` lists
//! its record names under and under which a field named `count` of the first
//! source is listed as well. Ensembl VEP 116.3 assigns the statistic as a
//! plain number on the input record's hash (`OutputFactory.pm:1307-1309`,
//! called at `:1004-1014`) and ends its run where a row pushes a record of
//! its allele under the same key (`:1301`, `:1303`, called at `:1426-1435`),
//! so the rows held for that shape are derived: the statistic, then the
//! values Ensembl VEP pushes after it. Where Ensembl VEP completes, under
//! records that name no allele, in JSON and on a `--most_severe` row, the
//! rows are the ones it writes, its second listing of the key left out (the
//! documented difference `duplicate_output_column_collapsed`).

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The 1-based line of the record in the custom corpus's `variants.vcf`.
const RECORD_LINE: usize = 9;

/// The key the statistic and the second source share.
const KEY: &str = "S_count";

/// The thirteen columns of the record's two rows, which no case changes.
const ROW_OF_THE_NAMED_ALLELE: &str =
    "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-";
const ROW_OF_THE_OTHER_ALLELE: &str =
    "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-";

/// The record Ensembl VEP 116.3 writes in JSON under a source `S` of
/// `type=exact` that asks for `count` and a source `S_count` of `type=exact`,
/// its `input` key aside: the statistic opens the list of `S` and the record
/// of the second source stands under `S_count`.
const BOTH_EXACT_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{"S":[{"count":1},{"allele":"C","name":"exact_2"}],"S_count":[{"allele":"C","name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

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

/// The custom corpus's site file, whose record at the site is named `exact_2`
/// and names the allele `C`.
fn site_file() -> PathBuf {
    corpus("GRCh38-custom").join("custom/sites.vcf.gz")
}

/// A site file written under `dir` whose one record is the corpus's record at
/// the site with an INFO field named `count`, the name of a statistic.
fn site_file_with_a_count_field(dir: &Path) -> PathBuf {
    common::indexed_vcf(
        dir,
        "count_field.vcf.gz",
        &[
            "##fileformat=VCFv4.2",
            "##INFO=<ID=count,Number=1,Type=Integer,Description=\"A field named as a statistic\">",
            "##contig=<ID=21>",
        ],
        "21\t5380280\texact_2\tA\tC\t.\tPASS\tcount=7",
    )
}

/// The source `S` over `file`, which asks for the statistic `count`.
fn statistic_source(file: &Path, source_type: &str) -> String {
    format!(
        "file={},short_name=S,format=vcf,type={source_type},summary_stats=count",
        file.display()
    )
}

/// The source `S_count` over `file`, whose record names are listed under the
/// key of the first source's statistic.
fn source_named_as_the_statistic(file: &Path, source_type: &str) -> String {
    format!(
        "file={},short_name=S_count,format=vcf,type={source_type}",
        file.display()
    )
}

/// The source `S` over `file` that asks for the statistic `count` and for the
/// record field of the same name.
fn source_with_a_count_field(file: &Path, source_type: &str) -> String {
    format!(
        "file={},short_name=S,format=vcf,type={source_type},fields=count,summary_stats=count",
        file.display()
    )
}

/// The output of a run over the record in `format` (`default`, `tab`, `vcf`,
/// `json` or `parquet`) under `sources`, in their order, and `flags`; `label`
/// names the run's output file.
fn annotate(dir: &Path, label: &str, format: &str, sources: &[String], flags: &[&str]) -> PathBuf {
    let out = dir.join(format!("{label}.{format}.out"));
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(record_file(dir))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus("GRCh38-everything").join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
        .args(flags);
    for source in sources {
        cmd.args(["--custom", source]);
    }
    match format {
        "default" => {}
        "parquet" => {
            cmd.args(["--output_format", "parquet"]);
        }
        other => {
            cmd.arg(format!("--{other}"));
        }
    }
    let run = cmd.output().unwrap();
    assert!(
        run.status.success(),
        "vep {format} {label} {flags:?}:\n{}",
        String::from_utf8_lossy(&run.stderr)
    );
    out
}

/// The cells the tab output of the run prints under [`KEY`], row by row, with
/// the number of columns its column line lists under that name.
fn tab_cells_under_the_key(dir: &Path, label: &str, sources: &[String]) -> (usize, Vec<String>) {
    tab_cells_under(dir, label, sources, KEY)
}

/// The cells the tab output of the run prints under `key`, row by row, with
/// the number of columns its column line lists under that name.
fn tab_cells_under(dir: &Path, label: &str, sources: &[String], key: &str) -> (usize, Vec<String>) {
    let text = fs::read_to_string(annotate(dir, label, "tab", sources, &[])).unwrap();
    let mut lines = text.lines().filter(|l| !l.starts_with("##"));
    let columns: Vec<&str> = lines
        .next()
        .expect("the column line")
        .trim_start_matches('#')
        .split('\t')
        .collect();
    let listed = columns.iter().filter(|c| **c == key).count();
    let place = columns
        .iter()
        .position(|c| *c == key)
        .expect("the key's column");
    let cells = lines
        .map(|row| {
            row.split('\t')
                .nth(place)
                .expect("the key's cell")
                .to_string()
        })
        .collect();
    (listed, cells)
}

/// The statistic is a plain number on the input record's hash before any row
/// pushes the records of its allele (`OutputFactory.pm:1004-1014`,
/// `:1307-1309`, then `:1426-1435`), which ends Ensembl VEP's run at the push
/// (`:1301` for a record's name, `:1303` for a field). The tab output holds
/// the statistic and the pushed value after it on the row of the allele the
/// record names, and the statistic alone on the other row, in one column.
#[test]
fn a_custom_key_named_as_a_statistic_prints_the_statistic_then_the_records_of_the_rows_allele() {
    let tmp = tempfile::tempdir().unwrap();
    let sites = site_file();
    let with_count = site_file_with_a_count_field(tmp.path());
    for (label, sources, cells) in [
        (
            "both_exact",
            vec![
                statistic_source(&sites, "exact"),
                source_named_as_the_statistic(&sites, "exact"),
            ],
            ["1,exact_2", "1"],
        ),
        (
            "statistic_overlap_names_exact",
            vec![
                statistic_source(&sites, "overlap"),
                source_named_as_the_statistic(&sites, "exact"),
            ],
            ["1,exact_2", "1"],
        ),
        (
            "count_field_exact",
            vec![source_with_a_count_field(&with_count, "exact")],
            ["1,7", "1"],
        ),
    ] {
        let (listed, printed) = tab_cells_under_the_key(tmp.path(), label, &sources);
        assert_eq!(listed, 1, "{label}: the key heads one column");
        assert_eq!(printed, cells, "{label}");
    }
}

/// Ensembl VEP assigns every statistic a source closed (`OutputFactory.pm:1307-1309`):
/// `count` for a source that asks for it, and `min`, `max`, `mean` and `sum`
/// where a matched record carries a score (`AnnotationSource/File.pm:424-434`).
/// The custom corpus's `custom/regions.bed.gz` holds a record at the site with
/// the score 4, so beside a source `R` over it that asks for all five, a source
/// named as any of the five keys ends Ensembl VEP 116.3's run, and the cell is
/// the statistic and then the record of the row's allele. The records of the
/// site file carry no score, so a source `S` over it closes no `min`: Ensembl
/// VEP 116.3 completes and the key `S_min` lists the second source's record
/// alone, the row held here.
#[test]
fn each_statistic_a_source_closes_is_followed_by_the_records_of_the_rows_allele() {
    let tmp = tempfile::tempdir().unwrap();
    let sites = site_file();
    let scored = format!(
        "file={},short_name=R,format=bed,type=overlap,summary_stats=min%max%mean%count%sum",
        corpus("GRCh38-custom")
            .join("custom/regions.bed.gz")
            .display()
    );
    for (statistic, value) in [
        ("min", "4"),
        ("max", "4"),
        ("mean", "4"),
        ("count", "1"),
        ("sum", "4"),
    ] {
        let key = format!("R_{statistic}");
        let named = format!(
            "file={},short_name={key},format=vcf,type=exact",
            sites.display()
        );
        let (listed, printed) = tab_cells_under(
            tmp.path(),
            &format!("scored_{statistic}"),
            &[scored.clone(), named],
            &key,
        );
        assert_eq!(listed, 1, "{key}: the key heads one column");
        assert_eq!(
            printed,
            [format!("{value},exact_2"), value.to_string()],
            "{key}"
        );
    }

    let unscored = format!(
        "file={},short_name=S,format=vcf,type=exact,summary_stats=min%max%mean%sum%count",
        sites.display()
    );
    let named = format!(
        "file={},short_name=S_min,format=vcf,type=exact",
        sites.display()
    );
    let (listed, printed) =
        tab_cells_under(tmp.path(), "unscored_min", &[unscored, named], "S_min");
    assert_eq!(listed, 1, "S_min: the key heads one column");
    assert_eq!(printed, ["exact_2", "-"], "S_min");
}

/// A record that names no allele is pushed on the input record's hash, where
/// the statistic is assigned over it (`OutputFactory.pm:1303` then
/// `:1307-1309` for a field of the statistic's own source; for another
/// source, in the order in which Ensembl VEP 116.3 completes). The cells are
/// the ones every completed run of Ensembl VEP writes: the statistic alone on
/// both rows, in either order of the two sources.
#[test]
fn records_of_the_whole_input_record_stay_under_the_statistic() {
    let tmp = tempfile::tempdir().unwrap();
    let sites = site_file();
    let with_count = site_file_with_a_count_field(tmp.path());
    for (label, sources) in [
        (
            "both_overlap",
            vec![
                statistic_source(&sites, "overlap"),
                source_named_as_the_statistic(&sites, "overlap"),
            ],
        ),
        (
            "both_overlap_names_first",
            vec![
                source_named_as_the_statistic(&sites, "overlap"),
                statistic_source(&sites, "overlap"),
            ],
        ),
        (
            "statistic_exact_names_overlap",
            vec![
                statistic_source(&sites, "exact"),
                source_named_as_the_statistic(&sites, "overlap"),
            ],
        ),
        (
            "count_field_overlap",
            vec![source_with_a_count_field(&with_count, "overlap")],
        ),
    ] {
        let (listed, printed) = tab_cells_under_the_key(tmp.path(), label, &sources);
        assert_eq!(listed, 1, "{label}: the key heads one column");
        assert_eq!(printed, ["1", "1"], "{label}");
    }
}

/// The whole rows of the default format and the record line of the VCF under
/// the two `exact` sources: the default format joins the statistic and the
/// record's name with a comma under the key's place among the Extra keys, and
/// a `CSQ` entry joins them with `&`.
#[test]
fn the_default_and_vcf_outputs_join_the_statistic_and_the_record() {
    let tmp = tempfile::tempdir().unwrap();
    let sites = site_file();
    let sources = [
        statistic_source(&sites, "exact"),
        source_named_as_the_statistic(&sites, "exact"),
    ];

    let text =
        fs::read_to_string(annotate(tmp.path(), "both_exact", "default", &sources, &[])).unwrap();
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        rows,
        [
            format!("{ROW_OF_THE_NAMED_ALLELE}\tIMPACT=MODIFIER;S=exact_2;S_count=1,exact_2"),
            format!("{ROW_OF_THE_OTHER_ALLELE}\tIMPACT=MODIFIER;S_count=1"),
        ]
    );

    let text =
        fs::read_to_string(annotate(tmp.path(), "both_exact", "vcf", &sources, &[])).unwrap();
    let record = text
        .lines()
        .find(|l| !l.starts_with('#'))
        .expect("the record's line");
    assert_eq!(
        record,
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||||||||||||||||||||||exact_2|1&exact_2,G|intergenic_variant|MODIFIER|||||||||||||||||||||||1"
    );
}

/// Ensembl VEP 116.3 completes the run of the two `exact` sources in JSON,
/// whose factory takes the custom annotations and their statistics off the
/// record before any hash is built (`OutputFactory/JSON.pm:249-258`), and on
/// a `--most_severe` row, which is the input record's hash alone and makes no
/// push of an allele's records (`OutputFactory.pm:387`): the record and the
/// row are the ones it writes.
#[test]
fn the_json_record_and_a_most_severe_row_hold_the_statistic_alone() {
    let tmp = tempfile::tempdir().unwrap();
    let sites = site_file();
    let sources = [
        statistic_source(&sites, "exact"),
        source_named_as_the_statistic(&sites, "exact"),
    ];

    let expected: serde_json::Value = serde_json::from_str(BOTH_EXACT_JSON).unwrap();
    let text =
        fs::read_to_string(annotate(tmp.path(), "both_exact", "json", &sources, &[])).unwrap();
    let mut record: serde_json::Value =
        serde_json::from_str(text.lines().next().expect("the record's line")).unwrap();
    let object = record.as_object_mut().expect("a JSON object");
    object.remove("input");
    object.remove("vep_rs_version");
    assert_eq!(record, expected);

    let text = fs::read_to_string(annotate(
        tmp.path(),
        "both_exact_most_severe",
        "tab",
        &sources,
        &["--most_severe"],
    ))
    .unwrap();
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        rows,
        ["rs1444031311\t21:5380280\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t1"]
    );
}

/// The Parquet output stores the key as a string that holds the `--tab` cell:
/// the statistic and the record's name after it on the row of the allele the
/// record names, the statistic on the other row.
#[cfg(feature = "parquet")]
#[test]
fn the_parquet_output_stores_the_tab_cell() {
    let tmp = tempfile::tempdir().unwrap();
    let sites = site_file();
    let sources = [
        statistic_source(&sites, "exact"),
        source_named_as_the_statistic(&sites, "exact"),
    ];
    let text =
        fs::read_to_string(annotate(tmp.path(), "both_exact", "tab", &sources, &[])).unwrap();
    let place = text
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .expect("the column line")
        .trim_start_matches('#')
        .split('\t')
        .position(|column| column == KEY)
        .expect("the key's column");
    let out = annotate(tmp.path(), "both_exact", "parquet", &sources, &[]);
    let mut cells: Vec<String> = vep_io::output::parquet::tab_rows(&out)
        .unwrap()
        .iter()
        .map(|row| {
            row.split('\t')
                .nth(place)
                .expect("the key's cell")
                .to_string()
        })
        .collect();
    cells.sort_unstable();
    assert_eq!(cells, ["1", "1,exact_2"]);
}
