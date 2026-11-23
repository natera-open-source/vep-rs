// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` key that spells `Uploaded_variation`, `Location` or `Allele`,
//! the three columns that name a row, on one record: line 9 of the custom
//! corpus's input, for whose site the corpus's `custom/sites.vcf.gz` holds a
//! record on one of the two alleles. Ensembl VEP 116.3 holds each of the three
//! as a plain value when it pushes a source's values under the same key
//! (`OutputFactory.pm:895`, `:931` and `:1347` set them, `:1301` and `:1303`
//! push), and its run ends there.
//!
//! The key holds both, so the rows held here are derived as those of the
//! entry's flag fields are: the column's value, then the values pushed after
//! it. `Uploaded_variation` and `Location` are values of the record's hash,
//! set before both pushes; a row assigns `Allele` between them, so the names
//! of the records that name its allele follow it alone. Where Ensembl VEP
//! completes, the rows are the ones it writes, its second listing of the key
//! left out (the documented difference `duplicate_output_column_collapsed`).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Output;

mod common;

/// The 1-based line of the record in the custom corpus's `variants.vcf`.
const RECORD_LINE: usize = 9;

/// The columns of the record's rows from `Gene` to `Existing_variation`.
const INTERGENIC: &str = "-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-";

/// The directory of the custom corpus.
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-custom")
}

/// The record alone, written under `dir`: the header lines of the custom
/// corpus's input and its line [`RECORD_LINE`].
fn record_file(dir: &Path) -> PathBuf {
    let input = fs::read_to_string(corpus().join("variants.vcf")).unwrap();
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

/// One command over the record: the source's short name and type, what
/// follows its `fields` in the `--custom` string, the output format
/// (`default`, `tab`, `vcf`, `json` or `parquet`) and the run's other flags.
struct Run {
    short_name: &'static str,
    source_type: &'static str,
    options: &'static str,
    format: &'static str,
    flags: &'static [&'static str],
}

impl Run {
    const fn of(short_name: &'static str, source_type: &'static str, format: &'static str) -> Self {
        Run {
            short_name,
            source_type,
            options: "",
            format,
            flags: &[],
        }
    }

    const fn under(mut self, flags: &'static [&'static str]) -> Self {
        self.flags = flags;
        self
    }

    fn name(&self) -> String {
        format!(
            "{} {}{} {} {:?}",
            self.short_name, self.source_type, self.options, self.format, self.flags
        )
    }

    fn out(&self, dir: &Path) -> PathBuf {
        dir.join(format!(
            "{}.{}{}.{}.{}.out",
            self.short_name,
            self.source_type,
            self.options.replace([',', '='], "."),
            self.format,
            self.flags.join("")
        ))
    }

    /// The run's result, with the custom corpus's site file as the source and
    /// the corpus's own cache.
    fn result(&self, dir: &Path) -> Output {
        let source = format!(
            "file={},short_name={},format=vcf,type={},fields=AF%CLASS{}",
            corpus().join("custom/sites.vcf.gz").display(),
            self.short_name,
            self.source_type,
            self.options
        );
        let mut cmd = common::isolated_vep(dir);
        cmd.arg("-i")
            .arg(record_file(dir))
            .arg("-o")
            .arg(self.out(dir))
            .arg("--json_cache")
            .arg(corpus().join("json_cache"))
            .args(["--offline", "--species", "homo_sapiens"])
            .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
            .args(["--custom", &source])
            .args(self.flags);
        match self.format {
            "default" => {}
            "parquet" => {
                cmd.args(["--output_format", "parquet"]);
            }
            other => {
                cmd.arg(format!("--{other}"));
            }
        }
        cmd.output().unwrap()
    }

    /// The output file of a run that completes.
    fn written(&self, dir: &Path) -> PathBuf {
        let run = self.result(dir);
        assert!(
            run.status.success(),
            "vep {}:\n{}",
            self.name(),
            String::from_utf8_lossy(&run.stderr)
        );
        self.out(dir)
    }

    /// The lines of the output a completed run writes that are no `##` line.
    fn lines(&self, dir: &Path) -> Vec<String> {
        fs::read_to_string(self.written(dir))
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with("##"))
            .map(str::to_owned)
            .collect()
    }

    /// The first three cells of each row of a completed tab or default run,
    /// with the number of columns the column line lists under `key`.
    fn naming_cells(&self, dir: &Path, key: &str) -> (usize, Vec<String>) {
        let lines = self.lines(dir);
        let listed = lines
            .iter()
            .filter(|line| line.starts_with("#Uploaded_variation"))
            .flat_map(|line| line.trim_start_matches('#').split('\t'))
            .filter(|column| *column == key)
            .count();
        let cells = lines
            .iter()
            .filter(|line| !line.starts_with('#'))
            .map(|row| row.split('\t').take(3).collect::<Vec<_>>().join("\t"))
            .collect();
        (listed, cells)
    }
}

/// The record's hash holds `Location` as a plain string (`OutputFactory.pm:931`)
/// before the values of a record that names no allele are pushed
/// (`:1004-1014`), and each row copies it (`:391`) before the values of a
/// record that names its allele are (`:1426-1435`). Ensembl VEP 116.3 ends its
/// run at the first such push; the key holds both, so the column prints the
/// location and then the names pushed after it, once.
#[test]
fn a_key_named_as_the_location_column_holds_the_location_then_the_sources_names() {
    let tmp = tempfile::tempdir().unwrap();
    for (run, cells) in [
        (
            Run::of("Location", "exact", "tab"),
            vec![
                "rs1444031311\t21:5380280,exact_2\tC",
                "rs1444031311\t21:5380280\tG",
            ],
        ),
        (
            Run::of("Location", "overlap", "tab"),
            vec![
                "rs1444031311\t21:5380280,exact_2\tC",
                "rs1444031311\t21:5380280,exact_2\tG",
            ],
        ),
        (
            Run::of("Location", "overlap", "tab").under(&["--most_severe"]),
            vec!["rs1444031311\t21:5380280,exact_2\t-"],
        ),
        (
            Run::of("Location", "overlap", "tab").under(&["--summary"]),
            vec!["rs1444031311\t21:5380280,exact_2\t-"],
        ),
    ] {
        let (listed, written) = run.naming_cells(tmp.path(), "Location");
        assert_eq!(listed, 1, "{}: the key heads one column", run.name());
        assert_eq!(written, cells, "{}", run.name());
    }
}

/// `Uploaded_variation` is set with `Location` (`OutputFactory.pm:895`) and
/// reads the same way: the identifier, then the names pushed after it.
#[test]
fn a_key_named_as_the_identifier_column_holds_the_identifier_then_the_sources_names() {
    let tmp = tempfile::tempdir().unwrap();
    for (run, cells) in [
        (
            Run::of("Uploaded_variation", "exact", "tab"),
            vec![
                "rs1444031311,exact_2\t21:5380280\tC",
                "rs1444031311\t21:5380280\tG",
            ],
        ),
        (
            Run::of("Uploaded_variation", "overlap", "tab"),
            vec![
                "rs1444031311,exact_2\t21:5380280\tC",
                "rs1444031311,exact_2\t21:5380280\tG",
            ],
        ),
        (
            Run::of("Uploaded_variation", "overlap", "tab").under(&["--most_severe"]),
            vec!["rs1444031311,exact_2\t21:5380280\t-"],
        ),
    ] {
        let (listed, written) = run.naming_cells(tmp.path(), "Uploaded_variation");
        assert_eq!(listed, 1, "{}: the key heads one column", run.name());
        assert_eq!(written, cells, "{}", run.name());
    }
}

/// A `--most_severe` or `--summary` row takes no record that names an allele
/// (`OutputFactory.pm:387` returns it before the rows of the alleles are
/// built), so under a VCF source of `type=exact` Ensembl VEP 116.3 completes
/// on this record with the column's own value, the row held here.
#[test]
fn a_row_of_the_whole_record_keeps_the_columns_own_value_under_an_exact_source() {
    let tmp = tempfile::tempdir().unwrap();
    for key in ["Uploaded_variation", "Location"] {
        for flag in [&["--most_severe"], &["--summary"]] {
            let run = Run::of(key, "exact", "tab").under(flag);
            let (listed, written) = run.naming_cells(tmp.path(), key);
            assert_eq!(listed, 1, "{}", run.name());
            assert_eq!(written, ["rs1444031311\t21:5380280\t-"], "{}", run.name());
        }
    }
}

/// The whole rows of the default format and the record of the VCF under a
/// source named `Location`: the default format prints the joined cell in the
/// column and the key in no `Extra` item, and each `CSQ` entry carries the
/// key as a field of its own, the two values joined with `&`.
#[test]
fn the_default_and_vcf_outputs_join_the_location_and_the_sources_name() {
    let tmp = tempfile::tempdir().unwrap();
    let rows: Vec<String> = Run::of("Location", "exact", "default")
        .lines(tmp.path())
        .into_iter()
        .filter(|line| !line.starts_with('#'))
        .collect();
    assert_eq!(
        rows,
        [
            format!(
                "rs1444031311\t21:5380280,exact_2\tC\t{INTERGENIC}\t\
                 IMPACT=MODIFIER;Location_AF=0.004;Location_CLASS=exact;Location_FILTER=PASS"
            ),
            format!("rs1444031311\t21:5380280\tG\t{INTERGENIC}\tIMPACT=MODIFIER"),
        ]
    );

    let text = fs::read_to_string(Run::of("Location", "exact", "vcf").written(tmp.path())).unwrap();
    let fields = text
        .lines()
        .find(|line| line.starts_with("##INFO=<ID=CSQ,"))
        .expect("the CSQ line");
    assert!(
        fields.ends_with("|SOURCE|Location|Location_AF|Location_CLASS\">"),
        "{fields}"
    );
    assert_eq!(
        text.lines().find(|line| !line.starts_with('#')),
        Some(
            "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER\
             ||||||||||||||||||||||21:5380280&exact_2|0.004|exact,G|intergenic_variant|MODIFIER\
             ||||||||||||||||||||||21:5380280||"
        )
    );
}

/// Lines 3 and 25 of the custom corpus's input are a deletion and an
/// insertion with no transcript in reach, each named by one record of the
/// source on its own allele. Under `--shift_3prime 1 --shift_genomic 1` an
/// intergenic row assigns `Location` again before it pushes the records of
/// its allele (`OutputFactory.pm:2001-2007`, then the call of `:2010`), so the
/// assigned location is followed by that record's name: the deletion moves
/// by no base and the insertion by eight.
#[test]
fn an_intergenic_row_that_assigns_its_location_again_keeps_the_names_of_its_alleles_records() {
    let tmp = tempfile::tempdir().unwrap();
    let input = fs::read_to_string(corpus().join("variants.vcf")).unwrap();
    let mut text: String = input
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    for line in [3, 25] {
        text.push_str(input.lines().nth(line - 1).expect("the record's line"));
        text.push('\n');
    }
    let records = tmp.path().join("indels.vcf");
    fs::write(&records, text).unwrap();
    let location_cells = |flags: &[&str]| {
        let out = tmp.path().join(format!("indels.{}.out", flags.len()));
        let source = format!(
            "file={},short_name=Location,format=vcf,type=exact",
            corpus().join("custom/sites.vcf.gz").display()
        );
        let run = common::isolated_vep(tmp.path())
            .arg("-i")
            .arg(&records)
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(corpus().join("json_cache"))
            .arg("--fasta")
            .arg(corpus().join("reference.fa.gz"))
            .args(["--offline", "--species", "homo_sapiens", "--tab"])
            .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
            .args(["--custom", &source])
            .args(flags)
            .output()
            .unwrap();
        assert!(
            run.status.success(),
            "vep {flags:?}:\n{}",
            String::from_utf8_lossy(&run.stderr)
        );
        fs::read_to_string(&out)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|row| {
                let cells: Vec<&str> = row.split('\t').collect();
                format!("{} {}", cells[1], cells[6])
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(
        location_cells(&[]),
        [
            "21:5028998-5028999,exact_1 intergenic_variant",
            "21:23254313-23254314,exact_5 intergenic_variant",
        ]
    );
    assert_eq!(
        location_cells(&["--shift_3prime", "1", "--shift_genomic", "1"]),
        [
            "21:5028998-5028999,exact_1 intergenic_variant",
            "21:23254321-23254322,exact_5 intergenic_variant",
        ]
    );
}

/// The Parquet output stores the tab output's cells: the joined cell under
/// `Location`.
#[cfg(feature = "parquet")]
#[test]
fn the_parquet_output_stores_the_joined_location() {
    let tmp = tempfile::tempdir().unwrap();
    let out = Run::of("Location", "exact", "parquet").written(tmp.path());
    let mut cells: Vec<String> = vep_io::output::parquet::tab_rows(&out)
        .unwrap()
        .iter()
        .map(|row| row.split('\t').take(3).collect::<Vec<_>>().join("\t"))
        .collect();
    cells.sort_unstable();
    assert_eq!(
        cells,
        [
            "rs1444031311\t21:5380280\tG",
            "rs1444031311\t21:5380280,exact_2\tC",
        ]
    );
}

/// A row assigns `Allele` (`OutputFactory.pm:1347`) over the values of a
/// record that names no allele and before the names of the records that name
/// its allele are pushed (`:1426-1435`), where Ensembl VEP 116.3 ends its run
/// under a VCF source of `type=exact`. The column prints the row's allele and
/// then those names, once, in the tab and default formats, and the first
/// field of a `CSQ` entry joins them with `&`.
#[test]
fn a_key_named_as_the_allele_column_holds_the_rows_allele_then_the_names_of_its_records() {
    let tmp = tempfile::tempdir().unwrap();
    let cells = [
        "rs1444031311\t21:5380280\tC,exact_2",
        "rs1444031311\t21:5380280\tG",
    ];
    for run in [
        Run::of("Allele", "exact", "tab"),
        Run::of("Allele", "exact", "tab").under(&["--flag_pick"]),
        Run::of("Allele", "exact", "default"),
    ] {
        let (listed, written) = run.naming_cells(tmp.path(), "Allele");
        assert_eq!(listed, 1, "{}: the key heads one column", run.name());
        assert_eq!(written, cells, "{}", run.name());
    }

    let rows: Vec<String> = Run::of("Allele", "exact", "default")
        .lines(tmp.path())
        .into_iter()
        .filter(|line| !line.starts_with('#'))
        .collect();
    assert_eq!(
        rows,
        [
            format!(
                "rs1444031311\t21:5380280\tC,exact_2\t{INTERGENIC}\t\
                 IMPACT=MODIFIER;Allele_AF=0.004;Allele_CLASS=exact;Allele_FILTER=PASS"
            ),
            format!("rs1444031311\t21:5380280\tG\t{INTERGENIC}\tIMPACT=MODIFIER"),
        ]
    );

    let text = fs::read_to_string(Run::of("Allele", "exact", "vcf").written(tmp.path())).unwrap();
    let fields = text
        .lines()
        .find(|line| line.starts_with("##INFO=<ID=CSQ,"))
        .expect("the CSQ line");
    assert!(
        fields.contains("Format: Allele|Consequence|")
            && fields.ends_with("|SOURCE|Allele_AF|Allele_CLASS\">"),
        "{fields}"
    );
    // The 26 fields: `Allele`, `Consequence`, `IMPACT`, 21 fields the rows
    // leave empty, then `Allele_AF` and `Allele_CLASS`.
    let unset = "|".repeat(22);
    assert_eq!(
        text.lines().find(|line| !line.starts_with('#')),
        Some(
            format!(
                "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=\
                 C&exact_2|intergenic_variant|MODIFIER{unset}0.004|exact,\
                 G|intergenic_variant|MODIFIER{unset}|"
            )
            .as_str()
        )
    );
}

/// The Parquet output stores the tab output's cells: the joined cell under
/// `Allele`.
#[cfg(feature = "parquet")]
#[test]
fn the_parquet_output_stores_the_joined_allele() {
    let tmp = tempfile::tempdir().unwrap();
    let out = Run::of("Allele", "exact", "parquet").written(tmp.path());
    let mut cells: Vec<String> = vep_io::output::parquet::tab_rows(&out)
        .unwrap()
        .iter()
        .map(|row| row.split('\t').take(3).collect::<Vec<_>>().join("\t"))
        .collect();
    cells.sort_unstable();
    assert_eq!(
        cells,
        [
            "rs1444031311\t21:5380280\tC,exact_2",
            "rs1444031311\t21:5380280\tG",
        ]
    );
}

/// The column line of a tab run over the record under a source named `Allele`.
const ALLELE_COLUMNS: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t\
     Consequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\t\
     Existing_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tAllele_AF\tAllele_CLASS";

/// A `--most_severe` or `--summary` row takes no record that names an allele
/// (`OutputFactory.pm:387`), so under a VCF source of `type=exact` named
/// `Allele` Ensembl VEP 116.3 completes and the row holds no value of the
/// source.
#[test]
fn an_allele_key_of_an_exact_source_runs_on_a_row_of_the_whole_record() {
    let tmp = tempfile::tempdir().unwrap();
    for flag in [&["--most_severe"], &["--summary"]] {
        assert_eq!(
            Run::of("Allele", "exact", "tab")
                .under(flag)
                .lines(tmp.path()),
            [
                ALLELE_COLUMNS,
                "rs1444031311\t21:5380280\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t-\t-\t\
                 -\t-\t-\t-\t-",
            ],
            "{flag:?}"
        );
    }
}

/// Under `num_records=0` an exact VCF source stores one record named `...`
/// that names no allele (`AnnotationSource/File.pm:436-439`): it is pushed
/// before the row assigns its allele, Ensembl VEP 116.3 completes, and the
/// rows are the ones it writes. A `--most_severe` row, which assigns no
/// allele, holds that name.
#[test]
fn an_allele_key_of_an_exact_source_that_lists_no_record_runs_as_ensembl_vep_completes_it() {
    let tmp = tempfile::tempdir().unwrap();
    let listing_none = |flags: &'static [&'static str]| {
        Run {
            options: ",num_records=0",
            ..Run::of("Allele", "exact", "tab")
        }
        .under(flags)
    };
    assert_eq!(
        listing_none(&[]).lines(tmp.path()),
        [
            ALLELE_COLUMNS,
            "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t\
             -\t-\t-\t-\t-\t-",
            "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t\
             -\t-\t-\t-\t-\t-",
        ]
    );
    assert_eq!(
        listing_none(&["--most_severe"]).lines(tmp.path()),
        [
            ALLELE_COLUMNS,
            "rs1444031311\t21:5380280\t...\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t-\t-\t\
             -\t-\t-\t-\t-",
        ]
    );
}

/// The JSON factory takes the custom annotations off the record before any
/// hash is built (`OutputFactory/JSON.pm:249-258`), so Ensembl VEP 116.3
/// completes in JSON under each of the three keys: the source's record stands
/// under `custom_annotations`, with its allele where the source is exact.
#[test]
fn the_json_output_runs_under_each_of_the_three_keys() {
    let tmp = tempfile::tempdir().unwrap();
    for key in ["Uploaded_variation", "Location", "Allele"] {
        for (source_type, allele) in [("exact", r#""allele":"C","#), ("overlap", "")] {
            let expected: serde_json::Value = serde_json::from_str(&format!(
                r#"{{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{{"{key}":[{{{allele}"fields":{{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]}},"name":"exact_2"}}]}},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"}},{{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}}"#
            ))
            .unwrap();
            let lines = Run::of(key, source_type, "json").lines(tmp.path());
            let mut record: serde_json::Value =
                serde_json::from_str(lines.first().expect("the record's line")).unwrap();
            let object = record.as_object_mut().expect("a JSON object");
            object.remove("input");
            object.remove("vep_rs_version");
            assert_eq!(record, expected, "{key} {source_type}");
        }
    }
}
