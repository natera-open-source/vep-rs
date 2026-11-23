// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `compare` program end to end on two-record fixtures, and the registry held to the
//! Ensembl VEP references of the golden corpora and to the support table of
//! `docs/plugins.md`: a field a plugin writes that no entry names is a field the comparison
//! never counts.

use flate2::read::MultiGzDecoder;
use regex::Regex;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::OnceLock;
use vep_plugin_pairs::registry::{
    all_plugin_names, core_named_columns, field_sets, plugin_alias, plugin_of_field, FieldMatcher,
};

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";

fn row(location: &str, allele: &str, feature: &str, extra: &str) -> String {
    format!(
        "v\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\tmissense_variant\t10\t10\t4\tA/T\tgcc/Acc\t-\t{extra}"
    )
}

/// A directory of this test's own, empty at the start.
fn scratch(test: &str) -> PathBuf {
    let directory = std::env::temp_dir().join(format!(
        "vep-plugin-pairs-compare-{}-{test}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&directory);
    fs::create_dir_all(&directory).unwrap();
    directory
}

fn write<S: AsRef<str>>(path: PathBuf, rows: &[S]) -> PathBuf {
    let mut text = format!("{HEADER}\n");
    for row in rows {
        text.push_str(row.as_ref());
        text.push('\n');
    }
    fs::write(&path, text).unwrap();
    path
}

fn run(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vep-plugin-pairs"))
        .arg("compare")
        .args(arguments)
        .output()
        .unwrap()
}

fn text(path: &Path) -> &str {
    path.to_str().unwrap()
}

/// The two-record pair whose sides agree.
fn pair(directory: &Path) -> (PathBuf, PathBuf) {
    let rows = [
        row("21:1000", "T", "ENST1", "IMPACT=MODERATE;REVEL=0.512"),
        row("21:2000", "G", "ENST2", "IMPACT=MODERATE;REVEL=0.100"),
    ];
    (
        write(directory.join("perl.txt"), &rows),
        write(directory.join("rust.txt"), &rows),
    )
}

/// Twelve records whose REVEL value differs between the two sides on every record: two
/// more mismatches than the default sample keeps.
fn twelve(directory: &Path, records: usize) -> (PathBuf, PathBuf) {
    let rows = |value: &str| -> Vec<String> {
        (0..records)
            .map(|i| {
                row(
                    &format!("21:{}", 1000 + i),
                    "T",
                    &format!("ENST{i}"),
                    &format!("REVEL={value}"),
                )
            })
            .collect()
    };
    (
        write(directory.join("perl.txt"), &rows("0.100")),
        write(directory.join("rust.txt"), &rows("0.900")),
    )
}

/// A JSON value as the report holds them.
#[derive(Debug, PartialEq)]
enum Json {
    Null,
    Bool(bool),
    Number(f64),
    Text(String),
    List(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl std::ops::Index<&str> for Json {
    type Output = Json;
    fn index(&self, name: &str) -> &Json {
        match self {
            Json::Object(members) => members
                .iter()
                .find(|(known, _)| known == name)
                .map(|(_, value)| value)
                .unwrap_or_else(|| panic!("no member {name}")),
            other => panic!("{other:?} has no member {name}"),
        }
    }
}

impl Json {
    fn items(&self) -> &[Json] {
        match self {
            Json::List(items) => items,
            other => panic!("{other:?} is no list"),
        }
    }

    fn texts(&self) -> Vec<&str> {
        self.items()
            .iter()
            .map(|item| match item {
                Json::Text(text) => text.as_str(),
                other => panic!("{other:?} is no string"),
            })
            .collect()
    }

    fn count(&self) -> u64 {
        match self {
            Json::Number(number) => *number as u64,
            other => panic!("{other:?} is no number"),
        }
    }
}

struct JsonReader<'a> {
    text: &'a [u8],
    at: usize,
}

impl JsonReader<'_> {
    fn peek(&mut self) -> u8 {
        while self.text[self.at].is_ascii_whitespace() {
            self.at += 1;
        }
        self.text[self.at]
    }

    fn word(&mut self, word: &str, value: Json) -> Json {
        assert_eq!(&self.text[self.at..self.at + word.len()], word.as_bytes());
        self.at += word.len();
        value
    }

    fn string(&mut self) -> String {
        assert_eq!(self.peek(), b'"');
        self.at += 1;
        let mut units: Vec<u16> = Vec::new();
        loop {
            let byte = self.text[self.at];
            self.at += 1;
            match byte {
                b'"' => break,
                b'\\' => {
                    let escape = self.text[self.at];
                    self.at += 1;
                    units.push(match escape {
                        b'n' => 10,
                        b't' => 9,
                        b'r' => 13,
                        b'b' => 8,
                        b'f' => 12,
                        b'u' => {
                            let digits = std::str::from_utf8(&self.text[self.at..self.at + 4]);
                            self.at += 4;
                            u16::from_str_radix(digits.unwrap(), 16).unwrap()
                        }
                        other => other as u16,
                    });
                }
                other => units.push(other as u16),
            }
        }
        String::from_utf16(&units).unwrap()
    }

    fn value(&mut self) -> Json {
        match self.peek() {
            b'n' => self.word("null", Json::Null),
            b't' => self.word("true", Json::Bool(true)),
            b'f' => self.word("false", Json::Bool(false)),
            b'"' => Json::Text(self.string()),
            b'[' => {
                self.at += 1;
                let mut items = Vec::new();
                while self.peek() != b']' {
                    items.push(self.value());
                    if self.peek() == b',' {
                        self.at += 1;
                    }
                }
                self.at += 1;
                Json::List(items)
            }
            b'{' => {
                self.at += 1;
                let mut members = Vec::new();
                while self.peek() != b'}' {
                    let name = self.string();
                    assert_eq!(self.peek(), b':');
                    self.at += 1;
                    members.push((name, self.value()));
                    if self.peek() == b',' {
                        self.at += 1;
                    }
                }
                self.at += 1;
                Json::Object(members)
            }
            _ => {
                let from = self.at;
                while self.at < self.text.len()
                    && matches!(
                        self.text[self.at],
                        b'0'..=b'9' | b'-' | b'+' | b'.' | b'e' | b'E'
                    )
                {
                    self.at += 1;
                }
                let number = std::str::from_utf8(&self.text[from..self.at]).unwrap();
                Json::Number(number.parse().unwrap())
            }
        }
    }
}

fn report(path: &Path) -> Json {
    let text = fs::read(path).unwrap();
    JsonReader { text: &text, at: 0 }.value()
}

fn stdout(output: &Output) -> String {
    String::from_utf8(output.stdout.clone()).unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8(output.stderr.clone()).unwrap()
}

#[test]
fn json_output_aggregates_carry_annotated_counts() {
    let directory = scratch("aggregates");
    let (perl, rust) = pair(&directory);
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL",
        "--json-output",
        text(&out),
    ]);
    assert_eq!(done.status.code(), Some(0), "{}", stderr(&done));
    let data = report(&out);
    let aggregate = &data["aggregates"]["REVEL"];
    assert_eq!(aggregate["perl_annotated"].count(), 2);
    assert_eq!(aggregate["rust_annotated"].count(), 2);
    assert_eq!(aggregate["match"].count(), 2);
    assert_eq!(aggregate["mismatch"].count(), 0);
    assert_eq!(aggregate["pass"], Json::Bool(true));
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn a_tolerance_half_way_between_two_texts_is_written_with_the_even_digit() {
    let directory = scratch("tolerance-tie");
    let (perl, rust) = pair(&directory);
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL",
        "--tolerance",
        "1125899906842624.25",
        "--json-output",
        text(&out),
    ]);
    assert_eq!(done.status.code(), Some(0), "{}", stderr(&done));
    assert!(
        stdout(&done).contains("\nTolerance: 1125899906842624.2\n"),
        "{}",
        stdout(&done)
    );
    let written = fs::read_to_string(&out).unwrap();
    assert!(
        written.contains("\"tolerance\": 1125899906842624.2,\n"),
        "{written}"
    );
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn json_output_flags_a_mismatch_and_a_rust_only_record() {
    let directory = scratch("flags");
    let perl = write(
        directory.join("perl.txt"),
        &[row("21:1000", "T", "ENST1", "REVEL=0.512")],
    );
    let rust = write(
        directory.join("rust.txt"),
        &[
            row("21:1000", "T", "ENST1", "REVEL=0.900"),
            row("21:3000", "C", "ENST3", "REVEL=0.300"),
        ],
    );
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL",
        "--json-output",
        text(&out),
    ]);
    assert_ne!(done.status.code(), Some(0));
    let data = report(&out);
    let aggregate = &data["aggregates"]["REVEL"];
    assert_eq!(aggregate["mismatch"].count(), 1);
    assert_eq!(aggregate["rust_only"].count(), 1);
    assert_eq!(aggregate["pass"], Json::Bool(false));
    let _ = fs::remove_dir_all(directory);
}

/// UTRAnnotator.pm:282-283 walk the upstream ORFs of a row and each one's pairs as hashes,
/// so Ensembl VEP writes them in another order from one process to the next; vep-rs writes
/// each group's pairs in the order of their keys.
#[test]
fn a_5utr_annotation_differing_in_pair_order_alone_is_a_match() {
    let directory = scratch("pair-order");
    let hash_order = "type=uORF:KozakStrength=Weak:DistanceToCDS=40,DistanceToCDS=7:KozakStrength=Strong:type=OutOfFrame_oORF";
    let key_order = "DistanceToCDS=40:KozakStrength=Weak:type=uORF,DistanceToCDS=7:KozakStrength=Strong:type=OutOfFrame_oORF";
    let groups_swapped = "type=OutOfFrame_oORF:DistanceToCDS=7:KozakStrength=Strong,KozakStrength=Weak:type=uORF:DistanceToCDS=40";
    let another_value = key_order.replace("KozakStrength=Weak", "KozakStrength=Moderate");
    assert_ne!(another_value, key_order);

    // The match and mismatch counts of one record a side, and the exit status.
    let compared = |case: &str, perl: PathBuf, rust: PathBuf| -> (u64, u64, Option<i32>) {
        let out = directory.join(format!("{case}.json"));
        let done = run(&[
            "--perl-output",
            text(&perl),
            "--rust-output",
            text(&rust),
            "--plugins",
            "UTRAnnotator",
            "--json-output",
            text(&out),
        ]);
        let data = report(&out);
        let aggregate = &data["aggregates"]["UTRAnnotator"];
        (
            aggregate["match"].count(),
            aggregate["mismatch"].count(),
            done.status.code(),
        )
    };
    let default_pair = |case: &str, field: &str, perl: &str, rust: &str| {
        let side = |name: &str, value: &str| {
            let extra = format!("IMPACT=MODIFIER;{field}={value}");
            write(
                directory.join(format!("{case}.{name}.txt")),
                &[row("21:1000", "T", "ENST1", &extra)],
            )
        };
        compared(case, side("perl", perl), side("rust", rust))
    };
    let annotation = "5UTR_annotation";
    assert_eq!(
        default_pair("pairs", annotation, hash_order, key_order),
        (1, 0, Some(0)),
        "the pairs of one group in another order"
    );
    assert_eq!(
        default_pair("groups", annotation, groups_swapped, key_order),
        (1, 0, Some(0)),
        "the groups, and the pairs of both, in another order"
    );
    assert_eq!(
        default_pair("value", annotation, hash_order, &another_value),
        (0, 1, Some(1)),
        "one value of one pair differs"
    );
    assert_eq!(
        default_pair("field", "5UTR_consequence", hash_order, key_order),
        (0, 1, Some(1)),
        "the same two texts under another field of the plugin"
    );
    assert_eq!(
        default_pair("candidate", annotation, key_order, hash_order),
        (0, 1, Some(1)),
        "vep-rs's own pairs out of the order of their keys"
    );

    // A VCF value joins its groups by `&`, since a comma parts the entries of a CSQ.
    let vcf_pair = |case: &str, perl: &str, rust: &str| {
        let side = |name: &str, value: &str| {
            let path = directory.join(format!("{case}.{name}.vcf"));
            let text = format!(
                "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations. Format: Allele|Feature|5UTR_annotation\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t1000\t.\tC\tT\t.\t.\tCSQ=T|ENST1|{}\n",
                value.replace(',', "&")
            );
            fs::write(&path, text).unwrap();
            path
        };
        compared(case, side("perl", perl), side("rust", rust))
    };
    assert_eq!(
        vcf_pair("vcf-groups", groups_swapped, key_order),
        (1, 0, Some(0)),
        "a VCF value with its groups and pairs in another order"
    );
    assert_eq!(
        vcf_pair("vcf-value", hash_order, &another_value),
        (0, 1, Some(1)),
        "a VCF value of which one pair differs"
    );
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn unknown_plugin_name_is_refused() {
    let directory = scratch("unknown");
    let (perl, rust) = pair(&directory);
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVELL",
    ]);
    assert_eq!(done.status.code(), Some(2));
    assert!(stderr(&done).contains("unknown plugin"));
    assert!(
        stderr(&done).starts_with("ERROR: unknown plugin(s): REVELL. Known: AVI, AlphaMissense, ")
    );
    let _ = fs::remove_dir_all(directory);
}

fn max_mismatches_caps_the_json_sample(case: &str, value: &str, kept: usize) {
    let directory = scratch(case);
    let (perl, rust) = twelve(&directory, 12);
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL",
        "--json-output",
        text(&out),
        "--max-mismatches",
        value,
    ]);
    assert_eq!(done.status.code(), Some(1));
    let data = report(&out);
    let entry = &data["results"]["perl"]["REVEL"];
    let locations: Vec<String> = entry["mismatches_sample"]
        .items()
        .iter()
        .map(|mismatch| mismatch["key"].texts()[0].to_string())
        .collect();
    let expected: Vec<String> = (0..kept).map(|i| format!("21:{}", 1000 + i)).collect();
    assert_eq!(locations, expected);
    assert_eq!(entry["total_mismatches"].count(), 12);
    assert_eq!(entry["stats"]["mismatch"].count(), 12);
    let aggregate = &data["aggregates"]["REVEL"];
    assert_eq!(aggregate["mismatch"].count(), 12);
    assert_eq!(aggregate["pass"], Json::Bool(false));
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn max_mismatches_caps_the_json_sample_under_the_default() {
    max_mismatches_caps_the_json_sample("under", "3", 3);
}

#[test]
fn max_mismatches_caps_the_json_sample_over_the_default() {
    max_mismatches_caps_the_json_sample("over", "11", 11);
}

#[test]
fn max_mismatches_caps_the_json_sample_zero_keeps_none() {
    max_mismatches_caps_the_json_sample("zero", "0", 0);
}

#[test]
fn negative_max_mismatches_is_refused_by_name() {
    let directory = scratch("negative");
    let (perl, rust) = pair(&directory);
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL",
        "--json-output",
        text(&out),
        "--max-mismatches",
        "-1",
    ]);
    assert_eq!(done.status.code(), Some(2));
    assert!(stderr(&done).contains("argument --max-mismatches"));
    assert!(
        stderr(&done).ends_with("error: argument --max-mismatches: must be 0 or greater, got -1\n")
    );
    assert!(!out.exists());
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn max_mismatches_caps_each_file_pair_in_directory_mode() {
    let directory = scratch("directories");
    let (perl_dir, rust_dir) = (directory.join("perl"), directory.join("rust"));
    fs::create_dir_all(&perl_dir).unwrap();
    fs::create_dir_all(&rust_dir).unwrap();
    let (perl, rust) = twelve(&directory, 12);
    fs::rename(&perl, perl_dir.join("first.txt")).unwrap();
    fs::rename(&rust, rust_dir.join("first.txt")).unwrap();
    let (perl, rust) = twelve(&directory, 5);
    fs::rename(&perl, perl_dir.join("second.txt")).unwrap();
    fs::rename(&rust, rust_dir.join("second.txt")).unwrap();
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-dir",
        text(&perl_dir),
        "--rust-dir",
        text(&rust_dir),
        "--plugins",
        "REVEL",
        "--json-output",
        text(&out),
        "--max-mismatches",
        "4",
    ]);
    assert_eq!(done.status.code(), Some(1));
    let data = report(&out);
    let first = &data["results"]["first"]["REVEL"];
    let second = &data["results"]["second"]["REVEL"];
    assert_eq!(first["mismatches_sample"].items().len(), 4);
    assert_eq!(first["total_mismatches"].count(), 12);
    assert_eq!(second["mismatches_sample"].items().len(), 4);
    assert_eq!(second["total_mismatches"].count(), 5);
    assert_eq!(data["aggregates"]["REVEL"]["mismatch"].count(), 17);
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn printed_report_shows_one_mismatch_under_max_mismatches_1() {
    let directory = scratch("printed");
    let (perl, rust) = twelve(&directory, 12);
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL",
        "--max-mismatches",
        "1",
    ]);
    assert_eq!(done.status.code(), Some(1));
    let printed = stdout(&done);
    assert!(printed.contains("Mismatch: 12"));
    let shown: Vec<&str> = printed
        .lines()
        .filter(|line| line.ends_with(": mismatch"))
        .map(str::trim)
        .collect();
    assert_eq!(shown, ["21:1000 T ENST0: mismatch"]);
    let _ = fs::remove_dir_all(directory);
}

/// Control: the default sample of ten holds whether or not the option is read.
#[test]
fn control_absent_max_mismatches_keeps_ten() {
    let directory = scratch("control");
    let (perl, rust) = twelve(&directory, 12);
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL",
        "--json-output",
        text(&out),
    ]);
    assert_eq!(done.status.code(), Some(1));
    let data = report(&out);
    let entry = &data["results"]["perl"]["REVEL"];
    assert_eq!(entry["mismatches_sample"].items().len(), 10);
    assert_eq!(entry["total_mismatches"].count(), 12);
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn gnomadc_fields_are_counted_under_the_names_ensembl_vep_writes() {
    let directory = scratch("gnomadc");
    let extra = "IMPACT=MODIFIER;gnomAD_10x_cov=0.0332;gnomAD_mean_cov=2.1875";
    let other = row("21:3000", "C", "ENST3", "IMPACT=MODERATE;REVEL=0.100");
    let perl = write(
        directory.join("perl.txt"),
        &[
            row("21:1000", "T", "ENST1", extra),
            row("21:2000", "G", "ENST2", extra),
            other.clone(),
        ],
    );
    let differing = format!(
        "{};gnomAD_total_DP_cov=7",
        extra.replace("0.0332", "0.5000")
    );
    let rust = write(
        directory.join("rust.txt"),
        &[
            row("21:1000", "T", "ENST1", extra),
            row("21:2000", "G", "ENST2", &differing),
            other,
        ],
    );
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "gnomADc",
        "--json-output",
        text(&out),
    ]);
    assert_ne!(done.status.code(), Some(0));
    let data = report(&out);
    let aggregate = &data["aggregates"]["gnomADc"];
    assert_eq!(aggregate["perl_annotated"].count(), 2);
    assert_eq!(aggregate["rust_annotated"].count(), 2);
    assert_eq!(aggregate["match"].count(), 1);
    assert_eq!(aggregate["mismatch"].count(), 1);
    assert_eq!(
        data["results"]["perl"]["gnomADc"]["fields"].texts(),
        ["gnomAD_mean_cov", "gnomAD_total_DP_cov", "gnomAD_10x_cov"]
    );
    let _ = fs::remove_dir_all(directory);
}

fn write_vcf(path: PathBuf, entries: &[(&str, &str)]) -> PathBuf {
    let mut text = String::from(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations. Format: Allele|Feature|REVEL|gnomAD_10x_cov\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for (position, csq) in entries {
        text.push_str(&format!("21\t{position}\t.\tC\tT\t.\t.\tCSQ={csq}\n"));
    }
    fs::write(&path, text).unwrap();
    path
}

#[test]
fn vcf_outputs_are_compared_on_the_csq_fields_a_plugin_writes() {
    let directory = scratch("vcf");
    let unannotated = [("3000", "T|ENST3||"), ("4000", "T|ENST4|.|.")];
    let perl = write_vcf(
        directory.join("perl.vcf"),
        &[
            ("1000", "T|ENST1|0.512|0.0332"),
            ("2000", "T|ENST2||0.0332"),
            unannotated[0],
            unannotated[1],
        ],
    );
    let rust = write_vcf(
        directory.join("rust.vcf"),
        &[
            ("1000", "T|ENST1|0.512|0.0332"),
            ("2000", "T|ENST2||0.5000"),
            unannotated[0],
            unannotated[1],
        ],
    );
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL,gnomADc",
        "--json-output",
        text(&out),
    ]);
    assert_ne!(done.status.code(), Some(0));
    let data = report(&out);
    let (revel, gnomadc) = (&data["aggregates"]["REVEL"], &data["aggregates"]["gnomADc"]);
    assert_eq!(revel["perl_annotated"].count(), 1);
    assert_eq!(revel["rust_annotated"].count(), 1);
    assert_eq!(revel["pass"], Json::Bool(true));
    assert_eq!(gnomadc["perl_annotated"].count(), 2);
    assert_eq!(gnomadc["rust_annotated"].count(), 2);
    assert_eq!(gnomadc["match"].count(), 1);
    assert_eq!(gnomadc["mismatch"].count(), 1);
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn a_name_two_requested_plugins_write_counts_for_one_of_them() {
    let directory = scratch("shared");
    let rows = [
        row(
            "21:1000",
            "T",
            "ENST1",
            "IMPACT=MODERATE;SIFT_pred=T;REVEL_score=0.512",
        ),
        row("21:2000", "G", "ENST2", "IMPACT=MODERATE;SIFT_pred=D"),
    ];
    let perl = write(directory.join("perl.txt"), &rows);
    let rust = write(directory.join("rust.txt"), &rows);
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "TranscriptAnnotator,dbNSFP",
        "--json-output",
        text(&out),
    ]);
    assert_ne!(done.status.code(), Some(0));
    let data = report(&out);
    let dbnsfp = &data["aggregates"]["dbNSFP"];
    let annotator = &data["aggregates"]["TranscriptAnnotator"];
    assert_eq!(dbnsfp["perl_annotated"].count(), 2);
    assert_eq!(dbnsfp["rust_annotated"].count(), 2);
    assert_eq!(dbnsfp["pass"], Json::Bool(true));
    assert_eq!(annotator["perl_annotated"].count(), 0);
    assert_eq!(annotator["pass"], Json::Bool(false));
    let results = &data["results"]["perl"];
    assert_eq!(
        results["dbNSFP"]["fields"].texts(),
        ["REVEL_score", "SIFT_pred"]
    );
    assert_eq!(results["dbNSFP"]["shared_fields"].texts(), ["SIFT_pred"]);
    assert_eq!(
        results["TranscriptAnnotator"]["fields"].texts(),
        ["transcript_annotator_<column>", "Comment"]
    );
    assert_eq!(
        results["TranscriptAnnotator"]["shared_fields"].texts(),
        ["SIFT_pred"]
    );
    assert_eq!(
        stdout(&done)
            .matches("Also entries of another requested plugin: SIFT_pred\n")
            .count(),
        2
    );

    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "TranscriptAnnotator",
        "--json-output",
        text(&out),
    ]);
    assert_eq!(done.status.code(), Some(0));
    let data = report(&out);
    assert_eq!(
        data["aggregates"]["TranscriptAnnotator"]["perl_annotated"].count(),
        2
    );
    assert!(
        data["results"]["perl"]["TranscriptAnnotator"]["shared_fields"]
            .items()
            .is_empty()
    );
    assert!(!stdout(&done).contains("Also entries of another requested plugin"));
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn the_comparison_leaves_a_core_named_column_to_the_core_key() {
    let directory = scratch("core-named");
    let rows = [row(
        "21:1000",
        "T",
        "ENST1",
        "IMPACT=MODERATE;APPRIS=principal;TSL=1",
    )];
    let perl = write(directory.join("perl.txt"), &rows);
    let rust = write(directory.join("rust.txt"), &rows);
    let out = directory.join("result.json");
    let done = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "dbNSFP",
        "--json-output",
        text(&out),
    ]);
    assert_ne!(done.status.code(), Some(0));
    let data = report(&out);
    let aggregate = &data["aggregates"]["dbNSFP"];
    assert_eq!(aggregate["perl_annotated"].count(), 0);
    assert_eq!(aggregate["rust_annotated"].count(), 0);
    assert_eq!(aggregate["pass"], Json::Bool(false));
    let _ = fs::remove_dir_all(directory);
}

/// Rows of four plugins over a few hundred keys, most keys written twice so that a later
/// row stands for an earlier one, with a value that differs between the sides on every
/// seventh key and keys that one side alone holds.
fn many_rows(side: usize) -> Vec<String> {
    let mut rows = Vec::new();
    for i in 0..1500usize {
        let key = (i * 7919) % 600;
        if key % 97 == side {
            continue;
        }
        let differing = if key % 7 == 0 { side } else { 0 };
        let extra = match i % 5 {
            0 => format!("IMPACT=LOW;REVEL=0.{key};CADD_RAW={i}"),
            1 => format!("REVEL=0.{key}{differing};gnomAD_{}x_cov=0.5", 10 + key % 3),
            2 => "IMPACT=MODIFIER".to_string(),
            3 => format!(
                "CADD_PHRED={key};CADD_RAW={differing}.25;dbNSFP_col{}=v{key}",
                key % 4
            ),
            _ => format!("REVEL={key}e-3;SIFT_pred=T{differing};REVEL=0.{key}"),
        };
        let location = format!("{}:{}", 1 + key % 22, 1000 + key);
        let allele = ["T", "G", "-", "\u{e9}"][key % 4];
        rows.push(row(&location, allele, &format!("ENST{}", key % 50), &extra));
    }
    rows
}

#[test]
fn the_same_pair_gives_the_same_outputs_on_any_thread_count_and_with_small_shards() {
    let directory = scratch("threads");
    let perl = write(directory.join("perl.txt"), &many_rows(0));
    let rust = write(directory.join("rust.txt"), &many_rows(1));
    let out = directory.join("result.json");
    let outputs = |more: &[&str]| {
        let _ = fs::remove_file(&out);
        let mut arguments = vec![
            "--perl-output",
            text(&perl),
            "--rust-output",
            text(&rust),
            "--plugins",
            "REVEL,CADD,gnomADc,dbNSFP",
            "--max-mismatches",
            "7",
            "--json-output",
            text(&out),
        ];
        arguments.extend_from_slice(more);
        let done = run(&arguments);
        assert_eq!(done.status.code(), Some(1), "{more:?}: {}", stderr(&done));
        assert_eq!(stderr(&done), "", "{more:?}");
        (stdout(&done), fs::read(&out).unwrap())
    };
    let one = outputs(&["--threads", "1"]);
    let data = report(&out);
    assert!(data["results"]["perl"]["REVEL"]["stats"]["mismatch"].count() > 10);
    assert!(data["results"]["perl"]["CADD"]["stats"]["match"].count() > 100);
    assert!(data["results"]["perl"]["dbNSFP"]["stats"]["perl_only"].count() > 0);
    assert!(data["results"]["perl"]["gnomADc"]["stats"]["rust_only"].count() > 0);
    assert_eq!(
        data["results"]["perl"]["REVEL"]["mismatches_sample"]
            .items()
            .len(),
        7
    );
    for threads in ["1", "2", "5", "16"] {
        for sizes in [
            &[][..],
            &["--block-bytes", "512", "--shard-bytes", "4096"],
            &["--block-bytes", "100", "--shard-bytes", "64"],
            &["--block-bytes", "1", "--shard-bytes", "1000000"],
            &["--shard-bytes", "1"],
        ] {
            let mut more = vec!["--threads", threads];
            more.extend_from_slice(sizes);
            assert!(outputs(&more) == one, "{more:?}");
        }
    }
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn the_same_directories_give_the_same_outputs_on_any_thread_count() {
    let directory = scratch("directories-threads");
    let (perl_dir, rust_dir) = (directory.join("perl"), directory.join("rust"));
    fs::create_dir_all(&perl_dir).unwrap();
    fs::create_dir_all(&rust_dir).unwrap();
    for (name, from) in [("a", 0), ("b", 500), ("c", 1000)] {
        write(
            perl_dir.join(format!("perl_{name}.txt")),
            &many_rows(0)[from..from + 400],
        );
        write(
            rust_dir.join(format!("rust_{name}.txt")),
            &many_rows(1)[from..from + 400],
        );
    }
    let out = directory.join("result.json");
    let outputs = |more: &[&str]| {
        let _ = fs::remove_file(&out);
        let mut arguments = vec![
            "--perl-dir",
            text(&perl_dir),
            "--rust-dir",
            text(&rust_dir),
            "--plugins",
            "CADD,REVEL",
            "--json-output",
            text(&out),
        ];
        arguments.extend_from_slice(more);
        let done = run(&arguments);
        assert_eq!(done.status.code(), Some(1), "{more:?}: {}", stderr(&done));
        (stdout(&done), fs::read(&out).unwrap())
    };
    let one = outputs(&["--threads", "1"]);
    assert_eq!(report(&out)["file_count"].count(), 3);
    for more in [
        &["--threads", "2"][..],
        &["--threads", "3", "--shard-bytes", "200"],
        &[
            "--threads",
            "16",
            "--block-bytes",
            "64",
            "--shard-bytes",
            "1",
        ],
    ] {
        assert!(outputs(more) == one, "{more:?}");
    }
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn the_scratch_is_made_under_the_work_directory_and_removed_when_the_run_ends() {
    let directory = scratch("work-dir");
    let work = directory.join("work");
    fs::create_dir_all(&work).unwrap();
    let (perl, rust) = pair(&directory);
    let empty = |work: &Path| fs::read_dir(work).unwrap().count() == 0;

    let passing = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL",
        "--work-dir",
        text(&work),
    ]);
    assert_eq!(passing.status.code(), Some(0), "{}", stderr(&passing));
    assert!(empty(&work));

    let (perl, rust) = twelve(&directory, 12);
    let failing = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL",
        "--work-dir",
        text(&work),
        "--shard-bytes",
        "16",
    ]);
    assert_eq!(failing.status.code(), Some(1));
    assert!(empty(&work));

    let unreadable = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&directory.join("missing.txt")),
        "--plugins",
        "REVEL",
        "--work-dir",
        text(&work),
    ]);
    assert_eq!(unreadable.status.code(), Some(1));
    assert!(stderr(&unreadable).starts_with("ERROR: [vep-plugin-pairs] "));
    assert!(stderr(&unreadable).contains("missing.txt: "));
    assert!(empty(&work));

    let nowhere = run(&[
        "--perl-output",
        text(&perl),
        "--rust-output",
        text(&rust),
        "--plugins",
        "REVEL",
        "--work-dir",
        text(&directory.join("no-such-directory")),
    ]);
    assert_eq!(nowhere.status.code(), Some(1));
    assert!(stderr(&nowhere).starts_with("ERROR: [vep-plugin-pairs] "));
    assert!(stdout(&nowhere).starts_with("Plugins: REVEL\nFile pairs: 1\n"));
    let _ = fs::remove_dir_all(directory);
}

#[test]
fn a_termination_signal_ends_the_run_with_its_scratch_removed() {
    use std::io::Write;
    use std::os::unix::ffi::OsStrExt;
    use std::process::Stdio;
    use std::time::{Duration, Instant};

    let directory = scratch("signal");
    let work = directory.join("work");
    fs::create_dir_all(&work).unwrap();
    let (_, rust) = pair(&directory);
    // A pipe that is never closed stands for an output that is still being read.
    let pipe = directory.join("piped.txt");
    let name = std::ffi::CString::new(pipe.as_os_str().as_bytes()).unwrap();
    // SAFETY: `name` is a path of this test's own directory, NUL-terminated.
    assert_eq!(unsafe { libc::mkfifo(name.as_ptr(), 0o600) }, 0);
    let child = Command::new(env!("CARGO_BIN_EXE_vep-plugin-pairs"))
        .args([
            "compare",
            "--perl-output",
            text(&pipe),
            "--rust-output",
            text(&rust),
            "--plugins",
            "REVEL",
            "--threads",
            "3",
            "--work-dir",
            text(&work),
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut writer = fs::OpenOptions::new().write(true).open(&pipe).unwrap();
    writeln!(writer, "{}", row("21:1000", "T", "ENST1", "REVEL=0.5")).unwrap();
    writer.flush().unwrap();
    let started = Instant::now();
    while fs::read_dir(&work).unwrap().count() == 0 {
        assert!(
            started.elapsed() < Duration::from_secs(30),
            "no scratch was made"
        );
        std::thread::sleep(Duration::from_millis(10));
    }
    // The run has read what the pipe holds once nothing of it is left to write.
    std::thread::sleep(Duration::from_millis(200));
    // SAFETY: the process is this test's own child, which has not been waited for.
    assert_eq!(
        unsafe { libc::kill(child.id() as libc::pid_t, libc::SIGTERM) },
        0
    );
    let done = child.wait_with_output().unwrap();
    drop(writer);
    assert_eq!(done.status.code(), Some(130), "{}", stderr(&done));
    assert_eq!(
        stderr(&done),
        "ERROR: [vep-plugin-pairs] interrupted before the report was written\n"
    );
    assert_eq!(
        stdout(&done),
        "Plugins: REVEL\nFile pairs: 1\nTolerance: 1e-06\n"
    );
    assert_eq!(fs::read_dir(&work).unwrap().count(), 0);
    let _ = fs::remove_dir_all(directory);
}

fn repository() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// A default-format reference of the golden corpora.
struct Reference {
    /// Its path under the corpus root.
    name: String,
    /// The plugins its command line names; `None` without a command line.
    plugins: Option<Vec<String>>,
    /// The words of its command line.
    flags: Vec<String>,
    /// Each Extra key its header declares, and its description.
    declared: Vec<(String, String)>,
    /// The declared keys and the keys its rows carry.
    keys: BTreeSet<String>,
}

fn references() -> &'static [Reference] {
    static REFERENCES: OnceLock<Vec<Reference>> = OnceLock::new();
    REFERENCES.get_or_init(|| {
        let golden = repository().join("tests/golden/116");
        let plugin_word = Regex::new(r"--plugin\s+([^,\s]+)").unwrap();
        let mut paths = Vec::new();
        for corpus in fs::read_dir(&golden).unwrap() {
            let corpus = corpus.unwrap().path();
            let Ok(entries) = fs::read_dir(&corpus) else {
                continue;
            };
            for expected in entries {
                let expected = expected.unwrap();
                let reference = expected.path().join("default.txt.gz");
                if expected
                    .file_name()
                    .to_string_lossy()
                    .starts_with("expected")
                    && reference.is_file()
                {
                    paths.push(reference);
                }
            }
        }
        paths.sort();
        let mut found = Vec::new();
        for path in paths {
            let mut text = String::new();
            MultiGzDecoder::new(fs::File::open(&path).unwrap())
                .read_to_string(&mut text)
                .unwrap();
            let (mut plugins, mut flags, mut declared) = (None, Vec::new(), Vec::new());
            let mut keys = BTreeSet::new();
            let mut declaring = false;
            for line in text.split_inclusive('\n') {
                if line.starts_with("## VEP command-line:") {
                    plugins = Some(
                        plugin_word
                            .captures_iter(line)
                            .map(|found| found[1].to_string())
                            .collect::<Vec<String>>(),
                    );
                    let words = line.split_once(':').unwrap().1;
                    flags = words.split_whitespace().map(str::to_string).collect();
                    declaring = false;
                } else if line.starts_with("## Extra column keys:") {
                    declaring = true;
                } else if line.starts_with("##") {
                    let described = line
                        .get(3..)
                        .unwrap_or("")
                        .trim_end_matches('\n')
                        .split_once(" : ");
                    if let (true, Some((key, description))) = (declaring, described) {
                        declared.push((key.to_string(), description.to_string()));
                        keys.insert(key.to_string());
                    }
                } else if !line.starts_with('#') {
                    let cells: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
                    if cells.len() >= 14 && cells[13] != "-" {
                        keys.extend(
                            cells[13]
                                .split(';')
                                .filter_map(|pair| pair.split_once('='))
                                .map(|(key, _)| key.to_string()),
                        );
                    }
                }
            }
            let name = path.strip_prefix(&golden).unwrap().to_string_lossy();
            found.push(Reference {
                name: name.into_owned(),
                plugins,
                flags,
                declared,
                keys,
            });
        }
        assert!(found.len() > 300, "{} references read", found.len());
        found
    })
}

fn loads_no_plugin(reference: &Reference) -> bool {
    reference.plugins.as_ref().is_some_and(Vec::is_empty)
}

/// The Extra keys of the references whose command line loads no plugin.
fn core_keys() -> BTreeSet<&'static str> {
    references()
        .iter()
        .filter(|reference| loads_no_plugin(reference))
        .flat_map(|reference| reference.keys.iter().map(String::as_str))
        .collect()
}

/// Each core key with the descriptions the references that load no plugin give it.
fn core_descriptions() -> BTreeMap<&'static str, BTreeSet<&'static str>> {
    let mut described: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    for reference in references()
        .iter()
        .filter(|reference| loads_no_plugin(reference))
    {
        for (key, description) in &reference.declared {
            described.entry(key).or_default().insert(description);
        }
    }
    described
}

fn registry_name(loaded: &str) -> &str {
    plugin_alias(loaded).unwrap_or(loaded)
}

fn plugin_matchers() -> Vec<FieldMatcher> {
    field_sets()
        .iter()
        .filter(|set| all_plugin_names().contains(&set.name))
        .map(|set| FieldMatcher::new(&set.entries))
        .collect()
}

#[test]
fn every_plugin_key_of_a_golden_reference_is_registered() {
    let core = core_keys();
    assert!(!core.is_empty(), "no reference without a plugin was read");
    let mut unregistered = Vec::new();
    for reference in references() {
        let Some(plugins) = &reference.plugins else {
            assert!(
                !reference.name.contains("-plugin-"),
                "{} has no command line naming its plugins",
                reference.name
            );
            continue;
        };
        for key in reference
            .keys
            .iter()
            .filter(|key| !core.contains(key.as_str()))
        {
            if plugin_of_field(key, plugins, Some(&reference.flags))
                .unwrap()
                .is_none()
            {
                unregistered.push(format!("{}: {key}", reference.name));
            }
        }
    }
    assert!(
        unregistered.is_empty(),
        "not one loaded plugin has an entry for: {}",
        unregistered.join("; ")
    );
    let matchers = plugin_matchers();
    let claimed: Vec<&str> = core
        .iter()
        .copied()
        .filter(|key| matchers.iter().any(|claims| claims.claims(key)))
        .collect();
    assert!(
        claimed.is_empty(),
        "an entry names a key a run without plugins writes: {}",
        claimed.join(", ")
    );
}

/// A core key that a plugin's reference declares under a description no reference without
/// a plugin gives it is a plugin's column under that name: the lookup answers the plugin
/// for it in that reference's run.
#[test]
fn core_named_columns_are_the_ones_the_golden_references_declare() {
    let described = core_descriptions();
    let mut columns: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for reference in references() {
        let Some(plugins) = reference
            .plugins
            .as_ref()
            .filter(|plugins| !plugins.is_empty())
        else {
            continue;
        };
        for (key, description) in &reference.declared {
            let known = described.get(key.as_str());
            if known.is_some_and(|known| !known.contains(description.as_str())) {
                assert_eq!(
                    plugins.len(),
                    1,
                    "{} declares {key} among several plugins",
                    reference.name
                );
                let plugin = registry_name(&plugins[0]).to_string();
                assert_eq!(
                    plugin_of_field(key, plugins, Some(&reference.flags)),
                    Ok(Some(plugin.clone()))
                );
                columns.entry(plugin).or_default().insert(key.clone());
            }
        }
    }
    let mut registered: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for column in core_named_columns() {
        registered
            .entry(column.plugin.clone())
            .or_default()
            .insert(column.column.clone());
    }
    assert_eq!(columns, registered);
}

#[test]
fn the_flags_of_a_core_named_column_are_the_ones_that_turn_its_core_key_on() {
    let described = core_descriptions();
    let with_command_line: Vec<&Reference> = references()
        .iter()
        .filter(|reference| reference.plugins.is_some())
        .collect();
    for column in core_named_columns() {
        let descriptions = &described[column.column.as_str()];
        for reference in &with_command_line {
            let is_core_key = reference.declared.iter().any(|(key, description)| {
                *key == column.column && descriptions.contains(description.as_str())
            });
            let sets_a_flag = column
                .flags
                .iter()
                .any(|flag| reference.flags.contains(flag));
            assert_eq!(is_core_key, sets_a_flag, "{}", reference.name);
        }
        for flag in &column.flags {
            assert!(
                with_command_line
                    .iter()
                    .any(|reference| reference.flags.contains(flag)),
                "{flag}"
            );
        }
    }
}

#[test]
fn every_registry_entry_is_read_from_a_golden_reference() {
    let mut written: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    for reference in references() {
        for plugin in reference.plugins.iter().flatten() {
            written
                .entry(registry_name(plugin))
                .or_default()
                .extend(reference.keys.iter().map(String::as_str));
        }
    }
    let mut unread = Vec::new();
    for set in field_sets()
        .iter()
        .filter(|set| all_plugin_names().contains(&set.name))
    {
        for entry in &set.entries {
            let claims = FieldMatcher::new(&[entry]);
            let read = written
                .get(set.name.as_str())
                .is_some_and(|keys| keys.iter().any(|key| claims.claims(key)));
            if !read {
                unread.push(format!("{}: {entry}", set.name));
            }
        }
    }
    assert!(
        unread.is_empty(),
        "no reference that loads the plugin writes: {}",
        unread.join("; ")
    );
}

#[test]
fn only_the_bare_transcript_annotator_columns_have_two_plugins() {
    let matchers = plugin_matchers();
    let mut candidates: BTreeSet<&str> = references()
        .iter()
        .flat_map(|reference| reference.keys.iter().map(String::as_str))
        .collect();
    for set in field_sets()
        .iter()
        .filter(|set| all_plugin_names().contains(&set.name))
    {
        candidates.extend(set.entries.iter().map(String::as_str));
    }
    let shared: Vec<&str> = candidates
        .into_iter()
        .filter(|key| matchers.iter().filter(|claims| claims.claims(key)).count() > 1)
        .collect();
    assert_eq!(shared, ["SIFT_pred", "SIFT_score"]);
}

/// The cells of a Markdown table row, split on the pipes not escaped with a backslash.
fn table_cells(line: &str) -> Vec<String> {
    let inner = line.trim();
    let inner = inner.strip_prefix('|').unwrap_or(inner);
    let inner = match inner.strip_suffix('|') {
        Some(rest) if !inner.ends_with("\\|") => rest,
        _ => inner,
    };
    let mut cells = vec![String::new()];
    let mut escaped = false;
    for c in inner.chars() {
        match c {
            '|' if escaped => {
                let cell = cells.last_mut().unwrap();
                cell.pop();
                cell.push('|');
            }
            '|' => cells.push(String::new()),
            _ => cells.last_mut().unwrap().push(c),
        }
        escaped = c == '\\';
    }
    cells.iter().map(|cell| cell.trim().to_string()).collect()
}

#[test]
fn the_registry_names_exactly_the_supported_plugins_of_the_support_table() {
    let page = fs::read_to_string(repository().join("docs/plugins.md")).unwrap();
    let lines: Vec<&str> = page.lines().collect();
    let heading = lines
        .iter()
        .position(|line| line.trim() == "## Support table")
        .expect("docs/plugins.md has the support table's heading");
    let header = (heading + 1..lines.len())
        .take_while(|at| !lines[*at].starts_with('#'))
        .find(|at| lines[*at].trim_start().starts_with('|'))
        .expect("a table stands under the heading");
    assert_eq!(table_cells(lines[header]).len(), 6);
    let link = Regex::new(r"^\[([^\]]+)\]\([^)]*\)$").unwrap();
    let mut supported = BTreeSet::new();
    let mut rows = 0;
    for line in lines[header + 2..]
        .iter()
        .take_while(|line| line.trim_start().starts_with('|'))
    {
        let cells = table_cells(line);
        assert_eq!(cells.len(), 6, "{line}");
        let name = match link.captures(&cells[0]) {
            Some(found) => found[1].to_string(),
            None => cells[0].clone(),
        };
        rows += 1;
        if cells[4].starts_with("supported") {
            supported.insert(name);
        }
    }
    assert!(rows > 0, "the support table has no row");
    let registered: BTreeSet<String> = all_plugin_names().iter().cloned().collect();
    assert_eq!(
        supported, registered,
        "the registry has one field set per supported plugin of the support table, naming the Extra keys the plugin's golden references carry"
    );
}
