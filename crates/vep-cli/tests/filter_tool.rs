// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `vep filter`, through the binary: the behaviour of Ensembl VEP's
//! `filter_vep` (release/116.3, whose line numbers the tests cite) on
//! hand-built inputs, and the `tests/golden/116/GRCh38-filter-tool` corpus,
//! whose `tool_cases` each name the arguments of one `filter_vep` run on an
//! annotated output and the file its standard output was recorded in.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-filter-tool")
}

/// Runs `vep filter` with `args`, feeding `stdin` when given.
fn filter_in(dir: &Path, args: &[&str], stdin: Option<&str>) -> Output {
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("filter").args(args).current_dir(dir);
    match stdin {
        Some(text) => {
            let mut child = cmd
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap();
            child
                .stdin
                .take()
                .unwrap()
                .write_all(text.as_bytes())
                .unwrap();
            child.wait_with_output().unwrap()
        }
        None => cmd.stdin(Stdio::null()).output().unwrap(),
    }
}

fn filter(args: &[&str]) -> Output {
    filter_in(&corpus_dir(), args, None)
}

fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The data lines of an output, header lines left out.
fn data_lines(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('#')).collect()
}

/// The output of a run that must succeed quietly.
fn ok(args: &[&str]) -> String {
    let output = filter(args);
    assert!(
        output.status.success(),
        "{args:?} exited {:?}: {}",
        output.status.code(),
        stderr(&output)
    );
    assert_eq!(stderr(&output), "", "{args:?}");
    stdout(&output)
}

/// The one `ERROR: [filter]` line of a refused run, exit 1.
fn refused(args: &[&str]) -> String {
    let output = filter(args);
    assert_eq!(
        output.status.code(),
        Some(1),
        "{args:?}: {}",
        stderr(&output)
    );
    let text = stderr(&output);
    assert!(
        text.starts_with("ERROR: [filter] ") && text.matches('\n').count() == 1,
        "{args:?}: {text:?}"
    );
    text.trim_end().to_string()
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|DISTANCE\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Three records: two entries of which the second passes `IMPACT is HIGH`,
/// one entry that fails, and one line without the consequence field.
const VCF_LINES: &str = "21\t100\tr1\tA\tG\t50\t.\tDP=13;CSQ=G|intron_variant|MODIFIER|TPTE|,G|stop_gained|HIGH|TPTE|\n\
21\t200\tr2\tC\tT\t.\tPASS\tDP=20;DB;CSQ=T|missense_variant|MODERATE|POTED|\n\
21\t300\tr3\tG\tA\t.\t.\tDP=-\n";

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    fs::write(&path, text).unwrap();
    path
}

fn vcf_fixture(dir: &Path) -> PathBuf {
    write(dir, "in.vcf", &format!("{VCF_HEADER}{VCF_LINES}"))
}

/// The header block of an output written with `--keep_csq` over an annotated
/// input: the input's lines, the earlier annotation's `##VEP`, `##INFO` and
/// `##VEP-command-line` among them, then the new annotation's lines. The two
/// Format lists differ, as a re-annotation under more flags leaves them; both
/// engines parse every item's entries under the last list.
const KEEP_CSQ_HEADER: &str = "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##VEP=\"v116\" API=\"v116\" cache=\"homo_sapiens/116_GRCh38\"\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT\">\n\
##VEP-command-line='vep --cache --vcf'\n\
##VEP=\"v116\" API=\"v116\" cache=\"homo_sapiens/116_GRCh38\"\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|DISTANCE\">\n\
##VEP-command-line='vep --cache --vcf --keep_csq --symbol --distance'\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Five records of a `--keep_csq` output, each INFO column holding `DP`, the
/// earlier annotation's `CSQ=` item of three fields, then the new one of five,
/// appended as the writer appends it. Under `IMPACT is HIGH`: r1 passes on its
/// new item alone, r2 on its earlier item alone, r3 on neither, r4 on one
/// entry of each item; r5 has no consequence field.
const KEEP_CSQ_LINES: &str = "21\t100\tr1\tA\tG\t.\t.\tDP=4;CSQ=G|intron_variant|MODIFIER;CSQ=G|stop_gained|HIGH|S1|\n\
21\t200\tr2\tC\tT\t.\tPASS\tDP=5;CSQ=T|stop_gained|HIGH;CSQ=T|missense_variant|MODERATE|S2|\n\
21\t300\tr3\tG\tA\t.\t.\tDP=6;CSQ=A|intron_variant|MODIFIER;CSQ=A|synonymous_variant|LOW|S3|\n\
21\t400\tr4\tT\tC\t.\t.\tDP=7;CSQ=C|stop_gained|HIGH,C|intron_variant|MODIFIER;CSQ=C|missense_variant|MODERATE|S4|,C|stop_gained|HIGH|S4|\n\
21\t500\tr5\tG\tC\t.\t.\tDP=8\n";

fn keep_csq_fixture(dir: &Path) -> PathBuf {
    write(
        dir,
        "keep_csq.vcf",
        &format!("{KEEP_CSQ_HEADER}{KEEP_CSQ_LINES}"),
    )
}

/// The header of a `--keep_csq` output re-annotated under `--af`: the input's
/// `AF` line, the earlier annotation's three-field `CSQ` line and the new
/// four-field one, under which both engines read every entry.
const SHORT_ENTRY_HEADER: &str = "##fileformat=VCFv4.2\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency in the cohort\">\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT\">\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|AF\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Four records under that header: s1 a three-field entry beside an INFO
/// `AF`, s2 an entry writing its own `AF`, s3 an entry writing an empty `AF`,
/// s4 a three-field entry on a line without an INFO `AF`.
const SHORT_ENTRY_LINES: &str = "21\t100\ts1\tA\tG\t.\t.\tAF=0.5;CSQ=G|missense_variant|MODERATE\n\
21\t200\ts2\tC\tT\t.\t.\tAF=0.5;CSQ=T|missense_variant|MODERATE|0.01\n\
21\t300\ts3\tG\tA\t.\t.\tAF=0.5;CSQ=A|missense_variant|MODERATE|\n\
21\t400\ts4\tT\tC\t.\t.\tDP=9;CSQ=C|missense_variant|MODERATE\n";

fn short_entry_fixture(dir: &Path) -> PathBuf {
    write(
        dir,
        "short_entry.vcf",
        &format!("{SHORT_ENTRY_HEADER}{SHORT_ENTRY_LINES}"),
    )
}

/// Every `tool_cases` entry of the corpus manifest reproduces the recorded
/// `filter_vep` output byte for byte, with exit 0 and nothing on standard
/// error; the manifest covers the two matrix rows and every case's expected
/// file is present.
#[test]
fn expected_outputs_match_filter_vep_on_the_corpus() {
    let dir = corpus_dir();
    let manifest: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join("manifest.json")).unwrap()).unwrap();
    let covers: Vec<&str> = manifest["covers"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap())
        .collect();
    assert_eq!(covers, ["filter-vep", "filter-vep-ontology"]);
    let cases = manifest["tool_cases"].as_array().unwrap();
    assert!(cases.len() >= 20, "{} cases", cases.len());
    let mut failures = Vec::new();
    for case in cases {
        let name = case["name"].as_str().unwrap();
        let input = case["input"].as_str().unwrap();
        let expected_path = dir.join(case["expected"].as_str().unwrap());
        let args: Vec<&str> = case["args"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| a.as_str().unwrap())
            .collect();
        let mut argv = vec!["-i", input];
        argv.extend(args);
        let output = filter_in(&dir, &argv, None);
        let Ok(expected) = fs::read(&expected_path) else {
            failures.push(format!(
                "{name}: no expected file {}",
                expected_path.display()
            ));
            continue;
        };
        if !output.status.success() || !output.stderr.is_empty() {
            failures.push(format!(
                "{name}: exit {:?}, stderr {:?}",
                output.status.code(),
                stderr(&output)
            ));
            continue;
        }
        if output.stdout != expected {
            let actual = stdout(&output);
            let expected = String::from_utf8_lossy(&expected);
            let first = actual
                .lines()
                .zip(expected.lines())
                .position(|(a, e)| a != e)
                .unwrap_or(actual.lines().count().min(expected.lines().count()));
            failures.push(format!(
                "{name}: {} lines written, {} expected; first difference at line {}:\n  written:  {:?}\n  expected: {:?}",
                actual.lines().count(),
                expected.lines().count(),
                first + 1,
                actual.lines().nth(first).unwrap_or(""),
                expected.lines().nth(first).unwrap_or("")
            ));
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

/// A VCF line passes when any entry of its consequence field passes
/// (282-291), is written unchanged, and the header block is written first;
/// `--only_matched` keeps the passing entries alone (294-297) and leaves a
/// line whose entries all pass as it was.
#[test]
fn vcf_line_passes_when_any_csq_block_passes() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    let out = ok(&[
        "-i",
        input,
        "--filter",
        "IMPACT is HIGH",
        "--no_missing_consequence_warning",
    ]);
    assert!(out.starts_with(VCF_HEADER), "{out}");
    assert_eq!(data_lines(&out), [VCF_LINES.lines().next().unwrap()]);
    let out = ok(&[
        "-i",
        input,
        "--only_matched",
        "--filter",
        "IMPACT is HIGH or IMPACT is MODERATE",
        "--no_missing_consequence_warning",
    ]);
    assert_eq!(
        data_lines(&out),
        [
            "21\t100\tr1\tA\tG\t50\t.\tDP=13;CSQ=G|stop_gained|HIGH|TPTE|",
            "21\t200\tr2\tC\tT\t.\tPASS\tDP=20;DB;CSQ=T|missense_variant|MODERATE|POTED|",
        ]
    );
}

/// Every INFO item under the field's name is read, so a `--keep_csq` output
/// passes when any entry of its earlier or its new annotation passes, as
/// `filter_vep` has it (261-266, 282-291): its pattern matches each item and
/// every match's entries join one list. The earlier item's entries are read
/// under the last `##INFO` line's Format list, so their fourth field, `SYMBOL`,
/// is undefined (363-366). `--soft_filter` marks the same lines, `--count`
/// counts them, and a line whose earlier item fails and whose new item passes
/// is rewritten under `--only_matched` as `filter_vep` writes it.
#[test]
fn keep_csq_output_is_filtered_on_every_consequence_item() {
    let dir = tempfile::tempdir().unwrap();
    let input = keep_csq_fixture(dir.path());
    let input = input.to_str().unwrap();
    let quiet = "--no_missing_consequence_warning";
    let lines: Vec<&str> = KEEP_CSQ_LINES.lines().collect();
    let out = ok(&["-i", input, "-f", "IMPACT is HIGH", quiet]);
    assert!(out.starts_with(KEEP_CSQ_HEADER), "{out}");
    assert_eq!(data_lines(&out), [lines[0], lines[1], lines[3]]);
    assert_eq!(
        ok(&["-i", input, "-c", "-f", "IMPACT is HIGH", quiet]),
        "3\n"
    );
    let out = ok(&["-i", input, "-f", "IMPACT is HIGH and SYMBOL nex", quiet]);
    assert_eq!(data_lines(&out), [lines[1], lines[3]]);
    let out = ok(&["-i", input, "--soft_filter", "-f", "IMPACT is HIGH", quiet]);
    assert_eq!(
        data_lines(&out),
        [
            "21\t100\tr1\tA\tG\t.\tfilter_vep_pass\tDP=4;CSQ=G|intron_variant|MODIFIER;CSQ=G|stop_gained|HIGH|S1|",
            "21\t200\tr2\tC\tT\t.\tPASS;filter_vep_pass\tDP=5;CSQ=T|stop_gained|HIGH;CSQ=T|missense_variant|MODERATE|S2|",
            "21\t300\tr3\tG\tA\t.\tfilter_vep_fail\tDP=6;CSQ=A|intron_variant|MODIFIER;CSQ=A|synonymous_variant|LOW|S3|",
            "21\t400\tr4\tT\tC\t.\tfilter_vep_pass\tDP=7;CSQ=C|stop_gained|HIGH,C|intron_variant|MODIFIER;CSQ=C|missense_variant|MODERATE|S4|,C|stop_gained|HIGH|S4|",
            "21\t500\tr5\tG\tC\t.\tfilter_vep_fail\tDP=8",
        ]
    );
    let out = ok(&[
        "-i",
        input,
        "--only_matched",
        "-f",
        "IMPACT is HIGH and SYMBOL is S1",
        quiet,
    ]);
    assert_eq!(
        data_lines(&out),
        ["21\t100\tr1\tA\tG\t.\t.\tDP=4;CSQ=;CSQ=G|stop_gained|HIGH|S1|"]
    );
    let output = filter(&["-i", input, "-f", "IMPACT is HIGH"]);
    assert_eq!(
        stderr(&output),
        "WARNING: the VEP annotations INFO field CSQ was not found in 1 line(s) of the input file\n"
    );
}

/// `--only_matched` rewrites each item under the field's name to its own
/// passing entries, an item none of whose entries pass to an empty value.
/// `filter_vep` pairs each passing hash with the chunk at its index in the one
/// chunk list, though with a second item its hashes hold the first item's
/// entries twice, and rewrites the first item alone, only when the passing
/// hashes and the chunks differ in number (263-265, 290, 294-296): on r2 the
/// two hashes of the earlier entry pass, as many as the chunks, so the line
/// stands with its failing new entry; on r4 it writes the first item as
/// `C|stop_gained|HIGH,C|missense_variant|MODERATE|S4|,` and the second
/// unchanged; on r3 under `--soft_filter` it empties the first item alone.
#[test]
fn only_matched_rewrites_each_item_to_its_own_passing_entries() {
    let dir = tempfile::tempdir().unwrap();
    let input = keep_csq_fixture(dir.path());
    let input = input.to_str().unwrap();
    let quiet = "--no_missing_consequence_warning";
    let out = ok(&["-i", input, "--only_matched", "-f", "IMPACT is HIGH", quiet]);
    assert_eq!(
        data_lines(&out),
        [
            "21\t100\tr1\tA\tG\t.\t.\tDP=4;CSQ=;CSQ=G|stop_gained|HIGH|S1|",
            "21\t200\tr2\tC\tT\t.\tPASS\tDP=5;CSQ=T|stop_gained|HIGH;CSQ=",
            "21\t400\tr4\tT\tC\t.\t.\tDP=7;CSQ=C|stop_gained|HIGH;CSQ=C|stop_gained|HIGH|S4|",
        ]
    );
    let out = ok(&[
        "-i",
        input,
        "--soft_filter",
        "--only_matched",
        "-f",
        "IMPACT is HIGH",
        quiet,
    ]);
    assert_eq!(
        data_lines(&out)[2],
        "21\t300\tr3\tG\tA\t.\tfilter_vep_fail\tDP=6;CSQ=;CSQ="
    );
    assert_eq!(
        data_lines(&out)[4],
        "21\t500\tr5\tG\tC\t.\tfilter_vep_fail\tDP=8"
    );
}

/// A consequence field an entry lacks, one past the entry's last `|`, is
/// undefined whatever the line's INFO column holds: `AF > 0.1`, `AF exists`
/// and `AF is 0.5` fail s1 and `AF nex` passes it, as they do s4, while s2
/// reads its own `0.01` and s3's empty `AF` is undefined. `filter_vep` starts
/// each entry's hash as a copy of the line's hash with its INFO keys (389),
/// writes the fields the entry has (395) and sets a header field undefined
/// only `unless exists` (413-416), so on s1 it reads the INFO `AF=0.5` as the
/// consequence field: `s1`; `s1 s2`; `s3 s4`; `s2`; `s1` for the five filters
/// in order, `2` under `--count` and `filter_vep_pass` under `--soft_filter`.
#[test]
fn short_entry_field_is_undefined_not_the_info_key() {
    let dir = tempfile::tempdir().unwrap();
    let input = short_entry_fixture(dir.path());
    let input = input.to_str().unwrap();
    let quiet = "--no_missing_consequence_warning";
    let lines: Vec<&str> = SHORT_ENTRY_LINES.lines().collect();
    let passing = |filter: &str| -> Vec<String> {
        let out = ok(&["-i", input, "-f", filter, quiet]);
        assert!(out.starts_with(SHORT_ENTRY_HEADER), "{filter}: {out}");
        data_lines(&out).iter().map(|l| l.to_string()).collect()
    };
    assert_eq!(passing("AF > 0.1"), Vec::<String>::new());
    assert_eq!(passing("AF exists"), [lines[1]]);
    assert_eq!(passing("AF nex"), [lines[0], lines[2], lines[3]]);
    assert_eq!(passing("AF < 0.1"), [lines[1]]);
    assert_eq!(passing("AF is 0.5"), Vec::<String>::new());
    assert_eq!(ok(&["-i", input, "-c", "-f", "AF exists", quiet]), "1\n");
    let out = ok(&["-i", input, "--soft_filter", "-f", "AF exists", quiet]);
    assert_eq!(
        data_lines(&out),
        [
            "21\t100\ts1\tA\tG\t.\tfilter_vep_fail\tAF=0.5;CSQ=G|missense_variant|MODERATE",
            "21\t200\ts2\tC\tT\t.\tfilter_vep_pass\tAF=0.5;CSQ=T|missense_variant|MODERATE|0.01",
            "21\t300\ts3\tG\tA\t.\tfilter_vep_fail\tAF=0.5;CSQ=A|missense_variant|MODERATE|",
            "21\t400\ts4\tT\tC\t.\tfilter_vep_fail\tDP=9;CSQ=C|missense_variant|MODERATE",
        ]
    );
}

/// `--soft_filter` writes every line, the FILTER column `filter_vep_pass` or
/// `filter_vep_fail`, `PASS` kept ahead of the verdict, adds the two `##FILTER`
/// lines before the column line (201-206, 305-317), writes every line whatever
/// `--start` says (236) and is refused on a tab input (232); with `--count` the
/// lines and the count are both written.
#[test]
fn soft_filter_marks_every_line_and_adds_the_two_filter_header_lines() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    let out = ok(&[
        "-i",
        input,
        "--soft_filter",
        "--start",
        "5",
        "--filter",
        "IMPACT is HIGH",
        "--no_missing_consequence_warning",
    ]);
    let headers: Vec<&str> = out.lines().filter(|l| l.starts_with('#')).collect();
    assert_eq!(
        &headers[headers.len() - 3..],
        [
            "##FILTER=<ID=filter_vep_pass,Description=\"Variant passes filter_vep\">",
            "##FILTER=<ID=filter_vep_fail,Description=\"Variant fails filter_vep\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        ]
    );
    assert_eq!(
        data_lines(&out),
        [
            "21\t100\tr1\tA\tG\t50\tfilter_vep_pass\tDP=13;CSQ=G|intron_variant|MODIFIER|TPTE|,G|stop_gained|HIGH|TPTE|",
            "21\t200\tr2\tC\tT\t.\tPASS;filter_vep_fail\tDP=20;DB;CSQ=T|missense_variant|MODERATE|POTED|",
            "21\t300\tr3\tG\tA\t.\tfilter_vep_fail\tDP=-",
        ]
    );
    let out = ok(&[
        "-i",
        input,
        "--soft_filter",
        "--count",
        "--filter",
        "IMPACT is HIGH",
        "--no_missing_consequence_warning",
    ]);
    assert_eq!(out.lines().count(), 4, "{out}");
    assert_eq!(out.lines().last(), Some("1"));
    assert!(!out.contains("##FILTER"), "{out}");
    let tab = write(
        dir.path(),
        "in.tab",
        "#Uploaded_variation\tLocation\tAllele\tConsequence\nr1\t21:100\tG\tintron_variant\n",
    );
    assert_eq!(
        refused(&[
            "-i",
            tab.to_str().unwrap(),
            "--soft_filter",
            "--filter",
            "Consequence exists"
        ]),
        "ERROR: [filter] --soft_filter is compatible only with VCF files"
    );
    assert_eq!(
        refused(&[
            "-i",
            tab.to_str().unwrap(),
            "--only_matched",
            "--filter",
            "Consequence exists"
        ]),
        "ERROR: [filter] --only_matched is compatible only with VCF files"
    );
}

/// `--count` prints the number of passing lines and nothing else (209, 343),
/// to the output file when one is named.
#[test]
fn count_prints_the_count_alone() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    assert_eq!(
        ok(&[
            "-i",
            input,
            "--count",
            "--filter",
            "IMPACT is HIGH or SYMBOL is POTED",
            "--no_missing_consequence_warning",
        ]),
        "2\n"
    );
    let out_path = dir.path().join("count.txt");
    let out = ok(&[
        "-i",
        input,
        "-c",
        "-f",
        "DP > 15",
        "-o",
        out_path.to_str().unwrap(),
        "--no_missing_consequence_warning",
    ]);
    assert_eq!(out, "");
    assert_eq!(fs::read_to_string(&out_path).unwrap(), "1\n");
}

/// `--list` prints `Available fields:`, a blank line and the sorted union of
/// the consequence fields, the columns and every other `##INFO` id (217-223),
/// on standard output whatever `--output_file` names, with no filter needed
/// and the data lines unread; on an input with header lines alone it prints
/// the same (333-341).
#[test]
fn list_prints_the_sorted_union_of_the_header_fields() {
    let dir = corpus_dir();
    let header = fs::read_to_string(dir.join("annotated.vcf")).unwrap();
    let mut fields = std::collections::BTreeSet::new();
    for line in header.lines().filter(|l| l.starts_with('#')) {
        if let Some(rest) = line.strip_prefix("##INFO=<ID=") {
            let id = rest.split(',').next().unwrap();
            if id == "CSQ" {
                let format = line.split("Format: ").nth(1).unwrap();
                fields.extend(
                    format
                        .trim_end_matches("\">")
                        .split('|')
                        .map(str::to_string),
                );
            } else {
                fields.insert(id.to_string());
            }
        } else if !line.starts_with("##") {
            fields.extend(line[1..].split('\t').map(str::to_string));
        }
    }
    let expected = format!(
        "Available fields:\n\n{}\n",
        fields.iter().cloned().collect::<Vec<_>>().join("\n")
    );
    assert_eq!(ok(&["-i", "annotated.vcf", "--list"]), expected);
    let tmp = tempfile::tempdir().unwrap();
    let named = tmp.path().join("out.vcf");
    assert_eq!(
        ok(&["-i", "annotated.vcf", "-l", "-o", named.to_str().unwrap()]),
        expected
    );
    assert_eq!(fs::read_to_string(&named).unwrap(), "");
    let headers_only = write(tmp.path(), "headers.vcf", VCF_HEADER);
    assert_eq!(
        ok(&["-i", headers_only.to_str().unwrap(), "--list"]),
        "Available fields:\n\nALT\nAllele\nCHROM\nConsequence\nDISTANCE\nDP\nFILTER\nID\nIMPACT\nINFO\nPOS\nQUAL\nREF\nSYMBOL\n"
    );
}

/// `--start N` skips the first N-1 passing lines, `--limit N` stops after N
/// passing lines are written and `--test N` reads N data lines at most
/// (195, 299-301, 322); an input with no data line still writes its headers
/// (326-331).
#[test]
fn start_limit_and_test_bound_the_lines() {
    let dir = tempfile::tempdir().unwrap();
    let lines: String = (1..=6)
        .map(|i| format!("21\t{i}00\tr{i}\tA\tG\t.\t.\tCSQ=G|intron_variant|MODIFIER|S{i}|\n"))
        .collect();
    let input = write(dir.path(), "six.vcf", &format!("{VCF_HEADER}{lines}"));
    let input = input.to_str().unwrap();
    let ids = |out: &str| -> Vec<String> {
        data_lines(out)
            .iter()
            .map(|l| l.split('\t').nth(2).unwrap().to_string())
            .collect()
    };
    assert_eq!(
        ids(&ok(&[
            "-i",
            input,
            "--start",
            "3",
            "--limit",
            "2",
            "-f",
            "IMPACT is MODIFIER"
        ])),
        ["r3", "r4"]
    );
    assert_eq!(
        ids(&ok(&["-i", input, "-s", "5", "-f", "IMPACT is MODIFIER"])),
        ["r5", "r6"]
    );
    assert_eq!(
        ok(&[
            "-i",
            input,
            "--test",
            "4",
            "--count",
            "-f",
            "IMPACT is MODIFIER"
        ]),
        "4\n"
    );
    assert_eq!(
        ids(&ok(&[
            "-i",
            input,
            "--test",
            "2",
            "--limit",
            "1",
            "-f",
            "IMPACT is MODIFIER"
        ])),
        ["r1"]
    );
    let headers_only = write(dir.path(), "headers.vcf", VCF_HEADER);
    assert_eq!(
        ok(&[
            "-i",
            headers_only.to_str().unwrap(),
            "-f",
            "IMPACT is MODIFIER"
        ]),
        VCF_HEADER
    );
}

/// A tab line's `Extra` column is exploded into fields (400-405); `-` and an
/// empty cell are undefined, so `nex` passes them and a unit with a value fails
/// them, while an `Allele` of `-` is defined (407-411); the fields a filter
/// reads are the columns, resolved in any case (216, 552-559).
#[test]
fn tab_extra_column_is_exploded_and_dash_reads_as_undefined() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "in.txt",
        "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n\
## Column descriptions:\n\
## Uploaded_variation : Identifier of uploaded variant\n\
## Extra column keys:\n\
## IMPACT : Subjective impact classification of consequence type\n\
#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tConsequence\tExtra\n\
r1\t21:100\t-\tENSG1\tENST1\tframeshift_variant\tIMPACT=HIGH;STRAND=1;SYMBOL=TPTE\n\
r2\t21:200\tG\t-\t-\tintergenic_variant\tIMPACT=MODIFIER\n\
r3\t21:300\tT\tENSG2\tENST2\tmissense_variant\tIMPACT=MODERATE;STRAND=-1;SYMBOL=\n",
    );
    let input = input.to_str().unwrap();
    let ids = |args: &[&str]| -> Vec<String> {
        let mut argv = vec!["-i", input];
        argv.extend_from_slice(args);
        data_lines(&ok(&argv))
            .iter()
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect()
    };
    assert_eq!(ids(&["-f", "IMPACT is HIGH"]), ["r1"]);
    assert_eq!(ids(&["-f", "STRAND is -1"]), ["r3"]);
    assert_eq!(ids(&["-f", "Allele is -"]), ["r1"]);
    assert_eq!(ids(&["-f", "Gene nex"]), ["r2"]);
    assert_eq!(ids(&["-f", "SYMBOL exists"]), ["r1"]);
    assert_eq!(ids(&["-f", "symbol is TPTE"]), ["r1"]);
    assert_eq!(
        ids(&["-f", "feature is ENST2 and impact is MODERATE"]),
        ["r3"]
    );
    assert_eq!(ids(&["-f", "Gene is -"]), Vec::<String>::new());
    assert_eq!(
        ok(&["-i", input, "--list"]),
        "Available fields:\n\nAllele\nConsequence\nExtra\nFeature\nGene\nIMPACT\nLocation\nUploaded_variation\n"
    );
}

/// A word without `=` inside an `Extra` column is a key with no value:
/// `filter_vep`'s split leaves it undefined (402), as the clean leaves
/// `FLAGX=` (407-411), so `exists` fails it, `nex` passes it, a unit with a
/// value fails it (FilterSet.pm:481-483) and `is` without a value passes it
/// (`eq` over two undefined, FilterSet.pm:704). The key is still the line's,
/// so the prefix `FLAG` resolves to it, and beside a second key it starts is
/// an ambiguous prefix (FilterSet.pm:574-580). Every expected set is
/// `filter_vep`'s on this input.
#[test]
fn tab_extra_bare_word_is_undefined_as_filter_vep_reads_it() {
    let dir = tempfile::tempdir().unwrap();
    let header = "## ENSEMBL VARIANT EFFECT PREDICTOR v116.3\n\
## Extra column keys:\n\
## IMPACT : Subjective impact classification of consequence type\n\
## STRAND : Strand of the feature\n\
#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tConsequence\tExtra\n";
    let ambiguous = write(
        dir.path(),
        "ambiguous.txt",
        &format!("{header}r1\t21:100\tG\tENSG1\tENST1\tmissense_variant\tIMPACT=MODERATE;FLAGX;FLAGY=1\n"),
    );
    assert_eq!(
        refused(&["-i", ambiguous.to_str().unwrap(), "-f", "FLAG"]),
        "ERROR: [filter] Multiple data fields found matching query field \"FLAG\": FLAGX, FLAGY"
    );
    let input = write(
        dir.path(),
        "in.txt",
        &format!(
            "{header}\
r1\t21:100\tG\tENSG1\tENST1\tmissense_variant\tIMPACT=MODERATE;FLAGX;STRAND=1\n\
r2\t21:200\tT\tENSG2\tENST2\tintron_variant\tIMPACT=MODIFIER;STRAND=-1\n\
r3\t21:300\tA\tENSG3\tENST3\tsynonymous_variant\tIMPACT=LOW;FLAGX=;STRAND=1\n\
r4\t21:400\tC\tENSG4\tENST4\tstop_gained\tIMPACT=HIGH;FLAGX=X;STRAND=-1\n\
r5\t21:500\tG\tENSG5\tENST5\tmissense_variant\tFLAGX\n\
r6\t21:600\tT\tENSG6\tENST6\tmissense_variant\tIMPACT=MODERATE;FLAGX=-;STRAND=1\n"
        ),
    );
    let input = input.to_str().unwrap();
    let undefined_on = "r1 r2 r3 r5 r6";
    for (filter, expected) in [
        ("FLAGX", "r4"),
        ("FLAGX exists", "r4"),
        ("FLAGX nex", undefined_on),
        ("FLAGX is", undefined_on),
        ("FLAGX is \"\"", ""),
        ("FLAGX is X", "r4"),
        ("FLAGX ne", undefined_on),
        ("FLAGX ne X", ""),
        ("FLAGX match X", "r4"),
        ("FLAGX in X,Y", "r4"),
        ("not FLAGX", undefined_on),
        ("FLAGX > 0", ""),
        ("FLAGX < 1", "r4"),
        ("FLAGX gte 0", "r4"),
        ("FLAG", "r4"),
        ("STRAND is -1", "r2 r4"),
        ("IMPACT nex", "r5"),
        ("IMPACT is MODERATE and FLAGX nex", "r1 r6"),
    ] {
        let output = ok(&["-i", input, "-f", filter]);
        let ids: Vec<&str> = data_lines(&output)
            .iter()
            .map(|l| l.split('\t').next().unwrap())
            .collect();
        assert_eq!(ids.join(" "), expected, "{filter}");
    }
}

/// The INFO keys of a VCF line are fields of every entry (252-258), read in
/// any case; a flag key without `=` is defined and empty, where `filter_vep`
/// leaves it undefined (253-257). An INFO value of `-` is undefined on an
/// entry and stands on a line without the field (its own hash is not
/// cleaned), where a numeric comparison reads it as 0, as in `filter_vep`.
#[test]
fn info_keys_are_fields_and_a_flag_is_defined_and_empty() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    let ids = |filter: &str| -> Vec<String> {
        data_lines(&ok(&[
            "-i",
            input,
            "-f",
            filter,
            "--no_missing_consequence_warning",
        ]))
        .iter()
        .map(|l| l.split('\t').nth(2).unwrap().to_string())
        .collect()
    };
    assert_eq!(ids("DP > 15"), ["r2"]);
    assert_eq!(ids("dp < 15"), ["r1", "r3"]);
    assert_eq!(ids("DP > 0"), ["r1", "r2"]);
    assert_eq!(ids("DB"), ["r2"]);
    assert_eq!(ids("DB exists and DB is"), ["r2"]);
    assert_eq!(ids("DB nex"), ["r1", "r3"]);
    assert_eq!(ids("DP nex"), Vec::<String>::new());
    assert_eq!(ids("FILTER is PASS"), ["r2"]);
    assert_eq!(ids("POS > 150 and POS < 250"), ["r2"]);
}

/// An input whose first data line arrives before any header line is refused
/// (199); so is a run with neither `--filter` nor `--list` (111), an input
/// file that does not exist (151), an output file that exists without
/// `--force_overwrite` (168-169) and a `--format` that is neither `vcf` nor
/// `tab` (276), each with `filter_vep`'s message.
#[test]
fn input_without_headers_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let bare = write(dir.path(), "bare.vcf", VCF_LINES);
    assert_eq!(
        refused(&["-i", bare.to_str().unwrap(), "-f", "DP > 1"]),
        "ERROR: [filter] No headers found in input file"
    );
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    assert_eq!(
        refused(&["-i", input]),
        "ERROR: [filter] No valid filters given"
    );
    assert_eq!(
        refused(&["-i", "/nonexistent/in.vcf", "-f", "DP > 1"]),
        "ERROR: [filter] Could not find input file /nonexistent/in.vcf"
    );
    let existing = write(dir.path(), "out.vcf", "");
    let existing = existing.to_str().unwrap();
    assert_eq!(
        refused(&["-i", input, "-f", "DP > 1", "-o", existing]),
        format!("ERROR: [filter] Output file {existing} already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite")
    );
    let out = filter(&[
        "-i",
        input,
        "-f",
        "DP > 15",
        "-o",
        existing,
        "--force_overwrite",
        "--no_missing_consequence_warning",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert_eq!(
        data_lines(&fs::read_to_string(existing).unwrap()),
        [VCF_LINES.lines().nth(1).unwrap()]
    );
    assert_eq!(
        refused(&["-i", input, "-f", "DP > 1", "--format", "json"]),
        "ERROR: [filter] Unable to parse data in format json"
    );
    // A filter Ensembl VEP's parser refuses is refused before the input is read.
    assert_eq!(
        refused(&["-i", "/nonexistent/in.vcf", "-f", "DP equals 1"]),
        "ERROR: [filter] No such operator \"equals\""
    );
}

/// A VCF line without the consequence field is evaluated on its columns and
/// INFO keys, whose values stand as written there (its `DP=-` is defined), and
/// the run ends with one warning naming the field and the count of such lines,
/// unless `--no_missing_consequence_warning` (267-271, 346-348); the field's
/// name follows `--vcf_info_field`.
#[test]
fn missing_csq_field_warns_at_the_end_unless_silenced() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    let output = filter(&["-i", input, "-f", "DP exists"]);
    assert!(output.status.success());
    assert_eq!(data_lines(&stdout(&output)).len(), 3);
    assert_eq!(
        stderr(&output),
        "WARNING: the VEP annotations INFO field CSQ was not found in 1 line(s) of the input file\n"
    );
    assert_eq!(
        stderr(&filter(&[
            "-i",
            input,
            "-f",
            "DP exists",
            "--no_missing_consequence_warning"
        ])),
        ""
    );
    let ann = write(
        dir.path(),
        "ann.vcf",
        &format!(
            "{}{}",
            VCF_HEADER.replace("ID=CSQ", "ID=ANN"),
            VCF_LINES.replace("CSQ=", "ANN=")
        ),
    );
    let output = filter(&[
        "-i",
        ann.to_str().unwrap(),
        "--vcf_info_field",
        "ANN",
        "-f",
        "IMPACT is HIGH",
    ]);
    assert_eq!(data_lines(&stdout(&output)).len(), 1);
    assert_eq!(
        stderr(&output),
        "WARNING: the VEP annotations INFO field ANN was not found in 1 line(s) of the input file\n"
    );
}

/// Standard error on a stream that takes no bytes never ends a run of `vep
/// filter`, a device that refuses every write or a pipe whose reader has gone.
/// A run whose one line there is the missing-field warning exits 0 with the
/// output of the same run on a readable stream, and a refused run keeps its
/// exit status of 1 though its `ERROR:` line is lost.
#[test]
fn vep_filter_keeps_its_run_on_an_unwritable_error_stream() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    let run = |args: &[&str], error_stream: Stdio| {
        common::isolated_vep_anywhere()
            .arg("filter")
            .args(args)
            .stdin(Stdio::null())
            .stderr(error_stream)
            .output()
            .unwrap()
    };
    let warns = ["-i", input, "-f", "DP exists"];
    let refusal = ["-i", "/nonexistent/in.vcf", "-f", "DP > 1"];
    let readable = run(&warns, Stdio::piped());
    assert!(readable.status.success(), "{}", stderr(&readable));
    assert_eq!(
        stderr(&readable),
        "WARNING: the VEP annotations INFO field CSQ was not found in 1 line(s) of the input file\n"
    );
    assert_eq!(data_lines(&stdout(&readable)).len(), 3);
    assert_eq!(run(&refusal, Stdio::piped()).status.code(), Some(1));
    for stream in common::unwritable_streams() {
        let output = run(&warns, common::unwritable(stream));
        assert!(output.status.success(), "{stream}: {}", output.status);
        assert_eq!(output.stdout, readable.stdout, "{stream}");
        let refused = run(&refusal, common::unwritable(stream));
        assert_eq!(
            refused.status.code(),
            Some(1),
            "{stream}: {}",
            refused.status
        );
        assert!(refused.stdout.is_empty(), "{stream}");
    }
}

/// `--ontology` rewrites `Consequence is <term>` to `Consequence is_child
/// <term>` (115-118), which passes an entry whose term is the parent or lies
/// below it in the Sequence Ontology and fails every other; a parent the
/// table lacks is an error before the input is read.
#[test]
fn ontology_rewrites_consequence_is_to_is_child() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    let ids = |args: &[&str]| -> Vec<String> {
        let mut argv = vec!["-i", input, "--no_missing_consequence_warning"];
        argv.extend_from_slice(args);
        data_lines(&ok(&argv))
            .iter()
            .map(|l| l.split('\t').nth(2).unwrap().to_string())
            .collect()
    };
    assert_eq!(
        ids(&["-f", "Consequence is coding_sequence_variant"]),
        Vec::<String>::new()
    );
    assert_eq!(
        ids(&["--ontology", "-f", "Consequence is coding_sequence_variant"]),
        ["r1", "r2"]
    );
    assert_eq!(
        ids(&["-y", "-f", "consequence is UTR_variant"]),
        Vec::<String>::new()
    );
    // `stop_gained` lies below `feature_truncation` as well as below
    // `nonsynonymous_variant` in the Sequence Ontology.
    assert_eq!(
        ids(&["-y", "-f", "Consequence is feature_truncation"]),
        ["r1"]
    );
    assert_eq!(ids(&["-y", "-f", "Consequence is stop_gained"]), ["r1"]);
    assert_eq!(
        refused(&[
            "-i",
            "/nonexistent/in.vcf",
            "--ontology",
            "-f",
            "Consequence is coding_variant"
        ]),
        "ERROR: [filter] no SO term named coding_variant in the consequence ontology"
    );
}

/// The flags that name an Ensembl database connection are refused by name,
/// before the input is read, `--version` among them (91-97).
#[test]
fn database_flags_are_refused_by_name() {
    for (flag, value) in [
        ("--host", "ensembldb.ensembl.org"),
        ("--user", "anonymous"),
        ("--pass", "x"),
        ("--password", "x"),
        ("--port", "3306"),
        ("--version", "116"),
        ("--registry", "reg.pm"),
    ] {
        let message = refused(&[
            "-i",
            "/nonexistent/in.vcf",
            "-y",
            "-f",
            "Consequence is coding_sequence_variant",
            flag,
            value,
        ]);
        let named = if flag == "--password" { "--pass" } else { flag };
        assert!(
            message.starts_with(&format!("ERROR: [filter] {named} is not supported: ")),
            "{flag}: {message}"
        );
    }
}

/// A field written as the prefix of one key reads that key's value on every
/// line, the first included, where `filter_vep` reads nothing on the first
/// (`FilterSet.pm:570-571`); a prefix of several keys ends the run with the
/// error of `FilterSet.pm:574-580`, after the header block.
#[test]
fn prefix_synonym_reads_the_value_on_the_first_line() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    let out = ok(&[
        "-i",
        input,
        "--only_matched",
        "-f",
        "Symb match TP",
        "--no_missing_consequence_warning",
    ]);
    assert_eq!(data_lines(&out), [VCF_LINES.lines().next().unwrap()]);
    assert_eq!(
        ok(&[
            "-i",
            input,
            "-c",
            "-f",
            "Dist nex",
            "--no_missing_consequence_warning"
        ]),
        "3\n"
    );
    let output = filter(&["-i", input, "-f", "I is HIGH"]);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        stderr(&output),
        "ERROR: [filter] Multiple data fields found matching query field \"I\": ID, INFO, IMPACT\n"
    );
    assert!(stdout(&output).starts_with(VCF_HEADER));
}

/// The command line is read as Getopt::Long reads `filter_vep`'s: one dash
/// or two, any case, a unique prefix, `=` for a value, and the short forms
/// `-i`, `-o`, `-f`, `-s`, `-c`, `-l`, `-y`; no argument at all prints the
/// usage and exits 0, as does `--help`.
#[test]
fn getopt_spellings_reach_the_tool() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    assert_eq!(
        ok(&[
            "-input_file",
            input,
            "-filter",
            "IMPACT is HIGH",
            "--COUNT",
            "--no_missing"
        ]),
        "1\n"
    );
    assert_eq!(
        ok(&[
            &format!("--input_file={input}"),
            "--fil=IMPACT is HIGH",
            "-c",
            "--no_missing_consequence_warning"
        ]),
        "1\n"
    );
    let usage = filter(&[]);
    assert!(usage.status.success());
    assert!(
        stdout(&usage).contains("--only_matched"),
        "{}",
        stdout(&usage)
    );
    let help = filter(&["--help"]);
    assert!(help.status.success());
    assert!(stdout(&help).contains("--soft_filter"));
    let unknown = filter(&["--nope", "-i", input, "-f", "DP"]);
    assert_eq!(unknown.status.code(), Some(2));
}

/// A name that is a prefix of two options is refused after one dash as after
/// two, where `filter_vep` ends on `Option li is ambiguous (limit, list)`
/// (72-103). Its letters are not read as short forms: `-li` lists no field of
/// the file after it, `-fo` takes no filter and `-on` writes no file named
/// `n`.
#[test]
fn a_one_dash_prefix_of_two_options_is_refused_as_after_two_dashes() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    let refusal = |name: &str, similar: &str, usage: &str| {
        format!(
            "error: unexpected argument '--{name}' found\n\n  tip: a similar argument exists: \
             '--{similar}'\n\nUsage: vep filter {usage}\n\nFor more information, try '--help'.\n"
        )
    };
    for (name, rest, expected) in [
        ("li", vec![input], refusal("li", "list", "--list")),
        (
            "fo",
            vec!["x", "-i", input],
            refusal("fo", "format", "--format <vcf|tab>"),
        ),
        (
            "on",
            vec!["-i", input, "-f", "IMPACT is HIGH"],
            refusal("on", "ontology", "--ontology"),
        ),
    ] {
        for dashes in ["-", "--"] {
            let option = format!("{dashes}{name}");
            let mut args = vec![option.as_str()];
            args.extend(&rest);
            let output = filter_in(dir.path(), &args, None);
            assert_eq!(output.status.code(), Some(2), "{args:?}");
            assert_eq!(stdout(&output), "", "{args:?}");
            assert_eq!(stderr(&output), expected, "{args:?}");
        }
    }
    assert!(
        !dir.path().join("n").exists(),
        "`-on` wrote its output to a file named `n`"
    );
}

/// A repeated option runs with its last value, as `filter_vep`'s
/// `GetOptions` (72-103) assigns each occurrence of a scalar or switch over
/// the last: `-s 1 -s 2` writes what `-s 2` writes, `--count --count` what
/// `--count` writes, `--limit 2 --limit 1` what `--limit 1` writes. A
/// repeated `--filter` keeps every filter (102), joined by `and`.
#[test]
fn repeated_options_run_with_the_last_value() {
    let dir = tempfile::tempdir().unwrap();
    let input = vcf_fixture(dir.path());
    let input = input.to_str().unwrap();
    let passing = "Consequence match variant";
    let quiet = "--no_missing_consequence_warning";
    let second = VCF_LINES.lines().nth(1).unwrap();
    let out = ok(&["-i", input, "-s", "1", "-s", "2", "-f", passing, quiet]);
    assert_eq!(data_lines(&out), [second]);
    assert_eq!(out, ok(&["-i", input, "-s", "2", "-f", passing, quiet]));
    assert_eq!(
        ok(&["-i", input, "--count", "--count", "-f", passing, quiet]),
        "2\n"
    );
    let first = VCF_LINES.lines().next().unwrap();
    let out = ok(&[
        "-i", input, "--limit", "2", "--limit", "1", "-f", passing, quiet,
    ]);
    assert_eq!(data_lines(&out), [first]);
    let both = |symbol: &str| {
        ok(&[
            "-i",
            input,
            "--count",
            "-f",
            "IMPACT is HIGH",
            "-f",
            &format!("SYMBOL is {symbol}"),
            quiet,
        ])
    };
    assert_eq!(both("TPTE"), "1\n");
    assert_eq!(both("POTED"), "0\n");
}

/// A `.gz` input, or any gzip stream, is decompressed (153-155), and an input
/// with no `--input_file` is read from standard input (161-163).
#[test]
fn gzip_and_standard_input_are_read() {
    let dir = tempfile::tempdir().unwrap();
    let text = format!("{VCF_HEADER}{VCF_LINES}");
    let gz = dir.path().join("in.vcf.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        fs::File::create(&gz).unwrap(),
        flate2::Compression::default(),
    );
    encoder.write_all(text.as_bytes()).unwrap();
    encoder.finish().unwrap();
    let args = [
        "-f",
        "IMPACT is HIGH",
        "-c",
        "--no_missing_consequence_warning",
    ];
    let mut with_file = vec!["-i", gz.to_str().unwrap()];
    with_file.extend_from_slice(&args);
    assert_eq!(ok(&with_file), "1\n");
    let renamed = dir.path().join("in.vcf.compressed");
    fs::copy(&gz, &renamed).unwrap();
    let mut with_renamed = vec!["-i", renamed.to_str().unwrap()];
    with_renamed.extend_from_slice(&args);
    assert_eq!(ok(&with_renamed), "1\n");
    let output = filter_in(dir.path(), &args, Some(&text));
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(stdout(&output), "1\n");
}

/// The standard output bytes of a run that must succeed quietly.
fn ok_bytes(args: &[&str]) -> Vec<u8> {
    let output = filter(args);
    assert!(
        output.status.success(),
        "{args:?} exited {:?}: {}",
        output.status.code(),
        stderr(&output)
    );
    assert_eq!(stderr(&output), "", "{args:?}");
    output.stdout
}

fn concat(parts: &[&[u8]]) -> Vec<u8> {
    parts.concat()
}

/// `text`, every character of which lies in Latin-1, as its Latin-1 bytes:
/// the one byte each character is in an input written under that encoding.
fn latin1(text: &str) -> Vec<u8> {
    text.chars()
        .map(|c| u8::try_from(u32::from(c)).unwrap())
        .collect()
}

/// A passing record whose ID holds the Latin-1 `é`, the byte 0xE9, which is
/// not UTF-8 where a tab follows it.
const LATIN1_PASS: &[u8] = b"21\t100\tr\xE9\tA\tG\t50\t.\tDP=13;CSQ=G|stop_gained|HIGH|TPTE|\n";
/// A failing record of ASCII alone.
const LATIN1_FAIL: &[u8] =
    b"21\t200\tr2\tC\tT\t50\tPASS\tDP=20;CSQ=T|missense_variant|MODERATE|POTED|\n";
/// A record with two Latin-1 bytes in its ID and two entries, of which the
/// second passes `IMPACT is HIGH`.
const LATIN1_TWO_ENTRIES: &[u8] = b"21\t300\tr\xE9\xE8\tG\tA\t50\t.\tDP=4;CSQ=G|intron_variant|MODIFIER|TPTE|,G|stop_gained|HIGH|TPTE|\n";

/// A line holding a byte outside UTF-8 is filtered and written as the bytes
/// it was read from, as `filter_vep` does with no decoding layer on its input
/// (157-158, 184-185): a passing line is written as read (319), `--count`
/// counts it, `--soft_filter` marks it (316), `--only_matched` rewrites its
/// consequence field around the bytes (296), a header line's byte is written
/// through, a gzip input reads the same, `--list` prints a field name as read,
/// an error quoting the keys of the input quotes their bytes, and a line that
/// is UTF-8 throughout is unchanged. Every expected output is the one Ensembl
/// VEP's `filter_vep` writes on the same bytes, the error's keys in `vep
/// filter`'s own order.
#[test]
fn non_utf8_lines_are_filtered_and_written_as_read() {
    let dir = tempfile::tempdir().unwrap();
    let quiet = "--no_missing_consequence_warning";
    let header = VCF_HEADER.as_bytes();
    let fixture = |name: &str, bytes: &[u8]| -> String {
        let path = dir.path().join(name);
        fs::write(&path, bytes).unwrap();
        path.to_str().unwrap().to_string()
    };

    let data = fixture(
        "latin1-data.vcf",
        &concat(&[header, LATIN1_PASS, LATIN1_FAIL]),
    );
    let passing = concat(&[header, LATIN1_PASS]);
    assert_eq!(
        ok_bytes(&["-i", &data, "-f", "IMPACT is HIGH", quiet]),
        passing
    );
    assert_eq!(
        ok_bytes(&["-i", &data, "-f", "IMPACT is MODERATE", quiet]),
        concat(&[header, LATIN1_FAIL])
    );
    assert_eq!(
        ok_bytes(&["-i", &data, "-c", "-f", "IMPACT is HIGH", quiet]),
        b"1\n"
    );
    let out_path = dir.path().join("out.vcf");
    let out = out_path.to_str().unwrap();
    assert_eq!(
        ok_bytes(&["-i", &data, "-o", out, "-f", "IMPACT is HIGH", quiet]),
        b""
    );
    assert_eq!(fs::read(&out_path).unwrap(), passing);

    let (meta, column_line) = VCF_HEADER.split_at(VCF_HEADER.find("#CHROM").unwrap());
    assert_eq!(
        ok_bytes(&["-i", &data, "--soft_filter", "-f", "IMPACT is HIGH", quiet]),
        concat(&[
            meta.as_bytes(),
            b"##FILTER=<ID=filter_vep_pass,Description=\"Variant passes filter_vep\">\n",
            b"##FILTER=<ID=filter_vep_fail,Description=\"Variant fails filter_vep\">\n",
            column_line.as_bytes(),
            b"21\t100\tr\xE9\tA\tG\t50\tfilter_vep_pass\tDP=13;CSQ=G|stop_gained|HIGH|TPTE|\n",
            b"21\t200\tr2\tC\tT\t50\tPASS;filter_vep_fail\tDP=20;CSQ=T|missense_variant|MODERATE|POTED|\n",
        ])
    );

    let rewrite = fixture(
        "latin1-rewrite.vcf",
        &concat(&[header, LATIN1_PASS, LATIN1_FAIL, LATIN1_TWO_ENTRIES]),
    );
    assert_eq!(
        ok_bytes(&[
            "-i",
            &rewrite,
            "--only_matched",
            "-f",
            "IMPACT is HIGH",
            quiet
        ]),
        concat(&[
            header,
            LATIN1_PASS,
            b"21\t300\tr\xE9\xE8\tG\tA\t50\t.\tDP=4;CSQ=G|stop_gained|HIGH|TPTE|\n",
        ])
    );

    let described = latin1(&VCF_HEADER.replace(
        "Description=\"Depth\"",
        "Description=\"Profondeur de s\u{e9}quen\u{e7}age\"",
    ));
    let ascii_pass = b"21\t100\tr1\tA\tG\t50\t.\tDP=13;CSQ=G|stop_gained|HIGH|TPTE|\n";
    let in_header = fixture(
        "latin1-header.vcf",
        &concat(&[&described, ascii_pass, LATIN1_FAIL]),
    );
    assert_eq!(
        ok_bytes(&["-i", &in_header, "-f", "IMPACT is HIGH", quiet]),
        concat(&[&described, ascii_pass])
    );

    let gz = dir.path().join("latin1-data.vcf.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        fs::File::create(&gz).unwrap(),
        flate2::Compression::default(),
    );
    encoder.write_all(&fs::read(&data).unwrap()).unwrap();
    encoder.finish().unwrap();
    assert_eq!(
        ok_bytes(&["-i", gz.to_str().unwrap(), "-f", "IMPACT is HIGH", quiet]),
        passing
    );

    let listed = fixture(
        "latin1-list.vcf",
        &concat(&[
            &latin1(&VCF_HEADER.replace("ID=DP,", "ID=D\u{e9}P,")),
            ascii_pass,
            LATIN1_FAIL,
        ]),
    );
    assert_eq!(
        ok_bytes(&["-i", &listed, "--list"]),
        b"Available fields:\n\nALT\nAllele\nCHROM\nConsequence\nDISTANCE\nD\xE9P\nFILTER\nID\nIMPACT\nINFO\nPOS\nQUAL\nREF\nSYMBOL\n"
    );
    let keys = fixture(
        "latin1-keys.vcf",
        &concat(&[
            &latin1(&VCF_HEADER.replace("ID=DP,", "ID=D\u{e9}P,")),
            &latin1(
                "21\t100\tr1\tA\tG\t50\t.\tD\u{e9}P=13;D\u{e9}Q=2;CSQ=G|stop_gained|HIGH|TPTE|\n",
            ),
        ]),
    );
    let output = filter(&["-i", &keys, "-f", "D is 13", quiet]);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        output.stderr,
        b"ERROR: [filter] Multiple data fields found matching query field \"D\": D\xE9P, D\xE9Q, DISTANCE\n"
    );

    let utf8_pass = "21\t100\tr\u{e9}\tA\tG\t50\t.\tDP=13;CSQ=G|stop_gained|HIGH|TPTE|\n";
    let utf8 = fixture(
        "utf8-data.vcf",
        &concat(&[header, utf8_pass.as_bytes(), LATIN1_FAIL]),
    );
    assert_eq!(
        ok_bytes(&["-i", &utf8, "-f", "IMPACT is HIGH", quiet]),
        concat(&[header, utf8_pass.as_bytes()])
    );
}

/// Writes `bytes` under `name` in `dir` and returns the path as text.
fn byte_fixture(dir: &Path, name: &str, bytes: &[u8]) -> String {
    let path = dir.join(name);
    fs::write(&path, bytes).unwrap();
    path.to_str().unwrap().to_string()
}

/// One record per ID of `ids`, at positions 100, 200, ..., each with one
/// passing entry.
fn records(ids: &[&[u8]]) -> Vec<u8> {
    let mut out = Vec::new();
    for (i, id) in ids.iter().enumerate() {
        out.extend_from_slice(format!("21\t{}\t", (i + 1) * 100).as_bytes());
        out.extend_from_slice(id);
        out.extend_from_slice(b"\tA\tG\t50\t.\tDP=13;CSQ=G|stop_gained|HIGH|TPTE|\n");
    }
    out
}

/// The positions of the data lines of a run's standard output bytes.
fn positions(output: &[u8]) -> Vec<u32> {
    output
        .split(|b| *b == b'\n')
        .filter(|line| !line.is_empty() && line[0] != b'#')
        .map(|line| {
            let pos = line.split(|b| *b == b'\t').nth(1).unwrap();
            std::str::from_utf8(pos).unwrap().parse().unwrap()
        })
        .collect()
}

/// Six IDs: `é` as UTF-8, `é` as Latin-1, `É` as UTF-8, `E`, `e`, and a
/// no-break space as UTF-8, at positions 100 to 600.
const CASED_IDS: &[&[u8]] = &[
    b"r\xC3\xA9",
    b"r\xE9",
    b"r\xC3\x89",
    b"rE",
    b"re",
    b"r\xC2\xA0",
];

/// `match` runs its pattern over the bytes of the value, as Perl's regex
/// engine runs `/$pattern/i` over a byte string (`FilterSet.pm:714`): `.` is
/// one byte, so `é` written as UTF-8 is two of them; a character class holds
/// bytes, so `[é]` is the class of its two bytes and `[\xC3][\xA9]` matches
/// the UTF-8 `é`; the case fold is ASCII alone, so `É` matches neither `é`;
/// `\w` and `\s` are ASCII, so no byte of a UTF-8 letter or a no-break space
/// is one; a quantifier binds the last byte of a multibyte character. Every
/// expected set of positions is the one `filter_vep` writes on the same
/// input.
#[test]
fn match_reads_the_bytes_of_the_value() {
    let dir = tempfile::tempdir().unwrap();
    let quiet = "--no_missing_consequence_warning";
    let input = byte_fixture(
        dir.path(),
        "cased.vcf",
        &concat(&[VCF_HEADER.as_bytes(), &records(CASED_IDS)]),
    );
    let passing =
        |filter: &str| -> Vec<u32> { positions(&ok_bytes(&["-i", &input, "-f", filter, quiet])) };
    assert_eq!(passing("ID match r.$"), [200, 400, 500]);
    assert_eq!(passing("ID match r.{2}$"), [100, 300, 600]);
    assert_eq!(passing("ID match r[\u{e9}]$"), [0u32; 0]);
    assert_eq!(passing(r"ID match r[\xC3][\xA9]"), [100]);
    assert_eq!(passing("ID match \u{c9}$"), [300]);
    assert_eq!(passing("ID match \u{e9}$"), [100]);
    assert_eq!(passing("ID match E$"), [400, 500]);
    assert_eq!(passing(r"ID match r\w$"), [400, 500]);
    assert_eq!(passing(r"ID match r\s$"), [0u32; 0]);
    assert_eq!(passing("ID match r\u{e9}+$"), [100]);
    assert_eq!(passing("ID nre r.$"), [100, 300, 600]);
}

/// An `in <file>` list is read as bytes, as Ensembl VEP reads the file with no
/// decoding layer (`FilterSet.pm:397-409`): a Latin-1 list matches the
/// Latin-1 value and a UTF-8 list the UTF-8 value, byte for byte. Each line
/// loses its carriage returns and every byte Perl's `\v` matches on a byte
/// string, the bytes 0x0B, 0x0C and 0x85 among them, so a UTF-8 next line
/// (`C2 85`) leaves its lead byte behind and a UTF-8 line separator (U+2028,
/// three bytes) stays whole. Every expected set is `filter_vep`'s on the same
/// files.
#[test]
fn in_list_file_is_read_as_bytes() {
    let dir = tempfile::tempdir().unwrap();
    let quiet = "--no_missing_consequence_warning";
    let header = VCF_HEADER.as_bytes();
    let cased = byte_fixture(
        dir.path(),
        "cased.vcf",
        &concat(&[header, &records(CASED_IDS)]),
    );
    let latin1 = byte_fixture(dir.path(), "latin1.txt", b"r\xE9\n");
    let utf8 = byte_fixture(dir.path(), "utf8.txt", b"r\xC3\xA9\n");
    let passing = |input: &str, filter: &str| -> Vec<u32> {
        positions(&ok_bytes(&["-i", input, "-f", filter, quiet]))
    };
    assert_eq!(passing(&cased, &format!("ID in {latin1}")), [200]);
    assert_eq!(passing(&cased, &format!("ID in {utf8}")), [100]);
    assert_eq!(
        passing(&cased, &format!("not ID in {latin1}")),
        [100, 300, 400, 500, 600]
    );

    let vertical = byte_fixture(
        dir.path(),
        "vertical.txt",
        b"a\x85\nb\xC2\x85\nc\xE2\x80\xA8\nd\x0Be\nf\x0Cg\nh\r\n",
    );
    let spaced = byte_fixture(
        dir.path(),
        "spaced.vcf",
        &concat(&[
            header,
            &records(&[
                b"a",
                b"a\x85",
                b"b",
                b"b\xC2",
                b"b\xC2\x85",
                b"c",
                b"c\xE2\x80\xA8",
                b"de",
                b"d\x0Be",
                b"fg",
                b"h",
            ]),
        ]),
    );
    assert_eq!(
        passing(&spaced, &format!("ID in {vertical}")),
        [100, 400, 700, 800, 1000, 1100]
    );
}

/// Whether a file system refused the creation of a file over its name, which
/// is outside UTF-8: the creation failed with `EILSEQ` or `EINVAL`, the two
/// errors a file system returns over the encoding of a name. A creation that
/// failed for any other reason panics.
#[cfg(unix)]
fn refused_for_its_name(created: std::io::Result<()>) -> bool {
    match created {
        Ok(()) => false,
        Err(error) if matches!(error.raw_os_error(), Some(libc::EILSEQ | libc::EINVAL)) => true,
        Err(error) => panic!("a file named outside UTF-8 was not made: {error}"),
    }
}

/// A creation that succeeded is no refusal, and each of the two errors over a
/// name's encoding is one.
#[cfg(unix)]
#[test]
fn a_refusal_of_a_file_name_is_told_from_a_creation() {
    use std::io::Error;
    assert!(!refused_for_its_name(Ok(())));
    for refusal in [libc::EILSEQ, libc::EINVAL] {
        assert!(refused_for_its_name(Err(Error::from_raw_os_error(refusal))));
    }
}

/// A creation that fails for any reason but the name fails the test.
#[cfg(unix)]
#[test]
#[should_panic(expected = "a file named outside UTF-8 was not made")]
fn a_creation_failure_other_than_the_name_fails_the_test() {
    refused_for_its_name(Err(std::io::Error::from_raw_os_error(libc::EACCES)));
}

/// Every command-line word is read as the bytes the shell passed, as Ensembl VEP
/// reads `@ARGV`: a filter value outside UTF-8 compares byte for byte under
/// `is`, `match` and a character class, with no case fold above ASCII, and an
/// input or output file whose name is outside UTF-8 is opened, where the file
/// system takes such a name. Each expected output is `filter_vep`'s on the
/// same bytes.
#[cfg(unix)]
#[test]
fn command_line_words_are_read_as_bytes() {
    use std::ffi::OsStr;
    use std::os::unix::ffi::OsStrExt;
    let dir = tempfile::tempdir().unwrap();
    let header = VCF_HEADER.as_bytes();
    let cased = concat(&[header, &records(CASED_IDS)]);
    let input = dir.path().join("cased.vcf");
    fs::write(&input, &cased).unwrap();
    let run = |args: &[&[u8]]| -> Output {
        let mut cmd = common::isolated_vep_anywhere();
        cmd.arg("filter")
            .arg("-i")
            .arg(&input)
            .arg("--no_missing_consequence_warning")
            .current_dir(dir.path());
        for arg in args {
            cmd.arg(OsStr::from_bytes(arg));
        }
        cmd.stdin(Stdio::null()).output().unwrap()
    };
    let passing = |filter: &[u8]| -> Vec<u32> {
        let output = run(&[b"-f", filter]);
        assert!(
            output.status.success() && output.stderr.is_empty(),
            "{filter:?}: exit {:?}, stderr {}",
            output.status.code(),
            stderr(&output)
        );
        positions(&output.stdout)
    };
    assert_eq!(passing(b"ID is r\xE9"), [200]);
    assert_eq!(passing(b"ID match r\xE9"), [200]);
    assert_eq!(passing(b"ID match r[\xE9]$"), [200]);
    assert_eq!(passing(b"ID match r\xC9"), [0u32; 0]);
    assert_eq!(passing(b"ID nre \xE9"), [100, 300, 400, 500, 600]);
    let output = run(&[b"-f", b"ID in list-\xE9.txt"]);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        output.stderr,
        b"ERROR: [filter] Could not find/parse list list-\xE9.txt\n"
    );

    let named = dir.path().join(OsStr::from_bytes(b"in-\xE9.vcf"));
    let created = fs::write(&named, concat(&[header, &records(&[b"r\xE9"])]));
    if refused_for_its_name(created) {
        eprintln!(
            "the scratch directory does not take a file name outside UTF-8: \
             the runs on such a file are left out"
        );
        return;
    }
    let out = dir.path().join(OsStr::from_bytes(b"out-\xE9.vcf"));
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("filter")
        .arg("-i")
        .arg(&named)
        .arg("-o")
        .arg(&out)
        .arg("-f")
        .arg("IMPACT is HIGH")
        .arg("--no_missing_consequence_warning")
        .current_dir(dir.path());
    let output = cmd.stdin(Stdio::null()).output().unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(output.stdout, b"");
    assert_eq!(
        fs::read(&out).unwrap(),
        concat(&[header, &records(&[b"r\xE9"])])
    );
}

/// `--list` prints the field names in the order of their bytes, as Perl's
/// `sort` orders byte strings (`filter_vep:221`): a raw byte 0x80 sorts before
/// the UTF-8 `é` (`C3 A9`), which sorts before a raw 0xE9, which sorts before
/// a four-byte UTF-8 character (`F0 ...`). The expected output is
/// `filter_vep`'s on the same header block.
#[test]
fn list_orders_the_field_names_by_byte() {
    let dir = tempfile::tempdir().unwrap();
    let info = |name: &[u8]| -> Vec<u8> {
        concat(&[
            b"##INFO=<ID=",
            name,
            b",Number=1,Type=String,Description=\"x\">\n",
        ])
    };
    let input = byte_fixture(
        dir.path(),
        "order.vcf",
        &concat(&[
            b"##fileformat=VCFv4.2\n",
            &info(b"D\x80"),
            &info(b"D\xC3\xA9"),
            &info(b"D\xE9"),
            &info(b"D\xF0\x9F\x98\x80"),
            &info(b"DZ"),
            &info(b"D"),
            b"##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence\">\n",
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            b"21\t100\tr1\tA\tG\t50\t.\tDZ=1;CSQ=G|stop_gained\n",
        ]),
    );
    assert_eq!(
        ok_bytes(&["-i", &input, "--list"]),
        concat(&[
            b"Available fields:\n\nALT\nAllele\nCHROM\nConsequence\nD\nDZ\nD\x80\nD\xC3\xA9\nD\xE9\nD\xF0\x9F\x98\x80\n",
            b"FILTER\nID\nINFO\nPOS\nQUAL\nREF\n",
        ])
    );
}

/// Every text rule of `vep filter` reads bytes as Perl's rules read a byte
/// string. The numeric extraction's `\w` is ASCII, so `é(0.5)` yields no
/// number, and a no-break space before `5` is no blank, so the value reads as
/// 0; `Inf`, `Infinity` and `NaN` are the numbers Perl reads them as. The lower- and
/// upper-case synonyms of a field name change its ASCII letters alone, and the
/// prefix search folds ASCII case alone, so `ék` names neither `ÉK` nor `ek`
/// and `é` is a prefix of neither. The filter string splits on ASCII blanks
/// alone, so a no-break space joins `ID` and `is` into one field name. A
/// first column of one UTF-8 letter has no word character, so the line is a
/// tab line. Every expected set is `filter_vep`'s on the same input.
#[test]
fn text_rules_read_bytes_as_perl_does() {
    let dir = tempfile::tempdir().unwrap();
    let quiet = "--no_missing_consequence_warning";
    let info = |name: &[u8]| -> Vec<u8> {
        concat(&[
            b"##INFO=<ID=",
            name,
            b",Number=1,Type=String,Description=\"x\">\n",
        ])
    };
    let csq = b"##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence\">\n";
    let column_line = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let record = |pos: u32, info: &[u8]| -> Vec<u8> {
        concat(&[
            format!("21\t{pos}\tr\tA\tG\t50\t.\t").as_bytes(),
            info,
            b";CSQ=G|stop_gained\n",
        ])
    };
    let values = byte_fixture(
        dir.path(),
        "values.vcf",
        &concat(&[
            b"##fileformat=VCFv4.2\n",
            &info(b"S"),
            &info(b"\xC3\x89K"),
            &info(b"ek"),
            csq,
            column_line,
            &record(100, b"S=\xC3\xA9(0.5);\xC3\x89K=1;ek=2"),
            &record(200, b"S=x(0.5);\xC3\x89K=1;ek=2"),
            &record(300, b"S=\xC2\xA05"),
            &record(400, b"S=\xD9\xA3"),
            &record(500, b"S=Inf"),
            &record(600, b"S=-inf"),
            &record(700, b"S=NaN"),
            &record(800, b"S=Infinity"),
        ]),
    );
    let passing = |input: &str, filter: &str| -> Vec<u32> {
        positions(&ok_bytes(&["-i", input, "-f", filter, quiet]))
    };
    assert_eq!(passing(&values, "S is 0.5"), [200]);
    assert_eq!(passing(&values, "S > 1"), [500, 800]);
    assert_eq!(passing(&values, "S < 1"), [100, 200, 300, 400, 600]);
    assert_eq!(passing(&values, "\u{e9}k is 1"), [0u32; 0]);
    assert_eq!(passing(&values, "\u{c9}k is 1"), [100, 200]);
    assert_eq!(passing(&values, "EK is 2"), [100, 200]);
    assert_eq!(passing(&values, "\u{e9} is 1"), [0u32; 0]);

    let cased = byte_fixture(
        dir.path(),
        "cased.vcf",
        &concat(&[VCF_HEADER.as_bytes(), &records(CASED_IDS)]),
    );
    assert_eq!(
        refused(&["-i", &cased, "-f", "ID\u{a0}is rE", quiet]),
        "ERROR: [filter] No such operator \"rE\""
    );

    let letter = byte_fixture(
        dir.path(),
        "letter.vcf",
        &concat(&[
            b"##fileformat=VCFv4.2\n",
            csq,
            column_line,
            b"\xC3\xA9\t100\tr1\tA\tG\t50\t.\tCSQ=G|stop_gained\n",
        ]),
    );
    assert_eq!(passing(&letter, "Consequence is stop_gained"), [0u32; 0]);
    assert_eq!(passing(&letter, "Consequence is 100"), [100]);
}

/// The 42 IDs of the case-folding runs: the fold partners of `ss`, `fi`,
/// `st`, `ff`, `ffi`, `ffl`, `fl`, `s`, `k`, `é` and `µ` in UTF-8 and in
/// Latin-1, with controls.
const FOLD_IDS: &[&[u8]] = &[
    br"ss",
    br"SS",
    br"sS",
    b"\xC3\x9F",
    b"\xDF",
    b"\xEF\xAC\x81",
    br"fi",
    br"FI",
    b"\xC5\xBF",
    br"s",
    br"S",
    b"\xE2\x84\xAA",
    br"k",
    br"K",
    b"\xE1\xBA\x9E",
    br"st",
    b"\xEF\xAC\x86",
    br"ass",
    b"a\xDFa",
    br"assa",
    b"\xB5",
    br"ff",
    b"\xEF\xAC\x80",
    br"ffi",
    b"\xEF\xAC\x83",
    b"a\xDF",
    b"\xDFs",
    br"sss",
    b"s\xDF",
    br"ssss",
    br"xss",
    b"x\xDF",
    br"xx",
    b"\xDF\xDF",
    b"\xC3\x9F\xC3\x9F",
    b"\xE9",
    b"\xC9",
    br"i",
    br"I",
    br"t",
    br"fl",
    br"ffl",
];

/// The 30 IDs of the brace and octal runs: the values `r{`, `r{x}`, `r{,3}`,
/// `rA`, `r` with the bytes 0x00, 0x01, 0x08, 0x10, 0xE9 and U+0100, and the
/// brace words.
const ESC_IDS: &[&[u8]] = &[
    br"r{",
    br"r{x}",
    br"r{,3}",
    br"a",
    br"aaa",
    br"aaaa",
    br"A",
    b"r\x00",
    br"r0",
    br"rA",
    br"r{2}",
    br"rr",
    br"r{1,2}",
    br"r{,}",
    br"r{1,",
    b"r\xE9",
    b"r\x01",
    br"{",
    br"r}",
    br"r",
    br"rrr",
    b"r\x081",
    b"r\x08",
    br"ra",
    br"r{1}",
    br"{x}",
    b"r\xC4\x80",
    b"r\x100",
    br"rA0",
    b"r\x001",
];

/// A fixture of one record per ID at positions 100, 200, ..., with a header
/// declaring `DP`, `PAT` and the consequence field, and a runner over it whose
/// filter word is bytes.
fn byte_runs<'a>(dir: &'a Path, ids: &[&[u8]]) -> impl Fn(&[u8]) -> Output + 'a {
    use std::ffi::OsStr;
    use std::os::unix::ffi::OsStrExt;
    let header = "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##INFO=<ID=PAT,Number=1,Type=String,Description=\"Pattern\">\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|DISTANCE\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let input = dir.join("shapes.vcf");
    fs::write(&input, concat(&[header.as_bytes(), &records(ids)])).unwrap();
    move |filter: &[u8]| -> Output {
        let mut cmd = common::isolated_vep_anywhere();
        cmd.arg("filter")
            .arg("-i")
            .arg(&input)
            .arg("--no_missing_consequence_warning")
            .arg("-f")
            .arg(OsStr::from_bytes(filter))
            .current_dir(dir);
        cmd.stdin(Stdio::null()).output().unwrap()
    }
}

/// The positions kept by a run of `filter` that must succeed quietly.
fn kept(run: &dyn Fn(&[u8]) -> Output, filter: &[u8]) -> Vec<u32> {
    let output = run(filter);
    assert!(
        output.status.success() && output.stderr.is_empty(),
        "{}: exit {:?}, stderr {}",
        String::from_utf8_lossy(filter),
        output.status.code(),
        stderr(&output)
    );
    positions(&output.stdout)
}

/// Asserts a failing run of `filter` whose error carries `message`.
fn refused_with(run: &dyn Fn(&[u8]) -> Output, filter: &[u8], message: &str) {
    let output = run(filter);
    assert!(
        !output.status.success() && stderr(&output).contains(message),
        "{}: exit {:?}, stderr {}",
        String::from_utf8_lossy(filter),
        output.status.code(),
        stderr(&output)
    );
}

/// `match` folds case as Perl's `/i` does over the same bytes
/// (`FilterSet.pm:714`). Under the default rules only ASCII folds, so `ss`
/// matches neither `ß` nor `ẞ`. A `\x{}` above 255 or a `\p{}` class switches
/// the pattern, and a `u` or `a` flag its own group, to Unicode rules over the
/// bytes read as Latin-1, under which folds are full: `ss` and `\x{DF}` match each other and
/// `ß`, `\x{1E9E}` folds to `ss`, the ligatures `\x{FB00}` to `\x{FB06}` fold
/// to `ff`, `fi`, `fl`, `ffi`, `ffl` and `st`, the long s `\x{17F}` and the
/// Kelvin sign `\x{212A}` fold to `s` and `k`, and the fold of a run of
/// literals is compared whole, so `s\x{DF}` matches `sss`, `ßs` and `sß`.
/// A quantifier, a quantified class and a negated class break or exclude the
/// fold as Perl's do; `aa` forbids a fold across ASCII; `(?-i)` and `(?^u)`
/// turn the fold off. A pattern whose parenthesis the filter string reads as
/// grouping is refused as `filter_vep` refuses it. The IDs are the 42 values
/// of `FOLD_IDS` and every expected set of positions is the one
/// `filter_vep` writes on the same input.
#[test]
fn match_folds_case_fully_under_unicode_rules() {
    let dir = tempfile::tempdir().unwrap();
    let run = byte_runs(dir.path(), FOLD_IDS);
    let passing = |filter: &[u8]| kept(&run, filter);
    let refused = |filter: &[u8], message: &str| refused_with(&run, filter, message);
    assert_eq!(
        passing(br"ID match ^ss$"),
        [100, 200, 300, 1800, 2800, 3000, 3100]
    );
    assert_eq!(passing(br"ID match ^\xDF$"), [500, 2600, 2900, 3200, 3400]);
    assert_eq!(
        passing(br"ID match ^\x{DF}$"),
        [500, 2600, 2900, 3200, 3400]
    );
    assert_eq!(passing(b"ID match ^\xC3\x9F$"), [400, 3500]);
    assert_eq!(passing(br"ID match ^fi$"), [700, 800, 2400]);
    assert_eq!(
        passing(br"ID match ^s$"),
        [100, 200, 300, 1000, 1100, 1800, 2700, 2800, 3000, 3100]
    );
    assert_eq!(passing(br"ID match ^k$"), [1300, 1400]);
    assert_eq!(passing(br"ID match ^\xE9$"), [3600]);
    assert_eq!(passing(b"ID match ^\xDF$"), [500, 2600, 2900, 3200, 3400]);
    assert_eq!(
        passing(br"ID match ^ss\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(
        passing(br"ID match ^\xDF\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(
        passing(br"ID match ^\x{DF}\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(
        passing(b"ID match ^\xDF\\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(passing(b"ID match ^\xC3\x9F\\x{100}?$"), [400, 3500]);
    assert_eq!(passing(br"ID match ^\x{FB01}\x{100}?$"), [700, 800, 2400]);
    assert_eq!(passing(br"ID match ^fi\x{100}?$"), [700, 800, 2400]);
    assert_eq!(
        passing(br"ID match ^\x{17F}\x{100}?$"),
        [100, 200, 300, 1000, 1100, 1800, 2700, 2800, 3000, 3100]
    );
    assert_eq!(
        passing(br"ID match ^s\x{100}?$"),
        [100, 200, 300, 1000, 1100, 1800, 2700, 2800, 3000, 3100]
    );
    assert_eq!(passing(br"ID match ^\x{212A}\x{100}?$"), [1300, 1400]);
    assert_eq!(passing(br"ID match ^k\x{100}?$"), [1300, 1400]);
    assert_eq!(
        passing(br"ID match ^\x{1E9E}\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(passing(br"ID match ^st\x{100}?$"), [1600]);
    assert_eq!(passing(br"ID match ^\x{FB06}\x{100}?$"), [1600]);
    assert_eq!(passing(br"ID match ^\x{FB05}\x{100}?$"), [1600]);
    assert_eq!(passing(br"ID match ^\x{FB00}\x{100}?$"), [2200]);
    assert_eq!(passing(br"ID match ^\x{FB03}\x{100}?$"), [2400]);
    assert_eq!(passing(br"ID match ^\x{FB02}\x{100}?$"), [4100, 4200]);
    assert_eq!(passing(br"ID match ^\x{FB04}\x{100}?$"), [4200]);
    assert_eq!(
        passing(br"ID match ^s+\x{100}?$"),
        [100, 200, 300, 1000, 1100, 1800, 2700, 2800, 3000, 3100]
    );
    assert_eq!(
        passing(br"ID match ^s{2}\x{100}?$"),
        [100, 200, 300, 1800, 2800, 3000, 3100]
    );
    assert_eq!(
        passing(br"ID match ^[\x{DF}]\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(passing(br"ID match ^[^\x{DF}]\x{100}?$"), [0u32; 0]);
    assert_eq!(
        passing(br"ID match ^[s]s\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(passing(br"ID match ^(?:s)s\x{100}?$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^(?:ss)\x{100}?$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^assa\x{100}?$"), [1900, 2000]);
    assert_eq!(passing(br"ID match ^a\x{DF}a\x{100}?$"), [1900, 2000]);
    assert_eq!(
        passing(br"ID match ^s\x{DF}\x{100}?$"),
        [2700, 2800, 2900, 3000]
    );
    assert_eq!(
        passing(br"ID match ^sss\x{100}?$"),
        [2700, 2800, 2900, 3000]
    );
    assert_eq!(
        passing(br"ID match ^\x{DF}s\x{100}?$"),
        [2700, 2800, 2900, 3000]
    );
    assert_eq!(passing(br"ID match ^ssss\x{100}?$"), [3000, 3400]);
    assert_eq!(passing(br"ID match ^\x{DF}\x{DF}\x{100}?$"), [3000, 3400]);
    assert_eq!(
        passing(br"ID match ^\x{DF}+\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(
        passing(br"ID match ^[\x{DF}x]{2}\x{100}?$"),
        [3000, 3100, 3200, 3300, 3400]
    );
    assert_eq!(
        passing(br"ID match ^[\x{DF}x]\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3300, 3400]
    );
    assert_eq!(passing(br"ID match ^x\x{DF}?\x{100}?$"), [3100, 3200, 3300]);
    assert_eq!(passing(br"ID match ^ass+\x{100}?$"), [1800]);
    assert_eq!(passing(br"ID match ^as{2}\x{100}?$"), [1800]);
    assert_eq!(passing(br"ID match ^(?:ss|zz)\x{100}?$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^\xE9\x{100}?$"), [3600, 3700]);
    assert_eq!(passing(br"ID match ^\x{3BC}\x{100}?$"), [2100]);
    assert_eq!(passing(br"ID match ^\xB5\x{100}?$"), [2100]);
    assert_eq!(passing(br"ID match ^\x{39C}\x{100}?$"), [2100]);
    assert_eq!(
        passing(br"ID match ^i\x{100}?$"),
        [700, 800, 2400, 3800, 3900]
    );
    assert_eq!(passing(br"ID match ^\x{130}\x{100}?$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^\x{131}\x{100}?$"), [0u32; 0]);
    assert_eq!(
        passing(br"ID match ^\x{17F}s\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(
        passing(br"ID match ^\x{17F}\x{17F}\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(
        passing(br"ID match ^S\x{DF}\x{100}?$"),
        [2700, 2800, 2900, 3000]
    );
    assert_eq!(
        passing(br"ID match ^\x{DF}.\x{100}?$"),
        [1900, 2000, 2700, 2800, 3000, 3400]
    );
    assert_eq!(
        passing(br"ID match ^.\x{DF}\x{100}?$"),
        [1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(
        passing(br"ID match ^s\w\x{100}?$"),
        [100, 200, 300, 1600, 1800, 2000, 2800, 2900, 3000, 3100]
    );
    assert_eq!(passing(br"ID match ^fl\x{100}?$"), [4100, 4200]);
    assert_eq!(passing(br"ID match ^ffl\x{100}?$"), [4200]);
    assert_eq!(
        passing(br"ID match ^[\x{DF}s]\x{100}?$"),
        [100, 200, 300, 500, 1000, 1100, 1800, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(
        passing(br"ID match ^[\x{DF}s]+\x{100}?$"),
        [100, 200, 300, 500, 1000, 1100, 1800, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    assert_eq!(
        passing(br"ID match ^[s]{2}\x{100}?$"),
        [100, 200, 300, 1800, 2800, 3000, 3100]
    );
    assert_eq!(
        passing(br"ID match ^\w\w\x{100}?$"),
        [
            100, 200, 300, 700, 800, 1600, 1800, 1900, 2000, 2200, 2400, 2600, 2700, 2800, 2900,
            3000, 3100, 3200, 3300, 3400, 4100, 4200
        ]
    );
    assert_eq!(passing(br"ID match ^(?:\x{DF}|zz)\x{100}?$"), [0u32; 0]);
    assert_eq!(
        passing(br"ID match ^\x{DF}z?\x{100}?$"),
        [100, 200, 300, 500, 1800, 2600, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    refused(
        br"ID match (?u)^ss$",
        "Multiple data fields found matching query field",
    );
    assert_eq!(
        passing(br"ID match ^ss\p{L}?$"),
        [100, 200, 300, 500, 1800, 1900, 2000, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3400]
    );
    refused(
        br"ID match (?a)^ss$",
        "Multiple data fields found matching query field",
    );
    refused(
        br"ID match (?a)^\xE9$",
        "Multiple data fields found matching query field",
    );
    refused(
        br"ID match (?d)^ss\x{100}?$",
        "Multiple data fields found matching query field",
    );
    refused(
        br"ID match (?l)^ss$",
        "Multiple data fields found matching query field",
    );
    assert_eq!(passing(br"ID match (?x)^s s$"), [0u32; 0]);
    assert_eq!(passing(br"ID match (?-i)^ss$"), [0u32; 0]);
    assert_eq!(passing(br"ID match (?u:^ss$)"), [0u32; 0]);
    assert_eq!(passing(br"ID match (?^u)^ss$"), [0u32; 0]);
    assert_eq!(passing(br"ID match (?iu)^ss$"), [0u32; 0]);
    refused(
        br"ID match (?n)^(ss)$",
        "Multiple data fields found matching query field",
    );
}

/// `match` reads braces and octal escapes as Perl does over the same bytes
/// (`FilterSet.pm:714`): `{` opens a quantifier only as `{n}`, `{n,}`,
/// `{n,m}` or `{,m}` after an atom and is a literal elsewhere, `{n,m}` with
/// `n` above `m` matches nothing, a leading zero is refused; `\o{}` and `\0`
/// to `\777` are characters by value, blanks and single underscores inside the
/// braces skipped and the first other character ending the value; `\1` to
/// `\9` are back-references, refused when no such group exists; `\10` and
/// above is octal when fewer groups exist; `\o{400}` and `\400` switch the
/// pattern to Unicode rules; `\N{U+}` reads hexadecimal numbers with single
/// underscores between digits and `.`-separated sequences, refusing a
/// malformed number with Perl's text and a sequence before a `-`. The
/// IDs are the 30 values of `ESC_IDS` and every expected set of
/// positions is the one `filter_vep` writes on the same input.
#[test]
fn match_reads_braces_and_octal_escapes_as_perl_does() {
    let dir = tempfile::tempdir().unwrap();
    let run = byte_runs(dir.path(), ESC_IDS);
    let passing = |filter: &[u8]| kept(&run, filter);
    let refused = |filter: &[u8], message: &str| refused_with(&run, filter, message);
    assert_eq!(passing(br"ID match ^r{$"), [100]);
    assert_eq!(passing(br"ID match ^r{x}$"), [200]);
    assert_eq!(
        passing(br"ID match ^r{,3}$"),
        [
            100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600,
            1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000
        ]
    );
    assert_eq!(
        passing(br"ID match ^a{,3}$"),
        [
            100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600,
            1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000
        ]
    );
    assert_eq!(passing(br"ID match ^r\o{101}$"), [1000, 2400]);
    assert_eq!(passing(br"ID match ^r\101$"), [1000, 2400]);
    assert_eq!(passing(br"ID match ^r\0$"), [800]);
    assert_eq!(passing(br"ID match ^r\00$"), [800]);
    assert_eq!(passing(br"ID match ^r\000$"), [800]);
    assert_eq!(passing(br"ID match ^r\0001$"), [3000]);
    assert_eq!(passing(br"ID match ^r\0101$"), [2200]);
    assert_eq!(passing(br"ID match ^r\351$"), [1600]);
    refused(br"ID match ^r\1$", "Reference to nonexistent group");
    assert_eq!(passing(br"ID match ^{$"), [100, 1800]);
    assert_eq!(passing(br"ID match ^r}$"), [1900]);
    assert_eq!(passing(br"ID match ^r{1,$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r{,}$"), [1400]);
    assert_eq!(passing(br"ID match ^r{1,2}$"), [1200, 2000, 2100]);
    assert_eq!(passing(br"ID match ^r\{1,2}$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r{2}$"), [1200, 2100]);
    assert_eq!(passing(br"ID match ^r{x$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r{2,1}$"), [0u32; 0]);
    refused(br"ID match ^r\o{ 101 }$", "Missing right brace on \\o{}");
    assert_eq!(passing(br"ID match ^r\o{1_01}$"), [1000, 2400]);
    refused(br"ID match ^r\o{}$", "Empty \\o{}");
    assert_eq!(passing(br"ID match ^r\o{8}$"), [800]);
    assert_eq!(passing(br"ID match ^r\o{101x}$"), [1000, 2400]);
    assert_eq!(passing(br"ID match ^r\o{400}$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r\400$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r\o{0}$"), [800]);
    refused(br"ID match ^r\8$", "Reference to nonexistent group");
    assert_eq!(passing(br"ID match ^r\07$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r\010$"), [2300]);
    assert_eq!(passing(br"ID match ^r\10$"), [2300]);
    assert_eq!(passing(br"ID match ^r\1010$"), [2900]);
    assert_eq!(passing(br"ID match ^r\0A$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r\2000$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r{1}$"), [1200, 2000, 2100]);
    assert_eq!(passing(br"ID match ^{x}$"), [200, 2600]);
    assert_eq!(passing(br"ID match ^r[{]$"), [100]);
    assert_eq!(passing(br"ID match ^r{3,2,1}$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r+{$"), [100]);
    assert_eq!(passing(br"ID match ^.{$"), [100]);
    assert_eq!(passing(br"ID match ^(?:r){$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r{1}{$"), [100]);
    assert_eq!(passing(br"ID match ^r\o{1,01}$"), [1700]);
    assert_eq!(passing(br"ID match ^r\o{101}+$"), [1000, 2400]);
    assert_eq!(passing(br"ID match ^r\o{101}0?$"), [1000, 2400, 2900]);
    assert_eq!(passing(br"ID match ^r\1010?$"), [1000, 2400, 2900]);
    assert_eq!(passing(br"ID match ^r{,x}$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r{1,x}$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r[\o{101}]$"), [1000, 2400]);
    assert_eq!(passing(br"ID match ^r[\101]$"), [1000, 2400]);
    assert_eq!(passing(br"ID match ^r[\0]$"), [800]);
    assert_eq!(passing(br"ID match ^r[\1]$"), [1700]);
    assert_eq!(
        passing(br"ID match ^r\o{101}?$"),
        [1000, 1200, 2000, 2100, 2400]
    );
    assert_eq!(passing(br"ID match ^r\o{400}?\xe9$"), [1600]);
    assert_eq!(passing(br"ID match ^r\400?\xe9$"), [1600]);
    assert_eq!(passing(br"ID match ^r\N{U+41}$"), [1000, 2400]);
    assert_eq!(passing(br"ID match ^r\N{U+4_1}$"), [1000, 2400]);
    assert_eq!(passing(br"ID match ^r\N{U+41.30}$"), [2900]);
    assert_eq!(passing(br"ID match ^r[\N{U+41.30}]$"), [2900]);
    assert_eq!(passing(br"ID match ^r\N{U+41}{2}$"), [0u32; 0]);
    refused(
        br"ID match ^r\N{U+}$",
        "Invalid hexadecimal number in \\N{U+...}",
    );
    refused(
        br"ID match ^r\N{U+41_}$",
        "Invalid hexadecimal number in \\N{U+...}",
    );
    refused(
        br"ID match ^r\N{U+41.}$",
        "Invalid hexadecimal number in \\N{U+...}",
    );
    refused(
        br"ID match ^r[\N{U+41.30}-]$",
        "a \\N{U+...} escape of several characters before a `-` is not supported",
    );
}

/// The filter string's tokenizer groups a parenthesis as `filter_vep`'s does
/// (`FilterSet.pm:223-306`), inside a pattern word too: `ID match ^(r)\1$` is
/// the unit `ID match \1$` with the child `r ex`, and a unit with children
/// evaluates its children alone, so the pattern is never compiled; `r`
/// resolves to `REF` through the prefix search, which `vep filter` reads on
/// the first line as well, where `filter_vep` drops it (the documented
/// divergence). An unbalanced parenthesis, a word in a group that is no
/// operator and an empty group are refused as `filter_vep` refuses them.
/// Every other expected set of positions and error is the one `filter_vep`
/// writes on the same input.
#[test]
fn filter_string_groups_a_parenthesis_as_filter_vep_does() {
    let dir = tempfile::tempdir().unwrap();
    let run = byte_runs(dir.path(), ESC_IDS);
    let passing = |filter: &[u8]| kept(&run, filter);
    let refused = |filter: &[u8], message: &str| refused_with(&run, filter, message);
    let all: Vec<u32> = (1..=30).map(|i| i * 100).collect();
    let with_r = [
        100, 200, 300, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1900, 2000, 2100,
        2200, 2300, 2400, 2500, 2700, 2800, 2900, 3000,
    ];
    let with_a = [400, 500, 600, 700, 1000, 2400, 2900];
    assert_eq!(passing(br"ID match ^(r)\1$"), all);
    assert_eq!(passing(br"ID match ^(r)\10$"), all);
    assert_eq!(passing(br"r ex"), all);
    assert_eq!(passing(br"REF ex"), all);
    assert_eq!(passing(br"ID match (r)"), all);
    assert_eq!(passing(br"ID match (r)x"), all);
    assert_eq!(passing(br"ID match #(r)"), all);
    assert_eq!(passing(br"ID match ((r))"), all);
    assert_eq!(passing(br"ID match ^(r){2}$"), all);
    assert_eq!(passing(br"ID match ^(ID)$"), all);
    assert_eq!(passing(br"ID match (DP is 13)"), all);
    assert_eq!(passing(br"ID match (ID match ^(r))"), all);
    assert_eq!(passing(br"ID match ^(r)\1$ or ID match a"), all);
    assert_eq!(passing(br"ID match ^(r)\1$ and ID match a"), with_a);
    assert_eq!(passing(br"(ID match r)"), with_r);
    assert_eq!(passing(br"ID match (ID match r) and (REF is A)"), with_r);
    assert_eq!(passing(br"ID match a(b)c"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^(?:r)$"), [0u32; 0]);
    assert_eq!(passing(br"ID match ^r(?#c)$"), [0u32; 0]);
    assert_eq!(passing(br"ID match r{(2)}"), [0u32; 0]);
    assert_eq!(passing(br"not ID match (r)"), [0u32; 0]);
    assert_eq!(passing(br"ID match (not r)"), [0u32; 0]);
    let incomplete = "Error parsing filter string - incomplete parentheses sets?";
    for filter in [
        &br"ID match ^r("[..],
        br"ID match r)",
        br"ID match (r",
        br"ID match r)x",
        br"ID match x(",
        br"ID match (r)(",
        br"ID match ^(r)\1$)",
    ] {
        refused(filter, incomplete);
    }
    refused(
        br"ID match ^(r)(a)$",
        "Multiple data fields found matching query field \"a\": ALT, Allele",
    );
    refused(
        br"ID match ^(r|a)$",
        "Multiple data fields found matching query field \"r|a\": REF, ALT, QUAL, Allele, IMPACT, DISTANCE",
    );
    refused(
        br"ID match ^(D)\1$",
        "Multiple data fields found matching query field \"D\": DP, DISTANCE",
    );
    refused(
        br"ID match ()",
        "Multiple data fields found matching query field \"\": CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, DP, Allele, Consequence, IMPACT, SYMBOL, DISTANCE",
    );
    refused(br"ID (match) r", "No such operator \"r\"");
    refused(br"ID in (r)", "No list/file given");
}
