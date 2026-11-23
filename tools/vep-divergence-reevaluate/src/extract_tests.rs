// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::engine::input_chromosomes;
use crate::pyjson::loads;
use crate::record::load_reproductions;
use crate::repo::repo_root;
use crate::testing::{repro, Scratch, DELETE};
use std::io::Write;

const TABLE: &str = "## ENSEMBL VARIANT EFFECT PREDICTOR v116.2\n## Output produced by a stand-in\n## Extra column keys:\n## HGVSp : HGVS protein sequence name\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\nrec\t21:100\tA\tENSG1\tENST1\tTranscript\tmissense_variant,splice_region_variant\t10\t10\t4\tL/P\tctA/ccA\t-\tIMPACT=MODERATE;HGVSp=ENSP1:p.Leu4Pro;DISTANCE=0\nrec\t21:100\tA\tENSG1\tENST2\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\nother\t21:300\tC\tENSG3\tENST3\tTranscript\tdownstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=4000\n";

fn engine_run(output: String, exit_code: i64, stderr: &str, stdout: &str) -> EngineRun {
    EngineRun {
        engine: "vep-rs",
        command: vec!["vep".to_string()],
        exit_code,
        output,
        stdout: stdout.to_string(),
        stderr: stderr.to_string(),
    }
}

/// What the compare edits read from a run.
fn read(edits: &[(&str, &str)], run: &EngineRun) -> Vec<String> {
    extract(&repro(edits), run).expect("a value")
}

#[test]
fn rows_field_and_count_read_the_table_by_where() {
    let scratch = Scratch::new("table");
    let run = engine_run(scratch.write("default.txt", TABLE.as_bytes()), 0, "", "");
    let rows = [("compare.kind", "'rows'"), ("compare.field", DELETE)];
    assert_eq!(
        read(&rows, &run),
        [
            "21:100 | A | ENST1 | missense_variant,splice_region_variant",
            "21:100 | A | ENST2 | intron_variant"
        ]
    );
    assert_eq!(
        read(
            &[
                rows[0],
                rows[1],
                ("compare.where", "{record = 'rec', feature = 'ENST1'}"),
                (
                    "compare.columns",
                    "['Feature', 'Codons', 'DISTANCE', 'SIFT']"
                ),
            ],
            &run
        ),
        ["ENST1 | ctA/ccA | 0 | -"]
    );
    assert_eq!(read(&[], &run), ["ENSP1:p.Leu4Pro", "-"]);
    assert_eq!(
        read(&[("compare.with_key", "true")], &run),
        ["ENST1=ENSP1:p.Leu4Pro", "ENST2=-"]
    );
    assert_eq!(
        read(
            &[
                ("compare.field", "'Protein_position'"),
                (
                    "compare.where",
                    "{consequence = 'splice_region_variant,missense_variant'}"
                ),
            ],
            &run
        ),
        ["4"]
    );
    let count = [("compare.kind", "'count'"), ("compare.field", DELETE)];
    assert_eq!(
        read(
            &[
                count[0],
                count[1],
                ("compare.where", "{feature = ['ENST1', 'ENST3']}")
            ],
            &run
        ),
        ["2"]
    );
    assert_eq!(
        read(
            &[
                count[0],
                count[1],
                ("compare.where", "{pattern = 'DISTANCE=4000'}")
            ],
            &run
        ),
        ["1"]
    );
    assert_eq!(
        read(
            &[rows[0], rows[1], ("compare.where", "{record = 'absent'}")],
            &run
        ),
        Vec::<String>::new()
    );
}

#[test]
fn header_lines_warnings_streams_exit_and_stats_read_their_files() {
    let scratch = Scratch::new("files");
    let output = scratch.write("default.txt", TABLE.as_bytes());
    scratch.write(
        "default.txt_warnings.txt",
        b"WARNING: line 3 : skipped /host/in.vcf\nWARNING: line 4 : kept\n",
    );
    scratch.write(
        "default.txt_summary.txt",
        b"[General statistics]\nLines of input read\t2\n\n[Position in protein]\n00-10%\t1\n10-20%\t2\n\n[Other]\nx\t1\n",
    );
    let run = engine_run(
        output,
        2,
        "WARNING: x\nERROR: [input] bad at /host/path/in.vcf\n",
        "done\n",
    );
    let kind = |kind: &'static str, key: &'static str, value: &'static str| {
        [
            ("compare.kind", kind),
            ("compare.field", DELETE),
            (key, value),
        ]
    };
    assert_eq!(
        read(&kind("'header'", "compare.pattern", "'^## HGVSp'"), &run),
        ["## HGVSp : HGVS protein sequence name"]
    );
    assert_eq!(
        read(&kind("'lines'", "compare.pattern", "'ENST3'"), &run),
        [TABLE.lines().nth(7).unwrap()]
    );
    assert_eq!(
        read(&kind("'warnings'", "compare.pattern", "'skipped'"), &run),
        ["WARNING: line 3 : skipped /host/in.vcf"]
    );
    assert_eq!(
        read(&kind("'stderr'", "compare.pattern", "'^ERROR'"), &run),
        ["ERROR: [input] bad at /host/path/in.vcf"]
    );
    assert_eq!(
        read(&kind("'stdout'", "compare.pattern", "'done'"), &run),
        ["done"]
    );
    assert_eq!(
        read(
            &[("compare.kind", "'exit'"), ("compare.field", DELETE)],
            &run
        ),
        ["2"]
    );
    assert_eq!(
        read(
            &kind("'stats'", "compare.section", "'[Position in protein]'"),
            &run
        ),
        ["00-10%\t1", "10-20%\t2"]
    );
    assert_eq!(
        read(&kind("'stats'", "compare.section", "'[Missing]'"), &run),
        Vec::<String>::new()
    );
    // A byte that is not UTF-8 reads as its escape.
    let raw = scratch.write(
        "raw.txt",
        b"#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\nr\xe9c\t21:1\tA\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t-\n",
    );
    assert_eq!(
        read(
            &[
                ("compare.field", "'Uploaded_variation'"),
                ("compare.where", DELETE)
            ],
            &engine_run(raw, 0, "", "")
        ),
        ["r\\xe9c"]
    );
    // A gzip file is read whatever its name.
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder
        .write_all(b"##fileformat=VCFv4.2\n#CHROM\tPOS\n22\t5\n21\t7\n")
        .unwrap();
    let zipped = scratch.write("zipped.bgz", &encoder.finish().unwrap());
    assert_eq!(input_chromosomes(&repro(&[]), &zipped), ["22", "21"]);
    let no_files = engine_run(scratch.path("missing.txt"), 1, "", "");
    assert_eq!(
        read(&kind("'warnings'", "compare.pattern", "'.'"), &no_files),
        Vec::<String>::new()
    );
    assert_eq!(
        read(
            &[
                ("compare.kind", "'rows'"),
                ("compare.field", DELETE),
                ("compare.where", DELETE)
            ],
            &no_files
        ),
        Vec::<String>::new()
    );
}

#[test]
fn json_reads_a_path_over_the_matching_records() {
    let scratch = Scratch::new("json");
    let records = concat!(
        r#"{"id": "rec", "input": "21\t100\trec\tG\tA", "seq_region_name": "21", "start": 100, "most_severe_consequence": "missense_variant", "transcript_consequences": [{"transcript_id": "ENST1", "hgvsp": "ENSP1:p.Leu4Pro", "consequence_terms": ["missense_variant"]}, {"transcript_id": "ENST2", "consequence_terms": ["intron_variant"]}]}"#,
        "\n",
        r#"{"id": "other", "input": "21\t300\tother\tG\tC", "seq_region_name": "21", "start": 300, "most_severe_consequence": "downstream_gene_variant"}"#,
        "\n"
    );
    let run = engine_run(scratch.write("json.jsonl", records.as_bytes()), 0, "", "");
    let json = |path: &'static str, where_: &'static str| {
        [
            ("compare.kind", "'json'"),
            ("compare.field", DELETE),
            ("compare.path", path),
            ("compare.where", where_),
        ]
    };
    assert_eq!(
        read(
            &json(
                "'transcript_consequences[transcript_id=ENST1].hgvsp'",
                "{id = 'rec'}"
            ),
            &run
        ),
        ["ENSP1:p.Leu4Pro"]
    );
    assert_eq!(
        read(
            &json(
                "'transcript_consequences[].consequence_terms'",
                "{record = 'rec'}"
            ),
            &run
        ),
        [r#"["missense_variant"]"#, r#"["intron_variant"]"#]
    );
    assert_eq!(
        read(&json("'most_severe_consequence'", DELETE), &run),
        ["missense_variant", "downstream_gene_variant"]
    );
    assert_eq!(
        read(
            &json("'most_severe_consequence'", "{location = '21:300'}"),
            &run
        ),
        ["downstream_gene_variant"]
    );
    assert_eq!(
        read(
            &json("'start'", r#"{input = "21\t300\tother\tG\tC"}"#),
            &run
        ),
        ["300"]
    );
    assert_eq!(
        read(
            &json("'transcript_consequences[].hgvsp'", "{id = 'other'}"),
            &run
        ),
        Vec::<String>::new()
    );
    let document = loads(r#"{"a": [{"k": "1", "v": 1}, {"k": "2", "v": 2}]}"#).unwrap();
    assert_eq!(
        json_path(&document, "a[k=2].v"),
        [&Value::Int("2".to_string())]
    );
    let document = loads(r#"{"a": {"b": null}}"#).unwrap();
    assert_eq!(json_path(&document, "a.b"), [&Value::None]);
    assert_eq!(json_value(&Value::None), "null");
    // A record that is no object stops the entry, in the words of the missing method.
    let list = engine_run(scratch.write("list.jsonl", b"[1, 2]\n"), 0, "", "");
    assert_eq!(
        extract(&repro(&json("'a'", DELETE)), &list),
        Err(RunError("'list' object has no attribute 'get'".to_string()))
    );
}

#[test]
fn vcf_reads_csq_fields_by_record_feature_and_allele() {
    let scratch = Scratch::new("vcf");
    let vcf = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|HGVSp\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t100\trec\tG\tA,C\t.\t.\tDP=4;CSQ=A|missense_variant|ENST1|ENSP1:p.Leu4Pro,A|intron_variant|ENST2|,C|missense_variant|ENST1|ENSP1:p.Leu4Arg\n21\t300\tother\tG\tC\t.\t.\t.\n";
    let run = engine_run(scratch.write("vcf.vcf", vcf.as_bytes()), 0, "", "");
    let field = |field: &'static str, where_: &'static str| {
        [
            ("compare.kind", "'vcf'"),
            ("compare.field", field),
            ("compare.where", where_),
        ]
    };
    assert_eq!(
        read(&field("'HGVSp'", "{id = 'rec', feature = 'ENST1'}"), &run),
        ["ENSP1:p.Leu4Pro", "ENSP1:p.Leu4Arg"]
    );
    assert_eq!(
        read(
            &field("'HGVSp'", "{id = 'rec', feature = 'ENST1', allele = 'C'}"),
            &run
        ),
        ["ENSP1:p.Leu4Arg"]
    );
    assert_eq!(
        read(
            &field("'Consequence'", "{location = '21:100', allele = 'A'}"),
            &run
        ),
        ["missense_variant", "intron_variant"]
    );
    assert_eq!(read(&field("'CSQ'", "{id = 'other'}"), &run), ["-"]);
    assert_eq!(
        read(
            &[
                ("compare.kind", "'vcf'"),
                ("compare.field", DELETE),
                ("compare.line", "true"),
                ("compare.where", "{pos = '300'}")
            ],
            &run
        ),
        ["21\t300\tother\tG\tC\t.\t.\t."]
    );
    assert_eq!(
        read(
            &[
                ("compare.kind", "'vcf'"),
                ("compare.field", "'DP'"),
                ("compare.info", "'DP'"),
                ("compare.where", "{id = 'rec'}")
            ],
            &run
        ),
        ["4"]
    );
}

fn compare(text: &str) -> Value {
    crate::record::parse_toml(text.as_bytes()).expect("a compare table")
}

fn list(values: &[&str]) -> Vec<String> {
    values.iter().map(|v| v.to_string()).collect()
}

#[test]
fn normalise_masks_paths_and_patterns_and_sorts_unless_told_otherwise() {
    let values = list(&["b at /opt/vep/src/x.pm line 3", "a pid 1234"]);
    let plain = compare("");
    let normal = |values: &[String], compare: &Value| normalise(values, compare).unwrap();
    assert_eq!(
        normal(&values, &plain),
        ["a pid 1234", "b at <path> line 3"]
    );
    assert_eq!(
        normal(&values, &compare("sorted = false\nmask_paths = false")),
        values
    );
    assert_eq!(
        normal(&values, &compare(r"mask = ['\d+']")),
        ["a pid <masked>", "b at <path> line <masked>"]
    );
    assert_eq!(
        normal(&list(&["ENSP1:p.Leu4Pro"]), &plain),
        ["ENSP1:p.Leu4Pro"]
    );
    assert_eq!(
        normal(&list(&["a description "]), &plain),
        ["a description"]
    );
    assert_eq!(
        normal(
            &list(&["-/A", "C//G", "*/Q", "21:100 | -/A | ENST1 | x"]),
            &plain
        ),
        ["*/Q", "-/A", "21:100 | -/A | ENST1 | x", "C//G"]
    );
    assert_eq!(
        normal(&list(&["gzip: /work/in/x.bgz: not in gzip format"]), &plain),
        ["gzip: <path> not in gzip format"]
    );
    assert_eq!(
        normal(&list(&["a description "]), &compare("strip = false")),
        ["a description "]
    );
    assert_eq!(
        normal(&list(&["21:100 | A | ENST1 | x"]), &plain),
        ["21:100 | A | ENST1 | x"]
    );
}

#[test]
fn verdict_is_closed_present_or_changed() {
    let expected = compare("ensembl_vep = ['0.00']\nvep_rs = ['-']");
    let plain = compare("");
    let judge = |expected: &Value, ensembl: &[&str], vep_rs: &[&str]| {
        verdict(expected, &plain, &list(ensembl), &list(vep_rs)).unwrap()
    };
    assert_eq!(judge(&expected, &["0.00"], &["-"]), "still-present");
    assert_eq!(judge(&expected, &["-"], &["-"]), "closed");
    assert_eq!(judge(&expected, &["0.00"], &["0.34"]), "changed");
    assert_eq!(judge(&expected, &["0.01"], &["-"]), "changed");
    let single = compare("ensembl_vep = 'x'\nvep_rs = 'y'");
    assert_eq!(judge(&single, &["x"], &["y"]), "still-present");
}

#[test]
fn a_failed_engine_run_is_an_error_unless_its_exit_is_recorded() {
    let scratch = Scratch::new("failure");
    let present = scratch.write("default.txt", b"#Uploaded_variation\n");
    let missing = scratch.path("missing.txt");
    let run = |engine: &'static str, exit_code: i64, output: &str, stderr: &str| EngineRun {
        engine,
        command: vec!["vep".to_string()],
        exit_code,
        output: output.to_string(),
        stdout: String::new(),
        stderr: stderr.to_string(),
    };
    let ok = run("ensembl-vep", 0, &present, "");
    let failed = run("ensembl-vep", 2, &missing, "ERROR: no cache\n");
    let silent = run("vep-rs", 0, &missing, "");
    let rows = [("compare.kind", "'rows'"), ("compare.field", DELETE)];
    let failure = |edits: &[(&str, &str)], engine: &str, run: &EngineRun| {
        run_failure(&repro(edits), engine, run)
    };
    assert_eq!(failure(&rows, "ensembl-vep", &ok), None);
    assert_eq!(
        failure(&rows, "ensembl-vep", &failed).as_deref(),
        Some("ensembl-vep exited 2: ERROR: no cache")
    );
    assert_eq!(
        failure(&rows, "vep-rs", &silent).as_deref(),
        Some("vep-rs wrote no output file")
    );
    let recorded = [
        rows[0],
        rows[1],
        ("expected.ensembl_vep", "[]"),
        ("expected.vep_rs", "['x']"),
        ("expected.ensembl_vep_exit", "2"),
    ];
    assert_eq!(failure(&recorded, "ensembl-vep", &failed), None);
    assert_eq!(
        failure(&recorded, "ensembl-vep", &ok).as_deref(),
        Some("ensembl-vep exited 0, the entry records 2: ")
    );
    let zero = [
        rows[0],
        rows[1],
        ("expected.ensembl_vep", "[]"),
        ("expected.vep_rs", "['x']"),
        ("expected.ensembl_vep_exit", "0"),
    ];
    assert_eq!(
        failure(&zero, "ensembl-vep", &run("ensembl-vep", 0, &missing, "")).as_deref(),
        Some("ensembl-vep wrote no output file")
    );
    for kind in ["'exit'", "'stderr'", "'stdout'"] {
        let edits = [
            ("compare.kind", kind),
            ("compare.field", DELETE),
            ("compare.pattern", "'x'"),
        ];
        assert_eq!(failure(&edits, "ensembl-vep", &failed), None, "{kind}");
    }
    let warnings = [
        ("compare.kind", "'warnings'"),
        ("compare.field", DELETE),
        ("compare.pattern", "'x'"),
    ];
    assert_eq!(failure(&warnings, "vep-rs", &silent), None);
    assert_eq!(
        failure(&warnings, "vep-rs", &failed).as_deref(),
        Some("vep-rs exited 2: ERROR: no cache")
    );
}

/// Ensembl VEP opens a compressed input once for each of its two parser objects
/// (`Parser.pm:133`, reached again from `Parser.pm:159`, through `Utils.pm:441`) and reads
/// only the second `gzip` to its end: the first reports an error that follows its output
/// only where that output is written before its pipe closes, so the count of a `gzip` line
/// on standard error differs from run to run. No committed run that reads standard error
/// records such a line for Ensembl VEP.
#[test]
fn no_committed_run_holds_ensembl_vep_to_a_count_of_gzip_lines() {
    let directory = format!("{}/scripts/divergences/reproductions", repo_root());
    let repros = load_reproductions(None, &directory).expect("the committed records");
    let mut read = 0;
    let mut held = Vec::new();
    for repro in &repros {
        for (place, run) in repro.runs().iter().enumerate() {
            if run.compare().get_str("kind") != Some("stderr") {
                continue;
            }
            read += 1;
            for line in as_list(run.expected().section("ensembl_vep")) {
                if line.starts_with("gzip:") {
                    held.push(format!("{} run {place}: {line}", repro.name));
                }
            }
        }
    }
    assert!(read > 0, "no committed run reads standard error");
    assert_eq!(held, Vec::<String>::new());
}
