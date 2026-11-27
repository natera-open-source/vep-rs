// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The structural-variant classes: the rows of a record the reference skipped, a breakend's
//! mate-side row read at the local coordinate, and the own-end rows the reference lost,
//! through the library and through the binary on planted files.

mod compare_common;

use std::collections::BTreeSet;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;

use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::cache::THREADS_VARIABLE;
use vep_formats_compare::classes::tables;
use vep_formats_compare::classes::INDIVIDUAL_ALLELE_ORDER;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::compare_files::compare_files;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Error;
use vep_formats_compare::input::read_text;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::json_objects;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::json::JSON_IDENTITY_KEY;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::remainder::hand_over;
use vep_formats_compare::remainder::FrontEnd;
use vep_formats_compare::remainder::STREAM;
use vep_formats_compare::remainder::WHOLE_VARIABLE;
use vep_formats_compare::rules::structural::StableId;
use vep_formats_compare::rules::structural::StructuralVariants;
use vep_formats_compare::rules::Side;

/// How one run of the binary ended.
struct Ran {
    code: i32,
    out: String,
    err: String,
}

/// The binary on `args`.
fn run(args: &[&str], whole: bool) -> Ran {
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-formats-compare"));
    command.args(args).env(THREADS_VARIABLE, "2");
    if whole {
        command.env(WHOLE_VARIABLE, "1");
    } else {
        command.env_remove(WHOLE_VARIABLE);
    }
    let done = command.output().expect("the tool runs");
    Ran {
        code: done.status.code().expect("the tool exits"),
        out: String::from_utf8(done.stdout).expect("the standard output is text"),
        err: String::from_utf8(done.stderr).expect("the standard error is text"),
    }
}

fn path_text(path: &Path) -> &str {
    path.to_str().expect("a path of the tests is text")
}

fn report_at(path: &Path) -> Value {
    json(&std::fs::read_to_string(path).expect("the report reads"))
}

/// An input VCF of `records` (tab-separated data lines).
fn structural_input(path: &Path, records: &[&str]) -> PathBuf {
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for record in records {
        text.push_str(record);
        text.push('\n');
    }
    std::fs::write(path, text).expect("the input is written");
    path.to_path_buf()
}

/// A candidate JSON cache of two transcripts near the breakends of the tests: ENST21 on
/// chromosome 21 at 100500-110000 and ENST22 on chromosome 22 at 5000100-5010000, both
/// forward-strand non-coding, ENST22 two exons; ENST21_FAR lies out of a breakend's reach.
fn breakend_cache(root: &Path) -> PathBuf {
    let shard = |chromosome: &str, name: &str, text: &str| {
        let directory = root.join("transcripts").join(chromosome);
        std::fs::create_dir_all(&directory).expect("the directory is made");
        std::fs::write(directory.join(name), text).expect("the shard is written");
    };
    shard(
        "21",
        "1-1000000.json",
        r#"[{"stable_id": "ENST21", "start": 100500, "end": 110000, "strand": 1, "biotype": "lncRNA", "exons": [{"start": 100500, "end": 110000}]},
            {"stable_id": "ENST21_FAR", "start": 200000, "end": 210000, "strand": 1, "biotype": "lncRNA", "exons": [{"start": 200000, "end": 210000}]}]"#,
    );
    shard(
        "22",
        "5000001-6000000.json",
        r#"[{"stable_id": "ENST22", "start": 5000100, "end": 5010000, "strand": 1, "biotype": "lncRNA",
             "exons": [{"start": 5000100, "end": 5000500}, {"start": 5005000, "end": 5010000}]}]"#,
    );
    root.to_path_buf()
}

/// The side tables of a comparison that names the structural-variant one alone.
fn side_of(structural: StructuralVariants) -> Side {
    Side {
        structural_variants: Some(structural),
        ..Side::default()
    }
}

fn structural(classes: &[&str], input: Option<&Path>, cache: Option<&Path>) -> StructuralVariants {
    StructuralVariants::new(classes, input, cache).expect("the side table reads")
}

fn refusal(classes: &[&str], input: Option<&Path>, cache: Option<&Path>) -> String {
    match StructuralVariants::new(classes, input, cache) {
        Err(Error::Comparator(text)) => text,
        Ok(_) => panic!("the side table read"),
    }
}

/// The report of two entry lists beside `side`.
fn compare_beside(
    side: &Side,
    expected: &[Entry],
    actual: &[Entry],
    documented: &mut Documented,
) -> Value {
    compare_entries(expected, actual, documented, &CompareOptions::new(side))
        .expect("the entries compare")
}

/// The report of two default-format texts under the built-in `classes` beside `side`, and
/// what was applied.
fn compare_structural(
    reference: &str,
    candidate: &str,
    side: &Side,
    classes: &[&str],
) -> (Value, Documented) {
    let mut documented = documented_with(classes);
    let expected = parse_default(reference).expect("the reference parses");
    let actual = parse_default(candidate).expect("the candidate parses");
    let report = compare_beside(side, &expected, &actual, &mut documented);
    (report, documented)
}

fn applied(documented: &Documented) -> std::collections::BTreeMap<String, i64> {
    counts(&documented.applied)
}

fn examples<'a>(report: &'a Value, class: &str) -> Vec<&'a str> {
    texts(at(report, &["mask_examples", class]))
}

const SKIPPED: &str = "reference_skipped_record";
const MATE_READ: &str = "breakend_mate_local_read";
const OWN_ROWS: &str = "bnd_own_chromosome_rows_lost";

/// The data lines of an input VCF whole. The comparison keeps of each no more than it reads
/// (`input_records`); a test that writes an output from a record, or reads a later column of
/// it, reads the lines itself.
fn whole_records(path: &Path) -> Vec<String> {
    let text = read_text(path).expect("the input reads");
    data_lines(&text)
}

/// Parser/VCF.pm:464 stores a QUAL of `0` as `.`, so the `input` of the committed
/// GRCh38-vcf-edge-cases JSON object for `sv_qual0` (record 40, QUAL `0`) is not the input
/// line. The object still keys to its record, as a VCF line does, so a candidate output
/// carrying the input-only object of the reference-skipped record 26 (the manifest's
/// documented `reference_skipped_record` key) before it compares green under the class and
/// shifts no later object onto its ordinal.
#[test]
fn a_json_object_names_its_record_by_the_first_five_columns() {
    let directory = scratch("structural-json-objects");
    let corpus = golden().join("GRCh38-vcf-edge-cases");
    let inputs = whole_records(&corpus.join("variants.vcf"));
    let manifest =
        json(&std::fs::read_to_string(corpus.join("manifest.json")).expect("the manifest reads"));
    let text = read_text(&corpus.join("expected").join("json.jsonl.gz")).expect("the output reads");
    let objects = json_objects(&text).expect("the objects parse");
    let columns: Vec<&str> = inputs[40].split('\t').collect();
    assert_eq!(columns[2..6], ["sv_qual0", "N", "<DEL>", "0"]);
    let echo = objects
        .iter()
        .find(|(_, object)| object.get("id").and_then(Value::as_str) == Some("sv_qual0"))
        .and_then(|(_, object)| object.get("input"))
        .and_then(Value::as_str)
        .expect("the object of the record");
    assert_eq!(echo, inputs[40].replace("\t0\tPASS\t", "\t.\tPASS\t"));
    assert_ne!(echo, inputs[40]);
    let entries = parse_json(&text, &inputs).expect("the output parses");
    let keys: Vec<&str> = entries
        .iter()
        .filter(|e| e.record == Some(40))
        .map(|e| e.key.as_str())
        .collect();
    assert_eq!(
        keys,
        ["rec40|record", "rec40|deletion|ENST00000286808|transcript"]
    );
    assert!(entries.iter().all(|e| e.record.is_some()));

    let extras = at(&manifest, &["vep_rs_only_tuples"])
        .as_list()
        .expect("a list");
    assert_eq!(extras.len(), 1);
    assert_eq!(
        extras[0].get("expected_divergence").and_then(Value::as_str),
        Some(SKIPPED)
    );
    assert_eq!(
        extras[0].get("record_indices"),
        Some(&Value::List(vec![Value::Int(26)]))
    );
    let records: Vec<&Entry> = entries
        .iter()
        .filter(|e| e.key.ends_with("|record"))
        .collect();
    let position = objects
        .iter()
        .zip(&records)
        .find(|(_, e)| e.record > Some(26))
        .map(|((n, _), _)| *n as usize)
        .expect("an object after the skipped record");
    let identified = |object: &Value| {
        let mut members = object.as_dict().expect("an object").clone();
        members.insert(
            JSON_IDENTITY_KEY.to_string(),
            Value::Str("0.5.0".to_string()),
        );
        pyjson::dumps(&Value::Dict(members), &pyjson::DumpOptions::default())
    };
    let mut lines: Vec<String> = objects
        .iter()
        .map(|(_, object)| identified(object))
        .collect();
    let input_only = Value::Dict(IndexMap::from([(
        "input".to_string(),
        Value::Str(inputs[26].clone()),
    )]));
    lines.insert(position, identified(&input_only));
    let expected_file = planted(&directory, "expected.jsonl", &text);
    let actual_file = planted(&directory, "actual.jsonl", &(lines.join("\n") + "\n"));

    let compared = |classes: &[&str]| {
        let masked = masks(classes);
        let mut documented = from_run(&manifest, "default", classes);
        compare_files(
            "json",
            &expected_file,
            &actual_file,
            Some(&inputs),
            &mut documented,
            &masked,
            &IndexMap::new(),
            &IndexMap::new(),
            4,
            false,
            None,
            None,
            None,
            &mut Side::default(),
        )
        .expect("the files compare")
    };
    let report = compared(&[SKIPPED]);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(number(&report, &["counts", "expected_entries"]), 107);
    assert_eq!(number(&report, &["counts", "actual_entries"]), 108);
    assert_eq!(number(&report, &["counts", "compared"]), 107);
    assert_eq!(
        number(
            &report,
            &["counts", "record_entries_on_reference_skipped_records"]
        ),
        1
    );
    assert_eq!(number(&report, &["one_sided", "expected_only"]), 0);
    assert_eq!(number(&report, &["one_sided", "actual_only"]), 0);
    let unmasked = compared(&[]);
    assert_eq!(verdict(&unmasked), "red");
    assert_eq!(number(&unmasked, &["one_sided", "expected_only"]), 0);
    assert_eq!(number(&unmasked, &["one_sided", "actual_only"]), 1);
    let unexpected: Vec<&str> = at(&unmasked, &["unexpected"])
        .as_dict()
        .expect("an object")
        .values()
        .flat_map(|group| texts(at(group, &["examples"])))
        .collect();
    assert_eq!(unexpected, ["rec26|record"]);
}

/// The manifest documents five vep-rs-only keys on records the reference skipped; the
/// candidate's consequence rows on them are set aside only under the class's mask, and its
/// record-level VCF line on such a record is counted, not compared.
#[test]
fn the_rows_of_a_documented_skipped_record_are_set_aside_under_the_class() {
    // The input's lines whole: the candidate's VCF takes one of them for an output line.
    let inputs = whole_records(&corpus().join("variants.vcf"));
    let (default, vcf, manifest) = (text_of("default"), text_of("vcf"), manifest());
    let skipped_row = "multi_ins_del\t21:30165801-30166000\tinsertion\tENSG00000156282\tENST00000286808\tTranscript\t\
                       feature_elongation,coding_sequence_variant,3_prime_UTR_variant\t-\t-\t-\t-\t-\t-\tIMPACT=HIGH;STRAND=1";
    let mut rows = data_lines(&default);
    rows.push(skipped_row.to_string());
    let expected = parse_default(&default).expect("the reference parses");
    let actual = parse_default(&rebuild(&default, &rows)).expect("the candidate parses");
    let unmasked = compare(&expected, &actual, &mut from_manifest(&manifest, &[]));
    assert_eq!(number(&unmasked, &["one_sided", "actual_only"]), 1);
    let mut masked = from_run(&manifest, "default", &[SKIPPED]);
    assert_eq!(verdict(&compare(&expected, &actual, &mut masked)), "green");

    let mut vcf_rows = data_lines(&vcf);
    vcf_rows.push(inputs[12].clone());
    let vcf_actual = parse_vcf(&rebuild(&vcf, &vcf_rows), &inputs).expect("the candidate parses");
    let vcf_expected = parse_vcf(&vcf, &inputs).expect("the reference parses");
    let report = compare(&vcf_expected, &vcf_actual, &mut masked);
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        number(
            &report,
            &["counts", "record_entries_on_reference_skipped_records"]
        ),
        1
    );
    let bare = compare(&vcf_expected, &vcf_actual, &mut no_documented());
    assert_eq!(number(&bare, &["one_sided", "actual_only"]), 1);
}

/// A run reads the input records for the vcf and json formats and for a structural-variant
/// class alone: a default or tab comparison beside `--manifest` names `--input` for the flag
/// set documented over that input, opens none of its lines (the file need not exist), and
/// leaves out the runs documented over another input.
#[test]
fn a_text_format_reads_the_input_name_for_the_manifest_alone() {
    let directory = scratch("structural-input-name");
    let reference = default_text(&[default_row(
        "rs1",
        "21:100",
        "A",
        "ENST1",
        "intron_variant",
        "IMPACT=MODIFIER",
    )]);
    let expected = planted(&directory, "ref.txt", &reference);
    let actual = planted(
        &directory,
        "cand.txt",
        &reference.replace("intron_variant", "missense_variant"),
    );
    let item = r#"{"location": "21:100", "allele": "A", "feature": "ENST1", "vep_rs_consequence_sets": ["missense_variant"], "expected_divergence": "some_class"}"#;
    let manifest = planted(
        &directory,
        "manifest.json",
        &format!(
            r#"{{"divergences": [], "flag_sets": [{{"name": "documented", "input": "documented.vcf", "divergences": [{item}]}}, {{"name": "other", "input": "other.vcf", "divergences": []}}]}}"#
        ),
    );
    let report = directory.join("report.json");
    let documented = directory.join("documented.vcf");
    let other = directory.join("other.vcf");
    let argv = [
        "--expected",
        path_text(&expected),
        "--actual",
        path_text(&actual),
        "--manifest",
        path_text(&manifest),
        "--report",
        path_text(&report),
    ];
    let with = |before: &[&str], input: &Path| -> Ran {
        let mut args: Vec<&str> = before.to_vec();
        args.extend(argv);
        args.extend(["--input", path_text(input)]);
        run(&args, false)
    };
    for fmt in ["default", "tab"] {
        let named = with(&["--format", fmt, "--masks", "some_class"], &documented);
        assert_eq!(named.code, 0, "{}", named.err);
        assert_eq!(
            at(&report_at(&report), &["masks", "applied"]),
            &json(r#"{"some_class": 1}"#)
        );
        let another = with(&["--format", fmt, "--masks", "some_class"], &other);
        assert_eq!(another.code, 1, "{}", another.err);
        assert_eq!(at(&report_at(&report), &["masks", "applied"]), &json("{}"));
    }
    for reading in [
        &["--format", "vcf"][..],
        &["--format", "default", "--masks", SKIPPED][..],
    ] {
        let unread = with(reading, &documented);
        assert_eq!(unread.code, 2, "{}", unread.err);
        let opening = format!(
            "ERROR: [compare_vep_formats] cannot read {}: ",
            path_text(&documented)
        );
        assert!(unread.err.starts_with(&opening), "{}", unread.err);
    }
}

/// Ensembl VEP drops a record of an unsupported type (Parser/VCF.pm:477-481) or a span over
/// --max_sv_size (Parser.pm:491-501) before annotation and writes nothing for it; the
/// candidate annotates it. Every candidate-only row of such a record is the class, decided by
/// the structural-variant comparator's skip routes over the input VCF, counted per row in
/// every format, the record's VCF line and JSON object counted and not compared; a
/// candidate-only row on a record the reference kept, a skip-route record the reference wrote
/// rows for (a run under --dont_skip) and a key a manifest documents are outside the
/// built-in's count.
#[test]
fn the_rows_of_a_record_the_reference_drops_are_set_aside() {
    let directory = scratch("structural-skipped-record");
    let records = [
        "21\t100\tdel1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=300",
        "21\t1000\tcpx1\tG\t<CPX>\t.\t.\tSVTYPE=CPX;END=1200",
    ];
    let vcf_input = structural_input(&directory.join("in.vcf"), &records);
    let kept = default_row(
        "del1",
        "21:101-300",
        "deletion",
        "ENST1",
        "intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    );
    let dropped: Vec<String> = ["ENST1", "ENST1b"]
        .iter()
        .map(|feature| {
            default_row(
                "cpx1",
                "21:1001-1200",
                "CPX",
                feature,
                "intron_variant",
                "IMPACT=MODIFIER;STRAND=1",
            )
        })
        .collect();
    let reference = default_text(&[&kept]);
    let candidate = default_text(&[&kept, &dropped[0], &dropped[1]]);
    let side = side_of(structural(&[SKIPPED], Some(&vcf_input), None));
    let bare = Side::default();
    let (report, documented) = compare_structural(&reference, &candidate, &side, &[SKIPPED]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(SKIPPED, 2)]));
    assert_eq!(
        examples(&report, SKIPPED),
        [
            "21:1001-1200|CPX|ENST1|Transcript: candidate-only row on a record the reference skipped, 'intron_variant'",
            "21:1001-1200|CPX|ENST1b|Transcript: candidate-only row on a record the reference skipped, 'intron_variant'",
        ]
    );
    let without = compare_structural(&reference, &candidate, &bare, &[SKIPPED]).0;
    assert_eq!(number(&without, &["one_sided", "actual_only"]), 2);
    assert_eq!(
        verdict(&compare_default(&reference, &candidate, &[]).0),
        "red"
    );
    // A candidate-only row on the kept record, and a skip-route record the reference wrote for.
    let extra_kept = default_row(
        "del1",
        "21:101-300",
        "deletion",
        "ENST1b",
        "intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    );
    let (red, documented) = compare_structural(
        &reference,
        &default_text(&[&kept, &extra_kept, &dropped[0], &dropped[1]]),
        &side,
        &[SKIPPED],
    );
    assert_eq!(verdict(&red), "red");
    assert_eq!(number(&red, &["one_sided", "actual_only"]), 1);
    assert_eq!(applied(&documented), count_map(&[(SKIPPED, 2)]));
    let (red, documented) = compare_structural(
        &default_text(&[&kept, &dropped[0]]),
        &candidate,
        &side,
        &[SKIPPED],
    );
    assert_eq!(verdict(&red), "red");
    assert_eq!(number(&red, &["one_sided", "actual_only"]), 1);
    assert_eq!(applied(&documented), count_map(&[]));
    // A manifest's documented tuple is the manifest's; the built-in counts the other row.
    let manifest = json(
        r#"{"vep_rs_only_tuples": [{"location": "21:1001-1200", "allele": "CPX", "feature": "ENST1", "feature_type": "Transcript",
            "record_indices": [1], "expected_divergence": "reference_skipped_record"}]}"#,
    );
    let mut documented = from_manifest(&manifest, &[SKIPPED]);
    let report = compare_beside(
        &side,
        &parse_default(&reference).expect("the reference parses"),
        &parse_default(&candidate).expect("the candidate parses"),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(SKIPPED, 2)]));
    // The VCF and JSON formats key the rows by record and count the record entry apart.
    let inputs = whole_records(&vcf_input);
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. \
                      Format: Allele|Consequence|Feature_type|Feature\">\n";
    let kept_line = format!(
        "{};CSQ=deletion|intron_variant|Transcript|ENST1\n",
        records[0]
    );
    let dropped_line = format!(
        "{};CSQ=CPX|intron_variant|Transcript|ENST1,CPX|intron_variant|Transcript|ENST1b\n",
        records[1]
    );
    let vcf_expected =
        parse_vcf(&format!("{csq_header}{kept_line}"), &inputs).expect("the reference parses");
    let vcf_actual = parse_vcf(&format!("{csq_header}{kept_line}{dropped_line}"), &inputs)
        .expect("the candidate parses");
    let mut documented = documented_with(&[SKIPPED]);
    let report = compare_beside(&side, &vcf_expected, &vcf_actual, &mut documented);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(SKIPPED, 2)]));
    assert_eq!(
        number(
            &report,
            &["counts", "record_entries_on_reference_skipped_records"]
        ),
        1
    );
    let without = compare_beside(
        &bare,
        &vcf_expected,
        &vcf_actual,
        &mut documented_with(&[SKIPPED]),
    );
    assert_eq!(number(&without, &["one_sided", "actual_only"]), 3);

    let object = |record: usize, allele: &str, features: &[&str]| -> String {
        let rows: Vec<String> = features
            .iter()
            .map(|feature| {
                format!(
                    r#"{{"variant_allele": "{allele}", "transcript_id": "{feature}", "consequence_terms": ["intron_variant"]}}"#
                )
            })
            .collect();
        format!(
            r#"{{"input": {}, "id": "{}", "most_severe_consequence": "intron_variant", "transcript_consequences": [{}]}}"#,
            pyjson::dumps_str(records[record], true),
            records[record].split('\t').nth(2).expect("an ID"),
            rows.join(", ")
        ) + "\n"
    };
    let json_expected =
        parse_json(&object(0, "deletion", &["ENST1"]), &inputs).expect("the reference parses");
    let json_actual = parse_json(
        &(object(0, "deletion", &["ENST1"]) + &object(1, "CPX", &["ENST1", "ENST1b"])),
        &inputs,
    )
    .expect("the candidate parses");
    let mut documented = documented_with(&[SKIPPED]);
    let report = compare_beside(&side, &json_expected, &json_actual, &mut documented);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(SKIPPED, 2)]));
    assert_eq!(
        number(
            &report,
            &["counts", "record_entries_on_reference_skipped_records"]
        ),
        1
    );
    // The class reads the input VCF in every format; the command line says so.
    assert!(refusal(&[SKIPPED], None, None).contains("needs --input"));
    let expected = planted(&directory, "ref.txt", &reference);
    let actual = planted(&directory, "cand.txt", &candidate);
    let argv = [
        "--format",
        "default",
        "--expected",
        path_text(&expected),
        "--actual",
        path_text(&actual),
        "--masks",
        SKIPPED,
    ];
    let refused = run(&argv, false);
    assert_eq!(refused.code, 2);
    assert_eq!(
        refused.err,
        "ERROR: [compare_vep_formats] --masks reference_skipped_record needs --input, the input VCF\n"
    );
    let report_path = directory.join("report.json");
    let mut args = argv.to_vec();
    args.extend([
        "--input",
        path_text(&vcf_input),
        "--report",
        path_text(&report_path),
    ]);
    let ran = run(&args, false);
    assert_eq!(ran.code, 0, "{}", ran.err);
    let report = report_at(&report_path);
    assert_eq!(
        at(&report, &["masks", "applied"]),
        &json(r#"{"reference_skipped_record": 2}"#)
    );
    assert_eq!(number(&report, &["json_cache_transcripts"]), 0);
}

/// Ensembl VEP evaluates a breakend's mate-side Transcript row at the record's own
/// coordinate (BaseVariationFeatureOverlapAllele.pm:257,273), the candidate at the mate's: a
/// bracket-mate pair whose two consequence sets are the local-coordinate and the
/// mate-coordinate derivations over the cache's transcript model is skipped whole, its
/// DISTANCE with it, counted once per row in every format; a candidate set other than the
/// mate derivation stays a difference, and the class reads the candidate JSON cache.
#[test]
fn a_mate_side_row_read_at_the_local_coordinate_is_set_aside_as_a_pair() {
    let directory = scratch("structural-mate-read");
    let vcf_input = structural_input(
        &directory.join("in.vcf"),
        &["21\t100000\tbnd1\tN\tN[22:5000000[\t.\t.\tSVTYPE=BND"],
    );
    let cache = breakend_cache(&directory.join("cache"));
    let mate = "N[22:5000000[";
    let own = default_row(
        "bnd1",
        "21:100001",
        "N.",
        "ENST21",
        "upstream_gene_variant",
        "IMPACT=MODIFIER;DISTANCE=499;STRAND=1",
    );
    let local_read = default_row(
        "bnd1",
        "21:100001",
        mate,
        "ENST22",
        "intergenic_variant",
        "IMPACT=MODIFIER;STRAND=1",
    );
    let mate_read = default_row(
        "bnd1",
        "21:100001",
        mate,
        "ENST22",
        "upstream_gene_variant",
        "IMPACT=MODIFIER;DISTANCE=100;STRAND=1",
    );
    let reference = default_text(&[&own, &local_read]);
    let candidate = default_text(&[&own, &mate_read]);
    let read = structural(&[MATE_READ], Some(&vcf_input), Some(&cache));
    let ids =
        |ids: &[&str]| -> BTreeSet<StableId> { ids.iter().map(|id| StableId::text(id)).collect() };
    assert_eq!(
        read.by_chr,
        IndexMap::from([
            ("21".to_string(), ids(&["ENST21", "ENST21_FAR"])),
            ("22".to_string(), ids(&["ENST22"])),
        ])
    );
    let side = side_of(read);
    let bare = Side::default();
    let (report, documented) = compare_structural(&reference, &candidate, &side, &[MATE_READ]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(MATE_READ, 1)]));
    assert_eq!(number(&report, &["counts", "registry_pairs_set_aside"]), 1);
    assert_eq!(number(&report, &["counts", "compared"]), 1);
    assert_eq!(
        examples(&report, MATE_READ),
        ["21:100001|N[22:5000000[|ENST22|Transcript: expected 'intergenic_variant', got 'upstream_gene_variant'"]
    );
    assert_eq!(
        verdict(&compare_structural(&reference, &candidate, &bare, &[MATE_READ]).0),
        "red"
    );
    let other_term = candidate.replace(
        "upstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=100",
        "downstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=100",
    );
    let (red, documented) = compare_structural(&reference, &other_term, &side, &[MATE_READ]);
    assert_eq!(verdict(&red), "red");
    assert_eq!(applied(&documented), count_map(&[]));
    assert!(refusal(&[MATE_READ], Some(&vcf_input), None).contains("needs --json-cache"));
    assert!(refusal(
        &[MATE_READ],
        Some(&vcf_input),
        Some(&directory.join("none"))
    )
    .contains("no transcripts/ directory"));
    // The VCF and JSON formats place the row at the record's Location through the input record.
    let inputs = whole_records(&vcf_input);
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. \
                      Format: Allele|Consequence|Feature_type|Feature|DISTANCE\">\n";
    let vcf = |mate_terms: &str, distance: &str| -> Vec<Entry> {
        let text = format!(
            "{csq_header}{};CSQ=N.|upstream_gene_variant|Transcript|ENST21|499,{mate}|{mate_terms}|Transcript|ENST22|{distance}\n",
            inputs[0]
        );
        parse_vcf(&text, &inputs).expect("the output parses")
    };
    let mut documented = documented_with(&[MATE_READ]);
    let report = compare_beside(
        &side,
        &vcf("intergenic_variant", ""),
        &vcf("upstream_gene_variant", "100"),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(MATE_READ, 1)]));

    let object = |severe: &str, mate_terms: &str, more: &str| -> Vec<Entry> {
        let text = format!(
            r#"{{"input": {}, "id": "bnd1", "most_severe_consequence": "{severe}", "transcript_consequences": [{{"variant_allele": "N.", "transcript_id": "ENST21", "consequence_terms": ["upstream_gene_variant"], "distance": 499}}, {{"variant_allele": "{mate}", "transcript_id": "ENST22", "consequence_terms": ["{mate_terms}"]{more}}}]}}"#,
            pyjson::dumps_str(&inputs[0], true)
        ) + "\n";
        parse_json(&text, &inputs).expect("the output parses")
    };
    let mut documented = documented_with(&[MATE_READ]);
    let report = compare_beside(
        &side,
        &object("upstream_gene_variant", "intergenic_variant", ""),
        &object(
            "upstream_gene_variant",
            "upstream_gene_variant",
            r#", "distance": 100"#,
        ),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(MATE_READ, 1)]));
    let expected = planted(&directory, "ref.txt", &reference);
    let actual = planted(&directory, "cand.txt", &candidate);
    let argv = [
        "--format",
        "default",
        "--expected",
        path_text(&expected),
        "--actual",
        path_text(&actual),
        "--masks",
        MATE_READ,
    ];
    let refused = run(&argv, false);
    assert_eq!(refused.code, 2);
    assert_eq!(
        refused.err,
        "ERROR: [compare_vep_formats] --masks breakend_mate_local_read needs --json-cache, the candidate JSON cache\n"
    );
    let report_path = directory.join("report.json");
    let mut args = argv.to_vec();
    args.extend([
        "--json-cache",
        path_text(&cache),
        "--report",
        path_text(&report_path),
    ]);
    let ran = run(&args, false);
    assert_eq!(ran.code, 0, "{}", ran.err);
    assert_eq!(
        number(&report_at(&report_path), &["json_cache_transcripts"]),
        3
    );
}

/// A mate spelled `CHR22` names no chromosome to the reference (StructuralVariationOverlap.pm
/// :86-103), so it writes no mate row, and where the own end reaches no transcript it writes
/// no CSQ at all, as for `chr13 31613634 C C[CHR17:7675000[`. A candidate-only mate row of
/// a transcript on the mate's chromosome within 5,000 bases of the mate is the class, and the
/// record's CSQ_present is left uncompared; a mate the reference resolves, a transcript
/// farther from the mate stay differences.
#[test]
fn the_mate_rows_of_a_record_the_reference_wrote_no_csq_for_are_set_aside() {
    let directory = scratch("structural-synonym-mate");
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. \
                      Format: Allele|Consequence|Feature_type|Feature|DISTANCE\">\n";
    let cache = breakend_cache(&directory.join("cache"));
    let judged = |mate: &str| -> (String, std::collections::BTreeMap<String, i64>) {
        let record = format!("21\t150000\tbnd3\tN\tN{mate}\t.\t.\tSVTYPE=BND");
        let vcf_input = structural_input(&directory.join("in.vcf"), &[record.as_str()]);
        let inputs = whole_records(&vcf_input);
        let side = side_of(structural(&[OWN_ROWS], Some(&vcf_input), Some(&cache)));
        let line = &inputs[0];
        let mut documented = documented_with(&[OWN_ROWS]);
        let report = compare_beside(
            &side,
            &parse_vcf(&format!("{csq_header}{line}\n"), &inputs).expect("the reference parses"),
            &parse_vcf(
                &format!(
                    "{csq_header}{line};CSQ=N{mate}|upstream_gene_variant|Transcript|ENST22|100\n"
                ),
                &inputs,
            )
            .expect("the candidate parses"),
            &mut documented,
        );
        (verdict(&report).to_string(), applied(&documented))
    };
    assert_eq!(
        judged("[CHR22:5000000["),
        ("green".to_string(), count_map(&[(OWN_ROWS, 1)]))
    );
    assert_eq!(judged("[chr22:5000000[").0, "red");
    assert_eq!(judged("[CHR22:4990000[").0, "red");
}

/// A mate-chromosome transcript names the record's slice first (AnnotationType/Transcript.pm
/// :122) and every own-chromosome transcript fails the slice test (:125-143), so the
/// reference keeps the mate rows and writes nothing under the own-end allele: a
/// candidate-only Transcript row under an own-end allele at such a Location, its transcript
/// on the Location's chromosome within 5,000 bases of the breakpoint by the cache's spans, is
/// the class, counted once per row in every format and leaving the JSON record's
/// most_severe_consequence uncompared; a farther transcript, a Location where the reference
/// kept an own-end row and a same-chromosome mate stay differences.
#[test]
fn the_own_end_rows_the_reference_lost_are_set_aside() {
    let directory = scratch("structural-own-rows");
    let vcf_input = structural_input(
        &directory.join("in.vcf"),
        &["21\t100000\tbnd1\tN\tN[22:5000000[\t.\t.\tSVTYPE=BND"],
    );
    let cache = breakend_cache(&directory.join("cache"));
    let mate = "N[22:5000000[";
    let mate_row = default_row(
        "bnd1",
        "21:100001",
        mate,
        "ENST22",
        "upstream_gene_variant",
        "IMPACT=MODIFIER;DISTANCE=100;STRAND=1",
    );
    let own = default_row(
        "bnd1",
        "21:100001",
        "N.",
        "ENST21",
        "upstream_gene_variant",
        "IMPACT=MODIFIER;DISTANCE=499;STRAND=1",
    );
    let far = default_row(
        "bnd1",
        "21:100001",
        "N.",
        "ENST21_FAR",
        "upstream_gene_variant",
        "IMPACT=MODIFIER;DISTANCE=99999;STRAND=1",
    );
    let reference = default_text(&[&mate_row]);
    let candidate = default_text(&[&mate_row, &own]);
    let side = side_of(structural(&[OWN_ROWS], None, Some(&cache)));
    let bare = Side::default();
    let (report, documented) = compare_structural(&reference, &candidate, &side, &[OWN_ROWS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(OWN_ROWS, 1)]));
    assert_eq!(
        examples(&report, OWN_ROWS),
        ["21:100001|N.|ENST21|Transcript: candidate-only own-end row, 'upstream_gene_variant'"]
    );
    assert_eq!(
        verdict(&compare_structural(&reference, &candidate, &bare, &[OWN_ROWS]).0),
        "red"
    );
    let (red, documented) = compare_structural(
        &reference,
        &default_text(&[&mate_row, &own, &far]),
        &side,
        &[OWN_ROWS],
    );
    assert_eq!(verdict(&red), "red");
    assert_eq!(number(&red, &["one_sided", "actual_only"]), 1);
    assert_eq!(applied(&documented), count_map(&[(OWN_ROWS, 1)]));
    let (red, documented) = compare_structural(
        &default_text(&[&mate_row, &far]),
        &default_text(&[&mate_row, &far, &own]),
        &side,
        &[OWN_ROWS],
    );
    assert_eq!(verdict(&red), "red");
    assert_eq!(applied(&documented), count_map(&[]));
    let same_chromosome = |text: &str| text.replace("N[22:5000000[", "N[21:5000000[");
    assert_eq!(
        verdict(
            &compare_structural(
                &same_chromosome(&reference),
                &same_chromosome(&candidate),
                &side,
                &[OWN_ROWS]
            )
            .0
        ),
        "red"
    );
    // The VCF and JSON formats place the rows through the input record; the JSON record's
    // most_severe_consequence follows the lost row and is left uncompared.
    let inputs = whole_records(&vcf_input);
    let side = side_of(structural(&[OWN_ROWS], Some(&vcf_input), Some(&cache)));
    let csq_header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. \
                      Format: Allele|Consequence|Feature_type|Feature|DISTANCE\">\n";
    let mate_csq = format!("{mate}|upstream_gene_variant|Transcript|ENST22|100");
    let mut documented = documented_with(&[OWN_ROWS]);
    let report = compare_beside(
        &side,
        &parse_vcf(
            &format!("{csq_header}{};CSQ={mate_csq}\n", inputs[0]),
            &inputs,
        )
        .expect("the reference parses"),
        &parse_vcf(
            &format!(
                "{csq_header}{};CSQ={mate_csq},N.|upstream_gene_variant|Transcript|ENST21|499\n",
                inputs[0]
            ),
            &inputs,
        )
        .expect("the candidate parses"),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(OWN_ROWS, 1)]));

    let mate_object = format!(
        r#"{{"variant_allele": "{mate}", "transcript_id": "ENST22", "consequence_terms": ["upstream_gene_variant"], "distance": 100}}"#
    );
    let own_object = r#"{"variant_allele": "N.", "transcript_id": "ENST21", "consequence_terms": ["intron_variant"]}"#;
    let object = |severe: &str, rows: &[&str]| -> Vec<Entry> {
        let text = format!(
            r#"{{"input": {}, "id": "bnd1", "most_severe_consequence": "{severe}", "transcript_consequences": [{}]}}"#,
            pyjson::dumps_str(&inputs[0], true),
            rows.join(", ")
        ) + "\n";
        parse_json(&text, &inputs).expect("the output parses")
    };
    let json_expected = object("upstream_gene_variant", &[&mate_object]);
    let json_actual = object("intron_variant", &[&mate_object, own_object]);
    let mut documented = documented_with(&[OWN_ROWS]);
    let report = compare_beside(&side, &json_expected, &json_actual, &mut documented);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(OWN_ROWS, 1)]));
    let without = compare_beside(
        &bare,
        &json_expected,
        &json_actual,
        &mut documented_with(&[OWN_ROWS]),
    );
    assert_eq!(verdict(&without), "red");
    assert_eq!(number(&without, &["one_sided", "actual_only"]), 1);
    assert_eq!(
        number(
            &without,
            &["columns", "most_severe_consequence", "differing"]
        ),
        1
    );
}

/// One transcript row of record `record`, at a Location of the record's own.
fn planted_row(record: usize, feature: &str) -> String {
    default_row(
        &format!("rec{record}"),
        &format!("21:{}", 1000 + record),
        &"ACGT"[record % 4..record % 4 + 1],
        feature,
        "intron_variant",
        &format!("IMPACT=MODIFIER;STRAND=1;DOMAINS=db:a{record},db:b{record}"),
    )
}

/// Two transcript rows for each of eight records.
fn planted_rows() -> Vec<String> {
    (0..8)
        .flat_map(|record| ["ENST1", "ENST2"].map(|feature| planted_row(record, feature)))
        .collect()
}

/// One comparison on the whole outputs and through the stream of record groups: the
/// standard output, the exit code, and what the stream did with the two outputs as one
/// line of rows kept back and groups handed over. The two runs' report,
/// standard output and exit code are the same bytes.
fn both_paths(
    directory: &Path,
    run_name: &str,
    fmt: &str,
    expected: &Path,
    actual: &Path,
    options: &[&str],
) -> (String, i32, String) {
    let mut seen: Vec<(Option<Vec<u8>>, String, i32, String)> = Vec::new();
    for (name, whole) in [("whole", true), ("handed", false)] {
        let report = directory.join(format!("{name}.{run_name}.json"));
        let mut args = vec![
            "--format",
            fmt,
            "--expected",
            path_text(expected),
            "--actual",
            path_text(actual),
            "--report",
            path_text(&report),
        ];
        args.extend(options);
        let ran = run(&args, whole);
        seen.push((std::fs::read(&report).ok(), ran.out, ran.code, ran.err));
    }
    let (handed, whole) = (seen.pop().expect("two runs"), seen.pop().expect("two runs"));
    assert_eq!(
        (&handed.0, &handed.1, handed.2),
        (&whole.0, &whole.1, whole.2),
        "{}{}",
        handed.3,
        whole.3
    );
    assert!(whole.3.is_empty() || whole.2 == 2, "{}", whole.3);
    let value_of = |name: &str| {
        options
            .iter()
            .position(|option| *option == name)
            .map(|at| options[at + 1])
    };
    let classes: Vec<&str> = value_of("--masks").into_iter().collect();
    let mut front_end = FrontEnd::new(STREAM);
    hand_over(
        &mut front_end,
        fmt,
        expected,
        actual,
        &documented_with(&classes),
        false,
        value_of("--keep-only"),
        &[],
        None,
    )
    .expect("the outputs are read");
    (
        handed.1,
        handed.2,
        format!("[compare_vep_formats] {}\n", front_end.note),
    )
}

/// A structural-variant class, `individual_allele_order`, `--keep-only` and a format the
/// stream does not read are compared on the whole outputs, and the stream's line names what
/// is outside it.
#[test]
fn the_comparisons_outside_the_stream_are_made_on_the_whole_outputs() {
    let directory = scratch("structural-outside-the-stream");
    let rows = planted_rows();
    let outside = |what: &str| {
        format!("[compare_vep_formats] the whole outputs were compared: {what} is outside the front end\n")
    };
    let input = corpus().join("variants.vcf");
    let rows_file = planted(&directory, "rows.txt", &default_text(&rows));
    assert_eq!(
        tables().structural_variant_classes(),
        [SKIPPED, MATE_READ, OWN_ROWS]
    );
    assert_eq!(tables().cache_reading_classes(), [MATE_READ, OWN_ROWS]);
    let mut classes = tables().structural_variant_classes().to_vec();
    classes.push(INDIVIDUAL_ALLELE_ORDER);
    for class in classes {
        if tables().in_cache_reading_classes(class) {
            continue;
        }
        let compared = both_paths(
            &directory,
            class,
            "default",
            &rows_file,
            &rows_file,
            &["--masks", class, "--input", path_text(&input)],
        );
        assert_eq!(
            compared,
            ("default: green\n".to_string(), 0, outside(class))
        );
    }
    let regulatory = default_row_of(
        "rec2",
        "21:1002",
        "G",
        "ENSR1",
        "regulatory_region_variant",
        "IMPACT=MODIFIER",
        "RegulatoryFeature",
        "ENSG1",
    );
    let mut with_regulatory = rows.clone();
    with_regulatory.push(regulatory);
    let regulatory_file = planted(
        &directory,
        "regulatory.txt",
        &default_text(&with_regulatory),
    );
    let compared = both_paths(
        &directory,
        "keep-only",
        "default",
        &regulatory_file,
        &regulatory_file,
        &["--keep-only", "regulatory_feature_rows"],
    );
    assert_eq!(
        compared,
        ("default: green\n".to_string(), 0, outside("--keep-only"))
    );
    let recoder = planted(
        &directory,
        "recoder.json",
        r#"[{"T": {"input": "rs1", "hgvsg": ["21:g.1A>T"]}}]"#,
    );
    let compared = both_paths(&directory, "recoder", "recoder", &recoder, &recoder, &[]);
    assert_eq!(
        compared,
        (
            "recoder: green\n".to_string(),
            0,
            outside("the recoder format")
        )
    );
    let warnings = planted(&directory, "warnings.txt", "");
    let ran = run(
        &[
            "--warnings",
            "--expected",
            path_text(&warnings),
            "--actual",
            path_text(&warnings),
        ],
        false,
    );
    assert_eq!(ran.code, 0, "{}", ran.err);
    assert_eq!(
        FrontEnd::new(STREAM).note,
        "the whole outputs were compared: only a default, tab, vcf, json or parquet comparison is handed over"
    );
}
