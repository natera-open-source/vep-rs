// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the SpliceAI class: the members of a row's value, the JSON object of a
//! row, the entries of the other formats, the line a question is, the fields the class
//! takes of those the format comparator accepts, and the line a run ends on where the
//! command leaves no answers.

use super::*;
use crate::classify::fields::FieldKey;
use crate::TempDir;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;
use std::path::PathBuf;

const GENEA: &str = "GENEA|0.15|0.08|0.06|0.17|-6|-41|-17|-9";
const GENEB: &str = "GENEB|0.00|1.00|0.10|0.05|-25|-23|37|17";

fn row(fields: &[(&str, &str)]) -> FieldRow {
    fields
        .iter()
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect()
}

/// A row that writes a value as its nine member fields.
fn split(value: &str) -> FieldRow {
    SPLICEAI_MEMBER_FIELDS
        .iter()
        .zip(value.split('|'))
        .map(|(field, part)| (field.to_string(), part.to_string()))
        .collect()
}

fn texts(parts: &[&str]) -> Vec<String> {
    parts.iter().map(|part| part.to_string()).collect()
}

#[test]
fn the_fields_of_a_row_are_the_twelve_the_class_reads() {
    assert_eq!(SPLICEAI_ROW_FIELDS.len(), 12);
    assert_eq!(SPLICEAI_ROW_FIELDS[1..10], SPLICEAI_MEMBER_FIELDS);
    assert_eq!(SPLICEAI_MULTI_ENTRY_CELL, "spliceai_multi_entry_cell");
    assert_eq!(SPLICEAI_NUMBER, r"-?\d+(\.\d+)?");
}

#[test]
fn the_members_of_a_row_are_nine_parts_or_nine_fields() {
    let nine = texts(&["a", "b", "c", "d", "e", "f", "g", "h", "i"]);
    let whole = row(&[("SpliceAI_pred", "a|b|c|d|e|f|g|h|i|j|k")]);
    assert_eq!(spliceai_members(&whole), Some(nine));
    assert_eq!(spliceai_members(&row(&[("SpliceAI_pred", "a|b")])), None);
    let unwritten = row(&[("SpliceAI_pred", "a|b|c|d|-|f|g|h|i")]);
    assert_eq!(spliceai_members(&unwritten), None);
    assert_eq!(spliceai_members(&row(&[("SpliceAI_pred", "")])), None);
    assert_eq!(spliceai_members(&FieldRow::new()), None);
    let mut fields = split("v|v|v|v|v|v|v|v|v");
    assert_eq!(spliceai_members(&fields), Some(texts(&["v"; 9])));
    fields.insert("SpliceAI_pred".to_string(), "-".to_string());
    assert_eq!(spliceai_members(&fields), Some(texts(&["v"; 9])));
    fields.remove("SpliceAI_pred_DP_DL");
    assert_eq!(spliceai_members(&fields), None);
}

#[test]
fn the_json_object_of_a_row_writes_a_number_by_its_value() {
    let object = |text: &str| {
        let members = texts(&["GENE", text, "1", "1", "1", "1", "1", "1", "1"]);
        let written = spliceai_json_object(Some(&members), "-");
        let rest = written
            .strip_prefix(r#"{"DP_AG":1,"DP_AL":1,"DP_DG":1,"DP_DL":1,"DS_AG":"#)
            .unwrap();
        rest.strip_suffix(r#","DS_AL":1,"DS_DG":1,"DS_DL":1,"SYMBOL":"GENE"}"#)
            .unwrap()
            .to_string()
    };
    for (text, written) in [
        ("0.00", "0"),
        ("1.00", "1"),
        ("0.10", "0.1"),
        ("-0.0", "0"),
        ("-0", "0"),
        ("007", "7"),
        ("-17", "-17"),
        ("0.30000000000000004", "0.30000000000000004"),
        ("123456789012345678", "123456789012345680"),
        ("12345678901234567890123", "12345678901234567741440"),
        ("0.000001", "1e-06"),
        ("0.0000001", "1e-07"),
        // What is no number is a text.
        ("1e3", r#""1e3""#),
        ("-9,G", r#""-9,G""#),
        ("0.1.2", r#""0.1.2""#),
        (".5", r#"".5""#),
        ("5.", r#""5.""#),
        ("+5", r#""+5""#),
        (" 5", r#"" 5""#),
        ("1_0", r#""1_0""#),
    ] {
        assert_eq!(object(text), written, "{text}");
    }
    assert_eq!(object(&format!("1{}", "0".repeat(400))), "Infinity");
    assert_eq!(object(&format!("-{}", "9".repeat(400))), "-Infinity");
    assert_eq!(spliceai_json_object(None, "PASS"), "-");
    // The symbol is a text whatever it holds, a text is written in ASCII, and an empty
    // flag is a flag.
    let members = texts(&["GEN\u{c9}", "x\"y", "1", "1", "1", "1", "1", "1", "1"]);
    let escaped = format!("GEN{}u00c9", '\\');
    assert_eq!(
        spliceai_json_object(Some(&members), "-"),
        format!(
            r#"{{"DP_AG":1,"DP_AL":1,"DP_DG":1,"DP_DL":1,"DS_AG":"x\"y","DS_AL":1,"DS_DG":1,"DS_DL":1,"SYMBOL":"{escaped}"}}"#
        )
    );
    let members = texts(&["0.5", "x", "1", "1", "1", "1", "1", "1", "1"]);
    assert_eq!(
        spliceai_json_object(Some(&members), ""),
        r#"{"DP_AG":1,"DP_AL":1,"DP_DG":1,"DP_DL":1,"DS_AG":"x","DS_AL":1,"DS_DG":1,"DS_DL":1,"SYMBOL":"0.5","SpliceAI_cutoff":""}"#
    );
}

/// The entry `classify` writes for a field of a row of the class.
fn written(feature: &str) -> Object {
    let mut entry = Object::new();
    entry.insert("location", "21:5001-5002");
    entry.insert("allele", "-");
    entry.insert("feature", feature);
    entry.insert("feature_type", "Transcript");
    entry.insert("record_indices", vec![Value::from(0usize)]);
    entry.insert("consequence_set", "intron_variant");
    entry.insert("field", "SpliceAI_pred");
    entry.insert("vep", "x");
    entry.insert("vep_rs", "y");
    entry.insert("expected_divergence", SPLICEAI_MULTI_ENTRY_CELL);
    entry
}

/// The field and the two values of each entry, with the place of `field` among its keys.
fn derived(entries: &[Object]) -> Vec<(String, String, String)> {
    let text =
        |entry: &Object, key: &str| entry.get(key).and_then(Value::as_str).unwrap().to_string();
    entries
        .iter()
        .map(|entry| {
            let keys: Vec<&str> = entry.keys().collect();
            assert_eq!(keys, written("ENST1").keys().collect::<Vec<_>>());
            assert_eq!(
                text(entry, "expected_divergence"),
                SPLICEAI_MULTI_ENTRY_CELL
            );
            (
                text(entry, "field"),
                text(entry, "vep"),
                text(entry, "vep_rs"),
            )
        })
        .collect()
}

#[test]
fn classify_writes_the_member_fields_of_a_row_of_the_spliceai_class_for_the_vcf_format() {
    let whole = format!("{GENEA},G|{GENEB}");
    let compared = ["SpliceAI_pred".to_string()];
    // The first entry's gene's row: Ensembl VEP's ninth part runs on into the second
    // entry's allele.
    let (perl, rust) = (
        row(&[("SpliceAI_pred", &whole)]),
        row(&[("SpliceAI_pred", GENEA)]),
    );
    let entries = derived(&spliceai_other_format_entries(
        &written("ENST1"),
        &perl,
        &rust,
        &compared,
    ));
    let own = r#"{"DP_AG":-6,"DP_AL":-41,"DP_DG":-17,"DP_DL":-9,"DS_AG":0.15,"DS_AL":0.08,"DS_DG":0.06,"DS_DL":0.17,"SYMBOL":"GENEA"}"#;
    let run_on = own.replace(r#""DP_DL":-9"#, r#""DP_DL":"-9,G""#);
    let entry = |field: &str, vep: &str, vep_rs: &str| {
        (field.to_string(), vep.to_string(), vep_rs.to_string())
    };
    assert_eq!(
        entries,
        [
            entry("SpliceAI_pred_DP_DL", "-9,G", "-9"),
            entry("spliceai", &run_on, own)
        ]
    );
    // The second entry's gene's row: all nine members, then the object.
    let (perl, rust) = (
        row(&[("SpliceAI_pred", "-")]),
        row(&[("SpliceAI_pred", GENEB)]),
    );
    let entries = derived(&spliceai_other_format_entries(
        &written("ENST2"),
        &perl,
        &rust,
        &compared,
    ));
    let mut members: Vec<_> = SPLICEAI_MEMBER_FIELDS
        .iter()
        .zip(GENEB.split('|'))
        .map(|(field, part)| entry(field, "-", part))
        .collect();
    members.push(entry(
        "spliceai",
        "-",
        r#"{"DP_AG":-25,"DP_AL":-23,"DP_DG":37,"DP_DL":17,"DS_AG":0,"DS_AL":1,"DS_DG":0.1,"DS_DL":0.05,"SYMBOL":"GENEB"}"#,
    ));
    assert_eq!(entries, members);
    // A run that compares the nine fields itself gains the object alone.
    let compared: Vec<String> = SPLICEAI_MEMBER_FIELDS
        .iter()
        .map(|field| field.to_string())
        .collect();
    let (perl, rust) = (split(&whole), split(GENEA));
    let entries = derived(&spliceai_other_format_entries(
        &written("ENST1"),
        &perl,
        &rust,
        &compared,
    ));
    assert_eq!(entries, [entry("spliceai", &run_on, own)]);
}

#[test]
fn the_members_of_a_vep_rs_row_of_several_entries_are_their_values_joined() {
    let other = "GENEA|0.31|0.00|0.44|0.09|3|-2|11|-30";
    let joined = row(&[("SpliceAI_pred", &format!("{GENEA},G|{other},G|{GENEA}"))]);
    assert_eq!(
        spliceai_entry_members(&joined),
        Some(texts(&[
            "GENEA,GENEA,GENEA",
            "0.15,0.31,0.15",
            "0.08,0.00,0.08",
            "0.06,0.44,0.06",
            "0.17,0.09,0.17",
            "-6,3,-6",
            "-41,-2,-41",
            "-17,11,-17",
            "-9,-30,-9",
        ]))
    );
    // One entry, the nine fields and a row without the value are read as
    // `spliceai_members` reads them; an entry of fewer than nine parts is no whole value.
    let single = row(&[("SpliceAI_pred", GENEA)]);
    assert_eq!(spliceai_entry_members(&single), spliceai_members(&single));
    let fields = split(GENEB);
    assert_eq!(spliceai_entry_members(&fields), spliceai_members(&fields));
    assert_eq!(
        spliceai_entry_members(&row(&[("SpliceAI_pred", "-")])),
        None
    );
    let short = row(&[("SpliceAI_pred", &format!("{GENEA},G|GENEA|0.31"))]);
    assert_eq!(spliceai_entry_members(&short), None);
}

#[test]
fn the_member_question_is_the_first_member_that_differs_on_rows_written_alike() {
    let other = "GENEA|0.31|0.00|0.44|0.09|3|-2|11|-30";
    let value = format!("{GENEA},G|{other}");
    let alike = row(&[("SpliceAI_pred", &value)]);
    assert_eq!(
        member_question(&alike, &alike),
        Some(MemberQuestion {
            name: "SpliceAI_pred_SYMBOL",
            vep: "GENEA".to_string(),
            vep_rs: "GENEA,GENEA".to_string(),
        })
    );
    // One entry gives both engines the same nine members; rows that differ in the one
    // field are a question of that field; a row without the value asks nothing.
    let single = row(&[("SpliceAI_pred", GENEA)]);
    assert_eq!(member_question(&single, &single), None);
    assert_eq!(member_question(&alike, &single), None);
    let none = row(&[("SpliceAI_pred", "-")]);
    assert_eq!(member_question(&none, &none), None);
    assert_eq!(member_question(&FieldRow::new(), &alike), None);
}

#[test]
fn classify_writes_the_other_member_fields_behind_the_entry_of_a_member_question() {
    let other = "GENEA|0.31|0.00|0.44|0.09|3|-2|11|-30";
    let value = format!("{GENEA},G|{other}");
    let alike = row(&[("SpliceAI_pred", &value)]);
    let mut asked = written("ENST1");
    asked.insert("field", "SpliceAI_pred_SYMBOL");
    let entries = derived(&spliceai_other_format_entries(
        &asked,
        &alike,
        &alike,
        &["SpliceAI_pred".to_string()],
    ));
    let fields: Vec<&str> = entries.iter().map(|entry| entry.0.as_str()).collect();
    let mut expected: Vec<&str> = SPLICEAI_MEMBER_FIELDS[1..].to_vec();
    expected.push("spliceai");
    assert_eq!(fields, expected);
    assert_eq!(
        entries[7],
        (
            "SpliceAI_pred_DP_DL".to_string(),
            "-9,G".to_string(),
            "-9,-30".to_string()
        )
    );
    assert_eq!(
        entries[8].2,
        r#"{"DP_AG":"-6,3","DP_AL":"-41,-2","DP_DG":"-17,11","DP_DL":"-9,-30","DS_AG":"0.15,0.31","DS_AL":"0.08,0.00","DS_DG":"0.06,0.44","DS_DL":"0.17,0.09","SYMBOL":"GENEA,GENEA"}"#
    );
}

#[test]
fn classify_writes_the_json_object_of_a_row_of_the_spliceai_class_as_the_harness_prints_one() {
    let whole = format!("{GENEA},G|{GENEB}");
    let compared = ["SpliceAI_pred".to_string(), "SpliceAI_cutoff".to_string()];
    let flagged =
        |value: &str, flag: &str| row(&[("SpliceAI_pred", value), ("SpliceAI_cutoff", flag)]);
    let object = |entries: Vec<Object>| derived(&entries).pop().unwrap();
    let own = r#"{"DP_AG":-6,"DP_AL":-41,"DP_DG":-17,"DP_DL":-9,"DS_AG":0.15,"DS_AL":0.08,"DS_DG":0.06,"DS_DL":0.17,"SYMBOL":"GENEA","SpliceAI_cutoff":"FAIL"}"#;
    let first = spliceai_other_format_entries(
        &written("ENST1"),
        &flagged(&whole, "FAIL"),
        &flagged(GENEA, "FAIL"),
        &compared,
    );
    assert_eq!(
        object(first),
        (
            "spliceai".to_string(),
            own.replace(r#""DP_DL":-9"#, r#""DP_DL":"-9,G""#),
            own.to_string()
        )
    );
    let second = spliceai_other_format_entries(
        &written("ENST2"),
        &flagged("-", "-"),
        &flagged(GENEB, "PASS"),
        &compared,
    );
    assert_eq!(
        object(second),
        (
            "spliceai".to_string(),
            "-".to_string(),
            r#"{"DP_AG":-25,"DP_AL":-23,"DP_DG":37,"DP_DL":17,"DS_AG":0,"DS_AL":1,"DS_DG":0.1,"DS_DL":0.05,"SYMBOL":"GENEB","SpliceAI_cutoff":"PASS"}"#.to_string()
        )
    );
}

fn key(feature: &str) -> FieldKey {
    FieldKey {
        location: "21:5001-5002".to_string(),
        allele: "-".to_string(),
        feature: feature.to_string(),
        feature_type: "Transcript".to_string(),
        consequences: "intron_variant".to_string(),
    }
}

fn question<'a>(
    name: &'a str,
    key: &'a FieldKey,
    values: (&'a str, &'a str),
    rows: (Option<&'a FieldRow>, Option<&'a FieldRow>),
) -> FieldQuestion<'a> {
    FieldQuestion {
        name,
        key,
        vep: values.0,
        vep_rs: values.1,
        strand: "1",
        flags: &[],
        perl_row: rows.0,
        rust_row: rows.1,
        minimal_split_record: false,
        genomic_span: None,
    }
}

#[test]
fn a_question_is_a_line_of_thirty_columns_with_a_dash_for_a_field_a_row_lacks() {
    let key = key("ENST2");
    let perl_row = row(&[
        ("SpliceAI_pred", "-"),
        ("SpliceAI_cutoff", ""),
        ("SYMBOL", "GENEB"),
        ("STRAND", "1"),
    ]);
    let rust_row = row(&[
        ("SYMBOL", "GENEB"),
        ("SpliceAI_cutoff", "PASS"),
        ("SpliceAI_pred", GENEB),
        ("HGVSc", "ENST2.1:c.5del"),
    ]);
    let asked = question(
        "SpliceAI_pred",
        &key,
        ("-", GENEB),
        (Some(&perl_row), Some(&rust_row)),
    );
    let line = question_line(&asked);
    let cells: Vec<&str> = line.strip_suffix('\n').unwrap().split('\t').collect();
    assert_eq!(cells.len(), 6 + 2 * SPLICEAI_ROW_FIELDS.len());
    assert_eq!(cells.len(), 30);
    assert_eq!(
        cells[..6],
        ["21:5001-5002", "-", "ENST2", "SpliceAI_pred", "-", GENEB]
    );
    // The reference's row: a value of `-` and a field the row lacks are one cell, an empty
    // value is its own, and a field the class does not read is no cell.
    let mut reference = ["-"; 12];
    (reference[10], reference[11]) = ("", "GENEB");
    assert_eq!(cells[6..18], reference);
    let mut candidate = ["-"; 12];
    (candidate[0], candidate[10], candidate[11]) = (GENEB, "PASS", "GENEB");
    assert_eq!(cells[18..], candidate);
    // A question without a row carries a dash for each of its fields, and the member
    // fields stand in their order behind the whole value.
    let members = split(GENEA);
    let asked = question(
        "SpliceAI_pred_DS_AG",
        &key,
        ("-", "0.15"),
        (None, Some(&members)),
    );
    let line = question_line(&asked);
    let cells: Vec<&str> = line.strip_suffix('\n').unwrap().split('\t').collect();
    assert_eq!(cells[6..18], ["-"; 12]);
    let parts: Vec<&str> = GENEA.split('|').collect();
    assert_eq!(cells[18], "-");
    assert_eq!(cells[19..28], parts[..]);
    assert_eq!(cells[28..], ["-", "-"]);
}

#[test]
fn the_class_takes_the_fields_the_format_comparator_accepts() {
    let manifest = Object::new();
    let run = FieldRun {
        corpus: "c",
        cache_dir: "c/json_cache",
        manifest: &manifest,
        fields: &["SpliceAI_pred".to_string()],
        spliceai_tables: &[],
        reference_fasta: None,
    };
    // A run without score files reads no field of a row beside the compared ones and
    // takes no field.
    let models = Models::of(&run);
    assert!(!models.score_files && models.accepted.is_empty());
    assert_eq!(models.row_fields(), &[] as &[&str]);
    let (first, second) = (key("ENST1"), key("ENST2"));
    let rust_row = row(&[("SpliceAI_pred", GENEB), ("SYMBOL", "GENEB")]);
    let asked = |name: &'static str, key: &FieldKey, models: Models| {
        let asked = question(name, key, ("-", GENEB), (None, Some(&rust_row)));
        let models = FieldModels {
            spliceai: models,
            ..FieldModels::default()
        };
        spliceai_multi_entry_cell(&asked, &models)
    };
    assert_eq!(asked("SpliceAI_pred", &second, models), Verdict::Pass);
    assert_eq!(
        asked("SpliceAI_pred", &second, Models::default()),
        Verdict::Pass
    );
    // A run with score files reads the twelve fields, and takes a field of a key the
    // comparator accepted and no other.
    let tables = ["c/plugin_data/scores.synthetic.vcf.gz".to_string()];
    let scored = || {
        Models::of(&FieldRun {
            spliceai_tables: &tables,
            ..run
        })
    };
    assert!(scored().score_files);
    assert_eq!(scored().row_fields(), &SPLICEAI_ROW_FIELDS);
    assert_eq!(asked("SpliceAI_pred", &second, scored()), Verdict::Pass);
    let accepting = || {
        let mut models = scored();
        let fields = models.accepted.entry(second.clone()).or_default();
        fields.insert("SpliceAI_pred".to_string());
        models
    };
    assert_eq!(
        asked("SpliceAI_pred", &second, accepting()),
        Verdict::Class(SPLICEAI_MULTI_ENTRY_CELL, None)
    );
    assert_eq!(
        asked("SpliceAI_cutoff", &second, accepting()),
        Verdict::Pass
    );
    assert_eq!(asked("SpliceAI_pred", &first, accepting()), Verdict::Pass);
}

fn ended(code: Option<i32>, error: &str) -> Ended {
    Ended {
        code,
        error: error.to_string(),
    }
}

#[test]
fn a_run_the_command_refuses_ends_on_its_detail_under_this_tools_name() {
    let scratch = TempDir::create(&std::env::temp_dir(), "golden_spliceai_failure_").unwrap();
    let table = scratch.path.join("snv.synthetic.vcf.gz");
    let table = table.to_str().unwrap().to_string();
    let other = scratch.path.join("indel.synthetic.vcf.gz");
    let tables = [other.to_str().unwrap().to_string(), table.clone()];
    let index = format!("{table}.tbi");
    let comparator = Comparator {
        program: PathBuf::from("tools/vep-formats-compare/target/release/vep-formats-compare"),
    };
    let line = |code: Option<i32>, error: &str| {
        let failure = failure(&comparator, &tables, &ended(code, error));
        assert_eq!(failure.code, 1);
        failure.line
    };
    let own = |detail: &str| format!("ERROR: [vep-golden-corpus] {detail}");
    let theirs = |detail: &str| format!("ERROR: [compare_vep_formats] {detail}\n");
    // An index the command cannot read: the reason is this tool's own reading of the file.
    let unread = theirs(&format!(
        "--spliceai-table: cannot read the tabix index {index}: No such file or directory"
    ));
    assert_eq!(
        line(Some(2), &unread),
        own(&format!(
            "--spliceai-table: {index} is not a readable tabix index ([Errno 2] No such file or directory: '{index}')"
        ))
    );
    std::fs::write(&index, b"no gzip").unwrap();
    assert_eq!(
        line(Some(2), &unread),
        own(&format!(
            "--spliceai-table: {index} is not a readable tabix index (Not a gzipped file (b'no'))"
        ))
    );
    // A gzip file in whose bytes the command reads no index.
    let mut gzip = GzEncoder::new(Vec::new(), Compression::default());
    gzip.write_all(b"TBI\x02").unwrap();
    std::fs::write(&index, gzip.finish().unwrap()).unwrap();
    assert_eq!(
        line(Some(2), &unread),
        own(&format!(
            "--spliceai-table: {index} is not a readable tabix index (not a tabix index)"
        ))
    );
    // A block of a score file that does not inflate.
    let block = theirs(&format!(
        "--spliceai-table: cannot read {table}: Error -3 while decompressing data: invalid block type"
    ));
    assert_eq!(
        line(Some(2), &block),
        own(&format!(
            "--spliceai-table: {table} holds a block that does not inflate"
        ))
    );
    // Any other refusal of the command is the run's, word for word, the last line of what
    // the command wrote.
    for detail in [
        format!("--spliceai-table: {table} is not a BGZF file"),
        format!("--spliceai-table: {table} is not a VCF (a header line of eight columns)"),
        "c/json_cache/transcripts/21/1-2.json is not JSON: Expecting value: line 1 column 1 (char 0)"
            .to_string(),
        "cannot read c/json_cache/transcripts/21/1-2.json.gz: Not a gzipped file (b'[1')"
            .to_string(),
        "--spliceai-table: cannot read the tabix index elsewhere.synthetic.vcf.gz.tbi: No such file or directory"
            .to_string(),
    ] {
        let written = format!("a line before\n{}", theirs(&detail));
        assert_eq!(line(Some(2), &written), own(&detail), "{detail}");
    }
    // A run the command ends in another way names the build and how it ended.
    let build =
        "the format comparator tools/vep-formats-compare/target/release/vep-formats-compare";
    assert_eq!(
        line(
            Some(1),
            &theirs("cannot use 'list' as a dict key (unhashable type: 'list')")
        ),
        own(&format!(
            "{build} ended with exit code 1: cannot use 'list' as a dict key (unhashable type: 'list')"
        ))
    );
    assert_eq!(
        line(
            Some(2),
            "usage: vep-formats-compare field-class [-h]\nvep-formats-compare field-class: error: unrecognized arguments: x\n"
        ),
        own(&format!(
            "{build} ended with exit code 2: vep-formats-compare field-class: error: unrecognized arguments: x"
        ))
    );
    assert_eq!(
        line(Some(3), &theirs("a line of another exit code")),
        own(&format!(
            "{build} ended with exit code 3: a line of another exit code"
        ))
    );
    assert_eq!(
        line(None, ""),
        own(&format!(
            "{build} was ended by a signal and wrote no error line"
        ))
    );
}
