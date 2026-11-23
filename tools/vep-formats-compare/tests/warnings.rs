// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The warnings mode on planted warnings files and skipped-variants files.

mod support;

use std::collections::BTreeSet;
use std::path::Path;
use std::path::PathBuf;
use support::member as at;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::warnings::compare_warnings;
use vep_formats_compare::warnings::parse_warning_lines;
use vep_formats_compare::warnings::without_forked_child_prefix;
use vep_formats_compare::warnings::without_plugin_failure;
use vep_formats_compare::warnings::LineNumber;
use vep_formats_compare::warnings::REGULATORY_NOTICE;

/// The file `name` of `directory`, holding `text`.
fn write(directory: &Path, name: &str, text: &str) -> PathBuf {
    let path = directory.join(name);
    std::fs::write(&path, text).expect("the file is written");
    path
}

fn json(text: &str) -> Value {
    pyjson::loads(text).expect("the test's JSON reads")
}

fn text(text: &str) -> Value {
    Value::Str(text.to_string())
}

fn compared(expected: &Path, actual: &Path, masks: &[&str]) -> Value {
    compare_warnings(expected, actual, None, None, masks).expect("the files are read")
}

fn descriptions(rows: &Value) -> Vec<&str> {
    let rows = rows.as_list().expect("the rows are a list");
    rows.iter()
        .map(|row| at(row, &["description"]).as_str().expect("a text"))
        .collect()
}

#[test]
fn warnings_agree_once_the_record_texts_and_the_paths_are_masked() {
    let scratch = support::scratch("warnings-masking");
    let expected_text = "WARNING: line 13 skipped (21 30165800 multi_ins_del N <INS>,<DEL> . . END=3016...): \
         Invalid allele string N/<INS>,<DEL> or possible parsing error\n\
         WARNING: line 16 skipped (21 30165850 del_without_end N <DEL> . . SVTYPE=DEL): start > end+1 : (START=30165851, END=30165850)\n\
         WARNING: Unable to open /data/ref/one.fa for reading\n";
    let expected = write(&scratch, "ref_warnings.txt", expected_text);
    let actual = write(
        &scratch,
        "cand_warnings.txt",
        "WARNING: line 16 skipped (21\t30165850\tdel_without_end\tN\t<DEL>): start > end+1 : (START=30165851, END=30165850)\n\
         WARNING: line 13 skipped (): Invalid allele string N/<INS>,<DEL> or possible parsing error\n\
         WARNING: Unable to open ./two.fa for reading\n",
    );
    let report = compared(&expected, &actual, &[]);
    assert_eq!(at(&report, &["verdict"]), &text("green"));
    assert_eq!(at(&report, &["warnings", "shared"]), &Value::Int(3));
    let entry = |number: Option<&str>, description: &str| {
        (number.map(LineNumber::of), description.to_string())
    };
    assert_eq!(
        parse_warning_lines(expected_text),
        BTreeSet::from([
            entry(
                Some("13"),
                "Invalid allele string N/<INS>,<DEL> or possible parsing error"
            ),
            entry(Some("16"), "start > end+1 : (START=30165851, END=30165850)"),
            entry(None, "Unable to open <path> for reading"),
        ])
    );
    std::fs::write(
        &actual,
        "WARNING: line 13 skipped (x): Invalid allele string N/<INS>,<DEL> or possible parsing error\n",
    )
    .expect("the file is written");
    let report = compared(&expected, &actual, &[]);
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    assert_eq!(
        at(&report, &["warnings", "expected_only"]),
        &json(
            r#"[{"line": 16, "description": "start > end+1 : (START=30165851, END=30165850)"},
                {"line": null, "description": "Unable to open <path> for reading"}]"#
        )
    );
    assert_eq!(
        at(&report, &["reasons"]),
        &json(r#"["warnings: 2 entries only in the expected file"]"#)
    );
}

#[test]
fn skipped_variants_files_compare_by_line_and_description() {
    let scratch = support::scratch("warnings-skipped");
    let expected = write(
        &scratch,
        "ref_skipped.txt",
        "[VEP skipped the following variants from /data/in/variants.vcf]\n\
         Line 13   \t21\t30165800\tmulti_ins_del\tN\t<INS>,<DEL>\tInvalid allele string N/<INS>,<DEL> or possible parsing error\n\
         Line 16   \t21 30165850 del_without_end N <DEL> . . SVTYPE=DEL   \tstart > end+1 : (START=30165851, END=30165850)\n",
    );
    let actual = write(
        &scratch,
        "cand_skipped.txt",
        "[VEP skipped the following variants from variants.vcf]\n\
         Line 13\t21 30165800 multi_ins_del N <INS>,<DEL> . . END=30166000\tInvalid allele string N/<INS>,<DEL> or possible parsing error\n\
         Line 16\tx\tstart > end+1 : (START=30165851, END=30165850)\n",
    );
    let warnings = write(&scratch, "w.txt", "");
    let report = compare_warnings(&warnings, &warnings, Some(&expected), Some(&actual), &[])
        .expect("the files are read");
    assert_eq!(at(&report, &["verdict"]), &text("green"));
    assert_eq!(at(&report, &["skipped", "shared"]), &Value::Int(2));
    let one_sided = compare_warnings(&warnings, &warnings, Some(&expected), None, &[]);
    assert_eq!(
        one_sided
            .expect_err("one file of the pair is refused")
            .to_string(),
        "--expected-skipped and --actual-skipped go together"
    );
}

#[test]
fn the_regulatory_mask_takes_the_notice_of_a_cache_without_regulatory_data() {
    let scratch = support::scratch("warnings-cache-notice");
    let expected = write(&scratch, "expected_warnings.txt", "");
    let actual = write(
        &scratch,
        "actual_warnings.txt",
        "WARNING: --regulatory: this cache holds no regulatory data; no RegulatoryFeature or MotifFeature row is written\n",
    );
    assert_eq!(
        at(&compared(&expected, &actual, &[]), &["verdict"]),
        &text("red")
    );
    let report = compared(&expected, &actual, &["regulatory_feature_rows"]);
    assert_eq!(at(&report, &["verdict"]), &text("green"));
    assert_eq!(
        at(&report, &["warnings", "regulatory_notice_missing"]),
        &Value::Bool(false)
    );
    assert_eq!(at(&report, &["warnings", "actual_entries"]), &Value::Int(0));
}

/// The warnings-mode half of the `regulatory_feature_rows` mask; `tests/manifest.rs` has the rows.
#[test]
fn the_regulatory_mask_requires_the_candidates_notice_and_leaves_it_out() {
    let scratch = support::scratch("warnings-notice");
    let expected = write(
        &scratch,
        "ref_warnings.txt",
        "WARNING: line 2 skipped (21 1 x A G . . .): Invalid allele string\n",
    );
    let skipped = "WARNING: line 2 skipped (): Invalid allele string\n";
    let notice = format!(
        "WARNING: --everything implies --regulatory, which has no effect: {REGULATORY_NOTICE}\n"
    );
    let actual = write(&scratch, "cand_warnings.txt", &format!("{notice}{skipped}"));
    assert_eq!(
        at(&compared(&expected, &actual, &[]), &["verdict"]),
        &text("red")
    );
    let masked = compared(&expected, &actual, &["regulatory_feature_rows"]);
    assert_eq!(at(&masked, &["verdict"]), &text("green"));
    assert_eq!(at(&masked, &["warnings", "shared"]), &Value::Int(1));
    std::fs::write(&actual, skipped).expect("the file is written");
    let missing = compared(&expected, &actual, &["regulatory_feature_rows"]);
    assert_eq!(at(&missing, &["verdict"]), &text("red"));
    assert_eq!(
        at(&missing, &["reasons"]),
        &json(r#"["warnings: the regulatory notice is missing from the actual file"]"#)
    );
}

#[test]
fn the_fork_class_removes_the_relayed_prefix_from_the_references_lines() {
    let scratch = support::scratch("warnings-fork");
    let expected = write(
        &scratch,
        "ref_warnings.txt",
        "WARNING: 11 : WARNING: Transcript-assembly mismatch in rs922296917\n\
         WARNING: Transcript-assembly mismatch in rs1247639677\n\
         WARNING: 12 : substr outside of string at /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm line 1063, <__ANONIO__> line 41.\n",
    );
    let actual = write(
        &scratch,
        "cand_warnings.txt",
        "WARNING: Transcript-assembly mismatch in rs922296917\n\
         WARNING: Transcript-assembly mismatch in rs1247639677\n",
    );
    let report = compared(&expected, &actual, &["fork_warning_pid_prefix"]);
    let masked = at(&report, &["warnings"]);
    assert_eq!(
        at(masked, &["masks_applied"]),
        &json(r#"{"fork_warning_pid_prefix": 2}"#)
    );
    assert_eq!(at(masked, &["shared"]), &Value::Int(2));
    assert_eq!(at(masked, &["actual_only"]), &json("[]"));
    // The interpreter's line is the child's own text, bare once the relay prefix is gone.
    assert_eq!(
        at(masked, &["expected_only"]),
        &json(
            r#"[{"line": null, "description": "substr outside of string at <path> line 1063, <__ANONIO__> line 41."}]"#
        )
    );
    let report = compared(&expected, &actual, &[]);
    let unmasked = at(&report, &["warnings"]);
    assert_eq!(at(unmasked, &["masks_applied"]), &json("{}"));
    assert_eq!(at(unmasked, &["shared"]), &Value::Int(1));
    assert_eq!(
        descriptions(at(unmasked, &["actual_only"])),
        ["Transcript-assembly mismatch in rs922296917"]
    );
    assert!(descriptions(at(unmasked, &["expected_only"]))
        .iter()
        .any(|description| description.starts_with("11 : WARNING: Transcript-assembly mismatch")));
    // A line with no relay prefix is returned as it is, and counted as unchanged.
    assert_eq!(
        without_forked_child_prefix("WARNING: 11 : WARNING: a\nWARNING: b\n"),
        ("WARNING: a\nWARNING: b\n".to_string(), 1)
    );
}

#[test]
fn the_plugin_classes_set_aside_the_references_failure_lines() {
    let scratch = support::scratch("warnings-plugins");
    let expected = write(
        &scratch,
        "ref_warnings.txt",
        "WARNING: 13 : Use of uninitialized value $orig_a_alt in hash element at /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/Sequence.pm line 1193, <__ANONIO__> line 4.\n\
         Use of uninitialized value $alt in numeric gt (>) at /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/Sequence.pm line 1252, <__ANONIO__> line 4.\n\
         WARNING: Plugin 'CADD' went wrong: \n\
         -------------------- EXCEPTION --------------------\n\
         MSG: Missing reference or alternate sequence\n\
         STACK Bio::EnsEMBL::Variation::Utils::Sequence::trim_sequences /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/Sequence.pm:968\n\
         STACK CADD::run /plugins/CADD.pm:251\n\
         Date (localtime)    = Fri Oct 30 20:39:30 2026\n\
         Ensembl API version = 116\n\
         ---------------------------------------------------\n\
         WARNING: Plugin 'DosageSensitivity' went wrong: Can't locate object method \"transcript\"\n\
         WARNING: Plugin 'pLI' went wrong: something else\n",
    );
    let actual = write(&scratch, "cand_warnings.txt", "");
    let report = compared(
        &expected,
        &actual,
        &[
            "cadd_structural_allele_match",
            "dosage_sensitivity_regulatory_row_warnings",
        ],
    );
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    // The two undefined-allele lines belong to the CADD failure and count under its class.
    assert_eq!(
        at(&report, &["warnings", "masks_applied"]),
        &json(
            r#"{"cadd_structural_allele_match": 3, "dosage_sensitivity_regulatory_row_warnings": 1}"#
        )
    );
    assert_eq!(
        at(&report, &["warnings", "expected_only"]),
        &json(r#"[{"line": null, "description": "Plugin 'pLI' went wrong: something else"}]"#)
    );
    let report = compared(&expected, &actual, &[]);
    let unmasked = at(&report, &["warnings"]);
    assert_eq!(at(unmasked, &["masks_applied"]), &json("{}"));
    assert_eq!(at(unmasked, &["expected_entries"]), &Value::Int(12));
    assert!(descriptions(at(unmasked, &["expected_only"]))
        .iter()
        .any(|description| description.starts_with("Plugin 'CADD'")));
    // The DosageSensitivity class alone leaves the CADD failure and its companion lines in place.
    let report = compared(
        &expected,
        &actual,
        &["dosage_sensitivity_regulatory_row_warnings"],
    );
    let dosage_only = at(&report, &["warnings"]);
    assert_eq!(
        at(dosage_only, &["masks_applied"]),
        &json(r#"{"dosage_sensitivity_regulatory_row_warnings": 1}"#)
    );
    assert_eq!(at(dosage_only, &["expected_entries"]), &Value::Int(11));
    let undefined = descriptions(at(dosage_only, &["expected_only"]))
        .iter()
        .filter(|description| description.contains("uninitialized value"))
        .count();
    assert_eq!(undefined, 2);
    // A block without its closing rule ends at the next warning line.
    assert_eq!(
        without_plugin_failure(
            "WARNING: Plugin 'CADD' went wrong: \n-------------------- EXCEPTION --------------------\nMSG: x\nWARNING: Plugin 'pLI' went wrong: kept\n",
            "Plugin 'CADD' went wrong",
            false,
        ),
        ("WARNING: Plugin 'pLI' went wrong: kept\n".to_string(), 1)
    );
}

/// The class of the mismatch warning a JSON run of Ensembl VEP writes with no record name.
const UNNAMED_MISMATCH_CLASS: &str = "json_mismatch_warning_names_no_record";

/// The interpreter's warning of the statement that joins the missing name, as Ensembl VEP
/// 116.3 words it, after `prefix` and with the statement's line `line`.
fn interpreter_line(prefix: &str, line: u32) -> String {
    format!(
        "{prefix}Use of uninitialized value in concatenation (.) or string at \
         /opt/vep/src/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm line {line}, \
         <__ANONIO__> line 4.\n"
    )
}

/// The warnings file of a forked JSON run of Ensembl VEP over two records that each draw the
/// warning on two rows: one child's relayed block, the unnamed line behind the first
/// interpreter warning.
fn forked_json_reference(line: u32) -> String {
    [
        interpreter_line("WARNING: 11 : ", line),
        "WARNING: Transcript-assembly mismatch in \n".to_string(),
        interpreter_line("", line),
        interpreter_line("", line),
        interpreter_line("", line),
    ]
    .concat()
}

const NAMED_LINES: &str = "WARNING: Transcript-assembly mismatch in rs922296917\n\
                           WARNING: Transcript-assembly mismatch in rs1247639677\n";

#[test]
fn the_json_mismatch_class_sets_aside_the_unnamed_line_beside_the_named_lines() {
    let scratch = support::scratch("warnings-json-mismatch");
    let expected = write(&scratch, "ref_warnings.txt", &forked_json_reference(1727));
    let actual = write(&scratch, "cand_warnings.txt", NAMED_LINES);
    let both = ["fork_warning_pid_prefix", UNNAMED_MISMATCH_CLASS];
    let report = compared(&expected, &actual, &both);
    assert_eq!(at(&report, &["verdict"]), &text("green"));
    let masked = at(&report, &["warnings"]);
    // The unnamed line, the four interpreter lines and the two named lines.
    assert_eq!(
        at(masked, &["masks_applied"]),
        &json(r#"{"fork_warning_pid_prefix": 1, "json_mismatch_warning_names_no_record": 7}"#)
    );
    assert_eq!(at(masked, &["expected_entries"]), &Value::Int(0));
    assert_eq!(at(masked, &["actual_entries"]), &Value::Int(0));

    // Without the fork class the relayed line keeps its prefix and stays the reference's
    // alone; the statement is still read through the prefix, so the bare lines go.
    let report = compared(&expected, &actual, &[UNNAMED_MISMATCH_CLASS]);
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    let relayed = at(&report, &["warnings"]);
    assert_eq!(
        at(relayed, &["masks_applied"]),
        &json(r#"{"json_mismatch_warning_names_no_record": 6}"#)
    );
    assert_eq!(
        descriptions(at(relayed, &["expected_only"])),
        [
            "11 : Use of uninitialized value in concatenation (.) or string at <path> line 1727, \
          <__ANONIO__> line 4."
        ]
    );
    assert_eq!(at(relayed, &["actual_only"]), &json("[]"));

    // Without the class the two files hold nothing in common.
    let report = compared(&expected, &actual, &["fork_warning_pid_prefix"]);
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    let unmasked = at(&report, &["warnings"]);
    assert_eq!(at(unmasked, &["shared"]), &Value::Int(0));
    assert_eq!(
        descriptions(at(unmasked, &["expected_only"])),
        [
            "Transcript-assembly mismatch in",
            "Use of uninitialized value in concatenation (.) or string at <path> line 1727, \
             <__ANONIO__> line 4."
        ]
    );
    assert_eq!(
        descriptions(at(unmasked, &["actual_only"])),
        [
            "Transcript-assembly mismatch in rs1247639677",
            "Transcript-assembly mismatch in rs922296917"
        ]
    );

    // A run without a fork writes the unnamed line alone, and the statement stands at
    // another line in another release.
    std::fs::write(&expected, "WARNING: Transcript-assembly mismatch in \n").expect("written");
    let report = compared(&expected, &actual, &[UNNAMED_MISMATCH_CLASS]);
    assert_eq!(at(&report, &["verdict"]), &text("green"));
    assert_eq!(
        at(&report, &["warnings", "masks_applied"]),
        &json(r#"{"json_mismatch_warning_names_no_record": 3}"#)
    );
    std::fs::write(&expected, forked_json_reference(1718)).expect("written");
    let report = compared(&expected, &actual, &both);
    assert_eq!(at(&report, &["verdict"]), &text("green"));
}

#[test]
fn the_json_mismatch_class_applies_to_its_own_shape_alone() {
    let scratch = support::scratch("warnings-json-mismatch-shape");
    let expected = write(&scratch, "ref_warnings.txt", "");
    let actual = write(&scratch, "cand_warnings.txt", "");
    let both = ["fork_warning_pid_prefix", UNNAMED_MISMATCH_CLASS];
    let nothing_applied = |reference: &str, candidate: &str| {
        std::fs::write(&expected, reference).expect("written");
        std::fs::write(&actual, candidate).expect("written");
        let report = compared(&expected, &actual, &[UNNAMED_MISMATCH_CLASS]);
        assert_eq!(
            at(&report, &["warnings", "masks_applied"]),
            &json("{}"),
            "{reference}"
        );
        at(&report, &["verdict"]).clone()
    };
    let unnamed = "WARNING: Transcript-assembly mismatch in \n";
    // A candidate that names no record keeps the reference's line a difference.
    assert_eq!(nothing_applied(unnamed, ""), text("red"));
    assert_eq!(nothing_applied(unnamed, "WARNING: other\n"), text("red"));
    // A candidate that writes the unnamed line itself agrees with nothing set aside.
    assert_eq!(nothing_applied(unnamed, unnamed), text("green"));
    // A reference without the unnamed line is no JSON run's: its named lines are compared.
    assert_eq!(nothing_applied(NAMED_LINES, NAMED_LINES), text("green"));
    assert_eq!(nothing_applied("", NAMED_LINES), text("red"));
    // A reference that holds a named line beside the unnamed one is not the class's shape.
    let mixed = format!("{unnamed}WARNING: Transcript-assembly mismatch in rs922296917\n");
    assert_eq!(nothing_applied(&mixed, NAMED_LINES), text("red"));

    // Another statement's interpreter warning, in another file or at another line of the
    // same file, and every other warning stay the reference's own.
    let other_line = interpreter_line("", 1800);
    let other_file = "substr outside of string at /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm line 1063, <__ANONIO__> line 41.\n";
    let reference = [
        forked_json_reference(1727).as_str(),
        other_line.as_str(),
        other_file,
        "WARNING: kept\n",
    ]
    .concat();
    std::fs::write(&expected, reference).expect("written");
    std::fs::write(&actual, format!("{NAMED_LINES}WARNING: kept\n")).expect("written");
    let report = compared(&expected, &actual, &both);
    let masked = at(&report, &["warnings"]);
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    assert_eq!(
        at(masked, &["masks_applied"]),
        &json(r#"{"fork_warning_pid_prefix": 1, "json_mismatch_warning_names_no_record": 7}"#)
    );
    assert_eq!(at(masked, &["shared"]), &Value::Int(1));
    assert_eq!(at(masked, &["actual_only"]), &json("[]"));
    assert_eq!(
        descriptions(at(masked, &["expected_only"])),
        [
            "Use of uninitialized value in concatenation (.) or string at <path> line 1800, \
             <__ANONIO__> line 4.",
            "substr outside of string at <path> line 1063, <__ANONIO__> line 41."
        ]
    );

    // An unnamed line the candidate writes beside its named ones stays the candidate's alone.
    std::fs::write(&expected, unnamed).expect("written");
    std::fs::write(&actual, format!("{NAMED_LINES}{unnamed}")).expect("written");
    let report = compared(&expected, &actual, &[UNNAMED_MISMATCH_CLASS]);
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    assert_eq!(
        descriptions(at(&report, &["warnings", "actual_only"])),
        ["Transcript-assembly mismatch in"]
    );
}

/// The six committed flag sets whose JSON reference holds the unnamed line, each with the
/// lines the class sets aside in the reference's relayed block and in the named lines.
const COMMITTED_JSON_REFERENCES: [(&str, i64); 6] = [
    ("GRCh38-refseq/expected", 7),
    ("GRCh38-refseq/expected_all-refseq", 7),
    ("GRCh38-refseq/expected_exclude-predicted", 7),
    ("GRCh38-refseq/expected_gene-symbol-id", 7),
    ("GRCh38-refseq/expected_hgvs", 7),
    ("GRCh37-refseq/expected", 4),
];

#[test]
fn the_committed_json_references_read_green_against_the_named_lines() {
    let scratch = support::scratch("warnings-json-mismatch-committed");
    let both = ["fork_warning_pid_prefix", UNNAMED_MISMATCH_CLASS];
    for (directory, set_aside) in COMMITTED_JSON_REFERENCES {
        let reference = support::golden().join(directory);
        let expected = reference.join("json.warnings.log");
        // The named lines are the ones the same directory's default-format reference holds
        // behind its relay prefix.
        let named = support::text_of(&reference.join("default.warnings.log"));
        let (named, _) = without_forked_child_prefix(&named);
        let actual = write(&scratch, "cand_warnings.txt", &named);
        let report = compared(&expected, &actual, &both);
        assert_eq!(at(&report, &["verdict"]), &text("green"), "{directory}");
        assert_eq!(
            at(&report, &["warnings", "masks_applied"]),
            &Value::Dict(
                [
                    ("fork_warning_pid_prefix".to_string(), Value::Int(1)),
                    (UNNAMED_MISMATCH_CLASS.to_string(), Value::Int(set_aside)),
                ]
                .into_iter()
                .collect()
            ),
            "{directory}"
        );
        let report = compared(&expected, &actual, &["fork_warning_pid_prefix"]);
        assert_eq!(at(&report, &["verdict"]), &text("red"), "{directory}");
    }
}
