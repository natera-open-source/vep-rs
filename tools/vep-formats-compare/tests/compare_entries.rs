// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison of two entry lists on the committed reference outputs and on planted rows.

mod compare_common;

use std::collections::BTreeSet;

use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::compare_files::compare_files_through;
use vep_formats_compare::consts::RECORD_LEVEL_COLUMNS;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::floors::json_floor_fields;
use vep_formats_compare::input::input_records;
use vep_formats_compare::input::read_text;
use vep_formats_compare::input::text_lines;
use vep_formats_compare::manifest::documented_form;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::json::JSON_IDENTITY_KEY;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::regulatory::only_regulatory_rows;
use vep_formats_compare::rules::Side;
use vep_formats_compare::shared::ABSENT;

fn default_entries(text: &str) -> Vec<Entry> {
    parse_default(text).expect("the text parses")
}

fn json_entries(text: &str, inputs: &[String]) -> Vec<Entry> {
    parse_json(text, inputs).expect("the text parses")
}

fn compare_with_floors(
    expected: &[Entry],
    actual: &[Entry],
    documented: &mut Documented,
    floors: &[(&str, i64)],
) -> Value {
    let floors: IndexMap<String, i64> = floors
        .iter()
        .map(|(name, minimum)| (name.to_string(), *minimum))
        .collect();
    let side = Side::default();
    let options = CompareOptions {
        floors: Some(&floors),
        ..CompareOptions::new(&side)
    };
    compare_entries(expected, actual, documented, &options).expect("the entries compare")
}

fn compare_without_order(
    expected: &[Entry],
    actual: &[Entry],
    documented: &mut Documented,
) -> Value {
    let side = Side::default();
    let options = CompareOptions {
        check_order: false,
        ..CompareOptions::new(&side)
    };
    compare_entries(expected, actual, documented, &options).expect("the entries compare")
}

#[test]
fn the_engine_identity_key_of_a_json_object_is_set_aside() {
    let text = text_of("json");
    let inputs = inputs();
    let reference = json_entries(&text, &inputs);
    let with_identity = |line: &str| -> String {
        let mut members = json(line).into_dict().expect("a line is an object");
        members.insert(
            JSON_IDENTITY_KEY.to_string(),
            Value::Str("0.5.0".to_string()),
        );
        pyjson::dumps(&Value::Dict(members), &pyjson::DumpOptions::default())
    };
    let mut json_lines: Vec<String> = text_lines(&text)
        .iter()
        .filter(|line| !line.is_empty())
        .map(|line| with_identity(line))
        .collect();
    json_lines.push(with_identity(
        r#"{"input": "21\t1\tskipped\tA\tG\t.\t.\t."}"#,
    ));
    let entries = json_entries(&(json_lines.join("\n") + "\n"), &inputs);
    let (annotated, skipped) = entries.split_at(reference.len());
    assert_eq!(annotated.len(), reference.len());
    for (ours, theirs) in annotated.iter().zip(&reference) {
        assert_eq!((&ours.key, &ours.fields), (&theirs.key, &theirs.fields));
    }
    assert!(entries
        .iter()
        .all(|e| !e.fields.contains_key(JSON_IDENTITY_KEY)));
    let keys: Vec<&str> = skipped.iter().map(|e| e.key.as_str()).collect();
    assert_eq!(keys, [format!("obj{RECORDS}|record")]);
    assert!(skipped[0].fields.is_empty());
    let report = compare(&reference, annotated, &mut no_documented());
    assert_eq!(verdict(&report), "green");
}

#[test]
fn the_verdict_is_red_on_any_difference() {
    let text = text_of("default");
    let expected = default_entries(&text);
    let rows = data_lines(&text);
    let same = compare(&expected, &default_entries(&text), &mut no_documented());
    assert_eq!(verdict(&same), "green");
    let mut mutated = rows.clone();
    mutated[0] = mutated[0].replace("IMPACT=MODIFIER", "IMPACT=LOW");
    let report = compare(
        &expected,
        &default_entries(&rebuild(&text, &mutated)),
        &mut no_documented(),
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(
        at(&report, &["one_sided"]),
        &json(r#"{"expected_only": 0, "actual_only": 0}"#)
    );
    assert_eq!(number(&report, &["columns", "IMPACT", "differing"]), 1);
    assert_eq!(number(&report, &["columns", "Extra", "differing"]), 1);
    assert_eq!(
        at(&report, &["columns", "IMPACT", "examples"]),
        &json(r#"{"intergenic_variant": ["21:11000001|N.|-|-: expected 'MODIFIER', got 'LOW'"]}"#)
    );
    assert!(reasons(&report).contains(&"field IMPACT: 1 differing entries"));
}

#[test]
fn a_row_one_side_alone_writes_is_counted_on_its_side() {
    let text = text_of("default");
    let expected = default_entries(&text);
    let rows = data_lines(&text);
    let dropped = compare(
        &expected,
        &default_entries(&rebuild(&text, &rows[1..])),
        &mut no_documented(),
    );
    assert_eq!(
        at(&dropped, &["one_sided"]),
        &json(r#"{"expected_only": 1, "actual_only": 0}"#)
    );
    assert_eq!(
        at(&dropped, &["missing"]),
        &json(r#"{"intergenic_variant": {"count": 1, "examples": ["21:11000001|N.|-|-"]}}"#)
    );
    let extra_row = rows[0].replace("21:11000001", "21:11000002");
    let mut more = rows.clone();
    more.push(extra_row);
    let added = compare(
        &expected,
        &default_entries(&rebuild(&text, &more)),
        &mut no_documented(),
    );
    assert_eq!(
        at(&added, &["one_sided"]),
        &json(r#"{"expected_only": 0, "actual_only": 1}"#)
    );
    assert_eq!(
        at(&added, &["unexpected", "intergenic_variant", "examples"]),
        &json(r#"["21:11000002|N.|-|-"]"#)
    );
    assert_eq!(verdict(&added), "red");
    assert_eq!(number(&added, &["counts", "expected_entries"]), ROWS);
    assert_eq!(number(&added, &["counts", "actual_entries"]), ROWS + 1);
}

/// The manifest documents the two mate-side rows of record 1 under
/// breakend_mate_local_read: ENST00000334068 as upstream_gene_variant where the reference
/// prints intergenic_variant, ENST00000390690 as coding_sequence_variant,feature_truncation
/// where it prints feature_truncation.
#[test]
fn a_mask_applies_the_keys_the_manifest_documents_under_it() {
    let text = text_of("default");
    let manifest = manifest();
    let expected = default_entries(&text);
    let rows = data_lines(&text);
    let mut mutated = rows.clone();
    let row_of = |feature: &str| {
        rows.iter()
            .position(|r| r.starts_with("bnd_intergenic_with_genic_mate\t") && r.contains(feature))
            .expect("the reference writes the row")
    };
    let target = row_of("ENST00000334068");
    mutated[target] =
        mutated[target].replace("\tintergenic_variant\t", "\tupstream_gene_variant\t");
    let other = row_of("ENST00000390690");
    mutated[other] = mutated[other].replace(
        "\tfeature_truncation\t",
        "\tfeature_truncation,coding_sequence_variant\t",
    );
    let actual = default_entries(&rebuild(&text, &mutated));
    let klass = "breakend_mate_local_read";

    let unmasked = compare(&expected, &actual, &mut from_manifest(&manifest, &[]));
    assert_eq!(verdict(&unmasked), "red");
    assert_eq!(
        number(&unmasked, &["columns", "Consequence", "differing"]),
        2
    );

    let mut masked = from_run(&manifest, "default", &[klass]);
    assert_eq!(counts(&masked.applied), count_map(&[(klass, 2)]));
    assert_eq!(
        counts(&masked.available),
        count_map(&[
            (klass, 2),
            ("reference_skipped_record", 4),
            ("stop_codon_overlap_insertion_range", 1)
        ])
    );
    let report = compare(&expected, &actual, &mut masked);
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        number(&report, &["counts", "documented_divergences_skipped"]),
        2
    );
    assert_eq!(number(&report, &["counts", "compared"]), ROWS - 2);

    mutated[target] =
        mutated[target].replace("\tupstream_gene_variant\t", "\tdownstream_gene_variant\t");
    let stale = compare(
        &expected,
        &default_entries(&rebuild(&text, &mutated)),
        &mut masked,
    );
    assert_eq!(verdict(&stale), "red");
    assert_eq!(texts(at(&stale, &["stale_divergences"])).len(), 1);
    assert!(reasons(&stale).contains(&"documented divergences not observed: 1"));

    // A flag set's documented keys apply alike once the pair names its run; unnamed, the
    // lists of the two runs are not read as one, and no class named applies nothing.
    let own: Vec<Value> = at(&manifest, &["divergences"])
        .as_list()
        .expect("a list")
        .iter()
        .filter(|d| d.get("expected_divergence").and_then(Value::as_str) == Some(klass))
        .cloned()
        .collect();
    let mut flag_set = json(r#"{"name": "set", "flags": ["--pick"]}"#);
    flag_set = with_member(&flag_set, "divergences", Value::List(own));
    let flag_sets = Value::List(vec![flag_set]);
    let nested = with_member(&manifest, "flag_sets", flag_sets.clone());
    let mut from_set = from_run(&nested, "set", &[klass]);
    assert_eq!(counts(&from_set.applied), count_map(&[(klass, 2)]));
    assert_eq!(
        verdict(&compare(&expected, &actual, &mut from_set)),
        "green"
    );
    let only_set = with_member(
        &with_member(&manifest, "divergences", Value::List(Vec::new())),
        "flag_sets",
        flag_sets,
    );
    assert_eq!(
        counts(&from_run(&only_set, "set", &[klass]).applied),
        count_map(&[(klass, 2)])
    );
    let report = compare(
        &expected,
        &actual,
        &mut from_run(&only_set, "set", &[klass]),
    );
    assert_eq!(verdict(&report), "green");
    assert!(counts(&from_run(&only_set, "default", &[klass]).applied).is_empty());
    let unnamed = Documented::from_manifest(Some(&only_set), &masks(&[klass]), None, None, None)
        .expect_err("two runs document the class");
    assert_eq!(
        unnamed.to_string(),
        "--run is not given, and the manifest lists 2 runs (default, set) that document \
         breakend_mate_local_read, which --masks names: name the pair's run"
    );
    assert!(counts(&from_manifest(&only_set, &[]).applied).is_empty());

    // The built-in class motif_score_reverse_strand_mnv: a differing MOTIF_SCORE_CHANGE on a
    // reverse-strand motif entry with a two-base allele is accepted under the class and no
    // manifest, in the Extra text too; a forward-strand or one-base allele stays a difference.
    let motif_class = "motif_score_reverse_strand_mnv";
    let motif = |strand: i64, score: &str| {
        format!("IMPACT=MODIFIER;STRAND={strand};MOTIF_NAME=ENSPFM0484;MOTIF_POS=9;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE={score}")
    };
    let row = |allele: &str, extra: &str| {
        format!("rev\t21:5090730-5090731\t{allele}\t-\tENSM00000574815\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\t{extra}")
    };
    let reference = [
        DEFAULT_HEADER.to_string(),
        row("AG", &motif(-1, "-0.092")),
        row("A", &motif(-1, "-0.070")).replace("5090730-5090731", "5090730"),
        row("AG", &motif(1, "-0.016")).replace("ENSM00000574815", "ENSM00000574814"),
    ]
    .join("\n")
        + "\n";
    let candidate = reference.replace("-0.092", "-0.136");
    let (report, built_in) = compare_default(&reference, &candidate, &[motif_class]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&built_in.applied), count_map(&[(motif_class, 2)]));
    assert_eq!(
        number(
            &report,
            &[
                "columns",
                "MOTIF_SCORE_CHANGE",
                "motif_score_reverse_strand_mnv"
            ]
        ),
        1
    );
    let red = |reference: &str, candidate: &str, classes: &[&str]| {
        verdict(&compare_default(reference, candidate, classes).0) == "red"
    };
    assert!(red(&reference, &candidate, &[]));
    for other in [
        candidate.replace("-0.070", "-0.071"),
        candidate.replace("-0.016", "-0.017"),
    ] {
        assert!(red(&reference, &other, &[motif_class]));
    }
    // A score the reference omits (the allele ran past the motif's end) and the candidate writes.
    let omitted = reference.replace(";MOTIF_SCORE_CHANGE=-0.092", "");
    assert!(!red(&omitted, &candidate, &[motif_class]));
    assert!(red(&omitted, &candidate, &[]));
    // The other way round, a score the candidate leaves out, is a difference in every format.
    let (unscored, _) = compare_default(&reference, &omitted, &[motif_class]);
    assert_eq!(verdict(&unscored), "red");
    assert_eq!(
        number(&unscored, &["columns", "MOTIF_SCORE_CHANGE", "differing"]),
        1
    );
    // Only a well-formed score is accepted: a difference of two similarity scores in [-1, 1] to
    // three decimals; the accepted entries are reported as the class's examples.
    for malformed in ["9.999", "abc", "-0.09200000", "1.5"] {
        let other = candidate.replace("-0.136", malformed);
        assert!(red(&reference, &other, &[motif_class]), "{malformed}");
        assert!(red(&omitted, &other, &[motif_class]), "{malformed}");
    }
    let (flipped, _) = compare_default(
        &reference,
        &candidate.replace("-0.136", "0.092"),
        &[motif_class],
    );
    assert_eq!(verdict(&flipped), "green");
    assert_eq!(
        at(
            &flipped,
            &[
                "columns",
                "MOTIF_SCORE_CHANGE",
                "motif_score_reverse_strand_mnv_examples"
            ]
        ),
        &json(
            r#"["21:5090730-5090731|AG|ENSM00000574815|MotifFeature: expected '-0.092', got '0.092'"]"#
        )
    );
    // The tab format carries the same columns without the Extra text.
    let tab_header = "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tSTRAND\tMOTIF_SCORE_CHANGE";
    let tab_row = |strand: i64, score: &str| {
        format!("rev\t21:5090730-5090731\tAG\tENSM00000574815\tMotifFeature\tTF_binding_site_variant\t{strand}\t{score}")
    };
    let tab_reference = [
        tab_header.to_string(),
        tab_row(-1, "-0.092"),
        tab_row(1, "-0.016").replace("ENSM00000574815", "ENSM00000574814"),
    ]
    .join("\n")
        + "\n";
    let tab_verdict = |candidate: &str| -> String {
        let expected = parse_tab(&tab_reference).expect("the reference parses");
        let actual = parse_tab(candidate).expect("the candidate parses");
        let report = compare(&expected, &actual, &mut documented_with(&[motif_class]));
        verdict(&report).to_string()
    };
    for (value, colour) in [("-0.136", "green"), ("abc", "red")] {
        assert_eq!(tab_verdict(&tab_reference.replace("-0.092", value)), colour);
    }
    assert_eq!(
        tab_verdict(&tab_reference.replace("-0.016", "-0.017")),
        "red"
    );

    // A manifest's flag sets document their own runs: with the compared input named, a set
    // over another input is left out (its record ordinals count that input's records), and
    // the lists of several runs that document a named class are not read as one.
    let skipped = r#"{"location": "21:100", "allele": "A", "feature": "-", "record_indices": [0], "vep_rs_consequence_sets": ["intergenic_variant"], "expected_divergence": "reference_skipped_record"}"#;
    let field_item = |value: &str| {
        format!(
            r#"{{"location": "21:200", "allele": "C", "feature": "ENST1", "field": "SIFT", "vep_rs": "{value}", "expected_divergence": "sift_rounding"}}"#
        )
    };
    let scoped = |more: &str| {
        json(&format!(
            r#"{{"field_divergences": [{}], "flag_sets": [
                {{"name": "same", "input": "variants.vcf", "field_divergences": [{}]}},
                {{"name": "other", "input": "structural.vcf", "divergences": [], "vep_rs_only_tuples": [{skipped}]}}{more}]}}"#,
            field_item("tolerated(0.1)"),
            field_item("tolerated(0.1)"),
        ))
    };
    let classes = masks(&["reference_skipped_record", "sift_rounding"]);
    let read = |manifest: &Value, input_name: Option<&str>| {
        Documented::from_manifest(Some(manifest), &classes, input_name, None, None)
    };
    let refusal = |manifest: &Value, input_name: Option<&str>| {
        read(manifest, input_name)
            .expect_err("several runs document a named class")
            .to_string()
    };
    let sift_key = (
        "21:200".to_string(),
        "C".to_string(),
        "ENST1".to_string(),
        String::new(),
    );
    let sift = |value: &str| fields(&[("SIFT", value)]);
    let first = scoped("");
    let all_records = |records: &[i64]| -> BTreeSet<i64> { records.iter().copied().collect() };
    assert_eq!(
        refusal(&first, None),
        "--run is not given, and the manifest lists 3 runs (default, same, other) that document \
         reference_skipped_record, sift_rounding, which --masks names: name the pair's run"
    );
    assert_eq!(
        refusal(&first, Some("variants.vcf")),
        "--run is not given, and the manifest lists 2 runs over variants.vcf (default, same) \
         that document sift_rounding, which --masks names: name the pair's run"
    );
    let named = Documented::from_manifest(
        Some(&first),
        &classes,
        Some("variants.vcf"),
        Some("same"),
        None,
    )
    .expect("it reads");
    assert!(named.skipped_records.is_empty());
    assert_eq!(
        named.fields_by_location.get(&sift_key),
        Some(&sift("tolerated(0.1)"))
    );
    let other = read(&first, Some("structural.vcf")).expect("it reads");
    assert_eq!(other.skipped_records, all_records(&[0]));
    assert!(other.fields_by_location.is_empty());
    let gz = r#", {"name": "gz", "input": "variants.vcf.gz"}"#;
    let second = scoped(gz);
    for unnamed_or_bare in ["variants.vcf.gz", "renamed.vcf"] {
        let inherited = read(&second, Some(unnamed_or_bare)).expect("it reads");
        assert_eq!(
            inherited.fields_by_location.get(&sift_key),
            Some(&sift("tolerated(0.1)"))
        );
        assert!(inherited.skipped_records.is_empty());
    }
    // A flag set without `divergences` of its own inherits the manifest's lists, and a field
    // it documents differently on the same key is refused.
    let third = scoped(&format!(
        r#", {{"name": "gz", "input": "variants.vcf.gz", "field_divergences": [{}]}}"#,
        field_item("tolerated(0.2)")
    ));
    let refused = read(&third, Some("variants.vcf.gz")).expect_err("two values of one field");
    assert!(refused
        .to_string()
        .contains("both 'tolerated(0.1)' and 'tolerated(0.2)'"));
    // Two runs that document one record ordinal differently are not read as one list, over
    // two inputs with none named or over the one named input.
    let by_ordinal = |value: &str| {
        format!(
            r#"{{"allele": "C", "feature": "ENST1", "field": "SIFT", "vep_rs": "{value}", "expected_divergence": "sift_rounding", "record_indices": [4]}}"#
        )
    };
    let two_inputs = |input: &str| {
        json(&format!(
            r#"{{"field_divergences": [{}], "flag_sets": [{{"name": "other", "input": "{input}", "divergences": [], "field_divergences": [{}]}}]}}"#,
            by_ordinal("tolerated(0.1)"),
            by_ordinal("tolerated(0.3)"),
        ))
    };
    assert_eq!(
        refusal(&two_inputs("structural.vcf"), None),
        "--run is not given, and the manifest lists 2 runs (default, other) that document \
         sift_rounding, which --masks names: name the pair's run"
    );
    assert_eq!(
        refusal(&two_inputs("variants.vcf"), Some("variants.vcf")),
        "--run is not given, and the manifest lists 2 runs over variants.vcf (default, other) \
         that document sift_rounding, which --masks names: name the pair's run"
    );

    // A documented value compares on the default format's spelling: the VCF's `&` is its comma.
    assert_eq!(
        documented_form("A549_(m&_58_y):ACTIVE&H1_(m):INACTIVE"),
        "A549_(m,_58_y):ACTIVE,H1_(m):INACTIVE"
    );
    assert_eq!(
        documented_form("A549_(m,_58_y):ACTIVE"),
        "A549_(m,_58_y):ACTIVE"
    );
    assert_eq!(documented_form(""), "-");
    assert_eq!(documented_form(ABSENT), "-");
    assert_eq!(documented_form("x%3Dy"), "x=y");
}

#[test]
fn the_regulatory_rows_alone_keep_their_records_checks_and_drop_other_records() {
    let inputs = vec![
        "1\t100\tv0\tA\tT\t.\t.\t.".to_string(),
        "1\t101\tv1\tA\tT\t.\t.\t.".to_string(),
    ];
    let motif = r#"{"variant_allele": "T", "motif_feature_id": "ENSM1", "consequence_terms": ["TF_binding_site_variant"]}"#;
    let object = |id: &str, severe: &str| {
        dumped(&format!(
            r#"{{"input": "1\t100\tv0\tA\tT\t.\t.\t.", "id": "{id}", "most_severe_consequence": "{severe}", "motif_feature_consequences": [{motif}]}}"#
        ))
    };
    let unrelated = dumped(
        r#"{"input": "1\t101\tv1\tA\tT\t.\t.\t.", "id": "v1", "most_severe_consequence": "intergenic_variant"}"#,
    );
    let regulatory = |text: &str| -> Vec<Entry> {
        only_regulatory_rows(&mut json_entries(text, &inputs)).expect("the rows are ranked")
    };
    let expected = regulatory(&format!(
        "{}\n{unrelated}",
        object("v0", "TF_binding_site_variant")
    ));
    let keys: Vec<&str> = expected.iter().map(|e| e.key.as_str()).collect();
    assert_eq!(keys, ["rec0|record", "rec0|T|ENSM1|motif_feature"]);
    let actual = regulatory(&object("wrong", "stop_gained"));
    let report = compare(&expected, &actual, &mut no_documented());
    assert_eq!(number(&report, &["columns", "id", "differing"]), 1);
    assert_eq!(
        number(
            &report,
            &["columns", "most_severe_consequence", "differing"]
        ),
        1
    );
    // A record without regulatory rows is not retained merely to match the other side.
    let actual = regulatory(&dumped(
        r#"{"input": "1\t100\tv0\tA\tT\t.\t.\t.", "id": "v0"}"#,
    ));
    assert!(actual.is_empty());
    let report = compare(&expected, &actual, &mut no_documented());
    assert_eq!(number(&report, &["one_sided", "expected_only"]), 2);
}

/// The committed manifest with one field divergence on the first row of the reference.
fn manifest_documenting(field: &str, vep: &str, vep_rs: &str, klass: &str) -> Value {
    let item = json(&format!(
        r#"[{{"location": "21:11000001", "allele": "N.", "feature": "-", "feature_type": "-",
            "record_indices": [0], "field": "{field}", "vep": "{vep}", "vep_rs": "{vep_rs}",
            "expected_divergence": "{klass}"}}]"#
    ));
    with_member(&manifest(), "field_divergences", item)
}

#[test]
fn a_field_divergence_is_held_to_the_documented_value() {
    let text = text_of("default");
    let documented = manifest_documenting("IMPACT", "MODIFIER", "LOW", "impact_example");
    let expected = default_entries(&text);
    let mut rows = data_lines(&text);
    rows[0] = rows[0].replace("IMPACT=MODIFIER", "IMPACT=LOW");
    let actual = default_entries(&rebuild(&text, &rows));
    let mut masked = from_run(&documented, "default", &["impact_example"]);
    assert_eq!(verdict(&compare(&expected, &actual, &mut masked)), "green");
    rows[0] = rows[0].replace("IMPACT=LOW", "IMPACT=HIGH");
    let stale = compare(
        &expected,
        &default_entries(&rebuild(&text, &rows)),
        &mut masked,
    );
    assert_eq!(verdict(&stale), "red");
    assert!(texts(at(&stale, &["stale_divergences"]))[0]
        .starts_with("21:11000001|N.|-|-: field IMPACT documented as 'LOW'"));
}

/// The class of the item `shared_key_manifest` documents.
const SHARED_KEY_CLASS: &str = "allele_number_repeated_or_reference_equal_alt";

/// A manifest of one `field_divergences` item on a key and consequence set three records
/// share: the allele number of the allele `T` at one position, 2 in the reference and 3 in
/// the candidate on the record that repeats an alternate. `rows_alike` is the item's count
/// of the rows on which both engines print one value, where it states one.
fn shared_key_manifest(rows_alike: Option<i64>) -> Value {
    shared_key_manifest_at("21:30289183", "[0, 1, 2]", rows_alike)
}

/// `shared_key_manifest` with its item at `location` and on the records of `records`, a
/// JSON list of their ordinals.
fn shared_key_manifest_at(location: &str, records: &str, rows_alike: Option<i64>) -> Value {
    let counted = rows_alike.map_or(String::new(), |count| format!(r#", "rows_alike": {count}"#));
    json(&format!(
        r#"{{"field_divergences": [{{"location": "{location}", "allele": "T", "feature": "ENST1",
            "feature_type": "Transcript", "record_indices": {records},
            "consequence_set": "stop_retained_variant", "field": "ALLELE_NUM", "vep": "2",
            "vep_rs": "3", "expected_divergence": "{SHARED_KEY_CLASS}"{counted}}}]}}"#
    ))
}

/// The stale line of an item of `shared_key_manifest_at` that documents `documented` rows
/// alike on the key at `location` where a comparison counts `alike`.
fn rows_alike_line(location: &str, documented: i64, alike: i64) -> String {
    format!("{location}|T|ENST1: field ALLELE_NUM documents {documented} rows alike, {alike} are")
}

/// The three records of the shared key: the first writes `T` alone, the second after one
/// other alternate and the third after an alternate it writes twice.
const SHARED_KEY_RECORDS: [(&str, &str); 3] =
    [("plain", "T"), ("second", "G,T"), ("repeat", "G,G,T")];

/// The rows of the three records on the shared key in `format`, numbered as `numbers`
/// gives them: found by their Location in the default and tab formats and by their records
/// in the VCF and JSON formats.
fn shared_key_rows(format: &str, numbers: [&str; 3]) -> Vec<Entry> {
    let inputs: Vec<String> = SHARED_KEY_RECORDS
        .iter()
        .map(|(name, alternates)| format!("21\t30289183\t{name}\tC\t{alternates}\t.\t.\t."))
        .collect();
    let rows = SHARED_KEY_RECORDS.iter().zip(numbers);
    let entries = match format {
        "default" => {
            let rows: Vec<String> = rows
                .map(|((name, _), number)| {
                    let extra = format!("ALLELE_NUM={number};IMPACT=LOW;STRAND=-1");
                    let terms = "stop_retained_variant";
                    default_row(name, "21:30289183", "T", "ENST1", terms, &extra)
                })
                .collect();
            parse_default(&default_text(&rows))
        }
        "tab" => {
            let mut text = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tALLELE_NUM\tIMPACT\n".to_string();
            for ((name, _), number) in rows {
                text += &format!(
                    "{name}\t21:30289183\tT\tENSG1\tENST1\tTranscript\tstop_retained_variant\t{number}\tLOW\n"
                );
            }
            parse_tab(&text)
        }
        "vcf" => {
            let mut text = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|ALLELE_NUM\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_string();
            for (input, number) in inputs.iter().zip(numbers) {
                let fixed = input.strip_suffix('.').expect("an INFO that is a point");
                text += &format!("{fixed}CSQ=T|stop_retained_variant|Transcript|ENST1|{number}\n");
            }
            parse_vcf(&text, &inputs)
        }
        _ => {
            let mut text = String::new();
            for (input, number) in inputs.iter().zip(numbers) {
                text += &format!(
                    "{{\"input\":\"{}\",\"transcript_consequences\":[{{\"variant_allele\":\"T\",\"transcript_id\":\"ENST1\",\"consequence_terms\":[\"stop_retained_variant\"],\"allele_num\":{number}}}]}}\n",
                    input.replace('\t', "\\t")
                );
            }
            parse_json(&text, &inputs)
        }
    };
    entries.expect("the rows parse")
}

/// A documented field on a key and consequence set several records share: an item that
/// states `rows_alike` holds the rows that differ to its value, and a row on which both
/// engines print one value is compared like any other, in every format. An item that
/// states no `rows_alike` holds every row of its key.
#[test]
fn a_documented_field_on_a_shared_key_holds_the_rows_that_differ() {
    for format in ["default", "tab", "vcf", "json"] {
        let compared = |candidate: [&str; 3], rows_alike: Option<i64>| {
            let manifest = shared_key_manifest(rows_alike);
            let mut documented = from_manifest(&manifest, &[SHARED_KEY_CLASS]);
            compare(
                &shared_key_rows(format, ["1", "2", "2"]),
                &shared_key_rows(format, candidate),
                &mut documented,
            )
        };
        let stale = |report: &Value| -> Vec<String> {
            texts(at(report, &["stale_divergences"]))
                .iter()
                .map(|line| line.to_string())
                .collect()
        };
        let field = if format == "json" {
            "allele_num"
        } else {
            "ALLELE_NUM"
        };
        let shown = |printed: &str| {
            format!("field {field} documented as '3', the candidate now prints '{printed}'")
        };
        // The row that differs prints the item's value, and the two others agree.
        let report = compared(["1", "2", "3"], Some(2));
        assert_eq!(verdict(&report), "green", "{format}");
        assert!(stale(&report).is_empty(), "{format}");
        // The row that differs is held to the item's value.
        let report = compared(["1", "2", "4"], Some(2));
        assert_eq!(verdict(&report), "red", "{format}");
        let lines = stale(&report);
        assert_eq!(lines.len(), 1, "{format} {lines:?}");
        assert!(lines[0].contains(&shown("4")), "{format} {lines:?}");
        // A row of another record that differs is held to it as well, and one row is left
        // alike where the item documents two.
        let report = compared(["2", "2", "3"], Some(2));
        assert_eq!(verdict(&report), "red", "{format}");
        let lines = stale(&report);
        assert_eq!(lines.len(), 2, "{format} {lines:?}");
        assert!(lines[0].contains(&shown("2")), "{format} {lines:?}");
        assert_eq!(lines[1], rows_alike_line("21:30289183", 2, 1), "{format}");
        // An item that states no `rows_alike` holds every row of its key.
        let report = compared(["1", "2", "3"], None);
        assert_eq!(verdict(&report), "red", "{format}");
        let lines = stale(&report);
        assert_eq!(lines.len(), 2, "{format} {lines:?}");
        assert!(lines[0].contains(&shown("1")), "{format} {lines:?}");
        assert!(lines[1].contains(&shown("2")), "{format} {lines:?}");
    }
}

/// The count an item states under `rows_alike` is held in every format: the rows of its
/// key on which both engines print one value are counted, over the rows of the Location in
/// the default and tab formats and over the entries of the item's records in the VCF and
/// JSON formats, and an item whose count is another number is one stale line. An item whose
/// key no compared row reaches with the field is held to no count, and an item that states
/// no `rows_alike` counts nothing.
#[test]
fn a_documented_field_on_a_shared_key_holds_its_count_of_rows_alike() {
    let here = "21:30289183";
    for format in ["default", "tab", "vcf", "json"] {
        let compared = |candidate: [&str; 3], manifest: &Value| {
            let mut documented = from_manifest(manifest, &[SHARED_KEY_CLASS]);
            let report = compare(
                &shared_key_rows(format, ["1", "2", "2"]),
                &shared_key_rows(format, candidate),
                &mut documented,
            );
            let stale: Vec<String> = texts(at(&report, &["stale_divergences"]))
                .iter()
                .map(|line| line.to_string())
                .collect();
            (verdict(&report).to_string(), stale)
        };
        let green = ("green".to_string(), Vec::new());
        let red = |line: String| ("red".to_string(), vec![line]);
        // Two rows agree beside the one that differs, the count the item states.
        let two = shared_key_manifest(Some(2));
        assert_eq!(compared(["1", "2", "3"], &two), green, "{format}");
        // Every row is alike: no row prints the item's value, and three agree.
        assert_eq!(
            compared(["1", "2", "2"], &two),
            red(rows_alike_line(here, 2, 3)),
            "{format}"
        );
        // Fewer rows, and more rows, are alike than the item states.
        for documented in [3, 1, 0] {
            assert_eq!(
                compared(["1", "2", "3"], &shared_key_manifest(Some(documented))),
                red(rows_alike_line(here, documented, 2)),
                "{format} {documented}"
            );
        }
        // An item whose key no row of the comparison reaches is held to no count: the
        // outputs of another run over the same manifest are not stale for it.
        let elsewhere = shared_key_manifest_at("21:30289999", "[7]", Some(2));
        assert_eq!(compared(["1", "2", "2"], &elsewhere), green, "{format}");
        // A count that is no whole number from 0 up reads as 0.
        for written in [r#""2""#, "2.0", "-2", "null"] {
            let odd = json(
                &pyjson::dumps_compact(&two)
                    .replace(r#""rows_alike":2"#, &format!(r#""rows_alike":{written}"#)),
            );
            assert_eq!(
                compared(["1", "2", "3"], &odd),
                red(rows_alike_line(here, 0, 2)),
                "{format} {written}"
            );
        }
        // An item that states no `rows_alike` holds every row to its value and counts none.
        let whole = shared_key_manifest(None);
        assert_eq!(compared(["3", "3", "3"], &whole), green, "{format}");
        // One `Documented` read for several comparisons counts each comparison's rows alone.
        let mut documented = from_manifest(&two, &[SHARED_KEY_CLASS]);
        for _ in 0..2 {
            let report = compare(
                &shared_key_rows(format, ["1", "2", "2"]),
                &shared_key_rows(format, ["1", "2", "3"]),
                &mut documented,
            );
            assert_eq!(verdict(&report), "green", "{format}");
        }
    }
}

/// A documented default-format value `1311092.0000` is the JSON number `1311092`, and
/// `0.2010` is `0.201`; another number, or another spelling of the same number, is stale.
#[test]
fn a_documented_numeric_value_is_read_in_the_json_number_form() {
    let text = text_of("default");
    let klass = "minimal_plugin_sees_the_rejoined_record";
    let documenting = |vep_rs: &str| -> Documented {
        let manifest = manifest_documenting("gnomAD_total_DP_cov", "1311160.0000", vep_rs, klass);
        from_run(&manifest, "default", &[klass])
    };
    let rows = data_lines(&text);
    let with_value = |value: &str| -> Vec<Entry> {
        let mut planted = rows.clone();
        planted[0] = planted[0].replace(
            "IMPACT=MODIFIER",
            &format!("IMPACT=MODIFIER;gnomAD_total_DP_cov={value}"),
        );
        default_entries(&rebuild(&text, &planted))
    };
    let expected = with_value("1311160.0000");
    let mut masked = documenting("1311092.0000");
    for same in ["1311092.0000", "1311092"] {
        let report = compare(&expected, &with_value(same), &mut masked);
        assert_eq!(verdict(&report), "green", "{same}");
    }
    for other in ["1311093", "1311092.0", "1311092.00000"] {
        let stale = compare(&expected, &with_value(other), &mut masked);
        assert_eq!(verdict(&stale), "red", "{other}");
        assert!(texts(at(&stale, &["stale_divergences"]))[0]
            .contains("gnomAD_total_DP_cov documented as '1311092.0000'"));
    }
    // A fraction keeps its significant digits; only a decimal number is shortened, never a
    // dotted name ending in a zero.
    for (recorded, printed, colour) in [
        ("0.2010", "0.201", "green"),
        ("0.2010", "0.21", "red"),
        ("p.Ter10", "p.Ter1", "red"),
        ("p.Ter10", "p.Ter10", "green"),
    ] {
        let report = compare(&expected, &with_value(printed), &mut documenting(recorded));
        assert_eq!(verdict(&report), colour, "{recorded} {printed}");
    }
}

/// Rows of the allele-flags reference carrying AMBIGUITY in the default format, and the record
/// objects carrying ambiguity in its JSON format: one per record with an ambiguity code.
const AMBIGUITY_ROWS: i64 = 115;
const AMBIGUITY_OBJECTS: usize = 13;

/// A floor is met, unmet or empty on the default format; a floor on a record-level column
/// is not applicable in the JSON format: `AMBIGUITY` is a key of
/// `VariationFeature_to_output_hash` (OutputFactory.pm:995-1000), as `SV`, `IND`,
/// `CHECK_REF`, `VARIANT_CLASS` and `MINIMISED` are, and the JSON factory copies every
/// such key once to the record object (JSON.pm:262), never into a consequence entry, so the
/// report lists the floor under `floors_not_applicable` for the text formats of the case to
/// check, and the default format meets the same floor.
#[test]
fn a_floor_is_met_unmet_or_not_applicable() {
    let text = text_of("default");
    let expected = default_entries(&text);
    let actual = default_entries(&text);
    let met = compare_with_floors(
        &expected,
        &actual,
        &mut no_documented(),
        &[("IMPACT", ROWS), ("OverlapBP", 1)],
    );
    assert_eq!(verdict(&met), "green");
    assert_eq!(
        at(&met, &["columns", "IMPACT"]),
        &json(&format!(
            r#"{{"expected_populated": {ROWS}, "actual_populated": {ROWS}, "differing": 0, "pairs_accepted": 0,
                "nearest_accepted": 0, "order_only": 0, "individual_allele_order": 0,
                "motif_score_reverse_strand_mnv": 0, "floor": {ROWS}, "floor_met": true}}"#
        ))
    );
    assert!(texts(at(&met, &["floors_not_applicable"])).is_empty());
    let empty = compare_with_floors(&expected, &actual, &mut no_documented(), &[("HGVSc", 1)]);
    assert_eq!(verdict(&empty), "red");
    assert_eq!(
        reasons(&empty),
        ["floor HGVSc=1: expected populated 0, actual populated 0"]
    );
    let short = compare_with_floors(
        &expected,
        &actual,
        &mut no_documented(),
        &[("IMPACT", ROWS + 1)],
    );
    assert_eq!(verdict(&short), "red");
    assert_eq!(
        at(&short, &["columns", "IMPACT", "floor_met"]),
        &Value::Bool(false)
    );

    let output_flags = golden().join("GRCh37-output-flags");
    let expected_dir = output_flags.join("expected_allele-flags");
    let inputs = input_records(&output_flags.join("variants.vcf")).expect("the input reads");
    for column in ["AMBIGUITY", "IND", "SV", "CHECK_REF"] {
        assert!(RECORD_LEVEL_COLUMNS.contains(&column));
        assert!(json_floor_fields(column).is_none());
    }
    let json_path = expected_dir.join("json.jsonl.gz");
    let json_text = read_text(&json_path).expect("the output reads");
    let objects = text_lines(&json_text)
        .iter()
        .filter(|line| line.contains("\"ambiguity\""))
        .count();
    assert_eq!(objects, AMBIGUITY_OBJECTS);
    let files = |fmt: &str,
                 path: &std::path::Path,
                 inputs: Option<&[String]>,
                 floors: &[(&str, i64)]|
     -> Value {
        let floors: IndexMap<String, i64> = floors
            .iter()
            .map(|(name, minimum)| (name.to_string(), *minimum))
            .collect();
        compare_files_through(
            None,
            fmt,
            path,
            path,
            inputs,
            &mut no_documented(),
            &BTreeSet::new(),
            &floors,
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
    let report = files(
        "json",
        &json_path,
        Some(&inputs),
        &[
            ("AMBIGUITY", AMBIGUITY_ROWS),
            ("VARIANT_CLASS", 1),
            ("IND", 1),
            ("SV", 1),
            ("CHECK_REF", 1),
            ("IMPACT", 1),
        ],
    );
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(
        texts(at(&report, &["floors_not_applicable"])),
        ["AMBIGUITY", "CHECK_REF", "IND", "SV", "VARIANT_CLASS"]
    );
    assert_eq!(
        at(&report, &["columns", "AMBIGUITY", "floor_met"]),
        &Value::Null
    );
    assert_eq!(
        at(&report, &["columns", "IMPACT", "floor_met"]),
        &Value::Bool(true)
    );
    let default_path = expected_dir.join("default.txt.gz");
    let text = files(
        "default",
        &default_path,
        None,
        &[("AMBIGUITY", AMBIGUITY_ROWS)],
    );
    assert_eq!(verdict(&text), "green", "{:?}", reasons(&text));
    assert!(texts(at(&text, &["floors_not_applicable"])).is_empty());
    assert_eq!(
        number(&text, &["columns", "AMBIGUITY", "expected_populated"]),
        AMBIGUITY_ROWS
    );
    assert_eq!(
        at(&text, &["columns", "AMBIGUITY", "floor_met"]),
        &Value::Bool(true)
    );
}

#[test]
fn an_order_difference_is_red_unless_the_order_is_not_read() {
    let text = text_of("default");
    let expected = default_entries(&text);
    let rows = data_lines(&text);
    let mut swapped = rows.clone();
    swapped.swap(0, 1);
    let actual = default_entries(&rebuild(&text, &swapped));
    let ordered = compare(&expected, &actual, &mut no_documented());
    assert_eq!(verdict(&ordered), "red");
    assert_eq!(
        at(&ordered, &["order"]),
        &json(
            r#"{"first_divergence": 0, "expected": "21:11000001|N.|-|-", "actual": "21:11000001|N[21:11500000[|-|-"}"#
        )
    );
    let unordered = compare_without_order(&expected, &actual, &mut no_documented());
    assert_eq!(verdict(&unordered), "green");
    assert_eq!(at(&unordered, &["order"]), &Value::Null);
}

/// start_cooccurrence_swap: Ensembl VEP's start_lost beside start_retained_variant on a
/// sequence allele against the candidate's set without start_lost (the structural allele drops
/// start_retained_variant instead); splice_family_swap: the candidate's splice_region_variant
/// the reference lost, with or without intron_variant. The pair is skipped whole, so IMPACT
/// and HGVSp on it never count, and every other shape stays a difference.
#[test]
fn the_registry_pair_classes_set_aside_the_excluding_rules_shapes_without_a_manifest() {
    let reference = default_text(&[
        default_row(
            "a",
            "21:100",
            "G",
            "ENST1",
            "start_lost,start_retained_variant,5_prime_UTR_variant",
            "IMPACT=HIGH;STRAND=1;HGVSp=ENSP1:p.Met1?",
        ),
        default_row(
            "b",
            "21:200-250",
            "deletion",
            "ENST2",
            "start_lost,start_retained_variant,coding_sequence_variant",
            "IMPACT=HIGH;STRAND=1",
        ),
        default_row(
            "c",
            "21:300",
            "T",
            "ENST3",
            "intron_variant",
            "IMPACT=MODIFIER;STRAND=1",
        ),
        default_row(
            "d",
            "21:400",
            "A",
            "ENST4",
            "missense_variant",
            "IMPACT=MODERATE;STRAND=1",
        ),
    ]);
    let candidate = default_text(&[
        default_row(
            "a",
            "21:100",
            "G",
            "ENST1",
            "start_retained_variant,5_prime_UTR_variant",
            "IMPACT=LOW;STRAND=1;HGVSp=ENSP1:p.Met1Val",
        ),
        default_row(
            "b",
            "21:200-250",
            "deletion",
            "ENST2",
            "start_lost,coding_sequence_variant",
            "IMPACT=HIGH;STRAND=1",
        ),
        default_row(
            "c",
            "21:300",
            "T",
            "ENST3",
            "splice_region_variant,intron_variant",
            "IMPACT=LOW;STRAND=1",
        ),
        default_row(
            "d",
            "21:400",
            "A",
            "ENST4",
            "missense_variant",
            "IMPACT=MODERATE;STRAND=1",
        ),
    ]);
    let both = ["start_cooccurrence_swap", "splice_family_swap"];
    let applied_both = count_map(&[("start_cooccurrence_swap", 2), ("splice_family_swap", 1)]);
    let (report, documented) = compare_default(&reference, &candidate, &both);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), applied_both);
    assert_eq!(number(&report, &["counts", "registry_pairs_set_aside"]), 3);
    assert_eq!(number(&report, &["counts", "compared"]), 1);
    assert!(
        texts(at(&report, &["mask_examples", "start_cooccurrence_swap"]))[0].starts_with(
            "21:100|G|ENST1|Transcript: expected 'start_lost,start_retained_variant,5_prime_UTR_variant'"
        )
    );
    // One class alone leaves the other shape a difference; no class leaves both.
    let (only_start, _) = compare_default(&reference, &candidate, &["start_cooccurrence_swap"]);
    assert_eq!(verdict(&only_start), "red");
    assert_eq!(
        number(&only_start, &["columns", "Consequence", "differing"]),
        1
    );
    let (neither, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(
        number(&neither, &["columns", "Consequence", "differing"]),
        3
    );
    // The reverse shapes are not the registry's: a candidate dropping the wrong start term for its
    // allele kind, or a reference-side loss of any other splice term.
    let wrong_term = candidate.replace(
        "start_retained_variant,5_prime_UTR_variant",
        "start_lost,5_prime_UTR_variant",
    );
    assert_eq!(
        verdict(&compare_default(&reference, &wrong_term, &both).0),
        "red"
    );
    let other_splice = candidate.replace(
        "splice_region_variant,intron_variant",
        "splice_donor_variant,intron_variant",
    );
    assert_eq!(
        verdict(&compare_default(&reference, &other_splice, &both).0),
        "red"
    );
    // A documented key keeps the manifest's reading and is not counted twice.
    let manifest = json(
        r#"{"divergences": [{"location": "21:100", "allele": "G", "feature": "ENST1", "expected_divergence": "start_cooccurrence_swap",
            "vep_rs_consequence_sets": ["5_prime_UTR_variant,start_retained_variant"]}]}"#,
    );
    let mut documented = from_manifest(&manifest, &both);
    let report = compare_texts(&reference, &candidate, &mut documented);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), applied_both);
    assert_eq!(
        number(&report, &["counts", "documented_divergences_skipped"]),
        1
    );
    assert_eq!(number(&report, &["counts", "registry_pairs_set_aside"]), 2);
}

/// The fsTer or extTer count of an HGVSp differs on a STRAND -1 row whose SOURCE names the
/// file that holds its transcript, and nothing else does; a forward-strand row, another
/// part of the notation or a missing value stays a difference.
#[test]
fn a_ter_count_on_a_reverse_strand_row_is_accepted_under_its_class() {
    let klass = "gxf_reverse_strand_utr_sequence";
    let reference = default_text(&[
        default_row(
            "a",
            "21:100",
            "-",
            "ENST1",
            "frameshift_variant",
            "IMPACT=HIGH;STRAND=-1;SOURCE=genes.gtf.gz;HGVSc=ENST1:c.100del;HGVSp=ENSP1:p.Lys34ArgfsTer41",
        ),
        default_row(
            "b",
            "21:200",
            "T",
            "ENST2",
            "stop_lost",
            "IMPACT=HIGH;STRAND=-1;SOURCE=genes.gtf.gz;HGVSp=ENSP2:p.Ter300GlnextTer20",
        ),
        default_row(
            "c",
            "21:300",
            "-",
            "ENST3",
            "frameshift_variant",
            "IMPACT=HIGH;STRAND=1;SOURCE=genes.gtf.gz;HGVSp=ENSP3:p.Lys10ArgfsTer5",
        ),
    ]);
    let candidate = reference
        .replace("fsTer41", "fsTer69")
        .replace("extTer20", "extTer18");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(klass, 2)]));
    let masked = json(r#"{"gxf_reverse_strand_utr_sequence": 2}"#);
    assert_eq!(at(&report, &["columns", "HGVSp", "masked"]), &masked);
    assert_eq!(at(&report, &["columns", "Extra", "masked"]), &masked);
    let colour = |reference: &str, candidate: &str, classes: &[&str]| -> String {
        verdict(&compare_default(reference, candidate, classes).0).to_string()
    };
    assert_eq!(colour(&reference, &candidate, &[]), "red");
    for other in [
        reference.replace("fsTer5", "fsTer7"),
        reference.replace("Lys34Arg", "Lys34Gly"),
        reference.replace(";HGVSp=ENSP1:p.Lys34ArgfsTer41", ""),
    ] {
        assert_eq!(colour(&reference, &other, &[klass]), "red");
    }
    // A count the reference leaves unknown is the same shape.
    let unknown = reference.replace("fsTer41", "fsTer?");
    assert_eq!(colour(&unknown, &candidate, &[klass]), "green");
    // The tab format carries the field as a column.
    let tab_header =
        "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tSTRAND\tSOURCE\tHGVSp";
    let tab = format!(
        "{tab_header}\na\t21:100\t-\tENST1\tTranscript\tframeshift_variant\t-1\tgenes.gtf.gz\tENSP1:p.Lys34ArgfsTer41\n"
    );
    let mut documented = documented_with(&[klass]);
    let report = compare(
        &parse_tab(&tab).expect("the reference parses"),
        &parse_tab(&tab.replace("fsTer41", "fsTer69")).expect("the candidate parses"),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(klass, 1)]));
}

/// Under MAX the reference's text sort reports the most negative covered score where the
/// candidate reports the numeric maximum, so the reference's value lies below the candidate's
/// and both are negative; any other pair of values is a difference.
#[test]
fn a_lower_negative_reference_conservation_is_accepted_under_its_class() {
    let klass = "conservation_max_text_order";
    let reference = default_text(&[default_row(
        "a",
        "21:44293093-44293094",
        "-",
        "ENST1",
        "intron_variant",
        "IMPACT=MODIFIER;STRAND=1;Conservation=-5.200",
    )]);
    let candidate = reference.replace("-5.200", "-2.170");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(klass, 1)]));
    let colour = |reference: &str, candidate: &str, classes: &[&str]| -> String {
        verdict(&compare_default(reference, candidate, classes).0).to_string()
    };
    for value in ["-6.100", "0.000", "1.250", "abc"] {
        let other = reference.replace("-5.200", value);
        assert_eq!(colour(&reference, &other, &[klass]), "red", "{value}");
    }
    assert_eq!(
        colour(
            &reference.replace("-5.200", "2.000"),
            &reference.replace("-5.200", "3.000"),
            &[klass]
        ),
        "red"
    );
    assert_eq!(colour(&reference, &candidate, &[]), "red");
    // The JSON key is conservation and its value a number.
    let inputs = vec!["21\t44293093\trs1\tTAC\tT\t.\t.\t.".to_string()];
    let reference_json = dumped(
        r#"{"input": "21\t44293093\trs1\tTAC\tT\t.\t.\t.", "id": "rs1", "most_severe_consequence": "intron_variant", "transcript_consequences": [
            {"variant_allele": "-", "transcript_id": "ENST1", "consequence_terms": ["intron_variant"], "conservation": -5.2}]}"#,
    ) + "\n";
    let mut documented = documented_with(&[klass]);
    let report = compare(
        &json_entries(&reference_json, &inputs),
        &json_entries(&reference_json.replace("-5.2", "-2.17"), &inputs),
        &mut documented,
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[(klass, 1)]));
}
