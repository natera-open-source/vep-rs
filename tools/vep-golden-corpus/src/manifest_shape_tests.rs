// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The shapes a manifest is held to: each refused shape with its reason, in the words every
//! subcommand refuses it with, and the shapes beside it that pass.

use super::*;
use crate::constants::{CORRECTED_BEHAVIOUR_CLASSES, HARNESS_ONLY_MASK_CLASSES};

fn object(text: &str) -> Object {
    match loads(text) {
        Ok(Value::Object(object)) => object,
        other => panic!("{text} is no JSON object: {other:?}"),
    }
}

fn shape(text: &str) -> Option<String> {
    manifest_shape_error(&object(text))
}

fn classes_text() -> String {
    mask_classes().join(", ")
}

#[test]
fn plugin_files_shape_is_held_at_every_read() {
    let good = r#"{"LoFTEE": ["plugin_data/human_ancestor.fa.gz", "plugin_data/gerp.bw"]}"#;
    assert_eq!(plugin_files_shape_error(&loads(good).unwrap()), None);
    assert_eq!(
        shape(&format!(r#"{{"flags": [], "plugin_files": {good}}}"#)),
        None
    );
    let cases = [
        (
            r#"["plugin_data/gerp.bw"]"#,
            "plugin_files is not a map of plugin names",
        ),
        (
            r#"{"": ["plugin_data/gerp.bw"]}"#,
            "plugin_files is not a map of plugin names",
        ),
        (
            r#"{"LoFTEE": "plugin_data/gerp.bw"}"#,
            "plugin_files for LoFTEE is not a list of paths",
        ),
        (
            r#"{"LoFTEE": ["plugin_data/gerp.bw", ""]}"#,
            "plugin_files for LoFTEE is not a list of paths",
        ),
        (
            r#"{"LoFTEE": ["/data/gerp.bw"]}"#,
            "plugin_files for LoFTEE names paths outside the corpus ['/data/gerp.bw']",
        ),
        (
            r#"{"LoFTEE": ["../gerp.bw"]}"#,
            "plugin_files for LoFTEE names paths outside the corpus ['../gerp.bw']",
        ),
        (
            r#"{"LoFTEE": ["plugin data/gerp.bw"]}"#,
            "plugin_files for LoFTEE names paths outside the corpus ['plugin data/gerp.bw']",
        ),
    ];
    for (value, reason) in cases {
        assert_eq!(
            plugin_files_shape_error(&loads(value).unwrap()).as_deref(),
            Some(reason),
            "{value}"
        );
        assert_eq!(
            shape(&format!(r#"{{"flags": [], "plugin_files": {value}}}"#)).as_deref(),
            Some(reason),
            "{value}"
        );
    }
}

#[test]
fn the_lists_of_names_and_the_modes_are_held_to_their_shapes() {
    let unknown = format!(
        "masks names unknown classes ['alpha', \"it's\", 'zeta'] (the classes are {})",
        classes_text()
    );
    let cases: [(&str, Option<&str>); 24] = [
        (r#"{"covers": []}"#, None),
        (r#"{"covers": null}"#, Some("covers is not a list of names")),
        (
            r#"{"covers": ["symbol", " "]}"#,
            Some("covers is not a list of names"),
        ),
        (r#"{"masks": []}"#, None),
        (r#"{"masks": ["order_only", "order_only"]}"#, None),
        (
            r#"{"masks": ["zeta", "alpha", "zeta", "it's"]}"#,
            Some(&unknown),
        ),
        (r#"{"masks": [" "]}"#, Some("masks is not a list of names")),
        (r#"{"reference_warnings": []}"#, None),
        (
            r#"{"reference_warnings": ["no FASTA", 5]}"#,
            Some("reference_warnings is not a list of texts"),
        ),
        (r#"{"plugin_files": {}}"#, None),
        (r#"{"plugin_files": {"P": []}}"#, None),
        (
            "{\"plugin_files\": {\"P\": [\"a\\tb\", \"ok\", \"/x\", \"c\u{a0}d\"]}}",
            Some("plugin_files for P names paths outside the corpus ['a\\tb', '/x', 'c\\xa0d']"),
        ),
        (
            r#"{"plugin_files": {" ": ["a"]}}"#,
            Some("plugin_files is not a map of plugin names"),
        ),
        (
            r#"{"plugin_files": null}"#,
            Some("plugin_files is not a map of plugin names"),
        ),
        (r#"{"perl_mode": "offline"}"#, None),
        (r#"{"perl_mode": "gff"}"#, None),
        (
            r#"{"perl_mode": "cache"}"#,
            Some("perl_mode is not one of offline, database, gtf, gff"),
        ),
        (
            r#"{"perl_mode": null}"#,
            Some("perl_mode is not one of offline, database, gtf, gff"),
        ),
        (
            r#"{"perl_mode": ["gtf"]}"#,
            Some("perl_mode is not one of offline, database, gtf, gff"),
        ),
        (r#"{"record_class_key": "regulatory_class"}"#, None),
        (
            r#"{"record_class_key": "other"}"#,
            Some("record_class_key is not one of regulatory_class"),
        ),
        (
            r#"{"record_class_key": null}"#,
            Some("record_class_key is not one of regulatory_class"),
        ),
        (r#"{"required_genotype_classes": []}"#, None),
        (
            r#"{"required_genotype_classes": ["het", "one_hom_ref_rest_missing"]}"#,
            None,
        ),
    ];
    for (text, reason) in cases {
        assert_eq!(shape(text).as_deref(), reason, "{text}");
    }
    // Every mask class passes, the harness's own among them.
    let every: Vec<&str> = mask_classes();
    assert_eq!(masks_shape_error(&Value::strings(&every), "masks"), None);
    assert_eq!(
        masks_shape_error(&Value::strings(&CORRECTED_BEHAVIOUR_CLASSES), "masks"),
        None
    );
    assert_eq!(
        masks_shape_error(
            &Value::strings(&HARNESS_ONLY_MASK_CLASSES),
            "flag set c masks"
        ),
        None
    );
}

#[test]
fn the_counts_the_species_and_the_cache_keys_are_held_to_their_shapes() {
    let cases: [(&str, Option<&str>); 37] = [
        (r#"{"keep_samples": 0}"#, None),
        (r#"{"keep_samples": -0}"#, None),
        (r#"{"keep_samples": 2, "samples": ["a", "b"]}"#, None),
        (
            r#"{"keep_samples": 2.0}"#,
            Some("keep_samples is not a count"),
        ),
        (
            r#"{"keep_samples": null}"#,
            Some("keep_samples is not a count"),
        ),
        (
            r#"{"keep_samples": "8"}"#,
            Some("keep_samples is not a count"),
        ),
        (
            r#"{"keep_samples": -1}"#,
            Some("keep_samples is not a count"),
        ),
        (
            r#"{"keep_samples": true}"#,
            Some("keep_samples is not a count"),
        ),
        (
            r#"{"keep_samples": 2}"#,
            Some("samples does not name 2 kept samples"),
        ),
        (
            r#"{"keep_samples": 2, "samples": ["a"]}"#,
            Some("samples does not name 2 kept samples"),
        ),
        (
            r#"{"keep_samples": 1, "samples": [""]}"#,
            Some("samples does not name 1 kept samples"),
        ),
        (
            r#"{"keep_samples": 1, "samples": "a"}"#,
            Some("samples does not name 1 kept samples"),
        ),
        (
            r#"{"keep_samples": 99999999999999999999, "samples": ["a"]}"#,
            Some("samples does not name 99999999999999999999 kept samples"),
        ),
        (
            r#"{"required_genotype_classes": "het"}"#,
            Some("required_genotype_classes is not a list of genotype classes"),
        ),
        (
            r#"{"required_genotype_classes": ["het", "ghost"]}"#,
            Some("required_genotype_classes is not a list of genotype classes"),
        ),
        (
            r#"{"required_genotype_classes": [["het"]]}"#,
            Some("required_genotype_classes is not a list of genotype classes"),
        ),
        (r#"{"species": "homo_sapiens"}"#, None),
        (r#"{"species": "mus_musculus_129s1svimj"}"#, None),
        (
            r#"{"species": "Homo sapiens"}"#,
            Some(r#"species is "Homo sapiens", not an Ensembl species name"#),
        ),
        (
            r#"{"species": 5}"#,
            Some("species is 5, not an Ensembl species name"),
        ),
        (
            r#"{"species": ""}"#,
            Some(r#"species is "", not an Ensembl species name"#),
        ),
        (
            r#"{"species": "9x"}"#,
            Some(r#"species is "9x", not an Ensembl species name"#),
        ),
        (
            "{\"species\": \"h\u{e9}\"}",
            Some("species is \"h\\u00e9\", not an Ensembl species name"),
        ),
        (
            r#"{"cache_layout": "ensembl", "release": "116", "assembly": "GRCh38"}"#,
            None,
        ),
        (
            r#"{"cache_layout": "ensembl", "release": 116, "assembly": "GRCh38"}"#,
            Some("cache_layout ensembl wants release and assembly"),
        ),
        (
            r#"{"cache_layout": "ensembl", "release": "116", "assembly": ""}"#,
            Some("cache_layout ensembl wants release and assembly"),
        ),
        (
            r#"{"cache_layout": "ensembl", "release": "116"}"#,
            Some("cache_layout ensembl wants release and assembly"),
        ),
        (
            r#"{"cache_layout": "perl"}"#,
            Some(r#"cache_layout is "perl", not "ensembl""#),
        ),
        (
            r#"{"cache_layout": null}"#,
            Some(r#"cache_layout is null, not "ensembl""#),
        ),
        (
            r#"{"cache_layout": ["ensembl"]}"#,
            Some(r#"cache_layout is ["ensembl"], not "ensembl""#),
        ),
        (
            r#"{"cache_layout": 1.50}"#,
            Some(r#"cache_layout is 1.5, not "ensembl""#),
        ),
        (
            r#"{"fasta_in_cache": true, "fasta": "reference.fa.gz"}"#,
            None,
        ),
        (r#"{"fasta_in_cache": false}"#, None),
        (
            r#"{"fasta_in_cache": true}"#,
            Some("fasta_in_cache wants fasta"),
        ),
        (
            r#"{"fasta_in_cache": true, "fasta": 5}"#,
            Some("fasta_in_cache wants fasta"),
        ),
        (
            r#"{"fasta_in_cache": "yes"}"#,
            Some(r#"fasta_in_cache is "yes", not a boolean"#),
        ),
        (
            r#"{"fasta_in_cache": 1}"#,
            Some("fasta_in_cache is 1, not a boolean"),
        ),
    ];
    for (text, reason) in cases {
        assert_eq!(shape(text).as_deref(), reason, "{text}");
    }
    let included = [
        (r#"{"s01": ["rs1"]}"#, true),
        ("{}", true),
        (r#"{"s01": []}"#, false),
        (r#"{" ": ["rs1"]}"#, false),
        (r#"["rs1"]"#, false),
        (r#"{"s01": ["rs1", ""]}"#, false),
        ("null", false),
    ];
    for (value, held) in included {
        let reason = shape(&format!(r#"{{"included_records": {value}}}"#));
        let expected =
            (!held).then(|| "included_records is not a map of dataset to record names".to_string());
        assert_eq!(reason, expected, "{value}");
    }
}

#[test]
fn a_flag_set_is_held_to_the_shapes_the_subcommand_writes() {
    let nameless = "has a flag_sets entry without a name (letters, digits, _ and -, not default)";
    let unknown = format!(
        "flag set a masks names unknown classes ['splice_rows'] (the classes are {})",
        classes_text()
    );
    let cases: [(&str, Option<&str>); 39] = [
        (r#""pick""#, Some("flag_sets is not a list")),
        ("null", Some("flag_sets is not a list")),
        (r#"["pick"]"#, Some(nameless)),
        (r#"[{"covers": []}]"#, Some(nameless)),
        (r#"[{"name": "default"}]"#, Some(nameless)),
        (r#"[{"name": "-x"}]"#, Some(nameless)),
        (r#"[{"name": 5}]"#, Some(nameless)),
        (r#"[{"name": "a b"}]"#, Some(nameless)),
        (r#"[{"name": "_a"}]"#, Some(nameless)),
        (r#"[{"name": "a\n"}]"#, Some(nameless)),
        ("[{\"name\": \"\u{e9}\"}]", Some(nameless)),
        (r#"[{"name": "a_b-1"}]"#, None),
        (
            r#"[{"name": "a"}, {"name": "a"}]"#,
            Some("names flag set a twice"),
        ),
        (
            r#"[{"name": "a", "covers": "x"}]"#,
            Some("flag set a covers is not a list of names"),
        ),
        (
            r#"[{"name": "a", "masks": "order_only"}]"#,
            Some("flag set a masks is not a list of names"),
        ),
        (
            r#"[{"name": "a", "masks": ["order_only", "splice_rows"]}]"#,
            Some(&unknown),
        ),
        (r#"[{"name": "a", "masks": []}]"#, None),
        (
            r#"[{"name": "a", "reference_warnings": "no FASTA"}]"#,
            Some("flag set a reference_warnings is not a list of texts"),
        ),
        (
            r#"[{"name": "a", "reference_warnings": [["no FASTA"]]}]"#,
            Some("flag set a reference_warnings is not a list of texts"),
        ),
        (
            r#"[{"name": "a", "reference_warnings": ["no FASTA", " "]}]"#,
            Some("flag set a reference_warnings is not a list of texts"),
        ),
        (
            r#"[{"name": "a", "flags": "--pick"}]"#,
            Some("flag set a flags is not a list of strings"),
        ),
        (
            r#"[{"name": "a", "flags": ["--pick", 1]}]"#,
            Some("flag set a flags is not a list of strings"),
        ),
        (
            r#"[{"name": "a", "input": 5}]"#,
            Some("flag set a input is not a path inside the corpus directory"),
        ),
        (
            r#"[{"name": "a", "input": "../x.vcf"}]"#,
            Some("flag set a input is not a path inside the corpus directory"),
        ),
        (
            r#"[{"name": "a", "input": "/abs"}]"#,
            Some("flag set a input is not a path inside the corpus directory"),
        ),
        (
            r#"[{"name": "a", "input": "a/../b"}]"#,
            Some("flag set a input is not a path inside the corpus directory"),
        ),
        (r#"[{"name": "a", "input": "a/./b"}]"#, None),
        (r#"[{"name": "a", "input": "..x/y"}]"#, None),
        (
            r#"[{"name": "a", "expected_dir": ""}]"#,
            Some("flag set a expected_dir is not a path inside the corpus directory"),
        ),
        (
            r#"[{"name": "a", "stdin": "true"}]"#,
            Some("flag set a stdin is not true or false"),
        ),
        (
            r#"[{"name": "a", "unordered": 1}]"#,
            Some("flag set a unordered is not true or false"),
        ),
        (
            r#"[{"name": "a", "formats": []}]"#,
            Some("flag set a formats is not a non-empty subset of default,tab,vcf,json"),
        ),
        (
            r#"[{"name": "a", "formats": "vcf"}]"#,
            Some("flag set a formats is not a non-empty subset of default,tab,vcf,json"),
        ),
        (
            r#"[{"name": "a", "formats": ["vcf", 1]}]"#,
            Some("flag set a formats is not a non-empty subset of default,tab,vcf,json"),
        ),
        (
            r#"[{"name": "a", "stdin": true, "flags": ["--pick"]}]"#,
            Some("flag set a runs on stdin without --format among its flags"),
        ),
        (
            r#"[{"name": "a", "stdin": true, "flags": ["--format=vcf"]}]"#,
            None,
        ),
        (r#"[{"name": "a", "stdin": false, "flags": []}]"#, None),
        (
            r#"[{"name": "db", "perl_mode": "cache"}]"#,
            Some("flag set db perl_mode is not one of offline, database, gtf, gff"),
        ),
        (r#"[{"name": "db", "perl_mode": "database"}]"#, None),
    ];
    for (sets, reason) in cases {
        let text = format!(r#"{{"covers": ["symbol"], "flag_sets": {sets}}}"#);
        assert_eq!(shape(&text).as_deref(), reason, "{sets}");
    }
}

#[test]
fn a_with_reference_key_is_held_to_its_shape_on_a_divergences_entry_alone() {
    let key = WITH_REFERENCE_SETS;
    let entry = |rest: &str| {
        format!(
            r#"{{"divergences": [{{"location": "21:100", "allele": "G", "feature": "ENST1", {rest}}}]}}"#
        )
    };
    let on = format!("the manifest's divergences entry on 21:100 G ENST1: {key}");
    let cases = [
        (format!(r#""{key}": "x""#), Some("is not a non-empty list")),
        (format!(r#""{key}": []"#), Some("is not a non-empty list")),
        (
            format!(r#""{key}": ["b,a"]"#),
            Some("holds an item that is not a sorted, comma-joined set of terms"),
        ),
        (
            format!(r#""{key}": ["a", ""]"#),
            Some("holds an item that is not a sorted, comma-joined set of terms"),
        ),
        (
            format!(r#""{key}": ["b", "a"]"#),
            Some("is not a sorted list of distinct sets"),
        ),
        (
            format!(r#""{key}": ["a", "a"]"#),
            Some("is not a sorted list of distinct sets"),
        ),
        (
            format!(r#""{key}": ["a", "b"], "vep_rs_consequence_sets": ["a", "b"]"#),
            Some("repeats vep_rs_consequence_sets and states no second outcome"),
        ),
        (
            format!(r#""{key}": ["a", "b"], "vep_consequence_set": "b""#),
            Some("names Ensembl VEP's set beside another"),
        ),
        (
            format!(r#""{key}": ["b"], "vep_consequence_set": "b""#),
            None,
        ),
        (
            format!(
                r#""{key}": ["a,b", "c"], "vep_rs_consequence_sets": ["c"], "vep_consequence_set": "d""#
            ),
            None,
        ),
    ];
    for (rest, reason) in cases {
        let expected = reason.map(|reason| format!("{on} {reason}"));
        assert_eq!(shape(&entry(&rest)), expected, "{rest}");
    }
    let elsewhere = [
        (
            format!(r#"{{"vep_rs_only_tuples": [{{"location": 5, "allele": null, "{key}": ["a"]}}]}}"#),
            format!("the manifest's vep_rs_only_tuples entry on 5 None : {key} belongs on a divergences entry"),
        ),
        (
            format!(
                r#"{{"field_divergences": [{{"location": ["x"], "allele": {{"k": 1.5}}, "feature": true, "{key}": ["a"]}}]}}"#
            ),
            format!(
                "the manifest's field_divergences entry on ['x'] {{'k': 1.5}} True: {key} belongs on a divergences entry"
            ),
        ),
        (
            format!(r#"{{"divergences": "x", "vep_rs_only_tuples": [3, {{"{key}": 1}}]}}"#),
            format!("the manifest's vep_rs_only_tuples entry on   : {key} belongs on a divergences entry"),
        ),
        (
            format!(
                r#"{{"flag_sets": [{{"name": "hgvs", "divergences": [{{"location": "1:5", "{key}": [1]}}]}}]}}"#
            ),
            format!(
                "flag set hgvs divergences entry on 1:5  : {key} holds an item that is not a sorted, comma-joined set of terms"
            ),
        ),
        (
            format!(
                r#"{{"flag_sets": [{{"name": "a"}}, {{"name": "b", "field_divergences": [{{"{key}": ["a"]}}]}}]}}"#
            ),
            format!("flag set b field_divergences entry on   : {key} belongs on a divergences entry"),
        ),
    ];
    for (text, reason) in elsewhere {
        assert_eq!(shape(&text), Some(reason), "{text}");
    }
}

#[test]
fn the_readings_of_a_manifest_follow_its_lists() {
    let manifest = object(
        r#"{"covers": ["symbol"], "flag_sets": [{"name": "af", "covers": ["max_af"]}, {"name": "pick"}]}"#,
    );
    assert_eq!(
        covers_lists(&manifest),
        [
            ("covers".to_string(), vec!["symbol".to_string()]),
            ("flag set af covers".to_string(), vec!["max_af".to_string()]),
            ("flag set pick covers".to_string(), Vec::new()),
        ]
    );
    assert_eq!(flag_set_index(&manifest, "pick"), Ok(1));
    assert_eq!(
        flag_set_index(&manifest, ""),
        Err(failure("the manifest has no flag set ''"))
    );
    assert_eq!(
        target_run(&manifest, Some("af")).map(|run| names_under(run, "covers")),
        Ok(vec!["max_af".to_string()])
    );
    assert_eq!(target_run(&manifest, None), Ok(&manifest));
    let items: HashSet<String> = ["symbol".to_string()].into_iter().collect();
    assert_eq!(
        names_outside_matrix(
            &["ghost".to_string(), "symbol".to_string()],
            &items,
            "m.tsv",
            &covers_owner(Some("af"))
        ),
        ["flag set af covers names ghost, which is not a row of m.tsv"]
    );
    assert_eq!(
        covers_names(&["b".to_string(), "a".to_string(), "b".to_string()]),
        Ok(vec!["b".to_string(), "a".to_string()])
    );
    assert_eq!(covers_names(&[]), Ok(Vec::new()));
    assert_eq!(
        covers_names(&["a".to_string(), " ".to_string()]),
        Err(failure("--covers wants a non-empty row name"))
    );
    assert_eq!(
        path_inside_corpus("../x", "--input"),
        Err(failure(
            "--input wants a path inside the corpus directory, got '../x'"
        ))
    );
    assert_eq!(
        path_inside_corpus("sub/variants.vcf", "--input"),
        Ok("sub/variants.vcf".to_string())
    );
    let options = run_options(&["--hgvs".to_string()], Some("reference.fa.gz"));
    assert_eq!(
        manifest_text(&options),
        "{\n \"flags\": [\n  \"--hgvs\"\n ],\n \"fasta\": \"reference.fa.gz\"\n}\n"
    );
    assert!(run_options(&[], None).is_empty() && run_options(&[], Some("")).is_empty());
    assert!(names_format(&["--pick", "--format=vcf"]) && names_format(&["--format"]));
    assert!(!names_format(&["--formats", "format"]));
}

#[test]
fn a_generated_block_is_a_command_and_the_digests_of_its_files() {
    let digest = "0123456789abcdef".repeat(4);
    let block = |command: &str, digests: &str| {
        let text = format!(r#"{{"command": {command}, "sha256": {digests}}}"#);
        generated_shape_error(&loads(&text).unwrap())
    };
    let one = format!(r#"{{"a": "{digest}"}}"#);
    // A word is any text but the empty one, a path any relative one without a `..` part,
    // and a key beside the two is no fault.
    assert_eq!(block(r#"[" "]"#, &one), None);
    assert_eq!(
        block(
            r#"["a", "{generated}"]"#,
            &format!(
                r#"{{"a b": "{digest}", "d/e/f": "{digest}", "./a": "{digest}", "a/..b": "{digest}"}}"#
            )
        ),
        None
    );
    let noted = format!(r#"{{"command": ["a"], "sha256": {one}, "note": "kept"}}"#);
    assert_eq!(generated_shape_error(&loads(&noted).unwrap()), None);
    assert_eq!(
        shape(&format!(
            r#"{{"covers": ["symbol"], "generated": {noted}}}"#
        )),
        None
    );
    let no_block = "generated is not a block of command and sha256";
    for text in ["null", r#""python3""#, "5", r#"["python3"]"#] {
        let reason = generated_shape_error(&loads(text).unwrap());
        assert_eq!(reason.as_deref(), Some(no_block), "{text}");
    }
    assert_eq!(
        shape(r#"{"generated": null}"#).as_deref(),
        Some(no_block),
        "a block is read wherever the key stands"
    );
    let no_command = "generated command is not a non-empty list of words";
    assert_eq!(
        generated_shape_error(&loads("{}").unwrap()).as_deref(),
        Some(no_command)
    );
    for command in ["null", "[]", r#""a b""#, r#"["a", 5]"#, r#"["a", ""]"#] {
        assert_eq!(
            block(command, &one).as_deref(),
            Some(no_command),
            "{command}"
        );
    }
    let no_map = "generated sha256 is not a map of paths to digests";
    for digests in ["null", "{}", r#"["a"]"#, r#""a""#] {
        assert_eq!(
            block(r#"["a"]"#, digests).as_deref(),
            Some(no_map),
            "{digests}"
        );
    }
    let outside = |digests: &str, listed: &str| {
        assert_eq!(
            block(r#"["a"]"#, digests),
            Some(format!(
                "generated sha256 names paths outside the generated directory {listed}"
            )),
            "{digests}"
        );
    };
    outside(&format!(r#"{{"": "{digest}"}}"#), "['']");
    outside(
        &format!(r#"{{"a/../b": "{digest}", "..": "{digest}", "./a": "{digest}"}}"#),
        "['..', 'a/../b']",
    );
    outside(
        &format!(r#"{{"/": "{digest}", "//a": "{digest}", "a//b": "{digest}"}}"#),
        "['/', '//a']",
    );
    // A path outside the directory is named before a digest of another shape.
    outside(r#"{"../a": 5, "b": 5}"#, "['../a']");
    let undigested = |digests: &str, listed: &str| {
        assert_eq!(
            block(r#"["a"]"#, digests),
            Some(format!(
                "generated sha256 gives no 64 hexadecimal digits for {listed}"
            )),
            "{digests}"
        );
    };
    let long = format!("{digest}0");
    let other = format!("{}G", &digest[..63]);
    undigested(
        &format!(r#"{{"b": "{long}", "a": "{other}", "c": "{digest}"}}"#),
        "['a', 'b']",
    );
    undigested(&format!(r#"{{"a": "{digest}\n"}}"#), "['a']");
    undigested(
        r#"{"a": null, "\u00e9": 1.5, "Z": true}"#,
        "['Z', 'a', '\u{e9}']",
    );
}

#[test]
fn a_hand_written_statement_stops_the_run_that_would_discard_it() {
    let own = object(
        r#"{"divergences": [{"location": "21:100", "allele": "G", "feature": "ENST1", "vep_rs_consequence_sets_with_reference_sequence": ["missense_variant"]}, {"location": "21:200", "allele": "T", "feature": "ENST2"}, "no entry"]}"#,
    );
    let set = object(
        r#"{"name": "pick", "divergences": [{"location": 5, "vep_rs_consequence_sets_with_reference_sequence": null}], "vep_rs_only_tuples": [{"location": "21:300", "vep_rs_consequence_sets_with_reference_sequence": []}]}"#,
    );
    let unlisted =
        object(r#"{"divergences": {"vep_rs_consequence_sets_with_reference_sequence": []}}"#);
    let runs = || {
        [
            ("the manifest's".to_string(), &own),
            ("flag set pick".to_string(), &set),
            ("flag set other".to_string(), &unlisted),
        ]
    };
    let names = [
        "the manifest's divergences entry on 21:100 G ENST1",
        "flag set pick divergences entry on 5  ",
    ];
    assert_eq!(
        reference_statements_to_discard(runs(), "select", true),
        Ok(names.map(str::to_string).to_vec())
    );
    assert_eq!(
        reference_statements_to_discard(runs(), "a flag-set rewrite", false),
        Err(failure(format!(
            "a flag-set rewrite would discard vep_rs_consequence_sets_with_reference_sequence, a key written by hand, on {}; {}: pass --drop-reference-statements to discard it",
            names[0], names[1]
        )))
    );
    // A run without such an entry is refused under neither.
    for allowed in [true, false] {
        assert_eq!(
            reference_statements_to_discard(
                [("flag set other".to_string(), &unlisted)],
                "select",
                allowed
            ),
            Ok(Vec::new())
        );
    }
}
