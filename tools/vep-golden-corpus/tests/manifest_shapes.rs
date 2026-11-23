// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs every subcommand that reads a manifest over manifests of each refused shape: the
//! run ends on the reason before anything is written, and the shapes beside them pass.

mod common;

use common::{
    clean, failed, fixture, manifest_file, manifest_text, matrix, object, run, scratch,
    two_record_corpus,
};
use vep_golden_corpus::constants::{mask_classes, HARNESS_ONLY_MASK_CLASSES};
use vep_golden_corpus::json::{dumps, dumps_object, loads, Style, Value};

/// The command lines that read the manifest of `corpus`, one a subcommand.
fn commands(corpus: &str, matrix: &str) -> Vec<Vec<String>> {
    let command = |args: &[&str]| args.iter().map(|arg| arg.to_string()).collect();
    let reference = fixture("two_record/ref.txt");
    let reference = reference.to_str().unwrap();
    let dataset = format!(
        "pool={reference}={}",
        fixture("two_record/in.vcf").to_str().unwrap()
    );
    vec![
        command(&[
            "select",
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--dataset",
            &dataset,
            "--out",
            corpus,
            "--k",
            "1",
            "--seed",
            "1",
            "--multi-allelic-per-dataset",
            "0",
        ]),
        command(&["check", "--corpus", corpus, "--dataset", &dataset]),
        command(&["covers", corpus, "--list", "--matrix", matrix]),
        command(&["flag-set", corpus, "--list", "--matrix", matrix]),
        command(&[
            "classify",
            "--corpus",
            corpus,
            "--vep-default",
            reference,
            "--vep-rs-output",
            reference,
        ]),
        command(&["generated", "--corpus", corpus, "--check"]),
    ]
}

#[test]
fn manifest_shapes_are_refused_by_every_subcommand() {
    let scratch = scratch("shapes");
    let out = two_record_corpus(&scratch.path);
    let matrix = matrix(&scratch.path, &["symbol"]);
    let manifest = manifest_file(&out);
    let named = manifest.to_str().unwrap();
    let good = manifest_text(&out);
    let commands = commands(&out, &matrix);
    let classes = mask_classes().join(", ");
    let unknown =
        format!("masks names unknown classes ['splice_rows'] (the classes are {classes})");
    let unknown_in_set = format!("flag set a {unknown}");
    let nameless = "has a flag_sets entry without a name (letters, digits, _ and -, not default)";
    let shapes: [(&str, &str); 46] = [
        (r#""covers": "symbol""#, "covers is not a list of names"),
        (r#""covers": [["symbol"]]"#, "covers is not a list of names"),
        (
            r#""covers": ["symbol", ""]"#,
            "covers is not a list of names",
        ),
        (
            r#""covers": ["symbol"], "masks": "regulatory_feature_rows""#,
            "masks is not a list of names",
        ),
        (
            r#""covers": ["symbol"], "masks": [["regulatory_feature_rows"]]"#,
            "masks is not a list of names",
        ),
        (
            r#""covers": ["symbol"], "masks": ["order_only", "splice_rows"]"#,
            &unknown,
        ),
        (
            r#""covers": ["symbol"], "reference_warnings": "no FASTA""#,
            "reference_warnings is not a list of texts",
        ),
        (
            r#""covers": ["symbol"], "reference_warnings": ["no FASTA", ""]"#,
            "reference_warnings is not a list of texts",
        ),
        (
            r#""covers": ["symbol"], "plugin_files": ["plugin_data/gerp.bw"]"#,
            "plugin_files is not a map of plugin names",
        ),
        (
            r#""covers": ["symbol"], "plugin_files": {"LoFTEE": ["../gerp.bw"]}"#,
            "plugin_files for LoFTEE names paths outside the corpus ['../gerp.bw']",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": "pick""#,
            "flag_sets is not a list",
        ),
        (r#""covers": ["symbol"], "flag_sets": ["pick"]"#, nameless),
        (
            r#""covers": ["symbol"], "flag_sets": [{"covers": []}]"#,
            nameless,
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "default"}]"#,
            nameless,
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "-x"}]"#,
            nameless,
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a"}, {"name": "a"}]"#,
            "names flag set a twice",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "covers": "x"}]"#,
            "flag set a covers is not a list of names",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "masks": "order_only"}]"#,
            "flag set a masks is not a list of names",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "masks": ["order_only", "splice_rows"]}]"#,
            &unknown_in_set,
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "reference_warnings": "no FASTA"}]"#,
            "flag set a reference_warnings is not a list of texts",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "reference_warnings": [["no FASTA"]]}]"#,
            "flag set a reference_warnings is not a list of texts",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "reference_warnings": ["no FASTA", " "]}]"#,
            "flag set a reference_warnings is not a list of texts",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "flags": "--pick"}]"#,
            "flag set a flags is not a list of strings",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "input": 5}]"#,
            "flag set a input is not a path inside the corpus directory",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "input": "../x.vcf"}]"#,
            "flag set a input is not a path inside the corpus directory",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "expected_dir": ""}]"#,
            "flag set a expected_dir is not a path inside the corpus directory",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "stdin": "true"}]"#,
            "flag set a stdin is not true or false",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "unordered": 1}]"#,
            "flag set a unordered is not true or false",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "formats": []}]"#,
            "flag set a formats is not a non-empty subset of default,tab,vcf,json",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "formats": "vcf"}]"#,
            "flag set a formats is not a non-empty subset of default,tab,vcf,json",
        ),
        (
            r#""covers": ["symbol"], "flag_sets": [{"name": "a", "stdin": true, "flags": ["--pick"]}]"#,
            "flag set a runs on stdin without --format among its flags",
        ),
        (
            r#""covers": ["symbol"], "cache_layout": "perl""#,
            r#"cache_layout is "perl", not "ensembl""#,
        ),
        (
            r#""covers": ["symbol"], "cache_layout": null"#,
            r#"cache_layout is null, not "ensembl""#,
        ),
        (
            r#""covers": ["symbol"], "cache_layout": "ensembl", "release": "116""#,
            "cache_layout ensembl wants release and assembly",
        ),
        (
            r#""covers": ["symbol"], "fasta_in_cache": "yes""#,
            r#"fasta_in_cache is "yes", not a boolean"#,
        ),
        (
            r#""covers": ["symbol"], "fasta_in_cache": true"#,
            "fasta_in_cache wants fasta",
        ),
        (
            r#""covers": ["symbol"], "keep_samples": "8""#,
            "keep_samples is not a count",
        ),
        (
            r#""covers": ["symbol"], "keep_samples": -1"#,
            "keep_samples is not a count",
        ),
        (
            r#""covers": ["symbol"], "keep_samples": true"#,
            "keep_samples is not a count",
        ),
        (
            r#""covers": ["symbol"], "keep_samples": 2"#,
            "samples does not name 2 kept samples",
        ),
        (
            r#""covers": ["symbol"], "keep_samples": 2, "samples": ["a"]"#,
            "samples does not name 2 kept samples",
        ),
        (
            r#""covers": ["symbol"], "keep_samples": 1, "samples": [""]"#,
            "samples does not name 1 kept samples",
        ),
        (
            r#""covers": ["symbol"], "required_genotype_classes": "het""#,
            "required_genotype_classes is not a list of genotype classes",
        ),
        (
            r#""covers": ["symbol"], "required_genotype_classes": ["het", "ghost"]"#,
            "required_genotype_classes is not a list of genotype classes",
        ),
        (
            r#""covers": ["symbol"], "species": "Homo sapiens""#,
            r#"species is "Homo sapiens", not an Ensembl species name"#,
        ),
        (
            r#""covers": ["symbol"], "record_class_key": "motif_class""#,
            "record_class_key is not one of regulatory_class",
        ),
    ];
    // A `generated` block is a command and the digests of its files.
    let digest = format!(r#""{}""#, "0".repeat(64));
    let short = "0".repeat(63);
    let upper = "A".repeat(64);
    let generated: [(String, &str); 10] = [
        (
            r#"["python3"]"#.to_string(),
            "generated is not a block of command and sha256",
        ),
        (
            format!(r#"{{"sha256": {{"a": {digest}}}}}"#),
            "generated command is not a non-empty list of words",
        ),
        (
            format!(r#"{{"command": [], "sha256": {{"a": {digest}}}}}"#),
            "generated command is not a non-empty list of words",
        ),
        (
            format!(r#"{{"command": "python3 write.py", "sha256": {{"a": {digest}}}}}"#),
            "generated command is not a non-empty list of words",
        ),
        (
            format!(r#"{{"command": ["python3", ""], "sha256": {{"a": {digest}}}}}"#),
            "generated command is not a non-empty list of words",
        ),
        (
            r#"{"command": ["python3"]}"#.to_string(),
            "generated sha256 is not a map of paths to digests",
        ),
        (
            r#"{"command": ["python3"], "sha256": {}}"#.to_string(),
            "generated sha256 is not a map of paths to digests",
        ),
        (
            r#"{"command": ["python3"], "sha256": ["a"]}"#.to_string(),
            "generated sha256 is not a map of paths to digests",
        ),
        (
            format!(
                r#"{{"command": ["python3"], "sha256": {{"../a": {digest}, "/b": {digest}, "c": {digest}}}}}"#
            ),
            "generated sha256 names paths outside the generated directory ['../a', '/b']",
        ),
        (
            format!(
                r#"{{"command": ["python3"], "sha256": {{"d": "{short}", "c": {digest}, "b": 5, "a": "{upper}"}}}}"#
            ),
            "generated sha256 gives no 64 hexadecimal digits for ['a', 'b', 'd']",
        ),
    ];
    let generated = generated.iter().map(|(block, reason)| {
        (
            format!(r#""covers": ["symbol"], "generated": {block}"#),
            *reason,
        )
    });
    let shapes = shapes
        .iter()
        .map(|(body, reason)| (body.to_string(), *reason));
    for (body, reason) in generated.chain(shapes) {
        std::fs::write(&manifest, format!(r#"{{"records": [], {body}}}"#)).unwrap();
        for command in &commands {
            assert_eq!(
                run(command),
                failed(&format!("{named} {reason}")),
                "{} on {body}",
                command[0]
            );
        }
    }
    // Every mask class passes at both levels, the harness's own among them, and so does an
    // empty list on a set, whose run masks nothing.
    let mut accepted = object(&good);
    accepted.insert(
        "masks",
        Value::strings(&["order_only", "regulatory_feature_rows"]),
    );
    let harness_classes = dumps(&Value::strings(&HARNESS_ONLY_MASK_CLASSES), &Style::PLAIN);
    let sets = format!(
        r#"[{{"name": "a", "masks": ["order_only"]}}, {{"name": "b", "masks": []}}, {{"name": "c", "masks": {harness_classes}}}]"#
    );
    accepted.insert("flag_sets", loads(&sets).unwrap());
    assert_eq!(
        HARNESS_ONLY_MASK_CLASSES,
        [
            "unloaded_module_plugin",
            "conservation_short_name_spelling",
            "pc_description_text",
            "file_cache_writes_no_source_field",
        ]
    );
    let on_one_line = dumps_object(&accepted, &Style::PLAIN);
    std::fs::write(&manifest, on_one_line).unwrap();
    assert_eq!(
        run(&commands[1]),
        clean(0, "observed 2, covered by corpus 2, missing 0\n")
    );
    assert_eq!(run(&commands[2]), clean(0, "symbol\n"));
    assert_eq!(
        run(&commands[3]),
        clean(
            0,
            concat!(
                "a\t\tvariants.vcf\tfalse\tfalse\tdefault,tab,vcf,json\texpected_a\t\n",
                "b\t\tvariants.vcf\tfalse\tfalse\tdefault,tab,vcf,json\texpected_b\t\n",
                "c\t\tvariants.vcf\tfalse\tfalse\tdefault,tab,vcf,json\texpected_c\t\n",
            )
        )
    );
    std::fs::write(&manifest, "[]").unwrap();
    for command in &commands {
        assert_eq!(
            run(command),
            failed(&format!("{named} is not a JSON object"))
        );
    }
    std::fs::write(&manifest, "{").unwrap();
    for command in &commands {
        assert_eq!(
            run(command),
            failed(&format!(
                "{named} is not readable JSON (Expecting property name enclosed in double quotes: line 1 column 2 (char 1))"
            ))
        );
    }
    std::fs::write(&manifest, &good).unwrap();
    // `--flag-set ''` names no set rather than the manifest's own list.
    assert_eq!(
        run(&[
            "covers",
            &out,
            "--flag-set",
            "",
            "--list",
            "--matrix",
            &matrix
        ]),
        failed("the manifest has no flag set ''")
    );
    let none = scratch.path.join("none.tsv");
    let none = none.to_str().unwrap();
    assert_eq!(
        run(&["covers", &out, "--add", "symbol", "--matrix", none]),
        failed(&format!(
            "{none} is not a file: name the parity matrix with --matrix"
        ))
    );
    assert_eq!(
        run(&["covers", &out, "--add", "symbol", "--matrix", ""]),
        failed(". is not a file: name the parity matrix with --matrix")
    );
    let mut check = commands[1].clone();
    check.extend(["--matrix".to_string(), String::new()]);
    assert_eq!(
        run(&check),
        failed(". is not a file: name the parity matrix with --matrix")
    );
    let binary = scratch.path.join("binary.tsv");
    std::fs::write(&binary, b"item\t\xff\xfe\n").unwrap();
    let binary = binary.to_str().unwrap();
    assert_eq!(
        run(&["covers", &out, "--add", "symbol", "--matrix", binary]),
        failed(&format!(
            "{binary} is not a readable text file ('utf-8' codec can't decode byte 0xff in position 5: invalid start byte)"
        ))
    );
    assert_eq!(manifest_text(&out), good);
}

#[test]
fn an_unreadable_manifest_or_matrix_is_refused_with_the_reason() {
    let scratch = scratch("unreadable");
    let out = two_record_corpus(&scratch.path);
    let absent = scratch.path.join("absent");
    let absent = absent.to_str().unwrap();
    assert_eq!(
        run(&["covers", absent, "--list"]),
        failed(&format!(
            "{absent}/manifest.json is not readable JSON ([Errno 2] No such file or directory: '{absent}/manifest.json')"
        ))
    );
    // A corpus named with a trailing separator is the same corpus, and the empty name is the
    // working directory.
    assert_eq!(run(&["covers", &format!("{out}/./")]), clean(0, "symbol\n"));
    assert_eq!(
        run(&["flag-set", "", "--list"]),
        failed(
            "manifest.json is not readable JSON ([Errno 2] No such file or directory: 'manifest.json')"
        )
    );
    let manifest = manifest_file(&out);
    let named = manifest.to_str().unwrap();
    std::fs::write(&manifest, b"{\"covers\": [\"caf\xc3\"]}").unwrap();
    assert_eq!(
        run(&["covers", &out]),
        failed(&format!(
            "{named} is not readable JSON ('utf-8' codec can't decode byte 0xc3 in position 16: invalid continuation byte)"
        ))
    );
    // A manifest with carriage returns reads as one with newlines, and its lines are counted so.
    std::fs::write(
        &manifest,
        "{\r\n \"covers\": [\r\n  \"symbol\",\r\n ]\r\n}\r\n",
    )
    .unwrap();
    assert_eq!(
        run(&["covers", &out]),
        failed(&format!(
            "{named} is not readable JSON (Illegal trailing comma before end of array: line 3 column 11 (char 25))"
        ))
    );
    std::fs::write(&manifest, "{\"covers\": [\"symbol\"]}").unwrap();
    let no_header = scratch.path.join("no_header.tsv");
    std::fs::write(&no_header, "symbol\tConfig.pm\n").unwrap();
    let no_header = no_header.to_str().unwrap();
    assert_eq!(
        run(&["covers", &out, "--add", "symbol", "--matrix", no_header]),
        failed(&format!(
            "{no_header} is not a parity matrix (first column must be item)"
        ))
    );
    let empty = scratch.path.join("empty.tsv");
    std::fs::write(&empty, "").unwrap();
    let empty = empty.to_str().unwrap();
    assert_eq!(
        run(&["covers", &out, "--add", "symbol", "--matrix", empty]),
        failed(&format!(
            "{empty} is not a parity matrix (first column must be item)"
        ))
    );
    // A row of the matrix is its first column up to a line boundary of `str.splitlines`; a
    // blank line names no row.
    let odd = scratch.path.join("odd.tsv");
    std::fs::write(
        &odd,
        "item\tx\r\nsymbol\r\n\r\n \t \r\nmax_af\tx\x0bpick\x0cper_gene",
    )
    .unwrap();
    let odd = odd.to_str().unwrap();
    assert_eq!(
        run(&[
            "covers", &out, "--add", "max_af", "--add", "pick", "--add", "per_gene", "--add", "",
            "--matrix", odd
        ]),
        failed("--covers wants a non-empty row name")
    );
    assert_eq!(
        run(&[
            "covers", &out, "--add", "max_af", "--add", "pick", "--add", "per_gene", "--add", " ",
            "--add", "ghost", "--matrix", odd
        ]),
        failed("--covers wants a non-empty row name")
    );
    assert_eq!(
        run(&[
            "covers", &out, "--add", "max_af", "--add", "pick", "--add", "per_gene", "--add",
            "ghost", "--list", "--matrix", odd
        ]),
        failed(&format!("covers names ghost, which is not a row of {odd}"))
    );
    assert_eq!(
        run(&[
            "covers", &out, "--add", "max_af", "--add", "pick", "--add", "per_gene", "--list",
            "--matrix", odd
        ]),
        clean(0, "symbol\nmax_af\npick\nper_gene\n")
    );
    assert_eq!(
        manifest_text(&out),
        "{\n \"covers\": [\n  \"symbol\",\n  \"max_af\",\n  \"pick\",\n  \"per_gene\"\n ]\n}\n"
    );
}
