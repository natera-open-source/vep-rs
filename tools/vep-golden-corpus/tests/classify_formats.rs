// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify --field` with the tab, VCF and JSON outputs of both engines on made-up
//! runs: on a key the default format documents, each field that differs in a further
//! format and that no entry of the key names is written under the name the format gives
//! it, a key without an entry is not read there, and a pair of outputs that cannot be read
//! ends the run with nothing written.

mod classify_support;
mod common;

use classify_support::{
    counted, field_divergence, fields_written, js, key, nothing_written, obj, record, records,
    row_with, transcript, Case,
};
use common::{clean, refused, Run};

const UNEXPLAINED: &str = "unexplained_residual";

const RS1: &str = "21\t100\trs1\tA\tG\t.\tPASS\t.";
const RS2: &str = "21\t200\trs2\tA\tT\t.\tPASS\t.";

/// The line of standard output of a run with `--field` and no key divergence.
fn fields_counted(count: usize) -> String {
    let line = counted((0, "{}"), (0, "{}"), "[]");
    let summary = match count {
        0 => "{}".to_string(),
        _ => format!("{{'{UNEXPLAINED}': {count}}}"),
    };
    format!(
        "{}; field divergences {count} ({summary})\n",
        line.strip_suffix('\n').unwrap()
    )
}

/// The default-format output of one engine: a row on each record, the first with the
/// `PLUGIN` value given.
fn default_rows(plugin: &str) -> String {
    let extra = |plugin: &str| format!("IMPACT=MODERATE;STRAND=1;PLUGIN={plugin}");
    row_with(
        "rs1",
        "21:100",
        "G",
        "ENST1",
        "missense_variant",
        &extra(plugin),
    ) + &row_with(
        "rs2",
        "21:200",
        "T",
        "ENST2",
        "intron_variant",
        &extra("x|1|2"),
    )
}

/// The tab-format output of one engine: the same two rows with the `OTHER` column given
/// on both, a column the run does not compare.
fn tab_rows(plugin: &str, other: &str) -> String {
    let header = "## a made-up tab output\n#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tPLUGIN\tOTHER\n";
    let row = |name: &str,
               location: &str,
               allele: &str,
               feature: &str,
               terms: &str,
               plugin: &str| {
        format!("{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t{terms}\t{plugin}\t{other}\n")
    };
    header.to_string()
        + &row("rs1", "21:100", "G", "ENST1", "missense_variant", plugin)
        + &row("rs2", "21:200", "T", "ENST2", "intron_variant", "x|1|2")
}

/// The VCF-format output of one engine: the three fields the format writes the plugin's
/// value as, on the first record's entry as given and on the second's as `second`.
fn vcf_lines(first: [&str; 3], second: [&str; 3]) -> String {
    let header = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|PLUGIN_NAME|PLUGIN_ONE|PLUGIN_TWO\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let line = |record: &str, kept: &str, key: &str, fields: [&str; 3]| {
        let columns = record.strip_suffix(".").unwrap();
        format!("{columns}{kept}CSQ={key}|{}\n", fields.join("|"))
    };
    header.to_string()
        + &line(RS1, "AC=1;", "G|missense_variant|ENST1", first)
        + &line(RS2, "", "T|intron_variant|ENST2", second)
}

/// The JSON-format output of one engine: a consequence object on each record with the
/// members given, each as JSON text.
fn json_lines(first: &str, second: &str) -> String {
    let object = |input: &str, allele: &str, feature: &str, terms: &str, members: &str| {
        format!(
            r#"{{"input": {}, "transcript_consequences": [{{"variant_allele": "{allele}", "transcript_id": "{feature}", "consequence_terms": ["{terms}"], {members}}}]}}"#,
            js(input)
        ) + "\n"
    };
    object(RS1, "G", "ENST1", "missense_variant", first)
        + &object(RS2, "T", "ENST2", "intron_variant", second)
}

/// A case of two records whose first row differs in `PLUGIN` in the default format, with
/// the outputs of every further format beside it.
fn four_format_case(tag: &str) -> (Case, String) {
    let case = Case::new(tag);
    let held = records(&[record(RS1), record(RS2)]);
    case.manifest(&obj(&[&held]));
    case.outputs(&default_rows("a,b|1|2"), &default_rows("a|1|3"));
    case.write("vep.tab", &tab_rows("a,b|1|2", "k"));
    case.write("rs.tab", &tab_rows("a|1|3", "m"));
    // The second record's entry differs in every further format and has no entry in the
    // default one.
    case.write("vep.vcf", &vcf_lines(["a&b", "", "2"], ["x", "1", "2"]));
    case.write("rs.vcf", &vcf_lines(["a", "1", "3"], ["y", "1", "2"]));
    case.write(
        "vep.json",
        &json_lines(
            r#""plugin_values": {"TWO": 2, "NAME": "a,b"}, "other": "k""#,
            r#""plugin_values": {"NAME": "x"}"#,
        ),
    );
    case.write(
        "rs.json",
        &json_lines(
            r#""plugin_values": {"TWO": 3, "NAME": "a"}, "other": "m", "score": 0.10"#,
            r#""plugin_values": {"NAME": "y"}"#,
        ),
    );
    (case, held)
}

/// The options that name the outputs of the formats given.
fn further(case: &Case, formats: &[&str]) -> Vec<String> {
    let mut options = Vec::new();
    for format in formats {
        options.push(format!("--vep-{format}"));
        options.push(case.path(&format!("vep.{format}")));
        options.push(format!("--vep-rs-{format}"));
        options.push(case.path(&format!("rs.{format}")));
    }
    options
}

fn classify(case: &Case, options: &[String]) -> Run {
    let mut arguments = vec!["--field", "PLUGIN"];
    arguments.extend(options.iter().map(String::as_str));
    case.classify(&arguments)
}

/// An entry of the first record's row.
fn entry(field: &str, vep: &str, vep_rs: &str) -> String {
    field_divergence(
        &key("21:100", "G", transcript("ENST1"), &[0]),
        "missense_variant",
        (field, vep, vep_rs),
        UNEXPLAINED,
    )
}

#[test]
fn classify_writes_a_field_under_the_name_each_format_gives_it() {
    let (case, held) = four_format_case("formats");
    let compared = entry("PLUGIN", "a,b|1|2", "a|1|3");
    let expected = |entries: &[String]| {
        let summary = match entries.len() {
            0 => "{}".to_string(),
            count => format!(r#"{{"{UNEXPLAINED}": {count}}}"#),
        };
        obj(&[
            &held,
            &nothing_written(),
            &fields_written(&["PLUGIN"], entries, &summary),
        ])
    };
    // The default format alone: the one compared field.
    assert_eq!(classify(&case, &[]), clean(0, &fields_counted(1)));
    assert_eq!(case.after(), expected(std::slice::from_ref(&compared)));

    // Every format: behind the compared field, the fields of each further format that
    // differ on the row and that no entry names yet, the formats in the order tab, VCF,
    // JSON and the fields of one in the byte order of their names. A value is spelt as the
    // default format spells it: `,` for the VCF format's `&`, `-` for its empty value and
    // for a member a JSON object lacks; a JSON object is written without blanks, its
    // members in the byte order of their names; `other` is the tab format's `OTHER`.
    let all = further(&case, &["tab", "vcf", "json"]);
    assert_eq!(classify(&case, &all), clean(0, &fields_counted(7)));
    let entries = [
        compared.clone(),
        entry("OTHER", "k", "m"),
        entry("PLUGIN_NAME", "a,b", "a"),
        entry("PLUGIN_ONE", "-", "1"),
        entry("PLUGIN_TWO", "2", "3"),
        entry(
            "plugin_values",
            r#"{"NAME":"a,b","TWO":2}"#,
            r#"{"NAME":"a","TWO":3}"#,
        ),
        entry("score", "-", "0.1"),
    ];
    assert_eq!(case.after(), expected(&entries));

    // The formats named are the formats read, and a field is written once: without the
    // tab format the JSON format's `other` is the first to name the field.
    let two = further(&case, &["json", "vcf"]);
    assert_eq!(classify(&case, &two), clean(0, &fields_counted(7)));
    let entries = [
        compared.clone(),
        entry("PLUGIN_NAME", "a,b", "a"),
        entry("PLUGIN_ONE", "-", "1"),
        entry("PLUGIN_TWO", "2", "3"),
        entry("other", "k", "m"),
        entry(
            "plugin_values",
            r#"{"NAME":"a,b","TWO":2}"#,
            r#"{"NAME":"a","TWO":3}"#,
        ),
        entry("score", "-", "0.1"),
    ];
    assert_eq!(case.after(), expected(&entries));

    // A row without a difference in the default format is read in no further format.
    case.outputs(&default_rows("a|1|3"), &default_rows("a|1|3"));
    assert_eq!(classify(&case, &all), clean(0, &fields_counted(0)));
    assert_eq!(case.after(), expected(&[]));
}

#[test]
fn classify_compares_a_multiset_column_of_a_further_format_as_one() {
    // Under `order_only` the items of `DOMAINS` compare as a multiset in the default
    // format, and so do the JSON format's `domains` on a row with an entry.
    let (case, held) = four_format_case("formats-order");
    let masks = r#""masks": ["order_only"]"#;
    case.manifest(&obj(&[&held, masks]));
    let domains = |first: &str, second: &str| {
        format!(
            r#""plugin_values": {{"NAME": "a"}}, "domains": [{{"db": "{first}", "name": "n"}}, {{"db": "{second}", "name": "n"}}]"#
        )
    };
    let unchanged = r#""plugin_values": {"NAME": "x"}"#;
    case.write("vep.json", &json_lines(&domains("p", "q"), unchanged));
    case.write("rs.json", &json_lines(&domains("q", "p"), unchanged));
    let json = further(&case, &["json"]);
    assert_eq!(classify(&case, &json), clean(0, &fields_counted(1)));
    case.manifest(&obj(&[&held]));
    assert_eq!(classify(&case, &json), clean(0, &fields_counted(2)));
    let entries = [
        entry("PLUGIN", "a,b|1|2", "a|1|3"),
        entry(
            "domains",
            r#"{"db":"p","name":"n"},{"db":"q","name":"n"}"#,
            r#"{"db":"q","name":"n"},{"db":"p","name":"n"}"#,
        ),
    ];
    let summary = format!(r#"{{"{UNEXPLAINED}": 2}}"#);
    assert_eq!(
        case.after(),
        obj(&[
            &held,
            &nothing_written(),
            &fields_written(&["PLUGIN"], &entries, &summary),
        ])
    );
}

#[test]
fn classify_pairs_the_outputs_of_a_reencoded_input_by_its_lines_and_in_turn() {
    // A flag set over `variants.vcf` re-encoded as Ensembl-format lines: a JSON object names
    // its record by the input's own line, and the lines of the VCF output, which carry no
    // input line, are the records in turn. The entries are the manifest's own run's, written
    // into the set.
    let (case, held) = four_format_case("formats-reencoded");
    let lines = ["21\t100\t100\tA/G\t+\trs1", "21\t200\t200\tA/T\t+\trs2"];
    case.write_in_corpus("variants.ensembl", (lines.join("\n") + "\n").as_bytes());
    let set = |more: &str| {
        format!(
            r#""flag_sets": [{{"name": "ensembl", "flags": ["--format", "ensembl"], "input": "variants.ensembl"{more}}}]"#
        )
    };
    case.manifest(&obj(&[&held, &set("")]));
    let reencode = |name: &str| {
        let text = std::fs::read_to_string(case.path(name)).unwrap();
        let text = text
            .replace(&js(RS1), &js(lines[0]))
            .replace(&js(RS2), &js(lines[1]));
        // No line of the VCF output is a line of the input.
        case.write(name, &text.replace("\trs1\tA\tG\t", "\t.\tN\tG\t"));
    };
    for name in ["vep.json", "rs.json", "vep.vcf", "rs.vcf"] {
        reencode(name);
    }
    let mut options = further(&case, &["vcf", "json"]);
    options.extend(["--flag-set".to_string(), "ensembl".to_string()]);
    assert_eq!(classify(&case, &options), clean(0, &fields_counted(7)));
    let entries = [
        entry("PLUGIN", "a,b|1|2", "a|1|3"),
        entry("PLUGIN_NAME", "a,b", "a"),
        entry("PLUGIN_ONE", "-", "1"),
        entry("PLUGIN_TWO", "2", "3"),
        entry("other", "k", "m"),
        entry(
            "plugin_values",
            r#"{"NAME":"a,b","TWO":2}"#,
            r#"{"NAME":"a","TWO":3}"#,
        ),
        entry("score", "-", "0.1"),
    ];
    let summary = format!(r#"{{"{UNEXPLAINED}": 7}}"#);
    let written = format!(
        ", {}, {}",
        nothing_written(),
        fields_written(&["PLUGIN"], &entries, &summary)
    );
    assert_eq!(case.after(), obj(&[&held, &set(&written)]));
    // A VCF output of another number of lines than the input has records is refused.
    let vcf = case.path("rs.vcf");
    let text = std::fs::read_to_string(&vcf).unwrap();
    case.write(
        "rs.vcf",
        text.trim_end_matches('\n').rsplit_once('\n').unwrap().0,
    );
    let before = case.text();
    assert_eq!(
        classify(&case, &options),
        refused(
            1,
            &format!(
                "ERROR: [vep-golden-corpus] {vcf} holds 1 lines for the 2 records of an input that is no VCF, whose lines are the records in turn"
            )
        )
    );
    assert_eq!(case.text(), before);
}

#[test]
fn classify_refuses_outputs_of_a_further_format_it_cannot_read_with_nothing_written() {
    let (case, _) = four_format_case("formats-refused");
    let before = case.text();
    let own = |message: &str| refused(1, &format!("ERROR: [vep-golden-corpus] {message}"));
    // One output of a format without the other.
    let mut alone = further(&case, &["vcf"]);
    alone.truncate(2);
    assert_eq!(
        classify(&case, &alone),
        own("--vep-vcf and --vep-rs-vcf name the two outputs of one format: give both or neither")
    );
    // The further formats are read on the keys `--field` documents.
    let options = further(&case, &["tab"]);
    let options: Vec<&str> = options.iter().map(String::as_str).collect();
    assert_eq!(
        case.classify(&options),
        own("the outputs of a further format are read on the keys --field documents: give --field")
    );
    // An entry of another width than the header names, a line that is no object, a row
    // before the column line and a file that is not there.
    let vcf = case.path("rs.vcf");
    let text = std::fs::read_to_string(&vcf).unwrap();
    case.write("rs.vcf", &text.replace("|a|1|3", "|a|1"));
    assert_eq!(
        classify(&case, &further(&case, &["vcf"])),
        own(&format!(
            "{vcf} holds an entry of 5 fields under a header that names 6 for CSQ"
        ))
    );
    case.write("rs.json", "[]\n");
    assert_eq!(
        classify(&case, &further(&case, &["json"])),
        own(&format!(
            "{} line 1 is no JSON object",
            case.path("rs.json")
        ))
    );
    case.write("vep.tab", "rs1\t21:100\n");
    assert_eq!(
        classify(&case, &further(&case, &["tab"])),
        own(&format!(
            "{} holds a row before the line that names its columns",
            case.path("vep.tab")
        ))
    );
    std::fs::remove_file(case.path("vep.tab")).unwrap();
    let missing = classify(&case, &further(&case, &["tab"]));
    let line = format!(
        "ERROR: [vep-golden-corpus] {} is not a readable text file (",
        case.path("vep.tab")
    );
    assert_eq!((missing.code, missing.stdout.as_str()), (1, ""));
    assert!(missing.stderr.starts_with(&line), "{}", missing.stderr);
    assert_eq!(case.text(), before);
}
