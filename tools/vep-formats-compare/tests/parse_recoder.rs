// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The recoder parser on two `variant_recoder` outputs and against stored answers.
//!
//! `tests/data/recoder/cases.tsv` holds `masks<TAB>text<TAB>answer` lines written by
//! CPython, a text as the lowercase hex of its UTF-8 bytes: the answer is `E` and the
//! comparator's error text, or `K` and `json.dumps` of the entries (key, allele, fields in
//! order) and the applied counts in the order they were first counted. The ignored test
//! reads the file of the same name, with many more lines, under the `recoder` folder of the
//! directory `VEP_FORMATS_COMPARE_PYDIFF_DIR` names.

use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::{Error, Uncaught};
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::recoder::{parse_recoder, recoder_identifier_indel, recoder_value};
use vep_formats_compare::pyjson::{self, DumpOptions, Value};

const BUILT_IN: &[&str] = &["refseq_transcript_entries", "lrg_transcript_entries"];
const ARCHIVE: &str = "archive_dbsnp_synonyms";
const IDENTIFIER_INDEL: &str = "identifier_indel_transcript_notations";

/// The reference's output: keys and values in hash order, with the RefSeq and LRG entries
/// of a merged database.
const EXPECTED: &str = r#"[
 {"T": {
   "vcf_string": ["21-25585733-C-T"],
   "hgvsg": ["NC_000021.9:g.25585733C>T"],
   "hgvsc": ["ENST00000352957.9:c.991G>A", "NM_017446.4:c.991G>A", "ENST00000307301.12:c.*18G>A"],
   "hgvsp": ["NP_059142.3:p.Ala331Thr", "ENSP00000284967.7:p.Ala331Thr", "LRG_1p1:p.Ala331Thr"],
   "ga4gh_vrs": [{
     "state": {"sequence": "T", "type": "SequenceState"},
     "location": {
       "interval": {"end": 25585733, "start": "25585732", "type": "SimpleInterval"},
       "sequence_id": "refseq:NC_000021.9",
       "type": "SequenceLocation"},
     "type": "Allele"}],
   "input": "rs142513484",
   "spdi": ["NC_000021.9:25585732:C:T"],
   "mane_select": [
     {"hgvsc": "ENST00000352957.9:c.991G>A", "hgvsg": "NC_000021.9:g.25585733C>T", "hgvsp": "ENSP00000284967.7:p.Ala331Thr"},
     {"hgvsg": "NC_000021.9:g.25585733C>T", "hgvsc": "NM_017446.4:c.991G>A", "hgvsp": "NP_059142.3:p.Ala331Thr"}],
   "var_synonyms": ["ClinVar: VCV000829541,RCV003317418", "OMIM: 606870.0005"],
   "id": ["rs142513484"]}},
 null,
 {"warnings": ["No variant found with ID 'rs_nope'"]}
]"#;

/// The candidate's output: sorted, from an Ensembl cache.
const ACTUAL: &str = r#"[
 {"T": {
   "ga4gh_vrs": [{
     "location": {
       "interval": {"end": 25585733, "start": "25585732", "type": "SimpleInterval"},
       "sequence_id": "refseq:NC_000021.9",
       "type": "SequenceLocation"},
     "state": {"sequence": "T", "type": "SequenceState"},
     "type": "Allele"}],
   "hgvsc": ["ENST00000307301.12:c.*18G>A", "ENST00000352957.9:c.991G>A"],
   "hgvsg": ["NC_000021.9:g.25585733C>T"],
   "hgvsp": ["ENSP00000284967.7:p.Ala331Thr"],
   "id": ["rs142513484"],
   "input": "rs142513484",
   "mane_select": [
     {"hgvsc": "ENST00000352957.9:c.991G>A", "hgvsg": "NC_000021.9:g.25585733C>T", "hgvsp": "ENSP00000284967.7:p.Ala331Thr"}],
   "spdi": ["NC_000021.9:25585732:C:T"],
   "var_synonyms": ["ClinVar: RCV003317418,VCV000829541", "OMIM: 606870.0005"],
   "vcf_string": ["21-25585733-C-T"]}},
 null,
 {"warnings": ["No variant found with ID 'rs_nope'"]}
]"#;

/// An identifier input's deletion allele, its transcript notations unshifted.
const DELETION: &str = r#"{
 "input": "rs2057997486",
 "hgvsg": ["NC_000021.9:g.34886862del"],
 "hgvsc": ["ENST00000300305.7:c.332del"],
 "hgvsp": ["ENSP00000300305.3:p.Leu112CysfsTer10"],
 "mane_select": [{"hgvsc": "ENST00000675419.1:c.332del", "hgvsg": "NC_000021.9:g.34886862del", "hgvsp": "ENSP00000501943.1:p.Leu112CysfsTer10"}]
}"#;

fn documented(masks: &[&str]) -> Documented {
    let masks: BTreeSet<String> = masks.iter().map(|class| class.to_string()).collect();
    Documented::from_manifest(None, &masks, None, None, None).unwrap()
}

fn json(text: &str) -> Value {
    pyjson::loads(text).unwrap_or_else(|error| panic!("{text}: {error}"))
}

fn compact(value: &Value) -> String {
    pyjson::dumps_compact(value)
}

/// One step into a value: a member of an object or an item of a list.
#[derive(Clone, Copy)]
enum Step<'a> {
    Key(&'a str),
    Item(usize),
}

/// `value` with the value at `path` replaced by `edit` of it; a member `edit` returns
/// nothing for is removed, and a member the object lacks is added.
fn edited(value: Value, path: &[Step], edit: &dyn Fn(Option<Value>) -> Option<Value>) -> Value {
    let Some((step, rest)) = path.split_first() else {
        return edit(Some(value)).expect("the whole value is kept");
    };
    match *step {
        Step::Key(name) => {
            let mut members = value.into_dict().expect("an object on the path");
            if rest.is_empty() {
                match edit(members.get(name).cloned()) {
                    Some(member) => {
                        members.insert(name.to_string(), member);
                    }
                    None => {
                        members.shift_remove(name);
                    }
                }
            } else {
                let member = members.get(name).cloned().expect("a member on the path");
                members.insert(name.to_string(), edited(member, rest, edit));
            }
            Value::Dict(members)
        }
        Step::Item(at) => {
            let mut items = value.into_list().expect("a list on the path");
            let item = std::mem::replace(&mut items[at], Value::Null);
            items[at] = edited(item, rest, edit);
            Value::List(items)
        }
    }
}

/// The text `text` with the list at `path` gaining `item`.
fn appended(text: &str, path: &[Step], item: &str) -> String {
    let value = edited(json(text), path, &|list| {
        let mut items = list.expect("a list").into_list().expect("a list");
        items.push(Value::Str(item.to_string()));
        Some(Value::List(items))
    });
    compact(&value)
}

/// The text `text` with the value at `path` replaced by the JSON `replacement`.
fn replaced(text: &str, path: &[Step], replacement: &str) -> String {
    compact(&edited(json(text), path, &|_| Some(json(replacement))))
}

/// The path of a key of the one allele of the first object.
fn of_t(key: &str) -> [Step<'_>; 3] {
    [Step::Item(0), Step::Key("T"), Step::Key(key)]
}

/// The names of the fields two entries hold with different values or on one side alone.
fn differing_fields(expected: &Entry, actual: &Entry) -> Vec<String> {
    let names: BTreeSet<&String> = expected.fields.keys().chain(actual.fields.keys()).collect();
    names
        .into_iter()
        .filter(|name| expected.fields.get(*name) != actual.fields.get(*name))
        .cloned()
        .collect()
}

fn keys(entries: &[Entry]) -> Vec<&str> {
    entries.iter().map(|entry| entry.key.as_str()).collect()
}

fn refusal(text: &str) -> String {
    match parse_recoder(text, &mut documented(BUILT_IN)) {
        Err(Error::Comparator(text)) => text,
        Ok(entries) => panic!("{text}: {} entries", entries.len()),
    }
}

#[test]
fn recoder_value_canonical_forms() {
    let mut documented = documented(BUILT_IN);
    let mut value = |key: &str, text: &str| recoder_value(key, &json(text), &mut documented);
    assert_eq!(
        value("hgvsc", r#"["b", "a", "b", "NM_1.1:c.1A>G"]"#).as_deref(),
        Some(r#"["a","b"]"#)
    );
    assert_eq!(
        value("var_synonyms", r#"["ClinVar: VCV2,RCV1"]"#).as_deref(),
        Some(r#"["ClinVar: RCV1,VCV2"]"#)
    );
    assert_eq!(value("input", r#""rs1""#).as_deref(), Some(r#""rs1""#));
    assert_eq!(value("strand", "[1, -1]").as_deref(), Some("[-1,1]"));
    assert_eq!(
        documented.applied.iter().collect::<Vec<_>>(),
        [("refseq_transcript_entries", 1)]
    );
}

#[test]
fn a_list_every_item_of_which_a_class_drops_is_no_value_and_an_empty_list_is_one() {
    let mut documented = documented(&[BUILT_IN[0], BUILT_IN[1], ARCHIVE]);
    let mut value = |key: &str, text: &str| recoder_value(key, &json(text), &mut documented);
    assert_eq!(value("hgvsp", r#"["NP_1:p.A1T", "LRG_1p1:p.A1T"]"#), None);
    assert_eq!(value("var_synonyms", r#"["Archive dbSNP: rs2,rs1"]"#), None);
    assert_eq!(value("hgvsp", "[]").as_deref(), Some("[]"));
    // A prefix is read under the keys that name a transcript, a protein or a region alone.
    assert_eq!(
        value("id", r#"["NM_1", "LRG_1"]"#).as_deref(),
        Some(r#"["LRG_1","NM_1"]"#)
    );
    // A triple is dropped at its first masked notation and counted once.
    assert_eq!(
        value(
            "mane_select",
            r#"[{"hgvsc": "LRG_1t1:c.1A>G", "hgvsp": "NP_1:p.A1T"}, {"hgvsp": "x", "hgvsc": "y"}]"#
        )
        .as_deref(),
        Some(r#"[{"hgvsc":"y","hgvsp":"x"}]"#)
    );
    assert_eq!(
        documented.applied.iter().collect::<Vec<_>>(),
        [
            ("refseq_transcript_entries", 1),
            ("lrg_transcript_entries", 2),
            ("archive_dbsnp_synonyms", 1)
        ]
    );
}

#[test]
fn recoder_format_compares_each_alleles_keys_as_sets_under_the_masks() {
    let mut documented = documented(BUILT_IN);
    let expected = parse_recoder(EXPECTED, &mut documented).unwrap();
    assert_eq!(
        documented.applied.iter().collect::<Vec<_>>(),
        [
            ("refseq_transcript_entries", 3),
            ("lrg_transcript_entries", 1)
        ]
    );
    let actual = parse_recoder(ACTUAL, &mut documented).unwrap();
    assert_eq!(expected, actual);
    assert_eq!(
        keys(&expected),
        ["rs142513484|T", "null", "warnings|warnings"]
    );
    assert_eq!(expected[0].allele, "T");
    assert_eq!(expected[0].fields.len(), 10);
    assert_eq!(expected[0].fields["input"], r#""rs142513484""#);
    assert_eq!(
        expected[0].fields["var_synonyms"],
        r#"["ClinVar: RCV003317418,VCV000829541","OMIM: 606870.0005"]"#
    );
    assert_eq!(
        (expected[1].allele.as_str(), expected[1].fields.len()),
        ("", 0)
    );
    assert_eq!(
        expected[2].fields.iter().collect::<Vec<_>>(),
        [(
            &"warnings".to_string(),
            &r#"["No variant found with ID 'rs_nope'"]"#.to_string()
        )]
    );
    assert_eq!(
        (
            documented.applied.get(BUILT_IN[0]),
            documented.applied.get(BUILT_IN[1])
        ),
        (3, 1)
    );
}

#[test]
fn recoder_format_names_a_differing_key_and_a_missing_object() {
    let expected = parse_recoder(EXPECTED, &mut documented(BUILT_IN)).unwrap();
    let changed = replaced(ACTUAL, &of_t("hgvsg"), r#"["NC_000021.9:g.25585733C>A"]"#);
    let changed = appended(&changed, &of_t("hgvsc"), "ENST00000999999.1:c.1A>G");
    let actual = parse_recoder(&changed, &mut documented(BUILT_IN)).unwrap();
    assert_eq!(
        differing_fields(&expected[0], &actual[0]),
        ["hgvsc", "hgvsg"]
    );
    assert_eq!(expected[1..], actual[1..]);

    // An object the candidate lacks is two entries fewer here: the null and the warnings.
    let first = compact(&Value::List(vec![json(ACTUAL)
        .into_list()
        .unwrap()
        .remove(0)]));
    let alone = parse_recoder(&first, &mut documented(BUILT_IN)).unwrap();
    assert_eq!(keys(&alone), ["rs142513484|T"]);

    // The objects come in the file's order.
    let mut objects = json(EXPECTED).into_list().unwrap();
    objects.swap(0, 1);
    let swapped = parse_recoder(&compact(&Value::List(objects)), &mut documented(BUILT_IN));
    assert_eq!(
        keys(&swapped.unwrap()),
        ["null", "rs142513484|T", "warnings|warnings"]
    );

    // The alleles of one object come in the order of their names, whatever the file's.
    let allele = compact(
        &json(EXPECTED).into_list().unwrap()[0]
            .get("T")
            .unwrap()
            .clone(),
    );
    let in_order = format!(r#"[{{"G": {allele}, "T": {allele}}}]"#);
    let reversed = format!(r#"[{{"T": {allele}, "G": {allele}}}]"#);
    let in_order = parse_recoder(&in_order, &mut documented(BUILT_IN)).unwrap();
    assert_eq!(keys(&in_order), ["rs142513484|G", "rs142513484|T"]);
    assert_eq!(
        in_order,
        parse_recoder(&reversed, &mut documented(BUILT_IN)).unwrap()
    );
}

#[test]
fn recoder_format_lrg_class_drops_the_genomic_entries_of_the_lrg_region() {
    let with_region = appended(EXPECTED, &of_t("hgvsg"), "LRG_1:g.12345C>T");
    let with_region = appended(&with_region, &of_t("spdi"), "LRG_1:12344:C:T");
    let with_region = appended(&with_region, &of_t("vcf_string"), "LRG_1-12345-C-T");
    let mut both = documented(BUILT_IN);
    let expected = parse_recoder(&with_region, &mut both).unwrap();
    let actual = parse_recoder(ACTUAL, &mut both).unwrap();
    assert_eq!(expected, actual);
    assert_eq!(both.applied.get("lrg_transcript_entries"), 4);

    let mut refseq_alone = documented(&BUILT_IN[..1]);
    let expected = parse_recoder(&with_region, &mut refseq_alone).unwrap();
    let actual = parse_recoder(ACTUAL, &mut refseq_alone).unwrap();
    assert_eq!(
        differing_fields(&expected[0], &actual[0]),
        ["hgvsg", "hgvsp", "spdi", "vcf_string"]
    );
    assert!(!refseq_alone.applied.contains("lrg_transcript_entries"));
}

#[test]
fn recoder_format_without_the_masks_sees_the_refseq_entries() {
    let mut none = documented(&[]);
    let expected = parse_recoder(EXPECTED, &mut none).unwrap();
    let actual = parse_recoder(ACTUAL, &mut none).unwrap();
    assert_eq!(
        differing_fields(&expected[0], &actual[0]),
        ["hgvsc", "hgvsp", "mane_select"]
    );
    assert!(none.applied.is_empty());
}

#[test]
fn recoder_format_archive_dbsnp_synonyms_class_drops_the_merged_identifiers() {
    let merged = "Archive dbSNP: rs60296503,rs17181160";
    let with_merged = appended(EXPECTED, &of_t("var_synonyms"), merged);
    let mut with_class = documented(&[BUILT_IN[0], BUILT_IN[1], ARCHIVE]);
    let expected = parse_recoder(&with_merged, &mut with_class).unwrap();
    let actual = parse_recoder(ACTUAL, &mut with_class).unwrap();
    assert_eq!(expected, actual);
    assert_eq!(with_class.applied.get(ARCHIVE), 1);

    let mut without = documented(BUILT_IN);
    let expected = parse_recoder(&with_merged, &mut without).unwrap();
    let actual = parse_recoder(ACTUAL, &mut without).unwrap();
    assert_eq!(differing_fields(&expected[0], &actual[0]), ["var_synonyms"]);

    // A list the class empties is the absent key the candidate writes.
    let merged_alone = replaced(EXPECTED, &of_t("var_synonyms"), &format!(r#"["{merged}"]"#));
    let no_key = compact(&edited(json(ACTUAL), &of_t("var_synonyms"), &|_| None));
    let mut with_class = documented(&[BUILT_IN[0], BUILT_IN[1], ARCHIVE]);
    let expected = parse_recoder(&merged_alone, &mut with_class).unwrap();
    let actual = parse_recoder(&no_key, &mut with_class).unwrap();
    assert!(!expected[0].fields.contains_key("var_synonyms"));
    assert_eq!(expected, actual);
}

#[test]
fn recoder_format_identifier_indel_class_drops_the_transcript_notations() {
    let shifted = replaced(
        DELETION,
        &[Step::Key("hgvsc")],
        r#"["ENST00000300305.7:c.334del"]"#,
    );
    let shifted = replaced(
        &shifted,
        &[Step::Key("mane_select"), Step::Item(0), Step::Key("hgvsc")],
        r#""ENST00000675419.1:c.334del""#,
    );
    let expected_text = format!(r#"[{{"-": {DELETION}}}]"#);
    let actual_text = format!(r#"[{{"-": {shifted}}}]"#);
    let masks = [BUILT_IN[0], BUILT_IN[1], IDENTIFIER_INDEL];

    let mut with_class = documented(&masks);
    let expected = parse_recoder(&expected_text, &mut with_class).unwrap();
    let actual = parse_recoder(&actual_text, &mut with_class).unwrap();
    assert_eq!(expected, actual);
    assert_eq!(keys(&expected), ["rs2057997486|-"]);
    assert_eq!(
        expected[0].fields.keys().collect::<Vec<_>>(),
        ["input", "hgvsg"]
    );
    assert_eq!(with_class.applied.get(IDENTIFIER_INDEL), 6);

    let mut without = documented(BUILT_IN);
    let expected = parse_recoder(&expected_text, &mut without).unwrap();
    let actual = parse_recoder(&actual_text, &mut without).unwrap();
    assert_eq!(
        differing_fields(&expected[0], &actual[0]),
        ["hgvsc", "mane_select"]
    );

    // The same alleles under an HGVS input are compared whole.
    let notation = r#""NC_000021.9:g.34886862del""#;
    let expected_text = replaced(
        &expected_text,
        &[Step::Item(0), Step::Key("-"), Step::Key("input")],
        notation,
    );
    let actual_text = replaced(
        &actual_text,
        &[Step::Item(0), Step::Key("-"), Step::Key("input")],
        notation,
    );
    let mut with_class = documented(&masks);
    let expected = parse_recoder(&expected_text, &mut with_class).unwrap();
    let actual = parse_recoder(&actual_text, &mut with_class).unwrap();
    assert_eq!(
        differing_fields(&expected[0], &actual[0]),
        ["hgvsc", "mane_select"]
    );
    assert!(!with_class.applied.contains(IDENTIFIER_INDEL));

    // An identifier's substitution is compared whole.
    let substitution = replaced(ACTUAL, &of_t("hgvsc"), r#"["ENST00000352957.9:c.991G>T"]"#);
    let mut with_class = documented(&masks);
    let expected = parse_recoder(EXPECTED, &mut with_class).unwrap();
    let actual = parse_recoder(&substitution, &mut with_class).unwrap();
    assert_eq!(differing_fields(&expected[0], &actual[0]), ["hgvsc"]);
    assert!(!with_class.applied.contains(IDENTIFIER_INDEL));
}

#[test]
fn an_identifier_indel_without_the_three_keys_is_counted_as_none_dropped() {
    let mut with_class = documented(&[IDENTIFIER_INDEL]);
    let entries = parse_recoder(
        r#"[{"-": {"input": "rs1", "id": ["rs1"]}}]"#,
        &mut with_class,
    );
    assert_eq!(entries.unwrap()[0].fields.len(), 2);
    assert_eq!(
        with_class.applied.iter().collect::<Vec<_>>(),
        [(IDENTIFIER_INDEL, 0)]
    );
}

#[test]
fn an_allele_is_an_identifier_indel_by_its_dash_its_notation_or_its_unequal_sequences() {
    let indel = |label: &str, allele: &str, fields: &str| {
        recoder_identifier_indel(label, allele, &json(fields))
    };
    assert!(indel("rs1", "-", "{}"));
    assert!(indel("rs1", "T", r#"{"hgvsg": ["NC_000021.9:g.100dup"]}"#));
    assert!(indel("rs1", "T", r#"{"spdi": ["NC_000021.9:99:A:AT"]}"#));
    assert!(indel("rs1", "T", r#"{"vcf_string": ["21-99-A-AT"]}"#));
    assert!(!indel(
        "rs1",
        "T",
        r#"{"hgvsg": ["NC_000021.9:g.100A>T"], "spdi": ["NC_000021.9:99:A:T"], "vcf_string": ["21-100-A-T"]}"#
    ));
    assert!(!indel(
        "rs1",
        "T",
        r#"{"hgvsg": "NC_000021.9:g.100dup", "spdi": 5, "vcf_string": "21-99-A-AT"}"#
    ));
    assert!(!indel(
        "rs1",
        "T",
        r#"{"spdi": [null, 7, "NC_000021.9:99:A"], "vcf_string": [["21-99-A-AT"]]}"#
    ));
    assert!(!indel("NC_000021.9:g.100dup", "-", "{}"));
    assert!(!indel("21\t100\trs1\tA\tAT\t.\t.\t.", "-", "{}"));
    // The white space of any script is a space, and a sequence is as long as its characters.
    assert!(!indel("rs\u{2003}1", "-", "{}"));
    assert!(indel("rs\u{200b}1", "-", "{}"));
    assert!(!indel("rs1", "T", r#"{"spdi": ["NC_1:99:é:T"]}"#));
}

#[test]
fn recoder_format_compares_ga4gh_vrs_start_as_written() {
    let start = [
        Step::Item(0),
        Step::Key("T"),
        Step::Key("ga4gh_vrs"),
        Step::Item(0),
        Step::Key("location"),
        Step::Key("interval"),
        Step::Key("start"),
    ];
    let numbered = replaced(ACTUAL, &start, "25585732");
    let expected = parse_recoder(EXPECTED, &mut documented(BUILT_IN)).unwrap();
    let actual = parse_recoder(&numbered, &mut documented(BUILT_IN)).unwrap();
    assert_eq!(differing_fields(&expected[0], &actual[0]), ["ga4gh_vrs"]);
}

#[test]
fn recoder_format_refuses_what_is_not_an_array_of_objects() {
    assert_eq!(
        refusal(r#"{"T": {}}"#),
        "the recoder output is not a JSON array"
    );
    assert_eq!(
        refusal("[1]"),
        "recoder object 1 is neither an object nor null"
    );
    assert_eq!(
        refusal("[null, [] ]"),
        "recoder object 2 is neither an object nor null"
    );
    assert_eq!(
        refusal(""),
        "the recoder output does not parse: Expecting value: line 1 column 1 (char 0)"
    );
    assert_eq!(
        refusal(r#"[{"it's": 5}]"#),
        r#"recoder object 1: allele "it's" is not an object"#
    );
    assert_eq!(
        refusal(r#"[null, {"T": {"input": "b"}, "G": {"input": "a'"}, "C": {"input": 5}}]"#),
        r#"recoder object 2 names more than one input: ["a'", 'b']"#
    );
}

#[test]
fn entries_of_one_key_are_numbered_from_the_second() {
    let text =
        r#"[{"T": {"input": "rs1"}}, null, {"T": {"input": "rs1"}}, null, {}, {"warnings": 1}]"#;
    let entries = parse_recoder(text, &mut documented(BUILT_IN)).unwrap();
    assert_eq!(
        keys(&entries),
        ["rs1|T", "null", "rs1|T#2", "null#2", "warnings|warnings"]
    );
}

/// A recoder output of one allele whose key `x` holds `value`.
fn holding(value: &str) -> String {
    format!(r#"[{{"T": {{"x": {value}}}}}]"#)
}

/// What ends the run on `text`, an input no reading of the comparison expects; nothing
/// where the text parses.
fn uncaught(text: &str) -> Option<String> {
    let parsed = std::panic::catch_unwind(|| parse_recoder(text, &mut documented(BUILT_IN)));
    match parsed {
        Ok(entries) => {
            entries.expect("the text parses");
            None
        }
        Err(payload) => Some(
            payload
                .downcast::<Uncaught>()
                .expect("an uncaught exception")
                .0,
        ),
    }
}

#[test]
fn a_value_nested_past_the_interpreters_frames_ends_the_run() {
    const RECURSION: &str = "maximum recursion depth exceeded";
    let lists = |depth: usize, inner: &str| "[".repeat(depth) + inner + &"]".repeat(depth);
    let objects = |depth: usize, inner: &str| r#"{"a":"#.repeat(depth) + inner + &"}".repeat(depth);
    // The deepest value read and the first refused; a list's items are written one by
    // one, an object whole.
    let cases = [
        (lists(497, "1"), lists(498, "1")),
        (lists(496, r#""s""#), lists(497, r#""s""#)),
        (lists(495, "null"), lists(496, "null")),
        (lists(497, ""), lists(498, "")),
        (objects(495, "1"), objects(496, "1")),
        (objects(495, "{}"), objects(496, "{}")),
        (
            format!("[{}]", objects(495, "1.5")),
            format!("[{}]", objects(496, "1.5")),
        ),
    ];
    for (read, refused) in cases {
        assert_eq!(uncaught(&holding(&read)), None);
        assert_eq!(uncaught(&holding(&refused)).as_deref(), Some(RECURSION));
    }
    let warnings = |depth: usize| format!(r#"[{{"warnings": {}}}]"#, lists(depth, "1"));
    assert_eq!(uncaught(&warnings(497)), None);
    assert_eq!(uncaught(&warnings(498)).as_deref(), Some(RECURSION));
}

fn unhex(hex: &str) -> String {
    let bytes: Vec<u8> = (0..hex.len() / 2)
        .map(|at| u8::from_str_radix(&hex[2 * at..2 * at + 2], 16).expect("two hex digits"))
        .collect();
    String::from_utf8(bytes).expect("the hex of UTF-8 bytes")
}

fn hex(text: &str) -> String {
    let mut out = String::with_capacity(2 * text.len());
    for byte in text.bytes() {
        let _ = write!(out, "{byte:02x}");
    }
    out
}

fn text_value(text: &str) -> Value {
    Value::Str(text.to_string())
}

/// What `parse_recoder` makes of a text under a mask set, spelt as the case files spell it.
fn answer(masks: &str, text: &str) -> String {
    let masks: Vec<&str> = masks.split(',').filter(|class| !class.is_empty()).collect();
    let mut documented = documented(&masks);
    let entries = match parse_recoder(text, &mut documented) {
        Ok(entries) => entries,
        Err(Error::Comparator(text)) => return format!("E{}", hex(&text)),
    };
    let shown = entries
        .iter()
        .map(|entry| {
            let bare = Entry {
                key: entry.key.clone(),
                allele: entry.allele.clone(),
                fields: entry.fields.clone(),
                ..Entry::default()
            };
            assert_eq!(
                *entry, bare,
                "an entry holds its key, allele and fields alone"
            );
            let fields = entry
                .fields
                .iter()
                .map(|(name, value)| Value::List(vec![text_value(name), text_value(value)]))
                .collect();
            Value::List(vec![
                text_value(&entry.key),
                text_value(&entry.allele),
                Value::List(fields),
            ])
        })
        .collect();
    let applied = documented
        .applied
        .iter()
        .map(|(class, count)| Value::List(vec![text_value(class), Value::Int(count)]))
        .collect();
    let both = Value::List(vec![Value::List(shown), Value::List(applied)]);
    format!("K{}", hex(&pyjson::dumps(&both, &DumpOptions::default())))
}

/// The count of the cases of `file`, every one of which `parse_recoder` answers as written.
fn check(file: &Path) -> usize {
    let cases =
        std::fs::read_to_string(file).unwrap_or_else(|error| panic!("{}: {error}", file.display()));
    let mut count = 0;
    for (at, line) in cases.lines().enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        let [masks, text, expected] = fields[..] else {
            panic!("{}:{}: not three fields", file.display(), at + 1);
        };
        let (masks, text) = (unhex(masks), unhex(text));
        let found = std::panic::catch_unwind(|| answer(&masks, &text))
            .unwrap_or_else(|_| panic!("{}:{}: the parser ends the run", file.display(), at + 1));
        assert!(
            found == expected,
            "{}:{}: masks {masks:?}\n found {}\n wants {}",
            file.display(),
            at + 1,
            unhex(&found[1..]),
            unhex(&expected[1..])
        );
        count += 1;
    }
    count
}

#[test]
fn the_parser_answers_the_committed_cases_as_they_are_written() {
    let file = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/recoder/cases.tsv");
    assert_eq!(check(&file), 159);
}

#[test]
#[ignore = "reads the cases under VEP_FORMATS_COMPARE_PYDIFF_DIR"]
fn the_parser_answers_every_generated_case_as_it_is_written() {
    let dir = std::env::var_os("VEP_FORMATS_COMPARE_PYDIFF_DIR")
        .expect("VEP_FORMATS_COMPARE_PYDIFF_DIR names the directory of the cases");
    let count = check(&PathBuf::from(dir).join("recoder/cases.tsv"));
    println!("{count} cases agree");
}
