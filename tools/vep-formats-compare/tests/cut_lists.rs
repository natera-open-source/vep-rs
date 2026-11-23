// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of a known variant's list the cache cut at its cap, on planted records in the
//! default, tab, VCF and JSON formats.

mod compare_common;
mod support;

use compare_common::at;
use compare_common::compare;
use compare_common::counts;
use compare_common::default_row;
use compare_common::default_text;
use compare_common::documented_with;
use compare_common::number;
use compare_common::texts;
use compare_common::verdict;
use indexmap::IndexMap;
use vep_formats_compare::cli::parse_masks;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;

const CUT_LIST: &str = "colocated_list_cut_at_cache_cap";
const CUT_LIST_INPUT: &str = "21\t500\trs9990000001\tA\tG\t.\t.\t.";
const FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];
const ROW_FORMATS: [&str; 3] = ["default", "tab", "vcf"];

/// The `pubmed` and the `var_synonyms` value a known variant stores, `None` for none.
type Known<'a> = (Option<&'a str>, Option<&'a str>);

/// Per format, the verdict under the class and the entries counted under it.
type Outcome = Vec<(&'static str, String, i64)>;

fn inputs() -> Vec<String> {
    vec![CUT_LIST_INPUT.to_string()]
}

/// What a record of two rows comes to under the class in each of `formats`: set aside, the
/// rows counted in the default, tab and VCF formats and the record's object in JSON.
fn set_aside(formats: &[&'static str]) -> Outcome {
    formats
        .iter()
        .map(|fmt| {
            (
                *fmt,
                "green".to_string(),
                if *fmt == "json" { 1 } else { 2 },
            )
        })
        .collect()
}

/// Refused, nothing counted.
fn refused(formats: &[&'static str]) -> Outcome {
    formats
        .iter()
        .map(|fmt| (*fmt, "red".to_string(), 0))
        .collect()
}

/// Nothing to set aside, the two outputs agreeing.
fn agree(formats: &[&'static str]) -> Outcome {
    formats
        .iter()
        .map(|fmt| (*fmt, "green".to_string(), 0))
        .collect()
}

/// The outcome of `rows` in the default, tab and VCF formats and of `json` in JSON.
fn rows_and_json(rows: (&str, i64), json: (&str, i64)) -> Outcome {
    FORMATS
        .iter()
        .map(|fmt| {
            let (verdict, counted) = if *fmt == "json" { json } else { rows };
            (*fmt, verdict.to_string(), counted)
        })
        .collect()
}

/// 113 invented eight-digit PubMed ids.
fn pubmed_whole() -> String {
    (1..114)
        .map(|n| (99_000_000 + n).to_string())
        .collect::<Vec<_>>()
        .join(",")
}

/// `pubmed_whole` as the cache cuts it at 1,024 bytes, inside a 114th id.
fn pubmed_cut() -> String {
    format!("{},9900011", pubmed_whole())
}

/// An invented `var_synonyms` value of exactly 1,024 bytes that ends in `tail`, and the
/// value before `tail`: `lead`, a UniProt group whose one id is as long as the length
/// needs, then a ClinVar group of 70 ids.
fn synonyms_ending(tail: &str, lead: &str) -> (String, String) {
    let ids: Vec<String> = (1..71)
        .map(|n| {
            format!(
                "{}{}",
                if n % 2 == 1 { "RCV" } else { "VCV" },
                999_000_000 + n
            )
        })
        .collect();
    let clinvar = format!("--ClinVar::{}", ids.join(","));
    let used = format!("{lead}UniProt::VAR_{clinvar}{tail}").len();
    let body = format!("{lead}UniProt::VAR_{}{clinvar}", "9".repeat(1024 - used));
    assert_eq!(body.len() + tail.len(), 1024);
    (format!("{body}{tail}"), body)
}

/// The fields of `text` without the trailing empty ones.
fn perl_split<'a>(text: &'a str, separator: &str) -> Vec<&'a str> {
    let mut fields: Vec<&str> = text.split(separator).collect();
    while fields.last() == Some(&"") {
        fields.pop();
    }
    fields
}

/// The ids of a stored list as Ensembl VEP's JSON prints them: split on `,`, an id that
/// reads as a number printed as one, digits that close on a point or on a decimal part of
/// zeros as the number of those digits (Utils.pm:241 adds 0 to the item).
fn json_ids(stored: &str) -> Value {
    let digits = |text: &str| !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit());
    Value::List(
        perl_split(stored, ",")
            .into_iter()
            .map(|item| match item.split_once('.') {
                None if digits(item) => Value::of_int(item),
                Some((whole, fraction))
                    if digits(whole) && fraction.bytes().all(|byte| byte == b'0') =>
                {
                    Value::of_int(whole)
                }
                Some((whole, fraction)) if digits(whole) && digits(fraction) => {
                    Value::Float(item.parse().expect("a number"))
                }
                _ => Value::Str(item.to_string()),
            })
            .collect(),
    )
}

fn object(members: Vec<(&str, Value)>) -> Value {
    Value::Dict(
        members
            .into_iter()
            .map(|(name, value)| (name.to_string(), value))
            .collect(),
    )
}

fn text(value: &str) -> Value {
    Value::Str(value.to_string())
}

/// The JSON object of a known variant storing the two values.
fn json_variant(n: usize, (pubmed, synonyms): Known) -> Value {
    let mut variant: IndexMap<String, Value> = IndexMap::new();
    variant.insert("id".to_string(), text(&format!("rs999000000{n}")));
    variant.insert("allele_string".to_string(), text("A/G"));
    if let Some(pubmed) = pubmed.filter(|pubmed| !pubmed.is_empty()) {
        variant.insert("pubmed".to_string(), json_ids(pubmed));
    }
    if let Some(synonyms) = synonyms.filter(|synonyms| !synonyms.is_empty()) {
        let mut listed: IndexMap<String, Value> = IndexMap::new();
        for group in perl_split(synonyms, "--") {
            let group = perl_split(group, "::");
            listed.insert(
                group[0].to_string(),
                json_ids(group.get(1).copied().unwrap_or("")),
            );
        }
        variant.insert("var_synonyms".to_string(), Value::Dict(listed));
    }
    Value::Dict(variant)
}

/// The record object of the JSON format over `variants`.
fn json_record(variants: Vec<Value>) -> Value {
    let consequence = object(vec![
        ("variant_allele", text("G")),
        ("transcript_id", text("ENST1")),
        (
            "consequence_terms",
            Value::List(vec![text("intron_variant")]),
        ),
    ]);
    object(vec![
        ("input", text(CUT_LIST_INPUT)),
        ("id", text("rs9990000001")),
        ("most_severe_consequence", text("intron_variant")),
        ("colocated_variants", Value::List(variants)),
        ("transcript_consequences", Value::List(vec![consequence])),
    ])
}

fn dumped(value: &Value) -> String {
    pyjson::dumps(value, &pyjson::DumpOptions::default())
}

/// One record's output in `fmt` as an engine writes it for known variants storing the
/// given `pubmed` and `var_synonyms` values: two rows in the default, tab and VCF formats,
/// `PUBMED` joining the variants' ids with `,` and `VAR_SYNONYMS` their values with `--`
/// (`;` as `%3B` in the default and VCF formats, `,` as `&` in VCF), and in JSON one
/// `colocated_variants` entry per variant.
fn colocated_output(fmt: &str, known: &[Known]) -> String {
    if fmt == "json" {
        let variants = known
            .iter()
            .enumerate()
            .map(|(n, known)| json_variant(n + 1, *known))
            .collect();
        return dumped(&json_record(variants)) + "\n";
    }
    let ids: Vec<&str> = known
        .iter()
        .filter_map(|(pubmed, _)| pubmed.filter(|pubmed| !pubmed.is_empty()))
        .flat_map(|pubmed| perl_split(pubmed, ","))
        .collect();
    let synonyms: Vec<&str> = known
        .iter()
        .filter_map(|(_, synonyms)| synonyms.filter(|synonyms| !synonyms.is_empty()))
        .collect();
    let mut lists = [ids.join(","), synonyms.join("--")];
    if fmt != "tab" {
        lists = lists.map(|value| value.replace(';', "%3B"));
    }
    match fmt {
        "default" => {
            let mut extra = "IMPACT=MODIFIER;STRAND=1".to_string();
            if !lists[0].is_empty() {
                extra += &format!(";PUBMED={}", lists[0]);
            }
            if !lists[1].is_empty() {
                extra += &format!(";VAR_SYNONYMS={}", lists[1]);
            }
            let rows: Vec<String> = (1..=2)
                .map(|n| {
                    default_row(
                        "rs9990000001",
                        "21:500",
                        "G",
                        &format!("ENST{n}"),
                        "intron_variant",
                        &extra,
                    )
                })
                .collect();
            default_text(&rows)
        }
        "tab" => {
            let header = "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tPUBMED\tVAR_SYNONYMS";
            let cells = lists
                .iter()
                .map(|value| if value.is_empty() { "-" } else { value })
                .collect::<Vec<_>>()
                .join("\t");
            let mut lines = vec![header.to_string()];
            lines.extend((1..=2).map(|n| {
                format!("rs9990000001\t21:500\tG\tENST{n}\tTranscript\tintron_variant\t{cells}")
            }));
            lines.join("\n") + "\n"
        }
        _ => {
            let fields = "Allele|Consequence|Feature_type|Feature|PUBMED|VAR_SYNONYMS";
            let header =
                format!("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: {fields}\">");
            let cells = lists.join("|").replace(',', "&");
            let csq = (1..=2)
                .map(|n| format!("G|intron_variant|Transcript|ENST{n}|{cells}"))
                .collect::<Vec<_>>()
                .join(",");
            let line = &CUT_LIST_INPUT[..CUT_LIST_INPUT.len() - 1];
            format!("{header}\n{line}CSQ={csq}\n")
        }
    }
}

/// The entries of an output of `fmt` as the format's parser reads them.
fn parsed(fmt: &str, text: &str) -> Vec<Entry> {
    match fmt {
        "default" | "tab" => parse_default(text),
        "vcf" => parse_vcf(text, &inputs()),
        _ => parse_json(text, &inputs()),
    }
    .expect("the output parses")
}

fn colocated_entries(fmt: &str, known: &[Known]) -> Vec<Entry> {
    parsed(fmt, &colocated_output(fmt, known))
}

/// The report and the `Documented` of one comparison in `fmt` under `classes`.
fn compare_colocated(
    fmt: &str,
    reference: &[Known],
    candidate: &[Known],
    classes: &[&str],
) -> (Value, Documented) {
    let mut documented = documented_with(classes);
    let report = compare(
        &colocated_entries(fmt, reference),
        &colocated_entries(fmt, candidate),
        &mut documented,
    );
    (report, documented)
}

fn outcome_in(formats: &[&'static str], reference: &[Known], candidate: &[Known]) -> Outcome {
    formats
        .iter()
        .map(|fmt| {
            let (report, documented) = compare_colocated(fmt, reference, candidate, &[CUT_LIST]);
            (
                *fmt,
                verdict(&report).to_string(),
                documented.applied.get(CUT_LIST),
            )
        })
        .collect()
}

fn cut_list_outcome(reference: &[Known], candidate: &[Known]) -> Outcome {
    outcome_in(&FORMATS, reference, candidate)
}

fn pubmed(value: &str) -> Known<'_> {
    (Some(value), None)
}

fn synonyms(value: &str) -> Known<'_> {
    (None, Some(value))
}

const NONE: Known = (None, None);

/// The counts under the names of `path`, as compact JSON.
fn counted(report: &Value, path: &[&str]) -> String {
    pyjson::dumps_compact(at(report, path))
}

#[test]
fn a_pubmed_list_cut_inside_its_last_id_is_set_aside() {
    let (cut, whole) = (pubmed_cut(), pubmed_whole());
    assert_eq!(cut.len(), 1024);
    let (reference, candidate) = ([pubmed(&cut)], [pubmed(&whole)]);
    assert_eq!(
        cut_list_outcome(&reference, &candidate),
        set_aside(&FORMATS)
    );
    for fmt in FORMATS {
        let (report, _) = compare_colocated(fmt, &reference, &candidate, &[]);
        assert_eq!(verdict(&report), "red", "{fmt}");
    }
    let (report, _) = compare_colocated("default", &reference, &candidate, &[CUT_LIST]);
    let masked = format!("{{\"{CUT_LIST}\":2}}");
    assert_eq!(counted(&report, &["columns", "PUBMED", "masked"]), masked);
    assert_eq!(counted(&report, &["columns", "Extra", "masked"]), masked);
    assert_eq!(
        texts(at(&report, &["mask_examples", CUT_LIST]))[..2],
        [
            "21:500|G|ENST1|Transcript: Extra expected 1056 characters, got 1048",
            "21:500|G|ENST1|Transcript: PUBMED expected 1024 characters, got 1016",
        ]
    );
    let (report, _) = compare_colocated("json", &reference, &candidate, &[CUT_LIST]);
    assert_eq!(
        counted(&report, &["columns", "colocated_variants", "masked"]),
        format!("{{\"{CUT_LIST}\":1}}")
    );
    // A last id as long as a whole one goes too.
    let kept = whole.replace(",99000113", ",9900113");
    let whole_length = format!("{kept},99000114");
    assert_eq!(whole_length.len(), 1024);
    assert_eq!(
        cut_list_outcome(&[pubmed(&whole_length)], &[pubmed(&kept)]),
        set_aside(&FORMATS)
    );
}

#[test]
fn a_synonyms_list_cut_inside_its_last_item_is_set_aside() {
    let kept_after = [
        (",RCV999", ""),
        ("--OMIM::9990", ""),
        ("--OMIM::9990,99", "--OMIM::9990"),
        ("--OMIM", ""),
        ("--OMI", ""),
        (",RCV999000999-", ""),
        ("--OMIM:", ""),
    ];
    for (tail, kept_tail) in kept_after {
        let (stored, body) = synonyms_ending(tail, "");
        let kept = format!("{body}{kept_tail}");
        let (reference, candidate) = ([synonyms(&stored)], [synonyms(&kept)]);
        assert_eq!(
            cut_list_outcome(&reference, &candidate),
            set_aside(&FORMATS),
            "{tail}"
        );
        for fmt in FORMATS {
            let (report, _) = compare_colocated(fmt, &reference, &candidate, &[]);
            assert_eq!(verdict(&report), "red", "{tail} {fmt}");
        }
    }
    // The cap counts bytes, and the `%3B` of the default and VCF formats is the one byte
    // of a `;`.
    let (stored, body) = synonyms_ending(",RCV999", "LSDB::VAR_9;99900\u{e9}--");
    assert_eq!(stored.chars().count(), 1023);
    assert_eq!(
        cut_list_outcome(&[synonyms(&stored)], &[synonyms(&body)]),
        set_aside(&FORMATS)
    );
}

#[test]
fn a_pubmed_list_off_the_cap_or_differing_earlier_is_refused() {
    let (cut, whole) = (pubmed_cut(), pubmed_whole());
    let kept = [pubmed(&whole)];
    assert_eq!(
        cut_list_outcome(&[pubmed(&cut)], &kept),
        set_aside(&FORMATS)
    );
    let (under, over) = (format!("{whole},990011"), format!("{whole},99000114"));
    assert_eq!((under.len(), over.len()), (1023, 1025));
    assert_eq!(
        cut_list_outcome(&[pubmed(&under)], &kept),
        refused(&FORMATS)
    );
    assert_eq!(cut_list_outcome(&[pubmed(&over)], &kept), refused(&FORMATS));
    let earlier = whole.replace("99000050", "99000950");
    assert_eq!(
        cut_list_outcome(&[pubmed(&cut)], &[pubmed(&earlier)]),
        refused(&FORMATS)
    );
    // A longer list whose first 1,024 bytes end inside an id holds no run of the class.
    let (longer, glued) = (format!("{cut}4,99000115"), format!("{whole}4,99000115"));
    assert_eq!(
        cut_list_outcome(&[pubmed(&longer)], &[pubmed(&glued)]),
        refused(&FORMATS)
    );
}

#[test]
fn a_synonyms_list_off_the_cap_or_differing_earlier_is_refused() {
    for tail in [",RCV999", "--OMIM::9990", "--OMIM"] {
        let (stored, body) = synonyms_ending(tail, "");
        let (reference, candidate) = ([synonyms(&stored)], [synonyms(&body)]);
        assert_eq!(
            cut_list_outcome(&reference, &candidate),
            set_aside(&FORMATS),
            "{tail}"
        );
        // One byte under the cap the value is refused, in a row always and in JSON unless its
        // last id is digits alone: JSON prints a piece that closes on a point as its digits, so
        // the object of such a value is also that of one cut at the cap after a point.
        let under = [stored.replacen('9', "", 1), body.replacen('9', "", 1)];
        assert_eq!(under[0].len(), 1023);
        let refusing: &[&'static str] = if tail == "--OMIM::9990" {
            &ROW_FORMATS
        } else {
            &FORMATS
        };
        assert_eq!(
            outcome_in(refusing, &[synonyms(&under[0])], &[synonyms(&under[1])]),
            refused(refusing),
            "{tail}"
        );
        let over = format!("{stored}0");
        assert_eq!(
            cut_list_outcome(&[synonyms(&over)], &candidate),
            refused(&FORMATS),
            "{tail}"
        );
        let earlier = body.replace("RCV999000011", "RCV999000911");
        assert_eq!(
            cut_list_outcome(&reference, &[synonyms(&earlier)]),
            refused(&FORMATS),
            "{tail}"
        );
    }
}

#[test]
fn a_missing_and_a_wrong_pubmed_value_stay_differences() {
    let (cut, whole) = (pubmed_cut(), pubmed_whole());
    let reference = [pubmed(&cut)];
    assert_eq!(
        cut_list_outcome(&reference, &[pubmed(&whole)]),
        set_aside(&FORMATS)
    );
    let missing = [NONE];
    assert_eq!(cut_list_outcome(&reference, &missing), refused(&FORMATS));
    let short = [pubmed("99000001,99000002")];
    assert_eq!(cut_list_outcome(&short, &missing), refused(&FORMATS));
    let (report, _) = compare_colocated("default", &reference, &missing, &[CUT_LIST]);
    assert_eq!(number(&report, &["columns", "PUBMED", "differing"]), 2);
    let wrong = [
        ("another fragment", format!("{whole},9900012")),
        (
            "two ids fewer",
            whole.rsplit_once(',').expect("a list").0.to_string(),
        ),
        ("an id more", format!("{whole},99000114")),
    ];
    for (planted, value) in &wrong {
        assert_eq!(
            cut_list_outcome(&reference, &[pubmed(value)]),
            refused(&FORMATS),
            "{planted}"
        );
    }
    // A JSON entry whose other member differs beside the cut list, or whose ids are
    // strings where the reference prints numbers.
    let known = json_variant(1, pubmed(&whole));
    let as_text = Value::List(
        known
            .get("pubmed")
            .and_then(Value::as_list)
            .expect("the ids")
            .iter()
            .map(|id| text(&pyjson::py_str(id)))
            .collect(),
    );
    for (member, value) in [("allele_string", text("A/T")), ("pubmed", as_text)] {
        let mut planted = known.as_dict().expect("an object").clone();
        planted.insert(member.to_string(), value);
        let record = dumped(&json_record(vec![Value::Dict(planted)]));
        let entries = parse_json(&record, &inputs()).expect("the record parses");
        let mut documented = documented_with(&[CUT_LIST]);
        let expected = colocated_entries("json", &reference);
        let report = compare(&expected, &entries, &mut documented);
        assert_eq!(verdict(&report), "red", "{member}");
        assert!(counts(&documented.applied).is_empty(), "{member}");
    }
}

#[test]
fn a_missing_and_a_wrong_synonyms_value_stay_differences() {
    let (stored, body) = synonyms_ending(",RCV999", "");
    let reference = [synonyms(&stored)];
    assert_eq!(
        cut_list_outcome(&reference, &[synonyms(&body)]),
        set_aside(&FORMATS)
    );
    let missing = [NONE];
    assert_eq!(cut_list_outcome(&reference, &missing), refused(&FORMATS));
    let short = [synonyms(
        "ClinVar::RCV999000001,VCV999000002--UniProt::VAR_999001",
    )];
    assert_eq!(cut_list_outcome(&short, &missing), refused(&FORMATS));
    let (report, _) = compare_colocated("tab", &reference, &missing, &[CUT_LIST]);
    assert_eq!(
        number(&report, &["columns", "VAR_SYNONYMS", "differing"]),
        2
    );
    let wrong = [
        ("another fragment", format!("{body},RCV998")),
        (
            "a source name the stored value lacks",
            format!("{body}--OMIM"),
        ),
        (
            "an id fewer",
            body.rsplit_once(',').expect("a list").0.to_string(),
        ),
        (
            "the value without its last group",
            body.split_once("--ClinVar").expect("a group").0.to_string(),
        ),
    ];
    for (planted, value) in &wrong {
        assert_eq!(
            cut_list_outcome(&reference, &[synonyms(value)]),
            refused(&FORMATS),
            "{planted}"
        );
    }
}

#[test]
fn a_source_name_kept_without_its_only_id_stays_a_difference() {
    let (stored, body) = synonyms_ending("--OMIM::9990", "");
    let kept_name = format!("{body}--OMIM");
    let (cut, whole, bare) = (synonyms(&stored), synonyms(&body), synonyms(&kept_name));
    let other = synonyms("UniProt::VAR_999002,VAR_999003");
    assert_eq!(cut_list_outcome(&[cut], &[whole]), set_aside(&FORMATS));
    assert_eq!(
        cut_list_outcome(&[cut, other], &[whole, other]),
        set_aside(&FORMATS)
    );
    assert_eq!(
        cut_list_outcome(&[other, cut], &[other, whole]),
        set_aside(&FORMATS)
    );
    assert_eq!(cut_list_outcome(&[cut], &[bare]), refused(&FORMATS));
    assert_eq!(
        cut_list_outcome(&[cut, other], &[bare, other]),
        refused(&FORMATS)
    );
    assert_eq!(
        cut_list_outcome(&[other, cut], &[other, bare]),
        refused(&FORMATS)
    );
    let listed: Vec<String> = [cut, whole, bare]
        .into_iter()
        .map(|known| {
            let entry = json_variant(1, known);
            let sources = entry.get("var_synonyms").expect("the sources");
            sources
                .get("OMIM")
                .map_or("None".to_string(), pyjson::dumps_compact)
        })
        .collect();
    assert_eq!(listed, ["[9990]", "None", "[]"]);
}

#[test]
fn the_joined_lists_of_two_known_variants_are_read() {
    let (cut_ids, whole_ids) = (pubmed_cut(), pubmed_whole());
    let (cut_synonyms, body) = synonyms_ending(",RCV999", "");
    let cut: Known = (Some(&cut_ids), Some(&cut_synonyms));
    let whole: Known = (Some(&whole_ids), Some(&body));
    let other: Known = (
        Some("9900201,99000202"),
        Some("UniProt::VAR_999002,VAR_999003"),
    );
    assert_eq!(
        cut_list_outcome(&[cut, other], &[whole, other]),
        set_aside(&FORMATS)
    );
    assert_eq!(
        cut_list_outcome(&[other, cut], &[other, whole]),
        set_aside(&FORMATS)
    );
    assert_eq!(
        cut_list_outcome(&[cut, cut], &[whole, whole]),
        set_aside(&FORMATS)
    );
    let output = colocated_output("tab", &[cut, other]);
    let joined: Vec<&str> = output.lines().nth(1).expect("a row").split('\t').collect();
    assert_eq!(joined[6], format!("{cut_ids},{}", other.0.unwrap()));
    assert_eq!(joined[7], format!("{cut_synonyms}--{}", other.1.unwrap()));
    // A first list one byte under the cap stays a difference.
    let under_ids = format!("{whole_ids},990011");
    let (under_synonyms, under_body) =
        (cut_synonyms.replacen('9', "", 1), body.replacen('9', "", 1));
    let under: Known = (Some(&under_ids), Some(&under_synonyms));
    let under_whole: Known = (Some(&whole_ids), Some(&under_body));
    assert_eq!(
        cut_list_outcome(&[under, other], &[under_whole, other]),
        refused(&FORMATS)
    );
    // So does a candidate that leaves out an id ending no run of 1,024 bytes.
    let shorter: Known = (Some("9900201"), Some("UniProt::VAR_999002"));
    assert_eq!(
        cut_list_outcome(&[cut, other], &[cut, shorter]),
        refused(&FORMATS)
    );
}

#[test]
fn a_run_across_two_known_variants_reads_as_one_list_in_a_row() {
    let cut = pubmed_cut();
    let reference = [pubmed(&cut), pubmed("99000201,99000202")];
    let candidate = [pubmed(&cut), pubmed("99000201")];
    let in_a_row = rows_and_json(("green", 2), ("red", 0));
    assert_eq!(cut_list_outcome(&reference, &candidate), in_a_row);
    let (joined, first) = synonyms_ending("--OMIM::9990", "");
    let reference = [synonyms(&first), synonyms("OMIM::9990")];
    let output = colocated_output("tab", &reference);
    assert!(output
        .lines()
        .nth(1)
        .expect("a row")
        .ends_with(&format!("\t{joined}")));
    assert_eq!(
        cut_list_outcome(&reference, &[synonyms(&first), NONE]),
        in_a_row
    );
}

#[test]
fn a_list_the_cache_cut_after_a_separator_holds_nothing_of_the_class() {
    let (cut, whole) = (pubmed_cut(), pubmed_whole());
    assert_eq!(
        cut_list_outcome(&[pubmed(&cut)], &[pubmed(&whole)]),
        set_aside(&FORMATS)
    );
    let ending = format!("{whole},990011,");
    assert_eq!(ending.len(), 1024);
    let ids = [pubmed(&ending)];
    assert_eq!(cut_list_outcome(&ids, &ids), agree(&FORMATS));
    assert_eq!(cut_list_outcome(&ids, &[pubmed(&whole)]), refused(&FORMATS));
    let rows_alone = rows_and_json(("red", 0), ("green", 0));
    for (tail, separator) in [(",", ","), ("--OMIM::", "::"), ("--", "--")] {
        let (stored, body) = synonyms_ending(tail, "");
        let value = [synonyms(&stored)];
        assert_eq!(cut_list_outcome(&value, &value), agree(&FORMATS), "{tail}");
        let shorter = body.rsplit_once(',').expect("a list").0;
        assert_eq!(
            cut_list_outcome(&value, &[synonyms(shorter)]),
            refused(&FORMATS),
            "{tail}"
        );
        let without_separator = stored.strip_suffix(separator).expect("the separator");
        assert_eq!(
            cut_list_outcome(&value, &[synonyms(without_separator)]),
            rows_alone,
            "{tail}"
        );
    }
}

#[test]
fn the_field_a_cut_value_with_no_id_does_not_write_is_set_aside() {
    let (long, short) = ("9".repeat(1024), "9".repeat(1023));
    let (item, other) = (pubmed(&long), pubmed("99000201"));
    let rows = set_aside(&ROW_FORMATS);
    assert_eq!(outcome_in(&ROW_FORMATS, &[item], &[NONE]), rows);
    assert_eq!(
        outcome_in(&ROW_FORMATS, &[item, other], &[NONE, other]),
        rows
    );
    assert_eq!(
        outcome_in(&ROW_FORMATS, &[other, item], &[other, NONE]),
        rows
    );
    assert_eq!(
        outcome_in(&ROW_FORMATS, &[pubmed(&short)], &[NONE]),
        refused(&ROW_FORMATS)
    );
    for stored in [
        format!("ClinVar{}", "9".repeat(1017)),
        format!("ClinVar::RCV{}", "9".repeat(1012)),
    ] {
        assert_eq!(stored.len(), 1024);
        assert_eq!(
            cut_list_outcome(&[synonyms(&stored)], &[NONE]),
            set_aside(&FORMATS)
        );
        let second = synonyms("UniProt::VAR_999002");
        assert_eq!(
            cut_list_outcome(&[synonyms(&stored), second], &[NONE, second]),
            set_aside(&FORMATS)
        );
        assert_eq!(
            cut_list_outcome(&[synonyms(&stored[..1023])], &[NONE]),
            refused(&FORMATS)
        );
    }
}

#[test]
fn the_json_shapes_that_hide_the_stored_length_stay_differences() {
    let hidden = rows_and_json(("green", 2), ("red", 0));
    for lead in [
        "OMIM--",
        "OMIM::--",
        "ClinVar::RCV999000901--",
        "OMIM::999001.0010--",
    ] {
        let (stored, body) = synonyms_ending(",RCV999", lead);
        assert_eq!(
            cut_list_outcome(&[synonyms(&stored)], &[synonyms(&body)]),
            hidden,
            "{lead}"
        );
    }
    // A piece that closes on a zero after its point prints two bytes shorter than it is
    // stored, and a whole id that closes on a zero one byte shorter: read with a point after
    // its last id the rebuilt value is still short in JSON, where the row formats read the
    // stored text.
    for (tail, lead) in [
        ("--OMIM::999001.0", ""),
        ("--OMIM::999001.", "LSDB::999003.0010--"),
    ] {
        let (stored, body) = synonyms_ending(tail, lead);
        assert_eq!(
            cut_list_outcome(&[synonyms(&stored)], &[synonyms(&body)]),
            hidden,
            "{tail} {lead}"
        );
    }
    // An id whose number prints as stored is set aside in JSON too.
    let (stored, body) = synonyms_ending(",RCV999", "OMIM::999001.0001--");
    assert_eq!(
        cut_list_outcome(&[synonyms(&stored)], &[synonyms(&body)]),
        set_aside(&FORMATS)
    );
}

#[test]
fn a_piece_cut_after_its_point_is_set_aside_where_json_prints_a_whole_number() {
    // The reference's object holds the piece without its point, a byte short of the stored
    // value: the only id of its source, and the last of two.
    let pieces = [
        ("--OMIM::999001.", "", "[999001]"),
        (
            "--OMIM::999001.0001,999002.",
            "--OMIM::999001.0001",
            "[999001.0001,999002]",
        ),
    ];
    for (tail, kept_tail, listed) in pieces {
        let (stored, body) = synonyms_ending(tail, "");
        let kept = format!("{body}{kept_tail}");
        let (reference, candidate) = ([synonyms(&stored)], [synonyms(&kept)]);
        let entry = json_variant(1, reference[0]);
        let sources = entry.get("var_synonyms").expect("the sources");
        assert_eq!(
            pyjson::dumps_compact(sources.get("OMIM").expect("the source")),
            listed,
            "{tail}"
        );
        assert_eq!(
            cut_list_outcome(&reference, &candidate),
            set_aside(&FORMATS),
            "{tail}"
        );
        for fmt in FORMATS {
            let (report, _) = compare_colocated(fmt, &reference, &candidate, &[]);
            assert_eq!(verdict(&report), "red", "{tail} {fmt}");
        }
        // A row refuses the same piece at any other stored length. The object does not hold
        // the point, so JSON reads a value one byte longer as its 1,024 bytes without it,
        // and refuses a value one byte shorter and one two bytes longer.
        let json_alone = rows_and_json(("red", 0), ("green", 1));
        let off_the_cap = [
            ("", 1023, refused(&FORMATS)),
            ("99", 1025, json_alone),
            ("999", 1026, refused(&FORMATS)),
        ];
        for (nines, length, outcome) in off_the_cap {
            let planted = [stored.replacen('9', nines, 1), kept.replacen('9', nines, 1)];
            assert_eq!(planted[0].len(), length);
            assert_eq!(
                cut_list_outcome(&[synonyms(&planted[0])], &[synonyms(&planted[1])]),
                outcome,
                "{tail} {length}"
            );
        }
        // So is a candidate that differs in an earlier id, or keeps the name of a source
        // whose only id is the piece.
        let mut wrong = vec![kept.replace("RCV999000011", "RCV999000911")];
        if kept_tail.is_empty() {
            wrong.push(format!("{kept}--OMIM"));
        }
        for value in &wrong {
            assert_eq!(
                cut_list_outcome(&reference, &[synonyms(value)]),
                refused(&FORMATS),
                "{tail}"
            );
        }
    }
    // Beside a piece of digits alone, a whole id that closes on a zero leaves the rebuilt
    // value one byte short, which the reading with the point makes 1,024 bytes, the length
    // the value is stored at.
    let (stored, body) = synonyms_ending("--OMIM::9990", "LSDB::999003.0010--");
    assert_eq!(
        cut_list_outcome(&[synonyms(&stored)], &[synonyms(&body)]),
        set_aside(&FORMATS)
    );
}

#[test]
fn a_value_of_whole_ids_json_prints_as_a_piece_cut_after_its_point_is_set_aside_in_json_alone() {
    // Two stored values whose every id is whole print the object of a value cut after the
    // point of its last id: 1,023 bytes that end in the digits of that id, and 1,024 bytes
    // that end in a `,` after them. JSON reads a candidate lacking that id as the shape,
    // where a row refuses the pair.
    let json_alone = rows_and_json(("red", 0), ("green", 1));
    let (stored, body) = synonyms_ending("--OMIM::9990", "");
    let under = [stored.replacen('9', "", 1), body.replacen('9', "", 1)];
    assert_eq!(under[0].len(), 1023);
    assert_eq!(
        cut_list_outcome(&[synonyms(&under[0])], &[synonyms(&under[1])]),
        json_alone
    );
    let (stored, body) = synonyms_ending("--OMIM::9990,", "");
    assert_eq!(
        cut_list_outcome(&[synonyms(&stored)], &[synonyms(&body)]),
        json_alone
    );
}

#[test]
fn any_group_of_a_json_object_reads_as_the_last() {
    let (stored, _) = synonyms_ending("--OMIM::9990", "LSDB::999001--");
    let groups: Vec<&str> = stored.split("--").collect();
    let [first, only, many, last] = groups[..] else {
        panic!("the value holds four groups");
    };
    let earlier = [
        ("a group of one id", vec![first, many, last]),
        (
            "the last id of a group",
            vec![first, only, many.rsplit_once(',').expect("a list").0, last],
        ),
    ];
    let json_alone = rows_and_json(("red", 0), ("green", 1));
    for (lacking, groups) in earlier {
        let candidate = groups.join("--");
        assert_eq!(
            cut_list_outcome(&[synonyms(&stored)], &[synonyms(&candidate)]),
            json_alone,
            "{lacking}"
        );
    }
}

#[test]
fn the_class_applies_only_when_the_masks_name_it() {
    let directory = support::scratch("cut-lists-masks");
    let (cut, whole) = (pubmed_cut(), pubmed_whole());
    let reference = directory.join("reference.txt");
    let candidate = directory.join("candidate.txt");
    std::fs::write(&reference, colocated_output("default", &[pubmed(&cut)])).unwrap();
    std::fs::write(&candidate, colocated_output("default", &[pubmed(&whole)])).unwrap();
    let report = directory.join("report.json");
    let run = |actual: &std::path::Path, masks: Option<&str>| {
        let mut args: Vec<&std::ffi::OsStr> = vec![
            "--format".as_ref(),
            "default".as_ref(),
            "--expected".as_ref(),
            reference.as_os_str(),
            "--report".as_ref(),
            report.as_os_str(),
            "--actual".as_ref(),
            actual.as_os_str(),
        ];
        if let Some(masks) = masks {
            args.push("--masks".as_ref());
            args.push(masks.as_ref());
        }
        support::run(&directory, &args)
    };
    let unmasked = run(&candidate, None);
    assert_eq!(unmasked.code, 1, "{}", unmasked.err());
    assert!(unmasked.out().contains("red"));
    assert!(!parse_masks(None).contains(CUT_LIST));
    assert_eq!(run(&candidate, Some(CUT_LIST)).code, 0);
    let masks = support::report(&report);
    assert_eq!(
        counted(&masks, &["masks", "applied"]),
        format!("{{\"{CUT_LIST}\":2}}")
    );
    assert_eq!(counted(&masks, &["masks", "unused"]), "[]");
    assert_eq!(run(&reference, Some(CUT_LIST)).code, 0);
    let masks = support::report(&report);
    assert_eq!(counted(&masks, &["masks", "applied"]), "{}");
    assert_eq!(
        counted(&masks, &["masks", "unused"]),
        format!("[\"{CUT_LIST}\"]")
    );
}
