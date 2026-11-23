// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison of a differing field: accepted value pairs, nearest sites, lists compared
//! as multisets and the classes that meet on one `Extra` text.

mod compare_common;

use std::collections::BTreeSet;

use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::classes::ORDER_ONLY;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::fields::extra_pairs_accepted;
use vep_formats_compare::input::input_records;
use vep_formats_compare::input::read_text;
use vep_formats_compare::input::text_lines;
use vep_formats_compare::lists::clin_sig_list_separator;
use vep_formats_compare::lists::order_only_fields;
use vep_formats_compare::lists::same_items;
use vep_formats_compare::lists::JSON_ARRAY;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::nearest::NearestSites;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::header_columns;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::Side;

type FieldPairs = IndexMap<String, BTreeSet<(String, String)>>;

fn default_entries(text: &str) -> Vec<Entry> {
    parse_default(text).expect("the text parses")
}

fn tab_entries(text: &str) -> Vec<Entry> {
    parse_tab(text).expect("the text parses")
}

fn pair_of(column: &str, reference: &str, candidate: &str) -> FieldPairs {
    let mut pairs = FieldPairs::new();
    pairs
        .entry(column.to_string())
        .or_default()
        .insert((reference.to_string(), candidate.to_string()));
    pairs
}

fn compare_with_pairs(
    expected: &[Entry],
    actual: &[Entry],
    documented: &mut Documented,
    pairs: &FieldPairs,
) -> Value {
    let side = Side::default();
    let options = CompareOptions {
        field_pairs: Some(pairs),
        ..CompareOptions::new(&side)
    };
    compare_entries(expected, actual, documented, &options).expect("the entries compare")
}

/// The report of two entry lists under the list fields `fmt` compares as multisets.
fn compare_as_lists(
    expected: &[Entry],
    actual: &[Entry],
    documented: &mut Documented,
    fmt: &str,
) -> Value {
    let side = Side::default();
    let lists = order_only_fields(fmt);
    let options = CompareOptions {
        order_only: Some(&lists),
        ..CompareOptions::new(&side)
    };
    compare_entries(expected, actual, documented, &options).expect("the entries compare")
}

/// A pair accepted on an Extra key also settles the raw Extra text that carries it; an
/// Extra whose keys are only reordered stays a difference, as every key still agrees.
#[test]
fn an_accepted_field_pair_reaches_the_extra_column() {
    let text = text_of("default");
    let expected = default_entries(&text);
    let reference_rows = data_lines(&text);
    let mut rows: Vec<String> = reference_rows
        .iter()
        .map(|r| r.replace("OverlapPC=100.00", "OverlapPC=100"))
        .collect();
    let actual = default_entries(&rebuild(&text, &rows));
    let changed = reference_rows
        .iter()
        .filter(|r| r.contains("OverlapPC=100.00"))
        .count();
    assert!(changed > 0);
    let changed = i64::try_from(changed).expect("a count");
    let pairs = pair_of("OverlapPC", "100.00", "100");
    let report = compare_with_pairs(&expected, &actual, &mut no_documented(), &pairs);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(
        number(&report, &["columns", "OverlapPC", "pairs_accepted"]),
        changed
    );
    assert_eq!(
        number(&report, &["columns", "Extra", "pairs_accepted"]),
        changed
    );
    rows[0] = rows[0].replace("IMPACT=MODIFIER", "IMPACT=LOW");
    let report = compare_with_pairs(
        &expected,
        &default_entries(&rebuild(&text, &rows)),
        &mut no_documented(),
        &pairs,
    );
    assert_eq!(number(&report, &["columns", "Extra", "differing"]), 1);
    assert_eq!(number(&report, &["columns", "IMPACT", "differing"]), 1);
    let ablation_extra = "IMPACT=HIGH;STRAND=-1;OverlapBP=23932;OverlapPC=100.00";
    let reordered: Vec<String> = reference_rows
        .iter()
        .map(|r| {
            r.replace(
                ablation_extra,
                "STRAND=-1;IMPACT=HIGH;OverlapBP=23932;OverlapPC=100.00",
            )
        })
        .collect();
    let ending = format!("\t{ablation_extra}");
    assert_eq!(
        reference_rows
            .iter()
            .filter(|r| r.ends_with(&ending))
            .count(),
        1
    );
    let report = compare(
        &expected,
        &default_entries(&rebuild(&text, &reordered)),
        &mut no_documented(),
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "Extra", "differing"]), 1);
    assert_eq!(number(&report, &["columns", "Extra", "pairs_accepted"]), 0);
    for column in ["IMPACT", "STRAND", "OverlapBP", "OverlapPC"] {
        assert_eq!(number(&report, &["columns", column, "differing"]), 0);
    }
    assert!(extra_pairs_accepted(
        "IMPACT=HIGH;OverlapPC=100.00",
        "IMPACT=HIGH;OverlapPC=100",
        &pairs
    ));
    assert!(!extra_pairs_accepted(
        "IMPACT=HIGH;OverlapPC=100.00",
        "OverlapPC=100;IMPACT=HIGH",
        &pairs
    ));
}

/// The nearest-transcript set of the identifiers corpus: every record where Ensembl VEP's
/// one transcript is one of several at the winning site differs by id and agrees by site;
/// the gene and symbol sets agree verbatim.
#[test]
fn the_nearest_sites_of_a_cache_settle_the_identifiers_corpus() {
    let nearest_corpus = golden().join("GRCh38-identifiers");
    let sites =
        NearestSites::from_cache(&nearest_corpus.join("json_cache")).expect("the cache reads");
    let site: BTreeSet<i64> = [29940022].into_iter().collect();
    assert_eq!(sites.sites.get("ENST00000389125"), Some(&site));
    assert_eq!(sites.sites.get("ENST00000399914"), Some(&site));
    for name in ["GRIK1", "ENSG00000171189"] {
        assert!(sites
            .sites
            .get(name)
            .is_some_and(|held| held.is_superset(&site)));
    }
    assert!(sites.agree(
        "ENST00000399914",
        "ENST00000389124,ENST00000389125,ENST00000399907,ENST00000399909,ENST00000399913,ENST00000399914,ENST00001131791"
    ));
    let text = read_text(
        &nearest_corpus
            .join("expected_nearest-transcript")
            .join("default.txt.gz"),
    )
    .expect("the output reads");
    let expected = default_entries(&text);
    let rows = data_lines(&text);
    let replaced = |to: &str| -> Vec<Entry> {
        let planted: Vec<String> = rows
            .iter()
            .map(|r| r.replace("NEAREST=ENST00000399913", to))
            .collect();
        default_entries(&rebuild(&text, &planted))
    };
    let changed = rows
        .iter()
        .filter(|r| r.contains("NEAREST=ENST00000399913"))
        .count();
    assert!(changed > 0);
    let changed = i64::try_from(changed).expect("a count");
    let actual = replaced("NEAREST=ENST00000389125,ENST00000399913");
    let plain = compare(&expected, &actual, &mut no_documented());
    assert_eq!(verdict(&plain), "red");
    assert_eq!(
        number(&plain, &["columns", "NEAREST", "differing"]),
        changed
    );
    let side = Side {
        nearest_sites: Some(sites),
        ..Side::default()
    };
    let with_sites = |actual: &[Entry]| -> Value {
        compare_entries(
            &expected,
            actual,
            &mut no_documented(),
            &CompareOptions::new(&side),
        )
        .expect("the entries compare")
    };
    let settled = with_sites(&actual);
    assert_eq!(verdict(&settled), "green", "{:?}", reasons(&settled));
    assert_eq!(
        number(&settled, &["columns", "NEAREST", "nearest_accepted"]),
        changed
    );
    assert_eq!(
        number(&settled, &["columns", "Extra", "nearest_accepted"]),
        changed
    );
    let report = with_sites(&replaced("NEAREST=ENST00000400546"));
    assert_eq!(verdict(&report), "red");
    assert_eq!(
        number(&report, &["columns", "NEAREST", "differing"]),
        changed
    );
    assert!(NearestSites::from_cache(&nearest_corpus.join("expected")).is_err());
}

/// Rows of the protein-features reference whose DOMAINS carries two or more items that differ,
/// so reversing the items changes the text: 118 of the 119 multi-item rows in the default
/// format, and 116 JSON consequence objects once the JSON re-split has mapped the items.
const REORDERABLE_DOMAINS_ROWS: i64 = 118;
const REORDERABLE_DOMAINS_OBJECTS: i64 = 116;

/// `text` with the items of every `column=` list of the Extra text reversed.
fn reversed_lists(text: &str, column: &str) -> String {
    let lists = regex::Regex::new(&format!(r"{column}=([^;\t\n]+)")).expect("the pattern compiles");
    lists
        .replace_all(text, |found: &regex::Captures| {
            let items: Vec<&str> = found[1].split(',').rev().collect();
            format!("{column}={}", items.join(","))
        })
        .into_owned()
}

/// The protein-features reference under `order_only`: a candidate printing the protein
/// features of a row in another order (a native cache orders them by feature id, the
/// Ensembl cache by its array) agrees as a multiset, counted per column as `order_only`;
/// without the class the same rows are differing fields, and a reorder of any other column,
/// or a repeated item, stays red under it.
#[test]
fn the_order_only_class_accepts_domains_in_any_order_in_every_format() {
    let protein_features = golden().join("GRCh38-protein-features");
    let reference_of = |name: &str| -> String {
        read_text(&protein_features.join("expected").join(name)).expect("the output reads")
    };
    let inputs = input_records(&protein_features.join("variants.vcf")).expect("the input reads");
    let default = reference_of("default.txt.gz");
    let expected = default_entries(&default);
    let actual = default_entries(&reversed_lists(&default, "DOMAINS"));
    let plain = compare(&expected, &actual, &mut no_documented());
    assert_eq!(verdict(&plain), "red");
    for column in ["DOMAINS", "Extra"] {
        assert_eq!(
            number(&plain, &["columns", column, "differing"]),
            REORDERABLE_DOMAINS_ROWS
        );
    }
    let mut documented = no_documented();
    let accepted = compare_as_lists(&expected, &actual, &mut documented, "default");
    assert_eq!(verdict(&accepted), "green");
    assert_eq!(number(&accepted, &["columns", "DOMAINS", "differing"]), 0);
    for column in ["DOMAINS", "Extra"] {
        assert_eq!(
            number(&accepted, &["columns", column, "order_only"]),
            REORDERABLE_DOMAINS_ROWS
        );
    }
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(ORDER_ONLY, 2 * REORDERABLE_DOMAINS_ROWS)])
    );
    let example = texts(at(
        &accepted,
        &["columns", "DOMAINS", "order_only_examples"],
    ))[0];
    assert!(example.starts_with(
        "21:14216193|G|ENST00000400577|Transcript: expected 'PDB-ENSP_mappings:2ywk.A,"
    ));
    assert!(example.contains("got 'AFDB-ENSP_mappings:AF-P57052-F1,PDB-ENSP_mappings:2ywk.A'"));

    // The tab format's DOMAINS column, the same rule.
    let tab = reference_of("tab.txt.gz");
    let domains_column = header_columns(&tab)
        .iter()
        .position(|name| name == "DOMAINS")
        .expect("the tab output has the column");
    let tab_rows: Vec<String> = data_lines(&tab)
        .iter()
        .map(|row| {
            let mut cells: Vec<String> = row.split('\t').map(str::to_string).collect();
            let items: Vec<&str> = cells[domains_column].split(',').rev().collect();
            cells[domains_column] = items.join(",");
            cells.join("\t")
        })
        .collect();
    let tab_report = compare_as_lists(
        &tab_entries(&tab),
        &tab_entries(&rebuild(&tab, &tab_rows)),
        &mut no_documented(),
        "tab",
    );
    assert_eq!(verdict(&tab_report), "green");
    assert_eq!(
        number(&tab_report, &["columns", "DOMAINS", "order_only"]),
        REORDERABLE_DOMAINS_ROWS
    );

    // The VCF format joins the items with & (VCF.pm:391).
    let vcf = reference_of("vcf.vcf.gz");
    let flipped: Vec<String> = data_lines(&vcf)
        .iter()
        .map(|row| {
            let (prefix, marker, csq) = match row.split_once(";CSQ=") {
                Some((prefix, csq)) => (prefix, ";CSQ=", csq),
                None => match row.split_once("CSQ=") {
                    Some((prefix, csq)) => (prefix, "CSQ=", csq),
                    None => (row.as_str(), "", ""),
                },
            };
            let chunks: Vec<String> = csq
                .split(',')
                .map(|chunk| {
                    let mut values: Vec<String> = chunk.split('|').map(str::to_string).collect();
                    if let Some(last) = values.last_mut() {
                        let items: Vec<&str> = last.split('&').rev().collect();
                        *last = items.join("&");
                    }
                    values.join("|")
                })
                .collect();
            format!("{prefix}{marker}{}", chunks.join(","))
        })
        .collect();
    let vcf_expected = parse_vcf(&vcf, &inputs).expect("the reference parses");
    let vcf_actual = parse_vcf(&rebuild(&vcf, &flipped), &inputs).expect("the candidate parses");
    let vcf_plain = compare(&vcf_expected, &vcf_actual, &mut no_documented());
    assert_eq!(
        number(&vcf_plain, &["columns", "DOMAINS", "differing"]),
        REORDERABLE_DOMAINS_ROWS
    );
    let vcf_report = compare_as_lists(&vcf_expected, &vcf_actual, &mut no_documented(), "vcf");
    assert_eq!(verdict(&vcf_report), "green");
    assert_eq!(
        number(&vcf_report, &["columns", "DOMAINS", "order_only"]),
        REORDERABLE_DOMAINS_ROWS
    );

    // The JSON format's domains array (JSON.pm:363-371), the objects compared as a multiset.
    let json_text = reference_of("json.jsonl.gz");
    let mut objects: Vec<String> = Vec::new();
    for line in text_lines(&json_text) {
        if line.is_empty() {
            continue;
        }
        let mut object = json(&line);
        if let Value::Dict(members) = &mut object {
            if let Some(Value::List(consequences)) = members.get_mut("transcript_consequences") {
                for consequence in consequences {
                    if let Value::Dict(held) = consequence {
                        if let Some(Value::List(domains)) = held.get_mut("domains") {
                            domains.reverse();
                        }
                    }
                }
            }
        }
        objects.push(pyjson::dumps(&object, &pyjson::DumpOptions::default()));
    }
    let json_expected = parse_json(&json_text, &inputs).expect("the reference parses");
    let json_actual =
        parse_json(&(objects.join("\n") + "\n"), &inputs).expect("the candidate parses");
    let json_plain = compare(&json_expected, &json_actual, &mut no_documented());
    assert_eq!(
        number(&json_plain, &["columns", "domains", "differing"]),
        REORDERABLE_DOMAINS_OBJECTS
    );
    let json_report = compare_as_lists(&json_expected, &json_actual, &mut no_documented(), "json");
    assert_eq!(verdict(&json_report), "green");
    assert_eq!(
        number(&json_report, &["columns", "domains", "order_only"]),
        REORDERABLE_DOMAINS_OBJECTS
    );
    assert_eq!(
        order_only_fields("json"),
        fields(&[
            ("domains", JSON_ARRAY),
            ("zyg", ","),
            ("nearestexonjb", ",")
        ])
    );

    // Another column's order, or a repeated item, is not an order-only agreement.
    let domains = "DOMAINS=PDB-ENSP_mappings:2ywk.A,AFDB-ENSP_mappings:AF-P57052-F1";
    let mut rows = data_lines(&default);
    let target = rows
        .iter()
        .position(|r| r.contains(domains))
        .expect("the reference writes the row");
    rows[target] = rows[target].replace(
        domains,
        "DOMAINS=PDB-ENSP_mappings:2ywk.A,PDB-ENSP_mappings:2ywk.A,AFDB-ENSP_mappings:AF-P57052-F1",
    );
    let repeated = compare_as_lists(
        &expected,
        &default_entries(&rebuild(&default, &rows)),
        &mut no_documented(),
        "default",
    );
    assert_eq!(verdict(&repeated), "red");
    assert_eq!(number(&repeated, &["columns", "DOMAINS", "differing"]), 1);
    let mut rows = data_lines(&default);
    rows[target] = rows[target].replacen(
        "\tmissense_variant\t",
        "\tsplice_region_variant,missense_variant\t",
        1,
    );
    let other = compare_as_lists(
        &expected,
        &default_entries(&rebuild(&default, &rows)),
        &mut no_documented(),
        "default",
    );
    assert_eq!(verdict(&other), "red");
    assert_eq!(number(&other, &["columns", "Consequence", "differing"]), 1);
    assert_eq!(number(&other, &["columns", "Consequence", "order_only"]), 0);
}

/// One Extra text may carry a pair order_only accepts beside a pair a field class accepts:
/// the two acceptances meet on the text, the class credited there and order_only through the
/// pair's own column.
#[test]
fn an_extra_text_differing_in_an_order_only_pair_and_a_class_pair_is_accepted() {
    let reference = default_text(&[default_row(
        "a",
        "21:100",
        "G",
        "ENST1",
        "5_prime_UTR_variant",
        "IMPACT=MODIFIER;STRAND=1;NearestExonJB=ENSE1+10+100,ENSE2-5+200;LOEUF=-",
    )]);
    let candidate = reference
        .replace("ENSE1+10+100,ENSE2-5+200", "ENSE2-5+200,ENSE1+10+100")
        .replace("LOEUF=-", "LOEUF=0");
    let entries = (default_entries(&reference), default_entries(&candidate));
    let mut documented = documented_with(&["order_only", "loeuf_zero_minimum"]);
    let report = compare_as_lists(&entries.0, &entries.1, &mut documented, "default");
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        counts(&documented.applied),
        count_map(&[("order_only", 1), ("loeuf_zero_minimum", 1)])
    );
    // Without order_only the reordered items stay a difference; without the class the LOEUF does.
    let report = compare(
        &entries.0,
        &entries.1,
        &mut documented_with(&["loeuf_zero_minimum"]),
    );
    assert_eq!(verdict(&report), "red");
    let report = compare_as_lists(
        &entries.0,
        &entries.1,
        &mut documented_with(&["order_only"]),
        "default",
    );
    assert_eq!(verdict(&report), "red");
}

/// NearestExonJB's exon items and the 5UTR_annotation groups, whose key=value pairs
/// Ensembl VEP writes in hash order, agree as multisets under order_only in every text format
/// and NearestExonJB in the JSON format; a changed pair stays a difference.
#[test]
fn the_order_only_class_covers_exon_items_and_annotation_pair_groups() {
    let utr = "5UTR_annotation=KozakStrength=Moderate:type=uORF:Evidence=False,type=uORF:Evidence=True:KozakStrength=Weak";
    let reference = default_text(&[default_row(
        "a",
        "21:100",
        "G",
        "ENST1",
        "5_prime_UTR_variant",
        &format!("IMPACT=MODIFIER;STRAND=1;NearestExonJB=ENSE1+10+100,ENSE2-5+200;{utr}"),
    )]);
    let candidate = reference
        .replace("ENSE1+10+100,ENSE2-5+200", "ENSE2-5+200,ENSE1+10+100")
        .replace(
            utr,
            "5UTR_annotation=Evidence=True:KozakStrength=Weak:type=uORF,Evidence=False:KozakStrength=Moderate:type=uORF",
        );
    let mut documented = documented_with(&["order_only"]);
    let report = compare_as_lists(
        &default_entries(&reference),
        &default_entries(&candidate),
        &mut documented,
        "default",
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[("order_only", 3)]));
    for column in ["NearestExonJB", "5UTR_annotation"] {
        assert_eq!(number(&report, &["columns", column, "order_only"]), 1);
    }
    assert_eq!(
        verdict(&compare_default(&reference, &candidate, &[]).0),
        "red"
    );
    let changed = candidate.replace("KozakStrength=Weak", "KozakStrength=Strong");
    let report = compare_as_lists(
        &default_entries(&reference),
        &default_entries(&changed),
        &mut documented_with(&["order_only"]),
        "default",
    );
    assert_eq!(verdict(&report), "red");
    // The tab and VCF spellings: a column, and & between items.
    let tab_header = "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tNearestExonJB\t5UTR_annotation";
    let tab = format!(
        "{tab_header}\na\t21:100\tG\tENST1\tTranscript\t5_prime_UTR_variant\tENSE1+10+100,ENSE2-5+200\tk=1:t=u,k=2:t=v\n"
    );
    let swapped = tab
        .replace("ENSE1+10+100,ENSE2-5+200", "ENSE2-5+200,ENSE1+10+100")
        .replace("k=1:t=u,k=2:t=v", "t=v:k=2,t=u:k=1");
    let mut documented = documented_with(&["order_only"]);
    let report = compare_as_lists(
        &tab_entries(&tab),
        &tab_entries(&swapped),
        &mut documented,
        "tab",
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[("order_only", 2)]));
    assert_eq!(
        order_only_fields("vcf")
            .get("5UTR_annotation")
            .map(String::as_str),
        Some("pairs:&")
    );
    assert!(same_items("k=1:t=u&k=2:t=v", "t=v:k=2&t=u:k=1", "pairs:&"));
    assert!(!same_items("k=1:t=u&k=2:t=v", "t=v:k=2&t=u:k=3", "pairs:&"));
    assert_eq!(
        order_only_fields("json")
            .get("nearestexonjb")
            .map(String::as_str),
        Some(",")
    );
}

/// `CLIN_SIG` under `--clin_sig_allele 1` holds each different term list of a row's known
/// variants once, joined by `;` in the key order of a hash (OutputFactory.pm:1120,
/// 1143-1144), so two lists written in the other order agree under order_only in every row
/// format (`clin_sig_allele_list_order`): the default format's `%3B` (VEP_output.pm:118),
/// the tab format's `;` and the VCF format's `%3B` between lists whose terms `&` joins
/// (VCF.pm:401-402). A list is compared whole, so a changed term, a list written twice and
/// the terms of one list in another order stay differences.
#[test]
fn the_order_only_class_covers_the_clin_sig_lists_of_two_known_variants() {
    let reference = default_text(&[default_row(
        "a",
        "21:500",
        "G",
        "ENST1",
        "missense_variant",
        "IMPACT=MODERATE;CLIN_SIG=likely_pathogenic,pathogenic%3Blikely_pathogenic;PHENO=1,1",
    )]);
    let in_variant_order = "CLIN_SIG=likely_pathogenic%3Blikely_pathogenic,pathogenic";
    let candidate = reference.replace(
        "CLIN_SIG=likely_pathogenic,pathogenic%3Blikely_pathogenic",
        in_variant_order,
    );
    let mut documented = documented_with(&["order_only"]);
    let report = compare_as_lists(
        &default_entries(&reference),
        &default_entries(&candidate),
        &mut documented,
        "default",
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[("order_only", 2)]));
    assert_eq!(number(&report, &["columns", "CLIN_SIG", "order_only"]), 1);
    assert_eq!(number(&report, &["columns", "Extra", "order_only"]), 1);
    assert_eq!(
        verdict(&compare_default(&reference, &candidate, &[]).0),
        "red"
    );
    for changed in [
        "CLIN_SIG=likely_pathogenic%3Blikely_pathogenic,benign",
        "CLIN_SIG=likely_pathogenic%3Blikely_pathogenic%3Blikely_pathogenic,pathogenic",
        "CLIN_SIG=likely_pathogenic%3Bpathogenic,likely_pathogenic",
        "CLIN_SIG=likely_pathogenic",
    ] {
        let report = compare_as_lists(
            &default_entries(&reference),
            &default_entries(&candidate.replace(in_variant_order, changed)),
            &mut documented_with(&["order_only"]),
            "default",
        );
        assert_eq!(verdict(&report), "red", "{changed}");
    }
    // The terms of one known variant, or of `--clin_sig_allele 0`, are no lists of their own.
    let terms = reference.replace(
        "CLIN_SIG=likely_pathogenic,pathogenic%3Blikely_pathogenic",
        "CLIN_SIG=likely_pathogenic,pathogenic",
    );
    let report = compare_as_lists(
        &default_entries(&terms),
        &default_entries(&terms.replace(
            "likely_pathogenic,pathogenic",
            "pathogenic,likely_pathogenic",
        )),
        &mut documented_with(&["order_only"]),
        "default",
    );
    assert_eq!(verdict(&report), "red");
    // The tab and VCF spellings: a column with `;` between the lists, and `%3B` between lists
    // whose terms `&` joins.
    let tab_header =
        "#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tCLIN_SIG";
    let tab = format!(
        "{tab_header}\na\t21:500\tG\tENST1\tTranscript\tmissense_variant\tlikely_pathogenic,pathogenic;likely_pathogenic\n"
    );
    let swapped = tab.replace(
        "likely_pathogenic,pathogenic;likely_pathogenic",
        "likely_pathogenic;likely_pathogenic,pathogenic",
    );
    let mut documented = documented_with(&["order_only"]);
    let report = compare_as_lists(
        &tab_entries(&tab),
        &tab_entries(&swapped),
        &mut documented,
        "tab",
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(counts(&documented.applied), count_map(&[("order_only", 1)]));
    assert_eq!(number(&report, &["columns", "CLIN_SIG", "order_only"]), 1);
    for (fmt, separator) in [("default", "%3B"), ("tab", ";"), ("vcf", "%3B")] {
        assert_eq!(clin_sig_list_separator(fmt), separator);
        assert_eq!(
            order_only_fields(fmt).get("CLIN_SIG").map(String::as_str),
            Some(separator)
        );
    }
    assert!(same_items(
        "likely_pathogenic&pathogenic%3Blikely_pathogenic",
        "likely_pathogenic%3Blikely_pathogenic&pathogenic",
        "%3B"
    ));
    assert!(!same_items(
        "likely_pathogenic&pathogenic%3Blikely_pathogenic",
        "likely_pathogenic%3Bpathogenic&likely_pathogenic",
        "%3B"
    ));
    // JSON carries each known variant's own terms and no joined cell.
    assert_eq!(order_only_fields("json").get("clin_sig"), None);
    assert_eq!(order_only_fields("json").get("CLIN_SIG"), None);
}

/// The default format's `Extra` text may carry a cut list beside a pair another class
/// accepts and a pair `order_only` accepts: the reference's text is compared as the
/// candidate reads its cut list, so the three acceptances meet on it. A pair no class
/// accepts stays a difference beside the cut list, which is still counted.
#[test]
fn a_cut_list_meets_another_class_on_one_extra_text() {
    let cut_list = "colocated_list_cut_at_cache_cap";
    // 113 made-up eight-digit PubMed ids, and that list as a cache cuts it at 1,024 bytes,
    // inside a 114th id.
    let pubmed_whole: Vec<String> = (1..114).map(|n| (99_000_000 + n).to_string()).collect();
    let pubmed_whole = pubmed_whole.join(",");
    let pubmed_cut = format!("{pubmed_whole},9900011");
    let extra = "IMPACT=MODIFIER;STRAND=1;SYMBOL=GENEX1P;SYMBOL_SOURCE=EntrezGene";
    let domains = "DOMAINS=Pfam:PF99001,PANTHER:PTHR99001";
    let row = |extra: &str| {
        default_row(
            "rs9990000001",
            "21:500",
            "G",
            "ENST1",
            "intron_variant",
            extra,
        )
    };
    let stored = format!("{extra};{domains};PUBMED={pubmed_cut}");
    let reference = default_text(&[row(&stored)]);
    let read = format!("{extra};HGNC_ID=HGNC:99001;DOMAINS=PANTHER:PTHR99001,Pfam:PF99001");
    let candidate = default_text(&[row(&format!("{read};PUBMED={pubmed_whole}"))]);
    let entries = (default_entries(&reference), default_entries(&candidate));
    let mut documented = documented_with(&[cut_list, "hgnc_id_copy", "order_only"]);
    let report = compare_as_lists(&entries.0, &entries.1, &mut documented, "default");
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(cut_list, 1), ("hgnc_id_copy", 1), ("order_only", 1)])
    );
    assert_eq!(
        at(&report, &["columns", "Extra", "masked"]),
        &json(r#"{"colocated_list_cut_at_cache_cap": 1, "hgnc_id_copy": 1}"#)
    );
    let mut documented = documented_with(&["hgnc_id_copy", "order_only"]);
    let report = compare_as_lists(&entries.0, &entries.1, &mut documented, "default");
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "PUBMED", "differing"]), 1);
    let planted = default_entries(&candidate.replace("IMPACT=MODIFIER", "IMPACT=LOW"));
    let mut documented = documented_with(&[cut_list, "hgnc_id_copy", "order_only"]);
    let report = compare_as_lists(&entries.0, &planted, &mut documented, "default");
    assert_eq!(verdict(&report), "red");
    assert_eq!(documented.applied.get(cut_list), 1);
    let columns = at(&report, &["columns"]).as_dict().expect("an object");
    let differing: BTreeSet<&str> = columns
        .iter()
        .filter(|(_, column)| column.get("differing") != Some(&Value::Int(0)))
        .map(|(name, _)| name.as_str())
        .collect();
    assert_eq!(differing, ["Extra", "IMPACT"].into_iter().collect());
}
