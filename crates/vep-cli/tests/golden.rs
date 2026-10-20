// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Golden-corpus tests: every consequence-set combination observed on the
//! datasets listed in `manuscript/data/dataset_inventory.csv`, per Ensembl
//! release and assembly, annotated with the corpus's own pruned cache and
//! compared column by column with the Ensembl VEP output committed beside it.
//!
//! One run per corpus; every mismatch is collected and reported in one failure
//! message, grouped by column and consequence set, so a change is measured
//! against the whole corpus rather than its first failing row.

mod common;

use std::collections::BTreeSet;

use common::{
    compare_entries, compare_tab_headers, corpora, parse_default, parse_json, parse_vcf,
    read_expected, run_vep, Documented, REFERENCE_SKIPPED_RECORD,
};

#[test]
fn default_output_matches_vep_on_every_corpus() {
    let corpora = corpora();
    assert!(!corpora.is_empty(), "no corpus found under tests/golden");
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora {
        let expected = read_expected(&corpus.dir, "default.txt");
        let actual = run_vep(corpus, "default", tmp.path(), &[]);
        let documented = Documented::from_manifest(&corpus.manifest);
        let header = compare_tab_headers(&expected, &actual);
        let body = compare_entries(
            &parse_default(&expected),
            &parse_default(&actual),
            &documented,
        );
        if !header.is_empty() || !body.is_empty() {
            failures.push(format!(
                "== corpus {}/{} (default format)\n{header}{body}",
                corpus.release, corpus.name
            ));
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

/// `--hgvsg` on the GRCh37-hgvs corpus: every row of a record carries the HGVSg
/// Ensembl VEP 116 prints for it (`vep --hgvsg --fasta reference.fa` on the
/// corpus input and reference), in the tab, VCF and JSON (`hgvsg`) formats alike:
/// an SNV, a deletion left where it is, a deletion and an insertion moved 3'
/// along the forward strand (one a `dup`), an insertion that duplicates nothing,
/// and a `delins`. The committed expected output was written without the flag,
/// so the values are listed here.
#[test]
fn hgvsg_matches_vep_in_every_format() {
    // (the record's first five VCF columns, Location, Allele, HGVSg)
    const EXPECTED: [(&str, &str, &str, &str); 6] = [
        (
            "21\t9937918\trs867932594\tG\tA",
            "21:9937918",
            "A",
            "21:g.9937918G>A",
        ),
        (
            "21\t9937918\trs867932594\tGGAGGGTCTGAATGTTTGTCCCTCACATAGGATTCCAGAACACTGCTACAAGGT\tG",
            "21:9937919-9937971",
            "-",
            "21:g.9937919_9937971del",
        ),
        (
            "21\t10910400\trs149363218\tCTA\tC",
            "21:10910401-10910402",
            "-",
            "21:g.10910402_10910403del",
        ),
        (
            "21\t18981461\trs768550367\tA\tATT",
            "21:18981461-18981462",
            "TT",
            "21:g.18981467_18981468dup",
        ),
        (
            "21\t19653357\trs62213417\tC\tCTT",
            "21:19653357-19653358",
            "TT",
            "21:g.19653358_19653359insTT",
        ),
        (
            "21\t27369678\t662051\tTAACAGGATC\tG",
            "21:27369678-27369687",
            "G",
            "21:g.27369678_27369687delinsG",
        ),
    ];
    let corpus = corpora()
        .into_iter()
        .find(|c| c.release == "116" && c.name == "GRCh37-hgvs")
        .expect("tests/golden/116/GRCh37-hgvs");
    let tmp = tempfile::tempdir().unwrap();
    let inputs = common::input_records(&corpus);
    let tab = parse_default(&run_vep(&corpus, "tab", tmp.path(), &["--hgvsg"]));
    let vcf = parse_vcf(&run_vep(&corpus, "vcf", tmp.path(), &["--hgvsg"]), &inputs);
    let json = parse_json(&run_vep(&corpus, "json", tmp.path(), &["--hgvsg"]), &inputs);
    for (line, location, allele, hgvsg) in EXPECTED {
        let record = inputs
            .iter()
            .position(|input| input.starts_with(&format!("{line}\t")))
            .unwrap_or_else(|| panic!("no input record {line}"));
        let rows = [
            (
                "tab",
                "HGVSg",
                tab.iter()
                    .filter(|e| e.location.as_deref() == Some(location) && e.allele == allele)
                    .collect::<Vec<_>>(),
            ),
            (
                "vcf",
                "HGVSg",
                vcf.iter()
                    .filter(|e| e.record == Some(record) && e.allele == allele)
                    .collect(),
            ),
            (
                "json",
                "hgvsg",
                json.iter()
                    .filter(|e| e.record == Some(record) && e.allele == allele)
                    .collect(),
            ),
        ];
        for (format, field, entries) in rows {
            assert!(
                !entries.is_empty(),
                "{format}: no row for {location} {allele}"
            );
            for e in entries {
                assert_eq!(
                    e.fields.get(field).map(String::as_str),
                    Some(hgvsg),
                    "{format}: {}",
                    e.key
                );
            }
        }
    }
}

/// The manifest is internally consistent with the committed inputs and
/// expected output: every record it lists is in `variants.vcf`, every
/// combination it claims to cover appears in VEP's output for the corpus, and
/// the generator left no combination without an exemplar.
#[test]
fn manifest_covers_every_combination_it_claims() {
    for corpus in corpora() {
        let m = &corpus.manifest;
        let records = m["records"].as_array().unwrap();
        let inputs = common::input_records(&corpus);
        assert_eq!(
            records.len(),
            inputs.len(),
            "{}/{}: manifest records vs variants.vcf data lines",
            corpus.release,
            corpus.name
        );
        for (i, (rec, line)) in records.iter().zip(inputs.iter()).enumerate() {
            assert_eq!(
                rec["vcf"].as_str().unwrap(),
                line,
                "{}/{}: record {i} differs between manifest and variants.vcf",
                corpus.release,
                corpus.name
            );
        }
        assert!(
            m["uncovered_combinations"].as_array().unwrap().is_empty(),
            "{}/{}: combinations without an exemplar: {:?}",
            corpus.release,
            corpus.name,
            m["uncovered_combinations"]
        );
        let expected = read_expected(&corpus.dir, "default.txt");
        let observed: BTreeSet<String> = parse_default(&expected)
            .iter()
            .map(|e| common::consequence_set(&e.consequence))
            .collect();
        let claimed: BTreeSet<String> = m["combinations"]
            .as_object()
            .unwrap()
            .keys()
            .map(|k| common::consequence_set(k))
            .collect();
        let unmet: Vec<&String> = claimed.difference(&observed).collect();
        assert!(
            unmet.is_empty(),
            "{}/{}: {} combination(s) in the manifest never appear in expected/default.txt: {:?}",
            corpus.release,
            corpus.name,
            unmet.len(),
            unmet
        );
        let listed: BTreeSet<String> = records
            .iter()
            .flat_map(|r| r["combinations"].as_array().unwrap().iter())
            .map(|c| common::consequence_set(c.as_str().unwrap()))
            .collect();
        assert_eq!(
            listed, claimed,
            "{}/{}: the records' combinations and the combination table disagree",
            corpus.release, corpus.name
        );
    }
}

/// The records the reference dropped before annotation are marked consistently
/// with its committed output: `reference_rows` 0 and a `reference_warning`
/// string go together; every `vep_rs_only_tuples` entry classed
/// `reference_skipped_record` names only marked records; the reference wrote no
/// VCF line and no JSON object for a marked record, and one of each for every
/// other record (Ensembl VEP writes one per record it annotates).
#[test]
fn skipped_record_markers_match_the_reference_output() {
    for corpus in corpora() {
        let label = format!("{}/{}", corpus.release, corpus.name);
        let m = &corpus.manifest;
        let records = m["records"].as_array().unwrap();
        let marked: BTreeSet<usize> = records
            .iter()
            .enumerate()
            .filter(|(_, r)| r["reference_rows"].as_u64() == Some(0))
            .map(|(i, _)| i)
            .collect();
        for (i, rec) in records.iter().enumerate() {
            assert_eq!(
                marked.contains(&i),
                rec["reference_warning"].is_string(),
                "{label}: record {i} has one of reference_rows 0 and reference_warning without the other"
            );
        }
        for item in m["vep_rs_only_tuples"].as_array().into_iter().flatten() {
            if item["expected_divergence"].as_str() != Some(REFERENCE_SKIPPED_RECORD) {
                continue;
            }
            for r in item["record_indices"].as_array().unwrap() {
                let r = r.as_u64().unwrap() as usize;
                assert!(
                    marked.contains(&r),
                    "{label}: a {REFERENCE_SKIPPED_RECORD} tuple names record {r}, which has reference rows"
                );
            }
        }
        let inputs = common::input_records(&corpus);
        let vcf = parse_vcf(&read_expected(&corpus.dir, "vcf.vcf"), &inputs);
        let json = parse_json(&read_expected(&corpus.dir, "json.jsonl"), &inputs);
        for e in vcf.iter().chain(json.iter()) {
            assert!(
                e.record.is_some(),
                "{label}: reference entry {} matches no input record",
                e.key
            );
        }
        let in_vcf: BTreeSet<usize> = vcf.iter().filter_map(|e| e.record).collect();
        let in_json: BTreeSet<usize> = json.iter().filter_map(|e| e.record).collect();
        for i in 0..records.len() {
            let (vcf_line, json_object) = (in_vcf.contains(&i), in_json.contains(&i));
            if marked.contains(&i) {
                assert!(
                    !vcf_line && !json_object,
                    "{label}: record {i} is marked reference_rows 0 but the reference wrote a VCF line ({vcf_line}) or a JSON object ({json_object}) for it"
                );
            } else {
                assert!(
                    vcf_line && json_object,
                    "{label}: record {i} is not marked reference_rows 0 but the reference wrote no VCF line ({vcf_line}) or no JSON object ({json_object}) for it"
                );
            }
        }
    }
}

/// Each corpus stays inside its in-repository size budget. A corpus that ships a
/// reference contig (`manifest["fasta"]`) has that file measured on its own
/// against a budget of its own, since a whole chromosome compresses to about
/// ten mebibytes and would otherwise consume the corpus budget by itself.
#[test]
fn each_corpus_fits_its_size_budget() {
    const BUDGET_BYTES: u64 = 15 * 1024 * 1024;
    const FASTA_BUDGET_BYTES: u64 = 11 * 1024 * 1024;
    for corpus in corpora() {
        let mut size = common::dir_size(&corpus.dir);
        if let Some(fasta) = corpus.manifest["fasta"].as_str() {
            let fasta_size = std::fs::metadata(corpus.dir.join(fasta)).unwrap().len();
            assert!(
                fasta_size <= FASTA_BUDGET_BYTES,
                "{}/{}: {fasta} is {fasta_size} bytes, over the {FASTA_BUDGET_BYTES} byte reference budget",
                corpus.release,
                corpus.name
            );
            size -= fasta_size;
        }
        assert!(
            size <= BUDGET_BYTES,
            "{}/{} is {} bytes on disk (reference excluded), over the {} byte budget",
            corpus.release,
            corpus.name,
            size,
            BUDGET_BYTES
        );
    }
}
