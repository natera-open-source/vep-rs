// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Golden-corpus tests: every consequence-set combination observed on the
//! datasets listed in `manuscript/data/dataset_inventory.csv`, per Ensembl
//! release and assembly, annotated with the corpus's own pruned cache and
//! compared column by column with the Ensembl VEP output committed beside it.
//!
//! One run per corpus and per `flag_sets` entry of its manifest; every mismatch
//! is collected and reported in one failure message, grouped by column and
//! consequence set, so a change is measured against the whole corpus rather
//! than its first failing row.

mod common;

use std::collections::BTreeSet;
use std::fs;
use std::io::Write;
use std::path::Path;

use common::{
    compare_run, compare_warnings, corpora, parse_default, parse_json, parse_vcf, read_expected,
    runs, Corpus, REFERENCE_SKIPPED_RECORD,
};

/// Every run's default-format output matches the reference, and the warnings
/// file the run writes beside it names the records the reference skipped, with
/// their reasons (the manifest's `reference_warning` strings), and no other.
#[test]
fn default_output_matches_vep_on_every_corpus() {
    let corpora = corpora();
    assert!(!corpora.is_empty(), "no corpus found under tests/golden");
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora {
        for run in runs(corpus).iter().filter(|r| r.compares("default")) {
            let (header, body) = compare_run(run, "default", tmp.path());
            if !header.is_empty() || !body.is_empty() {
                failures.push(format!(
                    "== corpus {} (default format)\n{header}{body}",
                    run.label()
                ));
            }
            let warnings = compare_warnings(run, tmp.path());
            if !warnings.is_empty() {
                failures.push(format!(
                    "== corpus {} (warnings file)\n{warnings}",
                    run.label()
                ));
            }
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

/// The manifest is internally consistent with the committed inputs and
/// expected output: every record it lists is in `variants.vcf`, every flag
/// set's input carries one data line per record, its `covers` names at least
/// one parity-matrix row, every combination it claims to cover appears in VEP's
/// output for the corpus, and the generator left no combination without an
/// exemplar.
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
        let covers = m["covers"].as_array();
        assert!(
            covers
                .is_some_and(|c| !c.is_empty()
                    && c.iter().all(|n| n.as_str().is_some_and(|s| !s.is_empty()))),
            "{}/{}: the manifest's covers must name at least one parity-matrix row",
            corpus.release,
            corpus.name
        );
        for run in runs(&corpus).iter().filter(|r| r.input != "variants.vcf") {
            assert_eq!(
                records.len(),
                run.input_records().len(),
                "{}: manifest records vs {} data lines",
                run.label(),
                run.input
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

fn copy_dir(src: &Path, dst: &Path) {
    fs::create_dir_all(dst).unwrap();
    for entry in fs::read_dir(src).unwrap().flatten() {
        let target = dst.join(entry.file_name());
        if entry.path().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), &target).unwrap();
        }
    }
}

/// A corpus runs once per `flag_sets` entry beside its own run. On a copy of
/// the smallest committed corpus, three entries point at the manifest's
/// expected outputs: one repeating its flags as written, one over a
/// bgzip-compressed copy of the input, one with the input piped on standard
/// input and `--format vcf` added to the flags, the input format a stdin run
/// must name; each compares the formats it names and no other, the compressed
/// input carries one data line per manifest record, and the entries inherit the
/// manifest's documented divergences. A fourth entry carrying an empty
/// `divergences` list of its own is compared without them, so the reference's
/// skipped records make it fail.
#[test]
fn flag_sets_run_once_per_entry_with_their_input_and_stdin() {
    let source = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-release116")
        .expect("the GRCh38-release116 corpus is committed");
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join("corpus");
    copy_dir(&source.dir, &dir);
    let plain = fs::read(dir.join("variants.vcf")).unwrap();
    // Two bgzf members, as a bgzip file over 64 KiB has: a single-member reader
    // would stop at the first and read half the records.
    let mut bgzip =
        noodles_bgzf::io::Writer::new(fs::File::create(dir.join("variants.vcf.gz")).unwrap());
    let (head, tail) = plain.split_at(plain.len() / 2);
    bgzip.write_all(head).unwrap();
    bgzip.flush().unwrap();
    bgzip.write_all(tail).unwrap();
    bgzip.finish().unwrap();
    let flags = source.manifest["flags"].clone();
    let mut stdin_flags = flags.as_array().unwrap().clone();
    stdin_flags.extend(["--format", "vcf"].map(serde_json::Value::from));
    let mut manifest = source.manifest.clone();
    manifest["flag_sets"] = serde_json::json!([
        {"name": "again", "flags": flags, "formats": ["default", "tab", "vcf", "json"], "expected_dir": "expected", "covers": []},
        {"name": "bgzip", "flags": flags, "input": "variants.vcf.gz", "formats": ["default", "json"], "expected_dir": "expected"},
        {"name": "stdin", "flags": stdin_flags, "stdin": true, "formats": ["vcf"], "expected_dir": "expected"},
        {"name": "strict", "flags": flags, "formats": ["default"], "expected_dir": "expected",
         "divergences": [], "vep_rs_only_tuples": []},
    ]);
    let corpus = Corpus {
        release: source.release.clone(),
        name: source.name.clone(),
        assembly: source.assembly.clone(),
        dir,
        manifest,
    };
    let runs = runs(&corpus);
    let names: Vec<&str> = runs.iter().map(|r| r.name.as_str()).collect();
    assert_eq!(names, ["default", "again", "bgzip", "stdin", "strict"]);
    assert_eq!(runs[0].label(), "116/GRCh38-release116");
    assert_eq!(runs[2].label(), "116/GRCh38-release116 [bgzip]");
    assert!(runs[2].compares("json") && !runs[2].compares("tab"));
    assert!(runs[3].stdin && !runs[2].stdin);
    let records = corpus.manifest["records"].as_array().unwrap().len();
    assert_eq!(runs[2].input_records().len(), records);
    assert_eq!(runs[2].input_records(), runs[0].input_records());

    let mut failures = Vec::new();
    for run in &runs[..4] {
        for format in ["default", "tab", "vcf", "json"] {
            if !run.compares(format) {
                continue;
            }
            let (header, body) = compare_run(run, format, tmp.path());
            if !header.is_empty() || !body.is_empty() {
                failures.push(format!("== {} ({format})\n{header}{body}", run.label()));
            }
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));

    let (_header, report) = compare_run(&runs[4], "default", tmp.path());
    assert!(
        report.contains("UNEXPECTED in vep-rs"),
        "an entry with its own empty divergence lists is compared without the manifest's:\n{report}"
    );
}

/// `runs` refuses the `flag_sets` entries the corpus tool's `flag-set` refuses,
/// naming the entry and the reason, so a manifest edited by hand cannot certify
/// a set's `covers` with a run the tool would not have written; the shapes the
/// tool writes (a name of letters, digits, `_` and `-` that starts with a letter
/// or digit) run.
#[test]
fn runs_refuses_the_flag_set_entries_the_corpus_tool_refuses() {
    let source = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-release116")
        .expect("the GRCh38-release116 corpus is committed");
    let label = format!("{}/{}", source.release, source.name);
    let with_flag_sets = |flag_sets: serde_json::Value| {
        let mut manifest = source.manifest.clone();
        manifest["flag_sets"] = flag_sets;
        Corpus {
            release: source.release.clone(),
            name: source.name.clone(),
            assembly: source.assembly.clone(),
            dir: source.dir.clone(),
            manifest,
        }
    };
    let not_a_name = |name: &str| {
        format!("{label}: flag set name {name:?} is not a plain word other than default")
    };
    let outside =
        |path: &str| format!("{label} flag set x: {path:?} is not a path inside the corpus");
    let refused = [
        (
            serde_json::json!([{"name": "default"}]),
            not_a_name("default"),
        ),
        (serde_json::json!([{"name": ""}]), not_a_name("")),
        (serde_json::json!([{"name": "-x"}]), not_a_name("-x")),
        (serde_json::json!([{"name": "_x"}]), not_a_name("_x")),
        (serde_json::json!([{"name": "a/b"}]), not_a_name("a/b")),
        (
            serde_json::json!([{"name": "x"}, {"name": "x"}]),
            format!("{label}: flag set x is named twice"),
        ),
        (
            serde_json::json!([{"name": "x", "formats": []}]),
            format!("{label} flag set x: no format to compare"),
        ),
        (
            serde_json::json!([{"name": "x", "formats": ["parquet"]}]),
            format!("{label} flag set x: unknown format parquet"),
        ),
        (serde_json::json!([{"name": "x", "input": ""}]), outside("")),
        (
            serde_json::json!([{"name": "x", "input": "../variants.vcf"}]),
            outside("../variants.vcf"),
        ),
        (
            serde_json::json!([{"name": "x", "input": "/variants.vcf"}]),
            outside("/variants.vcf"),
        ),
        (
            serde_json::json!([{"name": "x", "expected_dir": ""}]),
            outside(""),
        ),
        (
            serde_json::json!([{"name": "x", "expected_dir": "../expected"}]),
            outside("../expected"),
        ),
        (
            serde_json::json!([{"name": "x", "stdin": true, "flags": ["--symbol"]}]),
            format!("{label} flag set x: a stdin run names its input format with --format"),
        ),
        (
            serde_json::json!([{"name": "x", "stdin": "true"}]),
            format!("{label} flag set x: stdin is not of its type"),
        ),
        (
            serde_json::json!([{"name": "x", "unordered": 1}]),
            format!("{label} flag set x: unordered is not of its type"),
        ),
        (
            serde_json::json!([{"name": "x", "input": 5}]),
            format!("{label} flag set x: input is not of its type"),
        ),
        (
            serde_json::json!([{"name": "x", "formats": "vcf"}]),
            format!("{label} flag set x: formats is not of its type"),
        ),
        (
            serde_json::json!([{"name": "x", "flags": "--symbol"}]),
            format!("{label} flag set x: flags is not of its type"),
        ),
    ];
    let mut wrong = Vec::new();
    for (flag_sets, reason) in refused {
        let corpus = with_flag_sets(flag_sets.clone());
        let outcome =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| runs(&corpus).len()));
        let message = match &outcome {
            Ok(count) => format!("accepted, {count} runs"),
            Err(panic) => panic
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| panic.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_default(),
        };
        if message != reason {
            wrong.push(format!("{flag_sets}: expected {reason:?}, got {message:?}"));
        }
    }
    assert!(wrong.is_empty(), "\n{}", wrong.join("\n"));
    let corpus = with_flag_sets(serde_json::json!([
        {"name": "x1_-", "flags": ["--format", "vcf"], "stdin": true, "formats": ["vcf"], "expected_dir": "expected"},
        {"name": "eq", "flags": ["--format=vcf"], "stdin": true},
    ]));
    assert_eq!(runs(&corpus).len(), 3);
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
/// ten mebibytes and would otherwise consume the corpus budget by itself; a
/// corpus naming a sibling's reference (`../GRCh37-hgvs/reference.fa.gz`) ships
/// no bytes of it, so nothing is subtracted.
#[test]
fn each_corpus_fits_its_size_budget() {
    const BUDGET_BYTES: u64 = 15 * 1024 * 1024;
    const FASTA_BUDGET_BYTES: u64 = 11 * 1024 * 1024;
    for corpus in corpora() {
        let mut size = common::dir_size(&corpus.dir);
        if let Some(fasta) = corpus.manifest["fasta"].as_str() {
            let path = corpus.dir.join(fasta);
            let fasta_size = fs::metadata(&path).unwrap().len();
            assert!(
                fasta_size <= FASTA_BUDGET_BYTES,
                "{}/{}: {fasta} is {fasta_size} bytes, over the {FASTA_BUDGET_BYTES} byte reference budget",
                corpus.release,
                corpus.name
            );
            let inside_corpus = path
                .canonicalize()
                .unwrap()
                .starts_with(corpus.dir.canonicalize().unwrap());
            if inside_corpus {
                size -= fasta_size;
            }
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
