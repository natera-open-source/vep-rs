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
    apply_masks, compare_entries, compare_run, compare_warnings, corpora,
    list_columns_as_multisets, nearest_as_coordinates, parse_default, parse_json, parse_vcf,
    protein_coding_start_sites, read_expected, runs, Corpus, Documented, Entry, ORDER_ONLY,
    REFERENCE_SKIPPED_RECORD, REGULATORY_FEATURE_ROWS,
};

/// Every run's default-format output matches the reference, and the warnings
/// file the run writes beside it names the records the reference skipped, with
/// their reasons (the manifest's `reference_warning` strings), and the warnings
/// the reference wrote about the run as a whole (its `reference_warnings`
/// list), and no other.
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
/// a set's `covers` with a run the tool would not have written, and refuses a
/// manifest whose own `masks` or `reference_warnings` is not a list; the shapes
/// the tool writes (a name of letters, digits, `_` and `-` that starts with a
/// letter or digit) run.
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
        (
            serde_json::json!([{"name": "x", "perl_mode": true}]),
            format!("{label} flag set x: perl_mode is not of its type"),
        ),
        (
            serde_json::json!([{"name": "x", "perl_mode": "cache"}]),
            format!("{label} flag set x: perl_mode is offline or database, not \"cache\""),
        ),
        (
            serde_json::json!([{"name": "x", "reference_warnings": "text"}]),
            format!("{label} flag set x: reference_warnings is not of its type"),
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
    for (key, value) in [
        ("reference_warnings", serde_json::json!("text")),
        ("masks", serde_json::json!("order_only")),
    ] {
        let mut corpus = with_flag_sets(serde_json::json!([]));
        corpus.manifest[key] = value;
        let outcome =
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| runs(&corpus).len()));
        assert_eq!(
            outcome.unwrap_err().downcast_ref::<String>().cloned(),
            Some(format!("{label}: {key} is not of its type"))
        );
    }
    let corpus = with_flag_sets(serde_json::json!([
        {"name": "x1_-", "flags": ["--format", "vcf"], "stdin": true, "formats": ["vcf"], "expected_dir": "expected"},
        {"name": "eq", "flags": ["--format=vcf"], "stdin": true},
        {"name": "db", "flags": ["--ga4gh_vrs"], "formats": ["default", "tab", "json"], "perl_mode": "database"},
    ]));
    let runs = runs(&corpus);
    assert_eq!(runs.len(), 4);
    assert_eq!(runs[3].perl_mode, "database");
    assert_eq!(runs[1].perl_mode, "offline");
}

/// A `reference_warnings` text that is empty, or blank, is refused with the
/// corpus or the flag set named, as the corpus tool refuses it (a text of
/// nothing is a warning no run can write); a text that is not a string is
/// refused as the other typed keys are; a non-empty text runs.
#[test]
fn runs_refuses_an_empty_reference_warning_text() {
    let source = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-release116")
        .expect("the GRCh38-release116 corpus is committed");
    let label = format!("{}/{}", source.release, source.name);
    let with_warning_lists = |manifest_list: serde_json::Value, set_list: serde_json::Value| {
        let mut manifest = source.manifest.clone();
        manifest["reference_warnings"] = manifest_list;
        manifest["flag_sets"] = serde_json::json!([{"name": "x", "reference_warnings": set_list}]);
        Corpus {
            release: source.release.clone(),
            name: source.name.clone(),
            assembly: source.assembly.clone(),
            dir: source.dir.clone(),
            manifest,
        }
    };
    let panic_of = |corpus: &Corpus| -> Option<String> {
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| runs(corpus).len()))
            .err()
            .and_then(|p| p.downcast_ref::<String>().cloned())
    };
    let notice = "could not fetch sequence for tandem repeats; consequence calculation for tandem repeats is less precise when using --offline without --fasta";
    let empty = |owner: &str| Some(format!("{owner}: reference_warnings has an empty text"));
    let set_label = format!("{label} flag set x");
    assert_eq!(
        panic_of(&with_warning_lists(
            serde_json::json!([notice, ""]),
            serde_json::json!([notice])
        )),
        empty(&label)
    );
    assert_eq!(
        panic_of(&with_warning_lists(
            serde_json::json!([notice]),
            serde_json::json!([notice, " "])
        )),
        empty(&set_label)
    );
    assert_eq!(
        panic_of(&with_warning_lists(
            serde_json::json!([5]),
            serde_json::json!([notice])
        )),
        Some(format!("{label}: reference_warnings is not of its type"))
    );
    assert_eq!(
        panic_of(&with_warning_lists(
            serde_json::json!([notice]),
            serde_json::json!([notice, 5])
        )),
        Some(format!(
            "{set_label}: reference_warnings is not of its type"
        ))
    );
    assert_eq!(
        panic_of(&with_warning_lists(
            serde_json::json!([notice]),
            serde_json::json!([notice])
        )),
        None
    );
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

/// Under `--nearest transcript` a run's `NEAREST` values compare as the start
/// sites of the named transcripts: the reference names one transcript at a
/// winning site and the output every transcript there, so the ids differ where
/// the sites agree; a transcript starting at another site, and one the cache
/// lacks, still differ. The identifiers corpus's `nearest-transcript` set carries
/// such records (GRIK1's seven transcripts share the start site at 29940022).
#[test]
fn nearest_transcript_values_compare_by_start_site() {
    let corpus = corpora()
        .into_iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-identifiers")
        .expect("the identifiers corpus");
    let sites = protein_coding_start_sites(&corpus);
    assert_eq!(sites.get("ENST00000389125"), Some(&29_940_022));
    assert_eq!(sites.get("ENST00000399914"), Some(&29_940_022));
    let set = runs(&corpus)
        .into_iter()
        .find(|r| r.name == "nearest-transcript")
        .expect("the nearest-transcript set");
    assert!(set.nearest_by_coordinate());
    assert!(!common::default_run(&corpus).nearest_by_coordinate());

    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
    let row = |nearest: &str| {
        format!(
            "{header}del25\t21:30000001-30000025\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;NEAREST={nearest}\n"
        )
    };
    let resolved = |nearest: &str| {
        let mut entries = parse_default(&row(nearest));
        nearest_as_coordinates(&mut entries, &sites);
        (
            entries[0].fields["NEAREST"].clone(),
            entries[0].fields["Extra"].clone(),
        )
    };
    let one = resolved("ENST00000399914");
    let all = resolved(
        "ENST00000389124,ENST00000389125,ENST00000399907,ENST00000399909,ENST00000399913,ENST00000399914,ENST00001131791",
    );
    assert_eq!(one, all);
    assert_eq!(one.0, "29940022");
    assert_eq!(one.1, "IMPACT=MODIFIER;NEAREST=29940022");
    // The VCF format joins the ids with `&`.
    let (two, _) = resolved("ENST00000389125&ENST00000400546");
    assert_eq!(two, format!("29940022,{}", sites["ENST00000400546"]));
    assert_ne!(resolved("ENST00000400546").0, one.0);
    assert_eq!(resolved("ENST_UNKNOWN").0, "ENST_UNKNOWN");
}

/// Under `regulatory_feature_rows` the reference's RegulatoryFeature and
/// MotifFeature rows go in the text formats and the VCF, a JSON record that
/// carried a regulatory array loses it and has its most severe term recomputed
/// from the rows left (`intergenic_variant` for a record Ensembl VEP annotated
/// against a regulatory feature alone), and a record that lost nothing keeps the
/// reference's term as written; a corpus without the mask keeps every row.
#[test]
fn regulatory_feature_rows_mask_drops_the_reference_rows() {
    let masked = corpora()
        .into_iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-everything")
        .expect("the everything corpus");
    assert_eq!(masked.masks(), [REGULATORY_FEATURE_ROWS]);
    let unmasked = corpora()
        .into_iter()
        .find(|c| c.release == "116" && c.name == "GRCh38")
        .expect("the GRCh38 corpus");
    assert!(unmasked.masks().is_empty());

    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
    let text = format!(
        "{header}rs1\t21:100\tA\t-\tENSR00000001\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n\
         rs1\t21:100\tA\t-\tENSM00000001\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n\
         rs1\t21:100\tA\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n\
         rs2\t21:200\tG\tENSG1\tENST1\tTranscript\tmissense_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODERATE\n"
    );
    let mut entries = parse_default(&text);
    apply_masks(&unmasked.masks(), "default", &mut entries);
    assert_eq!(entries.len(), 4);
    apply_masks(&masked.masks(), "default", &mut entries);
    assert_eq!(
        entries.iter().map(|e| e.key.as_str()).collect::<Vec<_>>(),
        ["21:100|A|-|-", "21:200|G|ENST1|Transcript"]
    );

    let inputs = vec![
        "21\t100\trs1\tC\tA\t.\t.\t.".to_string(),
        "21\t200\trs2\tT\tG\t.\t.\t.".to_string(),
    ];
    let vcf = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|Feature_type|Feature\">\n\
               #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
               21\t100\trs1\tC\tA\t.\t.\tCSQ=A|regulatory_region_variant|MODIFIER|RegulatoryFeature|ENSR00000001,A|intergenic_variant|MODIFIER||\n\
               21\t200\trs2\tT\tG\t.\t.\tCSQ=G|missense_variant|MODERATE|Transcript|ENST1\n";
    let mut entries = parse_vcf(vcf, &inputs);
    apply_masks(&masked.masks(), "vcf", &mut entries);
    assert_eq!(
        entries.iter().map(|e| e.key.as_str()).collect::<Vec<_>>(),
        [
            "rec0|record",
            "rec0|A||",
            "rec1|record",
            "rec1|G|ENST1|Transcript"
        ]
    );

    let json = "{\"input\":\"21\\t100\\trs1\\tC\\tA\\t.\\t.\\t.\",\"most_severe_consequence\":\"regulatory_region_variant\",\
                \"regulatory_feature_consequences\":[{\"variant_allele\":\"A\",\"regulatory_feature_id\":\"ENSR00000001\",\"consequence_terms\":[\"regulatory_region_variant\"]}],\
                \"intergenic_consequences\":[{\"variant_allele\":\"A\",\"consequence_terms\":[\"intergenic_variant\"]}]}\n\
                {\"input\":\"21\\t200\\trs2\\tT\\tG\\t.\\t.\\t.\",\"most_severe_consequence\":\"splice_region_variant\",\
                \"transcript_consequences\":[{\"variant_allele\":\"G\",\"transcript_id\":\"ENST1\",\"consequence_terms\":[\"missense_variant\"]}]}\n";
    let mut entries = parse_json(json, &inputs);
    assert_eq!(
        entries[0].fields["regulatory_feature_consequences"],
        "{\"consequence_terms\":[\"regulatory_region_variant\"],\"regulatory_feature_id\":\"ENSR00000001\",\"variant_allele\":\"A\"}"
    );
    apply_masks(&masked.masks(), "json", &mut entries);
    assert!(!entries[0]
        .fields
        .contains_key("regulatory_feature_consequences"));
    assert_eq!(
        entries[0].fields["most_severe_consequence"],
        "intergenic_variant"
    );
    // The second record lost no row, so a reference term that disagrees with its
    // rows stays under comparison.
    assert_eq!(
        entries[2].fields["most_severe_consequence"],
        "splice_region_variant"
    );
    assert_eq!(entries.len(), 4);
}

/// Under `order_only` the items of the columns [`ORDER_ONLY`] names compare as
/// a multiset: the same items in another order agree in the default format (the
/// field and the raw `Extra` text), in the VCF's `&`-joined form and in the JSON
/// `domains` array, while a changed item, or one listed a different number of
/// times, still differs; an empty value stays as written. The native-builder
/// corpus's `domains` set declares the class for its run alone.
#[test]
fn order_only_mask_compares_domains_as_a_multiset() {
    let corpus = corpora()
        .into_iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-native-builder")
        .expect("the native-builder corpus");
    assert!(corpus.masks().is_empty());
    assert!(!common::default_run(&corpus).order_only());
    let set = runs(&corpus)
        .into_iter()
        .find(|r| r.name == "domains")
        .expect("the domains set");
    assert_eq!(set.masks, [ORDER_ONLY]);
    assert!(set.order_only() && !set.masks_regulatory_rows());

    let documented = Documented::from_manifest(&serde_json::json!({}));
    let inputs = vec!["21\t18256600\t.\tA\tG\t.\t.\t.".to_string()];
    let compare = |format: &str, expected: &str, actual: &str| -> String {
        let parse = |text: &str| -> Vec<Entry> {
            let mut entries = match format {
                "default" => parse_default(text),
                "vcf" => parse_vcf(text, &inputs),
                _ => parse_json(text, &inputs),
            };
            list_columns_as_multisets(&mut entries, format);
            entries
        };
        compare_entries(&parse(expected), &parse(actual), &documented, true)
    };

    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
    let row = |domains: &str| {
        format!(
            "{header}rs1\t21:18256600\tG\tENSG1\tENST1\tTranscript\tmissense_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODERATE;DOMAINS={domains};SYMBOL=CHODL\n"
        )
    };
    let reference = row("Gene3D:3.10.100.10,Pfam:PF00059,SMART:SM00034,Pfam:PF00059");
    assert_eq!(
        compare(
            "default",
            &reference,
            &row("SMART:SM00034,Pfam:PF00059,Pfam:PF00059,Gene3D:3.10.100.10")
        ),
        ""
    );
    let changed = compare(
        "default",
        &reference,
        &row("SMART:SM00034,Pfam:PF00060,Pfam:PF00059,Gene3D:3.10.100.10"),
    );
    assert!(
        changed.contains("FIELD DOMAINS: 1 mismatching"),
        "{changed}"
    );
    assert!(changed.contains("FIELD Extra: 1 mismatching"), "{changed}");
    let fewer = compare(
        "default",
        &reference,
        &row("SMART:SM00034,Pfam:PF00059,Gene3D:3.10.100.10"),
    );
    assert!(fewer.contains("FIELD DOMAINS: 1 mismatching"), "{fewer}");
    let mut empty = parse_default(&row("-"));
    list_columns_as_multisets(&mut empty, "default");
    assert_eq!(empty[0].fields["DOMAINS"], "-");
    assert_eq!(
        empty[0].fields["Extra"],
        "IMPACT=MODERATE;DOMAINS=-;SYMBOL=CHODL"
    );
    let freqs_row = |freqs: &str| {
        format!(
            "{header}rs1\t21:18256600\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;FREQS={freqs}\n"
        )
    };
    assert_eq!(
        compare(
            "default",
            &freqs_row("1KG_ALL:G:0.0456,1KG_ALL:T:0.0012"),
            &freqs_row("1KG_ALL:T:0.0012,1KG_ALL:G:0.0456")
        ),
        ""
    );
    let changed_freqs = compare(
        "default",
        &freqs_row("1KG_ALL:G:0.0456,1KG_ALL:T:0.0012"),
        &freqs_row("1KG_ALL:T:0.0012,1KG_ALL:G:0.0457"),
    );
    assert!(
        changed_freqs.contains("FIELD FREQS: 1 mismatching"),
        "{changed_freqs}"
    );

    let vcf = |domains: &str| {
        format!(
            "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|Feature_type|Feature|DOMAINS\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             21\t18256600\t.\tA\tG\t.\t.\tCSQ=G|missense_variant|MODERATE|Transcript|ENST1|{domains}\n"
        )
    };
    assert_eq!(
        compare(
            "vcf",
            &vcf("Gene3D:3.10.100.10&Pfam:PF00059&SMART:SM00034"),
            &vcf("SMART:SM00034&Gene3D:3.10.100.10&Pfam:PF00059")
        ),
        ""
    );
    let changed_vcf = compare(
        "vcf",
        &vcf("Gene3D:3.10.100.10&Pfam:PF00059&SMART:SM00034"),
        &vcf("SMART:SM00034&Gene3D:3.10.100.10&Pfam:PF00060"),
    );
    assert!(
        changed_vcf.contains("FIELD DOMAINS: 1 mismatching"),
        "{changed_vcf}"
    );

    let json = |domains: &str| {
        format!(
            "{{\"input\":\"21\\t18256600\\t.\\tA\\tG\\t.\\t.\\t.\",\"most_severe_consequence\":\"missense_variant\",\
             \"transcript_consequences\":[{{\"variant_allele\":\"G\",\"transcript_id\":\"ENST1\",\"consequence_terms\":[\"missense_variant\"],\"domains\":[{domains}]}}]}}\n"
        )
    };
    let objects = "{\"db\":\"Gene3D\",\"name\":3},{\"db\":\"Pfam\",\"name\":\"PF00059\"},{\"name\":\"SM00034\",\"db\":\"SMART\"}";
    assert_eq!(
        compare(
            "json",
            &json(objects),
            &json("{\"db\":\"SMART\",\"name\":\"SM00034\"},{\"db\":\"Gene3D\",\"name\":3},{\"db\":\"Pfam\",\"name\":\"PF00059\"}")
        ),
        ""
    );
    let changed_json = compare(
        "json",
        &json(objects),
        &json("{\"db\":\"SMART\",\"name\":\"SM00034\"},{\"db\":\"Gene3D\",\"name\":3},{\"db\":\"Pfam\",\"name\":\"PF00060\"}"),
    );
    assert!(
        changed_json.contains("FIELD domains: 1 mismatching"),
        "{changed_json}"
    );
}
