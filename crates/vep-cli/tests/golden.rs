// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Golden-corpus tests: every corpus under `tests/golden/` that holds a
//! `variants.vcf` (a directory `<release>/<name>`), annotated with the
//! corpus's own cache and compared column by column with the Ensembl VEP
//! output committed beside it.
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
    list_columns_as_multisets, masked_reference_text, nearest_as_coordinates, output_path,
    parse_default, parse_json, parse_vcf, parse_vcf_info_in_order, protein_coding_start_sites,
    read_expected, runs, vep_command, Corpus, Documented, Entry, CONSERVATION_SHORT_NAME_SPELLING,
    DEFAULT_SPECIES, ORDER_ONLY, PC_DESCRIPTION_TEXT, REFERENCE_SKIPPED_RECORD,
    REGULATORY_FEATURE_ROWS, UNLOADED_MODULE_PLUGIN, WITH_REFERENCE_SETS,
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
        species: source.species.clone(),
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

/// A manifest's `species` is what every run passes to `--species` and the
/// directory the Ensembl cache layout resolves under `cache_root/`; a manifest
/// that names none runs as `homo_sapiens`, which every such committed corpus is.
#[test]
fn manifest_species_names_the_run_and_the_ensembl_cache_directory() {
    let source = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-run-control")
        .expect("the GRCh38-run-control corpus is committed");
    assert_eq!(source.species, DEFAULT_SPECIES);
    assert!(source.ensembl_cache_layout());
    let tmp = tempfile::tempdir().unwrap();
    let species_arg = |corpus: &Corpus| -> Vec<String> {
        let runs = runs(corpus);
        let argv: Vec<String> = vep_command(&runs[0], tmp.path())
            .get_args()
            .map(|a| a.to_string_lossy().into_owned())
            .collect();
        let at = argv.iter().position(|a| a == "--species").unwrap();
        argv[at..at + 2].to_vec()
    };
    assert_eq!(species_arg(&source), ["--species", DEFAULT_SPECIES]);
    assert_eq!(
        source.cache_dir(),
        source
            .cache_root()
            .join(DEFAULT_SPECIES)
            .join(format!("{}_{}", source.release, source.assembly))
    );
    let mut manifest = source.manifest.clone();
    manifest["species"] = serde_json::Value::from("saccharomyces_cerevisiae");
    manifest["assembly"] = serde_json::Value::from("R64-1-1");
    let yeast = Corpus {
        release: source.release.clone(),
        name: source.name.clone(),
        assembly: "R64-1-1".to_string(),
        species: "saccharomyces_cerevisiae".to_string(),
        dir: source.dir.clone(),
        manifest,
    };
    assert_eq!(
        species_arg(&yeast),
        ["--species", "saccharomyces_cerevisiae"]
    );
    assert_eq!(
        yeast.cache_dir(),
        yeast
            .cache_root()
            .join("saccharomyces_cerevisiae")
            .join(format!("{}_R64-1-1", source.release))
    );
    for corpus in corpora() {
        assert_eq!(
            corpus.species,
            corpus.manifest["species"]
                .as_str()
                .unwrap_or(DEFAULT_SPECIES),
            "{}/{}",
            corpus.release,
            corpus.name
        );
    }
}

/// The VCF output of a run over an input that is not a VCF carries no copy of
/// the input line to match a record by, so its lines are the input records in
/// turn, one line per record, and a documented field divergence keyed on the
/// record applies to the line as it does in the other formats. A run reads a
/// VCF under `--format vcf`, or without a `--format` word when its input
/// carries a VCF header, and matches by the VCF line; a `--format`
/// word naming another format, or a header-less input under no `--format`
/// word (the format Ensembl VEP detects itself), pairs in turn.
#[test]
fn vcf_output_of_a_non_vcf_input_keys_its_lines_by_record() {
    let reads_vcf = |corpus_name: &str, set: &str| -> bool {
        let corpus = corpora()
            .into_iter()
            .find(|c| c.name == corpus_name)
            .unwrap_or_else(|| panic!("the {corpus_name} corpus is committed"));
        runs(&corpus)
            .into_iter()
            .find(|r| r.name == set)
            .unwrap_or_else(|| panic!("flag set {set} of {corpus_name}"))
            .reads_vcf()
    };
    assert!(reads_vcf("GRCh38-regulatory", "default"));
    assert!(reads_vcf("GRCh38-regulatory", "regulatory-shift"));
    assert!(reads_vcf("GRCh38-input-compressed", "stdin"));
    assert!(!reads_vcf("GRCh38-regulatory", "regulatory-ensembl-input"));
    assert!(!reads_vcf("GRCh38-input-ensembl", "ensembl-guess"));
    assert!(!reads_vcf("GRCh38-input-hgvs", "hgvs-guess"));
    let vcf = "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|Feature_type|Feature\">\n\
               #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
               21\t100\t21_100_G/A\tG\tA\t.\t.\tCSQ=A|TF_binding_site_variant|MODIFIER|MotifFeature|ENSM1\n\
               21\t200\trs2\tT\tG\t.\t.\tCSQ=G|missense_variant|MODERATE|Transcript|ENST1\n";
    let keys = |entries: &[Entry]| -> Vec<String> {
        entries.iter().map(|e| e.key.clone()).collect::<Vec<_>>()
    };
    let ensembl_inputs = vec![
        "21\t100\t100\tC/T\t-\t21_100_G/A".to_string(),
        "21\t200\t200\tT/G\t+\trs2".to_string(),
    ];
    assert_eq!(
        keys(&parse_vcf_info_in_order(vcf, &ensembl_inputs, "CSQ")),
        [
            "rec0|record",
            "rec0|A|ENSM1|MotifFeature",
            "rec1|record",
            "rec1|G|ENST1|Transcript"
        ]
    );
    // The same lines against the Ensembl-format records they came from match
    // nothing by VCF line: without the pairing in turn the entries carry no
    // record and a record-keyed divergence cannot reach them.
    assert_eq!(
        keys(&parse_vcf(vcf, &ensembl_inputs)),
        [
            "line0|record",
            "line0|A|ENSM1|MotifFeature",
            "line1|record",
            "line1|G|ENST1|Transcript"
        ]
    );
    // An output with a record missing or one too many is a misalignment the
    // pairing cannot tell from a match, so it is refused by count.
    let three_inputs = [
        &ensembl_inputs[..],
        &["21\t300\t300\tA/G\t+\trs3".to_string()],
    ]
    .concat();
    for inputs in [&ensembl_inputs[..1], &three_inputs[..]] {
        let refused = std::panic::catch_unwind(|| parse_vcf_info_in_order(vcf, inputs, "CSQ"));
        let message = refused
            .expect_err("a line count other than the record count is refused")
            .downcast::<String>()
            .unwrap();
        assert!(
            message.contains(&format!("2 data lines against {} records", inputs.len())),
            "{message}"
        );
    }
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
            species: source.species.clone(),
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
            format!("{label} flag set x: perl_mode is one of offline, database, gtf, gff, not \"cache\""),
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
        ("perl_mode", serde_json::json!(["gtf"])),
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
    let mut corpus = with_flag_sets(serde_json::json!([]));
    corpus.manifest["perl_mode"] = serde_json::json!("cache");
    let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| runs(&corpus).len()));
    assert_eq!(
        outcome.unwrap_err().downcast_ref::<String>().cloned(),
        Some(format!(
            "{label}: perl_mode is one of offline, database, gtf, gff, not \"cache\""
        ))
    );
    let corpus = with_flag_sets(serde_json::json!([
        {"name": "x1_-", "flags": ["--format", "vcf"], "stdin": true, "formats": ["vcf"], "expected_dir": "expected"},
        {"name": "eq", "flags": ["--format=vcf"], "stdin": true},
        {"name": "db", "flags": ["--ga4gh_vrs"], "formats": ["default", "tab", "json"], "perl_mode": "database"},
    ]));
    let corpus_runs = runs(&corpus);
    assert_eq!(corpus_runs.len(), 4);
    assert_eq!(corpus_runs[3].perl_mode, "database");
    assert_eq!(corpus_runs[1].perl_mode, "offline");
    assert!(
        corpus_runs[3].headers_by_descriptions() && !corpus_runs[3].transcript_source_reference()
    );
    // A manifest's own `perl_mode` is its default run's and the flag sets' unless
    // a set names its own.
    let mut corpus = with_flag_sets(serde_json::json!([
        {"name": "inherits"},
        {"name": "own", "perl_mode": "offline"},
    ]));
    corpus.manifest["perl_mode"] = serde_json::json!("gtf");
    let inheriting = runs(&corpus);
    let modes: Vec<&str> = inheriting.iter().map(|r| r.perl_mode.as_str()).collect();
    assert_eq!(modes, ["gtf", "gtf", "offline"]);
    assert!(inheriting[0].transcript_source_reference() && inheriting[0].headers_by_descriptions());
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
            species: source.species.clone(),
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
/// reference's term as written; a run without the mask keeps every row. The
/// mouse corpus's `everything` set declares the class, since its cache holds no
/// regulatory shards; the everything corpus of GRCh38, whose cache holds them,
/// declares none and compares the rows.
#[test]
fn regulatory_feature_rows_mask_drops_the_reference_rows() {
    let mouse = corpora()
        .into_iter()
        .find(|c| c.release == "116" && c.name == "GRCm39-mouse")
        .expect("the mouse corpus");
    let masked = runs(&mouse)
        .into_iter()
        .find(|r| r.name == "everything")
        .expect("the mouse everything set");
    assert_eq!(masked.masks, [REGULATORY_FEATURE_ROWS]);
    assert!(masked.masks_regulatory_rows());
    let unmasked = corpora()
        .into_iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-everything")
        .expect("the everything corpus");
    assert!(unmasked.masks().is_empty());
    assert!(!common::default_run(&unmasked).masks_regulatory_rows());

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
    apply_masks(&masked.masks, "default", &mut entries);
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
    apply_masks(&masked.masks, "vcf", &mut entries);
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
        entries.iter().map(|e| e.key.as_str()).collect::<Vec<_>>(),
        [
            "rec0|record",
            "rec0|A|ENSR00000001|regulatory",
            "rec0|A||intergenic",
            "rec1|record",
            "rec1|G|ENST1|transcript"
        ]
    );
    assert!(!entries[0]
        .fields
        .contains_key("regulatory_feature_consequences"));
    assert_eq!(
        entries[1].fields["consequence_terms"],
        "regulatory_region_variant"
    );
    apply_masks(&masked.masks, "json", &mut entries);
    assert_eq!(
        entries.iter().map(|e| e.key.as_str()).collect::<Vec<_>>(),
        [
            "rec0|record",
            "rec0|A||intergenic",
            "rec1|record",
            "rec1|G|ENST1|transcript"
        ]
    );
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
}

/// A documented field value is read in the form each format writes: the
/// manifest documents the default format's `1311092.0000`, the JSON format
/// writes `1311092`, and the two are one value; another number, or another
/// spelling of the same number, is stale.
#[test]
fn documented_numeric_field_value_compares_across_formats() {
    let documented = Documented::from_manifest(&serde_json::json!({
        "field_divergences": [{
            "location": "21:36170554", "allele": "CA", "feature": "ENST1",
            "feature_type": "Transcript", "record_indices": [0],
            "consequence_set": "intron_variant", "field": "gnomAD_total_DP_cov",
            "vep": "1311160.0000", "vep_rs": "1311092.0000",
            "expected_divergence": "minimal_plugin_sees_the_rejoined_record",
        }]
    }));
    let inputs = vec!["21\t36170554\trs1\tC\tCCA,T\t.\t.\t.".to_string()];
    let object = |total: &str| {
        format!(
            "{{\"input\":\"21\\t36170554\\trs1\\tC\\tCCA,T\\t.\\t.\\t.\",\"id\":\"rs1\",\
             \"transcript_consequences\":[{{\"variant_allele\":\"CA\",\"transcript_id\":\"ENST1\",\
             \"consequence_terms\":[\"intron_variant\"],\"gnomad_total_dp_cov\":{total}}}]}}\n"
        )
    };
    let compare = |actual: &str| {
        compare_entries(
            &parse_json(&object("1311160"), &inputs),
            &parse_json(&object(actual), &inputs),
            &documented,
            true,
        )
    };
    assert_eq!(compare("1311092"), "");
    let stale = compare("1311093");
    assert!(stale.contains("documented as"), "{stale}");
    // The same number spelt another way is not the documented value.
    for respelt in ["1311092.0", "1.311092e6"] {
        let stale = compare(respelt);
        assert!(stale.contains("documented as"), "{stale}");
    }
    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
    let row = |total: &str| {
        format!(
            "{header}rs1\t21:36170554\tCA\tENSG1\tENST1\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;gnomAD_total_DP_cov={total}\n"
        )
    };
    assert_eq!(
        compare_entries(
            &parse_default(&row("1311160.0000")),
            &parse_default(&row("1311092.0000")),
            &documented,
            true
        ),
        ""
    );
}

/// A record-level row a `--summary` run documents, which names no feature and
/// carries the terms of every feature, is found by its Location in the default
/// format, where a row printing another set is stale. The JSON format's record
/// object of the same record prints no consequence set, so the entry has
/// nothing to observe on it.
#[test]
fn documented_record_level_row_has_nothing_to_observe_on_a_json_record_object() {
    let documented = Documented::from_manifest(&serde_json::json!({
        "divergences": [{
            "location": "21:100-104", "allele": "-", "feature": "-", "feature_type": "-",
            "record_indices": [0],
            "vep_consequence_set": "intron_variant,stop_retained_variant",
            "vep_rs_consequence_sets": ["coding_sequence_variant,intron_variant"],
            "expected_divergence": "stop_retained_without_stop_codon_test",
        }]
    }));
    let inputs = vec!["21\t99\trs1\tGACCTG\tG\t.\t.\t.".to_string()];
    let object = |most_severe: &str| {
        format!(
            "{{\"input\":\"21\\t99\\trs1\\tGACCTG\\tG\\t.\\t.\\t.\",\"id\":\"rs1\",\
             \"most_severe_consequence\":\"{most_severe}\"}}\n"
        )
    };
    assert_eq!(
        compare_entries(
            &parse_json(&object("stop_retained_variant"), &inputs),
            &parse_json(&object("coding_sequence_variant"), &inputs),
            &documented,
            true
        ),
        ""
    );
    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
    let row = |terms: &str| {
        format!("{header}rs1\t21:100-104\t-\t-\t-\t-\t{terms}\t-\t-\t-\t-\t-\t-\tIMPACT=LOW\n")
    };
    let compare = |actual: &str| {
        compare_entries(
            &parse_default(&row("intron_variant,stop_retained_variant")),
            &parse_default(&row(actual)),
            &documented,
            true,
        )
    };
    assert_eq!(compare("intron_variant,coding_sequence_variant"), "");
    let stale = compare("intron_variant,frameshift_variant");
    assert!(stale.contains("vep-rs now emits"), "{stale}");
}

/// A `divergences` item that states vep-rs's sets for a run with a reference
/// sequence is read for the run's setting. A run without one is held to the
/// item's `vep_rs_consequence_sets`. In a run with one, stated sets that are
/// VEP's own leave the key compared like any other: a row that agrees with
/// VEP's passes, and one carrying the sets of a run without a reference
/// sequence is a difference of its Consequence. Each crossed run fails, stated
/// sets other than VEP's are the ones such a run is held to, and an item
/// without the key reads the same under both settings.
#[test]
fn documented_item_is_read_for_the_runs_reference_setting() {
    let kept = "3_prime_UTR_variant,stop_retained_variant";
    let lost = "3_prime_UTR_variant,stop_lost";
    let other = "3_prime_UTR_variant,coding_sequence_variant";
    let manifest = |with_reference: Option<&str>| {
        let mut item = serde_json::json!({
            "location": "21:100-104", "allele": "-", "feature": "ENST1",
            "feature_type": "Transcript", "record_indices": [0],
            "vep_consequence_set": kept,
            "vep_rs_consequence_sets": [lost],
            "expected_divergence": "stop_codon_removed_without_stop_codon_test",
        });
        if let Some(sets) = with_reference {
            item[WITH_REFERENCE_SETS] = serde_json::json!([sets]);
        }
        serde_json::json!({ "divergences": [item] })
    };
    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
    let row = |terms: &str| {
        format!(
            "{header}rs1\t21:100-104\t-\tENSG1\tENST1\tTranscript\t{terms}\t-\t-\t-\t-\t-\t-\tIMPACT=LOW\n"
        )
    };
    let compare = |documented: &Documented, actual: &str| {
        compare_entries(
            &parse_default(&row(kept)),
            &parse_default(&row(actual)),
            documented,
            true,
        )
    };
    let stated = manifest(Some(kept));
    for without_reference in [
        Documented::for_run(&stated, false),
        Documented::from_manifest(&stated),
    ] {
        assert_eq!(compare(&without_reference, lost), "");
        let crossed = compare(&without_reference, kept);
        assert!(crossed.contains("NOT OBSERVED IN THIS RUN: 1"), "{crossed}");
    }
    let with_reference = Documented::for_run(&stated, true);
    assert_eq!(compare(&with_reference, kept), "");
    let crossed = compare(&with_reference, lost);
    assert!(
        crossed.contains("FIELD Consequence: 1 mismatching entries")
            && !crossed.contains("NOT OBSERVED"),
        "{crossed}"
    );
    let stated_other = Documented::for_run(&manifest(Some(other)), true);
    assert_eq!(compare(&stated_other, other), "");
    let crossed = compare(&stated_other, lost);
    assert!(crossed.contains("NOT OBSERVED IN THIS RUN: 1"), "{crossed}");
    let plain = manifest(None);
    for setting in [false, true] {
        let documented = Documented::for_run(&plain, setting);
        assert_eq!(compare(&documented, lost), "");
        let stale = compare(&documented, kept);
        assert!(stale.contains("NOT OBSERVED IN THIS RUN: 1"), "{stale}");
    }
}

/// A `vep_rs_consequence_sets_with_reference_sequence` key of another shape
/// than a sorted list of distinct sets fails the run that reads it, under
/// either setting.
#[test]
#[should_panic(
    expected = "the item on 21:100-104 - ENST1: vep_rs_consequence_sets_with_reference_sequence names VEP's set beside another"
)]
fn documented_item_with_a_malformed_reference_key_is_refused() {
    Documented::for_run(
        &serde_json::json!({ "divergences": [{
            "location": "21:100-104", "allele": "-", "feature": "ENST1",
            "feature_type": "Transcript", "record_indices": [0],
            "vep_consequence_set": "3_prime_UTR_variant,stop_retained_variant",
            "vep_rs_consequence_sets": ["3_prime_UTR_variant,stop_lost"],
            "expected_divergence": "stop_codon_removed_without_stop_codon_test",
            "vep_rs_consequence_sets_with_reference_sequence": [
                "3_prime_UTR_variant,coding_sequence_variant",
                "3_prime_UTR_variant,stop_retained_variant",
            ],
        }] }),
        false,
    );
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

    // The UTRAnnotator plugin's `5UTR_annotation`: the pairs of a group and
    // the groups compare in any order; a changed value and a dropped group
    // still differ. In the VCF the groups are `&`-joined.
    let utr_row = |annotation: &str| {
        format!(
            "{header}rs1\t21:18256600\tG\tENSG1\tENST1\tTranscript\t5_prime_UTR_variant\t-\t-\t-\t-\t-\t-\t5UTR_annotation={annotation};5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;IMPACT=MODIFIER\n"
        )
    };
    let fixed = "Evidence=NA:KozakContext=CCCATGC:alt_type=uORF,Evidence=False:KozakContext=-:alt_type=InFrame_oORF";
    assert_eq!(
        compare(
            "default",
            &utr_row("alt_type=InFrame_oORF:KozakContext=-:Evidence=False,KozakContext=CCCATGC:Evidence=NA:alt_type=uORF"),
            &utr_row(fixed)
        ),
        ""
    );
    let changed_utr = compare(
        "default",
        &utr_row("alt_type=InFrame_oORF:KozakContext=-:Evidence=True,KozakContext=CCCATGC:Evidence=NA:alt_type=uORF"),
        &utr_row(fixed),
    );
    assert!(
        changed_utr.contains("FIELD 5UTR_annotation: 1 mismatching"),
        "{changed_utr}"
    );
    let fewer_utr = compare(
        "default",
        &utr_row("KozakContext=CCCATGC:Evidence=NA:alt_type=uORF"),
        &utr_row(fixed),
    );
    assert!(
        fewer_utr.contains("FIELD 5UTR_annotation: 1 mismatching"),
        "{fewer_utr}"
    );
    let utr_vcf = |annotation: &str| {
        format!(
            "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|Feature_type|Feature|5UTR_annotation\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             21\t18256600\t.\tA\tG\t.\t.\tCSQ=G|5_prime_UTR_variant|MODIFIER|Transcript|ENST1|{annotation}\n"
        )
    };
    assert_eq!(
        compare(
            "vcf",
            &utr_vcf("alt_type=InFrame_oORF:KozakContext=-&KozakContext=CCCATGC:alt_type=uORF"),
            &utr_vcf("KozakContext=CCCATGC:alt_type=uORF&KozakContext=-:alt_type=InFrame_oORF")
        ),
        ""
    );
}

/// Under [`UNLOADED_MODULE_PLUGIN`] the reference's run-time warnings about the
/// plugin module (`did not return a hashref`, `went wrong`) are set aside and
/// vep-rs's one `Failed to load plugin <Name>` line stands in for them: a
/// warnings file with that line alone agrees, one without it is reported, and
/// any other line is still a line the reference lacks. The faulty-plugin corpus
/// declares the class for its own run.
#[test]
fn unloaded_module_plugin_mask_stands_the_load_warning_in_for_the_reference() {
    let corpus = corpora()
        .into_iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-plugin-faulty")
        .expect("the faulty-plugin corpus");
    assert_eq!(corpus.masks(), [UNLOADED_MODULE_PLUGIN]);
    let run = common::default_run(&corpus);
    assert!(run.masks_unloaded_module_plugin() && !run.masks_regulatory_rows());
    let tmp = tempfile::tempdir().unwrap();
    let warnings = format!(
        "{}_warnings.txt",
        output_path(&run, "default", tmp.path()).display()
    );
    let load_failure = "WARNING: Failed to load plugin Faulty: Faulty is not a built-in plugin and no libFaulty.so or libFaulty.dylib is under /somewhere/plugin_data\n";

    fs::write(&warnings, load_failure).unwrap();
    assert_eq!(compare_warnings(&run, tmp.path()), "");

    fs::write(&warnings, "").unwrap();
    let report = compare_warnings(&run, tmp.path());
    assert!(
        report.contains("0 lines of") && report.contains("plugin module Faulty failed to load"),
        "{report}"
    );

    fs::write(
        &warnings,
        format!("{load_failure}WARNING: Plugin 'Faulty' went wrong: Faulty on purpose\n"),
    )
    .unwrap();
    let report = compare_warnings(&run, tmp.path());
    assert!(
        report.contains("warning not in the reference")
            && report.contains("Plugin 'Faulty' went wrong"),
        "{report}"
    );
}

/// Under [`CONSERVATION_SHORT_NAME_SPELLING`] the reference text is read with
/// every `phlyoP<N>way` key spelt `phyloP<N>way` (Ensembl VEP's `--phyloP` set
/// misspells the short name of the source it adds, `Config.pm:555`, while the
/// file it names and vep-rs's key are spelt `phyloP<N>way`): the header line,
/// the tab column, the `Extra` item, the CSQ field list, the `##INFO` line and
/// the JSON `custom_annotations` member alike, a key of any width, while the
/// file's name, a `phlyoP` outside a key and every other text stay as written.
/// The custom corpus's `phylop-phastcons` set declares the class.
#[test]
fn conservation_short_name_spelling_mask_respells_the_reference_key() {
    let corpus = corpora()
        .into_iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-custom")
        .expect("the custom corpus");
    let run = runs(&corpus)
        .into_iter()
        .find(|r| r.name == "phylop-phastcons")
        .expect("the phylop-phastcons set");
    assert_eq!(run.masks, [CONSERVATION_SHORT_NAME_SPELLING]);
    for name in ["default.txt", "tab.txt", "vcf.vcf", "json.jsonl"] {
        let reference = run.expected(name);
        assert!(reference.contains("phlyoP100way"), "{name}");
        let masked = masked_reference_text(&run.masks, reference.clone());
        assert!(!masked.contains("phlyoP"), "{name}");
        assert_eq!(
            masked,
            reference.replace("phlyoP100way", "phyloP100way"),
            "{name}"
        );
        assert_eq!(
            masked_reference_text(&[], reference.clone()),
            reference,
            "{name}"
        );
    }
    let masks = [CONSERVATION_SHORT_NAME_SPELLING.to_string()];
    let text = "\
## phlyoP100way : [PATH]/hg38.phyloP100way.bw
#Uploaded_variation\tLocation\tphlyoP100way\tphlyoP20way
rs1\t21:100\tIMPACT=MODIFIER;phlyoP100way=-0.98;phastCons100way=0.998\t-
##INFO=<ID=phlyoP470way,Number=.,Type=String,Description=\"[PATH]/hg38.phyloP470way.bw\">
{\"custom_annotations\":{\"phlyoP100way\":[{\"name\":-0.98}]}}
phlyoPway phlyoP phlyoP1 phlyoP1wa
";
    let expected = "\
## phyloP100way : [PATH]/hg38.phyloP100way.bw
#Uploaded_variation\tLocation\tphyloP100way\tphyloP20way
rs1\t21:100\tIMPACT=MODIFIER;phyloP100way=-0.98;phastCons100way=0.998\t-
##INFO=<ID=phyloP470way,Number=.,Type=String,Description=\"[PATH]/hg38.phyloP470way.bw\">
{\"custom_annotations\":{\"phyloP100way\":[{\"name\":-0.98}]}}
phlyoPway phlyoP phlyoP1 phlyoP1wa
";
    assert_eq!(masked_reference_text(&masks, text.to_string()), expected);
    assert_eq!(
        masked_reference_text(&[PC_DESCRIPTION_TEXT.to_string()], text.to_string()),
        text
    );
}

/// Under [`PC_DESCRIPTION_TEXT`] the reference's header line for a custom
/// source's `PC` field carries vep-rs's description of the share the field
/// reports in place of Ensembl VEP's: `AnnotationSourceAdaptor.pm:395-397`
/// describes `reciprocal=0` as the share of the input variant covered and
/// `reciprocal=1` as a minimum reciprocal overlap, while `File.pm:556-564`
/// computes the share of the annotation record covered without `reciprocal`
/// (`$ref_overlap_percentage`) and the share of the input variant covered with
/// it (`$overlap_percentage` over `$length`, the variant's). Both texts are
/// mapped, in the `## <key>_PC : ` shape and the `##INFO=<ID=<key>_PC` shape;
/// the same text on a line describing another key, a row carrying the field's
/// values and every other header line stay as written. The custom corpus's
/// `bed-overlap-stats` set (a `reciprocal=1` source) declares the class.
#[test]
fn pc_description_text_mask_maps_both_overlap_texts() {
    let corpus = corpora()
        .into_iter()
        .find(|c| c.release == "116" && c.name == "GRCh38-custom")
        .expect("the custom corpus");
    let run = runs(&corpus)
        .into_iter()
        .find(|r| r.name == "bed-overlap-stats")
        .expect("the bed-overlap-stats set");
    assert_eq!(run.masks, [PC_DESCRIPTION_TEXT]);
    let reciprocal =
        "Percentage of minimum reciprocal overlap between input variant and reference variant";
    let plain = "Percentage of input variant covered by reference variant";
    let reference = run.expected("default.txt");
    let masked = masked_reference_text(&run.masks, reference.clone());
    assert_eq!(
        masked,
        reference.replace(
            &format!("## regions_PC : {reciprocal} from"),
            "## regions_PC : Percentage of the input variant covered by the annotation record from"
        )
    );
    assert_ne!(masked, reference);
    let reference = run.expected("vcf.vcf");
    let masked = masked_reference_text(&run.masks, reference.clone());
    assert_eq!(
        masked,
        reference.replace(
            &format!("Description=\"{reciprocal} from"),
            "Description=\"Percentage of the input variant covered by the annotation record from"
        )
    );
    assert_ne!(masked, reference);

    let text = format!(
        "## regions_PC : {reciprocal} from [PATH]/regions.bed.gz
## plain_PC : {plain} from [PATH]/plain.bed.gz
##INFO=<ID=regions_PC,Number=.,Type=String,Description=\"{reciprocal} from [PATH]/regions.bed.gz\">
##INFO=<ID=plain_PC,Number=.,Type=String,Description=\"{plain} from [PATH]/plain.bed.gz\">
## regions_min : min data from [PATH]/regions.bed.gz
## other : {plain}
rs1\t21:5028998-5028999\tregions_PC=100.000,50.000;{plain}
"
    );
    let expected = format!(
        "## regions_PC : Percentage of the input variant covered by the annotation record from [PATH]/regions.bed.gz
## plain_PC : Percentage of the annotation record covered by the input variant from [PATH]/plain.bed.gz
##INFO=<ID=regions_PC,Number=.,Type=String,Description=\"Percentage of the input variant covered by the annotation record from [PATH]/regions.bed.gz\">
##INFO=<ID=plain_PC,Number=.,Type=String,Description=\"Percentage of the annotation record covered by the input variant from [PATH]/plain.bed.gz\">
## regions_min : min data from [PATH]/regions.bed.gz
## other : {plain}
rs1\t21:5028998-5028999\tregions_PC=100.000,50.000;{plain}
"
    );
    assert_eq!(
        masked_reference_text(&[PC_DESCRIPTION_TEXT.to_string()], text.clone()),
        expected
    );
    assert_eq!(masked_reference_text(&[], text.clone()), text);
    assert_eq!(
        masked_reference_text(
            &[CONSERVATION_SHORT_NAME_SPELLING.to_string()],
            text.clone()
        ),
        text
    );
}

/// A GTF or GFF reference loses its transcript source before the comparison:
/// the `SOURCE` column Ensembl VEP fills with the source's short name on every
/// transcript it reads (`AnnotationSource/File/BaseGXF.pm:635` sets
/// `_source_cache`, `OutputFactory.pm:1583` prints it, `Constants.pm:101` lists
/// the column once a custom source is loaded), the source's own annotation
/// column (empty on every row, since a transcript source annotates no record)
/// and the header lines describing both, in each format; every other column,
/// header line and value stays as written.
#[test]
fn transcript_source_reference_loses_its_source_column_and_header_lines() {
    let default_text = "\
## ENSEMBL VARIANT EFFECT PREDICTOR v116.2
## Using API version 116, DB version ?
## Extra column keys:
## IMPACT : Subjective impact classification of consequence type
## SYMBOL : Gene symbol (e.g. HGNC)
## SOURCE : Source of transcript
## genes.gtf.gz : [PATH]/genes.gtf.gz
## VEP command-line: vep --gtf [PATH]/genes.gtf.gz
#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra
rs1\t21:100\tG\tENSG1\tENST1\tTranscript\tmissense_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODERATE;SYMBOL=PWP2;SOURCE=genes.gtf.gz
rs2\t21:200\tA\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tSOURCE=genes.gtf.gz
";
    let default_masked = "\
## ENSEMBL VARIANT EFFECT PREDICTOR v116.2
## Using API version 116, DB version ?
## Extra column keys:
## IMPACT : Subjective impact classification of consequence type
## SYMBOL : Gene symbol (e.g. HGNC)
## VEP command-line: vep --gtf [PATH]/genes.gtf.gz
#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra
rs1\t21:100\tG\tENSG1\tENST1\tTranscript\tmissense_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODERATE;SYMBOL=PWP2
rs2\t21:200\tA\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t-
";
    assert_eq!(
        common::without_transcript_source("default", default_text, "CSQ"),
        default_masked
    );

    let tab_text = "\
## SOURCE : Source of transcript
## genes.gtf.gz : [PATH]/genes.gtf.gz
#Uploaded_variation\tLocation\tAllele\tIMPACT\tSOURCE\tgenes.gtf.gz\tSYMBOL
rs1\t21:100\tG\tMODERATE\tgenes.gtf.gz\t-\tPWP2
";
    assert_eq!(
        common::without_transcript_source("tab", tab_text, "CSQ"),
        "#Uploaded_variation\tLocation\tAllele\tIMPACT\tSYMBOL\nrs1\t21:100\tG\tMODERATE\tPWP2\n"
    );

    let vcf_text = "\
##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|SOURCE|SYMBOL|genes.gtf.gz\">
##INFO=<ID=genes.gtf.gz,Number=.,Type=String,Description=\"[PATH]/genes.gtf.gz\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
21\t100\trs1\tA\tG\t.\t.\tDP=5;CSQ=G|missense_variant|genes.gtf.gz|PWP2|,G|intron_variant|genes.gtf.gz|PWP2|
";
    let vcf_masked = "\
##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|SYMBOL\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
21\t100\trs1\tA\tG\t.\t.\tDP=5;CSQ=G|missense_variant|PWP2,G|intron_variant|PWP2
";
    assert_eq!(
        common::without_transcript_source("vcf", vcf_text, "CSQ"),
        vcf_masked
    );
    // Under --keep_csq the input's own CSQ header line, which names no format,
    // and its CSQ entry, the line's first, stay as written; the fields dropped
    // are the last header line's and the entry rewritten the line's last.
    let kept_csq_text = "\
##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"An earlier annotation\">
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|SOURCE|SYMBOL|genes.gtf.gz\">
##INFO=<ID=genes.gtf.gz,Number=.,Type=String,Description=\"[PATH]/genes.gtf.gz\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
21\t100\trs1\tA\tG\t.\t.\tCSQ=G|earlier|x;DP=5;CSQ=G|missense_variant|genes.gtf.gz|PWP2|
";
    let kept_csq_masked = "\
##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"An earlier annotation\">
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|SYMBOL\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
21\t100\trs1\tA\tG\t.\t.\tCSQ=G|earlier|x;DP=5;CSQ=G|missense_variant|PWP2
";
    assert_eq!(
        common::without_transcript_source("vcf", kept_csq_text, "CSQ"),
        kept_csq_masked
    );

    let json_text = "{\"input\":\"21\\t100\\trs1\\tA\\tG\",\"transcript_consequences\":[{\"transcript_id\":\"ENST1\",\"source\":\"genes.gtf.gz\",\"impact\":\"MODERATE\"}],\"intergenic_consequences\":[{\"impact\":\"MODIFIER\"}]}\n";
    let masked: serde_json::Value =
        serde_json::from_str(common::without_transcript_source("json", json_text, "CSQ").trim())
            .unwrap();
    assert_eq!(
        masked,
        serde_json::json!({
            "input": "21\t100\trs1\tA\tG",
            "transcript_consequences": [{"transcript_id": "ENST1", "impact": "MODERATE"}],
            "intergenic_consequences": [{"impact": "MODIFIER"}]
        })
    );
}

/// The custom-annotation corpus and the GTF-route corpus keep their size: every
/// run's reference holds at least the rows the committed references hold, the
/// corpora their records, and each custom source annotates at least the rows
/// it annotates in the committed references (the default-format rows carrying
/// the source's column), so a regenerated fixture that matched nothing, or a
/// cache built from a slice missing a gene, fails here before the comparison.
#[test]
fn custom_and_transcript_source_corpora_keep_their_floors() {
    /// A run's floors: the reference's rows and, for a custom source, the
    /// rows carrying its column.
    struct Floor(
        &'static str,
        &'static str,
        usize,
        Option<(&'static str, usize)>,
    );
    let floors = [
        Floor("GRCh38-custom", "default", 1586, None),
        Floor("GRCh38-custom", "vcf-exact", 1586, Some(("sites", 454))),
        Floor(
            "GRCh38-custom",
            "vcf-exact-positional",
            1586,
            Some(("sites", 454)),
        ),
        Floor(
            "GRCh38-custom",
            "vcf-overlap-coords",
            1586,
            Some(("sites", 944)),
        ),
        Floor(
            "GRCh38-custom",
            "vcf-multi-allelic",
            1586,
            Some(("sites", 454)),
        ),
        Floor(
            "GRCh38-custom",
            "vcf-clinvar-source",
            1586,
            Some(("ClinVar", 454)),
        ),
        Floor(
            "GRCh38-custom",
            "vcf-suppress-filter",
            1586,
            Some(("sites", 454)),
        ),
        Floor(
            "GRCh38-custom",
            "bed-overlap-stats",
            1586,
            Some(("regions", 944)),
        ),
        Floor(
            "GRCh38-custom",
            "bed-within-distance",
            1586,
            Some(("regions", 601)),
        ),
        Floor(
            "GRCh38-custom",
            "bed-chr-names",
            1586,
            Some(("regions", 944)),
        ),
        Floor(
            "GRCh38-custom",
            "gff-promoter",
            1586,
            Some(("GENCODE_PROMOTER", 944)),
        ),
        Floor(
            "GRCh38-custom",
            "bigwig-shortcut",
            1586,
            Some(("scores.bw", 197)),
        ),
        Floor(
            "GRCh38-custom",
            "bigwig-overlap-stats",
            1586,
            Some(("cons", 481)),
        ),
        Floor(
            "GRCh38-custom",
            "phastcons",
            1586,
            Some(("phastCons100way", 197)),
        ),
        Floor("GRCh38-gtf-builder", "default", 862, None),
        Floor("GRCh38-gtf-builder", "filter-biotype", 714, None),
    ];
    let records: [(&str, usize); 2] = [("GRCh38-custom", 73), ("GRCh38-gtf-builder", 74)];
    let corpora = corpora();
    let corpus_named = |name: &str| {
        corpora
            .iter()
            .find(|c| c.release == "116" && c.name == name)
            .unwrap_or_else(|| panic!("the {name} corpus"))
    };
    for (name, floor) in records {
        assert!(
            common::input_records(corpus_named(name)).len() >= floor,
            "{name}: fewer than {floor} records"
        );
    }
    for Floor(name, run_name, rows, column) in floors {
        let corpus = corpus_named(name);
        let run = runs(corpus)
            .into_iter()
            .find(|r| r.name == run_name)
            .unwrap_or_else(|| panic!("{name}: the run {run_name}"));
        let entries = parse_default(&run.expected("default.txt"));
        assert!(
            entries.len() >= rows,
            "{}: {} rows in the reference, fewer than {rows}",
            run.label(),
            entries.len()
        );
        if let Some((key, populated)) = column {
            let carrying = entries
                .iter()
                .filter(|e| e.fields.contains_key(key))
                .count();
            assert!(
                carrying >= populated,
                "{}: {carrying} rows carry {key}, fewer than {populated}",
                run.label()
            );
        }
    }
}
