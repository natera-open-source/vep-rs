// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `select` on small reference outputs and their inputs: the corpus and the manifest
//! it writes, the exemplars the options choose, the keys a regeneration carries, and the
//! runs that are refused before anything is written, a regeneration that would discard a
//! hand-written statement among them.

mod common;

use common::{
    at, clean, failed, key_place, lines_of, manifest, manifest_file, manifest_text, matrix,
    one_stratum_dataset, record_by_id, record_ids, records, reference_statement, row, row_with,
    run, scratch, select, text_of, two_record_dataset, under, vcf_of, write, write_manifest_line,
    Run, VCF_HEADER,
};
use std::collections::BTreeSet;
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::Path;
use vep_golden_corpus::cache_cost::CacheCost;
use vep_golden_corpus::classify::cache::corpus_cache_dir;
use vep_golden_corpus::constants::{mask_classes, WITH_REFERENCE_SETS};
use vep_golden_corpus::json::{loads, Object, Value};
use vep_golden_corpus::manifest::{generated_shape_error, manifest_shape_error};
use vep_golden_corpus::rows::{iter_rows, Row};
use vep_golden_corpus::sampling::{choose_exemplars, OrderedMap, Random, Sample};

/// The summary line of a run that covered everything it observed.
fn summary(combinations: usize, records: usize, classes: &str) -> String {
    let count = classes.split(", ").count();
    format!(
        "combinations observed {combinations}, covered {combinations}, records {records}, unresolved exemplars 0, variant classes {count} ({classes})\n"
    )
}

fn exists(path: &str) -> bool {
    Path::new(path).exists()
}

#[test]
fn select_builds_corpus_and_manifest() {
    let scratch = scratch("select-builds");
    let rows = row("rs1", "21:100", "G", "ENST1", "missense_variant")
        + &row("rs1", "21:100", "G", "ENST2", "intron_variant")
        + &row_with(
            "21_201_AT/A",
            "21:201",
            "-",
            "ENST1",
            "frameshift_variant",
            "IMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF",
        )
        + &row("21_300_A/G/T", "21:300", "G", "ENST3", "synonymous_variant")
        + &row("21_300_A/G/T", "21:300", "T", "ENST3", "missense_variant");
    let reference = write(&scratch.path, "ref.txt", &rows);
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!(
            "{VCF_HEADER}21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t150\trs9\tC\tT\t.\tPASS\t.\n21\t200\t.\tAT\tA\t.\tPASS\t.\n21\t300\t.\tA\tG,T\t.\tPASS\t.\n21\t400\trs4\tCA\tC,CGG\t.\tPASS\t.\n"
        ),
    );
    let dataset = format!("s01={reference}={input}");
    let out = under(&scratch.path, "corpus");
    let select_into = |out: &str, multi: &str| {
        run(&[
            "select",
            "--release",
            "115",
            "--assembly",
            "GRCh37",
            "--dataset",
            &dataset,
            "--out",
            out,
            "--k",
            "2",
            "--seed",
            "3",
            "--multi-allelic-per-dataset",
            multi,
            "--covers",
            "symbol",
        ])
    };
    assert_eq!(
        select_into(&out, "0"),
        clean(0, &summary(4, 3, "deletion 1, snv 4"))
    );
    let m = manifest(&out);
    let combinations: BTreeSet<&str> = at(&m, &["combinations"])
        .as_object()
        .unwrap()
        .keys()
        .collect();
    assert_eq!(
        combinations,
        BTreeSet::from([
            "missense_variant",
            "intron_variant",
            "frameshift_variant",
            "synonymous_variant"
        ])
    );
    assert_eq!(text_of(at(&m, &["uncovered_combinations"])), "[]");
    assert_eq!(text_of(at(&m, &["unresolved_exemplars"])), "[]");
    let lines: Vec<&str> = records(&m).iter().map(vcf_of).collect();
    assert_eq!(
        lines,
        [
            "21\t100\trs1\tA\tG\t.\tPASS\t.",
            "21\t200\t.\tAT\tA\t.\tPASS\t.",
            "21\t300\t.\tA\tG,T\t.\tPASS\t."
        ]
    );
    let variants = fs::read_to_string(under(Path::new(&out), "variants.vcf")).unwrap();
    let body = lines_of(&variants);
    assert_eq!(body[0], "##fileformat=VCFv4.2");
    assert!(body[1].starts_with("#CHROM") && body.len() == 5);
    // The record no row selected is not carried.
    assert!(!variants.contains("rs9"));

    // With multi-allelic sampling on, the SNV-only and the indel multi-allelic records
    // both join the corpus, tagged by stratum, without claiming combinations.
    let out2 = under(&scratch.path, "corpus2");
    assert_eq!(
        select_into(&out2, "2"),
        clean(0, &summary(4, 4, "deletion 1, snv 4"))
    );
    let m2 = manifest(&out2);
    let by_vcf = |line: &str| {
        let found = records(&m2).iter().find(|record| vcf_of(record) == line);
        found.expect("the record").clone()
    };
    let snv = by_vcf("21\t300\t.\tA\tG,T\t.\tPASS\t.");
    assert_eq!(
        text_of(snv.get("strata").unwrap()),
        r#"["multi_allelic_snv"]"#
    );
    let indel = by_vcf("21\t400\trs4\tCA\tC,CGG\t.\tPASS\t.");
    assert_eq!(
        text_of(indel.get("strata").unwrap()),
        r#"["multi_allelic_indel"]"#
    );
    assert_eq!(text_of(indel.get("combinations").unwrap()), "[]");
}

#[test]
fn select_writes_the_species_and_the_cache_dir_follows_it() {
    let scratch = scratch("select-species");
    let reference = write(
        &scratch.path,
        "ref.txt",
        &row("rs1", "IV:100", "G", "YDL001", "missense_variant"),
    );
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!("{VCF_HEADER}IV\t100\trs1\tA\tG\t.\tPASS\t.\n"),
    );
    let dataset = format!("s01={reference}={input}");
    let common = [
        "--dataset",
        &dataset,
        "--k",
        "1",
        "--multi-allelic-per-dataset",
        "0",
        "--covers",
        "symbol",
    ];
    let select_with = |first: &[&str]| {
        let mut args = vec!["select"];
        args.extend(first);
        args.extend(common);
        run(&args)
    };
    let one = summary(1, 1, "snv 1");
    let yeast = under(&scratch.path, "yeast");
    assert_eq!(
        select_with(&[
            "--release",
            "116",
            "--assembly",
            "R64-1-1",
            "--species",
            "saccharomyces_cerevisiae",
            "--out",
            &yeast
        ]),
        clean(0, &one)
    );
    let mut m = manifest(&yeast);
    assert_eq!(
        at(&m, &["species"]).as_str(),
        Some("saccharomyces_cerevisiae")
    );
    assert_eq!(at(&m, &["assembly"]).as_str(), Some("R64-1-1"));
    assert_eq!(corpus_cache_dir(&yeast, &m), format!("{yeast}/json_cache"));
    m.insert("cache_layout", "ensembl");
    assert_eq!(
        corpus_cache_dir(&yeast, &m),
        format!("{yeast}/cache_root/saccharomyces_cerevisiae/116_R64-1-1")
    );
    let human = under(&scratch.path, "human");
    assert_eq!(
        select_with(&["--release", "116", "--assembly", "GRCh38", "--out", &human]),
        clean(0, &one)
    );
    let mut m = manifest(&human);
    assert_eq!(at(&m, &["species"]).as_str(), Some("homo_sapiens"));
    m.remove("species");
    m.insert("cache_layout", "ensembl");
    assert_eq!(
        corpus_cache_dir(&human, &m),
        format!("{human}/cache_root/homo_sapiens/116_GRCh38")
    );
    assert_eq!(manifest_shape_error(&m), None);
    m.insert("species", "Saccharomyces cerevisiae");
    assert_eq!(
        manifest_shape_error(&m).as_deref(),
        Some("species is \"Saccharomyces cerevisiae\", not an Ensembl species name")
    );
    m.insert("species", 7i64);
    assert_eq!(
        manifest_shape_error(&m).as_deref(),
        Some("species is 7, not an Ensembl species name")
    );
    let bad = under(&scratch.path, "bad");
    assert_eq!(
        select_with(&[
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--species",
            "Homo sapiens",
            "--out",
            &bad
        ]),
        failed(
            "--species 'Homo sapiens' is not an Ensembl species name (lower-case letters, digits and underscores)"
        )
    );
    assert!(!exists(&bad));
    // A regeneration over the yeast manifest keeps its species without --species.
    assert_eq!(
        select_with(&["--release", "116", "--assembly", "R64-1-1", "--out", &yeast]),
        clean(0, &one)
    );
    assert_eq!(
        at(&manifest(&yeast), &["species"]).as_str(),
        Some("saccharomyces_cerevisiae")
    );
}

#[test]
fn select_focus_terms_take_more_exemplars() {
    let scratch = scratch("select-focus");
    let (mut rows, mut lines) = (String::new(), String::new());
    for i in 1..6 {
        let (name, location) = (format!("rs{i}"), format!("21:{}", i * 100));
        rows += &row(&name, &location, "G", "ENST1", "inframe_insertion");
        rows += &row(&name, &location, "G", "ENST2", "intron_variant");
        lines += &format!("21\t{}\trs{i}\tA\tG\t.\tPASS\t.\n", i * 100);
    }
    let reference = write(&scratch.path, "ref.txt", &rows);
    let input = write(&scratch.path, "in.vcf", &format!("{VCF_HEADER}{lines}"));
    let out = under(&scratch.path, "corpus");
    let dataset = format!("s01={reference}={input}");
    let done = run(&[
        "select",
        "--release",
        "115",
        "--assembly",
        "GRCh37",
        "--dataset",
        &dataset,
        "--out",
        &out,
        "--k",
        "1",
        "--focus-term",
        "inframe_insertion",
        "--focus-k",
        "4",
        "--seed",
        "1",
        "--multi-allelic-per-dataset",
        "0",
        "--covers",
        "symbol",
    ]);
    assert_eq!(done, clean(0, &summary(2, 4, "snv 5")));
    let m = manifest(&out);
    assert_eq!(
        text_of(at(&m, &["focus_terms"])),
        r#"["inframe_insertion"]"#
    );
    assert_eq!(
        at(&m, &["focus_exemplars_per_combination"]).as_i64(),
        Some(4)
    );
    let exemplars =
        |combination: &str| at(&m, &["combinations", combination, "exemplar_records"]).as_i64();
    assert_eq!(exemplars("inframe_insertion"), Some(4));
    // A record's combinations are the ones it was chosen for, so intron_variant keeps --k.
    assert_eq!(exemplars("intron_variant"), Some(1));
    assert!((4..=5).contains(&records(&m).len()));
    // Without --flag or --fasta-name the manifest carries neither key.
    assert!(!m.contains_key("flags") && !m.contains_key("fasta"));
}

/// A whole cache of one shard on chromosome 21 that holds `records`, a JSON list.
fn cache_of(directory: &Path, records: &str) -> String {
    let shards = directory.join("cache/transcripts/21");
    fs::create_dir_all(&shards).unwrap();
    fs::write(shards.join("1-1000000.json"), records).unwrap();
    under(directory, "cache")
}

fn sampled_row(uv: &str, location: &str, allele: &str, feature: &str, consequence: &str) -> Row {
    Row {
        dataset: "sel".to_string(),
        uploaded_variation: uv.to_string(),
        location: location.to_string(),
        allele: allele.to_string(),
        feature: feature.to_string(),
        feature_type: "Transcript".to_string(),
        consequence_set: consequence.to_string(),
        strand: "1".to_string(),
        flags: String::new(),
        positions: ["-", "-", "-"].map(str::to_string),
        family: Vec::new(),
        substitution: false,
    }
}

#[test]
fn select_with_a_cache_takes_the_exemplar_that_costs_the_fewest_bytes() {
    let scratch = scratch("select-cost");
    let peptide = "MKV".repeat(4000);
    let heavy = format!(
        r#"{{"stable_id": "ENST_HEAVY", "start": "100", "end": "900", "variation_effect_feature_cache": {{"peptide": "{peptide}"}}}}"#
    );
    let light = r#"{"stable_id": "ENST_LIGHT", "start": "500100", "end": "500900"}"#;
    let cache = cache_of(&scratch.path, &format!("[{heavy}, {light}]"));
    let rows = row("rs_heavy", "21:500", "G", "ENST_HEAVY", "missense_variant")
        + &row(
            "rs_light",
            "21:500500",
            "G",
            "ENST_LIGHT",
            "missense_variant",
        );
    let reference = write(&scratch.path, "ref.txt", &rows);
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!(
            "{VCF_HEADER}21\t500\trs_heavy\tA\tG\t.\tPASS\t.\n21\t500500\trs_light\tA\tG\t.\tPASS\t.\n"
        ),
    );
    let dataset = format!("sel={reference}={input}");
    let common = [
        "select",
        "--release",
        "116",
        "--assembly",
        "GRCh38",
        "--dataset",
        &dataset,
        "--k",
        "1",
        "--multi-allelic-per-dataset",
        "0",
        "--covers",
        "symbol",
    ];
    let one = summary(1, 1, "snv 1");
    let mut picked = BTreeSet::new();
    for seed in 1..7 {
        let seed = seed.to_string();
        let out = under(&scratch.path, &format!("plain{seed}"));
        let mut args = common.to_vec();
        args.extend(["--out", &out, "--seed", &seed]);
        assert_eq!(run(&args), clean(0, &one));
        picked.insert(record_ids(&out)[0].clone());
        let out = under(&scratch.path, &format!("costed{seed}"));
        let mut args = common.to_vec();
        args.extend(["--out", &out, "--seed", &seed, "--cache", &cache]);
        assert_eq!(run(&args), clean(0, &one));
        assert_eq!(record_ids(&out), ["rs_light"]);
        assert_eq!(
            at(&manifest(&out), &["cache_cost_flank"]).as_i64(),
            Some(5000)
        );
    }
    // Without a cache the tie is broken by the seed.
    assert_eq!(
        picked,
        BTreeSet::from(["rs_heavy".to_string(), "rs_light".to_string()])
    );
    let cost = CacheCost::new(&cache, 5000).unwrap();
    let far = sampled_row("x", "21:600000", "G", "ENST_LIGHT", "intron_variant");
    assert_eq!(cost.row_cost(&far), Ok(0), "nothing within the flank");
    let mate = sampled_row(
        "x",
        "21:600000",
        "N[21:400[",
        "ENST_HEAVY",
        "feature_truncation",
    );
    let whole = cost.span_cost("21", 400, 600000).unwrap();
    let light_alone = cost.span_cost("21", 500500, 500500).unwrap();
    assert_eq!(cost.row_cost(&mate), Ok(whole));
    assert!(whole > light_alone && light_alone > 0);
}

#[test]
fn cache_cost_skips_a_plain_file_under_transcripts() {
    // A plain file beside the sequence directories of `transcripts/` names no sequence, so
    // a cache that holds one costs and selects as the same cache without the file.
    let scratch = scratch("select-plain-file");
    let rows = row("rs_heavy", "21:500", "G", "ENST_HEAVY", "missense_variant")
        + &row(
            "rs_light",
            "21:500500",
            "G",
            "ENST_LIGHT",
            "missense_variant",
        );
    let reference = write(&scratch.path, "ref.txt", &rows);
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!(
            "{VCF_HEADER}21\t500\trs_heavy\tA\tG\t.\tPASS\t.\n21\t500500\trs_light\tA\tG\t.\tPASS\t.\n"
        ),
    );
    let dataset = format!("sel={reference}={input}");
    let note: Vec<String> = (0..4000).map(|number| number.to_string()).collect();
    let heavy = format!(
        r#"{{"stable_id": "ENST_HEAVY", "start": "100", "end": "900", "note": "{}"}}"#,
        note.join(" ")
    );
    let light = r#"{"stable_id": "ENST_LIGHT", "start": "500100", "end": "500900"}"#;
    let mut written = Vec::new();
    for name in ["directories_alone", "with_a_plain_file"] {
        let cache = scratch.path.join(name).join("cache");
        let sequence = cache.join("transcripts/21");
        fs::create_dir_all(&sequence).unwrap();
        fs::write(
            sequence.join("1-1000000.json"),
            format!("[{heavy}, {light}]"),
        )
        .unwrap();
        if name == "with_a_plain_file" {
            fs::write(cache.join("transcripts/notes.txt"), "no sequence\n").unwrap();
        }
        let cache = cache.to_str().unwrap();
        assert_eq!(CacheCost::new(cache, 5000).unwrap().sequences(), ["21"]);
        let out = under(&scratch.path.join(name), "corpus");
        let ran = run(&[
            "select",
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--dataset",
            &dataset,
            "--out",
            &out,
            "--k",
            "1",
            "--seed",
            "1",
            "--multi-allelic-per-dataset",
            "0",
            "--covers",
            "symbol",
            "--cache",
            cache,
        ]);
        assert_eq!(ran, clean(0, &summary(1, 1, "snv 1")), "{name}");
        assert_eq!(record_ids(&out), ["rs_light"], "{name}");
        written.push((
            fs::read(Path::new(&out).join("variants.vcf")).unwrap(),
            fs::read(manifest_file(&out)).unwrap(),
        ));
    }
    assert_eq!(written[0], written[1]);
}

#[test]
fn select_keeps_an_exemplar_of_every_variant_class() {
    let scratch = scratch("select-classes");
    let rows = row("rs1", "21:100", "G", "ENST1", "feature_truncation")
        + &row(
            "bnd1",
            "21:300000",
            "[21:300900[TTGC",
            "ENST2",
            "feature_truncation",
        );
    let reference = write(&scratch.path, "ref.txt", &rows);
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!(
            "{VCF_HEADER}21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t300000\tbnd1\tC\t[21:300900[TTGC\t.\tPASS\tSVTYPE=BND\n"
        ),
    );
    let peptide = "MKV".repeat(400);
    let cache = cache_of(
        &scratch.path,
        &format!(
            r#"[{{"stable_id": "ENST1", "start": "50", "end": "150"}}, {{"stable_id": "ENST2", "start": "299000", "end": "301000", "variation_effect_feature_cache": {{"peptide": "{peptide}"}}}}]"#
        ),
    );
    let out = under(&scratch.path, "corpus");
    let dataset = format!("sel={reference}={input}");
    let done = run(&[
        "select",
        "--release",
        "116",
        "--assembly",
        "GRCh38",
        "--dataset",
        &dataset,
        "--out",
        &out,
        "--k",
        "1",
        "--seed",
        "1",
        "--multi-allelic-per-dataset",
        "0",
        "--cache",
        &cache,
        "--covers",
        "symbol",
    ]);
    assert_eq!(
        done,
        clean(0, &summary(1, 2, "breakend_insertion 1, snv 1"))
    );
    let m = manifest(&out);
    for id in ["rs1", "bnd1"] {
        let combinations = record_by_id(&m, id).get("combinations").unwrap();
        assert_eq!(text_of(combinations), r#"["feature_truncation"]"#, "{id}");
    }
    assert_eq!(
        text_of(at(&m, &["variant_classes"])),
        r#"{"breakend_insertion": 1, "snv": 1}"#
    );
    let exemplars = at(
        &m,
        &["combinations", "feature_truncation", "exemplar_records"],
    );
    assert_eq!(exemplars.as_i64(), Some(2));
    // Without the class guarantee the cheaper SNV alone would carry the combination.
    let rows: Vec<Row> = iter_rows::<&str>(&reference, "sel", &[], None)
        .collect::<Result<_, _>>()
        .unwrap();
    let mut reservoir: OrderedMap<String, Sample> = OrderedMap::new();
    let seen = rows.len();
    reservoir.insert("feature_truncation".to_string(), Sample { seen, rows });
    let cost = CacheCost::new(&cache, 5000).unwrap();
    let chosen = choose_exemplars(
        &reservoir,
        1,
        &mut Random::new(1),
        0,
        &[],
        None,
        Some(&cost),
    );
    let chosen = chosen.unwrap();
    let names: Vec<&str> = chosen
        .get("feature_truncation")
        .unwrap()
        .iter()
        .map(|row| row.uploaded_variation.as_str())
        .collect();
    assert_eq!(names, ["rs1"]);
}

#[test]
fn select_writes_covers() {
    let scratch = scratch("select-covers");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    let extra = [
        "--covers",
        "symbol",
        "--covers",
        "max_af",
        "--covers",
        "symbol",
        "--flag=--symbol",
    ];
    assert_eq!(
        select(&reference, &input, &out, &extra),
        clean(0, &summary(2, 2, "snv 2"))
    );
    let m = manifest(&out);
    assert_eq!(text_of(at(&m, &["covers"])), r#"["symbol", "max_af"]"#);
    assert!(key_place(&m, "flags") < key_place(&m, "covers"));
    assert!(key_place(&m, "covers") < key_place(&m, "datasets"));
    assert!(!m.contains_key("flag_sets"));
}

#[test]
fn select_refuses_a_corpus_without_covers() {
    let scratch = scratch("select-no-covers");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    assert_eq!(
        select(&reference, &input, &out, &[]),
        failed(
            "select writes no manifest without --covers: name the parity-matrix rows the corpus exercises"
        )
    );
    assert!(!exists(&out));
}

#[test]
fn select_refuses_a_cache_without_transcripts() {
    let scratch = scratch("select-no-transcripts");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    let cache = under(&scratch.path, "cache");
    fs::create_dir(&cache).unwrap();
    for cache in [cache, under(&scratch.path, "absent")] {
        assert_eq!(
            select(
                &reference,
                &input,
                &out,
                &["--covers", "symbol", "--cache", &cache]
            ),
            failed(&format!("{cache} has no transcripts/ directory"))
        );
        assert!(!exists(&out));
    }
    // A reference that cannot be read is refused the same way.
    let absent = under(&scratch.path, "absent.txt");
    assert_eq!(
        select(&absent, &input, &out, &["--covers", "symbol"]),
        failed(&format!(
            "{absent} is not a readable text file ([Errno 2] No such file or directory: '{absent}')"
        ))
    );
}

#[test]
fn select_preserves_covers_on_regeneration() {
    let scratch = scratch("select-regeneration");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    let matrix = matrix(&scratch.path, &["symbol", "max_af", "pick"]);
    let two = summary(2, 2, "snv 2");
    assert_eq!(
        select(&reference, &input, &out, &["--covers", "symbol"]),
        clean(0, &two)
    );
    let set = run(&[
        "flag-set",
        &out,
        "--name",
        "pick",
        "--flag=--pick",
        "--covers",
        "pick",
        "--matrix",
        &matrix,
    ]);
    assert_eq!(set.code, 0);
    let mut m = manifest(&out);
    let entry = text_of(&at(&m, &["flag_sets"]).as_list().unwrap()[0]);
    let notice = "could not fetch sequence for tandem repeats; consequence calculation for tandem repeats is less precise when using --offline without --fasta";
    let plugin_files = r#"{"LoFTEE": ["plugin_data/human_ancestor.fa.gz", "plugin_data/gerp.bw"]}"#;
    let divergences = r#"[{"location": "21:100", "record_indices": [0]}]"#;
    let value = |text: &str| vep_golden_corpus::json::loads(text).unwrap();
    m.insert("masks", value(r#"["regulatory_feature_rows"]"#));
    m.insert("reference_warnings", Value::strings(&[notice]));
    m.insert("plugin_files", value(plugin_files));
    m.insert("divergences", value(divergences));
    let sets = m.get_mut("flag_sets").and_then(Value::as_list_mut).unwrap();
    let set = sets[0].as_object_mut().unwrap();
    set.insert("divergences", value(divergences));
    set.insert("divergence_summary", value("{}"));
    write_manifest_line(&out, &m);
    // The list the manifest holds satisfies the covers requirement.
    assert_eq!(select(&reference, &input, &out, &[]), clean(0, &two));
    let m = manifest(&out);
    let carried = |m: &vep_golden_corpus::json::Object| {
        assert_eq!(text_of(at(m, &["masks"])), r#"["regulatory_feature_rows"]"#);
        assert_eq!(at(m, &["reference_warnings"]), &Value::strings(&[notice]));
        assert_eq!(text_of(at(m, &["plugin_files"])), plugin_files);
        assert_eq!(text_of(at(m, &["flag_sets"])), format!("[{entry}]"));
    };
    assert_eq!(text_of(at(&m, &["covers"])), r#"["symbol"]"#);
    carried(&m);
    assert!(!m.contains_key("divergences"));
    assert_eq!(
        select(
            &reference,
            &input,
            &out,
            &["--covers", "max_af", "--covers", "symbol"]
        ),
        clean(0, &two)
    );
    let mut m = manifest(&out);
    assert_eq!(text_of(at(&m, &["covers"])), r#"["symbol", "max_af"]"#);
    carried(&m);
    let order = [
        "covers",
        "masks",
        "reference_warnings",
        "plugin_files",
        "flag_sets",
        "datasets",
    ];
    let places: Vec<usize> = order.iter().map(|key| key_place(&m, key)).collect();
    assert!(places.windows(2).all(|pair| pair[0] < pair[1]));
    m.insert("masks", value(r#"["splice_rows"]"#));
    write_manifest_line(&out, &m);
    let variants = under(Path::new(&out), "variants.vcf");
    fs::remove_file(&variants).unwrap();
    let written = manifest_text(&out);
    let named = manifest_file(&out);
    let named = named.to_str().unwrap();
    assert_eq!(
        select(&reference, &input, &out, &[]),
        failed(&format!(
            "{named} masks names unknown classes ['splice_rows'] (the classes are {})",
            mask_classes().join(", ")
        ))
    );
    assert!(
        manifest_text(&out) == written && !exists(&variants),
        "nothing is written"
    );
    m.insert("masks", value(r#"["regulatory_feature_rows"]"#));
    m.insert("reference_warnings", notice);
    write_manifest_line(&out, &m);
    let written = manifest_text(&out);
    assert_eq!(
        select(&reference, &input, &out, &[]),
        failed(&format!(
            "{named} reference_warnings is not a list of texts"
        ))
    );
    assert!(
        manifest_text(&out) == written && !exists(&variants),
        "nothing is written"
    );
}

#[test]
fn select_matrix_and_empty_names() {
    let scratch = scratch("select-matrix");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    let matrix = matrix(&scratch.path, &["symbol"]);
    assert_eq!(
        select(
            &reference,
            &input,
            &out,
            &["--covers", "symbol", "--covers", "ghost", "--matrix", &matrix]
        ),
        failed(&format!(
            "covers names ghost, which is not a row of {matrix}"
        ))
    );
    assert!(!exists(&out));
    assert_eq!(
        select(&reference, &input, &out, &["--covers", ""]),
        failed("--covers wants a non-empty row name")
    );
    assert!(!exists(&out));
    assert_eq!(
        select(
            &reference,
            &input,
            &out,
            &["--covers", "symbol", "--matrix", &matrix]
        ),
        clean(0, &summary(2, 2, "snv 2"))
    );
}

#[test]
fn select_include_carries_the_named_record_outside_the_sampling() {
    let scratch = scratch("select-include");
    let (reference, input) = one_stratum_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    assert_eq!(
        select(&reference, &input, &out, &["--covers", "symbol"]),
        clean(0, &summary(1, 1, "snv 1"))
    );
    // Seed 1 samples rs2, so the other record of the stratum is rs1.
    assert_eq!(record_ids(&out), ["rs2"]);
    let other = "rs1";
    let done = select(
        &reference,
        &input,
        &out,
        &["--covers", "symbol", "--include", "s01=rs1"],
    );
    let named = summary(1, 2, "snv 2").replace('\n', ", named records 1\n");
    assert_eq!(done, clean(0, &named));
    let m = manifest(&out);
    assert_eq!(record_ids(&out), ["rs1", "rs2"]);
    assert_eq!(
        text_of(at(&m, &["included_records"])),
        r#"{"s01": ["rs1"]}"#
    );
    let combinations = record_by_id(&m, other).get("combinations").unwrap();
    assert_eq!(text_of(combinations), r#"["missense_variant"]"#);
    assert!(!m.contains_key("kept_records"));
    assert!(key_place(&m, "exemplars_per_combination") < key_place(&m, "included_records"));
    assert!(key_place(&m, "included_records") < key_place(&m, "max_span"));
    let variants = under(Path::new(&out), "variants.vcf");
    fs::remove_file(&variants).unwrap();
    let written = manifest_text(&out);
    assert_eq!(
        select(
            &reference,
            &input,
            &out,
            &["--include", "s01=ghost", "--include", "s01=rs1"]
        ),
        failed("--include s01=ghost has no row in the s01 reference")
    );
    assert!(
        !exists(&variants) && manifest_text(&out) == written,
        "nothing is written"
    );
    assert_eq!(
        select(&reference, &input, &out, &["--include", "s02=rs1"]),
        failed("--include s02=rs1 names dataset 's02', which is not among --dataset")
    );
    assert_eq!(
        select(&reference, &input, &out, &["--include", "rs1"]),
        failed("--include wants sid=name, got 'rs1'")
    );
}

#[test]
fn select_keep_records_keeps_every_record_of_the_manifest() {
    let scratch = scratch("select-keep");
    let (reference, input) = one_stratum_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    let select = |extra: &[&str]| select(&reference, &input, &out, extra);
    let (one, two) = (summary(1, 1, "snv 1"), summary(1, 2, "snv 2"));
    let two_named = two.replace('\n', ", named records 2\n");
    let both = ["rs1", "rs2"];
    let included = |out: &str| manifest(out).get("included_records").map(text_of);
    assert_eq!(select(&["--covers", "symbol", "--k", "2"]), clean(0, &two));
    assert_eq!(record_ids(&out), both);
    assert_eq!(select(&[]), clean(0, &one));
    assert_eq!(
        record_ids(&out).len(),
        1,
        "the sampling alone drops one of the two"
    );
    assert_eq!(select(&["--k", "2"]), clean(0, &two));
    assert_eq!(select(&["--keep-records"]), clean(0, &two_named));
    let mut m = manifest(&out);
    assert_eq!(record_ids(&out), both);
    assert_eq!(at(&m, &["kept_records"]).as_i64(), Some(2));
    assert!(!m.contains_key("included_records"));
    assert_eq!(
        select(&["--keep-records", "--max-records", "1"]),
        clean(0, &two_named)
    );
    assert_eq!(
        record_ids(&out),
        both,
        "a kept record survives the record cap"
    );
    // A record both kept and included counts once.
    assert_eq!(
        select(&["--keep-records", "--include", "s01=rs1"]),
        clean(0, &two_named)
    );
    assert_eq!(included(&out).as_deref(), Some(r#"{"s01": ["rs1"]}"#));
    assert_eq!(select(&["--keep-records"]), clean(0, &two_named));
    // A regeneration that keeps the records keeps their included names.
    assert_eq!(included(&out).as_deref(), Some(r#"{"s01": ["rs1"]}"#));
    assert_eq!(
        select(&["--keep-records", "--include", "s01=rs2"]),
        clean(0, &two_named)
    );
    // Carried and new names merge.
    assert_eq!(
        included(&out).as_deref(),
        Some(r#"{"s01": ["rs1", "rs2"]}"#)
    );
    assert_eq!(select(&[]), clean(0, &one));
    // A sampling run carries no included names.
    assert_eq!(included(&out), None);
    let value = |text: &str| vep_golden_corpus::json::loads(text).unwrap();
    let variants = under(Path::new(&out), "variants.vcf");
    m.insert(
        "included_records",
        value(r#"{"s01": ["rs1"], "s02": ["rs9"]}"#),
    );
    write_manifest_line(&out, &m);
    assert_eq!(
        select(&["--keep-records"]),
        failed("carried include s02=rs9 names dataset 's02', which is not among --dataset")
    );
    m.insert("included_records", value(r#"{"s01": ["rs1", "rs9"]}"#));
    write_manifest_line(&out, &m);
    fs::remove_file(&variants).unwrap();
    let written = manifest_text(&out);
    assert_eq!(
        select(&["--keep-records"]),
        failed("carried include s01=rs9 has no row in the s01 reference")
    );
    assert!(
        !exists(&variants) && manifest_text(&out) == written,
        "nothing is written"
    );
    let named = manifest_file(&out);
    let named = named.to_str().unwrap();
    for malformed in [
        r#"["rs1"]"#,
        r#"{"s01": []}"#,
        r#"{"": ["rs1"]}"#,
        r#"{"s01": [""]}"#,
    ] {
        m.insert("included_records", value(malformed));
        write_manifest_line(&out, &m);
        assert_eq!(
            select(&["--keep-records"]),
            failed(&format!(
                "{named} included_records is not a map of dataset to record names"
            )),
            "{malformed}"
        );
    }
    m.insert("included_records", value(r#"{"s01": ["rs1"]}"#));
    let extra = |dataset: &str| {
        value(&format!(
            r#"{{"dataset": "{dataset}", "vcf": "21\t300\trs3\tA\tG\t.\tPASS\t.", "combinations": []}}"#
        ))
    };
    let kept = m.get_mut("records").and_then(Value::as_list_mut).unwrap();
    kept.push(extra("s01"));
    write_manifest_line(&out, &m);
    let written = manifest_text(&out);
    assert_eq!(
        select(&["--keep-records"]),
        failed("kept record 21 300 rs3 A G has no row in the s01 reference")
    );
    assert!(
        !exists(&variants) && manifest_text(&out) == written,
        "nothing is written"
    );
    let kept = m.get_mut("records").and_then(Value::as_list_mut).unwrap();
    kept.pop();
    kept.push(extra("s02"));
    write_manifest_line(&out, &m);
    assert_eq!(
        select(&["--keep-records"]),
        failed("kept record 21 300 rs3 A G names dataset 's02', which is not among --dataset")
    );
    let fresh = under(&scratch.path, "fresh");
    assert_eq!(
        common::select(
            &reference,
            &input,
            &fresh,
            &["--covers", "symbol", "--keep-records"]
        ),
        failed(&format!(
            "--keep-records wants a manifest to regenerate in {fresh}"
        ))
    );
    assert!(!exists(&fresh));
}

/// The first flag set of a manifest, open to change.
fn first_set(manifest: &mut Object) -> &mut Object {
    let sets = manifest.get_mut("flag_sets").and_then(Value::as_list_mut);
    sets.and_then(|sets| sets[0].as_object_mut())
        .expect("a flag set")
}

#[test]
fn select_refuses_to_drop_a_with_reference_key() {
    // `select` writes a manifest without the divergence keys of the one it regenerates, and
    // no `classify` writes the with-reference sets again: an entry that carries the key, on
    // the manifest's own list or a flag set's, stops the run before a reference is read,
    // each named with its owner, and nothing is written. `--drop-reference-statements`
    // names each on standard error and regenerates without them; a list whose entries carry
    // no such key goes with no option and no line.
    let scratch = scratch("select-statements");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    let matrix = matrix(&scratch.path, &["symbol", "pick"]);
    let two = summary(2, 2, "snv 2");
    assert_eq!(
        select(&reference, &input, &out, &["--covers", "symbol"]),
        clean(0, &two)
    );
    let set = [
        "flag-set",
        &out,
        "--name",
        "pick",
        "--flag=--pick",
        "--covers",
        "pick",
        "--matrix",
        &matrix,
    ];
    assert_eq!(run(&set).code, 0);
    let regenerated = manifest(&out);
    let plain = loads(
        r#"{"location": "21:200", "allele": "T", "feature": "ENST2", "record_indices": [1]}"#,
    )
    .unwrap();
    let mut stated = manifest(&out);
    stated.insert(
        "divergences",
        vec![
            reference_statement("21:100", "G", "ENST1", 0),
            plain.clone(),
        ],
    );
    first_set(&mut stated).insert(
        "divergences",
        vec![
            plain.clone(),
            reference_statement("21:200", "T", "ENST2", 1),
        ],
    );
    write_manifest_line(&out, &stated);
    let held = manifest_text(&out);
    let vcf = Path::new(&out).join("variants.vcf");
    let records = fs::read(&vcf).unwrap();
    let refusal = failed(&format!(
        "select would discard {WITH_REFERENCE_SETS}, a key written by hand, on the manifest's divergences entry on 21:100 G ENST1; flag set pick divergences entry on 21:200 T ENST2: pass --drop-reference-statements to discard it"
    ));
    // The dataset names a reference that is no file: the refusal comes before a reference
    // is read.
    let absent = under(&scratch.path, "absent.txt");
    assert_eq!(select(&absent, &input, &out, &[]), refusal);
    assert_eq!(select(&reference, &input, &out, &[]), refusal);
    assert_eq!(manifest_text(&out), held);
    assert_eq!(fs::read(&vcf).unwrap(), records);
    // A run the option lets past and another check refuses discards nothing and names no
    // entry.
    assert_eq!(
        select(&absent, &input, &out, &["--drop-reference-statements"]),
        failed(&format!(
            "{absent} is not a readable text file ([Errno 2] No such file or directory: '{absent}')"
        ))
    );
    assert_eq!(manifest_text(&out), held);
    // A manifest that cannot be written is found before an entry is named: the run ends on
    // the one line of the write. A user who may write any file skips the check.
    let file = manifest_file(&out);
    let mode = |bits: u32| {
        fs::set_permissions(&file, fs::Permissions::from_mode(bits)).unwrap();
    };
    mode(0o444);
    if fs::OpenOptions::new().write(true).open(&file).is_err() {
        let unwritten = select(&reference, &input, &out, &["--drop-reference-statements"]);
        let line = format!(
            "ERROR: [vep-golden-corpus] {} is not writable (",
            file.display()
        );
        assert_eq!((unwritten.code, unwritten.stdout.as_str()), (1, ""));
        assert!(unwritten.stderr.starts_with(&line), "{}", unwritten.stderr);
        assert_eq!(unwritten.stderr.lines().count(), 1, "{}", unwritten.stderr);
        assert_eq!(manifest_text(&out), held);
    }
    mode(0o644);
    let discards = |entry: &str| {
        format!("WARNING: [build_golden_corpus] select discards {WITH_REFERENCE_SETS} on {entry}\n")
    };
    assert_eq!(
        select(&reference, &input, &out, &["--drop-reference-statements"]),
        Run {
            code: 0,
            stdout: two.clone(),
            stderr: discards("the manifest's divergences entry on 21:100 G ENST1")
                + &discards("flag set pick divergences entry on 21:200 T ENST2"),
        }
    );
    assert_eq!(manifest(&out), regenerated);
    let mut unstated = manifest(&out);
    unstated.insert("divergences", vec![plain.clone()]);
    first_set(&mut unstated).insert("divergences", vec![plain]);
    write_manifest_line(&out, &unstated);
    assert_eq!(select(&reference, &input, &out, &[]), clean(0, &two));
    assert_eq!(manifest(&out), regenerated);
}

#[test]
fn select_keeps_the_generated_block() {
    // A regeneration carries the manifest's `generated` block (the command that writes the
    // files the corpus does not commit, and their digests) as it stands, after
    // `plugin_files` and before `datasets`; a manifest without a block gains none.
    let scratch = scratch("select-generated");
    let (reference, input) = two_record_dataset(&scratch.path);
    let out = under(&scratch.path, "corpus");
    assert_eq!(
        select(&reference, &input, &out, &["--covers", "symbol"]).code,
        0
    );
    let mut carried = manifest(&out);
    assert!(!carried.contains_key("generated"));
    assert_eq!(select(&reference, &input, &out, &[]).code, 0);
    assert!(!manifest(&out).contains_key("generated"));
    let block = loads(&format!(
        r#"{{"command": ["python3", "scripts/fixtures/write_set.py", "--out", "{{generated}}/set"], "sha256": {{"set/b": "{}", "set/a": "{}"}}}}"#,
        "1".repeat(64),
        "0".repeat(64)
    ))
    .unwrap();
    assert_eq!(generated_shape_error(&block), None);
    carried.insert("generated", block.clone());
    carried.insert(
        "plugin_files",
        loads(r#"{"MaxEntScan": ["plugin_data/notes.txt"]}"#).unwrap(),
    );
    write_manifest_line(&out, &carried);
    assert_eq!(select(&reference, &input, &out, &[]).code, 0);
    let written = manifest(&out);
    assert_eq!(written.get("generated"), Some(&block));
    let digests = at(&written, &["generated", "sha256"]).as_object().unwrap();
    assert_eq!(digests.keys().collect::<Vec<_>>(), ["set/b", "set/a"]);
    assert!(key_place(&written, "plugin_files") < key_place(&written, "generated"));
    assert!(key_place(&written, "generated") < key_place(&written, "datasets"));
}
