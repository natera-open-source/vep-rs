// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `select` and `check` with `--family-columns` and `--record-classes`: every family
//! stratum and every record class the references show keeps its exemplars, the manifest
//! records them, and `check` names the ones a corpus lacks.

mod common;

use common::{
    at, clean, failed, feature_row, manifest, matrix, record_by_id, record_ids, records, row,
    row_with, run, scratch, select, text_of, under, vcf_of, write, write_manifest_line, Run,
    VCF_HEADER,
};
use std::fs;
use vep_golden_corpus::json::Value;

const OBSERVED_ONE: &str = "observed 1, covered by corpus 1, missing 0\n";
const OBSERVED_TWO: &str = "observed 2, covered by corpus 2, missing 0\n";

/// The summary line of a run that covered everything it observed.
fn summary(combinations: usize, records: usize, classes: &str) -> String {
    let count = classes.split(", ").count();
    format!(
        "combinations observed {combinations}, covered {combinations}, records {records}, unresolved exemplars 0, variant classes {count} ({classes})\n"
    )
}

#[test]
fn a_two_base_substitution_is_classed_from_its_input_record() {
    // A two-base substitution and a two-base insertion print the same Location and Allele,
    // so such a row is classed with its input record.
    let scratch = scratch("strata-substitution");
    let coding =
        |name: &str, location: &str| row(name, location, "TG", "ENST1", "coding_sequence_variant");
    let reference = write(
        &scratch.path,
        "ref.txt",
        &(coding("mnv1", "21:100-101")
            + &coding("ins1", "21:200-201")
            + &coding("21_300_AC/TG", "21:300-301")),
    );
    let insertion = "21\t200\tins1\tA\tATG\t.\tPASS\t.\n";
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!(
            "{VCF_HEADER}21\t100\tmnv1\tAC\tTG\t.\tPASS\t.\n{insertion}21\t300\t.\tAC\tTG\t.\tPASS\t.\n"
        ),
    );
    let family = ["--covers", "symbol", "--family-columns", "IMPACT="];

    // One exemplar a combination: the second class enters as a stratum of its own.
    let out = under(&scratch.path, "corpus");
    assert_eq!(
        select(&reference, &input, &out, &family),
        clean(0, &summary(1, 2, "insertion 1, mnv_or_complex 1"))
    );
    let m = manifest(&out);
    assert_eq!(
        text_of(at(&m, &["variant_classes"])),
        r#"{"insertion": 1, "mnv_or_complex": 1}"#
    );
    assert_eq!(
        text_of(at(&m, &["family_strata"])),
        r#"{"IMPACT=MODIFIER;insertion": 1, "IMPACT=MODIFIER;mnv_or_complex": 1}"#
    );
    let refs: Vec<&str> = records(&m)
        .iter()
        .map(|record| vcf_of(record).split('\t').nth(3).unwrap())
        .collect();
    assert_eq!(refs, ["AC", "A"]);
    let dataset = format!("s01={reference}={input}");
    let check = |corpus: &str, dataset: &str| {
        run(&[
            "check",
            "--family-columns",
            "IMPACT=",
            "--corpus",
            corpus,
            "--dataset",
            dataset,
        ])
    };
    assert_eq!(
        check(&out, &dataset),
        clean(
            0,
            &format!("{OBSERVED_ONE}family strata observed 2, covered by corpus 2, missing 0\n")
        )
    );

    // A corpus of the insertion alone lacks the substitution's stratum.
    let ins_reference = write(&scratch.path, "ins.txt", &coding("ins1", "21:200-201"));
    let ins_input = write(
        &scratch.path,
        "ins.vcf",
        &format!("{VCF_HEADER}{insertion}"),
    );
    let lacking = under(&scratch.path, "lacking");
    assert_eq!(
        select(&ins_reference, &ins_input, &lacking, &family),
        clean(0, &summary(1, 1, "insertion 1"))
    );
    assert_eq!(
        text_of(at(&manifest(&lacking), &["variant_classes"])),
        r#"{"insertion": 1}"#
    );
    assert_eq!(
        check(&lacking, &dataset),
        clean(
            1,
            &format!(
                "{OBSERVED_ONE}family strata observed 2, covered by corpus 1, missing 1\n  MISSING family stratum IMPACT=MODIFIER;mnv_or_complex\n"
            )
        )
    );

    // With family columns `check` reads the input, and refuses an unreadable one after the
    // errors found before it.
    let absent = under(&scratch.path, "absent.vcf");
    let matrix = matrix(&scratch.path, &["other"]);
    let absent_dataset = format!("s01={reference}={absent}");
    assert_eq!(
        run(&[
            "check",
            "--family-columns",
            "IMPACT=",
            "--corpus",
            &out,
            "--dataset",
            &absent_dataset,
            "--matrix",
            &matrix
        ]),
        Run {
            code: 1,
            stdout: String::new(),
            stderr: format!(
                "ERROR: [build_golden_corpus] covers names symbol, which is not a row of {matrix}\nERROR: [build_golden_corpus] {absent} is not a readable text file ([Errno 2] No such file or directory: '{absent}')\n"
            ),
        }
    );

    // The record-class fill adds every row of the records it takes, read with the input too.
    let motif = |name: &str, location: &str, allele: &str, feature: &str| {
        feature_row(
            name,
            location,
            allele,
            feature,
            "MotifFeature",
            "TF_binding_site_variant",
            "IMPACT=MODIFIER;STRAND=1",
        )
    };
    let mut rows = motif("snv1", "21:300", "G", "ENSM1");
    for (name, location) in [("subA", "21:100-101"), ("subB", "21:500-501")] {
        for feature in ["ENSM1", "ENSM2"] {
            rows += &motif(name, location, "TG", feature);
        }
    }
    let fill_reference = write(&scratch.path, "fill.txt", &rows);
    let fill_input = write(
        &scratch.path,
        "fill.vcf",
        &format!(
            "{VCF_HEADER}21\t100\tsubA\tAC\tTG\t.\tPASS\t.\n21\t300\tsnv1\tA\tG\t.\tPASS\t.\n21\t500\tsubB\tAC\tTG\t.\tPASS\t.\n"
        ),
    );
    let fill = under(&scratch.path, "fill");
    let classes = [
        "--covers",
        "symbol",
        "--record-classes",
        "regulatory_class",
        "--record-class-exemplars",
        "4",
    ];
    assert_eq!(
        select(&fill_reference, &fill_input, &fill, &classes),
        clean(
            0,
            &summary(1, 3, "mnv_or_complex 3, snv 1")
                .replace('\n', ", record classes observed 4, covered 4\n")
        )
    );
    let filled = manifest(&fill);
    assert_eq!(
        text_of(at(&filled, &["variant_classes"])),
        r#"{"mnv_or_complex": 3, "snv": 1}"#
    );
    assert_eq!(
        at(&filled, &["record_classes", "motif_mnv"]).as_i64(),
        Some(2)
    );
}

#[test]
fn regulatory_class_key_covers_every_class() {
    let scratch = scratch("strata-regulatory");
    let rf = "regulatory_region_variant";
    let tf = "TF_binding_site_variant";
    let intron = |name: &str, location: &str, allele: &str| {
        row(name, location, allele, "ENST1", "intron_variant")
    };
    let regulatory = |name: &str, location: &str, allele: &str, feature: &str| {
        feature_row(
            name,
            location,
            allele,
            feature,
            "RegulatoryFeature",
            rf,
            "IMPACT=MODIFIER",
        )
    };
    let motif = |name: &str, location: &str, allele: &str, feature: &str, strand: &str| {
        let number = feature.trim_start_matches("ENSM");
        let extra = format!("IMPACT=MODIFIER;STRAND={strand};MOTIF_NAME=ENSPFM{number}");
        feature_row(name, location, allele, feature, "MotifFeature", tf, &extra)
    };
    let intergenic = |name: &str, location: &str, allele: &str| {
        feature_row(
            name,
            location,
            allele,
            "-",
            "-",
            "intergenic_variant",
            "IMPACT=MODIFIER",
        )
    };
    // Two records of one shape: the combination sampling takes one, and `rf_only` reaches
    // its third exemplar through the other. The multi-allelic record at 9000 is one the
    // stratum sampling carries: its classes are read by its names.
    let rows = [
        intron("rs_rf", "21:1000", "G"),
        regulatory("rs_rf", "21:1000", "G", "ENSR1"),
        intron("rs_rf2", "21:1100", "G"),
        regulatory("rs_rf2", "21:1100", "G", "ENSR1"),
        intron("rs_mf", "21:2000", "G"),
        motif("rs_mf", "21:2000", "G", "ENSM1", "1"),
        intron("rs_both", "21:3000", "G"),
        regulatory("rs_both", "21:3000", "G", "ENSR2"),
        motif("rs_both", "21:3000", "G", "ENSM2", "-1"),
        regulatory("rs_regonly", "21:4000", "G", "ENSR3"),
        intergenic("rs_regonly", "21:4000", "G"),
        intron("21_5000_AT/GC", "21:5000-5001", "GC"),
        motif("21_5000_AT/GC", "21:5000-5001", "GC", "ENSM3", "1"),
        intron("21_6000_AT/A", "21:6001", "-"),
        regulatory("21_6000_AT/A", "21:6001", "-", "ENSR4"),
        motif("21_6000_AT/A", "21:6001", "-", "ENSM4", "1"),
        intron("rs_plain", "21:7000", "G"),
        motif("21_9000_A/G/T", "21:9000", "G", "ENSM5", "1"),
        intergenic("21_9000_A/G/T", "21:9000", "G"),
        motif("21_9000_A/G/T", "21:9000", "T", "ENSM5", "1"),
        intergenic("21_9000_A/G/T", "21:9000", "T"),
    ]
    .concat();
    let reference = write(&scratch.path, "ref.txt", &rows);
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!(
            "{VCF_HEADER}21\t1000\trs_rf\tA\tG\t.\tPASS\t.\n21\t1100\trs_rf2\tA\tG\t.\tPASS\t.\n21\t2000\trs_mf\tA\tG\t.\tPASS\t.\n21\t3000\trs_both\tA\tG\t.\tPASS\t.\n21\t4000\trs_regonly\tA\tG\t.\tPASS\t.\n21\t5000\t.\tAT\tGC\t.\tPASS\t.\n21\t6000\t.\tAT\tA\t.\tPASS\t.\n21\t7000\trs_plain\tA\tG\t.\tPASS\t.\n21\t9000\t.\tA\tG,T\t.\tPASS\t.\n"
        ),
    );
    let out = under(&scratch.path, "corpus");
    let dataset = format!("reg={reference}={input}");
    let select_into = |out: &str, extra: &[&str]| {
        let mut args = vec![
            "select",
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--dataset",
            &dataset,
            "--out",
            out,
            "--k",
            "1",
            "--seed",
            "1",
            "--multi-allelic-per-dataset",
            "1",
            "--covers",
            "regulatory",
        ];
        args.extend(extra);
        run(&args)
    };
    let classes = |exemplars: &'static str| {
        [
            "--record-classes",
            "regulatory_class",
            "--record-class-exemplars",
            exemplars,
        ]
    };
    assert_eq!(
        select_into(&out, &classes("3")),
        clean(
            0,
            &summary(4, 8, "deletion 1, mnv_or_complex 1, snv 11")
                .replace('\n', ", record classes observed 8, covered 8\n")
        )
    );
    let mut m = manifest(&out);
    assert_eq!(
        at(&m, &["record_class_key"]).as_str(),
        Some("regulatory_class")
    );
    assert_eq!(at(&m, &["record_class_exemplars"]).as_i64(), Some(3));
    // A record by its ID, or by its POS where it has none.
    let record_classes = |m: &vep_golden_corpus::json::Object, name: &str| {
        let record = records(m).iter().find(|record| {
            let columns: Vec<&str> = vcf_of(record).split('\t').collect();
            columns[2] == name || (columns[2] == "." && columns[1] == name)
        });
        record.map(|record| text_of(record.get("record_classes").unwrap()))
    };
    let expected = [
        ("rs_rf", r#"["rf_only"]"#),
        ("rs_rf2", r#"["rf_only"]"#),
        ("rs_mf", r#"["mf_only", "motif_snv"]"#),
        ("rs_both", r#"["motif_reverse", "motif_snv", "rf_and_mf"]"#),
        ("rs_regonly", r#"["regulatory_only", "rf_only"]"#),
        ("5000", r#"["mf_only", "motif_mnv"]"#),
        ("6000", r#"["motif_indel", "rf_and_mf"]"#),
        ("9000", r#"["mf_only", "motif_snv", "regulatory_only"]"#),
    ];
    for (name, classes) in expected {
        assert_eq!(record_classes(&m, name).as_deref(), Some(classes), "{name}");
    }
    assert_eq!(record_classes(&m, "rs_plain"), None);
    assert_eq!(
        text_of(at(&m, &["record_classes"])),
        r#"{"mf_only": 3, "motif_indel": 1, "motif_mnv": 1, "motif_reverse": 1, "motif_snv": 3, "regulatory_only": 2, "rf_and_mf": 2, "rf_only": 3}"#
    );
    assert_eq!(text_of(at(&m, &["uncovered_record_classes"])), "[]");
    // With one exemplar a class the multi-allelic record enters through its stratum alone,
    // with no exemplar row, and its classes are read from the reference by its names.
    let one = under(&scratch.path, "one");
    assert_eq!(
        select_into(&one, &classes("1")),
        clean(
            0,
            &summary(4, 6, "deletion 1, mnv_or_complex 1, snv 3")
                .replace('\n', ", record classes observed 8, covered 8\n")
        )
    );
    let m_one = manifest(&one);
    let sampled = records(&m_one)
        .iter()
        .find(|record| vcf_of(record).starts_with("21\t9000\t"));
    let sampled = sampled.expect("the multi-allelic record");
    assert_eq!(text_of(sampled.get("combinations").unwrap()), "[]");
    assert_eq!(
        text_of(sampled.get("strata").unwrap()),
        r#"["multi_allelic_snv"]"#
    );
    assert_eq!(
        text_of(sampled.get("record_classes").unwrap()),
        r#"["mf_only", "motif_snv", "regulatory_only"]"#
    );
    let observed = "observed 4, covered by corpus 4, missing 0\n";
    assert_eq!(
        run(&["check", "--corpus", &out, "--dataset", &dataset]),
        clean(
            0,
            &format!("{observed}record classes observed 8, covered by corpus 8, missing 0\n")
        )
    );

    // A class the reference shows and no corpus record carries fails the check.
    let extra = feature_row(
        "rs_extra",
        "21:8000",
        "G",
        "ENSM9",
        "MotifFeature",
        tf,
        "IMPACT=MODIFIER;STRAND=-1",
    ) + &intergenic("rs_extra", "21:8000", "G");
    let wider = write(&scratch.path, "wider.txt", &(rows + &extra));
    let kept = m.get_mut("records").and_then(Value::as_list_mut).unwrap();
    kept.retain(|record| !vcf_of(record).contains("rs_both"));
    write_manifest_line(&out, &m);
    let wider_dataset = format!("reg={wider}={input}");
    assert_eq!(
        run(&["check", "--corpus", &out, "--dataset", &wider_dataset]),
        clean(
            1,
            &format!(
                "{observed}record classes observed 8, covered by corpus 7, missing 1\n  MISSING record class motif_reverse\n"
            )
        )
    );
    assert_eq!(
        select_into(&out, &["--record-classes", "splice_class"]),
        failed("--record-classes wants one of regulatory_class, got 'splice_class'")
    );
    m.insert("record_class_key", "splice_class");
    assert_eq!(
        vep_golden_corpus::manifest::manifest_shape_error(&m).as_deref(),
        Some("record_class_key is not one of regulatory_class")
    );
}

#[test]
fn family_columns_stratify_exemplars() {
    // Rows of one combination and one (class, strand) stratum differ only in the family
    // key: HGVSc present or absent, and the IMPACT value.
    let scratch = scratch("strata-family");
    let rows = row_with(
        "rs1",
        "21:100",
        "G",
        "ENST1",
        "missense_variant",
        "IMPACT=MODERATE;STRAND=1;HGVSc=ENST1.1:c.5C>A",
    ) + &row_with(
        "rs2",
        "21:200",
        "T",
        "ENST1",
        "missense_variant",
        "IMPACT=MODERATE;STRAND=1",
    ) + &row_with(
        "rs3",
        "21:300",
        "A",
        "ENST1",
        "missense_variant",
        "IMPACT=HIGH;STRAND=1;HGVSc=ENST1.1:c.9C>A",
    ) + &row_with(
        "rs4",
        "21:400",
        "-",
        "ENST1",
        "frameshift_variant",
        "IMPACT=HIGH;STRAND=1;HGVSc=ENST1.1:c.12del",
    );
    let reference = write(&scratch.path, "ref.txt", &rows);
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!(
            "{VCF_HEADER}21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t200\trs2\tC\tT\t.\tPASS\t.\n21\t300\trs3\tC\tA\t.\tPASS\t.\n21\t400\trs4\tCA\tC\t.\tPASS\t.\n"
        ),
    );
    let exemplars = |m: &vep_golden_corpus::json::Object| {
        at(m, &["combinations", "missense_variant", "exemplar_records"]).as_i64()
    };
    let plain = under(&scratch.path, "plain");
    assert_eq!(
        select(&reference, &input, &plain, &["--covers", "symbol"]),
        clean(0, &summary(2, 2, "deletion 1, snv 1"))
    );
    let m = manifest(&plain);
    assert!(exemplars(&m) == Some(1) && !m.contains_key("family_columns"));
    let family = under(&scratch.path, "family");
    let key = ["--covers", "symbol", "--family-columns", "HGVSc,IMPACT="];
    assert_eq!(
        select(&reference, &input, &family, &key),
        clean(0, &summary(2, 4, "deletion 1, snv 3"))
    );
    let m = manifest(&family);
    assert_eq!(
        text_of(at(&m, &["family_columns"])),
        r#"["HGVSc", "IMPACT="]"#
    );
    assert_eq!(
        text_of(at(&m, &["family_strata"])),
        r#"{"HGVSc+;IMPACT=HIGH;deletion": 1, "HGVSc+;IMPACT=HIGH;snv": 1, "HGVSc+;IMPACT=MODERATE;snv": 1, "HGVSc-;IMPACT=MODERATE;snv": 1}"#
    );
    assert_eq!(exemplars(&m), Some(3));
    let dataset = format!("s01={reference}={input}");
    let check = |columns: &str, corpus: &str| {
        run(&[
            "check",
            "--dataset",
            &dataset,
            "--family-columns",
            columns,
            "--corpus",
            corpus,
        ])
    };
    assert_eq!(
        check("HGVSc,IMPACT=", &family),
        clean(
            0,
            &format!("{OBSERVED_TWO}family strata observed 4, covered by corpus 4, missing 0\n")
        )
    );
    // The three missense rows tie on every ordering key, so the seed alone decides which
    // one the plain corpus keeps: rs2, whose stratum is the one snv stratum it covers.
    assert_eq!(record_ids(&plain), ["rs2", "rs4"]);
    assert_eq!(
        check("HGVSc,IMPACT=", &plain),
        clean(
            1,
            &format!(
                "{OBSERVED_TWO}family strata observed 4, covered by corpus 2, missing 2\n  MISSING family stratum HGVSc+;IMPACT=HIGH;snv\n  MISSING family stratum HGVSc+;IMPACT=MODERATE;snv\n"
            )
        )
    );
    assert_eq!(
        check("HGVSc;IMPACT", &plain),
        failed(
            "--family-columns wants NAME, NAME=, custom_match=EXACT/OVERLAP or A/B, got 'HGVSc;IMPACT'"
        )
    );
    // The record cap keeps a stratum's sole exemplar as it keeps a combination's: four
    // strata survive a cap of three, where the same cap over four records chosen without
    // the key leaves three.
    let capped = under(&scratch.path, "capped");
    let mut capped_key = key.to_vec();
    capped_key.extend(["--max-records", "3"]);
    assert_eq!(select(&reference, &input, &capped, &capped_key).code, 0);
    assert_eq!(records(&manifest(&capped)).len(), 4);
    let plain_capped = under(&scratch.path, "plain-capped");
    let cut = ["--covers", "symbol", "--k", "3", "--max-records", "3"];
    assert_eq!(select(&reference, &input, &plain_capped, &cut).code, 0);
    assert_eq!(records(&manifest(&plain_capped)).len(), 3);
}

#[test]
fn custom_match_key_covers_every_class() {
    let scratch = scratch("strata-custom");
    let base = "IMPACT=MODERATE;STRAND=1";
    let missense = |name: &str, location: &str, allele: &str, extra: &str| {
        row_with(
            name,
            location,
            allele,
            "ENST1",
            "missense_variant",
            &format!("{base}{extra}"),
        )
    };
    let rows = missense(
        "rs1",
        "21:100",
        "G",
        ";CE=exact_1;CO=exact_1;CE_AF=0.10;CO_AF=0.10",
    ) + &missense("rs2", "21:200", "T", ";CO=overlap_1;CO_AF=0.20")
        + &missense("rs3", "21:300", "A", "")
        + &missense(
            "rs4",
            "21:400",
            "C",
            ";CE=multi_1;CO=multi_1;CE_AF=0.30;CO_AF=0.30,0.40",
        )
        + &row_with(
            "rs5",
            "21:500-900",
            "deletion",
            "ENST1",
            "transcript_ablation",
            "IMPACT=HIGH;STRAND=1;CO=exact_1,overlap_1;CO_AF=0.10,0.20",
        );
    let reference = write(&scratch.path, "ref.txt", &rows);
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!(
            "{VCF_HEADER}21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t200\trs2\tC\tT\t.\tPASS\t.\n21\t300\trs3\tC\tA\t.\tPASS\t.\n21\t400\trs4\tA\tC\t.\tPASS\t.\n21\t499\trs5\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=900\n"
        ),
    );
    let every = under(&scratch.path, "every");
    let key = [
        "--covers",
        "custom",
        "--family-columns",
        "custom_match=CE/CO",
    ];
    assert_eq!(
        select(&reference, &input, &every, &key),
        clean(0, &summary(2, 5, "snv 4, symbolic_deletion 1"))
    );
    let m = manifest(&every);
    assert_eq!(
        text_of(at(&m, &["family_columns"])),
        r#"["custom_match=CE/CO"]"#
    );
    assert_eq!(
        text_of(at(&m, &["family_strata"])),
        r#"{"custom_match=exact;snv": 1, "custom_match=miss;snv": 1, "custom_match=multi_allelic;snv": 1, "custom_match=overlap;snv": 1, "custom_match=sv_span;symbolic_deletion": 1}"#
    );
    let dataset = format!("s01={reference}={input}");
    let check = |corpus: &str| {
        run(&[
            "check",
            "--dataset",
            &dataset,
            "--family-columns",
            "custom_match=CE/CO",
            "--corpus",
            corpus,
        ])
    };
    assert_eq!(
        check(&every),
        clean(
            0,
            &format!("{OBSERVED_TWO}family strata observed 5, covered by corpus 5, missing 0\n")
        )
    );
    let plain = under(&scratch.path, "plain");
    assert_eq!(
        select(&reference, &input, &plain, &["--covers", "custom"]),
        clean(0, &summary(2, 2, "snv 1, symbolic_deletion 1"))
    );
    // The four missense rows tie on every ordering key, so the seed alone decides which one
    // the plain corpus keeps: rs4, the multi-allelic match. The three classes it does not
    // keep are the missing strata.
    assert_eq!(record_ids(&plain), ["rs4", "rs5"]);
    assert_eq!(
        check(&plain),
        clean(
            1,
            &format!(
                "{OBSERVED_TWO}family strata observed 5, covered by corpus 2, missing 3\n  MISSING family stratum custom_match=exact;snv\n  MISSING family stratum custom_match=miss;snv\n  MISSING family stratum custom_match=overlap;snv\n"
            )
        )
    );
}

#[test]
fn refseq_family_columns_stratify_exemplars() {
    // Six missense rows of one (class, strand) stratum differ in the RefSeq and merged
    // cache columns alone.
    let scratch = scratch("strata-refseq");
    let matched = "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt";
    let missense = |name: &str, location: &str, allele: &str, feature: &str, extra: &str| {
        let extra = format!("IMPACT=MODERATE;STRAND=1;{extra}");
        row_with(name, location, allele, feature, "missense_variant", &extra)
    };
    let rows = missense(
        "rs1",
        "21:100",
        "G",
        "NM_1.1",
        &format!("{matched};BAM_EDIT=OK;GIVEN_REF=A;USED_REF=A"),
    ) + &missense(
        "rs2",
        "21:200",
        "T",
        "NM_1.1",
        "REFSEQ_MATCH=rseq_mrna_nonmatch,rseq_5p_mismatch;BAM_EDIT=OK;GIVEN_REF=C;USED_REF=T",
    ) + &missense(
        "rs3",
        "21:300",
        "A",
        "NM_1.1",
        &format!("{matched};BAM_EDIT=FAILED;GIVEN_REF=C;USED_REF=C"),
    ) + &missense(
        "rs4",
        "21:400",
        "C",
        "NM_1.1",
        &format!("{matched};BAM_EDIT=OK;GIVEN_REF=A;USED_REF=A;REFSEQ_OFFSET=-3"),
    ) + &missense(
        "rs5",
        "21:500",
        "G",
        "NM_1.1",
        &format!("{matched};SOURCE=RefSeq;BAM_EDIT=OK;GIVEN_REF=A;USED_REF=A"),
    ) + &missense(
        "rs6",
        "21:600",
        "T",
        "ENST1",
        "SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A",
    );
    let reference = write(&scratch.path, "ref.txt", &rows);
    let input = write(
        &scratch.path,
        "in.vcf",
        &format!(
            "{VCF_HEADER}21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t200\trs2\tC\tT\t.\tPASS\t.\n21\t300\trs3\tC\tA\t.\tPASS\t.\n21\t400\trs4\tA\tC\t.\tPASS\t.\n21\t500\trs5\tA\tG\t.\tPASS\t.\n21\t600\trs6\tA\tT\t.\tPASS\t.\n"
        ),
    );
    let key = "REFSEQ_MATCH=,SOURCE=,BAM_EDIT=,REFSEQ_OFFSET,GIVEN_REF/USED_REF";
    let every = under(&scratch.path, "every");
    assert_eq!(
        select(
            &reference,
            &input,
            &every,
            &["--covers", "refseq", "--family-columns", key]
        ),
        clean(0, &summary(1, 6, "snv 6"))
    );
    let m = manifest(&every);
    let columns: Vec<&str> = key.split(',').collect();
    assert_eq!(at(&m, &["family_columns"]), &Value::strings(&columns));
    let strata = [
        "REFSEQ_MATCH=-;SOURCE=Ensembl;BAM_EDIT=-;REFSEQ_OFFSET-;GIVEN_REF/USED_REF=same;snv",
        "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;SOURCE=-;BAM_EDIT=FAILED;REFSEQ_OFFSET-;GIVEN_REF/USED_REF=same;snv",
        "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;SOURCE=-;BAM_EDIT=OK;REFSEQ_OFFSET+;GIVEN_REF/USED_REF=same;snv",
        "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;SOURCE=-;BAM_EDIT=OK;REFSEQ_OFFSET-;GIVEN_REF/USED_REF=same;snv",
        "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;SOURCE=RefSeq;BAM_EDIT=OK;REFSEQ_OFFSET-;GIVEN_REF/USED_REF=same;snv",
        "REFSEQ_MATCH=rseq_mrna_nonmatch,rseq_5p_mismatch;SOURCE=-;BAM_EDIT=OK;REFSEQ_OFFSET-;GIVEN_REF/USED_REF=differ;snv",
    ];
    let counted: Vec<String> = strata
        .iter()
        .map(|stratum| format!("\"{stratum}\": 1"))
        .collect();
    assert_eq!(
        text_of(at(&m, &["family_strata"])),
        format!("{{{}}}", counted.join(", "))
    );
    let exemplars = at(
        &m,
        &["combinations", "missense_variant", "exemplar_records"],
    );
    assert_eq!(exemplars.as_i64(), Some(6));
    let dataset = format!("s01={reference}={input}");
    let check = |corpus: &str| {
        run(&[
            "check",
            "--dataset",
            &dataset,
            "--family-columns",
            key,
            "--corpus",
            corpus,
        ])
    };
    assert_eq!(
        check(&every),
        clean(
            0,
            &format!("{OBSERVED_ONE}family strata observed 6, covered by corpus 6, missing 0\n")
        )
    );
    let plain = under(&scratch.path, "plain");
    assert_eq!(
        select(&reference, &input, &plain, &["--covers", "refseq"]),
        clean(0, &summary(1, 1, "snv 1"))
    );
    // Seed 1 keeps rs3 of the six tied rows, so the other five strata are missing.
    assert_eq!(record_ids(&plain), ["rs3"]);
    let missing: Vec<String> = [0, 2, 3, 4, 5]
        .iter()
        .map(|place| format!("  MISSING family stratum {}\n", strata[*place]))
        .collect();
    assert_eq!(
        check(&plain),
        clean(
            1,
            &format!(
                "{OBSERVED_ONE}family strata observed 6, covered by corpus 1, missing 5\n{}",
                missing.concat()
            )
        )
    );
    assert!(record_by_id(&manifest(&plain), "rs3")
        .get("combinations")
        .is_some());
    assert!(fs::metadata(under(std::path::Path::new(&plain), "variants.vcf")).is_ok());
}
